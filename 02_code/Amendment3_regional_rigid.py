"""
Amendment3_regional_rigid.py
===========================================================================
PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - Amendment 3 (LOCKED specification v2.1, 1 Nov 2026).
Patient-specific regional rigid refinement of HEAD and PELVIS, composed on the current adopted transform T_A
(Amendment 2 for 97 patients; supplied rigid.tfm for train_0018, train_0049, train_0067). Tumour-blind.

Stage implemented here: numerical fit and PROVISIONAL selection (the blinded visual veto is a separate stage).
  T_A(p)  = FDG/rigid.tfm(p) - t   (adopted)   |  FDG/rigid.tfm(p) (fallback)       [point mapping FDG -> PSMA]
  T_r(p)  = dT_r(T_A(p)),  dT_r rigid, estimated in PSMA space by SYMMETRIC TRIMMED ICP:
            nearest neighbours in both directions; closest 80% kept PER DIRECTION; equal total weight (0.5) per
            direction; weighted Kabsch; stop < 0.01 mm and < 0.01 deg, or 100 iterations; start = identity (= T_A).
  one-way trimmed ICP (FDG -> PSMA) = bias diagnostic (difference > 1 deg or > 1 mm flagged).
Roles:   head   fit = complete eligible skull surface | select = brain (centroid, surface distance) |
                diagnostics = skull upper/lower spatial cross-validation (stability only), C1-C2 (supporting only)
         pelvis fit = hip bones only (labels 77/78; coded check: no femur/sacrum/S1 voxel) | select = sacrum + S1 |
                report = femoral heads (top 40 mm of femur label), iliac vessels / gluteal + iliopsoas surface distance
Coverage: common field of view (10 mm z-margin); skull >= 2000 points & >= 60 mm span; pelvis >= 2000 points per hip,
          both sides. Otherwise NOT ESTIMABLE -> Amendment 2.
Caps (reject, never clip): HEAD rotation > 20 deg or translation > 50 mm (Amendment 3a, provisional);
      PELVIS rotation > 10 deg or translation > 20 mm (unchanged). Translation = displacement of the PSMA fitting-surface
      centroid. Maximum 300 ICP iterations (Amendment 3a); convergence mandatory (more iterations improve the chance of
      convergence, they do not guarantee correct registration).
Selection (R4): adopt only if improvement of the selection centroid distance > max(1.0 mm, 2 x bootstrap SD) AND the
          selection surface distance is not worse by > 0.5 mm. Bootstrap: 200 replicates, each resampling the fitting
          points (both CTs, with replacement) and REPEATING the full symmetric ICP and the selection calculation
          (fitting uncertainty, not total measurement noise); computed only when the other conditions already pass.
Safety monitor (R7): among PROVISIONALLY adopting patients with evaluable reporting structures, > 10% worse by > 2 mm
          -> further adoption in that region suspended pending review (reported; nothing replaced automatically).
Thresholds are engineering choices, not clinically validated standards.

Usage (Anaconda Prompt, env totalseg):
  python Amendment3_regional_rigid.py --selftest
  python Amendment3_regional_rigid.py --preflight 3          implementation checks only (nothing written)
  python Amendment3_regional_rigid.py --overlay train_0001 train_0002 train_0003
                                                             diagnostic overlays: Amendment 2 vs proposed correction
  python Amendment3_regional_rigid.py                        definitive numerical stage (all 100 patients)
Outputs (06_output): Amendment3_fit_<run>.csv (one row per patient x region), Amendment3_transforms_<run>.csv
  (dT matrices), Amendment3_summary_<run>.csv, Amendment3_manifest_<run>.json, Amendment3_report_<run>.txt
Requires: Python 3.9+, numpy, pandas, scipy, SimpleITK.
===========================================================================
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
A2_PATIENT = "Audit1C_e_per_patient_20261031_141128.csv"
FALLBACKS = {"train_0018", "train_0049", "train_0067"}

LAB = dict(skull=91, brain=90, C1=50, C2=49, hip_L=77, hip_R=78, sacrum=25, S1=26, femur_L=75, femur_R=76,
           iliac_art_L=65, iliac_art_R=66, iliac_vein_L=67, iliac_vein_R=68,
           glut_max_L=80, glut_max_R=81, glut_med_L=82, glut_med_R=83, glut_min_L=84, glut_min_R=85,
           iliopsoas_L=88, iliopsoas_R=89)
TRIM = 0.80
MAX_ITER = 300          # Amendment 3a (1 Nov 2026): 100 -> 300; convergence remains mandatory
TOL_MM, TOL_DEG = 0.01, 0.01
FOV_MARGIN = 10.0
MIN_POINTS = 2000
MIN_SKULL_SPAN = 60.0
CAPS = {"head": (20.0, 50.0), "pelvis": (10.0, 20.0)}   # (deg, mm). Amendment 3a: head 10/20 -> 20/50 (provisional);
                                                        # pelvis unchanged. Engineering choices, not validated limits.
CAP_DEG, CAP_MM = CAPS["pelvis"]                        # default (synthetic tests)
MIN_IMPROVE = 1.0
SURF_WORSE = 0.5
N_BOOT = 200
SEED = 20261101
MAX_PTS = 4000
FLAG_ONEWAY_DEG, FLAG_ONEWAY_MM = 1.0, 1.0
FEMHEAD_MM = 40.0
SAFETY_FRAC, SAFETY_MM = 0.10, 2.0


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ rigid algebra ------------------------------
def kabsch_weighted(X, Y, w):
    """R, c minimising sum w ||R x + c - y||^2 (proper rotation)."""
    w = w / w.sum()
    mx, my = w @ X, w @ Y
    H = (X - mx).T @ ((Y - my) * w[:, None])
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, my - R @ mx


def rot_deg(R):
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))


def apply(R, c, P):
    return P @ R.T + c


def trimmed_pairs(src, tgt_tree, tgt):
    d, j = tgt_tree.query(src, k=1)
    keep = np.argsort(d)[: max(3, int(np.floor(TRIM * len(d))))]
    return src[keep], tgt[j[keep]], d[keep]


def icp(F0, P, symmetric=True):
    """F0: FDG points already mapped by T_A (PSMA space); P: PSMA points. Returns R, c, iterations, converged.
    Symmetric: pairs from both directions, trimmed separately, each direction total weight 0.5."""
    R, c = np.eye(3), np.zeros(3)
    treeP = cKDTree(P)
    for it in range(1, MAX_ITER + 1):
        F = apply(R, c, F0)
        xs, ys, _ = trimmed_pairs(F, treeP, P)                     # FDG -> PSMA
        X, Y = [xs], [ys]
        W = [np.full(len(xs), 0.5 / len(xs)) if symmetric else np.full(len(xs), 1.0 / len(xs))]
        if symmetric:
            treeF = cKDTree(F)
            qs, fs, _ = trimmed_pairs(P, treeF, F)                 # PSMA -> FDG: source = moving FDG point
            X.append(fs); Y.append(qs); W.append(np.full(len(qs), 0.5 / len(qs)))
        Ri, ci = kabsch_weighted(np.vstack(X), np.vstack(Y), np.concatenate(W))
        R, c = Ri @ R, Ri @ c + ci
        step_mm = float(np.max(np.linalg.norm(apply(Ri, ci, F) - F, axis=1))) if len(F) else 0.0
        if step_mm < TOL_MM and rot_deg(Ri) < TOL_DEG:
            return R, c, it, True
    return R, c, MAX_ITER, False


def sym_surface_distance(A, B):
    if len(A) == 0 or len(B) == 0:
        return np.nan
    return float((cKDTree(B).query(A)[0].mean() + cKDTree(A).query(B)[0].mean()) / 2)


def subsample(P, n, rng):
    return P if len(P) <= n else P[rng.choice(len(P), n, replace=False)]


# ------------------------------ geometry from images ------------------------------
class Geom:
    def __init__(self, img):
        self.o = np.array(img.GetOrigin(), float)
        self.s = np.array(img.GetSpacing(), float)
        self.D = np.array(img.GetDirection(), float).reshape(3, 3)
        self.size = np.array(img.GetSize(), int)

    def phys(self, idx_xyz):
        return self.o + (np.asarray(idx_xyz, float) * self.s) @ self.D.T

    def z_range(self):
        cn = np.array([[i, j, k] for i in (0, self.size[0] - 1) for j in (0, self.size[1] - 1) for k in (0, self.size[2] - 1)])
        z = self.phys(cn)[:, 2]
        return float(z.min()), float(z.max())

    def corners(self):
        cn = np.array([[i, j, k] for i in (0, self.size[0] - 1) for j in (0, self.size[1] - 1) for k in (0, self.size[2] - 1)])
        return self.phys(cn)


def label_sets(arr_zyx, labels, g):
    """surface points (6-neighbour boundary) and all voxel points for the union of labels."""
    m = np.isin(arr_zyx, labels)
    if not m.any():
        return np.zeros((0, 3)), np.zeros((0, 3))
    sl = ndimage.find_objects(m.astype(np.uint8))[0]
    sub = m[sl]
    surf = sub & ~ndimage.binary_erosion(sub, structure=ndimage.generate_binary_structure(3, 1), border_value=0)
    out = []
    for mm in (surf, sub):
        zz, yy, xx = np.nonzero(mm)
        idx = np.stack([xx + sl[2].start, yy + sl[1].start, zz + sl[0].start], axis=1)
        out.append(g.phys(idx))
    return out[0], out[1]


def affine_of(T):
    """Exact affine (A, b) of a SimpleITK rigid/affine transform from 4 points; verified on random points."""
    b = np.array(T.TransformPoint((0.0, 0.0, 0.0)))
    A = np.stack([np.array(T.TransformPoint(tuple(e))) - b for e in np.eye(3)], axis=1)
    rng = np.random.default_rng(1)
    for p in rng.uniform(-500, 500, (5, 3)):
        if np.linalg.norm(np.array(T.TransformPoint(tuple(p))) - (A @ p + b)) > 1e-6:
            raise RuntimeError("transform is not affine to 1e-6 mm")
    return A, b


# ------------------------------ per-region evaluation ------------------------------
def in_fov(P, zlo, zhi):
    return P[(P[:, 2] >= zlo) & (P[:, 2] <= zhi)] if len(P) else P


def centroid_dist(Pvox, Fvox_m, R=np.eye(3), c=np.zeros(3)):
    if len(Pvox) == 0 or len(Fvox_m) == 0:
        return np.nan
    return float(np.linalg.norm(apply(R, c, Fvox_m).mean(axis=0) - Pvox.mean(axis=0)))


def fit_region(fitP, fitF, sel, rng, boot=True, caps=None):
    """fitP/fitF: fitting surfaces (PSMA; FDG mapped by T_A). sel: list of (Pvox, Fvox_m, Psurf, Fsurf_m) selection
    structures. Returns dict of fit diagnostics and selection metrics."""
    out = {}
    R, c, it, conv = icp(fitF, fitP, symmetric=True)
    R1, c1, it1, conv1 = icp(fitF, fitP, symmetric=False)
    c0 = fitP.mean(axis=0)
    out.update(converged=conv, iterations=it, rot_deg=rot_deg(R), trans_mm=float(np.linalg.norm(apply(R, c, c0[None])[0] - c0)),
               oneway_converged=conv1, oneway_diff_deg=rot_deg(R1 @ R.T),
               oneway_diff_mm=float(np.linalg.norm(apply(R1, c1, c0[None])[0] - apply(R, c, c0[None])[0])))
    out["oneway_flag"] = bool(out["oneway_diff_deg"] > FLAG_ONEWAY_DEG or out["oneway_diff_mm"] > FLAG_ONEWAY_MM)
    cap_deg, cap_mm = caps if caps else (CAP_DEG, CAP_MM)
    out["cap_deg"], out["cap_mm"] = cap_deg, cap_mm
    out["within_caps"] = bool(out["rot_deg"] <= cap_deg and out["trans_mm"] <= cap_mm)

    def sel_metrics(Rx, cx):
        cds = [centroid_dist(Pv, Fv, Rx, cx) for Pv, Fv, _, _ in sel]
        sds = [sym_surface_distance(Ps, apply(Rx, cx, Fs)) for _, _, Ps, Fs in sel]
        return float(np.nanmean(cds)), float(np.nanmean(sds))
    cA, sA = sel_metrics(np.eye(3), np.zeros(3))
    cR, sR = sel_metrics(R, c)
    out.update(sel_centroid_A2=cA, sel_centroid_A3=cR, sel_improvement=cA - cR, sel_surface_A2=sA, sel_surface_A3=sR,
               sel_surface_worse=sR - sA)
    pre_ok = conv and out["within_caps"] and out["sel_improvement"] > MIN_IMPROVE and out["sel_surface_worse"] <= SURF_WORSE
    out["boot_sd"] = np.nan
    out["boot_n"] = 0
    if pre_ok and boot:
        imps = []
        for b in range(N_BOOT):
            iP = rng.integers(0, len(fitP), len(fitP))
            iF = rng.integers(0, len(fitF), len(fitF))
            Rb, cb, _, _ = icp(fitF[iF], fitP[iP], symmetric=True)          # full refit (author safeguard)
            imps.append(cA - sel_metrics(Rb, cb)[0])                       # selection recomputed with replicate
        out["boot_sd"] = float(np.std(imps, ddof=1))
        out["boot_n"] = len(imps)
    thr = max(MIN_IMPROVE, 2 * out["boot_sd"]) if np.isfinite(out["boot_sd"]) else np.inf
    out["threshold_mm"] = thr if np.isfinite(thr) else np.nan
    out["provisional_adopt"] = bool(pre_ok and boot and out["sel_improvement"] > thr)
    out["R"], out["c"] = R, c
    return out


def decision_reason(o, estimable, why_not):
    if not estimable:
        return f"NOT ESTIMABLE ({why_not}) -> Amendment 2"
    if not o["converged"]:
        return "not converged -> Amendment 2"
    if not o["within_caps"]:
        return "exceeds caps (rejected, not clipped) -> Amendment 2"
    if o["sel_improvement"] <= MIN_IMPROVE:
        return "selection improvement <= 1.0 mm -> Amendment 2"
    if o["sel_surface_worse"] > SURF_WORSE:
        return "selection surface distance worse > 0.5 mm -> Amendment 2"
    if not o["provisional_adopt"]:
        return "improvement <= 2 x bootstrap SD -> Amendment 2"
    return "PROVISIONALLY ADOPTED (visual veto pending)"


# ------------------------------ synthetic self-test ------------------------------
def rot_xyz(ax, ay, az):
    ax, ay, az = np.radians([ax, ay, az])
    Rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def synth_surface(rng, n=6000):
    """asymmetric closed surface (ellipsoid with bumps) ~ bone-like, mm."""
    u = rng.uniform(0, 2 * np.pi, n); v = np.arccos(rng.uniform(-1, 1, n))
    r = 1 + 0.25 * np.sin(3 * u) * np.sin(2 * v) + 0.15 * np.cos(5 * v)
    return np.stack([90 * r * np.sin(v) * np.cos(u), 60 * r * np.sin(v) * np.sin(u), 50 * r * np.cos(v)], axis=1)


def selftest():
    print("SELF-TEST (synthetic only; no study data)")
    rng = np.random.default_rng(0)
    res = []
    P = synth_surface(rng)
    Rt, ct = rot_xyz(4, -3, 5), np.array([6.0, -4.0, 3.0])
    # FDG = misaligned copy (inverse of truth) so ICP must recover (Rt, ct): Rt @ F + ct = P
    F = (P - ct) @ Rt                          # R^T (P - c)
    R, c, it, conv = icp(F, P)
    res.append(("symmetric ICP recovers known rotation/translation", conv and rot_deg(R @ Rt.T) < 0.05 and np.linalg.norm(c - ct) < 0.1))
    # 20% corrupted points + unequal partial coverage
    F2 = F.copy()
    bad = rng.choice(len(F2), int(0.2 * len(F2)), replace=False)
    F2[bad] += rng.normal(0, 15, (len(bad), 3))
    Pp = P[P[:, 2] > -20]                      # PSMA surface truncated (unequal coverage)
    R2, c2, _, conv2 = icp(F2, Pp)
    res.append(("robust to 20% corrupted points and partial coverage (< 0.5 deg, < 1 mm)",
                conv2 and rot_deg(R2 @ Rt.T) < 0.5 and np.linalg.norm(c2 - ct) < 1.0))
    R1, c1, _, _ = icp(F2, Pp, symmetric=False)
    res.append(("one-way solution computed for the diagnostic", np.isfinite(rot_deg(R1))))
    # equal weighting per direction: duplicating one set must not change the symmetric solution materially
    Rd, cd, _, _ = icp(np.vstack([F2, F2]), Pp)
    res.append(("equal direction weighting: duplicating FDG points leaves solution unchanged",
                rot_deg(Rd @ R2.T) < 0.05 and np.linalg.norm(cd - c2) < 0.1))
    X = rng.normal(size=(10, 3)); Y = X + 1
    w = np.concatenate([np.full(7, 0.5 / 7), np.full(3, 0.5 / 3)])
    res.append(("weights sum to 0.5 per direction", abs(w[:7].sum() - 0.5) < 1e-12 and abs(w[7:].sum() - 0.5) < 1e-12))
    Rk, ck = kabsch_weighted(X, Y, w)
    res.append(("weighted Kabsch exact on noiseless data", np.allclose(Rk, np.eye(3), atol=1e-10) and np.allclose(ck, 1)))
    # caps
    Rbig = rot_xyz(15, 0, 0)
    Fb = (P - 0) @ Rbig
    sel = [(P[:500], (P[:500] - 0) @ Rbig, P[:500], (P[:500]) @ Rbig)]
    ob = fit_region(P, Fb, sel, rng, boot=False)
    res.append(("rotation > 10 deg rejected (not clipped)", (not ob["within_caps"]) and ob["rot_deg"] > 10 and not ob["provisional_adopt"]))
    # selection rule with full-refit bootstrap: a genuine 6 mm / 5 deg misalignment with coupled selection structure
    selP = P[(P[:, 0] > 40)]                    # a coupled 'selection' structure (subset)
    selF = (selP - ct) @ Rt
    o = fit_region(subsample(P, 1500, rng), subsample(F, 1500, rng), [(selP, selF, selP, selF)], np.random.default_rng(5))
    res.append(("bootstrap repeats full refits (boot_n = 200, SD > 0) and adopts a clear improvement",
                o["boot_n"] == N_BOOT and o["boot_sd"] > 0 and o["provisional_adopt"]))
    o0 = fit_region(subsample(P, 1500, rng), subsample(P, 1500, rng), [(selP, selP, selP, selP)], rng)
    res.append(("already aligned: no adoption, no bootstrap", not o0["provisional_adopt"] and o0["boot_n"] == 0))
    res.append(("decision reason text", decision_reason(o0, True, "").endswith("Amendment 2")
                and decision_reason(o, True, "").startswith("PROVISIONALLY")
                and decision_reason(o, False, "too few points").startswith("NOT ESTIMABLE")))
    # label hygiene: the pelvic fitting set never contains femur / sacrum / S1 labels
    arr = np.zeros((20, 20, 20), np.uint8); arr[2:8, 2:8, 2:8] = LAB["hip_L"]; arr[10:16, 2:8, 2:8] = LAB["femur_L"]
    arr[2:8, 10:16, 2:8] = LAB["sacrum"]
    class G0:  # identity geometry
        def phys(self, idx):
            return np.asarray(idx, float)
    s_hip, v_hip = label_sets(arr, [LAB["hip_L"], LAB["hip_R"]], G0())
    res.append(("label hygiene: hip fitting set excludes femur and sacrum voxels",
                len(v_hip) == 216 and (v_hip[:, 2] < 8).all() and (v_hip[:, 1] < 8).all()))
    res.append(("surface extraction gives boundary voxels only", 0 < len(s_hip) < len(v_hip)))
    for n, ok in res:
        print(f"  [{'PASS' if ok else 'FAIL'}] {n}")
    allok = all(ok for _, ok in res)
    print(f"SELF-TEST {'PASSED' if allok else 'FAILED'} ({sum(ok for _, ok in res)}/{len(res)})")
    return allok


def render_overlay(cid, region, dP, dF, A, b, t, R, c, row, odir):
    """Diagnostic overlay (author validation of the proposed correction; not a selection step): PSMA CT with PSMA
    contours (green) and FDG contours resampled (magenta) under Amendment 2 (left) and under the proposed
    regional correction (right). Head: skull, brain, C1-C2; pelvis: hips, sacrum, S1, femora."""
    import SimpleITK as sitk
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labs = [LAB["skull"], LAB["brain"], LAB["C1"], LAB["C2"]] if region == "head" else \
        [LAB["hip_L"], LAB["hip_R"], LAB["sacrum"], LAB["S1"], LAB["femur_L"], LAB["femur_R"]]
    Rr = sitk.ReadTransform(str(dP / "rigid.tfm"))
    RA = sitk.CompositeTransform(3)
    RA.AddTransform(Rr)
    if cid not in FALLBACKS:
        RA.AddTransform(sitk.TranslationTransform(3, tuple(float(v) for v in t)))
    Ri = R.T
    inv = sitk.AffineTransform(3)
    inv.SetMatrix(tuple(Ri.ravel()))
    inv.SetTranslation(tuple(float(v) for v in (-Ri @ c)))
    RB = sitk.CompositeTransform(3)
    RB.AddTransform(RA)
    RB.AddTransform(inv)                       # applied first: dT^-1, then T_A^-1
    # composition check: RB(dT(T_A(p))) == p
    for p in np.random.default_rng(3).uniform(-200, 200, (3, 3)):
        q = R @ ((A @ p + b) - t) + c
        if np.linalg.norm(np.array(RB.TransformPoint(tuple(q))) - p) > 1e-4:
            raise RuntimeError(f"{cid}: overlay composition check failed")
    tsP = sitk.Cast(sitk.ReadImage(str(dP / "totseg_24.nii.gz")), sitk.sitkUInt8)
    tsF = sitk.Cast(sitk.ReadImage(str(dF / "totseg_24.nii.gz")), sitk.sitkUInt8)
    ct = sitk.Resample(sitk.ReadImage(str(dP / "CT.nii.gz")), tsP, sitk.Transform(), sitk.sitkLinear, -1000, sitk.sitkFloat32)
    rA = sitk.Resample(tsF, tsP, RA, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    rB = sitk.Resample(tsF, tsP, RB, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    # IMAGE-LEVEL check (author safeguard): the centroid of each FDG label as it appears in the RESAMPLED image must
    # equal the centroid of the same FDG voxels POINT-MAPPED (Amendment 2: T_A; candidate: dT o T_A), restricted to
    # voxels landing inside the PSMA grid. Tolerance = largest PSMA voxel dimension (nearest-neighbour resampling).
    def phys(im, kji):
        D = np.array(im.GetDirection()).reshape(3, 3); S = np.array(im.GetSpacing())
        return np.array(im.GetOrigin()) + (kji[:, ::-1] * S) @ D.T
    def inside(im, q):
        D = np.array(im.GetDirection()).reshape(3, 3); S = np.array(im.GetSpacing())
        ci = np.linalg.solve(D * S, (q - np.array(im.GetOrigin())).T).T
        return np.all((ci > 0.5) & (ci < np.array(im.GetSize()) - 1.5), axis=1)
    aF_, tol, checks = sitk.GetArrayFromImage(tsF), max(tsP.GetSpacing()), []
    for lab in labs:
        idx = np.argwhere(aF_ == lab)
        if len(idx) < 50:
            continue
        qA = (phys(tsF, idx.astype(float)) @ A.T + b) - t
        for nm, q, res_im in (("A2", qA, rA), ("A3a", qA @ R.T + c, rB)):
            k = inside(tsP, q)
            ra = sitk.GetArrayFromImage(res_im)
            ridx = np.argwhere(ra == lab)
            if k.sum() < 50 or len(ridx) < 50:
                continue
            err = float(np.linalg.norm(phys(tsP, ridx.astype(float)).mean(0) - q[k].mean(0)))
            checks.append((nm, lab, err))
    worst = max((e for _, _, e in checks), default=np.nan)
    status = "PASS" if checks and worst <= tol else "FAIL"
    print(f"    {cid} {region}: image-level transform check {status} (worst centroid discrepancy {worst:.2f} mm, "
          f"tolerance {tol:.2f} mm, {len(checks)} label/transform checks)")
    L = lambda im: sitk.GetArrayFromImage(sitk.DICOMOrient(im, "LPS"))
    ctO, pO, aO, bO = L(ct), L(tsP), L(rA), L(rB)
    sp = sitk.DICOMOrient(tsP, "LPS").GetSpacing()
    m = np.isin(pO, labs)
    if not m.any():
        return
    zz, yy, xx = np.nonzero(m)
    g = [int(np.ceil(20 / s)) for s in (sp[2], sp[1], sp[0])]
    lo = [max(0, v.min() - k) for v, k in zip((zz, yy, xx), g)]
    hi = [min(n, v.max() + 1 + k) for v, k, n in zip((zz, yy, xx), g, pO.shape)]
    z0, y0, x0 = int(np.median(zz)), int(np.median(yy)), int(np.median(xx))
    zs, ys, xs = slice(lo[0], hi[0]), slice(lo[1], hi[1]), slice(lo[2], hi[2])
    views = [("coronal", np.s_[zs, y0, xs], sp[0], sp[2], True), ("sagittal", np.s_[zs, ys, x0], sp[1], sp[2], True),
             ("axial", np.s_[z0, ys, xs], sp[0], sp[1], False)]
    fig, axs = plt.subplots(3, 2, figsize=(12, 16))
    for j, (name, res) in enumerate([("AMENDMENT 2", aO), (f"AMENDMENT 3a CANDIDATE "
                                     f"({row['rot_deg']:.1f} deg, {row['trans_mm']:.1f} mm, "
                                     f"{'converged' if row['converged'] else 'NOT converged'}, "
                                     f"{'within caps' if row['within_caps'] else 'EXCEEDS caps'})", bO)]):
        for i, (ttl, sl, w, h, flip) in enumerate(views):
            ax = axs[i, j]
            f = np.flipud if flip else (lambda a_: a_)
            img = f(ctO[sl])
            ext = (0, img.shape[1] * w, 0, img.shape[0] * h)
            ax.imshow(img, cmap="gray", vmin=-200, vmax=1000, extent=ext, origin="upper")
            for arr, col in ((pO, "lime"), (res, "magenta")):
                a1 = f(np.isin(arr[sl], labs)).astype(float)
                if a1.any():
                    ax.contour(a1, levels=[0.5], colors=col, linewidths=0.8, extent=ext, origin="upper")
            ax.set_title(f"{name} | {ttl}", fontsize=8.5)
    fig.suptitle(f"{cid} {region} | LEFT: Amendment 2 | RIGHT: Amendment 3a candidate | background PSMA CT (reference, "
                 f"not transformed) | green PSMA labels, magenta FDG labels resampled\nimage-level transform check {status} "
                 f"(worst {worst:.2f} mm, tol {tol:.2f} mm) | diagnostic only, no decision", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(odir / f"{cid}_{region}_A2_vs_A3a_candidate.png", dpi=100)
    plt.close(fig)


# ------------------------------ study run ------------------------------
def run(preflight=None, overlay_cases=None):
    import SimpleITK as sitk
    t0 = time.time()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    PP = pd.read_csv(OUT_DIR / A2_PATIENT).set_index("case_id")
    fb = set(PP.index[PP["would_fallback_if_adopted"].astype(str).str.lower() == "true"])
    if fb != FALLBACKS:
        sys.exit(f"fallback set mismatch {sorted(fb)}")
    cases = sorted(p.name for p in DATA_ROOT.iterdir() if p.is_dir() and p.name.startswith("train_"))
    if preflight:
        cases = cases[:preflight]
    if overlay_cases:
        cases = [c for c in cases if c in set(overlay_cases)]
        odir = OUT_DIR / f"Amendment3_preflight_overlays_{run_id}"
        odir.mkdir()
    rows, trows = [], []
    rng = np.random.default_rng(SEED)
    for i, cid in enumerate(cases, 1):
        dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
        T = sitk.ReadTransform(str(dF / "rigid.tfm"))
        A, b = affine_of(T)
        t = np.zeros(3) if cid in FALLBACKS else PP.loc[cid, ["t_x_LPS", "t_y_LPS", "t_z_LPS"]].to_numpy(float)
        TA = lambda Q: (Q @ A.T + b) - t if len(Q) else Q
        # composition check against the SimpleITK composite used elsewhere (T'(p) = T(p) - t)
        for p in np.random.default_rng(2).uniform(-300, 300, (3, 3)):
            if np.linalg.norm(TA(p[None])[0] - (np.array(T.TransformPoint(tuple(p))) - t)) > 1e-6:
                sys.exit(f"{cid}: composition check failed")
        imP = sitk.ReadImage(str(dP / "totseg_24.nii.gz")); imF = sitk.ReadImage(str(dF / "totseg_24.nii.gz"))
        gP, gF = Geom(imP), Geom(imF)
        aP, aF = sitk.GetArrayFromImage(imP), sitk.GetArrayFromImage(imF)
        zP = gP.z_range(); zF = TA(gF.corners())[:, 2]
        zlo, zhi = max(zP[0], zF.min()) + FOV_MARGIN, min(zP[1], zF.max()) - FOV_MARGIN

        def get(labels):
            sP, vP = label_sets(aP, labels, gP)
            sF, vF = label_sets(aF, labels, gF)
            return in_fov(sP, zlo, zhi), in_fov(vP, zlo, zhi), in_fov(TA(sF), zlo, zhi), in_fov(TA(vF), zlo, zhi)

        for region in ("head", "pelvis"):
            row = dict(case_id=cid, region=region, fallback_patient=cid in FALLBACKS, fov_zlo=zlo, fov_zhi=zhi)
            if region == "head":
                fsP, _, fsF, _ = get([LAB["skull"]])
                span = min(np.ptp(fsP[:, 2]) if len(fsP) else 0, np.ptp(fsF[:, 2]) if len(fsF) else 0)
                ok = len(fsP) >= MIN_POINTS and len(fsF) >= MIN_POINTS and span >= MIN_SKULL_SPAN
                why = f"skull points P {len(fsP)}, F {len(fsF)}, span {span:.0f} mm"
                bsP, bvP, bsF, bvF = get([LAB["brain"]])
                sel = [(bvP, bvF, subsample(bsP, MAX_PTS, rng), subsample(bsF, MAX_PTS, rng))] if len(bvP) and len(bvF) else []
                if not sel:
                    ok, why = False, why + "; brain (selection) unavailable"
                row.update(n_fit_P=len(fsP), n_fit_F=len(fsF), fit_span_mm=span)
                fitP, fitF = subsample(fsP, MAX_PTS, rng), subsample(fsF, MAX_PTS, rng)
            else:
                hL = get([LAB["hip_L"]]); hR = get([LAB["hip_R"]])
                ok = min(len(hL[0]), len(hL[2]), len(hR[0]), len(hR[2])) >= MIN_POINTS
                why = f"hip points L P{len(hL[0])}/F{len(hL[2])}, R P{len(hR[0])}/F{len(hR[2])} (both sides >= {MIN_POINTS})"
                fitP = subsample(np.vstack([hL[0], hR[0]]), MAX_PTS, rng)
                fitF = subsample(np.vstack([hL[2], hR[2]]), MAX_PTS, rng)
                sel = []
                for lab in ("sacrum", "S1"):
                    sP, vP, sF, vF = get([LAB[lab]])
                    if len(vP) and len(vF):
                        sel.append((vP, vF, subsample(sP, MAX_PTS, rng), subsample(sF, MAX_PTS, rng)))
                if not sel:
                    ok, why = False, why + "; sacrum/S1 (selection) unavailable"
                row.update(n_fit_P=min(len(hL[0]), len(hR[0])), n_fit_F=min(len(hL[2]), len(hR[2])), fit_span_mm=np.nan)
            o = None
            if ok:
                o = fit_region(fitP, fitF, sel, rng, boot=not overlay_cases, caps=CAPS[region])
                R, c = o.pop("R"), o.pop("c")
                row.update(o)
                trows.append(dict(case_id=cid, region=region, **{f"R{i}{j}": R[i, j] for i in range(3) for j in range(3)},
                                  **{f"c{k}": c[k] for k in range(3)}))
                if overlay_cases:
                    render_overlay(cid, region, dP, dF, A, b, t, R, c, row, odir)
                # diagnostics / reporting (never used for fitting or selection)
                if region == "head":
                    for lab in ("C1", "C2"):
                        _, vP, _, vF = get([LAB[lab]])
                        row[f"diag_{lab}_centroid_A2"] = centroid_dist(vP, vF)
                        row[f"diag_{lab}_centroid_A3"] = centroid_dist(vP, vF, R, c)
                    zm = np.median(fsP[:, 2])
                    for name, fpart, tpart in [("upper_to_lower", fsP[:, 2] >= zm, fsP[:, 2] < zm),
                                               ("lower_to_upper", fsP[:, 2] < zm, fsP[:, 2] >= zm)]:
                        fF = fsF[:, 2] >= zm if name == "upper_to_lower" else fsF[:, 2] < zm
                        tF = ~fF
                        if fpart.sum() > 100 and fF.sum() > 100 and tpart.sum() > 100 and tF.sum() > 100:
                            Rh, ch, _, _ = icp(subsample(fsF[fF], MAX_PTS, rng), subsample(fsP[fpart], MAX_PTS, rng))
                            row[f"diag_skullCV_{name}_surface_A2"] = sym_surface_distance(subsample(fsP[tpart], MAX_PTS, rng), subsample(fsF[tF], MAX_PTS, rng))
                            row[f"diag_skullCV_{name}_surface_half"] = sym_surface_distance(subsample(fsP[tpart], MAX_PTS, rng),
                                                                                            apply(Rh, ch, subsample(fsF[tF], MAX_PTS, rng)))
                            row[f"diag_skullCV_{name}_rot_diff_deg"] = rot_deg(Rh @ R.T)
                else:
                    worse = []
                    for side in ("L", "R"):
                        _, vP, _, vF = get([LAB[f"femur_{side}"]])
                        if len(vP) and len(vF):
                            hp = vP[vP[:, 2] >= vP[:, 2].max() - FEMHEAD_MM]; hf = vF[vF[:, 2] >= vF[:, 2].max() - FEMHEAD_MM]
                            a2, a3 = centroid_dist(hp, hf), centroid_dist(hp, hf, R, c)
                            row[f"rep_femhead_{side}_A2"], row[f"rep_femhead_{side}_A3"] = a2, a3
                            worse.append(a3 - a2)
                    row["rep_femhead_worst_change"] = max(worse) if worse else np.nan
                    for grp, labs in [("vessels", ["iliac_art_L", "iliac_art_R", "iliac_vein_L", "iliac_vein_R"]),
                                      ("muscles", ["glut_max_L", "glut_max_R", "glut_med_L", "glut_med_R", "glut_min_L",
                                                   "glut_min_R", "iliopsoas_L", "iliopsoas_R"])]:
                        sP, _, sF, _ = get([LAB[x] for x in labs])
                        sP, sF = subsample(sP, MAX_PTS, rng), subsample(sF, MAX_PTS, rng)
                        row[f"supp_{grp}_surface_A2"] = sym_surface_distance(sP, sF)
                        row[f"supp_{grp}_surface_A3"] = sym_surface_distance(sP, apply(R, c, sF))
            row["estimable"] = ok
            row["decision"] = decision_reason(row if o else {}, ok, why) if ok else f"NOT ESTIMABLE ({why}) -> Amendment 2"
            row["provisional_adopt"] = bool(ok and row.get("provisional_adopt", False))
            rows.append(row)
        print(f"[{i:3d}/{len(cases)}] {cid}  ({time.time() - t0:.0f} s)")
    F = pd.DataFrame(rows)
    if preflight or overlay_cases:
        pd.set_option("display.width", 250)
        cols = ["case_id", "region", "estimable", "n_fit_P", "n_fit_F", "converged", "iterations", "rot_deg", "trans_mm",
                "oneway_flag", "within_caps"]
        print(F[[c for c in cols if c in F]].round(2).to_string(index=False))
        print("OVERLAY CHECK: PNGs written to " + str(odir) + "; no tables written." if overlay_cases else
              "PREFLIGHT (implementation checks only): nothing written.")
        return
    # safety monitor (R7)
    summ = []
    for region in ("head", "pelvis"):
        r = F[F.region == region]
        ad = r[r.provisional_adopt]
        if region == "pelvis":
            ev = ad.dropna(subset=["rep_femhead_worst_change"]) if "rep_femhead_worst_change" in ad else ad.iloc[0:0]
            nworse = int((ev["rep_femhead_worst_change"] > SAFETY_MM).sum()) if len(ev) else 0
        else:
            cols = [c for c in ["diag_C1_centroid_A2", "diag_C2_centroid_A2"] if c in ad]
            ev = ad.dropna(subset=cols) if cols else ad.iloc[0:0]
            ch = pd.concat([ev.get("diag_C1_centroid_A3", np.nan) - ev.get("diag_C1_centroid_A2", np.nan),
                            ev.get("diag_C2_centroid_A3", np.nan) - ev.get("diag_C2_centroid_A2", np.nan)], axis=1).max(axis=1) if len(ev) else pd.Series(dtype=float)
            nworse = int((ch > SAFETY_MM).sum())
        frac = nworse / len(ev) if len(ev) else np.nan
        summ.append(dict(region=region, n_patients=len(r), n_estimable=int(r.estimable.sum()),
                         n_converged=int(r.get("converged", pd.Series(dtype=bool)).fillna(False).sum()),
                         n_within_caps=int(r.get("within_caps", pd.Series(dtype=bool)).fillna(False).sum()),
                         n_provisionally_adopted=len(ad), n_oneway_flag=int(r.get("oneway_flag", pd.Series(dtype=bool)).fillna(False).sum()),
                         safety_denominator=len(ev), safety_n_worse_gt_2mm=nworse, safety_fraction=frac,
                         safety_trigger=bool(np.isfinite(frac) and frac > SAFETY_FRAC)))
    S = pd.DataFrame(summ)
    out = []
    for name, df in [("fit", F), ("transforms", pd.DataFrame(trows)), ("summary", S)]:
        p = OUT_DIR / f"Amendment3_{name}_{run_id}.csv"
        df.to_csv(p, index=False)
        out.append(str(p))
    script = Path(__file__).resolve()
    man = dict(step="Amendment 3 numerical stage (provisional; visual veto pending)", run=run_id, script=str(script),
               script_MD5=md5(script), inputs={str(OUT_DIR / A2_PATIENT): md5(OUT_DIR / A2_PATIENT)}, outputs=out,
               parameters=dict(TRIM=TRIM, MAX_ITER=MAX_ITER, TOL_MM=TOL_MM, TOL_DEG=TOL_DEG, FOV_MARGIN=FOV_MARGIN,
                               MIN_POINTS=MIN_POINTS, MIN_SKULL_SPAN=MIN_SKULL_SPAN, CAPS=CAPS,
                               MIN_IMPROVE=MIN_IMPROVE, SURF_WORSE=SURF_WORSE, N_BOOT=N_BOOT, SEED=SEED, MAX_PTS=MAX_PTS,
                               FEMHEAD_MM=FEMHEAD_MM, SAFETY=(SAFETY_FRAC, SAFETY_MM)),
               runtime_s=round(time.time() - t0),
               software=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__,
                             SimpleITK=sitk.Version_VersionString(), platform=platform.platform()))
    (OUT_DIR / f"Amendment3_manifest_{run_id}.json").write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    rep = [f"Amendment 3 numerical stage - run {run_id} - {round(time.time() - t0)} s", S.round(3).to_string(index=False),
           "", F["decision"].groupby(F["region"]).value_counts().to_string()]
    (OUT_DIR / f"Amendment3_report_{run_id}.txt").write_text("\n".join(rep) + "\n", encoding="utf-8")
    print("\n".join(rep))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", type=int)
    ap.add_argument("--overlay", nargs="+", help="case IDs: diagnostic overlays of A2 vs proposed correction (no bootstrap, no tables)")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.preflight, a.overlay)


if __name__ == "__main__":
    main()
