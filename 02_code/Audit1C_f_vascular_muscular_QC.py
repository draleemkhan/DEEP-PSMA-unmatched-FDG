"""
Audit1C_f_vascular_muscular_QC.py
===========================================================================
PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - AUDIT 1C-f
Independent vascular and muscular registration QC (descriptive).
Implements protocol v2.9, Section 8.4 (specification; DRAFT until script and
synthetic-validation criteria are reviewed and the MD5 is locked).

Scope: CT-CT anatomical alignment of structures NOT used to estimate any
registration. Descriptive only. NO acceptance thresholds or categories, NO
tumour masks of any kind, NO lesion matching, correspondence, DMTVF or
discordance.

Transformations (no new optimisation)
  T_S  supplied (resampling PSMA/rigid.tfm), all patients
  T_A  adopted: translation-corrected R'(q) = PSMA/rigid.tfm(q + t) for 97 patients,
       supplied for the 3 Amendment 2 fallbacks (train_0018, train_0049, train_0067)
  T_C  translation-corrected for the 3 fallback patients only (descriptive, separate)
  t and fallback status are read from Audit 1C-e (per-patient table, run 20261031_141128).

Structures (Section 8.4.4; TotalSegmentator v2.4 task=total labels)
  vascular  : aorta 52, IVC 63, iliac vena L/R 67/68, iliac artery L/R 65/66
              -> PRESENT_COMPLETE on both CTs AND frac_largest_component >= 0.95 on both
  muscular  : iliopsoas L/R 88/89, gluteus medius L/R 82/83, gluteus minimus L/R 84/85,
              gluteus maximus L/R 80/81 -> PRESENT_COMPLETE on both CTs
  exploratory: prostate 22 (separate table; excluded from regional summaries)

Regions (Section 8.4.5; PSMA CT only): boundaries at midpoints between adjacent
vertebral centroids (C7-T1, T12-L1, L5-S1). Unreliable boundary vertebrae (not
PRESENT_COMPLETE on the PSMA CT, fragmented, adjacent to a fragmented vertebra, or
Audit 1B flagged levels as used in Amendment 1) are replaced by linear interpolation
of level position between the nearest reliable vertebrae above and below; every
substitution is logged. No extrapolation: a boundary without reliable vertebrae on
both sides is undefined and the adjacent regions are not evaluable (logged).

Metric (Section 8.4.6)
  - surface points = 6-neighbour contour of the FULL masks (no cropping), physical mm
  - region window W = z-interval common to both masks within the region
  - sampling window S = W shrunk by END_MARGIN_MM at each end; target window = W
  - BOTH directions: source surface points in S -> nearest target surface points of
    the other mask in W only (KD-tree, physical coordinates)
  - evaluable if len(S) >= MIN_WINDOW_MM and >= MIN_SURFACE_POINTS in each direction
  - symmetric MSD (pooled mean) = primary; pooled HD95, vessel slice-wise in-plane
    centroid displacement, muscle centroid distance (S-restricted) = supporting

Built-in synthetic validation (Section 8.4.7) runs FIRST; any failure stops the run.

Usage (Anaconda Prompt, env totalseg):
  python Audit1C_f_vascular_muscular_QC.py --synthetic-only   (synthetic validation only;
                                                               no patient data read)
  python Audit1C_f_vascular_muscular_QC.py --preflight 3      (first 3 patients + 1 fallback: eligibility,
                                                               boundaries, Audit 1C-e reproduction and evaluability
                                                               counts only; no residual values)
  python Audit1C_f_vascular_muscular_QC.py                    (definitive run; refused
                                                               until SPEC_APPROVED = True)
Requires: Python 3.9+, SimpleITK, numpy, pandas, scipy.
===========================================================================
"""
import argparse
import hashlib
import json
import platform
import sys
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import scipy
    from scipy import ndimage
    from scipy.spatial import cKDTree
    import SimpleITK as sitk
except ImportError as e:
    sys.exit(f"Missing package: {e}")

# ------------------------------ configuration (locked with the script MD5) ------------------------------
DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
AUDIT1B_LABELS = OUT_DIR / "Audit1B_label_level_20261029_134401.csv"
A1E_PATIENT = OUT_DIR / "Audit1C_e_per_patient_20261031_141128.csv"
A1E_STRUCT = OUT_DIR / "Audit1C_e_structure_level_20261031_141128.csv"

# The definitive run refuses to start until the specification and script are approved
# and the MD5 is recorded in protocol Section 8.4 (synthetic-only mode is always allowed).
SPEC_APPROVED = True                # approved 31 Oct 2026 after synthetic validation and pre-flight (protocol 8.4.10)

EXPECTED_PATIENTS = 100
VASCULAR = {52: "aorta", 63: "inferior_vena_cava", 67: "iliac_vena_left", 68: "iliac_vena_right",
            65: "iliac_artery_left", 66: "iliac_artery_right"}
MUSCULAR = {88: "iliopsoas_left", 89: "iliopsoas_right", 82: "gluteus_medius_left", 83: "gluteus_medius_right",
            84: "gluteus_minimus_left", 85: "gluteus_minimus_right", 80: "gluteus_maximus_left",
            81: "gluteus_maximus_right"}
EXPLORATORY = {22: "prostate"}
FRAG_MIN = 0.95                     # frac_largest_component >= 0.95 on both CTs (vessels only, 8.4.4)
END_MARGIN_MM = 10.0                # provisional (8.4.6)
MIN_WINDOW_MM = 20.0                # provisional; NOT a validated minimum vessel length
MIN_SURFACE_POINTS = 100            # per direction (8.4.6)
ROI_MARGIN_VOX = 3                  # union bounding box margin (as Amendment 1)
REPRO_TOL_MM = 0.01                 # aorta centroid distances must reproduce Audit 1C-e (hard stop)
COMPOSE_TOL_MM = 1e-6               # R'(q) == R(q + t) check (hard stop)
FALLBACKS_EXPECTED = ["train_0018", "train_0049", "train_0067"]

# vertebral level positions: S1 = 0, L5 = 1, ..., L1 = 5, T12 = 6, ..., T1 = 17, C7 = 18, ..., C1 = 24
VERT_LABELS = list(range(26, 51))                       # label = 26 + position
BOUNDARIES = {"C7T1": 17.5, "T12L1": 5.5, "L5S1": 0.5}
# Audit 1B flagged unreliable vertebrae (affected + adjacent), identical to Amendment 1 (FLAGGED)
FLAGGED = {"train_0088": [39, 40, 41, 42, 43], "train_0059": [26, 27, 28, 29, 30], "train_0029": [29, 30, 31, 32, 33]}

# synthetic validation (8.4.7)
SYN_GRIDS = {"anisotropic_0.98x0.98x3.0": (0.98, 0.98, 3.0), "isotropic_1.0": (1.0, 1.0, 1.0)}
SYN_SHIFTS_MM = [2.0, 5.0, 10.0]
SYN_ANGLES_DEG = [5.0, 10.0]
SYN_STABILITY_WINDOWS_MM = [10.0, 20.0, 40.0, 100.0]

# ------------------------------ mode, logging ------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--synthetic-only", action="store_true", help="run only the synthetic validation (no patient data)")
ap.add_argument("--preflight", type=int, default=0,
                help="implementation check on the first N patients: eligibility, boundaries, reproduction of Audit 1C-e "
                     "and evaluability counts only; NO residual (MSD/HD95/centroid) values are written or printed")
ARGS = ap.parse_args()
SYN_ONLY = ARGS.synthetic_only
PREFLIGHT = ARGS.preflight > 0
PREFIX = "Audit1C_f_SYNTHETIC_" if SYN_ONLY else ("Audit1C_f_PREFLIGHT_" if PREFLIGHT else "Audit1C_f_")
RUN = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG = []


def say(text=""):
    print(text, flush=True)
    LOG.append(str(text))


def out(name, ext):
    return OUT_DIR / f"{PREFIX}{name}_{RUN}.{ext}"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stop_run(reason):
    say(f"\n*** RUN STOPPED: {reason} ***")
    say("Implementation or input error; no result is reported. Document before any re-run.")
    out("RUN_STOPPED", "txt").write_text("\n".join(LOG), encoding="utf-8")
    sys.exit(1)


# ------------------------------ geometry and surfaces ------------------------------
class Geom:
    """Voxel index (x, y, z) -> physical (LPS, mm) for a SimpleITK image."""

    def __init__(self, img):
        self.o = np.array(img.GetOrigin(), float)
        self.s = np.array(img.GetSpacing(), float)
        self.D = np.array(img.GetDirection(), float).reshape(3, 3)
        self.size = np.array(img.GetSize(), int)

    def phys(self, idx_xyz):
        return self.o + (np.asarray(idx_xyz, float) * self.s) @ self.D.T


STRUCT6 = ndimage.generate_binary_structure(3, 1)       # 6-neighbourhood (fullyConnected = False)


def label_points(arr_zyx, lab, bbox_lo, bbox_hi, geom):
    """Full-mask voxels and 6-neighbour contour voxels of one label inside an ROI.
    Returns dict(vox=(N,3) phys, vox_k=(N,) z-index, surf=(M,3) phys, surf_k=(M,))."""
    z0, y0, x0 = bbox_lo[2], bbox_lo[1], bbox_lo[0]
    sub = arr_zyx[bbox_lo[2]:bbox_hi[2], bbox_lo[1]:bbox_hi[1], bbox_lo[0]:bbox_hi[0]] == lab
    cont = sub & ~ndimage.binary_erosion(sub, structure=STRUCT6, border_value=0)
    out_ = {}
    for key, m in (("vox", sub), ("surf", cont)):
        zz, yy, xx = np.nonzero(m)
        idx = np.stack([xx + x0, yy + y0, zz + z0], axis=1)
        out_[key] = geom.phys(idx) if len(idx) else np.zeros((0, 3))
        out_[key + "_k"] = zz + z0
    return out_


def nn_dist(src, tgt):
    if len(src) == 0 or len(tgt) == 0:
        return np.zeros(0)
    d, _ = cKDTree(tgt).query(src, k=1)
    return d


def measure(P, F, z_lo_region, z_hi_region, cls, margin=END_MARGIN_MM, min_len=MIN_WINDOW_MM,
            min_pts=MIN_SURFACE_POINTS, restrict_targets=True):
    """Section 8.4.6 measurement for one structure-segment. P, F = label_points() dicts of the PSMA
    mask and the resampled FDG mask (same grid). Region bounds in physical z (LPS, superior positive)."""
    r = dict(evaluable=False, reason="", W_lo_mm=np.nan, W_hi_mm=np.nan, S_len_mm=np.nan,
             n_src_P=0, n_src_F=0, n_tgt_P=0, n_tgt_F=0, MSD_mm=np.nan, HD95_mm=np.nan,
             inplane_dx_mm=np.nan, inplane_dy_mm=np.nan, inplane_mag_mm=np.nan, n_slices=0, centroid_mm=np.nan)
    if len(P["vox"]) == 0 or len(F["vox"]) == 0:
        r["reason"] = "label absent on one grid"
        return r
    W_lo = max(P["vox"][:, 2].min(), F["vox"][:, 2].min(), z_lo_region)
    W_hi = min(P["vox"][:, 2].max(), F["vox"][:, 2].max(), z_hi_region)
    r.update(W_lo_mm=W_lo, W_hi_mm=W_hi)
    if not W_hi > W_lo:
        r["reason"] = "no common extent in region"
        return r
    S_lo, S_hi = W_lo + margin, W_hi - margin
    r["S_len_mm"] = S_hi - S_lo
    if S_hi - S_lo < min_len:
        r["reason"] = f"sampling window < {min_len:g} mm"
        return r
    inS = lambda pts: (pts[:, 2] >= S_lo) & (pts[:, 2] <= S_hi)
    inW = lambda pts: (pts[:, 2] >= W_lo) & (pts[:, 2] <= W_hi)
    srcP, srcF = P["surf"][inS(P["surf"])], F["surf"][inS(F["surf"])]
    if restrict_targets:
        tgtP, tgtF = P["surf"][inW(P["surf"])], F["surf"][inW(F["surf"])]
    else:   # comparator used only in the synthetic region test
        tgtP, tgtF = P["surf"], F["surf"]
    r.update(n_src_P=len(srcP), n_src_F=len(srcF), n_tgt_P=len(tgtP), n_tgt_F=len(tgtF))
    if min(len(srcP), len(srcF)) < min_pts or min(len(tgtP), len(tgtF)) == 0:
        r["reason"] = f"< {min_pts} surface points in a direction"
        return r
    d = np.concatenate([nn_dist(srcP, tgtF), nn_dist(srcF, tgtP)])
    r.update(evaluable=True, MSD_mm=float(d.mean()), HD95_mm=float(np.percentile(d, 95)))
    vP, vF = inS(P["vox"]), inS(F["vox"])
    if cls == "vascular":
        kP, kF = P["vox_k"][vP], F["vox_k"][vF]
        xyP, xyF = P["vox"][vP][:, :2], F["vox"][vF][:, :2]
        common = np.intersect1d(np.unique(kP), np.unique(kF))
        if len(common):
            diffs = np.array([xyF[kF == k].mean(axis=0) - xyP[kP == k].mean(axis=0) for k in common])
            r.update(inplane_dx_mm=float(diffs[:, 0].mean()), inplane_dy_mm=float(diffs[:, 1].mean()),
                     inplane_mag_mm=float(np.linalg.norm(diffs, axis=1).mean()), n_slices=int(len(common)))
    else:
        if vP.any() and vF.any():
            r["centroid_mm"] = float(np.linalg.norm(F["vox"][vF].mean(axis=0) - P["vox"][vP].mean(axis=0)))
    return r


def naive_crop_msd(P, F, S_lo, S_hi):
    """Comparator ONLY for the synthetic cut-face check (never used for patients): both masks cropped
    to S, so their contours gain artificial cut faces (the end slices of the cropped masks); MSD over
    all resulting contour points."""
    surfs = []
    for X in (P, F):
        keep = (X["vox"][:, 2] >= S_lo) & (X["vox"][:, 2] <= S_hi)
        v, k = X["vox"][keep], X["vox_k"][keep]
        full_surf = set(map(tuple, np.round(X["surf"], 6)))
        lateral = np.array([tuple(np.round(q, 6)) in full_surf for q in v], dtype=bool)
        surfs.append(v[lateral | (k == k.min()) | (k == k.max())])
    d = np.concatenate([nn_dist(surfs[0], surfs[1]), nn_dist(surfs[1], surfs[0])])
    return float(d.mean())


# ------------------------------ synthetic validation ------------------------------
def syn_grid(spacing, half_extent=(80.0, 80.0, 130.0)):
    sp = np.array(spacing, float)
    size = [int(np.ceil(2 * h / s)) + 1 for h, s in zip(half_extent, sp)]
    img = sitk.Image(size, sitk.sitkUInt8)
    img.SetSpacing(tuple(sp))
    img.SetOrigin(tuple(-np.array(half_extent)))
    return img


def syn_shape(ref, kind):
    g = Geom(ref)
    nx, ny, nz = g.size
    zz, yy, xx = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    P = g.phys(np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1))
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    if kind == "cylinder":                                    # vessel-like, axis z, r = 8 mm
        m = (x ** 2 + y ** 2 <= 8.0 ** 2) & (np.abs(z) <= 100)
    elif kind == "tilted_cylinder":                            # axis tilted 30 deg in x-z
        a = np.array([np.sin(np.radians(30)), 0, np.cos(np.radians(30))])
        proj = P @ a
        perp = np.linalg.norm(P - np.outer(proj, a), axis=1)
        m = (perp <= 8.0) & (np.abs(proj) <= 100)
    elif kind == "ellipsoid":                                  # muscle-like
        m = (x / 30) ** 2 + (y / 20) ** 2 + (z / 60) ** 2 <= 1
    elif kind == "tapered_block":
        hw = 25 - (z + 60) * (15 / 120)
        m = (np.abs(z) <= 60) & (np.abs(x) <= hw) & (np.abs(y) <= 15)
    elif kind == "L_tube":                                     # region-crossing test (8.4.7)
        vert = (x ** 2 + y ** 2 <= 6.0 ** 2) & (z >= -8) & (z <= 100)
        horiz = ((y ** 2 + (z + 8) ** 2) <= 6.0 ** 2) & (x >= 0) & (x <= 60)
        m = vert | horiz
    else:
        raise ValueError(kind)
    arr = m.reshape(nz, ny, nx).astype(np.uint8)
    img = sitk.GetImageFromArray(arr)
    img.CopyInformation(ref)
    return img


def syn_move(img, transform):
    """Moved shape through the SAME path used for patients (sitk.Resample, nearest neighbour)."""
    return sitk.Resample(img, img, transform, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)


def syn_points(img):
    arr = sitk.GetArrayViewFromImage(img)
    g = Geom(img)
    return label_points(arr, 1, np.array([0, 0, 0]), g.size, g)


def expected_cyl_msd(r, d, n=200000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    p = np.stack([r * np.cos(th), r * np.sin(th)], axis=1)
    return float(np.abs(np.linalg.norm(p - np.array([d, 0.0]), axis=1) - r).mean())


def synthetic_validation():
    say("\n===== SYNTHETIC VALIDATION (Section 8.4.7; implementation check, not anatomical accuracy) =====")
    rows, fails = [], []
    inf = 1e9

    def check(name, ok, **kw):
        rows.append(dict(test=name, passed=bool(ok), **kw))
        if not ok:
            fails.append(name)
        say(f"  [{'PASS' if ok else 'FAIL'}] {name} " + " ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                                                         for k, v in kw.items()))

    tol_all = 0.5 * max(max(s) for s in SYN_GRIDS.values())
    grid_msd = {}
    for gname, sp in SYN_GRIDS.items():
        ref = syn_grid(sp)
        tol = 0.5 * max(sp)
        shapes = {k: syn_shape(ref, k) for k in ["cylinder", "tilted_cylinder", "ellipsoid", "tapered_block"]}
        pts = {k: syn_points(v) for k, v in shapes.items()}
        # T1 zero displacement
        for k, img in shapes.items():
            F = syn_points(syn_move(img, sitk.TranslationTransform(3, (0.0, 0.0, 0.0))))
            r = measure(pts[k], F, -inf, inf, "vascular" if "cylinder" in k else "muscular")
            check(f"{gname}|zero|{k}", r["evaluable"] and r["MSD_mm"] == 0.0, MSD=r["MSD_mm"])
        # T2 in-plane translations of the straight cylinder (x and y)
        for ax in (0, 1):
            msds = []
            for d in SYN_SHIFTS_MM:
                off = [0.0, 0.0, 0.0]; off[ax] = d
                F = syn_points(syn_move(shapes["cylinder"], sitk.TranslationTransform(3, tuple(off))))
                r = measure(pts["cylinder"], F, -inf, inf, "vascular")
                exp = expected_cyl_msd(8.0, d)
                msds.append(r["MSD_mm"])
                grid_msd[(gname, ax, d)] = r["MSD_mm"]
                check(f"{gname}|inplane_axis{ax}|{d:g}mm|centroid", abs(r["inplane_mag_mm"] - d) <= tol,
                      measured=r["inplane_mag_mm"], expected=d, tol=tol)
                check(f"{gname}|inplane_axis{ax}|{d:g}mm|MSD_vs_expected", abs(r["MSD_mm"] - exp) <= tol,
                      MSD=r["MSD_mm"], expected=exp, tol=tol)
                if d == 5.0 and ax == 0:   # T5 cut-face check
                    S_lo, S_hi = r["W_lo_mm"] + END_MARGIN_MM, r["W_hi_mm"] - END_MARGIN_MM
                    naive = naive_crop_msd(pts["cylinder"], F, S_lo, S_hi)
                    check(f"{gname}|cut_face|5mm", r["MSD_mm"] >= exp - tol and abs(r["MSD_mm"] - exp) <= tol,
                          windowed=r["MSD_mm"], naive_crop=naive, expected=exp, naive_bias=naive - exp)
            check(f"{gname}|inplane_axis{ax}|monotonic", all(np.diff(msds) > 0), MSDs=str([round(m, 3) for m in msds]))
        # T3 z-translation: cylinder insensitive (expected), ellipsoid detects
        for k, expect_detect in (("cylinder", False), ("ellipsoid", True)):
            msds = []
            for d in (5.0, 10.0):
                F = syn_points(syn_move(shapes[k], sitk.TranslationTransform(3, (0.0, 0.0, d))))
                r = measure(pts[k], F, -inf, inf, "vascular" if k == "cylinder" else "muscular")
                msds.append(r["MSD_mm"])
            if expect_detect:
                check(f"{gname}|z_shift|{k}|detected", msds[1] > msds[0] > 0, MSDs=str([round(m, 3) for m in msds]))
            else:
                check(f"{gname}|z_shift|{k}|axial_insensitivity", max(msds) <= tol, MSDs=str([round(m, 3) for m in msds]), tol=tol)
        # T4 rotations about x, y, z (shape centre = origin)
        for k in ("tilted_cylinder", "ellipsoid", "tapered_block"):
            for axis in range(3):
                msds = []
                for a in SYN_ANGLES_DEG:
                    T = sitk.Euler3DTransform()
                    ang = [0.0, 0.0, 0.0]; ang[axis] = np.radians(a)
                    T.SetRotation(*ang)
                    F = syn_points(syn_move(shapes[k], T))
                    r = measure(pts[k], F, -inf, inf, "muscular")
                    msds.append(r["MSD_mm"])
                check(f"{gname}|rotation|{k}|axis{axis}", msds[1] > msds[0] > 0, MSDs=str([round(m, 3) for m in msds]))
        # T6 region restriction (L-shaped tube crossing a boundary at z = 0)
        Limg = syn_shape(ref, "L_tube"); PL = syn_points(Limg)
        FL = syn_points(syn_move(Limg, sitk.TranslationTransform(3, (0.0, 5.0, 0.0))))
        rr = measure(PL, FL, 0.0, inf, "vascular", restrict_targets=True)
        ru = measure(PL, FL, 0.0, inf, "vascular", restrict_targets=False)
        exp = expected_cyl_msd(6.0, 5.0)
        check(f"{gname}|region_restriction|5mm", rr["MSD_mm"] >= ru["MSD_mm"] - 1e-9 and abs(rr["MSD_mm"] - exp) <= tol,
              restricted=rr["MSD_mm"], unrestricted=ru["MSD_mm"], expected=exp, tol=tol)
        # stress case: 20 mm shift along the limb axis; lateral partner distances (up to ~20 mm) exceed the
        # distance to the out-of-region limb just below the boundary (~11-15 mm on these grids), so an
        # unrestricted search returns falsely small distances; the restricted measurement must not
        FL = syn_points(syn_move(Limg, sitk.TranslationTransform(3, (20.0, 0.0, 0.0))))
        rr = measure(PL, FL, 0.0, inf, "vascular", restrict_targets=True)
        ru = measure(PL, FL, 0.0, inf, "vascular", restrict_targets=False)
        exp = expected_cyl_msd(6.0, 20.0)
        check(f"{gname}|region_restriction|20mm_stress", rr["MSD_mm"] > ru["MSD_mm"] and abs(rr["MSD_mm"] - exp) <= tol,
              restricted=rr["MSD_mm"], unrestricted=ru["MSD_mm"], expected=exp, tol=tol)
        # T7 stability (report only)
        F = syn_points(syn_move(shapes["cylinder"], sitk.TranslationTransform(3, (5.0, 0.0, 0.0))))
        exp = expected_cyl_msd(8.0, 5.0)
        for w in SYN_STABILITY_WINDOWS_MM:
            half = w / 2 + END_MARGIN_MM
            r = measure(pts["cylinder"], F, -half, half, "vascular", min_len=0.0, min_pts=1)
            rows.append(dict(test=f"{gname}|stability|window{w:g}mm", passed=True, report_only=True,
                             MSD=r["MSD_mm"], expected=exp, abs_error=abs(r["MSD_mm"] - exp),
                             n_src_P=r["n_src_P"], n_src_F=r["n_src_F"]))
            say(f"  [REPORT] {gname}|stability|window {w:g} mm: MSD {r['MSD_mm']:.3f} (expected {exp:.3f}), "
                f"surface points {r['n_src_P']}/{r['n_src_F']}")
    # grid agreement (anisotropic vs isotropic)
    for ax in (0, 1):
        for d in SYN_SHIFTS_MM:
            a = grid_msd[("anisotropic_0.98x0.98x3.0", ax, d)]; b = grid_msd[("isotropic_1.0", ax, d)]
            check(f"grid_agreement|axis{ax}|{d:g}mm", abs(a - b) <= tol_all, aniso=a, iso=b, tol=tol_all)
    df = pd.DataFrame(rows)
    return df, fails


# ------------------------------ start ------------------------------
script_path = Path(__file__).resolve()
SCRIPT_MD5 = md5(script_path)
say(f"Audit 1C-f (independent vascular and muscular registration QC) - run {RUN}")
say(f"Script: {script_path} (MD5 {SCRIPT_MD5})")
say("Protocol v2.9 Section 8.4. Descriptive only: no thresholds, no tumour data, no matching or discordance.")
if not SYN_ONLY and not PREFLIGHT and not SPEC_APPROVED:
    print("DEFINITIVE RUN REFUSED: SPEC_APPROVED is False. Review the specification, script and synthetic validation, "
          "set the flag, record the MD5 in protocol Section 8.4, then run. Nothing was processed.")
    sys.exit(2)
OUT_DIR.mkdir(parents=True, exist_ok=True)

SYN, syn_fails = synthetic_validation()
SYN.to_csv(out("synthetic_validation", "csv"), index=False)
if syn_fails:
    stop_run(f"synthetic validation failed: {syn_fails}")
say("Synthetic validation: ALL PASSED")
if SYN_ONLY:
    manifest = dict(step="Audit 1C-f synthetic validation only", run=RUN, script=str(script_path), script_MD5=SCRIPT_MD5,
                    passed=True, n_tests=int((SYN.get("report_only") != True).sum()),  # noqa: E712
                    software=dict(python=sys.version, SimpleITK=sitk.Version_VersionString(), numpy=np.__version__,
                                  pandas=pd.__version__, scipy=scipy.__version__, platform=platform.platform()))
    out("manifest", "json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    out("report", "txt").write_text("\n".join(LOG), encoding="utf-8")
    say("SYNTHETIC-ONLY run finished. No patient data were read.")
    sys.exit(0)

# ------------------------------ inputs ------------------------------
for p in (AUDIT1B_LABELS, A1E_PATIENT, A1E_STRUCT):
    if not p.is_file():
        stop_run(f"input not found: {p}")
L1B = pd.read_csv(AUDIT1B_LABELS)
L1B = L1B.set_index(["case_id", "tracer", "label_id"])
PE = pd.read_csv(A1E_PATIENT).set_index("case_id")
SE = pd.read_csv(A1E_STRUCT)
SE_AORTA = SE[SE["label_id"] == 52].set_index("case_id")
fallbacks = sorted(PE.index[PE["would_fallback_if_adopted"].astype(str).str.lower() == "true"])
if fallbacks != FALLBACKS_EXPECTED:
    stop_run(f"fallback patients {fallbacks} differ from the adoption decision {FALLBACKS_EXPECTED}")
if not (PE["status"] == "CORRECTION_COMPUTED").all():
    stop_run("Audit 1C-e per-patient table contains patients without a computed correction")
cases = sorted(p.name for p in DATA_ROOT.glob("train_*") if p.is_dir())
n_found = len(cases)
say(f"Case folders found: {n_found} (expected {EXPECTED_PATIENTS}); fallbacks {fallbacks}")
if PREFLIGHT:
    cases = cases[:ARGS.preflight] + [c for c in fallbacks if c not in cases[:ARGS.preflight]][:1]   # include one fallback (T_C path)
    say(f"*** PRE-FLIGHT: {cases}; implementation checks only; no residual values are written or printed ***")


def cat(cid, tracer, lab):
    try:
        r = L1B.loc[(cid, tracer, lab)]
        return str(r["category"]), float(r["frac_largest_component"])
    except KeyError:
        return "MISSING_IN_AUDIT1B", np.nan


def eligible(cid, lab, cls):
    cP, fP = cat(cid, "PSMA", lab); cF, fF = cat(cid, "FDG", lab)
    if cP != "PRESENT_COMPLETE" or cF != "PRESENT_COMPLETE":
        return False, f"not PRESENT_COMPLETE on both CTs (PSMA {cP}, FDG {cF})"
    if cls == "vascular" and not (fP >= FRAG_MIN and fF >= FRAG_MIN):
        return False, f"fragmentation > 5% (frac_largest PSMA {fP:.3f}, FDG {fF:.3f})"
    return True, ""


def boundaries(cid, cenP):
    """Region boundaries in physical z (PSMA CT) with interpolation for unreliable vertebrae."""
    bad = set(FLAGGED.get(cid, []))
    for lab in VERT_LABELS:
        c, f = cat(cid, "PSMA", lab)
        if f == f and f < FRAG_MIN:
            bad.update([lab - 1, lab, lab + 1])
    reliable = {}
    for lab in VERT_LABELS:
        c, f = cat(cid, "PSMA", lab)
        if c == "PRESENT_COMPLETE" and lab not in bad and lab in cenP:
            reliable[lab - 26] = cenP[lab][2]
    out_b, log = {}, []
    for name, pos in BOUNDARIES.items():
        below = [p for p in reliable if p < pos]; above = [p for p in reliable if p > pos]
        lo_direct, hi_direct = int(np.floor(pos)), int(np.ceil(pos))
        if not below or not above:
            out_b[name] = np.nan
            log.append(dict(case_id=cid, boundary=name, status="UNDEFINED",
                            detail="no reliable vertebra on one side (no extrapolation)"))
            continue
        pa, pb = max(below), min(above)
        za, zb = reliable[pa], reliable[pb]
        out_b[name] = za + (pos - pa) * (zb - za) / (pb - pa)
        sub = not (pa == lo_direct and pb == hi_direct)
        log.append(dict(case_id=cid, boundary=name, status="SUBSTITUTED" if sub else "DIRECT",
                        detail=f"vertebrae {26 + pa}/{26 + pb} (positions {pa}/{pb})", z_mm=out_b[name]))
    return out_b, log


def regions_from(b):
    """thorax: (T12L1, C7T1]; abdomen/retroperitoneum: (L5S1, T12L1]; pelvis: (-inf, L5S1].
    Regions needing an undefined boundary are skipped (logged). Defined boundaries must be ordered
    C7T1 > T12L1 > L5S1 (superior = larger physical z, LPS)."""
    inf = 1e9
    order = [b["C7T1"], b["T12L1"], b["L5S1"]]
    defined = [z for z in order if np.isfinite(z)]
    order_ok = all(defined[i] > defined[i + 1] for i in range(len(defined) - 1))
    regs = {}
    if np.isfinite(b["T12L1"]) and np.isfinite(b["C7T1"]):
        regs["thorax"] = (b["T12L1"], b["C7T1"])
    if np.isfinite(b["L5S1"]) and np.isfinite(b["T12L1"]):
        regs["abdomen_retroperitoneum"] = (b["L5S1"], b["T12L1"])
    if np.isfinite(b["L5S1"]):
        regs["pelvis"] = (-inf, b["L5S1"])
    return regs, order_ok


def tp(T, p):
    return np.array(T.TransformPoint(tuple(float(v) for v in p)))


def bbox_union(fP, fX, lab, size_xyz):
    b1, b2 = fP.GetBoundingBox(lab), fX.GetBoundingBox(lab)
    lo = np.minimum(np.array(b1[:3]), np.array(b2[:3])) - ROI_MARGIN_VOX
    hi = np.maximum(np.array(b1[:3]) + np.array(b1[3:]), np.array(b2[:3]) + np.array(b2[3:])) + ROI_MARGIN_VOX
    return np.maximum(lo, 0), np.minimum(hi, size_xyz)


# ------------------------------ main loop ------------------------------
seg_rows, elig_rows, bnd_rows, repro_rows, inventory, status_rows = [], [], [], [], [], []
LABELS = [(l, n, "vascular") for l, n in VASCULAR.items()] + [(l, n, "muscular") for l, n in MUSCULAR.items()] + \
         [(l, n, "exploratory") for l, n in EXPLORATORY.items()]

for ci, cid in enumerate(cases, 1):
    say(f"[{ci:3d}/{len(cases)}] {cid}")
    try:
        dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
        files = {"PSMA_totseg": dP / "totseg_24.nii.gz", "FDG_totseg": dF / "totseg_24.nii.gz",
                 "PSMA_tfm": dP / "rigid.tfm", "FDG_tfm": dF / "rigid.tfm"}
        for role, p in files.items():
            if not p.is_file():
                raise FileNotFoundError(f"{role}: {p}")
            st = p.stat()
            inventory.append(dict(case_id=cid, role=role, path=str(p), bytes=st.st_size,
                                  modified=datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"), md5=md5(p)))
        R = sitk.ReadTransform(str(files["PSMA_tfm"]))
        T = sitk.ReadTransform(str(files["FDG_tfm"]))
        segP = sitk.ReadImage(str(files["PSMA_totseg"]), sitk.sitkUInt8)
        segF = sitk.ReadImage(str(files["FDG_totseg"]), sitk.sitkUInt8)
        gP = Geom(segP)
        shpP = sitk.LabelShapeStatisticsImageFilter(); shpP.Execute(segP)
        shpF = sitk.LabelShapeStatisticsImageFilter(); shpF.Execute(segF)
        cenP = {int(l): np.array(shpP.GetCentroid(l)) for l in shpP.GetLabels()}
        cenF = {int(l): np.array(shpF.GetCentroid(l)) for l in shpF.GetLabels()}

        # adopted transformation (Audit 1C-e) and implementation checks (hard stops)
        t = PE.loc[cid, ["t_x_LPS", "t_y_LPS", "t_z_LPS"]].to_numpy(float)
        is_fb = cid in fallbacks
        Rc = sitk.CompositeTransform(3)
        Rc.AddTransform(R)
        Rc.AddTransform(sitk.TranslationTransform(3, tuple(float(v) for v in t)))     # applied first: R(q + t)
        test_q = [gP.phys(np.array(c)) for c in ([0, 0, 0], gP.size - 1, (gP.size - 1) / 2)]
        e_comp = max(np.linalg.norm(tp(Rc, q) - tp(R, q + t)) for q in test_q)
        if e_comp > COMPOSE_TOL_MM:
            stop_run(f"{cid}: corrected transform composition error {e_comp:.2e} mm")
        if 52 in cenP and 52 in cenF and cid in SE_AORTA.index:
            cs = float(np.linalg.norm(tp(T, cenF[52]) - cenP[52]))
            cc = float(np.linalg.norm(tp(T, cenF[52]) - t - cenP[52]))
            es, ec = abs(cs - SE_AORTA.loc[cid, "supplied_centroid_mm"]), abs(cc - SE_AORTA.loc[cid, "corrected_centroid_mm"])
            repro_rows.append(dict(case_id=cid, aorta_supplied_mm=cs, aorta_corrected_mm=cc,
                                   diff_supplied_vs_1Ce=es, diff_corrected_vs_1Ce=ec, composition_err_mm=e_comp))
            if es > REPRO_TOL_MM or ec > REPRO_TOL_MM:
                stop_run(f"{cid}: aorta centroid distances do not reproduce Audit 1C-e ({es:.4f}, {ec:.4f} mm)")
        else:
            repro_rows.append(dict(case_id=cid, note="aorta not available for reproduction check", composition_err_mm=e_comp))

        # region boundaries (PSMA CT)
        b, blog = boundaries(cid, cenP)
        bnd_rows += blog
        regs, order_ok = regions_from(b)
        if not order_ok:
            bnd_rows.append(dict(case_id=cid, boundary="ALL", status="ORDER_INVALID", detail=str(b)))
            regs = {}

        # transformations to evaluate
        trans = {"T_S": R, "T_A": (R if is_fb else Rc)}
        if is_fb:
            trans["T_C"] = Rc
        resampled, shapes = {}, {}
        for name, tr in trans.items():
            if name == "T_A" and is_fb:
                resampled[name], shapes[name] = resampled["T_S"], shapes["T_S"]
                continue
            img = sitk.Resample(segF, segP, tr, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
            f = sitk.LabelShapeStatisticsImageFilter(); f.Execute(img)
            resampled[name], shapes[name] = img, f
        arrP = sitk.GetArrayViewFromImage(segP)

        for lab, sname, cls in LABELS:
            ok, why = eligible(cid, lab, "vascular" if cls == "vascular" else "muscular")
            elig_rows.append(dict(case_id=cid, label_id=lab, structure=sname, cls=cls, eligible=ok, reason=why))
            if not ok:
                continue
            for name in trans:
                fX = shapes[name]
                if not (shpP.HasLabel(lab) and fX.HasLabel(lab)):
                    for rg in regs:
                        seg_rows.append(dict(case_id=cid, transform=name, label_id=lab, structure=sname, cls=cls, region=rg,
                                             evaluable=False, reason="label absent after resampling"))
                    continue
                lo, hi = bbox_union(shpP, fX, lab, gP.size)
                P = label_points(arrP, lab, lo, hi, gP)
                Fp = label_points(sitk.GetArrayViewFromImage(resampled[name]), lab, lo, hi, gP)
                for rg, (zlo, zhi) in regs.items():
                    r = measure(P, Fp, zlo, zhi, "vascular" if cls == "vascular" else "muscular")
                    seg_rows.append(dict(case_id=cid, transform=name, label_id=lab, structure=sname, cls=cls,
                                         region=rg, fallback_patient=is_fb, **r))
        status_rows.append(dict(case_id=cid, status="OK", regions=";".join(regs)))
    except SystemExit:
        raise
    except Exception as e:
        status_rows.append(dict(case_id=cid, status="FAILED", regions="", message=str(e)))
        say(f"   FAILED: {e}\n{traceback.format_exc(limit=2)}")

# ------------------------------ pre-flight output (no residual values) ------------------------------
if PREFLIGHT:
    SGp = pd.DataFrame(seg_rows)
    keep = ["case_id", "transform", "label_id", "structure", "cls", "region", "evaluable", "reason",
            "S_len_mm", "n_src_P", "n_src_F", "n_tgt_P", "n_tgt_F"]
    SGp = SGp[[c for c in keep if c in SGp.columns]]
    REp = pd.DataFrame(repro_rows)
    say("\n===== PRE-FLIGHT SUMMARY (implementation checks only) =====")
    say(pd.DataFrame(status_rows).to_string(index=False))
    if len(REp):
        say("Reproduction of Audit 1C-e (aorta centroid distances, hard stop 0.01 mm) and composition error:\n"
            + REp[[c for c in ["case_id", "diff_supplied_vs_1Ce", "diff_corrected_vs_1Ce", "composition_err_mm", "note"]
                   if c in REp.columns]].to_string(index=False))
    say("Boundaries:\n" + pd.DataFrame(bnd_rows).to_string(index=False))
    E = pd.DataFrame(elig_rows)
    say(f"Eligibility: {int(E['eligible'].sum())}/{len(E)} label-patient pairs eligible; excluded: "
        + "; ".join(f"{r.case_id} {r.structure}: {r.reason}" for r in E[~E['eligible']].itertuples()))
    if len(SGp):
        say("Evaluable segments by transform and region:\n"
            + SGp.groupby(["transform", "region"])["evaluable"].agg(["sum", "count"]).to_string())
        nev = SGp[~SGp["evaluable"]]
        if len(nev):
            say("Not evaluable (reasons): " + nev["reason"].value_counts().to_string())
    SGp.to_csv(out("segments_no_values", "csv"), index=False)
    E.to_csv(out("eligibility", "csv"), index=False)
    pd.DataFrame(bnd_rows).to_csv(out("boundaries_substitution_log", "csv"), index=False)
    REp.to_csv(out("reproduction_check", "csv"), index=False)
    out("report", "txt").write_text("\n".join(LOG), encoding="utf-8")
    say("PRE-FLIGHT finished. No residual values were computed for output.")
    sys.exit(0)

# ------------------------------ summaries (descriptive) ------------------------------
SG = pd.DataFrame(seg_rows)
ST = pd.DataFrame(status_rows)
main = SG[SG["transform"].isin(["T_S", "T_A"]) & (SG["cls"] != "exploratory")]
ev = main[main["evaluable"] == True]  # noqa: E712


def second_largest(s):
    s = np.sort(s.dropna().to_numpy())
    return float(s[-2]) if len(s) >= 2 else np.nan


PR = (ev.groupby(["case_id", "transform", "region", "cls"])["MSD_mm"]
      .agg(n_evaluable="count", median="median", P90=lambda s: s.quantile(.9), max="max",
           second_largest=second_largest).reset_index())
wide = PR.pivot_table(index=["case_id", "region", "cls"], columns="transform", values="median").reset_index()
if {"T_S", "T_A"}.issubset(wide.columns):
    wide["diff_T_A_minus_T_S"] = wide["T_A"] - wide["T_S"]
cohort = []
for (rg, cls), g in PR.groupby(["region", "cls"]):
    for tr, gg in g.groupby("transform"):
        s = gg["median"]
        cohort.append(dict(region=rg, cls=cls, transform=tr, n_patients=len(gg), median=s.median(),
                           Q1=s.quantile(.25), Q3=s.quantile(.75), P90=s.quantile(.9), max=s.max(),
                           median_n_evaluable=gg["n_evaluable"].median()))
    w = wide[(wide["region"] == rg) & (wide["cls"] == cls)]
    if "diff_T_A_minus_T_S" in w:
        dd = w["diff_T_A_minus_T_S"].dropna()
        cohort.append(dict(region=rg, cls=cls, transform="diff_T_A_minus_T_S", n_patients=len(dd), median=dd.median(),
                           Q1=dd.quantile(.25), Q3=dd.quantile(.75), P90=dd.quantile(.9), max=dd.max(),
                           n_worse_gt_1mm=int((dd > 1).sum()), n_worse_gt_2mm=int((dd > 2).sum()),
                           n_improved_gt_1mm=int((dd < -1).sum())))
CO = pd.DataFrame(cohort)
TC = SG[SG["transform"] == "T_C"]
EXP = SG[SG["cls"] == "exploratory"]

say("\n===== SUMMARY (descriptive; median of patient-level median MSD, mm) =====")
if len(CO):
    say(CO.round(2).to_string(index=False))
n_ok = int((ST["status"] == "OK").sum())
say(f"Patients processed OK: {n_ok}/{len(cases)}; segments evaluable (T_S/T_A, main set): {len(ev)}/{len(main)}")
say(f"Boundary substitutions: {sum(1 for r in bnd_rows if r.get('status') == 'SUBSTITUTED')}; "
    f"undefined: {sum(1 for r in bnd_rows if r.get('status') == 'UNDEFINED')}")

# ------------------------------ outputs ------------------------------
errors, written = [], []
for name, df in {"segments": SG, "patient_region_class": PR, "paired_differences": wide, "cohort_summary": CO,
                 "fallback_corrected_T_C": TC, "exploratory_prostate": EXP, "eligibility": pd.DataFrame(elig_rows),
                 "boundaries_substitution_log": pd.DataFrame(bnd_rows), "reproduction_check": pd.DataFrame(repro_rows),
                 "processing_status": ST, "synthetic_validation": SYN, "input_inventory": pd.DataFrame(inventory)}.items():
    try:
        f = out(name, "csv"); df.to_csv(f, index=False); written.append(str(f))
    except Exception as e:
        errors.append(f"{name}: {e}")
complete = (n_found == len(cases) == EXPECTED_PATIENTS and n_ok == EXPECTED_PATIENTS and not errors and not syn_fails)
completion = "COMPLETE" if complete else "INCOMPLETE"
manifest = dict(step="Audit 1C-f - independent vascular and muscular registration QC (protocol v2.9 Section 8.4)",
                run=RUN, completion_status=completion, script=str(script_path), script_MD5=SCRIPT_MD5,
                inputs={str(p): md5(p) for p in (AUDIT1B_LABELS, A1E_PATIENT, A1E_STRUCT)},
                parameters=dict(VASCULAR=VASCULAR, MUSCULAR=MUSCULAR, EXPLORATORY=EXPLORATORY, FRAG_MIN=FRAG_MIN,
                                END_MARGIN_MM=END_MARGIN_MM, MIN_WINDOW_MM=MIN_WINDOW_MM,
                                MIN_SURFACE_POINTS=MIN_SURFACE_POINTS, ROI_MARGIN_VOX=ROI_MARGIN_VOX,
                                REPRO_TOL_MM=REPRO_TOL_MM, BOUNDARIES=BOUNDARIES, FLAGGED=FLAGGED,
                                fallbacks=fallbacks),
                synthetic_validation_passed=not syn_fails, output_errors=errors, outputs=written,
                n_input_files=len(inventory),
                software=dict(python=sys.version, SimpleITK=sitk.Version_VersionString(), numpy=np.__version__,
                              pandas=pd.__version__, scipy=scipy.__version__, platform=platform.platform()),
                note="Descriptive only. No thresholds, tumour data, matching, correspondence or discordance.")
try:
    fM = out("manifest", "json"); tmp = Path(str(fM) + ".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    if json.loads(tmp.read_text(encoding="utf-8")).get("completion_status") != completion:
        raise RuntimeError("manifest read-back mismatch")
    tmp.replace(fM)
except Exception as e:
    completion = "INCOMPLETE"; say(f"MANIFEST ERROR: {e}")
say(f"\nAUDIT 1C-f RUN STATUS: {completion}")
out("report", "txt").write_text("\n".join(LOG), encoding="utf-8")
out(f"RUN_{completion}", "txt").write_text(f"Audit 1C-f run {RUN}: {completion}\n", encoding="utf-8")
