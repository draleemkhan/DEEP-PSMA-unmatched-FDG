"""
PU_psma_uptake_v1_0.py - PAPER 3, SAP section 8.2: PSMA uptake at FDG-lesion locations (exploratory, PROMISE-inspired).

Implements: locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d) section 8.2;
  00_protocol\\StepSAP_PSMAuptake_implementation_clarification_20261112.md (MD5 5d8deb6bce244aa480ce2aea0866a727);
  00_protocol\\StepSAP_PSMAuptake_implementation_clarification_addendum_20261112.md (MD5 88037af0077a32b5669d34f2402fcdd8).
Primary analysis, perturbation results and all locked classifications are read only, never modified.

Units: the 217 components with status V_U in the primary P3 run (P3_conn26_20261112_001915).
Per component:
  SUVmax  = max of the NATIVE PSMA PET (never resampled) over the component's unperturbed final mapped voxel set
            (P2b_conn26_20261112_000903, role 0 = final transform; linear_index_1based - 1 = C-order index (k, j, i)).
  extent  = maximum Euclidean distance (mm) between the voxel centres of the native FDG component
            (P1 label image, 26-connectivity, FDG PET grid, physical LPS coordinates). Not the P1 column max_extent_mm.
Per patient (all 100):
  references from the PSMA CT TotalSegmentator labels (totseg_24.nii.gz; aorta 52 = blood pool, liver 5):
  label -> PSMA PET grid by physical coordinates (nearest neighbour, identity transform) -> 1 x binary erosion with a
  3x3x3 element (border = background) -> remove PSMA TTB voxels -> SUVmean of native PSMA PET. Empty / absent ->
  reference unavailable (no substitute).
Classes (addendum section 3 precedence):
  1 extent <= 10 mm                     -> NC_partial_volume
  2 either reference unavailable        -> NC_reference_unavailable
  3 blood pool < liver: absent (<= bp), low (> bp and <= liver), positive (> liver)
    blood pool >= liver (flagged): absent (<= bp), positive (> bp), low unavailable
Exploratory count: SUVmax >= 3 (DEEP-PSMA segmentation threshold, not a clinical criterion).
Mixed-ligand cohort; ligand not inferred.
Outputs: 06_output\\PU_[PREFLIGHT_]<run>\\ PU_components_<run>.csv, PU_patients_<run>.csv, PU_summary_<run>.txt,
  PU_manifest_<run>.json.
Checks (stop on failure): MD5 of all inputs; per-case voxel-set CSV MD5 = P2b grid table; P1 label image MD5 =
  P1 patient table; PSMA PET grid = P2b grid record; PSMA TTB on PSMA PET grid; every V_U component has >= 1 final
  mapped voxel, all indices inside the PET array; native voxel count = P1 table.
Usage (folder 02_code, env totalseg):
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u PU_psma_uptake_v1_0.py --selftest
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u PU_psma_uptake_v1_0.py --preflight
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u PU_psma_uptake_v1_0.py
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
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
P1_RUN, P2B_RUN, P3_RUN = "P1_20261111_233936", "P2b_conn26_20261112_000903", "P3_conn26_20261112_001915"
CONN = 26
INPUTS = {
    (OUT_DIR, f"{P1_RUN}/P1_components_20261111_233936.csv"): "ffafbc7917e1b595eadacafab247fe88",
    (OUT_DIR, f"{P1_RUN}/P1_patients_20261111_233936.csv"): "f9b5b5908e523d1122547ef07d5c0afd",
    (OUT_DIR, f"{P1_RUN}/P1_manifest_20261111_233936.json"): "a27e4e4c5f44402c53487f9f8dd7dd8b",
    (OUT_DIR, f"{P2B_RUN}/P2b_manifest_20261112_000903.json"): "4320cad089bca8eecb128b9df5efd972",
    (OUT_DIR, f"{P2B_RUN}/P2b_psma_pet_grid_20261112_000903.csv"): "ad5504fe4fb4118b727c167b53019983",
    (OUT_DIR, f"{P3_RUN}/P3_components_20261112_001915.csv"): "ac0048cfe059c714b3c293fbe7e88279",
    (PROTO_DIR, "StepSAP_consolidated_LOCKED_20261111.md"): "0c64711db1207e31a70cf3f558ccbd7d",
    (PROTO_DIR, "StepSAP_PSMAuptake_implementation_clarification_20261112.md"): "5d8deb6bce244aa480ce2aea0866a727",
    (PROTO_DIR, "StepSAP_PSMAuptake_implementation_clarification_addendum_20261112.md"): "88037af0077a32b5669d34f2402fcdd8",
}
LABEL_AORTA, LABEL_LIVER = 52, 5
EXTENT_MM = 10.0
SUV_COUNT = 3.0
BRUTE_MAX = 3000
CLASSES = ["absent", "low", "positive", "NC_partial_volume", "NC_reference_unavailable"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ pure helpers (tested) ------------------------------
def phys_lps(origin, spacing, direction, kji):
    """physical LPS centres (mm) of voxels given as (k, j, i) rows (verbatim P1 convention)."""
    D = np.asarray(direction, float).reshape(3, 3); S = np.asarray(spacing, float)
    return np.asarray(origin, float) + (np.asarray(kji)[:, ::-1].astype(float) * S) @ D.T


def _brute_max(p):
    best = 0.0
    for a in range(0, len(p), 2000):
        blk = p[a:a + 2000]
        d2 = ((blk[:, None, :] - p[None, :, :]) ** 2).sum(-1)
        best = max(best, float(d2.max()))
    return float(np.sqrt(best))


def max_extent(pts):
    """maximum Euclidean distance between points (mm); exact. 0 for one point."""
    p = np.unique(np.asarray(pts, float), axis=0)
    if len(p) <= 1:
        return 0.0
    if len(p) <= BRUTE_MAX:
        return float(pdist(p).max())
    try:
        v = p[ConvexHull(p).vertices]               # the diameter is attained at hull vertices
    except Exception:                               # degenerate (collinear / coplanar): exact chunked brute force
        return _brute_max(p)
    return float(pdist(v).max()) if len(v) > 1 else 0.0


def erode1(mask):
    """1 x binary erosion, 3x3x3 element (26-neighbour), outside = background."""
    return ndimage.binary_erosion(mask, structure=np.ones((3, 3, 3), bool), iterations=1, border_value=0)


def reference(pet, label_on_pet, label_id, ttb):
    """returns dict(n_mapped, n_eroded, n_final, suvmean or nan)."""
    m = label_on_pet == label_id
    e = erode1(m)
    f = e & ~ttb
    n = int(f.sum())
    return dict(n_mapped=int(m.sum()), n_eroded=int(e.sum()), n_final=n,
                suvmean=float(pet[f].mean()) if n > 0 else np.nan)


def classify(suvmax, extent, bp, liver):
    """addendum section 3. Returns (class, bp_ge_liver_rule_applied)."""
    if extent <= EXTENT_MM:
        return "NC_partial_volume", False
    if not (np.isfinite(bp) and np.isfinite(liver)):
        return "NC_reference_unavailable", False
    if bp >= liver:
        return ("absent" if suvmax <= bp else "positive"), True
    if suvmax <= bp:
        return "absent", False
    if suvmax <= liver:
        return "low", False
    return "positive", False


def info_partial(suvmax, bp, liver):
    """transparency only for reference-unavailable components (never used as a class)."""
    s = []
    if np.isfinite(bp):
        s.append("<= blood pool" if suvmax <= bp else "> blood pool")
    if np.isfinite(liver):
        s.append("<= liver" if suvmax <= liver else "> liver")
    return "; ".join(s)


def lin1_to_kji(lin1, shape_kji):
    return np.column_stack(np.unravel_index(np.asarray(lin1, np.int64) - 1, shape_kji))


# ------------------------------ self-test ------------------------------
def selftest():
    res = []

    def ok(name, cond):
        res.append((name, bool(cond))); print(("PASS " if cond else "FAIL ") + name)

    I = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    ok("extent: single voxel = 0", max_extent(phys_lps((0, 0, 0), (2, 2, 2), I, np.array([[1, 1, 1]]))) == 0.0)
    ok("extent: two voxels along x (spacing 2.5) = 7.5",
       abs(max_extent(phys_lps((0, 0, 0), (2.5, 3, 4), I, np.array([[0, 0, 0], [0, 0, 3]]))) - 7.5) < 1e-12)
    g = np.array(np.meshgrid(range(3), range(3), range(3), indexing="ij")).reshape(3, -1).T
    ok("extent: 3x3x3 cube, spacing 2 -> 4*sqrt(3)",
       abs(max_extent(phys_lps((5, -3, 7), (2, 2, 2), I, g)) - 4 * np.sqrt(3)) < 1e-12)
    ok("extent: anisotropic spacing (4,4,3), cube 3x3x3 -> sqrt(8^2+8^2+6^2)",
       abs(max_extent(phys_lps((0, 0, 0), (4, 4, 3), I, g)) - np.sqrt(64 + 64 + 36)) < 1e-12)
    Dflip = (-1, 0, 0, 0, -1, 0, 0, 0, 1)
    ok("extent invariant to direction flip",
       abs(max_extent(phys_lps((0, 0, 0), (4, 4, 3), Dflip, g)) - np.sqrt(64 + 64 + 36)) < 1e-12)
    big = np.array(np.meshgrid(range(20), range(20), range(20), indexing="ij")).reshape(3, -1).T   # 8000 > BRUTE_MAX
    ok("extent: 20^3 cube via hull path = 19*sqrt(3)*s",
       abs(max_extent(phys_lps((0, 0, 0), (1.5, 1.5, 1.5), I, big)) - 19 * np.sqrt(3) * 1.5) < 1e-9)
    line = np.column_stack([np.zeros(4000, int), np.zeros(4000, int), np.arange(4000)])
    ok("extent: collinear 4000 points via degenerate fallback = 3999*s",
       abs(max_extent(phys_lps((0, 0, 0), (0.5, 1, 1), I, line)) - 3999 * 0.5) < 1e-9)
    plane = np.array(np.meshgrid([0], range(70), range(70), indexing="ij")).reshape(3, -1).T          # coplanar 4900
    ok("extent: coplanar 70x70 via fallback = 69*sqrt(2)",
       abs(max_extent(phys_lps((0, 0, 0), (1, 1, 1), I, plane)) - 69 * np.sqrt(2)) < 1e-9)
    rng = np.random.default_rng(1); pts = rng.normal(size=(3500, 3)) * [30, 10, 5]
    ok("extent: hull result = brute force on random cloud", abs(max_extent(pts) - _brute_max(np.unique(pts, axis=0))) < 1e-9)
    ok("10 mm rule: extent exactly 10 -> partial volume", classify(50, 10.0, 1, 5)[0] == "NC_partial_volume")
    ok("10 mm rule: 10.0001 -> classified", classify(50, 10.0001, 1, 5)[0] == "positive")

    m = np.zeros((9, 9, 9), bool); m[2:7, 2:7, 2:7] = True
    ok("erosion: 5^3 cube -> 3^3 = 27", int(erode1(m).sum()) == 27)
    m2 = np.zeros((6, 6, 6), bool); m2[0:4, 0:4, 0:4] = True
    ok("erosion: image border is background (4^3 at corner -> 2^3 at [1:3])",
       int(erode1(m2).sum()) == 8 and erode1(m2)[1:3, 1:3, 1:3].all())
    m3 = np.zeros((7, 7, 7), bool); m3[3, 3, 3] = True; m3[3, 3, 4] = True
    ok("erosion: thin structure vanishes", int(erode1(m3).sum()) == 0)
    m4 = np.zeros((9, 9, 9), bool); m4[2:7, 2:7, 2:7] = True; m4[2, 4, 4] = False   # face dent
    ok("erosion is 26-neighbour (dent on face removes adjacent-diagonal interior voxels)", int(erode1(m4).sum()) == 18)

    lab = np.zeros((9, 9, 9), np.int64); lab[2:7, 2:7, 2:7] = LABEL_LIVER
    pet = np.full((9, 9, 9), 4.0); pet[4, 4, 4] = 100.0
    ttb = np.zeros((9, 9, 9), bool); ttb[4, 4, 4] = True
    r = reference(pet, lab, LABEL_LIVER, ttb)
    ok("reference: TTB voxel excluded after erosion (27 -> 26), mean unaffected by the tumour voxel",
       r["n_mapped"] == 125 and r["n_eroded"] == 27 and r["n_final"] == 26 and r["suvmean"] == 4.0)
    ok("reference: absent label -> unavailable (nan), counts 0",
       np.isnan(reference(pet, lab, LABEL_AORTA, ttb)["suvmean"]) and reference(pet, lab, LABEL_AORTA, ttb)["n_mapped"] == 0)
    lab2 = np.zeros((9, 9, 9), np.int64); lab2[3:6, 3:6, 3:6] = LABEL_AORTA; t2 = np.zeros_like(ttb); t2[4, 4, 4] = True
    ok("reference: all eroded voxels tumour -> unavailable", np.isnan(reference(pet, lab2, LABEL_AORTA, t2)["suvmean"]))

    ok("classes: = blood pool -> absent", classify(1.5, 20, 1.5, 5)[0] == "absent")
    ok("classes: between -> low", classify(3.0, 20, 1.5, 5)[0] == "low")
    ok("classes: = liver -> low", classify(5.0, 20, 1.5, 5)[0] == "low")
    ok("classes: > liver -> positive", classify(5.0001, 20, 1.5, 5)[0] == "positive")
    ok("bp >= liver: <= bp absent, flagged", classify(4.0, 20, 4.0, 3.0) == ("absent", True))
    ok("bp >= liver: between liver and bp -> absent (no overlap)", classify(3.5, 20, 4.0, 3.0) == ("absent", True))
    ok("bp >= liver: > bp positive", classify(4.01, 20, 4.0, 3.0) == ("positive", True))
    ok("bp == liver: never low", all(classify(s, 20, 3.0, 3.0)[0] != "low" for s in (2.9, 3.0, 3.1)))
    ok("bp >= liver: low never assigned on a grid of values",
       all(classify(s, 20, 4.0, 3.0)[0] in ("absent", "positive") for s in np.linspace(0, 10, 101)))
    ok("ordinary: categories exhaustive and exclusive on a grid",
       all(classify(s, 20, 1.5, 5)[0] in ("absent", "low", "positive") for s in np.linspace(0, 10, 101)))
    ok("missing liver -> reference unavailable", classify(9, 20, 1.5, np.nan)[0] == "NC_reference_unavailable")
    ok("missing aorta -> reference unavailable", classify(9, 20, np.nan, 5)[0] == "NC_reference_unavailable")
    ok("partial volume takes precedence over missing reference", classify(9, 8, np.nan, np.nan)[0] == "NC_partial_volume")
    ok("partial info listed", info_partial(9, np.nan, 5) == "> liver" and info_partial(1, 1.5, np.nan) == "<= blood pool")
    ok("SUV >= 3 inclusive", (np.array([2.999, 3.0, 3.1]) >= SUV_COUNT).tolist() == [False, True, True])

    shape = (4, 5, 6); arr = np.arange(np.prod(shape)).reshape(shape)
    k, j, i = 2, 3, 4; lin1 = np.ravel_multi_index((k, j, i), shape) + 1
    kji = lin1_to_kji([lin1], shape)[0]
    ok("P2b linear_index_1based -> (k,j,i) C order", tuple(kji) == (k, j, i) and arr.ravel()[lin1 - 1] == arr[k, j, i])
    ok("1-based index equals MATLAB column-major index of (i,j,k)", lin1 == i + 6 * j + 6 * 5 * k + 1)

    try:
        import SimpleITK as sitk
        ct = sitk.Image([8, 8, 8], sitk.sitkUInt8); ct.SetSpacing((1, 1, 1)); ct.SetOrigin((0, 0, 0))
        a = sitk.GetArrayFromImage(ct); a[2:6, 2:6, 2:6] = LABEL_AORTA
        ct2 = sitk.GetImageFromArray(a); ct2.CopyInformation(ct)
        pet_g = sitk.Image([4, 4, 4], sitk.sitkFloat32); pet_g.SetSpacing((2, 2, 2)); pet_g.SetOrigin((0.6, 0.6, 0.6))
        on = sitk.GetArrayFromImage(sitk.Resample(ct2, pet_g, sitk.Transform(), sitk.sitkNearestNeighbor, 0, ct2.GetPixelID()))
        ok("CT label -> PET grid by physical coordinates (NN): 2x2x2 coarse voxels labelled", int((on == LABEL_AORTA).sum()) == 8
           and (on[1:3, 1:3, 1:3] == LABEL_AORTA).all())
        sh2 = sitk.Image([4, 4, 4], sitk.sitkFloat32); sh2.SetSpacing((2, 2, 2)); sh2.SetOrigin((2.6, 0.6, 0.6))
        on2 = sitk.GetArrayFromImage(sitk.Resample(ct2, sh2, sitk.Transform(), sitk.sitkNearestNeighbor, 0, ct2.GetPixelID()))
        ok("label mapping follows origin (shifted PET grid shifts labelled voxels in i)", (on2[1:3, 1:3, 0:2] == LABEL_AORTA).all()
           and int((on2 == LABEL_AORTA).sum()) == 8)
    except ImportError:
        ok("SimpleITK available for mapping tests", False)

    n_ok = sum(c for _, c in res)
    print(f"\nSELF-TEST {n_ok}/{len(res)} passed")
    return n_ok == len(res)


# ------------------------------ run ------------------------------
def grid_ok(img, row, tol=1e-3):
    g = (img.GetSize(), img.GetSpacing(), img.GetOrigin(), img.GetDirection())
    exp = ((int(row.size_x), int(row.size_y), int(row.size_z)), (row.spacing_x, row.spacing_y, row.spacing_z),
           (row.origin_x_LPS, row.origin_y_LPS, row.origin_z_LPS),
           tuple(float(row[f"dir_{a}{b}"]) for a in range(3) for b in range(3)))
    return g[0] == exp[0] and all(abs(x - y) <= tol for a, b in zip(g[1:], exp[1:]) for x, y in zip(a, b))


def same_grid(a, b, tol=1e-3):
    return a.GetSize() == b.GetSize() and max(abs(x - y) for x, y in zip(
        a.GetSpacing() + a.GetOrigin() + a.GetDirection(), b.GetSpacing() + b.GetOrigin() + b.GetDirection())) <= tol


def process(cid, VU, P1C, P1P, GR, imgs_md5):
    import SimpleITK as sitk
    dP = DATA_ROOT / cid / "PSMA"
    pet_img = sitk.ReadImage(str(dP / "PET.nii.gz"), sitk.sitkFloat64)
    ttb_img = sitk.ReadImage(str(dP / "TTB.nii.gz"))
    ts_img = sitk.ReadImage(str(dP / "totseg_24.nii.gz"))
    for f in ("PET.nii.gz", "TTB.nii.gz", "totseg_24.nii.gz"):
        imgs_md5[f"{cid}/PSMA/{f}"] = md5(dP / f)
    if not grid_ok(pet_img, GR.loc[cid]):
        raise RuntimeError(f"{cid}: PSMA PET grid differs from the P2b grid record")
    if not same_grid(pet_img, ttb_img):
        raise RuntimeError(f"{cid}: PSMA TTB not on the PSMA PET grid")
    pet = sitk.GetArrayFromImage(pet_img); ttb = sitk.GetArrayFromImage(ttb_img) > 0
    ts_on_pet = sitk.GetArrayFromImage(sitk.Resample(ts_img, pet_img, sitk.Transform(), sitk.sitkNearestNeighbor, 0,
                                                     ts_img.GetPixelID())).astype(np.int64)
    ts_arr = sitk.GetArrayFromImage(ts_img)
    rA = reference(pet, ts_on_pet, LABEL_AORTA, ttb); rL = reference(pet, ts_on_pet, LABEL_LIVER, ttb)
    thr = None
    tj = dP / "threshold.json"
    if tj.exists():
        try:
            thr = json.loads(tj.read_text()).get("suv_threshold")
        except Exception:
            thr = "unreadable"
    bp, li = rA["suvmean"], rL["suvmean"]
    prow = dict(case_id=cid, aorta_ct_voxels=int((ts_arr == LABEL_AORTA).sum()), aorta_mapped=rA["n_mapped"],
                aorta_eroded=rA["n_eroded"], aorta_final=rA["n_final"], blood_pool_SUVmean=bp,
                liver_ct_voxels=int((ts_arr == LABEL_LIVER).sum()), liver_mapped=rL["n_mapped"],
                liver_eroded=rL["n_eroded"], liver_final=rL["n_final"], liver_SUVmean=li,
                blood_pool_available=bool(np.isfinite(bp)), liver_available=bool(np.isfinite(li)),
                bp_ge_liver=bool(np.isfinite(bp) and np.isfinite(li) and bp >= li),
                psma_threshold_json=thr, psma_pet_max=float(pet.max()), n_VU=int(len(VU)))
    crows = []
    if len(VU):
        vs_path = OUT_DIR / P2B_RUN / f"P2b_voxelsets_{cid}.csv"
        if md5(vs_path) != GR.loc[cid, "csv_md5"]:
            raise RuntimeError(f"{cid}: voxel-set CSV MD5 differs from the P2b grid table")
        VS = pd.read_csv(vs_path); VS = VS[VS.role == 0]
        lab_path = OUT_DIR / P1_RUN / f"P1_{cid}_FDG_conn{CONN}.nii.gz"
        if md5(lab_path) != P1P.loc[cid, "label_image_md5"]:
            raise RuntimeError(f"{cid}: P1 label image MD5 differs from the P1 patient table")
        lab_img = sitk.ReadImage(str(lab_path)); lab = sitk.GetArrayFromImage(lab_img)
        idx = ndimage.value_indices(lab, ignore_value=0)
        for r in VU.itertuples():
            q = int(r.component_id)
            lin = VS.loc[VS.component_id == q, "linear_index_1based"].to_numpy()
            if len(lin) == 0:
                raise RuntimeError(f"{cid} c{q}: no final mapped voxel")
            if lin.min() < 1 or lin.max() > pet.size:
                raise RuntimeError(f"{cid} c{q}: mapped index outside the PSMA PET array")
            suvmax = float(pet.ravel()[lin - 1].max())
            kji = np.column_stack(idx[q])
            if len(kji) != int(P1C.loc[(cid, q), "voxels"]):
                raise RuntimeError(f"{cid} c{q}: native voxel count differs from P1")
            ext = max_extent(phys_lps(lab_img.GetOrigin(), lab_img.GetSpacing(), lab_img.GetDirection(), kji))
            cls, bpl = classify(suvmax, ext, bp, li)
            crows.append(dict(case_id=cid, component_id=q, volume_native_mL=float(r.volume_native_mL),
                              voxels_native=len(kji), voxels_mapped=len(lin), extent_mm=ext, SUVmax_PSMA=suvmax,
                              blood_pool_SUVmean=bp, liver_SUVmean=li, uptake_class=cls, bp_ge_liver_rule=bpl,
                              partial_reference_info=info_partial(suvmax, bp, li) if cls == "NC_reference_unavailable" else "",
                              SUVmax_ge_3=bool(suvmax >= SUV_COUNT)))
    return prow, crows


def summarise(C, P):
    L = [f"PSMA uptake (SAP 8.2), {len(C)} primary V_U components, {C.case_id.nunique() if len(C) else 0} patients "
         f"with V_U; references computed in {len(P)} patients", ""]
    L.append("Class                       components  native mL  patients")
    for c in CLASSES:
        s = C[C.uptake_class == c]
        L.append(f"  {c:26s}{len(s):10d}{s.volume_native_mL.sum():11.1f}{s.case_id.nunique():10d}")
    L.append(f"  {'total':26s}{len(C):10d}{C.volume_native_mL.sum():11.1f}{C.case_id.nunique():10d}")
    s3 = C[C.SUVmax_ge_3]
    L += ["", f"SUVmax >= 3 (DEEP-PSMA segmentation threshold, not a clinical criterion): {len(s3)} components, "
              f"{s3.volume_native_mL.sum():.1f} mL, {s3.case_id.nunique()} patients",
          f"Blood pool >= liver: {int(P.bp_ge_liver.sum())} patients (of whom with V_U: "
          f"{int(P[P.n_VU > 0].bp_ge_liver.sum())}); components classified under that rule: {int(C.bp_ge_liver_rule.sum())}",
          f"Blood pool unavailable: {int((~P.blood_pool_available).sum())} patients; liver unavailable: "
          f"{int((~P.liver_available).sum())} patients",
          f"Reference voxels after erosion and TTB exclusion: aorta median {P.aorta_final.median():.0f} "
          f"(min {P.aorta_final.min()}), liver median {P.liver_final.median():.0f} (min {P.liver_final.min()})",
          "Mixed-ligand cohort (68Ga-PSMA-617 / 18F-DCFPyL; per-patient ligand unavailable, not inferred). "
          "Classes are PROMISE-inspired, not PROMISE V2 scores; imaging evidence, not biological diagnoses."]
    return L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not selftest():
        sys.exit("self-test failed: no run")
    t0 = time.time()
    bad = [f"{d / n}: {md5(d / n)} != {m}" for (d, n), m in INPUTS.items() if md5(d / n) != m]
    if bad:
        sys.exit("INPUT MD5 MISMATCH:\n" + "\n".join(bad))
    P3 = pd.read_csv(OUT_DIR / P3_RUN / "P3_components_20261112_001915.csv")
    VU = P3[P3.status == "V_U"]
    if len(VU) != 217:
        sys.exit(f"expected 217 primary V_U components, found {len(VU)}")
    P1C = pd.read_csv(OUT_DIR / P1_RUN / "P1_components_20261111_233936.csv")
    P1C = P1C[P1C.connectivity == CONN].set_index(["case_id", "component_id"])
    P1P = pd.read_csv(OUT_DIR / P1_RUN / "P1_patients_20261111_233936.csv")
    P1P = P1P[P1P.connectivity == CONN].set_index("case_id")
    GR = pd.read_csv(OUT_DIR / P2B_RUN / "P2b_psma_pet_grid_20261112_000903.csv").set_index("case_id")
    cases = sorted(GR.index)
    if len(cases) != 100:
        sys.exit(f"expected 100 cases in the P2b grid table, found {len(cases)}")
    if a.preflight:
        cases = sorted(VU.case_id.unique())[:2]
    run = datetime.now().strftime("%Y%m%d_%H%M%S")
    odir = OUT_DIR / (f"PU_PREFLIGHT_{run}" if a.preflight else f"PU_{run}"); odir.mkdir()
    prows, crows, imgs = [], [], {}
    for n, cid in enumerate(cases, 1):
        pr, cr = process(cid, VU[VU.case_id == cid], P1C, P1P, GR, imgs)
        prows.append(pr); crows += cr
        print(f"[{n}/{len(cases)}] {cid}: V_U {pr['n_VU']}; aorta {pr['aorta_final']} vox, liver {pr['liver_final']} vox"
              + (f"; bp {pr['blood_pool_SUVmean']:.3f}, liver {pr['liver_SUVmean']:.3f}; threshold.json {pr['psma_threshold_json']}"
                 if a.preflight else ""), flush=True)
    P = pd.DataFrame(prows); C = pd.DataFrame(crows)
    if not a.preflight:
        if len(C) != 217 or abs(C.volume_native_mL.sum() - VU.volume_native_mL.sum()) > 1e-9:
            sys.exit("reconciliation failed: components or volume differ from the primary V_U set")
    out = {}
    for nm, df in (("components", C), ("patients", P)):
        p = odir / f"PU_{nm}_{run}.csv"; df.to_csv(p, index=False); out[p.name] = md5(p)
    S = summarise(C, P); p = odir / f"PU_summary_{run}.txt"; p.write_text("\n".join(S) + "\n", encoding="utf-8"); out[p.name] = md5(p)
    import SimpleITK as sitk, scipy
    man = dict(stage="PU PSMA uptake (SAP 8.2)", mode="PREFLIGHT" if a.preflight else "FULL", run=run,
               script=Path(__file__).name, script_md5=md5(Path(__file__)),
               inputs={f"{d / n}": m for (d, n), m in INPUTS.items()}, image_md5=imgs, cases=cases,
               params=dict(LABEL_AORTA=LABEL_AORTA, LABEL_LIVER=LABEL_LIVER, EXTENT_MM=EXTENT_MM, SUV_COUNT=SUV_COUNT,
                           erosion="3x3x3, 1 iteration, border background", mapping="NN, identity, physical coordinates"),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=scipy.__version__, SimpleITK=sitk.__version__),
               outputs=out, runtime_s=round(time.time() - t0), status="COMPLETE")
    (odir / f"PU_manifest_{run}.json").write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print("\n".join(S)); print(f"\nOutput: {odir}")


if __name__ == "__main__":
    main()
