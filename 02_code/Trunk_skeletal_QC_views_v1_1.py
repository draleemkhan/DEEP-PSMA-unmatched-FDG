"""
===========================================================================
Trunk_skeletal_QC_views_v1_1.py  -  PAPER 3: trunk skeletal registration QC (CT-to-CT review images), Amendment 2 only
===========================================================================
v1.1 (9 Nov 2026, display-only amendment; Trunk_QC_display_amendment_note_20261109.md): after preflight
  Trunk_QC_PREFLIGHT_review_20261109_112212 (v1.0, MD5 14c663249b02b67e5b99502fc95d437a) every view is cropped to a
  subregion-specific PSMA bone bounding box + 20 mm (identical array crop for PSMA and FDG; array slicing only: no
  additional resampling, smoothing or change of transforms, checks, views, slice positions or z-extent), and neck
  pages are drawn taller (row height x NECK_HEIGHT_SCALE). Crop label sets (TotalSegmentator total task):
    neck     C2..T1 vertebrae
    thorax   C7..L1 vertebrae, ribs, sternum, clavicles, scapulae   (humeri excluded: arm position is not graded)
    abdomen  T12..S1 vertebrae, ribs
  combined page: union of its parts. No crop label in the view -> full field of view, stated on the page.
Implements the LOCKED specification Trunk_skeletal_QC_spec_LOCKED_20261109.md (MD5 0151f8d44bb11661daa4b04a203d4185).
Tumour-blind: no PET, tumour mask, threshold or lesion information is read or shown. Nothing is decided here; no
transform, plane, cap or earlier QC result is changed. New file; it imports, UNCHANGED, the locked
Amendment3_regional_rigid.py and A3_visualQC_views_v1_2.py (transform reproduction, composition check, C4 image-level
check, checkerboard, display constants).

Pages: one per coded patient (TR001-TR100) and subregion, plus 15 hidden repeats (TR101-TR115), from the closed key
Trunk_QC_KEY_do_not_open_20261109.csv. Subregions (PSMA CT, physical LPS z; R6/Audit 1C-f vertebral midpoints):
  neck     C2/C3 .. C7/T1      thorax   C7/T1 .. T12/L1      abdomen  T12/L1 .. L5/S1 (L5/S1 undefined -> FOV end)
  undefined C7/T1 or T12/L1 -> the two adjacent subregions on one combined page (key: e.g. thorax+abdomen)
  C7T1/T12L1/L5S1 z from Audit1C_f_boundaries_substitution_log_20261031_164711.csv (DIRECT and SUBSTITUTED as logged);
  C2/C3 = midpoint of the PSMA z-centroids of C2 and C3 (same centroid function, verified on every DIRECT row, C5).
Condition: Amendment 2 only (T_A; supplied rigid.tfm for train_0018/0049/0067), resampling RA exactly as in the locked
  scripts. Columns: PSMA CT (reference) | FDG CT resampled | checkerboard (20 mm) | FDG bone edges (300 HU, cyan) on
  PSMA CT. Bone window [-200, 1300] HU. No segmentation contours, no case ID, no registration numbers.
Views (from PSMA anatomy only), z-extent = subregion +/- 20 mm:
  sagittal and coronal curved-slab MIP of the spine: +/-12 mm around the vertebral centroid line (centre line =
    linear interpolation in z of the PSMA C1..S1 label centroids, constant beyond the end vertebrae)
  axial slices at the PSMA centroids of C5 (neck), T4 and T8 (thorax), L1 and L3 (abdomen); full axial field of view
    (vertebra, ribs, sternum where present). A missing vertebra label -> slice omitted and stated on the page.
Integrity checks BEFORE any image; any failure STOPS the run:
  C1 MD5 of every input, imported script, specification, key and draw record
  C2 key: 100 patients, page sets as implied by the boundary log, 15 repeats (5 per type, single-subregion pages);
     draw re-run (reported; the key governs)
  C5 boundary reproduction: every DIRECT boundary in the Audit 1C-f log recomputed from the PSMA labels with the
     centroid function used for C2/C3, |difference| <= 0.01 mm
  C4 per patient and subregion: composition check (1e-4 mm) and image-level check on one vertebra (C5 / T8 / L3),
     tolerance = largest PSMA voxel dimension (as in A3_visualQC_views_v1_2.py; identity candidate, so A = B = RA)

Usage (PowerShell, folder 02_code; Amendment3_regional_rigid.py and A3_visualQC_views_v1_2.py in the same folder):
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" Trunk_skeletal_QC_views.py --selftest
  ... Trunk_skeletal_QC_views.py --preflight      TEST001/TEST002 = the same two patients as the earlier preflights
  ... Trunk_skeletal_QC_views.py                  verification of all 100, then 314 review pages
Outputs (06_output):
  Trunk_QC_[PREFLIGHT_]review_<run>\\TRxxx_<subregion>.png, RATING_TEMPLATE_Trunk_QC.csv, README_review.txt
  Trunk_QC_[PREFLIGHT_]verification_<run>.csv   (by TR code; no case IDs)
  Trunk_QC_[PREFLIGHT_]manifest_<run>.json
===========================================================================
"""
import argparse
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
_sp = importlib.util.spec_from_file_location("vqc", HERE / "A3_visualQC_views_v1_2.py")
v = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(v)
a3 = v.a3
OUT_DIR, DATA_ROOT = a3.OUT_DIR, a3.DATA_ROOT
PROTO_DIR = OUT_DIR.parent / "00_protocol"

SPEC = "Trunk_skeletal_QC_spec_LOCKED_20261109.md"
KEY = "Trunk_QC_KEY_do_not_open_20261109.csv"
DRAW = "Trunk_QC_draw_record_20261109.json"
BOUND = "Audit1C_f_boundaries_substitution_log_20261031_164711.csv"
EXPECTED_MD5 = {   # (folder, file): md5
    (HERE, "Amendment3_regional_rigid.py"): "b21a405145526d56a59bb5950ef53acf",
    (HERE, "A3_visualQC_views_v1_2.py"): "e8770c32865cf828ed27406a7fef6d9c",
    (PROTO_DIR, SPEC): "0151f8d44bb11661daa4b04a203d4185",
    (OUT_DIR, KEY): "4b2bc0868fb11378bd35ebe060bddbd0",
    (OUT_DIR, DRAW): "597589591ea5cd76b19219f4ecd91e89",
    (OUT_DIR, BOUND): "14b882a1ed3c0949c099a26291bc02c5",
    (OUT_DIR, v.FIT): "dc74686ed0b024d8524668f9b9e16c94",
    (OUT_DIR, v.TFM): "6dd27f366127d1941981489d57f0a99a",
    (OUT_DIR, v.A2P): "c66da788a6c73737163b4c09ed76f494",
    (OUT_DIR, v.KEY): "c9a4548d3daf13d3e3dadc96171a9a94",       # only to identify the two preflight patients
}
SEED = 20261109
N_REP_PER_TYPE = 5
TYPES = ("neck", "thorax", "abdomen")
TOL_BOUND_MM = 0.01
SLAB_MM = 12.0
MARGIN_MM = 20.0
# TotalSegmentator (total task) labels: C1=50 .. C7=44, T1=43 .. T12=32, L1=31 .. L5=27, S1=26 (C1/C2/S1 as in a3.LAB)
VERT = {**{f"C{i}": 51 - i for i in range(1, 8)}, **{f"T{i}": 44 - i for i in range(1, 13)},
        **{f"L{i}": 32 - i for i in range(1, 6)}, "S1": 26}
assert VERT["C1"] == a3.LAB["C1"] and VERT["C2"] == a3.LAB["C2"] and VERT["S1"] == a3.LAB["S1"]
PAIRS = {"C2C3": ("C2", "C3"), "C7T1": ("C7", "T1"), "T12L1": ("T12", "L1"), "L5S1": ("L5", "S1")}
AXIAL = {"neck": ["C5"], "thorax": ["T4", "T8"], "abdomen": ["L1", "L3"]}
C4_VERT = {"neck": "C5", "thorax": "T8", "abdomen": "L3"}
LAB_EXTRA = dict(scapula=(71, 72), clavicle=(73, 74), ribs=tuple(range(92, 116)), sternum=(116,))
assert a3.LAB["femur_L"] == 75 and a3.LAB["hip_L"] == 77 and a3.LAB["brain"] == 90   # same class map (total task)
CROP_MM = 20.0
NECK_HEIGHT_SCALE = 1.8


def _vr(a, b):
    names = list(VERT)
    return [VERT[n] for n in names[names.index(a):names.index(b) + 1]]


CROP_LABELS = {"neck": _vr("C2", "T1"),
               "thorax": _vr("C7", "L1") + list(LAB_EXTRA["ribs"] + LAB_EXTRA["sternum"] + LAB_EXTRA["clavicle"]
                                                 + LAB_EXTRA["scapula"]),
               "abdomen": _vr("T12", "S1") + list(LAB_EXTRA["ribs"])}
CENTRE_LINE = [f"C{i}" for i in range(1, 8)] + [f"T{i}" for i in range(1, 13)] + [f"L{i}" for i in range(1, 6)] + ["S1"]


# ------------------------------ pure helpers (tested) ------------------------------
def pages_for(undefined):
    """undefined: set of boundary names undefined for the patient -> list of page subregions."""
    c7, t12 = "C7T1" in undefined, "T12L1" in undefined
    if c7 and t12:
        return ["neck+thorax+abdomen"]
    if c7:
        return ["neck+thorax", "abdomen"]
    if t12:
        return ["neck", "thorax+abdomen"]
    return ["neck", "thorax", "abdomen"]


def redraw(cases_sorted, B):
    """the locked draw (spec section 3) -> key DataFrame (code, case_id, subregion, repeat_of)."""
    und = {(r.case_id, r.boundary) for r in B.itertuples() if r.status == "UNDEFINED"}
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(len(cases_sorted))
    rows = []
    for i, j in enumerate(perm):
        c = cases_sorted[j]
        for sub in pages_for({b for (cc, b) in und if cc == c}):
            rows.append(dict(code=f"TR{i + 1:03d}", case_id=c, subregion=sub, repeat_of=""))
    P = pd.DataFrame(rows)
    rep = []
    for typ in TYPES:
        el = sorted(f"{r.code}_{r.subregion}" for r in P.itertuples() if r.subregion == typ)
        rep += list(rng.choice(el, N_REP_PER_TYPE, replace=False))
    ordr = rng.permutation(len(rep))
    for k in range(len(rep)):
        code, sub = rep[ordr[k]].split("_", 1)
        cid = P[(P.code == code) & (P.subregion == sub)].case_id.iloc[0]
        rows.append(dict(code=f"TR{101 + k:03d}", case_id=cid, subregion=sub, repeat_of=code))
    return pd.DataFrame(rows)


def check_key(K, cases_sorted, B):
    """C2: structural checks (stop) + draw reproduction (reported)."""
    prob = []
    und = {(r.case_id, r.boundary) for r in B.itertuples() if r.status == "UNDEFINED"}
    prim, rep = K[K.repeat_of == ""], K[K.repeat_of != ""]
    pc = prim.groupby("code").case_id.agg(set)
    if sorted(pc.index) != [f"TR{i:03d}" for i in range(1, 101)]:
        prob.append("primary codes are not TR001-TR100")
    if any(len(s) != 1 for s in pc) or sorted(c for s in pc for c in s) != sorted(cases_sorted):
        prob.append("primary codes do not map one-to-one onto the 100 cases")
    for cid, g in prim.groupby("case_id"):
        if sorted(g.subregion) != sorted(pages_for({b for (cc, b) in und if cc == cid})):
            prob.append(f"page set wrong for a patient ({g.code.iloc[0]})")
    if sorted(rep.code) != [f"TR{i:03d}" for i in range(101, 116)]:
        prob.append("repeat codes are not TR101-TR115")
    if rep.subregion.value_counts().to_dict() != {t: N_REP_PER_TYPE for t in TYPES}:
        prob.append("repeats not 5 per subregion type")
    for r in rep.itertuples():
        src = prim[(prim.code == r.repeat_of) & (prim.subregion == r.subregion)]
        if len(src) != 1 or src.case_id.iloc[0] != r.case_id:
            prob.append(f"repeat {r.code} does not match its source page")
    if rep.duplicated(["repeat_of", "subregion"]).any():
        prob.append("a page is repeated twice")
    R = redraw(cases_sorted, B)
    same = R.reset_index(drop=True).equals(K[["code", "case_id", "subregion", "repeat_of"]].reset_index(drop=True))
    return prob, bool(same)


def subregion_z(sub, planes, z_fov_lo):
    """(z_lo, z_hi) in PSMA physical z for a (possibly combined) subregion. planes: dict C2C3/C7T1/T12L1/L5S1 -> z or
    nan. Undefined L5/S1 -> FOV bottom. Tie rule (SAP R.3): a point exactly on a plane belongs to the trunk, so each
    subregion is closed at its trunk-side ends (display only; V_I assignment is done elsewhere)."""
    parts = sub.split("+")
    top = {"neck": "C2C3", "thorax": "C7T1", "abdomen": "T12L1"}[parts[0]]
    bot = {"neck": "C7T1", "thorax": "T12L1", "abdomen": "L5S1"}[parts[-1]]
    z_hi = planes[top]
    z_lo = z_fov_lo if (bot == "L5S1" and not np.isfinite(planes[bot])) else planes[bot]   # FOV only for L5/S1
    if not (np.isfinite(z_hi) and np.isfinite(z_lo)) or z_lo >= z_hi:
        raise ValueError(f"subregion {sub}: invalid extent ({z_lo}, {z_hi})")
    return float(z_lo), float(z_hi)


def phys_kji(im, kji):
    """physical LPS coordinates of voxel centres given as (k, j, i) array indices of a SimpleITK image."""
    D = np.array(im.GetDirection()).reshape(3, 3); S = np.array(im.GetSpacing())
    return np.array(im.GetOrigin()) + (np.asarray(kji, float)[:, ::-1] * S) @ D.T


def label_centroid(im, arr, lab):
    """physical centroid (x, y, z) of all voxels with label lab; None if absent."""
    idx = np.argwhere(arr == lab)
    return None if len(idx) == 0 else phys_kji(im, idx).mean(0)


def midpoint_z(im, arr, a, b):
    ca, cb = label_centroid(im, arr, VERT[a]), label_centroid(im, arr, VERT[b])
    return np.nan if ca is None or cb is None else float((ca[2] + cb[2]) / 2)


def centre_line(lab_lps, nz):
    """vertebral centroid line on an LPS-oriented label array (z, y, x): per-slice (yc, xc) in voxel units, linear in z
    between label centroids, constant beyond the end ones. None if fewer than 2 vertebrae are labelled."""
    pts = []
    for nm in CENTRE_LINE:
        idx = np.argwhere(lab_lps == VERT[nm])
        if len(idx):
            pts.append(idx.mean(0))
    if len(pts) < 2:
        return None
    pts = np.array(sorted(pts, key=lambda p: p[0]))
    k = np.arange(nz)
    return np.interp(k, pts[:, 0], pts[:, 1]), np.interp(k, pts[:, 0], pts[:, 2])


def curved_mip(vol, line, sp_xyz, plane, slab_mm=SLAB_MM):
    """curved-slab MIP of vol (z, y, x, LPS) around the centre line: 'sag' -> max over x within +/-slab of xc(z),
    image (z, y); 'cor' -> max over y within +/-slab of yc(z), image (z, x). Superior up (flipud)."""
    yc, xc = line
    nz, ny, nx = vol.shape
    out = np.full((nz, ny if plane == "sag" else nx), -1000.0, np.float32)
    for k in range(nz):
        if plane == "sag":
            m = np.abs(np.arange(nx) - xc[k]) * sp_xyz[0] <= slab_mm
            if m.any():
                out[k] = vol[k][:, m].max(axis=1)
        else:
            m = np.abs(np.arange(ny) - yc[k]) * sp_xyz[1] <= slab_mm
            if m.any():
                out[k] = vol[k][m, :].max(axis=0)
    return np.flipud(out)


def crop_box(lab_zyx, labels, sp_xyz, margin_mm=CROP_MM):
    """(y0, y1, x0, x1) voxel slice bounds of the in-plane bounding box of the given labels (any z) + margin; None if
    none present. Pure array indices: cropping is slicing only."""
    m = np.isin(lab_zyx, labels)
    if not m.any():
        return None
    ys, xs = np.nonzero(m.any(axis=0))
    gy, gx = int(np.ceil(margin_mm / sp_xyz[1])), int(np.ceil(margin_mm / sp_xyz[0]))
    ny, nx = lab_zyx.shape[1:]
    return max(ys.min() - gy, 0), min(ys.max() + gy + 1, ny), max(xs.min() - gx, 0), min(xs.max() + gx + 1, nx)


def render_page4(title, rows, fn, height_scale=1.0):
    """rows: list of (view_name, sp_rc, P2, F2). 4 columns; no mask is ever passed here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["PSMA CT (reference)", "FDG CT (Amendment 2)", "checkerboard", "FDG bone edges on PSMA CT"]
    asp = [max(0.35, min(3.0, (P2.shape[0] * s[0]) / (P2.shape[1] * s[1]))) for _, s, P2, _ in rows]
    fig, axs = plt.subplots(len(rows), 4, figsize=(18, 0.6 + sum(4.2 * height_scale * a for a in asp)), squeeze=False,
                            gridspec_kw=dict(height_ratios=asp))
    for r, (name, sp_rc, P2, F2) in enumerate(rows):
        ext = (0, P2.shape[1] * sp_rc[1], 0, P2.shape[0] * sp_rc[0])
        panels = [P2, F2, v.checkerboard(P2, F2, sp_rc), P2]
        for j in range(4):
            ax = axs[r, j]
            ax.imshow(panels[j], cmap="gray", vmin=v.BONE_W[0], vmax=v.BONE_W[1], extent=ext, origin="upper",
                      aspect="equal")
            if j == 3 and (F2 > v.EDGE_HU).any() and (F2 <= v.EDGE_HU).any():
                ax.contour(F2, levels=[v.EDGE_HU], colors=["cyan"], linewidths=0.6, extent=ext, origin="upper")
            ax.set_title(f"{cols[j]} | {name}" if j == 0 else cols[j], fontsize=8)
            ax.tick_params(labelsize=5)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    return fn


# ------------------------------ patient-level work ------------------------------
def patient_setup(cid, PP, TF):
    """transforms (locked reproduction), RA with composition check, PSMA label image/array, planes."""
    import SimpleITK as sitk
    X = v.patient_transforms(cid, PP, TF); X["cid"] = cid
    RA, RB, comp = v.resampling_transforms(X["dP"], X["t"], np.eye(3), np.zeros(3), X["A"], X["b"],
                                           cid in a3.FALLBACKS)
    X.update(RA=RA, compose_err_mm=comp)
    imP = sitk.ReadImage(str(X["dP"] / "totseg_24.nii.gz"))
    return X, imP, sitk.GetArrayFromImage(imP)


def boundary_planes(cid, imP, aP, B):
    """C2C3 computed; C7T1/T12L1/L5S1 from the log; plus C5 reproduction rows for DIRECT boundaries."""
    planes, c5 = {"C2C3": midpoint_z(imP, aP, "C2", "C3")}, []
    for r in B[B.case_id == cid].itertuples():
        planes[r.boundary] = float(r.z_mm) if r.status != "UNDEFINED" else np.nan
        if r.status == "DIRECT":
            z = midpoint_z(imP, aP, *PAIRS[r.boundary])
            c5.append(dict(boundary=r.boundary, logged=float(r.z_mm), repro=z, diff=abs(z - float(r.z_mm))))
    return planes, c5


def fov_z_lo(im):
    S = np.array(im.GetSize())
    corners = np.array([[i, j, k] for i in (0, S[0] - 1) for j in (0, S[1] - 1) for k in (0, S[2] - 1)], float)
    return float(phys_kji(im, corners[:, ::-1])[:, 2].min())


def build_rows(sub, X, imP, planes):
    """PSMA sub-grid (LPS, full axial FOV, z = subregion +/- margin); FDG CT resampled with RA; views."""
    import SimpleITK as sitk
    z_lo, z_hi = subregion_z(sub, planes, fov_z_lo(imP))
    L = sitk.DICOMOrient(imP, "LPS")
    nz = L.GetSize()[2]
    zk = np.array([L.TransformIndexToPhysicalPoint((0, 0, int(k)))[2] for k in range(nz)])
    keep = np.nonzero((zk >= z_lo - MARGIN_MM) & (zk <= z_hi + MARGIN_MM))[0]
    if len(keep) < 2:
        return [], [f"subregion {sub} not covered by the PSMA CT"], (z_lo, z_hi)
    ref = sitk.RegionOfInterest(L, [L.GetSize()[0], L.GetSize()[1], int(keep.max() - keep.min() + 1)],
                                [0, 0, int(keep.min())])
    ctP = sitk.ReadImage(str(X["dP"] / "CT.nii.gz")); ctF = sitk.ReadImage(str(X["dF"] / "CT.nii.gz"))
    Pa = sitk.GetArrayFromImage(sitk.Resample(ctP, ref, sitk.Transform(), sitk.sitkLinear, -1000, sitk.sitkFloat32))
    Fa = sitk.GetArrayFromImage(sitk.Resample(ctF, ref, X["RA"], sitk.sitkLinear, -1000, sitk.sitkFloat32))
    lab = sitk.GetArrayFromImage(ref)
    labfull = sitk.GetArrayFromImage(L)
    sp = np.array(ref.GetSpacing())
    rows, notes = [], []
    box = crop_box(lab, sorted({l for part in sub.split("+") for l in CROP_LABELS[part]}), sp)
    if box is None:
        notes.append("no crop label in the view: full field of view shown")
        box = (0, lab.shape[1], 0, lab.shape[2])
    y0, y1, x0, x1 = box
    full = centre_line(labfull, nz)
    if full is None:
        notes.append("vertebral centre line not available: curved MIPs omitted")
    else:
        line = (full[0][keep.min():keep.max() + 1], full[1][keep.min():keep.max() + 1])
        for pl, nm in (("sag", "sagittal curved-slab MIP +/-12 mm"), ("cor", "coronal curved-slab MIP +/-12 mm")):
            sl = np.s_[:, y0:y1] if pl == "sag" else np.s_[:, x0:x1]
            rows.append((nm, v.plane_spacing(sp, pl), curved_mip(Pa, line, sp, pl)[sl], curved_mip(Fa, line, sp, pl)[sl]))
    for part in sub.split("+"):
        for vn in AXIAL[part]:
            idx = np.argwhere(lab == VERT[vn])
            if len(idx) == 0:
                notes.append(f"{vn} not labelled in the PSMA CT within the view: axial {vn} omitted")
                continue
            k = int(round(idx[:, 0].mean()))
            rows.append((f"axial {vn} centroid", v.plane_spacing(sp, "ax"), Pa[k, y0:y1, x0:x1], Fa[k, y0:y1, x0:x1]))
    return rows, notes, (z_lo, z_hi)


README = """Trunk skeletal registration QC - review instructions (spec {spec}; locked)
Skeletal alignment only. No PET, no lesion information. Organ position is NOT graded.
Each page = one coded patient (TRxxx) x one subregion (neck / thorax / abdomen; a combined page, e.g. thorax+abdomen,
is graded once and its grade applies to both parts). One condition only (Amendment 2): there is no X/Y comparison.
Columns: PSMA CT (reference) | FDG CT resampled | checkerboard ({sq:.0f} mm squares) | FDG bone edges ({edge:.0f} HU
contour, cyan) on PSMA CT. Bone window -200..1300 HU. Views: sagittal and coronal curved-slab MIP of the spine
(+/-12 mm around the vertebral centroid line) over the subregion +/- 20 mm; axial slices at C5 (neck), T4 and T8
(thorax), L1 and L3 (abdomen). Display v1.1: every view cropped to the subregion's PSMA bone bounding box + 20 mm
(same crop for PSMA and FDG; display only); neck pages drawn taller.
Grade each page (RATING_TEMPLATE_Trunk_QC.csv) on skeletal landmarks only:
  3  Excellent  - no discernible mismatch of corresponding skeletal landmarks
  2  Acceptable - minor residual mismatch, no meaningful misregistration
  1  Poor       - definite skeletal misalignment
  0  Failure    - gross misregistration or wrong correspondence (e.g. vertebral-level slip)
  NA Not assessable - insufficient coverage to judge
Landmarks: neck C3-C7 bodies, posterior elements, alignment; thorax T1-T12 bodies, posterior elements, alignment
(supporting: ribs, sternum); abdomen L1-L5 bodies, posterior elements, alignment (supporting: lower ribs).
Also record: vertebral_level_slip yes/no; rotation_evident yes/no; landmarks driving the grade; neck_posture_diff
yes/no (neck pages only); additional views requested (with the prompting landmark); comments.
Do not penalise arm position or organ position.
Do NOT open Trunk_QC_KEY_do_not_open_20261109.csv until grading of ALL pages is complete. Grade every page on its
own, without looking back at other pages.
"""
TEMPLATE_COLS = ["code", "subregion", "grade", "landmarks", "vertebral_level_slip", "rotation_evident",
                 "neck_posture_diff", "additional_views_requested", "comment"]


# ------------------------------ self-test (synthetic; no patient data) ------------------------------
def selftest():
    import SimpleITK as sitk
    res = []
    ok = lambda name, cond: res.append((name, bool(cond)))
    # 1 page sets
    ok("pages: none undefined", pages_for(set()) == ["neck", "thorax", "abdomen"])
    ok("pages: T12L1 undefined -> thorax+abdomen", pages_for({"T12L1", "L5S1"}) == ["neck", "thorax+abdomen"])
    ok("pages: C7T1 undefined -> neck+thorax", pages_for({"C7T1"}) == ["neck+thorax", "abdomen"])
    ok("pages: L5S1 alone does not combine", pages_for({"L5S1"}) == ["neck", "thorax", "abdomen"])
    # 2 draw: synthetic cohort, structure and determinism
    cs = [f"case_{i:03d}" for i in range(100)]
    Bs = pd.DataFrame([dict(case_id=c, boundary=b, status="DIRECT", detail="", z_mm=0.0) for c in cs
                       for b in ("C7T1", "T12L1", "L5S1")])
    Bs.loc[(Bs.case_id == "case_007") & Bs.boundary.isin(["T12L1", "L5S1"]), "status"] = "UNDEFINED"
    K1, K2 = redraw(cs, Bs), redraw(cs, Bs)
    ok("draw deterministic", K1.equals(K2))
    prob, same = check_key(K1, cs, Bs)
    ok("draw passes its own key check", not prob and same)
    ok("draw: 299 primary + 15 repeats", (K1.repeat_of == "").sum() == 299 and (K1.repeat_of != "").sum() == 15)
    ok("draw: combined page never repeated", not K1[K1.repeat_of != ""].subregion.str.contains(r"\+").any())
    Kt = K1.copy(); Kt.loc[Kt.index[-1], "case_id"] = "case_999"
    ok("key check detects a corrupted repeat", bool(check_key(Kt, cs, Bs)[0]))
    Kt = K1.copy(); Kt.loc[0, "subregion"] = "thorax"
    ok("key check detects a wrong page set", bool(check_key(Kt, cs, Bs)[0]))
    # 3 subregion extents
    pl = dict(C2C3=100.0, C7T1=50.0, T12L1=-100.0, L5S1=np.nan)
    ok("extent neck", subregion_z("neck", pl, -300) == (50.0, 100.0))
    ok("extent abdomen, L5S1 undefined -> FOV", subregion_z("abdomen", pl, -300) == (-300.0, -100.0))
    ok("extent combined", subregion_z("thorax+abdomen", dict(pl, T12L1=np.nan), -300) == (-300.0, 50.0))
    try:
        subregion_z("thorax", dict(pl, T12L1=np.nan), -300); ok("undefined T12L1 on a single page raises", False)
    except ValueError:
        ok("undefined T12L1 on a single page raises", True)
    # 4 centroid in physical space with a flipped, anisotropic grid
    arr = np.zeros((20, 10, 10), np.uint8); arr[4:6, 3:5, 3:5] = VERT["C2"]; arr[10:12, 3:5, 3:5] = VERT["C3"]
    im = sitk.GetImageFromArray(arr); im.SetSpacing((1.0, 1.0, 2.5)); im.SetOrigin((0.0, 0.0, 50.0))
    im.SetDirection((1, 0, 0, 0, 1, 0, 0, 0, -1))
    ok("midpoint z (flipped z, 2.5 mm)", abs(midpoint_z(im, arr, "C2", "C3") - (50 - 2.5 * 7.5)) < 1e-9)
    ok("absent label -> nan", np.isnan(midpoint_z(im, arr, "C2", "C4")))
    # 5 centre line and curved MIP
    lab = np.zeros((30, 40, 40), np.uint8); lab[2:4, 18:20, 10:12] = VERT["L3"]; lab[26:28, 18:20, 28:30] = VERT["T8"]
    line = centre_line(lab, 30)
    ok("centre line interpolates", abs(line[1][15] - (10.5 + (15 - 2.5) / 24 * 18)) < 1e-9)
    vol = np.full((30, 40, 40), -1000.0, np.float32)
    for k in range(30):
        vol[k, 20, int(round(line[1][k]))] = 1000.0          # bright rod on the curve
        vol[k, 20, 39] = 2000.0                              # decoy far outside the slab
    sag = curved_mip(vol, line, (1.0, 1.0, 1.0), "sag", slab_mm=3)
    ok("curved MIP finds the rod, excludes the decoy", np.all(sag[:, 20] == 1000.0) and sag.max() == 1000.0)
    ok("curved MIP superior up", curved_mip(np.arange(30, dtype=np.float32)[:, None, None] * np.ones((1, 4, 4)),
                                            (np.full(30, 2.0), np.full(30, 2.0)), (1, 1, 1), "cor")[0, 0] == 29)
    ok("centre line None if < 2 vertebrae", centre_line(np.zeros((5, 5, 5), np.uint8), 5) is None)
    # 6 composition with identity candidate and C4 on a synthetic pair (identity transforms)
    a = np.zeros((30, 30, 30), np.uint8); a[10:20, 10:20, 10:20] = VERT["T8"]
    imA = sitk.GetImageFromArray(a); imA.SetSpacing((2.0, 2.0, 3.0))
    T = sitk.CompositeTransform(3); T.AddTransform(sitk.TranslationTransform(3, (0.0, 0.0, 0.0)))
    w, tol, st = v.image_level_check(imA, imA, T, T, np.eye(3), np.zeros(3), np.zeros(3), np.eye(3), np.zeros(3),
                                     VERT["T8"])
    ok("C4 identity: error ~0, tol = max spacing", st == "ok" and w < 1e-6 and tol == 3.0)
    # 7 rendering (synthetic), no identifiers
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        P2 = np.random.default_rng(0).uniform(-200, 1300, (40, 60)).astype(np.float32)
        fn = render_page4("TR000 | thorax | selftest", [("sagittal", (2.0, 1.0), P2, P2 + 50)], Path(d) / "t.png")
        ok("render 4 columns", fn.exists() and fn.stat().st_size > 1000)
    ok("template has no case_id column", "case_id" not in TEMPLATE_COLS)
    # 8 v1.1 crop (display only)
    labc = np.zeros((4, 100, 80), np.uint8); labc[1, 40:50, 30:35] = VERT["C5"]; labc[2, 55:60, 20:25] = 92
    ok("crop box: union over z, margin in mm", crop_box(labc, [VERT["C5"], 92], (2.0, 2.0, 3.0)) == (30, 70, 10, 45))
    ok("crop box clipped to the grid", crop_box(labc, [VERT["C5"]], (1.0, 1.0, 1.0), 100) == (0, 100, 0, 80))
    ok("crop box None if no label", crop_box(labc, [VERT["L3"]], (1.0, 1.0, 1.0)) is None)
    volc = np.random.default_rng(1).normal(size=(4, 100, 80)).astype(np.float32)
    y0, y1, x0, x1 = crop_box(labc, [VERT["C5"], 92], (2.0, 2.0, 3.0))
    ok("crop is slicing only (values unchanged)", np.array_equal(volc[2, y0:y1, x0:x1], volc[2][30:70, 10:45]))
    ok("crop label sets by region", VERT["C3"] in CROP_LABELS["neck"] and 116 in CROP_LABELS["thorax"]
       and 69 not in CROP_LABELS["thorax"] and VERT["L5"] in CROP_LABELS["abdomen"] and 92 not in CROP_LABELS["neck"])
    for n, g in res:
        print(f"  [{'PASS' if g else 'FAIL'}] {n}")
    allok = all(g for _, g in res)
    print(f"selftest: {sum(g for _, g in res)}/{len(res)} passed")
    return allok


# ------------------------------ run ------------------------------
def run(preflight=False):
    import SimpleITK as sitk
    t0 = time.time()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    md = {name: ((a3.md5(d / name) if (d / name).exists() else "MISSING"), exp) for (d, name), exp in EXPECTED_MD5.items()}
    bad = {k: got for k, (got, exp) in md.items() if got != exp}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(md)}/{len(md)} match")
    FITT = pd.read_csv(OUT_DIR / v.FIT); TF = pd.read_csv(OUT_DIR / v.TFM)
    PP = pd.read_csv(OUT_DIR / v.A2P).set_index("case_id")
    B = pd.read_csv(OUT_DIR / BOUND)
    cases = sorted(FITT["case_id"].unique())
    K = pd.read_csv(OUT_DIR / KEY, keep_default_na=False)
    prob, redraw_ok = check_key(K, cases, B)
    if prob:
        sys.exit(f"STOP C2: key problems: {prob[:5]}")
    print(f"C2 key: 100 patients, {(K.repeat_of == '').sum()} pages + {(K.repeat_of != '').sum()} repeats OK; draw "
          f"re-run in numpy {np.__version__}: {'REPRODUCED' if redraw_ok else 'NOT reproduced (the locked key governs)'}")
    if preflight:
        oldK = pd.read_csv(OUT_DIR / v.KEY)
        others = [c for c in cases if c not in set(oldK["case_id"])][:2]          # same two as on 5 Nov
        und = {(r.case_id, r.boundary) for r in B.itertuples() if r.status == "UNDEFINED"}
        work = pd.DataFrame([dict(code=f"TEST{i + 1:03d}", case_id=c, subregion=s, repeat_of="")
                             for i, c in enumerate(others) for s in pages_for({b for (cc, b) in und if cc == c})])
        tag = "PREFLIGHT_"
    else:
        work, tag = K, ""
    pcode = {r.case_id: r.code for r in work.itertuples() if r.repeat_of == ""}
    # verification BEFORE any image: composition, C5 boundaries, C4 per subregion
    ver, c5rows, cache = [], [], {}
    for cid in sorted(pcode, key=lambda c: pcode[c]):
        X, imP, aP = patient_setup(cid, PP, TF)
        planes, c5 = boundary_planes(cid, imP, aP, B)
        for d in c5:
            c5rows.append(dict(code=pcode[cid], **d))
        imF = sitk.ReadImage(str(X["dF"] / "totseg_24.nii.gz"))
        for sub in work[(work.case_id == cid) & (work.repeat_of == "")].subregion:
            for part in sub.split("+"):
                w, tol, st = v.image_level_check(imP, imF, X["RA"], X["RA"], X["A"], X["b"], X["t"], np.eye(3),
                                                 np.zeros(3), VERT[C4_VERT[part]])
                ver.append(dict(code=pcode[cid], case=cid, subregion=sub, c4_vertebra=C4_VERT[part],
                                compose_err_mm=X["compose_err_mm"], image_check_worst_mm=w, image_check_tol_mm=tol,
                                image_check_status=st,
                                C4_pass=bool(X["compose_err_mm"] <= v.TOL_COMPOSE_MM
                                             and (st != "ok" or (np.isfinite(w) and w <= tol)))))
        cache[cid] = planes
        del imP, aP, imF
        print(f"  verified {pcode[cid]}: C4 {'PASS' if all(d['C4_pass'] for d in ver if d['case'] == cid) else 'FAIL'}"
              f", C5 {'PASS' if all(d['diff'] <= TOL_BOUND_MM for d in c5) else 'FAIL'}  ({time.time() - t0:.0f} s)")
    V = pd.DataFrame(ver)
    C5 = pd.DataFrame(c5rows)
    C5["C5_pass"] = C5["diff"] <= TOL_BOUND_MM
    vp = OUT_DIR / f"Trunk_QC_{tag}verification_{run_id}.csv"
    V.drop(columns="case").to_csv(vp, index=False)                       # by code only; no case IDs
    cp = OUT_DIR / f"Trunk_QC_{tag}boundary_check_{run_id}.csv"
    C5.to_csv(cp, index=False)
    print(f"C5 boundaries: {int(C5.C5_pass.sum())}/{len(C5)} DIRECT reproduced, max |diff| {C5['diff'].max():.2e} mm")
    stop = sorted(set(V.loc[~V.C4_pass, "code"]) | set(C5.loc[~C5.C5_pass, "code"]))
    if stop:
        sys.exit(f"STOP: verification failed for {stop}; no image generated. See {vp.name} and {cp.name}")
    print("C4 and C5: PASS -> generating images")
    odir = OUT_DIR / f"Trunk_QC_{tag}review_{run_id}"
    odir.mkdir()
    files, tmpl = [], []
    for cid in sorted(pcode, key=lambda c: pcode[c]):
        X, imP, _ = patient_setup(cid, PP, TF)
        for sub in sorted(set(work[work.case_id == cid].subregion)):
            rows, notes, (zlo, zhi) = build_rows(sub, X, imP, cache[cid])
            for kr in work[(work.case_id == cid) & (work.subregion == sub)].itertuples():
                if not rows:
                    print(f"  {kr.code} {sub}: no image ({'; '.join(notes)})"); continue
                title = (f"{kr.code} | {sub} | Amendment 2 | PSMA CT reference, FDG CT resampled | bone window "
                         f"{v.BONE_W[0]}..{v.BONE_W[1]} HU | no segmentation contours"
                         + (f"\nnote: {'; '.join(notes)}" if notes else ""))
                files.append(render_page4(title, rows, odir / f"{kr.code}_{sub}.png",
                                          NECK_HEIGHT_SCALE if sub == "neck" else 1.0))
                tmpl.append(dict(code=kr.code, subregion=sub, **{k: "" for k in TEMPLATE_COLS[2:]}))
        del imP
        print(f"  {pcode[cid]}: images written ({time.time() - t0:.0f} s)")
    tp = odir / "RATING_TEMPLATE_Trunk_QC.csv"
    pd.DataFrame(sorted(tmpl, key=lambda d: (d["code"], d["subregion"])), columns=TEMPLATE_COLS).to_csv(tp, index=False)
    rp = odir / "README_review.txt"
    rp.write_text(README.format(spec=SPEC, sq=v.CHECKER_MM, edge=v.EDGE_HU), encoding="utf-8")
    man = dict(step=f"Trunk skeletal registration QC images ({'PREFLIGHT' if preflight else 'FULL'})", run=run_id,
               spec=SPEC, script=Path(__file__).name, script_md5=a3.md5(Path(__file__)),
               inputs={k: got for k, (got, _) in md.items()}, draw_reproduced=bool(redraw_ok),
               verification=vp.name, verification_md5=a3.md5(vp), boundary_check=cp.name, boundary_check_md5=a3.md5(cp),
               n_pages=len(files), outputs={p.name: a3.md5(p) for p in files + [tp, rp]},
               params=dict(CROP_MM=CROP_MM, NECK_HEIGHT_SCALE=NECK_HEIGHT_SCALE, CROP_LABELS=CROP_LABELS, SEED=SEED, N_REP_PER_TYPE=N_REP_PER_TYPE, SLAB_MM=SLAB_MM, MARGIN_MM=MARGIN_MM,
                           TOL_BOUND_MM=TOL_BOUND_MM, TOL_COMPOSE_MM=v.TOL_COMPOSE_MM, BONE_W=v.BONE_W,
                           EDGE_HU=v.EDGE_HU, CHECKER_MM=v.CHECKER_MM, AXIAL=AXIAL, C4_VERT=C4_VERT),
               runtime_s=round(time.time() - t0),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             SimpleITK=sitk.Version_VersionString()), status="COMPLETE")
    (OUT_DIR / f"Trunk_QC_{tag}manifest_{run_id}.json").write_text(json.dumps(man, indent=2, default=str),
                                                                   encoding="utf-8")
    print(f"done: {len(files)} pages in {odir.name}; {time.time() - t0:.0f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(preflight=a.preflight)


if __name__ == "__main__":
    main()
