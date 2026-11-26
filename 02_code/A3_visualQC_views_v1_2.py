"""
===========================================================================
A3_visualQC_views.py  -  PAPER 3: Amendment 3a/3b targeted visual registration QC (CT-to-CT review images)
===========================================================================
Implements the LOCKED specification A3_visual_QC_spec_v1_0_20261105.md (MD5 10bea05a1792ffdc4401f5825b322c15).
Level 2 only (anatomical validity). No PET, tumour mask, threshold or lesion information is read or shown.
Nothing is decided here; no algorithm, parameter, cap, mask or pilot result is changed.

For each of the 30 coded patients (VR001-VR030, from the closed key) and each region (head, pelvis):
  condition A2  = Amendment 2 (supplied rigid.tfm for the fallback patients train_0018/0049/0067)
  condition A3a = Amendment 3a CANDIDATE from Amendment3_transforms_20261102_073040.csv (adopted OR rejected)
  shown as X / Y in the randomised order stored in the key. PSMA CT = reference; FDG CT resampled into PSMA space
  (linear). Bone window [-200, 1300] HU. Per slice: PSMA CT | FDG (X) | checkerboard (X) | FDG bone edges on PSMA (X)
  | FDG (Y) | checkerboard (Y) | edges (Y). NO TotalSegmentator contours are drawn.
  Slice positions are defined from PSMA anatomy only and are identical for X and Y:
    head:   axial skull-base candidates at the 5th and 12th percentiles of PSMA brain z | axial petrous level
            (densest bilateral skull-label bone > 1300 HU lateral to the midline, brain z P0-P40; display aid; if
            not found the page is FLAGGED for an additional view, nothing substituted) | axial upper vault (brain z
            P80) | coronal and sagittal through the PSMA skull centroid
            (revision after preflight 20261105_135552: the single P20 slice missed the skull base)
    pelvis: axial acetabular roofs (mean PSMA femoral-head centre z + 0.85 x radius, Amendment 3b table) | axial S1
            centroid | coronal through both PSMA femoral-head centres | sagittal through each femoral-head centre
  (fallback where a centre is unavailable: hip-bone label centroid; stated on the image)

Transform handling reproduces the locked Amendment3_regional_rigid.py exactly (imported unchanged):
  point mapping FDG -> PSMA:  T_A(p) = FDG/rigid.tfm(p) - t ;  candidate dT(T_A(p)) = R T_A(p) + c
  resampling (PSMA grid -> FDG): RA = PSMA/rigid.tfm o (+t) ;  RB = RA o dT^-1   (as in render_overlay)

Integrity checks BEFORE any image (spec section 9); any failure STOPS the run:
  C1 MD5 of every input, the locked script and the specification/key/draw files
  C2 key: 30 codes, 30 unique cases, X/Y assignments complete; control draw re-run (reported; the key governs)
  C3 per patient and region: stored fov_zlo/fov_zhi and sel_centroid_A2 / sel_centroid_A3 reproduced from the masks
     with the locked functions, |difference| <= 0.01 mm (TOL_MM)
  C4 per patient and region: composition check RB(dT(T_A(p))) = p (1e-4 mm) and image-level check (centroid of the
     selection label resampled with RA / RB vs point-mapped, tolerance = largest PSMA voxel dimension)

Usage (Anaconda PowerShell, env totalseg, folder 02_code, Amendment3_regional_rigid.py in the same folder):
  python A3_visualQC_views.py --selftest
  python A3_visualQC_views.py --preflight          2 patients OUTSIDE the 30 (TEST codes, separate folder)
  python A3_visualQC_views.py                      verification of all 30, then review images
  python A3_visualQC_views.py --overlap-survey     descriptive FOV overlap of skull/brain, all 100 patients (added
                                                   5 Nov 2026; no images; console = cohort distributions only)
  python A3_visualQC_views.py --diagnose VR009     read-only C4 diagnostic for one code (added 5 Nov 2026 after
                                                   the VR009 head C4 failure; no images, no files, no case ID;
                                                   C1-C4 unchanged)
Outputs (06_output):
  A3_visualQC_review_<run>\\VRxxx_head.png, VRxxx_pelvis.png, RATING_TEMPLATE_A3_visualQC.csv, README_review.txt
  A3_visualQC_verification_<run>.csv  (by VR code; no case IDs)
  A3_visualQC_manifest_<run>.json     (MD5s of inputs and outputs, parameters, software)
===========================================================================
"""
import argparse
import hashlib
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
_sp = importlib.util.spec_from_file_location("a3", HERE / "Amendment3_regional_rigid.py")
a3 = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(a3)
OUT_DIR, DATA_ROOT = a3.OUT_DIR, a3.DATA_ROOT

SPEC = "A3_visual_QC_spec_v1_0_20261105.md"
EXPECTED_MD5 = {
    "Amendment3_regional_rigid.py": "b21a405145526d56a59bb5950ef53acf",
    "Amendment3_fit_20261102_073040.csv": "dc74686ed0b024d8524668f9b9e16c94",
    "Amendment3_transforms_20261102_073040.csv": "6dd27f366127d1941981489d57f0a99a",
    "Audit1C_e_per_patient_20261031_141128.csv": "c66da788a6c73737163b4c09ed76f494",
    "Amendment3b_femhead_20261102_211242/Amendment3b_femhead_sides_20261102_211242.csv": "8aeff038bda56d60b4f396ca82f40374",
    SPEC: "10bea05a1792ffdc4401f5825b322c15",
    "A3_visualQC_KEY_do_not_open_20261105.csv": "c9a4548d3daf13d3e3dadc96171a9a94",
    "A3_visualQC_draw_record_20261105.json": "bae9d1ba252c48f1cf1d0bf1cc595419",
    "A3_visualQC_control_pool_20261105.csv": "7867e06a790c412a420c62a2a0c3e49a",
}
NOTE_VR009 = "A3_visualQC_note_20261105_VR009_C4_NA.md"
EXPECTED_MD5[NOTE_VR009] = "e25538129147e192c180af46b9c987bb"
# scoped C4 exception (locked note, 5 Nov 2026): image generation only; original C4 failure preserved; never a pass
C4_EXCEPTIONS = {("VR009", "head"): NOTE_VR009}
# ---- Amendment v1.1 (6 Nov 2026): full-cohort review, 100 patients + 20 hidden repeats ----
SPEC11 = "A3_visual_QC_spec_v1_1_amendment_20261106.md"
KEY100 = "A3_visualQC_KEY100_do_not_open_20261106.csv"
DRAW100 = "A3_visualQC_draw100_record_20261106.json"
SURVEY = "A3_overlap_survey_20261105_211213_CASE_IDS_DO_NOT_OPEN_until_unblinding.csv"
EXPECTED_MD5.update({SPEC11: "089e428056faeea5345a2f8602658166", KEY100: "aa72f6b91b2f7aa064eb6f2519033aa7",
                     DRAW100: "07bb3538feff0821b90124267be57528", SURVEY: "36853d5051537f19aae4bd6a382638cc"})
SEED100 = 20261106
# ---- second-reader plan (8 Nov 2026): 40 patients (26 critical + 14 random), codes SR01-SR40 ----
PLAN_SR = "A3_visualQC_second_reader_plan_20261108.md"
KEY_SR = "A3_visualQC_KEY_SR_do_not_open_20261108.csv"
DRAW_SR = "A3_visualQC_drawSR_record_20261108.json"
EXPECTED_MD5.update({PLAN_SR: "ff33c80f6b77dd758fee432a0f525f45", KEY_SR: "37404bc2316e9516a507fdd0df3311f4",
                     DRAW_SR: "a7d8259f74c1e8b0c2b2be5ebad7c10d"})
SEED_SR = 20261108
SR_CRITICAL = sorted(["train_0007", "train_0010", "train_0012", "train_0017", "train_0018", "train_0023", "train_0024",
                      "train_0028", "train_0033", "train_0042", "train_0052", "train_0053", "train_0056", "train_0066",
                      "train_0068", "train_0069", "train_0073", "train_0075", "train_0076", "train_0079", "train_0085",
                      "train_0087", "train_0089", "train_0091", "train_0095", "train_0100"])
N_REPEATS = 20
FIT = "Amendment3_fit_20261102_073040.csv"
TFM = "Amendment3_transforms_20261102_073040.csv"
A2P = "Audit1C_e_per_patient_20261031_141128.csv"
SIDES = "Amendment3b_femhead_20261102_211242/Amendment3b_femhead_sides_20261102_211242.csv"
KEY = "A3_visualQC_KEY_do_not_open_20261105.csv"
POOL = "A3_visualQC_control_pool_20261105.csv"
SEED = 20261105
TOL_REPRO_MM = a3.TOL_MM            # 0.01 mm (locked ICP tolerance)
TOL_COMPOSE_MM = 1e-4
BONE_W = (-200, 1300)
EDGE_HU = 300.0                     # display level for the FDG bone-edge contour (display only)
CHECKER_MM = 20.0
MARGIN_MM = 20.0
LAB = a3.LAB
REGION_LABELS = {"head": [LAB["skull"], LAB["brain"], LAB["C1"], LAB["C2"]],
                 "pelvis": [LAB["hip_L"], LAB["hip_R"], LAB["sacrum"], LAB["S1"], LAB["femur_L"], LAB["femur_R"]]}
SEL_LABEL = {"head": LAB["brain"], "pelvis": LAB["sacrum"]}   # label used for the image-level check


# ------------------------------ pure helpers (tested) ------------------------------
def checkerboard(a, b, spacing_rc, square_mm=CHECKER_MM):
    """2D checkerboard of images a and b (same shape); spacing_rc = (row mm, col mm)."""
    r = (np.arange(a.shape[0]) * spacing_rc[0] // square_mm).astype(int)
    c = (np.arange(a.shape[1]) * spacing_rc[1] // square_mm).astype(int)
    m = ((r[:, None] + c[None, :]) % 2).astype(bool)
    return np.where(m, b, a)


def xy_conditions(key_row, region):
    """returns ('A2' or 'A3a') for X and Y; checks the key is consistent."""
    x, y = key_row[f"{region}_X"], key_row[f"{region}_Y"]
    if {x, y} != {"A2", "A3a"}:
        raise RuntimeError(f"{key_row['code']} {region}: inconsistent X/Y assignment in key")
    return x, y


def redraw_controls(pool, targeted):
    """re-run the recorded draw (spec 3.2 / draw record). Returns controls, order, xy."""
    rng = np.random.default_rng(SEED)
    ctrl = sorted(rng.choice(sorted(pool), 10, replace=False).tolist())
    cases = sorted(targeted) + ctrl
    perm = rng.permutation(30)
    xy = rng.integers(0, 2, size=(30, 2))
    return ctrl, cases, perm, xy


def redraw100(cases_sorted):
    """re-run the v1.1 draw (draw record). Returns expected primary rows and repeat rows."""
    rng = np.random.default_rng(SEED100)
    perm = rng.permutation(len(cases_sorted)); xy = rng.integers(0, 2, size=(len(cases_sorted), 2))
    prim = [(f"VR{i + 1:03d}", cases_sorted[j], "A2" if xy[i, 0] == 0 else "A3a", "A2" if xy[i, 1] == 0 else "A3a")
            for i, j in enumerate(perm)]
    units = [(p[0], reg) for p in prim for reg in ("head", "pelvis")]
    rsel = rng.choice(len(units), N_REPEATS, replace=False); rord = rng.permutation(N_REPEATS)
    rxy = rng.integers(0, 2, N_REPEATS)
    reps = [(f"VR{101 + k:03d}",) + units[rsel[rord[k]]] + ("A2" if rxy[k] == 0 else "A3a",) for k in range(N_REPEATS)]
    return prim, reps          # prim: (code, case, head_X, pelvis_X); reps: (code, src_code, region, X)


def check_key100(K, cases_sorted):
    """structural check of the 100-patient key; returns (problems list, redraw_reproduced)."""
    K = K.fillna("")
    prob = []
    P = K[K["repeat_of"] == ""]; Rp = K[K["repeat_of"] != ""]
    if len(K) != 100 + N_REPEATS or K["code"].nunique() != len(K):
        prob.append("key must hold 120 unique codes")
    if sorted(P["code"]) != [f"VR{i:03d}" for i in range(1, 101)] or sorted(P["case_id"]) != list(cases_sorted):
        prob.append("primary codes VR001-VR100 must cover all 100 patients once")
    if sorted(Rp["code"]) != [f"VR{i:03d}" for i in range(101, 101 + N_REPEATS)]:
        prob.append("repeat codes must be VR101-VR120")
    pc = P.set_index("code")
    for _, r in P.iterrows():
        for reg in ("head", "pelvis"):
            try:
                xy_conditions(r, reg)
            except RuntimeError as e:
                prob.append(str(e))
    for _, r in Rp.iterrows():
        if r["repeat_of"] not in pc.index or pc.loc[r["repeat_of"], "case_id"] != r["case_id"] \
                or r["repeat_region"] not in ("head", "pelvis"):
            prob.append(f"{r['code']}: inconsistent repeat")
            continue
        try:
            xy_conditions(r, r["repeat_region"])
        except RuntimeError as e:
            prob.append(str(e))
    if Rp[["repeat_of", "repeat_region"]].duplicated().any():
        prob.append("a patient-region is repeated more than once")
    prim, reps = redraw100(list(cases_sorted))
    ok = not prob and all(pc.loc[c, "case_id"] == cid and pc.loc[c, "head_X"] == hx and pc.loc[c, "pelvis_X"] == px
                          for c, cid, hx, px in prim)
    rk = Rp.set_index("code") if len(Rp) else None
    ok = ok and rk is not None and all(rk.loc[c, "repeat_of"] == s and rk.loc[c, "repeat_region"] == reg
                                       and rk.loc[c, f"{reg}_X"] == x for c, s, reg, x in reps)
    return prob, bool(ok)


def redraw_sr(cases_sorted):
    pool = sorted(set(cases_sorted) - set(SR_CRITICAL))
    rng = np.random.default_rng(SEED_SR)
    rand = sorted(rng.choice(pool, 14, replace=False).tolist())
    allc = SR_CRITICAL + rand
    perm = rng.permutation(40); xy = rng.integers(0, 2, size=(40, 2))
    return [(f"SR{i + 1:02d}", allc[j], "A2" if xy[i, 0] == 0 else "A3a", "A2" if xy[i, 1] == 0 else "A3a")
            for i, j in enumerate(perm)]


def check_key_sr(SK, cases_sorted):
    """structural check of the second-reader key; returns (problems, redraw_reproduced)."""
    SK = SK.fillna("")
    prob = []
    if len(SK) != 40 or SK["code"].nunique() != 40 or SK["case_id"].nunique() != 40:
        prob.append("second-reader key must hold 40 unique codes and patients")
    if not set(SR_CRITICAL) <= set(SK["case_id"]):
        prob.append("critical stratum incomplete")
    if not set(SK["case_id"]) <= set(cases_sorted):
        prob.append("unknown case in key")
    for _, r in SK.iterrows():
        for reg in ("head", "pelvis"):
            try:
                xy_conditions(r, reg)
            except RuntimeError as e:
                prob.append(str(e))
    kk = SK.set_index("code")
    ok = not prob and all(c in kk.index and kk.loc[c, "case_id"] == cid and kk.loc[c, "head_X"] == hx
                          and kk.loc[c, "pelvis_X"] == px for c, cid, hx, px in redraw_sr(list(cases_sorted)))
    return prob, bool(ok)


def truncated_case(survey, oldkey):
    """the single patient with severe head truncation: brain overlap < 10% under A2 AND 3a in the locked survey, and
    identical to VR009 of the superseded 30-patient key. Returns case ID or None (-> STOP)."""
    s = survey[(survey["brain_overlap_A2"] < 0.1) & (survey["brain_overlap_A3a"] < 0.1)]
    if len(s) != 1:
        return None
    c = s["case_id"].iloc[0]
    old = oldkey.set_index("code").loc["VR009", "case_id"] if "VR009" in set(oldkey["code"]) else None
    return c if c == old else None


def within(a, b, tol):
    a, b = float(a), float(b)
    if np.isnan(a) and np.isnan(b):
        return True
    return bool(abs(a - b) <= tol)


def c4_gate(V, exceptions, md5_ok):
    """V: verification rows (dicts with code, region, C3_pass, C4_pass, compose_err_mm). Returns (stop_codes, allowed)
    where allowed = {(code, region): note} for C4 failures covered by a documented exception. An exception applies
    only if: listed, its note MD5 verified, C3 passed and composition passed. Any other failure stops."""
    stop, allowed = set(), {}
    for d in V:
        k = (d["code"], d["region"])
        if d["C3_pass"] and d["C4_pass"]:
            continue
        note = exceptions.get(k)
        if (note is not None and md5_ok.get(note, False) and d["C3_pass"]
                and d["compose_err_mm"] <= TOL_COMPOSE_MM and not d["C4_pass"]):
            allowed[k] = note
        else:
            stop.add(d["code"])
    return sorted(stop), allowed


def check_md5(paths):
    res = {}
    for name, exp in paths.items():
        p = (HERE / name) if name.endswith(".py") else (OUT_DIR / name)
        got = a3.md5(p) if p.exists() else "MISSING"
        res[name] = (got, got == exp)
    return res


# ------------------------------ transforms (reproduce the locked script) ------------------------------
def patient_transforms(cid, PP, TF):
    import SimpleITK as sitk
    dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
    T = sitk.ReadTransform(str(dF / "rigid.tfm"))
    A, b = a3.affine_of(T)
    t = np.zeros(3) if cid in a3.FALLBACKS else PP.loc[cid, ["t_x_LPS", "t_y_LPS", "t_z_LPS"]].to_numpy(float)
    out = dict(dP=dP, dF=dF, A=A, b=b, t=t, T=T)
    for region in ("head", "pelvis"):
        r = TF[(TF["case_id"] == cid) & (TF["region"] == region)]
        if len(r) != 1:
            raise RuntimeError(f"{cid} {region}: candidate transform missing in {TFM}")
        r = r.iloc[0]
        R = np.array([[r[f"R{i}{j}"] for j in range(3)] for i in range(3)], float)
        c = np.array([r["c0"], r["c1"], r["c2"]], float)
        out[region] = (R, c)
    return out


def resampling_transforms(dP, t, R, c, A, b, fallback):
    """exactly as Amendment3_regional_rigid.render_overlay: RA (A2) and RB (A3a candidate), with composition check."""
    import SimpleITK as sitk
    Rr = sitk.ReadTransform(str(dP / "rigid.tfm"))
    RA = sitk.CompositeTransform(3)
    RA.AddTransform(Rr)
    if not fallback:
        RA.AddTransform(sitk.TranslationTransform(3, tuple(float(v) for v in t)))
    Ri = R.T
    inv = sitk.AffineTransform(3)
    inv.SetMatrix(tuple(Ri.ravel()))
    inv.SetTranslation(tuple(float(v) for v in (-Ri @ c)))
    RB = sitk.CompositeTransform(3)
    RB.AddTransform(RA)
    RB.AddTransform(inv)
    worst = 0.0
    for p in np.random.default_rng(3).uniform(-200, 200, (5, 3)):
        qa = (A @ p + b) - t
        qb = R @ qa + c
        worst = max(worst, np.linalg.norm(np.array(RA.TransformPoint(tuple(qa))) - p),
                    np.linalg.norm(np.array(RB.TransformPoint(tuple(qb))) - p))
    return RA, RB, float(worst)


# ------------------------------ C3: reproduce stored selection metrics ------------------------------
def reproduce_selection(cid, X, fitrows):
    """recompute fov_zlo/zhi, sel_centroid_A2/A3 with the locked functions; compare with the fit table."""
    import SimpleITK as sitk
    imP = sitk.ReadImage(str(X["dP"] / "totseg_24.nii.gz")); imF = sitk.ReadImage(str(X["dF"] / "totseg_24.nii.gz"))
    gP, gF = a3.Geom(imP), a3.Geom(imF)
    aP, aF = sitk.GetArrayFromImage(imP), sitk.GetArrayFromImage(imF)
    A, b, t = X["A"], X["b"], X["t"]
    TA = lambda Q: (Q @ A.T + b) - t if len(Q) else Q
    zP = gP.z_range(); zF = TA(gF.corners())[:, 2]
    zlo, zhi = max(zP[0], zF.min()) + a3.FOV_MARGIN, min(zP[1], zF.max()) - a3.FOV_MARGIN

    def vox(labels):
        _, vP = a3.label_sets(aP, labels, gP)
        _, vF = a3.label_sets(aF, labels, gF)
        return a3.in_fov(vP, zlo, zhi), a3.in_fov(TA(vF), zlo, zhi)

    rows = []
    for region in ("head", "pelvis"):
        fr = fitrows[region]
        R, c = X[region]
        labs = [["brain"]] if region == "head" else [["sacrum"], ["S1"]]
        cA, cR = [], []
        for l in labs:
            vP, vF = vox([LAB[x] for x in l])
            if len(vP) and len(vF):
                cA.append(a3.centroid_dist(vP, vF)); cR.append(a3.centroid_dist(vP, vF, R, c))
        mA = float(np.nanmean(cA)) if cA else np.nan
        mR = float(np.nanmean(cR)) if cR else np.nan
        d = dict(region=region,
                 fov_zlo_diff=abs(zlo - fr["fov_zlo"]), fov_zhi_diff=abs(zhi - fr["fov_zhi"]),
                 sel_A2_stored=fr["sel_centroid_A2"], sel_A2_repro=mA, sel_A2_diff=abs(mA - fr["sel_centroid_A2"]),
                 sel_A3_stored=fr["sel_centroid_A3"], sel_A3_repro=mR, sel_A3_diff=abs(mR - fr["sel_centroid_A3"]))
        d["C3_pass"] = bool(within(zlo, fr["fov_zlo"], TOL_REPRO_MM) and within(zhi, fr["fov_zhi"], TOL_REPRO_MM)
                            and within(mA, fr["sel_centroid_A2"], TOL_REPRO_MM)
                            and within(mR, fr["sel_centroid_A3"], TOL_REPRO_MM))
        rows.append(d)
    return rows, imP, imF


# ------------------------------ C4: image-level check (as in render_overlay) ------------------------------
def image_level_check(imP, imF, RA, RB, A, b, t, R, c, lab):
    import SimpleITK as sitk
    tsP = sitk.Cast(imP, sitk.sitkUInt8); tsF = sitk.Cast(imF, sitk.sitkUInt8)

    def phys(im, kji):
        D = np.array(im.GetDirection()).reshape(3, 3); S = np.array(im.GetSpacing())
        return np.array(im.GetOrigin()) + (kji[:, ::-1] * S) @ D.T

    def inside(im, q):
        D = np.array(im.GetDirection()).reshape(3, 3); S = np.array(im.GetSpacing())
        ci = np.linalg.solve(D * S, (q - np.array(im.GetOrigin())).T).T
        return np.all((ci > 0.5) & (ci < np.array(im.GetSize()) - 1.5), axis=1)

    aF_ = sitk.GetArrayFromImage(tsF)
    idx = np.argwhere(aF_ == lab)
    tol = max(tsP.GetSpacing())
    if len(idx) < 50:
        return np.nan, tol, "label absent"
    qA = (phys(tsF, idx.astype(float)) @ A.T + b) - t
    errs = []
    for q, T_ in ((qA, RA), (qA @ R.T + c, RB)):
        k = inside(tsP, q)
        ra = sitk.GetArrayFromImage(sitk.Resample(tsF, tsP, T_, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8))
        ridx = np.argwhere(ra == lab)
        if k.sum() < 50 or len(ridx) < 50:
            errs.append(np.nan); continue
        errs.append(float(np.linalg.norm(phys(tsP, ridx.astype(float)).mean(0) - q[k].mean(0))))
    worst = float(np.nanmax(errs)) if np.isfinite(errs).any() else np.nan
    return worst, float(tol), "ok"


# ------------------------------ diagnostic (read-only; added 5 Nov 2026 for the VR009 C4 failure) ------------------------------
INTERIOR_MM = 10.0       # diagnostic only: 'interior' = at least this far from every face of the PSMA grid


def population_compare(q_ci, r_ci, size, q_mm, r_mm, border_vox):
    """q_ci: point-mapped FDG voxels as continuous PSMA indices (x, y, z); r_ci: resampled label voxels as PSMA indices
    (x, y, z); size: PSMA grid size (x, y, z); *_mm: the same points in mm. Returns centroid errors for three
    population definitions and the population sizes. Pure numpy (tested)."""
    size = np.asarray(size, float)
    def m_orig(ci):   # the C4 definition (unchanged)
        return np.all((ci > 0.5) & (ci < size - 1.5), axis=1)
    def m_nomargin(ci):
        return np.all((ci >= -0.5) & (ci <= size - 0.5), axis=1)
    def m_interior(ci):
        return np.all((ci >= border_vox) & (ci <= size - 1 - border_vox), axis=1)
    out = dict(n_mapped=int(len(q_ci)), frac_mapped_outside_grid=float(1 - m_nomargin(q_ci).mean()) if len(q_ci) else np.nan,
               n_resampled=int(len(r_ci)))
    for nm, fq, fr in (("C4_original", m_orig, None), ("no_margin", m_nomargin, None), ("interior", m_interior, m_interior)):
        kq = fq(q_ci)
        kr = fr(r_ci) if fr is not None else np.ones(len(r_ci), bool)
        if kq.sum() < 50 or kr.sum() < 50:
            out[f"err_{nm}_mm"] = np.nan; out[f"dxyz_{nm}_mm"] = (np.nan,) * 3
        else:
            d = r_mm[kr].mean(0) - q_mm[kq].mean(0)
            out[f"err_{nm}_mm"] = float(np.linalg.norm(d)); out[f"dxyz_{nm}_mm"] = tuple(round(float(v), 2) for v in d)
        out[f"n_mapped_{nm}"] = int(kq.sum()); out[f"n_resampled_{nm}"] = int(kr.sum())
    return out


def diagnose(code):
    """VR009-type C4 investigation: no images, no files written, no case ID printed; X/Y mapping withheld."""
    import SimpleITK as sitk
    md = check_md5(EXPECTED_MD5)
    bad = {k: v[0] for k, v in md.items() if not v[1]}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    K = pd.read_csv(OUT_DIR / KEY).set_index("code")
    if code not in K.index:
        sys.exit(f"{code} not in key")
    cid = K.loc[code, "case_id"]
    PP = pd.read_csv(OUT_DIR / A2P).set_index("case_id"); TF = pd.read_csv(OUT_DIR / TFM)
    X = patient_transforms(cid, PP, TF)
    imP = sitk.ReadImage(str(X["dP"] / "totseg_24.nii.gz")); imF = sitk.ReadImage(str(X["dF"] / "totseg_24.nii.gz"))
    tsP = sitk.Cast(imP, sitk.sitkUInt8); tsF = sitk.Cast(imF, sitk.sitkUInt8)
    D = np.array(tsP.GetDirection()).reshape(3, 3); S = np.array(tsP.GetSpacing()); O = np.array(tsP.GetOrigin())
    size = np.array(tsP.GetSize())
    to_ci = lambda q: np.linalg.solve(D * S, (q - O).T).T
    def phys(im, kji):
        Di = np.array(im.GetDirection()).reshape(3, 3); Si = np.array(im.GetSpacing())
        return np.array(im.GetOrigin()) + (kji[:, ::-1] * Si) @ Di.T
    aF_ = sitk.GetArrayFromImage(tsF)
    border = np.ceil(INTERIOR_MM / S).astype(int)
    print(f"DIAGNOSE {code} (case ID withheld; X/Y mapping withheld to preserve blinding)")
    print(f"PSMA label grid: size {tuple(int(v) for v in size)}, spacing {tuple(round(float(v), 3) for v in S)} mm; "
          f"C4 tolerance {max(S):.2f} mm; interior = >= {INTERIOR_MM:.0f} mm from every face")
    for region in ("head", "pelvis"):
        R, c = X[region]
        RA, RB, comp = resampling_transforms(X["dP"], X["t"], R, c, X["A"], X["b"], cid in a3.FALLBACKS)
        labs = [("brain", LAB["brain"]), ("skull", LAB["skull"])] if region == "head" else [("sacrum", LAB["sacrum"])]
        for lname, lab in labs:
            idx = np.argwhere(aF_ == lab)
            if len(idx) < 50:
                print(f"  {region} {lname}: label absent in FDG"); continue
            qA = (phys(tsF, idx.astype(float)) @ X["A"].T + X["b"]) - X["t"]
            for cond, q, T_ in (("A2", qA, RA), ("A3a", qA @ R.T + c, RB)):
                ra = sitk.GetArrayFromImage(sitk.Resample(tsF, tsP, T_, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8))
                ridx = np.argwhere(ra == lab).astype(float)
                r_ci = ridx[:, ::-1]
                o = population_compare(to_ci(q), r_ci, size, q, phys(tsP, ridx), border)
                ci = to_ci(q)
                faces = {f"{ax}{s}": float(((ci[:, k] < -0.5) if s == "-" else (ci[:, k] > size[k] - 0.5)).mean())
                         for k, ax in enumerate("xyz") for s in "-+"}
                faces = {k: round(v, 3) for k, v in faces.items() if v > 0}
                print(f"  {region:6s} {lname:6s} {cond:4s}: mapped outside grid {o['frac_mapped_outside_grid']:.3f} "
                      f"(by face {faces or 'none'}) | n mapped {o['n_mapped']}, n resampled {o['n_resampled']}")
                for nm in ("C4_original", "no_margin", "interior"):
                    print(f"      {nm:12s} err {o[f'err_{nm}_mm']:.2f} mm  dxyz {o[f'dxyz_{nm}_mm']}  "
                          f"(n mapped {o[f'n_mapped_{nm}']}, n resampled {o[f'n_resampled_{nm}']})")
    print("DIAGNOSE complete: nothing written.")


def overlap_fraction(q_ci, size):
    """fraction of points whose nearest PSMA voxel lies inside the grid (no-margin rule, as nearest-neighbour)."""
    if len(q_ci) == 0:
        return np.nan
    size = np.asarray(size, float)
    return float(np.all((q_ci >= -0.5) & (q_ci <= size - 0.5), axis=1).mean())


OVERLAP_BINS = [(0.9, 1.01, ">= 90%"), (0.5, 0.9, "50-90%"), (0.1, 0.5, "10-50%"), (-0.01, 0.1, "< 10%")]


def overlap_survey():
    """Descriptive field-of-view overlap audit (author decision 5 Nov 2026): for all 100 patients, the fraction of
    mapped FDG skull and brain voxels lying within the PSMA CT grid, under Amendment 2 and the 3a candidate.
    Measures overlap, NOT registration accuracy. No threshold, no exclusion, no images. Console: cohort-level
    distributions only. Per-patient table written to a file marked DO_NOT_OPEN until unblinding."""
    import SimpleITK as sitk
    t0 = time.time()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    md = check_md5(EXPECTED_MD5)
    bad = {k: v[0] for k, v in md.items() if not v[1]}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    PP = pd.read_csv(OUT_DIR / A2P).set_index("case_id"); TF = pd.read_csv(OUT_DIR / TFM)
    cases = sorted(TF["case_id"].unique())
    if len(cases) != 100:
        sys.exit(f"expected 100 patients in {TFM}, found {len(cases)}")
    rows = []
    for i, cid in enumerate(cases, 1):
        X = patient_transforms(cid, PP, TF)
        imP = sitk.ReadImage(str(X["dP"] / "totseg_24.nii.gz")); imF = sitk.ReadImage(str(X["dF"] / "totseg_24.nii.gz"))
        D = np.array(imP.GetDirection()).reshape(3, 3); S = np.array(imP.GetSpacing()); O = np.array(imP.GetOrigin())
        size = np.array(imP.GetSize())
        to_ci = lambda q: np.linalg.solve(D * S, (q - O).T).T
        DF = np.array(imF.GetDirection()).reshape(3, 3); SF = np.array(imF.GetSpacing()); OF = np.array(imF.GetOrigin())
        aF_ = sitk.GetArrayFromImage(imF)
        R, c = X["head"]
        row = dict(case_id=cid)
        for lname in ("skull", "brain"):
            idx = np.argwhere(aF_ == LAB[lname])
            row[f"{lname}_FDG_voxels"] = len(idx)
            if len(idx) == 0:
                row[f"{lname}_overlap_A2"] = row[f"{lname}_overlap_A3a"] = np.nan; continue
            p = OF + (idx[:, ::-1].astype(float) * SF) @ DF.T
            qA = (p @ X["A"].T + X["b"]) - X["t"]
            row[f"{lname}_overlap_A2"] = overlap_fraction(to_ci(qA), size)
            row[f"{lname}_overlap_A3a"] = overlap_fraction(to_ci(qA @ R.T + c), size)
        rows.append(row)
        del aF_, imP, imF
        if i % 10 == 0:
            print(f"  {i}/100 ({time.time() - t0:.0f} s)")
    T = pd.DataFrame(rows)
    out = OUT_DIR / f"A3_overlap_survey_{run_id}_CASE_IDS_DO_NOT_OPEN_until_unblinding.csv"
    T.to_csv(out, index=False)
    print("FIELD-OF-VIEW OVERLAP (mapped FDG voxels inside the PSMA CT grid) - descriptive only, not accuracy, "
          "no threshold applied")
    for col in ("skull_overlap_A2", "skull_overlap_A3a", "brain_overlap_A2", "brain_overlap_A3a"):
        v = T[col].dropna()
        bins = {lab: int(((v >= lo) & (v < hi)).sum()) for lo, hi, lab in OVERLAP_BINS}
        print(f"  {col:18s} n {len(v):3d} | min {v.min():.3f} P5 {v.quantile(.05):.3f} P25 {v.quantile(.25):.3f} "
              f"median {v.median():.3f} | patients {bins}")
    man = dict(step="A3 field-of-view overlap survey (descriptive)", run=run_id, script=Path(__file__).name,
               script_md5=a3.md5(Path(__file__)), inputs={k: v[0] for k, v in md.items()},
               output=out.name, output_md5=a3.md5(out), runtime_s=round(time.time() - t0),
               software=dict(python=platform.python_version(), numpy=np.__version__, SimpleITK=sitk.Version_VersionString()))
    (OUT_DIR / f"A3_overlap_survey_manifest_{run_id}.json").write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"per-patient table (do not open until unblinding): {out.name}; {time.time() - t0:.0f} s")


# ------------------------------ slice positions (PSMA anatomy only) ------------------------------
PETROUS_HU = 1300.0          # dense-bone level used only to LOCATE a candidate petrous slice (display only)
PETROUS_LAT_MM = 20.0        # lateral offset from the skull midline for the petrous search
PETROUS_MIN_VOX = 20         # minimum dense voxels on EACH side; otherwise no candidate -> flag, no substitution


def petrous_level(lab_zyx, ct_zyx, sp_xyz, z_lo, z_hi):
    """axial index in [z_lo, z_hi] with the most dense skull-label bone on BOTH sides lateral to the skull midline
    (petrous temporal bones). Returns (index or None, score). Display aid only; never a decision."""
    sk = lab_zyx == LAB["skull"]
    if not sk.any():
        return None, 0
    xc = np.argwhere(sk)[:, 2].mean()
    off = (np.arange(lab_zyx.shape[2]) - xc) * sp_xyz[0]
    best, score = None, 0
    for z in range(max(z_lo, 0), min(z_hi, lab_zyx.shape[0] - 1) + 1):
        d = sk[z] & (ct_zyx[z] > PETROUS_HU)
        l = int(d[:, off > PETROUS_LAT_MM].sum()); r = int(d[:, off < -PETROUS_LAT_MM].sum())
        if min(l, r) > score:
            best, score = z, min(l, r)
    return (best, score) if score >= PETROUS_MIN_VOX else (None, score)


def slice_positions(region, aP_lps, imP_lps, side_rows, ct_lps=None):
    """returns list of (name, plane, index) in the LPS-oriented PSMA label grid, plus notes."""
    import SimpleITK as sitk
    notes = []
    def cidx(phys_pt):
        return np.array(imP_lps.TransformPhysicalPointToContinuousIndex(tuple(float(v) for v in phys_pt)))  # (x, y, z)
    def lab_idx(l):
        return np.argwhere(aP_lps == l)          # (z, y, x)
    out = []
    if region == "head":
        br, sk = lab_idx(LAB["brain"]), lab_idx(LAB["skull"])
        if len(br) == 0 or len(sk) == 0:
            return [], ["brain or skull not labelled in PSMA CT"]
        zc, yc, xc = sk.mean(0)
        out += [("axial skull base candidate (brain z P5)", "ax", int(round(np.percentile(br[:, 0], 5)))),
                ("axial skull base candidate (brain z P12)", "ax", int(round(np.percentile(br[:, 0], 12))))]
        if ct_lps is not None:
            zp, sc = petrous_level(aP_lps, ct_lps, np.array(imP_lps.GetSpacing()),
                                   int(np.percentile(br[:, 0], 0)), int(np.percentile(br[:, 0], 40)))
            if zp is not None:
                out.append(("axial petrous level (densest bilateral skull base)", "ax", int(zp)))
            else:
                notes.append("SKULL-BASE LANDMARK SLICE NOT FOUND automatically: request an additional view")
        out += [("axial upper vault (brain z P80)", "ax", int(round(np.percentile(br[:, 0], 80)))),
                ("coronal skull centroid", "cor", int(round(yc))),
                ("sagittal skull centroid", "sag", int(round(xc)))]
    else:
        heads = {}
        for _, r in side_rows.iterrows():
            if bool(r.get("fitted_P", False)) and np.isfinite(r["cP_x"]):
                ci = cidx([r["cP_x"], r["cP_y"], r["cP_z"]])
                heads[r["side"]] = (ci, float(r["radius_P"]))
        sp = np.array(imP_lps.GetSpacing())
        if heads:
            zroof = np.mean([ci[2] + 0.85 * rad / sp[2] for ci, rad in heads.values()])
            ycor = np.mean([ci[1] for ci, _ in heads.values()])
        else:
            hp = np.vstack([lab_idx(LAB["hip_L"]), lab_idx(LAB["hip_R"])])
            if len(hp) == 0:
                return [], ["hip bones not labelled in PSMA CT"]
            zroof, ycor = hp[:, 0].mean(), hp[:, 1].mean()
            notes.append("femoral-head centres unavailable: hip-bone centroid used")
        out.append(("axial acetabular roofs", "ax", int(round(zroof))))
        s1 = lab_idx(LAB["S1"])
        if len(s1):
            out.append(("axial S1 centroid", "ax", int(round(s1[:, 0].mean()))))
        else:
            notes.append("S1 not labelled: S1 axial slice omitted")
        out.append(("coronal femoral-head centres", "cor", int(round(ycor))))
        for s, hl in (("R", LAB["hip_R"]), ("L", LAB["hip_L"])):
            if s in heads:
                out.append((f"sagittal {'right' if s == 'R' else 'left'} femoral head", "sag", int(round(heads[s][0][0]))))
            else:
                h = lab_idx(hl)
                if len(h):
                    out.append((f"sagittal {'right' if s == 'R' else 'left'} hip (label centroid)", "sag", int(round(h[:, 2].mean()))))
                    notes.append(f"{s} femoral-head centre unavailable: hip-bone centroid used")
    shp = aP_lps.shape
    lim = {"ax": shp[0], "cor": shp[1], "sag": shp[2]}
    return [(n, p, int(np.clip(i, 0, lim[p] - 1))) for n, p, i in out], notes


def cut(vol, plane, i):
    """vol (z, y, x) LPS -> 2D image with superior up for cor/sag; returns img, (row mm idx, col axis)."""
    if plane == "ax":
        return vol[i]
    if plane == "cor":
        return np.flipud(vol[:, i, :])
    return np.flipud(vol[:, :, i])


def plane_spacing(sp_xyz, plane):
    return {"ax": (sp_xyz[1], sp_xyz[0]), "cor": (sp_xyz[2], sp_xyz[0]), "sag": (sp_xyz[2], sp_xyz[1])}[plane]


# ------------------------------ rendering (no labels drawn) ------------------------------
def render_page(title, rows, fn):
    """rows: list of (slice_name, plane, sp_rc, P2, {'X': F2, 'Y': F2}). Writes PNG; no mask is ever passed here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["PSMA CT (reference)", "X: FDG CT", "X: checkerboard", "X: FDG bone edges on PSMA",
            "Y: FDG CT", "Y: checkerboard", "Y: FDG bone edges on PSMA"]
    # row heights follow each slice's physical aspect (removes empty space; image resolution unchanged)
    asp = [max(0.35, min(2.5, (P2.shape[0] * s[0]) / (P2.shape[1] * s[1]))) for _, _, s, P2, _ in rows]
    fig, axs = plt.subplots(len(rows), 7, figsize=(28, 0.6 + sum(3.6 * a for a in asp)), squeeze=False,
                            gridspec_kw=dict(height_ratios=asp))
    for r, (name, plane, sp_rc, P2, F) in enumerate(rows):
        ext = (0, P2.shape[1] * sp_rc[1], 0, P2.shape[0] * sp_rc[0])      # same pairing as the locked overlays
        panels = [P2, F["X"], checkerboard(P2, F["X"], sp_rc), None, F["Y"], checkerboard(P2, F["Y"], sp_rc), None]
        for j in range(7):
            ax = axs[r, j]
            base = P2 if panels[j] is None else panels[j]
            ax.imshow(base, cmap="gray", vmin=BONE_W[0], vmax=BONE_W[1], extent=ext, origin="upper", aspect="equal")
            if panels[j] is None:
                src = F["X"] if j == 3 else F["Y"]
                if (src > EDGE_HU).any() and (src <= EDGE_HU).any():
                    ax.contour(src, levels=[EDGE_HU], colors=["cyan"], linewidths=0.6, extent=ext, origin="upper")
            ax.set_title(f"{cols[j]} | {name}" if j == 0 else cols[j], fontsize=7.5)
            ax.tick_params(labelsize=5)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    return fn


def build_region_rows(region, X, imP_lab, imF_lab, side_rows):
    """resample FDG CT into a PSMA sub-grid for A2 and A3a; return page rows {'A2': F2, 'A3a': F2} and notes."""
    import SimpleITK as sitk
    R, c = X[region]
    RA, RB, comp = resampling_transforms(X["dP"], X["t"], R, c, X["A"], X["b"], X["cid"] in a3.FALLBACKS)
    if comp > TOL_COMPOSE_MM:
        raise RuntimeError(f"{X['cid']} {region}: composition check failed ({comp:.2e} mm)")
    ctP = sitk.ReadImage(str(X["dP"] / "CT.nii.gz"))
    ctF = sitk.ReadImage(str(X["dF"] / "CT.nii.gz"))
    # PSMA reference sub-grid: bounding box of the region labels (+ margin), on the PSMA label grid
    labP = sitk.GetArrayFromImage(imP_lab)
    m = np.isin(labP, REGION_LABELS[region])
    if not m.any():
        return [], ["region labels absent in PSMA CT"]
    sp = np.array(imP_lab.GetSpacing())
    zz, yy, xx = np.nonzero(m)
    g = np.ceil(MARGIN_MM / sp).astype(int)
    lo = np.maximum([xx.min() - g[0], yy.min() - g[1], zz.min() - g[2]], 0)
    hi = np.minimum([xx.max() + g[0] + 1, yy.max() + g[1] + 1, zz.max() + g[2] + 1], imP_lab.GetSize())
    ref = sitk.RegionOfInterest(imP_lab, [int(v) for v in hi - lo], [int(v) for v in lo])
    refO = sitk.DICOMOrient(ref, "LPS")
    P = sitk.DICOMOrient(sitk.Resample(ctP, ref, sitk.Transform(), sitk.sitkLinear, -1000, sitk.sitkFloat32), "LPS")
    Fv = {nm: sitk.DICOMOrient(sitk.Resample(ctF, ref, T_, sitk.sitkLinear, -1000, sitk.sitkFloat32), "LPS")
          for nm, T_ in (("A2", RA), ("A3a", RB))}
    Pa = sitk.GetArrayFromImage(P); Fa = {k: sitk.GetArrayFromImage(v) for k, v in Fv.items()}
    Lo = sitk.GetArrayFromImage(sitk.DICOMOrient(ref, "LPS"))
    pos, notes = slice_positions(region, Lo, refO, side_rows, ct_lps=Pa)
    spo = np.array(refO.GetSpacing())
    rows = []
    for name, plane, i in pos:
        rows.append((name, plane, plane_spacing(spo, plane), cut(Pa, plane, i), {k: cut(v, plane, i) for k, v in Fa.items()}))
    return rows, notes


README = """A3 visual registration QC - review instructions (spec {spec}; locked)
Level 2 only: anatomical validity of CT-to-CT alignment. No PET, no lesion information.
Each page = one coded patient (VRxxx) x one region. X and Y are two registration conditions in random order.
Columns: PSMA CT (reference) | X: FDG CT resampled | X: checkerboard ({sq:.0f} mm squares) | X: FDG bone edges
({edge:.0f} HU contour, cyan) on PSMA CT | the same three for Y. Slice positions are identical for X and Y.
Grade each region and condition (RATING_TEMPLATE_A3_visualQC.csv). Grade X and Y INDEPENDENTLY, each on its own
absolute scale, BEFORE recording the preference:
  3  Excellent  - no discernible mismatch of corresponding anatomical landmarks at the available image resolution
  2  Acceptable - minor residual mismatch without meaningful anatomical misregistration
  1  Poor       - definite anatomical misalignment requiring further assessment
  0  Failure    - gross misregistration or anatomically incorrect correspondence
  NA Not assessable - insufficient anatomical coverage to judge alignment (NOT a registration failure; record the
     available landmarks and, for the head, which parts of the vault are missing). A region graded NA does not
     remove the patient from the review of other regions.
Pages stamped "C4 EXCEPTION (documented)" are graded like any other page; the stamp means the image-level check
did not validate the registration (see A3_visualQC_note_20261105_VR009_C4_NA.md).
Head: grade the CRANIUM only (vault outer/inner table; clivus, sella, petrous ridges; orbits).
  Skull-base slices: two candidates (brain z P5, P12) and, where found, an automatically located petrous level.
  Record whether the petrous temporal bones AND the clivus are shown (skullbase_landmarks_shown yes / no). If no, or
  if the page states that the landmark slice was not found, request an additional view; nothing is substituted.
  Record separately (not graded): C1/C2 aligned / offset / not covered; mandible-jaw position same / different;
  neck posture difference yes / no. Do not penalise arm position.
Pelvis: grade the pelvic ring (iliac wings/crests, sacroiliac joints, sacrum promontory and foramina, pubic
  symphysis, ischial tuberosities), the acetabula and the corresponding femoral-head positions (native / prosthetic
  noted). Do NOT penalise femoral shaft position, lower-limb rotation or arm position: femoral orientation can
  change independently of the pelvis.
Preference per region: X better / Y better / equivalent. Name the landmark(s) driving each grade.
Additional views: if a landmark is not adequately shown, request extra slices (record the landmark); they are
generated identically for X and Y.
Do NOT open A3_visualQC_KEY100_do_not_open_20261106.csv (or any file marked DO_NOT_OPEN) until grading of ALL
pages is complete. Grade every page on its own, without looking back at other pages.
Category questions are answered after unblinding.
"""
TEMPLATE_COLS = ["code", "region", "grade_X", "grade_Y", "preference", "landmarks_X", "landmarks_Y",
                 "C1C2_X", "C1C2_Y", "mandible_position", "neck_posture_diff", "skullbase_landmarks_shown",
                 "additional_views_requested",
                 "comment"]


# ------------------------------ self-test ------------------------------
def selftest():
    import SimpleITK as sitk
    res = []
    a = np.zeros((40, 60)); b = np.ones((40, 60))
    cb = checkerboard(a, b, (1.0, 1.0), 20)
    res.append(("checkerboard alternates 20 mm squares from both images",
                cb[0, 0] == 0 and cb[0, 25] == 1 and cb[25, 0] == 1 and cb[25, 25] == 0 and cb.mean() == 0.5))
    res.append(("key X/Y consistency enforced", xy_conditions({"code": "VR001", "head_X": "A3a", "head_Y": "A2"}, "head")
                == ("A3a", "A2")))
    try:
        xy_conditions({"code": "VR001", "head_X": "A2", "head_Y": "A2"}, "head"); bad = False
    except RuntimeError:
        bad = True
    res.append(("inconsistent X/Y assignment rejected", bad))
    pool = [f"train_{i:04d}" for i in range(1, 74)]; T = [f"train_{i:04d}" for i in range(200, 220)]
    r1, r2 = redraw_controls(pool, T), redraw_controls(pool, T)
    res.append(("control draw deterministic (10 distinct, outside targeted)", r1[0] == r2[0] and len(set(r1[0])) == 10
                and not set(r1[0]) & set(T)))
    res.append(("reproduction tolerance: 0.005 mm passes, 0.02 mm fails",
                within(1.0, 1.005, TOL_REPRO_MM) and not within(1.0, 1.02, TOL_REPRO_MM)))
    # synthetic transform chain: inverse rigid pair + translation + candidate; composition and resampling
    Tf = sitk.Euler3DTransform(); Tf.SetRotation(0.05, -0.03, 0.08); Tf.SetTranslation((12.0, -7.0, 30.0))
    Tp = Tf.GetInverse()
    import tempfile
    _td = tempfile.TemporaryDirectory()                     # nothing is written to any study folder
    tmp = Path(_td.name)
    try:
        sitk.WriteTransform(Tp, str(tmp / "rigid.tfm"))
        A, bb = a3.affine_of(Tf)
        t = np.array([3.0, -2.0, 5.0]); R = a3.rot_xyz(2, -1, 3); c = np.array([1.5, 0.5, -2.0])
        RA, RB, comp = resampling_transforms(tmp, t, R, c, A, bb, False)
        res.append(("composition RA(T_A(p)) = p and RB(dT(T_A(p))) = p (< 1e-4 mm)", comp < TOL_COMPOSE_MM))
        # synthetic FDG volume with a bright block; resampled with RB it must land where dT(T_A) maps it
        F = sitk.Image(80, 80, 80, sitk.sitkFloat32); F.SetSpacing((2.0, 2.0, 2.0)); F.SetOrigin((-80.0, -80.0, -80.0))
        Fa = np.full((80, 80, 80), -1000, np.float32); Fa[30:45, 35:50, 25:40] = 1000
        Fi = sitk.GetImageFromArray(Fa); Fi.CopyInformation(F)
        Pgrid = sitk.Image(100, 100, 100, sitk.sitkFloat32); Pgrid.SetSpacing((2.0, 2.0, 2.0)); Pgrid.SetOrigin((-100.0, -100.0, -60.0))
        rs = sitk.GetArrayFromImage(sitk.Resample(Fi, Pgrid, RB, sitk.sitkLinear, -1000, sitk.sitkFloat32))
        kji = np.argwhere(rs > 0); got = np.array(Pgrid.GetOrigin()) + kji[:, ::-1].mean(0) * 2.0
        src = np.argwhere(Fa > 0); p = np.array(F.GetOrigin()) + src[:, ::-1].mean(0) * 2.0
        exp = R @ ((A @ p + bb) - t) + c
        res.append(("resampled FDG block lands at dT(T_A(p)) (< 1 voxel)", np.linalg.norm(got - exp) < 2.0))
    finally:
        _td.cleanup()
    # slice positions identical for X and Y; no labels reach the renderer; no case ID in titles
    vol = np.random.default_rng(0).normal(0, 300, (30, 40, 50)).astype(np.float32)
    rows = [("t", p, (1.0, 1.0), cut(vol, p, 10), {"X": cut(vol + 5, p, 10), "Y": cut(vol - 5, p, 10)})
            for p in ("ax", "cor", "sag")]
    res.append(("X and Y cut at identical positions (same shape per plane)",
                all(r[3].shape == r[4]["X"].shape == r[4]["Y"].shape for r in rows)))
    import inspect
    sig = inspect.signature(render_page).parameters
    res.append(("renderer accepts no label/mask input", list(sig) == ["title", "rows", "fn"]))
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        f = render_page("VR001 | head | X vs Y", rows, Path(td) / "t.png")
        res.append(("page renders (7 columns), file written", f.exists() and f.stat().st_size > 10000))
    res.append(("blinded title format contains no case ID", "train_" not in "VR001 | head | X vs Y"))
    # petrous finder: dense bilateral bone at z = 12 is found; unilateral-only volume gives no candidate (flag)
    lab = np.zeros((30, 40, 60), np.int16); ctv = np.zeros((30, 40, 60), np.float32)
    lab[:, 5:35, :] = LAB["skull"]
    ctv[12, 15:25, 2:8] = 2000; ctv[12, 15:25, 52:58] = 2000; ctv[20, 15:25, 2:8] = 2000
    zp, _ = petrous_level(lab, ctv, np.array([1.0, 1.0, 1.0]), 0, 29)
    ctu = np.zeros_like(ctv); ctu[12, 15:25, 2:8] = 2000
    zu, _ = petrous_level(lab, ctu, np.array([1.0, 1.0, 1.0]), 0, 29)
    res.append(("petrous finder: bilateral dense level found; unilateral only -> no candidate (flag, no substitution)",
                zp == 12 and zu is None))
    # diagnostic populations: an object truncated at the grid face. Point-mapped set extends beyond the grid;
    # resampled set is the in-grid part. Interior comparison must agree; no-margin comparison must agree; outside
    # fraction must be reported.
    g = np.stack(np.meshgrid(np.arange(-10, 30), np.arange(5, 25), np.arange(5, 25), indexing="ij"), -1).reshape(-1, 3).astype(float)
    size = np.array([40, 30, 30])
    r = g[(g[:, 0] >= 0)]
    o = population_compare(g, r, size, g, r, np.array([3, 3, 3]))
    res.append(("overlap fraction (no-margin rule): truncated object 0.75 inside",
                abs(overlap_fraction(g, size) - 0.75) < 1e-9))
    Vt = [dict(code="VR009", region="head", C3_pass=True, C4_pass=False, compose_err_mm=1e-12),
          dict(code="VR009", region="pelvis", C3_pass=True, C4_pass=True, compose_err_mm=1e-12),
          dict(code="VR010", region="head", C3_pass=True, C4_pass=True, compose_err_mm=1e-12)]
    ex = {("VR009", "head"): "note.md"}
    s1, a1 = c4_gate(Vt, ex, {"note.md": True})
    s2, _ = c4_gate(Vt, ex, {"note.md": False})
    s3, _ = c4_gate(Vt + [dict(code="VR011", region="head", C3_pass=True, C4_pass=False, compose_err_mm=1e-12)], ex,
                    {"note.md": True})
    s4, _ = c4_gate([dict(code="VR009", region="head", C3_pass=False, C4_pass=False, compose_err_mm=1e-12)], ex,
                    {"note.md": True})
    s5, _ = c4_gate([dict(code="VR009", region="pelvis", C3_pass=True, C4_pass=False, compose_err_mm=1e-12)], ex,
                    {"note.md": True})
    res.append(("C4 exception: only the listed code/region, only with verified note MD5, never with a C3 failure; "
                "any other failure stops", s1 == [] and a1 == {("VR009", "head"): "note.md"} and s2 == ["VR009"]
                and s3 == ["VR011"] and s4 == ["VR009"] and s5 == ["VR009"]))
    res.append(("diagnostic: truncated object -> outside fraction 0.25, no-margin and interior errors ~0",
                abs(o["frac_mapped_outside_grid"] - 0.25) < 1e-9 and o["err_no_margin_mm"] < 1e-9
                and o["err_interior_mm"] < 1e-9))
    selftest_v11(res)
    for n, ok in res:
        print(f"  [{'PASS' if ok else 'FAIL'}] {n}")
    allok = all(o for _, o in res)
    print(f"SELF-TEST {'PASSED' if allok else 'FAILED'} ({sum(o for _, o in res)}/{len(res)})")
    return allok


def selftest_v11(res):
    """Amendment v1.1 tests: key structure, redraw, repeats, truncated-case anchoring."""
    cases = [f"train_{i:04d}" for i in range(1, 101)]
    prim, reps = redraw100(cases)
    rows = [dict(code=c, case_id=cid, category="", head_X=hx, head_Y=("A3a" if hx == "A2" else "A2"), pelvis_X=px,
                 pelvis_Y=("A3a" if px == "A2" else "A2"), repeat_of="", repeat_region="") for c, cid, hx, px in prim]
    pc = {r["code"]: r["case_id"] for r in rows}
    for c, s, reg, x in reps:
        d = dict(code=c, case_id=pc[s], category="", head_X="", head_Y="", pelvis_X="", pelvis_Y="", repeat_of=s,
                 repeat_region=reg)
        d[f"{reg}_X"] = x; d[f"{reg}_Y"] = "A3a" if x == "A2" else "A2"
        rows.append(d)
    K = pd.DataFrame(rows)
    p1, ok1 = check_key100(K, cases)
    res.append(("v1.1 key: 100 primary + 20 repeats valid and draw reproduced", p1 == [] and ok1))
    K2 = K.copy(); K2.loc[K2.index[105], "case_id"] = "train_9999"
    p2, _ = check_key100(K2, cases)
    K3 = K.copy(); K3.loc[K3.index[3], "head_X"] = "A3a"; K3.loc[K3.index[3], "head_Y"] = "A3a"
    p3, _ = check_key100(K3, cases)
    res.append(("v1.1 key: inconsistent repeat and broken X/Y are detected", len(p2) > 0 and len(p3) > 0))
    res.append(("v1.1 repeats: 20 distinct patient-regions, none duplicated",
                len({(s, reg) for _, s, reg, _ in reps}) == N_REPEATS))
    sv = pd.DataFrame({"case_id": ["a", "b", "c"], "brain_overlap_A2": [1.0, 0.011, 0.95],
                       "brain_overlap_A3a": [1.0, 0.037, 0.99]})
    ok_old = pd.DataFrame({"code": ["VR009"], "case_id": ["b"]}); bad_old = pd.DataFrame({"code": ["VR009"], "case_id": ["c"]})
    sv2 = pd.concat([sv, pd.DataFrame({"case_id": ["d"], "brain_overlap_A2": [0.05], "brain_overlap_A3a": [0.05]})])
    res.append(("truncated case anchored only if survey and superseded key agree on exactly one patient",
                truncated_case(sv, ok_old) == "b" and truncated_case(sv, bad_old) is None
                and truncated_case(sv2, ok_old) is None))
    sr = redraw_sr(cases)
    SK = pd.DataFrame([dict(code=c, case_id=cid, stratum="", head_X=hx, head_Y=("A3a" if hx == "A2" else "A2"),
                            pelvis_X=px, pelvis_Y=("A3a" if px == "A2" else "A2")) for c, cid, hx, px in sr])
    q1, r1 = check_key_sr(SK, cases)
    SK2 = SK.copy(); SK2.loc[SK2.index[0], "head_X"] = SK2.loc[SK2.index[0], "head_Y"]
    q2, _ = check_key_sr(SK2, cases)
    res.append(("second-reader key: 40 patients incl. all 26 critical, draw reproduced; broken X/Y detected",
                q1 == [] and r1 and len(set(SR_CRITICAL) & set(SK.case_id)) == 26 and len(q2) > 0))


# ------------------------------ run ------------------------------
def run(preflight=False, second_reader=False):
    """Amendment v1.1: verification of all 100 patients, then 220 review pages (100 x head/pelvis + 20 hidden
    repeats). Preflight: the two patients already seen in the 5 Nov preflights (TEST codes, separate folder)."""
    import SimpleITK as sitk
    t0 = time.time()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    md = check_md5(EXPECTED_MD5)
    bad = {k: v[0] for k, v in md.items() if not v[1]}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(md)}/{len(md)} match")
    FITT = pd.read_csv(OUT_DIR / FIT); TF = pd.read_csv(OUT_DIR / TFM)
    PP = pd.read_csv(OUT_DIR / A2P).set_index("case_id")
    SD = pd.read_csv(OUT_DIR / SIDES)
    cases = sorted(FITT["case_id"].unique())
    oldK = pd.read_csv(OUT_DIR / KEY)
    K = pd.read_csv(OUT_DIR / KEY100).fillna("")
    prob, redraw_ok = check_key100(K, cases)
    if prob:
        sys.exit(f"STOP C2: key problems: {prob[:5]}")
    print(f"C2 key: 100 patients + {N_REPEATS} repeats OK; draw re-run in numpy {np.__version__}: "
          f"{'REPRODUCED' if redraw_ok else 'NOT reproduced in this numpy version (the locked key governs)'}")
    if preflight:
        others = [c for c in cases if c not in set(oldK["case_id"])][:2]          # same two as on 5 Nov
        work = pd.DataFrame({"code": ["TEST001", "TEST002"], "case_id": others, "head_X": ["A2", "A3a"],
                             "head_Y": ["A3a", "A2"], "pelvis_X": ["A2", "A3a"], "pelvis_Y": ["A3a", "A2"],
                             "repeat_of": ["", ""], "repeat_region": ["", ""]})
        tag, exceptions = "PREFLIGHT_", {}
    else:
        tc = truncated_case(pd.read_csv(OUT_DIR / SURVEY), oldK)
        if tc is None:
            sys.exit("STOP: the severely truncated patient could not be identified consistently from the overlap "
                     "survey and the superseded key (amendment v1.1 section 5)")
        work, tag, exceptions = K.sort_values("code"), "", {(tc, "head"): NOTE_VR009}
        if second_reader:
            SK = pd.read_csv(OUT_DIR / KEY_SR).fillna("")
            sprob, sredraw = check_key_sr(SK, cases)
            if sprob:
                sys.exit(f"STOP C2 (second reader): key problems: {sprob[:5]}")
            print(f"C2 second-reader key: 40 patients OK; draw re-run in numpy {np.__version__}: "
                  f"{'REPRODUCED' if sredraw else 'NOT reproduced in this numpy version (the locked key governs)'}")
            SK["repeat_of"] = ""; SK["repeat_region"] = ""
            work, tag = SK.sort_values("code"), "SECONDREADER_"
            redraw_ok = bool(redraw_ok and sredraw)
    pcode = {r["case_id"]: r["code"] for _, r in work.iterrows() if r["repeat_of"] == ""}   # case -> primary code
    odir = OUT_DIR / f"A3_visualQC_{tag}review{'40' if second_reader else '100'}_{run_id}"
    # C3 + C4 once per patient BEFORE any image
    ver, cache = [], {}
    for cid in sorted(pcode, key=lambda c: pcode[c]):
        X = patient_transforms(cid, PP, TF); X["cid"] = cid
        fr = {reg: FITT[(FITT["case_id"] == cid) & (FITT["region"] == reg)].iloc[0] for reg in ("head", "pelvis")}
        rows, imP, imF = reproduce_selection(cid, X, fr)
        for d in rows:
            reg = d["region"]
            R, c = X[reg]
            RA, RB, comp = resampling_transforms(X["dP"], X["t"], R, c, X["A"], X["b"], cid in a3.FALLBACKS)
            w, tol, st = image_level_check(imP, imF, RA, RB, X["A"], X["b"], X["t"], R, c, SEL_LABEL[reg])
            d.update(code=pcode[cid], case=cid, compose_err_mm=comp, image_check_worst_mm=w, image_check_tol_mm=tol,
                     image_check_status=st,
                     C4_pass=bool(comp <= TOL_COMPOSE_MM and (st != "ok" or (np.isfinite(w) and w <= tol))))
            ver.append(d)
        cache[cid] = X
        del imP, imF
        print(f"  verified {pcode[cid]}: C3 {'PASS' if all(d['C3_pass'] for d in rows) else 'FAIL'}, "
              f"C4 {'PASS' if all(d['C4_pass'] for d in ver[-2:]) else 'FAIL'}  ({time.time() - t0:.0f} s)")
    gate_in = [dict(d, code=d["case"]) for d in ver]                 # gate on case ID (exception is by case)
    stop_cases, allowed = c4_gate(gate_in, exceptions, {k: v[1] for k, v in md.items()})
    V = pd.DataFrame(ver)
    V["c4_exception"] = [f"DOCUMENTED EXCEPTION ({allowed[(c_, r_)]}, MD5 {md[allowed[(c_, r_)]][0]}); original C4 "
                         f"failure preserved" if (c_, r_) in allowed else "" for c_, r_ in zip(V["case"], V["region"])]
    V = V[["code", "region", "fov_zlo_diff", "fov_zhi_diff", "sel_A2_stored", "sel_A2_repro", "sel_A2_diff",
           "sel_A3_stored", "sel_A3_repro", "sel_A3_diff", "C3_pass", "compose_err_mm", "image_check_worst_mm",
           "image_check_tol_mm", "image_check_status", "C4_pass", "c4_exception"]].sort_values(["code", "region"])
    vp = OUT_DIR / f"A3_visualQC_{tag}verification100_{run_id}.csv"
    V.to_csv(vp, index=False)                                         # by code only; no case IDs
    print(f"C3 max |diff|: zlo {V.fov_zlo_diff.max():.2e}, zhi {V.fov_zhi_diff.max():.2e}, A2 {V.sel_A2_diff.max():.2e}, "
          f"A3 {V.sel_A3_diff.max():.2e} mm (tolerance {TOL_REPRO_MM} mm)")
    if stop_cases:
        sys.exit(f"STOP: verification failed for {sorted(pcode[c] for c in stop_cases)}; no image generated. "
                 f"See {vp.name}")
    allowed_codes = sorted(f"{pcode[c]} {r}" for c, r in allowed)
    if allowed:
        print(f"C4 documented exception(s) applied for image generation only: {allowed_codes}; original C4 failure "
              f"preserved in {vp.name}")
    print("C3 and C4: PASS (or documented exception) -> generating images")
    odir.mkdir()
    files, tmpl = [], []
    for cid in sorted(pcode, key=lambda c: pcode[c]):
        X = cache[cid]
        imP = sitk.ReadImage(str(X["dP"] / "totseg_24.nii.gz"))
        pages = [(r, reg) for _, r in work[work["case_id"] == cid].iterrows()
                 for reg in (("head", "pelvis") if r["repeat_of"] == "" else (r["repeat_region"],))]
        for reg in ("head", "pelvis"):
            todo = [r for r, g in pages if g == reg]
            if not todo:
                continue
            rows, notes = build_region_rows(reg, X, imP, None, SD[SD["case_id"] == cid])
            for kr in todo:
                code = kr["code"]
                if not rows:
                    print(f"  {code} {reg}: no image ({'; '.join(notes)})"); continue
                xc, yc = xy_conditions(kr, reg)
                prow = [(n, p, s, P2, {"X": F[xc], "Y": F[yc]}) for n, p, s, P2, F in rows]
                stamp = ("\nC4 EXCEPTION (documented): severe field-of-view truncation; image-level check failed under "
                         "the original definition; registration NOT validated by this check") if (cid, reg) in allowed else ""
                title = (f"{code} | {reg} | X vs Y (random order) | PSMA CT reference, FDG CT resampled | bone window "
                         f"{BONE_W[0]}..{BONE_W[1]} HU | no segmentation contours" + stamp
                         + (f"\nnote: {'; '.join(notes)}" if notes else ""))
                files.append(render_page(title, prow, odir / f"{code}_{reg}.png"))
                tmpl.append(dict(code=code, region=reg, **{k: "" for k in TEMPLATE_COLS[2:]}))
        del imP
        print(f"  {pcode[cid]}: images written ({time.time() - t0:.0f} s)")
    tp = odir / "RATING_TEMPLATE_A3_visualQC.csv"
    pd.DataFrame(sorted(tmpl, key=lambda d: (d["code"], d["region"])), columns=TEMPLATE_COLS).to_csv(tp, index=False)
    rp = odir / "README_review.txt"
    rp.write_text(README.format(spec=f"{SPEC} + {SPEC11}", sq=CHECKER_MM, edge=EDGE_HU), encoding="utf-8")
    man = dict(step=f"A3 visual registration QC images, full cohort ({'PREFLIGHT' if preflight else 'FULL'})",
               run=run_id, spec=[SPEC, SPEC11, NOTE_VR009], script=Path(__file__).name, script_md5=a3.md5(Path(__file__)),
               inputs={k: v[0] for k, v in md.items()}, draw_reproduced=bool(redraw_ok),
               c4_exceptions_applied=allowed_codes, verification=vp.name, verification_md5=a3.md5(vp),
               n_pages=len(files), outputs={p.name: a3.md5(p) for p in files + [tp, rp]},
               params=dict(TOL_REPRO_MM=TOL_REPRO_MM, TOL_COMPOSE_MM=TOL_COMPOSE_MM, BONE_W=BONE_W, EDGE_HU=EDGE_HU,
                           CHECKER_MM=CHECKER_MM, MARGIN_MM=MARGIN_MM, SEED100=SEED100, N_REPEATS=N_REPEATS),
               runtime_s=round(time.time() - t0),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             SimpleITK=sitk.Version_VersionString()), status="COMPLETE")
    (OUT_DIR / f"A3_visualQC_{tag}manifest100_{run_id}.json").write_text(json.dumps(man, indent=2, default=str),
                                                                         encoding="utf-8")
    print(f"done: {len(files)} pages in {odir.name}; verification {vp.name}; {time.time() - t0:.0f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--second-reader", action="store_true",
                    help="second-reader pages SR01-SR40 (plan 20261108): same images and checks, new codes and X/Y")
    ap.add_argument("--diagnose", metavar="VRxxx", help="read-only C4 diagnostic for one code (no images, no files)")
    ap.add_argument("--overlap-survey", action="store_true",
                    help="descriptive FOV overlap audit of skull/brain, all 100 patients (no images, no threshold)")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.overlap_survey:
        overlap_survey()
        return
    if a.diagnose:
        diagnose(a.diagnose)
        return
    run(preflight=a.preflight, second_reader=a.second_reader)


if __name__ == "__main__":
    main()
