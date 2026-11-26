"""
Audit1C_e_translation_refinement.py
===========================================================================
PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - AUDIT 1C-e
AMENDMENT 2: translation-only refinement of the supplied PSMA-FDG rigid
registration. Implements protocol v2.8, Section 15.2 (LOCKED text).

ONE-ATTEMPT RULE: this script is executed ONCE, after its MD5 has been recorded
in protocol Section 15.2.10. No parameter, structure set, criterion or code
may be changed after results are seen (any change = documented deviation).

Scope: anatomical CT-CT evaluation only. No tumour correspondence, no lesion
matching, no DMTVF, no discordance. TTB masks are used only to recompute the
field-of-view coverage measure (Section 8.1), never for matching.

Method (Section 15.2)
  d_s  = FDG/rigid.tfm(c_FDG,s) - c_PSMA,s      (LPS, mm; TotalSegmentator centroids)
  t    = component-wise median of d_s over the fitting structures (>= 5 required)
  T'(p) = FDG/rigid.tfm(p) - t                  point mapping FDG -> PSMA
  R'(q) = PSMA/rigid.tfm(q + t)                 resampling transform (sitk.Resample)
  Structure sets (fitting / heldout_bone / heldout_soft / head_neck) and their
  per-patient eligibility are taken VERBATIM from the Amendment 1 structure-level
  table (column evaluable_both == True), so they are identical to Amendment 1.

Mandatory verification BEFORE evaluation (every patient); any failure STOPS the
whole run (implementation error; never a fallback):
  (a) R'(T'(p)) = p and T'(R'(q)) = q for test points spanning both CT fields
  (b) T'(c_FDG,s) - c_PSMA,s = d_s - t for every fitting structure
  (c) composition order: T'(p) = T(p) - t and R'(q) = R(q + t) at the test points
  (d) the Amendment 1 change measure equals |t| at every fitting centroid (rigid R)
  tolerance: VERIFY_TOL_MM

Adoption criteria C1-C4 (Amendment 1, unchanged), computed BEFORE any
per-patient fallback: (i) over patients with a computed correction (decides),
(ii) over all patients (insufficient-landmarks patients enter with difference 0).
C4 is a visual review by the author of the overlays written here (CT only).

Implementation safeguards (author-approved corrections, 31 Oct 2026)
  - every Amendment 1 evaluable structure must be found again (else STOP)
  - supplied centroid distances must reproduce Amendment 1 within
    CENTROID_REPRO_TOL_MM for every structure (else STOP); Dice/MSD/HD95
    differences are reported and flagged against predefined thresholds
  - every patient, including insufficient-landmarks fallbacks, needs >= 1 valid
    held-out bone AND held-out soft measurement; otherwise C1-C3 are INVALID and
    the run is INCOMPLETE
  - Dice, MSD/HD95 (surface_metrics) and the C4 overlay() are copied VERBATIM from
    Amendment 1 (Audit1C_c_rigid_refinement.py, MD5 0c350ab57b59c1dea8fc729b434457b7)
  - pre-flight never prints or writes the estimated translation
  - the definitive run refuses to start until C4_PROCEDURE_VERIFIED = True

Inputs (read-only)
  D:\\DEEP PSMA\\DEEP PSMA EXTRACTED\\train_XXXX\\{PSMA,FDG}\\
      totseg_24.nii.gz, CT.nii.gz (overlays only), rigid.tfm,
      PET.nii.gz (PSMA, grid only), TTB.nii.gz (FDG, coverage only)
  06_output\\Audit1C_c_structure_level_20261029_212531.csv  (Amendment 1 sets)
Outputs: 06_output, single folder, prefix Audit1C_e_, run timestamp (Section 15.2.9).

PRE-FLIGHT (implementation check; NOT the experiment):
    python ...\\Audit1C_e_translation_refinement.py --preflight 3
  Processes the first N patients ONLY to: run the mandatory verification (a)-(c),
  and reproduce the SUPPLIED-transform metrics already published by Amendment 1
  and the supplied coverage of Audit 1C-b. No corrected-transform metric, no
  criterion, no overlay is computed or written, so no Amendment 2 result is seen.
  Outputs carry the prefix Audit1C_e_PREFLIGHT_. The flag is a command-line
  argument, so the script file (and its MD5) is identical for both modes.

DEFINITIVE RUN (once, after the MD5 is recorded in Section 15.2.10):
    conda activate totalseg
    python "D:\\DEEP PSMA\\PAPER3_Discordance\\02_code\\Audit1C_e_translation_refinement.py"
Requires: Python 3.9+, SimpleITK, numpy, pandas, matplotlib.
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
    import SimpleITK as sitk
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:  # mandatory outputs need all packages; fail before any processing
    sys.exit(f"Missing package: {e}")

# ------------------------------ configuration (LOCKED with the script MD5) ------------------------------
DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
A1_STRUCT = OUT_DIR / "Audit1C_c_structure_level_20261029_212531.csv"   # Amendment 1 structure sets
FILES = dict(totseg="totseg_24.nii.gz", ct="CT.nii.gz", tfm="rigid.tfm", pet="PET.nii.gz", ttb="TTB.nii.gz")

EXPECTED_PATIENTS = 100
MIN_FITTING = 5                    # >= 5 valid, anatomically distinct fitting structures (15.2.2)
VERIFY_TOL_MM = 1e-6               # inverse-pair / sign / composition verification (15.2.3)
N_TEST_GRID = 5                    # test points per axis in each CT field (5^3 per CT)
C1_DIFF_MAX = 0.0                  # median paired difference must be < 0 mm
C2_DIFF_MAX = 0.5                  # held-out soft median paired difference <= +0.5 mm
C3_WORSEN_MM = 2.0                 # patient worsened if held-out bone median centroid increases > 2 mm
C3_MAX_PATIENTS = 5                # ... in <= 5 patients
C4_N = 5                           # PRIMARY C4 set size (15.2.5)
# C4 selection (author decision 31 Oct 2026, Option B):
#  PRIMARY  = the C4_N patients with the largest refinement change, computed EXACTLY as Amendment 1's
#             refinement_change_median_mm: median over the patient's fitting-structure centroids (PSMA space)
#             of |R_supplied(q) - R_corrected(q)|, R = resampling transforms (Amendment 1 sorted(change, reverse)[:5]).
#  SAFETY EXTENSION (prospective; not part of the Amendment 1 replication) = the patient with the largest
#             held-out-bone deterioration and the patient with the largest improvement, if not already selected.
#  C4 is judged on ALL selected overlays; gross misregistration in any of them fails C4.
COVERAGE_MARGIN_MM = 10.0
VOL_RATIO_OK = (0.95, 1.05)            # Amendment 1 evaluability range; for the corrected transform REPORTED (sets are fixed)
VERT = list(range(25, 51))              # vertebral labels (Amendment 1 overlay centring)
WINDOW = (-200, 1000)                   # Amendment 1 overlay CT window
# C4 procedure must be confirmed against the Amendment 1 report (06_output\Audit1C_final_report_20261030.md)
# before the MD5 lock. The definitive run refuses to start while this is False (pre-flight is allowed).
C4_PROCEDURE_VERIFIED = True        # verified 31 Oct 2026 against Amendment 1 code (Option B; protocol 15.2.5)
# Reproduction of Amendment 1 SUPPLIED-transform values (Section 15.2; implementation check)
CENTROID_REPRO_TOL_MM = 0.01       # HARD STOP if any structure differs by more than this (C1-C3 use centroids)
REPRO_FLAG_DICE = 0.02             # flag only (not decisive for C1-C3): |Dice difference|
REPRO_FLAG_MSD_MM = 0.5            # flag only: |MSD difference| in mm
REPRO_FLAG_HD95_MM = 1.0           # flag only: |HD95 difference| in mm
# Evaluation completeness: EVERY patient (including insufficient-landmarks fallbacks) must have
# >= 1 valid measurement in each of these sets, otherwise the run is INCOMPLETE and C1-C3 are invalid.
REQUIRED_HELDOUT_SETS = ["heldout_bone", "heldout_soft"]
SETS = ["fitting", "heldout_bone", "heldout_soft", "head_neck"]

# ------------------------------ mode ------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--preflight", type=int, default=0,
                help="implementation check on the first N patients (supplied metrics + verification only)")
ARGS = ap.parse_args()
PREFLIGHT = ARGS.preflight > 0
PREFIX = "Audit1C_e_PREFLIGHT_" if PREFLIGHT else "Audit1C_e_"

# ------------------------------ helpers ------------------------------
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


def stop_run(reason: str):
    """Verification failure or fatal error: stop the whole experiment (Section 15.2.3)."""
    say(f"\n*** RUN STOPPED: {reason} ***")
    say("This is an implementation error, not a refinement result. Document as a protocol deviation before any re-run.")
    out("RUN_STOPPED", "txt").write_text("\n".join(LOG), encoding="utf-8")
    sys.exit(1)


def centroids(label_img):
    f = sitk.LabelShapeStatisticsImageFilter()
    f.Execute(label_img)
    return {int(k): np.array(f.GetCentroid(k)) for k in f.GetLabels()}, f


def corrected_transforms(T_point, R_resample, t):
    """T'(p) = T(p) - t (FDG->PSMA points); R'(q) = R(q + t) (resampling). ITK composite: last added applied first."""
    Tc = sitk.CompositeTransform(3)
    Tc.AddTransform(sitk.TranslationTransform(3, tuple(float(-x) for x in t)))
    Tc.AddTransform(T_point)                     # applied first
    Rc = sitk.CompositeTransform(3)
    Rc.AddTransform(R_resample)
    Rc.AddTransform(sitk.TranslationTransform(3, tuple(float(x) for x in t)))   # applied first
    return Tc, Rc


def grid_points(img, n):
    size = np.array(img.GetSize())
    axes = [np.linspace(0, s - 1, n) for s in size]
    pts = []
    for i in axes[0]:
        for j in axes[1]:
            for k in axes[2]:
                pts.append(np.array(img.TransformContinuousIndexToPhysicalPoint((float(i), float(j), float(k)))))
    return pts


def tp(T, p):
    return np.array(T.TransformPoint(tuple(float(x) for x in p)))


def surface_metrics(imgA, imgB, l, bbA, bbB, margin=3):
    """VERBATIM from Amendment 1 (02_code\\Audit1C_c_rigid_refinement.py, MD5 0c350ab57b59c1dea8fc729b434457b7):
    region of interest = union of the two label bounding boxes + 3-voxel margin (clipped to the image);
    pooled symmetric contour-to-boundary distances; MSD = mean, HD95 = 95th percentile."""
    sz = imgA.GetSize()
    lo = [int(max(0, min(bbA[i], bbB[i]) - margin)) for i in range(3)]
    hi = [int(min(sz[i], max(bbA[i] + bbA[i + 3], bbB[i] + bbB[i + 3]) + margin)) for i in range(3)]
    size = [hi[i] - lo[i] for i in range(3)]
    a = sitk.BinaryThreshold(sitk.RegionOfInterest(imgA, size, lo), l, l, 1, 0)
    b = sitk.BinaryThreshold(sitk.RegionOfInterest(imgB, size, lo), l, l, 1, 0)
    ca = sitk.GetArrayFromImage(sitk.BinaryContour(a, False, 0, 1)) > 0
    cb = sitk.GetArrayFromImage(sitk.BinaryContour(b, False, 0, 1)) > 0
    da = np.abs(sitk.GetArrayFromImage(sitk.SignedMaurerDistanceMap(a, insideIsPositive=False,
                                                                     squaredDistance=False, useImageSpacing=True)))
    db = np.abs(sitk.GetArrayFromImage(sitk.SignedMaurerDistanceMap(b, insideIsPositive=False,
                                                                     squaredDistance=False, useImageSpacing=True)))
    v = np.concatenate([db[ca], da[cb]])
    return (float(v.mean()), float(np.percentile(v, 95))) if v.size else (float("nan"), float("nan"))


def label_metrics(segP, res, ovf, shpP, shpR, shpF, lab):
    """Dice (LabelOverlapMeasuresImageFilter on the full grid, as Amendment 1), MSD/HD95 (Amendment 1
    surface_metrics), and the resampled/native volume ratio (Amendment 1 evaluability check, reported)."""
    if not shpR.HasLabel(lab):
        return np.nan, np.nan, np.nan, 0.0
    try:
        dice = float(ovf.GetDiceCoefficient(int(lab)))
    except Exception:
        dice = np.nan
    msd, hd = surface_metrics(segP, res, lab, shpP.GetBoundingBox(lab), shpR.GetBoundingBox(lab))
    vr = shpR.GetPhysicalSize(lab) / shpF.GetPhysicalSize(lab)
    return dice, msd, hd, vr


def fov_fraction(points, img, margin_mm=0.0):
    """Fraction of physical points inside the image field (voxel-edge bounds), optionally shrunk by margin_mm."""
    if len(points) == 0:
        return np.nan
    size = np.array(img.GetSize()); sp = np.array(img.GetSpacing())
    D = np.array(img.GetDirection()).reshape(3, 3); o = np.array(img.GetOrigin())
    ci = np.linalg.solve(D * sp, (np.asarray(points) - o).T).T          # continuous index (x, y, z)
    m = margin_mm / sp
    inside = np.all((ci >= -0.5 + m) & (ci <= size - 0.5 - m), axis=1)
    return float(inside.mean())


def affine_points(T, pts):
    """Apply a point transform to many points (exact; loop kept simple and transform-agnostic)."""
    return np.array([tp(T, p) for p in pts]) if len(pts) else np.zeros((0, 3))


# ------------------------------ start ------------------------------
script_path = Path(__file__).resolve()
SCRIPT_MD5 = md5(script_path)
say(f"Audit 1C-e (Amendment 2) - run {RUN}")
say(f"Script: {script_path} (MD5 {SCRIPT_MD5})")
say("Protocol v2.8 Section 15.2 (LOCKED). ONE execution only. No correspondence or discordance.")
if PREFLIGHT:
    say(f"*** PRE-FLIGHT MODE: first {ARGS.preflight} patients; verification + SUPPLIED metrics only; no Amendment 2 result is computed ***")
if not PREFLIGHT and not C4_PROCEDURE_VERIFIED:
    print("DEFINITIVE RUN REFUSED: C4_PROCEDURE_VERIFIED is False. Confirm the C4 selection procedure against "
          "Amendment 1, set the flag, record the script MD5 in protocol Section 15.2.10, then run. Nothing was processed.")
    sys.exit(2)
if not A1_STRUCT.is_file():
    stop_run(f"Amendment 1 structure table not found: {A1_STRUCT}")
A1 = pd.read_csv(A1_STRUCT)
A1 = A1[A1["evaluable_both"] == True]  # noqa: E712
say(f"Amendment 1 structure assignments: {len(A1)} evaluable structures, {A1['case_id'].nunique()} patients "
    f"(MD5 {md5(A1_STRUCT)})")

cases = sorted(p.name for p in DATA_ROOT.glob("train_*") if p.is_dir())
say(f"Case folders found: {len(cases)} (expected {EXPECTED_PATIENTS})")
n_found = len(cases)
if PREFLIGHT:
    cases = cases[:ARGS.preflight]

struct_rows, patient_rows, verify_rows, cover_rows, inventory = [], [], [], [], []

# ------------------------------ main loop ------------------------------
for ci_, cid in enumerate(cases, 1):
    say(f"[{ci_:3d}/{len(cases)}] {cid}")
    prow = dict(case_id=cid, status="", fallback_reason="", n_fitting=0, fitting_structures="",
                t_x_LPS=np.nan, t_y_LPS=np.nan, t_z_LPS=np.nan, t_magnitude_mm=np.nan, t_inplane_mm=np.nan)
    try:
        dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
        paths = {"PSMA_totseg": dP / FILES["totseg"], "FDG_totseg": dF / FILES["totseg"],
                 "PSMA_tfm": dP / FILES["tfm"], "FDG_tfm": dF / FILES["tfm"],
                 "PSMA_PET": dP / FILES["pet"], "FDG_TTB": dF / FILES["ttb"],
                 "PSMA_CT": dP / FILES["ct"], "FDG_CT": dF / FILES["ct"]}
        for role, p in paths.items():
            if not p.is_file():
                raise FileNotFoundError(f"{role}: {p}")
            st = p.stat()
            inventory.append(dict(case_id=cid, role=role, path=str(p), bytes=st.st_size,
                                  modified=datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
                                  md5=md5(p)))

        R = sitk.ReadTransform(str(paths["PSMA_tfm"]))   # resampling FDG -> PSMA grid (verified, Audit 1C-a)
        T = sitk.ReadTransform(str(paths["FDG_tfm"]))    # point mapping FDG -> PSMA (verified, Audit 1C-a)
        segP = sitk.ReadImage(str(paths["PSMA_totseg"]), sitk.sitkUInt8)
        segF = sitk.ReadImage(str(paths["FDG_totseg"]), sitk.sitkUInt8)
        cenP, shpP = centroids(segP)
        cenF, shpF = centroids(segF)

        A = A1[A1["case_id"] == cid]
        if len(A) == 0:
            stop_run(f"structure mismatch in {cid}: patient absent from the Amendment 1 structure table")
        missing = [f"{s} ({int(l)})" for l, s in zip(A["label_id"], A["structure"])
                   if int(l) not in cenP or int(l) not in cenF]
        if missing:
            stop_run(f"structure mismatch in {cid}: Amendment 1 evaluable structures not found now: {missing}")
        fit = A[A["set"] == "fitting"]
        fit = fit[fit["label_id"].isin(list(cenP.keys())) & fit["label_id"].isin(list(cenF.keys()))]
        dmat = np.array([tp(T, cenF[int(l)]) - cenP[int(l)] for l in fit["label_id"]]) if len(fit) else np.zeros((0, 3))
        prow["n_fitting"] = int(len(fit))
        prow["fitting_structures"] = ";".join(fit["structure"])

        if len(fit) >= MIN_FITTING:
            t = np.median(dmat, axis=0)
            prow.update(status="CORRECTION_COMPUTED", t_x_LPS=t[0], t_y_LPS=t[1], t_z_LPS=t[2],
                        t_magnitude_mm=float(np.linalg.norm(t)), t_inplane_mm=float(np.hypot(t[0], t[1])))
            Tc, Rc = corrected_transforms(T, R, t)

            # ---------- mandatory verification (stop on any failure) ----------
            ptsF = grid_points(segF, N_TEST_GRID); ptsP = grid_points(segP, N_TEST_GRID)
            e_a1 = max(np.linalg.norm(tp(Rc, tp(Tc, p)) - p) for p in ptsF)
            e_a2 = max(np.linalg.norm(tp(Tc, tp(Rc, q)) - q) for q in ptsP)
            e_b = max(np.linalg.norm((tp(Tc, cenF[int(l)]) - cenP[int(l)]) - (dmat[i] - t))
                      for i, l in enumerate(fit["label_id"]))
            e_c1 = max(np.linalg.norm(tp(Tc, p) - (tp(T, p) - t)) for p in ptsF)
            e_c2 = max(np.linalg.norm(tp(Rc, q) - tp(R, q + t)) for q in ptsP)
            # Amendment 1 change measure at the fitting-structure centroids (PSMA space), resampling transforms
            qfit = [cenP[int(l)] for l in fit["label_id"]]
            chg = np.array([np.linalg.norm(tp(R, q) - tp(Rc, q)) for q in qfit])
            prow["A1_change_median_mm"] = float(np.median(chg))
            # (d) interpretation check: for a rigid R, |R(q+t) - R(q)| = |t| at every point
            e_d = float(np.max(np.abs(chg - np.linalg.norm(t))))
            ok = all(e < VERIFY_TOL_MM for e in (e_a1, e_a2, e_b, e_c1, e_c2, e_d))
            verify_rows.append(dict(case_id=cid, roundtrip_FDG_max_mm=e_a1, roundtrip_PSMA_max_mm=e_a2,
                                    fitting_sign_check_max_mm=e_b, point_composition_max_mm=e_c1,
                                    resample_composition_max_mm=e_c2, change_equals_translation_max_mm=e_d,
                                    tolerance_mm=VERIFY_TOL_MM, passed=ok))
            if not ok:
                pd.DataFrame(verify_rows).to_csv(out("verification_log", "csv"), index=False)
                stop_run(f"verification failed in {cid}: roundtrip {e_a1:.2e}/{e_a2:.2e}, sign {e_b:.2e}, "
                         f"composition {e_c1:.2e}/{e_c2:.2e} mm")
        else:
            t = None
            prow["A1_change_median_mm"] = 0.0          # supplied retained: no change (as Amendment 1 fallbacks)
            prow.update(status="INSUFFICIENT_LANDMARKS_FALLBACK",
                        fallback_reason=f"{len(fit)} fitting structures (< {MIN_FITTING}); supplied retained")
            verify_rows.append(dict(case_id=cid, passed=np.nan, note="no correction computed"))

        # ---------- resample FDG labels onto the PSMA CT grid (nearest neighbour) ----------
        resS = sitk.Resample(segF, segP, R, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
        resC = sitk.Resample(segF, segP, Rc, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8) if (t is not None and not PREFLIGHT) else resS
        fS = sitk.LabelShapeStatisticsImageFilter(); fS.Execute(resS)
        oS = sitk.LabelOverlapMeasuresImageFilter(); oS.Execute(segP, resS)
        if t is not None and not PREFLIGHT:
            fC = sitk.LabelShapeStatisticsImageFilter(); fC.Execute(resC)
            oC = sitk.LabelOverlapMeasuresImageFilter(); oC.Execute(segP, resC)

        # ---------- structure-level evaluation (all four sets) ----------
        repro_fail = []
        for _, r in A.iterrows():
            lab = int(r["label_id"])
            if lab not in cenP or lab not in cenF:          # cannot happen after the check above
                stop_run(f"structure mismatch in {cid}: label {lab}")
            pS = tp(T, cenF[lab]); cS = float(np.linalg.norm(pS - cenP[lab]))
            a1c = r["supplied_centroid_mm"]
            if not np.isfinite(a1c) or abs(cS - a1c) > CENTROID_REPRO_TOL_MM:
                repro_fail.append(f"{r['structure']}: now {cS:.4f} vs Amendment 1 {a1c} mm")
            if t is not None and not PREFLIGHT:
                pC = tp(Tc, cenF[lab]); cC = float(np.linalg.norm(pC - cenP[lab]))
            else:
                cC = np.nan if PREFLIGHT else cS
            dS, mS, hS, vS = label_metrics(segP, resS, oS, shpP, fS, shpF, lab)
            if t is not None and not PREFLIGHT:
                dC, mC, hC, vC = label_metrics(segP, resC, oC, shpP, fC, shpF, lab)
            elif PREFLIGHT:
                dC, mC, hC, vC = np.nan, np.nan, np.nan, np.nan
            else:
                dC, mC, hC, vC = dS, mS, hS, vS
            struct_rows.append(dict(case_id=cid, label_id=lab, structure=r["structure"], set=r["set"],
                                    correction_status=prow["status"],
                                    supplied_centroid_mm=cS, corrected_centroid_mm=cC, diff_centroid_mm=cC - cS,
                                    supplied_dice=dS, corrected_dice=dC,
                                    supplied_MSD_mm=mS, corrected_MSD_mm=mC,
                                    supplied_HD95_mm=hS, corrected_HD95_mm=hC,
                                    supplied_volume_ratio=vS, corrected_volume_ratio=vC,
                                    corrected_volume_ratio_outside_A1_range=bool(np.isfinite(vC) and not (VOL_RATIO_OK[0] <= vC <= VOL_RATIO_OK[1])),
                                    A1_supplied_centroid_mm=r["supplied_centroid_mm"],
                                    A1_supplied_dice=r["supplied_dice"], A1_supplied_MSD_mm=r["supplied_MSD_mm"],
                                    A1_supplied_HD95_mm=r["supplied_HD95_mm"]))

        if repro_fail:
            stop_run(f"supplied centroid distances do not reproduce Amendment 1 in {cid} "
                     f"(tolerance {CENTROID_REPRO_TOL_MM} mm): {repro_fail[:5]}{' ...' if len(repro_fail) > 5 else ''}")

        # ---------- coverage: FDG TTB voxels inside the PSMA PET field of view ----------
        pet = sitk.ReadImage(str(paths["PSMA_PET"]))
        ttb = sitk.ReadImage(str(paths["FDG_TTB"]), sitk.sitkUInt8)
        idx = np.argwhere(sitk.GetArrayViewFromImage(ttb) > 0)[:, ::-1]   # (x, y, z)
        ptsT = [np.array(ttb.TransformIndexToPhysicalPoint(tuple(int(v) for v in i))) for i in idx]
        mapS = affine_points(T, ptsT)
        if PREFLIGHT:
            mapC = np.zeros((0, 3))
        else:
            mapC = affine_points(Tc, ptsT) if t is not None else mapS
        cover_rows.append(dict(case_id=cid, n_FDG_TTB_voxels=len(ptsT),
                               supplied_FDG_TTB_in_PSMA_PET_FOV=fov_fraction(mapS, pet),
                               corrected_FDG_TTB_in_PSMA_PET_FOV=fov_fraction(mapC, pet),
                               supplied_FDG_TTB_in_PSMA_PET_FOV_10mm=fov_fraction(mapS, pet, COVERAGE_MARGIN_MM),
                               corrected_FDG_TTB_in_PSMA_PET_FOV_10mm=fov_fraction(mapC, pet, COVERAGE_MARGIN_MM)))
        if PREFLIGHT:   # the estimated translation is never printed or written in pre-flight
            say(f"   {prow['status']} | fitting {prow['n_fitting']} | verification "
                + ("n/a" if t is None else "PASSED") + " | supplied centroids reproduce Amendment 1")
        else:
            say(f"   {prow['status']} | fitting {prow['n_fitting']} | t = "
                + ("n/a" if t is None else f"({t[0]:.1f}, {t[1]:.1f}, {t[2]:.1f}) mm, |t| {np.linalg.norm(t):.1f}"))
    except SystemExit:
        raise
    except Exception as e:  # processing failure (not a verification failure): recorded, run becomes INCOMPLETE
        prow.update(status="FAILED", fallback_reason=f"processing error: {e}")
        say(f"   FAILED: {e}\n{traceback.format_exc(limit=2)}")
    patient_rows.append(prow)

# ------------------------------ completeness and reproduction helpers ------------------------------
def evaluated_counts(SLx, cols):
    """Number of structures per patient and set with all `cols` finite."""
    if len(SLx) == 0:
        return pd.DataFrame(columns=["case_id"])
    ok = SLx[cols].notna().all(axis=1) & np.isfinite(SLx[cols]).all(axis=1)
    cnt = SLx[ok].groupby(["case_id", "set"]).size().unstack(fill_value=0)
    cnt.columns = [f"n_{c}_evaluated" for c in cnt.columns]
    return cnt.reset_index()


def completeness(counts, case_list):
    """Patients lacking >= 1 valid measurement in any REQUIRED_HELDOUT_SETS."""
    c = pd.DataFrame({"case_id": case_list}).merge(counts, on="case_id", how="left").fillna(0)
    bad = []
    for s in REQUIRED_HELDOUT_SETS:
        col = f"n_{s}_evaluated"
        vals = c[col] if col in c.columns else pd.Series(0, index=c.index)
        bad += c.loc[vals < 1, "case_id"].tolist()
    return sorted(set(bad)), c


def repro_flags(SLx):
    """Supplied-transform reproduction of Amendment 1: centroid (stop already enforced) + flagged metrics."""
    out_ = {}
    for m, a1, thr in [("supplied_centroid_mm", "A1_supplied_centroid_mm", CENTROID_REPRO_TOL_MM),
                       ("supplied_dice", "A1_supplied_dice", REPRO_FLAG_DICE),
                       ("supplied_MSD_mm", "A1_supplied_MSD_mm", REPRO_FLAG_MSD_MM),
                       ("supplied_HD95_mm", "A1_supplied_HD95_mm", REPRO_FLAG_HD95_MM)]:
        dd = (SLx[m] - SLx[a1]).abs()
        out_[m] = dict(threshold=thr, max_abs_diff=float(dd.max()), median_abs_diff=float(dd.median()),
                       n_structures_above_threshold=int((dd > thr).sum()), n_compared=int(dd.notna().sum()),
                       action="hard stop" if m == "supplied_centroid_mm" else "flag only")
    return out_


# ------------------------------ pre-flight: report implementation checks and stop ------------------------------
if PREFLIGHT:
    SLp = pd.DataFrame(struct_rows); Vp = pd.DataFrame(verify_rows); Cp = pd.DataFrame(cover_rows)
    keep = ["case_id", "label_id", "structure", "set", "supplied_centroid_mm", "A1_supplied_centroid_mm",
            "supplied_dice", "A1_supplied_dice", "supplied_MSD_mm", "A1_supplied_MSD_mm",
            "supplied_HD95_mm", "A1_supplied_HD95_mm"]
    SLp = SLp[keep] if len(SLp) else SLp
    say("\n===== PRE-FLIGHT SUMMARY (implementation checks only) =====")
    say(f"Patients processed: {len(patient_rows)}; statuses: {pd.Series([r['status'] for r in patient_rows]).value_counts().to_dict()}")
    if len(Vp):
        say("Verification (a)-(c):\n" + Vp.to_string(index=False))
    if len(SLp):
        for m, r in repro_flags(SLp).items():
            say(f"Reproduction of Amendment 1 {m} ({r['action']}, threshold {r['threshold']}): max |diff| "
                f"{r['max_abs_diff']:.4g}, median {r['median_abs_diff']:.4g}, above threshold {r['n_structures_above_threshold']}/{r['n_compared']}")
        badp, _ = completeness(evaluated_counts(SLp, ["supplied_centroid_mm"]), [r["case_id"] for r in patient_rows])
        say(f"Patients lacking a required held-out measurement (supplied): {badp if badp else 'none'}")
    if len(Cp):
        say("Supplied coverage (compare with Audit 1C-b FDG_TTB_in_PSMA_PET_FOV):\n"
            + Cp[["case_id", "supplied_FDG_TTB_in_PSMA_PET_FOV", "supplied_FDG_TTB_in_PSMA_PET_FOV_10mm"]].to_string(index=False))
    SLp.to_csv(out("supplied_reproduction", "csv"), index=False)
    Vp.to_csv(out("verification_log", "csv"), index=False)
    Cp[[c for c in Cp.columns if not c.startswith("corrected")]].to_csv(out("supplied_coverage", "csv"), index=False)
    out("report", "txt").write_text("\n".join(LOG), encoding="utf-8")
    say("PRE-FLIGHT finished. No Amendment 2 result was computed.")
    sys.exit(0)

# ------------------------------ patient-level summaries ------------------------------
SL = pd.DataFrame(struct_rows)
PT = pd.DataFrame(patient_rows)
for s in SETS:
    g = SL[SL["set"] == s].groupby("case_id")
    PT = PT.merge(g["supplied_centroid_mm"].median().rename(f"{s}_supplied_centroid_mm"), on="case_id", how="left")
    PT = PT.merge(g["corrected_centroid_mm"].median().rename(f"{s}_corrected_centroid_mm"), on="case_id", how="left")
    PT[f"{s}_diff_centroid_mm"] = PT[f"{s}_corrected_centroid_mm"] - PT[f"{s}_supplied_centroid_mm"]
    for m in ["dice", "MSD_mm", "HD95_mm"]:
        PT = PT.merge(g[f"supplied_{m}"].median().rename(f"{s}_supplied_{m}"), on="case_id", how="left")
        PT = PT.merge(g[f"corrected_{m}"].median().rename(f"{s}_corrected_{m}"), on="case_id", how="left")

cnt_all = evaluated_counts(SL, ["supplied_centroid_mm", "corrected_centroid_mm"])
PT = PT.merge(cnt_all, on="case_id", how="left")
for c in [c for c in PT.columns if c.startswith("n_") and c.endswith("_evaluated")]:
    PT[c] = PT[c].fillna(0).astype(int)
missing_heldout, _ = completeness(cnt_all, PT["case_id"].tolist())
n_failed = int((PT["status"] == "FAILED").sum())
evaluation_complete = (not missing_heldout) and n_failed == 0
if not evaluation_complete:
    say(f"*** EVALUATION INCOMPLETE: patients lacking a required held-out measurement: {missing_heldout}; "
        f"processing failures: {n_failed}. C1-C3 are INVALID; the run cannot be COMPLETE. ***")

# ------------------------------ criteria C1-C3 (before any fallback) ------------------------------
def criteria(df, label):
    hb_s, hb_c, hb_d = df["heldout_bone_supplied_centroid_mm"], df["heldout_bone_corrected_centroid_mm"], df["heldout_bone_diff_centroid_mm"]
    hs_d = df["heldout_soft_diff_centroid_mm"]
    n_worse = int((hb_d > C3_WORSEN_MM).sum())
    c1 = bool(hb_c.median() < hb_s.median() and hb_d.median() < C1_DIFF_MAX)
    c2 = bool(hs_d.median() <= C2_DIFF_MAX)
    c3 = bool(n_worse <= C3_MAX_PATIENTS)
    return dict(population=label, n_patients=len(df),
                heldout_bone_supplied_median=hb_s.median(), heldout_bone_corrected_median=hb_c.median(),
                heldout_bone_median_paired_diff=hb_d.median(),
                heldout_soft_median_paired_diff=hs_d.median(),
                n_patients_heldout_bone_worse_gt_2mm=n_worse,
                n_improved_gt_1mm=int((hb_d < -1).sum()), n_worsened_gt_1mm=int((hb_d > 1).sum()),
                C1_met=c1, C2_met=c2, C3_met=c3, C4="PENDING_AUTHOR_VISUAL_REVIEW")

ok_pts = PT[PT["status"] != "FAILED"]
CR = pd.DataFrame([criteria(ok_pts[ok_pts["status"] == "CORRECTION_COMPUTED"], "(i) patients with computed correction [decides]"),
                   criteria(ok_pts, "(ii) all patients (insufficient-landmarks enter with difference 0)")])
CR["n_insufficient_landmarks"] = int((PT["status"] == "INSUFFICIENT_LANDMARKS_FALLBACK").sum())
CR["criteria_valid"] = evaluation_complete
CR["patients_missing_required_heldout"] = ";".join(missing_heldout)

# per-patient fallback (applies ONLY if the author adopts the corrected transformation)
PT["would_fallback_if_adopted"] = (PT["status"] != "CORRECTION_COMPUTED") | (PT["heldout_bone_diff_centroid_mm"] > C3_WORSEN_MM)
PT.loc[(PT["status"] == "CORRECTION_COMPUTED") & (PT["heldout_bone_diff_centroid_mm"] > C3_WORSEN_MM),
       "fallback_reason"] = f"held-out bone worsened > {C3_WORSEN_MM} mm"

# ------------------------------ reproduction of Amendment 1 supplied values (centroids enforced above) ------------------------------
rep = repro_flags(SL) if len(SL) else {}

# ------------------------------ C4 overlays: Amendment 1 overlay() VERBATIM (labels "refined" -> "corrected") ------------------------------
ALL_LABELS = sorted(int(x) for x in pd.read_csv(A1_STRUCT)["label_id"].unique())


def overlay(cid, ct, tsP, res_s, res_r, path_png, txt):
    ctO = sitk.GetArrayFromImage(sitk.DICOMOrient(ct, "LPS"))
    pO = sitk.GetArrayFromImage(sitk.DICOMOrient(tsP, "LPS"))
    sp = sitk.DICOMOrient(ct, "LPS").GetSpacing()
    zz, yy, xx = np.nonzero(np.isin(pO, VERT))
    y0, x0 = int(np.median(yy)), int(np.median(xx))
    fig, axs = plt.subplots(2, 2, figsize=(12, 15))
    for j, (name, res) in enumerate([("supplied", res_s), ("corrected", res_r)]):
        rO = sitk.GetArrayFromImage(sitk.DICOMOrient(res, "LPS"))
        for i, (sl, w, ttl) in enumerate([(np.s_[:, y0, :], sp[0], "coronal (pt right = image left)"),
                                           (np.s_[:, :, x0], sp[1], "sagittal (anterior = left)")]):
            ax = axs[i, j]
            img = np.flipud(ctO[sl])
            ext = (0, img.shape[1] * w, 0, img.shape[0] * sp[2])
            ax.imshow(img, cmap="gray", vmin=WINDOW[0], vmax=WINDOW[1], extent=ext, origin="upper")
            a1 = np.flipud(np.isin(pO[sl], ALL_LABELS)).astype(float)
            a2 = np.flipud(np.isin(rO[sl], ALL_LABELS)).astype(float)
            if a1.any():
                ax.contour(a1, levels=[0.5], colors="lime", linewidths=0.7, extent=ext, origin="upper")
            if a2.any():
                ax.contour(a2, levels=[0.5], colors="magenta", linewidths=0.7, extent=ext, origin="upper")
            ax.set_title("%s | %s transform | %s" % (cid, name, ttl), fontsize=9)
            ax.set_xlabel("mm")
            ax.set_ylabel("mm")
    fig.suptitle("%s: PSMA CT; green = PSMA labels; magenta = FDG labels resampled\n%s" % (cid, txt), fontsize=10)
    fig.savefig(path_png, dpi=100)
    plt.close(fig)


# ------------------------------ C4 selection (PENDING author decision; see C4_PROCEDURE_VERIFIED) ------------------------------
primary = ok_pts.sort_values("A1_change_median_mm", ascending=False, kind="mergesort").head(C4_N)["case_id"].tolist()
cand = ok_pts[ok_pts["status"] == "CORRECTION_COMPUTED"].dropna(subset=["heldout_bone_diff_centroid_mm"])
worst = cand.nlargest(1, "heldout_bone_diff_centroid_mm")["case_id"].tolist()
best = cand.nsmallest(1, "heldout_bone_diff_centroid_mm")["case_id"].tolist()
safety = [c for c in dict.fromkeys(worst + best) if c not in primary]
c4_sel = primary + safety
c4_rows, overlay_files = [], []
for cid in c4_sel:
    reason = ("primary (Amendment 1 rule: largest change)" if cid in primary else
              "safety extension: largest held-out-bone deterioration" if cid in worst else
              "safety extension: largest held-out-bone improvement")
    c4_rows.append(dict(case_id=cid, selection_reason=reason,
                        heldout_bone_diff_centroid_mm=float(PT.loc[PT.case_id == cid, "heldout_bone_diff_centroid_mm"].iloc[0]),
                        t_magnitude_mm=float(PT.loc[PT.case_id == cid, "t_magnitude_mm"].iloc[0]),
                        A1_change_median_mm=float(PT.loc[PT.case_id == cid, "A1_change_median_mm"].iloc[0]),
                        visual_verdict="TO_BE_RECORDED_BY_AUTHOR"))
    try:
        dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
        ct = sitk.ReadImage(str(dP / FILES["ct"]))
        tsP = sitk.Cast(sitk.ReadImage(str(dP / FILES["totseg"])), sitk.sitkUInt8)
        tsF = sitk.Cast(sitk.ReadImage(str(dF / FILES["totseg"])), sitk.sitkUInt8)
        R = sitk.ReadTransform(str(dP / FILES["tfm"])); T = sitk.ReadTransform(str(dF / FILES["tfm"]))
        prow = PT[PT.case_id == cid].iloc[0]
        t = np.array([prow.t_x_LPS, prow.t_y_LPS, prow.t_z_LPS])
        _, Rc = corrected_transforms(T, R, t)
        rs = sitk.Resample(tsF, tsP, R, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
        rr = sitk.Resample(tsF, tsP, Rc, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
        f = out(f"C4_overlay_{cid}", "png")
        overlay(cid, ct, tsP, rs, rr, str(f), "held-out bone %.1f -> %.1f mm; translation |t| %.1f mm (criterion C4; %s)" % (
            prow.heldout_bone_supplied_centroid_mm, prow.heldout_bone_corrected_centroid_mm, prow.t_magnitude_mm, reason))
        overlay_files.append(str(f))
    except Exception as e:
        say(f"C4 overlay error for {cid}: {e}")

# ------------------------------ write outputs ------------------------------
errors = []
outputs = {"per_patient": PT, "structure_level": SL, "criteria_C1_C3": CR, "verification_log": pd.DataFrame(verify_rows),
           "coverage": pd.DataFrame(cover_rows), "C4_selection": pd.DataFrame(c4_rows),
           "input_inventory": pd.DataFrame(inventory)}
written = []
for name, df in outputs.items():
    try:
        f = out(name, "csv"); df.to_csv(f, index=False); written.append(str(f))
    except Exception as e:
        errors.append(f"{name}: {e}")
written += overlay_files

say("\n===== SUMMARY =====")
say(f"Patients: {len(PT)} | computed {int((PT.status == 'CORRECTION_COMPUTED').sum())} | "
    f"insufficient landmarks {int((PT.status == 'INSUFFICIENT_LANDMARKS_FALLBACK').sum())} | failed {int((PT.status == 'FAILED').sum())}")
say("Criteria (before any fallback):\n" + CR.to_string(index=False))
say(f"Reproduction of Amendment 1 supplied values: {rep}")
say(f"C4 selection: primary {primary}; safety extension {safety} (overlays written: {len(overlay_files)})")
say(f"Structures whose corrected-transform volume ratio is outside the Amendment 1 range {VOL_RATIO_OK} "
    f"(reported; structure sets are fixed): {int(SL['corrected_volume_ratio_outside_A1_range'].sum()) if len(SL) else 0}")

n_ok = int((PT["status"] != "FAILED").sum())
complete = (n_found == len(cases) == EXPECTED_PATIENTS and n_ok == EXPECTED_PATIENTS and not errors
            and evaluation_complete
            and len(primary) == C4_N and len(overlay_files) == len(c4_sel))
completion = "COMPLETE" if complete else "INCOMPLETE"

manifest = dict(step="Audit 1C-e - Amendment 2 translation-only refinement (protocol v2.8 Section 15.2)",
                run=RUN, completion_status=completion, script=str(script_path), script_MD5=SCRIPT_MD5,
                amendment1_structure_table=dict(path=str(A1_STRUCT), md5=md5(A1_STRUCT)),
                parameters=dict(MIN_FITTING=MIN_FITTING, VERIFY_TOL_MM=VERIFY_TOL_MM, N_TEST_GRID=N_TEST_GRID,
                                C1_DIFF_MAX=C1_DIFF_MAX, C2_DIFF_MAX=C2_DIFF_MAX, C3_WORSEN_MM=C3_WORSEN_MM,
                                C3_MAX_PATIENTS=C3_MAX_PATIENTS, C4_N=C4_N, C4_rule="primary: top C4_N by Amendment 1 change; safety: largest held-out deterioration and improvement",
                                COVERAGE_MARGIN_MM=COVERAGE_MARGIN_MM, CENTROID_REPRO_TOL_MM=CENTROID_REPRO_TOL_MM,
                                REPRO_FLAG_DICE=REPRO_FLAG_DICE, REPRO_FLAG_MSD_MM=REPRO_FLAG_MSD_MM,
                                REPRO_FLAG_HD95_MM=REPRO_FLAG_HD95_MM, REQUIRED_HELDOUT_SETS=REQUIRED_HELDOUT_SETS,
                                estimator="component-wise median of d_s over Amendment 1 fitting structures",
                                transforms="T'(p)=FDG_rigid(p)-t ; R'(q)=PSMA_rigid(q+t)",
                                surface_metrics="pooled symmetric contour-to-boundary distances; MSD mean; HD95 95th pct"),
                criteria=CR.to_dict("records"), reproduction_check_vs_amendment1=rep,
                evaluation_complete=evaluation_complete, patients_missing_required_heldout=missing_heldout,
                C4_procedure_verified=C4_PROCEDURE_VERIFIED,
                C4_selection=c4_sel, output_errors=errors,
                software=dict(python=sys.version, SimpleITK=sitk.Version_VersionString(),
                              numpy=np.__version__, pandas=pd.__version__, matplotlib=matplotlib.__version__,
                              platform=platform.platform()),
                n_input_files=len(inventory), outputs=written,
                note="Adoption decision is recorded by the author after C4 review. No correspondence or discordance.")
try:
    fM = out("manifest", "json"); tmp = Path(str(fM) + ".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    chk = json.loads(tmp.read_text(encoding="utf-8"))
    if chk.get("completion_status") != completion:
        raise RuntimeError("manifest read-back mismatch")
    tmp.replace(fM)
except Exception as e:
    completion = "INCOMPLETE"; say(f"MANIFEST ERROR: {e}")
say(f"\nAUDIT 1C-e RUN STATUS: {completion}  (C4 pending author visual review; adoption decision by author)")
out("report", "txt").write_text("\n".join(LOG), encoding="utf-8")
out(f"RUN_{completion}", "txt").write_text(f"Audit 1C-e run {RUN}: {completion}\n", encoding="utf-8")
