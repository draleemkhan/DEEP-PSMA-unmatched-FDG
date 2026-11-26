# -*- coding: utf-8 -*-
"""
Audit1C_c_rigid_refinement.py
===========================================================================
PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - AUDIT 1C-c (Protocol Amendment 1)
Bone-based rigid refinement of the supplied PSMA-FDG registration, with
held-out anatomical evaluation in all 100 patients.
Protocol: PAPER3_PSMA_FDG_discordance_protocol_v2.5, Section 15 (Amendment 1)

READ-ONLY on:   D:\\DEEP PSMA\\DEEP PSMA EXTRACTED
Reads (this study's output): Audit 1B label-level table
Writes ONLY to: D:\\DEEP PSMA\\PAPER3_Discordance\\06_output
  (refined transforms saved there as Audit1C_c_tfm_<case>_*.tfm; no subfolders)

Anatomy only. Tumour masks used ONLY for coverage fractions. NO tumour
correspondence, NO DMTVF, NO discordance results.

Pre-specified method (Amendment 1)
  fixed = PSMA CT, moving = FDG CT, Euler3D initialised with PSMA/rigid.tfm
  (SimpleITK convention: maps PSMA points -> FDG points)
  Mattes MI (50 bins), fixed mask = fitting structures dilated 5 mm,
  random sampling 10% (fixed seed), linear interpolation,
  regular-step gradient descent, physical-shift scales,
  shrink 8/4/2, smoothing 4/2/1 mm.
  Fitting structures: vertebrae T1-L5, S1, sacrum, hips (eligible only).
Evaluation (supplied vs refined; structures eligible under BOTH)
  held-out bone (primary): ribs, sternum, clavicles, scapulae
  held-out soft: liver, spleen, kidneys, aorta, heart, lung lobes, bladder
  fitting (in-sample) and head/neck (C1-C7, brain) reported separately
  centroid distance (point mapping), Dice, MSD, HD95
Pre-specified adoption criteria C1-C3 are evaluated and printed; C4
(visual) uses the overlays; the final decision is recorded by the author.

Usage:
  python Audit1C_c_rigid_refinement.py                 (all patients)
  python Audit1C_c_rigid_refinement.py train_0001 ...  (selected patients,
                                                        e.g. smoke test)
Requires: Python 3, SimpleITK >= 2.0, numpy, matplotlib.
===========================================================================
"""
import os
import sys
import csv
import json
import math
import time
import platform
import datetime

import numpy as np

try:
    import SimpleITK as sitk
except ImportError:
    sys.exit("SimpleITK not found in this Python environment.")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ config
DATA_ROOT = r"D:\DEEP PSMA\DEEP PSMA EXTRACTED"
OUT_DIR = r"D:\DEEP PSMA\PAPER3_Discordance\06_output"
AUDIT1B_LABEL_CSV = os.path.join(OUT_DIR, "Audit1B_label_level_20261029_134401.csv")
TS = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
TFM_DIR = OUT_DIR   # single locked output folder (no subfolders)

SEED = 20261029
MI_BINS = 50
SAMPLING = 0.10
SHRINK = [8, 4, 2]
SMOOTH = [4.0, 2.0, 1.0]
LR, MIN_STEP, N_ITER, RELAX = 2.0, 1e-3, 200, 0.5
MASK_DILATE_MM = 5.0

CAND_LIMITS = [5, 10, 15, 20]
FRAG_MIN = 0.95
VOL_RATIO_OK = (0.95, 1.05)
PET_BODY_SUV = 0.15
FOV_MARGIN_MM = 10.0
WORSEN_MM = 2.0
N_OVERLAYS = 5
WINDOW = (-200, 1000)

FIT = list(range(26, 44)) + [25, 77, 78]                      # S1..T1, sacrum, hips
SETS = {
    "heldout_bone": list(range(92, 116)) + [116, 73, 74, 71, 72],
    "heldout_soft": [5, 1, 2, 3, 52, 51, 10, 11, 12, 13, 14, 21],
    "fitting": FIT,
    "head_neck": list(range(44, 51)) + [90],
}
VERT = list(range(25, 51))
FLAGGED = {"train_0088": [39, 40, 41, 42, 43],
           "train_0059": [26, 27, 28, 29, 30],
           "train_0029": [29, 30, 31, 32, 33]}

NAMES = {1: "spleen", 2: "kidney_right", 3: "kidney_left", 5: "liver", 10: "lung_upper_lobe_left",
         11: "lung_lower_lobe_left", 12: "lung_upper_lobe_right", 13: "lung_middle_lobe_right",
         14: "lung_lower_lobe_right", 21: "urinary_bladder", 25: "sacrum", 51: "heart", 52: "aorta",
         71: "scapula_left", 72: "scapula_right", 73: "clavicula_left", 74: "clavicula_right",
         77: "hip_left", 78: "hip_right", 90: "brain", 116: "sternum"}
for i, v in enumerate(["S1", "L5", "L4", "L3", "L2", "L1", "T12", "T11", "T10", "T9", "T8", "T7", "T6",
                       "T5", "T4", "T3", "T2", "T1", "C7", "C6", "C5", "C4", "C3", "C2", "C1"]):
    NAMES[26 + i] = "vertebrae_" + v
for i in range(12):
    NAMES[92 + i] = "rib_left_%d" % (i + 1)
    NAMES[104 + i] = "rib_right_%d" % (i + 1)
SET_OF = {}
for sname, labs in SETS.items():
    for l in labs:
        SET_OF[l] = sname
ALL = sorted(SET_OF)


def out(name, ext):
    return os.path.join(OUT_DIR, "Audit1C_c_%s_%s.%s" % (name, TS, ext))


class Tee(object):
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.o = sys.__stdout__

    def write(self, s):
        self.o.write(s)
        self.f.write(s)

    def flush(self):
        self.o.flush()
        self.f.flush()


# ------------------------------------------------------------------ helpers
def dc(t):
    try:
        return t.Downcast()
    except Exception:
        return t


def load_tfm(path):
    t = sitk.ReadTransform(path)
    try:
        t = t.Downcast()
    except Exception:
        pass
    return t


class Grid(object):
    def __init__(self, src):
        self.o = np.array(src.GetOrigin(), float)
        self.s = np.array(src.GetSpacing(), float)
        self.D = np.array(src.GetDirection(), float).reshape(3, 3)
        self.size = np.array(src.GetSize(), float)
        self.Dinv = np.linalg.inv(self.D)

    def idx2phys(self, idx_xyz):
        return self.o + (idx_xyz * self.s).dot(self.D.T)

    def inside(self, p, margin_mm=0.0):
        ci = (p - self.o).dot(self.Dinv.T) / self.s
        m = margin_mm / self.s
        return np.all((ci >= -0.5 + m) & (ci <= self.size - 0.5 - m), axis=1)


def rigid_apply(T, p):
    R = np.array(T.GetMatrix()).reshape(3, 3)
    c = np.array(T.GetCenter())
    t = np.array(T.GetTranslation())
    return (p - c).dot(R.T) + c + t


def mask_points(img, grid):
    arr = sitk.GetArrayViewFromImage(img)
    z, y, x = np.nonzero(arr)
    return grid.idx2phys(np.stack([x, y, z], axis=1).astype(float))


def pet_body(pet):
    m = sitk.BinaryThreshold(pet, PET_BODY_SUV, 1e12, 1, 0)
    rel = sitk.RelabelComponent(sitk.ConnectedComponent(m, True), sortByObjectSize=True)
    return sitk.BinaryThreshold(rel, 1, 1, 1, 0)


def label_stats(img):
    f = sitk.LabelShapeStatisticsImageFilter()
    f.Execute(img)
    return {int(l): dict(vol=f.GetPhysicalSize(l), cen=np.array(f.GetCentroid(l)), bb=f.GetBoundingBox(l))
            for l in f.GetLabels()}


def surface_metrics(imgA, imgB, l, bbA, bbB, margin=3):
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


def med(v):
    v = [x for x in v if x == x]
    return float(np.median(v)) if v else float("nan")


def stats(v):
    v = np.array([x for x in v if x == x], float)
    if v.size == 0:
        return dict(n=0, median=float("nan"), Q1=float("nan"), Q3=float("nan"), P90=float("nan"),
                    min=float("nan"), max=float("nan"))
    return dict(n=int(v.size), median=float(np.median(v)), Q1=float(np.percentile(v, 25)),
                Q3=float(np.percentile(v, 75)), P90=float(np.percentile(v, 90)),
                min=float(v.min()), max=float(v.max()))


def write_csv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def read_audit1b(path):
    d = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                fr = float(r["frac_largest_component"])
            except (ValueError, KeyError):
                fr = float("nan")
            d[(r["case_id"], r["tracer"], int(r["label_id"]))] = (r["category"], fr)
    return d


def transform_change(T0, T1, pts):
    """Displacement (mm) between two resampling transforms at given PSMA points."""
    a = rigid_apply(T0, pts)
    b = rigid_apply(T1, pts)
    return np.linalg.norm(a - b, axis=1)


def refine(fixed, moving, mask, Tinit):
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMattesMutualInformation(MI_BINS)
    R.SetMetricSamplingStrategy(R.RANDOM)
    R.SetMetricSamplingPercentage(SAMPLING, SEED)
    R.SetMetricFixedMask(mask)
    R.SetInterpolator(sitk.sitkLinear)
    R.SetOptimizerAsRegularStepGradientDescent(learningRate=LR, minStep=MIN_STEP,
                                               numberOfIterations=N_ITER, relaxationFactor=RELAX)
    R.SetOptimizerScalesFromPhysicalShift()
    R.SetShrinkFactorsPerLevel(SHRINK)
    R.SetSmoothingSigmasPerLevel(SMOOTH)
    R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    init = sitk.Euler3DTransform()
    init.SetCenter(Tinit.GetCenter())
    init.SetMatrix(Tinit.GetMatrix())
    init.SetTranslation(Tinit.GetTranslation())
    R.SetInitialTransform(init, inPlace=False)
    T = R.Execute(fixed, moving)
    try:
        T = T.Downcast()
    except Exception:
        pass
    if not isinstance(T, sitk.Euler3DTransform):
        # composite wrapper -> extract the Euler transform
        try:
            T = sitk.CompositeTransform(T).GetNthTransform(0).Downcast()
        except Exception:
            pass
    if not hasattr(T, "GetMatrix"):
        raise RuntimeError("refined transform is not a rigid Euler3D transform")
    return T, R.GetMetricValue(), R.GetOptimizerIteration(), R.GetOptimizerStopConditionDescription()


def evaluate(tsP, tsF, stP, stF, Tres, Tpm, labels):
    """Metrics for one transform pair (resampling Tres, point mapping Tpm)."""
    res = sitk.Resample(tsF, tsP, Tres, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    stR = label_stats(res)
    ov = sitk.LabelOverlapMeasuresImageFilter()
    ov.Execute(tsP, res)
    outd = {}
    for l in labels:
        if l not in stR:
            outd[l] = None
            continue
        vr = stR[l]["vol"] / stF[l]["vol"]
        if not (VOL_RATIO_OK[0] <= vr <= VOL_RATIO_OK[1]):
            outd[l] = None
            continue
        m = np.array(Tpm.TransformPoint(tuple(float(x) for x in stF[l]["cen"])))
        try:
            dc = ov.GetDiceCoefficient(int(l))
        except Exception:
            dc = float("nan")
        msd, hd = surface_metrics(tsP, res, l, stP[l]["bb"], stR[l]["bb"])
        outd[l] = dict(cd=float(np.linalg.norm(stP[l]["cen"] - m)), dice=dc, msd=msd, hd95=hd, vr=vr)
    return outd, res


def overlay(cid, ct, tsP, res_s, res_r, path_png, txt):
    ctO = sitk.GetArrayFromImage(sitk.DICOMOrient(ct, "LPS"))
    pO = sitk.GetArrayFromImage(sitk.DICOMOrient(tsP, "LPS"))
    sp = sitk.DICOMOrient(ct, "LPS").GetSpacing()
    zz, yy, xx = np.nonzero(np.isin(pO, VERT))
    y0, x0 = int(np.median(yy)), int(np.median(xx))
    fig, axs = plt.subplots(2, 2, figsize=(12, 15))
    for j, (name, res) in enumerate([("supplied", res_s), ("refined", res_r)]):
        rO = sitk.GetArrayFromImage(sitk.DICOMOrient(res, "LPS"))
        for i, (sl, w, ttl) in enumerate([(np.s_[:, y0, :], sp[0], "coronal (pt right = image left)"),
                                           (np.s_[:, :, x0], sp[1], "sagittal (anterior = left)")]):
            ax = axs[i, j]
            img = np.flipud(ctO[sl])
            ext = (0, img.shape[1] * w, 0, img.shape[0] * sp[2])
            ax.imshow(img, cmap="gray", vmin=WINDOW[0], vmax=WINDOW[1], extent=ext, origin="upper")
            a1 = np.flipud(np.isin(pO[sl], ALL)).astype(float)
            a2 = np.flipud(np.isin(rO[sl], ALL)).astype(float)
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


# ------------------------------------------------------------------ main
def main():
    assert os.path.isdir(DATA_ROOT), "Data root not found"
    assert os.path.isfile(AUDIT1B_LABEL_CSV), "Audit 1B label-level table not found"
    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = out("run_log", "txt")
    sys.stdout = Tee(log_path)
    t0 = time.time()
    print("Audit 1C-c (Amendment 1) started %s" % datetime.datetime.now())
    print("Python %s | SimpleITK %s | numpy %s" % (platform.python_version(), sitk.Version.VersionString(), np.__version__))
    a1b = read_audit1b(AUDIT1B_LABEL_CSV)
    cases = sorted(d for d in os.listdir(DATA_ROOT)
                   if d.startswith("train_") and os.path.isdir(os.path.join(DATA_ROOT, d)))
    if len(sys.argv) > 1:
        cases = [c for c in cases if c in sys.argv[1:]]
    print("Patients: %d\n" % len(cases))

    s_rows, p_rows, change_rank = [], [], []
    for ci, cid in enumerate(cases):
        tc = time.time()
        dP = os.path.join(DATA_ROOT, cid, "PSMA")
        dF = os.path.join(DATA_ROOT, cid, "FDG")
        P = dict(case_id=cid, error="")
        try:
            Tp = load_tfm(os.path.join(dP, "rigid.tfm"))
            Tf = load_tfm(os.path.join(dF, "rigid.tfm"))
            ctP = sitk.Cast(sitk.ReadImage(os.path.join(dP, "CT.nii.gz")), sitk.sitkFloat32)
            ctF = sitk.Cast(sitk.ReadImage(os.path.join(dF, "CT.nii.gz")), sitk.sitkFloat32)
            tsP = sitk.Cast(sitk.ReadImage(os.path.join(dP, "totseg_24.nii.gz")), sitk.sitkUInt8)
            tsF = sitk.Cast(sitk.ReadImage(os.path.join(dF, "totseg_24.nii.gz")), sitk.sitkUInt8)
            stP, stF = label_stats(tsP), label_stats(tsF)

            # eligibility (Audit 1B completeness + vertebral rules)
            excl = set(FLAGGED.get(cid, []))
            for l in VERT:
                fp = a1b.get((cid, "PSMA", l), ("", float("nan")))[1]
                ff = a1b.get((cid, "FDG", l), ("", float("nan")))[1]
                if (fp == fp and fp < FRAG_MIN) or (ff == ff and ff < FRAG_MIN):
                    excl.update([l - 1, l, l + 1])
            elig = [l for l in ALL
                    if a1b.get((cid, "PSMA", l), ("NA", 0))[0] == "PRESENT_COMPLETE"
                    and a1b.get((cid, "FDG", l), ("NA", 0))[0] == "PRESENT_COMPLETE"
                    and l not in excl and l in stP and l in stF]
            fit = [l for l in elig if l in FIT]
            P["n_fitting_structures"] = len(fit)
            if len(fit) < 3:
                raise RuntimeError("fewer than 3 fitting structures")

            # fixed mask = fitting structures dilated
            arr = np.isin(sitk.GetArrayViewFromImage(tsP), fit).astype(np.uint8)
            mask = sitk.GetImageFromArray(arr)
            mask.CopyInformation(tsP)
            rad = [int(math.ceil(MASK_DILATE_MM / s)) for s in tsP.GetSpacing()]
            mask = sitk.BinaryDilate(mask, rad)
            del arr

            # refinement (fallback to supplied if it fails)
            fallback_reason = ""
            try:
                Tr, mval, nit, stop = refine(ctP, ctF, mask, Tp)
                Tr_inv = dc(Tr.GetInverse())
                P.update(final_metric=mval, iterations_last_level=nit, stop_condition=stop)
            except Exception as ex:
                Tr, Tr_inv = Tp, Tf
                fallback_reason = "registration_failed: %s" % ex
            # change produced by refinement over the fitting-structure centroids (PSMA space)
            pts = np.array([stP[l]["cen"] for l in fit])
            chg = transform_change(Tp, Tr, pts)
            P.update(refinement_change_median_mm=float(np.median(chg)), refinement_change_max_mm=float(np.max(chg)))
            if Tr is not Tp:
                R0 = np.array(Tp.GetMatrix()).reshape(3, 3)
                R1 = np.array(Tr.GetMatrix()).reshape(3, 3)
                Rd = R1.dot(R0.T)
                P["refinement_rotation_change_deg"] = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(Rd) - 1) / 2))))
            sitk.WriteTransform(Tr, os.path.join(TFM_DIR, "Audit1C_c_tfm_%s_refined_resampling_%s.tfm" % (cid, TS)))
            sitk.WriteTransform(Tr_inv, os.path.join(TFM_DIR, "Audit1C_c_tfm_%s_refined_pointmap_FDG_to_PSMA_%s.tfm" % (cid, TS)))

            # evaluation: supplied vs refined
            ev_s, res_s = evaluate(tsP, tsF, stP, stF, Tp, Tf, elig)
            ev_r, res_r = evaluate(tsP, tsF, stP, stF, Tr, Tr_inv, elig)
            for l in elig:
                a, b = ev_s[l], ev_r[l]
                row = dict(case_id=cid, label_id=l, structure=NAMES[l], set=SET_OF[l],
                           evaluable_both=(a is not None and b is not None))
                for tag, e in (("supplied", a), ("refined", b)):
                    if e is not None:
                        row.update({"%s_centroid_mm" % tag: e["cd"], "%s_dice" % tag: e["dice"],
                                    "%s_MSD_mm" % tag: e["msd"], "%s_HD95_mm" % tag: e["hd95"]})
                s_rows.append(row)
            both = [r for r in s_rows if r["case_id"] == cid and r["evaluable_both"]]
            for sname in SETS:
                g = [r for r in both if r["set"] == sname]
                P["n_%s" % sname] = len(g)
                for tag in ("supplied", "refined"):
                    P["%s_%s_centroid_mm" % (sname, tag)] = med([r["%s_centroid_mm" % tag] for r in g])
                    P["%s_%s_MSD_mm" % (sname, tag)] = med([r["%s_MSD_mm" % tag] for r in g])
                    P["%s_%s_dice" % (sname, tag)] = med([r["%s_dice" % tag] for r in g])
                P["%s_diff_centroid_mm" % sname] = P["%s_refined_centroid_mm" % sname] - P["%s_supplied_centroid_mm" % sname]
                P["%s_diff_MSD_mm" % sname] = P["%s_refined_MSD_mm" % sname] - P["%s_supplied_MSD_mm" % sname]

            # per-patient safeguard (C3): fallback if held-out bone worsens > 2 mm or registration failed
            d = P["heldout_bone_diff_centroid_mm"]
            if not fallback_reason and d == d and d > WORSEN_MM:
                fallback_reason = "heldout_bone_worse_than_%.0fmm" % WORSEN_MM
            P["fallback_to_supplied"] = bool(fallback_reason)
            P["fallback_reason"] = fallback_reason

            # coverage with supplied and refined point mapping (Section 8.1)
            petP = sitk.ReadImage(os.path.join(dP, "PET.nii.gz"))
            petF = sitk.ReadImage(os.path.join(dF, "PET.nii.gz"))
            gPP, gPF = Grid(petP), Grid(petF)
            ptF = mask_points(sitk.ReadImage(os.path.join(dF, "TTB.nii.gz")) > 0, gPF)
            bF = mask_points(pet_body(petF), gPF)
            frac = lambda x: float(np.mean(x)) if x.size else float("nan")
            for tag, M in (("supplied", Tf), ("refined", Tr_inv)):
                P["FDG_TTB_in_PSMA_PET_FOV_%s" % tag] = frac(gPP.inside(rigid_apply(M, ptF)))
                P["FDG_TTB_in_PSMA_PET_FOV_10mm_%s" % tag] = frac(gPP.inside(rigid_apply(M, ptF), FOV_MARGIN_MM))
                P["FDG_body_in_PSMA_PET_FOV_%s" % tag] = frac(gPP.inside(rigid_apply(M, bF)))

            change_rank.append((P["refinement_change_median_mm"], cid))
            print("[%3d/%3d] %s fit %d | held-out bone %.1f -> %.1f mm | held-out soft %.1f -> %.1f mm | "
                  "fitting %.1f -> %.1f | change %.1f mm | %s| %.0f s" % (
                      ci + 1, len(cases), cid, len(fit),
                      P["heldout_bone_supplied_centroid_mm"], P["heldout_bone_refined_centroid_mm"],
                      P["heldout_soft_supplied_centroid_mm"], P["heldout_soft_refined_centroid_mm"],
                      P["fitting_supplied_centroid_mm"], P["fitting_refined_centroid_mm"],
                      P["refinement_change_median_mm"],
                      ("FALLBACK(%s) " % fallback_reason) if fallback_reason else "", time.time() - tc))
            del ctP, ctF, tsP, tsF, res_s, res_r, mask, petP, petF
        except Exception as ex:
            P["error"] = str(ex)
            print("[%3d/%3d] %s ERROR: %s" % (ci + 1, len(cases), cid, ex))
        p_rows.append(P)

    # ---------------- cohort summary and pre-specified criteria ----------------
    ok = [p for p in p_rows if not p.get("error")]
    c_rows = []
    for sname in SETS:
        for metric in ("centroid_mm", "MSD_mm"):
            vs = [p.get("%s_supplied_%s" % (sname, metric), float("nan")) for p in ok]
            vr = [p.get("%s_refined_%s" % (sname, metric), float("nan")) for p in ok]
            dd = [p.get("%s_diff_%s" % (sname, metric), float("nan")) for p in ok]
            for tag, v in (("supplied", vs), ("refined", vr), ("diff_refined_minus_supplied", dd)):
                row = dict(set=sname, metric=metric, transform=tag, **stats(v))
                vv = np.array([x for x in v if x == x])
                if tag != "diff_refined_minus_supplied":
                    for L in CAND_LIMITS:
                        row["pct_patients_gt_%dmm" % L] = 100.0 * np.mean(vv > L) if vv.size else float("nan")
                else:
                    row["n_improved_gt_1mm"] = int(np.sum(vv < -1))
                    row["n_worsened_gt_1mm"] = int(np.sum(vv > 1))
                    row["n_worsened_gt_2mm"] = int(np.sum(vv > 2))
                c_rows.append(row)
    hb_s = stats([p.get("heldout_bone_supplied_centroid_mm", float("nan")) for p in ok])["median"]
    hb_r = stats([p.get("heldout_bone_refined_centroid_mm", float("nan")) for p in ok])["median"]
    hb_d = stats([p.get("heldout_bone_diff_centroid_mm", float("nan")) for p in ok])["median"]
    hs_d = stats([p.get("heldout_soft_diff_centroid_mm", float("nan")) for p in ok])["median"]
    n_w2 = sum(1 for p in ok if p.get("heldout_bone_diff_centroid_mm", float("nan")) > WORSEN_MM)
    n_fail = sum(1 for p in ok if str(p.get("fallback_reason", "")).startswith("registration_failed"))
    crit = dict(
        C1=dict(text="lower cohort median held-out-bone centroid with refinement and median paired difference < 0",
                supplied_median=hb_s, refined_median=hb_r, median_paired_diff=hb_d, met=bool(hb_r < hb_s and hb_d < 0)),
        C2=dict(text="held-out soft tissue median paired difference <= +0.5 mm", median_paired_diff=hs_d,
                met=bool(hs_d <= 0.5)),
        C3=dict(text="held-out bone worsens > 2 mm in <= 5 patients (per-patient fallback applied)",
                n_worsened_gt_2mm=n_w2, n_registration_failed=n_fail, met=bool(n_w2 <= 5)),
        C4=dict(text="visual verification of overlays (author)", met="AUTHOR_TO_CONFIRM"))
    files = {}
    for name, rows in [("structure_level", s_rows), ("patient_level", p_rows), ("cohort_summary", c_rows)]:
        pth = out(name, "csv")
        write_csv(pth, rows)
        files[name] = pth
    with open(out("adoption_criteria", "json"), "w", encoding="utf-8") as fh:
        json.dump(crit, fh, indent=2)

    # ---------------- figures and overlays ----------------
    figs = []
    try:
        fig, axs = plt.subplots(1, 3, figsize=(16, 5.5))
        for ax, sname in zip(axs, ("heldout_bone", "heldout_soft", "head_neck")):
            x = [p.get("%s_supplied_centroid_mm" % sname, np.nan) for p in ok]
            y = [p.get("%s_refined_centroid_mm" % sname, np.nan) for p in ok]
            ax.scatter(x, y, s=14)
            lim = max([v for v in x + y if v == v] + [1]) * 1.05
            ax.plot([0, lim], [0, lim], "k--", lw=0.8)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_xlabel("supplied transform (mm)")
            ax.set_ylabel("refined transform (mm)")
            ax.set_title("%s: patient median centroid distance" % sname, fontsize=9)
        p = out("fig_supplied_vs_refined", "png")
        fig.tight_layout()
        fig.savefig(p, dpi=120)
        plt.close(fig)
        figs.append(p)
    except Exception as ex:
        print("figure error: %s" % ex)
    for _, cid in sorted(change_rank, reverse=True)[:N_OVERLAYS]:
        try:
            dP = os.path.join(DATA_ROOT, cid, "PSMA")
            dF = os.path.join(DATA_ROOT, cid, "FDG")
            Tp = load_tfm(os.path.join(dP, "rigid.tfm"))
            Tr = load_tfm(os.path.join(TFM_DIR, "Audit1C_c_tfm_%s_refined_resampling_%s.tfm" % (cid, TS)))
            ct = sitk.ReadImage(os.path.join(dP, "CT.nii.gz"))
            tsP = sitk.Cast(sitk.ReadImage(os.path.join(dP, "totseg_24.nii.gz")), sitk.sitkUInt8)
            tsF = sitk.Cast(sitk.ReadImage(os.path.join(dF, "totseg_24.nii.gz")), sitk.sitkUInt8)
            rs = sitk.Resample(tsF, tsP, Tp, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
            rr = sitk.Resample(tsF, tsP, Tr, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
            pr = next(p for p in ok if p["case_id"] == cid)
            pth = out("overlay_%s" % cid, "png")
            overlay(cid, ct, tsP, rs, rr, pth, "held-out bone %.1f -> %.1f mm; refinement change %.1f mm (criterion C4)" % (
                pr["heldout_bone_supplied_centroid_mm"], pr["heldout_bone_refined_centroid_mm"],
                pr["refinement_change_median_mm"]))
            figs.append(pth)
        except Exception as ex:
            print("overlay error %s: %s" % (cid, ex))

    # ---------------- console summary ----------------
    print("\n===== SUMMARY =====")
    print("Patients processed: %d/%d (errors: %d); per-patient fallback to supplied: %d" % (
        len(ok), len(p_rows), len(p_rows) - len(ok), sum(1 for p in ok if p.get("fallback_to_supplied"))))
    for r in c_rows:
        if r["metric"] == "centroid_mm":
            extra = ("| >5/10/15/20 mm: " + "/".join("%.0f%%" % r["pct_patients_gt_%dmm" % L] for L in CAND_LIMITS)
                     if "pct_patients_gt_5mm" in r else
                     "| improved>1mm %d, worsened>1mm %d, worsened>2mm %d" % (
                         r["n_improved_gt_1mm"], r["n_worsened_gt_1mm"], r["n_worsened_gt_2mm"]))
            print("%-13s %-28s n=%3d median %6.2f (IQR %.2f to %.2f, P90 %.2f, max %.2f) %s" % (
                r["set"], r["transform"], r["n"], r["median"], r["Q1"], r["Q3"], r["P90"], r["max"], extra))
    for k, v in crit.items():
        print("%s: %s -> %s" % (k, v["text"], v["met"]))
    print("Refined transforms saved in %s" % TFM_DIR)
    print("Elapsed: %.1f min" % ((time.time() - t0) / 60.0))

    man = dict(audit="Audit 1C-c - bone-based rigid refinement (Protocol Amendment 1) with held-out evaluation",
               protocol="PAPER3_PSMA_FDG_discordance_protocol_v2.5, Section 15, Amendment 1",
               script=os.path.abspath(__file__), timestamp=TS, cases=cases,
               python=platform.python_version(), simpleitk=sitk.Version.VersionString(), numpy=np.__version__,
               matplotlib=matplotlib.__version__, platform=platform.platform(),
               registration=dict(metric="MattesMI", bins=MI_BINS, sampling=SAMPLING, seed=SEED, shrink=SHRINK,
                                 smoothing_mm=SMOOTH, optimizer="RegularStepGradientDescent", learning_rate=LR,
                                 min_step=MIN_STEP, iterations=N_ITER, relaxation=RELAX,
                                 mask="fitting structures dilated %.0f mm" % MASK_DILATE_MM,
                                 init="supplied PSMA/rigid.tfm"),
               structure_sets=SETS, fallback_rule="held-out bone worsens > %.0f mm or registration failure" % WORSEN_MM,
               criteria=crit, outputs=dict(files, figures=figs, transforms_dir=TFM_DIR, run_log=log_path),
               note="Read-only on original data. Anatomy and coverage only. No tumour correspondence, DMTVF or "
                    "discordance results. Final adoption decision recorded by the author.")
    with open(out("manifest", "json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=2, default=str)
    print("Outputs written to %s (timestamp %s)" % (OUT_DIR, TS))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
