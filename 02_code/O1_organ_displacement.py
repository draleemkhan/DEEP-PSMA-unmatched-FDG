"""
O1_organ_displacement.py - PAPER 3, stage O1: soft-tissue (organ) displacement under the frozen Amendment 2
registration. CT and TotalSegmentator labels only; no PET, tumour mask or correspondence data are read.

Implements: StepSAP_soft_tissue_ST_spec_LOCKED_20261111.md (MD5 055199c4fd5d46039dcc65b6ab4f2416), SP1-SP3, and
PET_correspondence_implementation_spec_LOCKED_20261111.md (MD5 a5a2ec19d284deb44bc9ceebf8636658), I1.

Structures: liver, spleen, right kidney, left kidney, right lung (lobes 12+13+14), left lung (lobes 10+11).
Transform: Amendment 2 = supplied rigid.tfm corrected by the per-patient translation t (Audit 1C-e);
  supplied transform for train_0018, train_0049, train_0067 (fallbacks).
    point mapping FDG -> PSMA   T'(p) = FDG/rigid.tfm(p) - t
    resampling (sitk.Resample)  R'(q) = PSMA/rigid.tfm(q + t)
Metrics (functions copied VERBATIM from Audit1C_e_translation_refinement.py / Amendment 1):
  centroid distance = |T'(c_FDG) - c_PSMA| (LabelShapeStatistics centroids of the native label maps)
  ASSD  = Amendment 1 surface_metrics MSD (pooled both-direction contour-to-boundary distances; mean)
  HD95  = pooled 95th percentile
  FDG labels resampled onto the PSMA CT grid by nearest neighbour with R'. Lungs: lobe labels merged into one
  label on each native label map before any computation.
Eligibility (SP3): PRESENT_COMPLETE and frac_largest_component >= 0.95 on BOTH CTs (Audit 1B); for a lung, every
  contributing lobe. The union need not be one connected component.
Validation (I1): liver, spleen, kidneys must reproduce Audit 1C-e (<= 0.01 mm) for every record with identical
  masks, transform, eligibility and method: corrected_* values for patients whose adopted transform is the
  correction, supplied_* values for the three fallbacks. Records not comparable are listed with the reason.
  Any unexplained difference STOPS the run (no output table is accepted).

Usage (folder 02_code; PowerShell):
  & "<totalseg python>" O1_organ_displacement.py --selftest
  & "<totalseg python>" O1_organ_displacement.py --preflight     2 patients (same TEST patients as before)
  & "<totalseg python>" O1_organ_displacement.py                 all 100
Outputs (06_output): O1_[PREFLIGHT_]organ_displacement_<run>.csv, O1_[PREFLIGHT_]validation_<run>.csv,
                     O1_[PREFLIGHT_]manifest_<run>.json
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

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
FALLBACKS = {"train_0018", "train_0049", "train_0067"}
INPUTS = {   # (folder, file): MD5
    (OUT_DIR, "Audit1C_e_per_patient_20261031_141128.csv"): "c66da788a6c73737163b4c09ed76f494",
    (OUT_DIR, "Audit1C_e_structure_level_20261031_141128.csv"): "289b57e69a6f119968b3ae3990db9fd9",
    (OUT_DIR, "Audit1B_label_level_20261029_134401.csv"): "896a1681e90da4fe793a3a449d737fde",
    (OUT_DIR, "A3_visualQC_KEY_do_not_open_20261105.csv"): "c9a4548d3daf13d3e3dadc96171a9a94",   # preflight IDs only
    (PROTO_DIR, "StepSAP_soft_tissue_ST_spec_LOCKED_20261111.md"): "055199c4fd5d46039dcc65b6ab4f2416",
    (PROTO_DIR, "PET_correspondence_implementation_spec_LOCKED_20261111.md"): "a5a2ec19d284deb44bc9ceebf8636658",
}
STRUCT = {"liver": [5], "spleen": [1], "kidney_right": [2], "kidney_left": [3],
          "lung_right": [12, 13, 14], "lung_left": [10, 11]}
NEW_LABEL = {"lung_right": 201, "lung_left": 202}          # merged lobe labels (not used by TotalSegmentator)
VALIDATE = ["liver", "spleen", "kidney_right", "kidney_left"]
FRAG_MIN = 0.95
TOL_MM = 0.01


# ------------------------------ verbatim helpers (Audit1C_e_translation_refinement.py) ------------------------------
def centroids(label_img):
    import SimpleITK as sitk
    f = sitk.LabelShapeStatisticsImageFilter()
    f.Execute(label_img)
    return {int(k): np.array(f.GetCentroid(k)) for k in f.GetLabels()}, f


def corrected_transforms(T_point, R_resample, t):
    """T'(p) = T(p) - t (FDG->PSMA points); R'(q) = R(q + t) (resampling). ITK composite: last added applied first."""
    import SimpleITK as sitk
    Tc = sitk.CompositeTransform(3)
    Tc.AddTransform(sitk.TranslationTransform(3, tuple(float(-x) for x in t)))
    Tc.AddTransform(T_point)                     # applied first
    Rc = sitk.CompositeTransform(3)
    Rc.AddTransform(R_resample)
    Rc.AddTransform(sitk.TranslationTransform(3, tuple(float(x) for x in t)))   # applied first
    return Tc, Rc


def tp(T, p):
    return np.array(T.TransformPoint(tuple(float(x) for x in p)))


def surface_metrics(imgA, imgB, l, bbA, bbB, margin=3):
    """VERBATIM from Amendment 1 (02_code\\Audit1C_c_rigid_refinement.py, MD5 0c350ab57b59c1dea8fc729b434457b7):
    region of interest = union of the two label bounding boxes + 3-voxel margin (clipped to the image);
    pooled symmetric contour-to-boundary distances; MSD = mean, HD95 = 95th percentile."""
    import SimpleITK as sitk
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


def surface_count(imgA, imgB, l, bbA, bbB, margin=3):
    """number of pooled surface distances used by surface_metrics (same ROI and contours; reporting only)."""
    import SimpleITK as sitk
    sz = imgA.GetSize()
    lo = [int(max(0, min(bbA[i], bbB[i]) - margin)) for i in range(3)]
    hi = [int(min(sz[i], max(bbA[i] + bbA[i + 3], bbB[i] + bbB[i + 3]) + margin)) for i in range(3)]
    size = [hi[i] - lo[i] for i in range(3)]
    n = 0
    for img in (imgA, imgB):
        x = sitk.BinaryThreshold(sitk.RegionOfInterest(img, size, lo), l, l, 1, 0)
        n += int((sitk.GetArrayFromImage(sitk.BinaryContour(x, False, 0, 1)) > 0).sum())
    return n


# ------------------------------ O1 helpers (tested) ------------------------------
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def eligibility(cid, labels, L1B):
    """SP3: PRESENT_COMPLETE and frac_largest >= 0.95 on both CTs for every contributing label."""
    why = []
    for tr in ("PSMA", "FDG"):
        for lab in labels:
            try:
                r = L1B.loc[(cid, tr, lab)]
                cat, fr = str(r["category"]), float(r["frac_largest_component"])
            except KeyError:
                cat, fr = "MISSING_IN_AUDIT1B", np.nan
            if cat != "PRESENT_COMPLETE":
                why.append(f"{tr} label {lab}: {cat}")
            elif not (fr >= FRAG_MIN):
                why.append(f"{tr} label {lab}: fragmented (frac_largest {fr:.3f})")
    return (len(why) == 0), "; ".join(why)


def relabel_lungs(arr):
    """merge lobe labels into one label per lung (native label array; other labels unchanged)."""
    out = arr.astype(np.uint8, copy=True)            # 201/202 fit in uint8 (as the original label maps)
    for name, new in NEW_LABEL.items():
        out[np.isin(arr, STRUCT[name])] = new
    return out


def label_of(name):
    return NEW_LABEL.get(name, STRUCT[name][0])


def preflight_cases(cases, oldkey):
    return [c for c in cases if c not in set(oldkey["case_id"])][:2]


def compare(new, ref, tol=TOL_MM):
    return bool(np.isfinite(new) and np.isfinite(ref) and abs(new - ref) <= tol)


# ------------------------------ per-patient processing ------------------------------
def process(cid, t, L1B):
    import SimpleITK as sitk
    dP, dF = DATA_ROOT / cid / "PSMA", DATA_ROOT / cid / "FDG"
    R = sitk.ReadTransform(str(dP / "rigid.tfm")); T = sitk.ReadTransform(str(dF / "rigid.tfm"))
    Tc, Rc = corrected_transforms(T, R, t)
    rawP = sitk.ReadImage(str(dP / "totseg_24.nii.gz"), sitk.sitkUInt8)
    rawF = sitk.ReadImage(str(dF / "totseg_24.nii.gz"), sitk.sitkUInt8)

    def merged(img):
        m = sitk.GetImageFromArray(relabel_lungs(sitk.GetArrayFromImage(img)))
        m.CopyInformation(img)
        return m
    segP, segF = merged(rawP), merged(rawF)
    cenP, shpP = centroids(segP); cenF, shpF = centroids(segF)
    res = sitk.Resample(segF, segP, Rc, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    shpR = sitk.LabelShapeStatisticsImageFilter(); shpR.Execute(res)
    # composition check: R'(T'(p)) = p at test points
    worst = 0.0
    for p in np.random.default_rng(3).uniform(-200, 200, (5, 3)):
        worst = max(worst, float(np.linalg.norm(tp(Rc, tp(Tc, p)) - p)))
    rows = []
    for name, labs in STRUCT.items():
        lab = label_of(name)
        ok, why = eligibility(cid, labs, L1B)
        row = dict(case_id=cid, structure=name, labels="+".join(map(str, labs)), eligible=ok, ineligible_reason=why,
                   centroid_mm=np.nan, ASSD_mm=np.nan, HD95_mm=np.nan, n_surface_distances=0)
        if ok:
            if lab not in cenP or lab not in cenF or not shpR.HasLabel(lab):
                row.update(eligible=False, ineligible_reason="label absent after merge/resampling")
            else:
                row["centroid_mm"] = float(np.linalg.norm(tp(Tc, cenF[lab]) - cenP[lab]))
                m, h = surface_metrics(segP, res, lab, shpP.GetBoundingBox(lab), shpR.GetBoundingBox(lab))
                n = surface_count(segP, res, lab, shpP.GetBoundingBox(lab), shpR.GetBoundingBox(lab))
                row.update(ASSD_mm=m, HD95_mm=h, n_surface_distances=n)
        rows.append(row)
    return rows, worst


# ------------------------------ self-test (synthetic) ------------------------------
def selftest():
    import SimpleITK as sitk
    res = []
    ok = lambda n, c: res.append((n, bool(c)))
    # relabel
    a = np.zeros((4, 4, 4), np.uint8); a[0, 0, 0] = 12; a[1, 1, 1] = 14; a[2, 2, 2] = 10; a[3, 3, 3] = 5
    r = relabel_lungs(a)
    ok("lobes merged per side; other labels kept", r[0, 0, 0] == 201 and r[1, 1, 1] == 201 and r[2, 2, 2] == 202
       and r[3, 3, 3] == 5)
    # eligibility
    L = pd.DataFrame([dict(case_id="c", tracer=tr, label_id=l, category="PRESENT_COMPLETE", frac_largest_component=1.0)
                      for tr in ("PSMA", "FDG") for l in (1, 2, 3, 5, 10, 11, 12, 13, 14)]).set_index(["case_id", "tracer", "label_id"])
    ok("eligible when complete on both", eligibility("c", [12, 13, 14], L)[0])
    L2 = L.copy(); L2.loc[("c", "FDG", 13), "frac_largest_component"] = 0.9
    ok("one fragmented lobe -> lung ineligible", not eligibility("c", [12, 13, 14], L2)[0])
    L3 = L.copy(); L3.loc[("c", "PSMA", 5), "category"] = "PRESENT_BOUNDARY_Z"
    ok("truncated on one CT -> ineligible", not eligibility("c", [5], L3)[0])
    ok("missing in Audit 1B -> ineligible", not eligibility("x", [5], L)[0])
    # transforms: identity T, R and translation t
    T = sitk.Euler3DTransform(); R = sitk.Euler3DTransform()
    t = np.array([2.0, -3.0, 4.0])
    Tc, Rc = corrected_transforms(T, R, t)
    p = np.array([10.0, 20.0, 30.0])
    ok("T'(p) = p - t", np.allclose(tp(Tc, p), p - t))
    ok("R'(q) = q + t", np.allclose(tp(Rc, p), p + t))
    ok("R'(T'(p)) = p", np.allclose(tp(Rc, tp(Tc, p)), p))
    # metrics: identical masks -> 0; shifted box -> known centroid and surface distances
    arr = np.zeros((40, 40, 40), np.uint8); arr[10:20, 10:20, 10:20] = 5
    img = sitk.GetImageFromArray(arr); img.SetSpacing((2.0, 2.0, 3.0))
    sh = sitk.LabelShapeStatisticsImageFilter(); sh.Execute(img)
    m, h = surface_metrics(img, img, 5, sh.GetBoundingBox(5), sh.GetBoundingBox(5))
    n = surface_count(img, img, 5, sh.GetBoundingBox(5), sh.GetBoundingBox(5))
    ok("identical masks: ASSD = HD95 = 0", m == 0 and h == 0 and n > 0)
    arr2 = np.zeros_like(arr); arr2[10:20, 10:20, 12:22] = 5          # shift 2 voxels in x = 4 mm
    img2 = sitk.GetImageFromArray(arr2); img2.CopyInformation(img)
    sh2 = sitk.LabelShapeStatisticsImageFilter(); sh2.Execute(img2)
    c1, _ = centroids(img); c2, _ = centroids(img2)
    ok("centroid shift = 4 mm", abs(np.linalg.norm(c2[5] - c1[5]) - 4.0) < 1e-9)
    m2, h2 = surface_metrics(img, img2, 5, sh.GetBoundingBox(5), sh2.GetBoundingBox(5))
    ok("shifted box: 0 < ASSD <= 4 mm, HD95 <= 4 mm", 0 < m2 <= 4.0 + 1e-9 and h2 <= 4.0 + 1e-9)
    # resampling with R'(q) = q + t moves a label by -t in PSMA space (NN; integer-voxel t)
    t2 = np.array([4.0, 0.0, 0.0])
    _, Rc2 = corrected_transforms(T, R, t2)
    rs = sitk.Resample(img, img, Rc2, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    cr, _ = centroids(rs)
    ok("resampled label shifted by -t", np.allclose(cr[5] - c1[5], -t2, atol=1e-9))
    # comparison rule
    ok("compare within 0.01 mm", compare(1.004, 1.0) and not compare(1.02, 1.0) and not compare(np.nan, 1.0))
    ok("fallback set", FALLBACKS == {"train_0018", "train_0049", "train_0067"})
    for n_, g in res:
        print(f"  [{'PASS' if g else 'FAIL'}] {n_}")
    print(f"selftest: {sum(g for _, g in res)}/{len(res)} passed")
    return all(g for _, g in res)


# ------------------------------ run ------------------------------
def run(preflight):
    import SimpleITK as sitk
    t0 = time.time(); run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = "O1_PREFLIGHT_" if preflight else "O1_"
    md = {name: md5(d / name) if (d / name).exists() else "MISSING" for (d, name) in INPUTS}
    bad = {n: md[n] for (d, n), e in INPUTS.items() if md[n] != e}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(INPUTS)}/{len(INPUTS)} match")
    PP = pd.read_csv(OUT_DIR / "Audit1C_e_per_patient_20261031_141128.csv").set_index("case_id")
    SL = pd.read_csv(OUT_DIR / "Audit1C_e_structure_level_20261031_141128.csv")
    L1B = pd.read_csv(OUT_DIR / "Audit1B_label_level_20261029_134401.csv").set_index(["case_id", "tracer", "label_id"])
    cases = sorted(PP.index)
    if len(cases) != 100:
        sys.exit("STOP: per-patient table does not list 100 patients")
    if preflight:
        cases = preflight_cases(cases, pd.read_csv(OUT_DIR / "A3_visualQC_KEY_do_not_open_20261105.csv"))
    rows, comp = [], {}
    for i, cid in enumerate(cases, 1):
        t = np.zeros(3) if cid in FALLBACKS else PP.loc[cid, ["t_x_LPS", "t_y_LPS", "t_z_LPS"]].to_numpy(float)
        r, w = process(cid, t, L1B)
        comp[cid] = w
        if w > 1e-6:
            sys.exit(f"STOP: composition check failed for {cid} ({w:.2e} mm)")
        rows += r
        print(f"  [{i:3d}/{len(cases)}] {cid}: eligible {sum(x['eligible'] for x in r)}/6  ({time.time() - t0:.0f} s)")
    D = pd.DataFrame(rows)
    # ---- validation against Audit 1C-e ----
    V = []
    for _, r in D[D.structure.isin(VALIDATE)].iterrows():
        ref = SL[(SL.case_id == r.case_id) & (SL.structure == r.structure)]
        pre = "supplied" if r.case_id in FALLBACKS else "corrected"
        if len(ref) != 1:
            V.append(dict(case_id=r.case_id, structure=r.structure, compared=False,
                          reason="no Audit 1C-e record (not evaluable in Amendment 1)", PASS=np.nan)); continue
        ref = ref.iloc[0]
        if not r.eligible:
            V.append(dict(case_id=r.case_id, structure=r.structure, compared=False,
                          reason=f"not eligible under SP3 ({r.ineligible_reason})", PASS=np.nan)); continue
        d = {m: abs(r[f"{m}_mm"] - ref[f"{pre}_{k}_mm"]) for m, k in (("centroid", "centroid"), ("ASSD", "MSD"), ("HD95", "HD95"))}
        V.append(dict(case_id=r.case_id, structure=r.structure, compared=True, reference=pre,
                      diff_centroid_mm=d["centroid"], diff_ASSD_mm=d["ASSD"], diff_HD95_mm=d["HD95"],
                      PASS=bool(all(compare(r[f"{m}_mm"], ref[f"{pre}_{k}_mm"]) for m, k in
                                    (("centroid", "centroid"), ("ASSD", "MSD"), ("HD95", "HD95"))))))
    V = pd.DataFrame(V)
    vp = OUT_DIR / f"{tag}validation_{run_id}.csv"; V.to_csv(vp, index=False)
    cmp_ = V[V.compared == True]
    nfail = int((cmp_.PASS == False).sum())
    print(f"Validation vs Audit 1C-e: {len(cmp_) - nfail}/{len(cmp_)} compared records reproduced (<= {TOL_MM} mm); "
          f"not compared: {int((V.compared == False).sum())} (reasons in {vp.name})")
    if nfail:
        sys.exit(f"STOP: {nfail} records do not reproduce Audit 1C-e; investigate before use (see {vp.name}); "
                 f"no displacement table written")
    dp = OUT_DIR / f"{tag}organ_displacement_{run_id}.csv"; D.to_csv(dp, index=False)
    elig = D.groupby("structure").eligible.sum().to_dict()
    print(f"Eligible per structure: {elig}")
    man = dict(stage="O1 organ displacement (CT only)", mode="PREFLIGHT" if preflight else "FULL", run=run_id,
               script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs=md, n_patients=len(cases),
               composition_max_mm=max(comp.values()), outputs={p.name: md5(p) for p in (dp, vp)},
               validation=dict(compared=int(len(cmp_)), reproduced=int(len(cmp_) - nfail),
                               not_compared=int((V.compared == False).sum())),
               params=dict(FRAG_MIN=FRAG_MIN, TOL_MM=TOL_MM, STRUCT=STRUCT, NEW_LABEL=NEW_LABEL),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             SimpleITK=sitk.Version_VersionString()),
               runtime_s=round(time.time() - t0), status="COMPLETE")
    mp = OUT_DIR / f"{tag}manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"done: {dp.name}; manifest {mp.name} MD5 {md5(mp)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.preflight)
