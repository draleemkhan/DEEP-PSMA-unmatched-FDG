"""
SAP_R1_region_rules.py - PAPER 3, Step 2 SAP section R.1 (draft MD5 44f039f2696394e8aaa95bf472c40aa0):
region rules for registration-QC V_I, their synthetic self-tests (R1.7 tests 1-6) and the plane table (R1.7 test 7).
No PET, lesion or correspondence data are read.

  python SAP_R1_region_rules.py --selftest     R1.7 tests 1-6 (synthetic; no patient data)
  python SAP_R1_region_rules.py --planes       R1.7 test 7: per-patient plane table (PSMA labels only)

--planes recomputes, for all 100 patients, the boundaries with the LOCKED Audit 1C-f rule (Audit1C_f_vascular_muscular_
QC.py, function `boundaries`: midpoint of adjacent vertebral z-centroids; unreliable vertebrae (not PRESENT_COMPLETE on
the PSMA CT, frac_largest < 0.95 or adjacent to such a vertebra, Audit 1B flagged levels) replaced by linear
interpolation in level position between the nearest reliable vertebrae; no extrapolation), extended to the C2/C3
boundary (level position 22.5). Stops unless every C7T1/T12L1/L5S1 value and status reproduces the Audit 1C-f log
(<= 0.01 mm), including the SUBSTITUTED and UNDEFINED rows. Reports the C2/C3 status for every patient and its
difference from the plain C2/C3 midpoint used by the trunk QC.
Output (06_output): SAP_R1_plane_table_<run>.csv (+ MD5 printed).
"""
import argparse
import hashlib
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

BAND_MM = 15.0
FAILED_NECK = {"train_0001", "train_0002", "train_0010", "train_0012", "train_0015", "train_0023", "train_0037",
               "train_0042", "train_0055", "train_0058", "train_0062", "train_0078", "train_0080", "train_0089"}
FAILED_HEAD = {"train_0007", "train_0012", "train_0023", "train_0028", "train_0066", "train_0068", "train_0089",
               "train_0010", "train_0033", "train_0042", "train_0053", "train_0079", "train_0095"}
FAILED_PELVIS = {"train_0069", "train_0089"}


# ------------------------------ region rules (R.3 + R1.3) ------------------------------
def region_of(z, P):
    """z: array of PSMA z (mm). P: dict z_head, z_c7t1, z_pelvis (nan = undefined). Returns array of region names.
    head z > z_head; neck z_c7t1 <= z <= z_head (both ends closed); pelvis z < z_pelvis; otherwise thorax/abdomen."""
    z = np.asarray(z, float)
    out = np.full(z.shape, "thorax_abdomen", dtype=object)
    if np.isfinite(P["z_pelvis"]):
        out[z < P["z_pelvis"]] = "pelvis"
    out[(z >= P["z_c7t1"]) & (z <= P["z_head"])] = "neck"
    out[z > P["z_head"]] = "head"
    return out


def r6_region(z, P):
    """R6 transform regions only (head / trunk / pelvis); tie on a plane -> trunk."""
    r = region_of(np.atleast_1d(z), P)
    return np.where(np.isin(r, ["neck", "thorax_abdomen"]), "trunk", r)


def relevant_transforms(ta_points, P, adopted):
    """SAP R.5 step 2. ta_points: (N,3) voxel centres mapped with T_A. adopted: {'head': bool, 'pelvis': bool} (final
    transform differs from T_A). Returns the list of transform names to apply: 'TA' always; the final transform of the
    R6 region of the T_A centroid; the adjacent region's transform if in the 15 mm band of an R6 plane or crossing it."""
    T = ["TA"]
    zc = float(ta_points[:, 2].mean())
    reg = r6_region(zc, P)[0]
    if reg in ("head", "pelvis") and adopted[reg]:
        T.append(reg)
    for plane, side in (("z_head", "head"), ("z_pelvis", "pelvis")):
        zp = P[plane]
        if not np.isfinite(zp):
            continue
        z = ta_points[:, 2]
        crosses = (z.min() < zp) and (z.max() > zp)
        if (abs(zc - zp) <= BAND_MM or crosses) and adopted[side] and side not in T:
            T.append(side)
    return T


def registration_vi(mapped, P, failed):
    """mapped: {transform name: (N,3) mapped points}. failed: set of failed region names for this patient
    ('head', 'neck', 'pelvis'). True -> the whole component is V_I (registration QC), counted once."""
    for pts in mapped.values():
        if np.isin(region_of(pts[:, 2], P), sorted(failed)).any():
            return True
    return False


def failed_regions(cid):
    return {r for r, s in (("head", FAILED_HEAD), ("neck", FAILED_NECK), ("pelvis", FAILED_PELVIS)) if cid in s}


# ------------------------------ R1.7 tests 1-6 ------------------------------
def selftest():
    res = []
    ok = lambda n, c: res.append((n, bool(c)))
    P = dict(z_head=100.0, z_c7t1=40.0, z_pelvis=-300.0)
    # 1 tie rules
    ok("1a point on z_head -> neck", region_of([100.0], P)[0] == "neck")
    ok("1b point on z_c7t1 -> neck", region_of([40.0], P)[0] == "neck")
    ok("1c just above z_head -> head", region_of([100.0 + 1e-9], P)[0] == "head")
    ok("1d just below z_c7t1 -> thorax/abdomen", region_of([40.0 - 1e-9], P)[0] == "thorax_abdomen")
    ok("1e point on z_pelvis -> trunk (R.3)", r6_region(-300.0, P)[0] == "trunk")
    box = lambda z0, z1: np.array([[x, y, z] for x in (0, 5) for y in (0, 5) for z in np.linspace(z0, z1, 9)])
    nofail = dict(head=False, pelvis=False)
    # 2 wholly inside a failed neck
    c = box(55, 85)
    ok("2 inside failed neck -> V_I", registration_vi({"TA": c}, P, {"neck"}))
    ok("2b inside accepted neck -> not V_I", not registration_vi({"TA": c}, P, set()))
    # 3 straddling C7/T1
    c = box(30, 50)
    ok("3a straddles C7/T1, neck failed -> V_I", registration_vi({"TA": c}, P, {"neck"}))
    ok("3b straddles C7/T1, neck accepted -> not V_I", not registration_vi({"TA": c}, P, set()))
    ok("3c no band at C7/T1 (TA only)", relevant_transforms(c, P, dict(head=True, pelvis=True)) == ["TA"])
    # 4 in the C2/C3 band; head transform maps it down into the failed neck
    c = box(104, 112)                                   # T_A: head region, centroid 108 mm (8 mm above z_head)
    T = relevant_transforms(c, P, dict(head=True, pelvis=False))
    ok("4a band at C2/C3 -> TA + head", T == ["TA", "head"])
    mapped = {"TA": c, "head": c - np.array([0, 0, 10.0])}   # head transform shifts 10 mm caudally
    ok("4b head-transform points land in failed neck -> V_I", registration_vi(mapped, P, {"neck"}))
    ok("4c TA alone would not touch the neck", not registration_vi({"TA": c}, P, {"neck"}))
    ok("4d band not used if head not adopted", relevant_transforms(c, P, nofail) == ["TA"])
    # 5 failed head and failed neck -> V_I once
    c = box(80, 130)
    ok("5 head + neck failed -> one V_I decision", registration_vi({"TA": c}, P, {"head", "neck"}) is True)
    ok("5b failed_regions(train_0089) = head, neck, pelvis", failed_regions("train_0089") == {"head", "neck", "pelvis"})
    ok("5c 29 failed regions in 23 patients",
       len(FAILED_HEAD) + len(FAILED_NECK) + len(FAILED_PELVIS) == 29
       and len(FAILED_HEAD | FAILED_NECK | FAILED_PELVIS) == 23)
    # 6 z_pelvis undefined -> thorax/abdomen to the end of the FOV; no pelvic band
    Q = dict(P, z_pelvis=np.nan)
    ok("6a undefined L5/S1 -> thorax/abdomen at -900 mm", region_of([-900.0], Q)[0] == "thorax_abdomen")
    ok("6b undefined L5/S1 -> no pelvic transform", relevant_transforms(box(-310, -290), Q, dict(head=True, pelvis=True)) == ["TA"])
    ok("6c defined L5/S1 band -> TA + pelvis", relevant_transforms(box(-310, -290), P, dict(head=False, pelvis=True)) == ["TA", "pelvis"])
    for n, g in res:
        print(f"  [{'PASS' if g else 'FAIL'}] {n}")
    print(f"selftest: {sum(g for _, g in res)}/{len(res)} passed")
    return all(g for _, g in res)


# ------------------------------ R1.7 test 7: plane table (PC) ------------------------------
POS = {"C2C3": 22.5, "C7T1": 17.5, "T12L1": 5.5, "L5S1": 0.5}       # Audit 1C-f level positions (S1 = 0 ... C1 = 24)
VERT_LABELS = list(range(26, 51))                                     # label = 26 + position
FLAGGED = {"train_0088": [39, 40, 41, 42, 43], "train_0059": [26, 27, 28, 29, 30], "train_0029": [29, 30, 31, 32, 33]}
FRAG_MIN = 0.95


def boundaries_1cf(cid, cenz, L1B):
    """identical to Audit1C_f_vascular_muscular_QC.boundaries (locked), with C2/C3 added."""
    def cat(lab):
        try:
            r = L1B.loc[(cid, "PSMA", lab)]
            return str(r["category"]), float(r["frac_largest_component"])
        except KeyError:
            return "MISSING_IN_AUDIT1B", np.nan
    bad = set(FLAGGED.get(cid, []))
    for lab in VERT_LABELS:
        c, f = cat(lab)
        if f == f and f < FRAG_MIN:
            bad.update([lab - 1, lab, lab + 1])
    reliable = {lab - 26: cenz[lab] for lab in VERT_LABELS
                if cat(lab)[0] == "PRESENT_COMPLETE" and lab not in bad and lab in cenz}
    out = {}
    for name, pos in POS.items():
        below = [p for p in reliable if p < pos]; above = [p for p in reliable if p > pos]
        if not below or not above:
            out[name] = (np.nan, "UNDEFINED", ""); continue
        pa, pb = max(below), min(above)
        za, zb = reliable[pa], reliable[pb]
        sub = not (pa == int(np.floor(pos)) and pb == int(np.ceil(pos)))
        out[name] = (za + (pos - pa) * (zb - za) / (pb - pa), "SUBSTITUTED" if sub else "DIRECT",
                     f"vertebrae {26 + pa}/{26 + pb} (positions {pa}/{pb})")
    return out


def planes():
    import SimpleITK as sitk
    import importlib.util
    here = Path(__file__).resolve().parent
    sp = importlib.util.spec_from_file_location("a3", here / "Amendment3_regional_rigid.py")
    a3 = importlib.util.module_from_spec(sp); sp.loader.exec_module(a3)
    O = a3.OUT_DIR
    exp = {here / "Amendment3_regional_rigid.py": "b21a405145526d56a59bb5950ef53acf",
           here / "Audit1C_f_vascular_muscular_QC.py": "e893c4ff12dd871c52d2c992f8c2b7b3",
           O / "Audit1C_f_boundaries_substitution_log_20261031_164711.csv": "14b882a1ed3c0949c099a26291bc02c5",
           O / "Amendment3_fit_20261102_073040.csv": "dc74686ed0b024d8524668f9b9e16c94"}
    bad = {p.name: a3.md5(p) for p, e in exp.items() if not p.exists() or a3.md5(p) != e}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    l1b = O / "Audit1B_label_level_20261029_134401.csv"
    print(f"C1 MD5: {len(exp)}/{len(exp)} match; Audit 1B label table MD5 {a3.md5(l1b)} (recorded)")
    L1B = pd.read_csv(l1b).set_index(["case_id", "tracer", "label_id"])
    LOG = pd.read_csv(O / "Audit1C_f_boundaries_substitution_log_20261031_164711.csv")
    cases = sorted(pd.read_csv(O / "Amendment3_fit_20261102_073040.csv").case_id.unique())
    rows, worst, stop = [], 0.0, []
    for cid in cases:
        im = sitk.ReadImage(str(a3.DATA_ROOT / cid / "PSMA" / "totseg_24.nii.gz"))
        sh = sitk.LabelShapeStatisticsImageFilter(); sh.Execute(sitk.Cast(im, sitk.sitkUInt8))
        cen = {int(l): np.array(sh.GetCentroid(l)) for l in sh.GetLabels()}
        B = boundaries_1cf(cid, {l: c[2] for l, c in cen.items()}, L1B)
        for name in ("C7T1", "T12L1", "L5S1"):
            lg = LOG[(LOG.case_id == cid) & (LOG.boundary == name)].iloc[0]
            z, st, _ = B[name]
            same = st == lg.status and ((st == "UNDEFINED" and not np.isfinite(lg.z_mm))
                                        or (st != "UNDEFINED" and abs(z - lg.z_mm) <= 0.01))
            if st != "UNDEFINED":
                worst = max(worst, abs(z - lg.z_mm))
            if not same:
                stop.append(f"{cid} {name}")
        plain = (cen[49][2] + cen[48][2]) / 2 if 49 in cen and 48 in cen else np.nan
        rows.append(dict(case_id=cid, z_head_C2C3=B["C2C3"][0], C2C3_status=B["C2C3"][1], C2C3_detail=B["C2C3"][2],
                         C2C3_plain_midpoint=plain, C2C3_rule_minus_plain_mm=B["C2C3"][0] - plain,
                         z_C7T1=B["C7T1"][0], C7T1_status=B["C7T1"][1], C7T1_detail=B["C7T1"][2],
                         z_T12L1=B["T12L1"][0], T12L1_status=B["T12L1"][1],
                         z_pelvis_L5S1=B["L5S1"][0], L5S1_status=B["L5S1"][1]))
    T = pd.DataFrame(rows)
    print(f"Audit 1C-f log reproduced: {300 - len(stop)}/300 boundaries (max |diff| {worst:.2e} mm)")
    print(f"C2/C3 status: {T.C2C3_status.value_counts().to_dict()}; max |rule - plain midpoint| "
          f"{np.nanmax(np.abs(T.C2C3_rule_minus_plain_mm)):.2e} mm")
    sub = T[T.C7T1_status == "SUBSTITUTED"][["case_id", "C7T1_detail", "z_C7T1"]]
    print("C7/T1 substituted:\n" + sub.to_string(index=False))
    if stop:
        sys.exit(f"STOP test 7: Audit 1C-f log not reproduced for {stop}; no table written")
    p = O / f"SAP_R1_plane_table_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    T.to_csv(p, index=False)
    print(f"written {p.name}  MD5 {hashlib.md5(p.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--planes", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.planes:
        planes()
