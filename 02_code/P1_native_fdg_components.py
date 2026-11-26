"""
P1_native_fdg_components.py - PAPER 3, stage P1: native FDG tumour components and eligibility.
Native FDG space ONLY: reads FDG TTB and FDG PET; no PSMA data, no transform, no matching, no discordance.

Implements: consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) section 2 (items 1-2) and
PET_correspondence_implementation_spec_LOCKED_20261111.md (MD5 a5a2ec19d284deb44bc9ceebf8636658) I0, I2.

Per patient (FDG):
  BW = TTB > 0 (as Audit 1D); connected components for 6, 18 and 26 connectivity (scipy.ndimage.label with the
  standard 3-D structuring elements); component IDs in raster order of first voxel, which equals MATLAB
  bwconncomp's column-major order on the NIfTI (i, j, k) array; voxel volume = |det(direction x spacing)| (mL).
  Per component: voxels, volume_mL, SUVmax, SUVmean (native FDG PET), centroid (LPS and RAS, mm), world-axis extent
  (mm, voxel width included, as Audit 1D), eligibility rule A (>= 0.1 mL AND >= 3 voxels) and rule B
  (>= 0.25 mL AND >= 3 voxels).
  Per patient: V_FDG,total (whole native TTB volume; identical for every connectivity), V_S under A and B.
Acceptance (any failure STOPS; no table accepted):
  - components per patient and connectivity identical to Audit 1D component level; voxels and volume per
    component ID identical (volume within 1e-6 mL); SUVmax/SUVmean within 1e-6; centroid within 1e-4 mm (RAS)
  - reconciliation: sum of component volumes = V_FDG,total (relative 1e-9), identical across connectivities
  - full run only: excluded components and pooled excluded volume reproduce Option D (26-conn A 871 / 0.24 %,
    B 1,389 / 0.61 %; 18-conn A 933 / 0.25 %, B 1,466 / 0.63 %; 6-conn A 1,428 / 0.36 %, B 1,997 / 0.78 %)
Outputs (06_output\\P1_[PREFLIGHT_]<run>\\): native component label images P1_<case>_FDG_conn<c>.nii.gz (int32,
  FDG PET grid), P1_components_<run>.csv, P1_patients_<run>.csv, P1_audit1D_comparison_<run>.csv, manifest.

Usage (folder 02_code):
  & "<totalseg python>" -u P1_native_fdg_components.py --selftest
  & "<totalseg python>" -u P1_native_fdg_components.py --preflight
  & "<totalseg python>" -u P1_native_fdg_components.py
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

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
INPUTS = {
    (OUT_DIR, "Audit1D_component_level_20261030_215034.csv"): "30f5f4614a345253a850c737bb37c1aa",
    (OUT_DIR, "Audit1D_scan_level_20261030_215034.csv"): "187366d01501ef79e7116205ef8e78b5",
    (OUT_DIR, "A3_visualQC_KEY_do_not_open_20261105.csv"): "c9a4548d3daf13d3e3dadc96171a9a94",   # preflight IDs only
    (PROTO_DIR, "StepSAP_consolidated_LOCKED_20261111.md"): "0c64711db1207e31a70cf3f558ccbd7d",
    (PROTO_DIR, "PET_correspondence_implementation_spec_LOCKED_20261111.md"): "a5a2ec19d284deb44bc9ceebf8636658",
}
CONNS = (6, 18, 26)
RULES = {"A": (0.1, 3), "B": (0.25, 3)}      # eligible if volume_mL >= v AND voxels >= n
OPTION_D = {(26, "A"): (871, 0.24), (26, "B"): (1389, 0.61), (18, "A"): (933, 0.25), (18, "B"): (1466, 0.63),
            (6, "A"): (1428, 0.36), (6, "B"): (1997, 0.78)}       # (excluded components, pooled % rounded to 0.01)
TOL_ML, TOL_SUV, TOL_MM, TOL_REL = 1e-6, 1e-6, 1e-4, 1e-9


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ core (tested) ------------------------------
def structure(conn):
    return ndimage.generate_binary_structure(3, {6: 1, 18: 2, 26: 3}[conn])


def label_components(bw_kji, conn):
    """bw_kji: boolean array in SimpleITK order (k, j, i). Raster (C) order over (k, j, i) equals MATLAB
    column-major order over (i, j, k), so labels are numbered as bwconncomp numbers its objects."""
    lab, n = ndimage.label(bw_kji, structure=structure(conn))
    return lab.astype(np.int32), int(n)


def eligible(volume_ml, voxels, rule):
    v, n = RULES[rule]
    return (volume_ml >= v) & (voxels >= n)


def phys_lps(img, kji):
    D = np.array(img.GetDirection()).reshape(3, 3); S = np.array(img.GetSpacing())
    return np.array(img.GetOrigin()) + (kji[:, ::-1].astype(float) * S) @ D.T


def preflight_cases(cases, oldkey):
    return [c for c in cases if c not in set(oldkey["case_id"])][:2]


# ------------------------------ per patient ------------------------------
def process(cid, odir):
    import SimpleITK as sitk
    d = DATA_ROOT / cid / "FDG"
    pet_img = sitk.ReadImage(str(d / "PET.nii.gz"), sitk.sitkFloat64)
    ttb_img = sitk.ReadImage(str(d / "TTB.nii.gz"))
    if (ttb_img.GetSize() != pet_img.GetSize()
            or max(abs(a - b) for a, b in zip(ttb_img.GetOrigin() + ttb_img.GetSpacing() + ttb_img.GetDirection(),
                                               pet_img.GetOrigin() + pet_img.GetSpacing() + pet_img.GetDirection())) > 1e-3):
        raise RuntimeError(f"{cid}: TTB not on the FDG PET grid")
    bw = sitk.GetArrayFromImage(ttb_img) > 0
    pet = sitk.GetArrayFromImage(pet_img)
    D = np.array(pet_img.GetDirection()).reshape(3, 3); S = np.array(pet_img.GetSpacing())
    vox_ml = abs(np.linalg.det(D * S)) / 1000.0
    half = 0.5 * np.abs(D * S).sum(axis=1)                    # world-axis half extent of one voxel
    tot_vox = int(bw.sum()); tot_ml = tot_vox * vox_ml
    comps, pats = [], []
    for c in CONNS:
        lab, n = label_components(bw, c)
        out = sitk.GetImageFromArray(lab); out.CopyInformation(ttb_img)
        fn = odir / f"P1_{cid}_FDG_conn{c}.nii.gz"
        sitk.WriteImage(out, str(fn), True)
        sum_ml = 0.0
        if n:
            idx = ndimage.value_indices(lab, ignore_value=0)
            for q in range(1, n + 1):
                kji = np.column_stack(idx[q]); pv = pet[idx[q]]
                w = phys_lps(pet_img, kji)
                ext = w.max(0) - w.min(0) + 2 * half
                ml = len(kji) * vox_ml; sum_ml += ml
                cl = w.mean(0)
                comps.append(dict(case_id=cid, connectivity=c, component_id=q, voxels=len(kji), volume_mL=ml,
                                  SUVmax=float(pv.max()), SUVmean=float(pv.mean()),
                                  centroid_x_LPS_mm=cl[0], centroid_y_LPS_mm=cl[1], centroid_z_LPS_mm=cl[2],
                                  centroid_x_RAS_mm=-cl[0], centroid_y_RAS_mm=-cl[1], centroid_z_RAS_mm=cl[2],
                                  max_extent_mm=float(ext.max()),
                                  eligible_A=bool(eligible(ml, len(kji), "A")),
                                  eligible_B=bool(eligible(ml, len(kji), "B"))))
        pats.append(dict(case_id=cid, connectivity=c, n_components=n, TTB_voxels=tot_vox, V_FDG_total_mL=tot_ml,
                         sum_component_mL=sum_ml, voxel_mL=vox_ml, label_image=fn.name, label_image_md5=md5(fn)))
    return comps, pats


# ------------------------------ self-test (synthetic) ------------------------------
def selftest():
    res = []
    ok = lambda n_, c_: res.append((n_, bool(c_)))
    bw = np.zeros((6, 6, 6), bool)
    bw[1, 1, 1] = bw[1, 1, 2] = True          # face pair
    bw[3, 3, 3] = bw[4, 4, 3] = True          # edge pair (differs in k and j)
    bw[0, 5, 5] = bw[1, 4, 4] = True          # corner pair
    ok("6-connectivity: 5 components", label_components(bw, 6)[1] == 5)
    ok("18-connectivity: 4 components", label_components(bw, 18)[1] == 4)
    ok("26-connectivity: 3 components", label_components(bw, 26)[1] == 3)
    # ordering: raster (k, j, i) order of first voxel == MATLAB column-major over (i, j, k)
    lab, n = label_components(bw, 26)
    firsts = [np.ravel_multi_index(np.argwhere(lab == q)[0], bw.shape) for q in range(1, n + 1)]
    ok("labels numbered by first voxel in raster order", firsts == sorted(firsts))
    mat_order = sorted(range(1, n + 1), key=lambda q: min(i + 6 * j + 36 * k for k, j, i in np.argwhere(lab == q)))
    ok("raster order (k,j,i) equals column-major order (i,j,k)", mat_order == list(range(1, n + 1)))
    # eligibility boundaries
    ok("rule A at 0.1 mL and 3 voxels -> eligible", eligible(0.1, 3, "A"))
    ok("rule A at 0.0999 mL -> excluded", not eligible(0.0999, 50, "A"))
    ok("rule A with 2 voxels -> excluded", not eligible(5.0, 2, "A"))
    ok("rule B at 0.25 mL -> eligible; 0.2499 -> excluded", eligible(0.25, 3, "B") and not eligible(0.2499, 9, "B"))
    # volume reconciliation
    lab6, n6 = label_components(bw, 6)
    ok("sum of component voxels = TTB voxels (all connectivities)",
       all(sum((label_components(bw, c)[0] == q).sum() for q in range(1, label_components(bw, c)[1] + 1)) == bw.sum()
           for c in CONNS))
    ok("empty mask -> 0 components", label_components(np.zeros((3, 3, 3), bool), 26)[1] == 0)
    for n_, g in res:
        print(f"  [{'PASS' if g else 'FAIL'}] {n_}")
    print(f"selftest: {sum(g for _, g in res)}/{len(res)} passed")
    return all(g for _, g in res)


# ------------------------------ run ------------------------------
def run(preflight):
    import SimpleITK as sitk
    t0 = time.time(); run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = "P1_PREFLIGHT_" if preflight else "P1_"
    md = {n: md5(d / n) if (d / n).exists() else "MISSING" for (d, n) in INPUTS}
    bad = {n: md[n] for (d, n), e in INPUTS.items() if md[n] != e}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(INPUTS)}/{len(INPUTS)} match")
    AC = pd.read_csv(OUT_DIR / "Audit1D_component_level_20261030_215034.csv")
    AC = AC[AC.tracer == "FDG"]
    cases = sorted(AC.case_id.unique())
    if len(cases) != 100:
        sys.exit("STOP: Audit 1D does not list 100 FDG patients")
    if preflight:
        cases = preflight_cases(cases, pd.read_csv(OUT_DIR / "A3_visualQC_KEY_do_not_open_20261105.csv"))
    odir = OUT_DIR / f"{tag}{run_id}"; odir.mkdir()
    comps, pats = [], []
    for i, cid in enumerate(cases, 1):
        c_, p_ = process(cid, odir)
        comps += c_; pats += p_
        print(f"  [{i:3d}/{len(cases)}] {cid}: components 6/18/26 = "
              f"{'/'.join(str(p['n_components']) for p in p_)}  ({time.time() - t0:.0f} s)")
    C = pd.DataFrame(comps); P = pd.DataFrame(pats)
    # ---- reconciliation ----
    rec = P.assign(rel=(P.sum_component_mL - P.V_FDG_total_mL).abs() / P.V_FDG_total_mL.where(P.V_FDG_total_mL > 0, 1))
    rec_ok = bool((rec.rel <= TOL_REL).all() and P.groupby("case_id").V_FDG_total_mL.nunique().eq(1).all())
    print(f"Reconciliation (sum of components = V_FDG,total; identical across connectivities): "
          f"{'PASS' if rec_ok else 'FAIL'}")
    # ---- Audit 1D comparison (per case, connectivity, component ID) ----
    A = AC[AC.case_id.isin(cases)].copy()
    M = C.merge(A, on=["case_id", "connectivity", "component_id"], how="outer", suffixes=("", "_1D"), indicator=True)
    M["count_ok"] = M["_merge"] == "both"
    both = M[M.count_ok]
    M.loc[M.count_ok, "voxels_ok"] = both.voxels == both.voxels_1D
    M.loc[M.count_ok, "volume_ok"] = (both.volume_mL - both.volume_mL_1D).abs() <= TOL_ML
    M.loc[M.count_ok, "SUV_ok"] = ((both.SUVmax - both.SUVmax_1D).abs() <= TOL_SUV) & ((both.SUVmean - both.SUVmean_1D).abs() <= TOL_SUV)
    M.loc[M.count_ok, "centroid_ok"] = (np.sqrt((both.centroid_x_RAS_mm - both.centroid_x_RAS_mm_1D) ** 2
                                                + (both.centroid_y_RAS_mm - both.centroid_y_RAS_mm_1D) ** 2
                                                + (both.centroid_z_RAS_mm - both.centroid_z_RAS_mm_1D) ** 2) <= TOL_MM)
    cols = ["voxels_ok", "volume_ok", "SUV_ok", "centroid_ok"]
    M["PASS"] = M.count_ok & M[cols].fillna(False).all(axis=1)
    cp = odir / f"P1_audit1D_comparison_{run_id}.csv"
    M.drop(columns=["_merge"]).to_csv(cp, index=False)
    n_cnt = C.groupby("connectivity").size().to_dict()
    print(f"Audit 1D comparison: {int(M.PASS.sum())}/{len(M)} components identical (count, voxels, volume, SUV, "
          f"centroid); components 6/18/26 = {n_cnt.get(6, 0)}/{n_cnt.get(18, 0)}/{n_cnt.get(26, 0)}")
    for k in cols:
        nb = int((M[k] == False).sum())
        if nb:
            print(f"   mismatches in {k}: {nb}")
    # ---- Option D (full run only) ----
    od_ok, od_rows = True, []
    if not preflight:
        for (c, r), (n_exp, pct_exp) in OPTION_D.items():
            sub = C[C.connectivity == c]; ex = sub[~sub[f"eligible_{r}"]]
            tot = P[P.connectivity == c].V_FDG_total_mL.sum()
            pct = 100 * ex.volume_mL.sum() / tot
            ok_ = (len(ex) == n_exp) and (round(pct, 2) == pct_exp)
            od_ok &= ok_
            od_rows.append(dict(connectivity=c, rule=r, excluded=len(ex), expected=n_exp, pooled_pct=pct,
                                expected_pct=pct_exp, PASS=ok_))
            print(f"Option D {c}-conn rule {r}: excluded {len(ex)} (expected {n_exp}), pooled {pct:.2f}% "
                  f"(expected {pct_exp:.2f}%) -> {'PASS' if ok_ else 'FAIL'}")
    allok = rec_ok and bool(M.PASS.all()) and od_ok
    if not allok:
        sys.exit(f"STOP: acceptance failed (reconciliation {rec_ok}, Audit 1D {bool(M.PASS.all())}, Option D {od_ok}); "
                 f"see {cp.name}; outputs in {odir.name} are NOT accepted")
    C["V_FDG_total_mL"] = C.case_id.map(P.groupby("case_id").V_FDG_total_mL.first())
    for r in RULES:
        P[f"V_S_{r}_mL"] = [C[(C.case_id == a) & (C.connectivity == b) & ~C[f"eligible_{r}"]].volume_mL.sum()
                            for a, b in zip(P.case_id, P.connectivity)]
    ccp = odir / f"P1_components_{run_id}.csv"; C.to_csv(ccp, index=False)
    pp = odir / f"P1_patients_{run_id}.csv"; P.to_csv(pp, index=False)
    extra = []
    if od_rows:
        op = odir / f"P1_optionD_check_{run_id}.csv"; pd.DataFrame(od_rows).to_csv(op, index=False); extra = [op]
    man = dict(stage="P1 native FDG components and eligibility", mode="PREFLIGHT" if preflight else "FULL",
               run=run_id, script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs=md,
               n_patients=len(cases), components_by_connectivity=n_cnt,
               outputs={p.name: md5(p) for p in [ccp, pp, cp] + extra},
               label_images={r.label_image: r.label_image_md5 for r in P.itertuples()},
               params=dict(CONNS=CONNS, RULES=RULES, TOL_ML=TOL_ML, TOL_SUV=TOL_SUV, TOL_MM=TOL_MM, TOL_REL=TOL_REL),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=__import__("scipy").__version__, SimpleITK=sitk.Version_VersionString()),
               runtime_s=round(time.time() - t0), status="COMPLETE")
    mp = odir / f"P1_manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"done: {odir.name}; manifest MD5 {md5(mp)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.preflight)
