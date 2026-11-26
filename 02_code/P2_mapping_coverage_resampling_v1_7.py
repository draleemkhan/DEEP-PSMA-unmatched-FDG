"""
P2_mapping_coverage_resampling.py - PAPER 3, stage P2: mapping of eligible native FDG components into PSMA space,
registration-QC flags, coverage flags and per-component resampling onto the native PSMA PET grid.
NO matching and NO DMTVF: PSMA PET and PSMA TTB are read for their GRID (header) only; no PSMA voxel value is read.

Implements: consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) sections 3-4 and Section R / R.1
(failed regions, R6 assignment, band, precedence inputs); implementation spec LOCKED (MD5
a5a2ec19d284deb44bc9ceebf8636658) I0, I3 (with the author's clarifications: per-component voxel sets, overlaps kept;
band alternative mapping with its own coverage and zero-survivor checks).

Per patient (connectivity --conn, default 26; components eligible under rule A; rule B is a subset):
  transforms  T_A(p) = FDG/rigid.tfm(p) - t (t = 0 for train_0018/0049/0067); region correction dT(q) = R q + c where
              the region's final transform is Amendment 3a (A3_final_acceptance_per_region.csv); resampling transforms
              RA / RB exactly as the locked scripts (A3_visualQC_views_v1_2.resampling_transforms), composition checked
  R6          region of the T_A-mapped centroid (planes from SAP_R1_plane_table; tie -> trunk); final transform =
              that region's transform if adopted, else T_A; relevant transforms and band from SAP_R1_region_rules
  flags       VI_registration: any native voxel centre, mapped with any relevant transform, in a failed region
              (head/neck/pelvis per SAP R.4 + R.1.4); coverage_fov: any mapped voxel centre (any relevant transform)
              outside the PSMA PET grid or closer than tau_i to its outer face; per-transform flags kept
  band        within 15 mm of an R6 plane (C2/C3, L5/S1) or crossing it AND the two sides use different transforms:
              alternative transform = the other side's; its own coverage and zero-survivor flags
  resampling  the native component label image is resampled (nearest neighbour) once per transform used; each
              component's voxel set is taken from the image of ITS transform (final; alt for band) and stored
              separately (overlaps between components with different transforms are kept); % volume change;
              zero-survivor flag; collision count
Outputs (06_output\\P2_[PREFLIGHT_]conn<c>_<run>\\): P2_components_<run>.csv, P2_voxelsets_<case>.npz
  (component_id, role, PSMA-PET linear indices in (k, j, i) C order), QC images P2_<case>_final_labels.nii.gz and
  P2_<case>_overlap_count.nii.gz, P2_grids_<run>.csv (fingerprints), manifest.
Acceptance (stop on failure): C1 MD5 of all inputs; transform composition <= 1e-4 mm; PSMA TTB header on the PSMA PET
  grid; component set (IDs, voxels, volumes) identical to P1; image-level check for components with >= 50 resampled
  voxels and no coverage flag: |centroid of resampled voxels - mapped centroid| <= max PSMA PET voxel dimension.

Revision (11 Nov 2026, 23:55): image-level acceptance filter fixed (explicit Boolean conversion of the coverage
flag; the first full run P2_conn26_20261111_235057 stopped because coverage-flagged components were not excluded;
that run is kept, not accepted). Three self-test cases added. No change to mapping, flags or resampling.

v1.2 (12 Nov 2026, 01:00; file P2_mapping_coverage_resampling_v1_2.py): registration sensitivity scenarios of SAP
section R.6 / R.1.6 and the unblinding note section 5 (--scenario; default none = identical to v1.1):
  S1  the 9 grade-1 regions use their Amendment 3a candidate and do not fail QC (heads train_0007, 0012, 0023, 0028,
      0066, 0068, 0089; pelvises train_0069, 0089)
  S2  the 6 NA heads do not fail QC (train_0010, 0033, 0042, 0053, 0079, 0095); transforms unchanged
  S4  failed necks released (neck removed from every failed set); transforms unchanged
Only the failed-region set (and, for S1, the final transform of the 9 regions) changes; everything else as v1.1.
Output folder P2_<scenario>_conn<c>_<run>; scenario recorded in the manifest.

v1.3 (12 Nov 2026, 01:20; file P2_mapping_coverage_resampling_v1_3.py): S3b and S3c (SAP section R.6, R.1.6), per the
dated implementation clarification 00_protocol\\StepSAP_S3_implementation_clarification_20261112.md
(MD5 513196d78278107851fcd4e70cdea376). Failed regions = the locked set (29 regions). Points = native FDG voxel centres
of each eligible component, mapped with T_A ONLY; a point is "failed" if its region (region_of, locked planes and tie
rules) is a failed region of the patient.
  S3b  VI_registration = (fraction of failed points > 0.5); everything else as v1.1
  S3c  all points failed -> VI_registration (whole); none failed -> not V_I by this rule; otherwise SPLIT: the failed
       points are removed; R6 final transform, relevant transforms, band and alternative transform are those of the
       WHOLE component (not recomputed); coverage (FOV / tau margin, every transform used), resampling, zero-survivor
       and the image-level check use the REMAINDER only. Resampling uses a native label image in which only the failed
       points of split components are set to 0 (nearest neighbour: other components' voxel sets are unchanged).
  New columns (every scenario): s3_n_failed_TA, s3_frac_failed_TA, s3_split. Volumes are NOT computed here: P3 splits
  the native component volume as n_failed x (native volume / native voxel count) (clarification rule 5).

v1.4 (12 Nov 2026, 08:40; file P2_mapping_coverage_resampling_v1_4.py): registration perturbation mode (SAP section 8.1)
per the dated implementation clarification 00_protocol\\StepSAP_perturbation_implementation_clarification_20261112.md
(MD5 5f280df40970eec504bb4f08d421d69d). Perturbations read from the MATLAB generator output (PERT_generate_v1_0.m,
MD5 fca5d75c04ec4b75bc5b7b23565137f3; CSV MD5 verified through its manifest). Per patient i and perturbation k:
c_i = mean of all eligible native voxel centres mapped with T_A (unperturbed); every mapped point of every transform
used is replaced by P_k(q) = R_k (q - c_i) + c_i + t_k; resampling uses x -> RT(P_k^-1(x)) (composite transform,
checked per patient, perturbation and transform: native points recovered within TOL_COMPOSE). Recomputed: coverage,
resampling, zero survivors, image-level check; kept from the unperturbed mapping: R6 final transform, relevant
transforms, band and alternative transform, registration-QC flag. QC images are not written in this mode.
  --pertdir PERT_perturbations_<run> --level L1|L2 [--preflight]   one parent folder per level, subfolders k01..k50
  --pertzero [--preflight]                                          gate: R = I, t = 0 must reproduce the primary P2 run

v1.5 (12 Nov 2026, 08:35; file P2_mapping_coverage_resampling_v1_5.py): zero-perturbation gate comparison fixed. In v1.4
the non-numeric columns were compared after astype(str); with the pandas version on the analysis PC, missing values
stay missing after that conversion and missing != missing, so identical empty cells (the 871 V_S candidates) were
reported as differences (run P2_PERTZERO_conn26_20261112_082256 stopped; kept, not accepted; an independent re-read
found 0 differences in all flagged columns). Missing values are now compared explicitly as equal (frames_identical,
self-tested). No change to mapping, perturbation, flags or resampling.

v1.6 (12 Nov 2026, 09:05; file P2_mapping_coverage_resampling_v1_6.py): image-level acceptance in PERTURBATION runs per
the dated note 00_protocol\\StepSAP_perturbation_imagecheck_note_20261112.md (MD5 ced70ff6ce6597d5e1e8cc5a730fc908;
author decision (a)). Both v1.5 level runs stopped on one small component each (offset 4.45 / 3.41 mm); a diagnostic
showed exact agreement with an independent nearest-neighbour derivation (sampling effect, not an error). Now, in
perturbation runs only, a component above the centroid tolerance is re-derived exactly (nn_rederive: each PSMA voxel
centre in a bounding box mapped back with the composed resampling transform, nearest native voxel by rounding half up,
native label read); the run stops only if the re-derived set differs. Such components are listed in
P2_imagecheck_exceptions_kNN.csv. Unperturbed runs unchanged. Composition, reconciliation and all other checks unchanged.

v1.7 (12 Nov 2026, 12:30; file P2_mapping_coverage_resampling_v1_7.py): supplied-registration comparator (POST-HOC) per
the dated note 00_protocol\\StepSAP_supplied_registration_comparator_note_20261112.md (MD5 ead22c1221bd993b2ce667d4d82ccadc).
  --supplied [--preflight]   one subfolder k01 (level SUPPLIED). Every mapped point of every transform used by a component
  (final and alternative) is replaced by the supplied rigid mapping S(p) = A p + b; every resampling transform is replaced by
  the supplied composite (supplied rigid only; as for the fallback patients). Kept from the primary mapping (frozen at
  component level): eligibility, R6 final-transform label, relevant transforms, band candidacy and alternative label,
  registration-QC flag. Recomputed: coverage, resampling, zero survivors, image-level check (exact re-derivation as v1.6).
  Final and alternative voxel sets are identical by construction (band V_I = 0 in P3). Gate: components of the three
  fallback patients whose primary final transform is TA and alternative is absent or TA must have voxel sets identical to
  the primary P2 run (subset size printed). All other modes are unchanged (zero-perturbation gate must reproduce the
  primary P2 run with this version).

Usage (folder 02_code; Amendment3_regional_rigid.py, A3_visualQC_views_v1_2.py, SAP_R1_region_rules.py alongside):
  & "<totalseg python>" -u P2_mapping_coverage_resampling.py --selftest
  & "<totalseg python>" -u P2_mapping_coverage_resampling.py --preflight
  & "<totalseg python>" -u P2_mapping_coverage_resampling.py            [--conn 18|6 for sensitivity partitions]
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
from scipy import ndimage

HERE = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("vqc", HERE / "A3_visualQC_views_v1_2.py")
v = importlib.util.module_from_spec(_s); _s.loader.exec_module(v)
a3 = v.a3
_s = importlib.util.spec_from_file_location("rr", HERE / "SAP_R1_region_rules.py")
rr = importlib.util.module_from_spec(_s); _s.loader.exec_module(rr)

DATA_ROOT, OUT_DIR = a3.DATA_ROOT, a3.OUT_DIR
PROTO_DIR = OUT_DIR.parent / "00_protocol"
P1_RUN = "P1_20261111_233936"
INPUTS = {
    (HERE, "Amendment3_regional_rigid.py"): "b21a405145526d56a59bb5950ef53acf",
    (HERE, "A3_visualQC_views_v1_2.py"): "e8770c32865cf828ed27406a7fef6d9c",
    (HERE, "SAP_R1_region_rules.py"): "beedabcfe956f3c9640a061b3855e867",
    (OUT_DIR, "Audit1C_e_per_patient_20261031_141128.csv"): "c66da788a6c73737163b4c09ed76f494",
    (OUT_DIR, "Amendment3_transforms_20261102_073040.csv"): "6dd27f366127d1941981489d57f0a99a",
    (OUT_DIR, "A3_visualQC_unblinded_20261108/A3_final_acceptance_per_region.csv"): "225e4e673b31783e5bf1fa210b9bf42c",
    (OUT_DIR, "SAP_R1_plane_table_20261111_214952.csv"): "db6fe2142e5c6fe7fe8ee6d8f9bc0f54",
    (OUT_DIR, f"{P1_RUN}/P1_components_20261111_233936.csv"): "ffafbc7917e1b595eadacafab247fe88",
    (OUT_DIR, f"{P1_RUN}/P1_patients_20261111_233936.csv"): "f9b5b5908e523d1122547ef07d5c0afd",
    (OUT_DIR, "A3_visualQC_KEY_do_not_open_20261105.csv"): "c9a4548d3daf13d3e3dadc96171a9a94",   # preflight IDs only
    (PROTO_DIR, "StepSAP_consolidated_LOCKED_20261111.md"): "0c64711db1207e31a70cf3f558ccbd7d",
    (PROTO_DIR, "PET_correspondence_implementation_spec_LOCKED_20261111.md"): "a5a2ec19d284deb44bc9ceebf8636658",
}
PRIMARY_P2 = "P2_conn26_20261111_235603"
PERT_GEN_MD5 = "fca5d75c04ec4b75bc5b7b23565137f3"
PERT_NOTE = (PROTO_DIR, "StepSAP_perturbation_implementation_clarification_20261112.md", "5f280df40970eec504bb4f08d421d69d")
SUP_NOTE = (PROTO_DIR, "StepSAP_supplied_registration_comparator_note_20261112.md", "ead22c1221bd993b2ce667d4d82ccadc")
FALLBACK_IDS = ("train_0018", "train_0049", "train_0067")
BAND_MM = rr.BAND_MM
TOL_COMPOSE = 1e-4
MIN_VOX_IMAGE_CHECK = 50


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ pure helpers (tested) ------------------------------
def fingerprint(size, spacing, origin, direction):
    """grid fingerprint: size; spacing and origin rounded to 1e-4 mm; direction cosines to 1e-6 (SAP item 5)."""
    s = "|".join([",".join(str(int(x)) for x in size),
                  ",".join(f"{round(float(x), 4):.4f}" for x in spacing),
                  ",".join(f"{round(float(x), 4):.4f}" for x in origin),
                  ",".join(f"{round(float(x), 6):.6f}" for x in direction)])
    return hashlib.md5(s.encode()).hexdigest()


def cont_index(pts, origin, spacing, direction):
    D = np.asarray(direction, float).reshape(3, 3)
    return np.linalg.solve(D * np.asarray(spacing, float), (np.atleast_2d(pts) - np.asarray(origin, float)).T).T


def coverage_fail(pts, origin, spacing, direction, size, tau):
    """True if any point lies outside the grid (outer voxel faces) or closer than tau (mm) to an outer face.
    Non-oblique grid assumed (verified: obliquity 0 in all PSMA PET grids, Audit 1D / tau4)."""
    if len(pts) == 0:
        return False
    ci = cont_index(pts, origin, spacing, direction)
    sp = np.asarray(spacing, float); sz = np.asarray(size, float)
    lo = (ci + 0.5) * sp                       # mm to the lower outer face along each axis
    hi = (sz - 0.5 - ci) * sp                  # mm to the upper outer face
    return bool(np.any(np.minimum(lo, hi) < tau))


def band_alternative(zc, zmin, zmax, P, final, adopted):
    """R6 band: within BAND_MM of an R6 plane or crossing it; alternative = the other side's transform if it differs.
    final in {'TA','head','pelvis'}; returns (band, alt or None)."""
    for plane, side in (("z_head", "head"), ("z_pelvis", "pelvis")):
        zp = P[plane]
        if not np.isfinite(zp):
            continue
        if abs(zc - zp) <= BAND_MM or (zmin < zp < zmax):
            other = "TA" if final == side else (side if adopted[side] else "TA")
            return True, (other if other != final else None)
    return False, None


def final_transform(zc, P, adopted):
    reg = rr.r6_region(zc, P)[0]
    return reg if reg in ("head", "pelvis") and adopted[reg] else "TA", reg


def as_bool(series):
    """explicit Boolean conversion of a flag column (object columns of Python bools/NaN would make ~ a bitwise int
    operation; fix after the stopped run P2_conn26_20261111_235057)."""
    return series.map(lambda x: x is True or (isinstance(x, (bool, np.bool_)) and bool(x))
                      or (isinstance(x, str) and x.strip().lower() == "true")).astype(bool)


def image_check_failures(E, maxsp):
    """components checked: >= MIN_VOX_IMAGE_CHECK resampled voxels (final) and NO coverage flag on the final mapping;
    failure: |resampled centroid - mapped centroid| > max PSMA PET voxel dimension. Returns (checked, failures)."""
    if len(E) == 0 or "voxels_resampled_final" not in E:
        return E.iloc[0:0], E.iloc[0:0]
    chk = E[(E["voxels_resampled_final"].fillna(0) >= MIN_VOX_IMAGE_CHECK) & ~as_bool(E["coverage_fov_final"])]
    return chk, chk[chk["image_check_mm"] > chk["case_id"].map(maxsp)]


S1_RECOVER = {("train_0007", "head"), ("train_0012", "head"), ("train_0023", "head"), ("train_0028", "head"),
              ("train_0066", "head"), ("train_0068", "head"), ("train_0089", "head"),
              ("train_0069", "pelvis"), ("train_0089", "pelvis")}
S2_NA_HEADS = {"train_0010", "train_0033", "train_0042", "train_0053", "train_0079", "train_0095"}


def scenario_failed(cid, scenario):
    """failed regions of a patient under a sensitivity scenario (SAP R.6, R.1.6, unblinding note 5)."""
    f = set(rr.failed_regions(cid))
    if scenario == "S1":
        f -= {reg for (c, reg) in S1_RECOVER if c == cid}
    elif scenario == "S2" and cid in S2_NA_HEADS:
        f.discard("head")
    elif scenario == "S4":
        f.discard("neck")
    return f


def scenario_adopted(cid, adopted, scenario):
    if scenario == "S1":
        return {reg: (adopted[reg] or (cid, reg) in S1_RECOVER) for reg in adopted}
    return dict(adopted)


def s3_failed_points(ta_z, P, failed):
    """S3 (clarification rules 1-2): Boolean per point, True if its T_A-mapped z lies in a failed region."""
    if not failed:
        return np.zeros(len(ta_z), bool)
    return np.isin(rr.region_of(np.asarray(ta_z, float), P), sorted(failed))


def s3_decision(inF, scenario):
    """returns (VI_registration of the whole component, keep mask for the assessed remainder or None, split flag).
    S3b: > 50 % failed -> V_I (exactly 50 % is not); S3c: all failed -> V_I; none -> unchanged; else split."""
    n, nf = len(inF), int(np.sum(inF))
    if scenario == "S3b":
        return nf > 0.5 * n, None, False
    if scenario == "S3c":
        if nf == n:
            return True, None, False
        if nf == 0:
            return False, None, False
        return False, ~np.asarray(inF, bool), True
    raise ValueError(scenario)


def perturb_points(q, R, t, c):
    """P(q) = R (q - c) + c + t for rows of q (PSMA physical LPS, mm)."""
    return (np.atleast_2d(q) - c) @ np.asarray(R, float).T + c + t


def supplied_points(pts, A, b):
    """v1.7: supplied rigid mapping S(p) = A p + b (no translation correction, no regional refinement)."""
    return np.atleast_2d(pts) @ np.asarray(A, float).T + np.asarray(b, float)


def fallback_subset(R0):
    """v1.7: primary P2 rows of the fallback validation subset (eligible, final TA, alternative absent or TA)."""
    alt = R0.alt_transform.fillna("").astype(str) if "alt_transform" in R0 else pd.Series("", index=R0.index)
    return R0[R0.case_id.isin(FALLBACK_IDS) & (R0.status_P2 == "eligible") & (R0.final_transform == "TA")
              & alt.isin(["", "TA"])]


def inverse_perturbation(R, t, c):
    """SimpleITK affine x -> R^T (x - c - t) + c (the inverse of perturb_points)."""
    import SimpleITK as sitk
    R = np.asarray(R, float); c = np.asarray(c, float); t = np.asarray(t, float)
    A = sitk.AffineTransform(3)
    A.SetMatrix(tuple(R.T.flatten())); A.SetTranslation(tuple(c - R.T @ (c + t)))
    return A


def compose(outer, inner):
    """x -> outer(inner(x)) (ITK composite: the transform added last is applied first)."""
    import SimpleITK as sitk
    C = sitk.CompositeTransform(3); C.AddTransform(outer); C.AddTransform(inner)
    return C


def nn_rederive(Tr, lab_img, lab, comp, mapped, g, margin=4):
    """independent nearest-neighbour derivation of a component's resampled voxel set (C-order linear indices on grid g):
    PSMA voxel centres in the bounding box of the mapped points (+ margin voxels) mapped with Tr (resampling direction),
    nearest native voxel by floor(x + 0.5), native label == comp."""
    Dp = np.array(g["direction"], float).reshape(3, 3); Sp = np.array(g["spacing"], float); Op = np.array(g["origin"], float)
    size = np.array(g["size"], int); shape = tuple(int(x) for x in size[::-1]); nsz = np.array(lab_img.GetSize(), int)
    ci = np.linalg.solve(Dp * Sp, (np.atleast_2d(mapped) - Op).T).T
    lo = np.maximum(np.floor(ci.min(0)).astype(int) - margin, 0); hi = np.minimum(np.ceil(ci.max(0)).astype(int) + margin, size - 1)
    out = set()
    for i in range(lo[0], hi[0] + 1):
        for j in range(lo[1], hi[1] + 1):
            for k in range(lo[2], hi[2] + 1):
                x = Op + Dp @ (np.array([i, j, k]) * Sp)
                y = Tr.TransformPoint(tuple(float(u) for u in x))
                n_ = np.floor(np.array(lab_img.TransformPhysicalPointToContinuousIndex(y)) + 0.5).astype(int)
                if np.all(n_ >= 0) and np.all(n_ < nsz) and lab[n_[2], n_[1], n_[0]] == comp:
                    out.add(int(np.ravel_multi_index((k, j, i), shape)))
    return out


def frames_identical(A0, A1, atol=1e-9):
    """columns of A0 that differ in A1: numeric with tolerance (NaN = NaN); others as text with missing = missing."""
    diffs = []
    if len(A0) != len(A1):
        return ["row count"]
    for col in A0.columns:
        if col not in A1.columns:
            diffs.append(f"missing {col}"); continue
        x, y = A0[col].reset_index(drop=True), A1[col].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            if not np.allclose(x.astype(float), y.astype(float), atol=atol, rtol=0, equal_nan=True):
                diffs.append(col)
        else:
            xs = [None if pd.isna(v) else str(v) for v in x]
            ys = [None if pd.isna(v) else str(v) for v in y]
            if xs != ys:
                diffs.append(col)
    return diffs


def load_perturbations(pertdir):
    pd_ = OUT_DIR / pertdir; rid = pertdir[-15:]
    m = json.loads((pd_ / f"PERT_manifest_{rid}.json").read_text(encoding="utf-8"))
    csv = pd_ / m["csv"]
    if m.get("script_md5") != PERT_GEN_MD5 or m.get("status") != "COMPLETE" or md5(csv) != m.get("csv_md5"):
        sys.exit("STOP C1: perturbation file not produced by the accepted generator or changed")
    T = pd.read_csv(csv)
    return T, {f"PERT_manifest_{rid}.json": md5(pd_ / f"PERT_manifest_{rid}.json"), m["csv"]: md5(csv)}


def preflight_cases(cases, oldkey):
    return [c for c in cases if c not in set(oldkey["case_id"])][:2]


# ------------------------------ per patient ------------------------------
def read_grid(path):
    import SimpleITK as sitk
    r = sitk.ImageFileReader(); r.SetFileName(str(path)); r.ReadImageInformation()
    return dict(size=r.GetSize(), spacing=r.GetSpacing(), origin=r.GetOrigin(), direction=r.GetDirection())


def grid_image(g):
    import SimpleITK as sitk
    im = sitk.Image([int(x) for x in g["size"]], sitk.sitkUInt8)
    im.SetSpacing(g["spacing"]); im.SetOrigin(g["origin"]); im.SetDirection(g["direction"])
    return im


def process(cid, conn, C1, PP, TF, FA, PL, odir, scenario="none", pert=None):
    import SimpleITK as sitk
    X = v.patient_transforms(cid, PP, TF); fb = cid in a3.FALLBACKS
    A, b, t = X["A"], X["b"], X["t"]
    adopted = scenario_adopted(cid, {reg: (FA.get((cid, reg)) == "3a") for reg in ("head", "pelvis")}, scenario)
    RT = {"TA": None}; comp = 0.0
    for reg in ("head", "pelvis"):
        R, c = X[reg]
        RA, RB, w = v.resampling_transforms(X["dP"], t, R, c, A, b, fb)
        comp = max(comp, w); RT["TA"] = RA
        if adopted[reg]:
            RT[reg] = RB
    if comp > TOL_COMPOSE:
        raise RuntimeError(f"{cid}: composition check failed ({comp:.2e} mm)")

    def mapf(name, pts):
        q = (pts @ A.T + b) - t
        if name == "TA":
            return q
        R, c = X[name]
        return q @ R.T + c
    g_pet = read_grid(X["dP"] / "PET.nii.gz"); g_ttb = read_grid(X["dP"] / "TTB.nii.gz")
    if g_pet["size"] != g_ttb["size"] or max(abs(a_ - b_) for a_, b_ in zip(
            g_pet["spacing"] + g_pet["origin"] + g_pet["direction"], g_ttb["spacing"] + g_ttb["origin"] + g_ttb["direction"])) > 1e-3:
        raise RuntimeError(f"{cid}: PSMA TTB not on the PSMA PET grid")
    ref = grid_image(g_pet)
    tau = float(np.sqrt(np.sum(np.square(g_pet["spacing"]))))
    P = dict(z_head=PL.loc[cid, "z_head_C2C3"], z_c7t1=PL.loc[cid, "z_C7T1"], z_pelvis=PL.loc[cid, "z_pelvis_L5S1"])
    failed = scenario_failed(cid, scenario)
    lab_img = sitk.ReadImage(str(OUT_DIR / P1_RUN / f"P1_{cid}_FDG_conn{conn}.nii.gz"))
    lab = sitk.GetArrayFromImage(lab_img)
    idx = ndimage.value_indices(lab, ignore_value=0)
    D = np.array(lab_img.GetDirection()).reshape(3, 3); S = np.array(lab_img.GetSpacing()); O = np.array(lab_img.GetOrigin())
    comps = C1[(C1.case_id == cid) & (C1.connectivity == conn)]
    if len(comps) != len(idx) or any(len(idx[q][0]) != n for q, n in zip(comps.component_id, comps.voxels)):
        raise RuntimeError(f"{cid}: native label image does not match the P1 component table")
    rows, need = [], set()
    lab_rs = lab.copy() if scenario == "S3c" else lab      # S3c: failed points of split components removed
    c_pert = None; RTp = {}; comp_pert = 0.0; mapped_final = {}
    sup = pert is not None and pert.get("supplied", False)
    if sup:                                                # v1.7: supplied composite for every transform used
        import SimpleITK as sitk
        RS = sitk.CompositeTransform(3); RS.AddTransform(sitk.ReadTransform(str(X["dP"] / "rigid.tfm")))
        RTp = {name: RS for name in list(RT) + ["head", "pelvis"]}
        c_pert = np.zeros(3)
    elif pert is not None:                                 # v1.4: patient rotation centre from T_A (unperturbed)
        acc_, n_ = np.zeros(3), 0
        for r in comps.itertuples():
            if r.eligible_A:
                kji_ = np.column_stack(idx[r.component_id])
                ta_ = mapf("TA", O + (kji_[:, ::-1].astype(float) * S) @ D.T)
                acc_ += ta_.sum(0); n_ += len(ta_)
        c_pert = acc_ / n_ if n_ else np.zeros(3)
        Pinv = inverse_perturbation(pert["R"], pert["t"], c_pert)
        RTp = {name: compose(T_, Pinv) for name, T_ in RT.items()}
    for r in comps.itertuples():
        row = dict(case_id=cid, connectivity=conn, component_id=r.component_id, voxels_native=r.voxels,
                   volume_native_mL=r.volume_mL, eligible_A=r.eligible_A, eligible_B=r.eligible_B)
        if not r.eligible_A:
            row.update(status_P2="V_S_candidate"); rows.append(row); continue
        kji = np.column_stack(idx[r.component_id])
        pts = O + (kji[:, ::-1].astype(float) * S) @ D.T
        ta = mapf("TA", pts); zc = float(ta[:, 2].mean())
        fin, reg = final_transform(zc, P, adopted)
        rel = rr.relevant_transforms(ta, P, adopted)
        band, alt = band_alternative(zc, ta[:, 2].min(), ta[:, 2].max(), P, fin, adopted)
        mapped = {n: (ta if n == "TA" else mapf(n, pts)) for n in set(rel) | {fin} | ({alt} if alt else set())}
        inF = s3_failed_points(ta[:, 2], P, failed)
        split = False
        if scenario in ("S3b", "S3c"):
            vi_reg, keep, split = s3_decision(inF, scenario)
            if split:
                mapped = {n: m[keep] for n, m in mapped.items()}
                fk = kji[~keep]
                lab_rs[fk[:, 0], fk[:, 1], fk[:, 2]] = 0
        else:
            vi_reg = rr.registration_vi({n: mapped[n] for n in rel}, P, failed) if failed else False
        if pert is not None:                               # v1.4: perturb every transform used; check composition
            if sup:                                        # v1.7: supplied mapping replaces every transform used
                mapped = {n: supplied_points(pts[keep] if split else pts, A, b) for n in mapped}
            else:
                mapped = {n: perturb_points(m, pert["R"], pert["t"], c_pert) for n, m in mapped.items()}
            sel = np.linspace(0, len(pts) - 1, min(len(pts), 20)).astype(int)
            for n in mapped:
                if n in RTp:
                    back = np.array([RTp[n].TransformPoint(tuple(float(v) for v in x)) for x in mapped[n][sel]])
                    comp_pert = max(comp_pert, float(np.abs(back - pts[sel]).max()))
        cov = {n: coverage_fail(mapped[n], g_pet["origin"], g_pet["spacing"], g_pet["direction"], g_pet["size"], tau)
               for n in mapped}
        need |= {fin} | ({alt} if alt else set())
        mc = mapped[fin].mean(0)
        if pert is not None:
            mapped_final[r.component_id] = mapped[fin]
        row.update(status_P2="eligible", TA_centroid_x=ta[:, 0].mean(), TA_centroid_y=ta[:, 1].mean(), TA_centroid_z=zc,
                   r6_region=reg, final_transform=fin, relevant_transforms="+".join(rel), band=band,
                   alt_transform=alt or "", VI_registration=bool(vi_reg),
                   coverage_fov_any=bool(any(cov[n] for n in set(rel) | {fin})),
                   coverage_fov_final=bool(cov[fin]), coverage_fov_alt=(bool(cov[alt]) if alt else np.nan),
                   mapped_centroid_final_x=mc[0], mapped_centroid_final_y=mc[1], mapped_centroid_final_z=mc[2], tau_mm=tau,
                   s3_n_failed_TA=int(inF.sum()), s3_frac_failed_TA=float(inF.mean()), s3_split=bool(split))
        rows.append(row)
    # ---- resampling: one nearest-neighbour image per transform used; voxel sets per component ----
    vs = {}
    sets_id, sets_role, sets_idx = [], [], []
    if scenario == "S3c":
        src_img = sitk.GetImageFromArray(lab_rs); src_img.CopyInformation(lab_img)
    else:
        src_img = lab_img
    if comp_pert > TOL_COMPOSE:
        raise RuntimeError(f"{cid}: perturbed composition check failed ({comp_pert:.2e} mm)")
    for name in sorted(need):
        T_ = RTp[name] if pert is not None else RT[name]
        res = sitk.GetArrayFromImage(sitk.Resample(src_img, ref, T_, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
        vs[name] = ndimage.value_indices(res, ignore_value=0)
    shape = tuple(int(x) for x in g_pet["size"][::-1])
    owner = np.zeros(shape, np.int32); count = np.zeros(shape, np.uint8)
    lin = lambda tup: np.ravel_multi_index(tup, shape) if len(tup[0]) else np.zeros(0, np.int64)
    for row in rows:
        if row["status_P2"] != "eligible":
            continue
        q = row["component_id"]
        for role, name in (("final", row["final_transform"]), ("alt", row["alt_transform"])):
            if not name:
                continue
            tup = vs[name].get(q, (np.zeros(0, int),) * 3)
            li = lin(tup)
            n_res = len(li)
            row[f"voxels_resampled_{role}"] = n_res
            row[f"zero_survivor_{role}"] = n_res == 0
            sets_id.append(np.full(n_res, q, np.int32)); sets_role.append(np.full(n_res, 0 if role == "final" else 1, np.int8))
            sets_idx.append(li.astype(np.int64))
            if role == "final":
                vnat = row["volume_native_mL"] * (1.0 - row["s3_frac_failed_TA"]) if row["s3_split"] else row["volume_native_mL"]
                row["pct_volume_change_final"] = 100.0 * (n_res * np.prod(g_pet["spacing"]) / 1000.0 - vnat) / vnat
                if n_res:
                    flat_o = owner.reshape(-1); flat_c = count.reshape(-1)
                    flat_c[li] = np.minimum(flat_c[li].astype(int) + 1, 255).astype(np.uint8)
                    flat_o[li] = q
                    kji = np.column_stack(np.unravel_index(li, shape))
                    rc = np.array(g_pet["origin"]) + (kji[:, ::-1] * np.array(g_pet["spacing"])) @ np.array(g_pet["direction"]).reshape(3, 3).T
                    row["image_check_mm"] = float(np.linalg.norm(rc.mean(0) - np.array([row["mapped_centroid_final_x"],
                                                  row["mapped_centroid_final_y"], row["mapped_centroid_final_z"]])))
                    if (pert is not None and n_res >= MIN_VOX_IMAGE_CHECK and not row["coverage_fov_final"]
                            and row["image_check_mm"] > max(g_pet["spacing"])):     # v1.6: exact re-derivation
                        s_ind = nn_rederive(RTp[name], lab_img, lab, q, mapped_final[q], g_pet)
                        row["nn_rederived"] = True; row["nn_rederive_equal"] = bool(s_ind == set(int(x) for x in li))
    # collisions (final voxel sets only)
    flat_c = count.reshape(-1)
    for row in rows:
        if row["status_P2"] == "eligible" and row.get("voxels_resampled_final", 0):
            li = lin(vs[row["final_transform"]][row["component_id"]])
            row["collision_voxels"] = int((flat_c[li] > 1).sum())
    np.savez_compressed(odir / f"P2_voxelsets_{cid}.npz",
                        component_id=np.concatenate(sets_id) if sets_id else np.zeros(0, np.int32),
                        role=np.concatenate(sets_role) if sets_role else np.zeros(0, np.int8),
                        linear_index=np.concatenate(sets_idx) if sets_idx else np.zeros(0, np.int64),
                        grid_shape_kji=np.array(shape))
    for arr, nm in ((owner, "final_labels"), (count, "overlap_count")):
        if pert is not None:
            break                                          # v1.4: QC images not written in perturbation mode
        im = sitk.GetImageFromArray(arr); im.CopyInformation(ref)
        sitk.WriteImage(im, str(odir / f"P2_{cid}_{nm}.nii.gz"), True)
    g_fdg = dict(size=lab_img.GetSize(), spacing=lab_img.GetSpacing(), origin=lab_img.GetOrigin(), direction=lab_img.GetDirection())
    grid = dict(case_id=cid, psma_pet_fingerprint=fingerprint(**g_pet), psma_ttb_fingerprint=fingerprint(**g_ttb),
                fdg_native_fingerprint=fingerprint(**g_fdg), tau_mm=tau, max_spacing_mm=float(max(g_pet["spacing"])),
                composition_max_mm=comp,
                head_final=("3a" if adopted["head"] else "A2"), pelvis_final=("3a" if adopted["pelvis"] else "A2"),
                failed_regions="+".join(sorted(failed)))
    if pert is not None:
        grid.update(pert_level=pert["level"], pert_k=pert["k"], pert_centre_x=c_pert[0], pert_centre_y=c_pert[1],
                    pert_centre_z=c_pert[2], pert_composition_max_mm=comp_pert)
    return rows, grid


# ------------------------------ self-test ------------------------------
def selftest():
    import SimpleITK as sitk
    res = []
    ok = lambda n_, c_: res.append((n_, bool(c_)))
    print("region rules (SAP R.1 module):"); ok("SAP_R1_region_rules self-test", rr.selftest())
    # fingerprint
    g = dict(size=(10, 10, 10), spacing=(4.0, 4.0, 3.0), origin=(1.0, 2.0, 3.0), direction=(1, 0, 0, 0, 1, 0, 0, 0, 1))
    f0 = fingerprint(**g)
    ok("fingerprint insensitive to 1e-6 mm", fingerprint(**dict(g, origin=(1.000001, 2.0, 3.0))) == f0)
    ok("fingerprint sensitive to 1e-3 mm", fingerprint(**dict(g, origin=(1.001, 2.0, 3.0))) != f0)
    # coverage: grid 10^3, 4 mm, origin 0 -> outer faces at -2 and 38 mm; tau 6.93
    G = dict(origin=(0, 0, 0), spacing=(4.0, 4.0, 4.0), direction=(1, 0, 0, 0, 1, 0, 0, 0, 1), size=(10, 10, 10))
    tau = float(np.sqrt(48))
    ok("centre point -> no coverage flag", not coverage_fail(np.array([[18.0, 18, 18]]), tau=tau, **G))
    ok("outside grid -> flag", coverage_fail(np.array([[18.0, 18, 40]]), tau=tau, **G))
    ok("inside but within tau of a face -> flag", coverage_fail(np.array([[18.0, 18, 3.0]]), tau=tau, **G))
    ok("exactly tau from the face -> no flag", not coverage_fail(np.array([[18.0, 18, -2 + tau]]), tau=tau, **G))
    Gf = dict(G, direction=(1, 0, 0, 0, 1, 0, 0, 0, -1), origin=(0, 0, 36.0))   # flipped z
    ok("flipped axis handled", coverage_fail(np.array([[18.0, 18, 37.0]]), tau=tau, **Gf)
       and not coverage_fail(np.array([[18.0, 18, 18.0]]), tau=tau, **Gf))
    # final transform and band alternative
    P = dict(z_head=100.0, z_c7t1=40.0, z_pelvis=-300.0)
    ad = dict(head=True, pelvis=True); na = dict(head=False, pelvis=False)
    ok("head region, adopted -> final head", final_transform(120, P, ad)[0] == "head")
    ok("head region, not adopted -> final TA", final_transform(120, P, na)[0] == "TA")
    ok("on z_head (tie) -> trunk, TA", final_transform(100, P, ad) == ("TA", "trunk"))
    ok("band head side: final head -> alt TA", band_alternative(108, 105, 111, P, "head", ad) == (True, "TA"))
    ok("band trunk side: final TA -> alt head", band_alternative(92, 90, 95, P, "TA", ad) == (True, "head"))
    ok("band but head not adopted -> no alternative", band_alternative(92, 90, 95, P, "TA", na) == (True, None))
    ok("crossing plane counts as band", band_alternative(60, 60, 130, P, "TA", ad)[0])
    ok("no pelvic plane -> no pelvic band", band_alternative(-300, -310, -290, dict(P, z_pelvis=np.nan), "TA", ad) == (False, None))
    # resampling: identity on the same grid reproduces voxel sets; translation shifts them; coarse grid loses tiny comps
    arr = np.zeros((20, 20, 20), np.int32); arr[5:8, 5:8, 5:8] = 1; arr[12, 12, 12] = 2
    img = sitk.GetImageFromArray(arr); img.SetSpacing((2.0, 2.0, 2.0))
    I = sitk.Transform(3, sitk.sitkIdentity)
    r0 = sitk.GetArrayFromImage(sitk.Resample(img, img, I, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("identity resampling reproduces labels", np.array_equal(r0, arr))
    Tr = sitk.TranslationTransform(3, (4.0, 0.0, 0.0))          # resampling: q -> q + 4 mm (x) => label moves -2 voxels
    r1 = sitk.GetArrayFromImage(sitk.Resample(img, img, Tr, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("translation shifts the voxel set by -2 voxels in x", np.array_equal(np.argwhere(r1 == 1)[:, 2].min(), 3))
    coarse = sitk.Image([4, 4, 4], sitk.sitkUInt8); coarse.SetSpacing((10.0, 10.0, 10.0)); coarse.SetOrigin((0.0, 0.0, 0.0))
    rc = sitk.GetArrayFromImage(sitk.Resample(img, coarse, I, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("single voxel onto a coarse grid can vanish (zero survivor)", (rc == 2).sum() == 0)
    # overlap of two components mapped with different transforms: both sets kept
    A_ = sitk.GetArrayFromImage(sitk.Resample(img, img, I, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    arr2 = np.zeros_like(arr); arr2[3:6, 5:8, 5:8] = 3
    img2 = sitk.GetImageFromArray(arr2); img2.CopyInformation(img)
    B_ = sitk.GetArrayFromImage(sitk.Resample(img2, img, I, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    s1 = set(np.flatnonzero(A_ == 1)); s3 = set(np.flatnonzero(B_ == 3))
    ok("overlapping components from different images both kept", len(s1 & s3) > 0 and len(s1) == 27 and len(s3) == 27)
    # acceptance filter (added after the stopped run): coverage-flagged components are excluded, object-dtype flags
    T = pd.DataFrame(dict(case_id=["c"] * 4, voxels_resampled_final=[100, 100, 10, 100],
                          coverage_fov_final=pd.Series([True, False, False, np.nan], dtype=object),
                          image_check_mm=[20.0, 1.0, 20.0, 1.0]))
    chk_, bad_ = image_check_failures(T.iloc[:3], pd.Series({"c": 4.0}))
    ok("coverage-flagged component excluded from the image-level check", len(chk_) == 1 and len(bad_) == 0)
    T2 = T.iloc[:3].copy(); T2.loc[1, "image_check_mm"] = 9.0
    ok("unflagged component beyond one voxel is a failure", len(image_check_failures(T2, pd.Series({"c": 4.0}))[1]) == 1)
    ok("as_bool on object/strings/NaN", list(as_bool(pd.Series([True, False, np.nan, "True", "false"], dtype=object)))
       == [True, False, False, True, False])
    # v1.2 scenarios
    ok("scenario none = locked failed set", scenario_failed("train_0089", "none") == {"head", "neck", "pelvis"})
    ok("S1: train_0089 head and pelvis recovered, neck still failed", scenario_failed("train_0089", "S1") == {"neck"})
    ok("S1: NA head (train_0010) still failed", "head" in scenario_failed("train_0010", "S1"))
    ok("S1: recovered regions get Amendment 3a", scenario_adopted("train_0069", dict(head=True, pelvis=False), "S1")
       == dict(head=True, pelvis=True))
    ok("S2: NA head released, its neck still failed", scenario_failed("train_0010", "S2") == {"neck"})
    ok("S2: grade-1 head (train_0012) still failed", "head" in scenario_failed("train_0012", "S2"))
    ok("S4: necks released, heads/pelvises unchanged", scenario_failed("train_0089", "S4") == {"head", "pelvis"}
       and scenario_failed("train_0015", "S4") == set())
    ok("S2/S4 do not change transforms", scenario_adopted("train_0069", dict(head=True, pelvis=False), "S4")
       == dict(head=True, pelvis=False))
    # v1.3: S3b / S3c (clarification rules 1-7)
    PS = dict(z_head=100.0, z_c7t1=40.0, z_pelvis=-300.0)
    zz = np.array([120.0, 110.0, 105.0, 90.0, 60.0])               # 3 head, 2 neck
    ok("S3 failed points: head failed -> 3 of 5", s3_failed_points(zz, PS, {"head"}).sum() == 3)
    ok("S3 failed points: neck included (R1.6) -> 5 of 5 with head+neck",
       s3_failed_points(zz, PS, {"head", "neck"}).all())
    ok("S3 failed points: tie on z_head -> neck, not head",
       not s3_failed_points([100.0], PS, {"head"})[0] and s3_failed_points([100.0], PS, {"neck"})[0])
    ok("S3 failed points: no failed region -> none", not s3_failed_points(zz, PS, set()).any())
    ok("S3b: 60 % failed -> V_I", s3_decision(np.array([1, 1, 1, 0, 0], bool), "S3b") == (True, None, False))
    ok("S3b: exactly 50 % -> not V_I (> 50 % required)", s3_decision(np.array([1, 1, 0, 0], bool), "S3b")[0] is False)
    ok("S3b: 40 % -> not V_I, no split", s3_decision(np.array([1, 1, 0, 0, 0], bool), "S3b") == (False, None, False))
    ok("S3c: all failed -> whole V_I, no remainder", s3_decision(np.ones(4, bool), "S3c") == (True, None, False))
    ok("S3c: none failed -> unchanged, no split", s3_decision(np.zeros(4, bool), "S3c") == (False, None, False))
    v_, k_, s_ = s3_decision(np.array([1, 1, 0, 0, 0], bool), "S3c")
    ok("S3c: partial -> split, remainder = the non-failed points", v_ is False and s_ is True and list(k_) == [False, False, True, True, True])
    ok("S3 scenarios use the locked failed set", scenario_failed("train_0089", "S3b") == scenario_failed("train_0089", "S3c")
       == {"head", "neck", "pelvis"})
    # S3c resampling: removing failed points of one component leaves every other component's voxel set unchanged
    arr3 = np.zeros((20, 20, 20), np.int32); arr3[4:10, 5:8, 5:8] = 1; arr3[4:10, 9:12, 5:8] = 2
    im3 = sitk.GetImageFromArray(arr3); im3.SetSpacing((2.0, 2.0, 2.0))
    cg = sitk.Image([14, 14, 14], sitk.sitkUInt8); cg.SetSpacing((3.0, 3.0, 3.0)); cg.SetOrigin((0.5, 0.5, 0.5))
    Tt = sitk.TranslationTransform(3, (1.3, -0.7, 2.1))
    rA = sitk.GetArrayFromImage(sitk.Resample(im3, cg, Tt, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    arr3b = arr3.copy(); arr3b[4:7][arr3b[4:7] == 1] = 0                 # failed half of component 1 removed
    im3b = sitk.GetImageFromArray(arr3b); im3b.CopyInformation(im3)
    rB = sitk.GetArrayFromImage(sitk.Resample(im3b, cg, Tt, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("S3c: other components' resampled voxel sets unchanged", np.array_equal(rA == 2, rB == 2))
    ok("S3c: remainder resampled voxels are a subset of the whole component's", (rB == 1).sum() > 0
       and not ((rB == 1) & (rA != 1)).any() and (rB == 1).sum() < (rA == 1).sum())
    arr3c = np.zeros_like(arr3); arr3c[12, 12, 12] = 1                 # one-voxel remainder on a coarse grid
    im3c = sitk.GetImageFromArray(arr3c); im3c.CopyInformation(im3)
    rc3 = sitk.GetArrayFromImage(sitk.Resample(im3c, coarse, I, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("S3c: a remainder can vanish (zero survivor -> V_I coverage in P3)", (rc3 == 1).sum() == 0)
    # v1.4: perturbation geometry
    cc = np.array([10.0, -5.0, 3.0]); th = np.radians(30.0)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]]); tt = np.array([1.0, 2.0, 3.0])
    qq = np.array([[14.0, -5.0, 3.0], [0.0, 0.0, 0.0]])
    ok("zero perturbation leaves points unchanged", np.array_equal(perturb_points(qq, np.eye(3), np.zeros(3), cc), qq))
    ok("rotation about the centre: centre moves only by t", np.allclose(perturb_points(cc, Rz, tt, cc), cc + tt))
    Pi = inverse_perturbation(Rz, tt, cc)
    back = np.array([Pi.TransformPoint(tuple(x)) for x in perturb_points(qq, Rz, tt, cc)])
    ok("inverse affine recovers the unperturbed points", np.abs(back - qq).max() < 1e-12)
    Ta = sitk.Euler3DTransform(); Ta.SetRotation(0.3, -0.2, 0.5); Ta.SetTranslation((5.0, -1.0, 2.0))
    Tb = sitk.TranslationTransform(3, (7.0, 0.0, -3.0))
    xx = (1.0, 2.0, 3.0)
    ok("composite order: compose(A, B)(x) = A(B(x))",
       np.allclose(compose(Ta, Tb).TransformPoint(xx), Ta.TransformPoint(Tb.TransformPoint(xx)), atol=1e-12)
       and not np.allclose(compose(Ta, Tb).TransformPoint(xx), Tb.TransformPoint(Ta.TransformPoint(xx)), atol=1e-6))
    Ip = sitk.Transform(3, sitk.sitkIdentity)
    rp = sitk.GetArrayFromImage(sitk.Resample(img, img, compose(Ip, inverse_perturbation(np.eye(3), np.array([4.0, 0, 0]),
                                                                                         np.zeros(3))),
                                              sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok("perturbed resampling moves labels in the point-mapping direction (+4 mm = +2 voxels)",
       np.argwhere(rp == 1)[:, 2].min() == 7 and (rp == 1).sum() == 27)
    # v1.5: gate comparison with missing values
    F0 = pd.DataFrame(dict(a=["x", np.nan, "y"], b=pd.Series([True, np.nan, False], dtype=object), c=[1.0, np.nan, 2.0]))
    F1 = F0.copy(); F2 = F0.copy(); F2.loc[2, "a"] = "z"; F3 = F0.copy(); F3.loc[1, "b"] = False
    ok("gate comparison: identical frames with missing values -> no difference", frames_identical(F0, F1) == [])
    ok("gate comparison: a changed text value and a missing -> value change are detected",
       frames_identical(F0, F2) == ["a"] and frames_identical(F0, F3) == ["b"])
    F4 = F0.astype({"a": "string"})
    ok("gate comparison: string dtype with <NA> equals object with NaN", frames_identical(F0, F4) == [])
    # v1.6: independent nearest-neighbour re-derivation equals SimpleITK resampling (thin component, rotation, other grid)
    arr4 = np.zeros((20, 20, 20), np.int32); arr4[3:16, 10, 10] = 5; arr4[10, 3:8, 4] = 6
    im4 = sitk.GetImageFromArray(arr4); im4.SetSpacing((3.4, 3.4, 3.0)); im4.SetOrigin((1.0, -2.0, 0.5))
    g4 = dict(size=(26, 26, 22), spacing=(2.86, 2.86, 3.27), origin=(-3.0, -6.0, -2.0), direction=(1, 0, 0, 0, 1, 0, 0, 0, 1))
    ref4 = grid_image(g4)
    T4 = sitk.Euler3DTransform(); T4.SetCenter((30.0, 30.0, 30.0)); T4.SetRotation(0.06, -0.04, 0.08); T4.SetTranslation((1.3, -0.8, 2.2))
    r4 = sitk.GetArrayFromImage(sitk.Resample(im4, ref4, T4, sitk.sitkNearestNeighbor, 0, sitk.sitkInt32))
    ok4 = True
    for q4 in (5, 6):
        kji4 = np.argwhere(arr4 == q4)
        nat = np.array([im4.TransformIndexToPhysicalPoint(tuple(int(u) for u in x[::-1])) for x in kji4])
        Ti = T4.GetInverse(); mp4 = np.array([Ti.TransformPoint(tuple(p)) for p in nat])
        ok4 = ok4 and nn_rederive(T4, im4, arr4, q4, mp4, g4) == set(np.flatnonzero(r4.reshape(-1) == q4).tolist())
    ok("v1.6: independent nearest-neighbour re-derivation = SimpleITK resampling (thin components, rotation, other grid)", ok4)
    # v1.7: supplied mapping and fallback validation subset
    rngS = np.random.default_rng(7); pS = rngS.uniform(-200, 200, (6, 3))
    AS = sitk.Euler3DTransform(); AS.SetRotation(0.03, -0.02, 0.05); AS = np.array(AS.GetMatrix()).reshape(3, 3)
    bS = np.array([4.0, -7.5, 12.0]); tS = np.array([2.0, 1.0, -3.0])
    ok("v1.7: supplied mapping = A p + b; T_A = supplied - t (equal when t = 0)",
       np.allclose(supplied_points(pS, AS, bS), pS @ AS.T + bS) and np.allclose(supplied_points(pS, AS, bS) - tS, (pS @ AS.T + bS) - tS)
       and np.allclose(supplied_points(pS, AS, bS) - np.zeros(3), (pS @ AS.T + bS) - np.zeros(3)))
    RF = pd.DataFrame(dict(case_id=["train_0018", "train_0018", "train_0049", "train_0067", "train_0001", "train_0067"],
                           component_id=[1, 2, 1, 1, 1, 2], status_P2=["eligible", "eligible", "eligible", "V_S_candidate", "eligible", "eligible"],
                           final_transform=["TA", "pelvis", "TA", np.nan, "TA", "TA"], alt_transform=["", "", "head", np.nan, "", "TA"]))
    FSt = fallback_subset(RF)
    ok("v1.7: fallback subset = fallback patients, eligible, final TA, alternative absent or TA",
       list(zip(FSt.case_id, FSt.component_id)) == [("train_0018", 1), ("train_0067", 2)])
    for n_, g_ in res:
        print(f"  [{'PASS' if g_ else 'FAIL'}] {n_}")
    print(f"selftest: {sum(g_ for _, g_ in res)}/{len(res)} passed")
    return all(g_ for _, g_ in res)


# ------------------------------ run ------------------------------
def run(preflight, conn, scenario="none"):
    import SimpleITK as sitk
    t0 = time.time(); run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"P2_{'PREFLIGHT_' if preflight else ''}{'' if scenario == 'none' else scenario + '_'}conn{conn}_"
    md = {n: md5(d / n) if (d / n).exists() else "MISSING" for (d, n) in INPUTS}
    bad = {n: md[n] for (d, n), e in INPUTS.items() if md[n] != e}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(INPUTS)}/{len(INPUTS)} match")
    PP = pd.read_csv(OUT_DIR / "Audit1C_e_per_patient_20261031_141128.csv").set_index("case_id")
    TF = pd.read_csv(OUT_DIR / "Amendment3_transforms_20261102_073040.csv")
    FAt = pd.read_csv(OUT_DIR / "A3_visualQC_unblinded_20261108/A3_final_acceptance_per_region.csv", keep_default_na=False)
    FAt = FAt[FAt.repeat_of == ""] if "repeat_of" in FAt.columns else FAt
    FA = {(r.case_id, r.region): r.final_transform for r in FAt.itertuples()}
    PL = pd.read_csv(OUT_DIR / "SAP_R1_plane_table_20261111_214952.csv").set_index("case_id")
    C1 = pd.read_csv(OUT_DIR / P1_RUN / "P1_components_20261111_233936.csv")
    cases = sorted(C1.case_id.unique())
    if len(cases) != 100 or len(FA) != 200:
        sys.exit("STOP: expected 100 patients and 200 final-transform records")
    if preflight:
        cases = preflight_cases(cases, pd.read_csv(OUT_DIR / "A3_visualQC_KEY_do_not_open_20261105.csv"))
    odir = OUT_DIR / f"{tag}{run_id}"; odir.mkdir()
    rows, grids = [], []
    for i, cid in enumerate(cases, 1):
        r_, g_ = process(cid, conn, C1, PP, TF, FA, PL, odir, scenario)
        rows += r_; grids.append(g_)
        el = [r for r in r_ if r["status_P2"] == "eligible"]
        print(f"  [{i:3d}/{len(cases)}] {cid}: components {len(r_)}, eligible {len(el)}  ({time.time() - t0:.0f} s)")
    R = pd.DataFrame(rows); G = pd.DataFrame(grids)
    # ---- acceptance ----
    c1 = C1[(C1.case_id.isin(cases)) & (C1.connectivity == conn)]
    rec = (len(R) == len(c1) and (R.sort_values(["case_id", "component_id"])[["component_id", "voxels_native"]].values
                                  == c1.sort_values(["case_id", "component_id"])[["component_id", "voxels"]].values).all())
    E = R[R.status_P2 == "eligible"]
    chk, img_bad = image_check_failures(E, G.set_index("case_id").max_spacing_mm)
    rp = odir / f"P2_components_{run_id}.csv"; R.to_csv(rp, index=False)
    gp = odir / f"P2_grids_{run_id}.csv"; G.to_csv(gp, index=False)
    print(f"Reconciliation with P1 (component IDs and native voxels): {'PASS' if rec else 'FAIL'}")
    print(f"Image-level check (>= {MIN_VOX_IMAGE_CHECK} resampled voxels, no coverage flag): "
          f"{len(chk) - len(img_bad)}/{len(chk)} within the max PSMA PET voxel dimension")
    print(f"Composition max {G.composition_max_mm.max():.2e} mm; eligible components {len(E)}; "
          f"flags: registration {int(as_bool(E.VI_registration).sum())}, coverage {int(as_bool(E.coverage_fov_any).sum())}, "
          f"zero-survivor {int(as_bool(E['zero_survivor_final']).sum()) if 'zero_survivor_final' in E else 0}, "
          f"band {int(as_bool(E.band).sum())} (with alternative {int((E.alt_transform != '').sum())}), "
          f"components with collisions {int((E.get('collision_voxels', pd.Series(dtype=float)).fillna(0) > 0).sum())}")
    if not rec or len(img_bad):
        sys.exit(f"STOP: acceptance failed (reconciliation {rec}, image-level failures {len(img_bad)}); "
                 f"see {rp.name}; outputs in {odir.name} NOT accepted")
    man = dict(stage="P2 mapping, registration-QC flags, coverage, resampling", mode="PREFLIGHT" if preflight else "FULL",
               connectivity=conn, scenario=scenario, run=run_id, script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs=md,
               n_patients=len(cases), outputs={p.name: md5(p) for p in sorted(odir.glob("*")) if p.is_file()},
               params=dict(BAND_MM=BAND_MM, TOL_COMPOSE=TOL_COMPOSE, MIN_VOX_IMAGE_CHECK=MIN_VOX_IMAGE_CHECK, P1_RUN=P1_RUN),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=__import__("scipy").__version__, SimpleITK=sitk.Version_VersionString()),
               runtime_s=round(time.time() - t0), status="COMPLETE")
    mp = odir / f"P2_manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"done: {odir.name}; manifest MD5 {md5(mp)}")


def run_pert(pertdir, level, preflight, zero, supplied=False):
    import SimpleITK as sitk
    t0 = time.time(); run_id = datetime.now().strftime("%Y%m%d_%H%M%S"); conn = 26
    md = {n: md5(d / n) if (d / n).exists() else "MISSING" for (d, n) in INPUTS}
    bad = {n: md[n] for (d, n), e in INPUTS.items() if md[n] != e}
    if md5(PERT_NOTE[0] / PERT_NOTE[1]) != PERT_NOTE[2]:
        bad[PERT_NOTE[1]] = md5(PERT_NOTE[0] / PERT_NOTE[1])
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    md[PERT_NOTE[1]] = PERT_NOTE[2]
    if supplied:
        if md5(SUP_NOTE[0] / SUP_NOTE[1]) != SUP_NOTE[2]:
            sys.exit(f"STOP C1: MD5 mismatch {SUP_NOTE[1]}")
        md[SUP_NOTE[1]] = SUP_NOTE[2]
    if supplied:
        ks = [dict(level="SUPPLIED", k=1, R=np.eye(3), t=np.zeros(3), supplied=True)]
    elif zero:
        ks = [dict(level="ZERO", k=0, R=np.eye(3), t=np.zeros(3))]
    else:
        PT_, mdp = load_perturbations(pertdir); md.update(mdp)
        L = PT_[PT_.level == level].sort_values("k")
        if len(L) != 50 or list(L.k) != list(range(1, 51)):
            sys.exit(f"STOP: level {level} does not hold perturbations 1..50")
        ks = [dict(level=level, k=int(r.k), R=np.array([[r.R11, r.R12, r.R13], [r.R21, r.R22, r.R23], [r.R31, r.R32, r.R33]]),
                   t=np.array([r.t_x, r.t_y, r.t_z])) for r in L.itertuples()]
        if preflight:
            ks = ks[:2]
    print(f"C1 MD5: {len(md)}/{len(md)} match; perturbations: {len(ks)} "
          f"({'supplied-registration comparator' if supplied else 'zero gate' if zero else level})")
    PP = pd.read_csv(OUT_DIR / "Audit1C_e_per_patient_20261031_141128.csv").set_index("case_id")
    TF = pd.read_csv(OUT_DIR / "Amendment3_transforms_20261102_073040.csv")
    FAt = pd.read_csv(OUT_DIR / "A3_visualQC_unblinded_20261108/A3_final_acceptance_per_region.csv", keep_default_na=False)
    FAt = FAt[FAt.repeat_of == ""] if "repeat_of" in FAt.columns else FAt
    FA = {(r.case_id, r.region): r.final_transform for r in FAt.itertuples()}
    PL = pd.read_csv(OUT_DIR / "SAP_R1_plane_table_20261111_214952.csv").set_index("case_id")
    C1 = pd.read_csv(OUT_DIR / P1_RUN / "P1_components_20261111_233936.csv")
    cases = sorted(C1.case_id.unique())
    if len(cases) != 100 or len(FA) != 200:
        sys.exit("STOP: expected 100 patients and 200 final-transform records")
    if preflight:
        cases = preflight_cases(cases, pd.read_csv(OUT_DIR / "A3_visualQC_KEY_do_not_open_20261105.csv"))
    tag = (f"P2_{'PREFLIGHT_' if preflight else ''}SUPPLIED_conn{conn}_{run_id}" if supplied else
           f"P2_{'PREFLIGHT_' if preflight else ''}PERT{'ZERO' if zero else '_' + level}_conn{conn}_{run_id}")
    parent = OUT_DIR / tag; parent.mkdir()
    subs = {}
    for pk in ks:
        sub = parent / f"k{pk['k']:02d}"; sub.mkdir()
        rows, grids = [], []
        for cid in cases:
            r_, g_ = process(cid, conn, C1, PP, TF, FA, PL, sub, "none", pert=pk)
            rows += r_; grids.append(g_)
        R = pd.DataFrame(rows); G = pd.DataFrame(grids)
        c1 = C1[(C1.case_id.isin(cases)) & (C1.connectivity == conn)]
        rec = (len(R) == len(c1) and (R.sort_values(["case_id", "component_id"])[["component_id", "voxels_native"]].values
                                      == c1.sort_values(["case_id", "component_id"])[["component_id", "voxels"]].values).all())
        E = R[R.status_P2 == "eligible"]
        chk, img_bad = image_check_failures(E, G.set_index("case_id").max_spacing_mm)
        sk = f"k{pk['k']:02d}"
        if len(img_bad) and "nn_rederive_equal" in img_bad.columns:           # v1.6: exact re-derivation decides
            verified = img_bad[as_bool(img_bad["nn_rederive_equal"])]
            img_bad = img_bad[~as_bool(img_bad["nn_rederive_equal"])]
        else:
            verified = img_bad.iloc[0:0]
        verified[["case_id", "component_id", "voxels_native", "voxels_resampled_final", "image_check_mm",
                  "final_transform", "pct_volume_change_final"]].to_csv(sub / f"P2_imagecheck_exceptions_{sk}.csv", index=False)
        rp = sub / f"P2_components_{sk}.csv"; R.to_csv(rp, index=False)
        gp = sub / f"P2_grids_{sk}.csv"; G.to_csv(gp, index=False)
        cmax = float(G.pert_composition_max_mm.max())
        print(f"  {pk['level']} k{pk['k']:02d}: reconciliation {'PASS' if rec else 'FAIL'}; image check "
              f"{len(chk) - len(img_bad) - len(verified)}/{len(chk)} within tolerance, {len(verified)} above it but "
              f"reproduced exactly (listed), {len(img_bad)} failed; perturbed composition max {cmax:.1e} mm; coverage flags "
              f"{int(as_bool(E.coverage_fov_any).sum())}; zero-survivor {int(as_bool(E['zero_survivor_final']).sum())}  "
              f"({time.time() - t0:.0f} s)")
        if not rec or len(img_bad) or cmax > TOL_COMPOSE:
            sys.exit(f"STOP: acceptance failed in {sk}; outputs in {tag} NOT accepted")
        man = dict(stage="P2 perturbation", mode="PREFLIGHT" if preflight else "FULL", connectivity=conn, scenario="none",
                   perturbation=dict(level=pk["level"], k=pk["k"], R=pk["R"].tolist(), t=pk["t"].tolist(),
                                     supplied=bool(pk.get("supplied", False))),
                   run=run_id, parent=tag, script=Path(__file__).name, script_md5=md5(Path(__file__)),
                   outputs={p.name: md5(p) for p in sorted(sub.glob("*")) if p.is_file()}, status="COMPLETE")
        mp = sub / f"P2_manifest_{sk}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
        subs[sk] = md5(mp)
    if supplied:                                           # v1.7 gate: fallback validation subset
        P0 = OUT_DIR / PRIMARY_P2; r0 = PRIMARY_P2[-15:]
        A0 = pd.read_csv(P0 / f"P2_components_{r0}.csv", keep_default_na=False, na_values=[""])
        A0["alt_transform"] = A0["alt_transform"].fillna("")
        FS = fallback_subset(A0[A0.case_id.isin(cases)])
        A1 = pd.read_csv(parent / "k01" / "P2_components_k01.csv").set_index(["case_id", "component_id"])
        diffs = []
        for cid in sorted(FS.case_id.unique()):
            z0 = np.load(P0 / f"P2_voxelsets_{cid}.npz"); z1 = np.load(parent / "k01" / f"P2_voxelsets_{cid}.npz")
            for q in FS[FS.case_id == cid].component_id:
                for role in (0, 1):
                    s0 = np.sort(z0["linear_index"][(z0["component_id"] == q) & (z0["role"] == role)])
                    s1 = np.sort(z1["linear_index"][(z1["component_id"] == q) & (z1["role"] == role)])
                    if not np.array_equal(s0, s1):
                        diffs.append(f"{cid} c{q} role {role}")
                r1 = A1.loc[(cid, q)]
                if bool(as_bool(pd.Series([r1.coverage_fov_final])).iloc[0]) != bool(as_bool(pd.Series(
                        [A0[(A0.case_id == cid) & (A0.component_id == q)].coverage_fov_final.iloc[0]])).iloc[0]):
                    diffs.append(f"{cid} c{q} coverage")
        print(f"GATE (fallback validation subset reproduces the primary P2 voxel sets): {len(FS)} components in "
              f"{FS.case_id.nunique()} patients ({FS.volume_native_mL.sum():.1f} mL): "
              f"{'PASS' if not diffs else 'FAIL ' + str(diffs[:10])}")
        if diffs:
            sys.exit("STOP: fallback validation gate failed")
        fsp = parent / f"P2_fallback_subset_{run_id}.csv"; FS[["case_id", "component_id", "volume_native_mL"]].to_csv(fsp, index=False)
    if zero:                                               # gate: identical to the primary P2 run
        P0 = OUT_DIR / PRIMARY_P2; r0 = PRIMARY_P2[-15:]
        A0 = pd.read_csv(P0 / f"P2_components_{r0}.csv"); A1 = pd.read_csv(parent / "k00" / "P2_components_k00.csv")
        A0 = A0[A0.case_id.isin(cases)].reset_index(drop=True)
        diffs = frames_identical(A0, A1)
        for cid in cases:
            z0 = np.load(P0 / f"P2_voxelsets_{cid}.npz"); z1 = np.load(parent / "k00" / f"P2_voxelsets_{cid}.npz")
            if any(not np.array_equal(z0[a], z1[a]) for a in ("component_id", "role", "linear_index", "grid_shape_kji")):
                diffs.append(f"voxel sets {cid}")
        print(f"GATE (zero perturbation reproduces the primary P2 run {PRIMARY_P2}: {len(A0.columns)} columns, "
              f"{len(cases)} patients' voxel sets): {'PASS' if not diffs else 'FAIL ' + str(diffs[:10])}")
        if diffs:
            sys.exit("STOP: zero-perturbation gate failed")
    man = dict(stage="P2 supplied-registration comparator (parent)" if supplied else "P2 perturbation (parent)",
               mode="PREFLIGHT" if preflight else "FULL", level="SUPPLIED" if supplied else "ZERO" if zero else level,
               n_perturbations=len(ks), n_patients=len(cases), run=run_id, script=Path(__file__).name,
               script_md5=md5(Path(__file__)), inputs=md, sub_manifests=subs,
               fallback_subset=({"file": fsp.name, "md5": md5(fsp), "n_components": int(len(FS))} if supplied else None),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=__import__("scipy").__version__, SimpleITK=sitk.Version_VersionString()),
               runtime_s=round(time.time() - t0), status="COMPLETE")
    mp = parent / f"P2_manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"done: {tag}; manifest MD5 {md5(mp)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--conn", type=int, default=26, choices=(6, 18, 26))
    ap.add_argument("--scenario", default="none", choices=("none", "S1", "S2", "S4", "S3b", "S3c"))
    ap.add_argument("--pertdir"); ap.add_argument("--level", choices=("L1", "L2")); ap.add_argument("--pertzero", action="store_true")
    ap.add_argument("--supplied", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.supplied:
        run_pert(None, None, a.preflight, False, supplied=True)
        sys.exit(0)
    if a.pertzero or a.pertdir:
        if not a.pertzero and not a.level:
            sys.exit("--level L1|L2 required with --pertdir")
        run_pert(a.pertdir, a.level, a.preflight, a.pertzero)
        sys.exit(0)
    run(a.preflight, a.conn, a.scenario)
