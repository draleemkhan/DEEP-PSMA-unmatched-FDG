# Paper 3 — Step 2 Statistical Analysis Plan (consolidated)
**Date:** 11 November 2026 (Asia/Kuwait, UTC+03:00)
**Status:** LOCKED (Section 12). Draft `StepSAP_consolidated_DRAFT_20261111.md` (MD5 a723582bb27ceda7b780214954589d15); this version adds the ligand outcome (8.2, 11), the rotation qualification (8.1) and the lock record. Written **tumour-blind**: no PSMA–FDG correspondence, DMTVF or discordance result has been computed or examined. After lock, any change requires a dated amendment stating whether any PET correspondence result had been seen.
**Governing documents:** protocol v2.10 (MD5 2be3d63a9f130dbc83b7b03e11c37b70). Where this SAP fixes an item the protocol left provisional or [LOCK-AFTER-AUDIT], this SAP governs; the protocol is not edited.

## 0. Locked components incorporated by reference (unchanged)
| Component | File (00_protocol unless stated) | MD5 |
|---|---|---|
| Section R — registration handling, failed-region rule, precedence, S1–S3 | `StepSAP_section_R_registration_handling_20261109.md` | eecce0280370aca026eb1b2edf922efe |
| Section R.1 — neck regions (C2/C3, C7/T1), 29 failed regions in 23 patients | `StepSAP_section_R1_cervical_region_implementation_LOCKED_20261111.md` | 8fe2d04b39a906594d0c3370c2859cb5 |
| Plane table (fixed input) | `06_output\SAP_R1_plane_table_20261111_214952.csv` | db6fe2142e5c6fe7fe8ee6d8f9bc0f54 |
| Matching distance and τ (D1–D3; τ sensitivity definitions) | `StepSAP_tau_matching_tolerance_spec_LOCKED_20261111.md` | 285debf473d92a32c765f6f7ee45d412 |
| Soft-tissue displacement, ST1A, ST1B, ST2 | `StepSAP_soft_tissue_ST_spec_LOCKED_20261111.md` | 055199c4fd5d46039dcc65b6ab4f2416 |
| Amendment 3 (R6, transforms) | `Amendment3_regional_rigid_refinement_spec_LOCKED_v2_1_20261101.md` | 0a3ea46a0d8907a4f18a53899f2f423a |
| Trunk QC results and decisions (S4, V_I by cause) | `06_output\Trunk_QC_unblinding_note_20261111.md` | eac3000180decd9ddafccb8480b1a23a |
| Audit 1D report (closes Audit 1D) | `06_output\Audit1D_report_20261111.md` | cea5c01ec20d63202b9fddf91d6d73e1 |

## 1. Endpoint and volume accounting (protocol §5, §16.2; unchanged)
- V_FDG,total = V_S + V_I + V_M + V_U (native FDG volumes), one status per eligible component, in fixed order.
- **Primary endpoint:** DMTVF_ref = V_U / V_FDG,total, per patient. **Primary uncertainty interval:** DMTVF_ref to (V_U + V_I) / V_FDG,total. f_small = V_S / V_FDG,total reported separately. All 100 patients retained; patient-level fractions are never grounds for exclusion. Cohort summaries as specified in protocol §5.

## 2. Lesion units and eligibility (items 1–2)
- **Connectivity (item 1):** 26-connectivity primary; 6 and 18 as sensitivity analyses (each with its own eligibility, no transfer between partitions).
- **FDG eligibility (item 2):** component eligible if **≥ 0.1 mL AND ≥ 3 native voxels** (rule A, primary); **≥ 0.25 mL AND ≥ 3 voxels** (rule B, sensitivity). Applied in native FDG space before any transform. Basis: Audit 1D report §4–5 (rule A excludes 0.24% of pooled FDG volume under 26-connectivity; no patient without an assessable component).
- **PSMA match targets:** all supplied PSMA TTB-positive voxels; no minimum volume.

## 3. Mapping, regions and registration QC (Section R, R.1; unchanged)
Component transform by R6 (T_A centroid); failed regions (head, pelvis, neck) → V_I (registration QC) by the any-voxel intersection rule; 15 mm band at C2/C3 and L5/S1; no band at C7/T1.

## 4. Resampling, coverage and grid integrity (items 3–5)
- **Resampling (protocol §14.2–14.3):** unique FDG component IDs assigned in native space; the integer label image is resampled onto the **native PSMA PET grid** by nearest neighbour with each component's transform.
- **Resampling loss (item 4, author decision 11 Nov 2026):** a component is **V_I (coverage)** only if **zero voxels** survive resampling. For every component the native and transformed voxel counts and the **percentage volume change** are recorded; no shrinkage threshold. Potential label collisions are recorded (§14.3).
- **Field-of-view coverage (item 3, author decision 11 Nov 2026):** a component is **V_I (coverage)** if **any** of its native FDG voxel centres, point-mapped with its relevant transform(s) (Section R.5), lies **outside the PSMA PET field of view or within τ_i of its boundary**. The boundary is the outer face of the PSMA PET grid (voxel centres ± half spacing); distances in physical mm. Rationale: a PSMA match within τ_i could lie outside the image, so unmatched status there is not assessable.
- **Grid fingerprint (item 5):** deterministic hash of grid size, spacing, origin (rounded to 1 × 10⁻⁴ mm) and direction cosines (rounded to 1 × 10⁻⁶); recorded for every image passed between Python and MATLAB; a mismatch stops processing.

## 5. Matching (τ specification; unchanged)
Matched if d ≤ τ_i, d = 0 when any voxel is shared, otherwise the minimum distance between boundary-voxel centres (mm) on the native PSMA PET grid; τ_i = PSMA PET voxel diagonal; any PSMA-labelled voxel counts; one-to-many recorded.

## 6. Classification precedence (one status per component)
1. **V_S** — below eligibility (native space).
2. **V_I (registration QC)** — intersects a failed region (Section R, R.1).
3. **V_I (boundary band)** — R6 band rule.
4. **V_I (coverage)** — field of view (item 3) or zero voxels after resampling (item 4).
5. **V_M / V_U** — matching rule (Section 5).

V_I is reported by mutually exclusive cause (first applicable). Anatomically implausible matches are flagged, not reclassified (Section 8.3).

## 7. Sensitivity analyses (each against the primary analysis, one change at a time)
| Code | Change | Source |
|---|---|---|
| C6, C18 | 6- or 18-connectivity | item 1 |
| EB | eligibility rule B | item 2 |
| T-strict, T-perm | τ = max voxel dimension; τ = 2 × max voxel dimension | τ spec §τ6 |
| S1 | 9 grade-1 head/pelvis regions recovered with Amendment 3a | Section R.6 |
| S2 | 6 NA heads not failed | Section R.6 |
| S3b, S3c | volume-weighted (> 50 %) / split handling of failed regions | Section R.6 |
| S4 | failed necks released (re-enter at the band step) | unblinding note §5 |
| ST2 | visceral V_I stress test | soft-tissue spec |
| OVL | overlap-only matching (shared voxel; τ = 0), agreement analysis: % consistent, kappa, burden difference, patient status change (agreement ≠ accuracy) | protocol §5 item 5 |
Nested tolerances are reported as sensitivity, not robustness. No combined sensitivity scenarios.

## 8. Secondary analyses (items 7–10)
### 8.1 Registration perturbation (item 7; principal robustness analysis)
- Applied on top of each component's final transform; everything else unchanged.
- **Translation magnitudes:** 5.16 mm and 8.50 mm = cohort median and P90 of the patient-level held-out-bone median centroid residual under the adopted registration (`Audit1C_e_per_patient_20261031_141128.csv`, MD5 c66da788a6c73737163b4c09ed76f494; corrected value, supplied value for the three fallback patients).
- **Rotation magnitudes:** 2.29° and 4.47° = median and P90 of the fitted rotation of the provisionally adopted **pelvic** Amendment 3a corrections (n = 89; `Amendment3_fit_20261102_073040.csv`, MD5 dc74686ed0b024d8524668f9b9e16c94). Pelvic rather than head values (median 6.86°, P90 13.32°) because head rotation reflects neck posture, not trunk rigid error. These are **empirically motivated perturbation magnitudes, not measured residual uncertainty and not validated estimates of residual trunk rotational error** (author approval, 11 Nov 2026).
- Two levels (median: 5.16 mm + 2.29°; P90: 8.50 mm + 4.47°), each with **50** perturbations: translation direction uniform on the sphere, rotation axis uniform on the sphere, rotation about the centroid of the patient's eligible FDG volume (PSMA space); seed `rng(20261113, 'twister')`, perturbations generated once and stored with MD5.
- Outputs: change in DMTVF_ref per patient; per lesion, the percentage of perturbations in which it remains unmatched; "stably unmatched" = ≥ 95 % (exploratory); patient discordance-status changes. Limitation: rigid perturbation does not simulate local deformation.

### 8.2 PSMA uptake at FDG-lesion locations (item 8; exploratory; PROMISE-inspired)
- Question B of protocol §5 item 3. Lesion value: SUVmax of the **native PSMA PET** within the FDG component mapped into PSMA space (PSMA values never resampled).
- References from the PSMA CT labels: **blood pool** = SUVmean of the aorta mask eroded by 1 voxel (PET grid); **liver** = SUVmean of the liver mask eroded by 1 voxel. Liver (not spleen) because neither DEEP-PSMA ligand (⁶⁸Ga-PSMA-617, ¹⁸F-DCFPyL) is predominantly hepatobiliary.
- Classes: **absent** ≤ blood pool; **low** > blood pool and ≤ liver; **positive** > liver. Components with maximum extent ≤ 10 mm: **not classifiable (partial volume)**.
- Described as **PROMISE-inspired**, not PROMISE V2 scores (Seifert et al., Eur Urol 2023): PROMISE is visual, uses the parotid (not labelled by TotalSegmentator total task) and separates intermediate from high.
- **Ligand:** the per-patient PSMA ligand is **not recorded** in the supplied files (threshold.json holds only `suv_threshold`, checked for train_0001 PSMA and FDG; NIfTI header description fields empty, protocol §16.5) or in the published dataset documentation (dataset and dataset-format pages: "a mix of 68Ga-PSMA-617 and 18F-DCFPyL", no per-case identifier). Uptake classes are reported for the **mixed-ligand cohort**, labelled exploratory; ligand heterogeneity is a stated limitation. The ligand is **not inferred** from image properties.

### 8.3 Anatomical distribution (item 9; exploratory)
- Bone: spinal regions (cervical, thoracic, lumbar, sacral) by the locked planes and labels; other skeletal groups as labelled. Visceral organs as labelled. Unlabelled soft tissue is never automatically nodal.
- Nodal regions (rule-based, PROMISE V2 nomenclature where reliable): **supradiaphragmatic**; **retroperitoneal** (above the aortic bifurcation, i.e. the inferior end of the aorta label); **iliac zone — region unassigned** (below the aortic bifurcation). No miN versus miM1a distinction is made for pelvic nodes (common iliac nodes are miM1a in PROMISE V2 but cannot be separated reliably).
- Site-indeterminate and unassignable soft tissue reported as their own categories; within-patient clustering accounted for (protocol §5 item 4).
- **Anatomical consistency flag:** a matched pair is flagged if the FDG and PSMA components fall in different coarse compartments (bone; visceral organ; nodal/soft tissue). Flags are reported, not used to reclassify.

### 8.4 Voxel-level sensitivity (item 10)
"Spatially unmatched FDG tumour-volume fraction": assessable FDG voxels (of V_M or V_U components) whose centre lies more than τ_i from every PSMA voxel, divided by V_FDG,total. Voxels of V_S and V_I components never enter the numerator but remain in the denominator. Secondary methodological sensitivity only; not intralesional heterogeneity.

## 9. Software and reproducibility (item 6)
- **Python** (geometry, resampling, point mapping, region and coverage tests): Python 3.11 env `totalseg`; locked transform code unchanged.
- **MATLAB R2025a (25.1.0.2943329)**, Windows 11 Enterprise (build 26200); Statistics and Machine Learning Toolbox 25.1; Image Processing Toolbox 25.1. **All statistics in MATLAB.**
- Every MATLAB analysis script starts with `restoredefaultpath`, then adds only the project folders and any required dependency explicitly; `which` is checked for every required function (e.g. `fitglm`, `corr`, `prctile`) and the run stops if any resolves outside the MathWorks installation.
- Every run writes a manifest: software and toolbox versions, function paths, input and output MD5s, random seeds, grid fingerprints. The MATLAB licence number is not recorded.
- Every script: synthetic self-test, preflight, then full run; identical-output rerun (MD5) for deterministic steps.

## 10. Run order and gates
1. Lock this SAP (MD5).
2. Organ displacement measurement (CT only; soft-tissue spec Option B).
3. PET correspondence: eligibility → resampling → V_I (registration, band, coverage) → matching; self-tests and preflight first.
4. Primary endpoint and uncertainty interval.
5. Sensitivity analyses (Section 7), ST1A/ST1B/ST2, then secondary analyses (Section 8).

## 11. Items resolved at lock
- **Ligand metadata verification (11 Nov 2026):** completed; per-patient ligand not available (Section 8.2). No method changed.

## 12. Lock record
- Items 1–10 and all choices above approved by the author on 11 Nov 2026 (items 3(b), 4(a), 9(b); translation 5.16 / 8.50 mm; rotation 2.29° / 4.47°).
- **Locked on: 11 November 2026, 23:14 (Asia/Kuwait, UTC+03:00)**, before any PET correspondence, DMTVF or discordance result.
- This document's MD5: recorded externally.
