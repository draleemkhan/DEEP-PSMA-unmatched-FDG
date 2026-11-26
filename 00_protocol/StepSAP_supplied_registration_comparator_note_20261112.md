# Dated note: supplied-registration comparator (POST-HOC, documented before results)
**Date:** 12 November 2026, 12:25 (Asia/Kuwait, UTC+03:00)
**Approved by the author:** 12 Nov 2026, 12:20 (decisions 1–6 and the conditional-comparison safeguard).
**Status:** **POST-HOC, DESCRIPTIVE.** Not part of the locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d). Proposed and documented after the primary, sensitivity, perturbation, uptake and cross-tabulation results had been seen; **no result under the supplied registration has been computed or viewed.** The locked primary pipeline and all results are unchanged. **Single comparator run; no preferred result will be selected after seeing it; no new biological claim.** All outputs below will be reported whatever they show.

## 1. Question
How does using the **supplied (dataset) rigid registration** instead of the **adopted registration strategy** (T_A = supplied rigid + translation correction, plus Amendment 3a regional head/pelvis corrections where adopted) change measured FDG–PSMA discordance, **conditional on the adopted pipeline's registration-QC decisions**?

This is a **conditional registration-choice comparison, not a comparison of two equally QC-qualified pipelines.** The comparison reflects the **complete registration strategy** (translation correction and regional refinement together), not the translation correction alone. A difference shows sensitivity of the endpoint to the registration choice; it does not establish which registration produces biologically correct correspondence.

## 2. Comparator definition
- Mapping of every eligible FDG component: the supplied rigid transform S(p) = A p + b (the transform whose residuals were audited in Audit 1C-e), for **every** transform the component uses (final and, where present, alternative). No translation correction, no regional refinement.
- Resampling: the supplied composite (the supplied rigid alone; the same resampling transform as used for the three fallback patients), nearest neighbour, label values only; composition checked per component (native points recovered within 1 × 10⁻⁴ mm).
- **Frozen at component level from the primary analysis:** eligibility, registration-QC status (V_I registration), R6 final-transform label, relevant transforms, band candidacy and alternative-transform label. Registration-QC flags are therefore **not** recomputed under the supplied mapping, so the withheld set is identical in both arms and the comparison isolates the registration effect (same convention as the perturbation analysis).
- **Recomputed under the supplied mapping:** coverage (FOV / τ margin), resampled voxel sets, zero survivors, matching (reference τ_i, any-voxel rule), classification and accounting; denominator = native V_FDG,total (unchanged).
- **Boundary band:** the anatomical planes and band candidacy are retained, but the supplied transform is applied on both sides, so final and alternative voxel sets are identical and **band V_I = 0 by construction.** A zero band volume is therefore not a finding; it is checked (run stops otherwise).
- Precedence unchanged: V_S → V_I registration → V_I band → V_I coverage → V_M / V_U.

## 3. Validation
- **Primary reproduction gate:** the classifier first reproduces the primary P3 tables from the primary inputs (generic loader), as in all previous comparisons.
- **P2 code change gate:** the revised P2 script must reproduce the primary P2 run exactly in zero-perturbation mode (components and voxel sets) before the comparator is run.
- **Fallback validation subset:** in train_0018, train_0049 and train_0067 the adopted trunk transform equals the supplied transform. Components of these patients whose primary final transform is T_A **and** whose alternative transform is absent or T_A must reproduce the primary analysis exactly (P2: identical final and alternative voxel sets; P3: identical status, cause and final distance). Components of these patients using a regional (head/pelvis) transform are excluded from the gate. **The size of the subset is reported**; if it is small, the gate is weak and this is stated.
- Reconciliation (V_S + V_I + V_M + V_U = V_FDG,total) and denominator identity per patient.

## 4. Outputs (both arms, reported regardless of result)
- Cohort V_M, V_U, V_I (registration, band, coverage) and V_S (mL).
- Median patient DMTVF_ref with IQR; patients with DMTVF_ref > 0, > 1 %, > 5 %, > 10 %; median upper bound.
- Per-patient change (supplied − adopted): number with increase, decrease, no change; median, IQR, minimum and maximum change; patient status changes (DMTVF_ref > 0) in both directions.
- Component-level agreement on components assessable in both arms (M/M, M/U, U/M, U/U), % consistent, Cohen's κ with patient-cluster bootstrap 95 % CI (2,000 replicates, seed 20261114; CI reported only if ≥ 90 % of replicates valid); V_I transitions reported separately.
- **QC-limitation output (made visible, not buried):** (a) number and native FDG volume of eligible components whose primary final transform is a regional Amendment 3a correction (head or pelvis), i.e. whose region was accepted only after regional correction, with their primary status distribution; (b) explicit statement that trunk registration QC was also assessed under T_A, not under the supplied registration, so **no region was QC-validated under the supplied registration.**
- Fallback-subset size and gate result.

## 5. Interpretation boundaries
- Held-out skeletal residuals (Audit 1C-e: median 8.00 mm supplied vs 5.16 mm adopted; chest-wall structures only) support better skeletal alignment of the adopted registration; they do not validate lesion-level correspondence of either arm.
- With any-voxel matching, misalignment can both remove true matches and create accidental ones; no direction of change is assumed.
- Results will be described as the effect of the complete registration strategy, conditional on the adopted QC decisions.

## 6. Implementation
- `P2_mapping_coverage_resampling_v1_7.py` (adds `--supplied`; all other modes unchanged), `P2b_export_voxelsets_csv_v1_7.py` (accepts P2 v1.7), `P3_matching_classification_v1_14.m` (adds mode `supplied`; accepts P2/P2b v1.7). New versioned files; earlier versions unchanged.
