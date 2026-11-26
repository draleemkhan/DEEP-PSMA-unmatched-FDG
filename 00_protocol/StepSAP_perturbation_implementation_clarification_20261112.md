# Dated implementation clarification: registration perturbation analysis (SAP §8.1)
**Date:** 12 November 2026, 08:25 (Asia/Kuwait, UTC+03:00)
**Relates to:** consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) §8.1 (magnitudes 5.16 / 8.50 mm and 2.29° / 4.47°, two levels × 50 perturbations, translation direction and rotation axis uniform on the sphere, rotation about the centroid of the patient's eligible FDG volume in PSMA space, `rng(20261113, 'twister')`, generated once and stored with MD5; outputs: change in DMTVF_ref per patient, per-lesion % of perturbations unmatched, "stably unmatched" ≥ 95 % (exploratory), patient status changes; limitation: rigid perturbation does not simulate local deformation); primary P2/P2b/P3 runs (`P2_conn26_20261111_235603`, `P2b_conn26_20261112_000903`, `P3_conn26_20261112_001915`).
**Status:** written after the primary and all sensitivity analyses, **before** any perturbation is generated or any perturbation result computed. The locked §8.1 is not changed; this note fixes implementation details it does not state (author decisions Q1–Q7, one correction and two safeguards, 12 Nov 2026, 08:15).

## Q1–Q2 Generation (MATLAB, once)
- `rng(20261113, 'twister')`; then, in this order, 50 translation directions = rows of `randn(50, 3)` normalised to unit length, and 50 rotation axes = rows of a second `randn(50, 3)` normalised to unit length (uniform on the sphere).
- **Paired levels:** perturbation k uses the same direction and axis at both levels; only the magnitudes differ (L1 median: 5.16 mm and 2.29°; L2 P90: 8.50 mm and 4.47°).
- **Shared across patients:** perturbation k is the same direction/axis/magnitude for every patient; only the rotation centre is patient-specific.
- Rotation matrix by the axis–angle (Rodrigues) formula; each perturbation stored with its rotation matrix (9 elements) and translation vector at full precision (%.17g) in a CSV with MD5 and manifest. Directions are in the PSMA physical (LPS) frame used by P2 (uniformity does not depend on the frame).

## Q3 Rotation centre
c_i = the mean of the centres of all eligible (26-connectivity, rule A) native FDG voxels of patient i mapped with **T_A** (unperturbed): one fixed point per patient, independent of the regional transform of each component.

## Q4 What is recomputed
- Perturbed mapping of any point already mapped by a transform used for a component: **P_k(q) = R_k (q − c_i) + c_i + t_k**, applied identically to **every transform used** (T_A, the final transform, and the band alternative).
- **Recomputed under each perturbation:** coverage (FOV and τ margin, margin at the reference τ), resampling (nearest neighbour, through the composed resampling transform x ↦ RT(P_k⁻¹(x))), zero survivors, both band categories, matching (reference τ), classification.
- **Kept from the primary analysis (unperturbed):** eligibility, the R6 transform assignment, the relevant transforms, band candidacy and the alternative transform, the registration-QC failed-region flag. These determine which registration applies, not where it lands.
- **Denominator (safeguard):** DMTVF_ref,ik = V_U,ik / V_FDG,total,i, with V_FDG,total,i the native FDG tumour volume, identical to the primary analysis for every perturbation (checked); only the classifications and the numerator change. Native-volume accounting unchanged.

## Q5 Lesion stability
Denominator = all 50 perturbations of a level; a component counts as unmatched in a perturbation only if it is V_U. Reported per component and level: n V_U, n V_M, n V_I (by cause), % unmatched. **Stably unmatched = V_U in ≥ 48 of 50** (≥ 95 %, exploratory). Reported for every primary V_U component; primary V_M components that become V_U in any perturbation are also reported.

## Q6 Patient summaries
Per patient and level: ΔDMTVF_ref,ik = perturbed − primary over the 50 perturbations (median, IQR, minimum, maximum; numbers of perturbations with an increase, a decrease and no change); status = DMTVF_ref > 0, % of perturbations whose status differs from the primary analysis, with 0 → positive and positive → 0 counted separately. Cohort: per perturbation, cohort V_U and median DMTVF_ref; summarised over the 50 perturbations.

## Correction (author): no direction assumed
Perturbations can **increase or decrease** DMTVF_ref (an unmatched component may become matched, a matched one unmatched, and components may enter or leave V_I). Positive and negative changes are reported without assuming a direction.

## Q7 Implementation and validation
- MATLAB `PERT_generate_v1_0.m` (generation; self-tests: orthonormality and det = 1, recovered angle and translation magnitude, axis invariance, pairing across levels, uniform-sphere moments).
- P2 v1.4 perturbation mode (one call per level, looping over the 50 perturbations; QC images not written), P2b v1.4, P3 v1.11 (loop over a level; per-perturbation tables and summaries).
- **Gates:** (1) a zero perturbation (R = I, t = 0) through the new P2 code must reproduce the primary P2 component table and voxel sets exactly, and through P3 the primary P3 tables; (2) per patient, perturbation and transform, composition check: the composed resampling transform applied to perturbed mapped points returns the native points within 10⁻⁴ mm; (3) synthetic self-tests; (4) two-patient preflight; then each level in full.

## Interpretation (author safeguards)
- **Stability ≠ accuracy:** a component unmatched in 48–50 perturbations is geometrically stable under this perturbation model; it is not thereby shown to be biologically FDG-positive / PSMA-negative.
- The two levels are **empirically motivated displacement magnitudes, not established probability distributions of registration error**; rigid perturbation does not simulate local deformation.
- Substantial instability would be a valid scientific result.
