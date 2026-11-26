# Dated note: exploratory cross-tabulation of perturbation stability and PSMA uptake (post-hoc, documented before results)
**Date:** 12 November 2026, 09:40 (Asia/Kuwait, UTC+03:00)
**Approved by the author:** 12 Nov 2026, 09:33.
**Status:** **POST-HOC and EXPLORATORY.** Not part of the locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d). Documented before results: at the time of writing, no full-level perturbation classification had been computed or viewed (only the two-patient L1 preflight summaries, mechanics checks, had been seen), and no PSMA uptake value had been computed. The primary endpoint, the locked classifications, thresholds and analyses are not modified.

## 1. Basis
- **Individually prespecified (locked SAP):** the component-level classification (primary analysis); the registration perturbation analysis (§8.1; implementation clarification MD5 5f280df40970eec504bb4f08d421d69d; stably unmatched = V_U in ≥ 48/50 perturbations, exploratory); the PSMA uptake classes (§8.2; absent ≤ blood pool, low > blood pool and ≤ liver, positive > liver; components with maximum extent ≤ 10 mm not classifiable).
- **Not prespecified:** their combined cross-tabulation, defined here. It is descriptive and exploratory; it does not modify the primary endpoint; **no new thresholds, exclusions or biological classifications** will be introduced after the results are examined.

## 2. Units and measurements
- Units: components classified **V_U in the primary analysis** (the locked units of §8.2).
- PSMA uptake class: as locked in §8.2, sampled at the **unperturbed (primary) mapped location** only. (Implementation details of §8.2 will be fixed in their own dated note before that analysis is run.)
- Perturbation counts, **per level (L1, L2) separately**, for each component: **N_U + N_M + N_I = 50**, with N_I split by cause (boundary band; coverage). Registration-QC status is kept from the primary analysis under perturbation, so V_I can arise only from the band or coverage steps. **The two levels are never combined and no composite stability score is created.**

## 3. Categories (per level)
- **Stably unmatched:** N_U ≥ 48 (actual N_U, N_M, N_I still reported).
- **Unstable:** N_U < 48, with two **non-exclusive** flags: **sometimes matched** (N_M > 0) and **sometimes indeterminate** (N_I > 0). A component with both is reported with both flags, not forced into one category.

## 4. Interpretation framework (imaging evidence, not validated biological diagnoses)
| Perturbation result (per level) | PSMA uptake class | Interpretation |
|---|---|---|
| Stably unmatched | Absent (≤ blood pool) | Imaging pattern **most consistent with** FDG+/PSMA-low disease (not confirmed biological discordance) |
| Stably unmatched | Low (> blood pool, ≤ liver) | Persistent spatial discordance with measurable low PSMA uptake |
| Stably unmatched | Positive (> liver) | Persistent annotation discordance despite appreciable PSMA signal (the signal may be unannotated tumour, physiological uptake, spill-in, registration error or another lesion) |
| Unstable: sometimes matched | Any | Correspondence depends on spatial alignment |
| Unstable: sometimes indeterminate | Any | Assessability depends on spatial alignment |
| Any | Not classifiable (maximum extent ≤ 10 mm) | Uptake-based interpretation withheld |

## 5. Reporting
- Each cell reported by number of components, native FDG volume and number of patients (components cluster within patients), separately for L1 and L2; flags and counts as above.
- **Stability ≠ accuracy:** stable means stable under this perturbation model (rotation about the patient's FDG centroid makes a component's displacement depend on its distance from that centre; 50 sampled directions cannot exclude PSMA disease in unsampled directions); it does not establish an empty PSMA neighbourhood of known radius.
- PSMA SUV ≥ 3 is the dataset's segmentation threshold, not a clinical criterion; the blood-pool/liver classes are clinically interpretable reference categories, not VISION or TheraP eligibility classifications; mixed-ligand cohort.
