# Step 2 SAP — Soft-tissue displacement and sensitivity analyses ST1A, ST1B, ST2: specification
**Date:** 11 November 2026 (Asia/Kuwait, UTC+03:00)
**Status:** LOCKED (SP8). Draft `StepSAP_soft_tissue_ST_spec_DRAFT_20261111.md` (MD5 8c3874628b9a995ff47f6c78abb0834c); this version adds two author clarifications (11 Nov 2026: lung-lobe fragmentation, SP3; warning-independent separation check, SP5) and the lock record. Tumour-blind: no PET correspondence, DMTVF or discordance result has been examined. The organ measurements (CT only) are run **after the full Step 2 SAP lock**.
**Builds on (unchanged):** trunk QC spec LOCKED §9; trunk QC unblinding note §5 (MD5 eac3000180decd9ddafccb8480b1a23a); SAP section R, R.1 (LOCKED); τ specification (LOCKED, MD5 285debf473d92a32c765f6f7ee45d412); protocol v2.10.

## SP0 Purpose
Skeletal QC establishes bony correspondence only. These analyses describe residual soft-tissue displacement under the frozen registration and test how strongly the results depend on correspondence in deformable organs. They are **associational and exploratory**; none changes registration acceptance, V_I in the primary analysis, or the primary endpoint.

## SP1 Structures
Six structures, TotalSegmentator total task (single label map, so mutually exclusive): **liver; spleen; right kidney; left kidney; right lung = union of the right upper, middle and lower lobes; left lung = union of the left upper and lower lobes.**

## SP2 Independence audit (completed, read-only, 11 Nov 2026)
- `Audit1C_e_translation_refinement.py` (MD5 ce038cb661670d673bdc7dbbb448f10b) computes each patient's Amendment 2 translation t as the component-wise median offset over the **fitting** set, taken verbatim from `Audit1C_c_structure_level_20261029_212531.csv` (MD5 77469a3a828691e41f43c2c92ed04eda): 21 bones (T1–T12, L1–L5, S1, sacrum, left and right hip). The union of all 100 per-patient fitting lists (`Audit1C_e_per_patient_20261031_141128.csv`, MD5 c66da788a6c73737163b4c09ed76f494) is the same 21 bones.
- **None of the six structures contributed to any patient's translation** (all are in the held-out soft set).
- **Disclosure (partial selection dependence):** the held-out soft set, including these organs, entered the cohort-level Amendment 2 adoption criterion C2 (held-out soft median paired difference ≤ +0.5 mm). The measurements are therefore not in-sample fitting residuals but are not wholly independent of the decision to adopt Amendment 2. Reported; no subset restriction.

## SP3 Displacement measurement (CT only; Option B: one new run for all six)
- **Transform:** the frozen Amendment 2 (T_A; supplied transform for train_0018, 0049, 0067). No organ-specific or new registration.
- **Masks:** whole structures, no cropping. FDG labels resampled onto the PSMA CT grid by nearest neighbour with the Amendment 2 resampling transform; lungs as unions of their lobe labels on each CT before any computation.
- **Metrics** (functions copied verbatim from Amendment 1 / Audit 1C-e):
  - centroid distance (mm): |T_A(c_FDG) − c_PSMA|, centroids of the native masks (FDG centroid point-mapped);
  - **ASSD (mm, principal descriptive metric)** = Amendment 1 `surface_metrics` MSD: mean of the pooled distances from each surface to the other (both directions);
  - HD95 (mm) = pooled 95th percentile of the same distances.
- **Eligibility for displacement:** on **both** CTs (Audit 1B, `Audit1B_label_level_20261029_134401.csv`, MD5 896a1681e90da4fe793a3a449d737fde), PRESENT_COMPLETE **and** not fragmented (frac_largest_component ≥ 0.95). For a lung, **every** contributing lobe must meet both conditions; completeness and fragmentation are assessed **per lobe**, and the lung union itself is **not** required to form a single connected component. Independent of PET.
- **Validation (required; any failure stops):** for liver, spleen, right and left kidney, the new values must reproduce Audit 1C-e (`Audit1C_e_structure_level_20261031_141128.csv`, corrected_centroid_mm, corrected_MSD_mm, corrected_HD95_mm) to ≤ 0.01 mm **for every record with identical masks, transform, eligibility and numerical method**. Records not expected to match (e.g. different eligibility because of the fragmentation filter) are listed with the reason; any unexplained difference is investigated before use, never silently accepted. Composition check of the transform (as in the earlier scripts).
- **Output:** one table, patient × structure: eligibility (and reason), centroid distance, ASSD, HD95, number of surface points; MD5 and manifest. Python (geometry); statistics in MATLAB.

## SP4 ST1A — patient-level association (exploratory)
- For each structure (6) and metric (3): Spearman ρ between displacement and the patient's **primary-analysis DMTVF_ref**, over patients with that structure eligible.
- Patient-level summary per metric: median over the patient's eligible structures, reported **with the number of contributing structures**; no composite across metrics.
- 95% CI: patient bootstrap, 2,000 resamples, percentile method, fixed seed (MATLAB `rng(20261111, 'twister')`). Reported with n evaluable patients. **No hypothesis tests, no p-values.**

## SP5 ST1B — lesion-level association (exploratory)
- **Units:** eligible FDG components (primary analysis) classified **V_M or V_U**; V_I and V_S excluded. Each component counts once (no volume weighting).
- **Organ assignment:** a component is assigned to a structure if **> 50 %** of its native FDG voxels lie in that structure's native FDG mask (for lungs: the union of lobes); otherwise unassigned and not analysed. The assigned structure must be eligible for displacement (SP3).
- **Outcome:** unmatched (V_U) = 1, matched (V_M) = 0.
- **Model:** logistic regression of the outcome on the assigned structure's **ASSD (per mm)**, adjusted for structure (categorical). MATLAB `fitglm` (binomial, logit).
- **Inference:** patient-cluster bootstrap (resample patients with all their components), 2,000 resamples, percentile 95% CI, fixed seed (`rng(20261112, 'twister')`). Report the odds ratio per mm and its 95% CI; **no p-values**.
- **Failed fit (objective; does not rely on MATLAB warnings alone):** a fit fails if any of the following holds: (i) MATLAB iteration-limit or separation warning; (ii) any non-finite coefficient or standard error; (iii) any structure category whose components all have the same outcome (quasi-complete separation by the adjustment term); (iv) any fitted probability < 1e-8 or > 1 − 1e-8 (numerical separation).
- **Objective fallback (author decision, 11 Nov 2026):** if (a) the full-data fit fails, **or** (b) more than **10 %** of the 2,000 bootstrap fits fail by the same criteria (i)–(iv), the model is not reported; descriptive counts and medians of ASSD by outcome and structure are reported instead. The 10 % threshold is a pragmatic prespecified rule, not a statistical standard. No minimum sample size is imposed.
- Organ displacement is a patient-and-structure-level exposure; components in the same structure share it (stated as a limitation).

## SP6 ST2 — visceral V_I stress test (exploratory)
- **Masks:** native FDG masks of the six structures (lungs: union of lobes), each lobe or organ included if **present** (any Audit 1B category other than absent) and **not fragmented** (frac_largest_component ≥ 0.95) on the **FDG CT**; truncated structures included. Independent of measured displacement. Structures available for ST2 recorded separately from those eligible for displacement.
- **Rule:** an eligible FDG component with **any voxel** in such a mask is V_I.
- **Precedence:** V_S → V_I (registration QC, attributed first) → V_I (ST2 organ) → V_I (boundary band) → V_I (coverage) → V_M / V_U. Denominator unchanged.
- **Reporting:** V_I by mutually exclusive cause; the additional V_I from ST2; DMTVF_ref, uncertainty interval and patient-level changes vs primary. Conservative only for components intersecting available masks.

## SP7 Interpretation
ST1A/ST1B assess association between measured displacement and correspondence status; ST2 estimates the consequence of withholding correspondence decisions in visceral organs. Neither establishes whether an individual FDG–PSMA mismatch is biological.

## SP8 Lock record
- Approved by the author on 11 Nov 2026 (decisions 1–5 and the two clarifications above).
- **Locked on: 11 November 2026, 22:16 (Asia/Kuwait, UTC+03:00)**, before any PET correspondence, DMTVF or discordance result, and before the organ measurements.
