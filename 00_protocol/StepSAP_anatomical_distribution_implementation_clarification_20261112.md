# Dated implementation clarification: anatomical distribution of unmatched FDG burden and anatomical consistency flag (SAP §8.3)
**Date:** 12 November 2026, 12:55 (Asia/Kuwait, UTC+03:00)
**Approved by the author:** 12 Nov 2026, 12:45–12:47 (A1–A8 with two safeguards).
**Relates to:** locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d) §8.3 (anatomical distribution; anatomical consistency flag); protocol v2.10 §4.4, §5 items 4–5 (site-indeterminate; unlabelled soft tissue never automatically nodal). **This completes an analysis required by the locked SAP; it is not a new exploratory experiment.** Implementation detail only; no locked definition changed.
**Disclosure:** written after the primary, sensitivity, perturbation, uptake, cross-tabulation and supplied-registration results had been seen; **no anatomical assignment of any FDG or PSMA component has been computed or viewed.**

## 1. Units (A1)
- **Distribution:** only the **217 components classified V_U in the primary analysis** (`P3_conn26_20261112_001915`). No matched-component distribution is produced.
- **Anatomical consistency flag (A7):** a **separate descriptive audit** of matched pairs (primary V_M components and their accepted PSMA components). It is reported separately and never added to the unmatched-burden distribution.

## 2. Anatomical labels (A2)
- Native **FDG CT** TotalSegmentator labels (`<case>\FDG\totseg_24.nii.gz`, official 117-class dictionary as verified in Audit 1B; label IDs taken by name from `Audit1B_label_level_20261029_134401.csv`, MD5 896a1681e90da4fe793a3a449d737fde).
- Points: physical centre of every native voxel of the component (P1 label image, 26-connectivity, FDG PET grid) → nearest FDG CT voxel by floor(x + 0.5) (as ST2a); outside the CT grid → "outside CT".
- Label groups:
  - **Bone (spinal):** vertebrae_C1 … vertebrae_S1 (26–50), sacrum (25).
  - **Bone (other skeletal groups):** skull (91); ribs and costal cartilages (92–115, 117); sternum (116); clavicles (73–74); scapulae (71–72); humeri (69–70); pelvic bones (hip_left/right, 77–78); femurs (75–76).
  - **Visceral organs (each reported separately):** spleen; kidneys (incl. kidney cysts); gallbladder; liver; stomach; pancreas; adrenal glands; lungs (all lobes); oesophagus; trachea; thyroid; small bowel; duodenum; colon; urinary bladder; prostate; heart; brain.
  - **Non-organ soft tissue:** every other CT voxel inside the grid (vessels, muscles, spinal cord and unlabelled tissue).

## 3. Assignment (A3; strictly > 50 % of native voxel centres)
In this order, fractions over all native voxel centres of the component:
1. Fraction inside the FDG CT grid ≤ 0.5 → **site-indeterminate: anatomy unavailable**.
2. Bone fraction > 0.5 → **bone**. Within bone: spinal fraction > 0.5 → spinal (region by §4); otherwise one other skeletal group > 0.5 → that group; otherwise **bone, mixed skeletal groups**.
3. One visceral organ > 0.5 → **that organ**.
4. Non-organ soft-tissue fraction > 0.5 → **soft tissue**, location by §5.
5. Otherwise (no tissue class > 0.5) → **site-indeterminate: anatomically ambiguous** (mixed tissue).

Safeguard: "anatomy unavailable" and "anatomically ambiguous" (site-indeterminate) are kept separate from each other and from "unassignable soft tissue" (§5).

## 4. Spinal regions: locked planes applied to native components (A4)
- The **locked plane table** `06_output\SAP_R1_plane_table_20261111_214952.csv` (MD5 db6fe2142e5c6fe7fe8ee6d8f9bc0f54) is used unchanged; no new plane system is introduced.
- The planes are defined in **PSMA space**. They are applied to native FDG components **exactly as the locked R6/S3 rules do**: every native voxel centre is mapped with **T_A** (the adopted trunk transform; the supplied transform for train_0018, 0049, 0067) and the component's craniocaudal coordinate is the **mean z of the T_A-mapped centres** (the R6 convention).
- Regions: **cervical** z ≥ z_C7T1; **thoracic** z_T12L1 ≤ z < z_C7T1; **lumbar** z_L5S1 ≤ z < z_T12L1; **sacral** z < z_L5S1. A component exactly on a plane goes to the superior region (consistent with the locked closed upper boundary of the neck region at C7/T1). Cervical includes C1–C2 (the spine is not divided at C2/C3).
- SUBSTITUTED planes are used as locked. If a plane needed for the decision is UNDEFINED (T12/L1: 1 patient; L5/S1: 3 patients), the component is reported as **spinal, region unassigned (plane undefined)**; no substitute plane is created.

## 5. Soft-tissue location (A5, A6; location only, never "nodal")
- **Safeguard:** TotalSegmentator does not establish that unlabelled soft tissue is a lymph node. Soft-tissue components are described by **location relative to the predefined nodal regions**, never classified as nodal metastases. Labels: "soft tissue — supradiaphragmatic location", "soft tissue — retroperitoneal level", "soft tissue — iliac zone (region unassigned)".
- **Supradiaphragmatic (operational proxy, not the anatomical diaphragm):** T_A-mapped mean z ≥ z_T12L1 (locked plane, as §4) and ≤ z_head_C2C3. Above z_head_C2C3 → **unassignable soft tissue (head)**. T12/L1 UNDEFINED → **soft tissue, location unassigned (plane undefined)**.
- **Below T12/L1:** aortic bifurcation = inferior end of the aorta label (52) in the native FDG CT (minimum z of aorta voxel centres). Native component mean z ≥ bifurcation → **retroperitoneal level**; < bifurcation → **iliac zone (region unassigned)**, unless below the inferior end of the pelvic bones (minimum z of hip_left/right voxel centres, native FDG CT) → **unassignable soft tissue (below pelvis)**. Aorta label absent → **soft tissue, infradiaphragmatic level unassigned**; pelvic-bone labels absent → the lower limit is not applied and this is recorded.
- The T12/L1 decision uses the locked PSMA-space plane (T_A mapping); the bifurcation and pelvic limits use native FDG CT labels. Both are stated in the Methods.
- No miN / miM1a distinction (SAP §8.3).

## 6. Anatomical consistency flag (A7; separate descriptive audit)
- Pairs: each primary V_M component and each PSMA component listed as matched to it in the primary run.
- Coarse compartment of each side (SAP: bone; visceral organ; nodal/soft tissue) by the same > 50 % rule (§3; site-indeterminate → **not evaluable**). FDG side: native FDG CT labels. PSMA side: PSMA TTB components (26-connectivity, component IDs in raster order of first voxel = MATLAB bwconncomp column-major order, as P1), native voxel centres → native **PSMA CT** labels (`<case>\PSMA\totseg_24.nii.gz`).
- A pair is **flagged** if both sides are evaluable and their compartments differ. Flags are reported (pairs; FDG components with ≥ 1 flagged pair: number, native mL, patients; not evaluable counts), **not used to reclassify**.
- ID check: for matched FDG components with final distance 0, at least one listed PSMA component must share a voxel with the FDG final mapped voxel set (verifies the PSMA component-ID convention; run stops otherwise).

## 7. Reporting (A8)
Per category: number of components, native FDG volume (mL) and number of patients (clustering acknowledged by reporting patients). Spinal regions, skeletal groups, organs and soft-tissue locations separately. No statistical tests. Reconciliation: 217 components and 144.9 mL; every component in exactly one category.

## 8. Validation
Synthetic self-tests (assignment order and > 50 % boundaries; anatomy unavailable vs ambiguous; spinal regions with plane ties and UNDEFINED planes; soft-tissue locations incl. head, bifurcation and pelvic limits and missing labels; compartment flag; PSMA component-ID ordering), two-patient preflight, full run, completion note. Script: `02_code\AD_anatomical_distribution_v1_0.py` (Python, env totalseg), MD5-checked inputs, outputs with manifest.
