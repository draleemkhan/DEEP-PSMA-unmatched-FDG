# Dated implementation clarification: PSMA uptake at FDG-lesion locations (SAP §8.2)
**Date:** 12 November 2026, 11:45 (Asia/Kuwait, UTC+03:00)
**Approved by the author:** 12 Nov 2026, 11:40 (U4, U6 approved; U7 not approved; reference-mask safeguard).
**Relates to:** locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d) §8.2. **Implementation detail only.** The locked lesion value (SUVmax of the native PSMA PET), reference definitions (eroded aorta and liver, 1 voxel, PET grid), the three uptake classes and the ≤ 10 mm rule are preserved. No change to the primary analysis, the perturbation results or any other locked item.
**Disclosure:** written after the primary classification and the perturbation results had been seen; **no PSMA uptake value, reference value or uptake class has been computed or viewed.** Related exploratory cross-tabulation: `PAPER3_exploratory_crosstab_note_20261112.md` (MD5 49b96844777d91a8f967256928d638b6).

## 1. Units
- The **217 components classified V_U in the primary analysis** (P3 primary run; reference matching rule). No other component is classified for uptake in this analysis.
- Location: the component's **unperturbed (primary) mapped voxel set on the PSMA PET grid**, as produced by P2 (final transform per R6/Amendment 3a; nearest-neighbour resampling of the FDG component label only). PSMA PET values are **never resampled**.

## 2. Lesion value (locked)
- **SUVmax** = maximum native PSMA PET value over the mapped voxel set. PSMA PET used as supplied (SUV units, no conversion); this is confirmed against the dataset documentation in the preflight and recorded.

## 3. Reference tissues (safeguard: documented before implementation)
For each patient, from the **PSMA CT** TotalSegmentator labels (task = total, official 117-class dictionary as verified in Audit B), in this order:
1. **Select** the label by name: `aorta` (blood pool) and `liver`.
2. **Map to the PSMA PET grid** by physical coordinates (nearest-neighbour; label values only; PSMA CT and PSMA PET share the scanner frame, so no registration is applied).
3. **Erode by 1 voxel on the PET grid**, 3 × 3 × 3 structuring element (26-connectivity, consistent with the conn26 component definition); voxels outside the image are treated as background.
4. **Exclude PSMA tumour-labelled voxels:** remove every voxel of the supplied PSMA TTB label (PET grid) from the eroded mask. Exclusion is applied after erosion, so it does not alter the erosion.
5. **Reference value** = SUVmean of the native PSMA PET over the remaining voxels.
- Recorded per patient: label present (yes/no), voxel counts after mapping, after erosion and after TTB exclusion, and the reference SUVmean.
- **Unavailable reference:** if a label is absent, or no voxel remains after erosion and exclusion, the reference is recorded as **unavailable** for that patient. **No substitute reference** (other organ, uneroded mask, other patient, cohort value) is used. Components in that patient whose class depends on the missing reference are reported as **not classifiable (reference unavailable)**, a category separate from the partial-volume category. (If only the liver is missing, "absent" can still be assigned when SUVmax ≤ blood pool; "low" vs "positive" cannot. If only the aorta is missing, "positive" can still be assigned when SUVmax > liver; "absent" vs "low" cannot. Such components are reported as reference-unavailable, with the partial information listed but not used as a class.)
- No minimum voxel count is introduced beyond "at least one voxel"; the voxel counts are reported so the reader can judge reference robustness.

## 4. Maximum extent (U4, approved)
- **Maximum Euclidean distance (mm) between the voxel centres of the native FDG component** (FDG grid, physical coordinates using the image spacing and direction). A single-voxel component has extent 0.
- Extent ≤ 10 mm → **not classifiable (partial volume)**. This is the geometric definition of the locked rule, not a new eligibility criterion. It is evaluated before the uptake class; such components still have SUVmax recorded and are included in the U6 count.

## 5. Classes (locked)
- **Absent:** SUVmax ≤ blood pool. **Low:** SUVmax > blood pool and ≤ liver. **Positive:** SUVmax > liver. Equality goes to the lower class, as the "≤" in the SAP implies.
- If, in a patient, the blood-pool reference is ≥ the liver reference, the classes are applied exactly as defined (the "low" interval is then empty) and the patient is flagged and reported; no reference is altered.
- Described as **PROMISE-inspired**, not PROMISE V2 scores.

## 6. Exploratory descriptive count (U6, approved)
- Number of the 217 components with **PSMA SUVmax ≥ 3**, with their native FDG volume and number of patients (all 217, including those not classifiable).
- Stated explicitly: **SUV ≥ 3 is the DEEP-PSMA segmentation threshold, not a clinical criterion for PSMA-positive disease.**

## 7. Ligand (U7, not approved)
- Per-patient ligand identity is unavailable (SAP §8.2). Results are reported for the **complete mixed-ligand cohort**; the ligand is **not inferred**. Stated as a limitation (⁶⁸Ga-PSMA-617 and ¹⁸F-DCFPyL differ in biodistribution and blood-pool/liver ratios).

## 8. Reporting
- Per class (absent, low, positive, not classifiable — partial volume, not classifiable — reference unavailable): number of components, native FDG volume (mL) and number of patients (components cluster within patients).
- Per patient: reference availability, voxel counts and reference values (aorta, liver); count of patients with blood pool ≥ liver.
- Per component: extent (mm), SUVmax, class.
- Interpretation follows the locked hierarchy: uptake classes are imaging evidence, not validated biological diagnoses ("most consistent with", never "confirmed").

## 9. Validation (before the full run)
- Self-tests on synthetic volumes: extent (single voxel = 0; known line and cube distances), erosion (known cube shrinks by one voxel per face), TTB exclusion, SUVmax over the mapped set, class boundaries including equality, missing-reference handling, blood pool ≥ liver handling.
- Preflight on two patients: label presence, voxel counts, grid identity (PSMA CT label mapped onto PET grid), and PSMA PET unit confirmation; reviewed before the full run.
- All outputs with MD5 manifest; the script reads the locked P2 primary mapped sets and P3 primary classification without modifying them.
