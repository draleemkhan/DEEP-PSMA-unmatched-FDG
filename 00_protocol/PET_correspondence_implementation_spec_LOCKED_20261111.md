# Paper 3 — Implementation specification: organ displacement (O1) and PET correspondence pipeline (P1–P3)
**Date:** 11 November 2026 (Asia/Kuwait, UTC+03:00)
**Status:** LOCKED (I7). Draft `PET_correspondence_implementation_spec_DRAFT_20261111.md` (MD5 43869af10ee1fe55b8d6007f696fd30c); this version adds the author's three clarifications of 11 Nov 2026 (band alternative mapping, label collisions, volume reconciliation) and the lock record.
**Governing document:** consolidated Step 2 SAP, LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) and the sections it incorporates (R, R.1, τ, soft tissue). This document implements them; it changes no rule. Where they conflict, the SAP governs.

## I0 Principles and invariants (apply to every stage)
1. **Native denominators.** Component volumes, V_FDG,total and every DMTVF denominator come from the **native FDG TTB** (voxel count × native voxel volume). Resampled volumes are used only for matching and are never summed into a denominator.
   - **V_FDG,total** = the **whole** native FDG TTB volume (all components, eligible or not) — the locked denominator (SAP §1). It is identical across connectivities and eligibility rules for a patient (checked).
   - **V_S** = native volume of components failing eligibility; counted **once**, in P1, and never re-entered later. **Eligible volume** = V_FDG,total − V_S = V_I + V_M + V_U, reported separately (with f_small = V_S / V_FDG,total); it is never used as the primary denominator.
2. **Stable IDs.** FDG component IDs are assigned once, in native space (P1), per connectivity, and carried unchanged through P2 and P3.
3. **Reconciliation** (checked at the end of each stage and again at the end of P3; any failure stops):
   - the component set (IDs, native voxel counts, native volumes) is identical across P1, P2 and P3;
   - per patient: Σ native component volumes = native FDG TTB volume (relative difference ≤ 1e-9);
   - per patient: V_S + V_I + V_M + V_U = V_FDG,total, every component has exactly one status, and status counts reconcile with component counts;
   - V_I is split by mutually exclusive cause (registration, boundary band, coverage) and the causes sum to V_I.
4. **Fixed inputs:** every input verified by MD5 (C1) or grid fingerprint before use; original data read-only.
5. **Testing:** each stage has its own synthetic self-test, then a preflight on two patients (TEST codes), then the full run. Preflight outputs are inspected for integrity only (reconciliation, counts, files), not for discordance results.
6. **Results stay closed until the full run completes:** no per-patient DMTVF is displayed during development; the first look is at the complete, reconciled P3 output.

## I1 Stage O1 — organ displacement (Python; CT only; soft-tissue spec SP3)
- **Inputs:** PSMA and FDG CT and totseg_24 (data root); frozen Amendment 2 transforms (as in the locked scripts); Audit 1B label table (MD5 896a1681…); Audit 1C-e structure-level table (for validation).
- **Processing:** six structures (lungs = lobe unions); FDG labels resampled onto the PSMA CT grid by nearest neighbour with the Amendment 2 resampling transform; centroid distance (point-mapped), ASSD and HD95 with the Amendment 1 `surface_metrics` function copied verbatim; eligibility per SP3 (per lobe for lungs).
- **Outputs:** `O1_organ_displacement_<run>.csv` (patient × structure: eligibility and reason, centroid, ASSD, HD95, surface-point counts); manifest.
- **Acceptance:** self-test (synthetic shapes with known shifts, as in Audit 1C-f); composition check; **liver, spleen and kidneys reproduce Audit 1C-e to ≤ 0.01 mm** for every record with identical masks, transform, eligibility and method; all other differences listed with reasons and investigated.

## I2 Stage P1 — native FDG components and eligibility (Python)
- **Inputs:** FDG TTB and FDG PET (native), per patient.
- **Processing:** connected components for 6, 18 and 26 connectivity (separately); per component: ID, voxel count, native volume (mL), native centroid (physical LPS), bounding box, maximum extent (mm), SUVmax and SUVmean from the native FDG PET (descriptors only); eligibility rule A and rule B per connectivity (SAP §2).
- **Outputs:** per patient and connectivity, an integer-labelled native component image (NIfTI, FDG PET grid) and `P1_components_<run>.csv`; patient table with native V_FDG,total; manifest with grid fingerprints.
- **Acceptance:** component counts per patient, tracer FDG and connectivity **identical** to Audit 1D (`Audit1D_component_level_20261030_215034.csv`, MD5 30f5f461…; cohort 4,588 / 3,959 / 3,838); volumes within 1e-6 mL; excluded counts and pooled excluded volume identical to the Option D re-run (rule A 871 / 0.24 %; rule B 1,389 / 0.61 % at 26-connectivity). Synthetic test: known blobs under 6/18/26 adjacency.

## I3 Stage P2 — mapping, registration QC, coverage and resampling (Python)
- **Inputs:** P1 outputs; native PSMA PET and PSMA TTB (headers and grid only for P2); frozen transforms (T_A, adopted head/pelvis transforms; fallbacks); `SAP_R1_plane_table_20261111_214952.csv` (MD5 db6fe214…); failed-region list (SAP R.4 + R.1.4).
- **Per eligible component** (primary connectivity and rule; repeated for sensitivity partitions):
  1. T_A-mapped centroid → R6 region; relevant transforms (SAP R.5 step 2); band flag (within 15 mm of an R6 plane or crossing it).
  2. **Registration QC:** all native voxel centres point-mapped with every relevant transform; any point in a failed region → flag `VI_registration` (functions from `SAP_R1_region_rules.py`, MD5 beedabcf…, imported unchanged).
  3. **Coverage (item 3):** any mapped voxel centre outside the PSMA PET field of view or within τ_i of its boundary (outer grid face) → flag `VI_coverage_fov`.
  4. **Resampling:** each component is resampled onto the **native PSMA PET grid** by nearest neighbour with its **final transform**, **component by component**, and stored as its own voxel set (sparse table: component ID, PSMA-grid voxel indices). Because components may use different transforms (R6), two components can map to the same PSMA voxel; **overlaps are allowed and both identities are kept** — no integer label overwrites another. Matching (P3) always uses each component's own voxel set. An integer label image is also written for visual QC only, with an overlap-count image and a per-component collision count (protocol §14.3). Per component: native and transformed voxel counts, % volume change, zero-survivor flag `VI_coverage_resampling`.
  5. **Band components — alternative transform:** for every band component, steps 3–4 are repeated with the **adjacent region's transform**: its own coverage test (FOV/τ) and its own resampled voxel set and zero-survivor check. If the alternative mapping is **not valid** (coverage fails or zero voxels survive), the alternative category is **"not assessable"**; it is never counted as unmatched (see P3).
- **Outputs:** resampled label image(s) on the PSMA PET grid (NIfTI, int32), `P2_components_<run>.csv` (all flags and counts), grid fingerprints of the PSMA PET grid and of every written image; manifest.
- **Acceptance (self-test):** known synthetic transforms (identity, translation, rotation) reproduce expected mapped positions (≤ 1e-6 mm) and resampled labels; R.1.7 region tests reused; coverage cases (inside, outside, within τ of the boundary, exactly at τ); zero-survivor case (sub-voxel component onto a coarser grid); collision case; fingerprint mismatch stops. Reconciliation (I0.3) with P1.

## I4 Stage P3 — matching, classification and DMTVF (MATLAB R2025a)
- **Start:** `restoredefaultpath`; add project folders only; `which` checks (`knnsearch`, `bwconncomp`, `niftiread`, `prctile`) must resolve to the MathWorks installation; else stop.
- **Inputs:** P2 resampled label images and table; native PSMA TTB (PSMA PET grid); P1 native volumes. Grid fingerprints recomputed and compared with P2; mismatch stops.
- **Matching (τ spec):** τ_i = voxel diagonal of the PSMA PET grid; for each component: d = 0 if any shared voxel, else the minimum boundary-voxel-centre distance (mm) to any PSMA-labelled voxel (boundary = 6-neighbour definition; `knnsearch`); matched if d ≤ τ_i. Correspondence graph to PSMA components (26-connectivity) recorded.
- **Boundary band (SAP R / R6):** for band components not already V_I (registration), the category under the final transform and under the adjacent transform are compared, each category being matched, unmatched or **not assessable** (invalid mapping, P2 step 5). Different categories → `VI_band`. A not-assessable alternative therefore makes the component V_I (band), never V_U. If the **final** mapping itself fails coverage, the component proceeds to V_I (coverage) under the precedence.
- **Precedence:** V_S → V_I (registration) → V_I (band) → V_I (coverage: FOV or resampling) → V_M / V_U; one status per component; cause recorded.
- **Outputs:** `P3_components_<run>.csv` (status, cause, d, matched PSMA IDs), `P3_patients_<run>.csv` (V_S, V_I by cause, V_M, V_U, V_FDG,total, DMTVF_ref, upper bound, f_small); manifest (versions, function paths, fingerprints, MD5s).
- **Acceptance:** τ spec tests T1–T8 re-run **in MATLAB on the identical synthetic inputs** exported by Python, with identical results; synthetic precedence cases (each cause, ties between causes; band component with an invalid alternative mapping → V_I band; two overlapping components from different transforms each matched on their own voxel set); reconciliation (I0.3), including V_FDG,total identical across partitions and V_S counted once.

## I5 Order and gates
O1 (self-test → preflight → full run) → P1 → P2 → P3, each gate passed before the next stage starts. Sensitivity partitions (C6, C18, EB) reuse P1–P3 unchanged with the alternative inputs; τ and S-type sensitivities re-run P3 (and P2 where the failed-region set changes). Perturbation (SAP §8.1), uptake classes, anatomical distribution and voxel-level analyses are separate later stages with their own short implementation notes.

## I6 Preflight patients
TEST001/TEST002 = the same two patients used in all earlier preflights (first two sorted case IDs not in the 5 Nov key). Their preflight outputs are checked for integrity only.

## I7 Lock record
- Approved by the author on 11 Nov 2026 with three clarifications (I0.1 volumes; P2 steps 4–5; P3 band rule).
- **Locked on: 11 November 2026, 23:21 (Asia/Kuwait, UTC+03:00)**, before any code for O1–P3 and before any PET correspondence.
