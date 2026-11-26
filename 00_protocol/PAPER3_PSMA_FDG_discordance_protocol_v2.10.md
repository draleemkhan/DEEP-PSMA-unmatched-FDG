# PAPER 3 — PSMA–FDG discordant tumour burden (DEEP-PSMA paired PET/CT): STUDY PROTOCOL (provisional, conceptually locked)
Version 2.10 — 31 October 2026 (v2.9 + Audit 1C-f pelvic overlay follow-up (train_0069, train_0037): regional rigid misalignment of the pelvis under the adopted transformation documented; adopted transformations unchanged; independence of pelvic QC evidence verified; conditional handling of pelvic FDG volume recorded (Section 8.4.11, follow-up); Audit 1C-g regional QC consolidation executed and accepted (Section 8.5); head/neck to be assessed separately (cervical vertebrae primary evidence; brain as posture warning); Step 3 physics note v2 accepted provisionally (Section 8.6; τ not locked). No new registration, correspondence or discordance analysis.)
Previous versions: v2.9 (31 October 2026), v2.8 (31 October 2026), v2.7 (31 October 2026), v2.6 (30 October 2026), v2.5 (29–30 October 2026), v2.2–v2.4 retained.
Status: conceptual design FINAL. Audits 1A, 1B and 1C completed and closed (technical and anatomical data only). Audit 1D executed (definitive run 20261030_215034, COMPLETE) and Option D descriptive summary done; Audit 1D report pending, audit not yet closed. No tumour correspondence, DMTVF or discordance result has been examined. A provisional framework (Section 16) has been agreed in discussion but is NOT locked. Amendment 2 (Section 15.2) executed once (locked script MD5 ce038cb661670d673bdc7dbbb448f10b; run 20261031_141128, COMPLETE); C1–C4 met; translation-corrected registration ADOPTED as the primary trunk registration (3 fallbacks); head/neck alignment not approved and subject to separate regional QC. Numerical parameters are set after the technical feasibility audits (1A–1D), independently justified, and locked at Step 2 BEFORE any discordance result is examined.
[VERIFY] = check against original source. [LOCK-AFTER-AUDIT] = numerical choice pending audits.

---

## 1. Working title
Robustness and Anatomical Characterization of Automated PSMA–FDG Discordant Tumour Burden Quantification in Metastatic Castration-Resistant Prostate Cancer: A Public Paired PET/CT Study

## 2. Positioning and novelty (after literature audit)
Anatomical distribution of discordance and discordant FDG volume fraction are recognised concepts. Our contribution is methodological: reproducible, automated, uncertainty-aware quantification of FDG-positive tumour burden without corresponding PSMA-labelled disease in an open dataset, with explicit separation of classification uncertainty and anatomical uncertainty.

Key literature (verified at published-abstract / full-text-web level; full PDF and supplement audit outstanding):
| Study | Relevance |
|---|---|
| Telli et al., 2026, Urologic Oncology | 70 patients; 29 (41.4%) discordant; 427 FDG+/PSMA− lesions (bone 236, nodes 95, visceral 88) → anatomical distribution not novel |
| Deng et al., 2026, Frontiers in Oncology (review) | Focal vs widespread discordance; discordant FDG volume fraction already discussed → we do not claim the concept |
| Lucas et al., 2026, EJNMMI | Automated registration-based lesion matching 95.1% (661/695), longitudinal PSMA in BCR → method established; accuracy not transferable to cross-tracer mCRPC |
| Jackson et al., 2025, J Nucl Med | Automated segmentation Dice PSMA 0.94, FDG 0.84, LuPSMA 0.97 → segmentation, not correspondence |
| Buteau et al., 2022, Lancet Oncology (TheraP) | PSMA SUVmean ≥10 predictive; FDG MTV ≥200 mL prognostic → clinical relevance of burden |
| Korean J Radiol 2024 (42 patients); Seifert et al. 2023 (89); Hartrampf et al. 2021 (32) | Site-based discordance, mismatch and eligibility, discordance during therapy |
Final pre-submission search still required (conference proceedings, DEEP-PSMA challenge papers, dataset-specific analyses).

## 3. Data
- DEEP-PSMA training set (Zenodo DOI 10.5281/zenodo.15281784; CC BY-NC 4.0; released 21 Nov 2025, annotations subsequently reviewed/updated): 100 mCRPC patients before 177Lu-PSMA, single centre (Peter MacCallum).
- Per patient: PSMA PET/CT (68Ga-PSMA-11 [VERIFY wording] or 18F-DCFPyL) + FDG PET/CT; EARL; SUVbw.
- Labels: PSMA TTB = SUV ≥3 minus manually removed physiological uptake (threshold.json = 2.999862); FDG TTB = case-specific liver-based threshold (threshold.json; typically 2.5–5).
- Official documentation: "Only voxels above this value should be considered for labelling of disease" → by design, PSMA uptake below SUV 3 can never be PSMA-labelled; label absence ≠ uptake absence (justifies secondary analysis 3).
- Supplied: rigid.tfm in each tracer folder (documented: PSMA folder = "FDG to PSMA"; FDG folder = "PSMA to FDG" — verified numerically in Audit 1C-a, Section 4.1); totseg_24.nii.gz = TotalSegmentator output on each CT, on the CT grid. Official dataset page: "TotalSegmentator (v2.4, task=total)" → "24" denotes version 2.4, not 24 classes. Label dictionary = official TotalSegmentator v2 `total` task, 117 classes (TotalSegmentator README "Class details" / `totalsegmentator/map_to_binary.py`); hard-coded by label ID in all scripts. Fast/full mode not documented. Generating version cannot be proven from the masks themselves → label–anatomy correspondence verified empirically in Audit 1B.
- Official documentation: the rigid registration "may be used to coarsely align PSMA and FDG images"; non-rigid registration was "beyond the scope of our data preparation to validate" for all cases (citable support for our QC-first approach and registration limitation).
- Official template code (GitHub Peter-MacCallum-Cancer-Centre/DEEP-PSMA-Algorithm, inference.py, convert_registration_json_to_transforms()): builds `fdg_fixed_transform` (Euler3D; corresponds to FDG/rigid.tfm) and `psma_fixed_transform = fdg_fixed_transform.GetInverse()` (corresponds to PSMA/rigid.tfm) → the two files are an inverse pair. The example does not resample images with them, so intended usage is not demonstrated.
- Not available: scan dates/interval, tracer per case, scanner per case, clinical outcomes.
- Selected population: all cases have measurable disease on both tracers → prevalence not generalisable.

### Cohort
- All 100 patients initially considered.
- All 100 original DEEP-PSMA cases enter a fresh eligibility audit (Audit 1A) using only `DEEP PSMA EXTRACTED`. No predefined exclusions. No previous scripts, eligibility decisions, exclusions or analysis outputs are reused.
- Eligibility is assessed solely on Paper 3 criteria.
- Coverage/registration problems handled at LESION level (indeterminate), not by excluding patients, unless a predefined case-level technical criterion fails. Incomplete CT coverage is not an automatic exclusion: case-level technical failure is distinguished from regional or lesion-level limitations.
- Status after Audits 1A–1C: all 100 patients retained; no case-level exclusions.

## 4. Core definitions

### 4.1 Coordinate system and transforms
- Reference space: PSMA PET/CT. FDG lesion masks transformed into PSMA space for correspondence.
- ITK resampling convention: resampling FDG onto the PSMA grid requires a PSMA→FDG point mapping; transforming FDG points (centroids, surface points, landmarks) into PSMA space requires FDG→PSMA. A .tfm file is not intrinsically labelled forward/inverse → direction verified numerically in Audit C for BOTH uses.
| Operation | Required mapping |
|---|---|
| Resample FDG onto PSMA grid | PSMA → FDG |
| Transform FDG centroid / surface points into PSMA space | FDG → PSMA |
| Compare anatomical landmarks in PSMA space | FDG → PSMA |
| Sample native FDG image at a PSMA location | PSMA → FDG |

**Working hypothesis (to be tested in Audit C, recorded before any data inspection):**
- `fdg_fixed_transform` = FDG is the fixed image → in ITK convention it maps FDG-space points → PSMA-space points and is used to resample PSMA onto the FDG grid; this matches the FDG/rigid.tfm label "PSMA to FDG" in its RESAMPLING sense.
- Therefore:
| Operation | File to use (hypothesis) |
|---|---|
| Resample FDG mask/PET onto PSMA grid | PSMA/rigid.tfm |
| Map FDG centroids / surface points / landmarks into PSMA space | FDG/rigid.tfm (despite its "PSMA to FDG" label) |
- File labels describe resampling; for point coordinates the "other" file applies. Audit C test: map FDG-CT vertebral (or other skeletal) centroids into PSMA space with each file and check which lands on the corresponding PSMA-CT structures; also check resampled FDG CT vs PSMA CT agreement.

**VERIFIED (Audit 1C-a, 29 Oct 2026; report 06_output\Audit1C_a_report_20261029.md):** hypothesis confirmed in 10/10 pilot cases using SimpleITK's actual conventions (resampling through `sitk.Resample`; point mapping tested separately with `TransformPoint`). The two files are an exact inverse pair in all 100 patients (max error 3.6 × 10⁻¹² mm).
| Operation | File (verified) |
|---|---|
| Resample FDG mask/PET/labels onto PSMA grid (`sitk.Resample`) | PSMA/rigid.tfm |
| Map FDG points into PSMA space (`TransformPoint`) | FDG/rigid.tfm (= PSMA/rigid.tfm inverse) |
| Sample native FDG image at a PSMA location | PSMA/rigid.tfm |
- Native PSMA and FDG world coordinates differ by up to ~2.5 m (different coordinate origins) → the transform is always required.

### 4.2 Lesion units
- Connected components of supplied TTB masks. Primary connectivity [LOCK-AFTER-AUDIT; 26 favoured initially]; 6, 18 and 26 compared in Audit 1D (not chosen by discordance results). Provisional (v2.7, Section 16): 26 primary; 6 and 18 sensitivity.
- Minimum lesion volume: one fixed primary threshold [LOCK-AFTER-AUDIT]; alternatives as supplementary sensitivity. Provisional (v2.7, Section 16): applies to FDG lesions only (assessment eligibility); all PSMA-labelled voxels remain match targets.
- Confluent skeletal disease: rule pending Audit 1D. Skeletal sub-components by individual vertebral level are not used (Audit 1B: apparent vertebral level misassignment in some patients); spinal regions (cervical, thoracic, lumbar, sacral) only.
- Lesion volume and TLG taken from the NATIVE FDG mask and FDG PET (before resampling); transformed mask used only for correspondence.

### 4.3 Lesion correspondence (primary method family)
- Distance-based: minimum surface-to-surface distance in mm, physical coordinates, common grid (FDG transformed into PSMA space).
- One prespecified reference tolerance [LOCK-AFTER-AUDIT]; no unjustified sub-voxel precision (PET voxels ~2.9–4 mm).
- Each FDG lesion assessed against all eligible PSMA components; one-to-many and many-to-one correspondences permitted and recorded (full correspondence graph kept).
- Anatomically implausible matches flagged by predefined anatomical consistency rules (only for structures actually labelled in totseg_24).

### 4.4 Indeterminate lesions (applied FIRST, before any robustness analysis)
Two distinct lesion-level classes (v2.5). They are assigned independently; a lesion may carry both.

| Classification | PSMA–FDG matching | DMTVF contribution | Anatomical site |
|---|---|---|---|
| Site-indeterminate | Assessable | Included according to matching result | Unknown |
| Correspondence-indeterminate | Not assessable | Excluded from numerator; retained in denominator | May be known |

**Site-indeterminate** (anatomy unavailable; correspondence assessable):
- Lesion (or part) outside CT coverage, or insufficient CT/TotalSegmentator information for site assignment.
- Matching uses the PET images, so missing CT coverage alone NEVER makes a lesion correspondence-indeterminate.
- Absence from a TotalSegmentator structure ≠ missing coverage: distinguish unavailable anatomy, segmentation failure, and unlabelled tissue.
- Unlabelled soft tissue is NEVER automatically classified as lymph node; nodal regions are assigned only by predefined anatomical rules (Section 9).

**Correspondence-indeterminate** (PSMA–FDG matching cannot be assessed):
- FDG lesion, after transformation into PSMA space, outside the PSMA PET field of view (criterion [LOCK-AFTER-AUDIT]).
- Inadequate local registration (failed predefined QC). Poor alignment NEVER increases the matching tolerance.
- Lesion disappearing during resampling (Section 14.3).
- Final correspondence-indeterminate status is established only at Step 4D, after transforming FDG lesions into PSMA space and checking coverage; it is never assigned from CT-coverage findings alone.

Rationale: separating the two classes prevents unnecessary underestimation of the primary endpoint.

### 4.5 Terminology
"FDG-positive lesion without a corresponding PSMA-positive lesion under the supplied annotation criteria" — never "PSMA-negative" in the biological sense.

## 5. Endpoints

### Primary endpoint — reference discordant metabolic tumour volume fraction
DMTVF_ref = V_FDG,unmatched / V_FDG,total × 100
- Single prespecified reference configuration.
- Denominator = total FDG-positive tumour volume INCLUDING correspondence-indeterminate lesions; numerator EXCLUDES correspondence-indeterminate lesions (conservative observed fraction).
- Site-indeterminate lesions are classified by their matching result like any other lesion (included in numerator if unmatched).
- Reporting: % patients with DMTVF_ref > 0; median (IQR) in all patients; median (IQR) among patients > 0; continuous values retained (no arbitrary dichotomisation).

### Secondary analyses (priority order)
1. **Indeterminate burden and methodological bounds (essential)**
   - IMTVF = V_correspondence-indeterminate / V_FDG,total × 100; proportion of patients with any correspondence-indeterminate burden.
   - Site-indeterminate burden reported separately (volume, % of FDG burden, number of patients); it does not enter the uncertainty interval.
   - Two separate intervals (anatomical and classification uncertainty are never mixed):
     - **Primary uncertainty interval (anatomical uncertainty):** from DMTVF_ref to (V_unmatched,ref + V_correspondence-indeterminate) / V_FDG,total × 100.
     - **Sensitivity envelope (classification uncertainty):** minimum to maximum DMTVF across the alternative matching configurations (supplementary analyses).
   - Methodological identification intervals, not confidence intervals.
2. **Registration perturbation (principal robustness analysis)**
   - Random multidirectional rigid perturbations (translation + rotation), fixed number per magnitude (provisional 50), locked random seed; magnitudes justified from Audit C residual alignment discrepancy [LOCK-AFTER-AUDIT].
   - Outputs: change in DMTVF_ref, lesion classification changes, patient-level discordance status changes.
   - Limitation: whole-image rigid perturbation does not simulate local respiratory deformation (liver dome, lung bases).
3. **PSMA uptake classes at unmatched FDG lesion locations**
   - Question A: does the FDG lesion correspond to a PSMA-labelled lesion? Question B: how much PSMA uptake is present at the corresponding location?
   - PSMA uptake sampled from the NATIVE PSMA PET at the transformed FDG lesion location (FDG lesion mapped into PSMA space); PSMA uptake values are never resampled.
   - Classes: absent / low / positive, relative to a prespecified reference tissue (e.g., liver or blood pool from TotalSegmentator on the PSMA CT) with prespecified handling of small lesions (partial volume) [LOCK-AFTER-AUDIT]. Thresholds not derived from observed separations. Quantitative classes, not visual interpretations.
4. **Anatomical distribution** of unmatched FDG burden: bone (spine reported by cervical, thoracic, lumbar and sacral regions, not individual vertebral levels; other skeletal sub-regions if labels allow), lymph nodes (nodal regions defined by predefined anatomical rules using labelled reference structures, e.g., aorta, IVC, iliac vessels, pelvic bones, vertebrae — unlabelled soft tissue is not automatically nodal), visceral organs; unassignable soft tissue reported as its own category; site-indeterminate burden reported as its own category; account for site-specific disease amount and within-patient clustering.
5. **Supplementary:** connectivity (6 vs 26, and 18 if relevant after Audit 1D); distance-tolerance sensitivity (nested thresholds reported as sensitivity, not robustness); alternative matching method (overlap-based) agreement (% consistent, kappa, burden difference, patient status change — agreement ≠ accuracy); anatomical consistency of matched pairs; patient-level PSMA vs FDG burden relationships (MTV, SUVmean, TL-PSMA, TLG).

### Exploratory robustness measure
- Perturbation-based stability: % of perturbations in which an FDG lesion remains unmatched at the reference tolerance; provisional stability criterion ≥95% [LOCK-AFTER-AUDIT].
- Rationale: nested distance tolerances make "unmatched across all tolerances" equal to the permissive-tolerance DMTVF (mathematical redundancy) → not used as a robustness definition.
- A connectivity-based robust fraction considered only if Audit 1D shows material partition effects; not manufactured otherwise.

### Uncertainty types (reported separately)
- Classification uncertainty: plausible analytic choices change correspondence.
- Anatomical uncertainty: registration or PET field-of-view coverage prevents assessment of correspondence (correspondence-indeterminate). Missing anatomical site information (site-indeterminate) is reported descriptively and is not part of this uncertainty.

## 6. Non-negotiable constraints / claims NOT made
- No expert adjudication.
- No claim of diagnostic accuracy of automated correspondence.
- No biological confirmation of PSMA-negative disease.
- No prediction of LuPSMA outcomes; DMTVF is a study-defined imaging measure, not a validated biomarker.
- No generalisable prevalence for unselected mCRPC.
- No registration accuracy claim against independent ground truth (TotalSegmentator-based QC = residual anatomical alignment discrepancy, which includes respiration, positioning and segmentation differences).

## 7. Registration QC framework (hierarchical)
| Level | Measurement | Purpose |
|---|---|---|
| Global | Skeletal landmark displacement | Overall rigid alignment |
| Regional | Vertebral / pelvic landmark displacement | Local alignment |
| Organ | Surface distances (liver incl. dome, kidneys, lungs, pelvic bones) | Regional discrepancies |
| Lesion neighbourhood | Local anatomical consistency | Plausibility of correspondence |
Decision rule: acceptable → eligible; borderline → predefined sensitivity analysis; inadequate → indeterminate. Limits, including the definition of the borderline category and its sensitivity analysis, [LOCK-AFTER-AUDIT], justified independently of discordance findings.
- No patient is classified as a registration failure or excluded before the QC criteria are locked at Step 2. Audit 1C-b reports regional alignment distributions (centroid and surface distances), anatomical coverage and the proportion of patients exceeding candidate limits (descriptive only); final limits are never selected from discordance outcomes.
- Registration refinement: Amendment 1 (MI rigid optimisation) NOT adopted (C3 and C4 not met). Amendment 2 (translation-only, Section 15.2) ADOPTED (31 Oct 2026): the corrected transformation is the primary trunk registration for 97 patients; the supplied transformation is retained for train_0018, train_0049 and train_0067 (per-patient fallback) and as a sensitivity analysis for all. Head/neck alignment is not approved by this decision (Section 15.2.11).
- A poorly aligned patient does not make every lesion indeterminate: regional (and, where feasible, local) QC determines which lesions are assessable. Final QC limits must be independently justified and locked at Step 2, before discordance is examined.
- Landmark/structure eligibility (Audit 1B): only structures apparently complete according to automated coverage criteria on BOTH CTs are scored. Femora (paired-complete in 10/100 patients) and skull (37/100) are not used; humeri, clavicles, scapulae and brain only where paired-complete. Vertebrae with apparent level misassignment or substantial fragmentation, and adjacent levels, are flagged as unreliable individual landmarks (final exclusion rule [LOCK-AFTER-AUDIT] at Step 2).
- Coverage safeguard: before computing any alignment metric, establish whether each structure is fully imaged on BOTH CTs (axial field of view, truncation at scan edges). Surface distances computed only for structures fully covered on both scans (or restricted to the common field of view); partially imaged organs flagged, not scored, so coverage differences are not mistaken for misregistration.

## 8. Data audits (sequential; each documented before the next)
| Audit | Task | Key questions |
|---|---|---|
| A (COMPLETED, CLOSED) | Cohort eligibility (fresh raw-data audit, all 100 cases) | (1) presence of all 12 expected files per case (per tracer: CT, PET, TTB, threshold.json, totseg_24, rigid.tfm; JPG previews / Thumbs.db ignored); (2) file readability and integrity; (3) PET, CT and TTB dimensions and voxel spacing; (4) within-tracer geometry (TTB on PET grid, totseg_24 on CT grid, PET–CT physical overlap); (5) PSMA and FDG threshold values; (6) non-empty PSMA and FDG tumour masks; (7) PET vs CT field-of-view coverage; (8) missing, corrupted or geometrically inconsistent data. Outputs: case-level QC table (100 cases), summary of technical abnormalities, eligibility assessment on Paper 3 criteria only, written Audit 1A report and reproducibility manifest |
| B (COMPLETED, CLOSED) | TotalSegmentator labels and anatomical coverage (all 100 cases, 200 CTs; MATLAB; no spatial transformation) | (1) labels vs official 117-class dictionary (v2.4, task=total), labels outside 1–117; (2) independent verification of the mapping: presence, volumes (robust outliers), centroid positions, laterality, vertebral ordering (sacrum→C1), craniocaudal organ order, fragmentation of single-piece structures; (3) per structure: present-complete / present touching a CT boundary (z, in-plane; suggests possible truncation, does not prove it) / absent with expected position outside CT, at CT boundary, unexplained (segmentation failure, anatomical absence e.g. post-surgical, or other — not interpreted clinically) or unknown; (4) completeness recorded separately for PSMA CT and FDG CT; paired candidates (complete on both) for registration QC — final suitability criteria locked at Step 2; (5) availability of reference structures for nodal-region rules and the PSMA uptake reference tissue. Outputs: label-, scan-, paired- and cohort-level tables, abnormalities, written Audit 1B report, manifest (dictionary version and source) |
| C (COMPLETED, CLOSED) | Transforms and registration QC (Python/SimpleITK) | 1C-a (COMPLETED): ten-case transform-direction and software-verification pilot (Section 4.1 verified). 1C-b (COMPLETED): all 100 patients with the verified convention — residual anatomical alignment by region (head/neck, thorax, upper abdomen, pelvis; centroid and surface distances), coverage measures (Section 8.1), coverage of each structure on both CTs, proportion of patients exceeding candidate limits (descriptive). 1C-c (COMPLETED; Amendment 1): bone-based rigid refinement, NOT adopted. No tumour correspondence, no DMTVF, no discordance results |
| D (EXECUTED — definitive run 20261030_215034 COMPLETE; report pending; not yet closed) | Lesion masks (MATLAB; native PET space; TTB masks unchanged) | (1) connected components with 6-, 18- and 26-connectivity: counts, merges, % TTB volume in merged components; (2) lesion volumes and per-component SUVmax/SUVmean from the native PET of the same tracer; (3) exploratory size thresholds 0.1/0.25/0.5/1/2 mL and 2/5/10 voxels, PSMA and FDG separately (descriptive; not eligibility criteria); (4) extent of confluent disease with physical-coordinate anatomical lookup in the same scan (spanning several labelled structures = automated observation, not proof of a single metastasis); (5) grid geometry and FDG/PSMA voxel-size compatibility. No intertracer registration, correspondence, overlap or discordance |
Original DEEP-PSMA files are never modified. Each audit is completed and documented before the next begins.

### 8.1 Revised coverage measures (v2.5; replace the Audit 1A PET-field criterion)
The Audit 1A criterion "CT covers ≥95% of the PET image extent on each axis" is retired: it measured the reconstructed PET field (largely air), not the patient. Patients are never excluded because reconstructed PET and CT fields of view differ. Coverage is instead measured as follows.

**Audit 1C — case level (per patient, per tracer; descriptive, no discordance outcomes examined)**
| Measure | Purpose |
|---|---|
| % of TTB volume inside the CT volume (in-plane vs z) | Tumour available for anatomical site assignment |
| % of PET body outline inside the CT volume (in-plane vs z); body extent beyond CT superiorly/inferiorly (mm) | Anatomically relevant CT coverage (body, not empty PET field) |
| CT in-plane edge truncation (slices with body at the image edge, by edge/direction) | Truncated anatomy (e.g., arms, lateral body wall) |
| % of FDG TTB volume (and FDG body) inside the PSMA PET field of view after FDG→PSMA point mapping (verified transform, Section 4.1); common z-range of both PET scans | Expected correspondence assessability |
| Coverage of each TotalSegmentator structure used for registration QC on both CTs | Registration QC only on fully imaged structures (Section 7) |

**Step 4D — lesion level (per FDG lesion)**
| Measure | Classification |
|---|---|
| Fraction of lesion volume inside the CT on which its site is assigned | Below threshold → site-indeterminate |
| Fraction of the transformed lesion inside the PSMA PET field of view, with a margin equal to the reference matching tolerance (a PSMA match must be detectable within the tolerance) | Below threshold → correspondence-indeterminate |
Thresholds [LOCK-AFTER-AUDIT], set at Step 2 from Audit 1C, independently of discordance findings. All measures in physical (mm) space using the image affines.

### 8.2 Audit 1D execution and descriptive findings (v2.7; single tracer; no intertracer analysis)
**Execution.** Script `02_code\Audit1D_lesion_mask_audit.m` (MATLAB; protocol v2.6 scope). Execution test on 2 patients (run 20261030_214350, status TEST_ONLY; not an audit result). Definitive run 20261030_215034: status COMPLETE — 100/100 patients, 200/200 tracer examinations OK, 0 PARTIAL, 0 FAILED, all mandatory outputs and verified manifest written; elapsed 25.5 min. Script MD5 of the definitive run recorded in its manifest. Outputs in 06_output.

**Grid geometry.** All 200 PET grids non-oblique, identical orientation; TTB on PET grid and totseg_24 on CT grid in all scans. Voxel volume 0.0076–0.098 mL (~13-fold); slice thickness 1.5–5 mm. PSMA and FDG voxel volumes differ in 21 patients (FDG/PSMA ratio 0.32–2.75). A fixed mL threshold therefore corresponds to different voxel counts (FDG: 0.1 mL = 1.0–4.8 voxels; 0.25 mL = 2.6–12.0 voxels).

**Connectivity (components, cohort total).** PSMA 7,902 / 6,999 / 6,841 and FDG 4,588 / 3,959 / 3,838 for 6 / 18 / 26-connectivity. Median % of TTB volume in merged components (6→26): PSMA 41.6%, FDG 15.1% (maximum ~98% both). 18→26: median ~0% but large in a minority of scans (whole merged-component volume counted; measures partition sensitivity, not segmentation error).

**Small FDG components (26-connectivity; 3,838 components; 452 single-voxel).**
| Below | Components | Pooled FDG volume | Per-patient volume, P90 / max | Median voxels | Median SUVmax |
|---|---|---|---|---|---|
| 0.1 mL | 20.2% | 0.18% | 1.3% / 27.8% | 1 | 3.2 |
| 0.25 mL | 36.2% | 0.61% | 5.9% / 50.0% | 2 | 3.3 |
| 0.5 mL | 50.1% | 1.41% | 12.1% / 70.0% | 3 | 3.4 |
| 1.0 mL | 62.5% | 2.86% | 22.6% / 100% | 5 | 3.6 |
- Combined rule < 0.1 mL OR < 3 voxels: 871 components, 0.24% pooled volume, per-patient P90 1.7%. At 0.25 mL a 3-voxel criterion adds nothing; a 5-voxel criterion adds 56 components (0.06%).
- Per-patient maxima arise in patients with very small total FDG TTB (e.g., train_0006 4.5 mL, train_0053 5.0 mL; 5 patients < 5 mL).
- Small components lie close to the labelling threshold (median SUVmax ~3.3 vs 6.4 for components ≥ 1 mL); descriptive only — per-patient FDG thresholds not compared.
- Present in 93 patients (< 0.25 mL), but the 10 highest-burden patients hold 36%.
- Largest-component share (f_largest, 26-connectivity): FDG median 0.43; > 0.5 in 36 patients, > 0.8 in 11. PSMA median 0.42; > 0.5 in 36, > 0.8 in 7.
- Large components examined anatomically: 1,234 (PSMA 839, FDG 395); spanning ≥ 2 labelled groups: 354 (automated observation, not proof of confluent metastasis).
- PSMA small components (context only; all retained as match targets): < 0.25 mL = 23.0% of components, 0.23% of pooled volume.

**Option D descriptive summary.** Python script `02_code\Audit1D_optionD_small_component_summary.py` (reads the two definitive Audit 1D CSVs; descriptive only; writes Audit1D_optionD_* outputs and manifest with input MD5s to 06_output).

### 8.3 Audit 1C-d — descriptive regional registration-QC review (v2.8; supplied transformation; no matching)
Script `02_code\Audit1C_d_regional_QC_review.py` (Python; re-analysis of the Audit 1C-b tables, run 20261029_194714; descriptive only; outputs Audit1C_d_* with manifest in 06_output). Groupings below are descriptive, not QC limits.
- Whole-body bony median centroid distance: ≤ 10 mm in 65 patients, 10–15 mm in 14, > 15 mm in 21.
- The residual is largely a per-patient rigid translation: in the > 15 mm group the median displacement vector is 0.98 of the median distance (≈ 20 mm in-plane, ≈ 1 mm craniocaudal); 10–15 mm group 0.93; ≤ 10 mm group 0.71.
- After removing each patient's median translation (in-sample, descriptive): residual median 3.8 mm; median per-patient P90 8.0 mm; 10–15 mm group P90 13.9 mm. Not used to select a tolerance.
- Trunk regions move together: Spearman ρ 0.85–0.90 between thorax, upper-abdomen and pelvis bony medians; within-patient SD across these regions median 1.1 mm (P90 3.0). Head/neck less concordant (ρ 0.59–0.74); excess over thorax P90 8.6 mm.
- Soft tissue minus bone (regional medians): thorax 0.5, upper abdomen 1.9, pelvis 2.4, head/neck 4.4 mm. Least stable organs: spleen, bladder, liver (median 10–12 mm; P90 ≈ 22–23 mm); most stable: aorta, heart.
- Implication (descriptive): with the supplied transformation, trunk regional QC acts in practice as a patient-level decision for the high-offset patients; this motivated drafting Amendment 2 (Section 15.2).

### 8.4 Audit 1C-f — Independent vascular and muscular registration QC (DRAFT specification, 31 October 2026)
**Status: specification and script APPROVED and LOCKED (author, 31 Oct 2026; Section 8.4.10). Single definitive run authorised once the script MD5 is recorded. No residual has been computed for any structure in this set before the definitive run (pre-flight reported implementation checks only).** Execution sequence (author): specification → review → script → synthetic validation → MD5 lock → single definitive run.

**8.4.1 Purpose and scope.** Descriptive, regional measurement of PSMA–FDG CT alignment on anatomical structures that were not used to estimate any registration (Amendments 1 and 2) and that are relevant to nodal and soft-tissue disease. Results are reported separately from the Amendment 2 held-out evaluation, as additional anatomical evidence; they are not proof of lesion-level registration accuracy. Outputs feed the later regional residual e_r (Section 16.5). **Not permitted in Audit 1C-f:** acceptance thresholds or categories, tumour masks (of any kind), lesion matching, correspondence, DMTVF or discordance.

**8.4.2 Inputs (read-only).** totseg_24 on both native CTs; supplied rigid.tfm files; adopted transformation per patient (translation t and fallback status from 06_output\Audit1C_e_per_patient_20261031_141128.csv); Audit 1B label-level table (06_output\Audit1B_label_level_20261029_134401.csv) for completeness and fragmentation; Audit 1B vertebral flags.

**8.4.3 Transformations evaluated (no new optimisation).**
- T_S: supplied (resampling PSMA/rigid.tfm; all 100 patients).
- T_A: adopted (translation-corrected for 97 patients; supplied for train_0018, train_0049, train_0067).
- T_C (descriptive only, reported separately): translation-corrected transformation for the 3 fallback patients (train_0018, train_0049, train_0067). T_A (97 + 3) remains the primary configuration.
- FDG label maps resampled onto the PSMA CT grid (nearest neighbour), as in Amendments 1 and 2.

**8.4.4 Structure set (fixed before any measurement).**
| Class | Labels (TotalSegmentator v2.4, task=total) | Inclusion rule |
|---|---|---|
| Vascular | aorta (52), inferior vena cava (63), iliac vena left/right (67/68), iliac artery left/right (65/66) | PRESENT_COMPLETE on both CTs (Audit 1B) AND frac_largest_component ≥ 0.95 on both CTs (fragmentation ≤ 5%) — applied identically to all vessels |
| Core muscular | iliopsoas left/right (88/89), gluteus medius left/right (82/83), gluteus minimus left/right (84/85) | PRESENT_COMPLETE on both CTs |
| Conditional muscular | gluteus maximus left/right (80/81) | PRESENT_COMPLETE on both CTs |
| Exploratory (reported separately; excluded from regional summaries) | prostate (22) | PRESENT_COMPLETE on both CTs |
- Expected effect of the fragmentation rule (Audit 1B): IVC excluded in up to 11 patients, iliac veins in 3/8, iliac arteries in 34/38 (left/right); exact numbers reported.
- Every exclusion is recorded with its reason; excluded structures are never replaced by others.

**8.4.5 Regions and boundaries (defined on the PSMA CT only, the reference space).**
- Boundaries are axial planes at the midpoint between the centroids (z) of adjacent vertebrae on the PSMA CT: B_thorax/abdomen = midpoint T12–L1; B_abdomen/pelvis = midpoint L5–S1; upper limit of the thorax = midpoint C7–T1.
- Regions: thorax (between C7–T1 and T12–L1), retroperitoneum/abdomen (between T12–L1 and L5–S1), pelvis (below L5–S1). Head/neck is not covered (no structure of this set lies there); it remains subject to separate head/neck QC.
- Landmark substitution: if a boundary-defining vertebra is not PRESENT_COMPLETE on the PSMA CT or is flagged unreliable (Audit 1B: train_0088 T2–T4 ±1 level; train_0059 L3–L5 ±1; train_0029 T12–L2 ±1, as recorded in Section 12), the boundary is estimated by linear interpolation of the vertebral-level positions from the nearest verified vertebrae above and below (e.g., T12–L1 from T11 and L2 centroids). Each substitution is written to a per-patient substitution log (patient, boundary, vertebrae used).
- Long structures are evaluated per region: each structure is split at the region boundaries, giving structure-segments (e.g., aorta → thoracic and abdominal segments; iliopsoas → abdominal and pelvic segments).

**8.4.6 Metrics (fixed before any measurement).**
- Grid: PSMA CT grid (FDG labels resampled with each transformation); all coordinates in physical (LPS, mm) space.
- **No cropping of masks.** Surface points: BinaryContour (fullyConnected = False) of the FULL PSMA mask and the FULL resampled FDG mask, converted to physical voxel-centre coordinates. No artificial cut surfaces are created.
- **Windows, applied identically in both directions** (author safeguard against nearest surfaces from outside the intended region):
  - Region window W_r: the z-interval common to both masks within region r.
  - Sampling window S = W_r shrunk by the end margin m = 10 mm at each end (provisional; tested synthetically).
  - Target window T = W_r (the full regional window). For each direction (PSMA → FDG and FDG → PSMA), distances are computed from the source surface points inside S to the nearest target surface points of the other mask inside T only (exact nearest-neighbour search on physical coordinates). Surface outside the region can therefore never supply a nearest point; restricting the target set can only increase, never decrease, a distance, and the margin m keeps sampled points away from the target window ends.
- **Evaluability (provisional; tested synthetically):** S ≥ 20 mm long AND ≥ 100 sampled surface points in EACH direction; otherwise the segment is "not evaluable" (reason and point counts recorded). The 20 mm length is not interpreted as a validated minimum for reliable vessel QC; the synthetic tests report MSD stability against window length and point count.
- **Primary metric — symmetric MSD:** mean of the pooled distances of both directions (Amendment 1's pooled definition, restricted as above).
- Supporting: pooled HD95 within the windows; for vessels, slice-wise in-plane centroid displacement (mean over S of per-slice cross-section centroid differences, x and y separately, LPS mm); for muscles, centroid distance of the S-restricted masks.
- Distance definition note: voxel-centre-to-voxel-centre surface distances (not Amendment 1's distance-map values); identical masks give exactly 0.
- Known limitation reported: vessels (approximately prismatic tubes) are insensitive to displacement along their own axis; they measure displacement perpendicular to their course. Muscles deform with posture; their residual includes posture-related deformation.

**8.4.7 Synthetic validation (built into the script; runs before any patient data; any failure stops the run).**
- Shapes: straight cylinder (vessel-like), cylinder tilted 30°, ellipsoid and tapered block (muscle-like), on an anisotropic PET-CT-like grid and an isotropic grid.
- Cases: zero displacement; known translations along x, y and z (e.g., 2, 5, 10 mm); rotations (e.g., 5° and 10° about each axis), applied through the same SimpleITK resampling path used for patients.
- Pass criteria (fixed in advance): zero displacement → MSD = 0 (identical masks); slice-wise in-plane centroid displacement recovers in-plane translations within 0.5 × the largest voxel spacing; MSD increases monotonically with translation magnitude and rotation angle; z-translation of the straight cylinder gives MSD ≈ 0 within the window (axial insensitivity demonstrated) while the ellipsoid detects it; results on anisotropic and isotropic grids agree within 0.5 × the largest voxel spacing.
- Cut-face check: for in-plane shifted cylinders, the windowed full-mask MSD must equal the analytically/numerically expected lateral-surface value within tolerance and must NOT be lower than the expected value; a naive "crop both masks, then MSD" comparator is computed only to demonstrate the bias it would introduce.
- Region-restriction check: a synthetic structure crossing a region boundary (e.g., a curved tube) with a known in-plane offset — the restricted measurement must not be lower than the expected value (no favourable nearest points from outside the region).
- Stability check: MSD for in-plane offsets reported against window length (10–100 mm) and sampled point count, to document the provisional 20 mm / 100-point evaluability rule.
- Qualification: synthetic tests validate the measurement implementation, not the anatomical accuracy of the registration.

**8.4.8 Outputs (06_output, single folder; run timestamp; descriptive only).**
- Structure-segment × transformation table: MSD, HD95, supporting metrics, window length (mm), number of sampled surface points, evaluability and exclusion reasons.
- Patient × region × class (vascular, muscular) table: number of evaluable segments, median, P90, maximum and second-largest MSD, per transformation; paired differences (T_A − T_S).
- Cohort summaries by region and class; exploratory prostate table (separate).
- Landmark substitution log; synthetic validation report; input inventory with checksums; manifest (parameters, software versions, script MD5, completion status as for Audits 1D and 1C-e).

**8.4.9 Author decisions (31 October 2026).** (1) Fragmentation ≤ 5% applied to all vessels, including IVC and iliac veins. (2) T_C evaluated descriptively for the 3 fallback patients, reported separately; T_A (97 + 3) primary. (3) End margin 10 mm accepted provisionally, subject to synthetic testing and a minimum number of genuine surface points; 20 mm is not a validated minimum vessel length. (4) Vertebral-midpoint boundaries approved, with predefined interpolation for unreliable levels and a complete substitution log. Additional safeguard: both directional measurements restricted consistently to the regional windows (8.4.6). Section 8.4 remains DRAFT until the script and synthetic-validation criteria are reviewed.

**8.4.10 Implementation check and lock record.**
- Script: `02_code\Audit1C_f_vascular_muscular_QC.py` (Python/SimpleITK/scipy). Modes: `--synthetic-only` (no patient data), `--preflight N` (implementation checks only; no residual values written or printed), definitive (refused until `SPEC_APPROVED = True`).
- Synthetic validation (author's PC, real SimpleITK; identical values to the developer test): ALL PASSED — zero displacement MSD = 0 for all shapes; in-plane translations 2/5/10 mm recovered (isotropic grid exactly 2.00/5.00/10.00 mm; anisotropic within 0.1 mm); MSD within tolerance of the analytical cylinder value and monotonic; straight-cylinder z-shift MSD = 0 (axial insensitivity), ellipsoid detects z-shift; rotations 5°/10° about x, y, z detected and monotonic; cut-face check: naive crop-then-measure biases MSD downward by 0.32–0.46 mm, the specified method does not; region restriction never lower than unrestricted, and in a 20 mm stress case the unrestricted search returns falsely smaller distances while the restricted one does not; grid agreement within tolerance. Disclosure: the first stress-case geometry (15 mm) did not bring out-of-region surface within reach on the 3 mm grid and therefore tested nothing; the test geometry (not the method or criteria) was corrected to 20 mm before the author's run.
- Pre-flight (run 20261031_164034; train_0001–0003 + fallback train_0018; script MD5 080d8dc11c271ca5250f96c06438aefe): 4/4 OK; Audit 1C-e aorta centroid distances reproduced with 0.0 mm difference (supplied and corrected); composition error 0.0 mm; all boundaries DIRECT and correctly ordered; eligibility 54/60 (exclusions exactly per rules: fragmented iliac arteries, gluteus maximus at the field-of-view boundary); T_C path exercised.
- Known consequence of the no-extrapolation rule (author decision: accepted as specified): in train_0059 (S1–L2 flagged unreliable) the L5–S1 boundary is undefined, so its pelvis and abdomen/retroperitoneum regions are not evaluated in Audit 1C-f (logged). Boundaries of train_0029 and train_0088 are expected to be interpolated (SUBSTITUTED, logged).
- Final code change: `SPEC_APPROVED = True`. **Locked script MD5: e893c4ff12dd871c52d2c992f8c2b7b3** (40,868 bytes; author's Get-FileHash on the copy in 02_code, 31 October 2026, matching the delivered file). No code changes after the MD5 is recorded.
- **LOCK COMPLETE (31 October 2026).** The single definitive Audit 1C-f run is AUTHORISED (no flag). Descriptive only; no thresholds, tumour data, matching or discordance.


**8.4.11 OUTCOME and decisions (run 20261031_164711; status COMPLETE; full report: 06_output\Audit1C_f_report_20261031.md).**
- Integrity: 100/100 processed; locked MD5 confirmed; synthetic validation passed; Audit 1C-e aorta centroid distances reproduced with 0.0 mm difference in all 100 (supplied and corrected); all evaluable segments ≥ 250 surface points per direction.
- Eligibility: 128 exclusions per rules (fragmentation 95: iliac arteries 72, IVC 11, iliac veins 11, aorta 1; not paired-complete 33: gluteus maximus 30, iliopsoas 1, prostate 2).
- Regional evaluability: 9 boundary substitutions (two across long gaps: train_0078 C7–T1 from T6/C6, train_0088 from T9/C7); 4 undefined boundaries → train_0007 no region (L5 and L2 fragmented; neighbour exclusion removes S1–L1), train_0059 and train_0091 thorax only. Limitations of the conservative boundary rules, not registration failures.
- Results (median of patient-level median MSD, T_S → T_A; P90 in brackets): thorax vascular 3.19 → 2.36 (9.92 → 3.76); abdomen vascular 3.29 → 1.75 (11.96 → 2.82), muscular 2.85 → 1.59 (9.49 → 3.17); pelvis vascular 2.57 → 2.09 (11.89 → 5.58), muscular 2.91 → 2.44 (8.60 → 4.89). Vessel cross-section in-plane displacement median 6.10 → 3.85 mm, no systematic direction.
- Deterioration > 2 mm (T_A − T_S): pelvis vascular 3 patients (train_0069 +10.1, train_0037 +7.1, train_0058 +2.3), pelvis muscular 4 (train_0069, train_0037, train_0064, train_0068), abdomen 1 (train_0069); thorax none. Pattern consistent with non-rigid differences (high Audit 1C-b vertebral gradient); not visible in Amendment 2's held-out set (no pelvic structure).
- Fallback patients (T_C, descriptive): corrected transformation would improve abdomen/pelvis but worsen thorax in train_0067 (+2.0) and train_0018 (+0.7).
- Prostate (exploratory): not evaluable in any patient (window rules require ≥ 40 mm common extent).
- Evidence density: median evaluable segments per patient — thorax vascular 2, abdomen vascular 2 / muscular 2, pelvis vascular 3 / muscular 8; only 107/487 patient-region-class combinations have ≥ 5; abdomen has no held-out bone.
- **Author decisions (31 Oct 2026):** Audit 1C-f accepted as technically complete with limitations documented; no supplementary alternative-landmark analysis now (3 patients, causes understood); Amendment 2 fallbacks preserved (T_C pattern reported only; reconsideration would need a separate documented decision); prostate non-evaluability documented, no further analysis; no five-structure requirement that most regions cannot satisfy — independent QC evidence (Amendment 2 held-out bone/soft tissue + Audit 1C-f vascular/muscular) to be combined by region before acceptance criteria are set; next: pelvic overlays of train_0069 and train_0037 under both transformations, then the regional acceptance framework. No tumour matching or discordance authorised.
**8.4.11 Follow-up — pelvic overlays under supplied and adopted transformations (v2.10, 31 October 2026; visual, descriptive).**
- Purpose: to distinguish residual registration error from anatomical deformation in the two patients with the largest pelvic deterioration under the adopted transformation (T_A). Visual only: no metrics, thresholds, tumour data, matching or discordance.
- Script: `02_code\Audit1C_f_pelvic_overlays.py`, MD5 a31cfc4f32538d03752ae940f4feefa8 [VERIFY against the script_MD5 field of `Audit1C_f_pelvic_overlay_manifest_<run>.json`]. PSMA CT with TotalSegmentator contours (PSMA labels green; FDG labels resampled magenta); coronal, sagittal and axial pelvic views; supplied (T_S) vs adopted (T_A); composition of the corrected transform re-verified (< 10⁻⁶ mm) before rendering. Run ID: [RUN ID from manifest].
- Patients with pelvic deterioration > 2 mm under T_A (Audit 1C-f, Section 8.4.11): train_0069, train_0037, train_0058, train_0064, train_0068. Overlays inspected for train_0069 and train_0037 only; the other three (smaller deterioration, ≈ 2.0–2.3 mm) were not inspected visually.
- Observations (visual estimate, not measured):
  - train_0037 (t = (−17.6, 0.7, 3.5) mm, |t| 18.0 mm): under T_S the hip bones, sacrum and femoral cortices are closely congruent. Under T_A all pelvic structures (bones, gluteal muscles, iliopsoas, iliac vessels) are displaced together towards the patient's left, by an amount comparable to |t_x|.
  - train_0069 (t = (−15.0, 7.4, −4.8) mm, |t| 17.4 mm): same pattern — near-congruent bone under T_S; under T_A a uniform displacement of the pelvic structures towards the patient's left and slightly anteriorly, with clear separation of the iliac vessels.
  - In both patients the pelvic structures move together under T_A: no visible tilt, rotation or structure-specific deformation within the pelvis, and no vertebral-level slip or gross rotational error. Bladder and prostate differ under both transformations (filling/physiological; not used as registration evidence).
- **Interpretation:** regional rigid misalignment of the pelvis under the adopted transformation. In these patients the supplied transformation was locally better aligned in the pelvis, whereas the translation t — estimated mainly from thoracolumbar vertebrae — improves the upper trunk. A difference in patient positioning between pelvis and trunk across the two examinations is a possible explanation; it is not proven. This refines the earlier wording "consistent with non-rigid differences" (Section 8.4.11 above): the difference is non-rigid between regions but appears rigid within the pelvis in these two patients.
- **Independence of pelvic QC evidence (verified against Amendment 1 (Section 15), Section 15.2.2 and the Audit 1C-f script):**
  - Amendment 2 fitting structures (in-sample; cannot independently support pelvic acceptance): vertebrae T1–L5, S1, sacrum, left and right hip.
  - Femora: not fitting structures, but excluded from all registration evaluations (paired-complete in 10/100, Audit 1B; Section 7); therefore no independent evidence either.
  - Amendment 2 held-out bone: ribs, sternum, clavicles, scapulae — no pelvic structure.
  - Amendment 2 held-out soft tissue in the pelvis: urinary bladder only (filling-dependent; limited value as alignment evidence).
  - Audit 1C-f (independent): iliac arteries and veins (65–68; iliac arteries frequently excluded for fragmentation) and gluteus maximus/medius/minimus and iliopsoas (80–85, 88, 89).
  - Hence independent pelvic evidence is vascular and muscular only; in-sample bony congruence (hips, sacrum) must not be counted towards regional acceptance.
- **Author decisions (31 Oct 2026):**
  - Adopted transformations unchanged; no further registration correction is introduced on the basis of these cases.
  - train_0069 and train_0037 documented as examples of regional deterioration under Amendment 2; all five patients with pelvic deterioration listed.
  - These findings support regional alignment-acceptance criteria rather than assuming that a globally improved transformation is satisfactory in every region.
  - Only genuinely held-out structures may independently support regional acceptance.
  - Conditional rule (in principle; not applied): if a region's alignment fails the subsequently locked QC criteria, the affected FDG tumour volume contributes to V_I, not to V_U. No patient or region is classified as indeterminate now, because acceptance criteria are not yet defined.
  - No new registration, correspondence or discordance analysis authorised.

### 8.5 Audit 1C-g — Regional consolidation of independent anatomical registration QC (v2.10, 31 October 2026; descriptive)
**8.5.1 Specification (approved 31 Oct 2026; D1–D5 as proposed; full text: `Audit1C_g_QC_consolidation_spec_DRAFT_20261031.md`).** Re-tabulation of existing Amendment 2 (Audit 1C-e) and Audit 1C-f outputs under the adopted registration (corrected for 97 patients; supplied for train_0018, train_0049, train_0067), with supplied values as companion columns. One row per patient × region (head/neck, thorax, abdomen, pelvis; 400 rows). Evidence classes kept separate: HB held-out bone (ribs, sternum, clavicles, scapulae; C1–C7 in head/neck), HR brain, V vessels and M muscles (Audit 1C-f, region-restricted), PV physiologically variable organs (lungs, heart, liver, spleen, kidneys, bladder; supporting only), F fitting structures (in-sample; context only). Aorta taken from Audit 1C-f only; prostate and femora excluded. Metrics never pooled: centroid distance and whole-structure MSD (wsMSD; Amendment 2), region-restricted MSD (rrMSD) and vessel in-plane displacement (Audit 1C-f). Per class and metric: n always reported; median, second-largest (n ≥ 2 only) and max. Fixed anatomical assignment of whole-structure measures to regions. No thresholds, categories, composite scores or tumour data.

**8.5.2 Lock and execution.** Script `02_code\Audit1C_g_QC_consolidation.py`, MD5 24f205e6b66c35d644794fd4d0489d36 (author Get-FileHash match). Definitive run 20261031_225909: status COMPLETE; 46/46 checks passed (100 patients; fallback set exact; every structure assigned or excluded; Audit 1C-f patient-region summaries and Amendment 2 per-patient set medians and counts reproduced with 0 differences; 400 rows, no duplicates; no pooled-MSD or classification columns). Input MD5s identical to the verified files (Section 12, v2.10 lock entry). Output table MD5 (as received) dff12596bc1d0b3b760ab5b4ed586291.

**8.5.3 Descriptive findings (patient-level medians across the cohort; adopted [supplied]; mm).**
- Evidence availability: no independent evidence in 7 patient-regions — abdomen and pelvis of train_0007, train_0059, train_0091 (undefined boundaries) and head/neck of train_0082; thin head/neck evidence in train_0098 (brain only), train_0028 and train_0042 (two cervical vertebrae). Thorax: vessels ≤ 2 per patient, no muscular evidence; independent evidence mainly held-out bone. Abdomen: no held-out bone; vessels and muscles only. Pelvis: muscles ≥ 3 in 97 patients. Region boundaries: thorax 90 direct / 9 substituted / 1 undefined; abdomen 90 / 7 / 3; pelvis 97 / 0 / 3.
- Trunk under the adopted registration: thorax HB centroid median 5.2 [8.0], P90 8.5 [20.2]; HB wsMSD 1.8 [3.1]; thorax V rrMSD 2.4 [3.2]; abdomen V rrMSD 1.8 [3.3], M rrMSD 1.6 [2.9]; pelvis V rrMSD 2.1 [2.6] (max 17.3), M rrMSD 2.4 [2.9] (max 9.7).
- Head/neck: C1–C7 centroid median 7.4 [9.4], P90 20.6, max 47.4 (> 15 mm in 21 patients); brain centroid median 15.3 [14.0], P90 45.0, max 117.3 (> 20 mm in 31 patients), with correspondingly large wsMSD; brain residuals exceed cervical residuals, compatible with head rotation/flexion between examinations (not proven). Head/neck is largely unaffected by the trunk correction.
- PV organs: pelvic (bladder) centroid median 9.2, max 40.7 — consistent with filling differences; supports supporting-only status.
- Metric scales differ (e.g., thorax HB centroid 5.2 vs wsMSD 1.8 mm in the same patients); tolerances must be metric-specific.
- train_0089 has the largest residuals in every region in which it is evaluable.
- `n_independent_structures` is not comparable across regions (ribs dominate the thorax count); evidence sufficiency must be judged per region and evidence class.

**8.5.4 Author decisions (31 Oct 2026).** Audit 1C-g accepted as complete. Head/neck is evaluated separately from the trunk with an independent acceptance strategy; head/neck residuals are not to be concealed in whole-body summaries. Head/neck lesions are NOT automatically classified as indeterminate. Step 3 literature review approved before any numerical criterion. No acceptance criteria, classification, tumour data or discordance analysis.

### 8.6 Step 3 preparatory physics note — sampling geometry and threshold-dependent boundary displacement (v2.10; note v2 ACCEPTED PROVISIONALLY 1 November 2026; τ NOT LOCKED)
**8.6.1 Scope and status.** Physics basis for a matching tolerance, established before any numerical registration-acceptance criterion. Independent of registration residuals, tumour masks and tumour correspondence. Results are **theoretical estimates from an idealised model under stated assumptions, not measured lesion-boundary errors**. Ranges are reported; no tolerance is proposed. Note: `PAPER3_PhysicsNote_sampling_boundary_v2_20261101.md` (supersedes v1 of 31 Oct 2026). Script `PhysicsNote_boundary_sampling_model_v2.py`, MD5 c015d76e3d8f0678834b10e6492df466 (v1 script MD5 cc08a8288ebb10ee9d58b60781a77ed8 retained, superseded). Inputs: Audit 1D grid geometry (run 20261030_215034); dataset SUV thresholds `SUV_thresholds_20261101_073432.csv` (MD5 bf7ce5d3174eab228ab992381506655b), extracted read-only by `Extract_SUV_thresholds.py` (MD5 943955e1176b49d6618b7230ad45d88a; 200/200 files OK; status COMPLETE). Self-tests passed (analytic profile vs brute-force 3D convolution, maximum error 0.008 on a lesion-to-background difference of 9; large-lesion limit to 0.001 mm).

**8.6.2 Literature context.** AAPM TG-132 (Brock et al. 2017; verified by the author): Table 3 — registration QA tolerances around the maximum voxel dimension (typically 2–3 mm); Table 4 — visual agreement within 1–2 voxels; radiotherapy registration-QA recommendations, not PET lesion-matching tolerances. EARL2 specifies recovery-coefficient bands, not spatial resolution or reconstruction filters. Threshold-segmented PET boundaries depend on resolution, size, contrast and partial volume; no universal boundary-displacement allowance exists. DEEP-PSMA: PSMA SUV ≥ 3; FDG liver-based threshold supplied per study; EARL-compliant reconstructions; no measured per-scanner resolution available.

**8.6.3 Dataset thresholds.** PSMA stored thresholds 2.9978–3.0037 (exactly 3.0 in 5; deviation negligible; model uses 3.0). FDG thresholds median 3.18 (P10–P90 2.38–3.96; range 1.96–6.02); 15 below 2.5 and 1 above 4.5. train_0064 has an unusually high FDG threshold (6.02): recorded as a dataset value only; no cause inferred and no further investigation on this basis.

**8.6.4 Component A — sampling geometry (matching grid = PSMA PET).** Maximum voxel dimension median 4.00 mm (P10–P90 3.27–4.20; range 2.79–5.47); 3D voxel diagonal median 6.10 mm (range 3.75–8.40); in-plane x = y; no obliquity; slice spacing 1.5–5.0 mm. Floor of about one maximum voxel dimension.

**8.6.5 Component B — Gaussian PSF boundary model.** Uniform sphere (D 6–50 mm) in uniform background, isotropic Gaussian PSF, absolute threshold. Measured inputs: PSMA T = 3.0; FDG T = the 100 per-study thresholds. Assumptions (not measurements): FWHM 4/6/8/10 mm; background PSMA 0.5/1.0/2.0, FDG 1.0/1.5/2.0; true uptake PSMA 5/10/20/40, FDG 4/6/10/20. 224 scenarios with background ≥ threshold (train_0034, train_0092 at FDG background 2.0) reported as boundary undefined. Findings (lesions ≥ 15 mm, FWHM 6 mm): δ median +2.4 mm PSMA (P10–P90 −0.6 to +3.9), +1.3 mm FDG (−1.8 to +3.3); extreme FDG shrinkage (to −9 mm) only when true uptake barely exceeds threshold. Partial-volume loss: 5,926/33,712 defined scenarios not segmented (mainly D ≤ 10 mm at FWHM ≥ 8 mm or high FDG thresholds). Identical co-located lesion, level-matched (primary): |δ_PSMA − δ_FDG| median 0.8 (FWHM 4) to 1.7–2.1 mm (FWHM 10), P90 1.6–4.4 mm; all-level combinations reported only as sensitivity analysis (median 1.2–3.2, P90 3.1–7.9 mm). Contact distance for 6–10 mm lesions: P10 4.2–8.7 mm (minima 1.2–4.8 mm in near-threshold cases); ≥ 15 mm lesions: P10 ≥ 10.5 mm.

**8.6.6 Limitations.** Spherical homogeneous lesions, uniform background, isotropic PSF; no noise, motion or voxelisation in Component B; absolute thresholds; identical true lesion geometry in both tracers by design; FWHM, background and uptake assumed; PSMA–FDG level-matching is a convention.

**8.6.7 Implications (for Step 3; nothing locked).** τ to be built from separately documented terms: sampling floor (about one maximum voxel dimension; median 4.0 mm, up to 5.5 mm) and a physics boundary term (typically 1–2 mm between tracers; up to about 3–4 mm P90 at poor resolution). Registration error stays outside τ. Candidate τ sensitivity analyses: sampling floor alone; floor + median tracer difference; floor + P90 tracer difference. **The way these terms combine (linear, root-sum-of-squares, or simulated percentile) requires mathematical justification before τ is defined; none adopted.** One reference matching tolerance with size-stratified sensitivity analyses; size alone does not justify size-dependent matching or assignment to V_I.

**8.6.8 Author decisions (31 Oct – 1 Nov 2026).** Physics first, criteria second; physics kept separate from registration residuals and tumour correspondence. Actual FDG thresholds replace the assumed range. Physics note v2 accepted provisionally, as theoretical estimates. No additional physics simulations at present. Head/neck: cervical vertebrae primary anatomical QC evidence; brain displacement a separate posture warning; neither automatically establishes acceptance. Matching tolerance (τ) and anatomical registration acceptance remain separate; neither may be adjusted to maximise the number of assessable patients. τ remains unlocked until its mathematical definition is finalised. Next: acceptance-framework specification using existing anatomical QC; descriptive evaluation of candidate thresholds across all 100 patients; rules selected and locked before any tumour correspondence is examined.

### 8.7 Audit 1C-h — Regional registration-acceptance framework (v2.10, 1 November 2026; design APPROVED, specification PROVISIONAL; nothing locked)
**8.7.1 Status.** (h1, h1b and h2z executed 1 Nov 2026; h3 specified, scripted and preflighted (100 patients); Stage 1 executed, then h3 abandoned before reading; registration investigation closed; Audit 1C-h v3.1 LOCKED; see 8.7.5–8.7.21.) Overall design approved by the author (1 Nov 2026); specification `Audit1C_h_acceptance_framework_spec_DRAFT_v2_20261101.md` remains provisional until D1–D8 are checked against their exact wording. No acceptance threshold locked; no classification, tumour data, correspondence or discordance analysis authorised. Approved now: h0 (this record) and drafting of the h1 script. Central principle: registration acceptance must reflect anatomical reliability, not the proportion of patients retained; no rule may be changed because many or few regions fail.

**8.7.2 Approved design elements.** Sequence h0 → h1 (known-translation test of metric behaviour) → h2 (segmentation-quality screen using existing pre-specified flags; PSMA/FDG volume ratio reported only, never an exclusion criterion, because 0.95–1.05 was a resampling check, not a validated segmentation criterion) → h3 (categories under all predefined candidates) → h4 (blinded face-validity review) → h5 (author review and lock). Categories: SATISFACTORY, BORDERLINE, UNACCEPTABLE, LIMITED_EVIDENCE (single-structure evidence only; never SATISFACTORY), NO_EVIDENCE; a class with n = 1 cannot alone make a region UNACCEPTABLE. Decision evidence: independent classes only (head/neck: cervical vertebrae; brain a posture warning only). If a region fails QC (UNACCEPTABLE or NO_EVIDENCE), its entire eligible FDG burden becomes V_I, decided before and independently of matching. Face-validity sample fixed before categories are computed (stratified random, ≤ 4 per main category and ≤ 2 LIMITED_EVIDENCE per region, ≤ 56 overlays, seed 20261101; blinded CT-only overlays; aligned / questionable / misaligned); it is a face-validity assessment only, not proof of registration accuracy, and a sample of 56 may miss important failures.

**8.7.3 h0 — Pre-specified rationale for the primary candidate (A1 + R1), recorded before any category or count is computed.**
- A1 (satisfactory e ≤ 1·v; borderline 1·v < e ≤ 2·v; unacceptable e > 2·v; v = the patient's maximum PSMA PET voxel dimension) is a **study-specific engineering criterion, not a physical guarantee** and not validated by TG-132 for PET/CT lesion matching.
- Consequence-based framing: registration error matters through its potential to displace lesion masks on the PET matching grid, so the PET voxel is used as the scale of consequence. A displacement smaller than one PET voxel can still change voxel membership, particularly near boundaries; A1 bounds the magnitude of this effect but does not exclude it.
- TG-132 (Table 3: tolerance around the maximum voxel dimension; Table 4: 1–2 voxels) supports voxel-scaled QA criteria by analogy only.
- Alignment is measured on CT-derived structures whereas v is a PET-grid quantity: an explicit study-specific adaptation.
- CT-derived residuals include segmentation noise; the in-sample fitting-structure residual under the adopted transformation (thorax centroid median 2.5 mm) is reported beside the anchors as a rough lower bound of that noise, not used to set them.
- R1 is primary for conservatism: a region is accepted only if all available independent evidence supports it (with the n = 1 and LIMITED_EVIDENCE rules).
- A2 (0.5·v / 1·v) and A3 (3 / 6 mm) and R2 (designated primary class per region) are predefined sensitivity candidates.
- **Commitment:** the primary candidate will not be changed because of the number of regions that fail; changes are permissible only through a dated amendment for a documented validity reason (h1 or h4).

**8.7.4 h1 clarifications (author, 1 Nov 2026).** The vessel P95 ≤ 1.0 mm allowance is a provisional implementation check, justified only because the vessel in-plane displacement is derived from voxelised masks (FDG labels resampled by nearest neighbour onto the PSMA grid; slice-wise centroids; the sampling window depends on the FDG z-extent). Where displacement is derived from transformed anatomical points, its change must equal the applied translation to numerical precision. The centroid bound |t| + 0.5 mm is a weak sanity check; vector changes are to be tested directly where centroid vectors are available.

**8.7.5 h1 — Known-translation test: execution and results (run 20261101_080246; status COMPLETE).** Script `02_code\Audit1C_h1_known_translation_test.py`, MD5 bc44924d1b7f85e12c458c857b1f8e9e (author hash confirmed); synthetic self-test 9/9; inputs Audit1C_f_segments (MD5 a8087e0ade5acb04ccf3ad0681d6eb3f), Audit1C_e_per_patient (c66da788a6c73737163b4c09ed76f494), Audit1C_e_structure_level (289b57e69a6f119968b3ae3990db9fd9); implementation checks 4/4 passed (fallback patients: T_A identical to T_S). A pandas deprecation warning (DataFrame.join lsuffix/rsuffix) was printed; it does not affect results. Expected change of every FDG-minus-PSMA measure between corrected and supplied transformation: −t.
- Bone (held-out bone and C1–C7; 3,468 pairs): 0 violations of the triangle inequality at 10⁻⁶ mm (maximum excess exactly 0), consistent with the exact −t vector change by construction (Amendment 2 check (b)).
- Vessel in-plane vector (740 pairs; 12 segments evaluable under one transformation only): deviation from −t_xy median 0.45 mm, P95 3.36 mm, max 13.69 mm → provisional criterion P95 ≤ 1.0 mm NOT met. Sign confirmed by comparators (opposite sign median 8.68 mm; no change 4.49 mm). P95 by region: thorax 1.07, abdomen 2.65, pelvis 4.23 mm; slice set unchanged 2.33, changed 5.33 mm.
- Muscle region-restricted centroid distance (951 pairs): 5 violations of |Δe| ≤ |t| + 0.5 mm (maximum excess over |t| 1.09 mm).

**8.7.6 h1 diagnostics (descriptive; existing outputs only; no images; script `Audit1C_h1_diagnostics.py`, MD5 9a68fd325e22bab96de32b5cbc43b094; results reproduced from MD5-identical inputs).**
- Vessels: residual deviation correlates with |t_z| (Spearman ρ 0.68) and not with |t_xy| (ρ 0.04); secondarily with sampling-window change (ρ 0.40). Per-vessel regression of the residual on t_z: iliac arteries and veins left–right slope −0.48 to −0.44 (left) and +0.30 to +0.42 (right), R² 0.78–0.90 (mirror-symmetric, consistent with laterally diverging oblique vessels); aorta and IVC slopes 0.03–0.08. P95 deviation by |t_z|: ≤ 1 mm 0.88; 1–3 mm 1.90; 3–6 mm 3.52; > 6 mm 7.74 mm. After removing the per-vessel t_z slope, P95 remains 1.86 mm; near-vertical aorta and IVC P95 1.24–1.37 mm; largest outliers in short segments (e.g., train_0041 iliac veins, 21 mm window, 7 slices). Mechanism: through-plane translation of oblique vessels appears as in-plane displacement, plus voxel-sampling and sampling-window effects.
- Muscles: violations (train_0045 gluteus medius R, train_0052 gluteus medius R and iliopsoas R, train_0069 gluteus medius L, train_0073 gluteus medius R) are not explained by |t_z| or window change (3/5 without window change); mechanism not identified (field-of-view truncation an untested hypothesis requiring images).

**8.7.7 Author conclusions (1 Nov 2026).**
1. Bone centroid metrics passed h1 (a test of metric behaviour; exact by construction) — retained as decision evidence.
2. Vessel in-plane displacement failed because vessel orientation and voxel sampling distort the displacement measurement — demoted to supporting evidence, as pre-specified. The 1.0 mm criterion is not relaxed.
3. Muscle centroid distance: five unexplained violations — remains demoted pending further evidence.
These are failures of the measurements for deciding registration acceptance under the current rules, not evidence that the registrations failed. h1 thresholds are not relaxed retrospectively; no region is classified; LIMITED_EVIDENCE is not applied by default to abdomen or pelvis.

**8.7.8 Next step (author).** Critically assess whether the existing region-restricted MSD (Audit 1C-f vascular and muscular surface measurements) can provide defensible anatomical acceptance evidence, before any new image-processing method is developed. No classification and no threshold changes.

**8.7.9 h1b — Descriptive assessment of region-restricted MSD (plan approved with P1–P4; run 20261101_082311; status COMPLETE).** Script `02_code\Audit1C_h1b_rrMSD_assessment.py`, MD5 36497eaf36397c4cc5f1467bd8ceef50; synthetic self-test 8/8; technical checks 3/3 (fallback set; fallback T_A identical to T_S for MSD/HD95; vessel pairing reproduces h1 exactly). Inputs: Audit1C_f_segments (a8087e0ade5acb04ccf3ad0681d6eb3f), Audit1C_e_per_patient (c66da788a6c73737163b4c09ed76f494), Audit1C_e_structure_level (289b57e69a6f119968b3ae3990db9fd9), Audit1C_h1_vessel_pairs_20261101_080246 (fb558fae4d755096261af0f3fb69e8eb). 1,691 pairs (740 vessel, 951 muscle); 368 aorta/IVC reference pairs; 99 patients for the thorax comparison. Safeguards: MSD measures surface agreement, not necessarily displacement; a low MSD never qualifies a region; a high MSD may reflect physiological deformation or segmentation differences.
- Q1 (aorta/IVC reference; imperfect, h1 P95 1.24–1.37 mm; not independent of MSD because both derive from the same masks): ΔrrMSD vs Δ in-plane displacement Spearman 0.98, sign concordance 0.95 (1.00 where |Δ reference| ≥ 1.4 mm); HD95 0.96 / 0.92. Across all vessels and muscles, |ΔrrMSD| correlates with |t_xy| (ρ 0.68–0.86) and not with |t_z| (ρ −0.06 to 0.14). Attenuated scale: ΔrrMSD about 0.2–0.4 mm per mm of shift (HD95 about 0.6). Floor (aorta/IVC states with in-plane displacement < 1 mm; n = 9): rrMSD median 0.65 mm.
- Q2: no material dependence of rrMSD on window length or surface-point count; short windows (20–30 mm) noisier (median |ΔrrMSD| 2.2–2.9 vs about 0.9 mm).
- Q3 (thorax): rrMSD vs held-out bone centroid — supplied ρ 0.86 (median) / 0.57 (second-largest); adopted ρ 0.40 / 0.24; change ρ 0.81 (sign concordance 0.80) / 0.56 (0.59). Pelvic fitting structures vs pelvic rrMSD ρ 0.72–0.76: IN-SAMPLE CONTEXT ONLY.
- Q4: aorta/IVC axial-dominant shifts (n = 8) |ΔrrMSD| / |t| median 0.08 vs in-plane-dominant (n = 112) 0.36; muscles 0.18.

**8.7.10 Author conclusions on h1b (1 Nov 2026).**
1. Region-restricted MSD is useful for detecting substantial in-plane disagreement but is not suitable as a standalone registration-acceptance metric (blind to cranio-caudal displacement; attenuated, non-displacement scale; agreement with bone weak at the small residuals of the adopted state).
2. Cranio-caudal alignment remains insufficiently assessed by independent evidence in the abdomen and pelvis.
3. No threshold changes: the failed metrics remain demoted; no region is classified. Option 1 (accepting the cranio-caudal gap as a limitation) is not accepted at this stage, because unrecognised cranio-caudal misalignment could undermine the PSMA–FDG discordance analysis.

**8.7.11 Next step (author).** Option 2: examine existing per-vertebra residuals for evidence of cranio-caudal displacement, distinguishing (i) absolute residual displacement under the adopted transformation and (ii) local variation around the fitted translation. Because these vertebrae contributed to fitting, neither constitutes independent validation. If existing data cannot resolve the gap, consider Option 3 (a targeted independent measurement, by dated amendment).

**8.7.12 h2z — Cranio-caudal evidence from existing bone centroid vectors (specification approved with Z1–Z5 and author safeguards; run 20261101_084243; status COMPLETE).** Script `02_code\Audit1C_h2z_craniocaudal.py`, MD5 c12377256f9413d2311ae06d0bf4986e (author hash confirmed); synthetic self-test 10/10. Inputs: Audit1C_b_structure_level_20261029_194714 (MD5 360d0a7af1edfb1f58d789a7037f807b), Audit1C_b_patient_level_20261029_194714 (71201ec32fbd12f922ddeebb54d06c80; context only), Audit1C_e_per_patient (c66da788a6c73737163b4c09ed76f494), Audit1C_e_structure_level (289b57e69a6f119968b3ae3990db9fd9). Adopted residual r = d − t (97 adopted patients; exact because T′(p) = T(p) − t); r = d for the 3 fallback patients (supplied = adopted); fallback corrected state and fallback leave-region-out reported separately as HYPOTHETICAL. Gradients: Theil–Sen slope against the PSMA centroid z (physical mm) over the patient's Amendment 2 fitting vertebrae T1–L5 and S1 (which exclude Audit 1B flagged levels, their neighbours and levels adjacent to fragmented vertebrae). Leave-region-out: same component-wise median estimator and the patient's original Amendment 2 fitting list; estimable only with ≥ 5 remaining structures.
- **Sign convention (found during technical validation):** Audit 1C-b stores centroid displacement as PSMA minus FDG-mapped, the opposite of the Amendment 2 convention. The stored vectors therefore reproduced −t exactly (errors up to 45 mm with the uncorrected sign). The script converts explicitly (d = −1 × stored) and check C2 verifies the converted convention in every patient.
- **Validation checks (9/9 passed):** (1) fallback set = train_0018, train_0049, train_0067; (2) one Audit 1C-b row per patient × structure; (3) |d| equals Audit 1C-b centroid_distance_mm (10⁻⁶ mm); (4) |d| equals the Amendment 2 supplied centroid distance for all shared structures (0.01 mm); (5) every Amendment 2 fitting structure has a d vector; (6) held-out bone / cervical structures without a d vector listed (341: scapulae 181, clavicles 160 — not in Audit 1C-b; ribs and sternum complete); (7) component-wise median of d over each patient's fitting structures reproduces t (10⁻⁶ mm); (8) fitting-structure names equal each patient's fitting list; (9) context only: Spearman ρ −0.04 between the Theil–Sen r_z slope and the Audit 1C-b vertebral gradient (definitions differ; not a check of agreement). 100 patients; 5,074 structures (2,467 independent). Leave-region-out estimable: pelvis 100; abdomen 99 (1 without lumbar fitting structures); thorax vertebrae 98 (2 with < 5 remaining structures; counted, not omitted).
- **Results (adopted patients; mm):**
  - Spine gradient (translation-invariant): Theil–Sen r_z median 0.12 mm/100 mm (P10–P90 −0.54 to 0.78; maximum 2.44, train_0078).
  - Thorax ribs and sternum (INDEPENDENT): signed median r_z median 0.19 (P10–P90 −2.79 to 2.62); median |r_z| 2.32; second-largest |r_z| 6.34; |signed median| > 3 mm in 17 patients and > 5 mm in 5.
  - Abdomen L1–L5 (IN-SAMPLE): median |r_z| 0.48; leave-abdomen-out median |r_z| 0.67 (8 patients > 2 mm; train_0091 4.57 mm on 2 vertebrae).
  - Pelvis S1, sacrum, hips (IN-SAMPLE): median |r_z| 1.25; leave-pelvis-out median |r_z| 1.38, with |signed median| > 2 mm in 25/97, > 3 mm in 11 and > 4 mm in 4: train_0078 (−5.13), train_0090 (−4.21), train_0087 (−4.17), train_0069 (+4.14). Leaving the pelvis out shifts t by a median 0.83 mm. Pelvic-deterioration patients (Audit 1C-f): train_0069 +4.14; train_0058 +2.99; train_0037 −0.79 (deterioration in-plane); train_0064 +1.21; train_0068 +1.12; train_0089 −2.29.
  - Leave-thorax-vertebrae-out: median |r_z| 0.91 (least stable estimate, as pre-stated).
- **Limitations:** lumbar, sacral and hip residuals are in-sample; leave-region-out is internal cross-validation, not independent anatomical validation; regional estimates rest on few structures (pelvis 3–4; abdomen ≤ 5), which are **correlated anatomical observations, not independent measurements**; rib and sternum z-residuals are independent but affected by respiratory motion and by rib segmentation extent; clavicles and scapulae unavailable; bone alignment cannot establish soft-tissue (nodal) alignment.

**8.7.13 Author conclusions on h2z (1 Nov 2026).**
1. Bone evidence: cranio-caudal alignment is generally reassuring, but several patients have meaningful regional offsets.
2. Pelvic exceptions are preserved individually, particularly train_0078, train_0090, train_0087 and train_0069.
3. The soft-tissue cranio-caudal gap remains unresolved; bone alignment cannot establish nodal alignment.
4. No regional classification; h2z is supporting evidence only.

**8.7.14 Next step (author).** Option 2: draft a narrowly focused specification for targeted independent soft-tissue landmarks (aortic bifurcation; iliac venous confluence) to provide cranio-caudal evidence in the abdomen and pelvis. Their detectability and reproducibility must be established before they are relied upon. New image processing requires a dated amendment.

**8.7.15 h3 — Targeted vascular landmarks: specification v2 and script (approved for script development 1 Nov 2026; nothing executed beyond preflight).** Specification `Audit1C_h3_vascular_landmarks_spec_DRAFT_v2_20261101.md` (concept and targeted scope approved; author decisions H1–H5): LB aortic bifurcation and LV iliac venous confluence; automated candidates (Def1 parent-vessel caudal end; Def2 branch junction; all connected components; agreement = internal consistency only, not independence; fragmentation recorded, never used to exclude); accuracy against blinded manual identification on native CT (one reader, the author; inter-reader reproducibility not assessed; own viewer with coded native CTs and MPR; both CTs of each sampled patient; identification criterion = anatomical junction tracked across consecutive axial slices and confirmed on MPR, not the first appearance of two outlines); sample 20 random patients (seed 20261101) + all with Def1–Def2 disagreement > 20 mm + 10 repeat CTs; landmark reliability established before any residual and never by agreement of residuals with bone or the transform. Provisional L-conditions (landmark reliability only, not registration-acceptance thresholds): L1 automated vs manual |Δz| median ≤ 1.5 mm and P90 ≤ 3 mm; L2 repeat-read 95% limits of agreement within ±2 mm; L3 CT z-spacing ≤ 3 mm on both CTs for interpretation (else precision-limited); L4 correct anatomical identification on blinded manual review (pass B; provisional ≥ 90%). Local evidence only (about L4–S1); landmark reproducibility does not guarantee registration accuracy elsewhere in the abdomen or pelvis. Script `02_code\Audit1C_h3_vascular_landmarks.py`, MD5 6afbd3409c16326393522413f6ba2a7d (stages: 1 candidates, sample and coded review sets; 1b pass-B images after locked readings; 2 reliability only; 3 residuals only for landmark types meeting all L-conditions); synthetic self-test 16/16 (author's PC).

**8.7.16 h3 preflight (2 patients; nothing written; 1 Nov 2026).** Image-reading code functional on real data (train_0001, train_0002; both CTs). LB: both definitions present on all 4 CTs, Def1–Def2 distance 4.9–5.2 mm, parent vessel in one component, both branches reaching the parent end. LV: Def1–Def2 distance 13.7 mm (train_0001) and 22.0–22.4 mm (train_0002); left iliac vein not reaching the IVC caudal end within 10 mm in both patients and neither branch in train_0002 — consistent with a gap between IVC and iliac-vein labels or with left common iliac vein anatomy; LV may therefore be weakly defined (to be decided by the L-conditions, not by tuning). CT z-spacing 3.27 mm on all 4 CTs, above the provisional L3 limit (≤ 3 mm): if typical, most patients would be precision-limited. Open before Stage 1 (which copies coded native CTs and requires explicit author permission): (i) sample size and reading load (number of patients with Def1–Def2 disagreement > 20 mm); (ii) CT slice-spacing distribution and the L3 decision (geometric, not result-driven); (iii) confirmation of L1–L4 values. A full candidate preflight (`--stage1 --preflight 100`, nothing written) is in progress to inform these decisions.

**8.7.17 h3 full candidate preflight (100 patients; both CTs; nothing written; 1 Nov 2026).** Landmark-quality descriptors only (no residuals).
- LB (aortic bifurcation): Def1–Def2 distance median 5.8 mm (P90 9.4; min 3.0); both iliac arteries reaching the aortic caudal end in 98% of CTs; > 20 mm in 2 patients (train_0051 about 20–21 mm; train_0024 234 mm).
- LV (iliac venous confluence): Def1–Def2 distance median 12.2 mm (P90 19.7); both iliac veins reaching the IVC caudal end in only 18% of CTs (a systematic gap between IVC and iliac-vein labels — a methodological finding about label continuity); > 20 mm in 11 patients.
- train_0024: Def1–Def2 about 231–234 mm (LB on the FDG CT; LV on both CTs) — almost certainly a stray label component captured by the all-components rule; retained unchanged; the predefined manual assessment will identify the cause.
- All landmarks ≥ 10 mm from the field-of-view ends.
- CT z-spacing (200 CTs): 1.0 mm 12; 1.5 mm 6; 2.0 mm 26; 2.5 mm 14; 2.79 mm 9; 3.0 mm 27; 3.27 mm 88; 4.0 mm 12; 5.0 mm 6. Patients with both CTs ≤ 3.0 mm: 41; ≤ 3.27 mm: 89; 11 with a 4–5 mm CT.
- Reading load: 12 patients with Def1–Def2 > 20 mm on either CT for either landmark; with the 20 random patients (seed 20261101; 1 overlapping) the manual sample is 31 patients = 62 CTs, plus 10 repeat CTs (72 CT reads, two landmarks each); about 124 pass-B images.
- Caution (author): agreement between automated definitions (LB median 5.8 mm) does not establish agreement with the true anatomical junction and is not reassuring against the 1.5 mm L1 target; manual comparison remains essential.

**8.7.18 Author decisions on h3 before Stage 1 (1 Nov 2026).**
1. L3 = CT z-spacing ≤ 3.3 mm on both CTs, documented as including the existing 3.27 mm acquisition protocol. L3 is an **eligibility criterion for evaluating landmark residuals, not evidence that a 3.3 mm CT can reliably resolve every 3 mm displacement**; actual measurement precision must still pass L1/L2.
2. Automated LV is retained for evaluation, without assuming it will pass; its poor label continuity is recorded as a methodological finding.
3. train_0024 is retained unchanged; the algorithm is not modified because of one observed case.
4. Stage 1 approved: creation of 72 coded CT copies exclusively within 06_output; original data untouched (read-only). Blinded reading approved.
5. Because the manual sample is enriched for automated disagreements, its accuracy statistics are not unbiased estimates for all 100 patients: random-sample and disagreement-enriched results are reported separately. Pre-specified decision basis: L1 and L4 from the random stratum; L2 from all repeat pairs (reader repeatability; the repeat set is too small to stratify); all strata reported.
6. No residual interpretation, regional acceptance classification or lesion correspondence.
Script revised before any Stage 1 execution to implement decisions 1 and 5 (L3 = 3.3 mm; stratified reliability with a decision row); no change to candidate definitions or sampling: `02_code\Audit1C_h3_vascular_landmarks.py`, MD5 a11b23f2d39352372282ccbc9172440f (supersedes 6afbd3409c16326393522413f6ba2a7d); synthetic self-test 17/17.

**8.7.19 h3 Stage 1 executed (run 20261101_104523; 1 Nov 2026).** Script MD5 a11b23f2d39352372282ccbc9172440f (author hash confirmed; self-test 17/17 on the author's PC). Automated candidates computed for all 100 patients (both CTs); fixed sample of 31 patients (20 random, seed 20261101, plus all with Def1–Def2 > 20 mm; 1 overlapping) → 62 coded native CTs for pass A (`Audit1C_h3_review_passA_20261101_104523`) and 10 re-coded repeat CTs (`Audit1C_h3_review_repeat_20261101_104523`), written exclusively within 06_output; original data read-only and untouched. Key, sample and candidate files are withheld from the reader until all readings are locked. Next: blinded pass A reading (author, 3D Slicer, RAS coordinates; identification criterion as fixed in the README), lock by MD5; then the repeat reading, locked by MD5; then Stage 1b (pass-B images). No residuals, classification or correspondence.

**8.7.20 Closure of the registration investigation and of h3 (author decision, 1 November 2026).**
- **Decision:** development of new registration measurements stops. h3 (manual vascular landmark assessment) is abandoned before any reading, because appropriately qualified independent reading is unavailable. This revises the earlier decision not to accept Option 1 (Section 8.7.10): the soft-tissue cranio-caudal gap in the abdomen and pelvis is now accepted as an **unresolved, explicitly acknowledged uncertainty**, handled by a conservative acceptance framework (Audit 1C-h v3) and a stricter sensitivity analysis, not by further measurement.
- **Rationale:** the completed audits (Amendments 1–2; Audits 1C-a to 1C-g; h1, h1b, h2z) provide substantial anatomical evidence; the remaining measurement would have resolved only local (about L4–S1) cranio-caudal uncertainty at a substantial reading cost; further sophisticated registration measurements are not expected to change the study materially.
- **Status of h3 outputs:** Stage 1 (run 20261101_104523) was executed but no reading was performed. The automated landmark candidates are **unvalidated and must not be used** as evidence of registration accuracy or in any analysis; they are archived only. The 72 coded CT copies remain in 06_output as an unused archived artefact (deletion only with separate author permission). Stages 1b–3 were not run.
- **Retained:** all completed registration audits and their findings, with their documented limitations (vessel and muscle displacement metrics demoted at h1; region-restricted MSD supporting only, blind to cranio-caudal displacement, h1b; bone cranio-caudal evidence in-sample or internal cross-validation only, h2z).
- **Methodological principle (author):** registration confidence and biological discordance are different uncertainties. A PSMA-negative/FDG-positive lesion is not interpreted as biologically discordant merely because its region passes anatomical QC; FDG specificity and lesion characterisation remain separate considerations.
- **Next:** Audit 1C-h specification v3 (three confidence levels; conservative primary analysis; stricter sensitivity analysis; limited gross visual QC), locked before any regional classification is computed; then regional QC; completion of the pending Audit 1D documentation; protocol lock; then lesion correspondence.

**8.7.21 Audit 1C-h v3.1 — LOCKED specification (author approval, 1 November 2026).** Specification `Audit1C_h_acceptance_framework_spec_DRAFT_v3_1_20261101.md` (supersedes v1–v3) is locked, with all V1–V6 decisions as approved:
- **Levels:** HIGHER (independent anatomical evidence supports acceptable regional alignment — higher confidence, not proven sub-voxel registration accuracy); QUALIFIED (supporting measurements reassuring, independent validation incomplete — never described as independently validated); INSUFFICIENT (inadequate evidence or concerning anatomical disagreement).
- **Fixed sequence (author safeguard):** evidence availability → bone limits → MSD warnings → confidence category → visual QC. Missing or unreliable evidence is never interpreted as a passing measurement.
- **V1:** region-restricted MSD > 0.5·v is a provisional warning, never INSUFFICIENT by itself; caps at QUALIFIED; triggers visual QC; INSUFFICIENT only on visual "grossly misaligned"; unreviewed or inconclusive warnings are retained; a reassuring review removes nothing.
- **V2:** abdomen/pelvis bone cranio-caudal (|leave-region-out signed median r_z|) or in-plane (leave-region-out median |r_xy|) > 1·v → INSUFFICIENT, with no borderline zone; fallback patients use in-sample supplied-state values. **V2a:** n ≥ 2 within a class is minimum evidence only; correlated observations; no independent confirmation implied.
- **V3:** borderline (1–2·v) or single-structure HB → QUALIFIED. **V4:** strict sensitivity — all QUALIFIED abdominal and pelvic eligible burden → V_I. **V5:** rule-defined visual QC sample (all V1 warnings; all single-structure/borderline-HB QUALIFIED head/neck and thorax regions; ≤ 4 HIGHER and ≤ 4 QUALIFIED per region type at random, seed 20261101); downgrade-only; workload reported before any review; impractical workload documented as a limitation applying to all unreviewed regions, never selective review. **V6:** A2/A3 vary thresholds only on identical measurements; no recomputation or optimisation.
- **Principle 7:** the h2z and Audit 1C-g distributions had been seen descriptively; all limits derive from the voxel anchor A1 recorded in h0 (8.7.3) before h1b and h2z.
- **Commitment (author):** no changes to registration transforms, thresholds or classifications based on observed category counts. If the predefined rules classify important pelvic regions as INSUFFICIENT, that result is accepted; the rules are not changed to rescue them.
- **Classification script** `02_code\Audit1C_h_classification.py`, MD5 e4bb1a335ce2e50427acff72677f293a (synthetic self-test 22/22; evidence-extraction check on uploaded copies of the inputs: 400 rows, v available for all, Audit 1C-g table MD5 matched, abdomen/pelvis bone evidence 193 leave-region-out, 6 in-sample fallback, 1 not usable; no classification computed). Stages: `--classify` (pre-visual levels A1/A2/A3, inclusion flags, visual QC workload report and rule-defined coded sample; no overlays) and `--apply-visual` (locked ratings; downward only). Script MD5 CONFIRMED on the author's PC (Get-FileHash E4BB1A335CE2E50427ACFF72677F293A, 1 Nov 2026) and LOCKED before any patient-level classification; no code changes after this point.

## 9. Parameters to lock after audits (before examining discordance)
- Reference matching tolerance (+ sensitivity values).
- Primary connectivity (6, 18 and 26 evaluated in Audit 1D).
- Minimum lesion volume.
- Confluent-disease rule.
- Registration QC acceptance limits, including the definition of the borderline category and its predefined sensitivity analysis (Section 7); structures suitable for paired registration QC (from Audit 1B paired completeness).
- Coverage thresholds for site- and correspondence-indeterminate status (Section 8.1).
- Fragmentation criterion distinguishing trivial fragments from substantial displaced components (provisional: >5% of label volume outside the largest 26-connected component).
- Exclusion rule for unreliable vertebral landmarks (affected and adjacent levels).
- Perturbation magnitudes, number per magnitude (provisional 50) and random seed.
- Perturbation-stability criterion for the exploratory robustness measure (provisional ≥95%).
- PSMA uptake reference tissue, class thresholds, and handling of small lesions (partial volume).
- Anatomical assignment rules, including predefined nodal-region rules (unlabelled soft tissue not automatically nodal).
- Grid-fingerprint rounding precision (spatial values, provisional 1e-4 mm; direction cosines) (Section 14.5).
- Statistical analysis plan.
- (v2.7) Lesion-level matching rule (provisional: any-contact within the reference tolerance; one-to-many permitted) and distance definition (boundary-voxel centres vs voxel boundaries).
- (v2.7) FDG eligibility rule: minimum physical volume plus minimum voxel count, applied in native FDG space, separately for each connectivity (candidates in Section 16).
- (v2.7) Order of volume-category assignment (below-threshold → correspondence-indeterminate → matched/unmatched).
- (v2.7) Regional registration-QC inheritance hierarchy for lesions without nearby landmarks (local → regional → indeterminate).
- (v2.7) Rule for FDG components spanning several QC regions (candidates in Section 16).
- (v2.7) Voxel-level spatially unmatched fraction: definition, handling of indeterminate voxels, inclusion in perturbation analysis.
Lock documented in a date-stamped protocol addendum (optional OSF pre-registration).

## 10. Limitations (stated in advance)
No scan dates (interval unknown); different labelling rules per tracer (PSMA fixed SUV 3, FDG liver-based; labels restricted to suprathreshold voxels by design); rigid registration between separate acquisitions, described by the data providers as coarse and not validated for lesion-level alignment; selected single-centre population; no outcomes; annotator-dependent physiological exclusion; differing head/neck posture between the two examinations, which a single trunk translation cannot correct (head/neck worsened > 2 mm in 28 patients and > 5 mm in 9 after the adopted correction; separate regional QC); supplied rigid PSMA–FDG registration with residual anatomical misalignment in some patients (whole-body bony centroid median 7.3 mm, P90 20.2 mm; ~20 patients with near-uniform 17–23 mm in-plane offset; Audit 1C-b); TotalSegmentator has no lymph-node class (nodal assignment relies on rule-based anatomical regions); TotalSegmentator generation mode (fast/full) undocumented; arms-up positioning places parts of the upper limbs outside the CT field of view in some patients (site-indeterminate tumour); automated TotalSegmentator labels show apparent vertebral level misassignment and fragmentation in some patients with abnormal spines (identified by automated screening and limited visual inspection, not by expert review); completeness of structures is assessed by automated coverage criteria only ("apparently complete"); absent structures (e.g., gallbladder, kidney, prostate) cannot be attributed to surgery or segmentation failure without clinical data; possible overlap with data owners' own work.

## 11. Target journals [VERIFY scope]
EJNMMI Research (if robustness/uncertainty analysis substantial); EJNMMI Reports; Nuclear Medicine Communications; Frontiers in Oncology (reviewers may expect clinical correlates).

## 12. Decision log (29–30 Oct 2026)
- Anatomical distribution demoted from primary aim (novelty: Telli 2026).
- Two-definition concept (lesion-level discordance vs intralesional mismatch) → intralesional mismatch not carried forward as an objective.
- Expert adjudication excluded (not feasible).
- RDMTVF across nested tolerances rejected (= permissive DMTVF); primary endpoint = DMTVF_ref; robustness = principal methodological objective via perturbation.
- One distance-based method family for the primary endpoint; overlap method only as supplementary agreement analysis (avoids circularity).
- Native-space FDG volume for burden; transformed masks only for correspondence.
- All 100 patients enter a fresh eligibility audit; lesion-level indeterminate handling.
- PSMA uptake classes retained as secondary analysis 3.
- (v2.1) Official documentation confirms suprathreshold-only labelling and coarse, unvalidated registration → added to Data and Limitations.
- (v2.1) Official template code shows the two rigid.tfm files are an inverse pair; working hypothesis recorded: PSMA/rigid.tfm for resampling FDG onto PSMA grid, FDG/rigid.tfm for mapping FDG points into PSMA space; Audit C tests it.
- (v2.1) Bounds split into a primary uncertainty interval (anatomical; indeterminate volume) and a separate sensitivity envelope (classification; alternative matching configurations).
- (v2.1) PSMA uptake sampled from native PSMA PET at the transformed FDG lesion location; no resampling of PSMA values.
- (v2.1) Registration QC preceded by a coverage check; partially imaged structures excluded from surface-distance metrics.
- (v2.2) Final analysis workflow (Steps 1A–13, with QC split 4A/4B/4C/4D and SAP locked at Step 2) and three operating rules adopted; working folder set to D:\DEEP PSMA\PAPER3_Discordance.
- (v2.3) Hybrid software architecture and data-exchange rules adopted (Section 14): Python/SimpleITK for all spatial operations; MATLAB for quantification, matching, statistics; common PSMA PET grid; lesion-preserving integer-label resampling; grid fingerprint; JSON manifests.
- (v2.4) Paper 3 declared a fully independent study starting from the original DEEP-PSMA data: all references to previous analyses, eligibility decisions and exclusions removed; Audit 1A redefined as a comprehensive raw-data audit of all 100 cases.
- (v2.5) Indeterminate lesions split into two classes (Section 4.4): site-indeterminate (anatomy unavailable; matching assessable; included in DMTVF according to matching result) and correspondence-indeterminate (matching not assessable; excluded from numerator, retained in denominator). A lesion outside CT coverage is not automatically correspondence-indeterminate, because matching uses PET. Final correspondence-indeterminate status assigned only at Step 4D. Sections 4.4, 5, 13 and 14.3 updated accordingly.
- (v2.5) Revised coverage measures added (Section 8.1): Audit 1A PET-field criterion retired; tumour-specific and anatomically relevant measures at case level (Audit 1C) and lesion level (Step 4D); thresholds locked at Step 2.
- (v2.5) Audit 1A (report: 01_audits\Audit1A_output.md): patients are not excluded because of differences between reconstructed PET and CT fields of view; 42 FLAGs accepted as reconstruction geometry; PET-field coverage criterion replaced by tumour-specific and anatomically relevant coverage measures (Section 8.1).
- (v2.5) Audit 1A COMPLETED: all 100 cases retained; 1,200/1,200 files present and readable; geometry, threshold and tumour-mask checks passed. train_0074 FAIL resolved by targeted coverage check (02_code\Audit1A_B_train0074_coverage.m): large PSMA PET reconstruction field (1.32 m), trunk fully covered, arms-up positioning (upper limbs partly outside CT); right upper-arm PSMA lesion site-indeterminate for the portion outside CT (2.65 mL = tumour volume outside CT, not necessarily whole-lesion volume); correspondence assessable. Six tracer entries in five patients have small amounts of tumour outside CT → site-indeterminate only.
- (v2.5) TotalSegmentator dictionary verified from official sources: v2.4, task=total, 117 classes; exact label-ID-to-name mapping hard-coded; version, task and source recorded in manifests.
- (v2.5) Audit 1B scope adopted (Section 8): dictionary verification by anatomical plausibility; absent vs possibly truncated structures distinguished without clinical interpretation; PSMA and FDG CT completeness recorded separately; paired-candidate structures identified for registration QC.
- (v2.5) Unlabelled soft tissue is not automatically classified as lymph node; predefined nodal-region rules required (Sections 4.4, 5, 9).
- (v2.5) Audit 1B COMPLETED and CLOSED (report: 06_output\Audit1B_report_20261029.md; decision log: 06_output\Audit1B_decision_log_20261029.md). 200/200 CTs read; all masks on CT grid; labels only 1–117; laterality 0 and craniocaudal 0 failures, vertebral order 1 (train_0088) → official 117-class mapping consistent with anatomy. No coverage-related absences. Completeness wording: "apparently complete according to automated coverage criteria". Kidney cysts treated as optional structures. Absent structures (gallbladder 14 patients; right kidney train_0082; prostate one CT in train_0025 and train_0087; right 12th rib train_0029 FDG) reported, not interpreted clinically; prostate label not evidence of an intact gland. Reference structures for nodal-region rules and liver (uptake reference) apparently complete in essentially all patients.
- (v2.5) Targeted vertebral checks (limited visual inspection): apparent level misassignment in train_0088 (T2–T4), train_0059 (L4, with L3 and L5) and train_0029 (T12, L1, L2), on both CTs. Affected and adjacent levels flagged as unreliable landmarks for registration QC (final rule at Step 2). Spinal disease reported by cervical, thoracic, lumbar and sacral regions. Fragmentation: trivial fragments distinguished from substantial displaced components; >5% criterion provisional.
- (v2.5) Cohort after Audits 1A and 1B: all 100 patients retained; no case-level exclusions.
- (v2.5) Audit 1C split into 1C-a (ten-case pilot) and 1C-b (full cohort). Safeguard: direction tested with SimpleITK's actual resampling convention; point mapping tested separately using the inverse as appropriate.
- (v2.5) Audit 1C-a COMPLETED (report: 06_output\Audit1C_a_report_20261029.md): Section 4.1 hypothesis verified in 10/10 cases (resampling: PSMA/rigid.tfm; point mapping: FDG/rigid.tfm = inverse); exact inverse pair in all 100; resampling preserves structure volumes (ratio 0.985–1.022). Residual bony centroid distance after the verified transform 3.7–22.3 mm (4/10 pilot cases ≈ 18–22 mm; near-constant and position-dependent patterns). No patient classified as a registration failure; no exclusions; QC limits and any registration refinement decided only from the full 100-patient anatomical distribution (Audit 1C-b), never from discordance outcomes.
- (v2.5) Audit 1C-b COMPLETED (outputs 06_output\Audit1C_b_*_20261029_194714.*; report to follow with Audit 1C closure): 100/100 patients; supplied transform — whole-body bony centroid median 7.3 mm (IQR 4.9–13.7; P90 20.2); bone mean surface distance median 2.7 mm; ~20 patients with near-uniform 17–23 mm in-plane bony offset; cervical increase (posture); optimal rigid re-fit of bony centroids would leave 2.8 mm median (in-sample). FDG TTB inside PSMA PET FOV after mapping: 100% in 93 patients; < 95% in 2 (train_0044, train_0047). No patient classified or excluded.
- (v2.5) AMENDMENT 1 (29 Oct 2026; Section 15): bone-based rigid refinement experiment recorded before results, with pre-specified method, held-out evaluation and adoption criteria; adoption conditional on independent anatomical QC.
- (v2.5) Amendment 1 LOCKED (29 Oct 2026, author approval before execution): C1–C4, fitting parameters, structure sets and fallback rule locked (script MD5 0c350ab57b59c1dea8fc729b434457b7); six-step decision sequence adopted; refinement described as anatomically evaluated registration refinement, not validated lesion correspondence.
- (v2.5) Amendment 1 OUTCOME (30 Oct 2026): locked run 20261029_212531 (100/100, no errors); C1 and C2 met, C3 not met (22 patients worse > 2 mm), C4 not met (gross misregistration, e.g. train_0070). Refinement NOT adopted.
- (v2.5) DECISION (author, 30 Oct 2026): the supplied rigid transformation (resampling PSMA/rigid.tfm; point mapping FDG/rigid.tfm) is the PRIMARY transformation for all 100 patients; per-patient table 06_output\Audit1C_per_patient_transformation_QC_table_20261030.csv. All 100 patients retained; no patient-level registration-failure classification. Lesion assessability to be determined by regional/local QC, with limits independently justified and locked at Step 2 before discordance is examined. No second refinement attempt at this stage.
- (v2.5) Audit 1C CLOSED (final report 06_output\Audit1C_final_report_20261030.md). Audit 1D may begin.
- (v2.6, 30 Oct 2026) Editorial consistency corrections only; no change to design, endpoints, definitions or locked decisions: (1) status line updated (Audits 1A–1C closed; no discordance result examined); (2) output folder 06_output added to Section 13 (single output folder, no subfolders; Audit 1A report remains in 01_audits); (3) Audit 1D connectivity scope corrected to 6/18/26 (Sections 4.2, 8, 9) and the approved Audit 1D scope recorded in Section 8; (4) Audit 1C-b and 1C-c marked completed in Section 8; (5) Section 9 lock list completed with parameters previously marked provisional elsewhere (borderline registration-QC category, perturbation-stability criterion, partial-volume handling for PSMA uptake classes, grid-fingerprint precision, perturbation seed, statistical analysis plan); (6) Section 14.1 row for 1C-c moved next to 1C; (7) Limitations updated with the full-cohort Audit 1C-b residual alignment figures in place of the pilot range.

- (v2.7, 30 Oct 2026) Audit 1D executed: execution test (2 patients, run 20261030_214350, TEST_ONLY) then definitive run 20261030_215034 — COMPLETE (100/100 patients, 200/200 scans OK, all outputs and verified manifest). Script revised before execution: v2.6 references, R2018b+, PET–TTB and CT–totseg grid verification, empty-mask handling, OK/PARTIAL/FAILED status, RAS-labelled coordinates, oblique-safe extent, full affine, separate prostate group, script MD5 and input inventory, completion status requiring all mandatory outputs, verified manifest write, execution-test switch.
- (v2.7, 31 Oct 2026) Option D descriptive summary of small components and grid geometry completed (Section 8.2). Descriptive only; no threshold selected.
- (v2.7, 31 Oct 2026) FDG eligibility: author approved two CANDIDATE rules (not locked) — primary ≥ 0.1 mL AND ≥ 3 native voxels; sensitivity ≥ 0.25 mL AND ≥ 3 native voxels (Section 16.3). Correction recorded: the primary candidate excludes 0.24% of pooled FDG volume (inclusive of the 0.18% removed by 0.1 mL alone); per-patient P90 1.7%, maximum 27.8%. Low-burden descriptors added to the Option D script.
- (v2.7, 31 Oct 2026) Provisional methodological framework agreed in discussion (Section 16).
- (v2.8, 31 Oct 2026) Audit 1C-d descriptive review (Section 8.3): residual misalignment is predominantly a per-patient rigid translation (≈ 20 mm in-plane in 21 patients); trunk regions concordant; head/neck separate.
- (v2.8, 31 Oct 2026) Author decision: the earlier decision "no second refinement attempt at this stage" is reconsidered in the light of Audit 1C-d. Amendment 2 (translation-only, closed-form, no optimiser) APPROVED FOR DRAFTING ONLY (Section 15.2): C1–C4 reused unchanged; C3 evaluated on the uncorrected refinement results before any per-patient fallback; inverse-pair verification and one-attempt rule mandatory. No execution and no lock until the author approves the final draft.
- (v2.10, 1 Nov 2026) AMENDMENT 3 LOCKED (specification v2.1): patient-specific regional rigid refinement of head and pelvis by symmetric trimmed ICP on bone surfaces, composed on the adopted transform; selection by independent-of-fit structures with a full-refit bootstrap; Amendment 2 as comparator and fallback; neck deferred; thresholds are engineering choices. Registration investigation reopened by the author for best patient-specific alignment (not a pipeline comparison), motivated by documented regional posture evidence. Implementation and synthetic testing authorised (Section 15.3).
- (v2.10, 1 Nov 2026) Audit 1C-h v3.1 LOCKED (author): three confidence levels; fixed sequence evidence → bone limits → MSD warnings → category → visual QC; V1–V6 as approved; commitment not to change rules or classifications because of counts. Classification script MD5 e4bb1a335ce2e50427acff72677f293a (self-test 22/22); MD5 confirmed on the author's PC and locked before any classification (Section 8.7.21).
- (v2.10, 1 Nov 2026) REGISTRATION INVESTIGATION CLOSED (author). h3 manual landmark assessment abandoned before any reading (qualified independent reading unavailable); automated landmark candidates unvalidated and not to be used; 72 coded CT copies archived, unused. Earlier rejection of Option 1 revised: the soft-tissue cranio-caudal gap is accepted as an acknowledged uncertainty, handled by Audit 1C-h v3 (three confidence levels, conservative primary analysis, stricter sensitivity analysis, limited gross visual QC) locked before classification. Registration confidence and biological discordance are distinct uncertainties (Section 8.7.20).
- (v2.10, 1 Nov 2026) Audit 1C-h3 Stage 1 EXECUTED (run 20261101_104523; script MD5 a11b23f2d39352372282ccbc9172440f): candidates for 100 patients; 31 patients sampled → 62 coded pass-A CTs + 10 repeat CTs, in 06_output only; originals untouched. Blinded reading next (Section 8.7.19).
- (v2.10, 1 Nov 2026) Audit 1C-h3 100-patient candidate preflight (nothing written): LB Def1–Def2 median 5.8 mm (2 patients > 20 mm); LV median 12.2 mm with iliac-vein labels reaching the IVC end in only 18% of CTs (11 patients > 20 mm); train_0024 outlier (about 231–234 mm) retained unchanged; CT z-spacing 1.0–5.0 mm (88/200 at 3.27 mm). Author decisions: L3 ≤ 3.3 mm as an eligibility criterion only; automated LV retained for evaluation; random-sample and disagreement-enriched reliability reported separately (decision: L1/L4 random stratum, L2 all repeats); Stage 1 approved (72 coded CT copies in 06_output only; originals untouched); blinded reading approved. Script MD5 a11b23f2d39352372282ccbc9172440f (self-test 17/17). No residual interpretation, classification or correspondence (Section 8.7.17–8.7.18).
- (v2.10, 1 Nov 2026) Audit 1C-h3 specification v2 approved for script development (H1–H5); script MD5 6afbd3409c16326393522413f6ba2a7d, self-test 16/16. Preflight on 2 patients (nothing written): LB well defined (Def1–Def2 about 5 mm); LV Def1–Def2 13.7–22.4 mm with iliac-vein labels not reaching the IVC end (possible weak definition); CT z-spacing 3.27 mm (above provisional L3 ≤ 3 mm). Before Stage 1 (copies coded CTs; needs author permission): reading load, CT spacing distribution and L3/L1–L4 confirmation; full 100-patient preflight in progress (Section 8.7.15–8.7.16).
- (v2.10, 1 Nov 2026) Audit 1C-h h2z COMPLETE (run 20261101_084243; script MD5 c12377256f9413d2311ae06d0bf4986e; 9/9 checks; Audit 1C-b sign convention PSMA minus FDG-mapped identified and converted, verified by reproduction of t in every patient). Bone cranio-caudal alignment generally reassuring (spine gradient median 0.12 mm/100 mm; leave-pelvis-out median |z| 1.38 mm) with individual pelvic exceptions preserved (train_0078 −5.1, train_0090 −4.2, train_0087 −4.2, train_0069 +4.1 mm). Regional bone estimates are correlated observations from few structures; leave-region-out is internal cross-validation only. Soft-tissue cranio-caudal gap unresolved. No classification; supporting evidence only. Next: landmark-assessment specification (aortic bifurcation, iliac venous confluence) (Section 8.7.12–8.7.14).
- (v2.10, 1 Nov 2026) Audit 1C-h h1b COMPLETE (run 20261101_082311; script MD5 36497eaf36397c4cc5f1467bd8ceef50; checks 3/3): region-restricted MSD tracks in-plane but not cranio-caudal displacement, on an attenuated scale, and agrees only weakly with bone at the adopted state. Author conclusions: MSD useful for detecting substantial in-plane disagreement, not a standalone acceptance metric; cranio-caudal alignment insufficiently assessed in abdomen and pelvis; failed metrics remain demoted; no classification; Option 1 not accepted. Next: Option 2 (per-vertebra residuals; in-sample, not independent validation) (Section 8.7.9–8.7.11).
- (v2.10, 1 Nov 2026) Audit 1C-h h1 EXECUTED (run 20261101_080246, COMPLETE; script MD5 bc44924d1b7f85e12c458c857b1f8e9e; implementation checks 4/4). Bone: 0 violations — retained as decision evidence. Vessel in-plane displacement: P95 deviation 3.36 mm (> 1.0 mm) — demoted to supporting; diagnostics (script MD5 9a68fd325e22bab96de32b5cbc43b094) attribute the failure to vessel orientation (through-plane shift of oblique iliac vessels; ρ 0.68 with |t_z|) and voxel sampling / window effects. Muscle centroid: 5 unexplained violations — remains demoted. Measurement failures, not registration failures. No threshold relaxed; no classification; LIMITED_EVIDENCE not applied by default. Next: critical assessment of region-restricted MSD as acceptance evidence before any new image processing (Section 8.7.5–8.7.8).
- (v2.10, 1 Nov 2026) Audit 1C-h design APPROVED; specification v2 PROVISIONAL (D1–D8 to be checked against exact wording). h0 recorded (Section 8.7.3: A1 + R1 rationale fixed before any category or count; A1 a study-specific engineering criterion, not a physical guarantee). Approved: LIMITED_EVIDENCE category; volume ratio report-only; fixed blinded face-validity sampling (face validity only); entire eligible regional FDG burden to V_I when correspondence is unassessable. h1 script drafting approved. No thresholds locked; no classification, correspondence or discordance analysis.
- (v2.10, 1 Nov 2026) Dataset SUV thresholds extracted read-only (Extract_SUV_thresholds.py MD5 943955e1176b49d6618b7230ad45d88a; run 20261101_073432, COMPLETE, 200/200; CSV MD5 bf7ce5d3174eab228ab992381506655b). Physics note v2 (script MD5 c015d76e3d8f0678834b10e6492df466; actual FDG thresholds; assumptions labelled; level-matched tracer comparison) ACCEPTED PROVISIONALLY as theoretical estimates (Section 8.6). train_0064 FDG threshold 6.02 recorded without inference. Combination of τ terms requires mathematical justification; τ NOT LOCKED. Matching tolerance and registration acceptance kept separate; neither adjusted to maximise assessable patients. No further physics simulations. Next: acceptance-framework specification. No correspondence or discordance analysis.
- (v2.10, 31 Oct 2026) Step 3 physics note DRAFTED (Section 8.6; script MD5 cc08a8288ebb10ee9d58b60781a77ed8): sampling floor about one maximum voxel dimension on the PSMA matching grid (median 4.0 mm, up to 5.5 mm); Gaussian-PSF boundary displacement reported as ranges under explicit background, contrast, size, resolution and threshold assumptions; small lesions most sensitive to registration error. TG-132 Tables 3–4 verified by the author (registration-QA context only). Author decisions: cervical vertebrae primary head/neck QC evidence, brain displacement a separate posture warning, neither automatically establishing acceptance; candidate criteria only after critical review of the note; τ NOT LOCKED. No registration, correspondence or discordance analysis.
- (v2.10, 31 Oct 2026) Audit 1C-g EXECUTED (run 20261031_225909, COMPLETE; 46/46 checks; script MD5 24f205e6b66c35d644794fd4d0489d36 and all six input MD5s confirmed) and accepted; descriptive findings and author decisions in Section 8.5 (head/neck to be evaluated separately with an independent acceptance strategy; no automatic indeterminate classification of head/neck lesions; Step 3 literature review approved; τ not to be fixed at 5, 8 or 10 mm before sampling geometry and threshold-dependent boundary uncertainty are examined separately). No acceptance criteria, classification or discordance analysis.
- (v2.10, 31 Oct 2026) Audit 1C-g LOCK: check-note wording corrected before the definitive run (author-approved; notes now recorded only for failed checks; validation logic, methodology and metrics unchanged; only the two note lines differ from MD5 2fb4baee…, which is superseded). Locked script MD5 24f205e6b66c35d644794fd4d0489d36. Re-verification: --selftest PASSED 17/17; --verify-only PASSED 46/46 (run 20261031_225535; 0 differences in all source reproductions). Verified input MD5s: Audit1C_e_per_patient c66da788a6c73737163b4c09ed76f494; Audit1C_e_structure_level 289b57e69a6f119968b3ae3990db9fd9 (both identical to the Audit 1C-f manifest); Audit1C_f_segments a8087e0ade5acb04ccf3ad0681d6eb3f; Audit1C_f_patient_region_class cbda77b01cb1037051639fbcc422d2e9; Audit1C_f_boundaries_substitution_log 14b882a1ed3c0949c099a26291bc02c5; Audit1C_f_eligibility 83db42bf716e1fc2faa209cff5aaf32b. Author AUTHORISED the single definitive Audit 1C-g run, conditional on the script MD5 and these input MD5s matching on the author's PC. No acceptance criteria, classification or discordance analysis.
- (v2.10, 31 Oct 2026) Audit 1C-g (Step 2; regional consolidation of independent anatomical QC; descriptive) specification APPROVED by the author with D1–D5 as proposed (supplied companion columns; C7 in head/neck; brain separate class HR; physiologically variable organs supporting only; second-largest descriptive only). Safeguards: whole-structure MSD (wsMSD) and region-restricted MSD (rrMSD) never pooled; second-largest only when n ≥ 2, structure count always reported. Script `02_code\Audit1C_g_QC_consolidation.py`, MD5 2fb4baee6c3710efccc113a0fb896b0d. Pre-execution checks: --selftest PASSED (17/17, synthetic, incl. 7 corrupted-input tests); --verify-only on the actual Audit 1C-e/1C-f outputs (run 20261031_225313; Audit 1C-e input MD5s c66da788a6c73737163b4c09ed76f494 and 289b57e69a6f119968b3ae3990db9fd9 identical to those recorded in the Audit 1C-f manifest) PASSED 46/46: fallback set = train_0018, train_0049, train_0067; every structure assigned or excluded; Audit 1C-f patient_region_class (n, median, max, second-largest) and Amendment 2 per-patient set medians (centroid and MSD; supplied and corrected; heldout_bone, heldout_soft, head_neck, fitting) and counts reproduced with 0 differences; 400 rows (100 × 4), no duplicates; no pooled MSD or classification columns. No values written. Definitive run PENDING author review. No acceptance criteria, classification or discordance analysis.
- (v2.10, 31 Oct 2026) Audit 1C-f pelvic overlay follow-up (train_0069, train_0037; script MD5 a31cfc4f32538d03752ae940f4feefa8 [VERIFY vs manifest]): regional rigid misalignment of the pelvis under the adopted transformation (supplied locally better; positioning difference a possible, unproven explanation). Adopted transformations UNCHANGED; no new correction. Independence verified: hips and sacrum are Amendment 2 fitting structures, femora unevaluated, so independent pelvic evidence is vascular/muscular (Audit 1C-f) only. Conditional rule recorded in principle: FDG volume in a region failing the later locked QC criteria contributes to V_I, not V_U; no region classified now. No registration, correspondence or discordance analysis authorised (Section 8.4.11 follow-up).
- (v2.8, 31 Oct 2026) Amendment 2 specifications completed (author): ≥ 5 valid, anatomically distinct fitting structures, otherwise insufficient-landmarks fallback, never omitted from evaluation accounting; C4 procedure and patient number identical to Amendment 1, including largest improvements and deteriorations, CT-only overlays, reviewers blind to tumour correspondence; required outputs defined (15.2.9). Section 15.2 remains DRAFT.
- (v2.9, 31 Oct 2026) Audit 1C-f executed (run 20261031_164711, COMPLETE) and accepted as technically complete; outcome and author decisions in Section 8.4.11 (regional-evaluability limitations documented; no supplementary landmark analysis; fallbacks preserved; prostate non-evaluable; evidence to be combined by region; next: pelvic overlays train_0069/train_0037, then regional acceptance framework).
- (v2.9, 31 Oct 2026) Audit 1C-f LOCK: script MD5 e893c4ff12dd871c52d2c992f8c2b7b3 recorded (author Get-FileHash match). No code changes after this point. Single definitive run authorised.
- (v2.9, 31 Oct 2026) Audit 1C-f synthetic validation (author's PC) and pre-flight passed; train_0059 regional limitation accepted (no extrapolation); script approved (`SPEC_APPROVED = True`); single definitive run authorised once the MD5 is recorded (Section 8.4.10).
- (v2.9, 31 Oct 2026) Audit 1C-f specification decisions (author): fragmentation ≤ 5% for all vessels; descriptive corrected transform for the 3 fallback patients; 10 mm margin provisional with ≥ 100 surface points per direction; vertebral-midpoint boundaries with substitution log; both directional distances restricted to the regional windows. Remains DRAFT.
- (v2.9, 31 Oct 2026) Audit 1C-f specification DRAFTED (Section 8.4; author-approved safeguards: verified regional landmarks with documented substitution; synthetic validation incl. zero displacement, three axes and rotations, cut-face check; separate reporting from Amendment 2's held-out evaluation). Not approved, not executed; no thresholds, tumour data or matching.
- (v2.9, 31 Oct 2026) Amendment 2 executed once (run 20261031_141128, COMPLETE): C1, C2, C3 met (3 patients worse > 2 mm); C4 PASS (author; 6 overlays; no gross trunk misregistration). Author ADOPTED the translation-corrected transformation as the primary trunk registration (97 patients), with supplied-transform fallback for train_0018, train_0049, train_0067 and the supplied transformation as sensitivity analysis. Head/neck alignment not approved (worse > 5 mm in 9 patients; posture) — separate regional QC required. Acceptance-framework discussion and Audit 1C-f direction recorded (Section 16.5; not locked).
- (v2.8, 31 Oct 2026) LOCK: script `02_code\Audit1C_e_translation_refinement.py` MD5 ce038cb661670d673bdc7dbbb448f10b recorded (author Get-FileHash match). No code changes after this point. Single definitive execution authorised; no lesion matching or discordance analysis authorised.
- (v2.8, 31 Oct 2026) Pre-execution clarifications after script review (author; before the MD5 lock; C1–C4 and estimator unchanged): (1) Amendment 1's actual C4 rule, read from its code, is "5 patients with the largest refinement change"; Section 15.2.5 now specifies this as the primary C4 set, with the largest held-out-bone deterioration and improvement added as a prospective safety extension (not an Amendment 1 replication); C4 judged on all overlays. (2) Equivalence of the Amendment 1 change measure with |t| verified analytically and checked numerically in every patient (15.2.3 (d)). (3) Dice, MSD/HD95 and the C4 overlay function copied verbatim from Amendment 1 after the pre-flight showed small MSD/HD95 differences (cause: missing 3-voxel margin in the region of interest, not a different distance definition).
- (v2.8, 31 Oct 2026) Final clarifications (author): C4 specified as 5 patients with Amendment 1's selection procedure; any inverse-verification failure stops the experiment (implementation error; never a fallback); head/neck deterioration reported and handled by regional QC, gross misregistration at C4 remains a safety finding. Author APPROVED the Amendment 2 methodology and LOCKED Section 15.2; single execution authorised once the script MD5 is recorded. C1–C4 and estimator unchanged. No correspondence or discordance analysis authorised. NOT LOCKED; to be formalised in a dated amendment and locked at Step 2 before any correspondence or discordance analysis. Audit 1D remains open until its written report is completed.

## 13. Analysis workflow and operating rules (v2.2, final; folder structure updated v2.6)
Working folder: `D:\DEEP PSMA\PAPER3_Discordance` (00_protocol, 01_audits, 02_code, 03_intermediate, 04_results, 05_logs, 06_output).
- 00_protocol: protocol versions and addenda.
- 01_audits: Audit 1A report (Audit1A_output.md).
- 02_code: all scripts.
- 06_output: single output folder (no subfolders) for all audit and analysis outputs from Audit 1B onward.
Original data in `D:\DEEP PSMA\DEEP PSMA EXTRACTED` is read-only and never modified. Paper 3 uses only the original dataset; no outputs from earlier projects are read or reused.

| Step | Task |
|---|---|
| 1A | Cohort eligibility audit |
| 1B | TotalSegmentator labels and anatomical coverage audit |
| 1C | Transformation direction and registration audit (develops/assesses QC rules; no discordance outcomes examined) |
| 1D | FDG/PSMA tumour-mask characteristics audit |
| 2 | Lock all methodological parameters AND the statistical analysis plan in a date-stamped protocol addendum (optional OSF) |
| 3 | Extract FDG and PSMA lesions; calculate native-space quantitative parameters (volume, SUVmean, TLG / TL-PSMA) |
| 4A | Case- and region-level registration QC (application of the rules locked at Step 2) |
| 4B | Transform FDG lesion masks into PSMA reference space (verified ITK convention) |
| 4C | Technical verification of the pipeline using selected CT and mask overlays (software check, not lesion adjudication) |
| 4D | Assign lesion-level indeterminate status (coverage measures, Section 8.1): site-indeterminate (CT/anatomical coverage) and correspondence-indeterminate (PSMA PET field of view after transformation, registration QC, resampling loss) |
| 5 | Distance-based lesion matching (locked reference tolerance); record one-to-many correspondences and anatomical-consistency flags |
| 6 | Calculate primary endpoint DMTVF_ref |
| 7 | Quantify indeterminate burden, primary uncertainty interval and sensitivity envelope |
| 8 | Registration perturbation analysis (locked seed, magnitudes, number) |
| 9 | PSMA uptake at unmatched FDG lesion locations (native PSMA PET; absent / low / positive) |
| 10 | Anatomical distribution of discordant burden |
| 11 | Supplementary sensitivity analyses (connectivity, tolerances, overlap vs distance, anatomical consistency, whole-body tracer burden relationships) |
| 12 | Execute the locked statistical analysis plan; generate tables and figures |
| 13 | Prepare manuscript and reproducibility materials |

### Mandatory operating rules
1. **Documentation:** every audit (1A–1D) has a completed, saved report before the next audit begins.
2. **Reproducibility:** software versions, parameters, random seeds, input file lists and output logs recorded throughout (05_logs), not only at manuscript preparation.
3. **Analysis integrity:** numerical thresholds and the statistical analysis plan are locked before discordance results are examined; any later change is documented as a protocol deviation with date and reason.

## 14. Software architecture and data-exchange rules (v2.3)

### 14.1 Division of work
| Step | Software |
|---|---|
| 1A Cohort eligibility audit | MATLAB |
| 1B TotalSegmentator label audit | MATLAB |
| 1C ITK transformation audit | Python / SimpleITK |
| 1C-c Bone-based rigid refinement (Amendment 1; not adopted) | Python / SimpleITK |
| 1C-d Descriptive regional QC review | Python (pandas) |
| 1C-e Translation-only refinement (Amendment 2; executed, ADOPTED) | Python / SimpleITK |
| 1C-f Independent vascular and muscular registration QC (Section 8.4; executed 20261031_164711, COMPLETE) | Python / SimpleITK |
| 1D Tumour-mask audit | MATLAB |
| 1D Option D descriptive summary | Python (pandas) |
| 2 Protocol and SAP lock | Documentation |
| 3 Lesion extraction and quantitative measurements | MATLAB |
| 4A Case/region registration QC | Python + MATLAB |
| 4B FDG → PSMA transformation | Python / SimpleITK |
| 4C Technical overlay verification | Python / SimpleITK |
| 4D Lesion-level indeterminate classification | MATLAB |
| 5 Distance-based correspondence | MATLAB |
| 6–7 DMTVF and uncertainty intervals | MATLAB |
| 8 Registration perturbations | Python / SimpleITK (MATLAB summaries) |
| 9 PSMA uptake classes | MATLAB |
| 10 Anatomical distribution | MATLAB |
| 11–12 Sensitivity analyses and statistics | MATLAB |
| 13 Manuscript outputs | MATLAB + documentation |
Principle: Python/SimpleITK handles ALL spatial transformations, orientation (origin, direction cosines, LPS/RAS), registration QC geometry, perturbations, resampling integrity and grid metadata. MATLAB handles lesion quantification, correspondence, discordance burden, statistics and manuscript outputs.

### 14.2 Coordinate safeguards
- SimpleITK uses LPS physical coordinates; MATLAB/NIfTI conventions vary by function and version → no axis convention assumed; orientation verified numerically with known landmarks (Audits 1C, step 4C) in Python only.
- Common-grid strategy: Python resamples FDG lesion labels onto the NATIVE PSMA PET grid (nearest-neighbour). PSMA PET and PSMA TTB unchanged. MATLAB works only on arrays sharing this grid and needs only the voxel spacings.
- Native FDG PET and masks remain the sole source of FDG volume, SUVmean and TLG.

### 14.3 Lesion preservation during resampling
- Unique FDG lesion IDs assigned in native FDG space BEFORE transformation; the integer-labelled lesion image is resampled (not separate binary masks).
- Per lesion recorded: native voxel count and volume; transformed voxel count; % volume change; disappearance or substantial shrinkage; potential collisions with other transformed lesions.
- A single resampled label image cannot reveal every collision → a separate diagnostic resampling check is used when necessary to identify competing labels.
- Disappeared lesions → correspondence-indeterminate, never discordant. Minimum-volume eligibility applied in native FDG space.

### 14.4 Physical distances
- Grid orthogonality (unskewed axes) verified from PSMA image metadata; only then distance (mm) = voxel index × spacing per axis.
- Metric defined explicitly: minimum Euclidean distance between boundary-voxel centres (in mm) of an FDG lesion and a PSMA lesion on the common grid — an approximation of continuous surface-to-surface distance.
- Implementation: surface point clouds with knnsearch (no large pairwise matrices); bwdist (isotropic voxel units) not used for physical distances.

### 14.5 Data exchange and manifests
- Python ↔ MATLAB exchange only via NIfTI images and structured CSV tables.
- Every processing step writes a JSON manifest:
| Field | Content |
|---|---|
| Case ID | train_XXXX |
| Input/output files | Exact paths |
| Reference grid | Size, spacing, origin, direction |
| Grid fingerprint | Deterministic hash of size + spacing + origin + direction; spatial values rounded to 1e-4 mm (provisional), direction cosines to an explicitly defined precision [LOCK-AFTER-AUDIT] |
| Transformation | File, mapping direction, inversion status |
| Interpolation | Nearest-neighbour for labels |
| Lesion IDs | Original and transformed identifiers |
| QC | Coverage, missing, shrunk or merged components, collisions |
| Software | MATLAB, Python, SimpleITK versions |
| Reproducibility | Parameters, random seeds, timestamps |
- MATLAB verifies identical array dimensions AND the expected grid fingerprint before every correspondence operation (identical dimensions alone do not guarantee identical physical grids).

## 15. Protocol amendments

### Amendment 1 — Bone-based rigid refinement of the supplied PSMA–FDG registration (29 October 2026)
Recorded BEFORE any refinement result is available. No tumour correspondence, DMTVF or discordance result has been examined.

**Rationale (Audit 1C-b, anatomy only).** With the supplied transform, whole-body bony centroid distance had a median of 7.3 mm (P90 20.2 mm); about 20 patients showed a near-uniform 17–23 mm, mainly in-plane, bony offset. An optimal rigid re-fit of bony centroids left a median residual of 2.8 mm (P90 5.0 mm; in-sample, optimistic), suggesting that most residual error is a global rigid error. The cervical spine showed additional posture-related, position-dependent differences.

**Method (pre-specified; Python/SimpleITK; script `02_code\Audit1C_c_rigid_refinement.py`).**
- Fixed image: PSMA CT. Moving image: FDG CT. Transform: Euler3D, initialised with the supplied PSMA/rigid.tfm (SimpleITK resampling convention: maps PSMA points to FDG points).
- Metric: Mattes mutual information (50 bins), restricted to a fixed-image mask = trunk-bone fitting structures on the PSMA CT dilated by 5 mm; random sampling 10% with a fixed seed; linear interpolation.
- Fitting structures: vertebrae T1–L5, S1, sacrum, hips — only those apparently complete according to automated coverage criteria on both CTs, excluding Audit 1B flagged levels (+ adjacent) and levels adjacent to fragmented vertebrae.
- Optimiser: regular-step gradient descent, physical-shift scaling; three resolution levels (shrink 8/4/2; smoothing 4/2/1 mm). Settings are fixed in advance and not tuned on evaluation results.
- Output per patient: refined resampling transform (PSMA → FDG) and its inverse for point mapping (FDG → PSMA), stored in 06_output.

**Evaluation (all 100 patients; supplied vs refined; same structures under both).**
- Held-out bone (primary): ribs 1–12 left/right, sternum, clavicles, scapulae (not used for fitting).
- Held-out soft tissue: liver, spleen, kidneys, aorta, heart, lung lobes, urinary bladder.
- Fitting structures: reported separately (in-sample).
- Head/neck (C1–C7, brain): reported separately; not used for fitting.
- Metrics: centroid distance (point mapping), Dice, mean surface distance, HD95; patient-level medians; paired differences (refined − supplied); proportion of patients improved/worsened; proportion exceeding candidate limits (5/10/15/20 mm, descriptive).
- Coverage (Section 8.1) recomputed with the refined transform.

**Pre-specified adoption criteria (anatomical QC only).**
- C1: lower cohort median of the patient-level held-out-bone median centroid distance with refinement, and median paired difference < 0 mm.
- C2: held-out soft tissue not meaningfully worse (median paired difference ≤ +0.5 mm).
- C3 (patient-level safeguard): held-out-bone median centroid distance worsens by > 2 mm in ≤ 5 patients. For any patient with such worsening, or a failed/non-converged registration, the supplied transform is retained for that patient (per-patient fallback based solely on anatomical QC).
- C4: visual verification of overlays for the patients with the largest changes shows no gross misregistration.
- If C1–C4 are met, the refined transform (with per-patient fallback) becomes the PRIMARY transformation and the supplied transform is used as a sensitivity analysis; otherwise the supplied transform remains primary. The final decision is recorded by the author before Audit 1D and before any discordance step.
- Head/neck alignment is not addressed by the global refinement; regional handling (QC limits / indeterminate status) is decided at Step 2.

**LOCKED (author approval, 29 October 2026, before execution).** Criteria C1–C4, the fixed fitting parameters, the structure sets and the per-patient fallback rule above are locked. Locked script: `02_code\Audit1C_c_rigid_refinement.py` (MD5 0c350ab57b59c1dea8fc729b434457b7). No parameter may be changed after results are seen; any deviation requires a new dated amendment.

**Decision sequence (locked):**
1. Lock C1–C4 and the fitting parameters (done).
2. Execute the same refinement procedure for all 100 patients.
3. Evaluate held-out bone and soft-tissue alignment.
4. Apply the predefined per-patient fallback rules.
5. Record the adoption decision and the final transformation (refined or supplied) for every patient.
6. Close Audit 1C before starting Audit 1D.

**Interpretation boundary.** Improved CT–CT anatomical alignment does not prove more accurate PSMA–FDG lesion matching. The refinement is described as an *anatomically evaluated registration refinement*, not as validated lesion correspondence.

**OUTCOME (30 October 2026; locked run 20261029_212531; report 06_output\Audit1C_final_report_20261030.md).**
- C1 met: held-out bone median 8.0 → 6.5 mm; median paired difference −0.59 mm.
- C2 met: held-out soft median paired difference −0.65 mm.
- C3 NOT met: 22 patients worsened by > 2 mm on held-out bone (per-patient fallback applied).
- C4 NOT met: gross misregistration on overlay review (e.g., train_0070, vertebral-level displacement; train_0066, head misalignment).
- **Decision (author): refinement not adopted; supplied rigid transformation remains PRIMARY for all 100 patients.** Refined transforms retained in 06_output for technical documentation only. No second refinement attempt at this stage. The experiment is kept in the technical documentation and is not a major manuscript result unless it becomes relevant to interpreting discordance uncertainty.

### 15.2 Amendment 2 — Translation-only refinement of the supplied PSMA–FDG registration (DRAFT, 31 October 2026)
**Status: LOCKED (text) — author approval 31 October 2026, after final clarifications; single execution AUTHORISED once the script is written and its MD5 recorded (15.2.10). No refinement result exists.** Recorded before any Amendment 2 result is available. No tumour correspondence, DMTVF or discordance result has been examined.

**15.2.1 Rationale (anatomy only).** Amendment 1 (mutual-information rigid optimisation) was not adopted: C3 failed (22 patients worse > 2 mm on held-out bone) and C4 failed (gross misregistration, e.g., vertebral level slip in train_0070). Audit 1C-d (Section 8.3) showed that the residual error of the supplied transformation is predominantly a per-patient rigid translation (≈ 20 mm in-plane in 21 patients; |median vector| / median distance 0.98), with trunk regions moving together. Amendment 2 tests one deterministic correction aimed at that error, without rotational freedom and without an optimiser. Amendment 2 is a subsequent methodological experiment motivated by anatomical QC findings; Amendment 1 and its failure are reported in full.

**15.2.2 Estimator (fixed in advance).**
- Displacement per structure, as defined in Audit 1C-b: d_s = FDG/rigid.tfm(c_FDG,s) − c_PSMA,s (LPS physical coordinates, mm), where c are TotalSegmentator structure centroids on the native CTs.
- Fitting structures (identical to Amendment 1): vertebrae T1–L5, S1, sacrum, left and right hip — only those apparently complete according to automated coverage criteria on both CTs, excluding Audit 1B flagged levels and their adjacent levels, and levels adjacent to fragmented vertebrae. Assignment fixed per patient before execution; landmarks are never selected according to performance.
- Translation t = component-wise median of d_s over the fitting structures (each axis separately).
- Minimum: ≥ 5 valid, anatomically distinct fitting structures (each vertebral level, the sacrum, and each hip count once). With fewer than 5, no correction is computed and the supplied transformation is retained for that patient, recorded as an **insufficient-landmarks fallback** — not as a successful refinement.
- Computed once per patient; no iteration, no parameter tuning, no re-runs.

**15.2.3 Corrected transformation (definition and verification).**
- Point mapping (FDG → PSMA): T′(p) = FDG/rigid.tfm(p) − t.
- Resampling (PSMA → FDG, `sitk.Resample` convention): R′(q) = PSMA/rigid.tfm(q + t).
- Then R′(T′(p)) = PSMA/rigid.tfm(FDG/rigid.tfm(p)) = p, because the supplied files are an exact inverse pair (Audit 1C-a).
- Mandatory numerical verification BEFORE any evaluation, every patient: (a) R′(T′(p)) = p and T′(R′(q)) = q for test points spanning the field of view, maximum error < 10⁻⁶ mm; (b) the fitting-structure displacements under T′ equal d_s − t to < 10⁻⁶ mm (confirms sign and composition order); (c) T′(p) = T(p) − t and R′(q) = R(q + t) at the test points; (d) the Amendment 1 change measure equals |t| at every fitting centroid. Any verification failure, in any patient, STOPS the whole experiment for investigation. It is an implementation error, never a refinement result and never a per-patient fallback. A corrected implementation may only be run after the error and the code change are documented as a protocol deviation (date, reason, new script MD5); C1–C4, the estimator and the structure sets may not change.
- Implementation in Python/SimpleITK (script `02_code\Audit1C_e_translation_refinement.py`, written only after approval of this draft; MD5 recorded at lock).

**15.2.4 Evaluation (all 100 patients; supplied vs corrected; identical structure sets to Amendment 1).**
- Held-out bone (primary): ribs 1–12 left/right, sternum, clavicles, scapulae.
- Held-out soft tissue: liver, spleen, kidneys, aorta, heart, lung lobes, urinary bladder.
- Fitting structures: reported separately (in-sample; not evidence of improvement).
- Head/neck (C1–C7, brain): reported separately; not used for fitting; not expected to be corrected by a trunk translation.
- Metrics as Amendment 1: centroid distance (point mapping), Dice, mean surface distance, HD95; patient-level medians; paired differences (corrected − supplied); patients improved/worsened; proportions above 5/10/15/20 mm (descriptive). Coverage (Section 8.1) recomputed with the corrected transformation.

**15.2.5 Adoption criteria — C1–C4 of Amendment 1, unchanged.**
- C1: lower cohort median of the patient-level held-out-bone median centroid distance, and median paired difference < 0 mm.
- C2: held-out soft tissue not meaningfully worse (median paired difference ≤ +0.5 mm).
- C3: held-out-bone median centroid distance worsens by > 2 mm in ≤ 5 patients.
- C4: visual verification of overlays shows no gross misregistration. Overlay type exactly as Amendment 1 (`overlay()` copied verbatim: PSMA CT, coronal and sagittal, PSMA labels green, resampled FDG labels magenta; CT and TotalSegmentator labels only, no tumour masks; reviewers blind to tumour correspondence).
  - **Primary C4 set (replicates Amendment 1's selection principle):** the 5 patients with the largest refinement change, computed exactly as Amendment 1's `refinement_change_median_mm` — median, over the patient's fitting-structure centroids in PSMA space, of the displacement between the supplied and the corrected resampling transforms. For a rigid supplied transform this equals |t| at every point (|R(q + t) − R(q)| = |Mt| = |t|, M orthogonal); the script computes the measure explicitly and verifies this equality (15.2.3 (d)).
  - **Prospective safety extension (NOT part of the Amendment 1 replication):** additionally, the patient with the largest held-out-bone deterioration and the patient with the largest held-out-bone improvement, if not already in the primary set (5–7 overlays in total).
  - **C4 verdict:** judged on ALL selected overlays; gross misregistration in any of them fails C4.
- **C1–C4 are evaluated on the corrected results for all patients BEFORE any per-patient fallback is applied.** The fallback cannot be used to satisfy C3 or any other criterion.
- Insufficient-landmarks patients (< 5 fitting structures) are never omitted: they appear in every evaluation and fallback table with their status and reason. C1–C4 are reported (i) over all patients with a computed correction and (ii) over all 100 patients, insufficient-landmarks patients entering (ii) with the supplied transformation (paired difference 0). Adoption requires C1–C4 to be met in (i); results in (ii) and the number of insufficient-landmarks patients are reported alongside, so that their exclusion from (i) cannot artificially improve any criterion.

**15.2.6 Decision and fallback.**
- If C1–C4 are all met: the corrected transformation becomes PRIMARY; per-patient fallback to the supplied transformation where held-out bone worsens by > 2 mm or the estimator has too few fitting structures (insufficient-landmarks fallback). Verification failure is not a fallback reason (15.2.3); the supplied transformation is retained as a sensitivity analysis.
- If any criterion fails: the supplied transformation remains PRIMARY for all patients; conservative regional QC applies (Section 7); no further refinement attempt.
- The decision is recorded by the author before any QC limit, tolerance or correspondence step.

**15.2.7 Head/neck.** The trunk translation is applied to the whole image, but cervical/head residuals are reported separately. Any head/neck deterioration (paired differences, patients worsened) is reported in full and addressed through separate head/neck regional QC at Step 2; failure to improve head/neck alignment is not by itself an Amendment 2 failure. However, gross misregistration identified during C4 visual review, in the head/neck or elsewhere, remains a safety finding: it counts against C4.

**15.2.8 One-attempt rule and interpretation boundary.**
- One execution of the locked script. No changes to estimator, structure sets, criteria or code after results are seen; any deviation requires a new dated amendment and is reported as such.
- If the residual is predominantly a rigid translation, held-out structures should improve; failure would indicate that the residual is not adequately described by a translation. Either outcome is reported.
- Improved CT–CT anatomical alignment does not establish accurate PSMA–FDG lesion correspondence; the correction is described as an anatomically evaluated translation refinement.
- QC limits and the reference matching tolerance are set only after the Amendment 2 decision, justified on held-out (not in-sample) residuals of the final transformation.

**15.2.9 Required outputs (06_output, single folder; run timestamp in every file name).**
- Per-patient table: translation vector t (x, y, z, magnitude, in-plane), number and list of fitting structures, correction status (computed / insufficient-landmarks fallback / verification failure), fallback reason.
- Fitting-structure residuals (in-sample) and held-out bone and soft-tissue residuals, supplied vs corrected, per structure and per patient; head/neck reported separately.
- C1–C4 results, computed (i) over patients with a computed correction and (ii) over all 100 patients, with the number of insufficient-landmarks patients.
- Inverse-transform verification log per patient (maximum round-trip errors, displacement-sign check).
- C4 QC overlays for the selected patients, with the selection list (primary / safety extension, change measure, translation magnitude, held-out-bone difference, verdict field for the author).
- Coverage (Section 8.1) recomputed with the corrected transformation.
- Reproducibility manifest: parameters, software versions (Python, SimpleITK, numpy, pandas), input file inventory with checksums, script MD5, timestamps, completion status (all patients processed and all mandatory outputs written, as for Audit 1D).

**15.2.10 Status and lock record.**
- Section 15.2 text LOCKED by the author on 31 October 2026 (C1–C4 and estimator unchanged from the approved draft; pre-execution clarifications recorded in Section 12).
- Implementation check (pre-flight, `--preflight 3`, train_0001–train_0003; supplied-transform metrics and verification only, no Amendment 2 result computed):
  - Pre-flight 1: verification (a)–(c) passed (≤ 5 × 10⁻¹³ mm); centroid distances and Dice reproduced Amendment 1 exactly (205/205); MSD/HD95 small differences (max 0.26 / 1.44 mm, urinary bladder train_0001), traced to a missing 3-voxel region-of-interest margin; Amendment 1's `surface_metrics()`, Dice filter and `overlay()` then copied verbatim.
  - Pre-flight 2 (final code): verification (a)–(d) passed in all 3 patients (≤ 3 × 10⁻¹³ mm; change measure = |t| to ≤ 3 × 10⁻¹³ mm); centroid distance, Dice, MSD and HD95 reproduce Amendment 1 with 0 difference in 205/205 structures; held-out measurements complete; supplied coverage identical to Audit 1C-b.
- C4 procedure verified against Amendment 1's code (primary: 5 largest refinement changes; safety extension: largest held-out-bone deterioration and improvement); `C4_PROCEDURE_VERIFIED = True`.
- **Locked script:** `02_code\Audit1C_e_translation_refinement.py`, MD5 **ce038cb661670d673bdc7dbbb448f10b** (40,868 bytes; confirmed by the author's Get-FileHash on the copy in 02_code, 31 October 2026, matching the delivered file).
- **No code changes after the MD5 is recorded.** Any change after that point is a protocol deviation (date, reason, new MD5) and requires a new author decision.
- **LOCK COMPLETE (31 October 2026).** Author approval: the single definitive Audit 1C-e execution is AUTHORISED (run without `--preflight`). No lesion matching, correspondence or discordance analysis is authorised.

**15.2.11 OUTCOME (31 October 2026; run 20261031_141128; status COMPLETE).**
- Execution: 100/100 patients; correction computed in 100 (0 insufficient-landmarks fallbacks, 0 failures); verification (a)–(d) passed; supplied centroid distance, Dice, MSD and HD95 reproduced Amendment 1 with 0 difference in 6,687/6,687 structures; held-out evaluation complete (criteria valid). Translation |t|: median 6.0 mm (IQR 3.2–13.1; max 24.5); fitting structures 11–21 per patient.
- Criteria (evaluated before any fallback; populations (i) and (ii) identical):
  | Criterion | Result | Met |
  |---|---|---|
  | C1 held-out bone | median 8.0 → 5.2 mm; median paired difference −1.76 mm | YES |
  | C2 held-out soft tissue | median paired difference −1.22 mm (≤ +0.5) | YES |
  | C3 held-out bone worse > 2 mm | 3 patients (≤ 5): train_0018 (+3.6), train_0049 (+3.0), train_0067 (+3.7) | YES |
  | C4 visual review (author) | 6 overlays (primary train_0012, train_0086, train_0066, train_0099, train_0021; safety extension train_0067): no vertebral-level slip, rotation/flip or gross trunk misregistration | YES (trunk) |
- Further descriptive results: held-out bone > 10 mm in 5 patients (36 supplied), > 15 mm in 0 (23 supplied); the 23 patients > 15 mm improved from median 20.1 to 6.5 mm; held-out Dice 0.32 → 0.55; held-out soft tissue worse > 2 mm in 6 patients (not a criterion); fitting structures 6.8 → 2.6 mm (in-sample). 26/6,687 structures (thin ribs/vertebrae/clavicle) had a corrected-transform resampled volume ratio outside 0.95–1.05 (all within range under the supplied transform); centroid-based C1–C3 unaffected (native centroids mapped).
- C4 per-patient notes (author): train_0086 PASS (minor head residual); train_0099 PASS trunk (head positional mismatch); train_0012 PASS trunk, head/neck clearly worsened (anterior displacement; head/neck 15.0 → 28.6 mm); train_0021 PASS trunk (mild head deterioration); train_0066 PASS trunk (head posture mismatch under both transforms); train_0067 PASS (supplied slightly better on trunk, consistent with its fallback; soft-tissue offset; no major new head deterioration). Limitation: selected coronal/sagittal views, not a full 3-D assessment.
- Head/neck (reported per 15.2.7): median paired difference +0.1 mm; worse > 2 mm in 28 patients, > 5 mm in 9 (max train_0012 +13.6 mm). Mechanism: differing head/neck posture between examinations; a trunk translation cannot correct it.
- **DECISION (author, 31 October 2026): Amendment 2 ADOPTED.** The translation-corrected transformation (point mapping T′(p) = FDG/rigid.tfm(p) − t; resampling R′(q) = PSMA/rigid.tfm(q + t)) is the PRIMARY trunk registration for 97 patients; the supplied transformation is retained for train_0018, train_0049 and train_0067 (held-out bone worse > 2 mm) and serves as a sensitivity analysis for all patients. Per-patient t and status: 06_output\Audit1C_e_per_patient_20261031_141128.csv.
- **Qualification:** this decision does NOT establish satisfactory head/neck alignment. Head/neck correspondence requires independent regional QC before lesion correspondence; regions failing it are correspondence-indeterminate. Anatomical CT–CT improvement does not establish accurate PSMA–FDG lesion correspondence.

### 15.3 Amendment 3 — Patient-specific regional rigid refinement of head and pelvis (LOCKED specification, 1 November 2026)
**15.3.1 Rationale.** The registration investigation was reopened by the author (1 Nov 2026) with a stated objective: the best anatomically supported transformation for each patient and region of this dataset (not a comparison of registration pipelines), because PSMA–FDG discordance can be stated confidently only under the best achievable alignment. Documented evidence before and independent of the Audit 1C-h counts motivates regional correction: head/neck posture (Amendment 2 C4; head/neck worse > 5 mm in 9 patients), pelvic rigid offsets (Audit 1C-f; train_0069/train_0037 overlays) and pelvic cranio-caudal offsets (h2z). Amendment 2 corrected only a whole-body translation; no rotational refinement has been adopted (Amendment 1, which included rotation, was not adopted). Tumour-blind throughout.

**15.3.2 Locked specification** (`Amendment3_regional_rigid_refinement_spec_LOCKED_v2_1_20261101.md`). Guiding principle (author): each patient receives the anatomically best-supported transformation, not the transformation that performs best on average; added complexity only where anatomical improvement is demonstrated.
- Targets: head and pelvis; neck (C3–C7) deferred to a separate, cautious investigation; thorax and abdomen stay on Amendment 2.
- Composition on the current adopted transform T_A (Amendment 2; supplied for the 3 fallback patients); exact inverse composition verified (< 10⁻⁶ mm).
- R1: symmetric trimmed ICP (80% per direction; equal total weight per direction; weighted Kabsch; < 0.01 mm / 0.01° or 100 iterations); one-way solution as a bias diagnostic (> 1° or > 1 mm flagged → visual review).
- R2: common field of view (10 mm margin); skull ≥ 2,000 points and ≥ 60 mm span; pelvis ≥ 2,000 points per hip, both sides; otherwise not estimable → Amendment 2.
- R3: caps 10° / 20 mm (experimental safeguards; reject, not clip).
- R4: adopt only if the selection-structure centroid improvement exceeds max(1.0 mm engineering minimum, 2 × bootstrap SD); the bootstrap (200 replicates) repeats the full fit on resampled fitting-surface points and recomputes the selection improvement (fitting uncertainty, not total measurement noise); selection surface distance not worse by > 0.5 mm; blinded visual veto only.
- R5: head — fit complete eligible skull / select brain / skull spatial cross-validation as stability diagnostic / C1–C2 supporting only; pelvis — fit hip bones only (labels 77/78; no femur, sacrum or S1 voxels: coded check) / select sacrum + S1 / report femoral heads (top 40 mm), supporting vessel/muscle surface distance and overlays. Coupled structures give internal evidence, not independent accuracy.
- R6: fixed anatomical assignment (pelvis below L5/S1; head above C2/C3; Amendment 2 elsewhere); 15 mm boundary band evaluated under both transforms at the correspondence stage → V_I if the category differs.
- R7: patient-specific safety monitor — among adopting patients with evaluable reporting structures, > 10% worse by > 2 mm suspends further adoption in that region pending review; no automatic replacement of selected transforms.
- Thresholds and caps are engineering choices, not clinically validated standards. No patient-specific result may change these rules without a documented, dated amendment. Amendment 2 remains the comparator and fallback; the original locked Audit 1C-h results are retained and reported alongside.
- The validity amendment of the Audit 1C-h classification rule (thorax second-largest statistic over about 25 structures; rib z; clavicle/scapula posture) remains a separate pending task.

**15.3.3 Implementation (numerical stage; not yet executed on study data).** Script `02_code\\Amendment3_regional_rigid.py`, MD5 d8cc08bfff82d8f7569f1469d91190f9; synthetic self-test 12/12 (known rotation/translation recovered; robustness to 20% corrupted points and unequal partial coverage; equal direction weighting verified by duplicating one point set; caps reject; full-refit bootstrap with 200 replicates; no adoption or bootstrap when already aligned; label hygiene). The visual-veto stage is separate; adoption decisions from this stage are PROVISIONAL. Next: author MD5 confirmation, self-test and 3-patient pre-flight (implementation checks only, nothing written). Author confirmed MD5 (Get-FileHash) and self-test 12/12 on the study PC.

**15.3.4 Pre-flight result (1 November 2026; train_0001–train_0003; implementation checks only, nothing written).** Runtime about 60–115 s per patient. Pelvis: all converged (28–36 iterations), 1.08–3.53°, 3.0–8.45 mm, within caps; no one-way flags. Head: train_0001 converged, 5.11°, 34.34 mm → exceeds the 20 mm cap; train_0002 not converged at 100 iterations (9.90°, 20.07 mm); train_0003 2.12°, 2.73 mm, within caps (FDG skull surface 48k points vs PSMA 107k).

**15.3.5 Amendment 3a — head-specific parameter change (dated 1 November 2026; PROVISIONAL).**
- **Disclosure:** these changes were made **after the pre-flight results of three patients (train_0001–train_0003) had been examined.** They are therefore informed by study data and are reported as such (Principle 7).
- **Rationale:** the original head caps (10°, 20 mm) were set before any head data were seen and are small relative to the head-posture differences already documented under Amendment 2 (brain centroid displacement median 15 mm, P90 45 mm, maximum 117 mm). One of three pre-flight patients was rejected by the cap and one did not converge within 100 iterations.
- **Changes:**

| Parameter | Amendment 3 (locked v2.1) | Amendment 3a |
|---|---|---|
| Head translation cap | 20 mm | **50 mm** (provisional) |
| Head rotation cap | 10° | **20°** (provisional) |
| Maximum ICP iterations (both regions) | 100 | **300**; convergence remains mandatory |
| Pelvic caps | 10°, 20 mm | unchanged |
| Brain selection rule, bootstrap, one-way diagnostic | — | unchanged |
| Visual veto | — | retained |

- **Caveats (author):** more iterations improve the opportunity to converge; they do not guarantee a correct registration. Larger caps increase the risk that ICP settles in a local minimum matching unrelated surfaces. Caps remain engineering safeguards, not validated limits.
- **Validation before any full run:** skull (and pelvic) overlays of train_0001–train_0003, comparing Amendment 2 with the proposed correction, are inspected by the author to confirm anatomical alignment rather than matching of unrelated surfaces. The full Amendment 3a run is **not authorised** until the author has validated these overlays.
- **Implementation (draft, pending author MD5 check):** `Amendment3_regional_rigid.py` modified: region-specific caps (`CAPS = {"head": (20°, 50 mm), "pelvis": (10°, 20 mm)}`), `MAX_ITER = 300`, new `--overlay <case IDs>` mode (diagnostic PNGs, left panel labelled Amendment 2, right panel Amendment 3a candidate with convergence and cap status; no bootstrap, no tables, no decisions). Two checks per patient: (i) transform composition verified numerically; (ii) image-level check: the centroid of each displayed FDG label in the resampled image must agree with the point-mapped centroid of the same FDG voxels (within the PSMA grid) to within the largest PSMA voxel dimension, for both Amendment 2 and the candidate; PASS/FAIL printed and shown on each figure. Background = PSMA CT (reference, not transformed). Script modification and overlay generation explicitly authorised by the author (1 Nov 2026); full-cohort run and adoption decisions not authorised. Self-test 12/12 unchanged. MD5 b21a405145526d56a59bb5950ef53acf (the locked Amendment 3 version, MD5 d8cc08bfff82d8f7569f1469d91190f9, remains the reproducibility reference).

**15.3.6 Overlay run and author's visual assessment (run 20261101_225414; 2 November 2026).** Self-test 12/12; image-level transform check PASS in all 6 patient × region overlays (worst centroid discrepancy 0.73–1.20 mm; tolerance 2.50–3.27 mm). Fits (no bootstrap): head train_0001 5.11°, 34.34 mm (89 iterations); train_0002 9.68°, 19.80 mm (80); train_0003 2.48°, 2.89 mm (58); pelvis train_0001 2.16°, 8.45 mm; train_0002 1.07°, 2.99 mm; train_0003 3.46°, 3.75 mm; all converged, within the Amendment 3a caps, no one-way flags.

*Point-sampling observation:* train_0002 head did not converge in the pre-flight (100 iterations) but converged here (80 iterations). train_0003 values differed slightly (≤ 0.4°, ≤ 0.2 mm). The cause is that the pre-flight bootstrap for train_0001 consumed random numbers, so later patients received different surface-point subsamples (≤ 4,000 points). Results with a fixed seed are reproducible; their sensitivity to subsampling is a stability issue. Author decision: implementation unchanged (Option A); the 200-refit bootstrap is the stability measure, noting that it does not guarantee convergence or detect systematic segmentation error.

*Author's visual assessment.* These are **preliminary anatomical plausibility checks, not independent validation of registration accuracy.** They are based on selected 2D views (coronal, sagittal, axial) of segmentation overlays, not a full examination of every anatomical level.

| Patient | Region | Visual finding (author) |
|---|---|---|
| train_0001 | Head | Major improvement: substantial skull displacement corrected across all three planes |
| train_0002 | Head | Major improvement: rotational and positional mismatch markedly reduced |
| train_0003 | Head | Clear improvement: previously modest misalignment further reduced |
| train_0001 | Pelvis | Moderate improvement: better pelvic contour correspondence |
| train_0002 | Pelvis | Minor improvement: baseline already reasonably aligned |
| train_0003 | Pelvis | Moderate improvement, particularly around the pelvic contours |

All six comparisons showed visually improved alignment, with no obvious gross deterioration in the displayed views. The 34 mm head correction in train_0001 corresponds to a substantial initial skull displacement visible on the overlays. **No correction was automatically adopted.** Adoption remains governed by the locked selection rule, bootstrap threshold and safety monitor, followed by the visual-veto stage.

*Authorisation (author, 2 November 2026):* full 100-patient Amendment 3a run authorised with the existing parameters; Amendment 2 preserved as comparator and fallback; all adoption decisions provisional; selection results, bootstrap behaviour and safety monitor examined before the visual-veto review; no parameter change based on the results without a documented amendment.

*Bootstrap limitation (identified on code review, 2 November 2026):* in the current implementation, each bootstrap replicate's convergence status is **not recorded**. Non-converged replicates are included in the bootstrap SD as if they were successful fits. The bootstrap failure rate is therefore not available from this implementation.

**15.3.7 Full Amendment 3a numerical run (run 20261102_073040; 2 November 2026).** Script MD5 b21a405145526d56a59bb5950ef53acf (Amendment 3a parameters; bootstrap 200 full refits; seed 20261101). 100/100 patients processed. **Runtime 21,346 s (about 5.9 h)** on the author's PC, recorded for reproducibility.

| Region | Patients | Estimable | Converged | Within caps | Provisionally adopted | One-way flag | Safety denominator | Worse > 2 mm | Fraction | Trigger |
|---|---|---|---|---|---|---|---|---|---|---|
| Head | 100 | 100 | 100 | 92 | 91 | 5 | 91 | 0 | 0.000 | no |
| Pelvis | 100 | 100 | 100 | 98 | 89 | 0 | 89 | 11 | 0.124 | **yes** |

Decisions: head: 91 provisionally adopted (visual veto pending); 8 exceed caps (rejected, not clipped) → Amendment 2; 1 selection improvement ≤ 1.0 mm → Amendment 2. Pelvis: 89 met the selection rule; 9 selection improvement ≤ 1.0 mm → Amendment 2; 2 exceed caps → Amendment 2.

**Pelvic safety-monitor trigger:** 11/89 (12.4%) of patients meeting the pelvic selection rule were worse by > 2 mm on the reporting structures (femoral heads), exceeding the predefined 10% limit. Per R7, **further pelvic adoption is suspended pending review**; the pelvis stays on Amendment 2 for all patients until a documented review and, if needed, a dated amendment. No transform is replaced automatically. **The safety rule is not changed.**

**Author's interpretation and review plan (2 November 2026):**
- Head: promising (91/100 provisional, no safety violations); final anatomical verification required.
- Pelvis: femoral-head deterioration may partly reflect a change in femoral orientation between examinations (the femora articulate with the pelvis), so femoral-head disagreement is not necessarily evidence of pelvic-bone misregistration. This cannot be concluded from the summary numbers.
- Sequence: (1) document the results (this section); (2) investigate the 11 pelvic cases first: sacrum/S1 improvement vs femoral-head deterioration, vessel/muscle measures and CT overlays, with the key question whether sacrum/S1 improved substantially while the femoral heads deteriorated; then decide whether the femoral-head metric is suitable as a pelvic registration safety measure; (3) review the head cases, including all five one-way flags; (4) complete the planned visual segmentation QC of the original CTs before regional registration decisions are finalised.
- Runtime is a reproducibility detail, not a scientific concern.

**15.3.8 Defective femoral-head reporting measure identified (2 November 2026; post-results).** Descriptive review of `Amendment3_fit_20261102_073040.csv` (MD5 dc74686ed0b024d8524668f9b9e16c94):
- Femoral-head centroid distances were physically implausible in several patients (e.g. 870–906 mm, 496 mm, 102 mm, 44 mm); across 100 patients, 8–13 per side exceeded 100 mm. 8 of the 11 triggering patients had an implausible value on at least one side.
- Probable cause (to be confirmed visually): the head was defined as the superior-most 40 mm of **all** femur-labelled voxels, so a stray superior femur-labelled fragment would be selected instead of the true femoral head. This is an implementation defect in the reporting measure (Claude), not evidence of pelvic misregistration.
- Triggering 11 vs other 78 adopted: sacrum/S1 centroid improvement median 3.1 mm in both; sacrum/S1 surface distance improved in all 11; vessel and muscle surface distances improved in all 11 (median −1.1 and −0.6 mm). No patient had both femoral heads worse; sides moved in opposite directions in 9/11. The 3 triggering patients with plausible values (train_0007, 0048, 0067) showed one head better and the other worse by about 2.4 mm.
- These findings do not establish that the 89 pelvic corrections are acceptable.

**Author decision (2 November 2026):** pelvic adoption remains **suspended**. (1) Visually inspect the 11 triggering patients to confirm that stray superior femoral fragments explain the extreme values; (2) correct the femoral-head definition to the proximal femoral component anatomically associated with the acetabulum, with connected-component checks (not simply the largest component); (3) recompute for all 100 patients from the **saved transforms** (no ICP refitting); (4) reapply the **unchanged** safety threshold (2 mm, 10%); (5) document this as a dated post-results amendment, keeping the original trigger result and identifying the defective measure. Replacing the femoral heads with another construct (vessels/muscles) was rejected, because it would change the safety-monitor construct after its results were seen.

**Author refinements to the proposed fragment rule (2 November 2026):**
1. Proximity to the hip bone (10 mm) is a starting value only. Proximity alone does not establish that a component belongs to the femoral head, so stray fragments near the pelvis could be retained.
2. No 5 mL minimum per component. A genuine femoral-head fragment separated from the shaft may be smaller. Volume plausibility is applied to the reconstructed proximal structure, and individual fragments are checked anatomically.
3. The superior-most 40 mm is not necessarily the femoral head, even after stray fragments are removed (it may include the femoral neck). The head should be defined relative to the acetabulum.
- **Not yet approved:** the 10 mm and 5 mL thresholds and the unchanged 40 mm definition. They are decided only after the diagnostic images show whether the rules identify the intended anatomy reliably.
- **Approved:** a diagnostic script (connected components, anatomical proximity, overlays, recorded exclusions); testing all 11 triggering patients plus randomly selected controls; use of the saved transforms without repeating ICP.
- The original femoral-head metric is kept as a documented historical result. If its anatomical definition changes, the revised metric is identified explicitly here.

**15.3.9 Femoral-head diagnostic script (2 November 2026; diagnostic only).** Script `02_code\Amendment3_femhead_diagnostic.py`, MD5 c639d9151dfc7427c15a924b2daae2f5 (author confirmed with Get-FileHash); synthetic self-test 7/7 (stray superior fragment found and not anchored; the original top-40 mm measure falls entirely in the stray fragment; sphere fit recovers the head centre and radius within 1.5 mm; laterality correct). It imports the unchanged functions of `Amendment3_regional_rigid.py` (MD5 b21a405145526d56a59bb5950ef53acf).
- **Scope:** no decision; no recomputed femoral-head distance or safety monitor; no refitting. Each CT is examined in its own native space (no transform).
- **Patients:** the 11 triggering patients; other patients with any femoral-head value > 20 mm (Amendment 2 or 3a); 6 random controls from the remaining provisionally adopted patients (seed 20261102).
- **Per patient, CT and side:**
  - femur connected components (26-connectivity) with volume, centroid, craniocaudal extent, minimum distance to the ipsilateral and contralateral hip-bone surface, and laterality;
  - the share of the ORIGINAL top-40 mm voxels lying in non-anchored components (tests the stray-fragment explanation);
  - candidate selections for judgement only: anchored set (components within 10 mm of the ipsilateral hip bone; not approved); acetabulum-based head = sphere fitted to anchored femur surface points within 6 mm of the ipsilateral hip-bone surface (articular contact; not approved), with radius, RMS residual and the distances from the original and anchored top-40 mm centroids to the sphere centre (tests whether "top 40 mm" captures the head or the neck);
  - reproduction of the original Amendment 2 femoral-head value exactly as in the main script, compared with the fit table (agreement required to 1e-3 mm).
- **Outputs:** `06_output\Amendment3_femhead_diagnostic_<run>\` with one PNG per patient (rows PSMA/FDG; whole-FOV coronal femur-component view, and left- and right-hip coronal slices through the sphere centre), components CSV, per-side CSV, report and manifest.
- **Status:** awaiting the author's run and visual review. Thresholds and the revised head definition are to be decided after the review, as a dated post-results amendment.

**15.3.10 Femoral-head diagnostic: results and author's visual review (run 20261102_144950; 2 November 2026).** Self-test 7/7; 27 patients (11 trigger, 10 implausible-other, 6 controls); runtime 1,482 s. Reproduction of the original Amendment 2 femoral-head values: **108/108 CT-sides agree** with the fit table. Per-side table `Amendment3_femhead_sides_20261102_144950.csv`, MD5 8898ceb70a65075a31edf10a4d4bb578.

*Numerical findings (descriptive):*
- The original top-40 mm region lay in femur components not near the ipsilateral hip bone in 31 CT-sides (20 patients): 8 of the 11 triggering patients, all 10 implausible-other patients and 2 controls (train_0041, train_0083; the error cancelled in the Amendment 3a − Amendment 2 difference). These fragments were often hundreds of mm from the pelvis. The 3 triggering patients without stray fragments (train_0007, 0048, 0067) show the genuine one-side-better/one-side-worse pattern (≈ 2.4 mm).
- The hip-bone label itself contains distant fragments in several patients (e.g. train_0031, 0040, 0065, 0089; in 6 CT-sides the ipsilateral hip label extends > 300 mm above the femur). Femur fragments adjacent to these are wrongly "anchored", and the sphere fit fails (radius 230–1,078 mm in 10 CT-sides, 6 patients: train_0017, 0022, 0031, 0040, 0065, 0089). The hip-bone craniocaudal range is therefore not a reliable plausibility check.
- Where selection was anatomically valid (98 CT-sides), sphere radius median 25.2 mm (IQR 24.5–26.8; RMS ≈ 1.1 mm), consistent with an adult femoral head; sphere-based head volume median 74 mL. However, 20 CT-sides had radius 30–40 mm with poorer fit (RMS ≈ 4 mm), **18 of 20 left-sided** (systematic; cause not established).
- Even after stray removal, the top-40 mm centroid lay a median 11.7 mm (4.1–22.0) from the sphere-fitted head centre, i.e. it includes femoral neck and trochanter.

*Author's visual review (train_0007, 0017, 0031, 0041):*

| Patient | Visual finding (author) | Assessment (author) |
|---|---|---|
| train_0007 | Both femoral heads identifiable; sphere radii ≈ 24.5–25 mm, good fits | Good reference case; the original top-40 mm region includes non-head anatomy |
| train_0017 | PSMA right and FDG left sphere radii ≈ 480–501 mm; selected slices far above the pelvis | Clear anatomical selection failure; stray fragments corrupt the measurements |
| train_0031 | Both PSMA sphere fits grossly wrong (≈ 489–1,078 mm); FDG right plausible; FDG left larger, less convincing | Severe selection failure, particularly on PSMA |
| train_0041 | Left femoral head identifiable but sphere radius ≈ 38–39 mm, including neck/trochanter; right-sided fits reasonable | Systematic left-sided fitting problem visually supported |

*Author's conclusions:*
1. The original top-40 mm method is **abandoned as a femoral-head definition** (it includes neck and trochanter even in clean segmentations).
2. Sphere fitting is promising but not yet reliable: it follows the femoral head where the anatomical selection is correct, but fits inappropriate structures when selection fails.
3. The hip-bone reference itself must be investigated; a contaminated reference mask can mislead the whole process. The images establish that anatomical selection fails, but not yet that every failure originates in the hip-bone mask rather than the sphere-fitting procedure.

*Decision (author, 2 November 2026):* Option C: a targeted diagnostic of hip-bone connected components (which component contains the acetabulum; distant fragments; locations on native CT), prioritising train_0031 and train_0017, then train_0041. No change to the ICP transforms. Pelvic adoption remains suspended until the measurement is corrected and validated.

**15.3.11 Hip-bone component diagnostic script (2 November 2026; diagnostic only).** Script `02_code\Amendment3_hipbone_diagnostic.py`, MD5 3abd7ea0560a693c5a4e6fc03064f0db; synthetic self-test 6/6 (distant hip-bone fragment not linked; largest hip component touches sacrum and femur; femur fragment adjacent to the distant hip fragment not anchored; trimmed sphere recovers both head centres and radii; femoral-neck contamination biases the plain fit while the trimmed fit recovers the radius). Imports the unchanged `Amendment3_regional_rigid.py` (b21a4051…) and `Amendment3_femhead_diagnostic.py` (c639d915…).
- **Scope:** no decision; no recomputed femoral-head distance or safety monitor; no ICP refitting; each CT in its own native space.
- **Patients:** the same 27 as run 20261102_144950; priority train_0031, train_0017, train_0041 (files prefixed PRIORITY_).
- **Per CT and side:** hip-bone components (26-connectivity) with volume, centroid, craniocaudal extent, distance to sacrum (sacroiliac joint), contralateral hip bone (symphysis) and ipsilateral femur, and height above the sacrum; candidate main component = largest; candidate linked set = components within 10 mm of the largest or of the sacrum (not approved). Femur anchoring repeated against the linked hip set only (10 mm, not approved). Sphere fits to femur surface points within 3 mm and 6 mm of the linked hip surface, plain and trimmed (iterative refit on the 80% of points with the smallest residual; probes the left-sided 30–40 mm radii). Comparison with the previous radius (uncleaned hip label).
- **Known limitation of the trimmed fit:** it removes moderate neck contamination, but a strongly biased starting fit can converge to a wrong sphere; the images remain the reference.
- **Outputs:** `06_output\Amendment3_hipbone_diagnostic_<run>\` with one PNG per patient (coronal MIP: linked hip blue, unlinked hip fragments orange, sacrum grey, anchored femur green, other femur red; left and right hip coronal slices through the trimmed sphere centre with contact points), hip-components CSV, per-side CSV, report, manifest.

**15.3.12 Hip-bone diagnostic: results and author's visual review (run 20261102_192737; 2 November 2026).** Self-test 6/6; 27 patients; runtime 1,889 s. Per-side table `Amendment3_hipbone_sides_20261102_192737.csv`, MD5 279560d875c3e6422306249f143a74ce.
- The largest hip-bone component touched the sacrum and the ipsilateral femur in **108/108** CT-sides. 39/108 CT-sides had > 1 hip-bone component; 37 had unlinked fragments (volume ≤ 8.4 mL, median 0; some ≈ 700 mm above the sacrum).
- With the cleaned (linked) hip reference, femur anchoring retained exactly one component in 106/108 CT-sides and two in 2 (train_0024 FDG right, train_0087 PSMA right).
- Sphere radius in 18–30 mm in **108/108** CT-sides (previously 78/108 with the uncleaned hip label; 10 > 100 mm, 20 in 30–100 mm): 6 mm plain median 25.0 (22.7–29.6), RMS 1.1; 6 mm trimmed 25.1 (22.7–27.6), RMS 0.7; 3 mm trimmed 24.9 (22.3–27.5), RMS 0.5. Left/right medians 25.0/25.2 mm: the left-sided 30–40 mm radii were also caused by hip-label fragments.
- Agreement of the same femoral head between CTs (|PSMA − FDG| radius, 6 mm trimmed): median 0.19, P95 0.60, max 0.73 mm (plain fit max 2.76 mm, train_0087 right).
- Even with clean labels, the top-40 mm centroid lay 3.8–18.8 mm (median 10.8) from the trimmed sphere centre.

*Author's visual review:* priority cases resolved; two-piece femurs:

| Patient | Finding (author) | Assessment (author) |
|---|---|---|
| train_0024 | FDG right femur: two retained components; fitted head correctly located at the hip joint, with marked abnormal CT density and contour irregularity; radius ≈ 23 mm; 6 mm trimmed RMS ≈ 1.8 mm (PSMA), 1.4 mm (FDG) | Localisation acceptable; right side flagged for additional review, not accepted automatically |
| train_0087 | PSMA right femur: two retained components; sphere follows the femoral head; radius ≈ 26.5 mm; RMS ≈ 0.7 (PSMA), 0.6 mm (FDG) | Localisation acceptable |

*Author's conclusions:* anatomically connected proximal components should be retained rather than requiring one continuous femur label. The revised procedure (cleaned hip reference + trimmed 6 mm sphere) is a **defensible candidate measurement method ready for full-cohort testing, not yet a fully validated one.** Abnormal fits are flagged automatically for visual review; a plausible radius alone does not establish acceptance.

**15.3.13 Amendment 3b — corrected femoral-head reporting measure (DRAFT, dated 2 November 2026; POST-RESULTS).**
- **Disclosure:** drafted after the Amendment 3a full-run results, the femoral-head diagnostic (27 patients) and the hip-bone diagnostic (27 patients) had been examined. The original femoral-head measure and its safety-monitor result (11/89, 12.4%, triggered) are retained and reported permanently as historical results, identified as defective.
- **Corrected definition (proposed):** per CT (native space) and side:
  1. hip reference = largest hip-bone component plus components within 10 mm of it or of the sacrum;
  2. proximal femur = femur components within 10 mm of the hip reference (several components allowed);
  3. contact points = proximal-femur surface voxels within 6 mm of the hip-reference surface;
  4. femoral-head centre = centre of a sphere fitted by iterative trimmed least squares (best 80% of points by residual, up to 10 rounds).
- **Measure:** distance between corresponding fitted femoral-head centres: under Amendment 2, ‖c_PSMA − T_A(c_FDG)‖; under Amendment 3a, ‖c_PSMA − ΔT(T_A(c_FDG))‖ (exact point mapping with the saved transforms; no resampling, no ICP refitting). Change = Amendment 3a − Amendment 2.
- **Safety monitor:** unchanged (worst evaluable side worse by > 2 mm; denominator = provisionally adopted pelvic patients with ≥ 1 evaluable side; trigger > 10%).
- **Evaluability and automatic flags (values to be approved):** not evaluable if no proximal femur, < 50 contact points, or either centre outside the common field of view; flagged for visual review if radius outside 20–30 mm, trimmed RMS > 1.2 mm, |radius PSMA − FDG| > 1.5 mm, or > 1 retained femur component. Flagged sides are not accepted automatically.
- **Status (2 November 2026):** APPROVED by the author for full-cohort evaluation (see 15.3.14).

**15.3.14 Amendment 3b — author approval and implementation (2 November 2026).**
- **Approved values:** hip-bone cleaning 10 mm; femur anchoring 10 mm; contact band 6 mm; trimmed sphere fit; minimum 50 points (provisional), applied to the points **retained in the final trimmed fit**, not to the initial candidate points; flags for review: radius outside 20–30 mm, RMS > 1.2 mm, |radius PSMA − FDG| > 1.5 mm, > 1 retained femur component. These are **study-specific provisional QC criteria**, identified as such.
- **Modification (author):** a femoral-head centre outside the common field of view is **not evaluable for the paired comparison**; it is not evidence of incorrect segmentation.
- **Same definition under both transforms:** per side, Δd = d(Amendment 3a) − d(Amendment 2); positive = deterioration, negative = improvement. The 2 mm deterioration and 10% patient-level thresholds are unchanged.
- **Aggregation preserved exactly:** patient worst change = maximum over evaluable sides (unilateral deterioration counts); patients without an evaluable side are excluded from the denominator; denominator = provisionally adopted pelvic patients with a worst change; trigger if > 10%.
- **Option A (author):** PRIMARY = all technically evaluable sides, flagged sides included; SENSITIVITY = flagged sides excluded. The safety-monitor decision remains **provisional** until flagged sides are visually adjudicated (train_0024 right is not accepted automatically).
- The original defective result (11/89, 12.4%) is preserved and reported alongside.
- **Implementation:** `02_code\Amendment3b_femhead_evaluation.py`, MD5 8bbd017bac4a6306fcfd31adfbeb0ad5; self-test 9/9 (known 6 mm offset gives Δd = −6 mm; wrong-direction correction gives positive Δd; centre outside FOV not evaluable; clean heads not flagged; RMS and multi-component flags; no fit not evaluable; max-over-sides aggregation; NaN excluded from denominator; exactly 10% does not trigger). Uses the saved Amendment 3a transforms (`Amendment3_transforms_20261102_073040.csv`) with exact point mapping; **no ICP refitting**. The femoral-head method is `analyse_ct()` of the hip-bone diagnostic, imported unchanged. Per patient it checks transform composition and that the common field of view equals the Amendment 3a fit table. Outputs: per-side, per-patient, summary (original / primary / sensitivity) and flagged-side tables; one review PNG per patient with a flagged side; report; manifest.
- **First run aborted (2 November 2026, about 4,267 s):** all 100 patients were processed and the CSV tables written, but the report step failed (`KeyError: 0`). The code read `ev.flags`, which in pandas is the DataFrame's built-in `.flags` attribute, not the "flags" column (the same attribute-shadowing error as earlier in h1/h1b/h2z). No report or manifest was written, so that output folder is **incomplete and is not used**. A first fix (MD5 6b9f71788c02e2c0b2bec454baaab84f) was **withdrawn before use**: it put `ev["flags"]` inside a double-quoted f-string, which is a syntax error on the author's Python 3.11 (verified with Python 3.11 here).
- **Author decision (2 November 2026):** do not repeat 70 minutes of image processing only to generate a report. Recover from the saved tables if their integrity can be verified; otherwise rerun the full pipeline.
- **Corrected script** MD5 f6df42596ef28af909d90f00d45fa48c; self-test 9/9; compiles under Python 3.11. The report code is in two functions, `flag_counts()` and `report_lines()`, with column access by key. New `--report-only RUN_FOLDER --processing-md5 MD5` mode: reads the saved tables, reads no image and recomputes nothing from images. It runs 12 integrity checks:
  1. required columns;
  2. 200 rows = 100 patients × 2 sides, no duplicates;
  3. patient IDs identical to the Amendment 3a fit table;
  4. Δd = d_A3a − d_A2 for every evaluable side;
  5. ≥ 50 retained points;
  6. flags re-derived from the stored radii, RMS and component counts;
  7. adoption status and original values agree with the fit table;
  8. patient worst changes re-aggregated from the sides;
  9. safety summary recomputed;
  10. ORIGINAL row reproduces 11/89;
  11. flagged-side list;
  12. one review PNG per flagged patient.
  
  Any failure means no report and a full rerun. The mode refuses to overwrite an existing report. Report and manifest record the processing-script MD5 (8bbd017bac4a6306fcfd31adfbeb0ad5, which produced the tables) and the reporting-script MD5 separately. Tested on synthetic tables: 12/12 pass, and a 0.5 mm alteration of one Δd is detected (integrity fail, no report).

**15.3.15 Amendment 3b results (run 20261102_211242; report-only recovery, 2 November 2026; PROVISIONAL).** Tables computed by script MD5 8bbd017bac4a6306fcfd31adfbeb0ad5 (processing about 4,267 s); report and manifest generated by MD5 f6df42596ef28af909d90f00d45fa48c after **12/12 integrity checks** (the ORIGINAL row reproduces 11/89). Per-side table MD5 8aeff038bda56d60b4f396ca82f40374; per-patient table MD5 a3b376472073072dd8c9fe7720d37d88.

| Analysis | Denominator | Worse > 2 mm | Fraction | Trigger |
|---|---|---|---|---|
| ORIGINAL (defective top-40 mm; historical) | 89 | 11 | 0.124 | yes |
| **Amendment 3b PRIMARY (all evaluable sides)** | 89 | **0** | 0.000 | no |
| Amendment 3b SENSITIVITY (flagged sides excluded) | 88 | 0 | 0.000 | no |

- Evaluable sides 200/200 (no fit failure, no centre outside the common FOV). Radius median 24.8 mm.
- Provisionally adopted pelvic patients (178 sides): head-centre distance d_A2 median 5.71 mm (P90 12.57, max 23.82); d_A3a median 0.49 mm (P90 0.82, max 2.67); Δd median −5.20 mm (P90 −1.87, max +1.58). One adopted side was worse: **train_0024 right, +1.58 mm**, the side already flagged by the author (RMS, two components). Patient worst change median −3.90 mm (max +1.58).
- Sensitivity denominator 88: train_0085 is excluded because both its sides are flagged.
- The 11 originally triggering patients: worst change now negative in 10, and +1.58 mm in train_0024.
- **Flagged sides for visual adjudication (7 sides, 6 patients, all provisionally adopted):** train_0024 R (rms; multi_component; Δd +1.58); train_0056 R (rms; multi_component; −1.12); train_0073 R (multi_component; −5.26); train_0075 L (rms; radius_diff, radius 24.15 vs 20.34 mm; −3.29); train_0085 L (radius 29.94/30.09 mm; −4.26); train_0085 R (multi_component; −2.96); train_0087 R (multi_component; −14.39; already judged acceptable by the author).
- Patients not adopted (descriptive only): Δd also negative in all 11, worst change −1.15 to −31.73 mm, including the two rejected for exceeding the pelvic caps (train_0069 −26.12 mm, train_0089 −31.73 mm). **Recorded as an observation; no change to caps or selection rules.**
- **Limitation (coupling):** the femoral head sits in the acetabulum, part of the hip bones that are the pelvic fitting surface. A small d_A3a therefore largely confirms hip-bone alignment and is **internal consistency evidence, not independent accuracy evidence**. The improvement in non-adopted patients is expected by construction for the same reason. Advantage: in a ball-and-socket joint the head centre is insensitive to leg posture, which affected the original measure. Independent evidence remains the vessel/muscle measures and the visual-veto stage.
- **Status:** the safety-monitor result is provisional until the 7 flagged sides are visually adjudicated. Lifting the pelvic suspension is the author's decision, to be made after that adjudication.

**15.3.16 Pelvic safety-monitor investigation CLOSED (author decision, dated 2 November 2026).**

*Visual adjudication of the 7 flagged femoral-head measurements (author):*

| Patient | Side | Decision |
|---|---|---|
| train_0024 | Right | Accept with qualification; **suspected hip prosthesis (unconfirmed)** |
| train_0056 | Right | Accept |
| train_0073 | Right | Accept |
| train_0075 | Left | Accept with qualification; **suspected hip prosthesis (unconfirmed)** |
| train_0085 | Left | Accept; radius flag retained |
| train_0085 | Right | Accept |
| train_0087 | Right | Accept |

*Decision:*
- Original defective result 11/89 (12.4%) is retained as a historical record; corrected result (Amendment 3b) 0/89 worse by > 2 mm; the original 10% threshold is passed; the seven flagged measurements have been visually reviewed.
- **Pelvic suspension LIFTED.** The 89 pelvic corrections remain **provisionally eligible**.
- **This closes the safety-monitor investigation only.** It does not declare the 89 pelvic registrations anatomically validated. Validation continues with the planned visual segmentation assessment of the 200 CTs (Step 1) and the visual CT-to-CT registration assessment of the final regional framework (Step 2: head and pelvis Amendment 3a where provisionally adopted; neck, thorax and abdomen Amendment 2), both with predefined scoring checklists. Within-examination PET–CT alignment is to be verified before the final transforms are applied to PET.
- The suspected prostheses in train_0024 (right) and train_0075 (left) are unconfirmed and documented as such. They are carried forward as items for Step 1 (segmentation near metal) and for the PET stages (possible metal artefact).

## 16. Provisional methodological framework (v2.7, 30–31 October 2026) — NOT LOCKED
Agreed in discussion after Audit 1D, using only within-tracer technical data (Audits 1A–1D). No correspondence or discordance result has been examined. To be formalised in a dated amendment and locked at Step 2; until then every item below may still change.

### 16.1 Agreed provisional decisions
| Item | Provisional decision |
|---|---|
| Primary connectivity | 26 (as originally favoured; not switched after the audit) |
| Connectivity sensitivity | 6 and 18 |
| FDG eligibility | Minimum physical volume plus minimum voxel count; determined in native FDG space; repeated independently for each connectivity (no transfer of eligibility between partitions). Approved candidates: primary ≥ 0.1 mL AND ≥ 3 voxels; sensitivity ≥ 0.25 mL AND ≥ 3 voxels (Section 16.3) |
| PSMA match targets | All supplied PSMA TTB-positive voxels; no minimum-volume exclusion from correspondence |
| Native masks and volumes | Unchanged; native FDG volume is the sole source of burden |
| Matching rule | Assessable FDG component matched if its minimum physical boundary distance to any PSMA-labelled region is ≤ the reference tolerance; one-to-many permitted |
| Registration | Supplied rigid transformation (Section 15 outcome) |
| Registration QC vs tolerance | Separate parameters; a region failing QC → correspondence-indeterminate; the tolerance is never widened to compensate for misregistration |
| Voxel-level sensitivity | "Spatially unmatched FDG tumour-volume fraction" (FDG voxels farther than the tolerance from any PSMA voxel); secondary methodological sensitivity only; not a measure of intralesional heterogeneity (Section 12, intralesional mismatch not an objective); indeterminate voxels never counted as unmatched |
| Perturbation analysis | Lesion-level DMTVF primary; voxel-level fraction secondary; identical predefined perturbations |
| Largest-component descriptor | f_largest = largest FDG component volume / total FDG TTB (26-connectivity), reported per patient |
| Watershed / anatomical splitting | Not introduced |
| Manuscript | Main text: primary endpoint, cohort burden, anatomical distribution, registration-related uncertainty, one robustness summary. Supplement: connectivity, minimum-volume, small-component burden, voxel-level, perturbations, detailed QC |

### 16.2 FDG volume accounting (four mutually exclusive categories)
V_FDG,total = V_S + V_I + V_M + V_U, assigned in fixed order:
1. V_S — below-threshold (native-space eligibility; unassessed, not presumed matched)
2. V_I — correspondence-indeterminate (eligible but correspondence not reliably assessable)
3. V_M — matched
4. V_U — unmatched

- Primary endpoint unchanged: DMTVF_ref = V_U / V_FDG,total.
- Primary uncertainty interval unchanged (Section 5): DMTVF_ref to (V_U + V_I) / V_FDG,total.
- f_small = V_S / V_FDG,total reported separately; not included in the uncertainty interval (eligibility choice ≠ anatomical uncertainty).
- Components shrinking or disappearing during resampling are handled by Section 14.3 (correspondence-indeterminate), never as below-threshold.
- Monotonicity note: without an FDG minimum volume, broader connectivity can only keep or reduce V_U under the any-contact rule; with a minimum volume, merged sub-threshold fragments can become eligible, so V_U can increase or decrease. Connectivity sensitivity results are interpreted accordingly.

### 16.3 Open decisions (to resolve before the Step 2 lock)
1. **FDG minimum volume — APPROVED CANDIDATE RULES (author, 31 Oct 2026; NOT LOCKED until the dated amendment and the Audit 1D report).** Technical justification by voxel sampling, not by expected discordance; these are minimum sampling requirements, not the minimum size of a biologically valid metastasis; the combined rule does not standardise effective spatial resolution across scanners (limitation). P90 excluded burden used as a feasibility diagnostic, not as the selection criterion. Each rule applied in native FDG space, independently for each connectivity.
   | Rule | Definition (component eligible if) | 26-conn: components excluded | Pooled FDG volume excluded | Per-patient median / P90 / max | Patients > 5 / 10 / 20% excluded | Patients with no assessable component |
   |---|---|---|---|---|---|---|
   | Primary candidate (A) | ≥ 0.1 mL AND ≥ 3 native voxels | 871 | 0.24% | 0.27% / 1.7% / 27.8% | 3 / 1 / 1 | 0 |
   | Sensitivity candidate (B) | ≥ 0.25 mL AND ≥ 3 native voxels | 1,389 | 0.61% | 0.89% / 5.9% / 50.0% | 16 / 2 / 2 | 0 |
   - Same rules under 6-connectivity: A 0.36% pooled (P90 2.9%, max 30.0%), B 0.78% (P90 6.7%, max 51.1%); under 18-connectivity: A 0.25% (P90 1.7%, max 27.8%), B 0.63% (P90 6.1%, max 50.0%). No patient loses all assessable FDG components under either rule or connectivity.
   - At 0.25 mL the 3-voxel criterion excludes no additional component; it is kept for consistency between the rules.
   - Patient-level excluded fractions are methodological descriptors only, never grounds for excluding patients; f_small reported separately from V_I. Maxima arise in patients with very small total FDG TTB (e.g., train_0006 4.5 mL).
   - 0.5 and 1 mL judged not defensible as primary (per-patient P90 12–23%).
   - Source: `02_code\Audit1D_optionD_small_component_summary.py` (low-burden descriptors section; outputs Audit1D_optionD_low_burden_*).
2. **Reference matching tolerance and regional registration-QC limits** — designed together from Audit 1C; the 7.3 mm median bony centroid distance is a registration metric, not a lesion-localisation error, and is not adopted directly as the tolerance. Indicative values discussed (trunk QC ≈ 10 mm; tolerance ≈ 8–10 mm) are NOT candidates for locking: they rest on in-sample residuals and must be justified independently after the registration method is settled (Amendment 2 outcome).
3. **QC inheritance** for lesions in unlabelled tissue: local landmarks → predefined anatomical region (compartment-aware; distance to the nearest landmark alone is insufficient) → correspondence-indeterminate.
4. **FDG components spanning several QC regions** — candidates: (a) strict: indeterminate if any part lies in a failing region; (b) volume-weighted: indeterminate if more than a fixed share lies in failing regions; (c) split: passing part assessed, failing part counted as V_I.
5. Distance definition (boundary-voxel centres vs voxel boundaries) and handling of partial PSMA field-of-view coverage.

### 16.4 Next steps
1. Audit 1C-f — COMPLETED (Section 8.4.11). Next: pelvic overlays of train_0069 and train_0037; combine regional QC evidence; [superseded sequence:] author approval of the specification → script with built-in synthetic-geometry test → MD5 lock → single run (supplied and adopted transformations; descriptive only).
2. Literature verification for the reference matching tolerance τ (boundary displacement of threshold-defined lesions, EARL performance, TG-132/MIRSIG applicability); τ fixed and documented before any acceptance ratio is computed.
3. Head/neck: separate regional QC; decide between correspondence-indeterminate status and a separately specified head/neck correction.
4. Regional acceptance rules (exploratory assessment, then lock).
5. Audit 1D written report (observations only) — closes Audit 1D.
6. Dated methodological amendment (decisions and rationale; separate from the audit report).
7. Step 2: lock parameters and statistical analysis plan.
8. Only then: PSMA–FDG correspondence processing.

### 16.5 Registration acceptance framework — discussion record (31 October 2026; NOT LOCKED)
- Three separate concepts: registration improvement; registration acceptance (reliability for correspondence assessment); lesion-matching tolerance. Improvement does not imply acceptance; the tolerance never compensates for misregistration.
- Proposed structure: R = e_r / τ (regional residual alignment error / reference matching tolerance). Candidate categories R ≤ 0.5 favourable, 0.5 < R ≤ 1 uncertain, R > 1 unreliable — engineering rules, exploratory only, not validated standards; not yet used to define the primary endpoint.
- τ: to be derived from the matching metric (boundary-voxel-centre distance; touching masks ≈ one voxel spacing of the coarser grid) and a published boundary-displacement term for threshold-defined lesions; registration error excluded from τ. The earlier proposal "FWHM + half voxel diagonal" was withdrawn. PET resolution ≠ voxel size; EARL harmonises recovery coefficients, not one effective FWHM.
- Acquisition metadata: scanner models across the cohort reported by the data owners [VERIFY exact source]; per-scan scanner, reconstruction and FWHM unavailable; NIfTI header description/auxiliary fields empty in the metadata audit (100 cases; tracer coverage to be confirmed). PET grid clusters (Audit 1D) usable only as a stratification proxy; in 21 patients PSMA and FDG grids differ (coarser grid relevant).
- e_r: from genuinely held-out structures; leave-one-out fitting-structure residuals withdrawn as substitutes (median of ~20 structures → nearly in-sample). Lesion-type-specific: bone lesions ← held-out bone; nodal/soft-tissue lesions ← vascular/muscular structures; solid organs descriptive only. ≥ 5 reliable structures per region; robust upper statistic (e.g., second-largest residual) rather than the maximum.
- Candidate independent landmarks (verified present in the 117-class dictionary and Audit 1B paired completeness): iliac veins (100/100), iliac arteries (100/100; fragmented in 34–38 patients), aorta, IVC, iliopsoas (99–100), gluteus medius/minimus (100), gluteus maximus (84–86), autochthon (100), prostate (98; exploratory). Pelvic vessels + muscles: median 12 paired-complete structures per patient (min 9). Residuals for these structures have not been computed.
- Audit 1C-f direction (author, 31 Oct 2026): core set iliac veins, iliopsoas, gluteus medius/minimus, aorta, IVC; iliac arteries only if fragmentation ≤ 5% on both CTs; gluteus maximus only if paired-complete; prostate exploratory. Vessels: symmetric MSD over the common anatomical extent (primary; centroid unsuitable for long vessels), with crop-plane surface points excluded; vessels are insensitive to displacement along their own axis (reported). Muscles: MSD with centroid as supporting measure; posture deformation acknowledged. Regions (pelvis, retroperitoneum, …) evaluated separately with predefined vertebral-level z-boundaries; long structures cropped per region. Both transformations (supplied; adopted with fallbacks); no new optimisation. Reporting: per-region median, P90, maximum, number of evaluable structures, patient-level differences. Metric and cropping verified on synthetic shapes with known offsets before any patient data are processed.
- Multi-region FDG components: primary endpoint remains lesion-level (four categories, one status per component); regional splitting only as a predefined sensitivity analysis. Primary multi-region rule pending (strict favoured in discussion; majority and split as sensitivity).

## 17. Anatomical validation plan (APPROVED by the author, 2 November 2026)

**Principle (author):** the aim is the most defensible anatomical correspondence, not the most impressive registration numbers. The next objective is to establish whether the anatomical data supporting the registration are trustworthy, not to improve an algorithm.

**Project tracker (2 November 2026):**

| Component | Status |
|---|---|
| Audit 1A dataset integrity; Audit 1B automated segmentation assessment | completed |
| Audit 1C registration assessment | numerical stages completed; **classification-rule validity amendment OPEN**; the 102 Audit 1C-h visual pages are folded into Phase 2 |
| Amendment 2 (translation) | completed |
| Amendment 3a (regional ICP) | computation completed; adoption provisional (visual veto in Phase 2) |
| Amendment 3b (femoral-head measure) and pelvic safety-monitor investigation | completed / closed (15.3.16) |
| Audit 1D written report | **OPEN** |
| Step 2 protocol lock (τ definition/combination, SAP) | **OPEN** |
| Phase 1 visual segmentation validation | pending (specification drafted) |
| Phase 2 visual CT-to-CT registration validation | pending |
| Phase 3 within-examination PET–CT alignment | pending |

**Phase 1: visual validation of TotalSegmentator masks (200 CTs).** Structures relevant to registration and later anatomical classification only (not all 117 labels); each CT against its own CT, plus PSMA–FDG consistency; five-point scale with error types recorded separately; predefined score-to-action rules; standardised automatically generated review sheets; machine-readable scores. Specification: `Phase1_segmentation_QC_spec_DRAFT_v2_20261102.md` (supersedes the v1 draft). Author decisions (2 Nov 2026): O1 organ list deferred to the SAP lock; O2 selection structures: 4–5 acceptable, 3 targeted assessment, 1–2 unreliable unless a documented check shows the error did not affect the selection measurement, and an unreliable selection measurement rejects the Amendment 3a correction; O3 ≥ 14 days; O4 quadratic-weighted kappa; O5 each vertebra scored, with individual-landmark (A) and regional (B) reliability distinguished, so that one unreliable vertebra does not invalidate a T1–L5 registration; O6 sessions ≤ 2 h. O7 approved as a screening flag, not an automatic rejection: > 3 unreliable fitting vertebrae or ≥ 2 adjacent unreliable vertebrae → flag for Phase 2; otherwise document; Phase 2 decides whether the remaining anatomy supports reliable registration. The decisions are approved; **Author conditional approval (3 November 2026), with three corrections:**
1. Amendment 3b cleaning applies only to the femoral-head measurement. **Verified against `Amendment3_regional_rigid.py` (b21a4051…): the Amendment 3a pelvic ICP fitted the uncleaned hip labels** (likewise skull fitting and brain, sacrum and S1 selection), so cleaning cannot justify fragments in fitting or selection surfaces.
2. G9 is judged on the reliability of the fitted femoral-head centre, not on the whole femur label.
3. Reproducibility: agreement and score distributions are primary; kappa only where estimable; pooled kappa descriptive.

**O8 decided (author, 3 November 2026): option (c), staged.**
1. Phase 1 flags fragmented fitting and selection structures (location, size, apparent effect on the surface used).
2. Quantify whether fragments actually contributed points to the ICP fit (sampled, and retained by trimming).
3. Refit with cleaned masks only for the affected cases, all other settings unchanged, and compare with the saved transforms (rotation, translation, alignment on the selection and supporting structures).
4. No automatic replacement; any change needs a dated amendment.

Amendment 2 centroid sensitivity is assessed separately. A separate specification will follow and must be completed before final registration acceptance.

**Phase 1 specification LOCKED (3 November 2026):** `Phase1_segmentation_QC_spec_LOCKED_v2_1_20261103.md`, MD5 528a2436fb87250813c6d7d2abd95903. (Superseded earlier status: the decisions had been approved before the full text was reviewed; the text has now been reviewed and is locked.)

**3 November 2026 (author approval of both deliverables):**
- **Phase 1 review-sheet script:** `02_code\Phase1_segmentation_QC_sheets.py`, MD5 4d4a04ca5f32cab53c964c4c488d74cd, written against the locked v2.1 specification with no methodological changes. Self-test 10/10; compiles under Python 3.11. Implementation notes:
  - spine views are curved slab MIPs (± 12 mm around the vertebral-centroid line at each axial level), so that scoliosis and lordosis do not take vertebrae out of plane;
  - the femoral-head sphere uses the unchanged Amendment 3b `analyse_ct()`;
  - label IDs are checked against TotalSegmentator's own class map;
  - the code key and the hidden pre-screen are written outside the sheet folder;
  - spine vertebrae are pre-filled 5 with a "level reviewed" confirmation per spine group.
  Status: awaiting author MD5 check, self-test and preflight.
- **O8 sensitivity specification:** `O8_fragment_sensitivity_spec_DRAFT_20261103.md`, DRAFT for review; no script. Stage A estimates the fragment contribution for all patients; Stage B refits selected patients with matched seeds, uncleaned vs cleaned; Stage C covers Amendment 2 centroids. It states explicitly that a negligible estimate cannot prove absence of effect, because the original subsamples were not saved. The six-hour exact replay is not planned (author). Phase 1 visual QC proceeds in parallel.
- **Phase 1 preflight (3 November 2026):** MD5 confirmed; self-test 10/10; preflight on train_0024 and train_0031 (TEST codes) rendered correctly (author).
- **Amendment v2.2 to the Phase 1 specification (dated 3 November 2026, author-approved): two-stage review.** Specification `Phase1_segmentation_QC_spec_v2_2_20261103.md` (first version MD5 6a931b305c4c859aec6d8d7c3677a6eb; corrected version below); §§1–11 unchanged as the Stage 2 criteria.
  - Stage 1 screens all 200 CTs: head, spine and pelvis scores 1–5 (4–5 pass = "no significant abnormality detected at screening", not individual validation), plus three essential items (fragment in a fitting/selection structure; metal/implant; femoral-head centre reliability, left and right separately).
  - Stage 2 applies the locked structure-level criteria to regions scored ≤ 3 or flagged by the hidden automated pre-screen. **Automated flags trigger targeted review, never automatic rejection.**
  - Both screening and detailed results are preserved. The repeat reading uses the screening scores.
  - Script `Phase1_segmentation_QC_sheets.py` updated: review images unchanged; new screening templates and README; hidden region flags. MD5 3c720cbeced4c9cb229018abf4c7d80f, self-test 12/12, compiles under Python 3.11 (supersedes 4d4a04ca…).
  - Open P1: automated flag thresholds (provisional).
- **Second preflight (3 November 2026, train_0024 and train_0031):** author checked the screening CSV (4 × 13), the consistency CSV (6 × 4), the region scores, the fragment and implant fields, separate left/right femoral-head fields, the Stage 2 referral rules and the absence of duplicates. Approved with two adjustments:
  1. the fragment item includes the T1–L5 vertebrae (Amendment 2 fitting);
  2. a doubtful or unreliable femoral-head centre sends the pelvis to Stage 2 regardless of the pelvis score.
  
  README and specification updated; review images unchanged (no regeneration needed). Script MD5 a9c72452fd08105734e5a1432b3f072c (self-test 12/12); specification v2.2 MD5 d9de62ef66970171442467894bf5588b. The full 100-patient run is authorised.

**Phase 2: visual validation of CT-to-CT registration.** Final framework: head and pelvis Amendment 3a where provisionally adopted, neck/thorax/abdomen Amendment 2. Head and pelvis: Amendment 2 vs Amendment 3a overlays in randomised left/right order. Assess actual CT anatomical correspondence, not only mask overlap, in axial, coronal and sagittal views. The 102 Audit 1C-h pages are kept as historical QC material, and their outstanding review requirement is incorporated so that previously flagged regions are not overlooked. Five-point registration scale; quantitative measures reported alongside, not overriding. **The locked visual-veto rule is preserved:** a visual finding can reject an Amendment 3a correction but never force one. Rejecting Amendment 3a does not establish that Amendment 2 is satisfactory; the retained transform must pass the applicable registration QC. Blinding is partial (the author has seen preflight and some pelvic overlays); this is acknowledged. Deliverable: a patient-by-region decision table.

**Phase 3: within-examination PET–CT alignment (200 PET/CT).** Reference regions: brain (where covered), kidneys, urinary bladder, liver, diaphragm. Movement- and respiration-related mismatch is scored and documented, with a predefined consequence rule (to be specified). Reportable diaphragmatic mismatch must be defined, because some mismatch is near-universal. Head motion between CT and PET acquisition is considered, given the head-specific correction.

**Refinements adopted (2 November 2026):**
1. Status table corrected (above).
2. Operational score definitions, each linked to an action.
3. Exact structures specified (T1–L5 = Amendment 2 fitting set; femoral head = Amendment 3b sphere-fitted head and retained femur components). The organ list for region assignment is fixed before review or deferred until after the SAP lock, never chosen during review.
4. Intra-reader reproducibility: 10% blinded repeat reading, selected with a fixed random seed, after a predefined interval, analysed with weighted Cohen's kappa (or another ordinal agreement statistic). A second independent reader would strengthen the validation but is not essential.
5. Visual-veto rule preserved.
6. Consequences of PET–CT misalignment predefined.
7. Standardised review sheets and machine-readable outputs.

**Author corrections:**
- **Registration uncertainty ≠ V_I.** Regional registration uncertainty is recorded first, as a property of the registration. Lesion-level correspondence status (including V_I) is assigned later by predefined lesion-level rules; for example, a large lesion with substantial overlap may still have reliable correspondence despite moderate registration uncertainty.
- **No unsupported millimetre or voxel tolerances in visual score definitions.** Visual scores use observable, operational criteria: visible boundary disagreement, missing anatomy, disconnected fragments, wrong labels, and whether an error affects a fitting surface. Any quantitative tolerances are computed separately.

**Sequence:** each phase is documented before the next begins. Phase 1 → confirm mask reliability; Phase 2 → finalise regional transforms; Phase 3 → establish whether the transforms can be transferred to PET; then finalise the framework for lesion correspondence.
