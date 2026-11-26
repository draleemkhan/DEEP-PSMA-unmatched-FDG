# Spatially unmatched FDG tumour burden in mCRPC (DEEP-PSMA) — analysis code

**Status:** archival reproduction package with automated verification; not a general-purpose PSMA–FDG discordance pipeline.

Code accompanying the manuscript:

> Khan MA, Ali MB, Fatima M. *Quantifying spatially unmatched FDG tumour burden in metastatic castration-resistant prostate cancer: dependence on image registration, matching criteria and analysis granularity.* [Journal, year, DOI to be added]

The study quantifies FDG tumour volume spatially unmatched to PSMA-labelled disease in the 100 patients of the public DEEP-PSMA training dataset. It examines how this volume depends on registration strategy, matching definition, analysis granularity and a geometric stress test. The measure is a pipeline-dependent measure of spatial correspondence, not a direct measure of biological FDG–PSMA discordance.

## What the repository contains

The repository holds:
- the exact archived analysis scripts used to generate the reported results;
- the locked protocol, statistical analysis plan and dated implementation specifications;
- the registration-QC records;
- a runner that reruns the 66 archived analysis runs and checks their outputs against reference checksums.

**No DEEP-PSMA image data or regenerated numerical analysis outputs are included.** The repository includes four categories of non-regenerable records required for verification or reproduction:
1. the registration-QC records (human judgements and QC records);
2. the ST2 completion-note addendum, whose checksum a script verifies (it quotes some results);
3. the original run manifests (run parameters, checksums and run-level summary values);
4. the reference checksums.

| Path | Contents |
|---|---|
| `02_code/` | 47 Python and MATLAB scripts, archived exactly as run |
| `00_protocol/` | Final protocol (v2.10), locked statistical analysis plan, and the locked specifications and dated clarifications whose checksums the scripts verify |
| `06_output/` | Registration-QC records (see below), the locked QC analysis plan and the ST2 completion-note addendum. The scripts expect these files here, in the same layout as the original project |
| `reference/reference_checksums_v1_0.csv` | MD5 checksums of the 24,977 original output files of the 66 reproduced runs |
| `reference/manifests/` | The 266 original run manifests, in their original folder layout |
| `tools/reproduce_runs_v1_0.csv` | The 66 runs, in original order, with the exact arguments used |
| `tools/build_reference_v1_0.py` | Script that built the reference checksum list from the original outputs |
| `reproduce_runner.py`, `reproduce_config.json` | Reproduction runner and its configuration (Python and MATLAB executables) |
| `reproduction_test_20261126/` | Reports and logs of the reproduction test of 26 November 2026 (see below) |

### Registration-QC records

The visual registration QC was carried out blind and has been unblinded and frozen. The records are provided as non-regenerable inputs:
- the readers' completed checklists (PDF);
- the grades transcribed from them: the first-reader table, whose transcription was checked against the PDF before unblinding, and the text extracted from the second-reader and trunk checklists;
- the Phase 1 register (PDF and parsed table);
- the blinding keys, the random draw records, the control pool and the overlap survey;
- the final registration-QC decision table (`06_output/A3_visualQC_unblinded_20261108/A3_final_acceptance_per_region.csv`).

The keys link random QC codes to public DEEP-PSMA case identifiers (`train_0001`–`train_0100`). The reader notes describe anatomical findings visible on the public images only. The records contain no personal information. Trunk (neck) QC decisions are also encoded in `SAP_R1_region_rules.py`.

## Data

The DEEP-PSMA dataset is not redistributed. Download it from Zenodo:

> Jackson P, Hofman M, Buteau JP, McIntosh L, Sun Y. *A repository of annotated PSMA and FDG PET/CT images for algorithm development in staging of mCRPC for treament with 177Lu-PSMA therapy* [data set]. Version v1. Zenodo; 2025. https://doi.org/10.5281/zenodo.15281784 (licence CC BY-NC 4.0)

MD5 checksums of the five archives used in this study (version v1):

| File | MD5 |
|---|---|
| 0001-0020.zip | ae8198db5e8fc975b192fcb955ac6aa8 |
| 0021-0040.zip | 66e985d54ca3b6d32c139c2e3ccb530a |
| 0041-0060.zip | 2668fe23ac06e12455c49009502a9496 |
| 0061-0080.zip | e95924233c1811d728983514aa4627f8 |
| 0081-0100.zip | f77e122734e77b90ab41d0570a1e9b3f |

Extract them so that the case folders (`train_0001` … `train_0100`) sit directly in `D:\DEEP PSMA\DEEP PSMA EXTRACTED`.

## Software and resources

**Original environment:** Windows 64-bit.

**Python 3.11.15 (Anaconda):**
- NumPy 2.4.6, pandas 3.1.0, SciPy 1.17.1, SimpleITK 2.5.6 and matplotlib 3.11.2 (see `requirements.txt`);
- scikit-learn 1.9.1, which is used only by `A3_visualQC_unblinded_analysis.py`. This is the version installed in the original analysis environment; it is not recorded in the run manifests.

**MATLAB R2025a (25.1):** base MATLAB, plus the Image Processing Toolbox and the Statistics and Machine Learning Toolbox for some scripts.

**Exceptions recorded in the original manifests:**
- `ST2_overlap_audit.py` was run in a Linux environment with Python 3.10.12, NumPy 2.2.6 and pandas 2.3.3.
- The random QC draws were made with NumPy 2.5.3. The draws are not rerun; their records are provided.

**Storage:** a full reproduction writes about 2 GB of outputs, in addition to the extracted dataset.

**Runtime:** not benchmarked. The longest single run, the regional registration refinement, took about 6 hours on the original workstation.

## Reproducing the analyses

The archived scripts use the absolute paths of the original workstation, and they verify the MD5 checksums of their inputs. The runner uses the original folder layout and does not alter the archived scripts.

1. **Folder layout.**
   - Place the repository at `D:\DEEP PSMA\PAPER3_Discordance`.
   - Place the extracted dataset at `D:\DEEP PSMA\DEEP PSMA EXTRACTED`.
   - If there is no D: drive, `subst D: <folder>` creates one.
2. **Configuration.** In `reproduce_config.json`, enter the paths of the Python executable (an environment with the versions above) and of `matlab.exe`.
3. **Checks.** Run:
   - `python reproduce_runner.py --selftest`, which tests the runner's own renaming and verification logic;
   - `python reproduce_runner.py --dry-run`, which checks the folders, the 100 cases, both executables, the package versions, and each archived script against the script MD5 recorded in its original manifest.
4. **Run.**
   - `python reproduce_runner.py` runs all 66 runs in their original order.
   - `--only 1-3` runs a subset.
   - Runs that have already passed are skipped, so an interrupted reproduction can be resumed.

**How the runner checks the results.** After each run it renames the new outputs to the original run names (time stamps), so that later scripts find the file names and checksums they expect. It then compares every output file with the reference list:
- **Result files** (`.csv`, `.npz`, `.nii.gz`) must be byte-identical after any required line-ending conversion. Any other difference stops the runner.
- **Run manifests** record run times and therefore always differ. Only when every result file of a run is identical is the new manifest kept as `*_rerun.json` and the original manifest restored in its place, so that later scripts can verify it.
- **Line endings.** A few original runs (the Step 5 QC analyses and the ST2 overlap audit) were made in a Linux environment and wrote tables with LF line endings, whereas Python on Windows writes CRLF. A table that differs only in line endings is converted to LF. It is accepted only if the converted file then has exactly the original MD5. Such files are reported as `PASS_EOL`.
- **Logs, images and text reports** are compared and reported but never stop the runner.

The runner never overwrites existing files. Its reports are written to `reproduction/`: one row per file, a per-run summary with run times, and a log for each run.

**What has and has not been demonstrated.** The scripts, inputs and reference checksums have been audited for completeness, and the runner's logic has been tested on synthetic data.

On 26 November 2026, five short runs were reproduced byte for byte in a fresh copy of the repository, using the original software environment. They required only the dataset and published repository files. Together, they covered MATLAB image processing, Python image processing, synthetic validation and QC-record analysis. All 321 result files matched the original checksums; 4 required conversion from CRLF to LF line endings.

| Run | Language | Result files identical |
|---|---|---|
| Audit1D | MATLAB | 9/9 |
| P1 | Python (SimpleITK) | 304/304 |
| Tau validation | Python | 4/4 |
| Second-reader comparison | Python | 3/3, after LF line-ending conversion |
| Trunk decode | Python | 1/1, after LF line-ending conversion |

The reports and logs are in `reproduction_test_20261126/`. They include the first attempt at the second-reader run, which failed before line-ending conversion was added (runner v1.1).

The runner runs every script, including scripts originally run on Linux, using the configured Windows Python environment. The reproduction test confirmed this for the second-reader comparison and the trunk decode. The ST2 overlap audit, originally run with Python 3.10.12, has not yet been reproduced under the configured Windows Python environment.

The dry run passed on an identical copy of the final repository. The only differences were runner v1.0 rather than v1.1, and the README; neither affects the checks performed by the dry run.

The full 66-run reproduction has not yet been run, and the runner reports any difference.

## Run order

The runner reproduces 66 runs (`tools/reproduce_runs_v1_0.csv`). The steps below group them.

| Step | Script(s) | Purpose |
|---|---|---|
| 1 | `Audit1B_totalsegmentator_labels.m`, `Audit1D_lesion_mask_audit.m` | CT label verification; tumour-mask components and grid geometry |
| 2 | `Audit1C_c_rigid_refinement.py`, `Audit1C_e_translation_refinement.py` | Structure-level residuals; adopted trunk translation correction (T_A) |
| 3 | `Audit1C_f_vascular_muscular_QC.py`, `Amendment3_regional_rigid.py`, `SAP_R1_region_rules.py` | Region boundaries; regional head and pelvic refinement; region rules and plane table |
| 4 | `A3_visualQC_views_v1_2.py`, `Trunk_skeletal_QC_views_v1_1.py` | Rendering of the blinded CT-only QC pages that the readers graded. **Not rerun:** the QC records are provided instead (see below) |
| 5 | `A3_visualQC_unblinded_analysis.py`, `A3_second_reader_comparison.py`, `Trunk_QC_unblinding_decode.py`, `A3_intrareader_kappa_v1_0.m` | QC decisions and reader agreement, recomputed from the published QC records |
| 6 | `P1_native_fdg_components.py` | Native FDG lesion components |
| 7 | `P2_mapping_coverage_resampling_v1_1.py` → `P2b_export_voxelsets_csv.py` → `P3_matching_classification_v1_4.m` | Primary analysis: mapping, coverage, resampling, matching, classification, UFVF (DMTVF in the code) |
| 8 | `P3V_voxel_level_v1_0.m`, `Tau_validation.py` | Voxel-level analysis; matching validation on synthetic grids |
| 9 | `P3_matching_classification_v1_5.m` (reference, strict, permissive), `_v1_10.m` (overlap only) | Matching-rule sensitivity |
| 10 | P2 `_v1_1` (6-, 18-connectivity) → P2b `_v1_1` → P3 `_v1_6` | Connectivity sensitivity |
| 11 | P2 `_v1_2` → P2b `_v1_2` → P3 `_v1_7` (eligibility rule B, S1, S2, S4); P2 `_v1_3` → P2b `_v1_3` → P3 `_v1_8` (S3b, S3c) | Eligibility and registration-QC handling |
| 12 | `ST2a_organ_intersection.py` → P3 `_v1_9` → `ST2_overlap_audit.py`; `O1_organ_displacement.py` → `ST1_association_v1_0.m` | Visceral organs and organ displacement |
| 13 | `PERT_generate_v1_0.m` → P2 `_v1_6` → P2b `_v1_6` → P3 `_v1_13` | Geometric stress test (imposed rigid perturbations) |
| 14 | `PU_psma_uptake_v1_0.py`, `PX_crosstab_stability_uptake_v1_0.py` | PSMA uptake; stability × uptake (post-hoc) |
| 15 | P2 `_v1_7` → P2b `_v1_7` → P3 `_v1_14` | Registration-strategy comparator (post-hoc) |
| 16 | `AD_anatomical_distribution_v1_0.py` | Anatomical distribution and cross-compartment flags |
| 17 | `FIG_paper3_figures_v1_1.m`, `FIG4_anatomy_v1_0.m` | Figures 1–4 (greyscale, 1200 dpi) |
| 18 | `TAB_S7_software_verification_v1_2.m` | Online Resource Table S7. **Verification only:** it checks the original archived runs against their manifests and completion notes, and cannot be rerun from the repository alone because it needs the original output records |

**Step 4 is not rerun.**
- Rendering the QC pages is not needed for any later step, because the reader decisions it led to are provided.
- `A3_visualQC_views_v1_2.py` checks the output of `Amendment3b_femhead_evaluation.py` (femoral-head safety evaluation; not included). That script's table-producing version (MD5 8bbd017bac4a6306fcfd31adfbeb0ad5) was not retained; only its later report-only version (MD5 f6df42596ef28af909d90f00d45fa48c) was. This is documented in protocol section 15.3.15.
- Later scripts load functions from `A3_visualQC_views_v1_2.py`. These functions are not affected, because that script's input checks run only when Step 4 itself is run.

The endpoint is named UFVF_ref in the manuscript and DMTVF_ref in the code and locked plan; only the name differs.

## Notes for running scripts individually

1. **Paths.** Each script uses absolute paths (`D:\DEEP PSMA\DEEP PSMA EXTRACTED` for data, `D:\DEEP PSMA\PAPER3_Discordance\06_output` for outputs and `…\00_protocol` for specifications).
   - Most scripts define these paths themselves.
   - `A3_visualQC_views_v1_2.py`, the P2 scripts and `Trunk_skeletal_QC_views_v1_1.py` inherit them from `Amendment3_regional_rigid.py`.
   - The MATLAB scripts define them in a `cfg` structure.
2. **Input checksums.** Each script verifies the MD5 of its inputs before running and stops on any mismatch. Inputs include scripts, specifications and the time-stamped outputs of earlier stages.
   - Regenerated outputs get new time stamps, so outside the runner the file names and MD5 constants in later scripts would have to be changed by hand.
   - The runner avoids this by renaming outputs to their original names and verifying them.
3. **Blinding key.** `A3_visualQC_KEY_do_not_open_20261105.csv` is checked by several scripts and is provided, so no script needs to be edited.
4. **Self-tests and gates.** The analysis scripts (P1–P3, sensitivity analyses, stress test, PSMA uptake, anatomy and figures) contain synthetic self-tests. The sensitivity analyses also include a reproduction gate: it must reproduce the primary results exactly before any sensitivity result is written.

## Licence

This repository is provided for **non-commercial use only**, under the licences specified below. Commercial use requires separate permission from the relevant copyright holders.

- **Analysis code, runner and tools** (`02_code/`, `reproduce_runner.py`, `tools/`): PolyForm Noncommercial License 1.0.0 (see `LICENSE`; https://polyformproject.org/licenses/noncommercial/1.0.0).
- **Documents and records**: CC BY-NC 4.0 (see `LICENSE-DOCUMENTS.md`). This covers:
  - the protocol, statistical analysis plan and dated implementation specifications (`00_protocol/`);
  - the registration-QC records and notes (`06_output/`);
  - the reference files (`reference/`).

  These may be reused and adapted for non-commercial purposes, with attribution to the authors and citation of this repository and the manuscript.
- **DEEP-PSMA dataset:** not redistributed, and remains subject to its own licence (CC BY-NC 4.0).

For permissions beyond these licences, please contact the corresponding author.

## Citation

See `CITATION.cff`. Please also cite the DEEP-PSMA dataset. [Repository DOI to be added after release]

## Contact

Muhammad Aleem Khan, Department of Nuclear Medicine, Ibn Sina Hospital, Kuwait City, Kuwait — draleemkhan@yahoo.com
