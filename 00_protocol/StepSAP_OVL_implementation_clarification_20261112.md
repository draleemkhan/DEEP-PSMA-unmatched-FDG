# Dated implementation clarification: OVL overlap-only agreement analysis
**Date:** 12 November 2026, 07:45 (Asia/Kuwait, UTC+03:00)
**Relates to:** protocol v2.10 (MD5 2be3d63a9f130dbc83b7b03e11c37b70) §5 secondary item 5 ("alternative matching method (overlap-based) agreement (% consistent, kappa, burden difference, patient status change — agreement ≠ accuracy)"); consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) §7 (OVL); τ specification LOCKED (MD5 285debf473d92a32c765f6f7ee45d412) §τ6 note; τ sensitivity clarification (MD5 f78a14882975a2f1bef79adcba5d1618, Option (a)); primary P3 run `P3_conn26_20261112_001915`; T-strict run `P3T_strict_conn26_20261112_003925`.
**Status:** written after the primary analysis and the sensitivity analyses (T-strict, T-perm, C6, C18, EB, S1, S2, S3b, S3c, S4, ST2), **before** any OVL result is computed. The locked SAP is not changed; this note fixes implementation details that it does not state. Author decisions O1–O6 and two reporting safeguards, 12 Nov 2026, 07:36.

## O1 Reclassification
- Matching tolerance **τ = 0** in every matching decision, including both categories of the boundary-band comparison. Since d = 0 only when the FDG component shares at least one voxel with the PSMA TTB mask (otherwise d ≥ the smallest voxel spacing > 0), a component is matched under OVL **only if at least one mapped FDG voxel overlaps a PSMA tumour-mask voxel**.
- Unchanged: eligibility (26-connectivity, rule A), registration step, coverage (margin fixed at the **reference τ**, as in τ Option (a)), resampling, precedence, native-volume accounting, denominator, component identifiers. V_I (band) may differ because both band categories use τ = 0.
- Gate: before any OVL result is written, the reference analysis recomputed by the same script must reproduce the primary P3 tables.

## O2 Agreement population (2 × 2)
- Units: eligible components that are **V_M or V_U under both methods**. Reported: number of eligible components; number V_I under both; number with a transition involving V_I (V_I under one method only), **excluded** from the 2 × 2 table; number assessable under both; the four 2 × 2 counts (M/M, M/U, U/M, U/U).
- % consistent = (M/M + U/U) / assessable; Cohen's κ = (p_o − p_e) / (1 − p_e) from the 2 × 2 counts. Component-level, unweighted.
- Comparisons: **OVL vs reference τ** (primary) and **OVL vs T-strict** (stored, MD5-checked T-strict tables).

## O3 Uncertainty of κ
- Patient-cluster bootstrap (all 100 patients resampled with replacement, each with all its assessable components), **2,000 resamples, percentile 95 % CI**, MATLAB `rng(20261114, 'twister')`, reset before each comparison.
- **Undefined κ (safeguard):** a replicate is invalid if it has no assessable component or 1 − p_e = 0 (no variation); it is recorded as invalid, never assigned a value. The number of valid replicates is reported.
- **CI estimability rule** (same pragmatic threshold as the locked ST1B fallback, SP5): the CI is reported only if ≥ 90 % of the 2,000 replicates are valid; otherwise the point estimate is reported with the statement that the CI was not estimable and why. If the point estimate itself is undefined, this is reported.

## O4 Patient status
Status = **DMTVF_ref > 0** (the dichotomy reported in protocol §5). Reported: the 2 × 2 patient table and the number of patients changing status; the > 1 / 5 / 10 % counts descriptively.

## O5 Burden difference
Per patient: DMTVF_ref (OVL) − DMTVF_ref (comparator), median and IQR; cohort V_U per method; V_I by cause per method.

## O6 Implementation
P3 v1.10 (MATLAB, R.8): τ mode `ovl` (τ = 0) and all agreement statistics; self-tests (face contact unmatched and shared voxel matched at τ = 0; κ on known tables; undefined κ flagged invalid); no P2 / P2b rerun.

## Interpretation
**Agreement ≠ accuracy.** High agreement between reference-τ and overlap-only matching would show methodological consistency, not biological validity.
