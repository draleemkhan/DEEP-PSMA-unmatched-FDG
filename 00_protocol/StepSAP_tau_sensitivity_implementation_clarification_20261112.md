# Dated implementation clarification: τ sensitivity analyses (T-strict, T-perm)
**Date:** 12 November 2026, 00:38 (Asia/Kuwait, UTC+03:00)
**Relates to:** consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) §7 (T-strict, T-perm); τ specification LOCKED (MD5 285debf473d92a32c765f6f7ee45d412) §τ6; primary P3 run `P3_conn26_20261112_001915` (completion note MD5 c02b3a52d0b172d2bfdafcd30d6401c8).
**Status:** written after the primary analysis and the voxel-level sensitivity analysis, **before** any τ sensitivity result is computed.

## The question
τ_i enters two locked rules: (1) **matching** — matched if d ≤ τ_i (τ spec D3); (2) the **coverage margin** — V_I (coverage) if any mapped voxel centre lies within τ_i of the PSMA PET field-of-view boundary (SAP item 3). The SAP defines the τ sensitivity analyses as alternative tolerances (strict = max voxel dimension; permissive = 2 × max voxel dimension) without stating whether the coverage margin follows them.

## Author decision (12 Nov 2026, 00:36): Option (a) — matching τ only
- **Varied:** the matching tolerance only: τ_strict,i = max(Δx_i, Δy_i, Δz_i); τ_i (reference) = voxel diagonal; τ_perm,i = 2 × max(Δx_i, Δy_i, Δz_i). This applies to every matching decision, including both categories in the boundary-band comparison.
- **Fixed:** the coverage margin stays at the locked reference τ_i (P2 flags unchanged); eligibility, regions, registration QC, band candidacy, resampling, precedence and native-volume accounting unchanged.
- **Interpretation:** these analyses measure sensitivity to the **matching rule**, not the sensitivity of the entire pipeline to τ (which would also vary the coverage margin; Option (b), not required for the current objective).
- **Gate:** before either sensitivity result is accepted, the same script run with the reference τ must reproduce the primary P3 component and patient tables exactly (statuses, causes, categories, distances, volumes).
- Nested tolerances are reported as sensitivity, not robustness (protocol §5).
