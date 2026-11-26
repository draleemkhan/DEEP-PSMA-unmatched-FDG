# Dated implementation clarification: sensitivity analyses S3b and S3c
**Date:** 12 November 2026, 01:15 (Asia/Kuwait, UTC+03:00)
**Relates to:** Step 2 SAP Section R LOCKED (`StepSAP_section_R_registration_handling_20261109.md`, MD5 eecce0280370aca026eb1b2edf922efe) §R.5 and §R.6 (S3); Section R.1 LOCKED (MD5 8fe2d04b39a906594d0c3370c2859cb5) §R1.6 (S3 applies to all failed regions, necks included); consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) §7; primary P3 run `P3_conn26_20261112_001915`.
**Status:** written after the primary analysis and the sensitivity analyses T-strict, T-perm, C6, C18, EB, S1, S2 and S4, **before** any S3 result is computed. The locked S3 definitions are not changed; this note fixes implementation details that the SAP does not state.

## Locked definitions (Section R.6, verbatim in substance)
- **S3b, volume-weighted:** V_I only if > 50 % of the component's mapped points (under T_A) lie in a failed region.
- **S3c, split:** the points in a failed region go to V_I and the remainder is assessed.
- Failed regions: the 29 locked regions in 23 patients (7 grade-1 heads, 6 NA heads, 2 grade-1 pelvises, 14 failed necks); region membership by the locked planes and tie rules (`SAP_R1_region_rules.region_of`).

## Author decisions (12 Nov 2026, 01:10)
1. **Points and weighting.** The points are the centres of all native FDG voxels of the eligible component, each weighted equally (equal native voxel volume, so the point fraction equals the volume fraction). In S3b, a component with ≤ 50 % of its points in a failed region is not V_I (registration); the whole component proceeds to the band, coverage and matching steps unchanged.
2. **Transform for the failed-region test.** **T_A only**, in both S3b and S3c (not "any relevant transform" as in the primary rule R.5). This keeps the point classification identical in S3b and S3c, so that any difference between them reflects only the handling rule.
3. **S3c remainder.** The R6 transform assignment (final transform), the relevant transforms and the band candidacy (including the alternative transform) are those of the **original whole component** and are not recomputed. Coverage (FOV and τ margin, every transform used), resampling, the zero-survivor check, both band categories and matching are recomputed on the **remainder's voxels only**.
4. **Eligibility.** No new eligibility test for the remainder (R.2: eligibility is native, fixed before any registration-dependent step; V_S is never reassigned).
5. **Accounting.** Each part's volume is computed from the **original component's native voxel volume**: V_failed part = n_failed × (native component volume / native voxel count); V_remainder = native component volume − V_failed part. **Resampled voxel counts are never used for volumes.** The failed part goes to V_I (registration); the remainder takes its own status (V_I band, V_I coverage, V_M or V_U). V_FDG,total (the denominator) is unchanged.
6. **Safeguard.** If the remainder has no surviving voxel after resampling (final mapping), it is V_I (coverage) under the existing zero-survivor rule.
7. **Edge cases.** All points in a failed region → whole component V_I (registration), no remainder. No point in a failed region → component unchanged by S3c (proceeds as in the primary analysis, apart from rule 2: points mapped into a failed region by a transform other than T_A no longer make it V_I).

## Implementation and verification
- P2 v1.3 (`--scenario S3b | S3c`): per-component T_A point count in failed regions and fraction; S3c builds a native label image with the failed points removed from split components only, and resamples it (nearest neighbour; the other components' voxel sets are unaffected). Synthetic self-tests added for every rule above.
- P2b v1.3: accepts P2 v1.3 runs (conversion unchanged).
- P3 v1.8 (`'scenario'` with S3b or S3c): splits the native volume per rule 5; same gate as v1.6/v1.7 (the generic loader must reproduce the primary tables before any S3 result is written); reconciliation and V_FDG,total identity checks.
- Workflow: self-tests → two-patient preflight (P2) → full cohort → P2b → P3.

## Interpretation (author)
S3b tests a **majority-overlap** rule; S3c tests **partial-volume attribution**. Both are sensitivity analyses of the handling of failed registration regions; they do not replace the primary result.
