# Dated addendum: PSMA uptake implementation (SAP §8.2) — erosion, blood pool ≥ liver, precedence
**Date:** 12 November 2026, 11:50 (Asia/Kuwait, UTC+03:00)
**Supplements:** `00_protocol\StepSAP_PSMAuptake_implementation_clarification_20261112.md` (MD5 5d8deb6bce244aa480ce2aea0866a727), kept unchanged. Approved by the author 12 Nov 2026, 11:43. Written before coding; no PSMA uptake value, reference value or class has been computed or viewed.

## 1. Erosion (corrects the justification in §3, step 3)
3 × 3 × 3 structuring element (26-neighbour binary erosion), 1 iteration, on the PET grid; voxels outside the image are background. **Justification:** a reproducible choice for **reference-mask erosion** (removing the mask rim to limit partial-volume and spill-in contamination of the reference SUVmean). It is a separate operation from FDG component connectivity; the earlier reference to "consistency with conn26" is withdrawn as a justification.

## 2. Blood pool ≥ liver (replaces the second bullet of §5)
Applied literally, the ordinary rules overlap when blood pool ≥ liver (a value with liver < SUVmax ≤ blood pool would be both "absent" and "positive"). To keep the categories **mutually exclusive**, in a patient with blood pool ≥ liver:
- **Absent:** SUVmax ≤ blood pool.
- **Positive:** SUVmax > blood pool.
- **Low:** unavailable (empty).

**This changes the usual positive-category definition (> liver) for these patients only.** Such patients are flagged; the number of affected patients and their components are reported. In all other patients the ordinary locked rules apply unchanged (absent ≤ blood pool; low > blood pool and ≤ liver; positive > liver; equality to the lower class).

## 3. Precedence of "not classifiable" categories (makes §3 and §4 explicit)
For each of the 217 primary V_U components, in this order:
1. Maximum extent ≤ 10 mm → **not classifiable (partial volume)**.
2. Otherwise, either reference unavailable in that patient → **not classifiable (reference unavailable)**. The comparison with the available reference is listed per component for transparency but is **not used as a class**.
3. Otherwise → absent / low / positive by §2 or the ordinary rules.
SUVmax is recorded for every component, and every component enters the SUVmax ≥ 3 count (U6).

## 4. Note on an existing column
The P1 table column `max_extent_mm` is a world-axis bounding-box extent including voxel width (Audit 1D convention). It is **not** the approved U4 definition and is not used; extent is computed afresh as the maximum Euclidean distance between voxel centres.
