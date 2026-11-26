# A3 visual registration QC — pre-unblinding analysis plan and acceptance rules
**Date:** 8 November 2026 (Asia/Kuwait, UTC+03:00)
**Relates to:** spec v1.0 (MD5 10bea05a1792ffdc4401f5825b322c15), amendment v1.1 (MD5 089e428056faeea5345a2f8602658166), note on VR009/C4 and NA (MD5 e25538129147e192c180af46b9c987bb).
**Status:** LOCKED (see §8). Written after blinded grading was complete and **before the key was opened**.

## 1. Frozen inputs
- **Graded record:** `A3_Visual_Registration_QC_Compact_120VR.pdf`, MD5 **15935a0f514768fa46c71a41d7f3699f** (saved 8 Nov 2026, 12:27).
  - 220 rows (VR001–VR100 head and pelvis; VR101–VR120 one region each).
  - Six head pages recorded as NA (insufficient coverage): VR007, VR011, VR032, VR038, VR041, VR047.
  - Corrected before freezing from an earlier version (MD5 93e7544f7361fe5d530eb4fdb65fbab4), in which those six were recorded as 0/0. The VR009 comment wrongly referred to the C4 exception and was removed. No other grade changed.
- **Key:** `A3_visualQC_KEY100_do_not_open_20261106.csv`, MD5 aa72f6b91b2f7aa064eb6f2519033aa7, **not yet opened**.
- Grades are extracted from the PDF into a table; the author confirms that the extraction matches the PDF before unblinding.

## 2. Two separate outcomes
1. **Algorithmic adoption:** whether Amendment 3a met the predefined numerical criteria (pilot run 20261102_073040). Unchanged by this review.
2. **Visual acceptability:** whether the registration finally retained for each patient-region is anatomically acceptable.

These are reported separately. A cap-rejected candidate that looks visually better is **not** counted as an adopted correction. An adopted correction graded ≤ 1 is a **visual failure despite numerical adoption**.

## 3. Final acceptance rule (per patient-region; primary pages VR001–VR100 only)
Acceptance threshold: **grade ≥ 2** (2 = minor residual mismatch without meaningful misregistration).

| Situation | Final transform | Classification |
|---|---|---|
| 3a adopted, 3a grade ≥ 2 | 3a | Accepted |
| 3a adopted, 3a grade ≤ 1 | Amendment 2 — **visual veto, recorded as such** | Accepted fallback if A2 ≥ 2; **registration uncertainty** if A2 ≤ 1 |
| 3a not adopted (cap, insufficient improvement), A2 ≥ 2 | Amendment 2 | Accepted |
| 3a not adopted, A2 ≤ 1 | Amendment 2 | **Registration uncertainty** (kept, not excluded) |
| NA (insufficient coverage) | As decided by the algorithm | **Not visually verified (coverage)**; kept and flagged |

- "Amendment 2" means the supplied transform for the fallback patients train_0018, 0049 and 0067.
- The visual review can reject an adopted correction; it **cannot adopt a rejected one**.
- Phase 1 segmentation QC findings (prosthesis, truncation, implants, femoral-head ratings, pre-screen flags) are mapped to each region as **context**. They do not override the visual grade.

## 4. Cap-rejected candidates: descriptive, region-separate (no decision threshold)
Skull (8 candidates) and pelvis (2 candidates) are analysed **separately and not pooled**. For each candidate the following is reported:
- the rotation and translation size, and which cap was exceeded;
- the candidate grade, the Amendment 2 grade and the recorded preference;
- a classification, using **both grades and preference**:
  - **Beneficial:** candidate grade ≥ 2 **and** strictly better than A2 (higher grade, consistent with the preference). Being acceptable alone is not enough.
  - **Better but not acceptable:** strictly better than A2, but candidate grade ≤ 1.
  - **Equivalent:** equal grades and preference "equivalent".
  - **Harmful:** candidate strictly worse than A2.
  - **Not assessable:** NA.

No cap is changed. All 10 pilot rejections stand. Any relaxation of the caps is **exploratory** and may only be evaluated prospectively, through a separate algorithmic rerun with a safety assessment and a dated amendment.

## 5. Other reporting (descriptive)
- **Targeted subgroup** (20 patients; patients in several categories counted once, with all categories listed) vs **routine group** (80 patients).
- **By adoption status** per region: adopted; cap-rejected; insufficient improvement; fallback patients.
- **Direction-disagreement cases** (T3) and **Amendment 3b flagged sides** (T4, prosthetic vs native).
- A2 and 3a grade distributions and preferences, by region.

## 6. Intra-reader agreement (VR101–VR120 vs their primary pages)
- **Primary measures:** exact grade agreement (%) and agreement within one grade (%), for the X grade, the Y grade and the preference, reported **separately for head and pelvis**.
- **Secondary:** weighted kappa for grades and kappa for preference, **pelvis only** (n = 17). Reported with the caveat that kappa is unstable when one category dominates (the prevalence paradox).
- **Head (n = 3):** no kappa. The three pairs are listed individually.
- Limitation reported: repeats are recognisable as repeats by code range (VR101–VR120), though not by source page; and the head/pelvis split (3/17) came from an unstratified random draw.

## 7. Sensitivity and limitations
- **TEST001/TEST002** (seen in preflight before grading): grades retained and flagged. The main results are repeated with these two patients excluded.
- Limitations reported: single reader; partial blinding where a very large correction is recognisable on sight; level 3 (suitability for lesion correspondence) not assessed.
- Any analysis not specified here is labelled exploratory.

## 8. Lock record
- Approved by the author on 8 November 2026, 12:32, including the beneficial-classification refinement (§4).
- Locked on: 8 November 2026, 12:35 (Asia/Kuwait, UTC+03:00), before the key was opened.
- This plan's MD5: recorded externally.
