"""
Trunk_QC_unblinding_decode.py - PAPER 3 trunk skeletal QC: decode the frozen blinded grades with the key (run once,
after Trunk_QC_grading_freeze_note_20261111.md). Reproduces Trunk_QC_decoded_grades_20261111.csv.
Usage: python Trunk_QC_unblinding_decode.py <checklist_pdftotext.txt> <Trunk_QC_KEY_do_not_open_20261109.csv> <out.csv>
  checklist_pdftotext.txt = `pdftotext -layout Trunk_Skeletal_QC_Blinded_Checklist_TR001-TR115.pdf`
                            (PDF MD5 d027cb90b3178a07c28f1561178ccdbd)
Final grades after reader-requested additional views (freeze note section 3): TR007 neck = 1, TR022 neck = 2.
"""
import re
import sys

import pandas as pd

FINAL = {("TR007", "neck"): ("1", "FINAL after additional views: definite cranio-caudal shift"),
         ("TR022", "neck"): ("2", "FINAL after additional views: mild cranio-caudal shift")}


def extract(txt):
    L = open(txt, encoding="utf-8").read().splitlines()
    rows, code = [], None
    for i, l in enumerate(L):
        m = re.match(r"^\s*(TR\d{3})?\s+(Abdomen|Neck|Thorax)\s+(\S+)\s*(.*)$", l)
        if m is None or m.group(3) == "+":
            continue
        code = m.group(1) or code
        reason = l[51:113].strip() or (L[i - 1][51:113].strip() + " " + L[i + 1][51:113].strip()).strip()
        rows.append(dict(code=code, subregion=m.group(2).lower(), grade=m.group(3), reason=reason,
                         flag_fields=" ".join(l[113:155].split())))
    i = [k for k, l in enumerate(L) if l.strip().startswith("Thorax +")][0]          # combined page (TR069)
    g = L[i + 1]
    rows.append(dict(code="TR069", subregion="thorax+abdomen", grade=g.split()[0], reason=g[51:113].strip(),
                     flag_fields=" ".join(g[113:155].split())))
    G = pd.DataFrame(rows)
    if len(G) != 314 or G.code.nunique() != 115:
        sys.exit(f"extraction incomplete: {len(G)} rows, {G.code.nunique()} codes")
    return G


def main(txt, key, out):
    G = extract(txt)
    G["grade_checklist"] = G["grade"]
    for (c, s), (gr, why) in FINAL.items():
        k = (G.code == c) & (G.subregion == s)
        G.loc[k, "grade"] = gr
        G.loc[k, "reason"] = G.loc[k, "reason"] + " | " + why
    K = pd.read_csv(key, keep_default_na=False)
    M = K.merge(G, on=["code", "subregion"], how="outer", indicator=True)
    if not (M["_merge"] == "both").all():
        sys.exit("key and checklist do not match one-to-one")
    M = M.rename(columns={"grade": "grade_final"})
    M = M[["code", "case_id", "subregion", "repeat_of", "grade_checklist", "grade_final", "reason", "flag_fields"]]
    M.sort_values(["code", "subregion"]).to_csv(out, index=False)


if __name__ == "__main__":
    main(*sys.argv[1:4])
