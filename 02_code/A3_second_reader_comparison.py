"""
A3_second_reader_comparison.py - PAPER 3: second-reader vs first-reader comparison (locked plan
A3_visualQC_second_reader_plan_20261108.md, MD5 ff33c80f6b77dd758fee432a0f525f45). First reader = primary result.
Usage: python A3_second_reader_comparison.py <sr_text_from_pdf> <sr_key> <first_reader_acceptance_csv> <out_dir>
  sr_text_from_pdf: `pdftotext -layout A3_Second_Reader_Blinded_Checklist_SR01-SR40.pdf` output
                    (PDF MD5 a40321a384ab3d6ad25c36c7ba62b77c)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PREF = {"x": "X", "y": "Y", "eq": "equivalent", "eqv": "equivalent", "na": "NA"}


def g(v):
    return np.nan if str(v).strip().upper() in ("NA", "", "NAN") else float(v)


def extract(txt):
    L = Path(txt).read_text().splitlines()
    codes = [(i, m.group(1)) for i, l in enumerate(L) for m in [re.search(r"\b(SR\d{2})\b", l)] if m]
    rows = []
    for i, l in enumerate(L):
        m = re.search(r"\b(Head|Pelvis)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(.*))?$", l)
        if not m:
            continue
        c = min(codes, key=lambda t: abs(t[0] - i))[1]
        rest = m.group(5) or ""
        sf = re.findall(r"(?:^|\s{2,})([a-z])(?=\s{2,}|\s*$)", rest) if m.group(1) == "Head" else []
        sf = sf[-5:] if len(sf) >= 5 else [""] * 5
        rows.append(dict(code=c, region=m.group(1).lower(), grade_X=m.group(2).upper(), grade_Y=m.group(3).upper(),
                         preference=PREF[m.group(4).lower()], C1C2_X=sf[0], C1C2_Y=sf[1], jaw=sf[2], neck=sf[3],
                         skullbase=sf[4], row_text=re.sub(r"\s{2,}", " | ", rest.strip())))
    G = pd.DataFrame(rows).sort_values(["code", "region"])
    if len(G) != 80 or G["code"].nunique() != 40:
        sys.exit(f"extraction incomplete: {len(G)} rows")
    return G


def acceptance(adopted, a2, a3):
    if np.isnan(a2) and np.isnan(a3):
        return ("3a" if adopted else "A2"), "not visually verified (coverage)"
    if adopted:
        return ("3a", "accepted") if a3 >= 2 else ("A2", "accepted fallback" if a2 >= 2 else "registration uncertainty")
    return "A2", ("accepted" if a2 >= 2 else "registration uncertainty")


def capcls(a2, a3, p):
    if np.isnan(a2) or np.isnan(a3):
        return "not assessable"
    if a3 > a2 and p == "A3a":
        return "beneficial" if a3 >= 2 else "better but not acceptable"
    if a3 == a2 and p == "equivalent":
        return "equivalent"
    if a3 < a2:
        return "harmful"
    return "inconsistent grade/preference"


def main(txt, key, r1csv, outd):
    outd = Path(outd); outd.mkdir(parents=True, exist_ok=True)
    G = extract(txt)
    G.to_csv(outd / "A3_second_reader_grades_extracted_20261109.csv", index=False)
    K = pd.read_csv(key).set_index("code")
    out = []
    for _, r in G.iterrows():
        k = K.loc[r["code"]]; xc, yc = k[f"{r['region']}_X"], k[f"{r['region']}_Y"]
        gr = {xc: g(r["grade_X"]), yc: g(r["grade_Y"])}
        out.append(dict(sr_code=r["code"], case_id=k["case_id"], stratum=k["stratum"], region=r["region"],
                        R2_A2=gr["A2"], R2_A3a=gr["A3a"],
                        R2_pref={"X": xc, "Y": yc, "equivalent": "equivalent", "NA": "NA"}[r["preference"]],
                        R2_notes=r["row_text"]))
    R2 = pd.DataFrame(out)
    R1 = pd.read_csv(r1csv, keep_default_na=False, na_values=[""])
    for c in ("grade_A2", "grade_A3a"):
        R1[c] = R1[c].map(g)
    R1["preference"] = R1["preference"].replace({"nan": "NA"}).fillna("NA")
    R1["adopted"] = R1["adopted"].astype(str) == "True"
    M = R2.merge(R1[["case_id", "region", "code", "status", "adopted", "grade_A2", "grade_A3a", "preference",
                     "classification", "final_transform"]], on=["case_id", "region"], how="left")
    M = M.rename(columns={"code": "vr_code", "grade_A2": "R1_A2", "grade_A3a": "R1_A3a", "preference": "R1_pref",
                          "classification": "R1_class", "final_transform": "R1_final"})
    M[["R2_final", "R2_class"]] = [acceptance(a, x, y) for a, x, y in zip(M.adopted, M.R2_A2, M.R2_A3a)]
    M["R1_cap"] = [capcls(a, b, p) if s == "cap-rejected" else "" for s, a, b, p in zip(M.status, M.R1_A2, M.R1_A3a, M.R1_pref)]
    M["R2_cap"] = [capcls(a, b, p) if s == "cap-rejected" else "" for s, a, b, p in zip(M.status, M.R2_A2, M.R2_A3a, M.R2_pref)]
    M.to_csv(outd / "A3_second_reader_vs_first_reader_20261109.csv", index=False)

    def agree(d):
        o = {}
        for m in ("A2", "A3a"):
            a, b = d[f"R1_{m}"], d[f"R2_{m}"]; ok = a.notna() & b.notna()
            o[f"{m}_exact"] = f"{int((a[ok] == b[ok]).sum())}/{int(ok.sum())}"
            o[f"{m}_within1"] = f"{int(((a[ok] - b[ok]).abs() <= 1).sum())}/{int(ok.sum())}"
        o["preference"] = f"{int((d.R1_pref == d.R2_pref).sum())}/{len(d)}"
        o["NA"] = f"{int((d.R1_A2.isna() == d.R2_A2.isna()).sum())}/{len(d)}"
        o["acceptance_class"] = f"{int((d.R1_class == d.R2_class).sum())}/{len(d)}"
        return o
    S = {f"{st}_{rg}": agree(M[(M.stratum == st) & (M.region == rg)]) for st in ("critical", "random") for rg in ("head", "pelvis")}
    S["all"] = agree(M)
    c = M[M.status == "cap-rejected"]
    S["cap_classification_same"] = f"{int((c.R1_cap == c.R2_cap).sum())}/{len(c)}"
    pd.DataFrame(S).T.to_csv(outd / "A3_second_reader_agreement_summary_20261109.csv")
    print(pd.DataFrame(S).T.to_string())


if __name__ == "__main__":
    main(*sys.argv[1:5])
