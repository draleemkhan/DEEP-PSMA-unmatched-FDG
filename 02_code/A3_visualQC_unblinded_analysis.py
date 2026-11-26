"""
A3_visualQC_unblinded_analysis.py - PAPER 3: unblinded analysis of the 100-patient visual registration QC.
Implements the LOCKED plan A3_visualQC_analysis_plan_20261108.md (MD5 d06b191c891b5e85d0e0d01fa9a61444).
Read-only on all inputs; writes tables to the output folder given as argument.
Usage: python A3_visualQC_unblinded_analysis.py <input_dir> <register_csv> <out_dir>
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXP = {"A3_visualQC_KEY100_do_not_open_20261106.csv": "aa72f6b91b2f7aa064eb6f2519033aa7",
       "A3_visualQC_grades_extracted_20261108.csv": "a0213a288301fdd9f24aa4cf13508df5",
       "A3_visualQC_analysis_plan_20261108.md": "d06b191c891b5e85d0e0d01fa9a61444",
       "Amendment3_fit_20261102_073040.csv": "dc74686ed0b024d8524668f9b9e16c94",
       "Phase1_QC_KEY_do_not_open_20261103_090259.csv": "e3f7d87ea76f8531a839bbcdea2ac26f",
       "A3_visualQC_KEY_do_not_open_20261105.csv": "c9a4548d3daf13d3e3dadc96171a9a94"}
FALLBACKS = {"train_0018", "train_0049", "train_0067"}
md5 = lambda p: hashlib.md5(Path(p).read_bytes()).hexdigest()


def g(v):
    return np.nan if str(v).strip().upper() in ("NA", "", "NAN") else float(v)


def main(ind, regcsv, outd):
    ind, outd = Path(ind), Path(outd)
    outd.mkdir(parents=True, exist_ok=True)
    for f, e in EXP.items():
        if md5(ind / f) != e:
            sys.exit(f"MD5 mismatch: {f}")
    K = pd.read_csv(ind / "A3_visualQC_KEY100_do_not_open_20261106.csv").fillna("")
    G = pd.read_csv(ind / "A3_visualQC_grades_extracted_20261108.csv", keep_default_na=False)
    F = pd.read_csv(ind / "Amendment3_fit_20261102_073040.csv")
    P1 = pd.read_csv(ind / "Phase1_QC_KEY_do_not_open_20261103_090259.csv").set_index("case_id")["code"]
    oldK = pd.read_csv(ind / "A3_visualQC_KEY_do_not_open_20261105.csv")
    REG = pd.read_csv(regcsv)
    cases = sorted(F["case_id"].unique())
    test_cases = [c for c in cases if c not in set(oldK["case_id"])][:2]

    # ---- unblind: one row per graded page ----
    kk = K.set_index("code")
    rows = []
    for _, r in G.iterrows():
        k = kk.loc[r["code"]]; reg = r["region"]
        xc = k[f"{reg}_X"]; yc = k[f"{reg}_Y"]
        gr = {xc: g(r["grade_X"]), yc: g(r["grade_Y"])}
        pref = {"X": xc, "Y": yc, "equivalent": "equivalent", "NA": "NA"}[r["preference"]]
        rows.append(dict(code=r["code"], case_id=k["case_id"], region=reg, category=k["category"],
                         repeat_of=k["repeat_of"], grade_A2=gr["A2"], grade_A3a=gr["A3a"], preference=pref,
                         C1C2_A2=r[f"C1C2_{'X' if xc == 'A2' else 'Y'}"] if reg == "head" else "",
                         C1C2_A3a=r[f"C1C2_{'X' if xc == 'A3a' else 'Y'}"] if reg == "head" else "",
                         jaw=r["jaw"], neck=r["neck"], skullbase=r["skullbase"], notes=r["row_text"]))
    U = pd.DataFrame(rows)
    fit = F.set_index(["case_id", "region"])
    U["adopted"] = [bool(fit.loc[(c, rg), "provisional_adopt"]) for c, rg in zip(U.case_id, U.region)]
    U["decision"] = [fit.loc[(c, rg), "decision"] for c, rg in zip(U.case_id, U.region)]
    U["rot_deg"] = [fit.loc[(c, rg), "rot_deg"] for c, rg in zip(U.case_id, U.region)]
    U["trans_mm"] = [fit.loc[(c, rg), "trans_mm"] for c, rg in zip(U.case_id, U.region)]
    U["cap"] = [f"{fit.loc[(c, rg), 'cap_deg']:.0f} deg / {fit.loc[(c, rg), 'cap_mm']:.0f} mm" for c, rg in zip(U.case_id, U.region)]
    U["status"] = np.where(U.adopted, "adopted",
                           np.where(U.decision.str.startswith("exceeds caps"), "cap-rejected", "insufficient improvement"))
    U["fallback_patient"] = U.case_id.isin(FALLBACKS)
    U["seen_before_grading"] = U.case_id.isin(test_cases)
    U["phase1_code"] = U.case_id.map(P1)
    U.to_csv(outd / "A3_unblinded_pages.csv", index=False)

    Pm = U[U.repeat_of == ""].copy()                        # primary pages decide

    # ---- acceptance rule (plan section 3) ----
    def accept(r):
        a2, a3 = r.grade_A2, r.grade_A3a
        if np.isnan(a2) and np.isnan(a3):
            return ("3a" if r.adopted else "A2"), "not visually verified (coverage)", False
        if r.adopted:
            if a3 >= 2:
                return "3a", "accepted", False
            return "A2", ("accepted fallback" if a2 >= 2 else "registration uncertainty"), True
        return "A2", ("accepted" if a2 >= 2 else "registration uncertainty"), False
    res = Pm.apply(accept, axis=1, result_type="expand"); res.columns = ["final_transform", "classification", "visual_veto"]
    Pm = pd.concat([Pm, res], axis=1)
    Pm.to_csv(outd / "A3_final_acceptance_per_region.csv", index=False)

    # ---- cap-rejected candidates (plan section 4) ----
    C = Pm[Pm.status == "cap-rejected"].copy()
    def cls(r):
        a2, a3 = r.grade_A2, r.grade_A3a
        if np.isnan(a2) or np.isnan(a3):
            return "not assessable"
        if a3 > a2 and r.preference == "A3a":
            return "beneficial" if a3 >= 2 else "better but not acceptable"
        if a3 == a2 and r.preference == "equivalent":
            return "equivalent"
        if a3 < a2:
            return "harmful"
        return "inconsistent grade/preference"
    C["classification"] = C.apply(cls, axis=1)
    C = C[["case_id", "code", "region", "category", "rot_deg", "trans_mm", "cap", "grade_A2", "grade_A3a", "preference",
           "classification", "seen_before_grading", "notes"]].sort_values(["region", "case_id"])
    C.to_csv(outd / "A3_cap_rejected_candidates.csv", index=False)

    # ---- summaries ----
    S = {}
    for reg in ("head", "pelvis"):
        d = Pm[Pm.region == reg]
        S[reg] = dict(
            n=len(d), adopted=int(d.adopted.sum()), cap_rejected=int((d.status == "cap-rejected").sum()),
            insufficient=int((d.status == "insufficient improvement").sum()),
            A2_grades=d.grade_A2.value_counts(dropna=False).sort_index().to_dict(),
            A3a_grades=d.grade_A3a.value_counts(dropna=False).sort_index().to_dict(),
            preference=d.preference.value_counts().to_dict(),
            adopted_3a_ge2=int(((d.adopted) & (d.grade_A3a >= 2)).sum()),
            adopted_3a_le1=int(((d.adopted) & (d.grade_A3a <= 1)).sum()),
            adopted_NA=int(((d.adopted) & d.grade_A3a.isna()).sum()),
            classification=d.classification.value_counts().to_dict(),
            final_transform=d.final_transform.value_counts().to_dict(),
            visual_vetoes=d[d.visual_veto][["case_id", "code", "grade_A2", "grade_A3a"]].to_dict("records"),
            uncertainty=d[d.classification == "registration uncertainty"][["case_id", "code", "status", "grade_A2", "grade_A3a"]].to_dict("records"))
        for grp, sel in (("targeted", d.category != "ROUTINE"), ("routine", d.category == "ROUTINE")):
            dd = d[sel]
            S[reg][f"{grp}_n"] = int(len(dd))
            S[reg][f"{grp}_classification"] = dd.classification.value_counts().to_dict()
            S[reg][f"{grp}_pref"] = dd.preference.value_counts().to_dict()
        ex = d[~d.seen_before_grading]
        S[reg]["sens_excl_TEST_classification"] = ex.classification.value_counts().to_dict()
        S[reg]["sens_excl_TEST_pref"] = ex.preference.value_counts().to_dict()
    # categories T3 / T4
    for cat in ("T3", "T4"):
        d = Pm[Pm.category.str.contains(cat)]
        S[cat] = d[["case_id", "region", "status", "grade_A2", "grade_A3a", "preference", "classification", "notes"]].to_dict("records")

    # ---- intra-reader agreement (plan section 6) ----
    R = U[U.repeat_of != ""]
    prim = U[U.repeat_of == ""].set_index(["code", "region"])
    ag = []
    for _, r in R.iterrows():
        p = prim.loc[(r.repeat_of, r.region)]
        ag.append(dict(repeat=r.code, source=r.repeat_of, region=r.region, case_id=r.case_id,
                       A2_first=p.grade_A2, A2_repeat=r.grade_A2, A3a_first=p.grade_A3a, A3a_repeat=r.grade_A3a,
                       pref_first=p.preference, pref_repeat=r.preference))
    A = pd.DataFrame(ag)
    A.to_csv(outd / "A3_repeat_pairs.csv", index=False)

    def agree(d):
        out = {}
        for m in ("A2", "A3a"):
            a, b = d[f"{m}_first"], d[f"{m}_repeat"]
            ok = a.notna() & b.notna()
            out[f"{m}_exact"] = f"{int((a[ok] == b[ok]).sum())}/{int(ok.sum())}"
            out[f"{m}_within1"] = f"{int(((a[ok] - b[ok]).abs() <= 1).sum())}/{int(ok.sum())}"
        out["pref_exact"] = f"{int((d.pref_first == d.pref_repeat).sum())}/{len(d)}"
        return out
    S["agreement"] = {reg: agree(A[A.region == reg]) for reg in ("head", "pelvis")}
    try:
        from sklearn.metrics import cohen_kappa_score
        pv = A[A.region == "pelvis"]
        kap = {}
        for m in ("A2", "A3a"):
            ok = pv[f"{m}_first"].notna() & pv[f"{m}_repeat"].notna()
            a, b = pv.loc[ok, f"{m}_first"], pv.loc[ok, f"{m}_repeat"]
            kap[f"{m}_weighted_kappa"] = (float(cohen_kappa_score(a, b, weights="linear", labels=[0, 1, 2, 3]))
                                          if len(set(a) | set(b)) > 1 else "not estimable (single category)")
        kap["pref_kappa"] = (float(cohen_kappa_score(pv.pref_first, pv.pref_repeat))
                             if len(set(pv.pref_first) | set(pv.pref_repeat)) > 1 else "not estimable")
        S["agreement"]["pelvis_kappa_secondary"] = kap
    except ImportError:
        S["agreement"]["pelvis_kappa_secondary"] = "scikit-learn not installed"
    S["head_repeat_pairs"] = A[A.region == "head"].to_dict("records")

    # ---- mapping to Phase 1 register (context only) ----
    REG["phase1_code"] = REG["code"]
    M = Pm.merge(REG.groupby("phase1_code").agg(
        p1_Frag=("Frag", lambda s: "/".join(s.astype(str))), p1_Mtl=("Mtl", lambda s: "/".join(s.astype(str))),
        p1_FHL=("FHL", lambda s: "/".join(s.astype(str))), p1_FHR=("FHR", lambda s: "/".join(s.astype(str))),
        p1_CH=("CH", "first"), p1_notes=("note", lambda s: " || ".join(str(x) for x in s.tolist() if isinstance(x, str) and x.strip()))),
        left_on="phase1_code", right_index=True, how="left")
    M.to_csv(outd / "A3_final_acceptance_with_phase1_context.csv", index=False)
    S["test_cases_seen_before_grading"] = test_cases
    S["inputs_md5"] = EXP
    (outd / "A3_unblinded_summary.json").write_text(json.dumps(S, indent=1, default=str), encoding="utf-8")
    print(json.dumps({k: v for k, v in S.items() if k in ("head", "pelvis", "agreement")}, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:4])
