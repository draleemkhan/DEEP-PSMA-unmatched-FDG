"""
PX_crosstab_stability_uptake_v1_0.py - PAPER 3: exploratory cross-tabulation of perturbation stability and PSMA uptake
class for the 217 components classified V_U in the primary analysis.

POST-HOC and EXPLORATORY: 00_protocol\\PAPER3_exploratory_crosstab_note_20261112.md (MD5 49b96844777d91a8f967256928d638b6),
documented 12 Nov 2026 09:40 before any full-level perturbation or uptake result. Joins existing outputs only:
no new thresholds, exclusions or biological categories; no new measurement.

Inputs: PSMA uptake run PU_20261112_115324; perturbation P3P_L1_20261112_101030 and P3P_L2_20261112_103235
(component tables). Levels are NEVER combined; no composite score.
Per level, per component: N_U + N_M + N_I = 50, N_I = N_I,band + N_I,coverage (registration-QC V_I and V_S must be 0
for primary V_U components, since registration-QC status is kept from the primary analysis; checked).
Stability (note section 3): stably unmatched = N_U >= 48; unstable = N_U < 48 with NON-EXCLUSIVE flags
  sometimes matched (N_M > 0) and sometimes indeterminate (N_I > 0). For presentation the unstable components are
  also shown in the three disjoint flag patterns (matched only / indeterminate only / both); the two flag totals are
  reported separately as well. Stable components' N_M and N_I (<= 2) are reported, not used to reclassify.
Uptake class: as in PU (absent, low, positive, NC_partial_volume, NC_reference_unavailable).
Every cell: components, native FDG mL, patients. Indeterminate causes reported separately (band, coverage).
Outputs: 06_output\\PX_<run>\\ PX_components_<run>.csv (joined table), PX_crosstab_<run>.csv (long format),
  PX_summary_<run>.txt, PX_manifest_<run>.json.
Usage (folder 02_code, env totalseg):
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u PX_crosstab_stability_uptake_v1_0.py --selftest
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u PX_crosstab_stability_uptake_v1_0.py
"""
import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
PU_RUN = "PU_20261112_115324"
LEVELS = {"L1": "P3P_L1_20261112_101030", "L2": "P3P_L2_20261112_103235"}
INPUTS = {
    (OUT_DIR, f"{PU_RUN}/PU_components_20261112_115324.csv"): "790db257935376785e5753ddd136d0eb",
    (OUT_DIR, f"{PU_RUN}/PU_manifest_20261112_115324.json"): "a7f16d8d9f8fdd0bedf5a4adfc9b4d06",
    (OUT_DIR, "P3P_L1_20261112_101030/P3P_L1_components_20261112_101030.csv"): "6d93c52c4617dd1ea44600e2f29431f3",
    (OUT_DIR, "P3P_L2_20261112_103235/P3P_L2_components_20261112_103235.csv"): "9df0305d6dfb6e747ccfe7f285838763",
    (PROTO_DIR, "PAPER3_exploratory_crosstab_note_20261112.md"): "49b96844777d91a8f967256928d638b6",
}
NPERT = 50
STABLE_MIN = 48
UPTAKE = ["absent", "low", "positive", "NC_partial_volume", "NC_reference_unavailable"]
PATTERN = ["stably unmatched", "unstable: sometimes matched only", "unstable: sometimes indeterminate only",
           "unstable: both flags"]
FLAGS = ["unstable, flag sometimes matched (any)", "unstable, flag sometimes indeterminate (any)"]
INTERP = {
    ("stably unmatched", "absent"): "imaging pattern most consistent with FDG+/PSMA-low disease (not confirmed)",
    ("stably unmatched", "low"): "persistent spatial discordance with measurable low PSMA uptake",
    ("stably unmatched", "positive"): "persistent annotation discordance despite appreciable PSMA signal "
                                      "(unannotated tumour, physiological uptake, spill-in, registration error or another lesion)",
    "matched": "correspondence depends on spatial alignment",
    "indeterminate": "assessability depends on spatial alignment",
    "NC": "uptake-based interpretation withheld",
}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ pure helpers (tested) ------------------------------
def level_counts(L):
    """L: P3P component table rows of primary V_U components. Returns N_U, N_M, N_I_band, N_I_cov, N_I; checks."""
    if (L.n_V_I_registration != 0).any() or (L.n_V_S != 0).any():
        raise RuntimeError("registration-QC V_I or V_S > 0 for a primary V_U component")
    N = pd.DataFrame(dict(N_U=L.n_V_U.astype(int), N_M=L.n_V_M.astype(int), N_I_band=L.n_V_I_band.astype(int),
                          N_I_cov=L.n_V_I_coverage.astype(int)), index=L.index)
    N["N_I"] = N.N_I_band + N.N_I_cov
    if not (N.N_U + N.N_M + N.N_I == NPERT).all():
        raise RuntimeError("N_U + N_M + N_I != 50")
    return N


def pattern(nu, nm, ni):
    if nu >= STABLE_MIN:
        return PATTERN[0]
    if nm > 0 and ni > 0:
        return PATTERN[3]
    if nm > 0:
        return PATTERN[1]
    if ni > 0:
        return PATTERN[2]
    raise RuntimeError("unstable component without N_M or N_I")   # impossible when the sum is 50


def cell(df):
    return dict(components=int(len(df)), native_mL=float(df.volume_native_mL.sum()), patients=int(df.case_id.nunique()))


def crosstab(J, lv):
    """long-format table for one level: rows = patterns + flag totals + total; columns = uptake classes + total."""
    out = []
    rows = [(p, J[f"{lv}_pattern"] == p) for p in PATTERN]
    rows += [(FLAGS[0], (J[f"{lv}_N_U"] < STABLE_MIN) & (J[f"{lv}_N_M"] > 0)),
             (FLAGS[1], (J[f"{lv}_N_U"] < STABLE_MIN) & (J[f"{lv}_N_I"] > 0)),
             ("all primary V_U", pd.Series(True, index=J.index))]
    for rname, rm in rows:
        for u in UPTAKE + ["total"]:
            m = rm & ((J.uptake_class == u) if u != "total" else True)
            sub = J[m]
            d = cell(sub)
            d.update(level=lv, row=rname, uptake_class=u,
                     n_with_band=int((sub[f"{lv}_N_I_band"] > 0).sum()), n_with_coverage=int((sub[f"{lv}_N_I_cov"] > 0).sum()),
                     sum_N_U=int(sub[f"{lv}_N_U"].sum()), sum_N_M=int(sub[f"{lv}_N_M"].sum()),
                     sum_N_I_band=int(sub[f"{lv}_N_I_band"].sum()), sum_N_I_cov=int(sub[f"{lv}_N_I_cov"].sum()))
            out.append(d)
    return pd.DataFrame(out)


def join(U, Ls):
    J = U[["case_id", "component_id", "volume_native_mL", "uptake_class", "SUVmax_PSMA", "extent_mm"]].copy()
    key = ["case_id", "component_id"]
    for lv, L in Ls.items():
        L = L[L.primary_status == "V_U"].set_index(key)
        N = level_counts(L).add_prefix(f"{lv}_")
        if len(N) != len(J) or set(N.index) != set(map(tuple, J[key].to_numpy())):
            raise RuntimeError(f"{lv}: component sets differ from the uptake table")
        J = J.join(N, on=key)
        if abs(L.volume_native_mL.sum() - J.volume_native_mL.sum()) > 1e-9:
            raise RuntimeError(f"{lv}: native volume differs")
        J[f"{lv}_pattern"] = [pattern(a, b, c) for a, b, c in zip(J[f"{lv}_N_U"], J[f"{lv}_N_M"], J[f"{lv}_N_I"])]
    return J


# ------------------------------ self-test ------------------------------
def selftest():
    res = []

    def ok(name, cond):
        res.append((name, bool(cond))); print(("PASS " if cond else "FAIL ") + name)

    ok("pattern: 48 unmatched -> stable", pattern(48, 2, 0) == PATTERN[0])
    ok("pattern: 47/3/0 -> matched only", pattern(47, 3, 0) == PATTERN[1])
    ok("pattern: 47/0/3 -> indeterminate only", pattern(47, 0, 3) == PATTERN[2])
    ok("pattern: 40/5/5 -> both", pattern(40, 5, 5) == PATTERN[3])
    ok("pattern: stable with N_M and N_I not reclassified", pattern(48, 1, 1) == PATTERN[0])
    mk = lambda u, m, b, c, r=0, s=0: dict(n_V_U=u, n_V_M=m, n_V_I_band=b, n_V_I_coverage=c, n_V_I_registration=r, n_V_S=s)
    L = pd.DataFrame([mk(50, 0, 0, 0), mk(40, 6, 3, 1)])
    N = level_counts(L)
    ok("N_I = band + coverage", list(N.N_I) == [0, 4])
    for bad, nm in ((mk(49, 0, 0, 0), "sum != 50 stops"), (mk(49, 0, 0, 0, r=1), "registration V_I > 0 stops"),
                    (mk(49, 0, 0, 0, s=1), "V_S > 0 stops")):
        try:
            level_counts(pd.DataFrame([bad])); ok(nm, False)
        except RuntimeError:
            ok(nm, True)
    U = pd.DataFrame(dict(case_id=["a", "a", "b", "c"], component_id=[1, 2, 1, 1], volume_native_mL=[1.0, 2.0, 0.5, 4.0],
                          uptake_class=["absent", "low", "NC_partial_volume", "absent"], SUVmax_PSMA=[1, 2, 3, 1],
                          extent_mm=[20, 20, 5, 30]))
    base = dict(case_id=["a", "a", "b", "c"], component_id=[1, 2, 1, 1], volume_native_mL=[1.0, 2.0, 0.5, 4.0],
                primary_status=["V_U"] * 4, n_V_I_registration=[0] * 4, n_V_S=[0] * 4)
    L1 = pd.DataFrame(dict(base, n_V_U=[50, 30, 47, 10], n_V_M=[0, 20, 0, 30], n_V_I_band=[0, 0, 3, 5], n_V_I_coverage=[0, 0, 0, 5]))
    extra = pd.DataFrame(dict(case_id=["z"], component_id=[9], volume_native_mL=[7.0], primary_status=["V_M"], n_V_I_registration=[0],
                              n_V_S=[0], n_V_U=[1], n_V_M=[49], n_V_I_band=[0], n_V_I_coverage=[0]))
    J = join(U, {"L1": pd.concat([L1, extra], ignore_index=True)})
    ok("join keeps primary V_U only and all components", len(J) == 4)
    X = crosstab(J, "L1")
    g = lambda r, u: X[(X.row == r) & (X.uptake_class == u)].iloc[0]
    ok("stable x absent: 1 comp, 1.0 mL, 1 patient", (g(PATTERN[0], "absent").components, g(PATTERN[0], "absent").native_mL) == (1, 1.0))
    ok("matched only x low: 1 comp, 2.0 mL", g(PATTERN[1], "low").components == 1 and g(PATTERN[1], "low").native_mL == 2.0)
    ok("indeterminate only x partial volume, band cause counted",
       g(PATTERN[2], "NC_partial_volume").components == 1 and g(PATTERN[2], "NC_partial_volume").n_with_band == 1
       and g(PATTERN[2], "NC_partial_volume").n_with_coverage == 0)
    ok("both x absent: band and coverage both counted", g(PATTERN[3], "absent").n_with_band == 1 and g(PATTERN[3], "absent").n_with_coverage == 1)
    ok("flag totals non-exclusive: matched any = 2, indeterminate any = 2",
       g(FLAGS[0], "total").components == 2 and g(FLAGS[1], "total").components == 2)
    ok("disjoint patterns sum to all", sum(g(p, "total").components for p in PATTERN) == g("all primary V_U", "total").components == 4)
    ok("volumes reconcile", abs(sum(g(p, "total").native_mL for p in PATTERN) - 7.5) < 1e-12)
    ok("patients counted per cell (absent total = 2 patients)", g("all primary V_U", "absent").patients == 2)
    ok("perturbation counts reconcile: sum N_U+N_M+N_I = 50 x components",
       g("all primary V_U", "total")[["sum_N_U", "sum_N_M", "sum_N_I_band", "sum_N_I_cov"]].sum() == 200)
    try:
        join(U.iloc[:3], {"L1": L1}); ok("component-set mismatch stops", False)
    except RuntimeError:
        ok("component-set mismatch stops", True)
    n_ok = sum(c for _, c in res)
    print(f"\nSELF-TEST {n_ok}/{len(res)} passed")
    return n_ok == len(res)


# ------------------------------ run ------------------------------
def fmt_table(X, lv):
    L = [f"Level {lv}: components / native mL / patients", "  " + "row".ljust(46) + "".join(u[:17].rjust(19) for u in UPTAKE + ["total"])]
    for r in PATTERN + FLAGS + ["all primary V_U"]:
        s = "  " + r.ljust(46)
        for u in UPTAKE + ["total"]:
            c = X[(X.row == r) & (X.uptake_class == u)].iloc[0]
            s += f"{c.components:5d} /{c.native_mL:6.1f} /{c.patients:3d}".rjust(19)
        L.append(s)
    L.append("  indeterminate causes (components with N_I,band > 0 / N_I,coverage > 0):")
    for r in PATTERN[2:]:
        c = X[(X.row == r) & (X.uptake_class == "total")].iloc[0]
        L.append(f"    {r}: band {c.n_with_band}, coverage {c.n_with_coverage}")
    return L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not selftest():
        sys.exit("self-test failed: no run")
    bad = [f"{d / n}: {md5(d / n)} != {m}" for (d, n), m in INPUTS.items() if md5(d / n) != m]
    if bad:
        sys.exit("INPUT MD5 MISMATCH:\n" + "\n".join(bad))
    U = pd.read_csv(OUT_DIR / PU_RUN / "PU_components_20261112_115324.csv")
    if len(U) != 217:
        sys.exit("expected 217 components in the uptake table")
    Ls = {lv: pd.read_csv(OUT_DIR / d / f"P3P_{lv}_components_{d[-15:]}.csv") for lv, d in LEVELS.items()}
    J = join(U, Ls)
    X = pd.concat([crosstab(J, lv) for lv in LEVELS], ignore_index=True)
    for lv in LEVELS:
        t = X[(X.level == lv) & (X.row == "all primary V_U") & (X.uptake_class == "total")].iloc[0]
        if t.components != 217 or abs(t.native_mL - U.volume_native_mL.sum()) > 1e-9:
            sys.exit(f"{lv}: reconciliation failed")
    run = datetime.now().strftime("%Y%m%d_%H%M%S"); odir = OUT_DIR / f"PX_{run}"; odir.mkdir()
    S = ["Exploratory cross-tabulation (POST-HOC, documented before results): perturbation stability x PSMA uptake class",
         "Units: 217 primary V_U components. Levels L1 (lower magnitude) and L2 (higher magnitude) never combined.",
         "Per level N_U + N_M + N_I = 50 (checked); stably unmatched = N_U >= 48; unstable flags non-exclusive.", ""]
    for lv in LEVELS:
        S += fmt_table(X[X.level == lv], lv) + [""]
    S += ["Interpretation framework (crosstab note section 4; imaging evidence, not validated biological diagnoses):"]
    S += [f"  {k[0]} x {k[1]}: {v}" for k, v in INTERP.items() if isinstance(k, tuple)]
    S += [f"  unstable, sometimes matched x any: {INTERP['matched']}", f"  unstable, sometimes indeterminate x any: {INTERP['indeterminate']}",
          f"  any x not classifiable: {INTERP['NC']}",
          "Stability != accuracy: stable under this perturbation model only (rotation about the patient FDG centroid; 50 sampled "
          "directions); perturbation magnitudes are empirically motivated, chest-wall-derived translations and pelvic-correction "
          "rotations, not anatomically calibrated lesion-level errors."]
    out = {}
    for nm, df in (("components", J), ("crosstab", X)):
        p = odir / f"PX_{nm}_{run}.csv"; df.to_csv(p, index=False); out[p.name] = md5(p)
    p = odir / f"PX_summary_{run}.txt"; p.write_text("\n".join(S) + "\n", encoding="utf-8"); out[p.name] = md5(p)
    man = dict(stage="PX exploratory cross-tabulation (post-hoc, documented before results)", run=run,
               script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs={f"{d / n}": m for (d, n), m in INPUTS.items()},
               params=dict(NPERT=NPERT, STABLE_MIN=STABLE_MIN), software=dict(python=platform.python_version(),
               numpy=np.__version__, pandas=pd.__version__), outputs=out, status="COMPLETE")
    (odir / f"PX_manifest_{run}.json").write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print("\n".join(S)); print(f"\nOutput: {odir}")


if __name__ == "__main__":
    main()
