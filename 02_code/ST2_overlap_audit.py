"""
ST2_overlap_audit.py - PAPER 3: descriptive organ-overlap audit of the ST2-withheld components (author approval,
12 Nov 2026, 07:31; ST2 completion note addendum MD5 d16a160add3f2fc9f5180857c486586b, item 4). Descriptive only:
no classification, no result of the locked SAP is changed. Uses the existing ST2a and P3S_ST2 tables (no imaging).

f_organ = (native FDG voxel centres inside the UNION of the available organ masks) / (native voxels of the component)
Union: the six structures come from ONE TotalSegmentator label map (mutually exclusive labels, SP1), and ST2a looks up
each voxel centre at exactly one CT voxel; each centre therefore carries at most one label and the per-structure counts
are disjoint, so their sum is the union count. This is verified (sum of counts <= native voxels - centres outside the
CT grid, for every component); otherwise the script stops rather than approximating.
Withheld components = P3S_ST2 cause 'ST2 organ' (checked identical to the ST2a flag, none previously V_I registration).
Outputs (06_output\\ST2_overlap_audit_<run>\\): per-component CSV, summary JSON, manifest (MD5s).
"""
import hashlib
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

OUT_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
INPUTS = {
    "ST2a_20261112_072150/ST2a_components_20261112_072150.csv": None,
    "ST2a_20261112_072150/ST2a_manifest_20261112_072150.json": "bca65313ea0d550286fbf9c38bcfc33b",
    "P3S_ST2_20261112_072438/P3S_ST2_components_20261112_072438.csv": "6511e142d7b88f64d405526800a867cb",
    "P3S_ST2_20261112_072438/P3S_manifest_20261112_072438.json": "641328dda1f30d59fd6ec985b20cf051",
    "P3S_ST2_completion_note_addendum_20261112.md": "d16a160add3f2fc9f5180857c486586b",
}
STRUCT = ["spleen", "kidney_right", "kidney_left", "liver", "lung_left", "lung_right"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def wquantile(x, w, p):
    """volume-weighted quantile: smallest x whose cumulative weight share >= p."""
    o = np.argsort(x, kind="stable"); x, w = np.asarray(x)[o], np.asarray(w)[o]
    c = np.cumsum(w) / w.sum()
    return float(x[np.searchsorted(c, p - 1e-12)])


def main():
    run = datetime.now().strftime("%Y%m%d_%H%M%S")
    m2a = json.loads((OUT_DIR / "ST2a_20261112_072150/ST2a_manifest_20261112_072150.json").read_text(encoding="utf-8"))
    INPUTS["ST2a_20261112_072150/ST2a_components_20261112_072150.csv"] = m2a["outputs"]["ST2a_components_20261112_072150.csv"]
    md = {k: md5(OUT_DIR / k) for k in INPUTS}
    bad = {k: v for k, v in md.items() if v != INPUTS[k]}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch {bad}")
    S = pd.read_csv(OUT_DIR / "ST2a_20261112_072150/ST2a_components_20261112_072150.csv")
    C = pd.read_csv(OUT_DIR / "P3S_ST2_20261112_072438/P3S_ST2_components_20261112_072438.csv")
    J = C.merge(S, on=["case_id", "component_id"], how="inner", validate="one_to_one", suffixes=("", "_st2a"))
    if len(J) != len(C) or len(C) != len(S):
        sys.exit("STOP: component tables do not align")
    flag = J.st2_flag.astype(str).str.lower().eq("true")
    W = J[J.cause == "ST2 organ"].copy()
    if not (flag == (J.cause == "ST2 organ")).all():
        sys.exit("STOP: ST2a flag and P3S ST2 status differ")
    cnt = W[[f"n_{s}" for s in STRUCT]].sum(axis=1)
    if not ((cnt <= W.voxels_native - W.n_outside_ct) & (cnt > 0)).all():
        sys.exit("STOP: per-structure counts cannot be the union count (union not reconstructable); not approximated")
    W["n_in_union"] = cnt.astype(int)
    W["f_organ"] = W.n_in_union / W.voxels_native
    W["volume_in_union_mL"] = W.volume_native_mL * W.f_organ
    v, f = W.volume_native_mL.values, W.f_organ.values
    summ = dict(
        n_components=int(len(W)), n_patients=int(W.case_id.nunique()), withheld_mL=float(v.sum()),
        withheld_inside_union_mL=float(W.volume_in_union_mL.sum()),
        f_organ_unweighted=dict(median=float(np.median(f)), q25=float(np.percentile(f, 25)), q75=float(np.percentile(f, 75)),
                                min=float(f.min()), max=float(f.max())),
        f_organ_volume_weighted=dict(median=wquantile(f, v, .5), q25=wquantile(f, v, .25), q75=wquantile(f, v, .75),
                                     mean=float(np.average(f, weights=v))),
        pct_withheld_volume_from_components_with_f_below={str(t): float(100 * v[f < t].sum() / v.sum()) for t in (.10, .25, .50)},
        n_components_with_f_below={str(t): int((f < t).sum()) for t in (.10, .25, .50)},
        by_structure={s: dict(n=int((W.st2_structure == s).sum()), withheld_mL=float(W.volume_native_mL[W.st2_structure == s].sum()),
                              inside_union_mL=float(W.volume_in_union_mL[W.st2_structure == s].sum()),
                              median_f=float(W.f_organ[W.st2_structure == s].median()) if (W.st2_structure == s).any() else None)
                      for s in STRUCT},
        prior_primary_status_note="see ST2 completion note: 145 matched (3,767.9 mL), 8 unmatched (3.1 mL)")
    od = OUT_DIR / f"ST2_overlap_audit_{run}"; od.mkdir()
    cp = od / f"ST2_overlap_components_{run}.csv"
    W[["case_id", "component_id", "voxels_native", "volume_native_mL", "n_outside_ct"] + [f"n_{s}" for s in STRUCT]
      + ["n_in_union", "f_organ", "volume_in_union_mL", "st2_structure", "n_structures_touched"]].to_csv(cp, index=False)
    sp = od / f"ST2_overlap_summary_{run}.json"; sp.write_text(json.dumps(summ, indent=2), encoding="utf-8")
    man = dict(stage="ST2 descriptive organ-overlap audit", run=run, script=Path(__file__).name, script_md5=md5(Path(__file__)),
               inputs=md, outputs={cp.name: md5(cp), sp.name: md5(sp)},
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__), status="COMPLETE")
    mp = od / f"ST2_overlap_manifest_{run}.json"; mp.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(json.dumps(summ, indent=1))
    print(f"done: {od.name}; manifest MD5 {md5(mp)}")


if __name__ == "__main__":
    main()
