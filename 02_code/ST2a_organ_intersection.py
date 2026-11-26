"""
ST2a_organ_intersection.py - PAPER 3, stage ST2a: intersection of eligible native FDG components with the native FDG
CT organ masks for the ST2 visceral V_I stress test. NATIVE FDG GEOMETRY ONLY: no registration, no PSMA data, no
matching, no DMTVF. Classification is done in MATLAB (P3_matching_classification_v1_9.m).

Implements: soft-tissue specification LOCKED (MD5 055199c4fd5d46039dcc65b6ab4f2416) SP1, SP6; dated implementation
clarification 00_protocol\\StepSAP_ST2_implementation_clarification_20261112.md (MD5 8cf2c01c306007ec00f2450f7e2bd993),
author decisions D1-D5 and the registration-first safeguard (applied in P3).

Per patient:
  availability  FDG rows of Audit 1B, per TotalSegmentator label (spleen 1, kidney R 2, kidney L 3, liver 5, lung lobes
                10, 11 (left), 12, 13, 14 (right)): available if category does not start with ABSENT (present = 1;
                disagreement stops) AND frac_largest_component >= 0.95. Lungs: union of the qualifying lobes only;
                included / excluded lobes and reasons recorded; a lung missing any lobe is flagged partial.
  intersection  eligible components (26-connectivity, rule A; P1 label image on the native FDG PET grid): physical
                centre of every native voxel -> continuous index on the FDG CT grid (totseg_24.nii.gz) -> nearest CT
                voxel by floor(x + 0.5); outside the CT grid -> no label. Counts per structure (available labels only).
                ST2 flag = any count > 0. Descriptive structure = largest count; tie -> lowest structure label ID
                (spleen 1, kidney_right 2, kidney_left 3, liver 5, lung_left 10, lung_right 12).
Outputs (06_output\\ST2a_[PREFLIGHT_]<run>\\): ST2a_components_<run>.csv, ST2a_availability_<run>.csv, manifest.
Checks (stop on failure): C1 MD5 of all inputs; P1 label image MD5 = P1 manifest record and voxel counts = P1 table;
  every O1 displacement-eligible structure ST2-available with all its labels; Audit 1B present/category agreement.
Usage (folder 02_code):
  & "<totalseg python>" -u ST2a_organ_intersection.py --selftest | --preflight | (full run)
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
P1_RUN = "P1_20261111_233936"
CONN = 26
INPUTS = {
    (OUT_DIR, f"{P1_RUN}/P1_components_20261111_233936.csv"): "ffafbc7917e1b595eadacafab247fe88",
    (OUT_DIR, f"{P1_RUN}/P1_manifest_20261111_233936.json"): "a27e4e4c5f44402c53487f9f8dd7dd8b",
    (OUT_DIR, "Audit1B_label_level_20261029_134401.csv"): "896a1681e90da4fe793a3a449d737fde",
    (OUT_DIR, "O1_organ_displacement_20261111_232613.csv"): "daaae3134932b6267807a128cb8f85dd",
    (OUT_DIR, "A3_visualQC_KEY_do_not_open_20261105.csv"): "c9a4548d3daf13d3e3dadc96171a9a94",   # preflight IDs only
    (PROTO_DIR, "StepSAP_soft_tissue_ST_spec_LOCKED_20261111.md"): "055199c4fd5d46039dcc65b6ab4f2416",
    (PROTO_DIR, "StepSAP_ST2_implementation_clarification_20261112.md"): "8cf2c01c306007ec00f2450f7e2bd993",
}
STRUCT = {"spleen": [1], "kidney_right": [2], "kidney_left": [3], "liver": [5],
          "lung_left": [10, 11], "lung_right": [12, 13, 14]}          # same structures and labels as O1
STRUCT_ID = {s: min(l) for s, l in STRUCT.items()}                   # tie rule: lowest label ID
LABEL_TO_STRUCT = {lab: s for s, ls in STRUCT.items() for lab in ls}
FRAG_MIN = 0.95


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------ pure helpers (tested) ------------------------------
def label_available(row):
    """row: Audit 1B FDG row (dict-like) or None. Returns (available, reason)."""
    if row is None:
        return False, "missing in Audit 1B"
    cat, pres, fr = str(row["category"]), int(row["present"]), float(row["frac_largest_component"])
    absent = cat.upper().startswith("ABSENT")
    if absent == bool(pres):
        raise RuntimeError(f"Audit 1B present={pres} disagrees with category {cat}")
    if absent:
        return False, f"absent ({cat})"
    if not (fr >= FRAG_MIN):
        return False, f"fragmented (frac_largest {fr:.3f})"
    return True, cat


def availability(cid, L1B):
    """per label: (available, reason); per structure: included / excluded labels, partial flag."""
    lab = {}
    for l in LABEL_TO_STRUCT:
        try:
            r = L1B.loc[(cid, "FDG", l)]
        except KeyError:
            r = None
        lab[l] = label_available(r)
    st = {}
    for s, ls in STRUCT.items():
        inc = [l for l in ls if lab[l][0]]
        exc = [f"{l}: {lab[l][1]}" for l in ls if not lab[l][0]]
        st[s] = dict(available=len(inc) > 0, labels_included="+".join(map(str, inc)), labels_excluded="; ".join(exc),
                     partial=(len(inc) > 0 and len(exc) > 0), complete=(len(exc) == 0),
                     truncated_included=[l for l in inc if lab[l][1] != "PRESENT_COMPLETE"])
    return lab, st


def ct_index(pts, origin, spacing, direction, size):
    """nearest CT voxel (k, j, i) of physical points (LPS, mm) by floor(x + 0.5); rows outside the grid -> -1."""
    D = np.asarray(direction, float).reshape(3, 3)
    ci = np.linalg.solve(D * np.asarray(spacing, float), (np.atleast_2d(pts) - np.asarray(origin, float)).T).T
    ijk = np.floor(ci + 0.5).astype(np.int64)
    inside = np.all((ijk >= 0) & (ijk < np.asarray(size, np.int64)), axis=1)
    kji = ijk[:, ::-1].copy(); kji[~inside] = -1
    return kji, inside


def labels_at(ct_arr, kji, inside):
    out = np.zeros(len(kji), np.int64)
    if inside.any():
        k = kji[inside]
        out[inside] = ct_arr[k[:, 0], k[:, 1], k[:, 2]]
    return out


def structure_counts(point_labels, avail_labels):
    """voxel-centre counts per structure, available labels only."""
    c = {s: 0 for s in STRUCT}
    for l in avail_labels:
        n = int(np.sum(point_labels == l))
        if n:
            c[LABEL_TO_STRUCT[l]] += n
    return c


def assign_structure(counts):
    """largest count; tie -> lowest structure label ID; '' if no count."""
    best = max(counts.values()) if counts else 0
    if best == 0:
        return ""
    return min((s for s, n in counts.items() if n == best), key=lambda s: STRUCT_ID[s])


def phys_lps(img, kji):            # verbatim P1
    D = np.array(img.GetDirection()).reshape(3, 3); S = np.array(img.GetSpacing())
    return np.array(img.GetOrigin()) + (kji[:, ::-1].astype(float) * S) @ D.T


def preflight_cases(cases, oldkey):
    return [c for c in cases if c not in set(oldkey["case_id"])][:2]


# ------------------------------ per patient ------------------------------
def process(cid, C1, L1B, O1, P1M):
    import SimpleITK as sitk
    lab_path = OUT_DIR / P1_RUN / f"P1_{cid}_FDG_conn{CONN}.nii.gz"
    if md5(lab_path) != P1M["label_images"].get(lab_path.name):
        raise RuntimeError(f"{cid}: P1 label image differs from the P1 manifest")
    lab_img = sitk.ReadImage(str(lab_path)); lab = sitk.GetArrayFromImage(lab_img)
    ct_img = sitk.ReadImage(str(DATA_ROOT / cid / "FDG" / "totseg_24.nii.gz")); ct = sitk.GetArrayFromImage(ct_img)
    labs, st = availability(cid, L1B)
    avail = [l for l, (a, _) in labs.items() if a]
    # O1 cross-check: displacement-eligible -> ST2-available with all labels
    for r in O1[O1.case_id == cid].itertuples():
        if str(r.eligible).lower() == "true" and not st[r.structure]["complete"]:
            raise RuntimeError(f"{cid}: {r.structure} displacement-eligible in O1 but not fully ST2-available")
    idx = ndimage.value_indices(lab, ignore_value=0)
    comps = C1[(C1.case_id == cid) & (C1.connectivity == CONN)]
    if len(comps) != len(idx) or any(len(idx[q][0]) != n for q, n in zip(comps.component_id, comps.voxels)):
        raise RuntimeError(f"{cid}: native label image does not match the P1 component table")
    g = dict(origin=ct_img.GetOrigin(), spacing=ct_img.GetSpacing(), direction=ct_img.GetDirection(), size=ct_img.GetSize())
    rows = []
    for r in comps.itertuples():
        row = dict(case_id=cid, component_id=r.component_id, voxels_native=r.voxels, volume_native_mL=r.volume_mL,
                   eligible_A=r.eligible_A)
        if r.eligible_A:
            pts = phys_lps(lab_img, np.column_stack(idx[r.component_id]))
            kji, inside = ct_index(pts, **g)
            pl = labels_at(ct, kji, inside)
            cnt = structure_counts(pl, avail)
            row.update(n_outside_ct=int((~inside).sum()), **{f"n_{s}": cnt[s] for s in STRUCT},
                       st2_flag=bool(sum(cnt.values()) > 0), st2_structure=assign_structure(cnt),
                       n_structures_touched=int(sum(v > 0 for v in cnt.values())))
        rows.append(row)
    av = [dict(case_id=cid, structure=s, structure_label_id=STRUCT_ID[s], available_ST2=v["available"],
               complete=v["complete"], partial=v["partial"], labels_included=v["labels_included"],
               labels_excluded=v["labels_excluded"], truncated_labels_included="+".join(map(str, v["truncated_included"])),
               displacement_eligible_O1=bool(O1[(O1.case_id == cid) & (O1.structure == s)].eligible.astype(str).str.lower()
                                            .eq("true").any()))
          for s, v in st.items()]
    return rows, av


# ------------------------------ self-test ------------------------------
def selftest():
    res = []
    ok = lambda n, c: res.append((n, bool(c)))
    R = lambda cat, pres, fr: dict(category=cat, present=pres, frac_largest_component=fr)
    ok("complete, unfragmented -> available", label_available(R("PRESENT_COMPLETE", 1, 1.0))[0])
    ok("truncated (boundary) -> available", label_available(R("PRESENT_BOUNDARY_INPLANE", 1, 0.99))[0])
    ok("absent -> not available", not label_available(R("ABSENT_UNEXPLAINED", 0, np.nan))[0])
    ok("frac 0.949 -> not available; 0.95 -> available", not label_available(R("PRESENT_COMPLETE", 1, 0.949))[0]
       and label_available(R("PRESENT_COMPLETE", 1, 0.95))[0])
    try:
        label_available(R("PRESENT_COMPLETE", 0, 1.0)); ok("present/category disagreement stops", False)
    except RuntimeError:
        ok("present/category disagreement stops", True)
    rows = [dict(case_id="c", tracer="FDG", label_id=l, category="PRESENT_COMPLETE", present=1,
                 frac_largest_component=1.0) for l in LABEL_TO_STRUCT]
    L = pd.DataFrame(rows).set_index(["case_id", "tracer", "label_id"])
    L.loc[("c", "FDG", 13), "frac_largest_component"] = 0.90
    L.loc[("c", "FDG", 1), ["category", "present", "frac_largest_component"]] = ["ABSENT_UNEXPLAINED", 0, np.nan]
    labs, st = availability("c", L)
    ok("lung: fragmented lobe excluded, union of qualifying lobes, flagged partial",
       st["lung_right"]["labels_included"] == "12+14" and st["lung_right"]["partial"] and not st["lung_right"]["complete"]
       and "13: fragmented" in st["lung_right"]["labels_excluded"])
    ok("complete lung not partial", st["lung_left"]["complete"] and not st["lung_left"]["partial"])
    ok("absent spleen -> structure not available", not st["spleen"]["available"])
    # centre-in-mask rule on synthetic grids: CT 1 mm, origin 0; PET centres
    G = dict(origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0), direction=(1, 0, 0, 0, 1, 0, 0, 0, 1), size=(10, 10, 10))
    ct = np.zeros((10, 10, 10), np.int64); ct[2:5, 2:5, 2:5] = 5                # liver voxels (k, j, i) 2..4 -> x,y,z 2..4
    P = np.array([[4.4, 3.0, 3.0], [4.6, 3.0, 3.0], [4.5, 3.0, 3.0], [-0.6, 3, 3], [9.4, 3, 3], [9.6, 3, 3]])
    kji, ins = ct_index(P, **G); pl = labels_at(ct, kji, ins)
    ok("centre 4.4 mm -> mask voxel 4 (in); 4.6 -> voxel 5 (out)", pl[0] == 5 and pl[1] == 0)
    ok("tie x.5 -> floor(x + 0.5) = upper voxel (4.5 -> 5, out)", pl[2] == 0)
    ok("outside CT grid -> no label; last voxel edge 9.4 in, 9.6 out", not ins[3] and ins[4] and not ins[5] and pl[3] == 0)
    Gf = dict(G, direction=(1, 0, 0, 0, 1, 0, 0, 0, -1), origin=(0.0, 0.0, 9.0))      # flipped z: z = 9 - k
    kf, inf_ = ct_index(np.array([[3.0, 3.0, 9.0 - 3.0]]), **Gf)
    ok("flipped CT axis handled", inf_[0] and labels_at(ct, kf, inf_)[0] == 5)
    ok("partial overlap without a centre in the mask is not counted (no supersampling)",
       structure_counts(labels_at(ct, *ct_index(np.array([[4.6, 3.0, 3.0]]), **G)), [5])["liver"] == 0)
    c = structure_counts(np.array([5, 5, 12, 12, 13, 0, 1]), [5, 12, 13])
    ok("counts use available labels only; lobes add to their lung", c == dict(spleen=0, kidney_right=0, kidney_left=0,
       liver=2, lung_left=0, lung_right=3))
    ok("largest count wins", assign_structure(c) == "lung_right")
    ok("tie -> lowest structure label ID (liver 5 < lung_right 12)",
       assign_structure(dict(c, lung_right=2)) == "liver")
    ok("no count -> unassigned", assign_structure({s: 0 for s in STRUCT}) == "")
    for n, g_ in res:
        print(f"  [{'PASS' if g_ else 'FAIL'}] {n}")
    print(f"selftest: {sum(g_ for _, g_ in res)}/{len(res)} passed")
    return all(g_ for _, g_ in res)


# ------------------------------ run ------------------------------
def run(preflight):
    import SimpleITK as sitk
    t0 = time.time(); run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    md = {n: md5(d / n) if (d / n).exists() else "MISSING" for (d, n) in INPUTS}
    bad = {n: md[n] for (d, n), e in INPUTS.items() if md[n] != e}
    if bad:
        sys.exit(f"STOP C1: MD5 mismatch or missing: {bad}")
    print(f"C1 MD5: {len(INPUTS)}/{len(INPUTS)} match")
    C1 = pd.read_csv(OUT_DIR / P1_RUN / "P1_components_20261111_233936.csv")
    P1M = json.loads((OUT_DIR / P1_RUN / "P1_manifest_20261111_233936.json").read_text(encoding="utf-8"))
    L1B = pd.read_csv(OUT_DIR / "Audit1B_label_level_20261029_134401.csv").set_index(["case_id", "tracer", "label_id"])
    O1 = pd.read_csv(OUT_DIR / "O1_organ_displacement_20261111_232613.csv")
    cases = sorted(C1.case_id.unique())
    if len(cases) != 100:
        sys.exit("STOP: expected 100 patients")
    if preflight:
        cases = preflight_cases(cases, pd.read_csv(OUT_DIR / "A3_visualQC_KEY_do_not_open_20261105.csv"))
    odir = OUT_DIR / f"ST2a_{'PREFLIGHT_' if preflight else ''}{run_id}"; odir.mkdir()
    rows, avs = [], []
    for i, cid in enumerate(cases, 1):
        r_, a_ = process(cid, C1, L1B, O1, P1M)
        rows += r_; avs += a_
        el = [r for r in r_ if r["eligible_A"]]
        print(f"  [{i:3d}/{len(cases)}] {cid}: eligible {len(el)}, ST2-flagged {sum(r['st2_flag'] for r in el)}  "
              f"({time.time() - t0:.0f} s)")
    R = pd.DataFrame(rows); A = pd.DataFrame(avs)
    E = R[R.eligible_A.astype(bool)]
    rp = odir / f"ST2a_components_{run_id}.csv"; R.to_csv(rp, index=False)
    ap_ = odir / f"ST2a_availability_{run_id}.csv"; A.to_csv(ap_, index=False)
    print(f"Eligible components {len(E)}; ST2-flagged {int(E.st2_flag.sum())} "
          f"(touching > 1 structure {int((E.n_structures_touched > 1).sum())}); voxel centres outside the CT grid "
          f"{int(E.n_outside_ct.sum())} in {int((E.n_outside_ct > 0).sum())} components")
    print("Availability (patients): " + "; ".join(
        f"{s} {int(A[A.structure == s].available_ST2.sum())} (partial {int(A[A.structure == s].partial.sum())})"
        for s in STRUCT))
    print("O1 cross-check (displacement-eligible -> fully ST2-available): PASS")
    man = dict(stage="ST2a organ intersection (native FDG geometry only)", mode="PREFLIGHT" if preflight else "FULL",
               run=run_id, script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs=md, n_patients=len(cases),
               connectivity=CONN, outputs={p.name: md5(p) for p in sorted(odir.glob("*")) if p.is_file()},
               params=dict(FRAG_MIN=FRAG_MIN, STRUCT=STRUCT, STRUCT_ID=STRUCT_ID, rounding="floor(x + 0.5)"),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=__import__("scipy").__version__, SimpleITK=sitk.Version_VersionString()),
               runtime_s=round(time.time() - t0), status="COMPLETE")
    mp = odir / f"ST2a_manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print(f"done: {odir.name}; manifest MD5 {md5(mp)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.preflight)
