"""
P2b_export_voxelsets_csv.py - PAPER 3: export the P2 per-component voxel sets and the PSMA PET grid geometry as CSV
for MATLAB (protocol v2.10 section 14.5: Python <-> MATLAB exchange only via NIfTI images and structured CSV tables;
P2 stored the voxel sets as .npz). Format conversion only: no value is computed, changed or filtered.

Inputs : 06_output\\P2_conn<c>_<run>\\P2_voxelsets_<case>.npz, P2_grids_<run>.csv (accepted P2 run, MD5-checked)
Outputs: 06_output\\P2b_conn<c>_<run>\\P2b_voxelsets_<case>.csv  (component_id, role[0 final, 1 alt], linear_index_1based)
         P2b_psma_pet_grid_<run>.csv (per case: size, spacing, origin and direction in SimpleITK LPS convention,
         fingerprint recomputed and checked against P2), manifest with MD5s.
linear_index_1based = P2 linear index (C order over (k, j, i)) + 1 = MATLAB column-major index of the (i, j, k) array
returned by niftiread for the same file.
Checks (stop on failure): input MD5s; round trip (CSV re-read == npz, every case); fingerprint recomputed == P2.
Usage: & "<totalseg python>" -u P2b_export_voxelsets_csv.py [--selftest]

v1.1 (12 Nov 2026, 00:50; file P2b_export_voxelsets_csv_v1_1.py): accepts any accepted P2 run (--p2dir, e.g. the
connectivity-sensitivity runs P2_conn18_<run>, P2_conn6_<run>). The P2 run is verified through its own manifest:
script = P2_mapping_coverage_resampling_v1_1.py (MD5 5022bbbd0e30945bdfcd2533f59e30d4), mode FULL, status COMPLETE,
and the MD5 of every .npz and of the grid table equal to the manifest record. Conversion unchanged.
  & "<totalseg python>" -u P2b_export_voxelsets_csv_v1_1.py --p2dir P2_conn18_<run>

v1.2 (12 Nov 2026, 01:05; file P2b_export_voxelsets_csv_v1_2.py): also accepts runs of P2 v1.2 (scenario runs S1/S2/S4;
P2_mapping_coverage_resampling_v1_2.py, MD5 d066195f1e5ffa3f1ece5a06333dff1a); output folder P2b_<scenario>_conn<c>_<run>. Conversion unchanged.

v1.3 (12 Nov 2026, 01:25; file P2b_export_voxelsets_csv_v1_3.py): also accepts runs of P2 v1.3 (scenarios S3b / S3c;
P2_mapping_coverage_resampling_v1_3.py, MD5 7c3e1b45b2f81e8d7e20289455ca32ca). Conversion unchanged.
"""
import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
P2_SCRIPT_MD5 = "5022bbbd0e30945bdfcd2533f59e30d4"
P2_ACCEPTED = {"P2_mapping_coverage_resampling_v1_1.py": "5022bbbd0e30945bdfcd2533f59e30d4",
               "P2_mapping_coverage_resampling_v1_2.py": "d066195f1e5ffa3f1ece5a06333dff1a",
               "P2_mapping_coverage_resampling_v1_3.py": "7c3e1b45b2f81e8d7e20289455ca32ca"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def to_frame(z):
    return pd.DataFrame(dict(component_id=z["component_id"].astype(np.int64), role=z["role"].astype(np.int64),
                             linear_index_1based=z["linear_index"].astype(np.int64) + 1))


def selftest():
    res = []
    ok = lambda n, c: res.append((n, bool(c)))
    shape = (4, 5, 6)                                   # (k, j, i) as in P2 (sitk array order)
    k, j, i = 2, 3, 4
    lin0 = np.ravel_multi_index((k, j, i), shape)
    nx, ny = shape[2], shape[1]
    ok("C-order (k,j,i) index + 1 == MATLAB column-major index of (i,j,k)", lin0 + 1 == i + nx * j + nx * ny * k + 1)
    z = dict(component_id=np.array([3, 3, 7]), role=np.array([0, 1, 0]), linear_index=np.array([0, 10, 119]))
    f = to_frame(z)
    ok("frame keeps ids, roles; index + 1", list(f.linear_index_1based) == [1, 11, 120] and list(f.role) == [0, 1, 0])
    for n, g in res:
        print(f"  [{'PASS' if g else 'FAIL'}] {n}")
    print(f"selftest: {sum(g for _, g in res)}/{len(res)} passed")
    return all(g for _, g in res)


def run(p2dir):
    import SimpleITK as sitk
    _s = importlib.util.spec_from_file_location("p2", HERE / "P2_mapping_coverage_resampling_v1_1.py")
    p2 = importlib.util.module_from_spec(_s); _s.loader.exec_module(p2)
    P2D = OUT_DIR / p2dir; P2_RUN = p2dir[-15:]
    mp2 = P2D / f"P2_manifest_{P2_RUN}.json"
    m2 = json.loads(mp2.read_text(encoding="utf-8"))
    sname = m2.get("script")
    if sname not in P2_ACCEPTED or m2.get("script_md5") != P2_ACCEPTED[sname] or md5(HERE / sname) != P2_ACCEPTED[sname]:
        sys.exit(f"STOP C1: {p2dir} was not produced by an accepted P2 script (or the script file changed)")
    if m2.get("mode") != "FULL" or m2.get("status") != "COMPLETE":
        sys.exit(f"STOP: {p2dir} is not a complete FULL run")
    scen = m2.get("scenario", "none")
    conn = int(m2["connectivity"])
    for name in (f"P2_grids_{P2_RUN}.csv", f"P2_components_{P2_RUN}.csv"):
        if md5(P2D / name) != m2["outputs"].get(name):
            sys.exit(f"STOP C1: {name} differs from its P2 manifest record")
    md = {mp2.name: md5(mp2), sname: P2_ACCEPTED[sname]}
    print(f"C1: P2 run {p2dir} (connectivity {conn}) verified via its manifest (MD5 {md5(mp2)})")
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    odir = OUT_DIR / f"P2b_{'' if scen == 'none' else scen + '_'}conn{conn}_{run_id}"; odir.mkdir()
    G = pd.read_csv(P2D / f"P2_grids_{P2_RUN}.csv")
    rows, outs = [], {}
    for cid in G.case_id:
        npz = P2D / f"P2_voxelsets_{cid}.npz"
        if md5(npz) != m2["outputs"].get(npz.name):
            sys.exit(f"STOP C1: {npz.name} differs from its P2 manifest record")
        z = np.load(npz)
        f = to_frame(z)
        fp = odir / f"P2b_voxelsets_{cid}.csv"; f.to_csv(fp, index=False)
        back = pd.read_csv(fp)
        if not (back.equals(f) and len(back) == len(z["linear_index"])):
            sys.exit(f"STOP: round trip failed for {cid}")
        g = p2.read_grid(DATA_ROOT / cid / "PSMA" / "PET.nii.gz")
        fpr = p2.fingerprint(**g)
        if fpr != G.set_index("case_id").loc[cid, "psma_pet_fingerprint"]:
            sys.exit(f"STOP: PSMA PET fingerprint differs from P2 for {cid}")
        shp = tuple(int(x) for x in z["grid_shape_kji"])
        if shp != tuple(int(x) for x in g["size"][::-1]):
            sys.exit(f"STOP: voxel-set grid shape differs from the PSMA PET header for {cid}")
        rows.append(dict(case_id=cid, size_x=g["size"][0], size_y=g["size"][1], size_z=g["size"][2],
                         spacing_x=g["spacing"][0], spacing_y=g["spacing"][1], spacing_z=g["spacing"][2],
                         origin_x_LPS=g["origin"][0], origin_y_LPS=g["origin"][1], origin_z_LPS=g["origin"][2],
                         **{f"dir_{a}{b}": g["direction"][3 * a + b] for a in range(3) for b in range(3)},
                         fingerprint=fpr, npz_md5=md5(npz), csv_md5=md5(fp), n_rows=len(f)))
        outs[fp.name] = md5(fp)
    gp = odir / f"P2b_psma_pet_grid_{run_id}.csv"; pd.DataFrame(rows).to_csv(gp, index=False)
    outs[gp.name] = md5(gp)
    man = dict(stage="P2b export of P2 voxel sets and PSMA PET grids to CSV (format conversion only)", run=run_id,
               source_P2_run=p2dir, connectivity=conn, scenario=scen, script=Path(__file__).name, script_md5=md5(Path(__file__)), inputs=md,
               outputs=outs, n_cases=len(rows), round_trip="identical for every case",
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             SimpleITK=sitk.Version_VersionString()), status="COMPLETE")
    mp = odir / f"P2b_manifest_{run_id}.json"; mp.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"done: {odir.name}; {len(rows)} cases; round trip and fingerprints OK; manifest MD5 {md5(mp)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--p2dir", default="P2_conn26_20261111_235603"); a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    run(a.p2dir)
