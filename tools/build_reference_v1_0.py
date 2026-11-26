"""
build_reference_v1_0.py - PAPER 3 repository: build the reference checksum list for reproduce_runner.

Reads the run list (tools/reproduce_runs_v1_0.csv) and, for every run, the ORIGINAL output files in the original
06_output folder (read-only). Writes reference/reference_checksums_v1_0.csv with one row per original output file:
  run_name, relpath (relative to 06_output, forward slashes), bytes, md5, cls
  cls = gating   (.csv .npz .nii .nii.gz .mat: must be byte-identical on reproduction)
        manifest (*manifest*.json: replaced by the original after all gating files of the run match)
        info     (everything else: logs, images, spreadsheets, text reports; compared and reported, never gating)
Run-list kinds: folder = output folder named run_name; flat = top-level entries whose name contains run_name (a
time stamp); argdir / argfile = folder named run_name, minus the published input files listed in 'exclude'.
Usage: python build_reference_v1_0.py <original 06_output> <run list csv> <out csv> [first_order last_order]
"""
import csv
import hashlib
import sys
from pathlib import Path

GATING = (".csv", ".npz", ".nii", ".nii.gz", ".mat")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def cls_of(name):
    n = name.lower()
    if "manifest" in n and n.endswith(".json"):
        return "manifest"
    return "gating" if n.endswith(GATING) else "info"


def files_of(out, r):
    name, kind = r["run_name"], r["kind"]
    excl = set(x for x in r["exclude"].split(";") if x)
    if kind == "flat":
        tops = [p for p in out.iterdir() if name in p.name]
    else:
        tops = [out / name]
    fs = []
    for t in tops:
        if not t.exists():
            sys.exit(f"STOP: {t} not found")
        fs += [t] if t.is_file() else [p for p in t.rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    return sorted(p for p in fs if p.name not in excl)


def selftest():
    assert cls_of("P2_manifest_k01.json") == "manifest" and cls_of("X_manifest_1.JSON") == "manifest"
    assert cls_of("a.csv") == "gating" and cls_of("b.nii.gz") == "gating" and cls_of("c.npz") == "gating"
    assert cls_of("summary.json") == "info" and cls_of("fig.png") == "info" and cls_of("log.txt") == "info"
    print("selftest 7/7 pass")


def main(out, runs, dst, lo=None, hi=None):
    out = Path(out)
    R = list(csv.DictReader(open(runs, encoding="utf-8")))
    if lo is not None:
        R = [r for r in R if int(lo) <= int(r["order"]) <= int(hi)]
    new = not Path(dst).exists()
    with open(dst, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["run_name", "relpath", "bytes", "md5", "cls"])
        for r in R:
            fs = files_of(out, r)
            for p in fs:
                w.writerow([r["run_name"], p.relative_to(out).as_posix(), p.stat().st_size, md5(p), cls_of(p.name)])
            print(f"{r['order']:>3} {r['run_name']}: {len(fs)} files")


if __name__ == "__main__":
    selftest()
    main(*sys.argv[1:])
