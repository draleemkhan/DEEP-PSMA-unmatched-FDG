"""
reproduce_runner.py - PAPER 3: rerun the archived analysis pipeline and verify byte-for-byte reproduction.

Version 1.0, 25 November 2026; 1.1, 26 November 2026 (see below). The archived scripts in 02_code are run UNCHANGED, in the original order
(tools/reproduce_runs_v1_0.csv, 66 runs). After each run the new outputs are renamed to the original run names
(time stamps), so that later scripts find the file names and checksums they expect, and every output file is compared
with the original (reference/reference_checksums_v1_0.csv):
  gating   (.csv .npz .nii .nii.gz .mat)  must be byte-identical (MD5); any difference stops the runner.
  manifest (*manifest*.json)               run manifests record run times and therefore always differ. Only when
                                           every gating file of the run matches is the new manifest kept as
                                           <name>_rerun.json and the original manifest (reference/manifests) put in
                                           its place, so that later scripts can verify it.
  info     (logs, images, text reports)    compared and reported; never stops the runner.
Nothing is ever overwritten: a run whose output names already exist (and is not recorded as passed) is refused.
Version 1.1, 26 November 2026 - line endings: a few original runs (Step 5 QC analyses, ST2 overlap audit) were made
in a Linux environment and wrote text tables with LF line endings; pandas on Windows writes CRLF. A text result file
(.csv) that differs ONLY in line endings is converted to LF; it is accepted only if the converted file then has exactly
the original MD5 (status PASS_EOL; the file on disk is then byte-identical to the original). No other difference is
tolerated.

Requirements: Windows; repository at D:\\DEEP PSMA\\PAPER3_Discordance; DEEP-PSMA v1 extracted at
D:\\DEEP PSMA\\DEEP PSMA EXTRACTED (the archived scripts use these absolute paths; `subst D: <folder>` creates a
D: drive if needed); reproduce_config.json with the Python and MATLAB executables.

Usage (from the repository root):
  python reproduce_runner.py --selftest            synthetic tests of renaming and verification (no data)
  python reproduce_runner.py --dry-run             check paths, executables, versions and archived scripts
  python reproduce_runner.py --list                print the run list
  python reproduce_runner.py [--only 1-3] [--continue]
Outputs: reproduction/state.json (passed runs; used to resume), reproduction/logs/<order>_<run>.log,
         reproduction/reproduction_report_<time>.csv (one row per file) and reproduction_summary_<time>.csv.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = Path(r"D:\DEEP PSMA\PAPER3_Discordance")
DATA = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
RUNS = ROOT / "tools" / "reproduce_runs_v1_0.csv"
REF = ROOT / "reference" / "reference_checksums_v1_0.csv"
REF_MAN = ROOT / "reference" / "manifests"
TOKEN = re.compile(r"\d{8}_\d{6}")
GATING = (".csv", ".npz", ".nii", ".nii.gz", ".mat")


# ------------------------------------------------------------------ helpers (tested) ----------
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


TEXT = (".csv",)


def lf(b):
    return b.replace(b"\r\n", b"\n")


def fix_eol(out, rows):
    """Convert to LF the text result files that differ from the original only in line endings (verified by MD5)."""
    fixed = []
    for x in rows:
        p = out / x["relpath"]
        if x["cls"] != "gating" or not p.name.lower().endswith(TEXT) or not p.exists() or md5(p) == x["md5"]:
            continue
        b = lf(p.read_bytes())
        if hashlib.md5(b).hexdigest() == x["md5"]:
            p.write_bytes(b)
            fixed.append(x["relpath"])
    return fixed


def ref_token(rows):
    """The single original time stamp of a run (from its reference file names), or None."""
    toks = {t for r in rows for t in TOKEN.findall(r["relpath"])}
    if len(toks) > 1:
        raise RuntimeError(f"more than one time stamp in reference: {sorted(toks)}")
    return toks.pop() if toks else None


def rename_tree(top, new, old):
    """Rename `top` and everything below it, replacing time stamp `new` by `old` in every name (deepest first)."""
    if top.is_dir():
        for p in sorted(top.rglob("*"), key=lambda q: len(q.parts), reverse=True):
            if new in p.name:
                p.rename(p.with_name(p.name.replace(new, old)))
    if new in top.name:
        dst = top.with_name(top.name.replace(new, old))
        if dst.exists():
            raise RuntimeError(f"refusing to overwrite {dst}")
        top.rename(dst)
        return dst
    return top


def output_files(out, r, token):
    name, kind = r["run_name"], r["kind"]
    excl = {x for x in r["exclude"].split(";") if x}
    tops = [p for p in out.iterdir() if token and token in p.name] if kind == "flat" else [out / name]
    fs = []
    for t in tops:
        if t.exists():
            fs += [t] if t.is_file() else [p for p in t.rglob("*") if p.is_file()]
    return sorted(p for p in fs if p.name not in excl and not p.name.endswith("_rerun.json"))


def verify(out, r, rows, token):
    """Compare the run's outputs with the reference. Returns (report rows, gating_ok)."""
    rep, ok, seen = [], True, set()
    for x in rows:
        p = out / x["relpath"]
        seen.add(p.resolve())
        got = md5(p) if p.exists() else "MISSING"
        same = got == x["md5"]
        st = "PASS" if same else ("MISSING" if got == "MISSING" else "DIFF")
        if x["cls"] == "gating" and not same:
            ok = False
        rep.append(dict(run=r["run_name"], relpath=x["relpath"], cls=x["cls"], expected=x["md5"], obtained=got,
                        status=st))
    for p in output_files(out, r, token):
        if p.resolve() not in seen:
            c = cls_of(p.name)
            if c == "gating":
                ok = False
            rep.append(dict(run=r["run_name"], relpath=p.relative_to(out).as_posix(), cls=c, expected="",
                            obtained=md5(p), status="EXTRA"))
    return rep, ok


def swap_manifests(out, rows, refman):
    """After a gating PASS: keep the new manifest as *_rerun.json and restore the original."""
    n = 0
    for x in rows:
        if x["cls"] != "manifest":
            continue
        p, src = out / x["relpath"], refman / x["relpath"]
        if p.exists() and md5(p) == x["md5"]:
            continue
        if p.exists():
            p.rename(p.with_name(p.stem + "_rerun.json"))
        shutil.copy2(src, p)
        if md5(p) != x["md5"]:
            raise RuntimeError(f"restored manifest MD5 mismatch: {p}")
        n += 1
    return n


def merge_dir(tmp, dst):
    """Move files from a temporary output folder into dst. An existing identical file is kept (the temporary copy is
    removed); an existing different file is left untouched and the new file kept as <name>.rerun for the report."""
    conflicts = []
    for p in sorted(tmp.rglob("*")):
        if not p.is_file():
            continue
        d = dst / p.relative_to(tmp)
        d.parent.mkdir(parents=True, exist_ok=True)
        if not d.exists():
            os.replace(p, d)
        elif md5(d) == md5(p) or (d.name.lower().endswith(TEXT) and lf(d.read_bytes()) == lf(p.read_bytes())):
            p.unlink()
        else:
            os.replace(p, d.with_name(d.name + ".rerun"))
            conflicts.append(d.relative_to(dst.parent).as_posix())
    shutil.rmtree(tmp)
    return conflicts


def selftest():
    res = []
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        # rename of a folder run
        f = out / "P1_20990101_000000"; (f / "k01").mkdir(parents=True)
        (f / "P1_components_20990101_000000.csv").write_text("a,b\n1,2\n")
        (f / "k01" / "P1_x_20990101_000000.npz").write_bytes(b"x")
        (f / "P1_manifest_20990101_000000.json").write_text('{"t": "new"}')
        top = rename_tree(f, "20990101_000000", "20261111_233936")
        res.append(("folder renamed", top.name == "P1_20261111_233936"))
        res.append(("inner files renamed", (top / "k01" / "P1_x_20261111_233936.npz").exists()))
        # verification: one gating match, manifest differs (allowed), info differs (allowed)
        (top / "log.txt").write_text("time 1")
        rows = [dict(relpath="P1_20261111_233936/P1_components_20261111_233936.csv", cls="gating",
                     md5=md5(top / "P1_components_20261111_233936.csv")),
                dict(relpath="P1_20261111_233936/k01/P1_x_20261111_233936.npz", cls="gating", md5=md5(top / "k01" / "P1_x_20261111_233936.npz")),
                dict(relpath="P1_20261111_233936/P1_manifest_20261111_233936.json", cls="manifest", md5="0" * 32),
                dict(relpath="P1_20261111_233936/log.txt", cls="info", md5="1" * 32)]
        r = dict(run_name="P1_20261111_233936", kind="folder", exclude="")
        rep, ok = verify(out, r, rows, "20261111_233936")
        res.append(("gating pass despite manifest/info differences", ok))
        # a changed result file must fail
        (top / "P1_components_20261111_233936.csv").write_text("a,b\n1,3\n")
        rep, ok = verify(out, r, rows, "20261111_233936")
        res.append(("changed result detected", not ok and any(x["status"] == "DIFF" for x in rep)))
        # an extra result file must fail
        (top / "P1_components_20261111_233936.csv").write_text("a,b\n1,2\n")
        (top / "unexpected.csv").write_text("z")
        rep, ok = verify(out, r, rows, "20261111_233936")
        res.append(("extra result file detected", not ok))
        (top / "unexpected.csv").unlink()
        # manifest swap
        refman = out / "_refman" / "P1_20261111_233936"; refman.mkdir(parents=True)
        (refman / "P1_manifest_20261111_233936.json").write_text('{"t": "orig"}')
        rows[2]["md5"] = md5(refman / "P1_manifest_20261111_233936.json")
        n = swap_manifests(out, rows, out / "_refman")
        res.append(("manifest restored, rerun kept", n == 1 and (top / "P1_manifest_20261111_233936_rerun.json").exists()
                     and md5(top / "P1_manifest_20261111_233936.json") == rows[2]["md5"]))
        # flat run rename and refusal to overwrite
        (out / "O1_x_20990101_000001.csv").write_text("1")
        rename_tree(out / "O1_x_20990101_000001.csv", "20990101_000001", "20261111_232613")
        res.append(("flat file renamed", (out / "O1_x_20261111_232613.csv").exists()))
        (out / "O1_x_20990101_000002.csv").write_text("2")
        try:
            rename_tree(out / "O1_x_20990101_000002.csv", "20990101_000002", "20261111_232613"); refused = False
        except RuntimeError:
            refused = True
        res.append(("refuses to overwrite", refused))
        # merge of an argument folder: identical kept, different kept aside
        dst = out / "A3_u"; dst.mkdir(); (dst / "same.csv").write_text("s"); (dst / "diff.csv").write_text("d1")
        tmp = out / "_rerun_A3_u"; tmp.mkdir(); (tmp / "same.csv").write_text("s"); (tmp / "diff.csv").write_text("d2")
        (tmp / "new.csv").write_text("n")
        c = merge_dir(tmp, dst)
        res.append(("merge", c == ["A3_u/diff.csv"] and (dst / "new.csv").exists() and (dst / "diff.csv.rerun").exists()
                    and (dst / "diff.csv").read_text() == "d1" and not tmp.exists()))
        res.append(("ref token", ref_token(rows) == "20261111_233936"))
        # line endings: CRLF-only difference accepted and converted; any other difference rejected
        e = out / "E"; e.mkdir(); (e / "t.csv").write_bytes(b"a,b\r\n1,2\r\n"); (e / "u.csv").write_bytes(b"a,b\r\n1,3\r\n")
        er = [dict(relpath="E/t.csv", cls="gating", md5=hashlib.md5(b"a,b\n1,2\n").hexdigest()),
              dict(relpath="E/u.csv", cls="gating", md5=hashlib.md5(b"a,b\n1,2\n").hexdigest())]
        fx = fix_eol(out, er)
        res.append(("CRLF-only difference converted to the original bytes", fx == ["E/t.csv"]
                    and md5(e / "t.csv") == er[0]["md5"] and (e / "u.csv").read_bytes() == b"a,b\r\n1,3\r\n"))
        d2 = out / "M"; d2.mkdir(); (d2 / "k.csv").write_bytes(b"x\n")
        t2 = out / "_rerun_M"; t2.mkdir(); (t2 / "k.csv").write_bytes(b"x\r\n")
        res.append(("merge treats CRLF-only difference as identical", merge_dir(t2, d2) == [] and (d2 / "k.csv").read_bytes() == b"x\n"))
    res.append(("classes", cls_of("a.nii.gz") == "gating" and cls_of("X_manifest.json") == "manifest"
                and cls_of("f.tif") == "info"))
    for n_, v in res:
        print(f"  {'PASS' if v else 'FAIL'}  {n_}")
    print(f"selftest {sum(v for _, v in res)}/{len(res)}")
    return all(v for _, v in res)


# ------------------------------------------------------------------ execution ----------
def load():
    runs = list(csv.DictReader(open(RUNS, encoding="utf-8")))
    ref = {}
    for x in csv.DictReader(open(REF, encoding="utf-8")):
        ref.setdefault(x["run_name"], []).append(x)
    return runs, ref


def command(cfg, r, newdir):
    code, out = PROJECT / "02_code", PROJECT / "06_output"
    if r["lang"] == "py":
        a = [t.replace("{OUT}", str(out)).replace("{NEWDIR}", str(newdir)) for t in r["args"].split()]
        return [cfg["python"], "-u", str(code / r["script"])] + a
    fn = r["script"][:-2]
    call = fn if r["lang"] == "m_script" else fn + (r["args"] or "")
    return [cfg["matlab"], "-batch", f"cd('{code}'); {call}"]


def run_one(cfg, r, rows, logdir):
    out = PROJECT / "06_output"
    token = ref_token(rows)
    kind = r["kind"]
    if kind in ("folder",) and (out / r["run_name"]).exists():
        raise RuntimeError(f"{r['run_name']} already exists in 06_output; remove it or mark it as passed")
    if kind == "flat" and any(token in p.name for p in out.iterdir()):
        raise RuntimeError(f"outputs with time stamp {token} already exist in 06_output")
    newdir = out / f"_rerun_{r['run_name']}" if kind in ("argdir", "argfile") else None
    if newdir is not None:
        if newdir.exists():
            raise RuntimeError(f"{newdir} exists from an interrupted run; remove it first")
        newdir.mkdir()
    before = {p.name for p in out.iterdir()}
    cmd = command(cfg, r, newdir)
    t0 = time.time()
    with open(logdir / f"{int(r['order']):02d}_{r['run_name']}.log", "w", encoding="utf-8") as lg:
        lg.write(" ".join(cmd) + "\n\n"); lg.flush()
        rc = subprocess.run(cmd, cwd=PROJECT / "02_code", stdout=lg, stderr=subprocess.STDOUT).returncode
    secs = time.time() - t0
    if rc != 0:
        return dict(exit_code=rc, seconds=secs), [], False, []
    new = sorted(n for n in {p.name for p in out.iterdir()} - before)
    if newdir is not None:
        if token:
            for p in list(newdir.iterdir()):
                for t in set(TOKEN.findall(p.name)) - {token}:
                    rename_tree(p, t, token)
        conflicts = merge_dir(newdir, out / r["run_name"])
    else:
        conflicts = []
        if not new:
            return dict(exit_code=rc, seconds=secs, note="no new output found"), [], False, []
        toks = {t for n in new for t in TOKEN.findall(n)}
        if len(toks) != 1:
            raise RuntimeError(f"cannot identify the new run time stamp from {new}")
        nt = toks.pop()
        for n in new:
            rename_tree(out / n, nt, token)
    eol = set(fix_eol(out, rows))
    rep, ok = verify(out, r, rows, token)
    for x in rep:
        if x["relpath"] in eol and x["status"] == "PASS":
            x["status"] = "PASS_EOL"
    ok = ok and not conflicts
    return dict(exit_code=rc, seconds=secs, conflicts=";".join(conflicts), eol=len(eol)), rep, ok, conflicts


def dry_run(cfg, runs, ref):
    ok = True

    def chk(name, v, info=""):
        nonlocal ok
        ok &= bool(v)
        print(f"  {'OK ' if v else 'FAIL'}  {name} {info}")

    chk("Windows", os.name == "nt")
    chk("repository location", str(ROOT).lower() == str(PROJECT).lower(), str(ROOT))
    chk("dataset folder", DATA.is_dir(), str(DATA))
    if DATA.is_dir():
        lv1 = [p for p in DATA.iterdir() if p.is_dir()]
        n = len({p.name for q in [DATA] + lv1 for p in q.iterdir() if p.is_dir() and p.name.startswith("train_")})
        chk("dataset cases", n == 100, f"({n} train_ folders found; 100 expected)")
    chk("python executable", Path(cfg.get("python", "")).exists(), cfg.get("python", ""))
    chk("matlab executable", Path(cfg.get("matlab", "")).exists(), cfg.get("matlab", ""))
    if Path(cfg.get("python", "")).exists():
        v = subprocess.run([cfg["python"], "-c", "import sys,numpy,pandas,scipy,SimpleITK,matplotlib,sklearn;"
                            "print(sys.version.split()[0],numpy.__version__,pandas.__version__,scipy.__version__,"
                            "SimpleITK.Version_VersionString(),matplotlib.__version__,sklearn.__version__)"],
                           capture_output=True, text=True)
        print(f"        python numpy pandas scipy SimpleITK matplotlib scikit-learn: {v.stdout.strip() or v.stderr.strip()[-200:]}")
        print("        original: 3.11.15 2.4.6 3.1.0 1.17.1 2.5.6 3.11.2 1.9.1")
    # archived scripts against the script MD5 recorded in each original manifest
    for r in runs:
        p = ROOT / "02_code" / r["script"]
        if not p.exists():
            chk(f"script {r['script']}", False, "missing"); continue
        rec = None
        for x in ref.get(r["run_name"], []):
            if x["cls"] == "manifest":
                m = json.load(open(REF_MAN / x["relpath"], encoding="utf-8"))
                rec = m.get("script_md5") or m.get("script_MD5")
                if rec:
                    break
        if rec:
            chk(f"{int(r['order']):02d} {r['script']}", md5(p) == rec, "(MD5 = original manifest)")
        else:
            print(f"  --   {int(r['order']):02d} {r['script']} (no script MD5 recorded in the original manifest)")
    chk("reference manifests present", all((REF_MAN / x["relpath"]).exists()
                                           for v in ref.values() for x in v if x["cls"] == "manifest"))
    print("dry run:", "all checks passed" if ok else "FAILED - fix the items above")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list", action="store_true"); ap.add_argument("--only")
    ap.add_argument("--continue", dest="cont", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    runs, ref = load()
    if a.list:
        for r in runs:
            print(f"{int(r['order']):>3}  {r['run_name']:<40} {r['script']} {r['args']}")
        return
    cfg = json.load(open(ROOT / "reproduce_config.json", encoding="utf-8"))
    if a.dry_run:
        sys.exit(0 if dry_run(cfg, runs, ref) else 1)
    if not selftest() or not dry_run(cfg, runs, ref):
        sys.exit("STOP: self-test or dry run failed")
    lo, hi = (map(int, a.only.split("-")) if a.only and "-" in a.only else
              ((int(a.only), int(a.only)) if a.only else (1, len(runs))))
    wd = ROOT / "reproduction"; (wd / "logs").mkdir(parents=True, exist_ok=True)
    sf = wd / "state.json"
    state = json.load(open(sf)) if sf.exists() else {}
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    rep_f, sum_f = wd / f"reproduction_report_{stamp}.csv", wd / f"reproduction_summary_{stamp}.csv"
    fields = ["run", "relpath", "cls", "expected", "obtained", "status"]
    with open(rep_f, "w", newline="", encoding="utf-8") as fr, open(sum_f, "w", newline="", encoding="utf-8") as fs:
        wr = csv.DictWriter(fr, fieldnames=fields); wr.writeheader()
        ws = csv.writer(fs); ws.writerow(["order", "run", "script", "result", "seconds", "gating_pass", "gating_total",
                                          "info_diff", "manifests_restored", "note"])
        for r in runs:
            o = int(r["order"])
            if not lo <= o <= hi:
                continue
            if state.get(r["run_name"]) == "PASS":
                print(f"{o:>3} {r['run_name']}: already passed (skipped)")
                continue
            rows = ref[r["run_name"]]
            print(f"{o:>3} {r['run_name']}: running {r['script']} ...", flush=True)
            info, rep, ok, conf = run_one(cfg, r, rows, wd / "logs")
            for x in rep:
                wr.writerow(x)
            fr.flush()
            g = [x for x in rep if x["cls"] == "gating"]
            nrest = swap_manifests(PROJECT / "06_output", rows, REF_MAN) if ok else 0
            res = "PASS" if ok else ("ERROR" if info.get("exit_code") else "FAIL")
            ws.writerow([o, r["run_name"], r["script"], res, round(info.get("seconds", 0)),
                         sum(x["status"] in ("PASS", "PASS_EOL") for x in g), len(g),
                         sum(x["cls"] == "info" and x["status"] != "PASS" for x in rep), nrest,
                         info.get("note", "") or info.get("conflicts", "") or (f"exit code {info['exit_code']}" if info.get("exit_code") else "")])
            fs.flush()
            print(f"      {res}: {sum(x['status'] in ('PASS', 'PASS_EOL') for x in g)}/{len(g)} result files identical"
                  f"{' (' + str(info['eol']) + ' after LF line-ending conversion)' if info.get('eol') else ''}; "
                  f"{round(info.get('seconds', 0))} s")
            if ok:
                state[r["run_name"]] = "PASS"
                json.dump(state, open(sf, "w"), indent=1)
            elif not a.cont:
                print(f"STOP at run {o}. See {rep_f.name} and logs/{o:02d}_{r['run_name']}.log")
                sys.exit(1)
    n = sum(1 for r in runs if state.get(r["run_name"]) == "PASS")
    print(f"{n}/{len(runs)} runs reproduced byte for byte. Report: {rep_f}")


if __name__ == "__main__":
    main()
