"""
Tau_validation.py - PAPER 3, Step 2 SAP: matching distance and reference tolerance tau
(spec StepSAP_tau_matching_tolerance_spec_DRAFT_v2_20261111.md): grid audit (tau4) and synthetic validation T1-T9
(tau5). Synthetic data and grid metadata only; no PET, tumour or correspondence data are read.

Usage:  python Tau_validation.py <Audit1D_grid_geometry_20261030_215034.csv> <out_dir>
Outputs (out_dir): Tau_grid_audit_<run>.csv, Tau_synthetic_tests_<run>.csv, Tau_T6_characterisation_<run>.csv,
                   Tau_T9_characterisation_<run>.csv, Tau_validation_manifest_<run>.json
Exit code 0 only if every pass/fail test (T1-T5, T7, T8) passes.

Geometry (spec tau2): a grid is (size, spacing, origin, direction); physical position of voxel index (i, j, k) =
origin + direction @ ((i, j, k) * spacing). Boundary voxel = labelled voxel with >= 1 of its 6 face neighbours
outside the label (grid edge counts as outside). d = 0 if the two labels share any voxel (D2); otherwise the minimum
Euclidean distance (mm) between boundary-voxel centres (D1), by k-d tree. tau_i = sqrt(dx^2 + dy^2 + dz^2) of the
matching (PSMA PET) grid; matched if d <= tau_i. Sensitivity tolerances: strict = max spacing; permissive =
2 x max spacing.
"""
import hashlib
import itertools
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

EXPECTED_GRID_MD5 = "8c8ac873aedca4114ccf033a15ab7efc"
SEED_T9 = 20261111
N_SHIFTS_T9 = 50
TOL = 1e-9


# ------------------------------ geometry ------------------------------
class Grid:
    def __init__(self, size, spacing, origin=(0.0, 0.0, 0.0), direction=np.eye(3)):
        self.size = tuple(int(s) for s in size)               # (nx, ny, nz); arrays are indexed [i, j, k]
        self.spacing = np.asarray(spacing, float)
        self.origin = np.asarray(origin, float)
        self.direction = np.asarray(direction, float)

    def phys(self, ijk):
        return self.origin + (np.asarray(ijk, float) * self.spacing) @ self.direction.T

    def index(self, p):
        """continuous index of physical points"""
        return np.linalg.solve(self.direction * self.spacing, (np.atleast_2d(p) - self.origin).T).T


def tau_ref(sp):
    return float(np.sqrt(np.sum(np.square(sp))))


def tau_strict(sp):
    return float(np.max(sp))


def tau_permissive(sp):
    return float(2 * np.max(sp))


def boundary(mask):
    """6-connected boundary voxels (grid edge counts as outside)."""
    m = np.pad(mask.astype(bool), 1, constant_values=False)
    inner = m[1:-1, 1:-1, 1:-1]
    nb_all = np.ones_like(inner)
    for ax in range(3):
        for s in (1, -1):
            nb_all &= np.roll(m, s, axis=ax)[1:-1, 1:-1, 1:-1]
    return inner & ~nb_all


def match_distance(fdg, psma, grid, brute=False):
    """D2 then D1. fdg, psma: boolean arrays on the same grid. Returns d in mm."""
    if np.any(fdg & psma):
        return 0.0
    a = grid.phys(np.argwhere(boundary(fdg)))
    b = grid.phys(np.argwhere(boundary(psma)))
    if len(a) == 0 or len(b) == 0:
        return np.inf
    if brute:
        return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())
    return float(cKDTree(b).query(a, k=1)[0].min())


def correspondences(fdg_labels, psma_labels, grid, tau):
    """all (fdg_id, psma_id, d) pairs with d <= tau (one-to-many and many-to-one)."""
    out = []
    for f in sorted(set(np.unique(fdg_labels)) - {0}):
        for p in sorted(set(np.unique(psma_labels)) - {0}):
            d = match_distance(fdg_labels == f, psma_labels == p, grid)
            if d <= tau + TOL:
                out.append((int(f), int(p), d))
    return out


# ------------------------------ synthetic scenes ------------------------------
def scene(grid, f_boxes, p_boxes):
    F = np.zeros(grid.size, bool); P = np.zeros(grid.size, bool)
    for (lo, hi) in f_boxes:
        F[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = True
    for (lo, hi) in p_boxes:
        P[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = True
    return F, P


def cases(sp):
    """(test, name, F boxes, P boxes, expected d) on a 30^3 grid; F is a 3^3 block at [10,13)."""
    f = [((10, 10, 10), (13, 13, 13))]
    dx, dy, dz = sp
    out = [("T1", "overlap (shared voxels)", f, [((12, 12, 12), (16, 16, 16))], 0.0),
           ("T2", "containment (F inside larger P)", f, [((6, 6, 6), (18, 18, 18))], 0.0)]
    for ax, nm, s in ((0, "x", dx), (1, "y", dy), (2, "z", dz)):
        lo = [10, 10, 10]; hi = [13, 13, 13]; lo[ax] = 13; hi[ax] = 16
        out.append(("T3", f"face contact {nm}", f, [(tuple(lo), tuple(hi))], s))
    for (a, b), nm in (((0, 1), "xy"), ((0, 2), "xz"), ((1, 2), "yz")):
        lo = [10, 10, 10]; hi = [13, 13, 13]
        for ax in (a, b):
            lo[ax] = 13; hi[ax] = 16
        out.append(("T4", f"edge contact {nm}", f, [(tuple(lo), tuple(hi))], float(np.hypot(sp[a], sp[b]))))
    out.append(("T5", "corner contact", f, [((13, 13, 13), (16, 16, 16))], tau_ref(sp)))
    for ax, nm, s in ((0, "x", dx), (1, "y", dy), (2, "z", dz)):
        lo = [10, 10, 10]; hi = [13, 13, 13]; lo[ax] = 14; hi[ax] = 17
        out.append(("T6", f"one-voxel gap {nm}", f, [(tuple(lo), tuple(hi))], 2 * s))
    out.append(("T7", "one-voxel diagonal gap", f, [((14, 14, 14), (17, 17, 17))], 2 * tau_ref(sp)))
    return out


def orientations(sp):
    """identity, axis flips and a shifted origin (LPS variants)."""
    yield "identity, origin 0", np.eye(3), (0.0, 0.0, 0.0)
    for flips in itertools.product((1, -1), repeat=3):
        if flips != (1, 1, 1):
            yield f"flip {flips}, origin shifted", np.diag(flips).astype(float), (-312.7, 145.3, -988.1)


# ------------------------------ T9: nearest-neighbour resampling from a coarser grid ------------------------------
def resample_nn(src_mask, src, dst):
    """nearest-neighbour resampling of a mask from grid src onto grid dst (identity transform)."""
    idx = np.indices(dst.size).reshape(3, -1).T
    ci = np.rint(src.index(dst.phys(idx))).astype(int)
    ok = np.all((ci >= 0) & (ci < np.array(src.size)), axis=1)
    out = np.zeros(int(np.prod(dst.size)), bool)
    out[ok] = src_mask[ci[ok, 0], ci[ok, 1], ci[ok, 2]]
    return out.reshape(dst.size)


def t9(coarse_sp, fine_sp, rng):
    """FDG block on the native FDG grid ("coarse"); PSMA region drawn on the PSMA grid ("fine") as the coarse block's face neighbour (contact)
    or one coarse voxel further (gap). Returns shares of shifts with d <= tau(fine)."""
    res = {"contact": [], "gap_one_coarse_voxel": []}
    for _ in range(N_SHIFTS_T9):
        shift = rng.uniform(0, 1, 3) * np.asarray(fine_sp)
        C = Grid((16, 16, 16), coarse_sp)
        Fc = np.zeros(C.size, bool); Fc[5:9, 5:9, 5:9] = True
        D = Grid(tuple(int(np.ceil(16 * c / f)) for c, f in zip(coarse_sp, fine_sp)), fine_sp, origin=shift)
        Ff = resample_nn(Fc, C, D)
        for kind, off in (("contact", 9), ("gap_one_coarse_voxel", 10)):
            # PSMA region = coarse voxels [off, off+3) along x, rasterised on the fine grid by voxel-centre membership
            ci = C.index(D.phys(np.indices(D.size).reshape(3, -1).T))
            r = np.rint(ci)
            P = ((r[:, 0] >= off) & (r[:, 0] < off + 3) & (r[:, 1] >= 5) & (r[:, 1] < 9)
                 & (r[:, 2] >= 5) & (r[:, 2] < 9)).reshape(D.size)
            res[kind].append(match_distance(Ff, P, D) <= tau_ref(fine_sp) + TOL)
    return {k: float(np.mean(v)) for k, v in res.items()}


# ------------------------------ run ------------------------------
def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main(grid_csv, out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    run = datetime.now().strftime("%Y%m%d_%H%M%S")
    gm = md5(grid_csv)
    if gm != EXPECTED_GRID_MD5:
        sys.exit(f"STOP C1: grid metadata MD5 {gm} != {EXPECTED_GRID_MD5}")
    G = pd.read_csv(grid_csv)
    P = G[G.tracer == "PSMA"].copy()
    if len(P) != 100 or P.case_id.nunique() != 100:
        sys.exit("STOP: grid metadata does not contain 100 PSMA rows")
    sp = P[["spacing_x_or_ratio", "spacing_y_or_ratio", "spacing_z_or_ratio"]].to_numpy(float)
    # ---- tau4 grid audit ----
    A = pd.DataFrame(dict(case_id=P.case_id, dx=sp[:, 0], dy=sp[:, 1], dz=sp[:, 2]))
    A["tau_strict_mm"] = sp.max(1); A["tau_ref_mm"] = np.sqrt((sp ** 2).sum(1)); A["tau_permissive_mm"] = 2 * sp.max(1)
    A["min_spacing_mm"] = sp.min(1); A["ratio_tau_over_2min"] = A.tau_ref_mm / (2 * A.min_spacing_mm)
    A["gap_along_min_axis_matched"] = A.tau_ref_mm >= 2 * A.min_spacing_mm - TOL
    A["obliquity_deg"] = P.obliquity_deg.to_numpy(); A["axis_code"] = P.axis_code_or_flag.to_numpy()
    A["PSMA_FDG_grid_differs"] = [bool(G[(G.case_id == c) & (G.tracer == "PAIR_FDG_over_PSMA")]
                                       [["spacing_x_or_ratio", "spacing_y_or_ratio", "spacing_z_or_ratio"]]
                                       .sub(1).abs().to_numpy().max() > 1e-4) for c in P.case_id]   # 1e-4: ignores float noise (train_0013 1.000002)
    if (A.obliquity_deg.abs() > 0).any():
        sys.exit("STOP: oblique PSMA PET grid present (protocol 14.4 condition not met)")
    ap = out / f"Tau_grid_audit_{run}.csv"; A.to_csv(ap, index=False)
    print(f"tau4 grid audit: 100 PSMA grids, {len(np.unique(sp.round(6), axis=0))} distinct spacings; obliquity 0 in all;"
          f" tau_ref median {A.tau_ref_mm.median():.2f} mm (P10-P90 {A.tau_ref_mm.quantile(.1):.2f}-"
          f"{A.tau_ref_mm.quantile(.9):.2f}; range {A.tau_ref_mm.min():.2f}-{A.tau_ref_mm.max():.2f}); "
          f"one-voxel gap along smallest axis matched in {int(A.gap_along_min_axis_matched.sum())}/100; "
          f"PSMA and FDG grids differ in {int(A.PSMA_FDG_grid_differs.sum())}")
    # ---- tau5 synthetic tests ----
    grids = [(4.0, 4.0, 4.0), (2.73, 2.73, 5.0), (4.07, 4.07, 1.5)] + [tuple(r) for r in np.unique(sp.round(6), axis=0)]
    rows, t6 = [], []
    for g in grids:
        T = tau_ref(g)
        for test, name, fb, pb, exp in cases(g):
            for oname, D, O in orientations(g):
                grid = Grid((30, 30, 30), g, O, D)
                F, Pm = scene(grid, fb, pb)
                d = match_distance(F, Pm, grid)
                db = match_distance(F, Pm, grid, brute=True)
                matched = d <= T + TOL
                if test in ("T1", "T2", "T3", "T4", "T5"):
                    ok = abs(d - exp) <= TOL and matched
                elif test == "T7":
                    ok = abs(d - exp) <= TOL and not matched
                else:                                         # T6 characterisation: distance must still be exact
                    ok = abs(d - exp) <= TOL
                rows.append(dict(grid=str(g), test=test, case=name, orientation=oname, expected_d_mm=exp, d_mm=d,
                                 d_brute_mm=db, tau_ref_mm=T, matched=matched,
                                 role="characterisation" if test == "T6" else "pass/fail",
                                 T8a_knn_equals_brute=abs(d - db) <= TOL, PASS=bool(ok and abs(d - db) <= TOL)))
                if test == "T6" and oname.startswith("identity"):
                    t6.append(dict(grid=str(g), axis=name[-1], gap_d_mm=d, tau_ref_mm=T, matched_at_ref=matched,
                                   matched_at_strict=d <= tau_strict(g) + TOL,
                                   matched_at_permissive=d <= tau_permissive(g) + TOL))
    R = pd.DataFrame(rows)
    # T8 (b) orientation/origin invariance; (c) one-to-many / many-to-one; (d) mm not voxel units
    inv = R.groupby(["grid", "test", "case"]).d_mm.agg(lambda s: s.max() - s.min()).max()
    g = Grid((30, 30, 30), (2.73, 2.73, 5.0))
    FL = np.zeros(g.size, int); PL = np.zeros(g.size, int)
    FL[10:13, 10:13, 10:13] = 1; FL[10:13, 18:21, 10:13] = 2
    PL[13:16, 10:13, 10:13] = 1; PL[10:13, 7:10, 10:13] = 2; PL[13:16, 18:21, 10:13] = 3
    pairs = {(f, p) for f, p, _ in correspondences(FL, PL, g, tau_ref(g.spacing))}
    t8c = pairs == {(1, 1), (1, 2), (2, 3)}
    Fz, Pz = scene(g, [((10, 10, 10), (13, 13, 13))], [((10, 10, 13), (13, 13, 16))])
    t8d = abs(match_distance(Fz, Pz, g) - 5.0) <= TOL
    extra = [dict(test="T8", case="(a) knn = brute force, all cases", PASS=bool(R.T8a_knn_equals_brute.all())),
             dict(test="T8", case=f"(b) orientation/origin invariance (max spread {inv:.1e} mm)", PASS=bool(inv <= TOL)),
             dict(test="T8", case="(c) one-to-many and many-to-one recorded", PASS=bool(t8c)),
             dict(test="T8", case="(d) distance in mm on anisotropic grid (z face = 5.0 mm)", PASS=bool(t8d))]
    R = pd.concat([R, pd.DataFrame(extra)], ignore_index=True)
    R["role"] = R["role"].fillna("pass/fail")
    rp = out / f"Tau_synthetic_tests_{run}.csv"; R.to_csv(rp, index=False)
    T6 = pd.DataFrame(t6); t6p = out / f"Tau_T6_characterisation_{run}.csv"; T6.to_csv(t6p, index=False)
    # ---- T9 ----
    rng = np.random.default_rng(SEED_T9)
    t9rows = []
    # FDG coarser than PSMA (as in most of the 21 patients) and, for completeness, FDG finer than PSMA
    for c, f in (((4.07, 4.07, 5.0), (2.73, 2.73, 3.0)), ((5.47, 5.47, 5.0), (4.07, 4.07, 3.0)),
                 ((2.73, 2.73, 3.0), (4.07, 4.07, 5.0))):
        r = t9(c, f, rng)
        t9rows.append(dict(coarse=str(c), fine=str(f), n_shifts=N_SHIFTS_T9, tau_fine_mm=tau_ref(f),
                           share_contact_kept=r["contact"], share_gap_matched=r["gap_one_coarse_voxel"]))
    T9 = pd.DataFrame(t9rows); t9p = out / f"Tau_T9_characterisation_{run}.csv"; T9.to_csv(t9p, index=False)
    # ---- summary ----
    pf = R[R.role == "pass/fail"]
    summ = pf.groupby("test").PASS.agg(["sum", "count"])
    print("pass/fail tests:\n" + summ.to_string())
    print(f"T6 (characterisation): one-voxel gaps matched at tau_ref in {int(T6.matched_at_ref.sum())}/{len(T6)} "
          f"axis x grid cases")
    print("T9 (characterisation):\n" + T9.to_string(index=False))
    allpass = bool(pf.PASS.all())
    man = dict(step="tau validation (grid audit + synthetic T1-T9)", run=run, script=Path(__file__).name,
               script_md5=md5(__file__), inputs={Path(grid_csv).name: gm},
               outputs={p.name: md5(p) for p in (ap, rp, t6p, t9p)}, all_pass_fail_tests_passed=allpass,
               params=dict(SEED_T9=SEED_T9, N_SHIFTS_T9=N_SHIFTS_T9, TOL=TOL),
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=__import__("scipy").__version__))
    mp = out / f"Tau_validation_manifest_{run}.json"; mp.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"ALL PASS/FAIL TESTS {'PASSED' if allpass else 'NOT PASSED'}; manifest {mp.name} MD5 {md5(mp)}")
    return allpass


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
