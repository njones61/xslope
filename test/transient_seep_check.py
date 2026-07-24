"""Phase-1 acceptance locks for the transient (time-dependent) seepage solver
(plan_transient_seep.md §8.1-2).  The transient solver adds a storage term and a
theta-scheme time integrator to the steady FE seepage machinery, reusing the same
batched build-once assembly.  These checks lock its correctness:

  1. ANALYTICAL — a 1-D sudden-head-change column against the erfc solution.
     A thin 2-D saturated column, initially at uniform head, has one end stepped
     to a higher head at t=0 (a repeated-time step series -> the initial condition
     falls out of the plan's "IC = steady solve at t=0 BCs" rule).  Head-vs-time at
     an interior point must track h0 + dH*erfc(x / (2*sqrt(D t))), D = k/Ss, and the
     discrete mass balance must close.

     Tolerance basis: the numeric error combines a spatial part (the wetting front
     is ~sqrt(D t) wide; near t=0 it spans only a few elements, which caps accuracy
     there) and a temporal part (backward Euler is O(dt)).  Away from the earliest
     sharp-front frames the agreement is <0.1% of dH; the peak over ALL saved frames
     stays under 1% of dH.  Mass balance closes to ~1e-14 (lumped backward Euler
     conserves the discrete storage exactly), locked at 1e-8.

  2. STEADY LIMIT — a constant-BC transient run to large t must equal the steady
     solver on the same mesh.  Confined column: (A) the default IC (a steady solve
     at the t=0 BCs) matches the steady solver and holds; (B) a non-steady (uniform)
     IC relaxes to the steady solution to machine precision.  Unconfined earth dam:
     the exit-face-per-step path relaxes toward the steady solution (smoke test of
     the seepage-face active set under time stepping).

  3. STEADY UNTOUCHED — build_seep_data on a numeric-BC (steady) file adds only
     empty transient bindings and leaves the baked bc arrays unchanged, so the
     steady solve is bit-identical.  (The full steady regression is run_tests.py
     --seep; this is a fast structural guard.)

Plus unit checks of the series evaluator and the storage functions.

Run directly:  PYTHONPATH=. python3 test/transient_seep_check.py
"""

import io
import contextlib
import os

import numpy as np
from scipy.special import erfc

from xslope.seep import (run_transient_seepage, solve_confined, build_seep_data,
                         build_tseep_data, _eval_series, storage_capacity_vec,
                         storage_potential_vec)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _quiet(fn, *a, **k):
    """Run fn silently (the solvers print progress); return its result."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _quad_grid(L, H, nx, ny):
    """Structured quad4 grid on [0,L]x[0,H]; returns (nodes, elements)."""
    xs = np.linspace(0, L, nx + 1)
    ys = np.linspace(0, H, ny + 1)
    nodes = np.array([[x, yv] for yv in ys for x in xs], dtype=float)

    def nid(i, j):
        return j * (nx + 1) + i
    elements = np.array([[nid(i, j), nid(i+1, j), nid(i+1, j+1), nid(i, j+1)]
                         for j in range(ny) for i in range(nx)], dtype=int)
    return nodes, elements, xs, ys


def _column_seep_data(nodes, elements, k, Ss, Sy=0.1, head_bindings=None,
                      const_dirichlet=None):
    """Minimal one-material seep_data for a saturated column."""
    n = len(nodes)
    bc_type = np.zeros(n, dtype=int)
    bc_values = np.zeros(n)
    for mask, val in (const_dirichlet or []):
        bc_type[mask] = 1
        bc_values[mask] = val
    ne = len(elements)
    return dict(
        nodes=nodes, elements=elements, element_types=np.full(ne, 4),
        element_materials=np.ones(ne, dtype=int), bc_type=bc_type,
        bc_values=bc_values, flux_nodal=np.zeros(n),
        k1_by_mat=np.array([k]), k2_by_mat=np.array([k]), angle_by_mat=np.array([0.0]),
        kr0_by_mat=np.array([0.001]), h0_by_mat=np.array([-1.0]),
        unsat_by_mat=np.array([0]), vg_a_by_mat=np.array([0.0]),
        vg_n_by_mat=np.array([0.0]), ss_by_mat=np.array([Ss]),
        sy_by_mat=np.array([Sy]), head_series_bindings=head_bindings or [],
        flux_series_bindings=[], unit_weight=9.81)


# ---------------------------------------------------------------------------
# Lock 1 — analytical erfc column
# ---------------------------------------------------------------------------
def check_analytical_erfc():
    fails = []
    L, H, nx, ny = 3.0, 0.2, 150, 2
    k, Ss = 1.0, 1.0
    D = k / Ss
    h0_far, h_step = 10.0, 11.0
    dH = h_step - h0_far
    x_meas = 0.2
    duration = 0.1

    nodes, elements, xs, ys = _quad_grid(L, H, nx, ny)
    tol = (xs[1] - xs[0]) * 0.25
    near = np.abs(nodes[:, 0] - 0.0) <= tol
    far = np.abs(nodes[:, 0] - L) <= tol

    seep_data = _column_seep_data(
        nodes, elements, k, Ss,
        head_bindings=[{"mask": near, "series": "res"}],
        const_dirichlet=[(far, h0_far)])

    # repeated-time step: value h0_far at t=0 (-> steady IC = uniform h0_far),
    # h_step for t>0.
    ta = np.array([0.0, 0.0])
    va = np.array([h0_far, h_step])
    tseep_data = dict(times=ta, series={"res": (ta, va)}, duration=duration,
                      save_interval=None, save_times=[0.02, 0.05, 0.1],
                      stage_1=None, stage_2=None, breakpoints=[0.0])

    sol = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0,
                 max_head_change_frac=0.02, verbose=False)
    if not sol["converged"]:
        fails.append("erfc: transient run reported non-convergence")

    # interior measurement node on the bottom row
    i = int(np.argmin(np.abs(xs - x_meas)))
    mnode = i                                  # j=0 row -> global index i
    xm = nodes[mnode, 0]

    max_rel = 0.0
    for fr in sol["frames"]:
        t = fr["time"]
        if t <= 0:
            continue
        num = fr["head"][mnode]
        ana = h0_far + dH * erfc(xm / (2.0 * np.sqrt(D * t)))
        max_rel = max(max_rel, abs(num - ana) / dH)
    if max_rel > 0.01:
        fails.append(f"erfc: max head error {max_rel:.3%} of dH exceeds 1% tolerance")

    closure = sol["mass_balance"]["final_closure"]
    if closure > 1e-8:
        fails.append(f"erfc: mass-balance closure {closure:.2e} exceeds 1e-8")
    return fails


# ---------------------------------------------------------------------------
# Lock 2a — confined steady limit
# ---------------------------------------------------------------------------
def check_steady_limit_confined():
    fails = []
    L, H, nx, ny = 4.0, 0.4, 40, 2
    k, Ss = 1.0, 0.5
    hL, hR = 12.0, 10.0
    rng = hL - hR

    nodes, elements, xs, ys = _quad_grid(L, H, nx, ny)
    n = len(nodes)
    tol = (xs[1] - xs[0]) * 0.25
    left = np.abs(nodes[:, 0] - 0.0) <= tol
    right = np.abs(nodes[:, 0] - L) <= tol
    seep_data = _column_seep_data(nodes, elements, k, Ss,
                                  const_dirichlet=[(left, hL), (right, hR)])

    bcs = [(idx, seep_data["bc_values"][idx]) for idx in range(n)
           if seep_data["bc_type"][idx] == 1]
    ne = len(elements)
    h_steady, *_ = _quiet(solve_confined, nodes, elements, seep_data["bc_type"], bcs,
                          np.full(ne, k), np.full(ne, k), np.zeros(ne),
                          np.full(ne, 4))

    # characteristic time ~ Ss L^2 / k = 8; run well past it.
    tseep_data = dict(times=np.array([]), series={}, duration=60.0,
                      save_interval=None, save_times=[], stage_1=None,
                      stage_2=None, breakpoints=[])

    # (A) default IC = steady solve at t=0.
    solA = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0, verbose=False)
    d_ic = np.abs(np.asarray(solA["frames"][0]["head"]) - h_steady).max()
    d_end = np.abs(np.asarray(solA["frames"][-1]["head"]) - h_steady).max()
    if d_ic > 1e-9:
        fails.append(f"steady-limit(A): IC-from-steady differs from steady by {d_ic:.2e}")
    if d_end > 1e-8 * rng:
        fails.append(f"steady-limit(A): t_end drifted from steady by {d_end:.2e}")

    # (B) asymmetric (uniform = hR) IC -> relax to steady, with real net storage gain.
    h_uniform = np.full(n, hR)
    solB = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0,
                  h_init=h_uniform, max_head_change_frac=0.1, verbose=False)
    d_endB = np.abs(np.asarray(solB["frames"][-1]["head"]) - h_steady).max()
    if d_endB > 1e-8 * rng:
        fails.append(f"steady-limit(B): uniform IC failed to relax to steady "
                     f"(residual {d_endB:.2e}, {d_endB/rng:.2e} of range)")
    closureB = solB["mass_balance"]["final_closure"]
    if closureB > 1e-6:
        fails.append(f"steady-limit(B): mass-balance closure {closureB:.2e} exceeds 1e-6")
    return fails


# ---------------------------------------------------------------------------
# Lock 2b — unconfined steady limit (exit-face-per-step smoke test)
# ---------------------------------------------------------------------------
def check_steady_limit_unconfined():
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import run_seepage_analysis

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/seep/files/xslope_earth_dam1.xlsx")
    d = load_slope_data(path)
    for m in d["materials"]:
        m["Ss"] = 1e-3
        m["Sy"] = 0.2
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    target = (max(xs) - min(xs)) / 18.0          # coarse -> ~60 nodes, fast

    mesh = _quiet(build_mesh_from_polygons, polys, target, "tri3")
    seep_data = _quiet(build_seep_data, mesh, d)
    steady = _quiet(run_seepage_analysis, seep_data, tol=1e-6, max_iter=400)
    if not steady.get("converged", True):
        fails.append("unconfined: steady reference did not converge")
        return fails
    n = len(seep_data["nodes"])
    h_steady = np.asarray(steady["head"])
    rng = h_steady.max() - h_steady.min()

    fixed = seep_data["bc_values"][seep_data["bc_type"] == 1]
    h_uniform = np.full(n, float(np.mean(fixed)))
    d_ic = np.abs(h_uniform - h_steady).max() / rng

    tseep_data = dict(times=np.array([]), series={}, duration=4.0,
                      save_interval=None, save_times=[], stage_1=None,
                      stage_2=None, breakpoints=[])
    sol = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0,
                 h_init=h_uniform, max_head_change_frac=0.15, verbose=False)
    if not sol["converged"]:
        fails.append("unconfined: transient reported non-convergence")
    h_end = np.asarray(sol["frames"][-1]["head"])
    d_end = np.abs(h_end - h_steady).max() / rng
    # the exit-face-per-step run must relax substantially toward steady...
    if d_end > 0.05:
        fails.append(f"unconfined: t_end still {d_end:.2%} of range from steady "
                     f"(exit-face relaxation too far off)")
    # ...and be much closer than the initial condition was.
    if d_end > 0.5 * d_ic:
        fails.append(f"unconfined: transient barely relaxed (IC {d_ic:.2%} -> "
                     f"end {d_end:.2%} of range)")
    return fails


# ---------------------------------------------------------------------------
# Lock 3 — steady build_seep_data untouched (structural guard)
# ---------------------------------------------------------------------------
def check_steady_untouched():
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/seep/files/xslope_earth_dam1.xlsx")
    d = load_slope_data(path)
    if "tseep" in d:
        fails.append("steady file unexpectedly carries a tseep key")
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 20.0, "tri3")
    seep_data = _quiet(build_seep_data, mesh, d)
    # New transient keys must be present but INERT on a steady file.
    if seep_data.get("head_series_bindings"):
        fails.append("steady file produced non-empty head_series_bindings")
    if seep_data.get("flux_series_bindings"):
        fails.append("steady file produced non-empty flux_series_bindings")
    if "ss_by_mat" not in seep_data or "sy_by_mat" not in seep_data:
        fails.append("ss_by_mat/sy_by_mat missing from seep_data")
    # A numeric-BC file must still bake real Dirichlet nodes.
    if not np.any(seep_data["bc_type"] > 0):
        fails.append("steady file baked no Dirichlet/exit-face nodes")
    # build_tseep_data returns None for a steady (no-tseep) model.
    if build_tseep_data(d) is not None:
        fails.append("build_tseep_data returned non-None for a steady model")
    return fails


# ---------------------------------------------------------------------------
# Unit checks — series evaluator and storage functions
# ---------------------------------------------------------------------------
def check_series_eval():
    fails = []
    # linear interpolation between anchors, hold outside the range
    lin = (np.array([0.0, 10.0, 30.0]), np.array([160.0, 155.0, 140.0]))
    cases = [(-5.0, 160.0), (0.0, 160.0), (5.0, 157.5), (10.0, 155.0),
             (20.0, 147.5), (30.0, 140.0), (99.0, 140.0)]
    for t, exp in cases:
        got = _eval_series(lin, t)
        if abs(got - exp) > 1e-9:
            fails.append(f"series linear: t={t} got {got}, expected {exp}")
    # step function: repeated time, right-continuous (new value applies from t on)
    step = (np.array([0.0, 5.0, 5.0, 10.0]), np.array([1.0, 1.0, 9.0, 9.0]))
    for t, exp in [(0.0, 1.0), (2.5, 1.0), (4.999, 1.0), (5.0, 9.0), (7.0, 9.0)]:
        got = _eval_series(step, t)
        if abs(got - exp) > 1e-6:
            fails.append(f"series step: t={t} got {got}, expected {exp}")
    return fails


def check_build_tseep_data():
    fails = []
    # steady model (no tseep) -> None
    if build_tseep_data({}) is not None:
        fails.append("build_tseep_data({}) should be None")
    # a tseep dict as fileio produces it: series carry None gaps that must be
    # reduced to each series' own non-blank anchors; breakpoints = union of anchors.
    slope_data = {"tseep": {
        "times": [0.0, 10.0, 20.0, 30.0],
        "series": {"res": [160.0, None, 140.0, None],   # anchors at t=0, 20
                   "q": [None, 1e-6, None, 5e-6]},        # anchors at t=10, 30
        "duration": 50.0, "save_interval": 5.0,
        "save_times": [12.0, 34.0], "stage_1": 0.0, "stage_2": 30.0}}
    td = build_tseep_data(slope_data)
    if td is None:
        return ["build_tseep_data returned None for a transient model"]
    ta_res, va_res = td["series"]["res"]
    if not (np.array_equal(ta_res, [0.0, 20.0]) and np.array_equal(va_res, [160.0, 140.0])):
        fails.append(f"res anchors wrong: {ta_res}, {va_res}")
    ta_q, va_q = td["series"]["q"]
    if not (np.array_equal(ta_q, [10.0, 30.0]) and np.array_equal(va_q, [1e-6, 5e-6])):
        fails.append(f"q anchors wrong: {ta_q}, {va_q}")
    if td["breakpoints"] != [0.0, 10.0, 20.0, 30.0]:
        fails.append(f"breakpoints wrong: {td['breakpoints']}")
    if td["duration"] != 50.0 or td["stage_2"] != 30.0:
        fails.append("controls not carried through")
    # interpolation through the None gap is linear (identical to filling it in)
    if abs(_eval_series(td["series"]["res"], 10.0) - 150.0) > 1e-9:
        fails.append("res interpolation through a blank gap is not linear")
    return fails


def check_storage_functions():
    fails = []
    Ss, Sy, h0 = 1e-4, 0.2, -2.0
    # saturated -> Ss; in-band -> Ss + Sy/|h0|; below band -> Ss.
    C = storage_capacity_vec(np.array([1.0, -1.0, -3.0]),
                             np.full(3, Ss), np.full(3, Sy), np.full(3, h0))
    if abs(C[0] - Ss) > 1e-15:
        fails.append(f"storage(sat) = {C[0]} != Ss")
    if abs(C[1] - (Ss + Sy / abs(h0))) > 1e-12:
        fails.append(f"storage(band) = {C[1]} != Ss+Sy/|h0|")
    if abs(C[2] - Ss) > 1e-15:
        fails.append(f"storage(dry) = {C[2]} != Ss (floor)")
    if np.any(C <= 0):
        fails.append("storage capacity is not strictly positive")
    # Phi is the integral of S: dPhi/dpsi ~ S by finite difference in the band.
    dp = 1e-6
    p = -1.0
    phi1 = storage_potential_vec(np.array([p + dp]), np.array([Ss]),
                                 np.array([Sy]), np.array([h0]))[0]
    phi0 = storage_potential_vec(np.array([p - dp]), np.array([Ss]),
                                 np.array([Sy]), np.array([h0]))[0]
    dPhi = (phi1 - phi0) / (2 * dp)
    Cmid = Ss + Sy / abs(h0)
    if abs(dPhi - Cmid) / Cmid > 1e-4:
        fails.append(f"dPhi/dpsi={dPhi} != S={Cmid} in the band")
    # vG capacity: positive and finite for n>1.
    Cvg = storage_capacity_vec(np.array([-1.0]), np.array([Ss]), np.array([Sy]),
                               np.array([0.0]), vg_a=np.array([0.5]),
                               vg_n=np.array([2.0]), model=np.array([1]))
    if not (np.isfinite(Cvg[0]) and Cvg[0] > Ss):
        fails.append(f"vG storage capacity {Cvg[0]} not finite/> Ss")
    return fails


def main():
    print("transient seepage phase-1 locks:")
    checks = [
        ("series evaluator", check_series_eval),
        ("build_tseep_data", check_build_tseep_data),
        ("storage functions", check_storage_functions),
        ("steady untouched", check_steady_untouched),
        ("analytical erfc", check_analytical_erfc),
        ("steady limit (confined)", check_steady_limit_confined),
        ("steady limit (unconfined)", check_steady_limit_unconfined),
    ]
    failures = []
    for name, fn in checks:
        fs = fn()
        print(f"  {name:28s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll transient seepage phase-1 locks passed.")


if __name__ == "__main__":
    main()
