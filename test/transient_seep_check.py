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
     the exit-face-per-step path relaxes toward the steady solution on BOTH a linear
     (tri3) and a quadratic (tri6) mesh — tri3 exercises the per-corner active-set
     rule, tri6 the ported quadratic all-or-nothing edge rule (midside exit nodes
     resolved per step).  Linear-vs-quadratic consistency (lock 2c): the SAME
     reservoir-drawdown problem solved on tri3 and tri6 agrees on the drawn-down
     saturated pore-pressure field (RMS ~1.9% of head range) and lands the downstream
     seepage-face exit point at the same elevation (0.00 ft delta), within coarse-mesh
     refinement expectations.

  3. STEADY UNTOUCHED — build_seep_data on a numeric-BC (steady) file adds only
     empty transient bindings and leaves the baked bc arrays unchanged, so the
     steady solve is bit-identical.  (The full steady regression is run_tests.py
     --seep; this is a fast structural guard.)

Plus unit checks of the series evaluator and the storage functions.

Phase-3 locks (plan §4 / §6 — outputs + drawdown):

  4. FRAMES ROUND TRIP — a small transient run exports {base}_tseep.csv (+ meta),
     which reloads to the exact head/u/phi it stored while DERIVING velocity and
     gradient on load, and one reloaded frame renders through the UNCHANGED
     plot_seep_solution (contours + flow lines + vectors) without raising. Every
     saved frame carries a real phi (solver product) and boundary inflow/outflow.

  5. DRAWDOWN STAGING — the in-memory §6 path (stage_transient_for_drawdown)
     reproduces the classic two-steady-file rapid-drawdown FS to machine zero when
     the two staged frames carry the same pore-pressure fields — proving the
     frame -> seep_u/seep_u2 plumbing is exact (physics consistency against a real
     transient model is the phase-5 corpus lock).

Phase-5 lock (plan §8.4 — the feature's headline integration test):

  6. DRAWDOWN CONSISTENCY — a SLOW (quasi-steady) transient reservoir drawdown on
     the earth-dam rapid-drawdown model must reproduce the classic
     two-independent-steady-solve rapid-drawdown FS. Unlike lock 5 (which fed the
     SAME field to both stages to prove the plumbing), this runs REAL transient
     physics: the upstream reservoir head is bound to a time series that holds
     high, ramps down to the low pool, then holds low; the stage-1 tag (t=0) and
     stage-2 tag (well after the drop) drive rapid_drawdown from the transient
     frames, and the FS is compared against two independent steady seepage solves
     (one at each pool level) fed through the identical rapid-drawdown machinery.

     Both paths share the stage-1 field to machine zero (the transient IC IS a
     steady solve at the t=0 pool level), so the ENTIRE FS gap is the stage-2
     residual: a quasi-steady drawdown only approaches the low-pool steady state
     asymptotically. The tolerance is derived, not guessed — the same slow run is
     tagged at three growing post-drawdown equilibration times and the FS gap is
     shown to shrink and plateau:

         equilibration (t units)      |FS_transient - FS_steady|
             2                              ~1.7e-4
             4                              ~1.5e-4
             8                              ~5e-5   (demonstrated plateau)

     The longest-equilibration agreement is locked at 5e-4 — an order of magnitude
     above the demonstrated ~5e-5 residual (absorbing mesh/step jitter) and ~30x
     below the discriminator: a FAST drawdown (short drop, no hold) leaves stage-2
     far from steady and misses the steady FS by ~1.5e-2, proving the lock is not
     vacuously satisfied. Kept fast with a coarse tri3 mesh (~60 nodes, linear per
     the exit-face caveat) and serial solves (~5 s).

Run directly:  PYTHONPATH=. python3 test/transient_seep_check.py
"""

import io
import contextlib
import os

import numpy as np
from scipy.special import erfc

from xslope.seep import (run_transient_seepage, solve_confined, build_seep_data,
                         build_tseep_data, _eval_series, storage_capacity_vec,
                         storage_potential_vec, export_transient_solution,
                         import_transient_solution, stage_transient_for_drawdown,
                         select_transient_frame_u, transient_frame_index)

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
#
# The per-step exit-face active set must relax an unconfined earth dam toward the
# steady solution on the SAME mesh, at any element order. tri3 exercises the
# per-corner rule; tri6 exercises the ported quadratic all-or-nothing edge rule
# (midside nodes resolved per step), so both orders are locked here.
# ---------------------------------------------------------------------------
def _unconfined_relax(element_type, target_div, duration=4.0, frac=0.15):
    """Build the earth-dam unconfined seepage problem at a given element order,
    relax a uniform initial head toward the steady solution on that same mesh, and
    return (d_ic, d_end, transient_converged, steady_converged) where d_ic/d_end
    are the initial and final max head deviations from steady as fractions of the
    head range."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import run_seepage_analysis

    path = os.path.join(_REPO_ROOT, "docs/seep/files/xslope_earth_dam1.xlsx")
    d = load_slope_data(path)
    for m in d["materials"]:
        m["Ss"] = 1e-3
        m["Sy"] = 0.2
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    target = (max(xs) - min(xs)) / target_div

    mesh = _quiet(build_mesh_from_polygons, polys, target, element_type)
    seep_data = _quiet(build_seep_data, mesh, d)
    steady = _quiet(run_seepage_analysis, seep_data, tol=1e-6, max_iter=600)
    if not steady.get("converged", True):
        return None, None, None, False
    n = len(seep_data["nodes"])
    h_steady = np.asarray(steady["head"])
    rng = h_steady.max() - h_steady.min()

    fixed = seep_data["bc_values"][seep_data["bc_type"] == 1]
    h_uniform = np.full(n, float(np.mean(fixed)))
    d_ic = np.abs(h_uniform - h_steady).max() / rng

    tseep_data = dict(times=np.array([]), series={}, duration=duration,
                      save_interval=None, save_times=[], stage_1=None,
                      stage_2=None, breakpoints=[])
    sol = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0,
                 h_init=h_uniform, max_head_change_frac=frac, verbose=False)
    h_end = np.asarray(sol["frames"][-1]["head"])
    d_end = np.abs(h_end - h_steady).max() / rng
    return d_ic, d_end, bool(sol["converged"]), True


def _check_unconfined_relax(element_type, target_div):
    fails = []
    d_ic, d_end, conv, steady_ok = _unconfined_relax(element_type, target_div)
    if not steady_ok:
        return [f"unconfined ({element_type}): steady reference did not converge"]
    if not conv:
        fails.append(f"unconfined ({element_type}): transient reported non-convergence")
    # the exit-face-per-step run must relax substantially toward steady...
    if d_end > 0.05:
        fails.append(f"unconfined ({element_type}): t_end still {d_end:.2%} of range "
                     f"from steady (exit-face relaxation too far off)")
    # ...and be much closer than the initial condition was.
    if d_end > 0.5 * d_ic:
        fails.append(f"unconfined ({element_type}): transient barely relaxed "
                     f"(IC {d_ic:.2%} -> end {d_end:.2%} of range)")
    return fails


def check_steady_limit_unconfined():
    # tri3: per-corner exit-face rule (linear elements).
    return _check_unconfined_relax("tri3", 18.0)


def check_steady_limit_unconfined_tri6():
    # tri6: quadratic all-or-nothing edge rule ported from the steady solver.
    # Coarser division (fewer curved elements) keeps ~the same node count as the
    # tri3 case and the run fast, while genuinely exercising midside exit nodes.
    return _check_unconfined_relax("tri6", 10.0)


# ---------------------------------------------------------------------------
# Lock 2c — linear-vs-quadratic drawdown consistency
#
# The SAME unconfined reservoir-drawdown problem solved on a linear (tri3) and a
# quadratic (tri6) mesh must agree on the drawn-down pore-pressure field and the
# downstream seepage-face exit point, within coarse-mesh refinement expectations.
# This is the integration test that the ported quadratic exit-face edge rule
# behaves like the trusted per-corner linear rule: both meshes are marched through
# the SAME reservoir series (hold high -> ramp to low -> hold), and the equilibrated
# low-pool frame is compared.
#
# Achieved on the coarse pair used here (tri3 ~55 nodes / tri6 ~105 nodes, ~10 s):
#   * downstream exit-point elevation:  identical to 0.00 ft (both land on 227.00,
#     the tailwater/toe elevation) — the seepage face terminates at the same place;
#   * saturated-zone head field:  RMS ~1.9% of the head range, max ~9% (the max sits
#     on the phreatic surface, where a small head error moves the h=y crossing the
#     most and where the two meshes differ by genuine refinement).
# The head field is compared only where BOTH solutions are saturated (pressure head
# >= 0): that is where the pore pressures driving stability live, and it excludes the
# unsaturated cap where head is legitimately mesh-sensitive. Tolerances carry real
# margin over the achieved deltas.
# ---------------------------------------------------------------------------
def check_drawdown_linear_vs_quadratic():
    import copy
    from scipy.interpolate import LinearNDInterpolator
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/lem/files/xslope_earth_dam_rapid.xlsx")
    if not os.path.exists(path):
        return [f"rapid-drawdown model not found: {path}"]
    H1, H2 = 302.0, 250.0          # high / low pool

    def _build(element_type, target_div):
        d = load_slope_data(path)
        for m in d["materials"]:
            m["Ss"], m["Sy"] = 1e-4, 0.2
        polys = get_material_polygons(d)
        xs = [x for x, _ in d["ground_surface"].coords]
        width = max(xs) - min(xs)
        mesh = _quiet(build_mesh_from_polygons, polys, width / target_div, element_type)
        d_ts = copy.deepcopy(d)
        for sh in d_ts["seepage_bc"]["specified_heads"]:
            sh["head"] = "res"
        seep_data = _quiet(build_seep_data, mesh, d_ts)
        return seep_data, list(d["ground_surface"].coords)

    def _drawdown(seep_data):
        t_pre, t_dd, hold = 0.3, 3.0, 8.0     # hold long enough to equilibrate
        t_end, dur = t_pre + t_dd, t_pre + t_dd + hold
        ta = np.array([0.0, t_pre, t_end, dur])
        tseep = dict(times=ta, series={"res": (ta, np.array([H1, H1, H2, H2]))},
                     duration=dur, save_interval=dur, save_times=[dur],
                     stage_1=0.0, stage_2=dur, breakpoints=list(ta))
        sol = _quiet(run_transient_seepage, seep_data, tseep, theta=1.0,
                     dt_max=0.5, max_head_change_frac=0.15, verbose=False)
        h = np.asarray(sol["frames"][transient_frame_index(sol, dur)]["head"])
        return sol, h

    sd3, ground = _build("tri3", 14.0)       # linear: per-corner rule
    sd6, _ = _build("tri6", 8.0)             # quadratic: ported edge rule
    if not sd6.get("head_series_bindings"):
        return ["reservoir head binding was not created (check seepage_bc mapping)"]

    sol3, h3 = _drawdown(sd3)
    sol6, h6 = _drawdown(sd6)
    if not sol3["converged"]:
        fails.append("linear (tri3) drawdown reported non-convergence")
    if not sol6["converged"]:
        fails.append("quadratic (tri6) drawdown reported non-convergence")

    rng = float(h3.max() - h3.min())

    # (1) saturated-zone head agreement (pressure head >= 0 in BOTH solutions).
    p3f = LinearNDInterpolator(sd3["nodes"], h3 - sd3["nodes"][:, 1])
    p6f = LinearNDInterpolator(sd6["nodes"], h6 - sd6["nodes"][:, 1])
    h3f = LinearNDInterpolator(sd3["nodes"], h3)
    h6f = LinearNDInterpolator(sd6["nodes"], h6)
    n = sd3["nodes"]
    gx = np.linspace(n[:, 0].min(), n[:, 0].max(), 60)
    gy = np.linspace(n[:, 1].min(), n[:, 1].max(), 45)
    GX, GY = np.meshgrid(gx, gy)
    pts = np.column_stack([GX.ravel(), GY.ravel()])
    p3, p6 = p3f(pts), p6f(pts)
    sat = np.isfinite(p3) & np.isfinite(p6) & (p3 >= 0) & (p6 >= 0)
    if sat.sum() < 100:
        fails.append(f"too few common saturated sample points ({int(sat.sum())})")
    diff = np.abs(h3f(pts) - h6f(pts))[sat]
    rms = float(np.sqrt((diff ** 2).mean()))
    dmax = float(diff.max())
    # RMS is the robust agreement measure; max is a loose gross-divergence guard.
    if rms > 0.05 * rng:
        fails.append(f"linear-vs-quadratic saturated head RMS {rms:.3f} "
                     f"({100*rms/rng:.1f}% of range {rng:.1f}) exceeds 5% lock")
    if dmax > 0.20 * rng:
        fails.append(f"linear-vs-quadratic saturated head max {dmax:.3f} "
                     f"({100*dmax/rng:.1f}% of range) exceeds 20% guard")

    # (2) downstream seepage-face exit-point elevation: the highest point on the
    # downstream ground line whose (interior) pressure head is >= 0.
    def _exit_elev(pf):
        gxy = np.asarray(ground)
        xmid = 0.5 * (gxy[:, 0].min() + gxy[:, 0].max())
        yspan = gxy[:, 1].max() - gxy[:, 1].min()
        eps = 1e-3 * yspan
        best = -np.inf
        for i in range(len(gxy) - 1):
            for s in np.linspace(0, 1, 25):
                x = gxy[i, 0] + s * (gxy[i + 1, 0] - gxy[i, 0])
                yb = gxy[i, 1] + s * (gxy[i + 1, 1] - gxy[i, 1])
                if x < xmid:
                    continue
                p = pf(x, yb - eps)          # nudge just inside the domain
                if np.isfinite(p) and p >= -1e-6 and yb > best:
                    best = yb
        return best

    e3, e6 = _exit_elev(p3f), _exit_elev(p6f)
    if not (np.isfinite(e3) and np.isfinite(e6)):
        fails.append(f"could not locate a downstream exit point (tri3={e3}, tri6={e6})")
    else:
        de = abs(e3 - e6)
        # a one-element shift on the downstream face is the coarse-mesh floor
        if de > 2.0:
            fails.append(f"exit-point elevation disagrees: tri3={e3:.2f}, "
                         f"tri6={e6:.2f} (delta {de:.2f} > 2.0 ft)")
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


# ---------------------------------------------------------------------------
# Lock 4 — frames round trip + render (phase 3, plan §4)
# ---------------------------------------------------------------------------
def check_frames_roundtrip():
    import tempfile
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/seep/files/xslope_earth_dam1.xlsx")
    d = load_slope_data(path)
    for m in d["materials"]:
        m["Ss"], m["Sy"] = 1e-3, 0.2
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 18.0, "tri3")
    seep_data = _quiet(build_seep_data, mesh, d)

    n = len(seep_data["nodes"])
    fixed = seep_data["bc_values"][seep_data["bc_type"] == 1]
    h_uniform = np.full(n, float(np.mean(fixed)))
    tseep_data = dict(times=np.array([]), series={}, duration=4.0,
                      save_interval=2.0, save_times=[], stage_1=None,
                      stage_2=None, breakpoints=[])
    sol = _quiet(run_transient_seepage, seep_data, tseep_data, theta=1.0,
                 h_init=h_uniform, max_head_change_frac=0.15, verbose=False)

    for fr in sol["frames"]:
        phi = fr.get("phi")
        if phi is None or not np.isfinite(np.asarray(phi)).any():
            fails.append(f"frame t={fr['time']}: phi missing/NaN (solver product)")
        if fr.get("inflow") is None or fr.get("outflow") is None:
            fails.append(f"frame t={fr['time']}: inflow/outflow missing")

    with tempfile.TemporaryDirectory() as td:
        base = os.path.join(td, "rt")
        csv_path, meta_path = _quiet(export_transient_solution, seep_data, sol,
                                     base, input_file=path)
        if not (os.path.exists(csv_path) and os.path.exists(meta_path)):
            return ["export did not write both {base}_tseep.csv and _meta.json"]
        loaded = _quiet(import_transient_solution, seep_data, base)
        if len(loaded["frames"]) != len(sol["frames"]):
            fails.append("reloaded frame count mismatch")
        for a, b in zip(sol["frames"], loaded["frames"]):
            if np.max(np.abs(np.asarray(a["head"]) - b["head"])) > 1e-9:
                fails.append(f"head round-trip drift at t={a['time']}")
            if np.max(np.abs(np.asarray(a["u"]) - b["u"])) > 1e-9:
                fails.append(f"u round-trip drift at t={a['time']}")
            if np.max(np.abs(np.asarray(a["phi"]) - b["phi"])) > 1e-9:
                fails.append(f"phi round-trip drift at t={a['time']}")
            v = b.get("velocity")
            if v is None or not np.isfinite(v).any():
                fails.append(f"derived-on-load velocity missing at t={a['time']}")

        # one reloaded frame must render through the UNCHANGED plotter
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from xslope.plot_seep import plot_seep_solution
        fig = plt.figure(figsize=(8, 4))
        try:
            _quiet(plot_seep_solution, seep_data, loaded["frames"][-1], fig=fig,
                   variable="head", flowlines=True, vectors=True, phreatic=True)
        except Exception as e:
            fails.append(f"plot_seep_solution raised on a reloaded frame: {e!r}")
        finally:
            plt.close(fig)
    return fails


# ---------------------------------------------------------------------------
# Lock 5 — in-memory drawdown staging reproduces the classic FS (phase 3, §6)
# ---------------------------------------------------------------------------
def check_drawdown_staging():
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope.advanced import rapid_drawdown

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/lem/files/xslope_earth_dam_rapid.xlsx")
    if not os.path.exists(path):
        return [f"rapid-drawdown model not found: {path}"]
    d = load_slope_data(path)
    if "seep_u" not in d or "seep_u2" not in d:
        return ["classic rapid model did not load seep_u/seep_u2 from files"]
    u1 = np.asarray(d["seep_u"], dtype=float).copy()
    u2 = np.asarray(d["seep_u2"], dtype=float).copy()
    circles = d.get("circles") or []
    if not circles:
        return ["rapid model carries no circles to slice"]
    circle = circles[0]

    okS, (dfA, _s) = _quiet(generate_slices, d, circle=circle, debug=False)
    if not okS:
        return ["generate_slices (classic) failed"]
    okA, resA = _quiet(rapid_drawdown, dfA, "spencer", debug_level=0)
    if not okA:
        return [f"classic drawdown failed: {resA}"]
    FS_A = resA["FS"]

    # synthetic transient solution: stage_1 frame == steady u1, stage_2 == u2.
    d["tseep"] = {"stage_1": 0.0, "stage_2": 10.0}
    ts_sol = {"times": [0.0, 5.0, 10.0],
              "frames": [{"time": 0.0, "u": u1},
                         {"time": 5.0, "u": 0.5 * (u1 + u2)},
                         {"time": 10.0, "u": u2}]}
    d["seep_u"] = d["seep_u2"] = None
    _quiet(stage_transient_for_drawdown, d, ts_sol)
    if np.any(np.asarray(d["seep_u"]) != u1) or np.any(np.asarray(d["seep_u2"]) != u2):
        fails.append("staged seep_u/seep_u2 differ from the stage frames")
    okS2, (dfB, _s2) = _quiet(generate_slices, d, circle=circle, debug=False)
    if not okS2:
        return fails + ["generate_slices (staged) failed"]
    okB, resB = _quiet(rapid_drawdown, dfB, "spencer", debug_level=0)
    if not okB:
        return fails + [f"staged drawdown failed: {resB}"]
    if abs(FS_A - resB["FS"]) > 1e-9:
        fails.append(f"drawdown FS mismatch: classic {FS_A:.6f} vs staged {resB['FS']:.6f}")

    # single-time selection picks the requested frame
    select_transient_frame_u(d, ts_sol, time=5.0)
    if np.max(np.abs(np.asarray(d["seep_u"]) - 0.5 * (u1 + u2))) > 0:
        fails.append("select_transient_frame_u(time=5) did not select the mid frame")
    # a time with no saved frame is a clear error, not a silent nearest-pull
    try:
        transient_frame_index(ts_sol, 7.3)
        fails.append("transient_frame_index accepted a time with no saved frame")
    except ValueError:
        pass
    return fails


# ---------------------------------------------------------------------------
# Lock 6 — slow-drawdown consistency vs the classic two-steady FS (phase 5, §8.4)
# ---------------------------------------------------------------------------
def check_drawdown_consistency():
    """A slow (quasi-steady) transient reservoir drawdown must reproduce the
    classic two-independent-steady-solve rapid-drawdown FS; a fast drawdown must
    not. See the module docstring (lock 6) for the tolerance derivation."""
    import copy
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.slice import generate_slices
    from xslope.advanced import rapid_drawdown

    fails = []
    path = os.path.join(_REPO_ROOT, "docs/lem/files/xslope_earth_dam_rapid.xlsx")
    if not os.path.exists(path):
        return [f"rapid-drawdown model not found: {path}"]
    d = load_slope_data(path)
    circles = d.get("circles") or []
    if not circles:
        return ["rapid model carries no circles to slice"]
    H1, H2 = 302.0, 250.0          # high / low pool levels (the model's two stages)

    # Storage params (the steady file leaves Ss/Sy blank); inject a modest storage
    # so the transient has a finite relaxation time.
    for m in d["materials"]:
        m["Ss"], m["Sy"] = 1e-4, 0.2

    # Coarse LINEAR mesh (~60 nodes) — the transient exit-face active set is
    # resolved per corner, so tri3 keeps the receding seepage face exact and fast.
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    width = max(xs) - min(xs)
    mesh = _quiet(build_mesh_from_polygons, polys, width / 16.0, "tri3")

    # Bind the upstream reservoir head to a series "res" (its value becomes a string
    # -> build_seep_data records a head_series_binding on the submerged face nodes;
    # run_transient_seepage holds each node at h(t) while submerged and turns it into
    # an exit face above the water line). The downstream exit face is unchanged.
    d_ts = copy.deepcopy(d)
    for sh in d_ts["seepage_bc"]["specified_heads"]:
        sh["head"] = "res"
    seep_data = _quiet(build_seep_data, mesh, d_ts)
    if not seep_data.get("head_series_bindings"):
        return ["reservoir head binding was not created (check seepage_bc mapping)"]

    def _run(series_ta, series_va, duration, save_times, dt_max, frac=0.1):
        tseep = dict(times=np.asarray(series_ta, dtype=float),
                     series={"res": (np.asarray(series_ta, dtype=float),
                                     np.asarray(series_va, dtype=float))},
                     duration=duration, save_interval=duration, save_times=save_times,
                     stage_1=0.0, stage_2=duration, breakpoints=list(series_ta))
        return _quiet(run_transient_seepage, seep_data, tseep, theta=1.0,
                      dt_max=dt_max, max_head_change_frac=frac, verbose=False)

    def _frame_u_at(sol, t):
        return np.asarray(sol["frames"][transient_frame_index(sol, t)]["u"], dtype=float)

    def _steady_u(level):
        # The transient IC IS a steady solve at the t=0 BC configuration, so frame 0
        # of a constant-series run is exactly the independent steady solution at that
        # pool level (same submerged-Dirichlet + seepage-face treatment as the
        # drawdown run, so the comparison isolates the time-marching residual).
        sol = _run([0.0], [level], 1.0, [], dt_max=1.0)
        return np.asarray(sol["frames"][0]["u"], dtype=float)

    def _fs(u1, u2):
        dd = copy.deepcopy(d)
        dd["mesh"], dd["seep_u"], dd["seep_u2"] = mesh, u1, u2
        okS, (df, _s) = _quiet(generate_slices, dd, circle=circles[0], debug=False)
        if not okS:
            raise RuntimeError("generate_slices failed")
        okA, res = _quiet(rapid_drawdown, df, "spencer", debug_level=0)
        if not okA:
            raise RuntimeError(f"rapid_drawdown failed: {res}")
        return res["FS"]

    # Classic path: two independent steady solves (high pool, low pool).
    u_steady_hi = _steady_u(H1)
    u_steady_lo = _steady_u(H2)
    FS_classic = _fs(u_steady_hi, u_steady_lo)

    # Slow drawdown: hold high, ramp to low over t_dd, hold low. Save the low-pool
    # state at three growing equilibration times so the FS gap vs equilibration can
    # be shown to plateau from ONE solve.
    t_pre, t_dd, holds = 0.5, 5.0, [2.0, 4.0, 8.0]
    t_end_dd = t_pre + t_dd
    save_ts = [t_end_dd + h for h in holds]
    dur = save_ts[-1]
    sol = _run([0.0, t_pre, t_end_dd, dur], [H1, H1, H2, H2], dur, save_ts, dt_max=0.25)
    if not sol["converged"]:
        fails.append("slow-drawdown transient reported non-convergence")

    # The stage-1 (t=0) transient field must equal the independent high-pool steady
    # solve to machine zero — that is what makes the FS gap purely the stage-2 residual.
    u_stage1 = _frame_u_at(sol, 0.0)
    d_ic = float(np.max(np.abs(u_stage1 - u_steady_hi)))
    if d_ic > 1e-9:
        fails.append(f"stage-1 transient field differs from high-pool steady by "
                     f"{d_ic:.2e} (expected machine zero — shared IC)")

    dfs = [abs(_fs(u_stage1, _frame_u_at(sol, st)) - FS_classic) for st in save_ts]

    # (a) convergence: the FS gap shrinks as equilibration grows (asymptotic approach).
    if not (dfs[-1] < dfs[0] and dfs[-1] < 0.6 * dfs[0]):
        fails.append(f"slow drawdown did not converge toward the steady FS as "
                     f"equilibration grew: gaps {[f'{x:.2e}' for x in dfs]}")
    # (b) plateau agreement: locked an order of magnitude above the demonstrated ~5e-5.
    TOL = 5e-4
    if dfs[-1] > TOL:
        fails.append(f"slow-drawdown FS gap {dfs[-1]:.2e} at the longest "
                     f"equilibration exceeds the {TOL:.0e} lock")

    # (c) discriminator: a FAST drawdown leaves stage-2 far from steady and misses.
    dur_f = t_pre + 0.4 + 0.1
    solf = _run([0.0, t_pre, t_pre + 0.4, dur_f], [H1, H1, H2, H2], dur_f, [],
                dt_max=0.05)
    dfs_fast = abs(_fs(_frame_u_at(solf, 0.0), _frame_u_at(solf, dur_f)) - FS_classic)
    if not (dfs_fast > 5e-3 and dfs_fast > 20.0 * dfs[-1]):
        fails.append(f"fast drawdown failed to discriminate: gap {dfs_fast:.2e} is "
                     f"not >> the slow-drawdown gap {dfs[-1]:.2e} (test may be vacuous)")
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
        ("steady limit (unconfined, tri3)", check_steady_limit_unconfined),
        ("steady limit (unconfined, tri6)", check_steady_limit_unconfined_tri6),
        ("drawdown tri3-vs-tri6 consistency", check_drawdown_linear_vs_quadratic),
        ("frames round trip + render", check_frames_roundtrip),
        ("drawdown staging (§6)", check_drawdown_staging),
        ("drawdown consistency (§8.4)", check_drawdown_consistency),
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
