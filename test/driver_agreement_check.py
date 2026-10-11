"""The three strength-reduction drivers must agree on problems with a known answer.

A strength-reduction trial can be run by three drivers (``fem_solver``):

  * ``'viscoplastic'`` — the plain Griffiths & Lane viscoplastic loop, with no
    corrector. It is the definition of the locked factors of safety.
  * ``'auto'`` — the default: the same loop, with a bounded Newton corrector
    offered the loop's state at a short checkpoint ladder and at every stopping
    rule. A corrector that reaches equilibrium inside the force, yield and
    displacement gates ends the trial as standing.
  * ``'newton'`` — the cold-start Newton-Raphson driver, with its own load walk.

They are three different routes to the same statement (the slope stands at F, or
it does not), so on a problem whose answer is known they must land on it, and on
each other. A defect in one driver's shared machinery — the corrector's
factorization cache once reused a wrong gather index and silently mis-solved
(LOOSE_ENDS item 2) — shows up here as a driver that disagrees.

Three small problems, each solved under every driver on the reference (NumPy)
kernel, on the same mesh and settings:

  1. **A homogeneous soil slope.** Griffiths & Lane (1999) Example 1 (2:1,
     c / gamma H = 0.05, phi = 20 deg) on a coarse tri6 mesh, against Bishop &
     Morgenstern's chart value 1.380.
  2. **A block on a plane.** A cohesionless slab on a jointed plane at 20 deg
     with phi_j = 30 deg: FS = tan phi_j / tan beta = 1.5865.
  3. **Two blocks on a plane.** The same slab cut in two by a vertical joint (a T
     on the plane and a crack to the surface): the same closed form.

What each problem must show:

  (a) every driver's bracket lands on the known answer, within the model's own
      allowance (the mesh and the interface's penalty stiffness, stated per
      problem) plus half a bisection step;
  (b) the three drivers' answers are within one bisection step of each other;
  (c) at nine tenths of the highest factor all three stood at, the three
      drivers' displacement fields agree to ``FIELD_TOL`` (1%) of the largest
      displacement; the reasoning for that number, and the measurement of the
      loop's own error that backs it, are at ``FIELD_TOL``.

A driver that raises on a model is a FAIL that names the driver and the error;
nothing is skipped silently.

Run directly:  PYTHONPATH=. python3 test/driver_agreement_check.py
"""
import contextlib
import io
import math
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402

FAILURES = []

#: The drivers, in the order they are reported.
DRIVERS = ('viscoplastic', 'auto', 'newton')

#: Bisection tolerance handed to solve_ssrm: the search stops once its bracket
#: is narrower than this. "One bisection step" in check (b) is the width of the
#: final bracket the drivers actually reached, which is at most this.
STEP = 0.01

#: The viscoplastic loop's stopping test (solve_fem's defaults): a trial has
#: converged when a sweep changes the field by less than TOLERANCE of max|u| and
#: every node's out-of-balance force is under FORCE_TOL of its own weight.
TOLERANCE = 1e-3
FORCE_TOL = 1e-3

#: How far two drivers' converged displacement fields may differ: the largest
#: nodal difference as a fraction of the largest nodal displacement.
#:
#: Why 1%. A converged viscoplastic state still moves by up to TOLERANCE of
#: max|u| per sweep, and the sweep is a contraction, so the distance left to the
#: fixed point is that last change times 1 / (1 - rho) for a contraction rate
#: rho. Allowing rho up to 0.9 puts each driver within 10 x TOLERANCE = 1% of the
#: equilibrium it is converging to. The Newton paths stop on the same FORCE_TOL
#: and are far tighter. The check measures the loop's own error on each problem
#: (the same solve re-run with both tolerances ten times tighter) and prints it,
#: and requires it to be under a tenth of this allowance: a model where the loop
#: itself is not converged that well cannot referee the drivers.
FIELD_TOL = 10.0 * TOLERANCE

#: The fields are compared at this fraction of the highest factor all three
#: drivers stood at. At the standing edge itself the slope is at its limit and
#: the plain loop's own convergence error grows without bound (on the soil slope
#: at F = 1.40 it is 6% of max|u| against a run ten times tighter), so the edge
#: cannot referee anything; nine tenths of it is well inside the stable range and
#: still carries plastic flow.
FIELD_AT = 0.9


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


class _reference_kernel:
    """Force every solve_fem call inside solve_ssrm onto the NumPy reference
    kernel (solve_ssrm has no fast_kernel argument; run_tests does the same)."""

    def __enter__(self):
        self._orig = fem.solve_fem
        orig = self._orig

        def _wrap(*a, **k):
            k['fast_kernel'] = False
            return orig(*a, **k)
        fem.solve_fem = _wrap
        return self

    def __exit__(self, *exc):
        fem.solve_fem = self._orig
        return False


# --------------------------------------------------------------------------
# The models
# --------------------------------------------------------------------------

def _soil_slope():
    """Griffiths & Lane Example 1 on the coarse tri6 mesh fem_guards uses."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    xlsx = os.path.join(_ROOT, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(xlsx)
        mesh = build_mesh_from_polygons(get_material_polygons(sd),
                                        target_size=8.0, element_type='tri6')
        return fem.build_fem_data(sd, mesh)


def _block():
    """joint_element_check's row 2: the slab on a plane at 20 deg, phi_j = 30."""
    import joint_element_check as jec
    with contextlib.redirect_stdout(io.StringIO()):
        _d, fd, _g = jec.slab_model(**jec.ROW2)
    return fd


def _two_blocks():
    """joint_junction_check's row i: the same slab cut by a vertical joint."""
    import joint_junction_check as jjc
    with contextlib.redirect_stdout(io.StringIO()):
        _d, _mesh, fd, _lines = jjc._slab(cut=True)
    return fd


_TAN = math.tan(math.radians(30.0)) / math.tan(math.radians(20.0))

#: name, builder, known answer, where it comes from, allowance (relative), and
#: the search settings. The jointed rows take joint_element_check's budget
#: (3,000 sweeps, ceiling 6,000) and its 3% allowance for the penalty
#: interface; the soil slope takes the tri6 row's budget from ssrm.md and a 4%
#: allowance for the coarse mesh (the tri6 row at target 6 locks 1.39 against the
#: chart's 1.380).
PROBLEMS = (
    dict(name='soil slope (Griffiths & Lane Ex. 1, coarse tri6)',
         build=_soil_slope, known=1.380,
         source="Bishop & Morgenstern chart", allow=0.04,
         ssrm=dict(F_min=1.3, F_max=1.5, max_iterations=4000)),
    dict(name='block on a plane (beta 20, phi_j 30)',
         build=_block, known=_TAN, source="tan phi_j / tan beta", allow=0.03,
         ssrm=dict(F_min=1.4, F_max=1.75, max_iterations=3000,
                   max_iterations_ceiling=6000)),
    dict(name='two blocks on a plane (vertical cut)',
         build=_two_blocks, known=_TAN, source="tan phi_j / tan beta",
         allow=0.03,
         ssrm=dict(F_min=1.4, F_max=1.75, max_iterations=3000,
                   max_iterations_ceiling=6000)),
)


# --------------------------------------------------------------------------
# One problem
# --------------------------------------------------------------------------

def _run_driver(fd, driver, ssrm_kw):
    t0 = time.time()
    try:
        with _reference_kernel(), contextlib.redirect_stdout(io.StringIO()):
            res = fem.solve_ssrm(fd, tolerance=STEP, fem_solver=driver,
                                 capture_failure_state=False, **ssrm_kw)
    except Exception as exc:                        # reported, never skipped
        return None, f"{type(exc).__name__}: {exc}", time.time() - t0
    return res, None, time.time() - t0


def _field(fd, driver, F, ssrm_kw, tight=False):
    """The displacement field one driver reaches at F, with the search's
    per-trial budget; None and the reason where it does not stand. ``tight``
    re-runs it with both stopping tolerances ten times tighter and a budget to
    match, which is the yardstick for the loop's own convergence error."""
    kw = dict(max_iterations=ssrm_kw['max_iterations'])
    if 'max_iterations_ceiling' in ssrm_kw:
        kw['max_iterations_ceiling'] = ssrm_kw['max_iterations_ceiling']
    if tight:
        kw = dict(max_iterations=60000, max_iterations_ceiling=60000,
                  tolerance=0.1 * TOLERANCE, force_tol=0.1 * FORCE_TOL)
    else:
        kw.update(tolerance=TOLERANCE, force_tol=FORCE_TOL)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            sol = fem.solve_fem(fd, F=F, fast_kernel=False, fem_solver=driver,
                                max_disp_factor=None, **kw)
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
    if not sol.get('converged'):
        return None, f"did not stand ({sol.get('exit_reason')})"
    return np.asarray(sol['displacements'], dtype=float), None


def run_problem(p):
    print(f"\n{p['name']}: known FS = {p['known']:.4f} ({p['source']})")
    fd = p['build']()
    runs = {}
    for drv in DRIVERS:
        res, err, wall = _run_driver(fd, drv, p['ssrm'])
        if err is not None:
            check(f"{p['name']}: the '{drv}' driver runs this model", False,
                  f"it raised {err}")
            continue
        trials = res.get('trials') or []
        iters = sum(int(t.get('iterations') or 0) for t in trials)
        n_corr = sum(1 for t in trials if t.get('corrector'))
        lo, hi = res['final_interval']
        runs[drv] = dict(FS=float(res['FS']), lo=float(lo), hi=float(hi),
                         lower=bool(res.get('fs_is_lower_bound')),
                         iters=iters, n_trials=len(trials), n_corr=n_corr,
                         wall=wall)
        print(f"    {drv:<12s} bracket [{lo:.4f}, {hi:.4f}]  FS = {res['FS']:.4f}"
              f"{' (lower bound)' if res.get('fs_is_lower_bound') else ''}  "
              f"{len(trials)} trials, {iters:,} iterations"
              f"{f', {n_corr} certified by the corrector' if n_corr else ''}"
              f", {wall:.1f} s")

    # (a) each driver on the known answer
    for drv, r in runs.items():
        err = (r['FS'] - p['known']) / p['known']
        ok = (not r['lower']) and abs(err) <= p['allow'] + 0.5 * STEP / p['known']
        check(f"{p['name']}: '{drv}' lands on the known answer", ok,
              f"FS {r['FS']:.4f} vs {p['known']:.4f} ({100 * err:+.2f}%, "
              f"allowance {100 * p['allow']:.0f}% + half a step)")

    # (b) the drivers against each other
    names = list(runs)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = runs[names[i]], runs[names[j]]
            d = abs(a['FS'] - b['FS'])
            step = max(a['hi'] - a['lo'], b['hi'] - b['lo'])
            same = (abs(a['lo'] - b['lo']) < 1e-12 and abs(a['hi'] - b['hi']) < 1e-12)
            check(f"{p['name']}: '{names[i]}' and '{names[j]}' within one "
                  f"bisection step", d <= step + 1e-12,
                  "the same bracket" if same else
                  f"answers differ by {d:.4f}, {d / step:.1f} steps of "
                  f"{step:.4f}")

    # (c) the converged fields, below the highest factor every driver stood at
    if len(runs) < 2:
        return
    F_c = FIELD_AT * min(r['lo'] for r in runs.values())
    fields = {}
    for drv in runs:
        u, why = _field(fd, drv, F_c, p['ssrm'])
        if u is None:
            check(f"{p['name']}: '{drv}' stands at F = {F_c:.4f} for the field "
                  f"comparison", False, why)
        else:
            fields[drv] = u
    if len(fields) < 2:
        return
    uv = {d: np.column_stack(fem._extract_nodal_uv(u, fd))
          for d, u in fields.items()}
    ref = 'viscoplastic' if 'viscoplastic' in uv else next(iter(uv))
    scale = float(np.max(np.linalg.norm(uv[ref], axis=1)))
    names = list(uv)
    # The yardstick: the plain loop against itself at ten times tighter.
    u_t, why = _field(fd, 'viscoplastic', F_c, p['ssrm'], tight=True)
    if u_t is None or 'viscoplastic' not in uv:
        check(f"{p['name']}: the viscoplastic loop's own error at F = "
              f"{F_c:.4f} can be measured", False,
              why or "the plain loop has no field here")
    else:
        uvt = np.column_stack(fem._extract_nodal_uv(u_t, fd))
        own = float(np.max(np.linalg.norm(uv['viscoplastic'] - uvt, axis=1)))
        own /= max(scale, 1e-300)
        check(f"{p['name']}: the viscoplastic loop's own error at F = "
              f"{F_c:.4f} is under a tenth of the field allowance",
              own <= 0.1 * FIELD_TOL,
              f"{100 * own:.3f}% of max|u| against a run ten times tighter")
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            diff = float(np.max(np.linalg.norm(uv[names[i]] - uv[names[j]], axis=1)))
            rel = diff / max(scale, 1e-300)
            check(f"{p['name']}: '{names[i]}' and '{names[j]}' fields agree "
                  f"at F = {F_c:.4f}", rel <= FIELD_TOL,
                  f"largest nodal difference {100 * rel:.3f}% of max|u| "
                  f"({scale:.4g}); allowance {100 * FIELD_TOL:.0f}%")
    for drv, r in runs.items():
        print(f"    summary  {drv:<12s} [{r['lo']:.4f}, {r['hi']:.4f}]  "
              f"{r['iters']:,} iterations  {r['wall']:.1f} s")


def main():
    print("=" * 72)
    print("Driver agreement: viscoplastic, auto (corrector) and newton")
    print("=" * 72)
    t0 = time.time()
    for p in PROBLEMS:
        run_problem(p)
    print("\n" + "=" * 72)
    print(f"({time.time() - t0:.0f} s)")
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All three drivers agree on every problem.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


if __name__ == '__main__':
    sys.exit(main())
