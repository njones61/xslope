"""The two SSRM drivers must agree on one locked case.

The strength-reduction bisection can be driven two ways (see SPIKE.md): the
viscoplastic initial-stiffness iteration, which is the default and the definition
of every locked and published factor of safety, and a Newton-Raphson iteration on
a consistent Mohr-Coulomb tangent, selected by ``fem_solver='newton'``. They are
different numerical routes to the same physical question — at what strength
reduction does the slope stop standing — so they have to reach the same answer.

This pins that on Griffiths & Lane Example 1 (docs/fem/files/xslope_griffiths1.xlsx)
at a deliberately coarse tri6 mesh, which the verification page locks at FS = 1.39.
Three things are asserted:

  * both drivers land within the bisection tolerance of the locked value;
  * they land within the bisection tolerance of EACH OTHER;
  * the Newton run leaves no trial undecided — every trial the bisection visits
    comes back converged or failed, never 'inconclusive'.

The last one is the point of the whole exercise. The viscoplastic path can stop a
trial at its iteration ceiling with the residual still falling and no verdict to
give; the Newton path decides a trial by whether the full gravity load is reached
in equilibrium, so it cannot.

The mesh is coarse on purpose — this check exists to catch a regression in either
driver, not to re-measure the benchmark, and it has to stay cheap enough to keep
in the suite.

Four more checks were added by the 2026-08-31 bug hunt (SPIKE.md, "Bug hunt").
That hunt tested each of the three suspects the spike had left open for the
benchmarks where the two drivers disagree, and cleared all three; what it found
instead was on the viscoplastic side. These lock the cleared properties, so a
future change to the Newton path cannot quietly reintroduce any of them:

  * the return map lands exactly on the yield surface on every branch, moves no
    mean stress with psi = 0, and never returns a state outside the cone;
  * a converged Newton trial certifies itself — its stress field satisfies the
    Mohr-Coulomb function to machine precision, measured in the INVARIANT form
    the viscoplastic path uses rather than the principal-stress form the return
    map is written on;
  * the verdict does not depend on the loading path (one gravity increment or
    sixteen give the same answer);
  * the verdict does not depend on the step-control knobs (the displacement
    runaway gate, the no-progress window, the per-increment iteration cap and the
    load-step floor can all be loosened without changing it).

Run directly:  PYTHONPATH=. python3 test/nr_ssrm_check.py
"""

import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

from xslope import fem as _fem
from xslope.fem import (build_fem_data, mc_return_map, solve_fem, solve_ssrm,
                        _NR_APEX, _NR_ELASTIC)
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons, get_material_polygons

MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files'
         / 'xslope_griffiths1.xlsx')
ELEMENT_TYPE = 'tri6'
TARGET_SIZE = 6.0
LOCKED_FS = 1.39         # docs/verification/ssrm.md, the coarse tri6 row
F_MIN, F_MAX = 1.30, 1.50
TOLERANCE = 0.05         # the bisection cell width, and the agreement tolerance
MAX_ITER = 4000
# Two trials on the same coarse mesh, one either side of its failure boundary, for
# the checks that ask what a verdict depends on rather than what it is.
F_STANDS, F_FAILS = 1.35, 1.45


def _fem_data():
    slope_data = load_slope_data(str(MODEL))
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=TARGET_SIZE,
                                    element_type=ELEMENT_TYPE)
    return build_fem_data(slope_data, mesh)


def run():
    failures = []
    fem_data = _fem_data()
    results = {}
    for solver in ('viscoplastic', 'newton'):
        results[solver] = solve_ssrm(
            fem_data, F_min=F_MIN, F_max=F_MAX, tolerance=TOLERANCE,
            max_iterations=MAX_ITER, fem_solver=solver,
            capture_failure_state=False)

    for solver, res in results.items():
        fs = res.get('FS')
        if fs is None:
            failures.append(f"{solver}: no factor of safety returned "
                            f"({res.get('message', 'no message')})")
            continue
        if abs(fs - LOCKED_FS) > TOLERANCE:
            failures.append(
                f"{solver}: FS = {fs:.4f} is {abs(fs - LOCKED_FS):.4f} from the "
                f"locked {LOCKED_FS:.2f}, outside the bisection tolerance "
                f"{TOLERANCE:g}")

    fs_vp = results['viscoplastic'].get('FS')
    fs_nr = results['newton'].get('FS')
    if fs_vp is not None and fs_nr is not None and abs(fs_vp - fs_nr) > TOLERANCE:
        failures.append(
            f"the two drivers disagree: viscoplastic {fs_vp:.4f} vs Newton "
            f"{fs_nr:.4f}, a gap of {abs(fs_vp - fs_nr):.4f} against the "
            f"bisection tolerance {TOLERANCE:g}")

    undecided = [t['F'] for t in results['newton'].get('trials', [])
                 if t.get('exit_reason') == 'inconclusive']
    if undecided:
        failures.append(
            "the Newton driver left trial(s) undecided at F = "
            + ", ".join(f"{f:.4f}" for f in undecided)
            + " — its verdict is meant to be binary")

    failures += check_return_map()
    failures += check_verdict_evidence(fem_data)
    failures += check_load_path_invariance(fem_data)
    failures += check_step_control_not_decisive(fem_data)

    return failures


def check_return_map():
    """The return map, audited against invariants it cannot satisfy by accident.

    Random trial stresses at four friction angles, checked for the properties the
    algorithm is supposed to have: a yielding point comes back ON the surface (not
    inside it, which would be extra strength, and not outside, which would be
    none); an elastic point comes back untouched; and with psi = 0 the correction
    is purely deviatoric, so no branch except the apex may move the mean stress.
    A branch that returned a point to the wrong place would show up here as a
    residual yield value, a moved mean stress, or a scrambled principal order.
    """
    fails = []
    rng = np.random.default_rng(7)
    n = 40000
    for phi_deg, c in ((0.0, 10.0), (20.0, 5.0), (35.0, 2.0), (45.0, 1.0)):
        snph, csph = np.sin(np.radians(phi_deg)), np.cos(np.radians(phi_deg))
        scale = 30.0
        sig = rng.normal(0.0, scale, size=(n, 4))
        sig[:, 2] *= 0.5
        out, branch = mc_return_map(sig, np.full(n, c), np.full(n, snph),
                                    np.full(n, csph), np.full(n, 4000.0))

        def principals(s):
            ctr = 0.5 * (s[:, 0] + s[:, 1])
            half = 0.5 * (s[:, 0] - s[:, 1])
            r = np.sqrt(half * half + s[:, 2] ** 2)
            return -np.sort(-np.stack([ctr + r, ctr - r, s[:, 3]], axis=1), axis=1)

        pt, pn = principals(sig), principals(out)
        a, bc = 0.5 * (1.0 + snph), 0.5 * (1.0 - snph)
        f_new = a * pn[:, 0] - bc * pn[:, 2] - c * csph
        tol = 1e-9 * max(c * csph, scale)
        tag = f"phi={phi_deg:g}"

        if f_new.max() > tol:
            fails.append(f"{tag}: a returned state is OUTSIDE the yield surface "
                         f"(max f = {f_new.max():.3e}, tolerance {tol:.1e})")
        on_surface = (branch != _NR_ELASTIC) & (branch != _NR_APEX)
        if on_surface.any() and np.abs(f_new[on_surface]).max() > tol:
            fails.append(f"{tag}: a yielded state did not land ON the surface "
                         f"(max |f| = {np.abs(f_new[on_surface]).max():.3e})")
        elastic = branch == _NR_ELASTIC
        if elastic.any() and np.abs(out[elastic] - sig[elastic]).max() > 0.0:
            fails.append(f"{tag}: an elastic point was modified by the return map")
        # psi = 0: every flow direction is deviatoric, so only the apex may move p.
        dp = ((out[:, 0] + out[:, 1] + out[:, 3])
              - (sig[:, 0] + sig[:, 1] + sig[:, 3])) / 3.0
        moved = (branch != _NR_APEX) & (np.abs(dp) > 1e-9 * scale)
        if moved.any():
            fails.append(f"{tag}: {int(moved.sum())} non-apex return(s) changed the "
                         f"mean stress, which psi = 0 flow cannot do "
                         f"(max |dp| = {np.abs(dp[moved]).max():.3e})")
        if int(((pn[:, 0] < pn[:, 1] - tol) | (pn[:, 1] < pn[:, 2] - tol)).sum()):
            fails.append(f"{tag}: a returned state has its principal stresses out "
                         f"of order")
        radius_grew = ((0.5 * (pn[:, 0] - pn[:, 2]) > 0.5 * (pt[:, 0] - pt[:, 2]) + tol)
                       & (branch != _NR_ELASTIC))
        if radius_grew.any():
            fails.append(f"{tag}: {int(radius_grew.sum())} return(s) GREW the "
                         f"deviatoric radius instead of shrinking it")
    return fails


def _newton(fem_data, F, **over):
    """One Newton trial, optionally with step-control globals overridden."""
    saved = {k: getattr(_fem, k) for k in over}
    for k, v in over.items():
        setattr(_fem, k, v)
    try:
        return solve_fem(fem_data, F=F, fem_solver='newton', max_disp_factor=None,
                         force_tol=1e-3)
    finally:
        for k, v in saved.items():
            setattr(_fem, k, v)


def check_verdict_evidence(fem_data):
    """A converged Newton trial has to carry the proof of its own verdict.

    'This slope stands at F' means a stress field that carries full gravity in
    equilibrium with nothing outside the yield surface. Both halves are measured
    and reported, so the check is on the numbers rather than on the label: the
    Dawson out-of-balance under the trial tolerance, and the largest Mohr-Coulomb
    violation at machine precision. The yield measure is computed in the invariant
    form the viscoplastic path uses, so it is an independent reading of the return
    map's own algebra.
    """
    fails = []
    sol = _newton(fem_data, F_STANDS)
    if not sol['converged']:
        fails.append(f"F = {F_STANDS} was expected to stand on this mesh and the "
                     f"Newton driver failed it ({sol.get('exit_reason')})")
        return fails
    # `unbalanced_force_ratio` is the like-for-like key: it is the Dawson
    # out-of-balance on BOTH drivers. (`residual` is the relative displacement
    # change on both, which is a different quantity and not the force gate.)
    if sol['unbalanced_force_ratio'] >= 1e-3:
        fails.append(f"the converged trial at F = {F_STANDS} reports an "
                     f"out-of-balance of {sol['unbalanced_force_ratio']:.3e}, at or "
                     f"above the 1e-3 force tolerance it is meant to have passed")
    viol = sol.get('nr_max_yield_violation')
    if viol is None:
        fails.append("the Newton result carries no 'nr_max_yield_violation': a "
                     "converged trial no longer certifies its own stress field")
    elif viol > 1e-6:
        fails.append(f"the converged trial at F = {F_STANDS} leaves Gauss points "
                     f"{viol:.3e} outside the yield surface (fraction of the local "
                     f"strength); the field is not statically admissible, so the "
                     f"verdict is not evidence")
    return fails


def check_load_path_invariance(fem_data):
    """The verdict must not depend on how the gravity load is applied.

    Mohr-Coulomb with non-associated flow is path-dependent, so this is a real
    question rather than a formality: if the answer moved with the increment size,
    a factor of safety would be reporting the solver's step schedule. Measured
    across one increment and sixteen at F either side of the boundary.
    """
    fails = []
    for F in (F_STANDS, F_FAILS):
        one = _newton(fem_data, F)
        many = _newton(fem_data, F, _NR_INIT_STEP=1.0 / 16, _NR_GROW=1.0)
        if bool(one['converged']) != bool(many['converged']):
            fails.append(
                f"F = {F}: the verdict depends on the loading path — one gravity "
                f"increment says {one['verdict']}, sixteen say {many['verdict']}")
    return fails


def check_step_control_not_decisive(fem_data):
    """Loosening every step-control knob must not change a verdict.

    These are solver controls with no physical meaning, and the spike's own worst
    finding was a configuration of them that moved a factor of safety by two
    thirds. So they are held to the standard that they may change how hard the
    solve works and never what it concludes: the no-progress window widened, the
    per-increment iteration cap raised, and the load-step floor dropped
    sixteenfold, all at once. The fourth knob that used to sit here, a displacement
    runaway gate, is gone from the driver — it decided a verdict, which is why it
    is gone (see SPIKE.md, "Bug hunt").
    """
    fails = []
    loose = dict(_NR_NO_PROGRESS=200, _NR_MAX_ITER=400, _NR_MIN_STEP=1.0 / 1024)
    for F in (F_STANDS, F_FAILS):
        shipped = _newton(fem_data, F)
        relaxed = _newton(fem_data, F, **loose)
        if bool(shipped['converged']) != bool(relaxed['converged']):
            fails.append(
                f"F = {F}: the verdict is decided by the step control — as shipped "
                f"it is {shipped['verdict']}, with every knob loosened it is "
                f"{relaxed['verdict']}")
    return fails


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nBoth SSRM drivers reproduce the locked factor of safety, agree with "
          "each other, and the Newton run decided every trial. The return map "
          "lands on the yield surface on every branch, a converged trial "
          "certifies its own stress field, and neither the loading path nor the "
          "step control changes a verdict.")


if __name__ == '__main__':
    main()
