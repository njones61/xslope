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
  * the verdict does not depend on the step-control knobs (the no-progress window,
    the per-increment iteration cap and the load-step floor can all be loosened
    without changing it), and the displacement gate that used to sit among them is
    not defined on the module at all.

Two more came from the 2026-08-31 corrections, and both hold behavior the first
four could not see:

  * a converged trial is admissible in DISPLACEMENT as well as in force — on the
    quad9 mesh, F = 1.396875 stands and F = 1.400 reaches machine-precision statics
    a sixth of the model height away and is refused for it, with the signal that
    says which of the two it failed;
  * the XSLOPE_FEM_SOLVER environment override announces itself, so a stale shell
    variable cannot recompute a session's factors of safety on the non-default
    driver in silence.

And one more for the monotonic strength-reduction ramp (SPIKE.md, "RAMP"), the
third SSRM route: it lands on the same limit the Newton bisection does, which is a
real question because the two carry different plastic histories to the same
strength, it reports the MIDPOINT of its final interval — the convention every
locked and published factor of safety here is defined on — and it never evaluates a
strength more than its initial increment above the highest one it has carried,
which is the whole of its cost advantage.

Run directly:  PYTHONPATH=. python3 test/nr_ssrm_check.py
"""

import os
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

# The same model on a FINER quad9 mesh, which is where the displacement bound
# actually decides something (see check_displacement_bound). One trial either side
# of the bound; the pair brackets the Newton factor of safety without paying for a
# whole bisection.
QUAD9_SIZE = 3.5
Q9_ADMISSIBLE, Q9_INADMISSIBLE = 1.396875, 1.400


def _mesh(element_type=ELEMENT_TYPE, target_size=TARGET_SIZE):
    slope_data = load_slope_data(str(MODEL))
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=target_size,
                                    element_type=element_type)
    return slope_data, mesh


def _fem_data(element_type=ELEMENT_TYPE, target_size=TARGET_SIZE):
    slope_data, mesh = _mesh(element_type, target_size)
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
    failures += check_displacement_bound()
    failures += check_env_override_announces_itself()
    failures += check_ramp(fem_data, results['newton'].get('FS'))

    return failures


# The ramp's first strength increment, written out rather than imported: the
# never-past-the-limit property is bounded BY that increment, so reading it from
# the module would make the assertion follow the code instead of holding it.
RAMP_DF_INIT = 0.05


def check_ramp(fem_data, fs_bisection):
    """The monotonic ramp lands on the bisection's answer and never overshoots it.

    Two properties, on the same coarse mesh the bisection ran above (see SPIKE.md,
    "RAMP"):

      * the answer. The ramp reduces strength along ONE warm-started history
        instead of solving each trial from zero, and non-associated Mohr-Coulomb is
        path-dependent, so this is a real question rather than a formality — the
        two routes have to reach the same limit.
      * it never goes past failure. The bisection's worst trials are the ones far
        above the limit, where the Newton driver costs 17x to 47x the viscoplastic
        driver's constitutive work to prove a load unreachable; the ramp stops one
        increment past the highest strength it has carried, by construction. That
        is asserted as an exact bound: no strength evaluated may exceed the last
        one carried by more than the ramp's initial increment.
      * it reports the MIDPOINT of its final interval. Every locked and published
        factor of safety in this repository is a bisection midpoint, and the ramp
        used to report the last strength it CARRIED, floored to 0.01 — the
        interval's lower edge, rounded down again. Across the eight spike
        benchmarks that read 0.0031 to 0.0119 low against the bisection, none of it
        a difference in what the two routes found. Asserted against the interval
        the same result reports, so the two can never drift apart.
    """
    fails = []
    res = solve_ssrm(fem_data, F_min=F_MIN, F_max=F_MAX, tolerance=TOLERANCE,
                     max_iterations=MAX_ITER, fem_solver='newton',
                     ssrm_driver='ramp', capture_failure_state=False)
    fs = res.get('FS')
    if fs is None:
        return [f"the ramp returned no factor of safety "
                f"({res.get('error', 'no message')})"]
    if abs(fs - LOCKED_FS) > TOLERANCE:
        fails.append(f"ramp: FS = {fs:.4f} is {abs(fs - LOCKED_FS):.4f} from the "
                     f"locked {LOCKED_FS:.2f}, outside the tolerance {TOLERANCE:g}")
    if fs_bisection is not None and abs(fs - fs_bisection) > TOLERANCE:
        fails.append(f"the two Newton SSRM drivers disagree: bisection "
                     f"{fs_bisection:.4f} vs ramp {fs:.4f}, a gap of "
                     f"{abs(fs - fs_bisection):.4f} against {TOLERANCE:g}")
    carried = res.get('ramp_last_carried')
    refused = res.get('ramp_first_refused')
    if carried is None or refused is None:
        fails.append("the ramp result no longer carries both raw edges of its "
                     "final interval (ramp_last_carried / ramp_first_refused), so "
                     "its reporting convention cannot be read")
    else:
        mid = 0.5 * (carried + refused)
        if abs(fs - mid) > 1e-12:
            fails.append(
                f"the ramp reports FS = {fs:.6f}, which is not the midpoint of its "
                f"final interval [{carried:.6f}, {refused:.6f}] (midpoint "
                f"{mid:.6f}). Every locked and published factor of safety here is a "
                f"bisection midpoint; reporting the lower edge, or flooring it to a "
                f"resolution, reads systematically low against them.")
        if not (carried <= fs <= refused):
            fails.append(
                f"the ramp reports FS = {fs:.6f}, outside its own final interval "
                f"[{carried:.6f}, {refused:.6f}]")
    evaluated = [t['F'] for t in res.get('trials', [])]
    if carried is None or not evaluated:
        fails.append("the ramp result carries no step record, so its "
                     "never-past-the-limit property cannot be read")
    else:
        top = max(evaluated)
        if top > carried + RAMP_DF_INIT + 1e-9:
            fails.append(
                f"the ramp evaluated F = {top:.4f}, which is "
                f"{top - carried:.4f} above the highest strength it carried "
                f"({carried:.4f}) and more than its {RAMP_DF_INIT:g} initial "
                f"increment. Its whole cost advantage is that it never solves far "
                f"past failure, which is where a load-controlled Newton trial is "
                f"at its most expensive.")
    return fails


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

    That last sentence used to be the whole guard against the gate coming back, and
    it guarded nothing: reintroducing ``_NR_RUNAWAY`` and its increment-abandoning
    test left this check passing, because loosening the OTHER three knobs does not
    disturb a gate that is not one of them. Two things close that. The constant may
    not exist on the module at all, and the behavior it used to break is measured
    directly in :func:`check_displacement_bound`.
    """
    fails = []
    if hasattr(_fem, '_NR_RUNAWAY'):
        fails.append(
            "xslope.fem defines _NR_RUNAWAY again. That constant was a displacement "
            "gate inside the load-increment loop: it abandoned an increment once "
            "max|u| passed a multiple of the elastic response, which reported trials "
            "as failed at strengths where the same driver reaches equilibrium with a "
            "statically admissible stress field. A displacement limit belongs on the "
            "FINAL state (see _NR_DISP_FACTOR), never on the solver's path.")
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


def check_displacement_bound():
    """The two halves of "the slope stands at F", measured on the trial that needs both.

    Equilibrium alone is not the question. Griffiths & Lane 1 on a quad9 mesh at
    size 3.5 is the case that separates the two halves: at F = 1.400 the Newton
    driver reaches machine-precision statics — an out-of-balance four orders of
    magnitude under the force tolerance — in a state that has moved about a sixth of
    the 50 m model height. That is outside the small-strain kinematics the whole
    model is written in, so it is not a slope standing; it is the displacement
    upturn, and the limit is the strength below it, which is where Griffiths & Lane
    read their own Example 1.

    Both directions are asserted, so the check is a bracket on the answer and not
    just a label:

      * F = 1.396875 converges AND stays under the bound (max|u| below a tenth of
        the model height), so the boundary is above it;
      * F = 1.400 is FAILED with the signal 'displacement_limit' AND with its
        out-of-balance UNDER the force tolerance, so it was rejected on
        displacement and not on forces, and the boundary is below it.

    The second assertion is what makes this mutation-sensitive. Under the removed
    ``_NR_RUNAWAY`` gate the same trial fails inside the increment loop and comes
    back at the load-step floor with a large out-of-balance — same FAILED label,
    different signal, different force reading — so a reintroduced gate cannot pass
    this by producing the same verdict for the wrong reason.
    """
    fails = []
    fd = _fem_data(element_type='quad9', target_size=QUAD9_SIZE)
    height = float(np.max(fd['nodes'][:, 1]) - np.min(fd['nodes'][:, 1]))
    bound = 0.1 * height

    stands = _newton(fd, Q9_ADMISSIBLE)
    if not stands['converged']:
        fails.append(
            f"F = {Q9_ADMISSIBLE} was expected to stand on the quad9 mesh and came "
            f"back {stands['verdict']} ({stands.get('exit_reason')}, "
            f"max|u| = {stands['max_displacement']:.4g} against a bound of "
            f"{bound:.4g}); the Newton factor of safety is no longer bracketed here")
    elif stands['max_displacement'] > bound:
        fails.append(
            f"F = {Q9_ADMISSIBLE} converged at max|u| = "
            f"{stands['max_displacement']:.4g}, past the {bound:.4g} bound it is "
            f"meant to sit under — the displacement bound is not being applied")

    over = _newton(fd, Q9_INADMISSIBLE)
    if over['converged']:
        fails.append(
            f"F = {Q9_INADMISSIBLE} came back CONVERGED at max|u| = "
            f"{over['max_displacement']:.4g} on a {height:.0f} m model "
            f"({100 * over['max_displacement'] / height:.1f}% of its height): a "
            f"state the small-strain model cannot represent is being reported as a "
            f"slope standing")
    else:
        if over.get('diverging_signal') != 'displacement_limit':
            fails.append(
                f"F = {Q9_INADMISSIBLE} failed for the wrong reason — signal "
                f"{over.get('diverging_signal')!r}, exit "
                f"{over.get('exit_reason')!r}, after {over['iterations']} "
                f"iterations. It is meant to reach equilibrium and be refused on "
                f"displacement; failing it inside the solve is the removed "
                f"_NR_RUNAWAY gate's signature.")
        if over['unbalanced_force_ratio'] >= 1e-3:
            fails.append(
                f"F = {Q9_INADMISSIBLE} was refused with an out-of-balance of "
                f"{over['unbalanced_force_ratio']:.3e}, at or above the 1e-3 force "
                f"tolerance: it never reached equilibrium, so this is a force-side "
                f"failure wearing a displacement label")
        if over['max_displacement'] <= bound:
            fails.append(
                f"F = {Q9_INADMISSIBLE} was refused on displacement at max|u| = "
                f"{over['max_displacement']:.4g}, which is UNDER the {bound:.4g} "
                f"bound")
    return fails


def check_env_override_announces_itself():
    """A stale XSLOPE_FEM_SOLVER may not redefine the default in silence.

    The variable swaps the driver for every solve in the process. Left set from an
    earlier session it would recompute every factor of safety on the non-default
    driver with nothing in the output to say so, and an explicit argument must stay
    silent because a caller that named the solver already knows.
    """
    import io
    import contextlib
    fails = []
    saved_env = os.environ.get('XSLOPE_FEM_SOLVER')
    saved_flag = _fem._FEM_SOLVER_ENV_ANNOUNCED
    try:
        os.environ['XSLOPE_FEM_SOLVER'] = 'newton'
        _fem._FEM_SOLVER_ENV_ANNOUNCED = False
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            got = _fem.resolve_fem_solver()
        if got != 'newton':
            fails.append(f"XSLOPE_FEM_SOLVER=newton resolved to {got!r}")
        if 'XSLOPE_FEM_SOLVER' not in buf.getvalue():
            fails.append(
                "XSLOPE_FEM_SOLVER=newton switched the default driver without "
                "printing a word: a stale shell variable would silently recompute "
                "every factor of safety on the non-default solver")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            _fem.resolve_fem_solver('newton')
        if buf.getvalue().strip():
            fails.append(
                "an EXPLICIT fem_solver='newton' printed a warning; only the "
                "environment override is meant to announce itself")
    finally:
        if saved_env is None:
            os.environ.pop('XSLOPE_FEM_SOLVER', None)
        else:
            os.environ['XSLOPE_FEM_SOLVER'] = saved_env
        _fem._FEM_SOLVER_ENV_ANNOUNCED = saved_flag
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
          "certifies its own stress field in force AND in displacement, neither "
          "the loading path nor the step control changes a verdict, and the "
          "environment override cannot swap the driver in silence. The monotonic "
          "ramp reaches the same limit along one warm-started history, reports it "
          "on the bisection's midpoint convention, and never solves past it.")


if __name__ == '__main__':
    main()
