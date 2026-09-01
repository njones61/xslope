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

Finally, the driver refuses what it cannot represent. Its return map is plain
Mohr-Coulomb, so a Hoek-Brown, power-curve or suction-bearing model must raise
rather than be solved on the wrong strength envelope; each is built as a real model
and the refusal message has to name the feature.

Reinforcement (SPIKE.md, "REINFORCEMENT") is the one feature the driver DOES carry,
and it is checked from both sides. What it carries: a reinforced bisection lands on
the strength measured for its mesh, its bar forces stay inside the capacity the
embedment develops, the bars it marks failed sit at that capacity, and a capped bar
holds the same force at two strength reductions — the reinforcement keeps its
capacity while the soil loses its. What it refuses: post-peak softening, which both
locked reinforced benchmarks declare, and piles.

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

# The reinforced case (SPIKE.md, "REINFORCEMENT"). The FEM reinforcement sample —
# six geogrid layers in an embankment — on a mesh coarse enough for the suite. It
# ships declaring post-peak softening (t_res 600 against t_max 800), which the
# Newton driver refuses, so the solved model has t_res unset: elastic-perfectly-
# plastic reinforcement, which is the model's own default and what the bar law the
# driver implements describes.
REINF_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files'
               / 'xslope_reinforce_fem.xlsx')
REINF_SIZE = 4.0
REINF_F_MIN, REINF_F_MAX = 1.2, 1.8
# Measured on this checkout, 2026-09-01: the Newton bisection reports 1.6078 on the
# interval [1.60313, 1.61250]. Locked to this mesh rather than to a published value,
# because the published reinforced values (1.497, 1.496) are defined on the SHIPPED
# models, which declare softening and which this driver refuses — see SPIKE.md,
# "REINFORCEMENT — results".
REINF_LOCKED_FS = 1.6078
REINF_TOL = 0.02


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
    failures += check_unsupported_features_refuse()
    failures += check_reinforcement_refusals()
    failures += check_reinforcement()
    failures += check_cohesionless_solve()
    failures += check_cohesionless_seed_depth()
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


def check_unsupported_features_refuse():
    """The Newton driver must REFUSE a model whose strength it cannot represent.

    The driver's return map is plain Mohr-Coulomb. A model built on a different
    strength envelope, or one that carries suction strength, would still solve on
    it — the guard is the only thing standing between such a model and a factor of
    safety computed with the wrong envelope, which is the failure mode that looks
    right. Three of the guarded features had no test:

      * Hoek-Brown (``option='hb'``), which the viscoplastic loop handles by
        inverting Balmer's curve for a tangent (c_i, phi_i) at each Gauss point;
      * the power-curve envelope (``option='pow'``), linearized the same way;
      * matric suction, which adds an apparent cohesion the return map knows
        nothing about.

    Each is built as a real model — the same coarse mesh the rest of this file
    uses, with the material re-declared on that envelope, or with a suction angle
    passed in — and the driver has to raise ``NotImplementedError`` naming the
    feature. The message is asserted, not just the exception type, so a guard that
    fires for some OTHER reason cannot pass this.

    The control matters as much as the guard: the same three models run on the
    default viscoplastic driver, which is where they belong.
    """
    import copy
    fails = []
    slope_data, mesh = _mesh()

    def _variant(**over):
        d = copy.deepcopy(slope_data)
        d['materials'][0].update(over)
        return build_fem_data(d, mesh)

    plain = build_fem_data(copy.deepcopy(slope_data), mesh)
    suction_name = plain['material_names'][0]

    cases = [
        ("Hoek-Brown",
         _variant(option='hb', hb_sci=30000.0, hb_gsi=60.0, hb_mi=10.0, hb_d=0.0),
         {}, 'hb_flag_by_elem'),
        ("power-curve",
         _variant(option='pow', pow_a=1.1, pow_b=0.9, pow_c=2.0, pow_d=4.0),
         {}, 'pow_flag_by_elem'),
        ("matric suction", plain, {'suction_phi_b': {suction_name: 15.0}}, None),
    ]

    for label, fd, kwargs, flag_key in cases:
        # The model really is on that feature, so a guard firing is a guard doing
        # its job rather than an unrelated refusal.
        if flag_key is not None and not np.any(fd[flag_key]):
            fails.append(
                f"the {label} test model carries no {flag_key} elements, so it "
                f"cannot exercise the guard — the model, not the guard, is broken")
            continue
        try:
            solve_fem(fd, F=1.0, fem_solver='newton', **kwargs)
        except NotImplementedError as exc:
            if label not in str(exc):
                fails.append(
                    f"the Newton driver refused the {label} model with a message "
                    f"that does not name it: {exc}")
        else:
            fails.append(
                f"the Newton driver ACCEPTED a {label} model. Its return map is "
                f"plain Mohr-Coulomb, so it would report a factor of safety "
                f"computed on the wrong strength envelope — the one failure mode "
                f"the guard exists for, because the answer looks right.")

    # Control: the viscoplastic driver, which does carry all three, must not be
    # refusing them. Two iterations is enough — the guard would raise before any.
    for label, fd, kwargs, _flag in cases:
        try:
            solve_fem(fd, F=1.0, max_iterations=2, **kwargs)
        except NotImplementedError as exc:
            fails.append(f"the DEFAULT viscoplastic driver now refuses the "
                         f"{label} model, which it is meant to handle: {exc}")
    return fails


def _reinf_fem_data(path=None, target_size=REINF_SIZE, soften=False,
                    keep_lines=None):
    """A reinforced model, meshed with its reinforcement lines as constraints.

    ``soften=False`` unsets ``t_res`` on every line, which is the model's own
    "no post-peak drop" default and the law the Newton bar element implements:
    tension-only, elastic-perfectly-plastic at the capacity the embedment develops.
    ``soften=True`` leaves the file as it ships, which is what the guard refuses.

    ``keep_lines`` is a slice of the reinforcement stack, 0-based, keeping only
    those layers — a leaner reinforcement layout over the same soil.
    """
    from xslope.mesh import extract_constraint_line_geometry
    slope_data = load_slope_data(str(path or REINF_MODEL))
    if keep_lines is not None:
        _all = slope_data.get('reinforcement_lines', []) or []
        slope_data['reinforcement_lines'] = [_all[i] for i in keep_lines]
    if not soften:
        for line in slope_data.get('reinforcement_lines', []) or []:
            line['t_res'] = float('nan')
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    mesh = build_mesh_from_polygons(
        get_material_polygons(slope_data, reinf_lines=lines),
        target_size=target_size, element_type='tri6', lines=lines,
        element_size_1d=slope_data.get('element_size_1d'))
    return build_fem_data(slope_data, mesh)


def check_reinforcement_refusals():
    """Bars are carried. Softening bars and piles are not, and they say so.

    The scope line is the whole point of the guard: the Newton driver's bar element
    is tension-only and elastic-perfectly-plastic, so a model whose bars declare a
    RUPTURE residual below the capacity their embedment develops is describing a law
    the driver does not implement. That law changes between solves — the post-peak
    set is a converged-state fixed point that re-opens the iteration each time it
    grows — which is where a Newton continuation is weakest.

    This matters more than an ordinary guard because the model it refuses is the
    shipped FEM reinforcement sample: both of the repository's locked reinforced FEM
    benchmarks declare t_res = 600 against t_max = 800. If the guard stopped firing,
    the driver would silently solve them with the bars holding their peak capacity
    forever and return a factor of safety that is too high and looks right.
    """
    fails = []
    shipped = _reinf_fem_data(soften=True)
    n_soft = int(np.count_nonzero(
        np.isfinite(np.asarray(shipped['t_res_by_1d_elem']))
        & (np.asarray(shipped['t_res_by_1d_elem'])
           < np.asarray(shipped['t_allow_by_1d_elem']) - 1e-12)))
    if n_soft == 0:
        fails.append(
            "the shipped reinforcement sample no longer declares post-peak "
            "softening on any bar element, so the softening refusal is untested "
            "here — the model, not the guard, has changed")
    else:
        try:
            solve_fem(shipped, F=1.0, fem_solver='newton')
        except NotImplementedError as exc:
            if 'softening' not in str(exc):
                fails.append(f"the Newton driver refused the softening model with a "
                             f"message that does not name softening: {exc}")
        else:
            fails.append(
                f"the Newton driver ACCEPTED a model declaring post-peak softening "
                f"on {n_soft} bar element(s). Its bar element is elastic-perfectly-"
                f"plastic, so those bars would hold their peak capacity forever and "
                f"the factor of safety would come back too high.")

    # Piles are beam elements with rotational degrees of freedom and a bending
    # capacity. Narrowing the guard to bars must not have narrowed it past them.
    pile_model = (Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files'
                  / 'xslope_piles_fem.xlsx')
    if pile_model.exists():
        fd = _reinf_fem_data(path=pile_model, target_size=4.0)
        if not fd.get('n_pile_elements', 0):
            fails.append("the pile test model carries no pile elements, so the pile "
                         "refusal is untested")
        else:
            try:
                solve_fem(fd, F=1.0, fem_solver='newton')
            except NotImplementedError as exc:
                if 'pile' not in str(exc):
                    fails.append(f"the Newton driver refused the pile model with a "
                                 f"message that does not name piles: {exc}")
            else:
                fails.append(
                    "the Newton driver ACCEPTED a pile model. Pile nodes carry a "
                    "rotational degree of freedom, which the driver's displacement "
                    "bound would then be comparing against a length, and the beam's "
                    "bending capacity is not in its element at all.")
    return fails


def check_reinforcement():
    """A reinforced Newton trial: the answer, and the bars' own evidence.

    The reinforcement sample on a coarse mesh, with softening unset so the Newton
    driver will take it. Four things, and only the first is about the number:

      * the Newton bisection lands on the strength measured for this mesh. Break the
        capacity cap — let a bar carry more than its embedment develops — and the
        slope holds at strengths it should not, so this moves.
      * the diagnostics are populated at full 1D-element length and physically
        consistent: no force outside [0, t_allow], every element the failed mask
        marks sitting AT its capacity, and the post-peak set empty, which is the only
        value reachable while softening is refused.
      * the bars agree with the viscoplastic driver's, which is the independent
        reading — the two solvers assemble the same bar from the same fem_data but
        drive it by completely different iterations.
      * the capacity is NOT reduced by F. Only the soil strength is, which is the
        vendor convention and what the limit-equilibrium engine does. A bar at
        capacity therefore reports the SAME force at two different strength
        reductions; dividing t_allow by F would make it fall.
    """
    fails = []
    fd = _reinf_fem_data()
    n_1d = len(fd['elements_1d'])
    t_allow = np.asarray(fd['t_allow_by_1d_elem'], dtype=float)
    if n_1d == 0 or not np.any(t_allow > 1e-9):
        return ["the reinforced test model carries no bar element with capacity, so "
                "nothing below tests the reinforcement"]

    res = solve_ssrm(fd, F_min=REINF_F_MIN, F_max=REINF_F_MAX, tolerance=0.01,
                     force_tol=1e-3, fem_solver='newton', max_iterations=16000,
                     capture_failure_state=False)
    fs = res.get('FS')
    if fs is None:
        return [f"the reinforced Newton bisection returned no factor of safety "
                f"({res.get('message') or res.get('error', 'no message')})"]
    if abs(fs - REINF_LOCKED_FS) > REINF_TOL:
        fails.append(
            f"reinforced: FS = {fs:.4f} is {abs(fs - REINF_LOCKED_FS):.4f} from the "
            f"{REINF_LOCKED_FS:.4f} measured for this mesh, outside {REINF_TOL:g}. "
            f"A reinforced factor of safety that moves is most likely the bar law: "
            f"a broken capacity cap lets the bars carry the slope past where they "
            f"can.")

    sol = solve_fem(fd, F=fs, fem_solver='newton', max_iterations=16000)
    forces = np.asarray(sol.get('forces_1d', []), dtype=float)
    at_cap = np.asarray(sol.get('failed_1d_elements', []), dtype=bool)
    softened = np.asarray(sol.get('softened_1d_elements', []), dtype=bool)
    if forces.shape != (n_1d,) or at_cap.shape != (n_1d,) or softened.shape != (n_1d,):
        return fails + [
            f"the reinforced Newton solution's diagnostics are not at 1D-element "
            f"length ({forces.shape}, {at_cap.shape}, {softened.shape} against "
            f"{n_1d}); every reader of forces_1d indexes them by element"]
    if np.any(forces < -1e-9):
        fails.append(f"{int(np.count_nonzero(forces < -1e-9))} bar(s) report a "
                     f"COMPRESSIVE force; the bar element is tension-only")
    over = forces > t_allow + 1e-6
    if np.any(over):
        fails.append(
            f"{int(over.sum())} bar(s) carry more than the capacity their embedment "
            f"develops (worst {np.max(forces[over] - t_allow[over]):.3g} over a cap "
            f"of {t_allow[over][np.argmax(forces[over] - t_allow[over])]:.3g}). The "
            f"capacity cap is what makes a reinforced slope failable at all.")
    if at_cap.any():
        off = np.abs(forces[at_cap] - t_allow[at_cap])
        if np.max(off) > 1e-6 * max(1.0, float(t_allow.max())):
            fails.append(f"a bar the failed mask marks is not at its capacity "
                         f"(worst gap {np.max(off):.3g})")
    if np.any(softened):
        fails.append(f"{int(softened.sum())} bar(s) are reported in the post-peak "
                     f"set on a driver that refuses softening")

    # The independent reading: the viscoplastic driver, same model, same strength.
    vp = solve_fem(fd, F=fs, max_iterations=16000)
    f_vp = np.asarray(vp.get('forces_1d', []), dtype=float)
    if f_vp.shape == forces.shape and f_vp.max() > 0:
        gap = float(np.max(np.abs(forces - f_vp)) / f_vp.max())
        if gap > 0.35:
            fails.append(
                f"the two drivers' bar forces differ by {gap:.2f} of the largest "
                f"force at F = {fs:.4f}. They assemble the SAME bar from the same "
                f"fem_data, so a gap this size is a difference in the bar law and "
                f"not in the soil solve.")

    # Capacity is not reduced by F: a bar at capacity holds it at every strength.
    lo = solve_fem(fd, F=max(1.0, fs - 0.2), fem_solver='newton', max_iterations=16000)
    f_lo = np.asarray(lo.get('forces_1d', []), dtype=float)
    both = at_cap & np.asarray(lo.get('failed_1d_elements', []), dtype=bool)
    if both.any():
        drift = np.max(np.abs(forces[both] - f_lo[both]))
        if drift > 1e-6 * max(1.0, float(t_allow.max())):
            fails.append(
                f"a bar at capacity reports {drift:.4g} more force at one strength "
                f"reduction than at another. Only the SOIL strength is reduced — the "
                f"reinforcement keeps its capacity, which is the vendor convention "
                f"and what the LEM applies — so a capped bar's force cannot depend "
                f"on F.")
    elif not at_cap.any():
        fails.append(
            "no bar reaches its capacity at the critical strength on this mesh, so "
            "the capped branch of the bar law — the branch that decides a reinforced "
            "factor of safety — is not exercised here")
    return fails


def check_cohesionless_solve():
    """A cohesionless base must not manufacture failures at strengths that stand.

    The reinforcement sample at its LOCKED tri6/2.0 mesh — 2,101 elements — over a
    lower material with c = 0, phi = 37. This is the case that broke the Newton
    driver (SPIKE.md, "THE COHESIONLESS SOLVE"): from a zero start the driver could
    not carry any part of gravity at F = 1.3, refusing the load at every increment
    down to the 1/64 floor, while converging at F = 1.2 below it and F = 1.45 above
    it. A verdict sequence that goes stands / fails / stands is wrong however the
    physics is read.

    Three assertions, and the middle one is the one the fix has to earn:

      * F = 1.45 converges. It is what makes the F = 1.3 verdict falsifiable without
        any reference to the other driver: strength reduction only ever RAISES c and
        tan(phi) as F falls, so an admissible stress field at F = 1.45 in equilibrium
        with full gravity is admissible at every lower F with the same load, and the
        lower-bound theorem then says the slope stands there. Its own yield-violation
        reading is asserted, because that is what makes the field admissible rather
        than merely converged.
      * F = 1.3 converges, to force equilibrium AND with a stress field on the yield
        surface. On the driver before the viscoplastic predictor this comes back
        FAILED at a load factor of 0.00 after 599 iterations, so the check fails
        there — which is the point of it.
      * ... and it converges to a state BELOW F = 1.45's, since a stronger soil under
        the same load cannot move further. A trial that passed the first two by
        finding some other equilibrium would not pass this.

    The mesh is the expensive one on purpose. The defect is invisible at tri6/4.0,
    marginal at tri6/3.0 and decisive here, so a cheaper mesh would lock nothing.
    """
    fails = []
    fd = _reinf_fem_data(target_size=2.0)
    if len(fd['elements']) < 1500:
        return [f"the cohesionless specimen meshed to {len(fd['elements'])} elements, "
                f"too coarse to carry the defect this check exists for"]

    above = _newton(fd, 1.45)
    if not above['converged']:
        return fails + [
            f"F = 1.45 on the cohesionless specimen came back {above['verdict']} "
            f"({above.get('exit_reason')}); it is the converged upper trial the whole "
            f"lower-bound argument below rests on, so nothing else here can be read"]
    if above['nr_max_yield_violation'] > 1e-6:
        fails.append(
            f"F = 1.45 converged with a worst yield violation of "
            f"{above['nr_max_yield_violation']:.3e} of the local strength; the "
            f"lower-bound argument needs an ADMISSIBLE field, not just an "
            f"equilibrated one")

    below = _newton(fd, 1.3)
    if not below['converged']:
        fails.append(
            f"F = 1.3 came back {below['verdict']} ({below.get('exit_reason')}, "
            f"signal {below.get('diverging_signal')!r}) after {below['iterations']} "
            f"iterations at a load factor of {below.get('nr_load_factor', 0.0):.2f}, "
            f"on a slope the SAME driver reports standing at F = 1.45 with a yield "
            f"violation of {above['nr_max_yield_violation']:.1e}. Strength reduction "
            f"only raises c and tan(phi) as F falls, so that field stands at 1.3 too: "
            f"this is the solver refusing a load the slope carries.")
    else:
        if below['unbalanced_force_ratio'] >= 1e-3:
            fails.append(
                f"F = 1.3 is reported CONVERGED at an out-of-balance of "
                f"{below['unbalanced_force_ratio']:.3e}, at or above the 1e-3 force "
                f"tolerance")
        if below['nr_max_yield_violation'] > 1e-6:
            fails.append(
                f"F = 1.3 is reported CONVERGED with a worst yield violation of "
                f"{below['nr_max_yield_violation']:.3e} of the local strength: the "
                f"stress field is outside the surface it is meant to lie on")
        if below['max_displacement'] > above['max_displacement']:
            fails.append(
                f"F = 1.3 converged to max|u| = {below['max_displacement']:.4g}, "
                f"further than F = 1.45's {above['max_displacement']:.4g}. A stronger "
                f"soil under the same gravity cannot move more; the two trials have "
                f"not found the same solution branch.")
    return fails


# The three-layer variant of the reinforcement sample — the same soil and the same
# six-layer geometry with only layers 1, 3 and 5 carrying capacity. It is the model
# that needed the ADAPTIVE predictor rung (SPIKE.md, "THE ADAPTIVE PREDICTOR"): the
# short fixed rungs rescue nothing on it, and the whole bisection sits 0.045 below
# the strength an admissible field proves it carries.
#
# The deciding trial, and where the numbers come from. An independent referee sharing
# no code with either driver certified this model standing at F = 1.20625 — a stress
# field in equilibrium with full gravity and nowhere outside the Mohr-Coulomb surface
# (SPIKE.md, "THE THREE-LAYER DISAGREEMENT"). It is the trial that decides the
# bisection on the bracket 1.0-1.6: converge here and the answer is 1.2109, refuse it
# and the answer is 1.1641.
THREE_LAYER_KEEP = (0, 2, 4)
THREE_LAYER_STANDS = 1.20625
THREE_LAYER_FAILS = 1.215625


def check_cohesionless_seed_depth():
    """A seed has to be deep enough, and the model that proves it.

    On the three-layer variant the plastic walk the corrector needs is twenty-odd
    THOUSAND viscoplastic iterations long, so a predictor on a fixed budget cannot
    reach it however the budget is chosen — six fixed budgets up to 32,000 were tried
    and each corrected to a different non-equilibrium state. What closes it is a rung
    that runs while the walk is still progressing and stops when it is not.

    Two assertions, and neither refers to the other driver:

      * F = 1.20625 converges to an ADMISSIBLE field — force equilibrium inside the
        trial tolerance AND a stress field on the yield surface. That is a statically
        admissible field in equilibrium with full gravity, so by the lower-bound
        theorem the slope stands there, and the bisection that refuses it is 0.045
        low. On the driver before the adaptive rung this comes back FAILED, which is
        the point of the check.
      * F = 1.215625, one bisection cell above it, still FAILS. The predictor is a
        seed generator and never a verdict: if a longer walk could talk the corrector
        into standing anywhere, this is where it would show.

    The mesh is the coarse tri6/4.0 one on purpose. Unlike the locked-mesh check
    above, the defect here is not a mesh-refinement effect — it is the length of the
    plastic walk — so the cheap mesh carries it.
    """
    fails = []
    fd = _reinf_fem_data(target_size=4.0, keep_lines=THREE_LAYER_KEEP)
    n_cap = int(np.count_nonzero(
        np.asarray(fd['t_allow_by_1d_elem']) > 1e-9))
    if n_cap == 0:
        return ["the three-layer specimen carries no bar with capacity, so nothing "
                "here is a reinforced measurement"]

    stands = _newton(fd, THREE_LAYER_STANDS)
    if not stands['converged']:
        fails.append(
            f"F = {THREE_LAYER_STANDS} came back {stands['verdict']} "
            f"({stands.get('exit_reason')}, signal "
            f"{stands.get('diverging_signal')!r}) after {stands['iterations']} "
            f"iterations at a load factor of "
            f"{stands.get('nr_load_factor', 0.0):.2f}, on {stands.get('nr_predictor_iterations', 0)} "
            f"predictor iterations. An admissible stress field in equilibrium with "
            f"full gravity exists at this strength, so the slope stands here and the "
            f"bisection that refuses it reads about 0.045 low: the seed the corrector "
            f"was handed is not deep enough.")
    else:
        if stands['unbalanced_force_ratio'] >= 1e-3:
            fails.append(
                f"F = {THREE_LAYER_STANDS} is reported CONVERGED at an out-of-balance "
                f"of {stands['unbalanced_force_ratio']:.3e}, at or above the 1e-3 "
                f"force tolerance")
        if stands['nr_max_yield_violation'] > 1e-6:
            fails.append(
                f"F = {THREE_LAYER_STANDS} is reported CONVERGED with a worst yield "
                f"violation of {stands['nr_max_yield_violation']:.3e} of the local "
                f"strength: the field is outside the surface it is meant to lie on, "
                f"so it proves nothing about whether the slope stands")

    above = _newton(fd, THREE_LAYER_FAILS)
    if above['converged']:
        fails.append(
            f"F = {THREE_LAYER_FAILS} came back CONVERGED. The predictor supplies a "
            f"starting state and must never supply a verdict; a converged trial one "
            f"bisection cell above the certified limit means a deeper seed is now "
            f"talking the corrector into standing where it did not.")
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
          "the loading path nor the step control changes a verdict, Hoek-Brown, "
          "power-curve, suction, softening-bar and pile models are refused by "
          "name, a reinforced bisection lands on its measured strength with every "
          "bar force inside the capacity its embedment develops, a cohesionless "
          "base does not manufacture a failure at a strength the driver's "
          "own converged answer above it proves standing, and the "
          "environment override cannot swap the driver in silence. The monotonic "
          "ramp reaches the same limit along one warm-started history, reports it "
          "on the bisection's midpoint convention, and never solves past it.")


if __name__ == '__main__':
    main()
