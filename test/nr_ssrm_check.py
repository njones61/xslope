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
Mohr-Coulomb, so a Hoek-Brown or power-curve model must raise rather than be solved
on the wrong strength envelope; each is built as a real model and the refusal
message has to name the feature.

Matric suction (SPIKE.md, "MATRIC SUCTION") is carried, not refused. Fredlund's
credit is an apparent cohesion on the same linear envelope, so it enters as a
per-Gauss-point cohesion and the check is written on that FIELD as well as on the
answer: the two drivers must build the same field, the cap must bound it, the trial
strength must reduce it, and the answer must move on the cap's value.

Reinforcement (SPIKE.md, "REINFORCEMENT") is the one feature the driver DOES carry,
and it is checked from both sides. What it carries: a reinforced bisection lands on
the strength measured for its mesh, its bar forces stay inside the capacity the
embedment develops, the bars it marks failed sit at that capacity, and a capped bar
holds the same force at two strength reductions — the reinforcement keeps its
capacity while the soil loses its. What it refuses: post-peak softening, which both
locked reinforced benchmarks declare, and piles.

The Rankine tension cutoff (SPIKE.md, "THE TENSION CUTOFF") is the second feature
the driver carries rather than refuses, and it matters because new materials are
authored with ``t_cut = 0``, so the old refusal covered the default case. It is
checked on four legs: the capped return map's invariants and its branch histogram,
which has to show the Mohr-Coulomb / Rankine intersection edge actually executing;
a cap above the Mohr-Coulomb apex reproducing the uncapped return bit for bit,
since no admissible state can reach that apex; the cap being reduced by the trial
strength when ``tension_srf`` says so; and the factor of safety on a capped model,
from both drivers.

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
from xslope.hoekbrown import hb_constants, hb_tangent
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
    failures += check_min_slip_depth()
    failures += check_reinforcement_refusals()
    failures += check_curved_envelopes()
    failures += check_reinforcement()
    failures += check_cohesionless_solve()
    failures += check_cohesionless_seed_depth()
    failures += check_tension_cutoff(fem_data)
    failures += check_matric_suction()
    failures += check_k0_initial_stress()
    failures += check_pile_element()
    failures += check_piles()
    failures += check_softening()
    failures += check_env_override_announces_itself()
    failures += check_ramp(fem_data, results['newton'].get('FS'))
    failures += check_rescue_cost(fem_data, results['newton'])
    failures += check_factorization(fem_data)
    failures += check_corrector(fem_data)

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


# The min_slip_depth filter (SPIKE.md, "THE MINIMUM SLIP DEPTH"). RS2-66b is the
# embankment-on-soft-ground model the filter exists for: a cohesionless fill face
# whose skin mechanism governs the unfiltered answer, over a soft band that carries
# the deep-seated one. It is run here on a mesh coarser than the locked 3.0 m so the
# rather than a whole bisection, with the same K0 = 1 and
# the same strength-reduced tensile cap.
DEPTH_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'verification'
               / 'files' / 'rocscience' / 'rs2_66b.xlsx')
DEPTH_SIZE = 3.0
DEPTH_FILTER = 4.0       # m below the ground surface: below the fill skin
DEPTH_F = 1.10           # stands with the filter, fails without it, on both drivers
DEPTH_MAX_ITER = 16000


def _depth_fem_data():
    slope_data = load_slope_data(str(DEPTH_MODEL))
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=DEPTH_SIZE,
                                    element_type='tri6')
    return build_fem_data(slope_data, mesh)


def _depth_trial(fem_data, solver, min_slip_depth):
    return solve_fem(fem_data, F=DEPTH_F, fem_solver=solver,
                     max_iterations=DEPTH_MAX_ITER, max_disp_factor=None,
                     min_slip_depth=min_slip_depth, k0=1.0, tension_srf=True,
                     debug_level=0)


def check_min_slip_depth():
    """A skin that fails does not decide a depth-filtered trial, on either driver.

    The filter's whole content is that a mechanism shallower than the cutoff is not
    the question being asked. On the viscoplastic driver that is one rule — the deep
    nodes alone are read for the out-of-balance — and it is enough there, because
    that driver's operator is the elastic stiffness, which no amount of skin flow
    makes singular, and its convergence test is a RATIO that a creeping skin's own
    displacement flatters rather than defeats. The Newton driver decides in three
    more places, and until this round none of them was filtered: the load increment
    it could not carry, the merit its line search minimized, and the displacement
    bound its converged state had to pass. All three belonged to the skin, so a
    depth-filtered Newton run answered the question the filter was meant to exclude.

    The check is a bracket, not a label. At one strength on one model:

      * WITH the filter, BOTH drivers converge, and the Newton state passes the same
        force gate and the same yield reading every other trial passes;
      * WITHOUT it, BOTH drivers fail at that same strength — so the skin really is
        governing there and the filter really is what changes the answer;
      * the skin has moved FURTHER than the deep field in the converged filtered
        state, which is what "the skin fails but the trial stands" looks like, and
        the deep field is inside the small-strain bound the verdict is read on.

    The second and third assertions are what make it mutation-sensitive. Drop the
    filter from the increment decision — hand ``_nr_equilibrate`` ``deep_free=None``
    — and the first assertion fails: the Newton trial dies at the load-step floor at
    a strength the viscoplastic driver carries.
    """
    fails = []
    fd = _depth_fem_data()
    height = float(np.max(fd['nodes'][:, 1]) - np.min(fd['nodes'][:, 1]))

    filtered = {s: _depth_trial(fd, s, DEPTH_FILTER)
                for s in ('viscoplastic', 'newton')}
    unfiltered = {s: _depth_trial(fd, s, None)
                  for s in ('viscoplastic', 'newton')}

    for solver, r in filtered.items():
        if not r['converged']:
            fails.append(
                f"{solver}: F = {DEPTH_F} with min_slip_depth = {DEPTH_FILTER} came "
                f"back {r['verdict']} ({r.get('exit_reason')}, "
                f"{r['iterations']} iterations, out-of-balance "
                f"{r['unbalanced_force_ratio']:.3e}, max|u| "
                f"{r['max_displacement']:.4g} on a {height:.0f} m model) — the "
                f"depth filter is not reaching this driver's verdict")

    for solver, r in unfiltered.items():
        if r['converged']:
            fails.append(
                f"{solver}: F = {DEPTH_F} stands with NO filter as well, so this "
                f"model no longer separates the skin from the deep mechanism and "
                f"the filtered result above proves nothing")

    nr = filtered['newton']
    if nr['converged']:
        if nr['unbalanced_force_ratio'] >= 1e-3:
            fails.append(
                f"the filtered Newton trial was called converged at an "
                f"out-of-balance of {nr['unbalanced_force_ratio']:.3e}, at or above "
                f"the 1e-3 force tolerance")
        if nr['nr_oob_global'] <= nr['unbalanced_force_ratio']:
            fails.append(
                f"the filtered Newton trial's whole-model out-of-balance "
                f"({nr['nr_oob_global']:.3e}) is no larger than the filtered one "
                f"({nr['unbalanced_force_ratio']:.3e}), so the filter is not "
                f"excluding anything and this trial says nothing about it")
        deep, skin = nr.get('nr_max_disp_deep'), nr.get('nr_max_disp_skin')
        if skin is None or deep is None:
            fails.append(
                "the filtered Newton trial reported no deep/skin displacement "
                "split, so nothing says which field the bound was read on")
        else:
            if skin <= deep:
                fails.append(
                    f"the skin moved {skin:.4g} against the deep field's "
                    f"{deep:.4g}: the excluded skin is not the part that is "
                    f"failing here, so this trial does not test the rule")
            if deep > 0.1 * height:
                fails.append(
                    f"the filtered Newton trial converged with a DEEP displacement "
                    f"of {deep:.4g} on a {height:.0f} m model, past the tenth of "
                    f"the height the bound allows — the bound is not being applied "
                    f"to the field the verdict rests on")
    return fails


def _curved_variant(option, mesh=None, slope_data=None, **params):
    """The check file's own coarse model with its material re-declared curved."""
    import copy
    if slope_data is None or mesh is None:
        slope_data, mesh = _mesh()
    d = copy.deepcopy(slope_data)
    d['materials'][0].update(option=option, **params)
    return build_fem_data(d, mesh)


HB_PARAMS = dict(hb_sci=30000.0, hb_gsi=60.0, hb_mi=10.0, hb_d=0.0)
POW_PARAMS = dict(pow_a=2.0, pow_b=0.85, pow_c=0.0, pow_d=5.0)

# The mixed model: four materials on one mesh, two Mohr-Coulomb, one Hoek-Brown,
# one power curve. Its numbers are written out rather than imported for the same
# reason every other lock in this file writes its own out — an assertion that
# reads the code it is checking asserts nothing.
MIXED_FILE = (Path(__file__).resolve().parents[1] / 'docs' / 'fem'
              / 'files' / 'xslope_noncircular_fem.xlsx')
MIXED_SIZE = 4.0
MIXED_HB = dict(hb_sci=886.6288783215639, hb_gsi=55.0, hb_mi=12.0, hb_d=0.0)
MIXED_POW = dict(pow_a=0.8826850799708558, pow_b=0.85, pow_c=0.0, pow_d=0.0)
MIXED_F_STANDS = 1.10
MIXED_F_FAILS = 2.10


def _mixed_fem_data():
    """Two Mohr-Coulomb materials, one Hoek-Brown and one power curve, one mesh.

    Constructed rather than found: every curved-envelope model in the corpus is
    SINGLE-material, so nothing shipped exercises per-element dispatch (see
    SPIKE.md, "CURVED ENVELOPES"). The two curved materials are fitted to the
    strengths the file's own Mohr-Coulomb materials carry at mid-height
    overburden, so the model is a slope and not a shape.
    """
    import copy
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    sd = load_slope_data(str(MIXED_FILE))
    sd = copy.deepcopy(sd)
    sd['materials'][1].update(option='hb', **MIXED_HB)
    sd['materials'][2].update(option='pow', **MIXED_POW)
    mesh = build_mesh_from_polygons(get_material_polygons(sd),
                                    target_size=MIXED_SIZE, element_type='tri6')
    return sd, build_fem_data(sd, mesh)


def check_curved_envelopes():
    """Hoek-Brown and the power curve, per Gauss point, on both drivers.

    Four legs, and each fails on a different defect (see SPIKE.md, "CURVED
    ENVELOPES"):

      * the RETURN MAP on both envelopes — every returned state on the linearized
        Mohr-Coulomb surface and under its tensile cap, the principal ordering
        intact, no elastic state modified, and a branch histogram that reaches
        every region, because a fuzz that never lands on a region proves nothing
        about it;
      * the LINEARIZATION is the one the viscoplastic driver takes — the tangent
        the return was computed on is the tangent of the F-reduced envelope at the
        abscissa the RETURNED stress produces, which is the fixed point the
        viscoplastic loop converges to;
      * the F-REDUCTION divides the tangent and not the constants: at F = 2 the
        cohesive intercept and tan(phi) are exactly half what they are at F = 1,
        on both envelopes;
      * and the MIXED model — four materials, two Mohr-Coulomb, one Hoek-Brown,
        one power curve, on one mesh — solves to the same verdicts on both
        drivers, with the per-element dispatch asserted from the solve rather than
        inferred: every Mohr-Coulomb element on the plain path and every curved one
        on its own, in the same solve.
    """
    fails = []
    slope_data, mesh = _mesh()

    # ---- the return map, per envelope ---------------------------------------
    rng = np.random.default_rng(20260901)
    n = 40000
    scale = 1000.0
    # mb / s / a are DERIVED from GSI, mi and D through the package's own
    # hb_constants, not written out, so a defect in that derivation reaches this
    # fuzz as well as the model legs below.
    _mb, _s, _a = hb_constants(HB_PARAMS['hb_gsi'], HB_PARAMS['hb_mi'],
                               HB_PARAMS['hb_d'])
    for label, code, params in (
            ('Hoek-Brown', _fem._NR_ENV_HB,
             dict(hb_sci=np.full(n, HB_PARAMS['hb_sci']),
                  hb_mb=np.full(n, float(_mb)), hb_s=np.full(n, float(_s)),
                  hb_a=np.full(n, float(_a)))),
            ('power curve', _fem._NR_ENV_POW,
             dict(pow_a=np.full(n, POW_PARAMS['pow_a']),
                  pow_b=np.full(n, POW_PARAMS['pow_b']),
                  pow_cp=np.full(n, POW_PARAMS['pow_c']),
                  pow_d=np.full(n, POW_PARAMS['pow_d'])))):
        for cap_label, t_cut in (('no cap', None), ('t_cut = 0', 0.0),
                                 ('t_cut inert', 1e9)):
            mu = rng.uniform(2e4, 2e5, n)
            env = {'code': np.full(n, code, dtype=np.int8), 'F': np.full(n, 1.3)}
            for k in ('pow_a', 'pow_b', 'pow_cp', 'pow_d',
                      'hb_sci', 'hb_mb', 'hb_s', 'hb_a'):
                env[k] = params.get(k, np.ones(n) if k == 'pow_b' else np.zeros(n))
            grp = {'env': env, 'mu': mu,
                   'c_r': np.full(n, 50.0), 'snph': np.full(n, np.sin(0.6)),
                   'csph': np.full(n, np.cos(0.6))}
            if t_cut is not None:
                grp['t_cap'] = np.full(n, t_cut)
                grp['lam'] = 1.5 * mu
            sig_tr = rng.normal(0.0, scale, (n, 4))
            sig_tr[:, 2] *= 0.6
            sig_tr[rng.integers(0, 4, n) == 0] *= 8.0
            out, branch, c, sn, cs, _ = _fem._nr_envelope_return(grp, sig_tr.copy())

            ctr, half = 0.5 * (out[:, 0] + out[:, 1]), 0.5 * (out[:, 0] - out[:, 1])
            R = np.hypot(half, out[:, 2])
            P = -np.sort(-np.stack([ctr + R, ctr - R, out[:, 3]], 1), axis=1)
            s1, s2, s3 = P[:, 0], P[:, 1], P[:, 2]
            A, Bc = 0.5 * (1 + sn), 0.5 * (1 - sn)
            sc = np.maximum.reduce([np.abs(s1), np.abs(s3), np.abs(c * cs),
                                    np.full(n, scale)])
            worst = float(np.max((A * s1 - Bc * s3 - c * cs) / sc))
            tag = f"{label} / {cap_label}"
            if worst > 1e-10:
                fails.append(
                    f"{tag}: a returned state sits {worst:.3e} of the stress "
                    f"scale OUTSIDE the linearized Mohr-Coulomb surface — "
                    f"strength the point has not got")
            if t_cut is not None and np.isfinite(t_cut):
                over = float(np.max((s1 - t_cut) / sc))
                if over > 1e-10:
                    fails.append(f"{tag}: a returned major principal stress is "
                                 f"{over:.3e} of the scale above its cap")
            nbad = int(np.count_nonzero((s1 < s2 - 1e-9 * sc)
                                        | (s2 < s3 - 1e-9 * sc)))
            if nbad:
                fails.append(f"{tag}: {nbad} returned state(s) came back with the "
                             f"principal stresses out of order")
            el = branch == 0
            if el.any() and np.any(out[el] != sig_tr[el]):
                fails.append(f"{tag}: an ELASTIC state was modified by the return")
            # the histogram: a region the fuzz never reaches is a region it says
            # nothing about
            got = set(np.unique(branch).tolist())
            want = {0, 1, 2, 3} | ({4} if t_cut is None or t_cut > 1e8
                                   else {5, 8, 9})
            missing = sorted(want - got)
            if missing:
                fails.append(
                    f"{tag}: the fuzz never reached branch(es) "
                    + ", ".join(_fem._NR_BRANCH_NAMES[b] for b in missing)
                    + " — it proves nothing about them")
            if _fem._NR_TENS_FALLBACK in got:
                fails.append(f"{tag}: state(s) came back on the unresolved "
                             f"fallback, so the candidate region list is incomplete")

            # the linearization IS the envelope's tangent at the returned state
            sp = -0.5 * (out[:, 0] + out[:, 1])
            c2, sn2, cs2, _ = _fem._nr_envelope_step(env['code'], c, sn, cs, sp,
                                                     env, None)
            cref = max(1.0, float(np.abs(c2 * cs2).mean()))
            resid = float(max(np.max(np.abs(c2 * cs2 - c * cs)) / cref,
                              np.max(np.abs(sn2 - sn))))
            if resid > 1e-3:
                fails.append(
                    f"{tag}: the tangent the return was taken on is not the "
                    f"tangent of the envelope at the abscissa the RETURNED stress "
                    f"produces — the self-consistency residual is {resid:.3e}, so "
                    f"the returned circle does not touch the curve")

    # ---- the tangent is the ENVELOPE's, from the material table -------------
    # The strongest leg, and the one a driver-against-driver comparison could not
    # supply: the tangent the Newton path derives is compared against one computed
    # from scratch out of the MATERIAL's own columns, through the public
    # hb_tangent (which re-derives mb, s and a from GSI, mi and D) and through the
    # power curve's own formula. A parameter wired from the wrong per-element
    # array is self-consistent with itself and only this can see it.
    for label, option, params in (('Hoek-Brown', 'hb', HB_PARAMS),
                                  ('power curve', 'pow', POW_PARAMS)):
        F_env = 1.7
        fd = _curved_variant(option, mesh, slope_data, **params)
        prep = _fem._prepare_fem_model(fd)
        n_el = len(fd['elements'])
        env = _fem._nr_envelope_by_elem(fd, np.full(n_el, F_env))
        groups = _fem._nr_build_groups(prep, np.zeros(n_el), np.zeros(n_el), None,
                                       env_by_elem=env)
        g = groups[0]
        m = g['n']
        sp = np.linspace(0.0, 4000.0, m)
        c0 = np.full(m, 40.0)
        sn0, cs0 = np.full(m, np.sin(0.45)), np.full(m, np.cos(0.45))
        c, sn, cs, _ = _fem._nr_envelope_step(g['env']['code'], c0, sn0, cs0, sp,
                                              g['env'], None)
        if option == 'hb':
            sig_n = np.maximum(sp * cs0 ** 2 - c0 * sn0 * cs0, 0.0)
            c_i, phi_i = hb_tangent(sig_n, params['hb_sci'], params['hb_gsi'],
                                    params['hb_mi'], params['hb_d'])
            c_ex = c_i / F_env
            t_ex = np.tan(np.radians(phi_i)) / F_env
        else:
            s_n = np.maximum(sp, 0.0)
            ref = max(1.0, float(s_n.mean()))
            s_ef = np.maximum(s_n + params['pow_d'], 1e-4 * ref)
            t_ex = (params['pow_a'] * params['pow_b']
                    * s_ef ** (params['pow_b'] - 1.0)) / F_env
            c_ex = ((params['pow_a'] * s_ef ** params['pow_b'] + params['pow_c'])
                    / F_env) - s_n * t_ex
        d_c = float(np.max(np.abs(c - c_ex)) / max(1.0, float(np.abs(c_ex).max())))
        d_t = float(np.max(np.abs(sn / cs - t_ex)) / max(1e-3, float(np.abs(t_ex).max())))
        if d_c > 1e-9 or d_t > 1e-9:
            fails.append(
                f"{label}: the tangent the Newton path derives is not the "
                f"envelope the material declares — the cohesive intercept differs "
                f"by {d_c:.3e} and tan(phi) by {d_t:.3e} relative, against an "
                f"expectation computed from the material's own columns")

    # ---- the F-reduction divides the tangent, not the constants -------------
    for label, option, params in (('Hoek-Brown', 'hb', HB_PARAMS),
                                  ('power curve', 'pow', POW_PARAMS)):
        got = {}
        for F in (1.0, 2.0):
            fd = _curved_variant(option, mesh, slope_data, **params)
            prep = _fem._prepare_fem_model(fd)
            n_el = len(fd['elements'])
            envF = np.full(n_el, float(F))
            env = _fem._nr_envelope_by_elem(fd, envF)
            groups = _fem._nr_build_groups(
                prep, np.zeros(n_el), np.zeros(n_el), None, env_by_elem=env)
            g = groups[0]
            sp = np.full(g['n'], 250.0)
            c, sn, cs, _ = _fem._nr_envelope_step(
                g['env']['code'], np.full(g['n'], 10.0),
                np.full(g['n'], np.sin(0.5)), np.full(g['n'], np.cos(0.5)),
                sp, g['env'], None)
            got[F] = (c.copy(), (sn / cs).copy())
        for what, i in (('cohesive intercept', 0), ('tan(phi)', 1)):
            r = got[2.0][i] / np.maximum(got[1.0][i], 1e-30)
            if not np.allclose(r, 0.5, rtol=1e-9):
                fails.append(
                    f"{label}: the {what} of the F-reduced tangent at F = 2 is "
                    f"{float(r.min()):.6f} to {float(r.max()):.6f} of its value "
                    f"at F = 1, where strength reduction divides it by exactly 2 "
                    f"— the reduction is on the TANGENT and not on the envelope's "
                    f"own constants")

    # ---- the mixed model ----------------------------------------------------
    sd_mix, fd_mix = _mixed_fem_data()
    n_el = len(fd_mix['elements'])
    hbm = np.asarray(fd_mix['hb_flag_by_elem'], bool)
    pwm = np.asarray(fd_mix['pow_flag_by_elem'], bool)
    n_mc = int(n_el - hbm.sum() - pwm.sum())
    if not (hbm.any() and pwm.any() and n_mc):
        fails.append(
            f"the mixed model does not mix: {int(hbm.sum())} Hoek-Brown, "
            f"{int(pwm.sum())} power-curve and {n_mc} Mohr-Coulomb elements of "
            f"{n_el} — it cannot exercise per-element dispatch")
    else:
        stash = []
        _orig = _fem._nr_build_groups
        def _grab(*a, **k):
            g = _orig(*a, **k)
            stash.append(g)
            return g
        _fem._nr_build_groups = _grab
        try:
            sol = solve_fem(fd_mix, F=MIXED_F_STANDS, fem_solver='newton',
                            force_tol=1e-3, k0=sd_mix.get('k0'))
        finally:
            _fem._nr_build_groups = _orig
        emat = np.asarray(fd_mix['element_materials'])
        seen = {}
        for g in stash[-1]:
            code = (g['env']['code'] if g.get('env') is not None
                    else np.zeros(g['n'], np.int8))
            for e, cd in zip(g['e_idx'], code):
                seen.setdefault(int(emat[e]), set()).add(int(cd))
        want = {1: _fem._NR_ENV_MC, 2: _fem._NR_ENV_HB,
                3: _fem._NR_ENV_POW, 4: _fem._NR_ENV_MC}
        for mat, code in want.items():
            got = seen.get(mat, set())
            if got != {code}:
                fails.append(
                    f"material {mat} ({sd_mix['materials'][mat-1].get('name')}) "
                    f"was dispatched to envelope code(s) {sorted(got)} where the "
                    f"model declares {code}: the dispatch is per ELEMENT and a "
                    f"model that mixes strength models in one mesh has to send "
                    f"each element to its own")
        if not sol['converged']:
            fails.append(
                f"the mixed model does not stand at F = {MIXED_F_STANDS} on the "
                f"Newton driver ({sol.get('exit_reason')}), so the dispatch above "
                f"was read off a state that means nothing")
        elif sol.get('nr_max_yield_violation', 1.0) > 1e-6:
            fails.append(
                f"the mixed model's converged field at F = {MIXED_F_STANDS} is "
                f"{sol['nr_max_yield_violation']:.2e} of the local strength "
                f"outside its own yield surface")
        for F, must_stand in ((MIXED_F_STANDS, True), (MIXED_F_FAILS, False)):
            for drv in ('viscoplastic', 'newton'):
                r = solve_fem(fd_mix, F=F, fem_solver=drv, force_tol=1e-3,
                              max_iterations=8000, k0=sd_mix.get('k0'))
                if bool(r['converged']) != must_stand:
                    fails.append(
                        f"the mixed model at F = {F} on the {drv} driver came "
                        f"back {'FAILED' if must_stand else 'CONVERGED'} where "
                        f"the other strength model in the same mesh is what "
                        f"decides it — the two drivers must agree on this bracket")
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
    """Nothing about reinforcement is refused any more, and that is the assertion.

    The driver used to refuse post-peak softening, which is what both of the
    repository's locked reinforced FEM benchmarks declare — so the guard was covering
    the showcase. It is carried now (see :func:`check_softening`), and what stands in
    the guard's place here is the opposite claim: the shipped reinforcement sample,
    with its softening as the file gives it, must SOLVE on the Newton driver rather
    than raise.
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
            "softening on any bar element, so nothing here is a softening "
            "measurement — the model, not the driver, has changed")
        return fails
    try:
        sol = solve_fem(shipped, F=1.0, fem_solver='newton', max_disp_factor=None)
    except NotImplementedError as exc:
        fails.append(
            f"the Newton driver still REFUSES the shipped reinforcement sample, "
            f"which declares post-peak softening on {n_soft} bar element(s): {exc}")
    else:
        if 'softened_1d_elements' not in sol:
            fails.append("a solved softening model returns no post-peak set")

    # Piles used to be asserted here too. They are CARRIED now (see
    # :func:`check_piles`), so a refusal expected here would be holding the driver
    # to a limitation it no longer has.
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
        # ... and the rescue chain is what carries it. The cost policy budgets the
        # chain by rungs (SPIKE.md, "THE COST OF THE RESCUE"), so the budget has to
        # be able to decide a verdict where the chain is decisive — otherwise the
        # verdict-identity assertion in check_rescue_cost is testing a knob attached
        # to nothing. Take the rungs away here and F = 1.3 must go back to what the
        # cold attempt alone reaches, which is the load-step floor.
        stripped = solve_fem(fd, F=1.3, fem_solver='newton', max_disp_factor=None,
                             force_tol=1e-3, _nr_rescue_rungs=0)
        if stripped['converged']:
            fails.append(
                f"F = 1.3 still converged with the rescue chain budgeted to zero "
                f"rungs ({stripped['iterations']} iterations, out-of-balance "
                f"{stripped['unbalanced_force_ratio']:.3e}). The predictor is meant "
                f"to be what carries this trial, so a rung budget that cannot take "
                f"it away is not reaching the chain at all.")
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


# The tension cutoff (SPIKE.md, "THE TENSION CUTOFF"). Measured on this checkout:
# the coarse tri6 model with t_cut = 0 on its one material reads 1.375 on BOTH
# drivers, one bisection cell below the uncapped 1.39 — the cap is doing real work
# and the two drivers see the same work.
TCUT_FS = 1.375
# The tension_srf leg. A cap of 30 at F = 1.2 binds either way, so the two settings
# are distinguishable: reduced with the trial it holds the major principal stress
# at 30/1.2 = 25, held at its authored value it holds it at 30.
TSRF_F, TSRF_CAP = 1.2, 30.0


def check_tension_cutoff(fem_data):
    """The Rankine tensile cap: the second surface, and the factor it is reduced by.

    The Newton driver used to refuse any model with a finite ``t_cut``. Since
    new materials are authored with ``t_cut = 0`` that refusal covered the default
    case, and the cap is also the documented remedy for the illegal tension the
    three-layer adjudication measured on a cohesionless soil (SPIKE.md, "THE
    THREE-LAYER DISAGREEMENT"). Four things are asserted, and each one fails on a
    different defect:

      * **the return map**, capped. Random trial states at four friction angles and
        three caps — none, zero, and one above the Mohr-Coulomb apex — must come
        back satisfying BOTH surfaces, ordered, with elastic states untouched. The
        BRANCH HISTOGRAM is asserted too: the Mohr-Coulomb / Rankine intersection
        edge and the hydrostatic-tension return must actually execute, because a
        fuzz that never lands on a region proves nothing about it, and no state may
        come back on the unresolved fallback.
      * **an inert cap is not a cap.** Every Mohr-Coulomb-admissible state has
        sigma_1 <= c cot(phi), so a cap at or above that apex can never bind, and
        the capped map must reproduce the uncapped one BIT FOR BIT. That is the
        cheapest possible statement that the second surface does not perturb the
        first.
      * **the cap is reduced with the trial** when ``tension_srf`` is on, which is
        what makes the factor of safety a factor on the whole envelope rather than
        on its shear half. At F = 1.2 with a cap of 30 the converged state's major
        principal stress sits at exactly 25 with the reduction and exactly 30
        without it.
      * **the factor of safety**, on the same coarse model the rest of this file
        uses, with the cap on: both drivers, agreeing with each other and with the
        strength measured for this mesh.
    """
    import numpy as _np
    from xslope.fem import (_NR_BRANCH_NAMES, _NR_TENS_FALLBACK, _NR_MC_TENS,
                            _NR_TENS3, build_constitutive_matrix_4)
    fails = []
    D4 = build_constitutive_matrix_4(30000.0, 0.3)
    mu_v, lam_v = D4[2, 2], D4[0, 1]
    rng = np.random.default_rng(17)
    n, scale = 40000, 30.0
    seen = {}
    for phi_deg, c in ((0.0, 10.0), (20.0, 5.0), (35.0, 2.0), (45.0, 1.0)):
        snph, csph = np.sin(np.radians(phi_deg)), np.cos(np.radians(phi_deg))
        a, bc = 0.5 * (1.0 + snph), 0.5 * (1.0 - snph)
        apex = np.inf if snph < 1e-12 else c * csph / snph
        sig = rng.normal(0.0, scale, size=(n, 4))
        sig[:, 2] *= 0.5
        args = (np.full(n, c), np.full(n, snph), np.full(n, csph), np.full(n, mu_v))
        plain, br_plain = mc_return_map(sig, *args)
        for tag, T in (("t_cut=0", 0.0), ("t_cut=%g" % (0.35 * c * csph),
                                          0.35 * c * csph),
                       ("t_cut inert", (1e6 if not np.isfinite(apex)
                                        else 10.0 * apex))):
            out, br = mc_return_map(sig, *args, t_cap=np.full(n, T),
                                    lam=np.full(n, lam_v))
            ctr = 0.5 * (out[:, 0] + out[:, 1])
            half = 0.5 * (out[:, 0] - out[:, 1])
            r = np.sqrt(half * half + out[:, 2] ** 2)
            pn = -np.sort(-np.stack([ctr + r, ctr - r, out[:, 3]], axis=1), axis=1)
            sc = np.maximum(np.abs(pn[:, 0]), scale)
            where = f"phi={phi_deg:g}, {tag}"
            f_new = (a * pn[:, 0] - bc * pn[:, 2] - c * csph) / sc
            if f_new.max() > 1e-12:
                fails.append(f"{where}: a returned state is OUTSIDE the "
                             f"Mohr-Coulomb surface by {f_new.max():.2e} of the "
                             f"stress scale")
            t_new = (pn[:, 0] - T) / sc
            if t_new.max() > 1e-12:
                fails.append(f"{where}: a returned state is above the tensile cap "
                             f"by {t_new.max():.2e} of the stress scale — the "
                             f"Rankine surface is not being enforced")
            if (np.minimum(pn[:, 0] - pn[:, 1], pn[:, 1] - pn[:, 2])
                    / sc).min() < -1e-12:
                fails.append(f"{where}: a returned state has its principal "
                             f"stresses out of order")
            el = br == 0
            if el.any() and np.abs(out[el] - sig[el]).max() > 0.0:
                fails.append(f"{where}: an elastic point was modified")
            nfb = int((br == _NR_TENS_FALLBACK).sum())
            if nfb:
                fails.append(f"{where}: {nfb} state(s) came back on the unresolved "
                             f"fallback — no candidate active set was consistent, "
                             f"which means the region list is incomplete")
            if tag == "t_cut inert":
                if not (np.array_equal(out, plain) and np.array_equal(br, br_plain)):
                    fails.append(
                        f"{where}: a cap ABOVE the Mohr-Coulomb apex changed the "
                        f"return. No admissible state can reach sigma_1 = c "
                        f"cot(phi), so an inert cap has to be bit-identical to no "
                        f"cap at all")
            for code, cnt in zip(*np.unique(br, return_counts=True)):
                seen[int(code)] = seen.get(int(code), 0) + int(cnt)
    for code in (_NR_MC_TENS, _NR_TENS3):
        if not seen.get(code):
            fails.append(
                f"the '{_NR_BRANCH_NAMES[code]}' branch never executed in "
                f"{4 * 3 * n} returns, so nothing here tests it. Either the fuzz "
                f"no longer reaches that region or the region is unreachable — "
                f"both are defects in this check")

    # --- the cap is reduced with the trial strength ---------------------------
    n_el = len(fem_data['elements'])
    for tsrf, want in ((True, TSRF_CAP / TSRF_F), (False, TSRF_CAP)):
        sol = solve_fem(fem_data, F=TSRF_F, fem_solver='newton',
                        max_iterations=MAX_ITER, max_disp_factor=None,
                        tension_cap_by_elem=_np.full(n_el, TSRF_CAP),
                        tension_srf=tsrf, force_tol=1e-3)
        if not sol['converged']:
            fails.append(f"tension_srf={tsrf}: the F = {TSRF_F} trial did not "
                         f"converge, so the cap in force cannot be read off it")
            continue
        st = sol['stresses']
        sx, sy, txy = -st[:, 0], -st[:, 1], st[:, 2]
        s1 = 0.5 * (sx + sy) + np.sqrt((0.5 * (sx - sy)) ** 2 + txy ** 2)
        # Two-sided, and the two sides say different things. Above the cap would be
        # a cap not enforced; far below it would be a cap that never binds, which
        # would make this leg vacuous — the point is that BOTH settings bind here
        # and land on their own value. (The stresses are element averages over the
        # Gauss points, so the reading sits a hair under the cap the points are at.)
        got = float(s1.max())
        if got > want * (1.0 + 1e-6):
            fails.append(
                f"tension_srf={tsrf}: the converged field reaches a major "
                f"principal stress of {got:.6f}, above the {want:.6f} cap in "
                f"force")
        elif got < want * 0.99:
            fails.append(
                f"tension_srf={tsrf}: the largest major principal stress in the "
                f"converged field is {got:.6f}, well under the {want:.6f} the cap "
                f"in force allows — the cap does not bind, so this leg tests "
                f"nothing. With the reduction on the cap is t_cut/F and without "
                f"it t_cut; a cap that is not divided by the trial makes the "
                f"factor of safety a factor on the shear envelope only")
        if sol.get('nr_max_tension_violation', 0.0) > 1e-8:
            fails.append(
                f"tension_srf={tsrf}: the converged field violates its own cap by "
                f"{sol['nr_max_tension_violation']:.2e} of the local strength")

    # --- and the factor of safety, both drivers -------------------------------
    name = fem_data['material_names'][0]
    got = {}
    for solver in ('viscoplastic', 'newton'):
        res = solve_ssrm(fem_data, F_min=F_MIN, F_max=F_MAX, tolerance=TOLERANCE,
                         max_iterations=MAX_ITER, fem_solver=solver,
                         tension_cutoff_by_material={name: 0.0},
                         capture_failure_state=False)
        got[solver] = res.get('FS')
        if got[solver] is None:
            fails.append(f"{solver}: no factor of safety on the capped model")
        elif abs(got[solver] - TCUT_FS) > TOLERANCE:
            fails.append(
                f"{solver}: with t_cut = 0 the factor of safety is "
                f"{got[solver]:.4f}, {abs(got[solver] - TCUT_FS):.4f} from the "
                f"{TCUT_FS} measured for this mesh")
    if None not in got.values() and abs(got['viscoplastic'] - got['newton']) > TOLERANCE:
        fails.append(
            f"the two drivers disagree on the CAPPED model: viscoplastic "
            f"{got['viscoplastic']:.4f} vs Newton {got['newton']:.4f}")
    return fails



# Matric suction (SPIKE.md, "MATRIC SUCTION"). The model is RS2 Part 4 VP102 at
# t = 60 s — docs/verification/rs2.md, benchmark RS2-P4-VP102-t-60-c3: its own
# transient seepage field on its own stored mesh (1,118 tri6 elements), an
# unsaturated friction angle of 37 deg supplied by the TEST TAG rather than the
# file, k0 = 1, tension_srf on, an SSR search zone, and a locked FS of
# 1.779 +/- 0.02. The tag's settings are written out rather than read through
# run_tests, so the assertion holds the mapping instead of following it.
SUC_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'verification' / 'files'
             / 'rocscience' / 'vp102t_60.xlsx')
SUC_PHI_B = {'Material 1': 37.0}
SUC_K0 = 1.0
SUC_MAX_ITER = 16000
# One trial either side of the two drivers' measured answers on this model
# (viscoplastic 1.7707, Newton 1.7926). The LOWER one is the load-bearing choice:
# the same model with the credit off is locked at 1.713 (its own c2 row), so a
# driver that loses the suction fails a trial it has to carry, and the leg cannot
# pass by accident.
SUC_F_STANDS, SUC_F_FAILS = 1.74, 1.82
# The cap legs. The suction reaches 91.4 stress units on this model and is
# positive at 474 Gauss points, so a cap of 5 binds at 421 of them — most of the
# credit — and a cap of 20 binds at 294 and still leaves enough to carry the
# trial. Both are asserted: a cap that changes every answer is as uninformative
# as one that changes none.
SUC_CAP_BINDS, SUC_CAP_LOOSE = 5.0, 20.0
# The pair of strengths the F-reduction is read across. BOTH must be trials the
# Newton driver CARRIES: a failing trial fires the viscoplastic predictor, whose
# own group build calls the same helper, and the capture would then hold two
# drivers' fields end to end rather than one driver's at two strengths.
SUC_F_LOW, SUC_F_HIGH = 1.0, 1.5


def _suction_fem_data():
    """The suction benchmark, on the mesh its stored seepage field is written on.

    ``u = seep`` means the nodal pore pressures belong to a particular mesh, so
    the model's own mesh is used rather than a fresh one — build_fem_data would
    silently zero the field if the node count differed.
    """
    slope_data = load_slope_data(str(SUC_MODEL))
    return build_fem_data(slope_data, slope_data['mesh'])


def _suction_field(fem_data, F, solver, phi_b=None, cap=None):
    """The apparent-cohesion field a DRIVER actually computes, in its own order.

    Captured at ``_suction_apparent_cohesion``, the one helper both drivers call,
    by running a single solver iteration through the driver's own entry point.
    Reading the field off a group the check built itself would prove only that the
    check can build a group; this proves what the solve is running on.
    """
    got = []
    original = _fem._suction_apparent_cohesion

    def spy(tanphib, s_capped, finv):
        out = original(tanphib, s_capped, finv)
        got.append(np.asarray(out).copy())
        return out

    _fem._suction_apparent_cohesion = spy
    try:
        solve_fem(fem_data, F=F, max_iterations=1, max_disp_factor=None,
                  tension_srf=True, k0=SUC_K0, fast_kernel=False,
                  suction_phi_b=(SUC_PHI_B if phi_b is None else phi_b),
                  suction_cap=cap, fem_solver=solver)
    finally:
        _fem._suction_apparent_cohesion = original
    return np.concatenate(got) if got else np.zeros(0)


def check_matric_suction():
    """Fredlund's apparent cohesion, on the Newton driver.

    The Newton driver used to refuse any model carrying a positive ``phi_b``,
    which put six locked vendor benchmarks out of reach for a credit that needs no
    new yield surface: ``min(max(0, -u_w), s_cap) tan(phi_b) / F`` is an APPARENT
    COHESION on the same linear envelope the return map already solves, and the
    return map already takes ``c`` per Gauss point.

    That makes the whole of this feature plumbing — which pore-pressure field,
    which cap, which F — so the check is written on the field and not only on a
    factor of safety. Five things are asserted, and each fails on a different
    defect:

      * **the two drivers compute the same field**, point for point, captured at
        the one helper they share by running a solver iteration through each
        driver's own entry point. A plumbing difference is the only thing that can
        show up here, which is the point.
      * **the credit is reduced by the trial strength**, as the viscoplastic
        driver reduces it and as docs/fem/overview.md states: the field at F = 1.5
        must be exactly two thirds of the field at F = 1 wherever the strength is
        actually reduced. The model carries an SSR search zone, so the points held
        at full strength must NOT scale, and both facts are asserted.
      * **``s_cap`` bounds the credit.** At F = 1 a cap of 5 must hold the largest
        credit to ``5 tan(phi_b)``, where the uncapped field on the same model
        reaches eighteen times that.
      * **the answer.** One trial either side of the measured limit, on BOTH
        drivers: F = 1.74 stands, F = 1.82 fails.
      * **the credit is load-bearing, and the cap moves it.** With the credit off,
        F = 1.74 must FAIL on both drivers — otherwise the trial above proves
        nothing about suction. With the binding cap it must fail too, and with the
        loose cap it must stand again, so what moved the answer is the cap's VALUE
        and not the presence of a cap keyword.
    """
    fails = []
    fd = _suction_fem_data()
    tanb = float(np.tan(np.radians(SUC_PHI_B['Material 1'])))

    # ---- (1) the two drivers' fields ---------------------------------------
    vp1 = _suction_field(fd, SUC_F_LOW, 'viscoplastic')
    nr1 = _suction_field(fd, SUC_F_LOW, 'newton')
    if nr1.size == 0:
        return ["the Newton driver computed no apparent cohesion at all on "
                f"{SUC_MODEL.name}, which declares an unsaturated friction angle "
                f"— the model, not the driver, may be the problem"]
    n_credited = int(np.count_nonzero(nr1 > 0.0))
    if n_credited == 0 or float(nr1.max()) <= 0.0:
        fails.append(
            f"{SUC_MODEL.name} credits suction at no Gauss point, so nothing "
            f"below is a measurement of matric suction")
    if vp1.shape != nr1.shape:
        fails.append(
            f"the two drivers built apparent-cohesion fields of different "
            f"lengths ({vp1.size} viscoplastic against {nr1.size} Newton), so "
            f"they are not walking the same Gauss points")
    else:
        d = float(np.max(np.abs(vp1 - nr1)))
        if d > 1e-12:
            fails.append(
                f"the two drivers' matric-suction apparent cohesion differs by "
                f"{d:.3e} at worst over {nr1.size} Gauss points. They share the "
                f"formula, so a difference here is a plumbing difference — the "
                f"signed pore-pressure field, the cap, or the trial factor")

    # ---- (2) the credit is reduced by F ------------------------------------
    nr_hi = _suction_field(fd, SUC_F_HIGH, 'newton')
    want = SUC_F_HIGH / SUC_F_LOW
    if nr_hi.shape != nr1.shape:
        fails.append(
            f"the apparent-cohesion field at F = {SUC_F_HIGH} has {nr_hi.size} "
            f"entries against {nr1.size} at F = {SUC_F_LOW}, so the two cannot be "
            f"compared point for point — both strengths must be trials the driver "
            f"CARRIES, since a failing one fires the viscoplastic predictor and "
            f"the capture then holds two drivers' fields")
    else:
        live = nr1 > 0.0
        ratio = nr1[live] / nr_hi[live]
        if not np.any(np.abs(ratio - want) < 1e-9):
            fails.append(
                f"the matric-suction apparent cohesion is NOT reduced by the "
                f"trial strength: raising F from {SUC_F_LOW} to {SUC_F_HIGH} left "
                f"the credit at a ratio of at most {float(ratio.max()):.6f} where "
                f"the viscoplastic driver divides it by {want:g}. The credit "
                f"enters the reduced envelope alongside c' and tan(phi'), so an "
                f"unreduced one reports a factor of safety on a stronger soil "
                f"than the trial actually has")
        _known = (np.abs(ratio - want) < 1e-9) | (np.abs(ratio - 1.0) < 1e-9)
        if not np.all(_known):
            fails.append(
                f"some credited Gauss points scale by neither 1 nor {want:g} "
                f"between F = {SUC_F_LOW} and F = {SUC_F_HIGH}; the only two "
                f"factors this model can produce are the SSR zone's full strength "
                f"and the reduced strength")

    # ---- (3) s_cap bounds it ------------------------------------------------
    capped = _suction_field(fd, SUC_F_LOW, 'newton', cap=SUC_CAP_BINDS)
    ceiling = SUC_CAP_BINDS * tanb
    if float(capped.max()) > ceiling * (1.0 + 1e-9):
        fails.append(
            f"s_cap does not bound the credit: with a cap of {SUC_CAP_BINDS:g} "
            f"the largest apparent cohesion is {float(capped.max()):.4f}, above "
            f"the {ceiling:.4f} that cap allows")
    if float(nr1.max()) <= ceiling:
        fails.append(
            f"the cap of {SUC_CAP_BINDS:g} does not bind on this model — the "
            f"uncapped credit peaks at {float(nr1.max()):.4f}, at or below the "
            f"{ceiling:.4f} the cap permits — so the leg above proves nothing")

    # ---- (4) and (5) the verdicts, on both drivers --------------------------
    def _verdict(F, phi_b=None, cap=None, solver='newton'):
        sol = solve_fem(fd, F=F, max_iterations=SUC_MAX_ITER, max_disp_factor=None,
                        tension_srf=True, k0=SUC_K0, fast_kernel=False,
                        failure_criterion='hybrid',
                        suction_phi_b=(SUC_PHI_B if phi_b is None else phi_b),
                        suction_cap=cap, fem_solver=solver)
        return bool(sol['converged'])

    for solver in ('viscoplastic', 'newton'):
        if not _verdict(SUC_F_STANDS, solver=solver):
            fails.append(
                f"{solver}: F = {SUC_F_STANDS} came back FAILED on "
                f"{SUC_MODEL.name}, where both drivers put the limit near 1.78 "
                f"with the suction credit and the vendor lock is 1.779")
        if _verdict(SUC_F_FAILS, solver=solver):
            fails.append(
                f"{solver}: F = {SUC_F_FAILS} came back CONVERGED, above the "
                f"limit both drivers measure on this model and above the vendor's "
                f"own 1.779")
        # The credit off: the same trial must fail, or the pair above is not a
        # measurement of suction.
        if _verdict(SUC_F_STANDS, phi_b={}, solver=solver):
            fails.append(
                f"{solver}: F = {SUC_F_STANDS} stands with the suction credit "
                f"switched OFF, so the trial above is not evidence that the "
                f"credit is carried — this model without suction is locked at "
                f"1.713")
        # The cap's VALUE, not its presence, is what moves the answer.
        if _verdict(SUC_F_STANDS, cap=SUC_CAP_BINDS, solver=solver):
            fails.append(
                f"{solver}: F = {SUC_F_STANDS} stands under a suction cap of "
                f"{SUC_CAP_BINDS:g}, which bounds the credit at 421 of the 474 "
                f"Gauss points that carry it — the cap is being ignored in the "
                f"solve even where the field says it is applied")
        if not _verdict(SUC_F_STANDS, cap=SUC_CAP_LOOSE, solver=solver):
            fails.append(
                f"{solver}: F = {SUC_F_STANDS} FAILS under a suction cap of "
                f"{SUC_CAP_LOOSE:g}, which leaves enough credit to carry it. "
                f"A cap keyword must move the answer by its VALUE, not by being "
                f"present at all")
    return fails


# The K0 initial stress (SPIKE.md, "K0 INITIAL STRESS"). The vendor lock this
# check reproduces is RS2-27 at its coarsest tagged mesh — docs/verification/rs2.md,
# `vp036.xlsx`, tri6 at target 1.5, k0 = 1, tension_srf off, FS 1.373 +/- 0.02. The
# tag's settings are written out rather than read through run_tests, so the
# assertion holds the mapping instead of following it.
K0_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'verification' / 'files'
            / 'rocscience' / 'vp036.xlsx')
K0_SIZE = 1.5
K0_LOCKED_FS = 1.373
K0_LOCK_TOL = 0.02
K0_F_MIN, K0_F_MAX = 1.1, 1.6
K0_MAX_ITER = 16000
# The level-ground block, from test/k0_level_ground_check.py, which is where the
# K0 field has a closed form. Machine-precision tolerances: the field is an exact
# equilibrium, so the only thing between the solver and zero is round-off.
K0_LG_MAX_U = 1e-10
K0_LG_REL_STRESS = 1e-9
K0_LG_MAX_ITER = 2


def _k0_level_ground():
    """The level-ground module's own model builders, loaded from beside this file."""
    import importlib.util
    path = Path(__file__).resolve().parent / 'k0_level_ground_check.py'
    spec = importlib.util.spec_from_file_location('k0_level_ground_check', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_k0_initial_stress():
    """The at-rest initial stress, on the driver that used to refuse it.

    K0 was the last refusal standing between this driver and the vendor corpus:
    every `fem_ssrm` lock in `docs/verification/` whose model carries a tensile cap
    also carries `k0=1`, so 142 locked benchmarks were unreachable for want of it.
    Four legs, and each fails on a different defect:

      * **level ground is exact.** On flat ground the field sigma'_v = -gamma z,
        sigma'_h = sigma'_z = K0 sigma'_v satisfies equilibrium identically for any
        K0, so there is nothing to redistribute: the Newton corrector must reach
        equilibrium on its first iteration, leave the mesh where it is, reproduce
        the imposed stresses and yield nowhere. Run dry and with the water table at
        the surface, where the K0 relation holds between EFFECTIVE stresses.
      * **the out-of-plane component is asserted, not assumed.** sigma_z does not
        appear in the in-plane equilibrium at all, so a wrong out-of-plane initial
        stress moves no node and changes no sigma_x: it would pass every leg above
        and still evaluate the yield surface on a state the model is not in. It is
        read here twice — directly off the groups the driver builds, against the
        analytic K0 sigma_v at each Gauss point, and through the von Mises stress
        of the K0 = 1 dry state, which is hydrostatic and therefore zero exactly
        when the third component is right.
      * **the two drivers build the same in-situ state.** The pre-equilibration is
        a full-strength solve whose (u, eps^p) every trial then starts from, so if
        the drivers disagree there, every trial is answering a different question.
        Compared on the reconstructed Gauss-point stress field rather than on
        either driver's reporting path.
      * **a locked vendor benchmark.** RS2-27, K0 = 1: both drivers inside the
        published tolerance and inside the bisection tolerance of each other.
    """
    import numpy as _np
    from xslope.fem import _prepare_fem_model, _nr_build_groups
    fails = []
    lg = _k0_level_ground()

    # ---- (1) and (2a): level ground, exactly ------------------------------
    fem_dry, fem_wet = lg._build(False), lg._build(True)
    for k0, wet in ((1.0, False), (1.0, True)):
        fd = fem_wet if wet else fem_dry
        tag = f"K0={k0:g} {'wet' if wet else 'dry'}"
        prep = _prepare_fem_model(fd, k0=k0)
        sv_exp, sh_exp = lg._expected(fd, prep, k0, wet)
        sol = solve_fem(fd, F=1.0, k0=k0, max_iterations=2000, fast_kernel=False,
                        _prepared=prep, fem_solver='newton')
        if not sol['converged']:
            fails.append(f"K0 level ground, {tag}: the Newton driver did not reach "
                         f"equilibrium on a field that is an exact one")
            continue
        ref = max(float(sv_exp.max()), 1e-30)
        st = sol['stresses']
        errs = {'sigma_y': float(_np.abs(st[:, 1] - sv_exp).max()) / ref,
                'sigma_x': float(_np.abs(st[:, 0] - sh_exp).max()) / ref,
                'tau_xy': float(_np.abs(st[:, 2]).max()) / ref}
        for name, err in errs.items():
            if err > K0_LG_REL_STRESS:
                fails.append(
                    f"K0 level ground, {tag}: {name} is off by {err:.3e} relative "
                    f"— the recovered stress is not the imposed K0 field")
        if float(sol['max_displacement']) > K0_LG_MAX_U:
            fails.append(
                f"K0 level ground, {tag}: max|u| = {sol['max_displacement']:.3e} "
                f"> {K0_LG_MAX_U:.0e} — the mesh moved under a field that is an "
                f"exact equilibrium")
        if int(sol['iterations']) > K0_LG_MAX_ITER:
            fails.append(
                f"K0 level ground, {tag}: equilibrium took {sol['iterations']} "
                f"iterations, not at most {K0_LG_MAX_ITER} — the K0 field is being "
                f"redistributed, so it is not in equilibrium as built")
        n_plastic = int(_np.count_nonzero(sol['plastic_elements']))
        if n_plastic:
            fails.append(
                f"K0 level ground, {tag}: {n_plastic} yielded element(s) under an "
                f"equilibrium state well inside the Mohr-Coulomb envelope")
        if not wet:
            # K0 = 1 dry is HYDROSTATIC, so the deviatoric stress is zero exactly
            # when sigma_z is right — and a wrong sigma_z shows up nowhere else on
            # level ground, because the out-of-plane component carries no in-plane
            # equilibrium and no in-plane stress.
            vm = float(_np.abs(st[:, 3]).max()) / ref
            if vm > K0_LG_REL_STRESS:
                fails.append(
                    f"K0 level ground, {tag}: the state is not hydrostatic — von "
                    f"Mises stress reaches {vm:.3e} of the overburden where K0 = 1 "
                    f"makes it zero. The out-of-plane initial stress is wrong.")

    # ---- (2b) the initial stress field itself, out-of-plane included -------
    k0_probe = 0.5
    prep = _prepare_fem_model(fem_dry, k0=k0_probe)
    n_el = len(fem_dry['elements'])
    groups = _nr_build_groups(prep, fem_dry['c_by_elem'],
                              _np.radians(fem_dry['phi_by_elem']),
                              None, None, k0=k0_probe)
    worst_h = worst_z = worst_v = 0.0
    n_gp = 0
    y_top = float(_np.max(fem_dry['nodes'][:, 1]))
    for grp, sg in zip(groups, prep['gp_groups_static']):
        s0 = grp.get('sig0')
        if s0 is None:
            fails.append("the K0 initial stress is not attached to the Newton "
                         "groups at all on a model that declares K0")
            break
        xy = lg._gauss_point_xy(fem_dry, prep)
        sv_an = _np.array([-lg.GAMMA * (y_top - xy[e][g][1]) for e, g in sg['pairs']])
        ref = max(float(_np.abs(sv_an).max()), 1e-30)
        n_gp += len(sv_an)
        worst_v = max(worst_v, float(_np.abs(s0[:, 1] - sv_an).max()) / ref)
        worst_h = max(worst_h, float(_np.abs(s0[:, 0] - k0_probe * sv_an).max()) / ref)
        worst_z = max(worst_z, float(_np.abs(s0[:, 3] - k0_probe * sv_an).max()) / ref)
    for name, err in (('vertical', worst_v), ('in-plane horizontal', worst_h),
                      ('OUT-OF-PLANE', worst_z)):
        if err > K0_LG_REL_STRESS:
            fails.append(
                f"the K0 initial stress the Newton driver builds is wrong in its "
                f"{name} component by {err:.3e} of the overburden, over {n_gp} "
                f"Gauss points, against the analytic K0 = {k0_probe:g} field")

    # ---- (3) the two drivers build the same in-situ state ------------------
    slope_data = load_slope_data(str(K0_MODEL))
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=K0_SIZE, element_type='tri6')
    fd = build_fem_data(slope_data, mesh)
    prep = _prepare_fem_model(fd, k0=1.0)

    def _gp_sigma(state):
        """sigma = sigma_0 + D (B u - eps^p), from first principles."""
        u = _np.asarray(state['u'], dtype=float)
        out = []
        for sg, ev in zip(prep['gp_groups_static'], state['evp']):
            B, D4, dof = sg['B'], sg['D4'], sg['dof']
            ep = _np.asarray(ev, dtype=float)
            eps = (B @ u[dof][:, :, None])[:, :, 0]
            eps4 = _np.empty((len(eps), 4))
            eps4[:, :3] = eps - ep[:, :3]
            eps4[:, 3] = -ep[:, 3]
            sig = (D4 @ eps4[:, :, None])[:, :, 0]
            sv0 = _np.array([prep['sv0_gp'][e][g] for e, g in sg['pairs']])
            ug = _np.array([prep['u_gp'][e][g] for e, g in sg['pairs']])
            sv = -sv0 + ug
            sh = 1.0 * sv          # K0 = 1 on this model
            z = _np.zeros_like(sv)
            out.append(sig + _np.stack([sh, sv, z, sh], axis=1))
        return _np.concatenate(out, axis=0)

    scale = max(float(max(max(row) for row in prep['sv0_gp'] if len(row))), 1e-30)
    states = {}
    for solver in ('viscoplastic', 'newton'):
        sol = solve_fem(fd, F=1.0, k0=1.0, max_iterations=K0_MAX_ITER,
                        max_disp_factor=None, force_tol=1e-3, early_exit=True,
                        failure_criterion='hybrid', fast_kernel=False,
                        _prepared=prep, fem_solver=solver)
        if not sol['converged'] or sol['_k0_state'] is None:
            fails.append(f"the {solver} driver could not establish the in-situ "
                         f"state on {K0_MODEL.name} at full strength")
        states[solver] = sol
    if all(s['_k0_state'] is not None for s in states.values()):
        d = _gp_sigma(states['viscoplastic']['_k0_state']) - \
            _gp_sigma(states['newton']['_k0_state'])
        rms = float(_np.sqrt(_np.mean(d ** 2))) / scale
        if rms > 1e-2:
            fails.append(
                f"the two drivers' K0 in-situ stress fields differ by {rms:.3e} "
                f"RMS relative to the overburden — every trial starts from this "
                f"state, so a disagreement here is a disagreement about the model")
        du = abs(states['viscoplastic']['max_displacement']
                 - states['newton']['max_displacement'])
        ref_u = max(states['viscoplastic']['max_displacement'], 1e-30)
        if du / ref_u > 0.05:
            fails.append(
                f"the two drivers' in-situ redistribution differs by "
                f"{du / ref_u:.1%} in max|u| "
                f"({states['viscoplastic']['max_displacement']:.4g} against "
                f"{states['newton']['max_displacement']:.4g})")

    # ---- (4) the locked vendor benchmark ----------------------------------
    fs = {}
    for solver in ('viscoplastic', 'newton'):
        res = solve_ssrm(fd, F_min=K0_F_MIN, F_max=K0_F_MAX, tolerance=0.01,
                         max_iterations=K0_MAX_ITER, k0=1.0, tension_srf=False,
                         fem_solver=solver, capture_failure_state=False)
        fs[solver] = res.get('FS')
        if fs[solver] is None:
            fails.append(f"{solver}: no factor of safety on {K0_MODEL.name} "
                         f"({res.get('error', 'no message')})")
        elif abs(fs[solver] - K0_LOCKED_FS) > K0_LOCK_TOL:
            fails.append(
                f"{solver}: FS = {fs[solver]:.4f} on {K0_MODEL.name} is "
                f"{abs(fs[solver] - K0_LOCKED_FS):.4f} from the locked "
                f"{K0_LOCKED_FS:g}, outside the published tolerance "
                f"{K0_LOCK_TOL:g}")
    if all(v is not None for v in fs.values()) and \
            abs(fs['viscoplastic'] - fs['newton']) > 0.01:
        fails.append(
            f"the two drivers disagree on {K0_MODEL.name} with K0 = 1: "
            f"viscoplastic {fs['viscoplastic']:.4f} against Newton "
            f"{fs['newton']:.4f}")
    return fails


# The pile cases (SPIKE.md, "PILES"), both at the mesh their own test tags carry.
#
# The FEM-3 tutorial's sheet pile wall is the one that holds an end — its tip is
# FIXED — so it is what the fixity and capacity legs run on.
PILE_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'tutorials' / 'files'
              / 'xslope_pile_wall.xlsx')
PILE_SIZE = 2.0
# The FEM pile sample: two pile rows, both capacities finite, both ends free. Its
# published lock is FS 1.380 (docs/fem/samples.md), and the two drivers bracket it
# identically — trial for trial, verdict for verdict, over the whole bisection. The
# pair below is that bracket, which holds the same fact as the bisection for four
# solves instead of two whole runs.
PILES_MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files'
               / 'xslope_piles_fem.xlsx')
PILES_SIZE = 2.0
PILES_LOCKED_FS = 1.380
PILES_STANDS, PILES_FAILS = 1.375, 1.384375


def _pile_fem_data(path=None, target_size=PILE_SIZE, element_size_1d=None):
    """A pile model, meshed with its pile axis as a constraint line.

    A reinforcement-only extraction drops the pile axes, and without them in the
    mesh there are no beam elements to extract — so the model would solve, silently,
    with no pile in it at all.
    """
    from xslope.mesh import (extract_constraint_line_geometry,
                             extract_point_constraints, extract_size_regions)
    slope_data = load_slope_data(str(path or PILE_MODEL))
    if element_size_1d is not None:
        slope_data['element_size_1d'] = element_size_1d
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    mesh = build_mesh_from_polygons(
        get_material_polygons(slope_data, reinf_lines=lines),
        target_size=target_size, element_type='tri6', lines=lines,
        element_size_1d=slope_data.get('element_size_1d'),
        point_constraints=extract_point_constraints(slope_data),
        size_regions=extract_size_regions(slope_data))
    return build_fem_data(slope_data, mesh)


def _beam_case(n_node, theta, L, EA, EI, V_cap=np.inf, M_cap=np.inf):
    """A one-element ``fem_data`` carrying exactly what _nr_build_piles reads."""
    from xslope.fem import _beam_local_stiffness, _beam_rotation
    nd = 3 * n_node
    Kl = _beam_local_stiffness(EA, EI, L, n_node == 3)
    T = _beam_rotation(np.cos(theta), np.sin(theta), n_node)
    dof = np.zeros((1, 9), dtype=np.int64)
    dof[0, :nd] = np.arange(nd)
    return {"n_pile_elements": 1, "K_global_pile_elems": [T.T @ Kl @ T],
            "K_local_by_pile_elem": [Kl], "dof_indices_pile": dof,
            "n_dof_by_pile_elem": np.array([nd]),
            "cos_theta_pile": np.array([np.cos(theta)]),
            "sin_theta_pile": np.array([np.sin(theta)]),
            "elem_length_by_pile_elem": np.array([L]),
            "EI_by_pile_elem": np.array([EI]), "EA_by_pile_elem": np.array([EA]),
            "V_cap_by_pile_elem": np.array([V_cap]),
            "M_cap_by_pile_elem": np.array([M_cap])}


def _nr_translational_max(u, tidx):
    """max|u| over the translational degrees of freedom, for the checks."""
    v = u if tidx is None else np.asarray(u)[tidx]
    return float(np.max(np.abs(v))) if np.size(v) else 0.0


def check_pile_element():
    """The beam element itself, before any soil: the tangent, and the closed forms.

    Two measurements, and neither needs a slope.

    **The consistent tangent is exact**, on the elastic branch, on the shear cap and
    across the plastic HINGE alike. Every action is linear in the element
    displacement, the shear clip is affine, and the hinge's plastic rotation
    ``p_A = S_AA^-1 (m_A - sign(m_A) M_cap)`` is affine too — so the analytic
    element tangent, with the rotational block condensed, must equal a central
    difference of the element's own internal force to round-off. 120 random
    elements, random orientation, length and rigidities, two-node and three-node,
    with the caps set from the actions the drawn displacement actually produces so
    that all eight combinations of (shear at its cap, end 1 hinged, end 2 hinged)
    are reached. A probe that straddles a branch boundary is skipped rather than
    averaged, and the histogram is asserted to have reached every branch — a check
    that never lands on a hinged branch would prove nothing about it.

    **An isolated beam reproduces the closed forms.** A six-element cantilever
    under a transverse end load must deflect PL^3/3EI and rotate PL^2/2EI at the
    tip; simply supported under a central load it must deflect PL^3/48EI. The
    assembly, the rotation into global coordinates and the residual are the
    round's own, so a corrupted bending stiffness cannot pass this.
    """
    from xslope.fem import _nr_build_piles, _nr_pile_force
    fails = []
    rng = np.random.default_rng(20260901)

    worst = 0.0
    branches = {}
    n_ok = 0
    for trial in range(120):
        n_node = 2 if trial % 2 == 0 else 3
        nd = 3 * n_node
        th = rng.uniform(0, 2 * np.pi)
        L = float(np.exp(rng.uniform(np.log(0.5), np.log(20.0))))
        EA = float(np.exp(rng.uniform(np.log(1e4), np.log(1e8))))
        EI = float(np.exp(rng.uniform(np.log(1e3), np.log(1e7))))
        u = rng.normal(0.0, 1.0, nd) * L * 1e-3
        pg0 = _nr_build_piles(_beam_case(n_node, th, L, EA, EI))[0]
        act = pg0['act'][0] @ u
        V_cap = abs(act[0]) * rng.choice([0.3, 0.8, 1.5, 3.0])
        # ONE moment capacity per element, as the file states it and as the
        # viscoplastic path reads it, drawn against one end or the other so that
        # both ends, either end and neither end are all reached.
        M_cap = abs(act[1 + rng.integers(2)]) * rng.choice([0.3, 0.8, 1.5, 3.0])
        if rng.random() < 0.25:
            V_cap = np.inf
        if rng.random() < 0.25:
            M_cap = np.inf
        pg = _nr_build_piles(_beam_case(n_node, th, L, EA, EI, V_cap, M_cap))[0]

        def _branch(g):
            # The branch the tangent is built on: the shear at its cap, and which
            # element ends actually took a plastic release.
            return (bool(g['_yV'][0]), bool(g['_p_rot'][0, 0]),
                    bool(g['_p_rot'][0, 1]))

        _nr_pile_force(pg, u, True)
        K = pg['_Ke'][0].copy()
        branch = _branch(pg)
        h = 1e-5 * float(np.max(np.abs(u)))
        D = np.zeros((nd, nd))
        straddled = False
        for j in range(nd):
            up, um = u.copy(), u.copy()
            up[j] += h
            um[j] -= h
            fp = _nr_pile_force(pg, up, False)[0].copy()
            b1 = _branch(pg)
            fm = _nr_pile_force(pg, um, False)[0].copy()
            b2 = _branch(pg)
            if b1 != branch or b2 != branch:
                straddled = True
                break
            D[:, j] = (fp - fm) / (2 * h)
        if straddled:
            continue
        worst = max(worst, float(np.max(np.abs(K - D)) / max(np.max(np.abs(K)), 1e-30)))
        branches[branch] = branches.get(branch, 0) + 1
        n_ok += 1

    if n_ok < 80:
        fails.append(f"only {n_ok} of 120 random beam elements could be verified; "
                     f"the tangent check is not exercising the element")
    if len(branches) < 8:
        fails.append(
            f"the pile tangent check reached only {len(branches)} of the 8 capacity "
            f"branches ({sorted(branches)}) — a branch it never lands on is a branch "
            f"it proves nothing about")
    if worst > 1e-6:
        fails.append(
            f"the pile beam element's consistent tangent disagrees with a central "
            f"difference of its own internal force by {worst:.2e} relative. Every "
            f"branch is affine in the element displacement — the hinge's plastic "
            f"rotation included — so the tangent is exact and this can only be an "
            f"error in the element.")

    # ---- the same actions the viscoplastic driver reads -------------------
    # The Newton element reads its axial force, shear and end moments off constant
    # ROWS, because it needs their gradients; the viscoplastic path computes the
    # same four quantities in closed form each iteration. They have to be the same
    # numbers, or the two drivers are capping different things.
    from xslope.fem import _pile_element_actions
    worst_act = 0.0
    for trial in range(60):
        n_node = 2 if trial % 2 == 0 else 3
        nd = 3 * n_node
        th = rng.uniform(0, 2 * np.pi)
        L = float(np.exp(rng.uniform(np.log(0.5), np.log(20.0))))
        EA = float(np.exp(rng.uniform(np.log(1e4), np.log(1e8))))
        EI = float(np.exp(rng.uniform(np.log(1e3), np.log(1e7))))
        fd = _beam_case(n_node, th, L, EA, EI)
        pg = _nr_build_piles(fd)[0]
        u = rng.normal(0.0, 1.0, nd) * L * 1e-3
        T_f, V, M1, M2, _ = _pile_element_actions(
            u, np.cos(th), np.sin(th), L, EA, EI,
            fd['K_local_by_pile_elem'][0], n_node)
        got = pg['act'][0] @ u
        for a, b in ((pg['ax'][0] @ u, T_f), (got[0], V), (got[1], M1), (got[2], M2)):
            worst_act = max(worst_act, abs(a - b) / max(abs(b), 1e-30))
    if worst_act > 1e-10:
        fails.append(
            f"the Newton element's action rows disagree with the viscoplastic "
            f"driver's own _pile_element_actions by {worst_act:.2e} relative — the "
            f"two drivers would be capping different quantities")

    # ---- the closed forms -------------------------------------------------
    E, I, A = 2.0e8, 3.0e-4, 1.2e-2
    EI, EA = E * I, E * A
    Ltot, P, n_elem = 6.0, 1000.0, 6

    def _chain(n_node, th):
        Le = Ltot / n_elem
        nd_e = 3 * n_node
        base = _beam_case(n_node, th, Le, EA, EI)
        dof = np.zeros((n_elem, 9), dtype=np.int64)
        for e in range(n_elem):
            row = [3 * e, 3 * e + 1, 3 * e + 2,
                   3 * (e + 1), 3 * (e + 1) + 1, 3 * (e + 1) + 2]
            if n_node == 3:
                m = n_elem + 1 + e
                row += [3 * m, 3 * m + 1, 3 * m + 2]
            dof[e, :nd_e] = row
        fd = dict(base)
        fd['n_pile_elements'] = n_elem
        for key in ('K_global_pile_elems', 'K_local_by_pile_elem'):
            fd[key] = [base[key][0]] * n_elem
        fd['dof_indices_pile'] = dof
        for key, val in (('n_dof_by_pile_elem', nd_e),
                         ('cos_theta_pile', np.cos(th)),
                         ('sin_theta_pile', np.sin(th)),
                         ('elem_length_by_pile_elem', Le),
                         ('EI_by_pile_elem', EI), ('EA_by_pile_elem', EA),
                         ('V_cap_by_pile_elem', np.inf),
                         ('M_cap_by_pile_elem', np.inf)):
            fd[key] = np.full(n_elem, val)
        n_nodes = n_elem + 1 + (n_elem if n_node == 3 else 0)
        return _nr_build_piles(fd), 3 * n_nodes

    def _solve(piles, n_dof, fixed, load):
        K = np.zeros((n_dof, n_dof))
        for pg in piles:
            for j in range(len(pg['idx'])):
                d = pg['dof'][j]
                K[np.ix_(d, d)] += pg['K'][j]
        free = np.array([i for i in range(n_dof) if i not in fixed])
        u = np.zeros(n_dof)
        u[free] = np.linalg.solve(K[np.ix_(free, free)], load[free])
        return u

    for n_node in (2, 3):
        for th in (0.0, np.pi / 2, 0.7):
            c, sn = np.cos(th), np.sin(th)
            piles, n_dof = _chain(n_node, th)
            tip = 3 * n_elem
            load = np.zeros(n_dof)
            load[tip], load[tip + 1] = -sn * P, c * P
            u = _solve(piles, n_dof, {0, 1, 2}, load)
            d_tip = -sn * u[tip] + c * u[tip + 1]
            e_tip = P * Ltot ** 3 / (3.0 * EI)
            e_rot = P * Ltot ** 2 / (2.0 * EI)
            piles2, n_dof2 = _chain(n_node, th)
            mid = 3 * (n_elem // 2)
            load2 = np.zeros(n_dof2)
            load2[mid], load2[mid + 1] = -sn * P, c * P
            u2 = _solve(piles2, n_dof2, {0, 1, 3 * n_elem, 3 * n_elem + 1}, load2)
            d_mid = -sn * u2[mid] + c * u2[mid + 1]
            e_mid = P * Ltot ** 3 / (48.0 * EI)
            for got, want, what in ((d_tip, e_tip, "cantilever tip deflection"),
                                    (u[tip + 2], e_rot, "cantilever tip rotation"),
                                    (d_mid, e_mid, "simply supported deflection")):
                rel = abs(got - want) / abs(want)
                if rel > 1e-8:
                    fails.append(
                        f"the {n_node}-node beam at {np.degrees(th):.0f} deg gets "
                        f"{what} {got:.8e} against the closed form {want:.8e} "
                        f"(relative {rel:.2e})")
    return fails


def check_piles():
    """The Newton driver carries pile beam elements. Four things beyond the element.

    **The displacement bound is a length.** A pile node carries a rotation as its
    third degree of freedom, and the bound the verdict is read on is ``max|u|``. It
    reads the TRANSLATIONAL degrees of freedom only — which is what the viscoplastic
    driver's own displacement-limit check reads — and the index set is asserted to be
    exactly the non-rotational ones. The behavioural half runs on a model meshed
    finely enough along the pile that its rotations exceed its displacements, so a
    bound that let the rotations in would read a different number: the trial has to
    stand under the translational reading and be refused under the raw one.

    **Fixity is already in force.** A held head rotation or a fixed tip is a
    constraint on ``free_dofs``, which this path reads from the prepared model, so
    nothing about it is built here — and that is asserted rather than assumed, on the
    model whose tip is fixed: the held degrees of freedom must be absent from
    ``free_dofs`` and exactly zero in the solution.

    **The capacity is enforced.** With the shear capacity tightened until it binds,
    the shear the element DELIVERS — read out of the residual the driver assembles,
    not out of the reported array — must be at the cap and not above it, and the
    state must differ from the uncapped one. Dropping the correction from the
    residual leaves the reported array clipped and the physics uncapped, which is
    the failure mode this leg exists for.

    **The moment capacity is a plastic hinge**, and the leg is written on what a
    hinge does rather than on what it reports. With ``M_cap`` tightened until it
    binds: an end moment READ OUT OF THE ASSEMBLED RESIDUAL — the rotational
    component of the element's own internal force, which is the moment the
    equilibrium carries — must sit at the capacity and not above it; the hinged ends
    must carry a plastic rotation; at most ONE element end per pile node may carry
    it, since a hinge is one release at one section; and the slope must MOVE, and
    move FARTHER, because a member that can deliver less moment cannot hold the soil
    back better. The last two are what a nodal moment applied at the shared
    rotational degree of freedom cannot do: it is equal and opposite on the two
    elements meeting there and cancels, so the state does not move at all while the
    reported moments read as capped. That is the form this round removed, and it is
    the mutation.

    **The factor of safety.** The FEM pile sample at its own tagged mesh stands at a
    strength below its published lock and fails at one above it, on BOTH drivers —
    the two agree trial for trial over the whole bisection there, so the bracket
    holds the same fact for four solves instead of two whole runs.
    """
    from xslope.fem import (_nr_build_piles, _nr_pile_force,
                            _nr_translational_dofs, _prepare_fem_model)
    fails = []
    fem_data = _pile_fem_data()
    n_pile = int(fem_data.get('n_pile_elements', 0))
    if n_pile == 0:
        return ["the pile test model carries no pile elements — the model, not the "
                "driver, is broken"]

    # ---- the translational index -------------------------------------------
    n_dof = int(fem_data['dof_offset'][len(fem_data['nodes'])])
    tidx = _nr_translational_dofs(fem_data, n_dof)
    is_pile = np.asarray(fem_data['is_pile_node'], dtype=bool)
    rot = np.asarray(fem_data['dof_offset'][:len(is_pile)])[is_pile] + 2
    if tidx is None or set(np.setdiff1d(np.arange(n_dof), tidx)) != set(rot.tolist()):
        fails.append(
            "the Newton driver's displacement bound does not read the translational "
            "degrees of freedom only: the set it excludes is not exactly the pile "
            "nodes' rotations, so the bound is comparing a length against a radian")

    # ---- fixity ------------------------------------------------------------
    prep = _prepare_fem_model(fem_data)
    free = set(int(d) for d in prep['free_dofs'])
    held = []
    for nodes_, rot_held, trans_held, end in (
            (fem_data['pile_tip_nodes'], fem_data['pile_tip_fixed'],
             fem_data['pile_tip_pinned'], 'tip'),
            (fem_data['pile_head_nodes'], fem_data['pile_head_fixed'],
             fem_data['pile_head_pinned'], 'head')):
        for k in range(len(nodes_)):
            base = int(fem_data['dof_offset'][int(nodes_[k])])
            if bool(rot_held[k]):
                held.append((base + 2, f"{end} rotation"))
            if bool(trans_held[k]):
                held += [(base, f"{end} x"), (base + 1, f"{end} y")]
    if not held:
        fails.append("the pile test model holds neither end, so the fixity leg is "
                     "not exercising anything")
    for dof, what in held:
        if dof in free:
            fails.append(f"the pile {what} is declared held and its degree of "
                         f"freedom is still free on the Newton path")

    # ---- the capacity, and what it delivers --------------------------------
    F_CAP = 1.2
    sol_free = solve_fem(fem_data, F=F_CAP, fem_solver='newton',
                         max_disp_factor=None)
    for dof, what in held:
        if abs(float(sol_free['displacements'][dof])) > 1e-12:
            fails.append(f"the pile {what} is held but moved "
                         f"{sol_free['displacements'][dof]:.3e} in the solution")
    V_free = float(np.max(np.abs(sol_free['forces_pile_lateral'])))
    u_free = _nr_translational_max(sol_free['displacements'], tidx)
    cap = 0.25 * V_free
    saved = fem_data['V_cap_by_pile_elem']
    try:
        fem_data['V_cap_by_pile_elem'] = np.full(n_pile, cap)
        sol_cap = solve_fem(fem_data, F=F_CAP, fem_solver='newton',
                            max_disp_factor=None)
    finally:
        fem_data['V_cap_by_pile_elem'] = saved
    u_cap = _nr_translational_max(sol_cap['displacements'], tidx)
    if int(np.count_nonzero(sol_cap['yielded_pile_V'])) == 0:
        fails.append(f"a shear cap at a quarter of the free shear ({cap:.1f} against "
                     f"{V_free:.1f}) binds on no element, so the capacity leg is not "
                     f"exercising the cap")
    if float(np.max(np.abs(sol_cap['forces_pile_lateral']))) > cap * (1 + 1e-9):
        fails.append(
            f"a pile element reports a shear of "
            f"{float(np.max(np.abs(sol_cap['forces_pile_lateral']))):.4f} above its "
            f"cap of {cap:.4f}")
    if np.allclose(sol_cap['displacements'], sol_free['displacements']):
        fails.append(
            "tightening the pile shear capacity to a quarter of the free shear left "
            "the displacement field unchanged, so the correction is not reaching the "
            "residual — the reported force would be clipped while the physics stayed "
            "uncapped, which is the failure mode this leg exists for")
    elif not u_cap > u_free:
        fails.append(
            f"capping the pile shear at {cap:.1f} where the uncapped pile carries "
            f"{V_free:.1f} made the slope move LESS, not more: max|u| went from "
            f"{u_free:.6g} to {u_cap:.6g}. A member that can deliver less force "
            f"cannot hold the soil back better — that is the signature of a "
            f"correction applied with the wrong sign (SPIKE.md, \"PILES\").")

    # ---- the moment capacity, which is a hinge -----------------------------
    M_free = float(np.max(np.abs(sol_free['forces_pile_moment'])))
    m_cap = 0.5 * M_free
    saved_M = fem_data['M_cap_by_pile_elem']
    try:
        fem_data['M_cap_by_pile_elem'] = np.full(n_pile, m_cap)
        sol_M = solve_fem(fem_data, F=F_CAP, fem_solver='newton',
                          max_disp_factor=None)
    finally:
        fem_data['M_cap_by_pile_elem'] = saved_M
    u_M = _nr_translational_max(sol_M['displacements'], tidx)
    p_rot = np.abs(np.asarray(sol_M.get('pile_plastic_rotation',
                                        np.zeros((n_pile, 2)))))
    if int(np.count_nonzero(sol_M['yielded_pile_M'])) == 0:
        fails.append(f"a moment cap at half the free demand ({m_cap:.1f} against "
                     f"{M_free:.1f}) binds on no element, so the moment leg is not "
                     f"exercising the cap")
    if not p_rot.any():
        fails.append(
            f"a moment cap at half the free demand binds but no element reports a "
            f"plastic rotation. A moment capacity is a RELEASE of rotational "
            f"continuity; a correction applied at the rotational degree of freedom "
            f"the two adjacent beam elements SHARE is equal and opposite there and "
            f"enforces nothing (SPIKE.md, \"THE PILE HINGE\").")
    # The moment the EQUILIBRIUM carries, read out of the internal force the driver
    # assembles rather than out of the reported array. A form that only clips what it
    # reports leaves this above the capacity.
    from xslope.fem import _nr_build_piles, _nr_pile_force
    fd_M = dict(fem_data)
    fd_M['M_cap_by_pile_elem'] = np.full(n_pile, m_cap)
    pgs = _nr_build_piles(fd_M)
    worst_del, owners = 0.0, {}
    node_of = fem_data.get('pile_elem_nodes', None)
    for pg in pgs:
        fint = _nr_pile_force(pg, np.asarray(sol_M['displacements']), False)
        for j, e in enumerate(pg['idx']):
            for end, slot in ((0, 2), (1, 5)):
                delivered = float(fint[j, slot])
                reported = float(sol_M['forces_pile_moment'][e][end])
                worst_del = max(worst_del,
                                abs(delivered - reported) / max(abs(reported), 1e-9))
                if abs(delivered) > m_cap * (1 + 1e-9):
                    fails.append(
                        f"pile element {int(e)} end {end + 1} DELIVERS a moment of "
                        f"{delivered:.6f} where its capacity is {m_cap:.6f} — the "
                        f"reported array is capped and the equilibrium is not")
            if node_of is not None and len(node_of) > int(e):
                for end in (0, 1):
                    if pg['_p_rot'][j, end]:
                        nd_ = int(node_of[int(e)][end])
                        owners.setdefault(nd_, []).append((int(e), end))
    if worst_del > 1e-9:
        fails.append(
            f"the pile end moment the Newton residual carries differs from the one "
            f"reported by {worst_del:.2e} relative — the report is not the physics")
    for nd_, who in owners.items():
        if len(who) > 1:
            fails.append(
                f"pile node {nd_} carries a plastic release on {len(who)} element "
                f"ends ({who}). A hinge is ONE release at one section; releasing "
                f"both leaves the node's rotation seeing equal and opposite "
                f"capacities whatever it does, undetermined by the beam.")
    if np.allclose(sol_M['displacements'], sol_free['displacements']):
        fails.append(
            f"tightening the pile MOMENT capacity to {m_cap:.1f} where the uncapped "
            f"pile carries {M_free:.1f} left the displacement field unchanged. That "
            f"is the signature of the cancelling nodal-moment form: the correction "
            f"is applied at the rotational degree of freedom the two adjacent "
            f"elements share, is equal and opposite there, and enforces nothing.")
    elif not u_M > u_free:
        fails.append(
            f"capping the pile moment at {m_cap:.1f} where the uncapped pile "
            f"carries {M_free:.1f} made the slope move LESS, not more: max|u| went "
            f"from {u_free:.6g} to {u_M:.6g}. A member that can deliver less moment "
            f"cannot hold the soil back better.")

    # ---- the factor of safety, both drivers --------------------------------
    # The FEM pile sample brackets its published lock on BOTH drivers: it stands at
    # a strength below the lock and fails at one above it, so the limit is inside
    # the bracket and the lock is inside it too. Two solves a driver rather than two
    # whole bisections, and the same statement.
    piles_fd = _pile_fem_data(path=PILES_MODEL, target_size=PILES_SIZE)
    if not (PILES_STANDS < PILES_LOCKED_FS < PILES_FAILS):
        fails.append(f"the pile bracket [{PILES_STANDS}, {PILES_FAILS}] does not "
                     f"contain the locked {PILES_LOCKED_FS}")
    for solver in ('viscoplastic', 'newton'):
        for F, want in ((PILES_STANDS, True), (PILES_FAILS, False)):
            sol = solve_fem(piles_fd, F=F, fem_solver=solver, max_iterations=16000,
                            max_disp_factor=None)
            got = bool(sol.get('stable', sol['converged']))
            if got != want:
                fails.append(
                    f"{solver}: the FEM pile sample at F = {F} came back "
                    f"{'CONVERGED' if got else 'FAILED'} where it must "
                    f"{'stand' if want else 'fail'} — the two trials bracket the "
                    f"published {PILES_LOCKED_FS} and both drivers agree on them "
                    f"(exit {sol.get('exit_reason')})")
    return fails


# Post-peak softening (SPIKE.md, "POST-PEAK SOFTENING"). The corpus cannot supply a
# specimen on which the latch fires on the NEWTON path: on both published reinforced
# models the Newton equilibrium at the deciding strengths leaves every softenable bar
# at about 0.97 of its peak, so nothing drops. The specimen is therefore constructed
# out of the same model, the same soil and the same six lines, with the capacities
# drawn down until the demand has to cross them — the same construction the
# reinforcement round used for its half-capacity row. Measured on this checkout:
# nothing softens at F = 1.1 on either driver, and at F = 1.2 the Newton latch drops
# seven bar elements where the viscoplastic one drops four.
SOFT_SCALE, SOFT_RES = 0.35, 0.5
SOFT_F_QUIET, SOFT_F_DROPS = 1.1, 1.2


def _soft_fem_data():
    """The constructed specimen: the geogrid model with capacities drawn down."""
    fd = _reinf_fem_data(target_size=4.0, soften=True)
    tal = np.asarray(fd['t_allow_by_1d_elem'], float) * SOFT_SCALE
    fd['t_allow_by_1d_elem'] = tal
    fd['t_res_by_1d_elem'] = np.where(
        np.isfinite(fd['t_res_by_1d_elem']), tal * SOFT_RES, np.nan)
    return fd


def check_softening():
    """The post-peak drop on the Newton path: it fires, it bites, and it is inert
    where nothing crosses.

    Three things, and the third is what makes the first two a measurement rather
    than a description.

    **It fires.** On a model whose bars are asked for more than they can hold, the
    Newton latch must drop some of them and no dropped bar may report a force above
    its residual. Remove the latch and this leg is the one that goes.

    **The drop actually lowers the cap.** A latch that records the set without moving
    the working capacity would report softening and enforce the peak, which is the
    same failure mode the pile round's moment cap had. So the reported force of every
    softened bar is checked against its RESIDUAL and not against its peak.

    **It is inert where nothing crosses.** At a strength where no softenable bar
    reaches its peak, the softened set must be empty on BOTH drivers, and the Newton
    solve must be bit-identical to the same model with no residual declared at all —
    same displacement field, same iteration count. That is what says the feature
    costs nothing on the models that do not use it.

    The driver no longer REFUSES softening, so the guard check cannot assert one; the
    first leg here is what replaces it.
    """
    fails = []
    fd = _soft_fem_data()
    tal = np.asarray(fd['t_allow_by_1d_elem'], float)
    tres = np.asarray(fd['t_res_by_1d_elem'], float)
    can = np.isfinite(tres) & (tres < tal - 1e-12)
    if not can.any():
        return ["the softening specimen declares no softenable bar, so nothing "
                "here is a softening measurement"]

    # --- it fires, and the drop lowers the cap -----------------------------
    try:
        sol = solve_fem(fd, F=SOFT_F_DROPS, fem_solver='newton',
                        max_disp_factor=None, max_iterations=16000)
    except NotImplementedError as exc:
        return [f"the Newton driver still REFUSES post-peak softening: {exc}"]
    soft = np.asarray(sol['softened_1d_elements'], bool)
    T = np.asarray(sol['forces_1d'], float)
    if not soft.any():
        fails.append(
            f"at F = {SOFT_F_DROPS} on a model whose bars are asked for more than "
            f"they can hold, the Newton latch dropped NO bar element. The latch reads "
            f"the uncapped demand against t_allow on a state in equilibrium with full "
            f"gravity, exactly as the viscoplastic driver does; if it never fires, "
            f"the post-peak law is not being applied at all.")
    else:
        over = np.flatnonzero(soft & (T > tres + 1e-6))
        if over.size:
            fails.append(
                f"{over.size} softened bar element(s) report a force above their "
                f"RESIDUAL — the largest is {float(T[over].max()):.3f} against a "
                f"residual of {float(tres[over].min()):.3f}. The latch recorded the "
                f"drop without lowering the working capacity, so the equilibrium is "
                f"still carrying the peak.")

    # --- inert where nothing crosses ---------------------------------------
    quiet = solve_fem(fd, F=SOFT_F_QUIET, fem_solver='newton',
                      max_disp_factor=None, max_iterations=16000)
    if np.any(np.asarray(quiet['softened_1d_elements'], bool)):
        fails.append(f"a bar softened at F = {SOFT_F_QUIET}, where no softenable bar "
                     f"reaches its peak — the trigger is not the demand")
    quiet_vp = solve_fem(fd, F=SOFT_F_QUIET, fem_solver='viscoplastic',
                         max_disp_factor=None, max_iterations=16000)
    if np.any(np.asarray(quiet_vp['softened_1d_elements'], bool)):
        fails.append(f"the VISCOPLASTIC driver softened a bar at F = {SOFT_F_QUIET} "
                     f"where the Newton driver did not, so the two are not reading "
                     f"the same trigger on the same model")
    fd_off = dict(fd)
    fd_off['t_res_by_1d_elem'] = np.full(len(tres), np.nan)
    off = solve_fem(fd_off, F=SOFT_F_QUIET, fem_solver='newton',
                    max_disp_factor=None, max_iterations=16000)
    if not np.array_equal(np.asarray(quiet['displacements']),
                          np.asarray(off['displacements'])):
        fails.append(
            "at a strength where nothing softens, the model that DECLARES a residual "
            "and the model that declares none do not return the same displacement "
            "field. The post-peak path must cost nothing where it does not fire.")
    if quiet['iterations'] != off['iterations']:
        fails.append(
            f"at a strength where nothing softens, declaring a residual changed the "
            f"iteration count from {off['iterations']} to {quiet['iterations']}")
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


# The rescue chain's cost policy (SPIKE.md, "THE COST OF THE RESCUE"). Phase 0
# attributed half the Newton driver's wall time to a chain that confirms failures the
# cold attempt had already reached, and the policy budgets that chain by how far a
# trial sits above the highest strength the search has shown to stand.
#
# THE POLICY SHIPS OFF, and these two numbers are NOT a safe setting. They are the
# thresholds a 46-row slice supported -- on that slice no rung converted a verdict
# more than half the initial bracket above a standing bound and the adaptive rung
# never more than a quarter -- and the 191-row corpus falsified them: on
# SSRM-TORGGLER the adaptive rung carries a trial at exactly half a bracket above the
# bound, which a 0.30 threshold takes away, and the factor of safety moves by a whole
# bisection cell. What this check holds is the MECHANISM: that the policy is off by
# default, that it is wired to the chain, and that on this fixture it changes what
# the driver spends and nothing else. It is not evidence that any threshold is safe.
COST_ADAPTIVE_FRAC = 0.30
COST_NONE_FRAC = 0.75


def _bisect_newton(fem_data, **over):
    """The check's own Newton bisection, with module globals optionally overridden."""
    saved = {k: getattr(_fem, k) for k in over}
    for k, v in over.items():
        setattr(_fem, k, v)
    try:
        return solve_ssrm(fem_data, F_min=F_MIN, F_max=F_MAX, tolerance=TOLERANCE,
                          max_iterations=MAX_ITER, fem_solver='newton',
                          capture_failure_state=False)
    finally:
        for k, v in saved.items():
            setattr(_fem, k, v)


def _verdict_fingerprint(res):
    """Everything a cost policy is forbidden to move."""
    return (round(float(res['FS']), 12),
            tuple(round(float(x), 12) for x in res['final_interval']),
            tuple((round(float(t['F']), 12), bool(t['stable']),
                   str(t['exit_reason'])) for t in res['trials']))


def _predictor_work(res):
    trials = res['trials']
    return (sum(int(t.get('nr_predictor_iterations', 0)) for t in trials),
            sum(1 for t in trials if t.get('nr_predictor_iterations')))


def check_rescue_cost(fem_data, base):
    """The chain's cost policy is off, is wired, and on this fixture moves no verdict.

    The viscoplastic predictor is what makes the Newton driver decide every trial --
    across the 191-row corpus it left none inconclusive where the viscoplastic driver
    left 37 -- and it is also two thirds of the driver's constitutive work, because
    it runs on every trial that dies at the load-step floor and half of every
    bisection does. Budgeting it by distance from the standing bound is the cheapest
    way to keep both, and the corpus says no threshold measured so far is safe (see
    COST_ADAPTIVE_FRAC above). So this locks the machinery rather than a setting.

    Four assertions, and the last two are what make the first two mean anything:

      * the policy SHIPS OFF, so the default Newton driver is the one every measured
        value in SPIKE.md is defined on;
      * this fixture still exercises the chain -- a mesh on which the predictor never
        fired would pass everything below without testing any of it;
      * turning the policy on changes what the driver SPENDS and nothing else here:
        the same factor of safety, the same final interval, the same verdict for
        every trial at the same strength, on strictly less predictor work;
      * and the knob reaches the chain: drive the threshold to zero and no trial
        above a standing bound may be charged a predictor iteration. Without that, a
        rule quietly doing nothing would pass the identity test perfectly.
    """
    fails = []

    for name, want in (('_NR_RESCUE_ADAPTIVE_FRAC', None),
                       ('_NR_RESCUE_NONE_FRAC', None),
                       ('_NR_SEED_MEMORY', False),
                       ('_NR_COLD_CHEAP', False)):
        got = getattr(_fem, name, '<not defined>')
        if got != want:
            fails.append(
                f"the rescue cost policy does not ship off: {name} is {got!r}, not "
                f"{want!r}. Every Newton factor of safety in SPIKE.md is measured "
                f"with the whole chain running on every trial that reaches it.")

    p0, fired0 = _predictor_work(base)
    if fired0 < 2:
        return fails + [
            f"the predictor fired on {fired0} of {len(base['trials'])} trials on this "
            f"fixture, so there is no rescue chain here to budget and nothing below "
            f"this line tests one"]

    cheap = _bisect_newton(fem_data,
                           _NR_RESCUE_ADAPTIVE_FRAC=COST_ADAPTIVE_FRAC,
                           _NR_RESCUE_NONE_FRAC=COST_NONE_FRAC)
    if _verdict_fingerprint(cheap) != _verdict_fingerprint(base):
        fails.append(
            f"budgeting the rescue chain by distance moved a verdict: FS "
            f"{base['FS']:.6f} on {base['final_interval']} became {cheap['FS']:.6f} "
            f"on {cheap['final_interval']}. A cost policy may change what the driver "
            f"spends and nothing else.")
    p1, _ = _predictor_work(cheap)
    if p1 >= p0:
        fails.append(
            f"the distance rule saved nothing: {p1} predictor iterations against "
            f"{p0} with the whole chain. The bracket-establishment trial alone sits a "
            f"full bracket above the standing bound, so it should never reach the "
            f"adaptive rung.")

    silent = _bisect_newton(fem_data, _NR_RESCUE_NONE_FRAC=0.0)
    charged = [t['F'] for t in silent['trials']
               if t.get('nr_predictor_iterations') and t.get('bracket')
               and float(t['F']) > float(t['bracket'][0])]
    if charged:
        fails.append(
            "the rung budget does not reach the chain: with the distance threshold at "
            "zero the predictor still ran on the trial(s) at F = "
            + ", ".join(f"{f:.4f}" for f in charged)
            + ". The identity result above cannot be read as evidence while the knob "
              "it turns is wired to nothing.")
    return fails


# The factorization is the Newton driver's single largest expense -- 35% of a
# trial's wall time on the 13-row representative set, against 20% for the whole
# line search (SPIKE.md, "THE FACTORIZATION"). Two things about it are held here.
#
# The FIRST is that the column ordering may be cached. splu's COLAMD ordering is a
# function of the sparsity pattern alone, a trial's tangent pattern is built once,
# and re-forming into that fixed pattern therefore re-derives the same permutation
# hundreds of times. _nr_factorize_tangent keeps it and hands SuperLU the already
# permuted matrix under permc_spec='NATURAL', which is bit-for-bit the same
# factorization -- not close, the same. That claim is the reason the reuse is
# allowed to ship, so it is asserted directly on a tangent taken out of a real
# trial rather than inferred from a factor of safety that matched.
#
# The SECOND is the REFRESH RULE, which is the lever this round measured and did not
# adopt. _nr_equilibrate holds its factorization while the residual keeps falling
# fourfold and re-forms when it does not, and on the representative set that rule
# fires on 95% to 97% of iterations: a line search that accepts a fractional step
# does not cut the residual by four, so the tangent is re-formed on nearly every
# iteration and the reuse is worth almost nothing. The count is locked so that a
# change to the rule shows up as a failing assertion rather than as a wall time
# nobody re-measures.
REFORM_FRACTION_MIN = 0.90   # the representative set measured 0.96, with room below
# A strength past this fixture's limit of about 1.39, so the trial walks its load
# increment down to the floor and the iteration count is large enough for the refresh
# rule to be measuring something. A converged trial at 1.35 takes 17 iterations.
FACTOR_PROBE = 1.45


# The Newton corrector (SPIKE.md, "THE CORRECTOR"). Every number below was
# measured on this checkout on 2026-09-03, on the same coarse tri6 fixture the rest
# of the file uses (387 elements), and each one is a property of the mechanism
# rather than of the benchmark:
#
#   * at F = 1.39 the viscoplastic loop cannot decide the trial inside ANY budget it
#     is given here — 200 iterations or 4,000, it comes back 'inconclusive' with the
#     residual still falling — and the corrector certifies the same trial from a
#     300-pass seed in 31 Newton iterations;
#   * at F = 1.4125 no corrector certifies anything, and the trial must come back
#     from the 'auto' path bit for bit as the viscoplastic loop left it;
#   * at F = 1.30 the loop settles on its own in 180 passes — before the first
#     checkpoint at 300 — so no attempt is made at all and the iteration count is
#     the shipped one.
CORR_F = 1.39
CORR_SMALL_BUDGET = 200    # too small for the ladder: only the rule exit can fire
CORR_BUDGET = 4000
CORR_REFUSED_F = 1.4125    # ... a strength no corrector certifies
CORR_SETTLES_F = 1.30      # ... and one the loop settles at before any checkpoint
CORR_SETTLES_ITERS = 180   # the shipped loop's iteration count there
CORR_LADDER_MAX_ITERS = 1000   # a certified trial ends far inside its budget


def check_corrector(fem_data):
    """Newton finishes the trials the viscoplastic loop cannot, and overrules none.

    The default driver is the viscoplastic iteration with a bounded Newton
    corrector called from inside it. Four properties make that safe, and this
    holds all four:

      * a trial the loop ends on a RULE — here the iteration ceiling, with the
        residual still coming down — is handed to the corrector, and a corrector
        that reaches equilibrium ends it as standing;
      * the checkpoint LADDER fires too, so a trial that would have spent its whole
        budget is finished from a 300-pass seed instead;
      * a certified trial carries its evidence, and the three gates are re-read here
        off the returned state rather than taken from the record: force balance to
        the Dawson tolerance at full gravity, a yield violation at or below
        `_CORRECTOR_YIELD_TOL` of the local strength, and a displacement inside the
        bound;
      * a corrector REFUSAL is not a verdict. Where no attempt is certified the
        trial comes back exactly as the viscoplastic loop left it, down to the
        displacement field, and `fem_solver='viscoplastic'` makes no attempt at all;
      * the LADDER IS TRIMMED. `_CORRECTOR_RULE_EXITS` names the loop exits worth an
        attempt, an exit outside it gets none and reports exactly what the shipped
        loop reports, and the rule route still works when an exit is put back;
      * and a CONVERGED viscoplastic state must be ADMISSIBLE. A state that settles
        in force and fails the yield reading is handed to the corrector; where the
        corrector certifies an admissible field the trial stands on that, and where
        it refuses too the trial is FAILED — `converged` False, `stable` False under
        every criterion, `exit_reason` 'yield_gate'.

    The refusal property is the safety property the corrector rests on: it can only
    ever convert a rule-decided refusal into a certified stand, so the corrector
    cannot read lower than the driver it wraps. The yield gate is the one thing here
    that CAN read lower, and it is meant to: a slope stands on an admissible field or
    it does not stand.
    """
    fails = []

    def _solve(F, solver, budget, **over):
        return solve_fem(fem_data, F=F, max_disp_factor=None, force_tol=1e-3,
                         max_iterations=budget, max_iterations_ceiling=budget,
                         early_failure=True, fem_solver=solver, **over)

    def _fingerprint(sol):
        u = np.ascontiguousarray(np.asarray(sol['displacements'], float))
        return (bool(sol['converged']), int(sol['iterations']),
                str(sol['exit_reason']), float(sol['unbalanced_force_ratio']),
                u.tobytes())

    def _gates(sol, where):
        """Re-read the three gates off the state, not off the evidence record."""
        out = []
        c = sol.get('corrector') or {}
        if c.get('driver_of_record') != 'corrector':
            out.append(f"{where}: the trial is certified but carries no "
                       f"driver_of_record: {c!r}")
            return out
        oob = float(sol.get('unbalanced_force_ratio', float('inf')))
        if not oob < 1e-3:
            out.append(f"{where}: certified with the out-of-balance at {oob:.3e}, "
                       f"which is not inside the force tolerance 1e-3")
        yv = sol.get('max_yield_violation')
        if yv is None or float(yv) > _fem._CORRECTOR_YIELD_TOL:
            out.append(f"{where}: certified with a yield violation of {yv!r} of the "
                       f"local strength, above the gate "
                       f"{_fem._CORRECTOR_YIELD_TOL:g}")
        dd, dl = c.get('max_disp_deep'), c.get('disp_limit')
        if dl is not None and (dd is None or float(dd) > float(dl)):
            out.append(f"{where}: certified with a displacement of {dd!r} against "
                       f"the bound {dl!r}")
        if int(c.get('vp_iterations', -1)) + int(c.get('nr_iterations', -1)) \
                != int(sol.get('iterations', 0)):
            out.append(f"{where}: the trial is charged {sol.get('iterations')!r} "
                       f"iterations, which is not the viscoplastic passes "
                       f"{c.get('vp_iterations')!r} plus the corrector's "
                       f"{c.get('nr_iterations')!r} — the seed's work is not being "
                       f"charged to the trial that used it")
        return out

    # --- the rule exit: a ceiling too small for the ladder to reach ------------
    vp_small = _solve(CORR_F, 'viscoplastic', CORR_SMALL_BUDGET)
    if vp_small['converged'] or vp_small['exit_reason'] != 'inconclusive':
        fails.append(
            f"the fixture has moved: at F = {CORR_F} on a {CORR_SMALL_BUDGET}-"
            f"iteration ceiling the viscoplastic loop is meant to come back "
            f"undecided, and it came back {vp_small['exit_reason']!r} "
            f"(converged={vp_small['converged']})")
    # `_CORRECTOR_RULE_EXITS` is the whole switch, and it is locked in both
    # directions. With the ceiling exit taken OUT the corrector must not be asked
    # there and the trial must come back exactly as the shipped loop leaves it; the
    # measurement behind keeping every exit is in SPIKE.md, "THE LADDER AND THE
    # YIELD GATE".
    _saved_exits = _fem._CORRECTOR_RULE_EXITS
    try:
        _fem._CORRECTOR_RULE_EXITS = tuple(
            x for x in _saved_exits if x != 'inconclusive')
        trimmed = _solve(CORR_F, 'auto', CORR_SMALL_BUDGET)
    finally:
        _fem._CORRECTOR_RULE_EXITS = _saved_exits
    if trimmed.get('corrector_attempts') or ():
        fails.append(
            f"with the ceiling exit out of _CORRECTOR_RULE_EXITS the corrector was "
            f"asked there anyway: {len(trimmed['corrector_attempts'])} attempt(s)")
    if _fingerprint(trimmed) != _fingerprint(vp_small):
        fails.append(
            f"with the ceiling exit trimmed out, the default path at F = {CORR_F} "
            f"reports {trimmed['exit_reason']!r} in {trimmed['iterations']} "
            f"iterations against the viscoplastic loop's {vp_small['exit_reason']!r} "
            f"in {vp_small['iterations']}")
    # ... and with it IN — the shipped tuple — the rule route decides the trial.
    if 'inconclusive' not in _fem._CORRECTOR_RULE_EXITS:
        fails.append(
            f"the shipped _CORRECTOR_RULE_EXITS is {_fem._CORRECTOR_RULE_EXITS!r} "
            f"and does not carry the ceiling exit, so the rule route below is not "
            f"the shipped one")
    auto_small = _solve(CORR_F, 'auto', CORR_SMALL_BUDGET)
    if not auto_small['converged']:
        fails.append(
            f"at F = {CORR_F} the corrector did not decide a trial the "
            f"viscoplastic loop left at its ceiling: {auto_small['exit_reason']!r}, "
            f"{len(auto_small.get('corrector_attempts') or [])} attempt(s)")
    else:
        fails += _gates(auto_small, f"F = {CORR_F} at the rule exit")
        at = (auto_small.get('corrector') or {}).get('checkpoint')
        if at != 'rule:inconclusive':
            fails.append(
                f"at F = {CORR_F} with a {CORR_SMALL_BUDGET}-iteration ceiling the "
                f"corrector is meant to be reached by the RULE exit, and the "
                f"evidence says {at!r}")

    # --- the checkpoint ladder ------------------------------------------------
    vp_full = _solve(CORR_F, 'viscoplastic', CORR_BUDGET)
    if vp_full['converged']:
        fails.append(
            f"the fixture has moved: at F = {CORR_F} the viscoplastic loop is meant "
            f"to spend its whole {CORR_BUDGET}-iteration budget undecided, and it "
            f"converged in {vp_full['iterations']}")
    auto_full = _solve(CORR_F, 'auto', CORR_BUDGET)
    if not auto_full['converged']:
        fails.append(
            f"at F = {CORR_F} the corrector did not decide the trial from any "
            f"checkpoint: {auto_full['exit_reason']!r}")
    else:
        fails += _gates(auto_full, f"F = {CORR_F} on the ladder")
        at = (auto_full.get('corrector') or {}).get('checkpoint')
        if at != f"vp{_fem._CORRECTOR_CHECKPOINTS[0]}":
            fails.append(
                f"at F = {CORR_F} the trial is meant to be finished from the FIRST "
                f"checkpoint seed, and the evidence says {at!r}")
        if auto_full['iterations'] > CORR_LADDER_MAX_ITERS:
            fails.append(
                f"at F = {CORR_F} the certified trial cost "
                f"{auto_full['iterations']} iterations against the viscoplastic "
                f"loop's {vp_full['iterations']}: the corrector is meant to finish "
                f"it from a short seed, not to ride the whole budget")

    # --- a refusal is not a verdict -------------------------------------------
    vp_ref = _solve(CORR_REFUSED_F, 'viscoplastic', CORR_BUDGET)
    auto_ref = _solve(CORR_REFUSED_F, 'auto', CORR_BUDGET)
    n_att = len(auto_ref.get('corrector_attempts') or [])
    if n_att == 0:
        fails.append(
            f"the fixture has moved: at F = {CORR_REFUSED_F} the corrector is meant "
            f"to be ASKED and to refuse, and it was never asked — the identity "
            f"below would then be measuring nothing")
    if any(a.get('certified') for a in (auto_ref.get('corrector_attempts') or ())):
        fails.append(
            f"the fixture has moved: at F = {CORR_REFUSED_F} an attempt was "
            f"certified, so this trial no longer measures what a refusal does")
    elif _fingerprint(vp_ref) != _fingerprint(auto_ref):
        fails.append(
            f"at F = {CORR_REFUSED_F} the {n_att} refused corrector attempt(s) "
            f"changed the trial: the viscoplastic loop reports "
            f"{vp_ref['exit_reason']!r} in {vp_ref['iterations']} iterations and "
            f"the default path {auto_ref['exit_reason']!r} in "
            f"{auto_ref['iterations']} — a refusal decides nothing and must leave "
            f"the trial exactly as it found it")

    # ... and the same, forced: with the yield gate made unreachable no state can
    # be certified, so the default path has to be the viscoplastic loop again.
    saved = _fem._CORRECTOR_YIELD_TOL
    try:
        _fem._CORRECTOR_YIELD_TOL = -1.0
        auto_gated = _solve(CORR_F, 'auto', CORR_BUDGET)
    finally:
        _fem._CORRECTOR_YIELD_TOL = saved
    if not (auto_gated.get('corrector_attempts') or ()):
        fails.append(
            "with the yield gate made unreachable the corrector was never asked, "
            "so the refusal path is untested")
    elif _fingerprint(vp_full) != _fingerprint(auto_gated):
        fails.append(
            f"with the yield gate made unreachable at F = {CORR_F} the default path "
            f"reports {auto_gated['exit_reason']!r} in {auto_gated['iterations']} "
            f"iterations against the viscoplastic loop's {vp_full['exit_reason']!r} "
            f"in {vp_full['iterations']}: a state the gate refuses may not end a "
            f"trial")

    # --- the escape hatch, and a trial that never needs a corrector ------------
    for sol, name in ((vp_small, f"F = {CORR_F} at the ceiling"),
                      (vp_full, f"F = {CORR_F} on the full budget"),
                      (vp_ref, f"F = {CORR_REFUSED_F}")):
        if sol.get('corrector') is not None or (sol.get('corrector_attempts') or ()):
            fails.append(
                f"fem_solver='viscoplastic' ran a corrector at {name}: the escape "
                f"hatch is meant to be the loop with nothing wrapped around it")
    settles = _solve(CORR_SETTLES_F, 'auto', CORR_BUDGET)
    if not settles['converged'] or settles['iterations'] != CORR_SETTLES_ITERS:
        fails.append(
            f"at F = {CORR_SETTLES_F} the loop settles before the first checkpoint, "
            f"so the default path must report the shipped "
            f"{CORR_SETTLES_ITERS} iterations, and it reports "
            f"{settles['iterations']} (converged={settles['converged']})")
    if settles.get('corrector_attempts'):
        fails.append(
            f"at F = {CORR_SETTLES_F} the loop converges in "
            f"{settles['iterations']} iterations, before the first checkpoint at "
            f"{_fem._CORRECTOR_CHECKPOINTS[0]}, and a corrector was asked anyway")

    # --- a CONVERGED viscoplastic state must be ADMISSIBLE ---------------------
    # The force gate cannot see a yield violation, so a trial that settles in force
    # is not yet a verdict: it must also pass the invariant-form yield reading
    # (_VP_YIELD_GATE). The gate is exercised by making it unreachable on the trial
    # the loop settles fastest on — the one that needs no corrector at all — so the
    # only thing that can change the answer is the gate itself.
    _saved_vp_gate = _fem._VP_YIELD_GATE
    try:
        _fem._VP_YIELD_GATE = -1.0
        gated = _solve(CORR_SETTLES_F, 'auto', CORR_BUDGET)
        # ... and with the corrector's own gate shut as well, nothing can be
        # certified, so there is no admissible field anywhere and the trial FAILS.
        _saved_c_gate = _fem._CORRECTOR_YIELD_TOL
        try:
            _fem._CORRECTOR_YIELD_TOL = -1.0
            failed = _solve(CORR_SETTLES_F, 'auto', CORR_BUDGET)
            hybrid = _solve(CORR_SETTLES_F, 'auto', CORR_BUDGET,
                            failure_criterion='hybrid')
        finally:
            _fem._CORRECTOR_YIELD_TOL = _saved_c_gate
    finally:
        _fem._VP_YIELD_GATE = _saved_vp_gate
    if not (gated.get('corrector_attempts') or ()):
        fails.append(
            f"with the yield gate unreachable the force-settled state at "
            f"F = {CORR_SETTLES_F} was never handed to the corrector: a state that "
            f"converges in FORCE and fails in YIELD is not a verdict")
    if gated['converged'] and (gated.get('corrector') or {}) \
            .get('driver_of_record') != 'corrector':
        fails.append(
            f"with the yield gate unreachable the trial at F = {CORR_SETTLES_F} "
            f"still came back CONVERGED on the viscoplastic state itself "
            f"({gated['exit_reason']!r}, {gated['iterations']} iterations): an "
            f"inadmissible converged state may not end a trial")
    if gated['converged']:
        fails += _gates(gated, f"F = {CORR_SETTLES_F} through the yield gate")
    for sol, name in ((failed, "the default criterion"),
                      (hybrid, "failure_criterion='hybrid'")):
        if sol['converged'] or sol.get('stable') or not sol.get('gate_failed'):
            fails.append(
                f"with both yield gates unreachable at F = {CORR_SETTLES_F} under "
                f"{name} the trial came back converged={sol['converged']}, "
                f"stable={sol.get('stable')}, gate_failed={sol.get('gate_failed')}, "
                f"exit {sol['exit_reason']!r}: with no admissible field anywhere the "
                f"trial has to be FAILED")
        elif sol['exit_reason'] != 'yield_gate' or sol.get('verdict') != 'FAILED':
            fails.append(
                f"the gate-failed trial at F = {CORR_SETTLES_F} under {name} reports "
                f"exit {sol['exit_reason']!r} / verdict {sol.get('verdict')!r} "
                f"instead of 'yield_gate' / 'FAILED', so nothing downstream can say "
                f"WHY it failed")

    # --- the admissibility reading is on every result, on both drivers ---------
    for sol, name in ((vp_full, 'the viscoplastic loop'), (auto_full, 'the default path')):
        for key in ('max_yield_violation', 'n_yield_above_1pct', 'yield_flagged',
                    'gate_failed'):
            if key not in sol:
                fails.append(f"{name} does not report {key!r}: the force gate "
                             f"cannot see a yield violation, so the reading has to "
                             f"be carried on every result")
    return fails


def check_factorization(fem_data):
    """The cached ordering is the same factorization, and the refresh rule is measured.

    Four assertions, and the last is what makes the third mean anything:

      * the measurement knobs ship as SPIKE.md's numbers were taken -- the refresh
        rule on its residual ratio, the non-symmetric factorization;
      * a cached-ordering factorization of a REAL tangent returns the same vector as
        splu's own call on the same matrix, bit for bit;
      * the driver re-forms its tangent on nearly every Newton iteration on this
        fixture, which is the measurement that says factorization reuse is not
        happening today;
      * and the refresh rule is wired: holding the tangent for three iterations must
        cut the factorization count. Without that the fraction above is a number the
        code could not move.
    """
    from scipy.sparse.linalg import splu as _splu
    fails = []

    for name, want in (('_NR_REFORM_EVERY', None),
                       ('_NR_REFORM_RATIO', 0.25),
                       ('_NR_FACTOR_SYMMETRIC', False),
                       ('_NR_FACTOR_CACHED_ORDER', True)):
        got = getattr(_fem, name, '<not defined>')
        if got != want:
            fails.append(
                f"the factorization controls do not ship as measured: {name} is "
                f"{got!r}, not {want!r}. Every Newton wall time in SPIKE.md is taken "
                f"on these settings.")

    # A tangent out of a real trial, not a synthetic matrix: the claim is about the
    # matrices this driver actually factorizes.
    grabbed = []
    real = _fem._nr_factorize_tangent

    def _spy(K, cache):
        if len(grabbed) < 1:
            grabbed.append(K.copy())
        return real(K, cache)

    _fem._nr_factorize_tangent = _spy
    try:
        solve_fem(fem_data, F=FACTOR_PROBE, fem_solver='newton', max_iterations=MAX_ITER)
    finally:
        _fem._nr_factorize_tangent = real

    if not grabbed:
        fails.append("no tangent was factorized on this fixture, so nothing below "
                     "this line tests the factorization at all")
    else:
        K = grabbed[0].tocsc()
        b = np.arange(1.0, K.shape[0] + 1.0)
        want = _splu(K).solve(b)
        cache = {}
        real(K, cache)                       # primes the cached ordering
        got = real(K, cache).solve(b)        # ... and this one uses it
        if not np.array_equal(want, got):
            d = float(np.max(np.abs(want - got)) / max(float(np.max(np.abs(want))),
                                                       1e-300))
            fails.append(
                f"the cached column ordering is not the same factorization: the "
                f"solution moved by {d:.3e} relative. It is allowed to ship only "
                f"because it is bit-identical, so any difference at all is a defect "
                f"and not a tolerance.")

    # The refresh rule, measured on the fixture and then moved.
    def _counts(**over):
        saved = {k: getattr(_fem, k) for k in over}
        prof = _fem._PROF_ON
        for k, v in over.items():
            setattr(_fem, k, v)
        _fem._PROF_ON = True
        _fem._prof_reset()
        try:
            solve_fem(fem_data, F=FACTOR_PROBE, fem_solver='newton',
                      max_iterations=MAX_ITER)
            p = _fem._prof_read()
        finally:
            _fem._PROF_ON = prof
            _fem._prof_reset()
            for k, v in saved.items():
                setattr(_fem, k, v)
        return (int(p.get('n_nr_factorize', 0)),
                int(p.get('n_nr_const_tan', 0)) + int(p.get('n_nr_const_res', 0)))

    n_fac, n_it = _counts()
    if n_it < 20:
        fails.append(f"the fixture ran {n_it} Newton iterations at F = {FACTOR_PROBE}, too few "
                     f"for the refresh rule below to be measuring anything")
    else:
        frac = n_fac / n_it
        if frac < REFORM_FRACTION_MIN:
            fails.append(
                f"the tangent refresh rule now fires on {frac:.0%} of iterations "
                f"({n_fac} factorizations in {n_it}), not the {REFORM_FRACTION_MIN:.0%} "
                f"or more this driver was measured at. Factorization is 35% of a "
                f"Newton trial, so a change here moves every wall time in SPIKE.md "
                f"and has to be re-measured rather than merely noticed.")
        n_fac3, _ = _counts(_NR_REFORM_EVERY=3)
        if n_fac3 >= n_fac:
            fails.append(
                f"the refresh rule is not wired: holding the tangent for three "
                f"iterations still took {n_fac3} factorizations against {n_fac}. The "
                f"fraction above cannot be read as evidence while the knob that "
                f"sets it reaches nothing.")
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
          "the loading path nor the step control changes a verdict, a bar's "
          "post-peak drop is latched on a converged state and walked down as a "
          "continuation while Hoek-Brown and "
          "power-curve materials are carried per Gauss point on the envelope the "
          "material declares — mixed with Mohr-Coulomb in one mesh, reduced on "
          "their tangent and not on their own constants — "
          "a reinforced bisection lands on its measured strength with every "
          "bar force inside the capacity its embedment develops, a cohesionless "
          "base does not manufacture a failure at a strength the driver's "
          "own converged answer above it proves standing, a Rankine tensile cap "
          "returns to both yield surfaces at once and is reduced with the trial "
          "strength while an inert cap changes nothing at all, the K0 initial "
          "stress is an exact equilibrium on level ground down to its "
          "out-of-plane component and reproduces a locked vendor benchmark on "
          "both drivers, the matric-suction apparent cohesion is the same field "
          "on both drivers and is bounded by its cap and reduced by the trial "
          "strength, the pile beam element reproduces the closed-form beam and its own "
          "tangent to round-off on every capacity branch while the displacement "
          "bound stays a length, a minimum slip depth keeps a failing surficial "
          "skin from ending a trial the deep mechanism carries, and the "
          "environment override cannot swap the driver in silence. The monotonic "
          "ramp reaches the same limit along one warm-started history, reports it "
          "on the bisection's midpoint convention, and never solves past it. The "
          "rescue chain's cost policy ships off, is wired to the chain it budgets, "
          "and where it is on it changes what the driver spends and no verdict. The "
          "tangent's cached column ordering is the same factorization bit for bit, "
          "and the refresh rule that decides how often that factorization is paid "
          "for is measured and wired. The Newton corrector finishes the trials the "
          "viscoplastic loop leaves at its ceiling — from the rule exit and from "
          "the first checkpoint seed alike — certifies each one in force, in yield "
          "and in displacement, charges the seed's work to the trial that used it, "
          "and where no state passes those gates it leaves the trial exactly as the "
          "loop left it, down to the displacement field.")


if __name__ == '__main__':
    main()
