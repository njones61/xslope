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
    failures += check_tension_cutoff(fem_data)
    failures += check_matric_suction()
    failures += check_k0_initial_stress()
    failures += check_pile_element()
    failures += check_piles()
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
      * the power-curve envelope (``option='pow'``), linearized the same way.

    Matric suction used to be the third. It is CARRIED now (see
    :func:`check_matric_suction`), so it is no longer asserted here — the guard
    that refused it is gone, and a check still expecting the refusal would be
    holding the driver to a limitation it no longer has.

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

    cases = [
        ("Hoek-Brown",
         _variant(option='hb', hb_sci=30000.0, hb_gsi=60.0, hb_mi=10.0, hb_d=0.0),
         {}, 'hb_flag_by_elem'),
        ("power-curve",
         _variant(option='pow', pow_a=1.1, pow_b=0.9, pow_c=2.0, pow_d=4.0),
         {}, 'pow_flag_by_elem'),
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
    """Bars are carried. Softening bars are not, and they say so.

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

    # Piles used to be asserted here. They are CARRIED now (see
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

    **The consistent tangent is exact**, because every action the capacity clips is
    LINEAR in the element displacement and every branch of the clip is affine — so
    the analytic element tangent must equal a central difference of the element's
    own internal force to round-off, on every branch. 120 random elements, random
    orientation, length and rigidities, two-node and three-node, with the caps set
    from the actions the drawn displacement actually produces so that all eight
    combinations of (shear, end moment 1, end moment 2) at their cap are reached.
    A probe that straddles a branch boundary is skipped rather than averaged, and
    the histogram is asserted to have reached every branch — a check that never
    lands on the capped branches would prove nothing about them.

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
        caps = np.abs(act) * rng.choice([0.3, 0.8, 1.5, 3.0], size=3)
        caps[rng.random(3) < 0.25] = np.inf
        pg = _nr_build_piles(_beam_case(n_node, th, L, EA, EI, caps[0], caps[1]))[0]
        pg['cap'][0, 2] = caps[2]
        f0 = _nr_pile_force(pg, u, True)
        K = pg['_Ke'][0].copy()
        branch = tuple(bool(b) for b in (np.abs(pg['_s'][0]) > pg['cap'][0]))
        h = 1e-5 * float(np.max(np.abs(u)))
        D = np.zeros((nd, nd))
        straddled = False
        for j in range(nd):
            up, um = u.copy(), u.copy()
            up[j] += h
            um[j] -= h
            fp = _nr_pile_force(pg, up, False)[0].copy()
            b1 = tuple(bool(b) for b in (np.abs(pg['_s'][0]) > pg['cap'][0]))
            fm = _nr_pile_force(pg, um, False)[0].copy()
            b2 = tuple(bool(b) for b in (np.abs(pg['_s'][0]) > pg['cap'][0]))
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
            f"branch is affine in the element displacement, so the tangent is exact "
            f"and this can only be an error in the element.")

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
          "power-curve, softening-bar and pile models are refused by "
          "name, a reinforced bisection lands on its measured strength with every "
          "bar force inside the capacity its embedment develops, a cohesionless "
          "base does not manufacture a failure at a strength the driver's "
          "own converged answer above it proves standing, a Rankine tensile cap "
          "returns to both yield surfaces at once and is reduced with the trial "
          "strength while an inert cap changes nothing at all, the K0 initial "
          "stress is an exact equilibrium on level ground down to its "
          "out-of-plane component and reproduces a locked vendor benchmark on "
          "both drivers, the matric-suction apparent cohesion is the same field "
          "on both drivers and is bounded by its cap and reduced by the trial "
          "strength, and the "
          "environment override cannot swap the driver in silence. The monotonic "
          "ramp reaches the same limit along one warm-started history, reports it "
          "on the bisection's midpoint convention, and never solves past it.")


if __name__ == '__main__':
    main()
