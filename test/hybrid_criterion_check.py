"""Unit checks for the SSRM's opt-in HYBRID failure criterion.

The hybrid criterion keeps non-convergence as the trigger for "failed" but requires
DISPLACEMENT EVIDENCE before a non-converged trial is counted as a failed slope:
max|u| beyond the trial's own elastic scale AND still growing. A trial frozen at
elastic scale is STABLE_STUCK — the solver failed, not the slope — and the bisection
treats it as standing. Anything in between keeps the legacy verdict, flagged
AMBIGUOUS. See xslope.fem.classify_nonconvergence and docs/fem/overview.md.

What this file locks:

  1. CLASSIFIER TABLE. The calibrated verdict for each corner of the (u_ratio,
     growth) plane, including the deliberate in-between band and the degenerate
     inputs (no elastic scale, too little history). Both signals are required for
     FAILED and both must be absent for STABLE_STUCK.

  2. BISECTION WIRING. Driving the real bisection driver (_ssrm_displacement_limit)
     with a stubbed solve_fem: a STABLE_STUCK trial moves the bracket UP under
     'hybrid' and DOWN under the default, an AMBIGUOUS trial moves it down under
     BOTH (legacy verdict stands), and the per-trial verdict table comes back on
     both settings.

  3. DEFAULT INVARIANCE. On a real solve, `stable` equals `converged` under the
     default criterion — the metadata is recorded but changes nothing — and the
     displacement history is actually populated (a converged trial reports
     verdict='CONVERGED').

Run directly:  PYTHONPATH=. python3 test/hybrid_criterion_check.py
Exits non-zero on any failure.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import xslope.fem as fem                                    # noqa: E402
from xslope.fem import classify_nonconvergence              # noqa: E402

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


# ===================== 1. classifier table =====================

def check_classifier():
    print("\n1. classify_nonconvergence — calibrated verdict table")

    # Measured stuck signature (RS2-62c, 2026-07-26): 1.0-1.1x elastic, frozen
    # across a 10k -> 80k budget increase. Both signals absent -> STABLE_STUCK.
    for ratio in (1.00, 1.05, 1.10):
        v, ur, g = classify_nonconvergence([ratio] * 40, 1.0, 'iteration_cap')
        check(f"frozen at {ratio:.2f}x elastic -> STABLE_STUCK",
              v == 'STABLE_STUCK' and abs(ur - ratio) < 1e-12 and abs(g) < 1e-12)

    # Measured failing signature: 4-21x elastic and still climbing.
    v, ur, g = classify_nonconvergence([1.0 + 0.5 * i for i in range(40)], 1.0,
                                       'iteration_cap')
    check("runaway (20x elastic, growing) -> FAILED",
          v == 'FAILED' and ur > 4.0 and g > 0.0, f"u_ratio={ur:.2f} growth={g:.2f}")

    # Griffiths & Lane Example 1's bisection-relevant band: ~1.5-3x and growing.
    v, ur, g = classify_nonconvergence([1.9 + 0.005 * i for i in range(40)], 1.0,
                                       'iteration_cap')
    check("G&L upturn band (2x elastic, growing) -> FAILED",
          v == 'FAILED', f"verdict={v} u_ratio={ur:.2f} growth={g:.3f}")

    # In-between: ONE signal without the other. Legacy verdict stands, flagged.
    v, _, _ = classify_nonconvergence([6.0] * 40, 1.0, 'iteration_cap')
    check("large but frozen -> AMBIGUOUS (legacy verdict stands)", v == 'AMBIGUOUS')
    v, ur, g = classify_nonconvergence([1.0 + 0.005 * i for i in range(40)], 1.0,
                                       'iteration_cap')
    check("elastic-scale but growing -> AMBIGUOUS", v == 'AMBIGUOUS',
          f"u_ratio={ur:.2f} growth={g:.3f}")
    # The undecided band between the two thresholds, held frozen: still AMBIGUOUS
    # because u_ratio clears the STUCK ceiling.
    v, _, _ = classify_nonconvergence([1.35] * 40, 1.0, 'iteration_cap')
    check("in the 1.25-1.5 band, frozen -> AMBIGUOUS", v == 'AMBIGUOUS')

    # The displacement cap is corroboration in its own right.
    v, _, _ = classify_nonconvergence([1.02] * 40, 1.0, 'disp_limit')
    check("displacement cap tripped -> FAILED regardless of history",
          v == 'FAILED')

    # Degenerate inputs must never claim STABLE_STUCK.
    v, ur, g = classify_nonconvergence([], 0.0, 'iteration_cap')
    check("no elastic scale -> AMBIGUOUS, no numbers",
          v == 'AMBIGUOUS' and ur is None and g is None)
    v, _, g = classify_nonconvergence([1.0, 1.0, 1.0], 1.0, 'iteration_cap')
    check("too little history -> AMBIGUOUS, no growth reported",
          v == 'AMBIGUOUS' and g is None)

    # A residual plateau is NOT itself a verdict: on a frozen state the classifier
    # still reads STABLE_STUCK. ('no_progress' is the legacy exit_reason string,
    # still accepted; solve_fem no longer ends a solve on a plateau.)
    v, _, _ = classify_nonconvergence([1.04] * 40, 1.0, 'no_progress')
    check("a plateau on a frozen state -> STABLE_STUCK", v == 'STABLE_STUCK')

    # Thresholds are where the calibration says they are.
    check("calibration constants unchanged",
          (fem._HYBRID_U_STUCK_MAX, fem._HYBRID_U_FAIL_MIN,
           fem._HYBRID_GROWTH_MIN, fem._HYBRID_WINDOW_FRAC)
          == (1.25, 1.5, 0.02, 0.25))


# ===================== 2. bisection wiring (stubbed solve_fem) =====================

def _stub_solution(F, kind):
    """A minimal solve_fem-shaped result. `kind` selects the verdict."""
    if kind == 'converged':
        return {"converged": True, "stable": True, "verdict": "CONVERGED",
                "u_ratio": None, "u_growth": None, "exit_reason": "converged",
                "iterations": 100}
    if kind == 'stuck':
        return {"converged": False, "stable": True, "verdict": "STABLE_STUCK",
                "u_ratio": 1.05, "u_growth": 0.0, "exit_reason": "no_progress",
                "iterations": 1600}
    if kind == 'ambiguous':
        return {"converged": False, "stable": False, "verdict": "AMBIGUOUS",
                "u_ratio": 6.0, "u_growth": 0.0, "exit_reason": "no_progress",
                "iterations": 1600}
    if kind == 'diverging':
        # Closed early by the early-failure rule: a failure verdict already made,
        # on evidence the trial produced before its budget ran out.
        return {"converged": False, "stable": False, "verdict": "FAILED",
                "u_ratio": 9.4, "u_growth": 1.1, "exit_reason": "diverging",
                "iterations": 330, "diverging_iteration": 330}
    if kind == 'inconclusive':
        # Reached the iteration ceiling with the residual still coming down:
        # neither settled nor failed.
        return {"converged": False, "stable": False, "verdict": "AMBIGUOUS",
                "u_ratio": 1.4, "u_growth": 0.0, "exit_reason": "inconclusive",
                "iterations": 50000}
    return {"converged": False, "stable": False, "verdict": "FAILED",
            "u_ratio": 12.0, "u_growth": 3.0, "exit_reason": "iteration_cap",
            "iterations": 1600}


class _stub_solve_fem:
    """Replace fem.solve_fem with a table lookup on F: `kind_of(F) -> kind`.

    _ssrm_displacement_limit calls solve_fem by bare module-global name, so the
    real bisection driver runs against the stub — the wiring under test is the
    driver's, not a re-implementation of it."""

    def __init__(self, kind_of):
        self._kind_of = kind_of
        self._orig = fem.solve_fem
        self.calls = []

    def __enter__(self):
        def _fake(fem_data, F=1.0, **kw):
            kind = self._kind_of(F)
            self.calls.append((F, kind, kw.get('failure_criterion')))
            sol = _stub_solution(F, kind)
            # solve_fem always reports `stable`; under the default criterion it
            # equals `converged`. Mirror that contract exactly.
            if kw.get('failure_criterion') != 'hybrid':
                sol = dict(sol, stable=sol["converged"])
            return sol
        fem.solve_fem = _fake
        return self

    def __exit__(self, *exc):
        fem.solve_fem = self._orig
        return False


def _bisect(kind_of, hybrid):
    with _stub_solve_fem(kind_of) as stub:
        res = fem._ssrm_displacement_limit(
            {}, F_min=1.0, F_max=2.0, tolerance=0.05, max_disp_factor=None,
            debug_level=0, hybrid=hybrid)
    return res, stub


def check_wiring():
    print("\n2. bisection wiring — the verdict a trial gets is the verdict the "
          "bracket uses")

    # A model that CONVERGES below 1.4, is STUCK from 1.4 to 1.8, and genuinely
    # FAILS above 1.8. The legacy criterion reads the stuck band as failure and
    # reports ~1.4; the hybrid walks through it and reports ~1.8.
    def kinds(F):
        if F < 1.4:
            return 'converged'
        if F < 1.8:
            return 'stuck'
        return 'failed'

    res_old, _ = _bisect(kinds, hybrid=False)
    res_new, _ = _bisect(kinds, hybrid=True)
    check("default criterion stops at the stuck band",
          abs(res_old['FS'] - 1.4) < 0.05, f"FS={res_old['FS']:.3f}")
    check("hybrid walks through the stuck band to the real boundary",
          abs(res_new['FS'] - 1.8) < 0.05, f"FS={res_new['FS']:.3f}")
    check("hybrid reports the method it used",
          'Hybrid' in res_new['method'] and 'Hybrid' not in res_old['method'])

    # AMBIGUOUS must NOT be rescued: the legacy verdict stands under both.
    def kinds_amb(F):
        return 'converged' if F < 1.4 else 'ambiguous'

    res_old, _ = _bisect(kinds_amb, hybrid=False)
    res_new, _ = _bisect(kinds_amb, hybrid=True)
    check("AMBIGUOUS keeps the legacy verdict on both criteria",
          abs(res_old['FS'] - res_new['FS']) < 1e-12,
          f"{res_old['FS']:.4f} vs {res_new['FS']:.4f}")

    # Per-trial metadata is recorded on BOTH settings (so an A/B costs no extra
    # solves), and every trial carries its verdict.
    for label, res in (('default', res_old), ('hybrid', res_new)):
        trials = res.get('trials')
        check(f"{label} returns a per-trial verdict table",
              bool(trials) and all(
                  {'F', 'role', 'converged', 'stable', 'verdict', 'u_ratio',
                   'growth', 'exit_reason', 'iterations'} <= set(t) for t in trials),
              f"{len(trials or [])} trials")

    # An INCONCLUSIVE trial — the iteration ceiling reached with the residual still
    # falling — is not a failure. The bisection carries on below it and reports the
    # final bracket's MIDPOINT, exactly as on any other run; what the trial changes is
    # what the bracket's upper edge means, which the result states in words.
    def kinds_inc(F):
        return 'converged' if F < 1.4 else ('inconclusive' if F < 1.8 else 'failed')

    res_inc, _ = _bisect(kinds_inc, hybrid=False)
    _lo, _hi = res_inc['final_interval']
    check("an inconclusive trial is still reported as the bracket midpoint",
          abs(res_inc['FS'] - 0.5 * (_lo + _hi)) < 1e-12 and _hi >= 1.4,
          f"FS={res_inc['FS']:.3f} interval={res_inc['final_interval']}")
    check("the run reports the inconclusive trial in words",
          bool(res_inc.get('note')) and 'inconclusive' in res_inc['note'].lower()
          and len(res_inc.get('inconclusive') or []) >= 1,
          f"note={res_inc.get('note')}")
    check("the factor of safety comes from a trial that reached equilibrium",
          res_inc['last_solution']['converged'] is True)

    # The criterion is threaded down to every trial, not just to the driver.
    _, stub = _bisect(kinds, hybrid=True)
    check("every trial is solved with failure_criterion='hybrid'",
          bool(stub.calls) and all(c[2] == 'hybrid' for c in stub.calls))
    _, stub = _bisect(kinds, hybrid=False)
    check("...and with 'non_convergence' on the default",
          bool(stub.calls) and all(c[2] == 'non_convergence' for c in stub.calls))

    # A trial closed early by the early-failure rule is a FAILURE on both criteria:
    # the rule fires only where the evidence is unambiguous, so there is nothing for
    # the hybrid to rescue.
    def kinds_div(F):
        return 'converged' if F < 1.4 else 'diverging'

    res_old, _ = _bisect(kinds_div, hybrid=False)
    res_new, _ = _bisect(kinds_div, hybrid=True)
    check("a 'diverging' trial moves the bracket DOWN on both criteria",
          abs(res_old['FS'] - 1.4) < 0.05 and abs(res_new['FS'] - 1.4) < 0.05,
          f"{res_old['FS']:.3f} vs {res_new['FS']:.3f}")


# ===================== 3. default invariance on a real solve =====================

def _real_model():
    """The Griffiths & Lane Example 1 slope on a coarse mesh, or None if absent."""
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    xlsx = os.path.join(here, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')
    if not os.path.exists(xlsx):
        return None
    slope_data = load_slope_data(xlsx)
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=8.0, element_type='tri6')
    return build_fem_data(slope_data, mesh)


def check_real_solve():
    print("\n3. real solve — metadata is recorded, default behavior is unchanged")

    from xslope.fem import solve_fem

    fem_data = _real_model()
    if fem_data is None:
        check("griffiths1 input present", False)
        return

    # F = 1.0 on this slope (FS ~ 1.4) settles quickly: a converged trial.
    sol = solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=4000,
                    max_disp_factor=None)
    check("converged trial reports verdict='CONVERGED'",
          sol['verdict'] == 'CONVERGED', f"converged={sol['converged']}")
    check("converged trial: stable == converged", sol['stable'] == sol['converged'])
    check("exit_reason recorded", sol['exit_reason'] == 'converged')
    check("elastic displacement scale is positive",
          sol['u_elastic_scale'] > 0.0, f"{sol['u_elastic_scale']:.4g}")

    # F well past failure: a non-converged trial with real displacement evidence.
    # `stable` must still equal `converged` because the criterion is the default.
    sol = solve_fem(fem_data, F=2.5, debug_level=0, max_iterations=1200,
                    max_disp_factor=None, early_exit=False)
    check("failing trial does not converge", sol['converged'] is False)
    check("failing trial: default criterion leaves stable == converged",
          sol['stable'] == sol['converged'],
          f"verdict={sol['verdict']} u_ratio={sol['u_ratio']}")
    check("failing trial carries displacement evidence",
          sol['u_ratio'] is not None and sol['u_ratio'] > 1.0,
          f"u_ratio={sol['u_ratio']:.2f}")
    check("failing trial is classified FAILED", sol['verdict'] == 'FAILED',
          f"verdict={sol['verdict']} growth={sol['u_growth']}")

    # Same solve under 'hybrid': a genuinely failing trial must NOT be rescued.
    sol_h = solve_fem(fem_data, F=2.5, debug_level=0, max_iterations=1200,
                      max_disp_factor=None, early_exit=False,
                      failure_criterion='hybrid')
    check("hybrid does not rescue a genuinely failing trial",
          sol_h['stable'] is False and sol_h['verdict'] == 'FAILED')

    # The no-progress watch cannot hand the classifier a truncated history,
    # because it no longer stops anything: a plateau in the out-of-balance residual
    # is recorded and the solve runs on. A trial ends on convergence, on its
    # iteration cap, or on the displacement cap, and on nothing else.
    check("the no-progress watch is reported",
          'plateau_iteration' in sol_h and 'early_exit_suppressed' in sol_h)
    # early_failure off: this trial IS running away, and the early-failure rule
    # (section 6) would close it long before the budget. What is under test here is
    # the plateau watch, which must not end a solve on its own.
    sol_d = solve_fem(fem_data, F=2.5, debug_level=0, max_iterations=2000,
                      max_disp_factor=None, early_exit=True, early_failure=False)
    check("a residual plateau never ends a solve",
          sol_d['exit_reason'] in ('converged', 'iteration_cap', 'disp_limit'),
          f"exit_reason={sol_d['exit_reason']}")
    check("a non-converged trial spends its whole budget",
          sol_d['converged'] is False and sol_d['iterations'] == 2000,
          f"iterations={sol_d['iterations']}")
    check("the early-failure fields are reported on every solve",
          'diverging_iteration' in sol_d and sol_d['diverging_iteration'] is None,
          f"diverging_iteration={sol_d['diverging_iteration']}")


def check_budget_extension():
    """The budget is extended while the residual is still falling, and the hard
    ceiling produces 'inconclusive' rather than a failure.

    The extension is the VISCOPLASTIC driver's contract: it marches a fixed
    number of steps, so running out of them is a statement about the budget and
    not about the slope. The Newton-Raphson driver reaches equilibrium by a
    different route and reports `budget_extensions` of zero on every trial, so
    the driver is named here rather than left to the default — on this model at
    F = 1.35 the default converges in 303 iterations and never extends, which
    tests nothing about the rule.
    """
    print("\n5. budget extension — reaching the budget is not a verdict either")
    from xslope.fem import solve_fem

    fem_data = _real_model()
    if fem_data is None:
        check("griffiths1 input present", False)
        return

    # Just below this slope's factor of safety (~1.37), where equilibrium takes
    # hundreds of iterations. A budget below what the trial needs: the residual
    # is still falling when it runs out, so the solve is extended rather than failed.
    sol = solve_fem(fem_data, F=1.35, debug_level=0, max_iterations=300,
                    max_iterations_ceiling=8000, max_disp_factor=None,
                    fem_solver='viscoplastic')
    check("a too-small budget is extended, not failed",
          sol['budget_extensions'] > 0 and sol['iterations'] > 300,
          f"extensions={sol['budget_extensions']} iterations={sol['iterations']}")
    check("and the extended trial reaches equilibrium", sol['converged'] is True,
          f"exit_reason={sol['exit_reason']}")

    # The same trial with the ceiling at the budget: nothing may be extended, and
    # stopping while still improving is INCONCLUSIVE, not failure.
    sol = solve_fem(fem_data, F=1.35, debug_level=0, max_iterations=300,
                    max_iterations_ceiling=300, max_disp_factor=None,
                    fem_solver='viscoplastic')
    check("the ceiling stops the extension", sol['budget_extensions'] == 0,
          f"iterations={sol['iterations']}")
    check("a solve still improving at the ceiling is inconclusive",
          sol['exit_reason'] == 'inconclusive' and sol['converged'] is False,
          f"exit_reason={sol['exit_reason']}")


def check_exit_suppression():
    """The suppression rule, on a synthetic solve whose history is controlled.

    RS2-62c at F = 0.800 is the measured case: the residual plateaus at iteration
    5,118 with max|u| at 1.03x elastic and no trailing growth — which the classifier
    alone reads as STABLE_STUCK — while the same trial run to 40,000 iterations
    reaches 1.72x and is growing. The rule under test is that a stuck-LOOKING state
    does not get its verdict taken on a short history."""
    print("\n4. a truncated history cannot produce STABLE_STUCK")

    # Truncated (what the exit would have seen) vs full budget (the truth).
    truncated = [1.03] * 40
    full = [1.03] * 40 + [1.03 + 0.7 * i / 100 for i in range(100)]
    v_trunc, _, _ = classify_nonconvergence(truncated, 1.0, 'no_progress')
    v_full, ur, g = classify_nonconvergence(full, 1.0, 'iteration_cap')
    check("the truncated history alone reads STABLE_STUCK",
          v_trunc == 'STABLE_STUCK')
    check("the full-budget history reads FAILED", v_full == 'FAILED',
          f"u_ratio={ur:.2f} growth={g:.2f}")
    check("so the two disagree — which is why a plateau may not stop the solve",
          v_trunc != v_full)


# ===================== 6. the early-failure rule =====================

def _series(txt):
    """A fixture series, one float per comma-separated token."""
    return [float(t) for t in txt.replace('\n', '').split(',')]


def _fires_at(disp, oob, **kw):
    """Replay a recorded trial exactly as solve_fem feeds it to the rule: prefix by
    prefix, one sample every _HYBRID_SAMPLE_EVERY iterations. Returns
    (iteration, signal) at the first firing, or (None, None)."""
    for j in range(len(disp)):
        sig = fem._early_failure(disp[:j + 1], oob[:j + 1], 1.0, **kw)
        if sig is not None:
            return j * fem._HYBRID_SAMPLE_EVERY, sig
    return None, None


def check_early_failure():
    """The rule that closes a running-away trial before its budget runs out.

    Its whole safety argument is a threshold placement, so that is what is locked
    here — against the recorded sampled series of three real strength-reduction
    trials (see the fixtures at the bottom of this file). Two are gross failures,
    one on each of the rule's two tests; the third is the near-critical trial that
    the reinforced slope's factor of safety turns on, which grows past FIVE times
    its elastic displacement, sits with a flat residual for thousands of
    iterations, and THEN converges at iteration 20,085. The rule must catch the
    first two and must not touch the third — and the mutation checks below show
    that the third is a live constraint, not a fixture that passes whatever the
    thresholds are: moving either threshold a little way down fires on it, which
    would move the factor of safety.
    """
    print("\n6. early failure — the thresholds, on recorded trials")

    check("early-failure constants unchanged",
          (fem._EARLY_FAIL_WINDOW, fem._EARLY_FAIL_WARMUP,
           fem._EARLY_FAIL_GAIN, fem._EARLY_FAIL_U_MAX) == (2000, 500, 1.0, 8.0),
          f"{fem._EARLY_FAIL_WINDOW}/{fem._EARLY_FAIL_WARMUP}/"
          f"{fem._EARLY_FAIL_GAIN}/{fem._EARLY_FAIL_U_MAX}")

    # (f) — the absolute-runaway test, on a trial that reaches 78 elastic
    # displacements. Recorded firing: iteration 1,160 of the 12,000 it used to spend.
    d, o = _series(RUNAWAY_DISP), _series(RUNAWAY_OOB)
    it, sig = _fires_at(d, o)
    check("gross runaway is caught, on the recorded iteration",
          (it, sig) == (1160, 'runaway'), f"{it} ({sig})")

    # (b) — the stalled-residual test, on a trial whose residual goes flat while the
    # field keeps moving. Recorded firing: iteration 3,140 of 36,000 (two extensions).
    d, o = _series(STALLED_DISP), _series(STALLED_OOB)
    it, sig = _fires_at(d, o)
    check("a stalled residual with the field still moving is caught",
          (it, sig) == (3140, 'stalled_residual'), f"{it} ({sig})")

    # The near-critical CONVERGED trial: silent over its whole 20,085 iterations.
    d, o = _series(NEARCRIT_DISP), _series(NEARCRIT_OOB)
    it, sig = _fires_at(d, o)
    check("the near-critical trial that CONVERGES is never touched",
          it is None, f"fired at {it} ({sig})")
    check("...and it is genuinely near the thresholds",
          4.5 < max(d) < 8.0, f"reaches {max(d):.3f} elastic displacements")

    # Mutation: the margins are 1.59x on the level and 3.0x on the window gain, so a
    # threshold moved inside them fires on a trial that reaches equilibrium.
    it_u, _ = _fires_at(d, o, u_max=5.0)
    check("a level threshold lowered to 5x elastic WOULD fire on it",
          it_u is not None, f"fired at {it_u}")
    it_g, _ = _fires_at(d, o, gain=0.3)
    check("a window-gain floor lowered to 0.3 WOULD fire on it",
          it_g is not None, f"fired at {it_g}")

    # Degenerate inputs: no elastic yardstick, no history.
    check("no elastic scale -> no verdict",
          fem._early_failure([9.0] * 40, [1.0] * 40, 0.0) is None)
    check("too little history -> no verdict",
          fem._early_failure([9.0, 9.5], [1.0, 1.0], 1.0) is None)


def check_early_failure_wiring():
    """The rule inside a real solve: the exit it produces, and the opt-out."""
    print("\n7. early failure — the exit it produces on a real solve")

    from xslope.fem import solve_fem

    fem_data = _real_model()
    if fem_data is None:
        check("griffiths1 input present", False)
        return

    # F = 2.5 on a slope whose FS is ~1.4: a gross runaway, the case the rule is for.
    sol = solve_fem(fem_data, F=2.5, debug_level=0, max_iterations=4000,
                    max_disp_factor=None, early_exit=False)
    check("a running-away trial exits 'diverging'",
          sol['exit_reason'] == 'diverging', f"exit_reason={sol['exit_reason']}")
    check("...well before its budget",
          sol['converged'] is False and sol['iterations'] < 4000,
          f"iterations={sol['iterations']} of 4000")
    check("...reporting the iteration and the test that fired",
          sol['diverging_iteration'] == sol['iterations'] - 1
          and sol['diverging_signal'] in ('runaway', 'stalled_residual'),
          f"{sol['diverging_signal']} at {sol['diverging_iteration']}")
    check("...and the bisection reads it as FAILED",
          sol['verdict'] == 'FAILED' and sol['stable'] is False,
          f"verdict={sol['verdict']}")

    # The opt-out. The at-failure capture solve uses it: the capture exists to let
    # the mechanism develop, which is exactly the runaway the rule stops.
    sol_off = solve_fem(fem_data, F=2.5, debug_level=0, max_iterations=1000,
                        max_disp_factor=None, early_exit=False,
                        early_failure=False)
    check("early_failure=False spends the whole budget",
          sol_off['iterations'] >= 1000 and sol_off['exit_reason'] != 'diverging',
          f"iterations={sol_off['iterations']} exit={sol_off['exit_reason']}")
    check("...and the field it leaves is the developed one",
          sol_off['max_displacement'] > sol['max_displacement'],
          f"{sol_off['max_displacement']:.3g} vs {sol['max_displacement']:.3g}")


def main():
    print("=" * 72)
    print("Hybrid failure-criterion checks")
    print("=" * 72)
    check_classifier()
    check_wiring()
    check_real_solve()
    check_exit_suppression()
    check_budget_extension()
    check_early_failure()
    check_early_failure_wiring()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All hybrid failure-criterion checks passed.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


# ===================== recorded fixtures =====================
#
# The sampled series (max|u| and the maximum nodal out-of-balance ratio, one sample
# every _HYBRID_SAMPLE_EVERY = 10 iterations) of three real strength-reduction
# trials, taken from two complete bisection walks on the tutorial models:
#
#   RUNAWAY   the FEM-1 embankment (3.5 ft tri6) at F = 1.500000 — FAILED, and
#             gross: 78 elastic displacements. Truncated at the recorded firing.
#   STALLED   the FEM-2 reinforced slope (2.0 ft tri6, elastic-perfectly-plastic
#             bars) at F = 1.625000 — FAILED after two budget extensions, on the
#             residual test. Truncated at the recorded firing.
#   NEARCRIT  the same model at F = 1.546875 — CONVERGED at iteration 20,085, and
#             it is the trial the reinforced slope's factor of safety turns on. Its
#             whole life is kept: the rule must stay silent across all of it.
#
# max|u| is written in units of that trial's own elastic displacement, which is
# exact — every comparison the rule makes is homogeneous of degree one in
# (max|u|, u_elastic) — so the fixtures replay with u_elastic = 1.0.

#RUNAWAY
RUNAWAY_DISP = """\
1.0664,1.2738,1.3507,1.4117,1.4677,1.5234,1.5778,1.6318,1.6857,1.7394,1.7933,1.847,1.9004,
1.9535,2.0063,2.0589,2.1114,2.1638,2.216,2.2682,2.3203,2.3722,2.4241,2.476,2.5277,2.5794,
2.631,2.6826,2.7341,2.7855,2.837,2.8884,2.9398,2.9912,3.0425,3.0937,3.145,3.1962,3.2473,
3.2984,3.3495,3.4005,3.4516,3.5025,3.5535,3.6044,3.6552,3.7061,3.7569,3.8077,3.8584,3.9092,
3.9599,4.0105,4.0612,4.1118,4.1624,4.213,4.2636,4.3141,4.3684,4.4341,4.4998,4.5657,4.6317,
4.6977,4.7636,4.8295,4.8955,4.9614,5.0273,5.0931,5.159,5.2249,5.2907,5.3565,5.4223,5.4881,
5.5539,5.6197,5.6855,5.7512,5.817,5.8827,5.9484,6.0141,6.0798,6.1455,6.2112,6.2769,6.3425,
6.4082,6.4739,6.5397,6.6056,6.6714,6.7372,6.803,6.8688,6.9346,7.0004,7.0662,7.132,7.1977,
7.2635,7.3292,7.3949,7.4607,7.5264,7.5921,7.6578,7.7235,7.7892,7.8549,7.9205,7.9862,8.0519"""

RUNAWAY_OOB = """\
1,1.25,1.92,1.96,1.88,1.83,1.79,1.76,1.74,1.72,1.71,1.7,1.7,1.69,1.68,1.68,1.68,1.67,1.67,
1.67,1.67,1.67,1.66,1.66,1.66,1.66,1.66,1.66,1.66,1.66,1.66,1.66,1.65,1.65,1.65,1.65,1.65,
1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,
1.65,1.65,1.65,1.65,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,
1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,
1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.64,1.63,1.64,1.64,1.64,1.64,1.63,
1.63,1.63,1.63,1.63,1.63,1.63,1.63,1.63"""

#STALLED
STALLED_DISP = """\
1.0209,1.1803,1.3117,1.4392,1.5629,1.6787,1.7876,1.8903,1.9873,2.079,2.1658,2.2479,2.3255,
2.3991,2.4689,2.535,2.5978,2.6574,2.7145,2.7692,2.8213,2.871,2.9185,2.9639,3.0074,3.0491,
3.0891,3.1274,3.1642,3.1997,3.2339,3.2666,3.2983,3.3289,3.3585,3.3871,3.4149,3.4418,3.4679,
3.4933,3.518,3.5422,3.5656,3.5885,3.6178,3.6481,3.6778,3.7069,3.7352,3.763,3.7902,3.8168,
3.8429,3.8685,3.8936,3.9182,3.9424,3.9662,3.9896,4.0127,4.0354,4.0578,4.0799,4.1017,4.1232,
4.1445,4.1655,4.1867,4.208,4.229,4.2499,4.2705,4.2909,4.3111,4.3311,4.3509,4.3705,4.3899,
4.4091,4.4281,4.4469,4.4656,4.484,4.5023,4.5204,4.5384,4.5562,4.5738,4.5912,4.6085,4.6257,
4.6426,4.6594,4.6761,4.6926,4.709,4.7252,4.7413,4.7573,4.7731,4.7887,4.8043,4.8197,4.835,
4.8501,4.8652,4.8801,4.8949,4.9096,4.9242,4.9388,4.9532,4.9676,4.9819,4.9961,5.0103,5.0244,
5.0384,5.0524,5.0664,5.0803,5.0942,5.1081,5.122,5.1358,5.1496,5.1634,5.1772,5.1909,5.2047,
5.2184,5.2322,5.2459,5.2596,5.2733,5.287,5.3006,5.3143,5.3279,5.3415,5.3551,5.3686,5.3821,
5.3956,5.4091,5.4225,5.4359,5.4492,5.4625,5.4758,5.489,5.5022,5.5153,5.5283,5.5413,5.5542,
5.5671,5.5799,5.5927,5.6054,5.618,5.6305,5.643,5.6554,5.6677,5.68,5.6922,5.7043,5.7163,
5.7283,5.7402,5.752,5.7638,5.7755,5.7871,5.7986,5.8101,5.8215,5.8328,5.8441,5.8553,5.8664,
5.8775,5.8885,5.8994,5.9103,5.9211,5.9318,5.9425,5.9531,5.9637,5.9742,5.9847,5.9951,6.0055,
6.0158,6.0261,6.0363,6.0466,6.0568,6.0669,6.0771,6.0872,6.0973,6.1074,6.1174,6.1275,6.1375,
6.1476,6.1577,6.1677,6.1778,6.1879,6.198,6.2081,6.2182,6.2284,6.2386,6.2488,6.259,6.2693,
6.2796,6.29,6.3003,6.3108,6.3213,6.3318,6.3424,6.353,6.3636,6.3743,6.3851,6.3958,6.4066,
6.4175,6.4284,6.4393,6.4503,6.4613,6.4723,6.4834,6.4945,6.5056,6.5168,6.528,6.5393,6.5506,
6.5619,6.5734,6.5848,6.5963,6.6078,6.6194,6.6309,6.6426,6.6542,6.6659,6.6777,6.6894,6.7012,
6.713,6.7249,6.7368,6.7487,6.7606,6.7725,6.7845,6.7965,6.8084,6.8204,6.8324,6.8444,6.8564,
6.8684,6.8804,6.8924,6.9044,6.9164,6.9284,6.9403,6.9523,6.9642,6.976,6.9879,6.9997,7.0115,
7.0233,7.035,7.0466,7.0583,7.0699,7.0814,7.0929,7.1044,7.1158,7.1271,7.1384,7.1497,7.1609,
7.1721,7.1832,7.1943,7.2053,7.2164,7.2274,7.2383,7.2493,7.2602,7.2711,7.282,7.2929,7.3038,
7.3146,7.3255,7.3363"""

STALLED_OOB = """\
5.39,6.55,5.25,4.65,4.38,4.18,4,3.82,3.66,3.5,3.34,3.19,3.05,2.91,2.77,2.64,2.52,2.4,2.28,
2.17,2.09,2.04,1.99,1.94,1.9,1.86,1.82,1.78,1.74,1.71,1.67,1.64,1.61,1.58,1.55,1.53,1.5,1.48,
1.46,1.43,1.41,1.4,1.38,1.36,1.35,1.32,1.29,1.26,1.24,1.22,1.19,1.18,1.16,1.14,1.13,1.11,1.1,
1.08,1.07,1.06,1.06,1.05,1.05,1.04,1.04,1.04,1.03,1.03,1.03,1.02,1.02,1.02,1.02,1.01,1.01,
1.01,1.01,1.01,1,1,0.999,0.996,0.994,0.992,0.989,0.986,0.984,0.981,0.978,0.975,0.972,0.969,
0.966,0.963,0.959,0.956,0.953,0.95,0.946,0.943,0.94,0.936,0.933,0.929,0.926,0.923,0.92,0.917,
0.914,0.911,0.908,0.905,0.902,0.9,0.898,0.896,0.894,0.892,0.891,0.89,0.888,0.887,0.887,0.886,
0.885,0.885,0.885,0.884,0.884,0.884,0.884,0.884,0.885,0.885,0.885,0.885,0.885,0.885,0.885,
0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.886,0.885,0.885,0.885,
0.884,0.883,0.883,0.882,0.881,0.88,0.879,0.878,0.877,0.876,0.874,0.873,0.872,0.87,0.869,
0.867,0.865,0.864,0.862,0.86,0.858,0.856,0.854,0.852,0.85,0.848,0.846,0.843,0.841,0.839,
0.837,0.835,0.833,0.831,0.829,0.827,0.825,0.823,0.821,0.819,0.818,0.816,0.814,0.813,0.811,
0.81,0.808,0.807,0.806,0.804,0.803,0.802,0.801,0.8,0.8,0.799,0.798,0.798,0.798,0.797,0.797,
0.797,0.797,0.797,0.797,0.798,0.798,0.799,0.8,0.801,0.802,0.803,0.805,0.806,0.808,0.809,
0.811,0.812,0.814,0.816,0.817,0.819,0.821,0.823,0.824,0.826,0.828,0.83,0.832,0.834,0.836,
0.838,0.84,0.843,0.845,0.847,0.849,0.852,0.854,0.856,0.858,0.861,0.863,0.865,0.867,0.869,
0.871,0.873,0.874,0.876,0.878,0.879,0.881,0.882,0.883,0.884,0.885,0.886,0.887,0.887,0.888,
0.888,0.889,0.889,0.889,0.889,0.889,0.889,0.889,0.889,0.888,0.887,0.886,0.885,0.883,0.881,
0.88,0.878,0.876,0.874,0.871,0.869,0.866,0.863,0.861,0.858,0.855,0.852,0.849,0.846,0.843,
0.84,0.837,0.835,0.832,0.829,0.827,0.824,0.822,0.82,0.818"""

#NEARCRIT
NEARCRIT_DISP = """\
1.013,1.1357,1.2387,1.3393,1.4394,1.5333,1.6215,1.7046,1.7831,1.8572,1.9274,1.9938,2.0567,
2.1163,2.1727,2.2269,2.2784,2.3273,2.3737,2.4177,2.4596,2.4993,2.537,2.573,2.6073,2.6399,
2.6711,2.7009,2.7295,2.7566,2.7827,2.8078,2.8317,2.8546,2.8767,2.8977,2.918,2.9373,2.956,
2.9737,2.991,3.0073,3.0232,3.0383,3.053,3.0669,3.0805,3.0934,3.106,3.118,3.1297,3.1409,
3.1518,3.1621,3.1723,3.1819,3.1914,3.2004,3.2093,3.2177,3.226,3.2339,3.2417,3.2491,3.2564,
3.2633,3.2702,3.2767,3.2832,3.2894,3.2955,3.3014,3.3072,3.3127,3.3182,3.3235,3.3288,3.3338,
3.3388,3.3436,3.3484,3.353,3.3576,3.362,3.3664,3.3706,3.3749,3.3789,3.383,3.3869,3.3909,
3.3947,3.3985,3.4021,3.4058,3.4093,3.4129,3.4163,3.4197,3.423,3.4263,3.4295,3.4328,3.4359,
3.439,3.4421,3.4451,3.4481,3.451,3.4539,3.4568,3.4596,3.4624,3.4651,3.4679,3.4706,3.4732,
3.4759,3.4785,3.4811,3.4836,3.4862,3.4887,3.4913,3.4937,3.4962,3.4987,3.5011,3.5035,3.5059,
3.5083,3.5107,3.513,3.5153,3.5176,3.5198,3.5221,3.5243,3.5265,3.5286,3.5308,3.5329,3.535,
3.5371,3.5391,3.5412,3.5432,3.5452,3.5471,3.549,3.5509,3.5528,3.5546,3.5564,3.5582,3.56,
3.5617,3.5634,3.565,3.5666,3.5682,3.5698,3.5713,3.5728,3.5743,3.5757,3.5771,3.5791,3.5811,
3.5831,3.585,3.5869,3.5887,3.5905,3.5922,3.5939,3.5955,3.5971,3.5987,3.6001,3.6016,3.6029,
3.6043,3.6056,3.6069,3.6081,3.6093,3.6104,3.6116,3.6126,3.6137,3.6147,3.6157,3.6167,3.6176,
3.6186,3.6194,3.6203,3.6211,3.6219,3.6227,3.6235,3.6242,3.625,3.6257,3.6264,3.6271,3.6277,
3.6284,3.629,3.6296,3.6302,3.6308,3.6314,3.632,3.6326,3.6331,3.6337,3.6342,3.6347,3.6352,
3.6357,3.6362,3.6367,3.6372,3.6376,3.6381,3.6386,3.639,3.6394,3.6399,3.6403,3.6407,3.6411,
3.6415,3.6418,3.6422,3.6426,3.6429,3.6433,3.6436,3.6439,3.6442,3.6445,3.6448,3.6451,3.6454,
3.6457,3.646,3.6462,3.6465,3.6467,3.6469,3.6472,3.6474,3.6476,3.6478,3.648,3.6482,3.6484,
3.6486,3.6487,3.6489,3.6491,3.6492,3.6494,3.6495,3.6497,3.6498,3.65,3.6501,3.6503,3.6504,
3.6505,3.6507,3.6508,3.651,3.6511,3.6513,3.6514,3.6516,3.6517,3.6519,3.6521,3.6522,3.6524,
3.6526,3.6528,3.653,3.6532,3.6534,3.6536,3.6538,3.654,3.6542,3.6545,3.6547,3.655,3.6552,
3.6555,3.6557,3.656,3.6563,3.6566,3.6569,3.6572,3.6575,3.6578,3.6582,3.6585,3.6589,3.6592,
3.6596,3.66,3.6604,3.6609,3.6613,3.6617,3.6622,3.6627,3.6632,3.6637,3.6642,3.6647,3.6653,
3.6658,3.6664,3.667,3.6676,3.6682,3.6689,3.6695,3.6702,3.6709,3.6716,3.6724,3.6731,3.6739,
3.6747,3.6755,3.6764,3.6772,3.6781,3.679,3.6799,3.6808,3.6818,3.6828,3.6838,3.6848,3.6858,
3.6869,3.6879,3.689,3.6901,3.6913,3.6924,3.6936,3.6948,3.696,3.6973,3.6985,3.6998,3.7011,
3.7024,3.7038,3.7051,3.7065,3.7079,3.7093,3.7108,3.7123,3.7138,3.7153,3.7168,3.7184,3.72,
3.7216,3.7233,3.7249,3.7266,3.7283,3.73,3.7318,3.7335,3.7353,3.7371,3.7389,3.7407,3.7426,
3.7445,3.7463,3.7482,3.7502,3.7521,3.754,3.756,3.7579,3.7599,3.7619,3.7639,3.7659,3.7679,
3.7699,3.7719,3.7739,3.7759,3.778,3.78,3.782,3.784,3.786,3.788,3.7899,3.7919,3.7939,3.7959,
3.7979,3.7999,3.8019,3.8038,3.8058,3.8078,3.8098,3.8117,3.8137,3.8157,3.8176,3.8196,3.8216,
3.8235,3.8255,3.8275,3.8295,3.8315,3.8334,3.8354,3.8374,3.8394,3.8414,3.8434,3.8453,3.8473,
3.8493,3.8512,3.8532,3.8552,3.8571,3.8591,3.861,3.8629,3.8649,3.8668,3.8687,3.8706,3.8725,
3.8744,3.8763,3.8782,3.88,3.8819,3.8837,3.8856,3.8874,3.8892,3.891,3.8928,3.8946,3.8964,
3.8981,3.8999,3.9016,3.9033,3.905,3.9067,3.9084,3.9101,3.9117,3.9133,3.9149,3.9165,3.9181,
3.9197,3.9212,3.9227,3.9242,3.9256,3.9271,3.9285,3.9299,3.9312,3.9326,3.9339,3.9352,3.9364,
3.9377,3.9389,3.9401,3.9412,3.9424,3.9435,3.9446,3.9457,3.9468,3.9479,3.9489,3.9499,3.9509,
3.9519,3.9529,3.9539,3.9548,3.9558,3.9567,3.9576,3.9585,3.9594,3.9602,3.9611,3.9619,3.9628,
3.9636,3.9644,3.9653,3.9661,3.9669,3.9677,3.9685,3.9692,3.97,3.9708,3.9716,3.9723,3.9731,
3.9738,3.9746,3.9753,3.9761,3.9768,3.9776,3.9783,3.979,3.9798,3.9805,3.9813,3.982,3.9828,
3.9835,3.9843,3.985,3.9858,3.9866,3.9874,3.9881,3.989,3.9898,3.9906,3.9914,3.9923,3.9931,
3.994,3.9949,3.9958,3.9968,3.9977,3.9986,3.9996,4.0006,4.0016,4.0026,4.0037,4.0047,4.0058,
4.0069,4.008,4.0091,4.0103,4.0115,4.0127,4.0139,4.0152,4.0165,4.0179,4.0193,4.0207,4.0222,
4.0237,4.0252,4.0268,4.0285,4.0302,4.032,4.0338,4.0357,4.0377,4.0397,4.0418,4.0439,4.0461,
4.0484,4.0507,4.0531,4.0555,4.058,4.0605,4.0631,4.0657,4.0684,4.0711,4.0739,4.0767,4.0796,
4.0825,4.0854,4.0884,4.0915,4.0945,4.0976,4.1007,4.1038,4.1069,4.11,4.1132,4.1163,4.1195,
4.1226,4.1258,4.1289,4.132,4.1352,4.1383,4.1414,4.1444,4.1474,4.1504,4.1534,4.1563,4.1592,
4.1621,4.1649,4.1677,4.1704,4.1732,4.1758,4.1785,4.181,4.1836,4.1861,4.1885,4.1909,4.1933,
4.1956,4.1979,4.2001,4.2023,4.2044,4.2064,4.2085,4.2104,4.2123,4.2142,4.216,4.2178,4.2195,
4.2212,4.2228,4.2243,4.2258,4.2273,4.2287,4.23,4.2314,4.2326,4.2339,4.2351,4.2362,4.2374,
4.2385,4.2395,4.2406,4.2416,4.2426,4.2435,4.2445,4.2454,4.2463,4.2471,4.248,4.2488,4.2497,
4.2505,4.2512,4.252,4.2528,4.2535,4.2542,4.2549,4.2556,4.2563,4.2569,4.2576,4.2582,4.2588,
4.2594,4.26,4.2606,4.2611,4.2617,4.2622,4.2628,4.2633,4.2638,4.2643,4.2648,4.2653,4.2658,
4.2662,4.2667,4.2671,4.2676,4.268,4.2685,4.2689,4.2693,4.2697,4.2701,4.2705,4.2709,4.2713,
4.2716,4.272,4.2724,4.2727,4.2731,4.2734,4.2738,4.2741,4.2745,4.2748,4.2751,4.2754,4.2758,
4.2761,4.2764,4.2767,4.277,4.2773,4.2776,4.2779,4.2782,4.2785,4.2788,4.2791,4.2794,4.2797,
4.28,4.2803,4.2806,4.2809,4.2812,4.2815,4.2818,4.2821,4.2824,4.2827,4.283,4.2833,4.2836,
4.2839,4.2842,4.2845,4.2848,4.2851,4.2854,4.2857,4.286,4.2863,4.2867,4.287,4.2873,4.2876,
4.2879,4.2883,4.2886,4.2889,4.2893,4.2896,4.29,4.2903,4.2907,4.291,4.2914,4.2917,4.2921,
4.2924,4.2928,4.2932,4.2936,4.294,4.2944,4.2948,4.2952,4.2956,4.2961,4.2965,4.2969,4.2974,
4.2978,4.2983,4.2988,4.2993,4.2998,4.3003,4.3008,4.3013,4.3019,4.3024,4.303,4.3036,4.3041,
4.3047,4.3054,4.306,4.3066,4.3073,4.308,4.3086,4.3093,4.31,4.3108,4.3115,4.3123,4.313,4.3138,
4.3146,4.3154,4.3162,4.3171,4.318,4.3188,4.3197,4.3206,4.3215,4.3225,4.3234,4.3244,4.3254,
4.3263,4.3273,4.3283,4.3294,4.3304,4.3314,4.3325,4.3336,4.3346,4.3358,4.3369,4.338,4.3392,
4.3404,4.3416,4.3428,4.3441,4.3454,4.3467,4.348,4.3494,4.3507,4.3521,4.3536,4.355,4.3565,
4.3579,4.3594,4.3609,4.3625,4.364,4.3655,4.3671,4.3687,4.3703,4.3719,4.3735,4.3751,4.3767,
4.3784,4.38,4.3817,4.3834,4.3851,4.3868,4.3885,4.3902,4.3919,4.3936,4.3953,4.3971,4.3988,
4.4005,4.4022,4.404,4.4057,4.4074,4.4092,4.4109,4.4127,4.4144,4.4162,4.4179,4.4197,4.4214,
4.4232,4.4249,4.4267,4.4284,4.4301,4.4318,4.4335,4.4352,4.4369,4.4386,4.4402,4.4419,4.4435,
4.4451,4.4467,4.4483,4.4499,4.4515,4.453,4.4545,4.456,4.4575,4.459,4.4604,4.4619,4.4633,
4.4647,4.466,4.4674,4.4687,4.47,4.4712,4.4725,4.4737,4.475,4.4762,4.4773,4.4785,4.4797,
4.4808,4.4819,4.483,4.4841,4.4852,4.4862,4.4873,4.4883,4.4893,4.4903,4.4912,4.4922,4.4931,
4.494,4.4949,4.4957,4.4966,4.4974,4.4982,4.499,4.4997,4.5005,4.5012,4.5019,4.5026,4.5033,
4.504,4.5046,4.5052,4.5059,4.5065,4.507,4.5076,4.5082,4.5087,4.5092,4.5097,4.5102,4.5107,
4.5112,4.5117,4.5121,4.5126,4.513,4.5135,4.5139,4.5143,4.5147,4.5151,4.5154,4.5158,4.5162,
4.5165,4.5169,4.5172,4.5175,4.5178,4.5182,4.5185,4.5188,4.5191,4.5193,4.5196,4.5199,4.5202,
4.5204,4.5207,4.5209,4.5212,4.5214,4.5217,4.5219,4.5221,4.5224,4.5226,4.5228,4.523,4.5232,
4.5234,4.5236,4.5238,4.524,4.5242,4.5244,4.5246,4.5248,4.5249,4.5251,4.5253,4.5254,4.5256,
4.5258,4.5259,4.5261,4.5262,4.5264,4.5265,4.5267,4.5268,4.527,4.5271,4.5272,4.5274,4.5275,
4.5276,4.5278,4.5279,4.528,4.5282,4.5283,4.5284,4.5285,4.5286,4.5288,4.5289,4.529,4.5291,
4.5292,4.5293,4.5294,4.5295,4.5296,4.5297,4.5298,4.5299,4.53,4.5301,4.5302,4.5303,4.5304,
4.5305,4.5306,4.5307,4.5308,4.5309,4.531,4.531,4.5311,4.5312,4.5313,4.5314,4.5315,4.5315,
4.5316,4.5317,4.5318,4.5319,4.5319,4.532,4.5321,4.5322,4.5322,4.5323,4.5324,4.5324,4.5325,
4.5326,4.5327,4.5327,4.5328,4.5329,4.5329,4.533,4.5331,4.5331,4.5332,4.5333,4.5333,4.5334,
4.5335,4.5335,4.5336,4.5337,4.5337,4.5338,4.5339,4.5339,4.534,4.5341,4.5341,4.5342,4.5342,
4.5343,4.5344,4.5344,4.5345,4.5346,4.5346,4.5347,4.5347,4.5348,4.5349,4.5349,4.535,4.535,
4.5351,4.5352,4.5352,4.5353,4.5353,4.5354,4.5354,4.5355,4.5356,4.5356,4.5357,4.5357,4.5358,
4.5359,4.5359,4.536,4.536,4.5361,4.5361,4.5362,4.5363,4.5363,4.5364,4.5364,4.5365,4.5365,
4.5366,4.5366,4.5367,4.5368,4.5368,4.5369,4.5369,4.537,4.537,4.5371,4.5372,4.5372,4.5373,
4.5373,4.5374,4.5374,4.5375,4.5375,4.5376,4.5377,4.5377,4.5378,4.5378,4.5379,4.538,4.538,
4.5381,4.5381,4.5382,4.5382,4.5383,4.5384,4.5384,4.5385,4.5386,4.5386,4.5387,4.5387,4.5388,
4.5389,4.5389,4.539,4.5391,4.5391,4.5392,4.5393,4.5393,4.5394,4.5395,4.5395,4.5396,4.5397,
4.5397,4.5398,4.5399,4.54,4.54,4.5401,4.5402,4.5403,4.5403,4.5404,4.5405,4.5406,4.5406,
4.5407,4.5408,4.5409,4.541,4.541,4.5411,4.5412,4.5413,4.5414,4.5415,4.5415,4.5416,4.5417,
4.5418,4.5419,4.542,4.5421,4.5422,4.5423,4.5424,4.5425,4.5426,4.5427,4.5428,4.5429,4.543,
4.5431,4.5432,4.5433,4.5434,4.5435,4.5436,4.5437,4.5438,4.5439,4.544,4.5441,4.5443,4.5444,
4.5445,4.5446,4.5447,4.5449,4.545,4.5451,4.5453,4.5454,4.5455,4.5457,4.5458,4.5459,4.5461,
4.5462,4.5464,4.5465,4.5467,4.5469,4.547,4.5472,4.5474,4.5476,4.5477,4.5479,4.5481,4.5483,
4.5485,4.5487,4.549,4.5492,4.5494,4.5496,4.5499,4.5501,4.5504,4.5506,4.5509,4.5512,4.5514,
4.5517,4.552,4.5523,4.5526,4.553,4.5533,4.5536,4.554,4.5543,4.5547,4.5551,4.5555,4.5559,
4.5563,4.5567,4.5571,4.5576,4.558,4.5585,4.559,4.5594,4.5599,4.5605,4.561,4.5615,4.5621,
4.5627,4.5633,4.5639,4.5645,4.5651,4.5657,4.5664,4.5671,4.5678,4.5685,4.5692,4.5699,4.5707,
4.5715,4.5723,4.5731,4.5739,4.5747,4.5756,4.5764,4.5773,4.5782,4.5791,4.58,4.5809,4.5819,
4.5829,4.5839,4.5849,4.5859,4.5869,4.5879,4.589,4.5901,4.5912,4.5922,4.5934,4.5945,4.5956,
4.5967,4.5979,4.599,4.6002,4.6013,4.6025,4.6037,4.6049,4.606,4.6072,4.6084,4.6097,4.6109,
4.6121,4.6134,4.6147,4.616,4.6173,4.6186,4.62,4.6214,4.6228,4.6242,4.6256,4.627,4.6285,4.63,
4.6315,4.633,4.6345,4.6361,4.6376,4.6392,4.6407,4.6423,4.6439,4.6455,4.6472,4.6488,4.6504,
4.6521,4.6538,4.6554,4.6571,4.6588,4.6605,4.6622,4.6639,4.6656,4.6674,4.6691,4.6708,4.6726,
4.6743,4.6761,4.6778,4.6796,4.6814,4.6832,4.685,4.6868,4.6886,4.6904,4.6922,4.694,4.6958,
4.6976,4.6994,4.7012,4.703,4.7047,4.7065,4.7083,4.7101,4.7119,4.7136,4.7154,4.7172,4.7189,
4.7207,4.7225,4.7242,4.726,4.7277,4.7294,4.7311,4.7328,4.7345,4.7361,4.7378,4.7394,4.741,
4.7426,4.7442,4.7457,4.7473,4.7488,4.7503,4.7518,4.7533,4.7547,4.7562,4.7576,4.7589,4.7603,
4.7616,4.7629,4.7641,4.7654,4.7666,4.7678,4.769,4.7701,4.7712,4.7723,4.7734,4.7744,4.7754,
4.7764,4.7774,4.7784,4.7793,4.7802,4.7811,4.7819,4.7828,4.7836,4.7844,4.7852,4.7859,4.7867,
4.7874,4.7881,4.7888,4.7894,4.7901,4.7907,4.7913,4.7919,4.7925,4.7931,4.7936,4.7942,4.7947,
4.7953,4.7958,4.7963,4.7968,4.7973,4.7977,4.7982,4.7987,4.7991,4.7996,4.8,4.8004,4.8008,
4.8013,4.8017,4.8021,4.8025,4.8028,4.8032,4.8036,4.804,4.8044,4.8047,4.8051,4.8054,4.8058,
4.8061,4.8065,4.8068,4.8071,4.8074,4.8077,4.8081,4.8084,4.8087,4.809,4.8093,4.8096,4.8098,
4.8101,4.8104,4.8107,4.811,4.8113,4.8115,4.8118,4.8121,4.8123,4.8126,4.8128,4.8131,4.8133,
4.8136,4.8138,4.8141,4.8143,4.8146,4.8148,4.8151,4.8153,4.8155,4.8158,4.816,4.8162,4.8164,
4.8167,4.8169,4.8171,4.8173,4.8175,4.8178,4.818,4.8182,4.8184,4.8186,4.8188,4.819,4.8193,
4.8195,4.8197,4.8199,4.8201,4.8203,4.8206,4.8208,4.821,4.8212,4.8214,4.8216,4.8219,4.8221,
4.8223,4.8225,4.8227,4.823,4.8232,4.8234,4.8236,4.8239,4.8241,4.8243,4.8246,4.8248,4.825,
4.8253,4.8255,4.8258,4.826,4.8262,4.8265,4.8268,4.827,4.8273,4.8275,4.8278,4.8281,4.8284,
4.8287,4.829,4.8293,4.8296,4.8299,4.8302,4.8305,4.8309,4.8312,4.8315,4.8319,4.8323,4.8326,
4.833,4.8334,4.8338,4.8342,4.8346,4.835,4.8355,4.8359,4.8364,4.8369,4.8373,4.8378,4.8384,
4.8389,4.8394,4.84,4.8405,4.8411,4.8417,4.8423,4.843,4.8436,4.8442,4.8449,4.8456,4.8463,
4.847,4.8477,4.8484,4.8492,4.85,4.8507,4.8515,4.8523,4.8531,4.854,4.8548,4.8556,4.8565,
4.8573,4.8582,4.8591,4.8599,4.8608,4.8617,4.8627,4.8636,4.8645,4.8654,4.8664,4.8674,4.8683,
4.8693,4.8703,4.8713,4.8724,4.8734,4.8745,4.8756,4.8766,4.8778,4.8789,4.88,4.8812,4.8824,
4.8836,4.8848,4.886,4.8872,4.8885,4.8898,4.8911,4.8924,4.8937,4.895,4.8964,4.8977,4.8991,
4.9005,4.9019,4.9033,4.9047,4.9061,4.9075,4.9089,4.9104,4.9118,4.9133,4.9147,4.9162,4.9177,
4.9192,4.9206,4.9221,4.9236,4.9251,4.9266,4.9281,4.9295,4.931,4.9325,4.934,4.9354,4.9369,
4.9384,4.9398,4.9413,4.9427,4.9442,4.9456,4.9471,4.9485,4.9499,4.9514,4.9528,4.9542,4.9556,
4.957,4.9585,4.9599,4.9613,4.9627,4.964,4.9654,4.9668,4.9681,4.9695,4.9708,4.9721,4.9734,
4.9747,4.976,4.9773,4.9785,4.9797,4.9809,4.9821,4.9833,4.9845,4.9856,4.9868,4.9879,4.989,
4.9901,4.9912,4.9923,4.9933,4.9944,4.9954,4.9964,4.9974,4.9984,4.9994,5.0004,5.0013,5.0022,
5.0031,5.004,5.0049,5.0058,5.0066,5.0074,5.0083,5.0091,5.0098,5.0106,5.0114,5.0121,5.0128,
5.0135,5.0142,5.0149,5.0155,5.0162,5.0168,5.0174,5.018,5.0186,5.0192,5.0197,5.0203,5.0208,
5.0213,5.0218,5.0223,5.0228,5.0232,5.0237,5.0241,5.0246,5.025,5.0254,5.0258,5.0262,5.0266,
5.0269,5.0273,5.0277,5.028,5.0284,5.0287,5.029,5.0294,5.0297,5.03,5.0303,5.0306,5.0309,
5.0312,5.0315,5.0318,5.0321,5.0324,5.0327,5.0329,5.0332,5.0335,5.0337,5.034,5.0342,5.0345,
5.0347,5.0349,5.0352,5.0354,5.0356,5.0358,5.0361,5.0363,5.0365,5.0367,5.0369,5.0371,5.0373,
5.0375,5.0376,5.0378,5.038,5.0382,5.0383,5.0385,5.0387,5.0388,5.039,5.0392,5.0393,5.0395,
5.0396,5.0397,5.0399,5.04,5.0402,5.0403,5.0404,5.0405,5.0407,5.0408,5.0409,5.041,5.0411,
5.0412,5.0413,5.0415,5.0416,5.0417,5.0418,5.0418,5.0419,5.042,5.0421,5.0422,5.0423,5.0424,
5.0424,5.0425,5.0426,5.0427,5.0427,5.0428,5.0429,5.0429,5.043,5.0431,5.0431,5.0432,5.0432,
5.0433,5.0433,5.0434,5.0435,5.0435,5.0435,5.0436,5.0436,5.0437,5.0437,5.0438,5.0438,5.0438,
5.0439,5.0439,5.0439,5.044,5.044,5.044,5.0441,5.0441,5.0441,5.0442,5.0442,5.0442,5.0442,
5.0443,5.0443,5.0443,5.0443,5.0443,5.0444,5.0444,5.0444,5.0444,5.0444,5.0445,5.0445,5.0445,
5.0445,5.0445,5.0445,5.0446,5.0446,5.0446,5.0446,5.0446,5.0446,5.0446,5.0446,5.0446,5.0447,
5.0447,5.0447,5.0447,5.0447,5.0447,5.0447,5.0447,5.0447,5.0447,5.0447,5.0448,5.0448,5.0448,
5.0448,5.0448,5.0448,5.0448"""

NEARCRIT_OOB = """\
5.16,5.3,4.29,3.76,3.52,3.35,3.19,3.04,2.89,2.75,2.62,2.49,2.36,2.24,2.12,2.04,1.98,1.93,
1.87,1.82,1.77,1.72,1.68,1.63,1.59,1.55,1.5,1.47,1.43,1.39,1.35,1.3,1.25,1.2,1.16,1.12,1.07,
1.04,0.997,0.961,0.927,0.9,0.878,0.856,0.834,0.813,0.793,0.773,0.754,0.735,0.717,0.7,0.683,
0.666,0.65,0.635,0.62,0.605,0.591,0.578,0.565,0.552,0.54,0.529,0.517,0.506,0.496,0.486,0.476,
0.466,0.457,0.448,0.44,0.432,0.424,0.416,0.41,0.403,0.397,0.391,0.385,0.38,0.376,0.371,0.368,
0.364,0.361,0.358,0.355,0.352,0.35,0.348,0.347,0.345,0.344,0.343,0.343,0.342,0.342,0.341,
0.341,0.341,0.341,0.341,0.341,0.341,0.341,0.342,0.342,0.342,0.342,0.343,0.343,0.343,0.343,
0.344,0.344,0.344,0.344,0.344,0.344,0.344,0.344,0.344,0.344,0.343,0.343,0.342,0.341,0.34,
0.339,0.338,0.337,0.335,0.333,0.332,0.33,0.328,0.326,0.324,0.321,0.319,0.316,0.314,0.311,
0.309,0.306,0.304,0.301,0.298,0.295,0.292,0.289,0.286,0.283,0.28,0.276,0.273,0.269,0.266,
0.262,0.258,0.254,0.25,0.246,0.242,0.237,0.233,0.228,0.224,0.219,0.214,0.209,0.204,0.199,
0.194,0.189,0.184,0.179,0.173,0.168,0.163,0.157,0.152,0.146,0.14,0.134,0.128,0.123,0.118,
0.113,0.108,0.103,0.0987,0.0945,0.0906,0.0869,0.0833,0.0799,0.0766,0.0735,0.0705,0.0676,
0.0649,0.0624,0.0599,0.0576,0.0554,0.0533,0.0513,0.0494,0.0475,0.0458,0.0441,0.0425,0.041,
0.0396,0.0382,0.037,0.0357,0.0346,0.0335,0.0324,0.0314,0.0304,0.0298,0.0295,0.0292,0.0289,
0.0285,0.0281,0.0277,0.0273,0.0269,0.0265,0.026,0.0255,0.0249,0.0243,0.0238,0.0232,0.0227,
0.0222,0.0217,0.0211,0.0205,0.0197,0.0192,0.0189,0.0186,0.0181,0.0175,0.0165,0.0158,0.0152,
0.0147,0.0142,0.0137,0.0132,0.0126,0.0121,0.0116,0.0107,0.0104,0.0103,0.0103,0.0103,0.0103,
0.0103,0.0104,0.0105,0.0107,0.0109,0.0111,0.0114,0.0119,0.0123,0.0127,0.0132,0.0136,0.0141,
0.0146,0.0152,0.0158,0.0164,0.017,0.0175,0.0181,0.0186,0.0193,0.02,0.0206,0.0213,0.0219,
0.0226,0.0234,0.0242,0.025,0.0258,0.0266,0.0276,0.0285,0.0295,0.0305,0.0315,0.0326,0.0338,
0.0349,0.036,0.0374,0.0387,0.04,0.0414,0.0428,0.0443,0.0459,0.0476,0.0493,0.0509,0.0527,
0.0545,0.0565,0.0585,0.0603,0.0622,0.0642,0.0661,0.0681,0.07,0.0719,0.0739,0.076,0.0783,
0.0806,0.0829,0.0854,0.0879,0.0904,0.0931,0.0955,0.0982,0.101,0.103,0.105,0.108,0.11,0.113,
0.115,0.117,0.119,0.122,0.124,0.126,0.127,0.129,0.131,0.133,0.135,0.137,0.139,0.141,0.143,
0.145,0.147,0.149,0.152,0.153,0.156,0.157,0.16,0.162,0.164,0.166,0.168,0.17,0.172,0.175,
0.176,0.178,0.179,0.181,0.183,0.185,0.186,0.188,0.19,0.192,0.193,0.195,0.197,0.198,0.2,0.201,
0.203,0.204,0.205,0.207,0.208,0.209,0.211,0.212,0.213,0.214,0.215,0.216,0.216,0.217,0.218,
0.218,0.219,0.219,0.219,0.22,0.22,0.22,0.221,0.221,0.221,0.221,0.221,0.222,0.222,0.222,0.222,
0.223,0.223,0.223,0.224,0.224,0.225,0.226,0.226,0.227,0.227,0.228,0.229,0.229,0.23,0.23,
0.231,0.231,0.232,0.232,0.232,0.232,0.232,0.232,0.232,0.232,0.231,0.231,0.231,0.23,0.23,
0.229,0.228,0.228,0.227,0.226,0.225,0.224,0.223,0.222,0.221,0.22,0.219,0.217,0.216,0.215,
0.214,0.212,0.211,0.21,0.208,0.207,0.205,0.204,0.202,0.2,0.198,0.196,0.194,0.193,0.191,0.189,
0.187,0.184,0.182,0.18,0.178,0.175,0.172,0.17,0.167,0.164,0.162,0.158,0.155,0.152,0.149,
0.146,0.143,0.14,0.137,0.134,0.131,0.129,0.126,0.123,0.121,0.118,0.116,0.114,0.112,0.109,
0.107,0.105,0.103,0.102,0.101,0.0992,0.0976,0.0961,0.0948,0.0941,0.0929,0.0917,0.0905,0.0893,
0.0888,0.088,0.0871,0.0861,0.0852,0.0844,0.0843,0.0838,0.0831,0.0825,0.0819,0.0814,0.0815,
0.0815,0.0812,0.081,0.0808,0.0806,0.0809,0.0815,0.0818,0.082,0.0823,0.0827,0.0832,0.0845,
0.0855,0.0864,0.0873,0.0882,0.0893,0.0908,0.0927,0.0942,0.0956,0.0971,0.0987,0.1,0.103,0.105,
0.106,0.108,0.11,0.112,0.115,0.117,0.12,0.122,0.125,0.127,0.13,0.133,0.136,0.139,0.143,0.147,
0.151,0.154,0.158,0.163,0.167,0.172,0.176,0.18,0.185,0.19,0.195,0.2,0.205,0.211,0.216,0.221,
0.226,0.232,0.237,0.242,0.248,0.253,0.259,0.264,0.269,0.274,0.28,0.284,0.289,0.293,0.298,
0.302,0.306,0.309,0.313,0.316,0.318,0.32,0.323,0.325,0.326,0.328,0.33,0.331,0.331,0.331,
0.332,0.332,0.331,0.33,0.329,0.328,0.326,0.324,0.323,0.321,0.319,0.316,0.315,0.312,0.309,
0.307,0.304,0.301,0.298,0.295,0.292,0.288,0.285,0.282,0.278,0.274,0.271,0.267,0.263,0.259,
0.255,0.251,0.247,0.242,0.238,0.233,0.228,0.223,0.218,0.213,0.207,0.201,0.195,0.189,0.184,
0.178,0.173,0.168,0.163,0.158,0.153,0.149,0.144,0.14,0.136,0.133,0.129,0.125,0.122,0.119,
0.115,0.112,0.109,0.106,0.103,0.101,0.0979,0.0953,0.0927,0.0903,0.0878,0.0854,0.0831,0.0809,
0.0786,0.0765,0.0745,0.0725,0.0706,0.0688,0.0669,0.0652,0.0635,0.0618,0.0602,0.0587,0.0573,
0.0559,0.0545,0.0533,0.052,0.0509,0.0498,0.0487,0.0477,0.0467,0.0458,0.0449,0.0441,0.0433,
0.0425,0.0418,0.0411,0.0404,0.0397,0.0391,0.0385,0.038,0.0374,0.0369,0.0364,0.036,0.0355,
0.0352,0.0348,0.0345,0.0341,0.0339,0.0336,0.0334,0.0331,0.033,0.0328,0.0326,0.0325,0.0324,
0.0323,0.0322,0.0321,0.0321,0.032,0.032,0.032,0.0321,0.0321,0.0322,0.0322,0.0323,0.0324,
0.0326,0.0327,0.0329,0.033,0.0332,0.0334,0.0336,0.0338,0.034,0.0342,0.0345,0.0347,0.035,
0.0353,0.0356,0.0358,0.0362,0.0365,0.0368,0.0372,0.0375,0.0379,0.0383,0.0386,0.039,0.0395,
0.0399,0.0404,0.041,0.0416,0.0422,0.0428,0.0435,0.0441,0.0449,0.0456,0.0464,0.0472,0.048,
0.0488,0.0497,0.0506,0.0515,0.0525,0.0535,0.0546,0.0557,0.0568,0.058,0.0593,0.0606,0.0618,
0.0632,0.0644,0.0658,0.0672,0.0688,0.0703,0.0719,0.0735,0.0751,0.0768,0.0785,0.0802,0.082,
0.0838,0.0857,0.0875,0.0895,0.0914,0.0934,0.0954,0.0974,0.0994,0.101,0.103,0.105,0.107,0.109,
0.111,0.113,0.115,0.117,0.119,0.121,0.122,0.124,0.126,0.129,0.131,0.133,0.135,0.138,0.14,
0.143,0.145,0.148,0.151,0.153,0.156,0.159,0.162,0.165,0.167,0.17,0.172,0.175,0.177,0.18,
0.182,0.184,0.187,0.189,0.191,0.193,0.195,0.197,0.199,0.201,0.202,0.204,0.206,0.208,0.209,
0.211,0.212,0.214,0.215,0.216,0.218,0.219,0.219,0.22,0.221,0.222,0.223,0.223,0.224,0.224,
0.225,0.225,0.226,0.226,0.226,0.226,0.226,0.226,0.225,0.225,0.224,0.224,0.223,0.222,0.222,
0.221,0.22,0.219,0.218,0.217,0.216,0.215,0.213,0.212,0.211,0.21,0.208,0.207,0.205,0.203,
0.202,0.2,0.198,0.196,0.193,0.191,0.189,0.187,0.184,0.182,0.18,0.178,0.175,0.173,0.171,0.169,
0.166,0.164,0.161,0.159,0.156,0.154,0.151,0.149,0.146,0.143,0.141,0.138,0.135,0.133,0.13,
0.127,0.125,0.122,0.119,0.117,0.114,0.112,0.109,0.107,0.104,0.102,0.0992,0.0968,0.0944,
0.0921,0.0898,0.0876,0.0854,0.0832,0.0811,0.0791,0.077,0.0751,0.0731,0.0713,0.0695,0.0677,
0.066,0.0643,0.0626,0.061,0.0595,0.0579,0.0564,0.055,0.0535,0.0522,0.0508,0.0495,0.0482,
0.0469,0.0457,0.0445,0.0433,0.0422,0.041,0.04,0.0389,0.0379,0.037,0.0361,0.0352,0.0343,
0.0334,0.0326,0.0318,0.0311,0.0303,0.0296,0.0289,0.0282,0.0275,0.0269,0.0262,0.0256,0.025,
0.0244,0.0239,0.0233,0.0228,0.0223,0.0218,0.0214,0.0209,0.0205,0.02,0.0196,0.0192,0.0189,
0.0185,0.0181,0.0178,0.0175,0.0172,0.0169,0.0166,0.0163,0.016,0.0157,0.0155,0.0152,0.015,
0.0148,0.0145,0.0143,0.0141,0.0138,0.0136,0.0134,0.0132,0.013,0.0128,0.0126,0.0124,0.0123,
0.0121,0.0119,0.0117,0.0116,0.0114,0.0113,0.0111,0.011,0.0108,0.0107,0.0106,0.0105,0.0103,
0.0102,0.0101,0.00996,0.00983,0.00972,0.0096,0.00949,0.00937,0.00927,0.00916,0.00906,0.00895,
0.00886,0.00875,0.00866,0.00855,0.00847,0.00837,0.00829,0.00821,0.00814,0.00806,0.00799,
0.00792,0.00786,0.0078,0.00775,0.00768,0.00764,0.00758,0.00753,0.00748,0.00744,0.00738,
0.00735,0.0073,0.00727,0.00722,0.00719,0.00715,0.00712,0.00708,0.00705,0.00701,0.00698,
0.00695,0.00692,0.00689,0.00686,0.00683,0.00681,0.00677,0.00675,0.00672,0.0067,0.00667,
0.00665,0.00661,0.00659,0.00655,0.00653,0.00649,0.00647,0.00643,0.00641,0.00637,0.00635,
0.00632,0.0063,0.00627,0.00625,0.00622,0.00621,0.00618,0.00616,0.00614,0.00613,0.0061,
0.00609,0.00607,0.00606,0.00604,0.00603,0.00601,0.00601,0.00599,0.00599,0.00597,0.00597,
0.00596,0.00596,0.00594,0.00594,0.00593,0.00594,0.00593,0.00593,0.00593,0.00594,0.00593,
0.00595,0.00595,0.00596,0.00597,0.00599,0.00599,0.00601,0.00602,0.00604,0.00605,0.00608,
0.00609,0.00612,0.00613,0.00616,0.00618,0.00621,0.00623,0.00627,0.00629,0.00633,0.00636,
0.0064,0.00643,0.00647,0.0065,0.00655,0.00658,0.00662,0.00666,0.0067,0.00674,0.00679,0.00682,
0.00688,0.00692,0.00698,0.00702,0.00708,0.00713,0.0072,0.00725,0.00732,0.00737,0.00745,
0.0075,0.00758,0.00764,0.00771,0.00777,0.00785,0.00792,0.008,0.00806,0.00814,0.00821,0.00829,
0.00836,0.00844,0.00849,0.00858,0.00865,0.00874,0.00882,0.00891,0.00899,0.00909,0.00916,
0.00926,0.00934,0.00944,0.00952,0.00962,0.00971,0.00981,0.0099,0.01,0.0101,0.0102,0.0103,
0.0104,0.0106,0.0107,0.0108,0.011,0.0111,0.0112,0.0113,0.0115,0.0116,0.0117,0.0119,0.012,
0.0122,0.0123,0.0125,0.0127,0.0129,0.0131,0.0133,0.0135,0.0137,0.0139,0.0142,0.0145,0.0147,
0.015,0.0153,0.0157,0.016,0.0164,0.0167,0.0172,0.0175,0.018,0.0184,0.0189,0.0193,0.0198,
0.0203,0.0208,0.0214,0.0219,0.0225,0.0231,0.0237,0.0243,0.0249,0.0256,0.0262,0.027,0.0277,
0.0284,0.0292,0.0299,0.0307,0.0316,0.0324,0.0332,0.0341,0.035,0.0359,0.0369,0.0378,0.0389,
0.0399,0.041,0.042,0.0432,0.0443,0.0455,0.0467,0.0479,0.0491,0.0504,0.0517,0.053,0.0542,
0.0555,0.0568,0.0583,0.0597,0.0612,0.0627,0.0641,0.0656,0.0672,0.0689,0.0706,0.0724,0.0742,
0.0759,0.0778,0.0796,0.0814,0.0833,0.0851,0.087,0.089,0.091,0.093,0.0949,0.0969,0.0989,0.101,
0.103,0.105,0.107,0.109,0.111,0.112,0.114,0.116,0.118,0.12,0.121,0.123,0.125,0.127,0.128,
0.13,0.132,0.134,0.135,0.137,0.139,0.141,0.143,0.145,0.147,0.15,0.152,0.154,0.157,0.159,
0.161,0.164,0.166,0.169,0.171,0.174,0.176,0.179,0.181,0.183,0.185,0.188,0.19,0.192,0.194,
0.196,0.198,0.2,0.202,0.204,0.205,0.207,0.209,0.21,0.212,0.213,0.215,0.216,0.217,0.219,0.22,
0.221,0.222,0.224,0.225,0.226,0.227,0.228,0.229,0.229,0.23,0.231,0.231,0.232,0.233,0.233,
0.234,0.234,0.235,0.235,0.236,0.236,0.236,0.236,0.236,0.237,0.237,0.237,0.237,0.237,0.237,
0.237,0.236,0.236,0.236,0.235,0.234,0.234,0.233,0.232,0.231,0.23,0.229,0.228,0.226,0.225,
0.224,0.222,0.22,0.219,0.217,0.214,0.212,0.209,0.206,0.203,0.199,0.196,0.193,0.189,0.186,
0.183,0.179,0.176,0.172,0.169,0.165,0.162,0.158,0.155,0.151,0.148,0.145,0.141,0.138,0.135,
0.131,0.128,0.125,0.122,0.118,0.115,0.112,0.11,0.107,0.104,0.101,0.0989,0.0964,0.0941,0.0917,
0.0895,0.0873,0.0852,0.0831,0.0812,0.0792,0.0773,0.0754,0.0737,0.0719,0.0702,0.0686,0.067,
0.0655,0.064,0.0625,0.0611,0.0598,0.0585,0.0572,0.056,0.0548,0.0536,0.0525,0.0514,0.0503,
0.0492,0.0482,0.0471,0.0461,0.0451,0.0442,0.0433,0.0424,0.0415,0.0407,0.0399,0.0391,0.0384,
0.0377,0.037,0.0364,0.0357,0.0351,0.0345,0.0339,0.0334,0.0329,0.0323,0.0318,0.0314,0.0309,
0.0304,0.03,0.0296,0.0292,0.0288,0.0284,0.028,0.0277,0.0273,0.027,0.0267,0.0263,0.026,0.0257,
0.0254,0.0251,0.0248,0.0246,0.0243,0.0241,0.0239,0.0237,0.0235,0.0234,0.0233,0.0231,0.023,
0.0229,0.0229,0.0228,0.0228,0.0227,0.0227,0.0227,0.0227,0.0226,0.0227,0.0227,0.0227,0.0227,
0.0227,0.0228,0.0228,0.0229,0.023,0.023,0.0231,0.0232,0.0233,0.0234,0.0235,0.0236,0.0238,
0.0239,0.0241,0.0243,0.0245,0.0247,0.0249,0.0252,0.0255,0.0259,0.0262,0.0266,0.027,0.0274,
0.0279,0.0283,0.0288,0.0293,0.0298,0.0303,0.0309,0.0315,0.0322,0.0328,0.0335,0.0342,0.035,
0.0358,0.0366,0.0375,0.0384,0.0393,0.0403,0.0413,0.0424,0.0434,0.0445,0.0456,0.0468,0.048,
0.0493,0.0505,0.0519,0.0532,0.0546,0.056,0.0574,0.0588,0.0602,0.0616,0.0631,0.0645,0.0661,
0.0676,0.0691,0.0705,0.072,0.0735,0.0751,0.0767,0.0784,0.08,0.0817,0.0833,0.085,0.0866,
0.0883,0.09,0.0917,0.0934,0.0951,0.0968,0.0985,0.1,0.102,0.104,0.106,0.107,0.109,0.111,0.113,
0.115,0.117,0.119,0.121,0.123,0.126,0.128,0.13,0.132,0.135,0.137,0.14,0.142,0.145,0.147,0.15,
0.152,0.155,0.157,0.16,0.162,0.164,0.167,0.169,0.171,0.173,0.176,0.178,0.179,0.181,0.183,
0.185,0.187,0.188,0.19,0.192,0.193,0.195,0.196,0.197,0.199,0.2,0.201,0.201,0.202,0.203,0.204,
0.204,0.205,0.205,0.206,0.206,0.206,0.207,0.207,0.207,0.207,0.207,0.206,0.206,0.206,0.206,
0.206,0.206,0.205,0.205,0.205,0.204,0.204,0.203,0.202,0.202,0.201,0.2,0.199,0.198,0.197,
0.195,0.194,0.193,0.191,0.19,0.188,0.186,0.185,0.183,0.181,0.18,0.178,0.176,0.174,0.173,
0.171,0.169,0.167,0.165,0.163,0.161,0.159,0.157,0.154,0.152,0.15,0.148,0.145,0.143,0.141,
0.139,0.136,0.134,0.132,0.129,0.127,0.124,0.122,0.12,0.117,0.115,0.113,0.11,0.108,0.106,
0.103,0.101,0.0986,0.0964,0.0941,0.0919,0.0896,0.0875,0.0853,0.0832,0.0811,0.0791,0.077,
0.0751,0.0732,0.0714,0.0696,0.0679,0.0663,0.0647,0.0632,0.0618,0.0603,0.059,0.0577,0.0564,
0.0552,0.054,0.0528,0.0517,0.0505,0.0494,0.0483,0.0472,0.0462,0.0452,0.0442,0.0432,0.0423,
0.0414,0.0404,0.0396,0.0387,0.0379,0.0371,0.0363,0.0354,0.0347,0.0339,0.0332,0.0324,0.0317,
0.031,0.0303,0.0296,0.029,0.0283,0.0276,0.027,0.0263,0.0257,0.0251,0.0245,0.0239,0.0233,
0.0228,0.0222,0.0217,0.0211,0.0206,0.0201,0.0196,0.0191,0.0186,0.0181,0.0177,0.0172,0.0168,
0.0164,0.016,0.0155,0.0151,0.0147,0.0144,0.014,0.0136,0.0132,0.0129,0.0125,0.0122,0.0119,
0.0115,0.0112,0.0109,0.0106,0.0103,0.00997,0.00968,0.00939,0.00911,0.00883,0.00857,0.00829,
0.00804,0.00778,0.00755,0.0073,0.00708,0.00684,0.00663,0.00641,0.00621,0.006,0.00581,0.00561,
0.00543,0.00524,0.00507,0.00489,0.00473,0.00456,0.0044,0.00425,0.0041,0.00395,0.0038,0.00365,
0.00351,0.00337,0.00324,0.00311,0.00299,0.00287,0.00275,0.00264,0.00253,0.00246,0.00238,
0.00233,0.00227,0.00222,0.00215,0.00211,0.00205,0.002,0.00194,0.0019,0.00184,0.0018,0.00175,
0.00171,0.00166,0.00163,0.00158,0.00154,0.0015,0.00146,0.00142,0.00138,0.00135,0.00133,
0.00128,0.00129,0.00122,0.00125,0.00116,0.00121,0.00111,0.00119,0.00109,0.00116,0.00107,
0.00114,0.00105,0.00111,0.00103,0.00108,0.001"""

if __name__ == '__main__':
    sys.exit(main())
