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

    # The early exit is NOT itself a verdict: a no-progress exit on a frozen state
    # is still STABLE_STUCK (this is the premise correction the hybrid exists for).
    v, _, _ = classify_nonconvergence([1.04] * 40, 1.0, 'no_progress')
    check("no-progress early exit on a frozen state -> STABLE_STUCK",
          v == 'STABLE_STUCK')

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

    # The criterion is threaded down to every trial, not just to the driver.
    _, stub = _bisect(kinds, hybrid=True)
    check("every trial is solved with failure_criterion='hybrid'",
          bool(stub.calls) and all(c[2] == 'hybrid' for c in stub.calls))
    _, stub = _bisect(kinds, hybrid=False)
    check("...and with 'non_convergence' on the default",
          bool(stub.calls) and all(c[2] == 'non_convergence' for c in stub.calls))


# ===================== 3. default invariance on a real solve =====================

def check_real_solve():
    print("\n3. real solve — metadata is recorded, default behavior is unchanged")

    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, solve_fem
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    xlsx = os.path.join(here, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')
    if not os.path.exists(xlsx):
        check("griffiths1 input present", False, xlsx)
        return

    slope_data = load_slope_data(xlsx)
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=8.0, element_type='tri6')
    fem_data = build_fem_data(slope_data, mesh)

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


def main():
    print("=" * 72)
    print("Hybrid failure-criterion checks")
    print("=" * 72)
    check_classifier()
    check_wiring()
    check_real_solve()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All hybrid failure-criterion checks passed.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
