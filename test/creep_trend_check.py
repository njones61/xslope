"""Checks for the TREND READING — how a trial still moving at its iteration limit
is read (`xslope.fem.creep_trend`, `creep_extrapolate`, the 'not_slowing' and
'slowing' endings and their sentences in the closing summary).

What this file locks:

  1. THE CLASSES. Synthetic block-end records: a movement shrinking at a steady
     ratio is 'dying'; one holding its pace is 'steady'; one speeding up is
     'growing'; one below the still line is 'still'; one that zigzags is
     'unclear'. On a jointed model a slip that does not slow makes the trial
     'steady' whatever the displacement increments say. Too few blocks, or no
     elastic yardstick, and the reading declines. The levels are the file's
     existing ones (the joint verdict's 0.9 rate ratio, the classifier's 0.02
     growth line, the joint verdict's 1e-4 "the field has stopped").

  2. A RATIO, NOT A COUNT. The same movement read on blocks of different length
     (a plain sweep and an accelerated one covering it in fewer iterations) is
     read the same way.

  3. THE EXTRAPOLATION. A field decaying geometrically is carried to its limit,
     by the field's one ratio and per degree of freedom (Aitken), and a movement
     that is not shrinking gives nothing to extrapolate.

  4. THE RECORDED TRACES (test/fixtures/creep_traces.json). The FEM-3 geogrid
     wall's trial at F = 1.25, at its page settings, is dying away at 100,000
     iterations; the FEM-3 wall alone at F = 1.140625 is holding steady at
     90,000 and is counted as sliding.

  5. THE SENTENCES. The closing summary quotes a standing edge counted standing
     by a slowing movement (from where it was, or from the estimated resting
     state), a failing edge counted as sliding by a movement that did not slow,
     and one still slowing at the limit that the corrector could not finish, in
     the Run dialog's words.

  6. THE SEEDS (minutes; skipped with --quick). The FEM-3 geogrid wall at
     F = 1.25, at its page settings: at 100,000 iterations the movement is dying
     away and the corrector refuses both seeds, the block-end state and the
     estimated resting state. At a 200,000 allowance each seed on its own
     certifies at 140,000 and the hold test holds it, at about 1.30 elastic
     displacements, where the plain sweep's own rest (465,581 iterations) is
     1.297. It is the corrector attempt at every block end after the warm-up that
     stands this trial; the extrapolation adds nothing measurable on it.

Run directly:  PYTHONPATH=. python3 test/creep_trend_check.py [--quick]
Exits non-zero on any failure.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402
from xslope.fem import creep_trend, creep_extrapolate      # noqa: E402

FAILURES = []
TRACES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'fixtures', 'creep_traces.json')


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


def _marks(first, ratio, n=6, start=1.0):
    """Block-end max|u| for a movement of ``first`` in the first block, each
    block moving ``ratio`` times the one before."""
    out = [start]
    step = first
    for _ in range(n - 1):
        out.append(out[-1] + step)
        step *= ratio
    return out


def check_classes():
    print("\n1. the classes")
    ue = 1.0
    r = creep_trend(_marks(0.02, 0.7), ue, 10000)
    check("a movement shrinking at 0.7 a block is dying away",
          r and r['trend'] == 'dying', f"{r and r['trend']}, ratio {r and r['ratio']:.3f}")
    check("its ratio is the geometric mean of the block ratios",
          r and abs(r['ratio'] - 0.7) < 1e-9)
    r = creep_trend(_marks(0.02, 1.0), ue, 10000)
    check("a movement holding its pace is steady", r and r['trend'] == 'steady',
          r and r['trend'])
    r = creep_trend(_marks(0.02, 0.95), ue, 10000)
    check("one slowing by only 5% a block is steady too (0.9 is the line)",
          r and r['trend'] == 'steady', r and r['trend'])
    r = creep_trend(_marks(0.02, 1.3), ue, 10000)
    check("a movement speeding up is growing", r and r['trend'] == 'growing',
          r and r['trend'])
    r = creep_trend(_marks(1e-6, 1.0), ue, 10000)
    check("a movement below the still line is still", r and r['trend'] == 'still',
          r and r['trend'])
    r = creep_trend(_marks(0.001, 1.0), ue, 10000)
    check("a steady movement under the 0.02 moving line is not called sliding",
          r and r['trend'] == 'unclear', r and r['trend'])
    zig = [1.0, 1.02, 1.021, 1.04, 1.041, 1.046]
    r = creep_trend(zig, ue, 10000)
    check("a zigzag that ends small is unclear (a block moved more than the one "
          "before)", r and r['trend'] == 'unclear', r and r['trend'])
    zig = [1.0, 1.02, 1.021, 1.04, 1.041, 1.06]
    r = creep_trend(zig, ue, 10000)
    check("a zigzag that keeps its pace is steady", r and r['trend'] == 'steady',
          r and r['trend'])
    back = [1.0, 1.02, 1.03, 1.035, 1.034, 1.0345]
    r = creep_trend(back, ue, 10000)
    check("a block that moved back is not dying away",
          r and r['trend'] != 'dying', r and r['trend'])
    # joints: the slip reading, 50,000 iterations at 10 a sample
    n = 5000
    steady_slip = list(np.linspace(1.0, 2.0, 2 * n))
    r = creep_trend(_marks(0.02, 0.7), ue, 10000, slip_hist=steady_slip)
    check("a slip that does not slow makes a dying displacement steady",
          r and r['trend'] == 'steady' and r['slip_ratio'] >= 0.9,
          f"{r and r['trend']}, slip ratio {r and r['slip_ratio']}")
    t = np.arange(2 * n, dtype=float)
    dying_slip = list(2.0 - np.exp(-t / 1500.0))
    r = creep_trend(_marks(0.02, 0.7), ue, 10000, slip_hist=dying_slip)
    check("a slip dying away with it leaves it dying away",
          r and r['trend'] == 'dying', f"{r and r['trend']}, {r and r['slip_ratio']}")
    check("too few blocks: no reading",
          creep_trend(_marks(0.02, 0.7, n=5), ue, 10000) is None)
    check("no elastic yardstick: no reading",
          creep_trend(_marks(0.02, 0.7), 0.0, 10000) is None)
    check("the levels are the file's own",
          fem._CREEP_DYING_MAX == fem._JOINT_MOVING_DECAY_MIN
          and fem._CREEP_MOVING == fem._HYBRID_GROWTH_MIN
          and fem._CREEP_STILL == fem._JOINT_SETTLED_GROWTH)


def check_ratio_not_count():
    print("\n2. a ratio of movements, not a count of iterations")
    # u(t) = 1 - exp(-t / T): a plain sweep reads it on 10,000-iteration blocks,
    # an accelerated one covering the same ground in a third of the iterations
    # reads it on 3,333-iteration blocks of the same movement.
    T = 40000.0
    plain = [1.0 - np.exp(-(50000 + 10000 * k) / T) for k in range(6)]
    fast = [1.0 - np.exp(-(50000 + 10000 * k) / T) for k in range(6)]
    a = creep_trend(plain, 0.01, 10000)
    b = creep_trend(fast, 0.01, 3333)
    check("the same movement on blocks of a different length reads the same",
          a['trend'] == b['trend'] == 'dying' and abs(a['ratio'] - b['ratio']) < 1e-12,
          f"{a['trend']} {a['ratio']:.4f} / {b['trend']} {b['ratio']:.4f}")


def check_extrapolation():
    print("\n3. the extrapolation")
    rng = np.random.default_rng(3)
    lim = rng.normal(size=200)
    a = rng.normal(size=200)
    snaps = [lim - a * 0.8 ** k for k in (4, 5, 6)]
    u, r = creep_extrapolate(snaps, per_component=False)
    check("one ratio: a geometric decay is carried to its limit",
          u is not None and np.abs(u - lim).max() < 1e-10 and abs(r - 0.8) < 1e-12,
          f"error {np.abs(u - lim).max():.1e}, r {r:.4f}")
    # two modes decaying at different rates: Aitken per component reads each
    b = rng.normal(size=200)
    mask = np.arange(200) < 100
    snaps = [lim - np.where(mask, a * 0.6 ** k, b * 0.9 ** k) for k in (4, 5, 6)]
    u1, _ = creep_extrapolate(snaps, per_component=False)
    u2, _ = creep_extrapolate(snaps, per_component=True)
    check("per component, two modes at two rates are each carried to their limit",
          np.abs(u2 - lim).max() < 1e-9 < np.abs(u1 - lim).max(),
          f"Aitken {np.abs(u2 - lim).max():.1e}, one ratio {np.abs(u1 - lim).max():.1e}")
    snaps = [lim + a * k for k in (0, 1, 2)]
    u, r = creep_extrapolate(snaps)
    check("a movement that is not shrinking has nowhere to be carried",
          u is None, f"r {r}")


def check_recorded():
    print("\n4. the recorded traces")
    with open(TRACES) as f:
        rec = json.load(f)
    for name, want, sweeps in (("grid_F1.25", ("dying",), 100000),
                               ("blocks_F1.140625", ("steady", "growing"), 90000)):
        t = rec[name]
        se = int(t["sample_every"])
        block = 10000
        idx = [(sweeps - block * k) // se for k in range(5, -1, -1)]
        marks = [t["disp"][min(i, len(t["disp"]) - 1)] for i in idx]
        slip = t["slip"][: sweeps // se]
        r = creep_trend(marks, t["u_elastic_scale"], block, slip_hist=slip,
                        sample_every=se)
        check(f"{name}: read {'/'.join(want)} at {sweeps:,} iterations",
              r is not None and r["trend"] in want,
              f"{r and r['trend']}, ratio {r and r['ratio']:.3f}, moved "
              f"{r and r['moved']:.4f} u_el, slip ratio {r and r['slip_ratio']}")
        check(f"{name}: the solver's own reading said the same",
              (t.get("reading") or {}).get("trend") in want,
              (t.get("reading") or {}).get("trend"))


def _run(trials, lo, hi):
    return {"converged": True, "FS": 0.5 * (lo + hi), "final_interval": [lo, hi],
            "elapsed_time": 600.0, "trials": trials}


def check_sentences():
    print("\n5. the sentences")
    slowing = dict(rule='slowing', window=50000, block=10000, iteration=100000,
                   increments=[0.02, 0.016, 0.0128, 0.0102, 0.0082], ratio=0.8,
                   extrapolated=0.0006, max_displacement=0.023,
                   seed='extrapolated',
                   corrector={'certified': True, 'hold': {'held': True}})
    sliding = dict(rule='not_slowing', window=50000, block=10000,
                   iteration=90000, ratio=0.96,
                   increments=[0.01, 0.0096, 0.0092, 0.0088, 0.0085])
    trials = [
        {"F": 1.25, "stable": True, "converged": True, "verdict": "CONVERGED",
         "iterations": 100011, "exit_reason": "converged",
         "corrector": {"checkpoint": "trend:100000"}, "stop_reading": slowing},
        {"F": 1.2578125, "stable": False, "converged": False, "verdict": "FAILED",
         "iterations": 90000, "exit_reason": "not_slowing", "stop_reading": sliding},
    ]
    summary = fem.ssrm_run_summary(_run(trials, 1.25, 1.2578125),
                                   fem_data={"unit_system": "si"})
    print("    " + summary)
    check("the standing edge is quoted by its slowing movement",
          "was still moving at iteration 100,000, but slowing" in summary
          and "had fallen 59% over the last 50,000" in summary)
    check("from the estimated resting state, with the hold test",
          "The solver checked whether it comes to rest and found that it does, "
          "0.0006 m further on, and that it stays there" in summary
          and "counted as standing at 0.023 m" in summary)
    s_as_is = fem.creep_sentence(dict(slowing, seed='as_is'), 1.25, "m")
    check("from where it was, when the block-end state itself certified",
          "found that it does, and that it stays there, so it was counted as "
          "standing at 0.023 m."
          in s_as_is, s_as_is)
    check("the failing edge is counted as sliding, with its ratio",
          "At F = 1.2578 it did not: over the last 50,000 iterations the "
          "movement did not slow (each block of 10,000 iterations moved the "
          "slope 96% as far as the one before), so the trial was counted as "
          "sliding." in summary)
    low = summary.lower()
    check("in the Run dialog's words",
          not any(w in low for w in ("budget", "sweep", "verdict", " edge")))
    refused = dict(rule='slowing_refused', window=50000, block=10000,
                   iteration=100000, ratio=0.8,
                   increments=[0.0185, 0.0140, 0.0111, 0.0091, 0.0075])
    trials = [
        {"F": 1.2421875, "stable": True, "converged": True,
         "verdict": "CONVERGED", "iterations": 1083, "exit_reason": "converged",
         "corrector": {"checkpoint": "vp1000"}},
        {"F": 1.25, "stable": False, "converged": False, "verdict": "AMBIGUOUS",
         "iterations": 100000, "exit_reason": "iteration_cap", "u_ratio": 1.25,
         "growth": 0.022, "stop_reading": refused},
    ]
    summary = fem.ssrm_run_summary(_run(trials, 1.2421875, 1.25))
    print("    " + summary)
    check("slowing at the limit with the corrector refused: its own sentence",
          "At F = 1.2500 the slope was still moving at the limit, but slowing "
          "(the movement per 10,000 iterations fell by 59% over the last "
          "50,000); the corrector could not find the balanced state from there, "
          "so the trial was counted as failed. The factor of safety depends on "
          "the iteration limit here. Raise Max iterations per trial and it may "
          "change." in summary)
    check("and in the Run dialog's words",
          not any(w in summary.lower()
                  for w in ("budget", "sweep", "verdict", " edge")))
    s = fem.creep_sentence(dict(sliding, ratio=1.0004), 1.375)
    check("a ratio that rounds to 1.00 reads as not slowing",
          "did not slow (each block of 10,000 iterations moved the slope 100% "
          "as far as the one before)" in s, s)
    acc_trials = [dict(t, acceleration={"on": True, "switched_off_at": None})
                  for t in trials]
    s_acc = fem.ssrm_run_summary(_run(acc_trials, 1.2421875, 1.25))
    check("an accelerated run says so, before the run time",
          "Convergence acceleration was on. The run took" in s_acc, s_acc[-90:])
    check("a plain run does not", "acceleration" not in summary)
    growing = dict(sliding, ratio=1.2)
    s = fem.creep_sentence(growing, 1.3)
    check("a growing movement says it grew", "the movement grew (each block of "
          "10,000 iterations moved the slope 120% as far as the one before)"
          in s, s)
    slip = dict(rule='not_slowing', window=50000, ratio=0.5, slip_frac=0.05,
                slip_ratio=0.99)
    s = fem.creep_sentence(slip, 1.3)
    check("a slip that did not slow is quoted by the slip",
          "the joint slip grew 5% and its rate did not slow (99% of the rate "
          "before)" in s, s)
    note = fem._verdict_note({"converged": False, "exit_reason": "not_slowing",
                              "verdict": "FAILED", "u_ratio": 1.74,
                              "stop_reading": sliding})
    check("the SSRM log line names it", "the movement did not slow (each block "
          "of 10,000 iterations moved the slope 96% as far as the one before)"
          in note, note)


def _grid_trial(budget, seeds):
    """The FEM-3 geogrid wall's trial at F = 1.25 at its page settings (the
    tutorial's lock tag), on the reference kernel, with Max iterations per trial
    set to ``budget`` and the trend reading's seeds limited to ``seeds``."""
    import contextlib
    import io
    import run_tests as rt
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    page = os.path.join(here, "docs", "tutorials", "fem03_block_wall_joints.md")
    tag = next(dict(t) for t in rt.parse_test_tags(page)
               if t.get("benchmark") == "FEM-3-grid-ssrm")
    tag["file"] = os.path.normpath(os.path.join(os.path.dirname(page), tag["file"]))
    with contextlib.redirect_stdout(io.StringIO()):
        fem_data, kw, f_min, f_max, tol = rt.build_fem_ssrm_case(tag)
    kw = dict(kw, capture_failure_state=False, max_iterations=int(budget),
              max_iterations_ceiling=int(budget))
    saved = fem.CREEP_SEEDS
    fem.CREEP_SEEDS = tuple(seeds)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            with rt._force_fast_kernel(fem, False):
                res = fem.solve_ssrm(fem_data, F_min=f_min, F_max=f_max,
                                     tolerance=tol, trial_factors=[1.25], **kw)
    finally:
        fem.CREEP_SEEDS = saved
    return next(t for t in res["trials"] if abs(t["F"] - 1.25) < 1e-12)


def check_seeds():
    print("\n6. the corrector's seeds on the geogrid wall at F = 1.25 (minutes)")
    ue = 0.017436613824979727          # the trial's elastic max|u|, m
    t = _grid_trial(100000, ('as_is', 'extrapolated'))
    rd = t.get("stop_reading") or {}
    seeds = [(a["seed"], a["certified"]) for a in rd.get("attempts") or []]
    check("at the page's 100,000 the trial is slowing and both seeds are refused",
          not t["stable"] and rd.get("rule") == "slowing_refused"
          and seeds == [("as_is", False), ("extrapolated", False)],
          f"{t['exit_reason']}, {rd.get('rule')}, {seeds}")
    for seeds_in, want in ((('as_is',), 'as_is'),
                           (('extrapolated',), 'extrapolated')):
        t = _grid_trial(200000, seeds_in)
        rd = t.get("stop_reading") or {}
        hold = (rd.get("corrector") or {}).get("hold") or {}
        u = float(rd.get("max_displacement") or 0.0)
        check(f"at a 200,000 allowance the {want} seed certifies at 140,000",
              t["stable"] and rd.get("rule") == "slowing"
              and rd.get("seed") == want and rd.get("iteration") == 140000,
              f"{t['exit_reason']}, {rd.get('rule')}, seed {rd.get('seed')}, "
              f"at {rd.get('iteration')}")
        check(f"  and the hold test holds it, at about 1.30 elastic "
              f"displacements", hold.get("held") is True
              and 1.29 <= u / ue <= 1.32, f"held {hold.get('held')}, "
              f"{u / ue:.3f} u_el")


def main():
    quick = "--quick" in sys.argv
    check_classes()
    check_ratio_not_count()
    check_extrapolation()
    check_recorded()
    check_sentences()
    if not quick:
        check_seeds()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"{len(FAILURES)} trend reading check(s) FAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("All trend reading checks passed.")


if __name__ == "__main__":
    main()
