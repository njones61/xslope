"""Checks for the TREND READING — how a trial still moving at its iteration limit
is read (`xslope.fem.creep_trend`, the 'not_slowing' and
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

  3. THE RECORDED TRACES (test/fixtures/creep_traces.json). The FEM-3 geogrid
     wall's trial at F = 1.25, at its page settings, is dying away at 100,000
     iterations; the FEM-3 wall alone at F = 1.140625 is holding steady at
     90,000 and is counted as sliding.

  4. THE SENTENCES. The closing summary quotes a standing edge counted standing
     by a slowing movement, a failing edge counted as sliding by a movement that did not slow,
     and one still slowing at the limit that the corrector could not finish, in
     the Run dialog's words.

  5. THE CORRECTOR ATTEMPT (about a minute; skipped with --quick). The Griffiths
     & Lane Example 2 slope (SSRM-G2) at F = 1.34375, built the way its figure
     producer builds it: the corrector refuses the state at every checkpoint and
     every block end before 11,200 iterations and certifies the block-end state
     as it stands at 11,200. The reading's record carries that one attempt and
     no extrapolated state, and the trial matches the one the row's stored
     result carries (docs/fem/files/xslope_griffiths2_fem_meta.json).

Run directly:  PYTHONPATH=. python3 test/creep_trend_check.py [--quick]
Exits non-zero on any failure.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402
from xslope.fem import creep_trend                         # noqa: E402

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


def check_recorded():
    print("\n3. the recorded traces")
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
    print("\n4. the sentences")
    slowing = dict(rule='slowing', window=50000, block=10000, iteration=100000,
                   increments=[0.02, 0.016, 0.0128, 0.0102, 0.0082], ratio=0.8,
                   max_displacement=0.023,
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
    check("the standing edge confirmed after slowing, in one sentence",
          "At F = 1.2500 the slope was still creeping at iteration 100,000 but "
          "slowing; the solver confirmed that it comes to rest, at 0.023 m, and "
          "counted it as standing." in summary, summary)
    standing = summary.split(" At F = 1.2578")[0]
    check("with no percentage, no distance further on, no 'stays there'",
          not any(w in standing for w in ("%", "further on", "stays there",
                                          "fallen")), standing)
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
    check("slowing at the limit with its rest unconfirmed: its own sentence",
          "At F = 1.2500 the slope was still creeping at the limit but slowing; "
          "the solver could not confirm that it comes to rest, so the trial was "
          "counted as failed. The factor of safety depends on the iteration "
          "limit here. Raise Max iterations per trial and it may change."
          in summary, summary)
    top = summary[summary.index("At F = 1.2500"):]
    check("and in the Run dialog's words, naming no mechanism",
          not any(w in top.lower()
                  for w in ("budget", "sweep", "verdict", " edge", "corrector",
                            "seed", "hold test", "out-of-balance", "%")), top)
    s = fem.creep_sentence(dict(sliding, ratio=1.0004), 1.375)
    check("a ratio that rounds to 1.00 reads as not slowing",
          "did not slow (each block of 10,000 iterations moved the slope 100% "
          "as far as the one before)" in s, s)
    acc_trials = [dict(t, acceleration={"on": True, "switched_off_at": None})
                  for t in trials]
    s_acc = fem.ssrm_run_summary(_run(acc_trials, 1.2421875, 1.25))
    check("the closing summary does not name the acceleration setting",
          "acceleration" not in s_acc and "acceleration" not in summary, s_acc[-90:])
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


#: The trial section 5 reads: the SSRM-G2 row's standing edge.
_G2_STEM = "xslope_griffiths2"
_G2_F = 1.34375
_G2_META = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "docs", "fem", "files", "xslope_griffiths2_fem_meta.json")


def _g2_trial():
    """The SSRM-G2 row's trial at F = 1.34375, on the model its figure producer
    (benchmarks/make_griffiths_figures.py) builds from the row's tag, solved on
    the reference kernel at the tag's settings."""
    import contextlib
    import io
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(here, "benchmarks"))
    import make_griffiths_figures as G
    tag = next(t for t in G.figure_tags() if G.stem(t) == _G2_STEM)
    with contextlib.redirect_stdout(io.StringIO()):
        sd = G.load_slope_data(tag["file"])
        lines, _n_reinf, _n_pile = G.extract_constraint_line_geometry(sd)
        polys = G.get_material_polygons(sd, reinf_lines=lines)
        mesh = G.build_mesh_from_polygons(
            polys, target_size=tag.get("target_size"),
            element_type=tag["element_type"], lines=lines,
            element_size_1d=sd.get("element_size_1d"),
            point_constraints=G.extract_point_constraints(sd),
            size_regions=G.extract_size_regions(sd), **G.RT._refine_kwargs(tag))
        fem_data = G.build_fem_data(sd, mesh)
        with G.RT._force_fast_kernel(fem, False):
            res = fem.solve_ssrm(fem_data, F_min=tag.get("f_min", 0.5),
                                 F_max=tag.get("f_max", 3.0),
                                 tolerance=tag.get("tolerance", 0.05),
                                 max_iterations=int(tag["max_iter"]),
                                 capture_failure_state=False,
                                 trial_factors=[_G2_F])
    return next(t for t in res["trials"] if abs(t["F"] - _G2_F) < 1e-12)


#: Fields a record of the trend reading's corrector attempt does not carry: the
#: attempt is made from the block-end state as it stands, and nothing else.
_NOT_RECORDED = ("seed", "extrapolated", "extrapolated_u_ratio", "field_ratio",
                 "extrapolation")


def check_corrector_attempt():
    print("\n5. the corrector attempt on SSRM-G2 at F = 1.34375 (about a minute)")
    t = _g2_trial()
    rd = t.get("stop_reading") or {}
    atts = t.get("corrector_attempts") or []
    labels = [str(a.get("at")) for a in atts]
    cert = [str(a.get("at")) for a in atts if a.get("certified")]
    check("the trial stands on the trend reading at 11,200 iterations",
          t["stable"] and rd.get("rule") == "slowing"
          and rd.get("iteration") == 11200
          and (t.get("corrector") or {}).get("checkpoint") == "trend:11200",
          f"{t['exit_reason']}, {rd.get('rule')}, at {rd.get('iteration')}, "
          f"checkpoint {(t.get('corrector') or {}).get('checkpoint')}")
    check("  every earlier checkpoint and block end was refused, and only "
          "trend:11200 certified",
          cert == ["trend:11200"] and labels[-1] == "trend:11200"
          and all(not a.get("certified") for a in atts[:-1])
          and {"vp300", "vp1000", "vp3000"} <= set(labels), f"{labels}")
    check("  every attempt was made from the state as it stood",
          not any(":" in lab.split("trend:", 1)[1] for lab in labels
                  if lab.startswith("trend:")), f"{labels}")
    ra = rd.get("attempts") or []
    check("  the reading's record carries one certified attempt and no "
          "extrapolation field",
          len(ra) == 1 and ra[0].get("certified") is True
          and not any(k in rd for k in _NOT_RECORDED)
          and not any(k in ra[0] for k in _NOT_RECORDED),
          f"{len(ra)} attempt(s), keys {sorted(rd)}")
    with open(_G2_META) as f:
        stored = next(x for x in json.load(f)["trials"]
                      if abs(x["F"] - _G2_F) < 1e-12)
    same = (t["iterations"] == stored["iterations"]
            and t["verdict"] == stored["verdict"]
            and t["exit_reason"] == stored["exit_reason"]
            and abs(t["u_ratio"] - stored["u_ratio"]) <= 1e-9 * stored["u_ratio"]
            and (t.get("corrector") or {}).get("checkpoint")
            == (stored.get("corrector") or {}).get("checkpoint"))
    check("  the trial is the one the row's stored result carries", same,
          f"{t['iterations']} / {stored['iterations']} iterations, u_ratio "
          f"{t['u_ratio']:.6f} / {stored['u_ratio']:.6f}")


def main():
    quick = "--quick" in sys.argv
    check_classes()
    check_ratio_not_count()
    check_recorded()
    check_sentences()
    if not quick:
        check_corrector_attempt()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"{len(FAILURES)} trend reading check(s) FAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("All trend reading checks passed.")


if __name__ == "__main__":
    main()
