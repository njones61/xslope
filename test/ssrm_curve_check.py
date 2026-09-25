"""Checks for the strength reduction run's displacement curve and closing summary.

What this file locks:

  1. THE TRIAL RECORD. On the viscoplastic path every trial solve_ssrm records
     carries ``max_displacement`` — the solve's own maximum displacement, the
     quantity a saved field's meta records under the same name — and the run ends
     with a closing summary on ``result['summary']``. One small solve of Griffiths
     & Lane Example 1 on a coarse mesh (seconds).

  2. THE ROUND TRIP. export_fem_solution writes the trials, displacements
     included, into ``{stem}_fem_meta.json``; import_fem_meta gives them back
     unchanged, and that restored record draws the curve. A meta file written
     before the record existed imports with no trials and no error, and the
     curve says why it cannot be drawn.

  3. THE PLOT. plot_ssrm_curve on a synthetic record holding a converged trial,
     one stopped by the sweep budget, one diverging and one past the displacement
     limit: every trial marked, each ending named in the key, the factor of
     safety ruled. plot_fem_results draws it alone and as the fourth of four
     stacked panels, and refuses a record it cannot draw with the reason.

  4. THE SUMMARY. ssrm_run_summary says how each bracket edge was decided, in
     each of the forms the failing edge can take, and says in so many words when
     the failing edge was the sweep budget's verdict.

Run directly:  PYTHONPATH=. python3 test/ssrm_curve_check.py
Exits non-zero on any failure.
"""
import contextlib
import io
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                           # noqa: E402
matplotlib.use("Agg")

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402

FAILURES = []
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


def _griffiths_coarse():
    """Griffiths & Lane Example 1 on the coarse tri6 mesh hybrid_criterion_check
    uses — the smallest shipped SSRM model."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    xlsx = os.path.join(_HERE, "docs", "fem", "files", "xslope_griffiths1.xlsx")
    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(xlsx)
        mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                        target_size=8.0, element_type="tri6")
        return fem.build_fem_data(slope_data, mesh)


_RUN = {}


def _solved():
    """One viscoplastic strength reduction run, shared by sections 1 and 2."""
    if "run" not in _RUN:
        fem_data = _griffiths_coarse()
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            result = fem.solve_ssrm(fem_data, F_min=1.0, F_max=2.0,
                                    tolerance=0.1, max_iterations=4000,
                                    fem_solver="viscoplastic",
                                    capture_failure_state=False)
        _RUN["run"] = (fem_data, result, out.getvalue())
    return _RUN["run"]


# ===================== 1. the trial record =====================

def check_trial_record():
    print("\n1. the viscoplastic trial record carries max_displacement")
    fem_data, result, log = _solved()
    trials = result.get("trials") or []
    check("the run bisected to a factor of safety",
          result.get("converged") and result.get("FS") is not None,
          f"FS={result.get('FS')}")
    check("it recorded trials", len(trials) >= 3, f"{len(trials)} trials")
    values = [t.get("max_displacement") for t in trials]
    check("every trial carries a finite max_displacement",
          all(v is not None and np.isfinite(v) and v > 0 for v in values),
          f"{values}")
    stood = [t["max_displacement"] for t in trials if t.get("stable")]
    fell = [t["max_displacement"] for t in trials if not t.get("stable")]
    check("a trial that failed moved further than every trial that stood",
          bool(stood) and bool(fell) and min(fell) > max(stood),
          f"stood {max(stood or [0]):.3g}, fell {min(fell or [0]):.3g}")
    # The last converged trial is the field the run hands back, so its recorded
    # displacement is that field's own.
    last = result["last_solution"]
    same = [t for t in trials if t.get("stable")
            and abs(t["F"] - float(last["F"])) < 1e-12]
    check("the recorded value is the solve's own max_displacement",
          bool(same) and abs(same[-1]["max_displacement"]
                             - float(last["max_displacement"])) < 1e-12)

    summary = result.get("summary") or ""
    check("the run carries a closing summary",
          summary.startswith(f"The factor of safety is {result['FS']:.3f}"),
          summary[:80])
    check("the summary is printed, once, as the run's last words",
          log.rstrip().endswith(summary) and log.count(summary) == 1)
    check("the summary names both edges and the wall time",
          "standing edge" in summary and "failing edge" in summary
          and "The run took" in summary, summary)


# ===================== 2. the round trip =====================

def check_round_trip():
    print("\n2. export then import round-trips the trials")
    from xslope.fem import (export_fem_solution, import_fem_meta,
                            import_fem_solution, ssrm_run_record)
    from xslope.plot_fem import ssrm_curve_unavailable
    fem_data, result, _log = _solved()
    tmp = tempfile.mkdtemp(prefix="xslope_ssrm_curve_")
    stem = os.path.join(tmp, "model")
    meta = dict(ssrm_run_record(result, fem_data, {"tolerance": 0.1}),
                FS=result["FS"], analysis="ssrm")
    with contextlib.redirect_stdout(io.StringIO()):
        export_fem_solution(fem_data, result["last_solution"], stem, meta=meta)
    check("no trial CSV is written beside the meta",
          not any(name.endswith("_fem_trials.csv") for name in os.listdir(tmp)),
          f"{sorted(os.listdir(tmp))}")
    read = import_fem_meta(stem) or {}
    back = read.get("trials") or []
    check("the meta sidecar gives the trials back",
          len(back) == len(result["trials"]), f"{len(back)} trials")
    check("with every displacement unchanged",
          [t.get("max_displacement") for t in back]
          == [float(t["max_displacement"]) for t in result["trials"]])
    check("the restored record draws the curve",
          ssrm_curve_unavailable(read) is None, ssrm_curve_unavailable(read))
    solution = import_fem_solution(fem_data, stem)
    check("the reloaded field does not carry the trials; the meta does",
          "trials" not in solution)

    # A file written before the record existed: no trials, no error.
    with open(f"{stem}_fem_meta.json") as f:
        old = json.load(f)
    old.pop("trials", None)
    with open(f"{stem}_fem_meta.json", "w") as f:
        json.dump(old, f)
    try:
        import_fem_solution(fem_data, stem)
        read = import_fem_meta(stem) or {}
        ok = "trials" not in read
    except Exception as exc:                                  # noqa: BLE001
        ok, read = False, {"error": repr(exc)}
    check("a pre-record file imports with no trials and no error", ok, str(read)[:80])
    why = ssrm_curve_unavailable(read)
    check("and the curve gives a one-line reason instead of drawing",
          bool(why) and "\n" not in why, why)
    # And one whose trials predate the displacement field.
    stale = {"trials": [{k: v for k, v in t.items() if k != "max_displacement"}
                        for t in result["trials"]]}
    why = ssrm_curve_unavailable(stale)
    check("trials saved before they carried a displacement are refused too",
          bool(why) and "saved before" in why, why)


# ===================== 3. the plot =====================

_RECORD = {
    "FS": 1.375, "final_interval": [1.35, 1.40],
    "trials": [
        {"F": 1.0, "stable": True, "converged": True, "exit_reason": "converged",
         "max_displacement": 0.05},
        {"F": 2.0, "stable": False, "converged": False, "exit_reason": "diverging",
         "max_displacement": 1.9},
        {"F": 1.5, "stable": False, "converged": False,
         "exit_reason": "disp_limit", "max_displacement": 1.2},
        {"F": 1.35, "stable": True, "converged": True, "exit_reason": "converged",
         "max_displacement": 0.1},
        {"F": 1.40, "stable": False, "converged": False,
         "exit_reason": "iteration_cap", "growth": 0.4, "u_ratio": 6.0,
         "iterations": 100000, "max_displacement": 0.6},
    ]}


def check_plot():
    print("\n3. the plot, on a synthetic record")
    import matplotlib.figure as mplfig
    from xslope.plot_fem import (SSRM_CURVE_TITLE, plot_fem_results,
                                 plot_ssrm_curve)

    fig = mplfig.Figure(figsize=(7.0, 4.0))
    ax = fig.add_subplot(111)
    plot_ssrm_curve(ax, _RECORD, fem_data={"unit_system": "SI"})
    marked = sum(len(line.get_xdata()) for line in ax.lines
                 if line.get_linestyle() == "None")
    check("every trial is marked", marked == len(_RECORD["trials"]),
          f"{marked} of {len(_RECORD['trials'])}")
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    for wanted in ("stood (converged)",
                   "did not converge within the sweep budget", "diverging",
                   "past the displacement limit", "final bracket", "FS = 1.375"):
        check(f"the key names {wanted!r}", wanted in labels, f"{labels}")
    title = ax.get_legend().get_title().get_text()
    check("the key says where an open marker sits", "when the trial stopped" in title,
          title)
    check("the axes are titled", ax.get_title() == SSRM_CURVE_TITLE)
    check("the displacement axis carries the declared length unit",
          ax.get_ylabel().endswith("(m)"), ax.get_ylabel())
    check("the open markers are open",
          all(line.get_markerfacecolor() == "white" for line in ax.lines
              if line.get_linestyle() == "None"
              and line.get_marker() in ("s", "v", "^")))
    fig.canvas.draw()
    box = ax.get_legend().get_window_extent()
    pts = ax.transData.transform([(t["F"], t["max_displacement"])
                                  for t in _RECORD["trials"]])
    check("the key covers no trial",
          not any(box.contains(x, y) for x, y in pts))

    fem_data = _griffiths_coarse()
    _f, result, _log = _solved()
    solution = result["last_solution"]
    fig = mplfig.Figure(figsize=(8.0, 5.0))
    with contextlib.redirect_stdout(io.StringIO()):
        plot_fem_results(fem_data, solution, plot_type=["ssrm_curve"], fig=fig,
                         ssrm_record=_RECORD, fs=1.375)
    check("plot_fem_results draws the curve alone on one axes",
          len(fig.axes) == 1 and fig.axes[0].get_aspect() != 1.0)
    fig = mplfig.Figure(figsize=(8.0, 12.0))
    four = ["deformation", "shear_strain", "displace_vector", "ssrm_curve"]
    with contextlib.redirect_stdout(io.StringIO()):
        _fig, axes = plot_fem_results(fem_data, solution, plot_type=four, fig=fig,
                                      ssrm_record=_RECORD, fs=1.375)
    check("and as the fourth of four stacked panels",
          len(axes) == 4 and axes[3].get_title() == SSRM_CURVE_TITLE)
    lefts = [round(a.get_position().x0, 3) for a in axes]
    check("the fourth panel keeps the stack's left edge",
          max(lefts) - min(lefts) < 0.05, f"{lefts}")
    try:
        plot_fem_results(fem_data, solution, plot_type=["ssrm_curve"],
                         fig=mplfig.Figure())
        refused = None
    except ValueError as exc:
        refused = str(exc)
    check("a request with no record is refused with the reason",
          bool(refused) and "\n" not in refused, refused)


# ===================== 4. the summary =====================

def _summary(failing, **extra):
    trials = [{"F": 1.35, "stable": True, "converged": True,
               "verdict": "CONVERGED", "iterations": 812,
               "exit_reason": "converged"},
              dict({"F": 1.40, "stable": False, "converged": False,
                    "iterations": 100000}, **failing)]
    result = dict({"converged": True, "FS": 1.375,
                   "final_interval": (1.35, 1.40), "trials": trials,
                   "failure_criterion": "hybrid", "elapsed_time": 664.0}, **extra)
    return fem.ssrm_run_summary(result, {"unit_system": "SI"})


def check_summary():
    print("\n4. the closing summary")
    # A trial still moving, but slowly, when its budget ran out (AMBIGUOUS).
    s = _summary({"exit_reason": "iteration_cap", "verdict": "AMBIGUOUS",
                  "growth": 0.5, "u_ratio": 5.0, "max_displacement": 0.25})
    check("it opens on the answer and its bracket",
          s.startswith("The factor of safety is 1.375, the midpoint of the "
                       "bracket from F = 1.3500 to F = 1.4000."), s)
    check("the standing edge is the trial that converged",
          "converged in 812 sweeps" in s, s)
    check("a slow budget stop is named in those words",
          "stopped at the 100,000-sweep budget with the section still moving, "
          "but slowly" in s, s)
    check("and the number is said to be the budget's",
          "the budget's, not the slope's, and a longer budget may move it" in s, s)
    check("and why it was counted as failed",
          "too slow to call it a failure, so it was counted as failed" in s, s)
    check("with the displacement, its elastic multiple and its growth",
          "0.25 m" in s and "5.0 times the elastic response" in s
          and "0.025 m" in s and "last 25,000 sweeps" in s, s)

    # A trial running away when its budget ran out (FAILED).
    s = _summary({"exit_reason": "iteration_cap", "verdict": "FAILED",
                  "growth": 0.74, "u_ratio": 4.88, "max_displacement": 0.277},
                 failure_criterion="non_convergence")
    check("a runaway at the budget is decided by the displacement evidence",
          "decided by the displacement evidence" in s
          and "ran out of sweeps at the 100,000-sweep budget while running away"
          in s, s)
    check("and is called a failure in progress",
          "That is a failure in progress, not a budget effect." in s, s)
    check("and the number is NOT called the budget's", "the budget's" not in s, s)
    check("with the displacement, its elastic multiple and its growth",
          "0.277 m, 4.9 times the elastic response, and still growing by 0.042 m"
          in s, s)
    check("and the wall time", s.endswith("The run took 11 min 4 s."), s)

    s = _summary({"exit_reason": "diverging", "verdict": "FAILED",
                  "iterations": 2501})
    check("a diverging edge says so, with its sweep",
          "diverged at sweep 2,501" in s and "budget" not in s, s)
    s = _summary({"exit_reason": "disp_limit", "verdict": "FAILED",
                  "iterations": 640})
    check("a displacement-limit edge says so, with its sweep",
          "passed it at sweep 640" in s and "displacement limit" in s, s)
    s = _summary({"exit_reason": "iteration_cap", "verdict": "STABLE_STUCK",
                  "growth": 0.0, "u_ratio": 1.02, "max_displacement": 0.05},
                 failure_criterion="non_convergence")
    check("a budget stop with the section no longer moving is not called moving",
          "still moving" not in s and "no longer moving" in s, s)
    s = fem.ssrm_run_summary({"converged": False, "FS": None, "trials": [],
                              "error": "SSRM: the slope still stands at F = 10.00, "
                                       "the highest factor the search may try.",
                              "elapsed_time": 42.0})
    check("a run with no answer says so", s.startswith(
        "No factor of safety was found. The slope still stands at F = 10.00"), s)
    s = fem.ssrm_run_summary({"probe": True, "trials": [
        {"F": 1.3, "stable": True, "iterations": 300},
        {"F": 1.4, "stable": False, "iterations": 4000}]})
    check("a probe names each trial and no factor of safety",
          "F = 1.3000 stood" in s and "F = 1.4000 did not stand" in s, s)


def run():
    """Every check; returns the list of failures (empty = pass), the entry point
    run_tests.py's module checks call."""
    del FAILURES[:]
    for fn in (check_trial_record, check_round_trip, check_plot, check_summary):
        try:
            fn()
        except Exception as exc:                              # noqa: BLE001
            import traceback
            traceback.print_exc()
            check(f"{fn.__name__} ran", False, repr(exc))
    return list(FAILURES)


def main():
    print("SSRM displacement curve and closing summary:")
    run()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURE(S):")
        for name in FAILURES:
            print(f"  - {name}")
        raise SystemExit(1)
    print("\nAll SSRM curve checks passed.")


if __name__ == "__main__":
    main()
