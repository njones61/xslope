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
     one stopped by the iteration limit, one diverging and one past the displacement
     limit: every trial marked, each ending named in the key, the factor of
     safety ruled. plot_fem_results draws it alone and as the fourth of four
     stacked panels, and refuses a record it cannot draw with the reason.

  4. THE SUMMARY. ssrm_run_summary says what happened at each end of the
     bracket, in each of the ways the trial at the top can end, in the Run
     dialog's words; it says when the answer depends on the iteration limit, and
     no summary uses "budget", "edge", "sweep" or "verdict".

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
    div = [t for t in trials if t.get("exit_reason") == "diverging"]
    check("a trial ended by a rule records the reading that fired it",
          bool(div) and all((t.get("stop_reading") or {}).get("rule")
                            == "diverging" and t["stop_reading"].get("u_ratio")
                            for t in div),
          f"{[t.get('stop_reading') for t in div]}")
    check("a converged trial records no reading",
          all(t.get("stop_reading") is None for t in trials if t.get("converged")))
    check("the summary says what happened at each end, and the wall time",
          "reached equilibrium" in summary and "did not:" in summary
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
    joined = [line for line in ax.lines if line.get_linestyle() == "-"]
    stood = sorted(t["F"] for t in _RECORD["trials"] if t["stable"])
    check("the line joins only the trials that reached equilibrium",
          len(joined) == 1 and list(joined[0].get_xdata()) == stood,
          f"{[list(l.get_xdata()) for l in joined]}")
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    for wanted in ("reached equilibrium",
                   "stopped at the iteration limit, still moving",
                   "displacements ran away",
                   "past the displacement limit", "final bracket", "FS = 1.375"):
        check(f"the key names {wanted!r}", wanted in labels, f"{labels}")
    title = ax.get_legend().get_title().get_text()
    check("the key says where an open marker sits",
          title == "open marker: where the trial was when it was stopped, "
                   "still moving", title)
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

_SUMMARIES = []


def _said(s):
    """Keep every summary this file produces, for the jargon check at the end."""
    _SUMMARIES.append(s)
    return s


def _summary(failing, **extra):
    trials = [{"F": 1.35, "stable": True, "converged": True,
               "verdict": "CONVERGED", "iterations": 812,
               "exit_reason": "converged"},
              dict({"F": 1.40, "stable": False, "converged": False,
                    "iterations": 100000}, **failing)]
    result = dict({"converged": True, "FS": 1.375,
                   "final_interval": (1.35, 1.40), "trials": trials,
                   "failure_criterion": "hybrid", "elapsed_time": 664.0}, **extra)
    return _said(fem.ssrm_run_summary(result, {"unit_system": "SI"}))


def check_summary():
    print("\n4. the closing summary")
    del _SUMMARIES[:]
    limit = ("The factor of safety depends on the iteration limit here. Raise "
             "Max iterations per trial and it may change.")

    # Still moving, but slowly, when the iteration limit stopped it (AMBIGUOUS).
    s = _summary({"exit_reason": "iteration_cap", "verdict": "AMBIGUOUS",
                  "growth": 0.5, "u_ratio": 5.0, "max_displacement": 0.25})
    check("it opens on the answer and its bracket",
          s.startswith("The factor of safety is 1.375, the midpoint of the "
                       "bracket F = 1.3500 to 1.4000."), s)
    check("the bottom of the bracket reached equilibrium",
          "At F = 1.3500 the slope reached equilibrium in 812 iterations." in s, s)
    check("a slow stop at the limit is named in those words",
          "At F = 1.4000 it did not: the trial hit the 100,000-iteration limit "
          "while still moving, but slowly" in s, s)
    check("and why it was counted as failed",
          "That is too much movement to call the slope settled and too little to "
          "call it a failure, so the trial was counted as failed." in s, s)
    check("and that the answer depends on the limit", limit in s, s)
    check("with the displacement, its elastic multiple and its growth",
          "its largest displacement was 0.25 m, 5.0 times the elastic value, and "
          "had grown by 0.025 m over the last 25,000 iterations" in s, s)
    check("and the wall time", s.endswith("The run took 11 min 4 s."), s)

    # Still moving fast when the iteration limit stopped it (FAILED).
    s = _summary({"exit_reason": "iteration_cap", "verdict": "FAILED",
                  "growth": 0.74, "u_ratio": 4.88, "max_displacement": 0.277},
                 failure_criterion="non_convergence")
    check("a fast-moving stop at the limit says the slope was failing",
          "the trial hit the 100,000-iteration limit while still moving fast — "
          "its largest displacement was 0.277 m, 4.9 times the elastic value, "
          "and had grown by 0.042 m" in s
          and "The slope was failing; more iterations would only have let it "
              "move further." in s, s)
    check("and does not say the answer depends on the limit",
          "depends on the iteration limit" not in s, s)

    s = _summary({"exit_reason": "diverging", "verdict": "FAILED",
                  "iterations": 2231})
    check("a diverging trial says the displacements ran away",
          "At F = 1.4000 it did not: the displacements ran away at iteration "
          "2,231." in s, s)
    s = _summary({"exit_reason": "disp_limit", "verdict": "FAILED",
                  "iterations": 640})
    check("a displacement-limit trial says so, with its iteration",
          "it passed the displacement limit at iteration 640." in s, s)
    s = _summary({"exit_reason": "inconclusive", "verdict": "AMBIGUOUS",
                  "growth": 0.01, "u_ratio": 1.1, "max_displacement": 0.05})
    check("an inconclusive trial is an open question",
          "the trial hit the 100,000-iteration limit with its out-of-balance "
          "force still falling. That is neither an equilibrium nor a failure, so "
          "the factor of safety carries it as an open question." in s, s)
    s = _summary({"exit_reason": "iteration_cap", "verdict": "STABLE_STUCK",
                  "growth": 0.0, "u_ratio": 1.02, "max_displacement": 0.05},
                 failure_criterion="non_convergence")
    check("a stop with the displacements stopped is not called moving",
          "moving" not in s and "with its displacements stopped" in s, s)
    stood = {"F": 1.35, "stable": True, "converged": False,
             "verdict": "STABLE_STUCK", "iterations": 100000,
             "exit_reason": "iteration_cap"}
    s = _said(fem.ssrm_run_summary({
        "converged": True, "FS": 1.375, "final_interval": (1.35, 1.40),
        "trials": [stood, {"F": 1.40, "stable": False, "iterations": 300,
                           "exit_reason": "diverging", "verdict": "FAILED"}]}))
    check("a trial that stood without meeting the tolerance is counted standing",
          "At F = 1.3500 the slope did not meet the force tolerance within "
          "100,000 iterations, but its displacements had stopped, so it was "
          "counted as standing." in s, s)
    s = _said(fem.ssrm_run_summary({
        "converged": False, "FS": None, "trials": [],
        "error": "SSRM: the slope still reaches equilibrium at F = 10.00, the top "
                 "of the range the search may try, so the factor of safety is "
                 "above it.", "elapsed_time": 42.0}))
    check("a run with no answer says so", s.startswith(
        "No factor of safety was found. The slope still reaches equilibrium at "
        "F = 10.00"), s)
    s = _said(fem.ssrm_run_summary({"probe": True, "trials": [
        {"F": 1.3, "stable": True, "iterations": 300},
        {"F": 1.4, "stable": False, "iterations": 4000}]}))
    check("a probe names each trial and no factor of safety",
          "At F = 1.3000 the slope reached equilibrium" in s
          and "At F = 1.4000 the slope did not reach equilibrium" in s, s)

    # Every rule's reading, quoted in one clause.
    s = _summary({"exit_reason": "steady_slip", "verdict": "FAILED",
                  "iterations": 90001,
                  "stop_reading": {"rule": "steady_slip", "window": 45000,
                                   "slip_frac": 0.06, "growth": 0.4,
                                   "rate_ratio": 1.02, "iteration": 90001}})
    check("sliding on the joints quotes the slip and its rate",
          "At F = 1.4000 it did not: over the last 45,000 iterations the joint "
          "slip grew 6% and its rate did not slow, so the run stopped waiting at "
          "iteration 90,001 and counted the slope as sliding." in s, s)
    s = _summary({"exit_reason": "steady_slip", "verdict": "FAILED",
                  "iterations": 90001,
                  "stop_reading": {"rule": "steady_slip", "window": 45000,
                                   "slip_frac": 0.0042, "growth": 0.4,
                                   "rate_ratio": 0.93, "iteration": 90001}})
    check("a slip under 1% and a slightly slower rate are quoted as they are",
          "grew 0.42% and its rate slowed by only 7%" in s, s)
    s = _summary({"exit_reason": "diverging", "verdict": "FAILED",
                  "iterations": 341,
                  "stop_reading": {"rule": "diverging", "signal": "runaway",
                                   "u_ratio": 15.12, "iteration": 341,
                                   "window": 2000, "gain": 3.0}})
    check("running away quotes the elastic multiple",
          "At F = 1.4000 it did not: the largest displacement reached 15.1 times "
          "the elastic value at iteration 341." in s, s)
    s = _summary({"exit_reason": "disp_limit", "verdict": "FAILED",
                  "iterations": 640,
                  "stop_reading": {"rule": "disp_limit", "displacement": 0.85,
                                   "limit": 0.5, "iteration": 640}})
    check("the displacement limit quotes the displacement and the limit",
          "At F = 1.4000 it did not: the displacement reached 0.85 m at iteration "
          "640, past the limit of 0.5 m." in s, s)
    s = _summary({"exit_reason": "inconclusive", "verdict": "AMBIGUOUS",
                  "stop_reading": {"rule": "inconclusive", "oob_from": 0.0031,
                                   "oob_to": 0.0024, "window": 1000,
                                   "force_tol": 0.001}})
    check("the iteration ceiling quotes the fall of the out-of-balance force",
          "still falling (from 0.0031 to 0.0024 over the last 1,000 iterations). "
          "That is neither an equilibrium nor a failure" in s, s)
    settled = {"F": 1.35, "stable": True, "converged": False,
               "verdict": "JOINT_SETTLED", "iterations": 60000,
               "exit_reason": "joint_settled",
               "stop_reading": {"rule": "joint_settled", "window": 15000,
                                "slip_frac": 0.00004, "growth": 0.001,
                                "soil_oob": 0.0004, "force_tol": 0.001}}
    s = _said(fem.ssrm_run_summary({
        "converged": True, "FS": 1.375, "final_interval": (1.35, 1.40),
        "trials": [settled, {"F": 1.40, "stable": False, "iterations": 300,
                             "exit_reason": "diverging", "verdict": "FAILED"}]}))
    check("settled joints quote the slip's growth",
          "At F = 1.3500 the joint slip and the displacements had stopped (the "
          "slip grew 0.004% over the last 15,000 iterations) while the joint "
          "forces kept flickering, so it was counted as standing." in s, s)
    old_file = dict(settled)
    old_file.pop("stop_reading")
    s = _said(fem.ssrm_run_summary({
        "converged": True, "FS": 1.375, "final_interval": (1.35, 1.40),
        "trials": [old_file, {"F": 1.40, "stable": False, "iterations": 300,
                              "exit_reason": "steady_slip", "verdict": "FAILED"}]}))
    check("a record saved before readings gets the earlier sentences",
          "its joints and displacements had stopped, so it was counted as "
          "standing" in s
          and "the slope was sliding steadily on its joints at iteration 300." in s,
          s)

    # The joint rule reports its reading without changing when it fires.
    traces = os.path.join(_HERE, "test", "fixtures", "joint_traces.json")
    with open(traces) as f:
        rows = json.load(f)
    same, filled = True, True
    for t in rows.values():
        args = (t["slip"], t["oob_soil"], t["disp"], t["u_elastic_scale"], 1e-3)
        kw = dict(joint_oob_hist=t["oob_joint"], budget=t["budget"],
                  sample_every=t["sample_every"])
        reading = {}
        plain = fem.joint_verdict(*args, **kw)
        read = fem.joint_verdict(*args, reading=reading, **kw)
        same &= plain == read
        filled &= (reading.get("rule") == read) if read else not reading
    check("joint_verdict fires the same with a reading asked for", same)
    check("and fills the reading exactly when it fires", filled)

    # The run's own summary from section 1 joins the synthetic ones.
    _f, result, _log = _solved()
    _SUMMARIES.append(result.get("summary") or "")
    import re
    jargon = [(w, s) for s in _SUMMARIES
              for w in ("budget", "edge", "sweep", "verdict")
              if re.search(w, s, re.IGNORECASE)]
    check("no summary says budget, edge, sweep or verdict", not jargon,
          f"{jargon[:1]}")
    contrast = [s for s in _SUMMARIES if re.search(r", not (by |the )", s)]
    check("no summary is built on an 'X, not Y' contrast", not contrast,
          f"{contrast[:1]}")


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
