"""A report never re-solves an engine that is already solved.

The Analysis Report documents what was run. Two of the three engines can take
minutes to run — a seepage march, a strength reduction bisection — and both save
their answers beside the model, which Studio restores on open and
``solutions_from_sidecars`` assembles from Python. A report that quietly solved
either of them again would cost the user those minutes for numbers they already
had, and could print an answer that is not the one their results view is showing.

The one engine a report may run is limit equilibrium, and only for a method the
report was asked to document that has no result yet
(:func:`xslope.report.run_requested_methods`). That run is a solve of its own,
announced in the progress bar, and it happens exactly once per method.

So this drives a real report over a model whose seepage AND finite element
solutions are already attached, with every seepage and FEM solve entry point in
the package replaced by a stub that raises:

  * ``xslope.seep.run_seepage_analysis``, ``solve_confined``,
    ``solve_unsaturated`` and ``run_transient_seepage`` — the steady, confined,
    unconfined and transient flow solves;
  * ``xslope.seep.solve_flow_function_confined`` and
    ``solve_flow_function_unsaturated`` — the stream function behind the flow
    net, which a report that recomputed it rather than reading the saved field
    would call while drawing the figure;
  * ``xslope.fem.solve_fem`` and ``solve_ssrm`` — the deformation solve and the
    strength reduction search.

Every name is confirmed to exist on its module BEFORE it is replaced, so a
renamed solver makes this check fail rather than silently guard nothing.

The report is built with every content option on, so the guard covers the
figures as well as the tables: a figure that re-derives a field it was handed is
the way this would come back. What must come out is a document on disk, a
Seepage Analysis section, a Deformation and Strength Reduction section, and
exactly one limit equilibrium run — counted by wrapping
``xslope.search.run_lem_analysis``, the one function every LEM run goes through.

Model: ``docs/seep/files/xslope_johnson_res.xlsx``, which ships with both a
seepage sidecar and a finite element sidecar beside it. No search: with no LEM
bundle in hand the report runs the surface the input specifies
(:func:`xslope.report.analysis_options`), which is one solve of about a second.

Run directly:  PYTHONPATH=. python3 test/report_no_resolve_check.py
"""

import os
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: The model reported on: a reservoir slope with a steady seepage solution and a
#: finite element solution saved beside it.
MODEL = os.path.join(_REPO, "docs", "seep", "files", "xslope_johnson_res.xlsx")

#: The one method the report is asked to document. Nothing has solved it, so the
#: report runs it — which is the run this check counts.
METHOD = "bishop"

#: Every way the package solves a flow problem or a finite element problem, as
#: ``(module path, attribute)``. A report may call none of them.
SOLVE_ENTRY_POINTS = [
    ("xslope.seep", "run_seepage_analysis"),
    ("xslope.seep", "solve_confined"),
    ("xslope.seep", "solve_unsaturated"),
    ("xslope.seep", "run_transient_seepage"),
    ("xslope.seep", "solve_flow_function_confined"),
    ("xslope.seep", "solve_flow_function_unsaturated"),
    ("xslope.fem", "solve_fem"),
    ("xslope.fem", "solve_ssrm"),
]


class Resolved(AssertionError):
    """Raised from a stub: the report called a solver it had an answer for."""


def _forbid_solves(failures):
    """Replace every solve entry point with a stub that raises, and return the
    callable that puts them all back.

    A name that is not there is reported rather than skipped: the guard is only
    worth anything while it names the functions the package actually solves
    with.
    """
    import importlib

    saved = []
    for module_name, attr in SOLVE_ENTRY_POINTS:
        module = importlib.import_module(module_name)
        if not callable(getattr(module, attr, None)):
            failures.append(f"{module_name}.{attr} is not a function any more; "
                            f"this check is guarding a name that has moved")
            continue

        def stub(*_args, _where=f"{module_name}.{attr}", **_kw):
            raise Resolved(f"the report called {_where}; the solution it was "
                           f"given was already solved")

        saved.append((module, attr, getattr(module, attr)))
        setattr(module, attr, stub)

    def restore():
        for module, attr, original in saved:
            setattr(module, attr, original)

    return restore


def check_report_reports_what_is_solved(failures):
    """A full report over an already-solved model, with the solvers booby-trapped."""
    from xslope.fileio import load_slope_data
    from xslope.report import (DEFAULT_OPTIONS, generate_report,
                               solutions_from_sidecars)

    if not os.path.exists(MODEL):
        failures.append(f"missing model {MODEL}")
        return

    slope_data = load_slope_data(MODEL)
    notes = []
    solutions = solutions_from_sidecars(MODEL, slope_data, notes)
    for note in notes:
        failures.append(f"a companion file could not be read: {note}")
    if not solutions.get("seep"):
        failures.append("the model's seepage sidecar produced no solution, so "
                        "this check would prove nothing")
        return
    if not solutions.get("fem"):
        failures.append("the model's finite element sidecar produced no "
                        "solution, so this check would prove nothing")
        return
    if solutions.get("lem"):
        failures.append("the sidecars carried a limit equilibrium bundle; the "
                        "report would then run no method and the count below "
                        "would prove nothing")
        return

    import xslope.search as search

    real_run = search.run_lem_analysis
    ran = []

    def counting_run(slope_data_, method, *args, **kwargs):
        ran.append(method)
        return real_run(slope_data_, method, *args, **kwargs)

    search.run_lem_analysis = counting_run
    restore = _forbid_solves(failures)
    options = dict(DEFAULT_OPTIONS, method=METHOD, input_path=MODEL,
                   title="Report Without Re-solving")
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "no_resolve_report.docx")
            try:
                ok, out = generate_report(slope_data, solutions, options, path)
            except Resolved as exc:
                failures.append(str(exc))
                return
            if not ok:
                failures.append(f"the report could not be generated: {out}")
                return
            if not os.path.exists(path) or os.path.getsize(path) < 10000:
                failures.append("no document was written")
                return
            titles = [title for _level, title in out["report"].section_titles()]
    finally:
        restore()
        search.run_lem_analysis = real_run

    if "Seepage Analysis" not in titles:
        failures.append(f"the seepage solution attached to the model got no "
                        f"section: {titles}")
    if not any(t.startswith("Deformation") for t in titles):
        failures.append(f"the finite element solution attached to the model got "
                        f"no section: {titles}")
    if ran != [METHOD]:
        failures.append(f"the report made {len(ran)} limit equilibrium run(s) "
                        f"({ran}); it was asked to document one method that had "
                        f"no result")


def check_no_solver_is_reachable(failures):
    """The report module names no seepage or FEM solver at all.

    The run above proves nothing was called on this model; this proves there is
    no branch to call one on another — a "solve it if it is missing" path that
    this model, being solved, would never reach.
    """
    import re

    source = open(os.path.join(_REPO, "xslope", "report.py"),
                  encoding="utf-8").read()
    # Prose about the solvers is what the module is full of; a call is the name
    # followed by an opening parenthesis, or the name imported.
    code = "\n".join(line for line in source.splitlines()
                     if not line.lstrip().startswith("#"))
    for _module, attr in SOLVE_ENTRY_POINTS:
        if re.search(rf"\b{attr}\s*\(", code) or re.search(
                rf"import[^\n]*\b{attr}\b", code):
            failures.append(f"xslope/report.py calls {attr}; a report documents "
                            f"the engines that were run and runs none of them")


def check_the_trap_is_armed(failures):
    """The stubs really do replace the solvers, and are really put back.

    A guard made of monkeypatches that missed its targets would pass this file
    on any code at all. So one of them is called on purpose and has to raise,
    and every original has to be back on its module afterwards — this check runs
    inside the same process as the rest of the suite.
    """
    import importlib

    before = {}
    for module_name, attr in SOLVE_ENTRY_POINTS:
        before[(module_name, attr)] = getattr(
            importlib.import_module(module_name), attr, None)

    restore = _forbid_solves(failures)
    try:
        module_name, attr = SOLVE_ENTRY_POINTS[0]
        try:
            getattr(importlib.import_module(module_name), attr)()
        except Resolved:
            pass
        else:
            failures.append(f"{module_name}.{attr} was replaced by a stub that "
                            f"does not raise; the guard would pass a report "
                            f"that re-solved everything")
    finally:
        restore()

    for (module_name, attr), original in before.items():
        if getattr(importlib.import_module(module_name), attr, None) is not original:
            failures.append(f"{module_name}.{attr} was left stubbed out, which "
                            f"would break every check that runs after this one")


CHECKS = [check_report_reports_what_is_solved, check_no_solver_is_reachable,
          check_the_trap_is_armed]


def run():
    """Returns a list of failure descriptions (empty on success)."""
    failures = []
    for chk in CHECKS:
        try:
            chk(failures)
        except Exception as e:                                 # noqa: BLE001
            failures.append(f"{chk.__name__} raised {type(e).__name__}: {e}")
    return failures


if __name__ == "__main__":
    bad = run()
    for f in bad:
        print(f"FAIL  {f}")
    print("report_no_resolve_check: "
          + (f"{len(bad)} failure(s)" if bad else f"{len(CHECKS)} checks passed"))
    sys.exit(1 if bad else 0)
