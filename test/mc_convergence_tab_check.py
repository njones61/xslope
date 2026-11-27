"""Standing check on the Convergence sub-tab of Studio's Reliability · MC view.

A Monte Carlo result's tab carries the FS histogram (Plot), the per-parameter
table (Table) and, between them, Convergence: the running probability of failure
against the number of realizations, drawn by ``xslope.plot.
plot_reliability_convergence``, the same function a script calls. The check runs
a small Monte Carlo campaign on a copy of the shipped reliability sample through
Studio's own runner and success handler, then asks:

* the sub-tab is there, between Plot and Table, and showing;
* its line is the running P_f of the result's own realizations, in order, one
  point per realization, ending at the result's P_f;
* it is the library's figure: the same lines as ``plot_reliability_convergence``
  drawn on its own;
* the Display panel's title and legend toggles reach it, and the run's stopping
  target draws its band;
* the panel grays what the showing sub-tab ignores: the fitted overlays off the
  histogram, the title on the Table;
* a response-surface result in the same tab hides it (its samples are a stride
  through the surrogate draws, not a sequence that converged), and a Monte Carlo
  result shows it again;
* an input edit and opening a file each close the tab and take the sub-tab with
  it, and the next run builds it afresh.

The window is built against a scratch settings store (``tools.studio_settings``)
and opens a copy of the workbook in a temporary folder. 300 realizations on a
fixed circle: a second or two.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from tools import studio_settings                      # noqa: E402

#: The shipped reliability sample (two parameters with standard deviations).
MODEL = os.path.join(_REPO, "docs/inputs/slope/xslope_reliability.xlsx")
N = 300
OPTS = {"app_mode": "lem", "engine": "mc", "method": "bishop", "search": False,
        "n_samples": N, "rng_seed": 20240117, "distribution": "normal",
        "sampling": "lhs", "converge_rel": None, "num_slices": 40}


def _flush():
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    QApplication.instance().processEvents()


def _run_mc(sd):
    """The Monte Carlo bundle exactly as Studio's runner emits it."""
    from studio.runners import ReliabilityRunner
    got = []
    runner = ReliabilityRunner(sd, dict(OPTS))
    runner.succeeded.connect(got.append)
    runner.failed.connect(lambda m: got.append(RuntimeError(m)))
    with contextlib.redirect_stdout(io.StringIO()):
        runner.run()                    # in this thread: the signal is direct
    if not got or isinstance(got[0], Exception):
        raise RuntimeError(f"the Monte Carlo run failed: {got}")
    return got[0]


def _succeed(mw, bundle):
    """Hand a result to the window as a finished run does, without its box."""
    from PySide6.QtWidgets import QMessageBox
    real = QMessageBox.information
    QMessageBox.information = staticmethod(lambda *a, **k: None)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mw._on_rel_succeeded(bundle)
    finally:
        QMessageBox.information = real


def _tab_state(view):
    """(labels, visible flags) of a result view's sub-tabs."""
    t = view.tabs
    return ([t.tabText(i) for i in range(t.count())],
            [t.isTabVisible(i) for i in range(t.count())])


def _lines(fig):
    """{label: (x, y)} of the lines on a figure's single axes."""
    import numpy as np
    ax = fig.axes[0]
    return {ln.get_label(): (np.asarray(ln.get_xdata(), float),
                             np.asarray(ln.get_ydata(), float))
            for ln in ax.get_lines()}


def _legend_texts(fig):
    leg = fig.axes[0].get_legend() if fig.axes else None
    return [t.get_text() for t in leg.get_texts()] if leg else []


def _check_drawn(mw, bundle, fails):
    import numpy as np
    from matplotlib.figure import Figure
    from xslope.plot import plot_reliability_convergence
    view = mw.reliability_hist_canvas
    if view is None:
        return fails.append("a Monte Carlo result opened no Reliability · MC tab")
    labels, shown = _tab_state(view)
    if labels != ["Plot", "Convergence", "Table"]:
        fails.append(f"the MC view's sub-tabs are {labels}, expected "
                     f"['Plot', 'Convergence', 'Table']")
        return
    if not shown[1]:
        fails.append("the Convergence sub-tab is hidden for a Monte Carlo result")
    conv = view.plot_tab("Convergence")
    view.tabs.setCurrentWidget(conv)
    _flush()
    conv.render_now()
    if not conv.figure.axes:
        return fails.append("the Convergence sub-tab drew nothing")
    rel = bundle["reliability"]
    fs = np.asarray(rel["fs_samples"], float)
    lines = _lines(conv.figure)
    run = lines.get("running $P_f$ (empirical)")
    if run is None:
        return fails.append(f"no running P_f line on the Convergence sub-tab "
                            f"(lines: {sorted(lines)})")
    x, y = run
    ok = np.isfinite(fs)
    want = np.cumsum((fs < 1.0) & ok) / np.maximum(np.cumsum(ok), 1) * 100
    if x.size != fs.size or fs.size != N:
        fails.append(f"the running P_f has {x.size} points for {fs.size} "
                     f"realizations ({N} requested)")
    elif not (np.array_equal(x, np.arange(1, N + 1))
              and np.allclose(y, want, rtol=0, atol=1e-12)):
        fails.append("the running P_f is not the result's realizations in order")
    if abs(y[-1] - rel["pf_empirical"] * 100) > 1e-9:
        fails.append(f"the running P_f ends at {y[-1]:.4f}%, the result's P_f is "
                     f"{rel['pf_empirical'] * 100:.4f}%")
    ref = _lines(plot_reliability_convergence(rel, fig=Figure(figsize=(9, 5.5))))
    if sorted(ref) != sorted(lines) or any(
            not (np.array_equal(ref[k][0], lines[k][0])
                 and np.array_equal(ref[k][1], lines[k][1], equal_nan=True))
            for k in ref):
        fails.append("Studio's Convergence lines differ from "
                     "plot_reliability_convergence's own")
    if conv.figure.axes[0].get_title() == "" or not _legend_texts(conv.figure):
        fails.append("the Convergence plot opened without its title or legend")

    # The Display panel's toggles reach the sub-tab; the stopping target's band.
    panel = mw._display_panels.get(view)
    real_opts = panel.options
    panel.options = lambda: {**real_opts(), "show_title": False,
                             "show_legend": False}
    try:
        mw._rerender_reliability_histogram()
        conv.render_now()
        if conv.figure.axes[0].get_title() or _legend_texts(conv.figure):
            fails.append("the Display panel's title/legend toggles do not reach "
                         "the Convergence sub-tab")
    finally:
        panel.options = real_opts
    mw.doc.results["reliability"] = dict(bundle, converge_rel=0.10)
    mw._rerender_reliability_histogram()
    conv.render_now()
    if not any("stop target" in t for t in _legend_texts(conv.figure)):
        fails.append("a run with a stopping target draws no target band")
    mw.doc.results["reliability"] = bundle


@studio_settings.isolated
def run():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    QApplication.instance() or QApplication([])
    import shiboken6
    from studio.main_window import MainWindow

    tmp = tempfile.mkdtemp(prefix="mc_conv_tab_")
    model = os.path.join(tmp, os.path.basename(MODEL))
    shutil.copy(MODEL, model)
    mw = MainWindow()
    mw.resize(1400, 900)
    mw.show()
    fails = []
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mw.open_path(model)
        bundle = _run_mc(mw.doc.slope_data)
        _succeed(mw, bundle)
        before = len(fails)
        _check_drawn(mw, bundle, fails)
        print(f"  {'Monte Carlo result: sub-tab drawn':<44s} "
              f"{'ok' if len(fails) == before else 'FAIL'}")

        # The Display panel grays what the showing sub-tab ignores: the fitted
        # overlays are the histogram's alone; the Table draws no title either.
        before = len(fails)
        view = mw.reliability_hist_canvas
        panel = mw._display_panels.get(view) if view is not None else None
        if panel is not None and view.plot_tab("Convergence") is not None:
            want = {"Plot": (True, True), "Convergence": (False, True),
                    "Table": (False, False)}
            for i in range(view.tabs.count()):
                label = view.tabs.tabText(i)
                view.tabs.setCurrentIndex(i)
                got = (panel._boxes["show_fits"].isEnabled(),
                       panel._show_title.isEnabled())
                if got != want[label]:
                    fails.append(f"on the {label} sub-tab the overlays/Title "
                                 f"toggles are enabled {got}, expected "
                                 f"{want[label]}")
            view.tabs.setCurrentIndex(0)
        print(f"  {'Display toggles grayed by sub-tab':<44s} "
              f"{'ok' if len(fails) == before else 'FAIL'}")

        # A response-surface result in the same tab hides it; MC shows it again.
        before = len(fails)
        view = mw.reliability_hist_canvas
        if view is not None and view.plot_tab("Convergence") is not None:
            # The view the success handler shows (its summary box reads the
            # surrogate's own keys, which this stand-in does not carry).
            mw.doc.results["reliability"] = dict(bundle, engine="rs")
            with contextlib.redirect_stdout(io.StringIO()):
                mw._show_reliability_histogram()
            labels, shown = _tab_state(view)
            if shown[labels.index("Convergence")]:
                fails.append("a response-surface result shows the Convergence "
                             "sub-tab")
            _succeed(mw, bundle)
            labels, shown = _tab_state(view)
            if not shown[labels.index("Convergence")]:
                fails.append("a Monte Carlo result after a response-surface one "
                             "left the Convergence sub-tab hidden")
        print(f"  {'response surface hides it, MC shows it':<44s} "
              f"{'ok' if len(fails) == before else 'FAIL'}")

        # Teardown: an edit, then opening a file, each take the sub-tab along;
        # the next run builds it again.
        for name, close in (("an input edit", lambda: mw.invalidate_results()),
                            ("opening a file", lambda: mw.open_path(model))):
            before = len(fails)
            view = mw.reliability_hist_canvas
            conv = view.plot_tab("Convergence") if view is not None else None
            with contextlib.redirect_stdout(io.StringIO()):
                close()
            _flush()
            if mw.reliability_hist_canvas is not None:
                fails.append(f"{name} left the Reliability · MC tab open")
            if conv is not None and shiboken6.isValid(conv):
                fails.append(f"{name} closed the tab but its Convergence canvas "
                             f"lives on")
            _succeed(mw, _run_mc(mw.doc.slope_data))
            _check_drawn(mw, mw.doc.results["reliability"], fails)
            print(f"  {name + ', then a new run':<44s} "
                  f"{'ok' if len(fails) == before else 'FAIL'}")
    finally:
        mw.doc._dirty = False
        mw.close()
        _flush()
        shutil.rmtree(tmp, ignore_errors=True)
    return fails


def main():
    print("Reliability · MC Convergence sub-tab:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll Convergence sub-tab checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
