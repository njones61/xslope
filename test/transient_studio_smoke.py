"""Acceptance smoke for the transient-seepage STUDIO layer. Headless (offscreen Qt)
so it needs no display.

Covers the transient Studio surface after the inputs-editor redesign:

  A. RunSeepDialog (simplified) — a Transient run-type choice appears only with a
     tseep sheet; the dialog carries NO rapid-drawdown stage widgets (stage times are
     model inputs, edited under Inputs → Transient) but shows a caption saying so;
     options() returns {mode, bc, tol} only; steady is byte-unchanged without tseep.
  B. SeepRunner transient path — builds seep + tseep data, runs the solver, and emits
     a bundle of plottable per-frame solution dicts; the run reports DETERMINATE
     progress (simulated-time fraction, monotonic 0→1, reaching 100%) and is
     CANCELLABLE (a cancel request stops the march and emits `cancelled` with no
     stored result).
  C. TransientSeepView + play bar — the bar has NO stage-tag buttons; playback ADVANCES
     frames on a real timer (distinct fields render at distinct times); a transient frame
     carries no stream function, so it renders WITHOUT flow lines (a storage-release state
     has no flow net — the plotter's guard draws none and shows the honest "no through-flow"
     note even when a caller requests them) and WITHOUT crashing.
  D. TransientEditor — opens on a tseep-bearing file and round-trips its data (no enable
     checkbox; the save-times list sits beside the series table); editing a stage time and
     a series value is reflected in result_tseep() and the live plot (series curve + stage
     reference lines); an all-blank editor on a no-tseep file → None (steady).
  E. MainWindow integration — a 'Transient' row is in the inputs tree; the bundle lands
     in a 'Seep · Transient' tab; the sidecar writes/restores; the analysis-time frame
     feeds seep_u / seep_u+seep_u2.

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

import numpy as np
import matplotlib
matplotlib.use("Agg")

try:
    from PySide6.QtWidgets import QApplication, QMainWindow, QMessageBox
    from PySide6.QtTest import QTest
except Exception:                       # engine-only install — nothing to test
    print("transient studio smoke: PySide6 not installed — skipped.")
    sys.exit(0)

_app = QApplication.instance() or QApplication([])

# Modal dialogs must never block a headless run.
QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import (build_seep_data, build_tseep_data, run_transient_seepage,
                         _transient_frame_solution)
from studio.dialogs import RunSeepDialog
from studio.runners import SeepRunner
from studio.transient import TransientSeepView
from studio.editors import TransientDialog, CATEGORY_EDITORS
from studio.display_panels import SeepDisplayPanel

# A real reservoir-drawdown transient fixture: a "pool" head series that draws down,
# so early frames have through-flow (real flow nets) and the drawdown exercises the
# save schedule + stage times.
DAM = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam_tseep.xlsx")
NO_TSEEP = os.path.join(_REPO, "docs/inputs/slope/xslope_dam.xlsx")

_SHARED = {}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _shared():
    """earth_dam_tseep meshed once and solved once (transient). Returns
    (slope_data, mesh, seep_data, tseep_data, solution, frames)."""
    if _SHARED:
        return (_SHARED["d"], _SHARED["mesh"], _SHARED["seep"], _SHARED["ts"],
                _SHARED["sol"], _SHARED["frames"])
    d = load_slope_data(DAM)
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 16.0, "tri3")
    d["mesh"] = mesh
    seep_data = _quiet(build_seep_data, mesh, d, seep_bc=1)
    tseep_data = _quiet(build_tseep_data, d)
    sol = _quiet(run_transient_seepage, seep_data, tseep_data, verbose=False)
    unconf = bool(sol.get("unconfined"))
    frames = [_transient_frame_solution(seep_data, fr["head"], fr["u"], fr.get("phi"),
                                        fr.get("inflow"), fr.get("outflow"), unconf,
                                        time=fr["time"])
              for fr in sol["frames"]]
    _SHARED.update(d=d, mesh=mesh, seep=seep_data, ts=tseep_data, sol=sol, frames=frames)
    return d, mesh, seep_data, tseep_data, sol, frames


# ------------------------------------------------------------------ A. dialog
def test_dialog():
    fails = []
    dlg0 = RunSeepDialog(has_bc2=False, has_tseep=False)
    o0 = dlg0.options()
    if o0.get("mode") != "steady" or set(o0) != {"mode", "bc", "tol", "max_iter"}:
        fails.append(f"steady (no tseep) options wrong: {o0}")
    if dlg0.run_type.findData("transient") != -1:
        fails.append("Transient choice offered when the file has no tseep sheet")

    dlg = RunSeepDialog(has_bc2=False, has_tseep=True)
    if dlg.run_type.findData("transient") < 0:
        fails.append("Transient choice missing when tseep present")
    if dlg.run_type.currentData() != "transient":
        fails.append("Run type did not default to transient")
    o = dlg.options()
    if o.get("mode") != "transient" or o.get("bc") != 1:
        fails.append(f"transient options wrong: {o}")
    if set(o) != {"mode", "bc", "tol"}:
        fails.append(f"transient options carry stage keys (should not): {o}")
    if hasattr(dlg, "seep_bc"):
        fails.append("run dialog exposes a BC selector; a steady run solves every set")
    # No stage widgets anywhere on the dialog.
    for attr in ("stage_1", "stage_2", "stage_enable", "transient_group"):
        if hasattr(dlg, attr):
            fails.append(f"run dialog still exposes stage widget '{attr}'")
    if not hasattr(dlg, "transient_caption"):
        fails.append("run dialog missing the 'edit under Inputs → Transient' caption")

    dlg.run_type.setCurrentIndex(dlg.run_type.findData("steady"))
    if dlg.options().get("mode") != "steady":
        fails.append("switching Run type to Steady did not change the mode")
    if dlg.options().get("bc") not in (1, "both"):
        fails.append(f"steady options carry an unexpected bc: {dlg.options()}")
    return fails


# ------------------------------------------------------------------ B. runner + progress + cancel
def test_runner():
    fails = []
    d, mesh, seep_data, tseep_data, sol, frames = _shared()

    # runner bundle shape
    runner = SeepRunner(d, {"mode": "transient"})
    got = {}
    runner.succeeded.connect(lambda b: got.update(b))
    err = {}
    runner.failed.connect(lambda m: err.setdefault("msg", m))
    _quiet(runner._run_transient, d, mesh)
    if err:
        fails.append(f"transient runner failed: {err['msg']}")
    if got.get("mode") != "transient":
        fails.append("runner bundle missing mode=transient")
    if len(got.get("frames") or []) < 2:
        fails.append("expected multiple frames")
    for fr in (got.get("frames") or []):
        if any(fr.get(k) is None for k in ("head", "u", "velocity", "gradient", "time")):
            fails.append("a frame is missing a plottable key")
            break

    # progress: monotonic 0→duration, reaches (duration, duration)
    dur = tseep_data["duration"]
    calls = []
    _quiet(run_transient_seepage, seep_data, tseep_data, verbose=False,
           progress_callback=lambda t, D: (calls.append((t, D)) or True))
    ts = [t for t, _ in calls]
    if not calls:
        fails.append("progress_callback was never called")
    else:
        if any(b < a - 1e-9 for a, b in zip(ts, ts[1:])):
            fails.append("progress times not monotonic")
        if abs(ts[-1] - dur) > 1e-9:
            fails.append(f"progress did not reach the duration ({ts[-1]} vs {dur})")
        if any(abs(D - dur) > 1e-9 for _, D in calls):
            fails.append("progress duration argument inconsistent")

    # cancel: a falsy return stops the march early with a cancelled result
    n = {"c": 0}

    def _cancel_cb(t, D):
        n["c"] += 1
        return n["c"] < 3
    solc = _quiet(run_transient_seepage, seep_data, tseep_data, verbose=False,
                  progress_callback=_cancel_cb)
    if not solc.get("cancelled"):
        fails.append("cancelled solve did not set cancelled=True")
    if n["c"] != 3:
        fails.append(f"cancel did not stop promptly (callback fired {n['c']}x)")
    if len(solc["frames"]) >= len(sol["frames"]):
        fails.append("cancelled run saved as many frames as a full run")

    # runner-level cancel: pre-set the flag → the worker emits `cancelled`, not succeeded
    r2 = SeepRunner(d, {"mode": "transient"})
    r2.cancel()
    outcome = {}
    r2.succeeded.connect(lambda b: outcome.setdefault("done", b))
    r2.cancelled.connect(lambda: outcome.setdefault("cancelled", True))
    _quiet(r2._run_transient, d, mesh)
    if not outcome.get("cancelled"):
        fails.append("pre-cancelled runner did not emit cancelled")
    if "done" in outcome:
        fails.append("pre-cancelled runner emitted a succeeded bundle")
    return fails


# ------------------------------------------------------------------ C. view + play + flow lines
def _force_draw(canvas):
    fig = canvas.figure
    fig.clear()
    canvas._draw_fn(fig)
    return fig


def _main_ax(fig):
    return max(fig.axes, key=lambda a: a.get_position().width * a.get_position().height)


def _gids(fig):
    return set(c.get_gid() for ax in fig.axes for c in ax.collections)


def _pmhash(canvas):
    import hashlib
    pm = canvas._pixitem.pixmap() if canvas._pixitem else None
    if pm is None:
        return None
    return hashlib.md5(bytes(pm.toImage().constBits())).hexdigest()[:12]


def test_view():
    fails = []
    d, mesh, seep_data, tseep_data, sol, frames = _shared()

    view = TransientSeepView()
    # No stage-tag buttons / signal on the play bar.
    for attr in ("_btn_tag1", "_btn_tag2", "tag_requested", "_emit_tag"):
        if hasattr(view, attr):
            fails.append(f"play bar still exposes removed stage-tag member '{attr}'")

    panel = SeepDisplayPanel(d.get("materials"))
    win = QMainWindow()
    win.setCentralWidget(view)
    win.resize(760, 520)
    win.show()
    QTest.qWait(80)
    view.set_frames(seep_data, frames, opts_getter=panel.options,
                    style_getter=lambda: None, keep_index=False)
    QTest.qWait(150)

    # Transient frames carry no stream function — a storage-release state has no flow
    # net — so no frame renders flow lines. Even when a caller requests them (the steady
    # display panel here does), the plotter's guard draws none and shows the honest
    # "no through-flow" note rather than crashing. The frame still renders and its title
    # carries the frame time.
    mid = len(frames) // 2
    view.set_index(mid)
    try:
        fig = _force_draw(view.canvas)
    except Exception as e:
        fig = None
        fails.append(f"transient frame crashed the renderer: {e!r}")
    if fig is not None:
        if "FLOWLINES" in _gids(fig):
            fails.append("transient frame drew flow lines (phi is not computed)")
        if "t =" not in _main_ax(fig).get_title():
            fails.append(f"frame title missing time: {_main_ax(fig).get_title()!r}")
        subs = " ".join(t.get_text() for t in _main_ax(fig).texts)
        if "no through-flow" not in subs:
            fails.append(f"transient frame missing the honest subtitle: {subs!r}")

    # Playback ADVANCES frames on the real timer and renders distinct fields.
    view.set_frames(seep_data, frames, opts_getter=panel.options,
                    style_getter=lambda: None, keep_index=False)
    view.set_index(0)
    QTest.qWait(120)
    h0 = _pmhash(view.canvas)
    view._btn_play.setChecked(True)     # start playback
    if not view._timer.isActive():
        fails.append("play did not start the playback timer")
    seen = set()
    idxs = []
    for _ in range(6):
        QTest.qWait(520)
        seen.add(_pmhash(view.canvas))
        idxs.append(view._idx)
    view._btn_play.setChecked(False)
    if max(idxs) < 2:
        fails.append(f"playback did not advance frames (idxs={idxs})")
    if len(seen) < 3:
        fails.append(f"playback rendered too few distinct frames ({len(seen)})")
    win.close()
    return fails


# ------------------------------------------------------------------ D. transient editor
def test_editor():
    fails = []
    d = load_slope_data(DAM)
    orig = d.get("tseep")
    if not orig:
        return ["earth_dam_tseep fixture carries no tseep data"]

    dlg = TransientDialog(orig, d, None)
    # The 'Enable transient analysis' checkbox was removed — blank fields alone mean
    # "no transient data", so a tseep-bearing file just round-trips its data.
    if hasattr(dlg, "_enable"):
        fails.append("editor still exposes the removed 'Enable transient analysis' checkbox")
    # New layout: the save-times list sits beside the series table (both present).
    if dlg._save_table is None or dlg._series_table is None:
        fails.append("editor missing the series or save-times table")
    rt = dlg.result_tseep()
    if rt is None:
        fails.append("editor produced None for a tseep-bearing file")

    def eq(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            return a.keys() == b.keys() and all(eq(a[k], b[k]) for k in a)
        if isinstance(a, list) and isinstance(b, list):
            return len(a) == len(b) and all(eq(x, y) for x, y in zip(a, b))
        if isinstance(a, float) or isinstance(b, float):
            if a is None or b is None:
                return a is b
            return abs(a - b) < 1e-9
        return a == b
    if not eq(orig, rt):
        fails.append(f"editor did not round-trip tseep: {orig} != {rt}")

    # Edit a stage time and a series value; both are reflected in result_tseep().
    dlg._stage_2.setText("120")
    dlg._series_table.item(0, 1).setText("17.5")
    out = dlg.result_tseep()
    if out.get("stage_2") != 120.0:
        fails.append(f"stage_2 edit not reflected: {out.get('stage_2')}")
    first_series = next(iter(out["series"].values()))
    if first_series[0] != 17.5:
        fails.append(f"series-value edit not reflected: {first_series[0]}")

    # Live plot draws the series curve and the stage reference lines.
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    dlg._draw_plot(ax)
    series_lines = [ln for ln in ax.lines if ln.get_marker() == "o"]
    stage_lines = [ln for ln in ax.lines if ln.get_linestyle() == "--"]
    if not series_lines:
        fails.append("editor plot drew no series curve")
    if len(stage_lines) != 2:
        fails.append(f"editor plot drew {len(stage_lines)} stage lines (expected 2)")
    plt.close(fig)

    # A no-tseep file: the all-blank editor produces None (round-trips a steady file),
    # with no checkbox to toggle — blank IS the steady state.
    d2 = load_slope_data(NO_TSEEP)
    dlg2 = TransientDialog(d2.get("tseep"), d2, None)
    if hasattr(dlg2, "_enable"):
        fails.append("no-tseep editor still exposes the removed enable checkbox")
    if dlg2.result_tseep() is not None:
        fails.append("blank editor produced a non-None tseep (should be None => steady)")
    return fails


# ------------------------------------------------------- E. MainWindow integration
def test_mainwindow():
    fails = []
    from studio.main_window import MainWindow
    d, mesh, seep_data, tseep_data, sol, frames = _shared()
    bundle = {"mode": "transient", "seep_data": seep_data, "transient": sol,
              "frames": frames, "options": {"mode": "transient"}}
    mw = MainWindow()
    try:
        mw.open_path(DAM)
        mw.doc.slope_data["mesh"] = mesh

        # The inputs tree carries a clickable 'Transient' row wired to the editor.
        labels = []
        root = mw.inputs_tree.invisibleRootItem()
        for i in range(root.childCount()):
            labels.append(root.child(i).text(0))
        if "Transient" not in labels:
            fails.append(f"'Transient' row missing from the inputs tree: {labels}")
        if "transient" not in CATEGORY_EDITORS:
            fails.append("no 'transient' category editor registered")

        tmpdir = tempfile.mkdtemp(prefix="tseep_studio_")
        stem = os.path.join(tmpdir, "damT")
        xlsx = stem + ".xlsx"
        shutil.copy(DAM, xlsx)
        mw.doc.path = xlsx

        _quiet(mw._on_transient_seep_succeeded, bundle)
        if mw.transient_seep_view is None:
            fails.append("Seep · Transient view not created by the succeeded handler")
        tabs = [mw.view_tabs.tabText(i) for i in range(mw.view_tabs.count())]
        if "Seep · Transient" not in tabs:
            fails.append(f"'Seep · Transient' tab missing; tabs={tabs}")
        if not os.path.exists(f"{stem}_tseep.csv"):
            fails.append("succeeded handler did not write the tseep.csv sidecar")

        # analysis-time frame (§6): single -> seep_u; rapid -> seep_u/seep_u2
        mw.transient_seep_view.set_index(mw.transient_seep_view.frame_count() - 1)
        mw.doc.slope_data["seep_u"] = None
        mw._apply_transient_analysis_frame(rapid=False)
        if mw.doc.slope_data.get("seep_u") is None:
            fails.append("analysis-time hook did not populate seep_u")
        mw.doc.slope_data["seep_u"] = None
        mw.doc.slope_data["seep_u2"] = None
        mw._apply_transient_analysis_frame(rapid=True)
        if (mw.doc.slope_data.get("seep_u") is None
                or mw.doc.slope_data.get("seep_u2") is None):
            fails.append("rapid analysis-time hook did not stage seep_u/seep_u2")

        # RESTORE round-trip from the on-disk sidecar
        view0 = mw.transient_seep_view
        idx = mw.view_tabs.indexOf(view0)
        if idx >= 0:
            mw.view_tabs.removeTab(idx)
        mw._display_panels.pop(view0, None)
        mw.transient_seep_view = None
        mw.doc.results.pop("transient_seep", None)
        _quiet(mw._restore_transient_sidecar, mesh, stem)
        if mw.transient_seep_view is None:
            fails.append("restore did not rebuild the transient view from the sidecar")

        shutil.rmtree(tmpdir, ignore_errors=True)
    finally:
        mw.doc._dirty = False
        mw.close()
    return fails


# ------------------------------------------------------------------ F. seep display panel + water levels
def test_display_panel():
    """A transient seep display panel OMITS the flow-net-only controls (Flow lines /
    Base material) and defaults the Water levels overlay ON; a steady panel keeps the
    controls and defaults Water levels OFF."""
    fails = []
    d, *_ = _shared()
    mats = d.get("materials")
    steady = SeepDisplayPanel(mats, transient=False)
    trans = SeepDisplayPanel(mats, transient=True)

    # Steady view keeps the flow-net controls (parented into its form layout).
    if steady.flowlines.parentWidget() is None:
        fails.append("steady seep panel is missing the Flow lines control")
    if steady.base_mat.parentWidget() is None:
        fails.append("steady seep panel is missing the Base material control")
    # Transient view OMITS them entirely (never added to the layout, so unparented).
    if trans.flowlines.parentWidget() is not None:
        fails.append("transient seep panel still exposes the Flow lines control")
    if trans.base_mat.parentWidget() is not None:
        fails.append("transient seep panel still exposes the Base material control")
    if trans.options().get("flowlines"):
        fails.append("transient seep panel requests flow lines (should never)")
    # Water-levels default: OFF steady (byte-stable), ON transient (playback aid).
    if steady.options().get("show_bc_levels"):
        fails.append("steady seep panel defaults Water levels ON (should be OFF)")
    if not trans.options().get("show_bc_levels"):
        fails.append("transient seep panel defaults Water levels OFF (should be ON)")
    return fails


def test_water_levels():
    """The show_bc_levels overlay draws each boundary's level, evaluated at the frame's
    time: the reservoir pool sits at h(0) on the first frame and has visibly dropped by
    a late frame; with the overlay off, no level line is drawn."""
    import matplotlib.pyplot as plt
    from xslope.plot_seep import plot_seep_solution
    from xslope.plot import _eval_bc_series_at
    fails = []
    d, mesh, seep_data, tseep_data, sol, frames = _shared()

    def _level_ys(frame, on):
        fig = plt.figure(figsize=(7, 3))
        _quiet(plot_seep_solution, seep_data, frame, fig=fig, levels=8,
               flowlines=False, mesh=False, show_legend=False, show_title=False,
               show_bc_levels=on)
        ax = _main_ax(fig)
        ys = [float(np.max(ln.get_ydata())) for ln in ax.lines
              if ln.get_gid() == "BC_LEVEL"]
        plt.close(fig)
        return ys

    first, late = frames[0], frames[-1]
    y_first = _level_ys(first, True)
    y_late = _level_ys(late, True)
    y_off = _level_ys(late, False)

    if not y_first:
        fails.append("show_bc_levels ON drew no water level on the first frame")
    if y_off:
        fails.append("show_bc_levels OFF still drew a water level line")
    # The reservoir is the highest level; it should equal the pool series at the frame
    # time and drop from the full-pool level to the drawn-down level.
    tseep = seep_data.get("tseep")
    h0 = _eval_bc_series_at(tseep, "pool", float(first["time"]))
    hL = _eval_bc_series_at(tseep, "pool", float(late["time"]))
    if y_first and abs(max(y_first) - h0) > 0.25:
        fails.append(f"first-frame reservoir level {max(y_first)} != h(0)={h0}")
    if y_late and abs(max(y_late) - hL) > 0.25:
        fails.append(f"late-frame reservoir level {max(y_late)} != h(t)={hL}")
    if y_first and y_late and not (max(y_late) < max(y_first) - 1e-6):
        fails.append(f"pool level did not drop: first {max(y_first)}, late {max(y_late)}")
    return fails


def main():
    print("transient seepage studio smoke:")
    checks = [("run dialog (simplified, no stage widgets)", test_dialog),
              ("seep runner + progress + cancel", test_runner),
              ("transient view + playback + flow lines", test_view),
              ("transient inputs editor + live plot", test_editor),
              ("mainwindow integration + sidecar", test_mainwindow),
              ("seep display panel (transient omits flow-net controls)", test_display_panel),
              ("transient water-level overlay", test_water_levels)]
    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll transient studio smoke checks passed.")


if __name__ == "__main__":
    main()
