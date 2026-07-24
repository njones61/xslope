"""Phase-4 acceptance smoke for the transient-seepage STUDIO layer
(plan_transient_seep.md §7).  Headless (offscreen Qt) so it needs no display.

Phase 4 wires the committed transient solver into XSlope Studio: a Transient
choice on the Run Seepage dialog with editable rapid-drawdown stage times, a
background run that produces a frame sequence, a Seep · Transient results tab with
a play bar (transport buttons, frame slider, jump-to-time, playback, stage tags),
per-frame rendering through the UNCHANGED plot_seep_solution (so the ``t = {value}
{unit}`` title annotation appears for free), sidecar save/restore of the frame
bundle, and the §6 analysis-time frame feeding an LEM/FEM ``u = seep`` run.

Checks:
  A. RunSeepDialog — the Transient choice appears only with a tseep sheet; the
     stage fields init from the sheet; options() returns the right mode; steady is
     byte-unchanged without tseep; stage_1 >= stage_2 is rejected.
  B. SeepRunner transient path — builds seep + tseep data, runs the solver, and
     emits a bundle of plottable per-frame solution dicts.
  C. TransientSeepView + play bar — set_frames, slider scrub, jump-to-time,
     play/pause timer step, per-frame render (a real flow net), the time-annotated
     title (bare when no unit declared, unit shown when declared), and the tag signal.
  D. MainWindow integration — the bundle lands in a 'Seep · Transient' tab; the
     sidecar (tseep.csv/meta) writes on the succeeded handler and restores via
     _restore_transient_sidecar; a play-bar tag writes stage times back onto the
     tseep controls (and dirties the doc); the analysis-time frame lands in
     slope_data['seep_u'] (single frame) / seep_u + seep_u2 (rapid).

Skips cleanly (exit 0) when PySide6 is not installed, like the editor round-trip
guard — engine-only installs have no studio layer to exercise.
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
    from PySide6.QtWidgets import QApplication, QMessageBox
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
from studio.dialogs import RunSeepDialog
from studio.runners import SeepRunner
from studio.transient import TransientSeepView

DAM = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam1.xlsx")

# Shared, built once (mesh + one transient solve) and reused across the view and
# MainWindow checks — a full transient solve per check would triple the runtime.
_SHARED = {}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _shared():
    """earth_dam1 + Ss/Sy + a constant-BC tseep dict (relaxes from the steady IC —
    a valid transient with real flow-net frames at every save time), meshed once,
    solved once. Returns (slope_data, mesh, bundle)."""
    if _SHARED:
        return _SHARED["d"], _SHARED["mesh"], _SHARED["bundle"]
    d = load_slope_data(DAM)
    for m in d["materials"]:
        m["Ss"] = 1e-3
        m["Sy"] = 0.2
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 16.0, "tri3")
    d["mesh"] = mesh
    d["tseep"] = {"times": [], "series": {}, "duration": 3.0, "save_interval": 1.0,
                  "save_times": [], "stage_1": 0.0, "stage_2": 2.0}
    runner = SeepRunner(d, {"mode": "transient", "stage_1": 0.0, "stage_2": 2.0})
    got = {}
    runner.succeeded.connect(lambda b: got.update(b))
    err = {}
    runner.failed.connect(lambda msg: err.setdefault("msg", msg))
    _quiet(runner._run_transient, d, mesh)
    if err:
        raise RuntimeError(f"transient runner failed: {err['msg']}")
    _SHARED.update(d=d, mesh=mesh, bundle=got)
    return d, mesh, got


# ------------------------------------------------------------------ A. dialog
def test_dialog():
    fails = []
    dlg0 = RunSeepDialog(has_bc2=False, has_tseep=False)
    o0 = dlg0.options()
    if o0.get("mode") != "steady" or set(o0) != {"mode", "bc", "tol"}:
        fails.append(f"steady (no tseep) options wrong: {o0}")
    if dlg0.run_type.findData("transient") != -1:
        fails.append("Transient choice offered when the file has no tseep sheet")

    dlg = RunSeepDialog(has_bc2=False, has_tseep=True,
                        tseep={"stage_1": 0.0, "stage_2": 7.0})
    if dlg.run_type.findData("transient") < 0:
        fails.append("Transient choice missing when tseep present")
    if dlg.run_type.currentData() != "transient":
        fails.append("Run type did not default to transient")
    o = dlg.options()
    if o.get("mode") != "transient" or o.get("bc") != 1:
        fails.append(f"transient options wrong: {o}")
    if o.get("stage_1") != 0.0 or o.get("stage_2") != 7.0:
        fails.append(f"stage fields did not init from the sheet: {o}")
    if not dlg.stage_enable.isChecked():
        fails.append("stage_enable should be on when the sheet carries stages")
    if dlg.bc.isEnabled():
        fails.append("BC selector should be disabled in transient mode")

    dlg.stage_enable.setChecked(False)
    o2 = dlg.options()
    if o2.get("stage_1") is not None or o2.get("stage_2") is not None:
        fails.append(f"stages not cleared when disabled: {o2}")
    dlg.run_type.setCurrentIndex(dlg.run_type.findData("steady"))
    if dlg.options().get("mode") != "steady":
        fails.append("switching Run type to Steady did not change the mode")
    if not dlg.bc.isEnabled():
        fails.append("BC selector should re-enable in steady mode")

    dlgv = RunSeepDialog(has_bc2=False, has_tseep=True,
                         tseep={"stage_1": 5.0, "stage_2": 2.0})
    dlgv.stage_enable.setChecked(True)
    dlgv.accept()
    if dlgv.result() != 0:
        fails.append("accept() allowed stage_1 >= stage_2")

    dlgb = RunSeepDialog(has_bc2=False, has_tseep=True, tseep={})
    if dlgb.stage_enable.isChecked():
        fails.append("stage_enable on for a sheet with no stage times")
    if dlgb.options().get("stage_1") is not None:
        fails.append("blank-stage sheet produced a stage time")
    return fails


# ------------------------------------------------------------------ B. runner
def test_runner():
    fails = []
    _d, _mesh, got = _shared()
    if got.get("mode") != "transient":
        fails.append("runner bundle missing mode=transient")
    frames = got.get("frames") or []
    if len(frames) < 2:
        fails.append(f"expected multiple frames, got {len(frames)}")
    for fr in frames:
        if any(fr.get(k) is None for k in ("head", "u", "velocity", "gradient", "time")):
            fails.append("a frame is missing a plottable key")
            break
    sol = got.get("transient") or {}
    if "mass_balance" not in sol or "frames" not in sol:
        fails.append("runner bundle 'transient' missing solver ledger/frames")
    return fails


# ------------------------------------------------------------------ C. view
def _force_draw(canvas):
    fig = canvas.figure
    fig.clear()
    canvas._draw_fn(fig)
    return fig


def _main_ax(fig):
    return max(fig.axes, key=lambda a: a.get_position().width * a.get_position().height)


def test_view():
    fails = []
    _d, _mesh, bundle = _shared()
    seep_data, frames = bundle["seep_data"], bundle["frames"]
    view = TransientSeepView()
    tagged = {}
    view.tag_requested.connect(lambda n, t: tagged.__setitem__(n, t))
    opts = {"variable": "head", "flowlines": True, "vectors": True, "phreatic": True}
    view.set_frames(seep_data, frames, opts_getter=lambda: opts,
                    style_getter=lambda: None, keep_index=False)
    n = view.frame_count()
    if n != len(frames):
        fails.append("frame_count != number of frames")
    if view._slider.maximum() != n - 1:
        fails.append("slider max not aligned to frame count")
    if view._idx != n - 1:
        fails.append("set_frames(keep_index=False) did not land on the last frame")

    view._slider.setValue(0)
    if view._idx != 0:
        fails.append("slider scrub to 0 did not update the index")

    tmid = view._times[len(view._times) // 2]
    view._readout.setText(f"{tmid:g}")
    view._on_readout()
    if abs(view.current_time - tmid) > 1e-9:
        fails.append("jump-to-time did not select the nearest frame")

    view.set_index(0)
    view._btn_play.setChecked(True)
    before = view._idx
    view._advance_playing()
    if view._idx != before + 1:
        fails.append("playback did not advance the frame")
    view._btn_play.setChecked(False)
    if view._timer.isActive():
        fails.append("pausing did not stop the timer")

    view.set_index(n - 1)
    fig = _force_draw(view.canvas)
    title = _main_ax(fig).get_title()
    if "t =" not in title:
        fails.append(f"frame title missing the time annotation: {title!r}")
    with tempfile.TemporaryDirectory() as td:
        png = os.path.join(td, "frame.png")
        fig.savefig(png, dpi=100)
        img = matplotlib.pyplot.imread(png)
        nonwhite = float(np.mean(np.any(img[..., :3] < 0.92, axis=-1)))
    if nonwhite < 0.05:
        fails.append(f"rendered frame nearly blank (nonwhite {nonwhite:.3f})")

    # declared time unit shows in the annotation; undeclared stays a bare number
    seep_data["unit_system"] = "SI"
    seep_data["time_unit"] = "day"
    view.render_current()
    if "day" not in _main_ax(_force_draw(view.canvas)).get_title():
        fails.append("declared time unit not shown in the frame title")
    seep_data.pop("unit_system", None)
    seep_data.pop("time_unit", None)

    view.set_index(0)
    view._emit_tag(1)
    if tagged.get(1) != view.current_time:
        fails.append("Set Stage 1 did not emit the current time")
    return fails


# ------------------------------------------------------- D. MainWindow integration
def test_mainwindow():
    fails = []
    from studio.main_window import MainWindow
    d, mesh, bundle = _shared()
    mw = MainWindow()
    try:
        mw.open_path(DAM)
        mw.doc.slope_data["mesh"] = mesh              # same object the bundle was built on
        mw.doc.slope_data["materials"] = d["materials"]
        mw.doc.slope_data["tseep"] = dict(d["tseep"])

        tmpdir = tempfile.mkdtemp(prefix="tseep_studio_")
        stem = os.path.join(tmpdir, "damT")
        xlsx = stem + ".xlsx"
        shutil.copy(DAM, xlsx)
        mw.doc.path = xlsx

        _quiet(mw._on_transient_seep_succeeded, bundle)
        if mw.transient_seep_view is None:
            fails.append("Seep · Transient view not created by the succeeded handler")
        labels = [mw.view_tabs.tabText(i) for i in range(mw.view_tabs.count())]
        if "Seep · Transient" not in labels:
            fails.append(f"'Seep · Transient' tab missing; tabs={labels}")
        if mw._display_panels.get(mw.transient_seep_view) is None:
            fails.append("no Display panel registered for the transient view")
        if not os.path.exists(f"{stem}_tseep.csv"):
            fails.append("succeeded handler did not write the tseep.csv sidecar")
        if not os.path.exists(f"{stem}_tseep_meta.json"):
            fails.append("succeeded handler did not write the meta json")

        # analysis-time frame (§6): single frame -> seep_u; rapid -> seep_u/seep_u2
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

        # play-bar tag writes the stage time back onto the tseep controls + dirties
        mw.doc._dirty = False
        t_tag = mw.transient_seep_view._times[0]
        mw._on_transient_tag(1, t_tag)
        if mw.doc.slope_data["tseep"].get("stage_1") != t_tag:
            fails.append("tag did not update tseep stage_1")
        if not mw.doc.dirty:
            fails.append("tagging a stage did not mark the document dirty")

        # RESTORE round-trip: drop the view, rebuild from the on-disk sidecar
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
        else:
            rb = mw.doc.results.get("transient_seep")
            if not rb or len(rb["frames"]) != len(bundle["frames"]):
                fails.append("restored frame count mismatch")
            if not _force_draw(mw.transient_seep_view.canvas).axes:
                fails.append("restored frame failed to render")

        # _sync_sidecars removes the files when the transient result is gone
        mw.doc.results.pop("transient_seep", None)
        mw.transient_seep_view = None
        _quiet(mw._sync_sidecars, stem)
        if os.path.exists(f"{stem}_tseep.csv"):
            fails.append("_sync_sidecars did not remove the stale tseep.csv")

        shutil.rmtree(tmpdir, ignore_errors=True)
    finally:
        mw.doc._dirty = False        # skip the close-time "save changes?" modal
        mw.close()
    return fails


def main():
    print("transient seepage phase-4 studio smoke:")
    checks = [("run dialog (transient choice + stages)", test_dialog),
              ("seep runner transient path", test_runner),
              ("transient view + play bar", test_view),
              ("mainwindow integration + sidecar", test_mainwindow)]
    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:42s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll phase-4 studio smoke checks passed.")


if __name__ == "__main__":
    main()
