"""Studio-capture pipeline: regenerate the transient-seepage Studio screenshots
that ``docs/studio/analysis.md`` and ``docs/studio/editing.md`` embed, headlessly
(offscreen Qt) so it needs no display and produces byte-stable, native-quality grabs.

This is the same offscreen-grab practice the reliability / parametric dialog
images were captured with: build the real Studio widget, ``.show()`` it under the
``offscreen`` QPA platform, let the layout settle, then ``QWidget.grab()`` the
widget to a PNG. Nothing here touches a live display or the user's window.

Three captures land in ``docs/studio/images/``:

  * ``analysis_run_seep_transient.png`` — the Run Seepage dialog in **Transient**
    mode. Built from :class:`studio.dialogs.RunSeepDialog` with a tseep sheet present
    (so the Transient choice appears and defaults on). The dialog carries NO
    rapid-drawdown stage fields — those are model inputs, edited under Inputs →
    Transient — just a caption pointing there. Sized to its ``sizeHint`` to match the
    neighboring ``analysis_*_dialog.png`` images.
  * ``analysis_seep_transient_playbar.png`` — the **Seep · Transient** results view
    with the play bar (transport buttons, frame slider, ``t =`` readout, Speed
    selector — no stage-tag buttons) under a rendered flow-net frame, grabbed on a
    through-flow frame so the flow net is fully developed.
  * ``editing_transient_editor.png`` — the **Transient** inputs editor
    (:class:`studio.editors.TransientDialog`): run controls, the extra-save-times
    list, the named time-series table, and the live series-vs-time plot with the
    rapid-drawdown stage reference lines. Captured on the earth-dam reservoir fixture
    (its ``pool`` drawdown series makes a clear plot).

Both transient captures solve the SMALLEST viable transient model — the earth-dam
reservoir-drawdown fixture on a coarse tri3 mesh, a handful of frames (a
seconds-long solve) — through the real :class:`studio.runners.SeepRunner` transient
path. The heavy solver is owned elsewhere; this stays deliberately trivial.

Run:  python tools/capture_studio_screenshots.py     # regenerate all three PNGs

Exits 0 with a note if PySide6 is not installed (engine-only install — no Studio
layer to capture), mirroring the transient studio smoke test.
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")

try:
    from PySide6.QtWidgets import QApplication, QMessageBox
except Exception:                       # engine-only install — no studio layer
    print("capture_studio_screenshots: PySide6 not installed — skipped.")
    sys.exit(0)

OUT_DIR = os.path.join(REPO_ROOT, "docs", "studio", "images")
# A real reservoir-drawdown transient fixture (materials carry Ss/Sy, the tseep sheet
# carries the pool series + stage times), so nothing needs to be synthesized here.
DAM = os.path.join(REPO_ROOT, "docs/seep/files/xslope_earth_dam_tseep.xlsx")

_app = QApplication.instance() or QApplication([])
# Modal dialogs must never block a headless run.
QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)


def _settle(cycles=12):
    """Pump the event loop so deferred layout + the canvas's debounced render fire."""
    for _ in range(cycles):
        _app.processEvents()
        time.sleep(0.02)


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def capture_run_dialog():
    """Run Seepage dialog, Transient mode — no stage fields, just the caption."""
    from studio.dialogs import RunSeepDialog

    dlg = RunSeepDialog(has_bc2=True, has_tseep=True, defaults={"mode": "transient"})
    dlg.resize(dlg.sizeHint())
    dlg.show()
    _settle()
    out = os.path.join(OUT_DIR, "analysis_run_seep_transient.png")
    dlg.grab().save(out)
    dlg.close()
    return out


def _solve_transient():
    """Smallest viable transient solve: the earth-dam reservoir-drawdown fixture on a
    coarse tri3 mesh. Returns (slope_data, runner bundle {seep_data, frames, ...})."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from studio.runners import SeepRunner

    d = load_slope_data(DAM)
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 16.0, "tri3")
    d["mesh"] = mesh
    runner = SeepRunner(d, {"mode": "transient"})
    bundle, err = {}, {}
    runner.succeeded.connect(lambda b: bundle.update(b))
    runner.failed.connect(lambda msg: err.setdefault("msg", msg))
    _quiet(runner._run_transient, d, mesh)
    if err:
        raise RuntimeError(f"transient solve failed: {err['msg']}")
    return d, bundle


def capture_playbar():
    """Seep · Transient results view with the play bar under a through-flow frame."""
    from studio.transient import TransientSeepView

    _d, bundle = _solve_transient()
    seep_data, frames = bundle["seep_data"], bundle["frames"]
    opts = {"variable": "head", "levels": 20, "flowlines": True,
            "vectors": True, "phreatic": True}
    view = TransientSeepView()
    view.resize(1000, 780)               # ~ the neighboring result-view footprint
    view.set_frames(seep_data, frames, opts_getter=lambda: opts,
                    style_getter=lambda: None, keep_index=False)
    view.show()
    _settle()
    # A mid frame: the reservoir has dropped but through-flow is still strong, so the
    # flow net is fully developed (flow lines present).
    view.set_index(max(view.frame_count() // 3, 1))
    _settle()
    view.canvas._render_current()            # force the raster into the scene
    _settle()
    out = os.path.join(OUT_DIR, "analysis_seep_transient_playbar.png")
    view.grab().save(out)
    view.close()
    return out


def capture_transient_editor():
    """The Transient inputs editor with its live series-vs-time plot populated."""
    from xslope.fileio import load_slope_data
    from studio.editors import TransientDialog

    d = load_slope_data(DAM)
    dlg = TransientDialog(d.get("tseep"), d, None)
    dlg.resize(1080, 620)
    dlg.show()
    _settle()
    dlg._preview.refresh_now()
    _settle()
    dlg._preview.canvas._render_current()     # force the plot raster into the scene
    _settle()
    out = os.path.join(OUT_DIR, "editing_transient_editor.png")
    dlg.grab().save(out)
    dlg.close()
    return out


def main():
    print("capture_studio_screenshots: regenerating transient Studio images")
    for fn in (capture_run_dialog, capture_playbar, capture_transient_editor):
        path = fn()
        print(f"  wrote {os.path.relpath(path, REPO_ROOT)}")
    print("done")


if __name__ == "__main__":
    main()
