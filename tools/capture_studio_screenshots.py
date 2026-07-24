"""Studio-capture pipeline: regenerate the transient-seepage Studio screenshots
that ``docs/studio/analysis.md`` embeds, headlessly (offscreen Qt) so it needs no
display and produces byte-stable, native-quality grabs.

This is the same offscreen-grab practice the reliability / parametric dialog
images were captured with: build the real Studio widget, ``.show()`` it under the
``offscreen`` QPA platform, let the layout settle, then ``QWidget.grab()`` the
widget to a PNG. Nothing here touches a live display or the user's window.

Two captures land in ``docs/studio/images/``:

  * ``analysis_run_seep_transient.png`` — the Run Seepage dialog in **Transient**
    mode with the rapid-drawdown stage fields visible. Built from
    :class:`studio.dialogs.RunSeepDialog` with a tseep sheet present (so the
    Transient choice appears and defaults on) and seeded stage times, sized to its
    ``sizeHint`` to match the neighboring ``analysis_*_dialog.png`` images.
  * ``analysis_seep_transient_playbar.png`` — the **Seep · Transient** results view
    with the play bar (transport buttons, frame slider, ``t =`` readout, Speed
    selector, Set Stage buttons) under a rendered flow-net frame. Built by solving
    the SMALLEST viable transient model (earth_dam1 on a coarse tri3 mesh, a handful
    of one-time-unit frames — a seconds-long solve) through the real
    :class:`studio.runners.SeepRunner` transient path, loading the frames into a
    :class:`studio.transient.TransientSeepView`, and grabbing the whole widget on
    its last (fully developed) frame. The heavy transient solver is owned elsewhere;
    this stays deliberately trivial.

Run:  python tools/capture_studio_screenshots.py     # regenerate both PNGs

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
DAM = os.path.join(REPO_ROOT, "docs/seep/files/xslope_earth_dam1.xlsx")

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
    """Run Seepage dialog, Transient mode, rapid-drawdown stage fields shown."""
    from studio.dialogs import RunSeepDialog

    tseep = {"stage_1": 10.0, "stage_2": 30.0, "duration": 100.0}
    dlg = RunSeepDialog(has_bc2=True, has_tseep=True, tseep=tseep,
                        defaults={"mode": "transient"})
    dlg.resize(dlg.sizeHint())
    dlg.show()
    _settle()
    out = os.path.join(OUT_DIR, "analysis_run_seep_transient.png")
    dlg.grab().save(out)
    dlg.close()
    return out


def _solve_transient():
    """Smallest viable transient solve: earth_dam1 + Ss/Sy on a coarse tri3 mesh,
    a constant-BC series that relaxes from the steady IC so every saved frame is a
    real flow net. Returns the runner bundle {seep_data, frames, ...}."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from studio.runners import SeepRunner

    d = load_slope_data(DAM)
    for m in d["materials"]:
        m["Ss"] = 1e-3
        m["Sy"] = 0.2
    polys = get_material_polygons(d)
    xs = [x for x, _ in d["ground_surface"].coords]
    mesh = _quiet(build_mesh_from_polygons, polys, (max(xs) - min(xs)) / 16.0, "tri3")
    d["mesh"] = mesh
    d["tseep"] = {"times": [], "series": {}, "duration": 4.0, "save_interval": 1.0,
                  "save_times": [], "stage_1": 1.0, "stage_2": 3.0}
    runner = SeepRunner(d, {"mode": "transient", "stage_1": 1.0, "stage_2": 3.0})
    bundle, err = {}, {}
    runner.succeeded.connect(lambda b: bundle.update(b))
    runner.failed.connect(lambda msg: err.setdefault("msg", msg))
    _quiet(runner._run_transient, d, mesh)
    if err:
        raise RuntimeError(f"transient solve failed: {err['msg']}")
    return bundle


def capture_playbar():
    """Seep · Transient results view with the play bar under a rendered frame."""
    from studio.transient import TransientSeepView

    bundle = _solve_transient()
    seep_data, frames = bundle["seep_data"], bundle["frames"]
    opts = {"variable": "head", "levels": 20, "flowlines": True,
            "vectors": True, "phreatic": True}
    view = TransientSeepView()
    view.resize(1000, 780)               # ~ the neighboring result-view footprint
    view.set_frames(seep_data, frames, opts_getter=lambda: opts,
                    style_getter=lambda: None, keep_index=False)
    view.show()
    _settle()
    view.set_index(view.frame_count() - 1)   # last, fully developed frame
    _settle()
    view.canvas._render_current()            # force the raster into the scene
    _settle()
    out = os.path.join(OUT_DIR, "analysis_seep_transient_playbar.png")
    view.grab().save(out)
    view.close()
    return out


def main():
    print("capture_studio_screenshots: regenerating transient Studio images")
    for fn in (capture_run_dialog, capture_playbar):
        path = fn()
        print(f"  wrote {os.path.relpath(path, REPO_ROOT)}")
    print("done")


if __name__ == "__main__":
    main()
