"""Studio-capture pipeline: regenerate the Studio dialog/editor screenshots that
``docs/studio/analysis.md`` and ``docs/studio/editing.md`` embed, headlessly
(offscreen Qt) so it needs no display and produces byte-stable, native-quality grabs.

This is the same offscreen-grab practice the reliability / parametric dialog
images were captured with: build the real Studio widget, ``.show()`` it under the
``offscreen`` QPA platform, let the layout settle, then ``QWidget.grab()`` the
widget to a PNG. Nothing here touches a live display or the user's window. Only the
full main-window shots are captured by hand; every dialog image is generated here.

The v21 input captures — the polygon Type/Size row, the profile Size, the
distributed-load Direction and the Run FEM Side BC — run on real sample models with
the v21 field set IN MEMORY for the shot. The corpus files stay untouched: a
screenshot fixture is not a reason to edit a sample model.

Transient captures:

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

v21 input captures:

  * ``editing_polygon_dialog.png`` / ``editing_polygon_dialog_refine.png`` — the
    Polygons editor on the levee model, once on a material zone carrying a local
    mesh Size and once on the ``refine`` region, where the Type dropdown greys out
    Mat ID and the material name.
  * ``editing_geometry_dialog.png`` — the Profile lines editor, with a Size on the
    selected line.
  * ``editing_dloads_editor.png`` — the Distributed loads editor showing two loads
    with different Directions, one normal and one vertical (the preview draws the
    vertical one's arrows straight up).
  * ``analysis_run_fem_dialog.png`` — the Run FEM dialog, now carrying Side BC.

Run:  python tools/capture_studio_screenshots.py     # regenerate every PNG

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
    """Run Seepage dialog, Transient mode — no stage fields, just the caption.

    Built on the earth-dam fixture rather than an empty model, so the model-checks
    box carries what that model actually reports (two anisotropy notes) instead of
    an empty pane."""
    from xslope.fileio import load_slope_data
    from studio.dialogs import RunSeepDialog

    dlg = RunSeepDialog(has_bc2=True, has_tseep=True, defaults={"mode": "transient"},
                        slope_data=load_slope_data(DAM))
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
    # Match the transient display panel's defaults: NO flow lines (a transient state
    # has no flow net), the instantaneous water levels on, and velocity vectors to read
    # the flow direction.
    opts = {"variable": "head", "levels": 20, "flowlines": False,
            "vectors": True, "phreatic": True, "show_bc_levels": True}
    view = TransientSeepView()
    view.resize(1000, 780)               # ~ the neighboring result-view footprint
    view.set_frames(seep_data, frames, opts_getter=lambda: opts,
                    style_getter=lambda: None, keep_index=False)
    view.show()
    _settle()
    # A mid frame: the reservoir has drawn down, so the pool water level sits well below
    # the crest and the velocity vectors show the drainage direction.
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
    dlg.resize(1220, 640)
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


# --------------------------------------------------------------------------- #
# v21 input fields: polygon Type / Size, profile Size, dload Direction, Side BC
# --------------------------------------------------------------------------- #
# A polygon-based model (no profile sheet, so the Polygons editor is the geometry
# editor) and a profile-based one, both small enough to read at dialog size.
LEVEE_POLY = os.path.join(REPO_ROOT, "docs/seep/files/xslope_levee_poly.xlsx")
LAYERS = os.path.join(REPO_ROOT, "docs/lem/files/xslope_eight_layers.xlsx")
EMBANKMENT = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_embankment_mods.xlsx")
# A finite element sample (E/nu set), an LEM sample, and one carrying BOTH surface
# families — the three models the preflight/Run-dialog captures are taken on.
GRIFFITHS = os.path.join(REPO_ROOT, "docs/fem/files/xslope_griffiths2.xlsx")
DAM_LEM = os.path.join(REPO_ROOT, "docs/inputs/slope/xslope_dam.xlsx")
BOTH_FAMILIES = os.path.join(REPO_ROOT,
                             "docs/verification/files/rocscience/vp042.xlsx")


def _grab(dlg, name, settle=True):
    """Show, settle and grab a dialog to ``docs/studio/images/<name>``."""
    dlg.show()
    if settle:
        _settle()
    out = os.path.join(OUT_DIR, name)
    dlg.grab().save(out)
    dlg.close()
    return out


def _polygon_dialog(select):
    """The Polygons editor on the levee model, with the v21 fields populated for the
    shot: a local mesh Size on the levee zone and a refine region over the grout
    curtain's tip (where a real model would want finer elements). Both are set on the
    loaded dict, never written back to the sample file."""
    from xslope.fileio import load_slope_data
    from studio.editors import PolygonEditor

    d = load_slope_data(LEVEE_POLY)
    d["polygons"][0]["size"] = 1.5
    # The grout curtain's TIP — where the flow field turns hardest and a mesh most
    # needs resolving. Straddles the curtain and the foundation on both sides, which
    # is the case a material Size cannot express and a refine region can.
    grout = next(p["polygon"] for p in d["polygons"]
                 if (d["materials"][p["mat_id"]].get("name") or "") == "grout")
    gx0, gy0, gx1, _gy1 = grout.bounds
    pad = 3.0
    d["refine_zones"] = [{"polygon": [(gx0 - pad, gy0), (gx1 + pad, gy0),
                                      (gx1 + pad, gy0 + 3.5), (gx0 - pad, gy0 + 3.5)],
                          "size": 0.6}]
    return PolygonEditor().build(d, None, select=select)


def capture_polygon_editor():
    """Polygons editor on a material zone: Type 'material', Mat ID live, Size set."""
    return _grab(_polygon_dialog(select=0), "editing_polygon_dialog.png")


def capture_polygon_editor_refine():
    """Polygons editor on the refine region: Mat ID and the material name greyed."""
    dlg = _polygon_dialog(select=None)
    dlg.list.setCurrentRow(dlg.list.count() - 1)      # the refine row, listed last
    return _grab(dlg, "editing_polygon_dialog_refine.png")


def capture_profile_editor():
    """Profile lines editor with a Size on the selected line."""
    from xslope.fileio import load_slope_data
    from studio.editors import ProfileEditor

    d = load_slope_data(LAYERS)
    d["profile_lines"][3]["size"] = 2.5
    return _grab(ProfileEditor().build(d, None, select=3),
                 "editing_geometry_dialog.png")


def capture_dloads_editor():
    """Distributed loads editor with one normal load and one vertical."""
    from xslope.fileio import load_slope_data
    from studio.editors import DloadsEditor

    d = load_slope_data(EMBANKMENT)
    # Two loads on this model; make the second a gravity surcharge so the dialog
    # shows both Directions at once and the preview shows the vertical arrows.
    d["dload_dirs"] = ["normal", "vertical"]
    dlg = DloadsEditor().build(d, None, select=(0, 1))
    dlg.show()
    _settle()
    dlg._preview.refresh_now()
    _settle()
    dlg._preview.canvas._render_current()
    _settle()
    out = os.path.join(OUT_DIR, "editing_dloads_editor.png")
    dlg.grab().save(out)
    dlg.close()
    return out


def capture_run_fem_dialog():
    """Run FEM dialog: the v21 Side BC selector and the model checks.

    On a real finite element sample (Griffiths & Lane's slope), so the checks box
    shows a finding the model genuinely carries — the material has no tensile cap,
    which grants it unbounded tension and raises the strength-reduction factor of
    safety with nothing else on screen to show it."""
    from xslope.fileio import load_slope_data
    from studio.dialogs import RunFemDialog

    d = load_slope_data(GRIFFITHS)
    dlg = RunFemDialog(defaults={}, slope_data=d,
                       material_names=[m.get("name", "") for m in d["materials"]])
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_run_fem_dialog.png")


# --------------------------------------------------------------------------- #
# Preflight: the findings a Run dialog shows, and the methods it dims
# --------------------------------------------------------------------------- #
# Both run on real sample models with ONE thing changed in memory for the shot, and
# both changes are the ones a user actually makes: a model whose starting circles
# have not been entered yet, and a model carrying both surface families.

def capture_run_lem_preflight():
    """Run LEM with an error standing: Run refused, and the fix offered as a button.

    The dam sample with its circles sheet emptied — the state every model passes
    through while it is being built. The error is the registry's own sentence, and
    the remedy button beside it is the starting-circle generator."""
    from xslope.fileio import load_slope_data
    from studio.dialogs import RunLemDialog

    d = load_slope_data(DAM_LEM)
    d["circles"] = []
    dlg = RunLemDialog(defaults={}, slope_data=d)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_run_lem_preflight.png")


def capture_run_lem_methods():
    """The method list on a non-circular surface: OMS and Bishop dimmed, with why.

    VP42 carries both a circular and a non-circular surface, so the dialog offers
    the family choice; picking the non-circular one re-filters the method list on
    the spot. The popup is a window of its own, which is what is grabbed here."""
    from xslope.fileio import load_slope_data
    from studio.dialogs import RunLemDialog

    d = load_slope_data(BOTH_FAMILIES)
    dlg = RunLemDialog(defaults={}, slope_data=d)
    dlg.surface.setCurrentIndex(1)                 # non-circular
    dlg.resize(dlg.sizeHint())
    dlg.show()
    _settle()
    dlg.method.showPopup()
    _settle()
    out = os.path.join(OUT_DIR, "analysis_run_lem_methods.png")
    dlg.method.view().window().grab().save(out)
    dlg.close()
    return out


def main():
    print("capture_studio_screenshots: regenerating Studio dialog images")
    for fn in (capture_run_dialog, capture_playbar, capture_transient_editor,
               capture_polygon_editor, capture_polygon_editor_refine,
               capture_profile_editor, capture_dloads_editor,
               capture_run_fem_dialog, capture_run_lem_preflight,
               capture_run_lem_methods):
        path = fn()
        print(f"  wrote {os.path.relpath(path, REPO_ROOT)}")
    print("done")


if __name__ == "__main__":
    main()
