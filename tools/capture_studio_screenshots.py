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

Stability-time captures (the Run dialogs' transient controls):

  * ``analysis_run_lem_dialog.png`` — the Run LEM dialog on a model that carries a
    transient seepage solution, so the **Seepage time** group renders live rather
    than disabled. The solution is a real march of the Johnson reservoir-drawdown
    fixture, so the saved-frame dropdown lists the instants that solve actually
    saved.

The remaining dialogs and dock panels — Build mesh, the steady Run Seepage dialog,
the welcome window, the DXF import wizard, the Global parameters form, the two
Assistant dialogs, the Inputs tree and the two Display panels — are captured the
same way. The Display docks
and the Inputs tree are grabbed off a real :class:`studio.main_window.MainWindow`
built offscreen, so the dock frame, its title bar and the pinned **Styles…** button
are the window's own rather than a stand-in assembled here.

The Assistant settings dialog runs against a THROWAWAY ``QSettings`` file and a
config whose ``api_key`` returns a fixed dummy: the shot must not depend on (or
show anything about) whoever runs it, and it must never touch the real keychain.

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


def _solve_transient(path=None):
    """Smallest viable transient solve: a reservoir-drawdown fixture on a coarse tri3
    mesh. Returns (slope_data, runner bundle {seep_data, frames, transient, ...})."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from studio.runners import SeepRunner

    d = load_slope_data(path or DAM)
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


# --------------------------------------------------------------------------- #
# The Run LEM dialog against a transient seepage solution
# --------------------------------------------------------------------------- #
# The Johnson reservoir-drawdown fixture: a rapid-drawdown LEM model that also
# carries a tseep sheet, so it is the one sample where the Run LEM dialog shows
# every control it has — the seepage-time selector included.
JOHNSON_TSEEP = os.path.join(REPO_ROOT, "docs/lem/files/xslope_johnson_res_rapid.xlsx")


def capture_run_lem_dialog():
    """The Run LEM dialog with the Seepage time group live.

    A stability run against a transient seepage march reads ONE instant of it, so
    the dialog names which. The group is only meaningful with a solution in hand —
    the saved-frame dropdown lists the instants a march actually saved — so this is
    taken after really marching the fixture, on the same coarse tri3 mesh and
    through the same runner the transient captures use.

    The model handed to the dialog is the one the solve worked on — mesh included —
    because that is the state a user is in when they open Run LEM: the seepage
    analysis has just run, and its frames are what the stability run will read. One
    frame is staged into ``seep_u`` first, which is what a run leaves behind, so the
    model checks read a model that is ready rather than one still missing its
    pore-pressure field."""
    from xslope.seep import select_transient_frame_u
    from studio.dialogs import RunLemDialog

    d, bundle = _solve_transient(JOHNSON_TSEEP)
    select_transient_frame_u(d, bundle["transient"])
    dlg = RunLemDialog(defaults={}, slope_data=d, transient=bundle["transient"])
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_run_lem_dialog.png")


# --------------------------------------------------------------------------- #
# Build mesh, steady Run Seepage, DXF import, Global parameters
# --------------------------------------------------------------------------- #
# A steady seepage sample with a second BC set, so the Run Seepage dialog shows the
# BC-set choice the docs describe rather than a lone "Set 1".
SEEP_TWO_BC = os.path.join(REPO_ROOT, "docs/lem/files/xslope_earth_dam_rapid.xlsx")


def capture_build_mesh_dialog():
    """Build mesh, as it opens: auto-size on, feature refinement off."""
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults={})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_build_mesh_dialog.png")


def capture_build_mesh_dialog_refine():
    """Build mesh with an explicit target size and feature refinement on, so the
    Refinement factor is live and Size divisions is greyed."""
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults={"auto_size": False, "target_size": 1.0,
                                    "refine_near_features": True})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_build_mesh_dialog_refine.png")


def capture_unpack_package_dialog():
    """Opening a project package: where its files will be unpacked to.

    The package path is a stand-in that does not exist on this machine, which is
    what the shot wants — a destination that is still free, so the dialog shows the
    default folder and the plain Unpack and Open button rather than the
    already-exists question."""
    from studio.dialogs import UnpackPackageDialog

    dlg = UnpackPackageDialog("/Users/you/Projects/xslope_earth_dam.xslz")
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "usage_unpack_package_dialog.png")


def capture_welcome_dialog():
    """The welcome window as a first launch shows it.

    On a throwaway settings file, for two reasons: the shot must not read (or
    write) the preferences of whoever runs the capture, and a store with nothing in
    it is what a first launch is — so **Don't show this again** appears ticked, the
    way the person in the screenshot would find it."""
    import tempfile
    from PySide6.QtCore import QSettings
    from studio.welcome import WelcomeDialog

    with tempfile.TemporaryDirectory() as tmp:
        settings = QSettings(os.path.join(tmp, "welcome.ini"), QSettings.IniFormat)
        dlg = WelcomeDialog(settings=settings)
        dlg.resize(dlg.sizeHint())
        return _grab(dlg, "overview_welcome_window.png")


def capture_run_seep_dialog():
    """The Run Seepage dialog in Steady mode, on a model with two BC sets.

    Two sets is what makes the BC-set selector a real choice (set 1, set 2, or
    both), which is the thing the surrounding prose describes; a one-set model
    shows a dropdown with nothing to pick."""
    from xslope.fileio import load_slope_data
    from studio.dialogs import RunSeepDialog

    d = load_slope_data(SEEP_TWO_BC)
    dlg = RunSeepDialog(has_bc2=True, has_tseep=False, defaults={"mode": "steady"},
                        slope_data=d)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "analysis_run_seep_dialog.png")


def capture_report_dialog():
    """File → Generate Report…, composed the way a real submittal would be.

    The solutions are stated rather than solved: the Methods list marks which
    bundles were run, so a synthetic pair shows the control doing its job without
    this script running an analysis. Two methods are ticked, which is the picture
    the surrounding prose describes — a report that documents both in full. The
    metadata fields are filled in because empty boxes photograph as a broken
    dialog.

    The contents tree is opened to its full height before the grab. The dialog's
    own size hint sizes the tree to the controls beside it, which on a tree of
    four engine branches leaves the last of them below a scrollbar — a picture of
    a control with most of it cut off. Measured from the rows themselves, so a row
    added to :data:`studio.report_dialog.CONTENT_TREE` grows the shot with it."""
    from PySide6.QtCore import Qt

    from xslope.fileio import load_slope_data
    from studio.report_dialog import ReportDialog

    d = _quiet(load_slope_data, DAM_LEM)
    solutions = {"lem": [{"results": {"FS": 1.532, "method": "spencer"},
                          "method": "spencer"},
                         {"results": {"FS": 1.548, "method": "bishop"},
                          "method": "bishop"}]}
    dlg = ReportDialog(slope_data=d, solutions=solutions, model_path=DAM_LEM,
                       default_method="spencer")
    for item in dlg._method_items():
        if item.data(Qt.UserRole) == "bishop":
            item.setCheckState(Qt.Checked)
    # The output path is overwritten with an illustrative one: the real default
    # is wherever this script happens to run, which is nobody's project folder.
    dlg.path.setText("/Users/you/projects/riverside/xslope_dam_report.docx")
    dlg.title.setText("Riverside Dam — Downstream Slope")
    dlg.project_number.setText("2026-114")
    dlg.organization.setText("Example Engineering")
    dlg.author.setText("N. L. Jones")
    _open_tree(dlg.contents)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "reports_dialog.png")


def _open_tree(tree):
    """Size ``tree`` to its own contents: every expanded row showing at once, and
    every label printed rather than elided."""
    from PySide6.QtWidgets import QTreeWidgetItemIterator

    rows = 0
    it = QTreeWidgetItemIterator(tree, QTreeWidgetItemIterator.NotHidden)
    while it.value():
        rows += tree.visualItemRect(it.value()).height()
        it += 1
    tree.setMinimumHeight(rows + tree.header().height() + 2 * tree.frameWidth())
    tree.setMinimumWidth(tree.sizeHintForColumn(0) + 2 * tree.frameWidth())


def capture_dxf_wizard():
    """The DXF import wizard on a drawing xslope wrote itself.

    The dam sample is exported through :func:`xslope.cad.export_dxf` to a temporary
    file and read straight back with the reader the importer uses, so the layer list
    and the suggested targets are the real ones — material zones, per-material
    profile lines, the search circles and the piezometric line — rather than a
    hand-written fixture."""
    import tempfile
    from xslope.cad import export_dxf, read_dxf_layers, suggest_dxf_target
    from xslope.fileio import load_slope_data
    from studio.dialogs import DxfImportDialog

    d = load_slope_data(DAM_LEM)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "dam.dxf")
        _quiet(export_dxf, d, path)
        layers, _warnings = read_dxf_layers(path)
    dlg = DxfImportDialog(layers, suggest_dxf_target)
    return _grab(dlg, "analysis_dxf_import_wizard.png")


# --------------------------------------------------------------------------- #
# The two-view line editors (reinforcement, piles)
#
# Both views of each, on the rich samples: the LIST view is what the editor opens
# on and what the reference page describes, the TABLE view is the bulk-entry path
# and is shot with EVERY column on show (both usage toggles) and the preview
# collapsed, because its subject is the column run itself — which is the sheet's,
# and which is what makes a block copied out of a worksheet or a tutorial land.
# --------------------------------------------------------------------------- #
REINFORCE = os.path.join(REPO_ROOT, "docs/lem/files/xslope_reinforce.xlsx")
PILES = os.path.join(REPO_ROOT, "docs/lem/files/xslope_piles.xlsx")


def _line_editor(editor_cls, path, mode, size):
    """A two-view line editor on ``path``, switched to ``mode`` and sized for a shot."""
    from xslope.fileio import load_slope_data

    dlg = editor_cls().build(_quiet(load_slope_data, path), None)
    dlg.set_view_mode(mode)
    for cb in dlg._toggles.values():        # every column, LEM and FEM
        if not cb.isChecked():
            cb.setChecked(True)
    dlg.resize(*size)
    if mode == "table":
        # Collapse the preview: the table's columns are the subject here, and the
        # section beside them is already on the list-view shot. The dialog is then
        # widened to whatever the columns need, measured — every column has to be on
        # the shot, and which columns there are is the point of it.
        from PySide6.QtWidgets import QHeaderView

        dlg._table_split.setSizes([size[0], 0])
        table = dlg._table.table
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        _settle()
        chrome = dlg.width() - table.viewport().width()   # everything but the columns
        dlg.resize(table.horizontalHeader().length() + chrome, size[1])
        dlg._table_split.setSizes([dlg.width(), 0])
    else:
        dlg._list_view.list.setCurrentRow(0)
        # Tall enough for the whole form: the groups are the shot's subject, and the
        # last of them (Type, with the Dir and Appl its preset fills) must not be
        # below the scroll. Measured from the form itself rather than guessed.
        _settle()
        scroll = dlg._list_view._form_scroll
        chrome = dlg.height() - scroll.viewport().height()
        dlg.resize(size[0], scroll.widget().sizeHint().height() + chrome)
    return dlg


def capture_reinforcement_editor():
    """Reinforcement editor, list view: the five form groups on the first geogrid
    layer, with the Type group's filled Dir/Appl and the preview's pullout markers."""
    from studio.editors import ReinforcementEditor

    return _grab(_line_editor(ReinforcementEditor, REINFORCE, "list", (1240, 660)),
                 "editing_reinforcement_editor.png")


def capture_reinforcement_table():
    """Reinforcement editor, table view: every column, in the reinforce sheet's
    order, Label first."""
    from studio.editors import ReinforcementEditor

    return _grab(_line_editor(ReinforcementEditor, REINFORCE, "table", (930, 460)),
                 "editing_reinforcement_table.png")


def capture_piles_editor():
    """Piles editor, list view: the four form groups on the first pile."""
    from studio.editors import PilesEditor

    return _grab(_line_editor(PilesEditor, PILES, "list", (1240, 660)),
                 "editing_piles_editor.png")


def capture_piles_table():
    """Piles editor, table view: every column, in the piles sheet's order."""
    from studio.editors import PilesEditor

    return _grab(_line_editor(PilesEditor, PILES, "table", (830, 460)),
                 "editing_piles_table.png")


def capture_seep_bc_editor():
    """Seep BC editor on the earth dam's boundary set, with a demonstration
    infiltration flux added across the crest (the no-flow strip between the pool
    head and the exit face) and selected — so the shot shows the flux value field,
    its points, and the preview's emphasis on the selected boundary."""
    from xslope.fileio import load_slope_data
    from studio.editors import SeepBcEditor

    data = _quiet(load_slope_data, DAM)
    bc = data["seepage_bc"]
    x0 = max(x for x, _ in bc["specified_heads"][0]["coords"])
    x1 = min(x for x, _ in bc["exit_face"])
    crest = [(x, y) for x, y in data["ground_surface"].coords if x0 <= x <= x1]
    bc["specified_fluxes"] = [{"flux": 1.5, "coords": crest}]
    dlg = SeepBcEditor().build(data, None,
                               select=(0, len(bc["specified_heads"])))
    dlg.resize(1080, 560)
    return _grab(dlg, "editing_seep_bc_editor.png")


def capture_global_form():
    """The Global parameters form: the Units / Time / Tension SRF selectors above
    the numeric fields, on a sample that declares a unit system (so the unit-weight
    row carries its unit suffix)."""
    from xslope.fileio import load_slope_data
    from studio.editors import GlobalEditor

    dlg = GlobalEditor().build(load_slope_data(DAM_LEM), None)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "editing_global_form.png")


# --------------------------------------------------------------------------- #
# Assistant dialogs
# --------------------------------------------------------------------------- #

def capture_assistant_settings():
    """The Assistant settings dialog on a throwaway settings file.

    Never the real one: the shot must not depend on which provider the person
    running it happens to use, and the key field must not be filled from the
    keychain. A stub config returns a fixed dummy key, which renders as the same
    row of dots a stored key does.

    ``auto_refresh=False`` for the same reason: the dialog normally enumerates the
    provider's models when it opens, and a screenshot must not depend on a network
    round trip (or make one). The shot therefore shows the offline fallback list,
    which is what a first open looks like anyway."""
    import tempfile
    from PySide6.QtCore import QSettings
    from studio.ai.config import AssistantConfig
    from studio.ai.settings_dialog import AssistantSettingsDialog

    class _CaptureConfig(AssistantConfig):
        def api_key(self, provider):
            return "x" * 48            # a stored key, without reading anyone's

    with tempfile.TemporaryDirectory() as tmp:
        settings = QSettings(os.path.join(tmp, "capture.ini"), QSettings.IniFormat)
        dlg = AssistantSettingsDialog(_CaptureConfig(settings), auto_refresh=False)
        dlg.resize(dlg.sizeHint())
        return _grab(dlg, "assistant_settings.png")


def capture_assistant_confirm():
    """The confirm-before-run dialog, built by the assistant's own ``_confirm_run``.

    That method ends in ``exec()``, so the grab happens from inside a patched
    ``QDialog.exec`` — the dialog on screen is the one the assistant would show, not
    a copy of it assembled here."""
    from PySide6.QtWidgets import QDialog
    from studio.ai.assistant import Assistant

    code = (
        "# Set material properties (undrained clay, phi=0)\n"
        "slope_data['materials'] = [{\n"
        "    'name': 'Clay', 'gamma': 130.0, 'option': 'mc', 'c': 400.0, 'phi': 0.0,\n"
        "    'u': 'none', 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,\n"
        "}]\n"
        "\n"
        "# Geometry: toe (0,10), crest (40,30), base at elev 0\n"
        "slope_data['profile_lines'] = [\n"
        "    {'coords': [(-60.0, 10.0), (0.0, 10.0), (40.0, 30.0), (120.0, 30.0)],\n"
        "     'mat_id': 0}]\n"
        "slope_data['max_depth'] = 0.0   # rigid base 10 ft below toe\n"
        "slope_data['gamma_water'] = 62.4\n"
        "print('material & profile set')\n"
    )
    out = os.path.join(OUT_DIR, "assistant_confirm.png")
    real_exec = QDialog.exec

    def _grab_instead(dlg):
        dlg.show()
        _settle()
        dlg.grab().save(out)
        dlg.close()
        return QDialog.Rejected

    class _Host:                       # _confirm_run only parents the dialog
        _mw = None

    QDialog.exec = _grab_instead
    try:
        Assistant._confirm_run(_Host(), code)
    finally:
        QDialog.exec = real_exec
    return out


# --------------------------------------------------------------------------- #
# Dock panels: the Inputs tree and the two Display panels
# --------------------------------------------------------------------------- #
# These live in a QDockWidget, and the dock's own title bar is part of what the
# docs point at — so they are grabbed off a real MainWindow rather than a
# stand-in dock assembled here. One window serves all three.

_MW = None


def _main_window(path=DAM_LEM):
    """A real Studio window, built offscreen, with ``path`` loaded.

    ``doc.load`` rather than ``open_path``: opening also writes the file into the
    recent-files list, and a screenshot run must not edit anyone's settings."""
    global _MW
    from studio.main_window import MainWindow

    if _MW is None:
        _MW = MainWindow()
        _MW.resize(1400, 900)
        _MW.show()
        _settle()
    _quiet(_MW.doc.load, path)
    _MW._populate_inputs_tree()
    _settle()
    return _MW


DOCK_WIDTH = 290                       # what MainWindow._arrange_docks gives the column


def _grab_dock(dock, name, content_height):
    """Grab a dock sized to show ``content_height`` pixels of its widget.

    The dock's own chrome — title bar, frame — is MEASURED (the difference between
    the dock and its widget once laid out) rather than guessed at, so the same call
    fits a tree and a form panel without a tuned constant per image."""
    dock.resize(DOCK_WIDTH, 400)
    _settle()
    chrome = dock.height() - dock.widget().height()
    dock.resize(DOCK_WIDTH, chrome + int(content_height))
    _settle()
    out = os.path.join(OUT_DIR, name)
    dock.grab().save(out)
    return out


def capture_inputs_tree():
    """The Inputs dock: every input category with its count, the editable ones
    underlined. On the dam sample, which fills most of the list.

    Sized to the rows it actually has: the header plus the bottom of the last item,
    so adding a category grows the image instead of scrolling out of it."""
    mw = _main_window()
    tree = mw.inputs_tree
    last = tree.topLevelItem(tree.topLevelItemCount() - 1)
    height = (tree.header().height() + tree.visualItemRect(last).bottom()
              + 2 * tree.frameWidth() + 4)
    return _grab_dock(mw.inputs_dock, "editing_inputs_tree.png", height)


def capture_mode_selector():
    """The toolbar's analysis-mode strip, shown on FEM.

    Grabbed off the real toolbar rather than a stand-in: the strip's whole point is
    how it sits next to its neighbours, so the crop runs from the **Mode:** label to
    the right edge of the Run button and MEASURES those two positions instead of
    carrying a tuned rectangle."""
    from PySide6.QtCore import QRect
    from PySide6.QtWidgets import QToolBar

    mw = _main_window()
    bar = next(tb for tb in mw.findChildren(QToolBar)
               if tb.objectName() == "main_toolbar")
    mw.set_mode_index(2)               # FEM: Build Mesh and Run FEM both showing
    _settle()
    left = mw.mode_label.geometry().left()
    right = bar.widgetForAction(mw.act_run).geometry().right()
    rect = QRect(left, 0, right - left + 1, bar.height())
    out = os.path.join(OUT_DIR, "interface_mode_selector.png")
    bar.grab(rect).save(out)
    mw.set_mode_index(0)               # leave the shared window as it was found
    return out


def _grab_display_dock(mw, panel, name):
    """Show ``panel`` in the window's own Display dock and grab the dock.

    The panel is pushed onto the dock's stack directly instead of being reached by
    solving and switching tabs: the widget is the real one either way, and a
    screenshot is not a reason to run an SSRM."""
    mw.display_stack.addWidget(panel)
    mw.display_stack.setCurrentWidget(panel)
    mw.styles_btn.setEnabled(True)
    _settle()
    try:
        return _grab_dock(mw.display_dock, name,
                          mw.display_dock.widget().sizeHint().height())
    finally:
        mw.display_stack.removeWidget(panel)
        panel.deleteLater()


def capture_display_panel_seep():
    """The Display panel for a steady seepage solution: the plotted variable,
    contour levels, the flow-net controls (steady-only), and the legend layout."""
    from studio.display_panels import SeepDisplayPanel

    mw = _main_window()
    panel = SeepDisplayPanel(mw.doc.slope_data.get("materials") or [])
    return _grab_display_dock(mw, panel, "analysis_display_panel.png")


def capture_display_dock_fem():
    """The Display panel for FEM results: plot type, the deformation controls, the
    field-state switch, and the mesh/reinforcement overlays."""
    from studio.display_panels import FemResultsDisplayPanel

    mw = _main_window()
    panel = FemResultsDisplayPanel()
    return _grab_display_dock(mw, panel, "interface_display_dock.png")



def analysis_reliability_dialog_lem():
    """The LEM Reliability dialog for the Studio docs — engine, MC controls
    including the convergence stop, and the read-only σ summary. Generated so
    the image tracks the dialog; the FEM variant's image predates this producer
    and its dialog is unchanged by the LEM-only rows."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    import contextlib, io
    from xslope.fileio import load_slope_data
    from studio.dialogs import ReliabilityDialog
    with contextlib.redirect_stdout(io.StringIO()):
        d = load_slope_data(os.path.join(REPO_ROOT, "docs/lem/files/xslope_prob_submerged_KEY.xlsx"))
    dlg = ReliabilityDialog(defaults={"engine": "mc", "method": "spencer",
                                      "num_slices": 40, "search": True},
                            slope_data=d, app_mode="lem")
    dlg.resize(dlg.sizeHint())
    dlg.show()
    app.processEvents()
    pix = dlg.grab()
    out = os.path.join(REPO_ROOT, "docs/studio/images/analysis_reliability_dialog_lem.png")
    pix.save(out)
    print("  wrote", os.path.relpath(out, REPO_ROOT))


def main():
    print("capture_studio_screenshots: regenerating Studio dialog images")
    analysis_reliability_dialog_lem()
    for fn in (capture_run_dialog, capture_playbar, capture_transient_editor,
               capture_polygon_editor, capture_polygon_editor_refine,
               capture_profile_editor, capture_dloads_editor,
               capture_run_fem_dialog, capture_run_lem_preflight,
               capture_run_lem_methods, capture_run_lem_dialog,
               capture_build_mesh_dialog, capture_build_mesh_dialog_refine,
               capture_run_seep_dialog, capture_unpack_package_dialog,
               capture_welcome_dialog,
               capture_report_dialog,
               capture_reinforcement_editor, capture_reinforcement_table,
               capture_piles_editor, capture_piles_table,
               capture_seep_bc_editor,
               capture_dxf_wizard, capture_global_form,
               capture_assistant_settings, capture_assistant_confirm,
               capture_inputs_tree, capture_display_panel_seep,
               capture_display_dock_fem):
        path = fn()
        print(f"  wrote {os.path.relpath(path, REPO_ROOT)}")
    global _MW
    if _MW is not None:                 # closes the assistant's worker thread too
        _MW.close()
        _MW = None
    print("done")


if __name__ == "__main__":
    main()
