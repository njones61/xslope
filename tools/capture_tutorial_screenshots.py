"""Studio dialog captures for the tutorial pages, headlessly (offscreen Qt).

Same practice as ``tools/capture_studio_screenshots.py``, which captures the
dialogs the Studio reference pages embed: build the real Studio widget, ``.show()``
it under the ``offscreen`` QPA platform, let the layout settle, then
``QWidget.grab()`` it to a PNG. Nothing here touches a live display.

The difference is the subject. The reference pages photograph a dialog to show what
the dialog *is*, and pick whichever sample makes its controls render live. A
tutorial photographs a dialog to show the reader **their own model at a named point
in the build**, so every shot here runs on that tutorial's model and in the state
the step before it leaves behind — the materials editor with this problem's one
material in it, the Run LEM dialog before the circles exist and after.

Full main-window captures are the owner's and are not attempted here; the tutorial
carries a generated placeholder for each (see ``tools/make_tutorial_figures.py``).

Run:  python3 tools/capture_tutorial_screenshots.py           # every shot
      python3 tools/capture_tutorial_screenshots.py lem01     # by name

Exits 0 with a note if PySide6 is not installed (engine-only install — no Studio
layer to capture), mirroring the Studio capture pipeline.
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
    print("capture_tutorial_screenshots: PySide6 not installed — skipped.")
    sys.exit(0)

OUT_DIR = os.path.join(REPO_ROOT, "docs", "tutorials", "images")
LEM01 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_embankment.xlsx")

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


def _grab(dlg, name, settle=True):
    dlg.show()
    if settle:
        _settle()
    out = os.path.join(OUT_DIR, name)
    dlg.grab().save(out)
    dlg.close()
    print("-> %s" % name)
    return out


def _list_view(dlg, width):
    """Put a two-view line editor in its LIST view and size it to the whole form.

    The view is set explicitly because the editor opens on the last one USED (a
    session setting), and the height is measured off the form rather than guessed:
    the shot's subject is the groups, and the last of them must not fall below the
    scroll when a group is added.
    """
    dlg.set_view_mode("list")
    dlg.resize(width, dlg.height())
    dlg.show()
    _settle()
    scroll = dlg._list_view._form_scroll
    chrome = dlg.height() - scroll.viewport().height()
    dlg.resize(width, scroll.widget().sizeHint().height() + chrome)
    _settle()
    return dlg


def _mat_table(dlg, through="u"):
    """Put the materials editor in its TABLE view, wide enough to reach ``through``
    and tall enough for every row.

    Both dimensions are measured off the built table rather than guessed. The width
    is taken to the right edge of the column named ``through`` — the last one the
    tutorial's own material rows fill — because a shot cut off at ``c`` photographs
    a table whose pore-pressure column the reader is being told to set. The height
    is the rows' own, so the last material never falls below the scroll.

    The view is set explicitly for the reason ``_list_view`` documents: the editor
    opens on the last one used.
    """
    from studio.editors import MaterialsEditor

    dlg._set_mode("table")
    dlg.show()
    _settle()
    tbl = dlg._table.table
    hh = tbl.horizontalHeader()
    col = [f.key for f in MaterialsEditor.FIELDS].index(through)
    rows = sum(tbl.rowHeight(r) for r in range(tbl.rowCount()))
    # Twice: the first resize takes the slack out of the layout's stretches, and
    # the chrome measured against a roomy dialog is larger than the chrome of the
    # tight one it produces. The second pass measures the tight one.
    for _ in range(2):
        w = (hh.sectionPosition(col) + hh.sectionSize(col)
             + (dlg.width() - tbl.viewport().width()))
        dlg.resize(w, rows + (dlg.height() - tbl.viewport().height()))
        _settle()
    return dlg


def _line_table(dlg, through):
    """Put a two-view line editor in its TABLE view, sized to reach ``through``.

    The line editors' own default is the LIST view and the last view used is a
    session setting, so the mode is set explicitly — the same reason ``_mat_table``
    gives. Width and height are measured off the built table for the same reason
    too: a reinforcement table cut off before its last column photographs a row
    whose columns the page is telling the reader to fill.
    """
    dlg._set_mode("table")
    dlg.show()
    _settle()
    tbl = dlg._table.table
    hh = tbl.horizontalHeader()
    keys = [f.key for f in dlg._fields]
    col = keys.index(through)
    rows = sum(tbl.rowHeight(r) for r in range(tbl.rowCount()))
    want = rows + tbl.horizontalHeader().height()
    # Height is grown by the measured deficit rather than computed once: the table
    # shares a splitter with the section preview, so a taller dialog does not hand
    # the table all of the extra — and a table one row short photographs a file
    # whose last line the reader cannot see.
    from PySide6.QtWidgets import QSplitter

    for _ in range(4):
        w = (hh.sectionPosition(col) + hh.sectionSize(col)
             + (dlg.width() - tbl.viewport().width()))
        deficit = max(0, want - tbl.viewport().height())
        dlg.resize(w, dlg.height() + deficit)
        # The splitter keeps its own proportions through a resize, so the extra
        # height lands on the preview unless the table's pane is asked for it.
        for sp in dlg.findChildren(QSplitter):
            sizes = sp.sizes()
            if len(sizes) == 2 and sp.isAncestorOf(tbl) and deficit:
                sp.setSizes([sizes[0] + deficit, max(1, sizes[1])])
        _settle()
    return dlg


@contextlib.contextmanager
def _app_defaults():
    """Run a capture against the app's OWN defaults, not this machine's stored ones.

    The usage toggles above an editor's table — LEM / FEM / Seepage / Reliability —
    decide WHICH COLUMNS the table view shows, and each is remembered in QSettings
    the moment anyone clicks it. A capture that simply opened the editor would
    therefore photograph whatever the person running the producer last had ticked,
    and the material rows a tutorial teaches would gain or lose columns between one
    machine and the next. Every ``editor_toggles`` key is removed for the length of
    the run — each dialog then falls back to the default coded beside it — and the
    machine's own settings are put back afterwards, whatever happens.
    """
    from PySide6.QtCore import QSettings

    settings = QSettings("XSlope", "XSlope Studio")
    settings.beginGroup("editor_toggles")
    stashed = {k: settings.value(k) for k in settings.allKeys()}
    for key in stashed:
        settings.remove(key)
    settings.endGroup()
    settings.sync()
    try:
        yield
    finally:
        settings.beginGroup("editor_toggles")
        for key, value in stashed.items():
            settings.setValue(key, value)
        settings.endGroup()
        settings.sync()


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


# --------------------------------------------------------------------------- #
# LEM-1 — Simple Embankment
# --------------------------------------------------------------------------- #
def lem01_materials():
    """The materials editor, **list view**, on this problem's one material.

    List view rather than the table it opens on: with a single material the table
    is one row of a very wide sheet, while the list view puts the strength envelope
    beside the numbers — and on a φ = 0 material that plot is a horizontal line at
    τ = c, which is the point the tutorial's step is making.
    """
    from studio.editors import MaterialsEditor

    dlg = _lem_only(MaterialsEditor().build(_load(LEM01), None))
    dlg._set_mode("list")
    dlg.resize(1180, 720)
    return _grab(dlg, "lem01_studio_materials.png")


def lem01_profile():
    """The profile-lines editor on the finished geometry: three vertices, and the
    preview drawing the section they make."""
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM01), None, select=0)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem01_studio_profile.png")


def lem01_run_lem_no_surface():
    """Run LEM on the model as the geometry steps leave it — no failure surface yet.

    The circles are dropped from the loaded model (never from the file): the model
    checks then report the missing surface as an error, **Run** is disabled, and the
    remedy button beside the finding is the starting-circle generator the tutorial's
    next step presses.
    """
    from studio.dialogs import RunLemDialog

    d = _load(LEM01)
    d["circles"] = []
    d["circular"] = False
    dlg = RunLemDialog(defaults={}, slope_data=d)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem01_studio_run_lem_no_surface.png")


def lem01_circles():
    """The circles editor just after **Generate starting circles…** — what the
    reader audits.

    The button is *pressed* here rather than the rows being pre-loaded, because the
    summary line the step teaches the reader to read ("one candidate dropped for not
    daylighting inside the model") only exists as the answer to a press. The model
    arrives the way the geometry steps leave it — no circles — so the generator's
    table is empty and it fills it without asking, which is also the flow the page
    describes.
    """
    from studio.editors import CirclesEditor

    d = _load(LEM01)
    d["circles"] = []
    d["circular"] = False
    dlg = CirclesEditor().build(d, None)
    with contextlib.redirect_stdout(io.StringIO()):
        dlg._run_generate()
    return _grab(dlg, "lem01_studio_circles.png")


def lem01_run_lem():
    """The Run LEM dialog on the complete model, filled in the way the step beside
    it dictates: Spencer, an auto search, 40 slices — and a clean model-checks
    column.

    The dialog is photographed in the state the reader's own dialog is in when
    they press **Run**, not in the state it opens in: the page's numbered list is
    read against this figure, and a capture showing the defaults would have the
    reader checking their choices against a picture of the choices not yet made.
    Surface takes no default — this model defines circles and no non-circular
    surface, so the row is the fixed label the step describes.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM01))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem01_studio_run_lem.png")


SHOTS = {
    "lem01_materials": lem01_materials,
    "lem01_profile": lem01_profile,
    "lem01_run_lem": lem01_run_lem,
    "lem01_circles": lem01_circles,
    # lem01_run_lem_no_surface retired with the preflight-remedy flow the page
    # no longer teaches.
}



def lem01_canvas():
    """The main window at the end of the Studio path's step 3: geometry entered,
    no circles yet. A full-window offscreen grab — the deferred canvas render
    needs a forced synchronous kick (no real screen, no paint events), and the
    only loss vs a hand capture is the macOS title bar. The assistant dock is
    hidden (the Studio path does not use it) and the log cleared."""
    from studio.main_window import MainWindow
    win = MainWindow()
    win.resize(1600, 1000)
    win.open_path(os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_embankment.xlsx"))
    win.doc.slope_data["circles"] = []
    win._populate_inputs_tree()
    for name in ("assistant_dock", "chat_dock", "ai_dock"):
        d = getattr(win, name, None)
        if d is not None:
            d.hide()
    win.show()
    _settle()
    win.canvas.render_inputs(win.doc.slope_data)
    _settle()
    win.canvas._render_timer.stop()
    win.canvas._render_current()
    win.log.clear()
    _settle()
    pix = win.grab()
    out = os.path.join(OUT_DIR, "lem01_studio_canvas.png")
    pix.save(out)
    print("-> lem01_studio_canvas.png  (%dx%d, offscreen main window)"
          % (pix.width(), pix.height()))
    win.close()


SHOTS["lem01_canvas"] = lem01_canvas


# --------------------------------------------------------------------------- #
# Tutorial 0 — Building Models Three Ways
# --------------------------------------------------------------------------- #
def t0_studio_window():
    """The whole Studio window, with both of the paths it hosts on screen at once.

    Tutorial 0's claim is that the editors and the assistant are two ways into one
    open project, so this shot keeps the Assistant dock rather than hiding it the
    way LEM-1's canvas figure does — the Inputs tree on the left, the section in the
    middle, the assistant on the right, all of one document. Same offscreen
    main-window technique as ``lem01_canvas`` (forced synchronous canvas render, no
    real screen), on LEM-1's model, so the reader recognises the section when they
    reach their first build.

    The dock's caption names the active provider and model, which it reads from the
    machine's stored settings — so an uncaptioned capture photographs whatever the
    person running it happens to have selected, which for this figure meant a model
    a reader could neither find in the list nor reach. The capture therefore pins
    the selection to the shipped default (Claude, and the first entry of the
    provider's static list) for the length of the grab and puts the machine's own
    settings back afterwards, whatever happens. The pinned model is read from the
    provider table rather than spelled here, so the figure follows the default if it
    changes.

    Pinned rather than deleted: with no model stored, the dock falls back to the
    recommendations manifest when one has been fetched into these same settings,
    which is again machine-dependent.
    """
    from PySide6.QtCore import QSettings

    from studio.ai.config import PROVIDERS
    from studio.main_window import MainWindow

    #: The settings the app itself uses — the real ones, which is why every write
    #: below is undone in the finally block.
    settings = QSettings("XSlope", "XSlope Studio")
    pinned = {"ai/provider": "anthropic",
              "ai/model/anthropic": PROVIDERS["anthropic"]["models"][0]}
    # `None` records "this key was not set", which must be restored by REMOVING it
    # again rather than by writing an empty string.
    stashed = {k: (settings.value(k) if settings.contains(k) else None)
               for k in pinned}
    win = None
    try:
        for key, value in pinned.items():
            settings.setValue(key, value)
        settings.sync()

        win = MainWindow()
        win.resize(1600, 1000)
        win.open_path(LEM01)
        win.show()
        _settle()
        win.canvas.render_inputs(win.doc.slope_data)
        _settle()
        win.canvas._render_timer.stop()
        win.canvas._render_current()
        win.log.clear()
        _settle()
        pix = win.grab()
        out = os.path.join(OUT_DIR, "t0_studio_window.png")
        pix.save(out)
        print("-> t0_studio_window.png  (%dx%d, offscreen main window, %s)"
              % (pix.width(), pix.height(), pinned["ai/model/anthropic"]))
        return out
    finally:
        for key, value in stashed.items():
            if value is None:
                settings.remove(key)
            else:
                settings.setValue(key, value)
        settings.sync()
        if win is not None:
            win.close()


SHOTS["t0_studio_window"] = t0_studio_window


def t0_unpack_exists():
    """The package-open dialog on a destination that is ALREADY THERE.

    The two states are mutually exclusive — ``UnpackPackageDialog._refresh`` hides
    **Unpack and Open** exactly when **Open Existing** and **Extract Fresh** appear
    — so the Studio reference page's capture (a free destination) cannot illustrate
    the question this state asks, and the tutorial shows both.

    The path in the figure is the same stand-in its sibling uses, so the two read as
    one dialog in two states rather than two dialogs. Making that path *exist* is
    the whole difficulty: the dialog asks the filesystem. So the filesystem's answer
    is stubbed for the one path, rather than the dialog being driven into the state
    by hand — the visibility rule, the status line and the numbered fresh folder are
    then all the dialog's own work, which is what the tutorial's prose describes.
    The stub answers True for that path ONLY: ``_fresh_dest`` probes ``-2``, ``-3``
    … for the first free name and would not terminate against a blanket True.
    """
    from studio import dialogs
    from studio.dialogs import UnpackPackageDialog

    dest = "/Users/you/Projects/xslope_earth_dam"
    dlg = UnpackPackageDialog(dest + ".xslz")
    real_exists = dialogs.os.path.exists
    dialogs.os.path.exists = lambda p: p == dest or real_exists(p)
    try:
        dlg._refresh()
        dlg.resize(dlg.sizeHint())
        return _grab(dlg, "t0_unpack_exists.png")
    finally:
        dialogs.os.path.exists = real_exists


SHOTS["t0_unpack_exists"] = t0_unpack_exists


# --------------------------------------------------------------------------- #
# LEM-2 — Loads on the Crest
# --------------------------------------------------------------------------- #
LEM02 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_crest_surcharge.xlsx")


def lem02_dloads():
    """The distributed-loads editor holding the crest surcharge.

    The model arrives with the load already in it, because this is the state the
    Studio path's step *ends* in and the figure is what the reader checks their
    own screen against: the two points, the intensity in both, and the Direction
    the level crest makes irrelevant. The preview beside the table draws the load
    on the section, which is where a mistyped X shows up.
    """
    from studio.editors import DloadsEditor

    dlg = DloadsEditor().build(_load(LEM02), None)
    dlg.resize(1000, 620)
    return _grab(dlg, "lem02_studio_dloads.png")


def lem02_lloads():
    """The line-loads editor with the surcharge re-stated as a single force.

    Built on the model with its distributed load REMOVED, since the page's line
    load replaces the surcharge rather than joining it — the preview would
    otherwise draw a section carrying twice the force the step describes.
    """
    from studio.editors import LineLoadsEditor

    d = _load(LEM02)
    d["dloads"], d["dload_dirs"] = [], []
    d["line_loads"] = [{"x": 30.0, "y": 20.0, "P": 7500.0, "angle": -90.0,
                        "label": "footing"}]
    dlg = LineLoadsEditor().build(d, None)
    dlg.resize(1000, 620)
    return _grab(dlg, "lem02_studio_lloads.png")


def lem02_parametric():
    """The Parametric dialog in **Design** mode, set up for the page's sweep.

    Every control is driven the way a reader drives it rather than pre-filled:
    the mode is chosen, then the material and the property, and only then the
    bounds — because the dialog seeds From/To from the property the picker lands
    on, and a dialog whose fields were written before its parameter was picked
    photographs a range the reader could not have got to. The **Sweeping** echo
    reads the picked parameter for the same reason: it is the picker's answer,
    not a caption.
    """
    from studio.dialogs import SensitivityDialog

    dlg = SensitivityDialog(defaults={"mode": "design", "method": "spencer",
                                      "num_slices": 40, "search": True},
                            slope_data=_load(LEM02))
    for combo, text in ((dlg.mode, "Design"), (dlg.material, "soil")):
        i = next(n for n in range(combo.count()) if text in combo.itemText(n))
        combo.setCurrentIndex(i)
    i = next(n for n in range(dlg.prop.count()) if dlg.prop.itemText(n) == "c")
    dlg.prop.setCurrentIndex(i)
    dlg.d_from.setValue(500.0)
    dlg.d_to.setValue(1200.0)
    dlg.d_steps.setValue(8)
    dlg.d_target.setValue(1.5)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem02_studio_parametric.png")


# --------------------------------------------------------------------------- #
# LEM-3 — A Layered Slope
# --------------------------------------------------------------------------- #
LEM03 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_mult_layers.xlsx")


def lem03_materials():
    """The materials editor, **table view**, on this problem's two materials.

    Table view rather than the list view LEM-1 photographs: with two materials the
    subject is the pair — the order that fixes the Mat IDs the profile lines
    reference, and one row's numbers read against the other's. The list view shows
    one material at a time, which is the wrong shape for that comparison.
    """
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_lem_only(MaterialsEditor().build(_load(LEM03), None))),
                 "lem03_studio_materials.png")


def lem03_circles():
    """The circles editor just after **Generate starting circles…** on layered ground.

    The button is pressed rather than the rows pre-loaded, for the reason LEM-1's
    shot documents: the summary line is the answer to a press. The model arrives
    with its circles dropped, the state the geometry steps leave — so what the
    table holds is the generator's own set, three circles on a two-layer section,
    and the tutorial's audit is read against it.
    """
    from studio.editors import CirclesEditor

    d = _load(LEM03)
    d["circles"] = []
    d["circular"] = False
    dlg = CirclesEditor().build(d, None)
    with contextlib.redirect_stdout(io.StringIO()):
        dlg._run_generate()
    return _grab(dlg, "lem03_studio_circles.png")


SHOTS["lem03_materials"] = lem03_materials
SHOTS["lem03_circles"] = lem03_circles


# --------------------------------------------------------------------------- #
# LEM-5 — A Weak Layer, Non-Circular
# --------------------------------------------------------------------------- #
LEM05 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_noncircular.xlsx")


def lem05_noncirc():
    """The non-circular editor holding the four points, as the step leaves it.

    The model arrives with the surface already in it, because this is the state
    the Studio path's step ends in: the vertex table, the Movement column, and the
    preview drawing the polyline across the section so the reader can see the
    track lying inside the clay seam before anything is run.
    """
    from studio.editors import NonCircEditor

    dlg = NonCircEditor().build(_load(LEM05), None)
    dlg.resize(1100, 640)
    return _grab(dlg, "lem05_studio_noncirc.png")


def lem05_generate():
    """The same editor just after **Generate from the weak zone…** is pressed.

    The table arrives empty — the state a reader who has entered no surface is in
    — so the generator fills it without asking, and the summary line under the
    button is the answer to the press rather than a caption: which zone it chose,
    the strengths it compared to choose it, and the ramp angles it built. Pressing
    the button is the whole point of the shot; pre-loading the rows would lose
    the line the page teaches the reader to read.
    """
    from studio.editors import NonCircEditor

    d = _load(LEM05)
    d["non_circ"] = []
    dlg = NonCircEditor().build(d, None)
    with contextlib.redirect_stdout(io.StringIO()):
        dlg._run_generate()
    dlg.resize(1100, 640)
    return _grab(dlg, "lem05_studio_generate.png")


def lem05_run_lem():
    """Run LEM on the finished model, filled the way the step beside it dictates.

    **Surface** is the row that differs from every earlier tutorial: this model
    defines a non-circular surface and no circles, so the selector is a fixed
    label reading *Non-circular*, the mirror of the fixed *Circular* label LEM-2's
    model produces.

    **Analysis** is *Auto search*, which is the first analysis the page runs: the
    four points the model carries are where a search starts, not what it answers.
    The page's held-surface comparisons switch the same row to *Single surface*
    and say so.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM05))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem05_studio_run_lem.png")


SHOTS["lem05_noncirc"] = lem05_noncirc
SHOTS["lem05_generate"] = lem05_generate
SHOTS["lem05_run_lem"] = lem05_run_lem

SHOTS["lem02_dloads"] = lem02_dloads
SHOTS["lem02_lloads"] = lem02_lloads
SHOTS["lem02_parametric"] = lem02_parametric


# --------------------------------------------------------------------------- #
# LEM-4 — Water in the Slope
# --------------------------------------------------------------------------- #
LEM04 = os.path.join(REPO_ROOT,
                     "docs/lem/files/xslope_method_slices_problem.xlsx")


def lem04_piezo():
    """The piezometric-lines editor holding the eight points, as the step leaves it.

    The model arrives with the line already in it, because this is the state the
    Studio path's step ends in: the point table on the **Line 1** tab, the empty
    **Line 2 (rapid drawdown)** tab beside it, and the preview drawing the line
    falling across the section — which is where a mistyped elevation shows up,
    as a point off the smooth descent from 80 to 40.
    """
    from studio.editors import PiezoEditor

    dlg = PiezoEditor().build(_load(LEM04), None)
    dlg.resize(1100, 640)
    return _grab(dlg, "lem04_studio_piezo.png")


SHOTS["lem04_piezo"] = lem04_piezo


# --------------------------------------------------------------------------- #
# LEM-6 — Polygon Geometry
# --------------------------------------------------------------------------- #
LEM06 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_sloping_bottom.xlsx")


def lem06_polygons():
    """The polygons editor on the foundation zone, as LEM-6's Studio step leaves it.

    The second zone is selected rather than the first, because the foundation's
    ring is the one this page is about: its four vertices are where the dipping
    base is entered, and the preview draws that base under the section while the
    other zone is dimmed behind it. The editor is built on a model that already
    carries both zones — Studio edits polygons on a polygon-based project, so
    this is the state a reader reaches by opening one, which is what the step
    says.
    """
    from studio.editors import PolygonEditor

    dlg = PolygonEditor().build(_load(LEM06), None, select=1)
    dlg.resize(1400, 560)
    return _grab(dlg, "lem06_studio_polygons.png")


SHOTS["lem06_polygons"] = lem06_polygons


# --------------------------------------------------------------------------- #
# LEM-8 — A Reinforced Slope
# --------------------------------------------------------------------------- #
LEM08 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_reinforce.xlsx")




def _all_usages(dlg):
    """Tick every usage toggle except Reliability — the app's own default.

    Explicit for the same reason the paste tests pin it: toggle state persists
    per editor, so a dry page's _lem_only capture earlier in the same batch
    would otherwise leak into a wet page's shot and hide the columns its
    14-column paste needs."""
    for tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(tag != "rel")
    return dlg


def _lem_only(dlg):
    """Untick every non-LEM usage toggle so the capture shows the LEM problem.

    The owner's rule for the LEM tutorials: an LEM page's editor capture carries
    no FEM or seepage columns. The toggles are real checkboxes, so unticking
    them is the same action a reader takes, and the toggle bar in the capture
    shows the state that produced it.
    """
    toggles = getattr(dlg, "_toggles", None) or {}
    for tag, cb in toggles.items():
        cb.setChecked(tag == "lem")
    return dlg

def lem08_reinforcement():
    """The reinforcement editor on the six geogrid layers, in its list view.

    The list view rather than the table, because it is the view the editor opens
    on and the one that groups a line the way the page teaches it — Identity,
    Geometry, Capacity, Anchorage, Type — beside a preview that draws the pullout
    breakpoints the capacity envelope is built from. The first line is selected,
    which is the line the page reads the envelope on: the one the critical
    surface passes beneath.
    """
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(LEM08), None)
    return _grab(_list_view(_lem_only(dlg), 1400), "lem08_studio_reinforcement.png")


SHOTS["lem08_reinforcement"] = lem08_reinforcement

def lem08_reinforcement_table():
    """The reinforcement editor's table view on the six geogrid layers — the
    view the page teaches first, because the two Problem-section blocks paste
    straight into it. LEM columns only, per the page's toggles."""
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(LEM08), None)
    _lem_only(dlg)
    # The preview stacks below the table, so the width is just the columns:
    # measured out to Spacing, the last LEM column.
    return _grab(_line_table(dlg, through="spacing"),
                 "lem08_studio_reinforcement_table.png")


SHOTS["lem08_reinforcement_table"] = lem08_reinforcement_table



# --------------------------------------------------------------------------- #
# LEM-9 — A Tieback Wall
# --------------------------------------------------------------------------- #
LEM09 = os.path.join(REPO_ROOT,
                     "docs/verification/files/rocscience/vp049.xlsx")


def lem09_reinforcement_table():
    """The reinforcement editor's table view on the two tiebacks — the view the
    page teaches first, because the two Problem-section blocks paste straight
    into it. LEM columns only, per the page's toggles."""
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(LEM09), None)
    _lem_only(dlg)
    return _grab(_line_table(dlg, through="spacing"),
                 "lem09_studio_reinforcement_table.png")


SHOTS["lem09_reinforcement_table"] = lem09_reinforcement_table


def lem09_reinforcement():
    """The reinforcement editor on the upper tieback, in its list view.

    The same editor LEM-8 photographs, on a row whose Type is ``anchor`` rather
    than blank: the shot's subject is the Type group at the bottom of the form,
    where the preset has filled Dir with ``axial`` and Appl with ``active``, and
    the Anchorage group above it, where the two pullout lengths are 0 at the wall
    and the bond length at the far end. The first line is selected because it is
    the one the page reads its capacity off.
    """
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(LEM09), None)
    return _grab(_list_view(_lem_only(dlg), 1400), "lem09_studio_reinforcement.png")


def lem09_piles_table():
    """The piles editor's table view on the soldier pile — the view the page
    teaches first: the worksheet's columns minus the derived θp, taking the
    Problem-section row as two pieces. LEM columns only, per the page's
    toggles."""
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(LEM09), None)
    _lem_only(dlg)
    return _grab(_line_table(dlg, through="M_cap"),
                 "lem09_studio_piles_table.png")


SHOTS["lem09_piles_table"] = lem09_piles_table


def lem09_piles():
    """The piles editor on the soldier pile, in its list view.

    Built on LEM-9's model, so the form holds the pile the page enters — the
    vertical axis at the wall face, the shear force per foot of wall, and the
    diameter and spacing beside it — rather than an empty row.
    """
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(LEM09), None)
    return _grab(_list_view(_lem_only(dlg), 1400), "lem09_studio_piles.png")


SHOTS["lem09_reinforcement"] = lem09_reinforcement
SHOTS["lem09_piles"] = lem09_piles


# --------------------------------------------------------------------------- #
# Studio-path parity — a completed dialog for every build step
#
# The Excel path shows the finished worksheet under each of its steps; the Studio
# path shows the finished dialog under each of its. The shots below fill in the
# steps that carried none, and every one of them photographs the editor in the
# state its own step ends in, on that tutorial's own model.
#
# The Run LEM dialog is captured per page for the same reason: the numbered choices
# under "Running the analysis" are read against a picture of the choices MADE, so
# each page's shot carries that page's method and that page's surface row — which
# is a fixed label on a model defining one surface family and a selector on a model
# defining both.
# --------------------------------------------------------------------------- #
def _run_lem(model, name, method="spencer", analysis="auto_search"):
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": method, "analysis": analysis,
                                 "num_slices": 40},
                       slope_data=_load(model))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, name)


# --------------------------------------------------------------------------- #
# LEM-10 — Finding the Global Minimum (open-and-run: the file is the start)
# --------------------------------------------------------------------------- #
LEM10 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_mult_min_KEY.xlsx")


def lem10_run_lem():
    """The Run LEM dialog on the loaded file, in the state the first run is made
    in: Spencer, auto search, 40 slices — nothing else touched, because the
    whole point of the first run is the file exactly as downloaded."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM10))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem10_studio_run_lem.png")


def lem10_run_lem_filtered():
    """The same dialog with the surficial filter on at 5 ft — the state Part A's
    filter step tells the reader to put it in. Grid search stays unticked; that
    tool belongs to Part B."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40,
                                 "min_slip_depth": 5.0},
                       slope_data=_load(LEM10))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem10_studio_run_lem_filtered.png")


def lem10_circles():
    """The circles editor after the tutorial's one edit: the file's deep circle
    replaced by the embankment circle (16.875, 30, Depth 0) — tangent to the top
    of the foundation, the seed the second search runs from."""
    from studio.editors import CirclesEditor

    d = _load(LEM10)
    shallow = dict(d["circles"][0])
    shallow.update({"Xo": 16.875, "Yo": 30.0, "Depth": 0.0, "R": None})
    dlg = CirclesEditor().build(dict(d, circles=[shallow]), None)
    return _grab(dlg, "lem10_studio_circles.png")


LEM10_JB = os.path.join(REPO_ROOT, "docs/lem/files/xslope_james_bay.xlsx")


def lem10_jb_run_lem():
    """Part B's first run: the dyke's single seed, Spencer auto search with the
    surficial filter at the metric 2 m — the state the step dictates."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40, "min_slip_depth": 2.0},
                       slope_data=_load(LEM10_JB))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem10_jb_run_lem.png")


def lem10_jb_run_lem_grid():
    """Part B's second run: the same dialog with Grid search ticked."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40, "min_slip_depth": 2.0,
                                 "grid_seed": True},
                       slope_data=_load(LEM10_JB))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem10_jb_run_lem_grid.png")


SHOTS["lem10_run_lem"] = lem10_run_lem
SHOTS["lem10_run_lem_filtered"] = lem10_run_lem_filtered
SHOTS["lem10_circles"] = lem10_circles
SHOTS["lem10_jb_run_lem"] = lem10_jb_run_lem
SHOTS["lem10_jb_run_lem_grid"] = lem10_jb_run_lem_grid


def lem01_global():
    """The global-parameters form as LEM-1's first Studio step leaves it.

    The model carries the answers the step dictates — Imperial, water at 62.4, and
    the tension crack and seismic fields at 0 — so the shot is the form filled,
    which is what the reader checks their own against. The unit-weight field is the
    one Units autofills, and it is photographed holding the autofilled value rather
    than a typed one.
    """
    from studio.editors import GlobalEditor

    dlg = GlobalEditor().build(_load(LEM01), None)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem01_studio_global.png")


def lem02_run_lem():
    return _run_lem(LEM02, "lem02_studio_run_lem.png")


def lem03_profile():
    """The profile-lines editor with both of LEM-3's lines entered, the CONTACT
    selected.

    The second line rather than the first: it is the one the step's paragraph is
    about — the top of the foundation, and by the top-of-a-layer rule everything
    beneath it down to the maximum depth — and selecting it puts its two vertices
    in the table while the preview draws it running the full width of the section
    under the embankment.
    """
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM03), None, select=1)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem03_studio_profile.png")


def lem03_run_lem():
    return _run_lem(LEM03, "lem03_studio_run_lem.png")


def lem04_materials():
    """The materials editor, table view, on LEM-4's three soils.

    Table view for the reason LEM-3's shot documents — the row order is what fixes
    the Mat IDs, and three rows are read against each other. This is the first
    tutorial whose materials fill BOTH unit-weight columns and set ``u`` to
    ``piezo``, so the frame runs through the ``u`` column rather than stopping at
    the strength values.
    """
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_all_usages(MaterialsEditor().build(_load(LEM04), None))),
                 "lem04_studio_materials.png")


def lem04_profile():
    """The profile-lines editor on LEM-4's three lines, the first selected — the
    ground surface, whose vertices carry the crest, the face and the toe."""
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM04), None, select=0)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem04_studio_profile.png")


def lem04_circles():
    """The circles editor holding LEM-4's single starting circle.

    Typed rather than generated: this page's step enters one circle by hand, and
    the preview below the table is where its Depth is read — an elevation, drawn
    bottoming out below the toe.
    """
    from studio.editors import CirclesEditor

    return _grab(CirclesEditor().build(_load(LEM04), None),
                 "lem04_studio_circles.png")


def lem04_run_lem():
    return _run_lem(LEM04, "lem04_studio_run_lem.png")


def lem05_materials():
    """The materials editor, table view, on LEM-5's four soils — the seam of soft
    clay among them, told from the sand fill by its cohesion."""
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_all_usages(MaterialsEditor().build(_load(LEM05), None))),
                 "lem05_studio_materials.png")


def lem05_profile():
    """The profile-lines editor on LEM-5's four lines, the clay seam selected.

    Line 3 is the top of the soft clay — the layer the whole page is about — so it
    is the one selected: its two vertices are in the table, and the preview draws
    the seam inside the section with the other three lines behind it.
    """
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM05), None, select=2)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem05_studio_profile.png")


def lem05_piezo():
    """The piezometric-lines editor holding LEM-5's two-point water table.

    The first tab is selected explicitly: the second is the rapid-drawdown line,
    which this model leaves empty, and the shot's subject is the pair — a flat
    line at elevation −2 on Line 1, nothing on Line 2.
    """
    from studio.editors import PiezoEditor

    dlg = PiezoEditor().build(_load(LEM05), None)
    dlg._tabs.setCurrentIndex(0)
    dlg.resize(1100, 560)
    return _grab(dlg, "lem05_studio_piezo.png")


def lem06_materials():
    """The materials editor, table view, on LEM-6's two soils — the pair the
    polygon rings reference by the Mat ID their row order fixes."""
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_lem_only(MaterialsEditor().build(_load(LEM06), None))),
                 "lem06_studio_materials.png")


def lem06_circles():
    """The circles editor as LEM-6's step leaves it: the generator's pair with the
    toe circle's Depth typed over.

    The generator cannot propose the deep seed on a dipping base — there is no
    single base elevation to target — so the step presses it and then edits one
    cell. The model as delivered carries exactly that result, which is why the shot
    is built on the model rather than on a press: the table holds the edited value
    the step ends with, not the value the button produced.
    """
    from studio.editors import CirclesEditor

    return _grab(CirclesEditor().build(_load(LEM06), None),
                 "lem06_studio_circles.png")


def lem06_run_lem():
    return _run_lem(LEM06, "lem06_studio_run_lem.png")


def lem08_materials():
    """The materials editor, table view, on LEM-8's two soils — the same γ and φ on
    both rows, told apart by the cohesion of the face band."""
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_lem_only(MaterialsEditor().build(_load(LEM08), None))),
                 "lem08_studio_materials.png")


def lem08_profile():
    """The profile-lines editor on LEM-8's two lines, the face band selected.

    Line 1 is the band of cohesive fill along the face — three vertices, the toe,
    the top of the 1.25:1 face and 2 ft of crest — and the preview draws it against
    line 2 running out in front of the toe and back along the crest, which is the
    pair the step describes.
    """
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM08), None, select=0)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem08_studio_profile.png")


def lem08_dloads():
    """The distributed-loads editor holding LEM-8's crest surcharge — two points at
    240 psf, and the Direction the level crest makes irrelevant."""
    from studio.editors import DloadsEditor

    dlg = DloadsEditor().build(_load(LEM08), None)
    dlg.resize(1000, 620)
    return _grab(dlg, "lem08_studio_dloads.png")


def lem08_circles():
    """The circles editor holding LEM-8's two starting circles, typed rather than
    generated: one tangent to the toe elevation, one to the bottom of the model."""
    from studio.editors import CirclesEditor

    return _grab(CirclesEditor().build(_load(LEM08), None),
                 "lem08_studio_circles.png")


def lem08_run_lem():
    return _run_lem(LEM08, "lem08_studio_run_lem.png")


def lem09_materials():
    """The materials editor, table view, on LEM-9's two soils — both unit-weight
    columns filled at the same value, and no water in either row."""
    from studio.editors import MaterialsEditor

    return _grab(_mat_table(_lem_only(MaterialsEditor().build(_load(LEM09), None))),
                 "lem09_studio_materials.png")


def lem09_profile():
    """The profile-lines editor on LEM-9's two lines, the lower layer selected.

    Line 2 carries the wall face itself — the vertical run at x = 0 — so it is the
    one selected, and the preview draws the face against the ground surface behind
    it, which is the shape a reader most needs to check on this model.
    """
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(LEM09), None, select=1)
    dlg.resize(1180, 720)
    return _grab(dlg, "lem09_studio_profile.png")


def lem09_noncirc():
    """The non-circular editor holding LEM-9's failure surface as entered — the
    wedge running back from the wall toe, one vertex per row with its Movement."""
    from studio.editors import NonCircEditor

    dlg = NonCircEditor().build(_load(LEM09), None)
    dlg.resize(1100, 640)
    return _grab(dlg, "lem09_studio_noncirc.png")


def lem09_run_lem():
    """Run LEM on the wall, filled as LEM-9's step dictates — Janbu (Corrected)
    rather than Spencer, on a model whose Surface row is the fixed non-circular
    label."""
    return _run_lem(LEM09, "lem09_studio_run_lem.png", method="janbu")


SHOTS.update({
    "lem01_global": lem01_global,
    "lem02_run_lem": lem02_run_lem,
    "lem03_profile": lem03_profile,
    "lem03_run_lem": lem03_run_lem,
    "lem04_materials": lem04_materials,
    "lem04_profile": lem04_profile,
    "lem04_circles": lem04_circles,
    "lem04_run_lem": lem04_run_lem,
    "lem05_materials": lem05_materials,
    "lem05_profile": lem05_profile,
    "lem05_piezo": lem05_piezo,
    "lem06_materials": lem06_materials,
    "lem06_circles": lem06_circles,
    "lem06_run_lem": lem06_run_lem,
    "lem08_materials": lem08_materials,
    "lem08_profile": lem08_profile,
    "lem08_dloads": lem08_dloads,
    "lem08_circles": lem08_circles,
    "lem08_run_lem": lem08_run_lem,
    "lem09_materials": lem09_materials,
    "lem09_profile": lem09_profile,
    "lem09_noncirc": lem09_noncirc,
    "lem09_run_lem": lem09_run_lem,
})


# --------------------------------------------------------------------------- #
# LEM-7 — Strength Options Beyond Mohr-Coulomb (open-and-run, one edit per part)
#
# Every materials shot is the LIST view: the subject of this page is the strength
# option, and only the list view puts the option's own fields and the strength plot
# that confirms them beside each other. Each part is photographed twice — the row
# as the file carries it, and the row after the page's one edit — so the reader
# checks both states against a picture.
# --------------------------------------------------------------------------- #
LEM07_A = os.path.join(REPO_ROOT, "docs/lem/files/xslope_baker_clay.xlsx")
LEM07_B = os.path.join(REPO_ROOT, "docs/lem/files/xslope_low_clay.xlsx")


def _lem07_material(model, name, row=0, edit=None, width=1180, height=760):
    """The materials editor in list view, on one material of ``model``.

    ``edit`` is applied to that material's row before the editor is built, so the
    shot carries the state the tutorial's step ends in rather than a pre-edit row
    dressed up in prose. The LEM-only toggles are the owner's rule for an LEM
    page's editor captures.
    """
    from studio.editors import MaterialsEditor

    d = _load(model)
    if edit:
        mats = [dict(m) for m in d["materials"]]
        mats[row].update(edit)
        d = dict(d, materials=mats)
    dlg = _lem_only(MaterialsEditor().build(d, None))
    dlg.set_view_mode("list")
    dlg._list_view.list.setCurrentRow(row)
    dlg.resize(width, height)
    return _grab(dlg, name)


def lem07_materials_pow():
    """Part A's material as the file carries it: option `pow`, with the power-curve
    coefficients and the strength plot that draws the curve they define."""
    return _lem07_material(LEM07_A, "lem07_studio_materials_pow.png")


def lem07_materials_mc():
    """The same material after Part A's one edit — option `mc` with Baker's fitted
    c′ = 11.64 kPa and φ′ = 24.7°. The power-curve fields are gone from the form
    and the plot is a straight line."""
    return _lem07_material(LEM07_A, "lem07_studio_materials_mc.png",
                           edit={"option": "mc", "c": 11.64, "phi": 24.7})


def lem07_run_lem():
    """Part A's run: Spencer, auto search, the dialog's own 40 slices — the state
    the first run is made in, on the file exactly as downloaded."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM07_A))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem07_studio_run_lem.png")


def lem07_low_materials_cp():
    """Part B's lower layer as the file carries it: option `cp`, with the c/p rate
    and the reference elevation the profile is measured from."""
    return _lem07_material(LEM07_B, "lem07_studio_low_materials_cp.png", row=2)


def lem07_low_materials_const():
    """The same layer after Part B's edit — option `mc` at a constant 22.5 kN/m²,
    φ = 0. The c/p and r-elev fields leave the form with the option that read
    them."""
    return _lem07_material(LEM07_B, "lem07_studio_low_materials_const.png", row=2,
                           edit={"option": "mc", "c": 22.5})


def lem07_low_run_lem():
    """Part B's run: Bishop's Simplified, auto search, 50 slices — the slice count
    the step raises from the dialog's default."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "bishop", "analysis": "auto_search",
                                 "num_slices": 50},
                       slope_data=_load(LEM07_B))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem07_studio_low_run_lem.png")


SHOTS.update({
    "lem07_materials_pow": lem07_materials_pow,
    "lem07_materials_mc": lem07_materials_mc,
    "lem07_run_lem": lem07_run_lem,
    "lem07_low_materials_cp": lem07_low_materials_cp,
    "lem07_low_materials_const": lem07_low_materials_const,
    "lem07_low_run_lem": lem07_low_run_lem,
})


# --------------------------------------------------------------------------- #
# LEM-11 — Reliability (open-and-run, one edit)
#
# The only tutorial whose subject is a column set that is HIDDEN by default: the
# mat sheet's standard deviations appear when the materials editor's Reliability
# toggle is ticked. Every materials shot here therefore pins LEM *and* Reliability
# rather than reusing _lem_only (which unticks Reliability) or _all_usages (which
# unticks it too) — a page whose sentence is "tick Reliability and the σ columns
# appear" cannot be illustrated by an editor with the box clear.
# --------------------------------------------------------------------------- #
LEM11 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_reliability.xlsx")


def _lem_and_rel(dlg):
    """Tick LEM and Reliability, untick FEM and Seepage.

    Pinned rather than defaulted for the reason ``_all_usages`` documents: the
    toggles persist per editor, and Reliability's coded default is OFF, so this
    page's subject would vanish from the shot on any machine that had not just
    ticked it by hand.
    """
    for tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(tag in ("lem", "rel"))
    return dlg


def _lem11_materials(name, sigma_c=None):
    """The materials editor in list view with the Reliability toggle on.

    List view rather than the table: the σ of a parameter is a property OF that
    parameter, and only the list view puts the two in one field — ``γ`` and its
    ``± σ`` box side by side — where the table separates them by twenty columns.
    ``sigma_c`` applies the page's edit before the editor is built, so the second
    shot carries the state its step ends in.
    """
    from studio.editors import MaterialsEditor

    d = _load(LEM11)
    if sigma_c is not None:
        mats = [dict(m) for m in d["materials"]]
        mats[0]["sigma_c"] = sigma_c
        d = dict(d, materials=mats)
    dlg = _lem_and_rel(MaterialsEditor().build(d, None))
    dlg.set_view_mode("list")
    dlg._list_view.list.setCurrentRow(0)
    dlg.resize(1180, 620)
    return _grab(dlg, name)


def lem11_materials():
    """The clay as the file carries it — γ = 120 ± 8, c = 400 ± 100 — with the
    Reliability toggle ticked so both σ boxes are on the form."""
    return _lem11_materials("lem11_studio_materials.png")


def lem11_materials_edit():
    """The clay after the page's one edit: s(c) halved to 50 psf, everything else
    — γ, its σ, c itself — untouched."""
    return _lem11_materials("lem11_studio_materials_edit.png", sigma_c=50.0)


def lem11_run_lem():
    """The deterministic run the page makes first: Spencer, auto search, the
    dialog's own 40 slices, on the file exactly as downloaded."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM11))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem11_studio_run_lem.png")


def _lem11_reliability_dialog(engine, name):
    from studio.dialogs import ReliabilityDialog

    dlg = ReliabilityDialog(defaults={"engine": engine, "method": "spencer",
                                      "num_slices": 40, "search": True},
                            slope_data=_load(LEM11), app_mode="lem")
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, name)


def lem11_reliability_taylor():
    """The Reliability dialog on the Taylor series, Spencer — including the
    read-only σ summary, which is where the run states what it will vary."""
    return _lem11_reliability_dialog("taylor",
                                     "lem11_studio_reliability_taylor.png")


def lem11_reliability_mc():
    """The same dialog switched to Monte Carlo: the sample count, seed and
    distribution come live, which is the whole difference in the inputs."""
    return _lem11_reliability_dialog("mc", "lem11_studio_reliability_mc.png")


def lem11_reliability_rs():
    """The dialog on the response surface: seed, distribution and sampling stay
    live (they describe the draws); the sample count and the convergence stop
    gray out, because a surrogate's realizations cost arithmetic."""
    return _lem11_reliability_dialog("rs", "lem11_studio_reliability_rs.png")


def lem11_parametric_variance():
    """The Parametric dialog set to the variance Pareto — the plot type that is
    offered only because this model carries standard deviations.

    The sweep table is left empty, which is the state the step's two changes
    (Method, Plot type) actually leave behind: the Pareto reads every σ-carrying
    material and ignores the table, as the dialog's own note under it says.
    """
    from studio.dialogs import SensitivityDialog

    dlg = SensitivityDialog(slope_data=_load(LEM11), app_mode="lem",
                            defaults={"method": "spencer", "plot_type": "variance"})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem11_studio_parametric_variance.png")


SHOTS.update({
    "lem11_materials": lem11_materials,
    "lem11_materials_edit": lem11_materials_edit,
    "lem11_run_lem": lem11_run_lem,
    "lem11_reliability_taylor": lem11_reliability_taylor,
    "lem11_reliability_mc": lem11_reliability_mc,
    "lem11_reliability_rs": lem11_reliability_rs,
    "lem11_parametric_variance": lem11_parametric_variance,
})


# --------------------------------------------------------------------------- #
# LEM-12 — Piles (open-and-explore)
#
# The subject is a column that is EMPTY in the shipped file: H is blank on both
# pile rows, which is what puts the run on the Ito & Matsui path. Every piles shot
# here pins LEM only, so the FEM tail (E, I, Area, Head, Tip) is out of the way and
# the blank H sits among the columns that feed it — D, S, Vcap and Mcap.
# --------------------------------------------------------------------------- #
LEM12 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_piles.xlsx")
#: The page's two edits, each applied before its editor is built so the shot
#: carries the state its own step ends in.
LEM12_H = (2540.7, 1827.0)
LEM12_S = 12.0


def _lem12_model(H=None, S=None):
    d = _load(LEM12)
    rows = []
    for i, p in enumerate(d["pile_lines"]):
        row = dict(p)
        if H is not None:
            row["H"] = H[i]
        if S is not None:
            row["S"] = S
        rows.append(row)
    return dict(d, pile_lines=rows)


def _lem12_piles_table(name, **edit):
    """The piles editor's table view on both rows, LEM columns only."""
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_lem12_model(**edit), None)
    _lem_only(dlg)
    return _grab(_line_table(dlg, through="M_cap"), name)


def lem12_piles_table():
    """The two pile rows as the file carries them: H empty on both, D = 2 and
    S = 6 beside it, and the two structural capacities at the end of the LEM
    columns. The empty H is the shot's subject."""
    return _lem12_piles_table("lem12_studio_piles_table.png")


def lem12_piles_h():
    """The same table after the page's specified-force edit — one number typed
    into each row's H — with every other cell as it was."""
    return _lem12_piles_table("lem12_studio_piles_h.png", H=LEM12_H)


def lem12_piles_spacing():
    """The same table after the spacing edit: S = 12 on both rows, H still
    blank, so the force is recomputed at the wider spacing."""
    return _lem12_piles_table("lem12_studio_piles_spacing.png", S=LEM12_S)


def lem12_piles():
    """The upper pile in the editor's list view, where the Capacity / design
    group puts the blank H above the D and S it is computed from and the two
    capacities that bound it."""
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(LEM12), None)
    return _grab(_list_view(_lem_only(dlg), 1400), "lem12_studio_piles.png")


def lem12_run_lem():
    """The first run: Spencer, auto search, 40 slices, on the file as
    downloaded."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM12))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem12_studio_run_lem.png")


def lem12_run_lem_grid():
    """The same dialog with Grid search ticked — the state the bypass step
    leaves behind, everything else unchanged."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40, "grid_seed": True},
                       slope_data=_load(LEM12))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem12_studio_run_lem_grid.png")


SHOTS.update({
    "lem12_piles_table": lem12_piles_table,
    "lem12_piles": lem12_piles,
    "lem12_piles_h": lem12_piles_h,
    "lem12_piles_spacing": lem12_piles_spacing,
    "lem12_run_lem": lem12_run_lem,
    "lem12_run_lem_grid": lem12_run_lem_grid,
})


# --------------------------------------------------------------------------- #
# SEEP-1 — Seepage Under a Sheetpile (built from scratch)
#
# The first seepage tutorial, so its captures follow the build from an empty
# project: the global parameters that declare the unit system and the time unit
# every conductivity and discharge is quoted in, the one soil, the profile line
# whose notch IS the sheetpile, and the boundary conditions — the step the page
# says is the one that decides the answer.
#
# Every editor shot pins the Seepage toggle and unticks the rest, the mirror of
# the LEM pages' _lem_only: a seepage page's material row has no business showing
# a cohesion, and the k columns it does teach are hidden unless Seepage is on.
#
# Then the two dialogs that have no LEM counterpart at all — Build Mesh (twice:
# the uniform mesh, and the same dialog with the tip refinement turned on) and
# Run Seepage — and the Parametric dialog in its seepage mode, whose output
# quantity is the discharge rather than a factor of safety.
# --------------------------------------------------------------------------- #
SEEP01 = os.path.join(REPO_ROOT, "docs/seep/files/xslope_clay_blanket.xlsx")
#: The mesh the page builds first, in the Build mesh dialog's own controls:
#: tri3 (the seepage default), auto-sized at 100 divisions across the 50 m
#: section, which is the 0.5 m target size the figures are computed at.
SEEP01_MESH = {"element_type": "tri3", "auto_size": True, "size_divisions": 100,
               "target_size": 0.5}
#: The refinement the page turns on afterwards, everything else unchanged.
SEEP01_REFINE = 4.0


def _seep_only(dlg):
    """Untick every usage toggle but Seepage, so the capture shows the seepage
    problem.

    The counterpart of ``_lem_only``, and pinned for the same reason: the toggles
    are remembered per editor in QSettings, so a materials shot that simply opened
    the editor would carry whatever the last person to click had ticked — and on
    this page the k columns being visible is the whole point of the shot.
    """
    for tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(tag == "seep")
    return dlg


def seep01_global():
    """Global parameters as the seepage build leaves them: SI, and a declared
    Time of days.

    The time unit is why this shot exists on a page whose LEM siblings skip it —
    conductivity, flux and discharge all carry a time dimension, and the declared
    unit is what puts `m/day` on the input forms and `m³/day per m` on the flow
    net's title.
    """
    from studio.editors import GlobalEditor

    dlg = GlobalEditor().build(_load(SEEP01), None)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep01_studio_global.png")


def seep01_materials():
    """The one soil in table view, with Seepage ticked and nothing else.

    The table rather than the list view, which is what the LEM pages photograph:
    with only the Seepage toggle on, the table IS the seepage row — name, the two
    principal conductivities, the rotation, the unsaturated model and its
    parameters, in the mat sheet's own order — while the list view puts those
    fields in a group below three strength groups a confined seepage problem never
    fills.

    Widened through the last seepage column rather than through ``unsat``, because
    the page's claim is that this view shows the seepage band AND NOTHING ELSE: a
    shot cut off inside the band shows a column with no header, which is neither
    the band nor a table.
    """
    from studio.editors import MaterialsEditor

    dlg = _seep_only(MaterialsEditor().build(_load(SEEP01), None))
    return _grab(_mat_table(dlg, through="vg_n"),
                 "seep01_studio_materials.png")


def seep01_profile():
    """The profile line, whose six vertices carry the sheetpile: the pair either
    side of x = 30 and the vertex at (30, 7) between them make the slot the mesh
    opens as a no-flow face."""
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(SEEP01), None, select=0)
    dlg.resize(1180, 720)
    return _grab(dlg, "seep01_studio_profile.png")


def seep01_seep_bc():
    """The seepage boundary conditions on Set 1, upstream boundary selected — the
    head value and the two points it is held along."""
    from studio.editors import SeepBcEditor

    dlg = SeepBcEditor().build(_load(SEEP01), None, select=(0, 0))
    dlg.resize(1080, 560)
    return _grab(dlg, "seep01_studio_seep_bc.png")


def _seep01_build_mesh(name, **over):
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults=dict(SEEP01_MESH, **over))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, name)


def seep01_build_mesh():
    """Build Mesh as the page's first run leaves it: tri3, auto-sized at 100
    divisions, no refinement. The quadrilateral style is dimmed, because a
    triangular element type has none."""
    return _seep01_build_mesh("seep01_studio_build_mesh.png")


def seep01_build_mesh_refine():
    """The same dialog with **Refine near features** ticked and the factor at 4 —
    the one change the tip-refinement step makes, with the element type and the
    target size where they were."""
    return _seep01_build_mesh("seep01_studio_build_mesh_refine.png",
                              refine_near_features=True,
                              refine_factor=SEEP01_REFINE)


def seep01_run_seep():
    """Run Seepage on the meshed model: one BC set, the solve tolerance, and the
    model checks beside them. This file defines no transient sheet and no second
    BC set, so neither the Run type selector nor the extra BC choices appear."""
    from studio.dialogs import RunSeepDialog

    data = _load(SEEP01)
    dlg = RunSeepDialog(defaults={"tol": 1e-4}, slope_data=data,
                        has_bc2=bool((data.get("seepage_bc2") or {})
                                     .get("specified_heads")),
                        has_tseep=bool(data.get("tseep")))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep01_studio_run_seep.png")


def seep01_parametric():
    """The Parametric dialog in seepage mode, filled in as the page's conductivity
    sweep: Design mode, the soil's k₁ swept from 30 to 300 m/day in ten steps, and
    a target discharge of 100.

    The bounds are pinned rather than left at the dialog's own ±50% default,
    because the numbers the page's sweep quotes are measured at exactly these
    ones, and a shot showing 15 to 45 would be a picture of a different study.
    """
    from studio.dialogs import SensitivityDialog

    dlg = SensitivityDialog(slope_data=_load(SEEP01), app_mode="seep",
                            defaults={"mode": "design", "low": 30.0, "high": 300.0,
                                      "steps": 10, "target_fs": 100.0})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep01_studio_parametric.png")


SHOTS.update({
    "seep01_global": seep01_global,
    "seep01_materials": seep01_materials,
    "seep01_profile": seep01_profile,
    "seep01_seep_bc": seep01_seep_bc,
    "seep01_build_mesh": seep01_build_mesh,
    "seep01_build_mesh_refine": seep01_build_mesh_refine,
    "seep01_run_seep": seep01_run_seep,
    "seep01_parametric": seep01_parametric,
})


# --------------------------------------------------------------------------- #
# SEEP-2 — Unconfined Seepage Through a Zoned Dam (open and explore)
#
# The reader opens a finished file rather than building one, so these shots
# photograph what is already in it: three materials with their conductivities and
# their unsaturated model, the boundary set whose exit face is what makes the
# problem unconfined, and the two dialogs the run goes through.
#
# The materials editor is photographed twice on purpose. The table is where the
# three zones are read against each other — three conductivities three orders of
# magnitude apart, on one screen. The list view is where a single zone's
# unsaturated model is changed, and it draws the kr curve of whichever model is
# selected, which is the control this page spends its length on.
# --------------------------------------------------------------------------- #
SEEP02 = os.path.join(REPO_ROOT, "docs/seep/files/xslope_johnson_res.xlsx")
#: The mesh the page builds, in the Build mesh dialog's own controls: tri3,
#: auto-sized at 120 divisions across the 750 ft section — the 6.25 ft target the
#: sample page's discharge and its SEEP2D cross-check were both computed on.
SEEP02_MESH = {"element_type": "tri3", "auto_size": True, "size_divisions": 120,
               "target_size": 6.25}


def seep02_materials():
    """The three zones in table view with only the Seepage columns showing.

    The whole seepage material is on one screen this way: a shell at 1 ft/day, a
    core a thousand times tighter, a foundation between them, and the unsaturated
    model and its parameters carried per row rather than once for the model.
    """
    from studio.editors import MaterialsEditor

    dlg = _seep_only(MaterialsEditor().build(_load(SEEP02), None))
    return _grab(_mat_table(dlg, through="vg_n"), "seep02_studio_materials.png")


def seep02_materials_unsat():
    """The shell in list view, where the unsaturated model is chosen.

    List view rather than the table, because this is the one control the page
    changes and the form draws the selected model's own conductivity curve beside
    it — so the reader sees what picking `vg` instead of `lf` does to the curve
    before running anything.
    """
    from studio.editors import MaterialsEditor

    dlg = _seep_only(MaterialsEditor().build(_load(SEEP02), None))
    dlg._set_mode("list")
    dlg.resize(1180, 780)
    dlg.show()
    _settle()
    # Scrolled to the foot of the form rather than left at its head: the group
    # this page is about is the last one, and a shot of the Identity fields with
    # the conductivity model below the fold photographs the wrong end of it.
    from PySide6.QtWidgets import QScrollArea

    scroll = dlg._list_view.findChildren(QScrollArea)[0]
    scroll.verticalScrollBar().setValue(scroll.verticalScrollBar().maximum())
    _settle()
    return _grab(dlg, "seep02_studio_materials_unsat.png")


def seep02_seep_bc():
    """The boundary set with the **exit face** selected — the entry that makes this
    problem unconfined, and the one SEEP-1's confined model left empty."""
    from studio.editors import SeepBcEditor

    dlg = SeepBcEditor().build(_load(SEEP02), None, select=(0, 2))
    dlg.resize(1080, 560)
    return _grab(dlg, "seep02_studio_seep_bc.png")


def seep02_build_mesh():
    """Build Mesh as this page sets it: tri3 again, and 120 divisions rather than
    the dialog's own 100, which is what puts the target size at the 6.25 ft the
    sample page's discharge was computed on."""
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults=dict(SEEP02_MESH))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep02_studio_build_mesh.png")


def seep02_run_seep():
    """Run Seepage on the meshed dam. The same dialog SEEP-1 photographs, on a
    model where its **Convergence tol** is not inert: this problem iterates, and
    the page's convergence section is about what that field is compared against."""
    from studio.dialogs import RunSeepDialog

    data = _load(SEEP02)
    dlg = RunSeepDialog(defaults={"tol": 1e-4}, slope_data=data,
                        has_bc2=bool((data.get("seepage_bc2") or {})
                                     .get("specified_heads")),
                        has_tseep=bool(data.get("tseep")))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep02_studio_run_seep.png")


def seep02_profile():
    """The profile-lines editor on the foundation line, whose dive to elevation
    60 carves the cutoff key out of the foundation. The build step's featured
    lesson — three stacked lines, the selected one emphasized in the preview."""
    from studio.editors import ProfileEditor

    dlg = ProfileEditor().build(_load(SEEP02), None, select=2)
    dlg.resize(1180, 720)
    return _grab(dlg, "seep02_studio_profile.png")


SHOTS.update({
    "seep02_materials": seep02_materials,
    "seep02_materials_unsat": seep02_materials_unsat,
    "seep02_profile": seep02_profile,
    "seep02_seep_bc": seep02_seep_bc,
    "seep02_build_mesh": seep02_build_mesh,
    "seep02_run_seep": seep02_run_seep,
})


# --------------------------------------------------------------------------- #
# SEEP-3 — Transient Seepage: Reservoir Drawdown Through a Cored Earth Dam
#
# The page's build has two halves, and the captures follow them in order. The
# reader first builds an ordinary STEADY seepage model — three boundary entries on
# one set — and solves it at full pool; the BC shot photographs that set. Then the
# model is made transient, and the last three shots are the three controls that do
# it: the Transient editor where the schedule and the pool series are entered, the
# Run Seepage dialog with the Run type selector that only a file with a tseep sheet
# has, and the results view a transient run lands in, which is a canvas with a play
# bar under it rather than a single picture.
# --------------------------------------------------------------------------- #
SEEP03 = os.path.join(REPO_ROOT,
                      "docs/tutorials/files/xslope_earth_dam_drawdown.xlsx")
#: The mesh the page builds, in the Build Mesh dialog's own controls: tri3,
#: auto-sized at 64 divisions across the 110 m section — the 1.71875 m target
#: every SEEP-3 figure is computed at.
SEEP03_DIVISIONS = 64
#: The instant the play bar is parked on for its shot: the pool is a third of the
#: way through its fall and the interior surface is visibly above it, which is what
#: the page is pointing at when it tells the reader to scrub.
SEEP03_PLAYBAR_TIME = 30.0


def seep03_seep_bc():
    """The boundary set as the page's STEADY half leaves it: the upstream face held
    at full pool, the tailwater, and the exit face on the downstream slope, with the
    upstream entry selected.

    The head is shown as the number 18, not as the ``pool`` series the finished file
    binds there. That binding is the next section's work — the reader types a series
    name into this same cell after the schedule exists — so the value is set back to
    the number here, on the loaded dict only. The file on disk is untouched.
    """
    from studio.editors import SeepBcEditor

    data = _load(SEEP03)
    head = data["seepage_bc"]["specified_heads"][0]
    head["kind"] = "head"
    head["head"] = float(data["tseep"]["series"]["pool"][0])
    dlg = SeepBcEditor().build(data, None, select=(0, 0))
    dlg.resize(1080, 560)
    return _grab(dlg, "seep03_studio_seep_bc.png")


def seep03_transient_editor():
    """The Transient editor holding this model's whole schedule.

    Everything the page's transient half asks for is on this one form: the run
    duration and the save interval, the extra save times, and the series table with
    its first column renamed from the template's default ``t1`` to ``pool`` — the
    name the upstream boundary's value cell refers to. The stage-time fields are
    deliberately EMPTY: they flag the rapid-drawdown states a stability analysis
    reads, and SEEP-3 stops at the pore-pressure field.

    The preview is rendered explicitly rather than waited for: the pane's redraw is
    debounced, so a grab taken straight after ``show()`` catches a blank canvas.
    """
    from studio.editors import TransientDialog

    data = _load(SEEP03)
    dlg = TransientDialog(data.get("tseep"), data, None)
    dlg.resize(1220, 640)
    dlg.show()
    _settle()
    dlg._preview.refresh_now()
    _settle()
    dlg._preview.canvas._render_current()
    _settle()
    out = os.path.join(OUT_DIR, "seep03_studio_transient.png")
    dlg.grab().save(out)
    dlg.close()
    print("-> seep03_studio_transient.png")
    return out


def seep03_run_transient():
    """Run Seepage in **Transient** mode — the Run type selector the file's tseep
    sheet adds, and the model checks a transient run is held to.

    Transient checks are stricter than steady ones (a declared time base, storage on
    every material), so the panel beside the controls is reporting on the transient
    rules here, not the steady ones.
    """
    from studio.dialogs import RunSeepDialog

    data = _load(SEEP03)
    dlg = RunSeepDialog(defaults={"mode": "transient", "tol": 1e-4},
                        slope_data=data,
                        has_bc2=bool((data.get("seepage_bc2") or {})
                                     .get("specified_heads")),
                        has_tseep=bool(data.get("tseep")))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep03_studio_run_seep.png")


def seep03_playbar():
    """The **Seep · Transient** results tab, parked mid-drawdown.

    The run is the page's own: the tutorial's file on the page's own mesh, marched
    through ``SeepRunner`` — the same call the Run button makes — so the frame under
    the play bar is a frame the reader will have. The display options are the
    transient panel's defaults: no flow lines (a storage-release state has no flow
    net), the instantaneous water levels on, and velocity vectors, which is how a
    transient frame's flow direction is read.
    """
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    from studio.runners import SeepRunner
    from studio.transient import TransientSeepView

    data = _load(SEEP03)
    xs = [x for x, _ in data["ground_surface"].coords]
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(get_material_polygons(data),
                                        (max(xs) - min(xs)) / SEEP03_DIVISIONS,
                                        "tri3")
    data["mesh"] = mesh
    runner = SeepRunner(data, {"mode": "transient"})
    bundle, err = {}, {}
    runner.succeeded.connect(lambda b: bundle.update(b))
    runner.failed.connect(lambda m: err.setdefault("msg", m))
    with contextlib.redirect_stdout(io.StringIO()):
        runner._run_transient(data, mesh)
    if err:
        raise RuntimeError("seep03_playbar: transient run failed: %s" % err["msg"])

    opts = {"variable": "head", "levels": 12, "flowlines": False,
            "vectors": True, "phreatic": True, "show_bc_levels": True}
    view = TransientSeepView()
    # The canvas sizes its figure to the viewport, and this dam is five times as
    # wide as it is tall, so a tall view leaves the frame stranded in white space.
    view.resize(1000, 540)
    view.set_frames(bundle["seep_data"], bundle["frames"],
                    opts_getter=lambda: opts, style_getter=lambda: None,
                    keep_index=False)
    view.show()
    _settle()
    times = [float(f["time"]) for f in bundle["frames"]]
    view.set_index(min(range(len(times)),
                       key=lambda i: abs(times[i] - SEEP03_PLAYBAR_TIME)))
    _settle()
    view.canvas._render_current()          # force the raster into the scene
    _settle()
    out = os.path.join(OUT_DIR, "seep03_studio_playbar.png")
    view.grab().save(out)
    view.close()
    print("-> seep03_studio_playbar.png")
    return out


SHOTS.update({
    "seep03_seep_bc": seep03_seep_bc,
    "seep03_transient_editor": seep03_transient_editor,
    "seep03_run_transient": seep03_run_transient,
    "seep03_playbar": seep03_playbar,
})


# --------------------------------------------------------------------------- #
# SEEP-4 — Infiltration & Flux Boundaries
#
# The page's whole subject is one editor: the boundary set, with three flux
# entries in it beside the reservoir head and the toe drain. So the first shot is
# the Seep BC dialog with a flux entry open — the one on the upstream face, whose
# first point is the same (20, 10) the head boundary ends on, which is the shared
# node the page's collision rule is about. The second is the Run Seepage dialog in
# steady mode, the state the reader solves from.
#
# There is no results-view shot: the two solved fields are drawn by
# ``tools/make_tutorial_figures.py`` (seep04_dry.png / seep04_wet.png) on one
# pinned color scale so the pair compares, and a Studio grab of the same plot in
# its panel chrome adds nothing the earlier seepage tutorials have not already
# shown of that view.
# --------------------------------------------------------------------------- #
SEEP04 = os.path.join(REPO_ROOT,
                      "docs/tutorials/files/xslope_dam_infiltration.xlsx")
#: The boundary-list row the BC shot opens: 0 is the reservoir head, 1-3 the three
#: infiltration blocks, 4 the toe drain. Row 1 is the upstream-face block, which
#: starts at the waterline the head boundary ends on.
SEEP04_BC_ROW = 1


def seep04_seep_bc():
    """The boundary set complete: the reservoir head, the three infiltration
    blocks and the toe drain, with the upstream-face block open.

    Its points table reads (20, 10) → (24, 12): the block starts exactly where the
    reservoir head boundary stops, so the two share a node, and its value is the
    projected 8.94427191e-09 rather than the vertical rain rate.
    """
    from studio.editors import SeepBcEditor

    data = _load(SEEP04)
    dlg = SeepBcEditor().build(data, None, select=(0, SEEP04_BC_ROW))
    dlg.resize(1080, 560)
    return _grab(dlg, "seep04_studio_seep_bc.png")


def seep04_run_seep():
    """Run Seepage in **Steady** mode — the model checks a steady unconfined run is
    held to, on a model whose boundary set is a head, three fluxes and an exit face.
    """
    from studio.dialogs import RunSeepDialog

    data = _load(SEEP04)
    dlg = RunSeepDialog(defaults={"mode": "steady", "tol": 1e-4},
                        slope_data=data,
                        has_bc2=bool((data.get("seepage_bc2") or {})
                                     .get("specified_heads")),
                        has_tseep=bool(data.get("tseep")))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "seep04_studio_run_seep.png")


SHOTS.update({
    "seep04_seep_bc": seep04_seep_bc,
    "seep04_run_seep": seep04_run_seep,
})


# --------------------------------------------------------------------------- #
# FEM-1 — Strength Reduction Basics
# --------------------------------------------------------------------------- #
FEM01_START = os.path.join(REPO_ROOT,
                           "docs/tutorials/files/xslope_ssrm_embankment_start.xlsx")
FEM01_DONE = os.path.join(REPO_ROOT,
                          "docs/tutorials/files/xslope_ssrm_embankment.xlsx")


def _fem_only(dlg):
    """Untick every usage toggle but FEM, the mirror of ``_lem_only``.

    The FEM band is not a band of its own: unit weight, the strength option and
    c/φ are shared with the limit-equilibrium set, so ticking FEM alone still
    shows the whole material — with E and ν, which is what the page's step adds,
    at the right of the row instead of forty columns away.
    """
    toggles = getattr(dlg, "_toggles", None) or {}
    for tag, cb in toggles.items():
        cb.setChecked(tag == "fem")
    return dlg


def _fem01_meshed(path=FEM01_DONE):
    """The completed model with the tutorial's own mesh attached.

    The mesh is BUILT here rather than read from a sidecar: neither tutorial file
    ships one, because building it is a step the page teaches. Studio's Run FEM
    action is unreachable without a mesh, so the dialog is photographed in the
    state the Build Mesh step leaves behind.

    The element type and size are the COMPLETED file's, whichever model is being
    meshed: the starter declares neither (the reader types them into Build Mesh),
    and photographing its Run FEM dialog on a different discretization would put
    two meshes in one tutorial.
    """
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    done = _load(FEM01_DONE)
    data = done if path == FEM01_DONE else _load(path)
    with contextlib.redirect_stdout(io.StringIO()):
        data["mesh"] = build_mesh_from_polygons(
            get_material_polygons(data), done["target_size"], done["element_type"],
            size_regions=extract_size_regions(data))
    return data


def fem01_materials():
    """The materials editor, table view, FEM columns: the row the page's second
    step completes.

    Photographed on the COMPLETED file, so the figure is the state the reader is
    working toward — E = 2,088,500 psf and ν = 0.3 filled in beside the strength
    the limit-equilibrium run already used. The window reaches through ``nu``,
    the last column this problem fills.
    """
    from studio.editors import MaterialsEditor

    dlg = _fem_only(MaterialsEditor().build(_load(FEM01_DONE), None))
    return _grab(_mat_table(dlg, through="nu"), "fem01_studio_materials.png")


def fem01_run_lem():
    """Run LEM on the STARTER — the model with no elastic constants at all.

    The page opens on a factor of safety by slices, and this is the figure that
    says why it can: the checks column beside the controls is clean and **Run** is
    enabled on a material whose E and ν cells are empty, because neither is a
    limit-equilibrium input. The same file will not run the finite element engine
    until they are filled, which is the next shot's subject.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(FEM01_START))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem01_studio_run_lem.png")


def fem01_build_mesh():
    """Build Mesh set to the tutorial's mesh: tri6, auto-sizing off, 3.5 ft.

    Auto-size is unticked so the target size box is live — the completed file
    declares 3.5, and a declared size is what turns auto-sizing off (see
    ``MainWindow._file_defaults``), so this is the dialog the reader's own file
    opens. tri6 is the element type the dialog already defaults to; the page's
    point is that the default is quadratic and must stay that way.
    """
    from studio.dialogs import BuildMeshDialog

    data = _load(FEM01_DONE)
    dlg = BuildMeshDialog(defaults={"element_type": data["element_type"],
                                    "target_size": float(data["target_size"]),
                                    "auto_size": False})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem01_studio_build_mesh.png")


def fem01_run_fem():
    """The Run FEM dialog on the meshed model, as the reader presses **Run**.

    SSRM with the bracket the page walks — F min 1.0, F max 2.0 — the 0.01
    bisection tolerance, and the rest at the dialog's own defaults: rollers on the
    sides, no K0 initialization, non-convergence as the failure criterion, and the
    at-failure capture on with its 0.15 margin, which is what produces the
    mechanism the results figures draw.
    """
    from studio.dialogs import RunFemDialog

    data = _fem01_meshed()
    dlg = RunFemDialog(defaults={"analysis": "ssrm", "F_min": 1.0, "F_max": 2.0,
                                 "tolerance": 0.01},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem01_studio_run_fem.png")


def fem01_run_fem_no_elastic():
    """The same dialog on the STARTER, meshed but with E and ν still blank.

    **Run** is disabled and the checks column names both missing constants. This
    is the page's transition figure: the limit-equilibrium run above needed
    neither, and the finite element run refuses without them.
    """
    from studio.dialogs import RunFemDialog

    data = _fem01_meshed(FEM01_START)
    dlg = RunFemDialog(defaults={"analysis": "ssrm", "F_min": 1.0, "F_max": 2.0,
                                 "tolerance": 0.01},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem01_studio_run_fem_no_elastic.png")


SHOTS.update({
    "fem01_materials": fem01_materials,
    "fem01_run_lem": fem01_run_lem,
    "fem01_build_mesh": fem01_build_mesh,
    "fem01_run_fem": fem01_run_fem,
    "fem01_run_fem_no_elastic": fem01_run_fem_no_elastic,
})


# --------------------------------------------------------------------------- #
# FEM-2 — Reinforcement: LEM against FEM
# --------------------------------------------------------------------------- #
FEM02_START = os.path.join(REPO_ROOT,
                           "docs/tutorials/files/xslope_reinforced_slope_start.xlsx")
FEM02_DONE = os.path.join(REPO_ROOT,
                          "docs/tutorials/files/xslope_reinforced_slope.xlsx")
#: The line whose detail panel is photographed: the most heavily loaded of the six.
FEM02_DETAIL_LINE = "Line 5"


def _fem02_meshed(path=FEM02_DONE):
    """The model with the tutorial's own mesh attached — the state the Build Mesh
    step leaves behind, and the only state Studio's Run FEM action is reachable in.

    The element type and target size are the COMPLETED file's whichever model is
    meshed, and the reinforcement lines go in as constraint lines so the bars land
    on mesh edges, exactly as ``studio.runners.MeshWorker`` does it. **Refine thin
    zones is off**, which is what the page instructs: the dialog's own default is
    on, and on this section it moves the mesh from 2,101 elements to 5,096 and the
    peak-residual factor of safety from 1.5117 to 1.4414.
    """
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_size_regions, get_material_polygons)

    done = _load(FEM02_DONE)
    data = done if path == FEM02_DONE else _load(path)
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(data)
    with contextlib.redirect_stdout(io.StringIO()):
        data["mesh"] = build_mesh_from_polygons(
            get_material_polygons(data, reinf_lines=lines), done["target_size"],
            done["element_type"], lines=lines or None,
            element_size_1d=done.get("element_size_1d"),
            size_regions=extract_size_regions(data))
    return data


def fem02_materials():
    """The materials editor, table view, FEM columns, AFTER the elastic pair is
    entered — the state the page's Materials step leaves behind.

    The starter delivers both soils with E and nu blank, so this is the filled
    state the reader is working toward, photographed on the COMPLETED file the
    way FEM-1's ``fem01_studio_materials`` is. The blank state is not shot: the
    Run FEM refusal above it is what shows the reader the cells are empty.
    """
    from studio.editors import MaterialsEditor

    dlg = _fem_only(MaterialsEditor().build(_load(FEM02_DONE), None))
    return _grab(_mat_table(dlg, through="nu"), "fem02_studio_materials.png")


def fem02_reinforce_blank():
    """The reinforcement table as the STARTER delivers it: six lines with their
    geometry, tensile capacity and pullout lengths, and Tres, E and Area empty.

    Only the FEM toggle is ticked: the LEM inputs are complete, so the table
    shows just the geometry and the three columns the reader is about to fill.
    """
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(FEM02_START), None)
    for tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(tag == "fem")
    return _grab(_line_table(dlg, through="area"),
                 "fem02_studio_reinforce_blank.png")


def fem02_reinforce_filled():
    """The same table on the COMPLETED file — Tres, E and Area entered on all six
    lines. The state the reader is working toward, beside the blank one."""
    from studio.editors import ReinforcementEditor

    dlg = ReinforcementEditor().build(_load(FEM02_DONE), None)
    for tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(tag == "fem")
    return _grab(_line_table(dlg, through="area"),
                 "fem02_studio_reinforce_filled.png")


def fem02_build_mesh():
    """Build Mesh set to the tutorial's mesh: tri6, auto-sizing off, 2 ft, and
    **Refine thin zones unticked**.

    The last of those is the only box on this dialog the page asks the reader to
    change from its default. The shell is a 1.19 ft facing band, thinner than one
    element at 2 ft, so the refinement fires on it — and the mesh it produces is a
    different answer, not a better-resolved version of this one.
    """
    from studio.dialogs import BuildMeshDialog

    data = _load(FEM02_DONE)
    dlg = BuildMeshDialog(defaults={"element_type": data["element_type"],
                                    "target_size": float(data["target_size"]),
                                    "auto_size": False,
                                    "refine_thin_zones": False})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem02_studio_build_mesh.png")


def fem02_run_fem():
    """Run FEM on the meshed completed model, with the checks column beside it.

    SSRM over the dialog's own bracket, and everything else at the dialog's own
    defaults — including **Max iterations per trial**, which opens at 12,000. The
    budget is where the automatic extension starts rather than where a trial dies,
    so the reader changes nothing here.
    """
    from studio.dialogs import RunFemDialog

    data = _fem02_meshed()
    dlg = RunFemDialog(defaults={"analysis": "ssrm", "F_min": 1.0, "F_max": 2.0,
                                 "tolerance": 0.01},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem02_studio_run_fem.png")



def fem02_details():
    """The 1D Details panel on the most heavily loaded line of the peak-residual
    run — the axial force along the bar over its capacity envelope, and the bond
    transfer rate under it.

    The run is made HERE rather than read from a sidecar: neither tutorial file
    ships companions, so the panel is photographed on a solve at the page's own
    settings (tri6 at 2 ft with the thin-zone refinement off, the dialog's
    bracket, 12,000 iterations a trial). That costs about two minutes.
    """
    from studio.fem_details_dialog import FemDetailsDialog
    from PySide6.QtCore import Qt
    from xslope.fem import build_fem_data, solve_ssrm

    data = _fem02_meshed()
    fem_data = build_fem_data(data, data["mesh"])
    with contextlib.redirect_stdout(io.StringIO()):
        result = solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.01,
                            debug_level=0, failure_criterion="non_convergence",
                            max_iterations=12000)
    dlg = FemDetailsDialog(fem_data, result["last_solution"], data,
                           model_path=FEM02_DONE,
                           failure_solution=result.get("failure_solution"))
    dlg.resize(1140, 660)
    for row in range(dlg.list.count()):
        entry = dlg.list.item(row).data(Qt.UserRole)
        if entry and entry["label"] == FEM02_DETAIL_LINE:
            dlg.list.setCurrentRow(row)
            break
    return _grab(dlg, "fem02_studio_1d_details.png")


SHOTS.update({
    "fem02_materials": fem02_materials,
    "fem02_reinforce_blank": fem02_reinforce_blank,
    "fem02_reinforce_filled": fem02_reinforce_filled,
    "fem02_build_mesh": fem02_build_mesh,
    "fem02_run_fem": fem02_run_fem,
    "fem02_details": fem02_details,
})


# --------------------------------------------------------------------------- #
# FEM-3 — Piles: LEM against FEM
# --------------------------------------------------------------------------- #
#: The discrete row: one model, two sample pages.  The page links the
#: limit-equilibrium copy and enters the mesh and the bracket by hand, so every
#: dialog below is photographed on a mesh built here at the page's own settings
#: rather than on a committed companion.
FEM03_PILES = os.path.join(REPO_ROOT, "docs/lem/files/xslope_piles.xlsx")
#: The continuous wall: the completed file of the page's second half, written by
#: ``tools/build_pile_wall_tutorial.py`` on the same slope.
FEM03_WALL_DONE = os.path.join(REPO_ROOT,
                               "docs/tutorials/files/xslope_pile_wall.xlsx")
#: The mesh the page builds for both models, and the finer 1D element size its
#: optional refinement step enters.
FEM03_ELEMENT_TYPE = "tri6"
FEM03_TARGET_SIZE = 2.0
FEM03_REFINED_1D = 0.5
#: The strength-reduction settings the page runs every trial at.  The bracket
#: reaches 2.0 because the socketed runs stand at 1.6 and above; everything else
#: is the dialog as it opens.
FEM03_BRACKET = (1.0, 2.0)
FEM03_TOLERANCE = 0.01
FEM03_MAX_ITERATIONS = 12000
FEM03_CRITERION = "non_convergence"


def _fem03_meshed(path, element_size_1d=None):
    """The model with a mesh attached — the state Build Mesh leaves behind, and
    the only state Studio's Run FEM action is reachable in.

    The mesh is built at the page's own settings, tri6 at 2 ft, with the pile or
    wall lines carried in as constraint lines.
    """
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_point_constraints,
                             extract_size_regions, get_material_polygons)

    data = _load(path)
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(data)
    with contextlib.redirect_stdout(io.StringIO()):
        data["mesh"] = build_mesh_from_polygons(
            get_material_polygons(data, reinf_lines=lines),
            FEM03_TARGET_SIZE, FEM03_ELEMENT_TYPE, lines=lines or None,
            element_size_1d=element_size_1d,
            point_constraints=extract_point_constraints(data),
            size_regions=extract_size_regions(data))
    return data


def fem03_piles_table():
    """The two pile rows with BOTH usage bands shown — every column the two
    engines read, side by side.

    E is given on both rows and I and Area are blank: the beam formulation
    computes the section constants from the diameter when those cells are empty,
    and then divides EA and EI by the spacing.  So the pair of cells the finite
    element engine ultimately reads its stiffness from is D and S, which the
    editor shows in both bands, uncolored, the way it shows the reinforcement
    editor's Spacing.  Both bands are photographed rather than FEM alone because
    the subject is one row read two ways: the same D and S feed the Ito & Matsui
    force on the left and the smeared beam section on the right.

    The table is taken out to Tip, the last column, because the page's spacing
    section turns on what Head and Tip are set to.
    """
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(FEM03_PILES), None)
    for _tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(True)
    return _grab(_line_table(dlg, through="tip_fixity"),
                 "fem03_studio_piles_table.png")


def fem03_run_fem_piles():
    """Run FEM on the meshed pile model: strength reduction over the page's
    bracket, with the checks column beside it."""
    from studio.dialogs import RunFemDialog

    data = _fem03_meshed(FEM03_PILES)
    dlg = RunFemDialog(defaults={"analysis": "ssrm",
                                 "F_min": FEM03_BRACKET[0],
                                 "F_max": FEM03_BRACKET[1],
                                 "tolerance": FEM03_TOLERANCE},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem03_studio_run_fem_piles.png")


def fem03_wall_row():
    """The one row the reader adds, in the editor's table view: the wall's two
    endpoints, its axial and bending section constants entered directly, spacing
    1, a moment capacity, and D and Vcap left empty.

    Both usage bands are shown rather than the FEM band alone, because the row's
    subject is which cells a continuous member fills and which it leaves blank,
    and those cells sit in both bands.  The table reaches Tip, which is the cell
    the page's second run changes.
    """
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(FEM03_WALL_DONE), None)
    # Both bands ticked, explicitly: which band is shown is a session setting, so
    # a shot that does not set it photographs whatever the last dialog left.
    for _tag, cb in (getattr(dlg, "_toggles", None) or {}).items():
        cb.setChecked(True)
    return _grab(_line_table(dlg, through="tip_fixity"),
                 "fem03_studio_wall_row.png")


def fem03_build_mesh_1d():
    """Build Mesh with the page's optional refinement entered: tri6 at 2 ft, and
    a 1D element size of 0.5 ft.

    The 1D size box is the only one the refinement step changes.  Left blank it
    follows the mesh element size, which is what every other run on the page
    uses.  Every other control is left at the dialog's own default, Refine thin
    zones included: this section carries one material, so the thin-zone plan is
    empty and the box changes nothing.
    """
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults={"element_type": FEM03_ELEMENT_TYPE,
                                    "target_size": FEM03_TARGET_SIZE,
                                    "auto_size": False,
                                    "element_size_1d": FEM03_REFINED_1D})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem03_studio_build_mesh_1d.png")


def fem03_run_fem_wall():
    """Run FEM on the meshed wall model, at the bracket the page's numbers are
    measured over."""
    from studio.dialogs import RunFemDialog

    data = _fem03_meshed(FEM03_WALL_DONE)
    dlg = RunFemDialog(defaults={"analysis": "ssrm",
                                 "F_min": FEM03_BRACKET[0],
                                 "F_max": FEM03_BRACKET[1],
                                 "tolerance": FEM03_TOLERANCE},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem03_studio_run_fem_wall.png")


def fem03_wall_details():
    """The 1D Details panel on the wall with its tip fixed: lateral displacement,
    shear, bending moment and soil reaction down its length.

    The moment capacity the row declares draws a dashed line on the moment panel;
    no Ito & Matsui envelope is drawn, because the row declares no diameter — which
    is correct for a member with no gaps for soil to arch across, and is the point
    the page reads off this shot.

    The run is made HERE rather than read from a sidecar: the tutorial files ship
    no companions, so the panel is photographed on a solve at the page's own
    settings, which takes about a minute.
    """
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope.fem import build_fem_data, solve_ssrm

    data = _fem03_meshed(FEM03_WALL_DONE)
    data["pile_lines"][0]["tip_fixity"] = "fixed"
    with contextlib.redirect_stdout(io.StringIO()):
        fem_data = build_fem_data(data, data["mesh"])
        result = solve_ssrm(fem_data, F_min=FEM03_BRACKET[0],
                            F_max=FEM03_BRACKET[1], tolerance=FEM03_TOLERANCE,
                            debug_level=0, failure_criterion=FEM03_CRITERION,
                            max_iterations=FEM03_MAX_ITERATIONS)
    dlg = FemDetailsDialog(fem_data, result["last_solution"], data,
                           model_path=FEM03_WALL_DONE,
                           failure_solution=result.get("failure_solution"))
    dlg.resize(1140, 660)
    dlg.list.setCurrentRow(dlg.list.count() - 1)
    return _grab(dlg, "fem03_studio_wall_1d_details.png")


SHOTS.update({
    "fem03_piles_table": fem03_piles_table,
    "fem03_run_fem_piles": fem03_run_fem_piles,
    "fem03_wall_row": fem03_wall_row,
    "fem03_build_mesh_1d": fem03_build_mesh_1d,
    "fem03_run_fem_wall": fem03_run_fem_wall,
    "fem03_wall_details": fem03_wall_details,
})


# --------------------------------------------------------------------------- #
# COMBO-1 — Seepage into Stability
#
# The reader opens a finished file and runs it three ways, so these shots
# photograph the three run dialogs and the one table that connects them: the
# materials sheet's pore-pressure column, which is what makes a seepage solution
# reach a stability run at all.
#
# Every dialog is photographed on a model that has been meshed AND solved for
# seepage, because that is the state the reader is in when each one is opened —
# and because the checks column is the subject of half of them: a material set to
# u = seep with no solved field is an error that blocks the run.
# --------------------------------------------------------------------------- #
COMBO01 = os.path.join(REPO_ROOT, "docs/tutorials/files/xslope_johnson_res.xlsx")
#: Build Mesh at the dialog's own defaults — quadratic triangles, auto-sizing on,
#: 100 divisions across the 750 ft section.
COMBO01_MESH = {"element_type": "tri6", "auto_size": True, "size_divisions": 100,
                "target_size": 7.5}
COMBO01_METHOD = "spencer"
COMBO01_SLICES = 40
COMBO01_BRACKET = (1.0, 2.0)
COMBO01_TOLERANCE = 0.01

_combo01_cache = {}


def _combo01_solved():
    """The workbook with the tutorial's mesh built on it and the seepage solution
    attached, exactly as the reader's session carries it after the seepage run.

    Built here rather than read from a companion file: the tutorial copy of the
    workbook ships no sidecars, because building the mesh and solving the seepage
    are the first two steps the page teaches. Cached, because three dialogs are
    photographed on the same state and the solve is not free.
    """
    if "data" not in _combo01_cache:
        from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                                 get_material_polygons)
        from xslope.seep import (apply_steady_stability_field, build_seep_data,
                                 run_seepage_analysis)

        data = _load(COMBO01)
        with contextlib.redirect_stdout(io.StringIO()):
            data["mesh"] = build_mesh_from_polygons(
                get_material_polygons(data), COMBO01_MESH["target_size"],
                COMBO01_MESH["element_type"],
                size_regions=extract_size_regions(data))
            seep_data = build_seep_data(data["mesh"], data)
            solution = run_seepage_analysis(seep_data, tol=1e-4, max_iter=400)
            apply_steady_stability_field(data, solution, bc=1)
        _combo01_cache["data"] = data
    return _combo01_cache["data"]


def combo01_build_mesh():
    """Build Mesh with nothing changed: quadratic triangles, auto-sizing on, 100
    divisions.

    The element type is the shot's subject. It is the dialog's own default, and on
    this file it is also a requirement rather than a preference — the same mesh
    goes on to carry the strength reduction, which linear elements cannot run
    honestly.
    """
    from studio.dialogs import BuildMeshDialog

    dlg = BuildMeshDialog(defaults=dict(COMBO01_MESH))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo01_studio_build_mesh.png")


def combo01_run_seep():
    """Run Seepage on the meshed dam, at the dialog's own two defaults."""
    from studio.dialogs import RunSeepDialog

    data = _load(COMBO01)
    dlg = RunSeepDialog(defaults={"tol": 1e-4, "max_iter": 400}, slope_data=data,
                        has_bc2=bool((data.get("seepage_bc2") or {})
                                     .get("specified_heads")),
                        has_tseep=bool(data.get("tseep")))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo01_studio_run_seep.png")


def combo01_materials():
    """The three zones in table view with the LEM columns showing — the shot the
    page's hinge is read off.

    The **u** column is the subject: all three rows read `seep`, which is what
    sends the solved seepage field to every slice base and every Gauss point.
    The LEM band rather than the seepage band, because u is a stability input —
    it says what the stability engines do with a field the seepage run produced,
    and a row left on `none` runs the same seepage solve and then ignores it.
    """
    from studio.editors import MaterialsEditor

    dlg = _lem_only(MaterialsEditor().build(_load(COMBO01), None))
    return _grab(_mat_table(dlg, through="u"), "combo01_studio_materials.png")


def combo01_run_lem():
    """Run LEM on the solved model, with Method set to Spencer — the one field
    the page changes from the dialog's defaults.

    Photographed with the seepage solution attached, so the checks column is
    clean. The same dialog on the same file before the seepage run carries an
    error instead: three materials read u = seep and there is no field for them
    to read.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": COMBO01_METHOD,
                                 "analysis": "auto_search",
                                 "num_slices": COMBO01_SLICES},
                       slope_data=_combo01_solved())
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo01_studio_run_lem.png")


def combo01_run_fem():
    """Run FEM on the same solved model, at the dialog's own defaults: strength
    reduction over [1.0, 2.0], a 0.01 bisection tolerance, 12,000 iterations a
    trial, rollers on the sides and non-convergence as the failure criterion.

    Nothing here names the seepage solution, which is the point: the pore
    pressures reach the Gauss points through the materials' own u column, and the
    run dialog has no water control of its own.
    """
    from studio.dialogs import RunFemDialog

    data = _combo01_solved()
    dlg = RunFemDialog(defaults={"analysis": "ssrm", "F_min": COMBO01_BRACKET[0],
                                 "F_max": COMBO01_BRACKET[1],
                                 "tolerance": COMBO01_TOLERANCE},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo01_studio_run_fem.png")


SHOTS.update({
    "combo01_build_mesh": combo01_build_mesh,
    "combo01_run_seep": combo01_run_seep,
    "combo01_materials": combo01_materials,
    "combo01_run_lem": combo01_run_lem,
    "combo01_run_fem": combo01_run_fem,
})


# --------------------------------------------------------------------------- #
# COMBO-2 — Rapid Drawdown
#
# Two workbooks, photographed at the point in the page where each dialog is
# opened. The starter's shots carry the piezometric pair and u = piezo; the
# completed model's carry the second boundary set, the schedule and the Run LEM
# dialog with the drawdown checkbox on, where the two stage-time fields appear.
# --------------------------------------------------------------------------- #
COMBO02_START = os.path.join(REPO_ROOT,
                             "docs/tutorials/files/xslope_johnson_rapid_start.xlsx")
COMBO02 = os.path.join(REPO_ROOT,
                       "docs/tutorials/files/xslope_johnson_rapid.xlsx")
#: The mesh the page builds: linear triangles, auto-sizing on, 100 divisions
#: across the 750 ft section.
COMBO02_MESH = {"element_type": "tri3", "auto_size": True, "size_divisions": 100,
                "target_size": 7.5}
COMBO02_METHOD = "spencer"
COMBO02_SLICES = 40

_combo02_cache = {}


def _combo02_solved():
    """The completed workbook with the page's mesh built and BOTH steady fields
    attached — the state the reader is in when Run LEM is opened for the
    two-steady run, and the state the drawdown checks are reported against.

    Cached: two solves and a mesh for one dialog is not free, and the shipped file
    carries no sidecars because solving is what the page teaches.
    """
    if "data" not in _combo02_cache:
        from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                                 get_material_polygons)
        from xslope.seep import (apply_steady_stability_field, build_seep_data,
                                 run_seepage_analysis)

        data = _load(COMBO02)
        with contextlib.redirect_stdout(io.StringIO()):
            data["mesh"] = build_mesh_from_polygons(
                get_material_polygons(data), COMBO02_MESH["target_size"],
                COMBO02_MESH["element_type"],
                size_regions=extract_size_regions(data))
            for bc in (1, 2):
                seep_data = build_seep_data(data["mesh"], data, seep_bc=bc)
                solution = run_seepage_analysis(seep_data, tol=1e-4, max_iter=400)
                apply_steady_stability_field(data, solution, bc=bc)
        _combo02_cache["data"] = data
    return _combo02_cache["data"]


#: Widest a preview-bearing editor is photographed at. A docs figure is rendered
#: at 1000 px, so a dialog wider than this is downscaled past the point where its
#: field labels stay legible.
_PREVIEW_MAX_WIDTH = 1400


def _size_preview(dlg, slope_data, width=_PREVIEW_MAX_WIDTH):
    """Size a preview-bearing editor so its canvas carries the section's own aspect.

    The preview figure is drawn to the canvas viewport, so a tall pane on a section
    four times as wide as it is high strands the drawing in white space with the
    padding all in y. Both dimensions are measured off the built dialog rather than
    guessed: the height is the form's own ``sizeHint`` — the table is then never
    clipped and the preview is no taller than it has to be — and the width is
    solved for from two measurements of how much of it the layout gives the
    canvas, so the viewport comes out as close to the section's width-to-height
    ratio as the cap allows.
    """
    dlg.show()
    _settle()
    dlg.resize(dlg.width(), dlg.sizeHint().height())
    _settle()
    samples = []
    for w in (900, width):
        dlg.resize(w, dlg.height())
        _settle()
        samples.append((w, dlg._preview.canvas.view.viewport().width()))
    (w0, v0), (w1, v1) = samples
    xs = [x for x, _ in slope_data["ground_surface"].coords]
    ys = [y for _, y in slope_data["ground_surface"].coords]
    aspect = (max(xs) - min(xs)) / max(1e-9, max(ys) - min(ys))
    want = aspect * dlg._preview.canvas.view.viewport().height()
    slope = (v1 - v0) / float(w1 - w0) if w1 != w0 else 0.0
    target = w1 if slope <= 0 else w0 + (want - v0) / slope
    dlg.resize(int(min(max(target, w0), width)), dlg.height())
    _settle()
    dlg._preview.refresh_now()
    _settle()
    dlg._preview.canvas._render_current()
    _settle()
    return dlg


def combo02_materials():
    """The three zones in table view with the LEM columns showing.

    The **d** and **psi** columns are the subject: only the core carries the pair,
    which is what makes it the one material the drawdown procedure treats as
    undrained. The LEM band, because both columns are read by the limit
    equilibrium engine alone.
    """
    from studio.editors import MaterialsEditor

    dlg = _lem_only(MaterialsEditor().build(_load(COMBO02_START), None))
    return _grab(_mat_table(dlg, through="u"), "combo02_studio_materials.png")


def combo02_piezo():
    """The piezometric editor with **Line 2** active.

    Both lines are on the preview — the active one bold with its vertices, the
    other dimmed — so the pair reads as what it is: one surface at full pool and
    one at the end of the drawdown, six points and five, meeting again downstream
    of the core where the drawdown has not reached.
    """
    from studio.editors import PiezoEditor

    data = _load(COMBO02_START)
    dlg = PiezoEditor().build(data, None)
    dlg._tabs.setCurrentIndex(1)
    _size_preview(dlg, data)
    return _grab(dlg, "combo02_studio_piezo.png", settle=False)


def combo02_seep_bc2():
    """The boundary editor on its **Set 2 (rapid drawdown)** tab, with the
    upstream entry selected.

    Set 2 is the drawn-down steady problem: two heads at the elevation-100
    tailwater datum and nothing else. It takes plain heads only — a reservoir
    boundary and a series-bound value both belong on Set 1 — which is what the
    help line above the tabs says.
    """
    from studio.editors import SeepBcEditor

    dlg = SeepBcEditor().build(_load(COMBO02), None, select=(1, 0))
    dlg.resize(1080, 560)
    return _grab(dlg, "combo02_studio_seep_bc2.png")


def combo02_transient():
    """The Transient editor with the pool schedule and BOTH stage times filled.

    The stage fields are the difference between this shot and SEEP-3's, where they
    are deliberately blank: they name the two saved frames the drawdown analysis
    reads, stage 1 at t = 0 with the reservoir still full and stage 2 at t = 50
    with it down. Both are forced into the saved-frame schedule, so each is a
    computed frame.

    The preview is rendered explicitly rather than waited for: the pane's redraw
    is debounced, so a grab taken straight after ``show()`` catches a blank canvas.
    """
    from studio.editors import TransientDialog

    data = _load(COMBO02)
    dlg = TransientDialog(data.get("tseep"), data, None)
    dlg.resize(1220, 640)
    dlg.show()
    _settle()
    dlg._preview.refresh_now()
    _settle()
    dlg._preview.canvas._render_current()
    _settle()
    out = os.path.join(OUT_DIR, "combo02_studio_run_transient.png")
    dlg.grab().save(out)
    dlg.close()
    print("-> combo02_studio_run_transient.png")
    return out


def _combo02_march():
    """The solved transient march, as the Run LEM dialog receives it from the
    results document once Part 3's seepage run has finished.

    The dialog needs it for two things the drawdown page photographs: the stage
    fields offer instants out of a solved march, and the model checks read the two
    stage times as the statement that this run takes both of its states from the
    march rather than from boundary set 2. Cached beside the mesh, since the march
    is the long run of the page.
    """
    if "march" not in _combo02_cache:
        from xslope.seep import (build_seep_data, build_tseep_data,
                                 run_seepage_analysis, run_transient_seepage)

        data = _combo02_solved()
        with contextlib.redirect_stdout(io.StringIO()):
            seep_data = build_seep_data(data["mesh"], data, seep_bc=1)
            run_seepage_analysis(seep_data, tol=1e-4, max_iter=400)
            _combo02_cache["march"] = run_transient_seepage(
                seep_data, build_tseep_data(data))
    return _combo02_cache["march"]


def _combo02_cleared_bc2(slope_data):
    """``slope_data`` with boundary set 2 emptied, as Part 3 has the reader clear it."""
    return dict(slope_data,
                seepage_bc2={"specified_heads": [], "specified_fluxes": [],
                             "exit_face": []})


def _combo02_run_lem(slope_data, name, transient=None):
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": COMBO02_METHOD,
                                 "analysis": "auto_search",
                                 "num_slices": COMBO02_SLICES, "rapid": True},
                       slope_data=slope_data, transient=transient)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, name)


def combo02_run_lem():
    """Run LEM with **Rapid drawdown** ticked, on the starter file.

    Two fields differ from the dialog's defaults: Method (Spencer) and the
    drawdown checkbox itself, which runs the three-stage procedure and puts the
    run through the drawdown checks as well as the ordinary limit equilibrium
    ones. Analysis stays on Auto search, so the run finds its own critical circle
    from the starting circle the file carries. The starter carries no transient
    schedule, so no stage-time fields are on the form yet.
    """
    return _combo02_run_lem(_load(COMBO02_START), "combo02_studio_run_lem.png")


def combo02_run_lem_staged():
    """The same dialog on the completed model, which carries the schedule.

    The **Rapid-drawdown stage times** group is what the schedule adds: two
    instants rather than the single seepage time an ordinary run on a transient
    model selects, holding the tseep sheet's own stage_1 and stage_2.

    The model is Part 3's: the march solved, and boundary set 2 cleared, which is
    what leaves the checks column carrying the page's one standing warning.
    """
    return _combo02_run_lem(_combo02_cleared_bc2(_combo02_solved()),
                            "combo02_studio_run_lem_staged.png",
                            transient=_combo02_march())


def combo02_run_lem_staged_bc2():
    """The same dialog with boundary set 2 still on the file.

    Two stage times and a solved march say the drawn-down state comes from the
    march, so the second boundary set states an analysis this run neither solves
    nor reads, and the checks say so rather than letting it look consulted.
    """
    return _combo02_run_lem(_combo02_solved(),
                            "combo02_studio_run_lem_staged_bc2.png",
                            transient=_combo02_march())


SHOTS.update({
    "combo02_materials": combo02_materials,
    "combo02_piezo": combo02_piezo,
    "combo02_seep_bc2": combo02_seep_bc2,
    "combo02_transient": combo02_transient,
    "combo02_run_lem": combo02_run_lem,
    "combo02_run_lem_staged": combo02_run_lem_staged,
    "combo02_run_lem_staged_bc2": combo02_run_lem_staged_bc2,
})


# --------------------------------------------------------------------------- #
#  COMBO-3 — factor of safety versus time
# --------------------------------------------------------------------------- #
COMBO03 = os.path.join(REPO_ROOT,
                       "docs/tutorials/files/xslope_earth_dam_fs_time.xlsx")
COMBO03_METHOD = "spencer"
COMBO03_SLICES = 40
#: The instants the shipped march saves, which the Run LEM and Parametric dialogs
#: offer as a list. Named here rather than read off the march: both selectors read
#: the times off the solution and nothing else, and importing nineteen frames of
#: nodal head for a dialog shot is work the shot does not need. The producer in
#: ``make_tutorial_figures.py`` reads the real march and prints the same nineteen.
COMBO03_TIMES = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 47.0,
                 55.0, 65.0, 80.0, 100.0, 130.0, 180.0, 240.0, 300.0)
#: The frame the Run LEM shot is set to: the page's first single-instant run, the
#: full reservoir the march starts from. The dialog opens on the LAST saved frame,
#: so a shot on this one is the reader having changed it, which is what the page is
#: showing.
COMBO03_DIALOG_TIME = 0.0
#: The frame the play bar is parked on: three quarters of the way down the fall,
#: where the pool has left the upper face and the core still holds its head.
COMBO03_PLAYBAR_TIME = 35.0

_combo03_cache = {}


def _combo03_solved():
    """The shipped model with the mesh and the march it carries, one frame staged.

    The state the reader is in at the first stability run: the loader attaches the
    companion mesh, the march is read off ``_tseep.csv``, and one instant of it is
    placed on the model — which a file whose materials read ``u = seep`` needs
    before Run LEM will report anything but an error panel. Cached, because the
    import is the same work for every shot.
    """
    if "data" not in _combo03_cache:
        from xslope.seep import (build_seep_data, import_transient_solution,
                                 select_transient_frame_u)

        data = _load(COMBO03)
        with contextlib.redirect_stdout(io.StringIO()):
            seep_data = build_seep_data(data["mesh"], data, seep_bc=1)
            solution = import_transient_solution(seep_data,
                                                 os.path.splitext(COMBO03)[0])
            select_transient_frame_u(data, solution, time=COMBO03_DIALOG_TIME)
        _combo03_cache["data"] = data
        _combo03_cache["solution"] = solution
        _combo03_cache["seep_data"] = seep_data
    return _combo03_cache["data"]


def combo03_materials():
    """The two zones in table view with the LEM columns showing.

    The band this page reads out: a unit weight above the water table and a
    saturated one below it, a drained effective-stress envelope on each zone, and
    ``u`` on ``seep`` so every slice base reads its pore pressure from the solved
    field.
    """
    from studio.editors import MaterialsEditor

    dlg = _lem_only(MaterialsEditor().build(_load(COMBO03), None))
    return _grab(_mat_table(dlg, through="u"), "combo03_studio_materials.png")


def combo03_circles():
    """The circles editor on the dam: a starting circle on each face.

    The page's starting-circle step is read off this shot — two rows in the table,
    the two arcs the preview draws over the two slopes, and a **Search window**
    group holding one value, the minimum slip depth. Which face governs is what the
    curve measures, so the entry and exit ranges that would confine the search to
    one of them are left blank.
    """
    from studio.editors import CirclesEditor

    return _grab(CirclesEditor().build(_load(COMBO03), None),
                 "combo03_studio_circles.png")


def combo03_run_lem():
    """Run LEM on the shipped model, with the **Seepage time** group set to the
    instant the page runs.

    The group is what a march adds: an ordinary run reads ONE instant, and this is
    where that instant is named. Its **Saved frame** list holds the nineteen the
    march stored, the note under it restates which field the run will read, and
    the checkbox writes the choice to the model so a scripted re-run reads the
    same frame. The shot is taken on the frame the page's first run uses rather
    than on the dialog's opening one, because the page is showing the reader how to
    change it.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": COMBO03_METHOD,
                                 "analysis": "auto_search",
                                 "num_slices": COMBO03_SLICES},
                       slope_data=_combo03_solved(),
                       transient={"times": list(COMBO03_TIMES)})
    dlg.seep_time.frame.setCurrentIndex(COMBO03_TIMES.index(COMBO03_DIALOG_TIME))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo03_studio_run_lem.png")


def combo03_playbar():
    """The **Seep · Transient** results tab on the march the file ships.

    Every pore pressure the page reads comes from one of the frames under this play
    bar, and the reader gets them without solving anything: the workbook arrives
    with ``_tseep.csv`` beside it. Parked mid-drawdown, at the instant the page
    runs singly, where the pool has left the upper face and the core still holds
    its head. The display options are the transient panel's defaults — no flow
    lines, since a storage-release state has no flow net.
    """
    from studio.transient import TransientSeepView

    data = _combo03_solved()
    seep_data = _combo03_cache["seep_data"]
    solution = _combo03_cache["solution"]

    opts = {"variable": "head", "levels": 12, "flowlines": False,
            "vectors": True, "phreatic": True, "show_bc_levels": True}
    view = TransientSeepView()
    # The canvas sizes its figure to the viewport, and this dam is five times as
    # wide as it is tall, so a tall view leaves the frame stranded in white space.
    view.resize(1000, 540)
    view.set_frames(seep_data, solution["frames"],
                    opts_getter=lambda: opts, style_getter=lambda: None,
                    keep_index=False)
    view.show()
    _settle()
    times = [float(f["time"]) for f in solution["frames"]]
    view.set_index(min(range(len(times)),
                       key=lambda i: abs(times[i] - COMBO03_PLAYBAR_TIME)))
    _settle()
    view.canvas._render_current()          # force the raster into the scene
    _settle()
    out = os.path.join(OUT_DIR, "combo03_studio_playbar.png")
    view.grab().save(out)
    view.close()
    print("-> combo03_studio_playbar.png")
    return out


def combo03_parametric():
    """Run → Parametric… in Factor-of-safety-vs-time mode.

    The whole curve in one run: the parameter picker steps aside (no input is
    substituted at any point) and the march's nineteen saved frames take its place,
    every one ticked. Method and slice count are the single-instant run's, and
    every solver control is at its default — the two circles the file carries reach
    both faces of the dam, so nothing here has to be turned on.
    """
    from studio.dialogs import SensitivityDialog

    dlg = SensitivityDialog(defaults={"method": COMBO03_METHOD,
                                      "num_slices": COMBO03_SLICES},
                            slope_data=_combo03_solved(), app_mode="lem",
                            transient={"times": list(COMBO03_TIMES)})
    dlg.mode.setCurrentIndex(dlg.mode.findData("fs_vs_time"))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo03_studio_parametric.png")


def combo03_fs_time_result():
    """The **FS vs Time** result tab the run opens.

    The tab draws the factors of safety the run computed, so this shot costs the
    whole sweep — nineteen searches. The frames are read off the shipped
    ``_tseep.csv`` rather than marched again, which is what the reader's copy does
    too, and the sweep is the page's own: Spencer, searching at each instant.
    """
    from studio.main_window import SweepCanvas
    from studio.runners import SensitivityRunner

    data = _combo03_solved()
    solution = _combo03_cache["solution"]
    opts = {"mode": "fs_vs_time", "engine_mode": "lem", "method": COMBO03_METHOD,
            "num_slices": COMBO03_SLICES, "search": True,
            "times": [float(t) for t in solution["times"]]}
    runner = SensitivityRunner(data, opts, transient=solution)
    bundle, err = {}, {}
    runner.succeeded.connect(lambda b: bundle.update(b))
    runner.failed.connect(lambda m: err.setdefault("msg", m))
    with contextlib.redirect_stdout(io.StringIO()):
        runner._run_fs_vs_time()
    if err or not bundle:
        raise RuntimeError("combo03_fs_time_result: the sweep failed: %s"
                           % err.get("msg", "no result"))
    canvas = SweepCanvas()
    canvas.resize(1000, 620)
    canvas.show()
    _settle()
    canvas.render_fs_vs_time(bundle, slope_data=data)
    _settle()
    canvas._render_current()
    _settle()
    out = os.path.join(OUT_DIR, "combo03_studio_fs_time.png")
    canvas.grab().save(out)
    canvas.close()
    print("-> combo03_studio_fs_time.png")
    return out


#: Part 2's shipped model. It carries the mesh and the march as companions, so
#: loading it is the whole of the state the reader opens the dialog in, and the
#: loader attaches the mesh without a build.
COMBO03R_FILE = os.path.join(REPO_ROOT,
                             "docs/tutorials/files/xslope_johnson_fs_time.xlsx")

#: The instants the shipped march saves, which the Parametric dialog offers as a
#: checklist. Named here rather than read off the march: the frames list reads the
#: times off the solution and nothing else, and importing 12 frames of nodal head
#: for one dialog shot is work the shot does not need. The producer in
#: ``make_tutorial_figures.py`` reads the real march and prints the same twelve.
COMBO03R_TIMES = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
                  60.0, 70.0, 80.0, 100.0, 130.0, 170.0, 220.0, 300.0, 400.0, 500.0)


def combo03_rapid_parametric():
    """Run → Parametric… in Factor-of-safety-vs-time mode with **Rapid
    drawdown at each time** ticked, on Part 2's shipped model.

    The checkbox is Part 2's whole control: it sits under the frames list because
    it changes what a ticked instant means — stage 2 of a fall that began at the
    march's initial pool rather than a state on its own — and ticking it holds
    **Re-search the critical surface at each step** on and greys it, which the shot
    has to show along with the note the dialog rewrites underneath.

    Shot on the shipped model rather than on COMBO-2's, because the Model checks
    panel beside the controls is part of what the page reads: that file carries no
    boundary set 2, so the drawdown's one remaining warning is the free-draining
    one.
    """
    from studio.dialogs import SensitivityDialog

    dlg = SensitivityDialog(defaults={"method": COMBO03_METHOD,
                                      "num_slices": COMBO03_SLICES},
                            slope_data=_load(COMBO03R_FILE), app_mode="lem",
                            transient={"times": list(COMBO03R_TIMES)})
    dlg.mode.setCurrentIndex(dlg.mode.findData("fs_vs_time"))
    dlg.rapid.setChecked(True)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "combo03_rapid_parametric.png")


SHOTS.update({
    "combo03_materials": combo03_materials,
    "combo03_circles": combo03_circles,
    "combo03_run_lem": combo03_run_lem,
    "combo03_playbar": combo03_playbar,
    "combo03_parametric": combo03_parametric,
    "combo03_fs_time_result": combo03_fs_time_result,
    "combo03_rapid_parametric": combo03_rapid_parametric,
})


# --------------------------------------------------------------------------- #
# LEM-13 — Rock Slope (Hoek-Brown)
#
# The materials shot is LIST view: only the list view puts the four Hoek-Brown
# field inputs beside the derived mb/s/a readout and the envelope those constants
# draw, which is the whole subject of Part A. Part B photographs the Parametric
# dialog in its Design mode, set to the GSI sweep the page runs.
# --------------------------------------------------------------------------- #
LEM13_A = os.path.join(REPO_ROOT, "docs/tutorials/files/xslope_rock_slope.xlsx")


def _lem13_material(model, name, edit=None, width=1180, height=780):
    """The materials editor in list view, on the one rock material.

    ``edit`` is applied before the editor is built, so the shot carries the state
    the step ends in rather than a pre-edit row described in prose.
    """
    from studio.editors import MaterialsEditor

    d = _load(model)
    if edit:
        d = dict(d, materials=[dict(d["materials"][0], **edit)])
    dlg = _lem_only(MaterialsEditor().build(d, None))
    dlg.set_view_mode("list")
    dlg._list_view.list.setCurrentRow(0)
    dlg.resize(width, height)
    return _grab(dlg, name)


def lem13_materials():
    """Part A's rock as the file carries it: option `hb`, the four field inputs, and
    the mb / s / a the editor derives from them beside the envelope they define."""
    return _lem13_material(LEM13_A, "lem13_studio_materials.png")


def lem13_run_lem():
    """Part A's run: Spencer, auto search, the dialog's own 40 slices."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM13_A))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem13_studio_run_lem.png")


def lem13_run_fem():
    """Part A's strength reduction, on the meshed model, at the bracket the
    completed file declares.

    The mesh is built here because the file ships no sidecar — Studio's Run FEM
    action is unreachable without one, so the dialog is photographed in the state
    the Build Mesh step leaves behind.
    """
    from studio.dialogs import RunFemDialog
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    data = _load(LEM13_A)
    with contextlib.redirect_stdout(io.StringIO()):
        data["mesh"] = build_mesh_from_polygons(
            get_material_polygons(data), data["target_size"], data["element_type"],
            size_regions=extract_size_regions(data))
    dlg = RunFemDialog(defaults={"analysis": "ssrm",
                                 "F_min": float(data["ssrm_f_min"]),
                                 "F_max": float(data["ssrm_f_max"]),
                                 "tolerance": 0.01,
                                 "k0": float(data["k0"])},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem13_studio_run_fem.png")


def lem13_build_mesh():
    """Build Mesh on the size the completed file declares: tri6 at 0.9 m, with
    auto-sizing off because a declared size is what turns it off."""
    from studio.dialogs import BuildMeshDialog

    data = _load(LEM13_A)
    dlg = BuildMeshDialog(defaults={"element_type": data["element_type"],
                                    "target_size": float(data["target_size"]),
                                    "auto_size": False})
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem13_studio_build_mesh.png")


def lem13_run_lem_corps():
    """Run LEM with the Corps of Engineers method selected on the Hoek-Brown rock.

    The method-and-material pairing raises its own check: a fixed-inclination
    force-equilibrium method against an envelope whose instantaneous friction angle
    can pass 55 degrees at low confinement. Nothing else on the dialog changes.
    """
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "corps", "analysis": "auto_search",
                                 "num_slices": 40},
                       slope_data=_load(LEM13_A))
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem13_studio_run_lem_corps.png")


def lem13_parametric():
    """The Parametric dialog set up for Part B's GSI sweep.

    Design mode rather than Sensitivity: the page wants the whole curve across an
    explicit range, not a percentage band about the current value. The parameter
    picker is driven to ``mat:rock:hb_gsi`` the same way a reader drives it — pick
    the material, then the property — so the shot carries the reference the run
    actually sweeps, echoed on the **Sweeping** row.
    """
    from studio.dialogs import SensitivityDialog

    data = _load(LEM13_A)
    dlg = SensitivityDialog(defaults={"mode": "design", "method": "spencer",
                                     "num_slices": 40, "low": 5.0, "high": 20.0,
                                     "steps": 6, "target_fs": 1.5,
                                     "search": True},
                           slope_data=data, app_mode="lem")
    i = dlg.prop.findData("mat:rock:hb_gsi")
    if i < 0:
        raise SystemExit("the Parametric dialog cannot reference mat:rock:hb_gsi")
    dlg.prop.setCurrentIndex(i)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "lem13_studio_parametric.png")


SHOTS.update({
    "lem13_materials": lem13_materials,
    "lem13_run_lem": lem13_run_lem,
    "lem13_build_mesh": lem13_build_mesh,
    "lem13_run_fem": lem13_run_fem,
    "lem13_run_lem_corps": lem13_run_lem_corps,
    "lem13_parametric": lem13_parametric,
})


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    names = [n for n in SHOTS if not argv or any(a in n for a in argv)]
    if not names:
        print("no shot matching %s; known shots: %s" % (argv, ", ".join(sorted(SHOTS))))
        return 1
    with _app_defaults():
        for name in names:
            SHOTS[name]()
    print("\nwrote %d screenshot(s) to docs/tutorials/images/" % len(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
