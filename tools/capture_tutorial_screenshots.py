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
    dlg._set_mode("table")
    # The preview stacks below the table, so the width is just the columns:
    # wide enough that every LEM column (Label through Spacing) shows.
    dlg.resize(1460, 760)
    return _grab(dlg, "lem08_studio_reinforcement_table.png")


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
    dlg._set_mode("table")
    dlg.resize(1460, 760)
    return _grab(dlg, "lem09_studio_reinforcement_table.png")


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
    dlg._set_mode("table")
    dlg.resize(1460, 700)
    return _grab(dlg, "lem09_studio_piles_table.png")


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
    """The same dialog with the section's two tools on: grid seeding and the
    surficial filter at 5 ft — the state the tools section tells the reader to
    put it in."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "spencer", "analysis": "auto_search",
                                 "num_slices": 40, "grid_seed": True,
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


SHOTS["lem10_run_lem"] = lem10_run_lem
SHOTS["lem10_run_lem_filtered"] = lem10_run_lem_filtered
SHOTS["lem10_circles"] = lem10_circles


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

    return _grab(_mat_table(MaterialsEditor().build(_load(LEM04), None)),
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

    return _grab(_mat_table(MaterialsEditor().build(_load(LEM05), None)),
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
