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

    dlg = MaterialsEditor().build(_load(LEM01), None)
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

    dlg = MaterialsEditor().build(_load(LEM03), None)
    dlg._set_mode("table")
    dlg.resize(1180, 380)
    return _grab(dlg, "lem03_studio_materials.png")


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
    return _grab(_list_view(dlg, 1400), "lem08_studio_reinforcement.png")


SHOTS["lem08_reinforcement"] = lem08_reinforcement


# --------------------------------------------------------------------------- #
# LEM-9 — A Tieback Wall
# --------------------------------------------------------------------------- #
LEM09 = os.path.join(REPO_ROOT,
                     "docs/verification/files/rocscience/vp049.xlsx")


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
    return _grab(_list_view(dlg, 1400), "lem09_studio_reinforcement.png")


def lem09_piles():
    """The piles editor on the soldier pile, in its list view.

    Built on LEM-9's model, so the form holds the pile the page enters — the
    vertical axis at the wall face, the shear force per foot of wall, and the
    diameter and spacing beside it — rather than an empty row.
    """
    from studio.editors import PilesEditor

    dlg = PilesEditor().build(_load(LEM09), None)
    return _grab(_list_view(dlg, 1400), "lem09_studio_piles.png")


SHOTS["lem09_reinforcement"] = lem09_reinforcement
SHOTS["lem09_piles"] = lem09_piles


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    names = [n for n in SHOTS if not argv or any(a in n for a in argv)]
    if not names:
        print("no shot matching %s; known shots: %s" % (argv, ", ".join(sorted(SHOTS))))
        return 1
    for name in names:
        SHOTS[name]()
    print("\nwrote %d screenshot(s) to docs/tutorials/images/" % len(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
