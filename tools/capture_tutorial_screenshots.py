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
    """The Run LEM dialog on the complete model, set up for the tutorial's first run:
    Bishop, a single surface, 40 slices — and a clean model-checks column."""
    from studio.dialogs import RunLemDialog

    dlg = RunLemDialog(defaults={"method": "bishop", "analysis": "single_surface",
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

    The dock's caption reads the provider and model stored on the machine that runs
    this, which is the dock's own behaviour ("the dock shows the active provider ·
    model") rather than anything the shot arranges — so a re-capture elsewhere shows
    that machine's choice, and the page never quotes the caption's contents.
    """
    from studio.main_window import MainWindow
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
    print("-> t0_studio_window.png  (%dx%d, offscreen main window)"
          % (pix.width(), pix.height()))
    win.close()
    return out


SHOTS["t0_studio_window"] = t0_studio_window


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
