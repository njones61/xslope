"""Standing checks on Studio's analysis-mode switch — the segmented strip in the
toolbar that selects LEM / Seepage / FEM.

The mode is the most-used control in the window: it decides what the Inputs view
draws and emphasizes, what the Run button runs, whether Build Mesh is there at
all, and which study dialogs apply. It is chosen dozens of times in a session and
is stored in no input file, so a strip that looked right but reached nothing would
leave no trace anywhere else — the window would simply stop following the click.
What is guarded here:

  A. THE STRIP — one checkable segment per entry in ``MODES``, labelled with that
     entry's label, in one exclusive group with exactly one segment checked, LEM on
     a fresh window. The segments sit in a zero-spacing container under an explicit
     padding rule rather than relying on the platform style, which is the whole
     point of the control: Windows and macOS space toolbar widgets differently.
  B. THE SWITCH — a click on a segment leaves the window in exactly the state a
     direct call to the mode handler leaves it in: the mode itself, the Run action's
     label / enablement, Build Mesh's visibility, the Parametric and Reliability
     gates, the mode the canvas is re-rendered with, and the inputs tree being
     rebuilt. This is the check that fails if the group is wired to nothing — the
     highlight would still move and the window would not. Re-picking the mode already
     in force stays a no-op, as it was when this was a drop-down that signalled only
     on a change of index.
  C. THE SHORTCUTS — Ctrl+1..N in ``MODES`` order (⌘ on macOS, which is Qt's own
     mapping of Ctrl), delivered as real key events, each moving both the mode and
     the highlight; and every segment names its shortcut in its tooltip, which is
     the only place a user finds out the shortcut exists.
  D. THE RELOAD — opening a file selects the mode that fits it and the highlight
     follows, silently: the file's own mode must not run the mode-changed handler a
     second time on top of the load's own render.

Nothing here solves anything: the models are opened and rendered, and no analysis
is run.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A plain limit-equilibrium model — opens in LEM mode.
DAM_LEM = os.path.join(_REPO, "docs/inputs/slope/xslope_dam.xlsx")
#: Materials with hydraulic properties and no strength — opens in Seepage mode.
LEVEE_SEEP = os.path.join(_REPO, "docs/seep/files/xslope_levee_poly.xlsx")


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        return fn(*args, **kwargs)


def _window(path=DAM_LEM):
    """A real Studio window with ``path`` loaded.

    ``doc.load`` rather than ``open_path``: opening also writes the file into the
    user's recent-files list, and a test must not edit anyone's settings."""
    from studio.main_window import MainWindow

    mw = MainWindow()
    if path is not None:
        _quiet(mw.doc.load, path)
    return mw


def _close(mw):
    mw.doc._dirty = False               # never raise the save prompt on a test
    mw.close()


def _instrument(mw):
    """Record the side effects a mode change is supposed to have.

    The canvas render and the inputs-tree rebuild are the two the window performs
    on itself, so they are recorded by standing in for them; everything else is
    read back off the window afterwards."""
    seen = {"render_modes": [], "tree_rebuilds": 0}
    real_render = mw.canvas.render_inputs

    def render_inputs(sd, mode=None, **kwargs):
        seen["render_modes"].append(mode)
        return real_render(sd, mode=mode, **kwargs)

    def populate():
        seen["tree_rebuilds"] += 1

    mw.canvas.render_inputs = render_inputs
    mw._populate_inputs_tree = populate
    return seen


def _state(mw, seen):
    """Everything a mode change is observable through, as one comparable value."""
    return {
        "mode": mw._mode,
        "checked": mw.mode_buttons.checkedId(),
        "run_text": mw.act_run.text(),
        "run_enabled": mw.act_run.isEnabled(),
        "build_mesh_visible": mw.act_build_mesh.isVisible(),
        "build_mesh_enabled": mw.act_build_mesh.isEnabled(),
        "sensitivity_enabled": mw.act_sensitivity.isEnabled(),
        "reliability_visible": mw.act_reliability.isVisible(),
        "render_modes": list(seen["render_modes"]),
        "tree_rebuilds": seen["tree_rebuilds"],
    }


# ------------------------------------------------------------------- A. strip
def test_strip():
    fails = []
    from studio.main_window import MODES

    mw = _window(path=None)
    try:
        buttons = mw.mode_buttons.buttons()
        if len(buttons) != len(MODES):
            fails.append(f"the mode strip has {len(buttons)} segments, "
                         f"want {len(MODES)} (one per MODES entry)")
        labels = [mw.mode_buttons.button(i).text() for i in range(len(MODES))
                  if mw.mode_buttons.button(i) is not None]
        want = [label for label, _ in MODES]
        if labels != want:
            fails.append(f"the segments read {labels}, want {want} in MODES order")
        if not mw.mode_buttons.exclusive():
            fails.append("the mode group is not exclusive — two modes could be lit")
        for b in buttons:
            if not b.isCheckable():
                fails.append(f"the {b.text()!r} segment is not checkable — it can "
                             "never show which mode is in force")
        checked = [b.text() for b in buttons if b.isChecked()]
        if checked != ["LEM"]:
            fails.append(f"a fresh window lights {checked}, want ['LEM'] "
                         "(the mode a new window starts in)")

        # A drop-down is what this replaced; leaving one behind would be a second
        # source of truth for the same state.
        if hasattr(mw, "mode_combo"):
            fails.append("the window still carries a mode_combo — two controls "
                         "for one mode")

        # The strip has to look like one strip on every platform, so the segments
        # are laid out with no gap and padded by an explicit rule rather than by
        # whatever the native style would choose.
        box = buttons[0].parent()
        layout = box.layout()
        if layout is None or layout.spacing() != 0:
            got = "no layout" if layout is None else layout.spacing()
            fails.append(f"the segments are spaced {got} px apart, want 0 — they "
                         "would read as separate buttons")
        if layout is not None:
            m = layout.contentsMargins()
            if (m.left(), m.top(), m.right(), m.bottom()) != (0, 0, 0, 0):
                fails.append("the segment container carries margins, so the strip "
                             "sits loose in the toolbar")
    finally:
        _close(mw)
    return fails


# ------------------------------------------------------------------ A2. styling
def test_styling():
    """The strip's metrics are spelled out, not inherited from the platform."""
    fails = []
    from PySide6.QtWidgets import QToolBar

    mw = _window(path=None)
    try:
        bar = None
        for tb in mw.findChildren(QToolBar):
            if tb.objectName() == "main_toolbar":
                bar = tb
        if bar is None:
            fails.append("no main_toolbar to carry the segment styling")
            return fails
        css = bar.styleSheet()
        if "mode_segments" not in css:
            fails.append("the toolbar stylesheet says nothing about the mode "
                         "segments — their look is the platform's to choose")
        if "padding" not in css:
            fails.append("no explicit padding on the segments; macOS and Windows "
                         "pad a text tool button differently")
        if "checked" not in css:
            fails.append("no :checked rule — the current mode would not stand out")
        # Every segment must be reachable by the stylesheet's own selector.
        for b in mw.mode_buttons.buttons():
            if b.property("segment") not in ("first", "middle", "last"):
                fails.append(f"the {b.text()!r} segment carries no position "
                             "property, so its end of the strip is not rounded")
    finally:
        _close(mw)
    return fails


# ------------------------------------------------------------------ B. switch
def test_switch():
    fails = []
    from studio.main_window import MODES

    # The handler is the specification: what a click does must be what calling it
    # does. Two identical windows, one driven each way.
    for index, (label, key) in enumerate(MODES):
        if key == "lem":
            continue                     # already the mode a fresh window is in
        by_hand = _window()
        by_click = _window()
        try:
            hand_seen = _instrument(by_hand)
            click_seen = _instrument(by_click)
            by_hand._on_mode_changed(index)
            by_click.mode_buttons.button(index).click()
            want = _state(by_hand, hand_seen)
            got = _state(by_click, click_seen)
            # The direct call does not move the highlight (the drop-down set that
            # itself); everything else must match exactly.
            want["checked"] = index
            for field in sorted(want):
                if got[field] != want[field]:
                    fails.append(
                        f"clicking {label!r} left {field}={got[field]!r}, but the "
                        f"mode handler leaves {field}={want[field]!r}")
            if got["render_modes"] != [key]:
                fails.append(f"clicking {label!r} re-rendered the canvas with "
                             f"{got['render_modes']}, want ['{key}']")
            if got["tree_rebuilds"] != 1:
                fails.append(f"clicking {label!r} rebuilt the inputs tree "
                             f"{got['tree_rebuilds']} times, want 1")
        finally:
            _close(by_hand)
            _close(by_click)

    # Exclusivity, walked across every mode: one lit segment, and it is the mode.
    mw = _window()
    try:
        for index, (label, key) in enumerate(MODES):
            mw.mode_buttons.button(index).click()
            lit = [b.text() for b in mw.mode_buttons.buttons() if b.isChecked()]
            if lit != [label]:
                fails.append(f"after picking {label!r} the strip lights {lit}")
            if mw._mode != key:
                fails.append(f"clicking {label!r} left the window in "
                             f"{mw._mode!r} mode")

        # Re-picking the mode already in force changes nothing — the drop-down
        # signalled on a change of index, not on every pick.
        seen = _instrument(mw)
        mw.mode_buttons.button(len(MODES) - 1).click()
        if seen["render_modes"] or seen["tree_rebuilds"]:
            fails.append("re-picking the current mode re-rendered the window")
        if not mw.mode_buttons.button(len(MODES) - 1).isChecked():
            fails.append("re-picking the current mode dropped its highlight")
    finally:
        _close(mw)
    return fails


# --------------------------------------------------------------- C. shortcuts
def test_shortcuts():
    fails = []
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QKeySequence
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication
    from studio.main_window import MODES

    mw = _window()
    try:
        mw.show()
        mw.activateWindow()
        for _ in range(5):
            QApplication.processEvents()

        for index, (label, key) in enumerate(MODES):
            if index >= 9:
                continue                 # Ctrl+1..9 is all the digits there are
            seq = QKeySequence(f"Ctrl+{index + 1}")
            native = seq.toString(QKeySequence.NativeText)
            tip = mw.mode_buttons.button(index).toolTip()
            if native not in tip:
                fails.append(f"the {label!r} segment's tooltip ({tip!r}) does not "
                             f"name its {native} shortcut")
            # Come from a different mode every time, so a key that does nothing
            # cannot be mistaken for one that worked.
            other = (index + 1) % len(MODES)
            mw.set_mode_index(other)
            QTest.keyClick(mw, getattr(Qt, f"Key_{index + 1}"), Qt.ControlModifier)
            for _ in range(5):
                QApplication.processEvents()
            if mw._mode != key:
                fails.append(f"{native} left the window in {mw._mode!r} mode, "
                             f"want {key!r}")
            if mw.mode_buttons.checkedId() != index:
                fails.append(f"{native} switched the mode without moving the "
                             "highlight onto that segment")
    finally:
        _close(mw)
    return fails


# ------------------------------------------------------------------ D. reload
def test_reload():
    fails = []
    from studio.main_window import MODES

    keys = [m for _, m in MODES]
    mw = _window(path=None)
    try:
        calls = []
        real = mw._on_mode_changed
        mw._on_mode_changed = lambda i: (calls.append(i), real(i))[1]

        for path, want in ((LEVEE_SEEP, "seep"), (DAM_LEM, "lem")):
            del calls[:]
            _quiet(mw.doc.load, path)
            name = os.path.basename(path)
            if mw._mode != want:
                fails.append(f"{name} opened in {mw._mode!r} mode, want {want!r}")
            if mw.mode_buttons.checkedId() != keys.index(mw._mode):
                lit = [b.text() for b in mw.mode_buttons.buttons() if b.isChecked()]
                fails.append(f"{name} opened in {mw._mode!r} mode with {lit} lit — "
                             "the strip disagrees with the window")
            if calls:
                fails.append(f"opening {name} ran the mode handler {len(calls)} "
                             "extra time(s) on top of the load's own render")
    finally:
        del mw._on_mode_changed         # back to the window's own bound method
        _close(mw)
    return fails


CHECKS = [("the segments + the fresh-window mode", test_strip),
          ("explicit, cross-platform strip metrics", test_styling),
          ("a click drives the mode handler", test_switch),
          ("Ctrl+N shortcuts + their tooltips", test_shortcuts),
          ("the mode a reopened file selects", test_reload)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("mode segments: PySide6 not installed — checks skipped.")
        return []
    failures = []
    for name, fn in CHECKS:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("analysis-mode switch (Studio toolbar):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll analysis-mode switch checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
