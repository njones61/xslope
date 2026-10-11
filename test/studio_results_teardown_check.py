"""Standing check on Studio's FEM results view being closed and the session going on.

The FEM results view carries its own toolbar ("1D Details…" and "Continue with a
higher limit…"), and the main window keeps a handle on each button so it can dim or
hide it. Three ordinary actions close the view because its results are stale: OK on
an input editor, a mesh rebuild, and opening another file. The buttons go with the
view. A handle left pointing at a deleted button raises the next time the Run menu
updates, which happens after every one of those actions, and the update stops part
way: Run FEM stays disabled for the rest of the session. Nothing about the model is
wrong, so the only symptom is a reader who cannot make a second run.

Each leg opens a results view, closes it the way the action does, lets Qt delete the
widgets, and then updates the Run menu exactly as the window does.

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

#: FEM-5's toppling model: a FEM model with joint lines, the case that met this.
MODEL = os.path.join(_REPO, "docs/tutorials/files/xslope_rock_joints.xlsx")
#: The file the third leg opens next.
OTHER = os.path.join(_REPO, "docs/tutorials/files/xslope_pile_wall.xlsx")


def _flush():
    """Run Qt's pending deletions, so a deleted button is really gone."""
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    QApplication.instance().processEvents()


def _leg(mw, name, close):
    """Open a results view, close it with ``close``, then update the Run menu."""
    fails = []
    with contextlib.redirect_stdout(io.StringIO()):
        mw._show_fem_results()
    if mw.fem_continue_btn is None or mw.fem_details_btn is None:
        return [f"{name}: the FEM results view opened without its toolbar buttons"]
    with contextlib.redirect_stdout(io.StringIO()):
        close()
    _flush()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mw._update_run_actions()
            mw._update_fem_details_action()
    except RuntimeError as exc:
        fails.append(f"{name}: the Run menu update raised after the results view "
                     f"closed ({exc}); Run FEM stays disabled for the session")
    if mw.fem_continue_btn is not None or mw.fem_details_btn is not None:
        fails.append(f"{name}: the window still holds the closed view's toolbar "
                     f"buttons")
    return fails


def run():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    QApplication.instance() or QApplication([])
    from studio.main_window import MainWindow

    mw = MainWindow()
    fails = []
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mw.open_path(MODEL)
        legs = [
            ("OK on an input editor", lambda: mw.invalidate_results()),
            ("a mesh rebuild", lambda: mw.invalidate_results(clear_mesh=True)),
            ("opening another file", lambda: mw.open_path(OTHER)),
        ]
        for name, close in legs:
            got = _leg(mw, name, close)
            print(f"  {name:<40s} {'ok' if not got else 'FAIL'}")
            fails += got
    finally:
        mw.doc._dirty = False
        mw.close()
    return fails


def main():
    print("FEM results view closed mid-session:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll results-view checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
