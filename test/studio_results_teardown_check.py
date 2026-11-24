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

The fourth check is the other side of the first: OK on an editor that changed
nothing closes nothing. Every editor used to drop every result on OK, so opening
Materials just to read it took a loaded transient seepage solution with it, and a
Save afterwards deleted that solution's files beside the workbook. Each editor is
opened and accepted unchanged on the transient earth dam: the results stay, no undo
step is added, the file stays clean, and a redo branch survives. A real edit still
drops the results.

Every window here is built against a scratch settings store
(``tools.studio_settings``): opening a file appends it to Studio's recent-files list,
and a check must never reorder the person's own.

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

from tools import studio_settings                      # noqa: E402

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


@studio_settings.isolated
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
        got = _noop_ok(mw)
        print(f"  {'OK on an editor that changed nothing':<40s} "
              f"{'ok' if not got else 'FAIL'}")
        fails += got
    finally:
        mw.doc._dirty = False
        mw.close()
    return fails


#: The transient earth dam: circles, profile lines, seepage blocks, a transient
#: schedule and a saved transient solution, so most editors have something to show.
NOOP_MODEL = os.path.join(_REPO, "docs/tutorials/files/xslope_earth_dam_fs_time.xlsx")


def _noop_ok(mw):
    """OK on every editor with nothing changed keeps every result and the history;
    a real edit still drops them."""
    import shutil
    import tempfile
    from PySide6.QtWidgets import QDialog
    from studio.main_window import CATEGORY_EDITORS
    fails = []
    tmp = tempfile.mkdtemp(prefix="noop_ok_")
    base = os.path.splitext(os.path.basename(NOOP_MODEL))[0]
    for f in os.listdir(os.path.dirname(NOOP_MODEL)):
        if f.startswith(base + ".") or f.startswith(base + "_"):
            shutil.copy(os.path.join(os.path.dirname(NOOP_MODEL), f), tmp)
    model = os.path.join(tmp, os.path.basename(NOOP_MODEL))
    real_exec = QDialog.exec
    QDialog.exec = lambda self: 1          # every editor accepted as opened
    try:
        for cat in CATEGORY_EDITORS:
            mw.doc._dirty = False
            with contextlib.redirect_stdout(io.StringIO()):
                mw.open_path(model)
            keep = object()
            mw.doc.results["transient_seep"] = keep
            n_undo = len(mw.doc._undo)
            with contextlib.redirect_stdout(io.StringIO()):
                mw.edit_category(cat)
            if mw.doc.results.get("transient_seep") is not keep:
                fails.append(f"OK on the unchanged {cat} editor dropped the results")
            if len(mw.doc._undo) != n_undo or mw.doc.dirty:
                fails.append(f"OK on the unchanged {cat} editor added an undo step "
                             f"or marked the file changed")
        # A redo branch survives an unchanged OK; a real edit still drops results.
        mw.doc._dirty = False
        with contextlib.redirect_stdout(io.StringIO()):
            mw.open_path(model)
        mats = mw.doc.slope_data["materials"]
        mw.doc.begin_edit("test edit")
        mats[0]["c"] = float(mats[0].get("c") or 0.0) + 1.0
        mw.doc.commit_edit()
        mw.doc.undo()
        if not mw.doc._redo:
            fails.append("the test could not set up a redo branch")
        with contextlib.redirect_stdout(io.StringIO()):
            mw.edit_category("materials")
        if not mw.doc._redo:
            fails.append("OK on the unchanged materials editor lost the redo branch")
        keep = object()
        mw.doc.results["transient_seep"] = keep
        editor = CATEGORY_EDITORS["materials"]
        real_apply = editor.apply

        def bump(sd, dlg):
            real_apply(sd, dlg)
            sd["materials"][0]["c"] = float(sd["materials"][0].get("c") or 0.0) + 1.0
        editor.apply = bump
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mw.edit_category("materials")
        finally:
            editor.apply = real_apply
        if mw.doc.results.get("transient_seep") is keep:
            fails.append("a real Materials edit kept a result it should have dropped")
    finally:
        QDialog.exec = real_exec
        mw.doc._dirty = False
        shutil.rmtree(tmp, ignore_errors=True)
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
