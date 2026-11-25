"""Standing check: opening a file starts Studio's run dialogs on the run options it
declares (Norm 646a).

A workbook's main sheet may declare how it is meant to be run: the LEM method, the
mesh target size, the strength-reduction bracket, and so on. Studio's Run LEM, Build
Mesh and Run FEM dialogs open on those values, and a choice made in a dialog then
sticks for the rest of the work on that file. Before 646a the remembered choices
outlived the file: after a 0.8 m mesh on FEM-3's block wall, the geotextile files
that declare 1.2 m opened Build Mesh at 0.8 m, and FEM-3 Part 3 came out at 1.535 /
2.128 / 1.277 where the page's 1.2 m mesh gives 1.566 / 2.167 / 1.285. The same held
for a declared LEM method and SSRM bracket. Now opening a file (or starting a new
one, or importing) empties the remembered choices, so each dialog opens on what the
new file declares.

Each case: open file A, make a choice in the dialog's memory as a run does, check the
choice sticks on A (the dialog must not snap back to the file mid-session), then open
file B and check the dialog opens on B's declaration, not on A's choice.

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

TUT = os.path.join(_REPO, "docs/tutorials/files")
BLOCK_WALL = os.path.join(TUT, "xslope_block_wall.xlsx")          # declares 0.8 m
LINER = os.path.join(TUT, "xslope_liner_bonded.xlsx")             # declares 1.2 m
W1 = os.path.join(TUT, "w1_build_from_image_after.xlsx")          # declares spencer
PILE_WALL = os.path.join(TUT, "xslope_pile_wall.xlsx")            # declares F max 2.0
ROCK = os.path.join(TUT, "xslope_rock_slope.xlsx")                # declares F max 1.6


def _open(mw, path):
    mw.doc._dirty = False
    with contextlib.redirect_stdout(io.StringIO()):
        mw.open_path(path)


def _case(mw, which, file_a, choice, file_b, key, want_b):
    fails = []
    _open(mw, file_a)
    cache = {"mesh": "_last_mesh_opts", "lem": "_last_lem_opts", "fem": "_last_fem_opts"}[which]
    setattr(mw, cache, dict(choice))               # what a run on file A leaves behind
    got_a = mw._file_defaults(which).get(key)
    if got_a != choice[key]:
        fails.append(f"{which}: on the same file the dialog opens on {key} = {got_a!r}, "
                     f"not the choice just made ({choice[key]!r})")
    _open(mw, file_b)
    got_b = mw._file_defaults(which).get(key)
    if got_b != want_b:
        fails.append(f"{which}: after opening {os.path.basename(file_b)} the dialog opens "
                     f"on {key} = {got_b!r}, not the {want_b!r} the file declares "
                     f"(the last file's choice was {choice[key]!r})")
    return fails


def test_mesh():
    from studio.main_window import MainWindow
    mw = MainWindow()
    try:
        return _case(mw, "mesh", BLOCK_WALL,
                     {"element_type": "tri6", "auto_size": False, "target_size": 0.5},
                     LINER, "target_size", 1.2)
    finally:
        mw.close()


def test_lem():
    from studio.main_window import MainWindow
    mw = MainWindow()
    try:
        return _case(mw, "lem", BLOCK_WALL, {"method": "janbu", "num_slices": 30},
                     W1, "method", "spencer")
    finally:
        mw.close()


def test_fem():
    from studio.main_window import MainWindow
    mw = MainWindow()
    try:
        return _case(mw, "fem", PILE_WALL, {"F_min": 0.5, "F_max": 3.0},
                     ROCK, "F_max", 1.6)
    finally:
        mw.close()


CHECKS = [("Build Mesh opens on the new file's size", test_mesh),
          ("Run LEM opens on the new file's method", test_lem),
          ("Run FEM opens on the new file's bracket", test_fem)]


@studio_settings.isolated
def run():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    QApplication.instance() or QApplication([])
    fails = []
    for name, fn in CHECKS:
        try:
            got = fn()
        except Exception as e:                        # noqa: BLE001
            import traceback
            traceback.print_exc()
            got = [f"{name}: raised {type(e).__name__}: {e}"]
        print(f"  {name:<46s} {'ok' if not got else 'FAIL'}")
        fails += got
    return fails


def main():
    print("Run dialogs open on a newly opened file's run options:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll dialog-default checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
