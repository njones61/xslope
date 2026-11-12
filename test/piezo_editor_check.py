"""Standing checks on the Studio Piezometric lines editor's Type, one per line.

Each piezometric line carries a Type, the piezo sheet's Type cell (B3 for line 1, E3 for
line 2): piezo, where u is the full static head below the line, or phreatic, where the
static head is reduced by cos^2 of the line's slope. Until 2026-11 the editor showed only
the points, so a phreatic line could be set only in the workbook, and the CE 544
seepage deck told students so. Now each tab has a Type box above its table.

  A. THE EDITOR SHOWS EACH LINE'S TYPE. Opened on a model whose line 1 is phreatic and
     whose line 2 is piezo, the two boxes read phreatic and piezo, and the other way
     round.
  B. OK WRITES THE CHOICE. Flipping both boxes and pressing OK sets piezo_phreatic and
     piezo_phreatic2 to match.
  C. OPEN AND OK CHANGES NOTHING. Untouched, the editor writes back the Types it was
     opened on, in all four combinations.
  D. THE TYPE SURVIVES A SAVE. The model saved and read back carries the Types the
     editor set, and the piezo sheet's B3 and E3 read 'phreatic' or 'piezo'.
  E. THE TABS KEEP THEIR NAMES. Tutorial LEM-4 names them: "Line 1" and
     "Line 2 (rapid drawdown)".

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A shipped model with a piezometric line on line 1.
MODEL = os.path.join(_REPO, "docs/lem/files/xslope_gsat_piezo.xlsx")
TABS = ("Line 1", "Line 2 (rapid drawdown)")


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _boxes(dlg):
    """The Type each tab's box shows, as the bool it stands for."""
    return tuple(bool(dlg.tab_prop(i)) for i in range(2))


def _set(dlg, i, phreatic):
    combo, values = dlg._tab_combos[i]
    combo.setCurrentIndex(values.index(phreatic))


def run():
    try:
        from PySide6.QtWidgets import QApplication, QTabWidget
    except ImportError:
        print("  PySide6 not installed: skipped")
        return []
    QApplication.instance() or QApplication([])
    from studio.editors import PiezoEditor
    fails = []
    base = _load(MODEL)

    # A and C: what the boxes show, and that OK untouched writes it back
    for types in ((True, False), (False, True), (False, False), (True, True)):
        data = dict(base, piezo_phreatic=types[0], piezo_phreatic2=types[1])
        ed = PiezoEditor()
        dlg = ed.build(data, None)
        shown = _boxes(dlg)
        if shown != types:
            fails.append(f"A. opened on Types {types}, the boxes read {shown}")
        out = dict(data)
        ed.apply(out, dlg)
        if (out["piezo_phreatic"], out["piezo_phreatic2"]) != types:
            fails.append(f"C. opened on {types} and OK'd untouched, the editor wrote "
                         f"{(out['piezo_phreatic'], out['piezo_phreatic2'])}")
        dlg.deleteLater()
        print(f"  opened on {types}: boxes {shown}, OK untouched writes "
              f"{(out['piezo_phreatic'], out['piezo_phreatic2'])}")

    # B and D: flip both, OK, save, read back
    from xslope.fileio import save_slope_data_to_xlsx
    from openpyxl import load_workbook
    data = dict(base, piezo_phreatic=False, piezo_phreatic2=True)
    ed = PiezoEditor()
    dlg = ed.build(data, None)
    _set(dlg, 0, True)
    _set(dlg, 1, False)
    out = dict(data)
    ed.apply(out, dlg)
    tabs = [t.tabText(i) for t in dlg.findChildren(QTabWidget) for i in range(t.count())]
    dlg.deleteLater()
    if (out["piezo_phreatic"], out["piezo_phreatic2"]) != (True, False):
        fails.append(f"B. flipped to (True, False), OK wrote "
                     f"{(out['piezo_phreatic'], out['piezo_phreatic2'])}")
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "piezo_type.xlsx")
        with contextlib.redirect_stdout(io.StringIO()):
            save_slope_data_to_xlsx(out, path)
        back = _load(path)
        ws = load_workbook(path)["piezo"]
        cells = (ws["B3"].value, ws["E3"].value)
    got = (bool(back.get("piezo_phreatic")), bool(back.get("piezo_phreatic2")))
    print(f"  flipped to (True, False): OK wrote {(out['piezo_phreatic'], out['piezo_phreatic2'])}, "
          f"saved B3/E3 {cells}, read back {got}")
    if got != (True, False):
        fails.append(f"D. saved and read back, the Types are {got}")
    if cells != ("phreatic", "piezo"):
        fails.append(f"D. the piezo sheet's B3 and E3 read {cells}")

    # E: the tab names
    for name in TABS:
        if name not in tabs:
            fails.append(f"E. the editor has no {name!r} tab; its tabs read {tabs}")
    return fails


if __name__ == "__main__":
    f = run()
    print("\n".join(["FAILED:"] + f) if f else "\nAll piezometric-line editor checks passed.")
    sys.exit(1 if f else 0)
