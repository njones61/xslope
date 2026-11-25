"""Standing check on four Studio editor defects the 1.1.0 tutorial pass met.

1. **A new profile line or polygon starts with no vertex rows.** It used to start
   with one (0, 0) row "ready to type into", so a tutorial step that pressed Add row
   once per point left a fourth row, and entering the points from the top left that
   (0, 0) at the END of the line: LEM-1's embankment became (0,0) (20,20) (60,20)
   (0,0). Every other point table in Studio starts empty, and a pasted block grows
   the table to its own length either way.

2. **A profile line whose x values turn back is refused, with the line and the
   point named.** The zones are built from the lines left to right, so such a line
   silently builds a different section: LEM-3's section with a stray (0, 0) at the
   end of the foundation line searches to FS 1.007 against 1.245, and LEM-1's
   gives no surface at all. The Profile lines editor refuses OK on it, and the
   model checks' ``geometry.profile_line_turns_back`` refuses a run of a file that
   carries one (a file still opens with it: anything Studio can save reopens).
   A vertical step (equal x, as at a sheetpile) is not a turn back.

3. **A blank bottom elevation stays blank.** The Profile lines editor showed a
   model with no Max depth as 0 and wrote 0 back on OK — an elevation nobody
   entered, and a real change on an OK that changed nothing. Blank now reads,
   shows and saves as blank; a typed number is that number; text that is not a
   number is refused.

4. **The seepage boundary list names the value being typed.** A boundary's row
   ("Flux 1  (q = …)") was rebuilt only when a row was added or removed, so a value
   typed into the field showed in the list only after OK and reopening the editor.

5. **The Materials table heads Poisson's ratio ν**, as the list view does, not
   "n" — which is also the van Genuchten exponent's letter.

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

#: LEM-1's embankment (one line) and LEM-3's two layers.
LEM01 = os.path.join(_REPO, "docs/lem/files/xslope_simple_embankment.xlsx")
LEM03 = os.path.join(_REPO, "docs/lem/files/xslope_simple_mult_layers.xlsx")
#: SEEP-4's completed dam: one head boundary, three flux blocks, the exit face.
SEEP04 = os.path.join(_REPO, "docs/tutorials/files/xslope_dam_infiltration.xlsx")
#: RS2-51's wall: a polygon model with no bottom elevation.
NO_BOTTOM = os.path.join(_REPO, "docs/verification/files/rocscience/vp090_fem.xlsx")


def _fail(cond, message):
    return [] if cond else [message]


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _paste(table, rows):
    """Paste ``rows`` into an _EditableTable at its first cell with a real
    Ctrl/Cmd+V, as a reader pastes a table copied off the page."""
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QKeyEvent
    from PySide6.QtWidgets import QApplication
    QApplication.clipboard().setText("\n".join("\t".join(str(c) for c in r)
                                               for r in rows))
    table.table.setCurrentCell(0, 0)
    QApplication.sendEvent(table.table,
                           QKeyEvent(QEvent.KeyPress, Qt.Key_V, Qt.ControlModifier, "v"))


def _type_rows(table, rows):
    """Type ``rows`` into the table's rows from the top, cell by cell."""
    for i, r in enumerate(rows):
        for j, v in enumerate(r):
            table.table.item(i, j).setText(str(v))


@contextlib.contextmanager
def _warnings_caught():
    """Every QMessageBox.warning shown inside the block, as (title, text)."""
    from PySide6.QtWidgets import QMessageBox
    seen = []
    real = QMessageBox.warning
    QMessageBox.warning = staticmethod(
        lambda parent, title, text, *a, **k: seen.append((title, text))
        or QMessageBox.Ok)
    try:
        yield seen
    finally:
        QMessageBox.warning = real


def _accepted(dlg):
    from PySide6.QtWidgets import QDialog
    dlg.setResult(QDialog.Rejected)
    dlg.accept()
    return dlg.result() == QDialog.Accepted


# --------------------------------------------------------------------------- #
# 1. A new item's vertex table starts empty
# --------------------------------------------------------------------------- #
def test_new_line_starts_empty():
    from studio.editors import PolygonEditor, ProfileEditor
    out = []
    model = _load(LEM01)
    blanks = [("Profile lines", ProfileEditor(), dict(model, profile_lines=[])),
              ("Polygons", PolygonEditor(),
               dict(model, profile_lines=[], polygons=[], ssr_zones=[],
                    refine_zones=[], joint_zones=[]))]
    for name, editor, blank in blanks:
        dlg = editor.build(blank, None)
        dlg._add_line()
        n = dlg.table.table.rowCount()
        out += _fail(n == 0, f"{name}: a new item starts with {n} vertex row(s); "
                             f"Add row is then one short of a stray row")
        for _ in range(3):
            dlg.table._add_row()
        n = dlg.table.table.rowCount()
        out += _fail(n == 3, f"{name}: Add row three times on a new item gives {n} rows")
        dlg.deleteLater()

    # LEM-1's step as the 1.1.0 page first wrote it: Add line, Add row three times,
    # the three points typed from the top. The line is the three points, no more.
    dlg = ProfileEditor().build(dict(model, profile_lines=[]), None)
    dlg._add_line()
    for _ in range(3):
        dlg.table._add_row()
    pts = [(0, 0), (20, 20), (60, 20)]
    _type_rows(dlg.table, pts)
    got = dlg.result_lines()[0]["coords"]
    out += _fail([tuple(map(float, p)) for p in got] == [(0.0, 0.0), (20.0, 20.0),
                                                        (60.0, 20.0)],
                 f"LEM-1's line typed into three added rows reads {got}")
    dlg.deleteLater()

    # A pasted block brings its own rows into the empty table.
    dlg = ProfileEditor().build(dict(model, profile_lines=[]), None)
    dlg._add_line()
    _paste(dlg.table, pts)
    got = dlg.result_lines()[0]["coords"]
    out += _fail(len(got) == 3, f"LEM-1's three points pasted into a new line gave "
                                f"{len(got)} vertices: {got}")
    dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# 2. A profile line that turns back
# --------------------------------------------------------------------------- #
def test_turned_back_line_refused():
    from studio.editors import PolygonEditor, ProfileEditor
    out = []
    model = _load(LEM03)
    lines = [dict(ln, coords=list(ln["coords"])) for ln in model["profile_lines"]]
    lines[1]["coords"].append((0.0, 0.0))           # the stray row, at the end
    dlg = ProfileEditor().build(dict(model, profile_lines=lines), None)
    with _warnings_caught() as seen:
        ok = _accepted(dlg)
    out += _fail(not ok, "the Profile lines editor accepted a line that turns back "
                         "to (0, 0) at its end")
    text = seen[0][1] if seen else ""
    out += _fail("Line 2 turns back: point 3 (0, 0)" in text,
                 f"the refusal does not name the line and the point: {text!r}")
    out += _fail(dlg.list.currentRow() == 1,
                 f"the refusal left line {dlg.list.currentRow() + 1} selected, not "
                 f"the line it names")
    dlg.deleteLater()

    # Right to left is a turn back at its second point; a vertical step is not.
    for coords, want_ok, what in (
            ([(90.0, 20.0), (40.0, 20.0), (0.0, 0.0)], False, "entered right to left"),
            ([(0.0, 0.0), (40.0, 20.0), (40.0, 15.0), (90.0, 15.0)], True,
             "with a vertical step")):
        one = [dict(lines[0], coords=coords), dict(model["profile_lines"][1])]
        dlg = ProfileEditor().build(dict(model, profile_lines=one), None)
        with _warnings_caught():
            ok = _accepted(dlg)
        out += _fail(ok == want_ok, f"the Profile lines editor "
                                    f"{'refused' if want_ok else 'accepted'} a line "
                                    f"{what}")
        dlg.deleteLater()

    # A polygon ring runs either way round: the Polygons editor does not ask.
    poly = _load(os.path.join(_REPO, "docs/inputs/slope/xslope_simple1.xlsx"))
    dlg = PolygonEditor().build(dict(poly, profile_lines=[], ssr_zones=[],
                                     refine_zones=[], joint_zones=[]), None)
    with _warnings_caught():
        ok = _accepted(dlg)
    out += _fail(ok, "the Polygons editor refused a model's own rings")
    dlg.deleteLater()
    return out


def test_turned_back_line_in_the_model_checks():
    import copy
    from xslope.preflight import preflight, rules
    out = []
    rid = "geometry.profile_line_turns_back"
    if rid not in {r.id for r in rules()}:
        return [f"the model checks have no {rid} rule"]
    model = _load(LEM03)
    bad = copy.deepcopy(model)
    bad["profile_lines"][1]["coords"].append((0.0, 0.0))
    for analysis in ("lem", "seep", "fem"):
        rep = preflight(bad, analysis, ids=[rid])
        msgs = [f.message for f in rep.errors if f.rule_id == rid]
        out += _fail(msgs and "Profile line 2 turns back: point 3 (0, 0) lies to the "
                     "left of point 2 (90, 0)" in msgs[0],
                     f"{analysis}: a stray (0, 0) at the end of profile line 2 is "
                     f"reported as {msgs}")
    for name, sd in (("LEM-3 as shipped", model),
                     ("a vertical step", dict(copy.deepcopy(model), profile_lines=[
                         {"coords": [(0.0, 0.0), (40.0, 20.0), (40.0, 15.0),
                                     (90.0, 15.0)], "mat_id": 0, "size": None},
                         copy.deepcopy(model["profile_lines"][1])]))):
        rep = preflight(sd, "lem", ids=[rid])
        out += _fail(not [f for f in rep.findings if f.rule_id == rid],
                     f"{name}: the turn-back rule fired")
    return out


# --------------------------------------------------------------------------- #
# 3. A blank bottom elevation
# --------------------------------------------------------------------------- #
def test_blank_bottom_stays_blank():
    import copy
    from studio.editors import ProfileEditor
    out = []
    editor = ProfileEditor()

    # A profile model with no bottom yet (profile!B2 empty reads as None).
    model = _load(LEM01)
    model["max_depth"] = None
    dlg = editor.build(model, None)
    shown = dlg._max_depth_edit.text() if dlg._max_depth_edit is not None else None
    out += _fail(shown == "", f"a model with no bottom shows Max depth as {shown!r}")
    sd = copy.deepcopy(model)
    try:
        editor.apply(sd, dlg)
        out += _fail(sd.get("max_depth") is None,
                     f"OK on a model with no bottom wrote Max depth "
                     f"{sd.get('max_depth')!r}")
        out += _fail(not sd.get("polygons"),
                     "a model with no bottom got zones on OK (the loader builds none)")
    except Exception as e:                            # noqa: BLE001
        out.append(f"OK on a model with no bottom raised {type(e).__name__}: {e}")
    dlg.deleteLater()

    # RS2-51's polygon wall carries no bottom either: unchanged, it stays None.
    wall = _load(NO_BOTTOM)
    dlg = editor.build(wall, None)
    out += _fail(dlg.result_max_depth() is None,
                 f"the wall's blank bottom reads back as {dlg.result_max_depth()!r}")
    dlg.deleteLater()

    # Typed: a number is that number, a cleared field is blank, text is refused.
    model = _load(LEM01)
    for text, want in (("-5", -5.0), ("", None)):
        dlg = editor.build(model, None)
        dlg._max_depth_edit.setText(text)
        sd = copy.deepcopy(model)
        editor.apply(sd, dlg)
        out += _fail(sd.get("max_depth") == want,
                     f"Max depth typed {text!r} saved {sd.get('max_depth')!r}")
        dlg.deleteLater()
    dlg = editor.build(model, None)
    dlg._max_depth_edit.setText("ten")
    with _warnings_caught() as seen:
        ok = _accepted(dlg)
    out += _fail(not ok and seen, "Max depth 'ten' was accepted")
    dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# 4. The seepage boundary list
# --------------------------------------------------------------------------- #
def test_boundary_list_names_the_value():
    from studio.editors import SeepBcEditor
    out = []
    model = _load(SEEP04)
    dlg = SeepBcEditor().build(model, None)
    w = dlg._sets[0]
    n_heads, n_flux = len(w._heads), len(w._fluxes)
    w._add_flux()
    row = n_heads + n_flux
    w.flux_edit.setText("1e-08")
    label = w.list.item(row).text()
    out += _fail(label == f"Flux {n_flux + 1}  (q = 1e-08)",
                 f"a new flux block's row reads {label!r} while 1e-08 is typed")
    w.list.setCurrentRow(w.list.count() - 1)            # Exit face
    label = w.list.item(row).text()
    out += _fail(label == f"Flux {n_flux + 1}  (q = 1e-08)",
                 f"after moving to the exit face the new block's row reads {label!r}")
    # SEEP-4's binding step: a rate replaced by a series name.
    w.list.setCurrentRow(n_heads)
    w.flux_edit.setText("storm_face")
    label = w.list.item(n_heads).text()
    out += _fail(label == "Flux 1  (q = storm_face)",
                 f"Flux 1 bound to storm_face reads {label!r}")
    # A head on set 1 and a cleared value.
    w.list.setCurrentRow(0)
    w.head_edit.setText("12.5")
    out += _fail(w.list.item(0).text() == "Head 1  (h = 12.5)",
                 f"Head 1 typed 12.5 reads {w.list.item(0).text()!r}")
    w.head_edit.setText("")
    out += _fail(w.list.item(0).text() == "Head 1  (h = blank)",
                 f"Head 1 cleared reads {w.list.item(0).text()!r}")
    # What OK stores is what the list said.
    w.head_edit.setText("10")
    res = w.result()
    out += _fail(res["specified_fluxes"][0]["flux"] == "storm_face"
                 and res["specified_fluxes"][-1]["flux"] == 1e-08
                 and res["specified_heads"][0]["head"] == 10.0,
                 f"the set saved {res['specified_heads'][0]['head']!r}, "
                 f"{[f['flux'] for f in res['specified_fluxes']]}")
    # Set 2 is constant: a word typed there is 0, and the list says so.
    w2 = dlg._sets[1]
    w2._add_head()
    w2.head_edit.setText("pool")
    out += _fail(w2.list.item(0).text() == "Head 1  (h = 0)",
                 f"a word typed on set 2 reads {w2.list.item(0).text()!r}; it saves 0")
    dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# 5. Poisson's ratio
# --------------------------------------------------------------------------- #
def test_poissons_ratio_header():
    from studio.editors import MaterialsEditor
    model = _load(LEM01)
    dlg = MaterialsEditor().build(model, None)
    t = dlg._table
    j = [f.key for f in t._fields].index("nu")
    head = t.table.horizontalHeaderItem(j).text()
    dlg.deleteLater()
    return _fail(head == "ν", f"the Materials table heads Poisson's ratio {head!r}")


CHECKS = [
    ("a new profile line or polygon starts empty", test_new_line_starts_empty),
    ("a profile line that turns back is refused", test_turned_back_line_refused),
    ("the model checks name a turned-back line", test_turned_back_line_in_the_model_checks),
    ("a blank bottom elevation stays blank", test_blank_bottom_stays_blank),
    ("the boundary list names the typed value", test_boundary_list_names_the_value),
    ("the Materials table heads Poisson's ratio", test_poissons_ratio_header),
]


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
    print("Studio model editors:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll model-editor checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
