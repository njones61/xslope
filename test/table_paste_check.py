"""Standing checks on the table clipboard shared by every Studio editor -- what a
pasted block fills, what it refuses to fill, and whether the tutorials' own tables
are blocks this table can actually take.

Every editable input in Studio is a table of records over the same widget, and the
numbers in those tables nearly always exist somewhere else first: a worksheet, a
report, one of the tutorial pages, whose tables are laid out column for column as
their destination and which tell the reader to copy a block and paste it rather
than retype it. That sentence is only true while the paste is, so this file checks
both halves of it.

  A. THE MECHANICS. A block anchored at the current cell fills right and down;
     rows GROW to fit it (a table that would not lengthen could only be pasted into
     after the rows had been added by hand, which is most of the typing the paste
     removes -- and all of it when the table starts empty); columns do NOT (the
     columns are the fields the input has, so a wider block is a block from
     somewhere else, and its extra columns are dropped and counted rather than
     shifted onto fields they do not name). Values land the way typed ones do,
     which is what makes a pasted rendering of a stored number keep the stored
     number. A choice column takes any spelling of one of its own entries and
     nothing else. A cell the row's own state holds read-only is left alone. Each
     of those refusals is COUNTED and said out loud, because a paste that quietly
     filled 22 of 24 cells is worse than one that filled none.
  B. THE TUTORIALS' TABLES. Every taught table named below is read out of the page
     itself, pasted into the real editor for that input, and the resulting records
     compared against the published completed model for that tutorial -- to the
     precision the page prints. A page whose table gains a column, loses a column,
     or drifts from the model it builds fails here, on the page, rather than in a
     reader's model.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import itertools
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_TUTORIALS = os.path.join(_REPO, "docs", "tutorials")
_MODELS = os.path.join(_REPO, "docs", "lem", "files")
#: LEM-9's completed model — the tieback wall lives with the verification corpus it
#: comes from (Rocscience VP49), not in docs/lem/files.
LEM09_MODEL = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                           "vp049.xlsx")


def _fail(cond, message):
    return [] if cond else [message]


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


# --------------------------------------------------------------------------- #
# Driving the clipboard
# --------------------------------------------------------------------------- #
def _tsv(rows):
    return "\n".join("\t".join(str(c) for c in r) for r in rows)


def _paste(table, text, row=0, col=0):
    """Put ``text`` on the clipboard and paste it into ``table`` at (row, col) with
    a real Ctrl/Cmd+V key event -- through the key handling a user's keystroke goes
    through, not by calling the paste method."""
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QKeyEvent
    from PySide6.QtWidgets import QApplication

    QApplication.clipboard().setText(text)
    table.table.setCurrentCell(row, col)
    QApplication.sendEvent(
        table.table,
        QKeyEvent(QEvent.KeyPress, Qt.Key_V, Qt.ControlModifier, "v"))


def _copy(table, rows):
    """Select ``rows`` and copy them with a real Ctrl/Cmd+C key event; returns the
    clipboard text."""
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QKeyEvent
    from PySide6.QtWidgets import QApplication

    QApplication.clipboard().setText("")
    table.table.clearSelection()
    for r in rows:
        table.table.selectRow(r)
        table.table.setSelectionMode(table.table.SelectionMode.MultiSelection)
    table.table.setSelectionMode(table.table.SelectionMode.ExtendedSelection)
    QApplication.sendEvent(
        table.table,
        QKeyEvent(QEvent.KeyPress, Qt.Key_C, Qt.ControlModifier, "c"))
    return QApplication.clipboard().text()


def _summary(table):
    return table.paste_summary.text()


def _new_table(fields, rows=(), new_row=None, **kw):
    from studio.editors import _EditableTable
    if new_row is None:
        def new_row():
            return {f.key: f.default for f in fields}
    return _EditableTable(fields, [dict(r) for r in rows], new_row, **kw)


# --------------------------------------------------------------------------- #
# A. The mechanics
# --------------------------------------------------------------------------- #
def test_paste_fills_an_empty_table_and_grows_rows():
    """The taught flow: open an empty editor, paste the page's block, done. No Add
    row, no Generate — the rows arrive with the values."""
    from studio.editors import PiezoEditor

    t = _new_table(PiezoEditor.FIELDS)
    block = [(0, 80), (75, 79), (112, 76), (140, 70),
             (170, 61), (189.5, 52), (204.3, 40), (320, 40)]
    _paste(t, _tsv(block) + "\n")          # trailing newline: what a sheet copies
    rows = t.result_rows()
    out = _fail(len(rows) == 8, f"pasting 8 rows into an empty table left {len(rows)}")
    if len(rows) == 8:
        got = [(r["x"], r["y"]) for r in rows]
        out += _fail(got == [(float(x), float(y)) for x, y in block],
                     f"the pasted points came back as {got}")
    out += _fail(_summary(t) == "Pasted 8 rows × 2 columns.",
                 f"the status line read {_summary(t)!r}")
    return out


def test_paste_lands_inside_an_existing_table():
    """A block pasted at a cell partway down fills from there and grows only by what
    it needs — the rows above and the columns left of the anchor are untouched."""
    from studio.editors import PiezoEditor

    t = _new_table(PiezoEditor.FIELDS,
                   [{"x": 1.0, "y": 2.0}, {"x": 3.0, "y": 4.0}])
    _paste(t, _tsv([(9, 9), (10, 10)]), row=1, col=0)
    rows = [(r["x"], r["y"]) for r in t.result_rows()]
    return _fail(rows == [(1.0, 2.0), (9.0, 9.0), (10.0, 10.0)],
                 f"a paste at row 2 of a 2-row table gave {rows}")


def test_choice_columns_match_without_regard_to_case():
    """LEM-5's Movement column, pasted as the page prints it and as a spreadsheet
    might have lower-cased it — plus one word that is not one of the three, which
    must leave its cell alone and be counted rather than blanking the setting."""
    from studio.editors import NonCircEditor

    t = _new_table(NonCircEditor.FIELDS)
    _paste(t, _tsv([(0, 0, "Free"), (10, -5.8, "horiz"), (30, -5.8, "HORIZ"),
                    (39, 10, "sideways")]))
    got = [r["Movement"] for r in t.result_rows()]
    out = _fail(got == ["Free", "Horiz", "Horiz", "Free"],
                f"the Movement column came back as {got}")
    out += _fail("1 cell skipped (1 not a listed choice)" in _summary(t),
                 f"an unmatched choice was not reported: {_summary(t)!r}")
    # The rest of that row still landed: one bad cell is not a bad row.
    last = t.result_rows()[3]
    out += _fail((last["X"], last["Y"]) == (39.0, 10.0),
                 f"the row carrying the unmatched choice lost its numbers: {last}")
    return out


def test_extra_columns_are_dropped_and_counted():
    """A block wider than the table: the columns past the last one are dropped, and
    the status line says how many rather than letting them land on other fields."""
    from studio.editors import PiezoEditor

    t = _new_table(PiezoEditor.FIELDS)
    _paste(t, _tsv([(1, 2, 3, 4), (5, 6, 7, 8)]))
    rows = [(r["x"], r["y"]) for r in t.result_rows()]
    out = _fail(rows == [(1.0, 2.0), (5.0, 6.0)],
                f"a 4-column block into a 2-column table gave {rows}")
    out += _fail(_summary(t)
                 == "Pasted 2 rows × 2 columns; 4 cells skipped (4 past the last column).",
                 f"the status line read {_summary(t)!r}")
    return out


def test_hidden_columns_are_not_filled():
    """A column folded away by a usage toggle is not a column the block fills: the
    block belongs to the columns the user can see."""
    from studio.editors import ReinforcementEditor

    t = _new_table(ReinforcementEditor.FIELDS, new_row=_new_reinf_row)
    t.apply_usage_filter({"lem"})          # the FEM tail (Tres, E, Area) folds away
    _paste(t, _tsv([["layer 1", 0, 0, 20, 0, "geosynthetic", "tangent", "active",
                     800, 4, 4, 0, 0, 1, 77, 77, 77]]))
    row = t.result_rows()[0]
    out = _fail(row["spacing"] == 1.0, f"the last visible column took {row['spacing']}")
    out += _fail((row["t_res"], row["E"], row["area"]) == (0.0, 0.0, 0.0),
                 "a hidden column was filled from the block")
    out += _fail("3 cells skipped (3 past the last column)" in _summary(t),
                 f"the status line read {_summary(t)!r}")
    return out


def test_read_only_cells_are_skipped():
    """The matric-suction pair is inert on a material that reads no signed pore
    pressure, and the editor grays it out. A block pasted across it leaves those
    cells as they were and says so."""
    from studio.editors import MaterialsEditor, _mat_dim_keys

    fields = MaterialsEditor.FIELDS
    keys = [f.key for f in fields]
    row = {f.key: f.default for f in fields}
    row.update(name="clay", option="mc", u="none", ru=None, phi_b=None, s_cap=None)
    t = _new_table(fields, [row], dim_rule=_mat_dim_keys)
    # Anchored on the live cell to their left and running across them — the only way
    # a block reaches a grayed cell, since a grayed cell cannot be made the current
    # one to anchor on.
    _paste(t, _tsv([[0.25, 15, 300]]), row=0, col=keys.index("ru"))
    got = t.result_rows()[0]
    out = _fail(got["ru"] == 0.25, f"the live cell of the block took {got['ru']!r}")
    out += _fail((got["phi_b"], got["s_cap"]) == (None, None),
                 f"a grayed cell took a pasted value: {got['phi_b']}, {got['s_cap']}")
    out += _fail("2 cells skipped (2 read-only)" in _summary(t),
                 f"the status line read {_summary(t)!r}")
    return out


def test_a_pasted_rounding_keeps_the_stored_value():
    """A pasted cell goes in the way a typed one does, which is the whole of this:
    a circle's R is stored at full precision and SHOWN rounded, so pasting back the
    number the cell shows must keep the stored value — exactly as retyping it
    does — instead of moving the model by the rounding."""
    from studio.editors import CirclesEditor

    stored = 41.23105625617661
    t = _new_table(CirclesEditor.FIELDS,
                   [{"Xo": 10.0, "Yo": 40.0, "Option": "Intercept", "Depth": 0.0,
                     "Xi": 0.0, "Yi": 0.0, "R": stored}])
    shown = t.table.item(0, 6).text()
    _paste(t, shown, row=0, col=6)
    got = t.result_rows()[0]["R"]
    return (_fail(shown == "41.2311", f"R was displayed as {shown!r}")
            + _fail(got == stored, f"pasting the displayed R stored {got!r}"))


def test_paste_notifies_once():
    """A paste is one edit: the preview redraws for the finished table, not for
    every cell of a half-filled one — and it does redraw, which is the parity with
    a hand edit that the live previews depend on."""
    from studio.editors import PiezoEditor

    fired = []
    t = _new_table(PiezoEditor.FIELDS, on_change=lambda: fired.append(1))
    _paste(t, _tsv([(1, 2), (3, 4), (5, 6)]))
    return _fail(len(fired) == 1,
                 f"a 3-row paste fired {len(fired)} change notifications, not 1")


def test_copy_round_trips():
    """Copy is the natural other half: a selection out as the same tab-separated
    block, which this table's own paste reads back into identical rows."""
    from studio.editors import NonCircEditor

    t = _new_table(NonCircEditor.FIELDS,
                   [{"X": 0.0, "Y": 0.0, "Movement": "Free"},
                    {"X": 10.5, "Y": -5.8, "Movement": "Horiz"}])
    before = t.result_rows()
    text = _copy(t, [0, 1])
    out = _fail(text == "0\t0\tFree\n10.5\t-5.8\tHoriz",
                f"the copied block was {text!r}")
    # Anchored on the last row: the block runs from there and the table grows under
    # it, which is how a copied block is appended to a table already holding rows.
    _paste(t, text, row=1, col=0)
    rows = t.result_rows()
    out += _fail(len(rows) == 3, f"pasting a copied 2-row block gave {len(rows)} rows")
    if len(rows) == 3:
        out += _fail(rows[1:] == before,
                     f"the round trip changed the rows: {rows}")
    return out


def test_cancel_discards_a_paste():
    """The dialog's own OK/Cancel boundary governs a paste like any other edit: the
    document is written by the editor's apply(), and Studio runs apply() only when
    the dialog was accepted."""
    from studio.editors import CirclesEditor

    model = _load(os.path.join(_MODELS, "xslope_simple_mult_layers.xlsx"))
    before = [dict(c) for c in model["circles"]]
    editor = CirclesEditor()
    dlg = editor.build(model, None)
    _paste(dlg._editable, _tsv([(99, 99, "Depth", 9)]))
    dlg.reject()
    out = _fail(model["circles"] == before,
                "a cancelled dialog left the paste in the model")
    # ...and the same paste, accepted, does land — or the check above proves nothing.
    dlg2 = editor.build(model, None)
    _paste(dlg2._editable, _tsv([(99, 99, "Depth", 9)]))
    dlg2.accept()
    editor.apply(model, dlg2)
    out += _fail(model["circles"][0]["Xo"] == 99.0,
                 "an accepted paste did not reach the model")
    return out


def _new_reinf_row():
    from studio.editors import _new_reinf
    return _new_reinf()


def test_the_reinforcement_two_block_paste():
    """LEM-8's reinforcement table goes in as TWO blocks, because the sheet's Dir and
    Appl sit between them and are filled from Type. The first block starts at Label —
    the editor's first column, as it is the sheet's — and the capacity values anchor
    at Tmax; the capacity block only lands contiguously while the editor's columns
    stay in the worksheet's order."""
    from studio.editors import ReinforcementEditor

    t = _new_table(ReinforcementEditor.FIELDS, new_row=_new_reinf_row)
    ends = [("layer 1", 0, 0, 20, 0), ("layer 2", 5, 4, 25, 4),
            ("layer 3", 10, 8, 30, 8), ("layer 4", 15, 12, 35, 12),
            ("layer 5", 20, 16, 40, 16), ("layer 6", 25, 20, 45, 20)]
    _paste(t, _tsv(ends))
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    _paste(t, _tsv([(800, 4, 4, 0, 0, 1)] * 6), row=0, col=keys.index("t_max"))
    rows = t.result_rows()
    out = _fail(len(rows) == 6, f"the two-block paste left {len(rows)} lines")
    out += _fail(_summary(t) == "Pasted 6 rows × 6 columns.",
                 f"the second block reported {_summary(t)!r}")
    for i, r in enumerate(rows[:6]):
        out += _fail(r["label"] == ends[i][0],
                     f"line {i + 1} label came back as {r['label']!r}")
        got = (r["x1"], r["y1"], r["x2"], r["y2"])
        out += _fail(got == tuple(float(v) for v in ends[i][1:]),
                     f"line {i + 1} endpoints came back as {got}")
        cap = (r["t_max"], r["lp1"], r["lp2"], r["tend1"], r["tend2"], r["spacing"])
        out += _fail(cap == (800.0, 4.0, 4.0, 0.0, 0.0, 1.0),
                     f"line {i + 1} capacities came back as {cap}")
    return out


def test_the_support_editors_columns_are_the_worksheets():
    """The reinforcement and piles editors' columns, in order, against the sheets
    they mirror — read out of the template rather than restated here.

    This is what makes every block on this page's tutorials land: the tables are
    printed in the worksheet's column order and pasted into the editor, so the two
    orders are one fact. A column added, dropped or moved on either side shifts a
    reader's values one field to the left with nothing on screen saying so.

    Only ``#`` (the sheet's row number) and the piles sheet's ``qp`` are dropped: qp
    is the pile force angle, which the editor derives from the pile's own axis on
    save, and which is why LEM-9's pile goes in as two blocks."""
    import openpyxl

    from studio.editors import PilesEditor, ReinforcementEditor

    template = os.path.join(_REPO, "docs", "inputs", "input_template.xlsx")
    wb = openpyxl.load_workbook(template, read_only=True)
    out = []
    for sheet, editor, derived in (("reinforce", ReinforcementEditor, ()),
                                   ("piles", PilesEditor, ("qp",))):
        ws = wb[sheet]
        # The header run from column A, stopping at the first empty cell: both sheets
        # keep their drop-down source lists in the same row, far off to the right.
        headers = [str(c.value).strip() for c in itertools.takewhile(
            lambda c: c.value is not None,
            next(ws.iter_rows(min_row=2, max_row=2)))]
        want = [h for h in headers if h != "#" and h not in derived]
        got = [f.header for f in editor.FIELDS]
        out += _fail(got == want,
                     f"the {editor.label} editor's columns are {got}, the {sheet!r} "
                     f"sheet's are {want}")
    wb.close()
    return out


# --------------------------------------------------------------------------- #
# A2. The Type preset — the sheet's Dir/Appl formulas, in the editor
#
# reinforce!H and !I are VLOOKUPs on the Type column, so on the worksheet a support
# type is one pick that answers both. The editor has to do the same thing, or the
# Type a reader picks in Studio is a label with no physics behind it: an Anchor
# whose Dir stayed 'tangent' is solved as a geosynthetic, at a different force
# direction, with nothing on screen saying so.
# --------------------------------------------------------------------------- #
def _preset_combos(t):
    """The (type, dir, appl) combo widgets of a reinforcement table's first row."""
    from studio.editors import ReinforcementEditor
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    return tuple(t.table.cellWidget(0, keys.index(k))
                 for k in ("type", "dir", "appl"))


def test_type_preset_fills_dir_and_appl():
    """Picking a Type fills Dir and Appl from the preset table; changing either
    afterwards is an override that STAYS (nothing rewrites it); picking the same Type
    again re-asserts the preset over that override, which is what the sheet's formula
    does to a cell that has not been overtyped."""
    from studio.editors import ReinforcementEditor

    t = _new_table(ReinforcementEditor.FIELDS, new_row=_new_reinf_row,
                   preset_spec=ReinforcementEditor.PRESET_SPEC)
    t._add_row()
    typ, direction, appl = _preset_combos(t)
    typ.setCurrentText("anchor")
    out = _fail((direction.currentText(), appl.currentText()) == ("axial", "active"),
                f"picking anchor left Dir/Appl at "
                f"({direction.currentText()!r}, {appl.currentText()!r})")
    # nail is the preset that differs in BOTH columns, so a fill that only ever wrote
    # one of them cannot pass this.
    typ.setCurrentText("nail")
    out += _fail((direction.currentText(), appl.currentText()) == ("axial", "passive"),
                 f"picking nail left Dir/Appl at "
                 f"({direction.currentText()!r}, {appl.currentText()!r})")
    # An override survives: the user's own answer is not overwritten while the Type
    # stands.
    direction.setCurrentText("tangent")
    row = t.result_rows()[0]
    out += _fail((row["dir"], row["appl"]) == ("tangent", "passive"),
                 f"an overridden Dir came back as {row['dir']!r}/{row['appl']!r}")
    # ...until the Type is picked again. Emitting `activated` is what a user's pick of
    # the entry already showing does; it changes no index, so only that signal carries
    # it.
    typ.activated.emit(typ.currentIndex())
    row = t.result_rows()[0]
    out += _fail((row["dir"], row["appl"]) == ("axial", "passive"),
                 f"re-picking nail left Dir/Appl at {row['dir']!r}/{row['appl']!r}")
    # A blank Type fills nothing — the sheet leaves those cells empty for a generic
    # line, and the values standing in them are the loader's defaults for one.
    typ.setCurrentText("")
    row = t.result_rows()[0]
    out += _fail((row["dir"], row["appl"]) == ("axial", "passive"),
                 f"a blank Type rewrote Dir/Appl to {row['dir']!r}/{row['appl']!r}")
    return out


def test_a_pasted_type_fills_dir_and_appl():
    """A pasted Type fills exactly as a picked one does — the paths a reader actually
    uses are the drop-down and the block, and the tutorials teach the block. A block
    that carries Dir or Appl itself still wins, because those columns are written
    after the Type that filled them."""
    from studio.editors import ReinforcementEditor

    keys = [f.key for f in ReinforcementEditor.FIELDS]
    # Two lines already entered (the endpoints block is what puts them there), then
    # the Type column pasted down them -- LEM-8's and LEM-9's Studio flow.
    t = _new_table(ReinforcementEditor.FIELDS,
                   [_new_reinf_row(), _new_reinf_row()], new_row=_new_reinf_row,
                   preset_spec=ReinforcementEditor.PRESET_SPEC)
    _paste(t, _tsv([["Anchor"], ["Nail"]]), row=0, col=keys.index("type"))
    rows = t.result_rows()
    out = _fail(len(rows) == 2, f"a 2-row Type block left {len(rows)} rows")
    out += _fail([(r["type"], r["dir"], r["appl"]) for r in rows]
                 == [("anchor", "axial", "active"), ("nail", "axial", "passive")],
                 f"the pasted Types came back as "
                 f"{[(r['type'], r['dir'], r['appl']) for r in rows]}")
    # Pasting the SAME type over a row that already holds it re-asserts the preset:
    # the combo's index does not change, so this only works if the paste fills
    # explicitly rather than riding the index-changed signal.
    t.table.cellWidget(0, keys.index("dir")).setCurrentText("tangent")
    _paste(t, _tsv([["Anchor"]]), row=0, col=keys.index("type"))
    out += _fail(t.result_rows()[0]["dir"] == "axial",
                 f"re-pasting the same Type left Dir at "
                 f"{t.result_rows()[0]['dir']!r}")
    # A block spanning Type, Dir and Appl: the block's own values are the answer.
    t2 = _new_table(ReinforcementEditor.FIELDS, [_new_reinf_row()],
                    new_row=_new_reinf_row,
                    preset_spec=ReinforcementEditor.PRESET_SPEC)
    _paste(t2, _tsv([["Anchor", "Tangent", "Passive"]]), row=0, col=keys.index("type"))
    got = t2.result_rows()[0]
    out += _fail((got["dir"], got["appl"]) == ("tangent", "passive"),
                 f"a block carrying Dir/Appl came back as "
                 f"{got['dir']!r}/{got['appl']!r} — the preset overwrote it")
    return out


def test_the_list_view_fills_dir_and_appl_too():
    """The same fill in the list view, where a reader entering one anchor at a time
    meets it: the two views edit the same rows, so a Type picked in either has to mean
    the same thing."""
    from studio.editors import ReinforcementEditor

    model = _load(os.path.join(_MODELS, "xslope_reinforce.xlsx"))
    editor = ReinforcementEditor()
    dlg = editor.build(dict(model, reinforcement_lines=[_new_reinf_row()]), None)
    dlg.set_view_mode("list")
    edits = dlg._list_view._edits
    edits["type"].setCurrentText("anchor")
    out = _fail((edits["dir"].currentText(), edits["appl"].currentText())
                == ("axial", "active"),
                f"picking anchor in the list view left Dir/Appl at "
                f"({edits['dir'].currentText()!r}, {edits['appl'].currentText()!r})")
    landed = dict(model, reinforcement_lines=[])
    dlg.accept()
    editor.apply(landed, dlg)
    row = landed["reinforcement_lines"][0]
    out += _fail((row["type"], row["dir"], row["appl"])
                 == ("anchor", "axial", "active"),
                 f"the list view's Anchor reached the model as "
                 f"{(row['type'], row['dir'], row['appl'])}")
    # A line loaded with an override is NOT re-filled by being looked at: selecting it
    # must not rewrite what the file says.
    over = dict(_new_reinf_row(), type="anchor", dir="tangent", appl="passive")
    dlg2 = editor.build(dict(model, reinforcement_lines=[over]), None)
    dlg2.set_view_mode("list")
    landed2 = dict(model, reinforcement_lines=[])
    dlg2.accept()
    editor.apply(landed2, dlg2)
    got = landed2["reinforcement_lines"][0]
    out += _fail((got["dir"], got["appl"]) == ("tangent", "passive"),
                 f"opening the editor on an overridden line rewrote it to "
                 f"{got['dir']!r}/{got['appl']!r}")
    return out


# --------------------------------------------------------------------------- #
# B. The tutorials' own tables
# --------------------------------------------------------------------------- #
_RULE = re.compile(r":?-{2,}:?")


def _page_tables(page):
    """Every markdown table on a tutorial page, as lists of rows of cell strings —
    the first row being the headers.

    Cells are read as the RENDERED page hands them to the clipboard: a reader copies
    from the page, not from its source, so ```none``` arrives as ``none``
    and ``**1.247**`` as ``1.247``."""
    out, cur = [], None
    with open(os.path.join(_TUTORIALS, page), encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("|") and s.endswith("|"):
                cells = [c.strip().strip("`").replace("**", "")
                         for c in s.strip("|").split("|")]
                if all(_RULE.fullmatch(c) for c in cells):
                    continue                       # the |---|---| separator row
                cur = (cur or []) + [cells]
            else:
                if cur:
                    out.append(cur)
                cur = None
    if cur:
        out.append(cur)
    return out


def _taught(page, headers, occurrence=0):
    """The ``occurrence``-th table on ``page`` whose headers are exactly ``headers``,
    as (rows, headers). Raises when the page no longer carries it — a page that
    renamed or reshaped a taught table is the failure, not a skip."""
    found = [t for t in _page_tables(page) if t[0] == list(headers)]
    if len(found) <= occurrence:
        raise AssertionError(
            f"{page}: no table #{occurrence + 1} with headers {list(headers)} "
            f"(found {len(found)})")
    return found[occurrence][1:]


def _matches(pasted, printed, published):
    """A pasted number against the published model's own value.

    The page prints what a reader can reasonably type; the model file may carry more
    digits than that (a generated vertex at 10.682271139 is taught as 10.6823). So
    the tolerance is half a unit in the last place the PAGE prints — a page that
    drops a digit or names the wrong column fails, a page that rounds does not."""
    if isinstance(published, str) or isinstance(pasted, str):
        return str(pasted).lower() == str(published).lower()
    dec = len(printed.split(".")[1]) if "." in printed else 0
    return abs(float(pasted) - float(published)) <= 0.5 * 10 ** (-dec) + 1e-9


def _compare(label, pasted_rows, printed_rows, published_rows, keys):
    """Field for field: what the page's block produced against the model the page
    ships, on the keys the block carries.

    A key the loaded model does not carry (a circle's Option, which the loader
    consumes into R rather than storing) is compared against the page's own cell —
    the block still has to have landed in that column."""
    out = _fail(len(pasted_rows) == len(published_rows),
                f"{label}: the page's table built {len(pasted_rows)} records, "
                f"the completed model has {len(published_rows)}")
    for i, (got, printed, want) in enumerate(
            zip(pasted_rows, printed_rows, published_rows)):
        for k, cell in zip(keys, printed):
            ref = want[k] if k in want else cell
            out += _fail(_matches(got[k], cell, ref),
                         f"{label}: row {i + 1} {k} pasted as {got[k]!r} from "
                         f"{cell!r}; expected {ref!r}")
    return out


def _xy_lines(dlg, blocks):
    """Paste one x/y block per line/polygon into a MatGeometryDialog, adding an item
    for each — the profile-line and polygon editors' own flow."""
    for rows in blocks:
        dlg._add_line()
        _paste(dlg.table, _tsv(rows))
    return dlg.result_lines()


def test_lem03_profile_lines_and_circles():
    """LEM-3's two profile lines and its two starting circles."""
    from studio.editors import CirclesEditor, ProfileEditor

    page = "lem03_layered_slope.md"
    model = _load(os.path.join(_MODELS, "xslope_simple_mult_layers.xlsx"))
    line1 = _taught(page, ["x (ft)", "y (ft)"], 0)
    line2 = _taught(page, ["x (ft)", "y (ft)"], 1)

    blank = dict(model, profile_lines=[])
    dlg = ProfileEditor().build(blank, None)
    lines = _xy_lines(dlg, [line1, line2])
    out = _fail(len(lines) == 2, f"LEM-3 profile lines: pasted {len(lines)} lines")
    for n, (got, printed, want) in enumerate(
            zip(lines, [line1, line2], model["profile_lines"])):
        out += _fail(len(got["coords"]) == len(want["coords"]),
                     f"LEM-3 line {n + 1}: {len(got['coords'])} vertices pasted, "
                     f"{len(want['coords'])} in the completed model")
        for (gx, gy), cell, (wx, wy) in zip(got["coords"], printed, want["coords"]):
            out += _fail(_matches(gx, cell[0], wx) and _matches(gy, cell[1], wy),
                         f"LEM-3 line {n + 1}: pasted ({gx}, {gy}) against "
                         f"({wx}, {wy})")

    circles = _taught(page, ["Xo", "Yo", "Option", "Depth"])
    editor = CirclesEditor()
    cdlg = editor.build(dict(model, circles=[]), None)
    _paste(cdlg._editable, _tsv(circles))
    landed = dict(model, circles=[])
    cdlg.accept()
    editor.apply(landed, cdlg)
    out += _compare("LEM-3 circles", landed["circles"], circles, model["circles"],
                    ["Xo", "Yo", "Option", "Depth"])
    # Option is not a stored key — what it does is set R, which is what the model
    # carries. A block that pasted "Depth" into the wrong column would not.
    for i, (got, want) in enumerate(zip(landed["circles"], model["circles"])):
        out += _fail(abs(got["R"] - want["R"]) < 1e-6,
                     f"LEM-3 circle {i + 1}: R came out {got['R']}, "
                     f"the model has {want['R']}")
    return out


def test_lem04_piezometric_line():
    """LEM-4's eight-point piezometric line — the page's own "the one most worth
    pasting" — into the piezo editor's first tab."""
    from studio.editors import PiezoEditor

    page = "lem04_water_in_the_slope.md"
    model = _load(os.path.join(_MODELS, "xslope_method_slices_problem.xlsx"))
    pts = _taught(page, ["x (ft)", "y (ft)"], 3)      # after the three profile lines
    editor = PiezoEditor()
    dlg = editor.build(dict(model, piezo_line=[], piezo_line2=[]), None)
    _paste(dlg._editables[0], _tsv(pts))
    landed = dict(model, piezo_line=[], piezo_line2=[])
    editor.apply(landed, dlg)
    got = landed["piezo_line"]
    want = model["piezo_line"]
    out = _fail(len(got) == len(want) == 8,
                f"LEM-4 piezo: {len(got)} points pasted, {len(want)} in the model")
    for i, ((gx, gy), cell, (wx, wy)) in enumerate(zip(got, pts, want)):
        out += _fail(_matches(gx, cell[0], wx) and _matches(gy, cell[1], wy),
                     f"LEM-4 piezo point {i + 1}: pasted ({gx}, {gy}) against "
                     f"({wx}, {wy})")
    return out


def test_lem05_noncircular_surface():
    """LEM-5's four-point surface, Movement column live: the taught block carries a
    choice column, and the two seam points are Horiz."""
    from studio.editors import NonCircEditor

    page = "lem05_weak_layer_noncircular.md"
    model = _load(os.path.join(_MODELS, "xslope_noncircular.xlsx"))
    pts = _taught(page, ["X (ft)", "Y (ft)", "Movement"])
    editor = NonCircEditor()
    dlg = editor.build(dict(model, non_circ=[]), None)
    _paste(dlg._editable, _tsv(pts))
    landed = dict(model, non_circ=[])
    editor.apply(landed, dlg)
    out = _compare("LEM-5 non-circular surface", landed["non_circ"], pts,
                   model["non_circ"], ["X", "Y", "Movement"])
    out += _fail(_summary(dlg._editable) == "Pasted 4 rows × 3 columns.",
                 f"LEM-5's block reported {_summary(dlg._editable)!r}")
    return out


def test_lem06_polygon_rings():
    """LEM-6's two polygon rings — the vertices listed once each, the ring closed by
    the editor."""
    from studio.editors import PolygonEditor

    page = "lem06_polygon_geometry.md"
    model = _load(os.path.join(_MODELS, "xslope_sloping_bottom.xlsx"))
    rings = [_taught(page, ["x (ft)", "y (ft)"], 0),
             _taught(page, ["x (ft)", "y (ft)"], 1)]
    editor = PolygonEditor()
    dlg = editor.build(dict(model, polygons=[], ssr_zones=[], refine_zones=[]), None)
    lines = _xy_lines(dlg, rings)
    out = _fail(len(lines) == 2, f"LEM-6: pasted {len(lines)} polygons")
    for n, (got, printed, want) in enumerate(zip(lines, rings, model["polygons"])):
        coords = list(want["polygon"].exterior.coords)[:-1]     # drop the closing dup
        out += _fail(len(got["coords"]) == len(coords),
                     f"LEM-6 polygon {n + 1}: {len(got['coords'])} vertices pasted, "
                     f"{len(coords)} in the completed model")
        for (gx, gy), cell, (wx, wy) in zip(got["coords"], printed, coords):
            out += _fail(_matches(gx, cell[0], wx) and _matches(gy, cell[1], wy),
                         f"LEM-6 polygon {n + 1}: pasted ({gx}, {gy}) against "
                         f"({wx}, {wy})")
    return out


def test_lem08_reinforcement_blocks():
    """LEM-8's two reinforcement blocks pasted into the real editor's table view: the
    Label..y2 table whole (Label is the editor's first column, as it is the sheet's),
    then Type set per line, then the capacity block at Tmax."""
    from studio.editors import ReinforcementEditor

    page = "lem08_reinforced_slope.md"
    model = _load(os.path.join(_MODELS, "xslope_reinforce.xlsx"))
    ends = _taught(page, ["Label", "x1", "y1", "x2", "y2"])
    caps = _taught(page, ["Tmax", "Lp1", "Lp2", "Tend1", "Tend2", "Spacing"])
    editor = ReinforcementEditor()
    dlg = editor.build(dict(model, reinforcement_lines=[]), None)
    dlg.set_view_mode("table")
    dlg._table.apply_usage_filter({"lem", "fem"})
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    _paste(dlg._table, _tsv(ends))                             # Label..y2, whole
    out = _fail(_summary(dlg._table) == "Pasted 6 rows × 5 columns.",
                f"LEM-8's endpoint block reported {_summary(dlg._table)!r}")
    _paste(dlg._table, _tsv([["Geosynthetic"]] * len(ends)),
           row=0, col=keys.index("type"))
    _paste(dlg._table, _tsv(caps), row=0, col=keys.index("t_max"))
    out += _fail(_summary(dlg._table) == "Pasted 6 rows × 6 columns.",
                 f"LEM-8's capacity block reported {_summary(dlg._table)!r}")
    landed = dict(model, reinforcement_lines=[])
    dlg.accept()
    editor.apply(landed, dlg)
    rows = landed["reinforcement_lines"]
    out += _compare("LEM-8 labels and endpoints", rows, ends,
                    model["reinforcement_lines"],
                    ["label", "x1", "y1", "x2", "y2"])
    out += _compare("LEM-8 capacities", rows, caps, model["reinforcement_lines"],
                    ["t_max", "lp1", "lp2", "tend1", "tend2", "spacing"])
    for i, (got, want) in enumerate(zip(rows, model["reinforcement_lines"])):
        for k in ("type", "dir", "appl"):
            out += _fail(str(got[k]).lower() == str(want[k]).lower(),
                         f"LEM-8 line {i + 1} {k}: {got[k]!r} against {want[k]!r}")
    return out


def test_lem09_anchor_blocks():
    """LEM-9's two anchors, built the way its Studio step now tells the reader to: the
    Label..y2 block, the Type set to Anchor, the capacity block at Tmax — and then
    every field measured against the model the page ships.

    This is the leg the preset is FOR. The page says picking Anchor answers Dir and
    Appl; vp049 carries dir='axial', appl='active' on both lines; and a reader who
    never types either has to land on exactly that, or the page teaches a model
    different from the one it publishes results for."""
    from studio.editors import ReinforcementEditor

    page = "lem09_tieback_wall.md"
    model = _load(LEM09_MODEL)
    ends = _taught(page, ["Label", "x1", "y1", "x2", "y2"])
    caps = _taught(page, ["Tmax", "Lp1", "Lp2", "Tend1", "Tend2", "Spacing"])
    editor = ReinforcementEditor()
    dlg = editor.build(dict(model, reinforcement_lines=[]), None)
    dlg.set_view_mode("table")
    dlg._table.apply_usage_filter({"lem", "fem"})
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    _paste(dlg._table, _tsv(ends))                             # Label..y2, whole
    _paste(dlg._table, _tsv([["Anchor"]] * len(ends)), row=0, col=keys.index("type"))
    _paste(dlg._table, _tsv(caps), row=0, col=keys.index("t_max"))
    out = _fail(_summary(dlg._table) == "Pasted 2 rows × 6 columns.",
                f"LEM-9's capacity block reported {_summary(dlg._table)!r}")
    landed = dict(model, reinforcement_lines=[])
    dlg.accept()
    editor.apply(landed, dlg)
    rows = landed["reinforcement_lines"]
    want = model["reinforcement_lines"]
    out += _compare("LEM-9 labels and endpoints", rows, ends, want,
                    ["label", "x1", "y1", "x2", "y2"])
    out += _compare("LEM-9 capacities", rows, caps, want,
                    ["t_max", "lp1", "lp2", "tend1", "tend2", "spacing"])
    for i, (got, w) in enumerate(zip(rows, want)):
        for k in ("type", "dir", "appl"):
            out += _fail(str(got[k]).lower() == str(w[k]).lower(),
                         f"LEM-9 anchor {i + 1} {k}: pasted {got[k]!r}, the model "
                         f"has {w[k]!r} — neither Dir nor Appl is typed by the reader")
    return out


def test_lem09_pile_row():
    """LEM-9's soldier pile, in two blocks either side of the sheet's qp and Appl —
    the columns the page leaves empty and the editor derives or defaults. The second
    block lands contiguously only while the piles editor's columns are the piles
    sheet's."""
    from studio.editors import PilesEditor

    page = "lem09_tieback_wall.md"
    model = _load(LEM09_MODEL)
    head = _taught(page, ["Label", "x1", "y1", "x2", "y2", "H"])
    size = _taught(page, ["D", "S"])
    editor = PilesEditor()
    dlg = editor.build(dict(model, pile_lines=[]), None)
    dlg.set_view_mode("table")
    dlg._table.apply_usage_filter({"lem", "fem"})
    keys = [f.key for f in PilesEditor.FIELDS]
    _paste(dlg._table, _tsv(head))
    _paste(dlg._table, _tsv(size), row=0, col=keys.index("D_pile"))
    out = _fail(_summary(dlg._table) == "Pasted 1 row × 2 columns.",
                f"LEM-9's pile size block reported {_summary(dlg._table)!r}")
    landed = dict(model, pile_lines=[])
    dlg.accept()
    editor.apply(landed, dlg)
    out += _compare("LEM-9 pile", landed["pile_lines"], head, model["pile_lines"],
                    ["label", "x1", "y1", "x2", "y2", "H"])
    out += _compare("LEM-9 pile size", landed["pile_lines"], size,
                    model["pile_lines"], ["D_pile", "S"])
    got = landed["pile_lines"][0]
    # The two columns the blocks step over: Appl defaults to active (a blank cell on
    # the sheet means the same), and θ is derived from the axis rather than pasted.
    out += _fail(got["appl"] == "active",
                 f"the pile's Appl came back {got['appl']!r}; the page leaves the "
                 f"column empty and a blank Appl is active")
    out += _fail(abs(got["theta_p"] - float(model["pile_lines"][0]["theta_p"])) < 1e-9,
                 f"the pile's θ came back {got['theta_p']!r}, the model has "
                 f"{model['pile_lines'][0]['theta_p']!r}")
    return out


CHECKS = [
    ("a block fills an empty table, growing rows",
     test_paste_fills_an_empty_table_and_grows_rows),
    ("a block lands at the anchor cell", test_paste_lands_inside_an_existing_table),
    ("choice columns match case-insensitively",
     test_choice_columns_match_without_regard_to_case),
    ("extra columns are dropped and counted",
     test_extra_columns_are_dropped_and_counted),
    ("hidden columns are not filled", test_hidden_columns_are_not_filled),
    ("read-only cells are skipped", test_read_only_cells_are_skipped),
    ("a pasted rounding keeps the stored value",
     test_a_pasted_rounding_keeps_the_stored_value),
    ("a paste is one edit", test_paste_notifies_once),
    ("copy round-trips through paste", test_copy_round_trips),
    ("cancel discards a paste", test_cancel_discards_a_paste),
    ("the reinforcement two-block paste", test_the_reinforcement_two_block_paste),
    ("the support editors' columns are the sheets'",
     test_the_support_editors_columns_are_the_worksheets),
    ("a picked Type fills Dir and Appl", test_type_preset_fills_dir_and_appl),
    ("a pasted Type fills Dir and Appl", test_a_pasted_type_fills_dir_and_appl),
    ("the list view fills Dir and Appl too",
     test_the_list_view_fills_dir_and_appl_too),
    ("LEM-3 profile lines and circles", test_lem03_profile_lines_and_circles),
    ("LEM-4 piezometric line", test_lem04_piezometric_line),
    ("LEM-5 non-circular surface", test_lem05_noncircular_surface),
    ("LEM-6 polygon rings", test_lem06_polygon_rings),
    ("LEM-8 reinforcement blocks", test_lem08_reinforcement_blocks),
    ("LEM-9 anchor blocks", test_lem09_anchor_blocks),
    ("LEM-9 pile row", test_lem09_pile_row),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("table paste: PySide6 not installed — skipped.")
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
    print("table paste:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll table paste checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
