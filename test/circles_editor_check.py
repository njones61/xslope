"""Standing checks on the Studio Circles editor -- what it shows, what it saves, and
what its Generate button does with work already in the table.

Three ways this dialog can be wrong, all of them quiet:

  A. A HIDDEN COLUMN. Seven columns beside a preview do not fit, and the one that
     falls off the right edge is R. Nothing errors: the user simply never sees a
     field, and a circle sized by an option they could not read is a circle they
     cannot check. Asserted as a property of the fitted widths -- every column's
     total fits the viewport at the size the dialog OPENS at, with no horizontal
     scroll bar -- so it holds for any font, any style and any row content rather
     than for the one screenshot someone looked at.
  B. A ROUNDED VALUE SAVED. Numbers are shown rounded for reading, which is only
     safe while the stored value survives. A circle whose Option is Intercept comes
     back from apply() with R = 41.23105625617661; if opening the editor and
     pressing OK wrote back the 41.2311 the cell displayed, every open would move
     the model a little, and nothing would report it. The check drives a real
     editor round trip on values that cannot survive six-digit rounding.
  C. A SILENT OVERWRITE. The generator replaces the table's contents. Rows already
     in it are the user's -- typed circles, a set from an earlier press -- and a
     button that discards them without asking is a button nobody can risk pressing.
     The check asserts the question is asked and that each of its three answers does
     what it says, including that Cancel leaves the table exactly as it was.

Also pins what the LEM-1 tutorial teaches word for word: the button's label, that
the generator fills an empty table without asking, and that its summary line -- the
audit the tutorial walks the reader through, dropped candidate included -- is
visible in the dialog afterwards.

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

#: The LEM-1 tutorial's model: one material, one profile line, a rigid base at the
#: toe. Its geometry is what the tutorial's Studio path generates a circle from.
LEM01 = os.path.join(_REPO, "docs/lem/files/xslope_simple_embankment.xlsx")

#: The button's label, quoted verbatim by the tutorial ("press **Generate starting
#: circles…**"). A rewording here is a tutorial that no longer matches the app.
GENERATE_LABEL = "Generate starting circles…"


def _fail(cond, message):
    return [] if cond else [message]


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _bare_model():
    """The tutorial's model as its geometry steps leave it: no failure surface."""
    d = _load(LEM01)
    d["circles"] = []
    d["circular"] = False
    return d


def _many_circles(n):
    """The same model carrying ``n`` circles — a table taller than any screen, where
    the row-number gutter is wider and the two panes cannot both have what they ask
    for. Both are ways a fit measured on a one-row table stops being a fit."""
    d = _load(LEM01)
    base = dict(d["circles"][0])
    d["circles"] = [dict(base, Xo=float(i), Yo=40.0 + i * 0.123456789)
                    for i in range(n)]
    return d


def _shown(dlg):
    """Show ``dlg`` offscreen and let the layout settle, so widths and heights are
    the ones a user would meet rather than pre-layout hints."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance()
    dlg.show()
    for _ in range(8):
        app.processEvents()
    return dlg


def _generate(dlg, answer=None):
    """Press Generate, answering any confirmation with ``answer`` (a QMessageBox
    standard button). Returns the list of confirmations that were raised."""
    from PySide6.QtWidgets import QMessageBox

    raised, real = [], QMessageBox.exec

    def fake(box):
        raised.append((box.text(), box.informativeText(),
                       [b.text() for b in box.buttons()]))
        return QMessageBox.Cancel if answer is None else answer

    QMessageBox.exec = fake
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            dlg._run_generate()
    finally:
        QMessageBox.exec = real
    return raised


# ---------------------------------------------------------------------------
# A. Every column is visible at the size the dialog opens at
# ---------------------------------------------------------------------------

def test_all_columns_visible():
    """No column falls off the right edge of the table the dialog opens with."""
    from studio.editors import CATEGORY_EDITORS

    editor = CATEGORY_EDITORS["circles"]
    out = []
    for label, model in (("with a circle", _load(LEM01)),
                         ("empty", _bare_model()),
                         ("forty circles", _many_circles(40))):
        dlg = _shown(editor.build(model, None))
        table = dlg._editable.table
        header = table.horizontalHeader()
        widths = [header.sectionSize(c) for c in range(table.columnCount())]
        out += _fail(table.columnCount() == len(editor.FIELDS),
                     f"{label}: table shows {table.columnCount()} of "
                     f"{len(editor.FIELDS)} columns")
        out += _fail(sum(widths) <= table.viewport().width(),
                     f"{label}: the columns total {sum(widths)} px in a "
                     f"{table.viewport().width()} px viewport -- "
                     f"{[f.header for f in editor.FIELDS][-1]} is off the edge")
        out += _fail(not table.horizontalScrollBar().isVisible(),
                     f"{label}: the table opens with a horizontal scroll bar")
        out += _fail(all(w > 0 for w in widths),
                     f"{label}: a column opened at zero width: {widths}")
        # Every header label must be legible, not elided to an ellipsis: a section
        # narrower than its own text is a column whose name the user cannot read.
        fm = header.fontMetrics()
        for c, f in enumerate(editor.FIELDS):
            out += _fail(widths[c] >= fm.horizontalAdvance(f.header),
                         f"{label}: column {f.header!r} is narrower than its header")
        dlg.deleteLater()
    return out


def test_preview_is_below_the_table():
    """The table is the full width of the dialog and the section sits under it."""
    from PySide6.QtCore import Qt
    from studio.editors import CATEGORY_EDITORS

    dlg = _shown(CATEGORY_EDITORS["circles"].build(_load(LEM01), None))
    out = _fail(dlg._split.orientation() == Qt.Vertical,
                "the circles preview is still beside the table, not below it")
    table_rect = dlg._editable.geometry()
    preview_rect = dlg._preview.geometry()
    out += _fail(preview_rect.top() >= table_rect.bottom(),
                 "the preview does not start below the table")
    out += _fail(abs(table_rect.width() - preview_rect.width()) <= 1,
                 f"the table ({table_rect.width()}) and the preview "
                 f"({preview_rect.width()}) are not the same width")
    # The caption under the preview must be fully inside the pane; a wrapped line
    # cut off along the bottom edge is text that is present and unreadable.
    caption = dlg._preview._caption
    out += _fail(caption is not None and
                 caption.geometry().bottom() <= dlg._preview.rect().bottom(),
                 "the preview caption is clipped by the bottom of its pane")
    floor = dlg._preview.minimumHeight()
    out += _fail(floor > 0 and dlg._preview.height() >= floor,
                 f"the preview opened {dlg._preview.height()} px tall against its own "
                 f"{floor} px minimum")
    dlg.deleteLater()

    # A table taller than the screen must not squeeze the preview out of the dialog:
    # the table is the pane with a scroll bar, so it is the pane that gives way.
    dlg = _shown(CATEGORY_EDITORS["circles"].build(_many_circles(40), None))
    # The splitter's allocation, not the widget's own height: a squeezed pane keeps
    # reporting the height it asked for while the splitter draws it at nothing.
    given = dlg._split.sizes()[1]
    out += _fail(given >= dlg._preview.minimumHeight(),
                 f"forty circles left the preview {given} px of the splitter, under "
                 f"its {dlg._preview.minimumHeight()} px minimum")
    out += _fail(dlg._editable.table.verticalScrollBar().isVisible(),
                 "forty circles did not put a scroll bar on the table, so the rows "
                 "beyond the pane cannot be reached")
    dlg.deleteLater()
    return out


# ---------------------------------------------------------------------------
# B. Display rounding never reaches the stored value
# ---------------------------------------------------------------------------

def test_display_rounds_but_stores_the_truth():
    """Cells read as rounded numbers; OK writes back the full-precision ones."""
    import math

    from studio.editors import CATEGORY_EDITORS

    editor = CATEGORY_EDITORS["circles"]
    cols = {f.key: j for j, f in enumerate(editor.FIELDS)}

    # (1) The display. An Intercept circle as CirclesEditor.apply leaves it: R from
    #     the intercept point, Depth from R. Neither survives six digits, and the
    #     raw repr of both is what the editor used to put on screen.
    r = math.hypot(10.0, 40.0)
    sd = _load(LEM01)
    sd["circles"] = [{"Xo": 10.0, "Yo": 40.0, "Option": "Intercept",
                      "Depth": 40.0 - r, "Xi": 0.0, "Yi": 0.0, "R": r}]
    dlg = editor.build(sd, None)
    shown_r = dlg._editable.table.item(0, cols["R"]).text()
    shown_d = dlg._editable.table.item(0, cols["Depth"]).text()
    out = _fail(shown_r == "41.2311",
                f"R is displayed as {shown_r!r}, not rounded for reading")
    out += _fail(shown_d == "-1.23106",
                 f"Depth is displayed as {shown_d!r}, not rounded for reading")
    out += _fail(len(shown_r) < len(repr(r)),
                 f"the displayed R {shown_r!r} is no shorter than repr({r})")
    dlg.deleteLater()

    # (2) The storage. Every value here is carried straight through apply() -- Xo and
    #     Yo are never recomputed, and with Option = Depth it is R that is derived
    #     from Depth rather than the other way round -- so what comes back out is
    #     what the editor decided to save, and every one of them is a number the
    #     display cannot hold.
    exact = {"Xo": 10.123456789012345, "Yo": 40.98765432109876,
             "Depth": 0.12345678901234567}
    sd = _load(LEM01)
    sd["circles"] = [dict(exact, Option="Depth", Xi=0.0, Yi=0.0,
                          R=exact["Yo"] - exact["Depth"])]
    dlg = editor.build(sd, None)
    for key, value in exact.items():
        shown = dlg._editable.table.item(0, cols[key]).text()
        out += _fail(shown != repr(value),
                     f"{key} is displayed as its full repr {shown!r}")
    editor.apply(sd, dlg)
    saved = sd["circles"][0]
    for key, value in exact.items():
        out += _fail(saved[key] == value,
                     f"OK stored {key} = {saved[key]!r}; the model held {value!r} "
                     f"and nothing in this dialog was edited")
    dlg.deleteLater()

    # (3) A real edit still wins: typing over a cell replaces the stored value.
    sd2 = _load(LEM01)
    sd2["circles"] = [{"Xo": 10.0, "Yo": 40.0, "Option": "Radius",
                       "Depth": 40.0 - r, "Xi": 0.0, "Yi": 0.0, "R": r}]
    dlg2 = editor.build(sd2, None)
    dlg2._editable.table.item(0, cols["R"]).setText("38.5")
    editor.apply(sd2, dlg2)
    out += _fail(sd2["circles"][0]["R"] == 38.5,
                 f"a typed R came back as {sd2['circles'][0]['R']!r}, not 38.5")
    dlg2.deleteLater()
    return out


# ---------------------------------------------------------------------------
# C. The Generate button
# ---------------------------------------------------------------------------

def test_generate_exists_and_is_labelled():
    """The button is there, is live on a model with geometry, and says what the
    tutorial says it says."""
    from studio.editors import CATEGORY_EDITORS

    dlg = CATEGORY_EDITORS["circles"].build(_bare_model(), None)
    button = getattr(dlg, "generate_button", None)
    if button is None:
        return ["the circles editor has no Generate button"]
    out = _fail(button.text() == GENERATE_LABEL,
                f"the button reads {button.text()!r}, not {GENERATE_LABEL!r} -- the "
                f"label the LEM-1 tutorial tells the reader to press")
    out += _fail(button.isEnabled(),
                 "the button is dimmed on a model whose geometry yields circles")
    out += _fail(bool(button.toolTip()),
                 "the live button carries no tooltip saying what it will build")
    dlg.deleteLater()
    return out


def test_generate_dims_without_geometry():
    """No geometry, no circles: the button dims and its tooltip says why."""
    from studio.editors import CATEGORY_EDITORS

    sd = _bare_model()
    sd["profile_lines"] = []
    sd["polygons"] = []
    sd["ground_surface"] = None
    dlg = CATEGORY_EDITORS["circles"].build(sd, None)
    out = _fail(not dlg.generate_button.isEnabled(),
                "the button is live on a model with no geometry at all")
    tip = dlg.generate_button.toolTip()
    out += _fail(len(tip) > 20,
                 f"the dimmed button's tooltip does not say why: {tip!r}")
    dlg.deleteLater()
    return out


def test_generate_fills_an_empty_table_without_asking():
    """The tutorial's own flow: an empty table, one press, rows and a summary.

    Nothing is at risk when the table is empty, so nothing is asked -- a confirmation
    here would be a dialog raised over a decision the user cannot get wrong. What
    the press must leave behind is the generator's summary, in the dialog: it is
    what the reader audits the row against.
    """
    from studio.editors import CATEGORY_EDITORS

    dlg = CATEGORY_EDITORS["circles"].build(_bare_model(), None)
    raised = _generate(dlg)
    out = _fail(not raised,
                f"pressing Generate on an EMPTY table asked {len(raised)} question(s)")
    rows = dlg.result_rows()
    out += _fail(bool(rows), "pressing Generate left the table empty")
    out += _fail(all(r.get("Option") == "Depth" for r in rows),
                 "a generated circle arrived without Option = Depth, so its stored "
                 "R and Depth disagree with the option that sizes it")
    out += _fail(all(isinstance(r.get("R"), float) and r["R"] > 0 for r in rows),
                 f"a generated circle has no positive radius: {rows}")

    summary = dlg.generate_summary
    out += _fail(summary.isVisibleTo(dlg), "the generator's summary is not shown")
    text = summary.text()
    out += _fail("circle" in text,
                 f"the summary does not say what was made: {text!r}")
    out += _fail("daylighting" in text,
                 f"the summary drops the candidate the generator rejected, which is "
                 f"the half of the audit that says what was NOT built: {text!r}")
    out += _fail(len(dlg._editable.result_rows()) == len(rows),
                 "the rows the summary describes are not the rows in the table")

    # What lands in the table is what the model gets on OK.
    sd = _bare_model()
    dlg2 = CATEGORY_EDITORS["circles"].build(sd, None)
    _generate(dlg2)
    proposed = dlg2.result_rows()
    CATEGORY_EDITORS["circles"].apply(sd, dlg2)
    out += _fail(len(sd["circles"]) == len(proposed) and sd["circular"],
                 "OK did not write the generated circles into the model")
    dlg.deleteLater()
    dlg2.deleteLater()
    return out


def test_generate_asks_before_touching_existing_rows():
    """With rows in the table the button asks, and each answer does what it says."""
    from PySide6.QtWidgets import QMessageBox
    from studio.editors import CATEGORY_EDITORS

    editor = CATEGORY_EDITORS["circles"]
    out = []

    # -- Cancel: the table is left exactly as it was. This is the mutation that
    #    matters -- an implementation that replaced the rows and asked nothing (or
    #    asked and ignored the answer) fails here and nowhere else.
    dlg = editor.build(_load(LEM01), None)
    before = dlg.result_rows()
    raised = _generate(dlg, QMessageBox.Cancel)
    out += _fail(len(raised) == 1,
                 f"pressing Generate over {len(before)} existing row(s) raised "
                 f"{len(raised)} confirmations; it must ask before touching them")
    if raised:
        text, informative, buttons = raised[0]
        out += _fail("Replace" in buttons,
                     f"the question offers no Replace: {buttons}")
        out += _fail("Append" in buttons,
                     f"the question offers no Append: {buttons}")
        out += _fail("Cancel" in buttons,
                     f"the question offers no way out: {buttons}")
        out += _fail(str(len(before)) in text,
                     f"the question does not say how much is at stake: {text!r}")
        out += _fail(bool(informative),
                     "the question does not say what the generator built")
    out += _fail(dlg.result_rows() == before,
                 "answering Cancel changed the table anyway")
    out += _fail(not dlg.generate_summary.isVisibleTo(dlg),
                 "a cancelled generate still reported that it did something")
    dlg.deleteLater()

    # -- Replace: the table holds the generated rows and nothing else.
    dlg = editor.build(_load(LEM01), None)
    _generate(dlg, QMessageBox.Apply)
    replaced = dlg.result_rows()
    out += _fail(len(replaced) == 1,
                 f"Replace left {len(replaced)} rows; the generator builds 1 on this "
                 f"model")
    dlg.deleteLater()

    # -- Append: the existing rows are still there, with the generated ones after.
    dlg = editor.build(_load(LEM01), None)
    before = dlg.result_rows()
    _generate(dlg, QMessageBox.Yes)
    after = dlg.result_rows()
    out += _fail(len(after) == len(before) + len(replaced),
                 f"Append left {len(after)} rows, not {len(before)} + "
                 f"{len(replaced)}")
    out += _fail(after[:len(before)] == before,
                 "Append did not keep the rows that were already in the table")
    out += _fail(dlg.generate_summary.isVisibleTo(dlg)
                 and "Added" in dlg.generate_summary.text(),
                 f"the summary does not say rows were added: "
                 f"{dlg.generate_summary.text()!r}")
    dlg.deleteLater()
    return out


def test_generate_is_the_remedys_generator():
    """The button and the Run LEM remedy call one function, so they cannot disagree.

    Not an import-path assertion: the circles the button puts in the table are
    compared against the circles the preflight remedy writes into the model, on the
    same model, field for field.
    """
    from xslope.preflight import preflight
    from studio.editors import CATEGORY_EDITORS

    sd = _bare_model()
    dlg = CATEGORY_EDITORS["circles"].build(sd, None)
    _generate(dlg)
    from_button = dlg.result_rows()

    with contextlib.redirect_stdout(io.StringIO()):
        report = preflight(_bare_model(), "lem",
                           remedies=["generate_starting_circles"])
    from_remedy = (report.model or {}).get("circles") or []
    out = _fail(bool(report.applied),
                f"the Run LEM preflight applied no circle remedy: {report.applied}")
    out += _fail(len(from_button) == len(from_remedy),
                f"the button built {len(from_button)} circle(s), the Run LEM remedy "
                f"{len(from_remedy)}")
    for i, (a, b) in enumerate(zip(from_button, from_remedy)):
        for key in ("Xo", "Yo", "Depth", "R"):
            out += _fail(abs(float(a[key]) - float(b[key])) < 1e-12,
                         f"circle {i} {key}: button {a[key]!r} vs remedy {b[key]!r}")
    dlg.deleteLater()
    return out


CHECKS = [
    ("every column is visible at the opening size", test_all_columns_visible),
    ("the preview sits below the table", test_preview_is_below_the_table),
    ("display rounds, storage does not", test_display_rounds_but_stores_the_truth),
    ("generate button exists and is labelled", test_generate_exists_and_is_labelled),
    ("generate dims with no geometry", test_generate_dims_without_geometry),
    ("generate fills an empty table", test_generate_fills_an_empty_table_without_asking),
    ("generate asks before overwriting", test_generate_asks_before_touching_existing_rows),
    ("button and remedy share a generator", test_generate_is_the_remedys_generator),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("circles editor: PySide6 not installed — skipped.")
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
    print("circles editor:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll circles editor checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
