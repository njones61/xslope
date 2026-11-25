"""LEM-15, Yield Coefficient and Newmark Displacement, in Studio.

Walked, on copies of the linked model and record: the material read against the
page's row; Run LEM with Find the yield coefficient k_y and the two Log lines the
page quotes; Run → Newmark… with the record chosen through Browse… (the file
dialog answered with the copied record), its header details, the k_y the dialog
carries over, the as-recorded run and its Log line; the record scaled by 2 and
by 3 and the scale table (peak, displacement both ways, slip episodes). Not
walked: the third candidate surface (the page says xslope does not compute it).
"""
import contextlib
import re

from tutorial_walks.harness import (StepFailed, editable, field, headers, page_rows, press,
                                    rows, set_field, set_value, show_parameters_for,
                                    tick, value)

PAGE = "lem15_yield_newmark.md"
TIMEOUT = 1800


@contextlib.contextmanager
def _file_dialog_answers(path):
    """QFileDialog.getOpenFileName answered with ``path`` (the reader's pick)."""
    from PySide6.QtWidgets import QFileDialog
    real = QFileDialog.getOpenFileName
    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (path, ""))
    try:
        yield
    finally:
        QFileDialog.getOpenFileName = real


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_pseudostatic_start.xlsx")))
        want = page_rows(p.table(r"set the \*\*Show parameters for:\*\* toggles to \*\*LEM\*\* alone:"))[0]
        seen = {}

        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            t = editable(dlg)
            i = headers(t)[0].index("name")
            seen["row"] = rows(t)[0][i:i + len(want)]
            return "reject"
        st.edit("Materials", act)
        for a, b in zip(want, seen["row"]):
            w.check(a == b or (a == "" and b in ("", "0")) or (_f(a) is not None and _f(a) == _f(b)),
                    "the material reads %s, not the page's %s" % (seen["row"], want))

    with w.step("Finding the yield coefficient"):
        def act(dlg):
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
            w.check(int(value(field(dlg, "Number of slices"))) == 40, "slices are not 40")
            tick(dlg, "Find the yield coefficient k_y", True)
        _s, log = st.run_lem(act)
        w.check_log(p.code_block(r"It finds k<sub>y</sub> on both surfaces\. Click \*\*Run\*\*\. The Log reports:"),
                    log, "yield coefficient")
        b = st.results.get("lem_solution") or {}
        w.check(abs(float((b.get("results") or {}).get("FS") or 0) - 1.0) < 0.001,
                "the LEM Solution view's circle is not solved to FS = 1 (%s)" % (b.get("results") or {}).get("FS"))
        ky = p.number(r"\*\*Yield coefficient k_y\*\* is filled from the run in Part 1, (\d\.\d+)")[0]

    record = w.copy_model(w.page_link("anderson_dam_1989_243.txt"), "record")

    def newmark(scaling, val=None, browse=False, seen=None):
        def act(dlg):
            if browse:
                with _file_dialog_answers(record):
                    press(dlg, "Browse")
            w.require(dlg.path.text() == record, "the record field reads %r" % dlg.path.text())
            if seen is not None:
                seen["info"] = dlg.info.text()
                seen["ky"] = value(field(dlg, "Yield coefficient k_y"))
            tick(dlg, scaling, True)
            if scaling == "Scale factor":
                set_value(dlg.factor, val)
        mark = st.mark()
        st.trigger("Run", "Newmark", act)
        r = st.results.get("newmark")
        if r is None:
            raise StepFailed("the Newmark run left no result")
        return r, st.log_since(mark)

    with w.step("The record"):
        head = p.find(r"horizontal component 243°: ([\d,]+) points at (\d\.\d+) s, with a peak of (\d\.\d+) g").groups()
        txt = open(record, encoding="utf-8", errors="replace").read().splitlines()
        data = [l for l in txt if re.match(r"^\s*-?\d", l)]
        w.check(len(data) == int(head[0].replace(",", "")), "the record has %d points, not %s" % (len(data), head[0]))

    with w.step("Running the record as recorded"):
        seen = {}
        r, log = newmark("As recorded", browse=True, seen=seen)
        for part in (head[0], head[1], head[2]):
            w.check(part in seen["info"], "the dialog's record line does not show %s: %r" % (part, seen["info"]))
        w.check(abs(float(seen["ky"]) - float(ky)) < 5e-5, "k_y is filled with %s, not %s" % (seen["ky"], ky))
        w.check_log(p.code_block(r"each one adds a step to the displacement\. The Log reports:"), log, "as recorded")

    with w.step("Scaling the record"):
        tbl = p.table(r"It scales the record so that its peak equals the value entered\.")
        for row in tbl[1:]:
            if row[0] == "1":
                r1 = r
            else:
                r1, _log = newmark("Scale factor", float(row[0]))
            dirs = r1["directions"]
            gov, other = dirs[r1["governing"]], dirs[-r1["governing"]]
            at = r"\| %s \| " % row[0]
            w.expect("peak at %s×" % row[0], r1["pga"], at + r"(\d\.\d+) \|")
            w.expect("D at %s×" % row[0], r1["D"], at + r"\d\.\d+ \| (\d+\.\d+) \|")
            w.expect("other direction at %s×" % row[0], other["D"],
                     at + r"\d\.\d+ \| \d+\.\d+ \| (\d+\.\d+) \|")
            w.expect("slip episodes at %s×" % row[0], gov["episodes"],
                     at + r"\d\.\d+ \| \d+\.\d+ \| \d+\.\d+ \| (\d+) \|")

def _f(v):
    try:
        return float(str(v).replace("−", "-"))
    except (TypeError, ValueError):
        return None
