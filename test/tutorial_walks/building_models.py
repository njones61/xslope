"""Tutorial 0, Building Models Three Ways: the machinery the page describes.

The page is a reference page with no model of its own to build and no number to
compute; it tells the reader what the template holds and what Studio does on an
edit, an undo, a save and a package. Walked: the template's sheets and version;
the Inputs tree's categories on the completed LEM-1 model (a copy); an edit
marked unsaved, labeled on the undo history, and undone back to clean; File →
New then File → Save opening the Save As dialog; Export Project Package… (its
dialog, and the offer to save first after an edit); opening the package (the
Open project package dialog, then Open Existing / Extract Fresh on the second
open); and the Assistant dock's Settings… dialog offering Ollama. Not walked:
the assistant's conversation (a provider call; W-1 records real sessions).
"""
import contextlib
import os

from tutorial_walks.harness import StepFailed, button, repo_path

PAGE = "building_models.md"
TIMEOUT = 900


@contextlib.contextmanager
def _file_dialogs(calls, save_path=""):
    """The static file pickers answered (open: nothing; save: ``save_path``),
    each call recorded as (kind, title, filter)."""
    from PySide6.QtWidgets import QFileDialog
    real = (QFileDialog.getOpenFileName, QFileDialog.getSaveFileName)
    QFileDialog.getOpenFileName = staticmethod(
        lambda *a, **k: (calls.append(("open",) + tuple(a[1:4])), ("", ""))[1])
    QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (calls.append(("save",) + tuple(a[1:4])), (save_path, ""))[1])
    try:
        yield
    finally:
        QFileDialog.getOpenFileName, QFileDialog.getSaveFileName = real


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("The Excel template"):
        import openpyxl
        wb = openpyxl.load_workbook(repo_path("docs/inputs/input_template.xlsx"), read_only=True)
        sheets = set(wb.sheetnames)
        for name in ("main", "mat", "profile", "polygon", "piezo", "dloads", "reinforce", "joints",
                     "piles", "circles", "non-circ", "seep bc", "tseep"):
            w.check(name in sheets, "the template has no %r sheet" % name)
        w.check(w.page_link("input_template.xlsx") == "docs/inputs/input_template.xlsx",
                "the page does not link the master template")

    model = w.copy_model(w.page_link("xslope_simple_embankment.xlsx"))

    with w.step("Studio's editors: the Inputs tree"):
        st.open(model)
        items = st.tree_items()
        for cat in ("Global parameters", "Materials", "Profile lines", "Circles"):
            w.check(cat in items, "the Inputs tree has no %r" % cat)

    with w.step("Every edit is undoable; unsaved changes are marked"):
        def act(dlg):
            from tutorial_walks.harness import column, editable, press
            press(dlg, "Table view")
            t = editable(dlg)
            it = t.table.item(0, column(t, "γ"))
            if it is None:
                raise StepFailed("the γ cell cannot be typed into")
            it.setText(str(float(it.text() or 0) + 1.0))
        st.edit("Materials", act)
        w.check("*" in st.mw.windowTitle(), "an edited project shows no asterisk: %r" % st.mw.windowTitle())
        label = p.find(r"labeled steps — \*\"(Edit Materials)\"\*").group(1)
        history = _undo_history(st)
        w.check(any(label in h for h in history),
                "the undo history does not list %r: %s" % (label, history))
        st.trigger("Edit", "Undo")
        w.check("*" not in st.mw.windowTitle(), "undoing to the saved state does not clear the asterisk")

    with w.step("Saving a project: Save on a new project"):
        calls = []
        with _file_dialogs(calls):
            st.trigger("File", "New")
            st.trigger("File", "Save")
        w.check(calls and calls[-1][0] == "save",
                "File → Save on a never-saved project does not open the Save As dialog (%s)" % calls)

    with w.step("Sending a project: Export Project Package"):
        st.open(model)
        out = os.path.splitext(model)[0] + ".xslz"
        seen = {}

        def pkg(dlg):
            seen["title"] = dlg.windowTitle()
            seen["suffix"] = dlg.defaultSuffix()
            seen["files"] = list(dlg.selectedFiles())
            dlg.selectFile(out)
        st.trigger("File", "Export Project Package", pkg)
        w.check(seen.get("suffix") == "xslz", "the save dialog's suffix is %r" % seen.get("suffix"))
        w.check(any(f.endswith(os.path.basename(out)) for f in seen.get("files") or []),
                "the save dialog does not fill in %s" % os.path.basename(out))
        w.require(os.path.exists(out), "no package was written")

        def edit(dlg):
            from tutorial_walks.harness import column, editable, press
            press(dlg, "Table view")
            t = editable(dlg)
            it = t.table.item(0, column(t, "γ"))
            it.setText(str(float(it.text() or 0) + 1.0))
        st.edit("Materials", edit)
        n = len(st.dialogs.boxes)
        st.trigger("File", "Export Project Package")
        asked = [b for b in st.dialogs.boxes[n:] if b[0] == "question"]
        w.check(asked and "Save" in asked[-1][1], "after an edit Export does not offer to save first (%s)" % asked)
        st.trigger("Edit", "Undo")

    with w.step("Opening a package"):
        seen = {}

        def first(dlg):
            seen["title"] = dlg.windowTitle()
            seen["dest"] = dlg.dest.text()
            button(dlg, "Change")
            dlg.btn_unpack.click()
            return "keep"
        st.run_queued([first], lambda: st.mw.open_path(out), "open the package")
        st.settle(10)
        w.check(seen.get("dest") == os.path.splitext(out)[0],
                "the dialog does not offer the folder named for the package: %r" % seen.get("dest"))
        w.check(os.path.dirname(st.doc.path) == os.path.splitext(out)[0],
                "the opened workbook is not the extracted one: %s" % st.doc.path)

        def again(dlg):
            seen["buttons"] = [b.text() for b in (dlg.btn_existing, dlg.btn_fresh, dlg.btn_unpack)
                               if b.isVisibleTo(dlg)]
            dlg.btn_fresh.click()
            return "keep"
        st.run_queued([again], lambda: st.mw.open_path(out), "open the package again")
        st.settle(10)
        w.check(seen.get("buttons") == ["Open Existing", "Extract Fresh"],
                "the second open offers %s" % seen.get("buttons"))
        w.check(os.path.dirname(st.doc.path) != os.path.splitext(out)[0],
                "Extract Fresh did not unpack into a numbered folder beside it")

    with w.step("The AI assistant: Settings…"):
        from PySide6.QtWidgets import QComboBox, QPushButton
        btns = [b for b in st.mw.findChildren(QPushButton) if b.text() == "Settings…"]
        w.require(btns, "the Assistant dock has no Settings… button")
        seen = {}

        def settings(dlg):
            seen["providers"] = [c.itemText(i) for c in dlg.findChildren(QComboBox) for i in range(c.count())]
            return "reject"
        st.run_queued([settings], btns[0].click, "Assistant Settings…")
        w.check(any(t.startswith("Ollama") for t in seen.get("providers") or []),
                "the provider list offers no Ollama: %s" % seen.get("providers"))


def _undo_history(st):
    """The labels the toolbar Undo button's dropdown lists."""
    from PySide6.QtWidgets import QToolButton
    out = []
    for b in st.mw.findChildren(QToolButton):
        if b.text().replace("&", "") == "Undo" and b.menu() is not None:
            b.menu().aboutToShow.emit()
            st.settle(2)
            out += [a.text() for a in b.menu().actions()]
    return out
