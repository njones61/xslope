"""W-3, Generating the Analysis Report: both reports built in Studio.

Walked: the solved workbook opened with its companions (the solutions Studio
restores), the quick Spencer run, File → Generate Report… with the page's
choices, the generated document read against the page (the factors of safety
table, the attached seepage and strength reduction numbers, Spencer's closing
F and θ, Bishop's two sums, the traceability digest), and the second report on
the example template (its letterhead in the headers and footers, the title page
from the Organization field). Not walked: the page counts and contents-page
numbers (they come from the Word or LibreOffice finishing step, which this walk
does not depend on), and the .xslz download (the package is built with the
docs, not shipped in the repository).
"""
import hashlib
import os
import re
import zipfile

from tutorial_walks.harness import StepFailed, field, repo_path, set_field, value

PAGE = "w03_report.md"
TIMEOUT = 1800


def _docx_text(path, part="word/document.xml"):
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if re.fullmatch(part, n)]
        xml = " ".join(z.read(n).decode("utf-8", "replace") for n in names)
    return " ".join(re.sub(r"<[^>]+>", " ", xml).split())


def _tick_method(dlg, prefix, on):
    from PySide6.QtCore import Qt
    lw = dlg.methods
    hits = [i for i in range(lw.count()) if lw.item(i).text().startswith(prefix)]
    if len(hits) != 1:
        raise StepFailed("Methods lists no single %r (%s)" % (prefix, [lw.item(i).text() for i in range(lw.count())]))
    lw.item(hits[0]).setCheckState(Qt.Checked if on else Qt.Unchecked)


def _generate(st, action, timeout=900):
    st.trigger("File", "Generate Report", action)
    st.wait(lambda: st.mw._report_runner is None, timeout, "the report")


def walk(w):
    p = w.page
    st = w.studio()
    title = p.find(r"Fill in \*\*Project title\*\* `([^`]+)`").group(1)
    number = p.find(r"\*\*Project number\*\* `([^`]+)`").group(1)

    with w.step("The model"):
        path = w.copy_model(w.page_link("xslope_johnson_res_solved.xlsx"))
        mark = st.mark()
        st.open(path)
        log = st.log_since(mark)
        w.check("Restored saved seepage solution" in log, "the seepage solution is not restored on open")
        m = re.search(r"Restored saved FEM solution .*?\(SSRM FS = (\d\.\d+)\)", log)
        w.check(m is not None, "the FEM solution is not restored on open")
        w.check(st.mode == "fem", "the file opens in %s mode" % st.mode)
        st.set_mode("LEM")
        st.run_lem(lambda dlg: None)
        w.expect("quick Spencer run", st.results["lem_solution"]["results"]["FS"],
                 r"The search returns \*\*FS = (1\.\d+)\*\*, COMBO-1's own answer")

    out = os.path.join(os.path.dirname(path), "xslope_johnson_res_solved_report.docx")

    with w.step("Generating the report"):
        def act(dlg):
            w.check(value(field(dlg, "Format")) == "Word document (.docx)", "Format opens on %s" % value(field(dlg, "Format")))
            w.check(dlg.path.text().endswith("_report.docx") and
                    os.path.dirname(dlg.path.text()) == os.path.dirname(path),
                    "Save as does not default to <model>_report.docx beside the model: %s" % dlg.path.text())
            w.check(dlg.template.text() in ("", "Shipped template") or "Shipped" in dlg.template.placeholderText(),
                    "Template does not read Shipped template: %r" % dlg.template.text())
            _tick_method(dlg, "Bishop", True)
            dlg.title.setText(title)
            dlg.project_number.setText(number)
            dlg.organization.setText(p.find(r"\*\*Organization\*\* `([^`]+)`").group(1))
            dlg.author.setText("A. Reader")
        _generate(st, act)
        w.require(os.path.exists(out), "no report was written at %s" % out)

    with w.step("What the document contains"):
        text = _docx_text(out)
        tbl = p.table(r"records that each searched for a surface of its own\.")
        for r in tbl[1:]:
            w.check(r[2].split()[0] in text, "the report does not state %s (%s)" % (r[2], r[0]))
        for label, pat in (("Spencer's closing F", r"they close at F = (\d\.\d+) and θ"),
                           ("Spencer's closing θ", r"and θ = (−\d+\.\d+)°"),
                           ("Bishop's first sum", r"Bishop divides its own two sums, (\d+) and"),
                           ("Bishop's second sum", r"two sums, \d+ and (\d+),")):
            want = p.find(pat).group(1).replace("−", "-")
            w.check(want in text.replace("−", "-"), "the report does not print %s %s" % (label, want))
        digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
        w.check(digest in text, "the traceability stamp is not the opened workbook's SHA-256")

    with w.step("A company template"):
        tpl = os.path.join(w.workdir, "report_template_example.docx")
        import shutil
        shutil.copy2(repo_path(w.page_link("report_template_example.docx")), tpl)
        out2 = os.path.join(os.path.dirname(path), p.find(r"Change \*\*Save as\*\* to `([^`]+)`").group(1))
        org = p.find(r"set \*\*Organization\*\* to `([^`]+)`").group(1)

        def act(dlg):
            from PySide6.QtWidgets import QFileDialog, QPushButton
            browse = [b for b in dlg.findChildren(QPushButton) if b.text() == "Browse…"]
            w.require(len(browse) == 2, "the dialog has %d Browse… buttons" % len(browse))
            real = QFileDialog.getOpenFileName
            QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (tpl, ""))
            try:
                browse[1].click()
            finally:
                QFileDialog.getOpenFileName = real
            w.check(dlg.template.text() == tpl, "the Template field reads %r after Browse…" % dlg.template.text())
            dlg.path.setText(out2)
            w.check(dlg.project_number.text() == "", "Project number does not reopen empty")
            _tick_method(dlg, "Bishop", False)
            dlg.title.setText(title)
            dlg.project_number.setText(number)
            dlg.organization.setText(org)
        _generate(st, act)
        w.require(os.path.exists(out2), "no report was written at %s" % out2)
        heads = _docx_text(out2, r"word/header\d+\.xml")
        foots = _docx_text(out2, r"word/footer\d+\.xml")
        w.check("ACME Geotechnical" in heads, "the firm's name is not in the report's headers")
        w.check("ACME Geotechnical" in foots, "the firm's line is not in the report's footers")
        w.check(org.upper() in _docx_text(out2).upper(), "the title page does not read %s" % org)
