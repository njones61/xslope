"""W-1, The AI Assistant: the dock, the settings dialog, and every "Check its
work" item a reader can do without the assistant.

Not scriptable here: the conversations themselves. Each is a call to a paid
model provider whose replies differ from run to run, and the page says so — it
shows one recorded session per use case and asks the reader to check the result,
not to match the wording. So the walk does what each "Check its work" list asks
of the saved workbooks the page links: reload, search, read the editors, run the
strength reduction, open the report. The sweep and ad hoc sections' re-runs, the
cost table and the "Other models" snapshots are not walked (the first two repeat
LEM-3's numbers, which LEM-3's walk holds; the last two are records of past
sessions).
"""
import hashlib
import re

from tutorial_walks.harness import StepFailed, field, repo_path, set_field

PAGE = "w01_ai_assistant.md"
TIMEOUT = 3600


def _search(st, path, w):
    st.open(w.copy_model(w.page_link(path), path[:-5]))
    st.set_mode("LEM")
    st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Auto search")))
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("%s: the search left no solution" % path)
    return b


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Setting it up: the dock"):
        from PySide6.QtWidgets import QDockWidget, QPushButton
        docks = {d.windowTitle(): d for d in st.mw.findChildren(QDockWidget)}
        w.require("Assistant" in docks, "there is no Assistant dock (docks: %s)" % sorted(docks))
        texts = {b.text() for b in docks["Assistant"].findChildren(QPushButton)}
        for t in ("New chat", "Settings…", "Send"):
            w.check(t in texts, "the Assistant dock has no %r button (%s)" % (t, sorted(texts)))

    with w.step("Building a model from a drawing: check its work"):
        b = _search(st, "w1_build_from_image_after.xlsx", w)
        w.expect("search", b["results"]["FS"], r"\*\*Reopen the workbook and search\.\*\* (1\.\d+), on a circle")
        w.check(abs(b["search"]["fs_cache"][0]["Depth"]) < 1e-6, "the circle is not tangent to the contact")
        mats = [(m["name"], float(m["gamma"]), float(m["c"]), float(m["phi"]), m["u"]) for m in st.data["materials"]]
        w.check([x[1:] for x in mats] == [(130.0, 400.0, 0.0, "none"), (135.0, 800.0, 0.0, "none")],
                "the materials read %s" % mats)
        lines = st.data["profile_lines"]
        xs = [float(x) for x, _y in lines[-1]["coords"]]
        w.expect("foundation line left end", min(xs), r"The foundation line runs from x = (−\d+) to")
        w.expect("foundation line right end", max(xs), r"runs from x = −\d+ to (\d+) and the")

    with w.step("Modifying the model: check its work"):
        fs = p.find(r"\*\*Reload each workbook and search\.\*\* (\d\.\d+), (\d\.\d+) and (\d\.\d+)\.").groups()
        lows = p.find(r"Elevation (\d+), then (\d+), then (−\d+)\.").groups()
        for i, name in enumerate(("w1_modify_after_1.xlsx", "w1_modify_after_2.xlsx", "w1_modify_after_3.xlsx")):
            b = _search(st, name, w)
            w.expect_written("%s search" % name, b["results"]["FS"], fs[i])
            w.expect_written("%s lowest point" % name, b["search"]["fs_cache"][0]["Depth"], lows[i])
            if i == 0:
                c = st.data["circles"]
                w.check(len(c) == 3 and sorted(float(x["Depth"]) for x in c) == [-10.0, -10.0, 0.0],
                        "the first workbook's circles are %s" % [x["Depth"] for x in c])
            if i == 1:
                w.check(all(m["u"] == "none" for m in st.data["materials"]), "u is not none on both soils")

    with w.step("Stiffnesses and a strength reduction run: check its work"):
        st.open(w.copy_model(w.page_link("w1_elastic_fem_after.xlsx"), "fem"))
        st.set_mode("FEM")
        st.build_mesh(lambda dlg: None)
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"tri6 mesh at target size 1\.2 — ([\d,]+) nodes")
        st.run_fem(lambda dlg: None)
        r = (st.results.get("fem_solution") or {}).get("ssrm_result")
        w.require(r is not None, "the run left no result")
        w.expect("FS", r["FS"], r"\*\*Read (1\.\d+) as a bracket\.\*\*")
        lo, hi = r["final_interval"]
        w.expect("bracket low", lo, r"\[(1\.\d+), 1\.\d+\] is")
        w.expect("bracket high", hi, r"\[1\.\d+, (1\.\d+)\] is")

    with w.step("A broken file: check its work"):
        st.open(w.copy_model(w.page_link("w1_diagnose_start.xlsx"), "diag0"))
        w.check(float(st.data.get("max_depth")) == float(p.number(r"the file that came in read (−\d+)")[0]),
                "the start file's maximum depth is %s" % st.data.get("max_depth"))
        b = _search(st, "w1_diagnose_after.xlsx", w)
        w.expect("repaired search", b["results"]["FS"], r"\*\*Reload the repaired workbook and search\.\*\* (1\.\d+) on")
        w.expect("repaired maximum depth", float(st.data.get("max_depth")), r"It reads (−\d+) in the repaired file")
        w.check(st.data["materials"][0]["name"].lower().startswith("foundation"),
                "the foundation is not listed first")

    with w.step("Generating a report: check its work"):
        import zipfile
        path = repo_path(w.page_link("w1_report_after.docx"))
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8", "replace")
        text = " ".join(re.sub(r"<[^>]+>", " ", xml).split())
        quoted = p.find(r"\*\*Open the report at §3\.4\.2\.\*\* \"(.+?)\" —").group(1)
        w.check(" ".join(quoted.split()) in text.replace(" ,", ","), "the report does not say %r" % quoted)
        wb = re.search(r"\[[^\]]*\]\(((?:\.\./[a-z]+/)?files/xslope_simple_mult_layers\.xlsx)\)", p.text)
        if wb:
            import os
            from tutorial_walks.harness import TUTORIALS
            digest = hashlib.sha256(open(os.path.normpath(os.path.join(TUTORIALS, wb.group(1))), "rb").read()).hexdigest()
            w.check(digest in text, "the report's traceability digest is not the linked workbook's")
        else:
            w.fail("the page does not link the workbook the report was written from")
