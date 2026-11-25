"""LEM-2, Loads on the Crest: the surcharge added in Studio to a copy of LEM-1's
completed model, then the page's runs.

Walked: path C (Studio) — the distributed load entered with the page's Add row
count and Direction left at Normal — the Spencer search and the loaded-against-
unloaded table, the other three moment methods on the same circle, the
interslice forces the warning strip quotes, the load replaced by a line load
(Remove, then Line loads) and its run, the two face loads of "Direction of a
distributed load" (through the engine: the page gives their result, not the
clicks that enter them), Edit → Undo twice back to the distributed load, and the
Parametric design sweep on c. Not walked: path A (the assistant), path B (the
dloads sheet cell by cell; its workbook is ../lem/files/xslope_crest_surcharge.xlsx,
which the page's tags run), and the factor-of-safety-by-method table (tags).
"""
import re

from tutorial_walks.harness import (StepFailed, button, editable, field, page_rows,
                                    paste, press, set_field, type_cells,
                                    value)

PAGE = "lem02_loads_on_the_crest.md"
TIMEOUT = 3600


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _crit(b):
    return b["search"]["fs_cache"][0]


def _radius(c):
    return c["R"] if c.get("R") is not None else c["Yo"] - c["Depth"]


def _lem(method, analysis="Auto search"):
    def act(dlg):
        set_field(dlg, "Method", method)
        set_field(dlg, "Analysis", analysis)
    return act


def _interslice(b):
    """(most tensile, largest compression) the warning strip prints."""
    text = " ".join(str(x) for x in (b["results"].get("warnings") or [])).replace("−", "-")
    m = re.search(r"min Z = (-?[\d.]+) vs max compression ([\d.]+)", text)
    if not m:
        raise StepFailed("the warnings give no interslice forces: %r" % text)
    return float(m.group(1)), float(m.group(2))


def walk(w):
    p = w.page
    st = w.studio()
    load = page_rows(p.table(r"Two points define the load here"))

    with w.step("Open a copy of LEM-1's completed model"):
        st.open(w.copy_model(w.page_link("xslope_simple_embankment.xlsx")))
        unloaded = w.lem(dict(st.data), "spencer")

    with w.step("C: the distributed load"):
        n = p.count(r"1\. Press \*\*Add row\*\* (\w+) and enter")
        seen = {}

        def act(dlg):
            from PySide6.QtWidgets import QTabWidget
            tabs = dlg.findChildren(QTabWidget)[0]
            seen["tabs"] = [tabs.tabText(i) for i in range(tabs.count())]
            press(dlg, "Add load")
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(load), "Add row %d times gives %d rows for %d "
                      "points" % (n, t.table.rowCount(), len(load)))
            type_cells(t, load)
            seen["dir"] = value(field(dlg, "Direction"))
        st.edit("Distributed loads", act)
        w.check(seen["tabs"] == ["Set 1", "Set 2 (rapid drawdown)"],
                "the editor's tabs read %s" % seen["tabs"])
        w.check(seen["dir"] == "Normal (perpendicular to the line)",
                "Direction opens on %r" % seen["dir"])
        dl = st.data.get("dloads") or []
        w.require(len(dl) == 1 and [(pt["X"], pt["Y"], pt["Normal"]) for pt in dl[0]] ==
                  [tuple(float(c) for c in r) for r in load], "the load saved as %s" % dl)

    with w.step("Running the analysis"):
        st.run_lem(_lem("Spencer"))
        b = _bundle(st)
        w.expect("loaded Spencer FS", b["results"]["FS"], r"\*\*FS = (0\.\d+)\.\*\* The unloaded slope")
        drop = 100.0 * (1 - b["results"]["FS"] / unloaded["results"]["FS"])
        w.expect("drop from the unloaded slope (%)", drop, r"takes (\d+)% off it")
        rows_ = [("Factor of safety (Spencer)", lambda x: x["results"]["FS"]),
                 ("Circle radius (ft)", lambda x: _radius(_crit(x))),
                 ("Surface exits the crest at x =", lambda x: max(c[0] for c in x["failure_surface"].coords)),
                 ("Length of the surface (ft)", lambda x: float(x["slice_df"]["dl"].sum())),
                 ("Weight of the sliding mass (lb/ft)", lambda x: float(x["slice_df"]["w"].sum()))]
        for label, get in rows_:
            at = r"\| %s \| " % re.escape(label)
            w.expect("unloaded: %s" % label, get(unloaded), at + r"([\d,.]+) \|")
            w.expect("loaded: %s" % label, get(b), at + r"[\d,.]+ \| ([\d,.]+) \|")
        loaded = b

    with w.step("The warnings from Tutorial LEM-1"):
        xo = _crit(loaded)
        for meth in ("Bishop's Simplified", "Ordinary Method of Slices", "Morgenstern-Price"):
            st.run_lem(_lem(meth))
            bb = _bundle(st)
            w.expect("%s FS" % meth, bb["results"]["FS"], r"all three return \*\*(0\.\d+) on the same")
            c = _crit(bb)
            w.check(abs(c["Xo"] - xo["Xo"]) < 1e-6 and abs(c["Yo"] - xo["Yo"]) < 1e-6,
                    "%s is not on Spencer's circle" % meth)
        t1, c1 = _interslice(unloaded)
        w.expect("LEM-1 most tensile interslice force", t1, r"most tensile interslice force was (−\d+) lb/ft")
        w.expect("LEM-1 largest compression", c1, r"against a largest compression of (\d+) — ")
        w.expect("LEM-1 tension share (%)", 100 * abs(t1) / c1, r"of (?:\d+) — (\d+)% of it")
        t2, c2 = _interslice(loaded)
        w.expect("loaded most tensile interslice force", t2, r"Under the surcharge it is (−\d+) lb/ft")
        w.expect("loaded largest compression", c2, r"lb/ft against (\d+), or")
        w.expect("loaded tension share (%)", 100 * abs(t2) / c2, r"against \d+, or (\d+)%")

    with w.step("The same force as a line load"):
        def rm(dlg):
            from PySide6.QtWidgets import QTabWidget
            s1 = dlg.findChildren(QTabWidget)[0].widget(0)
            s1.list.setCurrentRow(0)
            press(s1, "Remove")
        st.edit("Distributed loads", rm)
        w.require(not st.data.get("dloads"), "the distributed load is still there")
        row = page_rows(p.table(r"one row, in the columns the `lloads` worksheet"))

        def ll(dlg):
            press(dlg, "Add row")
            paste(editable(dlg), row)
        st.edit("Line loads", ll)
        lls = st.data.get("line_loads") or []
        w.require(len(lls) == 1, "%d line loads saved" % len(lls))
        st.run_lem(_lem("Spencer"))
        b = _bundle(st)
        w.expect("line-load FS", b["results"]["FS"], r"\*\*FS = (0\.\d+)\*\*, against 0\.919")
        xs = [c[0] for c in b["failure_surface"].coords]
        w.expect("line-load crest exit", max(xs), r"It exits the crest at x = (\d+\.\d) — just past")
        arc = float(b["slice_df"]["dl"].sum())
        arc0 = float(loaded["slice_df"]["dl"].sum())
        w.expect("shorter arc", arc0 - arc, r"That is (\d+\.\d) ft less arc")
        w.expect("c × arc, line load", 500.0 * arc, r"c times the arc length: ([\d,]+) lb/ft against")
        w.expect("c × arc, strip load", 500.0 * arc0, r"lb/ft against ([\d,]+)\. The two loads")

    with w.step("Direction of a distributed load (engine)"):
        base = dict(st.data, line_loads=[])
        face = [{"X": 5.0, "Y": 5.0, "Normal": 750.0}, {"X": 15.0, "Y": 15.0, "Normal": 750.0}]
        bn = w.lem(dict(base, dloads=[face], dload_dirs=["normal"]), "spencer")
        w.expect("face load, normal", bn["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* with \*\*Direction\*\* = `normal`")
        bv = w.lem(dict(base, dloads=[face], dload_dirs=["vertical"]), "spencer")
        w.expect("face load, vertical", bv["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* with \*\*Direction\*\* = `vertical`")
        cv = _crit(bv)
        w.check((round(cv["Xo"], 6), round(cv["Yo"], 6)) == (10.0, 40.0),
                "the vertical face load's circle is not the unloaded (10, 40)")
        w.expect("normal against vertical (%)", 100 * (bn["results"]["FS"] / bv["results"]["FS"] - 1),
                 r"The two settings differ by (\d+)% here")

    with w.step("Finding the cohesion for a target factor of safety: Edit → Undo twice"):
        st.trigger("Edit", "Undo")
        st.trigger("Edit", "Undo")
        w.require(len(st.data.get("dloads") or []) == 1 and not st.data.get("line_loads"),
                  "two Undos leave dloads %s and line loads %s" % (st.data.get("dloads"),
                                                                    st.data.get("line_loads")))

    with w.step("Finding the cohesion for a target factor of safety: the design sweep"):
        vals = {k: p.number(r"\*\*%s\*\* `([\d.]+)`" % re.escape(k))[0]
                for k in ("From", "To", "Steps", "Target FS")}
        seen = {}

        def act(dlg):
            set_field(dlg, "Mode", "Design (FS target)")
            set_field(dlg, "Material", "soil")
            set_field(dlg, "Property", "c")
            seen["echo"] = value(field(dlg, "Sweeping"))
            for k, v in vals.items():
                set_field(dlg, k, v)
            seen["research"] = button(dlg, "Re-search the critical surface at each step").isChecked()
        st.run_parametric(act)
        w.check("mat:soil:c" in str(seen["echo"]), "Sweeping reads %r, not mat:soil:c" % seen["echo"])
        w.check(seen["research"], "Re-search the critical surface at each step is not ticked")
        d = st.results.get("design") or {}
        w.require(d.get("bracketed"), "the sweep did not bracket the target: %s" % d.get("message"))
        w.expect("cohesion for FS 1.5", d["crossing"], r"\*\*c = ([\d.]+) psf\.\*\*")
        w.expect("stronger than 500 psf (%)", 100 * (d["crossing"] / 500.0 - 1),
                 r"The soil would have to be (\d+)% stronger")
