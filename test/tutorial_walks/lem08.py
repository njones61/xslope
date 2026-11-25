"""LEM-8, A Reinforced Slope (Geogrids): built in Studio from File → New, then the
page's runs.

Walked: path C (Studio) — the two soils pasted, the two profile lines, the
crest surcharge, the six reinforcement lines (endpoint block pasted, Type set to
geosynthetic row by row and what it fills, capacity block pasted at Tmax), the
two circles — the Spencer search, the sliding mass and the tension it takes, the
crossing table under the solution plot, the critical circle pasted into Circles
and solved alone at the file's, zero and 10 ft pullout lengths (with the force
and the limiting end of each crossing), and everything put back as the page
says. Through the engine (the page reports the results, not the clicks): the
reinforcement removed, the length study and its 10 ft crossings, and Dir and
Appl switched. Not walked: path A (the assistant), path B (the template cell
by cell), the methods table and the UTEXASED circle (both held by the page's
tags), and the List view's five groups.
"""
import re

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, press, set_field, set_value,
                                    show_parameters_for, value)

PAGE = "lem08_reinforced_slope.md"
TIMEOUT = 3600


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _run(analysis="Auto search"):
    def act(dlg):
        set_field(dlg, "Method", "Spencer")
        set_field(dlg, "Analysis", analysis)
    return act


def _crossings(b, sd):
    from xslope.crossings import crossings
    found = crossings(b["slice_df"], sd, b["results"]["FS"])
    if found is None:
        raise StepFailed("the solution has no crossing record")
    return found


def _set_lp(st, value):
    def act(dlg):
        press(dlg, "Table view")
        t = editable(dlg)
        n = t.table.rowCount()
        paste(t, [[value, value]] * n, 0, column(t, "Lp1"))
    st.edit("Reinforcement lines", act)
    got = {(float(r["lp1"]), float(r["lp2"])) for r in st.data["reinforcement_lines"]}
    if got != {(float(value), float(value))}:
        raise StepFailed("Lp1 / Lp2 saved as %s, not %s" % (got, value))


def _first_circle(st, row):
    def act(dlg):
        t = editable(dlg)
        paste(t, [row], 0, column(t, "Xo"))
    st.edit("Circles", act)


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"so the table ends at φ:")
    line1 = page_rows(p.table(r"\*\*Profile Line 1 — material 1"))
    line2 = page_rows(p.table(r"\*\*Profile Line 2 — material 2"))
    load = page_rows(p.table(r"\*\*Surcharge\*\* — 240 psf over the crest"))
    ends = page_rows(p.table(r"Each line here starts on the face and runs 20 ft back"))
    caps = page_rows(p.table(r"with the same eight values on every line:"))
    circles = page_rows(p.table(r"\*\*Starting circles\*\* — two, with centers 40 ft up"))

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials"):
        n = p.count(r"press \*\*Add row\*\* (\w+), and enter \(or copy-paste\) the two materials")

        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(mats) - 1, "Add row %d times gives %d rows"
                      % (n, t.table.rowCount()))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", act)
        got = [(m["name"], float(m["gamma"]), float(m["c"]), float(m["phi"])) for m in st.data["materials"]]
        w.require(got == [(r[0], float(r[1]), float(r[4]), float(r[5])) for r in mats[1:]],
                  "the materials saved as %s" % got)

    with w.step("C2: Profile lines"):
        depth = p.number(r"Maximum depth = `(-?\d+)`, the elevation of the bottom")[0]

        def act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", depth)
            for pts, mat in ((line1, "1: shell"), (line2, "2: base")):
                press(dlg, "Add line")
                set_field(dlg, "Material", mat)
                paste(editable(dlg), pts)
        st.edit("Profile lines", act)
        got = [[tuple(map(float, xy)) for xy in l["coords"]] for l in st.data["profile_lines"]]
        w.require(got == [[tuple(map(float, r)) for r in pts] for pts in (line1, line2)],
                  "the profile lines saved as %s" % got)

    with w.step("C3: The surcharge"):
        seen = {}

        def act(dlg):
            press(dlg, "Add load")
            paste(editable(dlg), load)
            seen["dir"] = value(field(dlg, "Direction"))
        st.edit("Distributed loads", act)
        w.check(seen["dir"] == "Normal (perpendicular to the line)",
                "Direction reads %r" % seen["dir"])
        dl = st.data.get("dloads") or []
        w.require(len(dl) == 1 and [(q["X"], q["Y"], q["Normal"]) for q in dl[0]] ==
                  [tuple(map(float, r)) for r in load], "the surcharge saved as %s" % dl)
        w.require((st.data.get("dload_dirs") or ["normal"])[0] == "normal",
                  "the surcharge's direction is %s" % st.data.get("dload_dirs"))

    with w.step("C4: Reinforcement lines"):
        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            t = editable(dlg)
            paste(t, ends, 0, column(t, "Label"))
            w.require(t.table.rowCount() == len(ends), "the endpoint paste gives %d rows"
                      % t.table.rowCount())
            from PySide6.QtWidgets import QComboBox
            c = column(t, "Type")
            for r in range(len(ends)):
                cw = t.table.cellWidget(r, c)
                cb = cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox)
                set_value(cb, "geosynthetic")
            paste(t, caps, 0, column(t, "Tmax"))
        st.edit("Reinforcement lines", act)
        rl = st.data["reinforcement_lines"]
        w.require(len(rl) == len(ends), "%d lines saved" % len(rl))
        w.check(all((r.get("type"), r.get("dir"), r.get("appl")) == ("geosynthetic", "tangent", "active")
                    for r in rl), "Type geosynthetic does not fill Dir tangent and Appl active: %s"
                % [(r.get("type"), r.get("dir"), r.get("appl")) for r in rl])
        w.require(all(float(r["t_max"]) == 800 and float(r["lp1"]) == 4 and float(r["lp2"]) == 4
                      for r in rl), "the capacity block did not land on every line")

    with w.step("C5: Starting circles"):
        def act(dlg):
            paste(editable(dlg), circles)
        st.edit("Circles", act)
        got = [(float(c["Xo"]), float(c["Yo"]), float(c["Depth"])) for c in st.data["circles"]]
        w.require(got == [(float(r[0]), float(r[1]), float(r[3])) for r in circles],
                  "the circles saved as %s" % got)

    with w.step("What the search finds"):
        st.run_lem(_run())
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered")
        w.expect("center x", crit["Xo"], r"centered at \((−\d+\.\d+), \d+\.\d+\) with a radius")
        w.expect("center y", crit["Yo"], r"centered at \(−\d+\.\d+, (\d+\.\d+)\) with a radius")
        w.expect("radius", crit["Yo"] - crit["Depth"], r"with a radius of (\d+\.\d+)\. The solution")
        xs = [c[0] for c in b["failure_surface"].coords]
        w.expect("crest exit", max(xs), r"comes out on the crest at x ≈ (\d+)")
        w.expect("sliding mass", float(b["slice_df"]["w"].sum()), r"It weighs ([\d,]+) lb/ft")
        found = _crossings(b, st.data)
        tension = sum(float(r["applied"]) for r in found["rows"])
        w.expect("tension", tension, r"contributes ([\d,]+) lb/ft of tension")
        w.check(len(found["rows"]) == p.count(r"— (\w+) lines at the full 800"),
                "the circle crosses %d lines" % len(found["rows"]))
        searched = b

    with w.step("Where the surface meets the lines"):
        view = st.mw.solution_canvas
        title = view._crossings_title.text()
        w.check("Not crossed: Line 1." in title, "the note above the table does not list Line 1: %r" % title)
        table = view._crossings
        cells = [[table.item(r, c).text() if table.item(r, c) else "" for c in range(table.columnCount())]
                 for r in range(table.rowCount())]
        w.check(len(cells) == 5 and all(row[-1] == "Tmax" for row in cells),
                "the crossing table's last column does not read Tmax on five rows: %s" % cells)
        line2_row = next((r for r in found["rows"] if r["name"] == "Line 2"), None)
        w.require(line2_row is not None, "Line 2 is not crossed")
        w.expect("Line 2 crossing x", line2_row["x"], r"and is crossed at x = (\d+\.\d+), \d")
        w.expect("Line 2 distance from end 1", line2_row["x"] - 5.0, r"x = \d+\.\d+, (\d+\.\d+) ft from end 1")
        w.check(all(min(abs(r["x"] - ln["x1"]), abs(ln["x2"] - r["x"])) > 8.0
                    for r in found["rows"] for ln in [st.data["reinforcement_lines"][r["index"]]]),
                "a crossing lands within 8 ft of an end of its line")

    with w.step("Where the surface meets the lines: pullout lengths on the critical circle"):
        held = page_rows(p.table(r"paste the critical circle over the first row:"))[0]
        _first_circle(st, held)
        st.run_lem(_run("Single surface"))
        w.expect("single surface, file's lengths", _bundle(st)["results"]["FS"],
                 r"run it once as it stands: it returns the search's (1\.\d+)")
        _set_lp(st, "0")
        st.run_lem(_run("Single surface"))
        w.expect("single surface, Lp 0", _bundle(st)["results"]["FS"],
                 r"fully anchored, the answer stays \*\*(1\.\d+)\*\*")
        _set_lp(st, "10")
        st.run_lem(_run("Single surface"))
        b = _bundle(st)
        w.expect("single surface, Lp 10", b["results"]["FS"], r"drops to \*\*(1\.\d+)\*\*\. The table's")
        found = _crossings(b, st.data)
        w.expect("tension at Lp 10", sum(float(r["applied"]) for r in found["rows"]),
                 r"mobilize ([\d,]+) lb/ft between them")
        ends1 = {r["name"] for r in found["rows"] if "end 1" in r["limit"]}
        ends2 = {r["name"] for r in found["rows"] if "end 2" in r["limit"]}
        w.check(ends1 == {"Line 2", "Line 6"} and ends2 == {"Line 3", "Line 4", "Line 5"},
                "the limiting ends read end 1 %s, end 2 %s" % (sorted(ends1), sorted(ends2)))
        _set_lp(st, "4")
        _first_circle(st, circles[0])
        w.require((float(st.data["circles"][0]["Xo"]), float(st.data["circles"][0]["Depth"])) ==
                  (float(circles[0][0]), float(circles[0][3])), "the first circle did not go back")

    sd = dict(st.data)

    with w.step("Effect of the reinforcement (engine)"):
        bare = w.lem(dict(sd, reinforcement_lines=[]), "spencer")
        w.expect("unreinforced FS", bare["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* against 1\.589")
        w.expect("unreinforced mass", float(bare["slice_df"]["w"].sum()), r"sliver behind the face, ([\d,]+) lb/ft")
        crit = searched["search"]["fs_cache"][0]
        held_c = dict(sd["circles"][0], Xo=crit["Xo"], Yo=crit["Yo"], Depth=crit["Depth"],
                      R=crit["Yo"] - crit["Depth"])
        one = w.lem(dict(sd, reinforcement_lines=[], circles=[held_c]), "spencer",
                    analysis="single_surface")
        w.expect("unreinforced on the reinforced circle", one["results"]["FS"],
                 r"the unreinforced section returns (1\.\d+)\.")

    with w.step("Effect of the reinforcement length (engine)"):
        tbl = p.table(r"Every length gets its own search and its own critical circle:")
        for r in tbl[1:]:
            L = float(r[0])
            lines = [dict(ln, x2=ln["x1"] + L) for ln in sd["reinforcement_lines"]]
            model = dict(sd, reinforcement_lines=lines)
            b = w.lem(model, "spencer")
            found = _crossings(b, model)
            at = r"\| %s \| " % re.escape(r[0])
            w.expect("FS at %s ft" % r[0], b["results"]["FS"], at + r"(\d\.\d+) \|")
            w.expect("tension at %s ft" % r[0], sum(float(x["applied"]) for x in found["rows"]),
                     at + r"\d\.\d+ \| ([\d,]+) \|")
            w.expect("lines crossed at %s ft" % r[0], len(found["rows"]),
                     at + r"\d\.\d+ \| [\d,]+ \| (\d) of 6 \|")
            if L == 10:
                got = sorted(round(float(x["applied"])) for x in found["rows"])
                want = sorted(int(v) for v in p.find(r"where they mobilize (\d+) and (\d+) lb/ft").groups())
                w.check(got == want, "the 10 ft crossings mobilize %s, not %s" % (got, want))

    with w.step("Effect of Dir and Appl (engine)"):
        ax = w.lem(dict(sd, reinforcement_lines=[dict(ln, dir="axial") for ln in sd["reinforcement_lines"]]),
                   "spencer")
        w.expect("all Axial", ax["results"]["FS"], r"the search returns \*\*(1\.\d+)\*\*: the horizontal")
        pa = w.lem(dict(sd, reinforcement_lines=[dict(ln, appl="passive") for ln in sd["reinforcement_lines"]]),
                   "spencer")
        w.expect("all Passive", pa["results"]["FS"], r"the search returns \*\*(1\.\d+)\*\*\. Neither")
