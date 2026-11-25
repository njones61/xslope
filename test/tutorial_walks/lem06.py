"""LEM-6, Polygon Geometry: built in Studio from File → New, then the page's runs.

Walked: path C (Studio) — the materials pasted, the two polygon rings pasted
from the page's tables, the generated circles and the deep seed typed in — the
Studio search, the solution's slices and weight, the deep circle alone, the
circle that will not fit (refused, then truncated as a composite surface), the
composite search, and the weak foundation with and without composite surfaces.
Not walked: path A (the assistant) and path B (the template cell by cell, whose
workbook is the completed file the page's tags run).
"""
from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, press, rows, set_field, show_parameters_for,
                                    tick, type_cells)

PAGE = "lem06_polygon_geometry.md"
TIMEOUT = 1800


def _radius(c):
    return c["R"] if c.get("R") is not None else c["Yo"] - c["Depth"]


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"so the table ends at φ:")
    ring1 = page_rows(p.table(r"\*\*Polygon 1 — material 1"))
    ring2 = page_rows(p.table(r"\*\*Polygon 2 — material 2"))
    circles = p.table(r"\*\*Starting circles\*\* — two, sharing a center")

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials"):
        n = p.count(r"press \*\*Add row\*\* (\w+) and enter \(or copy-paste\) the two materials")

        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            for _ in range(n):
                press(dlg, "Add row")
            w.require(t.table.rowCount() == len(mats) - 1, "Add row %d time(s) gives %d rows"
                      % (n, t.table.rowCount()))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", act)
        got = [(m["name"], float(m["gamma"]), float(m["c"])) for m in st.data["materials"]]
        w.require(got == [(r[0], float(r[1]), float(r[4])) for r in mats[1:]],
                  "the materials saved as %s" % got)

    with w.step("C2: Polygons"):
        def act(dlg):
            for ring, mat in ((ring1, "1: embankment"), (ring2, "2: foundation")):
                press(dlg, "Add polygon")
                t = editable(dlg)
                w.require(t.table.rowCount() == 0, "a new polygon starts with %d row(s)" % t.table.rowCount())
                paste(t, ring)
                set_field(dlg, "Material", mat)
                w.check(field(dlg, "Type").currentText() == "material", "Type does not read material")
        st.edit("Polygons", act)
        polys = st.data.get("polygons") or []
        w.require(len(polys) == 2, "%d polygons saved" % len(polys))
        for pg, ring in zip(polys, (ring1, ring2)):
            got = {tuple(round(float(v), 9) for v in xy) for xy in pg["polygon"].exterior.coords}
            w.require(got == {tuple(float(v) for v in r) for r in ring},
                      "a ring saved as %s" % sorted(got))
        w.check(not st.data.get("profile_lines"), "the model has profile lines")

    with w.step("C3: Starting circles"):
        deep = p.number(r"change the toe circle's \*\*Depth\*\* to `(-?[\d.]+)`")[0]

        def act(dlg):
            press(dlg, "Generate starting circles")
            st.settle(10)
            t = editable(dlg)
            got = rows(t)
            w.require(len(got) == 2, "the generator proposes %d circles, not two" % len(got))
            type_cells(t, [[deep]], 0, column(t, "Depth"))
        st.edit("Circles", act)
        c = st.data["circles"]
        want = [(float(r[0]), float(r[1]), float(r[3])) for r in circles[1:]]
        got = [(float(x["Xo"]), float(x["Yo"]), float(x["Depth"])) for x in c]
        w.check(got == want, "the circles table reads %s, not the page's %s" % (got, want))

    with w.step("What the search finds"):
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"),
                                set_field(dlg, "Analysis", "Auto search")))
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered")
        w.expect("center x", crit["Xo"], r"centered at \((\d+\.\d+), \d+\.\d+\) with a radius")
        w.expect("center y", crit["Yo"], r"centered at \(\d+\.\d+, (\d+\.\d+)\) with a radius")
        w.expect("radius", _radius(crit), r"with a radius of (\d+\.\d+)\. Notice")
        w.expect("lowest point", crit["Depth"], r"Its lowest point is at elevation (\d\.\d+)")
        df = b["slice_df"]
        w.check(len(df) == 40 and all(abs(float(v) - 400) < 1e-9 for v in df["c"]),
                "not all 40 slice bases carry c = 400")
        w.expect("sliding mass", float(df["w"].sum()), r"the sliding mass weighs ([\d,]+) lb/ft")
        sd = dict(st.data)
        one = w.lem(dict(sd, circles=[sd["circles"][0]]), "spencer", analysis="single_surface")
        w.expect("deep circle alone", one["results"]["FS"], r"comes to (1\.\d+), carrying")
        w.expect("deep circle weight", float(one["slice_df"]["w"].sum()), r"carrying ([\d,]+) lb/ft against")

    with w.step("Circles that will not fit"):
        new_d = p.number(r"\*\*Depth\*\* `(-?\d+)` in place of")[0]
        st.edit("Circles", lambda dlg: type_cells(editable(dlg), [[new_d]], 0, column(editable(dlg), "Depth")))
        mark = st.mark()
        st.run_lem(lambda dlg: set_field(dlg, "Analysis", "Single surface"))
        quoted = p.find(r"> \*\*LEM run failed\*\* — (.+?) The run is refused").group(1).strip()
        said = st.log_since(mark) + " ".join(t for _k, _ti, t in st.dialogs.boxes) + st.status()
        w.check(quoted in " ".join(said.split()), "the refusal does not read %r" % quoted)
        st.run_lem(lambda dlg: (set_field(dlg, "Analysis", "Single surface"),
                                tick(dlg, "Composite surfaces (truncate circles at the base)", True)))
        b = st.results.get("lem_solution")
        w.require(b and b.get("results"), "the composite run left no solution")
        w.expect("composite FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\.\*\* The surface follows")
        df = b["slice_df"]
        on_base = [i for i in range(len(df)) if abs(abs(float(df["alpha"].iloc[i])) - 3.366) < 0.05]
        w.expect("slices on the base", len(on_base), r"— (\d+) of its \d+ slice bases sitting")
        w.expect("slices", len(df), r"— \d+ of its (\d+) slice bases sitting")
        if on_base:
            x0 = float(df["x_l"].iloc[on_base[0]]) if "x_l" in df else None
            x1 = float(df["x_r"].iloc[on_base[-1]]) if "x_r" in df else None
            if x0 is not None:
                w.expect("composite meets the bedrock", x0, r"meets the bedrock at x = (\d+\.\d+)")
                w.expect("composite leaves the bedrock", x1, r"for [\d.]+ ft to x = (\d+\.\d+)")
        st.run_lem(lambda dlg: (set_field(dlg, "Analysis", "Auto search"),
                                tick(dlg, "Composite surfaces (truncate circles at the base)", True)))
        w.expect("composite search", st.results["lem_solution"]["results"]["FS"],
                 r"box ticked returns \*\*(1\.\d+)\*\*")

    with w.step("A weak foundation"):
        c_weak = p.number(r"give the foundation c = (\d+) psf")[0]

        def mats_act(dlg):
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            type_cells(t, [[c_weak]], 1, column(t, "c"))
        st.edit("Materials", mats_act)
        st.edit("Circles", lambda dlg: type_cells(editable(dlg), [[deep]], 0, column(editable(dlg), "Depth")))
        st.run_lem(lambda dlg: (set_field(dlg, "Analysis", "Auto search"),
                                tick(dlg, "Composite surfaces (truncate circles at the base)", False)))
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("weak FS", b["results"]["FS"], r"\*\*FS = (0\.\d+)\*\*, and the mechanism")
        w.expect("weak lowest point", crit["Depth"], r"bottoms out at elevation (−\d+\.\d+)")
        w.expect("weak weight", float(b["slice_df"]["w"].sum()), r"with ([\d,]+) lb/ft of soil")
        st.run_lem(lambda dlg: (set_field(dlg, "Analysis", "Auto search"),
                                tick(dlg, "Composite surfaces (truncate circles at the base)", True)))
        w.expect("weak composite FS", st.results["lem_solution"]["results"]["FS"],
                 r"composite surfaces on: it returns (0\.\d+)")
