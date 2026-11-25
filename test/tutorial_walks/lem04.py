"""LEM-4, Water in the Slope: built in Studio from File → New, then the page's runs.

Walked: path C (Studio) — the three soils pasted with LEM and FEM ticked, the
three profile lines, the piezometric line with the page's Add row count, the
circle — the search and the Log line that gives its circle, the circle held
and solved alone, the slice-base pore pressures, u set to none and the dry run
with the strength the water removes, both rows of the methods table (each a
single-surface run in Studio), the dry auto search, and the γsat cells cleared
and put back with Edit → Undo. Not walked: path A (the assistant), path B (the
template cell by cell; its workbook is the completed file the page's tags run),
the "around 1.37 or higher" of circles kept above soil 3 (the page gives no
run for it), the r_u section (prose only), and the factor-of-safety-by-method
tables (tags).
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, press, set_field, show_parameters_for,
                                    value)

PAGE = "lem04_water_in_the_slope.md"
TIMEOUT = 3600

METHODS = (("OMS", "Ordinary Method of Slices"), ("Bishop", "Bishop's Simplified"),
           ("Janbu", "Janbu"), ("Corps", "Corps of Engineers"), ("Lowe", "Lowe & Karafiath"),
           ("Spencer", "Spencer"), ("M-P", "Morgenstern-Price"))


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _run(method="Spencer", analysis="Auto search"):
    def act(dlg):
        set_field(dlg, "Method", method)
        set_field(dlg, "Analysis", analysis)
    return act


def _set_u(st, word):
    """Every material row's u, set as the page says (the u column, row by row)."""
    def act(dlg):
        press(dlg, "Table view")
        t = editable(dlg)
        paste(t, [[word]] * t.table.rowCount(), 0, column(t, "u"))
    st.edit("Materials", act)
    got = [m.get("u") for m in st.data["materials"]]
    if got != [word] * len(got):
        raise StepFailed("u reads %s, not %s on every row" % (got, word))


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"both reach the solver as zero:")
    lines = [page_rows(p.table(r"\*\*Profile Line %d — material %d" % (i, i))) for i in (1, 2, 3)]
    piezo = page_rows(p.table(r"\*\*Water\*\* — a piezometric line of eight points"))
    circle = page_rows(p.table(r"\*\*The starting circle\*\* — one seed"))

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials"):
        n = p.count(r"press \*\*Add row\*\* (\w+) times and enter \(or copy-paste\) the three soils")

        def act(dlg):
            show_parameters_for(dlg, {"LEM", "FEM"})
            press(dlg, "Table view")
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(mats) - 1, "Add row %d times gives %d rows"
                      % (n, t.table.rowCount()))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", act)
        got = [(m["name"], float(m["gamma"]), float(m["gamma_sat"]), float(m["c"]), float(m["phi"]), m["u"])
               for m in st.data["materials"]]
        want = [(r[0], float(r[1]), float(r[2]), float(r[4]), float(r[5]), r[13]) for r in mats[1:]]
        w.require(got == want, "the materials saved as %s" % got)

    with w.step("C1: Profile lines"):
        def act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", "0")
            for i, pts in enumerate(lines, 1):
                press(dlg, "Add line")
                set_field(dlg, "Material", "%d: soil %d" % (i, i))
                paste(editable(dlg), pts)
        st.edit("Profile lines", act)
        got = [[tuple(map(float, xy)) for xy in l["coords"]] for l in st.data["profile_lines"]]
        w.require(got == [[tuple(map(float, r)) for r in pts] for pts in lines],
                  "the profile lines saved as %s" % got)

    with w.step("C2: The piezometric line"):
        n = p.count(r"Press \*\*Add row\*\* (\w+) times and enter the points from the table above")
        seen = {}

        def act(dlg):
            seen["type"] = value(field(dlg, "Type"))
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(piezo), "Add row %d times gives %d rows for %d "
                      "points" % (n, t.table.rowCount(), len(piezo)))
            paste(t, piezo)
        st.edit("Piezometric lines", act)
        w.check(str(seen["type"]).startswith("piezo"), "Type reads %r, not piezo" % seen["type"])
        got = [tuple(map(float, xy)) for xy in st.data["piezo_line"]]
        w.require(got == [tuple(map(float, r)) for r in piezo], "the line saved as %s" % got)

    with w.step("C3: The circle"):
        def act(dlg):
            press(dlg, "Add row")
            paste(editable(dlg), circle)
        st.edit("Circles", act)
        c = st.data["circles"]
        w.require(len(c) == 1 and (float(c[0]["Xo"]), float(c[0]["Yo"]), float(c[0]["Depth"])) ==
                  (float(circle[0][0]), float(circle[0][1]), float(circle[0][3])),
                  "the circle saved as %s" % c)
        below = 40.0 - float(c[0]["Depth"])
        w.expect("circle bottom below the toe", below, r"bottoming out (\d+\.\d) ft below the toe")

    with w.step("What the search finds"):
        _s, log = st.run_lem(_run())
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, and the critical circle is not")
        for k, pat in (("Xo", r"\*\*Xo = (\d+\.\d+),"), ("Yo", r"Yo = (\d+\.\d+), Depth"),
                       ("Depth", r"Depth = (\d+\.\d+)\*\*")):
            w.expect("critical %s" % k, crit[k], pat)
        quoted = p.find(r"in Studio's \*\*Log\*\* pane, as `([^`]+)`").group(1)
        w.check(quoted in log, "the Log does not print %r" % quoted)
        held = re.search(r"x=([\d.]+), y=([\d.]+), depth=([\d.]+)", quoted).groups()

    with w.step("Holding the circle still"):
        def act(dlg):
            paste(editable(dlg), [[held[0], held[1], "Depth", held[2]]], 0,
                  column(editable(dlg), "Xo"))
        st.edit("Circles", act)
        st.run_lem(_run(analysis="Single surface"))
        b = _bundle(st)
        w.expect("held-circle FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* — the search's answer")
        df = b["slice_df"]
        u = df["u"].values
        w.expect("bases with pore pressure", int((u > 1e-9).sum()), r"on (\d+) of the 40 slice bases")
        k = int(np.argmax(u))
        w.expect("largest pore pressure", float(u[k]), r"the largest, \*\*(\d+) psf at x")
        w.expect("its x", float(df["x_c"].values[k]), r"psf at x = (\d+\.\d)\*\*")
        w.expect("its depth below the line", float(u[k]) / 62.4, r"a base sitting (\d+\.\d+) ft below")
        w.check(all(v == 0 for v in u[:3]), "the first three slices do not read zero")
        wet = b

    with w.step("Effect of the pore pressure"):
        _set_u(st, "none")
        st.run_lem(_run(analysis="Single surface"))
        dry = _bundle(st)
        w.expect("dry FS", dry["results"]["FS"], r"\*\*FS = (1\.\d+)\.\*\* Same circle")
        w.expect("base length", float(dry["slice_df"]["dl"].sum()), r"same (\d+\.\d) ft of base")
        w.expect("sliding mass", float(dry["slice_df"]["w"].sum()), r"same ([\d,]+) lb/ft of soil")
        w.expect("cost of the water (%)", 100 * (1 - wet["results"]["FS"] / dry["results"]["FS"]),
                 r"the pore pressure costs \*\*(\d+)% of the factor of safety\*\*")
        for label, bb in (("wet", wet), ("dry", dry)):
            dfx = bb["slice_df"]
            w.expect("c′Δℓ (%s)" % label, float((dfx["c"] * dfx["dl"]).sum()),
                     r"c′ Δℓ summed over the base is \*\*([\d,]+) lb/ft\*\*")
        dfw = wet["slice_df"]
        lost = float((dfw["u"] * dfw["dl"] * np.tan(np.radians(dfw["phi"].astype(float)))).sum())
        w.expect("friction the water removes", lost, r"\*\*([\d,]+) lb/ft\*\* of frictional resistance")

    with w.step("Effect of the pore pressure: the methods table"):
        at = p.find(r"the pattern holds everywhere:").end()
        for state, row_label in (("none", r"u = `none`"), ("piezo", "Piezometric line")):
            _set_u(st, state)
            for i, (head, meth) in enumerate(METHODS):
                st.run_lem(_run(meth, "Single surface"))
                fs = _bundle(st)["results"]["FS"]
                w.expect("%s, %s" % (head, state), fs,
                         r"\| %s \|" % row_label + r" [\d.]+ \|" * i + r" ([\d.]+) \|",
                         nth=0, start=at)

    with w.step("Effect of the pore pressure: the dry search"):
        _set_u(st, "none")
        st.run_lem(_run())
        w.expect("dry search", _bundle(st)["results"]["FS"], r"foundation — (1\.\d+)\. Isolating")
        _set_u(st, "piezo")

    with w.step("The two unit-weight columns"):
        def clear(dlg):
            press(dlg, "Table view")
            t = editable(dlg)
            c = column(t, "γsat")
            for r in range(t.table.rowCount()):
                it = t.table.item(r, c)
                if it is None:
                    raise StepFailed("the γsat cell of row %d cannot be cleared" % (r + 1))
                it.setText("")
        st.edit("Materials", clear)
        w.require(all(m.get("gamma_sat") in (None, "", 0, 0.0) for m in st.data["materials"]),
                  "γsat did not clear: %s" % [m.get("gamma_sat") for m in st.data["materials"]])
        st.run_lem(_run(analysis="Single surface"))
        b = _bundle(st)
        w.expect("γ only, weight", float(b["slice_df"]["w"].sum()), r"\| γ only \(gsat blank\) \| ([\d,]+) \|")
        w.expect("γ only, FS", b["results"]["FS"], r"\| γ only \(gsat blank\) \| [\d,]+ \| (1\.\d+) \|")
        w.expect("split, weight", float(wet["slice_df"]["w"].sum()), r"γ_sat below \| ([\d,]+) \|")
        w.expect("split, FS", wet["results"]["FS"], r"γ_sat below \| [\d,]+ \| (1\.\d+) \|")
        w.check(np.allclose(b["slice_df"]["u"].values, wet["slice_df"]["u"].values),
                "clearing γsat changed u")
        st.trigger("Edit", "Undo")
        w.require([float(m.get("gamma_sat") or 0) for m in st.data["materials"]] ==
                  [float(r[2]) for r in mats[1:]], "Undo did not put the γsat cells back")
