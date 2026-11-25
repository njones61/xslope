"""LEM-10, Finding the Global Minimum: both parts in Studio, on copies of the
page's files.

Walked: Part A — the file's circle searched and solved alone, the embankment
circle pasted over it and searched, the infinite-slope check, the surficial
filter at 5 and at 8 ft, and both circles on the sheet in either order. Part B —
the James Bay dike from its single seed with the 2 m filter, with grid search,
and from the generated starting circles. Not walked: the factor-of-safety-by-
method table (the page's tags), and the trial-circle fan of the search plot.
"""
import math

from tutorial_walks.harness import (StepFailed, column, editable, page_rows, paste,
                                    press, set_field, tick)

PAGE = "lem10_global_minimum.md"
TIMEOUT = 3600


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _run(analysis="Auto search", skin=None, grid=None):
    def act(dlg):
        set_field(dlg, "Method", "Spencer")
        set_field(dlg, "Analysis", analysis)
        if skin is not None:
            tick(dlg, "Ignore surficial (skin) failures", skin is not False)
            if skin is not False:
                set_field(dlg, "Min slip depth", skin)
        if grid is not None:
            tick(dlg, "Grid search (auto-seed the circular search)", grid)
    return act


def _circles(st, block):
    def act(dlg):
        t = editable(dlg)
        t.table.selectAll()
        press(dlg, "Remove selected")
        paste(t, block, 0, column(t, "Xo"))
    st.edit("Circles", act)
    got = [(float(c["Xo"]), float(c["Yo"]), float(c["Depth"])) for c in st.data["circles"]]
    want = [(float(r[0]), float(r[1]), float(r[3])) for r in block]
    if got != want:
        raise StepFailed("the circles read %s, not %s" % (got, want))


def walk(w):
    p = w.page
    st = w.studio()

    # ------------------------------------------------------------------ Part A
    with w.step("Part A: opening the model"):
        st.open(w.copy_model(w.page_link("xslope_mult_min_KEY.xlsx")))
        deep = [[str(c["Xo"]), str(c["Yo"]), "Depth", str(c["Depth"])] for c in st.data["circles"]]
        w.require(len(deep) == 1, "the file has %d circles, not one" % len(deep))
        w.check(float(deep[0][3]) == float(st.data["max_depth"]),
                "the file's circle does not reach the bottom of the clay")

    with w.step("Part A: running the analysis"):
        st.run_lem(_run())
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("deep FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered at \(16")
        w.expect("center x", crit["Xo"], r"centered at \((\d+\.\d+), \d+\.\d+\) with a radius of \d+\.\d+ —")
        w.expect("center y", crit["Yo"], r"centered at \(\d+\.\d+, (\d+\.\d+)\) with a radius of \d+\.\d+ —")
        w.expect("radius", crit["Yo"] - crit["Depth"], r"with a radius of (\d+\.\d+) —")
        w.expect("tangent depth", crit["Depth"], r"tangent to the limiting depth at elevation (−\d+)")
        # The toe is the embankment line's lowest vertex; the circle leaves the
        # ground on the far side of it from the crest.
        line1 = st.data["profile_lines"][0]["coords"]
        toe_x = min(line1, key=lambda xy: xy[1])[0]
        xs = [c[0] for c in b["failure_surface"].coords]
        w.expect("exit beyond the toe", abs(toe_x - min(xs)), r"exiting (\d+) ft beyond the toe")
        df = b["slice_df"]
        w.expect("surface length", float(df["dl"].sum()), r"(\d+\.\d) ft of failure surface carrying")
        w.expect("sliding mass", float(df["w"].sum()), r"carrying ([\d,]+) lb/ft of soil")
        one = w.lem(dict(st.data), "spencer", analysis="single_surface")
        w.expect("the starting circle alone", one["results"]["FS"], r"solved as entered it gives (1\.\d+)")

    emb = page_rows(p.table(r"tangent to the top of the foundation instead of the bottom"))

    with w.step("Part A: changing the starting circle"):
        _circles(st, emb)
        st.run_lem(_run())
        b = _bundle(st)
        w.expect("sliver FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\.\*\* The gray arcs")
        df = b["slice_df"]
        w.expect("sliver length", float(df["dl"].sum()), r"high on the slope — (\d\.\d) ft long")
        w.expect("sliver mass", float(df["w"].sum()), r"moving (\d\.\d) lb/ft of sand")
        sand = st.data["materials"][0]
        beta = math.degrees(math.atan(1 / 2.25))
        inf = math.tan(math.radians(float(sand["phi"]))) / math.tan(math.radians(beta))
        w.expect("infinite-slope factor of safety", inf, r"= \*\*(1\.\d+)\*\* for this 2\.25:1 face")
        w.expect("face angle", beta, r"tan 30° / tan (\d+\.\d+)° =")

    with w.step("Part A: the surficial filter"):
        d5 = p.number(r"tick the box and set the depth to `(\d+)`")[0]
        st.run_lem(_run(skin=d5))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("filtered FS", b["results"]["FS"], r"the search now returns \*\*(1\.\d+)\*\*, on a circle")
        w.expect("filtered radius", crit["Yo"] - crit["Depth"], r"on a circle of radius (\d+\.\d) ft")
        d8 = p.number(r"Only at (\d+) ft does the embankment seed reach")[0]
        st.run_lem(_run(skin=d8))
        w.expect("filtered at 8 ft", _bundle(st)["results"]["FS"], r"does the embankment seed reach (1\.\d+)\.")

    with w.step("Part A: searching from both circles at once"):
        for block in (deep + emb, emb + deep):
            _circles(st, block)
            st.run_lem(_run(skin=False))
            w.expect("both circles, %s first" % ("deep" if block[0] == deep[0] else "embankment"),
                     _bundle(st)["results"]["FS"], r"follows the deep one and returns \*\*(1\.\d+)\*\*")

    # ------------------------------------------------------------------ Part B
    with w.step("Part B: running the analysis from the single seed"):
        st.open(w.copy_model(w.page_link("xslope_james_bay.xlsx")))
        w.require(len(st.data["circles"]) == 1, "the dike file has %d circles" % len(st.data["circles"]))
        d2 = p.number(r"\*\*Min slip depth\*\* = `(\d+)` — the metric")[0]
        st.run_lem(_run(skin=d2))
        b = _bundle(st)
        w.expect("single seed", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, with no warnings")
        w.check(not b["results"].get("warnings"), "the run has warnings: %s" % b["results"].get("warnings"))
        w.check(abs(b["search"]["fs_cache"][0]["Depth"] - float(st.data["max_depth"])) < 1e-6,
                "the circle does not bottom on the base of the model")
        single = b["results"]["FS"]

    with w.step("Part B: grid search"):
        st.run_lem(_run(skin=d2, grid=True))
        b = _bundle(st)
        w.expect("grid search", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a wider circle")
        w.expect("below the single seed (%)", 100 * (1 - b["results"]["FS"] / single),
                 r"daylights beyond it — (\d+)% below")
        w.expect("single seed above (%)", 100 * (single / b["results"]["FS"] - 1),
                 r"single reasonable starting circle converged (\d+)% high")

        def gen(dlg):
            press(dlg, "Generate starting circles")
            st.settle(10)

        def replace(box):
            # "Replace what is already here?": the page reads the full generated set.
            from PySide6.QtWidgets import QMessageBox
            box.button(QMessageBox.Apply).click()
            return "keep"
        st.edit("Circles", gen, replace)
        st.run_lem(_run(skin=d2, grid=False))
        w.expect("generated circles", _bundle(st)["results"]["FS"],
                 r"set also finds (1\.\d+) with the")
