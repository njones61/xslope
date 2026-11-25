"""LEM-3, A Layered Slope: built in Studio, then the page's runs on the file.

Walked: path C (Studio) from File → New — the two materials pasted from the
page's table, the two profile lines with the page's Add row counts, the
generated circles and their summary — and the Studio search from those three
circles. Then, as the page says, its figures and numbers come from the
completed file's two circles: the three starting circles solved alone, the
search and its warnings, the tension crack, the starts that all agree, the weak
foundation and the three force-equilibrium radii. Not walked: path A (the
assistant), path B (the template cell by cell; its workbook is the completed
file), and the far-out starting circle of "Guarding against local minima",
whose coordinates the page gives only roughly.
"""
import re

from tutorial_walks.harness import (StepFailed, column, editable, page_rows, paste,
                                    press, rows, set_field, show_parameters_for,
                                    type_cells)

PAGE = "lem03_layered_slope.md"
TIMEOUT = 1800
DONE = "docs/lem/files/xslope_simple_mult_layers.xlsx"


def _radius(c):
    return c["R"] if c.get("R") is not None else c["Yo"] - c["Depth"]


def walk(w):
    p = w.page
    st = w.studio()
    lit = p.lit

    mats = p.table(r"## The problem")
    line1 = p.table(lit("**Profile Line 1 — material 1"))
    line2 = p.table(lit("**Profile Line 2 — material 2"))

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials, Add row and paste the two materials"):
        n = p.count(r"Press \*\*Add row\*\* (\w+) and enter \(or copy-paste\) the two materials")

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            for _ in range(n):
                press(dlg, "Add row")
            w.require(t.table.rowCount() == len(mats) - 1,
                      "Add row %d time(s) gives %d rows for %d materials"
                      % (n, t.table.rowCount(), len(mats) - 1))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", act)
        got = [(m["name"], float(m["gamma"]), m["option"], float(m["c"]), float(m["phi"]))
               for m in st.data["materials"]]
        want = [(r[0], float(r[1]), r[3], float(r[4]), float(r[5])) for r in mats[1:]]
        w.require(got == want, "the materials saved as %s, not the page's %s" % (got, want))

    with w.step("C2: Profile lines, two lines with the page's Add row counts"):
        depth = p.number(r"set \*\*Max depth \(bottom boundary elevation\):\*\* to `(-?\d+)`")[0]
        n1 = p.count(r"`1: embankment`, and \*\*Add row\*\* (\w+)(?: times)? for")
        n2 = p.count(r"`2: foundation`, and \*\*Add row\*\* (\w+) for")

        def act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", depth)
            for mat, n_rows, tbl in (("1: embankment", n1, line1), ("2: foundation", n2, line2)):
                press(dlg, "Add line")
                set_field(dlg, "Material", mat)
                for _ in range(n_rows):
                    press(dlg, "Add row")
                t = editable(dlg)
                w.require(t.table.rowCount() == len(tbl) - 1,
                          "line %s: Add row %d time(s) gives %d rows for %d points"
                          % (mat, n_rows, t.table.rowCount(), len(tbl) - 1))
                type_cells(t, page_rows(tbl))
        st.edit("Profile lines", act)
        lines = st.data["profile_lines"]
        w.require([[tuple(map(float, xy)) for xy in l["coords"]] for l in lines] ==
                  [[tuple(map(float, r)) for r in page_rows(t)] for t in (line1, line2)],
                  "the profile lines saved as %s" % [l["coords"] for l in lines])

    with w.step("C3: Circles, Generate starting circles"):
        seen = {}

        def act(dlg):
            from PySide6.QtWidgets import QLabel
            press(dlg, "Generate starting circles")
            st.settle(10)
            seen["summary"] = [l.text() for l in dlg.findChildren(QLabel)
                               if l.text().startswith("Generated")]
        st.edit("Circles", act)
        quoted = p.find(r"says so: \*\"(Generated [^\"]+)\"\*").group(1)
        w.check(any(quoted in s for s in seen["summary"]),
                "the generator's summary reads %s, not the page's %r" % (seen["summary"], quoted))
        depths = p.table(r"Read the \*\*Depth\*\* column")
        got = [float(c["Depth"]) for c in st.data["circles"]]
        for i, r in enumerate(depths[1:]):
            if i < len(got):
                w.expect("generated circle %d depth" % (i + 1), got[i],
                         r"\| %d \| (−?[\d.]+) \|" % (i + 1))
        w.check(len(got) == len(depths) - 1, "the generator made %d circles" % len(got))

    with w.step("Running the analysis from Studio's three circles"):
        def act(dlg):
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
        st.run_lem(act)
        b = st.results.get("lem_solution")
        w.require(b and b.get("results"), "the run left no solution")
        crit = b["search"]["fs_cache"][0]
        w.expect("FS from Studio's circles", b["results"]["FS"],
                 r"reports the same \*\*(1\.\d+)\*\*")
        w.expect("x of the center from Studio's circles", crit["Xo"],
                 r"the search stops at \((\d+\.\d), \d+\.\d\)")
        w.expect("y of the center from Studio's circles", crit["Yo"],
                 r"the search stops at \(\d+\.\d, (\d+\.\d)\)")
        w.expect("radius from Studio's circles", _radius(crit),
                 r"with a radius of (\d+\.\d), and the warnings")

    sd = w.load(w.copy_model(DONE))

    with w.step("Solving the starting circles"):
        circles = list(sd["circles"])
        toe = dict(st.data["circles"][0])
        tbl = p.table(r"The factor of safety is the number on that run's solution plot")
        runs = {"0": circles[0], "−4.72": toe, "−10": circles[1]}
        for r in tbl[1:]:
            c = runs[r[1]]
            b = w.lem(dict(sd, circles=[c]), "spencer", analysis="single_surface")
            w.expect("%s FS" % r[0], b["results"]["FS"],
                     r"\| %s \| %s \| \**(1\.\d+)\**" % (re.escape(r[0]), re.escape(r[1])))
            w.expect("%s weight" % r[0], float(b["slice_df"]["w"].sum()),
                     r"\| %s \| %s \| \**1\.\d+\** \| ([\d,]+) \|" % (re.escape(r[0]), re.escape(r[1])))

    with w.step("The search result"):
        b = w.lem(sd, "spencer")
        crit = b["search"]["fs_cache"][0]
        w.expect("FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered")
        w.expect("x of the center", crit["Xo"], r"on a circle centered at \((\d+\.\d), ")
        w.expect("y of the center", crit["Yo"], r"on a circle centered at \(\d+\.\d, (\d+\.\d)\)")
        w.expect("radius", _radius(crit), r"\) with a radius of (\d+\.\d)\. Its")
        w.expect("lowest point", crit["Depth"], r"lowest point is at elevation (\d\.\d+)")
        warn = " ".join(str(x) for x in (b["results"].get("warnings") or []))
        w.say("warnings:", warn)
        df = b["slice_df"]
        w.check(len(df) == 40 and all(abs(float(c) - 400) < 1e-9 for c in df["c"]),
                "not every one of the 40 slice bases carries c = 400")
        for label, pat in (("most tensile interslice force", r"most tensile (−\d+) lb/ft"),
                           ("largest interslice compression", r"largest compression of (\d+)"),
                           ("thrust outside (%)", r"outside the slice on (\d+)% of the boundaries"),
                           ("slices with base tension", r"tension on the base of (\d+) slices near")):
            got = _from_warnings(label, warn, b)
            w.expect(label, got, pat)
        uns = b["search"]["unsolved"]
        w.check(uns["unsolved"] == 0, "Spencer left %d trial surfaces unsolved" % uns["unsolved"])
        w.expect("trial surfaces solved", uns["attempted"], r"Spencer answered on all (\d+) trial")
        zc = p.number(r"z<sub>c</sub> = 2c/γ = (\d+\.\d+) ft")[0]
        bc = w.lem(dict(sd, tcrack_depth=float(zc), tcrack_water=0.0), "spencer")
        w.expect("cracked FS", bc["results"]["FS"], r"clean solution at \*\*FS = (1\.\d+)\*\*")
        w.check(not bc["results"].get("warnings"), "the cracked run is not clean: %s"
                % bc["results"].get("warnings"))
        w.check(abs(bc["search"]["fs_cache"][0]["Depth"]) < 1e-6,
                "the cracked circle is not tangent to the contact (lowest point %.3f)"
                % bc["search"]["fs_cache"][0]["Depth"])

    with w.step("Guarding against local minima"):
        for c in list(sd["circles"]) + [toe]:
            bb = w.lem(dict(sd, circles=[c]), "spencer")
            w.expect("search from Depth %s" % c["Depth"], bb["results"]["FS"],
                     r"tangent to the contact, FS = (1\.\d+)\.")

    with w.step("A weak foundation"):
        weak = dict(sd, materials=[dict(m) for m in sd["materials"]])
        weak["materials"][1]["c"] = float(p.number(r"give the foundation c = (\d+) psf")[0])
        b = w.lem(weak, "spencer")
        w.expect("weak-foundation FS", b["results"]["FS"], r"\*\*FS = (0\.\d+)\*\*, and it is a different")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.expect("depth tangent to the rock", b["search"]["fs_cache"][0]["Depth"],
                 r"now tangent to the rock at elevation (−\d+)")
        w.expect("exit beyond the toe", -min(xs), r"it exits (\d+) ft beyond the toe")
        w.expect("reach into the crest", max(xs), r"reaches back to x = (\d+) in the crest")
        w.expect("weight", float(b["slice_df"]["w"].sum()), r"It weighs ([\d,]+) lb/ft")
        b0 = w.lem(sd, "spencer")
        w.expect("weight of the answer above", float(b0["slice_df"]["w"].sum()),
                 r"against the ([\d,]+) of the answer above")
        b1 = w.lem(dict(weak, circles=[sd["circles"][0]]), "spencer", analysis="single_surface")
        w.expect("contact circle solved alone", b1["results"]["FS"],
                 r"solved alone it returns the same (1\.\d+)")

    with w.step("Factor of safety by method: the force-equilibrium radii"):
        for meth, label in (("janbu", "Janbu"), ("lowe", "Lowe"), ("corps", "Corps")):
            b = w.lem(sd, meth)
            w.expect("%s radius" % label, _radius(b["search"]["fs_cache"][0]),
                     r"%s(?: R =)? (\d+\.\d)" % label)


def _from_warnings(label, text, bundle):
    """A number the warning strip prints."""
    t = text.replace("−", "-")
    if label.startswith("most tensile"):
        m = re.search(r"(-\d[\d,]*(?:\.\d+)?)", t)
    elif label.startswith("largest interslice"):
        m = re.search(r"compression[^\d-]*(\d[\d,]*(?:\.\d+)?)", t)
    elif label.startswith("thrust"):
        m = re.search(r"thrust outside the slice on (\d+(?:\.\d+)?)%", t)
    else:
        m = re.search(r"(\d+) slices? near the crest", t) or re.search(r"base of (\d+) slice", t)
    if not m:
        raise StepFailed("the warnings do not give the %s: %r" % (label, text))
    return float(m.group(1).replace(",", ""))
