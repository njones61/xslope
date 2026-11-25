"""LEM-9, A Tieback Wall: built in Studio from File → New, then the page's runs.

Walked: path C (Studio) — the two soils, the two profile lines, the wedge's
three points, the two tiebacks (endpoint block pasted, Type set to tieback and
what it fills, capacity block pasted at Tmax), the soldier pile — the Janbu
search and where it left the wedge, the sliding mass, the two crossings and the
force they deliver, the pile's line in the Log, the circular comparison from a
generated starting circle (Janbu, Spencer and Bishop), the manual's wedge solved
alone at 60 slices, and the anchors removed and searched. Through the engine (the
page reports the results, not the clicks): the two other rows of the
support table, and VP85's active and passive cases. Not walked: path A (the
assistant), path B (the template cell by cell), the five-method table (tags),
and the List view's groups.
"""
import math

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, press, set_field, set_value,
                                    show_parameters_for, value)

PAGE = "lem09_tieback_wall.md"
TIMEOUT = 3600


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _run(method="Janbu", analysis="Auto search", surface=None, slices=None):
    def act(dlg):
        set_field(dlg, "Method", method)
        set_field(dlg, "Analysis", analysis)
        if surface is not None:
            set_field(dlg, "Surface", surface)
        if slices is not None:
            set_field(dlg, "Number of slices", slices)
    return act


def _crossings(b, sd):
    from xslope.crossings import crossings
    found = crossings(b["slice_df"], sd, b["results"]["FS"])
    if found is None:
        raise StepFailed("the solution has no crossing record")
    return found


def _set_combo(t, row, header, word):
    from PySide6.QtWidgets import QComboBox
    cw = t.table.cellWidget(row, column(t, header))
    cb = cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox)
    set_value(cb, word)


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"so the table ends at φ:")
    line1 = page_rows(p.table(r"\*\*Profile Line 1 — material 1"))
    line2 = page_rows(p.table(r"\*\*Profile Line 2 — material 2"))
    ends = page_rows(p.table(r"Both are declined 25°:"))
    caps = page_rows(p.table(r"capacity values below are their own table"))
    pile = page_rows(p.table(r"because a pile resists in shear and bending"))
    wedge = page_rows(p.table(r"entered as three points:"))

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials and profile lines"):
        def m_act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            for _ in range(len(mats) - 1):
                press(dlg, "Add row")
            t = editable(dlg)
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", m_act)
        got = [(m["name"], float(m["gamma"]), float(m["gamma_sat"]), float(m["c"]), float(m["phi"]))
               for m in st.data["materials"]]
        w.require(got == [(r[0], float(r[1]), float(r[2]), float(r[4]), float(r[5])) for r in mats[1:]],
                  "the materials saved as %s" % got)
        depth = p.number(r"Maximum depth = `(-?\d+)`, the elevation of the bottom")[0]

        def p_act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", depth)
            for pts, mat in ((line1, "1: Layer 1"), (line2, "2: Layer 2")):
                press(dlg, "Add line")
                set_field(dlg, "Material", mat)
                paste(editable(dlg), pts)
        st.edit("Profile lines", p_act)
        got = [[tuple(map(float, xy)) for xy in l["coords"]] for l in st.data["profile_lines"]]
        w.require(got == [[tuple(map(float, r)) for r in pts] for pts in (line1, line2)],
                  "the profile lines saved as %s" % got)

    with w.step("C2: The failure surface"):
        st.edit("Non-circular pts", lambda dlg: paste(editable(dlg), wedge))
        nc = st.data.get("non_circ") or []
        w.require([(q["X"], q["Y"], q["Movement"]) for q in nc] ==
                  [(float(r[0]), float(r[1]), r[2]) for r in wedge], "the wedge saved as %s" % nc)

    with w.step("C3: Reinforcement lines"):
        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            t = editable(dlg)
            paste(t, ends, 0, column(t, "Label"))
            w.require(t.table.rowCount() == len(ends), "the endpoint paste gives %d rows"
                      % t.table.rowCount())
            for r in range(len(ends)):
                _set_combo(t, r, "Type", "tieback")
            paste(t, caps, 0, column(t, "Tmax"))
        st.edit("Reinforcement lines", act)
        rl = st.data["reinforcement_lines"]
        w.check(all((r.get("type"), r.get("dir"), r.get("appl")) == ("tieback", "axial", "active")
                    for r in rl), "Type tieback does not fill Dir axial and Appl active: %s"
                % [(r.get("type"), r.get("dir"), r.get("appl")) for r in rl])
        w.require([(float(r["t_max"]), float(r["lp1"]), float(r["lp2"])) for r in rl] ==
                  [(float(c[0]), float(c[1]), float(c[2])) for c in caps],
                  "the capacities saved as %s" % [(r["t_max"], r["lp1"], r["lp2"]) for r in rl])

    with w.step("C4: The soldier pile"):
        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            t = editable(dlg)
            paste(t, pile, 0, column(t, "Label"))
        st.edit("Piles", act)
        pl = st.data.get("pile_lines") or []
        w.require(len(pl) == 1 and float(pl[0]["H"]) == float(pile[0][5]) and
                  str(pl[0].get("appl")).lower() == pile[0][6].lower(),
                  "the pile saved as %s" % pl)

    with w.step("What the search finds"):
        seen = {}

        def act(dlg):
            seen["surface"] = value(field(dlg, "Surface"))
            _run()(dlg)
        _s, log = st.run_lem(act)
        w.check(seen["surface"] == "Non-circular", "Surface reads %r" % seen["surface"])
        b = _bundle(st)
        w.expect("Janbu FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a surface that runs from")
        pts = list(b["failure_surface"].coords)
        if len(pts) != 3:
            raise StepFailed("the searched surface has %d vertices, not the wedge's three" % len(pts))
        w.expect("middle x", pts[1][0], r"through \((\d+\.\d+), \d+\.\d+\) and out")
        w.expect("middle y", pts[1][1], r"through \(\d+\.\d+, (\d+\.\d+)\) and out")
        w.expect("exit x", pts[2][0], r"and out at \((\d+\.\d+), \d+\.\d+\)\.")
        w.expect("exit y", pts[2][1], r"and out at \(\d+\.\d+, (\d+\.\d+)\)\.")
        low = math.degrees(math.atan2(pts[1][1] - pts[0][1], pts[1][0] - pts[0][0]))
        up = math.degrees(math.atan2(pts[2][1] - pts[1][1], pts[2][0] - pts[1][0]))
        w.expect("lower leg", low, r"the lower leg rises at (\d+\.\d)°")
        w.expect("upper leg", up, r"and the upper at (\d+\.\d)°")
        w.expect("exit pulled closer", float(wedge[-1][0]) - pts[2][0], r"pulled the exit point (\d+) ft closer")
        df = b["slice_df"]
        w.expect("sliding mass", float(df["w"].sum()), r"The sliding mass weighs ([\d,]+) lb/ft")
        w.expect("base", float(df["dl"].sum()), r"lb/ft on an (\d+\.\d) ft base")
        found = _crossings(b, st.data)
        by_y = {round(st.data["reinforcement_lines"][r["index"]]["y1"]): r for r in found["rows"]}
        w.require(set(by_y) == {8, 20}, "the surface does not cross both anchors: %s" % sorted(by_y))
        w.expect("lower anchor crossing", by_y[8]["x"], r"crosses the lower anchor (\d+\.\d) ft out")
        w.expect("upper anchor crossing", by_y[20]["x"], r"and the upper one (\d+\.\d) ft out")
        total = sum(float(r["applied"]) for r in found["rows"])
        w.expect("anchor force", total, r"deliver their full capacity — ([\d,]+) lb/ft between them")
        rad = math.radians(25.0)
        w.expect("holding back", total * math.cos(rad), r"resolves into ([\d,]+) lb/ft holding")
        w.expect("pressing down", total * math.sin(rad), r"lb/ft holding the mass back and ([\d,]+) lb/ft pressing")
        for r in found["rows"]:
            ln = st.data["reinforcement_lines"][r["index"]]
            behind = math.hypot(ln["x2"] - r["x"], ln["y2"] - r["y"])
            w.check(behind > 20.0, "only %.1f ft of %s lies behind the surface" % (behind, r["name"]))
        quoted = "H stated = 5,900.0 lb/ft"
        w.check(quoted in log, "the Log does not give the pile's %r" % quoted)
        searched = b

    with w.step("Comparison with a circular search"):
        def act(dlg):
            press(dlg, "Generate starting circles")
            st.settle(10)
        st.edit("Circles", act)
        w.check(len(st.data["circles"]) == 1, "the generator made %d circles, not one"
                % len(st.data["circles"]))
        for meth, pat in (("Janbu", r"It gives \*\*FS = (1\.\d+), 3% above"),
                          ("Spencer", r"Spencer (1\.\d+) against 1\.414"),
                          ("Bishop's Simplified", r"a circle restores its moment arm — (1\.\d+)\.")):
            st.run_lem(_run(meth, surface="Circular"))
            b = _bundle(st)
            w.expect("circular %s" % meth, b["results"]["FS"], pat)
            if meth == "Janbu":
                pts = list(b["failure_surface"].coords)
                entry = min(pts, key=lambda q: q[0])
                w.expect("circle leaves the wall face at", entry[1], r"It leaves the wall face at elevation (\d+),")
                crossed = {round(st.data["reinforcement_lines"][r["index"]]["y1"])
                           for r in _crossings(b, st.data)["rows"]}
                w.check(crossed == {20}, "the circle crosses the anchors at %s, not only the upper"
                        % sorted(crossed))

    with w.step("The wedge the reference manual gives"):
        n = p.number(r"and raise the slice count to (\d+):")[0]
        _s, log = st.run_lem(_run(analysis="Single surface", surface="Non-circular", slices=n))
        b = _bundle(st)
        r = b["results"]
        w.expect("wedge, Janbu", r["FS"], r"\*\*FS = (1\.\d+)\*\* with Janbu's correction factor")
        w.expect("correction factor", r.get("fo"), r"correction factor f₀ = (\d\.\d+) applied")
        w.expect("wedge, uncorrected", r.get("FS_base"), r"applied, (1\.\d+) without it")
        df = b["slice_df"]
        w.expect("wedge mass", float(df["w"].sum()), r"The wedge carries ([\d,]+) lb/ft on a")
        w.expect("wedge base", float(df["dl"].sum()), r"lb/ft on a (\d+\.\d) ft base, \d+%")
        w.expect("more soil (%)", 100 * (float(df["w"].sum()) / float(searched["slice_df"]["w"].sum()) - 1),
                 r"ft base, (\d+)% more soil than")
        pts = list(b["failure_surface"].coords)
        w.expect("wedge lower leg", math.degrees(math.atan2(pts[1][1] - pts[0][1], pts[1][0] - pts[0][0])),
                 r"leaves the toe at (\d+)° against")
        st.run_lem(_run("Spencer", analysis="Single surface", surface="Non-circular", slices=n))
        w.expect("wedge, Spencer", _bundle(st)["results"]["FS"], r"and (1\.\d+) with Spencer's method")

    with w.step("Effect of the tiebacks: the anchors removed"):
        def act(dlg):
            press(dlg, "Table view")
            t = editable(dlg)
            t.table.selectAll()
            press(dlg, "Remove selected")
        st.edit("Reinforcement lines", act)
        w.require(not st.data.get("reinforcement_lines"), "the anchors are still in the model")
        st.run_lem(_run(surface="Non-circular", slices="40"))
        b = _bundle(st)
        w.expect("no anchors", b["results"]["FS"], r"\*\*FS = (0\.\d+)\*\* against 1\.431")
        ex = list(b["failure_surface"].coords)[-1]
        w.expect("exit x", ex[0], r"exits at \((\d+\.\d+), \d+\.\d+\) on the cut slope")
        w.expect("exit y", ex[1], r"exits at \(\d+\.\d+, (\d+\.\d+)\) on the cut slope")
        w.expect("mass", float(b["slice_df"]["w"].sum()), r"cuts ([\d,]+) lb/ft of soil against")
        w.expect("pile only (table)", b["results"]["FS"], r"\| Pile only \| (\d\.\d+) \|")

    with w.step("Effect of the tiebacks: the support table (engine)"):
        st.trigger("Edit", "Undo")
        sd = dict(st.data)
        w.require(len(sd.get("reinforcement_lines") or []) == 2, "Undo did not put the anchors back")
        only = w.lem(dict(sd, pile_lines=[]), "janbu", surface="noncircular")
        w.expect("anchors only", only["results"]["FS"], r"\| Anchors only \| (\d\.\d+) \|")
        neither = w.lem(dict(sd, pile_lines=[], reinforcement_lines=[]), "janbu", surface="noncircular")
        w.expect("neither", neither["results"]["FS"], r"\| Neither \| (\d\.\d+) \|")

    with w.step("Active or passive (engine)"):
        a = w.load(w.copy_model("docs/verification/files/rocscience/vp085a.xlsx", "vp"))
        bb = w.load(w.copy_model("docs/verification/files/rocscience/vp085b.xlsx", "vp"))
        ra = w.lem(a, "bishop", analysis="single_surface", num_slices=60)
        w.expect("vp085a, active", ra["results"]["FS"], r"\| `vp085a` \| Active \| (\d\.\d+) \|")
        rb = w.lem(bb, "bishop", analysis="single_surface", num_slices=60)
        w.expect("vp085b, passive", rb["results"]["FS"], r"\| `vp085b` \| Passive \| (\d\.\d+) \|")
        lines = a["reinforcement_lines"]
        pas = w.lem(dict(a, reinforcement_lines=[dict(ln, appl="passive") for ln in lines]),
                    "bishop", analysis="single_surface", num_slices=60)
        w.expect("vp085a circle, passive", pas["results"]["FS"], r"active gives 1\.\d+ and passive (\d\.\d+), against")
        none = w.lem(dict(a, reinforcement_lines=[]), "bishop", analysis="single_surface", num_slices=60)
        w.expect("vp085a circle, no anchor", none["results"]["FS"], r"against (\d\.\d+) with no\s?anchor")
