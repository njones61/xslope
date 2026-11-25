"""LEM-12, Piles: the pile rows read in Studio, the searches run there, and the
page's fixed-surface tables through the engine.

Walked in Studio, on a copy of the linked workbook: the Piles editor (Table view
on LEM alone, its eleven columns and two rows, then List view's four groups and
the per-shaft capacity labels); Run LEM's Spencer search, its Log block and the
critical surface's numbers; the file's circle as a single surface; the spacing
set to 12 and searched; the spacing set to 3 (the Model checks warning and the
Log warning); the force typed into H and searched; H cleared and the grid search,
then the grid search with a minimum slip depth. Through the engine (the page
reports these as fixed-surface studies a reader makes by editing cells on a held
circle, or as searches on a changed model): the slope with no piles, each row
removed on the critical circle, the Ito & Matsui coefficients and forces, the
capacities removed, the spacing sweep, and the two-circle table. Not walked: the
Analysis Report (File → Generate Report… writes a Word file whose tables need a
document reader to check), the assistant, and the method table (held by the
page's tag).
"""
import copy
import re

from tutorial_walks.harness import (StepFailed, column, editable, field, headers,
                                    labels, preflight_items, press, rows, set_field,
                                    set_value, show_parameters_for, tick)

PAGE = "lem12_piles.md"
TIMEOUT = 3600


def _radius(c):
    return c["R"] if c.get("R") is not None else c["Yo"] - c["Depth"]


def _held(sd, crit):
    """``sd`` with the critical circle as its only circle."""
    out = copy.deepcopy(sd)
    out["circles"] = [dict(sd["circles"][0], Xo=crit["Xo"], Yo=crit["Yo"], Depth=crit["Depth"],
                           R=crit["Yo"] - crit["Depth"])]
    return out


def _piles(sd, **change):
    out = copy.deepcopy(sd)
    for pl in out["pile_lines"]:
        pl.update(change)
    return out


def _report(bundle):
    """The pile report of a solve, by row label."""
    recs = bundle["slice_df"].attrs.get("pile_report") or []
    return {r["label"]: r for r in recs}


def _set_cells(st, **cells):
    """Piles, Table view: every row's ``cells`` typed."""
    def act(dlg):
        press(dlg, "Table view")
        t = editable(dlg)
        for h, v in cells.items():
            c = column(t, h)
            vals = v if isinstance(v, (list, tuple)) else [v] * t.table.rowCount()
            for r, x in enumerate(vals):
                it = t.table.item(r, c)
                if it is None:
                    raise StepFailed("the %s cell of row %d cannot be typed into" % (h, r + 1))
                it.setText(str(x))
    st.edit("Piles", act)


def _spencer(st, analysis="Auto search", grid=None, min_slip=None, seen=None):
    def act(dlg):
        set_field(dlg, "Method", "Spencer")
        set_field(dlg, "Analysis", analysis)
        if grid is not None:
            tick(dlg, "Grid search (auto-seed the circular search)", grid)
        if min_slip is not None:
            w_ = field(dlg, "Min slip depth")
            if seen is not None:
                seen["min_slip_enabled_before"] = w_.isEnabled()
            tick(dlg, "Ignore surficial (skin) failures", True)
            if seen is not None:
                seen["min_slip_enabled_after"] = w_.isEnabled()
            set_value(w_, min_slip)
        if seen is not None:
            seen["items"] = preflight_items(dlg) or []
            seen["log_at_open"] = st.log_text()
    return st.run_lem(act)


def walk(w):
    p = w.page
    st = w.studio()
    tbl = p.table(r"The columns after the endpoints control how the pile force")

    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_piles.xlsx")))

    with w.step("The pile rows"):
        seen = {}

        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "Table view")
            t = editable(dlg)
            seen["hdr"] = [h.split(" (")[0] for h in headers(t)[0]]
            seen["rows"] = rows(t)
            press(dlg, "List view")
            from PySide6.QtWidgets import QGroupBox
            seen["groups"] = [g.title() for g in dlg.findChildren(QGroupBox) if g.isVisibleTo(dlg)]
            seen["labels"] = [l.text() for l in labels(dlg)]
            return "reject"
        st.edit("Piles", act)
        w.check(seen["hdr"] == tbl[0], "Table view on LEM alone heads %s, not the page's %s"
                % (seen["hdr"], tbl[0]))
        for want, got in zip(tbl[1:], seen["rows"]):
            same = all(a == b or (_f(a) is not None and _f(a) == _f(b)) for a, b in zip(want, got))
            w.check(same, "a pile row reads %s, not the page's %s" % (got, want))
        groups = p.find(r"as a form in four groups — \*\*(\w+)\*\*, \*\*(\w+)\*\*, \*\*([\w /]+)\*\*, \*\*(\w+)\*\*").groups()
        for g in groups:
            w.check(g in seen["groups"], "List view has no %r group (it has %s)" % (g, seen["groups"]))
        for lab in re.findall(r"\*\*(Vcap \(per element, lb\))\*\* and \*\*(Mcap \(per element, lb·ft\))\*\*",
                              p.flat)[0]:
            w.check(any(l.rstrip(":") == lab for l in seen["labels"]),
                    "List view has no label %r" % lab)

    with w.step("Running the analysis"):
        _s, log = _spencer(st)
        b = st.results["lem_solution"]
        w.check_log(p.code_block(r"the force each row applies:"), log, "Spencer search")
        crit = b["search"]["fs_cache"][0]
        w.expect("FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered at")
        w.expect("center x", crit["Xo"], r"centered at \((\d+\.\d+), \d+\.\d+\) with a radius of \d+\.\d+, tangent")
        w.expect("center y", crit["Yo"], r"centered at \(\d+\.\d+, (\d+\.\d+)\) with a radius of \d+\.\d+, tangent")
        w.expect("radius", _radius(crit), r"with a radius of (\d+\.\d+), tangent")
        w.expect("tangent elevation", crit["Depth"], r"tangent at elevation (−\d+\.\d+)\. The solution plot")
        w.expect("height above the base", crit["Depth"] - float(st.data["max_depth"]),
                 r"whose lowest point is (\d+\.\d+) ft above the base")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.expect("entry in front of the toe", -min(xs), r"The surface enters (\d+\.\d) ft out in front of the toe")
        w.expect("exit behind the crest", max(xs) - 20.0, r"exits (\d+\.\d) ft behind the crest and carries")
        rep = _report(b)
        for lab, k in (("lower row", 0), ("upper row", 1)):
            r = rep.get(lab)
            w.require(r is not None, "the run reports no crossing of the %s" % lab)
            w.expect("%s crossing elevation" % lab, r["y"],
                     r"at elevation (−\d+\.\d+) on the lower shaft and (−\d+\.\d+) on the upper", group=k + 1)
            w.expect("%s depth to the surface" % lab, r["depth"],
                     r"— (\d+\.\d+) ft and (\d+\.\d+) ft of soil above the surface", group=k + 1)
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Single surface")))
        w.expect("the file's circle alone", st.results["lem_solution"]["results"]["FS"],
                 r"runs along the rigid base and gives (\d\.\d+);")

    sd = w.load(w.copy_model(w.page_link("xslope_piles.xlsx"), "engine"))
    base = w.lem(sd, "spencer")
    crit = base["search"]["fs_cache"][0]
    held = _held(sd, crit)

    with w.step("Contribution of each pile row"):
        b0 = w.lem(dict(copy.deepcopy(sd), pile_lines=[]), "spencer")
        c0 = b0["search"]["fs_cache"][0]
        w.expect("no-pile FS", b0["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* on a shallower circle from the toe, tangent")
        w.expect("no-pile tangent", c0["Depth"], r"from the toe, tangent at elevation (−\d\.\d+), moving")
        w.expect("no-pile mass", float(b0["slice_df"]["w"].sum()), r"moving ([\d,]+) lb/ft against the")
        tab = p.table(r"we remove each row on its own:")
        keep = {"Neither": [], "Lower row only": ["lower row"], "Upper row only": ["upper row"],
                "Both": ["lower row", "upper row"]}
        for r in tab[1:]:
            m = copy.deepcopy(held)
            m["pile_lines"] = [pl for pl in m["pile_lines"] if pl["label"] in keep[r[0]]]
            bb = w.lem(m, "spencer", analysis="single_surface")
            w.expect("%s on the critical circle" % r[0], bb["results"]["FS"],
                     r"\| %s \| (\d\.\d+) \|" % re.escape(r[0]))

    with w.step("How the force is computed"):
        from xslope.ito_matsui import ito_matsui_coefficients
        mat = sd["materials"][0]
        pl = sd["pile_lines"][0]
        a1, a2 = ito_matsui_coefficients(float(pl["S"]) - float(pl["D_pile"]), float(pl["D_pile"]),
                                         float(mat["phi"]))
        w.expect("A1", a1, r"they are \*\*A₁ = (\d+\.\d+) ft\*\*")
        w.expect("A2", a2, r"and \*\*A₂ = (\d+\.\d+) ft\*\*\.")
        rep = _report(base)
        tbl2 = p.table(r"On the critical surface above:")
        for i, lab in enumerate(("lower row", "upper row")):
            r = rep[lab]
            row = lambda name: re.escape(name) + r"[^|]*\| ([\d,.]+) \| ([\d,.]+) \|"
            w.expect("%s z_f" % lab, r["depth"], row("Depth to the surface, z<sub>f</sub> (ft)"), group=i + 1)
            w.expect("%s F_pile" % lab, r["F_soil"], row("Soil force per pile, F<sub>pile</sub> (lb)"), group=i + 1)
            w.expect("%s F_pile/S" % lab, r["F_soil"] / float(pl["S"]),
                     row("Soil force per foot of slope, F<sub>pile</sub>/S (lb/ft)"), group=i + 1)

    with w.step("Effect of the structural capacity"):
        rep = _report(base)
        for i, lab in enumerate(("lower row", "upper row")):
            r = rep[lab]
            row = lambda name: re.escape(name) + r"[^|]*\| ([\d,.]+) \| ([\d,.]+) \|"
            w.expect("%s L_m" % lab, r["L_m"], row("Moment arm, L<sub>m</sub> (ft)"), group=i + 1)
            w.expect("%s force used per pile" % lab, r["F_used"], row("Force used, per pile (lb)"), group=i + 1)
            w.expect("%s force per foot" % lab, r["H_width"], row("Force used, per foot of slope (lb/ft)"), group=i + 1)
        caps = p.table(r"with only the two capacity cells changing")
        caps_at = p.find(r"with only the two capacity cells changing").end()
        for r in caps[1:]:
            ch = {"Vcap and Mcap (as shipped)": {}, "Mcap only": {"V_cap": None},
                  "Vcap only": {"M_cap": None}, "Neither": {"V_cap": None, "M_cap": None}}[r[0]]
            bb = w.lem(_piles(held, **ch), "spencer", analysis="single_surface")
            tot = sum(x["H_width"] for x in _report(bb).values())
            at = r"\| %s \| " % re.escape(r[0])
            w.expect("ΣH, %s" % r[0], tot, at + r"([\d,.]+) \|", nth=0, start=caps_at)
            w.expect("FS, %s" % r[0], bb["results"]["FS"], at + r"[\d,.]+ \| (\d+\.\d+) \|",
                     nth=0, start=caps_at)

    with w.step("Effect of pile spacing: the fixed-surface sweep"):
        from xslope.ito_matsui import ito_matsui_coefficients
        sw = p.table(r"ΣH is what survives the capacities and reaches the slices:")
        for r in sw[1:]:
            S = float(r[0])
            bb = w.lem(_piles(held, S=S), "spencer", analysis="single_surface")
            rep = _report(bb)
            a1, a2 = ito_matsui_coefficients(S - 2.0, 2.0, float(sd["materials"][0]["phi"]))
            at = r"\| %s \| [\d.]+ \| " % re.escape(r[0])
            w.expect("A1 at S = %s" % r[0], a1, at + r"(\d+\.\d+) \|")
            w.expect("A2 at S = %s" % r[0], a2, at + r"\d+\.\d+ \| (\d+\.\d+) \|")
            w.expect("lower-row soil force at S = %s" % r[0], rep["lower row"]["F_soil"],
                     at + r"\d+\.\d+ \| \d+\.\d+ \| ([\d,]+) / [\d,]+ \|")
            w.expect("upper-row soil force at S = %s" % r[0], rep["upper row"]["F_soil"],
                     at + r"\d+\.\d+ \| \d+\.\d+ \| [\d,]+ / ([\d,]+) \|")
            w.expect("ΣH at S = %s" % r[0], sum(x["H_width"] for x in rep.values()),
                     at + r"\d+\.\d+ \| \d+\.\d+ \| [\d,]+ / [\d,]+ \| ([\d,.]+) \|")
            w.expect("FS at S = %s" % r[0], bb["results"]["FS"],
                     at + r"\d+\.\d+ \| \d+\.\d+ \| [\d,]+ / [\d,]+ \| [\d,.]+ \| (\d\.\d+) \|")

    with w.step("Effect of pile spacing: 12 ft in Studio"):
        _set_cells(st, S=p.number(r"\*\*S\*\* from `6` to `(\d+)`")[0])
        _spencer(st)
        b = st.results["lem_solution"]
        w.expect("FS at 12 ft", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, and the surface is not the deep one")
        xs = [x for x, _y in b["failure_surface"].coords]
        ys = [y for _x, y in b["failure_surface"].coords]
        w.check(abs(min(xs)) < 0.05, "the 12 ft surface does not leave the toe (x = %.3f)" % min(xs))
        w.expect("exit behind the crest at 12 ft", max(xs) - 20.0, r"exits (\d+\.\d) ft behind the crest, and carries")
        w.expect("mass at 12 ft", float(b["slice_df"]["w"].sum()), r"and carries ([\d,]+) lb/ft against the deep")

    with w.step("Effect of pile spacing: 3 ft and the warnings"):
        # The engine prints its very-large-force warning once per process. The
        # walk's own engine runs above (the fixed-surface tables) are not a
        # reader's Studio session and must not use it up.
        import xslope.slice as _slice
        _slice._ito_matsui_warned = False
        _set_cells(st, S="3")
        seen = {}
        mark = st.mark()
        _spencer(st, seen=seen)
        w.check(any("1 warning" in t for t in seen["items"]) and any("S/D" in t for t in seen["items"]),
                "Model checks does not warn about S/D = 1.5: %s" % seen["items"][:4])
        w.check("very large" in seen["log_at_open"][mark:],
                "opening Run LEM… does not print the warning that the computed force is very large")
        _set_cells(st, S="6")

    with w.step("Entering the pile force directly"):
        hs = p.table(r"two values the report printed:")
        _set_cells(st, H=[r[0] for r in hs[1:]])
        _s, log = _spencer(st)
        b = st.results["lem_solution"]
        c = b["search"]["fs_cache"][0]
        w.expect("stated-force FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a shallower circle from the toe — tangent")
        w.expect("stated-force tangent", c["Depth"], r"from the toe — tangent at elevation (−\d\.\d+), [\d.]+ ft below")
        w.expect("stated-force mass", float(b["slice_df"]["w"].sum()), r"at its deepest, carrying ([\d,]+) lb/ft\. Solved")
        tot = re.findall(r"Total pile resistance = ([\d,.]+) lb/ft", log)
        w.expect("stated ΣH", float(tot[-1].replace(",", "")) if tot else None,
                 r"the stated numbers remain at the ([\d,.]+) that was right")
        shallow = _held(sd, c)
        stated = [float(r[0]) for r in hs[1:]]
        sd_stated = copy.deepcopy(sd)
        for pl, h in zip(sd_stated["pile_lines"], stated):
            pl["H"] = h
        two = p.table(r"Solved on that circle alone, the two ways of entering the pile no longer agree:")
        cases = [(held, "The auto search's critical circle"), (shallow, "The shallower circle above")]
        for circ, name in cases:
            at = r"\| %s \| " % re.escape(name)
            comp = w.lem(dict(copy.deepcopy(sd), circles=circ["circles"]), "spencer", analysis="single_surface")
            stat = w.lem(dict(copy.deepcopy(sd_stated), circles=circ["circles"]), "spencer", analysis="single_surface")
            none = w.lem(dict(copy.deepcopy(sd), circles=circ["circles"], pile_lines=[]), "spencer",
                         analysis="single_surface")
            w.expect("computed, %s" % name, comp["results"]["FS"], at + r"(\d\.\d+) \|")
            w.expect("stated, %s" % name, stat["results"]["FS"], at + r"\d\.\d+ \| (\d\.\d+) \|")
            w.expect("no piles, %s" % name, none["results"]["FS"], at + r"\d\.\d+ \| \d\.\d+ \| (\d\.\d+) \|")
            if circ is shallow:
                w.expect("computed ΣH on the shallow circle", sum(x["H_width"] for x in _report(comp).values()),
                         r"the computed path develops \*\*([\d,.]+) lb/ft\*\*")
        _set_cells(st, H="")
        w.require(all(pl.get("H") in (None, "") for pl in st.data["pile_lines"]), "H did not clear")

    with w.step("Looking above the pile row"):
        seen = {}
        _s, log = _spencer(st, grid=True)
        b = st.results["lem_solution"]
        w.expect("grid FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, below the 1\.845")
        coords = list(b["failure_surface"].coords)
        low = min(coords, key=lambda c: c[0])
        w.expect("grid entry x", low[0], r"at \((\d+\.\d+), \d+\.\d+\), the head of the upper shaft")
        w.expect("grid crest exit", max(c[0] for c in coords), r"runs back to the crest at x = (\d+\.\d)")
        w.expect("grid mass", float(b["slice_df"]["w"].sum()), r"moving ([\d,]+) lb/ft\. The depth")
        tot = re.findall(r"Total pile resistance = ([\d,.]+) lb/ft", log)
        w.expect("force on the bypass surface", float(tot[-1].replace(",", "")) if tot else 0.0,
                 r"the Log reports \*\*(\d\.\d) lb/ft\*\*")
        fs_grid = b["results"]["FS"]
        depth = p.number(r"Min slip depth\*\* box below it, and (\d+) ft there")[0]
        _spencer(st, grid=True, min_slip=depth, seen=seen)
        w.check(seen["min_slip_enabled_before"] is False and seen["min_slip_enabled_after"] is True,
                "Min slip depth is not grayed until Ignore surficial (skin) failures is ticked")
        w.check(abs(st.results["lem_solution"]["results"]["FS"] - fs_grid) < 5e-4,
                "a %s ft minimum slip depth changes the grid answer" % depth)
        b1 = w.lem(dict(copy.deepcopy(sd), pile_lines=[]), "spencer", grid_seed=True)
        w.expect("no piles, from the grid", b1["results"]["FS"],
                 r"returns 1\.149 from the file's circle and (1\.\d+) from the grid")


def _f(v):
    try:
        return float(str(v).replace("−", "-"))
    except (TypeError, ValueError):
        return None
