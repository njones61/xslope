"""FEM-4, Piles: LEM vs FEM, in Studio, then the spacing table and the benchmark.

Walked: the bare slope both ways; the pile model's Spencer search and the pile
forces its Log prints; the mesh, the Tip pinned in Piles, the strength
reduction and its 1D Details at failure; the tips fixed, the run, and 1D
Details at both field states; the pinned rows at 3 and 12 ft read at their
captured mechanisms for the spacing table (through the engine and
``fem_details.pile_profile``, the reading 1D Details makes); the VP106 Bishop
pair; and the wall entered in Piles from the starter, meshed, run and read in
1D Details. Not walked: the controls table of "What decides the finite element
answer", the tip-fixed spacing line, the limit equilibrium spacing curve (all
three are in the page's figure, from tools/make_tutorial_figures.py), and the
0.5 ft beam run of "A finer beam".
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, column, editable, field, list_rows,
                                    preflight_items, press, rows, run_button_enabled,
                                    set_field, set_value, show_parameters_for, tick,
                                    type_cells, value)

PAGE = "fem04_piles.md"
TIMEOUT = 3 * 3600


def _fem(st):
    b = st.results.get("fem_solution") or {}
    r = b.get("ssrm_result")
    if r is None:
        raise StepFailed("the run left no strength-reduction result")
    return b, r


def _bracket(dlg):
    set_field(dlg, "F min (SSRM)", "1.00")
    set_field(dlg, "F max (SSRM)", "2.00")


def _set_piles(st, **cells):
    def act(dlg):
        press(dlg, "Table view")
        t = editable(dlg)
        for h, v in cells.items():
            c = column(t, h)
            for r in range(t.table.rowCount()):
                cw = t.table.cellWidget(r, c)
                if cw is not None:
                    from PySide6.QtWidgets import QComboBox
                    cb = cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox)
                    set_value(cb, v)
                else:
                    t.table.item(r, c).setText(v)
    st.edit("Piles", act)


def _details(st):
    btn = st.mw.fem_details_btn
    if btn is None or not btn.isEnabled():
        raise StepFailed("there is no 1D Details… button on the results toolbar")
    btn.click()
    st.settle()
    dlg = st.mw.fem_details_dlg
    if dlg is None:
        raise StepFailed("1D Details did not open")
    return dlg


def _select(st, dlg, text, state=None):
    if state is not None:
        set_value(dlg.field_state, state)
    hits = [i for i in range(dlg.list.count()) if text in dlg.list.item(i).text()]
    if not hits:
        raise StepFailed("1D Details lists no %r (%s)" % (text, list_rows(dlg.list)))
    dlg.list.setCurrentRow(hits[0])
    st.settle(5)
    return dlg.current_profile()


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("The slope on its own"):
        st.open(w.copy_model(w.page_link("xslope_pile_wall_start.xlsx")))
        w.check(st.mode == "lem", "the file opens in %s mode" % st.mode)
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"),
                                set_field(dlg, "Analysis", "Auto search")))
        w.expect("bare-slope Spencer FS", st.results["lem_solution"]["results"]["FS"],
                 r"40 slices, \*\*Run\*\*: \*\*FS = (1\.\d+)\*\*")
        st.set_mode("FEM")
        st.build_mesh(lambda dlg: w.check(value(field(dlg, "Element type")) == "Quadratic triangles (tri6)",
                                          "Build Mesh does not open on tri6"))
        m = st.data["mesh"]
        w.expect("bare mesh nodes", len(m["nodes"]), r"mesh comes out at \*\*([\d,]+) nodes and [\d,]+ triangles\*\* with no")
        w.expect("bare mesh triangles", len(m["elements"]), r"\*\*[\d,]+ nodes and ([\d,]+) triangles\*\* with no")

        def act(dlg):
            w.check(abs(float(value(field(dlg, "F min (SSRM)"))) - 1.0) < 1e-9 and
                    abs(float(value(field(dlg, "F max (SSRM)"))) - 2.0) < 1e-9,
                    "the bracket does not open on 1.00 / 2.00")
        st.run_fem(act)
        _b, r = _fem(st)
        w.expect("bare-slope SSRM FS", r["FS"], r"and \*\*Run\*\*\. \*\*FS = (1\.\d+)\*\*\.")

    with w.step("The discrete pile row: the limit equilibrium answer"):
        st.doc.mark_clean() if hasattr(st.doc, "mark_clean") else None
        st.open(w.copy_model(w.page_link("xslope_piles.xlsx")))
        st.set_mode("LEM")
        _s, log = st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"),
                                          set_field(dlg, "Analysis", "Auto search")))
        b = st.results["lem_solution"]
        w.expect("piles Spencer FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a deep circle")
        tot = re.findall(r"Total pile resistance = ([\d,.]+) lb/ft", log)
        w.expect("force the rows deliver", float(tot[-1].replace(",", "")) if tot else None,
                 r"The rows deliver \*\*([\d,.]+) lb/ft\*\*")
        per = [float(x.replace(",", "")) for x in re.findall(r"-> ([\d,]+) lb per pile", log)]
        w.check(len(per) == 2, "the Log does not give the two rows' force per pile: %s" % per)
        if len(per) == 2:
            w.expect("lower row per shaft", per[0], r"— ([\d,]+) lb from\s+each shaft of the lower row")
            w.expect("upper row per shaft", per[1], r"lower row and ([\d,]+) lb from each of the upper")

    with w.step("The finite element answer"):
        st.set_mode("FEM")
        w.check(not st.menu_action("Run", "Run FEM").isEnabled(), "Run FEM… is enabled before a mesh")
        size = p.number(r"Untick \*\*Auto-size from geometry\*\* and set \*\*Target element size\*\* to `(\d+)`")[0]

        def mesh(dlg):
            set_field(dlg, "Element type", "Quadratic triangles (tri6)")
            tick(dlg, "Auto-size from geometry", False)
            set_field(dlg, "Target element size", size)
        st.build_mesh(mesh)
        m = st.data["mesh"]
        w.expect("pile mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and [\d,]+ triangles\*\*, with the two")
        w.expect("pile mesh triangles", len(m["elements"]), r"\*\*[\d,]+ nodes and ([\d,]+) triangles\*\*, with the two")
        w.expect("beam elements", _n1d(m), r"discretized into \*\*(\d+) beam elements\*\*", nth=0)
        _set_piles(st, Tip="pinned")
        w.require(all(pl.get("tip_fixity") == "pinned" for pl in st.data["pile_lines"]), "Tip is not pinned")
        w.require(st.data.get("mesh") is not None, "setting Tip dropped the mesh; the page says it is unaffected")
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
            seen["run"] = run_button_enabled(dlg)
            _bracket(dlg)
        st.run_fem(act)
        w.check(seen["run"] and any("No problems found" in t for t in seen["items"]),
                "Model checks does not find no problems: %s" % seen["items"][:3])
        w.check(any("Area" in t and "I" in t and ("solid" in t or "circular" in t) for t in seen["items"]),
                "Model checks shows no note about the derived solid section: %s" % seen["items"])
        _b, r = _fem(st)
        w.expect("pinned SSRM FS", r["FS"], r"\*\*FS = (1\.\d+)\*\*\. Spencer's method gave")

    with w.step("Head and tip fixity: 1D Details on the pinned rows"):
        dlg = _details(st)
        w.check(dlg.field_state.currentText() == "At failure", "1D Details does not open on At failure")
        prof = _select(st, dlg, "upper row")
        w.expect("upper row peak moment", abs(prof["max_moment"]), r"peaks 16 ft below the head at ([\d,]+) lb·ft per foot")
        w.expect("upper row peak moment depth", prof["max_moment_depth"], r"peaks (\d+) ft below the head at [\d,]+ lb·ft per foot")
        w.expect("upper row per shaft", abs(prof["max_moment"]) * 6.0, r"of slope, ([\d,]+) lb·ft per\s?shaft")
        w.expect("upper row share of capacity (%)", prof["peak_utilization"], r"per shaft, (\d+)% of the 60,000", scale=100.0)
        dlg.close()

    with w.step("Head and tip fixity: the tips fixed"):
        _set_piles(st, Tip="fixed")
        w.require(st.data.get("mesh") is not None, "setting Tip dropped the mesh")
        st.run_fem(_bracket)
        b, r = _fem(st)
        w.expect("tip-fixed SSRM FS", r["FS"], r"Run the same bracket again\. \*\*FS = (1\.\d+)\*\*")
        from xslope import fem_details
        # Read at the last converged state: this run's capture is a runaway a few
        # iterations long, where every element reads as yielded (the page's 1D
        # Details step below goes to Last converged for the same reason).
        fail = fem_details.field_solution(b["solution"], "converged", b.get("failure_solution")) or {}
        ym = np.asarray(fail.get("yielded_pile_M") if fail.get("yielded_pile_M") is not None else [], bool)
        yv = np.asarray(fail.get("yielded_pile_V") if fail.get("yielded_pile_V") is not None else [], bool)
        w.check(int(ym.sum()) == p.count(r"limits the gain: (\w+) of the 18 beam elements"),
                "%d beam elements stand at their moment capacity" % int(ym.sum()))
        w.check(int(yv.sum()) == 0, "%d beam elements stand at their shear capacity" % int(yv.sum()))
        dlg = _details(st)
        at_fail = _select(st, dlg, "upper row")
        w.check(at_fail["peak_utilization"] is not None and at_fail["peak_utilization"] <= 1.0 + 1e-6,
                "1D Details at failure reads %.0f%% of capacity on the tip-fixed upper row"
                % (100 * (at_fail["peak_utilization"] or 0)))
        prof = _select(st, dlg, "upper row", "Last converged")
        w.expect("tip-fixed head movement", abs(prof["u_lateral"][0]), r"The head moves (\d\.\d+) ft and the displacement")
        w.expect("tip-fixed moment", abs(prof["max_moment"]), r"reaches\s+the ([\d,]+) lb·ft per foot that is")
        w.expect("tip-fixed largest shear", abs(prof["max_shear"]), r"The largest shear is ([\d,]+) lb/ft at")
        w.expect("tip-fixed largest shear depth", prof["max_shear_depth"], r"The largest shear is [\d,]+ lb/ft at (\d+) ft")
        dlg.close()
        _set_piles(st, Tip="pinned")

    with w.step("What spacing does to each engine: the table"):
        from xslope import fem_details
        from xslope.fem import build_fem_data, solve_ssrm
        from xslope.mesh import (build_mesh_from_polygons, extract_constraint_line_geometry,
                                 extract_point_constraints, extract_size_regions,
                                 get_material_polygons)
        tbl = p.table(r"Each row below is read at that spacing's captured mechanism:")
        base = w.load(w.copy_model(w.page_link("xslope_piles.xlsx"), "engine"))
        for row in tbl[1:]:
            S = float(row[0])
            sd = dict(base, pile_lines=[dict(pl, S=S, tip_fixity="pinned") for pl in base["pile_lines"]])
            lines, _a, _b2 = extract_constraint_line_geometry(sd)
            mesh = w.quiet(build_mesh_from_polygons, get_material_polygons(sd, reinf_lines=lines),
                           float(size), "tri6", lines=lines or None,
                           point_constraints=extract_point_constraints(sd),
                           size_regions=extract_size_regions(sd))
            fd = build_fem_data(sd, mesh)
            res = w.quiet(solve_ssrm, fd, F_min=1.0, F_max=2.0, tolerance=0.01, debug_level=0,
                          capture_failure_state=True, failure_criterion="non_convergence",
                          max_iterations=12000)
            sol = dict(res["last_solution"])
            fail = res.get("failure_solution")
            peak = max(abs(fem_details.pile_profile(fd, sol, i, slope_data=sd, field_state="failure",
                                                    failure_solution=fail)["max_moment"] or 0)
                       for i in range(len(sd["pile_lines"])))
            at = r"\| %s \| [^|]+ \| " % re.escape(row[0])
            w.expect("S = %s peak moment per unit width" % row[0], peak, at + r"([\d,]+) \|")
            w.expect("S = %s peak moment per shaft" % row[0], peak * S, at + r"[\d,]+ \| ([\d,]+) \|")
            w.expect("S = %s fraction of Mcap (%%)" % row[0], 100 * peak * S / 60000.0,
                     at + r"[\d,]+ \| [\d,]+ \| (\d+)% \|")

    with w.step("A benchmark with a three-dimensional answer: the Bishop pair"):
        a = w.lem(w.load(w.copy_model(w.page_link("vp106a.xlsx"), "vp")), "bishop")
        c = w.lem(w.load(w.copy_model(w.page_link("vp106c.xlsx"), "vp")), "bishop")
        w.expect("Bishop with no pile", a["results"]["FS"], r"a Bishop search reads (1\.\d+) with no pile")
        w.expect("Bishop with the pile", c["results"]["FS"], r"with no pile and (1\.\d+) with it")
        w.expect("Bishop credit", c["results"]["FS"] / a["results"]["FS"], r"a credit of ×(1\.\d+)\.")

    with w.step("The continuous wall: adding the wall"):
        st.open(w.copy_model(w.page_link("xslope_pile_wall_start.xlsx"), "wall"))
        vals = {}
        m = p.find(r"Enter the geometry — `x1` (\S+), `y1` (\S+), `x2` (\S+), `y2` (\S+) — then `E` `(\d+)`")
        vals.update(dict(zip(("x1", "y1", "x2", "y2", "E"), m.groups())))
        vals["y2"] = vals["y2"].replace("−", "-")
        vals["I"] = p.find(r"`I` `([\d.]+)` and `Area`").group(1)
        vals["Area"] = p.find(r"`Area` `([\d.]+)` from the table above").group(1)
        vals["Mcap"] = p.find(r"`Mcap` `(\d+)`").group(1)
        vals["S"] = p.find(r"and `S` `(\d+)`\. Leave").group(1)
        seen = {}

        def act(dlg):
            press(dlg, "Table view")
            press(dlg, "Add row")
            t = editable(dlg)
            r = t.table.rowCount() - 1
            got = dict(zip([h for h in _hdr(t)], rows(t)[r]))
            w.check(got.get("Label") == "Pile", "the new row's Label reads %r, not 'Pile'" % got.get("Label"))
            w.check(got.get("Head") == "free" and got.get("Tip") == "free",
                    "the new row's Head / Tip read %s / %s" % (got.get("Head"), got.get("Tip")))
            t.table.item(r, column(t, "Label")).setText("sheet pile wall")
            for h, v in vals.items():
                it = t.table.item(r, column(t, h))
                if it is None:
                    raise StepFailed("the %s cell cannot be typed into" % h)
                it.setText(v)
            from PySide6.QtWidgets import QComboBox
            cw = t.table.cellWidget(r, column(t, "Tip"))
            set_value(cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox), "fixed")
        st.edit("Piles", act)

        def reopen(dlg):
            press(dlg, "Table view")
            t = editable(dlg)
            seen["E"] = rows(t)[0][_hdr(t).index("E")]
            return "reject"
        st.edit("Piles", reopen)
        shown = p.find(r"which the table shows as `([^`]+)` the next time it opens").group(1)
        w.check(seen["E"] == shown, "the reopened table shows E as %r, not %r" % (seen["E"], shown))
        done = w.load(w.copy_model(w.page_link("xslope_pile_wall.xlsx"), "done"))
        mine = st.data["pile_lines"][0]
        for k in ("x1", "y1", "x2", "y2", "E", "I", "Area", "M_cap", "S", "head_fixity", "tip_fixity"):
            w.check(mine.get(k) == done["pile_lines"][0].get(k) or
                    (_f(mine.get(k)) is not None and _f(mine.get(k)) == _f(done["pile_lines"][0].get(k))),
                    "the typed row's %s is %r; the completed file's is %r"
                    % (k, mine.get(k), done["pile_lines"][0].get(k)))

    with w.step("Re-meshing and running again"):
        st.set_mode("FEM")
        st.build_mesh(lambda dlg: None)
        m = st.data["mesh"]
        w.expect("wall mesh nodes", len(m["nodes"]), r"The mesh comes out at\s?\*\*([\d,]+) nodes and [\d,]+ triangles\*\*, one")
        w.expect("wall mesh triangles", len(m["elements"]), r"\*\*[\d,]+ nodes and ([\d,]+) triangles\*\*, one")
        w.expect("wall beam elements", _n1d(m), r"the wall discretized into \*\*(\d+) beam elements\*\*")
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
        st.run_fem(act)
        w.check(any("No problems found" in t for t in seen["items"]), "Model checks: %s" % seen["items"][:3])
        _b, r = _fem(st)
        w.expect("wall SSRM FS", r["FS"], r"\*\*FS = (1\.\d+)\*\*\. The wall takes the slope")

    with w.step("What the wall carries"):
        dlg = _details(st)
        w.check(dlg.field_state.currentText() == "At failure", "1D Details does not open on At failure")
        prof = _select(st, dlg, "sheet pile wall", "Last converged")
        w.expect("wall toe moment, converged", abs(prof["max_moment"]), r"toe — \*\*([\d,]+) lb·ft/ft\*\*, (?:\d+)% of")
        w.expect("wall share, converged (%)", prof["peak_utilization"], r"lb·ft/ft\*\*, (\d+)% of what the section", scale=100.0)
        w.expect("wall shear, converged", abs(prof["max_shear"]), r"peaking at ([\d,]+) lb/ft at \d+ ft\. The \*\*soil")
        w.expect("wall head movement", abs(prof["u_lateral"][0]), r"\*\*lateral displacement\*\* is (\d\.\d+) ft at the head")
        prof = _select(st, dlg, "sheet pile wall", "At failure")
        w.expect("wall toe moment, at failure", abs(prof["max_moment"]), r"raises the toe moment to ([\d,]+) lb·ft/ft")
        w.expect("wall share, at failure (%)", prof["peak_utilization"], r"lb·ft/ft, (\d+)% of the\s?capacity", scale=100.0)
        w.expect("wall shear, at failure", abs(prof["max_shear"]), r"from [\d,]+ to ([\d,]+) lb/ft, and its peak")
        dlg.close()


def _hdr(t):
    from tutorial_walks.harness import headers
    return headers(t)[0]


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _n1d(mesh):
    e = mesh.get("elements_1d")
    return 0 if e is None else len(e)
