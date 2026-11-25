"""FEM-2, Reinforcement: LEM vs FEM, in Studio from the starter.

Walked: Spencer on the starter and its five crossings; the mesh at the page's
settings; the soils' elastic properties typed into the columns the page names;
the layers' E and Area; the elastic-perfectly-plastic run and 1D Details at both
field states; the residual entered and its run; the overburden pullout law and
its run; Spencer again with the law. Not walked: the stiffness and residual
sweeps of the last two sections (the figure producer's runs, held by the page's
tags where it has them).
"""
import re

from tutorial_walks.harness import (StepFailed, column, editable, field, list_rows,
                                    preflight_items, press, run_button_enabled,
                                    set_field, set_value, show_parameters_for, tick,
                                    type_cells, value)

PAGE = "fem02_reinforcement.md"
TIMEOUT = 3 * 3600


def _fem_fs(st):
    r = (st.results.get("fem_solution") or {}).get("ssrm_result")
    if r is None:
        raise StepFailed("the run left no strength-reduction result")
    return r["FS"]


def _reinf(st, table_view=False, **cols):
    def act(dlg):
        if table_view:
            press(dlg, "Table view")
        show_parameters_for(dlg, {"FEM"})
        t = editable(dlg)
        for h, v in cols.items():
            c = column(t, h)
            for r in range(t.table.rowCount()):
                it = t.table.item(r, c)
                if it is None:
                    raise StepFailed("the %s cell of row %d cannot be typed into" % (h, r + 1))
                it.setText(v)
    st.edit("Reinforcement lines", act)


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the starter file"):
        st.open(w.copy_model(w.page_link("xslope_reinforced_slope_start.xlsx")))
        w.check(st.mode == "lem", "the starter opens in %s mode" % st.mode)

    with w.step("The limit equilibrium answer"):
        seen = {}

        def act(dlg):
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
            seen["items"] = preflight_items(dlg) or []
        st.run_lem(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]),
                "Model checks: %s" % seen["items"][:3])
        w.check(any("stiffness" in t.lower() for t in seen["items"]),
                "Model checks has no note about the lines' missing axial stiffness: %s" % seen["items"])
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("Spencer FS", b["results"]["FS"], r"\*\*Spencer's method gives FS = (1\.\d+)\*\*")
        w.expect("center x", crit["Xo"], r"circle centered at \((−\d+\.\d+), \d+\.\d+\)")
        w.expect("center y", crit["Yo"], r"circle centered at \(−\d+\.\d+, (\d+\.\d+)\)")
        w.expect("radius", crit["Yo"] - crit["Depth"], r"with a radius of (\d+\.\d+) ft\. The surface")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.expect("crest exit", max(xs), r"comes out on the crest at \*\*x = (\d+\.\d)\*\*")
        from xslope.crossings import crossings
        found = crossings(b["slice_df"], st.data, b["results"]["FS"]) or {"rows": []}
        crossed = found["rows"]
        w.check(len(crossed) == 5, "the circle crosses %d layers, not five" % len(crossed))
        w.check(all(abs(float(r["applied"]) - 800.0) < 1e-6 for r in crossed),
                "not every crossing takes 800 lb/ft: %s" % [r["applied"] for r in crossed])
        w.expect("total tension", sum(float(r["applied"]) for r in crossed),
                 r"for ([\d,]+) lb/ft of tension in total")

    with w.step("Building the mesh"):
        st.set_mode("FEM")
        size = p.number(r"Untick \*\*Auto-size from geometry\*\* and set \*\*Target element size\*\* to `(\d+)`")[0]

        def act(dlg):
            set_field(dlg, "Element type", "Quadratic triangles (tri6)")
            tick(dlg, "Auto-size from geometry", False)
            set_field(dlg, "Target element size", size)
            tick(dlg, "Refine thin zones", False)
        st.build_mesh(act)
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(m["elements"]), r"nodes and ([\d,]+) triangles\*\*, with the six")
        w.expect("bar elements", (0 if m.get("elements_1d") is None else len(m["elements_1d"])), r"discretized into \*\*(\d+) bar elements\*\*")

    with w.step("The stiffnesses the run needs: the soils"):
        named = p.find(r"in the\s?columns headed `([^`]+)` and `([^`]+)`;").groups()
        tbl = p.table(r"in the\s?columns headed `[^`]+` and `[^`]+`;")
        w.check(list(named) == tbl[0], "the sentence names %s, the table heads %s" % (named, tbl[0]))

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"FEM"})
            t = editable(dlg)
            for r in range(t.table.rowCount()):
                for h, v in zip(tbl[0], tbl[1]):
                    type_cells(t, [[v]], r, column(t, h))
        st.edit("Materials", act)
        for m in st.data["materials"]:
            w.require((float(m["E"]), float(m["nu"])) == (float(tbl[1][0]), float(tbl[1][1])),
                      "%s saved E / ν as %s / %s" % (m["name"], m["E"], m["nu"]))

    with w.step("The stiffnesses the run needs: the layers"):
        tbl = p.table(r"Fill `E \(psf\)` and `Area` on all six rows")
        _reinf(st, table_view=bool(p.finditer(r"Open \*\*Reinforcement\*\* in the Inputs dock, press \*\*Table view\*\*")),
               **dict(zip(tbl[0], tbl[1])))
        for rl in st.data["reinforcement_lines"]:
            w.require(float(rl.get("E") or 0) == float(tbl[1][0]) and float(rl.get("area") or rl.get("Area") or 0)
                      == float(tbl[1][1]), "a layer saved E / Area as %s / %s"
                      % (rl.get("E"), rl.get("area", rl.get("Area"))))
        w.require(st.data.get("mesh") is not None, "entering the stiffnesses dropped the mesh")

    with w.step("The elastic-perfectly-plastic run"):
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
            seen["run"] = run_button_enabled(dlg)
        st.run_fem(act)
        w.check(seen["run"] and any("1 warning" in t for t in seen["items"]),
                "Model checks does not read 1 warning with Run enabled: %s" % seen["items"][:2])
        w.expect("elastic-perfectly-plastic FS", _fem_fs(st), r"\*\*FS = (1\.\d+)\*\*, \*\*1\.4% below\*\*")

    with w.step("Reading what each layer is doing"):
        btn = st.mw.fem_details_btn
        w.require(btn is not None and btn.isEnabled(), "there is no 1D Details… button")
        btn.click()
        st.settle()
        dlg = st.mw.fem_details_dlg
        w.check(dlg.field_state.currentText() == "At failure", "1D Details does not open on At failure")
        lines = [t for t in list_rows(dlg.list) if t.strip().lower().startswith("line")]
        w.check(len(lines) == 6 and all("100%" in t for t in lines),
                "at failure the six layers do not all read 100%%: %s" % lines)
        set_value(dlg.field_state, "Last converged")
        st.settle(5)
        lines = [t for t in list_rows(dlg.list) if t.strip().lower().startswith("line")]
        states = {t.split()[0] + " " + t.split()[1]: t for t in lines}
        w.check("pullout" in states.get("line 1", "").lower() or "pullout" in states.get("Line 1", "").lower(),
                "line 1 does not read pullout at the last converged state: %s" % lines)
        w.check(sum("yielded" in t.lower() for t in lines) == 5,
                "lines 2 to 6 do not all read yielded: %s" % lines)
        dlg.close()

    with w.step("The peak-residual run"):
        tres = p.number(r"enter `(\d+)` in the `Tres` column on all six")[0]
        _reinf(st, Tres=tres)
        st.run_fem(lambda dlg: None)
        w.expect("peak-residual FS", _fem_fs(st), r"\*\*FS = (1\.\d+)\*\*\. Post-peak behavior cost")

    with w.step("Pullout from the overburden"):
        ad = p.number(r"enter `(\d+)` under `Adhesion \(psf\)`")[0]
        de = p.number(r"and `(\d+)` under `Delta` on all six")[0]
        _reinf(st, **{"Tres": "", "Adhesion (psf)": ad, "Delta": de})
        st.run_fem(lambda dlg: None)
        w.expect("overburden-law FS", _fem_fs(st), r"\*\*FS = (1\.\d+)\*\*, against 1\.566 with the stated")
        st.set_mode("LEM")
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Auto search")))
        b = st.results["lem_solution"]
        w.expect("Spencer with the law", b["results"]["FS"],
                 r"the limit equilibrium one by \d\.\d%, to \*\*(1\.\d+)\*\* when the search")

