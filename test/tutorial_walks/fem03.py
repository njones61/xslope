"""FEM-3, A Block Wall on Slip Joints: the three parts in Studio.

Walked: Part 1 — the seven contacts pasted from the page's table, the mesh at the
page's settings and its Log line, the run's dialog defaults at the page's
iteration budget and the Log's opening line, the factor of safety; Part 2 — the
three geogrid rows and their properties pasted into the FEM view of the
reinforcement editor, the mesh rebuilt and its counts, the run, the bound it
reports and the Log's paragraph, and the Continue button the bound brings;
Part 3 — the four sheet models opened, meshed as declared and run, and the
model checks on the two bonded sheets. Not walked: the optional
500,000-iteration run of "Results with a higher iteration limit" (about an
hour; the page reads it from a stored run) and the 1D Details readings (the
page's figures are drawn from tools/make_block_wall_figures.py).
"""
import re

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, preflight_items, press, set_field,
                                    show_parameters_for, tick, value)

PAGE = "fem03_block_wall_joints.md"
TIMEOUT = 4 * 3600


def _mesh_counts(log):
    m = re.findall(r"Mesh built: (\d+) nodes, (\d+) elements(?:, \d+ 1D elements)?, (\d+) joint "
                   r"elements on (\d+) jointed line", log)
    if not m:
        raise StepFailed("the Log reports no mesh with joint elements")
    return tuple(int(x) for x in m[-1])


def _wall_mesh(dlg):
    set_field(dlg, "Element type", "Quadratic triangles (tri6)")
    tick(dlg, "Auto-size from geometry", False)
    set_field(dlg, "Target element size", "0.8")
    tick(dlg, "Refine thin zones", False)


def _ssrm(budget, seen=None):
    def act(dlg):
        if seen is not None:
            seen["f"] = {k: value(field(dlg, k)) for k in
                         ("Analysis", "F min (SSRM)", "F max (SSRM)", "Tolerance (SSRM)",
                          "Failure criterion")}
            from tutorial_walks.harness import button
            seen["accel"] = button(dlg, "Accelerate convergence").isChecked()
        set_field(dlg, "Max iterations per trial", budget)
        if seen is not None:
            seen["items"] = preflight_items(dlg) or []
    return act


def _result(st):
    r = (st.results.get("fem_solution") or {}).get("ssrm_result")
    if r is None:
        raise StepFailed("the run left no strength-reduction result")
    return r


def walk(w):
    p = w.page
    st = w.studio()
    budget = p.number(r"\*\*Max iterations per trial — change it to ([\d,]+)\.\*\*")[0]

    with w.step("Entering the wall's contacts"):
        st.open(w.copy_model(w.page_link("xslope_block_wall_start.xlsx")))
        w.check(str(st.data.get("unit_system")).lower() == "si", "the starter is not metric")
        w.check(not st.data.get("joint_lines"), "the joints worksheet is not empty")
        say = r"enter the seven rows below, or paste them straight into the (?:worksheet|table)\."
        rows = page_rows(p.table(r"(?:E|press \*\*Table view\*\*, and e)" + say[1:]))
        table_view = bool(p.finditer(r"press \*\*Table view\*\*, and " + say))

        def act(dlg):
            if table_view:
                press(dlg, "Table view")
            t = editable(dlg)
            paste(t, rows, 0, column(t, "Label"))
        st.edit("Joints", act)
        jl = st.data.get("joint_lines") or []
        w.require([j["label"] for j in jl] == [r[0] for r in rows], "the joints saved as %s"
                  % [j["label"] for j in jl])
        done = w.load(w.copy_model(w.page_link("xslope_block_wall.xlsx"), "done"))
        for a, b in zip(jl, done["joint_lines"]):
            w.check(all(float(a[k]) == float(b[k]) for k in ("x1", "y1", "x2", "y2", "c", "phi")),
                    "joint %s differs from the completed file's" % a["label"])

    with w.step("Building the mesh"):
        st.set_mode("FEM")
        _s, log = st.build_mesh(_wall_mesh)
        n, e, j, lines = _mesh_counts(log)
        w.expect("nodes", n, r"The Log pane reports \*\*([\d,]+) nodes, \d+ elements, and \d+ joint elements on 7")
        w.expect("elements", e, r"\*\*[\d,]+ nodes, (\d+) elements, and \d+ joint elements on 7")
        w.expect("joint elements", j, r"elements, and (\d+) joint elements on 7")

    with w.step("Running the wall alone"):
        seen = {}
        _s, log = st.run_fem(_ssrm(budget, seen))
        f = seen["f"]
        w.check(f["Analysis"].startswith("SSRM") and abs(float(f["F min (SSRM)"]) - 1) < 1e-9 and
                abs(float(f["F max (SSRM)"]) - 2) < 1e-9 and abs(float(f["Tolerance (SSRM)"]) - 0.01) < 1e-9,
                "the run does not open on SSRM, 1.0 to 2.0, 0.01: %s" % f)
        w.check(f["Failure criterion"].startswith("Hybrid"), "Failure criterion opens on %s" % f["Failure criterion"])
        w.check(seen["accel"], "Accelerate convergence is not ticked")
        opening = p.find(r"The Log's opening lines for the run report whether it was on: `([^`]+)`").group(1)
        w.check(opening in log, "the Log does not open with %r" % opening)
        w.expect("wall FS", _result(st)["FS"], r">>\*\*FS = (1\.\d+)\*\*", nth=0)

    with w.step("Entering the three layers"):
        geo = page_rows(p.table(r"Then enter the three rows below, or paste them into the worksheet\."))
        props = p.table(r"Tres stays blank; the three layers are alike\.")

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"FEM"})
            t = editable(dlg)
            for _ in geo:
                press(dlg, "Add row")
            paste(t, geo, 0, column(t, "Label"))
            hdr = props[0]
            first = {"Adhesion": "Adhesion (kPa)", "E": "E (kPa)"}.get(hdr[0], hdr[0])
            paste(t, page_rows(props), 0, column(t, first))
        st.edit("Reinforcement lines", act)
        rl = st.data.get("reinforcement_lines") or []
        w.require([r["label"] for r in rl] == [g[0] for g in geo], "the layers saved as %s"
                  % [r["label"] for r in rl])
        done = w.load(w.copy_model(w.page_link("xslope_block_wall_grid.xlsx"), "done_grid"))
        keys = ("x1", "y1", "x2", "y2", "t_max", "tend1", "tend2", "adhesion", "delta", "E")
        for a, b in zip(rl, done["reinforcement_lines"]):
            diff = [k for k in keys if (_f(a.get(k)) or 0) != (_f(b.get(k)) or 0)]
            w.check(not diff and str(a.get("joint")).lower() == str(b.get("joint")).lower(),
                    "layer %s differs from the completed file's in %s" % (a["label"], diff or "Joint"))

    with w.step("What the geogrid adds"):
        if st.data.get("mesh") is not None:
            w.fail("adding the layers kept the old mesh; the page rebuilds it")
        _s, log = st.build_mesh(_wall_mesh)
        n, e, j, lines = _mesh_counts(log)
        w.expect("grid nodes", n, r"The mesh grows to \*\*([\d,]+) nodes, \d+ elements and \d+ joint")
        w.expect("grid elements", e, r"nodes, (\d+) elements and \d+ joint elements on 10")
        w.expect("grid joint elements", j, r"elements and (\d+) joint elements on 10 jointed")
        _s, log = st.run_fem(_ssrm(budget))
        r = _result(st)
        w.check(bool(r.get("fs_is_lower_bound")), "the run does not report a bound")
        # A lower bound is printed rounded down, as Studio's title and status bar print it.
        w.expect("grid FS bound", int(r["FS"] * 100) / 100.0, r">>\*\*FS ≥ (1\.\d+)\*\*")
        para = " ".join(p.find(r"and the Log gives the reason: > (.+?) The wall came to rest at every")
                        .group(1).replace(">", " ").split())
        w.check(para in " ".join(log.split()), "the Log does not give the page's paragraph: %r" % para[:80])
        btn = getattr(st.mw, "fem_continue_btn", None)
        w.check(btn is not None and not btn.isHidden() and btn.isEnabled(),
                "there is no Continue with a higher limit… button on the results toolbar")
        w.check("FS ≥ 1.50" in st.status(), "the status bar does not read FS ≥ 1.50: %r" % st.status())

    with w.step("Part 3: the four sheet models"):
        for name, pat in (("xslope_base_geotextile_bonded.xlsx", r"\(files/xslope_base_geotextile_bonded\.xlsx\):(?: <!--.*?-->)? >>\*\*FS = (\d\.\d+)\*\*"),
                          ("xslope_base_geotextile_jointed.xlsx", r"\(files/xslope_base_geotextile_jointed\.xlsx\):(?: <!--.*?-->)? >>\*\*FS = (\d\.\d+)\*\*"),
                          ("xslope_liner_bonded.xlsx", r"\(files/xslope_liner_bonded\.xlsx\):(?: <!--.*?-->)? >>\*\*FS = (\d\.\d+)\*\*"),
                          ("xslope_liner_jointed.xlsx", r"\(files/xslope_liner_jointed\.xlsx\):(?: <!--.*?-->)? >>\*\*FS = (\d\.\d+)\*\*")):
            st.open(w.copy_model(w.page_link(name), name[:-5]))
            st.set_mode("FEM")
            said = p.finditer(r"Mesh each file at \*\*tri6, (\d+\.\d+) m\*\*")

            def mesh(dlg, said=said):
                if said:
                    set_field(dlg, "Element type", "Quadratic triangles (tri6)")
                    tick(dlg, "Auto-size from geometry", False)
                    set_field(dlg, "Target element size", said[0].group(1))
            st.build_mesh(mesh)
            seen = {}
            st.run_fem(_ssrm(budget, seen))
            n_warn = sum(1 for t in seen["items"] if t.startswith("warning:"))
            if name == "xslope_liner_bonded.xlsx":
                w.check(n_warn == 3, "the bonded liner's Model checks show %d warnings, not three" % n_warn)
            elif name == "xslope_base_geotextile_bonded.xlsx":
                w.check(n_warn == 2, "the bonded base geotextile shows %d warnings, not two" % n_warn)
            else:
                w.check(n_warn == 0, "%s shows %d warnings" % (name, n_warn))
            w.expect(name, _result(st)["FS"], pat)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
