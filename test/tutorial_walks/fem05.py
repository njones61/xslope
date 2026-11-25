"""FEM-5, A Rock Slope on Its Joints: both parts in Studio from the starter.

Walked: Part 1 — the two joint rows pasted from the page's table, the mesh and
its Log line, the run's dialog defaults at the page's iteration budget, the
factor of safety, and the slip along the bedding plane at the last converged
state. Part 2 — the joint region entered in the polygons editor, the column set
from Build network… with the page's fields, the base plane row, the mesh and its
counts, the run; the page's claim about entering the coordinates at three
decimals, checked by meshing them; the Voronoi generator's line count in the
same dialog; and the completed Voronoi model meshed and run. Not walked: the
spacing table of "Solution time for a joint network" (three more meshes of the
same region, printed by the figure producer) and the 200,000-iteration rerun.
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    paste, press, set_field, set_value, tick, value)

PAGE = "fem05_rock_slope_joints.md"
TIMEOUT = 3 * 3600


def _mesh_counts(log):
    m = re.findall(r"Mesh built: (\d+) nodes, (\d+) elements(?:, \d+ 1D elements)?, (\d+) joint "
                   r"elements on (\d+) jointed line", log)
    if not m:
        raise StepFailed("the Log reports no mesh with joint elements")
    return tuple(int(x) for x in m[-1])


def _mesh(dlg):
    set_field(dlg, "Element type", "Quadratic triangles (tri6)")
    tick(dlg, "Auto-size from geometry", False)
    set_field(dlg, "Target element size", "1.5")


def _run(w, budget):
    seen = {}

    def act(dlg):
        seen["f"] = {k: value(field(dlg, k)) for k in
                     ("Analysis", "F min (SSRM)", "F max (SSRM)", "Tolerance (SSRM)", "Failure criterion")}
        from tutorial_walks.harness import button
        seen["accel"] = button(dlg, "Accelerate convergence").isChecked()
        set_field(dlg, "Max iterations per trial", budget)
    return act, seen


def _fs(st):
    r = (st.results.get("fem_solution") or {}).get("ssrm_result")
    if r is None:
        raise StepFailed("the run left no strength-reduction result")
    return r["FS"]


def walk(w):
    p = w.page
    st = w.studio()
    budget = p.number(r"Set \*\*Max iterations per trial\*\* to \*\*([\d,]+)\*\*")[0]

    with w.step("Opening the starter file"):
        st.open(w.copy_model(w.page_link("xslope_rock_joints_start.xlsx"), "part1"))
        w.check(str(st.data.get("unit_system")).lower() == "si", "the starter is not metric")
        w.check(len(st.data.get("polygons") or []) == 1 and not st.data.get("profile_lines"),
                "the section is not a single polygon")
        w.check(not st.data.get("joint_lines"), "the joints worksheet is not empty")

    with w.step("Entering the two joint lines"):
        say = r"the two rows below, or paste them straight into the (?:worksheet|table)\."
        rows = page_rows(p.table(r"[Ee]nter " + say))
        table_view = bool(p.finditer(r"press \*\*Table view\*\*, and enter " + say))

        def act(dlg):
            if table_view:
                press(dlg, "Table view")
            t = editable(dlg)
            paste(t, rows, 0, column(t, "Label"))
        st.edit("Joints", act)
        jl = st.data.get("joint_lines") or []
        w.require(len(jl) == 2, "%d joint rows saved" % len(jl))
        for j, r in zip(jl, rows):
            w.require(j["label"] == r[0] and [float(j[k]) for k in ("x1", "y1", "x2", "y2", "c", "phi")]
                      == [float(x) for x in r[1:7]], "a joint row saved as %s" % j)

    with w.step("Building the mesh"):
        st.set_mode("FEM")
        _s, log = st.build_mesh(_mesh)
        n, e, j, lines = _mesh_counts(log)
        w.expect("nodes", n, r"The Log pane reports what it built: \*\*(\d+) nodes, \d+ elements")
        w.expect("elements", e, r"built: \*\*\d+ nodes, (\d+) elements, and")
        w.expect("joint elements", j, r"elements, and (\d+) joint elements on 2 jointed lines")

    with w.step("Running the slab"):
        act, seen = _run(w, budget)
        st.run_fem(act)
        f = seen["f"]
        w.check(f["Analysis"].startswith("SSRM") and abs(float(f["F min (SSRM)"]) - 1) < 1e-9 and
                abs(float(f["F max (SSRM)"]) - 2) < 1e-9 and abs(float(f["Tolerance (SSRM)"]) - 0.01) < 1e-9,
                "the run does not open on SSRM, 1.0 to 2.0, 0.01: %s" % f)
        w.check(f["Failure criterion"].startswith("Hybrid"), "Failure criterion opens on %s, not Hybrid" % f["Failure criterion"])
        w.check(seen["accel"], "Accelerate convergence is not ticked")
        w.expect("slab FS", _fs(st), r">>\*\*FS = (1\.\d+)\*\* The results view")
        b = st.results["fem_solution"]
        from xslope import fem_details as FD
        sol, fsol = b["solution"], b.get("failure_solution")
        bed = next((r for r in FD.list_lines(b["fem_data"], sol, st.data, field_state="converged",
                                             failure_solution=fsol)
                    if r["kind"] == "joint" and "bedding" in str(r.get("label", ""))), None)
        if bed is None:
            w.fail("the bedding joint is not in the results")
        else:
            pr = FD.joint_profile(b["fem_data"], sol, bed["index"], st.data, field_state="converged",
                                  failure_solution=fsol)
            slip = np.abs(np.asarray(pr["slip"], float)) * 1000.0
            w.expect("largest bedding slip (mm)", float(slip.max()), r"grows from \d\.\d+ mm at the release joint to (\d\.\d+) mm")
            w.expect("smallest bedding slip (mm)", float(slip.min()), r"grows from (\d\.\d+) mm at the release joint")

    region = page_rows(p.table(r"up the face to the crest at \(14, 12\), and back along the ground surface"))
    base = page_rows(p.table(r"Enter the row below (?:in the joints editor|on the joints editor's \*\*Table view\*\*),? or paste it in\."))

    def part2(sub, region_rows, base_rows):
        st.open(w.copy_model(w.page_link("xslope_rock_joints_start.xlsx"), sub))

        def poly(dlg):
            press(dlg, "Add polygon")
            set_field(dlg, "Type", "joints")
            set_field(dlg, "Name", "Toppling zone")
            paste(editable(dlg), region_rows)
        st.edit("Polygons", poly)
        seen = {}

        def joints(dlg):
            press(dlg, "Build network")
            if p.finditer(r"Enter the row below on the joints editor's \*\*Table view\*\*"):
                press(dlg, "Table view")
            t = editable(dlg)
            press(dlg, "Add row")
            paste(t, base_rows, t.table.rowCount() - 1, column(t, "Label"))

        def network(dlg):
            set_field(dlg, "Name", "col")
            set_field(dlg, "Kind", "Parallel set")
            set_field(dlg, "Dip (deg)", p.number(r"\*\*(\d+)°\*\* here is a set standing")[0])
            set_field(dlg, "Spacing", p.number(r"between one trace and the next, (\d+\.\d+) m")[0])
            set_field(dlg, "Offset", p.number(r"(\d+\.\d+) m here moves the traces")[0])
            set_field(dlg, "Within", "Joint region: Toppling zone")
            set_field(dlg, "c", "0")
            set_field(dlg, "phi", "35")
            st.settle(40)
            seen["status"] = dlg._status.text()
        st.edit("Joints", joints, network)
        return seen

    with w.step("Part 2: the joint region, the column set and the base plane"):
        seen = part2("part2", region, base)
        n_lines = p.count(r"that set comes to \*\*(\w+)\*\* traces")
        w.check("%d joint lines" % n_lines in seen["status"],
                "Build network… reports %r, not %d joint lines" % (seen["status"], n_lines))
        jl = st.data.get("joint_lines") or []
        names = [j["label"] for j in jl]
        w.require(names[:n_lines] == ["col-%02d" % (i + 1) for i in range(n_lines)] and names[-1] == base[0][0],
                  "the joints worksheet holds %s" % names)
        w.check(st.tree_count("Joints") == str(n_lines + 1), "the Joints row reads %s" % st.tree_count("Joints"))

    with w.step("Running the toppling stack"):
        st.set_mode("FEM")
        _s, log = st.build_mesh(_mesh)
        counts = _mesh_counts(log)
        n, e, j, lines = counts
        w.expect("toppling nodes", n, r"The Log pane reports \*\*(\d+) nodes, \d+ elements, and \d+ joint elements on 6")
        w.expect("toppling elements", e, r"reports \*\*\d+ nodes, (\d+) elements, and \d+ joint elements on 6")
        w.expect("toppling joint elements", j, r"elements, and (\d+) joint elements on 6 jointed lines")
        act, _seen = _run(w, budget)
        st.run_fem(act)
        w.expect("toppling FS", _fs(st), r">>\*\*FS = (1\.\d+)\*\* !\[The deformed blocks: the columns lean")

    with w.step("Part 2 at three decimals"):
        def three(rows_):
            return [[("%.3f" % float(c)) if re.fullmatch(r"-?\d+\.\d+", c) else c for c in r] for r in rows_]
        part2("part2_3dp", three(region), three(base))
        st.set_mode("FEM")
        fails = bool(p.finditer(r"rounded to three decimals, the mesh build fails"))
        mark = st.mark()
        try:
            _s, log = st.build_mesh(_mesh)
            built = st.data.get("mesh") is not None and "Mesh built" in log
        except StepFailed:
            built, log = False, st.log_since(mark)
        if fails:
            w.check(not built, "the page says the mesh build fails on coordinates rounded to "
                               "three decimals; it builds (%s)" % (str(_mesh_counts(log)) if built else ""))
        else:
            w.check(built, "the three-decimal entry does not mesh")
        if built:
            w.check(_mesh_counts(log) == counts, "the three-decimal entry meshes to %s, not the "
                                                 "four-decimal %s" % (str(_mesh_counts(log)), str(counts)))
            if not fails:
                # Nothing on the page asks for four decimals: the rounded entry has to
                # give the page's answer too.
                act, _seen = _run(w, budget)
                st.run_fem(act)
                w.expect("toppling FS at three decimals", _fs(st),
                         r">>\*\*FS = (1\.\d+)\*\* !\[The deformed blocks: the columns lean")

    with w.step("Solution time for a joint network: the Voronoi model"):
        seen = {}

        def joints(dlg):
            press(dlg, "Build network")
            return "reject"

        def network(dlg):
            set_field(dlg, "Kind", "Voronoi (blocky mass)")
            set_field(dlg, "Block size", p.number(r"at a block size of (\d+) m")[0])
            set_field(dlg, "Seed", p.number(r"and a seed of (\d+) —")[0])
            set_field(dlg, "Within", "Joint region: Toppling zone")
            st.settle(60)
            seen["status"] = dlg._status.text()
            return "reject"
        st.edit("Joints", joints, network)
        n_vor = p.number(r"it writes \*\*(\d+) joint lines\*\*")[0]
        w.check("%s joint lines" % n_vor in seen["status"],
                "the Voronoi generator reports %r, not %s joint lines" % (seen["status"], n_vor))
        st.open(w.copy_model(w.page_link("xslope_rock_voronoi.xlsx"), "voronoi"))
        st.set_mode("FEM")
        _s, log = st.build_mesh(_mesh)
        n, e, j, lines = _mesh_counts(log)
        w.expect("Voronoi nodes", n, r"it comes to \*\*([\d,]+) nodes, \d+ elements and \d+ joint")
        w.expect("Voronoi elements", e, r"nodes, (\d+) elements and \d+ joint elements\*\*")
        w.expect("Voronoi joint elements", j, r"elements and (\d+) joint elements\*\*")
        act, _seen = _run(w, budget)
        st.run_fem(act)
        w.expect("Voronoi FS", _fs(st), r">>\*\*FS = (1\.\d+)\*\* The mass comes")
