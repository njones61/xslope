"""COMBO-2, Rapid Drawdown: the three routes, built up in Studio from the starter.

Walked: Part 1 (the materials, the two piezometric lines, the drawdown run and
its Log), Part 2 (the lines cleared, both boundary sets from the page's tables,
the mesh, the steady solves of both sets, u set to seep, the drawdown run),
and Part 3 (Set 2 cleared and what that does to the mesh, the mesh rebuilt, the
schedule and stage times typed into the Transient editor, the reservoir boundary
bound to the series, the transient run's frames, the two stage frames, and the
drawdown run on them). Not walked: the drained comparison runs at the end of
Part 2, the sweep of the core's d, and the verification problems "Which source
to use" quotes (VP96 to VP99, whose numbers the verification tags and the
tutorial restatement check hold).
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, button, column, editable, field,
                                    list_rows, page_rows, paste, press, rows,
                                    seep_bc_set, select_row, set_field, set_value,
                                    show_parameters_for, tick, type_cells, value)

PAGE = "combo02_rapid_drawdown.md"
TIMEOUT = 2 * 3600


def _stages(log):
    m = re.findall(r"Stage 1 FS = ([\d.]+), Stage 2 FS = ([\d.]+), Stage 3 FS = ([\w. ]+?)\s*$",
                   log, re.M)
    if not m:
        raise StepFailed("the Log reports no rapid drawdown stages")
    return m[-1]


def _rapid(dlg):
    set_field(dlg, "Method", "Spencer")
    set_field(dlg, "Analysis", "Auto search")
    tick(dlg, "Rapid drawdown", True)


def _check_stage_table(w, log, after, label):
    """The page's Stage | | FS table against the Log's stage line."""
    tbl = w.page.table(after)
    got = _stages(log)
    for r, g in zip(tbl[1:], got):
        if r[2].startswith("not"):
            w.check(g.startswith("not"), "%s: stage %s reads %s, not 'not required'" % (label, r[0], g))
        else:
            w.check(g == r[2], "%s: stage %s reads %s; the page prints %s" % (label, r[0], g, r[2]))


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("The dam: open the starter"):
        st.open(w.copy_model(w.page_link("xslope_johnson_rapid_start.xlsx")))
        w.check(st.mode == "lem", "the starter opens in %s mode" % st.mode)

    with w.step("The materials, and the one undrained zone"):
        tbl = p.table(r"narrows the table to the columns the drawdown run reads")

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM"})
            return "reject"
        st.edit("Materials", act)
        by = {m["name"]: m for m in st.data["materials"]}
        for r in tbl[1:]:
            m = by.get(r[1], {})
            got = [float(m.get(k) or 0) for k in ("gamma", "gamma_sat", "c", "phi", "d", "psi")]
            w.check(got == [float(x) for x in r[2:8]] and m.get("u") == r[8],
                    "%s reads %s / u %s, not the page's %s" % (r[1], got, m.get("u"), r[2:]))

    with w.step("The two lines"):
        l1 = page_rows(p.table(r"\*\*Line 1 — the full-pool surface\.\*\*"))
        l2 = page_rows(p.table(r"\*\*Line 2 — the surface at the end of the drawdown\.\*\*"))
        got1 = [tuple(map(float, c)) for c in st.data["piezo_line"]]
        got2 = [tuple(map(float, c)) for c in st.data["piezo_line2"]]
        w.check(got1 == [tuple(map(float, r)) for r in l1], "Line 1 is %s" % got1)
        w.check(got2 == [tuple(map(float, r)) for r in l2], "Line 2 is %s" % got2)
        st.edit("Piezometric lines", lambda dlg: None)

    with w.step("Running the drawdown on the pair"):
        _s, log = st.run_lem(_rapid)
        w.check_log(p.code_block(r"Click \*\*Run\*\*. The Log pane follows the search:"), log, "Part 1")
        _check_stage_table(w, log, r"the last lines of the log give the three stages", "Part 1")
        df = st.results["lem_solution"]["slice_df"]
        n_core = p.count(r"\*\*(\w+) cut the core\*\*")
        w.check(int((df["mat"] == 2).sum()) == n_core,
                "%d slices cut the core, not %d" % (int((df["mat"] == 2).sum()), n_core))
        w.expect("peak stage-2 pore pressure", float(df["u"].max()), r"reaching (\d+) psf on the slice table's `u` column, and the green")

    with w.step("Clearing the piezometric lines"):
        def act(dlg):
            from PySide6.QtWidgets import QTabWidget
            tabs = dlg.findChildren(QTabWidget)[0]
            for i in range(tabs.count()):
                tabs.setCurrentIndex(i)
                t = editable(dlg)
                t.table.selectAll()
                press(dlg, "Remove selected")
        st.edit("Piezometric lines", act)
        w.check(st.tree_count("Piezometric lines") == "0",
                "the Piezometric lines row reads %s, not 0" % st.tree_count("Piezometric lines"))

    def bc_set(s, heads, exit_pts):
        for _ in heads:
            press(s, "Add head")
        for i, (val, pts) in enumerate(heads, 1):
            select_row(s.list, "Head %d" % i)
            set_field(s, "Head value (ft)", val)
            paste(editable(s), pts)
        select_row(s.list, "Exit face")
        paste(editable(s), exit_pts)

    with w.step("The second boundary set"):
        at2 = p.find(r"\*\*Set 2 \(rapid drawdown\)\*\* is a second").start()
        h1 = (p.number(r"\*Head 1 — the reservoir\*, at \*\*Head value \(ft\)\*\* `(\d+)`")[0],
              page_rows(p.table(r"\*Head 1 — the reservoir\*")))
        h2 = (p.number(r"\*Head 2 — the tailwater\*, at `(\d+)`", nth=0)[0],
              page_rows(p.table(r"\*Head 2 — the tailwater\*, at")))
        ex = page_rows(p.table(r"\*Exit face\*, over the whole downstream slope, from"))
        r1 = (p.number(r"\*Head 1 — the residual pool\*, at \*\*Head value \(ft\)\*\* `(\d+)`")[0],
              page_rows(p.table(r"\*Head 1 — the residual pool\*")))
        r2 = (h2[0], page_rows(p.table(r"\*Head 2 — the tailwater\*, at `\d+`, unchanged")))
        ex2 = page_rows(p.table(r"\*Exit face\*, over the whole downstream slope, unchanged"))

        def act(dlg):
            bc_set(seep_bc_set(dlg, 0), [h1, h2], ex)
            bc_set(seep_bc_set(dlg, 1), [r1, r2], ex2)
        st.edit("Seep BC", act)
        sb, sb2 = st.data["seepage_bc"], st.data.get("seepage_bc2") or {}
        w.require(len(sb.get("specified_heads") or []) == 2 and len(sb2.get("specified_heads") or []) == 2,
                  "the boundary sets saved with %d and %d heads"
                  % (len(sb.get("specified_heads") or []), len(sb2.get("specified_heads") or [])))

    def mesh_auto(dlg):
        set_field(dlg, "Element type", "Linear triangles (tri3)")
        tick(dlg, "Auto-size from geometry", True)
        w.check(int(float(value(field(dlg, "Size divisions")))) == 100,
                "Auto-size is not at 100 size divisions")

    with w.step("Meshing, and both solves"):
        st.set_mode("Seepage")
        st.build_mesh(mesh_auto)
        mesh = st.data.get("mesh")
        w.require(mesh is not None, "no mesh was built")
        w.expect("mesh nodes", len(mesh["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(mesh["elements"]), r"nodes and ([\d,]+) triangles\*\*:")
        _s, log = st.run_seep(lambda dlg: None)
        its = [int(x) for x in re.findall(r"Converged in (\d+) iterations", log)]
        sols = st.results["seep_solutions"]
        w.require(1 in sols and 2 in sols, "one steady run did not solve both sets (%s)" % list(sols))
        w.expect("Set 1 iterations", its[0] if its else None, r"\*\*Set 1\*\* settles in \*\*(\d+) iterations\*\*")
        w.expect("Set 1 discharge", sols[1]["solution"]["flowrate"], r"total discharge of \*\*(\d+\.\d+) ft³/day per ft\*\*")
        w.expect("Set 2 iterations", its[1] if len(its) > 1 else None, r"\*\*Set 2\*\* settles in \*\*(\d+) iterations\*\*")
        w.expect("Set 2 discharge", sols[2]["solution"]["flowrate"], r"iterations\*\* to \*\*(\d+\.\d+) ft³/day per ft\*\*")
        for k, pat in ((1, r"The head runs from (\d+\.\d+) to 160\.000 ft"),
                       (2, r"The head runs from (\d+\.\d+) to 110\.000 ft")):
            w.expect("Set %d lowest head" % k, float(np.min(sols[k]["solution"]["head"])), pat)

    with w.step("Setting the materials' pore-pressure option"):
        def act(dlg):
            from PySide6.QtWidgets import QComboBox
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            c = column(t, "u")
            for r in range(t.table.rowCount()):
                cw = t.table.cellWidget(r, c)
                cb = cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox)
                set_value(cb, "seep")
        st.set_mode("LEM")
        st.edit("Materials", act)
        w.require(all(m["u"] == "seep" for m in st.data["materials"]), "u is not seep on every row")
        _s, log = st.run_lem(_rapid)
        _check_stage_table(w, log, r"settles a little upstream of where the pair's did", "Part 2")
        b = st.results["lem_solution"]
        w.expect("Part 2 FS", b["results"]["FS"], r"The rapid-drawdown factor of safety is (1\.\d+)\*\*, against")
        crit = b["search"]["fs_cache"][0]
        w.expect("Part 2 center x", crit["Xo"], r"at \((\d+\.\d+), \d+\.\d+\) with a radius of \d+\.\d+ ft:")
        w.expect("Part 2 center y", crit["Yo"], r"at \(\d+\.\d+, (\d+\.\d+)\) with a radius of \d+\.\d+ ft:")
        w.expect("Part 2 radius", crit.get("R") or crit["Yo"] - crit["Depth"],
                 r"with a radius of (\d+\.\d+) ft:")
        w.expect("Part 2 peak u", float(b["slice_df"]["u"].max()), r"peaks at (\d+) psf on its `u`")

    with w.step("Clearing boundary set 2"):
        def act(dlg):
            s2 = seep_bc_set(dlg, 1)
            s2.list.setCurrentRow(0)
            n = p.count(r"click \*\*Remove head\*\* (\w+), which takes out both heads")
            for _ in range(n):
                press(s2, "Remove head")
            select_row(s2.list, "Exit face")
            t = editable(s2)
            t.table.selectAll()
            press(s2, "Remove selected")
        st.edit("Seep BC", act)
        sb, sb2 = st.data["seepage_bc"], st.data.get("seepage_bc2") or {}
        w.require(not (sb2.get("specified_heads") or []) and not (sb2.get("exit_face") or []),
                  "Set 2 is not empty: %s" % sb2)
        w.require(len(sb.get("specified_heads") or []) == 2, "Set 1 was touched")
        says_cleared = bool(p.finditer(r"With that boundary gone the mesh no longer matches the model, "
                                       r"so Studio clears it"))
        if says_cleared:
            w.check(st.data.get("mesh") is None, "the page says clearing Set 2 clears the mesh; "
                                                 "Studio keeps it")
            pt = p.find(r"Set 2's residual pool ended at \((\d+), (\d+)\)")
            want = tuple(float(x) for x in pt.groups())
            w.check(want == tuple(map(float, r1[1][-1])),
                    "the page's (%g, %g) is not where Set 2's head ended" % want)
        if st.data.get("mesh") is None:
            st.set_mode("Seepage")
            st.build_mesh(mesh_auto)
        mesh = st.data["mesh"]
        w.expect("rebuilt mesh nodes", len(mesh["nodes"]), r"It comes out at the same ([\d,]+) nodes")
        w.expect("rebuilt mesh triangles", len(mesh["elements"]), r"the same [\d,]+ nodes and ([\d,]+) triangles")

    with w.step("The pool schedule and the stage times"):
        sched = p.table(r"\*\*The breakpoints\.\*\*")
        extra = p.find(r"\*\*Extra save times\*\* ((?:`\d+`,? )+)").group(1)
        extra = re.findall(r"`(\d+)`", extra)
        dur = p.number(r"\*\*Duration \(day\)\*\* `(\d+)`")[0]
        save = p.number(r"\*\*Save interval \(day\)\*\* `(\d+)`")[0]
        s1 = p.number(r"\*\*Stage 1 time \(day\)\*\* `(\d+)`")[0]
        s2 = p.number(r"\*\*Stage 2 time \(day\)\*\* `(\d+)`")[0]

        def act(dlg):
            from PySide6.QtWidgets import QTableWidgetItem
            set_value(dlg._name_edits[0], sched[0][1])
            tw = dlg._series_table
            for r, row in enumerate(sched[1:]):
                for c, v in enumerate(row):
                    it = tw.item(r, c)
                    if it is None:
                        tw.setItem(r, c, QTableWidgetItem(v))
                    else:
                        it.setText(v)
            sv = dlg._save_table
            for r, v in enumerate(extra):
                it = sv.item(r, 0)
                if it is None:
                    sv.setItem(r, 0, QTableWidgetItem(v))
                else:
                    it.setText(v)
            set_field(dlg, "Duration (day)", dur)
            set_field(dlg, "Save interval (day)", save)
            set_field(dlg, "Stage 1 time (day)", s1)
            set_field(dlg, "Stage 2 time (day)", s2)
        st.set_mode("Seepage")
        st.edit("Transient", act)
        ts = st.data.get("tseep") or {}
        w.require((float(ts.get("stage_1")), float(ts.get("stage_2"))) == (float(s1), float(s2)),
                  "the stage times saved as %s / %s" % (ts.get("stage_1"), ts.get("stage_2")))
        w.require(sorted(float(x) for x in (ts.get("save_times") or [])) ==
                  sorted(float(x) for x in extra),
                  "the extra save times saved as %s" % ts.get("save_times"))

    with w.step("Binding the reservoir boundary"):
        def act(dlg):
            s1w = seep_bc_set(dlg, 0)
            select_row(s1w.list, "Head 1")
            set_field(s1w, "Type", "reservoir")
            set_field(s1w, "Head value (ft)", "pool")
        st.edit("Seep BC", act)
        h = st.data["seepage_bc"]["specified_heads"][0]
        w.require(str(h.get("head")) == "pool", "Head 1 saved as %r" % h.get("head"))

    with w.step("Running the transient seepage analysis"):
        def act(dlg):
            set_field(dlg, "Run type", "Transient (time-dependent)")
            w.check(not field(dlg, "Convergence tol").isEnabled(), "Convergence tol is not grayed out")
            w.check(int(float(value(field(dlg, "Max iterations")))) == 2000,
                    "Max iterations does not move to 2000")
        _s, log = st.run_seep(act, timeout=3600)
        block = p.code_block(r"then prints a line per saved frame:")
        want_t = [float(m) for l in block for m in re.findall(r"t=(\d+): frame saved", l)]
        got_t = [float(m) for m in re.findall(r"t=(\d+(?:\.\d+)?): frame saved", log)]
        w.check(got_t == want_t, "the run saves frames at %s, not the page's %s" % (got_t, want_t))
        steps = re.findall(r"steps so far=(\d+)", log)
        w.expect_about("transient steps", int(steps[-1]) if steps else None,
                       r"The run took about \*\*(\d+) steps\*\*", 5)
        frames = (st.results.get("transient_seep") or {}).get("frames") or []
        w.check(len(frames) == p.count(r"\*\*(\w+) frames\*\* come out of the run"),
                "the run saves %d frames" % len(frames))
        f0 = next(f for f in frames if float(f["time"]) == 0.0)
        f50 = next(f for f in frames if float(f["time"]) == 50.0)
        w.expect("t = 0 inflow", f0["inflow"], r"inflow and outflow both read (\d+\.\d+) ft³/day per ft")
        w.expect("t = 0 outflow", f0["outflow"], r"inflow and outflow both read (\d+\.\d+) ft³/day per ft")
        w.check(abs(float(f50["inflow"])) < 5e-4, "inflow at t = 50 is %.4g, not zero" % float(f50["inflow"]))
        w.expect("t = 50 outflow", f50["outflow"], r"while (\d+\.\d+) ft³/day per ft still leaves")

    with w.step("The drawdown run on the transient frames"):
        seen = {}

        def act(dlg):
            _rapid(dlg)
            seen["s1"] = value(field(dlg, "Stage 1 time"))
            seen["s2"] = value(field(dlg, "Stage 2 time"))
        st.set_mode("LEM")
        _s, log = st.run_lem(act)
        w.check(float(seen["s1"]) == 0 and float(seen["s2"]) == 50,
                "the stage times read %s / %s" % (seen["s1"], seen["s2"]))
        _check_stage_table(w, log, r"pool the schedule has reached at t = 50", "Part 3")
        b = st.results["lem_solution"]
        w.expect("Part 3 FS", b["results"]["FS"], r"The rapid-drawdown factor of safety is (1\.\d+)\*\*, on a circle at")
        crit = b["search"]["fs_cache"][0]
        w.expect("Part 3 center x", crit["Xo"], r"on a circle at \((\d+\.\d+), \d+\.\d+\)")
        w.expect("Part 3 center y", crit["Yo"], r"on a circle at \(\d+\.\d+, (\d+\.\d+)\)")
        w.expect("Part 3 radius", crit.get("R") or crit["Yo"] - crit["Depth"],
                 r"\) with a radius of (\d+\.\d+) ft, its lowest point")
        df = b["slice_df"]
        w.expect("Part 3 peak u", float(df["u"].max()), r"reaches \*\*(\d+) psf\*\* on the slice table's")
        w.expect("Part 3 core slices", int((df["mat"] == 2).sum()), r"Three of the (\d+) core slices")


def _has(dlg, label):
    try:
        field(dlg, label)
        return True
    except StepFailed:
        return False
