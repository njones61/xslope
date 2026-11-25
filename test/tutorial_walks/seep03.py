"""SEEP-3, Transient Seepage: Reservoir Drawdown, built in Studio from the starter.

Walked, on a copy of the starter: the units and the material table read against
the page's row; the head and exit-face boundaries from the page's tables; the
mesh and its Log line; the steady run, its Log block, iteration count, Base
material, discharge, head and pore-pressure ranges and the phreatic surface over
the core; the schedule typed into the Transient editor; the reservoir boundary
bound to the series; the transient run (its dialog, notes, Log block, frames and
step counts); the frames' readings of the lag over the core; the head-history
table at its three nodes; and the inflow-outflow table. Not walked: the two
residual lines the page says differ between computers (the iteration count and
the flow are held instead), the four-panel figure and the assistant-drawn
history plot (made outside Studio's own views), and the mass-balance ledger's
cumulative volumes (not printed per frame).
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, editable, field, page_rows, paste,
                                    preflight_items, press, rows, seep_bc_set,
                                    select_row, set_field, set_value,
                                    show_parameters_for, tick, value)

PAGE = "seep03_reservoir_drawdown.md"
TIMEOUT = 3600


def _water_table(mesh, head, x, top=22.0):
    """The phreatic surface's elevation on the vertical at ``x``."""
    import matplotlib.tri as mtri
    nn = np.asarray(mesh["nodes"], float)
    el = np.asarray(mesh["elements"])[:, :3]
    f = mtri.LinearTriInterpolator(mtri.Triangulation(nn[:, 0], nn[:, 1], el),
                                   np.asarray(head, float) - nn[:, 1])
    ys = np.linspace(0.0, top, 22001)
    ps = np.asarray(f(np.full_like(ys, x), ys).filled(np.nan))
    ok = ~np.isnan(ps)
    ys, ps = ys[ok], ps[ok]
    cr = np.where((ps[:-1] > 0) & (ps[1:] <= 0))[0]
    if not len(cr):
        return None
    i = cr[-1]
    return float(ys[i] + (ys[i + 1] - ys[i]) * ps[i] / (ps[i] - ps[i + 1]))


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the starter file"):
        st.open(w.copy_model(w.page_link("xslope_earth_dam_drawdown_start.xlsx")))
        st.set_mode("Seepage")
        w.check(str(st.data.get("unit_system")).lower() == "si" and st.data.get("time_unit") == "day",
                "the starter is not SI / day")
        w.check(not (st.data.get("seepage_bc") or {}).get("specified_heads") and not st.data.get("tseep"),
                "the starter already has boundary conditions or a schedule")
        want = p.table(r"the rightmost of which need a scroll to the right to reach:")
        seen = {}

        def act(dlg):
            from tutorial_walks.harness import headers
            press(dlg, "Table view")
            show_parameters_for(dlg, {"Seepage"})
            t = editable(dlg)
            hdr = headers(t)[0]
            i = hdr.index("name")
            seen["hdr"] = hdr[i:i + len(want[0])]
            seen["rows"] = [r[i:i + len(want[0])] for r in rows(t)]
        st.edit("Materials", act)
        w.check(seen["hdr"] == want[0], "Seepage alone heads %s, not the page's %s" % (seen["hdr"], want[0]))
        for a, b in zip(want[1:], seen["rows"]):
            w.check(all(x == y or (_f(x) is not None and _f(x) == _f(y)) for x, y in zip(a, b)),
                    "a material row reads %s, not the page's %s" % (b, a))

    with w.step("Building the boundary conditions"):
        h1 = (p.number(r"\*\*Head 1 — the reservoir\.\*\* Leave \*\*Type:\*\* at `head` and set \*\*Head value \(m\):\*\* to `(\d+)`")[0],
              page_rows(p.table(r"\*\*Head 1 — the reservoir\.\*\*")))
        h2 = (p.number(r"\*\*Head 2 — the tailwater\.\*\* Press \*\*Add head\*\* again, set \*\*Head value \(m\):\*\* to `(\d+)`")[0],
              page_rows(p.table(r"\*\*Head 2 — the tailwater\.\*\*")))
        ex = page_rows(p.table(r"\*\*Exit face\.\*\* An exit face is the boundary"))

        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            press(s1, "Add head")
            press(s1, "Add head")
            for i, (val, pts) in enumerate((h1, h2), 1):
                select_row(s1.list, "Head %d" % i)
                set_field(s1, "Head value (m)", val)
                paste(editable(s1), pts)
            select_row(s1.list, "Exit face")
            paste(editable(s1), ex)
        st.edit("Seep BC", act)
        sb = st.data["seepage_bc"]
        w.require([float(h["head"]) for h in sb["specified_heads"]] == [float(h1[0]), float(h2[0])],
                  "the heads saved as %s" % sb["specified_heads"])

    with w.step("Building the mesh"):
        div = p.number(r"set \*\*Size divisions\*\* to `(\d+)`")[0]

        def act(dlg):
            set_field(dlg, "Element type", "Linear triangles (tri3)")
            tick(dlg, "Auto-size from geometry", True)
            set_field(dlg, "Size divisions", div)
        _s, log = st.build_mesh(act)
        line = p.find(r"the one the Log states when the mesh is built, `([^`]+)`").group(1)
        w.check(line in log, "the Log does not print %r" % line)
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(m["elements"]), r"nodes and ([\d,]+) triangles\*\*:")

    with w.step("The steady solution at full pool"):
        seen = {}

        def act(dlg):
            seen["run_type"] = _has(dlg, "Run type")
            w.check(abs(float(value(field(dlg, "Convergence tol"))) - 0.0001) < 1e-12, "Convergence tol is not 0.0001")
            w.check(int(float(value(field(dlg, "Max iterations")))) == 400, "Max iterations is not 400")
        _s, log = st.run_seep(act)
        w.check(not seen["run_type"], "the steady dialog has a Run type row before the model has a schedule")
        block = p.code_block(r"The Log pane opens with the path the solver took\. The head and the tail of it:")
        keep = [l for l in block if not l.startswith(("Iteration 22: Newton step, residual", "Converged in"))]
        w.check_log(keep, log, "steady run")
        its = re.findall(r"Converged in (\d+) iterations", log)
        w.expect("iterations", int(its[-1]) if its else None, r"Converged in (\d+) iterations \(residual")
        canvas = st.mw.seep_solution_canvas.get(1)
        panel = st.mw._display_panels.get(canvas)
        w.check(panel is not None and panel.base_mat.currentText() == "2: core", "Base material does not read 2: core")
        s = st.results["seep_solutions"][1]["solution"]
        w.expect("discharge", s["flowrate"], r"\*\*The total discharge is (\d\.\d+) m³/day per m\*\*")
        w.expect("lowest head", float(np.min(s["head"])), r"The head ranges from \*\*(\d+\.\d+) m to")
        w.expect("highest head", float(np.max(s["head"])), r"m to (\d+\.\d+) m\*\*, the two boundary")
        u = np.asarray(s["u"], float)
        w.expect("lowest pore pressure", u.min(), r"pore pressure runs from \*\*(−\d+\.\d) kPa to")
        w.expect("highest pore pressure", u.max(), r"kPa to (\d+\.\d) kPa\*\*;")
        mesh = st.data["mesh"]
        w.expect("water table at x = 50", _water_table(mesh, s["head"], 50.0), r"the surface stands at \*\*(\d+\.\d+) m at x = 50\*\*")
        w.expect("water table at x = 63", _water_table(mesh, s["head"], 63.0), r"and at \*\*(\d+\.\d+) m at x = 63\*\*")

    with w.step("Building the transient schedule"):
        w.check(st.tree_count("Transient") == "off", "the Transient row does not read off")
        sched = p.table(r"\*\*The breakpoints\.\*\* Fill the first three rows")
        dur = p.number(r"\*\*Duration \(day\)\.\*\* Set it to `(\d+)`")[0]
        save = p.number(r"\*\*Save interval \(day\)\.\*\* Set it to `(\d+)`")[0]
        extra = re.findall(r"`(\d+)`", p.find(r"Press \*\*Add\*\* under the list and enter (.+?)\. These are").group(1))

        def act(dlg):
            from PySide6.QtWidgets import QTableWidgetItem
            names = [e.text() for e in dlg._name_edits]
            w.check(names == ["t1", "t2", "t3", "t4", "t5"], "the series names open as %s" % names)
            set_value(dlg._name_edits[0], sched[0][1])
            tw = dlg._series_table
            for r, row in enumerate(sched[1:]):
                for c, v in enumerate(row):
                    it = tw.item(r, c)
                    if it is None:
                        tw.setItem(r, c, QTableWidgetItem(v))
                    else:
                        it.setText(v)
            set_field(dlg, "Duration (day)", dur)
            set_field(dlg, "Save interval (day)", save)
            sv = dlg._save_table
            for r, v in enumerate(extra):
                it = sv.item(r, 0)
                if it is None:
                    sv.setItem(r, 0, QTableWidgetItem(v))
                else:
                    it.setText(v)
        st.edit("Transient", act)
        ts = st.data.get("tseep") or {}
        w.require(float(ts.get("duration") or 0) == float(dur), "Duration saved as %s" % ts.get("duration"))
        w.require(sorted(float(x) for x in ts.get("save_times") or []) == sorted(float(x) for x in extra),
                  "the extra save times saved as %s" % ts.get("save_times"))
        w.check(ts.get("stage_1") in (None, "") and ts.get("stage_2") in (None, ""), "a stage time is set")
        w.check(st.tree_count("Transient") == "on", "the Transient row does not read on")
        w.require(st.data.get("mesh") is not None, "the schedule dropped the mesh")

    with w.step("Binding the reservoir boundary to the series"):
        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            select_row(s1.list, "Head 1")
            set_field(s1, "Type", "reservoir")
            set_field(s1, "Head value (m)", sched[0][1])
        st.edit("Seep BC", act)
        h = st.data["seepage_bc"]["specified_heads"][0]
        w.require(str(h.get("head")) == sched[0][1], "Head 1 saved as %r" % h.get("head"))
        w.require(st.data.get("mesh") is not None, "binding the series dropped the mesh")

    with w.step("Running the transient seepage analysis"):
        seen = {}

        def act(dlg):
            set_field(dlg, "Run type", "Transient (time-dependent)")
            w.check(not field(dlg, "Convergence tol").isEnabled(), "Convergence tol is not grayed out")
            w.check(int(float(value(field(dlg, "Max iterations")))) == 2000, "Max iterations does not jump to 2000")
            seen["items"] = preflight_items(dlg) or []
        _s, log = st.run_seep(act, timeout=1800)
        n_notes = sum(1 for t in seen["items"] if t.startswith("info:"))
        w.check(any("No problems found for this run" in t for t in seen["items"]) and n_notes == 2,
                "Model checks does not read No problems found with 2 notes: %s" % seen["items"][:4])
        w.check_log(p.code_block(r"and then prints a line per saved frame:"), log, "transient run")
        frames = (st.results.get("transient_seep") or {}).get("frames") or []
        times = [float(f["time"]) for f in frames]
        listed = [float(x) for x in re.findall(r"\d+", p.find(r"Overlaps are dropped, leaving ([\d, and]+)\.").group(1))]
        w.check(times == listed, "the frames are at %s, not %s" % (times, listed))
        steps = {float(t): int(n) for t, n in re.findall(r"t=(\d+): frame saved \(steps so far=(\d+)", log)}
        w.expect("steps", steps.get(360.0), r"\*\*The run took (\d+) steps to produce")
        w.expect("steps to day 47", steps.get(47.0), r"It reaches day 47 on step (\d+)")
        w.expect("share of steps by day 47 (%)", 100.0 * steps.get(47.0, 0) / max(steps.get(360.0, 1), 1),
                 r"take (\d+)% of its steps\*\*")
        w.expect("steps over the last 180 days", steps.get(360.0, 0) - steps.get(180.0, 0),
                 r"covers the last 180 days in \*\*(\d+) steps\*\*")

    def frame(t):
        f = min(frames, key=lambda fr: abs(float(fr["time"]) - t))
        if abs(float(f["time"]) - t) > 1e-6:
            raise StepFailed("no frame at t = %g" % t)
        return f

    with w.step("Reading the frames"):
        f30 = frame(30.0)
        w.check(abs(float(f30["inflow"])) < 5e-4, "day 30's inflow is not 0")
        w.expect("day 30 outflow", f30["outflow"], r"\*\*Inflow 0 / Outflow (\d\.\d+) m³/day per m\*\*")
        mesh = st.data["mesh"]
        lag = {}
        for t in (47.0, 80.0, 120.0, 180.0):
            wt = _water_table(mesh, frame(t)["head"], 54.5)
            lag[t] = None if wt is None else wt - 2.0
        w.expect("water table over the core, day 47", (lag[47.0] or 0) + 2.0,
                 r"the water table inside the core stands at \*\*elevation (\d+\.\d)\*\*")
        w.expect("lag, day 47", lag[47.0], r"a lag of \*\*(\d+\.\d) m, more than half")
        w.expect("lag, day 80", lag[80.0], r"from 8\.6 m on day 47 to (\d\.\d) m on day 80")
        w.expect("lag, day 120", lag[120.0], r"(\d\.\d) m on this frame")
        w.expect("lag, day 180", lag[180.0], r"and (\d\.\d) m on day 180")

    with w.step("The head history inside the core and the shell"):
        nn = np.asarray(st.data["mesh"]["nodes"], float)
        pts = {"core": (54.5, 8.9), "shell": (30.7, 8.9)}
        idx = {k: int(np.argmin(np.hypot(nn[:, 0] - x, nn[:, 1] - y))) for k, (x, y) in pts.items()}
        anchor = r"The table gives the values at four instants:"
        tbl = p.table(anchor)
        start = p.find(anchor).end()
        for r in tbl[1:]:
            f = frame(float(r[0]))
            at = r"\| %s \| [\d.]+ \| " % re.escape(r[0])
            hc = float(np.asarray(f["head"])[idx["core"]])
            hs = float(np.asarray(f["head"])[idx["shell"]])
            kw = dict(nth=0, start=start)
            w.expect("core h at t = %s" % r[0], hc, at + r"(\d+\.\d+) \|", **kw)
            w.expect("shell h at t = %s" % r[0], hs, at + r"\d+\.\d+ \| (\d+\.\d+) \|", **kw)
            w.expect("core ψ at t = %s" % r[0], hc - nn[idx["core"], 1],
                     at + r"\d+\.\d+ \| \d+\.\d+ \| ([+−]\d+\.\d+) \|", **kw)
            w.expect("shell ψ at t = %s" % r[0], hs - nn[idx["shell"], 1],
                     at + r"\d+\.\d+ \| \d+\.\d+ \| [+−]\d+\.\d+ \| ([+−]\d+\.\d+) \|", **kw)

    with w.step("Inflow, outflow, and the water released from storage"):
        tbl = p.table(r"Each saved frame reports the flow across the boundary in each direction:")
        start = p.find(r"Each saved frame reports the flow across the boundary in each direction:").end()
        for r in tbl[1:]:
            f = frame(float(r[0]))
            at = r"\| %s \| [\d.]+ \| " % re.escape(r[0])
            w.expect("inflow at t = %s" % r[0], f["inflow"], at + r"(\d\.\d+) \|", nth=0, start=start)
            w.expect("outflow at t = %s" % r[0], f["outflow"], at + r"\d\.\d+ \| (\d\.\d+) \|", nth=0, start=start)
        peak = max(frames, key=lambda fr: float(fr["outflow"]))
        w.expect("peak outflow", peak["outflow"], r"\*\*The outflow peaks at (\d\.\d+) on day 30\*\*")
        w.check(float(peak["time"]) == 30.0, "the outflow peaks on day %g, not 30" % float(peak["time"]))
        w.expect("last outflow as a share of the peak (%)", 100.0 * float(frames[-1]["outflow"]) / float(peak["outflow"]),
                 r"finishes at \*\*(\d\.\d+)% of the peak\*\*")


def _has(dlg, label):
    try:
        field(dlg, label)
        return True
    except Exception:
        return False


def _f(v):
    try:
        return float(str(v).replace("−", "-"))
    except (TypeError, ValueError):
        return None
