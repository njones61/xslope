"""SEEP-4, Infiltration and Flux Boundaries: the whole page in Studio, from the
starter file.

Walked: the starter's material row, the head and exit-face boundaries from the
page's tables, the mesh, the dry run, the three rain blocks (their list labels
as typed, and the mesh Studio keeps), the wet run, the doubled rain and the
other rates of the six-rate table, the storage properties, the storm typed into
the Transient editor, the blocks bound to the series (their labels), the
transient run, and the frames the page reads (days 0, 200, 300, 430, 600).
Not read: the per-boundary split of "Where the water comes from" (Studio does
not print it per boundary) and the storage balance of "Checking the run".
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, column, editable, field, list_rows,
                                    page_rows, paste, press, seep_bc_set, select_row,
                                    set_field, set_value, show_parameters_for,
                                    type_cells, value)

PAGE = "seep04_dam_infiltration.md"
TIMEOUT = 3600
START = "docs/tutorials/files/xslope_dam_infiltration_start.xlsx"
CREST_X = 26.0


def _iterations(log):
    m = re.findall(r"Converged in (\d+) iterations", log)
    return int(m[-1]) if m else None


def _exit_active(log):
    m = re.findall(r"(\d+)/(\d+) exit face active", log)
    return (int(m[-1][0]), int(m[-1][1])) if m else (None, None)


def _water_table(mesh, head, x=CREST_X, top=12.0):
    """The phreatic surface's elevation on the vertical at ``x`` and the pressure
    head at the ground surface there, interpolated on the mesh."""
    import matplotlib.tri as mtri
    nn = np.asarray(mesh["nodes"], float)
    el = np.asarray(mesh["elements"])[:, :3]
    f = mtri.LinearTriInterpolator(mtri.Triangulation(nn[:, 0], nn[:, 1], el),
                                   np.asarray(head, float) - nn[:, 1])
    ys = np.linspace(0.0, top, 12001)
    ps = np.asarray(f(np.full_like(ys, x), ys).filled(np.nan))
    ok = ~np.isnan(ps)
    ys, ps = ys[ok], ps[ok]
    cr = np.where((ps[:-1] > 0) & (ps[1:] <= 0))[0]
    wt = None
    if len(cr):
        i = cr[-1]
        wt = float(ys[i] + (ys[i + 1] - ys[i]) * ps[i] / (ps[i] - ps[i + 1]))
    return wt, float(ps[-1])


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the starter file"):
        st.open(w.copy_model(w.page_link("xslope_dam_infiltration_start.xlsx")))
        st.set_mode("Seepage")
        w.check(str(st.data.get("unit_system")).lower() == "si"
                and str(st.data.get("time_unit")).lower().startswith("sec"),
                "the starter is not SI / sec: %s / %s" % (st.data.get("unit_system"),
                                                         st.data.get("time_unit")))
        want = p.table(r"to \*\*Seepage\*\* alone. One row")
        got = {}

        def mats(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"Seepage"})
            t = editable(dlg)
            from tutorial_walks.harness import rows as _rows, headers
            names, _cols = headers(t)
            vals = _rows(t)[0]
            got.update(dict(zip(names, vals)))
            return "reject"
        st.edit("Materials", mats)
        for h, v in zip(want[0], want[1]):
            have = got.get(h)
            if have is None:
                w.fail("the Materials table has no column %r (it has %s)" % (h, sorted(got)))
                continue
            same = (have == v) or (_num(have) is not None and _num(have) == _num(v))
            w.check(same, "Materials: %s reads %r, the page says %r" % (h, have, v))

    head_pts = page_rows(p.table(r"\*\*Head 1 — the reservoir\.\*\*"))
    exit_pts = page_rows(p.table(r"\*\*Exit face — the toe drain\.\*\*"))
    head_val = p.number(r"set \*\*Head value \(m\):\*\* to `(\d+)`")[0]

    with w.step("Building the boundary conditions"):
        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            press(s1, "Add head")
            set_field(s1, "Head value (m)", head_val)
            paste(editable(s1), head_pts)
            select_row(s1.list, "Exit face")
            paste(editable(s1), exit_pts)
        st.edit("Seep BC", act)
        bc = st.data["seepage_bc"]
        w.require(len(bc["specified_heads"]) == 1 and
                  [tuple(map(float, c)) for c in bc["specified_heads"][0]["coords"]] ==
                  [tuple(map(float, r)) for r in head_pts],
                  "the head boundary saved as %s" % bc["specified_heads"])
        w.require([tuple(map(float, c)) for c in bc["exit_face"]] ==
                  [tuple(map(float, r)) for r in exit_pts],
                  "the exit face saved as %s" % bc["exit_face"])

    def mesh_dialog(dlg):
        set_field(dlg, "Element type", "Linear triangles (tri3)")
        from tutorial_walks.harness import tick
        tick(dlg, "Auto-size from geometry", False)
        w.check(abs(float(value(field(dlg, "Target element size"))) - 1.0) < 1e-9,
                "Target element size does not read 1.0")

    with w.step("Building the mesh"):
        st.build_mesh(mesh_dialog)
        mesh = st.data.get("mesh")
        w.require(mesh is not None, "no mesh was built")
        w.expect("mesh nodes", len(mesh["nodes"]), r"\*\*(\d+) nodes and \d+ triangles\*\*")
        w.expect("mesh triangles", len(mesh["elements"]), r"\*\*\d+ nodes and (\d+) triangles\*\*")

    def run_steady(dlg):
        from tutorial_walks.harness import preflight_items
        items = preflight_items(dlg) or []
        w.check(any("No problems found for this run" in t for t in items),
                "the Model checks panel does not read 'No problems found for this run.': %s" % items)

    with w.step("The dry-weather solution"):
        _s, log = st.run_seep(run_steady)
        sol = st.results["seep_solutions"][1]["solution"]
        w.expect("dry discharge", sol["flowrate"],
                 r"\*\*The discharge is (2\.\d+ × 10<sup>−7</sup>) m³/s per m\*\*")
        n_exit = p.count(r"the exit-face nodes — (\w+) of them")
        w.check("Number of exit face nodes: %d" % n_exit in log,
                "the solver does not report %d exit-face nodes" % n_exit)
        n_sweeps = p.count(r"this model closes in about (\w+) sweeps")
        w.check(_iterations(log) is not None and abs(_iterations(log) - n_sweeps) <= 1,
                "the dry run closes in %s sweeps, not about %d" % (_iterations(log), n_sweeps))
        a, n = _exit_active(log)
        w.check((a, n) == (13, 13), "the dry run ends with %s/%s exit-face nodes draining" % (a, n))
        w.expect("dry head max", float(np.max(sol["head"])), r"head ranges from \*\*0 to (\d+) m\*\*", nth=0)
        w.expect("dry head min", float(np.min(sol["head"])), r"head ranges from \*\*(\d+) to \d+ m\*\*", nth=0)
        wt, psi = _water_table(st.data["mesh"], sol["head"])
        w.expect("dry water table under the crest", wt, r"under the crest at about elevation (\d+\.\d)")
        w.expect("dry pressure head at the crest", psi, r"reaches (−\d+\.\d) m at the crest itself")
        dry_q = sol["flowrate"]

    fluxes = []
    for i, name in ((1, "the upstream face above the waterline"), (2, "the crest"),
                    (3, "the downstream face")):
        v = p.find(r"\*\*Flux %d — %s\.\*\* Flux value `([^`]+)`" % (i, re.escape(name))).group(1)
        pts = page_rows(p.table(r"\*\*Flux %d — %s\.\*\*" % (i, re.escape(name))))
        fluxes.append((v, pts))

    with w.step("Entering the flux boundary conditions"):
        seen = {}

        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            for v, pts in fluxes:
                press(s1, "Add flux")
                set_field(s1, "Flux value (m/sec)", v)
                paste(editable(s1), pts)
            seen["list"] = list_rows(s1.list)
        st.edit("Seep BC", act)
        lst = [re.sub(r"\s+", " ", t) for t in seen["list"]]
        if p.finditer(r"label in the list shows its rate the next time the editor opens"):
            # What the page says the list does while the editor is still open.
            w.check(not any("q = " in t and "blank" not in t for t in lst[1:4]),
                    "the page says a new block's label shows its rate only when the editor "
                    "reopens; the list shows the rates as typed: %s" % lst)
        quoted = p.find(r"The list now (?:reads|holds) `(Head 1 \(h = \d+\))`").group(1)
        w.check(lst and lst[0] == quoted, "the list's first row reads %r, not %r" % (lst[:1], quoted))
        for i, (v, _pts) in enumerate(fluxes, 1):
            row = next((t for t in lst if t.startswith("Flux %d" % i)), None)
            w.check(row is not None and _num(re.sub(r".*q = ([^)]+)\).*", r"\1", row)) is not None
                    and abs(_num(re.sub(r".*q = ([^)]+)\).*", r"\1", row)) - float(v)) <= 1e-6 * float(v),
                    "the list does not show Flux %d's rate as typed: %r" % (i, row))
        w.check(lst and lst[-1] == "Exit face", "the list does not end with 'Exit face': %s" % lst)
        got = [str(f["flux"]) for f in st.data["seepage_bc"]["specified_fluxes"]]
        w.require([_num(g) for g in got] == [float(v) for v, _ in fluxes],
                  "the fluxes saved as %s" % got)

    with w.step("Running it again"):
        mesh = st.data.get("mesh")
        if p.finditer(r"Adding boundary blocks clears the mesh"):
            # The page has the reader rebuild a mesh it says the blocks dropped.
            w.check(mesh is None, "the page says adding the boundary blocks clears the "
                                  "mesh; Studio keeps it")
            if mesh is None:
                st.build_mesh(mesh_dialog)
                mesh = st.data.get("mesh")
        w.require(mesh is not None, "adding the rain blocks dropped the mesh, and the page "
                                    "says there is no need to rebuild it")
        w.require(st.menu_action("Run", "Run Seep").isEnabled(), "Run → Run Seep… is grayed out")
        _s, log = st.run_seep(lambda dlg: None)
        sol = st.results["seep_solutions"][1]["solution"]
        w.expect("wet discharge", sol["flowrate"],
                 r"\*\*The discharge is (4\.\d+ × 10<sup>−7</sup>) m³/s per m\*\*")
        same = p.count(r"in the same (\w+) sweeps it took in dry weather")
        w.check(_iterations(log) == same,
                "the wet run takes %s sweeps, not the same %d" % (_iterations(log), same))
        a, n = _exit_active(log)
        w.check((a, n) == (13, 13), "the wet run ends with %s/%s exit-face nodes draining" % (a, n))
        w.expect("more water than dry (%)", 100.0 * (sol["flowrate"] / dry_q - 1.0),
                 r"the drain is passing \*\*(\d+)% more water\*\*")
        wet_wt, _psi = _water_table(st.data["mesh"], sol["head"])

    def rates(vals):
        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            for i, v in enumerate(vals, 1):
                select_row(s1.list, "Flux %d" % i)
                set_field(s1, "Flux value (m/sec)", v)
        return act

    with w.step("Doubling the rain rate"):
        dbl = p.table(r"double all three flux values")
        st.edit("Seep BC", rates([r[2] for r in dbl[1:]]))
        w.require(st.data.get("mesh") is not None, "changing the flux values dropped the mesh")
        st.run_seep(lambda dlg: None)
        sol = st.results["seep_solutions"][1]["solution"]
        w.expect("doubled discharge", sol["flowrate"],
                 r"The discharge comes out at \*\*(6\.\d+ × 10<sup>−7</sup>) m³/s per m\*\*")
        wt, _psi = _water_table(st.data["mesh"], sol["head"])
        w.expect("doubled water table", wt, r"stands at \*\*(\d+\.\d+) m\*\*, up from")
        w.expect("water table at the file's rate", wet_wt, r"m\*\*, up from (\d+\.\d+) m\.")

    with w.step("Six rates on one section"):
        tbl = p.table(r"draws the whole family on one section")
        at = p.find(r"draws the whole family on one section").end()
        cos = 2.0 / 5 ** 0.5
        for r in tbl[1:]:
            qk = r[0]
            if qk in ("0.1", "0.2"):
                continue
            q = float(qk) * 1e-7
            vals = ["%.10g" % (q * cos), "%.10g" % q, "%.10g" % (q * cos)] if q else ["0", "0", "0"]
            st.edit("Seep BC", rates(vals))
            _s, log = st.run_seep(lambda dlg: None)
            sol = st.results["seep_solutions"][1]["solution"]
            w.expect("discharge at q/k = %s" % qk, sol["flowrate"],
                     r"\| %s \| [^|]+ \| ([^|]+) \|" % re.escape(qk), nth=0, start=at)
            wt, _psi = _water_table(st.data["mesh"], sol["head"])
            if r[3].startswith("none"):
                w.check(wt is None, "at q/k = %s the dam has a water table at %.2f" % (qk, wt or 0))
                m = p.find(r"> (\d+) specified-flux node\(s\) finished with positive pore "
                           r"pressure \(max u = (\d+\.\d+)\)")
                w.check("%s specified-flux node(s) finished with positive pore pressure (max u = %s)"
                        % (m.group(1), m.group(2)) in log,
                        "the q/k = %s run does not warn as the page quotes" % qk)
            else:
                w.expect("water table at q/k = %s" % qk, wt,
                         r"\| %s \| [^|]+ \| [^|]+ \| (\d+\.\d+) \|" % re.escape(qk),
                         nth=0, start=at)
            a, n = _exit_active(log)
            w.check(a == n == 13, "at q/k = %s %s/%s drain nodes drain" % (qk, a, n))
        st.edit("Seep BC", rates([v for v, _ in fluxes]))
        w.require([float(f["flux"]) for f in st.data["seepage_bc"]["specified_fluxes"]] ==
                  [float(v) for v, _ in fluxes], "the three rates did not go back")

    with w.step("Storage properties"):
        tbl = p.table(r"fill the two columns at the right of the row")

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"Seepage"})
            t = editable(dlg)
            for h, v in zip(tbl[0][1:], tbl[1][1:]):
                type_cells(t, [[v]], 0, column(t, h))
        st.edit("Materials", act)
        m = st.data["materials"][0]
        w.require((float(m.get("Ss") or 0), float(m.get("Sy") or 0)) ==
                  (float(tbl[1][1]), float(tbl[1][2])),
                  "Ss / Sy saved as %s / %s" % (m.get("Ss"), m.get("Sy")))

    with w.step("The storm"):
        sched = p.table(r"type the schedule into the series table cell by cell")
        dur = p.number(r"Set \*\*Duration \(sec\)\*\* to `(\d+)`")[0]
        save = p.number(r"\*\*Save interval \(sec\)\*\* to `(\d+)`")[0]

        def act(dlg):
            from PySide6.QtWidgets import QTableWidgetItem
            names = dlg._name_edits
            w.require(names[0].text() == "t1" and names[1].text() == "t2",
                      "the Series names boxes do not read t1 / t2: %s / %s"
                      % (names[0].text(), names[1].text()))
            set_value(names[0], sched[0][1])
            set_value(names[1], sched[0][2])
            tw = dlg._series_table
            for r, row in enumerate(sched[1:]):
                if r >= tw.rowCount():
                    raise StepFailed("the series table has %d rows" % tw.rowCount())
                for c, v in enumerate(row):
                    it = tw.item(r, c)
                    if it is None:
                        tw.setItem(r, c, QTableWidgetItem(v))
                    else:
                        it.setText(v)
            set_field(dlg, "Duration (sec)", dur)
            set_field(dlg, "Save interval (sec)", save)
        st.edit("Transient", act)
        ts = st.data.get("tseep") or {}
        w.require(float(ts.get("duration") or 0) == float(dur), "Duration saved as %s" % ts.get("duration"))
        w.require(sorted((ts.get("series") or {}).keys()) == sorted(sched[0][1:]),
                  "the series saved as %s" % sorted((ts.get("series") or {}).keys()))

    with w.step("Binding the rain blocks to the schedule"):
        bind = p.table(r"replace each of the three flux values with the name")
        seen = {}

        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            for i, r in enumerate(bind[1:], 1):
                select_row(s1.list, "Flux %d" % i)
                set_field(s1, "Flux value (m/sec)", r[2])
            seen["list"] = [re.sub(r"\s+", " ", t) for t in list_rows(s1.list)]
        st.edit("Seep BC", act)
        quoted = [m.group(1) for m in p.finditer(r"`(Flux \d \(q = \w+\))`")]
        for q in quoted:
            w.check(q in seen["list"], "the list does not read %r once the names are typed "
                                       "(it reads %s)" % (q, seen["list"]))
        w.check(seen["list"][0] == "Head 1 (h = 10)", "the head boundary is not still 10")
        w.require(st.data.get("mesh") is not None, "binding the series dropped the mesh")

    with w.step("Running the transient analysis"):
        def act(dlg):
            set_field(dlg, "Run type", "Transient (time-dependent)")
        _s, log = st.run_seep(act, timeout=1800)
        tb = st.results.get("transient_seep") or {}
        frames = tb.get("frames") or []
        w.expect("saved frames", len(frames), r"The run reports \*\*(\d+) saved frames\*\*")
        m = re.findall(r"steps so far=(\d+)", log)
        w.expect_about("transient steps", int(m[-1]) if m else None,
                       r"reached in about \*\*(\d+) steps\*\*", 10)

    with w.step("Reading the frames"):
        day = 86400.0
        nn = np.asarray(st.data["mesh"]["nodes"], float)

        def frame(t):
            f = min(frames, key=lambda fr: abs(float(fr["time"]) - t))
            w.require(abs(float(f["time"]) - t) < 1.0, "no frame at t = %g" % t)
            return f
        f200 = frame(200 * day)
        w.expect("day 200 inflow", f200["inflow"], r"\*\*Inflow (4\.\d+ × 10<sup>−7</sup>) / Outflow")
        w.expect("day 200 outflow", f200["outflow"], r"Inflow 4\.\d+ × 10<sup>−7</sup> / Outflow (4\.\d+ × 10<sup>−7</sup>) m³/sec")
        _wt, psi200 = _water_table(st.data["mesh"], f200["head"])
        w.expect("day 200 pressure head at the crest", psi200, r"has eased to \*\*(−\d+\.\d) m\*\*")
        ic = int(np.argmin(np.hypot(nn[:, 0] - CREST_X, nn[:, 1] - 12.0)))
        if f200.get("u") is not None:
            w.expect("day 200 pore pressure at the crest", float(np.asarray(f200["u"])[ic]),
                     r"−\d+ kPa dry against (−\d+) kPa here")
        f300 = frame(300 * day)
        w.expect("day 300 inflow", f300["inflow"], r"reads \*\*Inflow (2\.\d+ × 10<sup>−7</sup>) / Outflow")
        w.expect("day 300 outflow", f300["outflow"], r"Inflow 2\.\d+ × 10<sup>−7</sup> / Outflow (3\.\d+ × 10<sup>−7</sup>)\*\*")
        f0 = frame(0.0)
        w.expect("day 0 drain", f0["outflow"], r"The drain reads (2\.\d+ × 10<sup>−7</sup>) m³/s per m at day 0")
        f600 = frame(600 * day)
        w.check(abs(float(f600["outflow"]) / dry_q - 1) < 0.01,
                "by day 600 the drain is %.2f%% from dry weather" % (100 * (float(f600["outflow"]) / dry_q - 1)))


def _num(s):
    try:
        return float(str(s).replace("−", "-"))
    except (TypeError, ValueError):
        return None


def _exit_nodes(st):
    return len((st.data.get("seepage_bc") or {}).get("exit_face") or [])
