"""COMBO-3, Factor of Safety vs Time: both parts in Studio.

Walked, Part 1: the starter opened with its mesh and transient solution, the
nineteen frames and the mesh counts, the strengths typed into Materials (and
what that OK does to the loaded solution), the two starting circles and the
search limit, the Run LEM dialog's Seepage time group and the two single-frame
runs with the Log lines the page quotes, and the Parametric sweep over every
frame with its table. Part 2: the Johnson dam opened in FEM mode, the three
editors the page opens and cancels, and both sweeps (rapid drawdown, then
drained) with their tables. Not read: the figures' own readings of where the
circles break out (they are the plotted circles the tables already give).
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, button, column, editable, field,
                                    preflight_items, rows, set_field, set_value,
                                    show_parameters_for, value)

PAGE = "combo03_fs_vs_time.md"
TIMEOUT = 3 * 3600


def _tabs(st):
    vt = st.mw.view_tabs
    return [vt.tabText(i) for i in range(vt.count())]


def _frames(st):
    return (st.results.get("transient_seep") or {}).get("frames") or []


def _fmt(t):
    return ("%g" % t)


def walk(w):
    p = w.page
    st = w.studio()

    # ------------------------------------------------------------------ Part 1
    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_earth_dam_fs_time_start.xlsx")))
        w.check(st.mode == "seep", "Studio opens the starter in %s mode, not Seepage" % st.mode)

    with w.step("The mesh and the seepage solution shipped with the file"):
        w.require(any("Transient" in t for t in _tabs(st)),
                  "there is no Seep · Transient tab (tabs: %s)" % _tabs(st))
        times = [float(f["time"]) for f in _frames(st)]
        listed = p.find(r"steps through \*\*\w+\*\* instants — t = ([\d, and]+) — and every run").group(1)
        want = [float(x) for x in re.findall(r"\d+", listed)]
        w.check(times == want, "the saved frames are at %s, not the page's %s" % (times, want))
        w.check(len(times) == p.count(r"steps through \*\*(\w+)\*\* instants"),
                "the play bar steps through %d instants" % len(times))
        mesh = st.data["mesh"]
        w.expect("mesh nodes", len(mesh["nodes"]), r"solved on \*\*([\d,]+) nodes and", nth=0)
        w.expect("mesh triangles", len(mesh["elements"]), r"nodes and ([\d,]+) triangles\*\*", nth=0)
        shipped = [np.asarray(f["head"], float).copy() for f in _frames(st)]

    with w.step("Entering the soil properties"):
        st.set_mode("LEM")
        typed = {}
        for name in ("shell", "core"):
            m = p.find(r"\*\*γ (\d+)\*\*, \*\*γsat (\d+)\*\*, \*\*c (\d+)\*\*, \*\*φ (\d+)\*\* and "
                       r"\*\*u = (\w+)\*\* for the `%s`" % name)
            typed[name] = dict(zip(("γ", "γsat", "c", "φ", "u"), m.groups()))

        def act(dlg):
            from PySide6.QtWidgets import QComboBox
            press_table = button(dlg, "Table view")
            press_table.click()
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            hdr_rows = rows(t)
            for r, row in enumerate(hdr_rows):
                nm = next((n for n in typed if n in row), None)
                if nm is None:
                    continue
                for key, val in typed[nm].items():
                    c = column(t, key)
                    cw = t.table.cellWidget(r, c)
                    cb = cw if isinstance(cw, QComboBox) else (cw.findChild(QComboBox) if cw else None)
                    if cb is not None:
                        set_value(cb, val)
                    else:
                        it = t.table.item(r, c)
                        if it is None:
                            raise StepFailed("cell %s of %s cannot be typed into" % (key, nm))
                        it.setText(val)
        st.edit("Materials", act)
        mats = {m["name"]: m for m in st.data["materials"]}
        for nm, vals in typed.items():
            m = mats.get(nm) or {}
            got = (float(m.get("gamma") or 0), float(m.get("gamma_sat") or 0), float(m.get("c") or 0),
                   float(m.get("phi") or 0), m.get("u"))
            want = (float(vals["γ"]), float(vals["γsat"]), float(vals["c"]), float(vals["φ"]), vals["u"])
            w.require(got == want, "%s saved as %s, not %s" % (nm, got, want))
        kept = bool(_frames(st)) and any("Transient" in t for t in _tabs(st))
        if p.finditer(r"the seepage solution the file loaded is gone from the session"):
            w.check(not kept, "the page says this OK clears the loaded seepage solution and "
                              "closes the Seep · Transient tab; Studio keeps both")
        else:
            w.require(kept, "the Materials OK dropped the loaded seepage solution (tabs %s)" % _tabs(st))
            w.check(p.says("stays loaded"), "the page no longer says the solution stays loaded")

    with w.step("Starting circles"):
        tbl = p.table(r"Click \*\*Circles\*\*. The file has two")
        seen = {}

        def act(dlg):
            seen["rows"] = rows(editable(dlg))
            return "reject"
        st.edit("Circles", act)
        for i, r in enumerate(tbl[1:]):
            c = st.data["circles"][i]
            w.check((float(c["Xo"]), float(c["Yo"]), float(c["Depth"])) ==
                    (float(r[0]), float(r[1]), float(r[3])),
                    "circle %d is (%s, %s, %s), not the page's %s" % (i + 1, c["Xo"], c["Yo"], c["Depth"], r))
        win = p.table(r"The file sets one:")
        for r in win[1:]:
            key = {"Min slip depth": "min_slip_depth"}.get(r[0])
            w.check(key is not None and float((st.data.get("search_window") or {}).get(key) or 0)
                    == float(r[1]), "the search window does not set %s = %s" % (r[0], r[1]))
        w.check(bool(_frames(st)), "Cancel on Circles dropped the seepage solution")

    if p.finditer(r"### Solving the seepage again"):
        with w.step("Solving the seepage again"):
            st.set_mode("Seepage")
            st.run_seep(lambda dlg: w.check(value(field(dlg, "Run type")).startswith("Transient"),
                                            "Run type does not read Transient"))
            st.set_mode("LEM")

    def at_frame(label, notes=None):
        def act(dlg):
            w.check(value(field(dlg, "Method")) == "Spencer", "Run LEM does not open on Spencer")
            w.check(value(field(dlg, "Analysis")) == "Auto search", "Analysis does not open on Auto search")
            w.check(int(value(field(dlg, "Number of slices"))) == 40, "slices do not open on 40")
            w.check(value(field(dlg, "Time")) == "Saved frame", "Time does not open on Saved frame")
            if notes is not None:
                notes.append(value(field(dlg, "Saved frame")))
            set_field(dlg, "Saved frame", label)
            if notes is not None:
                notes.append(preflight_items(dlg) or [])
        return act

    with w.step("Stability at one time step: t = 0"):
        notes = []
        _s, log = st.run_lem(at_frame("t = 0 day", notes))
        w.check(notes[0] == "t = %s day" % _fmt(_frames(st)[-1]["time"]),
                "the frame list opens on %r, not the last saved time step" % notes[0])
        quoted = " ".join(p.find(r"shows one note, which restates the dependency: > (.+?) Click \*\*Run\*\*")
                          .group(1).replace(">", " ").split())
        w.check(any(quoted in " ".join(t.split()) for t in notes[1]),
                "Model checks does not show the page's note: %s" % notes[1])
        w.check_log(p.code_block(r"Click \*\*Run\*\*. The Log reports the frame it read"), log, "t = 0")
        loads = st.data.get("dloads_derived") or []
        peak = max((float(pt["Normal"]) for blk in loads for pt in blk), default=None)
        if peak is not None:
            w.expect("water load at the heel", peak, r"blue arrows reach \*\*(\d+\.\d) kPa\*\*")

    with w.step("Stability at one time step: t = 5"):
        _s, log = st.run_lem(at_frame("t = 5 day"))
        w.check_log(p.code_block(r"set \*\*Saved frame\*\* to `t = 5 day` and click \*\*Run\*\*"), log, "t = 5")

    def sweep(rapid):
        def act(dlg):
            set_field(dlg, "Mode", "Factor of safety vs time")
            w.check(value(field(dlg, "Method")) == "Spencer", "Parametric does not open on Spencer")
            lw = dlg.times
            from PySide6.QtCore import Qt
            rb = button(dlg, "Rapid drawdown at each time")
            if rapid is None:
                w.check(not rb.isEnabled(), "Rapid drawdown at each time is not grayed out")
            else:
                w.check(rb.isEnabled(), "Rapid drawdown at each time is grayed out")
                set_value(rb, rapid)
            ticked = [lw.item(i).checkState() == Qt.Checked for i in range(lw.count())]
            if rapid:
                # t = 0 is stage 1 of every drawdown: it unticks and dims.
                w.check(not ticked[0] and all(ticked[1:]),
                        "with Rapid drawdown ticked, t = 0 is not the one frame unticked")
            else:
                w.check(all(ticked), "not every saved frame is ticked")
            rs = button(dlg, "Re-search the critical surface at each step")
            w.check(rs.isChecked(), "Re-search the critical surface at each step is not ticked")
        return act

    with w.step("Sweeping all time steps"):
        _s, log = st.run_parametric(sweep(None))
        w.check_log(p.code_block(r"a row per instant giving the factor of safety and the circle"),
                    log, "Part 1 sweep")

    # ------------------------------------------------------------------ Part 2
    with w.step("Part 2: opening the model"):
        st.doc.mark_clean() if hasattr(st.doc, "mark_clean") else setattr(st.doc, "_dirty", False)
        st.open(w.copy_model(w.page_link("xslope_johnson_fs_time.xlsx")))
        w.check(st.mode == "fem", "Studio opens the file in %s mode, not FEM" % st.mode)
        st.set_mode("LEM")
        w.require(bool(_frames(st)), "the file opened without its transient solution")

    with w.step("Part 2: what the file contains"):
        tbl = p.table(r"free-draining through the drawdown")

        def mats(dlg):
            button(dlg, "Table view").click()
            show_parameters_for(dlg, {"LEM"})
            return "reject"
        st.edit("Materials", mats)
        by = {m["name"]: m for m in st.data["materials"]}
        for r in tbl[1:]:
            m = by.get(r[1], {})
            got = [m.get("gamma"), m.get("gamma_sat"), m.get("c"), m.get("phi"), m.get("d"), m.get("psi")]
            w.check([float(x or 0) for x in got] == [float(x) for x in r[2:8]] and m.get("u") == r[8],
                    "%s reads %s / u %s, not the page's %s" % (r[1], got, m.get("u"), r[2:]))
        st.edit("Circles", lambda dlg: "reject")
        c = st.data["circles"][0]
        w.check((float(c["Xo"]), float(c["Yo"])) == (275.0, 235.0) and len(st.data["circles"]) == 1,
                "the file's circle is not the page's single one at (275, 235)")
        st.set_mode("Seepage")
        st.edit("Transient", lambda dlg: "reject")
        ts = st.data.get("tseep") or {}
        w.check(float(ts.get("stage_1") or 0) == float(p.number(r"\*\*Stage 1 time\*\* is `(\d+)`")[0]),
                "Stage 1 time is %s" % ts.get("stage_1"))
        w.check(float(ts.get("stage_2") or 0) == float(p.number(r"\*\*Stage 2 time\*\* is `(\d+)`")[0]),
                "Stage 2 time is %s" % ts.get("stage_2"))
        st.set_mode("LEM")
        w.require(bool(_frames(st)), "Cancel on the three editors dropped the transient solution")

    with w.step("Running the rapid drawdown sweep"):
        _s, log = st.run_parametric(sweep(True))
        w.check_log(p.code_block(r"table appears when the last search finishes"),
                    log, "rapid drawdown sweep")

    with w.step("The same time steps with drained strengths"):
        _s, log = st.run_parametric(sweep(False))
        w.check_log(p.code_block(r"so the sweep takes a fraction of the time"), log, "drained sweep")
