"""SEEP-2, Unconfined Seepage Through a Zoned Dam: built in Studio from File → New.

Walked: the units and time unit; the three material rows pasted after the page's
Add row count; the three profile lines; the two heads and the exit face from the
page's tables; the built model held against the completed workbook; the mesh and
its Log size line; the run, its quoted Log block, the iteration count, the Base
material the Display panel opens on, the discharge, the head and pore-pressure
ranges and the largest suction; the convergence-tolerance table; the van
Genuchten and Gardner runs of the three-model table; the run that does not
converge (its Log lines, status bar and plot title) and the same run at 400
iterations; the core-conductivity table. Not walked: the section integrals of
"Where the water goes" and the flow-net counts of "Scaling the flow net" (read
off the solved field by the figure producer, not by anything Studio prints).
"""
import re

import numpy as np

from tutorial_walks.harness import (column, editable, field, page_rows, paste,
                                    preflight_items, press, seep_bc_set, select_row,
                                    set_field, set_value, show_parameters_for, tick,
                                    type_cells, value)

PAGE = "seep02_johnson_dam.md"
TIMEOUT = 3600


def _sol(st):
    return st.results["seep_solutions"][1]["solution"]


def _iters(log):
    m = re.findall(r"Converged in (\d+) iterations", log)
    return int(m[-1]) if m else None


def _wet(log):
    m = re.findall(r"(\d+)/(\d+) exit face active", log)
    return (int(m[-1][0]), int(m[-1][1])) if m else (None, None)


def _set_cols(st, cells):
    """Materials on Seepage alone: ``cells`` is {column: [value per row]}."""
    def act(dlg):
        from PySide6.QtWidgets import QComboBox
        press(dlg, "Table view")
        show_parameters_for(dlg, {"Seepage"})
        t = editable(dlg)
        for h, vals in cells.items():
            c = column(t, h)
            for r, v in enumerate(vals):
                cw = t.table.cellWidget(r, c)
                if cw is not None:
                    set_value(cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox), v)
                else:
                    type_cells(t, [[v]], r, c)
    st.edit("Materials", act)


def _run(st, tol=None, max_iter=None):
    def act(dlg):
        if tol is not None:
            set_field(dlg, "Convergence tol", tol)
        if max_iter is not None:
            set_field(dlg, "Max iterations", max_iter)
    return st.run_seep(act)


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"top of the section first — and fill the three rows:")
    lines = [page_rows(p.table(r"\*\*Line %d — `%s`\*\*" % (i, n))) for i, n in
             ((1, "shell"), (2, "core"), (3, "foundation"))]

    with w.step("Building the model: File → New, Seepage"):
        st.trigger("File", "New")
        st.set_mode("Seepage")

    with w.step("1. Global parameters"):
        def act(dlg):
            set_field(dlg, "Units", "Imperial")
            w.check(abs(float(value(field(dlg, "Unit weight of water"))) - 62.4) < 1e-9,
                    "the unit weight of water does not fill itself with 62.4")
            set_field(dlg, "Time", "day")
        st.edit("Global parameters", act)

    with w.step("2. Materials"):
        n = p.count(r"Press \*\*Add row\*\* (\w+) times — the zones go in the order")

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"Seepage"})
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(mats) - 1, "Add row %d times gives %d rows" % (n, t.table.rowCount()))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", act)
        got = [(m["name"], float(m["k1"]), float(m["k2"]), m["unsat"], float(m["kr0"]), float(m["h0"]))
               for m in st.data["materials"]]
        want = [(r[0], float(r[1]), float(r[2]), r[4], float(r[5]), float(r[6])) for r in mats[1:]]
        w.require(got == want, "the materials saved as %s" % got)

    with w.step("3. Profile lines"):
        depth = p.number(r"set \*\*Max depth \(bottom boundary elevation\):\*\* to `(\d+)`")[0]
        names = [m["name"] for m in st.data["materials"]]

        def act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", depth)
            for i, pts in enumerate(lines):
                press(dlg, "Add line")
                set_field(dlg, "Material", "%d: %s" % (i + 1, names[i]))
                paste(editable(dlg), pts)
        st.edit("Profile lines", act)
        got = [[tuple(map(float, c)) for c in l["coords"]] for l in st.data["profile_lines"]]
        w.require(got == [[tuple(map(float, r)) for r in pts] for pts in lines],
                  "the profile lines saved as %s" % got)

    with w.step("4. Boundary conditions"):
        h1 = (p.number(r"\*\*Head 1 — the reservoir\.\*\* Press \*\*Add head\*\*, leave \*\*Type:\*\* at `head`, set \*\*Head value \(ft\):\*\* to `(\d+)`")[0],
              page_rows(p.table(r"\*\*Head 1 — the reservoir\.\*\*")))
        h2 = (p.number(r"\*\*Head value \(ft\):\*\* = `(\d+)`, applied along the downstream")[0],
              page_rows(p.table(r"\*\*Head 2 — the tailwater\.\*\*")))
        ex = page_rows(p.table(r"\*\*Exit face\.\*\* Select the \*\*Exit face\*\* entry"))

        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            for val, pts in (h1, h2):
                press(s1, "Add head")
                set_field(s1, "Head value (ft)", val)
                paste(editable(s1), pts)
            select_row(s1.list, "Exit face")
            paste(editable(s1), ex)
        st.edit("Seep BC", act)
        done = w.load(w.copy_model(w.page_link("xslope_johnson_res.xlsx"), "done"))
        for k in ("seepage_bc", "profile_lines", "max_depth"):
            w.check(_same(st.data.get(k), done.get(k)), "the built model's %s differs from the completed file's" % k)
        for a, b in zip(st.data["materials"], done["materials"]):
            w.check(all(_same(a.get(k), b.get(k)) for k in ("k1", "k2", "alpha", "unsat", "kr0", "h0")),
                    "material %s differs from the completed file's" % a["name"])

    with w.step("Building the mesh"):
        div = p.number(r"Set \*\*Size divisions\*\* to `(\d+)`")[0]

        def act(dlg):
            set_field(dlg, "Element type", "Linear triangles (tri3)")
            tick(dlg, "Auto-size from geometry", True)
            set_field(dlg, "Size divisions", div)
        _s, log = st.build_mesh(act)
        w.check("Auto element size: 6.250" in log, "the Log does not state the 6.250 ft size")
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(m["elements"]), r"nodes and ([\d,]+) triangles\*\*, with")

    with w.step("Running the analysis"):
        seen = {}

        def act(dlg):
            w.check(abs(float(value(field(dlg, "Convergence tol"))) - 0.0001) < 1e-12, "Convergence tol is not 0.0001")
            seen["items"] = preflight_items(dlg) or []
        _s, log = st.run_seep(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]), "Model checks: %s" % seen["items"][:3])
        w.check_log(p.code_block(r"The Log pane opens with the path the solver took:"), log, "unconfined run")
        w.expect("iterations", _iters(log), r"The run finishes in \*\*(\d+) iterations\*\*")
        canvas = st.mw.seep_solution_canvas.get(1)
        panel = st.mw._display_panels.get(canvas)
        w.check(panel is not None and panel.base_mat.currentText() == "3: foundation",
                "Base material does not read 3: foundation")
        s = _sol(st)
        w.expect("discharge", s["flowrate"], r"\*\*The total discharge is (\d\.\d+) ft³/day per ft\*\*")
        w.expect("lowest head", float(np.min(s["head"])), r"The head ranges from \*\*(\d+\.\d+) ft to")
        w.expect("highest head", float(np.max(s["head"])), r"ft to (\d+\.\d+) ft\*\*, which is the two")
        u = np.asarray(s["u"], float)
        w.expect("lowest pore pressure", u.min(), r"The pore pressure runs from \*\*(−[\d,]+\.\d) psf to")
        w.expect("highest pore pressure", u.max(), r"psf to ([\d,]+\.\d) psf\*\*\.")
        nn = np.asarray(st.data["mesh"]["nodes"], float)
        psi = np.asarray(s["head"], float) - nn[:, 1]
        k = int(np.argmin(psi))
        w.expect("largest suction", -psi[k], r"The largest suction is \*\*(\d+\.\d+) ft\*\*")
        w.expect("suction x", nn[k, 0], r"−\d+\.\d+ ft — at \((\d+\.\d), \d+\.\d\)")
        w.expect("suction y", nn[k, 1], r"−\d+\.\d+ ft — at \(\d+\.\d, (\d+\.\d)\)")

    with w.step("The three conditions: the tolerance table"):
        anchor = r"it turns out to change the answer hardly at all:"
        tbl = p.table(anchor)
        start = p.find(anchor).end()
        for r in tbl[1:]:
            _s, log = _run(st, tol=r[0])
            at = r"\| %s \| [\d.]+ \| " % re.escape(r[0])
            w.expect("q at tol %s" % r[0], _sol(st)["flowrate"], at + r"(\d\.\d+) \|", nth=0, start=start)
            w.expect("iterations at tol %s" % r[0], _iters(log), at + r"\d\.\d+ \| (\d+) \|", nth=0, start=start)
        _run(st, tol="0.0001")

    with w.step("Running all three unsaturated models"):
        par = p.table(r"means the three-way comparison compares models rather than soils\.")
        n_rows = len(st.data["materials"])
        for model, a_col, n_col in (("vg", 3, 4), ("gard", 5, 6)):
            _set_cols(st, {"unsat": [model] * n_rows, "vg_a": [r[a_col] for r in par[1:]],
                           "vg_n": [r[n_col] for r in par[1:]]})
            _s, log = _run(st)
            at = r"\| `%s` \| " % model
            w.expect("q, %s" % model, _sol(st)["flowrate"], at + r"(\d\.\d+) \|")
            w.expect("iterations, %s" % model, _iters(log), at + r"\d\.\d+ \| (\d+) \|")

    with w.step("A run that does not converge"):
        kr0 = p.number(r"with `kr0` = `([\d.]+)` and `h0`")[0]
        h0 = p.find(r"and `h0` = `(-?\d+)` on all three materials").group(1)
        _set_cols(st, {"unsat": ["lf"] * n_rows, "kr0": [kr0] * n_rows, "h0": [h0] * n_rows})
        w.require(st.data.get("mesh") is not None, "the materials edit dropped the mesh")
        cap = p.number(r"Set \*\*Max iterations\*\* on the \*\*Run Seepage\*\* dialog to `(\d+)`")[0]
        _s, log = _run(st, max_iter=cap)
        w.check(_sol(st).get("converged") is False, "the %s-iteration run converged" % cap)
        w.check_log(p.code_block(r"uses all 50 iterations and then stops:"), log, "the run that stops")
        w.expect("discharge at the limit", _sol(st)["flowrate"], r"returns one, (\d\.\d+) here")
        quoted = " ".join(p.find(r"the status bar reads \*(Seepage \(BC set 1\) did not converge: .+?)\*, and the plot").group(1).split())
        w.check(quoted in " ".join(st.status().split()), "the status bar reads %r" % st.status())
        title = p.find(r"the plot title reads \*\*(Seepage Solution \(did not converge\))\*\*").group(1)
        w.check(title in _titles(st), "the plot title is not %r (%s)" % (title, _titles(st)))
        _s, log = _run(st, max_iter="400")
        w.expect("iterations at 400", _iters(log), r"the same run converges, in \*\*(\d+) iterations\*\*")
        w.expect("q at 400", _sol(st)["flowrate"], r"iterations\*\*, at \*q\* = \*\*(\d\.\d+)\*\*")
        _set_cols(st, {"kr0": [r[5] for r in mats[1:]], "h0": [r[6] for r in mats[1:]]})

    with w.step("Raising the core's conductivity"):
        anchor = r"two cells on one row — \*\*k1\*\* and \*\*k2\*\*"
        tbl = p.table(anchor)
        start = p.find(anchor).end()
        for r in tbl[1:]:
            def act(dlg, k=r[0]):
                press(dlg, "Table view")
                show_parameters_for(dlg, {"Seepage"})
                t = editable(dlg)
                for h in ("k1 (ft/day)", "k2 (ft/day)"):
                    type_cells(t, [[k]], 1, column(t, h))
            st.edit("Materials", act)
            _s, log = _run(st)
            at = r"\| %s \| " % re.escape(r[0])
            w.expect("q at core k %s" % r[0], _sol(st)["flowrate"], at + r"(\d+\.\d+) \|", nth=0, start=start)
            wet, of = _wet(log)
            w.check("%s of %s" % (wet, of) == r[2], "at core k %s %s/%s exit-face nodes are wet, not %s"
                    % (r[0], wet, of, r[2]))
            w.expect("iterations at core k %s" % r[0], _iters(log),
                     at + r"[\d.]+ \| [^|]+ \| [^|]+ \| (\d+) \|", nth=0, start=start)


def _titles(st):
    c = st.mw.seep_solution_canvas.get(1)
    fig = getattr(c, "figure", None) or getattr(getattr(c, "canvas", None), "figure", None)
    if fig is None:
        return []
    out = [a.get_title() for a in fig.axes]
    if getattr(fig, "_suptitle", None) is not None:
        out.append(fig._suptitle.get_text())
    return out


def _f(v):
    try:
        return float(str(v).replace("−", "-"))
    except (TypeError, ValueError):
        return None


def _same(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        keys = {k for k in set(a) | set(b) if a.get(k) not in (None, [], "") or b.get(k) not in (None, [], "")}
        return all(_same(a.get(k), b.get(k)) for k in keys)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    fa, fb = _f(a), _f(b)
    return a == b or (fa is not None and fb is not None and abs(fa - fb) < 1e-9)
