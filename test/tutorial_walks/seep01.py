"""SEEP-1, Confined Seepage Under a Sheetpile: built in Studio from File → New.

Walked: path C (Studio) — the units and time unit, the material row, the
profile line with the page's Add row count, the two head blocks with theirs —
the mesh and its Log lines, the run and its Log lines, the discharge and the
head and pore-pressure ranges, the mesh study at five sizes, the five element
types at 1 m, the refinement at the sheetpile toe, and the Parametric design
sweep on k1. Not walked: path A (the assistant), path B (the template cell by
cell; its workbook is the completed file the page's tags run), the flow-net
counts and the gradient readings of the two reading sections.
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    preflight_items, press, seep_bc_set, select_row,
                                    set_field, show_parameters_for, tick, type_cells,
                                    value)

PAGE = "seep01_sheetpile.md"
TIMEOUT = 3600


def _q(st):
    return st.results["seep_solutions"][1]["solution"]["flowrate"]


def walk(w):
    p = w.page
    st = w.studio()
    mat = p.table(r"\*\*Material\*\* — one soil, isotropic")
    verts = page_rows(p.table(r"two base corners the maximum depth supplies"))

    with w.step("C: File → New, Seepage mode"):
        st.trigger("File", "New")
        st.set_mode("Seepage")
        w.check(st.menu_action("Run", "Run Seep") is not None, "the Run menu has no Run Seep…")

    with w.step("C1: Global parameters"):
        def act(dlg):
            set_field(dlg, "Units", "SI")
            w.check(abs(float(value(field(dlg, "Unit weight of water"))) - 9.81) < 1e-9,
                    "the unit weight of water does not fill itself with 9.81")
            set_field(dlg, "Time", "yr")
        st.edit("Global parameters", act)
        w.require(str(st.data.get("time_unit")) == "yr", "Time saved as %r" % st.data.get("time_unit"))

    with w.step("C2: Materials"):
        def act(dlg):
            show_parameters_for(dlg, {"Seepage"})
            press(dlg, "Add row")
            t = editable(dlg)
            for h, v in zip(mat[0][:4], mat[1][:4]):
                hh = {"k1": "k1 (m/yr)", "k2": "k2 (m/yr)"}.get(h, h)
                type_cells(t, [[v]], 0, column(t, hh))
        st.edit("Materials", act)
        m = st.data["materials"][0]
        w.require((m["name"], float(m["k1"]), float(m["k2"]), float(m.get("alpha") or 0)) ==
                  ("soil", 30.0, 30.0, 0.0), "the material saved as %s" % m)
        w.check(str(m.get("unsat") or "lf") == mat[1][4], "unsat reads %r" % m.get("unsat"))

    with w.step("C3: Profile lines"):
        n = p.count(r"or \*\*Add row\*\* (\w+) times and type them")

        def act(dlg):
            press(dlg, "Add line")
            set_field(dlg, "Max depth (bottom boundary elevation)", "0")
            set_field(dlg, "Material", "1: soil")
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(verts), "Add row %d times gives %d rows for %d vertices"
                      % (n, t.table.rowCount(), len(verts)))
            type_cells(t, verts)
        st.edit("Profile lines", act)
        got = [tuple(map(float, xy)) for xy in st.data["profile_lines"][0]["coords"]]
        w.require(got == [tuple(map(float, r)) for r in verts], "the line saved as %s" % got)

    with w.step("C4: Boundary conditions"):
        bcs = p.table(r"\*\*Boundary conditions\*\* — two specified-head blocks")
        n = p.count(r"3\. \*\*Add row\*\* (\w+), and enter")

        def act(dlg):
            s1 = seep_bc_set(dlg, 0)
            for r in bcs[1:]:
                press(s1, "Add head")
                set_field(s1, "Head value (m)", r[2])
                for _ in range(n):
                    press(s1, "Add row")
                t = editable(s1)
                w.require(t.table.rowCount() == 2, "Add row %d times on a new head gives %d rows"
                          % (n, t.table.rowCount()))
                pts = [r[3].strip("()").split(","), r[4].strip("()").split(",")]
                type_cells(t, [[a.strip() for a in pt] for pt in pts])
        st.edit("Seep BC", act)
        heads = st.data["seepage_bc"]["specified_heads"]
        w.require(len(heads) == 2 and [float(h["head"]) for h in heads] == [13.0, 10.0],
                  "the head blocks saved as %s" % heads)
        w.require(not st.data["seepage_bc"].get("exit_face"), "the exit face is not empty")

    def mesh(etype="Linear triangles (tri3)", auto=None, size=None, divisions=None, refine=None,
             factor=None):
        def act(dlg):
            set_field(dlg, "Element type", etype)
            if auto is not None:
                tick(dlg, "Auto-size from geometry", auto)
            if size is not None:
                set_field(dlg, "Target element size", size)
            if divisions is not None:
                set_field(dlg, "Size divisions", divisions)
            if refine is not None:
                tick(dlg, "Refine near features", refine)
            if factor is not None:
                set_field(dlg, "Refinement factor", factor)
        return act

    with w.step("Building the mesh"):
        def act(dlg):
            w.check(value(field(dlg, "Element type")) == "Quadratic triangles (tri6)",
                    "Build Mesh does not open on tri6")
            mesh()(dlg)
            w.check(int(float(value(field(dlg, "Size divisions")))) == 100, "Size divisions is not 100")
        _s, log = st.build_mesh(act)
        w.check_log(p.code_block(r"The Log pane reports the size of what it built"), log, "mesh")
        m = st.data["mesh"]
        w.expect("nodes", len(m["nodes"]), r"\*\*([\d,]+) nodes and [\d,]+ triangles\*\* — the element")
        w.expect("triangles", len(m["elements"]), r"\*\*[\d,]+ nodes and ([\d,]+) triangles\*\* — the element")

    with w.step("Running the analysis"):
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
        _s, log = st.run_seep(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]),
                "Model checks: %s" % seen["items"][:3])
        w.check_log(p.code_block(r"The Log pane reports which path the solver took"), log, "run")
        s = st.results["seep_solutions"][1]["solution"]
        w.expect("discharge", s["flowrate"], r"\*\*The total discharge is (\d+\.\d+) m³/yr per m\.\*\*")
        w.expect("lowest head", float(np.min(s["head"])), r"The head ranges from \*\*(\d+\.\d+) m to")
        w.expect("highest head", float(np.max(s["head"])), r"m to (\d+\.\d+) m\*\*, which is the two")
        w.expect("largest pore pressure", float(np.max(s["u"])), r"Pore pressure runs from \*\*0 to (\d+\.\d) kPa\*\*")
        w.expect("smallest pore pressure", float(np.min(s["u"])), r"Pore pressure runs from \*\*(\d+) to")

    with w.step("How fine the mesh has to be"):
        tbl = p.table(r"type each size into \*\*Target element size\*\* — and run each one:")
        for r in tbl[1:]:
            st.build_mesh(mesh(auto=False, size=r[0], refine=False))
            m = st.data["mesh"]
            st.run_seep(lambda dlg: None)
            at = r"\| %s \| " % re.escape(r[0])
            w.expect("nodes at %s m" % r[0], len(m["nodes"]), at + r"([\d,]+) \|")
            w.expect("triangles at %s m" % r[0], len(m["elements"]), at + r"[\d,]+ \| ([\d,]+) \|")
            w.expect("q at %s m" % r[0], _q(st), at + r"[\d,]+ \| [\d,]+ \| (\d+\.\d+) \|")
        st.build_mesh(mesh(auto=True, divisions="120"))
        st.run_seep(lambda dlg: None)
        w.expect("q at width/120", _q(st), r"the discharge is \*\*(\d+\.\d+)\*\*, and it lands")
        tbl = p.table(r"with each of the five element types gives:")
        names = {"tri3": "Linear triangles (tri3)", "tri6": "Quadratic triangles (tri6)",
                 "quad4": "Linear quads (quad4)", "quad8": "Quadratic quads (quad8)",
                 "quad9": "Quadratic quads (quad9)"}
        for r in tbl[1:]:
            st.build_mesh(mesh(names[r[0]], auto=False, size="1.0"))
            m = st.data["mesh"]
            st.run_seep(lambda dlg: None)
            at = r"\| `%s` \| " % r[0]
            w.expect("%s nodes" % r[0], len(m["nodes"]), at + r"([\d,]+) \|")
            w.expect("%s elements" % r[0], len(m["elements"]), at + r"[\d,]+ \| ([\d,]+) \|")
            w.expect("%s q" % r[0], _q(st), at + r"[\d,]+ \| [\d,]+ \| (\d+\.\d+) \|")

    with w.step("Refining at the sheetpile toe"):
        st.build_mesh(mesh(auto=True, divisions="100", refine=True, factor="4"))
        m = st.data["mesh"]
        st.run_seep(lambda dlg: None)
        w.expect("refined nodes", len(m["nodes"]), r"The whole mesh grew from [\d,]+ nodes to ([\d,]+)")
        w.expect("refined q", _q(st), r"\*\*(\d+\.\d+) m³/yr per m\*\*, against [\d.]+ on the uniform")
        nn = np.asarray(m["nodes"], float)
        d = np.sort(np.hypot(nn[:, 0] - 30.0, nn[:, 1] - 7.0))
        w.expect("toe's nearest neighbor", d[1], r"from \*\*0\.500 m away to (\d\.\d+) m\*\*")

    with w.step("Discharge against conductivity"):
        st.build_mesh(mesh(auto=True, divisions="100", refine=False))

        def act(dlg):
            set_field(dlg, "Mode", "Design (q target)")
            set_field(dlg, "Material / BC", "soil")
            set_field(dlg, "Property", "k1")
            set_field(dlg, "From", "30")
            set_field(dlg, "To", "300")
            set_field(dlg, "Steps", "10")
            set_field(dlg, "Target q", "100")
        _s, log = st.run_parametric(act)
        w.check_log(p.code_block(r"On this mesh that is under a second\. It reports:"), log, "sweep")
