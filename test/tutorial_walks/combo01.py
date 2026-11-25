"""COMBO-1, Seepage → LEM → FEM: the three runs in Studio on a copy of the file.

Walked: the file opened with no mesh (Run LEM refused for the missing mesh); the
mesh at the dialog's defaults; Run LEM refused for the missing field; the steady
seepage run, its iteration count, discharge and head range, and the two files it
writes beside the workbook; the u column on LEM alone and an OK that changes
nothing; u set to none on all three rows (the warning, and the search's answer),
then Edit → Undo; the Spencer search, its Log lines and the critical circle's
geometry and pore pressures; the strength reduction at the dialog's defaults.
Not walked: the iteration lines of the quoted seepage Log block, which the page
says differ by a few between computers (the iteration count is held "about").
"""
import os
import re

import numpy as np

from tutorial_walks.harness import (button, column, editable, field, headers,
                                    preflight_items, press, rows, run_button_enabled,
                                    set_field, set_value, show_parameters_for, value)

PAGE = "combo01_seepage_stability.md"
TIMEOUT = 2 * 3600


def _set_u(st, u):
    def act(dlg):
        from PySide6.QtWidgets import QComboBox
        press(dlg, "Table view")
        show_parameters_for(dlg, {"LEM"})
        t = editable(dlg)
        c = column(t, "u")
        for r in range(t.table.rowCount()):
            cw = t.table.cellWidget(r, c)
            set_value(cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox), u)
    st.edit("Materials", act)


def _checks(st, menu_text="Run LEM"):
    seen = {}

    def act(dlg):
        seen["items"] = preflight_items(dlg) or []
        seen["run"] = run_button_enabled(dlg)
        return "reject"
    st.trigger("Run", menu_text, act)
    return seen


def walk(w):
    p = w.page
    st = w.studio()
    path = w.copy_model(w.page_link("xslope_johnson_res.xlsx"))

    with w.step("The dam: the file as first opened"):
        st.open(path)
        w.check(st.data.get("mesh") is None, "the file opens with a mesh")
        st.set_mode("LEM")
        seen = _checks(st)
        quoted = p.find(r"\*\"(takes pore pressure from a seepage solution \(u = seep\), but this model carries no mesh)\"\*").group(1)
        w.check(any(quoted in " ".join(t.split()) for t in seen["items"]) and seen["run"] is False,
                "Run LEM on the file as opened is not refused for the missing mesh: %s" % seen["items"][:4])

    with w.step("One mesh for all three analyses"):
        st.set_mode("Seepage")

        def act(dlg):
            w.check(value(field(dlg, "Element type")) == "Quadratic triangles (tri6)", "Element type is not tri6")
            w.check(button(dlg, "Auto-size from geometry").isChecked(), "Auto-size from geometry is not ticked")
            w.check(int(float(value(field(dlg, "Size divisions")))) == 100, "Size divisions is not 100")
            tes = field(dlg, "Target element size")
            w.check(not tes.isEnabled() and abs(float(value(tes)) - 1.0) < 1e-9,
                    "Target element size is not the grayed 1.000")
        _s, log = st.build_mesh(act)
        w.check("Auto element size: 7.500" in log, "the Log does not state the 7.500 ft size the build uses")
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(m["elements"]), r"nodes and ([\d,]+) triangles\*\*, with the seepage")
        st.set_mode("LEM")
        seen = _checks(st)
        w.check(any("carries no pore-pressure field" in t for t in seen["items"]) and seen["run"] is False,
                "with a mesh and no solution Run LEM is not refused for the missing field: %s" % seen["items"][:4])
        st.set_mode("Seepage")

    with w.step("Solving it"):
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
        _s, log = st.run_seep(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]), "Model checks: %s" % seen["items"][:3])
        its = re.findall(r"Converged in (\d+) iterations", log)
        w.expect_about("iterations", int(its[-1]) if its else None, r"settles in about \*\*(\d+) iterations\*\*", 3)
        s = st.results["seep_solutions"][1]["solution"]
        w.expect("discharge", s["flowrate"], r"\*\*The total discharge is (\d\.\d+) ft³/day per ft\*\*")
        w.expect("lowest head", float(np.min(s["head"])), r"The head runs from (\d+\.\d+) ft to")
        w.expect("highest head", float(np.max(s["head"])), r"ft to (\d+\.\d+) ft, the two boundary")
        here = os.path.dirname(path)
        for name in re.findall(r"`(xslope_johnson_res_(?:mesh\.json|seep\.csv))`", p.flat):
            w.check(os.path.exists(os.path.join(here, name)), "%s is not written beside the workbook" % name)

    with w.step("The pore-pressure option (u column)"):
        st.set_mode("LEM")
        seen = {}

        def look(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            c = column(t, "u")
            i = headers(t)[0].index("u")
            seen["u"] = [r[i] for r in rows(t)]
            from PySide6.QtWidgets import QComboBox
            cw = t.table.cellWidget(0, c)
            cb = cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox)
            seen["items"] = [cb.itemText(k) for k in range(cb.count())]
        st.edit("Materials", look)
        want = [r[0] for r in p.table(r"so a section can mix them:")[1:]]
        w.check(seen["items"] == want, "the u column offers %s, not %s" % (seen["items"], want))
        w.check(seen["u"] == ["seep"] * 3, "the three materials read %s, not seep" % seen["u"])
        w.check(st.results.get("seep_solutions"), "an OK that changed nothing dropped the seepage solution")
        _set_u(st, "none")
        seen = {}

        def run_none(dlg):
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
            seen["items"] = preflight_items(dlg) or []
        st.run_lem(run_none)
        quoted = p.find(r"raises a \*\*warning\*\* — \*\"(This model carries a solved seepage field, but no material takes its pore pressure from it)\"\*").group(1)
        w.check(any(quoted in " ".join(t.split()) for t in seen["items"]),
                "u = none draws no warning that no material reads the field: %s" % seen["items"][:4])
        w.expect("FS with u = none", st.results["lem_solution"]["results"]["FS"],
                 r"the Spencer search below returns \*\*FS = (1\.\d+)\*\*")
        st.trigger("Edit", "Undo")
        w.require([m["u"] for m in st.data["materials"]] == ["seep"] * 3, "Undo did not put u back to seep")
        w.require(st.results.get("seep_solutions"), "the u edits dropped the seepage solution")

    with w.step("The limit equilibrium run"):
        seen = {}

        def act(dlg):
            w.check(value(field(dlg, "Method")) == "Spencer", "Run LEM does not open on Spencer")
            seen["items"] = preflight_items(dlg) or []
            set_field(dlg, "Analysis", "Auto search")
        _s, log = st.run_lem(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]), "Model checks: %s" % seen["items"][:3])
        w.check_log(p.code_block(r"Its last two refinement steps and its closing lines read:"), log, "Spencer search")
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("Spencer FS", b["results"]["FS"], r"\*\*Spencer's method gives FS = (1\.\d+)\*\*")
        w.expect("radius", crit["Yo"] - crit["Depth"], r"with a radius of (\d+\.\d+) ft\. It enters")
        coords = list(b["failure_surface"].coords)
        left = min(coords, key=lambda c: c[0])
        w.expect("entry x", left[0], r"enters the upstream face at x = (\d+\.\d), elevation")
        w.expect("entry elevation", left[1], r"x = \d+\.\d, elevation (\d+\.\d) — above")
        w.expect("bottom elevation", min(c[1] for c in coords), r"bottoms at elevation (\d+\.\d)")
        w.expect("exit x", max(c[0] for c in coords), r"comes out on the foundation surface at x = (\d+\.\d)")
        w.expect("candidate circles", len(b["search"]["fs_cache"]), r"The search evaluated (\d+) candidate circles")
        u = b["slice_df"]["u"].to_numpy(dtype=float)
        w.expect("largest slice-base pore pressure", u.max(), r"it runs from 0 to ([\d,]+) psf")
        w.check(int((u <= 0).sum()) == 2, "%d slices carry no pore pressure, not two" % int((u <= 0).sum()))

    with w.step("The finite element run"):
        st.set_mode("FEM")
        w.require(st.menu_action("Run", "Run FEM").isEnabled(), "Run FEM… is not available at once")
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
            seen["run"] = run_button_enabled(dlg)
            seen["f"] = [value(field(dlg, k)) for k in ("Analysis", "F min (SSRM)", "F max (SSRM)", "Tolerance (SSRM)")]
        st.run_fem(act, timeout=2 * 3600)
        w.check(seen["run"] and any("No problems found for this run" in t for t in seen["items"]),
                "Model checks: %s" % seen["items"][:3])
        w.check(str(seen["f"][0]).startswith("SSRM") and [float(x) for x in seen["f"][1:]] == [1.0, 2.0, 0.01],
                "the dialog opens on %s" % seen["f"])
        r = (st.results.get("fem_solution") or {}).get("ssrm_result")
        w.require(r is not None, "the run left no strength-reduction result")
        w.expect("SSRM FS", r["FS"], r"\*\*Strength reduction gives FS = (1\.\d+)\*\*")
