"""FEM-1, Strength Reduction: the starter run in Studio, limit equilibrium first.

Walked: the starter opened, the Run LEM dialog and the three methods' answers,
FEM mode and the grayed Run FEM, the mesh at the page's settings, the refused
run, the elastic properties typed into the columns the page names, the run's
dialog defaults, and the strength reduction and its trial table. Not walked: the sweeps the page reports but
does not ask the reader to make (E, tolerance, iteration budget, Hybrid, fixed
sides), which the page's tags and the verification pages hold.
"""
import re


from tutorial_walks.harness import (StepFailed, column, editable, field, page_rows,
                                    preflight_items, press, rounds_to, run_button_enabled,
                                    set_field, show_parameters_for, tick, type_cells,
                                    value)

PAGE = "fem01_strength_reduction.md"
TIMEOUT = 3600


def _crit(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b, b["search"]["fs_cache"][0]


def _radius(c):
    return c["R"] if c.get("R") is not None else c["Yo"] - c["Depth"]


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the starter file"):
        st.open(w.copy_model(w.page_link("xslope_ssrm_embankment_start.xlsx")))
        w.check(st.mode == "lem", "the starter opens in %s mode, not LEM" % st.mode)
        c = st.data["circles"][0]
        w.expect("starting circle x", c["Xo"], r"centered at \((\d+), \d+\) with a radius of \d+ ft, a")
        w.expect("starting circle radius", _radius(c), r"with a radius of (\d+) ft, a base-tangent")

    with w.step("The limit equilibrium answer"):
        def act(dlg):
            w.check(value(field(dlg, "Method")) == "Spencer", "Run LEM does not open on Spencer")
            w.check(value(field(dlg, "Analysis")) == "Auto search", "Analysis does not open on Auto search")
            w.check(int(value(field(dlg, "Number of slices"))) == 40, "slices do not open on 40")
            items = preflight_items(dlg) or []
            w.check(any("No problems found for this run" in t for t in items),
                    "Model checks does not read 'No problems found for this run.': %s" % items)
        st.run_lem(act)
        b, sp = _crit(st)
        w.expect("Spencer FS", b["results"]["FS"], r"\*\*Spencer's method gives FS = (1\.\d+)\*\*")
        w.expect("Spencer center x", sp["Xo"], r"toe circle centered at \((\d+\.\d+), \d+\.\d+\)")
        w.expect("Spencer center y", sp["Yo"], r"toe circle centered at \(\d+\.\d+, (\d+\.\d+)\)")
        w.expect("Spencer radius", _radius(sp), r"with a radius of (\d+\.\d+) ft, entering")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.expect("entry x", min(xs), r"entering the crest at x = (\d+\.\d)")
        w.expect("exit x", max(xs), r"leaving at x = (\d+\.\d)")
        for meth, pat in (("Morgenstern-Price", r"Morgenstern-Price \((1\.\d+)\) settles"),
                          ("Bishop's Simplified", r"Bishop's simplified method \((1\.\d+)\)")):
            st.run_lem(lambda dlg, m=meth: set_field(dlg, "Method", m))
            bb, cc = _crit(st)
            w.expect("%s FS" % meth, bb["results"]["FS"], pat)
            if meth == "Morgenstern-Price":
                w.check(abs(cc["Xo"] - sp["Xo"]) < 1e-6 and abs(_radius(cc) - _radius(sp)) < 1e-6,
                        "Morgenstern-Price is not on Spencer's circle")
            else:
                w.check(abs(cc["Xo"] - sp["Xo"]) < 1e-6, "Bishop's center x differs from Spencer's")
                w.expect("Bishop's radius minus Spencer's", _radius(cc) - _radius(sp),
                         r"a radius (\d\.\d+) ft larger")

    with w.step("Building the mesh"):
        st.set_mode("FEM")
        run = st.menu_action("Run", "Run FEM")
        w.check(not run.isEnabled(), "Run FEM… is not disabled before a mesh exists")
        tip = p.find(r"It stays disabled until a mesh exists — \*\"(.+?)\"\*").group(1)
        w.check(tip in (run.toolTip() or ""), "Run FEM…'s tooltip reads %r, not %r" % (run.toolTip(), tip))
        size = p.number(r"set that box to `(\d+\.\d+)`")[0]

        def act(dlg):
            set_field(dlg, "Element type", "Quadratic triangles (tri6)")
            tick(dlg, "Auto-size from geometry", False)
            set_field(dlg, "Target element size", size)
        st.build_mesh(act)
        mesh = st.data.get("mesh")
        w.require(mesh is not None, "no mesh was built")
        w.expect("mesh nodes", len(mesh["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(mesh["elements"]), r"nodes and ([\d,]+) triangles\*\*")

    with w.step("Young's modulus and Poisson's ratio: the refused run"):
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
            seen["run"] = run_button_enabled(dlg)
            return "reject"
        st.trigger("Run", "Run FEM", act)
        head = p.find(r"\*\*Model checks — (\d+ errors · \d+ warning)\*\*, and the \*\*Run\*\* button is disabled").group(1)
        w.check(any(head in t for t in seen["items"]), "Model checks does not read %r: %s" % (head, seen["items"][:2]))
        w.check(seen["run"] is False, "Run is not disabled")

    with w.step("Young's modulus and Poisson's ratio: entering them"):
        tbl = p.table(r"Enter the two values, in the columns headed")
        named = p.find(r"in the columns headed `([^`]+)` and `([^`]+)`").groups()
        w.check(list(named) == tbl[0], "the sentence names %s, the table heads %s" % (named, tbl[0]))

        def act(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"FEM"})
            t = editable(dlg)
            for h, v in zip(tbl[0], tbl[1]):
                type_cells(t, [[v]], 0, column(t, h))
        st.edit("Materials", act)
        m = st.data["materials"][0]
        w.require((float(m["E"]), float(m["nu"])) == (float(tbl[1][0]), float(tbl[1][1])),
                  "E / ν saved as %s / %s" % (m["E"], m["nu"]))
        w.require(st.data.get("mesh") is not None, "entering E and ν dropped the mesh")

    with w.step("Running the strength reduction"):
        seen = {}

        def act(dlg):
            seen["items"] = preflight_items(dlg) or []
            seen["run"] = run_button_enabled(dlg)
            seen["fields"] = {k: value(field(dlg, k)) for k in
                              ("Analysis", "F min (SSRM)", "F max (SSRM)", "Tolerance (SSRM)",
                               "Max iterations per trial", "Iteration ceiling", "Failure criterion")}
            seen["k0"] = value(field(dlg, "K0 initial stress")) if _has(dlg, "K0 initial stress") else None
        _s, log = st.run_fem(act)
        w.check(any("1 warning" in t for t in seen["items"]) and seen["run"],
                "Model checks does not read 1 warning with Run enabled: %s" % seen["items"][:2])
        cut = p.find(r"312\.5 / tan 20° = \*\*(\d+\.\d) psf\*\*").group(1)
        w.check(any(cut in t for t in seen["items"]), "the warning does not name %s psf" % cut)
        f = seen["fields"]
        w.check(f["Analysis"].startswith("SSRM"), "Analysis opens on %s" % f["Analysis"])
        for k, want in (("F min (SSRM)", 1.0), ("F max (SSRM)", 2.0), ("Tolerance (SSRM)", 0.01),
                        ("Max iterations per trial", 12000), ("Iteration ceiling", 50000)):
            w.check(abs(float(str(f[k]).replace(",", "")) - want) < 1e-9, "%s opens on %s" % (k, f[k]))
        w.check(f["Failure criterion"] == "Non-convergence", "Failure criterion opens on %s" % f["Failure criterion"])
        r = (st.results.get("fem_solution") or {}).get("ssrm_result")
        w.require(r is not None, "the run left no strength-reduction result")
        w.expect("FS", r["FS"], r"\*\*FS = (1\.\d+)\*\*, the midpoint of the final bracket")
        lo, hi = r["final_interval"]
        w.expect("bracket low", lo, r"final bracket \[(1\.\d+), 1\.\d+\]")
        w.expect("bracket high", hi, r"final bracket \[1\.\d+, (1\.\d+)\]")
        tbl = p.table(r"The Log reports each trial:")
        trials = r.get("trials") or []
        w.check(len(trials) == len(tbl) - 1, "the run made %d trials, the table has %d" % (len(trials), len(tbl) - 1))
        for t, row in zip(trials, tbl[1:]):
            ok_f = rounds_to(float(t["F"]), row[1])[0]
            conv = bool(t.get("converged"))
            ok_v = conv == (row[2] == "converged")
            ok_i = int(t.get("iterations") or 0) == int(row[3].replace(",", ""))
            w.check(ok_f and ok_v and ok_i, "trial %s: the run gives F %.4f %s in %s iterations; the "
                    "page prints %s" % (row[0], float(t["F"]), "converged" if conv else "failed",
                                        t.get("iterations"), row[1:4]))


def _has(dlg, label):
    try:
        field(dlg, label)
        return True
    except StepFailed:
        return False
