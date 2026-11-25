"""LEM-5, A Weak Layer, Non-Circular: built in Studio from File → New, then the
page's runs.

Walked: path C (Studio) — the four soils pasted with LEM and FEM ticked, the
four profile lines, the piezometric line with the page's Add row count, the
four non-circular points typed with the page's Add row count and then cancelled,
and the weak-zone generator on the empty table (its summary and its rows) — the
Run LEM dialog's fixed Surface label and dimmed controls, the Spencer search and
where it left each point, the base split by material, the dimmed circle methods
and the tooltip the page quotes, the crest point pulled in to a 72.4° exit and
the refusal, the circular search from the generated starting circles, and the
generated surface solved alone. Through the engine (the page reports their
results, not the clicks): the two other Movement patterns, the three fixed
shapes of the end-ramp table, and the track elevations. Not walked: path A (the
assistant), path B (the template cell by cell), the five-method search table
(tags), and the ramp angles quoted for each method's surface.
"""
import re

from tutorial_walks.harness import (StepFailed, button, column, editable, field,
                                    labels, page_rows, paste, press, rows,
                                    set_field, show_parameters_for, value)

PAGE = "lem05_weak_layer_noncircular.md"
TIMEOUT = 3600

ENGINE = (("Janbu", "janbu"), ("Corps of Engineers", "corps"), ("Lowe & Karafiath", "lowe"),
          ("Spencer", "spencer"), ("Morgenstern-Price", "mprice"))


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _num(s):
    return float(str(s).replace("−", "-"))


def _single(w, sd, method):
    """A single-surface solve of the model's non-circular surface, or None where
    the method finds no solution on it."""
    try:
        b = w.lem(sd, method, analysis="single_surface", surface="noncircular")
    except Exception:
        return None
    fs = (b.get("results") or {}).get("FS")
    return fs if fs is not None and fs == fs else None


def _base_by_mat(df, mat):
    return float(df.loc[df["mat"] == mat, "dl"].sum())


def walk(w):
    p = w.page
    st = w.studio()
    mats = p.table(r"each soil weighs γ both above and below the water table:")
    lines = [page_rows(p.table(r"\*\*Profile Line %d — material %d" % (i, i))) for i in (1, 2, 3, 4)]
    piezo = page_rows(p.table(r"\*\*Water\*\* — a piezometric line, level at elevation"))
    surf = page_rows(p.table(r"\*\*The failure surface\*\* — four points"))

    with w.step("C: File → New"):
        st.trigger("File", "New")

    with w.step("C1: Materials, profile lines and water"):
        n = p.count(r"press \*\*Add row\*\* (\w+) times and enter \(or copy-paste\) the four soils")

        def m_act(dlg):
            show_parameters_for(dlg, {"LEM", "FEM"})
            press(dlg, "Table view")
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(mats) - 1, "Add row %d times gives %d rows"
                      % (n, t.table.rowCount()))
            paste(t, page_rows(mats), 0, column(t, "name"))
        st.edit("Materials", m_act)
        clay = [m for m in st.data["materials"] if float(m.get("c") or 0) > 0]
        w.check([m["name"] for m in clay] == ["Soft Clay"] and clay[0]["u"] == "none",
                "Soft Clay is not the only row with cohesion and u none")
        depth = p.number(r"set \*\*Max depth \(bottom boundary elevation\):\*\* to `(-?\d+)`")[0]
        names = [r[0] for r in mats[1:]]

        def p_act(dlg):
            set_field(dlg, "Max depth (bottom boundary elevation)", depth)
            for i, pts in enumerate(lines, 1):
                press(dlg, "Add line")
                set_field(dlg, "Material", "%d: %s" % (i, names[i - 1]))
                paste(editable(dlg), pts)
        st.edit("Profile lines", p_act)
        got = [[tuple(map(float, xy)) for xy in l["coords"]] for l in st.data["profile_lines"]]
        w.require(got == [[tuple(map(float, r)) for r in pts] for pts in lines],
                  "the profile lines saved as %s" % got)
        n_pz = p.count(r"on the \*\*Line 1\*\* tab, \*\*Add row\*\* (\w+) for")

        def z_act(dlg):
            for _ in range(n_pz):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(piezo), "Add row %d times gives %d rows"
                      % (n_pz, t.table.rowCount()))
            paste(t, piezo)
        st.edit("Piezometric lines", z_act)
        w.require([tuple(map(float, xy)) for xy in st.data["piezo_line"]] ==
                  [tuple(map(float, r)) for r in piezo], "the piezometric line saved wrong")

    with w.step("C2: The failure surface, typed"):
        n = p.count(r"Press \*\*Add row\*\* (\w+) times and fill them from the table above")
        quoted = " ".join(p.find(r"what a blank Movement means: \*\"(Points ordered .+?)\"\*").group(1).split())
        seen = {}

        def act(dlg):
            seen["help"] = [" ".join(l.text().split()) for l in labels(dlg)]
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(surf), "Add row %d times gives %d rows"
                      % (n, t.table.rowCount()))
            paste(t, surf)
            seen["rows"] = rows(t)
            return "reject"
        st.edit("Non-circular pts", act)
        w.check(any(quoted in h for h in seen["help"]), "the help line does not read %r" % quoted)
        w.check([[_num(r[0]), _num(r[1]), r[2]] for r in seen["rows"]] ==
                [[_num(r[0]), _num(r[1]), r[2]] for r in surf],
                "the typed rows read %s" % seen["rows"])

    with w.step("C3: Generating the surface"):
        quoted = " ".join(p.find(r"under the button: > (Generated .+?) Each part").group(1)
                          .replace(">", " ").split())
        seen = {}

        def act(dlg):
            w.require(editable(dlg).table.rowCount() == 0, "the table is not empty")
            press(dlg, "Generate from the weak zone")
            st.settle(5)
            seen["summary"] = [" ".join(l.text().split()) for l in labels(dlg)
                               if l.text().startswith("Generated")]
        st.edit("Non-circular pts", act)
        w.check(any(quoted == s or quoted in s for s in seen["summary"]),
                "the generator reports %s, not %r" % (seen["summary"], quoted))
        nc = st.data.get("non_circ") or []
        w.require([(round(q["X"], 4), round(q["Y"], 4), q["Movement"]) for q in nc] ==
                  [(round(_num(r[0]), 4), round(_num(r[1]), 4), r[2]) for r in surf],
                  "the generated rows are %s, not the page's table" % nc)

    with w.step("Running the analysis"):
        seen = {}

        def act(dlg):
            seen["surface"] = value(field(dlg, "Surface"))
            for name in ("Composite surfaces (truncate circles at the base)",
                         "Grid search (auto-seed the circular search)"):
                seen[name] = button(dlg, name, visible=False).isEnabled()
            seen["tol"] = dlg.tol.isEnabled()      # the Geometric tol box
            seen["notes"] = [" ".join(l.text().split()) for l in labels(dlg)]
            m = field(dlg, "Method")
            items = {m.itemText(i): m.model().item(i) for i in range(m.count())}
            seen["methods"] = {k: (it.isEnabled(), " ".join(it.toolTip().split())) for k, it in items.items()}
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
        st.run_lem(act)
        w.check(seen["surface"] == "Non-circular", "Surface reads %r, not Non-circular" % seen["surface"])
        w.check(not any(seen[k] for k in ("Composite surfaces (truncate circles at the base)",
                                          "Grid search (auto-seed the circular search)", "tol")),
                "composite surfaces, grid seeding or the geometric tolerance is not dimmed")
        note = " ".join(p.find(r"The dialog's own note says what the two analyses do — \*\"(.+?)\"\*")
                        .group(1).split())
        w.check(any(note in n_ for n_ in seen["notes"]), "the dialog has no note %r" % note)
        tip = " ".join(p.find(r"Hold the pointer over the OMS entry to read why: > (.+?) \*\*OMS and Bishop")
                       .group(1).replace(">", " ").split())
        oms = seen["methods"].get("Ordinary Method of Slices (OMS)", (True, ""))
        bishop = seen["methods"].get("Bishop's Simplified", (True, ""))
        w.check(not oms[0] and not bishop[0], "OMS and Bishop are not dimmed")
        w.check(oms[1] == tip, "the OMS tooltip reads %r, not the page's %r" % (oms[1], tip))
        b = _bundle(st)
        w.expect("search FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\.\*\* Gray is every surface")
        coords = list(b["failure_surface"].coords)
        seam = [c for c in coords if abs(c[1] + 5.8) < 1e-6]
        moved = {"toe end": coords[0], "crest end": coords[-1]}
        if seam:
            moved["seam 1"] = min(seam)
            moved["seam 2"] = max(seam)
        for key, label, nth in (("toe end", "toe end", None), ("seam 1", "seam", 0),
                                ("seam 2", "seam", 1), ("crest end", "crest end", None)):
            if key not in moved:
                w.fail("the searched surface has no vertex on the seam at y = −5.8")
                continue
            x, y = moved[key]
            pat = r"\| %s \| \([^)]*\) \| \((−?\d+\.\d+), (−?\d+\.\d+)\) \|" % label
            w.expect("%s x" % key, x, pat, nth=nth)
            w.expect("%s y" % key, y, pat, nth=nth, group=2)
        df = b["slice_df"]
        w.expect("sliding mass", float(df["w"].sum()), r"surface carries ([\d,]+) lb/ft of soil")
        w.expect("base length", float(df["dl"].sum()), r"lb/ft of soil on a (\d+\.\d) ft base")
        w.expect("base in the clay", _base_by_mat(df, 3), r"\*\*(\d+\.\d) ft in the clay\*\*")
        w.expect("base in the sand", _base_by_mat(df, 2), r"clay\*\*, (\d+\.\d) ft in the sand")
        w.expect("base in the fill", _base_by_mat(df, 1), r"in the sand and (\d+\.\d) ft in the fill")
        w.check(not b["results"].get("warnings"), "Spencer returns warnings: %s" % b["results"].get("warnings"))

    sd = dict(st.data)

    with w.step("What the search found: other Movement patterns (engine)"):
        def moves(pattern):
            nc = [dict(q) for q in sd["non_circ"]]
            for q, mv in zip(nc, pattern):
                q["Movement"] = mv
            return dict(sd, non_circ=nc)
        b = w.lem(moves(["Free"] * 4), "spencer", surface="noncircular")
        w.expect("all four Free", b["results"]["FS"], r"leave behind returns \*\*(1\.\d+)\*\*")
        b = w.lem(moves(["Free", "Fixed", "Fixed", "Free"]), "spencer", surface="noncircular")
        w.expect("interior pair Fixed", b["results"]["FS"], r"the best the search reaches is \*\*(1\.\d+)\*\*")

    with w.step("Limit on the end-ramp angle: the crest point pulled in"):
        x_in = p.number(r"pull the crest point in from x = 39\.1 to x = (\d+)")[0]
        before = float(st.data["non_circ"][-1]["X"])

        def act(dlg):
            t = editable(dlg)
            paste(t, [[x_in]], len(rows(t)) - 1, column(t, "X"))
        st.edit("Non-circular pts", act)
        mark = st.mark()
        n_boxes = len(st.dialogs.boxes)
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Auto search")))
        boxes = st.dialogs.boxes[n_boxes:]
        refusal = p.find(r"a box titled \*\*LEM run failed\*\*: > (.+?) and the Log pane").group(1).strip()
        w.check(any(t == "LEM run failed" and refusal in msg for _k, t, msg in boxes),
                "no 'LEM run failed' box reading %r (boxes: %s)" % (refusal, boxes))
        why = " ".join(p.find(r"and the Log pane shows the reason: > (.+?) Nothing was searched")
                       .group(1).replace(">", " ").split())
        w.check(why in " ".join(st.log_since(mark).split()), "the Log does not print %r" % why)
        st.trigger("Edit", "Undo")
        w.require(abs(float(st.data["non_circ"][-1]["X"]) - before) < 1e-9,
                  "Undo did not put the crest point back (x = %s, was %s)"
                  % (st.data["non_circ"][-1]["X"], before))

    with w.step("Limit on the end-ramp angle: three fixed shapes (engine)"):
        tbl = p.table(r"an 82\.4° exit:")
        xs = (_num(surf[-1][0]), float(x_in),
              float(p.number(r"one pulled in further to x = (\d+\.\d+)")[0]))
        for j, x in enumerate(xs, 1):
            nc = [dict(q) for q in sd["non_circ"]]
            nc[-1]["X"] = x
            shape = dict(sd, non_circ=nc)
            for label, meth in ENGINE:
                fs = _single(w, shape, meth)
                row = next(r for r in tbl[1:] if r[0] == label)
                printed = row[j]
                if printed == "no solution":
                    w.check(fs is None, "%s at the %s exit gives %s where the page has no solution"
                            % (label, tbl[0][j], fs))
                else:
                    pat = r"\| %s \|" % re.escape(label) + r" [^|]+ \|" * (j - 1) + r" ([\d.]+) \|"
                    w.expect("%s at the %s exit" % (label, tbl[0][j]), fs, pat)

    with w.step("Where the track belongs inside the seam (engine)"):
        tbl = p.table(r"what the depth alone changes\.")
        for j in range(1, len(tbl[0])):
            y = _num(tbl[0][j])
            nc = [dict(q) for q in sd["non_circ"]]
            nc[1]["Y"] = nc[2]["Y"] = y
            fs = _single(w, dict(sd, non_circ=nc), "spencer")
            pat = r"\| Factor of safety \(Spencer\) \|" + r" [^|]+ \|" * (j - 1) + r" \**([\d.]+)\** \|"
            w.expect("track at %s" % tbl[0][j], fs, pat)

    with w.step("Comparison with a circular search"):
        def act(dlg):
            press(dlg, "Generate starting circles")
            st.settle(10)
        st.edit("Circles", act)
        w.check(len(st.data["circles"]) == 6, "the generator made %d circles, not six"
                % len(st.data["circles"]))
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Auto search"),
                                set_field(dlg, "Surface", "Circular")))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("circular FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered")
        w.expect("circle x", crit["Xo"], r"centered at \((\d+\.\d+), \d+\.\d+\) with a radius")
        w.expect("circle y", crit["Yo"], r"centered at \(\d+\.\d+, (\d+\.\d+)\) with a radius")
        w.expect("circle radius", crit["Yo"] - crit["Depth"], r"with a radius of (\d+\.\d+), whose")
        w.expect("circle lowest point", crit["Depth"], r"lowest point is at elevation \*\*(−\d\.\d+)\*\*")
        df = b["slice_df"]
        w.expect("circle base", float(df["dl"].sum()), r"of its (\d+\.\d) ft base lies inside")
        w.expect("circle base in the clay", _base_by_mat(df, 3), r"and (\d+\.\d) ft of its")

    with w.step("The generated surface as a starting point"):
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Single surface"),
                                set_field(dlg, "Surface", "Non-circular")))
        w.expect("generated surface alone", _bundle(st)["results"]["FS"],
                 r"\*\*FS = (2\.\d+)\*\* on the fixed surface")
