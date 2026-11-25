"""LEM-1, Simple Embankment: built in Studio from File → New, then run.

Walked: path C (Studio) step by step, the Spencer and Bishop runs and the Log
block the page quotes, the dry and water-filled tension crack, and the
minimum-crack-depth sweep (through the engine, as the page leaves the tool to
the reader). Not walked: path A (the assistant: a provider call, see W-1) and
path B (the Excel template cell by cell; the completed workbook it produces is
the one the page's tags run).
"""
from tutorial_walks.harness import (StepFailed, editable, field, press, rows,
                                    set_field, show_parameters_for, type_cells,
                                    value)

PAGE = "lem01_simple_embankment.md"
TIMEOUT = 1800


def _fs(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b["results"]["FS"], b


def _strip(st):
    sv = st.mw.solution_canvas
    if sv is None:
        return ""
    lab = sv._warning
    return lab.text() if lab.isVisibleTo(sv) else ""


def _lem(method):
    def act(dlg):
        set_field(dlg, "Method", method)
        set_field(dlg, "Analysis", "Auto search")
    return act


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("C: File → New"):
        st.trigger("File", "New")
        w.require(st.doc.is_open, "File → New opened no project")

    with w.step("C1: Global parameters, Units Imperial"):
        def glob(dlg):
            set_field(dlg, "Units", "Imperial")
            w.check(abs(float(value(field(dlg, "Unit weight of water"))) - 62.4) < 1e-9,
                    "the unit weight of water does not fill itself with 62.4")
        st.edit("Global parameters", glob)
        w.require(str(st.data.get("unit_system")).lower().startswith("imp"),
                  "Units did not become Imperial (%r)" % st.data.get("unit_system"))

    with w.step("C2: Materials, one soil in List view"):
        def mats(dlg):
            w.check(press is not None and _checked(dlg, "Table view"),
                    "the Materials editor does not open on Table view")
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "List view")
            press(dlg, "Add")
            set_field(dlg, "Name", "soil")
            set_field(dlg, "γ", "125")
            set_field(dlg, "Model (option)", "mc")
            set_field(dlg, "c (psf)", "500")
            set_field(dlg, "φ", "0")
            set_field(dlg, "Model (u)", "none")
        st.edit("Materials", mats)
        m = st.data["materials"]
        w.require(len(m) == 1, "Materials holds %d rows, not one" % len(m))
        m = m[0]
        w.require((m["name"], float(m["gamma"]), m["option"], float(m["c"]),
                   float(m["phi"]), m["u"]) == ("soil", 125.0, "mc", 500.0, 0.0, "none"),
                  "the material saved as %r" % ({k: m.get(k) for k in
                                                 ("name", "gamma", "option", "c", "phi", "u")},))

    with w.step("C3: Profile lines, Add line, Add row, three points"):
        n_rows = p.count(r"\*\*Add row\*\* (\w+)(?: times?)?(?:, for \w+ rows,)? and enter `0, 0`")
        pts = [(0, 0), (20, 20), (60, 20)]

        def prof(dlg):
            press(dlg, "Add line")
            set_field(dlg, "Max depth (bottom boundary elevation)", "0")
            set_field(dlg, "Material", "1: soil")
            for _ in range(n_rows):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(pts),
                      "Add row %d time(s) on a new line gives %d rows for the page's %d "
                      "points" % (n_rows, t.table.rowCount(), len(pts)))
            type_cells(t, pts)
        st.edit("Profile lines", prof)
        lines = st.data.get("profile_lines") or []
        w.require(len(lines) == 1, "%d profile lines saved" % len(lines))
        got = [tuple(map(float, xy)) for xy in lines[0]["coords"]]
        w.require(got == [tuple(map(float, xy)) for xy in pts],
                  "the profile line saved as %s" % (got,))
        w.require(float(st.data.get("max_depth")) == 0.0, "Max depth did not save as 0")

    with w.step("C4: Circles, Generate starting circles"):
        seen = {}

        def circ(dlg):
            press(dlg, "Generate starting circles")
            st.settle(10)
            from PySide6.QtWidgets import QLabel
            seen["summary"] = [l.text() for l in dlg.findChildren(QLabel)
                               if l.text().startswith("Generated")]
        st.edit("Circles", circ)
        w.check(any("not daylighting inside the model" in s for s in seen["summary"]),
                "the generator does not report the candidate it dropped for not "
                "daylighting inside the model: %s" % seen["summary"])
        c = st.data.get("circles") or []
        w.require(len(c) == 1, "the generator left %d circles" % len(c))
        c = c[0]
        w.check((float(c["Xo"]), float(c["Yo"])) == (10.0, 40.0),
                "the generated center is (%s, %s), not (10, 40)" % (c["Xo"], c["Yo"]))
        w.check(abs(float(c.get("R") or (float(c["Yo"]) - float(c["Depth"]))) - 40.0) < 1e-6,
                "the generated circle is not R = 40")

    with w.step("Running the analysis: Spencer, Auto search"):
        def run_dlg(dlg):
            w.check(value(field(dlg, "Method")) == "Spencer",
                    "Run LEM does not open on Spencer")
            w.check(int(value(field(dlg, "Number of slices"))) == 40,
                    "the slice count does not default to 40")
            set_field(dlg, "Analysis", "Auto search")
        _secs, log = st.run_lem(run_dlg)
        fs, b = _fs(st)
        w.expect("Spencer's FS", fs, r"\*\*FS = (1\.\d+)\*\*\. In Studio")
        strip = _strip(st)
        w.check(strip != "", "no amber warning strip after the Spencer run")
        block = p.code_block(r"Studio's \*\*Log\*\* pane show this")
        for line in block:
            w.check(line.strip() in log, "the Log does not print the page's line %r" % line.strip())
        crit = b["search"]["fs_cache"][0]
        w.expect("Spencer's radius", crit["R"] if "R" in crit else crit["Yo"] - crit["Depth"],
                 r"radius 33\.4 against Spencer's (\d+\.\d)")
        xs = [x for x, _ in b["failure_surface"].coords]
        w.expect("Spencer's exit x", max(xs), r"rather than x = (\d+)")
        w.check(abs(crit["Xo"] - 10) < 1e-6 and abs(crit["Yo"] - 40) < 1e-6,
                "the center left (10, 40): (%.2f, %.2f)" % (crit["Xo"], crit["Yo"]))

    with w.step("Comparing with Bishop's method"):
        st.run_lem(_lem("Bishop's Simplified"))
        fs, b = _fs(st)
        w.expect("Bishop's FS", fs, r"Bishop finds \*\*FS = (1\.\d+)\*\*")
        crit = b["search"]["fs_cache"][0]
        w.expect("Bishop's radius", crit["R"] if "R" in crit else crit["Yo"] - crit["Depth"],
                 r"radius (\d+\.\d) against Spencer's")
        xs = [x for x, _ in b["failure_surface"].coords]
        w.expect("Bishop's exit x", max(xs), r"crest at x = (\d+) rather")

    with w.step("Adding a tension crack (dry, 8 ft)"):
        depth = p.number(r"set \*\*Tension crack depth\*\* = `(\d+)`")[0]

        def crack(dlg):
            set_field(dlg, "Tension crack depth", depth)
            set_field(dlg, "Water in crack", "0")
        st.edit("Global parameters", crack)
        _s, log = st.run_lem(_lem("Spencer"))
        fs, _b = _fs(st)
        w.expect("cracked Spencer FS", fs, r"\*\*FS = (1\.\d+)\*\*, and this time the solution is clean")
        w.check(_strip(st) == "", "the cracked run still shows the amber strip: %r" % _strip(st))
        n_bad, n_all = p.number(r"Log still reports (\d+) of the \d+ trial")[0], \
            p.number(r"Log still reports \d+ of the (\d+) trial")[0]
        w.check("could not solve %s of %s trial surfaces" % (n_bad, n_all) in log,
                "the Log does not report %s of %s unsolved" % (n_bad, n_all))
        for meth in ("Bishop's Simplified", "Morgenstern-Price"):
            st.run_lem(_lem(meth))
            fs2, _b = _fs(st)
            w.expect("cracked %s FS" % meth, fs2, r"\*\*FS = (1\.\d+)\*\*, and this time the solution is clean")

    with w.step("Filling the crack with water"):
        wd = p.number(r"set \*\*Water in crack\*\* = `(\d+)`")[0]
        st.edit("Global parameters", lambda dlg: set_field(dlg, "Water in crack", wd))
        st.run_lem(_lem("Spencer"))
        fs, b = _fs(st)
        w.expect("water-filled Spencer FS", fs, r"\*\*FS = (1\.\d+)\*\*, 6\.5% below")
        df = b["slice_df"]
        sig = (df["n_eff"] / df["dl"]).values
        neg = [s for s in sig if s < 0]
        w.check(len(neg) == 1, "%d base-stress bars are red, not one" % len(neg))
        if neg:
            w.expect("tension on the slice against the crack", neg[0],
                     r"about (−\d+) psf of\s+tension")
        w.expect("P_w on the crack face", float(df["t"].values[-1]),
                 r"times 8\^2 = (\d\{,\}\d+) \\text")
        for meth in ("Bishop's Simplified", "Morgenstern-Price"):
            st.run_lem(_lem(meth))
            fs2, _b = _fs(st)
            w.expect("water-filled %s FS" % meth, fs2, r"\*\*FS = (1\.\d+)\*\*, 6\.5% below")
        st.edit("Global parameters", lambda dlg: set_field(dlg, "Water in crack", "0"))
        w.require(float(st.data.get("tcrack_water") or 0) == 0.0, "Water in crack did not go back to 0")

    with w.step("Finding the minimum crack depth"):
        sd = dict(st.data)
        sd.pop("mesh", None)
        for depth, pattern in (("4", r"(\d+)% of the interior boundaries are in tension at a 4 ft"),
                               ("6", r"(\d+)% at 6 ft"), ("7", r"(\d+)% at 7 ft"),
                               ("7.4", r"still (\d+)% at 7\.4 ft")):
            b = w.lem(dict(sd, tcrack_depth=float(depth), tcrack_water=0.0), "spencer")
            pct = _tension_pct(b)
            w.expect("tension at a %s ft crack (%%)" % depth, pct, pattern)
        b = w.lem(dict(sd, tcrack_depth=7.5, tcrack_water=0.0), "spencer")
        w.check(_tension_pct(b) in (None, 0), "the 7.5 ft crack still has tension")
        w.expect("FS at a 7.5 ft crack", b["results"]["FS"], r"factor of safety at 7½ ft is \*\*(1\.\d+)\*\*")


def _checked(dlg, text):
    from tutorial_walks.harness import button
    return button(dlg, text).isChecked()


def _tension_pct(bundle):
    """The share of interior boundaries in tension, as the warning strip words it."""
    import re
    for wtxt in (bundle["results"].get("warnings") or []):
        m = re.search(r"(\d+(?:\.\d+)?)% of", str(wtxt))
        if m and "tension" in str(wtxt).lower() and "base" not in str(wtxt).lower():
            return float(m.group(1))
    return None
