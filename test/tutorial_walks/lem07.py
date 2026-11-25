"""LEM-7, Strength Options Beyond Mohr-Coulomb: the three parts in Studio, each
on a copy of the page's file.

Walked: Part A — the power-curve material read in List view, the Spencer search
and the Log line the page quotes, the linear fit entered and searched, the two
fits' crossing and their strengths at 10 kPa, and the London clay pair. Part B —
the cp material read, the Bishop search at 50 slices and the share of its
surface in the lowest clay, the constant strength entered and searched, and the
four constants of the comparison table (through the engine: the page reports
them, not the clicks). Part C — both clays entered as ratios in List view and
searched. Not walked: the equation in the materials plot's title (a picture),
and the hand-worked strength table of Part C (arithmetic on the page's inputs).
"""
import math
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, _process, field, labels, press,
                                    set_field, show_parameters_for, value)

PAGE = "lem07_strength_envelopes.md"
TIMEOUT = 3600


def _bundle(st):
    b = st.results.get("lem_solution")
    if not b or not b.get("results"):
        raise StepFailed("the run left no solution")
    return b


def _run(method, slices=None):
    def act(dlg):
        set_field(dlg, "Method", method)
        set_field(dlg, "Analysis", "Auto search")
        if slices is not None:
            set_field(dlg, "Number of slices", slices)
    return act


def _list_view(dlg, row=None):
    """The materials editor on LEM alone, in List view, on material ``row``."""
    from studio.editors import _MaterialListView
    show_parameters_for(dlg, {"LEM"})
    press(dlg, "List view")
    lv = dlg.findChild(_MaterialListView)
    if lv is None:
        raise StepFailed("the materials editor has no List view")
    if row is not None:
        lv.list.setCurrentRow(row)
        _process()
    return lv


def _mean_sigma(df):
    """The length-weighted effective normal stress along the slice bases."""
    return float(df["n_eff"].sum() / df["dl"].sum())


def walk(w):
    p = w.page
    st = w.studio()

    # ------------------------------------------------------------------ Part A
    with w.step("Part A: opening the model"):
        st.open(w.copy_model(w.page_link("xslope_baker_clay.xlsx")))
        seen = {}

        def act(dlg):
            _list_view(dlg)
            seen["option"] = value(field(dlg, "Model (option)"))
            seen["a"] = value(field(dlg, "pow_a"))
            seen["b"] = value(field(dlg, "pow_b"))
            return "reject"
        st.edit("Materials", act)
        w.check(seen["option"] == "pow", "Model (option) reads %r, not pow" % seen["option"])
        w.expect("pow_a", float(seen["a"]), r"here `pow_a` = (\d\.\d+)")
        w.expect("pow_b", float(seen["b"]), r"and `pow_b` = (\d\.\d+) with")

    with w.step("Part A: running the analysis"):
        _s, log = st.run_lem(_run("Spencer"))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("power-curve FS", b["results"]["FS"], r"\*\*FS = (0\.\d+)\*\* — below one")
        w.expect("center x", crit["Xo"], r"on a circle centered at \((−\d+\.\d+), \d+\.\d+\) tangent")
        w.expect("center y", crit["Yo"], r"on a circle centered at \(−\d+\.\d+, (\d+\.\d+)\) tangent")
        w.expect("tangent elevation", crit["Depth"], r"tangent at elevation (−\d+\.\d+), \d")
        quoted = p.find(r"— `(Sliding mass = [^`]+)` in Studio's Log").group(1)
        w.check(quoted in log, "the Log does not print %r" % quoted)
        w.expect("mean effective normal stress", _mean_sigma(b["slice_df"]),
                 r"the effective normal stress along its slice bases averages (\d+\.\d) kPa")
        pow_sigma = _mean_sigma(b["slice_df"])

    with w.step("Part A: entering the linear fit"):
        fit = p.table(r"In the `mat` worksheet these are three adjacent cells:")

        def act(dlg):
            _list_view(dlg)
            set_field(dlg, "Model (option)", fit[1][0])
            set_field(dlg, "c (kPa)", fit[1][1])
            set_field(dlg, "φ", fit[1][2])
        st.edit("Materials", act)
        m = st.data["materials"][0]
        w.require((m["option"], float(m["c"]), float(m["phi"])) == ("mc", 11.64, 24.7),
                  "the material saved as %s / %s / %s" % (m["option"], m["c"], m["phi"]))
        st.run_lem(_run("Spencer"))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        w.expect("linear-fit FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\* — Slide reports 1\.536")
        w.expect("center x", crit["Xo"], r"on a circle centered at\s?\((−\d+\.\d+), \d+\.\d+\), tangent")
        w.expect("center y", crit["Yo"], r"on a circle centered at\s?\(−\d+\.\d+, (\d+\.\d+)\), tangent")
        w.expect("tangent elevation", crit["Depth"], r"tangent to elevation (\d+\.\d+) and")
        w.expect("surface length", float(b["slice_df"]["dl"].sum()), r"and (\d+\.\d+) m long")
        w.expect("sliding mass", float(b["slice_df"]["w"].sum()), r"carrying (\d+\.\d) kN/m of soil\. That")
        w.expect("mean effective normal stress", _mean_sigma(b["slice_df"]),
                 r"(\d+\.\d) along the Mohr-Coulomb one")

    with w.step("Part A: where the two fits cross"):
        a, n = 1.107, 0.86
        c, phi = 11.64, math.radians(24.7)
        grid = np.linspace(1.0, 300.0, 299001)
        diff = a * grid ** n - (c + grid * math.tan(phi))
        # The first crossing: the straight line starts above the curve (its
        # cohesion) and the curve overtakes it there.
        i = int(np.argmax(np.sign(diff[1:]) != np.sign(diff[:-1])))
        cross = float(grid[i])
        w.expect("crossing stress", cross, r"The two envelopes cross at σ′ = (\d+\.\d) kPa")
        w.expect("power curve at 10 kPa", a * 10 ** n, r"At 10 kPa the curve gives (\d+\.\d+) kPa")
        w.expect("line at 10 kPa", c + 10 * math.tan(phi), r"where the line gives (\d+\.\d+), and")
        w.check(pow_sigma < cross, "the power curve's surface is not in the low-stress range")

    with w.step("Part A: London clay, where the two fits agree"):
        for name, pat in (("vp061a.xlsx", r"\*\*FS = (1\.\d+)\*\*, against Slide's 1\.468"),
                          ("vp061b.xlsx", r"\*\*FS = (1\.\d+)\*\*, against Slide's 1\.366")):
            st.open(w.copy_model(w.page_link(name)))
            if st.mode != "lem":
                # The page has to say so: Run LEM is not on the Run menu in another mode.
                w.require(p.finditer(r"switch the mode strip to \*\*LEM\*\*"),
                          "%s opens in %s mode and the page does not say to switch to LEM" % (name, st.mode))
                st.set_mode("LEM")
            st.run_lem(_run("Spencer"))
            w.expect("London clay %s" % name, _bundle(st)["results"]["FS"], pat)

    # ------------------------------------------------------------------ Part B
    with w.step("Part B: opening the model"):
        st.open(w.copy_model(w.page_link("xslope_low_clay.xlsx"), "partB"))
        seen = {}

        def act(dlg):
            _list_view(dlg, 2)
            for k in ("Model (option)", "c (kPa)", "r-elev", "c/p"):
                seen[k] = value(field(dlg, k))
            return "reject"
        st.edit("Materials", act)
        w.check(seen["Model (option)"] == "cp", "the third row's option reads %r" % seen["Model (option)"])
        w.expect("c", float(seen["c (kPa)"]), r"`c` = (\d+) is the strength at the reference")
        w.expect("r-elev", float(seen["r-elev"]), r"`r-elev` = (\d+) is that")
        w.expect("c/p", float(seen["c/p"]), r"`c/p` = (\d+\.\d+) is the rate")

    with w.step("Part B: running the analysis"):
        slices = p.number(r"raising \*\*Number of slices\*\* to `(\d+)`")[0]
        st.run_lem(_run("Bishop's Simplified", slices))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        df = b["slice_df"]
        w.expect("cp FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, against Low's published")
        w.expect("center x", crit["Xo"], r"The circle is centered at \((\d+\.\d+), \d+\.\d+\)")
        w.expect("center y", crit["Yo"], r"The circle is centered at \(\d+\.\d+, (\d+\.\d+)\)")
        w.expect("surface length", float(df["dl"].sum()), r"\), (\d+\.\d+) m of surface carrying")
        w.expect("sliding mass", float(df["w"].sum()), r"m of surface carrying (\d+\.\d) kN/m")
        w.expect("lowest point", crit["Depth"], r"its lowest point is at \*\*elevation (\d\.\d+)\*\*")
        low = df["mat"] == 3
        in_low = float(df.loc[low, "dl"].sum())
        w.expect("length in the lowest clay", in_low, r"Of the [\d.]+ m of surface, (\d+\.\d+) m lies")
        mobilized = float((df.loc[low, "c"] * df.loc[low, "dl"]).sum() / in_low)
        w.expect("mobilized strength", mobilized, r"that stretch averages (\d+\.\d+) kPa")
        cp_run = b

    with w.step("Part B: replacing the profile with a constant strength"):
        const = p.table(r"a typical choice for a single value:")

        def act(dlg):
            _list_view(dlg, 2)
            set_field(dlg, "Model (option)", const[1][0])
            set_field(dlg, "c (kPa)", const[1][1])
            set_field(dlg, "φ", const[1][2])
        st.edit("Materials", act)
        st.run_lem(_run("Bishop's Simplified", slices))
        b = _bundle(st)
        df = b["slice_df"]
        w.expect("constant FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, and the critical circle now sits")
        w.expect("constant lowest point", b["search"]["fs_cache"][0]["Depth"],
                 r"sits \*\*tangent to elevation (\d)\*\*")
        w.expect("constant surface length", float(df["dl"].sum()), r"It is (\d+\.\d+) m long against")
        w.expect("constant length in the lowest clay", float(df.loc[df["mat"] == 3, "dl"].sum()),
                 r"and (\d+\.\d+) m of it runs through the lowest clay")

    with w.step("Part B: four constants (engine)"):
        tbl = p.table(r"we run the same edit at four different constants:")
        base = dict(st.data)
        for r in tbl[1:]:
            if r[1] in ("—", "-"):
                continue
            mats = [dict(m) for m in base["materials"]]
            mats[2].update(option="mc", c=float(r[1]), phi=0.0)
            b = w.lem(dict(base, materials=mats), "bishop", num_slices=int(slices))
            at = r"\| %s \| %s \| " % (re.escape(r[0]), re.escape(r[1]))
            w.expect("Bishop at %s kPa" % r[1], b["results"]["FS"], at + r"\**(\d\.\d+)\** \|")
            w.expect("tangent at %s kPa" % r[1], b["search"]["fs_cache"][0]["Depth"],
                     at + r"\**\d\.\d+\** \| \**(\d\.\d+)\** \|")
            if r[1] == "22.88":
                longer = float(b["slice_df"]["dl"].sum()) - float(cp_run["slice_df"]["dl"].sum())
                more = float(b["slice_df"]["w"].sum()) - float(cp_run["slice_df"]["w"].sum())
                w.expect("lengthens by", longer, r"drops to the base, lengthens by (\d\.\d) m")
                w.expect("more soil", more, r"and picks up (\d+) kN/m more soil")
        mats = [dict(m) for m in base["materials"]]
        mats[2].update(option="mc", c=30.0, phi=0.0)
        b = w.lem(dict(base, materials=mats), "bishop", num_slices=int(slices))
        w.expect("constant from the bottom", b["results"]["FS"], r"one taken from the bottom gives (1\.\d+)")

    # ------------------------------------------------------------------ Part C
    with w.step("Part C: entering both clays as ratios"):
        st.open(w.copy_model(w.page_link("xslope_low_clay.xlsx"), "partC"))
        ratio = p.table(r"Change \*\*Model \(option\)\*\* to `ratio` and enter these values:")
        seen = {}

        def act(dlg):
            for row in (2, 1):
                _list_view(dlg, row)
                set_field(dlg, "Model (option)", ratio[1][0])
                set_field(dlg, "c (kPa)", ratio[1][1])
                set_field(dlg, "c/p", ratio[1][3])
                seen.setdefault("phi_shown", []).append(
                    any(l.text() == "φ" for l in labels(dlg)))
        st.edit("Materials", act)
        w.check(not any(seen["phi_shown"]), "the ratio form still shows φ")
        got = [(m["option"], float(m["c"]), float(m.get("cp") or 0)) for m in st.data["materials"][1:3]]
        w.require(got == [("ratio", 15.0, 0.1875)] * 2, "the two clays saved as %s" % got)

    with w.step("Part C: searching with both clays as ratios"):
        _s, log = st.run_lem(_run("Bishop's Simplified", slices))
        b = _bundle(st)
        crit = b["search"]["fs_cache"][0]
        df = b["slice_df"]
        w.expect("ratio FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, up from 1\.142")
        w.expect("ratio lowest point", crit["Depth"], r"its lowest point is at elevation (\d\.\d+),")
        w.expect("ratio surface length", float(df["dl"].sum()), r"and it is (\d+\.\d+) m long")
        quoted = p.find(r"The Log prints that elevation as `([^`]+)`").group(1)
        w.check(quoted in log, "the Log does not print %r" % quoted)
        low = df["mat"] == 3
        in_low = float(df.loc[low, "dl"].sum())
        w.expect("ratio length in the lower clay", in_low, r"Of its length, (\d+\.\d+) m runs")
        w.expect("ratio mobilized strength", float((df.loc[low, "c"] * df.loc[low, "dl"]).sum() / in_low),
                 r"where the strength averages (\d+\.\d+) kPa")
