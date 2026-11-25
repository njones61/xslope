"""LEM-11, Reliability: the whole page in Studio, on a copy of the page's file.

Walked: the derived water loads, the two standard deviations read in List
view, the deterministic search and the starting circle alone, the Reliability
dialog on the Taylor series (its defaults and the σ list it shows) and the
per-parameter table, Monte Carlo at the dialog's settings and again with random
sampling, the response surface and its checks, the β-against-P_f table, the
Parametric variance Pareto and rank-correlation runs, σ(c) halved in List view,
and both engines again. Not walked: the rank correlations' values and the
Pareto's bars (the plots; the shares are checked from the Taylor table they are
computed from), and the "fuller test" of the response surface against thirty
thousand real solves (the page reports it; the reader does not run it).
"""
import re

from tutorial_walks.harness import (StepFailed, button, field, labels, press,
                                    set_field, set_value, show_parameters_for,
                                    value)

PAGE = "lem11_reliability.md"
TIMEOUT = 2 * 3600


def _rel(st):
    b = st.results.get("reliability")
    if not b or not b.get("reliability"):
        raise StepFailed("the run left no reliability result")
    return b["reliability"]


def _sigma_edit(dlg, label):
    """The ± σ box beside the list-view field ``label``."""
    from PySide6.QtWidgets import QLabel
    from studio.editors import _MaterialListView
    lv = dlg.findChild(_MaterialListView)
    for key, cell in lv._cell_widgets.items():
        labs = [l.text() for l in cell.findChildren(QLabel)]
        if labs and labs[0] == label and key in lv._sigma_widgets:
            return lv._edits[key], lv._sigma_widgets[key][1]
    raise StepFailed("the list view has no ± σ box beside %r" % label)


def _engine(method, **fields):
    def act(dlg):
        set_field(dlg, "Method", method)
        for k, v in fields.items():
            set_field(dlg, k, v)
    return act


def _row(label):
    return r"\| %s \| " % re.escape(label)


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_reliability.xlsx")))
        from xslope.water import with_water_loads
        derived = with_water_loads(st.data).get("dloads_derived") or []
        loads = [float(q["Normal"]) for blk in derived for q in blk]
        w.require(loads, "no water load was derived from the line")
        w.expect("deepest water load", max(loads), r"water is 40 ft deep and presses (\d+) psf")
        crest = [float(q["Normal"]) for blk in derived for q in blk
                 if abs(float(q["Y"]) - 30.0) < 1e-9]
        w.expect("water load over the crest", min(crest) if crest else None,
                 r"10 ft of water presses (\d+) psf")

    with w.step("The standard deviations"):
        seen = {}

        def act(dlg):
            show_parameters_for(dlg, {"LEM", "Reliability"})
            press(dlg, "List view")
            for lab in ("γ", "c (psf)"):
                v, s = _sigma_edit(dlg, lab)
                seen[lab] = (v.text(), s.text())
            return "reject"
        st.edit("Materials", act)
        g = p.find(r"The clay reads \*\*γ\*\* = `(\d+)` \*\*± σ\*\* `(\d+)` and \*\*c\*\* = `(\d+)` \*\*± σ\*\* `(\d+)`")
        w.check([float(x) for x in seen["γ"]] == [float(g.group(1)), float(g.group(2))],
                "γ reads %s ± %s" % seen["γ"])
        w.check([float(x) for x in seen["c (psf)"]] == [float(g.group(3)), float(g.group(4))],
                "c reads %s ± %s" % seen["c (psf)"])

    with w.step("The deterministic run"):
        start = dict(st.data["circles"][0])
        one = w.lem(dict(st.data), "spencer", analysis="single_surface")
        w.expect("the file's circle alone", one["results"]["FS"], r"solved on its own it gives (1\.\d+), so")
        st.run_lem(lambda dlg: (set_field(dlg, "Method", "Spencer"), set_field(dlg, "Analysis", "Auto search")))
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("deterministic FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, on a circle centered at")
        w.expect("center x", crit["Xo"], r"on a circle centered at \((\d+\.\d+), \d+\.\d+\), still")
        w.expect("center y", crit["Yo"], r"on a circle centered at \(\d+\.\d+, (\d+\.\d+)\), still")
        w.expect("radius before", start.get("R") or start["Yo"] - start["Depth"], r"radius (\d+) ft to \d+\.\d+")
        w.expect("radius after", crit["Yo"] - crit["Depth"], r"radius \d+ ft to (\d+\.\d+),")
        w.expect("surface before", float(one["slice_df"]["dl"].sum()), r"failure surface (\d+\.\d+) ft to")
        w.expect("surface after", float(b["slice_df"]["dl"].sum()), r"failure surface \d+\.\d+ ft to (\d+\.\d+),")
        w.expect("sliding mass", float(b["slice_df"]["w"].sum()), r"carrying ([\d,]+) lb/ft of clay")

    with w.step("The Taylor series"):
        seen = {}

        def act(dlg):
            seen["method"] = value(field(dlg, "Method"))
            seen["lem"] = value(field(dlg, "LEM method"))
            seen["slices"] = value(field(dlg, "Number of slices"))
            seen["search"] = button(dlg, "Search for the critical surface at each solve").isChecked()
            seen["sigmas"] = [" ".join(l.text().split()) for l in labels(dlg) if "±" in l.text()]
        st.run_reliability(act)
        w.check(seen["method"] == "Taylor series (TSPM)" and seen["lem"] == "Spencer" and
                int(seen["slices"]) == 40, "the dialog opens on %s / %s / %s" %
                (seen["method"], seen["lem"], seen["slices"]))
        w.check(seen["search"], "the search box is not ticked")
        for q in re.findall(r"`(mat:soil:\w+ = \d+ ± \d+ \(COV \d+%\))`", p.text):
            w.check(any(q in s for s in seen["sigmas"]), "the σ list does not show %r" % q)
        r = _rel(st)
        tbl = p.table(r"The per-parameter table is written to the Log pane:")
        by = {"Mat %s %s" % (q["material_id"], q["param"]): q for q in r["param_info"]}
        for row in tbl[1:]:
            q = by.get(row[0])
            if q is None:
                w.fail("the Taylor series has no row %r (it has %s)" % (row[0], sorted(by)))
                continue
            for j, key in ((5, "F_plus"), (6, "F_minus"), (7, "delta_F")):
                w.expect("%s %s" % (row[0], key), q[key], _row(row[0]) + r"[^|]+ \| " * (j - 1) + r"(\d\.\d+) \|")
        for label, key, pat, scale in (
                ("σ_F", "sigma_F", r"\*\*σ<sub>F</sub> = (\d\.\d+)\*\*", 1),
                ("COV_F", "COV_F", r"\*\*COV<sub>F</sub> = (\d\.\d+)\*\*", 1),
                ("β_LN", "beta_ln", r"\*\*β<sub>LN</sub> = (\d\.\d+)\*\*, a", 1),
                ("reliability (%)", "reliability", r"\*\*reliability of (\d+\.\d+)%\*\*", 100),
                ("P_f (%)", "prob_failure", r"\*\*probability of failure of (\d+\.\d+)%\*\*", 100)):
            w.expect("Taylor %s" % label, r[key], pat, scale=scale)
        w.expect("Taylor β_N", r.get("beta_normal"),
                 r"β<sub>N</sub> = \(F<sub>MLV</sub> − 1\)/σ<sub>F</sub> = \*\*(\d\.\d+)\*\*")
        w.expect("Taylor P_f, normal (%)", r.get("pf_normal"),
                 r"gives a probability of failure of \*\*(\d+\.\d+)%\*\*; the run summary", scale=100)
        taylor = r

    def mc(sampling=None):
        fields = {}
        if sampling is not None:
            fields["MC sampling"] = sampling

        def act(dlg):
            set_field(dlg, "Method", "Monte Carlo")
            seen["mc"] = {k: value(field(dlg, k)) for k in ("MC samples", "MC seed", "MC distribution",
                                                             "MC sampling")}
            from PySide6.QtWidgets import QCheckBox
            seen["mc_search"] = [b_.text() for b_ in dlg.findChildren(QCheckBox)
                                 if b_.text().startswith("Search for")]
            seen["converge"] = button(dlg, "Stop when P_f converges").isChecked()
            for k, v in fields.items():
                set_field(dlg, k, v)
        return act

    with w.step("Monte Carlo"):
        seen = {}
        st.run_reliability(mc())
        want = {"MC samples": p.number(r"\*\*MC samples\*\* `(\d+)`")[0],
                "MC seed": p.number(r"\*\*MC seed\*\* `(\d+)`")[0],
                "MC distribution": p.find(r"\*\*MC distribution\*\* `(\w+)`").group(1),
                "MC sampling": p.find(r"\*\*MC sampling\*\* `([\w ]+)`").group(1)}
        for k, v in want.items():
            got = str(seen["mc"][k]).replace(",", "")
            w.check(got == v or (_f(got) is not None and _f(got) == _f(v)),
                    "%s opens on %r, not %r" % (k, seen["mc"][k], v))
        w.check(not seen["converge"], "Stop when P_f converges is ticked")
        w.check("Search for the critical surface at the mean values" in seen["mc_search"],
                "the search box does not read 'at the mean values': %s" % seen["mc_search"])
        r = _rel(st)
        n = int(r.get("n_used") or float(want["MC samples"]))
        pf = float(r.get("pf_empirical", r["prob_failure"]))
        w.expect("MC mean", r["mean_FS"], r"\*\*Mean FS = (1\.\d+), σ")
        w.expect("MC σ_F", r["sigma_F"], r"σ<sub>F</sub> = (0\.\d+)\*\*, and the shaded")
        w.expect("MC failures", round(pf * n), r"holds \*\*([\d,]+) of the [\d,]+ realizations")
        w.expect("MC P_f (%)", pf, r"P<sub>f</sub> = (\d+\.\d+)%\*\*, counted directly", scale=100)
        w.expect("MC β normal", r["beta_normal"], r"= \*\*(\d\.\d+)\*\*, and the lognormal")
        w.expect("MC β LN", r["beta_ln"], r"β<sub>LN</sub> = \*\*(\d\.\d+)\*\* computed from")
        by = {q["param"]: q for q in taylor["param_info"]}
        fm = taylor["F_MLV"]
        g = by.get("gamma")
        c = by.get("c")
        if g and c:
            w.expect("γ pair above F_MLV", (g["F_plus"] + g["F_minus"]) / 2 - fm,
                     r"The unit weight's pair averages \*\*\+(\d\.\d+) above\*\*")
            w.check(abs((c["F_plus"] + c["F_minus"]) / 2 - fm) < 5e-8,
                    "the cohesion's pair does not average to F_MLV to seven digits")
        w.expect("MC mean above F_MLV", r["mean_FS"] - fm, r"the Monte Carlo mean sits \*\*\+(\d\.\d+) above\*\*")
        mc_lhs = r

    with w.step("Trying random sampling"):
        st.run_reliability(mc("Random"))
        # The dialog reopens on the last run's settings; the runs after this one
        # are on Latin hypercube only where the page says to set it back.
        back = p.finditer(r"set \*\*MC sampling\*\* back to `(Latin hypercube)`")
        after = {"MC sampling": back[0].group(1)} if back else {}
        r = _rel(st)
        w.expect("random P_f (%)", float(r.get("pf_empirical", r["prob_failure"])),
                 r"The result is \*\*P<sub>f</sub> = (\d+\.\d+)%\*\*", scale=100)
        mc_rand = r

    with w.step("The response surface"):
        st.run_reliability(_engine("Response surface (RS)", **after))
        r = _rel(st)
        w.expect("RS P_f (%)", float(r.get("pf_empirical", r["prob_failure"])),
                 r"\*\*P<sub>f</sub> = (\d+\.\d)%\*\* and", scale=100)
        w.expect("RS β LN", r["beta_ln"], r"\*\*β<sub>LN</sub> = (\d\.\d+)\*\*, beside Monte")
        w.expect("design solves", r.get("n_design"), r"The engine solves the real model (\d+) times to fit")
        w.expect("gate solves", r.get("n_gate"), r"to check the fit: (\d+) solves spread")
        w.expect("failure-region solves", r.get("tail_n"), r"and (\d+) concentrated where the")
        w.expect("typical error", r.get("gate_rms"), r"with a typical error of (\d\.\d+) \(R²")
        w.expect("R²", r.get("gate_r2"), r"\(R² = (\d\.\d+)\)")
        w.check(int(r.get("gate_pf_disagree_n") or 0) == 0,
                "%s of the checking solves landed across FS = 1" % r.get("gate_pf_disagree_n"))
        rs = r

    with w.step("Reading β against P_f"):
        rows_ = (("Taylor series (TSPM)", taylor["F_MLV"], taylor["sigma_F"], taylor["beta_ln"],
                  taylor["prob_failure"]),
                 ("Monte Carlo, 10,000 samples (Latin hypercube)", mc_lhs["mean_FS"], mc_lhs["sigma_F"],
                  mc_lhs["beta_ln"], mc_lhs.get("pf_empirical", mc_lhs["prob_failure"])),
                 ("Monte Carlo, 10,000 samples (random)", mc_rand["mean_FS"], mc_rand["sigma_F"],
                  mc_rand["beta_ln"], mc_rand.get("pf_empirical", mc_rand["prob_failure"])),
                 ("Response surface, 10⁷ formula realizations", rs["mean_FS"], rs["sigma_F"],
                  rs["beta_ln"], rs.get("pf_empirical", rs["prob_failure"])))
        for label, f, s, beta, pf in rows_:
            at = _row(label)
            w.expect("%s F" % label, f, at + r"(\d\.\d+) \|")
            w.expect("%s σ_F" % label, s, at + r"[\d.]+ \| (\d\.\d+) \|")
            w.expect("%s β" % label, beta, at + r"[\d.]+ \| [\d.]+ \| (\d\.\d+) ")
            w.expect("%s P_f" % label, float(pf), at + r"[\d.]+ \| [\d.]+ \| [^|]+ \| (\d+\.\d+)%", scale=100)

    with w.step("Which uncertainty to reduce"):
        terms = {q["param"]: (q["delta_F"] / 2) ** 2 for q in taylor["param_info"]}
        total = sum(terms.values())
        w.expect("cohesion's share (%)", terms.get("c", 0) / total,
                 r"the cohesion accounts for \*\*(\d+\.\d)%\*\* of the variance", scale=100)
        w.expect("unit weight's share (%)", terms.get("gamma", 0) / total,
                 r"the unit weight \*\*(\d+\.\d)%\*\*\.", scale=100)
        st.run_parametric(lambda dlg: set_field(dlg, "Plot type", "Variance Pareto (σ)"))
        mcn = p.number(r"set \*\*Plot type\*\* to `Monte Carlo rank correlation \(σ\)` and \*\*MC samples\*\* to `(\d+)`")[0]
        st.run_parametric(lambda dlg: (set_field(dlg, "Plot type", "Monte Carlo rank correlation (σ)"),
                                       set_field(dlg, "MC samples", mcn)))

    with w.step("Halving σ(c)"):
        page_sigma = p.table(r"the first two cells then read:")[1]

        def act(dlg):
            show_parameters_for(dlg, {"LEM", "Reliability"})
            press(dlg, "List view")
            _v, s = _sigma_edit(dlg, "c (psf)")
            set_value(s, page_sigma[1])
        st.edit("Materials", act)
        m = st.data["materials"][0]
        w.require(float(m.get("sigma_c") or 0) == float(page_sigma[1]) and
                  float(m.get("sigma_gamma") or 0) == float(page_sigma[0]) and float(m["c"]) == 400,
                  "σ(c) / σ(γ) / c saved as %s / %s / %s" % (m.get("sigma_c"), m.get("sigma_gamma"), m["c"]))
        st.run_reliability(_engine("Taylor series (TSPM)"))
        t50 = _rel(st)
        st.run_reliability(_engine("Monte Carlo", **after))
        m50 = _rel(st)
        for label, f, s, beta, pf in (("TSPM, σ(c) = 50", t50["F_MLV"], t50["sigma_F"], t50["beta_ln"],
                                       t50["prob_failure"]),
                                      ("Monte Carlo, σ(c) = 50", m50["mean_FS"], m50["sigma_F"], m50["beta_ln"],
                                       m50.get("pf_empirical", m50["prob_failure"]))):
            at = _row(label)
            w.expect("%s F" % label, f, at + r"(\d\.\d+) \|")
            w.expect("%s σ_F" % label, s, at + r"[\d.]+ \| (\d\.\d+) \|")
            w.expect("%s β" % label, beta, at + r"[\d.]+ \| [\d.]+ \| (\d\.\d+) \|")
            w.expect("%s P_f" % label, float(pf), at + r"[\d.]+ \| [\d.]+ \| [\d.]+ \| (\d+\.\d+)%", scale=100)
        n = int(m50.get("n_used") or 10000)
        w.expect("failed count at σ(c) = 50", round(float(m50.get("pf_empirical", m50["prob_failure"])) * n),
                 r"fell from [\d,]+ realizations to (\d+)\.")
        terms = {q["param"]: (q["delta_F"] / 2) ** 2 for q in t50["param_info"]}
        w.expect("unit weight's share after (%)", terms.get("gamma", 0) / sum(terms.values()),
                 r"it now makes up (\d+)% of", scale=100)


def _f(s):
    try:
        return float(str(s).replace(",", ""))
    except (TypeError, ValueError):
        return None
