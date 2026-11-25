"""LEM-13, A Rock Slope (Hoek-Brown): Part A and Part B in Studio.

Walked, on a copy of the linked workbook: the Hoek-Brown material in List view
on LEM alone (the four observations, and σci typed as 30 drawing the units
warning before it is put back); the Spencer search and Bishop; the critical
surface's slice table read for the stresses and tangents the page quotes; the
mesh the file declares and the strength reduction at the dialog's defaults; the
Parametric design sweeps over GSI and over D with the rock-mass constants each
step solved with; the Corps of Engineers pairing note in Run LEM's checks and
the methods that converge. Not walked: the envelope figures (drawn by the figure
producer; their tangent at the mean stress is checked here through the engine's
own Hoek-Brown tangent), and the RS2-60 normalized models the last paragraph
mentions.
"""
import re

import numpy as np

from tutorial_walks.harness import (StepFailed, button, field, preflight_items, press,
                                    set_field, show_parameters_for, value)

PAGE = "lem13_rock_slope.md"
TIMEOUT = 2 * 3600


def _sigma_n(df):
    return (df["n_eff"] / df["dl"]).to_numpy(dtype=float)


def walk(w):
    p = w.page
    st = w.studio()

    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_rock_slope.xlsx")))

    with w.step("The Hoek-Brown material"):
        seen = {}
        obs = dict(zip(("σci", "GSI", "mi", "D"),
                       p.find(r"\*\*σci\*\* = (\d+), \*\*GSI\*\* = (\d+), \*\*mi\*\* = (\d+), \*\*D\*\* = (\d+)").groups()))

        def act(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "List view")
            seen["option"] = value(field(dlg, "Model (option)"))
            for k in obs:
                seen[k] = value(field(dlg, k))
            return "reject"
        st.edit("Materials", act)
        w.check(seen["option"] == "hb", "Model (option) reads %r, not hb" % seen["option"])
        for k, v in obs.items():
            w.check(float(str(seen[k]).replace(",", "")) == float(v), "%s reads %r, not %s" % (k, seen[k], v))
        from xslope.hoekbrown import hb_constants
        mb, s, a = hb_constants(float(obs["GSI"]), float(obs["mi"]), float(obs["D"]))
        w.expect("m_b", mb, r"m<sub>b</sub> = (\d\.\d+), \*s\* = ")
        w.expect("s", s, r"\*s\* = (\d\.\d+ × 10<sup>−5</sup>), \*a\* = ")
        w.expect("a", a, r"\*a\* = (\d\.\d+) — reproducing")
        w.expect("unconfined mass strength", float(obs["σci"]) * s ** a,
                 r"σ<sub>ci</sub>·s<sup>a</sup>, is (\d+\.\d) kPa at this GSI")

    with w.step("σci typed as 30: the units warning"):
        def low(dlg):
            show_parameters_for(dlg, {"LEM"})
            press(dlg, "List view")
            set_field(dlg, "σci", "30")
        st.edit("Materials", low)
        seen = {}

        def look(dlg):
            seen["items"] = preflight_items(dlg) or []
            return "reject"
        st.trigger("Run", "Run LEM", look)
        w.check(any("σci" in t and "1000" in t for t in seen["items"]),
                "typed as 30, σci draws no units warning: %s" % seen["items"][:4])
        st.trigger("Edit", "Undo")
        w.require(float(st.data["materials"][0].get("hb_sci") or 0) == float(obs["σci"]),
                  "Undo did not put σci back to %s" % obs["σci"])

    with w.step("Running the analysis"):
        seen = {}

        def act(dlg):
            set_field(dlg, "Method", "Spencer")
            set_field(dlg, "Analysis", "Auto search")
            seen["items"] = preflight_items(dlg) or []
        st.run_lem(act)
        w.check(any("No problems found for this run" in t for t in seen["items"]),
                "Model checks: %s" % seen["items"][:3])
        b = st.results["lem_solution"]
        crit = b["search"]["fs_cache"][0]
        w.expect("Spencer FS", b["results"]["FS"], r"\*\*FS = (1\.\d+)\*\*, against Hammah's published")
        w.expect("center x", crit["Xo"], r"on a circle centered at \((\d+\.\d+), \d+\.\d+\) tangent")
        w.expect("center y", crit["Yo"], r"on a circle centered at \(\d+\.\d+, (\d+\.\d+)\) tangent")
        w.expect("tangent elevation", crit["Depth"], r"tangent at elevation (\d+\.\d+), [\d.]+ m of failure")
        df = b["slice_df"]
        w.expect("surface length", float(df["dl"].sum()), r"elevation \d+\.\d+, (\d+\.\d+) m of failure surface")
        w.expect("mass", float(df["w"].sum()), r"surface carrying (\d+) kN/m of rock")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.check(abs(min(xs) - 20.0) < 0.01, "the surface does not exit at the toe (x = %.3f)" % min(xs))
        st.run_lem(lambda dlg: set_field(dlg, "Method", "Bishop's Simplified"))
        w.expect("Bishop FS", st.results["lem_solution"]["results"]["FS"],
                 r"Bishop's simplified method on the same file gives (1\.\d+) against")

    with w.step("Strength along the critical surface"):
        sig = _sigma_n(df)
        c = df["c"].to_numpy(dtype=float)
        phi = df["phi"].to_numpy(dtype=float)
        dl = df["dl"].to_numpy(dtype=float)
        w.expect("largest normal stress", sig.max(), r"effective normal stress reaches (\d+\.\d) kPa")
        w.expect("smallest normal stress", sig.min(), r"falls to (−\d+\.\d) kPa at the crest")
        w.expect("smallest c", c.min(), r"runs from c<sub>i</sub> = (\d+\.\d) to")
        w.expect("largest c", c.max(), r"c<sub>i</sub> = \d+\.\d to (\d+\.\d) kPa")
        w.expect("smallest φ", phi.min(), r"φ<sub>i</sub> from (\d+\.\d)° to")
        w.expect("largest φ", phi.max(), r"φ<sub>i</sub> from \d+\.\d° to (\d+\.\d)°")
        w.expect("weighted c", float((c * dl).sum() / dl.sum()), r"averages c<sub>i</sub> = (\d+\.\d) kPa")
        w.expect("weighted φ", float((phi * dl).sum() / dl.sum()), r"and φ<sub>i</sub> = (\d+\.\d)°, at a mean")
        mean_sig = float((sig * dl).sum() / dl.sum())
        w.expect("mean normal stress", mean_sig, r"at a mean normal stress of (\d+\.\d) kPa")
        from xslope.hoekbrown import hb_tangent
        ci, phii = hb_tangent(mean_sig, float(obs["σci"]), float(obs["GSI"]), float(obs["mi"]), float(obs["D"]))
        w.expect("tangent cohesion at the mean stress", float(ci), r"cohesion c<sub>i</sub> = (\d+\.\d) kPa \(the line's intercept\)")
        w.expect("tangent friction angle at the mean stress", float(phii), r"friction angle φ<sub>i</sub> = (\d+\.\d)° \(its slope\)")

    with w.step("The same file through strength reduction"):
        st.set_mode("FEM")
        w.check(not st.menu_action("Run", "Run FEM").isEnabled(), "Run FEM… is enabled before a mesh")
        size = p.number(r"\*\*Target element size\*\* reads `(\d\.\d+)`")[0]

        def mesh(dlg):
            w.check(value(field(dlg, "Element type")) == "Quadratic triangles (tri6)", "Element type is not tri6")
            w.check(not button(dlg, "Auto-size from geometry").isChecked(), "Auto-size is ticked")
            w.check(abs(float(value(field(dlg, "Target element size"))) - float(size)) < 1e-9,
                    "Target element size does not read %s" % size)
        st.build_mesh(mesh)
        m = st.data["mesh"]
        w.expect("mesh nodes", len(m["nodes"]), r"The mesh comes out at \*\*([\d,]+) nodes and")
        w.expect("mesh triangles", len(m["elements"]), r"nodes and ([\d,]+) triangles\*\*:")
        seen = {}

        def run(dlg):
            seen.update({k: value(field(dlg, k)) for k in
                         ("Analysis", "F min (SSRM)", "F max (SSRM)", "Tolerance (SSRM)",
                          "Max iterations per trial", "Iteration ceiling", "Failure criterion")})
        st.run_fem(run, timeout=2 * 3600)
        w.check(str(seen["Analysis"]).startswith("SSRM"), "Analysis opens on %s" % seen["Analysis"])
        for k, want in (("F min (SSRM)", 0.8), ("F max (SSRM)", 1.6), ("Tolerance (SSRM)", 0.01),
                        ("Max iterations per trial", 12000), ("Iteration ceiling", 50000)):
            w.check(abs(float(str(seen[k]).replace(",", "")) - want) < 1e-9, "%s opens on %s" % (k, seen[k]))
        w.check(seen["Failure criterion"] == "Non-convergence", "Failure criterion opens on %s" % seen["Failure criterion"])
        r = (st.results.get("fem_solution") or {}).get("ssrm_result")
        w.require(r is not None, "the run left no strength-reduction result")
        w.expect("SSRM FS", r["FS"], r"\*\*FS = (1\.\d+)\*\*\. Hammah et al\. report")
        st.set_mode("LEM")

    def design(prop, lo, hi, steps, target=None):
        def act(dlg):
            set_field(dlg, "Mode", "Design (FS target)")
            w.check(value(field(dlg, "Method")) == "Spencer", "Method is not Spencer")
            w.check(int(value(field(dlg, "Number of slices"))) == 40, "slices are not 40")
            set_field(dlg, "Material", "rock")
            set_field(dlg, "Property", prop)
            w.check(value(field(dlg, "Sweeping")) == "mat:rock:%s" % prop,
                    "the Sweeping row reads %r" % value(field(dlg, "Sweeping")))
            set_field(dlg, "From", lo)
            set_field(dlg, "To", hi)
            set_field(dlg, "Steps", steps)
            if target is not None:
                set_field(dlg, "Target FS", target)
            w.check(button(dlg, "Re-search the critical surface at each step").isChecked(),
                    "Re-search the critical surface at each step is not ticked")
        st.run_parametric(act, timeout=2 * 3600)
        bnd = st.results.get("design")
        if bnd is None:
            raise StepFailed("the sweep left no design result")
        df = bnd["df"]
        sw = df.loc[~df["is_base"] & df["success"]].sort_values("value")
        return bnd, list(zip(sw["value"].astype(float), sw["fs"].astype(float)))

    from xslope.hoekbrown import hb_constants

    with w.step("Sweeping the Geological Strength Index"):
        lo, hi, n, tgt = p.find(r"set \*\*From\*\* `(\d+)`, \*\*To\*\* `(\d+)`, \*\*Steps\*\* `(\d+)` and "
                                r"\*\*Target FS\*\* `(\d\.\d)`").groups()
        bnd, pts = design("hb_gsi", lo, hi, n, tgt)
        w.expect("GSI for FS = 1.5", bnd.get("crossing"), r"\*\*GSI = (\d+\.\d+) for FS = 1\.5\*\*")
        tbl = p.table(r"Each step solved with its own rock-mass constants:")
        for (gsi, fs), row in zip(pts, tbl[1:]):
            at = r"\| %s \| " % row[0]
            mb, s, a = hb_constants(gsi, float(obs["mi"]), 0.0)
            w.check(abs(gsi - float(row[0])) < 1e-9, "the sweep solved GSI %g, the table lists %s" % (gsi, row[0]))
            w.expect("m_b at GSI %s" % row[0], float(mb), at + r"(\d\.\d+) \|")
            w.expect("s at GSI %s" % row[0], float(s), at + r"\d\.\d+ \| ([^|]+) \|")
            w.expect("a at GSI %s" % row[0], float(a), at + r"\d\.\d+ \| [^|]+ \| (\d\.\d+) \|")
            w.expect("σci·s^a at GSI %s" % row[0], float(obs["σci"]) * float(s) ** float(a),
                     at + r"\d\.\d+ \| [^|]+ \| \d\.\d+ \| ([\d.]+) \|")
            w.expect("FS at GSI %s" % row[0], fs, at + r"\d\.\d+ \| [^|]+ \| \d\.\d+ \| [\d.]+ \| (\d\.\d+) \|")

    with w.step("The disturbance factor"):
        lo, hi, n = p.find(r"\*\*Property\*\* to `hb_d`, \*\*From\*\* `(\d+)`, \*\*To\*\* `(\d+)`, \*\*Steps\*\* `(\d+)`").groups()
        bnd, pts = design("hb_d", lo, hi, n)
        w.check(not bnd.get("bracketed"), "the D sweep brackets the target")
        spans = p.find(r"the factor of safety spans \[(\d\.\d+), (\d\.\d+)\]").groups()
        w.check("spans [%s, %s]" % spans in (bnd.get("message") or ""),
                "the banner does not read spans [%s, %s]: %r" % (spans + (bnd.get("message"),)))
        tbl = p.table(r"The curve crosses FS = 1 instead")
        for (d, fs), row in zip(pts, tbl[1:]):
            at = r"\| %s \| " % re.escape(row[0])
            mb, s, a = hb_constants(float(obs["GSI"]), float(obs["mi"]), d)
            w.expect("m_b at D = %s" % row[0], float(mb), at + r"(\d\.\d+) \|")
            w.expect("s at D = %s" % row[0], float(s), at + r"\d\.\d+ \| ([^|]+) \|")
            w.expect("σci·s^a at D = %s" % row[0], float(obs["σci"]) * float(s) ** float(a),
                     at + r"\d\.\d+ \| [^|]+ \| ([\d.]+) \|")
            w.expect("FS at D = %s" % row[0], fs, at + r"\d\.\d+ \| [^|]+ \| [\d.]+ \| (\d\.\d+) \|")

    with w.step("Steep friction angles and method choice"):
        from xslope.hoekbrown import hb_tangent
        sci, mi = float(obs["σci"]), float(obs["mi"])
        phi = lambda g, m: float(hb_tangent(0.0, sci, g, m, 0.0)[1])
        w.expect("φ at zero confinement, GSI 5", phi(5.0, mi), r"zero confinement is (\d+\.\d)° at GSI = 5")
        w.expect("φ at zero confinement, GSI 70", phi(70.0, mi), r"and (\d+\.\d)° at GSI = 70")
        peak = max(phi(g, mi) for g in np.arange(5.0, 70.01, 0.25))
        w.expect("peak φ between", peak, r"peaking at only (\d+\.\d)° in between")
        w.expect("m_i 25, GSI 5", phi(5.0, 25.0), r"the same two ends give (\d+\.\d)° and")
        w.expect("m_i 25, GSI 70", phi(70.0, 25.0), r"give \d+\.\d° and (\d+\.\d)°\.")
        for meth, warns in (("Corps of Engineers", True), ("Lowe & Karafiath", True), ("Morgenstern-Price", False)):
            seen = {}

            def look(dlg, m=meth):
                set_field(dlg, "Method", m)
                seen["items"] = preflight_items(dlg) or []
                return "reject"
            st.trigger("Run", "Run LEM", look)
            said = any("Hoek-Brown" in t for t in seen["items"][1:])
            w.check(said == warns, "%s: the checks column %s the Hoek-Brown pairing: %s"
                    % (meth, "does not name" if warns else "names", seen["items"]))
        sd = w.load(w.copy_model(w.page_link("xslope_rock_slope.xlsx"), "engine"))
        for meth in ("corps", "lowe"):
            bb = w.lem(sd, meth)
            w.check(bb["results"].get("FS") is not None and not bb["search"]["unsolved"]["unsolved"],
                    "%s does not converge on this section" % meth)
