"""LEM-14, Pseudostatic Analysis: all three parts in Studio, from the starter.

Walked, on copies of the linked workbooks: k_h typed into Global parameters;
the material read against the page's row; the one-stage search; the page's
two-stage row pasted into Materials and the two-stage runs (K_c = 1 envelope,
interpolated, strength factor 0.8) with the Log lines the page quotes; the water
table pasted into Piezometric lines after the page's Add row count; u set to
piezo and the wet runs; the R-envelope row pasted back and the wet one-stage
search; the strength reduction on the completed wet model (mesh and run at the
file's declared settings, the k_h direction note). Not walked: the comparison
paragraphs' readings of Duncan et al. (they quote the book, not a run).
"""
import re

from tutorial_walks.harness import (column, editable, field, page_rows, paste,
                                    preflight_items, press, set_field, set_value,
                                    show_parameters_for, tick, value)

PAGE = "lem14_pseudostatic.md"
TIMEOUT = 2 * 3600


def _paste_row(st, row):
    def act(dlg):
        press(dlg, "Table view")
        show_parameters_for(dlg, {"LEM"})
        t = editable(dlg)
        paste(t, [row], 0, column(t, "name"))
    st.edit("Materials", act)


def _run(st, two_stage=False, strength=None, factor=None):
    def act(dlg):
        set_field(dlg, "Method", "Spencer")
        set_field(dlg, "Analysis", "Auto search")
        tick(dlg, "Two-stage seismic", two_stage)
        if strength is not None:
            set_field(dlg, "Stage 2 strength", strength)
        if factor is not None:
            set_field(dlg, "Strength factor", factor)
    return st.run_lem(act)


def walk(w):
    p = w.page
    st = w.studio()
    rows = [page_rows(t)[0] for _line, t in p.tables() if t and t[0][:2] == ["name", "γ (pcf)"]]
    w.require(len(rows) == 3, "the page has %d material rows, not three" % len(rows))
    r_dry, two_stage, r_wet = rows

    with w.step("Opening the model"):
        st.open(w.copy_model(w.page_link("xslope_pseudostatic_start.xlsx")))

    with w.step("Entering the seismic coefficient"):
        kh = p.number(r"set \*\*Horizontal seismic coefficient k_h\*\* to `(\d\.\d+)`")[0]

        def act(dlg):
            set_field(dlg, "Horizontal seismic coefficient k_h", kh)
            w.check(float(value(field(dlg, "Vertical seismic coefficient k_v (+ down)")) or 0) == 0.0,
                    "k_v does not read 0")
        st.edit("Global parameters", act)
        w.require(float(st.data.get("k_seismic") or 0) == float(kh), "k_h saved as %s" % st.data.get("k_seismic"))

    with w.step("One stage with the R envelope"):
        m = st.data["materials"][0]
        got = [m["name"], m["gamma"], m.get("gamma_sat"), m["option"], m["c"], m["phi"], m.get("cp"),
               m.get("r_elev"), m.get("d"), m.get("psi"), m["u"]]
        for want, have in zip(r_dry, got):
            w.check(want == str(have) or (want == "" and have in (None, "", 0)) or
                    (_f(want) is not None and _f(want) == _f(have)),
                    "the starter's material reads %s, not the page's %s" % (got, r_dry))
        _run(st)
        w.expect("dry one-stage FS", st.results["lem_solution"]["results"]["FS"],
                 r"to the toe: .*?\*\*FS = (0\.\d+)\.\*\*")

    blocks = {}
    for m_ in re.finditer(r"```text\n(Critical FS = [\d.]+\n[^`]+?)```", p.text):
        lines = m_.group(1).strip().splitlines()
        blocks[lines[0]] = lines

    def check_block(log, fs_line, label):
        if fs_line not in blocks:
            w.fail("the page quotes no Log block starting %r" % fs_line)
            return
        w.check_log(blocks[fs_line], log, label)

    with w.step("Two stages"):
        _paste_row(st, two_stage)
        m = st.data["materials"][0]
        w.require((float(m["c"]), float(m["phi"]), float(m["d"]), float(m["psi"])) ==
                  tuple(float(x) for x in (two_stage[4], two_stage[5], two_stage[8], two_stage[9])),
                  "the pasted row saved as c %s φ %s d %s ψ %s" % (m["c"], m["phi"], m["d"], m["psi"]))
        _s, log = _run(st, True, "Kc = 1 envelope", "1")
        check_block(log, "Critical FS = 1.068", "two stages, Kc = 1")
        _s, log = _run(st, True, "Interpolated (as in rapid drawdown)", "1")
        check_block(log, "Critical FS = 1.282", "two stages, interpolated")

    with w.step("Reducing the undrained strength"):
        f = p.number(r"set \*\*Strength factor\*\* to `(\d\.\d+)`")[0]
        _s, log = _run(st, True, "Kc = 1 envelope", f)
        check_block(log, "Critical FS = 0.855", "strength factor %s" % f)

    with w.step("Entering the water table"):
        n = p.count(r"leave \*\*Type\*\* at piezo, add (\d+) rows")
        pts = page_rows(p.table(r"add \d+ rows, and enter \(or paste\) the points"))

        def act(dlg):
            for _ in range(n):
                press(dlg, "Add row")
            t = editable(dlg)
            w.require(t.table.rowCount() == len(pts), "%d Add row clicks give %d rows for %d points"
                      % (n, t.table.rowCount(), len(pts)))
            paste(t, pts)
        st.edit("Piezometric lines", act)
        got = [tuple(map(float, xy)) for xy in st.data.get("piezo_line") or []]
        w.require(got == [tuple(map(float, r)) for r in pts], "the water table saved as %s" % got[:3])

        def u_piezo(dlg):
            from PySide6.QtWidgets import QComboBox
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM"})
            t = editable(dlg)
            cw = t.table.cellWidget(0, column(t, "u"))
            set_value(cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox), "piezo")
        st.edit("Materials", u_piezo)
        w.require(st.data["materials"][0]["u"] == "piezo", "u is not piezo")

    with w.step("Two stages with the water table"):
        _s, log = _run(st, True, "Kc = 1 envelope", "1")
        check_block(log, "Critical FS = 0.934", "wet, Kc = 1")
        _s, log = _run(st, True, "Interpolated (as in rapid drawdown)", "1")
        check_block(log, "Critical FS = 1.161", "wet, interpolated")

    with w.step("One stage with the water table"):
        _paste_row(st, r_wet)
        w.require(st.data["materials"][0]["u"] == "piezo" and float(st.data["materials"][0]["d"] or 0) == 0.0,
                  "the R-envelope row did not paste back")
        _run(st, False)
        b = st.results["lem_solution"]
        w.expect("wet one-stage FS", b["results"]["FS"], r"the pore pressure on each slice: .*?\*\*FS = (0\.\d+)\.\*\*")
        xs = [x for x, _y in b["failure_surface"].coords]
        w.check(min(xs) < 0.0 or max(xs) > 30.0 + 1e-6,
                "the wet one-stage circle does not come out on the level ground in front of the toe")
        done = w.load(w.copy_model(w.page_link("xslope_pseudostatic.xlsx"), "done"))
        for k in ("k_seismic", "piezo_line"):
            w.check(_same(st.data.get(k), done.get(k)), "the walked model's %s differs from the completed file's" % k)

    with w.step("Strength reduction at the same coefficient"):
        st.open(w.copy_model(w.page_link("xslope_pseudostatic.xlsx"), "fem"))
        st.set_mode("FEM")

        def mesh(dlg):
            w.check(value(field(dlg, "Element type")) == "Quadratic triangles (tri6)", "Element type is not tri6")
            w.check(abs(float(value(field(dlg, "Target element size"))) - 1.5) < 1e-9,
                    "Target element size does not read 1.500")
        st.build_mesh(mesh)
        seen = {}

        def run(dlg):
            seen["fmin"] = value(field(dlg, "F min (SSRM)"))
            seen["fmax"] = value(field(dlg, "F max (SSRM)"))
            seen["items"] = preflight_items(dlg) or []
        st.run_fem(run, timeout=2 * 3600)
        w.check((float(seen["fmin"]), float(seen["fmax"])) == (0.5, 1.5),
                "the bracket opens on %s / %s, not 0.50 / 1.50" % (seen["fmin"], seen["fmax"]))
        note = p.find(r"shows one note, \*(Horizontal seismic coefficient k_h = [\d.]+ pushes in \+x)\*").group(1)
        w.check(any(note in t for t in seen["items"]), "Model checks does not show %r: %s" % (note, seen["items"]))
        r = (st.results.get("fem_solution") or {}).get("ssrm_result")
        w.require(r is not None, "the run left no strength-reduction result")
        w.expect("SSRM FS", r["FS"], r"\*\*FS = (0\.\d+)\.\*\* The shear strain at failure")


def _f(v):
    try:
        return float(str(v).replace("−", "-"))
    except (TypeError, ValueError):
        return None


def _same(a, b):
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    fa, fb = _f(a), _f(b)
    return a == b or (fa is not None and fb is not None and abs(fa - fb) < 1e-9)
