"""Standing checks on what Studio's editors gray, show and save, held to what the
engines read — the Studio half of the input-graying audit.

A grayed cell is a promise that nothing reads it, and a cell shown live is a promise
that something does. A dialog that shows a value the file does not hold, or saves
one the user did not type, changes the model with nothing on screen saying so. Each
check below pins one of those promises to the engine's reading:

  A. MATERIALS. The graying rule (``_mat_dim_keys``) over every strength option,
     pore-pressure model and unsaturated model, blanks included: a blank u reads as
     none and a blank unsat as lf, as the loader reads them; the suction pair is
     live only for u = piezo or seep on an option that takes suction; gard keeps h0
     (transient seepage reads it); s(γ) is gray for elastic, pow and hb. The list
     view shows and grays by the same rule.
  B. REINFORCEMENT. Lp1/Lp2 gray exactly when Adhesion and Delta are both filled;
     Adhesion and Delta are never grayed and nothing in the editor clears them; kn,
     ks and Jred gray unless Joint is Yes; the table and list views agree row for
     row; a new line's development lengths are blank, not 0 (fully anchored).
  C. THE MESH SIGNATURE. An edit to anything the mesher reads (a line's Joint, a
     jointed line's Tend, a load or seepage-BC vertex, a derived water load, an SSR
     zone's Size) marks the mesh stale; a property the mesher does not read does not.
  D. NON-CIRCULAR MOVEMENT. A blank or an unrecognized word is shown as the search
     reads it and saved as it was: never turned into Free on OK.
  E. CIRCLES. The Option grays the size fields it does not read; a blank or
     unknown Option is shown and saved as it is, unsized, and grays nothing; any
     unsized circle (that, or an Option whose cell is blank) keeps its Option and
     what was typed in its size cells until it is sized.
  F. POLYGONS. Size stays live for every Type.
  G. SEEPAGE BOUNDARIES. A value never entered stays blank through the editor, on
     both sets, so the model checks can name it; it never becomes a 0.
  H. TENSION SRF. Run FEM's box is live when a joints-sheet line's cutoff is the
     only one the setting divides.
  I. LINE LOADS. Studio's mesh puts a node at each line load, as the scripts'
     meshes do; moving one makes the mesh stale, changing its magnitude does not.
  J. RELIABILITY. The dialog's model checks judge the run it will make: its
     engine and whether it is a rapid drawdown.

Every Studio setting these dialogs read or write goes to a scratch ini (QSettings
is replaced for the duration), and every model is loaded from a copy.

Skips cleanly when PySide6 is not installed.
"""
import contextlib
import copy
import io
import itertools
import math
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_SIMPLE = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_simple1.xlsx")
_DAM = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_dam.xlsx")
_POLYGONS = os.path.join(_REPO, "docs", "tutorials", "files",
                         "xslope_rock_toppling.xlsx")
_LINER = os.path.join(_REPO, "docs", "tutorials", "files", "xslope_liner_jointed.xlsx")
_ROCK_JOINTS = os.path.join(_REPO, "docs", "tutorials", "files",
                            "xslope_rock_joints.xlsx")


def _fail(cond, message):
    return [] if cond else [message]


def _blank(v):
    return v is None or (isinstance(v, float) and math.isnan(v)) or (
        isinstance(v, str) and not v.strip())


@contextlib.contextmanager
def _sandbox():
    """Scratch settings and scratch model copies for one check.

    Studio's dialogs construct ``QSettings("XSlope", "XSlope Studio")`` with a local
    import, and on macOS that store cannot be redirected (welcome_window_check.py
    measured it), so the class itself is replaced while the check runs: every store
    a dialog opens is one ini in the scratch folder. Yields a function that copies a
    shipped model into the scratch folder and loads the copy."""
    import PySide6.QtCore as QtCore
    from xslope.fileio import load_slope_data

    real = QtCore.QSettings
    tmp = tempfile.mkdtemp(prefix="audit_studio_")
    ini = os.path.join(tmp, "studio.ini")

    class _Scratch(real):
        def __init__(self, *args, **kwargs):
            super().__init__(ini, real.IniFormat)

    def load(path):
        dst = os.path.join(tmp, os.path.basename(path))
        if not os.path.exists(dst):
            shutil.copy(path, dst)
        with contextlib.redirect_stdout(io.StringIO()):
            return load_slope_data(dst)

    QtCore.QSettings = _Scratch
    try:
        yield load
    finally:
        QtCore.QSettings = real
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# A. Materials
# --------------------------------------------------------------------------- #
_OPTIONS = ["", "mc", "cp", "ratio", "pow", "hb", "elastic"]
_U = ["", "none", "piezo", "seep", "ru"]
_UNSAT = ["", "lf", "vg", "gard"]


def _engine_gray(opt, u, um):
    """The keys of the audit's rules, grayed or not, from what the engines read —
    written out here rather than taken from the editor, so the two are compared.

    u is read for every option, elastic included (Norm, 533a: an elastic material
    reads its own u); ru only when u = ru, on any option; the suction pair only
    when u is piezo or seep and the option takes suction (not cp, ratio or
    elastic). A blank u is none and a blank unsat is lf (the loader's readings).
    kr0 is read by lf; h0 by lf and gard; a and n by vg and gard; l by vg alone.
    Reliability refuses s(γ) on elastic, pow and hb. (s(d) and s(ψ) are left to
    the sheet comparison in table_paste_check.py: whether reliability should read
    them is an open question.)"""
    u = u or "none"
    um = um or "lf"
    gray = {
        "u": False,
        "ru": u != "ru",
        "phi_b": opt in ("cp", "ratio", "elastic") or u not in ("piezo", "seep"),
        "kr0": um in ("vg", "gard"),
        "h0": um == "vg",
        "vg_a": um == "lf",
        "vg_n": um == "lf",
        "vg_l": um in ("lf", "gard"),
        "sigma_gamma": opt in ("elastic", "pow", "hb"),
    }
    gray["s_cap"] = gray["phi_b"]
    return gray


def test_material_graying_follows_the_engine():
    """The table's rule against the engines' reading, for every option, u and unsat
    the editor can meet, the blanks a script or a pasted block can leave included."""
    from studio.editors import _mat_dim_keys
    out = []
    for opt, u, um in itertools.product(_OPTIONS, _U, _UNSAT):
        got = _mat_dim_keys({"option": opt, "u": u, "unsat": um})
        for key, gray in _engine_gray(opt, u, um).items():
            if (key in got) != gray:
                out.append(f"option={opt!r}, u={u!r}, unsat={um!r}: {key} is "
                           f"{'gray' if key in got else 'live'}; the engines "
                           f"{'do not read' if gray else 'read'} it")
    return out[:12] + ([f"... and {len(out) - 12} more"] if len(out) > 12 else [])


def _material_view(row, reliability=True):
    from studio.editors import MaterialsEditor, _MaterialListView, _new_material
    lv = _MaterialListView(MaterialsEditor.FIELDS, [dict(_new_material(), **row)],
                           _new_material, reliability)
    lv._load(0)
    return lv


def test_material_form_follows_the_rule():
    """The list view: gard shows h0 live (and hides kr0); a σ the rule grays is
    grayed by itself while its value stays live — s(γ) on pow — and one it leaves
    live stays live — s(d)/s(ψ) beside a live d/ψ on mc; a blank option still hides
    c and φ (the FEM refuses a blank option)."""
    out = []
    lv = _material_view({"name": "sand", "option": "mc", "unsat": "gard"})
    out += _fail("h0" not in lv._model_hidden and lv._cell_widgets["h0"].isEnabled(),
                 "unsat = gard: the form hides or grays h0, which transient seepage "
                 "reads for Gardner")
    out += _fail("kr0" in lv._model_hidden,
                 "unsat = gard: the form shows kr0, which Gardner does not read")
    out += _fail(lv._cell_widgets["d"].isEnabled()
                 and lv._edits["sigma_d"].isEnabled()
                 and lv._edits["sigma_psi"].isEnabled(),
                 "option = mc: d and s(d)/s(ψ) should be live; d live="
                 f"{lv._cell_widgets['d'].isEnabled()}, s(d) live="
                 f"{lv._edits['sigma_d'].isEnabled()}")
    out += _fail(lv._edits["sigma_gamma"].isEnabled(),
                 "option = mc: s(γ) is gray, but reliability reads it")
    lv._opt_combo.setCurrentText("pow")
    out += _fail(lv._cell_widgets["gamma"].isEnabled()
                 and not lv._edits["sigma_gamma"].isEnabled(),
                 "option = pow: γ should stay live and s(γ) gray")
    # An elastic material reads its own u (533a): u stays live, ru follows u as
    # on any row, and the suction pair stays gray (no strength to add it to).
    lv._opt_combo.setCurrentText("elastic")
    lv._u_combo.setCurrentText("ru")
    out += _fail(lv._cell_widgets["u"].isEnabled()
                 and lv._cell_widgets["ru"].isEnabled(),
                 f"option = elastic, u = ru: u live={lv._cell_widgets['u'].isEnabled()}, "
                 f"ru live={lv._cell_widgets['ru'].isEnabled()}; an elastic material "
                 f"reads both")
    lv._u_combo.setCurrentText("seep")
    out += _fail(not lv._cell_widgets["phi_b"].isEnabled(),
                 "option = elastic, u = seep: φ_b is live, but an elastic material "
                 "takes no suction strength")
    blank = _material_view({"name": "fill", "option": ""})
    out += _fail({"c", "phi"} <= blank._model_hidden,
                 "a blank option no longer hides c and φ in the form")
    return out


# --------------------------------------------------------------------------- #
# B. Reinforcement
# --------------------------------------------------------------------------- #
_NAN = float("nan")
_JOINTS = ["", "No", "Yes", "yes"]
_PULLOUT = {"neither": (_NAN, _NAN), "adhesion only": (5.0, _NAN),
            "delta only": (_NAN, 30.0), "both": (5.0, 30.0)}


def _reinf_rows():
    from studio.editors import _new_reinf
    rows, cases = [], []
    for k, (joint, (law, (adh, dlt))) in enumerate(
            itertools.product(_JOINTS, _PULLOUT.items())):
        rows.append(dict(_new_reinf(), label=f"L{k + 1}", x1=0.0, y1=-2.0 - k * 0.2,
                         x2=10.0, y2=-2.0 - k * 0.2, t_max=50.0, lp1=2.0, lp2=3.0,
                         adhesion=adh, delta=dlt, joint=joint, kn=1e5, ks=1e4,
                         jred="No", spacing=1.0))
        cases.append((joint, law))
    return rows, cases


def _reinf_want(joint, law):
    """Gray per the engines: Lp unread when the pair is filled; Adhesion/Delta read
    (by the pullout law or the FEM interface) or about to be; kn/ks/Jred read only
    on a jointed line."""
    jointed = joint.lower() == "yes"
    return {"lp1": law == "both", "lp2": law == "both", "adhesion": False,
            "delta": False, "kn": not jointed, "ks": not jointed,
            "jred": not jointed}


def test_reinforcement_rule():
    from studio.editors import ReinforcementEditor
    rows, cases = _reinf_rows()
    out = []
    for row, (joint, law) in zip(rows, cases):
        got = ReinforcementEditor.dim_keys(row)
        for key, gray in _reinf_want(joint, law).items():
            out += _fail((key in got) == gray,
                         f"Joint={joint!r}, {law}: {key} is "
                         f"{'gray' if key in got else 'live'}")
    return out


def test_reinforcement_views_agree_and_clear_nothing():
    """Both views, on one line per case: the table's grayed cells and the list's
    dimmed ones are the rule's, row for row; the list's Pullout read-out names the
    law the row's values put in force; and no control in the form other than the
    fields themselves moves an Adhesion or Delta value."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QComboBox
    from studio.editors import ReinforcementEditor

    out = []
    rows, cases = _reinf_rows()
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    with _sandbox() as load:
        sd = load(_LINER)
        sd["reinforcement_lines"] = rows
        dlg = ReinforcementEditor().build(sd, None)
        dlg.set_view_mode("table")
        t = dlg._table.table
        for i, (joint, law) in enumerate(cases):
            for key, gray in _reinf_want(joint, law).items():
                j = keys.index(key)
                w = t.cellWidget(i, j)
                live = (w.isEnabled() if isinstance(w, QComboBox)
                        else bool(t.item(i, j).flags() & Qt.ItemIsEditable))
                out += _fail(live != gray, f"table view, Joint={joint!r}, {law}: "
                                           f"{key} is {'live' if live else 'gray'}")
        dlg.set_view_mode("list")
        lv = dlg._list_view
        for i, (joint, law) in enumerate(cases):
            lv.list.setCurrentRow(i)
            for key, gray in _reinf_want(joint, law).items():
                live = lv._cells[key].isEnabled()
                out += _fail(live != gray, f"list view, Joint={joint!r}, {law}: "
                                           f"{key} is {'live' if live else 'gray'}")
            want = ("Overburden (Adhesion, Delta)" if law == "both"
                    else "Development length (Lp1, Lp2)")
            shown = getattr(lv, "_readout_label", None)
            out += _fail(shown is not None and shown.text() == want,
                         f"list view, {law}: the Pullout read-out shows "
                         f"{shown.text() if shown is not None else None!r}, "
                         f"not {want!r}")
            # Every combo in the form that is not a field's own editor (the old
            # Pullout switch was one) is driven through all its entries; Adhesion
            # and Delta must come through untouched.
            # Checked after every entry: the old switch cleared the pair on one
            # entry and put it back on the next.
            before = (lv._edits["adhesion"].text(), lv._edits["delta"].text())
            fields = set(map(id, lv._edits.values()))
            moved = None
            for combo in lv._form_scroll.findChildren(QComboBox):
                if id(combo) in fields:
                    continue
                for k in range(combo.count()):
                    combo.setCurrentIndex(k)
                    now = (lv._edits["adhesion"].text(), lv._edits["delta"].text())
                    if now != before and moved is None:
                        moved = (combo.itemText(k), now)
            out += _fail(moved is None,
                         f"list view, Joint={joint!r}, {law}: picking "
                         f"{moved[0] if moved else ''!r} changed Adhesion/Delta from "
                         f"{before} to {moved[1] if moved else ''}")
        # Typing the second of the pair moves the line onto the overburden law
        # in the list view, at once.
        i = cases.index(("Yes", "adhesion only"))
        lv.list.setCurrentRow(i)
        lv._edits["delta"].setText("28")
        lv._edits["delta"].editingFinished.emit()
        out += _fail(not lv._cells["lp1"].isEnabled()
                     and lv._readout_label.text().startswith("Overburden"),
                     "list view: filling Delta beside Adhesion left Lp1 live or the "
                     "read-out on the development lengths")
        lv._edits["joint"].setCurrentText("No")
        out += _fail(not lv._cells["kn"].isEnabled(),
                     "list view: kn stayed live after Joint was set to No")
        dlg.deleteLater()
    return out


def test_a_new_line_has_blank_development_lengths():
    """A new line's Lp1/Lp2 are blank, as an empty sheet cell is: 0 would be a fully
    anchored end that passes every check. Through the table's Add row and OK."""
    from studio.editors import ReinforcementEditor, _new_reinf
    out = []
    fresh = _new_reinf()
    out += _fail(_blank(fresh["lp1"]) and _blank(fresh["lp2"]),
                 f"a new line opens with Lp1={fresh['lp1']!r}, Lp2={fresh['lp2']!r}")
    with _sandbox() as load:
        sd = load(_LINER)
        sd["reinforcement_lines"] = []
        editor = ReinforcementEditor()
        dlg = editor.build(sd, None)
        dlg.set_view_mode("table")
        dlg._table._add_row()
        editor.apply(sd, dlg)
        got = sd["reinforcement_lines"][0]
        out += _fail(_blank(got["lp1"]) and _blank(got["lp2"]),
                     f"a line added in the table saves Lp1={got['lp1']!r}, "
                     f"Lp2={got['lp2']!r}")
        dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# C. The mesh signature
# --------------------------------------------------------------------------- #
def _signature_model():
    from shapely.geometry import Polygon
    return {
        "profile_lines": [], "max_depth": None,
        "polygons": [{"polygon": Polygon([(0, 0), (40, 0), (40, 10), (0, 10)]),
                      "mat_id": 0, "size": None}],
        "reinforcement_lines": [
            {"x1": 5.0, "y1": 5.0, "x2": 20.0, "y2": 5.0, "t_max": 30.0,
             "tend1": 0.0, "tend2": 0.0, "joint": "", "spacing": 1.0},
            {"x1": 5.0, "y1": 7.0, "x2": 20.0, "y2": 7.0, "t_max": 30.0,
             "tend1": 0.0, "tend2": 0.0, "joint": "Yes", "spacing": 1.0}],
        "pile_lines": [], "joint_lines": [], "refine_zones": [], "joint_zones": [],
        "ssr_zones": [{"kind": "ssr_reduce", "label": "reduce",
                       "polygon": [(10, 2), (30, 2), (30, 8), (10, 8)], "size": 1.0}],
        "dloads": [[{"X": 0.0, "Y": 10.0, "Normal": 100.0},
                    {"X": 10.0, "Y": 10.0, "Normal": 100.0}]],
        "dloads2": [],
        "seepage_bc": {"specified_heads": [{"head": 9.0,
                                            "coords": [(0, 0), (0, 9)]}],
                       "specified_fluxes": [], "exit_face": [(40, 0), (40, 4)]},
        "seepage_bc2": {"specified_heads": [], "specified_fluxes": [],
                        "exit_face": []},
    }


def test_mesh_signature_sees_what_the_mesher_reads():
    from studio.main_window import MainWindow
    sig = MainWindow.mesh_signature
    out = []

    def edit(fn):
        sd = _signature_model()
        fn(sd)
        return sd

    base = sig(_signature_model())
    out += _fail(base is not None, "the mesh signature failed on a plain model")
    stale = {
        "a line's Joint set to Yes":
            lambda sd: sd["reinforcement_lines"][0].update(joint="Yes"),
        "a jointed line's Tend1": lambda sd: sd["reinforcement_lines"][1].update(tend1=20.0),
        "a distributed-load vertex": lambda sd: sd["dloads"][0][1].update(X=12.0),
        "a seepage head BC vertex":
            lambda sd: sd["seepage_bc"]["specified_heads"][0].update(
                coords=[(0, 0), (0, 8)]),
        "the exit face": lambda sd: sd["seepage_bc"].update(exit_face=[(40, 0), (40, 5)]),
        "an SSR zone's Size": lambda sd: sd["ssr_zones"][0].update(size=0.5),
        "a polygon's Size": lambda sd: sd["polygons"][0].update(size=2.0),
    }
    for what, fn in stale.items():
        out += _fail(sig(edit(fn)) != base,
                     f"the mesh is kept when {what} changes, which the mesher reads")
    kept = {
        "a bonded line's Tend1": lambda sd: sd["reinforcement_lines"][0].update(tend1=20.0),
        "a Tmax": lambda sd: sd["reinforcement_lines"][1].update(t_max=45.0),
        "a load's magnitude": lambda sd: sd["dloads"][0][0].update(Normal=250.0),
        "a head value": lambda sd: sd["seepage_bc"]["specified_heads"][0].update(head=7.0),
        "an SSR zone's kind": lambda sd: sd["ssr_zones"][0].update(kind="ssr_hold"),
        "Joint 'Yes' read back as 'yes'":
            lambda sd: sd["reinforcement_lines"][1].update(joint="yes"),
    }
    for what, fn in kept.items():
        out += _fail(sig(edit(fn)) == base,
                     f"the mesh is thrown away when {what} changes, which the "
                     f"mesher does not read")

    # Automatic water loads: the derived pool's vertices are pinned by the mesher,
    # so moving the water that derives them makes the mesh stale. On the dam the
    # pool stands at the reservoir's seepage head, so a head VALUE — not mesh
    # input on its own (above) — moves the mesh's pinned points here.
    from xslope.water import DERIVED_KEYS, water_loads_mode, with_water_loads
    with _sandbox() as load:
        dam = load(_DAM)
    if water_loads_mode(dam) != "auto" or not with_water_loads(dam).get(DERIVED_KEYS[1]):
        out.append("xslope_dam.xlsx no longer derives a water load; pick another "
                   "model for the derived-load leg")
        return out
    before = sig(dam)
    raised = copy.deepcopy(dam)
    head = raised["seepage_bc"]["specified_heads"][0]
    head["head"] = float(head["head"]) - 10.0
    moved = with_water_loads(raised).get(DERIVED_KEYS[1])
    out += _fail(moved != with_water_loads(dam).get(DERIVED_KEYS[1]),
                 "lowering the reservoir head did not move the derived water load")
    out += _fail(sig(raised) != before,
                 "the mesh is kept when the derived water load moves")
    return out


# --------------------------------------------------------------------------- #
# D. Non-circular Movement
# --------------------------------------------------------------------------- #
def test_noncirc_movement_is_shown_and_saved_as_read():
    from studio.editors import NonCircEditor
    out = []
    pts = [(0.0, _NAN), (10.0, _NAN), (20.0, "free"), (30.0, "slide"),
           (40.0, "Horiz"), (50.0, None)]
    rows = [{"X": x, "Y": 5.0 - 0.1 * x, "Movement": mv} for x, mv in pts]
    with _sandbox() as load:
        sd = load(_SIMPLE)
        sd["non_circ"] = copy.deepcopy(rows)
        editor = NonCircEditor()
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        col = [f.key for f in NonCircEditor.FIELDS].index("Movement")
        shown = [t.cellWidget(i, col).currentText() for i in range(len(rows))]
        holder = [t.cellWidget(i, col).placeholderText() for i in range(len(rows))]
        out += _fail(shown == ["", "", "Free", "slide", "Horiz", ""],
                     f"the Movement cells show {shown}: a blank must stay blank and "
                     f"an unknown word as typed")
        out += _fail(holder[0] == "Free (blank)" and holder[1] == "Fixed (blank)"
                     and holder[5] == "Free (blank)",
                     f"the blank cells' placeholders read {holder[0]!r}, "
                     f"{holder[1]!r}, {holder[5]!r}; the search slides a blank end "
                     f"point and holds a blank interior one")
        editor.apply(sd, dlg)
        saved = [p.get("Movement") for p in sd["non_circ"]]
        same = all((_blank(a) and _blank(b)) or a == b
                   for a, b in zip(saved, [mv for _, mv in pts]))
        out += _fail(same, f"OK on an untouched table saved Movement {saved}; the "
                           f"file held {[mv for _, mv in pts]}")
        dlg.deleteLater()

        # A real pick is saved; an added row makes the old last row interior.
        sd["non_circ"] = copy.deepcopy(rows)
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        t.cellWidget(3, col).setCurrentText("Fixed")
        dlg._editable._add_row()
        out += _fail(t.cellWidget(5, col).placeholderText() == "Fixed (blank)",
                     "after a row was added, the former end point's blank still "
                     "reads as an end point")
        editor.apply(sd, dlg)
        out += _fail(sd["non_circ"][3]["Movement"] == "Fixed",
                     f"picking Fixed over 'slide' saved "
                     f"{sd['non_circ'][3]['Movement']!r}")
        dlg.deleteLater()

    # The preview's reading against the search's own, where the search has one.
    try:
        from studio.editors import _noncirc_point_movement
    except ImportError:
        return out + ["studio.editors has no _noncirc_point_movement, the preview's "
                      "reading of Movement"]
    try:
        from xslope.search import noncirc_point_movements
    except ImportError:
        return out + ["xslope.search.noncirc_point_movements is gone; the preview "
                      "has no engine reading to be held to"]
    known = [r for r in rows if r["Movement"] != "slide"]
    engine = noncirc_point_movements(known)
    ours = [_noncirc_point_movement(r["Movement"], i, len(known))
            for i, r in enumerate(known)]
    out += _fail(ours == engine, f"the preview reads Movement {ours}; the search "
                                 f"reads {engine}")
    out += _fail(_noncirc_point_movement("slide", 2, 5) is None,
                 "the preview reads an unknown Movement word as a movement")
    return out


# --------------------------------------------------------------------------- #
# E. Circles
# --------------------------------------------------------------------------- #
_CIRCLE_GRAY = {"Depth": {"Xi", "Yi", "R"}, "Radius": {"Depth", "Xi", "Yi"},
                "Intercept": {"Depth", "R"}}


def test_circles_option_grays_the_unread_fields():
    from PySide6.QtCore import Qt
    from studio.editors import CirclesEditor
    out = []
    keys = [f.key for f in CirclesEditor.FIELDS]
    rows = [{"Xo": 10.0, "Yo": 40.0, "Option": opt, "Depth": 0.0, "Xi": 0.0,
             "Yi": 0.0, "R": 40.0} for opt in _CIRCLE_GRAY]
    with _sandbox() as load:
        sd = load(_SIMPLE)
        sd["circles"] = rows
        dlg = CirclesEditor().build(sd, None)
        t = dlg._editable.table

        def gray(i):
            return {k for k in ("Depth", "Xi", "Yi", "R")
                    if not t.item(i, keys.index(k)).flags() & Qt.ItemIsEditable}

        for i, opt in enumerate(_CIRCLE_GRAY):
            out += _fail(gray(i) == _CIRCLE_GRAY[opt],
                         f"Option={opt}: grays {sorted(gray(i))}, the sheet grays "
                         f"{sorted(_CIRCLE_GRAY[opt])}")
        t.cellWidget(0, keys.index("Option")).setCurrentText("Radius")
        out += _fail(gray(0) == _CIRCLE_GRAY["Radius"],
                     f"switching a Depth circle to Radius left {sorted(gray(0))} gray")
        dlg.deleteLater()
    return out


def _same_cells(a, b):
    """Two size-cell dicts hold the same numbers, a blank matching a blank."""
    if not isinstance(a, dict) or set(a) != set(b):
        return False
    return all((_blank(a[k]) and _blank(b[k])) or a[k] == b[k] for k in b)


def test_circles_blank_option_stays_blank():
    """A circle with a blank Option (the loader keeps it, unsized) is shown blank,
    grays nothing — the sheet grays nothing for it — and is saved blank and unsized
    unless an Option is picked: a Depth typed beside a blank Option must not size it.
    A word that is none of the three is shown and saved as typed."""
    from PySide6.QtCore import Qt
    from studio.editors import CirclesEditor
    out = []
    keys = [f.key for f in CirclesEditor.FIELDS]
    col = keys.index("Option")
    rows = [{"Xo": 10.0, "Yo": 40.0, "Option": None, "Depth": _NAN, "R": _NAN},
            {"Xo": 12.0, "Yo": 40.0, "Option": "dpth", "Depth": _NAN, "R": _NAN},
            {"Xo": 14.0, "Yo": 40.0, "Depth": 0.0, "R": 40.0}]
    with _sandbox() as load:
        sd = load(_SIMPLE)
        sd["circles"] = copy.deepcopy(rows)
        editor = CirclesEditor()
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        shown = [t.cellWidget(i, col).currentText() for i in range(3)]
        out += _fail(shown == ["", "dpth", "Depth"],
                     f"the Option cells show {shown}: a blank Option must stay blank "
                     f"and an unknown word as typed")
        for i in (0, 1):
            gray = [k for k in ("Depth", "Xi", "Yi", "R")
                    if not t.item(i, keys.index(k)).flags() & Qt.ItemIsEditable]
            out += _fail(not gray, f"row {i + 1} (Option {rows[i]['Option']!r}) "
                                   f"grays {gray}; the sheet grays nothing for it")
        # A Depth typed beside the blank Option is not a size, but it is kept.
        t.item(0, keys.index("Depth")).setText("5")
        editor.apply(sd, dlg)
        got = sd["circles"]
        out += _fail(_blank(got[0].get("Option")) and _blank(got[0]["R"])
                     and _blank(got[0]["Depth"]),
                     f"the blank-Option circle saved Option={got[0].get('Option')!r}, "
                     f"R={got[0]['R']!r}, Depth={got[0]['Depth']!r}; it states no size")
        cells = got[0].get("size_cells") or {}
        out += _fail(cells.get("Depth") == 5.0
                     and all(_blank(cells.get(k)) for k in ("Xi", "Yi", "R")),
                     f"the Depth typed beside a blank Option was not kept: size cells "
                     f"{cells!r}; want Depth 5 and the untyped cells blank")
        out += _fail(got[1].get("Option") == "dpth" and _blank(got[1]["R"]),
                     f"the unknown Option saved as {got[1].get('Option')!r}, "
                     f"R={got[1]['R']!r}")
        out += _fail(got[2]["R"] == 40.0, f"the Depth circle saved R={got[2]['R']!r}")
        dlg.deleteLater()

        # A loaded blank-Option circle shows the size cells it was saved with, and
        # picking an Option sizes it from the one that Option names.
        kept = {"Xo": 10.0, "Yo": 40.0, "Option": None, "Depth": _NAN, "R": _NAN,
                "size_cells": {"Depth": 8.0, "Xi": 3.0, "Yi": _NAN, "R": 31.0}}
        sd["circles"] = [copy.deepcopy(kept)]
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        shown = [t.item(0, keys.index(k)).text() for k in ("Depth", "Xi", "Yi", "R")]
        out += _fail(shown == ["8", "3", "", "31"],
                     f"a blank-Option circle's kept size cells show as {shown}")
        editor.apply(sd, dlg)
        out += _fail(_same_cells(sd["circles"][0].get("size_cells"),
                                 kept["size_cells"]),
                     f"OK on an untouched blank-Option circle changed its size cells "
                     f"to {sd['circles'][0].get('size_cells')!r}")
        dlg.deleteLater()
        sd["circles"] = [copy.deepcopy(kept)]
        dlg = editor.build(sd, None)
        dlg._editable.table.cellWidget(0, col).setCurrentText("Depth")
        editor.apply(sd, dlg)
        got = sd["circles"][0]
        out += _fail((got["Depth"], got["R"]) == (8.0, 32.0)
                     and "size_cells" not in got,
                     f"picking Depth on a circle carrying Depth = 8 saved Depth="
                     f"{got['Depth']!r}, R={got['R']!r}, size cells "
                     f"{got.get('size_cells')!r}")
        dlg.deleteLater()

        # A circle whose Option names a blank cell (Radius, R empty) is unsized
        # too: its Option and typed cells survive OK untouched, and filling R
        # sizes it.
        radius = {"Xo": 10.0, "Yo": 40.0, "Option": "Radius", "Depth": _NAN,
                  "R": _NAN,
                  "size_cells": {"Depth": 8.0, "Xi": 3.0, "Yi": _NAN, "R": _NAN}}
        sd["circles"] = [copy.deepcopy(radius)]
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        shown = ([t.cellWidget(0, col).currentText()]
                 + [t.item(0, keys.index(k)).text() for k in ("Depth", "Xi", "Yi", "R")])
        out += _fail(shown == ["Radius", "8", "3", "", ""],
                     f"a Radius circle with R blank shows as {shown}")
        editor.apply(sd, dlg)
        got = sd["circles"][0]
        out += _fail(got.get("Option") == "Radius" and _blank(got["R"])
                     and _same_cells(got.get("size_cells"), radius["size_cells"]),
                     f"OK on an untouched Radius circle with R blank saved Option="
                     f"{got.get('Option')!r}, R={got['R']!r}, size cells "
                     f"{got.get('size_cells')!r}")
        dlg.deleteLater()
        sd["circles"] = [copy.deepcopy(radius)]
        dlg = editor.build(sd, None)
        dlg._editable.table.item(0, keys.index("R")).setText("30")
        editor.apply(sd, dlg)
        got = sd["circles"][0]
        out += _fail((got["R"], got["Depth"]) == (30.0, 10.0)
                     and "size_cells" not in got,
                     f"filling R on that circle saved R={got['R']!r}, Depth="
                     f"{got['Depth']!r}, size cells {got.get('size_cells')!r}")
        dlg.deleteLater()

        # Picking Intercept with only Xi typed leaves it unsized, Xi kept.
        sd["circles"] = copy.deepcopy(rows[:1])
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        t.cellWidget(0, col).setCurrentText("Intercept")
        t.item(0, keys.index("Xi")).setText("6")
        editor.apply(sd, dlg)
        got = sd["circles"][0]
        out += _fail(got.get("Option") == "Intercept" and _blank(got["R"])
                     and (got.get("size_cells") or {}).get("Xi") == 6.0,
                     f"Intercept picked with Yi blank saved Option="
                     f"{got.get('Option')!r}, R={got['R']!r}, size cells "
                     f"{got.get('size_cells')!r}")
        dlg.deleteLater()

        # Picking an Option sizes the circle from the cell it names.
        sd["circles"] = copy.deepcopy(rows)
        dlg = editor.build(sd, None)
        t = dlg._editable.table
        t.cellWidget(0, col).setCurrentText("Radius")
        t.item(0, keys.index("R")).setText("30")
        editor.apply(sd, dlg)
        got = sd["circles"][0]
        out += _fail((got.get("Option"), got["R"], got["Depth"]) == ("Radius", 30.0, 10.0),
                     f"picking Radius with R = 30 saved {got.get('Option')!r}, "
                     f"R={got['R']!r}, Depth={got['Depth']!r}")
        dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# F. Polygons
# --------------------------------------------------------------------------- #
def test_polygon_size_is_live_for_every_type():
    from studio.editors import PolygonEditor
    out = []
    with _sandbox() as load:
        sd = load(_POLYGONS)
        ring = [tuple(c) for c in sd["polygons"][0]["polygon"].exterior.coords][:-1]
        sd["ssr_zones"] = [{"kind": k, "label": k, "polygon": ring, "size": 1.0}
                           for k in ("ssr_reduce", "ssr_hold", "ssr_elastic")]
        sd["refine_zones"] = [{"polygon": ring, "size": 1.0}]
        dlg = PolygonEditor().build(sd, None)
        for i in range(dlg.list.count()):
            dlg.list.setCurrentRow(i)
            out += _fail(dlg.size_edit.isEnabled(),
                         f"polygon {i + 1} ({dlg._kind_at()}): Size is grayed")
        dlg.list.setCurrentRow(0)
        for k in range(dlg.type_combo.count()):
            dlg.type_combo.setCurrentIndex(k)
            out += _fail(dlg.size_edit.isEnabled(),
                         f"Type {dlg.type_combo.currentText()!r}: Size is grayed")
        dlg.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# G. Seepage boundary values
# --------------------------------------------------------------------------- #
def test_seep_bc_blank_value_stays_blank():
    """A seepage boundary whose value was never entered keeps a blank value through
    the editor, on both boundary sets, head and flux alike — so the model checks
    can name it (seep.bc_value_missing) instead of a head or flux of 0 standing in
    for it. A value emptied by hand saves blank too, and a new boundary starts
    blank."""
    from studio.editors import _SeepBcSetWidget
    out = []
    for set_no, constant in ((1, False), (2, True)):
        bc = {"specified_heads": [{"head": None, "kind": "head",
                                   "coords": [(0.0, 0.0), (0.0, 9.0)]},
                                  {"head": 9.0, "kind": "head",
                                   "coords": [(40.0, 0.0), (40.0, 9.0)]}],
              "specified_fluxes": [{"flux": None, "coords": [(10.0, 10.0), (20.0, 10.0)]}],
              "exit_face": [(40.0, 0.0), (40.0, 4.0)]}
        w = _SeepBcSetWidget(copy.deepcopy(bc), constant_only=constant, set_no=set_no)
        n = w.list.count()
        for row in range(n):              # visit every boundary: each one commits
            w.list.setCurrentRow(row)
        w.list.setCurrentRow(1)           # empty the filled head by hand
        w.head_edit.setText("")
        w.list.setCurrentRow(n - 1)
        got = w.result()
        heads = [h["head"] for h in got["specified_heads"]]
        flux = got["specified_fluxes"][0]["flux"]
        out += _fail(_blank(heads[0]) and _blank(flux),
                     f"set {set_no}: an untouched blank head/flux saved as "
                     f"{heads[0]!r}/{flux!r}, not blank")
        out += _fail(_blank(heads[1]),
                     f"set {set_no}: a head emptied by hand saved as {heads[1]!r}, "
                     f"not blank")
        w._add_head()
        w._add_flux()
        w.list.setCurrentRow(w.list.count() - 1)
        got = w.result()
        out += _fail(_blank(got["specified_heads"][-1]["head"])
                     and _blank(got["specified_fluxes"][-1]["flux"]),
                     f"set {set_no}: a new boundary starts at head "
                     f"{got['specified_heads'][-1]['head']!r}, flux "
                     f"{got['specified_fluxes'][-1]['flux']!r}, not blank")
        w.deleteLater()
    return out


# --------------------------------------------------------------------------- #
# H. Tension SRF
# --------------------------------------------------------------------------- #
def test_tension_srf_counts_joint_lines():
    """Run FEM's Tension SRF box is live on a strength reduction whenever the
    setting changes a number: a material OR a joints-sheet line with a tensile
    cutoff above zero (a joint's t_cut is divided by the trial F under the same
    setting). Dimmed when every cutoff is blank or 0, and the tooltip says why."""
    from studio.dialogs import RunFemDialog
    out = []
    with _sandbox() as load:
        sd = load(_ROCK_JOINTS)
        if not sd.get("joint_lines"):
            return ["xslope_rock_joints.xlsx no longer carries joints-sheet lines"]
        for m in sd.get("materials") or []:
            m["t_cut"] = 0.0
        for j in sd["joint_lines"]:
            j["t_cut"] = 0.0

        def box(model):
            dlg = RunFemDialog(defaults={"analysis": "ssrm"}, slope_data=model)
            dlg._sync_enabled()
            state = (dlg.tension_srf.isEnabled(), dlg.tension_srf.toolTip())
            dlg.deleteLater()
            return state

        live, tip = box(sd)
        out += _fail(not live, "every t_cut is 0, yet the Tension SRF box is live")
        out += _fail("joint" in tip.lower(),
                     "the dimmed box's tooltip names materials only")
        sd["joint_lines"][0]["t_cut"] = 5.0
        live, _ = box(sd)
        out += _fail(live, "a joints-sheet line has t_cut = 5, which the Tension SRF "
                           "setting divides by F, yet the box is dimmed")
    return out


# --------------------------------------------------------------------------- #
# I. Line loads in Studio's mesh
# --------------------------------------------------------------------------- #
def test_studio_mesh_pins_line_loads():
    """Studio's Build mesh passes each line load's point of application to the
    mesher as a point constraint, as every scripted build does — so the FEM finds
    a node at the load rather than moving it to the nearest one with a warning —
    and passes none on a model without line loads, whose mesh is built exactly as
    before. Moving a line load makes the mesh stale; changing its magnitude or
    angle does not. The builder is stood in for, so nothing is meshed here."""
    import xslope.mesh as mesh_mod
    from studio.main_window import MainWindow
    from studio.runners import MeshWorker
    out = []
    with _sandbox() as load:
        sd = load(_SIMPLE)
    gs = list(sd["ground_surface"].coords)
    k = max(range(len(gs) - 1), key=lambda i: abs(gs[i + 1][0] - gs[i][0]))
    (x0, y0), (x1, y1) = gs[k], gs[k + 1]
    point = (x0 + 0.37 * (x1 - x0), y0 + 0.37 * (y1 - y0))   # mid-segment, no vertex

    seen = []
    real = mesh_mod.build_mesh_from_polygons

    def recorder(polygons, **kwargs):
        seen.append(kwargs)
        return {"nodes": [], "elements": [], "element_types": [],
                "element_materials": []}

    def build(model):
        worker = MeshWorker()
        errors = []
        worker.failed.connect(errors.append)
        with contextlib.redirect_stdout(io.StringIO()):
            worker.build(model, {"element_type": "tri6"})
        return errors

    mesh_mod.build_mesh_from_polygons = recorder
    try:
        plain = copy.deepcopy(sd)
        plain["line_loads"] = []
        errors = build(plain)
        out += _fail(not errors and seen and seen[-1].get("point_constraints") is None,
                     f"a model with no line load reached the mesher with "
                     f"point_constraints={seen[-1].get('point_constraints') if seen else '?'!r}"
                     f"{' (' + errors[0] + ')' if errors else ''}; it must mesh as before")
        loaded = copy.deepcopy(sd)
        loaded["line_loads"] = [{"x": point[0], "y": point[1], "P": 100.0,
                                 "angle": -90.0, "label": "L1"}]
        errors = build(loaded)
        got = seen[-1].get("point_constraints") if seen else None
        out += _fail(not errors and got and [tuple(p) for p in got] == [point],
                     f"Studio's mesh build passed point_constraints={got!r} for a line "
                     f"load at {point}; the mesher puts no node there without it")
    finally:
        mesh_mod.build_mesh_from_polygons = real

    sig = MainWindow.mesh_signature
    base = sig(loaded)
    moved = copy.deepcopy(loaded)
    moved["line_loads"][0]["x"] += 1.0
    out += _fail(sig(moved) != base, "moving a line load keeps the mesh, which has a "
                                     "node at its old point and none at the new one")
    pushed = copy.deepcopy(loaded)
    pushed["line_loads"][0].update(P=250.0, angle=-60.0)
    out += _fail(sig(pushed) == base, "changing a line load's magnitude or angle "
                                      "throws the mesh away")
    return out


# --------------------------------------------------------------------------- #
# J. The Reliability dialog's checks judge the run it will make
# --------------------------------------------------------------------------- #
def test_reliability_checks_see_engine_and_rapid():
    """The Reliability dialog's model checks are run for the run it will make:
    its engine and whether it is a rapid drawdown. An mc material carrying s(d)
    beside s(φ): a Monte Carlo rapid drawdown run reads s(d), so the panel shows no
    reliability.dead_sigma_columns finding; unchecking Rapid drawdown makes it the
    INFO that s(d) is left out; and the panel follows both controls as they move.
    The engine reaches the checks too (the Taylor series refuses a σ above its
    mean, the sampling engines truncate their draws instead)."""
    from studio.dialogs import ReliabilityDialog
    out = []
    with _sandbox() as load:
        sd = load(_SIMPLE)
        mats = sd["materials"]
        mats[0].update(option="mc", sigma_phi=2.0, sigma_d=50.0, sigma_psi=0.0)

        def dead(dlg):
            rep = dlg.preflight.report
            return [f.severity for f in (rep.findings if rep is not None else ())
                    if f.rule_id == "reliability.dead_sigma_columns"]

        dlg = ReliabilityDialog(defaults={"engine": "mc", "rapid": True},
                                slope_data=sd, app_mode="lem")
        out += _fail(dead(dlg) == [],
                     f"Monte Carlo with rapid drawdown: the panel reports "
                     f"dead_sigma_columns {dead(dlg)}, but that run reads s(d)")
        dlg.rapid.setChecked(False)
        out += _fail(dead(dlg) == ["info"],
                     f"Monte Carlo without rapid drawdown: dead_sigma_columns "
                     f"{dead(dlg)}, expected the INFO that s(d) is left out")
        dlg.rapid.setChecked(True)
        out += _fail(dead(dlg) == [], "checking Rapid drawdown again did not clear "
                                      "the finding: the panel does not follow it")

        # The engine: a σ above its mean is refused by the Taylor series only.
        mats[0]["sigma_phi"] = float(mats[0].get("phi") or 0.0) + 5.0
        dlg2 = ReliabilityDialog(defaults={"engine": "taylor"}, slope_data=sd,
                                 app_mode="lem")

        def exceeds(d):
            rep = d.preflight.report
            return any(f.rule_id == "reliability.sigma_exceeds_mean"
                       for f in (rep.findings if rep is not None else ()))
        out += _fail(exceeds(dlg2), "the Taylor series run shows no "
                                    "sigma_exceeds_mean for s(φ) above φ")
        dlg2.engine.setCurrentIndex(dlg2.engine.findData("mc"))
        out += _fail(not exceeds(dlg2), "switching the engine to Monte Carlo left "
                                        "sigma_exceeds_mean standing; the sampling "
                                        "engines truncate their draws instead")
        dlg.deleteLater()
        dlg2.deleteLater()
    return out


CHECKS = [
    ("material graying = engine", test_material_graying_follows_the_engine),
    ("material form = the rule", test_material_form_follows_the_rule),
    ("reinforcement rule", test_reinforcement_rule),
    ("reinforcement views agree", test_reinforcement_views_agree_and_clear_nothing),
    ("new line: Lp blank", test_a_new_line_has_blank_development_lengths),
    ("mesh signature", test_mesh_signature_sees_what_the_mesher_reads),
    ("non-circ Movement as read", test_noncirc_movement_is_shown_and_saved_as_read),
    ("circles Option graying", test_circles_option_grays_the_unread_fields),
    ("circles blank Option", test_circles_blank_option_stays_blank),
    ("polygon Size live", test_polygon_size_is_live_for_every_type),
    ("seep BC blank value", test_seep_bc_blank_value_stays_blank),
    ("Tension SRF counts joints", test_tension_srf_counts_joint_lines),
    ("Studio mesh pins line loads", test_studio_mesh_pins_line_loads),
    ("Reliability checks see the run", test_reliability_checks_see_engine_and_rapid),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("audit studio: PySide6 not installed — skipped.")
        return []
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    failures = []
    for name, fn in CHECKS:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:34s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += [f"{name}: {f}" for f in fs]
    return failures


def main():
    print("audit studio:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll audit studio checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
