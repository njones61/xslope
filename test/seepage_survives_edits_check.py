"""Standing check: an input edit keeps the mesh and the seepage solutions it cannot
change, and drops the ones it can.

Every editor OK used to drop every result, and every seepage-boundary edit the
mesh. So a strength edit took a loaded transient seepage solution with it (COMBO-3:
Run LEM then had no frame to read), and SEEP-4's three rain blocks, every vertex of
which the mesh already had, left Run Seep asking for a new mesh.

**The mesh.** A boundary or load vertex reaches the mesh only one way: the mesher
inserts it into the zone polygon whose edge it falls on, so the mesh gets a node
there (``mesh.get_material_polygons``). The window's mesh signature now reads those
polygons rather than the raw vertex lists. A block drawn over vertices the zones
already have changes nothing and keeps the mesh; one that lands part way along an
edge, or the removal of one that did, changes the zones and drops it. Each case is
held to an independent build: the mesh built after the edit is node for node the
kept one, or it is not.

**The seepage solutions.** Each steady boundary set's solution and the transient
march record what they were solved from — the mesh object and every model input
the solve reads (``MainWindow.seepage_inputs``: each material's conductivities and
unsaturated parameters, the storage for the march, the boundary set, the tseep
sheet as the solve reads it, the unit weight of water, the units). An edit keeps a
solution whose record still matches, and drops one whose record does not. Driven
through the real edit path (``MainWindow.edit_category``) on the transient earth
dam (COMBO-3) and the solved Johnson dam (SEEP-2), with a real edit each: a
strength, a name, a circle, a seismic coefficient, a stage time and an undo keep
the solution; a conductivity, a storage, the unit weight of water, a head, the
duration and a new mesh drop it. Every other result (the stability and finite
element runs, and a Monte Carlo reliability result) is still dropped by every edit.

**The field a dropped solution leaves.** A solve puts its pore pressures into the
model (``seep_u``, ``seep_u2``; a run on a transient solution, the frame it reads
and its instant), and they used to outlive the solution: after a k1 edit the model
checks passed ``u = seep`` on the old field and Save rewrote ``_seep.csv`` from
it. A dropped solution now takes its field with it, Save deletes the companion,
and the file reopens with no seepage solution; a field no Studio solve put there
stays, and a kept steady solution keeps its field in the model. With the field
gone, Run FEM's model checks refuse the run (seep_field.missing, 602a) instead of
letting it go on at zero pore pressure.

Every window is built against a scratch settings store. Skips cleanly (exit 0)
when PySide6 is not installed.
"""
import contextlib
import copy
import io
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from tools import studio_settings                      # noqa: E402

_FILES = os.path.join(_REPO, "docs/tutorials/files")
#: SEEP-4: the completed dam (a head on the upstream face, three rain blocks over
#: the faces and crest, an exit face), meshed as the page meshes it.
SEEP04 = os.path.join(_FILES, "xslope_dam_infiltration.xlsx")
SEEP04_MESH = {"element_type": "tri3", "auto_size": False, "target_size": 1.0}
#: COMBO-2: set 2's head ends at (220, 110), part way up the shell face.
COMBO02 = os.path.join(_FILES, "xslope_johnson_rapid.xlsx")
COMBO02_MESH = {"element_type": "tri3", "auto_size": True, "size_divisions": 100}
#: COMBO-3's transient earth dam, saved with its mesh and transient solution.
TRANSIENT = os.path.join(_FILES, "xslope_earth_dam_fs_time.xlsx")
#: Its starter: the same dam and solution, the strengths still to be entered.
TRANSIENT_START = os.path.join(_FILES, "xslope_earth_dam_fs_time_start.xlsx")
#: SEEP-2's Johnson dam, saved with its mesh and steady solution (set 1).
STEADY = os.path.join(_FILES, "xslope_johnson_res_solved.xlsx")


def _fail(cond, message):
    return [] if cond else [message]


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _load(path):
    from xslope.fileio import load_slope_data
    return _quiet(load_slope_data, path)


def _copy_model(src, tmp):
    """Copy a model and every companion file beside it into ``tmp``."""
    base = os.path.splitext(os.path.basename(src))[0]
    for f in os.listdir(os.path.dirname(src)):
        if f.startswith(base + ".") or f.startswith(base + "_"):
            shutil.copy(os.path.join(os.path.dirname(src), f), tmp)
    return os.path.join(tmp, os.path.basename(src))


def _build_mesh(sd, opts):
    """The mesh Studio's Build mesh would make of ``sd`` (runners.MeshWorker)."""
    from studio.runners import MeshWorker
    got = {}
    w = MeshWorker()
    w.succeeded.connect(lambda m: got.setdefault("mesh", m))
    w.failed.connect(lambda e: got.setdefault("error", e))
    _quiet(w.build, sd, dict(opts))
    if "mesh" not in got:
        raise RuntimeError(f"mesh build failed: {got.get('error')}")
    return got["mesh"]


def _same_mesh(a, b):
    import numpy as np
    return (a["nodes"].shape == b["nodes"].shape
            and a["elements"].shape == b["elements"].shape
            and np.array_equal(a["nodes"], b["nodes"])
            and np.array_equal(a["elements"], b["elements"]))


def _set_bc(sd, bc, bc2=None):
    sd["seepage_bc"] = copy.deepcopy(bc)
    if bc2 is not None:
        sd["seepage_bc2"] = copy.deepcopy(bc2)
        sd["has_seepage_bc2"] = bool(bc2.get("specified_heads") or bc2.get("exit_face"))


# --------------------------------------------------------------------------- #
# The mesh
# --------------------------------------------------------------------------- #
def test_mesh_kept_for_blocks_on_its_vertices(mw, tmp):
    """SEEP-4's rain: the mesh built with the head and the exit face is kept when
    the three flux blocks are added, and it is the mesh the rain model builds."""
    from studio.main_window import MainWindow
    out = []
    path = _copy_model(SEEP04, tmp)
    full = _load(path)
    dry = copy.deepcopy(full)
    dry["seepage_bc"]["specified_fluxes"] = []
    sig = MainWindow.mesh_signature
    out += _fail(sig(dry) == sig(full),
                 "SEEP-4: adding the three rain blocks changes the mesh signature, "
                 "though every vertex of them is already a zone vertex")
    m_dry, m_full = _build_mesh(dry, SEEP04_MESH), _build_mesh(full, SEEP04_MESH)
    out += _fail(_same_mesh(m_dry, m_full),
                 f"SEEP-4: the mesh built with the rain blocks differs from the one "
                 f"built without ({len(m_dry['nodes'])} vs {len(m_full['nodes'])} "
                 f"nodes), so keeping it would be wrong")

    # Through the editor: the window holds the dry model and its mesh; OK on Seep
    # BC with the three blocks keeps that mesh.
    mw.doc._dirty = False
    _quiet(mw.open_path, path)
    _set_bc(mw.doc.slope_data, dry["seepage_bc"])
    mw.doc.slope_data["mesh"] = m_dry
    mw.doc.results["mesh"] = m_dry
    _edit(mw, "seep_bc", lambda sd: _set_bc(sd, full["seepage_bc"]))
    out += _fail(mw.doc.slope_data.get("mesh") is m_dry,
                 "SEEP-4: OK on Seep BC with the three rain blocks dropped the mesh")
    return out


def test_mesh_dropped_for_a_new_point(mw, tmp):
    """A block whose end lands part way along a zone edge pins a new node there, and
    removing one that did unpins it: both change the mesh a build makes, so both
    drop the mesh. COMBO-2's set 2 head ends at (220, 110) on the shell face;
    clearing set 2 is the second case."""
    from studio.main_window import MainWindow
    out = []
    sig = MainWindow.mesh_signature
    full = _load(_copy_model(SEEP04, tmp))
    moved = copy.deepcopy(full)
    moved["seepage_bc"]["specified_fluxes"].append(
        {"flux": 1e-08, "coords": [(24.0, 12.0), (26.0, 12.0)]})   # mid-crest
    out += _fail(sig(moved) != sig(full),
                 "a flux block ending mid-crest at (26, 12) keeps the mesh signature")
    a, b = _build_mesh(full, SEEP04_MESH), _build_mesh(moved, SEEP04_MESH)
    out += _fail(not _same_mesh(a, b),
                 "a flux block ending mid-crest builds the same mesh; the signature "
                 "would be dropping a mesh that holds")

    path = _copy_model(COMBO02, tmp)
    with2 = _load(path)
    without2 = copy.deepcopy(with2)
    empty = {"specified_heads": [], "specified_fluxes": [], "exit_face": []}
    _set_bc(without2, with2["seepage_bc"], empty)
    out += _fail(sig(with2) != sig(without2),
                 "COMBO-2: clearing set 2 keeps the mesh signature")
    a, b = _build_mesh(with2, COMBO02_MESH), _build_mesh(without2, COMBO02_MESH)
    out += _fail(not _same_mesh(a, b),
                 "COMBO-2: the mesh built without set 2 is the one built with it")
    mw.doc._dirty = False
    _quiet(mw.open_path, path)
    mw.doc.slope_data["mesh"] = a
    mw.doc.results["mesh"] = a
    _edit(mw, "seep_bc", lambda sd: _set_bc(sd, sd["seepage_bc"], empty))
    out += _fail(mw.doc.slope_data.get("mesh") is None,
                 "COMBO-2: clearing set 2 kept a mesh with its pinned node")
    return out


# --------------------------------------------------------------------------- #
# The seepage solutions
# --------------------------------------------------------------------------- #
def _edit(mw, category, change):
    """OK on the ``category`` editor after ``change(slope_data)`` on top of what the
    editor writes back: the real edit path, with the dialog accepted as opened."""
    from PySide6.QtWidgets import QDialog
    from studio.main_window import CATEGORY_EDITORS
    editor = CATEGORY_EDITORS[category]
    real_apply, real_exec = editor.apply, QDialog.exec

    def apply(sd, dlg):
        real_apply(sd, dlg)
        change(sd)
    editor.apply = apply
    QDialog.exec = lambda self: 1
    try:
        _quiet(mw.edit_category, category)
    finally:
        editor.apply = real_apply
        QDialog.exec = real_exec


def _bump(key, by=1.0, mat=0):
    def change(sd):
        m = sd["materials"][mat]
        m[key] = float(m.get(key) or 0.0) + by
    return change


def _tseep(key, value):
    def change(sd):
        sd["tseep"] = dict(sd["tseep"], **{key: value})
    return change


def test_transient_solution(mw, tmp):
    out = []
    path = _copy_model(TRANSIENT, tmp)

    def opened():
        mw.doc._dirty = False
        _quiet(mw.open_path, path)
        bundle = mw.doc.results.get("transient_seep")
        if not bundle:
            raise RuntimeError("the transient earth dam opened with no transient "
                               "solution")
        return bundle

    def first_circle_depth(sd):
        sd["circles"][0]["Depth"] = float(sd["circles"][0].get("Depth") or 0.0) - 1.0

    def head_value(sd):
        sd["seepage_bc"]["specified_heads"][1]["head"] = 2.5

    def series_value(sd):
        ts = copy.deepcopy(sd["tseep"])
        name = next(iter(ts["series"]))
        ts["series"][name] = [None if v is None else float(v) + 0.5
                              for v in ts["series"][name]]
        sd["tseep"] = ts

    def rename(sd):
        sd["materials"][0]["name"] = "outer shell"

    keeps = [
        ("materials", "a cohesion", _bump("c")),
        ("materials", "a friction angle", _bump("phi")),
        ("materials", "a unit weight", _bump("gamma")),
        ("materials", "a material's name", rename),
        ("circles", "a circle's depth", first_circle_depth),
        ("global", "the horizontal seismic coefficient",
         lambda sd: sd.update(k_seismic=0.1)),
        ("transient", "the stability time", _tseep("stability_time", 7.0)),
    ]
    drops = [
        ("materials", "k1", _bump("k1", 0.5)),
        ("materials", "the specific storage", _bump("Ss", 1e-4)),
        ("materials", "h0", _bump("h0", -0.5)),
        ("global", "the unit weight of water",
         lambda sd: sd.update(gamma_water=float(sd["gamma_water"]) + 0.1)),
        ("seep_bc", "the tailwater head", head_value),
        ("transient", "the duration",
         lambda sd: _tseep("duration", float(sd["tseep"]["duration"]) * 2)(sd)),
        ("transient", "a series value", series_value),
    ]
    for category, what, change in keeps:
        bundle = opened()
        mesh = mw.doc.slope_data.get("mesh")
        _edit(mw, category, change)
        kept = mw.doc.results.get("transient_seep") is bundle
        out += _fail(kept and mw.transient_seep_view is not None
                     and mw.transient_solution() is not None,
                     f"an edit of {what} dropped the transient seepage solution "
                     f"(Run LEM then has no frame)")
        out += _fail(mw.doc.slope_data.get("mesh") is mesh,
                     f"an edit of {what} dropped the mesh")
        out += _fail(mw.doc.can_undo(), f"the edit of {what} made no undo step: the "
                                        f"test did not edit")
        if what == "a material's name" and kept:
            names = bundle["seep_data"].get("material_names") or []
            out += _fail(names[:1] == ["outer shell"],
                         f"the kept solution still names the shell {names[:1]}")
    for category, what, change in drops:
        opened()
        _edit(mw, category, change)
        out += _fail(mw.doc.results.get("transient_seep") is None
                     and mw.transient_seep_view is None,
                     f"an edit of {what} kept the transient seepage solution, "
                     f"which reads it")

    # COMBO-3's own step: the starter file opens with its transient solution, and
    # the reader types γ, γsat, c, φ and u = seep on both zones.
    start = _copy_model(TRANSIENT_START, tmp)
    mw.doc._dirty = False
    _quiet(mw.open_path, start)
    bundle = mw.doc.results.get("transient_seep")
    typed = {"shell": dict(gamma=20.0, gamma_sat=21.0, c=0.0, phi=32.0, u="seep"),
             "core": dict(gamma=19.0, gamma_sat=20.0, c=10.0, phi=25.0, u="seep")}
    _edit(mw, "materials", lambda sd: [m.update(typed[m["name"]])
                                       for m in sd["materials"]])
    out += _fail(bundle is not None and mw.doc.results.get("transient_seep") is bundle
                 and mw.transient_solution() is not None,
                 "COMBO-3: entering the strengths on the starter dropped its transient "
                 "seepage solution")
    out += _fail([m.get("u") for m in mw.doc.slope_data["materials"]] == ["seep"] * 2,
                 "COMBO-3: the strengths did not reach the model; the test did not edit")

    # An undo of a strength edit keeps it too; a new mesh drops it.
    bundle = opened()
    _edit(mw, "materials", _bump("c"))
    _quiet(mw._undo)
    out += _fail(mw.doc.results.get("transient_seep") is bundle,
                 "undoing a cohesion edit dropped the transient seepage solution")
    mesh = copy.copy(mw.doc.slope_data["mesh"])
    _quiet(mw._on_mesh_succeeded, mesh)
    out += _fail(mw.doc.results.get("transient_seep") is None,
                 "a new mesh kept the transient seepage solution solved on the old one")
    # Every other result still goes with any edit.
    opened()
    mw.doc.results["lem_solution"] = {"FS": 1.0}
    # A Monte Carlo result and its histogram tab, which an edit used to leave open
    # (only opening another file closed them).
    mw.doc.results["reliability"] = {"engine": "mc", "reliability": None}
    _quiet(mw._show_reliability_histogram)
    _edit(mw, "materials", _bump("c"))
    out += _fail("lem_solution" not in mw.doc.results,
                 "a cohesion edit kept the stability result")
    out += _fail("reliability" not in mw.doc.results
                 and mw.reliability_hist_canvas is None,
                 "a cohesion edit kept the Monte Carlo reliability result and its tab")
    return out


def test_steady_solution(mw, tmp):
    out = []
    path = _copy_model(STEADY, tmp)

    def opened():
        mw.doc._dirty = False
        _quiet(mw.open_path, path)
        bundle = (mw.doc.results.get("seep_solutions") or {}).get(1)
        if not bundle:
            raise RuntimeError("the solved Johnson dam opened with no steady solution")
        return bundle

    def add_set2(sd):
        # An exit face on set 2 over two zone vertices: no new mesh point, and (a
        # head on set 2 would be different) no water load derived from it either.
        sd["seepage_bc2"] = {"specified_heads": [], "specified_fluxes": [],
                             "exit_face": [(380.0, 180.0), (550.0, 100.0)]}
        sd["has_seepage_bc2"] = True

    for category, what, change, keep in (
            ("materials", "a cohesion", _bump("c"), True),
            ("seep_bc", "boundary set 2, which set 1 does not read", add_set2, True),
            ("materials", "k1", _bump("k1", 0.5), False),
            ("materials", "the van Genuchten exponent", _bump("vg_n", 0.1), False),
            ("seep_bc", "the reservoir head",
             lambda sd: sd["seepage_bc"]["specified_heads"][0].update(head=150.0),
             False)):
        bundle = opened()
        _edit(mw, category, change)
        got = (mw.doc.results.get("seep_solutions") or {}).get(1)
        tab = mw.seep_solution_canvas.get(1)
        if keep:
            out += _fail(got is bundle and tab is not None,
                         f"an edit of {what} dropped the steady seepage solution")
        else:
            out += _fail(got is None and tab is None,
                         f"an edit of {what} kept the steady seepage solution, "
                         f"which reads it")
    return out


def _field(mw, key="seep_u"):
    return mw.doc.slope_data.get(key)


def _equal(a, b):
    import numpy as np
    return (a is not None and b is not None
            and np.asarray(a).shape == np.asarray(b).shape
            and bool(np.array_equal(np.asarray(a, dtype=float),
                                    np.asarray(b, dtype=float))))


def _missing_field_reported(sd):
    from xslope.preflight import preflight
    rep = preflight(sd, "lem", {"surface": "circular", "search": True},
                    ids=["seep_field.missing"])
    return bool(rep.errors)


def test_stale_field_cleared(mw, tmp):
    """A dropped seepage solution takes its pore-pressure field out of the model.

    It used to stay: after a k1 edit on SEEP-2's solved dam the solution tab closed
    but ``seep_u`` kept the old field, the model checks passed it, a stability run
    read it, and Save rewrote ``_seep.csv`` from it, so the file reopened with the
    new k1 and the old field."""
    out = []
    sub = os.path.join(tmp, "stale")
    os.makedirs(sub, exist_ok=True)
    path = _copy_model(STEADY, sub)

    def opened(p=path):
        mw.doc._dirty = False
        _quiet(mw.open_path, p)
        return (mw.doc.results.get("seep_solutions") or {}).get(1)

    bundle = opened()
    out += _fail(bundle is not None and _equal(_field(mw), bundle["solution"]["u"]),
                 "SEEP-2: the solved dam opened without its solution's field in the "
                 "model")
    _edit(mw, "materials", _bump("c"))
    out += _fail(bundle is not None and _equal(_field(mw), bundle["solution"]["u"]),
                 "SEEP-2: a cohesion edit took the solution's field out of the model")
    _edit(mw, "materials", _bump("k1", 0.5))
    out += _fail(_field(mw) is None,
                 "SEEP-2: after a k1 edit dropped the solution, its field stayed in "
                 "the model")
    out += _fail(_missing_field_reported(mw.doc.slope_data),
                 "SEEP-2: after a k1 edit the model checks still pass u = seep on "
                 "the old field")
    csv = os.path.splitext(path)[0] + "_seep.csv"
    _quiet(mw.save)
    out += _fail(not os.path.exists(csv),
                 "SEEP-2: Save after the k1 edit rewrote _seep.csv from the old field")
    back = _load(path)
    out += _fail(back.get("seep_u") is None and _missing_field_reported(back),
                 "SEEP-2: the saved file reopens with a seepage field for a solution "
                 "it no longer has")
    out += _fail(opened() is None and _field(mw) is None,
                 "SEEP-2: reopened in Studio, the saved file restores a solution")

    # A field with no Studio solve behind it stays: the saved field as read at
    # open, with no solution restored beside it (as when the restore fails, or a
    # field imported from SEEP/W).
    sub2 = os.path.join(tmp, "stale_file")
    os.makedirs(sub2, exist_ok=True)
    path2 = _copy_model(STEADY, sub2)
    opened(path2)
    mw.doc.results.pop("seep_solutions", None)
    field = copy.deepcopy(_field(mw))
    _edit(mw, "materials", _bump("k1", 0.5))
    out += _fail(field is not None and _equal(_field(mw), field),
                 "a field with no solution behind it was cleared by a k1 edit")

    # A transient frame put into the model for a run goes with the transient
    # solution, and its instant with it.
    sub3 = os.path.join(tmp, "stale_t")
    os.makedirs(sub3, exist_ok=True)
    tpath = _copy_model(TRANSIENT, sub3)
    for what, change, stays in (("a cohesion", _bump("c"), True),
                                ("k1", _bump("k1", 0.5), False)):
        mw.doc._dirty = False
        _quiet(mw.open_path, tpath)
        ok = _quiet(mw._apply_transient_analysis_frame)
        frame = copy.deepcopy(_field(mw))
        t = mw.doc.slope_data.get("seep_u_time")
        if not ok or frame is None or t is None:
            out.append("COMBO-3: no transient frame reached the model for a run")
            continue
        _edit(mw, "materials", change)
        if stays:
            out += _fail(_equal(_field(mw), frame)
                         and mw.doc.slope_data.get("seep_u_time") == t,
                         f"COMBO-3: {what} edit took the transient frame's field "
                         f"out of the model")
        else:
            out += _fail(_field(mw) is None
                         and mw.doc.slope_data.get("seep_u_time") is None,
                         f"COMBO-3: after a {what} edit dropped the transient "
                         f"solution, its frame (t = {t:g}) stayed in the model")
    return out


FEM_REFUSAL = "carries no pore-pressure field on its mesh"


def _run_fem_dialog(mw):
    """Run → Run FEM… as the menu opens it; returns the dialog, closed unrun."""
    from PySide6.QtWidgets import QDialog
    seen = []
    real = QDialog.exec

    def exec_(self):
        seen.append(self)
        return 0                                   # Cancel: nothing is run
    QDialog.exec = exec_
    try:
        _quiet(mw.run_fem)
    finally:
        QDialog.exec = real
    return seen[0] if seen else None


def test_fem_refused_without_its_field(mw, tmp):
    """A finite element run with u = seep and no field on its mesh is refused
    (602a). After a k1 edit drops SEEP-2's solution and its field, Run FEM's model
    checks used to pass, and the run went on at zero pore pressure in every
    element; the Run button now stays off, with the sentence saying why. With the
    solution loaded the run is not refused for it."""
    out = []
    sub = os.path.join(tmp, "fem_refused")
    os.makedirs(sub, exist_ok=True)
    path = _copy_model(STEADY, sub)
    for what, change, refused in (("as opened", None, False),
                                  ("after k1 x10", _bump("k1", 9.0), True)):
        mw.doc._dirty = False
        _quiet(mw.open_path, path)
        if change is not None:
            _edit(mw, "materials", change)
        dlg = _run_fem_dialog(mw)
        if dlg is None:
            out.append(f"SEEP-2 {what}: Run FEM opened no dialog")
            continue
        rep = dlg.preflight.report
        mine = [f for f in (rep.findings if rep is not None else ())
                if f.rule_id == "seep_field.missing"]
        if refused:
            out += _fail(mine and FEM_REFUSAL in mine[0].message
                         and dlg.preflight.blocked and not dlg._ok.isEnabled(),
                         f"SEEP-2 {what}: Run FEM is not refused for the missing "
                         f"field (field in model: {_field(mw) is not None}); the run "
                         f"would go on at zero pore pressure")
        else:
            out += _fail(not mine, f"SEEP-2 {what}: Run FEM is refused for a field "
                                   f"the model carries")
        out += _fail(mw._fem_runner is None, f"SEEP-2 {what}: an FEM run started")
        dlg.deleteLater()

    # A transient model carries no stored field: the run stages a frame of its
    # loaded solution before the solver starts, and the dialog says so. COMBO-3's
    # dam, opened with its solution, is not refused for the field.
    sub = os.path.join(tmp, "fem_transient")
    os.makedirs(sub, exist_ok=True)
    mw.doc._dirty = False
    _quiet(mw.open_path, _copy_model(TRANSIENT, sub))
    dlg = _run_fem_dialog(mw)
    rep = dlg.preflight.report if dlg is not None else None
    mine = [f for f in (rep.findings if rep is not None else ())
            if f.rule_id == "seep_field.missing"]
    out += _fail(dlg is not None and mw.transient_solution() is not None and not mine,
                 "COMBO-3's dam with its transient solution loaded: Run FEM is refused "
                 "for a field the run stages from the solution")
    if dlg is not None:
        dlg.deleteLater()
    return out


def test_kept_steady_field_restored(mw, tmp):
    """A steady solution that stays puts its own field back where a dropped
    transient frame had replaced it."""
    import numpy as np
    out = []
    sub = os.path.join(tmp, "restored")
    os.makedirs(sub, exist_ok=True)
    path = _copy_model(TRANSIENT, sub)
    mw.doc._dirty = False
    _quiet(mw.open_path, path)
    tb = mw.doc.results["transient_seep"]
    # A steady set-1 solution on the same mesh and inputs (its values made
    # distinct from every frame), applied as a solve applies it.
    u0 = np.asarray(tb["transient"]["frames"][0]["u"], dtype=float) + 1.0
    steady = {"solution": {"u": u0}, "options": {"bc": 1},
              "solved_from": mw._seepage_record(1)}
    mw.doc.results["seep_solutions"] = {1: steady}
    _quiet(mw._apply_seep_field, steady, 1)
    _quiet(mw._apply_transient_analysis_frame)
    out += _fail(not _equal(_field(mw), u0),
                 "the test could not stage a transient frame over the steady field")
    _edit(mw, "materials", _bump("Ss", 1e-4))          # the march reads it; set 1 not
    out += _fail((mw.doc.results.get("seep_solutions") or {}).get(1) is steady
                 and mw.doc.results.get("transient_seep") is None,
                 "a storage edit did not keep the steady solution and drop the "
                 "transient one")
    out += _fail(_equal(_field(mw), u0),
                 "with the transient solution dropped, the kept steady solution's "
                 "field is not the model's")
    return out


CHECKS = [
    ("blocks on existing vertices keep the mesh", test_mesh_kept_for_blocks_on_its_vertices),
    ("a newly pinned point drops the mesh", test_mesh_dropped_for_a_new_point),
    ("a transient solution survives what it does not read",
     test_transient_solution),
    ("a steady solution survives what it does not read", test_steady_solution),
    ("a dropped solution takes its field with it", test_stale_field_cleared),
    ("a kept steady solution keeps its field", test_kept_steady_field_restored),
    ("Run FEM is refused without its field", test_fem_refused_without_its_field),
]


@studio_settings.isolated
def run():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    QApplication.instance() or QApplication([])
    from studio.main_window import MainWindow

    mw = MainWindow()
    fails = []
    tmp = tempfile.mkdtemp(prefix="seep_survives_")
    try:
        for name, fn in CHECKS:
            try:
                got = fn(mw, tmp)
            except Exception as e:                    # noqa: BLE001
                import traceback
                traceback.print_exc()
                got = [f"{name}: raised {type(e).__name__}: {e}"]
            print(f"  {name:<54s} {'ok' if not got else 'FAIL'}")
            fails += got
    finally:
        mw.doc._dirty = False
        mw.close()
        shutil.rmtree(tmp, ignore_errors=True)
    return fails


def main():
    print("Edits keep the mesh and seepage solutions they cannot change:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll invalidation checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
