"""Standing checks on the engine fixes the input-graying audit led to (522-526,
the FEM's blank strength option, and the second piezometric line).

Each check is held to a number worked out here by hand, or to a contrasting model
whose answer is already known, and each one fails on the engine as it stood
before the fix:

  1. 522 FEM pore pressure is per material. Every element takes its pore pressure
     from its OWN material's u, as the limit equilibrium slicer takes each slice
     base's: none is zero, piezo is gamma_w times the depth below the piezometric
     line, ru is that material's own ru times the soil column. A row no zone uses,
     and an elastic row, switch nothing for anyone else; mixed sources run.
  2. 523 No suction strength for cp and ratio. phi_b adds the apparent cohesion
     tan(phi_b) s only on an effective-stress envelope (mc, pow, hb), in both
     engines.
  3. 524 Seepage ignores joints. A mesh split along a joint line seeps exactly as
     the unsplit mesh does, the field is defined on both faces of the joint, and
     the finite element run reading it finds pore pressure on both faces.
  4. 525 A Size on an SSR zone or a joints polygon refines the mesh inside it, as
     on a refine polygon.
  5. 526 A joint's tension cutoff follows the Tension SRF setting: divided by the
     trial factor when it is on, held when it is off.
  6. The FEM refuses a blank strength option on a material the mesh uses.
  7. The second piezometric line shapes and refuses only a rapid drawdown run; a
     single-stage factor of safety is the same with or without it.
  8. The report's finite element section names each material's pore-pressure
     source: "at every node" only where every meshed material shares one.
  9. The Taylor-series reliability takes a slice count (Studio's Reliability and
     Variance Pareto dialogs offer one, and the run used to ignore it and always
     take 40): none given is exactly 40 on a specified surface, and 60 slices is a
     different answer.

Every model is small and nothing is strength-reduced: the seepage check makes two
unconfined solves of one earth dam, and the joint check runs two viscoplastic
iterations. Needs gmsh for the meshes; skips the mesh checks without it.
"""
import contextlib
import copy
import io
import math
import os
import sys
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np

LFACE = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_lface.xlsx")
EARTH_DAM = os.path.join(_REPO, "docs", "inputs", "seep", "xslope_earth_dam1.xlsx")
ROCK_JOINTS = os.path.join(_REPO, "docs", "tutorials", "files", "xslope_rock_joints.xlsx")

_CACHE = {}


@contextlib.contextmanager
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            yield


def _gmsh():
    try:
        import gmsh  # noqa: F401
        return True
    except Exception:
        return False


def _load(path):
    from xslope.fileio import load_slope_data
    if path not in _CACHE:
        with _quiet():
            _CACHE[path] = load_slope_data(path)
    return copy.deepcopy(_CACHE[path])


def _lface_mesh():
    """lface (three mc materials on u = piezo) on a coarse tri6 mesh."""
    if "lface_mesh" not in _CACHE:
        from xslope.mesh import get_material_polygons, build_mesh_from_polygons
        sd = _load(LFACE)
        with _quiet():
            polys = get_material_polygons(sd)
            xs = [x for x, _ in sd["ground_surface"].coords]
            _CACHE["lface_mesh"] = build_mesh_from_polygons(
                polys, target_size=(max(xs) - min(xs)) / 25, element_type="tri6")
    return _CACHE["lface_mesh"]


def _lface_fem_model():
    sd = _load(LFACE)
    sd["piezo_phreatic"] = False          # the hand value below is plain hydrostatic
    for m in sd["materials"]:
        m.update(E=1e5, nu=0.3)
    return sd


def _fem_prepared(sd, mesh):
    from xslope.fem import build_fem_data, _prepare_fem_model
    with _quiet():
        fd = build_fem_data(sd, mesh)
        prep = _prepare_fem_model(fd)
    return fd, prep


def _gauss_points(fd, prep, e):
    """(x, y) of every Gauss point of element e, from its shape functions."""
    idx = fd["elements"][e][:fd["element_types"][e]]
    xy = np.asarray(fd["nodes"])[idx]
    return [(float(g["N"] @ xy[:, 0]), float(g["N"] @ xy[:, 1]))
            for g in prep["elem_gp_data"][e]]


def _hydrostatic(sd, x, y):
    """gamma_w times the depth below the piezometric line, worked out here."""
    pl = sorted(sd["piezo_line"])
    yp = float(np.interp(x, [p[0] for p in pl], [p[1] for p in pl]))
    return float(sd["gamma_water"]) * max(0.0, yp - y)


# ---------------------------------------------------------------------------
# 1. 522: the FEM honors each material's own u
# ---------------------------------------------------------------------------

def check_fem_pore_pressure_per_material():
    if not _gmsh():
        return []
    fails = []
    mesh = _lface_mesh()
    em = np.asarray(mesh["element_materials"], dtype=int)

    # (a) Silt u = none beside Sand and Clay on piezo: Silt is dry, the others
    # hydrostatic to the line. The pooled source put the line's pressure into Silt.
    sd = _lface_fem_model()
    sd["materials"][0]["u"] = "none"
    fd, prep = _fem_prepared(sd, mesh)
    worst_dry, worst_wet = 0.0, 0.0
    for e in range(len(em)):
        for (x, y), u in zip(_gauss_points(fd, prep, e), prep["u_gp"][e]):
            if em[e] == 1:
                worst_dry = max(worst_dry, abs(u))
            else:
                worst_wet = max(worst_wet, abs(u - _hydrostatic(sd, x, y)))
    if worst_dry != 0.0:
        fails.append(f"522: a u = none material carries pore pressure in the FEM "
                     f"(max {worst_dry:.4g} at a Gauss point; want 0)")
    if worst_wet > 1e-6:
        fails.append(f"522: the piezo materials beside a dry one depart from "
                     f"gamma_w * depth by up to {worst_wet:.4g}")

    # (b) Every zone dry, a row no zone uses on piezo, and Clay elastic on piezo:
    # the soils stay dry; the elastic elements take their own row's line.
    sd = _lface_fem_model()
    for m in sd["materials"]:
        m["u"] = "none"
    sd["materials"][2].update(option="elastic", u="piezo")
    sd["materials"].append(dict(sd["materials"][0], name="Unused", u="piezo"))
    fd, prep = _fem_prepared(sd, mesh)
    soil = max((abs(u) for e in range(len(em)) if em[e] != 3 for u in prep["u_gp"][e]),
               default=0.0)
    el_err = max((abs(u - _hydrostatic(sd, x, y))
                  for e in range(len(em)) if em[e] == 3
                  for (x, y), u in zip(_gauss_points(fd, prep, e), prep["u_gp"][e])),
                 default=0.0)
    if soil != 0.0:
        fd_src = fd.get("pp_option")
        fails.append(f"522: an unused row or an elastic row on u = piezo put pore "
                     f"pressure into u = none soils (max {soil:.4g}; source "
                     f"{fd_src!r})")
    if el_err > 1e-6:
        fails.append(f"522: the elastic row's own elements depart from its own "
                     f"piezometric line by {el_err:.4g}")

    # (c) Mixed sources run: Silt u = none carrying a stale ru, Sand u = ru, Clay
    # on piezo. Sand's pore pressure is the one the all-ru model gives it (its own
    # ru on its own column, a field the single-source path always computed), Silt
    # stays dry, Clay stays hydrostatic.
    sd = _lface_fem_model()
    sd["materials"][0].update(u="none", ru=0.4)
    sd["materials"][1].update(u="ru", ru=0.2)
    try:
        fd, prep = _fem_prepared(sd, mesh)
    except ValueError as exc:
        fails.append(f"522: mixed none / ru / piezo sources refused: {exc}")
    else:
        ref = _lface_fem_model()
        for m in ref["materials"]:
            m.update(u="ru", ru=0.2)
        _, prep_ref = _fem_prepared(ref, mesh)
        d_sand = max(abs(a - b) for e in range(len(em)) if em[e] == 2
                     for a, b in zip(prep["u_gp"][e], prep_ref["u_gp"][e]))
        silt = max(abs(u) for e in range(len(em)) if em[e] == 1 for u in prep["u_gp"][e])
        clay = max(abs(u - _hydrostatic(sd, x, y)) for e in range(len(em)) if em[e] == 3
                   for (x, y), u in zip(_gauss_points(fd, prep, e), prep["u_gp"][e]))
        if d_sand > 1e-9 * max(1.0, max(max(g) for g in prep_ref["u_gp"])):
            fails.append(f"522: the ru material's pore pressure differs from the "
                         f"all-ru model's by {d_sand:.4g}")
        if silt != 0.0:
            fails.append(f"522: a u = none row's stale ru was applied (max u {silt:.4g})")
        if clay > 1e-6:
            fails.append(f"522: the piezo material in a mixed model departs from "
                         f"gamma_w * depth by {clay:.4g}")
    return fails


# ---------------------------------------------------------------------------
# 2. 523: no suction strength for cp and ratio
# ---------------------------------------------------------------------------

def _lface_slices(sd):
    from xslope.slice import generate_slices
    with _quiet():
        ok, out = generate_slices(sd, circle=sd["circles"][0], num_slices=30,
                                  debug=False)
    if not ok:
        raise RuntimeError(out)
    return out[0]


def check_suction_effective_stress_only():
    fails = []
    tan25 = math.tan(math.radians(25.0))
    # LEM, cp: every row undrained on u = piezo, phi_b = 25. No slice gains a
    # suction cohesion, and the factor of safety is the one with phi_b blank.
    from xslope.solve import bishop
    base = _load(LFACE)
    for m in base["materials"]:
        m.update(option="cp", c=400.0, cp=0.0, phi=0.0, u="piezo",
                 phi_b=None, s_cap=None)
    with_pb = copy.deepcopy(base)
    for m in with_pb["materials"]:
        m["phi_b"] = 25.0
    df0, df1 = _lface_slices(base), _lface_slices(with_pb)
    if float(np.abs(df1["c_suction"]).max()) != 0.0:
        fails.append(f"523: cp slices carry suction cohesion up to "
                     f"{float(df1['c_suction'].max()):.4g} (want 0)")
    with _quiet():
        fs0, fs1 = bishop(df0)[1]["FS"], bishop(df1)[1]["FS"]
    if fs0 != fs1:
        fails.append(f"523: phi_b on cp rows moved Bishop's FS {fs0:.6f} -> {fs1:.6f}")
    # LEM, ratio: the same.
    rat = copy.deepcopy(with_pb)
    for m in rat["materials"]:
        m.update(option="ratio", c=50.0, cp=0.25)
    dfr = _lface_slices(rat)
    if float(np.abs(dfr["c_suction"]).max()) != 0.0:
        fails.append("523: ratio slices carry suction cohesion (want 0)")
    # LEM, mc contrast: the credit is still there on an effective-stress
    # material, and it is tan(phi_b) times the height above the line, by hand.
    mc = _load(LFACE)
    mc["piezo_phreatic"] = False
    for m in mc["materials"]:
        m.update(u="piezo", phi_b=25.0, s_cap=None)
    dfm = _lface_slices(mc)
    s_hand = (np.maximum(0.0, dfm["y_u"] - dfm["piezo_y"]) * float(mc["gamma_water"])
              * tan25)
    err = float(np.max(np.abs(dfm["c_suction"] - s_hand)))
    if not (float(s_hand.max()) > 0.0 and err < 1e-6 * max(1.0, float(s_hand.max()))):
        fails.append(f"523: the mc contrast lost its suction cohesion (hand max "
                     f"{float(s_hand.max()):.4g}, error {err:.3g})")

    # FEM: Silt cp, Sand ratio, Clay mc, all on piezo with phi_b = 25. Only the mc
    # elements carry tan(phi_b), whether phi_b comes from the file or a kwarg.
    if _gmsh():
        from xslope.fem import build_fem_data, _resolve_suction_by_elem
        mesh = _lface_mesh()
        sd = _lface_fem_model()
        sd["materials"][0].update(option="cp", c=400.0, cp=0.0)
        sd["materials"][1].update(option="ratio", c=50.0, cp=0.25)
        for m in sd["materials"]:
            m.update(u="piezo", phi_b=25.0)
        with _quiet():
            fd = build_fem_data(sd, mesh)
        em = np.asarray(fd["element_materials"], dtype=int)
        for label, kw in (("file", None),
                          ("kwarg", {m["name"]: 25.0 for m in sd["materials"]})):
            tpb, _, active = _resolve_suction_by_elem(fd, kw, None, em)
            want = np.where(em == 3, tan25, 0.0)
            if not active or float(np.max(np.abs(tpb - want))) > 1e-12:
                fails.append(f"523: FEM tan(phi_b) per element ({label}) is not "
                             f"tan 25 on the mc material and 0 on cp / ratio: "
                             f"cp {sorted(set(np.round(tpb[em == 1], 6)))}, ratio "
                             f"{sorted(set(np.round(tpb[em == 2], 6)))}")
    return fails


# ---------------------------------------------------------------------------
# 3. 524: seepage ignores joints
# ---------------------------------------------------------------------------

def check_seepage_through_joint():
    if not _gmsh():
        return []
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_joint_options)
    from xslope.seep import build_seep_data, run_seepage_analysis
    from xslope.fem import build_fem_data
    fails = []
    sd = _load(EARTH_DAM)
    sd["joint_lines"] = [{"label": "J1", "x1": 80.0, "y1": 0.5, "x2": 80.0, "y2": 12.0,
                          "c": 0.0, "phi": 30.0, "c_res": float("nan"),
                          "phi_res": float("nan"), "dil": float("nan"), "t_cut": 0.0,
                          "kn": float("nan"), "ks": float("nan"), "jred": ""}]
    lines, _, _ = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    out = {}
    for tag, jo in (("unsplit", None), ("split", extract_joint_options(sd))):
        with _quiet():
            mesh = build_mesh_from_polygons(polys, 3.0, "tri6", lines=lines,
                                            joint_lines=jo)
            sol = run_seepage_analysis(build_seep_data(mesh, sd, seep_bc=1))
        out[tag] = (mesh, sol)
    m_s, s_s = out["split"]
    q_u, q_s = out["unsplit"][1]["flowrate"], s_s["flowrate"]
    if len(m_s.get("elements_joint", [])) == 0:
        return ["524: the joint did not split the mesh (no joint elements)"]
    # The contrasting model is the same mesh unsplit: the flow through the dam
    # must not see the joint.
    if abs(q_s - q_u) > 1e-6 * abs(q_u):
        fails.append(f"524: the flow through the dam is {q_s:.6g} with the joint "
                     f"and {q_u:.6g} without (a joint is not a no-flow wall)")
    nodes = np.asarray(m_s["nodes"])
    _, inv, cnt = np.unique(nodes, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    copies = np.flatnonzero(cnt[inv] > 1)
    spread = max(float(np.ptp(s_s["head"][inv == g])) for g in set(inv[copies]))
    if spread > 1e-9:
        fails.append(f"524: the two faces of the joint carry heads {spread:.4g} apart")
    # The finite element run reading that field: pore pressure on both faces.
    fsd = copy.deepcopy(sd)
    for m in fsd["materials"]:
        m.update(option="mc", c=5.0, phi=30.0, gamma=120.0, E=1e5, nu=0.3, u="seep")
    fsd["seep_u"] = np.asarray(s_s["u"], dtype=float)
    with _quiet():
        fd = build_fem_data(fsd, m_s)
    us = np.asarray(fd["u_seep"] if fd.get("u_seep") is not None else fd["u"])
    wet = [g for g in set(inv[copies]) if np.all(us[inv == g] > 0.0)]
    gap = max(float(np.ptp(us[inv == g])) for g in set(inv[copies]))
    if not wet or gap > 1e-9:
        fails.append(f"524: the FEM's seepage pore pressure on the joint's two faces: "
                     f"{len(wet)} wet stations, largest difference {gap:.4g}")
    return fails


# ---------------------------------------------------------------------------
# 4. 525: a Size on an SSR zone or a joints polygon reaches the mesher
# ---------------------------------------------------------------------------

def check_size_on_overlays():
    from xslope.mesh import extract_size_regions
    fails = []
    sd = _load(LFACE)
    xs = [x for x, _ in sd["ground_surface"].coords]
    ys = [y for _, y in sd["ground_surface"].coords]
    x0, x1 = min(xs), max(xs)
    target = (x1 - x0) / 15.0
    small = target / 4.0
    # A box well inside the section, below the lowest ground point.
    yb = min(ys)
    box = [(x0 + 0.40 * (x1 - x0), yb - 3.0 * target), (x0 + 0.60 * (x1 - x0), yb - 3.0 * target),
           (x0 + 0.60 * (x1 - x0), yb - 1.0 * target), (x0 + 0.40 * (x1 - x0), yb - 1.0 * target)]
    sd["ssr_zones"] = [{"kind": "reduce", "polygon": box, "label": "SSR reduce",
                        "size": small}]
    sd["joint_zones"] = [{"polygon": box, "label": "set A", "size": small * 1.5,
                          "mat_id": None}]
    regs = extract_size_regions(sd)
    sizes = sorted(r["size"] for r in regs)
    if sizes != sorted([small, small * 1.5]):
        fails.append(f"525: extract_size_regions returned sizes {sizes}; want the SSR "
                     f"zone's {small:g} and the joints polygon's {small * 1.5:g}")
    if not _gmsh():
        return fails
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from shapely.geometry import Polygon, Point
    ring = Polygon(box)

    def _inside_edge(size_regions):
        with _quiet():
            mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size=target,
                                            element_type="tri3", size_regions=size_regions)
        nd = np.asarray(mesh["nodes"])
        lens = []
        for e, t in zip(mesh["elements"], mesh["element_types"]):
            idx = list(e[:t])
            c = nd[idx].mean(axis=0)
            if ring.contains(Point(c)):
                for a, b in zip(idx, idx[1:] + idx[:1]):
                    lens.append(float(np.hypot(*(nd[a] - nd[b]))))
        return float(np.mean(lens)) if lens else float("nan")

    plain = _inside_edge(None)
    sized = _inside_edge(regs)
    if not (sized <= 1.5 * small and plain > 2.0 * small):
        fails.append(f"525: inside the sized SSR zone the mean edge is {sized:.3g} "
                     f"(asked {small:.3g}; {plain:.3g} without the Size)")
    return fails


# ---------------------------------------------------------------------------
# 5. 526: a joint's tension cutoff follows Tension SRF
# ---------------------------------------------------------------------------

def _one_joint(tcut, jred=True):
    """One straight joint element, its faces 1 apart in x, normal +y: three node
    pairs, the upper face's dofs 0..5 and the lower face's 6..11."""
    n = 1
    return {
        "n": n, "dof": np.arange(12).reshape(1, 12), "w": np.ones((n, 3)),
        "tx": np.ones(n), "ty": np.zeros(n), "nx": np.zeros(n), "ny": np.ones(n),
        "kn": np.full(n, 100.0), "ks": np.full(n, 100.0),
        "cj": np.zeros(n), "tanphi": np.full(n, 0.5), "tcut": np.full(n, tcut),
        "jred": np.full(n, jred), "tip": np.zeros((n, 3), dtype=bool),
    }


def check_joint_tension_cutoff():
    fails = []
    try:
        from xslope.joint import joint_reduced_tcut, joint_state
    except ImportError as exc:
        return [f"526: {exc}"]
    jd = _one_joint(4.0)
    if not np.allclose(joint_reduced_tcut(jd, 2.0, True), [2.0]):
        fails.append("526: t_cut is not divided by F = 2 with Tension SRF on")
    if not np.allclose(joint_reduced_tcut(jd, 2.0, False), [4.0]):
        fails.append("526: t_cut is reduced with Tension SRF off")
    if not np.allclose(joint_reduced_tcut(_one_joint(4.0, jred=False), 2.0, True), [4.0]):
        fails.append("526: a joint that opts out of the reduction had its t_cut reduced")
    # Pull the faces 0.03 apart: t_n = -k_n * 0.03 = -3. With t_cut = 4 the joint
    # holds at full strength; at F = 2 (cutoff 2) it opens; with Tension SRF off
    # (cutoff 4) it holds.
    u = np.zeros(12)
    u[1:6:2] = 0.03                         # upper face moves +y (away)
    cj, tp = jd["cj"], jd["tanphi"]
    held = joint_state(jd, u, cj, tp)["open"]
    opened = joint_state(jd, u, cj, tp, tcut_r=joint_reduced_tcut(jd, 2.0, True))["open"]
    off = joint_state(jd, u, cj, tp, tcut_r=joint_reduced_tcut(jd, 2.0, False))["open"]
    if held.any() or not opened.all() or off.any():
        fails.append(f"526: t_n = -3 against t_cut 4: open at F = 1 {held.any()}, at "
                     f"F = 2 with SRF {opened.all()}, at F = 2 without {off.any()} "
                     f"(want False, True, False)")

    # The finite element solve hands the reduced cutoff to the interface sweep.
    if not (_gmsh() and os.path.exists(ROCK_JOINTS)):
        return fails
    import xslope.fem as fem
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_joint_options)
    sd = _load(ROCK_JOINTS)
    for j in sd["joint_lines"]:
        j["t_cut"] = 3.0
    lines, _, _ = extract_constraint_line_geometry(sd)
    xs = [x for x, _ in sd["ground_surface"].coords]
    with _quiet():
        mesh = build_mesh_from_polygons(get_material_polygons(sd, reinf_lines=lines),
                                        (max(xs) - min(xs)) / 12, "tri6", lines=lines,
                                        joint_lines=extract_joint_options(sd))
        fd = fem.build_fem_data(sd, mesh)
    jd_m = fd.get("joint_data")
    if jd_m is None or jd_m["n"] == 0:
        return fails + ["526: the rock-joints model built no joint elements"]
    seen = []
    real = fem.joint_vp_sweep

    def spy(*a, **kw):
        seen.append(kw.get("tcut_r"))
        return real(*a, **kw)

    want_on = np.where(jd_m["jred"], 3.0 / 2.0, 3.0)
    for srf, want in ((True, want_on), (False, np.full(jd_m["n"], 3.0))):
        seen.clear()
        fem.joint_vp_sweep = spy
        try:
            with _quiet():
                fem.solve_fem(fd, F=2.0, max_iterations=2, tension_srf=srf,
                              fem_solver="viscoplastic", early_exit=False)
        except Exception:
            pass                            # two iterations decide nothing
        finally:
            fem.joint_vp_sweep = real
        got = seen[0] if seen else None
        if got is None or not np.allclose(got, want):
            fails.append(f"526: solve_fem at F = 2, Tension SRF {srf}: the interface "
                         f"sweep read t_cut {None if got is None else np.unique(got)}, "
                         f"want {np.unique(want)}")
    groups = fem._nr_build_joints(fd, 2.0, tension_srf=True)
    if not np.allclose(groups[0].get("tcut_r", np.nan), want_on):
        fails.append("526: the Newton driver's joint group does not carry t_cut / F")
    return fails


# ---------------------------------------------------------------------------
# 6. The FEM refuses a blank strength option on a meshed material
# ---------------------------------------------------------------------------

def check_fem_blank_option():
    if not _gmsh():
        return []
    from xslope.fem import build_fem_data
    fails = []
    mesh = _lface_mesh()
    sd = _lface_fem_model()
    sd["materials"][0].update(option="", c=123.0, phi=21.0)
    try:
        with _quiet():
            build_fem_data(sd, mesh)
    except ValueError as exc:
        if "Silt" not in str(exc):
            fails.append(f"blank option: refused without naming the material: {exc}")
    else:
        fails.append("blank option: a meshed material with no strength option ran "
                     "(as Mohr-Coulomb on c = 123, phi = 21)")
    # A blank row no zone uses is fine.
    sd = _lface_fem_model()
    sd["materials"].append(dict(sd["materials"][0], name="Seep only", option=""))
    try:
        with _quiet():
            build_fem_data(sd, mesh)
    except ValueError as exc:
        fails.append(f"blank option: an unused blank row was refused: {exc}")
    return fails


# ---------------------------------------------------------------------------
# 7. The second piezometric line in a single-stage run
# ---------------------------------------------------------------------------

def check_piezo_line2_single_stage():
    from xslope.slice import generate_slices
    from xslope.solve import bishop
    fails = []
    base = _load(LFACE)
    circle = base["circles"][0]
    pl = sorted(base["piezo_line"])

    def fs(sd, **kw):
        with _quiet():
            ok, out = generate_slices(sd, circle=circle, num_slices=30, debug=False, **kw)
            if not ok:
                raise RuntimeError(out)
            return bishop(out[0])[1]["FS"], out[0]

    fs0, df0 = fs(base)
    short = copy.deepcopy(base)
    short["piezo_line2"] = [(pl[0][0], pl[0][1]), (pl[0][0] + 1.0, pl[0][1])]
    full = copy.deepcopy(base)
    full["piezo_line2"] = [(pl[0][0], pl[0][1] - 15.0), (pl[-1][0], pl[-1][1] - 15.0)]
    for label, sd in (("a short", short), ("a full-width", full)):
        try:
            f, df = fs(sd)
        except Exception as exc:
            fails.append(f"line 2: {label} second line stopped a single-stage run: {exc}")
            continue
        if f != fs0 or len(df) != len(df0):
            fails.append(f"line 2: {label} second line moved a single-stage run: FS "
                         f"{fs0:.6f} -> {f:.6f}, {len(df0)} -> {len(df)} slices")
    # Rapid drawdown still reads it: a short line refuses, a full one shapes the
    # slices (a boundary where the surface crosses it).
    try:
        fs(short, rapid=True)
    except ValueError as exc:
        if "second piezometric line" not in str(exc):
            fails.append(f"line 2: a rapid run refused a short line for another "
                         f"reason: {exc}")
    except TypeError as exc:
        fails.append(f"line 2: generate_slices has no rapid switch: {exc}")
    else:
        fails.append("line 2: a rapid run read past the end of a short second line")
    return fails


# ---------------------------------------------------------------------------
# 8. The report names each material's pore-pressure source
# ---------------------------------------------------------------------------

SHEETS_MAT = os.path.join(_REPO, "docs", "usage", "sample_sheets", "sheets_mat.xlsx")


def check_report_pore_basis():
    from xslope.report import _fem_pore_basis, resolve_options
    from xslope.fem import _normalize_pp_option
    fails = []
    opts = resolve_options({})

    def said(fd):
        block = _fem_pore_basis(fd, {}, opts)
        return block.text if block is not None else ""

    # sheets_mat.xlsx, every row meshed: two rows on the piezometric line, one on
    # ru, three on none (the sheet carries no geometry, so the record is the one
    # build_fem_data writes: each row's source, its name, the meshed rows).
    m = _load(SHEETS_MAT)["materials"]
    text = said({"pp_option_by_mat": [_normalize_pp_option(x.get("u")) for x in m],
                 "material_names": [x["name"] for x in m],
                 "element_materials": np.arange(1, len(m) + 1)})
    for part in ("in Embankment and Soft Clay is taken from the piezometric line",
                 "in Drain Sand from a pore pressure ratio r_u",
                 "Rockfill, Weathered Rock and Bedrock carry none"):
        if part not in text:
            fails.append(f"report: the mixed sheets_mat model is described as "
                         f"{text!r}, which does not say {part!r}")
    if "every node" in text:
        fails.append(f"report: a mixed model is said to take pore pressure at "
                     f"every node: {text!r}")
    if not _gmsh():
        return fails
    from xslope.fem import build_fem_data
    mesh = _lface_mesh()
    # lface as shipped: one source, the sentence it always had.
    sd = _lface_fem_model()
    with _quiet():
        fd = build_fem_data(sd, mesh)
    want = "Pore pressure at every node is taken from the piezometric line."
    if said(fd) != want:
        fails.append(f"report: the single-source model now reads {said(fd)!r}; "
                     f"want {want!r}")
    # Silt dry beside the piezo rows.
    sd["materials"][0]["u"] = "none"
    with _quiet():
        fd = build_fem_data(sd, mesh)
    want = ("Pore pressure in Sand and Clay is taken from the piezometric line; "
            "Silt carries none.")
    if said(fd) != want:
        fails.append(f"report: a dry Silt beside piezo rows reads {said(fd)!r}; "
                     f"want {want!r}")
    return fails


# ---------------------------------------------------------------------------
# 9. The Taylor-series reliability takes a slice count
# ---------------------------------------------------------------------------

RELIABILITY = os.path.join(_REPO, "docs", "lem", "files", "xslope_reliability.xlsx")


def check_taylor_num_slices():
    from xslope.reliability import reliability
    from xslope.sensitivity import variance_contribution
    fails = []
    got = {}
    for n in (None, 40, 60):
        sd = _load(RELIABILITY)
        with _quiet():
            ok, res = reliability(sd, "spencer", engine="taylor", search=False,
                                  debug_level=-1, num_slices=n)
        if not ok:
            return [f"Taylor reliability failed at num_slices={n}: {res}"]
        got[n] = res["F_MLV"]
    if got[None] != got[40]:
        fails.append(f"no slice count gives F_MLV {got[None]!r}, 40 slices "
                     f"{got[40]!r}: the default must stay 40")
    if got[60] == got[40]:
        fails.append(f"60 slices gives the 40-slice F_MLV {got[60]!r}: the slice "
                     f"count does not reach the run")
    sd = _load(RELIABILITY)
    with _quiet():
        ok, var = variance_contribution(sd, "spencer", search=False, num_slices=60)
    if not ok or var["F_MLV"] != got[60]:
        fails.append(f"variance contribution at 60 slices gives "
                     f"{var.get('F_MLV') if ok else var!r}, not the Taylor "
                     f"run's {got[60]!r}")
    return fails


CHECKS = [
    ("522 FEM pore pressure per material", check_fem_pore_pressure_per_material),
    ("523 suction only on mc/pow/hb", check_suction_effective_stress_only),
    ("524 seepage through a joint", check_seepage_through_joint),
    ("525 Size on SSR and joints polygons", check_size_on_overlays),
    ("526 joint t_cut follows Tension SRF", check_joint_tension_cutoff),
    ("FEM refuses a blank strength option", check_fem_blank_option),
    ("piezo line 2 only in rapid drawdown", check_piezo_line2_single_stage),
    ("report names each material's pore source", check_report_pore_basis),
    ("Taylor reliability takes a slice count", check_taylor_num_slices),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    failures = []
    for name, fn in CHECKS:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("audit engine checks:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("all passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
