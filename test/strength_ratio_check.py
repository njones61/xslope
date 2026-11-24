"""The strength ratio material (option = ratio): s_u = max(c, c/p * sigma'_v).

sigma'_v is the soil column above the slice base's pore-pressure point (gamma_sat
below the water table, as the 'ru' option weighs it), plus the vertical part of the
distributed load on the slice top, ponded water included, minus the pore pressure.
Checks against closed forms and identities, no tolerances to tune:

* **Planar wedge, exact.** On one plane, every method that satisfies force
  equilibrium along it gives F = sum(r (W_i - u_i b_i)) / (sin a cos a sum W_i),
  and with no water F = r / (sin a cos a), whatever the wedge's shape. Duncan,
  Wright & Brandon's Fig. 10.1 wedge, its material switched to a ratio.
* **The strength on every slice.** On the same wedge, each slice's c is
  max(c, r (W/b - u)), read back from the slice table.
* **Level ground, identical to c/p.** Under level ground with one soil and the water
  table at the surface, a ratio r and c/p = r gamma' with r-elev at the ground give
  the same strength at every slice base.
* **Ponded water.** A pond over that ground, entered as water standing above it,
  leaves every ratio strength unchanged: the pond's weight and the pore pressure it
  adds cancel in sigma'_v.
* **The minimum.** Where r sigma'_v < c, the strength is c.
* **Loader and stage 2.** The loader reads option = ratio; the two-stage seismic
  analysis's stage 2 strength factor reduces ratio slices (ratio_flag).
* **The finite element strengths.** Each element of a ratio material takes the same
  strength from its centroid: r gamma (ground - y) on the dry wedge; the c/p
  option's strength under level ground with the water at the surface; and the same
  again under a pond, whose load and pore pressure cancel.
* **The Gauss-point strengths the solver yields at.** The same three identities at
  every Gauss point of the prepared model, and both solution paths' working groups
  carry them reduced by the trial F.

Runs in a few seconds.
"""
import copy
import math
import os
import sys
import warnings

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from xslope import solve  # noqa: E402
from xslope.fileio import load_slope_data  # noqa: E402
from xslope.slice import generate_slices  # noqa: E402

WEDGE = os.path.join(ROOT, "docs", "verification", "files", "published", "dwb_fig10_1_wedge.xlsx")
FORCE_METHODS = ("spencer", "mprice", "corps", "lowe")
R = 0.25


def _ratio(sd, r=R, c=0.0, u="none"):
    sd = copy.deepcopy(sd)
    for m in sd["materials"]:
        m.update(option="ratio", c=c, cp=r, phi=0.0, d=0.0, psi=0.0, u=u)
    return sd


def _slices(sd, **kw):
    surf = {"non_circ": sd["non_circ"]} if sd.get("non_circ") else {"circle": sd["circles"][0]}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ok, out = generate_slices(sd, num_slices=kw.get("n", 30), **surf)
    if not ok:
        raise RuntimeError(out)
    return out[0]


def _fs(df, method):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ok, res = getattr(solve, method)(df.copy())
    return float(res["FS"]) if ok else None


def check_wedge_closed_form():
    fails = []
    sd = _ratio(load_slope_data(WEDGE))
    df = _slices(sd)
    a = math.radians(float(df["alpha"].iloc[0]))
    if not np.allclose(df["alpha"], df["alpha"].iloc[0]):
        return ["the wedge's slices are not on one plane"]
    want = R / (math.sin(a) * math.cos(a))
    for m in FORCE_METHODS:
        f = _fs(df, m)
        if f is None or abs(f - want) > 1e-6 * want:
            fails.append(f"wedge, {m}: F = {f}, closed form r/(sin a cos a) = {want:.6f}")
    # each slice's strength: c/p times W/b (dry, no load)
    got = df["c"].to_numpy()
    exp = R * df["w"].to_numpy() / df["dx"].to_numpy()
    if not np.allclose(got, exp, rtol=1e-9):
        fails.append(f"wedge slice strengths differ from r W/b by up to "
                     f"{np.max(np.abs(got - exp)):.3g}")
    if not (df["ratio_flag"].all() and (df["phi"] == 0).all()):
        fails.append("wedge slices: ratio_flag not set, or phi not 0")
    return fails


def check_wedge_minimum():
    sd = load_slope_data(WEDGE)
    df = _slices(_ratio(sd, c=0.0))
    cmin = float(np.median(df["c"]))
    df2 = _slices(_ratio(sd, c=cmin))
    want = np.maximum(cmin, df["c"].to_numpy())
    if not np.allclose(df2["c"].to_numpy(), want):
        return ["the minimum strength c is not applied as max(c, r sigma'_v)"]
    if not ((df2["c"] == cmin).any() and (df2["c"] > cmin).any()):
        return ["the minimum test did not exercise both sides of max()"]
    return []


def _level(sd_wedge, water=False, pond=0.0):
    """Level ground at y = 10 over the whole model, from the wedge model's own
    geometry and sheets: one soil, gamma 18, gamma_sat 20, a circle through it."""
    from shapely.geometry import LineString, Polygon
    sd = copy.deepcopy(sd_wedge)
    x0, x1, ytop, ybot = -20.0, 40.0, 10.0, -20.0
    sd["ground_surface"] = LineString([(x0, ytop), (x1, ytop)])
    sd["profile_lines"] = [{"coords": [(x0, ytop), (x1, ytop)], "mat_id": 0}]
    sd["max_depth"] = ybot
    sd["polygons"] = [{"polygon": Polygon([(x0, ytop), (x1, ytop), (x1, ybot), (x0, ybot)]), "mat_id": 0}]
    sd["materials"] = sd["materials"][:1]
    sd["materials"][0].update(gamma=18.0, gamma_sat=20.0)
    sd["non_circ"] = []
    sd["circles"] = [{"Xo": 10.0, "Yo": 20.0, "Depth": -2.0, "R": 22.0}]
    if water or pond:
        level = ytop + pond
        sd["piezo_line"] = [(x0, level), (x1, level)]
        for m in sd["materials"]:
            m["u"] = "piezo"
    return sd


def check_level_ground_matches_cp():
    fails = []
    base = load_slope_data(WEDGE)
    try:
        sd = _level(base, water=True)
        gw = float(sd["gamma_water"])
        ratio = _ratio(sd, u="piezo")
        cp = copy.deepcopy(ratio)
        for m in cp["materials"]:
            m.update(option="cp", c=0.0, cp=R * (20.0 - gw), r_elev=10.0)
        dr, dc = _slices(ratio), _slices(cp)
    except Exception as exc:              # the level model is built here, by hand
        return [f"level-ground model failed to slice: {exc!r}"]
    if len(dr) != len(dc) or not np.allclose(dr["c"], dc["c"], rtol=1e-9, atol=1e-9):
        fails.append(f"level ground: ratio strengths differ from c/p = r gamma' by up to "
                     f"{np.max(np.abs(dr['c'].to_numpy() - dc['c'].to_numpy())):.3g}")
    return fails


def check_pond_cancels():
    base = load_slope_data(WEDGE)
    try:
        dry_top = _slices(_ratio(_level(base, water=True), u="piezo"))
        ponded = _ratio(_level(base, pond=3.0), u="piezo")
        ponded["water_loads"] = "auto"
        dp = _slices(ponded)
    except Exception as exc:
        return [f"ponded level-ground model failed to slice: {exc!r}"]
    if dp["dload"].abs().sum() <= 0:
        return ["the pond put no load on the slices (automatic water loads off?)"]
    if not np.allclose(dp["c"], dry_top["c"], rtol=1e-9, atol=1e-6):
        return [f"a 3 m pond changed the ratio strength by up to "
                f"{np.max(np.abs(dp['c'].to_numpy() - dry_top['c'].to_numpy())):.3g}; "
                f"its weight and its pore pressure should cancel in sigma'_v"]
    return []


def check_two_stage_factor_reduces_ratio():
    from xslope.advanced import seismic_two_stage
    sd = _ratio(load_slope_data(WEDGE), r=0.6)   # F = 1.39 before the earthquake
    sd["k_seismic"] = 0.05
    df = _slices(sd)
    one, red = df.copy(), df.copy()
    ok1, r1 = seismic_two_stage(one, "spencer", debug_level=0, strength="kc1")
    ok2, r2 = seismic_two_stage(red, "spencer", debug_level=0,
                                strength={"strength": "kc1", "factor": 0.8})
    if not (ok1 and ok2):
        # a wedge with no d/psi material is refused by preflight but solves here
        return [f"two-stage on the ratio wedge did not solve: {r1 if not ok1 else r2}"]
    if not np.allclose(red["c"], 0.8 * one["c"]):
        return ["the stage 2 strength factor did not reduce the ratio slices"]
    if abs(r2["FS"] - 0.8 * r1["FS"]) > 1e-6 * r1["FS"]:
        return [f"stage 2 F at factor 0.8 is {r2['FS']:.6f}, not 0.8 x {r1['FS']:.6f}"]
    return []


def _fem(sd):
    from xslope.fem import build_fem_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    from xslope.water import with_water_loads
    sd = copy.deepcopy(sd)
    for m in sd["materials"]:
        m.update(E=1.0e5, nu=0.3)
    sd = with_water_loads(sd)
    mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size=2.5,
                                    element_type="tri6", debug=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_fem_data(sd, mesh)


def _centroids(fd):
    nodes, el, et = fd["nodes"], fd["elements"], fd["element_types"]
    return np.array([nodes[el[e][:et[e]]].mean(axis=0) for e in range(len(el))])


def check_fem_wedge_dry():
    sd = _ratio(load_slope_data(WEDGE))
    fd = _fem(sd)
    xy = _centroids(fd)
    g = float(sd["materials"][0]["gamma"])
    gx, gy = np.array(sd["ground_surface"].coords).T
    want = R * g * np.maximum(0.0, np.interp(xy[:, 0], gx, gy) - xy[:, 1])
    if not np.allclose(fd["c_by_elem"], want, rtol=1e-9, atol=1e-9):
        return [f"FEM, dry wedge: element strengths differ from r gamma (ground - y) by up to "
                f"{np.max(np.abs(fd['c_by_elem'] - want)):.3g}"]
    if np.any(fd["phi_by_elem"] != 0.0):
        return ["FEM, dry wedge: a ratio element has phi != 0"]
    return []


def _gauss_points(fd):
    """The prepared model, and each static group's Gauss-point x, y and
    strength-ratio s_u (NaN off a ratio material)."""
    from xslope.fem import _prepare_fem_model
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prep = _prepare_fem_model(fd)
    nodes, el, et = fd["nodes"], fd["elements"], fd["element_types"]
    out = []
    for grp in prep["gp_groups_static"]:
        xy = np.array([prep["elem_gp_data"][e][g]["N"] @ nodes[el[e][:et[e]]]
                       for e, g in grp["pairs"]])
        out.append((xy[:, 0], xy[:, 1], grp.get("c_ratio")))
    return prep, out


def check_fem_gauss_points():
    fails = []
    sd = _ratio(load_slope_data(WEDGE))
    g = float(sd["materials"][0]["gamma"])
    gx, gy = np.array(sd["ground_surface"].coords).T
    fd = _fem(sd)
    prep, groups = _gauss_points(fd)
    for x, y, cr in groups:
        want = R * g * np.maximum(0.0, np.interp(x, gx, gy) - y)
        if cr is None or not np.allclose(cr, want, rtol=1e-9, atol=1e-9):
            fails.append("FEM, dry wedge: Gauss-point strengths differ from r gamma "
                         "(ground - y) at the point")
    # Both paths' working groups reduce them by the trial F.
    from xslope.fem import _nr_build_groups
    F = 1.7
    n = len(fd["elements"])
    c_red = fd["c_by_elem"] / F
    groups_nr = _nr_build_groups(prep, c_red, np.zeros(n), None,
                                 finv_by_elem=np.full(n, 1.0 / F))
    for grp, (_x, _y, cr) in zip(groups_nr, groups):
        if not np.allclose(grp["c_r"], cr / F, rtol=1e-12, atol=0.0):
            fails.append("FEM, Newton groups: Gauss-point strengths are not s_u / F")
    base = load_slope_data(WEDGE)
    wet = _level(base, water=True)
    gw = float(wet["gamma_water"])
    pond = _ratio(_level(base, pond=3.0), u="piezo")
    pond["water_loads"] = "auto"
    for label, model in (("level ground", _ratio(wet, u="piezo")), ("a 3 m pond", pond)):
        _p, gps = _gauss_points(_fem(model))
        for x, y, cr in gps:
            want = R * (20.0 - gw) * np.maximum(0.0, 10.0 - y)
            if cr is None or not np.allclose(cr, want, rtol=1e-9, atol=1e-6):
                fails.append(f"FEM, {label}: Gauss-point strengths differ from "
                             "r gamma' (ground - y)")
                break
    return fails


def check_fem_level_ground_and_pond():
    fails = []
    base = load_slope_data(WEDGE)
    wet = _level(base, water=True)
    gw = float(wet["gamma_water"])
    ratio = _ratio(wet, u="piezo")
    cp = copy.deepcopy(ratio)
    for m in cp["materials"]:
        m.update(option="cp", c=0.0, cp=R * (20.0 - gw), r_elev=10.0)
    pond = _ratio(_level(base, pond=3.0), u="piezo")
    pond["water_loads"] = "auto"
    fr, fc, fp = _fem(ratio), _fem(cp), _fem(pond)
    if not np.allclose(fr["c_by_elem"], fc["c_by_elem"], rtol=1e-9, atol=1e-9):
        fails.append(f"FEM, level ground: ratio strengths differ from c/p = r gamma' by up to "
                     f"{np.max(np.abs(fr['c_by_elem'] - fc['c_by_elem'])):.3g}")
    if len(fp["c_by_elem"]) != len(fr["c_by_elem"]) or \
            not np.allclose(fp["c_by_elem"], fr["c_by_elem"], rtol=1e-9, atol=1e-6):
        fails.append("FEM, a 3 m pond changed the ratio strengths; its load and pore "
                     "pressure should cancel")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_wedge_closed_form, check_wedge_minimum, check_level_ground_matches_cp,
                  check_pond_cancels, check_two_stage_factor_reduces_ratio,
                  check_fem_wedge_dry, check_fem_level_ground_and_pond,
                  check_fem_gauss_points):
        try:
            fails += check()
        except Exception as exc:          # a crash is a failure, not a pass
            fails.append(f"{check.__name__} raised {exc!r}")
    return fails


if __name__ == "__main__":
    failures = run()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("Strength ratio checks passed.")
