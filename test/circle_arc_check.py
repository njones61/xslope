"""Slices on a circle: equal angles, and each base as long as its arc.

The slicer used to cut a circle into slices of equal width and give each base
its width over the cosine of the circle's slope at the slice's centerline.
Near a steep exit one such slice spans many degrees of arc, and both parts
fail there. The base length is short, because the base is curved, and the
centerline angle is too flat for the whole base. On Duncan, Wright and
Brandon's Figure 8.1 problem (a 20 ft, 2V:1H clay slope, phi = 0, c = 500 psf,
with a 10,000 lb/ft horizontal force at mid-height), the critical circle meets
the crest vertically. Its last slice was 0.58 ft wide and spanned 14 degrees
of arc, and its base read 3.4 ft for an arc of 4.8. The bases summed to
32.6 ft for an arc of 34.1, Bishop read 2.124 on that circle where the exact
answer is 2.2235, and the Bishop search came in 4% under the exact minimum.
Giving each base its arc length alone fixed the clay but made a sand slope
with pore pressure worse (OMS 1.4% and Janbu 1.1% low at ru = 0.3), because
the old short base had cancelled the too-flat angle in the water force.

A circle is now cut into equal angles, so slices narrow where it steepens,
and a base on the circle takes the length of its arc, 2 R asin(chord / 2R).
The slicer also keeps the surface's own end points when it merges slice
boundaries closer than 0.01; it used to drop the end for an interior point.

Each base on the circle takes its inclination from the chord between its
corners, which is the circle's tangent at the middle of the arc; its forces act
at the middle of that chord, on the slice's centerline and on the radius through
the arc's middle; and it reads its pore pressure at the middle of the arc, where
a smoothly varying u equals its average over the base to second order. The
point on the circle directly below the slice's center, where both were read
before, sits off the middle of the arc by up to a fifth of the base length on
the last slice of the vertical-exit circle, and read u there up to 5.6% high.

The checks use the same slope without the force, so the exact factor of safety
of the clay on a circle is c R L / (W x), with L the arc in the soil, W the
weight of the sliding mass and x its lever arm about the center (1.3453 on the
vertical-exit circle, where the old slicer read 1.2864):

1. The clay: OMS on the vertical-exit circle, and OMS and Bishop on a circle
   exiting the crest at 85 deg, match c R L / (W x) at 40 and 160 slices.
   (Bishop has no admissible solution at an exactly vertical exit: the last
   slice's base normal force goes negative.)
2. A sand slope with ru = 0.3 on the vertical-exit circle: OMS, Bishop, Janbu
   and Spencer at 40 slices match their own 400-slice values.
3. The bases sum to the length of the failure surface: on the vertical-exit
   circle; on it under a crest load starting 0.0066 ft inside the exit, where
   the old merge dropped the end and stopped 0.5 ft below the crest; and on a
   composite surface (a circle cut off by a firm layer 5 ft below the toe,
   exiting the crest at 88 deg), where the run along the layer keeps equal
   widths and straight bases.
4. A c-phi slope with ru = 0.3 and its mirror image give the same bases and
   factors of safety in OMS, Bishop, Janbu and Spencer, on each of those
   surfaces, and every one of those methods has a solution on both.
5. On the vertical-exit circle at 40 slices, every base's inclination is its
   chord's, its normal passes through the circle's center, and its u is within
   0.5% of the exact average along its arc, under a sloping water table and
   under ru = 0.3.

Run directly:  PYTHONPATH=. python3 test/circle_arc_check.py
"""

import contextlib
import copy
import io
import os
import tempfile

import numpy as np

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(_REPO, "docs", "lem", "files", "xslope_acads_simple.xlsx")
GAMMA = 100.0
CLAY = dict(c=500.0, phi=0.0, ru=0.0)
SAND = dict(c=0.0, phi=30.0, ru=0.3)
CPHI = dict(c=200.0, phi=30.0, ru=0.3)
_SLOPE = [(0.0, 0.0), (10.0, 20.0), (90.0, 20.0)]   # 20 ft at 2V:1H on a firm layer
# (ground, firm layer elevation, (circle center x, center y, radius))
STEEP = (_SLOPE, 0.0, (2.7317, 20.0049, 20.0049))   # tangent, meets the crest vertically
EXIT85 = (_SLOPE, 0.0, (3.0, 22.0, 22.0))           # tangent, meets the crest at 85 deg
COMPOSITE = ([(-30.0, 0.0)] + _SLOPE, -5.0, (3.0, 21.0, 28.0))   # dips 2 ft below the layer
LOAD_FROM = 22.73                      # 0.0066 ft inside the steep circle's exit
METHODS = ("oms", "bishop", "janbu", "spencer")
TOL = 2e-3                             # relative, on F: the old slicer missed by 1 to 4%
WATER = [(-40.0, 8.0), (-5.0, 8.0), (25.0, 19.5), (100.0, 19.5)]   # a sloping water table
GAMMA_W = 62.4
U_TOL = 5e-3          # relative, per base: u read below the slice's center missed by 5.6%


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _ground(model, mirror):
    pts = model[0]
    return [(-x, y) for x, y in reversed(pts)] if mirror else list(pts)


def _slope(model, soil=CLAY, mirror=False, load_from=None):
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
    sd = copy.deepcopy(_quiet(load_slope_data, BASE))
    sd["unit_system"] = "imperial"
    sd["gamma_water"] = 62.4
    sd["profile_lines"] = [{"coords": _ground(model, mirror), "mat_id": 0, "size": None}]
    sd["materials"] = [dict(sd["materials"][0], name="Soil", gamma=GAMMA, gamma_sat=None,
                            c=soil["c"], phi=soil["phi"], option="mc",
                            u="ru" if soil["ru"] else "none", ru=soil["ru"])]
    sd["max_depth"] = model[1]
    for key in ("dloads", "dloads2", "line_loads", "reinforcement_lines"):
        if key in sd:
            sd[key] = []
    if load_from is not None:
        sgn = -1.0 if mirror else 1.0
        load = [{"X": sgn * load_from, "Y": 20.0, "Normal": 100.0},
                {"X": sgn * 60.0, "Y": 20.0, "Normal": 100.0}]
        sd["dloads"] = [sorted(load, key=lambda p: p["X"])]
        sd["dload_dirs"] = ["normal"]
    sd["circles"] = [{"Xo": 0.0, "Yo": 30.0, "Depth": 0.0, "R": 30.0}]
    sd["circular"] = True
    path = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False).name
    try:
        _quiet(save_slope_data_to_xlsx, sd, path)
        return _quiet(load_slope_data, path)
    finally:
        os.unlink(path)


def _with_water(sd, line):
    sd = copy.deepcopy(sd)
    sd["piezo_line"] = line
    for m in sd["materials"]:
        m["u"] = "piezo"
    return sd


def _circle(model, mirror=False):
    xc, yc, r = model[2]
    return {"Xo": -xc if mirror else xc, "Yo": yc, "Depth": yc - r, "R": r}


def _slices(sd, circle, n, composite=False):
    from xslope.slice import generate_slices
    ok, out = _quiet(generate_slices, sd, circle=circle, num_slices=n,
                     composite=composite)
    if not ok:
        raise RuntimeError(f"generate_slices failed: {out}")
    return out[0]


def _fs(df, methods=METHODS):
    from xslope import solve
    out = {}
    for m in methods:
        ok, res = _quiet(getattr(solve, m), df.copy())
        out[m] = res["FS"] if ok else None
    return out


def _exact(model):
    """Surface length, and the clay's c R L / (W x), from shapely geometry with
    no slices: the sliding mass is the soil inside the circle, above the layer."""
    from shapely.geometry import LineString, Point, Polygon
    pts = _ground(model, False)
    floor = model[1]
    soil = Polygon(pts + [(pts[-1][0], floor), (pts[0][0], floor)])
    circle = _circle(model)
    disk = Point(circle["Xo"], circle["Yo"]).buffer(circle["R"], resolution=4096)
    mass = soil.intersection(disk)
    surface = mass.boundary.length - LineString(pts).intersection(disk).length
    arc = disk.exterior.intersection(soil).length
    lever = abs(mass.centroid.x - circle["Xo"])
    fs = CLAY["c"] * circle["R"] * arc / (GAMMA * mass.area * lever)
    return surface, fs


def check_clay_matches_the_exact_answer():
    fails = []
    for label, model, methods in (("the vertical-exit circle", STEEP, ("oms",)),
                                  ("the 85 deg exit", EXIT85, ("oms", "bishop"))):
        _surface, want = _exact(model)
        sd = _slope(model)
        for n in (40, 160):
            got = _fs(_slices(sd, _circle(model), n), methods)
            for m in methods:
                if got[m] is None:
                    fails.append(f"clay, {label}, {m}, {n} slices: no solution")
                elif abs(got[m] - want) > TOL * want:
                    fails.append(f"clay, {label}, {m}, {n} slices: F = {got[m]:.4f}; the "
                                 f"exact c R L / (W x) on this circle is {want:.4f}")
    return fails


def check_sand_with_water_converges():
    fails = []
    sd = _slope(STEEP, SAND)
    coarse = _fs(_slices(sd, _circle(STEEP), 40))
    fine = _fs(_slices(sd, _circle(STEEP), 400))
    for m in METHODS:
        a, b = coarse[m], fine[m]
        if a is None or b is None:
            fails.append(f"sand, ru = 0.3, {m}: no solution (40 slices {a}, 400 slices {b})")
        elif abs(a - b) > TOL * b:
            fails.append(f"sand, ru = 0.3, {m}: 40 slices give {a:.4f}, 400 give {b:.4f}")
    return fails


_SURFACES = (("the vertical-exit circle", STEEP, False, None),
             ("the vertical-exit circle under a crest load", STEEP, False, LOAD_FROM),
             ("the composite surface", COMPOSITE, True, None))


def check_bases_sum_to_the_surface():
    fails = []
    for label, model, composite, load in _SURFACES:
        want, _fs_exact = _exact(model)
        got = float(_slices(_slope(model, load_from=load), _circle(model), 40,
                            composite=composite)["dl"].sum())
        if abs(got - want) > 1e-4 * want:
            fails.append(f"{label}: the bases sum to {got:.4f} ft; the surface is "
                         f"{want:.4f} ft long")
    return fails


def check_mirror():
    fails = []
    for label, model, composite, load in _SURFACES:
        left = _slices(_slope(model, CPHI, load_from=load), _circle(model), 40,
                       composite=composite)
        right = _slices(_slope(model, CPHI, mirror=True, load_from=load),
                        _circle(model, mirror=True), 40, composite=composite)
        a, b = float(left["dl"].sum()), float(right["dl"].sum())
        if abs(a - b) > 1e-9 * a:
            fails.append(f"{label}: bases sum to {a:.6f} left-facing, {b:.6f} mirrored")
        fl, fr = _fs(left), _fs(right)
        for m in METHODS:
            if fl[m] is None or fr[m] is None:
                fails.append(f"{label}, c-phi with ru = 0.3, {m}: no solution "
                             f"(left-facing {fl[m]}, mirrored {fr[m]})")
            elif abs(fl[m] - fr[m]) > 1e-6 * fl[m]:
                fails.append(f"{label}, c-phi with ru = 0.3, {m}: left-facing "
                             f"{fl[m]:.6f}, mirrored {fr[m]:.6f}")
    return fails


def check_base_points():
    fails = []
    circle = _circle(STEEP)
    xo, yo, r = circle["Xo"], circle["Yo"], circle["R"]
    gx, gy = (np.array(v) for v in zip(*STEEP[0]))
    wx, wy = (np.array(v) for v in zip(*WATER))
    for label, sd, u_of in (
            ("a sloping water table", _with_water(_slope(STEEP, dict(CPHI, ru=0.0)), WATER),
             lambda x, y: GAMMA_W * (np.interp(x, wx, wy) - y)),
            ("ru = 0.3", _slope(STEEP, CPHI),
             lambda x, y: CPHI["ru"] * GAMMA * (np.interp(x, gx, gy) - y))):
        df = _slices(sd, circle, 40)
        scored = 0
        for s in df.itertuples():
            chord = np.degrees(np.arctan2(s.y_rb - s.y_lb, s.x_r - s.x_l))
            if abs(s.alpha - chord) > 1e-9:
                fails.append(f"{label}: the base at x = {s.x_c:.3f} is inclined "
                             f"{s.alpha:.4f} deg, its chord {chord:.4f}")
            a = np.radians(s.alpha)
            arm = (s.x_c - xo) * np.cos(a) + (s.y_cb - yo) * np.sin(a)
            if abs(arm) > 1e-9 * r:
                fails.append(f"{label}: the normal on the base at x = {s.x_c:.3f} "
                             f"passes {arm:.2e} from the circle's center")
            t = np.linspace(np.arcsin((s.x_l - xo) / r), np.arcsin((s.x_r - xo) / r), 2001)
            u = u_of(xo + r * np.sin(t), yo - r * np.cos(t))
            if u.min() <= 0.0:
                continue                    # dry, or wet over part of the base only
            avg = float(np.trapz(u, t) / (t[-1] - t[0]))
            scored += 1
            if abs(s.u - avg) > U_TOL * avg:
                fails.append(f"{label}: the base at x = {s.x_c:.3f} reads u = {s.u:.2f}; "
                             f"its average along the arc is {avg:.2f}")
        if scored < 30:
            fails.append(f"{label}: only {scored} bases are wet along their whole length")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_clay_matches_the_exact_answer, check_sand_with_water_converges,
                  check_bases_sum_to_the_surface, check_mirror, check_base_points):
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
    print("Circle slicing checks passed.")
