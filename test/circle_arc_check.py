"""The base of a slice on a circle is as long as the arc it follows.

The slicer cuts slices of equal width. Each base's length used to be its width
over the cosine of the circle's slope at the slice's centerline, which is exact
for a straight base and short for a curved one. Near a steep exit the circle
turns through many degrees within one slice, so the error grows: on a 20 ft
clay slope (phi = 0, c = 500 psf) with a circle meeting the crest vertically,
the last slice was 0.57 ft wide and its base read 3.4 ft for an arc of 4.8.
The bases summed to 32.6 ft for an arc of 34.1, and OMS and Bishop read 2.124
where the exact answer is 2.2235. A search seeks those circles out, because
the error lowers F there: on Duncan, Wright and Brandon's Figure 8.1 problem
the Bishop search came in 4% under the exact minimum (found 2026-11-05).

A base on the circle now takes the length of its arc, 2 R asin(chord / 2R).

Checks, on that slope with no reinforcement, so the exact factor of safety on
a circle is c R L / (W x), with L the arc in the soil, W the weight of the
sliding mass and x its lever arm about the center:

1. OMS and Bishop on the vertical-exit circle match c R L / (W x) at 40 and
   160 slices.
2. The bases sum to the length of the failure surface: on that circle, and on
   a composite surface (a circle cut off by a firm layer 5 ft below the toe,
   exiting the crest at 88 deg), where the run along the layer keeps its
   straight-base length.
3. The slope and its mirror image give the same bases and factors of safety.

Run directly:  PYTHONPATH=. python3 test/circle_arc_check.py
"""

import contextlib
import copy
import io
import os
import tempfile

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(_REPO, "docs", "lem", "files", "xslope_acads_simple.xlsx")
C, GAMMA = 500.0, 100.0
# (ground, firm layer elevation, circle center x, center y, radius)
STEEP = ([(0.0, 0.0), (10.0, 20.0), (90.0, 20.0)], 0.0,    # 20 ft at 2V:1H on the layer;
         (2.7317, 20.0049, 20.0049))                       # tangent, meets the crest vertically
COMPOSITE = ([(-30.0, 0.0), (0.0, 0.0), (10.0, 20.0), (90.0, 20.0)], -5.0,
             (3.0, 21.0, 28.0))                            # dips 2 ft below the layer
TOL = 3e-3                             # relative, on F; the old bases missed by 4.5%


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _ground(model, mirror):
    pts = model[0]
    return [(-x, y) for x, y in reversed(pts)] if mirror else list(pts)


def _slope(model, mirror=False):
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
    sd = copy.deepcopy(_quiet(load_slope_data, BASE))
    pts = _ground(model, mirror)
    sd["unit_system"] = "imperial"
    sd["gamma_water"] = 62.4
    sd["profile_lines"] = [{"coords": pts, "mat_id": 0, "size": None}]
    sd["materials"] = [dict(sd["materials"][0], name="Clay", gamma=GAMMA,
                            gamma_sat=None, c=C, phi=0.0, option="mc", u="none")]
    sd["max_depth"] = model[1]
    for key in ("dloads", "dloads2", "line_loads", "reinforcement_lines"):
        if key in sd:
            sd[key] = []
    sd["circles"] = [{"Xo": 0.0, "Yo": 30.0, "Depth": 0.0, "R": 30.0}]
    sd["circular"] = True
    path = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False).name
    try:
        _quiet(save_slope_data_to_xlsx, sd, path)
        return _quiet(load_slope_data, path)
    finally:
        os.unlink(path)


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


def _exact(model, mirror=False):
    """Surface length, and c R L / (W x) for a circle in the soil, from shapely
    geometry with no slices: the sliding mass is the soil inside the circle,
    above the firm layer."""
    from shapely.geometry import LineString, Point, Polygon
    pts = _ground(model, mirror)
    floor = model[1]
    soil = Polygon(pts + [(pts[-1][0], floor), (pts[0][0], floor)])
    circle = _circle(model, mirror)
    disk = Point(circle["Xo"], circle["Yo"]).buffer(circle["R"], resolution=4096)
    mass = soil.intersection(disk)
    ground_in = LineString(pts).intersection(disk).length
    surface = mass.boundary.length - ground_in
    arc = disk.exterior.intersection(soil).length
    lever = abs(mass.centroid.x - circle["Xo"])
    fs = C * circle["R"] * arc / (GAMMA * mass.area * lever)
    return surface, fs


def check_steep_exit_matches_the_exact_answer():
    from xslope import solve
    fails = []
    sd = _slope(STEEP)
    circle = _circle(STEEP)
    _surface, want = _exact(STEEP)
    for n in (40, 160):
        df = _slices(sd, circle, n)
        for m in ("oms", "bishop"):
            ok, res = _quiet(getattr(solve, m), df.copy())
            if not ok:
                fails.append(f"{m}, {n} slices: no solution ({res})")
            elif abs(res["FS"] - want) > TOL * want:
                fails.append(f"{m}, {n} slices: F = {res['FS']:.4f}; the exact answer "
                             f"c R L / (W x) on this circle is {want:.4f}")
    return fails


def check_bases_sum_to_the_surface():
    fails = []
    for label, model, composite in (("the vertical-exit circle", STEEP, False),
                                    ("the composite surface", COMPOSITE, True)):
        want, _fs = _exact(model)
        got = float(_slices(_slope(model), _circle(model), 40,
                            composite=composite)["dl"].sum())
        if abs(got - want) > 1e-4 * want:
            fails.append(f"{label}: the bases sum to {got:.4f} ft; the surface is "
                         f"{want:.4f} ft long")
    return fails


def check_mirror():
    from xslope import solve
    fails = []
    for label, model, composite in (("the vertical-exit circle", STEEP, False),
                                    ("the composite surface", COMPOSITE, True)):
        left = _slices(_slope(model), _circle(model), 40, composite=composite)
        right = _slices(_slope(model, mirror=True), _circle(model, mirror=True), 40,
                        composite=composite)
        a, b = float(left["dl"].sum()), float(right["dl"].sum())
        if abs(a - b) > 1e-9 * a:
            fails.append(f"{label}: bases sum to {a:.6f} left-facing, {b:.6f} mirrored")
        for m in ("oms", "bishop", "spencer"):
            ok1, r1 = _quiet(getattr(solve, m), left.copy())
            ok2, r2 = _quiet(getattr(solve, m), right.copy())
            f1 = r1["FS"] if ok1 else None
            f2 = r2["FS"] if ok2 else None
            if (f1 is None) != (f2 is None) or (f1 is not None and abs(f1 - f2) > 1e-6 * f1):
                fails.append(f"{label}, {m}: left-facing {f1}, mirrored {f2}")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_steep_exit_matches_the_exact_answer,
                  check_bases_sum_to_the_surface, check_mirror):
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
    print("Circle arc-length checks passed.")
