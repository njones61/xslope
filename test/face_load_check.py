"""A distributed load that starts at the crest of a wall stays on the crest.

A wall face is often drawn leaning by a hair, (-0.001, 0) to (0, 25), so the
ground stays a function of x. The slicer merges slice boundaries closer than
0.01, which drops both of the face's vertices and the load's first point, so
one slice runs from just before the face to past the crest. A slice top is the
chord between its two top corners, and the load was laid along that chord with
its intensity read by x: the 150 psf surcharge that starts at the crest ramped
up along the 25 ft face and was applied as a horizontal push into the wall,
1,876 lb/ft on that slice, where 92 lb/ft acts downward on the crest. It raised
the factor of safety of a nailed wall from 1.310 to 1.400 (found reconciling
the FHWA GEC 7 Appendix C wall, 2026-11-05).

Such a slice now takes its load from the ground under it
(:func:`xslope.slice.load_on_ground_pieces`), with the load placed by distance
along the ground, so the face, which lies before the load's first point,
carries none. Every slice whose top is one straight piece of ground keeps the
chord calculation.

Checks, on a 25 ft wall in a c-phi soil with the surcharge from the crest to
25 ft back and a circle that passes under the toe:

1. No slice carries more load than its plan width times the surcharge, and
   the slices carry the load on the crest the circle encloses, no more.
2. The leaning face and an exactly vertical face give the same factor of
   safety, as do the leaning face with the load starting 0.01 ft behind the
   crest. Before the fix the leaning face read 2 to 4% high (OMS 1.567 against
   1.532, Spencer 1.757 against 1.717).
3. The leaning model and its mirror image give the same factor of safety, so
   the inclination of the new resultant follows the right-facing convention.

Checks 2 and 3 cover OMS, Bishop, Janbu, Spencer and Morgenstern-Price. The
Corps of Engineers and Lowe-Karafiath methods are left out: they incline the
side forces by the slope of each slice top, and at a wall face one slice's top
runs up the face whichever way the face is drawn, so those two read the face's
geometry as well as its load (Corps 1.409 leaning, 1.439 vertical; Lowe-
Karafiath finds no solution on the leaning face, before the fix as after).

Run directly:  PYTHONPATH=. python3 test/face_load_check.py
"""

import contextlib
import copy
import io
import math
import os
import tempfile

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(_REPO, "docs", "lem", "files", "xslope_acads_simple.xlsx")
H, Q, LOAD_END = 25.0, 150.0, 25.0
CIRCLE = (-5.0, 35.0, 40.0)            # center x, center y, radius: under the toe
METHODS = ("oms", "bishop", "janbu", "spencer", "mprice")
TOL = 1e-3                             # relative: a 0.001 ft lean moves nothing


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _wall(lean=0.001, load_x0=0.0, mirror=False):
    """The wall as a loaded model: written to a workbook and read back, so the
    ground surface and polygons are the loader's own."""
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
    sd = copy.deepcopy(_quiet(load_slope_data, BASE))
    sgn = -1.0 if mirror else 1.0
    pts = [(-60.0, 0.0), (-lean, 0.0), (0.0, H), (60.0, H)]
    if mirror:
        pts = [(-x, y) for x, y in reversed(pts)]
    sd["unit_system"] = "imperial"
    sd["gamma_water"] = 62.4
    sd["profile_lines"] = [{"coords": pts, "mat_id": 0, "size": None}]
    sd["materials"] = [dict(sd["materials"][0], name="Silty sand", gamma=125.0,
                            gamma_sat=None, c=100.0, phi=32.0, option="mc", u="none")]
    sd["max_depth"] = -20.0
    load = [{"X": sgn * load_x0, "Y": H, "Normal": Q},
            {"X": sgn * LOAD_END, "Y": H, "Normal": Q}]
    sd["dloads"] = [sorted(load, key=lambda p: p["X"])]
    sd["dload_dirs"] = ["normal"]
    xc, yc, r = CIRCLE
    sd["circles"] = [{"Xo": sgn * xc, "Yo": yc, "Depth": yc - r, "R": r}]
    sd["circular"] = True
    path = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False).name
    try:
        _quiet(save_slope_data_to_xlsx, sd, path)
        return _quiet(load_slope_data, path)
    finally:
        os.unlink(path)


def _solve(sd):
    from xslope import solve
    from xslope.slice import generate_slices
    ok, out = _quiet(generate_slices, sd, circle=sd["circles"][0], num_slices=40)
    if not ok:
        raise RuntimeError(f"generate_slices failed: {out}")
    df = out[0]
    fs = {}
    for m in METHODS:
        ok, res = _quiet(getattr(solve, m), df.copy())
        fs[m] = res["FS"] if ok else None
    return df, fs


def check_load_stays_on_the_crest():
    fails = []
    df, _fs = _solve(_wall())
    # The crest the circle encloses: from the face to where it leaves the top.
    xc, yc, r = CIRCLE
    x_exit = xc + math.sqrt(r * r - (H - yc) ** 2)
    want = Q * (min(x_exit, LOAD_END) - 0.0)
    total = float(df["dload"].sum())
    if abs(total - want) > 1e-6 * want:
        fails.append(f"the slices carry {total:,.1f} lb/ft of surcharge; the crest "
                     f"inside the circle carries {want:,.1f}")
    over = df[df["dload"] > Q * df["dx"] * (1 + 1e-9) + 1e-9]
    for _i, row in over.iterrows():
        fails.append(f"slice {int(row['slice #'])} ({row['x_l']:.4f} to {row['x_r']:.4f}) "
                     f"carries {row['dload']:,.1f} lb/ft, more than its "
                     f"{Q * row['dx']:,.1f} lb/ft of crest")
    return fails


def check_lean_changes_nothing():
    fails = []
    _df, lean = _solve(_wall())
    for label, sd in (("an exactly vertical face", _wall(lean=0.0)),
                      ("the load starting 0.01 behind the crest", _wall(load_x0=0.01))):
        _df2, other = _solve(sd)
        for m in METHODS:
            a, b = lean[m], other[m]
            if a is None or b is None:
                fails.append(f"{m}: no solution (leaning {a}, {label} {b})")
            elif abs(a - b) > TOL * b:
                fails.append(f"{m}: the leaning face gives {a:.4f}, {label} {b:.4f}")
    return fails


def check_mirror():
    fails = []
    _df, left = _solve(_wall())
    _df2, right = _solve(_wall(mirror=True))
    for m in METHODS:
        a, b = left[m], right[m]
        if a is None or b is None or abs(a - b) > 1e-4 * a:
            fails.append(f"{m}: left-facing {a}, mirror {b}")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_load_stays_on_the_crest, check_lean_changes_nothing,
                  check_mirror):
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
    print("Face-load checks passed.")
