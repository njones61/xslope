"""Tension-crack mirror-symmetry regression test.

A tension crack — and the water thrust inside it when the crack is full — must
act the same way regardless of which way the slope faces, so reflecting a model
left-to-right must leave every method's factor of safety unchanged. This guards
the right-facing sign of the tension-crack water force across all six methods.

Background (audit finding, June 2026): the dry crack (geometry truncation only)
was already symmetric, but the water-filled crack was asymmetric by 5-8% for
oms/bishop/janbu/corps/lowe — only Spencer was correct. The cause was a sign
mismatch: slice.py pre-negated the crack water force on right-facing slopes for
Spencer's real-coordinate force balance, but the orientation-normalized methods
consume it as a driving-positive magnitude and did not want that negation. The
fix stores the force as a positive magnitude (like the seismic kw) and lets
Spencer flip it in its own right-facing swap block. Two checks are pinned here:

  1. mirror symmetry — a model and its left-right reflection give the same FS for
     every method, both dry and water-filled; and
  2. physical sign — a water-filled crack must LOWER the FS relative to a dry
     crack of the same depth (the water thrust is destabilizing).

Run directly:  PYTHONPATH=. python3 test/tension_crack_symmetry_check.py
"""

import contextlib
import copy
import os
import io

from shapely.affinity import scale
from shapely.geometry import LineString

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import oms, bishop, spencer, janbu, corps, lowe

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(_REPO, "docs/lem/files/xslope_acads_simple.xlsx")
TCRACK_DEPTH = 4.0
TCRACK_WATER = 3.0
METHODS = [("oms", oms), ("bishop", bishop), ("janbu", janbu),
           ("corps", corps), ("lowe", lowe), ("spencer", spencer)]
TOL_PCT = 0.05  # mirror asymmetry must be below this (well under arc-discretization noise)


def _mirror_x(geom):
    return scale(geom, xfact=-1.0, yfact=1.0, origin=(0, 0))


def mirror_slope_data(d):
    """Reflect a model about x=0, reproducing what load_slope_data would build for
    the mirror-image input (ground_surface re-sorted ascending in x)."""
    m = copy.deepcopy(d)
    gs = _mirror_x(d["ground_surface"])
    m["ground_surface"] = LineString(sorted(list(gs.coords), key=lambda c: c[0]))
    if d.get("domain_polygon") is not None:
        m["domain_polygon"] = _mirror_x(d["domain_polygon"])
    m["polygons"] = [{"polygon": _mirror_x(p["polygon"]), "mat_id": p["mat_id"]}
                     for p in d["polygons"]]
    if d.get("profile_lines"):
        m["profile_lines"] = [{"coords": [(-x, y) for (x, y) in p["coords"]],
                               "mat_id": p["mat_id"]} for p in d["profile_lines"]]
    m["circles"] = [dict(c, Xo=-c["Xo"]) for c in d["circles"]]
    return m


def _solve_all(d, ns=40):
    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(d, circle=d["circles"][0], num_slices=ns)
    assert ok, f"generate_slices failed: {res}"
    slice_df, _ = res
    out = {}
    for name, fn in METHODS:
        with contextlib.redirect_stdout(io.StringIO()):
            ok2, r = fn(slice_df.copy())
        out[name] = r["FS"] if ok2 else None
    return out, slice_df


def _load(water):
    d = load_slope_data(MODEL)
    d["tcrack_depth"] = TCRACK_DEPTH
    d["tcrack_water"] = water
    return d


def check_symmetry(water, label):
    base = _load(water)
    left, sdf_l = _solve_all(base)
    right, sdf_r = _solve_all(mirror_slope_data(base))
    rf_l = bool(sdf_l["y_lt"].iat[0] > sdf_l["y_rt"].iat[-1])
    rf_r = bool(sdf_r["y_lt"].iat[0] > sdf_r["y_rt"].iat[-1])
    assert not rf_l and rf_r, f"expected left right_facing=False, mirror=True (got {rf_l}, {rf_r})"

    print(f"  {label} crack (depth={TCRACK_DEPTH}, water={water}):")
    failures = []
    for name, _ in METHODS:
        fl, fr = left[name], right[name]
        if fl is None or fr is None:
            failures.append(f"{name}: solve failed (left={fl}, right={fr})")
            continue
        asym = abs(fl - fr) / fl * 100
        flag = "ok" if asym < TOL_PCT else "FAIL"
        print(f"    {name:8s} left={fl:.6f} mirror={fr:.6f}  asym={asym:.4f}%  {flag}")
        if asym >= TOL_PCT:
            failures.append(f"{name}: {asym:.3f}% mirror asymmetry (left={fl:.5f}, mirror={fr:.5f})")
    return failures


def check_water_lowers_fs():
    """A water-filled crack must lower the FS relative to a dry crack of the same depth."""
    dry, _ = _solve_all(_load(0.0))
    wet, _ = _solve_all(_load(TCRACK_WATER))
    print("  water must LOWER FS vs dry crack:")
    failures = []
    for name, _ in METHODS:
        ok = wet[name] < dry[name]
        print(f"    {name:8s} dry={dry[name]:.5f} -> wet={wet[name]:.5f}  {'ok' if ok else 'WRONG'}")
        if not ok:
            failures.append(f"{name}: water did not lower FS (dry={dry[name]:.5f}, wet={wet[name]:.5f})")
    return failures


def main():
    print("Tension-crack mirror-symmetry + physical-sign check:")
    failures = []
    failures += check_symmetry(0.0, "DRY")
    failures += check_symmetry(TCRACK_WATER, "WATER")
    failures += check_water_lowers_fs()
    return failures


def run():
    """Failures as a list, for run_tests.py."""
    return main()


def _cli():
    failures = main()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll tension-crack symmetry and sign checks passed.")


if __name__ == "__main__":
    _cli()
