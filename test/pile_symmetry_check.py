"""Pile mirror-symmetry regression test.

A stabilizing pile must resist the sliding mass regardless of which way the
slope faces, so reflecting a model left-to-right (and reflecting the pile with
it) must leave every method's factor of safety unchanged. This guards the
right-facing pile-moment sign fix in oms()/bishop()/spencer().

Background: an earlier bug applied the pile's vertical-component moment with the
wrong sign on right-facing slopes (a real ~3-5% asymmetry for battered piles).
Two separate red herrings also bit during diagnosis and are pinned here too:
  - the convention is theta_p RELATIVE to the resisting direction (keep theta_p
    the same when mirroring, do NOT map theta_p -> 180 - theta_p);
  - the mirrored ground_surface must be sorted ascending in x (as load_slope_data
    produces) or the Ito & Matsui ground lookup (np.interp) returns garbage.

Run directly:  PYTHONPATH=. python3 test/pile_symmetry_check.py
"""

import copy

import numpy as np
from shapely.affinity import scale
from shapely.geometry import LineString

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath

PILE_MODEL = "docs/lem/files/xslope_piles.xlsx"
METHODS = [("oms", oms), ("bishop", bishop), ("spencer", spencer),
           ("janbu", janbu), ("corps", corps_engineers), ("lowe", lowe_karafiath)]
TOL_PCT = 0.05  # mirror asymmetry must be below this (well under arc-discretization noise)


def _mirror_x(geom):
    return scale(geom, xfact=-1.0, yfact=1.0, origin=(0, 0))


def mirror_slope_data(d):
    """Reflect a model about x=0, reproducing what load_slope_data would build
    for the mirror-image input (ground_surface re-sorted ascending in x)."""
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
    # Reflect the pile geometry; keep theta_p unchanged (relative-to-resisting convention).
    m["pile_lines"] = [dict(p, x1=-p["x1"], x2=-p["x2"]) for p in d["pile_lines"]]
    return m


def _solve_all(d, ns=40):
    ok, res = generate_slices(d, circle=d["circles"][0], num_slices=ns)
    assert ok, f"generate_slices failed: {res}"
    slice_df, _ = res
    out = {}
    for name, fn in METHODS:
        ok2, r = fn(slice_df.copy())
        out[name] = r["FS"] if ok2 else None
    return out, slice_df


def check(theta_p, expect_pile=True):
    base = load_slope_data(PILE_MODEL)
    for p in base["pile_lines"]:
        p["theta_p"] = theta_p
        if not expect_pile:
            p["H"] = 0.0  # control: no pile force
    left, sdf_l = _solve_all(base)
    right, sdf_r = _solve_all(mirror_slope_data(base))
    rf_l = bool(sdf_l["y_lt"].iat[0] > sdf_l["y_rt"].iat[-1])
    rf_r = bool(sdf_r["y_lt"].iat[0] > sdf_r["y_rt"].iat[-1])
    assert not rf_l and rf_r, f"expected left right_facing=False, mirror=True (got {rf_l}, {rf_r})"

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
            failures.append(f"{name}: asym {asym:.3f}% >= {TOL_PCT}% (left={fl:.5f}, mirror={fr:.5f})")
    return failures


def test_pile_mirror_symmetry():
    all_failures = []
    print("  control (no pile force), theta_p=0:")
    all_failures += check(0.0, expect_pile=False)
    for theta in (0.0, 30.0, -20.0):
        print(f"  pile, theta_p={theta:+.0f}:")
        all_failures += check(theta, expect_pile=True)
    assert not all_failures, "Pile mirror-symmetry violations:\n  " + "\n  ".join(all_failures)


if __name__ == "__main__":
    try:
        test_pile_mirror_symmetry()
    except AssertionError as e:
        print(f"\nFAILED:\n{e}")
        raise SystemExit(1)
    print("\nPASS: pile forces are mirror-symmetric across all methods.")
