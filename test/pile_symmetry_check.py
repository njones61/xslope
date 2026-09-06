"""A slope and its mirror image must give the same factor of safety.


A stabilizing pile must resist the sliding mass regardless of which way the
slope faces, so reflecting a model left-to-right (and reflecting the pile with
it) must leave every method's factor of safety unchanged. The pile model is
where a reflection touches the most machinery — the pile geometry, its
horizontal force, the arms every method takes it on, and the ground lookup the
Ito & Matsui capacity reads — so it is the model this is checked on.

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
import os

import numpy as np
from shapely.affinity import scale
from shapely.geometry import LineString

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import oms, bishop, spencer, janbu, corps, lowe

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PILE_MODEL = os.path.join(_REPO, "docs/lem/files/xslope_piles.xlsx")
METHODS = [("oms", oms), ("bishop", bishop), ("spencer", spencer),
           ("janbu", janbu), ("corps", corps), ("lowe", lowe)]
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
    ok, res = generate_slices(d, circle=d["circles"][0], num_slices=ns, debug=False)
    assert ok, f"generate_slices failed: {res}"
    slice_df, _ = res
    out = {}
    for name, fn in METHODS:
        ok2, r = fn(slice_df.copy())
        out[name] = r["FS"] if ok2 else None
    return out, slice_df


def check(theta_p, expect_pile=True, ns=40, quiet=False):
    base = load_slope_data(PILE_MODEL)
    if expect_pile:
        for p in base["pile_lines"]:
            p["theta_p"] = theta_p
    else:
        # The control is the model with no pile at all. It used to be the model
        # with H = 0, which preflight now rejects — H is a capacity, and a
        # capacity of zero is an input error rather than a pile that does
        # nothing — so the pile lines are removed instead. Same control, stated
        # honestly.
        base["pile_lines"] = []
    left, sdf_l = _solve_all(base, ns)
    right, sdf_r = _solve_all(mirror_slope_data(base), ns)
    rf_l = bool(sdf_l["y_lt"].iat[0] > sdf_l["y_rt"].iat[-1])
    rf_r = bool(sdf_r["y_lt"].iat[0] > sdf_r["y_rt"].iat[-1])
    assert not rf_l and rf_r, f"expected left right_facing=False, mirror=True (got {rf_l}, {rf_r})"

    failures = []
    worst = (0.0, None)
    for name, _ in METHODS:
        fl, fr = left[name], right[name]
        if fl is None or fr is None:
            failures.append(f"{name}: solve failed (left={fl}, right={fr})")
            continue
        asym = abs(fl - fr) / fl * 100
        if asym > worst[0]:
            worst = (asym, name)
        if asym >= TOL_PCT:
            failures.append(f"{name} at {ns} slices: asym {asym:.3f}% >= {TOL_PCT}% "
                            f"(left={fl:.5f}, mirror={fr:.5f})")
    if not quiet and worst[1] is not None:
        print(f"    {ns:3d} slices: worst {worst[1]} {worst[0]:.4f}%")
    return failures


#: Slice counts the pair is checked at. Which slice a pile is credited to used to
#: depend on the order the slices were built in, so the asymmetry appeared only
#: where a pile landed exactly on a slice boundary — at 30 and 60 slices on this
#: model, where Corps of Engineers moved 1.44% and 0.72% between the model and
#: its mirror, and nowhere else. A single slice count would have missed it.
SLICE_COUNTS = (30, 40, 50, 60, 80)


def leg_the_mirror_pair_agrees():
    """Every method, both facings, at every slice count."""
    failures = []
    print("  control (no pile):")
    for ns in SLICE_COUNTS:
        failures += check(0.0, expect_pile=False, ns=ns)
    for theta in (0.0, 30.0, -20.0):
        print(f"  pile, theta_p = {theta:+.0f}:")
        for ns in SLICE_COUNTS:
            failures += check(theta, expect_pile=True, ns=ns)
    return failures


def leg_a_pile_lands_on_a_boundary():
    """The case only bites where a pile sits exactly on a slice corner.

    If no slice count in SLICE_COUNTS puts a pile on a boundary any more, the
    leg above is no longer testing what it was written for and says so.
    """
    base = load_slope_data(PILE_MODEL)
    on_boundary = []
    for ns in SLICE_COUNTS:
        ok, res = generate_slices(base, circle=base["circles"][0], num_slices=ns,
                                  debug=False)
        if not ok:
            continue
        df = res[0]
        rows = df[df["h_pile"] != 0]
        for i in rows.index:
            x = float(df["x_pile"][i])
            tol = 1e-9 * max(1.0, abs(x))
            if abs(x - float(df["x_r"][i])) <= tol or abs(x - float(df["x_l"][i])) <= tol:
                on_boundary.append(ns)
                break
    if not on_boundary:
        return ["no slice count in SLICE_COUNTS puts a pile on a slice boundary, "
                "so the corner-claim case is no longer exercised"]
    print(f"  a pile lands on a slice boundary at {sorted(set(on_boundary))} slices")
    return []


def _mutation(label, apply, restore, leg, fails):
    apply()
    try:
        caught = leg()
    finally:
        restore()
    if not caught:
        fails.append(f"{label}: the leg passed with the defect in place")
    else:
        print(f"  mutation  {label} -> caught ({len(caught)} failure(s))")


def leg_mutations():
    """Claiming a corner crossing by build order must put the asymmetry back."""
    from xslope import slice as xslice
    fails = []
    original = xslice._corner_claim_is_this_slice

    def first_seen_wins(x_cross, x_l, x_r, i, n_slices, right_facing):
        """The rule as it stood: whichever base is built first keeps it."""
        return True

    _mutation("the corner claimed by build order",
              lambda: setattr(xslice, '_corner_claim_is_this_slice', first_seen_wins),
              lambda: setattr(xslice, '_corner_claim_is_this_slice', original),
              leg_the_mirror_pair_agrees, fails)
    return fails


LEGS = [
    ("a pile lands on a slice boundary", leg_a_pile_lands_on_a_boundary),
    ("the mirror pair agrees", leg_the_mirror_pair_agrees),
    ("mutations", leg_mutations),
]


def run():
    failures = []
    for label, fn in LEGS:
        print(f"[{label}]")
        try:
            failures.extend(fn())
        except Exception as e:
            failures.append(f"{label}: raised {type(e).__name__}: {e}")
    return failures


if __name__ == "__main__":
    import sys
    fails = run()
    if fails:
        print("\nFAILURES:")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("\nPASS: a pile model and its mirror image read the same, "
          "at every slice count.")
