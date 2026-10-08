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

A piles sheet carries no force angle, so the angle a solver applies is derived
from the pile's end points: the force is perpendicular to the pile, the higher
end is the head whichever end is entered first, and the angle is measured in
the resisting frame — positive when the tip lies upslope of the head — so it
depends on which way the slope faces. The legs that state no angle send the
model through a workbook and back, so the angle is the one the loader and the
slice builder derive, as for any file a user writes:
  - a battered pile and its mirror image read the same, with the force tilted
    upward when the tip lies upslope of the head and downward when it lies
    downslope, on both facings;
  - a pile entered tip first reads the same as one entered head first;
  - a vertical pile's force is horizontal (angle 0) in either entry order and on
    either facing.

Run directly:  PYTHONPATH=. python3 test/pile_symmetry_check.py
"""

import copy
import os
import tempfile
from math import atan, degrees

import numpy as np
from shapely.affinity import scale
from shapely.geometry import LineString

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
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


#: The force a battered pile carries is stated: the Ito & Matsui computation is
#: for vertical piles only and refuses a battered one.
H_STATED = 5000.0
#: How far the tip of a battered pile lies from below its head, horizontally.
BATTER = 4.0
#: Two solves of identical slices, one with its pile entered tip first: the same
#: arithmetic in a different order.
ORDER_TOL_REL = 1e-9


def _through_a_workbook(d, folder, name):
    """Write a model to a workbook and read it back. The piles sheet has no force
    angle, so the angle the solvers apply is the one derived from the end points,
    as it is for any file a user writes."""
    return load_slope_data(save_slope_data_to_xlsx(d, os.path.join(folder,
                                                                   name + ".xlsx")))


def _with_piles(d, offset, tip_first=False):
    """The model with every pile's tip ``offset`` UPSLOPE of its head (negative =
    downslope), its force stated, entered head first or tip first. The sample
    descends to the left, so upslope is +x; its mirror image carries the same
    piles reflected, so their tips lie upslope there too."""
    m = copy.deepcopy(d)
    for p in m["pile_lines"]:
        p["H"] = H_STATED
        p["x2"] = p["x1"] + offset
        if tip_first:
            p["x1"], p["y1"], p["x2"], p["y2"] = p["x2"], p["y2"], p["x1"], p["y1"]
    return m


def _facing(sdf):
    return bool(sdf["y_lt"].iat[0] > sdf["y_rt"].iat[-1])


def _applied_angles(sdf):
    """The force angle (degrees) the solvers read on each slice a pile crosses."""
    rows = sdf[sdf["h_pile"] != 0]
    return sorted(round(degrees(t), 9) for t in rows["theta_p"])


def _expected_angles(d, offset):
    """atan(d_u / (y_head - y_tip)) for each of the model's piles."""
    return sorted(round(degrees(atan(offset / (max(p["y1"], p["y2"])
                                               - min(p["y1"], p["y2"])))), 9)
                  for p in d["pile_lines"])


def leg_a_battered_pile_reads_the_same_mirrored():
    """A battered pile, its angle derived from its end points, and its mirror
    image: every method, at every slice count, with the force tilted upward when
    the tip lies upslope of the head and downward when it lies downslope."""
    failures = []
    base = load_slope_data(PILE_MODEL)
    with tempfile.TemporaryDirectory() as tmp:
        for offset, word in ((BATTER, "upslope"), (-BATTER, "downslope")):
            m = _with_piles(base, offset)
            left = _through_a_workbook(m, tmp, "left")
            right = _through_a_workbook(mirror_slope_data(m), tmp, "right")
            want = _expected_angles(m, offset)
            print(f"  battered pile, tip {BATTER:g} {word} of the head "
                  f"(theta_p {', '.join(f'{a:+.2f}' for a in want)}):")
            for ns in SLICE_COUNTS:
                fl, sdf_l = _solve_all(left, ns)
                fr, sdf_r = _solve_all(right, ns)
                if not (not _facing(sdf_l) and _facing(sdf_r)):
                    failures.append(f"tip {word}: expected the model left-facing and "
                                    f"its mirror right-facing")
                for side, sdf in (("model", sdf_l), ("mirror", sdf_r)):
                    got = _applied_angles(sdf)
                    if got != want:
                        failures.append(f"tip {word}, {ns} slices, {side}: the "
                                        f"solvers read theta_p {got}, the pile's "
                                        f"inclination is {want}")
                worst = (0.0, None)
                for name, _ in METHODS:
                    a, b = fl[name], fr[name]
                    if a is None or b is None:
                        failures.append(f"tip {word}, {name}: solve failed "
                                        f"(model={a}, mirror={b})")
                        continue
                    asym = abs(a - b) / a * 100
                    if asym > worst[0]:
                        worst = (asym, name)
                    if asym >= TOL_PCT:
                        failures.append(f"tip {word}, {name} at {ns} slices: asym "
                                        f"{asym:.3f}% >= {TOL_PCT}% (model={a:.5f}, "
                                        f"mirror={b:.5f})")
                if worst[1] is not None:
                    print(f"    {ns:3d} slices: worst {worst[1]} {worst[0]:.4f}%")
    return failures


def leg_the_entry_order_does_not_matter():
    """A pile entered tip first reads the same as the same pile entered head
    first: vertical and battered both ways, on the model and its mirror image."""
    failures = []
    base = load_slope_data(PILE_MODEL)
    ns = 40
    with tempfile.TemporaryDirectory() as tmp:
        for offset in (0.0, BATTER, -BATTER):
            for side in ("model", "mirror"):
                runs = []
                for tip_first in (False, True):
                    m = _with_piles(base, offset, tip_first)
                    if side == "mirror":
                        m = mirror_slope_data(m)
                    runs.append(_solve_all(_through_a_workbook(m, tmp, "order"), ns)[0])
                head, tip = runs
                worst = 0.0
                for name, _ in METHODS:
                    a, b = head[name], tip[name]
                    if a is None or b is None:
                        failures.append(f"offset {offset:+g}, {side}, {name}: solve "
                                        f"failed (head first={a}, tip first={b})")
                        continue
                    rel = abs(a - b) / a
                    worst = max(worst, rel)
                    if rel > ORDER_TOL_REL:
                        failures.append(f"offset {offset:+g}, {side}, {name}: entered "
                                        f"tip first {b:.5f}, head first {a:.5f}")
                print(f"  tip offset {offset:+g}, {side:6s}: largest difference "
                      f"{worst:.1e} (spencer {head['spencer']:.4f} head first, "
                      f"{tip['spencer']:.4f} tip first)")
    return failures


def leg_a_vertical_pile_pushes_horizontally():
    """A vertical pile's force angle is 0 as the file is read, in either entry
    order, and the solvers apply 0 on both facings: its force is the horizontal
    force it has always been."""
    failures = []
    base = load_slope_data(PILE_MODEL)
    with tempfile.TemporaryDirectory() as tmp:
        for tip_first in (False, True):
            m = _with_piles(base, 0.0, tip_first)
            for side, d in (("model", m), ("mirror", mirror_slope_data(m))):
                loaded = _through_a_workbook(d, tmp, "vertical")
                order = "tip first" if tip_first else "head first"
                stored = [p["theta_p"] for p in loaded["pile_lines"]]
                if stored != [0.0] * len(stored):
                    failures.append(f"vertical, {order}, {side}: the loader stored "
                                    f"theta_p {stored}, not 0")
                _fs, sdf = _solve_all(loaded, 40)
                got = _applied_angles(sdf)
                if not got or any(a != 0.0 for a in got):
                    failures.append(f"vertical, {order}, {side}: the solvers read "
                                    f"theta_p {got}, not 0")
                print(f"  vertical, {order:10s}, {side:6s}: stored {stored}, "
                      f"applied {got}")
    return failures


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

    # The pile force angle: the loader and the slice builder both take it from
    # slice.pile_force_angle, so replacing that function replaces the rule in
    # both places.
    angle = xslice.pile_force_angle

    def upslope_is_always_plus_x(x1, y1, x2, y2, right_facing=None):
        """Perpendicular to the pile with the higher end as the head, but the
        tip's offset measured along +x whichever way the slope faces."""
        if y2 > y1:
            x1, y1, x2, y2 = x2, y2, x1, y1
        return degrees(np.arctan2(x2 - x1, y1 - y2))

    def the_first_end_is_the_head(x1, y1, x2, y2, right_facing=None):
        """Measured in the resisting frame, but with the first end entered taken
        as the head, whichever end is higher."""
        if x2 == x1:
            return 0.0 if y1 > y2 else 180.0
        if right_facing is None:
            return None
        d_u = -(x2 - x1) if right_facing else (x2 - x1)
        return degrees(np.arctan2(d_u, y1 - y2))

    for label, mutant, leg in (
            ("the pile's offset read along +x on either facing",
             upslope_is_always_plus_x, leg_a_battered_pile_reads_the_same_mirrored),
            ("the first end entered taken as the head",
             the_first_end_is_the_head, leg_the_entry_order_does_not_matter),
            ("the first end entered taken as the head (vertical)",
             the_first_end_is_the_head, leg_a_vertical_pile_pushes_horizontally)):
        _mutation(label,
                  lambda m=mutant: setattr(xslice, 'pile_force_angle', m),
                  lambda: setattr(xslice, 'pile_force_angle', angle),
                  leg, fails)
    return fails


LEGS = [
    ("a pile lands on a slice boundary", leg_a_pile_lands_on_a_boundary),
    ("the mirror pair agrees", leg_the_mirror_pair_agrees),
    ("a battered pile reads the same mirrored",
     leg_a_battered_pile_reads_the_same_mirrored),
    ("the entry order does not matter", leg_the_entry_order_does_not_matter),
    ("a vertical pile pushes horizontally", leg_a_vertical_pile_pushes_horizontally),
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
