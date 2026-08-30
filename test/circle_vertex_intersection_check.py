"""A circle that daylights exactly on a ground-surface vertex.

Ground crossings are solved segment by segment, so a circle passing precisely
through a vertex of the ground surface reports that vertex TWICE — once as the
end root of the segment arriving at it, once as the start root of the segment
leaving it. ``get_sorted_intersections`` then saw three points where there are
two crossings, took its "the circle re-enters beyond the toe" branch, and kept
the last two for a left-facing slope: the two copies of the vertex. The clipped
failure surface came back as a zero-length two-point segment standing at that
one x — drawn as nothing (only the center arrow appeared) and, in a solve, a
degenerate surface that was scored rather than refused.

The trigger on the earth dam is exact arithmetic, not round-off: center (7, 55),
R = 55, crest edge vertex (51, 22), and 44^2 + 33^2 = 55^2. Nudging the center a
tenth of a unit gives the arc the circle should always have had.

The fix collapses coincident crossings to one representative before anything
counts them (``_dedupe_coincident``, tolerance 1e-6 in model units — the
module's geometric tolerance). What is pinned here:

  1. the vertex circle yields exactly two distinct crossings, the full arc from
     the toe at x ~ 0.814 to the vertex at x = 51;
  2. the clipped surface spans that arc instead of collapsing to a point;
  3. FS is continuous through the vertex — the vertex circle and circles a
     thousandth of a unit either side agree to 1e-3, and the tenth-of-a-unit
     nudge from the bug report agrees to well under a percent;
  4. a TANGENT circle, whose single touch point was likewise doubled, is refused
     by the existing reason instead of passing on a zero-length surface;
  5. the six shipped reinforcement circles, which all daylight on a vertex and
     whose duplicate points the old pruning happened to survive, keep the exact
     same two endpoints — the fix moves no locked result;
  6. a mutation: with the dedupe disabled, leg 1 reproduces the original bug.

Run directly:  PYTHONPATH=. python3 test/circle_vertex_intersection_check.py
"""

import contextlib
import copy
import io
import math
import warnings

from xslope import slice as sl
from xslope.fileio import load_slope_data
from xslope.slice import (circle_polyline_intersections, generate_slices,
                          get_sorted_intersections)
from xslope.solve import spencer

DAM = "docs/tutorials/files/xslope_earth_dam_fs_time.xlsx"

# The exact-vertex circle and the vertex it lands on.
VERTEX_CIRCLE = (7.0, 55.0, 55.0)
VERTEX = (51.0, 22.0)
TOE_X = 0.814237          # the other daylight point, on the upstream face
TOE_Y = 0.348959

# Nudges of the center along the R = Yo family (the toe stays on the same face).
FINE_NUDGE = 0.001        # FS must match the vertex circle to FS_TOL
COARSE_NUDGE = 0.1        # the bug report's nudge; a real 0.05-unit longer arc
FS_TOL = 1e-3
FS_TOL_COARSE = 0.02      # 0.4% of FS ~ 5.4; the coarse nudge moves both ends

XY_TOL = 1e-4             # daylight points, model units
COARSE_XY_TOL = 0.06      # the coarse nudge's own geometric offset

# Reinforcement circles that also daylight on a vertex (duplicated toe point).
REINFORCE = [
    ("docs/fem/files/xslope_reinforce_fem.xlsx", 0.0, 40.0, 40.0),
    ("docs/lem/files/xslope_reinforce.xlsx", 0.0, 40.0, 40.0),
    ("docs/lem/files/xslope_reinforce_hw.xlsx", 0.0, 40.0, 40.0),
    ("docs/lem/files/xslope_reinforce_rface.xlsx", 100.0, 40.0, 40.0),
    ("docs/tutorials/files/xslope_reinforced_slope.xlsx", 0.0, 40.0, 40.0),
    ("docs/tutorials/files/xslope_reinforced_slope_start.xlsx", 0.0, 40.0, 40.0),
]


def _quiet_load(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()):
            return load_slope_data(path)


def _dry(path):
    """The dam with every pore-pressure source switched off — the arc geometry is
    the subject here, and u = none needs no mesh or piezometric line."""
    d = _quiet_load(path)
    for m in d["materials"]:
        m["u"] = "none"
    d["piezo_line"] = []
    d["piezo_line2"] = []
    d["seep_u"] = None
    return d


def _circle(Xo, Yo, R):
    return {"Xo": Xo, "Yo": Yo, "R": R, "Depth": None, "Xi": None, "Yi": None,
            "Option": "Circle", "Movement": "Fixed"}


def _xy(points):
    return [(p.x, p.y) for p in points]


def _fmt(points):
    return "[" + ", ".join(f"({x:.6f}, {y:.6f})" for x, y in _xy(points)) + "]"


def _solve(base, Xo, Yo, R, num_slices=40):
    """Single-surface Spencer on one circle. Returns (FS, clipped_surface, msg)."""
    d = copy.deepcopy(base)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()):
            ok, res = generate_slices(d, circle=_circle(Xo, Yo, R),
                                      num_slices=num_slices)
    if not ok:
        return None, None, str(res)
    slice_df, surface = res
    with contextlib.redirect_stdout(io.StringIO()):
        ok2, r = spencer(slice_df.copy())
    if not ok2:
        return None, surface, str(r)
    return r["FS"], surface, ""


def check_vertex_arc(ground):
    """Leg 1: two distinct crossings, the full arc, not three points."""
    failures = []
    Xo, Yo, R = VERTEX_CIRCLE
    raw = circle_polyline_intersections(Xo, Yo, R, ground)
    print(f"  vertex circle Xo={Xo} Yo={Yo} R={R}")
    print(f"    crossings: {_fmt(raw)}")
    if len(raw) != 2:
        failures.append(f"vertex circle: {len(raw)} crossings reported, expected 2 "
                        f"(a doubled vertex is one crossing) -> {_fmt(raw)}")
        return failures

    ok, msg, pruned = get_sorted_intersections(None, ground,
                                               circle_params={"Xo": Xo, "Yo": Yo, "R": R})
    if not ok:
        failures.append(f"vertex circle: get_sorted_intersections refused it ({msg})")
        return failures
    print(f"    kept:      {_fmt(pruned)}")
    want = [(TOE_X, TOE_Y), VERTEX]
    for (gx, gy), (wx, wy), end in zip(_xy(pruned), want, ("left", "right")):
        if abs(gx - wx) > XY_TOL or abs(gy - wy) > XY_TOL:
            failures.append(f"vertex circle {end} end at ({gx:.6f}, {gy:.6f}), "
                            f"expected ({wx:.6f}, {wy:.6f})")
    return failures


def check_nudged_arc(ground):
    """Leg 1b: the nudged circle gives the same arc, offset only by the nudge."""
    failures = []
    Xo, Yo, R = VERTEX_CIRCLE
    e = COARSE_NUDGE
    nudged = circle_polyline_intersections(Xo, Yo + e, R + e, ground)
    print(f"  nudged circle Xo={Xo} Yo={Yo + e} R={R + e}")
    print(f"    crossings: {_fmt(nudged)}")
    if len(nudged) != 2:
        failures.append(f"nudged circle: {len(nudged)} crossings, expected 2")
        return failures
    for (gx, gy), (wx, wy), end in zip(_xy(nudged), [(TOE_X, TOE_Y), VERTEX],
                                       ("left", "right")):
        if abs(gx - wx) > COARSE_XY_TOL:
            failures.append(f"nudged circle {end} end x={gx:.6f} is more than "
                            f"{COARSE_XY_TOL} from the vertex circle's {wx:.6f}")
    return failures


def check_surface_not_collapsed(base):
    """Leg 2: the clipped surface spans the arc instead of standing at one x."""
    failures = []
    Xo, Yo, R = VERTEX_CIRCLE
    fs, surface, msg = _solve(base, Xo, Yo, R)
    if surface is None:
        failures.append(f"vertex circle produced no surface ({msg})")
        return failures
    xmin, _, xmax, _ = surface.bounds
    span = xmax - xmin
    print(f"    clipped surface: x {xmin:.4f} -> {xmax:.4f} (span {span:.4f}), "
          f"length {surface.length:.4f}, {len(surface.coords)} points")
    if span < 45.0:
        failures.append(f"clipped surface spans only {span:.6f} in x — the arc from "
                        f"{TOE_X:.3f} to {VERTEX[0]} is ~50 units wide")
    if len(surface.coords) < 10:
        failures.append(f"clipped surface has {len(surface.coords)} points — a "
                        f"discretized arc has many")
    return failures


def check_fs_continuity(base):
    """Leg 3: FS is continuous through the vertex circle."""
    failures = []
    Xo, Yo, R = VERTEX_CIRCLE
    fs0, _, msg0 = _solve(base, Xo, Yo, R)
    if fs0 is None:
        failures.append(f"vertex circle did not solve ({msg0})")
        return failures
    print(f"    vertex circle          FS = {fs0:.6f}")
    for e, tol in ((FINE_NUDGE, FS_TOL), (-FINE_NUDGE, FS_TOL),
                   (COARSE_NUDGE, FS_TOL_COARSE)):
        fs, _, msg = _solve(base, Xo, Yo + e, R + e)
        if fs is None:
            failures.append(f"nudge {e:+g} did not solve ({msg})")
            continue
        d = abs(fs - fs0)
        flag = "ok" if d <= tol else "OFF"
        print(f"    nudge {e:+7.4f}          FS = {fs:.6f}  |d| = {d:.6f} "
              f"(tol {tol}) {flag}")
        if d > tol:
            failures.append(f"nudge {e:+g}: FS {fs:.6f} differs from the vertex "
                            f"circle's {fs0:.6f} by {d:.6f} > {tol}")
    return failures


def check_tangent_refused(ground):
    """Leg 4: a tangent circle is refused, not passed on as a zero-length surface."""
    failures = []
    # Tangent from above to the upstream face, which runs (0, 0) -> (42, 18).
    Xo, Yo = 30.0, 20.0
    R = abs(18 * Xo - 42 * Yo) / math.hypot(18, 42)
    pts = circle_polyline_intersections(Xo, Yo, R, ground)
    ok, msg, _ = get_sorted_intersections(None, ground,
                                          circle_params={"Xo": Xo, "Yo": Yo, "R": R})
    print(f"  tangent circle Xo={Xo} Yo={Yo} R={R:.6f}")
    print(f"    crossings: {_fmt(pts)}  ->  success={ok} msg={msg!r}")
    if len(pts) != 1:
        failures.append(f"tangent circle: {len(pts)} crossings, expected 1 touch point")
    if ok:
        failures.append("tangent circle was accepted — a single touch point is not "
                        "a failure surface and must be refused")
    elif "intersection points" not in msg:
        failures.append(f"tangent circle refused without the usual reason: {msg!r}")
    return failures


def check_reinforce_unmoved():
    """Leg 5: the shipped vertex-daylighting circles keep their exact endpoints."""
    failures = []
    for path, Xo, Yo, R in REINFORCE:
        d = _quiet_load(path)
        ground = d["ground_surface"]
        ok, msg, pruned = get_sorted_intersections(
            None, ground, circle_params={"Xo": Xo, "Yo": Yo, "R": R})
        if not ok:
            failures.append(f"{path}: its circle is now refused ({msg})")
            continue
        # Reproduce what the pruning returned before the dedupe existed.
        old = _legacy_pruned(Xo, Yo, R, ground)
        same = (old is not None and len(old) == len(pruned) and
                all(abs(a[0] - b[0]) <= 1e-12 and abs(a[1] - b[1]) <= 1e-12
                    for a, b in zip(_xy(old), _xy(pruned))))
        print(f"    {path.split('/')[-1]:38s} {_fmt(pruned)} "
              f"{'unmoved' if same else 'MOVED'}")
        if not same:
            failures.append(f"{path}: endpoints moved, was {_fmt(old) if old else None} "
                            f"now {_fmt(pruned)}")
    return failures


def _legacy_pruned(Xo, Yo, R, ground):
    """The pre-fix result: crossings with duplicates, then the old count-and-prune."""
    with _dedupe_disabled():
        pts = sorted(circle_polyline_intersections(Xo, Yo, R, ground), key=lambda p: p.x)
    if len(pts) < 2:
        return None
    if len(pts) == 2:
        return pts
    keep = pts[:2] if pts[0].y > pts[-1].y else pts[-2:]
    return sorted(keep, key=lambda p: p.x)


@contextlib.contextmanager
def _dedupe_disabled():
    """Mutation harness: put the module back the way it was before the fix."""
    original = sl._dedupe_coincident
    sl._dedupe_coincident = lambda points, tol=1e-6: list(points)
    try:
        yield
    finally:
        sl._dedupe_coincident = original


def check_mutation(ground):
    """Leg 6: without the dedupe, the vertex circle reproduces the original bug."""
    failures = []
    Xo, Yo, R = VERTEX_CIRCLE
    with _dedupe_disabled():
        raw = circle_polyline_intersections(Xo, Yo, R, ground)
        ok, msg, pruned = get_sorted_intersections(
            None, ground, circle_params={"Xo": Xo, "Yo": Yo, "R": R})
    print(f"    dedupe off: crossings {_fmt(raw)}")
    print(f"    dedupe off: kept      {_fmt(pruned) if pruned else None}")
    if len(raw) != 3:
        failures.append(f"mutation: expected the doubled vertex (3 crossings), got "
                        f"{len(raw)} — this check no longer exercises the bug")
    if pruned is None or len(pruned) != 2:
        failures.append("mutation: expected two points kept")
    elif abs(pruned[0].x - pruned[1].x) > 1e-9:
        failures.append(f"mutation: expected the collapsed pair at one x, got "
                        f"{_fmt(pruned)} — this check no longer exercises the bug")
    return failures


def run():
    failures = []
    base = _dry(DAM)
    ground = base["ground_surface"]

    print("Ground-surface vertex circle (earth dam, crest edge at (51, 22)):")
    failures += check_vertex_arc(ground)
    failures += check_nudged_arc(ground)

    print("  clipped surface must span the arc:")
    failures += check_surface_not_collapsed(base)

    print("  FS continuity through the vertex (single-surface Spencer, u = none):")
    failures += check_fs_continuity(base)

    print("  tangent circle must be refused:")
    failures += check_tangent_refused(ground)

    print("  shipped vertex-daylighting circles must not move:")
    failures += check_reinforce_unmoved()

    print("  mutation — the bug returns with the dedupe disabled:")
    failures += check_mutation(ground)

    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll vertex-circle intersection checks passed.")


if __name__ == "__main__":
    main()
