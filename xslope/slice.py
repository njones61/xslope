# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from math import sin, cos, tan, radians, atan, atan2, degrees, sqrt

import warnings

import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection

from .mesh import find_element_containing_point, interpolate_at_point
from .water import seep_water_table_profile
from .hoekbrown import hb_tangent

_ito_matsui_warned = False  # module-level flag: warn once about large H


def _support_choice(row, key, vocabulary, default, index, what):
    """One structural-row vocabulary field (Dir, Appl), or a loud refusal.

    A blank field -- absent, None, NaN or empty -- takes the documented default,
    which is what a blank cell means on the sheet. Anything else must be a word
    from the vocabulary.

    The engine used to test these with ``== "tangent"`` and ``== "active"`` and
    take the other branch for everything else, so an entry that is not a
    direction at all was applied AXIALLY and an entry that is not an application
    was applied as PASSIVE capacity. A ``Dir = 0.0`` written into the sheet came
    back as a different, entirely valid model -- a different force, at a
    different inclination, with no sign anywhere that the input had been
    reinterpreted. The loader and the preflight both refuse these values; when
    they are bypassed (a model assembled in memory, an importer, a sweep) the
    engine states the problem rather than guessing at it.
    """
    raw = row.get(key, None)
    if raw is None or (isinstance(raw, float) and raw != raw):
        return default
    value = str(raw).strip().lower()
    if not value:
        return default
    if value not in vocabulary:
        label = row.get("label") or f"{what} {index + 1}"
        raise ValueError(
            f"{what} '{label}' (row {index + 1}) has {key.capitalize()} = "
            f"{raw!r}, which is not one of: {', '.join(vocabulary)}. Fix the "
            f"value on the sheet; the engine will not guess which one was meant.")
    return value


def get_circular_y_coordinates(x_coords, Xo, Yo, R):
    """
    Calculate y-coordinates on a circular failure surface for given x-coordinates.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        Xo, Yo (float): Center coordinates of the circle
        R (float): Radius of the circle
        
    Returns:
        numpy.ndarray: Y-coordinates on the circle (bottom half)
    """
    x_coords = np.asarray(x_coords)
    # Calculate y-coordinates for the bottom half of the circle
    # y = Yo - sqrt(R^2 - (x - Xo)^2)
    dx_squared = (x_coords - Xo) ** 2
    # Handle cases where x is outside the circle
    valid_mask = dx_squared <= R ** 2
    y_coords = np.full_like(x_coords, np.nan)
    y_coords[valid_mask] = Yo - np.sqrt(R ** 2 - dx_squared[valid_mask])
    return y_coords


def get_circular_intersection_points(ground_surface, Xo, Yo, R, x_min, x_max):
    """
    Find intersection points between a circular failure surface and ground surface.
    
    Parameters:
        ground_surface (LineString): Ground surface geometry
        Xo, Yo, R (float): Circle parameters
        x_min, x_max (float): X-range to search
        
    Returns:
        tuple: (x_min, x_max, y_left, y_right, success)
    """
    # Create a dense set of points on the circle for intersection testing
    x_test = np.linspace(x_min, x_max, 1000)
    y_circle = get_circular_y_coordinates(x_test, Xo, Yo, R)
    
    # Create circle line for intersection
    valid_mask = ~np.isnan(y_circle)
    if not np.any(valid_mask):
        return None, None, None, None, False
    
    circle_coords = list(zip(x_test[valid_mask], y_circle[valid_mask]))
    circle_line = LineString(circle_coords)
    
    # Find intersections
    intersections = circle_line.intersection(ground_surface)
    
    if isinstance(intersections, Point):
        points = [intersections]
    elif isinstance(intersections, MultiPoint):
        points = list(intersections.geoms)
    elif isinstance(intersections, GeometryCollection):
        points = [g for g in intersections.geoms if isinstance(g, Point)]
    else:
        points = []
    
    if len(points) < 2:
        return None, None, None, None, False
    
    # Sort by x and take the two endpoints
    points = sorted(points, key=lambda p: p.x)
    x_min, x_max = points[0].x, points[-1].x
    y_left, y_right = points[0].y, points[-1].y
    
    return x_min, x_max, y_left, y_right, True


def get_ground_surface_y_coordinates(x_coords, ground_surface):
    """
    Get y-coordinates on the ground surface for given x-coordinates using interpolation.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        ground_surface (LineString): Ground surface geometry
        
    Returns:
        numpy.ndarray: Y-coordinates on the ground surface
    """
    x_coords = np.asarray(x_coords)
    ground_coords = np.array(ground_surface.coords)
    ground_x = ground_coords[:, 0]
    ground_y = ground_coords[:, 1]
    
    # Sort by x to ensure proper interpolation
    sort_idx = np.argsort(ground_x)
    ground_x = ground_x[sort_idx]
    ground_y = ground_y[sort_idx]
    
    # Interpolate y-coordinates
    y_coords = np.interp(x_coords, ground_x, ground_y, left=np.nan, right=np.nan)
    return y_coords


def get_piezometric_y_coordinates(x_coords, piezo_line):
    """
    Get y-coordinates on the piezometric surface for given x-coordinates.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        piezo_line (list): Piezometric line coordinates
        
    Returns:
        numpy.ndarray: Y-coordinates on the piezometric surface
    """
    if not piezo_line:
        return np.full_like(x_coords, np.nan)
    
    x_coords = np.asarray(x_coords)
    piezo_coords = np.array(piezo_line)
    piezo_x = piezo_coords[:, 0]
    piezo_y = piezo_coords[:, 1]
    
    # Sort by x to ensure proper interpolation
    sort_idx = np.argsort(piezo_x)
    piezo_x = piezo_x[sort_idx]
    piezo_y = piezo_y[sort_idx]
    
    # Interpolate y-coordinates
    y_coords = np.interp(x_coords, piezo_x, piezo_y, left=np.nan, right=np.nan)
    return y_coords


def _dedupe_coincident(points, tol=1e-6):
    """
    Collapse runs of coincident points to one representative, preserving order.

    Crossings are found segment by segment, so a circle that passes exactly
    through a polyline VERTEX produces that vertex twice — once as the end root
    of the segment arriving at it, once as the start root of the segment leaving
    it. A tangent circle likewise yields its single touch point twice, from the
    two equal roots of one segment. Both are one geometric crossing, and callers
    that count crossings must not see them as two.

    ``tol`` is the module's geometric tolerance in model units (see
    _resolve_right_facing and the clip-endpoint test in generate_failure_surface);
    two genuinely distinct daylight points a micron apart would be a degenerate
    sliver, not a surface worth slicing.
    """
    kept = []
    for p in points:
        if any(abs(p.x - q.x) <= tol and abs(p.y - q.y) <= tol for q in kept):
            continue
        kept.append(p)
    return kept


def _circle_polyline_all(Xo, Yo, R, polyline):
    """
    Every crossing of a full circle with a polyline, one point per crossing.

    Coincident roots from adjoining segments (a vertex hit) or from one segment
    (a tangency) are collapsed by _dedupe_coincident, so the length of the
    returned list is a crossing count callers can trust.

    Callers almost always want circle_polyline_intersections instead — see its
    docstring for why crossings above the equator must not reach the slicer.
    """
    intersections = []
    coords = list(polyline.coords)
    for i in range(len(coords) - 1):
        x1, y1 = coords[i]
        x2, y2 = coords[i+1]
        dx = x2 - x1
        dy = y2 - y1

        # Quadratic coefficients for t
        a = dx**2 + dy**2
        if a == 0:
            continue  # zero-length segment
        b = 2 * (dx * (x1 - Xo) + dy * (y1 - Yo))
        c = (x1 - Xo)**2 + (y1 - Yo)**2 - R**2

        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            continue  # No intersection

        sqrt_disc = np.sqrt(discriminant)
        for sign in [-1, 1]:
            t = (-b + sign * sqrt_disc) / (2 * a)
            if 0 <= t <= 1:
                intersections.append(Point(x1 + t * dx, y1 + t * dy))
    return _dedupe_coincident(intersections)


def circle_polyline_intersections(Xo, Yo, R, polyline):
    """
    Find intersection points between the bottom half of a circle and a polyline (LineString).
    Returns a list of shapely Point objects.

    The `yi < Yo` test is load-bearing, not an optimization. A crossing above the
    equator means the center sits below the ground surface, so the daylight points
    bound an arc longer than a semicircle — reverse curvature. The arc built by
    generate_failure_surface is the bottom semicircle, so such a circle cannot be
    represented: clipping to the daylight x-range would splice a vertical face from
    the daylight point down to the bottom arc, inventing an arbitrary-depth tension
    crack tied to no input and artificially lowering FS. Dropping these points makes
    the circle read as "never reaches the ground" and it is rejected instead, which
    is the correct outcome for a search. Removing this test drops vp023 bishop
    1.130 -> 0.820 and five other locked searches.

    A circle rejected here can still be scored if the input specifies a tension
    crack — see _recover_ends_via_tcrack.
    """
    return [p for p in _circle_polyline_all(Xo, Yo, R, polyline) if p.y < Yo]


def crossings_above_center(Xo, Yo, R, polyline):
    """The ground crossings ``circle_polyline_intersections`` throws away.

    The complement of that function's filter, and the reason a circle can be
    reported as "never reaches the ground" while a plot of it plainly shows it
    cutting the section. Kept here so the slicer's error message, the preflight
    rule and the search all name the same points from one rule rather than three
    re-derivations of it.
    """
    return [p for p in _circle_polyline_all(Xo, Yo, R, polyline) if p.y >= Yo]


def _recover_ends_via_tcrack(ground_surface, circle, tcrack_depth):
    """
    Recover the two daylight points for a reverse-curvature circle whose uphill end
    is resolved by an explicit tension crack. Returns [left, right] or None.

    A tension crack lowers the effective ground surface, but only upstream: the arc
    exits at the crack on the uphill side, while the toe still daylights on the true
    ground (soil cracks in tension behind the crest, not at the toe in compression).
    So each end is taken from its own surface. Lowering the whole surface instead
    drags the toe down with it and invents a crack there — worth ~8% of FS on VP30,
    and drifting with crack depth.

    The principled depth is y_uphill - Yo, which puts the uphill exit exactly on the
    equator; the tolerance below admits that case, which a strict `<` would reject
    while a hair deeper succeeded.

    Only reachable when tcrack_depth > 0, i.e. when the user asked for a crack.
    Searches set no crack, so they still reject these circles outright.
    """
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    crack_surface = LineString([(x, y - tcrack_depth) for x, y in ground_surface.coords])
    cpts = sorted((p for p in _circle_polyline_all(Xo, Yo, R, crack_surface)
                   if p.y <= Yo + 1e-9), key=lambda p: p.x)
    gpts = sorted(circle_polyline_intersections(Xo, Yo, R, ground_surface), key=lambda p: p.x)
    if len(cpts) < 2 or not gpts:
        return None
    right_facing = cpts[0].y > cpts[-1].y
    if right_facing:
        return [cpts[0], gpts[-1]]   # uphill exit at the crack, toe on the true ground
    return [gpts[0], cpts[-1]]


def get_sorted_intersections(failure_surface, ground_surface, circle_params=None):
    """
    Find and sort the intersection points between the failure and ground surfaces,
    pruning extras if the circle exits and re-enters the ground beyond the toe.
    If circle_params is provided, use analytic intersection.
    Returns:
        success (bool), msg (str), points (list of shapely Point)
    """
    if circle_params is not None:
        Xo, Yo, R = circle_params['Xo'], circle_params['Yo'], circle_params['R']
        points = circle_polyline_intersections(Xo, Yo, R, ground_surface)
    else:
        intersections = failure_surface.intersection(ground_surface)
        if isinstance(intersections, MultiPoint):
            points = list(intersections.geoms)
        elif isinstance(intersections, Point):
            points = [intersections]
        elif isinstance(intersections, GeometryCollection):
            points = [g for g in intersections.geoms if isinstance(g, Point)]
        else:
            points = []

    # Collapse coincident crossings BEFORE counting or pruning. A circle through a
    # ground-surface vertex reports it once per adjoining segment; unmerged, that
    # third point sent the pair-picking below to two copies of the same vertex and
    # the clipped surface came back as a zero-length segment instead of the arc.
    # Shapely's intersection on the non-circular path can double a vertex the same
    # way. After this, a surface that still has fewer than two DISTINCT crossings is
    # rejected below by the usual reason, as a tangent or grazing circle should be.
    points = _dedupe_coincident(points)

    # need at least two
    if len(points) < 2:
        msg = f"Expected at least 2 intersection points, but got {len(points)}."
        # A circle can lose a crossing two ways, and only one of them is about the
        # circle being too small or too shallow. The other is the equator filter in
        # circle_polyline_intersections: a crossing ABOVE the center is discarded,
        # because the arc joining it to the toe is longer than a semicircle and the
        # bottom-semicircle surface built below cannot represent it. Left unnamed,
        # that arrives as a count and sends the reader to the ground surface, which
        # is not where the fault is.
        if circle_params is not None:
            above = crossings_above_center(Xo, Yo, R, ground_surface)
            if above:
                p = max(above, key=lambda q: q.y)
                msg += (f" The circle (Xo={Xo:g}, Yo={Yo:g}, R={R:g}) meets the "
                        f"ground at ({p.x:.4g}, {p.y:.4g}), above its own center at "
                        f"y={Yo:g}: the arc would rise above its center on that "
                        f"side, which this slicer does not build. Lower the center "
                        f"or enlarge the radius so both crossings sit below it.")
        return False, msg, None

    # sort by x
    points = sorted(points, key=lambda p: p.x)

    # if exactly two, we're done. A flat arc (both crossings at the same
    # elevation, e.g. a level-ground bearing-capacity prism) is NO LONGER
    # rejected here: facing for it is resolved downstream from the surface's
    # asymmetry (see _resolve_right_facing). For a normal slope this branch is
    # unchanged — it returns the same two points as before.
    if len(points) == 2:
        return True, "", points

    # more than two: decide which pair to keep. On a normal slope the higher end
    # marks the crest side; on a flat arc (y_first == y_last) this falls to the
    # else branch (last two), and the facing is settled downstream as above.
    y_first, y_last = points[0].y, points[-1].y
    if y_first > y_last:
        # right-facing: keep first two
        pruned = points[:2]
    else:
        # left-facing: keep last two
        pruned = points[-2:]

    # sort those two again by x (just in case)
    pruned = sorted(pruned, key=lambda p: p.x)
    return True, "", pruned


def _resolve_right_facing(y_left, y_right, surface_coords, x_min, x_max,
                          override=None, tol=1e-6):
    """Decide whether a failure surface is right-facing.

    ``right_facing`` is the master sliding-direction flag consumed all over
    ``generate_slices``/``solve``. The historical, and still normal-slope,
    definition is ``right_facing = y_left > y_right``: when the LEFT ground
    crossing sits higher than the right one, the crest is on the left, the toe
    on the right, and the mass slides to the RIGHT (``True``). Downstream this
    flips the sign of the slice base angle ``alpha`` (and ``beta``, the tension-
    crack side, and line-load horizontals), so getting it wrong mirrors the
    whole mechanism.

    A *flat arc* — both ground crossings at the same elevation, within ``tol``
    — makes ``y_left > y_right`` undecidable and was previously rejected
    outright. That is collateral damage for legitimate level-ground bearing-
    capacity mechanisms (e.g. a Prandtl prism under a footing). For that case
    only, facing is inferred from the surface's ASYMMETRY: the mass slides
    toward the crossing nearest the deepest point of the surface.

        deepest point LEFT  of the crossings' midpoint -> slides left  -> left-facing  (right_facing=False)
        deepest point RIGHT of the crossings' midpoint -> slides right -> right-facing (right_facing=True)

    A perfectly symmetric flat arc (deepest point at the midpoint) carries no
    asymmetry to read; with no override it deterministically defaults to
    left-facing (``right_facing=False``, which also matches what ``solve``'s own
    ``y_lb[0] > y_rb[-1]`` recomputation yields for equal-elevation ends) and
    returns a note so the caller can surface it.

    ``override`` (a truthy/falsy value, or ``None`` for "auto") wins whenever it
    is not ``None`` and suppresses the note — the caller asked for a specific
    facing and should not also be warned.

    IMPORTANT: for a normal slope (``abs(y_left - y_right) >= tol``) this returns
    exactly ``y_left > y_right`` — byte-identical to the prior inline logic. The
    asymmetry branch fires ONLY in the equal-elevation case that used to be
    rejected. ``tol`` matches the flat-arc tolerance in ``get_sorted_intersections``.

    Returns ``(right_facing: bool, note: str or None)``.
    """
    if override is not None:
        return bool(override), None

    if abs(y_left - y_right) >= tol:
        # NORMAL slope: unchanged historical behavior.
        return (y_left > y_right), None

    # FLAT ARC: infer facing from surface asymmetry.
    x_mid = 0.5 * (x_min + x_max)
    interior = [(x, y) for (x, y) in surface_coords if x_min < x < x_max]
    pts = interior if interior else list(surface_coords)
    x_deep = min(pts, key=lambda p: p[1])[0]
    if x_deep < x_mid - tol:
        return False, None
    if x_deep > x_mid + tol:
        return True, None
    note = ("Flat-arc failure surface is symmetric about its crossings' "
            "midpoint; facing is indeterminate. Defaulted to left-facing "
            "(right_facing=False). Pass slope_data['right_facing'] or the "
            "generate_slices right_facing= kwarg to override.")
    return False, note


def generate_failure_surface(ground_surface, circular, circle=None, non_circ=None, tcrack_depth=0,
                             right_facing=None):
    """
    Generates a failure surface based on either a circular or non-circular definition.

    Parameters:
        ground_surface (LineString): The ground surface geometry.
        circular (bool): Whether to use circular failure surface.
        circle (dict, optional): Dictionary with keys 'Xo', 'Yo', 'Depth', and 'R'.
        non_circ (list, optional): List of dicts with keys 'X', 'Y', and 'Movement'.
        tcrack_depth (float, optional): Tension crack depth.

    Returns:
        tuple: (success, result)
            - If success is True:
                result = (x_min, x_max, y_left, y_right, clipped_surface)
            - If success is False:
                result = error message string
    """
    # --- Step 1: Build failure surface ---
    if circular and circle:
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
        theta_range = np.linspace(np.pi, 2 * np.pi, 100)
        arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        failure_coords = arc
        failure_surface = LineString(arc)
    elif non_circ:
        failure_coords = [(pt['X'], pt['Y']) for pt in non_circ]
        failure_surface = LineString(failure_coords)
    else:
        return False, "Either a circular or non-circular failure surface must be provided."

    # --- Step 2: Intersect with original ground surface to determine slope facing and toe ---
    if circular and circle:
        success, msg, points = get_sorted_intersections(failure_surface, ground_surface, circle_params=circle)
    else:
        success, msg, points = get_sorted_intersections(failure_surface, ground_surface)

    # A reverse-curvature circle daylights above its equator, which Step 2 cannot see,
    # so it reads as never reaching the ground. An explicit tension crack lowers the
    # effective ground on the uphill side and can give that end a valid exit; recover
    # the ends from their own surfaces and skip Step 3, which has then already been
    # applied. Gated on tcrack_depth, so searches still reject these circles.
    tcrack_recovered = False
    if not success and circular and circle and tcrack_depth > 0:
        recovered = _recover_ends_via_tcrack(ground_surface, circle, tcrack_depth)
        if recovered is not None:
            points, success, tcrack_recovered = recovered, True, True

    if not success:
        return False, msg

    x_min, x_max = points[0].x, points[1].x
    y_left, y_right = points[0].y, points[1].y
    # Byte-identical to `y_left > y_right` for a normal slope; the asymmetry rule
    # only engages on a flat arc (equal-elevation crossings, formerly rejected).
    right_facing, _facing_note = _resolve_right_facing(
        y_left, y_right, failure_coords, x_min, x_max, override=right_facing)

    # --- Step 3: If tension crack exists, find intersection with tension crack surface ---
    if tcrack_depth > 0 and not tcrack_recovered:
        # Create tension crack surface as parallel offset of entire ground surface
        tcrack_surface = LineString([(x, y - tcrack_depth) for x, y in ground_surface.coords])

        # Find intersection of failure surface with tension crack surface
        if circular and circle:
            tcrack_points = circle_polyline_intersections(Xo, Yo, R, tcrack_surface)
        else:
            tcrack_intersection = failure_surface.intersection(tcrack_surface)
            if isinstance(tcrack_intersection, Point):
                tcrack_points = [tcrack_intersection]
            elif isinstance(tcrack_intersection, MultiPoint):
                tcrack_points = list(tcrack_intersection.geoms)
            elif isinstance(tcrack_intersection, GeometryCollection):
                tcrack_points = [g for g in tcrack_intersection.geoms if isinstance(g, Point)]
            else:
                tcrack_points = []

        if len(tcrack_points) >= 1:
            # Sort tension crack intersection points by x
            tcrack_points = sorted(tcrack_points, key=lambda p: p.x)

            if right_facing:
                # Right-facing slope: tension crack is on the left (upslope)
                # Use leftmost tcrack intersection as new x_min
                new_left = tcrack_points[0]
                x_min = new_left.x
                y_left = new_left.y
            else:
                # Left-facing slope: tension crack is on the right (upslope)
                # Use rightmost tcrack intersection as new x_max
                new_right = tcrack_points[-1]
                x_max = new_right.x
                y_right = new_right.y

    # --- Step 4: Clip the failure surface between intersection x-range ---
    # Filter coordinates within the x-range
    filtered_coords = [pt for pt in failure_coords if x_min <= pt[0] <= x_max]

    # Add the exact intersection points if they're not already in the filtered list
    left_intersection = (x_min, y_left)
    right_intersection = (x_max, y_right)

    # Check if intersection points are already in the filtered list (with tolerance)
    tol = 1e-6
    has_left = any(abs(pt[0] - x_min) < tol and abs(pt[1] - y_left) < tol for pt in filtered_coords)
    has_right = any(abs(pt[0] - x_max) < tol and abs(pt[1] - y_right) < tol for pt in filtered_coords)

    if not has_left:
        filtered_coords.insert(0, left_intersection)
    if not has_right:
        filtered_coords.append(right_intersection)

    # Sort by x-coordinate to ensure proper ordering
    filtered_coords.sort(key=lambda pt: pt[0])

    clipped_surface = LineString(filtered_coords)

    return True, (x_min, x_max, y_left, y_right, clipped_surface)


def update_slice_materials(slice_df, materials):
    """Rewrite the material-derived slice columns for a new materials list, on a
    slice table generated with ``material_breakdown=True`` — the fast path a
    Monte Carlo campaign takes instead of re-slicing an unchanged geometry.

    Updates w, kw, y_cg (from the stored per-band heights and moments, in the
    same band order and arithmetic the generator used) and the base-strength
    columns c/phi/c1/phi1/d/psi (mc and cp options). Raises ValueError when the
    table cannot honestly take the shortcut — no breakdown stored, a base
    material on a stress-dependent envelope (pow/hb), a weight-dependent pore
    pressure (ru), or suction strength — and the caller falls back to a full
    rebuild.
    """
    mb = slice_df.attrs.get('mat_breakdown')
    if mb is None:
        raise ValueError("slice table carries no material breakdown")
    if mb.get('suction'):
        raise ValueError("suction strength depends on u; rebuild")
    base = mb['base_idx']
    if (base < 0).any():
        raise ValueError("a slice has no base material; rebuild")
    for bi in np.unique(base):
        opt = materials[bi].get('option')
        if opt not in ('mc', 'cp'):
            raise ValueError(f"base material option {opt!r} is stress-dependent; rebuild")
    for m in materials:
        if (m.get('u') or 'none') == 'ru':
            raise ValueError("ru pore pressure depends on weight; rebuild")

    gam = [m['gamma'] for m in materials]
    gsat = [m.get('gamma_sat') if m.get('gamma_sat') is not None else m['gamma']
            for m in materials]
    Hm, Hs, HYm, HYs = mb['Hm'], mb['Hs'], mb['HYm'], mb['HYs']
    dx = slice_df['dx'].to_numpy()
    n = len(slice_df)
    w = np.zeros(n)
    sgh = np.zeros(n)
    sghy = np.zeros(n)
    # Accumulate per band in band order — the generator's own summation order,
    # so a realization at the most-likely values reproduces it exactly.
    for p, mi in enumerate(mb['mat_of_poly']):
        g, gs = gam[mi], gsat[mi]
        band = gs * Hs[:, p] + g * Hm[:, p]
        w += band * dx
        sgh += band
        sghy += gs * HYs[:, p] + g * HYm[:, p]
    slice_df['w'] = w
    slice_df['kw'] = abs(mb['k_seismic']) * w
    with np.errstate(invalid='ignore', divide='ignore'):
        slice_df['y_cg'] = np.where(sgh > 0, sghy / np.maximum(sgh, 1e-300), np.nan)

    y_cb = slice_df['y_cb'].to_numpy()
    c = np.empty(n); phi = np.empty(n)
    c1 = np.zeros(n); phi1 = np.zeros(n)
    dcol = np.zeros(n); psi = np.zeros(n)
    for bi in np.unique(base):
        mask = base == bi
        m = materials[bi]
        if m['option'] == 'mc':
            c[mask] = m['c']; phi[mask] = m['phi']
            c1[mask] = m['c']; phi1[mask] = m['phi']
            dcol[mask] = m['d']; psi[mask] = m['psi']
        else:  # cp
            c[mask] = m['c'] + np.maximum(0.0, m['r_elev'] - y_cb[mask]) * m['cp']
            phi[mask] = 0.0
    slice_df['c'] = c; slice_df['phi'] = phi
    slice_df['c1'] = c1; slice_df['phi1'] = phi1
    slice_df['d'] = dcol; slice_df['psi'] = psi
    return slice_df


def get_y_from_intersection(geom):
    """
    Extracts the maximum Y-coordinate from a geometric intersection result.

    This function handles different geometric types resulting from intersections,
    including Point, MultiPoint, LineString, and GeometryCollection. If the input
    geometry is not one of these or is empty, the function returns None.

    Parameters:
        geom (shapely.geometry.base.BaseGeometry): The geometry object from which
            to extract the Y-coordinate(s).

    Returns:
        float or None: The maximum Y-coordinate found in the geometry, or None if not found.
    """
    if isinstance(geom, Point):
        return geom.y
    elif isinstance(geom, MultiPoint):
        return max(pt.y for pt in geom.geoms)
    elif isinstance(geom, LineString):
        return max(y for _, y in geom.coords) if geom.coords else None
    elif isinstance(geom, GeometryCollection):
        pts = [g for g in geom.geoms if isinstance(g, Point)]
        return max(pt.y for pt in pts) if pts else None
    return None


def calc_dload_resultant(x_l, y_lt, x_r, y_rt, qL, qR, dl):
    """
    Compute:
      - D    : total resultant force from a trapezoidal load varying
               linearly from intensity qL at (x_l,y_lt) to qR at (x_r,y_rt).
      - d_x  : x‐coordinate of the resultant's centroid on the top edge
      - d_y  : y‐coordinate of the resultant's centroid on the top edge

    Parameters
    ----------
    x_l, y_lt : float
        Coordinates of the left‐end of the top edge.
    x_r, y_rt : float
        Coordinates of the right‐end of the top edge.
    qL : float
        Load intensity (force per unit length) at (x_l, y_lt).
    qR : float
        Load intensity (force per unit length) at (x_r, y_rt).
    dl : float
        Actual length along the inclined surface.

    Returns
    -------
    D    : float
           Total resultant (area of trapezoid) = ½ (qL + qR) * dl
    d_x  : float
           Global x‐coordinate of the centroid of that trapezoid
    d_y  : float
           Global y‐coordinate of the centroid (lies on the line segment
           between (x_l,y_lt) and (x_r,y_rt))

    Notes
    -----
    1.  If x_r == x_l (zero‐width slice), this will return D=0 and place
        the "centroid" at (x_l, y_lt).
    2.  For a nonzero‐width trapezoid, the horizontal centroid‐offset from
        x_l is:
             x_offset = (x_r – x_l) * ( qL + 2 qR ) / [3 (qL + qR) ]
        provided (qL + qR) ≠ 0.  If qL + qR ≈ 0, it simply places the
        centroid at the midpoint in x.
    3.  The vertical coordinate d_y is found by linear‐interpolation:
          t = x_offset / (x_r – x_l)
          d_y = y_lt + t ·(y_rt – y_lt)

    """
    dx = x_r - x_l

    # 1) Total resultant force (area under trapezoid) using actual length
    D = 0.5 * (qL + qR) * dl

    # 2) Horizontal centroid offset from left end
    sum_q = qL + qR
    if abs(sum_q) < 1e-12:
        # nearly zero trapezoid => centroid at geometric midpoint
        x_offset = dx * 0.5
    else:
        x_offset = dx * (qL + 2.0 * qR) / (3.0 * sum_q)

    # 3) Global x‐coordinate of centroid
    d_x = x_l + x_offset

    # 4) Corresponding y‐coordinate by linear interpolation along top edge
    t = x_offset / dx
    d_y = y_lt + t * (y_rt - y_lt)

    return D, d_x, d_y


def _polygon_bands(polygon):
    """Decompose a material-zone polygon into y-simple bands.

    A vertical line may enter and leave a zone polygon more than once (an
    inclined dam core, a folded seam). Reducing such a polygon to a single
    (top, bottom) envelope pair counts the gap between its bands as if it were
    solid material, double-counting slice weight against the neighboring zones
    that actually fill the gap. Instead, sweep the polygon's vertex x-values
    and pair boundary crossings per strip into trapezoids, then chain
    trapezoids that share a y-interval at strip boundaries into bands. Each
    band IS y-simple, so np.interp on its (top, bottom) tables is exact.

    Returns a list of (txs, tys, bys) numpy arrays, one per band; y-simple
    polygons yield exactly one band, identical to the former envelope tables.
    """
    ring = list(polygon.exterior.coords)
    if len(ring) > 1 and ring[0] == ring[-1]:
        ring = ring[:-1]
    n = len(ring)
    xs = sorted(set(p[0] for p in ring))
    bands = []       # each: {'xs': [...], 'top': [...], 'bot': [...]}
    open_bands = []  # bands whose last sample sits at the current strip's left edge
    for xa, xb in zip(xs[:-1], xs[1:]):
        if xb - xa <= 1e-12:
            continue
        xm = 0.5 * (xa + xb)
        crossings = []
        for i in range(n):
            x1, y1 = ring[i]
            x2, y2 = ring[(i + 1) % n]
            if x1 == x2 or not (min(x1, x2) <= xm <= max(x1, x2)):
                continue
            t_a = (xa - x1) / (x2 - x1)
            t_b = (xb - x1) / (x2 - x1)
            crossings.append((y1 + t_a * (y2 - y1), y1 + t_b * (y2 - y1)))
        # Strips are split at every vertex x, so edges cannot cross inside a
        # strip; sorting by midpoint height orders them bottom-to-top even
        # when two edges touch at a strip boundary (a wedge vertex).
        crossings.sort(key=lambda ab: ab[0] + ab[1])
        new_open = []
        used = [False] * len(open_bands)
        for j in range(0, len(crossings) - 1, 2):
            bot_a, bot_b = crossings[j]
            top_a, top_b = crossings[j + 1]
            band = None
            for k, ob in enumerate(open_bands):
                if used[k]:
                    continue
                if min(ob['top'][-1], top_a) - max(ob['bot'][-1], bot_a) > 1e-9:
                    band = ob
                    used[k] = True
                    break
            if band is None:
                band = {'xs': [xa], 'top': [top_a], 'bot': [bot_a]}
                bands.append(band)
            elif band['top'][-1] != top_a or band['bot'][-1] != bot_a:
                # A vertical edge steps the boundary at xa: keep both samples
                # (duplicate x) so the tables carry the jump.
                band['xs'].append(xa)
                band['top'].append(top_a)
                band['bot'].append(bot_a)
            band['xs'].append(xb)
            band['top'].append(top_b)
            band['bot'].append(bot_b)
            new_open.append(band)
        open_bands = new_open
    return [(np.array(b['xs']), np.array(b['top']), np.array(b['bot']))
            for b in bands if len(b['xs']) >= 2]


def _build_polygon_edges(polygons):
    """Precompute per-band top/bottom boundary tables and x-range once per
    generate_slices call, so slice weights can be computed with fast numpy
    interpolation instead of per-slice geometric intersection.

    Emits one entry per y-simple band (see _polygon_bands), so a polygon that
    folds back in x contributes several entries. 'poly_index' ties a band to
    its source polygon (stable layer-height columns, material fallback);
    'is_primary' marks one band per polygon for the exterior-based loops
    (slice breakpoints, failure-surface crossings) that must visit each
    polygon exactly once."""
    edges = []
    for p_idx, poly in enumerate(polygons):
        geom = poly['polygon']
        for b_idx, (txs, tys, bys) in enumerate(_polygon_bands(geom)):
            edges.append({
                'mat_id': poly['mat_id'],
                'poly_index': p_idx,
                'is_primary': b_idx == 0,
                'txs': txs, 'tys': tys, 'bxs': txs, 'bys': bys,
                'xmin': float(txs[0]), 'xmax': float(txs[-1]),
                'exterior': geom.exterior,
            })
    return edges


def _domain_bottom_is_flat(domain, n=24, tol=1e-6):
    """True if the lower boundary of the domain polygon is horizontal.

    Sampled by vertical lines: the domain is flat-bottomed when its minimum y at
    every x equals the global minimum y. Flat bottoms (profile-line inputs and
    polygon inputs with a horizontal base) can use the cheap scalar depth clamp
    instead of a geometric containment check (plan_polygons.md §5.4).
    """
    minx, miny, maxx, maxy = domain.bounds
    span = maxx - minx
    if span <= 0:
        return True
    for x in np.linspace(minx + 0.01 * span, maxx - 0.01 * span, n):
        inter = domain.intersection(LineString([(x, miny - 1.0), (x, maxy + 1.0)]))
        if inter.is_empty or abs(inter.bounds[1] - miny) > tol:
            return False
    return True


def _domain_containment(slope_data):
    """Return a prepared domain polygon for fast repeated containment checks, or
    None when no check is needed (flat bottom, or no domain). Cached on slope_data
    so the prep() and flatness test run once per search, not per trial surface."""
    if '_domain_check' not in slope_data:
        domain = slope_data.get('domain_polygon')
        if domain is None or _domain_bottom_is_flat(domain):
            slope_data['_domain_check'] = None
        else:
            from shapely.prepared import prep
            # small tolerance: clipped surfaces legitimately touch the domain
            # boundary (ground-vertex entries, bedrock-grazing arcs) and exact
            # covers() fails on float noise there
            slope_data['_domain_check'] = prep(domain.buffer(1e-6))
    return slope_data['_domain_check']


def domain_lower_envelope(domain):
    """Return the lower boundary of the domain polygon as (x, y) points, left to
    right, sampled at the exterior vertex x-values (piecewise-linear between)."""
    ring = list(domain.exterior.coords)
    if len(ring) > 1 and ring[0] == ring[-1]:
        ring = ring[:-1]
    n = len(ring)
    pts = []
    for x in sorted(set(p[0] for p in ring)):
        ys = []
        for i in range(n):
            x1, y1 = ring[i]
            x2, y2 = ring[(i + 1) % n]
            if (x1 <= x <= x2) or (x2 <= x <= x1):
                if x1 == x2:
                    ys.extend([y1, y2])
                else:
                    ys.append(y1 + (x - x1) / (x2 - x1) * (y2 - y1))
        if ys:
            pts.append((x, min(ys)))
    return pts


def _elastic_cores(slope_data):
    """Cached (material name, shrunken-core Polygon) for every ``option='elastic'``
    zone. An elastic material is impenetrable (v16): a failure surface may RIDE
    ALONG its boundary (a slide daylighting on the top of bedrock) but may not cut
    into its interior. The core is the zone buffered inward by a small tolerance, so
    a surface on the boundary misses it while one dipping inside hits it.

    Empty list for the overwhelmingly common no-elastic model, so the crossing
    check that uses it is a zero-cost early-out there.
    """
    if '_elastic_cores' not in slope_data:
        mats = slope_data.get('materials') or []
        domain = slope_data.get('domain_polygon')
        if domain is not None:
            minx, miny, maxx, maxy = domain.bounds
            tol = max(1e-6, 1e-4 * max(maxx - minx, maxy - miny))
        else:
            tol = 1e-6
        cores = []
        for p in slope_data.get('polygons') or []:
            mid = p.get('mat_id')
            if mid is None or not (0 <= mid < len(mats)):
                continue
            if str(mats[mid].get('option', '')).strip().lower() != 'elastic':
                continue
            core = p['polygon'].buffer(-tol)
            if not core.is_empty:
                cores.append((mats[mid].get('name', f'Material {mid + 1}'), core))
        slope_data['_elastic_cores'] = cores
    return slope_data['_elastic_cores']


def _domain_floor(slope_data):
    """Cached (xs, ys) arrays of the domain's lower boundary, for np.interp.

    None when the problem has no domain polygon. The floor is the impenetrable
    boundary of the model — bedrock for a polygon input, the max_depth line for a
    profile-line input.
    """
    if '_domain_floor' not in slope_data:
        domain = slope_data.get('domain_polygon')
        pts = domain_lower_envelope(domain) if domain is not None else []
        slope_data['_domain_floor'] = (
            (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
            if len(pts) >= 2 else None)
    return slope_data['_domain_floor']


class CompositeSurface:
    """A circular arc truncated at the bottom of the domain.

    A trial circle that dips below the impenetrable boundary of the model (bedrock,
    or the user's max_depth line) is not an admissible slip surface. Rather than
    reject it, every LEM code runs the surface ALONG the boundary between the two
    points where the arc crosses it. That is a COMPOSITE surface, and because the
    floor is single-valued in x it is simply the upper envelope of the two:

        y(x) = max(arc_y(x), floor_y(x))

    The arc portion keeps its exact circular geometry (y and alpha come from the
    circle equation, not from a polyline chord), and the floor portion takes the
    slope of the floor segment it lies on. The two crossing points are returned so
    the caller can force slice boundaries there, which keeps every slice base
    wholly on one branch or the other — no slice straddles the kink.

    The circle center is still carried into slice_df, because the moment methods
    (OMS, Bishop) take moments about it. Their moment arms are no longer the
    constant R, though; see solve.oms/solve.bishop.
    """

    def __init__(self, circle, fx, fy, crossings, x_min, x_max):
        self.circle = circle
        self.fx, self.fy = fx, fy
        self.crossings = crossings
        self.x_min, self.x_max = x_min, x_max

    def floor_y(self, x):
        return np.interp(x, self.fx, self.fy)

    def y(self, x):
        """Elevation of the composite surface at x (vectorized)."""
        x = np.asarray(x, dtype=float)
        c = self.circle
        arc = get_circular_y_coordinates(x, c['Xo'], c['Yo'], c['R'])
        # nan outside the circle's x-span: those x lie on the floor by definition
        arc = np.where(np.isnan(arc), -np.inf, arc)
        return np.maximum(arc, self.floor_y(x))

    def on_floor(self, x, tol=1e-9):
        c = self.circle
        arc = get_circular_y_coordinates(np.atleast_1d(float(x)), c['Xo'], c['Yo'], c['R'])[0]
        return np.isnan(arc) or arc <= self.floor_y(x) + tol

    def alpha_deg(self, x):
        """Base inclination at x, degrees, positive when the surface rises to the
        right (the same left-to-right convention as the non-circular path; the
        caller applies the right-facing sign flip)."""
        if self.on_floor(x):
            i = int(np.clip(np.searchsorted(self.fx, x) - 1, 0, len(self.fx) - 2))
            dx = self.fx[i + 1] - self.fx[i]
            return degrees(atan2(self.fy[i + 1] - self.fy[i], dx)) if dx > 0 else 0.0
        c = self.circle
        dx_circle = x - c['Xo']
        R = c['R']
        if abs(dx_circle) >= R:
            return 0.0
        return degrees(atan(dx_circle / sqrt(R ** 2 - dx_circle ** 2)))

    def line_string(self, y_left, y_right, n=400):
        """Dense polyline of the composite surface, for plotting and for the
        material-boundary / piezo-line intersections that build slice breakpoints."""
        xs = set(np.linspace(self.x_min, self.x_max, n))
        xs.update(self.crossings)
        xs.update(x for x in self.fx if self.x_min < x < self.x_max)
        xs = np.array(sorted(xs))
        ys = self.y(xs)
        ys[0], ys[-1] = y_left, y_right      # pin to the exact ground intersections
        return LineString(list(zip(xs, ys)))


def build_composite_surface(slope_data, circle, x_min, x_max, n=2000):
    """Build the CompositeSurface for a circle that dips below the domain floor.

    Returns None when the arc clears the floor over its whole span — the ordinary
    circular case, which must stay on its exact-arc fast path.
    """
    floor = _domain_floor(slope_data)
    if floor is None:
        return None
    fx, fy = floor
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']

    xs = np.linspace(x_min, x_max, n)
    arc = get_circular_y_coordinates(xs, Xo, Yo, R)
    d = arc - np.interp(xs, fx, fy)          # < 0 where the arc is below the floor
    below = np.isnan(d) | (d < -1e-9)
    if not below.any():
        return None

    # Exact crossings by bisection on each sign change of d(x). An irregular floor
    # can be crossed more than twice; every crossing becomes a slice breakpoint.
    def d_at(x):
        a = get_circular_y_coordinates(np.array([x]), Xo, Yo, R)[0]
        return -1e9 if np.isnan(a) else a - float(np.interp(x, fx, fy))

    crossings = []
    for i in range(len(xs) - 1):
        if below[i] != below[i + 1]:
            lo, hi = xs[i], xs[i + 1]
            lo_below = below[i]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if (d_at(mid) < 0) == lo_below:
                    lo = mid
                else:
                    hi = mid
            crossings.append(0.5 * (lo + hi))

    return CompositeSurface(circle, fx, fy, crossings, x_min, x_max)


def generate_slices(slope_data, circle=None, non_circ=None, num_slices=40, debug=True,
                    composite=False, right_facing=None,
                    suction_phi_b=None, suction_cap=None, check_inputs=True,
                    material_breakdown=False):

    """
    Generates vertical slices between the ground surface and a failure surface for slope stability analysis.

    This function supports both circular and non-circular failure surfaces and computes
    geometric and mechanical properties for each slice, including weight, base geometry,
    water pressures, distributed loads, and reinforcement effects.

    Parameters:
        data (dict): Dictionary containing all input data
        circle (dict, optional): Dictionary with keys 'Xo', 'Yo', 'Depth', and 'R' defining the circular failure surface.
        non_circ (list, optional): List of dicts defining a non-circular failure surface with keys 'X', 'Y', and 'Movement'.
        num_slices (int, optional): Desired number of slices to generate (default is 40).
        debug (bool, optional): Whether to print debug information (default is True).
        right_facing (bool or None, optional): Facing override. None (default)
            auto-detects facing (``y_left > y_right`` for a normal slope, and the
            surface-asymmetry rule for a flat arc). A bool forces the facing and
            wins over any auto-detection; it also accepts an in-memory
            ``slope_data['right_facing']`` key (the kwarg takes precedence).
        suction_phi_b (dict or None, optional): Opt-in matric-suction strength
            (Fredlund extended Mohr-Coulomb), LEM only. Maps ``{material name:
            phi_b degrees}``. Default None auto-wires from the materials' template
            (v17) ``phi_b`` values; an explicit dict overrides the file, and an empty
            dict forces suction off. For a slice whose base material is named in the dict,
            the base matric suction ``s = max(0, -u)`` (from the UNCLAMPED pore
            pressure: a piezometric line's hydrostatic negative head above the
            line, or an unsaturated seepage solution's negative u) contributes an
            apparent cohesion ``c_suction = s * tan(phi_b)``. This realises
            ``tau = c' + (sigma - u_a) tan(phi') + (u_a - u_w) tan(phi_b)`` with
            ``u_a = 0``: the effective-normal term keeps u clamped at 0 (water-only,
            numerically safe) while suction is carried entirely as apparent cohesion.
            Default None => c_suction = 0.0 for every slice, bit-identical to the
            clamped baseline.
        suction_cap (float, dict, or None, optional): Optional upper bound (stress
            units) on the suction ``s`` before it is converted to apparent cohesion,
            so a deep piezometric surface cannot grow unbounded hydrostatic suction.
            A scalar caps every material identically; a dict ``{material name: cap}``
            caps per material (a material absent from the dict is uncapped). Default
            None auto-wires from the materials' template (v17) ``s_cap`` values (an
            explicit value overrides the file). Ignored when ``suction_phi_b`` is None.
        check_inputs (bool, optional): Run :func:`xslope.preflight.preflight` on
            ``slope_data`` first and refuse with a :class:`~xslope.preflight.PreflightError`
            when it finds an error -- an input that would crash the run or make its
            answer provably wrong. Default True. Pass False to bypass the checks:
            the searches do this after gating once at their own entry (the model
            does not change between trial surfaces), and the sweep engines do it
            because a deliberately perturbed value is not a user mistake.

    Returns:
        tuple:
            - pd.DataFrame: Slice table where each row includes geometry, strength, and external force values.
            - shapely.geometry.LineString: The clipped failure surface between the ground surface intersections.

    Notes:
        - Supports Method A interpretation of reinforcement: T reduces driving forces.
        - Handles pore pressure and distributed loads using linear interpolation at slice centers.
        - Automatically includes all geometry breakpoints in slice generation.
        - Must specify exactly one of 'circle' or 'non_circ'.
    """

    # === PREFLIGHT ===
    # The run gate. Everything the analysis needs and this model does not carry is
    # reported here, in the template's own vocabulary, before any geometry is built
    # -- rather than as a shapely exception naming no field, or as a number that is
    # quietly wrong. Errors refuse; warnings and infos are carried on the report and
    # do not block.
    if check_inputs:
        from .preflight import preflight as _preflight
        # `surface_supplied` says the surface came in as an ARGUMENT, so the
        # sheet-emptiness rule must not fire: a caller that builds its own circle
        # (a search's trial surface, a sweep step, a plot) is not missing one.
        _preflight(slope_data, 'lem',
                   {'surface': 'noncircular' if circle is None and non_circ is not None
                               else 'circular',
                    'surface_supplied': circle is not None or non_circ is not None},
                   ).raise_for_errors()

    # Validate material properties
    materials = slope_data["materials"]
    if all(m.get('c', 0) == 0 and m.get('phi', 0) == 0 and m.get('gamma', 0) == 0 for m in materials):
        return False, "All materials have empty strength properties (c, phi, gamma). Check your input template."

    # Unpack data
    # Geometry is represented internally as material-zone polygons. Slice weights,
    # base material, layer heights, and slice-boundary breakpoints are all computed
    # from these polygons.
    polygons = slope_data["polygons"]
    poly_edges = _build_polygon_edges(polygons)
    n_polygons = len(polygons)
    ground_surface = slope_data["ground_surface"]
    piezo_line = slope_data["piezo_line"]
    piezo_line2 = slope_data.get("piezo_line2", [])  # Second piezometric line
    gamma_w = slope_data["gamma_water"]
    tcrack_depth = slope_data["tcrack_depth"]
    tcrack_water = slope_data["tcrack_water"]
    k_seismic = slope_data['k_seismic']
    # === WATER LOADS (v22 main!D23) ===
    # In automatic mode the weight of standing water is not something the user
    # typed: it is derived here from the model's own water definition -- the
    # piezometric line, or the seepage head boundaries where a seepage analysis is
    # defined -- for the stage and moment being solved. The derived blocks arrive
    # in keys of their own and are APPENDED to the user's, never merged into them,
    # so the slice forces carry both while the plots and the double-count rule can
    # still tell which is which. Water always acts normal to the surface, so a
    # derived block takes the default direction: the dirs lists stay as long as the
    # user's own lists, and every index past their end reads 'normal'.
    #
    # In manual mode with_water_loads returns the caller's model by identity and
    # both derived lists are empty, so every list below is the one this function
    # has always built and the analysis is bit-identical.
    # Resolve any overburden-dependent pullout law BEFORE with_water_loads, which
    # in auto mode hands back a shallow copy: the refreshed point list has to land
    # on the caller's model, not on a copy this function then drops.
    if slope_data.get("reinforcement_lines"):
        from .fileio import ensure_reinforce_pullout
        ensure_reinforce_pullout(slope_data)
    from .water import with_water_loads, DERIVED_KEYS
    slope_data = with_water_loads(slope_data)
    dloads = list(slope_data["dloads"]) + list(slope_data.get(DERIVED_KEYS[1]) or [])
    dloads2 = list(slope_data.get("dloads2") or []) + list(
        slope_data.get(DERIVED_KEYS[2]) or [])
    dload_dirs = slope_data.get("dload_dirs") or []
    dload2_dirs = slope_data.get("dload2_dirs") or []

    # Auto-wire the per-material matric-suction parameters from the template (v17)
    # when the caller passes no explicit kwarg -- t_cut override semantics: an
    # explicit kwarg WINS over the file. load_slope_data carries phi_b/s_cap on each
    # material dict; a material with phi_b set contributes suction strength, capped
    # at its own s_cap. When no material carries phi_b the derived dict is empty ->
    # None, bit-identical to the pre-v17 default-off path. (None means "unset, read
    # the file"; pass an empty dict to force suction off regardless of the file.)
    if suction_phi_b is None:
        _file_phi_b = {m.get('name'): m.get('phi_b') for m in materials
                       if m.get('phi_b') is not None}
        suction_phi_b = _file_phi_b or None
    if suction_cap is None:
        _file_cap = {m.get('name'): m.get('s_cap') for m in materials
                     if m.get('s_cap') is not None}
        suction_cap = _file_cap or None

    # Opt-in matric-suction strength (Fredlund extended Mohr-Coulomb). Warn on a
    # phi_b keyed to a material name that does not exist, so a typo silently
    # produces zero suction rather than an error.
    if suction_phi_b:
        _mat_names = {m.get('name') for m in materials}
        for _nm in suction_phi_b:
            if _nm not in _mat_names:
                warnings.warn(
                    f"suction_phi_b names material '{_nm}', which is not in the model "
                    f"(materials: {sorted(n for n in _mat_names if n)}). No suction "
                    "strength will be applied for it.")

    # Warn once if seep pore pressure is selected but seep data is missing. This is
    # only ever a heads-up now: any slice whose BASE material takes u='seep' raises
    # from the per-slice branch below (_no_seep_field_error) rather than reading
    # zero, so the once-per-process latch can no longer hide a real dropout — it
    # only suppresses repeat notices during a search, which calls this thousands of
    # times. The case that reaches here and does NOT raise is the benign one: a
    # seep material exists in the model but no slice base lands in it, so there is
    # no pore pressure to lose.
    has_seep_materials = any(m["u"] == "seep" for m in materials)
    has_seep_data = (slope_data.get('mesh') is not None
                     and slope_data.get('seep_u') is not None)
    if has_seep_materials and not has_seep_data and not getattr(generate_slices, '_seep_warned', False):
        print("WARNING: Seep pore pressure option selected but required seep files were not "
              "found. Any slice base in a seep material will be refused rather than read as dry.")
        generate_slices._seep_warned = True

    # Determine failure surface type
    if circle is not None:
        circular = True
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
    else:
        circular = False

    # Prepare reinforcement lines data. Preferred source is the raw
    # 'reinforcement_lines' dicts, which carry the full capacity envelope
    # (t_max/lp/tend) plus the v12 Dir/Appl settings; the available tension at a
    # crossing is evaluated exactly via fileio.reinforce_available_tension.
    #
    # An EXPLICIT 'reinforcement_lines' list is the model's own statement about
    # its reinforcement, the empty list included: a model that says it has no
    # reinforcement lines has none. 'reinforce_lines' is DERIVED from that list
    # (fileio.build_reinforce_lines), so falling back to it whenever the source
    # was empty resurrected the reinforcement a user had just deleted -- "remove
    # the reinforcement" silently did nothing, because the point lists built
    # before the edit were still there and still crossed the surface. The legacy
    # point-list path (tangent/active, interpolated on X) is therefore reached
    # only when the model carries no 'reinforcement_lines' key at all --
    # slope_data assembled by hand.
    reinf_lines_data = []
    _reinf_source = slope_data.get("reinforcement_lines")
    if _reinf_source:
        from .fileio import reinforce_available_tension  # shared envelope
        for _ri, r in enumerate(_reinf_source):
            dxl = r["x2"] - r["x1"]
            dyl = r["y2"] - r["y1"]
            length = np.hypot(dxl, dyl)
            if length == 0:
                continue
            # Line inclination in the same angular convention as the slice base
            # angle alpha: measured from +x with the direction normalized so
            # cos(psi) > 0 (left-to-right), and mirrored for right-facing slopes
            # exactly as alpha/beta are.
            if dxl < 0:
                dxl, dyl = -dxl, -dyl
            psi_line = atan2(dyl, dxl)   # radians
            reinf_lines_data.append({
                "geom": LineString([(r["x1"], r["y1"]), (r["x2"], r["y2"])]),
                "envelope": (r["t_max"], r["lp1"], r["lp2"],
                             r.get("tend1") or 0.0, r.get("tend2") or 0.0),
                "pullout": r.get("_pullout_profile"),
                "length": length,
                "dir": _support_choice(r, "dir", ("tangent", "axial"),
                                       "tangent", _ri, "Reinforcement line"),
                "appl": _support_choice(r, "appl", ("active", "passive"),
                                        "active", _ri, "Reinforcement line"),
                "psi": psi_line,
                "avail": reinforce_available_tension,
                "claimed": set(),
            })
    elif _reinf_source is None and slope_data.get("reinforce_lines"):
        for line in slope_data["reinforce_lines"]:
            xs = [pt["X"] for pt in line]
            ts = [pt["T"] for pt in line]
            geom = LineString([(pt["X"], pt["Y"]) for pt in line])
            reinf_lines_data.append({"xs": xs, "ts": ts, "geom": geom,
                                     "dir": "tangent", "appl": "active",
                                     "claimed": set()})

    # Prepare pile lines data
    pile_lines_data = []
    if slope_data.get("pile_lines"):
        for _pi, pile in enumerate(slope_data["pile_lines"]):
            geom = LineString([(pile["x1"], pile["y1"]), (pile["x2"], pile["y2"])])
            pile_lines_data.append({
                "geom": geom,
                "H": pile["H"] if pile["H"] is not None else None,
                "theta_p": pile["theta_p"],
                "D_pile": pile.get("D_pile"),
                "S": pile.get("S"),
                "V_cap": pile.get("V_cap"),
                "M_cap": pile.get("M_cap"),
                # Same guard as the reinforcement Dir/Appl above: a value that is
                # not an application was applied as PASSIVE capacity.
                "appl": _support_choice(pile, "appl", ("active", "passive"),
                                        "active", _pi, "Pile"),
                "label": pile.get("label", ""),
                "claimed": set(),
            })
    # One record per pile-row crossing of the failure surface, carried out on
    # df.attrs['pile_report'] so the run summary can print the Ito & Matsui
    # chain (soil force -> governing capacity -> applied force) without a
    # report round-trip.
    pile_report = []

    # Line loads (v12 'lloads' sheet): concentrated forces per unit width on the
    # ground surface, assigned below to the slice whose top contains each point.
    line_loads_data = list(slope_data.get("line_loads") or [])

    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])

    # Resolve the facing override once: the explicit kwarg wins, else an
    # in-memory slope_data['right_facing'] key, else None (auto-detect). Only
    # ever non-None for the rare flat-arc / caller-forced case; for a normal
    # slope this stays None and facing is detected exactly as before.
    facing_override = right_facing if right_facing is not None else slope_data.get('right_facing')

    # Generate failure surface
    success, result = generate_failure_surface(ground_surface, circular, circle=circle, non_circ=non_circ,
                                               tcrack_depth=tcrack_depth, right_facing=facing_override)
    if success:
        x_min, x_max, y_left, y_right, clipped_surface = result
    else:
        return False, "Failed to generate surface:" + result

    # === Composite surfaces (opt-in) ===
    # With composite=True, a circle that dips below the bottom of the model is
    # truncated at it and follows it between the crossings (see CompositeSurface) —
    # the standard way an LEM code handles a slip surface that would otherwise cut
    # through bedrock. It is OFF by default because the floor of a profile-line
    # model is `max_depth`, a search bound the user chose rather than a real
    # impenetrable boundary; truncating there would be meaningless. With
    # composite=False a below-floor circle is rejected outright, as before.
    #
    # The exact-arc fast paths below stay in force whenever the arc clears the
    # floor, which is the overwhelmingly common case; `use_arc` gates them.
    comp = None
    if circular and composite:
        comp = build_composite_surface(slope_data, circle, x_min, x_max)
        if comp is not None:
            clipped_surface = comp.line_string(y_left, y_right)
    use_arc = circular and comp is None

    def surf_y(xs):
        """Failure-surface elevation at each x. Exact for the arc and composite
        cases; a vertical-line intersection for a hand-entered polyline."""
        if use_arc:
            return get_circular_y_coordinates(xs, Xo, Yo, R)
        if comp is not None:
            return comp.y(xs)
        return np.array([get_y_from_intersection(
            clipped_surface.intersection(LineString([(x, -1e6), (x, 1e6)]))) for x in xs])

    # Reject failure surfaces that still leave the domain polygon (plan_polygons.md
    # §5.4) — a hand-entered non-circular surface below the floor, or a composite
    # surface that failed to close. Only irregular bottoms need the geometric
    # covers() test; for a flat bottom the same check is a scalar depth comparison.
    prepared_domain = _domain_containment(slope_data)
    if prepared_domain is not None:
        if not prepared_domain.covers(clipped_surface):
            return False, "Failure surface extends outside the domain polygon"
    else:
        domain = slope_data.get('domain_polygon')
        if domain is not None:
            y_bot = domain.bounds[1]
            surf_min_y = min(y for _, y in clipped_surface.coords)
            if surf_min_y < y_bot - 1e-6:
                return False, (
                    f"Failure surface reaches y={surf_min_y:.3f}, below the bottom of "
                    f"the domain (y={y_bot:.3f}). Raise the surface, lower max_depth, or "
                    f"pass composite=True to truncate the circle at the bottom of the "
                    f"model and run it along the base.")

    # An elastic material is impenetrable (v16): the surface may ride along its
    # boundary but may not cut into it. Reject a surface that penetrates one, the
    # same way a below-floor surface is rejected above; inside circular_search this
    # scores the trial circle fs_fail and drops it (the search-side rejection).
    for _elastic_name, _core in _elastic_cores(slope_data):
        if _core.intersects(clipped_surface):
            return False, (
                f"Failure surface crosses elastic (impenetrable) material "
                f"'{_elastic_name}'. An elastic zone cannot fail; keep the surface "
                f"above it (it may run along its boundary).")

    # Determine if the failure surface is right-facing. Byte-identical to the
    # historical `y_left > y_right` for a normal slope; the surface-asymmetry
    # rule (and the facing override) engage ONLY on a flat arc — the equal-
    # elevation case that used to be rejected in get_sorted_intersections.
    right_facing, facing_note = _resolve_right_facing(
        y_left, y_right, list(clipped_surface.coords), x_min, x_max,
        override=facing_override)

    # === BEGIN : Find set of points that should correspond to slice boundaries. ===

    # Find set of points that are on the polygon boundaries if the points are above
    # the failure surface (these become slice-boundary breakpoints).
    fixed_xs = set()

    # Vectorized approach for polygon boundary vertices
    for pe in poly_edges:
        if not pe['is_primary']:
            continue
        poly_coords = np.array(pe['exterior'].coords)
        x_coords = poly_coords[:, 0]
        y_coords = poly_coords[:, 1]

        # Filter points within x-range
        mask = (x_coords >= x_min) & (x_coords <= x_max)
        x_filtered = x_coords[mask]
        y_filtered = y_coords[mask]

        if len(x_filtered) == 0:
            continue

        # Check if points are above the failure surface
        if circular:
            failure_y = surf_y(x_filtered)
            above_mask = y_filtered > failure_y
        else:
            # For non-circular, use geometric intersection (slower but necessary)
            above_mask = np.zeros(len(x_filtered), dtype=bool)
            for i, (x, y) in enumerate(zip(x_filtered, y_filtered)):
                vertical_line = LineString([(x, -1e6), (x, 1e6)])
                failure_y = get_y_from_intersection(clipped_surface.intersection(vertical_line))
                if failure_y is not None and y > failure_y:
                    above_mask[i] = True

        # Add points that are above the failure surface
        fixed_xs.update(x_filtered[above_mask])

    fixed_xs.update([x_min, x_max])

    # The arc/floor crossings are kinks in the base: force a slice boundary at each
    # so that no slice base straddles one.
    if comp is not None:
        fixed_xs.update(x for x in comp.crossings if x_min < x < x_max)

    # Add transition points from dloads
    if dloads:
        fixed_xs.update(
            pt['X'] for line in dloads for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from dloads2
    if dloads2:
        fixed_xs.update(
            pt['X'] for line in dloads2 for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from non_circ
    if non_circ:
        fixed_xs.update(
            pt['X'] for pt in non_circ
            if x_min <= pt['X'] <= x_max
        )

    # Find intersections between material-zone boundaries and the failure surface
    if use_arc:
        # For circular failure surfaces, we can use a more efficient approach
        # by creating a dense circle representation and finding intersections
        theta_range = np.linspace(np.pi, 2 * np.pi, 200)
        circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        circle_line = LineString(circle_coords)
        for pe in poly_edges:
            if not pe['is_primary']:
                continue
            intersection = pe['exterior'].intersection(circle_line)
            if not intersection.is_empty:
                if hasattr(intersection, 'x'):
                    if x_min <= intersection.x <= x_max:
                        fixed_xs.add(intersection.x)
                elif hasattr(intersection, 'geoms'):
                    for geom in intersection.geoms:
                        if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                            fixed_xs.add(geom.x)
    else:
        # For non-circular failure surfaces, use the original approach
        for pe in poly_edges:
            if not pe['is_primary']:
                continue
            intersection = pe['exterior'].intersection(clipped_surface)
            if not intersection.is_empty:
                if hasattr(intersection, 'x'):
                    if x_min <= intersection.x <= x_max:
                        fixed_xs.add(intersection.x)
                elif hasattr(intersection, 'geoms'):
                    for geom in intersection.geoms:
                        if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                            fixed_xs.add(geom.x)

    # Find intersections with piezometric lines. A single point is not a line: it
    # cannot become a LineString and it cannot be interpolated, so it is treated as
    # no line here and reported by the u='piezo' guard below.
    if piezo_line and len(piezo_line) >= 2:
        piezo_geom1 = LineString(piezo_line)
        if use_arc:
            # Use dense circle representation for intersection
            theta_range = np.linspace(np.pi, 2 * np.pi, 200)
            circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
            circle_line = LineString(circle_coords)
            intersection1 = piezo_geom1.intersection(circle_line)
        else:
            intersection1 = piezo_geom1.intersection(clipped_surface)
            
        if not intersection1.is_empty:
            if hasattr(intersection1, 'x'):
                # Single point intersection
                if x_min <= intersection1.x <= x_max:
                    fixed_xs.add(intersection1.x)
            elif hasattr(intersection1, 'geoms'):
                # Multiple points or line intersection
                for geom in intersection1.geoms:
                    if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                        fixed_xs.add(geom.x)

    if piezo_line2 and len(piezo_line2) >= 2:
        piezo_geom2 = LineString(piezo_line2)
        if use_arc:
            theta_range = np.linspace(np.pi, 2 * np.pi, 200)
            circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
            circle_line = LineString(circle_coords)
            intersection2 = piezo_geom2.intersection(circle_line)
        else:
            intersection2 = piezo_geom2.intersection(clipped_surface)
            
        if not intersection2.is_empty:
            if hasattr(intersection2, 'x'):
                # Single point intersection
                if x_min <= intersection2.x <= x_max:
                    fixed_xs.add(intersection2.x)
            elif hasattr(intersection2, 'geoms'):
                # Multiple points or line intersection
                for geom in intersection2.geoms:
                    if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                        fixed_xs.add(geom.x)

    # Remove duplicate points that are very close to each other
    tolerance = 1e-6
    cleaned_xs = []
    for x in sorted(fixed_xs):
        if not cleaned_xs or abs(x - cleaned_xs[-1]) > tolerance:
            cleaned_xs.append(x)
    
    fixed_xs = cleaned_xs

    # Generate slice boundaries
    segment_lengths = [fixed_xs[i + 1] - fixed_xs[i] for i in range(len(fixed_xs) - 1)]
    total_length = sum(segment_lengths)
    all_xs = [fixed_xs[0]]
    for i in range(len(fixed_xs) - 1):
        x_start = fixed_xs[i]
        x_end = fixed_xs[i + 1]
        segment_length = x_end - x_start
        n_subdiv = max(1, int(round((segment_length / total_length) * num_slices)))
        xs = np.linspace(x_start, x_end, n_subdiv + 1).tolist()
        all_xs.extend(xs[1:])

    # === END : Find set of points that should correspond to slice boundaries. ===

    # Remove thin slices (width < 1e-2), including at the ends
    min_width = 1e-2
    cleaned_xs = [all_xs[0]]
    for x in all_xs[1:]:
        if abs(x - cleaned_xs[-1]) >= min_width:
            cleaned_xs.append(x)
    
    # If the last slice is thin, merge it with the previous one
    if len(cleaned_xs) > 2 and abs(cleaned_xs[-1] - cleaned_xs[-2]) < min_width:
        cleaned_xs.pop(-2)
    
    all_xs = cleaned_xs

    # Pre-compute all y-coordinates for efficiency
    slice_x_coords = np.array(all_xs)
    slice_centers = (slice_x_coords[:-1] + slice_x_coords[1:]) / 2
    
    # Get failure surface y-coordinates
    y_lb_all = surf_y(slice_x_coords[:-1])
    y_rb_all = surf_y(slice_x_coords[1:])
    y_cb_all = surf_y(slice_centers)

    # Get ground surface y-coordinates
    y_lt_all = get_ground_surface_y_coordinates(slice_x_coords[:-1], ground_surface)
    y_rt_all = get_ground_surface_y_coordinates(slice_x_coords[1:], ground_surface)
    y_ct_all = get_ground_surface_y_coordinates(slice_centers, ground_surface)

    # Precompute each polygon's top/bottom boundary at every slice center, plus an
    # in-range mask. This is the polygon-based replacement for per-layer profile
    # heights — fast numpy interpolation, no per-slice geometric intersection.
    centers_arr = np.asarray(slice_centers)
    poly_top_all = []
    poly_bot_all = []
    poly_inrange = []
    for pe in poly_edges:
        poly_top_all.append(np.interp(centers_arr, pe['txs'], pe['tys']))
        poly_bot_all.append(np.interp(centers_arr, pe['bxs'], pe['bys']))
        poly_inrange.append((centers_arr >= pe['xmin'] - 1e-9) &
                            (centers_arr <= pe['xmax'] + 1e-9))

    # Get piezometric y-coordinates. A slice center sitting ON a line's own end is
    # INSIDE it, so the query is first snapped across floating-point noise: the
    # corpus is full of piezometric lines drawn exactly to the domain edge, and a
    # last-ulp disagreement there must not read as "outside the line" below.
    def _snap_to_extent(x, line):
        if not line:
            return x
        xs = [p[0] for p in line]
        lo, hi = min(xs), max(xs)
        tol = 1e-9 * max(1.0, hi - lo)
        x = np.asarray(x, dtype=float).copy()
        x[(x < lo) & (x >= lo - tol)] = lo
        x[(x > hi) & (x <= hi + tol)] = hi
        return x

    piezo_y_all = get_piezometric_y_coordinates(
        _snap_to_extent(slice_centers, piezo_line), piezo_line)
    piezo_y2_all = get_piezometric_y_coordinates(
        _snap_to_extent(slice_centers, piezo_line2), piezo_line2)

    # === Piezometric-line extent guard ===
    # get_piezometric_y_coordinates interpolates with left/right = NaN, so a slice
    # center beyond the line's own x-range reads NaN. A material with u='piezo' at
    # such a slice used to fall through to u = 0 -- no pore pressure, silently, and
    # a non-conservative factor of safety. A piezometric line that legitimately
    # stops short (an upstream pool with no downstream pond, say) is fine as long as
    # the failure surface stays inside it; the moment a slice that ACTUALLY samples
    # it falls outside, that is a modelling error and is raised below, per slice.
    # No line AT ALL under u='piezo' is the same silent-zero class -- there is no
    # reading of "take the pore pressure from the piezometric line, but there isn't
    # one" that means a dry slope (a dry slope is u='none') -- so it is refused too,
    # below, at the same layer. A missing SECOND line is different: only rapid
    # drawdown reads it, and its absence is reported by the drawdown driver.
    _pz_have = bool(piezo_line) and len(piezo_line) >= 2
    _pz_ext = ((min(p[0] for p in piezo_line), max(p[0] for p in piezo_line))
               if _pz_have else None)
    _pz2_ext = ((min(p[0] for p in piezo_line2), max(p[0] for p in piezo_line2))
                if piezo_line2 and len(piezo_line2) >= 2 else None)
    if circular:
        _surf_desc = f"circle Xo={Xo:.5g}, Yo={Yo:.5g}, R={R:.5g}"
    else:
        _surf_desc = "the non-circular failure surface"

    def _piezo_extent_error(idx, xc, mat_name, ext, which):
        return ValueError(
            f"Pore pressure: slice {idx + 1} of {_surf_desc} has its base center at "
            f"x = {xc:.6g}, outside the {which}'s x-extent "
            f"[{ext[0]:.6g}, {ext[1]:.6g}]. Its base material '{mat_name}' takes pore "
            f"pressure from that line (u='piezo'), so there is no piezometric surface "
            f"to read and the pore pressure would silently be zero -- an unsafe, "
            f"non-conservative factor of safety. Extend the piezometric line to span "
            f"the failure surface, or give the material outside it a pore-pressure "
            f"option that applies there.")

    def _no_piezo_line_error(mat_name):
        return ValueError(
            f"Pore pressure: material '{mat_name}' takes pore pressure from a "
            f"piezometric line (u='piezo'), but the file defines no piezometric "
            f"line -- the piezo sheet carries fewer than two points. Every slice "
            f"base in that material would read zero pore pressure -- an unsafe, "
            f"non-conservative factor of safety. Draw the piezometric line, or set "
            f"u='none' for a dry model (or 'seep' / 'ru' for the other "
            f"pore-pressure sources).")

    def _no_seep_field_error(mat_name, missing):
        return ValueError(
            f"Pore pressure: material '{mat_name}' takes pore pressure from a "
            f"seepage solution (u='seep'), but this model carries no {missing} -- "
            f"the '{{base}}_mesh.json' / '{{base}}_seep.csv' companion files did not "
            f"load. Every slice base in that material would read zero pore pressure "
            f"-- an unsafe, non-conservative factor of safety. Re-run the seepage "
            f"analysis, or set u='none' for a dry model (or 'piezo' / 'ru' for the "
            f"other pore-pressure sources).")

    # === Water table for the gamma/gamma_sat weight split (template v12) ===
    # The water table is GLOBAL — one phreatic surface per problem, independent
    # of any material's pore-pressure option (the "sidecar" model). It governs
    # slice WEIGHT only; base pore pressure stays with the per-material u option.
    # Source precedence: a seepage solution's u = 0 contour (root-found on the
    # UNCLAMPED signed field) beats a hand-drawn piezo line; rapid drawdown
    # deliberately keys weight to the PRE-drawdown (stage-1) surfaces, never the
    # staged ones — the premise of rapid drawdown is that the soil stays
    # saturated while the pore pressures fall.
    any_gamma_sat = any(m.get('gamma_sat') is not None for m in materials)
    water_table_y_all = np.full(len(slice_centers), np.nan)
    if any_gamma_sat:
        if has_seep_data:
            had_cache = slope_data.get('_water_table_profile') is not None
            cache = seep_water_table_profile(slope_data)
            if not had_cache and piezo_line:
                print("gamma_sat weight split: using the seepage solution's u = 0 "
                      "contour as the water table; the piezometric line is NOT "
                      "used for unit weight (it may still supply pore pressure "
                      "and plotting).")
            xs_grid, wt = cache
            water_table_y_all = np.interp(slice_centers, xs_grid, wt)
        elif piezo_line:
            water_table_y_all = np.asarray(piezo_y_all, dtype=float)
        else:
            warnings.warn(
                "One or more materials specify gamma_sat, but the model has no "
                "water table (no piezometric line or seepage solution) — the "
                "saturated unit weight can never apply and gamma is used throughout.")

    # Interpolation functions for distributed loads. np.interp requires the
    # sample points to be in ascending-X order; sort each load line so a load
    # entered right-to-left is not silently interpolated to zero everywhere.
    #
    # Blocks are split by DIRECTION (v21): a 'normal' block acts perpendicular to the
    # loaded surface, a 'vertical' block resolves the same resultant straight down.
    # The two are summed as separate resultants below, so a model may mix them.
    def _interp_funcs(lines, dirs, want):
        out = []
        for i, line in enumerate(lines):
            if str(dirs[i] if i < len(dirs) else 'normal').lower() != want:
                continue
            pts = sorted(line, key=lambda pt: pt['X'])
            xs = [pt['X'] for pt in pts]
            normals = [pt['Normal'] for pt in pts]
            out.append(lambda x, xs=xs, normals=normals:
                       np.interp(x, xs, normals, left=0, right=0))
        return out

    dload_interp_funcs = _interp_funcs(dloads, dload_dirs, 'normal')
    dload_vert_funcs = _interp_funcs(dloads, dload_dirs, 'vertical')
    dload2_interp_funcs = _interp_funcs(dloads2, dload2_dirs, 'normal')
    dload2_vert_funcs = _interp_funcs(dloads2, dload2_dirs, 'vertical')

    # Generate slices
    slices = []
    # Per-slice, per-band geometry the material-derived columns are built from,
    # captured only on request (material_breakdown=True): the moist/saturated
    # heights and their first moments, per polygon band, in band order. They are
    # what update_slice_materials() needs to rebuild w, y_cg, kw and the base
    # strengths for a new materials list without re-slicing — the geometry never
    # changes, only the properties multiplied onto it.
    mb_Hm, mb_Hs, mb_HYm, mb_HYs, mb_base = ([], [], [], [], []) \
        if material_breakdown else (None, None, None, None, None)
    # Coverage tally for u='seep' base points. A single base outside the seepage
    # mesh is a warning (a mesh that stops short of one end of a long surface is a
    # modelling judgement); EVERY base outside it is not a judgement, it is a
    # broken read, and it is silent — the whole surface then prices as dry. See the
    # post-loop check.
    _seep_pts = 0
    _seep_outside = 0
    for i in range(len(all_xs) - 1):
        x_l, x_r = slice_x_coords[i], slice_x_coords[i + 1]
        x_c = slice_centers[i]
        dx = x_r - x_l

        # Get pre-computed y-coordinates
        y_lb = y_lb_all[i]
        y_rb = y_rb_all[i]
        y_cb = y_cb_all[i]
        y_lt = y_lt_all[i]
        y_rt = y_rt_all[i]
        y_ct = y_ct_all[i]

        # Skip if any coordinates are invalid
        if any(np.isnan([y_lb, y_rb, y_cb, y_lt, y_rt, y_ct])):
            continue

        # Calculate beta (slope angle of the top edge) in degrees
        beta = degrees(atan2(y_rt - y_lt, x_r - x_l))
        if right_facing:
            beta = -beta

        # Calculate dl for the top surface (for distributed loads)
        dl_top = sqrt((x_r - x_l)**2 + (y_rt - y_lt)**2)

        # Calculate layer heights and weight from the material-zone polygons.
        # Each y-simple band's vertical extent at the slice center gives a layer
        # band; overlapping it with [failure surface, ground surface] gives the
        # height. Heights are accumulated per source polygon so the h1..hN
        # columns stay one-per-zone even when a folded zone has several bands.
        heights = [0] * n_polygons
        soil_weight = 0
        base_material_idx = None
        if material_breakdown:
            _hm = [0.0] * len(poly_edges)
            _hs = [0.0] * len(poly_edges)
            _hym = [0.0] * len(poly_edges)
            _hys = [0.0] * len(poly_edges)
        base_overlap_bot = float('inf')  # elevation of the deepest present layer's base
        sum_gam_h_y = 0  # for calculating center of gravity of slice
        sum_gam_h = 0    # ditto

        for p_idx, pe in enumerate(poly_edges):
            # Material index (0-based); fall back to polygon order if unset/out of range
            mat_id = pe['mat_id']
            if mat_id is not None and 0 <= mat_id < len(materials):
                mat_index = mat_id
            else:
                mat_index = pe['poly_index']

            poly_top = poly_top_all[p_idx][i]
            poly_bot = poly_bot_all[p_idx][i]
            if (not poly_inrange[p_idx][i]) or np.isnan(poly_top) or np.isnan(poly_bot):
                continue

            # Layer band [poly_bot, poly_top] clipped to [failure surface, ground]
            overlap_top = min(y_ct, poly_top)
            overlap_bot = max(y_cb, poly_bot)
            h = max(0, overlap_top - overlap_bot)
            heights[pe['poly_index']] += h

            if h > 0:
                gamma = materials[mat_index]['gamma']
                g_sat = materials[mat_index].get('gamma_sat')
                y_w = water_table_y_all[i] if any_gamma_sat else np.nan
                if g_sat is not None and not np.isnan(y_w):
                    # Split the band at the water table: gamma_sat below, gamma
                    # (moist) above. Degenerate cases fall out of the clamp: a
                    # water table below the band gives an all-moist band, above
                    # it an all-saturated one (the right answer under ponded
                    # water, whose free water arrives as a distributed load).
                    y_split = min(max(y_w, overlap_bot), overlap_top)
                    h_sat = y_split - overlap_bot
                    h_moist = overlap_top - y_split
                    sum_gam_h_y += (g_sat * h_sat * (y_split + overlap_bot) / 2
                                    + gamma * h_moist * (overlap_top + y_split) / 2)
                    sum_gam_h += g_sat * h_sat + gamma * h_moist
                    soil_weight += (g_sat * h_sat + gamma * h_moist) * dx
                    if material_breakdown:
                        _hs[p_idx] = h_sat
                        _hm[p_idx] = h_moist
                        _hys[p_idx] = h_sat * (y_split + overlap_bot) / 2
                        _hym[p_idx] = h_moist * (overlap_top + y_split) / 2
                else:
                    sum_gam_h_y += h * gamma * (overlap_top + overlap_bot) / 2
                    sum_gam_h += h * gamma
                    soil_weight += h * gamma * dx
                    if material_breakdown:
                        _hm[p_idx] = h
                        _hym[p_idx] = h * (overlap_top + overlap_bot) / 2
                # Base material = the DEEPEST present layer (smallest base
                # elevation), independent of polygon iteration order. The
                # previous "last h>0 layer wins" relied on polygons being
                # stored top-to-bottom; an out-of-order block silently bound
                # the wrong base strength.
                if overlap_bot < base_overlap_bot:
                    base_overlap_bot = overlap_bot
                    base_material_idx = mat_index

        # Center of gravity
        y_cg = (sum_gam_h_y) / sum_gam_h if sum_gam_h > 0 else None

        # Distributed load
        def _resultant(funcs):
            """(D, d_x, d_y, qL, qR) for one direction's worth of load blocks on this
            slice top. Nonzero center: distinguish a real load from the case where the
            load starts/ends on a side (allows negative uplift/suction loads)."""
            if not funcs:
                qL_ = qR_ = 0
            else:
                qC_ = sum(f(x_c) for f in funcs)          # intensity at center
                if qC_ == 0:
                    qL_ = qR_ = 0
                else:
                    qL_ = sum(f(x_l) for f in funcs)      # intensity at left-top corner
                    qR_ = sum(f(x_r) for f in funcs)      # intensity at right-top corner
            return calc_dload_resultant(x_l, y_lt, x_r, y_rt, qL_, qR_, dl_top) \
                + (qL_, qR_)

        def _combine(normal_funcs, vert_funcs):
            """One equivalent resultant (D, d_x, d_y, beta_eff) for this slice top.

            The surface-normal blocks act at the top-edge inclination `beta`; the
            vertical blocks act straight down (beta = 0) with the SAME resultant
            magnitude — a gravity surcharge, the vendor's 'type: vertical'. Summing
            the two force vectors and re-expressing the sum as one magnitude and one
            inclination is exact, because every solver consumes the load only as
            (D sin beta, -D cos beta) at the point (d_x, d_y).

            With no vertical block the normal result is returned UNTOUCHED (not
            round-tripped through atan2), so the historical path is bit-identical.
            """
            Dn, dxn, dyn, qLn, qRn = _resultant(normal_funcs)
            if not vert_funcs:
                return Dn, dxn, dyn, beta, qLn, qRn
            Dv, dxv, dyv, qLv, qRv = _resultant(vert_funcs)
            qL_, qR_ = qLn + qLv, qRn + qRv
            if Dv == 0.0:
                return Dn, dxn, dyn, beta, qL_, qR_
            if Dn == 0.0:
                return Dv, dxv, dyv, 0.0, qL_, qR_
            b = radians(beta)
            Fx = Dn * sin(b)
            Fy = -Dn * cos(b) - Dv
            D = sqrt(Fx * Fx + Fy * Fy)
            w = abs(Dn) + abs(Dv)
            return (D,
                    (abs(Dn) * dxn + abs(Dv) * dxv) / w,
                    (abs(Dn) * dyn + abs(Dv) * dyv) / w,
                    degrees(atan2(Fx, -Fy)), qL_, qR_)

        dload, d_x, d_y, beta_d, qL, qR = _combine(dload_interp_funcs,
                                                   dload_vert_funcs)

        # Second distributed load (rapid-drawdown stage 2)
        dload2, d_x2, d_y2, beta_d2, qL2, qR2 = _combine(dload2_interp_funcs,
                                                         dload2_vert_funcs)

        # Seismic force
        kw = abs(k_seismic) * soil_weight

        # === BEGIN : "Tension crack water force" ===

        # By default, zero out t and its line‐of‐action:
        t_force = 0.0
        y_t_loc  = 0.0

        # Only nonzero for the appropriate end‐slice:
        if tcrack_water is not None and tcrack_water > 0:
            # Horizontal resultant of water in tension crack (triangular distribution):
            #    t = (1/2) * γ_w * (d_tc)^2
            # Here, gamma_w is the unit weight of the crack‐water (y_w),
            # and tcrack_water is the depth of water in the crack (d_tc).
            t_force = 0.5 * gamma_w * (tcrack_water ** 2)

            if right_facing:
                # Right‐facing slope → water pushes on left side of the first slice (i == 0)
                if i == 0:
                    # line of action is d_tc/3 above the bottom left corner y_lb.
                    # t_force is stored as a positive magnitude (like the seismic kw);
                    # each method resolves its own sliding-direction sign. The orientation-
                    # normalized methods (oms/bishop/janbu/corps/lowe) treat it as driving;
                    # spencer, which works in real coordinates, flips it for right-facing.
                    y_t_loc = y_lb + (tcrack_water / 3.0)
                else:
                    # other slices = no tension‐crack force
                    t_force = 0.0
                    y_t_loc = 0.0

            else:
                # Left‐facing slope → water pushes on right side of the last slice (i == n-1)
                if i == (len(all_xs) - 2):  # last slice index = (number_of_slices − 1)
                    # line of action is d_tc/3 above the bottom right corner y_rb
                    y_t_loc = y_rb + (tcrack_water / 3.0)
                else:
                    t_force = 0.0
                    y_t_loc = y_rb
        # === END: "Tension crack water force" ===

        # === BEGIN : "Reinforcement lines" ===

        # 1) Build this slice's base as a LineString from (x_l, y_lb) to (x_r, y_rb):
        slice_base = LineString([(x_l, y_lb), (x_r, y_rb)])

        # 2) For each reinforcement line crossing this base, evaluate the available
        # tension and route it by the line's Dir/Appl settings:
        #   tangent + active  -> the classical scalar p (the pre-v12 behavior)
        #   tangent + passive -> scalar p_pt (methods apply it with psi = alpha,
        #                        factored by FS)
        #   axial             -> component sums at angle psi (the line's own
        #                        inclination) with the crossing point folded into
        #                        first-moment sums, so any moment center works:
        #                        M about (Xo,Yo) = Yo*cx - my + mx - Xo*cy
        p_sum = 0.0        # tangent, active (legacy column)
        p_pt_sum = 0.0     # tangent, passive
        pa_cx = pa_cy = pa_mx = pa_my = 0.0   # axial, active
        pp_cx = pp_cy = pp_mx = pp_my = 0.0   # axial, passive
        for rl in reinf_lines_data:
            intersec = slice_base.intersection(rl["geom"])
            if intersec.is_empty:
                continue

            # Since we guarantee only one intersection point, it must be a Point:
            if not isinstance(intersec, Point):
                # (In the extremely unlikely case that intersection is not a Point,
                #  skip it. Our assumption is only one Point per slice-base.)
                continue

            # A crossing that lands exactly on a shared slice corner is returned
            # by shapely to BOTH adjacent bases; without a claim check the line's
            # tension is counted twice (measured: Bishop 1.679 -> 1.998 on vp030a
            # with the geosynthetic moved onto the corner). First base to see the
            # point keeps it.
            pt_key = (round(intersec.x, 9), round(intersec.y, 9))
            if pt_key in rl["claimed"]:
                continue
            rl["claimed"].add(pt_key)

            if "envelope" in rl:
                s_along = rl["geom"].project(intersec)
                t_mx, lp1, lp2, te1, te2 = rl["envelope"]
                t_i = rl["avail"](s_along, rl["length"] - s_along,
                                  t_mx, lp1, lp2, te1, te2,
                                  pullout=rl.get("pullout"))
            else:
                # legacy point-list path (hand-built slope_data): interp on X
                t_i = np.interp(intersec.x, rl["xs"], rl["ts"], left=0.0, right=0.0)
            if t_i <= 0.0:
                continue

            if rl["dir"] == "tangent":
                if rl["appl"] == "active":
                    p_sum += t_i
                else:
                    p_pt_sum += t_i
            else:
                # Axial: force at the line's own inclination, applied at the
                # crossing point. Components are stored in the REAL frame
                # (cos(psi) > 0); each solution method applies its own
                # right-facing flip, exactly as it does for the tangent force.
                psi_i = rl["psi"]
                ci, si = cos(psi_i), sin(psi_i)
                if rl["appl"] == "active":
                    pa_cx += t_i * ci
                    pa_cy += t_i * si
                    pa_mx += t_i * si * intersec.x
                    pa_my += t_i * ci * intersec.y
                else:
                    pp_cx += t_i * ci
                    pp_cy += t_i * si
                    pp_mx += t_i * si * intersec.x
                    pp_my += t_i * ci * intersec.y

        # Now p_sum is the TOTAL T‐pull acting at this slice's base.
        # === END: "Reinforcement lines" ===

        # === BEGIN: "Pile lines" ===
        h_pile = 0.0
        h_pile_pas = 0.0   # passive (Appl) portion of h_pile, factored by FS
        theta_p_val = 0.0
        x_pile = 0.0
        y_pile = 0.0
        for pl in pile_lines_data:
            intersec = slice_base.intersection(pl["geom"])
            if intersec.is_empty:
                continue
            if isinstance(intersec, Point):
                # Same shared-corner claim check as the reinforcement lines
                # above: a crossing exactly on a slice corner is returned to
                # both adjacent bases and would count the pile force twice.
                pt_key = (round(intersec.x, 9), round(intersec.y, 9))
                if pt_key in pl["claimed"]:
                    continue
                pl["claimed"].add(pt_key)
                pile_H = pl["H"]
                pile_H_was_auto = False
                F_pile_single = None
                ito_segments = None

                # Auto-compute H using Ito & Matsui if H is not specified but D and S are
                if pile_H is None and pl["D_pile"] is not None and pl["S"] is not None:
                    # Ito & Matsui is only valid for vertical piles
                    pile_coords = pl["geom"].coords
                    if abs(pile_coords[0][0] - pile_coords[1][0]) > 1e-6:
                        raise ValueError(
                            f'Ito & Matsui auto-computation requires vertical piles. '
                            f'Pile "{pl["label"]}" is battered '
                            f'(x1={pile_coords[0][0]}, x2={pile_coords[1][0]}). '
                            f'Specify H directly for battered piles.')
                    from .ito_matsui import intersect_pile_with_materials, compute_ito_matsui_force
                    gs_coords = np.array(ground_surface.coords)
                    gs_coords = gs_coords[np.argsort(gs_coords[:, 0])]  # np.interp needs ascending x
                    y_ground_at_pile = np.interp(intersec.x, gs_coords[:, 0], gs_coords[:, 1])
                    ito_segments = intersect_pile_with_materials(
                        intersec.x, y_ground_at_pile, intersec.y,
                        polygons, materials
                    )
                    pile_H, F_pile_single = compute_ito_matsui_force(pl["D_pile"], pl["S"], ito_segments)
                    pile_H_was_auto = True
                    # Warn once if Ito & Matsui produces very large H
                    global _ito_matsui_warned
                    if not _ito_matsui_warned and pile_H > 50000:
                        depth = y_ground_at_pile - intersec.y
                        print(f'[WARNING] Ito & Matsui computed very large H={pile_H:.0f} '
                              f'(D={pl["D_pile"]}, S={pl["S"]}, depth={depth:.1f}). '
                              f'This may exceed the pile structural capacity. '
                              f'Consider specifying H directly or increasing pile spacing.')
                        _ito_matsui_warned = True
                elif pile_H is None:
                    pile_H = 0.0

                # Structural capacity check (V_cap and M_cap are per-pile values)
                V_cap = pl.get("V_cap")
                M_cap = pl.get("M_cap")
                if (V_cap is not None or M_cap is not None) and pile_H > 0 and pl["S"] is not None:
                    S_pile = pl["S"]
                    if F_pile_single is None:
                        F_pile_single = pile_H * S_pile  # user-specified H
                    F_capped = F_pile_single

                    if V_cap is not None:
                        F_capped = min(F_capped, V_cap)

                    if M_cap is not None:
                        if pile_H_was_auto and ito_segments:
                            from .ito_matsui import compute_ito_matsui_force_and_moment_arm
                            _, _, L_m = compute_ito_matsui_force_and_moment_arm(
                                pl["D_pile"], S_pile, ito_segments)
                        else:
                            gs_coords = np.array(ground_surface.coords)
                            gs_coords = gs_coords[np.argsort(gs_coords[:, 0])]  # np.interp needs ascending x
                            y_gnd = np.interp(intersec.x, gs_coords[:, 0], gs_coords[:, 1])
                            depth = y_gnd - intersec.y
                            L_m = depth / 3.0 if depth > 0 else 0.0
                        if L_m > 0:
                            F_capped = min(F_capped, M_cap / L_m)

                    pile_H = F_capped / S_pile

                # The crossing's diagnostics, for the run summary. 'governed'
                # names which bound produced the applied force.
                _rec = {"label": pl["label"], "x": float(intersec.x),
                        "y": float(intersec.y), "computed": bool(pile_H_was_auto),
                        "S": pl.get("S"), "appl": pl["appl"],
                        "H_width": float(pile_H), "F_soil": None, "F_used": None,
                        "governed": None, "L_m": None, "depth": None}
                try:
                    _gs = np.array(ground_surface.coords)
                    _gs = _gs[np.argsort(_gs[:, 0])]
                    _rec["depth"] = float(np.interp(intersec.x, _gs[:, 0],
                                                    _gs[:, 1]) - intersec.y)
                except Exception:
                    pass
                if F_pile_single is not None:
                    _rec["F_soil"] = float(F_pile_single)
                    _F_used = pile_H * pl["S"] if pl.get("S") else None
                    _rec["F_used"] = None if _F_used is None else float(_F_used)
                    if _F_used is not None:
                        if _F_used >= F_pile_single - 1e-9:
                            _rec["governed"] = "soil force"
                        else:
                            # L_m was computed in the capacity block above
                            # whenever this row carries an M_cap.
                            _Vb = V_cap if V_cap is not None else float("inf")
                            _Mb = (M_cap / L_m) if (M_cap is not None and L_m > 0) \
                                else float("inf")
                            if _Mb <= _Vb:
                                _rec["governed"] = "bending (Mcap/Lm)"
                                _rec["L_m"] = float(L_m)
                            else:
                                _rec["governed"] = "shear (Vcap)"
                pile_report.append(_rec)

                h_pile += pile_H
                if pl["appl"] == "passive":
                    h_pile_pas += pile_H
                theta_p_val = pl["theta_p"]  # last pile's angle if multiple (unusual)
                x_pile = intersec.x
                y_pile = intersec.y
        # === END: "Pile lines" ===

        # === BEGIN: "Line loads" ===
        # Fold each line load acting on this slice's top into an equivalent
        # distributed-load resultant: magnitude ll_res at an equivalent
        # inclination ll_beta (the dload convention, force = (F sin b, -F cos b))
        # acting at (ll_x_pt, ll_y_pt). The solvers then reuse the proven D-term
        # sign machinery verbatim. The horizontal component is mirrored for
        # right-facing slopes, exactly as the top-edge beta is above. A straight-
        # down load (angle = -90) maps to ll_beta = 0.
        ll_h = ll_v = 0.0
        ll_wx = ll_wy = ll_wsum = 0.0
        _is_last_slice = (i == len(all_xs) - 2)
        for _ll in line_loads_data:
            _in = (x_l - 1e-9) <= _ll["x"] < x_r or \
                  (_is_last_slice and abs(_ll["x"] - x_r) <= 1e-9)
            if not _in:
                continue
            _delta = radians(_ll["angle"])
            _h = _ll["P"] * cos(_delta)
            if right_facing:
                _h = -_h
            ll_h += _h
            ll_v += _ll["P"] * sin(_delta)
            ll_wx += _ll["P"] * _ll["x"]
            ll_wy += _ll["P"] * _ll["y"]
            ll_wsum += _ll["P"]
        if ll_wsum > 0.0:
            ll_res = sqrt(ll_h * ll_h + ll_v * ll_v)
            ll_beta = degrees(atan2(ll_h, -ll_v))
            ll_x_pt = ll_wx / ll_wsum
            ll_y_pt = ll_wy / ll_wsum
        else:
            ll_res = 0.0
            ll_beta = 0.0
            ll_x_pt = 0.0
            ll_y_pt = 0.0
        # === END: "Line loads" ===

        # Process piezometric line and pore pressures using pre-computed coordinates
        piezo_y = piezo_y_all[i]
        piezo_y2 = piezo_y2_all[i]
        
        hw = 0
        hw2 = 0
        u = 0
        u2 = 0
        # Signed base pore pressure BEFORE the u >= 0 clamp. Stays 0 unless the
        # pore-pressure source produces a negative (suction) value; consumed only
        # by the opt-in apparent-cohesion suction option, never by the effective-
        # normal term (which always sees the clamped u below).
        u_unclamped = 0.0
        # Determine pore pressure method from material property
        mat_u = materials[base_material_idx]['u'] if base_material_idx is not None else 'none'
        if mat_u == 'none':
            u = 0
            u2 = 0
        elif mat_u == 'piezo':
            # This slice really does sample the line -- so there must BE one, and
            # this slice must lie inside it.
            if not _pz_have:
                raise _no_piezo_line_error(materials[base_material_idx]['name'])
            if _pz_ext is not None and np.isnan(piezo_y):
                raise _piezo_extent_error(
                    i, x_c, materials[base_material_idx]['name'], _pz_ext,
                    "piezometric line")
            if _pz2_ext is not None and np.isnan(piezo_y2):
                raise _piezo_extent_error(
                    i, x_c, materials[base_material_idx]['name'], _pz2_ext,
                    "second piezometric line")
            if not np.isnan(piezo_y) and piezo_y > y_cb:
                hw = piezo_y - y_cb
            if not np.isnan(piezo_y2) and piezo_y2 > y_cb:
                hw2 = piezo_y2 - y_cb
            u = hw * gamma_w if not np.isnan(piezo_y) else 0
            u2 = hw2 * gamma_w if not np.isnan(piezo_y2) else 0
            # Signed head to the piezometric line: negative above the line, where
            # its magnitude is the hydrostatic matric suction. No phreatic cos^2
            # correction (that models steady parallel seepage below the line).
            if not np.isnan(piezo_y):
                u_unclamped = (piezo_y - y_cb) * gamma_w
            # Lines declared Type='phreatic' on the piezo sheet get the
            # phreatic-inclination correction (XSTABL / Slide "Hu: auto"):
            # for steady seepage roughly parallel to an inclined phreatic
            # surface, the head at the base is reduced by cos^2(theta) of
            # the LOCAL phreatic-line slope above the slice.
            def _cos2(line, xc):
                if not line or len(line) < 2:
                    return 1.0
                xs = [p[0] for p in line]; ys = [p[1] for p in line]
                for k in range(len(xs) - 1):
                    if xs[k] <= xc <= xs[k + 1] and xs[k + 1] > xs[k]:
                        m = (ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k])
                        return 1.0 / (1.0 + m * m)
                return 1.0
            if slope_data.get('piezo_phreatic'):
                u *= _cos2(slope_data.get('piezo_line'), x_c)
            if slope_data.get('piezo_phreatic2'):
                u2 *= _cos2(slope_data.get('piezo_line2'), x_c)
        elif mat_u == 'seep':
            # Seepage-based pore pressure calculation using mesh interpolation.
            # A missing mesh or field is REFUSED here, not read as u = 0 -- the same
            # doctrine as the u='piezo' no-line branch above. preflight's
            # seep_field.missing rule catches this for a normal run, but preflight is
            # an INPUT check callers are allowed to skip (sensitivity sweeps and the
            # GeoStudio import probe both pass check_inputs=False), and load_slope_data
            # degrades to mesh=None / seep_u=None on a companion-file read error with
            # only a printed WARNING. That combination let a seep model solve
            # completely dry: on the right-facing seep sample it returns FS 3.18
            # against the wet 2.08, a +53% non-conservative error with nothing raised.
            if slope_data.get('mesh') is not None and slope_data.get('seep_u') is not None:
                mesh = slope_data['mesh']
                seep_u = slope_data['seep_u']

                # Interpolate pore pressure at the slice center base point.
                # signed=True delivers the RAW field including negative (suction)
                # values above the water table; this branch then applies its OWN
                # clamp for the effective-normal u, mirroring the piezo branch.
                point = (x_c, y_cb)
                u_val, found = interpolate_at_point(
                    mesh['nodes'],
                    mesh['elements'],
                    mesh['element_types'],
                    seep_u,
                    point,
                    return_found=True,
                    signed=True
                )
                _seep_pts += 1
                if not found:
                    _seep_outside += 1
                    warnings.warn(
                        "Seepage pore pressure: a slice base point fell outside the "
                        "seepage mesh and was assigned u = 0. Check that the mesh "
                        "spans the full depth of the failure surface (a value of 0 "
                        "below the phreatic surface over-predicts the factor of safety).")
                # Unsaturated seepage solutions carry negative u (suction) above the
                # water table; keep the signed value for the suction option before
                # clamping the effective-normal pore pressure at 0. max(0, signed)
                # is identical to the old max(0, clamped), so the effective-normal
                # u is byte-unchanged for every non-suction caller.
                u_unclamped = u_val
                u = max(0.0, u_val)
            else:
                raise _no_seep_field_error(
                    materials[base_material_idx]['name'],
                    'mesh' if slope_data.get('mesh') is None
                    else 'nodal pore-pressure field')

            # Check for second seep solution (rapid drawdown)
            if slope_data.get('seep_u2') is not None:
                mesh = slope_data['mesh']
                seep_u2 = slope_data['seep_u2']

                # Interpolate pore pressure at the slice center base point
                point = (x_c, y_cb)
                u2_val, found2 = interpolate_at_point(
                    mesh['nodes'],
                    mesh['elements'],
                    mesh['element_types'],
                    seep_u2,
                    point,
                    return_found=True
                )
                if not found2:
                    warnings.warn(
                        "Seepage pore pressure (second solution): a slice base point "
                        "fell outside the seepage mesh and was assigned u = 0.")
                u2 = max(0.0, u2_val)
            else:
                u2 = 0
        elif mat_u == 'ru':
            # Pore pressure ratio (template v12): u = ru * sigma_v, where sigma_v is
            # the soil-column vertical stress at the base center. sum_gam_h is
            # sum(gamma_i * h_i) over the material bands of this slice, identically
            # W/dx, so the quantity already exists. By definition (Bishop &
            # Morgenstern) sigma_v is the SOIL column only -- distributed loads and
            # crack water are excluded. There is no staged variant: u2 = u.
            ru_ratio = materials[base_material_idx].get('ru', 0.0) if base_material_idx is not None else 0.0
            u = ru_ratio * sum_gam_h
            u2 = u
        else:
            # Unreachable: load_slope_data validates u against this same set. A bare
            # `u = 0` here silently deleted pore pressure whenever the two lists drifted.
            raise ValueError(
                f"Material '{materials[base_material_idx]['name']}' has an unrecognized "
                f"pore pressure option u='{mat_u}'. Expected one of: none, piezo, seep, ru."
            )

        # Calculate alpha (slope angle of the failure surface) more efficiently
        delta = 0.01
        if use_arc:
            # For circular failure surface, use parametric equation for derivative
            # The slope at any point on the circle is: dy/dx = (x - Xo) / sqrt(R^2 - (x - Xo)^2)
            dx_circle = x_c - Xo
            if abs(dx_circle) < R:  # Check if point is on the circle
                alpha = degrees(atan(dx_circle / sqrt(R**2 - dx_circle**2)))
            else:
                # Fallback to numerical method
                y1 = get_circular_y_coordinates([x_c - delta], Xo, Yo, R)[0]
                y2 = get_circular_y_coordinates([x_c + delta], Xo, Yo, R)[0]
                alpha = degrees(atan2(y2 - y1, 2 * delta))
        elif comp is not None:
            # Exact on both branches: the circle equation on the arc, the floor
            # segment's own slope on the floor. Every slice lies wholly on one
            # branch (the crossings are slice boundaries), so this never averages
            # across the kink the way a +/-delta finite difference would.
            alpha = comp.alpha_deg(x_c)
        else:
            # For non-circular failure surface, use geometric intersection
            # Clamp the probe inside the surface span: an end slice narrower than
            # delta otherwise sends a probe past the surface, the lookup returns
            # None, and the fallback scores a near-vertical end segment as flat
            # with dl = its sliver width — the base resistance vanishes and a
            # search chases the artifact (vp103d measured 6% low; the same bug
            # let acads' 89 deg scarp under the max_base_angle guard).
            failure_line = clipped_surface
            xs_min, _, xs_max, _ = failure_line.bounds
            xl = max(x_c - delta, xs_min)
            xr = min(x_c + delta, xs_max)
            y1 = get_y_from_intersection(
                failure_line.intersection(LineString([(xl, -1e6), (xl, 1e6)])))
            y2 = get_y_from_intersection(
                failure_line.intersection(LineString([(xr, -1e6), (xr, 1e6)])))
            if y1 is not None and y2 is not None and xr > xl:
                alpha = degrees(atan2(y2 - y1, xr - xl))
            else:
                alpha = 0

        if right_facing:
            alpha = -alpha
        dl = dx / cos(radians(alpha))

        if base_material_idx is None:
            phi = 0
            c = 0
            c1 = 0      # not used in rapid drawdown, but must be defined
            phi1 = 0    # not used in rapid drawdown, but must be defined
            d = 0       # not used in rapid drawdown, but must be defined
            psi = 0     # not used in rapid drawdown, but must be defined
            pow_flag = False
            hb_flag = False
        else:
            pow_flag = False
            hb_flag = False
            mat_option = materials[base_material_idx]['option']
            if mat_option == 'hb':
                # Generalized Hoek-Brown. Like 'pow' this is a curved envelope, so
                # the strength depends on the base normal stress, which depends on
                # FS: solve._with_nonlinear_strength re-linearizes it into an
                # instantaneous tangent (c_i, phi_i) each outer iteration. Seed
                # that iteration here with a no-FS normal-stress estimate.
                mat = materials[base_material_idx]
                sigma0 = max(0.0, sum_gam_h - u)
                c_arr, phi_arr = hb_tangent(sigma0, mat['hb_sci'], mat['hb_gsi'],
                                            mat['hb_mi'], mat['hb_d'])
                c = float(c_arr)
                phi = float(phi_arr)
                c1 = 0; phi1 = 0; d = 0; psi = 0
                hb_flag = True
            elif mat_option == 'pow':
                # Power-curve envelope tau = a*(sigma'_n + d)^b + c_p. Strength
                # depends on the base normal stress, which depends on FS, so the
                # solvers iterate: each method carries an outer loop (see
                # solve._with_nonlinear_strength) that re-linearizes the curve at
                # the current sigma'_n into an instantaneous-tangent (c_i, phi_i).
                # Here we seed that iteration with a no-FS estimate of the normal
                # stress (infinite-slope form on the soil column).
                mat = materials[base_material_idx]
                pw_a, pw_b = mat['pow_a'], mat['pow_b']
                pw_c, pw_d = mat['pow_c'], mat['pow_d']
                sigma0 = max(0.0, sum_gam_h - u)
                s_eff0 = max(sigma0 + pw_d, 1e-4 * max(1.0, sigma0))
                slope0 = pw_a * pw_b * s_eff0 ** (pw_b - 1.0)
                tau0 = pw_a * s_eff0 ** pw_b + pw_c
                c = tau0 - sigma0 * slope0
                phi = degrees(atan(slope0))
                c1 = 0; phi1 = 0; d = 0; psi = 0
                pow_flag = True
            elif mat_option not in ('mc', 'cp'):
                # Was silently treated as 'cp'. A blank option is legal on seep-only
                # material rows, but not on one a failure surface passes through.
                raise ValueError(
                    f"Material '{materials[base_material_idx]['name']}' lies on the failure "
                    f"surface but has no valid strength option (option='{mat_option}'). "
                    "Expected one of: mc, cp."
                )
            elif mat_option == 'mc':
                c = materials[base_material_idx]['c']
                phi = materials[base_material_idx]['phi']
                c1 = c       # make a copy for use in rapid drawdown
                phi1 = phi   # make a copy for use in rapid drawdown
                d = materials[base_material_idx]['d']
                psi = materials[base_material_idx]['psi']
            else:
                # 'cp' option: undrained strength c at the reference elevation r_elev,
                # increasing by the rate cp per unit elevation below it (clamped to c
                # at/above r_elev): Su = c + cp * max(0, r_elev - y).
                mat = materials[base_material_idx]
                c = mat['c'] + max(0.0, mat['r_elev'] - y_cb) * mat['cp']
                phi = 0
                c1 = 0      # not used in rapid drawdown, but must be defined
                phi1 = 0    # not used in rapid drawdown, but must be defined
                d = 0       # not used in rapid drawdown, but must be defined
                psi = 0     # not used in rapid drawdown, but must be defined

        # Apparent cohesion from matric suction (Fredlund extended Mohr-Coulomb),
        # opt-in. For a base material named in suction_phi_b, convert the base
        # suction s = max(0, -u_unclamped) into an apparent cohesion
        # c_suction = s * tan(phi_b), optionally capping s at suction_cap. The
        # solvers add this to c in the resisting term c*dl; the effective-normal
        # term keeps the clamped u, so this is exactly the (u_a - u_w) tan(phi_b)
        # term with u_a = 0. Default (suction_phi_b None) => 0.0, bit-identical.
        c_suction = 0.0
        if suction_phi_b and base_material_idx is not None:
            _bm_name = materials[base_material_idx]['name']
            phi_b_deg = suction_phi_b.get(_bm_name)
            if phi_b_deg:
                s = max(0.0, -u_unclamped)
                # suction_cap may be a single scalar (one cap for every material) or
                # a per-material dict {name: cap}; a material absent from the dict --
                # or paired with None -- is uncapped.
                if suction_cap is not None:
                    _cap = (suction_cap.get(_bm_name)
                            if isinstance(suction_cap, dict) else suction_cap)
                    if _cap is not None:
                        s = min(s, _cap)
                c_suction = s * tan(radians(phi_b_deg))

        # Prepare slice data with conditional circle parameters
        slice_data = {
            'slice #': i + 1, # Slice numbering starts at 1
            'x_l': x_l,     # left x-coordinate of the slice
            'y_lb': y_lb,   # left y-coordinate of the slice base
            'y_lt': y_lt,   # left y-coordinate of the slice top
            'x_r': x_r,     # right x-coordinate of the slice
            'y_rb': y_rb,   # right y-coordinate of the slice base
            'y_rt': y_rt,   # right y-coordinate of the slice top
            'x_c': x_c,     # center x-coordinate of the slice
            'y_cb': y_cb,   # center y-coordinate of the slice base
            'y_ct': y_ct,   # center y-coordinate of the slice top
            'y_cg': y_cg,   # center of gravity y-coordinate of the slice
            'dx': dx,       # width of the slice
            'alpha': alpha,  # slope angle of the bottom of the slice in degrees
            'dl': dl,        # length of the slice along the failure surface
            **{f'h{j+1}': h for j, h in enumerate(heights)},  # heights of each layer in the slice
            'w': soil_weight,  # weight of the slice
            'qL': qL,  # distributed load intensity at left edge
            'qR': qR,  # distributed load intensity at right edge
            'dload': dload,     # distributed load resultant (area of trapezoid)
            'd_x': d_x, # dist load resultant x-coordinate (point d)
            'd_y': d_y, # dist load resultant y-coordinate (point d)
            'qL2': qL2,  # second distributed load intensity at left edge
            'qR2': qR2,  # second distributed load intensity at right edge
            'dload2': dload2,     # second distributed load resultant (area of trapezoid)
            'd_x2': d_x2, # second dist load resultant x-coordinate (point d)
            'd_y2': d_y2, # second dist load resultant y-coordinate (point d)
            # Inclination of the distributed-load resultant, in degrees. Equal to the
            # top-edge slope angle for a surface-normal load (the historical meaning
            # and still the default); 0 for a vertical (gravity-surcharge) load; and
            # the inclination of the vector sum when a slice carries both.
            'beta': beta_d,
            'beta2': beta_d2,  # ditto for the stage-2 load (see advanced.rapid_drawdown)
            'kw': kw,   # seismic force
            't': t_force,  # tension crack water force
            'y_t': y_t_loc,  # y-coordinate of the tension crack water force line of action
            'p': p_sum,   # sum of reinforcement line T values that intersect base of slice (tangent, active).
            'p_pt': p_pt_sum,   # tangent, passive reinforcement (factored by FS)
            'pa_cx': pa_cx,   # axial active reinforcement: sum T*cos(psi)
            'pa_cy': pa_cy,   # axial active: sum T*sin(psi)
            'pa_mx': pa_mx,   # axial active: sum T*sin(psi)*x_r (moment linearization)
            'pa_my': pa_my,   # axial active: sum T*cos(psi)*y_r
            'pp_cx': pp_cx,   # axial passive reinforcement components (factored by FS)
            'pp_cy': pp_cy,
            'pp_mx': pp_mx,
            'pp_my': pp_my,
            'lload': ll_res,   # line-load resultant on slice top (dload convention)
            'll_beta': ll_beta,   # equivalent inclination of the line-load resultant (deg)
            'll_x': ll_x_pt,   # line-load application point
            'll_y': ll_y_pt,
            'h_pile': h_pile,    # pile force magnitude per unit width (0 if no pile)
            'h_pile_pas': h_pile_pas,  # passive portion of h_pile (factored by FS)
            'theta_p': radians(theta_p_val),  # pile force angle from horizontal in radians
            'x_pile': x_pile,    # x-coordinate of pile-failure surface intersection
            'y_pile': y_pile,    # y-coordinate of pile-failure surface intersection
            'n_eff': 0, # Placeholder for effective normal force
            'z': 0,     # Placeholder for interslice side forces
            'theta': 0, # Placeholder for interslice angles
            'piezo_y': piezo_y,  # y-coordinate of the piezometric surface at x_c
            'piezo_y2': piezo_y2,  # y-coordinate of the piezometric surface at x_c for second piezometric line (rapid drawdown)
            'hw': hw,   # height of water at x_c
            'u': u,     # pore pressure at x_c
            'hw2': hw2, # height of water at x_c for second piezometric line (rapid drawdown)
            'u2': u2,   # pore pressure at x_c for second piezometric line (rapid drawdown)
            'mat': base_material_idx + 1 if base_material_idx is not None else None,  # index of the base material (1-indexed)
            'c': c,      # cohesion of the base material
            'c_suction': c_suction,  # apparent cohesion from matric suction (0.0 unless suction_phi_b opt-in)
            'phi': phi,
            'pow_flag': pow_flag,  # power-curve strength: iterate tangent (c,phi) on sigma'_n
            'pow_a': (materials[base_material_idx]['pow_a'] if pow_flag else 0.0),
            'pow_b': (materials[base_material_idx]['pow_b'] if pow_flag else 0.0),
            'pow_c': (materials[base_material_idx]['pow_c'] if pow_flag else 0.0),
            'pow_d': (materials[base_material_idx]['pow_d'] if pow_flag else 0.0),
            'hb_flag': hb_flag,  # Hoek-Brown strength: iterate tangent (c,phi) on sigma'_n
            'hb_sci': (materials[base_material_idx]['hb_sci'] if hb_flag else 0.0),
            'hb_gsi': (materials[base_material_idx]['hb_gsi'] if hb_flag else 0.0),
            'hb_mi': (materials[base_material_idx]['hb_mi'] if hb_flag else 0.0),
            'hb_d': (materials[base_material_idx]['hb_d'] if hb_flag else 0.0),
            'c1': c1,    # cohesion of the base material for rapid drawdown
            'phi1': phi1,  # friction angle of the base material for rapid drawdown
            'd': d,       # d cohesion of the base material for rapid drawdown
            'psi': psi,   # psi friction angle of the base material for rapid drawdown
        }
        
        # Add circle parameters only for circular failure surfaces
        if circular:
            slice_data.update({
                'r': R,      # radius of the circular failure surface
                'xo': Xo,    # x-coordinate of the center of the circular failure surface
                'yo': Yo,    # y-coordinate of the center of the circular failure surface
            })
        else:
            slice_data.update({
                'r': None,   # not applicable for non-circular failure surface
                'xo': None,  # not applicable for non-circular failure surface
                'yo': None,  # not applicable for non-circular failure surface
            })
        slices.append(slice_data)
        if material_breakdown:
            mb_Hm.append(_hm)
            mb_Hs.append(_hs)
            mb_HYm.append(_hym)
            mb_HYs.append(_hys)
            mb_base.append(-1 if base_material_idx is None else base_material_idx)

    # Every u='seep' base point outside the seepage mesh: the field was not read at
    # all, so the whole surface priced as dry. The per-point warning above cannot
    # carry this — warnings.warn de-duplicates by location, and batch runners
    # routinely filter warnings entirely, so the loudest possible signal for a
    # wholly-unread field was nothing at all. Refuse it instead.
    if _seep_pts and _seep_outside == _seep_pts:
        raise ValueError(
            f"Pore pressure: all {_seep_pts} slice base point(s) in a u='seep' "
            f"material fell outside the seepage mesh, so every one of them read "
            f"u = 0 and the surface was priced as if dry -- an unsafe, "
            f"non-conservative factor of safety. The mesh and the failure surface "
            f"do not overlap at all: check that the seepage solution belongs to "
            f"this model and that its mesh spans the failure surface.")

    df = pd.DataFrame(slices)

    # Surface a flat-arc facing note in the returned data (rather than guessing
    # silently). Only set on a symmetric flat arc with no override; None otherwise.
    df.attrs['facing_note'] = facing_note
    if pile_lines_data:
        df.attrs['pile_report'] = pile_report
    if material_breakdown:
        _mats_of_poly = []
        for pe in poly_edges:
            mat_id = pe['mat_id']
            _mats_of_poly.append(mat_id if (mat_id is not None
                                            and 0 <= mat_id < len(materials))
                                 else pe['poly_index'])
        df.attrs['mat_breakdown'] = {
            'Hm': np.asarray(mb_Hm), 'Hs': np.asarray(mb_Hs),
            'HYm': np.asarray(mb_HYm), 'HYs': np.asarray(mb_HYs),
            'base_idx': np.asarray(mb_base, dtype=int),
            'mat_of_poly': _mats_of_poly,
            'k_seismic': k_seismic,
            'suction': bool(suction_phi_b),
        }
    if facing_note and debug:
        print("WARNING: " + facing_note)

    # Hand solve() an explicit facing ONLY when it cannot be trusted to re-derive
    # the same one from geometry: a caller override, or a flat arc whose equal-
    # elevation ends make each solver's own y_lb/y_cb facing test degenerate
    # (and, worse, method-dependent). For a normal slope with no override this
    # key is left UNSET, so every solver falls back to its historical geometric
    # test and stays byte-identical — the flat-arc rule is inert everywhere else.
    if facing_override is not None or abs(y_left - y_right) < 1e-6:
        df.attrs['right_facing'] = bool(right_facing)

    # Slice data were built by iterating from left to right. Flip the order slice data for right-facing slopes.
    # Slice 1 should be at the bottom and slice n at the top. This makes the slice data consistent with the
    # sign convention for alpha and the free-body diagram used to calculate forces.
    # if right_facing:
    #     df = df.iloc[::-1].reset_index(drop=True)

    return True, (df, clipped_surface)