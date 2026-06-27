"""Map an Inputs-canvas double-click back to the slope_data object under it.

The Inputs canvas renders a rasterised pixmap (not a live Matplotlib canvas), so
Matplotlib pick events aren't available. MplCanvas maps the click to axes *data*
coordinates and this module hit-tests those against the input geometry, returning
the CATEGORY_EDITORS key to open (plan §6, §8) — or None when nothing is near.

Lines and points win over polygon interiors; clicking inside a material zone (no
line nearby) falls back to the materials editor for that zone.
"""

from __future__ import annotations

from shapely.geometry import LineString, Point


def _xy(p):
    """Coerce a geometry point to an (x, y) tuple. Points come as (x, y) tuples
    (profile/piezo) or {'X','Y'} dicts (dloads/non_circ)."""
    if isinstance(p, dict):
        return (p.get("X"), p.get("Y"))
    return (p[0], p[1])


def _line_dist(pt, points):
    """Distance from `pt` to the polyline through `points`, or inf if degenerate."""
    coords = [_xy(p) for p in points]
    coords = [c for c in coords if c[0] is not None and c[1] is not None]
    if len(coords) < 2:
        return float("inf")
    try:
        return LineString(coords).distance(pt)
    except Exception:
        return float("inf")


def pick_category(slope_data, x, y, tol):
    """Return the editor category for the input feature nearest to (x, y) within
    `tol` data units, or None. Geometry that isn't editable for this project type
    (e.g. polygons on a profile-based file) is simply not offered."""
    d = slope_data
    pt = Point(x, y)
    cands = []  # (distance, category)

    for p in d.get("pile_lines") or []:
        cands.append((_line_dist(pt, [(p["x1"], p["y1"]), (p["x2"], p["y2"])]), "piles"))
    for r in d.get("reinforcement_lines") or []:
        cands.append((_line_dist(pt, [(r["x1"], r["y1"]), (r["x2"], r["y2"])]), "reinforce"))
    for line in d.get("dloads") or []:
        cands.append((_line_dist(pt, line), "dloads"))
    for key in ("piezo_line", "piezo_line2"):
        cands.append((_line_dist(pt, d.get(key) or []), "piezo"))
    for c in d.get("circles") or []:
        try:
            r = ((x - c["Xo"]) ** 2 + (y - c["Yo"]) ** 2) ** 0.5
            cands.append((abs(r - c["R"]), "circles"))
        except Exception:
            pass
    cands.append((_line_dist(pt, d.get("non_circ") or []), "non_circ"))

    # Geometry: profile lines (profile-based files) or polygon edges (polygon-based).
    profile_lines = d.get("profile_lines") or []
    polygons = d.get("polygons") or []
    if profile_lines:
        for pl in profile_lines:
            cands.append((_line_dist(pt, pl.get("coords") or []), "profile"))
    elif polygons:
        for poly in polygons:
            try:
                cands.append((poly["polygon"].exterior.distance(pt), "polygons"))
            except Exception:
                pass

    near = sorted((dist, cat) for dist, cat in cands if dist <= tol)
    if near:
        return near[0][1]

    # Fallback: a click inside a material zone edits that zone's material.
    for poly in polygons:
        try:
            if poly["polygon"].contains(pt):
                return "materials"
        except Exception:
            pass
    return None
