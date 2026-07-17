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


def pick_category(slope_data, x, y, tol, mode=None):
    """Return ``(category, index)`` for the input feature nearest to (x, y) within
    `tol` data units, or None. `index` is the feature's position in its list (so
    the editor can pre-highlight it); it's None when not meaningful. Geometry that
    isn't editable for this project type (e.g. polygons on a profile-based file) is
    simply not offered.

    `mode` ('lem'/'seep'/'fem') matches the picks to what the Inputs plot actually
    draws: seepage BCs only in seep mode, distributed loads only outside it (the
    two otherwise overlap along the surface). None includes everything."""
    d = slope_data
    pt = Point(x, y)
    cands = []  # (distance, category, index)

    for i, p in enumerate(d.get("pile_lines") or []):
        cands.append((_line_dist(pt, [(p["x1"], p["y1"]), (p["x2"], p["y2"])]), "piles", i))
    for i, r in enumerate(d.get("reinforcement_lines") or []):
        cands.append((_line_dist(pt, [(r["x1"], r["y1"]), (r["x2"], r["y2"])]), "reinforce", i))
    # Line loads are single points on the ground surface (drawn as an arrow); a
    # point hit-test picks the one nearest the click.
    for i, ll in enumerate(d.get("line_loads") or []):
        try:
            cands.append((pt.distance(Point(ll["x"], ll["y"])), "line_loads", i))
        except Exception:
            pass
    if mode != "seep":
        # Both sets are drawn (purple / orange); index is (set, block) so the
        # editor can jump to the right tab + load.
        for s, key in enumerate(("dloads", "dloads2")):
            for i, line in enumerate(d.get(key) or []):
                cands.append((_line_dist(pt, line), "dloads", (s, i)))
    if mode is None or mode == "seep":
        # index is (set, row), matching the editor's list order: heads (0…), then
        # fluxes, then the exit face as the trailing row.
        for s, key in enumerate(("seepage_bc", "seepage_bc2")):
            bc = d.get(key) or {}
            heads = bc.get("specified_heads") or []
            fluxes = bc.get("specified_fluxes") or []
            for j, sh in enumerate(heads):
                cands.append((_line_dist(pt, sh.get("coords") or []), "seep_bc", (s, j)))
            for j, sf in enumerate(fluxes):
                cands.append((_line_dist(pt, sf.get("coords") or []),
                              "seep_bc", (s, len(heads) + j)))
            cands.append((_line_dist(pt, bc.get("exit_face") or []),
                          "seep_bc", (s, len(heads) + len(fluxes))))
    for key in ("piezo_line", "piezo_line2"):
        cands.append((_line_dist(pt, d.get(key) or []), "piezo", None))
    for i, c in enumerate(d.get("circles") or []):
        try:
            r = ((x - c["Xo"]) ** 2 + (y - c["Yo"]) ** 2) ** 0.5
            cands.append((abs(r - c["R"]), "circles", i))
        except Exception:
            pass
    cands.append((_line_dist(pt, d.get("non_circ") or []), "non_circ", None))

    # Geometry: profile lines (profile-based files) or polygon edges (polygon-based).
    profile_lines = d.get("profile_lines") or []
    polygons = d.get("polygons") or []
    if profile_lines:
        for i, pl in enumerate(profile_lines):
            cands.append((_line_dist(pt, pl.get("coords") or []), "profile", i))
        # The max-depth line (the model's bottom boundary, drawn across the profile
        # x-extent) also opens the profile editor.
        md = d.get("max_depth")
        if md is not None:
            xs = [x for pl in profile_lines for (x, _) in (pl.get("coords") or [])]
            if xs:
                cands.append((_line_dist(pt, [(min(xs), md), (max(xs), md)]), "profile", None))
    elif polygons:
        for i, poly in enumerate(polygons):
            try:
                cands.append((poly["polygon"].exterior.distance(pt), "polygons", i))
            except Exception:
                pass

    near = sorted((c for c in cands if c[0] <= tol), key=lambda c: c[0])
    if near:
        return near[0][1], near[0][2]

    # Fallback: a click inside a material zone edits that zone's material.
    for i, poly in enumerate(polygons):
        try:
            if poly["polygon"].contains(pt):
                return "materials", poly.get("mat_id")
        except Exception:
            pass
    return None
