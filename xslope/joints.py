"""Joint NETWORKS: generators that write the rows of the ``joints`` sheet.

A jointed rock mass is rarely described one line at a time. What a section
carries is a SET — bedding at a spacing and a dip, a conjugate set crossing it,
a blocky mass with no preferred orientation — and what the mesher needs is the
individual traces those sets resolve into, clipped to the ground they run
through. These three functions do that conversion, and nothing else: each
returns a list of joint-line dicts in exactly the form
``slope_data['joint_lines']`` holds and :func:`xslope.fileio.save_slope_data_to_xlsx`
writes to the sheet, so a generated network can be inspected, edited by hand,
saved and reloaded like any other input.

    :func:`parallel_set`   one set: a dip, a spacing, an optional persistence
    :func:`cross_jointed`  two sets that cross, checked for the things the
                           mesher refuses
    :func:`voronoi`        a blocky mass: a seeded Voronoi tessellation

**The set.** :class:`JointSet` is one description — a kind, its parameters, a
region — and :meth:`JointSet.generate` resolves it into the rows. Nothing keeps
the description: what a file holds is the joint lines themselves, which is what
the mesher, the solver and the reader all work on. A set is changed by removing
its rows and building another.

**The region.** Every generator clips its traces to a region, which is the whole
section by default, a named material (or several), a polygon of Type ``joints``
on the polygon sheet by name or number, or a polygon given outright — optionally
cut to an elevation band. A trace that leaves the region is cut at the boundary;
a trace that would lie ALONG the boundary is dropped rather than emitted, because
the mesh split needs material on both sides of every joint and there is none
outside the section. That is the one rule the caller cannot see coming, so it is
applied here rather than left to the mesher's refusal. The elevation band is the
exception to it: a band's edges are lines drawn through material, so a trace
along one of them is kept.

**Dip.** ``dip_deg`` is the trace's inclination measured counter-clockwise from
the positive x axis: 0 is horizontal, 90 is vertical, and a set that descends to
the right — dipping right, in the field sense — is a NEGATIVE angle. A set and
its negative are mirror images about the vertical, which is what makes a
conjugate pair ``cross_jointed(parallel_set(..., 60), parallel_set(..., -60))``.

**Labels.** Rows are labeled ``<label>-01``, ``<label>-02``, … in the order the
traces come out, which is along the set's own normal from the low side. The
label reaches the plots, the 1D details view and the report, so a network keeps
one name per set rather than a hundred anonymous lines — and :func:`set_name`
reads the set back off it, which is what lets the rows of one network be found
and removed together.

**Properties.** ``props`` is applied to every row of the set: any of ``c``,
``phi``, ``c_res``, ``phi_res``, ``dil``, ``t_cut``, ``kn``, ``ks`` and ``jred``,
with the sheet's own defaults for the rest. ``phi`` has no default on the sheet
and none here; a set generated without one is refused by preflight by name, as a
hand-entered row without one is.

Run a generated network through the mesher as it stands: the traces of two sets
meet, and the mesh split copies each shared node once per wedge of material
around it (see :mod:`xslope.joint`). Nothing about a generated line differs from
a typed one.
"""

from __future__ import annotations

import math

import numpy as np
from shapely.geometry import (GeometryCollection, LineString, MultiLineString,
                              MultiPolygon, Point, Polygon)
from shapely.ops import unary_union

#: Length below which a clipped trace is dropped as a sliver rather than emitted.
#: A segment shorter than this cannot be meshed into even one 1D element at any
#: reasonable size, and a zero-length one is a mesher failure.
MIN_TRACE_FRAC = 1.0e-3

#: How close to a vertex of the region a clipped endpoint may land before it is
#: SNAPPED onto it, as a fraction of the region's diameter.
#:
#: Clipping a trace against a polygon is floating-point arithmetic, and it can
#: put an endpoint a part in 10^15 off the corner it is meant to land on. That
#: leaves a sliver of boundary edge between the two, and gmsh's edge recovery
#: cannot mesh one: it splits the offending edges and tries again, over and over.
#: Snapping the endpoint onto the vertex it is already at removes the sliver; the
#: tolerance is six orders of magnitude above the rounding it repairs and many
#: below any geometry a user states.
SNAP_TOL_FRAC = 1.0e-9

#: How close to the region's boundary a trace may lie along its whole length
#: before it is read as running ALONG the boundary rather than through the
#: material, as a fraction of the region's diameter. The mesh split needs
#: material on BOTH sides of a joint, and a line on the outside of the section
#: has it on one; the mesher refuses such a line by name and this keeps a
#: generator from producing one.
BOUNDARY_TOL_FRAC = 1.0e-6

#: Every key a generated row may carry, with the value that means "the sheet's
#: own default". ``phi`` is NaN because the sheet requires it: a set generated
#: without one reaches preflight as a joint with no friction angle, named.
_ROW_DEFAULTS = {
    "c": 0.0,
    "phi": float("nan"),
    "c_res": float("nan"),
    "phi_res": float("nan"),
    "dil": float("nan"),
    "t_cut": 0.0,
    "kn": float("nan"),
    "ks": float("nan"),
    "jred": "",
}


# ---------------------------------------------------------------------------
# Region and row helpers
# ---------------------------------------------------------------------------

def resolve_region(slope_data, region=None, band=None):
    """The polygon a generator clips its traces to.

    ``region`` may be

    ``None``
        the whole section, ``slope_data['domain_polygon']``;
    a material name, or an iterable of them
        the union of the polygons carrying those materials — the way a bedding
        set confined to one rock unit is asked for;
    ``'mat:<name>'`` or ``'mat:<name>+<name>'``
        the same thing, said explicitly. This is the form the set record uses,
        where a bare word would be ambiguous;
    ``'poly:<name>'`` or ``'poly:<n>'``
        a polygon of Type ``joints`` on the polygon sheet, by the name in its
        block header or by its number among the joint regions (1-based). This is
        how a set is confined to ground no material boundary draws — one block of
        an outcrop, the part of a unit behind a wall;
    a shapely polygon, or an iterable of ``(x, y)``
        that polygon outright.

    ``band`` is an optional ``(y_min, y_max)`` elevation band, either end
    ``None`` for open: the region is cut to it, which is what "the sandstone
    above elevation 40" asks for. It intersects whatever the region resolved to.

    Returns a shapely ``Polygon`` or ``MultiPolygon``.
    """
    poly = _resolve_region_only(slope_data, region)
    if band is None:
        return poly
    return _clip_to_band(poly, band)


def _clip_to_band(poly, band):
    """A region cut to an elevation band, ``(y_min, y_max)`` with either end
    ``None`` for open."""
    from shapely.geometry import box as _box
    lo, hi = band
    minx, miny, maxx, maxy = poly.bounds
    pad = max(maxx - minx, maxy - miny, 1.0)
    y0 = miny - pad if lo is None else float(lo)
    y1 = maxy + pad if hi is None else float(hi)
    if y1 <= y0:
        raise ValueError(
            f"The elevation band {y0:g} to {y1:g} is empty: its upper elevation "
            f"is at or below its lower one.")
    cut = poly.intersection(_box(minx - pad, y0, maxx + pad, y1))
    keep = [g for g in (cut.geoms if hasattr(cut, "geoms") else [cut])
            if isinstance(g, (Polygon, MultiPolygon)) and not g.is_empty]
    out = unary_union(keep) if keep else Polygon()
    if out.is_empty:
        raise ValueError(
            f"No part of the region lies in the elevation band "
            f"{'open' if lo is None else f'{float(lo):g}'} to "
            f"{'open' if hi is None else f'{float(hi):g}'}. The region spans "
            f"elevations {miny:g} to {maxy:g}.")
    return out


def _joint_zone_polygon(slope_data, ref):
    """The polygon of a Type ``joints`` region, by name or by number."""
    zones = list((slope_data or {}).get("joint_zones") or [])
    if not zones:
        raise ValueError(
            "This model has no joint region: no polygon on the polygon sheet "
            "declares Type 'joints'. Draw one (Polygons editor, or the polygon "
            "sheet's Type cell) and name it in its block header, or give the "
            "region as a material name or a polygon.")
    ref = str(ref).strip()
    names = [str(z.get("label") or "").strip() for z in zones]
    if ref.isdigit():
        n = int(ref)
        if not 1 <= n <= len(zones):
            raise ValueError(
                f"This model has {len(zones)} joint region(s), so there is no "
                f"joint region {n}.")
        z = zones[n - 1]
    else:
        hits = [i for i, nm in enumerate(names) if nm.lower() == ref.lower()]
        if not hits:
            listed = ", ".join(repr(nm) if nm else f"#{i + 1} (unnamed)"
                               for i, nm in enumerate(names))
            raise ValueError(
                f"No joint region named {ref!r}. The model's joint regions are: "
                f"{listed}. A region's name is the polygon sheet's block header "
                f"over its column.")
        if len(hits) > 1:
            raise ValueError(
                f"{len(hits)} joint regions are named {ref!r}, so the name does "
                f"not say which one. Rename one of them in its block header on "
                f"the polygon sheet, or name the region by its number.")
        z = zones[hits[0]]
    coords = [(float(x), float(y)) for x, y in (z.get("polygon") or [])]
    if len(coords) < 3:
        raise ValueError(
            f"Joint region {ref!r} has {len(coords)} vertices, which is not a "
            f"polygon.")
    return Polygon(coords)


def _resolve_region_only(slope_data, region=None):
    """:func:`resolve_region` without the elevation band."""
    if region is None:
        dom = (slope_data or {}).get("domain_polygon")
        if dom is None or dom.is_empty:
            raise ValueError(
                "The model has no domain polygon to generate joints in. Pass "
                "region= a material name or a polygon, or load a model whose "
                "geometry has been built.")
        return dom
    if isinstance(region, (Polygon, MultiPolygon)):
        return region
    if isinstance(region, str):
        text = region.strip()
        if text.lower().startswith("poly:"):
            return _joint_zone_polygon(slope_data, text[5:])
        if text.lower().startswith("mat:"):
            region = [p.strip() for p in text[4:].split("+") if p.strip()]
        else:
            region = [region]
    region = list(region)
    if not region:
        raise ValueError("region names no material and no polygon.")
    if all(isinstance(r, str) for r in region):
        names = {r.strip().lower() for r in region}
        mats = (slope_data or {}).get("materials") or []
        ids = {i for i, m in enumerate(mats)
               if str(m.get("name", "")).strip().lower() in names}
        unknown = names - {str(m.get("name", "")).strip().lower() for m in mats}
        if unknown:
            raise ValueError(
                "No material named " + ", ".join(sorted(repr(u) for u in unknown))
                + ". The model's materials are: "
                + ", ".join(repr(str(m.get("name", ""))) for m in mats) + ".")
        polys = [p["polygon"] for p in (slope_data.get("polygons") or [])
                 if int(p.get("mat_id", -1)) in ids]
        if not polys:
            raise ValueError(
                "The named material(s) carry no polygon in this model, so there "
                "is no region to generate joints in.")
        merged = unary_union(polys)
        return merged
    try:
        return Polygon(region)
    except Exception as exc:                      # pragma: no cover - input guard
        raise ValueError(f"region is not a polygon, a material name or a list "
                         f"of (x, y) points: {exc}")


def _row_props(props):
    """The property dict every row of a set carries, defaults filled in."""
    out = dict(_ROW_DEFAULTS)
    for k, v in dict(props or {}).items():
        key = str(k).strip().lower()
        if key not in _ROW_DEFAULTS:
            raise ValueError(
                f"{k!r} is not a joint property. The joints sheet's columns are: "
                + ", ".join(sorted(_ROW_DEFAULTS)) + ".")
        out[key] = v
    return out


def _row(label, p1, p2, props):
    row = {"label": label,
           "x1": float(p1[0]), "y1": float(p1[1]),
           "x2": float(p2[0]), "y2": float(p2[1])}
    row.update(props)
    return row


def _boundary_of(region):
    """Every boundary ring of the region, as one geometry to measure against."""
    if isinstance(region, MultiPolygon):
        return unary_union([g.boundary for g in region.geoms])
    return region.boundary


def _region_vertices(region):
    """Every vertex of the region's boundary, as an (n, 2) array."""
    polys = region.geoms if isinstance(region, MultiPolygon) else [region]
    pts = []
    for g in polys:
        pts.extend(list(g.exterior.coords))
        for ring in g.interiors:
            pts.extend(list(ring.coords))
    return np.asarray(pts, dtype=float) if pts else np.zeros((0, 2))


def _snap(pt, verts, tol):
    """A clipped endpoint, moved onto the region vertex it is already at.

    See :data:`SNAP_TOL_FRAC`: an endpoint a rounding error away from a corner
    leaves a sliver of boundary edge the mesher cannot recover.
    """
    if len(verts) == 0 or tol <= 0.0:
        return pt
    d = np.hypot(verts[:, 0] - pt[0], verts[:, 1] - pt[1])
    i = int(np.argmin(d))
    return (float(verts[i, 0]), float(verts[i, 1])) if d[i] <= tol else pt


def _segments(clipped):
    """The ``LineString`` pieces of a shapely clip result, in order."""
    if clipped.is_empty:
        return []
    if isinstance(clipped, LineString):
        return [clipped]
    if isinstance(clipped, (MultiLineString, GeometryCollection)):
        return [g for g in clipped.geoms if isinstance(g, LineString)
                and not g.is_empty and g.length > 0.0]
    return []


def _keep(seg, boundary, min_len, bnd_tol):
    """Whether a clipped trace is a joint the mesher can split along.

    Two things disqualify it: it is shorter than one element, or it lies ALONG
    the region's boundary, where there is material on one side only. The second
    is measured on interior sample points rather than on the endpoints, because
    a trace ENDING on the boundary is perfectly ordinary — that is a crack tip
    reaching the face — and only a trace whose whole length hugs the boundary is
    the case the mesher refuses.
    """
    if seg.length <= min_len:
        return False
    for f in (0.1, 0.3, 0.5, 0.7, 0.9):
        p = seg.interpolate(f, normalized=True)
        if boundary.distance(Point(p)) > bnd_tol:
            return True
    return False


def _chop(seg, persistence):
    """One clipped trace cut into its persistent pieces.

    ``persistence`` is ``(trace_len, gap)``: a joint that exists for
    ``trace_len`` then stops for ``gap`` — a rock bridge — then resumes.
    ``None`` leaves the trace continuous.
    """
    if persistence is None:
        return [seg]
    trace_len, gap = persistence
    trace_len = float(trace_len)
    gap = float(gap)
    if trace_len <= 0.0:
        raise ValueError("persistence trace length must be positive.")
    if gap < 0.0:
        raise ValueError("persistence gap cannot be negative.")
    out = []
    s = 0.0
    L = seg.length
    while s < L:
        e = min(s + trace_len, L)
        if e - s > 0.0:
            out.append(LineString([seg.interpolate(s), seg.interpolate(e)]))
        s = e + gap
    return out


# ---------------------------------------------------------------------------
# The generators
# ---------------------------------------------------------------------------

def parallel_set(slope_data, dip_deg, spacing, offset=0.0, persistence=None,
                 region=None, label="set1", props=None, band=None):
    """One set of parallel joints, as the rows of the ``joints`` sheet.

    Parameters
    ----------
    slope_data : dict
        The model, read for its domain polygon and material polygons.
    dip_deg : float
        The traces' inclination, counter-clockwise from the positive x axis:
        0 horizontal, 90 vertical, negative for a set descending to the right.
    spacing : float
        The perpendicular distance between neighbouring traces.
    offset : float, optional
        Where the set sits across its own normal. One trace passes through the
        origin offset by this much, so ``offset=0`` puts a trace through
        ``(0, 0)`` extended and ``offset=spacing/2`` shifts the whole set half a
        spacing. Changing it moves the set; it does not change how many traces
        the region gets, except where the shift takes one across a corner.
    persistence : (trace_len, gap), optional
        A discontinuous set: each trace exists for ``trace_len``, stops for
        ``gap`` — a rock bridge — and resumes. ``None`` is a fully persistent
        set, every trace continuous across the region.
    region : optional
        Where the set exists; see :func:`resolve_region`.
    label : str, optional
        The set's name. Rows come out ``label-01``, ``label-02``, …
    props : dict, optional
        Joint properties applied to every row: ``c``, ``phi``, ``c_res``,
        ``phi_res``, ``dil``, ``t_cut``, ``kn``, ``ks``, ``jred``.
    band : (y_min, y_max), optional
        An elevation band the region is cut to, either end ``None`` for open.
        A trace lying along the band's own edge is KEPT, unlike one lying along
        the region's boundary: a band is a line drawn through material, and the
        mesh split has material on both sides of it.

    Returns
    -------
    list of dict
        Joint-line rows, ready for ``slope_data['joint_lines']``.
    """
    spacing = float(spacing)
    if spacing <= 0.0:
        raise ValueError("Joint spacing must be positive.")
    whole = _resolve_region_only(slope_data, region)
    poly = whole if band is None else _clip_to_band(whole, band)
    if poly.is_empty:
        raise ValueError("The region for this joint set is empty.")
    p = _row_props(props)

    th = math.radians(float(dip_deg))
    dx, dy = math.cos(th), math.sin(th)
    nx, ny = -dy, dx                              # the set's own normal

    minx, miny, maxx, maxy = poly.bounds
    diag = math.hypot(maxx - minx, maxy - miny)
    if diag <= 0.0:
        raise ValueError("The region for this joint set has no extent.")
    min_len = MIN_TRACE_FRAC * diag
    bnd_tol = BOUNDARY_TOL_FRAC * diag
    snap_tol = SNAP_TOL_FRAC * diag
    # "On the boundary" is measured against the region BEFORE the band cut it: the
    # band's own edges are lines drawn through material, and a trace along one of
    # them splits a mesh perfectly well.
    boundary = _boundary_of(whole)
    verts = _region_vertices(poly)

    # Where the region sits along the set's normal, so the traces cover it and
    # nothing beyond it.
    corners = [(minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)]
    s_vals = [cx * nx + cy * ny for cx, cy in corners]
    s_lo, s_hi = min(s_vals), max(s_vals)
    k_lo = math.floor((s_lo - float(offset)) / spacing)
    k_hi = math.ceil((s_hi - float(offset)) / spacing)

    rows = []
    for k in range(int(k_lo), int(k_hi) + 1):
        s = float(offset) + k * spacing
        # The infinite trace at this offset, drawn long enough to cross the region.
        mid = (s * nx, s * ny)
        half = diag
        line = LineString([(mid[0] - half * dx, mid[1] - half * dy),
                           (mid[0] + half * dx, mid[1] + half * dy)])
        for seg in _segments(line.intersection(poly)):
            for piece in _chop(seg, persistence):
                if not _keep(piece, boundary, min_len, bnd_tol):
                    continue
                a = _snap(piece.coords[0], verts, snap_tol)
                b = _snap(piece.coords[-1], verts, snap_tol)
                if math.hypot(b[0] - a[0], b[1] - a[1]) <= min_len:
                    continue
                rows.append(_row(f"{label}-{len(rows) + 1:02d}", a, b, p))
    return rows


def cross_jointed(set1, set2):
    """Two sets in one network, checked for what the mesher refuses.

    Concatenates the two lists and relabels nothing — each set keeps its own
    name, which is what tells the two apart in the plots and the report. What it
    adds is the check: two joint lines may MEET at a point, which is exactly what
    a cross-jointed mass is, but they may not lie ON one another over a stretch,
    and the mesher refuses a pair that does. A conjugate pair at two different
    dips never overlaps; one generated twice at the same dip does, on every
    trace, and that is the mistake this catches.

    Returns the combined row list.
    """
    set1 = list(set1)
    set2 = list(set2)
    if not set1 or not set2:
        raise ValueError("cross_jointed needs two non-empty joint sets.")
    names1 = {str(r.get("label")) for r in set1}
    clash = names1 & {str(r.get("label")) for r in set2}
    if clash:
        raise ValueError(
            "Both sets carry the label(s) " + ", ".join(sorted(clash))
            + ". Give the two sets different label= names, so each joint line "
              "names the set it belongs to.")
    tol = _overlap_tol(set1 + set2)
    for a in set1:
        la = LineString([(a["x1"], a["y1"]), (a["x2"], a["y2"])])
        for b in set2:
            lb = LineString([(b["x1"], b["y1"]), (b["x2"], b["y2"])])
            shared = la.intersection(lb)
            if isinstance(shared, LineString) and shared.length > tol:
                raise ValueError(
                    f"Joint lines {a['label']!r} and {b['label']!r} lie on one "
                    f"another over {shared.length:g} of length. Two joint lines "
                    f"may meet at a point but not share a stretch: the mesh "
                    f"split has no material between them there. The two sets "
                    f"need different dips.")
    return set1 + set2


def _overlap_tol(rows):
    """How much shared length reads as a crossing rather than an overlap."""
    xs = [r["x1"] for r in rows] + [r["x2"] for r in rows]
    ys = [r["y1"] for r in rows] + [r["y2"] for r in rows]
    diag = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
    return BOUNDARY_TOL_FRAC * max(diag, 1.0)


def voronoi(slope_data, block_size, seed, region=None, label="vor", props=None,
            band=None):
    """A blocky mass: a seeded Voronoi tessellation, as joint lines.

    A rock mass with no through-going set is described by its BLOCK SIZE rather
    than by a dip and a spacing, and the standard idealisation of it is a Voronoi
    tessellation — cells of about the stated size, in no preferred orientation,
    every cell wall a joint. The seed points are a jittered grid at that spacing,
    drawn from ``numpy``'s own generator with the stated ``seed``, so the same
    two numbers reproduce the same network exactly.

    Only the walls BETWEEN cells become joint lines. A cell wall that falls on
    the region's boundary is the outside of the section, which has material on
    one side only, and is dropped.

    Parameters
    ----------
    slope_data : dict
        The model, read for its region polygons.
    block_size : float
        The characteristic cell size: the seed spacing.
    seed : int
        The random seed. Required, not optional — a tessellation nobody can
        reproduce is not an input.
    region, label, props, band
        As :func:`parallel_set`.

    Returns
    -------
    list of dict
    """
    from scipy.spatial import Voronoi          # local: scipy is a heavy import

    block_size = float(block_size)
    if block_size <= 0.0:
        raise ValueError("Voronoi block size must be positive.")
    whole = _resolve_region_only(slope_data, region)
    poly = whole if band is None else _clip_to_band(whole, band)
    if poly.is_empty:
        raise ValueError("The region for this Voronoi network is empty.")
    p = _row_props(props)

    minx, miny, maxx, maxy = poly.bounds
    diag = math.hypot(maxx - minx, maxy - miny)
    min_len = MIN_TRACE_FRAC * diag
    bnd_tol = BOUNDARY_TOL_FRAC * diag
    snap_tol = SNAP_TOL_FRAC * diag
    boundary = _boundary_of(whole)             # see parallel_set: the band's edges
    verts = _region_vertices(poly)             # are not the section's boundary

    # A jittered grid, padded one cell beyond the region so the cells at its edge
    # are bounded by neighbours rather than by Voronoi's rays to infinity.
    rng = np.random.default_rng(int(seed))
    pad = 2.0 * block_size
    xs = np.arange(minx - pad, maxx + pad + block_size, block_size)
    ys = np.arange(miny - pad, maxy + pad + block_size, block_size)
    gx, gy = np.meshgrid(xs, ys)
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    pts = pts + rng.uniform(-0.35, 0.35, size=pts.shape) * block_size
    if len(pts) < 4:
        raise ValueError(
            "The block size is larger than the region, so there is no "
            "tessellation to make. Use a smaller block_size.")

    vor = Voronoi(pts)
    rows = []
    for (a, b) in vor.ridge_vertices:
        if a < 0 or b < 0:
            continue                              # a ridge running to infinity
        wall = LineString([vor.vertices[a], vor.vertices[b]])
        for seg in _segments(wall.intersection(poly)):
            if not _keep(seg, boundary, min_len, bnd_tol):
                continue
            a = _snap(seg.coords[0], verts, snap_tol)
            b = _snap(seg.coords[-1], verts, snap_tol)
            if math.hypot(b[0] - a[0], b[1] - a[1]) <= min_len:
                continue
            rows.append(_row(f"{label}-{len(rows) + 1:02d}", a, b, p))
    if not rows:
        raise ValueError(
            "The tessellation produced no joint line inside the region. A block "
            "size close to the region's own size leaves every cell wall on its "
            "boundary; use a smaller block_size.")
    return rows


# ---------------------------------------------------------------------------
# The set, and the names its rows carry
# ---------------------------------------------------------------------------
#
# A set is one description — a dip, a spacing, a region — that RESOLVES into many
# joint lines. What a file keeps is the lines: the rows of the joints sheet are
# the input, and the mesher, the solver, the report and the reader all work on
# them. The description itself is not written down, and a set is not a thing the
# file remembers (the owner's ruling: "why not just store the joints? if you want
# to recreate them, make another set").
#
# What a generated row carries beyond its geometry and its properties is its
# NAME: the set's name and the row's position in it — "bed-01", "bed-02", … So a
# network reads as one thing in a list, a message or a table, and the rows of one
# set can be found together, which is how a set is changed: remove it (Remove
# set, in the joints editor) and build another.

#: Every parameter a kind of set takes, with the value that means "not stated".
_SET_PARAMS = {
    "parallel": {"dip": None, "spacing": None, "offset": 0.0,
                 "trace_len": None, "gap": None},
    "cross": {"dip": None, "spacing": None, "offset": 0.0,
              "dip2": None, "spacing2": None, "offset2": 0.0,
              "trace_len": None, "gap": None},
    "voronoi": {"block_size": None, "seed": None},
}


def _canonical_region(region):
    """A set's region in one spelling.

    ``'Sandstone'``, ``['Sandstone', 'Shale']`` and ``'mat:Sandstone'`` all name
    a material, and the generators accept every one of them; a set keeps the
    explicit form, so two sets described the same way compare equal.
    """
    if region is None:
        return None
    if isinstance(region, str):
        text = region.strip()
        low = text.lower()
        if low.startswith("poly:"):
            return "poly:" + text[5:].strip()
        if low.startswith("mat:"):
            return "mat:" + "+".join(p.strip() for p in text[4:].split("+")
                                     if p.strip())
        return f"mat:{text}"
    if isinstance(region, (Polygon, MultiPolygon)):
        return region
    items = list(region)
    if not all(isinstance(r, str) for r in items):
        return Polygon(items)                # a coordinate list, given outright
    names = [r.strip() for r in items]
    if not all(names):
        raise ValueError("A joint set's region names an empty material.")
    return "mat:" + "+".join(names)


def set_name(label):
    """The set a row label belongs to, or ``''`` for a row that is in none.

    A set's rows are named ``bed-01``, ``bed-02``, …, so the set is whatever
    stands before the last ``-NN``. A label of any other shape — ``'base
    joint'``, ``'bed'``, ``''`` — names no set, which is what a hand-entered line
    carries and what keeps it out of anything done to a set as a whole.
    """
    name, _, index = str(label or "").strip().rpartition("-")
    return name.strip() if name.strip() and index.isdigit() else ""


class JointSet:
    """One joint set: its kind, its parameters, its region and the properties
    every one of its rows carries.

    What the Build network dialog fills in and :meth:`generate` resolves into the
    rows of the ``joints`` sheet. Nothing keeps it afterwards — the rows are the
    input — so a set is changed by removing its rows and building another.

    Attributes
    ----------
    name : str
        The set's name. Its rows are named ``<name>-01``, ``<name>-02``, …
    kind : {'parallel', 'cross', 'voronoi'}
    params : dict
        The kind's own parameters; :data:`_SET_PARAMS` lists which each takes.
    region : optional
        As :func:`resolve_region`.
    band : (y_min, y_max), optional
        An elevation band the region is cut to, either end ``None`` for open.
    props : dict, optional
        The joint properties every row of the set gets. They land in the rows'
        own columns, where one row's strength can still be edited by hand.
    """

    def __init__(self, name, kind, params=None, region=None, band=None,
                 props=None):
        kind = str(kind).strip().lower()
        if kind not in _SET_PARAMS:
            raise ValueError(
                f"{kind!r} is not a joint set kind. The kinds are: "
                + ", ".join(sorted(_SET_PARAMS)) + ".")
        name = str(name).strip()
        if not name:
            raise ValueError(
                "A joint set needs a name; its rows are named with it.")
        self.name = name
        self.kind = kind
        self.params = dict(_SET_PARAMS[kind])
        for k, v in dict(params or {}).items():
            if k not in self.params:
                raise ValueError(
                    f"{k!r} is not a parameter of a {kind} joint set. It takes: "
                    + ", ".join(sorted(_SET_PARAMS[kind])) + ".")
            self.params[k] = v
        self.region = _canonical_region(region)
        self.band = tuple(band) if band is not None else None
        self.props = dict(props or {})

    def __eq__(self, other):
        if not isinstance(other, JointSet):
            return NotImplemented
        return (self.name == other.name and self.kind == other.kind
                and self.params == other.params and self.region == other.region
                and self.band == other.band)

    def __repr__(self):                     # pragma: no cover - debugging aid
        return f"<JointSet {self.name!r} ({self.kind})>"

    def row_label(self, index):
        """The Label cell of row ``index`` (1-based) of this set."""
        return f"{self.name}-{int(index):02d}"

    # -- generating ---------------------------------------------------------
    def generate(self, slope_data):
        """The rows this set resolves to, each named for the set.

        The same list the generators return, so a generated set is inspected,
        edited, saved and reloaded like a typed one.
        """
        p = self.params
        persistence = (None if p.get("trace_len") is None
                       else (p["trace_len"], p.get("gap") or 0.0))
        if self.kind == "voronoi":
            if p.get("block_size") is None or p.get("seed") is None:
                raise ValueError(
                    "A Voronoi joint set needs a block size and a seed; the seed "
                    "is what makes the network reproducible.")
            rows = voronoi(slope_data, p["block_size"], int(p["seed"]),
                           region=self.region, label=self.name,
                           props=self.props, band=self.band)
        elif self.kind == "cross":
            for key in ("dip", "dip2", "spacing", "spacing2"):
                if p.get(key) is None:
                    raise ValueError(
                        "A cross-jointed set needs a dip and a spacing for both "
                        "of its sets.")
            # Generated under throwaway names so cross_jointed's own check — two
            # sets may MEET but may not lie on one another — reads them as a
            # pair; the set's own names go on below.
            first = parallel_set(slope_data, p["dip"], p["spacing"],
                                 offset=p.get("offset") or 0.0,
                                 persistence=persistence, region=self.region,
                                 label="a", props=self.props, band=self.band)
            second = parallel_set(slope_data, p["dip2"], p["spacing2"],
                                  offset=p.get("offset2") or 0.0,
                                  persistence=persistence, region=self.region,
                                  label="b", props=self.props, band=self.band)
            rows = cross_jointed(first, second)
        else:
            if p.get("dip") is None or p.get("spacing") is None:
                raise ValueError(
                    "A parallel joint set needs a dip and a spacing.")
            rows = parallel_set(slope_data, p["dip"], p["spacing"],
                                offset=p.get("offset") or 0.0,
                                persistence=persistence, region=self.region,
                                label=self.name, props=self.props,
                                band=self.band)
        for i, row in enumerate(rows):
            row["label"] = self.row_label(i + 1)
        return rows


def sets_in(slope_data):
    """Every named set in the model, in the order its rows first appear.

    Returns a list of ``(name, [row indices])`` — the grouping "Remove set" and
    anything else that acts on a whole network uses. A row whose label names no
    set (:func:`set_name`) belongs to none of them and appears in none.
    """
    rows = (slope_data or {}).get("joint_lines") or []
    order, found = [], {}
    for i, row in enumerate(rows):
        name = set_name(row.get("label"))
        if not name:
            continue
        if name not in found:
            found[name] = []
            order.append(name)
        found[name].append(i)
    return [(name, found[name]) for name in order]
