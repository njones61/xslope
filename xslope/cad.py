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

"""CAD (DXF) import/export for xslope geometry.

Converts between DXF files and xslope input templates (see plans/plan_io.md §3):

  - export_dxf(slope_data, path): write a complete model to a layered DXF —
    material-zone polygons on per-material layers, plus profile lines, search
    circles, reinforcement, distributed loads, and piezo lines on reserved layers.

  - import_dxf(dxf_path, template, out_path): read material-zone polygons from a
    DXF and write them to the `polygons` sheet of an input template. Import is
    polygons-only; reserved feature layers are ignored. Material identity comes
    from the layer name.

The reader is robust to real-world DXFs: LWPOLYLINE and heavyweight POLYLINE,
arc bulges (flattened), unclosed rings (auto-closed), and loose LINE segments
(stitched into rings). The poly_test/ fixtures exercise each case.
"""

import re

import ezdxf
from ezdxf import colors as ezcolors
from ezdxf import path as ezpath

from .fileio import write_cells_to_xlsx, cell_ref

# Standard AutoCAD Color Index palette (index -> RGB), for mapping plot colors to
# an ACI color. ACI (code 62) is honored by essentially every CAD viewer, whereas
# layer true-color (code 420) is not.
_ACI_RGB = [ezcolors.int2rgb(ezcolors.DXF_DEFAULT_COLORS[i]) for i in range(256)]


def _nearest_aci(rgb):
    """Nearest AutoCAD Color Index to an (r, g, b) 0-255 color. Near-black maps to
    ACI 250 (pure black in the palette) rather than ACI 7 — ACI 7 is *white* in the
    palette and many viewers render it gray, not black, on a light background."""
    r, g, b = rgb
    if max(r, g, b) < 40:
        return 250
    best, best_d = 7, None
    for aci in range(1, 256):
        if aci == 7:
            continue
        cr, cg, cb = _ACI_RGB[aci]
        d = (cr - r) ** 2 + (cg - g) ** 2 + (cb - b) ** 2
        if best_d is None or d < best_d:
            best, best_d = aci, d
    return best


# Reserved layer names used by export for non-polygon features. On import these
# are ignored (import is polygons-only); every other layer is a material zone.
RESERVED_LAYERS = {
    'FAILURE_SURFACE', 'SEARCH_CIRCLES', 'REINFORCEMENT', 'DLOADS', 'PIEZO',
}
RESERVED_PREFIXES = ('PROFILE_',)

_DEFAULT_VERSION = 'R2010'
# AutoCAD Color Index palette, cycled by material id for visual distinction.
_ACI_PALETTE = [40, 30, 8, 5, 3, 1, 6, 2, 4, 7]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _layer_name(material_name, mat_id):
    """Material name -> a valid, CAD-conventional DXF layer name.

    Uppercase, spaces/illegal chars -> underscores (e.g. 'Silty Clay' ->
    'SILTY_CLAY'); falls back to MAT_<id> (1-based) when unnamed. DXF disallows
    < > / \\ " : ; ? * | = ` and control chars.
    """
    name = (material_name or '').strip()
    if not name:
        return f'MAT_{mat_id + 1}'
    name = re.sub(r'[<>/\\":;?*|=`\s]+', '_', name.upper()).strip('_')
    return name or f'MAT_{mat_id + 1}'


def _is_reserved_layer(name):
    up = name.strip().upper()
    return up in RESERVED_LAYERS or any(up.startswith(p) for p in RESERVED_PREFIXES)


def _dedupe_closing_vertex(coords):
    """Drop a trailing vertex equal to the first (a closed ring closes implicitly)."""
    if len(coords) >= 2 and coords[0] == coords[-1]:
        return coords[:-1]
    return coords


def _dedupe_consecutive(coords):
    """Remove consecutive duplicate vertices (zero-length segments). Some CAD
    importers choke on degenerate polylines and drop them and their neighbours."""
    out = []
    for c in coords:
        p = (float(c[0]), float(c[1]))
        if out and abs(out[-1][0] - p[0]) < 1e-9 and abs(out[-1][1] - p[1]) < 1e-9:
            continue
        out.append(p)
    return out


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def _entity_coords(e, arc_sag):
    """Flatten a LWPOLYLINE or POLYLINE to a list of (x, y), tessellating any arc
    bulges to chord sag `arc_sag`. Uses ezdxf.path, which handles both entity
    types and arc segments (LWPolyline has no .flattening() of its own)."""
    try:
        p = ezpath.make_path(e)
        return [(v.x, v.y) for v in p.flattening(arc_sag)]
    except Exception:
        # Last-resort fallback (no arc tessellation).
        try:
            return [(pt[0], pt[1]) for pt in e.get_points('xy')]
        except (AttributeError, TypeError):
            return [(v.dxf.location.x, v.dxf.location.y) for v in e.vertices]


def _entity_is_closed(e):
    for attr in ('is_closed', 'closed'):
        if hasattr(e, attr):
            try:
                return bool(getattr(e, attr))
            except TypeError:
                pass
    return False


def _stitch_lines(segments, tol=1e-6):
    """Stitch loose LINE segments into closed rings by matching shared endpoints.

    `segments` is a list of ((x1, y1), (x2, y2)). Returns (rings, leftover) where
    rings is a list of coordinate lists (closed, without the duplicate end vertex)
    and leftover is the count of segments that could not be chained into a ring.
    """
    def key(p):
        return (round(p[0] / tol), round(p[1] / tol))

    used = [False] * len(segments)
    rings = []
    leftover = 0

    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        a, b = segments[start]
        chain = [a, b]
        seg_count = 1
        extended = True
        while extended and key(chain[0]) != key(chain[-1]):
            extended = False
            tail = key(chain[-1])
            for j in range(len(segments)):
                if used[j]:
                    continue
                p, q = segments[j]
                if key(p) == tail:
                    chain.append(q); used[j] = True; seg_count += 1; extended = True; break
                if key(q) == tail:
                    chain.append(p); used[j] = True; seg_count += 1; extended = True; break

        if len(chain) >= 4 and key(chain[0]) == key(chain[-1]):
            rings.append(_dedupe_closing_vertex(chain))
        else:
            leftover += seg_count  # segments consumed in a chain that didn't close

    return rings, leftover


def read_dxf_polygons(dxf_path, arc_sag=0.05):
    """Read closed zones from a DXF as {'coords', 'layer'} dicts, with warnings.

    Handles LWPOLYLINE and POLYLINE (arc bulges flattened, unclosed rings
    auto-closed), and stitches loose LINE segments (grouped by layer) into rings.

    Returns:
        (polygons, warnings): polygons is a list of {'coords': [(x, y), ...],
        'layer': str}; warnings is a list of human-readable strings.
    """
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    polygons = []
    warnings = []

    for e in msp.query('LWPOLYLINE POLYLINE'):
        coords = _dedupe_closing_vertex(_entity_coords(e, arc_sag))
        layer = e.dxf.layer
        if len(coords) < 3:
            warnings.append(f"skipped polyline on layer '{layer}' with < 3 vertices")
            continue
        if not _entity_is_closed(e):
            warnings.append(f"open polyline on layer '{layer}' auto-closed")
        polygons.append({'coords': coords, 'layer': layer})

    # Loose LINE entities, grouped by layer, stitched into rings.
    lines_by_layer = {}
    for e in msp.query('LINE'):
        s = (e.dxf.start.x, e.dxf.start.y)
        t = (e.dxf.end.x, e.dxf.end.y)
        lines_by_layer.setdefault(e.dxf.layer, []).append((s, t))
    for layer, segs in lines_by_layer.items():
        rings, leftover = _stitch_lines(segs)
        for ring in rings:
            if len(ring) >= 3:
                polygons.append({'coords': ring, 'layer': layer})
        if rings:
            warnings.append(f"stitched {len(rings)} ring(s) from loose LINE segments on layer '{layer}'")
        if leftover:
            warnings.append(f"{leftover} LINE segment(s) on layer '{layer}' could not be stitched into a closed ring")

    return polygons, warnings


def dxf_to_polygons(dxf_path, layers=None, arc_sag=0.05):
    """Read material-zone polygons from a DXF (the import primitive).

    Reserved feature layers (PROFILE_*, FAILURE_SURFACE, SEARCH_CIRCLES,
    REINFORCEMENT, DLOADS, PIEZO) are ignored — import is polygons-only. If
    `layers` is given (an iterable of layer names), only those layers are kept.

    Returns:
        (polygons, warnings): polygons is a list of {'coords', 'layer'}.
    """
    polygons, warnings = read_dxf_polygons(dxf_path, arc_sag=arc_sag)
    keep = set(layers) if layers is not None else None
    out = []
    for p in polygons:
        if _is_reserved_layer(p['layer']):
            continue
        if keep is not None and p['layer'] not in keep:
            continue
        out.append(p)
    return out, warnings


def summarize_dxf(dxf_path, arc_sag=0.05):
    """Print a validation summary of the material-zone polygons in a DXF and
    return the polygon list. Use before import to confirm the extraction."""
    from shapely.geometry import Polygon

    polygons, warnings = dxf_to_polygons(dxf_path, arc_sag=arc_sag)
    print(f"{'Poly ID':>7} │ {'Layer Name':<18} │ {'Vertices':>8} │ {'Area':>10}")
    print(f"{'─'*7}─┼─{'─'*18}─┼─{'─'*8}─┼─{'─'*10}")
    for i, p in enumerate(polygons, start=1):
        area = Polygon(p['coords']).area if len(p['coords']) >= 3 else 0.0
        print(f"{i:>7} │ {p['layer']:<18} │ {len(p['coords']):>8} │ {area:>10.2f}")
    layers = sorted({p['layer'] for p in polygons})
    print(f"\nUnique material layers ({len(layers)}): {', '.join(layers)}")
    for w in warnings:
        print(f"  warning: {w}")
    return polygons


# ---------------------------------------------------------------------------
# Import: DXF -> polygons sheet
# ---------------------------------------------------------------------------

def import_dxf(dxf_path, template, out_path, material_map=None, arc_sag=0.05,
               seed_materials=True):
    """Read material-zone polygons from a DXF and write them to a template's
    `polygons` sheet.

    Parameters
    ----------
    dxf_path : str
        Source DXF file.
    template : str
        Path to an xlsx input template to copy and populate.
    out_path : str
        Output xlsx path (template is copied here, then the polygon/mat sheets
        are populated in place).
    material_map : dict, optional
        {layer_name: mat_id} with 1-based mat_id. If omitted, each unique layer
        (in first-appearance order) is assigned 1, 2, 3, ...
    arc_sag : float
        Chord sag for flattening arc bulges.
    seed_materials : bool
        If True, write each layer's name into the `mat` sheet as the material
        name (placeholder properties left for the user to fill in).

    Returns
    -------
    dict: {'polygons', 'layer_to_mat', 'warnings', 'out_path'}.
    """
    import shutil

    polygons, warnings = dxf_to_polygons(dxf_path, arc_sag=arc_sag)
    if not polygons:
        raise ValueError(f"No material-zone polygons found in {dxf_path}")

    # Build layer -> 1-based mat_id mapping.
    unique_layers = []
    for p in polygons:
        if p['layer'] not in unique_layers:
            unique_layers.append(p['layer'])
    if material_map is None:
        layer_to_mat = {lyr: i + 1 for i, lyr in enumerate(unique_layers)}
    else:
        layer_to_mat = dict(material_map)
        missing = [lyr for lyr in unique_layers if lyr not in layer_to_mat]
        if missing:
            raise ValueError(f"material_map is missing layers: {missing}")

    shutil.copy(template, out_path)
    updates = {}

    # polygon sheet: block p -> x_col = 1 + 3*(p-1), y_col = x_col + 1.
    # Mat ID at (row 5, y_col); coordinates from row 8 down.
    poly_cells = {}
    for p_idx, poly in enumerate(polygons, start=1):
        x_col = 1 + 3 * (p_idx - 1)
        y_col = x_col + 1
        poly_cells[cell_ref(5, y_col)] = layer_to_mat[poly['layer']]
        coords = _dedupe_closing_vertex(poly['coords'])
        for i, (x, y) in enumerate(coords):
            poly_cells[cell_ref(8 + i, x_col)] = float(x)
            poly_cells[cell_ref(8 + i, y_col)] = float(y)
    # Clear leftover Mat IDs in unused template blocks (template ships 1..15).
    for p_idx in range(len(polygons) + 1, 16):
        y_col = 2 + 3 * (p_idx - 1)
        poly_cells[cell_ref(5, y_col)] = ''
    updates['polygon'] = poly_cells

    # mat sheet: seed material names (row 9 = mat 1; name in column 2).
    if seed_materials:
        mat_cells = {}
        for lyr, mat_id in sorted(layer_to_mat.items(), key=lambda kv: kv[1]):
            mat_cells[cell_ref(9 + (mat_id - 1), 2)] = lyr
        updates['mat'] = mat_cells

    write_cells_to_xlsx(out_path, updates)
    return {'polygons': polygons, 'layer_to_mat': layer_to_mat,
            'warnings': warnings, 'out_path': out_path}


# ---------------------------------------------------------------------------
# Export: model -> layered DXF
# ---------------------------------------------------------------------------

def _ensure_layer(doc, name, color):
    if name not in doc.layers:
        doc.layers.add(name=name, color=color)


# ---------------------------------------------------------------------------
# Generic matplotlib Axes -> DXF (backs the save_dxf option on plot functions)
# ---------------------------------------------------------------------------

def axes_to_dxf(ax, dxf_path, version=_DEFAULT_VERSION):
    """Write the geometry drawn on a matplotlib Axes to a layered DXF.

    Walks the Axes artists (Line2D, collections, patches) and emits each as DXF
    geometry, **layered by artist gid, then label, then a type-based default**.
    Plot functions tag their draw calls with gid='LAYER_NAME' so e.g. slices,
    stress bars, the failure surface, and material zones each land on their own
    layer. Axes chrome (spines, ticks, gridlines, legend, title) is not an artist
    in ax.lines/collections/patches and is therefore skipped.

    Curve geometry that matplotlib already computed (e.g. contour LineCollections
    from tricontour) is captured directly, so contour/flow-net plots export
    without re-deriving anything.
    """
    import numpy as np
    import matplotlib.colors as mcolors
    from matplotlib.collections import PathCollection
    from matplotlib.patches import Circle as MplCircle

    doc = ezdxf.new(version)
    msp = doc.modelspace()

    def finite_pts(xs, ys):
        """Finite (x, y) pairs with consecutive duplicates removed. Adjacent slices
        share boundaries, so polylines like the line of thrust, failure surface, and
        slice outlines repeat points — producing zero-length segments. Some CAD
        importers choke on such degenerate polylines, dropping not only that entity
        but neighboring ones (observed: the thrust line's duplicate vertices caused
        the entire max-depth line to vanish on import)."""
        out = []
        for x, y in zip(xs, ys):
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            p = (float(x), float(y))
            if out and abs(out[-1][0] - p[0]) < 1e-9 and abs(out[-1][1] - p[1]) < 1e-9:
                continue
            out.append(p)
        return out

    def add_open_path(pts, att):
        """Add an open polyline. A 2-point segment is written as a LINE entity
        rather than a 2-vertex LWPOLYLINE: some CAD viewers cull an axis-aligned
        LWPOLYLINE with a zero-height/zero-width bounding box (e.g. the perfectly
        horizontal max-depth line), but a LINE always renders."""
        if len(pts) == 2:
            msp.add_line(pts[0], pts[1], dxfattribs=att)
        elif len(pts) >= 2:
            msp.add_lwpolyline(pts, dxfattribs=att)

    def artist_rgb(artist):
        """The artist's display color as an (r, g, b) 0-255 tuple, or None. Prefers
        the fill color, then the edge, then the line color — so a blue-filled patch
        with a black outline (e.g. pore pressure) reads as blue, while an unfilled
        patch (e.g. the hatched stress bars) reads as its edge color."""
        def rgb_of(getter):
            f = getattr(artist, getter, None)
            if f is None:
                return None
            try:
                rgba = mcolors.to_rgba_array(f())
            except (ValueError, TypeError):
                return None
            if len(rgba) == 0 or rgba[0][3] <= 0:   # missing or fully transparent
                return None
            r, g, b, _ = rgba[0]
            return (int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))

        return rgb_of('get_facecolor') or rgb_of('get_edgecolor') or rgb_of('get_color')

    def attribs(artist, default):
        """Return dxfattribs (layer + per-entity color) for an artist. The layer is
        gid -> label -> default; the color is the artist's own color set on *each
        entity* (ACI for universal CAD support, plus exact true-color) so artists
        that share a layer but differ in color (e.g. the profile lines) each keep
        their own color."""
        name = artist.get_gid()
        if not name:
            lbl = artist.get_label()
            name = lbl if (lbl and not lbl.startswith('_')) else default
        name = _layer_name(name, 0)
        rgb = artist_rgb(artist)
        aci = _nearest_aci(rgb) if rgb is not None else None
        if name not in doc.layers:
            if aci is not None:
                lyr = doc.layers.add(name=name, color=aci)
                lyr.rgb = rgb            # exact color for viewers that honor it
            else:
                doc.layers.add(name=name)
        att = {'layer': name}
        if aci is not None:
            att['color'] = aci           # ACI for universal support
            att['true_color'] = ezcolors.rgb2int(rgb)   # exact RGB (e.g. true black)
        return att

    # Line2D artists: failure surface, thrust line, slices, profiles, piezo, etc.
    # A marker-only line (linestyle 'None', e.g. 'k.' node markers, 'bs'/'ro' BC
    # markers) is emitted as points, not a connecting polyline.
    for ln in ax.get_lines():
        pts = finite_pts(ln.get_xdata(), ln.get_ydata())
        if not pts:
            continue
        att = attribs(ln, 'LINES')
        ls = ln.get_linestyle()
        has_line = ls not in ('None', 'none', ' ', '', None) and (ln.get_linewidth() or 0) > 0
        if has_line and len(pts) >= 2:
            add_open_path(pts, att)
        else:
            for pt in pts:
                msp.add_point(pt, dxfattribs=att)

    # Collections: line collections (tricontour lines, mesh edges), scatter
    # offsets (points), and polygon collections (filled contours / element fills).
    for coll in ax.collections:
        att = attribs(coll, 'COLLECTION')

        # Open polylines: anything that exposes segments (LineCollection and the
        # line ContourSets from tricontour) — never force-closed.
        segs = None
        if hasattr(coll, 'get_segments'):
            try:
                segs = coll.get_segments()
            except Exception:
                segs = None
        if segs:
            for seg in segs:
                pts = finite_pts(seg[:, 0], seg[:, 1])
                if len(pts) >= 2:
                    add_open_path(pts, att)
            continue

        if isinstance(coll, PathCollection):
            # Scatter: the offsets are the data points (the single path is just the
            # marker shape).
            for xy in np.atleast_2d(coll.get_offsets()):
                if len(xy) == 2 and np.isfinite(xy[0]) and np.isfinite(xy[1]):
                    msp.add_point((float(xy[0]), float(xy[1])), dxfattribs=att)
            continue

        # PatchCollection / PolyCollection / contour fills: the paths are the
        # geometry (in data coords). Close a ring only if its vertices actually
        # close, so open contour lines stay open.
        for path in getattr(coll, 'get_paths', lambda: [])():
            v = path.vertices
            pts = finite_pts(v[:, 0], v[:, 1])
            if len(pts) >= 2:
                closed = (abs(pts[0][0] - pts[-1][0]) < 1e-9 and
                          abs(pts[0][1] - pts[-1][1]) < 1e-9)
                msp.add_lwpolyline(pts, close=closed, dxfattribs=att)

    # Patches: filled material zones / stress bars (ax.fill -> Polygon), circles.
    for p in ax.patches:
        att = attribs(p, 'PATCHES')
        if isinstance(p, MplCircle):
            c = p.center
            msp.add_circle((float(c[0]), float(c[1])), float(p.radius), dxfattribs=att)
            continue
        try:
            verts = p.get_path().transformed(p.get_patch_transform()).vertices
        except Exception:
            continue
        pts = finite_pts(verts[:, 0], verts[:, 1])
        if len(pts) >= 2:
            msp.add_lwpolyline(pts, close=True, dxfattribs=att)

    doc.saveas(dxf_path)
    return dxf_path


def export_dxf(slope_data, dxf_path, version=_DEFAULT_VERSION):
    """Write an xslope model to a layered DXF (see plans/plan_io.md §3.3).

    Material-zone polygons go on per-material layers (named after the material);
    profile lines, search circles, reinforcement, distributed loads, and piezo
    lines go on their reserved layers. Whatever is present in slope_data is
    written; absent features are skipped.
    """
    materials = slope_data.get('materials', [])
    doc = ezdxf.new(version)
    msp = doc.modelspace()

    def mat_layer(mat_id):
        name = materials[mat_id]['name'] if (mat_id is not None and 0 <= mat_id < len(materials)) else None
        return _layer_name(name, mat_id if mat_id is not None else 0)

    # Material-zone polygons (closed LWPOLYLINE, one layer per material).
    for poly in slope_data.get('polygons', []) or []:
        mat_id = poly['mat_id']
        lname = mat_layer(mat_id)
        _ensure_layer(doc, lname, _ACI_PALETTE[(mat_id or 0) % len(_ACI_PALETTE)])
        coords = _dedupe_consecutive(_dedupe_closing_vertex(list(poly['polygon'].exterior.coords)))
        msp.add_lwpolyline(coords, close=True, dxfattribs={'layer': lname})

    # Profile lines (open LWPOLYLINE) on PROFILE_<mat>.
    for line in slope_data.get('profile_lines', []) or []:
        mat_id = line.get('mat_id')
        lname = f"PROFILE_{mat_layer(mat_id)}"
        _ensure_layer(doc, lname, 7)
        msp.add_lwpolyline(_dedupe_consecutive(line['coords']), close=False, dxfattribs={'layer': lname})

    # Search circles on SEARCH_CIRCLES: write the clipped failure-surface arc (the
    # portion within the slope), as plot_inputs does — not the full circle — plus a
    # center marker. Falls back to a full CIRCLE if the arc can't be formed.
    circles = slope_data.get('circles') or []
    ground_surface = slope_data.get('ground_surface')
    if circles:
        from .slice import generate_failure_surface
        _ensure_layer(doc, 'SEARCH_CIRCLES', 1)
        tcrack_depth = slope_data.get('tcrack_depth', 0) or 0
        for c in circles:
            arc = None
            if ground_surface is not None:
                ok, res = generate_failure_surface(ground_surface, circular=True,
                                                   circle=c, tcrack_depth=tcrack_depth)
                if ok:
                    clipped = res[4]
                    arc = _dedupe_consecutive(clipped.coords)
            if arc and len(arc) >= 2:
                msp.add_lwpolyline(arc, dxfattribs={'layer': 'SEARCH_CIRCLES'})
            elif c.get('R'):
                msp.add_circle((c['Xo'], c['Yo']), c['R'], dxfattribs={'layer': 'SEARCH_CIRCLES'})
            msp.add_point((c['Xo'], c['Yo']), dxfattribs={'layer': 'SEARCH_CIRCLES'})

    # Reinforcement lines (LINE) on REINFORCEMENT.
    reinforce = slope_data.get('reinforce_lines') or []
    if reinforce:
        _ensure_layer(doc, 'REINFORCEMENT', 2)
        for line in reinforce:
            pts = [(p['X'], p['Y']) for p in line]
            for a, b in zip(pts[:-1], pts[1:]):
                msp.add_line(a, b, dxfattribs={'layer': 'REINFORCEMENT'})

    # Distributed loads (LWPOLYLINE) on DLOADS.
    for key in ('dloads', 'dloads2'):
        for line in slope_data.get(key) or []:
            pts = _dedupe_consecutive([(p['X'], p['Y']) for p in line])
            if len(pts) >= 2:
                _ensure_layer(doc, 'DLOADS', 3)
                msp.add_lwpolyline(pts, close=False, dxfattribs={'layer': 'DLOADS'})

    # Piezometric lines (LWPOLYLINE) on PIEZO.
    for key in ('piezo_line', 'piezo_line2'):
        pl = slope_data.get(key) or []
        pts = _dedupe_consecutive(pl)
        if len(pts) >= 2:
            _ensure_layer(doc, 'PIEZO', 5)
            msp.add_lwpolyline(pts, close=False, dxfattribs={'layer': 'PIEZO'})

    doc.saveas(dxf_path)
    return dxf_path
