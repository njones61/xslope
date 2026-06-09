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
from ezdxf import path as ezpath

from .fileio import write_cells_to_xlsx, cell_ref


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
        coords = _dedupe_closing_vertex(list(poly['polygon'].exterior.coords))
        msp.add_lwpolyline(coords, close=True, dxfattribs={'layer': lname})

    # Profile lines (open LWPOLYLINE) on PROFILE_<mat>.
    for line in slope_data.get('profile_lines', []) or []:
        mat_id = line.get('mat_id')
        lname = f"PROFILE_{mat_layer(mat_id)}"
        _ensure_layer(doc, lname, 7)
        msp.add_lwpolyline(list(line['coords']), close=False, dxfattribs={'layer': lname})

    # Search circles (CIRCLE + center POINT) on SEARCH_CIRCLES.
    circles = slope_data.get('circles') or []
    if circles:
        _ensure_layer(doc, 'SEARCH_CIRCLES', 1)
        for c in circles:
            if c.get('R'):
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
            pts = [(p['X'], p['Y']) for p in line]
            if len(pts) >= 2:
                _ensure_layer(doc, 'DLOADS', 3)
                msp.add_lwpolyline(pts, close=False, dxfattribs={'layer': 'DLOADS'})

    # Piezometric lines (LWPOLYLINE) on PIEZO.
    for key in ('piezo_line', 'piezo_line2'):
        pl = slope_data.get(key) or []
        pts = [(x, y) for x, y in pl]
        if len(pts) >= 2:
            _ensure_layer(doc, 'PIEZO', 5)
            msp.add_lwpolyline(pts, close=False, dxfattribs={'layer': 'PIEZO'})

    doc.saveas(dxf_path)
    return dxf_path
