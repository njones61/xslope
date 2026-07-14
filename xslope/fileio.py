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

import os
import pickle
import re
import shutil
import subprocess
import tempfile
import zipfile
from xml.sax.saxutils import escape as xml_escape

import numpy as np
import pandas as pd
from lxml import etree
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from .mesh import import_mesh_from_json, build_polygons

def build_ground_surface(profile_lines):
    """
    Constructs the topmost ground surface LineString from a set of profile lines.

    The function finds the highest elevation at each x-coordinate across all profile lines,
    which represents the true ground surface.

    Parameters:
        profile_lines (list of dict): A list of profile lines, each represented
            as a dict with 'coords' key containing a list of (x, y) coordinate tuples.

    Returns:
        shapely.geometry.LineString: A LineString of the top surface, or an empty LineString
        if fewer than two valid points are found.
    """
    
    if not profile_lines:
        return LineString([])
    
    # Extract coordinate lists from profile line dicts
    coord_lists = [line['coords'] for line in profile_lines]
    
    # Step 1: Gather all points from all profile lines
    all_points = []
    for line in coord_lists:
        all_points.extend(line)
    
    # Step 2: Group points by x-coordinate and find the highest y for each x
    x_groups = {}
    for x, y in all_points:
        if x not in x_groups:
            x_groups[x] = y
        else:
            x_groups[x] = max(x_groups[x], y)
    
    # Step 3: For each candidate point, check if any profile line is above it
    ground_surface_points = []
    for x, y in sorted(x_groups.items()):
        # Create a vertical line at this x-coordinate
        vertical_line = LineString([(x, y - 1000), (x, y + 1000)])
        
        # Check intersections with all profile lines
        is_topmost = True
        for profile_line in coord_lists:
            line = LineString(profile_line)
            if line.length == 0:
                continue
            
            # Find intersection with this profile line
            intersection = line.intersection(vertical_line)
            if not intersection.is_empty:
                # Get the y-coordinate of the intersection
                if hasattr(intersection, 'y'):
                    # Single point intersection
                    if intersection.y > y + 1e-6:  # Allow small numerical tolerance
                        is_topmost = False
                        break
                elif hasattr(intersection, 'geoms'):
                    # Multiple points or line intersection
                    for geom in intersection.geoms:
                        if hasattr(geom, 'y') and geom.y > y + 1e-6:
                            is_topmost = False
                            break
                    if not is_topmost:
                        break
        
        if is_topmost:
            ground_surface_points.append((x, y))
    
    # Ensure we have at least 2 points
    if len(ground_surface_points) < 2:
        return LineString([])

    return LineString(ground_surface_points)


def build_ground_surface_from_polygons(polygons):
    """
    Derive the ground surface and domain polygon from material-zone polygons.

    The domain polygon is the union of all input polygons (the total extent of the
    model). The ground surface is the upper boundary of that union, traced from the
    top-left corner to the top-right corner.

    Parameters:
        polygons (list of dict): each dict has a 'polygon' key (shapely Polygon).

    Returns:
        tuple(LineString, shapely.geometry.Polygon): the ground surface (empty
        LineString if it cannot be formed) and the domain polygon.
    """
    domain = unary_union([p['polygon'] for p in polygons])
    if domain.geom_type == 'MultiPolygon':
        # Disjoint zones: use the largest connected component for the surface.
        domain = max(domain.geoms, key=lambda g: g.area)

    ring = list(domain.exterior.coords)[:-1]  # drop the closing duplicate vertex
    n = len(ring)
    if n < 3:
        return LineString([]), domain

    # The leftmost and rightmost vertices (ties broken by highest elevation) split
    # the exterior ring into two arcs; the one with the higher mean elevation is
    # the top boundary.
    imin = min(range(n), key=lambda i: (ring[i][0], -ring[i][1]))
    imax = max(range(n), key=lambda i: (ring[i][0], ring[i][1]))

    def arc(a, b):
        path, i = [], a
        while True:
            path.append(ring[i])
            if i == b:
                break
            i = (i + 1) % n
        return path

    arc1, arc2 = arc(imin, imax), arc(imax, imin)
    mean_y = lambda a: sum(p[1] for p in a) / len(a)
    upper = arc1 if mean_y(arc1) >= mean_y(arc2) else arc2

    # Order the surface left-to-right.
    if upper[0][0] > upper[-1][0]:
        upper = upper[::-1]

    return LineString(upper), domain


def _parse_polygon_sheet(xls, materials):
    """
    Parse the optional 'polygon' sheet into a list of material-zone polygons.

    The sheet uses the same block layout as the 'profile' sheet: polygon p occupies
    columns (x_col, y_col) with x_col = 3*p, a "Polygon #N" header in row 4, the
    Mat ID in row 5 of the y column, and (x, y) vertices from row 8 down. mat_id is
    stored 0-based to match the profile_lines / materials indexing convention.

    Returns:
        list of dict: [{'polygon': shapely Polygon, 'mat_id': int}, ...], empty if
        the sheet is absent or contains no polygons.
    """
    if 'polygon' not in xls.sheet_names:
        return []

    df = xls.parse('polygon', header=None)
    header_row, mat_id_row, coords_start_row = 3, 4, 7  # Excel rows 4, 5, 8
    polygons = []

    col = 0
    while col + 1 < df.shape[1]:
        x_col, y_col = col, col + 1

        # Stop when the block header is empty (no further polygon blocks).
        try:
            header_val = str(df.iloc[header_row, x_col]).strip()
        except Exception:
            break
        if not header_val or header_val.lower() == 'nan':
            break

        # Read vertices (row 8 down) until the first fully empty row.
        coords = []
        row = coords_start_row
        while row < df.shape[0]:
            try:
                x_val, y_val = df.iloc[row, x_col], df.iloc[row, y_col]
            except Exception:
                break
            if pd.isna(x_val) and pd.isna(y_val):
                break
            if pd.notna(x_val) and pd.notna(y_val):
                coords.append((float(x_val), float(y_val)))
            row += 1

        # A header with no vertices is an unused block (template default) — skip.
        if coords:
            if len(coords) >= 2 and coords[0] == coords[-1]:
                coords = coords[:-1]  # drop explicit closing vertex (Shapely closes)
            col_letter = chr(65 + col) if col < 26 else f"col {col + 1}"
            if len(coords) < 3:
                raise ValueError(
                    f"Each polygon must have at least 3 vertices. Polygon in block "
                    f"starting at column {col_letter} has {len(coords)}.")

            mat_id_val = df.iloc[mat_id_row, y_col]
            try:
                mat_id = int(float(mat_id_val)) - 1  # 1-based sheet -> 0-based
            except (ValueError, TypeError):
                mat_id = None
            if mat_id is None or mat_id < 0 or mat_id >= len(materials):
                raise ValueError(
                    f"Polygon in block starting at column {col_letter} has an invalid "
                    f"Mat ID ({mat_id_val!r}); it must reference a material in the "
                    f"'mat' sheet (1..{len(materials)}).")

            poly = Polygon(coords)
            if not poly.is_valid:
                raise ValueError(
                    f"Polygon in block starting at column {col_letter} is not a valid "
                    f"polygon (self-intersecting or degenerate).")

            polygons.append({'polygon': poly, 'mat_id': mat_id})

        col += 3  # next block (A->D->G->...)

    return polygons


def _validate_polygons_no_overlap(polygons):
    """Material zones must tile the section without overlapping — adjacent zones
    share matching edges; a zone is never drawn on top of another. Overlapping
    zones mesh incorrectly (a high-conductivity zone bridges over a low-
    conductivity barrier and can inflate the seepage flowrate several-fold), so
    reject them at load time instead of silently producing wrong results.

    Touching at shared edges/vertices is fine (zero overlap area); only a positive
    intersection area is an error.
    """
    tol = 1e-6
    for i in range(len(polygons)):
        pi = polygons[i]['polygon']
        for j in range(i + 1, len(polygons)):
            area = pi.intersection(polygons[j]['polygon']).area
            if area > tol:
                mi, mj = polygons[i].get('mat_id'), polygons[j].get('mat_id')
                m1 = mi + 1 if mi is not None else '?'
                m2 = mj + 1 if mj is not None else '?'
                raise ValueError(
                    f"Material zones overlap: polygon #{i + 1} (Mat ID {m1}) and "
                    f"polygon #{j + 1} (Mat ID {m2}) overlap by {area:.4g} sq units. "
                    f"Material zones must tile the section without overlapping — adjacent "
                    f"zones share matching edges, and a zone that sits inside or cuts "
                    f"through another (a lens, or a dam core) must be carved out of its "
                    f"neighbor (e.g. a shell around a core is one concave polygon with a "
                    f"notch), not drawn on top of it.")


def reinforce_available_tension(d1, d2, t_max, lp1, lp2, tend1=0.0, tend2=0.0):
    """Available tensile force at a point along a reinforcement line — the
    capacity envelope shared by the LEM point list and the FEM element taper:

        T = min( Tmax,
                 Tend1 + Tmax*d1/Lp1,     (Lp1 = 0 -> end 1 fully anchored)
                 Tend2 + Tmax*d2/Lp2 )    (Lp2 = 0 -> end 2 fully anchored)

    d1/d2 are the distances from the point to end 1 / end 2. Tend* are end
    anchorage capacities (plate/connection/anchor); 0 reproduces the classical
    friction-only taper exactly. One implementation for both engines, so the
    two can never drift.
    """
    cap1 = t_max if lp1 <= 0 else min(t_max, tend1 + t_max * d1 / lp1)
    cap2 = t_max if lp2 <= 0 else min(t_max, tend2 + t_max * d2 / lp2)
    return max(0.0, min(cap1, cap2))


def _reinforce_line_points(x1, y1, x2, y2, Tmax, Tres, Lp1, Lp2, E, Area,
                           Tend1=0.0, Tend2=0.0):
    """Build the LEM tension-distribution point list for ONE reinforcement line
    from its raw endpoints, pullout lengths, and end anchorage. Returns [] for a
    zero-length line.

    The available tension is the piecewise-linear capacity envelope of
    :func:`reinforce_available_tension`; the point list holds its breakpoints so
    linear interpolation between points reproduces the envelope exactly.

    Extracted from load_slope_data so the same derivation can be reused (e.g. by
    the GUI) to rebuild the display format after editing the raw line data."""
    import math
    line_length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    if line_length == 0:
        return []

    # Unit vector from (x1,y1) to (x2,y2)
    dx = (x2 - x1) / line_length
    dy = (y2 - y1) / line_length

    def T_at(s):
        return reinforce_available_tension(s, line_length - s, Tmax, Lp1, Lp2,
                                           Tend1, Tend2)

    # Candidate breakpoints: the endpoints plus every kink of the envelope —
    # where each end's ramp reaches Tmax, and where the two ramps cross.
    cands = {0.0, line_length}
    if Lp1 > 0 and Tend1 < Tmax:
        cands.add(min(line_length, (Tmax - Tend1) * Lp1 / Tmax))
    if Lp2 > 0 and Tend2 < Tmax:
        cands.add(max(0.0, line_length - (Tmax - Tend2) * Lp2 / Tmax))
    if Lp1 > 0 and Lp2 > 0:
        m1, m2 = Tmax / Lp1, Tmax / Lp2
        s_x = (Tend2 + m2 * line_length - Tend1) / (m1 + m2)
        # only a breakpoint when the ramps cross BELOW the Tmax plateau
        if 0.0 < s_x < line_length and (Tend1 + m1 * s_x) < Tmax:
            cands.add(s_x)

    # Emit ordered, deduplicated breakpoints. The stored Tres mirrors the old
    # semantics: the material residual where the full capacity is available,
    # zero inside a friction ramp (sudden pullout) unless the governing end is
    # anchored, in which case the anchor hardware survives up to min(Tres, Tend).
    # An UNSET Tres (NaN) stays unset everywhere along the line — no post-peak
    # drop anywhere — rather than collapsing to a number through min().
    tol = max(1e-9, 1e-9 * line_length)
    line_points = []
    for s in sorted(cands):
        if line_points and abs(s - line_points[-1]["_s"]) < tol:
            continue
        T = T_at(s)
        if Tres != Tres:                      # NaN: unset, no post-peak drop
            tres_pt = float('nan')
        elif T >= Tmax - 1e-12:
            tres_pt = Tres
        else:
            governing_tend = tend_g = 0.0
            cap1 = Tmax if Lp1 <= 0 else Tend1 + Tmax * s / Lp1
            cap2 = Tmax if Lp2 <= 0 else Tend2 + Tmax * (line_length - s) / Lp2
            tend_g = Tend1 if cap1 <= cap2 else Tend2
            tres_pt = min(Tres, tend_g)
        line_points.append({"X": x1 + s * dx, "Y": y1 + s * dy, "T": T,
                            "Tres": tres_pt, "E": E, "Area": Area, "_s": s})
    for p in line_points:
        del p["_s"]

    return line_points


def build_reinforce_lines(reinforcement_lines):
    """Map raw FEM-format reinforcement dicts (x1,y1,x2,y2,t_max,t_res,lp1,lp2,E,
    area, and optionally tend1/tend2) to the LEM display/analysis format (a list
    of tension point-lists). Lines resolving to fewer than 2 points are dropped,
    matching load_slope_data."""
    lines = []
    for r in reinforcement_lines:
        pts = _reinforce_line_points(
            r["x1"], r["y1"], r["x2"], r["y2"], r["t_max"], r["t_res"],
            r["lp1"], r["lp2"], r["E"], r["area"],
            r.get("tend1", 0.0), r.get("tend2", 0.0))
        if len(pts) >= 2:
            lines.append(pts)
    return lines


# Highest input-template version this build can read. Bump together with the
# template (docs/inputs/input_template.xlsx, main!D5) and its reader support.
SUPPORTED_TEMPLATE_VERSION = 15


def _read_seep_bc_sheet(seep_df, sheet_name):
    """Parse one ``seep bc`` sheet into specified heads, fluxes, and the exit face.

    Layout (v15): exit face in B/C from row 5 down; then up to 5 BC blocks in
    E/F, H/I, K/L, N/O, Q/R. Row 3 of a block holds the TYPE cell (x column) and
    the VALUE (y column); coordinates start at row 5. A block ends the scan when
    its VALUE cell is empty.

    Backward compatibility: templates v14 and earlier put the literal label
    "Head:" in the type cell, so a block is a flux BC only when the type text
    starts with "flux" — every older file therefore still loads as a head BC.
    """
    bc = {"specified_heads": [], "specified_fluxes": [], "exit_face": []}

    def _read_coords(x_col, y_col, start_row=4):
        coords = []
        row = start_row
        while row < seep_df.shape[0]:
            try:
                x_val = seep_df.iloc[row, x_col]
                y_val = seep_df.iloc[row, y_col]
                if pd.isna(x_val):
                    break
                if pd.notna(x_val) and pd.notna(y_val):
                    coords.append((float(x_val), float(y_val)))
            except Exception:
                break
            row += 1
        return coords

    bc["exit_face"] = _read_coords(1, 2)

    type_row = 2  # Excel row 3 (0-indexed)
    col = 4       # column E
    while col + 1 < seep_df.shape[1]:
        if seep_df.shape[0] <= type_row:
            break
        value = seep_df.iloc[type_row, col + 1]
        if pd.isna(value):
            break  # empty VALUE cell ends the block scan

        type_cell = seep_df.iloc[type_row, col]
        is_flux = (pd.notna(type_cell)
                   and str(type_cell).strip().lower().startswith("flux"))

        coords = _read_coords(col, col + 1)
        if is_flux:
            block = (col - 4) // 3 + 1
            if len(coords) < 2:
                raise ValueError(
                    f"Flux BC #{block} on sheet '{sheet_name}' has "
                    f"{len(coords)} coordinate(s). A flux BC is applied over the "
                    "edges of a polyline and needs at least 2 points."
                )
            bc["specified_fluxes"].append({"flux": float(value), "coords": coords})
        elif coords:
            bc["specified_heads"].append({"head": float(value), "coords": coords})

        col += 3  # E -> H -> K -> ...

    return bc


def load_slope_data(filepath):
    """
    This function reads input data from various Excel sheets and parses it into
    structured components used throughout the slope stability analysis framework.
    It handles circular and non-circular failure surface data, reinforcement, piezometric
    lines, and distributed loads.

    Validation is enforced to ensure required geometry and material information is present:
    - Circular failure surface: must contain at least one valid row with Xo and Yo
    - Non-circular failure surface: required if no circular data is provided
    - Profile lines: must contain at least one valid set, and each line must have ≥ 2 points
    - Materials: must match the number of profile lines
    - Piezometric line: only included if it contains ≥ 2 valid rows
    - Distributed loads and reinforcement: each block must contain ≥ 2 valid entries

    Raises:
        ValueError: if required inputs are missing or inconsistent.

    Returns:
        dict: Parsed and validated global data structure for analysis
    """

    xls = pd.ExcelFile(filepath)
    globals_data = {}

    # === STATIC GLOBALS ===
    main_df = xls.parse('main', header=None)

    try:
        template_version = main_df.iloc[4, 3]  # Excel row 5, column D
        gamma_water = float(main_df.iloc[7, 3])  # Excel row 8, column D
        tcrack_depth = float(main_df.iloc[8, 3])  # Excel row 9, column D
        tcrack_water = float(main_df.iloc[9, 3])  # Excel row 10, column D
        k_seismic = float(main_df.iloc[10, 3])  # Excel row 11, column D
    except Exception as e:
        raise ValueError(f"Error reading static global values from 'main' tab: {e}")

    # Template version gate. Refuse files NEWER than this build understands, so a
    # newer template can never be silently mis-read by an older install (new
    # columns ignored, options like u='ru' silently zeroed). Shipped models carry
    # versions 8-13; all load through the header-name-driven readers below.
    try:
        _tv = int(float(template_version))
    except (TypeError, ValueError):
        raise ValueError(
            f"Unrecognized template version {template_version!r} in cell D5 of the "
            f"'main' sheet. Expected a number (current version is "
            f"{SUPPORTED_TEMPLATE_VERSION}).")
    if _tv > SUPPORTED_TEMPLATE_VERSION:
        raise ValueError(
            f"This input file is template version {_tv}, but this installation of "
            f"xslope supports versions up to {SUPPORTED_TEMPLATE_VERSION}. "
            "Update xslope to read this file.")

    # === PROFILE LINES ===
    profile_df = xls.parse('profile', header=None)

    max_depth = float(profile_df.iloc[1, 1])  # Excel B2 = row 1, column 1

    profile_lines = []
    
    # New format: single data block, profile lines arranged horizontally
    # First profile line: columns A:B, second: D:E, third: G:H, etc.
    # Header row is row 4 (index 3), mat_id is in B5 (row 4, column 1)
    # XY coordinates start in row 7 (index 6)
    header_row = 3  # Excel row 4 (0-indexed)
    mat_id_row = 4  # Excel row 5 (0-indexed)
    coords_start_row = 7  # Excel row 8 (0-indexed)
    
    col = 0  # Start with column A (index 0)
    while col < profile_df.shape[1]:
        x_col = col
        y_col = col + 1
        
        # Check if header row is empty (stop reading if empty)
        try:
            header_val = str(profile_df.iloc[header_row, x_col]).strip()
            if not header_val or header_val.lower() == 'nan':
                break  # No more profile lines
        except:
            break  # No more profile lines
        
        # Read mat_id from B5 (row 4, column 1) for this profile line
        # Convert from 1-based to 0-based for internal use
        try:
            mat_id_val = profile_df.iloc[mat_id_row, y_col]
            if pd.isna(mat_id_val):
                mat_id = None
            else:
                # Convert to integer and subtract 1 to make it 0-based.
                # Range validation happens later, once materials are parsed.
                mat_id = int(float(mat_id_val)) - 1
                if mat_id < 0:
                    mat_id = None  # Invalid mat_id
        except (ValueError, TypeError):
            mat_id = None
        
        # Read XY coordinates starting from row 7, stop at first empty row
        coords = []
        row = coords_start_row
        while row < profile_df.shape[0]:
            try:
                x_val = profile_df.iloc[row, x_col]
                y_val = profile_df.iloc[row, y_col]
                
                # Stop at first empty row (both x and y are empty)
                if pd.isna(x_val) and pd.isna(y_val):
                    break
                
                # If at least one coordinate is present, try to convert
                if pd.notna(x_val) and pd.notna(y_val):
                    coords.append((float(x_val), float(y_val)))
            except:
                break
            row += 1
        
        # Validate that we have at least 2 points
        if len(coords) == 1:
            raise ValueError(f"Each profile line must contain at least two points. Profile line starting at column {chr(65 + col)} has only one point.")
        
        if len(coords) >= 2:
            # Store as dict with coords and mat_id
            profile_lines.append({
                'coords': coords,
                'mat_id': mat_id
            })
        
        # Move to next profile line (skip 3 columns: A->D, D->G, etc.)
        col += 3

    # Ground surface, domain polygon, and tensile-crack line are built from the
    # unified polygon representation after materials are parsed — see below.

    # === MATERIALS (Optimized Parsing) ===
    # Locate the header row rather than assuming row 8 -- it has moved before (pre-v8
    # templates put it on row 3) and a legend row above the table moves it again.
    _mat_raw = xls.parse('mat', header=None)
    _mat_hdr = _find_mat_header_row(
        lambda r, c: (_mat_raw.iloc[r - 1, c - 1]
                      if r <= _mat_raw.shape[0] and c <= _mat_raw.shape[1] else None))
    mat_df = _mat_raw.iloc[_mat_hdr:].reset_index(drop=True)
    # Header names are matched underscore-insensitively: the master template's
    # styled headers concatenate to 'pow_a'/'vg_a' while older sheets carry
    # 'powa'/'vga' - both must resolve to the same column.
    mat_df.columns = [('' if pd.isna(v) else str(v).strip().replace('_', ''))
                      for v in _mat_raw.iloc[_mat_hdr - 1]]
    materials = []

    def _num(x):
        v = pd.to_numeric(x, errors="coerce")
        return float(v) if pd.notna(v) else 0.0

    def _pick(row, *names):
        """First of ``names`` that is actually a COLUMN on this sheet.

        v14 renamed the two unsaturated-curve parameters from the van-Genuchten-
        specific 'vga'/'vgn' to the law-agnostic 'a'/'n' (they now serve van
        Genuchten and Gardner alike), and renamed Poisson's ratio from 'n' to
        'nu' — which had collided with the new 'n'. Order matters: a v14 sheet
        has BOTH 'n' (the curve exponent) and 'nu', so nu must prefer 'nu' and
        only fall back to 'n' on a pre-v14 sheet that has no 'nu' column."""
        for nm in names:
            if nm in row.index:
                return row.get(nm, 0)
        return 0

    def _choice(x, default):
        """Normalize a free-text option cell. An empty cell reaches here as the
        float NaN, which ``str()`` renders as 'nan' -- treat both as unset."""
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return default
        s = str(x).strip().lower()
        return default if s in ('', 'nan') else s

    # Read materials row by row until we encounter an empty material name (Column B)
    # Data starts at Excel row 9 (0-indexed row 0 after header=7)
    for i in range(len(mat_df)):
        row = mat_df.iloc[i]
        
        # Check if material name (Column B) is empty - stop reading if empty
        material_name = row.get('name', '')
        if pd.isna(material_name) or str(material_name).strip() == '':
            break  # Stop reading when we encounter an empty material name
        
        # For seep workflows, 'g' (unit weight) and shear strength properties are not required.
        # A material row is considered "missing" only if EVERY property column after
        # 'name' is empty. (Expressed position-free: a column insert can never
        # silently narrow this check the way the old hardcoded C:X window could.)
        props_empty = row.iloc[2:].isna().all() if mat_df.shape[1] > 2 else True
        if props_empty:
            # Excel row number: first data row sits just below the located header
            excel_row = _mat_hdr + 1 + i
            raise ValueError(
                "CRITICAL ERROR: Material row has empty property fields. "
                f"Material '{material_name}' (Excel row {excel_row}) has no property "
                "values after the name column."
            )

        # Unsaturated relative-permeability model (template v11+): 'lf' (linear
        # front, kr0/h0 apply), 'vg' (van Genuchten) or 'gard' (Gardner, v14+).
        # The last two share the a/n parameter pair. Older templates lack the
        # column -> default 'lf' (current behavior). Read by header name, so
        # column position is version-independent.
        unsat_raw = row.get('unsat', 'lf')
        unsat_val = str(unsat_raw).strip().lower() if (
            pd.notna(unsat_raw) and str(unsat_raw).strip().lower() not in ('', 'nan')
        ) else 'lf'
        if unsat_val not in ('lf', 'vg', 'gard'):
            unsat_val = 'lf'

        # Excel row number: first data row sits just below the located header
        excel_row = _mat_hdr + 1 + i

        # Pore pressure option. An unrecognized value used to fall through to
        # u = 0 in slice.py, silently deleting pore pressure and inflating FS.
        u_val = _choice(row.get('u'), 'none')
        if u_val not in ('none', 'piezo', 'seep', 'ru'):
            raise ValueError(
                f"Material '{material_name}' (mat sheet, Excel row {excel_row}) has an "
                f"unrecognized pore pressure option u='{u_val}'. "
                "Expected one of: none, piezo, seep, ru (or leave blank for none)."
            )

        # Strength model. Blank is allowed -- seep-only material rows carry no
        # strength -- but slice.py raises if a blank one reaches a failure surface.
        option_val = _choice(row.get('option'), '')
        if option_val not in ('', 'mc', 'cp', 'pow', 'hb'):
            raise ValueError(
                f"Material '{material_name}' (mat sheet, Excel row {excel_row}) has an "
                f"unrecognized strength option option='{option_val}'. "
                "Expected one of: mc, cp, pow, hb."
            )

        # v12 columns, read by header name. Older templates lack them entirely:
        # row.get() returns None and the defaults preserve pre-v12 behavior exactly.
        gamma_val = _num(row.get("g", 0))
        _gsat_num = pd.to_numeric(row.get('gsat'), errors='coerce')
        gamma_sat_val = float(_gsat_num) if pd.notna(_gsat_num) else None
        if gamma_sat_val is not None and gamma_sat_val < gamma_val:
            raise ValueError(
                f"Material '{material_name}' (mat sheet, Excel row {excel_row}) has "
                f"gsat = {gamma_sat_val} < g = {gamma_val}. The saturated unit weight "
                "cannot be less than the moist unit weight; leave gsat blank to use "
                "g throughout.")

        ru_val = _num(row.get('ru', 0))
        if u_val == 'ru' and ru_val < 0:
            raise ValueError(
                f"Material '{material_name}' (mat sheet, Excel row {excel_row}) selects "
                f"u='ru' but has a negative pore pressure ratio ru = {ru_val}.")

        pow_a_val = _num(row.get('powa', 0))
        pow_b_val = _num(row.get('powb', 0))
        pow_c_val = _num(row.get('powc', 0))
        pow_d_val = _num(row.get('powd', 0))
        if option_val == 'pow' and (pow_a_val <= 0 or pow_b_val <= 0):
            raise ValueError(
                f"Material '{material_name}' (mat sheet, Excel row {excel_row}) selects "
                f"option='pow' but pow_a ({pow_a_val}) and pow_b ({pow_b_val}) must both "
                "be positive for the envelope tau = pow_a*(sigma_n + pow_d)^pow_b + pow_c.")

        # Generalized Hoek-Brown (v14). mb/s/a are DERIVED from GSI/mi/D at use
        # time (xslope.hoekbrown), so only the four field-observable inputs are
        # entered. GSI is defined on (0, 100] and D on [0, 1]; sigma_ci and mi
        # must be positive or the envelope collapses.
        hb_sci_val = _num(row.get('hbsci', 0))
        hb_gsi_val = _num(row.get('hbgsi', 0))
        hb_mi_val = _num(row.get('hbmi', 0))
        hb_d_val = _num(row.get('hbd', 0))
        if option_val == 'hb':
            if hb_sci_val <= 0 or hb_mi_val <= 0:
                raise ValueError(
                    f"Material '{material_name}' (mat sheet, Excel row {excel_row}) selects "
                    f"option='hb' but hb_sci ({hb_sci_val}) and hb_mi ({hb_mi_val}) must "
                    "both be positive (intact strength and the intact Hoek-Brown constant).")
            if not (0 < hb_gsi_val <= 100):
                raise ValueError(
                    f"Material '{material_name}' (mat sheet, Excel row {excel_row}) selects "
                    f"option='hb' but has hb_gsi = {hb_gsi_val}. The Geological Strength "
                    "Index must lie in (0, 100].")
            if not (0 <= hb_d_val <= 1):
                raise ValueError(
                    f"Material '{material_name}' (mat sheet, Excel row {excel_row}) selects "
                    f"option='hb' but has hb_d = {hb_d_val}. The disturbance factor must "
                    "lie in [0, 1] (0 = undisturbed, 1 = heavily blast-damaged).")

        materials.append({
            "name": str(material_name).strip(),
            "gamma": gamma_val,
            "gamma_sat": gamma_sat_val,
            "option": option_val,
            "c": _num(row.get('c', 0)),
            "phi": _num(row.get('f', 0)),
            "cp": _num(row.get('c/p', 0)),
            "r_elev": _num(row.get('r-elev', 0)),
            "d": _num(row.get('d', 0)) if pd.notna(row.get('d')) else 0,
            "psi": _num(row.get('psi', 0)) if pd.notna(row.get('psi')) else 0,
            "pow_a": pow_a_val,
            "pow_b": pow_b_val,
            "pow_c": pow_c_val,
            "pow_d": pow_d_val,
            "u": u_val,
            "ru": ru_val,
            "sigma_gamma": _num(row.get('s(g)', 0)),
            "sigma_c": _num(row.get('s(c)', 0)),
            "sigma_phi": _num(row.get('s(f)', 0)),
            "sigma_cp": _num(row.get('s(c/p)', 0)),
            "sigma_d": _num(row.get('s(d)', 0)),
            "sigma_psi": _num(row.get('s(psi)', 0)),
            "k1": _num(row.get('k1', 0)),
            "k2": _num(row.get('k2', 0)),
            "alpha": _num(row.get('alpha', 0)),
            "unsat": unsat_val,
            "kr0" : _num(row.get('kr0', 0)),
            "h0" : _num(row.get('h0', 0)),
            # v14: 'a'/'n' (law-agnostic); pre-v14: 'vga'/'vgn'. Same two slots
            # feed van Genuchten and Gardner, exactly as SEEP2D's uspar(1..2) do.
            # The OLD names are checked first because they are unambiguous: on a
            # pre-v14 sheet the bare 'n' column is POISSON'S RATIO, not the curve
            # exponent, so preferring 'n' there would silently read nu into vg_n.
            "vg_a": _num(_pick(row, 'vga', 'a')),
            "vg_n": _num(_pick(row, 'vgn', 'n')),
            "E": _num(row.get('E', 0)),
            # v14 renamed this to 'nu'; pre-v14 sheets call it 'n'
            "nu": _num(_pick(row, 'nu', 'n')),
            # Generalized Hoek-Brown (v14). mb/s/a are DERIVED from GSI/mi/D by
            # the Hoek-Brown 2002 relations, not entered — nobody has them to hand.
            "hb_sci": hb_sci_val,
            "hb_gsi": hb_gsi_val,
            "hb_mi": hb_mi_val,
            "hb_d": hb_d_val,
        })

    # === UNIFIED POLYGON REPRESENTATION ===
    # Geometry is always represented internally as material-zone polygons. If the
    # 'polygon' sheet is populated it is used directly; otherwise the profile lines
    # are converted to polygons via build_polygons(). All downstream code (slicing,
    # search) works from polygons.
    polygons_from_sheet = _parse_polygon_sheet(xls, materials)

    if polygons_from_sheet:
        if profile_lines:
            raise ValueError(
                "Both the 'profile' and 'polygon' sheets contain data. Use one "
                "geometry method, not both.")
        _validate_polygons_no_overlap(polygons_from_sheet)
        polygons = polygons_from_sheet
        # max_depth is a profile-sheet concept; it has no meaning for polygon input.
        max_depth = None
    elif profile_lines:
        # Convert profile lines -> polygons. max_depth is used ONLY here, as the
        # bottom boundary for build_polygons (mat_id is 0-based in both).
        polygons = [
            {'polygon': Polygon(p['coords']), 'mat_id': p['mat_id']}
            for p in build_polygons(slope_data={'profile_lines': profile_lines,
                                                'max_depth': max_depth})
        ]
    else:
        polygons = []

    # Derive the ground surface and domain polygon from the polygons. The domain
    # polygon (not max_depth) defines the bottom/side boundaries for all downstream
    # use (slice generation, search containment).
    domain_polygon = None
    if polygons:
        ground_surface, domain_polygon = build_ground_surface_from_polygons(polygons)
    else:
        ground_surface = LineString([])

    # === BUILD TENSILE CRACK LINE ===
    tcrack_surface = None
    if tcrack_depth > 0 and not ground_surface.is_empty:
        tcrack_surface = LineString(
            [(x, y - tcrack_depth) for (x, y) in ground_surface.coords])

    # === MESH AND SEEPAGE ANALYSIS FILES ===
    base, _ = os.path.splitext(filepath)
    mesh_filename = f"{base}_mesh.json"

    # Load mesh if it exists (used by both seep and fem workflows)
    mesh = None
    if os.path.exists(mesh_filename):
        try:
            mesh = import_mesh_from_json(mesh_filename)
        except Exception as e:
            print(f"WARNING: Error reading mesh file: {e}. Continuing without mesh.")

    # Load seepage solution files if any materials use seep pore pressure
    has_seep_materials = any(material["u"] == "seep" for material in materials)
    seep_u = None
    seep_u2 = None

    if has_seep_materials:
        try:
            solution1_filename = f"{base}_seep.csv"
            solution2_filename = f"{base}_seep2.csv"

            if mesh is not None and os.path.exists(solution1_filename):
                solution1_df = pd.read_csv(solution1_filename)
                solution1_df = solution1_df.iloc[:-1]
                seep_u = solution1_df["u"].to_numpy()

                if os.path.exists(solution2_filename):
                    solution2_df = pd.read_csv(solution2_filename)
                    solution2_df = solution2_df.iloc[:-1]
                    seep_u2 = solution2_df["u"].to_numpy()

        except Exception as e:
            print(f"WARNING: Error reading seepage files: {e}. Continuing without seep data.")

    # === PIEZOMETRIC LINE ===
    piezo_df = xls.parse('piezo', header=None)
    piezo_line = []
    piezo_line2 = []

    # Locate the x/y header row by content so both layouts load: v13 adds a
    # 'Type:' row (piezo | phreatic) above the headers, shifting data down one
    # row. The Type value sits one row above the header in the y column.
    _hdr = 2  # v12 default (headers in Excel row 3)
    for _r in range(min(8, piezo_df.shape[0])):
        if str(piezo_df.iloc[_r, 0]).strip().lower() == 'x':
            _hdr = _r
            break
    def _line_type(col):
        if _hdr >= 1:
            v = piezo_df.iloc[_hdr - 1, col]
            if pd.notna(v) and str(v).strip().lower() == 'phreatic':
                return True
            if pd.notna(v) and str(v).strip().lower() not in ('', 'piezo', 'type:'):
                raise ValueError(
                    f"Unrecognized piezometric line Type {v!r} on the 'piezo' "
                    f"sheet. Expected 'piezo' (static head, default) or "
                    f"'phreatic' (cos^2 inclination correction).")
        return False
    piezo_phreatic = _line_type(1)
    piezo_phreatic2 = _line_type(4)

    start_row = _hdr + 1
    x_col = 0  # Column A
    y_col = 1  # Column B
    
    row = start_row
    while row < piezo_df.shape[0]:
        try:
            x_val = piezo_df.iloc[row, x_col]
            y_val = piezo_df.iloc[row, y_col]
            
            # Stop at first empty row (both x and y are empty)
            if pd.isna(x_val) and pd.isna(y_val):
                break
            
            # If at least one coordinate is present, try to convert
            if pd.notna(x_val) and pd.notna(y_val):
                piezo_line.append((float(x_val), float(y_val)))
        except:
            break
        row += 1
    
    # Validate first piezometric line
    if len(piezo_line) == 1:
        raise ValueError("First piezometric line must contain at least two points.")
    
    # Read second piezometric line (columns D:E, starting at row 4, Excel row 4 = index 3)
    # Keep reading until we encounter an empty row
    x_col2 = 3  # Column D
    y_col2 = 4  # Column E
    
    row = start_row
    while row < piezo_df.shape[0]:
        try:
            x_val = piezo_df.iloc[row, x_col2]
            y_val = piezo_df.iloc[row, y_col2]
            
            # Stop at first empty row (both x and y are empty)
            if pd.isna(x_val) and pd.isna(y_val):
                break
            
            # If at least one coordinate is present, try to convert
            if pd.notna(x_val) and pd.notna(y_val):
                piezo_line2.append((float(x_val), float(y_val)))
        except:
            break
        row += 1
    
    # Validate second piezometric line (only if it has data)
    if len(piezo_line2) == 1:
        raise ValueError("Second piezometric line must contain at least two points if provided.")

    # === DISTRIBUTED LOADS ===
    # Read first set from "dloads" tab
    dload_df = xls.parse('dloads', header=None)
    dloads = []
    
    # Start reading from column B (index 1), each distributed load uses 3 columns (X, Y, Normal)
    # Keep reading to the right until we encounter an empty distributed load
    start_row = 3  # Excel row 4 (0-indexed row 3)
    col = 1  # Start with column B (index 1)
    
    while col < dload_df.shape[1]:
        x_col = col
        y_col = col + 1
        normal_col = col + 2
        
        # Check if dataframe has enough rows before accessing start_row
        if dload_df.shape[0] <= start_row:
            break  # Not enough rows, stop reading

        # Check if this distributed load block is empty (check first row for X coordinate)
        if pd.isna(dload_df.iloc[start_row, x_col]):
            break  # Stop reading when we encounter an empty distributed load
        
        # Read points for this distributed load, keep reading down until empty row
        block_points = []
        row = start_row
        while row < dload_df.shape[0]:
            try:
                x_val = dload_df.iloc[row, x_col]
                y_val = dload_df.iloc[row, y_col]
                normal_val = dload_df.iloc[row, normal_col]
                
                # Stop at first empty row (all three values are empty)
                if pd.isna(x_val) and pd.isna(y_val) and pd.isna(normal_val):
                    break
                
                # If at least coordinates are present, try to convert
                if pd.notna(x_val) and pd.notna(y_val):
                    normal = float(normal_val) if pd.notna(normal_val) else 0.0
                    block_points.append({
                        "X": float(x_val),
                        "Y": float(y_val),
                        "Normal": normal
                    })
            except:
                break
            row += 1
        
        # Validate that we have at least 2 points
        if len(block_points) == 1:
            raise ValueError(f"Each distributed load must contain at least two points. Distributed load starting at column {chr(65 + col)} has only one point.")
        
        if len(block_points) >= 2:
            dloads.append(block_points)
        
        # Move to next distributed load (skip 4 columns: 3 for the dload + 1 empty column)
        col += 4
    
    # Read second set from "dloads (2)" tab
    dloads2 = []
    try:
        dload_df2 = xls.parse('dloads (2)', header=None)
        
        # Start reading from column B (index 1), each distributed load uses 3 columns (X, Y, Normal)
        # Keep reading to the right until we encounter an empty distributed load
        col = 1  # Start with column B (index 1)
        
        while col < dload_df2.shape[1]:
            x_col = col
            y_col = col + 1
            normal_col = col + 2
            
            # Check if dataframe has enough rows before accessing start_row
            if dload_df2.shape[0] <= start_row:
                break  # Not enough rows, stop reading
            
            # Check if this distributed load block is empty (check first row for X coordinate)
            if pd.isna(dload_df2.iloc[start_row, x_col]):
                break  # Stop reading when we encounter an empty distributed load
            
            # Read points for this distributed load, keep reading down until empty row
            block_points = []
            row = start_row
            while row < dload_df2.shape[0]:
                try:
                    x_val = dload_df2.iloc[row, x_col]
                    y_val = dload_df2.iloc[row, y_col]
                    normal_val = dload_df2.iloc[row, normal_col]
                    
                    # Stop at first empty row (all three values are empty)
                    if pd.isna(x_val) and pd.isna(y_val) and pd.isna(normal_val):
                        break
                    
                    # If at least coordinates are present, try to convert
                    if pd.notna(x_val) and pd.notna(y_val):
                        normal = float(normal_val) if pd.notna(normal_val) else 0.0
                        block_points.append({
                            "X": float(x_val),
                            "Y": float(y_val),
                            "Normal": normal
                        })
                except:
                    break
                row += 1
            
            # Validate that we have at least 2 points
            if len(block_points) == 1:
                raise ValueError(f"Each distributed load must contain at least two points. Distributed load starting at column {chr(65 + col)} has only one point.")
            
            if len(block_points) >= 2:
                dloads2.append(block_points)
            
            # Move to next distributed load (skip 4 columns: 3 for the dload + 1 empty column)
            col += 4
    except (ValueError, KeyError):
        # If "dloads (2)" tab doesn't exist, just leave dloads2 as empty list
        pass

    # === CIRCLES ===

    # Read the first 3 rows to get the max depth
    raw_df = xls.parse('circles', header=None)  # No header, get full sheet

    # Read the circles data starting from row 2 (index 1)
    circles_df = xls.parse('circles', header=1)
    raw = circles_df.dropna(subset=['Xo', 'Yo'], how='any')
    circles = []
    for _, row in raw.iterrows():
        Xo = row['Xo']
        Yo = row['Yo']
        Option = row.get('Option', None)
        Depth = row.get('Depth', None)
        Xi = row.get('Xi', None)
        Yi = row.get('Yi', None)
        R = row.get('R', None)
        # For each circle, fill in the radius and depth values depending on the circle option
        if Option == 'Depth':
            R = Yo - Depth
        elif Option == 'Intercept':
            R = ((Xi - Xo) ** 2 + (Yi - Yo) ** 2) ** 0.5
            Depth = Yo - R
        elif Option == 'Radius':
            Depth = Yo - R
        else:
            raise ValueError(f"Unknown option '{Option}' for circles.")
        circle = {
            "Xo": Xo,
            "Yo": Yo,
            "Depth": Depth,
            "R": R,
        }
        circles.append(circle)

    # === NON-CIRCULAR SURFACES ===
    noncirc_df = xls.parse('non-circ')
    non_circ = list(noncirc_df.iloc[1:].dropna(subset=['Unnamed: 0']).apply(
        lambda row: {
            "X": float(row['Unnamed: 0']),
            "Y": float(row['Unnamed: 1']),
            "Movement": row['Unnamed: 2']
        }, axis=1))

    # === REINFORCEMENT LINES ===
    # Header-name-driven, so both template layouts load: v11 is
    # (# | x1 y1 x2 y2 Tmax Tres Lp1 Lp2 E Area); v12 inserts Label at column B,
    # regroups by analysis type, and adds Type/Dir/Appl/Tend1/Tend2/Spacing.
    reinforce_df = xls.parse('reinforce', header=1)  # Header in row 2 (0-indexed row 1)
    reinforce_df.columns = [str(c).strip().lower() for c in reinforce_df.columns]
    reinforcement_lines = []    # FEM format: list of dicts with raw line endpoints and properties

    # Type presets fill Dir/Appl when those cells are blank (mirrors the in-sheet
    # default formulas); explicit values win. Same table as the template legend.
    _TYPE_PRESETS = {
        'geosynthetic': ('tangent', 'active'),
        'nail':         ('axial',   'passive'),
        'tieback':      ('axial',   'active'),
        'anchor':       ('axial',   'active'),
    }

    for i, row in reinforce_df.iterrows():
        excel_row = i + 3
        # Stop reading when x1 is empty (NOT column B positionally -- in v12,
        # column B is the optional Label and may be blank on a data row)
        if pd.isna(row.get('x1')):
            break

        label = (str(row['label']).strip()
                 if 'label' in reinforce_df.columns and pd.notna(row.get('label'))
                 else f"Line {i + 1}")

        # Check if other required coordinates are present
        if pd.isna(row.get('y1')) or pd.isna(row.get('x2')) or pd.isna(row.get('y2')):
            continue  # Skip rows with incomplete coordinate data

        # If coordinates are present, check for required parameters (Tmax, Lp1, Lp2)
        if pd.isna(row.get('tmax')) or pd.isna(row.get('lp1')) or pd.isna(row.get('lp2')):
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                "coordinates but missing required parameters (Tmax, Lp1, Lp2). "
                "All three must be specified.")

        # v12 support-type columns; defaults reproduce pre-v12 behavior exactly
        # (generic tensile line: tangent direction, active application).
        rtype = _choice(row.get('type'), '')
        if rtype not in ('',) + tuple(_TYPE_PRESETS):
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                f"an unrecognized Type='{rtype}'. Expected one of: "
                f"{', '.join(_TYPE_PRESETS)} (or leave blank for a generic line).")
        _dir_def, _appl_def = _TYPE_PRESETS.get(rtype, ('tangent', 'active'))

        direction = _choice(row.get('dir'), _dir_def)
        if direction not in ('tangent', 'axial'):
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                f"an unrecognized Dir='{direction}'. Expected: tangent or axial.")

        appl = _choice(row.get('appl'), _appl_def)
        if appl not in ('active', 'passive'):
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                f"an unrecognized Appl='{appl}'. Expected: active or passive.")

        _sp_num = pd.to_numeric(row.get('spacing'), errors='coerce')
        spacing = float(_sp_num) if pd.notna(_sp_num) else 1.0
        if spacing <= 0:
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                f"Spacing = {spacing}; it must be positive (blank or 1 for geosynthetics).")

        tend1 = _num(row.get('tend1', 0))
        tend2 = _num(row.get('tend2', 0))
        if tend1 < 0 or tend2 < 0:
            raise ValueError(
                f"Reinforcement line '{label}' (reinforce sheet, Excel row {excel_row}) has "
                f"a negative end anchorage capacity (Tend1 = {tend1}, Tend2 = {tend2}).")

        try:
            # Extract coordinates and parameters into the raw (FEM) format. All
            # capacity terms and the axial stiffness are per unit width: discrete
            # supports enter per-element values plus Spacing and are divided here,
            # once, for both engines. (spacing defaults to 1 -> no-op for v11 files
            # and geosynthetics.)
            reinforcement_lines.append({
                "x1": float(row['x1']), "y1": float(row['y1']),
                "x2": float(row['x2']), "y2": float(row['y2']),
                "t_max": float(row['tmax']) / spacing,
                # A BLANK Tres means "no post-peak drop" — the bar is elastic-
                # perfectly-plastic and holds its capacity once it yields. It does
                # NOT mean zero. Zero is a legitimate, and very aggressive, entry:
                # it says the bar ruptures brittly and carries nothing afterwards.
                # Defaulting a blank cell to 0.0 silently made brittle rupture the
                # behaviour of every file that never mentions Tres. NaN carries the
                # "unset" sense through to the FEM, which softens only where t_res
                # is finite.
                "t_res": (float(row['tres']) / spacing
                          if pd.notna(row.get('tres')) else float('nan')),
                "lp1": float(row['lp1']) if not pd.isna(row['lp1']) else 0.0,
                "lp2": float(row['lp2']) if not pd.isna(row['lp2']) else 0.0,
                "E": float(row['e']) if pd.notna(row.get('e')) else float('nan'),
                "area": (float(row['area']) if pd.notna(row.get('area')) else float('nan')) / spacing,
                "label": label,
                "type": rtype,
                "dir": direction,
                "appl": appl,
                "tend1": tend1 / spacing,
                "tend2": tend2 / spacing,
                "spacing": spacing,
            })
        except Exception as e:
            raise ValueError(f"Error processing reinforcement line '{label}' in row {excel_row}: {e}")

    # LEM tension-distribution format, derived from the raw endpoints/pullout data.
    reinforce_lines = build_reinforce_lines(reinforcement_lines)


    # === PILE LINES ===
    pile_lines = []
    if 'piles' in xls.sheet_names:
        # Header-name-driven with case-normalized headers, so the v12 column
        # regrouping (and the v11 layout) both load. The force-angle column is
        # 'qp' in the template (theta_p rendered via Symbol font) but 'theta' in
        # some older files -- accept either.
        piles_df = xls.parse('piles', header=1)
        piles_df.columns = [str(c).strip().lower() for c in piles_df.columns]
        _theta_col = 'qp' if 'qp' in piles_df.columns else 'theta'
        for i, row in piles_df.iterrows():
            # Stop reading when column x1 is empty
            if pd.isna(row.get('x1')):
                break
            # Check required coordinates
            if pd.isna(row.get('y1')) or pd.isna(row.get('x2')) or pd.isna(row.get('y2')):
                continue
            try:
                x1, y1 = float(row['x1']), float(row['y1'])
                x2, y2 = float(row['x2']), float(row['y2'])
                H = float(row['h']) if pd.notna(row.get('h')) else None
                if pd.notna(row.get(_theta_col)):
                    theta_p = float(row[_theta_col])
                else:
                    # Auto-compute: perpendicular to pile axis (0 for vertical)
                    dx = x2 - x1
                    dy = y2 - y1
                    theta_p = np.degrees(np.arctan2(dx, -dy))
                D_pile = float(row['d']) if pd.notna(row.get('d')) else None
                S = float(row['s']) if pd.notna(row.get('s')) else None
                E_pile = float(row['e']) if pd.notna(row.get('e')) else None
                I_pile = float(row['i']) if pd.notna(row.get('i')) else None
                area = float(row['area']) if pd.notna(row.get('area')) else None
                V_cap = float(row['vcap']) if pd.notna(row.get('vcap')) else None
                M_cap = float(row['mcap']) if pd.notna(row.get('mcap')) else None
                fixity_raw = str(row['fixity']).strip().lower() if pd.notna(row.get('fixity')) else 'free'
                if fixity_raw not in ('free', 'fixed'):
                    raise ValueError(f"Fixity must be 'free' or 'fixed', got '{fixity_raw}'")
                fixity = fixity_raw
                # Force application (v12, LEM only): Active = allowable force, not
                # divided by FS (default, pre-v12 behavior); Passive = ultimate
                # capacity divided by FS.
                appl_raw = str(row['appl']).strip().lower() if pd.notna(row.get('appl')) else 'active'
                if appl_raw not in ('active', 'passive'):
                    raise ValueError(f"Appl must be 'active' or 'passive', got '{appl_raw}'")
                label = str(row['label']) if pd.notna(row.get('label')) else f"Pile {i+1}"

                # Validate
                line_length = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
                if line_length == 0:
                    continue
                if H is not None and H <= 0:
                    raise ValueError(f"H must be positive, got {H}")
                if V_cap is not None and V_cap <= 0:
                    raise ValueError(f"Vcap must be positive, got {V_cap}")
                if M_cap is not None and M_cap <= 0:
                    raise ValueError(f"Mcap must be positive, got {M_cap}")
                if (V_cap is not None or M_cap is not None) and S is None:
                    raise ValueError(f"S (pile spacing) is required when Vcap or Mcap are specified")

                pile_lines.append({
                    "x1": x1, "y1": y1,
                    "x2": x2, "y2": y2,
                    "H": H,
                    "theta_p": theta_p,
                    "D_pile": D_pile,
                    "S": S,
                    "E": E_pile,
                    "I": I_pile,
                    "area": area,
                    "V_cap": V_cap,
                    "M_cap": M_cap,
                    "fixity": fixity,
                    "appl": appl_raw,
                    "label": label,
                })
            except Exception as e:
                raise ValueError(f"Error processing pile in row {i + 3}: {e}")

    # === LINE LOADS (v12 'lloads' sheet) ===
    # A concentrated force per unit width applied at a point on the ground
    # surface (e.g. the weight of a shotcrete facing plate). Absent in pre-v12
    # templates -> empty list.
    line_loads = []
    if 'lloads' in xls.sheet_names:
        lloads_df = xls.parse('lloads', header=1)
        lloads_df.columns = [str(c).strip().lower() for c in lloads_df.columns]
        # Snap tolerance for "on the ground surface": 0.5% of the model height,
        # floored to a small absolute value for degenerate geometries.
        if not ground_surface.is_empty:
            _gs_ys = [p[1] for p in ground_surface.coords]
            _ll_tol = max(1e-6, 0.005 * (max(_gs_ys) - min(_gs_ys)))
        for i, row in lloads_df.iterrows():
            excel_row = i + 3
            if pd.isna(row.get('x')):
                break
            ll_label = (str(row['label']).strip()
                        if 'label' in lloads_df.columns and pd.notna(row.get('label'))
                        else f"Load {i + 1}")
            if pd.isna(row.get('y')) or pd.isna(row.get('p')):
                raise ValueError(
                    f"Line load '{ll_label}' (lloads sheet, Excel row {excel_row}) needs "
                    "x, y, and P; one or more are blank.")
            ll_x, ll_y = float(row['x']), float(row['y'])
            ll_p = float(row['p'])
            if ll_p <= 0:
                raise ValueError(
                    f"Line load '{ll_label}' (lloads sheet, Excel row {excel_row}) has "
                    f"P = {ll_p}; the magnitude must be positive (use Angle for direction).")
            ll_angle = float(row['angle']) if pd.notna(row.get('angle')) else -90.0
            # The load must act on the ground surface: snap small mismatches from
            # rounded coordinates, refuse anything farther than the tolerance.
            if not ground_surface.is_empty:
                _pt = Point(ll_x, ll_y)
                _d = ground_surface.distance(_pt)
                if _d > _ll_tol:
                    raise ValueError(
                        f"Line load '{ll_label}' (lloads sheet, Excel row {excel_row}) at "
                        f"({ll_x}, {ll_y}) is {_d:.3g} away from the ground surface "
                        f"(tolerance {_ll_tol:.3g}). Line loads must act on the ground surface.")
                _snapped = ground_surface.interpolate(ground_surface.project(_pt))
                ll_x, ll_y = float(_snapped.x), float(_snapped.y)
            line_loads.append({
                "x": ll_x, "y": ll_y,
                "P": ll_p,
                "angle": ll_angle,
                "label": ll_label,
            })


    # === SEEPAGE ANALYSIS BOUNDARY CONDITIONS ===
    seepage_bc = _read_seep_bc_sheet(xls.parse('seep bc', header=None), 'seep bc')
    try:
        seepage_bc2 = _read_seep_bc_sheet(xls.parse('seep bc (2)', header=None), 'seep bc (2)')
    except (ValueError, KeyError) as e:
        if isinstance(e, ValueError) and 'Flux BC' in str(e):
            raise
        # Sheet absent (older workbook) -> no second BC set.
        seepage_bc2 = {"specified_heads": [], "specified_fluxes": [], "exit_face": []}

    # === VALIDATION ===
 
    circular = len(circles) > 0
    # Check if this is a seep-only analysis (has seep BCs but no slope stability surfaces)
    has_seepage_bc = (len(seepage_bc.get("specified_heads", [])) > 0 or
                     len(seepage_bc.get("specified_fluxes", [])) > 0 or
                     len(seepage_bc.get("exit_face", [])) > 0)
    is_seepage_only = has_seepage_bc and not circular and len(non_circ) == 0
    # A mesh-based run with no LEM surfaces is a seepage or FEM (SSRM) analysis;
    # neither needs circular/non-circular failure surfaces.
    is_mesh_analysis = mesh is not None and not circular and len(non_circ) == 0

    # Only require circular/non-circular data for a pure LEM run (no seep BCs, no mesh)
    if not is_seepage_only and not is_mesh_analysis and not circular and len(non_circ) == 0:
        raise ValueError("Input must include either circular or non-circular surface data.")
    if not polygons:
        raise ValueError("Geometry is missing: provide either the 'profile' sheet or the 'polygon' sheet.")
    if not materials:
        raise ValueError("Materials sheet is empty.")

    # Every polygon must reference a material that exists in the 'mat' sheet. The
    # polygon-sheet path validates this at parse time; profile-derived polygons are
    # validated here (materials are not yet parsed when profile lines are read).
    for poly in polygons:
        mid = poly.get('mat_id')
        if mid is not None and (mid < 0 or mid >= len(materials)):
            raise ValueError(
                f"A geometry zone references an invalid Mat ID ({mid + 1}); it must "
                f"reference a material in the 'mat' sheet (1..{len(materials)}).")

    # For slope-stability (non seep-only) runs, every material referenced by the
    # geometry must have a positive unit weight. A gamma of 0 silently produces
    # zero slice weights and meaningless factors of safety; seep-only runs do not
    # use unit weight, so they are exempt.
    if not is_seepage_only:
        for poly in polygons:
            mid = poly.get('mat_id')
            if mid is not None and 0 <= mid < len(materials):
                if materials[mid].get('gamma', 0) <= 0:
                    raise ValueError(
                        f"Material '{materials[mid]['name']}' (Mat ID {mid + 1}) has a "
                        f"non-positive unit weight (gamma = {materials[mid]['gamma']}). "
                        f"A positive unit weight is required for slope-stability analysis.")


    # Add everything to globals_data
    globals_data["template_version"] = template_version
    globals_data["gamma_water"] = gamma_water
    globals_data["tcrack_depth"] = tcrack_depth
    globals_data["tcrack_water"] = tcrack_water
    globals_data["k_seismic"] = k_seismic
    globals_data["max_depth"] = max_depth
    globals_data["profile_lines"] = profile_lines
    globals_data["polygons"] = polygons
    globals_data["domain_polygon"] = domain_polygon
    globals_data["ground_surface"] = ground_surface
    globals_data["tcrack_surface"] = tcrack_surface
    globals_data["materials"] = materials
    globals_data["piezo_line"] = piezo_line
    globals_data["piezo_phreatic"] = piezo_phreatic
    globals_data["piezo_phreatic2"] = piezo_phreatic2
    globals_data["piezo_line2"] = piezo_line2
    globals_data["circular"] = circular # True if circles are present
    globals_data["circles"] = circles
    globals_data["non_circ"] = non_circ
    globals_data["dloads"] = dloads
    globals_data["dloads2"] = dloads2
    globals_data["reinforce_lines"] = reinforce_lines
    globals_data["reinforcement_lines"] = reinforcement_lines
    globals_data["pile_lines"] = pile_lines
    globals_data["line_loads"] = line_loads
    globals_data["seepage_bc"] = seepage_bc
    globals_data["seepage_bc2"] = seepage_bc2
    globals_data["has_seepage_bc2"] = bool(seepage_bc2.get("specified_heads")
                                           or seepage_bc2.get("specified_fluxes")
                                           or seepage_bc2.get("exit_face"))
    
    # Add mesh if available (used by both seep and fem workflows)
    globals_data["mesh"] = mesh

    # Add seep solution data if available
    if has_seep_materials:
        globals_data["seep_u"] = seep_u
        if seep_u2 is not None:
            globals_data["seep_u2"] = seep_u2

    return globals_data


def default_template_path():
    """Filesystem path to the blank XSLOPE input template bundled with the
    package (``xslope/resources/input_template.xlsx``).

    Pass it as the ``template=`` argument to :func:`save_slope_data_to_xlsx` to
    create a new file from the standard template — useful for library-only
    (no-GUI) users, since the template ships inside the installed package and so
    is available regardless of where pip put it. The docs copy
    (``docs/inputs/input_template.xlsx``) is the editable master; the two are kept
    byte-identical by a check in ``run_tests.py``.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "resources", "input_template.xlsx")


# in-memory material key -> 'mat' sheet header text. Shared by the writer and by the
# input-building scripts, so a new column is declared in exactly one place.
MAT_NUM_HEADERS = [
    ('gamma', 'g'), ('c', 'c'), ('phi', 'f'), ('cp', 'c/p'), ('r_elev', 'r-elev'),
    ('d', 'd'), ('psi', 'psi'),
    ('pow_a', 'powa'), ('pow_b', 'powb'), ('pow_c', 'powc'), ('pow_d', 'powd'),
    ('ru', 'ru'),
    ('sigma_gamma', 's(g)'), ('sigma_c', 's(c)'), ('sigma_phi', 's(f)'),
    ('sigma_cp', 's(c/p)'), ('sigma_d', 's(d)'), ('sigma_psi', 's(psi)'),
    ('hb_sci', 'hb_sci'), ('hb_gsi', 'hb_gsi'),
    ('hb_mi', 'hb_mi'), ('hb_d', 'hb_d'),
    ('k1', 'k1'), ('k2', 'k2'), ('alpha', 'alpha'),
    # v14: the curve-parameter pair is law-agnostic ('a'/'n'), serving both van
    # Genuchten and Gardner; Poisson's ratio moved off 'n' to 'nu' to break the
    # header collision that created.
    ('kr0', 'kr0'), ('h0', 'h0'), ('vg_a', 'a'), ('vg_n', 'n'),
    ('E', 'E'), ('nu', 'nu'),
]
# Optional numerics: written only when set (None must stay a blank cell -- e.g.
# a blank gsat means "fall back to gamma", which 0.0 would not).
MAT_OPT_NUM_HEADERS = [('gamma_sat', 'gsat')]
MAT_STR_HEADERS = [('option', 'option'), ('u', 'u'), ('unsat', 'unsat')]

# The 'mat' header row is located by scanning for its sentinel cells rather than
# assumed to be row 8: the header has moved before (pre-v8 templates put it on row 3,
# and none of those load today), and a legend row inserted above the table moves it
# again. Bounded scan -- the legend above the table is short.
_MAT_HEADER_SCAN_ROWS = 20


def _find_mat_header_row(cell):
    """Row index (1-based) of the 'mat' sheet header, found by its sentinel cells:
    column A == 'mat' and column B == 'name'.

    ``cell(row, col) -> value`` is a lookup callable, so this works against either an
    openpyxl worksheet or a header-less DataFrame.

    The row analogue of :data:`MAT_NUM_HEADERS`' by-name column lookup: neither the
    header's row nor its columns may be hardcoded, so inserting a legend row or a new
    property column cannot silently shift what gets read or written.

    Raises:
        ValueError: if no header row is found in the first 20 rows.
    """
    def _s(v):
        return '' if v is None else str(v).strip().lower()

    for r in range(1, _MAT_HEADER_SCAN_ROWS + 1):
        if _s(cell(r, 1)) == 'mat' and _s(cell(r, 2)) == 'name':
            return r
    raise ValueError(
        "The 'mat' sheet has no header row: expected a row whose first two cells are "
        f"'mat' and 'name' within the first {_MAT_HEADER_SCAN_ROWS} rows."
    )


def _read_mat_header_cols(filepath):
    """(header_row, cols) for the 'mat' sheet, read from the destination workbook.

    ``cols`` maps each header's text to its 1-based column index, so
    :func:`save_slope_data_to_xlsx` locates material columns BY NAME rather than by
    hardcoded position and adapts to the template version automatically (v11 inserted
    ``unsat`` and ``vg_a``/``vg_n``, shifting ``kr0``/``h0``/``E``/``nu``). Headers
    absent in an older template are simply skipped by the caller. ``header_row`` is
    found the same way -- see :func:`_find_mat_header_row`.
    """
    import openpyxl
    wb = openpyxl.load_workbook(filepath, read_only=True, data_only=False)
    try:
        ws = wb['mat']
        header_row = _find_mat_header_row(lambda r, c: ws.cell(row=r, column=c).value)
        cols = {}
        for c in range(1, ws.max_column + 1):
            v = ws.cell(row=header_row, column=c).value
            if v is not None and str(v).strip() != '':
                # underscore-insensitive, matching load_slope_data's mat parse
                cols[str(v).strip().replace('_', '')] = c
        return header_row, cols
    finally:
        wb.close()


def mat_header_cols(filepath):
    """Public accessor for the 'mat' sheet's ``(header_row, {header: column})`` map.

    Input-building scripts must write material cells through this rather than by
    hardcoded column number, or a future column insert shifts them silently.
    """
    return _read_mat_header_cols(filepath)


def _inplace_save_would_drop(filepath, materials):
    """True if editing ``filepath`` in place would silently drop material data
    because its 'mat' sheet lacks a column the data needs — e.g. a pre-v11 file
    (no ``unsat``/``vg_a``/``vg_n`` columns) being saved with a van Genuchten
    material. The caller then upgrades the file through the current template instead
    of dropping the data. Only flags when a material carries a NON-default value for
    a missing column, so old files that don't use the new feature still save in
    place (preserving their formatting)."""
    if not materials:
        return False
    have = set(_read_mat_header_cols(filepath))
    # (header, key, default) for columns that may be absent in older templates.
    checks = [('unsat', 'unsat', 'lf'), ('vg_a', 'vg_a', 0.0), ('vg_n', 'vg_n', 0.0)]
    for header, key, default in checks:
        if header in have:
            continue
        for m in materials:
            val = m.get(key)
            if val is None:
                continue
            if isinstance(default, str):
                if str(val).strip().lower() not in ('', 'nan', default):
                    return True
            else:
                try:
                    if float(val) != default:
                        return True
                except (TypeError, ValueError):
                    pass
    return False


def save_slope_data_to_xlsx(slope_data, filepath, template=None):
    """
    Write an in-memory ``slope_data`` dict back to an XSLOPE Excel input file.

    This is the inverse of :func:`load_slope_data`: it maps every editable input
    category (global parameters, materials, geometry, piezometric lines, failure
    surfaces, distributed loads, reinforcement, piles, and seepage boundary
    conditions) back into the template's sheet/cell layout. It builds on
    :func:`write_cells_to_xlsx`, which edits the workbook at the XML level and so
    preserves all formatting, formulas, charts, and drawings.

    Round-trip guarantee: for every input category,
    ``load_slope_data(<file written by this function>)`` reproduces the source
    ``slope_data``. Derived geometry (``ground_surface``, ``domain_polygon``,
    ``tcrack_surface``, and — for profile input — ``polygons``) is recomputed by
    the loader and therefore not written here.

    Parameters
    ----------
    slope_data : dict
        A dict in the form returned by :func:`load_slope_data`.
    filepath : str
        Destination ``.xlsx`` path.
    template : str, optional
        Path to a blank XSLOPE template to copy before writing (e.g. "Save As" from
        a chosen file). If ``None`` (the default, and the normal "Save"), the current
        standard template from :func:`default_template_path` is used. Either way the
        destination is written from a fresh, data-free template, so ``filepath`` need
        not exist and no stale rows can survive from a previous version of the file.

    Returns
    -------
    str
        ``filepath``.

    Notes
    -----
    Geometry is written to either the ``profile`` sheet (when ``profile_lines`` is
    present) or the ``polygon`` sheet (otherwise) — never both, matching the
    loader's mutual-exclusivity rule. Formula cells (e.g. the profile/polygon
    row-6 material-name XLOOKUP) are never written; the material IDs that feed
    them are, and :func:`write_cells_to_xlsx` flags the workbook for a full
    recalculation on open so the dependent formulas refresh.

    Circle surfaces are always written with ``Option = "Depth"`` (the loader
    collapses every circle to ``Xo/Yo/Depth/R`` with ``R = Yo - Depth``, so this
    reproduces the radius regardless of how the original was specified). The pile
    ``theta`` (``qp``) column is left blank because the loader auto-derives it
    from the pile endpoints.
    """
    # Always write into a FRESH copy of the standard template rather than editing the
    # destination in place. The template carries no data — only structural helpers the
    # loader ignores (index numbers in column A, dropdown reference lists) — so every
    # save rebuilds each input category from slope_data on a clean sheet. This is
    # robust by construction: a shortened list (a deleted reinforcement line, profile
    # point, material, …) leaves no stale/orphaned rows behind, and it also auto-
    # upgrades older files to the current template format. The trade-off is that
    # user-added custom formulas/formatting in the destination are not preserved.
    # Callers may pass an explicit `template` (e.g. Save As from a chosen file).
    if template is None:
        template = default_template_path()
    shutil.copy(template, filepath)

    def _f(v):
        return float(v)

    def _isnan(v):
        """True for an unset numeric field (NaN or None) — written as a blank cell."""
        if v is None:
            return True
        try:
            return float(v) != float(v)
        except (TypeError, ValueError):
            return False

    updates = {}

    # === main ===
    updates['main'] = {
        'D8': _f(slope_data['gamma_water']),
        'D9': _f(slope_data['tcrack_depth']),
        'D10': _f(slope_data['tcrack_water']),
        'D11': _f(slope_data['k_seismic']),
    }

    # === mat ===  Both the header ROW and its COLUMNS are located by name in the
    # destination file, never hardcoded, so the writer adapts to the template version
    # automatically (v11 inserted 'unsat' and 'vg_a'/'vg_n', shifting kr0/h0/E/nu; a
    # legend row above the table would shift the header row). A header absent in an
    # older template is skipped.
    mat_header_row, mat_cols = _read_mat_header_cols(filepath)
    mat = {}
    # mat_cols is keyed underscore-insensitively (see _read_mat_header_cols), so every
    # lookup must normalize the declared header the same way -- otherwise a header
    # written with an underscore ('hb_sci') misses its column and the writer SILENTLY
    # skips it, zeroing the property on the next load.
    def _col(header):
        return mat_cols.get(header.replace('_', ''))

    for idx, material in enumerate(slope_data.get('materials', [])):
        row = mat_header_row + 1 + idx
        if 'mat' in mat_cols:
            mat[cell_ref(row, mat_cols['mat'])] = idx + 1     # 1-based mat number
        if 'name' in mat_cols:
            mat[cell_ref(row, mat_cols['name'])] = str(material.get('name', ''))
        # option / u / unsat are strings; leave the cell blank when unset (the loader
        # reads an empty cell back as a default, so writing literal text is noise).
        for key, header in MAT_STR_HEADERS:
            col = _col(header)
            if col is None:
                continue
            val = material.get(key)
            if val is not None and str(val).strip().lower() not in ('', 'nan'):
                mat[cell_ref(row, col)] = str(val)
        for key, header in MAT_NUM_HEADERS:
            col = _col(header)
            if col is None:
                continue
            mat[cell_ref(row, col)] = _f(material.get(key, 0) or 0)
        for key, header in MAT_OPT_NUM_HEADERS:
            col = _col(header)
            if col is None:
                continue
            val = material.get(key)
            if val is not None:
                mat[cell_ref(row, col)] = _f(val)
    if mat:
        updates['mat'] = mat

    # === geometry: profile OR polygon (mutually exclusive, matching the loader) ===
    profile_lines = slope_data.get('profile_lines') or []
    if profile_lines:
        prof = {}
        md = slope_data.get('max_depth')
        prof['B2'] = _f(md) if md is not None else 0.0
        for n, line in enumerate(profile_lines):              # n is 0-based
            x_col = 1 + n * 3                                  # A, D, G, ...
            y_col = x_col + 1
            # Write the row-4 block header. load_slope_data() detects how many
            # profile lines exist by scanning for a non-empty header here, so this
            # makes the file self-describing rather than relying on the template's
            # pre-labeled blocks.
            prof[cell_ref(4, x_col)] = f"Profile Line #{n + 1}"
            mat_id = line.get('mat_id')
            if mat_id is not None:
                prof[cell_ref(5, y_col)] = int(mat_id) + 1    # 0-based -> 1-based
            for i, (x, y) in enumerate(line['coords']):
                prof[cell_ref(8 + i, x_col)] = _f(x)
                prof[cell_ref(8 + i, y_col)] = _f(y)
        updates['profile'] = prof
    else:
        polygons = slope_data.get('polygons') or []
        if polygons:
            poly_u = {}
            for n, pdict in enumerate(polygons):
                x_col = 1 + n * 3
                y_col = x_col + 1
                poly_u[cell_ref(4, x_col)] = f"Polygon #{n + 1}"   # block header (see profile)
                mat_id = pdict.get('mat_id')
                if mat_id is not None:
                    poly_u[cell_ref(5, y_col)] = int(mat_id) + 1
                coords = list(pdict['polygon'].exterior.coords)
                if len(coords) >= 2 and coords[0] == coords[-1]:
                    coords = coords[:-1]                       # loader closes implicitly
                for i, (x, y) in enumerate(coords):
                    poly_u[cell_ref(8 + i, x_col)] = _f(x)
                    poly_u[cell_ref(8 + i, y_col)] = _f(y)
            updates['polygon'] = poly_u

    # === piezo ===  (v13 layout: Type row at Excel row 3 — 'piezo' static head
    # or 'phreatic' cos^2 correction — x/y headers row 4, data from row 5)
    piezo = {}
    piezo[cell_ref(3, 2)] = 'phreatic' if slope_data.get('piezo_phreatic') else 'piezo'
    piezo[cell_ref(3, 5)] = 'phreatic' if slope_data.get('piezo_phreatic2') else 'piezo'
    for i, (x, y) in enumerate(slope_data.get('piezo_line') or []):
        piezo[cell_ref(5 + i, 1)] = _f(x)
        piezo[cell_ref(5 + i, 2)] = _f(y)
    for i, (x, y) in enumerate(slope_data.get('piezo_line2') or []):
        piezo[cell_ref(5 + i, 4)] = _f(x)
        piezo[cell_ref(5 + i, 5)] = _f(y)
    if piezo:
        updates['piezo'] = piezo

    # === circles ===  (header row 2, data rows 3+; always written as Option="Depth")
    circ = {}
    for n, c in enumerate(slope_data.get('circles') or []):
        row = 3 + n
        circ[cell_ref(row, 1)] = n + 1
        circ[cell_ref(row, 2)] = _f(c['Xo'])
        circ[cell_ref(row, 3)] = _f(c['Yo'])
        circ[cell_ref(row, 4)] = 'Depth'
        circ[cell_ref(row, 5)] = _f(c['Depth'])
    if circ:
        updates['circles'] = circ

    # === non-circ ===  (data rows 3+; cols A=X, B=Y, C=Movement)
    nonc = {}
    for i, p in enumerate(slope_data.get('non_circ') or []):
        row = 3 + i
        nonc[cell_ref(row, 1)] = _f(p['X'])
        nonc[cell_ref(row, 2)] = _f(p['Y'])
        mv = p.get('Movement')
        if mv is not None and not (isinstance(mv, float) and pd.isna(mv)):
            nonc[cell_ref(row, 3)] = str(mv)
    if nonc:
        updates['non-circ'] = nonc

    # === dloads / dloads (2) ===  (3-col blocks from col B, +1 gap; data from row 4)
    def _dload_updates(blocks):
        u = {}
        for n, block in enumerate(blocks):
            x_col = 2 + n * 4                                  # B, F, J, ...
            for i, pt in enumerate(block):
                u[cell_ref(4 + i, x_col)] = _f(pt['X'])
                u[cell_ref(4 + i, x_col + 1)] = _f(pt['Y'])
                u[cell_ref(4 + i, x_col + 2)] = _f(pt['Normal'])
        return u
    d1 = _dload_updates(slope_data.get('dloads') or [])
    if d1:
        updates['dloads'] = d1
    d2 = _dload_updates(slope_data.get('dloads2') or [])
    if d2:
        updates['dloads (2)'] = d2

    # === reinforce ===  (raw endpoint form in 'reinforcement_lines' round-trips)
    # v12 layout: # | Label | x1 y1 x2 y2 | Type Dir Appl | Tmax Lp1 Lp2 Tend1
    # Tend2 Spacing | Tres E Area. Capacity terms were divided by Spacing at load,
    # so they are multiplied back here -- the file carries per-element values.
    # Dir/Appl carry in-sheet default formulas driven by Type: those cells are
    # written ONLY when the value differs from what the Type preset (or the
    # generic default) would produce, so the formulas survive a round-trip and
    # the preset behavior stays live for hand editing.
    _REINF_PRESETS = {
        'geosynthetic': ('tangent', 'active'),
        'nail':         ('axial',   'passive'),
        'tieback':      ('axial',   'active'),
        'anchor':       ('axial',   'active'),
    }
    reinf = {}
    for n, r in enumerate(slope_data.get('reinforcement_lines') or []):
        row = 3 + n
        sp = float(r.get('spacing', 1.0) or 1.0)
        reinf.update({
            cell_ref(row, 1): n + 1,
            cell_ref(row, 2): str(r.get('label', f"Line {n + 1}")),
            cell_ref(row, 3): _f(r['x1']), cell_ref(row, 4): _f(r['y1']),
            cell_ref(row, 5): _f(r['x2']), cell_ref(row, 6): _f(r['y2']),
            cell_ref(row, 10): _f(r['t_max']) * sp,
            cell_ref(row, 11): _f(r['lp1']), cell_ref(row, 12): _f(r['lp2']),
            cell_ref(row, 13): _f(r.get('tend1', 0.0)) * sp,
            cell_ref(row, 14): _f(r.get('tend2', 0.0)) * sp,
            cell_ref(row, 15): sp,
            # unset Tres round-trips as a BLANK cell, not a literal NaN
            cell_ref(row, 16): (None if _isnan(r.get('t_res'))
                                else _f(r.get('t_res', 0.0)) * sp),
            cell_ref(row, 17): _f(r['E']), cell_ref(row, 18): _f(r['area']) * sp,
        })
        rtype = str(r.get('type', '') or '')
        d_def, a_def = _REINF_PRESETS.get(rtype, ('tangent', 'active'))
        if rtype:
            reinf[cell_ref(row, 7)] = rtype.capitalize()
        if str(r.get('dir', d_def)) != d_def:
            reinf[cell_ref(row, 8)] = str(r['dir']).capitalize()
        if str(r.get('appl', a_def)) != a_def:
            reinf[cell_ref(row, 9)] = str(r['appl']).capitalize()
    if reinf:
        updates['reinforce'] = reinf

    # === piles ===  (v12 layout: # | Label | x1 y1 x2 y2 | H qp Appl | D S Vcap
    # Mcap | E I Area Fixity; header row 2, data rows 3+; qp left blank — auto-
    # derived by the loader from the pile endpoints)
    piles_u = {}
    for n, p in enumerate(slope_data.get('pile_lines') or []):
        row = 3 + n
        piles_u[cell_ref(row, 1)] = n + 1
        piles_u[cell_ref(row, 2)] = str(p.get('label', f"Pile {n + 1}"))
        piles_u[cell_ref(row, 3)] = _f(p['x1'])
        piles_u[cell_ref(row, 4)] = _f(p['y1'])
        piles_u[cell_ref(row, 5)] = _f(p['x2'])
        piles_u[cell_ref(row, 6)] = _f(p['y2'])
        for key, col in [('H', 7), ('D_pile', 10), ('S', 11), ('V_cap', 12),
                         ('M_cap', 13), ('E', 14), ('I', 15), ('area', 16)]:
            val = p.get(key)
            if val is not None:
                piles_u[cell_ref(row, col)] = _f(val)
        if str(p.get('appl', 'active')) == 'passive':
            piles_u[cell_ref(row, 9)] = 'Passive'
        piles_u[cell_ref(row, 17)] = str(p.get('fixity', 'free'))
    if piles_u:
        updates['piles'] = piles_u

    # === lloads ===  (v12: # | Label | x | y | P | Angle; data rows 3+)
    lloads_u = {}
    for n, ll in enumerate(slope_data.get('line_loads') or []):
        row = 3 + n
        lloads_u[cell_ref(row, 1)] = n + 1
        lloads_u[cell_ref(row, 2)] = str(ll.get('label', f"Load {n + 1}"))
        lloads_u[cell_ref(row, 3)] = _f(ll['x'])
        lloads_u[cell_ref(row, 4)] = _f(ll['y'])
        lloads_u[cell_ref(row, 5)] = _f(ll['P'])
        lloads_u[cell_ref(row, 6)] = _f(ll.get('angle', -90.0))
    if lloads_u:
        updates['lloads'] = lloads_u

    # === seep bc / seep bc (2) ===
    def _seep_updates(bc):
        u = {}
        for i, (x, y) in enumerate(bc.get('exit_face') or []):
            u[cell_ref(5 + i, 2)] = _f(x)                     # B
            u[cell_ref(5 + i, 3)] = _f(y)                     # C
        blocks = [('head', b['head'], b['coords'])
                  for b in bc.get('specified_heads') or []]
        blocks += [('flux', b['flux'], b['coords'])
                   for b in bc.get('specified_fluxes') or []]
        for k, (kind, value, coords) in enumerate(blocks):
            x_col = 5 + k * 3                                 # E, H, K, ...
            y_col = x_col + 1
            u[cell_ref(3, x_col)] = kind                      # type cell (head/flux)
            u[cell_ref(3, y_col)] = _f(value)                 # head or flux value
            for i, (x, y) in enumerate(coords):
                u[cell_ref(5 + i, x_col)] = _f(x)
                u[cell_ref(5 + i, y_col)] = _f(y)
        return u
    s1 = _seep_updates(slope_data.get('seepage_bc') or {})
    if s1:
        updates['seep bc'] = s1
    s2 = _seep_updates(slope_data.get('seepage_bc2') or {})
    if s2:
        updates['seep bc (2)'] = s2

    updates = {k: v for k, v in updates.items() if v}
    write_cells_to_xlsx(filepath, updates)
    return filepath


def save_data_to_pickle(data, filepath):
    """
    Save a data object to a pickle file.
    
    This function serializes the data object and saves it to the specified filepath.
    Useful for saving processed data from Excel templates for later use.
    
    Parameters:
        data: The data object to save (typically a dictionary from load_slope_data)
        filepath (str): The file path where the pickle file should be saved
        
    Returns:
        None
        
    Raises:
        IOError: If the file cannot be written
        PickleError: If the data cannot be serialized
    """
    try:
        with open(filepath, 'wb') as f:
            pickle.dump(data, f)
    except Exception as e:
        raise IOError(f"Failed to save data to pickle file '{filepath}': {e}")


def load_data_from_pickle(filepath):
    """
    Load a data object from a pickle file.
    
    This function deserializes a data object from the specified pickle file.
    Useful for loading previously saved data without re-processing Excel templates.
    
    Parameters:
        filepath (str): The file path of the pickle file to load
        
    Returns:
        The deserialized data object (typically a dictionary)
        
    Raises:
        FileNotFoundError: If the pickle file doesn't exist
        IOError: If the file cannot be read
        PickleError: If the data cannot be deserialized
    """
    try:
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        return data
    except FileNotFoundError:
        raise FileNotFoundError(f"Pickle file not found: '{filepath}'")
    except Exception as e:
        raise IOError(f"Failed to load data from pickle file '{filepath}': {e}")
    

def print_dictionary(dictionary):
    """
    Print the contents of a dictionary to the console.
    This can be used for slope_data, seep_data, or any other dictionary.
    """
    for key, value in dictionary.items():
        print(f"\n=== {key} ===")
        if isinstance(value, list):
            for item in value:
                print(item)
        else:
            print(value)


# ---------------------------------------------------------------------------
# Surgical xlsx writer
# ---------------------------------------------------------------------------
# Writes individual cell values into an existing xlsx template by editing the
# sheet XML in place, rather than re-saving the whole workbook (which would
# discard formulas, formatting, and named ranges). Used to populate templates
# programmatically — e.g. from the polygon-sheet builders and the planned CAD
# importer.
#
# Two pitfalls this handles (do not "simplify" them away):
#   * Writing a value into a cell that some formula depends on invalidates the
#     cached calcChain. Excel then "recovers" the sheet and silently discards
#     the edits. We force fullCalcOnLoad and delete calcChain.xml so Excel
#     rebuilds it. (See the xlsx calcChain/recalc note in project memory.)
#   * Never write a plain value into a cell that itself holds a formula — that
#     breaks calcChain.xml. Callers must target precedent/value cells only.
#
# Requires the `zip` CLI (used to replace single entries inside the archive in
# place; Python's zipfile cannot update an existing member without rewriting
# the whole file).

_XLSX_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
_XLSX_R_NS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
_XLSX_PKG_NS = 'http://schemas.openxmlformats.org/package/2006/relationships'


def _col_num_to_letter(n):
    result = ''
    while n > 0:
        n, rem = divmod(n - 1, 26)
        result = chr(65 + rem) + result
    return result


def cell_ref(row, col):
    """Return an A1-style cell reference for 1-based (row, col)."""
    return f'{_col_num_to_letter(col)}{row}'


def _parse_cell_ref(ref):
    m = re.match(r'^([A-Z]+)(\d+)$', ref)
    col_str, row = m.group(1), int(m.group(2))
    col = 0
    for c in col_str:
        col = col * 26 + (ord(c) - ord('A') + 1)
    return row, col


def _modify_existing_cell(cell_xml, value):
    open_match = re.match(r'(<c\s[^>]*?)(/?>)', cell_xml)
    if not open_match:
        return cell_xml
    open_tag_attrs = re.sub(r'\s+t="[^"]*"', '', open_match.group(1))
    # Blank the cell for None / NaN / inf. A numeric cell holding "nan" in particular
    # makes openpyxl's reader raise int('nan') on the next load (silent corruption).
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return f'{open_tag_attrs}/>'
    if isinstance(value, float):
        value = round(value, 10)
    if isinstance(value, str):
        return (f'{open_tag_attrs} t="inlineStr"><is><t>'
                f'{xml_escape(value)}</t></is></c>')
    else:
        return f'{open_tag_attrs}><v>{value}</v></c>'


def _build_new_cell(ref, value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return f'<c r="{ref}"/>'          # blank cell (None / NaN / inf)
    if isinstance(value, float):
        value = round(value, 10)
    if isinstance(value, str):
        return f'<c r="{ref}" t="inlineStr"><is><t>{xml_escape(value)}</t></is></c>'
    else:
        return f'<c r="{ref}"><v>{value}</v></c>'


def _modify_sheet_xml(xml_bytes, cells):
    xml_text = xml_bytes.decode('utf-8')
    rows_data = {}
    for ref, value in cells.items():
        row_num, col_num = _parse_cell_ref(ref)
        rows_data.setdefault(row_num, []).append((ref, col_num, value))
    for row_num in rows_data:
        rows_data[row_num].sort(key=lambda x: x[1])
    for row_num, row_cells in sorted(rows_data.items()):
        row_pattern = re.compile(
            r'(<row\s+r="%d"[^>]*>)(.*?)(</row>)' % row_num, re.DOTALL)
        row_match = row_pattern.search(xml_text)
        if row_match:
            row_open = row_match.group(1)
            row_content = row_match.group(2)
            row_close = row_match.group(3)
            for ref, col_num, value in row_cells:
                cell_pattern = re.compile(
                    r'<c\s+r="%s"(?:\s[^>]*)*/>' % re.escape(ref) +
                    r'|<c\s+r="%s"(?:\s[^>]*)?>.*?</c>' % re.escape(ref), re.DOTALL)
                cell_match = cell_pattern.search(row_content)
                if cell_match:
                    new_cell = _modify_existing_cell(cell_match.group(0), value)
                    row_content = (row_content[:cell_match.start()] +
                                   new_cell + row_content[cell_match.end():])
                else:
                    row_content = row_content + _build_new_cell(ref, value)
            new_row = row_open + row_content + row_close
            xml_text = xml_text[:row_match.start()] + new_row + xml_text[row_match.end():]
        else:
            cells_xml = ''.join(_build_new_cell(ref, value) for ref, _, value in row_cells)
            new_row = f'<row r="{row_num}">{cells_xml}</row>'
            sd_close = xml_text.rfind('</sheetData>')
            xml_text = xml_text[:sd_close] + new_row + xml_text[sd_close:]
    return xml_text.encode('utf-8')


def _reset_view(xml_bytes):
    """Reset a sheet's saved scroll position/selection to A1 so populated cells are
    visible on open. The template ships some sheets (e.g. polygon) scrolled far right
    (topLeftCell="AA1", activeCell="AA8") — Excel then opens to an empty-looking region
    and the populated left-hand columns appear blank until you scroll back."""
    t = xml_bytes.decode('utf-8')
    t = re.sub(r'\s+topLeftCell="[^"]*"', '', t)
    t = re.sub(r'<selection\b[^>]*/>', '<selection activeCell="A1" sqref="A1"/>', t)
    return t.encode('utf-8')


def _force_full_recalc(xml_text):
    m = re.search(r'<calcPr\b[^>]*?/?>', xml_text)
    if m:
        tag = m.group(0)
        if 'fullCalcOnLoad=' in tag:
            new_tag = re.sub(r'fullCalcOnLoad="[^"]*"', 'fullCalcOnLoad="1"', tag)
        elif tag.endswith('/>'):
            new_tag = tag[:-2] + ' fullCalcOnLoad="1"/>'
        else:
            new_tag = tag[:-1] + ' fullCalcOnLoad="1">'
        return xml_text[:m.start()] + new_tag + xml_text[m.end():]
    insert = '<calcPr calcId="191029" fullCalcOnLoad="1"/>'
    idx = xml_text.find('</sheets>')
    pos = (idx + len('</sheets>')) if idx != -1 else xml_text.find('</workbook>')
    return xml_text[:pos] + insert + xml_text[pos:]


def _drop_calcchain(tmpdir, paths_to_zip, zf_read):
    """Stage edits that remove xl/calcChain.xml and its two references.

    The cached calcChain becomes stale the moment we change a formula's precedent
    cell at the XML level. Excel then "recovers" the affected sheet and discards
    our edits (symptom: a populated sheet opens blank). Deleting calcChain.xml and
    its references makes Excel rebuild it from scratch — safe because we also set
    fullCalcOnLoad="1", so every formula is recomputed on open. Returns True if a
    calcChain part was present and needs deleting from the archive."""
    ct = zf_read('[Content_Types].xml').decode('utf-8')
    if 'calcChain' not in ct:
        return False
    ct = re.sub(r'<Override\s+PartName="/xl/calcChain\.xml"[^>]*/>', '', ct)
    out = os.path.join(tmpdir, '[Content_Types].xml')
    with open(out, 'wb') as f:
        f.write(ct.encode('utf-8'))
    paths_to_zip.append('[Content_Types].xml')

    rels = zf_read('xl/_rels/workbook.xml.rels').decode('utf-8')
    rels = re.sub(r'<Relationship\s+[^>]*Target="calcChain\.xml"[^>]*/>', '', rels)
    out = os.path.join(tmpdir, 'xl/_rels/workbook.xml.rels')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'wb') as f:
        f.write(rels.encode('utf-8'))
    paths_to_zip.append('xl/_rels/workbook.xml.rels')
    return True


def write_cells_to_xlsx(filepath, updates):
    """Write cell values into an existing xlsx template in place.

    Parameters
    ----------
    filepath : str
        Path to the xlsx file to modify (edited in place).
    updates : dict
        Mapping of {sheet_name: {cell_ref: value}}. Use cell_ref(row, col) to
        build the A1-style references. Values may be str, int, or float. Target
        only value/precedent cells — never cells that hold formulas.
    """
    with zipfile.ZipFile(filepath) as zf:
        wb_xml = etree.fromstring(zf.read('xl/workbook.xml'))
        rels_xml = etree.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
        rid_map = {r.get('Id'): r.get('Target')
                   for r in rels_xml.iter('{%s}Relationship' % _XLSX_PKG_NS)}
        sheet_paths = {}
        for s in wb_xml.iter('{%s}sheet' % _XLSX_NS):
            rid = s.get('{%s}id' % _XLSX_R_NS)
            if rid and rid in rid_map:
                sheet_paths[s.get('name')] = f'xl/{rid_map[rid]}'
    tmpdir = tempfile.mkdtemp()
    abs_filepath = os.path.abspath(filepath)
    try:
        paths_to_zip = []
        for sheet_name, cells in updates.items():
            path = sheet_paths[sheet_name]
            with zipfile.ZipFile(filepath) as zf:
                orig_xml = zf.read(path)
            modified_xml = _reset_view(_modify_sheet_xml(orig_xml, cells))
            out_path = os.path.join(tmpdir, path)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path, 'wb') as f:
                f.write(modified_xml)
            paths_to_zip.append(path)
        with zipfile.ZipFile(filepath) as zf:
            wb_text = zf.read('xl/workbook.xml').decode('utf-8')
        wb_out = os.path.join(tmpdir, 'xl/workbook.xml')
        os.makedirs(os.path.dirname(wb_out), exist_ok=True)
        with open(wb_out, 'wb') as f:
            f.write(_force_full_recalc(wb_text).encode('utf-8'))
        paths_to_zip.append('xl/workbook.xml')
        with zipfile.ZipFile(filepath) as zf:
            drop_cc = _drop_calcchain(tmpdir, paths_to_zip, zf.read)
        for path in paths_to_zip:
            subprocess.run(['zip', abs_filepath, path],
                           cwd=tmpdir, capture_output=True, text=True)
        if drop_cc:
            subprocess.run(['zip', '-d', abs_filepath, 'xl/calcChain.xml'],
                           capture_output=True, text=True)
    finally:
        shutil.rmtree(tmpdir)