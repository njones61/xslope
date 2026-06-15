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
    top-left corner to the top-right corner (see plans/plan_polygons.md §4).

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
    mat_df = xls.parse('mat', header=7)  # header=7 because the header row is row 8 in Excel (0-indexed row 7)
    materials = []

    def _num(x):
        v = pd.to_numeric(x, errors="coerce")
        return float(v) if pd.notna(v) else 0.0

    # Read materials row by row until we encounter an empty material name (Column B)
    # Data starts at Excel row 9 (0-indexed row 0 after header=7)
    for i in range(len(mat_df)):
        row = mat_df.iloc[i]
        
        # Check if material name (Column B) is empty - stop reading if empty
        material_name = row.get('name', '')
        if pd.isna(material_name) or str(material_name).strip() == '':
            break  # Stop reading when we encounter an empty material name
        
        # For seep workflows, 'g' (unit weight) and shear strength properties are not required.
        # A material row is considered "missing" only if Excel columns C:X are empty.
        # (Excel A:B are number and name; C:X contain the actual property fields.)
        start_col = 2  # C
        end_col = min(mat_df.shape[1], 24)  # X is column 24 (1-based) -> index 23, so slice end is 24
        c_to_x_empty = True if start_col >= end_col else row.iloc[start_col:end_col].isna().all()
        if c_to_x_empty:
            # Excel row number: header is on row 8, first data row is row 9
            excel_row = i + 9
            raise ValueError(
                "CRITICAL ERROR: Material row has empty property fields. "
                f"Material '{material_name}' (Excel row {excel_row}) is blank in columns C:X."
            )

        materials.append({
            "name": str(material_name).strip(),
            "gamma": _num(row.get("g", 0)),
            "option": str(row.get('option', '')).strip().lower(),
            "c": _num(row.get('c', 0)),
            "phi": _num(row.get('f', 0)),
            "cp": _num(row.get('c/p', 0)),
            "r_elev": _num(row.get('r-elev', 0)),
            "d": _num(row.get('d', 0)) if pd.notna(row.get('d')) else 0,
            "psi": _num(row.get('psi', 0)) if pd.notna(row.get('psi')) else 0,
            "u": str(row.get('u', 'none')).strip().lower(),
            "sigma_gamma": _num(row.get('s(g)', 0)),
            "sigma_c": _num(row.get('s(c)', 0)),
            "sigma_phi": _num(row.get('s(f)', 0)),
            "sigma_cp": _num(row.get('s(c/p)', 0)),
            "sigma_d": _num(row.get('s(d)', 0)),
            "sigma_psi": _num(row.get('s(psi)', 0)),
            "k1": _num(row.get('k1', 0)),
            "k2": _num(row.get('k2', 0)),
            "alpha": _num(row.get('alpha', 0)),
            "kr0" : _num(row.get('kr0', 0)),
            "h0" : _num(row.get('h0', 0)),
            "E": _num(row.get('E', 0)),
            "nu": _num(row.get('n', 0))
        })

    # === UNIFIED POLYGON REPRESENTATION ===
    # Geometry is always represented internally as material-zone polygons (see
    # plans/plan_polygons.md §5.1, §9.2). If the 'polygon' sheet is populated it is
    # used directly; otherwise the profile lines are converted to polygons via
    # build_polygons(). All downstream code (slicing, search) works from polygons.
    polygons_from_sheet = _parse_polygon_sheet(xls, materials)

    if polygons_from_sheet:
        if profile_lines:
            raise ValueError(
                "Both the 'profile' and 'polygon' sheets contain data. Use one "
                "geometry method, not both.")
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

    # Read first piezometric line (columns A:B, starting at row 4, Excel row 4 = index 3)
    # Keep reading until we encounter an empty row
    start_row = 3  # Excel row 4 (0-indexed row 3)
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
    reinforce_df = xls.parse('reinforce', header=1)  # Header in row 2 (0-indexed row 1)
    reinforce_lines = []        # LEM format: list of point lists with tension distributions
    reinforcement_lines = []    # FEM format: list of dicts with raw line endpoints and properties

    # Process rows starting from row 3 (Excel) which is 0-indexed row 0 in pandas after header=1
    # Keep reading until we encounter an empty value in column B
    for i, row in reinforce_df.iterrows():
        # Check if column B (x1 coordinate) is empty - stop reading if empty
        if pd.isna(row.iloc[1]):
            break  # Stop reading when column B is empty
        
        # Check if other required coordinates are present
        if pd.isna(row.iloc[2]) or pd.isna(row.iloc[3]) or pd.isna(row.iloc[4]):
            continue  # Skip rows with incomplete coordinate data
            
        # If coordinates are present, check for required parameters (Tmax, Lp1, Lp2)
        if pd.isna(row.iloc[5]) or pd.isna(row.iloc[7]) or pd.isna(row.iloc[8]):
            raise ValueError(f"Reinforcement line in row {i + 3} has coordinates but missing required parameters (Tmax, Lp1, Lp2). All three must be specified.")
            
        try:
            # Extract coordinates and parameters
            x1, y1 = float(row.iloc[1]), float(row.iloc[2])  # Columns B, C
            x2, y2 = float(row.iloc[3]), float(row.iloc[4])  # Columns D, E
            Tmax = float(row.iloc[5])  # Column F
            Tres = float(row.iloc[6])  # Column G
            Lp1 = float(row.iloc[7]) if not pd.isna(row.iloc[7]) else 0.0   # Column H
            Lp2 = float(row.iloc[8]) if not pd.isna(row.iloc[8]) else 0.0   # Column I
            E = float(row.iloc[9])     # Column J
            Area = float(row.iloc[10]) # Column K

            # Store raw line data for FEM (endpoints + properties)
            reinforcement_lines.append({
                "x1": x1, "y1": y1, "x2": x2, "y2": y2,
                "t_max": Tmax, "t_res": Tres,
                "lp1": Lp1, "lp2": Lp2,
                "E": E, "area": Area,
            })

            # Calculate line length and direction
            import math
            line_length = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if line_length == 0:
                continue  # Skip zero-length lines
                
            # Unit vector from (x1,y1) to (x2,y2)
            dx = (x2 - x1) / line_length
            dy = (y2 - y1) / line_length
            
            line_points = []
            
            # Handle different cases based on pullout lengths
            if Lp1 + Lp2 >= line_length:
                # Line too short - create single interior point
                if Lp1 == 0 and Lp2 == 0:
                    # Both ends anchored - uniform tension
                    line_points = [
                        {"X": x1, "Y": y1, "T": Tmax, "Tres": Tres, "E": E, "Area": Area},
                        {"X": x2, "Y": y2, "T": Tmax, "Tres": Tres, "E": E, "Area": Area}
                    ]
                else:
                    # Find equilibrium point where tensions are equal
                    # T1 = Tmax * d1/Lp1, T2 = Tmax * d2/Lp2
                    # At equilibrium: d1/Lp1 = d2/Lp2 and d1 + d2 = line_length
                    if Lp1 == 0:
                        # End 1 anchored, all tension at end 1
                        line_points = [
                            {"X": x1, "Y": y1, "T": Tmax, "Tres": Tres, "E": E, "Area": Area},
                            {"X": x2, "Y": y2, "T": 0.0, "Tres": 0, "E": E, "Area": Area}
                        ]
                    elif Lp2 == 0:
                        # End 2 anchored, all tension at end 2
                        line_points = [
                            {"X": x1, "Y": y1, "T": 0.0, "Tres": 0, "E": E, "Area": Area},
                            {"X": x2, "Y": y2, "T": Tmax, "Tres": Tres, "E": E, "Area": Area}
                        ]
                    else:
                        # Both ends have pullout - find equilibrium point
                        ratio_sum = 1.0/Lp1 + 1.0/Lp2
                        d1 = line_length / (Lp2 * ratio_sum)
                        d2 = line_length / (Lp1 * ratio_sum)
                        T_eq = Tmax * d1 / Lp1  # = Tmax * d2 / Lp2
                        
                        # Interior point location
                        x_int = x1 + d1 * dx
                        y_int = y1 + d1 * dy
                        
                        line_points = [
                            {"X": x1, "Y": y1, "T": 0.0, "Tres": 0, "E": E, "Area": Area},
                            {"X": x_int, "Y": y_int, "T": T_eq, "Tres": Tres, "E": E, "Area": Area},
                            {"X": x2, "Y": y2, "T": 0.0, "Tres": 0, "E": E, "Area": Area}
                        ]
            else:
                # Normal case - line long enough for 4 points
                points_to_add = []
                
                # Point 1: Start point
                points_to_add.append((x1, y1, 0.0, 0.0))
                
                # Point 2: At distance Lp1 from start (if Lp1 > 0)
                if Lp1 > 0:
                    x_p2 = x1 + Lp1 * dx
                    y_p2 = y1 + Lp1 * dy
                    points_to_add.append((x_p2, y_p2, Tmax, Tres))
                else:
                    # Lp1 = 0, so start point gets Tmax tension
                    points_to_add[0] = (x1, y1, Tmax, Tres)
                
                # Point 3: At distance Lp2 back from end (if Lp2 > 0)
                if Lp2 > 0:
                    x_p3 = x2 - Lp2 * dx
                    y_p3 = y2 - Lp2 * dy
                    points_to_add.append((x_p3, y_p3, Tmax, Tres))
                else:
                    # Lp2 = 0, so end point gets Tmax tension
                    pass  # Will be handled when adding end point
                
                # Point 4: End point
                if Lp2 > 0:
                    points_to_add.append((x2, y2, 0.0, 0.0))
                else:
                    points_to_add.append((x2, y2, Tmax, Tres))
                
                # Remove duplicate points (same x,y coordinates)
                unique_points = []
                tolerance = 1e-6
                for x, y, T, Tres in points_to_add:
                    is_duplicate = False
                    for ux, uy, uT, uTres in unique_points:
                        if abs(x - ux) < tolerance and abs(y - uy) < tolerance:
                            # Update tension to maximum value at this location
                            for i, (px, py, pT, pTres) in enumerate(unique_points):
                                if abs(x - px) < tolerance and abs(y - py) < tolerance:
                                    unique_points[i] = (px, py, max(pT, T), max(pTres, Tres))
                            is_duplicate = True
                            break
                    if not is_duplicate:
                        unique_points.append((x, y, T, Tres))
                
                # Convert to required format
                line_points = [{"X": x, "Y": y, "T": T, "Tres": Tres, "E": E, "Area": Area} for x, y, T, Tres in unique_points]
            
            if len(line_points) >= 2:
                reinforce_lines.append(line_points)
                
        except Exception as e:
            raise ValueError(f"Error processing reinforcement line in row {row.name + 3}: {e}")


    # === PILE LINES ===
    pile_lines = []
    if 'piles' in xls.sheet_names:
        piles_df = xls.parse('piles', header=1)
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
                H = float(row['H']) if pd.notna(row.get('H')) else None
                if pd.notna(row.get('theta')):
                    theta_p = float(row['theta'])
                else:
                    # Auto-compute: perpendicular to pile axis (0 for vertical)
                    dx = x2 - x1
                    dy = y2 - y1
                    theta_p = np.degrees(np.arctan2(dx, -dy))
                D_pile = float(row['D']) if pd.notna(row.get('D')) else None
                S = float(row['S']) if pd.notna(row.get('S')) else None
                E_pile = float(row['E']) if pd.notna(row.get('E')) else None
                I_pile = float(row['I']) if pd.notna(row.get('I')) else None
                area = float(row['Area']) if pd.notna(row.get('Area')) else None
                V_cap = float(row['Vcap']) if pd.notna(row.get('Vcap')) else None
                M_cap = float(row['Mcap']) if pd.notna(row.get('Mcap')) else None
                fixity_raw = str(row['Fixity']).strip().lower() if pd.notna(row.get('Fixity')) else 'free'
                if fixity_raw not in ('free', 'fixed'):
                    raise ValueError(f"Fixity must be 'free' or 'fixed', got '{fixity_raw}'")
                fixity = fixity_raw
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
                    "label": label,
                })
            except Exception as e:
                raise ValueError(f"Error processing pile in row {i + 3}: {e}")


    # === SEEPAGE ANALYSIS BOUNDARY CONDITIONS ===
    # Read first set from "seep bc" sheet
    seep_df = xls.parse('seep bc', header=None)
    seepage_bc = {"specified_heads": [], "exit_face": []}
    
    # Exit Face BC: starts at B5 (row 4, columns 1 and 2), continues down until empty x value
    exit_coords = []
    exit_start_row = 4  # Excel row 5 (0-indexed row 4)
    exit_x_col = 1  # Column B
    exit_y_col = 2  # Column C
    
    row = exit_start_row
    while row < seep_df.shape[0]:
        try:
            x_val = seep_df.iloc[row, exit_x_col]
            y_val = seep_df.iloc[row, exit_y_col]
            
            # Stop at first empty x value
            if pd.isna(x_val):
                break
            
            # If x is present, try to convert (y can be empty but we'll still add the point)
            if pd.notna(x_val) and pd.notna(y_val):
                exit_coords.append((float(x_val), float(y_val)))
        except:
            break
        row += 1
    seepage_bc["exit_face"] = exit_coords
    
    # Specified Head BCs: start at columns E:F, then H:I, etc.
    # Head value is in row 3 (index 2), XY values start at row 5 (index 4)
    # Keep reading to the right until head value in row 3 is empty
    head_row = 2  # Excel row 3 (0-indexed row 2)
    data_start_row = 4  # Excel row 5 (0-indexed row 4)
    col = 4  # Start with column E (index 4)
    
    while col < seep_df.shape[1]:
        x_col = col
        y_col = col + 1
        head_col = col + 1  # Head value is in the Y column (F, I, L, etc.)
        
        # Check if head value in row 3 is empty - stop reading if empty
        if seep_df.shape[0] <= head_row:
            break
        head_val = seep_df.iloc[head_row, head_col]
        if pd.isna(head_val):
            break  # Stop reading when head value is empty
        
        # Read XY coordinates starting from row 5, continue down until empty
        coords = []
        row = data_start_row
        while row < seep_df.shape[0]:
            try:
                x_val = seep_df.iloc[row, x_col]
                y_val = seep_df.iloc[row, y_col]
                
                # Stop at first empty x value
                if pd.isna(x_val):
                    break
                
                # If x is present, try to convert
                if pd.notna(x_val) and pd.notna(y_val):
                    coords.append((float(x_val), float(y_val)))
            except:
                break
            row += 1
        
        if coords:  # Only add if we have coordinates
            seepage_bc["specified_heads"].append({"head": float(head_val), "coords": coords})
        
        # Move to next specified head BC (skip 3 columns: E->H, H->K, etc.)
        col += 3
    
    # Read second set from "seep bc (2)" sheet
    seepage_bc2 = {"specified_heads": [], "exit_face": []}
    try:
        seep_df2 = xls.parse('seep bc (2)', header=None)
        
        # Exit Face BC: starts at B5 (row 4, columns 1 and 2), continues down until empty x value
        exit_coords2 = []
        row = exit_start_row
        while row < seep_df2.shape[0]:
            try:
                x_val = seep_df2.iloc[row, exit_x_col]
                y_val = seep_df2.iloc[row, exit_y_col]
                
                # Stop at first empty x value
                if pd.isna(x_val):
                    break
                
                # If x is present, try to convert
                if pd.notna(x_val) and pd.notna(y_val):
                    exit_coords2.append((float(x_val), float(y_val)))
            except:
                break
            row += 1
        seepage_bc2["exit_face"] = exit_coords2
        
        # Specified Head BCs: same structure as first sheet
        col = 4  # Start with column E (index 4)
        while col < seep_df2.shape[1]:
            x_col = col
            y_col = col + 1
            head_col = col + 1  # Head value is in the Y column
            
            # Check if head value in row 3 is empty - stop reading if empty
            if seep_df2.shape[0] <= head_row:
                break
            head_val = seep_df2.iloc[head_row, head_col]
            if pd.isna(head_val):
                break  # Stop reading when head value is empty
            
            # Read XY coordinates starting from row 5, continue down until empty
            coords = []
            row = data_start_row
            while row < seep_df2.shape[0]:
                try:
                    x_val = seep_df2.iloc[row, x_col]
                    y_val = seep_df2.iloc[row, y_col]
                    
                    # Stop at first empty x value
                    if pd.isna(x_val):
                        break
                    
                    # If x is present, try to convert
                    if pd.notna(x_val) and pd.notna(y_val):
                        coords.append((float(x_val), float(y_val)))
                except:
                    break
                row += 1
            
            if coords:  # Only add if we have coordinates
                seepage_bc2["specified_heads"].append({"head": float(head_val), "coords": coords})
            
            # Move to next specified head BC (skip 3 columns: E->H, H->K, etc.)
            col += 3
    except (ValueError, KeyError):
        # If "seep bc (2)" sheet doesn't exist, just leave seepage_bc2 as empty
        pass

    # === VALIDATION ===
 
    circular = len(circles) > 0
    # Check if this is a seep-only analysis (has seep BCs but no slope stability surfaces)
    has_seepage_bc = (len(seepage_bc.get("specified_heads", [])) > 0 or 
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
    globals_data["piezo_line2"] = piezo_line2
    globals_data["circular"] = circular # True if circles are present
    globals_data["circles"] = circles
    globals_data["non_circ"] = non_circ
    globals_data["dloads"] = dloads
    globals_data["dloads2"] = dloads2
    globals_data["reinforce_lines"] = reinforce_lines
    globals_data["reinforcement_lines"] = reinforcement_lines
    globals_data["pile_lines"] = pile_lines
    globals_data["seepage_bc"] = seepage_bc
    globals_data["seepage_bc2"] = seepage_bc2
    globals_data["has_seepage_bc2"] = bool(seepage_bc2.get("specified_heads") or seepage_bc2.get("exit_face"))
    
    # Add mesh if available (used by both seep and fem workflows)
    globals_data["mesh"] = mesh

    # Add seep solution data if available
    if has_seep_materials:
        globals_data["seep_u"] = seep_u
        if seep_u2 is not None:
            globals_data["seep_u2"] = seep_u2

    return globals_data

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
# importer (see plans/plan_io.md).
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
    if isinstance(value, float):
        value = round(value, 10)
    open_match = re.match(r'(<c\s[^>]*?)(/?>)', cell_xml)
    if not open_match:
        return cell_xml
    open_tag_attrs = open_match.group(1)
    open_tag_attrs = re.sub(r'\s+t="[^"]*"', '', open_tag_attrs)
    if isinstance(value, str):
        return f'{open_tag_attrs} t="inlineStr"><is><t>{value}</t></is></c>'
    else:
        return f'{open_tag_attrs}><v>{value}</v></c>'


def _build_new_cell(ref, value):
    if isinstance(value, float):
        value = round(value, 10)
    if isinstance(value, str):
        return f'<c r="{ref}" t="inlineStr"><is><t>{value}</t></is></c>'
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