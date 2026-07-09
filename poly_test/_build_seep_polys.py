"""Build polygon-sheet versions of the four BYU CE544 seepage problems.

For each source (profile-based) input file, this:
  1. copies the v10 template,
  2. faithfully copies the main globals, the mat sheet, and the seep bc
     (and seep bc (2)) sheets from the source,
  3. converts the profile-line geometry to polygons via build_polygons() and
     writes them to the polygon sheet,
  4. forces a full recalc on open so the XLOOKUP material-name cells refresh.

Run from the repo root:  PYTHONPATH=. python3 poly_test/_build_seep_polys.py
"""
import os
import shutil

from xslope.fileio import (load_slope_data, write_cells_to_xlsx, cell_ref,
                           mat_header_cols, MAT_NUM_HEADERS, MAT_STR_HEADERS)
from xslope.mesh import build_polygons

# The surgical xlsx writer now lives in xslope.fileio
# (write_cells_to_xlsx, cell_ref) and is imported above.


# --- Populate helpers --------------------------------------------------------

import math


def _is_num(v):
    """True only for a real, non-zero, non-NaN number worth writing."""
    if v is None or isinstance(v, str):
        return False
    try:
        f = float(v)
    except (TypeError, ValueError):
        return False
    return not math.isnan(f) and f != 0


def _clean_str(v):
    """Return a real string to write, or None to skip (drops pandas 'nan'/blank)."""
    if v is None:
        return None
    s = str(v).strip()
    if s == '' or s.lower() == 'nan':
        return None
    return s


def build_polygon_version(src, dst, template="docs/inputs/input_template.xlsx"):
    shutil.copy(template, dst)
    data = load_slope_data(src)
    polys = build_polygons(data)
    updates = {}

    # --- main: global options ---
    updates['main'] = {
        cell_ref(8, 4): float(data['gamma_water']),
        cell_ref(9, 4): float(data['tcrack_depth']),
        cell_ref(10, 4): float(data['tcrack_water']),
        cell_ref(11, 4): float(data['k_seismic']),
    }

    # --- mat: copy every populated field from the source materials ---
    # Header row and columns are located BY NAME in the destination file. The previous
    # hardcoded v10 map (kr0=21, h0=22, E=23, nu=24) silently shifted when v11 inserted
    # 'unsat' at column 21.
    mat_hdr, mat_cols = mat_header_cols(dst)
    updates['mat'] = {}
    for i, mat in enumerate(data['materials']):
        row = mat_hdr + 1 + i
        for key, header in [('name', 'name')] + MAT_STR_HEADERS:
            col = mat_cols.get(header)
            s = _clean_str(mat.get(key))
            if col is not None and s is not None:
                updates['mat'][cell_ref(row, col)] = s
        for key, header in MAT_NUM_HEADERS:
            col = mat_cols.get(header)
            v = mat.get(key)
            if col is not None and _is_num(v):
                updates['mat'][cell_ref(row, col)] = float(v)

    # --- polygon sheet ---
    # Block p (1-based): x_col = 1 + 3*(p-1), y_col = x_col + 1.
    # Row 5: Mat ID at (5, y_col).  Row 6: XLOOKUP -> DO NOT WRITE.  Coords row 8+.
    # build_polygons returns 0-based mat_id (sheet = +1).
    updates['polygon'] = {}
    max_pts = 0
    for p, poly in enumerate(polys, start=1):
        x_col = 1 + 3 * (p - 1)
        y_col = x_col + 1
        updates['polygon'][cell_ref(5, y_col)] = poly['mat_id'] + 1
        for i, (x, y) in enumerate(poly['coords']):
            updates['polygon'][cell_ref(8 + i, x_col)] = float(x)
            updates['polygon'][cell_ref(8 + i, y_col)] = float(y)
        max_pts = max(max_pts, len(poly['coords']))
    # Clear leftover template defaults in unused blocks (the template ships with
    # Mat IDs 1..15 pre-filled in row 5 for every block — blank the extras).
    for p in range(len(polys) + 1, 16):
        y_col = 2 + 3 * (p - 1)
        updates['polygon'][cell_ref(5, y_col)] = ''

    # --- seep bc (and seep bc (2)) ---
    def write_bc(sheet_name, bc):
        cells = {}
        for i, (x, y) in enumerate(bc.get('exit_face', [])):
            cells[cell_ref(5 + i, 2)] = float(x)  # B
            cells[cell_ref(5 + i, 3)] = float(y)  # C
        for h, head in enumerate(bc.get('specified_heads', []), start=1):
            x_col = 5 + 3 * (h - 1)   # E, H, K, ...
            y_col = x_col + 1
            cells[cell_ref(3, y_col)] = float(head['head'])  # head value in row 3
            for i, (x, y) in enumerate(head['coords']):
                cells[cell_ref(5 + i, x_col)] = float(x)
                cells[cell_ref(5 + i, y_col)] = float(y)
        if cells:
            updates[sheet_name] = cells

    write_bc('seep bc', data['seepage_bc'])
    if data.get('has_seepage_bc2'):
        write_bc('seep bc (2)', data['seepage_bc2'])

    updates = {k: v for k, v in updates.items() if v}
    write_cells_to_xlsx(dst, updates)
    return polys, data


PROBLEMS = [
    ("SEA TRENCH",  "docs/seep/files/xslope_sea_trench.xlsx", "poly_test/xslope_sea_trench_poly.xlsx"),
    ("LEVEE",       "docs/seep/files/xslope_levee_full.xlsx", "poly_test/xslope_levee_poly.xlsx"),
    ("JOHNSON RES", "docs/seep/files/xslope_johnson_res.xlsx", "poly_test/xslope_johnson_res_poly.xlsx"),
    ("EARTH DAM",   "docs/seep/files/xslope_earth_dam1.xlsx", "poly_test/xslope_earth_dam_poly.xlsx"),
]

if __name__ == "__main__":
    for name, src, dst in PROBLEMS:
        polys, data = build_polygon_version(src, dst)
        nheads = len(data['seepage_bc'].get('specified_heads', []))
        nexit = len(data['seepage_bc'].get('exit_face', []))
        bc2 = " +bc2" if data.get('has_seepage_bc2') else ""
        print(f"{name:12s} -> {dst}")
        print(f"             {len(data['materials'])} materials, {len(polys)} polygons, "
              f"{nheads} specified head(s), {nexit}-pt exit face{bc2}")
