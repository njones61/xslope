"""Build poly_test/xslope_lost_lake_poly.xlsx by copying the current template and
populating it with the Lost Lake inputs, using the polygon sheet (derived from
the profile-line geometry via build_polygons) instead of the profile sheet."""
import os
import shutil

from xslope.fileio import load_slope_data, write_cells_to_xlsx, cell_ref, mat_header_cols
from xslope.mesh import build_polygons

# The surgical xlsx writer now lives in xslope.fileio
# (write_cells_to_xlsx, cell_ref) and is imported above.


# --- Copy template ---
src = "docs/inputs/input_template.xlsx"
dst = "poly_test/xslope_lost_lake_poly.xlsx"
shutil.copy(src, dst)

# --- Load source data + build polygons from profile lines ---
data = load_slope_data("docs/inputs/seep/xslope_lost_lake.xlsx")
polys = build_polygons(data)
mat_names = [m['name'] for m in data['materials']]

updates = {}

# --- main: global options (mirror source) ---
updates['main'] = {
    cell_ref(8, 4): 62.4,   # unit weight of water
    cell_ref(9, 4): 0,      # tension crack depth
    cell_ref(10, 4): 0,     # depth of water in crack
    cell_ref(11, 4): 0,     # seismic coefficient
}

# --- mat: names + seepage conductivity params (seepage-only problem) ---
# Source rows 9..16 -> materials 1..8. Columns: B name, R k1, S k2, T alpha,
# U kr0, V h0  (cols 2,18,19,20,21,22). Strength cols left blank as in source.
seep_props = [
    # name,           k1,     k2,     alpha, kr0,    h0
    ("Rip Rap",       100000, 100000, 0,     0.0001, -1),
    ("Shell",         250,    50,     0,     0.0001, -1),
    ("Clay Core",     0.1,    0.1,    0,     0.0001, -1),
    ("Drain",         10000,  10000,  0,     0.0001, -1),
    ("Filter",        1000,   1000,   0,     0.0001, -1),
    ("Glacial Till",  4000,   2000,   0,     0.0001, -1),
    ("Grouted Bedrock", 250,  250,    0,     0.0001, -1),
    ("Bedrock",       2000,   1000,   0,     0.0001, -1),
]
# Locate the mat header row and its columns BY NAME in the destination file. Writing by
# hardcoded column number silently shifts when a column is inserted -- v11 added 'unsat'
# at column 21, which would have sent kr0 into it and h0 into kr0.
mat_hdr, mat_cols = mat_header_cols(dst)
updates['mat'] = {}
for i, (name, k1, k2, alpha, kr0, h0) in enumerate(seep_props):
    row = mat_hdr + 1 + i
    for header, val in (('name', name), ('k1', k1), ('k2', k2),
                        ('alpha', alpha), ('kr0', kr0), ('h0', h0)):
        col = mat_cols.get(header)
        if col is not None:
            updates['mat'][cell_ref(row, col)] = val

# --- polygon sheet ---
# Block layout: polygon p (1-based) -> x_col = 1 + 3*(p-1), y_col = x_col + 1.
# Row 5: Mat ID value at (5, y_col).  Row 6: XLOOKUP formula -> DO NOT WRITE.
# Coords start at row 8.  build_polygons returns 0-based mat_id (sheet = +1).
updates['polygon'] = {}
for p, poly in enumerate(polys, start=1):
    x_col = 1 + 3 * (p - 1)
    y_col = x_col + 1
    updates['polygon'][cell_ref(5, y_col)] = poly['mat_id'] + 1
    for i, (x, y) in enumerate(poly['coords']):
        updates['polygon'][cell_ref(8 + i, x_col)] = float(x)
        updates['polygon'][cell_ref(8 + i, y_col)] = float(y)

# --- seep bc: exit face + specified head #1 (mirror source) ---
updates['seep bc'] = {}
exit_face = [(260, 230), (316, 202), (320, 200), (324, 198), (475, 198)]
for i, (x, y) in enumerate(exit_face):
    updates['seep bc'][cell_ref(5 + i, 2)] = float(x)  # B
    updates['seep bc'][cell_ref(5 + i, 3)] = float(y)  # C
# Specified head #1: head value F3 = 225, XY in E/F starting row 5
head1 = [(0, 200), (150, 200), (225, 225)]
updates['seep bc'][cell_ref(3, 6)] = 225.0  # F3
for i, (x, y) in enumerate(head1):
    updates['seep bc'][cell_ref(5 + i, 5)] = float(x)  # E
    updates['seep bc'][cell_ref(5 + i, 6)] = float(y)  # F

# --- write ---
updates = {k: v for k, v in updates.items() if v}
write_cells_to_xlsx(dst, updates)
print(f"Wrote {dst}")
print(f"Polygons written: {len(polys)}")
