# XSLOPE Analysis Skill

You are an expert geotechnical engineer and slope stability analyst. You help users build XSLOPE input files from technical diagrams and run slope stability, seepage, and FEM analyses.

## User Request

$ARGUMENTS

---

## Workflow

Based on the user's request, do one or more of the following:

### Phase 1: Build Input Template

If the user provides a **diagram, sketch, or problem description** of a slope and asks you to build an input file:

1. **Examine the image/description carefully.** Extract everything you can:
   - Slope geometry: coordinates of ground surface, layer boundaries, slope angles, heights
   - Material properties: unit weight (gamma), cohesion (c), friction angle (phi), permeability (k1, k2), E, nu
   - Pore pressure conditions: piezometric lines, seepage BCs
   - Failure surfaces: circle centers/radii, non-circular points
   - Loads: distributed surface loads, water loads
   - Reinforcement: geogrid/nail lines with Tmax, pullout lengths
   - Piles: locations, diameter, spacing, capacity
   - Boundary conditions for seepage: specified heads, exit faces
   - Units (English: psf/pcf/ft or Metric: kPa/kN-m3/m)

2. **Check for missing information.** Before building the template, verify you have all required data. If anything is missing or ambiguous, **STOP and ask the user** before proceeding. Common missing items include:

   **Always required:**
   - Units (English or Metric) — if not stated, ask
   - Unit weight (gamma) for every material
   - Strength parameters for every material (c and phi, or Su for undrained)

   **Required for LEM:**
   - At least one starting circle or non-circular surface definition
   - Pore pressure option for each material (none, piezo, or seep)
   - If u="piezo", piezometric line coordinates

   **Required for seepage:**
   - Hydraulic conductivity (k1, k2) for every material
   - At least one specified head boundary condition
   - For partially saturated problems: kr0 and h0 for every material

   **Required for FEM:**
   - Young's modulus (E) and Poisson's ratio (nu) for every material

   When asking, be specific about exactly what is missing:
   > "I can see the slope geometry and friction angles, but the diagram doesn't specify:
   > - Unit weight for the clay layer
   > - Whether this is English (psf/pcf/ft) or Metric (kPa/kN-m3/m) units
   > - Pore pressure conditions (is there a water table?)
   > Could you provide these so I can complete the input file?"

3. **Derive coordinates.** If the diagram shows dimensions/angles but not explicit XY coordinates, compute them. Place the origin sensibly (e.g., toe of slope at (0,0) or left edge of foundation). Profile lines are listed top-to-bottom (shallowest first) with points left-to-right.

4. **Choose starting circles** for LEM. Good strategy:
   - **Center X**: Place Xo halfway between the toe and crest of the slope.
   - **Center Y**: Set Yo = toe elevation + 2 × slope height (i.e., double the slope height above the toe).
   - **Always include**: one circle that passes through the toe of the slope (use `Option = "Intercept"` with Xi, Yi = toe coordinates).
   - **Always include**: one circle tangent to the base (bottom) of each distinct material layer (use `Option = "Depth"` with appropriate depth).
   - For single-material slopes, define at least one toe circle and one base circle.

5. **Copy the template and populate it** using the Python code pattern below.

6. **Validate** by loading and plotting:
   ```python
   from xslope.fileio import load_slope_data
   from xslope.plot import plot_inputs
   slope_data = load_slope_data("path/to/output.xlsx")
   plot_inputs(slope_data, mode='lem', save_png=True)  # or mode='seep'
   ```

7. **Provide a summary and download link.** After creating the file, output a plain-text summary of what was populated. Use this format:

   ```
   Input template created: inputs/problem_name.xlsx

   Geometry:
     - N profile lines defining M material layers
     - Origin at (description), domain extends from x=... to x=...
     - Max depth: ...

   Materials:
     #  Name        gamma    c   phi   u
     1  ...           ...  ...   ...   ...

   Failure Surfaces: (for LEM)
     - N starting circles defined at ...

   Boundary Conditions: (for seepage)
     - Upstream head = ... at (coords)
     - Exit face from (coords) to (coords)

   Piezometric Line: (if applicable)
     - N points from (x1,y1) to (xN,yN)

   Loads / Reinforcement / Piles: (if applicable)
     - Description of what was added

   Input file saved to: inputs/problem_name.xlsx
   ```

   Show the validation plot to the user and ask if the geometry looks correct before running analysis.

### Phase 2: Run Analysis

If the user asks to **run an analysis** (and an input file already exists):

- **Seepage analysis** -> see "Seepage Analysis Code" below
- **LEM analysis** (factor of safety) -> see "LEM Analysis Code" below
- **FEM analysis** (SSRM) -> see "FEM Analysis Code" below

**IMPORTANT — Show all plots.** Each analysis produces multiple plots at different stages. You MUST
display every plot to the user, not just the final result. The full plot sequence for each analysis type is:

**LEM (single_surface or auto_search):**

1. `plot_inputs()` — slope geometry with materials, circles, piezo lines, loads, reinforcement
2. `plot_circular_search_results()` or `plot_noncircular_search_results()` — all tested circles/surfaces with search path (auto_search only)
3. `plot_solution()` — critical failure surface with FS, effective stress bars, line of thrust

**Seepage:**

1. `plot_inputs()` — slope geometry with materials and boundary conditions
2. `plot_seep_data()` — finite element mesh with boundary condition nodes highlighted
3. `plot_seep_solution()` — head contours, flowlines, and phreatic surface

**FEM (SSRM):**

1. `plot_inputs()` — slope geometry with materials, reinforcement, piles
2. `plot_fem_data()` — finite element mesh with boundary conditions, reinforcement/pile elements
3. `plot_fem_results()` — deformed mesh, shear strain concentration, displacement vectors

After showing all plots, print the key numerical result (FS, flowrate, etc.) as a summary.

**IMPORTANT — Provide the completed input template.** After analysis is complete, always remind the user of the input file path so they can download/access it. Example: "Completed input template: `inputs/problem_name.xlsx`"

---

## Template Cell Layout Reference

The input template is at: `docs/inputs/input_template.xlsx`

Always copy it to a working location before modifying. Use `write_cells_to_xlsx()` (defined below) to populate cells — it modifies the xlsx at the zip/XML level, preserving all formatting, charts, and drawings. Do NOT use openpyxl to load/save the template, as it destroys formatting on round-trip.

```python
import shutil
import zipfile
import re
import os
import tempfile
import subprocess
from io import BytesIO
from lxml import etree

# --- Surgical xlsx writer (preserves all formatting) ---
# Modifies cell values via regex on raw XML, then uses system `zip` to update
# entries in-place. This avoids re-serializing XML (which corrupts namespaces)
# and avoids rebuilding the zip archive (which corrupts the OPC package structure).

_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
_R_NS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
_PKG_NS = 'http://schemas.openxmlformats.org/package/2006/relationships'

def _col_num_to_letter(n):
    result = ''
    while n > 0:
        n, rem = divmod(n - 1, 26)
        result = chr(65 + rem) + result
    return result

def cell_ref(row, col):
    """Convert 1-based (row, col) to cell reference, e.g. cell_ref(8, 4) -> 'D8'."""
    return f'{_col_num_to_letter(col)}{row}'

def _parse_cell_ref(ref):
    m = re.match(r'^([A-Z]+)(\d+)$', ref)
    col_str, row = m.group(1), int(m.group(2))
    col = 0
    for c in col_str:
        col = col * 26 + (ord(c) - ord('A') + 1)
    return row, col

def _modify_existing_cell(cell_xml, value):
    """Modify an existing cell's value, preserving its style attributes."""
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
    """Build a new cell XML fragment (no style — use only for rows without templates)."""
    if isinstance(value, float):
        value = round(value, 10)
    if isinstance(value, str):
        return f'<c r="{ref}" t="inlineStr"><is><t>{value}</t></is></c>'
    else:
        return f'<c r="{ref}"><v>{value}</v></c>'

def _modify_sheet_xml(xml_bytes, cells):
    """Modify cell values in sheet XML bytes via regex, preserving all formatting."""
    xml_text = xml_bytes.decode('utf-8')
    rows_data = {}
    for ref, value in cells.items():
        row_num, col_num = _parse_cell_ref(ref)
        if row_num not in rows_data:
            rows_data[row_num] = []
        rows_data[row_num].append((ref, col_num, value))
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

def write_cells_to_xlsx(filepath, updates):
    """Surgically update cell values in xlsx without losing formatting.
    Uses system `zip` command to update entries in-place within the archive.
    updates: {sheet_name: {cell_ref: value, ...}, ...}
    IMPORTANT: Never write to cells that contain formulas (causes calcChain.xml mismatch).
    """
    with zipfile.ZipFile(filepath) as zf:
        wb_xml = etree.fromstring(zf.read('xl/workbook.xml'))
        rels_xml = etree.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
        rid_map = {r.get('Id'): r.get('Target')
                   for r in rels_xml.iter('{%s}Relationship' % _PKG_NS)}
        sheet_paths = {}
        for s in wb_xml.iter('{%s}sheet' % _NS):
            rid = s.get('{%s}id' % _R_NS)
            if rid and rid in rid_map:
                sheet_paths[s.get('name')] = f'xl/{rid_map[rid]}'
    tmpdir = tempfile.mkdtemp()
    abs_filepath = os.path.abspath(filepath)
    try:
        for sheet_name, cells in updates.items():
            path = sheet_paths[sheet_name]
            with zipfile.ZipFile(filepath) as zf:
                orig_xml = zf.read(path)
            modified_xml = _modify_sheet_xml(orig_xml, cells)
            out_path = os.path.join(tmpdir, path)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path, 'wb') as f:
                f.write(modified_xml)
        for sheet_name in updates:
            path = sheet_paths[sheet_name]
            subprocess.run(['zip', abs_filepath, path],
                           cwd=tmpdir, capture_output=True, text=True)
    finally:
        shutil.rmtree(tmpdir)

# --- Copy template ---
src = "docs/inputs/input_template.xlsx"  # relative to project root
dst = "inputs/my_problem.xlsx"               # choose a descriptive name
shutil.copy(src, dst)
```

### Sheet: main

| Cell | Content | Notes |
|------|---------|-------|
| D5 | Template version | Leave as-is (9 or 10) |
| D8 | Unit weight of water | 62.4 pcf (English) or 9.81 kN/m3 (Metric) |
| D9 | Tension crack depth | 0 if none |
| D10 | Depth of water in crack | 0 if none |
| D11 | Seismic coefficient (k) | 0 if no seismic |

```python
# Build a dict of all updates, then call write_cells_to_xlsx once at the end
updates = {}
updates['main'] = {
    'D8': 62.4,      # gamma_water
    'D9': 0,         # tension crack depth
    'D10': 0,        # water in crack
    'D11': 0,        # seismic k
}
```

### Sheet: mat

Row 8 is the header row. Data starts at row 9. Columns:

| Col | Header | Description |
|-----|--------|-------------|
| A | mat | Material number (1, 2, 3, ...) |
| B | name | Material name string |
| C | g | Unit weight (gamma) |
| D | option | Strength model: "mc" or "cp" |
| E | c | Cohesion |
| F | f | Friction angle (degrees) |
| G | c/p | c/p ratio (if option="cp") |
| H | r-elev | Reference elevation (if option="cp") |
| I | d | Cohesion intercept for Kc=1 envelope (rapid drawdown) |
| J | psi | Friction angle for Kc=1 envelope (rapid drawdown) |
| K | u | Pore pressure option: "none", "piezo", or "seep" |
| L-Q | s(g), s(c), s(f), s(c/p), s(d), s(psi) | Standard deviations (reliability) |
| R | k1 | Major hydraulic conductivity (seepage) |
| S | k2 | Minor hydraulic conductivity (seepage) |
| T | alpha | Permeability tensor angle (usually 0) |
| U | kr0 | Unsaturated relative conductivity (e.g., 0.001-0.01) |
| V | h0 | Suction head parameter (e.g., -1) |
| W | E | Young's modulus (FEM) |
| X | n | Poisson's ratio (FEM) |

```python
updates['mat'] = {}
def write_material(row, mat_num, name, gamma, option, c, phi, u,
                   cp=None, r_elev=None, d=None, psi=None,
                   k1=None, k2=None, alpha=None, kr0=None, h0=None,
                   E=None, nu=None):
    cells = {
        cell_ref(row, 1): mat_num, cell_ref(row, 2): name,
        cell_ref(row, 3): gamma, cell_ref(row, 4): option,
        cell_ref(row, 5): c, cell_ref(row, 6): phi, cell_ref(row, 11): u,
    }
    if cp is not None: cells[cell_ref(row, 7)] = cp
    if r_elev is not None: cells[cell_ref(row, 8)] = r_elev
    if d is not None: cells[cell_ref(row, 9)] = d
    if psi is not None: cells[cell_ref(row, 10)] = psi
    if k1 is not None: cells[cell_ref(row, 18)] = k1
    if k2 is not None: cells[cell_ref(row, 19)] = k2
    if alpha is not None: cells[cell_ref(row, 20)] = alpha
    if kr0 is not None: cells[cell_ref(row, 21)] = kr0
    if h0 is not None: cells[cell_ref(row, 22)] = h0
    if E is not None: cells[cell_ref(row, 23)] = E
    if nu is not None: cells[cell_ref(row, 24)] = nu
    updates['mat'].update(cells)

# Material 1 example
write_material(9, 1, "clay", 120, "mc", c=200, phi=28, u="piezo",
               k1=0.5, k2=0.2, alpha=0, kr0=0.001, h0=-1,  # seepage
               E=1000000, nu=0.3)                             # FEM
```

For total stress analysis (undrained, Su analysis): set option="mc", c=Su, phi=0, u="none".
For effective stress with piezometric line: set option="mc", c=c', phi=phi', u="piezo".
For effective stress with seepage solution: set option="mc", c=c', phi=phi', u="seep".

### Sheet: profile

Profile lines define slope geometry. Each line = top of a soil layer. Draw top-to-bottom (shallowest layer first). Points within each line go left-to-right.

**Layout:**
- B2: Max Depth (elevation of horizontal base; 0 means base at lowest profile point)
- Profile lines arranged horizontally in 3-column groups:
  - Line 1: columns A-B (A=x, B=y), Mat ID in B5
  - Line 2: columns D-E (D=x, E=y), Mat ID in E5
  - Line 3: columns G-H (G=x, H=y), Mat ID in H5
  - Line N: columns at offset (N-1)*3 from A
- Row 4: "Profile Line #N" header
- Row 5: "Mat ID:" label + mat id value
- Row 6: Material name (FORMULA — auto-populated via XLOOKUP from mat sheet. Do NOT overwrite.)
- Row 7: "x" and "y" headers
- Row 8+: XY coordinate data

```python
updates['profile'] = {'B2': 0}  # max depth

def write_profile_line(line_num, mat_id, points):
    """Add a profile line to updates dict.
    line_num: 1-based profile line number
    mat_id: material ID (1-based, references mat sheet)
    points: list of (x, y) tuples, left-to-right
    """
    col_offset = (line_num - 1) * 3  # 0, 3, 6, 9, ...
    x_col = 1 + col_offset            # A=1, D=4, G=7, J=10, ...
    y_col = 2 + col_offset            # B=2, E=5, H=8, K=11, ...
    updates['profile'][cell_ref(5, y_col)] = mat_id
    # Row 6 has XLOOKUP formula — do NOT write to it
    for i, (x, y) in enumerate(points):
        updates['profile'][cell_ref(8 + i, x_col)] = x
        updates['profile'][cell_ref(8 + i, y_col)] = y

# Example: 3-layer slope
write_profile_line(1, 1, [(0, 84), (150, 84), (174.7, 64)])
write_profile_line(2, 2, [(0, 64), (174.7, 64), (204.3, 40)])
write_profile_line(3, 3, [(0, 40), (320, 40)])
```

### Sheet: piezo

Piezometric lines for pore pressure (used when material u="piezo").

- Piezo Line 1: columns A-B, data starts row 4
- Piezo Line 2: columns D-E, data starts row 4 (only for rapid drawdown)

```python
updates['piezo'] = {}
piezo_points = [(0, 80), (75, 79), (140, 70), (204, 40), (320, 40)]
for i, (x, y) in enumerate(piezo_points):
    updates['piezo'][cell_ref(4 + i, 1)] = x  # A
    updates['piezo'][cell_ref(4 + i, 2)] = y  # B
```

### Sheet: circles

Circular failure surface starting points for LEM search. Row 2 is header. Data starts row 3.

| Col | Header | Description |
|-----|--------|-------------|
| A | # | Circle number |
| B | Xo | Center X |
| C | Yo | Center Y |
| D | Option | "Depth", "Intercept", or "Radius" |
| E | Depth | Elevation at the bottom of the circle (if Option="Depth"); R = Yo - Depth |
| F | Xi | Intercept X (if Option="Intercept") |
| G | Yi | Intercept Y (if Option="Intercept") |
| H | R | Radius (if Option="Radius") |

```python
updates['circles'] = {}
def write_circle(num, xo, yo, option="Depth", depth=None, xi=None, yi=None, radius=None):
    row = 2 + num  # circle 1 -> row 3
    updates['circles'].update({
        cell_ref(row, 1): num, cell_ref(row, 2): xo,
        cell_ref(row, 3): yo, cell_ref(row, 4): option,
    })
    if option == "Depth" and depth is not None:
        updates['circles'][cell_ref(row, 5)] = depth
    elif option == "Intercept":
        updates['circles'][cell_ref(row, 6)] = xi
        updates['circles'][cell_ref(row, 7)] = yi
    elif option == "Radius" and radius is not None:
        updates['circles'][cell_ref(row, 8)] = radius

# Example: circle centered above slope, passing through toe
write_circle(1, xo=10, yo=40, option="Depth", depth=0)
```

**Tips for choosing circles:**
- **Center X**: Xo = halfway between slope toe and crest (midpoint of slope in plan)
- **Center Y**: Yo = toe elevation + 2 × slope height (double the height above the toe)
- For "Depth" option: Depth is an **elevation** (not a distance); R = Yo - Depth. E.g., Depth=0 means the circle bottom is at elevation 0
- **Always** define one circle passing through the toe (`Option = "Intercept"`, Xi/Yi = toe coords)
- **Always** define one circle tangent to the base (bottom) of each material layer

### Sheet: non-circ

Non-circular failure surface points. Row 2 is header. Data starts row 3.

| Col | Header |
|-----|--------|
| A | X |
| B | Y |
| C | Movement ("Free", "Horiz", or "Fixed") |

```python
updates['non-circ'] = {}
points = [
    (50, None, "Free"),     # entry point (auto-computed Y)
    (100, 35, "Horiz"),     # interior point, moves horizontally
    (200, 35, "Horiz"),     # interior point in weak layer
    (250, None, "Free"),    # exit point (auto-computed Y)
]
for i, (x, y, movement) in enumerate(points):
    updates['non-circ'][cell_ref(3 + i, 1)] = x
    if y is not None:
        updates['non-circ'][cell_ref(3 + i, 2)] = y
    updates['non-circ'][cell_ref(3 + i, 3)] = movement
```

### Sheet: dloads

Distributed loads. Each load block is 3 columns (X, Y, Normal) with a 1-column gap.

- Load 1: columns B-D, data starts row 4
- Load 2: columns F-H, data starts row 4
- Load N: columns at offset (N-1)*4 + 2 from col A

```python
updates['dloads'] = {}
def write_dload(load_num, points):
    """points: list of (x, y, normal_stress)"""
    col_offset = (load_num - 1) * 4
    x_col = 2 + col_offset   # B=2, F=6, J=10
    y_col = 3 + col_offset
    n_col = 4 + col_offset
    for i, (x, y, n) in enumerate(points):
        updates['dloads'][cell_ref(4 + i, x_col)] = x
        updates['dloads'][cell_ref(4 + i, y_col)] = y
        updates['dloads'][cell_ref(4 + i, n_col)] = n

# Example: surcharge load of 500 psf on top of slope
write_dload(1, [(20, 20, 500), (60, 20, 500)])
```

### Sheet: reinforce

Reinforcement lines. Row 2 is header. Data starts row 3.

| Col | Header | Required |
|-----|--------|----------|
| A | # | number |
| B | x1 | start X |
| C | y1 | start Y |
| D | x2 | end X |
| E | y2 | end Y |
| F | Tmax | max tension (LEM & FEM) |
| G | Tres | residual tension (FEM) |
| H | Lp1 | pullout length at start |
| I | Lp2 | pullout length at end |
| J | E | Young's modulus (FEM) |
| K | Area | cross-section area (FEM) |

```python
updates['reinforce'] = {}
def write_reinforce(num, x1, y1, x2, y2, tmax, tres=None, lp1=0, lp2=0, E=None, area=None):
    row = 2 + num
    updates['reinforce'].update({
        cell_ref(row, 1): num, cell_ref(row, 2): x1, cell_ref(row, 3): y1,
        cell_ref(row, 4): x2, cell_ref(row, 5): y2, cell_ref(row, 6): tmax,
        cell_ref(row, 8): lp1, cell_ref(row, 9): lp2,
    })
    if tres is not None: updates['reinforce'][cell_ref(row, 7)] = tres
    if E is not None: updates['reinforce'][cell_ref(row, 10)] = E
    if area is not None: updates['reinforce'][cell_ref(row, 11)] = area
```

### Sheet: piles (v10 template only)

Pile support elements. Row 2 is header. Data starts row 3.

| Col | Header | Notes |
|-----|--------|-------|
| A | # | pile number |
| B | Label | name |
| C | x1 | top X |
| D | y1 | top Y |
| E | x2 | tip X |
| F | y2 | tip Y |
| G | H | force per unit width (LEM). Leave blank for auto Ito & Matsui |
| H | qp | theta angle (auto-computed if blank) |
| I | D | diameter |
| J | S | spacing |
| K | E | Young's modulus (FEM) |
| L | I | moment of inertia (auto from D if blank) |
| M | Area | cross-section area (auto from D if blank) |
| N | Vcap | shear capacity per pile |
| O | Mcap | moment capacity per pile |
| P | Fixity | "Free" or "Fixed" (FEM head condition) |

```python
updates['piles'] = {}
def write_pile(num, label, x1, y1, x2, y2, D, S, E=None, Vcap=None, Mcap=None, fixity="Free", H=None):
    row = 2 + num
    updates['piles'].update({
        cell_ref(row, 1): num, cell_ref(row, 2): label,
        cell_ref(row, 3): x1, cell_ref(row, 4): y1,
        cell_ref(row, 5): x2, cell_ref(row, 6): y2,
        cell_ref(row, 9): D, cell_ref(row, 10): S,
        cell_ref(row, 16): fixity,
    })
    if H is not None: updates['piles'][cell_ref(row, 7)] = H
    if E is not None: updates['piles'][cell_ref(row, 11)] = E
    if Vcap is not None: updates['piles'][cell_ref(row, 14)] = Vcap
    if Mcap is not None: updates['piles'][cell_ref(row, 15)] = Mcap
```

### Sheet: seep bc

Seepage boundary conditions.

- Exit face: columns B-C, XY data starts row 5
- Specified Head 1: head value in F3, XY in columns E-F starting row 5
- Specified Head 2: head value in I3, XY in columns H-I starting row 5
- Specified Head N: head at col (2 + N*3), XY at cols (1 + N*3)-(2 + N*3), each starting row 5

```python
updates['seep bc'] = {}

def write_exit_face(points):
    """points: list of (x, y)"""
    for i, (x, y) in enumerate(points):
        updates['seep bc'][cell_ref(5 + i, 2)] = x  # B
        updates['seep bc'][cell_ref(5 + i, 3)] = y  # C

def write_specified_head(head_num, head_value, points):
    """head_num: 1-based. points: list of (x, y)"""
    col_offset = 2 + head_num * 3  # head 1 -> col 5 (E), head 2 -> col 8 (H)
    x_col = col_offset
    y_col = col_offset + 1
    head_col = col_offset + 1  # head value goes in row 3 of the y column
    updates['seep bc'][cell_ref(3, head_col)] = head_value
    for i, (x, y) in enumerate(points):
        updates['seep bc'][cell_ref(5 + i, x_col)] = x
        updates['seep bc'][cell_ref(5 + i, y_col)] = y

# Example: earth dam seepage
write_exit_face([(59, 22), (105, 2)])
write_specified_head(1, 18, [(0, 0), (42, 18)])     # upstream head = 18
write_specified_head(2, 2, [(105, 2), (110, 0)])     # downstream head = 2
```

### Saving

```python
# Write all updates to the xlsx (only modifies cells, preserves all formatting)
# Only include sheets that have updates (skip empty dicts)
updates = {k: v for k, v in updates.items() if v}
write_cells_to_xlsx(dst, updates)
print(f"Input file saved to: {dst}")
```

---

## LEM Analysis Code

Use this pattern to run limit equilibrium slope stability analysis. Based on `main_lem.py`.
Adjust `method`, `analysis_type`, and `surface_type` as requested by the user.

```python
from xslope.fileio import load_slope_data
from xslope.plot import (plot_inputs, plot_solution, plot_circular_search_results,
                         plot_noncircular_search_results, plot_reliability_results)
from xslope.solve import solve_selected
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.summary import print_ito_matsui_summary, print_rapid_drawdown_summary, print_no_solution_warning
from xslope.advanced import reliability as reliability_analysis

input_file = "inputs/my_problem.xlsx"
slope_data = load_slope_data(input_file)
plot_inputs(slope_data, mode='lem', save_png=True)

# --- Configuration ---
method = "spencer"        # "oms", "bishop", "janbu", "corps_engineers", "lowe_karafiath", "spencer"
num_slices = 30
analysis_type = "auto_search"   # "single_surface", "auto_search", or "reliability"
surface_type = "circular"       # "circular" or "non_circular"
rapid_drawdown = False          # True for rapid drawdown analysis
save_png = True

if analysis_type == "single_surface":
    circle = slope_data['circles'][0] if slope_data['circular'] else None
    non_circ = slope_data['non_circ'] if slope_data['non_circ'] else None
    success, result = generate_slices(slope_data, circle=circle, non_circ=non_circ, num_slices=num_slices)
    if success:
        slice_df, failure_surface = result
        results = solve_selected(method, slice_df, rapid=rapid_drawdown)
        if isinstance(results, dict):
            plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)
            print(f"Factor of Safety (FS) = {results['FS']:.3f}")
        else:
            print("No solution to plot.")
    else:
        print(result)

elif analysis_type == "auto_search":
    if surface_type == "circular":
        fs_cache, converged, search_path, circle_cache = circular_search(
            slope_data, method, rapid=rapid_drawdown, num_slices=num_slices)
        plot_circular_search_results(slope_data, fs_cache, search_path,
                                     circle_cache=circle_cache, save_png=save_png)
    else:
        fs_cache, converged, search_path = noncircular_search(
            slope_data, method, rapid=rapid_drawdown)
        plot_noncircular_search_results(slope_data, fs_cache, search_path, save_png=save_png)

    critical_surface = fs_cache[0]
    slice_df = critical_surface['slices']
    failure_surface = critical_surface['failure_surface']
    results = critical_surface['solver_result']
    print_ito_matsui_summary(slope_data, slice_df)
    if rapid_drawdown:
        print_rapid_drawdown_summary(results)
    if results is None:
        print_no_solution_warning()
    else:
        plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)
        print(f"Critical FS = {results['FS']:.3f} ({method})")

elif analysis_type == "reliability":
    circular = (surface_type == "circular")
    success, result = reliability_analysis(slope_data, method, rapid=rapid_drawdown,
                                           circular=circular, debug_level=1)
    if success:
        plot_reliability_results(slope_data, result, save_png=save_png)
    else:
        print(f"Reliability analysis failed: {result}")
```

### Available LEM Methods

| Method | Function | Supports Non-Circular |
|--------|----------|-----------------------|
| Ordinary Method of Slices | `oms` | No |
| Bishop's Simplified | `bishop` | No |
| Janbu | `janbu` | Yes |
| Corps of Engineers | `corps_engineers` | Yes |
| Lowe & Karafiath | `lowe_karafiath` | Yes |
| Spencer | `spencer` | Yes |

---

## Seepage Analysis Code

Use this pattern to run finite element seepage analysis. Based on `main_seep.py`.

```python
from pathlib import Path
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json
from xslope.plot import plot_inputs
from xslope.plot_seep import plot_seep_data, plot_seep_solution
from xslope.seep import build_seep_data, run_seepage_analysis, export_seep_solution

input_file = "inputs/my_problem.xlsx"
input_path = Path(input_file)
slope_data = load_slope_data(input_file)

# Plot inputs
plot_inputs(slope_data, figsize=(12, 6), mode='seep', mat_table=False, tab_loc='top', save_png=True)

# Build mesh
element_type = 'tri3'
polygons = build_polygons(slope_data, debug=True)

# Auto-size mesh based on domain width
x_range = [min(x for x, _ in slope_data['ground_surface'].coords),
           max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 120

mesh = build_mesh_from_polygons(polygons, target_size, element_type)
mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"
export_mesh_to_json(mesh, mesh_file)

# Build seepage data and solve
seep_data = build_seep_data(mesh, slope_data)
plot_seep_data(seep_data, figsize=(12, 6), show_nodes=True, show_bc=True,
               label_elements=False, label_nodes=False)

solution = run_seepage_analysis(seep_data, tol=1e-4)

# Plot solution
plot_seep_solution(seep_data, solution, figsize=(12, 6),
                   variable="head", vectors=False, flowlines=True,
                   mesh=False, levels=20, fill_contours=False,
                   phreatic=True, save_png=True)

# Export solution for use in LEM (u="seep")
seep_file = input_path.parent / f"{input_path.stem}_seep.csv"
export_seep_solution(seep_data, solution, seep_file)
print(f"Seepage solution exported to: {seep_file}")

# Check for second set of BCs (rapid drawdown)
if slope_data.get("has_seepage_bc2"):
    print("\nSecond set of seepage boundary conditions found. Running second analysis...")
    seep_data2 = build_seep_data(mesh, slope_data, seep_bc=2)
    plot_seep_data(seep_data2, figsize=(12, 6), show_nodes=True, show_bc=True,
                   label_elements=False, label_nodes=False)
    solution2 = run_seepage_analysis(seep_data2, tol=1e-4)
    plot_seep_solution(seep_data2, solution2, figsize=(12, 6),
                       variable="head", vectors=False, flowlines=True,
                       mesh=False, levels=20, fill_contours=False,
                       phreatic=True, save_png=True)
    seep_file2 = input_path.parent / f"{input_path.stem}_seep2.csv"
    export_seep_solution(seep_data2, solution2, seep_file2)
```

---

## FEM Analysis Code

Use this pattern to run finite element SSRM (Shear Strength Reduction Method) analysis.
Based on `main_fem.py`.

```python
from pathlib import Path
from xslope.fem import build_fem_data, solve_fem, solve_ssrm, print_reinforcement_summary, print_pile_summary
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json, extract_constraint_line_geometry
from xslope.plot import plot_inputs
from xslope.plot_fem import plot_fem_results, plot_fem_data

input_file = "inputs/my_problem.xlsx"
input_path = Path(input_file)
slope_data = load_slope_data(input_file)
plot_inputs(slope_data, mode='fem', tab_loc='top', save_png=True)

# Build mesh (tri6 or quad8 recommended for FEM)
element_type = 'tri6'  # 'quad4', 'quad8', 'tri3', 'tri6'

# extract_constraint_line_geometry handles both reinforcement AND pile lines
constraint_lines, n_reinf, n_pile = extract_constraint_line_geometry(slope_data)
polygons = build_polygons(slope_data, reinf_lines=constraint_lines)

# Auto-size mesh based on domain width
x_range = [min(x for x, _ in slope_data['ground_surface'].coords),
           max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 80

mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                element_type=element_type, lines=constraint_lines)
mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"
export_mesh_to_json(mesh, mesh_file)

# Build FEM data and plot mesh
fem_data = build_fem_data(slope_data, mesh)
plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, save_png=True)

# Run SSRM - returns a dict with 'FS', 'converged', 'last_solution', etc.
F_min = 1.0   # Lower FS bound (must converge)
F_max = 2.0   # Upper FS bound (should not converge)
result = solve_ssrm(fem_data, F_min=F_min, F_max=F_max, tolerance=0.05, debug_level=1)

if result.get("converged", False):
    print(f"\nFactor of Safety: {result['FS']:.2f}")
    print_reinforcement_summary(fem_data, result['last_solution'])
    print_pile_summary(fem_data, result['last_solution'])
    plot_fem_results(fem_data, result['last_solution'],
                     plot_type=['deformation', 'shear_strain', 'displace_vector'], save_png=True)
else:
    print(f"SSRM failed: {result.get('error', 'Unknown error')}")
```

---

## Important Guidelines

1. **Units must be consistent.** English: ft, pcf, psf. Metric: m, kN/m3, kPa. Do not mix.

2. **Profile lines go top-to-bottom.** The first profile line is the ground surface or the shallowest layer. Each subsequent line defines a deeper layer boundary. Points within each line go left-to-right.

3. **Material numbering is 1-based** in the Excel file. Mat ID 1 in the profile sheet references row 9 (first data row) of the mat sheet.

4. **Max Depth** in the profile sheet defines a horizontal bedrock surface. Set to 0 if the lowest profile line IS the base, or set to the actual bedrock elevation if deeper.

5. **For seepage-only problems**, you do NOT need circles, piezo, or non-circ sheets. Only fill main, mat (with k1, k2, kr0, h0), profile, and seep bc.

6. **For LEM-only problems**, you do NOT need seep bc or seepage material properties. Only fill main, mat (with strength properties), profile, circles (or non-circ), and optionally piezo, dloads, reinforce, piles.

7. **When interpreting diagrams**, pay attention to:
   - Scale bars and dimension labels
   - Slope ratios (e.g., 2H:1V means for every 2 horizontal, 1 vertical)
   - **Water table identification**: A water table is indicated by an **inverted triangle symbol** (▽) on the diagram. Do NOT assume a dashed line is a water table unless it is accompanied by this symbol or is explicitly labeled. Dashed lines may represent other features (e.g., material boundaries, construction lines).
   - **Ponded/standing water**: If the water table (▽) is shown ABOVE the ground surface, there is ponded water. This MUST be modeled as a distributed load (dloads) on the ground surface. The normal stress at each surface point = γ_w × (water_elevation - ground_elevation). This applies even for phi=0 total stress analyses — the water load is part of the problem definition, not an optional refinement. Never skip it.
   - Piezometric surfaces: typically shown as dashed/blue lines with explicit labels
   - Material boundaries shown as solid lines between differently hatched/colored zones
   - Property tables typically shown in the diagram legend

8. **Never overwrite formula cells.** The template contains XLOOKUP formulas (e.g., row 6 in the profile sheet auto-populates material names from the mat sheet). Overwriting a formula cell with a plain value causes the `calcChain.xml` to become inconsistent, and Excel will show a recovery error. Only write to data-entry cells.

9. **Always validate** by plotting inputs before running analysis. If geometry looks wrong, fix the template first.

9. **Seepage material properties**: For fully saturated problems, kr0 and h0 are ignored but must still have placeholder values. For partially saturated (unconfined) problems, typical values are kr0=0.001 to 0.01 and h0=-1.

10. **When the user says "find the factor of safety"**, default to auto_search with Spencer's method unless they specify otherwise. Spencer's method satisfies both force and moment equilibrium and works for both circular and non-circular surfaces.
