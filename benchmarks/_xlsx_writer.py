"""Surgical xlsx cell writer (preserves all template formatting/charts/drawings).

Lifted verbatim from the xslope Claude skill (docs/usage/claude/xslope.md). Modifies
cell values via regex on the raw sheet XML, forces a full recalc on open, and drops the
stale calcChain so dependent XLOOKUP formulas (profile/polygon material names) refresh
instead of opening blank. Used by the benchmark build_*.py scripts.
"""

import shutil
import zipfile
import re
import os
import tempfile
import subprocess
from lxml import etree

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


def _is_blank(value):
    """True for a value that must be written as an EMPTY cell rather than text.

    ``None``/NaN/inf mean "unset". They must clear the cell, never land in it as the
    literal string 'None' or 'nan' — a numeric cell holding "nan" makes openpyxl's
    reader raise on the next load. Mirrors fileio._modify_existing_cell so both
    writers blank identically; a builder needs this to CLEAR a template default
    (e.g. the stock 'day' time unit on a model that declares no time unit).
    """
    if value is None:
        return True
    return isinstance(value, float) and (value != value or value in (float('inf'), float('-inf')))


def _modify_existing_cell(cell_xml, value):
    open_match = re.match(r'(<c\s[^>]*?)(/?>)', cell_xml)
    if not open_match:
        return cell_xml
    open_tag_attrs = open_match.group(1)
    open_tag_attrs = re.sub(r'\s+t="[^"]*"', '', open_tag_attrs)
    if _is_blank(value):
        return f'{open_tag_attrs}/>'
    if isinstance(value, float):
        value = round(value, 10)
    if isinstance(value, str):
        return f'{open_tag_attrs} t="inlineStr"><is><t>{value}</t></is></c>'
    else:
        return f'{open_tag_attrs}><v>{value}</v></c>'


def _build_new_cell(ref, value):
    if _is_blank(value):
        return f'<c r="{ref}"/>'
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
    """Surgically update cell values in xlsx without losing formatting.

    updates: {sheet_name: {cell_ref: value, ...}, ...}
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
                _t = rid_map[rid]
                sheet_paths[s.get('name')] = _t.lstrip('/') if _t.startswith('/') else f'xl/{_t}'
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


# --- Convenience writers for the data-entry sheets ---

def new_file(dst, src="docs/inputs/input_template.xlsx"):
    shutil.copy(src, dst)
    return dst


_TEMPLATE_VERSION_CACHE = {}


def _template_version(path):
    """Template version from the destination workbook's main!D5, or 0 if unreadable.

    Read from the FILE, not assumed, for the same reason material_cells locates its
    columns by header name: the main-sheet global block has moved across template
    versions (v18 inserted the Units and Time selectors at D8:D9 and pushed
    gamma_w/tcrack/crack-water/seismic down two rows). A builder that hardcodes the
    <=v17 positions writes gamma_w into the Units selector and the tension-crack
    depth into gamma_w — which is exactly how the LEM sample files came to read back
    with gamma_water = 0 and no unit system.
    """
    key = os.path.abspath(path)
    if key not in _TEMPLATE_VERSION_CACHE:
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=False)
        try:
            raw = wb['main'].cell(row=5, column=4).value       # main!D5
        finally:
            wb.close()
        try:
            _TEMPLATE_VERSION_CACHE[key] = int(float(raw))
        except (TypeError, ValueError):
            _TEMPLATE_VERSION_CACHE[key] = 0
    return _TEMPLATE_VERSION_CACHE[key]


def _unit_label(unit_system, gamma_w):
    """The template's Units-selector text ('SI' / 'Imperial'), or None if unlabeled.

    When the caller does not declare a system, it is inferred from the unit weight of
    water — the one quantity physics pins — through xslope.units, the same inference
    the loader applies to a pre-selector file. Off-band (e.g. a seawater override)
    stays None: unlabeled, never mislabeled.
    """
    from xslope.units import normalize_unit_system, infer_system_from_gamma_water
    us = normalize_unit_system(unit_system) or infer_system_from_gamma_water(gamma_w)
    return 'SI' if us == 'si' else 'Imperial' if us == 'imperial' else None


def main_cells(gamma_w=9.81, tcrack_depth=0, tcrack_water=0, seismic=0,
               unit_system=None, time_unit=None,
               lem_method=None, num_slices=None, k0=None, tension_srf=None,
               element_type=None, target_size=None, ssrm_f_min=None, ssrm_f_max=None,
               side_bc=None, water_loads='auto',
               template="docs/inputs/input_template.xlsx"):
    """Cells for the 'main' sheet's global block, placed BY TEMPLATE VERSION.

    v18 declares the unit system (D8) and time unit (D9) above the numeric globals,
    so gamma_w/tcrack/crack-water/seismic sit at D10:D13; <=v17 has no selectors and
    they sit at D8:D11. Both selectors are written unconditionally on v18 — blank
    when undeclared — so the template's own pre-filled defaults ('Imperial', 'day')
    can never leak into a built file. This mirrors fileio.save_slope_data_to_xlsx.

    v19 adds eight optional run options at D14:D21 (LEM method, number of slices, K0,
    Tension SRF, mesh element type, mesh target size, SSRM F min/F max). They are
    written unconditionally too, for the leak reason above — the blank template ships
    D17 pre-filled 'YES', so a builder that declares nothing must still clear it.
    v21 adds the Side BC selector at D22, which ships pre-filled 'rollers' and is
    cleared the same way.

    v22 adds the water-load mode at D23, and it is the one cell here whose default
    is a value rather than a blank: 'auto' means the engine derives the standing
    water from the model's own water definition, which is what a built file should
    say. It is written EXPLICITLY for the same reason the others are cleared --
    inheriting the template's pre-fill would make the mode an accident of which
    template the build happened to run against, and this cell decides whether a
    reservoir is counted once or twice. A builder whose model must keep its water
    typed in passes water_loads='manual' with its reason.
    """
    ver = _template_version(template)
    if ver >= 18:
        cells = {'D8': _unit_label(unit_system, gamma_w),
                 'D9': str(time_unit) if time_unit else None,
                 'D10': gamma_w, 'D11': tcrack_depth,
                 'D12': tcrack_water, 'D13': seismic}
        if ver >= 19:
            cells.update({
                'D14': str(lem_method) if lem_method else None,
                'D15': int(num_slices) if num_slices is not None else None,
                'D16': k0,
                'D17': None if tension_srf is None else ('YES' if tension_srf else 'NO'),
                'D18': str(element_type) if element_type else None,
                'D19': target_size,
                'D20': ssrm_f_min,
                'D21': ssrm_f_max,
            })
        if ver >= 21:
            # v21 side BC (D22). Written unconditionally for the leak reason above:
            # the blank template ships D22 pre-filled 'rollers', so a builder that
            # declares nothing must still clear it.
            cells['D22'] = str(side_bc).lower() if side_bc else None
        if ver >= 22:
            wl = str(water_loads or 'auto').strip().lower()
            if wl not in ('auto', 'manual'):
                raise ValueError(f"water_loads must be 'auto' or 'manual', not "
                                 f"{water_loads!r}")
            cells['D23'] = wl
        return cells
    return {'D8': gamma_w, 'D9': tcrack_depth, 'D10': tcrack_water, 'D11': seismic}


def search_window_cells(entry_x_min=None, entry_x_max=None,
                        exit_x_min=None, exit_x_max=None,
                        center_box_x_min=None, center_box_x_max=None,
                        center_box_y_min=None, center_box_y_max=None,
                        max_tangent_depth=None, min_slip_depth=None,
                        template="docs/inputs/input_template.xlsx"):
    """Cells for the v19 'circles' sheet search window (K8:K17).

    The ten optional circular-search limits, in the sheet's own order (see
    fileio.SEARCH_WINDOW_KEYS). Returns {} on a pre-v19 template, which has no such
    block. Every value is written, so an unset limit leaves a blank cell.
    """
    if _template_version(template) < 19:
        return {}
    from xslope.fileio import SEARCH_WINDOW_KEYS
    vals = dict(entry_x_min=entry_x_min, entry_x_max=entry_x_max,
                exit_x_min=exit_x_min, exit_x_max=exit_x_max,
                center_box_x_min=center_box_x_min, center_box_x_max=center_box_x_max,
                center_box_y_min=center_box_y_min, center_box_y_max=center_box_y_max,
                max_tangent_depth=max_tangent_depth, min_slip_depth=min_slip_depth)
    return {cell_ref(8 + i, 11): vals[k] for i, k in enumerate(SEARCH_WINDOW_KEYS)}


_MAT_HEADER_CACHE = {}


def _mat_header(template):
    """(header_row, {stripped_header: col}) for ``template``'s 'mat' sheet, cached per
    path. Delegates to fileio.mat_header_cols so the benchmark writer resolves columns
    exactly as load_slope_data reads them back — see material_cells."""
    from xslope.fileio import mat_header_cols
    key = os.path.abspath(template)
    if key not in _MAT_HEADER_CACHE:
        _MAT_HEADER_CACHE[key] = mat_header_cols(template)
    return _MAT_HEADER_CACHE[key]


def material_cells(mat_num, name, gamma, option, c, phi, u,
                   k1=None, k2=None, alpha=None, kr0=None, h0=None, E=None, nu=None,
                   phi_b=None, s_cap=None, t_cut=None,
                   sigma_gamma=None, sigma_c=None, sigma_phi=None,
                   template="docs/inputs/input_template.xlsx"):
    """Cells for one 'mat' sheet material row, located BY HEADER NAME.

    The mat-sheet header row AND its column order have both shifted across template
    versions (v11 inserted unsat/vg_a/vg_n; v16 t_cut/E/nu; v17 relocated phi_b/s_cap
    to Q:R, right of ru), so every field is placed through fileio.mat_header_cols() —
    the same name->column map load_slope_data reads back — never by hardcoded index.
    The data row is the located header row + ``mat_num`` (mat_num is 1-based). Header
    lookup is underscore-insensitive, matching the loader ('phi_b' -> 'phib' etc.).

    ``template`` MUST name the same workbook ``new_file`` copied to create the
    destination — the refs returned here are resolved against ITS column layout, so
    a destination built from a different template takes every field one column off
    per inserted header, silently, and only visible on the next load."""
    header_row, cols = _mat_header(template)
    row = header_row + mat_num
    cells = {}

    def put(header, val):
        if val is None:
            return
        col = cols.get(header.replace('_', ''))
        if col is not None:
            cells[cell_ref(row, col)] = val

    put('mat', mat_num)
    put('name', name)
    put('g', gamma)
    put('option', option)
    put('c', c)
    put('f', phi)          # phi's mat-sheet header is 'f'
    put('u', u)
    put('k1', k1)
    put('k2', k2)
    put('alpha', alpha)
    put('kr0', kr0)
    put('h0', h0)
    put('E', E)
    put('nu', nu)
    put('phi_b', phi_b)    # v17 matric-suction pair (LEM & FEM), now at cols Q:R
    put('s_cap', s_cap)
    put('t_cut', t_cut)
    # Reliability/probabilistic standard deviations. Their headers are parenthesized
    # ('s(g)', 's(c)', 's(f)') and carry no underscore, so they resolve unchanged
    # through the same by-name lookup.
    put('s(g)', sigma_gamma)
    put('s(c)', sigma_c)
    put('s(f)', sigma_phi)
    return cells


def profile_line_cells(line_num, mat_id, points, size=None,
                       template="docs/inputs/input_template.xlsx"):
    """Cells for one 'profile' sheet block, placed BY TEMPLATE VERSION.

    v21 inserts an optional 'Size:' row (Excel row 7) between the material-name echo
    and the coordinates, so the vertices start at row 9 instead of row 8. Reading the
    destination's own version, rather than hardcoding, is the same discipline
    material_cells and main_cells already follow — a hardcoded v20 map writes the first
    vertex into the Size cell and silently loses a point.
    """
    ver = _template_version(template)
    coord_row = 9 if ver >= 21 else 8
    col_offset = (line_num - 1) * 3
    x_col = 1 + col_offset
    y_col = 2 + col_offset
    cells = {cell_ref(5, y_col): mat_id}
    if ver >= 21:
        cells[cell_ref(7, y_col)] = size          # blank when unset
    for i, (x, y) in enumerate(points):
        cells[cell_ref(coord_row + i, x_col)] = x
        cells[cell_ref(coord_row + i, y_col)] = y
    return cells


def polygon_cells(poly_num, mat_id, points, poly_type=None, size=None,
                  template="docs/inputs/input_template.xlsx"):
    """Cells for one 'polygon' sheet block, placed BY TEMPLATE VERSION.

    Points may be CW or CCW; repeat the first point at the end to close explicitly.
    x/y go in cols (1,2) for polygon #1, (4,5) for #2, etc.

      <= v20   Mat ID row 5 (negative = an SSR sentinel),        vertices from row 8
      v21      Type row 5, Mat ID row 6, Size row 8,             vertices from row 10

    ``poly_type`` is the v21 Type word ('material' — the default — 'ssr reduce',
    'ssr hold', 'ssr elastic', 'refine'); on a pre-v21 template an SSR type is written
    as its sentinel Mat ID instead, so a builder can name the kind either way.
    """
    from xslope.fileio import (POLYGON_TYPE_WORDS, SSR_ZONE_SENTINELS,
                               SSR_ZONE_LABELS)
    ver = _template_version(template)
    word = (poly_type or 'material').strip().lower()
    if word not in POLYGON_TYPE_WORDS:
        raise ValueError(f"unknown polygon type {poly_type!r}; expected one of "
                         f"{', '.join(sorted(POLYGON_TYPE_WORDS))}")
    kind = POLYGON_TYPE_WORDS[word]
    col_offset = (poly_num - 1) * 3
    x_col = 1 + col_offset
    y_col = 2 + col_offset
    if ver >= 21:
        coord_row = 10
        cells = {cell_ref(5, y_col): word,
                 cell_ref(6, y_col): mat_id if kind == 'material' else None,
                 cell_ref(8, y_col): size}
    else:
        coord_row = 8
        if kind == 'refine':
            raise ValueError("a 'refine' polygon has no pre-v21 representation; "
                             f"template {template} is version {ver}")
        _sentinel = {v: k for k, v in SSR_ZONE_SENTINELS.items()}
        cells = {cell_ref(5, y_col):
                 mat_id if kind == 'material' else _sentinel[kind]}
    for i, (x, y) in enumerate(points):
        cells[cell_ref(coord_row + i, x_col)] = x
        cells[cell_ref(coord_row + i, y_col)] = y
    return cells


def circle_cells(num, xo, yo, option="Depth", depth=None, xi=None, yi=None, radius=None):
    row = 2 + num
    cells = {cell_ref(row, 1): num, cell_ref(row, 2): xo,
             cell_ref(row, 3): yo, cell_ref(row, 4): option}
    if option == "Depth" and depth is not None:
        cells[cell_ref(row, 5)] = depth
    elif option == "Intercept":
        cells[cell_ref(row, 6)] = xi
        cells[cell_ref(row, 7)] = yi
    elif option == "Radius" and radius is not None:
        cells[cell_ref(row, 8)] = radius
    return cells


def noncirc_cells(points):
    """points: list of (x, y, movement)."""
    cells = {}
    for i, (x, y, movement) in enumerate(points):
        cells[cell_ref(3 + i, 1)] = x
        cells[cell_ref(3 + i, 2)] = y
        cells[cell_ref(3 + i, 3)] = movement
    return cells


def piezo_cells(points):
    cells = {}
    for i, (x, y) in enumerate(points):
        cells[cell_ref(4 + i, 1)] = x
        cells[cell_ref(4 + i, 2)] = y
    return cells


def dload_cells(load_num, points, direction=None,
                template="docs/inputs/input_template.xlsx"):
    """Cells for one 'dloads' (or 'dloads (2)') block, placed BY TEMPLATE VERSION.

    points: list of (x, y, n) — load intensity n along the line, linearly interpolated
    between points. Blocks are 4 columns apart starting at col 2 (X, Y, N).

    v21 adds a Direction cell in the block's N column and moves the data down one row
    (4 -> 5); v22 puts a conditional note above the blocks and pushes everything down
    two more rows (Direction at Excel row 5, data from row 7). ``direction`` is 'normal'
    (the default, written blank so the sheet reads the way it always has) or 'vertical'.
    """
    ver = _template_version(template)
    if ver >= 22:
        data_row, dir_row = 7, 5
    elif ver >= 21:
        data_row, dir_row = 5, 3
    else:
        data_row, dir_row = 4, None
    col0 = 2 + (load_num - 1) * 4
    cells = {}
    if dir_row is not None:
        word = (direction or 'normal').strip().lower()
        if word not in ('normal', 'vertical'):
            raise ValueError(f"unknown dload direction {direction!r}; expected "
                             "'normal' or 'vertical'")
        cells[cell_ref(dir_row, col0 + 2)] = None if word == 'normal' else word
    elif direction and str(direction).strip().lower() != 'normal':
        raise ValueError(f"a '{direction}' distributed load has no pre-v21 "
                         f"representation; template {template} is version {ver}")
    for i, (x, y, n) in enumerate(points):
        cells[cell_ref(data_row + i, col0)] = x
        cells[cell_ref(data_row + i, col0 + 1)] = y
        cells[cell_ref(data_row + i, col0 + 2)] = n
    return cells


def seep_bc_cells(exit_face=None, head1=None, head1_pts=None,
                  head2=None, head2_pts=None):
    """Cells for the 'seep bc' sheet.

    Layout (1-based cols): Exit Face x/y = cols 2/3; Specified Head #1 value at
    (row 3, col 6), points x/y = cols 5/6; Specified Head #2 value at (row 3,
    col 9), points x/y = cols 8/9. All polyline points start at row 5.
    """
    cells = {}
    if exit_face:
        for i, (x, y) in enumerate(exit_face):
            cells[cell_ref(5 + i, 2)] = x
            cells[cell_ref(5 + i, 3)] = y
    if head1 is not None:
        cells[cell_ref(3, 6)] = head1
    if head1_pts:
        for i, (x, y) in enumerate(head1_pts):
            cells[cell_ref(5 + i, 5)] = x
            cells[cell_ref(5 + i, 6)] = y
    if head2 is not None:
        cells[cell_ref(3, 9)] = head2
    if head2_pts:
        for i, (x, y) in enumerate(head2_pts):
            cells[cell_ref(5 + i, 8)] = x
            cells[cell_ref(5 + i, 9)] = y
    return cells


# --- v22 water loads: the model states the reservoir, the file expresses it ---

def _ground_surface(slope_data):
    """The ground surface of a model a BUILDER assembled, derived from its geometry.

    A builder starts from a donor file and replaces its geometry, so the
    ``ground_surface`` still sitting on the dict belongs to the donor -- the loader
    derives that key, and nothing has re-derived it since. Measuring a water line
    against it would compare this model's pool to somebody else's hillside, which
    is a silent wrong answer rather than an error. So it is rebuilt here from the
    model's OWN profile lines or polygons, exactly as load_slope_data does.
    """
    from shapely.geometry import Polygon, LineString
    from xslope.fileio import build_polygons, build_ground_surface_from_polygons
    polys = []
    if slope_data.get('profile_lines'):
        polys = [{'polygon': Polygon(p['coords']), 'mat_id': p['mat_id']}
                 for p in build_polygons(slope_data={
                     'profile_lines': slope_data['profile_lines'],
                     'max_depth': slope_data.get('max_depth')})]
    elif slope_data.get('polygons'):
        polys = [{'polygon': p['polygon'], 'mat_id': p['mat_id']}
                 for p in slope_data['polygons']]
    if not polys:
        return slope_data.get('ground_surface') or LineString([])
    return build_ground_surface_from_polygons(polys)[0]


def emit_water_mode(slope_data, label=''):
    """Express a model's standing water the way a v22 file states it.

    A builder transcribes the source problem, and published problems draw ponded
    water as a load: a hydrostatic block on the submerged face, read off the
    figure. That transcription is evidence — it is how we show the figure was read
    correctly — so it stays in the builder. What changes is the FILE: in automatic
    mode the reservoir is stated once, as the water definition, and the dloads
    sheet carries the loads a user typed for their own reasons.

    So this removes exactly the blocks the derivation reproduces, per stage, and
    keeps everything else verbatim. Nothing is removed on the strength of looking
    like water: a block is removed only when the derived load matches it on
    pressure, resultant AND reach, which makes every rebuild a re-proof that the
    derivation still reproduces what the builder transcribed.

    A transcribed block over the pool's reach that the derivation does NOT
    reproduce raises, rather than being dropped or silently kept. That case is a
    disagreement between the figure we read and the water definition we wrote, and
    it is the finding the corpus migration exists to surface.

    ``manual`` models are returned unchanged: their water is typed in on purpose,
    each with a recorded reason (see the migration manifest).

    Returns a new dict; the caller's model is not mutated.
    """
    from xslope.water import derive_water_loads, match_water_blocks
    mode = str(slope_data.get('water_loads') or 'auto').strip().lower()
    if mode == 'manual':
        return dict(slope_data)
    out = dict(slope_data)
    out['water_loads'] = 'auto'
    out['ground_surface'] = _ground_surface(out)
    for stage, key, dirs_key, sheet in ((1, 'dloads', 'dload_dirs', 'dloads'),
                                        (2, 'dloads2', 'dload2_dirs', 'dloads (2)')):
        blocks = list(out.get(key) or [])
        if not blocks:
            continue
        derived = derive_water_loads(out, stage=stage)
        if not derived['blocks']:
            continue
        found = match_water_blocks(blocks, derived['blocks'])
        if found['missing']:
            a, b = derived['reach']
            over = []
            for i in found['unmatched']:
                xs = [float(p['X']) for p in blocks[i]]
                if max(xs) >= a - 1e-9 and min(xs) <= b + 1e-9:
                    over.append(i + 1)
            if over:
                raise ValueError(
                    f"{label or 'this model'}: the {sheet} sheet carries block(s) "
                    f"{over} over the pool's reach (x = {a:.6g} to {b:.6g}) that the "
                    f"load derived from {derived['source']} does not reproduce. "
                    f"Either the transcription or the water definition is wrong; "
                    f"set water_loads='manual' with a recorded reason if the model "
                    f"is meant to keep its water typed in.")
        drop = {i for i, _j, _m in found['pairs']}
        if not drop:
            continue
        dirs = list(out.get(dirs_key) or [])
        while len(dirs) < len(blocks):
            dirs.append('normal')
        out[key] = [b for i, b in enumerate(blocks) if i not in drop]
        out[dirs_key] = [d for i, d in enumerate(dirs) if i not in drop]
    return out
