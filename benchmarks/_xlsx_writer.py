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


# --- Convenience writers for the data-entry sheets ---

def new_file(dst, src="docs/inputs/input_template.xlsx"):
    shutil.copy(src, dst)
    return dst


def main_cells(gamma_w=9.81, tcrack_depth=0, tcrack_water=0, seismic=0):
    return {'D8': gamma_w, 'D9': tcrack_depth, 'D10': tcrack_water, 'D11': seismic}


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
                   template="docs/inputs/input_template.xlsx"):
    """Cells for one 'mat' sheet material row, located BY HEADER NAME.

    The mat-sheet header row AND its column order have both shifted across template
    versions (v11 inserted unsat/vg_a/vg_n; v16 t_cut/E/nu; v17 relocated phi_b/s_cap
    to Q:R, right of ru), so every field is placed through fileio.mat_header_cols() —
    the same name->column map load_slope_data reads back — never by hardcoded index.
    The data row is the located header row + ``mat_num`` (mat_num is 1-based). Header
    lookup is underscore-insensitive, matching the loader ('phi_b' -> 'phib' etc.)."""
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
    return cells


def profile_line_cells(line_num, mat_id, points):
    col_offset = (line_num - 1) * 3
    x_col = 1 + col_offset
    y_col = 2 + col_offset
    cells = {cell_ref(5, y_col): mat_id}
    for i, (x, y) in enumerate(points):
        cells[cell_ref(8 + i, x_col)] = x
        cells[cell_ref(8 + i, y_col)] = y
    return cells


def polygon_cells(poly_num, mat_id, points):
    """Cells for the 'polygon' sheet: explicit closed boundary ring.

    Mat ID value sits at row 5 (col = 2 + 3*(poly_num-1)); points start at row 8
    with x/y in cols (1,2) for polygon #1, (4,5) for #2, etc. Points may be CW or
    CCW; repeat the first point at the end to close explicitly.
    """
    col_offset = (poly_num - 1) * 3
    x_col = 1 + col_offset
    y_col = 2 + col_offset
    cells = {cell_ref(5, y_col): mat_id}
    for i, (x, y) in enumerate(points):
        cells[cell_ref(8 + i, x_col)] = x
        cells[cell_ref(8 + i, y_col)] = y
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


def dload_cells(load_num, points):
    """Cells for the 'dloads' sheet. points: list of (x, y, n) — normal stress n
    applied to the ground surface, linearly interpolated between points.
    Layout: 4-column blocks starting at col 2 (X, Y, N), points from row 4."""
    col0 = 2 + (load_num - 1) * 4
    cells = {}
    for i, (x, y, n) in enumerate(points):
        cells[cell_ref(4 + i, col0)] = x
        cells[cell_ref(4 + i, col0 + 1)] = y
        cells[cell_ref(4 + i, col0 + 2)] = n
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
