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
import shutil
import zipfile
import re
import os
import tempfile
import subprocess
from lxml import etree

from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons

# --- Surgical xlsx writer (from xslope skill, with recalc-on-open fix) ---
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
    # Column map (1-based): name=2, gamma=3, option=4, c=5, phi=6, cp=7,
    # r_elev=8, d=9, psi=10, u=11, sigmas s(g..psi)=12..17,
    # k1=18, k2=19, alpha=20, kr0=21, h0=22, E=23, nu=24
    num_cols = {
        'gamma': 3, 'c': 5, 'phi': 6, 'cp': 7, 'r_elev': 8, 'd': 9, 'psi': 10,
        'sigma_gamma': 12, 'sigma_c': 13, 'sigma_phi': 14, 'sigma_cp': 15,
        'sigma_d': 16, 'sigma_psi': 17, 'k1': 18, 'k2': 19, 'alpha': 20,
        'kr0': 21, 'h0': 22, 'E': 23, 'nu': 24,
    }
    updates['mat'] = {}
    for i, mat in enumerate(data['materials']):
        row = 9 + i
        for key, col in (('name', 2), ('option', 4), ('u', 11)):
            s = _clean_str(mat.get(key))
            if s is not None:
                updates['mat'][cell_ref(row, col)] = s
        for key, col in num_cols.items():
            v = mat.get(key)
            if _is_num(v):
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
