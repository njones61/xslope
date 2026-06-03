"""Build poly_test/xslope_lost_lake_poly.xlsx by copying the v10 template and
populating it with the Lost Lake inputs, using the polygon sheet (derived from
the profile-line geometry via build_polygons) instead of the profile sheet."""
import shutil
import zipfile
import re
import os
import tempfile
import subprocess
from lxml import etree

from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons

# --- Surgical xlsx writer (from xslope skill) ---
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


def _force_full_recalc(xml_text):
    """Set fullCalcOnLoad="1" on <calcPr> so Excel recomputes formula caches on open.
    Writing cell values at the XML level can't fire a recalc event, so dependent
    formulas (e.g. the XLOOKUP material-name lookups) would otherwise stay blank."""
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


def _reset_view(xml_bytes):
    """Reset a sheet's saved scroll/selection to A1 so populated cells are visible on
    open (template ships polygon scrolled to topLeftCell="AA1" -> looks empty)."""
    t = xml_bytes.decode('utf-8')
    t = re.sub(r'\s+topLeftCell="[^"]*"', '', t)
    t = re.sub(r'<selection\b[^>]*/>', '<selection activeCell="A1" sqref="A1"/>', t)
    return t.encode('utf-8')


def _drop_calcchain(tmpdir, paths_to_zip, zf_read):
    """Remove xl/calcChain.xml + its references so Excel rebuilds it on open.
    Without this, changing a formula's precedent cell leaves a stale calcChain and
    Excel "recovers" the sheet, discarding our edits (sheet opens blank)."""
    ct = zf_read('[Content_Types].xml').decode('utf-8')
    if 'calcChain' not in ct:
        return False
    ct = re.sub(r'<Override\s+PartName="/xl/calcChain\.xml"[^>]*/>', '', ct)
    with open(os.path.join(tmpdir, '[Content_Types].xml'), 'wb') as f:
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
updates['mat'] = {}
for i, (name, k1, k2, alpha, kr0, h0) in enumerate(seep_props):
    row = 9 + i
    updates['mat'][cell_ref(row, 2)] = name
    updates['mat'][cell_ref(row, 18)] = k1
    updates['mat'][cell_ref(row, 19)] = k2
    updates['mat'][cell_ref(row, 20)] = alpha
    updates['mat'][cell_ref(row, 21)] = kr0
    updates['mat'][cell_ref(row, 22)] = h0

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
