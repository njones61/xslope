"""Build the v12 input template from the v11 master by direct XML surgery.

Edits the xlsx XML inside the zip (same approach as fix_cf_rules.py) because an
openpyxl load/save cycle is destructive for this workbook: it guts the plot
sheet's 87-series scatter chart and rewrites the rich-text headers. Untouched
zip members are copied byte-for-byte.

v12 changes (the full manifest, reproduced here):
  main       - template_version 11 -> 12 (cell D5)
  mat        - insert g-sat after g; pow_a..pow_d after psi; ru after u
               (adjacent placement per plan; legend, CF, DV, merges, widths,
               chartless) + pow/ru legend entries + fix J2 "presssue" typo
  reinforce  - insert Label at col B; append type/dir/appl/Tend1/Tend2/spacing
               at M-R with dropdown validation; chart series refs shifted
  piles      - append appl column with dropdown validation
  lloads     - new sheet: # | Label | x | y | P | angle

Usage:  python make_template_v12.py [--out PATH]
Reads docs/inputs/input_template_v11.xlsx (archived master), writes the v12 master.
⚠ The shipped master contains hand-edited legend text boxes this script cannot
reproduce — rebuild only to a scratch --out path unless you mean to discard them.
"""

import re
import sys
import zipfile

SRC = 'docs/inputs/input_template_v11.xlsx'
OUT = 'docs/inputs/input_template.xlsx'

MAT_SHEET = 'xl/worksheets/sheet3.xml'
REINF_SHEET = 'xl/worksheets/sheet11.xml'
PILES_SHEET = 'xl/worksheets/sheet12.xml'
MAIN_SHEET = 'xl/worksheets/sheet1.xml'
CHART = 'xl/charts/chart1.xml'
WORKBOOK = 'xl/workbook.xml'
WB_RELS = 'xl/_rels/workbook.xml.rels'
CONTENT_TYPES = '[Content_Types].xml'
LLOADS_SHEET = 'xl/worksheets/sheet15.xml'


def col_to_idx(col):
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - 64)
    return n


def idx_to_col(n):
    s = ''
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


# --- mat column mapping: g-sat after C(g), pow_a..d after J(psi), ru after K(u)
def mat_map(idx):
    if idx <= 3:
        return idx
    if idx <= 10:          # D(option)..J(psi)
        return idx + 1
    if idx == 11:          # K(u)
        return idx + 5
    return idx + 6         # L(s(g)) onward


# reinforce v12 (Norm, 2026-07-09): Label at B, then columns grouped by
# analysis type like piles —
#   # Label x1 y1 x2 y2 | Type Dir Appl (LEM, green) |
#   Tmax Lp1 Lp2 Tend1 Tend2 Spacing (both, red) | Tres E Area (FEM, blue)
# v11 source order: # x1 y1 x2 y2 Tmax Tres Lp1 Lp2 E Area
REINF_PERM = {2: 3, 3: 4, 4: 5, 5: 6,      # geometry right one (Label at B)
              6: 10,                        # Tmax
              7: 16,                        # Tres -> blue block
              8: 11, 9: 12,                 # Lp1, Lp2
              10: 17, 11: 18}               # E, Area -> blue block


def reinf_map(idx):
    return REINF_PERM.get(idx, idx)


def shift_cell_refs(xml, mapping):
    """Shift the column of every <c r="XX9"> cell reference."""
    def sub(m):
        return f'<c r="{idx_to_col(mapping(col_to_idx(m.group(1))))}{m.group(2)}"'
    return re.sub(r'<c r="([A-Z]{1,3})(\d+)"', sub, xml)


def shift_sqref(ref, mapping):
    """Shift an A1:B2-style range (or single ref), preserving any $."""
    def sub(m):
        return f'{m.group(1)}{idx_to_col(mapping(col_to_idx(m.group(2))))}{m.group(3)}{m.group(4)}'
    return re.sub(r'(\$?)([A-Z]{1,3})(\$?)(\d+)', sub, ref)


def shift_formula_cols(xml, mapping, tags=('formula', 'formula1')):
    """Shift $COL refs inside <formula>/<formula1> elements."""
    def fix(m):
        inner = re.sub(r'\$([A-Z]{1,3})(?=\$?\d)',
                       lambda mm: '$' + idx_to_col(mapping(col_to_idx(mm.group(1)))),
                       m.group(2))
        return m.group(1) + inner + m.group(3)
    for tag in tags:
        xml = re.sub(rf'(<{tag}>)(.*?)(</{tag}>)', fix, xml)
    return xml


def shift_block_refs(xml, mapping, tag_attr_patterns):
    """Shift ranges in attributes like sqref="..." or ref="..."."""
    for pat in tag_attr_patterns:
        xml = re.sub(pat + r'"([^"]+)"',
                     lambda m: m.group(0).replace(m.group(1), ' '.join(
                         shift_sqref(part, mapping) for part in m.group(1).split())),
                     xml)
    return xml


def shift_cols_block(xml, mapping, new_cols_from_left):
    """Rewrite <cols> per mapping; new columns clone the left neighbor's def."""
    m = re.search(r'<cols>(.*?)</cols>', xml)
    if not m:
        return xml
    per_col = {}
    for cm in re.finditer(r'<col ([^/]*?)/>', m.group(1)):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', cm.group(1)))
        for i in range(int(attrs['min']), int(attrs['max']) + 1):
            a = dict(attrs)
            per_col[mapping(i)] = a
    for new_idx, left_idx in new_cols_from_left.items():
        if mapping(left_idx) in per_col:
            per_col[new_idx] = dict(per_col[mapping(left_idx)])
    parts = []
    for i in sorted(per_col):
        a = per_col[i]
        a['min'] = a['max'] = str(i)
        parts.append('<col ' + ' '.join(f'{k}="{v}"' for k, v in a.items()) + '/>')
    return xml[:m.start()] + '<cols>' + ''.join(parts) + '</cols>' + xml[m.end():]


def insert_cell(xml, row, col, cell_xml):
    """Insert cell into <row r="row"> keeping ascending column order."""
    rm = re.search(rf'(<row r="{row}"[^>]*>)(.*?)(</row>)', xml)
    if not rm:
        raise RuntimeError(f'row {row} not found')
    body = rm.group(2)
    target = col_to_idx(col)
    out, inserted = [], False
    # lazy attr match: a self-closing cell must terminate at its own '/>' and
    # never swallow following cells to reach a later '</c>'
    for cm in re.finditer(r'<c r="([A-Z]{1,3})\d+"[^>]*?(?:/>|>.*?</c>)', body):
        existing = col_to_idx(cm.group(1))
        if not inserted and existing == target:      # replace in place
            out.append(cell_xml)
            inserted = True
            continue
        if not inserted and existing > target:
            out.append(cell_xml)
            inserted = True
        out.append(cm.group(0))
    if not inserted:
        out.append(cell_xml)
    return xml[:rm.start()] + rm.group(1) + ''.join(out) + rm.group(3) + xml[rm.end():]


def move_cell(xml, old_ref, new_ref):
    return xml.replace(f'<c r="{old_ref}"', f'<c r="{new_ref}"', 1)


def istr(ref, text, style):
    return f'<c r="{ref}" s="{style}" t="inlineStr"><is><t>{text}</t></is></c>'


def styled(ref, style):
    return f'<c r="{ref}" s="{style}"/>'


def get_style(xml, ref):
    m = re.search(rf'<c r="{ref}" s="(\d+)"', xml)
    return m.group(1) if m else '0'


def shift_rows(xml, first, delta):
    """Shift every row >= first down by delta: row elements, cell refs,
    ref/sqref ranges, and row numbers inside <formula>/<formula1>."""
    xml = re.sub(r'<row r="(\d+)"',
                 lambda m: f'<row r="{int(m.group(1)) + (delta if int(m.group(1)) >= first else 0)}"',
                 xml)
    xml = re.sub(r'<c r="([A-Z]{1,3})(\d+)"',
                 lambda m: f'<c r="{m.group(1)}{int(m.group(2)) + (delta if int(m.group(2)) >= first else 0)}"',
                 xml)

    def ref_part(p):
        return re.sub(r'(\$?[A-Z]{1,3}\$?)(\d+)',
                      lambda m: m.group(1) + str(int(m.group(2)) +
                                                 (delta if int(m.group(2)) >= first else 0)),
                      p)
    xml = re.sub(r'((?:sqref|ref)=)"([^"]+)"',
                 lambda m: f'{m.group(1)}"{ref_part(m.group(2))}"', xml)
    xml = re.sub(r'(<formula1?>)(.*?)(</formula1?>)',
                 lambda m: m.group(1) + ref_part(m.group(2)) + m.group(3), xml)
    return xml


# ---------------------------------------------------------------- mat sheet
def build_mat(xml):
    # 0. blank spacer row between the legend (rows 2-6) and the table:
    #    group titles 7->8, headers 8->9, data 9->10. NOTE: requires the
    #    _find_mat_header_row() code work — header row is no longer row 8.
    xml = shift_rows(xml, 7, 1)

    # 1. shift everything by the mapping
    xml = shift_cell_refs(xml, mat_map)
    xml = shift_block_refs(xml, mat_map, [r'<mergeCell ref=', r'<dimension ref='])
    # CF and DV ranges live in sqref="..." attributes:
    xml = re.sub(r'(sqref=)"([^"]+)"',
                 lambda m: f'{m.group(1)}"{" ".join(shift_sqref(p, mat_map) for p in m.group(2).split())}"',
                 xml)
    xml = shift_formula_cols(xml, mat_map)
    xml = shift_cols_block(xml, mat_map,
                           {4: 3, 12: 11, 13: 11, 14: 11, 15: 11, 17: 11})

    # 2. repair legend blocks that the shift split apart
    #    strength desc E3/E4 -> back to D3/D4; pore desc P3..P5 -> L3..L5
    for old, new in [('E3', 'D3'), ('E4', 'D4'), ('P3', 'L3'), ('P4', 'L4'), ('P5', 'L5')]:
        xml = move_cell(xml, old, new)

    # 3. legend additions (plain text, style cloned from existing legend cells)
    s_code = get_style(xml, 'C3')
    s_desc = get_style(xml, 'D3')
    xml = insert_cell(xml, 5, 'C', istr('C5', 'pow', s_code))
    xml = insert_cell(xml, 5, 'D',
                      istr('D5', 'Power curve envelope: t = pow_a*(sn + pow_d)^pow_b + pow_c', s_desc))
    xml = insert_cell(xml, 6, 'K', istr('K6', 'ru', s_code))
    xml = insert_cell(xml, 6, 'L', istr('L6', 'Pore pressure ratio (u = ru x sigma_v)', s_desc))

    # 3b. fix the J2 typo -> now at K2 after the shift ("Pore presssue (u) options")
    #     K2 is a shared string; point it at an inline replacement instead.
    xml = re.sub(r'<c r="K2"( s="\d+")? t="s"><v>\d+</v></c>',
                 lambda m: f'<c r="K2"{m.group(1) or ""} t="inlineStr">'
                           f'<is><t>Pore pressure (u) options</t></is></c>', xml)

    # 4. new header cells (row 9); plain style like 'option' (s=21 family).
    #    g-sat renders as gamma_sat via rich-text runs (Symbol-font "g" +
    #    subscript "sat"), matching Norm's hand styling of 2026-07-09.
    #    NOTE the machine-readable header text is therefore "gsat" — the
    #    header-name-driven reader must match "gsat", not "g-sat".
    s_hdr = get_style(xml, 'E9')     # shifted 'option' header style

    def rich_hdr(ref, base, sub, base_font='Aptos Narrow (Body)', charset=''):
        """Header as rich text: base run + subscript run. The machine-readable
        value is the concatenation (gsat, powa, ..., ru)."""
        return (f'<c r="{ref}" s="{s_hdr}" t="inlineStr"><is>'
                f'<r><rPr><sz val="12"/><color theme="0"/><rFont val="{base_font}"/>'
                f'{charset}</rPr><t>{base}</t></r>'
                f'<r><rPr><vertAlign val="subscript"/><sz val="12"/><color theme="0"/>'
                f'<rFont val="Aptos Narrow (Body)"/></rPr><t>{sub}</t></r></is></c>')

    xml = insert_cell(xml, 9, 'D',
                      rich_hdr('D9', 'g', 'sat', 'Symbol', '<charset val="2"/>'))
    for ref, base, sub in [('L9', 'pow', 'a'), ('M9', 'pow', 'b'),
                           ('N9', 'pow', 'c'), ('O9', 'pow', 'd'), ('Q9', 'r', 'u')]:
        xml = insert_cell(xml, 9, ref[:-1], rich_hdr(ref, base, sub))

    # 5. styled empty cells for new columns in every data row (>= 10)
    new_from_left = {4: 3, 12: 11, 13: 11, 14: 11, 15: 11, 17: 16}
    for rm in list(re.finditer(r'<row r="(\d+)"', xml)):
        rnum = int(rm.group(1))
        if rnum < 10:
            continue
        for new_idx, left_idx in sorted(new_from_left.items()):
            left_ref = f'{idx_to_col(left_idx)}{rnum}'
            sm = re.search(rf'<c r="{left_ref}" s="(\d+)"', xml)
            if sm:
                xml = insert_cell(xml, rnum, idx_to_col(new_idx),
                                  styled(f'{idx_to_col(new_idx)}{rnum}', sm.group(1)))

    # 6. widen the Base Values merge to cover pow_a..pow_d (C7:K7 -> C7:O7)
    xml = xml.replace('<mergeCell ref="C8:K8"/>', '<mergeCell ref="C8:O8"/>')

    # 7. new conditional-formatting rules (reuse the existing gray dxf)
    dm = re.search(r'<cfRule type="expression" dxfId="(\d+)" priority="(\d+)"', xml)
    dxf = dm.group(1)
    prios = [int(p) for p in re.findall(r'priority="(\d+)"', xml)] or [0]
    p = max(prios)

    def cf(sqref, formula):
        nonlocal p
        p += 1
        return (f'<conditionalFormatting sqref="{sqref}">'
                f'<cfRule type="expression" dxfId="{dxf}" priority="{p}">'
                f'<formula>{formula}</formula></cfRule></conditionalFormatting>')

    new_rules = (
        cf('F10:F51', '$E10="pow"') +        # c unused under pow
        cf('G10:G51', '$E10="pow"') +        # f unused under pow
        cf('H10:I51', '$E10="pow"') +        # c/p, r-elev unused under pow
        cf('L10:O51', 'AND($E10&lt;&gt;"",$E10&lt;&gt;"pow")') +  # pow_* gray once a non-pow option is chosen; blank rows stay white
        cf('Q10:Q51', 'AND($P10&lt;&gt;"",$P10&lt;&gt;"ru")') +   # ru gray once a non-ru u option is chosen; blank rows stay white
        cf('S10:S51', '$E10="pow"') +        # s(c)
        cf('T10:T51', '$E10="pow"') +        # s(f)
        cf('U10:U51', '$E10="pow"')          # s(c/p)
    )
    last_cf_end = xml.rindex('</conditionalFormatting>') + len('</conditionalFormatting>')
    xml = xml[:last_cf_end] + new_rules + xml[last_cf_end:]

    # 8. extend the dropdown lists: option gains pow, u gains ru
    xml = xml.replace('<formula1>$C$3:$C$4</formula1>', '<formula1>$C$3:$C$5</formula1>')
    xml = xml.replace('<formula1>$K$3:$K$5</formula1>', '<formula1>$K$3:$K$6</formula1>')
    return xml


# ---------------------------------------------------------- reinforce sheet
# type -> (dir, appl) preset map; also written to hidden cells Z8:AB11 that
# feed the default VLOOKUP formulas in Dir/Appl
TYPE_PRESETS = [('Geosynthetic', 'Tangent', 'Active'),
                ('Nail', 'Axial', 'Passive'),
                ('Tieback', 'Axial', 'Active'),
                ('Anchor', 'Axial', 'Active')]


def build_reinforce(xml, s_green):
    xml = shift_cell_refs(xml, reinf_map)
    xml = sort_row_cells(xml)               # REINF_PERM is non-monotonic
    xml = re.sub(r'(<dimension ref=)"([^"]+)"', r'\1"A2:AB25"', xml)
    xml = shift_cols_block(xml, reinf_map, {2: 2})  # Label clones old B width

    s_hdr_geom = get_style(xml, 'C2')   # shifted x1 header style
    s_red = get_style(xml, 'J2')        # shifted Tmax header style (LEM & FEM)
    xml = insert_cell(xml, 2, 'B', istr('B2', 'Label', s_hdr_geom))
    for ref, text, style in [('G2', 'Type', s_green), ('H2', 'Dir', s_green),
                             ('I2', 'Appl', s_green), ('M2', 'Tend1', s_red),
                             ('N2', 'Tend2', s_red), ('O2', 'Spacing', s_red)]:
        xml = insert_cell(xml, 2, ref[:-1], istr(ref, text, style))

    # dropdown source lists + the preset lookup table (hidden off to the right)
    for ref, text in [('Z2', 'Geosynthetic'), ('Z3', 'Nail'), ('Z4', 'Tieback'), ('Z5', 'Anchor'),
                      ('AA2', 'Tangent'), ('AA3', 'Axial'),
                      ('AB2', 'Active'), ('AB3', 'Passive')]:
        row = int(ref[2:]) if ref[1].isalpha() else int(ref[1:])
        col = ref[:2] if ref[1].isalpha() else ref[0]
        xml = insert_cell(xml, row, col, istr(ref, text, '0'))
    for j, (t, d, a) in enumerate(TYPE_PRESETS):
        r = 8 + j
        cells = (istr(f'Z{r}', t, '0') + istr(f'AA{r}', d, '0') + istr(f'AB{r}', a, '0'))
        prev = re.search(rf'<row r="{r - 1}"[^>]*>.*?</row>', xml)
        rm = re.search(rf'<row r="{r}"([^>]*)>(.*?)</row>', xml)
        if rm:
            xml = (xml[:rm.start()] + f'<row r="{r}"{rm.group(1)}>{rm.group(2)}{cells}</row>'
                   + xml[rm.end():])
        else:
            xml = xml[:prev.end()] + f'<row r="{r}" spans="26:28">{cells}</row>' + xml[prev.end():]

    # styled data cells for the new columns, rows 3..22
    for rnum in range(3, 23):
        c_style = re.search(rf'<c r="C{rnum}" s="(\d+)"', xml)
        if not c_style:
            continue
        s = c_style.group(1)
        for col in ['B', 'G', 'H', 'I', 'M', 'N', 'O']:
            xml = insert_cell(xml, rnum, col, styled(f'{col}{rnum}', s))
        # Dir/Appl carry default formulas driven by Type (blank when Type blank)
        for col, li in [('H', '2'), ('I', '3')]:
            f = (f'IF($G{rnum}="","",IFERROR('
                 f'VLOOKUP($G{rnum},$Z$8:$AB$11,{li},0),""))')
            xml = insert_cell(xml, rnum, col,
                              f'<c r="{col}{rnum}" s="{s}"><f>{f}</f></c>')

    dv = ('<dataValidations count="3">'
          '<dataValidation type="list" allowBlank="1" showInputMessage="1" '
          'showErrorMessage="1" sqref="G3:G22"><formula1>$Z$2:$Z$5</formula1></dataValidation>'
          '<dataValidation type="list" allowBlank="1" showInputMessage="1" '
          'showErrorMessage="1" sqref="H3:H22"><formula1>$AA$2:$AA$3</formula1></dataValidation>'
          '<dataValidation type="list" allowBlank="1" showInputMessage="1" '
          'showErrorMessage="1" sqref="I3:I22"><formula1>$AB$2:$AB$3</formula1></dataValidation>'
          '</dataValidations>')
    xml = xml.replace('<pageMargins', dv + '<pageMargins', 1)
    return xml


# -------------------------------------------------------------- piles sheet
# v12 groups columns by analysis type (Norm, 2026-07-09):
#   LEM-only (green):  H, qp, appl        FEM-only (blue): E, I, Area, Fixity
#   LEM & FEM (red):   D, S, Vcap, Mcap
# Final order: # Label x1 y1 x2 y2 H qp appl D S Vcap Mcap E I Area Fixity
PILES_PERM = {9: 10, 10: 11,            # D, S right one (appl lands at I)
              14: 12, 15: 13,           # Vcap, Mcap after S
              11: 14, 12: 15, 13: 16,   # E, I, Area after Mcap
              16: 17}                   # Fixity last


def piles_map(idx):
    return PILES_PERM.get(idx, idx)


def sort_row_cells(xml):
    """Re-sort cells within each row into ascending column order (required
    after a non-monotonic column permutation)."""
    def fix(m):
        toks = re.findall(r'<c r="[A-Z]{1,3}\d+"[^>]*?(?:/>|>.*?</c>)', m.group(2))
        toks.sort(key=lambda t: col_to_idx(re.match(r'<c r="([A-Z]{1,3})', t).group(1)))
        return m.group(1) + ''.join(toks) + m.group(3)
    return re.sub(r'(<row r="\d+"[^>]*>)(.*?)(</row>)', fix, xml)


def build_piles(xml):
    xml = shift_cell_refs(xml, piles_map)
    xml = sort_row_cells(xml)
    xml = re.sub(r'(sqref=)"([^"]+)"',
                 lambda m: f'{m.group(1)}"{" ".join(shift_sqref(p, piles_map) for p in m.group(2).split())}"',
                 xml)
    xml = shift_cols_block(xml, piles_map, {})

    # appl is LEM-only (active/passive factoring of H) -> clone the green
    # LEM-only header style from H (G2), NOT Fixity's blue FEM-only style
    s_hdr = get_style(xml, 'G2')
    xml = insert_cell(xml, 2, 'I', istr('I2', 'Appl', s_hdr))
    for ref, text in [('AA2', 'Active'), ('AA3', 'Passive')]:
        xml = insert_cell(xml, int(ref[2:]), 'AA', istr(ref, text, '0'))
    for rnum in range(3, 13):
        pm = re.search(rf'<c r="H{rnum}" s="(\d+)"', xml)
        if pm:
            xml = insert_cell(xml, rnum, 'I', styled(f'I{rnum}', pm.group(1)))
    dv = ('<dataValidation type="list" allowBlank="1" showInputMessage="1" '
          'showErrorMessage="1" sqref="I3:I12"><formula1>$AA$2:$AA$3</formula1>'
          '</dataValidation>')
    xml = re.sub(r'<dataValidations count="1">', '<dataValidations count="2">', xml)
    xml = xml.replace('</dataValidations>', dv + '</dataValidations>', 1)
    xml = re.sub(r'(<dimension ref=)"([^"]+)"', r'\1"A2:AA15"', xml)
    return xml


# -------------------------------------------------------------- new lloads
def build_lloads(reinf_xml):
    root = re.match(r'(<\?xml[^>]*\?>\s*)?<worksheet[^>]*>', reinf_xml).group(0)
    s_hdr = get_style(reinf_xml, 'C2')
    s_num = get_style(reinf_xml, 'A3')
    headers = ''.join(istr(f'{c}2', t, s_hdr) for c, t in
                      zip('ABCDEF', ['#', 'Label', 'x', 'y', 'P', 'Angle']))
    rows = [f'<row r="2" spans="1:6">{headers}</row>']
    for i in range(1, 21):
        r = i + 2
        cells = f'<c r="A{r}" s="{s_num}"><v>{i}</v></c>' + \
                ''.join(styled(f'{c}{r}', s_num) for c in 'BCDEF')
        rows.append(f'<row r="{r}" spans="1:6">{cells}</row>')
    return (root
            + '<dimension ref="A2:F22"/>'
            + '<sheetViews><sheetView workbookViewId="0"/></sheetViews>'
            + '<sheetFormatPr baseColWidth="10" defaultRowHeight="16"/>'
            + '<cols><col min="1" max="1" width="4.83203125" customWidth="1"/>'
            + '<col min="2" max="2" width="17.33203125" customWidth="1"/>'
            + '<col min="3" max="6" width="10.83203125" customWidth="1"/></cols>'
            + '<sheetData>' + ''.join(rows) + '</sheetData>'
            + '<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" '
              'header="0.3" footer="0.3"/></worksheet>')


def assert_cell_order(xml, name):
    """Excel refuses rows whose cells are not in ascending column order."""
    for rm in re.finditer(r'<row r="(\d+)"[^>]*>(.*?)</row>', xml):
        cols = [col_to_idx(c) for c in
                re.findall(r'<c r="([A-Z]{1,3})\d+"', rm.group(2))]
        if cols != sorted(cols) or len(cols) != len(set(cols)):
            raise RuntimeError(f'{name} row {rm.group(1)}: cell order broken: {cols}')


# ------------------------------------------------------------------- main
def build(src=SRC, out=OUT):
    zin = zipfile.ZipFile(src)
    parts = {}

    mat = zin.read(MAT_SHEET).decode('utf-8')
    parts[MAT_SHEET] = build_mat(mat)

    piles = zin.read(PILES_SHEET).decode('utf-8')
    s_green = get_style(piles, 'G2')          # LEM-only header green (H column)
    parts[PILES_SHEET] = build_piles(piles)

    reinf = zin.read(REINF_SHEET).decode('utf-8')
    parts[REINF_SHEET] = build_reinforce(reinf, s_green)

    main = zin.read(MAIN_SHEET).decode('utf-8')
    main2 = re.sub(r'(<c r="D5"[^>]*><v>)11(</v>)', r'\g<1>12\g<2>', main)
    if main2 == main:
        raise RuntimeError('main!D5 version cell (11) not found')
    parts[MAIN_SHEET] = main2

    chart = zin.read(CHART).decode('utf-8')
    parts[CHART] = re.sub(
        r'(reinforce!\$)([A-Z]{1,2})(\$)',
        lambda m: m.group(1) + idx_to_col(reinf_map(col_to_idx(m.group(2)))) + m.group(3),
        chart)

    wb = zin.read(WORKBOOK).decode('utf-8')
    # the Dir/Appl default formulas are new cells with no cached values or
    # calcChain entries — force a full recalc and drop calcChain (Excel
    # rebuilds it); see memory note project-xlsx-calcchain-recalc
    if 'fullCalcOnLoad' not in wb:
        if '<calcPr' in wb:
            wb = re.sub(r'<calcPr ', '<calcPr fullCalcOnLoad="1" ', wb, count=1)
        else:
            wb = wb.replace('</workbook>', '<calcPr fullCalcOnLoad="1"/></workbook>')
    sheet_ids = [int(s) for s in re.findall(r'sheetId="(\d+)"', wb)]
    new_sid = max(sheet_ids) + 1
    piles_entry = re.search(r'<sheet name="piles"[^/]*/>', wb).group(0)
    wb = wb.replace(piles_entry,
                    piles_entry + f'<sheet name="lloads" sheetId="{new_sid}" r:id="rId19"/>')
    parts[WORKBOOK] = wb

    rels = zin.read(WB_RELS).decode('utf-8')
    rels = re.sub(r'<Relationship[^>]*calcChain[^>]*/>', '', rels)  # part is dropped below
    rels = rels.replace('</Relationships>',
                        '<Relationship Id="rId19" '
                        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                        'Target="worksheets/sheet15.xml"/></Relationships>')
    parts[WB_RELS] = rels

    ct = zin.read(CONTENT_TYPES).decode('utf-8')
    ct = re.sub(r'<Override PartName="/xl/calcChain\.xml"[^/]*/>', '', ct)
    ct = ct.replace('</Types>',
                    '<Override PartName="/xl/worksheets/sheet15.xml" '
                    'ContentType="application/vnd.openxmlformats-officedocument.'
                    'spreadsheetml.worksheet+xml"/></Types>')
    parts[CONTENT_TYPES] = ct

    parts[LLOADS_SHEET] = build_lloads(reinf)

    for name in (MAT_SHEET, REINF_SHEET, PILES_SHEET, LLOADS_SHEET):
        assert_cell_order(parts[name], name)

    with zipfile.ZipFile(out, 'w', compression=zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename == 'xl/calcChain.xml':
                continue
            data = parts.pop(item.filename, None)
            zout.writestr(item, data.encode('utf-8') if isinstance(data, str)
                          else zin.read(item.filename))
        for name, data in parts.items():
            zout.writestr(name, data.encode('utf-8'))
    print(f'wrote {out}')


if __name__ == '__main__':
    out = OUT
    if '--out' in sys.argv:
        out = sys.argv[sys.argv.index('--out') + 1]
    build(out=out)
