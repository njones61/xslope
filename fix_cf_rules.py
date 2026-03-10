"""Fix conditional formatting rules on 'mat' sheet and update version on 'main' sheet.

Edits the xlsx XML directly (zip archive) to avoid openpyxl's lossy load/save cycle
which corrupts rich text formatting (e.g. sigma characters).

Changes:
  - =$D9="cp" range E9:F50 -> F9:F50
  - =$D9="cp" range M9:N50 -> N9:N50
  - Remove =$D9="cp" rule on R9:V9 (if present)
  - Update version from 8 to 9 on 'main' sheet cell D5
  - Replace psi (ψ) in J8 with plain text "psi"
  - Replace s(ψ) in Q8 with rich text "σ(psi)" using mixed Symbol/Aptos fonts
"""

import sys
import zipfile
import shutil
import tempfile
import os
import re
from xml.etree import ElementTree as ET

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
NS_R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'


def get_sheet_file_map(zf):
    """Map sheet names to their xml file paths inside the xlsx zip."""
    wb_xml = ET.parse(zf.open('xl/workbook.xml'))
    rels_xml = ET.parse(zf.open('xl/_rels/workbook.xml.rels'))

    # Build rId -> target map
    rid_map = {}
    for rel in rels_xml.iter('{http://schemas.openxmlformats.org/package/2006/relationships}Relationship'):
        rid_map[rel.get('Id')] = 'xl/' + rel.get('Target')

    # Build sheet name -> file path map
    sheet_map = {}
    for sheet in wb_xml.iter(f'{{{NS}}}sheet'):
        name = sheet.get('name')
        rid = sheet.get(f'{{{NS_R}}}id')
        if rid in rid_map:
            sheet_map[name] = rid_map[rid]
    return sheet_map


# New shared string entries for psi fix
# "psi" - plain text for J8
SS_PSI = '<si><t>psi</t></si>'
# "σ(psi)" - rich text with Symbol font for σ( and ), Aptos for "psi" - for Q8
SS_SIGMA_PSI = (
    '<si>'
    '<r><rPr><sz val="12"/><color theme="0"/><rFont val="Symbol"/><charset val="2"/></rPr><t>s(</t></r>'
    '<r><rPr><sz val="12"/><color theme="0"/><rFont val="Aptos Narrow (Body)"/></rPr><t>psi</t></r>'
    '<r><rPr><sz val="12"/><color theme="0"/><rFont val="Symbol"/><charset val="2"/></rPr><t>)</t></r>'
    '</si>'
)


def fix_shared_strings(xml_bytes):
    """Add new shared string entries for psi and return new indices."""
    xml_str = xml_bytes.decode('utf-8')

    # Count existing strings to determine new indices
    existing_count = len(re.findall(r'<si[ >]', xml_str))
    psi_idx = existing_count
    sigma_psi_idx = existing_count + 1

    # Insert new <si> entries before closing </sst>
    new_entries = SS_PSI + SS_SIGMA_PSI
    xml_str = xml_str.replace('</sst>', new_entries + '</sst>')

    # Update count and uniqueCount attributes
    new_count = existing_count + 2
    xml_str = re.sub(r'(<sst[^>]*\bcount=")(\d+)(")', lambda m: f'{m.group(1)}{new_count}{m.group(3)}', xml_str)
    xml_str = re.sub(r'(<sst[^>]*\buniqueCount=")(\d+)(")', lambda m: f'{m.group(1)}{new_count}{m.group(3)}', xml_str)

    return xml_str.encode('utf-8'), psi_idx, sigma_psi_idx


def fix_mat_xml(xml_bytes, psi_idx=None, sigma_psi_idx=None):
    """Fix conditional formatting sqref values and psi cells in the mat sheet XML."""
    xml_str = xml_bytes.decode('utf-8')
    changes = []

    # Fix 1: E9:F50 -> F9:F50 (for cp rule)
    pattern1 = r'(<conditionalFormatting sqref=")E9:F50("><cfRule[^>]*><formula>\$D9="cp"</formula>)'
    if re.search(pattern1, xml_str):
        xml_str = re.sub(pattern1, r'\g<1>F9:F50\2', xml_str)
        changes.append('  [mat] E9:F50 -> F9:F50')

    # Fix 2: M9:N50 -> N9:N50 (for cp rule)
    pattern2 = r'(<conditionalFormatting sqref=")M9:N50("><cfRule[^>]*><formula>\$D9="cp"</formula>)'
    if re.search(pattern2, xml_str):
        xml_str = re.sub(pattern2, r'\g<1>N9:N50\2', xml_str)
        changes.append('  [mat] M9:N50 -> N9:N50')

    # Fix 3: Remove entire R9:V9 conditionalFormatting block (for cp rule)
    pattern3 = r'<conditionalFormatting sqref="R9:V9"><cfRule[^>]*><formula>\$D9="cp"</formula></cfRule></conditionalFormatting>'
    if re.search(pattern3, xml_str):
        xml_str = re.sub(pattern3, '', xml_str)
        changes.append('  [mat] removed R9:V9 rule')

    # Fix 4: Replace J8 cell to point to new "psi" shared string and change style
    # Style 26 uses Symbol font (renders "psi" as "πσι"), style 21 uses normal font
    if psi_idx is not None:
        pattern_j8 = r'<c r="J8" s="\d+" t="s"><v>\d+</v></c>'
        replacement_j8 = f'<c r="J8" s="21" t="s"><v>{psi_idx}</v></c>'
        if re.search(pattern_j8, xml_str):
            xml_str = re.sub(pattern_j8, replacement_j8, xml_str)
            changes.append(f'  [mat] J8 -> "psi" style=21 (shared string {psi_idx})')

    # Fix 5: Replace Q8 cell to point to new "σ(psi)" shared string
    if sigma_psi_idx is not None:
        pattern_q8 = r'(<c r="Q8"[^>]*t="s"><v>)\d+(</v></c>)'
        if re.search(pattern_q8, xml_str):
            xml_str = re.sub(pattern_q8, rf'\g<1>{sigma_psi_idx}\2', xml_str)
            changes.append(f'  [mat] Q8 -> "σ(psi)" (shared string {sigma_psi_idx})')

    return xml_str.encode('utf-8'), changes


def fix_main_xml(xml_bytes):
    """Update version from 8 to 9 in cell D5 of the main sheet."""
    xml_str = xml_bytes.decode('utf-8')
    changes = []

    # D5 cell with numeric value 8 -> 9
    pattern = r'(<c r="D5"[^>]*><v>)8(</v></c>)'
    if re.search(pattern, xml_str):
        xml_str = re.sub(pattern, r'\g<1>9\2', xml_str)
        changes.append('  [main] version 8 -> 9')
    elif '<c r="D5"' in xml_str and '>9<' in xml_str:
        changes.append('  [main] version already 9')
    else:
        changes.append('  [main] version cell not found or unexpected value')

    return xml_str.encode('utf-8'), changes


def fix_workbook(filepath, dry_run=False):
    with zipfile.ZipFile(filepath, 'r') as zf:
        sheet_map = get_sheet_file_map(zf)

        if 'mat' not in sheet_map:
            print('  [mat] sheet not found, skipping')
            return
        if 'main' not in sheet_map:
            print('  [main] sheet not found, skipping')
            return

        mat_path = sheet_map['mat']
        main_path = sheet_map['main']

        mat_xml = zf.read(mat_path)
        main_xml = zf.read(main_path)
        ss_xml = zf.read('xl/sharedStrings.xml')

        new_ss_xml, psi_idx, sigma_psi_idx = fix_shared_strings(ss_xml)
        new_mat_xml, mat_changes = fix_mat_xml(mat_xml, psi_idx, sigma_psi_idx)
        new_main_xml, main_changes = fix_main_xml(main_xml)

    all_changes = main_changes + mat_changes
    for c in all_changes:
        print(c)

    if not mat_changes and all(('already' in c or 'not found' in c) for c in main_changes):
        print('  No changes needed')
        return

    if dry_run:
        print(f'  [DRY RUN] Would save: {filepath}')
        return

    # Rewrite the zip, replacing modified sheets and shared strings
    tmp_fd, tmp_path = tempfile.mkstemp(suffix='.xlsx')
    os.close(tmp_fd)
    try:
        with zipfile.ZipFile(filepath, 'r') as zf_in:
            with zipfile.ZipFile(tmp_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf_out:
                for item in zf_in.infolist():
                    if item.filename == mat_path:
                        zf_out.writestr(item, new_mat_xml)
                    elif item.filename == main_path:
                        zf_out.writestr(item, new_main_xml)
                    elif item.filename == 'xl/sharedStrings.xml':
                        zf_out.writestr(item, new_ss_xml)
                    else:
                        zf_out.writestr(item, zf_in.read(item.filename))
        shutil.move(tmp_path, filepath)
        print(f'  Saved: {filepath}')
    except Exception:
        os.unlink(tmp_path)
        raise


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python fix_cf_rules.py <file1.xlsx> [file2.xlsx ...]")
        print("       python fix_cf_rules.py --dry-run <file1.xlsx> [...]")
        sys.exit(1)

    dry_run = '--dry-run' in sys.argv
    files = [f for f in sys.argv[1:] if f != '--dry-run']

    for f in files:
        print(f"\nProcessing: {f}")
        try:
            fix_workbook(f, dry_run=dry_run)
        except Exception as e:
            print(f"  ERROR: {e}")
