"""Put every worksheet's rows and cells back in order, changing nothing else.

    python3 tools/fix_xlsx_cell_order.py FILE_OR_DIR [...]          # report what it would fix
    python3 tools/fix_xlsx_cell_order.py --write FILE_OR_DIR [...]  # fix in place

Excel requires each sheet's rows in ascending order and each row's cells in column order; a file that breaks
either opens with "We found a problem with some content" and an offer to repair it. XSLOPE's writer broke the
second from template v19 (release 0.2.1) until it was fixed in 4c88cb8b: the main sheet's run options
(D14:D22) followed the cells in columns F and G. load_slope_data reads such files correctly, so no result was
affected.

This tool sorts the cells of any row that is out of order (and the rows of any sheet that is), and copies every
other part of the workbook unchanged: the same entries, in the same order, with their timestamps and
compression method (the data is compressed again, so the compressed bytes can differ). A row holding anything
besides cells stops it with an error instead of a guess. Re-saving a file in XSLOPE Studio repairs it too.
"""
import re
import sys
import zipfile
from pathlib import Path

CELL = re.compile(r'<c\s[^>]*?/>|<c\s[^>]*?>.*?</c>', re.DOTALL)
ROW = re.compile(r'<row\s+r="(\d+)"[^>]*?(?:/>|>.*?</row>)', re.DOTALL)


def col(cell):
    n = 0
    for ch in re.match(r'<c\s+r="([A-Z]+)', cell).group(1):
        n = n * 26 + ord(ch) - 64
    return n


def fix_row(row_xml):
    m = re.match(r'(<row\b[^>]*>)(.*)(</row>)$', row_xml, re.DOTALL)
    if not m:
        return row_xml
    head, body, tail = m.groups()
    cells = CELL.findall(body)
    if [col(c) for c in cells] == sorted(col(c) for c in cells):
        return row_xml
    rest = CELL.sub("", body)
    if rest.strip():
        raise ValueError("row holds something besides cells: %r" % rest[:80])
    return head + "".join(sorted(cells, key=col)) + tail


def fix_sheet(xml_text):
    rows = list(ROW.finditer(xml_text))
    if not rows:
        return xml_text
    fixed = [fix_row(r.group(0)) for r in rows]
    nums = [int(r.group(1)) for r in rows]
    if nums != sorted(nums):
        fixed = [f for _, f in sorted(zip(nums, fixed), key=lambda p: p[0])]
    start, end = rows[0].start(), rows[-1].end()
    between = xml_text[start:end]
    if ROW.sub("", between).strip():
        raise ValueError("sheetData holds something between rows")
    return xml_text[:start] + "".join(fixed) + xml_text[end:]


def fix(path, write):
    with zipfile.ZipFile(path) as z:
        infos = z.infolist()
        data = {i.filename: z.read(i.filename) for i in infos}
    changed = []
    for name in data:
        if re.match(r"xl/worksheets/sheet\d+\.xml$", name):
            text = data[name].decode("utf-8")
            new = fix_sheet(text)
            if new != text:
                data[name] = new.encode("utf-8")
                changed.append(name)
    if changed and write:
        tmp = path.with_suffix(".tmp_fix.xlsx")
        with zipfile.ZipFile(tmp, "w") as out:
            for i in infos:
                out.writestr(i, data[i.filename], compress_type=i.compress_type)
        tmp.replace(path)
    return changed


def main():
    write = "--write" in sys.argv
    files = []
    for a in (x for x in sys.argv[1:] if x != "--write"):
        p = Path(a)
        files += [p] if p.is_file() else sorted(f for f in p.rglob("*.xlsx") if not f.name.startswith("~$"))
    n = 0
    for f in files:
        ch = fix(f, write)
        if ch:
            n += 1
            print(("fixed " if write else "would fix ") + str(f), [c.split("/")[-1] for c in ch])
    print(f"{n} of {len(files)} workbooks {'fixed' if write else 'to fix'}")


if __name__ == "__main__":
    main()
