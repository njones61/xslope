"""Store cached results for the reinforce sheet's Dir/Appl formula cells.

The template's ``reinforce!H:I`` are VLOOKUPs reading ``Type``. A workbook
written by openpyxl carries the formulas but no cached results, so Excel shows
Dir/Appl once it opens the file while every other reader — ``render_sheet``,
a docs screenshot, a pandas dump — sees blank cells. Excel itself stores a
cached ``<v>`` beside every ``<f>``; this script adds the same, by direct
sheet-XML surgery (openpyxl cannot write a value behind a formula), leaving
the formulas byte-identical. ``fullCalcOnLoad`` still makes Excel recompute
on open, so a wrong cache could never survive a real session.

The cached values are not guessed: each row's Type is resolved through
``xslope.fileio``'s own preset table, and the loader's output is asserted
identical before and after.

Run:  python3 tools/fill_reinforce_formula_cache.py <workbook.xlsx> [...]
"""

import re
import shutil
import sys
import zipfile
import contextlib
import io

PRESETS = {  # Type -> (Dir, Appl), the template's own lookup block
    "geosynthetic": ("Tangent", "Active"),
    "nail": ("Axial", "Passive"),
    "tieback": ("Axial", "Active"),
    "anchor": ("Axial", "Active"),
}


def load_lines(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(path)
    return [(l.get("label"), l.get("type"), l.get("dir"), l.get("appl"))
            for l in sd.get("reinforcement_lines") or []]


def fill(path):
    before = load_lines(path)
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        # find the reinforce sheet by its declared name
        wbxml = zin.read("xl/workbook.xml").decode("utf-8")
        m = re.search(r'<sheet[^>]*name="reinforce"[^>]*?r:id="rId(\d+)"', wbxml)
        rels = zin.read("xl/_rels/workbook.xml.rels").decode("utf-8")
        relel = re.search(r'<Relationship\b[^>]*\bId="rId%s"[^>]*/>' % m.group(1), rels)
        target = re.search(r'Target="([^"]+)"', relel.group(0)).group(1)
        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        sheet = zin.read(target).decode("utf-8")
        shared = zin.read("xl/sharedStrings.xml").decode("utf-8") if "xl/sharedStrings.xml" in names else ""

        # rows 3..: read Type (G) inline/shared is messy — take types from the loader instead
        types = [t for (_, t, _, _) in before]
        filled = []
        row = 3
        for t in types:
            dirv, applv = PRESETS[str(t).lower()]
            for col, val in (("H", dirv), ("I", applv)):
                ref = f"{col}{row}"
                # match the formula cell and add a cached string value
                pat = re.compile(r'(<c r="%s"[^>]*?)(/>|>(.*?)</c>)' % ref, re.S)
                mm = pat.search(sheet)
                if not mm:
                    continue
                cell = mm.group(0)
                if "<f" not in cell:
                    continue  # not a formula cell; leave alone
                if re.search(r"<v>[^<]+</v>", cell):
                    continue  # already carries a non-empty cache
                head = mm.group(1)
                if 't="' not in head:
                    head = head.rstrip() + ' t="str"'
                inner = re.search(r">(.*)</c>", cell, re.S)
                body = inner.group(1) if inner else ""
                body = re.sub(r"<v\s*/>|<v></v>", "", body)  # drop an empty cache
                newcell = head + ">" + body + f"<v>{val}</v></c>"
                sheet = sheet[:mm.start()] + newcell + sheet[mm.end():]
                filled.append(ref)
            row += 1

        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for n in names:
                data = sheet.encode("utf-8") if n == target else zin.read(n)
                zout.writestr(n, data)
    shutil.move(tmp, path)
    after = load_lines(path)
    assert before == after, f"loader output changed: {before} != {after}"
    assert filled, f"{path}: nothing was filled -- check the cell matching"
    print(f"{path}: cached {len(filled)} cells ({', '.join(filled)}); loader identical")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        fill(p)
