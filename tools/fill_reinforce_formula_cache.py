"""Store cached results for the reinforce sheet's Dir/Appl formula cells.

The template's ``reinforce!H:I`` are VLOOKUPs reading ``Type``. A workbook
written outside Excel (``save_slope_data_to_xlsx``, openpyxl) carries the
formulas but no cached results, so Excel shows Dir/Appl once it opens the file
while every other reader — ``render_sheet``, a docs screenshot, a pandas dump —
sees blank cells. Excel itself stores a cached ``<v>`` beside every ``<f>``;
this script adds the same, by direct sheet-XML surgery (openpyxl cannot write a
value behind a formula), leaving the formulas byte-identical.
``fullCalcOnLoad`` still makes Excel recompute on open, so a wrong cache could
never survive a real session.

The cached values are not guessed: each row's Type is resolved through
``xslope.fileio``'s own preset table, and the loader's output is asserted
identical before and after.

The corpus builders call :func:`fill_reinforce_cache` on every workbook they
write, so a rebuilt file carries the results; on a file with nothing to fill
(no reinforce sheet, no typed row, every result already stored) it returns
without touching the file.

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


def _cell(sheet, ref):
    """The match for cell ``ref`` in the sheet XML, or None."""
    return re.search(r'(<c r="%s"[^>]*?)(/>|>(.*?)</c>)' % ref, sheet, re.S)


def _unfilled(sheet):
    """H/I formula cells with no stored result on rows whose Type (G) is set."""
    out = []
    for m in re.finditer(r'<c r="([HI])(\d+)"[^>]*?(?:/>|>(.*?)</c>)', sheet, re.S):
        row, body = int(m.group(2)), m.group(3) or ""
        if row < 3 or "<f" not in body or re.search(r"<v>[^<]+</v>", body):
            continue
        g = _cell(sheet, f"G{row}")
        gbody = (g.group(3) or "") if g else ""
        if re.search(r"<v>[^<]+</v>|<t[^>]*>[^<]+</t>", gbody):
            out.append(f"{m.group(1)}{row}")
    return out


def _reinforce_sheet(zin):
    """Archive name of the reinforce sheet, or None when the workbook has none."""
    wbxml = zin.read("xl/workbook.xml").decode("utf-8")
    m = re.search(r'<sheet[^>]*name="reinforce"[^>]*?r:id="rId(\d+)"', wbxml)
    if not m:
        return None
    rels = zin.read("xl/_rels/workbook.xml.rels").decode("utf-8")
    relel = re.search(r'<Relationship\b[^>]*\bId="rId%s"[^>]*/>' % m.group(1), rels)
    target = re.search(r'Target="([^"]+)"', relel.group(0)).group(1)
    target = target.lstrip("/")
    if not target.startswith("xl/"):
        target = "xl/" + target
    return target


def _shared_strings(zin):
    if "xl/sharedStrings.xml" not in zin.namelist():
        return []
    xml = zin.read("xl/sharedStrings.xml").decode("utf-8")
    return ["".join(re.findall(r"<t[^>]*>([^<]*)</t>", si))
            for si in re.findall(r"<si>(.*?)</si>", xml, re.S)]


def _type_of(sheet, row, shared):
    """The Type (G) text of ``row``: inline, formula-string or shared string."""
    g = _cell(sheet, f"G{row}")
    if not g:
        return ""
    head, body = g.group(1), g.group(3) or ""
    if 't="s"' in head:
        v = re.search(r"<v>(\d+)</v>", body)
        return shared[int(v.group(1))] if v else ""
    t = re.findall(r"<t[^>]*>([^<]*)</t>", body)
    if t:
        return "".join(t)
    v = re.search(r"<v>([^<]*)</v>", body)
    return v.group(1) if v else ""


def fill_reinforce_cache(path):
    """Cache the Dir/Appl results of ``path``'s reinforce sheet in place.

    Returns the cells filled. An empty list means the file needed nothing and
    was not rewritten, so a builder can call this after every write.
    """
    with zipfile.ZipFile(path) as zin:
        target = _reinforce_sheet(zin)
        if target is None or not _unfilled(zin.read(target).decode("utf-8")):
            return []

    before = load_lines(path)
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        shared = _shared_strings(zin)
        sheet = zin.read(target).decode("utf-8")
        filled = []
        for ref in _unfilled(sheet):
            col, row = ref[0], int(ref[1:])
            preset = PRESETS.get(_type_of(sheet, row, shared).strip().lower())
            if preset is None:  # a Type outside the lookup block: IFERROR's ""
                continue
            val = preset[0] if col == "H" else preset[1]
            mm = _cell(sheet, ref)
            cell = mm.group(0)
            head = mm.group(1)
            if 't="' not in head:
                head = head.rstrip() + ' t="str"'
            inner = re.search(r">(.*)</c>", cell, re.S)
            body = inner.group(1) if inner else ""
            body = re.sub(r"<v\s*/>|<v></v>", "", body)  # drop an empty cache
            sheet = sheet[:mm.start()] + head + ">" + body + f"<v>{val}</v></c>" \
                + sheet[mm.end():]
            filled.append(ref)

        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
            for n in names:
                data = sheet.encode("utf-8") if n == target else zin.read(n)
                zout.writestr(n, data)
    shutil.move(tmp, path)
    # The loader reads a stored result in place of the Type's preset, so an
    # unchanged Dir/Appl on every line is the proof that each cache is right.
    after = load_lines(path)
    assert before == after, f"loader output changed: {before} != {after}"
    return filled


def fill(path):
    filled = fill_reinforce_cache(path)
    if filled:
        print(f"{path}: cached {len(filled)} cells ({', '.join(filled)}); loader identical")
    else:
        print(f"{path}: nothing to cache")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        fill(p)
