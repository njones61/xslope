"""The shipped documentation index is current.

The Studio assistant answers "does xslope support …" by searching
``xslope/resources/docs_index.json`` — the whole documentation, folded into one
resource so a pip install can read the pages it was built from. That resource is
generated (``tools/build_docs_index.py``) and committed, which means it can go
stale in exactly the way that matters: a page renamed, a heading retitled, a new
worksheet documented, a capability added — and the assistant keeps answering out
of the last build.

The defect this whole path exists to prevent was an answer given from memory: an
assistant asked whether xslope supports line loads said no, and offered a
workaround, on an input the template, the loader, the FEM and six documentation
pages all carry. An index a release forgot to rebuild puts it back in that
position on whatever was added since, so the guard is not optional bookkeeping.

This check regenerates the index in memory and compares it with the committed
copy. The generator is deterministic — everything sorted, nothing timestamped —
so a pass means the resource is exactly what a rebuild would write. A failure
names the pages and sections that moved, not just "the file differs", because
the fix (``python tools/build_docs_index.py``) is trivial and the useful
information is what changed.

The template is a source of this index and is never written by it: the check
reads ``docs/inputs/input_template.xlsx`` and compares, nothing more.
"""

from __future__ import annotations

import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

#: How many differing names a failure line spells out before it says "and N more".
_NAME_CAP = 6


def _names(items):
    items = sorted(items)
    if len(items) <= _NAME_CAP:
        return ", ".join(items)
    return ", ".join(items[:_NAME_CAP]) + f", and {len(items) - _NAME_CAP} more"


def _section_map(index):
    """``{'<url>#<anchor>': (heading, length)}`` for every section in an index."""
    out = {}
    for page in index.get("pages") or []:
        for sec in page.get("sections") or []:
            key = "%s#%s" % (page.get("url", ""), sec.get("anchor", ""))
            out[key] = (sec.get("heading", ""), int(sec.get("length") or 0))
    return out


def _diff(built, shipped):
    """The specific staleness, as failure lines."""
    out = []

    built_pages = {p["url"] for p in built.get("pages") or []}
    ship_pages = {p["url"] for p in shipped.get("pages") or []}
    if built_pages - ship_pages:
        out.append("pages missing from the shipped index: %s"
                   % _names(built_pages - ship_pages))
    if ship_pages - built_pages:
        out.append("pages in the shipped index that the documentation no longer "
                   "has: %s" % _names(ship_pages - built_pages))

    a, b = _section_map(built), _section_map(shipped)
    common_pages = built_pages & ship_pages
    added = {k for k in a if k not in b and k.split("#")[0] in common_pages}
    gone = {k for k in b if k not in a and k.split("#")[0] in common_pages}
    if added:
        out.append("sections missing from the shipped index: %s" % _names(added))
    if gone:
        out.append("sections in the shipped index that no page has: %s"
                   % _names(gone))
    changed = {k for k in set(a) & set(b) if a[k] != b[k]}
    if changed:
        out.append("sections whose heading or text changed: %s" % _names(changed))

    built_sheets = {s["name"]: s for s in
                    (built.get("template") or {}).get("sheets") or []}
    ship_sheets = {s["name"]: s for s in
                   (shipped.get("template") or {}).get("sheets") or []}
    if set(built_sheets) != set(ship_sheets):
        out.append("template sheets differ: index has %s, the template has %s"
                   % (_names(set(ship_sheets)), _names(set(built_sheets))))
    else:
        moved = {n for n in built_sheets if built_sheets[n] != ship_sheets[n]}
        if moved:
            out.append("template sheets whose columns or help text changed: %s"
                       % _names(moved))

    if not out:
        # Everything the diff knows how to name agrees, so what differs is the
        # packed body blob or a top-level field. Say that rather than pass.
        out.append("the index differs from a rebuild in its packed section text "
                   "or a top-level field")
    return out


def run():
    """Failures as a list of strings; empty when the shipped index is current."""
    import json

    from tools import build_docs_index as bdi

    if not bdi.INDEX_PATH.exists():
        return ["xslope/resources/docs_index.json is missing — run: "
                "python tools/build_docs_index.py"]

    buf = io.StringIO()
    with redirect_stdout(buf):
        built = bdi.build_index(verbose=False)
    text = bdi.serialize(built)

    shipped_text = bdi.INDEX_PATH.read_text(encoding="utf-8")
    if shipped_text == text:
        return []

    try:
        shipped = json.loads(shipped_text)
    except Exception as exc:
        return ["xslope/resources/docs_index.json does not parse (%r) — run: "
                "python tools/build_docs_index.py" % (exc,)]

    return [f + " — run: python tools/build_docs_index.py"
            for f in _diff(built, shipped)]


def main():
    print("documentation index sync:")
    failures = run()
    for f in failures:
        print("  - " + f)
    if failures:
        raise SystemExit(1)
    print("  the shipped index matches a rebuild.")


if __name__ == "__main__":
    main()
