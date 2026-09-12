"""Section order: the rows run down the page in the summary table's order.

A corpus page is a summary table followed by one section per problem, and the
table is what a reader navigates by. If the sections are in a different order
from the table, then following the table down the page means jumping about —
and the order drifts the moment a row is written last and appended at the end,
which is exactly how it happens.

So: **the anchored sections a page's summary table links must appear in the page
in the order the table lists them.** Nothing here says what that order should be
— the table decides, as it decides the dots and the locked values — this only
says the page agrees with it.

What the check reads
--------------------
* **Rows.** Every row of every ``corpus-summary`` block (or the marker the page
  config names) whose first cell links a LOCAL anchor. A row pointing at another
  page names a section that page owns, and a row with no link names no section;
  both are skipped.
* **Sections.** Only headings carrying an explicit ``{#anchor}`` that a row
  names. A heading the table never links — a methodology section, a shared
  discussion — may sit anywhere, so it is not read. An anchor the page defines
  only as an inline ``<a id=...>`` is skipped too: ``dots`` reports those, and
  an inline anchor has no section of its own to order.
* **Several rows naming one section.** A section covering three problems of a
  manual is ordered by the FIRST row that names it, which is where a reader
  meets it.
* **A row marked *covered*** is a cross-reference: the status term means another
  row owns the build, and that row is where the section belongs. So a covered row
  does not order anything. Without this, one page's row 35 — covered by its own
  Part IV problem 70 — would demand that problem 70's section be moved out of
  the Part IV block and into the middle of the first one.

The report names the pair that is out of order and both line numbers, so the fix
is a move rather than a hunt. ``--fix`` is deliberately not offered: moving a
section means moving its prose, its tag and its figure together, and a mechanical
reshuffle of a markdown file is how a caption ends up under the wrong image.

Usage
-----
    python -m tools.verification_checks.order
    python -m tools.verification_checks.order docs/verification/rs2.md
"""
import os
import re
import sys

from .dots import HEADING, INLINE_ANCHOR, ROW_LINK, SEPARATOR

#: The status term for a row whose problem is built under another row or page
#: (docs/verification/index.md#status-terms). Such a row links a section it does
#: not own, so it says nothing about where that section belongs.
COVERED = re.compile(r"\*covered\*", re.I)


def summary_anchors(lines, marker="corpus-summary"):
    """(line number, label, anchor) per summary row that links a local anchor,
    in the order the table lists them."""
    out, in_block, seen_header = [], False, False
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if not in_block:
            if s.startswith("<div") and marker in s:
                in_block, seen_header = True, False
            continue
        if s.startswith("</div"):
            in_block = False
            continue
        if not s.startswith("|"):
            continue
        if set(s) <= SEPARATOR:               # the |---|:-:| rule
            continue
        if not seen_header:                   # the header row
            seen_header = True
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        link = ROW_LINK.match(cells[0]) if cells else None
        if not link:
            continue
        target = link.group("target")
        if not target.startswith("#"):
            continue                          # a section another page owns
        if COVERED.search(cells[-1] if cells else ""):
            continue                          # a cross-reference, not the owner
        out.append((i, link.group("label"), target[1:]))
    return out


def section_lines(lines):
    """anchor -> line number, for every anchored heading."""
    out = {}
    for i, line in enumerate(lines, 1):
        m = HEADING.match(line)
        if m and m.group("anchor") not in out:
            out[m.group("anchor")] = i
    return out


def scan(path, cfg=None):
    """Read one page. Returns (problems, notes)."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    marker = getattr(cfg, "summary_marker", None) or "corpus-summary"
    rows = summary_anchors(lines, marker)
    if not rows:
        return [], ["no summary table with local anchors"]
    heads = section_lines(lines)
    inline = {m.group("anchor") for m in INLINE_ANCHOR.finditer("\n".join(lines))}

    wanted, first_row = [], {}
    for lineno, label, anchor in rows:
        if anchor in first_row:
            continue                          # several rows, one section
        first_row[anchor] = (lineno, label)
        if anchor in heads:
            wanted.append(anchor)

    notes = []
    skipped = [a for a in first_row if a not in heads]
    if skipped:
        notes.append(f"{len(skipped)} row anchor(s) with no heading of their own "
                     f"(inline or undefined): {', '.join(sorted(skipped)[:6])}")

    problems = []
    for k in range(1, len(wanted)):
        prev, cur = wanted[k - 1], wanted[k]
        if heads[cur] < heads[prev]:
            problems.append(
                f"section {{#{cur}}} (line {heads[cur]}, row "
                f"'{first_row[cur][1]}') comes before {{#{prev}}} "
                f"(line {heads[prev]}, row '{first_row[prev][1]}'), but the "
                f"summary table lists {first_row[prev][1]} first")
    return problems, notes


def run(path, cfg, report=print):
    """Check one page. Returns the failure count."""
    problems, notes = scan(path, cfg)
    for note in notes:
        report(f"  order     : note — {note}")
    if not problems:
        report("  order     : clean")
        return 0
    for p in problems:
        report(f"  order     : {os.path.basename(path)}: {p}")
    return len(problems)


def main(argv):
    from .pages import ORDER, PAGES
    here = os.path.dirname(os.path.abspath(__file__))
    pagedir = os.path.join(os.path.dirname(os.path.dirname(here)),
                           "docs", "verification")
    names = [a for a in argv[1:] if not a.startswith("-")] or list(ORDER)
    total = 0
    for n in names:
        p = n if n.endswith(".md") else os.path.join(pagedir, n + ".md")
        key = os.path.basename(p)[:-3]
        print(f"{key}:")
        total += run(p, PAGES.get(key))
    return 1 if total else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
