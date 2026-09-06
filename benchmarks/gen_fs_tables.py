"""Regenerate the per-method factor-of-safety tables and test tags in
docs/lem/samples.md.

A sample problem is *managed* by this script when its ``<!-- test: ... -->`` tag
is preceded by an ``<!-- fs-table -->`` / ``<!-- /fs-table -->`` block with
nothing but whitespace between the two. For a managed problem this runs every
applicable method and rewrites two things in place:

  * the visible "Factor of safety by method" table — only its data row; the
    caption, the header, the column alignment and any footnote are the page's
    and stay untouched; and
  * the ``fs_<method>=`` values inside the test tag, so run_tests.py checks
    every method. Every other key in the tag (``rapid``, ``benchmark``, the
    search-window keys, anything added later) is round-tripped verbatim.

Everything else on the page is passed through byte for byte: a tag with no table
above it is not a managed problem and is never touched, which is how the γ_sat
section, the right-facing reinforcement guard, the unreinforced Hassiotis case
and the second pile row keep their hand-written tags and their deliberate
absence of a table. To bring a new section under management, add an empty

    <!-- fs-table -->
    <!-- /fs-table -->

marker pair above its tag and rerun; the block is filled in.

**Held values.** run_tests.py checks a locked factor of safety to ±0.01, so any
recomputed value within 0.01 of the one on the page is the same lock. Those are
held at the page's printed digits rather than rewritten, and the drift is
reported at the end of the run. This keeps a rerun byte-identical when nothing
has meaningfully moved, so a real solver change shows up as a real diff instead
of being buried in third-decimal churn. ``--relock`` adopts the fresh values
regardless; that is the deliberate act of re-locking the page.

    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py            # rewrite samples.md
    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py --dry-run  # report, don't write
    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py --relock   # adopt drifted values
    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py --only piles  # one problem

Each method runs its own critical-surface search (the surfaces differ by method),
which is the meaningful comparison.
"""

import argparse
import contextlib
import io
import os
import re
import sys

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.search import (circular_search, noncircular_search, file_search_window,
                           noncircular_search_opts)
from xslope.solve import solve_selected

SAMPLES_MD = "docs/lem/samples.md"

#: run_tests.py's default absolute FS tolerance. A recomputed value inside this
#: band of the page's is the same lock, so the page's digits are held.
HOLD_TOL = 0.01

# Search-window keys a tag may carry. They say which mechanism the problem is
# about, so they are run WITH the search. (They are also written back out
# unchanged, along with every other key in the tag.)
WINDOW_KEYS = ('entry_range', 'exit_range', 'tangent_depth', 'center_box',
               'min_slip_depth')

# Tag types this script computes a table for. Anything else (reliability,
# gsat_pair, roundtrip, ...) is another checker's tag and is left alone.
FS_TYPES = ("single_circle", "circular_search", "noncircular_search")

# Display order; short tag name -> (solver name, column header).
METHODS = [
    ("oms", "oms", "OMS"),
    ("bishop", "bishop", "Bishop"),
    ("janbu", "janbu", "Janbu"),
    ("corps", "corps", "Corps"),
    ("lowe", "lowe", "Lowe"),
    ("spencer", "spencer", "Spencer"),
    ("mprice", "mprice", "M-P"),
]
# Methods that require a circular surface (skipped for non-circular searches).
CIRCULAR_ONLY = {"oms", "bishop"}

# (file, short-method) combinations that are method limitations rather than valid
# results: reported as "n/a" in the table (with a footnote) and omitted from the
# test tag. OMS breaks on the fully-submerged upstream face of earth_dam_up: its
# Fellenius base normal goes negative on the deepest slices under a full
# reservoir, so the strength it computes there is meaningless (see the OMS
# method note).
NA = "n/a"
EXCLUDED = {
    ("files/xslope_earth_dam_up.xlsx", "oms"),
}
EXCLUDED_NOTE = (
    "\\* OMS is not reported for this problem. Its base normal force is the "
    "Fellenius value $W\\cos\\alpha + D\\cos(\\alpha-\\beta) - u\\,\\Delta\\ell$, "
    "which under a full reservoir goes negative on the deepest slices — a "
    "quarter of them on the surface it settles on here — so the shear "
    "resistance it computes there is meaningless and the factor of safety it "
    "reports (0.886) sits far below every other method's. See the "
    "[OMS method note](oms.md)."
)
EXCLUDED_NOTE_MARK = "OMS is not reported for this problem"

TEST_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')
TABLE_RE = re.compile(r'[ \t]*<!--\s*fs-table\s*-->.*?<!--\s*/fs-table\s*-->\n?', re.DOTALL)
PAIR_RE = re.compile(r'(?P<key>[A-Za-z_][A-Za-z0-9_]*)\s*=\s*(?P<val>[^,]*)')
SEP_RE = re.compile(r'^\s*\|(?:\s*:?-{2,}:?\s*\|)+\s*$')
ROW_RE = re.compile(r'^\s*\|.*\|\s*$')


def parse_raw(tag_str):
    """Parse a tag body into a dict, keeping the file path relative (as written)."""
    params = {}
    for pair in tag_str.split(','):
        key, _, value = pair.partition('=')
        params[key.strip()] = value.strip()
    return params


def applicable(short, test_type):
    if test_type == "noncircular_search" and short in CIRCULAR_ONLY:
        return False
    return True


def compute_fs(params):
    """Run every applicable method for one problem; return {short: FS or None}."""
    file_rel = params["file"]
    test_type = params.get("type", "circular_search")
    num_slices = int(params.get("num_slices", 30))
    # rapid=true -> 3-stage rapid-drawdown analysis (works with any test type).
    rapid = str(params.get("rapid", "false")).strip().lower() in ("true", "1", "yes")
    path = os.path.join(os.path.dirname(SAMPLES_MD), file_rel)

    slope_data = load_slope_data(path)
    # The tag's own window wins; anything it leaves open comes from the model's
    # circles sheet, which is the same reading every other search path takes.
    win = {}
    for key in WINDOW_KEYS:
        raw = params.get(key)
        if raw is None or str(raw).strip() == '':
            continue
        vals = tuple(float(v) for v in str(raw).split(';') if v.strip() != '')
        win[key] = vals if len(vals) > 1 else vals[0]
    win.update(file_search_window(slope_data, already=win))
    results = {}
    for short, solver, _ in METHODS:
        if (file_rel, short) in EXCLUDED:
            results[short] = NA
            continue
        if not applicable(short, test_type):
            results[short] = None
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                if test_type == "single_circle":
                    circle = slope_data["circles"][0]
                    ok, res = generate_slices(slope_data, circle=circle, num_slices=num_slices)
                    if not ok:
                        results[short] = None
                        continue
                    slice_df, _ = res
                    r = solve_selected(solver, slice_df, rapid=rapid)
                    results[short] = r["FS"] if isinstance(r, dict) else None
                elif test_type == "noncircular_search":
                    fs_cache, _, _ = noncircular_search(slope_data, solver, num_slices=num_slices, rapid=rapid,
                                                        **noncircular_search_opts(win))
                    results[short] = fs_cache[0]["FS"] if fs_cache and fs_cache[0]["FS"] < 9999 else None
                else:  # circular_search
                    fs_cache, _, _, _ = circular_search(slope_data, solver, num_slices=num_slices, rapid=rapid,
                                                        **win)
                    results[short] = fs_cache[0]["FS"] if fs_cache and fs_cache[0]["FS"] < 9999 else None
        except Exception as e:
            print(f"    ! {solver} failed on {file_rel}: {e}", file=sys.stderr)
            results[short] = None
    return results


def resolve_display(fs, params, relock, drift, failures):
    """Decide the printed text for each method.

    Returns {short: text or None or NA}. A method whose freshly computed value
    is within HOLD_TOL of the value already in the tag keeps the tag's own text,
    so a rerun that has not moved a lock produces no diff; the difference is
    appended to `drift` for the end-of-run report.
    """
    display = {}
    for short, _, _ in METHODS:
        v = fs.get(short)
        if v is None and params.get(f"fs_{short}") is not None:
            # The method failed or found nothing on this run, but the page holds
            # a locked value for it. Keep the page's value and say so: blanking
            # the cell would leave the table contradicting its own tag and would
            # read as "this method does not apply here", which is a different
            # statement from "this run did not produce a number".
            failures.append((params["file"], short, params[f"fs_{short}"]))
            display[short] = params[f"fs_{short}"]
            continue
        if v is None or v == NA:
            display[short] = v
            continue
        fresh = f"{v:.3f}"
        page = params.get(f"fs_{short}")
        if page is None:
            display[short] = fresh
            continue
        try:
            old = float(page)
        except ValueError:
            display[short] = fresh
            continue
        if fresh != page:
            drift.append((params["file"], short, page, fresh, abs(v - old)))
        if not relock and abs(v - old) <= HOLD_TOL:
            display[short] = page          # same lock; hold the page's digits
        else:
            display[short] = fresh
    return display


def cell(display, short):
    """One table cell for a method."""
    v = display.get(short)
    if v == NA:
        return "n/a\\*"
    return v if v is not None else "—"


def data_row(display):
    return "| " + " | ".join(cell(display, s) for s, _, _ in METHODS) + " |"


def build_table(display):
    """A whole fs-table block, for a section that has only the empty markers."""
    headers = [hdr for _, _, hdr in METHODS]
    lines = [
        "<!-- fs-table -->",
        "**Factor of safety by method** (each method's own critical surface):",
        "",
        "| " + " | ".join(headers) + " |",
        "|" + "---:|" * len(headers),
        data_row(display),
    ]
    if any(display.get(s) == NA for s, _, _ in METHODS):
        lines += ["", EXCLUDED_NOTE]
    lines.append("<!-- /fs-table -->")
    return "\n".join(lines)


def rewrite_table(block, display):
    """Update an existing fs-table block in place.

    Only the data row is regenerated. The caption, the header row, the column
    alignment row and any footnote are the page's own text and are copied
    through unchanged — they carry hand-written statements ("each method on the
    same specified circle", "13.7 m row") that no solver run can reconstruct.
    A block with markers but no table yet is generated from scratch.
    """
    trailing = "\n" if block.endswith("\n") else ""
    lines = block[:len(block) - len(trailing)].split("\n")
    sep = next((i for i, ln in enumerate(lines)
                if SEP_RE.match(ln) and i and ROW_RE.match(lines[i - 1])), None)
    if sep is None or sep + 1 >= len(lines) or not ROW_RE.match(lines[sep + 1]):
        return build_table(display) + trailing
    lines[sep + 1] = data_row(display)
    if any(display.get(s) == NA for s, _, _ in METHODS) and \
            not any(EXCLUDED_NOTE_MARK in ln for ln in lines):
        lines[-1:-1] = ["", EXCLUDED_NOTE]
    return "\n".join(lines) + trailing


def rewrite_tag_body(body, display):
    """Update the fs_* values inside a tag body, leaving every other key alone.

    The body is edited by substitution rather than rebuilt from a whitelist of
    keys: a key this script does not know about (``rapid=true`` on the Johnson
    rapid-drawdown lock, say) has to survive the round trip. Rebuilding would
    silently drop it, and with it the analysis the lock is checking.
    """
    have = set()

    def sub(m):
        key = m.group("key")
        if not key.startswith("fs_"):
            return m.group(0)
        short = key[3:]
        text = display.get(short)
        if text is None or text == NA:
            return m.group(0)
        have.add(short)
        return m.group(0)[:m.start("val") - m.start()] + text

    new = PAIR_RE.sub(sub, body)
    # A method that now has a value but no key yet (a method added to METHODS
    # since the page was written) is appended after the last fs_ pair.
    missing = [(s, display[s]) for s, _, _ in METHODS
               if s not in have and isinstance(display.get(s), str) and display[s] != NA]
    if missing:
        ends = [m.end() for m in PAIR_RE.finditer(new) if m.group("key").startswith("fs_")]
        at = ends[-1] if ends else len(new)
        new = new[:at] + "".join(f", fs_{s}={t}" for s, t in missing) + new[at:]
    return new


def table_before(tables, pos, content):
    """The fs-table block that immediately precedes `pos` (whitespace only in
    between), or None. That adjacency is what marks a problem as managed."""
    for t in tables:
        if t.end() <= pos and content[t.end():pos].strip() == "":
            return t
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="report results, do not write")
    ap.add_argument("--relock", action="store_true",
                    help="adopt recomputed values even when they are within "
                         f"{HOLD_TOL} of the page's (a deliberate re-lock)")
    ap.add_argument("--only", action="append", default=[], metavar="SUBSTR",
                    help="only problems whose file name contains SUBSTR "
                         "(repeatable); every other section is passed through")
    args = ap.parse_args()

    with open(SAMPLES_MD) as f:
        content = f.read()

    tables = list(TABLE_RE.finditer(content))
    out = []
    last = 0
    drift = []
    failures = []
    n = 0
    for m in TEST_RE.finditer(content):
        params = parse_raw(m.group(1))
        if params.get("type") not in FS_TYPES:
            continue
        tbl = table_before(tables, m.start(), content)
        if tbl is None:
            # No table above it: a hand-written tag, not a managed problem.
            continue
        if args.only and not any(s in params.get("file", "") for s in args.only):
            continue
        # num_slices policy for a tag that does not state one: 50 for
        # verification benchmarks, 40 for ordinary samples. A tag that states
        # one is run at the number it states, and keeps it.
        if "num_slices" not in params:
            params["num_slices"] = "50" if "benchmark" in params else "40"
        print(f"  {params.get('file')}  ({params.get('type')}, ns={params['num_slices']})")
        fs = compute_fs(params)
        display = resolve_display(fs, params, args.relock, drift, failures)
        cells = "  ".join(f"{s}={fs[s]:.3f}" if isinstance(fs[s], float)
                          else f"{s}={fs[s] if fs[s] is not None else '—'}"
                          for s, _, _ in METHODS)
        print(f"      {cells}")
        n += 1

        out.append(content[last:tbl.start()])
        out.append(rewrite_table(tbl.group(0), display))
        out.append(content[tbl.end():m.start(1)])
        out.append(rewrite_tag_body(m.group(1), display))
        last = m.end(1)
    out.append(content[last:])
    new_content = "".join(out)

    print(f"\n{n} managed problem(s) recomputed")
    if failures:
        print("\nProduced no value this run; the page's locked value was kept:")
        for f_rel, short, page in failures:
            print(f"  {f_rel:48s} {short:8s} kept {page}")
    if drift:
        verb = "adopted" if args.relock else "held at the page's value"
        print(f"\nDrift against the page ({verb}; run_tests tolerance is {HOLD_TOL}):")
        for f_rel, short, page, fresh, d in drift:
            flag = "  <-- OUTSIDE TOLERANCE" if d > HOLD_TOL else ""
            print(f"  {f_rel:48s} {short:8s} page {page} -> {fresh}  (Δ {d:.4f}){flag}")
    else:
        print("No drift: every recomputed value prints the digits already on the page.")

    if args.dry_run:
        print("\n--- DRY RUN (not written) ---")
        return
    with open(SAMPLES_MD, "w") as f:
        f.write(new_content)
    print(f"\nWrote {SAMPLES_MD}")


if __name__ == "__main__":
    main()
