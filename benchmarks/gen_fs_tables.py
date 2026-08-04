"""Regenerate the per-method factor-of-safety tables and test tags in
docs/lem/samples.md.

For every LEM sample problem (identified by its ``<!-- test: ... -->`` tag) this
runs all applicable methods and rewrites two things in place:

  * the test tag, as a compact multi-method form
    ``<!-- test: file=..., type=..., num_slices=..., fs_oms=..., fs_bishop=..., ... -->``
    so run_tests.py checks every method, not just one; and
  * a visible "Factor of safety by method" table (wrapped in
    ``<!-- fs-table -->`` / ``<!-- /fs-table -->`` markers) immediately above the
    tag, so readers can compare methods at a glance.

Each method runs its own critical-surface search (the surfaces differ by method),
which is the meaningful comparison. The operation is idempotent: rerun any time a
solver changes.

    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py            # rewrite samples.md
    PYTHONPATH=. python3 benchmarks/gen_fs_tables.py --dry-run  # print, don't write
"""

import argparse
import contextlib
import io
import os
import re
import sys

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.search import circular_search, noncircular_search
from xslope.solve import solve_selected

SAMPLES_MD = "docs/lem/samples.md"

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

TEST_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')
TABLE_RE = re.compile(r'[ \t]*<!--\s*fs-table\s*-->.*?<!--\s*/fs-table\s*-->\n?', re.DOTALL)


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
                    fs_cache, _, _ = noncircular_search(slope_data, solver, num_slices=num_slices, rapid=rapid)
                    results[short] = fs_cache[0]["FS"] if fs_cache and fs_cache[0]["FS"] < 9999 else None
                else:  # circular_search
                    fs_cache, _, _, _ = circular_search(slope_data, solver, num_slices=num_slices, rapid=rapid)
                    results[short] = fs_cache[0]["FS"] if fs_cache and fs_cache[0]["FS"] < 9999 else None
        except Exception as e:
            print(f"    ! {solver} failed on {file_rel}: {e}", file=sys.stderr)
            results[short] = None
    return results


def build_table(fs):
    headers = [hdr for _, _, hdr in METHODS]
    def cell(short):
        v = fs.get(short)
        if v == NA:
            return "n/a\\*"
        return f"{v:.3f}" if v is not None else "—"
    row = [cell(short) for short, _, _ in METHODS]
    lines = [
        "<!-- fs-table -->",
        "**Factor of safety by method** (each method's own critical surface):",
        "",
        "| " + " | ".join(headers) + " |",
        "|" + "---:|" * len(headers),
        "| " + " | ".join(row) + " |",
    ]
    if any(fs.get(s) == NA for s, _, _ in METHODS):
        lines += ["", EXCLUDED_NOTE]
    lines.append("<!-- /fs-table -->")
    return "\n".join(lines)


def build_tag(params, fs):
    # Preserve identity keys in a stable order; replace any method/expected_fs/fs_* with fresh fs_*.
    parts = [f"file={params['file']}", f"type={params.get('type', 'circular_search')}"]
    if "num_slices" in params:
        parts.append(f"num_slices={params['num_slices']}")
    for short, _, _ in METHODS:
        v = fs.get(short)
        if v is not None and v != NA:
            parts.append(f"fs_{short}={v:.3f}")
    if "benchmark" in params:
        parts.append(f"benchmark={params['benchmark']}")
    return "<!-- test: " + ", ".join(parts) + " -->"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="print results, do not write")
    args = ap.parse_args()

    with open(SAMPLES_MD) as f:
        content = f.read()

    # 1) Strip any existing fs-table blocks so we regenerate cleanly.
    content = TABLE_RE.sub("", content)

    # 2) Walk the test tags in order; rewrite the first tag for each (file,type)
    #    and drop later duplicates (e.g. a file that had one tag per method).
    out = []
    last = 0
    seen = {}
    for m in TEST_RE.finditer(content):
        params = parse_raw(m.group(1))
        key = (params.get("file"), params.get("type"))
        out.append(content[last:m.start()])
        # Only FS-by-method tags get a generated table; pass everything else
        # through untouched (e.g. type=reliability carries expected_beta, not fs_*).
        if params.get("type") not in ("single_circle", "circular_search", "noncircular_search"):
            out.append(m.group(0))
            last = m.end()
            continue
        if key in seen:
            print(f"  (collapsed duplicate tag for {key[0]})")
            last = m.end()
            # also swallow a trailing newline left by the removed tag
            if content[last:last + 1] == "\n":
                last += 1
            continue
        # num_slices policy: 50 for verification benchmarks, 40 for ordinary samples.
        params['num_slices'] = 50 if 'benchmark' in params else 40
        print(f"  {params.get('file')}  ({params.get('type', 'circular_search')}, ns={params['num_slices']})")
        fs = compute_fs(params)
        seen[key] = fs
        cells = "  ".join(f"{s}={fs[s]:.3f}" if isinstance(fs[s], float)
                          else f"{s}={fs[s] if fs[s] is not None else '—'}"
                          for s, _, _ in METHODS)
        print(f"      {cells}")
        out.append(build_table(fs) + "\n\n" + build_tag(params, fs))
        last = m.end()
    out.append(content[last:])
    new_content = "".join(out)

    # Tidy: collapse any 3+ blank lines left by removals down to 2.
    new_content = re.sub(r"\n{3,}", "\n\n", new_content)

    if args.dry_run:
        print("\n--- DRY RUN (not written) ---")
        return
    with open(SAMPLES_MD, "w") as f:
        f.write(new_content)
    print(f"\nWrote {SAMPLES_MD}")


if __name__ == "__main__":
    main()
