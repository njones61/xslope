#!/usr/bin/env python3
"""Write the mesh each finite element lock was cut on into its test tag.

Every ``fem_ssrm`` and ``fem_tilt`` tag in the docs gets ``expected_elements`` and
``expected_nodes``: the counts of the mesh ``run_tests._fem_case_model`` builds for
that tag. ``run_tests.py`` then checks them on every run (``_check_mesh_counts``),
so a mesher change (a gmsh upgrade, a meshing fix) is reported as a moved mesh
before anything is solved, instead of reaching the lock only as a small shift in
the factor of safety that reads exactly like a solver change.

Only the two count fields are touched: an existing pair is replaced, a missing
one is appended, and nothing else on the line is re-rendered. Write the counts
when a lock is cut or re-cut, never to make a moved mesh pass: a count that
changed means the lock beside it was cut on a different mesh.

Usage:
    python tools/lock_mesh_counts.py                    # report, write nothing
    python tools/lock_mesh_counts.py --write            # write every tag
    python tools/lock_mesh_counts.py --page docs/verification/rs2.md --write
    python tools/lock_mesh_counts.py --jobs 3
"""
import argparse
import contextlib
import glob
import io
import os
import re
import sys
from multiprocessing import get_context

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

TYPES = ("fem_ssrm", "fem_tilt")
_TAG = re.compile(r"<!--\s*test:\s*(.*?)\s*-->")
_COUNT = re.compile(r",\s*expected_(elements|nodes)=[\d;]+")


def tags(pages):
    """Every FEM lock tag on ``pages`` as ``(page, line_no, test)``, where ``test``
    is the dict ``run_tests.parse_test_tags`` makes of that tag (an FEM tag is
    never expanded into several rows, so the page's FEM rows and its FEM tag lines
    pair up in order)."""
    import run_tests as RT
    out = []
    for page in pages:
        rows = [t for t in RT.parse_test_tags(page)
                if t.get("type") in TYPES and t.get("file")]
        with open(page) as fh:
            line_nos = [i for i, line in enumerate(fh, 1)
                        if (m := _TAG.search(line)) and _is_fem(m.group(1))]
        if len(rows) != len(line_nos):
            raise RuntimeError(f"{page}: {len(rows)} FEM rows but {len(line_nos)} tag lines")
        out += [(page, n, t) for n, t in zip(line_nos, rows)]
    return out


def _is_fem(tag_text):
    kv = dict(p.partition("=")[::2] for p in tag_text.split(","))
    kv = {k.strip(): v.strip() for k, v in kv.items()}
    return kv.get("type") in TYPES and bool(kv.get("file"))


def count(item):
    """``(page, line_no, elements, nodes, error)`` for one tag, meshed with its
    current counts removed so a stale pair cannot refuse the build."""
    import run_tests as RT
    page, line_no, test = item
    test = {k: v for k, v in test.items()
            if k not in ("expected_elements", "expected_nodes")}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            _sd, mesh = RT._fem_case_model(test)
    except Exception as exc:                       # reported, never written
        return page, line_no, None, None, f"{type(exc).__name__}: {exc}"
    return page, line_no, len(mesh["elements"]), len(mesh["nodes"]), None


def write_counts(page, line_no, elements, nodes):
    """Replace or append the two fields on one tag line, in place."""
    with open(page) as fh:
        lines = fh.readlines()
    line = _COUNT.sub("", lines[line_no - 1].rstrip("\n"))
    head, sep, tail = line.rpartition("-->")
    if not sep:
        raise ValueError(f"{page}:{line_no} is not a test tag line")
    lines[line_no - 1] = (f"{head.rstrip().rstrip(',')}, expected_elements={elements}, "
                          f"expected_nodes={nodes} {sep}{tail}\n")
    with open(page, "w") as fh:
        fh.writelines(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--page", action="append",
                    help="a docs page (repeatable); default: every page under docs/")
    ap.add_argument("--write", action="store_true", help="write the tags")
    ap.add_argument("--jobs", type=int, default=3)
    args = ap.parse_args(argv)

    os.chdir(_ROOT)
    pages = args.page or sorted(glob.glob("docs/**/*.md", recursive=True))
    items = tags(pages)
    print(f"{len(items)} FEM lock tag(s) on {len({p for p, _, _ in items})} page(s)")
    with get_context("spawn").Pool(max(1, args.jobs)) as pool:
        results = sorted(pool.map(count, items))

    changed = errors = 0
    for (page, line_no, el, nd, err), (_p, _l, test) in zip(
            results, sorted(items, key=lambda x: (x[0], x[1]))):
        if err:
            errors += 1
            print(f"  ERROR {page}:{line_no}  {err}")
            continue
        old = (test.get("expected_elements"), test.get("expected_nodes"))
        if all(got in (want if isinstance(want, tuple) else (want,))
               for got, want in zip((el, nd), old)):
            continue
        changed += 1
        was = "none" if old == (None, None) else f"{old[0]} elements, {old[1]} nodes"
        print(f"  {page}:{line_no}  {el} elements, {nd} nodes (was {was})")
        if args.write:
            write_counts(page, line_no, el, nd)
    print(f"{changed} tag(s) {'written' if args.write else 'to write'}, "
          f"{errors} could not be meshed")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
