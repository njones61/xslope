#!/usr/bin/env python3
"""Write the two trials that define a strength-reduction lock into its test tag.

A ``solve_ssrm`` bisection ends on a bracket: the highest trial factor the model
stood at, and the lowest it failed at. The lock is that bracket's midpoint, so
the pair is what the lock IS — and a lock that is still reproducible is one whose
model still stands at the first and still fails at the second. ``run_tests.py``
can check it that way (``check=edges``) for two solves instead of nine, which is
what this tool supplies the input for: it reads each lock's persisted trial
record and appends ``f_stand``, ``f_fail`` and ``check=edges`` to the tag.

The record is the one the figure producers write beside each model
(``ssrm_run_record`` -> ``*_fem_meta.json``), and it is the ONLY source. A pair
invented from ``expected_fs`` and the tolerance would be a guess at where the
bisection closed, and a guess one bracket step out passes while the lock moves
underneath it. So a lock with no record stays in bracket mode and is listed here
as needing one; it gains its edges the next time it is re-cut, not by being
re-run for this.

Four things must hold before a pair is written. Both edges must be trials the
solve DECIDED — a trial that ran out of iteration budget (STABLE_STUCK,
AMBIGUOUS, INCONCLUSIVE, or simply at the ceiling) has not shown the model
standing or failing anywhere, and a bracket closed on one is a statement about
the budget. The pair must be the bracket the run ended on, not two trials from
the middle of it. It must straddle the lock. And it must be no wider than twice
the tag's tolerance, since the bisection stops inside the tolerance. A record
that fails any of these is reported with the reason and the tag is left alone —
which is what happens to a sidecar written before its lock was re-cut, where the
trials are a different model's.

``run_tests.py``'s ``lock_edges`` row re-checks the last three properties on
every ``check=edges`` tag in the docs, so a pair that stops matching its lock
fails the suite rather than quietly checking the wrong thing.

Usage:
    python tools/lock_edges.py                       # report, write nothing
    python tools/lock_edges.py --write               # write the tags
    python tools/lock_edges.py --page docs/verification/rs2.md --write
    python tools/lock_edges.py --benchmark RS2-48,RJ-18 --write
    python tools/lock_edges.py --missing              # list the locks with no record
    python tools/lock_edges.py --benchmark FEM-1-ssrm --recut   # cut one, then write it

``--recut`` is for the locks no figure producer draws: it re-solves the tag's own
bracket on the reference kernel, writes the record the run closed on, and then
writes the pair. Where a producer DOES draw the row, run the producer instead —
its run writes the figure and the record together, so the picture on the page and
the pair in the tag come off one bracket.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "tools"))

import ssrm_trial_audit as audit          # noqa: E402  (path set above)

#: Verdicts that mean a trial ran out of budget rather than deciding.
UNDECIDED = audit.UNDECIDED


def _fmt(value):
    """A trial factor as a tag writes it: enough digits to name the trial, none
    of the binary-fraction tail a JSON round trip leaves behind
    (1.2531250000000003 -> 1.253125)."""
    return f"{float(value):.10g}"


def _decided(trial, ceiling):
    """Did this trial answer the standing/failing question?

    ``ssrm_trial_audit.trial_decided`` is the single implementation: a verdict in
    the decided set, not inconclusive, and either inside its budget or carrying a
    corrector certification — an independent driver having reached equilibrium
    from that trial's own state answers the question whatever the sweep count.
    The two tools must agree, or a pair this one writes is a pair the audit reads
    as budget-bound.
    """
    return audit.trial_decided(trial, ceiling)


def edges_from_record(kv, meta):
    """``(f_stand, f_fail)`` for one lock, or ``(None, reason)``.

    ``kv`` is the parsed test tag and ``meta`` the path to its meta sidecar."""
    try:
        record = json.load(open(meta))
    except (OSError, ValueError):
        return None, "the sidecar is not readable"
    trials = record.get("trials")
    if not trials:
        return None, "the sidecar carries no trial record"
    interval = record.get("final_interval")
    if not interval or len(interval) != 2:
        return None, "the record carries no final bracket"

    try:
        stated = int(float(kv.get("max_iter", 12000)))
    except (TypeError, ValueError):
        stated = 12000
    try:
        declared = int(float(kv.get("max_iter_ceiling", audit.DEFAULT_CEILING)))
    except (TypeError, ValueError):
        declared = audit.DEFAULT_CEILING
    ceiling = max(declared, stated)

    decided = [t for t in trials if _decided(t, ceiling)]
    stands = [float(t["F"]) for t in decided
              if t.get("verdict") in ("CONVERGED", "JOINT_SETTLED")]
    fails = [float(t["F"]) for t in decided if t.get("verdict") == "FAILED"]
    if not stands:
        return None, "no trial was decided standing"
    if not fails:
        return None, "no trial was decided failing"
    f_stand, f_fail = max(stands), min(fails)

    # The pair must BE the bracket the run closed on. Where it is not, one of the
    # closing edges was a trial nothing could rule on, and the lock it produced is
    # a reading of the iteration budget (tools/ssrm_trial_audit.py names those).
    if (abs(f_stand - float(interval[0])) > 1e-9
            or abs(f_fail - float(interval[1])) > 1e-9):
        return None, (f"the final bracket [{_fmt(interval[0])}, {_fmt(interval[1])}] "
                      f"is not the decided pair [{_fmt(f_stand)}, {_fmt(f_fail)}] — "
                      f"an edge of it was never decided")

    expected = kv.get("expected_fs")
    if expected is None:
        return None, "the tag locks no factor of safety"
    expected = float(expected)
    if not (f_stand < expected <= f_fail):
        return None, (f"the record's bracket [{_fmt(f_stand)}, {_fmt(f_fail)}] does "
                      f"not straddle the lock {expected:g} — it is a different run "
                      f"of this model")
    tol = float(kv.get("tolerance", 0.05))
    if (f_fail - f_stand) > 2.0 * tol + 1e-12:
        return None, (f"the bracket is {f_fail - f_stand:.4g} wide, more than twice "
                      f"the tag's tolerance {tol:g}")
    return (f_stand, f_fail), None


def recut(page, kv, overrides, capture=False):
    """Re-cut one lock's bracket and persist the trial record it closed on.

    For the rows whose figures a producer draws, the producer's own run is the
    re-cut: it solves the bracket, writes the fields and the record together, and
    the figure on the page is then the mechanism that record describes. This is
    the path for the rest — a sweep point the page plots rather than figures, a
    row with no figure at all — where there is a lock to check and nothing that
    would redraw anything.

    The bracket is built by ``run_tests.build_fem_ssrm_case``, the suite's own
    tag-to-model mapping, so the trials recorded here are the trials the suite
    runs; and it is solved on the pure-NumPy reference kernel, the path a locked
    factor is defined by. The at-failure capture is skipped — it is an extra solve
    past the bracket, and nothing here draws a field.

    Returns ``(FS, record_path)``.
    """
    import run_tests as RT
    import xslope.fem as _fem
    from xslope.fem import solve_ssrm, ssrm_run_record

    def _key(tag):
        """What identifies one row among the rows of a page: the workbook, the
        discretization, the bracket and the lock. Two tags agreeing on all four
        are the same run asked for twice."""
        out = [os.path.basename(str(tag.get('file'))),
               str(tag.get('benchmark', '')), str(tag.get('element_type', ''))]
        for key in ('target_size', 'f_min', 'f_max', 'expected_fs'):
            value = tag.get(key)
            out.append('' if value in (None, '') else f"{float(value):.10g}")
        return tuple(out)

    tests = [t for t in RT.parse_test_tags(page) if t.get('type') == 'fem_ssrm']
    want = [t for t in tests if _key(t) == _key(kv)]
    if len(want) != 1:
        raise SystemExit(f"{kv.get('benchmark') or kv['file']}: the tag matched "
                         f"{len(want)} runnable rows on {os.path.basename(page)}")
    test = want[0]

    fem_data, kwargs, f_min, f_max, tol = RT.build_fem_ssrm_case(test)
    kwargs.setdefault('capture_failure_state', capture)
    with RT._force_fast_kernel(_fem, False):
        result = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=tol,
                            debug_level=0, **kwargs)
    if not result.get('converged'):
        raise SystemExit(f"{kv.get('benchmark') or kv['file']}: the bracket did "
                         f"not close: {result.get('error')}")

    record = ssrm_run_record(result, fem_data=fem_data, options={
        'tolerance': tol, 'F_min': f_min, 'F_max': f_max,
        'ssr_exclude': kwargs.get('ssr_exclude'),
    })
    meta = {'benchmark': kv.get('benchmark'), 'analysis': 'ssrm',
            'FS': float(result['FS']), 'expected_fs': kv.get('expected_fs'),
            'file': kv.get('file'), 'max_iter': int(float(kv.get('max_iter', 4000))),
            'record_only': True}
    meta.update(record)
    path = audit.row_meta_name(audit.stem_path(page, kv, overrides), kv)
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=2)
    return float(result['FS']), path


def write_edges(page, line_no, f_stand, f_fail):
    """Append the two fields and ``check=edges`` to one tag, in place.

    The tag's existing text is never re-parsed or re-rendered — the fields are
    appended to the line as it stands, so no other value on it can be rewritten
    by this tool. Returns the new line."""
    with open(page) as fh:
        lines = fh.readlines()
    line = lines[line_no - 1]
    head, sep, tail = line.rstrip("\n").rpartition("-->")
    if not sep:
        raise ValueError(f"{page}:{line_no} is not a test tag line")
    new = (f"{head.rstrip().rstrip(',')}, f_stand={_fmt(f_stand)}, "
           f"f_fail={_fmt(f_fail)}, check=edges {sep}{tail}\n")
    lines[line_no - 1] = new
    with open(page, "w") as fh:
        fh.writelines(lines)
    return new


def tags_with_lines(pages):
    """Every ``fem_ssrm`` tag on ``pages`` as ``(page, line_no, kv)``."""
    out = []
    for page in pages:
        with open(page) as fh:
            for i, line in enumerate(fh, 1):
                kv = audit._kv(line)
                if kv and kv.get("type") == "fem_ssrm" and kv.get("file"):
                    out.append((page, i, kv))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--page", action="append",
                    help="restrict to this page (repeatable)")
    ap.add_argument("--benchmark",
                    help="restrict to these benchmark ids (comma-separated)")
    ap.add_argument("--write", action="store_true",
                    help="write the tags (default: report only)")
    ap.add_argument("--missing", action="store_true",
                    help="list the locks with no usable trial record and exit")
    ap.add_argument("--recut", action="store_true",
                    help="re-solve the selected locks' brackets and persist the "
                         "trial record each closed on (for rows no figure "
                         "producer draws); implies --write")
    args = ap.parse_args(argv)
    if args.recut and not (args.benchmark or args.page):
        ap.error("--recut needs --benchmark or --page: a re-cut is minutes of "
                 "solving per row, and naming the rows is asking for them")

    pages = args.page or sorted(
        set(glob.glob(os.path.join(_ROOT, "docs", "**", "*.md"), recursive=True)))
    pages = [p if os.path.isabs(p) else os.path.join(_ROOT, p) for p in pages]
    wanted = None
    if args.benchmark:
        wanted = {b.strip().lower() for b in args.benchmark.split(",") if b.strip()}
    overrides = audit._sidecar_overrides()

    rows = tags_with_lines(pages)
    ready, already, blocked = [], [], []
    for page, line_no, kv in rows:
        name = kv.get("benchmark") or os.path.basename(kv["file"])
        if wanted is not None and str(name).lower() not in wanted:
            continue
        if str(kv.get("check", "")).strip().lower() == "edges":
            already.append((name, page, line_no))
            continue
        meta = audit.meta_path(page, kv, overrides)
        if args.recut and (wanted is not None or args.page):
            import time
            t0 = time.time()
            fs, meta = recut(page, kv, overrides)
            exp = float(kv["expected_fs"]) if kv.get("expected_fs") else None
            moved = f"  lock={exp:g} d={fs - exp:+.4f}" if exp is not None else ""
            print(f"recut {name:24s} FS={fs:.4f}{moved}  "
                  f"({time.time() - t0:.0f}s)  {os.path.relpath(meta, _ROOT)}",
                  flush=True)
        if not meta:
            blocked.append((name, page, "no meta sidecar beside the model"))
            continue
        pair, why = edges_from_record(kv, meta)
        if pair is None:
            blocked.append((name, page, why))
            continue
        ready.append((name, page, line_no, pair, kv))

    print(f"{len(rows)} fem_ssrm lock(s) on {len(set(p for p, _, _ in rows))} page(s)")
    if args.missing:
        print(f"\n{len(blocked)} lock(s) cannot be converted from what is on disk. "
              f"Each stays in bracket mode until its next re-cut, when the run that "
              f"cuts it writes a trial record:")
        for name, page, why in blocked:
            print(f"  {name:24s} {os.path.relpath(page, _ROOT):42s} {why}")
        return 0

    print(f"\n{len(ready)} lock(s) can be checked on their bracket edges:")
    for name, page, line_no, (f_stand, f_fail), kv in ready:
        print(f"  {name:24s} {os.path.relpath(page, _ROOT):42s} "
              f"stands at {_fmt(f_stand):>10s}, fails at {_fmt(f_fail):>10s} "
              f"(lock {float(kv['expected_fs']):g}, tolerance "
              f"{float(kv.get('tolerance', 0.05)):g})")
    if already:
        print(f"\n{len(already)} lock(s) already carry a pair:")
        for name, page, line_no in already:
            print(f"  {name:24s} {os.path.relpath(page, _ROOT)}:{line_no}")
    if blocked:
        print(f"\n{len(blocked)} lock(s) stay in bracket mode "
              f"(--missing lists them with reasons)")

    if (args.write or args.recut) and ready:
        # Highest line first, so a rewrite cannot move a line this pass has yet to
        # read. (It cannot today — the fields are appended in place — but a tool
        # that edits files by line number should not depend on that.)
        for name, page, line_no, pair, _kv in sorted(ready, key=lambda r: -r[2]):
            write_edges(page, line_no, *pair)
        print(f"\nwrote {len(ready)} tag(s). The pages that changed need "
              f"recertifying (tools/verification_checks/certify.py), and "
              f"run_tests.py's lock_edges row re-checks every pair.")
    elif ready:
        print("\n(--write to apply)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
