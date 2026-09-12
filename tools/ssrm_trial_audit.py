#!/usr/bin/env python3
"""Which strength-reduction locks were cut at the sweep ceiling rather than at a
mechanism.

A ``solve_ssrm`` bracket halves an interval by asking, of each trial factor,
whether the model stands there. A trial that reaches its iteration budget without
deciding is recorded ``STABLE_STUCK`` or ``AMBIGUOUS`` and the bracket treats it
as not standing — so a model whose equilibria are slow to converge is cut LOW,
and the factor of safety the row locks is a statement about the budget rather
than about the slope. The geotextile wall family was found that way: four of
RS2-52's nine trials ended at the ceiling, both edges of its final bracket among
them, and lifting the budget moved the row by ten percent.

This reads the trial record every lock's meta sidecar carries
(``ssrm_run_record``, persisted by ``export_fem_solution``) and reports, per lock:
how many trials decided, how many ended at the ceiling, and how close the longest
DECIDED trial came to it. A lock with no trial at the ceiling was decided by the
model; one with a trial at the ceiling on either edge of its final bracket was
decided by the budget, and its value is not yet a measurement.

**The effective ceiling is not the tag's number.** ``solve_fem`` extends a budget
that is still making progress, up to ``max_iterations_ceiling`` (default 50 000),
and takes ``ceiling = max(ceiling, max_iterations)`` — so a tag saying
``max_iter=30000`` actually runs to 50 000, and a tag above 50 000 with no ceiling
of its own has a ceiling EQUAL to its budget and cannot extend at all. A tag
states its headroom with ``max_iter_ceiling``. That is what this reports as the
effective ceiling.

Most committed sidecars predate the trial record and carry none; those rows are
listed as NOT CAPTURED, which is a thing to fix by re-running them, not a pass.

It also reads the other thing a tag can say about the same bracket. A row tagged
``check=edges`` is verified by re-solving at ``f_stand`` and ``f_fail`` rather
than by bisecting, so those two factors have to BE the bracket its record ended
on; a pair that has drifted from the record is reported here, because such a row
passes while checking two factors that belong to a different run.
``tools/lock_edges.py`` is what writes them.

Usage:
    python tools/ssrm_trial_audit.py [--page docs/verification/rs2.md] [--all]

Without ``--all`` only the rows that have something to report are printed.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TAG = re.compile(r"<!--\s*test:\s*(.*?)\s*-->")

#: The solver's own default hard stop on budget extension
#: (``fem.solve_fem(max_iterations_ceiling=50000)``).
DEFAULT_CEILING = 50000

#: Verdicts that mean the trial ran out of budget rather than deciding.
UNDECIDED = ("STABLE_STUCK", "AMBIGUOUS", "INCONCLUSIVE")

#: Verdicts that ANSWER the standing/failing question, which is not the same list
#: as "converged or failed". ``JOINT_SETTLED`` is a decision on a jointed trial
#: whose slip, displacement field and soil residual have all settled and whose
#: joint residual has stopped falling into a limit cycle no budget brings down
#: (fem.joint_verdict); it never met the force tolerance and never claims to, but
#: it stands, and an edge that reads it is answered. run_tests.DECIDED_VERDICTS
#: is the same list, and the two must not disagree.
DECIDED = ("CONVERGED", "FAILED", "JOINT_SETTLED")


def _kv(line):
    m = TAG.search(line)
    if not m:
        return None
    out = {}
    for part in m.group(1).split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _sidecar_overrides():
    """``{benchmark: stem}`` from the RS2 figure producer, where a row's field is
    written under a stem that is not its workbook's."""
    try:
        sys.path.insert(0, os.path.join(_ROOT, "benchmarks", "rocscience"))
        import make_rs2_figures as F
        return dict(F.SIDECAR_STEM)
    except Exception:
        return {}


def tags(pages):
    """Every ``fem_ssrm`` tag on ``pages``, as ``(page, kv)``."""
    out = []
    for page in pages:
        for line in open(page):
            kv = _kv(line)
            if kv and kv.get("type") == "fem_ssrm" and kv.get("file"):
                out.append((page, kv))
    return out


def meta_path(page, kv, overrides):
    """The meta sidecar a row's run writes, or None where the path does not
    resolve."""
    stem = overrides.get(kv.get("benchmark", ""))
    book = os.path.normpath(os.path.join(os.path.dirname(page), kv["file"]))
    if stem:
        path = os.path.join(os.path.dirname(book), stem)
    else:
        path = os.path.splitext(book)[0]
    meta = path + "_fem_meta.json"
    return meta if os.path.exists(meta) else None


def audit_one(kv, meta):
    """What one lock's trial record says. Returns a dict, or None with no record."""
    try:
        rec = json.load(open(meta))
    except Exception:
        return None
    trials = rec.get("trials")
    if not trials:
        return None
    try:
        stated = int(float(kv.get("max_iter", 12000)))
    except (TypeError, ValueError):
        stated = 12000
    try:
        declared = int(float(kv.get("max_iter_ceiling", DEFAULT_CEILING)))
    except (TypeError, ValueError):
        declared = DEFAULT_CEILING
    ceiling = max(declared, stated)
    rows = []
    for t in trials:
        try:
            rows.append((float(t.get("F")), str(t.get("verdict")),
                         int(t.get("iterations", 0))))
        except (TypeError, ValueError):
            continue
    if not rows:
        return None
    undecided = [r for r in rows if r[1] in UNDECIDED or r[2] >= ceiling]
    decided = [r for r in rows if r not in undecided]
    # The final bracket's two edges: the highest standing trial and the lowest
    # refused one. A lock is budget-bound when either edge never decided.
    stand = [r for r in rows if r[1] == "CONVERGED"]
    refuse = [r for r in rows if r[1] not in ("CONVERGED",) ]
    edge_lo = max((r[0] for r in stand), default=None)
    edge_hi = min((r[0] for r in refuse), default=None)
    edge_undecided = sum(1 for r in undecided
                         if edge_lo is not None and abs(r[0] - edge_lo) < 1e-12
                         or edge_hi is not None and abs(r[0] - edge_hi) < 1e-12)
    return {
        "edges": (edge_lo, edge_hi),
        "trials": len(rows),
        "decided": len(decided),
        "at_ceiling": len(undecided),
        "ceiling": ceiling,
        "stated": stated,
        "max_decided": max((r[2] for r in decided), default=0),
        "max_any": max(r[2] for r in rows),
        "edge_undecided": edge_undecided,
        "interval": rec.get("final_interval"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--page", action="append",
                    help="restrict to this page (repeatable)")
    ap.add_argument("--all", action="store_true",
                    help="print every lock, not only the ones with a finding")
    args = ap.parse_args(argv)

    pages = args.page or sorted(
        set(glob.glob(os.path.join(_ROOT, "docs", "**", "*.md"), recursive=True)))
    pages = [p if os.path.isabs(p) else os.path.join(_ROOT, p) for p in pages]
    overrides = _sidecar_overrides()

    rows = tags(pages)
    print(f"{len(rows)} fem_ssrm locks on {len(set(p for p, _ in rows))} page(s)")
    flagged, captured, missing = [], 0, []
    checked_on_edges, stale_edges = 0, []
    for page, kv in rows:
        meta = meta_path(page, kv, overrides)
        a = audit_one(kv, meta) if meta else None
        name = kv.get("benchmark") or os.path.basename(kv["file"])
        # What the TAG says about the same bracket. A row checked on its edges is
        # re-solved at exactly these two factors instead of bisected, so a pair
        # that has drifted from the record it was written off is a row checking
        # something other than its own lock.
        tag_edges = None
        if str(kv.get("check", "")).strip().lower() == "edges":
            checked_on_edges += 1
            try:
                tag_edges = (float(kv["f_stand"]), float(kv["f_fail"]))
            except (KeyError, TypeError, ValueError):
                tag_edges = None
        if a is None:
            missing.append((name, os.path.relpath(page, _ROOT)))
            continue
        captured += 1
        finding = a["at_ceiling"] > 0
        if finding:
            flagged.append((name, a))
        rec_lo, rec_hi = a["edges"]
        if (tag_edges and rec_lo is not None and rec_hi is not None
                and (abs(tag_edges[0] - rec_lo) > 1e-9
                     or abs(tag_edges[1] - rec_hi) > 1e-9)):
            stale_edges.append((name, tag_edges, (rec_lo, rec_hi)))
        if args.all or finding:
            head = "BUDGET-BOUND" if a["edge_undecided"] else \
                   ("at-ceiling trials" if finding else "decided")
            edges = (f", checked on edges [{tag_edges[0]:g}, {tag_edges[1]:g}]"
                     if tag_edges else "")
            print(f"  {name:22s} {head:18s} "
                  f"{a['decided']}/{a['trials']} decided, "
                  f"{a['at_ceiling']} at the {a['ceiling']} ceiling "
                  f"(tag max_iter={a['stated']}), "
                  f"longest decided {a['max_decided']}{edges}")

    print(f"\n{captured} lock(s) carry a trial record; {len(flagged)} have a trial "
          f"that ended at the ceiling")
    print(f"{checked_on_edges} lock(s) are checked on their bracket edges rather "
          f"than re-bisected (tools/lock_edges.py writes the pair)")
    if stale_edges:
        print(f"{len(stale_edges)} of them carry a pair that is NOT the bracket "
              f"their own record ended on — the row is checking two factors that "
              f"belong to some other run, and the pair has to be rewritten:")
        for name, tag_pair, rec_pair in stale_edges:
            print(f"  {name:22s} tag [{tag_pair[0]:g}, {tag_pair[1]:g}] vs "
                  f"record [{rec_pair[0]:g}, {rec_pair[1]:g}]")
    if missing:
        print(f"{len(missing)} lock(s) carry NO trial record — the sidecar predates "
              f"it, so the budget question cannot be answered from the corpus and "
              f"has to be captured by re-running the row:")
        for name, page in missing[:200]:
            print(f"  {name:22s} {page}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
