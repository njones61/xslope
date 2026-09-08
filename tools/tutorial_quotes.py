#!/usr/bin/env python3
"""Which tutorial lines quote a lock that is about to move.

A re-lock round re-records a tag and re-measures the verification section that
carries it.  The tutorials are the pages that quietly go stale: they walk the
reader through the same model, print the same answer in prose, in a results
table and in a console transcript, and nothing in the round's own scope looks at
them.  The 2026-09-07 strength-reduction re-lock moved 56 locks and needed a
separate sweep of every other page afterwards to find the four that restated
one.

This is that sweep, run BEFORE the round instead of after it.  Name what is
moving — a benchmark id, a lock value, or the tags a diff already changed — and
it lists every tutorial line that quotes it, so the pages are fixed in the same
round rather than found later.

    # every tutorial line that quotes this benchmark's locks, or names it
    python tools/tutorial_quotes.py --benchmark FEM-1-ssrm

    # a value on its way out, before the tag is edited
    python tools/tutorial_quotes.py --value 1.3633 --value 1.587

    # every lock a diff has already moved, old value by old value
    python tools/tutorial_quotes.py --since HEAD~1
    python tools/tutorial_quotes.py --since b44de48e --root /tmp/old/tutorials

A quoted value is matched at the precision the tutorials write it: the lock's
own digits, and the lock correctly rounded down to ``--dp`` places (three by
default, which is what a page means when it prints FEM-1's 1.3633 as 1.363).
Matching is bounded on both sides, so 1.363 is not found inside 11.3634.

Exit status is 0 whether or not anything is found: this reports, it does not
gate.
"""
import argparse
import os
import re
import subprocess
import sys
from decimal import Decimal, InvalidOperation

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.verification_checks.tags import _forms, _tag_kv, _wanted  # noqa: E402

#: The tag keys that hold a locked factor of safety.  Same set the tutorial
#: restatement sweep scores against, so the two tools agree on what a lock is.
LOCK_KEYS = ("expected_fs*", "fs_*", "expected", "expected_first")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: A test tag on a line of a diff, with the +/- stripped.
DIFF_TAG = re.compile(r"^([+-])(?!\+\+|--)(.*<!-- test:.*)$")


def _lock_values(kv):
    """{tag key: [value, ...]} for the locks one tag carries."""
    out = {}
    for k, v in kv.items():
        if not _wanted(k, LOCK_KEYS):
            continue
        vals = []
        for part in str(v).split(";"):
            part = part.split(":")[-1].strip()
            try:
                Decimal(part)
            except InvalidOperation:
                continue
            vals.append(part)
        if vals:
            out[k] = vals
    return out


def doc_tags(repo=REPO):
    """Every ``<!-- test: -->`` tag under docs/, as (page, line, kv)."""
    import glob
    out = []
    for page in sorted(glob.glob(os.path.join(repo, "docs", "**", "*.md"),
                                 recursive=True)):
        rel = os.path.relpath(page, repo)
        for i, line in enumerate(open(page).read().split("\n")):
            kv = _tag_kv(line)
            if kv:
                out.append((rel, i + 1, kv))
    return out


def by_benchmark(repo=REPO):
    """{benchmark id: [(page, line, kv), ...]}."""
    out = {}
    for page, line, kv in doc_tags(repo):
        if "benchmark" in kv:
            out.setdefault(kv["benchmark"], []).append((page, line, kv))
    return out


def moved_locks(since, repo=REPO):
    """Locks a diff has already moved: [(old, new, key, where, model)].

    Read from ``git diff <since>`` over docs/, pairing each removed tag with the
    added tag that carries the same input file and benchmark id.  A key whose
    value differs between the two is a lock the round moved, and its OLD value
    is what the tutorials may still be printing.
    """
    diff = subprocess.run(
        ["git", "-C", repo, "diff", "-U0", since, "--", "docs"],
        capture_output=True, text=True, check=True).stdout
    removed, added = {}, {}
    for line in diff.split("\n"):
        m = DIFF_TAG.match(line)
        if not m:
            continue
        kv = _tag_kv(m.group(2))
        if not kv:
            continue
        key = (os.path.basename(kv.get("file", "")), kv.get("benchmark", ""),
               kv.get("type", ""), kv.get("element_type", ""),
               kv.get("method", ""), kv.get("f_min", ""))
        (removed if m.group(1) == "-" else added).setdefault(key, []).append(kv)
    out = []
    for key, olds in removed.items():
        news = added.get(key, [])
        for n, old in enumerate(olds):
            new = news[n] if n < len(news) else {}
            ovals, nvals = _lock_values(old), _lock_values(new)
            for k, vals in ovals.items():
                for j, v in enumerate(vals):
                    nv = nvals.get(k, [])
                    if j < len(nv) and nv[j] == v:
                        continue
                    where = old.get("benchmark") or key[0] or "?"
                    out.append((v, nv[j] if j < len(nv) else None, k, where,
                                key[0]))
    return sorted(set(out))


def _needles(value, dp):
    """The strings a page may print `value` as, longest first."""
    forms = [f for f in _forms(value, dp)
             if "." not in f or len(f.split(".")[1]) >= dp]
    return sorted(set(forms), key=len, reverse=True)


def pages_under(root):
    out = []
    for base, _dirs, names in os.walk(root):
        for name in sorted(names):
            if name.endswith(".md") and name != "index.md":
                out.append(os.path.join(base, name))
    return sorted(out)


def quotes(root, needles, repo=REPO, model=None):
    """[(page, line, text, needle, same model?)] for lines printing a needle.

    ``model`` is the base name of the workbook the lock is recorded on.  A page
    that links it is walking the very model whose answer is moving, so its hit
    is almost certainly a restatement; a page that does not may be printing the
    same three digits about something else entirely.  The distinction is
    marked, never used to drop a line — the Cai & Ugai table on FEM-3 restates a
    verification lock on a model FEM-3 never links.
    """
    pats = [(n, re.compile(r"(?<![\w.,\-−])" + re.escape(n) + r"(?![\d\w.])"))
            for n in needles]
    hits = []
    for path in pages_under(root):
        rel = os.path.relpath(path, repo) if path.startswith(repo) else path
        body = open(path).read()
        same = bool(model) and model in body
        for i, line in enumerate(body.split("\n")):
            if line.lstrip().startswith("<!-- test:"):
                continue          # the tag is what moves; it is not a quote
            for needle, pat in pats:
                if pat.search(line.replace("−", "-")):
                    hits.append((rel, i + 1, line.strip(), needle, same))
                    break
    return sorted(hits, key=lambda h: (not h[4], h[0], h[1]))


def names(root, ids, repo=REPO):
    """[(page, line no, line text, id)] for every line naming a benchmark id."""
    hits = []
    for path in pages_under(root):
        rel = os.path.relpath(path, repo) if path.startswith(repo) else path
        for i, line in enumerate(open(path).read().split("\n")):
            if line.lstrip().startswith("<!-- test:"):
                continue
            for ident in ids:
                if ident in line:
                    hits.append((rel, i + 1, line.strip(), ident))
                    break
    return hits


def run(targets, ids, root, dp=3, repo=REPO, report=print):
    """Report every tutorial line quoting one of `targets` or naming an id.

    ``targets`` is [(old, new or None, tag key, where it is locked, model)].
    Returns the number of lines found.
    """
    total = 0
    for old, new, key, where, model in targets:
        hits = quotes(root, _needles(old, dp), repo, model)
        moving = f"{old} -> {new}" if new else old
        n_same = sum(1 for h in hits if h[4])
        report(f"\n{key}={moving}  [{where}]  — {len(hits)} tutorial line"
               f"{'' if len(hits) == 1 else 's'}"
               f"{f', {n_same} on a page that links {model}' if n_same else ''}")
        for page, line, text, needle, same in hits:
            report(f"  {'*' if same else ' '} {page}:{line}  ({needle})  "
                   f"{text[:110]}")
        total += len(hits)
    if ids:
        hits = names(root, ids, repo)
        report(f"\nbenchmark id named in prose — {len(hits)} line"
               f"{'' if len(hits) == 1 else 's'}")
        for page, line, text, ident in hits:
            report(f"  {page}:{line}  ({ident})  {text[:110]}")
        total += len(hits)
    report(f"\n{total} tutorial line{'' if total == 1 else 's'} quote a moving "
           f"lock across {len(targets)} value{'' if len(targets) == 1 else 's'}")
    return total


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="List tutorial lines that quote a lock about to move.")
    ap.add_argument("--benchmark", action="append", default=[],
                    help="benchmark id whose locks are moving (repeatable)")
    ap.add_argument("--value", action="append", default=[],
                    help="a lock value on its way out (repeatable)")
    ap.add_argument("--since", default=None,
                    help="git ref: take every lock the diff since it moved")
    ap.add_argument("--root", default=None,
                    help="directory of pages to search "
                         "(default docs/tutorials)")
    ap.add_argument("--dp", type=int, default=3,
                    help="fewest decimal places a page may print a lock to "
                         "(default 3)")
    ap.add_argument("--repo", default=REPO, help="repository root")
    args = ap.parse_args(argv)

    root = args.root or os.path.join(args.repo, "docs", "tutorials")
    targets, ids = [], list(args.benchmark)

    for v in args.value:
        targets.append((v, None, "value", "given", None))

    if args.benchmark:
        index = by_benchmark(args.repo)
        for ident in args.benchmark:
            tagged = index.get(ident)
            if not tagged:
                print(f"no tag carries benchmark={ident}", file=sys.stderr)
                continue
            for page, line, kv in tagged:
                model = os.path.basename(kv.get("file", ""))
                for key, vals in _lock_values(kv).items():
                    for v in vals:
                        targets.append((v, None, key, f"{page}:{line}", model))

    if args.since:
        targets += moved_locks(args.since, args.repo)

    if not targets and not ids:
        ap.error("name something that is moving: --benchmark, --value "
                 "or --since")

    seen, uniq = set(), []
    for t in targets:
        if t[:2] not in seen:
            seen.add(t[:2])
            uniq.append(t)
    run(uniq, ids, root, args.dp, args.repo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
