#!/usr/bin/env python3
"""Restatement guard for the tutorial pages under docs/tutorials.

A tutorial page is a walkthrough of runs the reader repeats, so nearly every
factor of safety it prints is a number XSLOPE produced.  Those numbers are what
drift: a solver round moves a locked answer, the verification page and the
sample table are re-measured because a tag guards them, and the tutorial keeps
printing what it printed the day it was written.

This sweep reads every factor-of-safety-shaped number a tutorial ATTRIBUTES TO
ITSELF and asks whether a test tag stands behind it.  The locks in scope for a
page are

* the tags on the page itself, and
* the tags anywhere under ``docs/`` on a model file the page links —
  ``docs/lem/samples.md`` locks the seven methods of the same workbook LEM-3
  walks the reader through, and ``docs/verification`` locks the vendor models
  the LEM pages borrow.

Each number lands in one of three buckets:

``guarded``
    it restates a lock in scope — the lock verbatim, or the lock correctly
    rounded to fewer places.  Agreement is by PRINTED FORM, not by the tag's
    tolerance: a tutorial that prints 1.313 where the lock reads 1.314 is
    restating the lock wrongly however small the gap is, and that one digit is
    the whole failure mode this guard exists for.

``disagreeing``
    the number's own context names a lock — a column headed ``Janbu`` in a
    method table, a row label or a sentence that names the method or the
    strength-reduction run — and it agrees with none of the locks it names.

``unguarded``
    nothing in scope locks the number at all.  A sweep row, a variant run, a
    reading taken off a solved field: reproducible only by hand, and nothing
    tells anyone when it goes stale.

The sweep REPORTS.  It returns a count, never a failure: a tutorial legitimately
re-runs a sample under settings the sample's own tag does not use, and every
finding is a sentence someone has to read before it is a defect.

Usage: python -m tools.verification_checks.tutorials [page.md ...]
"""
import glob
import os
import re
import sys
from collections import namedtuple
from decimal import Decimal, InvalidOperation

from .deltas import AUTH_HDR_BASE, XCOL, mask_sci
from .tags import (METHOD_WORDS, PAREN, SENT_END, SOURCE_EXTRA, _forms,
                   _identity, _tag_kv, _wanted)
from .untagged import (COORD_FIRST, COORD_SECOND, FS_HI, FS_LO, FS_SHAPED,
                       LABEL_BEFORE, MINUS, SEEPAGE_TYPES, TABLE_SEP,
                       UNIT_AFTER, _cells, _mask, quantity_spans, sections)

#: The tag keys whose value is a factor of safety this sweep can match against.
#: ``expected_first`` is the first factor of safety of an FS-vs-time march and
#: ``expected`` its whole series, one value per step.
LOCK_KEYS = ("expected_fs*", "fs_*", "expected", "expected_first")

#: Decimal places a tutorial may restate a lock to.  A lock of 1.3711 is
#: satisfied by a printed 1.371 or 1.37, never by 1.372.
ROUND_DP = 2

#: A page whose model files carry only these lock types publishes heads,
#: pressures and flow rates; nothing shaped like a factor of safety in it is
#: one.  Inherited from the untagged sweep, and applied per section for the same
#: reason: one stability lock anywhere in a section puts it back in the running.
SEEP_ONLY = SEEPAGE_TYPES

#: Recognises the model files a page links, so the page inherits their locks.
MODEL_FILE = re.compile(r"[\w./-]+\.xlsx")

#: A table column header, or a row label, that names an authority rather than
#: this program: the cells under it are the source's numbers, not a run's.
AUTH = re.compile(AUTH_HDR_BASE, re.I)

#: A source-side attributor in prose.  What follows belongs to the source until
#: the sentence ends.
SOURCE = re.compile(SOURCE_EXTRA + "|" + AUTH_HDR_BASE, re.I)

Lock = namedtuple("Lock", "value tol slots origin benchmark key type file series")

Finding = namedtuple("Finding", "line token verdict context lock")


def _locks_of(kv, origin, base):
    """The factor-of-safety locks one test tag carries, with their identity."""
    try:
        tol = Decimal(str(kv.get("tolerance", "0.01")))
    except InvalidOperation:
        tol = Decimal("0.01")
    out = []
    for k, v in kv.items():
        if not _wanted(k, LOCK_KEYS):
            continue
        slots = _identity(k, kv)
        # A series lock names one instant of a march, not the page's answer for
        # that method: ``expected=`` lists a factor of safety per step, and
        # ``expected_first`` locks the first step of one.
        series = ";" in str(v) or k == "expected_first"
        for part in str(v).split(";"):
            part = part.split(":")[-1].strip()
            try:
                Decimal(part)
            except InvalidOperation:
                continue
            out.append(Lock(part, tol, slots, origin, kv.get("benchmark", "?"),
                            k, kv.get("type", ""), base, series))
    return out


def doc_locks(repo=None):
    """Every test tag under ``docs/``, keyed by the model file it names.

    The map is what lets a tutorial inherit a guard it does not carry itself:
    LEM-3 walks ``xslope_simple_mult_layers.xlsx``, whose seven method locks
    live on ``docs/lem/samples.md``, and LEM-9 walks a Rocscience model locked
    on a verification page.  Keyed by base name, because a tag writes the path
    relative to its own page.
    """
    if repo is None:
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(os.path.dirname(here))
    by_file = {}
    pattern = os.path.join(repo, "docs", "**", "*.md")
    for page in sorted(glob.glob(pattern, recursive=True)):
        rel = os.path.relpath(page, repo)
        for i, line in enumerate(open(page).read().split("\n")):
            kv = _tag_kv(line)
            if not kv or "file" not in kv:
                continue
            base = os.path.basename(kv["file"])
            by_file.setdefault(base, []).extend(
                _locks_of(kv, f"{rel}:{i + 1}", base))
    return by_file


def _agrees(token, lock):
    """True when a printed token restates `lock` as the lock is written.

    The token must be the lock's own digits, or the lock correctly rounded to
    fewer places (never more, and never past ``ROUND_DP``).  The tag's
    tolerance is deliberately NOT consulted: a tolerance says how far a rerun
    may land from the locked answer, not how far the prose may.
    """
    return token in _forms(lock.value, ROUND_DP)


def _names(slots, text):
    """True when `text` names every slot of a lock's identity.

    Whole words only.  A method name is short enough to hide inside an ordinary
    one — "Lowe" inside "lower row only", "the shallower circle" — and a row
    label that happens to contain one is not a row about that method.
    """
    return bool(slots) and all(
        any(re.search(r"\b" + re.escape(a) + r"\b", text) for a in slot)
        for slot in slots)


def _method_slots(text):
    """Identity slots a header cell, a row label or a sentence declares."""
    low = text.lower()
    slots = []
    for words in METHOD_WORDS.values():
        if any(w in low for w in words):
            slots.append(words)
    if re.search(r"\bssrm?\b|\bsrf\b|strength reduction", low):
        slots.append(("ssrm", "srf", "strength reduction"))
    return slots


#: A run SETTING, not a result: the ends of a strength-reduction bracket, a
#: convergence tolerance, a relaxation factor, a residual, a mesh size.  These
#: are typed by the reader, so a tutorial printing one is quoting its own input.
#: Read against the 80 characters before the token, because the setting is
#: normally named at the head of the sentence that gives its value.
INPUT_BEFORE = re.compile(
    r"\b(?:F ?min|F ?max|f_min|f_max|(?:lower|upper) bound|bracket|"
    r"relax(?:ation)?|residual|closure|"
    r"tolerance|grid|target (?:element )?size|element size|convergence)\b"
    r"(?:[^.]|\.(?=\d)){0,60}$", re.I)

#: A console line reporting an INTERMEDIATE search iteration.  The factor of
#: safety on it is a step on the way to the answer, not the answer.
ITERATION_LINE = re.compile(r"iteration\s*\d|^\s*Iteration\s", re.I)

#: What labels a number in a console transcript as a factor of safety.  A code
#: fence is a verbatim log — grid spacings, coordinates, residuals and factors
#: of safety in one stream — so only the numbers the log itself names are read.
FS_LABEL = re.compile(r"\b(?:FS|FOS|SRF|factor of safety)\s*[=:]\s*$", re.I)


def _fenced(lines):
    """Line indices inside a ``` code fence."""
    out, on = set(), False
    for i, line in enumerate(lines):
        if line.lstrip().startswith("```"):
            on = not on
            out.add(i)
            continue
        if on:
            out.add(i)
    return out


def _candidate(line, token, start, end):
    """The standard reasons a factor-of-safety-shaped token is not one."""
    try:
        value = Decimal(token)
    except InvalidOperation:
        return False
    if not (FS_LO <= value <= FS_HI):
        return False
    after = line[end:end + 12]
    before = line[max(0, start - 80):start]
    if INPUT_BEFORE.search(before):
        return False
    before = before[-24:]
    if UNIT_AFTER.match(after) or LABEL_BEFORE.search(before):
        return False
    if COORD_SECOND.search(before) or (
            before.rstrip().endswith("(") and COORD_FIRST.match(after)):
        return False
    return True


def _tables(lines, sec):
    """Per table-body line: (row label, {column -> header}, auth columns, width).

    A tutorial's method table puts the methods in the HEADER and one row of
    answers under it — ``| OMS | Bishop | Janbu | ... |`` — where a verification
    table puts the method in the row label and the answer in an XSLOPE column.
    Both are read, and a cell's identity is its column header and its row label
    together, so either orientation binds a number to the lock it restates.
    """
    out, i = {}, sec[0]
    while i < sec[1] - 1:
        if "|" not in lines[i] or not TABLE_SEP.match(lines[i + 1]):
            i += 1
            continue
        hdr = [t for t, _, _ in _cells(lines[i])]
        auth = {k for k, t in enumerate(hdr)
                if AUTH.search(t) and not XCOL.search(t)
                and not _method_slots(t)}
        j = i + 2
        body = []
        while j < sec[1] and lines[j].lstrip().startswith("|"):
            body.append(j)
            j += 1
        nxt = j
        # A column header binds its cells to a lock only where the table has ONE
        # body row.  With several, the rows are variants of each other — a
        # pore-pressure option per row, a surcharge per row — and a lock, if one
        # exists, belongs to a single row; the header would otherwise bind every
        # variant to the same locked answer.  Row labels still bind, so a table
        # written the other way round (method per row) reads normally.
        single = len(body) == 1
        # A table whose every column after the first is headed by a NUMBER is a
        # sweep — a load, a dip angle, a GSI per column — so its cells are a
        # curve of results and the method its row label names holds for the
        # whole curve, not for any one column.  Nothing in it restates a lock.
        sweep = len(hdr) > 1 and not any(re.search(r"[A-Za-z]", h)
                                         for h in hdr[1:])
        for j in body:
            cells = _cells(lines[j])
            # The first cell is a LABEL only where it is words.  A method table
            # written with the methods in the header has a number there, and
            # reading it as the row's name would bind every cell of the row to
            # whatever the first answer happens to be.
            label = cells[0][0] if cells else ""
            if not re.search(r"[A-Za-z]", label):
                label = ""
            headers = {k: hdr[k] for k in range(min(len(hdr), len(cells)))
                       if re.search(r"[A-Za-z]", hdr[k])} if single else {}
            if sweep:
                label = ""
            out[j] = (label, headers, auth, len(hdr))
        i = max(nxt, i + 1)
    return out


def _classify(token, context, locks):
    """(verdict, the lock the finding names) for one attributed number.

    A number is GUARDED by any lock in scope it restates, named or not: the
    page prints the locked answer, which is all the guard asks.  It is only
    DISAGREEING where its own context names a lock — a method, or the
    strength-reduction run — and none of the locks it names is what it prints.

    A lock that is one element of a series (an FS-vs-time march, a stage set)
    never narrows.  The series says what the answer is at each of its own
    instants and claims nothing about the other numbers the page prints under
    the same method, so treating it as an identity would make every other
    reading on the page disagree with it.
    """
    named = [l for l in locks if not l.series and _names(l.slots, context)]
    hit = next((l for l in (named or locks) if _agrees(token, l)), None)
    if hit is None and named:
        hit = next((l for l in locks if _agrees(token, l)), None)
    if hit is not None:
        return "guarded", hit
    if not named:
        return "unguarded", None
    near = min(named, key=lambda l: abs(Decimal(token) - Decimal(l.value)))
    return "disagreeing", near


def scan(path, by_file=None, repo=None):
    """Every attributed factor of safety on one tutorial, with its verdict."""
    if by_file is None:
        by_file = doc_locks(repo)
    raw = open(path).read().replace(MINUS, "-")
    lines = raw.split("\n")
    files = {os.path.basename(m) for m in MODEL_FILE.findall(raw)}
    locks = [l for f in sorted(files) for l in by_file.get(f, ())]
    types = {l.type for l in locks}
    fence = _fenced(lines)

    findings = []
    if not locks:
        return findings, files
    for sec in sections(lines):
        body = "\n".join(lines[sec[0]:sec[1]])
        # What quantity the section's numbers ARE is said by the tags on the
        # files it names; a section that names no file inherits the page's.
        named_here = {f for f in files if f in body}
        own = {kv.get("type", "") for i in range(sec[0], sec[1])
               if (kv := _tag_kv(lines[i]))}
        stypes = own or {l.type for l in locks
                         if l.file in named_here} or types
        if stypes and stypes <= SEEP_ONLY:
            continue
        qspans = quantity_spans(lines, sec)
        tables = _tables(lines, sec)
        blocked, sentence = False, ""
        for i in range(sec[0], sec[1]):
            line = lines[i]
            stripped = line.strip()
            if not stripped or stripped.startswith("<!--") \
                    or line.startswith("#"):
                blocked, sentence = False, ""
                continue
            if i in fence:
                # A console transcript: read the numbers the log names as a
                # factor of safety, and nothing else on the line.
                blocked, sentence = False, ""
                if ITERATION_LINE.search(line):
                    continue
                for m in FS_SHAPED.finditer(mask_sci(_mask(line))):
                    tok = m.group(1)
                    if not FS_LABEL.search(line[max(0, m.start(1) - 24):
                                                m.start(1)]):
                        continue
                    if not _candidate(line, tok, m.start(1), m.end(1)):
                        continue
                    verdict, lock = _classify(tok, line.lower(), locks)
                    findings.append(Finding(i + 1, tok, verdict,
                                            "console log", lock))
                continue
            masked = mask_sci(_mask(line))
            row = tables.get(i)
            if row is not None:
                label, headers, auth, ncol = row
                if TABLE_SEP.match(line):
                    continue
                blocked, sentence = False, ""
                for k, (text, a, b) in enumerate(_cells(line)):
                    # The first cell of a multi-column table is the row's name.
                    if k in auth or (k == 0 and ncol > 1):
                        continue
                    ctx = (label + " " + headers.get(k, "")).lower()
                    cell = PAREN.sub(lambda p: " " * len(p.group(0)),
                                     mask_sci(_mask(text)))
                    for m in FS_SHAPED.finditer(cell):
                        tok = m.group(1)
                        s, e = a + m.start(1), a + m.end(1)
                        if not _candidate(line, tok, s, e):
                            continue
                        if any(x <= s and e <= y
                               for x, y in qspans.get(i, ())):
                            continue
                        verdict, lock = _classify(tok, ctx, locks)
                        findings.append(Finding(i + 1, tok, verdict,
                                                ctx.strip() or "(table cell)",
                                                lock))
                continue
            if stripped.startswith("|"):
                continue
            pos = 0
            scanner = re.compile(
                "|".join(f"(?:{p})" for p in (SOURCE.pattern,
                                              SENT_END.pattern,
                                              FS_SHAPED.pattern)), re.I)
            for m in scanner.finditer(masked):
                sentence += masked[pos:m.start()]
                pos = m.end()
                tok = m.group(0)
                if SENT_END.fullmatch(tok):
                    blocked, sentence = False, ""
                    continue
                if m.group(1) is None:
                    blocked, sentence = True, sentence + tok
                    continue
                sentence += tok
                if blocked or not _candidate(line, tok, m.start(), m.end()):
                    continue
                if any(x <= m.start() and m.end() <= y
                       for x, y in qspans.get(i, ())):
                    continue
                verdict, lock = _classify(tok, sentence.lower(), locks)
                findings.append(Finding(i + 1, tok, verdict,
                                        "prose", lock))
            sentence += masked[pos:]
    return findings, files


def tutorial_pages(repo=None):
    if repo is None:
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(os.path.dirname(here))
    return [p for p in sorted(glob.glob(os.path.join(repo, "docs", "tutorials",
                                                     "*.md")))
            if os.path.basename(p) != "index.md"]


def run(paths=None, report=print, verbose=True, repo=None):
    """Sweep the tutorials and report the per-page tally.

    Returns the number of disagreeing restatements, which the suite prints and
    does not fail on: the guard reports until every page's residue has been
    read.
    """
    by_file = doc_locks(repo)
    pages = paths or tutorial_pages(repo)
    total = {"guarded": 0, "disagreeing": 0, "unguarded": 0}
    rows = []
    for path in pages:
        findings, _files = scan(path, by_file)
        tally = {k: sum(1 for f in findings if f.verdict == k)
                 for k in total}
        for k in total:
            total[k] += tally[k]
        rows.append((os.path.basename(path), tally, findings))
    width = max((len(r[0]) for r in rows), default=10)
    report(f"tutorial restatements: {sum(total.values())} attributed numbers "
           f"read across {len(rows)} pages — guarded {total['guarded']}, "
           f"disagreeing {total['disagreeing']}, unguarded "
           f"{total['unguarded']}  (reported, not enforced)")
    report(f"  {'page'.ljust(width)}  guarded  disagreeing  unguarded")
    for name, tally, _f in rows:
        report(f"  {name.ljust(width)}  {tally['guarded']:7d}  "
               f"{tally['disagreeing']:11d}  {tally['unguarded']:9d}")
    if verbose:
        for name, _t, findings in rows:
            for f in findings:
                if f.verdict == "guarded":
                    continue
                if f.verdict == "disagreeing":
                    report(f"  {name} L{f.line} prints {f.token} "
                           f"({f.context[:48]}) but the lock is "
                           f"{f.lock.key}={f.lock.value} [{f.lock.benchmark}], "
                           f"tagged at {f.lock.origin}")
                else:
                    report(f"  {name} L{f.line} {f.token} unguarded "
                           f"({f.context[:48]})")
    return total["disagreeing"]


def _cli():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    quiet = "--quiet" in sys.argv[1:]
    run(args or None, verbose=not quiet)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
