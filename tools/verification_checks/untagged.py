#!/usr/bin/env python3
"""Untagged factor-of-safety numbers.

A verification section prints three kinds of number that carry an argument: the
factor of safety a ``<!-- test: -->`` tag locks, the value the source or vendor
published, and the comparison between them.  Anything else that LOOKS like a
factor of safety — a mesh-sweep row, a depth-cutoff row, a with/without variant,
a reading taken off a field at the critical strength — is a companion
measurement that nothing regenerates.  It cannot be defended when a lock moves,
and the section starts contradicting itself the moment one does.

So this check reads every section and reports each factor-of-safety-shaped
number that is neither

* within the tolerance of a tag that section carries, nor
* printed in a column whose header names the source or the vendor,

after the section's own inputs, geometry, dimensions and percentages are taken
out of the running.  A companion measurement that is worth keeping earns its
own tag; the rest is trimmed with the prose that exists to carry it.

Usage: python -m tools.verification_checks.untagged <page.md>
"""
import re
import sys
from decimal import Decimal, InvalidOperation

from .deltas import AUTH_HDR_BASE, XCOL, _hdr_cols, mask_sci

#: The pages set a minus sign as U+2212.
MINUS = "−"

#: A factor of safety, as these pages print one: two to four decimals, between
#: 0.1 and 10.  Tighter than "a decimal number" on purpose — an element size, a
#: coordinate and a unit weight are all written with one decimal or with none,
#: and a stress or a modulus runs to four figures or more.
FS_SHAPED = re.compile(r"(?<![\w.,\-−])(\d\.\d{2,4})(?![\d\w.])")
FS_LO, FS_HI = Decimal("0.10"), Decimal("10")

#: What a number is NOT a factor of safety when it is followed by: a unit, a
#: dimension, a scale.  Matched against the text immediately after the token.
UNIT_AFTER = re.compile(
    r"^\s*(?:%|°|″|'|\"|:|/|×|x\b|m\b|mm\b|m²|m³|ft\b|ft²|ft³|in\b|s\b|"
    r"kPa|MPa|psf|psi|pcf|ksi|kN|kNm|lb|kip|H\b|D\b|g\b|degrees?\b|"
    r"times\b|per\b|elements?\b|nodes?\b)")

#: ... and what it is not one when it is preceded by: a figure or table number,
#: a version, an equation.  Matched against the 24 characters before the token.
LABEL_BEFORE = re.compile(
    r"(?:Fig\.?|Figure|Table|Eq\.?|Equation|Problem|§|v|version|"
    r"Part|item|no\.?)\s*(?:\d+(?:\.\d+)?\s*[-–]\s*)?$", re.I)

#: Spans whose numbers are not prose: inline code, math, link targets, image
#: paths, HTML comments (a test tag is one).  Masked to spaces so every offset
#: stays valid.
SPANS = re.compile(
    r"`[^`]*`|\$[^$]*\$|\]\([^)]*\)|!\[[^\]]*\]|<!--.*?-->|"
    r"https?://\S+|\{[^}]*\}", re.S)

TAG = re.compile(r"<!-- test:(.*?)-->")

#: Tag keys whose value is a locked factor of safety.
LOCK_KEYS = ("expected_fs", "expected", "fs_", "points")

#: The check REPORTS; it does not yet fail a page.  Every finding is a sentence
#: someone has to read — a companion measurement to trim, a number to tag, or a
#: quantity that only looks like a factor of safety and belongs in
#: ``untagged_allow``.  Turning this on before that reading is done would push
#: the pages toward blanket allowances, which is the opposite of what the check
#: is for.  Flip it to True once each page's residue has been adjudicated.
ENFORCING = False


def _mask(text):
    return SPANS.sub(lambda m: " " * len(m.group(0)), text)


def _tag_kv(line):
    m = TAG.search(line)
    if not m:
        return None
    kv = {}
    for part in m.group(1).split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            kv[k.strip()] = v.strip()
    return kv


def _locked(kv):
    """(value, tolerance) pairs a tag locks, as Decimals."""
    tol = Decimal(str(kv.get("tolerance", "0.01")))
    out = []
    for k, v in kv.items():
        if not k.startswith(LOCK_KEYS):
            continue
        for part in str(v).split(";"):
            part = part.split(":")[-1].strip()
            try:
                out.append((Decimal(part), tol))
            except InvalidOperation:
                pass
    return out


def sections(lines):
    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{2,4} ", l)]
    out = []
    for k, (i, l) in enumerate(heads):
        j = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        out.append((i, j, l.strip()))
    if heads and heads[0][0] > 0:
        out.insert(0, (0, heads[0][0], "(preamble)"))
    elif not heads:
        out.append((0, len(lines), "(whole page)"))
    return out


def source_values(lines, sec, auth_hdr):
    """Every number printed in a source/vendor column of this section's tables.

    A table with an XSLOPE column pairs each of its other columns against it, so
    the authority columns are all of them but the XSLOPE one; a table with no
    XSLOPE column is read for the columns whose header names an authority.
    """
    body = lines[sec[0]:sec[1]]

    def pick(hdr):
        if any(XCOL.search(c) for c in hdr):
            return [k for k, c in enumerate(hdr) if not XCOL.search(c)]
        return [k for k, c in enumerate(hdr) if auth_hdr.search(c)]

    vals = set()
    cols = _hdr_cols(body, pick)
    for j, ks in cols.items():
        cells = [c.strip() for c in body[j].strip().strip("|").split("|")]
        for k in ks:
            if k < len(cells):
                for m in FS_SHAPED.finditer(mask_sci(cells[k])):
                    try:
                        vals.add(Decimal(m.group(1)))
                    except InvalidOperation:
                        pass
    return vals


def corpus_locks():
    """Every verification page's locks, keyed by the input file they name.

    A page routinely restates a value another page locks — geostudio's ACADS
    rows are Rocscience's models, and the same workbook is tagged there.  The
    number is tag-guarded either way, so a section that names the file gets the
    benefit of its lock wherever the tag lives.
    """
    import glob
    import os
    out = {}
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(here))
    for page in glob.glob(os.path.join(repo, "docs", "verification", "*.md")):
        for line in open(page):
            kv = _tag_kv(line)
            if not kv or "file" not in kv:
                continue
            out.setdefault(os.path.basename(kv["file"]), []).extend(_locked(kv))
    return out


def run(path, cfg, report=print):
    raw = open(path).read().replace(MINUS, "-")
    lines = raw.split("\n")
    corpus = corpus_locks()
    auth_hdr = re.compile(
        AUTH_HDR_BASE + ("|" + "|".join(cfg.auth_hdr_extra)
                         if cfg.auth_hdr_extra else ""), re.I)
    secs = sections(lines)

    # The dotted corpus-summary table restates each section's headline numbers
    # under a link to it.  Its locks are guarded from the other end, by the tags
    # check's reverse sweep, and reading it here would report every section's
    # source value twice over.
    skip = set()
    insum = False
    for i, line in enumerate(lines):
        if cfg.summary_marker and cfg.summary_marker in line:
            insum = True
        elif insum and line.strip() == "</div>":
            insum = False
        elif insum:
            skip.add(i)

    # A tag placed before the first heading is a page-level bank; its locks
    # belong to whichever section links its input file.
    page_locks = []
    for i, line in enumerate(lines):
        kv = _tag_kv(line)
        if kv is None or any(s[0] <= i < s[1] and s[2] != "(preamble)"
                             for s in secs):
            continue
        page_locks.append((kv.get("file", ""), _locked(kv)))

    allow = list(cfg.untagged_allow)
    fired, flags, scanned = set(), [], 0

    for sec in secs:
        locks = []
        for i in range(sec[0], sec[1]):
            kv = _tag_kv(lines[i])
            if kv:
                locks += _locked(kv)
        body_text = "\n".join(lines[sec[0]:sec[1]])
        for f, ls in page_locks:
            if f and f in body_text:
                locks += ls
        for f, ls in corpus.items():
            if f in body_text:
                locks += ls
        srcs = source_values(lines, sec, auth_hdr)

        for i in range(sec[0], sec[1]):
            line = lines[i]
            if i in skip or line.lstrip().startswith("<!-- test:"):
                continue
            if line.startswith("#"):
                continue          # a heading numbers the section, not a run
            for m in FS_SHAPED.finditer(mask_sci(_mask(line))):
                tok = m.group(1)
                try:
                    v = Decimal(tok)
                except InvalidOperation:
                    continue
                if not (FS_LO <= v <= FS_HI):
                    continue
                if UNIT_AFTER.match(line[m.end():m.end() + 12]):
                    continue
                if LABEL_BEFORE.search(line[max(0, m.start() - 24):m.start()]):
                    continue
                scanned += 1
                if any(abs(v - lv) <= tol for lv, tol in locks):
                    continue
                if v in srcs:
                    continue
                ex = next((e for e in allow
                           if e[0] == tok and e[1] in line), None)
                if ex:
                    fired.add(ex)
                    continue
                flags.append((i + 1, tok, sec[2], line.strip()[:120]))

    dead = [e for e in allow if e not in fired]
    report(f"untagged FS numbers: {len(flags)} flagged of {scanned} "
           f"factor-of-safety-shaped numbers read; dead allowances: "
           f"{len(dead)}{'' if ENFORCING else '  (reported, not enforced)'}")
    for ln, tok, sec, text in flags:
        report(f"  L{ln} {tok}  in {sec[:60]!r} :: {text}")
    for e in dead:
        report(f"  DEAD untagged allowance (never fires): {e[0]!r}, {e[1]!r}")
    return (len(flags) + len(dead)) if ENFORCING else 0


def _cli():
    from .pages import config_for
    path = sys.argv[1]
    return run(path, config_for(path))


if __name__ == "__main__":
    sys.exit(_cli())
