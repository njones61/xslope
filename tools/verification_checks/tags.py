#!/usr/bin/env python3
"""Tag-vs-text audit, BIDIRECTIONAL.

Forward : every locked value in a test tag (`expected_fs`, `expected_kc`,
          `expected_flowrate`, `fs_*`, ...) must appear verbatim in the section
          that carries the tag.  A tag value that is a semicolon LIST —
          ``expected=1.686;0.941;...`` (a factor of safety per march step),
          ``points=20:5:7.166;...`` (a solved head per station) — locks one
          value per element, and each element is checked the same way, so a
          wrong digit in one printed cell of an eleven-row table fails.
Restate: the forward pass is existential — it asks whether the lock is printed,
          not whether every printing of it is right.  A section that states its
          locked factor of safety twice (a results-table cell and a note under
          it) could therefore drift in one place with the other left standing.
          So every number the section ATTRIBUTES TO ITSELF — a results-table
          cell in an XSLOPE column of a row whose label names the locked method,
          or a number introduced in prose by "XSLOPE's" / "this model's" — must
          agree with one of the locks that attribution reaches.
Reverse : every value the page PRESENTS AS LOCKED must have a tag behind it.
          A value is presented as locked when it appears in the Results column
          of a summary row whose Match cell carries a colour dot (a dot scores
          "the match quality of what is locked").

The tag is truth; a disagreement is a text defect — the text conforms, never
the tag.

Usage: python -m tools.verification_checks.tags <page.md>
"""
import re
import sys
from decimal import Decimal, ROUND_HALF_UP

TAG = re.compile(r"<!-- test:(.*?)-->")
DOT = re.compile(r"🟢|🟡|🔴|🟣")

#: The pages set a minus sign as U+2212; tag values are ASCII.  Normalised in
#: the body before matching so a locked −0.0932 is found where it is printed.
MINUS = "−"

#: A number printed in the pages' scientific notation, `2.500×10⁻⁵`.
SCI_PRINTED = re.compile(r"(\d+(?:\.\d+)?)\s?×\s?10\s?(⁻|⁺)?([⁰¹²³⁴⁵⁶⁷⁸⁹]+)")
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")


def _sci_match(value, body):
    """True when `body` prints `value` in scientific notation.

    A discharge tag reads ``expected_flowrate=2.3069e-05`` while the page
    prints ``2.307×10⁻⁵``: the same number, at the precision the page reports
    it to.  The mantissa is therefore compared numerically, rescaled to the
    printed exponent and rounded to the printed mantissa's own precision, so a
    different exponent or a different value still fails.
    """
    try:
        want = Decimal(value)
    except Exception:
        return False
    for mant, sign, exp in SCI_PRINTED.findall(body):
        e = int(exp.translate(SUP)) * (-1 if sign == "⁻" else 1)
        dp = len(mant.split(".")[1]) if "." in mant else 0
        try:
            got = want.scaleb(-e).quantize(Decimal(1).scaleb(-dp),
                                           rounding=ROUND_HALF_UP)
        except Exception:
            continue
        if got == Decimal(mant):
            return True
    return False


def _forms(v, dp):
    """The strings that satisfy a tag value.

    The value itself, plus — where the page config allows it — the same value
    correctly rounded, at every precision from one place shy of the tag's own
    down to ``dp``.  Rounding only ever coarsens, so a genuinely different
    value still fails, and a value in exponential notation is left alone
    (quantizing ``2.500e-05`` to two places would offer a meaningless ``0.00``
    as a match; exponential values are matched by `_sci_match` instead).
    """
    out = [v]
    # A whole number large enough to be digit-grouped: the pages print element
    # and node counts as "3 166" (and occasionally "3,166"), which is the same
    # integer the mesh-size tag locks.
    if v.isdigit() and len(v) > 3:
        out.append(f"{int(v):,}".replace(",", " "))
        out.append(f"{int(v):,}")
    if dp is None or "." not in v or "e" in v.lower():
        return out
    nd = len(v.split(".")[1])
    for d in range(nd - 1, dp - 1, -1):
        try:
            out.append(str(Decimal(v).quantize(Decimal(1).scaleb(-d),
                                               rounding=ROUND_HALF_UP)))
        except Exception:
            pass
    return out


def _printed(v, body, dp):
    """True when `body` prints the tag value `v` at some accepted precision."""
    if any(re.search(r"(?<![\w.])" + re.escape(f) + r"(?![\w.])", body)
           for f in _forms(v, dp)):
        return True
    return bool("e" in v or "E" in v) and _sci_match(v, body)


def _elements(v):
    """The locked values of a semicolon list, one per element.

    ``expected=1.686;0.941`` locks the element itself.  ``points=20:5:7.166``
    locks the head at a station: the coordinates are inputs, like `time=` or
    `target_size=`, which no check reads back off the page, so the locked value
    is the last colon field.
    """
    return [p.split(":")[-1].strip() for p in v.split(";")]


def _element_forms(part):
    """The values that satisfy one list element.

    The locked value, and — for a ``points`` element, whose coordinates travel
    with it — the same lock read as a PRESSURE head, ψ = h − y at that station.
    Sections routinely tabulate a solved field in whichever of the two the
    comparison is made in (SEEPW-T02 prints ψ against SEEP/W's ψ, GW19 prints ψ
    along the top boundary against the vendor's figure), and ψ is computed from
    the element's own y, so accepting it is exactly as strong as matching h.
    """
    f = [p.strip() for p in part.split(":")]
    out = [f[-1]]
    if len(f) >= 3:
        try:
            out.append(str(Decimal(f[-1]) - Decimal(f[-2])))
        except Exception:
            pass
    return out


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


def _wanted(key, patterns):
    for p in patterns:
        if p.endswith("*"):
            if key.startswith(p[:-1]):
                return True
        elif key == p:
            return True
    return False


def forward(path, cfg, report=print):
    """Every tagged value must be printed in the section carrying the tag."""
    lines = [l.replace(MINUS, "-") for l in open(path).read().split("\n")]
    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{2,3} ", l)]
    bounds = []
    for k, (i, l) in enumerate(heads):
        j = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        bounds.append((i, j, l))

    def body_of(b):
        return "\n".join(l for l in lines[b[0]:b[1]]
                         if not l.lstrip().startswith("<!-- test:"))

    def owning_section(kv):
        """The section that publishes a page-level tag's results: the one
        section that links the tag's input file.  A list lock is counted, not
        merely searched for, so its scope has to be the section that presents
        it — over a whole-page body an unrelated number could stand in for a
        missing element.  ``None`` when no single section claims the file."""
        f = kv.get("file")
        if not f:
            return None
        hits = [b for b in bounds if f in "\n".join(lines[b[0]:b[1]])]
        return hits[0] if len(hits) == 1 else None

    problems, rows, fired = [], [], set()
    lists, elements, decl_fired = 0, 0, set()

    for idx, line in enumerate(lines):
        kv = _tag_kv(line)
        if kv is None:
            continue
        # A tag placed before the first heading is a page-level tag bank, not a
        # section tag: its scope is the whole page.
        sec = next((b for b in bounds if b[0] <= idx < b[1]), None)
        page_level = sec is None
        if page_level:
            sec = (0, len(lines), "(page-level tag bank)")
        body = body_of(sec)
        vals = {k: v for k, v in kv.items() if _wanted(k, cfg.tag_value_keys)}
        for k, v in vals.items():
            if ";" in v:
                lists += 1
                parts = v.split(";")
                els = _elements(v)
                elements += len(els)
                lsec = sec
                if page_level:
                    lsec = owning_section(kv) or sec
                lbody = body_of(lsec)
                ok = [any(_printed(f, lbody, cfg.tag_round_dp)
                          for f in _element_forms(p)) for p in parts]
                hit = [e for e, o in zip(els, ok) if o]
                rows += [dict(line=idx + 1, bench=kv.get("benchmark", "?"),
                              key=k, value=e, section=lsec[2].strip(),
                              in_section=o) for e, o in zip(els, ok)]
                # A page may publish only part of a probe set — a table of the
                # stations that carry the argument, a regression set it prints
                # nothing of.  The config declares how many elements the page
                # prints, and the count is exact, so a mistyped printed value
                # (one fewer found) fails just as an undeclared one does.
                dec = next(((s, n) for s, n in cfg.tag_list_published
                            if s in line), None)
                want = len(els) if dec is None else dec[1]
                if dec is not None:
                    decl_fired.add(dec)
                if len(hit) == want:
                    continue
                if dec is None:
                    problems += [(idx + 1, kv.get("benchmark", "?"),
                                  f"{k}[]", e, lsec[2].strip())
                                 for e in els if e not in hit]
                else:
                    problems.append(
                        (idx + 1, kv.get("benchmark", "?"), f"{k}[]",
                         f"{len(hit)} of {len(els)} elements printed, "
                         f"declared {want}", lsec[2].strip()))
                continue
            found = _printed(v, body, cfg.tag_round_dp)
            rows.append(dict(line=idx + 1, bench=kv.get("benchmark", "?"),
                             key=k, value=v, section=sec[2].strip(),
                             in_section=bool(found)))
            if found:
                continue
            ex = next(((val, key) for val, key in cfg.tag_exempt
                       if val == v and key in line), None)
            if ex:
                fired.add(ex)
                continue
            problems.append((idx + 1, kv.get("benchmark", "?"), k, v,
                             sec[2].strip()))
    dead = [e for e in cfg.tag_exempt if e not in fired]
    dead_decl = [d for d in cfg.tag_list_published if d not in decl_fired]
    report(f"tags scanned: {len({r['line'] for r in rows})}  "
           f"values checked: {len(rows)} ({lists} list locks, {elements} "
           f"elements)  missing from their section: {len(problems)}  "
           f"dead tag exemptions: {len(dead) + len(dead_decl)}")
    for p in problems:
        report(f"  L{p[0]} [{p[1]}] {p[2]}={p[3]}  -- not in {p[4]!r}")
    for val, key in dead:
        report(f"  DEAD tag exemption (never fires): {val!r}, {key!r}")
    for sub, n in dead_decl:
        report(f"  DEAD list declaration (matches no tag): {sub!r}, {n}")
    return len(problems) + len(dead) + len(dead_decl)


#: An XSLOPE-side attributor: the sentence says the number that follows is this
#: program's own result.  A number with no attributor before it in its sentence
#: is not claimed by the page as a restatement of anything, and belongs to the
#: untagged sweep rather than to this one.
XSIDE = re.compile(
    r"\bXSLOPE(?:'s|’s)?\b|\bthis model(?:'s|’s)?\b|"
    r"\bthe model(?:'s|’s)?\b|\bthis run(?:'s|’s)?\b", re.I)

#: ... and a source-side attributor, which withdraws it: what follows belongs to
#: the source until an XSLOPE-side marker takes the sentence back.  The page's
#: own authority vocabulary (the column headers its tables use for the source)
#: plus the generic ways prose names one.
SOURCE_EXTRA = (
    r"\bthe (?:paper|source|authors?|chart|charts|vendor|program|reference)"
    r"(?:'s|’s)?\b|\btheir\b|\bpublished\b|\breported\b|\bquotes?\b|"
    r"\bFE (?:FOS|FS)\b|\bSLOPE/W\b|\bSEEP/W\b|\bGMS\b|\bSEEP2D\b")

#: Ends a sentence, and with it the attribution the sentence carried.  Matched
#: only where a period is followed by space, so a decimal point never resets;
#: an abbreviation ("p. 394", "Fig. 5") resets early, which only ever makes the
#: sweep read fewer numbers.
SENT_END = re.compile(r"[.;:!?](?=\s|$)")

#: The words a tag's own fields put into the label of the row that publishes it.
#: Each slot is a tuple of alternatives; a row label (or a sentence) names a
#: lock when every slot the lock declares is named in it.
METHOD_WORDS = {
    "spencer": ("spencer",),
    "bishop": ("bishop",),
    "janbu": ("janbu",),
    "oms": ("oms", "ordinary method", "fellenius"),
    "corps": ("corps",),
    "lowe": ("lowe", "karafiath"),
    "mprice": ("morgenstern", "m-p", "m–p", "mprice"),
}

#: A parenthetical aside inside a table cell — `1.415 *(search 1.411)*`, `1.34
#: (quad8)`.  The row's result is what the cell states outright; a number the
#: cell puts in parentheses is a second quantity the page names there, and the
#: row label does not describe it.
PAREN = re.compile(r"\([^)]*\)")


def _identity(key, kv):
    """Identity slots of one locked value: what its row label must name.

    The method is normally in the KEY — a LEM tag locks one value per method
    as ``fs_bishop=0.985, fs_spencer=…`` — and otherwise in ``method=``.  The
    engine, the element type and a constrained search add their own slot, so
    two SSRM locks on the same model are told apart by the mesh they name and a
    tangency-constrained search by the word "tangent".
    """
    slots = []
    t = kv.get("type", "").lower()
    if "ssrm" in t or t == "fem_elements":
        slots.append(("ssrm", "srf", "strength reduction"))
    m = kv.get("method", "").lower()
    km = re.fullmatch(r"fs_(\w+)", key)
    if km:
        m = km.group(1).lower()
    if m:
        slots.append(METHOD_WORDS.get(m, (m,)))
    et = kv.get("element_type", "").lower()
    if et:
        slots.append((et,))
    if "tangent_depth" in kv:
        slots.append(("tangent",))
    return slots


def _names(slots, text):
    """True when `text` names every slot of a lock's identity."""
    return bool(slots) and all(any(a in text for a in slot) for slot in slots)


def _agrees(tok, lock, dp):
    """True when a printed number restates `lock` at a precision it allows.

    Either the printed string is one of the forms the lock may be rounded to
    (the same rule the forward pass matches with), or it is inside the tag's own
    tolerance — which is the tag's own statement of what counts as this value.
    """
    value, tol = lock[0], lock[1]
    if tok in _forms(value, dp):
        return True
    try:
        return abs(Decimal(tok) - Decimal(value)) <= tol
    except Exception:
        return False


#: A cross-page reference to another verification section: `(rocscience.md#vp52)`.
XREF = re.compile(r"\(([a-z_0-9]+\.md)#([\w.:-]+)\)")


def _corpus_tags():
    """Every verification page's tags, keyed by input file and by anchor.

    A page routinely restates a value another page locks — geostudio's ACADS
    rows are Rocscience's models, tagged there — so a section that names the
    file gets the benefit of the lock wherever the tag lives.  Same scope the
    untagged sweep gives a corpus lock, and the same reason: the number is
    tag-guarded either way, and it is the restatement that can drift.

    A section also gets the locks of the section it cross-references by anchor
    ("**Rocscience detail:** [VP52](rocscience.md#vp52)"), because that link is
    the page saying the two present the same problem.  A pair of sections that
    tabulate the same model under different variants normally links only one of
    the input files, and the variant it publishes is locked in the other's.
    """
    import glob
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(here))
    by_file, by_anchor = {}, {}
    for page in sorted(glob.glob(os.path.join(repo, "docs", "verification",
                                              "*.md"))):
        name = os.path.basename(page)
        lines = open(page).read().split("\n")
        heads = [i for i, l in enumerate(lines) if re.match(r"^#{2,3} ", l)]
        secs = [(i, heads[k + 1] if k + 1 < len(heads) else len(lines))
                for k, i in enumerate(heads)]
        tagged = [(i, kv) for i, l in enumerate(lines) if (kv := _tag_kv(l))
                  and "file" in kv]
        bank = [(i, kv) for i, kv in tagged
                if not any(a <= i < b for a, b in secs)]
        for i, kv in tagged:
            by_file.setdefault(os.path.basename(kv["file"]), []).append(
                (kv, f"{name}:{i + 1}"))
        for a, b in secs:
            m = re.search(r"\{#([\w.:-]+)\}", lines[a])
            if not m:
                continue
            body = "\n".join(lines[a:b])
            here_tags = [(kv, f"{name}:{i + 1}") for i, kv in tagged
                         if a <= i < b]
            here_tags += [(kv, f"{name}:{i + 1}") for i, kv in bank
                          if kv["file"] in body]
            by_anchor[f"{name}#{m.group(1)}"] = here_tags
    return by_file, by_anchor


def _section_locks(lines, sec, page_bank, corpus):
    """Every lock in scope for one section, with its identity and its tag line.

    A section's own tags, a page-level bank tag whose input file the section
    links, any other page's tag on a file this section links, and the tags of a
    section this one cross-references by anchor.
    """
    by_file, by_anchor = corpus
    out = []
    seen = set()

    def add(kv, ln):
        if id(kv) in seen:
            return
        seen.add(id(kv))
        tol = kv.get("tolerance", "0.01")
        try:
            tol = Decimal(str(tol))
        except Exception:
            tol = Decimal("0.01")
        for k, v in kv.items():
            if not _wanted(k, ("expected_fs*", "fs_*", "expected")):
                continue
            slots = _identity(k, kv)
            for part in str(v).split(";"):
                part = part.split(":")[-1].strip()
                try:
                    Decimal(part)
                except Exception:
                    continue
                out.append((part, tol, slots, ln, kv.get("benchmark", "?"), k))

    body = "\n".join(lines[sec[0]:sec[1]])
    for i in range(sec[0], sec[1]):
        kv = _tag_kv(lines[i])
        if kv:
            add(kv, i + 1)
    for kv, ln in page_bank:
        if kv.get("file") and kv["file"] in body:
            add(kv, ln)
    for base, tagged in by_file.items():
        if base in body:
            for kv, ln in tagged:
                add(kv, ln)
    for page, anchor in XREF.findall(body):
        for kv, ln in by_anchor.get(f"{page}#{anchor}", ()):
            add(kv, ln)
    return out


def _table_rows(lines, sec):
    """(row line, label, [(cell text, start, end)]) for XSLOPE columns.

    A results table is one with a column headed XSLOPE; those are the cells that
    hold what the page computed, and the row's first cell says what quantity.

    Read only where the table has ONE such column.  A table with several —
    `XSLOPE composite` beside `XSLOPE circles-only`, `XSLOPE FE seepage` beside
    `XSLOPE piezo line` — publishes a row's method under several variants, and
    the row label names none of them; which lock belongs to which column is a
    reading of the header, not something the row states.  Those tables are left
    to the delta check, which pairs each column against the source beside it.
    """
    from .deltas import XCOL
    out, i = [], sec[0]
    while i < sec[1] - 1:
        if not lines[i].lstrip().startswith("|") or \
                not re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
            i += 1
            continue
        hdr = [c.strip() for c in lines[i].strip().strip("|").split("|")]
        xk = [k for k, c in enumerate(hdr) if XCOL.search(c)]
        j = i + 2
        while j < sec[1] and lines[j].lstrip().startswith("|"):
            if len(xk) == 1:
                bars = [p for p, ch in enumerate(lines[j]) if ch == "|"]
                cells = [(lines[j][a + 1:b], a + 1, b)
                         for a, b in zip(bars, bars[1:])]
                if len(cells) > xk[0]:
                    out.append((j, cells[0][0].lower(), [cells[xk[0]]]))
            j += 1
        i = max(j, i + 1)
    return out


def restatements(path, cfg, report=print):
    """Every number a section attributes to itself must agree with a lock.

    The forward pass is existential.  This one is universal over the numbers the
    section CLAIMS as its own results, which is a narrower set than "every
    number shaped like a factor of safety": a cell in the XSLOPE column of a row
    whose label names a locked method, or a number a sentence introduces with
    "XSLOPE's" / "this model's".  Anything the page attributes to a source, and
    anything it attributes to nobody, is left where it already is — to the
    delta check and to the untagged sweep.
    """
    from .deltas import AUTH_HDR_BASE, mask_sci
    from .untagged import (COORD_FIRST, COORD_SECOND, FS_HI, FS_LO, FS_SHAPED,
                           LABEL_BEFORE, SEEPAGE_TYPES, UNIT_AFTER, _mask)

    source = re.compile(
        SOURCE_EXTRA + "|" + AUTH_HDR_BASE
        + ("|" + "|".join(cfg.auth_hdr_extra) if cfg.auth_hdr_extra else ""),
        re.I)
    raw = open(path).read().replace(MINUS, "-")
    lines = raw.split("\n")
    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{2,3} ", l)]
    bounds = [(i, heads[k + 1][0] if k + 1 < len(heads) else len(lines), l)
              for k, (i, l) in enumerate(heads)]
    page_bank = [(kv, i + 1) for i, l in enumerate(lines)
                 if (kv := _tag_kv(l)) and not any(b[0] <= i < b[1]
                                                   for b in bounds)]
    corpus = _corpus_tags()

    #: One scan for the four things the prose arm reacts to, in priority order:
    #: an XSLOPE-side attributor, a source, the end of the sentence carrying
    #: them, and a number.  Each alternative is wrapped, so the combined pattern
    #: has exactly one group whichever alternative fired.
    prose = re.compile("|".join(f"(?:{p})" for p in (
        XSIDE.pattern, source.pattern, SENT_END.pattern, FS_SHAPED.pattern)),
        re.I)

    def candidate(line, tok, start, end):
        """The standard reasons a factor-of-safety-shaped token is not one."""
        try:
            v = Decimal(tok)
        except Exception:
            return False
        if not (FS_LO <= v <= FS_HI):
            return False
        after = line[end:end + 12]
        before = line[max(0, start - 24):start]
        if UNIT_AFTER.match(after) or LABEL_BEFORE.search(before):
            return False
        if COORD_SECOND.search(before) or (
                before.rstrip().endswith("(") and COORD_FIRST.match(after)):
            return False
        return True

    problems, checked = [], 0
    for sec in bounds:
        locks = _section_locks(lines, sec, page_bank, corpus)
        if not locks:
            continue
        types = {kv.get("type", "") for i in range(sec[0], sec[1])
                 if (kv := _tag_kv(lines[i]))}
        types |= {kv.get("type", "") for kv, _ in page_bank
                  if kv.get("file", "\x00") in "\n".join(lines[sec[0]:sec[1]])}
        if types and types <= SEEPAGE_TYPES:
            continue          # heads and flow rates, not factors of safety

        def flag(i, tok, why, bound):
            near = min(bound, key=lambda l: abs(Decimal(tok) - Decimal(l[0])))
            problems.append((i + 1, tok, why, near[0], near[3], near[4],
                             near[5]))

        # -- results-table cells: the XSLOPE column of a row the locks name ---
        table_rows = _table_rows(lines, sec)
        for j, label, cells in table_rows:
            bound = [l for l in locks if _names(l[2], label)]
            if not bound:
                continue
            for text, a, _b in cells:
                masked = PAREN.sub(lambda p: " " * len(p.group(0)),
                                   mask_sci(_mask(text)))
                for m in FS_SHAPED.finditer(masked):
                    tok = m.group(1)
                    if not candidate(text, tok, m.start(1), m.end(1)):
                        continue
                    checked += 1
                    if not any(_agrees(tok, l, cfg.tag_round_dp)
                               for l in bound):
                        flag(j, tok, f"row {label.strip()!r}", bound)

        # -- prose: a number the sentence attributes to XSLOPE ---------------
        rows = {j for j, _, _ in table_rows}
        side, sent = None, ""
        for i in range(sec[0], sec[1]):
            line = lines[i]
            if not line.strip() or line.lstrip().startswith("<!-- ") \
                    or line.startswith("#") or line.lstrip().startswith("|"):
                side, sent = None, ""
                continue
            if i in rows:
                continue
            masked = mask_sci(_mask(line))
            pos = 0
            for m in prose.finditer(masked):
                sent += masked[pos:m.start()]
                pos = m.end()
                tok = m.group(0)
                if SENT_END.fullmatch(tok):
                    side, sent = None, ""
                    continue
                if XSIDE.fullmatch(tok):
                    side, sent = "x", sent + tok
                    continue
                if source.fullmatch(tok):
                    side, sent = None, sent + tok
                    continue
                sent += tok
                if side != "x" or m.group(1) is None \
                        or not candidate(line, tok, m.start(), m.end()):
                    continue
                low = sent.lower()
                narrowed = [l for l in locks if _names(l[2], low)]
                bound = narrowed or locks
                checked += 1
                if not any(_agrees(tok, l, cfg.tag_round_dp) for l in bound):
                    flag(i, tok, "prose attributed to XSLOPE", bound)
            sent += masked[pos:]

    report(f"restatements: {checked} attributed numbers read; "
           f"disagreeing with their lock: {len(problems)}")
    for ln, tok, why, val, tagln, bench, key in problems:
        report(f"  L{ln} prints {tok} ({why}) but the lock is {key}={val} "
               f"[{bench}], tagged at L{tagln}")
    return len(problems)


def reverse(path, cfg, report=print):
    """Values presented as locked that no tag backs."""
    if not cfg.locked_value_re:
        report("reverse sweep: page carries no dotted summary table (skipped)")
        return 0
    lines = open(path).read().split("\n")
    tagvals = set()
    for l in lines:
        kv = _tag_kv(l)
        if kv is None:
            continue
        for k, v in kv.items():
            if not _wanted(k, cfg.tag_value_keys):
                continue
            # a list lock backs each of its elements
            tagvals.update(_elements(v) if ";" in v else [v])
    XV = re.compile(cfg.locked_value_re)
    insum = False
    bad = []
    for i, l in enumerate(lines):
        if cfg.summary_marker in l:
            insum = True
            continue
        if insum and l.strip() == "</div>":
            insum = False
            continue
        if not insum or not l.startswith("| ["):
            continue
        c = l.strip().strip("|").split("|")
        if len(c) < 4 or not DOT.search(c[1]):
            continue                       # nodata rows carry no lock claim
        # only the XSLOPE side of each ` · ` case: the part before " vs ",
        # since everything after it is the source being compared against.
        xside = " ".join(seg.split(" vs ")[0] for seg in c[3].split("·"))
        for v in XV.findall(xside):
            if v not in tagvals:
                bad.append((i + 1, v, c[0].strip(), c[3].strip()[:110]))
    report(f"reverse sweep: dotted summary rows checked; "
           f"values presented as locked with no tag: {len(bad)}")
    for b in bad:
        report(f"  L{b[0]} {b[2]} value {b[1]} has no tag :: {b[3]}")
    return len(bad)


def run(path, cfg, report=print):
    return (forward(path, cfg, report) + restatements(path, cfg, report)
            + reverse(path, cfg, report))


def _cli():
    from .pages import config_for
    path = sys.argv[1]
    return run(path, config_for(path))


if __name__ == "__main__":
    sys.exit(_cli())
