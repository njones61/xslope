#!/usr/bin/env python3
"""Tag-vs-text audit, BIDIRECTIONAL.

Forward : every locked value in a test tag (`expected_fs`, `expected_kc`,
          `expected_flowrate`, `fs_*`, ...) must appear verbatim in the section
          that carries the tag.  A tag value that is a semicolon LIST —
          ``expected=1.686;0.941;...`` (a factor of safety per march step),
          ``points=20:5:7.166;...`` (a solved head per station) — locks one
          value per element, and each element is checked the same way, so a
          wrong digit in one printed cell of an eleven-row table fails.
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
    return forward(path, cfg, report) + reverse(path, cfg, report)


def _cli():
    from .pages import config_for
    path = sys.argv[1]
    return run(path, config_for(path))


if __name__ == "__main__":
    sys.exit(_cli())
