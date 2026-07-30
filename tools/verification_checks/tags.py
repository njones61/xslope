#!/usr/bin/env python3
"""Tag-vs-text audit, BIDIRECTIONAL.

Forward : every locked value in a test tag (`expected_fs`, `expected_kc`,
          `expected_flowrate`, `fs_*`, ...) must appear verbatim in the section
          that carries the tag.
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

#: A number printed in the pages' scientific notation, `2.500×10⁻⁵`.
SCI_PRINTED = re.compile(r"(\d+(?:\.\d+)?)\s?×\s?10\s?(⁻|⁺)?([⁰¹²³⁴⁵⁶⁷⁸⁹]+)")
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")


def _sci_match(value, body):
    """True when `body` prints `value` in scientific notation.

    A discharge tag reads ``expected_flowrate=2.2985e-05`` while the page
    prints ``2.299×10⁻⁵``: the same number, at the precision the page reports
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
    lines = open(path).read().split("\n")
    heads = [(i, l) for i, l in enumerate(lines) if re.match(r"^#{2,3} ", l)]
    bounds = []
    for k, (i, l) in enumerate(heads):
        j = heads[k + 1][0] if k + 1 < len(heads) else len(lines)
        bounds.append((i, j, l))
    problems, rows, fired = [], [], set()

    def forms(v):
        """The strings that satisfy a tag value: the value itself, plus — where
        the page config allows it — the same value correctly rounded to the
        precision the page's prose uses."""
        out = [v]
        dp = cfg.tag_round_dp
        if dp is not None and "." in v and len(v.split(".")[1]) > dp:
            try:
                out.append(str(Decimal(v).quantize(Decimal(1).scaleb(-dp),
                                                   rounding=ROUND_HALF_UP)))
            except Exception:
                pass
        return out

    for idx, line in enumerate(lines):
        kv = _tag_kv(line)
        if kv is None:
            continue
        # A tag placed before the first heading is a page-level tag bank, not a
        # section tag: its scope is the whole page.
        sec = next((b for b in bounds if b[0] <= idx < b[1]),
                   (0, len(lines), "(page-level tag bank)"))
        body = "\n".join(l for l in lines[sec[0]:sec[1]]
                         if not l.lstrip().startswith("<!-- test:"))
        vals = {k: v for k, v in kv.items() if _wanted(k, cfg.tag_value_keys)}
        for k, v in vals.items():
            found = any(re.search(r"(?<![\w.])" + re.escape(f) + r"(?![\w.])",
                                  body) for f in forms(v))
            if not found and ("e" in v or "E" in v):
                found = _sci_match(v, body)
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
    report(f"tags scanned: {len({r['line'] for r in rows})}  "
           f"values checked: {len(rows)}  missing from their section: "
           f"{len(problems)}  dead tag exemptions: {len(dead)}")
    for p in problems:
        report(f"  L{p[0]} [{p[1]}] {p[2]}={p[3]}  -- not in {p[4]!r}")
    for val, key in dead:
        report(f"  DEAD tag exemption (never fires): {val!r}, {key!r}")
    return len(problems) + len(dead)


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
            if _wanted(k, cfg.tag_value_keys):
                tagvals.add(v)
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
