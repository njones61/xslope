#!/usr/bin/env python3
"""Authoritative delta re-derivation for a verification page.

Every printed comparison must satisfy

    delta = (XSLOPE - source) / source * 100,  Decimal ROUND_HALF_UP, 1 dp

from two numbers the page prints in the SAME sentence or the SAME table row.

What is checked
---------------
* percentages in every printed form: `(+1.2%)`, `+1.2%`, `(-7%)`, `1.23%`,
  `n %` (spaced), `**+1.2%**`;
* absolute factor-of-safety differences (`worth +0.013`, `costs 0.18`);
* percentages spelled in words (`two and a half percent`) — reported for
  adjudication, since no arithmetic can be attached to them;
* the precision rule: a signed percentage printed with other than one decimal
  is reported unless it is explicitly hedged (`about`, `~`, `roughly`, a
  range), and percentages that are not comparisons at all (domain shares, mesh
  fractions) are excluded by an explicit pattern rather than by a bare word.

How a pair is chosen — SIGN-AWARE, nearest first
------------------------------------------------
Scopes are tried in order and the FIRST that reproduces the printed value wins:

  cell     inside one table cell: source = the last number before the delta,
           XSLOPE = each earlier number in the cell, nearest first.
  matrix   the table header names an XSLOPE column: source = the number in the
           delta's own cell, XSLOPE = that row's XSLOPE cell.
  row      other cells of the same table row, nearest column first.
  sentence in prose: source = the last number before the delta in the sentence,
           XSLOPE = each earlier number in that sentence, nearest first.

For a SIGNED percentage only (first - second)/second is computed in the order
the page prints the two numbers, so a sign flip cannot be certified.  (An
UNSIGNED percentage states a magnitude and the direction is in the words, so
both orientations are admissible there.)  No scope enumerates all ordered pairs
of a paragraph.  Cross-scope references must be named in the page config's
whitelist with both numbers spelled out; the script still does the arithmetic,
so a wrong whitelist entry fails like anything else.

Usage: python -m tools.verification_checks.deltas <page.md> [out.json]
"""
import json
import re
import sys
from decimal import Decimal, ROUND_HALF_UP

MINUS = "−"
#: A number as the pages print one, including thousands separators
#: (`56,020` is one value, not "56" followed by "020").
NUM = re.compile(r"(?<![\w.,])(\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)(?![\w%]|\.\d|,\d)")
PCT = re.compile(r"(\(\s*)?(?:\*\*)?([+−-]?)(\d+(?:\.\d+)?)\s?%")
XCOL = re.compile(r"XSLOPE|X-?SLOPE", re.I)
HEDGE = re.compile(r"(?:about|around|~|roughly|approximately|within|under|over|"
                   r"below|above|≈|to|and|or|by|–|-)\s*$")
WORDED = re.compile(
    r"\b(?:one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|"
    r"twenty|thirty|forty|fifty|a|two and a half|one and a half|half|quarter|"
    r"third)\s+(?:and a half\s+)?(?:a\s+)?per\s?cent(?:age point)?s?\b", re.I)
ABS = re.compile(r"\b(?:worth|costs?|adds?|buys?|lifts? it by|moves? it by|"
                 r"higher by|lower by)\s+\**([+−-]?\d+\.\d+)\**(?![%\d]|\s*→)")
WITHIN = re.compile(r"\bwithin\s+(?:≈\s*)?$")
LT = re.compile(r"[<≤]\s*$")

#: Percentages that are shares of something rather than comparisons of two
#: printed factors.  Phrases only, never a bare word — a bare word would let a
#: real comparison hide behind an incidental noun.
NOT_A_COMPARISON_BASE = (
    r"of the domain|of the mesh|of the model domain|"
    r"of the model|by area|by count|by the vendor|by materials|by the polygon|"
    r"by the polygons|of the elements|Mohr-Coulomb element fraction|"
    r"element-area fraction|element fraction|of the section|of its pixels|"
    r"and draw a polygon|is held at full strength|is reduced|by element count|"
    r"of `?#?\w+`?'s (?:own )?(?:elements|domain|mesh|section|area)")

#: Column headers naming an authority: a delta in another cell of the same row
#: pairs that cell's value against them.
AUTH_HDR_BASE = (
    r"RS2|Slide2|referee|published|D&W|Duncan|USACE|Teoman|Loukidis|Cheng|"
    r"Li \(|Perry|Low|Nakamura|Hammah|PLAXIS|Z-Soil|GEO FEM|Baker|Greco|Zhu|"
    r"L&H|ACADS|Giam|Pilot|Prandtl|closed form|Fredlund|Tandjiria|Borges|"
    r"Wolff|Hassan|Ng & Shi|Huang")

SHARE_HDR = re.compile(r"held elastic|share|fraction|% of|Domain", re.I)
BREAK = re.compile(r"(?<=[a-z\)\]%*’])[.;]\s|\s·\s|\s—\s")
LINK = re.compile(r"\]\(#([\w:.-]+)\)")


def fmt(v, dp=1):
    q = v.quantize(Decimal(1).scaleb(-dp), rounding=ROUND_HALF_UP)
    zero = ("0." + "0" * dp) if dp else "0"
    return zero if q == 0 else ("+" if q > 0 else MINUS) + str(abs(q))


def pct(x, s):
    return (Decimal(x) - Decimal(s)) / Decimal(s) * 100


#: The exponent of a number in scientific notation — `×10⁻⁵`, `x 10^-5`.
#: Its "10" is not a value the page is comparing anything against, and left in
#: place it becomes "the last number before the delta" and poisons every pair on
#: a page that reports discharges.  Masking preserves string length so every
#: match offset stays valid.
SCI = re.compile(r"[×x]\s?10\s?(?:[⁻⁺][⁰¹²³⁴⁵⁶⁷⁸⁹]+|[⁰¹²³⁴⁵⁶⁷⁸⁹]+|\^?[-+]?\d+)")


def mask_sci(s):
    return SCI.sub(lambda m: " " * len(m.group(0)), s)


def nums(s):
    return [m.group(1).replace(",", "")
            for m in NUM.finditer(mask_sci(s))]


def _hdr_cols(lines, pick):
    """Map table-body line index -> column indices chosen by `pick(header)`."""
    out, i = {}, 0
    while i < len(lines):
        if lines[i].lstrip().startswith("|") and i + 1 < len(lines) \
                and re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
            hdr = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            cols = pick(hdr)
            j = i + 2
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                out[j] = cols
                j += 1
            i = j
        else:
            i += 1
    return out


def share_cols(lines, share_hdr=None):
    """Columns whose header declares a share/fraction, not a comparison."""
    share_hdr = share_hdr or SHARE_HDR
    return _hdr_cols(lines, lambda hdr: [k for k, c in enumerate(hdr)
                                         if share_hdr.search(c)])


def auth_map(lines, auth_hdr):
    """For tables with NO XSLOPE column: the columns naming an authority.
    A delta in another cell pairs that cell's value against them."""
    def pick(hdr):
        if any(XCOL.search(c) for c in hdr):
            return []
        return [k for k, c in enumerate(hdr) if auth_hdr.search(c)]
    return _hdr_cols(lines, pick)


def xcol_map(lines):
    return _hdr_cols(lines, lambda hdr: [k for k, c in enumerate(hdr)
                                         if XCOL.search(c)])


def sentence_bounds(text, pos):
    starts = [0] + [m.end() for m in BREAK.finditer(text) if m.end() <= pos]
    ends = [m.start() for m in BREAK.finditer(text) if m.start() > pos] + [len(text)]
    return max(starts), min(ends)


def table_pairs(lines, xmap, upto):
    """Structural (XSLOPE, source) pairs from the nearest matrix table above
    `upto` in the same section: every XSLOPE cell against every other cell of
    its OWN row.  Not an all-pairs enumeration of the block."""
    start = 0
    for k in range(upto - 1, -1, -1):
        if re.match(r"^#{2,4} ", lines[k]):
            start = k
            break
    pairs = []
    for j in range(start, upto):
        if j not in xmap or not xmap[j]:
            continue
        cells = lines[j].strip().strip("|").split("|")
        for xc in xmap[j]:
            if xc >= len(cells):
                continue
            for x in nums(cells[xc]):
                for k, c in enumerate(cells):
                    if k == xc:
                        continue
                    for s in nums(c):
                        pairs.append((x, s))
    return pairs


def paragraphs(lines):
    """(text, [(char_offset, line_index)]) for each prose paragraph."""
    out, buf, map_ = [], [], []
    for i, l in enumerate(lines):
        prose = l.strip() and not l.lstrip().startswith(("|", "<!-- test:",
                                                         "!["))
        if prose:
            map_.append((sum(len(b) + 1 for b in buf), i))
            buf.append(l)
        else:
            if buf:
                out.append((" ".join(buf), map_))
            buf, map_ = [], []
    if buf:
        out.append((" ".join(buf), map_))
    return out


def _pct_units(text):
    """Yield (match, sentence_start, sentence_end) for percentages in `text`."""
    for m in PCT.finditer(text):
        a, b = sentence_bounds(text, m.start())
        yield m, a, b


def _cands_prose(sent, rel):
    """Sign-aware candidates inside one sentence.

    forward  : source = last number before the delta, XSLOPE = each earlier
               number, nearest first        ("X vs S (d)")
    backward : XSLOPE = last number before the delta, source = each following
               number, nearest first        ("X, d on S")
    """
    sent = mask_sci(sent)
    pre = [(mm.group(1).replace(",", ""), mm.start())
           for mm in NUM.finditer(sent[:rel])]
    post = [mm.group(1).replace(",", "") for mm in NUM.finditer(sent[rel:])]
    out = []
    for s, _ in list(reversed(pre))[:3]:
        for x, _ in reversed(pre):
            if x is not s:
                out.append((x, s, "sentence-fwd"))
    for x, _ in list(reversed(pre))[:2]:
        for s in post:
            out.append((x, s, "sentence-back"))
    if not pre:
        # the claim leads its block ("worth 2.8% here ... 1.577 -> 1.534"):
        # the nearest pair after it is the one it is about
        for i in range(min(len(post), 4)):
            for j in range(min(len(post), 4)):
                if i != j:
                    out.append((post[i], post[j], "sentence-lead"))
    return out


def anchor_spans(lines):
    """anchor -> (start, end) line range of the section it names."""
    marks = []
    for i, l in enumerate(lines):
        m = re.match(r"^#{2,4} .*\{#([\w:.-]+)\}\s*$", l) or \
            re.match(r'^<a id="([^"]+)"></a>', l)
        if m:
            marks.append((i, m.group(1)))
    out = {}
    for k, (i, a) in enumerate(marks):
        j = marks[k + 1][0] if k + 1 < len(marks) else len(lines)
        out[a] = (i, j)
    return out


def reach(lines, spans, own, ctx):
    """The text a claim may draw its operands from: the block it sits in (its
    own paragraph or table row) plus every section it explicitly links to.
    Test tags are excluded — they are machine metadata, not printed content, so
    an operand surviving only in a tag does not keep a claim checkable."""
    def visible(a, b):
        return "\n".join(l for l in lines[a:b]
                         if not l.lstrip().startswith("<!-- test:"))
    txt = own
    for a in LINK.findall(ctx + " " + own):
        if a in spans:
            s, e = spans[a]
            txt += "\n" + visible(s, e)
    return txt


def section_of(lines, idx):
    for k in range(idx, -1, -1):
        if re.match(r"^#{2,4} ", lines[k]):
            return k
    return 0


def run(path, cfg, out_json=None, report=print):
    """Re-derive every printed comparison on `path` under `cfg`.

    Returns (failures, summary_dict).  `report` receives one line at a time.
    """
    not_a_comparison = re.compile(
        r"\d+(?:\.\d+)?\s?%\**\s+(?:" + NOT_A_COMPARISON_BASE
        + ("|" + "|".join(cfg.not_a_comparison_extra)
           if cfg.not_a_comparison_extra else "") + r")")
    # A list of shares sharing one "of ..." tail — "0.13%, 0.31% and 0.14% of
    # u0".  Only the last member is followed directly by the phrase, so the
    # earlier ones are recognised here, anchored immediately after their own
    # percent sign so an unrelated share later in the line cannot excuse them.
    share_list = re.compile(
        r"^\**(?:\s*[,;]?\s*(?:and\s+)?\d+(?:\.\d+)?\s?%\**)+\s+(?:"
        + NOT_A_COMPARISON_BASE
        + ("|" + "|".join(cfg.not_a_comparison_extra)
           if cfg.not_a_comparison_extra else "") + r")")

    prefix_re = (re.compile(r"(?:" + "|".join(cfg.not_a_comparison_prefix)
                            + r")$")
                 if cfg.not_a_comparison_prefix else None)

    def is_share(txt, start, end):
        if prefix_re is not None and prefix_re.search(txt[max(0, start - 48):start]):
            return True
        return bool(not_a_comparison.search(txt[start:start + 70])
                    or share_list.match(txt[end:end + 110]))

    auth_hdr = re.compile(
        AUTH_HDR_BASE + ("|" + "|".join(cfg.auth_hdr_extra)
                         if cfg.auth_hdr_extra else ""), re.I)
    WHITELIST, BOUNDS = cfg.whitelist, cfg.bounds
    ABS_BOUNDS, WORDED_OK = cfg.abs_bounds, cfg.worded_ok

    lines = open(path).read().split("\n")
    xmap = xcol_map(lines)
    scmap = share_cols(lines, re.compile(
        SHARE_HDR.pattern + ("|" + "|".join(cfg.share_hdr_extra)
                             if cfg.share_hdr_extra else ""), re.I))
    amap = auth_map(lines, auth_hdr)
    results, notes = [], []
    fired = set()          # whitelist / bounds entries that actually execute
    orphaned = []          # whitelist entries whose operands left the page
    spans = anchor_spans(lines)
    heads = [i for i, l in enumerate(lines) if re.match(r"^#{2,3} ", l)] + [len(lines)]

    def reachable(ctx, sec, para=""):
        own = para or ctx
        if sec is not None:
            s = max([h for h in heads if h <= sec] or [0])
            e = next((h for h in heads if h > s), len(lines))
            own += "\n" + "\n".join(l for l in lines[s:e]
                                    if not l.lstrip().startswith("<!-- test:"))
        return reach(lines, spans, own, ctx)

    def num_in(v, text):
        return re.search(r"(?<![\w.,])" + re.escape(v) + r"(?![\d.,])", text) is not None

    derived = {}          # section start -> {printed delta already derived}

    def match(printed, x, s, unsigned, dp=1):
        got = fmt(pct(x, s), dp)
        if got == printed:
            return True
        return unsigned and got.lstrip("+" + MINUS) == printed.lstrip("+" + MINUS)

    def record_pct(ln, printed, dp, hedged, cand, ctx, structural=None,
                   sec=None, row_repeats=None, unsigned=False, in_table=False,
                   para_text="", lt=False):
        # every token is reconciled at the precision it is printed to
        quote_ok = True
        if lt:
            # `<0.01%` / `≤2%` states an UPPER BOUND on a difference the page
            # does print.  Where the page prints the pair, that is checkable
            # arithmetic rather than an exemption: the nearest non-degenerate
            # pair — the one the sentence is about — must differ by no more than
            # the printed bound.  If it does, the claim is certified here.  If
            # it does not, pair reconciliation is switched off (a bound is not a
            # delta, so no other pair may certify it) and the claim must be
            # named in BOUNDS, like any other bound over a set the page does not
            # print.
            want = Decimal(printed.lstrip("+" + MINUS))
            for x, s, how in cand:
                if Decimal(x) == Decimal(s):
                    continue                  # a value against itself bounds nothing
                try:
                    got = pct(x, s).copy_abs()
                except Exception:
                    continue
                if got <= want:
                    results.append(dict(
                        line=ln, printed="<" + printed.lstrip("+" + MINUS),
                        kind="upper-bound", ctx=ctx[:170], status="OK",
                        xslope=x, source=s, via=how + "/bound",
                        recomputed=str(got.quantize(Decimal("0.0001")))))
                    return
                break
            cand, structural, row_repeats = [], None, None
            quote_ok = False
        if dp != 1:
            notes.append(dict(line=ln, kind="precision", printed=printed,
                              signed=not unsigned, hedged=hedged, text=ctx[:170]))
        # A "within X%" bound is a claim about the WHOLE operand set, not about
        # one pair, so pair reconciliation cannot check it — some pair in a rich
        # section will always match.  Such tokens must be pinned by value in
        # BOUNDS, adjudicated once against the set they summarise.
        if WITHIN.search(ctx[:ctx.find(printed.lstrip("+" + MINUS) + "%")]
                         if printed.lstrip("+" + MINUS) + "%" in ctx else ""):
            cand, structural, row_repeats = [], None, None
            quote_ok = False
        hit = None
        for x, s, how in cand:
            try:
                if match(printed, x, s, unsigned, dp):
                    hit = (x, s, how)
                    break
            except Exception:
                pass
        if hit is None and structural:
            for x, s in structural:
                try:
                    if match(printed, x, s, unsigned, dp):
                        hit = (x, s, "table-above")
                        break
                except Exception:
                    pass
        if hit is None and row_repeats and printed in row_repeats:
            hit = ("-", "-", "row-repeat")
        if hit is None and sec is not None and not in_table and quote_ok \
                and printed in derived.get(sec, ()) and not cand:
            # only a sentence that prints no operands of its own may quote a
            # figure derived elsewhere in the section; a sentence with numbers
            # must reconcile from them
            hit = ("-", "-", "quotes-section-delta")
        if hit is None:
            for pd, key, xw, sw in WHITELIST:
                if pd == printed and key in ctx and match(printed, xw, sw, unsigned, dp):
                    # The tuple's arithmetic is not enough.  Both operands must
                    # STILL be printed where the claim can reach them — its own
                    # section, or a section the sentence links to — so that
                    # editing an operand cannot silently orphan a certified
                    # comparison.
                    scope = reachable(ctx, sec, para_text or ctx)
                    if not (num_in(xw, scope) and num_in(sw, scope)):
                        orphaned.append((pd, key, xw, sw, ln))
                        continue
                    fired.add((pd, key))
                    hit = (xw, sw, "whitelist")
                    break
        if hit is None:
            for pd, key in BOUNDS:
                if pd == printed and key in ctx:
                    fired.add((pd, key))
                    hit = ("-", "-", "bound/share")
                    break
        e = dict(line=ln, printed=printed, kind="percent", ctx=ctx[:170])
        if hit:
            e.update(status="OK", xslope=hit[0], source=hit[1], via=hit[2],
                     recomputed=(fmt(pct(hit[0], hit[1]), dp) if hit[0] != "-"
                                 else printed))
            if hit[0] != "-" and sec is not None:
                derived.setdefault(sec, set()).add(printed)
        else:
            near, seen = [], set()
            for x, s, how in cand[:16]:
                if (x, s) in seen:
                    continue
                seen.add((x, s))
                try:
                    near.append((x, s, fmt(pct(x, s), dp)))
                except Exception:
                    pass
            e.update(status="UNRECONCILED", candidates=near[:6])
        results.append(e)

    # ---------------- table rows ----------------
    for idx, line in enumerate(lines):
        if not line.lstrip().startswith("|"):
            continue
        if any(sub in line for sub in cfg.share_rows):
            continue          # a declared value row, not a row of comparisons
        ln = idx + 1
        cells = line.strip().strip("|").split("|")
        base = line.index("|") + 1
        for m in PCT.finditer(line):
            if is_share(line, m.start(), m.end()):
                continue
            acc, ci = base, None
            for k, c in enumerate(cells):
                if acc <= m.start() < acc + len(c):
                    ci = k
                    break
                acc += len(c) + 1
            if ci is not None and ci in (scmap.get(idx) or []):
                continue                      # a declared share column
            sign, mag = m.group(2), m.group(3)
            dp = len(mag.split(".")[1]) if "." in mag else 0
            printed = (fmt(Decimal(0), dp) if float(mag) == 0 else
                       (MINUS if sign in ("-", MINUS) else "+") + mag)
            own = cells[ci] if ci is not None else ""
            pre = own[: m.start() - acc] if ci is not None else ""
            srcs = nums(pre) or nums(own)
            cand = []
            # tables laid out the other way round: the source is a column of
            # its own and each value cell carries its delta against it
            for ac in (amap.get(idx) or []):
                if ci is not None and ac != ci and ac < len(cells):
                    for s2 in nums(cells[ac]):
                        for x2 in srcs:
                            cand.append((x2, s2, "authority-column"))
            if srcs:
                # against a printed range, either bound may be the one quoted.
                # The source is identified by POSITION, not by value: a token
                # paired against itself yields 0.0 and would certify a printed
                # "0%" no matter what the row actually says.
                pn = nums(pre)
                rng = len(pn) >= 2 and re.search(r"\d\s?[–-]\s?\d", pre)
                bidx = ([len(pn) - 1] + ([len(pn) - 2] if rng else [])
                        if pn else [None])
                for bi in bidx:
                    s = pn[bi] if bi is not None else srcs[-1]
                    for k2 in range(len(pn) - 1, -1, -1):
                        if k2 != bi:
                            cand.append((pn[k2], s, "cell"))
                    for xc in (xmap.get(idx) or []):
                        if ci is not None and xc != ci and xc < len(cells):
                            for x in nums(cells[xc]):
                                cand.append((x, s, "matrix"))
                s = srcs[-1]
                for k in range(len(cells) - 1, -1, -1):
                    if k != ci:
                        for x in nums(cells[k]):
                            cand.append((x, s, "row"))
            allp = [(fmt(Decimal(0), len(g3.split(".")[1]) if "." in g3 else 0)
                     if float(g3) == 0 else
                     (MINUS if g2 in ("-", MINUS) else "+") + g3)
                    for g1, g2, g3 in
                    ((mm.group(1), mm.group(2), mm.group(3))
                     for mm in PCT.finditer(line)) if g1]
            # a repeat only counts when the value is printed more than once in
            # the row, otherwise the token would certify itself
            rr = {v for v in allp if allp.count(v) > 1}
            record_pct(ln, printed, dp,
                       bool(HEDGE.search(line[:m.start()][-16:])),
                       cand, line.strip(), None, section_of(lines, idx), rr,
                       unsigned=not sign, in_table=True,
                       lt=bool(LT.search(line[:m.start()])))

    # ---------------- prose paragraphs ----------------
    for text, cmap in paragraphs(lines):
        def line_of(pos, cmap=cmap):
            best = cmap[0][1]
            for off, li in cmap:
                if off <= pos:
                    best = li
            return best + 1
        for m in WORDED.finditer(text):
            ctxw = text[max(0, m.start() - 110):m.start() + 90]
            if any(k in ctxw for k in WORDED_OK):
                continue
            notes.append(dict(line=line_of(m.start()), kind="worded-percent",
                              text=ctxw))
        for m in ABS.finditer(text):
            claimed = m.group(1)
            a, b = sentence_bounds(text, m.start())
            # sentence scope, nearest pair, sign-aware: the difference is taken
            # in the order the sentence prints the two numbers, exactly as the
            # percent path does.  No all-ordered-pairs pool.
            pre = nums(text[a:m.start()])
            post = nums(text[m.end():b])
            want = Decimal(claimed.replace(MINUS, "-"))
            signed = claimed[0] in "+-" + MINUS
            got = None
            pairs = []
            if len(pre) >= 2:
                pairs += [(pre[-1], pre[-2]), (pre[-2], pre[-1])][:1 if signed else 2]
            for s in post[:2]:
                if pre:
                    pairs.append((pre[-1], s))
                    if not signed:
                        pairs.append((s, pre[-1]))
            if not pre and len(post) >= 2:
                pairs += [(post[0], post[1]), (post[1], post[0])]
            for x, s in pairs:
                d = Decimal(x) - Decimal(s)
                if (d if signed else d.copy_abs()).quantize(Decimal("0.0001")) == \
                        (want if signed else want.copy_abs()):
                    got = (x, s)
                    break
            if got is None:
                for pd, key in ABS_BOUNDS:
                    if pd == claimed and key in text[a:b]:
                        fired.add((pd, key))
                        got = ("-", "-")
                        break
            if got:
                results.append(dict(line=line_of(m.start()), printed=claimed,
                                    kind="absolute", status="OK", pair=got,
                                    ctx=text[a:b].strip()[:170]))
            else:
                results.append(dict(line=line_of(m.start()), printed=claimed,
                                    kind="absolute", status="UNRECONCILED",
                                    candidates=[(x, s, str(Decimal(x) - Decimal(s)))
                                                for x, s in pairs[:4]],
                                    ctx=text[a:b].strip()[:170]))

        struct = table_pairs(lines, xmap, cmap[0][1])
        sent_starts = sorted({sentence_bounds(text, mm.start())[0]
                              for mm in PCT.finditer(text)})
        for m, a, b in _pct_units(text):
            if is_share(text, m.start(), m.end()):
                continue
            sign, mag = m.group(2), m.group(3)
            dp = len(mag.split(".")[1]) if "." in mag else 0
            printed = (fmt(Decimal(0), dp) if float(mag) == 0 else
                       (MINUS if sign in ("-", MINUS) else "+") + mag)
            hedged = bool(HEDGE.search(text[:m.start()][-16:]))
            cand = _cands_prose(text[a:b], m.start() - a)
            if not hedged:
                cand += [(x, s, "paragraph") for x, s, _ in
                         _cands_prose(text, m.start())]
            prior = [s for s in sent_starts if s < a]
            if prior:
                pa = prior[-1]
                ptxt = text[pa:a]
                cand += [(x, s, "prior-sentence")
                         for x, s, _ in _cands_prose(ptxt, len(ptxt))]
            record_pct(line_of(m.start()), printed, dp, hedged, cand,
                       text[a:b].strip(), struct, section_of(lines, cmap[0][1]),
                       unsigned=not sign, para_text=text,
                       lt=bool(LT.search(text[:m.start()])))

    bad = [r for r in results if r["status"] == "UNRECONCILED"]
    prec_all = [n for n in notes if n["kind"] == "precision"]
    for n in prec_all:
        n.setdefault("signed", False)
    # A signed delta must be printed to one decimal, unless it is
    # explicitly hedged ("to within −0.27%", "about 2%"), where the
    # extra digit is part of the hedge rather than a second convention.
    prec = [n for n in prec_all if n["signed"] and not n["hedged"]]
    word = [n for n in notes if n["kind"] == "worded-percent"]
    unp = [n for n in notes if n["kind"] == "absolute-unpaired"]
    dead = [(pd, key) for pd, key, *_ in WHITELIST if (pd, key) not in fired] \
        + [(pd, key) for pd, key in BOUNDS if (pd, key) not in fired] \
        + [(pd, key) for pd, key in ABS_BOUNDS if (pd, key) not in fired]
    ok = sum(1 for r in results if r["status"] == "OK")
    report(f"{path}: comparisons={len(results)} clean={ok} unreconciled={len(bad)}"
           f" | precision violations={len(prec)}/{len(prec_all)} worded={len(word)}"
           f" unpaired-absolute={len(unp)} dead-exemptions={len(dead)}"
           f" orphaned-operands={len(orphaned)}")
    for pd, key in dead:
        report(f"  DEAD exemption (never fires): {pd!r}, {key!r}")
    for pd, key, xw, sw, ln in orphaned:
        report(f"  L{ln} ORPHANED whitelist operand: {pd} needs {xw} and {sw}, "
               f"but they are no longer printed where the claim can reach them "
               f"({key!r})")
    for r in sorted(bad, key=lambda r: r["line"]):
        report(f"  L{r['line']} [{r['kind']}] printed {r['printed']} "
               f"cands={r.get('candidates') or r.get('pair')}")
        report(f"      {r['ctx']}")
    for n in prec_all:
        tag = ("precision" if n["signed"] and not n["hedged"] else
               "precision-hedged" if n["signed"] else "precision-note")
        report(f"  L{n['line']} [{tag}] {n['printed']} :: {n['text'][:130]}")
    for n in word:
        report(f"  L{n['line']} [worded] {n['text'][:150]}")
    for n in unp:
        report(f"  L{n['line']} [absolute-unpaired] {n['printed']} :: {n['text'][:130]}")
    if out_json:
        json.dump(dict(results=results, notes=notes), open(out_json, "w"),
                  indent=1, ensure_ascii=False)
    return (len(bad) + len(prec) + len(word) + len(dead) + len(orphaned),
            dict(comparisons=len(results), clean=ok, unreconciled=len(bad),
                 precision=len(prec), worded=len(word), dead=len(dead),
                 orphaned=len(orphaned)))


def _cli():
    from .pages import config_for
    path = sys.argv[1]
    fails, _ = run(path, config_for(path),
                   sys.argv[2] if len(sys.argv) > 2 else None)
    return fails


if __name__ == "__main__":
    sys.exit(_cli())
