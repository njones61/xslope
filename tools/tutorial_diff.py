"""Track-changes diff of a docs page between two git revisions, as HTML and PDF.

Usage: python3 tools/tutorial_diff.py <base-rev> <head-rev> <page.md> [<page.md> ...]
       --changes-only   print only changed blocks, with the heading above each
                        and one block of context on either side
       --out DIR        (default: xslope_private/reports/tutorial_diffs)

The page is rendered as prose -- headings, bold, code blocks and tables come out
as themselves, images as a gray "[figure: name]" line, and the lock comments are
dropped. Unchanged text is plain, removed words are red and struck through,
added words are green, bold and underlined. The PDF is produced by LibreOffice
headless (never Word).
"""
import datetime, difflib, html, os, re, subprocess, sys

import markdown

SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
OUT_DIR = "/Users/njones/python_projects/xslope_private/reports/tutorial_diffs"

# LibreOffice's HTML import ignores <style> element rules (and leaks del/ins
# decoration into the rest of the document when it half-reads them), so every
# rule below is applied inline instead. Only @page survives in a <style> block.
BODY = "font-family:Georgia,'Times New Roman',serif;font-size:11pt;color:#1a1a1a;line-height:1.45"
DEL_S = "color:#b00020;text-decoration:line-through"
INS_S = "color:#0a6d1f;text-decoration:underline;font-weight:bold"
FIG_S = "color:#808080;font-style:italic"
MONO = "'DejaVu Sans Mono','Courier New',monospace"

TAG_STYLE = {
    "h1": f"font-family:Georgia,serif;font-size:16pt;color:#000;margin:18pt 0 6pt 0",
    "h2": f"font-family:Georgia,serif;font-size:13.5pt;color:#000;margin:16pt 0 5pt 0",
    "h3": f"font-family:Georgia,serif;font-size:12pt;color:#000;margin:14pt 0 4pt 0",
    "h4": f"font-family:Georgia,serif;font-size:11pt;color:#000;margin:12pt 0 4pt 0",
    "h5": f"font-family:Georgia,serif;font-size:11pt;color:#000;margin:12pt 0 4pt 0",
    "h6": f"font-family:Georgia,serif;font-size:11pt;color:#000;margin:12pt 0 4pt 0",
    "p": BODY + ";margin:0 0 8pt 0",
    "li": BODY + ";margin:0 0 3pt 0",
    "blockquote": BODY + ";margin:0 0 8pt 18pt",
    "pre": f"font-family:{MONO};font-size:8.5pt;color:#222;margin:0",
    "code": f"font-family:{MONO};font-size:9.5pt",
    "th": "border:1px solid #999;padding:3px;font-family:Georgia,serif;font-size:9pt;text-align:left",
    "td": "border:1px solid #999;padding:3px;font-family:Georgia,serif;font-size:9pt",
}

MD_EXT = ["tables", "fenced_code", "attr_list", "sane_lists"]

# ---------------------------------------------------------------- source prep

FRONT = re.compile(r"\A---\n(.*?)\n---\n", re.S)
COMMENT = re.compile(r"<!--.*?-->", re.S)
IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+)[^)]*\)(\{[^}]*\})?")
PILL = re.compile(r'<span class="tg-pill">([^<]*)</span>')
TAGS = re.compile(r"<[^>]+>")
FIG_SPAN = re.compile(r'<span style="' + re.escape(FIG_S) + r'">\[figure: [^<]*\]</span>')


def prep(text):
    """Strip the front matter and lock tags; turn images into figure lines."""
    fm = ""
    m = FRONT.match(text)
    if m:
        fm, text = m.group(1), text[m.end():]
    text = COMMENT.sub("", text)
    text = IMAGE.sub(lambda m: f'<span style="{FIG_S}">[figure: '
                               f'{os.path.basename(m.group(1))}]</span>', text)
    text = PILL.sub(r"\1 · ", text)
    return fm, text


def split_blocks(text):
    """Blank-line blocks, with fenced code kept whole. Raw HTML furniture is
    flattened to its text so the words inside it still diff."""
    blocks, buf, fence = [], [], None
    for line in text.split("\n"):
        f = re.match(r"^\s*(```+|~~~+)", line)
        if fence:
            buf.append(line)
            if f and line.strip().startswith(fence):
                blocks.append("\n".join(buf)); buf, fence = [], None
            continue
        if f:
            if buf:
                blocks.append("\n".join(buf)); buf = []
            fence = f.group(1); buf = [line]
            continue
        if line.strip() == "":
            if buf:
                blocks.append("\n".join(buf)); buf = []
        else:
            buf.append(line)
    if buf:
        blocks.append("\n".join(buf))

    out = []
    for b in blocks:
        if is_fence(b):
            out.append(b); continue
        if b.lstrip().startswith("<") and not FIG_SPAN.match(b.lstrip()):
            out.extend(emphasis(x) for x in flatten_html(b))
            continue
        if b.strip() in ("---", "***", "___"):
            continue
        out.append(emphasis(b))
    return out


BOLD = re.compile(r"\*\*(?=\S)(.+?)(?<=\S)\*\*", re.S)
ITAL = re.compile(r"(?<![\w*])\*(?=\S)([^*\n]+?)(?<=\S)\*(?![\w*])")


def emphasis(b):
    """Resolve bold and italic to tags before the diff runs. A rewritten bold
    phrase cannot keep its ** markers balanced across both revisions, and an
    orphaned marker prints as literal asterisks; an orphaned tag prints as
    nothing."""
    return ITAL.sub(r"<i>\1</i>", BOLD.sub(r"<b>\1</b>", b))


HTML_LINE = re.compile(r"^\s*</?(div|p|table|tr|td|th|ul|ol|li|section)\b")


def flatten_html(b):
    """A glance box or pill row: keep the words, drop the markup, and let each
    tile or paragraph of it stand as its own block."""
    segs, cur = [], []
    for ln in b.split("\n"):
        if HTML_LINE.match(ln):
            if cur:
                segs.append("\n".join(cur)); cur = []
            segs.append(ln)
        else:
            cur.append(ln)
    if cur:
        segs.append("\n".join(cur))
    out = []
    for s in segs:
        t = TAGS.sub(" ", s)
        t = "\n".join(re.sub(r"[ \t]{2,}", " ", ln).strip() for ln in t.split("\n"))
        t = re.sub(r"\n{2,}", "\n", t).strip()
        if t:
            out.append(html.unescape(t))
    return out


def is_fence(b):
    return bool(re.match(r"^\s*(```|~~~)", b))


def is_heading(b):
    return bool(re.match(r"^\s*#{1,6}\s", b))


# --------------------------------------------------------------- word diffing

# Inline spans that must never be split by a <del>/<ins> boundary.
TOKEN = re.compile(
    r"\$\$.*?\$\$"          # display math
    r"|\$[^$\n]+\$"          # inline math
    r"|`[^`\n]*`"             # code span
    r"|<span[^>]*>.*?</span>"  # a figure line
    r"|\[[^\]\n]*\]\([^)\n]*\)"  # link
    r"|</?[bi]>"               # emphasis, resolved above
    r"|\*+|_{1,2}(?![\w])"    # any marker that survived, on its own
    r"|\s+"
    r"|[^\s*<]+"
    r"|\S+", re.S)

# Markdown structure that must stay outside a styled span. The list and heading
# markers only count at the start of a line -- "84337150." is a number, not an
# ordered-list marker.
BARE = re.compile(r"^(\||\*+|_{1,2}|</?[bi]>)$")   # never inside a styled span
STRUCT = re.compile(r"^(\||#{1,6}|[-*+]|\d{1,3}\.|>|:?-{3,}:?)$")


def toks(text):
    return TOKEN.findall(text)


def wrap(run, style, line_start=True):
    """Wrap a token run, letting newlines and structural markers through bare so
    lists, tables and headings keep their markdown meaning."""
    out, seg = [], []
    bol = line_start

    def flush():
        if seg:
            s = "".join(seg)
            out.append(f'<span style="{style}">{s}</span>' if s.strip() else s)
            seg.clear()

    for t in run:
        if "\n" in t:
            flush(); out.append(t); bol = True
        elif BARE.match(t) or (bol and STRUCT.match(t)):
            flush(); out.append(t)
        else:
            if t.strip():
                bol = False
            seg.append(t)
    flush()
    return "".join(out)


def nwords(run):
    return sum(1 for t in run if t.strip())


def coalesce(ops, a, b, limit=6):
    """A rewritten sentence comes back from SequenceMatcher as a dozen one-word
    changes threaded on the few words that happened to survive, which reads as
    salad. Swallow a short surviving run when the change around it is at least
    as long, so the sentence prints once struck and once whole."""
    merged, i = [], 0
    while i < len(ops):
        op = ops[i]
        if op[0] == "equal":
            merged.append(op); i += 1; continue
        _, i1, i2, j1, j2 = op
        changed = nwords(a[i1:i2]) + nwords(b[j1:j2])
        k = i + 1
        while k + 1 < len(ops) and ops[k][0] == "equal" and ops[k + 1][0] != "equal":
            e = nwords(a[ops[k][1]:ops[k][2]])
            if e > limit or changed < e:
                break
            nx = ops[k + 1]
            changed += e + nwords(a[nx[1]:nx[2]]) + nwords(b[nx[3]:nx[4]])
            i2, j2 = nx[2], nx[4]
            k += 2
        merged.append(("replace", i1, i2, j1, j2))
        i = k
    return merged


def merge(old, new):
    """Merged markdown for a modified block: word runs wrapped in del/ins."""
    a, b = toks(old), toks(new)
    out = []
    ops = difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes()
    if not old.lstrip().startswith("|"):     # a table diffs cell by cell
        ops = coalesce(ops, a, b)
    for op, i1, i2, j1, j2 in ops:
        if op in ("equal",):
            out.append("".join(a[i1:i2]))
        else:
            bol_a = i1 == 0 or "\n" in a[i1 - 1]
            bol_b = j1 == 0 or "\n" in b[j1 - 1]
            d = wrap(a[i1:i2], DEL_S, bol_a) if op in ("delete", "replace") else ""
            n = wrap(b[j1:j2], INS_S, bol_b) if op in ("insert", "replace") else ""
            # A struck word butted against its replacement is unreadable.
            tight = (d and n and a[i2 - 1].strip() and b[j1].strip())
            gap = " " if tight else ""
            out.append(d + gap + n)
    return "".join(out)


def fence_diff(old, new):
    """Code fences are literal, so the del/ins spans go into the HTML directly."""
    def body(b):
        lines = b.split("\n")
        return "\n".join(lines[1:-1] if len(lines) > 2 and re.match(r"^\s*(```|~~~)", lines[-1])
                         else lines[1:])
    def span(run, style):
        parts = []
        for chunk in re.split(r"(\n)", "".join(run)):
            parts.append(chunk if chunk == "\n" or not chunk.strip()
                         else f'<span style="{style}">{html.escape(chunk)}</span>')
        return "".join(parts)

    a, b_ = toks(body(old)), toks(body(new))
    out = []
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b_, autojunk=False).get_opcodes():
        if op == "equal":
            out.append(html.escape("".join(a[i1:i2])))
        else:
            d = span(a[i1:i2], DEL_S) if op in ("delete", "replace") else ""
            n = span(b_[j1:j2], INS_S) if op in ("insert", "replace") else ""
            gap = " " if d and n and a[i2 - 1].strip() and b_[j1].strip() else ""
            out.append(d + gap + n)
    return f'<pre style="{TAG_STYLE["pre"]}">' + "".join(out) + "</pre>"


# ----------------------------------------------------------------- rendering

def render(md_text):
    return markdown.Markdown(extensions=MD_EXT).convert(md_text)


def style_html(h):
    """Push the inline styles onto every tag LibreOffice will meet."""
    def rep(m):
        tag, rest = m.group(1).lower(), m.group(2)
        s = TAG_STYLE.get(tag)
        if not s or "style=" in rest:
            return m.group(0)
        return f"<{m.group(1)}{rest} style=\"{s}\">"
    h = re.sub(r"<(\w+)((?:\s[^>]*?)?)>", rep, h)
    # LibreOffice repeats a <thead> across a page break and prints the repeat blank.
    h = re.sub(r"</?(thead|tbody)>", "", h)
    h = h.replace("<table>", '<table border="1" cellspacing="0" cellpadding="3" '
                             'style="border-collapse:collapse;width:100%;margin:0 0 8pt 0">')
    return h


TEXT_SPLIT = re.compile(r"(<[^>]+>)")


def colorize(h, style):
    """Paint every text node of an added or removed block."""
    out = []
    for piece in TEXT_SPLIT.split(h):
        if piece.startswith("<") or not piece.strip():
            out.append(piece)
        else:
            out.append(f'<span style="{style}">{piece}</span>')
    return "".join(out)


def block_html(md_text, kind=None):
    if is_fence(md_text):
        lines = md_text.split("\n")
        inner = lines[1:-1] if len(lines) > 2 and re.match(r"^\s*(```|~~~)", lines[-1]) else lines[1:]
        h = f'<pre style="{TAG_STYLE["pre"]}">' + html.escape("\n".join(inner)) + "</pre>"
    else:
        h = style_html(render(md_text))
    if kind == "del":
        h = colorize(h, DEL_S)
    elif kind == "ins":
        h = colorize(h, INS_S)
    return h


# ------------------------------------------------------------------ pairing

def norm(b):
    return re.sub(r"\s+", " ", b).strip()


def ratio(o, n):
    return difflib.SequenceMatcher(None, norm(o).split(), norm(n).split(),
                                   autojunk=False).ratio()


def pair(olds, news):
    """Greedy best-similarity pairing inside a replace run."""
    pairs, used = {}, set()
    for i, o in enumerate(olds):
        best, score = None, 0.45
        for j, n in enumerate(news):
            if j in used or is_fence(o) != is_fence(n) or is_heading(o) != is_heading(n):
                continue
            r = ratio(o, n)
            if r > score:
                best, score = j, r
        if best is not None:
            pairs[i] = best; used.add(best)
    return pairs


def build(old_blocks, new_blocks):
    """A list of (kind, html, md_for_heading_tracking) in reading order."""
    items = []
    sm = difflib.SequenceMatcher(None, [norm(b) for b in old_blocks],
                                 [norm(b) for b in new_blocks], autojunk=False)
    ops = []
    for o in sm.get_opcodes():
        # A rewritten block often lands as a bare delete beside a bare insert;
        # fuse those so the pairing below can diff them word by word.
        if ops and o[0] != "equal" and ops[-1][0] != "equal":
            p = ops.pop()
            ops.append(("replace", p[1], o[2], p[3], o[4]))
        else:
            ops.append(o)
    for op, i1, i2, j1, j2 in ops:
        if op == "equal":
            for j in range(j1, j2):
                items.append(("same", new_blocks[j]))
        elif op == "delete":
            for i in range(i1, i2):
                items.append(("del", old_blocks[i]))
        elif op == "insert":
            for j in range(j1, j2):
                items.append(("ins", new_blocks[j]))
        else:
            olds, news = old_blocks[i1:i2], new_blocks[j1:j2]
            p = pair(olds, news)
            rev = {v: k for k, v in p.items()}
            emitted = set()
            for jj, n in enumerate(news):
                if jj in rev:
                    ii = rev[jj]
                    for k in range(len(olds)):
                        if k < ii and k not in p and k not in emitted:
                            items.append(("del", olds[k])); emitted.add(k)
                    items.append(("mod", (olds[ii], n))); emitted.add(ii)
                else:
                    items.append(("ins", n))
            for k in range(len(olds)):
                if k not in p and k not in emitted:
                    items.append(("del", olds[k]))
    return items


def item_html(kind, payload):
    if kind == "same":
        return block_html(payload)
    if kind == "del":
        return block_html(payload, "del")
    if kind == "ins":
        return block_html(payload, "ins")
    old, new = payload
    # A block rewritten this far interleaves into word salad; the old and the
    # new read as themselves instead.
    if ratio(old, new) < 0.55:
        return block_html(old, "del") + block_html(new, "ins")
    if is_fence(old) or is_fence(new):
        return fence_diff(old, new)
    return style_html(render(merge(old, new)))


# ------------------------------------------------------------------- page

def header(path, base, head, changed, total, changes_only):
    scope = "changed blocks only" if changes_only else "whole page"
    return (
        f'<p style="font-family:Georgia,serif;font-size:15pt;color:#000;margin:0 0 4pt 0">'
        f'<b>{html.escape(path)}</b></p>'
        f'<p style="font-family:Georgia,serif;font-size:10pt;color:#444;margin:0 0 2pt 0">'
        f'{html.escape(base)} &#8594; {html.escape(head)} &#183; '
        f'{datetime.date.today().isoformat()} &#183; {changed} of {total} '
        f'paragraphs changed &#183; {scope}</p>'
        f'<p style="font-family:Georgia,serif;font-size:10pt;margin:0 0 14pt 0">'
        f'<span style="{DEL_S}">red struck = removed since {html.escape(base)}</span>'
        f'<span style="color:#444"> &#183; </span>'
        f'<span style="{INS_S}">green = added</span>'
        f'<span style="color:#444"> &#183; plain black = unchanged</span></p>'
        f'<hr>'
    )


def page_diff(base, head, path, changes_only=False):
    fm_o, old = prep(git_show(base, path))
    fm_n, new = prep(git_show(head, path))
    ob, nb = split_blocks(old), split_blocks(new)
    items = build(ob, nb)
    changed = sum(1 for k, _ in items if k != "same")

    if changes_only:
        keep, n = set(), len(items)
        for i, (k, _) in enumerate(items):
            if k == "same":
                continue
            keep.add(i)
            if i > 0:
                keep.add(i - 1)
            if i + 1 < n:
                keep.add(i + 1)
            for j in range(i - 1, -1, -1):          # nearest heading above
                k2, p2 = items[j]
                md = p2[1] if k2 == "mod" else p2
                if is_heading(md):
                    keep.add(j); break
        sel, prev = [], -1
        for i in sorted(keep):
            if prev >= 0 and i > prev + 1:
                sel.append(("gap", None))
            sel.append(items[i]); prev = i
        items = sel

    parts = [header(path, base, head, changed, len(nb), changes_only)]
    if fm_o.strip() != fm_n.strip():
        parts.append(f'<p style="{TAG_STYLE["p"]};color:#666;font-size:9.5pt">'
                     f'<b>front matter</b></p>')
        parts.append(style_html(render(merge(fm_o, fm_n))))
    for kind, payload in items:
        if kind == "gap":
            parts.append('<p style="color:#999;font-size:10pt;margin:6pt 0">'
                         '· · ·</p>')
        else:
            parts.append(item_html(kind, payload))
    if changed == 0:
        parts.append(f'<p style="{TAG_STYLE["p"]}"><i>(no changes)</i></p>')

    doc = ("<!doctype html><html><head><meta charset='utf-8'>"
           "<style>@page{size:8.5in 11in;margin:0.85in 0.8in 0.9in 0.8in}</style>"
           f'</head><body style="{BODY}">' + "\n".join(parts) + "</body></html>")
    return doc, changed


def git_show(rev, path):
    r = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ""


def main(argv):
    out_dir = OUT_DIR
    changes_only = False
    if "--changes-only" in argv:
        argv.remove("--changes-only"); changes_only = True
    if "--out" in argv:
        i = argv.index("--out"); out_dir = argv[i + 1]; del argv[i:i + 2]
    if len(argv) < 3:
        print(__doc__); return 1
    base, head, pages = argv[0], argv[1], argv[2:]
    os.makedirs(out_dir, exist_ok=True)
    for page in pages:
        doc, n = page_diff(base, head, page, changes_only)
        stem = (os.path.splitext(os.path.basename(page))[0]
                + f"_{base}_to_{head}".replace("/", "-"))
        hp = os.path.join(out_dir, stem + ".html")
        with open(hp, "w") as f:
            f.write(doc)
        pdf = os.path.join(out_dir, stem + ".pdf")
        if os.path.exists(pdf):
            os.remove(pdf)
        subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf",
                        "--outdir", out_dir, hp], capture_output=True, text=True)
        print(f"{page}: {n} changed block(s) -> {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
