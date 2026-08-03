"""Voice check: the verification pages are documentation, not status reports.

These pages are read by a stranger who wants to know whether XSLOPE reproduces a
published problem and what the caveats are.  They are not addressed to the
maintainer, and they are not a record of how the work went.  Prose that says
"the obvious suspect was...", "their factors are withdrawn", "four candidate
causes have been measured" is written to whoever commissioned the work; a reader
who was not there gets nothing from it.  The replacement voice states what IS —
the comparison, the number, the physical reason, the caveat a user needs.

This is a lint, so it is deliberately narrow: it flags phrasings that are
campaign voice in essentially any context, not everything that could be written
better.  Judgement cases stay with the writer.  Each hit names the page, the
line, the phrase and the line's text, so the fix is never a guessing game.

Prose only.  Test tags and other HTML comments, fenced code blocks, inline code
spans, link targets and image paths are removed before matching, so a material
named `us` in a code span or a URL containing "our" cannot fire.

Genuine usage prose sometimes needs a banned word — a caveat addressed to the
reader ("we recommend" in a how-to, "your model") — so each page config carries
``voice_allow``: a list of (phrase, distinctive substring of the line) pairs,
each of which must fire, exactly like the other exemption lists.  A dead
allowance fails, which is what stops the list silently accumulating.
"""
import re
import sys

#: Phrases that read as campaign voice wherever they appear.  Grouped by what is
#: wrong with them, because the fix differs: a first-person sentence needs a
#: subject change, a process sentence usually needs deleting.
BANNED = [
    # -- first person: these pages have no "we" -----------------------------
    (r"\bwe\b", "first person"),
    (r"\bwe['’](?:re|ve|ll|d)\b", "first person"),
    (r"\bour\b", "first person"),
    (r"\bours\b", "first person"),
    (r"\bus\b", "first person"),
    (r"\bI (?:have|had|will|would|ran|measured|checked|found|built|note|think|see)\b",
     "first person"),
    (r"\blet['’]s\b", "first person"),

    # -- project process rather than physics ---------------------------------
    (r"\bwithdraw(?:n|s|ing)?\b", "project process"),
    (r"\bwithdrew\b", "project process"),
    (r"\bre-?lock(?:ed|s|ing)?\b", "project process"),
    (r"\bheld pending\b", "project process"),
    (r"\bpending (?:a|the) (?:decision|review|run|rebuild)\b", "project process"),
    (r"\bfresh solve\b", "project process"),
    (r"\b(?:has|have|had) been (?:measured|probed|re-?run|re-?checked)\b",
     "project process"),
    (r"\b(?:was|were) probed\b", "project process"),
    (r"\brecorded as measured\b", "project process"),
    (r"\bit was built\b", "project process"),
    (r"\bthis (?:session|round|pass|sweep run)\b", "project process"),
    (r"\bsub-?agents?\b", "project process"),
    (r"\bTODO\b", "project process"),
    (r"\bFIXME\b", "project process"),

    # -- investigation narrative addressed at whoever ordered it -------------
    (r"\bthe obvious (?:suspect|candidate)\b", "investigation narrative"),
    (r"\bthe first place to look\b", "investigation narrative"),
    (r"\bexonerat(?:e|ed|es|ing)\b", "investigation narrative"),
    (r"\bcandidate causes?\b", "investigation narrative"),
    (r"\bworth checking\b", "investigation narrative"),
    (r"\bthe question was never\b", "investigation narrative"),
    (r"\bforcing agreement would\b", "investigation narrative"),
    (r"\btun(?:ing|ed) [^.]{0,40}? to the answer\b", "investigation narrative"),

    # -- time measured from the project, not from the problem ----------------
    (r"\bthe earlier (?:reading|assumption|value|offset|number|figure)\b",
     "project-relative time"),
    (r"\bcorrecting the earlier\b", "project-relative time"),
    (r"\bused to (?:fall|be|sit|read|give)\b", "project-relative time"),
    (r"\bno longer (?:do|does)\b", "project-relative time"),
    (r"\bit now sits\b", "project-relative time"),
    # NOT banned: a bare "now sits"/"now reads".  Both are ordinary descriptive
    # English for a variant of the model under discussion ("the firm base now
    # sits at depth D = 1.5H", "with the cap removed the file now reads ..."),
    # and banning them cost more in false positives than the campaign-voice
    # cases they caught, which "it now sits" and "used to ..." already reach.
]

FENCE = re.compile(r"^\s*(```|~~~)")
COMMENT = re.compile(r"<!--.*?-->", re.S)
CODESPAN = re.compile(r"`[^`]*`")
LINKTARGET = re.compile(r"\]\([^)]*\)")
IMGPATH = re.compile(r"!\[[^\]]*\]\([^)]*\)")
URL = re.compile(r"https?://\S+")


def _prose_lines(path):
    """(line number, prose) for every line that carries prose.

    Fenced blocks are dropped whole; HTML comments (which is what a test tag is)
    are dropped even when they span lines; code spans, link targets and bare
    URLs are blanked in place so their offsets do not shift and a banned word
    inside a file name or a citation URL cannot fire.
    """
    raw = open(path, encoding="utf-8").read()
    # Blank multi-line HTML comments in place, keeping the newlines so the line
    # numbers a reader is given still point at the page.
    raw = COMMENT.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), raw)
    out, in_fence = [], False
    for i, line in enumerate(raw.split("\n"), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        text = line
        # Keep image ALT text (a reader sees it); drop only the path.
        text = IMGPATH.sub(lambda m: m.group(0).split("](")[0] + "]", text)
        text = LINKTARGET.sub("]", text)
        text = URL.sub(" ", text)
        text = CODESPAN.sub(lambda m: " " * len(m.group(0)), text)
        out.append((i, text))
    return out


def scan(path, cfg=None):
    """Every banned phrase on one page: (line, phrase, kind, text)."""
    allow = list(getattr(cfg, "voice_allow", []) or [])
    fired = set()
    hits = []
    for lineno, text in _prose_lines(path):
        for pattern, kind in BANNED:
            for m in re.finditer(pattern, text, re.I):
                phrase = m.group(0)
                ex = next(((p, s) for p, s in allow
                           if p.lower() == phrase.lower() and s in text), None)
                if ex:
                    fired.add(ex)
                    continue
                hits.append((lineno, phrase, kind, text.strip()))
    dead = [a for a in allow if a not in fired]
    return hits, dead


def run(path, cfg, report=print):
    """Check one page.  Returns the failure count."""
    hits, dead = scan(path, cfg)
    name = path.split("/")[-1]
    if not hits and not dead:
        report(f"  voice     : clean")
        return 0
    for lineno, phrase, kind, text in hits:
        report(f"  voice     : {name}:{lineno} {kind} — \"{phrase}\"")
        report(f"              {text[:150]}")
    for phrase, sub in dead:
        report(f"  voice     : dead allowance {phrase!r} / {sub!r} — "
               f"no line matches it any more")
    return len(hits) + len(dead)


def main(argv):
    from .pages import PAGES
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    pagedir = os.path.join(os.path.dirname(os.path.dirname(here)),
                           "docs", "verification")
    names = argv[1:] or sorted(PAGES)
    total = 0
    for n in names:
        p = n if n.endswith(".md") else os.path.join(pagedir, n + ".md")
        key = os.path.basename(p)[:-3]
        print(f"{key}:")
        total += run(p, PAGES.get(key), report=print)
    print(f"\n{total} voice problem(s)")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
