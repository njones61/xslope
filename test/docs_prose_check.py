"""The docs read like a colleague wrote them: a scan for the phrases, sentence
shapes, spellings and vocabulary that the documentation does not use.

What it scans
-------------
* every ``docs/**/*.md`` page, outside fenced code, inline code, HTML comments
  and link targets, and outside a page's References section (a cited title is
  quoted as printed);
* every user-facing string literal under ``xslope/`` and ``studio/``: a literal
  with a space in it, docstrings excluded. Log lines, dialog text, preflight
  messages and the strength-reduction closing summary all pass through here.
  ``studio/ai`` is left out: its strings are instructions to a model.

The rules
---------
Each rule is a name, a scope and a pattern. The groups:

``flourish``   the phrases that mark a page as machine-written: "worth noting",
               "delve", "crucially", "a testament to", "the key takeaway"...
``signpost``   a paragraph that announces itself instead of saying its point:
               "Two things follow.", "One more thing to know:", "The last
               lesson is patience.", "and this is why."
``reader``     imagined readers: "a reader who tries it will meet the refusal".
``contrast``   "not a limitation but a feature" in prose; in program strings
               also the comma contrast, "a failure, not a budget effect".
``coined``     phrases coined in a session that mean nothing to a stranger:
               "meshed and marched", "ships beside the workbook", a bare
               "march" for a transient run.
``mechanism``  solver internals in a tutorial: corrector, hold test, budget,
               verdict, out-of-balance, continuum, Gauss point.
``british``    the spellings the docs do not use: centre, metre, colour, grey,
               labelled, behaviour, mobilised, modelled, whilst...

A hit is a failure. The fix is to rewrite the sentence, never to widen a rule
or add an exemption; a rule that fires on a sentence a colleague would write
is narrowed, and the narrowing is the commit. (A rule for "X is what makes Y"
was tried and dropped: on pages Norm had reviewed it fired on sentences he had
let stand.)

Run it alone::

    python test/docs_prose_check.py            # every hit, grouped by rule
    python test/docs_prose_check.py --summary  # counts by rule and by file

or through ``run_tests.py``, where it is the ``docs_prose`` row.
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
PY_DIRS = (ROOT / "xslope", ROOT / "studio")

# Scopes: "docs" = every docs page, "tutorials" = docs/tutorials only,
# "strings" = user-facing Python string literals.
ALL = ("docs", "strings")

I = re.IGNORECASE

RULES: list[tuple[str, tuple[str, ...], re.Pattern]] = [
    # ---- flourish -------------------------------------------------------
    ("flourish", ALL, re.compile(
        r"\bworth (?:(?:a|an|the) )?(?:\w+ing|a look|a pause|a mention)\b", I)),
    ("flourish", ALL, re.compile(r"\bit bears (?:noting|repeating|mention)", I)),
    ("flourish", ALL, re.compile(
        r"\b(?:needless to say|it goes without saying|at the end of the day|"
        r"in a nutshell|the beauty of|a testament to|rest assured|"
        r"game[- ]changer|deep dive|dive into|in essence)\b", I)),
    ("flourish", ALL, re.compile(
        r"\b(?:delve[sd]?|delving|tapestry|seamless(?:ly)?|leverag(?:e|es|ed|ing)|"
        r"showcas(?:e|es|ed|ing)|myriad|plethora|"
        r"empower(?:s|ed|ing)?|holistic(?:ally)?|pivotal|nuanced?|"
        r"crucial(?:ly)?|importantly|notably|arguably|elegant(?:ly)?|"
        r"robustly|paradigm|realm)\b", I)),
    ("flourish", ALL, re.compile(
        r"\b(?:it is|it's) (?:important|worth|useful|helpful) to "
        r"(?:note|remember|mention|keep in mind)\b", I)),
    ("flourish", ALL, re.compile(
        r"\b(?:the key (?:insight|takeaway|point|idea|thing|lesson)|takeaways?|"
        r"upshot|bottom line|the short version|"
        r"the (?:last|first|final|real|main|big|one) lesson|"
        r"lesson (?:here|of this|for practice))\b", I)),
    # ---- signpost -------------------------------------------------------
    # A count noun that announces a list ("Two things follow.") rather than
    # counting something ("Two points, both at y = 20"; "Four notes come back").
    ("signpost", ALL, re.compile(
        r"(?:^|(?<=[.!?:] ))(?:One|Two|Three|Four|Five|Several|A few|A couple of) "
        r"(?:more )?(?:things?|habits|points|notes|rules|lessons|consequences|"
        r"observations|caveats|remarks|takeaways)\b[^.:;]{0,40}?"
        r"\b(?:follow|stand out|come out of|to carry|to note|to know|to watch|"
        r"to take away|to remember|to keep in mind|are worth|for practice)\b")),
    ("signpost", ALL, re.compile(
        r"\b(?:things?|habits|lessons|consequences) "
        r"(?:follow|to note|stand out|to take away|to know|to watch|to carry)\b", I)),
    ("signpost", ALL, re.compile(r"(?:^|(?<=[.!?:] ))One more (?:thing|point|note)\b")),
    ("signpost", ALL, re.compile(r"\b(?:for practice|this is why)\b", I)),
    ("signpost", ALL, re.compile(
        r"(?:^|(?<=[.!?] ))(?:Put|Said|Stated) (?:simply|another way|differently|plainly)\b")),
    ("signpost", ALL, re.compile(r"(?:^|(?<=[.!?] ))In short\b")),
    # ---- reader ---------------------------------------------------------
    ("reader", ALL, re.compile(r"\ba reader who\b", I)),
    ("reader", ALL, re.compile(
        r"\breaders? (?:will|would|might|may|who) "
        r"(?:meet|find|notice|wonder|ask|see|expect|try|reach|land)\b", I)),
    ("reader", ALL, re.compile(r"\bmeets? the (?:refusal|error|message|warning|dialog)\b", I)),
    # ---- contrast -------------------------------------------------------
    ("contrast", ALL, re.compile(
        r"\bnot (?:a|an|the|its|their) \w+ but (?:a|an|the|its|their)\b", I)),
    ("contrast", ("strings",), re.compile(
        r"\w, not (?:(?:a|an|the|its|their|just|one|only) )?\w+")),
    # ---- coined ---------------------------------------------------------
    ("coined", ALL, re.compile(
        r"\b(?:meshed and marched|ships? (?:beside|alongside)|shipped (?:beside|alongside)|"
        r"(?:a|the|each|every|per|re)[- ]?march(?:es)?)\b", I)),
    # ---- mechanism (tutorials only) --------------------------------------
    ("mechanism", ("tutorials",), re.compile(
        r"\b(?:corrector|hold test|out-of-balance|verdicts?|continuum|"
        r"constitutive|Gauss points?|stress points?|line nets?|residual forensics)\b", I)),
    # ---- british --------------------------------------------------------
    ("british", ALL, re.compile(
        r"\b(?:centre[sd]?|centring|metres?|colour(?:s|ed|ing)?|greys?|greyed|"
        r"labell(?:ed|ing)|behaviours?|mobilis(?:e|es|ed|ing|ation|able)|"
        r"analys(?:ed|ing)|normalis(?:e|es|ed|ing|ation)|modell(?:ed|ing)|"
        r"(?:optim|minim|maxim|initial|stabil|linear|discret|ideal|penal|"
        r"recogn|visual|summar|emphas|character|parametr|organ|util)is"
        r"(?:e|es|ed|ing|ation)|favour(?:s|ed|ing|able|ite)?|programme|catalogue|"
        r"defence|licence|artefacts?|aluminium|whilst|amongst)\b", I)),
]

_FENCE = re.compile(r"^\s*(```|~~~)")
_COMMENT = re.compile(r"<!--.*?-->", re.S)
_INLINE_CODE = re.compile(r"`[^`\n]*`")
_LINK_TARGET = re.compile(r"\]\([^)]*\)")
_MARKER = re.compile(r"^\s*(?:>\s?|[-*+]\s+|\d+[.)]\s+|#{1,6}\s+|\|)?")
_REFS = re.compile(r"^(#{1,6})\s*references\b", I)
_HEADING = re.compile(r"^(#{1,6})\s")


def _paragraphs(md: Path):
    """Yield (start_line, text) for each prose paragraph of a page, with code,
    comments, link targets and the References section removed and the lines
    of a paragraph joined so sentence-start patterns see whole sentences."""
    raw = md.read_text(encoding="utf-8")
    # Blank out HTML comments (the test tags) but keep line numbers.
    raw = _COMMENT.sub(lambda m: "\n" * m.group(0).count("\n"), raw)
    lines = raw.split("\n")
    in_fence = False
    refs_level = None
    buf, start = [], None
    for i, line in enumerate(lines, 1):
        if _FENCE.match(line):
            in_fence = not in_fence
            line = ""
        elif in_fence:
            line = ""
        h = _HEADING.match(line)
        if h:
            if refs_level is not None and len(h.group(1)) <= refs_level:
                refs_level = None
            r = _REFS.match(line)
            if r:
                refs_level = len(r.group(1))
        if refs_level is not None and not h:
            line = ""
        line = _INLINE_CODE.sub("", line)
        line = _LINK_TARGET.sub("]", line)
        stripped = _MARKER.sub("", line).strip()
        if stripped:
            if start is None:
                start = i
            buf.append(stripped)
        elif buf:
            yield start, " ".join(buf)
            buf, start = [], None
    if buf:
        yield start, " ".join(buf)


def _strings(py: Path):
    """Yield (line, text) for each user-facing string literal in a module: a
    constant with a space in it, twelve characters or longer, that is not a
    docstring."""
    try:
        tree = ast.parse(py.read_text(encoding="utf-8"))
    except SyntaxError:
        return
    doc_ids = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                doc_ids.add(id(body[0].value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in doc_ids:
            s = node.value
            if " " in s and len(s) >= 12:
                yield node.lineno, s


def _excerpt(text: str, m: re.Match, width: int = 90) -> str:
    a = max(0, m.start() - width // 2)
    b = min(len(text), m.end() + width // 2)
    out = text[a:b].replace("\n", " ")
    return ("…" if a else "") + out + ("…" if b < len(text) else "")


def scan():
    """Return the hits: (rule, scope_label, path, line, match, excerpt)."""
    hits = []
    for md in sorted(DOCS.rglob("*.md")):
        rel = md.relative_to(ROOT)
        scopes = {"docs"} | ({"tutorials"} if "tutorials" in md.parts else set())
        for line, text in _paragraphs(md):
            for rule, scope, pat in RULES:
                if not scopes & set(scope):
                    continue
                for m in pat.finditer(text):
                    hits.append((rule, str(rel), line, m.group(0), _excerpt(text, m)))
    for d in PY_DIRS:
        for py in sorted(d.rglob("*.py")):
            if "__pycache__" in py.parts or "ai" in py.parts:
                continue
            rel = py.relative_to(ROOT)
            for line, text in _strings(py):
                for rule, scope, pat in RULES:
                    if "strings" not in scope:
                        continue
                    for m in pat.finditer(text):
                        hits.append((rule, str(rel), line, m.group(0), _excerpt(text, m)))
    return hits


def run():
    """The battery entry point: a list of failure strings, empty when clean."""
    hits = scan()
    return [f"{p}:{ln}: [{rule}] {mt!r} — {ex}" for rule, p, ln, mt, ex in hits]


def main(argv):
    hits = scan()
    if "--summary" in argv:
        from collections import Counter
        by_rule = Counter(h[0] for h in hits)
        by_file = Counter(h[1] for h in hits)
        print("by rule:")
        for k, v in by_rule.most_common():
            print(f"  {v:5d}  {k}")
        print("by file:")
        for k, v in by_file.most_common():
            print(f"  {v:5d}  {k}")
        print(f"{len(hits)} hits")
        return 1 if hits else 0
    last = None
    for rule, p, ln, mt, ex in sorted(hits):
        if rule != last:
            print(f"\n[{rule}]")
            last = rule
        print(f"  {p}:{ln}: {mt!r}\n      {ex}")
    if hits:
        print(f"\n{len(hits)} hits")
        return 1
    print("docs prose: clean")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
