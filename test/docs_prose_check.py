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
``punchline``  the rhetorical identification in place of a plain sentence:
               "the zones are the story", "that is the point.", "That is the
               whole of what makes...", "This pressure is precisely why...",
               a bold "**What decides this row is...**". Plain definitions
               ("phi is the friction angle", "That gap is the clay blanket")
               do not match.
``slogan``     a bold lead-in or note title with a dash or colon tail that
               editorializes: "**What to observe — the zones are the story:**",
               "The result depends on the mesh — but not on the bracket",
               "**The field approaches a new steady state — slowly, and
               unevenly.**"; and a dash tail that draws a moral, "— a
               reminder that...". Labels such as "**Head 1 — the reservoir.**"
               do not match.
``figurative`` figures of speech standing in for the physics: "holds its head
               up", "a hot island of trapped head", "walks down the slope",
               "earns its keep", "on display".
``intensifier`` drama: "precisely why", "the real problem", "the story",
               "strikingly", "the crux", "a reminder that".
``teaser``     a bold lead-in that withholds its point: "**What makes it
               transient.**", "**Why it matters:**".
``entity``     a tutorial as the subject of a verb: "LEM-1 left", "SEEP-3
               builds" (say "in Tutorial SEEP-3").
``process``    project bookkeeping in public prose: "regression guard",
               "held pending", "recorded as open".
``voice``      the program given a voice or a will: "said out loud", "the
               search says", "the model admits", "the plot tells", "the run
               says so in the Log". The program reports, prints, shows,
               returns, lists. A text artifact may say things ("the warning
               says", "the Log says"), and results may agree ("the methods
               agree", "XSLOPE agrees with Slide"): those do not match.
``heading``    headings that are not plain labels of what the section
               contains: a contrast ("The fix is in the ground, not the
               settings"), a question ("Why not a circle?"), an "X is the Y"
               punchline in any clause not opened by a question word ("the
               zones are the story"; "What is in a report" is a label), a
               program voice ("What the other methods say", "What it
               knows"), a rhetorical noun ("the story", "the lesson").

The rules from ``punchline`` down are docs-only, except ``voice``, which also
covers program strings; program strings keep the first group. ``SELFTEST`` plants a violation for every rule and a near-miss
sentence beside it; ``run()`` fails if a planted violation is missed or a
near miss is flagged, so a rule cannot be loosened into uselessness or
widened onto ordinary sentences without the battery saying so.

A hit is a failure. The fix is to rewrite the sentence, never to widen a rule
or add an exemption; a rule that fires on a sentence a colleague would write
is narrowed, and the narrowing is the commit. (A rule for "X is what makes Y"
was tried and dropped: on pages Norm had reviewed it fired on sentences he had
let stand.)

Run it alone::

    python test/docs_prose_check.py            # every hit, grouped by rule
    python test/docs_prose_check.py --summary  # counts by rule and by file
    python test/docs_prose_check.py docs/lem/samples.md   # named pages only
    python test/docs_prose_check.py --selftest # the planted-violation test

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

# Sentence start inside a joined paragraph.
_SS = r"(?:^|(?<=[.!?:] ))"


class _Bold:
    """A rule over the bold spans of a paragraph (and a note's quoted title):
    a span matches when its text satisfies ``cond``. ``lead_only`` keeps the
    spans that open a sentence. ``finditer`` yields the span matches, so the
    scanner treats it like a compiled pattern."""

    _SPAN = re.compile(r"(?<![\w*])\*\*(?=\S)([^*]+?)(?<=\S)\*\*|^!!! ?\w+ \"([^\"]+)\"")
    _LEAD = re.compile(r"[.!?:]\s$")

    def __init__(self, cond: re.Pattern, lead_only: bool = False):
        self.cond = cond
        self.lead_only = lead_only

    def finditer(self, text: str):
        for m in self._SPAN.finditer(text):
            if self.lead_only and m.start() and not self._LEAD.search(text[:m.start()]):
                continue
            if self.cond.search(m.group(1) or m.group(2)):
                yield m


DOCS_ONLY = ("docs",)

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
        r"recogn|visual|summar|emphas|character|parametr|parameter|organ|util)is"
        r"(?:e|es|ed|ing|ation)|favour(?:s|ed|ing|able|ite)?|programme|catalogue|"
        r"defence|licence|artefacts?|aluminium|whilst|amongst)\b", I)),
    # ---- punchline ------------------------------------------------------
    # "The X is the Y" where Y is a rhetorical noun, not a definition.
    ("punchline", DOCS_ONLY, re.compile(
        r"\b(?:is|are|was|were) (?:the|its|their) (?:real |entire |)"
        r"(?:story|crux|moral|catch|twist|punchline|lesson|heart of the matter)\b", I)),
    ("punchline", DOCS_ONLY, re.compile(
        r"\b(?:is|are|was|were) the whole (?:of|reason|question|answer|story|point|"
        r"difference|problem|mechanism|game|trick|[\d.]+\b)", I)),
    ("punchline", DOCS_ONLY, re.compile(
        r"\b(?:is|are|was) (?:the|its|their) (?:real |entire |)point\s*[.;:!)]", I)),
    ("punchline", DOCS_ONLY, re.compile(
        _SS + r"(?:The|That|This|Here) (?:point|lesson|story|crux|catch|trick|moral|upshot|takeaway)"
        r"(?: here)? (?:is|was) (?:that|the|to|not|this|simple|clear)\b")),
    ("punchline", DOCS_ONLY, re.compile(
        _SS + r"(?:That|This) is (?:the|what|where|why|how) (?:lesson|point|story|statement|"
        r"crux|catch|trick|moral)\b")),
    ("punchline", DOCS_ONLY, re.compile(
        _SS + r"(?:That|This|These|Those) (?!is\b|are\b|was\b)[\w'-]+(?: [\w'-]+){0,2} "
        r"(?:is|are|was) (?:(?:precisely|exactly|really|simply) (?:why|the reason)|"
        r"the source of everything)\b")),
    ("punchline", DOCS_ONLY, _Bold(re.compile(
        r"^What (?!to\b)(?:[\w'-]+ ){1,5}(?:is|are) (?:the|a|an)\b"), lead_only=True)),
    # ---- slogan ---------------------------------------------------------
    # A bold lead-in or note title whose dash/colon tail editorializes: a
    # "but/not/never" turn, an "is the" verdict, a comma contrast, or a short
    # verbless flourish after a full clause. "**Head 1 — the reservoir.**" is
    # a label (fewer than four words before the dash) and does not match.
    ("slogan", DOCS_ONLY, _Bold(re.compile(
        r" [—–] (?:but|not|and not|yet|never)\b|"
        r"[—–:] [^—–:]*\b(?:is|are|was) (?:the|what|where|why)\b|"
        r" [—–] [^—–]*\w, not \w|"
        r"^(?:[^\s—–]+ ){4,}[—–] (?:(?:and|but|yet|or|[a-z]+ly),? ?){1,4}[.!]?$"))),
    ("slogan", DOCS_ONLY, re.compile(
        r"[—–] (?:a|an|the) (?:(?:useful|stark|sobering|good|clear|simple|classic|textbook|"
        r"telltale) )?(?:reminder|hallmark|signature|lesson|moral|testament|recipe|"
        r"essence)\b", I)),
    # ---- figurative -----------------------------------------------------
    ("figurative", DOCS_ONLY, re.compile(
        r"\b(?:holds? (?:its|their) (?:head|heads|breath)|keeps? (?:its|their) head|"
        r"earns? (?:its|their) keep|speaks? for itself|front and center|on (?:full )?display|"
        r"wins? out|comes? alive|steals? the|hot island|islands? of (?:trapped|high|low)|"
        r"pockets? of trapped|smoking gun|sweet spot|silver bullet|double-edged|"
        r"walks? (?:down|up|along|across|back)|"
        r"(?:on|by) (?:its|each zone's|their) own [\w/]+ clock)\b", I)),
    # ---- intensifier ----------------------------------------------------
    ("intensifier", DOCS_ONLY, re.compile(
        r"\b(?:precisely (?:why|what|where|how|because)|exactly why)\b", I)),
    ("intensifier", DOCS_ONLY, re.compile(
        r"\bthe real (?:story|question|problem|issue|lesson|answer|reason|culprit|danger|"
        r"test|work|point|cause|difference|risk)\b", I)),
    ("intensifier", DOCS_ONLY, re.compile(
        r"\b(?:the|a|whole|same|different|full) story\b|\btells? (?:the|a) story\b", I)),
    ("intensifier", DOCS_ONLY, re.compile(
        r"\b(?:strikingly|remarkabl[ey]|tellingly|starkly|"
        r"staggering(?:ly)?|unmistakabl[ey]|the crux|a (?:\w+ )?reminder that)\b", I)),
    # ---- teaser ---------------------------------------------------------
    ("teaser", DOCS_ONLY, _Bold(re.compile(
        r"^(?:What|Why) (?:makes|made|it matters|this matters|this means|that means|"
        r"it buys|changes)\b"), lead_only=True)),
    # ---- entity ---------------------------------------------------------
    ("entity", DOCS_ONLY, re.compile(
        r"(?<![Tt]utorial )(?<![Tt]utorials )\b(?:LEM|SEEP|FEM|COMBO|W)-\d+(?:'s)? "
        r"(?:builds?|built|adds?|added|shows?|showed|teaches|taught|runs?|ran|uses?|used|"
        r"covers?|takes?|took|leaves?|left|sets?|walks?|introduces?|introduced|drew|draws?)\b")),
    # ---- process --------------------------------------------------------
    ("process", DOCS_ONLY, re.compile(
        r"\b(?:regression guard|held pending|recorded as open)\b", I)),
    # ---- voice ----------------------------------------------------------
    # The program as a speaker or a mind. The subject is an agent that has
    # no voice (search, model, solver, plot, run...), never a text artifact
    # (warning, Log, message, output, answer), which may "say" what it says;
    # and "agrees" is caught only with such a subject, so numeric agreement
    # ("the two methods agree", "XSLOPE agrees with Slide") does not match.
    ("voice", ALL, re.compile(r"\bout loud\b", I)),
    ("voice", ALL, re.compile(
        r"\b(?:the|this|that|its|each|every|a) (?:[\w'-]+ )?"
        r"(?:search(?:es)?|model|solver|program|engine|analysis|software|plot|figure|run|"
        r"XSLOPE|Studio) (?:says|said|is saying|tells|told|admits|admitted|agrees|agreed|"
        r"knows|wants|insists|believes|thinks|complains)\b", I)),
]


class _HeadingPunch:
    """An "X is the Y" pronouncement in a heading clause (split at dashes,
    colons, semicolons and commas) that a question word does not open:
    "The fix is in the ground", "the zones are the story" match; "What is in
    a report" is a label and does not."""

    _CLAUSE = re.compile(r"[^—–:;,]+")
    _IS = re.compile(r"\b(?:is|are|was|were) (?:the|in|what|where|why|how)\b", I)
    _WH = re.compile(r"^\s*(?:What|How|Which|When|Where|Why|Who)\b", I)

    def finditer(self, text: str):
        for c in self._CLAUSE.finditer(text):
            if self._WH.match(c.group(0)):
                continue
            for m in self._IS.finditer(c.group(0)):
                yield m


# Heading rules run on each heading's text alone (its {#anchor} and inline
# code removed); docs only. A heading is a plain label of what the section
# contains.
HEADING_RULES: list[tuple[str, tuple[str, ...], object]] = [
    ("heading", DOCS_ONLY, re.compile(r",\s*not\b|\s[—–]\s*(?:but\s+)?not\b", I)),
    ("heading", DOCS_ONLY, re.compile(r"\?\s*$")),
    ("heading", DOCS_ONLY, _HeadingPunch()),
    ("heading", DOCS_ONLY, re.compile(
        r"\b(?:it|search(?:es)?|model|solver|program|engine|analysis|plot|figure|methods?|"
        r"XSLOPE|Studio|assistant) (?:says?|said|tells?|knows?|wants?|admits?|agrees?|"
        r"thinks?|believes?)\b", I)),
    ("heading", DOCS_ONLY, re.compile(
        r"\bthe (?:story|lesson|moral|catch|trick|secret|punchline|crux)\b", I)),
]

# Each rule's planted violation and the nearest ordinary sentence it must
# leave alone: (rule, text, must_hit).
SELFTEST: list[tuple[str, str, bool]] = [
    ("punchline", "**What to observe:** the zones are the story.", True),
    ("punchline", "It is small enough to check by hand, and that is the point.", True),
    ("punchline", "That is the whole of what makes a boundary time-varying.", True),
    ("punchline", "That cascade is the whole 0.031.", True),
    ("punchline", "The point is the answer that looks fine.", True),
    ("punchline", "That is the statement the table makes.", True),
    ("punchline", "This retained core pressure is precisely why drawdown is dangerous.", True),
    ("punchline", "**What decides this row is the rock bridges.** The bridges fail first.", True),
    ("punchline", "That opposition is the source of everything we measure here.", True),
    ("punchline", "That tension is why the two searches disagree.", False),
    ("punchline", "That amplification is why the yield acceleration is the compared quantity.", False),
    ("punchline", "Here phi is the friction angle and c is the cohesion.", False),
    ("punchline", "That gap is the clay blanket. That last edge is the dipping bedrock.", False),
    ("punchline", "The point is selected automatically after the coarse sweep.", False),
    ("punchline", "The seepage solution covers the whole model. This is the main sheet.", False),
    ("punchline", "Each line's capacity is the whole mat, and the core is the whole material.", False),
    ("slogan", "**What to observe — the zones are the story:**", True),
    ('slogan', '!!! note "The result depends on the mesh — but not on the bracket"', True),
    ("slogan", "**The field approaches a new steady state — slowly, and unevenly.** By day 400 it is.", True),
    ("slogan", "**Van Genuchten discharge vs SEEP2D — a reporting difference, not a solver difference**", True),
    ("slogan", "The probability is 11.6% — a reminder that a high factor of safety is not enough.", True),
    ("slogan", "**Head 1 — the reservoir.** Select the upstream face.", False),
    ("slogan", "**Profile Line 1 — material 1:** enter the points below.", False),
    ("slogan", "**Stage 1 — full pool, drained.** The first stage uses the drained strengths.", False),
    ("slogan", "**Surficial skin slides are filtered — in grid mode only.** On a steep face...", False),
    ("slogan", "Read every warning — a warning passed over silently is a defect.", False),
    ("figurative", "The low-permeability core holds its head up.", True),
    ("figurative", "The core is a hot island of trapped total head.", True),
    ("figurative", "The exit point walks down the slope.", True),
    ("figurative", "The pile holds its own line, and the search walks the entry along the ground.", False),
    ("figurative", "The core keeps a high head long after the shell has drained.", False),
    ("intensifier", "This is precisely why rapid drawdown is dangerous.", True),
    ("intensifier", "The real problem is the core.", True),
    ("intensifier", "The story is in the colors on the layers.", True),
    ("intensifier", "The two results differ strikingly.", True),
    ("intensifier", "The factor of safety drops dramatically.", False),
    ("intensifier", "The total is reproduced exactly. Studio shows exactly what the log reports.", False),
    ("teaser", "**What makes it transient.** Two things are added.", True),
    ("teaser", "**Why it matters:** the core drains slowly.", True),
    ("teaser", "**What to observe:** the core drains slowly. **How it is set.** On the mat sheet.", False),
    ("entity", "The warnings LEM-1 left behind are cleared here.", True),
    ("entity", "Tutorial SEEP-3 builds the dam. In Tutorial LEM-1 the slope is homogeneous.", False),
    ("process", "This file is included as a regression guard on the support mechanics.", True),
    ("process", "The model is included to check the support mechanics on a mirrored slope.", False),
    ("flourish", "Crucially, the core drains slowly. Notably, the shell does not.", True),
    ("signpost", "Two things follow. The core drains slowly.", True),
    ("voice", "The search says this out loud. Its run output includes the line:", True),
    ("voice", "Being uncertain out loud is not a failure.", True),
    ("voice", "Fifty-six trial circles admit no solution, and the model admits it.", True),
    ("voice", "The run says so in the Log rather than counting it as a failure.", True),
    ("voice", "If it never crosses, the plot says which way to widen the range.", True),
    ("voice", "XSLOPE's circular search on the weak-layer model agrees.", True),
    ("voice", "The plot tells the whole history of the search.", True),
    ("voice", "The run output reports this. The search prints the circle when it converges.", False),
    ("voice", "The warning says the table has no tensile cutoff, and the Log says the same.", False),
    ("voice", "Spencer and Morgenstern-Price agree to 0.1%, and XSLOPE agrees with Slide.", False),
    ("voice", "With no tension in the model, the methods agree.", False),
    ("voice", "Nothing in the output says so, and the answer says which of the two it computed.", False),
    ("voice", "A system that admits no admissible root is reported as such.", False),
    ("heading", "The fix is in the ground, not the settings", True),
    ("heading", "Why not a circle?", True),
    ("heading", "How deep does the crack really need to be?", True),
    ("heading", "What to observe — the zones are the story", True),
    ("heading", "What the other methods say", True),
    ("heading", "What it knows before it runs anything", True),
    ("heading", "The lesson of the second search", True),
    ("heading", "Factor of safety by method", False),
    ("heading", "Adding a tension crack", False),
    ("heading", "What is in a report", False),
    ("heading", "Part 3 — When a sheet is a slip surface, and when it is bonded", False),
    ("heading", "London clay, where the two fits agree", False),
    ("heading", "Choose how you want to build it", False),
    ("heading", "How the search finds the answer", False),
    ("heading", "A flux is a rate normal to the boundary", False),
]


def selftest() -> list[str]:
    """Run every SELFTEST case through its rule group; return the failures."""
    fails = []
    for rule, text, must in SELFTEST:
        hit = any(True for name, scope, pat in RULES + HEADING_RULES if name == rule
                  for _ in pat.finditer(text))
        if hit != must:
            fails.append(f"selftest [{rule}] {'missed' if must else 'false hit'}: {text!r}")
    return fails

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


_HEADING_TEXT = re.compile(r"^#{1,6}\s+(.*?)\s*(?:\{[^}]*\})?\s*#*\s*$")


def _headings(md: Path):
    """Yield (line, text) for each Markdown heading of a page outside fenced
    code and HTML comments, with its {#anchor}/attribute list and inline code
    removed."""
    raw = md.read_text(encoding="utf-8")
    raw = _COMMENT.sub(lambda m: "\n" * m.group(0).count("\n"), raw)
    in_fence = False
    for i, line in enumerate(raw.split("\n"), 1):
        if _FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = _HEADING_TEXT.match(line)
        if m:
            yield i, _INLINE_CODE.sub("", m.group(1)).strip()


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


def scan_page(md: Path):
    """The hits on one docs page: (rule, path, line, match, excerpt)."""
    md = Path(md).resolve()
    try:
        rel = str(md.relative_to(ROOT))
    except ValueError:
        rel = str(md)
    scopes = {"docs"} | ({"tutorials"} if "tutorials" in md.parts else set())
    hits = []
    for line, text in _paragraphs(md):
        for rule, scope, pat in RULES:
            if not scopes & set(scope):
                continue
            for m in pat.finditer(text):
                hits.append((rule, rel, line, m.group(0), _excerpt(text, m)))
    for line, text in _headings(md):
        for rule, scope, pat in HEADING_RULES:
            if not scopes & set(scope):
                continue
            for m in pat.finditer(text):
                hits.append((rule, rel, line, m.group(0), text))
    return hits


def scan(pages=None):
    """Return the hits: (rule, path, line, match, excerpt). With ``pages``,
    only those docs pages are scanned and the program strings are skipped."""
    if pages:
        return [h for p in pages for h in scan_page(p)]
    hits = []
    for md in sorted(DOCS.rglob("*.md")):
        hits.extend(scan_page(md))
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
    return selftest() + [f"{p}:{ln}: [{rule}] {mt!r} — {ex}" for rule, p, ln, mt, ex in hits]


def main(argv):
    if "--selftest" in argv:
        fails = selftest()
        for f in fails:
            print(f)
        print(f"selftest: {len(SELFTEST) - len(fails)}/{len(SELFTEST)} cases pass")
        return 1 if fails else 0
    pages = [a for a in argv if not a.startswith("--")]
    hits = scan(pages or None)
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
