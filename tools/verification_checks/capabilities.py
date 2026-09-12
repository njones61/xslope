#!/usr/bin/env python3
"""Capability negations: a page may not say XSLOPE lacks something without a check.

A verification row that cannot be built has to say why, and the cheapest thing to
write is that XSLOPE has no input for it.  Written from memory that sentence is
wrong about as often as it is right: the ``lloads`` sheet (x, y, P, Angle) has
carried a concentrated line load since v14 and ``fem.py`` applies it, and two
public texts still said the loads sheet could not carry one.  A reader takes
such a sentence as the capability inventory it looks like, and a blocked row
built on it never gets revisited.

So every sentence on these pages that NEGATES a capability — *does not carry*,
*cannot model*, *has no counterpart*, *not supported*, *not implemented* — must
name a capability that is on ``ABSENT`` below, the one list of things XSLOPE
genuinely does not have.  Each entry there carries the grep that establishes it.
A sentence that reaches no entry is the failure this check exists for: either
the capability is real and the entry (with its evidence) is missing, or the page
is wrong about XSLOPE and the prose is what has to change.  An entry that
matches nothing is reported dead, like every other exemption list here, so the
list cannot quietly accumulate claims the pages stopped making.

What fires
----------
A negation fires only when it is ABOUT XSLOPE.  Three conditions, all checked
against the sentence the negation sits in:

* the sentence names XSLOPE or one of its artifacts — the package, a sheet of
  the input template, the loader, the importer, an engine (LEM/FEM/SSRM/seepage),
  or a named part of the model it builds;
* the negation's subject is not the SOURCE.  "a quantity the example does not
  carry", "storage the vendor's model does not have" and "the manual prints no
  radius" are statements about the published problem, and say nothing about what
  XSLOPE can do;
* the negation's object is not the verification work itself.  A value that
  "cannot carry the comparison", a surface that "cannot be locked" and a row
  with "no published counterpart" are about the scoring, not about a capability.

Prose and table cells alike: a blocked row's reason lives in a cell, which is
exactly where these sentences are written.

Usage: python -m tools.verification_checks.capabilities [page ...]
"""
import os
import re
import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class Absent:
    """One capability XSLOPE genuinely does not have.

    ``capability`` names it in the terms a reader would search for; ``page`` and
    ``marker`` say which page's rows cite it (``marker`` is a distinctive
    substring of the sentence, and may match several rows that cite the same
    absence); ``evidence`` records where it was looked for and not found, so the
    next person to doubt the entry re-runs the grep instead of the reasoning.
    """
    capability: str
    page: str
    marker: str
    evidence: str


#: The list.  One entry per capability per page that cites it.  ADD AN ENTRY ONLY
#: AFTER GREPPING: xslope/fileio.py, xslope/solve.py, xslope/fem.py,
#: xslope/seep.py, xslope/joints.py, xslope/sensitivity.py, the template's sheet
#: list and docs/usage/input_template.md decide the question, not recollection.
ABSENT = [
    # ---------------------------------------------------------------- strength
    Absent("orientation-dependent (anisotropic) shear strength",
           "geostudio", "no such term",
           "no anisotropic/orientation/dip-dependent strength in fileio.py or "
           "solve.py; the only 'anisotropic' in the docs is hydraulic "
           "conductivity k1/k2 (docs/usage/input_template.md:382)"),
    Absent("orientation-dependent (anisotropic) shear strength",
           "geostudio", "no XSLOPE counterpart",
           "as above; the .gsz importer flags the anisotropic-function material "
           "rather than dropping the orientation term"),
    Absent("orientation-dependent (anisotropic) shear strength",
           "geostudio", "Compound strength",
           "as above; SLOPE/W's compound-strength model is the same "
           "orientation-dependent idea in its second formulation"),
    Absent("orientation-dependent (anisotropic) shear strength",
           "rocscience", "no orientation-dependent strength model",
           "as above"),

    # --------------------------------------------------------------- standards
    Absent("design-standard partial factors (Eurocode 7 and the like)",
           "geostudio", "no native",
           "slide2.py:691 declines a Slide2 model that selects a design standard "
           "— 'xslope applies no partial factors'; no partial-factor input on any "
           "sheet of docs/usage/input_template.md"),

    # ----------------------------------------------------------------- methods
    Absent("an interslice force function the user supplies (Slide2's GLE)",
           "rocscience", "no XSLOPE counterpart",
           "solve.py mprice() takes f_type 'constant' or 'half_sine' only "
           "(_mp_f_vals); f_type is not an input — no hit in fileio.py or "
           "docs/usage/input_template.md"),
    Absent("the Corps of Engineers 2-stage drawdown procedure, and staged "
           "Lowe & Karafiath",
           "rocscience", "not implemented",
           "advanced.py rapid_drawdown() is the three-stage Duncan-Wright-Wong "
           "procedure and the only staged strength path; no 2-stage variant in "
           "solve.py or advanced.py"),
    Absent("the Corps of Engineers 2-stage drawdown procedure",
           "rocscience", "*not supported* — Corps 2-stage",
           "as above"),
    Absent("permanent-displacement (Newmark) analysis from a yield acceleration",
           "rocscience", "reads no XSLOPE model",
           "no newmark / yield-acceleration / displacement-analysis hit anywhere "
           "in xslope/*.py; the seismic input is a coefficient, not a record"),

    # ------------------------------------------------------------------- water
    Absent("an exponential unsaturated conductivity law",
           "rocscience_groundwater", "does not implement",
           "seep.py offers three kr laws: linear frontal (kr_frontal), van "
           "Genuchten (KR_VG) and Gardner (KR_GARD), and its Gardner is the "
           "POWER form kr = 1/(1 + a*psi^n) (kr_gardner_vec, seep.py:1553-1556), "
           "not k = ks*exp(-alpha*psi)"),
    Absent("RS2's built-in 'Simple' conductivity and water-content functions",
           "rocscience", "which XSLOPE does not implement",
           "seep.py's material laws are linear frontal / van Genuchten / Gardner "
           "(KR_LF, KR_VG, KR_GARD); no vendor-preset function table"),
    Absent("RS2's built-in 'Simple' conductivity and water-content functions",
           "rocscience_groundwater", "RS2's built-in Simple curve",
           "as above; the vendor file does not store the curve's parameters "
           "either, so a linear front stands in for it"),

    # -------------------------------------------------------------- FEM inputs
    Absent("a tension crack in the finite-element model",
           "rs2", "no FEM representation",
           "the main sheet's tcrack depth/water (fileio.py:1303) is read by LEM "
           "slice generation only; no tension-crack input reaches fem.py, and "
           "docs/usage/input_template.md:256 states t_cut is the FEM's tensile "
           "control and the tension-crack parameters the LEM's"),
    Absent("staged construction: a stage that excavates part of the section",
           "rs2_joints", "no staged construction",
           "no 'excavat' anywhere in fem.py, fileio.py, mesh.py, solve.py or "
           "joints.py; fem.py's stage_list is a LOAD stage list (one entry, "
           "fem.py:5988), never a change of geometry or a removal of elements. "
           "The opening ITSELF is buildable — the domain is the union of the "
           "material zones (preflight.py:2583), so a region no zone covers is "
           "simply not meshed"),

    # ------------------------------------------------------------------ joints
    Absent("a hyperbolic displacement- and work-softening joint law",
           "rs2_joints", "interface element does not have",
           "no hyperbolic / softening term in joints.py or joint.py; the "
           "interface carries peak and residual Mohr-Coulomb strengths with "
           "dilation, and nothing between them"),

    # ---------------------------------------------------------- reinforcement
    Absent("block shear and facing flexure in a segmental wall",
           "published", "no element for block shear or facing flexure",
           "reinforcement is a line with Tmax, a bond law and end anchorage "
           "(fileio.py reinforce_available_tension); no facing element in fem.py "
           "and no facing input on the reinforce sheet. The face CONNECTION is "
           "not in this class: Tend1/Tend2 are the connection/plate capacity "
           "(docs/usage/input_template.md:863)"),
    Absent("a pullout resistance factor that varies along one line",
           "published", "no input to carry that variation",
           "Delta is one value per reinforcement line (fileio.py "
           "reinforce_pullout_profile), so F* may differ between lines and not "
           "along one"),
    Absent("surface loads in the reinforcement pullout overburden",
           "published", "overburden law does not read",
           "reinforce_effective_overburden (fileio.py:728) is the soil column's "
           "own weight less the declared pore pressure; the dloads sheet's "
           "pressures do not enter it"),

    # ------------------------------------------------------------ reliability
    Absent("a pore-pressure ratio ru as an uncertain input",
           "rocscience", "does not perturb ru",
           "the mat sheet's standard-deviation columns are s(g), s(c), s(f), "
           "s(c/p), s(d), s(psi) (fileio.py:1871-1876, 2961-2962); ru carries "
           "no sigma, so no driver can sample it"),
    Absent("a spatial-averaging (autocorrelation) length for material properties",
           "geostudio", "no autocorrelation-length input",
           "reliability.py samples each material property independently — no "
           "correlation / autocorrelation / spatial-averaging term in "
           "reliability.py or sensitivity.py, and geostudio.py:990 declines a "
           ".gsz that declares one"),
]


# --------------------------------------------------------------------- matching

#: The negations, built from a grep of docs/verification/*.md for every way these
#: pages phrase an absence.  Each is ``(pattern, ours)``, matched
#: case-insensitively; ``ours`` says the phrasing is a claim about XSLOPE
#: wherever it appears on a verification page (a row's *not supported* verdict,
#: an "XSLOPE has no ...", a "has no input"), so it fires without needing the
#: sentence to name XSLOPE again.  The rest need the sentence to name XSLOPE or
#: one of its artifacts — see ``XSLOPE_TOKENS``.
#:
#: NOT here, deliberately: the bare word "unsupported", which on ssrm.md means a
#: slope with no plate in it ("SSRM FS, unsupported"), and "no lock possible",
#: which is a scoring verdict rather than a claim about the software.
NEGATIONS = [
    (r"\b(?:do(?:es)?|did)\s+not\s+(?:carry|carries|support|have|has|model|"
     r"implement|offer|provide|expose|read|apply|include|accept|recognize|"
     r"recognise|know|perturb|sample|vary)\b", False),
    (r"\bcannot\s+(?:carry|model|represent|apply|reproduce|express|take|accept|"
     r"be\s+(?:built|modell?ed|represented|posed|carried|entered))\b", False),
    (r"\bcan(?:'|’)?t\s+(?:carry|model|represent|reproduce|do)\b", False),
    (r"\bis\s+not\s+(?:available|modell?ed)\b", False),
    (r"\bhas\s+no\s+way\b", False),
    (r"\b(?:has|have)\s+no\s+(?:\w+[- ]){0,3}(?:counterpart|equivalent|"
     r"representation|term|model|feature|option|analogue|analog|element|"
     r"parameter|setting|sheet|column|way)\b", False),
    (r"\blacks\b", False),
    # ... and the phrasings that are about XSLOPE by construction
    (r"\bnot\s+supported\b", True),
    (r"\b(?:is\s+)?not\s+implemented\b", True),
    (r"\b(?:has|have)\s+no\s+input\b", True),
    (r"\bno\s+XSLOPE\b", True),
    (r"\bXSLOPE\s+has\s+no\b", True),
    (r"\bnothing\s+in\s+XSLOPE\b(?!\s+is\s+missing)", True),
    (r"\bno\s+native\b", True),
    (r"\bno\s+staged\s+construction\b", True),
    (r"\bno\s+support\s+for\b", True),
]

#: A matched phrase that is not a claim about the software at all.  "no
#: published counterpart" says the SOURCE prints nothing to compare against,
#: which is the commonest thing these rows have to say about a searched value.
NOT_A_CLAIM = re.compile(r"published\s+counterpart", re.I)

#: The status vocabulary the corpus tables use.  A sentence that recites three
#: or more of these words is the legend that defines them, not a verdict on a
#: problem, and names no capability.
STATUS_WORDS = re.compile(
    r"\*{0,2}(?:built|covered|partial|planned|blocked|reported|"
    r"no lock possible|not supported|not implemented)\*{0,2}", re.I)

#: XSLOPE, and the parts of it these pages name.  A negation is a capability
#: claim only when the sentence naming it names one of these.
XSLOPE_TOKENS = re.compile(
    r"\bXSLOPE(?:'|’)?s?\b|"
    r"\b(?:LEM|FEM|SSRM)\b|"
    r"\b(?:lloads|dloads|reinforce|circles|non-?circ|profile|polygon|mat|piles|"
    r"joints?|seep)\s+sheet\b|\bloads\s+sheet\b|\binput\s+template\b|"
    r"\bthe\s+(?:loader|importer|mesher|preflight|template)\b|"
    r"\b(?:interface\s+element|material\s+model|overburden\s+law|"
    r"capacity\s+envelope)\b", re.I)

#: A negation whose SUBJECT is the source says nothing about XSLOPE.  Matched
#: against the text immediately before the negation.
SOURCE_SUBJECT = re.compile(
    r"\b(?:the\s+)?(?:vendor(?:'|’)?s?|manual(?:'|’)?s?|source(?:'|’)?s?|"
    r"example(?:'|’)?s?|published|author(?:'|’)?s?|reference(?:'|’)?s?|"
    r"problem(?:'|’)?s?|paper(?:'|’)?s?|his|their)"
    r"(?:\s+\w+){0,2}\s*$", re.I)

#: The vendors run engines with the same names XSLOPE's carry.  "RS2's own SSRM"
#: and "Slide2's Spencer" are the vendor's, and a token inside one of these does
#: not make the sentence a statement about XSLOPE.
VENDOR_ENGINE = re.compile(
    r"\b(?:RS2|Slide2?|SLOPE/W|SEEP/W|GeoStudio|PLAXIS|UDEC|FLAC|XSTABL)"
    r"(?:'|’)?s?\s+(?:own\s+)?(?:SSRM|SSR|FEM|LEM)\b", re.I)


def _xslope_spans(text):
    """Where the text names XSLOPE or one of its artifacts, the vendors' own
    engines excluded."""
    vendor = [m.span() for m in VENDOR_ENGINE.finditer(text)]
    return [m.span() for m in XSLOPE_TOKENS.finditer(text)
            if not any(a <= m.start() and m.end() <= b for a, b in vendor)]


#: ``it``/``they``/``which`` carries the subject of the clause before it, so a
#: pronoun subject is resolved by looking back through the sentence: where what
#: stands there is a vendor or a source, the negation is about the vendor.
PRONOUN_SUBJECT = re.compile(r"\b(?:it|they|which|that|this)\s*(?:\*{1,2})?\s*$", re.I)
VENDOR_NAME = re.compile(
    r"\b(?:RS2|Slide2?|SLOPE/W|SEEP/W|GeoStudio|PLAXIS|UDEC|FLAC|XSTABL|"
    r"Rocscience|Seequent|SEEP2D|the\s+vendor|the\s+manual|the\s+example|"
    r"the\s+source|the\s+author)", re.I)

#: ... and one whose OBJECT is the verification work — a comparison, a lock, a
#: dot, a published counterpart — is about the scoring, not about a capability.
WORK_OBJECT = re.compile(
    r"^\s*(?:it|them|this|that|these|those)?\s*"
    r"(?:the\s+)?(?:\w+\s+){0,2}(?:comparison|comparisons|lock|locks|dot|dots|"
    r"leg|legs|tally|pairing|pairings|scoring|verdict|column|row|table|section|"
    r"figure)\b|"
    r"^\s*(?:be\s+)?(?:locked|compared|scored|paired|regenerated|tabulated)\b",
    re.I)

FENCE = re.compile(r"^\s*(```|~~~)")
COMMENT = re.compile(r"<!--.*?-->", re.S)
CODESPAN = re.compile(r"`[^`]*`")
LINKTARGET = re.compile(r"\]\([^)]*\)")
IMGPATH = re.compile(r"!\[[^\]]*\]\([^)]*\)")
URL = re.compile(r"https?://\S+")
SENTENCE_END = re.compile(r"(?<=[.!?])\s+")


def _clean(line):
    """One line with the non-prose blanked, exactly as the voice check does it."""
    text = IMGPATH.sub(lambda m: m.group(0).split("](")[0] + "]", line)
    text = LINKTARGET.sub("]", text)
    text = URL.sub(" ", text)
    return CODESPAN.sub(lambda m: " " * len(m.group(0)), text)


def _blocks(path):
    """(first line number, text) for every block of prose on the page.

    A table row is one block: the row states one problem's verdict and the cell
    carrying the reason is the one the check is written for, while the cells
    beside it name the problem the reason is about.  Every other run of
    consecutive non-blank lines is joined into a paragraph, because these pages
    wrap mid-sentence and the subject of a negation is regularly on the line
    above it.  Fenced blocks, HTML comments (test tags), code spans, link
    targets and bare URLs are removed first, so a sheet name in a code span
    cannot by itself make a sentence look like a capability claim.
    """
    raw = open(path, encoding="utf-8").read()
    raw = COMMENT.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), raw)
    blocks, para, start, in_fence = [], [], None, False
    for lineno, line in enumerate(raw.split("\n"), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        text = _clean(line)
        is_row = "|" in text
        if not text.strip() or is_row or text.lstrip().startswith("#"):
            if para:
                blocks.append((start, " ".join(para)))
                para, start = [], None
            if is_row:
                blocks.append((lineno, text))
            continue
        if not para:
            start = lineno
        para.append(text.strip())
    if para:
        blocks.append((start, " ".join(para)))
    return blocks


def _units(path):
    """(line number, sentence) for every sentence on the page."""
    out = []
    for lineno, block in _blocks(path):
        for sentence in SENTENCE_END.split(block):
            if sentence.strip():
                out.append((lineno, sentence.strip()))
    return out


def scan(path, cfg=None):
    """(hits, dead) for one page.

    A hit is ``(line, phrase, sentence)``: a capability negation no ``ABSENT``
    entry covers.  ``dead`` lists the entries for this page that no sentence
    cites any more.
    """
    page = os.path.splitext(os.path.basename(path))[0]
    entries = [e for e in ABSENT if e.page == page]
    fired, hits, seen = set(), [], set()
    for lineno, sentence in _units(path):
        spans = _xslope_spans(sentence)
        ours = bool(spans)
        if len(STATUS_WORDS.findall(sentence)) >= 3:
            continue                                   # the status legend itself
        for pattern, always in NEGATIONS:
            if not (always or ours):
                continue
            for m in re.finditer(pattern, sentence, re.I):
                if NOT_A_CLAIM.search(m.group(0)):
                    continue
                before = sentence[:m.start()]
                if SOURCE_SUBJECT.search(before[-48:]):
                    continue
                # "RS2's own SSRM ... and it does not model ...": a pronoun
                # subject belongs to whichever of the two was named last.
                last_vendor = max((v.end() for v in VENDOR_NAME.finditer(before)),
                                  default=-1)
                last_ours = max((e for _, e in spans if e <= m.start()), default=-1)
                if PRONOUN_SUBJECT.search(before[-24:]) and last_vendor > last_ours:
                    continue
                if WORK_OBJECT.match(sentence[m.end():m.end() + 60]):
                    continue
                covered = [e for e in entries if e.marker.lower() in sentence.lower()]
                if covered:
                    fired.update(covered)
                    continue
                if (lineno, sentence) in seen:         # one report per sentence
                    continue
                seen.add((lineno, sentence))
                hits.append((lineno, m.group(0), sentence))
    dead = [e for e in entries if e not in fired]
    return hits, dead


def run(path, cfg, report=print):
    """Check one page.  Returns the failure count."""
    hits, dead = scan(path, cfg)
    name = os.path.basename(path)
    if not hits and not dead:
        report("  capability: clean")
        return 0
    for lineno, phrase, sentence in hits:
        report(f"  capability: {name}:{lineno} unlisted absence — \"{phrase}\"")
        report(f"              {sentence[:170]}")
        report( "              grep the package before this sentence stands: if "
                "the capability is real, add it to ABSENT in "
                "tools/verification_checks/capabilities.py with the evidence; if "
                "XSLOPE has it, the sentence is the defect")
    for e in dead:
        report(f"  capability: dead entry {e.capability!r} / {e.marker!r} — "
               f"no sentence on {e.page}.md cites it any more")
    return len(hits) + len(dead)


def main(argv):
    from .pages import PAGES
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
    print(f"\n{total} capability problem(s)")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
