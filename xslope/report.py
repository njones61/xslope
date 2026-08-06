# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The Analysis Report: one document model, rendered per format.

An analysis that ends in a screenshot ends in an argument. This module turns a
model and its solutions into a formal document — the thing that gets submitted,
reviewed and filed.

**The content tree comes first.** :func:`build_report` produces a
:class:`Report`: an ordered list of :class:`Section` objects, each holding
:class:`Prose`, :class:`Figure`, :class:`Table`, :class:`KeyValues` and
:class:`Bullets` blocks. It knows nothing about Word, PDF or LaTeX — a figure is
a PNG on disk with a caption and a number, a table is headers and rows of
strings. Every question about *what the report says* is answered here, and is
testable without a renderer.

**Renderers consume the tree.** :mod:`xslope.report_docx` is the first
(``python-docx`` over the styles template in ``xslope/resources``). PDF and
LaTeX renderers read the same tree.

:func:`generate_report` is the public entry point and works headless::

    ok, out = generate_report(slope_data, {"lem": bundle},
                              {"title": "North Levee"}, "report.docx")

``solutions`` maps an engine name to its solved bundle, or to a list of bundles
when several methods were run::

    {"lem": {"slice_df": df, "failure_surface": fs, "results": res,
             "search": search_or_None, "method": "spencer"}}

Figures are rendered by the same plotting functions Studio draws with, at 300 dpi,
so what is in the report is what was on screen. They are embedded in the document
and their files are thrown away with the temporary directory they were drawn in,
unless the caller names a ``figure_dir`` to keep them in. The exception is the
slice force diagrams (:data:`FORCE_DIAGRAMS`), which are drawn by hand rather
than plotted: they ship inside the package and are embedded from there.

Every figure and every table is cited by number in the prose, as a technical
report is read: :func:`cite` writes the phrase and the cross-reference together,
always from the number the counter assigned, and the sentence that carries it is
built where the block is built so that a block the options leave out takes its
citation with it.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from datetime import datetime

#: Figure size, in inches, for the report's full-width plots. Chosen to sit
#: inside a 1 in margin on US Letter portrait without being scaled down.
FIGURE_SIZE = (9.0, 5.5)

#: Rendering resolution for every figure the report embeds.
FIGURE_DPI = 300

#: The method a report falls back to when nobody names one.
DEFAULT_METHOD = "spencer"

#: Display names for the limit-equilibrium methods, in the order the factor of
#: safety summary lists them. Which methods EXIST is read from the solver (see
#: :func:`supported_methods`); this dict only says how they are spelled and in
#: what order, and a method missing from it is still reported, under its key.
METHOD_NAMES = {
    "oms": "Ordinary Method of Slices",
    "bishop": "Bishop's Simplified Method",
    "janbu": "Janbu's Simplified Method",
    "corps": "Corps of Engineers Method",
    "lowe": "Lowe and Karafiath Method",
    "spencer": "Spencer's Method",
    "mprice": "Morgenstern-Price Method",
}

#: What each method assumes and which equilibrium conditions it satisfies, in
#: one paragraph. It opens the method's own Results section so that section
#: stands on its own: a reader who arrives at a slice table or a factor of
#: safety learns what produced it without paging back to the Calculations
#: section, which may be switched off. Each is written from the method's own
#: documentation page (:data:`METHOD_DOC_PAGES`) and says the same thing in the
#: same terms; the Calculations section carries the link, so these do not.
METHOD_SUMMARIES = {
    "oms": (
        "The Ordinary Method of Slices satisfies moment equilibrium of the "
        "sliding mass about the center of the circle, and no other condition. "
        "Interslice forces are neglected entirely and the base normal force is "
        "obtained by resolving the slice forces perpendicular to its own base, "
        "so the factor of safety follows in one pass with no iteration. The "
        "method takes moments and therefore needs a center of rotation: it "
        "applies to circular and composite surfaces only."),
    "bishop": (
        "Bishop's Simplified Method satisfies vertical force equilibrium of "
        "each slice and moment equilibrium of the sliding mass about the center "
        "of the circle. The interslice forces are assumed horizontal, which "
        "neglects the interslice shear. The base normal force that follows from "
        "vertical equilibrium depends on the factor of safety itself, so the "
        "two are solved together by iteration. The method applies to circular "
        "and composite surfaces."),
    "janbu": (
        "Janbu's Simplified Method satisfies vertical force equilibrium of each "
        "slice and horizontal force equilibrium of the sliding mass as a whole; "
        "it does not satisfy moment equilibrium. The interslice shear is "
        "neglected, and the base normal force is the same iterative relation "
        "Bishop's method uses. The factor of safety from that balance is then "
        "multiplied by Janbu's empirical correction factor f_o, which "
        "compensates for the neglected interslice shear. The method applies to "
        "circular and non-circular surfaces."),
    "corps": (
        "The Corps of Engineers Method satisfies force equilibrium — both "
        "horizontal and vertical — of every slice, and does not satisfy moment "
        "equilibrium. The interslice forces are assumed parallel at an "
        "inclination the method's own convention fixes from the geometry, which "
        "here is the ground surface at each slice boundary. The method applies "
        "to circular and non-circular surfaces; because the factor of safety "
        "depends on that assumed inclination, it is best read as a comparison "
        "against a complete-equilibrium method."),
    "lowe": (
        "The Lowe and Karafiath Method satisfies force equilibrium — both "
        "horizontal and vertical — of every slice, and does not satisfy moment "
        "equilibrium. The interslice force at each boundary is inclined at the "
        "average of the ground-surface slope and the slip-surface slope there. "
        "The method applies to circular and non-circular surfaces; tying the "
        "inclination to the base slope makes the result sensitive to the shape "
        "of a non-circular surface, so it too is best read as a comparison "
        "against a complete-equilibrium method."),
    "spencer": (
        "Spencer's Method is a complete-equilibrium procedure: it satisfies "
        "force equilibrium — both horizontal and vertical — and moment "
        "equilibrium of the sliding mass. Its single assumption is that the "
        "interslice forces are parallel, at one constant inclination θ, and the "
        "factor of safety and θ are solved together from those two conditions. "
        "The method applies to circular and non-circular surfaces and is the "
        "recommended basis for design."),
    "mprice": (
        "The Morgenstern-Price Method is a complete-equilibrium procedure: it "
        "satisfies force equilibrium — both horizontal and vertical — and "
        "moment equilibrium of the sliding mass. The interslice inclination is "
        "allowed to vary along the surface as tan θ = λ·f(x) for a chosen "
        "function f, with λ and the factor of safety solved together; Spencer's "
        "method is the special case f(x) = 1. The method applies to circular "
        "and non-circular surfaces."),
}

#: Where the published documentation lives. The same base URL ``mkdocs.yml``
#: declares as ``site_url`` and ``tools/make_corpus_index.py`` writes into the
#: corpus index; that tool is not shipped in the wheel, so the constant is here
#: too, and the two are checked against ``mkdocs.yml`` by the report checks.
DOCS_BASE_URL = "https://xslope.readthedocs.io/en/latest/"

#: The documentation page each limit-equilibrium method is derived on, as a path
#: under ``docs/``. The force-equilibrium methods share one page.
METHOD_DOC_PAGES = {
    "oms": "lem/oms.md",
    "bishop": "lem/bishop.md",
    "janbu": "lem/janbu.md",
    "corps": "lem/force_eq.md",
    "lowe": "lem/force_eq.md",
    "spencer": "lem/spencer.md",
    "mprice": "lem/mprice.md",
}


#: The slice free-body diagram each method's derivation displays as its COMPLETE
#: force set — the figure under "Complete Formulation" on the page
#: :data:`METHOD_DOC_PAGES` names, drawn by hand and kept in ``docs/lem/images``.
#: Methods that share a derivation share a diagram: Janbu's takes the Ordinary
#: Method's, the two force-equilibrium methods take the one their common page
#: draws, and Morgenstern-Price takes Spencer's, which its own page displays and
#: refers every symbol to.
#:
#: The diagram carries every force a slice can take. A method's equation prints
#: the subset the model populates, which is the same relation the derivation
#: pages hold between their figure and their equations.
FORCE_DIAGRAMS = {
    "oms": "oms_complete.png",
    "janbu": "oms_complete.png",
    "bishop": "bishop_complete.png",
    "corps": "slice_fe_complete.png",
    "lowe": "slice_fe_complete.png",
    "spencer": "spencer3_forces.png",
    "mprice": "spencer3_forces.png",
}

#: Printed width of a force diagram, in inches. A hand-drawn slice is about 500
#: pixels across: at a third of the text width it prints at better than 150 dpi
#: and leaves the page it heads room for the equation it explains, where the same
#: drawing stretched to the full 6.5 in would be a page of its own at 90 dpi.
FORCE_DIAGRAM_WIDTH_IN = 3.25


def force_diagram(method):
    """Filesystem path to the slice force diagram for ``method``, or ``""``.

    The PNG ships inside the package (``xslope/resources``), as the input
    template and the report template do: ``docs/`` is not in the wheel, so the
    copy the report embeds cannot be the docs master. The two are kept
    byte-identical by a check in ``run_tests.py``.
    """
    name = FORCE_DIAGRAMS.get(str(method or "").lower())
    if not name:
        return ""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "resources", name)
    return path if os.path.exists(path) else ""


def docs_url(page):
    """``"lem/oms.md"`` -> the published page's URL.

    mkdocs serves a page at its path with the suffix dropped and a trailing
    slash, and an ``index`` page at its directory. One helper, so a report
    section that cites a documentation page — this one now, seepage and finite
    element later — does not each invent the mapping.
    """
    parts = str(page or "").strip("/").split("/")
    if not parts or not parts[0]:
        return DOCS_BASE_URL
    parts[-1] = parts[-1].rsplit(".", 1)[0] if "." in parts[-1] else parts[-1]
    if parts[-1] == "index":
        parts = parts[:-1]
    tail = "/".join(parts)
    return DOCS_BASE_URL + (f"{tail}/" if tail else "")


def method_doc_url(method):
    """The published derivation page for a solver method, or ``""``."""
    page = METHOD_DOC_PAGES.get(str(method).lower())
    return docs_url(page) if page else ""


#: Roughly how many characters of 8.5 pt bold header text fit across the 6.5 in of
#: table a portrait page allows. Divided by the column count it gives the longest
#: word a header may carry: a longer one does not wrap between columns, because
#: Word breaks a word that will not fit rather than widening past the page —
#: "Application" in a ten-column table prints as "Applicatio n".
HEADER_CHARS_PER_ROW = 95


def max_header_word(n_cols):
    """The longest word a header may carry in an ``n_cols`` portrait table."""
    return max(4, int(HEADER_CHARS_PER_ROW / max(1, n_cols)))


# ---------------------------------------------------------------------------
# The content tree
# ---------------------------------------------------------------------------

@dataclass
class Block:
    """Base for everything a section can contain. ``kind`` is what a renderer
    dispatches on, so a renderer never needs ``isinstance``."""

    kind: str = field(init=False, default="block")


@dataclass
class Prose(Block):
    """A paragraph of running text.

    ``links`` is ``[(display text, target), ...]``. A target beginning with ``#``
    is a bookmark elsewhere in the document — the cross-reference from the
    calculations to the slice table — and anything else is a URL. The renderer
    finds each display text in ``text`` and turns that phrase into a link, so the
    sentence is written once, as a sentence.

    ``bold`` is a list of phrases in ``text`` to set in bold, by the same
    find-the-phrase rule: an answer a reader is looking for — the factor of
    safety in the sentence that states it — should be findable without reading
    the sentence.
    """

    text: str
    links: list = field(default_factory=list)
    bold: list = field(default_factory=list)

    def __post_init__(self):
        self.kind = "prose"


@dataclass
class Bullets(Block):
    """An unordered list of short statements."""

    items: list

    def __post_init__(self):
        self.kind = "bullets"


@dataclass
class KeyValues(Block):
    """A label/value block — the shape of every "what was this run" statement.

    Rendered as a two-column borderless table, which reads as a definition list
    in Word and stays a table everywhere else.
    """

    items: list                       # [(label, value), ...], both strings
    title: str = ""

    def __post_init__(self):
        self.kind = "keyvalues"


@dataclass
class Figure(Block):
    """A rendered PNG, its caption, and where it came from.

    ``source`` records which analysis produced it (``"bishop critical surface"``)
    so two reports that differ only in the method they document differ visibly in
    their tree, not only in their pixels.

    ``landscape`` asks the renderer for a rotated page, as a table does — the
    slice-key figure takes one so that it stands beside the slice table it keys
    rather than a page away from it. ``width_in`` of zero means "as wide as the
    page allows", which is how a figure asks to be large without naming a number
    that only holds for one paper size.
    """

    path: str
    caption: str
    number: int = 0
    source: str = ""
    width_in: float = 6.5
    landscape: bool = False

    def __post_init__(self):
        self.kind = "figure"


@dataclass
class Table(Block):
    """Headers, rows of strings, and a caption.

    ``landscape`` asks the renderer for a rotated page — the slice table's whole
    reason for existing. ``legend`` is ``[(term, definition), ...]``, printed
    beneath the table. ``bookmark`` names a Word bookmark placed on the table, so
    a paragraph elsewhere can link to it. ``bold_rows`` holds the indices of the
    rows a reader is meant to find first — the methods a report documents in
    detail, among all the ones it lists.

    ``align`` is per-column justification — ``"l"``, ``"c"`` or ``"r"``, one per
    column, or a single letter for the whole table. It is normally left alone: a
    table that names none takes the report's one policy, read off its own cells
    by :func:`xslope.columns.infer_alignment` — a column of numbers centered, a
    column of anything else ranged left, and every header justified with the
    column it heads. A builder passes it only where the cells cannot say what the
    column is, and says in a comment why.

    ``totals`` is a final row, set in bold and ruled off from the body: the sums
    of the columns that HAVE a sum. A reader who is asked to believe a quotient
    of two sums should be able to find those sums at the foot of the table the
    terms came from.

    ``fit`` is ``"content"`` — size each column to what it holds and let the
    table end where its content ends — or ``"page"``, which stretches the set to
    the full text width. Content is the default: a three-column table ruled
    across a seven-inch page is a table pretending to be a page.
    """

    headers: list
    rows: list
    caption: str = ""
    number: int = 0
    landscape: bool = False
    legend: list = field(default_factory=list)
    bookmark: str = ""
    bold_rows: list = field(default_factory=list)
    align: object = None
    totals: list = field(default_factory=list)
    fit: str = "content"

    def __post_init__(self):
        self.kind = "table"
        if self.align is None:
            from .columns import infer_alignment
            self.align = infer_alignment(self.headers, self.rows)


@dataclass
class Math(Block):
    """A displayed equation, in the notation :mod:`xslope.report_docx` compiles
    to Word's own math (see :func:`xslope.report_docx.omath`).

    Nothing is printed beside the equation. A report numbers its own equations
    consecutively; these are transcribed from several derivations, which number
    theirs independently, and printing those numbers here would run (1), (2),
    (7), (1) down the page with two different equations under the same number.
    So the number lives in the sentence above the equation — "equation (7) of the
    derivation" — where it reads as the reference it is, and the link in the
    section's opening statement is what makes it resolvable.
    """

    notation: str

    def __post_init__(self):
        self.kind = "math"


#: How deep the headings are numbered and styled. A section deeper than this is
#: written in the deepest heading style there is, and numbered with it — the one
#: rule :func:`Report.section_numbers` and :mod:`xslope.report_docx` both follow.
HEADING_LEVELS = 3


@dataclass
class Section:
    """One numbered part of the document: a heading, its blocks, its children.

    ``anchor`` is the bookmark a sentence elsewhere cross-references the section
    by — see :func:`cite_section`. A section that nothing cites leaves it blank
    and is given one when the report is finished, so every heading is reachable.
    """

    title: str
    blocks: list = field(default_factory=list)
    children: list = field(default_factory=list)
    anchor: str = ""

    def walk(self, level=1):
        """This section and every descendant, as ``(level, section)`` pairs."""
        yield level, self
        for child in self.children:
            yield from child.walk(level + 1)


@dataclass
class Report:
    """The whole document: its title-page metadata and its sections."""

    meta: dict = field(default_factory=dict)
    sections: list = field(default_factory=list)

    def section_titles(self):
        """Every heading, in document order, as ``(level, title)`` — what a test
        asserts section order against."""
        out = []
        for s in self.sections:
            out.extend((lvl, sec.title) for lvl, sec in s.walk())
        return out

    def section_numbers(self):
        """Every heading, in document order, as ``(number, level, section)``.

        The number is the one Word will print in front of it — "2", "2.1",
        "2.1.3" — from the multilevel list :mod:`xslope.report_docx` binds to the
        heading styles. Word owns the numbering; this walk exists only so that
        what is written into the document before Word ever opens it — the
        contents page, and the cached result of every cross-reference to a
        section — already says what Word will say. One walk, so the two cannot
        disagree.

        A section deeper than :data:`HEADING_LEVELS` is written in the deepest
        heading style, and is counted at that depth, which is what Word does with
        it.
        """
        out = []
        counts = []
        for root in self.sections:
            for level, sec in root.walk():
                depth = min(level, HEADING_LEVELS)
                del counts[depth:]
                while len(counts) < depth:
                    counts.append(0)
                counts[depth - 1] += 1
                out.append((".".join(str(c) for c in counts), level, sec))
        return out

    def blocks(self, kind=None):
        """Every block in the report, optionally of one ``kind``."""
        out = []
        for s in self.sections:
            for _lvl, sec in s.walk():
                out.extend(b for b in sec.blocks if kind is None or b.kind == kind)
        return out

    def figures(self):
        return self.blocks("figure")

    def tables(self):
        return self.blocks("table")


class _Counter:
    """Sequential figure and table numbers, assigned as blocks are built."""

    def __init__(self):
        self.figure = 0
        self.table = 0

    def next_figure(self):
        self.figure += 1
        return self.figure

    def next_table(self):
        self.table += 1
        return self.table


# ---------------------------------------------------------------------------
# Citations
#
# Every figure and every table is cited by number, in prose that prints whenever
# the block does — the convention a technical report is read under. The number
# always comes from the counter that assigned it, never from a literal, and the
# citation is a live cross-reference: :mod:`xslope.report_docx` bookmarks every
# numbered caption under the name :func:`cite_anchor` builds, so "Figure 4" in a
# sentence jumps to Figure 4.
#
# A block whose number is zero was never assigned one — the options switched it
# off — and cites as nothing, so the sentence that would have named it is written
# without the reference rather than naming Figure 0.
# ---------------------------------------------------------------------------

def cite_anchor(kind, number):
    """The bookmark name a numbered caption carries. ``kind`` is ``"Figure"`` or
    ``"Table"``."""
    return f"xslope_{str(kind).lower()}_{int(number)}"


def cite(kind, number):
    """``(phrase, links)`` for a citation of one numbered block.

    ``phrase`` is ``"Figure 4"`` — what a sentence says — and ``links`` is the
    list a :class:`Prose` carries, so the phrase renders as a cross-reference.
    Both are empty for a number of zero.
    """
    if not number:
        return "", []
    phrase = f"{kind} {int(number)}"
    return phrase, [(phrase, f"#{cite_anchor(kind, number)}")]


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

#: Every option :func:`build_report` reads, with its default. A caller passes a
#: partial dict; anything absent takes the value here. The content toggles are
#: the dialog's checkbox tree, one key per box.
DEFAULT_OPTIONS = {
    # --- title page metadata ---
    "title": "Slope Stability Analysis",
    "project_number": "",
    "organization": "",
    "author": "",
    "date": "",                       # blank = today, formatted long
    "signature_lines": False,         # prepared-by / checked-by (Norm: off by default)

    # --- content ---
    "traceability": True,
    "project_definition": True,
    "pd_figure": True,
    "pd_coords": True,                # label the model figure's geometry points
                                      # with their (x, y); read only when
                                      # pd_figure draws the figure
    "pd_reinforcement": True,
    "pd_units": True,
    "lem": True,
    "lem_materials": True,            # strengths and pore pressures: what the
                                      # method of slices reads off each material
    "lem_loads": True,                # the loads the stability analysis applies
    "lem_search": True,
    "lem_search_figure": True,
    "lem_solution_figure": True,
    "lem_slice_table": True,
    "lem_slice_key": True,
    "lem_calculations": True,
    "lem_rapid": True,
    "seep": True,
    "seep_inputs_figure": True,       # the model as the flow solver reads it
    "seep_materials": True,
    "seep_kr_figure": True,           # every material's unsaturated curve, on
                                      # one axes; drawn only where a material
                                      # carries an unsaturated model
    "seep_mesh_figure": True,         # the mesh with every boundary condition
    "seep_flownet": True,
    "fem": True,
    "fem_inputs_figure": True,        # the model as the FEM solver reads it
    "fem_materials": True,
    "fem_loads": True,                # the loads the deformation analysis carries
    "fem_mesh_figure": True,          # the mesh the section was solved on
    "fem_figure": True,
    "fem_reinforcement": True,        # what the solution put in the bars, and
    "fem_reinforcement_figure": True, # the profiles along the governing ones
    "fem_piles": True,                # the same for the piles; each prints only
    "fem_piles_figure": True,         # where the model carries that member
    "model_checks": False,            # opt-in (Norm: off by default)

    # --- what the report documents ---
    "method": None,                   # which method(s) the detail follows; a name
                                      # or a list of them
    "input_path": None,               # the .xlsx, for the traceability stamp
    "solved_at": None,                # datetime of the solve; None = now
    "style": None,                    # Studio's live display style
    "preflight": None,                # a PreflightReport captured at solve time
    "progress": None,                 # called (done, total, label) per figure
    "dpi": FIGURE_DPI,
    "figsize": FIGURE_SIZE,
}


def resolve_options(options=None):
    """A full options dict: the caller's values over :data:`DEFAULT_OPTIONS`."""
    out = dict(DEFAULT_OPTIONS)
    out.update(options or {})
    return out


# ---------------------------------------------------------------------------
# Reading the model
# ---------------------------------------------------------------------------

def engine_bundles(solutions, engine):
    """The bundles a ``solutions`` mapping carries for one engine, as a list.

    Accepts one bundle or a list of them, so a single run and a sweep are the
    same argument. Every engine is read the same way: the seepage runner emits
    one bundle per boundary condition set and the LEM runner one per method, and
    both arrive here as a list.
    """
    got = (solutions or {}).get(engine)
    if got is None:
        return []
    if isinstance(got, dict):
        return [got]
    return list(got)


def lem_bundles(solutions):
    """The LEM bundles in a ``solutions`` mapping, as a list."""
    return engine_bundles(solutions, "lem")


def seep_bundles(solutions):
    """The seepage bundles, as a list. A bundle is ``{"seep_data", "solution",
    "options"}`` — what Studio's seepage runner emits, once per boundary
    condition set."""
    return engine_bundles(solutions, "seep")


def fem_bundles(solutions):
    """The finite element bundles, as a list. A bundle is ``{"fem_data",
    "solution", "FS", "analysis", "failure_solution"}`` — what Studio's FEM
    runner emits."""
    return engine_bundles(solutions, "fem")


def bundle_method(bundle):
    """The method name a bundle documents. The bundle's own ``method`` wins; the
    solver's ``results['method']`` is the fallback, because not every method
    records its name in its result dict."""
    name = bundle.get("method") or (bundle.get("results") or {}).get("method")
    return str(name).lower() if name else ""


def solved_methods(solutions):
    """Every LEM method present in ``solutions``, in the order they were run."""
    out = []
    for b in lem_bundles(solutions):
        name = bundle_method(b)
        if name and name not in out:
            out.append(name)
    return out


def method_label(name):
    """``"Spencer's Method"`` — the display name, or the bare key if unknown."""
    return METHOD_NAMES.get(str(name).lower(), str(name).title())


def method_summary(name):
    """One paragraph on what a method assumes and what it satisfies, or ``""``.

    See :data:`METHOD_SUMMARIES`. A method the dict does not name gets nothing
    rather than a guess — an unnamed method is still reported, just without a
    sentence claiming things about it.
    """
    return METHOD_SUMMARIES.get(str(name).lower(), "")


def supported_methods():
    """Every limit-equilibrium method the solver offers, in reporting order.

    Read from the template's own method list and confirmed against
    :mod:`xslope.solve`: a name is offered only if the module really defines a
    solver for it, so this list cannot promise a method the engine does not have.
    A method the solver gains appears in the report's summary table without a line
    here changing. :data:`METHOD_NAMES` fixes the order; anything it does not name
    is listed last, under its own key.
    """
    from . import solve
    from .fileio import LEM_METHODS

    names = [n for n in LEM_METHODS
             if n != "all" and callable(getattr(solve, n, None))]
    order = list(METHOD_NAMES)
    names.sort(key=lambda n: order.index(n) if n in order else len(order))
    return names


def method_list(method):
    """The methods a ``method`` option names, as a list of lower-case names.

    A report used to document one method and the option was its name; it now
    documents as many as the dialog's list is ticked for. A bare string is still
    one method — every caller written against the old option keeps working, and
    the list is what everything downstream reads.
    """
    if method is None:
        return []
    if isinstance(method, str):
        return [method.strip().lower()] if method.strip() else []
    out = []
    for name in method:
        name = str(name).strip().lower()
        if name and name not in out:
            out.append(name)
    return out


def featured_methods(solutions, opts=None):
    """The methods the report documents in DETAIL, in the order it documents them.

    The caller's ``method`` option decides, in the order it names them. With
    nothing asked for, the first method that was actually run is featured; with
    nothing run either, the default method is, so a report always has one detail
    block rather than none.
    """
    wanted = method_list((opts or {}).get("method"))
    if wanted:
        return wanted
    run = solved_methods(solutions)
    return [run[0]] if run else [DEFAULT_METHOD]


def select_bundle(solutions, method=None):
    """The bundle whose critical surface the report is built around.

    The named method wins — the first of them, where several are named, since
    they all describe the same surface — and otherwise the first bundle. Returns
    None when there is nothing to select from.
    """
    bundles = lem_bundles(solutions)
    if not bundles:
        return None
    for want in method_list(method):
        for b in bundles:
            if bundle_method(b) == want:
                return b
    return bundles[0]


def search_bundle(solutions):
    """The bundle that carries a search, or None.

    The search belongs to the search, not to a method: it located one surface and
    every method the report features is reported on that one surface, so the
    search is documented once wherever it was run.
    """
    for b in lem_bundles(solutions):
        if b.get("search"):
            return b
    return None


def surface_family(slope_data, solutions=None):
    """``"circular"`` or ``"noncircular"`` for the surface a report documents.

    The dialog offers a method only where the method can run on this surface, and
    the report says so where it cannot; both ask here, so the two answers cannot
    differ.
    """
    bundle = select_bundle(solutions) or {}
    return _surface_family(bundle.get("slice_df"), slope_data)


def file_digest(path):
    """SHA-256 of a file, or ``""`` when it cannot be read."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 16), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return ""


def referenced_materials(slope_data):
    """``[(index, material), ...]`` for the materials the geometry actually uses.

    A materials sheet routinely carries rows left over from an earlier version of
    the model. Geometry is what makes a material real: profile lines and material
    polygons are both numbered in materials order, so a row beyond the last zone
    describes nothing and is left out of the report.
    """
    materials = slope_data.get("materials") or []
    zones = max(len(slope_data.get("profile_lines") or []),
                len(slope_data.get("polygons") or []))
    used = zones if zones else len(materials)
    return [(i + 1, m) for i, m in enumerate(materials[:used])]


def _unit_labels(slope_data):
    from .plot import declared_unit_labels
    return declared_unit_labels(slope_data)


#: Element names by node count, as the mesher and the solvers name them.
ELEMENT_NAMES = {3: "tri3", 6: "tri6", 4: "quad4", 8: "quad8", 9: "quad9"}


def mesh_summary(container):
    """``"1,510 nodes, 706 elements (tri6)"`` for anything carrying a mesh — a
    ``slope_data['mesh']``, a ``seep_data`` or a ``fem_data``, which all hold the
    same three arrays. Empty when it holds no readable mesh."""
    try:
        n_nodes = len(container["nodes"])
        n_elems = len(container["elements"])
        types = sorted({int(t) for t in container["element_types"]})
    except Exception:
        return ""
    kinds = ", ".join(ELEMENT_NAMES.get(t, str(t)) for t in types)
    return f"{n_nodes:,} nodes, {n_elems:,} elements ({kinds})"


def _num(value):
    """A float, or None for anything that is not one."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _populated(rows, key):
    """True when any row carries a meaningful value under ``key``."""
    for r in rows:
        v = r.get(key)
        if isinstance(v, str):
            if v.strip():
                return True
            continue
        n = _num(v)
        if n is not None and abs(n) > 1e-12:
            return True
    return False


def _join(items):
    """``"a, b and c"`` — a list of things, read out."""
    items = [i for i in items if i]
    if len(items) <= 1:
        return items[0] if items else ""
    return f"{', '.join(items[:-1])} and {items[-1]}"


def water_features(slope_data):
    """What the model actually says about water.

    A report must never describe a feature the model does not have: telling the
    reviewer of a bone-dry section how its water loads are derived is not a
    harmless extra sentence, it is a statement about the analysis that is not
    true. Every water sentence and every water row in the report is written from
    this dict, so there is one answer to "is there water in this model?".

    Keys — ``piezo`` (the piezometric lines with a real polyline), ``heads`` (the
    stages carrying specified-head seepage boundaries), ``surfaces`` (the stages a
    water surface can be located for), ``pore`` (the pore-pressure methods the
    referenced materials name, ``"none"`` excluded), and ``any``, true when the
    model carries any of them.
    """
    from .water import water_line_for_stage

    sd = slope_data or {}
    piezo = [name for key, name in (("piezo_line", "Piezometric Line 1"),
                                    ("piezo_line2", "Piezometric Line 2"))
             if len(sd.get(key) or []) >= 2]
    heads = [stage for stage, key in ((1, "seepage_bc"), (2, "seepage_bc2"))
             if ((sd.get(key) or {}).get("specified_heads") or [])]
    surfaces = []
    for stage in (1, 2):
        try:
            line = water_line_for_stage(sd, stage=stage)
        except Exception:
            line = {}
        if line.get("points"):
            surfaces.append(stage)
    pore = sorted({str(m.get("u") or "none").strip().lower()
                   for _i, m in referenced_materials(sd)} - {"none", ""})
    out = {"piezo": piezo, "heads": heads, "surfaces": surfaces, "pore": pore}
    out["any"] = bool(piezo or heads or surfaces or pore)
    return out


def report_analyses(solutions, opts):
    """The analysis types this report documents, in :mod:`xslope.preflight`'s own
    vocabulary — what the model checks are filtered against."""
    out = []
    if opts.get("lem"):
        bundle = select_bundle(solutions, opts.get("method"))
        if bundle is not None:
            results = bundle.get("results") or {}
            out.append("rapid" if "stage1_FS" in results else "lem")
    if opts.get("seep") and seep_bundles(solutions):
        out.append("seep")
    for bundle in (fem_bundles(solutions) if opts.get("fem") else []):
        name = "ssrm" if str(bundle.get("analysis")) == "ssrm" else "fem"
        if name not in out:
            out.append(name)
    return out


def relevant_findings(findings, analyses):
    """The findings that concern the analyses a report contains.

    A preflight report is captured for the run it gated, and a model checked for
    one engine can carry findings about another — a mesh the finite element engine
    would refuse says nothing about a limit equilibrium report. Each rule already
    declares the analyses it applies to (``@rule(..., analyses=("fem",))``), so
    relevance is read straight off the registry rather than guessed from the
    message: a finding is kept when its rule applies to any analysis in the report,
    and a rule tagged ``("*",)`` applies to all of them. A finding whose rule is not
    in the registry is kept — an id this build cannot resolve is not evidence of
    irrelevance.
    """
    from .preflight import expand_analysis, rules

    findings = list(findings or [])
    if not analyses:
        return findings
    wanted = set()
    for name in analyses:
        try:
            wanted |= set(expand_analysis(name))
        except ValueError:
            continue
    if not wanted:
        return findings
    by_id = {r.id: r for r in rules()}
    return [f for f in findings
            if by_id.get(f.rule_id) is None or by_id[f.rule_id].applies_to(wanted)]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _new_figure(figsize):
    """A standalone Agg figure — no pyplot, so the builder is safe off the main
    thread and leaks no global figure manager state."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MplFigure
    fig = MplFigure(figsize=figsize)
    FigureCanvasAgg(fig)
    return fig


def _render(draw, path, opts):
    """Draw into a fresh figure and write it to ``path``. Returns the path, or
    None when the plot could not be produced (a missing optional input must not
    take the whole report down).

    The directory is made here rather than up front: a report with every figure
    switched off should leave no folder behind to explain.
    """
    fig = _new_figure(opts["figsize"])
    try:
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        draw(fig)
        fig.savefig(path, dpi=opts["dpi"], bbox_inches="tight")
    except Exception:
        import traceback
        traceback.print_exc()
        return None
    return path


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def _traceability_section(slope_data, solutions, opts):
    """The "could someone reproduce this?" stamp."""
    from ._version import __version__

    items = [("xslope version", __version__)]
    path = opts.get("input_path")
    if path:
        items.append(("Input file", os.path.basename(path)))
        digest = file_digest(path)
        if digest:
            items.append(("Input file SHA-256", digest))
    else:
        items.append(("Input file", "not saved to a file"))

    solved = opts.get("solved_at") or datetime.now()
    items.append(("Analysis run", solved.strftime("%Y-%m-%d %H:%M")))

    mesh = slope_data.get("mesh")
    if mesh is not None:
        items.append(("Mesh", mesh_summary(mesh) or "present"))

    items.append(("Report generated", datetime.now().strftime("%Y-%m-%d %H:%M")))

    sec = Section("Traceability")
    sec.blocks.append(Prose(
        "The results in this report were computed by xslope from the input file "
        "below. Its SHA-256 digest identifies that file exactly."))
    sec.blocks.append(KeyValues(items))
    return sec


def _property_table(slope_data, fields, caption, counter):
    """One row per referenced material, one column per property that is real.

    ``fields`` is ``(key, header, formatter, always-shown)``. A column marked
    always-shown is printed whether or not the materials carry a value for it —
    it is part of what the analysis is — and every other column is printed only
    where some material populates it, so a table never rules off a column of
    blanks. Each engine's table is this function with its own field list.
    """
    rows_src = referenced_materials(slope_data)
    if not rows_src:
        return None
    mats = [m for _i, m in rows_src]
    keep = [f for f in fields if f[3] or _populated(mats, f[0])]
    headers = [f[1] for f in keep]
    rows = []
    for idx, m in rows_src:
        rows.append([str(idx) if key == "__index" else fmt(m)
                     for key, _h, fmt, _always in keep])
    return Table(headers, rows, caption, counter.next_table())


def _unit_suffix(slope_data):
    """``unit("stress")`` -> ``" (kPa)"``, and ``""`` for a model that declares no
    unit system — the suffix every property header is built with."""
    lbl = _unit_labels(slope_data)

    def unit(key):
        return f" ({lbl[key]})" if lbl and lbl.get(key) else ""

    return unit


def _fmt(value, spec="{:.2f}"):
    n = _num(value)
    return "" if n is None else spec.format(n)


def _materials_table(slope_data, counter):
    """The materials table: referenced rows, populated columns."""
    unit = _unit_suffix(slope_data)

    # (key, header, formatter, always-shown)
    fields = [
        ("__index", "#", lambda m: "", True),
        ("name", "Material", lambda m: str(m.get("name") or ""), True),
        ("gamma", f"γ{unit('unit_weight')}", lambda m: _fmt(m.get("gamma"), "{:.1f}"), True),
        ("gamma_sat", f"γ_sat{unit('unit_weight')}", lambda m: _fmt(m.get("gamma_sat"), "{:.1f}"), False),
        ("option", "Strength", lambda m: str(m.get("option") or "").upper(), True),
        ("c", f"c{unit('stress')}", lambda m: _fmt(m.get("c"), "{:.1f}"), True),
        ("phi", "φ (deg)", lambda m: _fmt(m.get("phi"), "{:.1f}"), True),
        ("cp", "c/p", lambda m: _fmt(m.get("cp"), "{:.3f}"), False),
        ("r_elev", f"r elev.{unit('length')}", lambda m: _fmt(m.get("r_elev"), "{:.1f}"), False),
        ("d", f"d{unit('stress')}", lambda m: _fmt(m.get("d"), "{:.1f}"), False),
        ("psi", "ψ (deg)", lambda m: _fmt(m.get("psi"), "{:.1f}"), False),
        ("t_cut", f"σ_t{unit('stress')}", lambda m: _fmt(m.get("t_cut"), "{:.1f}"), False),
        ("ru", "r_u", lambda m: _fmt(m.get("ru"), "{:.3f}"), False),
        ("u", "Pore pressure", lambda m: str(m.get("u") or "none"), True),
    ]
    return _property_table(slope_data, fields, "Material properties", counter)


def _seep_materials_table(slope_data, counter):
    """The conductivities and unsaturated parameters the flow domain is solved
    with. A material's strength has no bearing on its flow, so this is a table of
    its own rather than more columns on the strength one."""
    unit = _unit_suffix(slope_data)
    models = {"lf": "linear front", "vg": "van Genuchten", "gard": "Gardner"}
    fields = [
        ("__index", "#", lambda m: "", True),
        ("name", "Material", lambda m: str(m.get("name") or ""), True),
        ("k1", f"k₁{unit('k')}", lambda m: _fmt(m.get("k1"), "{:.3g}"), True),
        ("k2", f"k₂{unit('k')}", lambda m: _fmt(m.get("k2"), "{:.3g}"), True),
        ("alpha", "α (deg)", lambda m: _fmt(m.get("alpha"), "{:.1f}"), True),
        ("unsat", "Unsaturated model",
         lambda m: models.get(str(m.get("unsat") or "").strip().lower(), ""), False),
        ("kr0", "k_r0", lambda m: _fmt(m.get("kr0"), "{:.4g}"), False),
        ("h0", f"h₀{unit('length')}", lambda m: _fmt(m.get("h0"), "{:.2f}"), False),
        ("vg_a", "a", lambda m: _fmt(m.get("vg_a"), "{:.4g}"), False),
        ("vg_n", "n", lambda m: _fmt(m.get("vg_n"), "{:.3f}"), False),
    ]
    return _property_table(slope_data, fields, "Seepage material properties",
                           counter)


def _kr_materials(slope_data):
    """The materials the unsaturated conductivity figure is drawn from — empty
    for a model analyzed saturated throughout, which has no curve to draw.

    Read by the section that draws the figure and by the count that promises it,
    so the two cannot disagree about whether there is one.
    """
    from .plot import material_kr_curves
    mats = [m for _i, m in referenced_materials(slope_data)]
    return mats if material_kr_curves(mats) else []


def _fem_materials_table(slope_data, counter):
    """The strength and stiffness the finite element analysis is solved with. The
    stiffness columns are what separate it from the limit equilibrium table: a
    limit equilibrium analysis never reads E or ν, and this one does."""
    unit = _unit_suffix(slope_data)
    fields = [
        ("__index", "#", lambda m: "", True),
        ("name", "Material", lambda m: str(m.get("name") or ""), True),
        ("gamma", f"γ{unit('unit_weight')}", lambda m: _fmt(m.get("gamma"), "{:.1f}"), True),
        ("c", f"c{unit('stress')}", lambda m: _fmt(m.get("c"), "{:.1f}"), True),
        ("phi", "φ (deg)", lambda m: _fmt(m.get("phi"), "{:.1f}"), True),
        ("E", f"E{unit('stress')}", lambda m: _fmt(m.get("E"), "{:,.0f}"), True),
        ("nu", "ν", lambda m: _fmt(m.get("nu"), "{:.2f}"), True),
        ("t_cut", f"σ_t{unit('stress')}", lambda m: _fmt(m.get("t_cut"), "{:.1f}"), False),
    ]
    return _property_table(slope_data, fields, "Finite element material properties",
                           counter)


#: What each pore-pressure option on the materials sheet MEANS, in the words the
#: prose uses. The sheet's own codes are printed in the table, against the
#: material each belongs to; a sentence that repeated the code would be telling
#: the reader what they can already see instead of what it is.
PORE_SOURCES = {
    "piezo": "the piezometric line",
    "seep": "the computed seepage field",
    "ru": "a pore pressure ratio r_u applied to the overburden",
    "cons": "the consolidation ratio entered for the material",
}


def _water_items(slope_data, feats):
    """Where the pore pressures a stability analysis reads come from, as rows.

    The seepage boundary conditions are NOT here: they are an input to the flow
    problem and are stated in the section that solves it. What a limit
    equilibrium analysis needs to know is which water surface it stands under,
    how water above the ground surface is loaded onto it, and how each material's
    pore pressure is taken. Empty for a model that carries none of that, and the
    caller says so in a sentence instead.
    """
    from .water import water_line_for_stage, water_loads_mode

    items = []
    if feats["surfaces"]:
        mode = water_loads_mode(slope_data)
        items.append(("Water loads", (
            "derived by the engine from the model's own water surface (automatic)"
            if mode == "auto" else
            "taken from the distributed loads as entered (manual)")))

    for key, name in (("piezo_line", "Piezometric Line 1"),
                      ("piezo_line2", "Piezometric Line 2")):
        pts = slope_data.get(key) or []
        if len(pts) >= 2:
            items.append((name, f"{len(pts)} points, elevation "
                                f"{min(p[1] for p in pts):g} to "
                                f"{max(p[1] for p in pts):g}"))

    for stage in (1, 2):
        if stage in feats["surfaces"]:
            line = water_line_for_stage(slope_data, stage=stage)
            items.append((f"Water surface (stage {stage})",
                          f"from {line['source']}"))
    return items


def _loads_table(slope_data, counter):
    """The distributed loads as entered, one row per DEFINING POINT.

    A load block is a polyline of (x, y, pressure) points, and all three vary
    along it: a load on a sloping face changes elevation point by point, and a
    trapezoidal load changes pressure. Reporting a block as its x range and its
    largest pressure described a rectangle standing on level ground — one special
    case of what the sheet can hold, and not the interesting one. The points are
    what the user entered and what the engine integrates, so the points are what
    the report prints.
    """
    blocks = slope_data.get("dloads") or []
    if not blocks:
        return None
    lbl = _unit_labels(slope_data)
    su = f" ({lbl['stress']})" if lbl and lbl.get("stress") else ""
    lu = f" ({lbl['length']})" if lbl and lbl.get("length") else ""
    dirs = slope_data.get("dload_dirs") or []
    headers = ["Load", "Point", f"x{lu}", f"y{lu}", f"Pressure{su}", "Direction"]
    rows = []
    for i, blk in enumerate(blocks):
        # The load number and its direction are properties of the BLOCK, printed
        # once against its first point rather than repeated down every row.
        for j, point in enumerate(blk):
            rows.append([
                str(i + 1) if j == 0 else "",
                str(j + 1),
                _fmt(point.get("X"), "{:g}"),
                _fmt(point.get("Y"), "{:g}"),
                _fmt(point.get("Normal"), "{:g}"),
                str(dirs[i] if i < len(dirs) else "normal") if j == 0 else "",
            ])
    return Table(headers, rows, "Distributed loads", counter.next_table())


def _loads_section(slope_data, feats, counter, seismic=True, already=0):
    """The loads an engine applies, as a section of that engine's own inputs.

    A distributed load is not a property of the section: it is something an
    analysis puts on it, and the two engines put it on differently — a limit
    equilibrium analysis integrates it along the ground surface of each slice, a
    finite element analysis carries it as a traction on the boundary. So each
    engine presents the loads it applies, from this one builder, rather than a
    general table that belongs to neither.

    ``seismic`` prints the pseudo-static coefficient, which only the limit
    equilibrium analysis applies. ``already`` is the number of a loads table an
    earlier section has printed: the blocks are identical — same points, same
    pressures — so the second engine points at the first table rather than
    setting the same numbers twice, where two copies could disagree.
    """
    sub = Section("Loads")
    if already:
        where, links = cite("Table", already)
        sub.blocks.append(Prose(
            f"The analysis carries the distributed loads of {where}, applied as "
            f"tractions on the boundary of the mesh.", links=links))
        return sub

    table = _loads_table(slope_data, counter)
    if table is not None:
        where, links = cite("Table", table.number)
        sub.blocks.append(Prose(
            f"Each distributed load is entered as a polyline whose points carry "
            f"a position and a pressure, and {where} gives those points as "
            f"entered. The load is integrated along the ground surface between "
            f"them.", links=links))
        sub.blocks.append(table)
    elif feats["surfaces"]:
        sub.blocks.append(Prose(
            "The model carries no distributed loads entered by hand. Any water "
            "standing on the section is measured by the engine from the water "
            "surface and applied as a distributed load."))
    else:
        sub.blocks.append(Prose("The model carries no distributed loads."))
    k = _num(slope_data.get("k_seismic"))
    if seismic and k:
        sub.blocks.append(Prose(
            f"A pseudo-static seismic coefficient of k = {k:g} is applied, "
            f"acting horizontally toward the toe on every slice."))
    return sub


def _reinforcement_table(slope_data, counter):
    lines = slope_data.get("reinforcement_lines") or []
    if not lines:
        return None
    lbl = _unit_labels(slope_data)
    fu = f" ({lbl['force_per_len']})" if lbl and lbl.get("force_per_len") else ""
    lu = f" ({lbl['length']})" if lbl and lbl.get("length") else ""
    headers = ["Label", "Start (x, y)", "End (x, y)", f"T_max{fu}", f"T_res{fu}",
               f"L_p1{lu}", f"L_p2{lu}", f"Spacing{lu}", "Direction", "Applied"]
    rows = []
    for ln in lines:
        rows.append([
            str(ln.get("label") or ""),
            f"({_fmt(ln.get('x1'), '{:g}')}, {_fmt(ln.get('y1'), '{:g}')})",
            f"({_fmt(ln.get('x2'), '{:g}')}, {_fmt(ln.get('y2'), '{:g}')})",
            _fmt(ln.get("t_max"), "{:g}"), _fmt(ln.get("t_res"), "{:g}"),
            _fmt(ln.get("lp1"), "{:g}"), _fmt(ln.get("lp2"), "{:g}"),
            _fmt(ln.get("spacing"), "{:g}"),
            str(ln.get("dir") or "tangent"), str(ln.get("appl") or "active"),
        ])
    return Table(headers, rows, "Reinforcement lines", counter.next_table())


def _piles_table(slope_data, counter):
    piles = slope_data.get("pile_lines") or []
    if not piles:
        return None
    lbl = _unit_labels(slope_data)
    fu = f" ({lbl['force_per_len']})" if lbl and lbl.get("force_per_len") else ""
    lu = f" ({lbl['length']})" if lbl and lbl.get("length") else ""
    headers = ["Label", "Top (x, y)", "Bottom (x, y)", f"H{fu}", "θ (deg)",
               f"D{lu}", f"Spacing{lu}", "Applied"]
    rows = []
    for p in piles:
        rows.append([
            str(p.get("label") or ""),
            f"({_fmt(p.get('x1'), '{:g}')}, {_fmt(p.get('y1'), '{:g}')})",
            f"({_fmt(p.get('x2'), '{:g}')}, {_fmt(p.get('y2'), '{:g}')})",
            _fmt(p.get("H"), "{:g}") or "computed (Ito and Matsui)",
            _fmt(p.get("theta_p"), "{:g}"),
            _fmt(p.get("D_pile"), "{:g}"), _fmt(p.get("S"), "{:g}"),
            str(p.get("appl") or "active"),
        ])
    return Table(headers, rows, "Piles", counter.next_table())


def _units_prose(slope_data):
    """What the report's numbers are in — including the unit weight of water,
    which is a constant of the units and not a description of any water the model
    carries: a dry section is still analyzed in a system that has one."""
    from .units import normalize_unit_system
    system = normalize_unit_system(slope_data.get("unit_system"))
    lbl = _unit_labels(slope_data) or {}
    gw = _num(slope_data.get("gamma_water"))
    water = ""
    if gw is not None:
        suffix = f" {lbl['unit_weight']}" if lbl.get("unit_weight") else ""
        # A constant of the unit system, not a description of any water the
        # model carries: a dry section is still analyzed in a system that has
        # one, and a sentence here about water loads would be a statement about
        # an analysis that has none.
        water = f" The unit weight of water is {gw:g}{suffix}."
    if system is None:
        return Prose(
            "The model declares no unit system. Every quantity is in the units "
            "the inputs were entered in, and results are consistent with those "
            "units throughout." + water)
    name = "SI" if system == "si" else "US customary"
    return Prose(
        f"All quantities are in {name} units: lengths in {lbl.get('length', '')}, "
        f"stresses and pressures in {lbl.get('stress', '')}, unit weights in "
        f"{lbl.get('unit_weight', '')}, and forces per unit thickness of section "
        f"in {lbl.get('force_per_len', '')}. Angles are in degrees. The analysis "
        f"is two-dimensional: every force is per unit thickness normal to the "
        f"section." + water)


def _project_definition_section(slope_data, opts, counter, figure_dir,
                                progress=None):
    """The model itself: the section, what it is made of, and what is built into
    it.

    Only what more than one engine reads belongs here. A material's shear
    strength is a limit equilibrium and finite element input and its conductivity
    is a seepage input, so each engine gives its own materials table; the loads
    are applied by the analysis rather than owned by the section; the pore
    pressures a stability analysis reads are stated where that analysis is
    documented. What is left is the section, the units, and the reinforcement and
    piles, which are structure the section carries whichever engine looks at it.
    """
    sec = Section("Project Definition")
    feats = water_features(slope_data)

    # The figure is drawn first so the prose that introduces it can be written
    # from what was actually produced — a report never announces a figure that
    # the option switched off or the plot failed to render.
    figure = None
    # Point coordinates are a property of the figure, so the option is read only
    # where the figure is drawn: switching the labels off never switches the
    # figure off, and switching the figure off never consults the labels.
    coords = bool(opts.get("pd_coords", DEFAULT_OPTIONS["pd_coords"]))
    if opts["pd_figure"]:
        path = os.path.join(figure_dir, "model.png")

        def draw(fig):
            from .plot import plot_inputs
            plot_inputs(slope_data, fig=fig, mode="shared", show_title=False,
                        frame="content", style=opts.get("style"),
                        label_coordinates=coords)

        if progress:
            progress("the analysis model")
        if _render(draw, path, opts):
            figure = Figure(path, "Analysis model", counter.next_figure(),
                            source="shared model")

    geometry = ("material zone polygons"
                if slope_data.get("polygons") and not slope_data.get("profile_lines")
                else "profile lines")
    mats = referenced_materials(slope_data)
    text = (f"The section is defined by {len(mats)} material "
            f"{'zone' if len(mats) == 1 else 'zones'} described with {geometry}.")
    # Named, in the order the model lists them, so the reader meets the words the
    # legend and every properties table use before either appears.
    names = [str(m.get("name") or "").strip() for _i, m in mats]
    names = [n for n in names if n]
    if names:
        text += (f" The {'zone is' if len(names) == 1 else 'zones are'} "
                 f"{_join(names)}.")
    links = []
    if figure is not None:
        # Only what the model carries: a figure caption that lists water surfaces
        # on a dry section describes a different model.
        shows = ["the geometry and materials"]
        if feats["surfaces"]:
            shows.append("the water surface"
                         if len(feats["surfaces"]) == 1 else "the water surfaces")
        if slope_data.get("dloads"):
            shows.append("the distributed loads")
        if slope_data.get("reinforcement_lines"):
            shows.append("the reinforcement lines")
        if slope_data.get("pile_lines"):
            shows.append("the piles")
        if coords:
            shows.append("every geometry point labeled with its coordinates")
        where, links = cite("Figure", figure.number)
        text += (f" Every analysis in this report is run on the model of "
                 f"{where}: {_join(shows)}.")
    sec.blocks.append(Prose(text, links=links))

    # The units statement leads: a reader meets the numbers knowing what they are
    # in, rather than finding out at the end of the section.
    if opts["pd_units"]:
        sec.blocks.append(_units_prose(slope_data))

    if figure is not None:
        sec.blocks.append(figure)

    # Reinforcement and piles are separate features with separate tables, and a
    # section for one is not a heading the other hides under.
    if opts["pd_reinforcement"]:
        reinf = _reinforcement_table(slope_data, counter)
        if reinf is not None:
            where, links = cite("Table", reinf.number)
            sec.children.append(Section("Reinforcement", [
                Prose(f"{where} gives each reinforcement line: its endpoints, "
                      f"the tensile and residual capacities it can develop, the "
                      f"pullout lengths at either end, its out-of-plane spacing, "
                      f"and whether it acts as an active or a passive force.",
                      links=links),
                reinf]))
        piles = _piles_table(slope_data, counter)
        if piles is not None:
            where, links = cite("Table", piles.number)
            sec.children.append(Section("Piles", [
                Prose(f"{where} gives each pile: its head and tip, the lateral "
                      f"force it carries, its inclination, its diameter and "
                      f"out-of-plane spacing, and whether it acts as an active "
                      f"or a passive force. A pile that states no force has one "
                      f"computed for it by the Ito and Matsui method.",
                      links=links),
                piles]))

    return sec


def _search_section(slope_data, bundle, opts, counter, figure_dir, method,
                    progress=None):
    """What the search did, and the surfaces it tried."""
    search = bundle.get("search")
    if not search:
        return None
    sub = Section("Search for the Critical Surface")
    fs_cache = search.get("fs_cache") or []
    path_pts = search.get("search_path") or []
    kind = search.get("kind", "circular")

    items = [("Surface family", "circular" if kind == "circular" else "non-circular"),
             ("Method", method_label(method)),
             ("Trial surfaces evaluated", f"{len(fs_cache):,}")]
    if path_pts:
        items.append(("Refinement stages", str(len(path_pts))))
    valid = [c for c in fs_cache if _num(c.get("FS")) is not None
             and _num(c.get("FS")) < 9000]
    if valid:
        items.append(("Factor of safety range",
                      f"{min(_num(c['FS']) for c in valid):.3f} to "
                      f"{max(_num(c['FS']) for c in valid):.3f}"))

    window = (slope_data.get("search_window") or {})
    stated = {k: v for k, v in window.items() if v is not None}
    if stated:
        items.append(("Search window",
                      ", ".join(f"{k} = {v:g}" if isinstance(v, (int, float))
                                else f"{k} = {v}" for k, v in sorted(stated.items()))))
    else:
        items.append(("Search window", "none declared; the search was unconstrained"))

    # The figure is drawn before the paragraph that introduces it, so the
    # sentence naming the surface can name the figure it is highlighted on.
    figure = None
    if opts["lem_search_figure"]:
        fpath = os.path.join(figure_dir, "search.png")

        def draw(fig):
            from .plot import (plot_circular_search_results,
                               plot_noncircular_search_results)
            if kind == "circular":
                plot_circular_search_results(
                    slope_data, fs_cache, search_path=path_pts,
                    circle_cache=search.get("circle_cache"), fig=fig,
                    show_title=False, style=opts.get("style"))
            else:
                plot_noncircular_search_results(
                    slope_data, fs_cache, search_path=path_pts, fig=fig,
                    show_title=False, style=opts.get("style"))

        if progress:
            progress("the trial surfaces the search evaluated")
        if _render(draw, fpath, opts):
            figure = Figure(
                fpath, "Trial surfaces evaluated by the search, with the critical "
                       "surface highlighted", counter.next_figure(),
                source=f"{method} search")

    where, links = cite("Figure", figure.number if figure is not None else 0)
    highlighted = f", highlighted in {where}," if where else ""
    sub.blocks.append(Prose(
        f"The critical surface was located by automated search, which refines a "
        f"grid of trial surfaces until the factor of safety stops improving. The "
        f"reported surface{highlighted} is the lowest the search reached.",
        links=links))
    sub.blocks.append(KeyValues(items))
    if figure is not None:
        sub.blocks.append(figure)
    return sub


def _rapid_section(results, counter):
    """The three-stage rapid drawdown detail, when the run was one."""
    if not results or "stage1_FS" not in results:
        return None
    sub = Section("Rapid Drawdown")
    sub.blocks.append(Prose(
        "The surface was analyzed for rapid drawdown by the three-stage procedure "
        "of Duncan, Wright and Wong. Stage 1 establishes the consolidation "
        "stresses under the full pool, stage 2 applies undrained strengths to the "
        "drawn-down section, and stage 3 re-checks the same section with drained "
        "strengths where those are lower. The reported factor of safety is the "
        "lower of stages 2 and 3."))
    items = [("Stage 1 — full pool, drained", f"{results['stage1_FS']:.3f}"),
             ("Stage 2 — drawn down, undrained", f"{results['stage2_FS']:.3f}"),
             ("Stage 3 — drawn down, drained", f"{results['stage3_FS']:.3f}"),
             ("Governing factor of safety", f"{results['FS']:.3f}")]
    sub.blocks.append(KeyValues(items))
    return sub


def _surface_family(slice_df, slope_data):
    """``"circular"`` or ``"noncircular"`` for the surface a slice table describes.

    Read from the slice table's own circle center where there is one, because that
    is what the moment methods actually test; the model's declared family is the
    fallback for a bundle that carries no slices.
    """
    if slice_df is not None and len(slice_df) and "r" in slice_df.columns:
        import pandas as pd
        r = slice_df["r"].iloc[0]
        return "noncircular" if r is None or pd.isna(r) else "circular"
    return "circular" if (slope_data or {}).get("circular") else "noncircular"


def _solve_on(name, slice_df, rapid=False):
    """``(slice_df, results)`` for one method on an already-built slice table, or
    ``(None, None)``.

    ``solve_selected`` is the solver's own entry point, prints its result and
    returns the error string rather than a dict when a method does not converge;
    the printing is swallowed here because a report is not a console session, and
    a non-dict return is what "did not converge" means. The slice table is copied:
    the solvers write their working columns into it, and one method's own table
    must not pick up another method's arithmetic. The copy is handed back because
    those working columns — ``n_eff`` above all — are exactly what a slice table
    and a calculations section print.
    """
    import contextlib
    import io

    from .solve import solve_selected
    work = slice_df.copy()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = solve_selected(name, work, rapid=rapid)
    except Exception:
        return None, None
    if not isinstance(out, dict):
        return None, None
    return work, dict(out, method=out.get("method") or name)


def _solve_for_summary(name, slice_df, rapid=False):
    """One method's answer on an already-built slice table, or None."""
    _df, res = _solve_on(name, slice_df, rapid)
    return res


def detail_bundle(slope_data, solutions, method):
    """``(bundle, note)`` — what the report documents in detail for one method.

    A method that was RUN is documented from its own bundle, and ``note`` is
    empty. A method that was not is solved here, on the critical surface the
    report documents, exactly as the factor of safety summary solves it: the same
    surface, the same slice geometry, a different method. ``note`` is then the
    sentence that says so, because a reader must never have to work out whether a
    block came from the run or from the report.

    ``(None, reason)`` where the method cannot be documented at all — it does not
    apply to this surface family, or it does not converge on it. Either way the
    note is a whole sentence, printed as it stands: the heading above it already
    names the method, and the checker's own refusal names it too, so a wrapper
    sentence would only say the name a third time.
    """
    from .preflight import method_surface_reason

    method = str(method or "").lower()
    label = method_label(method)
    for b in lem_bundles(solutions):
        if bundle_method(b) == method:
            return b, ""
    if method not in supported_methods():
        return None, f"{label} is not a method xslope offers."

    base = select_bundle(solutions)
    base_df = (base or {}).get("slice_df")
    if base_df is None or not len(base_df):
        return None, (f"There is no slice table to work {label} through on this "
                      f"surface.")
    family = _surface_family(base_df, slope_data)
    reason = method_surface_reason(method, family)
    if reason:
        return None, reason[0].upper() + reason[1:]
    rapid = "stage1_FS" in (base.get("results") or {})
    df, res = _solve_on(method, base_df, rapid)
    if df is None:
        return None, (f"{label} did not converge on this surface, so no detail "
                      f"is reported for it.")
    return ({"slice_df": df, "failure_surface": base.get("failure_surface"),
             "results": res, "search": None, "method": method},
            "It was not run in the analysis; the report solved it on the same "
            "critical surface, so the comparison is between methods rather than "
            "surfaces.")


def _solution_parameters(res):
    """The method-specific numbers that go beside a factor of safety."""
    extra = []
    if _num(res.get("theta")) is not None:
        extra.append(f"θ = {_num(res['theta']):.2f} deg")
    if _num(res.get("fo")) is not None:
        extra.append(f"f₀ = {_num(res['fo']):.3f}")
    if _num(res.get("lambda")) is not None:
        extra.append(f"λ = {_num(res['lambda']):.3f}")
    return ", ".join(extra)


def _fs_table(slope_data, solutions, opts, counter):
    """The factor of safety each documented method reported, on its own surface.

    Only the methods this report documents, and only their own answers. The
    summary used to list every method xslope offers, filling in the ones that had
    not been run by solving them on the surface the featured method had found —
    which made a column of numbers that were not comparable in the way a column
    of numbers reads as: each method has its own critical surface, and a method
    solved on somebody else's is not reporting its factor of safety. Which
    methods a reader wants compared is a decision the reader makes, by choosing
    what the report documents.

    None for a report of a single method: two rows are a comparison, one row is
    the same number the method's own section states, and a table that restates it
    is a table that can disagree with it.
    """
    methods = featured_methods(solutions, opts)
    if len(methods) < 2:
        return None

    solved = {}
    for b in lem_bundles(solutions):
        name = bundle_method(b)
        if name and name not in solved:
            solved[name] = b

    rows = []
    for name in methods:
        bundle = solved.get(name) or {}
        res = bundle.get("results") or {}
        fs = _num(res.get("FS"))
        rows.append([method_label(name),
                     "did not converge" if fs is None else f"{fs:.3f}",
                     _surface_provenance(bundle).capitalize(),
                     "" if fs is None else _solution_parameters(res)])

    if not rows:
        return None
    return Table(["Method", "Factor of safety", "Surface", "Solution parameters"],
                 rows, "Computed factors of safety", counter.next_table())


def _surface_provenance(bundle):
    """How the surface a factor of safety belongs to was arrived at.

    "Critical" is a word a search earns. A factor of safety computed on a surface
    the user entered is the factor of safety of THAT surface and not the minimum
    over any family, and calling it critical invites a reader to take it for one.
    """
    return ("the critical surface this method searched for"
            if (bundle or {}).get("search")
            else "the surface specified in the input")


# ---------------------------------------------------------------------------
# The calculation
#
# The section prints the equation the SOLVER evaluates, with the converged
# numbers in it. Two things keep it from drifting into a textbook:
#
# 1. Every per-slice term is built here from the columns the solver wrote —
#    ``n_eff`` above all, which is the base normal at the converged factor of
#    safety — and the sums are formed from those arrays. Nothing is re-derived.
# 2. The quotient is evaluated and compared with the solver's own factor of
#    safety before anything is printed (:func:`_closes`). A model whose terms
#    this module does not carry fails that comparison and gets no calculation
#    rather than a wrong one. The comparison is made against the tolerance the
#    solver itself converged to, so a solution that came as close as it was
#    asked to is never refused.
#
# The equations themselves are the ones on the method's documentation page, in
# that page's own symbols, so the report and the derivation can be read side by
# side. Where the code and the page differ the code wins and the page is the bug.
# The three methods whose page publishes a march rather than a quotient are the
# exception, and their section says what its equation is instead of citing the
# page for it — see :data:`WHOLE_MASS_BALANCE_METHODS`.
# ---------------------------------------------------------------------------

#: The Word bookmark placed on the slice table, and linked to from the
#: calculations — the per-slice terms of every sum are columns of that table.
SLICE_TABLE_BOOKMARK = "xslope_slice_table"

# ---------------------------------------------------------------------------
# THE FORCE REGISTRY
#
# One table of the forces a slice can carry, and what each of them contributes
# to every equation this module prints. The equations are many — two moment
# sums, a horizontal balance, Spencer's two force sums, the base-normal equation
# — and each carries a different subset of the forces in a different notation, so
# the term set used to be five parallel hand-kept lists plus the arrays, the
# omission sentence, the passive gate and the nomenclature. A force added to some
# of them and not the rest is printed in one equation and denied in the next,
# which is a defect the reader can see and the code could not.
#
# Every force is one :class:`ForceTerm`, and every equation is assembled from the
# same entries. A contribution is either a :class:`Term` — the sign, how it is
# written, and the per-slice values it is formed from — or a
# :class:`NotApplicable` carrying the reason, which is a statement about the
# mechanics rather than a gap: Spencer's F_h and F_v exclude the base normal
# because equations (1) and (2) are defined that way, and the shear on the top of
# the slice is in those published equations and in no xslope model.
# ---------------------------------------------------------------------------

#: The equations the registry feeds, by the name every entry gives its
#: contribution to each. :class:`ForceTerm` gives none of them a default, so a
#: force added to the moment sums and forgotten in the horizontal balance is a
#: TypeError at import rather than a silent omission on the page.
CONSUMERS = ("moment_res", "moment_drv", "force_res", "force_drv",
             "spencer_h", "spencer_v", "march_x", "march_y", "normal")


@dataclass(frozen=True)
class NotApplicable:
    """Why a force contributes nothing to one of the equations.

    Not a blank: a reason, in the mechanics of the equation it is absent from.

    ``published`` is the terms the equation carries AS PUBLISHED for a force the
    solver does not carry — the shear on the top of the slice, which is in
    Spencer's equations (1) and (2) and in no xslope model. A section that
    transcribes the published equation before reducing it to this model needs
    both: the term, so the transcription is complete, and the reason, so the
    reduction can state why it went. Empty for every absence that is a statement
    about the equation itself, where the published form carries no such term.
    """

    reason: str
    published: tuple = ()


@dataclass(frozen=True)
class Term:
    """One printed term of one equation.

    ``values`` is called with the calculation context (:class:`_Calc`) and
    returns the term's value on every slice — what decides whether the model
    exercises it, and what the sum is formed from.

    ``rank`` places the term inside its own equation where that is not the order
    of the registry. Each equation is printed in the order its own derivation
    writes it, and they do not agree: the moment sum opens with the weight, the
    horizontal balance carries the surface loads after the body forces, and the
    base-normal equation ends with the two terms that come off the base.

    ``always`` prints the term whether or not this model exercises it. Only the
    base-normal equation has such terms — it is a statement of vertical
    equilibrium rather than a sum over the forces this model applies, and the
    weight and the mobilized cohesion are in it however small they are.
    """

    sign: int
    symbol: str
    values: object
    rank: float = None
    always: bool = False


@dataclass(frozen=True)
class Symbol:
    """A symbol one force introduces into the equations, and what it means.

    ``group`` is where the nomenclature lists it: the angles, then the moment
    arms and the quantities that live in the equation alone, then the letters
    that rename a slice-table column, then the letters Spencer's page takes from
    UTEXAS. Within a group the order is the registry's — which is the order the
    equations carry the forces in — unless ``rank`` says otherwise, as it does
    for the UTEXAS letters: those are listed in the order Spencer's own symbol
    list introduces them.

    Symbols that are slice-table COLUMNS are not declared here. Their definitions
    come from the column registry (:mod:`xslope.columns`), which is what the
    table's legend is written from, so the equation and the table cannot describe
    the same quantity two different ways.

    A symbol in a group of :data:`SECTION_SYMBOL_GROUPS` is defined for one
    section rather than for the report: T is the shear on the top of the slice on
    Spencer's page and the tension-crack water force on every other, and a letter
    that means two things cannot go in one nomenclature.
    """

    name: str
    group: str
    meaning: str
    rank: float = None


#: The nomenclature's groups, in the order it prints them.
SYMBOL_GROUPS = ("angle", "arm", "letter", "spencer")

#: Groups that do NOT go into :data:`EQUATION_SYMBOLS`: a letter one derivation
#: gives a force the others give another, which the section that prints it hands
#: to :func:`equation_symbols` as an override.
SECTION_SYMBOL_GROUPS = ("spencer_only",)


@dataclass(frozen=True)
class ForceTerm:
    """One force a slice can carry, and its part in every printed equation.

    ``columns`` are the slice-table columns the force's magnitude arrives in:
    what makes it present in a model, what the omission sentence is decided by,
    and — for the passive entries — what the passive gate is built from.
    ``arrays`` are the named per-slice arrays its contributions read, as
    ``(name, column, degrees)``; ``degrees`` converts on the way in.

    ``feature`` is what the omission sentence calls the force. The empty string
    is a force every model has, which is never reported absent. Active and
    passive support share a feature name: a model reinforced entirely by passive
    capacity prints its P_p terms, and testing the tangent column alone denied
    the reinforcement in the equation directly above the sentence.

    ``passive`` marks capacity that mobilizes with the soil and so carries 1/F.
    """

    key: str
    columns: tuple
    arrays: tuple
    symbols: tuple
    feature: str
    passive: bool
    moment_res: object
    moment_drv: object
    force_res: object
    force_drv: object
    spencer_h: object
    spencer_v: object
    march_x: object
    march_y: object
    normal: object


#: Excluded from Spencer's F_h and F_v by the definition of equations (1) and
#: (2): they are the forces on the slice OTHER than the base normal, the base
#: shear and the interslice forces, which is what makes Q solvable from them.
_NOT_IN_SPENCER_SUMS = NotApplicable(
    "equations (1) and (2) sum the forces on the slice other than the base "
    "normal, the base shear and the interslice forces")

#: A passive force divides by the factor of safety, so it stands on both sides
#: of every equation it enters. The moment methods can still show it — it makes
#: a resisting moment of its own — but no quotient and no force sum can.
_PASSIVE_CARRIES_F = NotApplicable(
    "passive capacity mobilizes with the soil and enters divided by F, which "
    "the compact form cannot show")

#: Horizontal forces have no vertical component to give the base-normal equation
#: or Spencer's F_v, and vertical forces none to give the horizontal balance.
_NO_VERTICAL_COMPONENT = NotApplicable("the force is horizontal")
_NO_HORIZONTAL_COMPONENT = NotApplicable("the force is vertical")

#: A force xslope's solver does not simulate, in an equation that publishes no
#: term for it either.
_NOT_SIMULATED = NotApplicable("xslope does not simulate it")


def _unexercised(C):
    """A published term no model carries: zero on every slice, every time.

    The transcription of a published equation prints it; nothing sums it. Written
    as a per-slice array like every other contribution so that the registry's own
    check can evaluate it alongside them.
    """
    return C.A["W"] * 0.0

#: Every force of the general equations, in the order the equations carry them:
#: what comes off the base first, then the body forces, then the forces applied
#: to the slice, then the support that mobilizes with the soil.
FORCE_TERMS = (
    ForceTerm(
        key="strength",
        columns=("c", "phi", "dl", "n_eff"),
        arrays=(("dl", "dl", False), ("N", "n_eff", False)),
        symbols=(Symbol("a_S", "arm",
                        "moment arm of the base shear about the center of "
                        "rotation"),),
        feature="", passive=False,
        moment_res=(Term(+1, "(c·Δl + N'·tan φ)·a_S",
                         lambda C: (C.A["c"] * C.A["dl"]
                                    + C.A["N"] * C.A["tan_phi"])
                         * C.arms["a_S"]),),
        moment_drv=NotApplicable("the mobilized strength resists"),
        force_res=(Term(+1, "(c·Δl + N'·tan φ)·cos α",
                        lambda C: (C.A["c"] * C.A["dl"]
                                   + C.A["N"] * C.A["tan_phi"])
                        * C.A["cos_a"]),),
        force_drv=NotApplicable("the mobilized strength resists"),
        spencer_h=_NOT_IN_SPENCER_SUMS,
        spencer_v=_NOT_IN_SPENCER_SUMS,
        # The march writes the strength mobilized on the base, c_m = c/F and
        # tan φ_m = tan φ/F: the cohesion on the right of each equation and the
        # friction on the left, with the base normal it multiplies. Both are in
        # the equilibrium of every slice however small this model's cohesion is.
        march_x=(Term(-1, "c_m·Δl·cos α",
                      lambda C: C.A["c"] * C.A["dl"] * C.A["cos_a"],
                      rank=0, always=True),),
        march_y=(Term(-1, "c_m·Δl·sin α",
                      lambda C: C.A["c"] * C.A["dl"] * C.A["sin_a"],
                      rank=0, always=True),),
        # The mobilized cohesion is the term that puts F in the base-normal
        # equation and makes the solution iterative, and it is printed whether
        # or not this model's soils have any.
        normal=(Term(-1, "frac{c·Δl·sin α}{F}",
                     lambda C: C.A["c"] * C.A["dl"] * C.A["sin_a"],
                     rank=99, always=True),),
    ),
    ForceTerm(
        key="normal",
        columns=("n_eff", "u", "dl"),
        arrays=(("u", "u", False),),
        symbols=(Symbol("a_N", "arm",
                        "moment arm of the total base normal force about the "
                        "center of rotation"),),
        feature="", passive=False,
        moment_res=NotApplicable("the total base normal force drives"),
        moment_drv=(Term(-1, "(N' + u·Δl)·a_N",
                         lambda C: (C.A["N"] + C.A["u"] * C.A["dl"])
                         * C.arms["a_N"]),),
        force_res=NotApplicable("the total base normal force drives"),
        force_drv=(Term(+1, "(N' + u·Δl)·sin α",
                        lambda C: (C.A["N"] + C.A["u"] * C.A["dl"])
                        * C.A["sin_a"]),),
        spencer_h=_NOT_IN_SPENCER_SUMS,
        spencer_v=_NOT_IN_SPENCER_SUMS,
        # The pore water on the base is known ahead of the march, so it is on
        # the right of each equation with the other knowns; the effective
        # normal, which is not, is on the left.
        march_x=(Term(+1, "u·Δl·sin α",
                      lambda C: C.A["u"] * C.A["dl"] * C.A["sin_a"],
                      rank=4, always=True),),
        march_y=(Term(-1, "u·Δl·cos α",
                      lambda C: C.A["u"] * C.A["dl"] * C.A["cos_a"],
                      rank=4, always=True),),
        # The water on the base is the one part of the total normal force that
        # is known ahead of the solution, so it is on the numerator's side.
        normal=(Term(-1, "u·Δl·cos α",
                     lambda C: C.A["u"] * C.A["dl"] * C.A["cos_a"], rank=98),),
    ),
    ForceTerm(
        key="W",
        columns=("w",),
        arrays=(("W", "w", False),),
        symbols=(Symbol("x_r", "arm",
                        "horizontal moment arm of the slice weight about the "
                        "center of rotation"),),
        feature="", passive=False,
        moment_res=NotApplicable("the weight drives"),
        # Equation (8a) opens with the weight, ahead of the base terms.
        moment_drv=(Term(+1, "W·x_r", lambda C: C.A["W"] * C.arms["x_r"],
                         rank=-1),),
        force_res=_NO_HORIZONTAL_COMPONENT,
        force_drv=_NO_HORIZONTAL_COMPONENT,
        spencer_h=_NO_HORIZONTAL_COMPONENT,
        spencer_v=(Term(-1, "W", lambda C: C.A["W"]),),
        march_x=_NO_HORIZONTAL_COMPONENT,
        march_y=(Term(+1, "W", lambda C: C.A["W"], rank=5),),
        normal=(Term(+1, "W", lambda C: C.A["W"], always=True),),
    ),
    ForceTerm(
        key="D",
        columns=("dload",),
        arrays=(("D", "dload", False), ("beta", "beta", True),
                ("d_x", "d_x", False), ("d_y", "d_y", False)),
        symbols=(Symbol("β", "angle",
                        "inclination of the distributed load from vertical "
                        "(perpendicular to the slope)"),
                 # The transcribed equations carry the surface load whether or
                 # not this model applies one, and a model that applies none has
                 # no column for the letter to be defined from.
                 Symbol("D", "letter",
                        "resultant of the distributed load acting on the top "
                        "of the slice, per unit thickness"),
                 Symbol("a_dx", "arm",
                        "horizontal moment arm of the distributed load"),
                 Symbol("a_dy", "arm",
                        "vertical moment arm of the distributed load")),
        feature="distributed load", passive=False,
        moment_res=NotApplicable(
            "a surface load is carried on the driving side, where its two "
            "components enter with their own signs"),
        moment_drv=(Term(+1, "D·cos β·a_dx",
                         lambda C: C.A["D"] * C.cos("beta") * C.arms["a_dx"]),
                    Term(-1, "D·sin β·a_dy",
                         lambda C: C.A["D"] * C.sin("beta") * C.arms["a_dy"])),
        force_res=NotApplicable(
            "a surface load is carried on the driving side, where the sign of "
            "its horizontal component decides which way it acts"),
        # The horizontal balance and Spencer's equation (1) both write the
        # surface loads after the body forces.
        force_drv=(Term(-1, "D·sin β",
                        lambda C: C.A["D"] * C.sin("beta"), rank=5.5),),
        spencer_h=(Term(+1, "P sin β",
                        lambda C: C.A["D"] * C.sin("beta"), rank=5.5),),
        spencer_v=(Term(-1, "P cos β", lambda C: C.A["D"] * C.cos("beta")),),
        march_x=(Term(-1, "D sin β",
                      lambda C: C.A["D"] * C.sin("beta"), rank=6),),
        march_y=(Term(+1, "D cos β",
                      lambda C: C.A["D"] * C.cos("beta"), rank=7),),
        normal=(Term(+1, "D cos β", lambda C: C.A["D"] * C.cos("beta")),),
    ),
    ForceTerm(
        key="kW",
        columns=("kw",),
        arrays=(("kW", "kw", False), ("y_cg", "y_cg", False)),
        symbols=(Symbol("a_s", "arm",
                        "moment arm of the seismic force, taken at the slice "
                        "center of gravity"),
                 # Spencer's section transcribes the published equations before
                 # reducing them, so kW is printed on a model with no seismic
                 # load, where the column it would be read from is not.
                 Symbol("kW", "letter",
                        "horizontal seismic force on the slice, per unit "
                        "thickness — column kW of the slice table")),
        feature="seismic load", passive=False,
        moment_res=NotApplicable("the seismic force drives"),
        moment_drv=(Term(+1, "kW·a_s", lambda C: C.A["kW"] * C.arms["a_s"]),),
        force_res=NotApplicable("the seismic force drives"),
        force_drv=(Term(+1, "kW", lambda C: C.A["kW"]),),
        spencer_h=(Term(-1, "kW", lambda C: C.A["kW"]),),
        spencer_v=_NO_VERTICAL_COMPONENT,
        march_x=(Term(+1, "kW", lambda C: C.A["kW"], rank=7),),
        march_y=_NO_VERTICAL_COMPONENT,
        normal=_NO_VERTICAL_COMPONENT,
    ),
    ForceTerm(
        key="T",
        columns=("t",),
        arrays=(("T", "t", False), ("y_t", "y_t", False)),
        symbols=(Symbol("a_t", "arm",
                        "moment arm of the tension-crack water force"),
                 Symbol("T", "letter",
                        "resultant force of the water in a tension crack — "
                        "column T_c of the slice table"),
                 # "written T below" while T stands two rows above it as the
                 # shear on the top of the slice is a nomenclature contradicting
                 # itself; the other derivations are named instead.
                 Symbol("V", "spencer",
                        "resultant force of the water in a tension crack — "
                        "column T_c of the slice table, which the other "
                        "methods' equations write T", rank=2)),
        feature="tension-crack water force", passive=False,
        moment_res=NotApplicable("the water in the crack drives"),
        moment_drv=(Term(+1, "T·a_t", lambda C: C.A["T"] * C.arms["a_t"]),),
        force_res=NotApplicable("the water in the crack drives"),
        force_drv=(Term(+1, "T", lambda C: C.A["T"]),),
        spencer_h=(Term(-1, "V", lambda C: C.A["T"]),),
        spencer_v=_NO_VERTICAL_COMPONENT,
        march_x=(Term(+1, "T", lambda C: C.A["T"], rank=8),),
        march_y=_NO_VERTICAL_COMPONENT,
        normal=_NO_VERTICAL_COMPONENT,
    ),
    ForceTerm(
        key="P",
        # Tangent reinforcement arrives in P and axial in the pa_ columns; the
        # equations write both with the one letter, at the angle ψ.
        columns=("p", "pa_cx", "pa_cy"),
        arrays=(("P", "p", False), ("pa_cx", "pa_cx", False),
                ("pa_cy", "pa_cy", False), ("pa_mx", "pa_mx", False),
                ("pa_my", "pa_my", False)),
        symbols=(Symbol("ψ", "angle",
                        "angle of the reinforcement force from horizontal"),
                 Symbol("a_rx", "arm",
                        "horizontal moment arm of the reinforcement force"),
                 Symbol("a_ry", "arm",
                        "vertical moment arm of the reinforcement force"),
                 # On an axially reinforced model the tangent column P is zero
                 # on every slice and is not printed, so the letter is defined
                 # here rather than left to the column registry.
                 Symbol("P", "letter",
                        "reinforcement force crossing the failure surface, at "
                        "the angle ψ from horizontal"),
                 Symbol("R", "spencer",
                        "reinforcement force on the slice base, at the angle ψ "
                        "from horizontal — column P of the slice table, and "
                        "written P below", rank=1)),
        feature="reinforcement crossing the failure surface", passive=False,
        moment_res=NotApplicable(
            "active reinforcement is carried on the driving side, where its "
            "negative sign makes it resist"),
        moment_drv=(Term(-1, "P·a_S", lambda C: C.A["P"] * C.arms["a_S"]),
                    Term(-1, "(P cos ψ·a_ry + P sin ψ·a_rx)",
                         lambda C: C.A["pa_cx"] * C.arms["Yo"] - C.A["pa_my"]
                         + C.A["pa_mx"] - C.arms["Xo"] * C.A["pa_cy"])),
        force_res=NotApplicable(
            "active reinforcement is carried on the driving side, where its "
            "negative sign makes it resist"),
        force_drv=(Term(-1, "P cos ψ",
                        lambda C: C.A["P"] * C.A["cos_a"] + C.A["pa_cx"]),),
        spencer_h=(Term(+1, "R cos ψ",
                        lambda C: C.A["P"] * C.A["cos_a"] + C.A["pa_cx"]),),
        spencer_v=(Term(+1, "R sin ψ",
                        lambda C: C.A["P"] * C.A["sin_a"] + C.A["pa_cy"]),),
        march_x=(Term(-1, "P cos ψ",
                      lambda C: C.A["P"] * C.A["cos_a"] + C.A["pa_cx"],
                      rank=1),),
        march_y=(Term(-1, "P sin ψ",
                      lambda C: C.A["P"] * C.A["sin_a"] + C.A["pa_cy"],
                      rank=1),),
        normal=(Term(-1, "P sin ψ",
                     lambda C: C.A["P"] * C.A["sin_a"] + C.A["pa_cy"]),),
    ),
    ForceTerm(
        key="H",
        columns=("h_pile",),
        arrays=(("H", "h_pile", False), ("theta_p", "theta_p", False),
                ("x_pile", "x_pile", False), ("y_pile", "y_pile", False)),
        symbols=(Symbol("θ_p", "angle",
                        "angle of the pile force from horizontal "
                        "(positive counterclockwise, upward)"),
                 Symbol("a_ex", "arm",
                        "horizontal moment arm of the pile force"),
                 Symbol("a_ey", "arm",
                        "vertical moment arm of the pile force"),
                 Symbol("H", "spencer",
                        "pile or pier force on the slice at the failure "
                        "surface — column H_p of the slice table", rank=3)),
        feature="pile force", passive=False,
        moment_res=NotApplicable(
            "an active pile force is carried on the driving side, where its "
            "negative sign makes it resist"),
        # The column carries the passive share too, and that share is on the
        # resisting side under its own entry.
        moment_drv=(Term(-1, "(H cos θ_p·a_ey + H sin θ_p·a_ex)",
                         lambda C: (C.A["H"] - C.A["H_p"]) * C.cos("theta_p")
                         * C.arms["a_ey"]
                         + (C.A["H"] - C.A["H_p"]) * C.sin("theta_p")
                         * C.arms["a_ex"]),),
        force_res=NotApplicable(
            "an active pile force is carried on the driving side, where its "
            "negative sign makes it resist"),
        force_drv=(Term(-1, "H cos θ_p",
                        lambda C: C.A["H"] * C.cos("theta_p")),),
        spencer_h=(Term(+1, "H cos θ_p",
                        lambda C: C.A["H"] * C.cos("theta_p")),),
        spencer_v=(Term(+1, "H sin θ_p",
                        lambda C: C.A["H"] * C.sin("theta_p")),),
        march_x=(Term(-1, "H cos θ_p",
                      lambda C: C.A["H"] * C.cos("theta_p"), rank=2),),
        march_y=(Term(-1, "H sin θ_p",
                      lambda C: C.A["H"] * C.sin("theta_p"), rank=2),),
        normal=(Term(-1, "H sin θ_p",
                     lambda C: C.A["H"] * C.sin("theta_p")),),
    ),
    ForceTerm(
        key="L",
        columns=("lload",),
        arrays=(("L", "lload", False), ("ll_b", "ll_beta", True),
                ("ll_x", "ll_x", False), ("ll_y", "ll_y", False)),
        symbols=(Symbol("δ", "angle",
                        "angle of the line load from horizontal "
                        "(−90° is straight down)"),
                 Symbol("a_fx", "arm",
                        "horizontal moment arm of the line load"),
                 Symbol("a_fy", "arm", "vertical moment arm of the line load"),
                 Symbol("L", "arm",
                        "line load applied on top of the slice, per unit "
                        "thickness")),
        feature="line load", passive=False,
        moment_res=NotApplicable(
            "a surface load is carried on the driving side, where its two "
            "components enter with their own signs"),
        moment_drv=(Term(+1, "L·a_fx",
                         lambda C: C.A["L"] * C.cos("ll_b") * C.arms["a_fx"]),
                    Term(-1, "L·a_fy",
                         lambda C: C.A["L"] * C.sin("ll_b") * C.arms["a_fy"])),
        force_res=NotApplicable(
            "a surface load is carried on the driving side, where the sign of "
            "its horizontal component decides which way it acts"),
        force_drv=(Term(-1, "L cos δ",
                        lambda C: C.A["L"] * C.sin("ll_b")),),
        spencer_h=(Term(+1, "L cos δ", lambda C: C.A["L"] * C.sin("ll_b")),),
        spencer_v=(Term(+1, "L sin δ", lambda C: C.A["L"] * C.cos("ll_b")),),
        march_x=(Term(-1, "L cos δ",
                      lambda C: C.A["L"] * C.sin("ll_b"), rank=3),),
        march_y=(Term(-1, "L sin δ",
                      lambda C: C.A["L"] * C.cos("ll_b"), rank=3),),
        normal=(Term(-1, "L sin δ", lambda C: C.A["L"] * C.cos("ll_b")),),
    ),
    ForceTerm(
        # In the published equations (1) and (2) and in no xslope model: the
        # solver does not simulate a shear force on the top of the slice, and a
        # force no model can carry is not an omission the report can report.
        key="T_top",
        columns=(), arrays=(),
        # Spencer's page writes the tension-crack water force V and the shear on
        # the top of the slice T; every other page writes that same crack force
        # T. So this letter is defined for Spencer's section alone.
        symbols=(Symbol("T", "spencer_only",
                        "shear force on the top of the slice"),),
        feature="", passive=False,
        moment_res=NotApplicable("xslope does not simulate it"),
        moment_drv=NotApplicable("xslope does not simulate it"),
        force_res=NotApplicable("xslope does not simulate it"),
        force_drv=NotApplicable("xslope does not simulate it"),
        # In the published equations and in no model: the terms are transcribed,
        # and the reduction below them says they are not simulated. The ranks put
        # each one where its own equation writes it — after the distributed load
        # in (1), and after it again in (2), where that load is ranked earlier.
        spencer_h=NotApplicable(
            "the shear on the top of the slice is in equation (1) as published "
            "and is not simulated",
            published=(Term(+1, "T cos β", _unexercised, rank=5.75),)),
        spencer_v=NotApplicable(
            "the shear on the top of the slice is in equation (2) as published "
            "and is not simulated",
            published=(Term(+1, "T sin β", _unexercised, rank=3.5),)),
        march_x=_NOT_SIMULATED,
        march_y=_NOT_SIMULATED,
        normal=_NOT_SIMULATED,
    ),
    ForceTerm(
        key="P_p",
        columns=("p_pt", "pp_cx", "pp_cy", "pp_mx", "pp_my"),
        arrays=(("P_pt", "p_pt", False), ("pp_cx", "pp_cx", False),
                ("pp_cy", "pp_cy", False), ("pp_mx", "pp_mx", False),
                ("pp_my", "pp_my", False)),
        symbols=(Symbol("P_p", "letter",
                        "passive reinforcement force, which mobilizes with the "
                        "soil and so carries 1/F"),),
        feature="reinforcement crossing the failure surface", passive=True,
        moment_res=(Term(+1, "P_p·a_S",
                         lambda C: C.A["P_pt"] * C.arms["a_S"]),
                    Term(+1, "(P_p cos ψ·a_ry + P_p sin ψ·a_rx)",
                         lambda C: C.A["pp_cx"] * C.arms["Yo"] - C.A["pp_my"]
                         + C.A["pp_mx"] - C.arms["Xo"] * C.A["pp_cy"])),
        moment_drv=NotApplicable("passive capacity resists"),
        force_res=_PASSIVE_CARRIES_F,
        force_drv=_PASSIVE_CARRIES_F,
        spencer_h=_PASSIVE_CARRIES_F,
        spencer_v=_PASSIVE_CARRIES_F,
        march_x=_PASSIVE_CARRIES_F,
        march_y=_PASSIVE_CARRIES_F,
        normal=_PASSIVE_CARRIES_F,
    ),
    ForceTerm(
        key="H_p",
        columns=("h_pile_pas",),
        arrays=(("H_p", "h_pile_pas", False),),
        symbols=(Symbol("H_p", "letter",
                        "passive pile force, which mobilizes with the soil and "
                        "so carries 1/F"),),
        feature="pile force", passive=True,
        moment_res=(Term(+1, "(H_p cos θ_p·a_ey + H_p sin θ_p·a_ex)",
                         lambda C: C.A["H_p"] * C.cos("theta_p")
                         * C.arms["a_ey"]
                         + C.A["H_p"] * C.sin("theta_p") * C.arms["a_ex"]),),
        moment_drv=NotApplicable("passive capacity resists"),
        force_res=_PASSIVE_CARRIES_F,
        force_drv=_PASSIVE_CARRIES_F,
        spencer_h=_PASSIVE_CARRIES_F,
        spencer_v=_PASSIVE_CARRIES_F,
        march_x=_PASSIVE_CARRIES_F,
        march_y=_PASSIVE_CARRIES_F,
        normal=_PASSIVE_CARRIES_F,
    ),
)


def _group_symbols(group):
    """The registry's symbols of one nomenclature group, in printed order."""
    rows = [(index if s.rank is None else s.rank, place, s)
            for index, term in enumerate(FORCE_TERMS)
            for place, s in enumerate(term.symbols) if s.group == group]
    return {s.name: s.meaning
            for _rank, _place, s in sorted(rows, key=lambda r: (r[0], r[1]))}


#: What every symbol in a printed equation means, in the words the derivation
#: pages use. The calculations section prints the ones its own equation carries,
#: and nothing else, so a reader never has to guess what a letter stands for and
#: never has to read past ten definitions to find the one they wanted.
#:
#: The forces' own symbols — their angles, their moment arms, and the letters
#: that rename a column — come from :data:`FORCE_TERMS`, so a force arrives
#: defined. What is written out here is what belongs to no one force: the factor
#: of safety, the geometry of the base, the interslice inclination, and the
#: quantities the iterations are judged by.
EQUATION_SYMBOLS = {
    "F": "factor of safety",
    "α": "inclination of the slice base from horizontal",
    **_group_symbols("angle"),
    "θ": "inclination of the interslice forces from horizontal",
    "θ_i": "inclination of the interslice resultant on the left side of the "
           "slice",
    "θ_{i+1}": "inclination of the interslice resultant on the right side of "
               "the slice",
    "λ": "scaling factor on the interslice force function f(x)",
    "f(x)": "the interslice force function; tan θ = λ·f(x)",
    "c_m": "cohesion mobilized on the slice base, c/F",
    "φ_m": "friction mobilized on the slice base, tan φ_m = tan φ/F",
    "Z_i": "interslice resultant on the left side of the slice, carried in from "
           "the slice before it",
    "Z_{i+1}": "interslice resultant on the right side of the slice, which the "
               "march solves for and carries into the next slice",
    "m_α": "the base-normal denominator, which carries the factor of safety and "
           "is what makes the solution iterative",
    "f_o": "Janbu's empirical correction factor for the neglected interslice shear",
    "F_corr": "factor of safety after Janbu's correction",
    **_group_symbols("arm"),
    **_group_symbols("letter"),
    **_group_symbols("spencer"),
    "Q": "resultant of the interslice forces on the slice, at the inclination θ "
         "— column Q_s of the slice table",
    "y_Q": "elevation at which the interslice resultant Q acts",
    "x_b": "horizontal coordinate of the slice base mid-point",
    "R_1": "force imbalance of the whole sliding mass at a trial (F, θ) — zero "
           "at the solution",
    "R_2": "moment imbalance of the whole sliding mass at a trial (F, θ) — zero "
           "at the solution",
    "Z_n": "interslice force left over at the far end of the march — zero at the "
           "solution",
    "M_o": "moment of the whole sliding mass about the coordinate origin",
}


def equation_symbols(notation, slice_df=None, unit_labels=None, overrides=None):
    """``[(symbol, meaning), ...]`` for the symbols a printed equation uses.

    The equation is scanned for every symbol the report can define, and only
    those are returned, in the order they are defined here — a nomenclature that
    listed every symbol xslope knows would be a glossary, and the reader is
    holding one equation.

    A symbol that is a column of the slice table is defined from the column
    registry and marked as such, so a reader who wants the number goes to the
    table rather than looking for it in the prose.

    ``overrides`` is a ``{symbol: meaning}`` mapping that wins over both, for a
    section whose equations are published in a notation of their own. Spencer's
    page follows UTEXAS, in which P is the distributed-load resultant and R the
    reinforcement force — the two letters every other derivation here writes D
    and P — and a row that says so is the only way the printed equations and the
    slice table can be read side by side.
    """
    from .columns import header, selected_columns

    text = str(notation or "")
    out = []
    seen = set()

    by_label = {}
    for c in (selected_columns(slice_df) if slice_df is not None else ()):
        by_label.setdefault(c.label, c)
    present = _present_symbols(
        text, list(overrides or ()) + list(by_label) + list(EQUATION_SYMBOLS))

    for symbol, meaning in (overrides or {}).items():
        if symbol in present and symbol not in seen:
            seen.add(symbol)
            out.append((symbol, meaning))

    # The slice table's own columns first: those are the ones a reader can look
    # up a value for. Only the columns the table PRINTS — one dropped for being
    # zero on every slice is not where anybody will find a number, and pointing
    # at it sends the reader looking for a column that is not there.
    for label, c in by_label.items():
        if label in present and label not in seen:
            seen.add(label)
            out.append((header(c, unit_labels), c.description.rstrip(".")))

    for symbol, meaning in EQUATION_SYMBOLS.items():
        if symbol in present and symbol not in seen:
            seen.add(symbol)
            out.append((symbol, meaning))
    return out


def _present_symbols(text, candidates):
    """Which of ``candidates`` an equation really uses — longest match wins.

    A plain substring test reads the residual R_1 as the reinforcement force R,
    and defines, under an equation about equilibrium, a force that equation never
    mentions. Each symbol claims its characters longest first, and a shorter one
    is present only where a longer has left something for it: R_1 takes both of
    its characters, so R is present only if it stands somewhere on its own.
    """
    masked = text
    present = set()
    for symbol in sorted(set(candidates), key=len, reverse=True):
        if symbol and symbol in masked:
            present.add(symbol)
            masked = masked.replace(symbol, "\x00" * len(symbol))
    return present

#: How closely the printed quotient must reproduce the solver's factor of safety
#: before the calculation is shown at all. This is the arithmetic identity, in
#: exact floats — not the printed-precision reproduction, which is looser and is
#: what the checks assert. A relative 1e-6 catches a missing term, which moves a
#: factor of safety in the third digit or worse.
#:
#: It is one of the two statements :func:`_closes` will take; the other is the
#: solver's own tolerance, because a report that demands more than the solver
#: was asked to deliver refuses converged solutions.
CALC_TOLERANCE = 1e-6

#: How much looser than the solver the report's own check is.
#:
#: Bishop's and Janbu's iterations stop when the STEP in F falls below their
#: tolerance, and the force-equilibrium march and Morgenstern-Price when the
#: root-finder's step does. A step of that size does not put the answer that
#: close to the exact root: the remaining distance is the step over the rate the
#: iteration contracts at, which on a slowly contracting surface is tens of times
#: larger. Bishop's solution on the deep-seated vp061a sits 15 tolerances from
#: its own fixed point, and its quotient was refused for it.
#:
#: A hundred covers that, and still refuses an error a tenth the size of the
#: smallest a missing term can make.
CALC_SAFETY_FACTOR = 100

#: Which solver's convergence a method's printed numbers come out of. Corps of
#: Engineers and Lowe & Karafiath are two conventions for the interslice
#: inclination, and both are marched by ``force_equilibrium``.
SOLVER_FUNCTIONS = {"corps": "force_equilibrium", "lowe": "force_equilibrium"}


def _solver_tolerance(method):
    """The tolerance the solver behind ``method`` converges to, read off its own
    signature.

    Read rather than restated, so the two cannot drift: if a solver is retuned,
    the gate that judges its answers moves with it.

    Zero for the Ordinary Method of Slices, which is closed form and does not
    iterate at all; its arithmetic identity is exact to the last bit.
    """
    import inspect

    from . import solve

    fn = getattr(solve, SOLVER_FUNCTIONS.get(method, method), None)
    if fn is None:
        return 0.0
    tol = inspect.signature(fn).parameters.get("tol")
    if tol is None or tol.default is inspect.Parameter.empty:
        return 0.0
    return float(tol.default)


def _closes(residual, scale, method):
    """Is a residual small enough for the calculation to be printed?

    Two statements, and either one will do.

    The first is relative: the residual is nothing against the magnitudes it
    cancels within. That is what catches a term this module does not carry.

    The second is absolute and is the solver's own. Spencer's Newton pair stops
    when both imbalances fall below 1e-4 in force and moment units; the iterative
    quotient methods stop at a step in F below 1e-6. A solution that converged
    exactly as far as it was asked to has to be printable, and a relative 1e-6
    was refusing several: Spencer's on gs2_46 and vp037, Janbu's on vp040, and
    Bishop's and Janbu's on vp061a.

    The absolute statement is also what a collapsed scale needs. On a surface
    where every slice's Q acts along a line through the coordinate origin, its
    moment about that origin vanishes on every slice and the sum of the
    magnitudes comes to a few parts in a billion of one force-length unit. A
    relative test there compares rounding with rounding — 5e-10 against 2e-9
    reads as a quarter — and refused a solution whose moment equation is
    satisfied identically.
    """
    return abs(float(residual)) <= max(
        CALC_TOLERANCE * abs(float(scale)),
        CALC_SAFETY_FACTOR * _solver_tolerance(method))

#: Methods that take moments about the center of rotation, and methods that
#: balance horizontal forces over the whole sliding mass. Every method xslope
#: offers is in one list or the other.
MOMENT_METHODS = ("oms", "bishop")
FORCE_METHODS = ("janbu", "corps", "lowe", "spencer", "mprice")

#: Methods whose printed equation is not the equation their own page publishes.
#:
#: The force-equilibrium page derives a slice-by-slice march — its equations (6)
#: and (10), each slice's interslice force from the one before it — and
#: Morgenstern-Price's page derives no quotient at all. What the section prints
#: for these three is the horizontal force balance of the whole sliding mass at
#: the converged solution: true BECAUSE the march closed, since the interslice
#: forces cancel in the sum over the slices, but not the equation the derivation
#: publishes. Citing the page for it was citing it for something it does not
#: contain.
#:
#: Janbu's page writes that balance directly, as its equation (7), so Janbu's is
#: the published equation and is cited as one.
WHOLE_MASS_BALANCE_METHODS = ("corps", "lowe", "mprice")

#: Support features whose forces mobilize with the soil and so carry 1/F. They
#: put the factor of safety on both sides of every term they touch, which the
#: compact form cannot show; a model that uses one gets no worked calculation
#: from the methods that close on a quotient.
PASSIVE_COLUMNS = tuple(column for term in FORCE_TERMS if term.passive
                        for column in term.columns)


def _column(df, name, n=None):
    """One column as a float array, or zeros when the model does not carry it."""
    import numpy as np
    if name in df.columns:
        return df[name].values.astype(float)
    return np.zeros(len(df) if n is None else n)


def _any(values):
    import numpy as np
    return bool(np.any(np.abs(np.asarray(values, dtype=float)) > 1e-12))


def _calc_state(bundle):
    """The per-slice state the reported factor of safety was computed from.

    For an ordinary run that is the slice table as solved. For a rapid drawdown
    run the reported factor of safety belongs to a drawn-down stage: the solver
    hands back that stage's strengths, pore pressures and base normals in the
    caller's table, but the drawn-down loads live in the stage-2 columns, so they
    are read from there. Returns ``(DataFrame, stage label or "")``.
    """
    df = bundle.get("slice_df")
    results = bundle.get("results") or {}
    if df is None or "stage1_FS" not in results:
        return df, ""
    stage2 = _num(results.get("stage2_FS"))
    stage3 = _num(results.get("stage3_FS"))
    if stage2 is None or stage3 is None:
        return df, ""
    stage = 2 if stage2 <= stage3 else 3
    df = df.copy()
    for base, stage2 in (("u", "u2"), ("dload", "dload2"), ("d_x", "d_x2"),
                         ("d_y", "d_y2"), ("beta", "beta2")):
        if stage2 in df.columns:
            df[base] = df[stage2].values
    return df, (f"Stage {stage} — the drawn-down section with "
                f"{'undrained' if stage == 2 else 'drained'} strengths")


def _calc_arrays(df):
    """Every per-slice quantity the factor of safety equation needs, by the
    symbol the documentation gives it.

    The base geometry and the mobilized strength are formed here, because they
    are read by every equation and by the registry's own arithmetic. Everything
    else is the array a force declares in :data:`FORCE_TERMS`, so a force added
    there arrives with its columns and needs nothing added here.
    """
    import numpy as np

    n = len(df)
    alpha = np.radians(df["alpha"].values.astype(float))
    c = df["c"].values.astype(float)
    if "c_suction" in df.columns:
        c = c + df["c_suction"].values.astype(float)
    A = {
        "n": n,
        "alpha": alpha, "sin_a": np.sin(alpha), "cos_a": np.cos(alpha),
        "tan_phi": np.tan(np.radians(df["phi"].values.astype(float))),
        "c": c,
    }
    for term in FORCE_TERMS:
        for name, column, degrees in term.arrays:
            values = _column(df, column, n)
            A[name] = np.radians(values) if degrees else values
    return A


class _Calc:
    """The state one model's terms are evaluated on.

    Handed to every :class:`Term`'s ``values``. ``A`` is the named arrays; the
    moment arms need the center of rotation and the facing, so they are formed
    only when a moment equation asks for them — the force methods have no center
    of rotation to form them about.
    """

    def __init__(self, df, A, right_facing=False):
        self.df = df
        self.A = A
        self.right_facing = right_facing
        self._arms = None

    def sin(self, name):
        import numpy as np
        return np.sin(self.A[name])

    def cos(self, name):
        import numpy as np
        return np.cos(self.A[name])

    @property
    def arms(self):
        if self._arms is None:
            self._arms = _moment_arms_table(self.df, self.A, self.right_facing)
        return self._arms


def _moment_arms_table(df, A, right_facing):
    """Every moment arm of the equations, about the center of rotation.

    The x-arms are mirrored on a right-facing slope, which is the frame the
    solver's own arms and the slice inclinations live in.
    """
    from .solve import _moment_arms

    Xo = float(df["xo"].iloc[0])
    Yo = float(df["yo"].iloc[0])
    mirror = -1.0 if right_facing else 1.0
    x_r, a_S, a_N = _moment_arms(df, Xo, Yo, A["alpha"], right_facing)
    return {
        "Xo": Xo, "Yo": Yo, "x_r": x_r, "a_S": a_S, "a_N": a_N,
        "a_dx": (A["d_x"] - Xo) * mirror, "a_dy": Yo - A["d_y"],
        "a_s": Yo - A["y_cg"], "a_t": Yo - A["y_t"],
        "a_ex": (A["x_pile"] - Xo) * mirror, "a_ey": Yo - A["y_pile"],
        "a_fx": (A["ll_x"] - Xo) * mirror, "a_fy": Yo - A["ll_y"],
    }


def _equation_terms(consumer, C):
    """The terms of one equation, in the order that equation is published in.

    Every entry of :data:`FORCE_TERMS` contributes terms or says why it does
    not, and the terms are gathered in the registry's order except where one
    carries a rank of its own.
    """
    return _ordered_terms(consumer, published=False)


def _published_terms(consumer):
    """The terms of one equation AS PUBLISHED — including the ones no xslope
    model carries.

    :func:`_equation_terms` gives the equation the solver evaluates; this gives
    the equation the derivation prints, which is the same table plus whatever a
    :class:`NotApplicable` records as published. The two are assembled from the
    one registry so that a section can transcribe the published form and then
    reduce it without the two forms drifting apart.
    """
    return _ordered_terms(consumer, published=True)


def _ordered_terms(consumer, published):
    """One equation's terms, in the order that equation is published in."""
    rows = []
    for index, term in enumerate(FORCE_TERMS):
        got = getattr(term, consumer)
        if isinstance(got, NotApplicable):
            if not published:
                continue
            got = got.published
        rows.extend((index if t.rank is None else t.rank, place, t)
                    for place, t in enumerate(got))
    return [t for _rank, _place, t in sorted(rows, key=lambda r: (r[0], r[1]))]


def _evaluate(terms, C):
    """``[(sign, symbol, values), ...]`` — one equation's terms on this model."""
    return [(t.sign, t.symbol, t.values(C)) for t in terms]


def _keep(terms, scale):
    """The terms this model actually exercises.

    A term that is zero on every slice describes something the model does not do,
    and a column of zeros in a submittal is a question rather than an answer. The
    test is against the size of the largest term, not against zero: on a true
    circular arc the base normal's moment is zero by construction and arrives as
    float noise a few parts in 1e17 of the driving moment, which is an absence
    however it rounds.
    """
    import numpy as np

    floor = 1e-9 * abs(scale)
    return [t for t in terms
            if float(np.max(np.abs(np.asarray(t[2], dtype=float)))) > floor]


def _features():
    """``[(plain name, columns), ...]`` — every force the omission sentence can
    name, with every column it can arrive in, in the order the sentence lists
    them.

    A feature is present if ANY of its columns carries a value. Reinforcement is
    tangent, axial or passive and piles are active or passive, so each of those
    is two registry entries under one name: a model reinforced entirely by
    passive capacity has zero in the tangent and axial columns while the printed
    equation carries its P_p terms, and testing one entry alone denied a force
    the equation directly above the sentence was printing.
    """
    out = {}
    for term in FORCE_TERMS:
        if term.feature:
            out[term.feature] = out.get(term.feature, ()) + term.columns
    return list(out.items())


def _absent_features(A, df):
    """The forces of the general equation this model does not carry at all."""
    return [name for name, columns in _features()
            if not any(_any(_column(df, col, A["n"])) for col in columns)]


def _moment_terms(df, A, right_facing):
    """``(resisting_terms, driving_terms)`` for a moment method, as equation (8a)
    of the Ordinary Method of Slices page writes them.

    Each list is ``[(sign, symbol, values), ...]`` and each moment is the signed
    sum of its terms' values. Passive support — capacity that mobilizes with the
    soil — contributes a resisting moment of its own, which is why the resisting
    side is a list rather than one term.
    """
    C = _Calc(df, A, right_facing)
    return (_evaluate(_equation_terms("moment_res", C), C),
            _evaluate(_equation_terms("moment_drv", C), C))


def _force_terms(A):
    """``(resisting, driving_terms)`` for a force-equilibrium method.

    Janbu's page writes the balance directly (its equation 7); the march the
    other four run is the force-equilibrium page's equations (6) and (10), and
    the same march under Spencer's and Morgenstern-Price's converged interslice
    forces. All five carry the horizontal component of the surface loads on the
    driving side: a load written as (D·sin β, −D·cos β) has a SIGNED horizontal
    component in the sliding-direction frame, and a load normal to a slope face
    (β > 0) pushes into the slope.
    """
    C = _Calc(None, A)
    return (_evaluate(_equation_terms("force_res", C), C),
            _evaluate(_equation_terms("force_drv", C), C))


def _signed_notation(kept):
    """A plain signed sum of the live terms — ``−W − P cos β``. The sibling of
    :func:`_sum_notation` for an equation whose terms are per-slice forces rather
    than sums over the slices."""
    out = ""
    for sign, symbol, _values in kept:
        if not out:
            out = symbol if sign > 0 else f"−{symbol}"
        else:
            out += f" + {symbol}" if sign > 0 else f" − {symbol}"
    return out or "0"


def _reduction(printed, consumers, absent):
    """The clause that takes a published equation down to this model.

    Read off the same registry the published form is assembled from, so that what
    the sentence says went is what went. Three things can take a term out of the
    published form, and each is stated as what it is: a force the model does not
    carry at all, named as a force; the shear on the top of the slice, which is
    in Spencer's published equations and in no xslope model; and a term of a
    force the model DOES carry that happens to be zero on every slice — a
    distributed load applied vertically has no horizontal component, and its
    P sin β goes while its P cos β stays.

    ``absent`` is what the model carries none of at all, by feature, so that the
    sentence names a FORCE where the force is missing and a TERM where the force
    is there and one of its components is not. A model with active reinforcement
    and no passive capacity carries reinforcement: its P_p terms are dropped as
    terms, not as a reinforcement the equation directly above prints.

    Empty when nothing went, which is when the published equation is already this
    model's and stands alone.
    """
    wholly, partly, top = [], [], False
    for term in FORCE_TERMS:
        published = []
        for consumer in consumers:
            got = getattr(term, consumer)
            published.extend(got.published if isinstance(got, NotApplicable)
                             else got)
        gone = [t.symbol for t in published if t.symbol not in printed]
        if not gone:
            continue
        if term.key == "T_top":
            top = True
        elif term.feature in absent:
            if term.feature not in wholly:
                wholly.append(term.feature)
        else:
            partly.extend(gone)

    parts = []
    if wholly:
        parts.append(f"the model carries no {_join(wholly)}")
    if top:
        parts.append("no shear on the top of the slice is simulated")
    if partly:
        parts.append(f"{_join(partly)} "
                     f"{'is' if len(partly) == 1 else 'are'} zero on every "
                     f"slice")
    if not parts:
        return ""
    lead = (parts[0] if len(parts) == 1
            else ", ".join(parts[:-1] + [f"and {parts[-1]}"]))
    return f"{lead[0].upper()}{lead[1:]}"


@dataclass(frozen=True)
class Transcription:
    """The published equations one method's section prints in full before it
    prints this model's, and the sentences that name them.

    ``consumers`` are the registry's contributions the full form is assembled
    from — the same ones the model's own equation is assembled from, so the two
    cannot drift apart or away from the page.

    ``build`` is the shape the derivation publishes the equation in: a quotient
    for F, a balance of two sides, or the pair of per-slice equations a march
    solves. ``lead`` introduces the transcription and names the equation the
    derivation numbers it; ``reduces`` finishes the sentence that follows it,
    after the clause saying what this model does not carry; ``solved`` is that
    sentence where nothing was dropped and the model's own equation still has to
    follow, which is so wherever the published form and the printed one are
    different shapes.

    ``link`` is a page other than the method's own that the lead names, as
    ``(phrase, page)``: Bishop's moment arms come from the Ordinary Method's
    derivation and Morgenstern-Price's march from the force-equilibrium page,
    and a number cited from another document has to say which.
    """

    consumers: tuple
    build: str
    lead: str
    reduces: str
    solved: str = ""
    link: tuple = ()


#: What each method transcribes, and from which of its derivation's equations.
#:
#: The two moment methods print the moment equilibrium their pages publish and
#: then that equation solved for F, because the solved form is a single quotient
#: whose driving side, written out for every force a slice can carry, is wider
#: than any page. Janbu's page publishes the solved form itself. The three
#: methods that march print the per-slice equilibrium the march solves: the
#: force-equilibrium page's equations (6) and (7), which is the derivation Corps
#: of Engineers and Lowe & Karafiath are documented on and the one
#: Morgenstern-Price's own page sends the reader to for its march. What those
#: three print BELOW their transcription is not a published equation at all — it
#: is the horizontal balance of the whole mass at the converged solution, and the
#: section says so (:data:`WHOLE_MASS_BALANCE_METHODS`).
TRANSCRIPTIONS = {
    "oms": Transcription(
        consumers=("moment_res", "moment_drv"), build="balance",
        lead="Equation (7) of the derivation takes moments about the center of "
             "rotation, with the strength mobilized on the base and the support "
             "that mobilizes with it on the left. It is written here in the "
             "general moment arms of equation (8a), which carry the base "
             "normal's own moment — zero on a true circle:",
        reduces="so equation (8a) — the same equilibrium solved for F — is:",
        solved="Equation (8a) is the same equilibrium solved for F:"),
    "bishop": Transcription(
        consumers=("moment_res", "moment_drv"), build="balance",
        lead="Equation (9) of the derivation takes moments about the center of "
             "rotation, with the strength mobilized on the base and the support "
             "that mobilizes with it on the left. It is written here in the "
             "general moment arms of equation (8a) of the Ordinary Method of "
             "Slices derivation, which carry the base normal's own moment — "
             "zero on a true circle:",
        reduces="so equation (9), solved for F, is:",
        solved="Solved for F, equation (9) is:",
        link=("Ordinary Method of Slices", "lem/oms.md")),
    "janbu": Transcription(
        consumers=("force_res", "force_drv"), build="quotient",
        lead="Equation (7) of the derivation balances the horizontal forces on "
             "the whole sliding mass, with the total base normal "
             "N = N' + u·Δl:",
        reduces="so equation (7) reduces to:"),
    "corps": Transcription(
        consumers=("march_x", "march_y"), build="march",
        lead="Equations (6) and (7) of the derivation are the horizontal and "
             "vertical equilibrium of one slice. The march solves them on each "
             "slice in turn for the base normal N' and the interslice resultant "
             "Z_{i+1} on its right, given the Z_i carried in from the slice "
             "before it:",
        reduces="so equations (6) and (7) reduce to:"),
    "lowe": Transcription(
        consumers=("march_x", "march_y"), build="march",
        lead="Equations (6) and (7) of the derivation are the horizontal and "
             "vertical equilibrium of one slice. The march solves them on each "
             "slice in turn for the base normal N' and the interslice resultant "
             "Z_{i+1} on its right, given the Z_i carried in from the slice "
             "before it:",
        reduces="so equations (6) and (7) reduce to:"),
    "mprice": Transcription(
        consumers=("march_x", "march_y"), build="march",
        lead="The march solves the same per-slice system as the "
             "force-equilibrium methods — equations (6) and (7) of that "
             "derivation, the horizontal and vertical equilibrium of one slice "
             "— for the base normal N' and the interslice resultant Z_{i+1} on "
             "its right, given the Z_i carried in from the slice before it:",
        reduces="so equations (6) and (7) reduce to:",
        link=("force-equilibrium methods", "lem/force_eq.md")),
    "spencer": Transcription(
        consumers=("spencer_h", "spencer_v"), build="spencer",
        lead="", reduces="so equations (1) and (2) reduce to:"),
}


#: The unknowns of the per-slice march, which are not forces the model applies
#: but what the march computes: the effective base normal and the interslice
#: resultant on the right of the slice, both on the left of each equation, and
#: the resultant carried in from the slice before, on the right among the knowns
#: at the rank the page writes it.
_MARCH_FRAME = {
    "march_x": ("N'·(tan φ_m·cos α − sin α) − Z_{i+1}·cos θ_{i+1}",
                Term(-1, "Z_i·cos θ_i", _unexercised, rank=5, always=True)),
    "march_y": ("N'·(tan φ_m·sin α + cos α) − Z_{i+1}·sin θ_{i+1}",
                Term(-1, "Z_i·sin θ_i", _unexercised, rank=6, always=True)),
}


def _march_equations(A):
    """``(full, reduced, printed)`` — equations (6) and (7) of the
    force-equilibrium derivation as published, this model's reduction of them,
    and the terms that reduction keeps.

    Both forms come off the one registry, as Spencer's two do: the published pair
    from every contribution the registry declares, the model's from the ones it
    exercises. The strength, the pore pressure and the interslice resultants are
    in the equilibrium of every slice whatever their size, so they are kept
    however small they come out.
    """
    C = _Calc(None, A)
    full, reduced, printed = [], [], []
    for consumer, (frame, carried) in _MARCH_FRAME.items():
        terms = sorted(list(_published_terms(consumer)) + [carried],
                       key=lambda t: 0 if t.rank is None else t.rank)
        kept = [t for t in terms if t.always or _any(t.values(C))]
        printed += [t.symbol for t in kept]
        full.append(f"{frame} = "
                    f"{_signed_notation([(t.sign, t.symbol, None) for t in terms])}")
        reduced.append(f"{frame} = "
                       f"{_signed_notation([(t.sign, t.symbol, None) for t in kept])}")
    return full, reduced, printed


def _moment_balance():
    """Moment equilibrium about the center of rotation — equation (7) of the
    Ordinary Method of Slices derivation and equation (9) of Bishop's — carrying
    every force a slice can take, in the general moment arms.

    The mobilized side is what divides by F: the strength on the base and the
    support that mobilizes with it. Everything else stands at its full value, on
    the side its sign puts it — which is why the driving terms the registry
    signs negative are written here as resisting moments, the way both pages
    write them.
    """
    left = [f"sum{{frac{{{t.symbol}}}{{F}}}}"
            for t in _published_terms("moment_res")]
    left += [f"sum{{{t.symbol}}}" for t in _published_terms("moment_drv")
             if t.sign < 0]
    right = [f"sum{{{t.symbol}}}" for t in _published_terms("moment_drv")
             if t.sign > 0]
    return f"{' + '.join(left)} = {' + '.join(right)}"


def _quotient_full(res_consumer, drv_consumer):
    """One method's factor of safety as its derivation publishes it: every force
    a slice can take, in the quotient the page solves for F."""
    res = [(t.sign, t.symbol, None) for t in _published_terms(res_consumer)]
    drv = [(t.sign, t.symbol, None) for t in _published_terms(drv_consumer)]
    return f"F = frac{{{_sum_notation(res)}}}{{{_sum_notation(drv)}}}"


def _transcription(method, A, printed, absent):
    """``(full, reduced, sentence)`` for one method's Calculations section.

    ``full`` is what the derivation publishes, ``reduced`` what is left of it on
    this model where the two are the same shape, and ``sentence`` says what went
    and what is below it. An empty sentence is a model that drops nothing: the
    published equation is this model's.
    """
    spec = TRANSCRIPTIONS[method]
    reduced = []
    if spec.build == "march":
        full, reduced, printed = _march_equations(A)
    elif spec.build == "balance":
        full = [_moment_balance()]
    else:
        full = [_quotient_full(*spec.consumers)]
    clause = _reduction(printed, spec.consumers, absent)
    if clause:
        sentence = f"{clause}, {spec.reduces}"
    else:
        sentence = spec.solved
    return full, reduced, sentence


def _spencer_force_sums(A, absent):
    """Equations (1) and (2) of the Spencer page, this model's reduction of them,
    and the symbols they need defining.

    The published equations are transcribed in full — every force the derivation
    carries, in the order that page writes them — and then reduced to the terms
    this model has, with a sentence between saying what went and why. Both forms
    are assembled from :data:`FORCE_TERMS`: the full one from
    :func:`_published_terms`, the reduced one from :func:`_equation_terms` with
    the zero terms dropped, so the transcription and the model's own arithmetic
    cannot drift from each other or from the page.

    That page follows UTEXAS, in which P is the distributed-load resultant, R the
    reinforcement force and V the tension-crack water force — three letters the
    other derivations here write D, P and T — and writes T for the shear on the
    top of the slice, which every other page writes for the crack force. The
    returned symbol definitions are what let the two be read side by side; they
    are handed to :func:`equation_symbols`, which prefers them to the column
    registry.

    Returns ``(blocks, symbols)``.
    """
    import numpy as np

    from .columns import BY_KEY

    C = _Calc(None, A)
    horizontal = _evaluate(_equation_terms("spencer_h", C), C)
    vertical = _evaluate(_equation_terms("spencer_v", C), C)
    scale = max(float(np.max(np.abs(np.asarray(v, dtype=float))))
                for _s, _sym, v in horizontal + vertical)

    kept_h = _keep(horizontal, scale)
    kept_v = _keep(vertical, scale)
    printed = [sym for _s, sym, _v in kept_h + kept_v]

    full = {name: [(t.sign, t.symbol, None) for t in _published_terms(name)]
            for name in ("spencer_h", "spencer_v")}
    blocks = [Math(f"F_h = {_signed_notation(full['spencer_h'])}"),
              Math(f"F_v = {_signed_notation(full['spencer_v'])}")]
    clause = _reduction(printed, ("spencer_h", "spencer_v"), absent)
    if clause:
        blocks.append(Prose(f"{clause}, {TRANSCRIPTIONS['spencer'].reduces}"))
        blocks.append(Math(f"F_h = {_signed_notation(kept_h)}"))
        blocks.append(Math(f"F_v = {_signed_notation(kept_v)}"))

    # T is the shear on the top of the slice here and the tension-crack water
    # force everywhere else, so its definition is the section's and not the
    # report's.
    symbols = dict(_group_symbols("spencer_only"))
    transcribed = [sym for _s, sym, _v in full["spencer_h"] + full["spencer_v"]]
    if any(sym.startswith("P ") for sym in transcribed + printed):
        # P is the one letter the section would otherwise use for two different
        # forces: the load on top of the slice here, the reinforcement below.
        dload, reinf = BY_KEY["dload"].label, BY_KEY["p"].label
        meaning = (f"resultant of the distributed load on the top of the slice, "
                   f"in the symbols of equations (1) and (2) — column {dload} of "
                   f"the slice table")
        if any(sym.startswith("R ") for sym in transcribed + printed):
            meaning += (f", and written {dload} below, where {reinf} is instead "
                        f"the reinforcement force")
        symbols["P"] = meaning
    return blocks, symbols


def _sum_notation(kept):
    """One side of the printed equation: one Σ per live term, signed."""
    out = ""
    for sign, symbol, _values in kept:
        joiner = "" if not out else (" + " if sign > 0 else " − ")
        if not out and sign < 0:
            joiner = "−"
        out += f"{joiner}sum{{{symbol}}}"
    return out or "0"


#: Equations (27) and (28) of the Spencer page — the two equilibrium equations
#: whose common root IS the solution. The page writes them for an assumed pair
#: (F_0, θ_0), because that is what a Newton step is taken from; the report
#: prints them at the pair that was found, so they are written in F and θ.
SPENCER_EQUILIBRIUM = ("R_1 = sum{Q}",
                       "R_2 = sum{Q·(x_b·sin θ − y_Q·cos θ)}")


def _spencer_state(df):
    """Spencer's per-slice ``F_h``, ``F_v``, ``Q`` and ``y_Q`` at the converged
    solution.

    Taken from the solver rather than rebuilt: its debug level writes ``F_h``,
    ``F_v``, ``Q``, ``y_q`` and the slice moment ``Mo`` into the table it is
    given, and the residuals the section prints are sums of exactly those
    numbers. The solve is repeated on a copy for that reason alone, and its
    console output — the whole point of a debug level — is swallowed.

    The force sums are lifted for the same reason the rest is: a second
    implementation of equations (1) and (2) here would be free to drift from the
    one the factor of safety came out of, which is the drift this whole section
    exists to make impossible. They are the solver's EFFECTIVE values — passive
    support at its mobilized 1/F share — so that equation (23) reproduces the Q
    beside them row by row.
    """
    import contextlib
    import io

    from .solve import spencer
    work = df.copy()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            ok, res = spencer(work, debug_level=2)
    except Exception:
        return None
    if not ok or not {"Q", "F_h", "F_v"} <= set(work.columns):
        return None
    import numpy as np
    Q = work["Q"].values.astype(float)
    y_q = work["y_q"].values.astype(float)
    theta = np.radians(res["theta"])
    x_b = work["x_c"].values.astype(float)
    return {
        "Q": Q, "y_q": y_q, "FS": res["FS"], "theta": res["theta"],
        "F_h": work["F_h"].values.astype(float),
        "F_v": work["F_v"].values.astype(float),
        # The solver's own facing test, not a second opinion on it: on a
        # right-facing surface every value above belongs to the mirrored slice,
        # and the section has to say so for a reader to check a row.
        "right_facing": bool(work.attrs.get(
            "right_facing", work["y_cb"].values[0] > work["y_cb"].values[-1])),
        # Equations (27) and (28): the force and moment imbalances the solver
        # drove to zero, recomputed from the same per-slice values the table
        # prints so that the section's numbers and its columns are one thing.
        "R1": float(np.sum(Q)),
        "R2": float(np.sum(Q * (x_b * np.sin(theta) - y_q * np.cos(theta)))),
        # What each residual is small COMPARED WITH. A number near zero says
        # nothing on its own; these are the sums of the magnitudes the
        # cancellation happens in, and they are what the section is gated on.
        "scale": float(np.sum(np.abs(Q))),
        "m_scale": float(np.sum(np.abs(
            Q * (x_b * np.sin(theta) - y_q * np.cos(theta))))),
    }


def _mp_residuals_for(df, results):
    """Morgenstern-Price's two equilibrium residuals at its converged solution,
    from the solver's own march — the force closure Z_n and the moment of the
    whole mass about the origin. ``(force, moment)``, or None."""
    from .solve import _mp_f_vals, _mp_march
    try:
        lam = float(results["lambda"])
        FS = float(results["FS"])
        f_vals = _mp_f_vals(df, results.get("f_type", "half_sine"))
        right_facing = bool(df.attrs.get(
            "right_facing", df["y_cb"].values[0] > df["y_cb"].values[-1]))
        _N, _Z, force, moment = _mp_march(df, lam, f_vals, FS, right_facing)
    except Exception:
        return None
    return float(force), float(moment)


#: What is true of a solution the equations cannot be worked through for. Each
#: is a fact about the analysis, printed where the working would have been, so
#: that a method block never goes quiet without saying why.
#:
#: Passive support is the common one. Its capacity mobilizes with the soil, so it
#: enters divided by the factor of safety and stands on both sides of the balance
#: the method solves; there is no quotient to evaluate and no pair of sums to
#: print. The moment methods CAN show it — it makes a resisting moment of its own
#: — so this is said only of the methods that close on a force balance.
PASSIVE_NOTE = (
    "The support in this model is passive: it mobilizes with the soil, so its "
    "force enters the equilibrium of every slice divided by the factor of "
    "safety and stands on both sides of the balance this method solves. There "
    "is no quotient of two sums behind this factor of safety to work through.")

NO_BASE_NORMAL_NOTE = (
    "The base normal force at the converged factor of safety is not among this "
    "solution's per-slice values, and every term of the equilibrium is formed "
    "from it, so no working is printed for it.")

NO_CENTER_NOTE = (
    "This surface has no center of rotation, and the equilibrium this method "
    "solves is taken in moments about one, so no working is printed for it.")

NO_QUOTIENT_NOTE = (
    "The driving terms of this model sum to zero over the slices, so the "
    "quotient this method closes on has no value to print.")

NO_INTERSLICE_NOTE = (
    "The interslice resultants this method's equations are written in did not "
    "re-form on this surface, so no working is printed for it.")

#: The range a number has to fall in to be quoted as a factor of safety.
#:
#: A slope at 0.01 carries a hundredth of the strength it needs and one at 100 a
#: hundred times, and neither is a slope: no model of a real section solves to
#: either, and no design decision turns on the difference between them and the
#: numbers further out. What lies outside is the arithmetic of a frame that was
#: never solved — an equation evaluated on initial values divides by a driving
#: sum near zero and comes out at 1e14 — and quoting that to six decimals dresses
#: a nonsense magnitude as a measurement. So the mismatch is stated instead, and
#: the only factor of safety printed is the one the solution reports.
CREDIBLE_FS = (0.01, 100.0)


def _mismatch_note(computed, FS):
    """What is true when the equation, evaluated on this solution's own values,
    does not return the solution.

    Both numbers are the statement when both are factors of safety — 1.752
    against 1.760 is a mismatch that can be sized by eye. A value outside
    :data:`CREDIBLE_FS` is not a factor of safety to compare, and what it is
    outside is stated in its place.
    """
    low, high = CREDIBLE_FS
    if low <= computed <= high:
        what = (f"gives a factor of safety of {computed:.6f}, and the solution "
                f"reports {FS:.6f}")
    else:
        what = (f"comes out outside the range a factor of safety can take "
                f"({low:g} to {high:g}), and the solution reports {FS:.6f}")
    return (f"Evaluated on this model's converged values the equation {what}. A "
            f"quotient that does not return the solution is not the working "
            f"behind it, and none is printed.")


def calculation(slope_data, bundle, method):
    """The factor of safety calculation the report prints, or None.

    Returns a dict carrying the slice table with the per-slice contribution
    columns added, the sums, the equation notation, and everything the section
    needs to write itself. None means the calculation could not be shown for this
    model — see :func:`_closes` and :data:`PASSIVE_COLUMNS`.
    """
    return _calculation(slope_data, bundle, method)[0]


def _calculation(slope_data, bundle, method):
    """``(calculation, note)`` — the working, or what is true instead.

    A solution whose equilibrium this module cannot work through is not a
    solution to go quiet over. Every refusal below carries the fact that
    produced it, and that fact is printed where the working would have been, so
    the two can never disagree: the note comes out of the same test that
    withheld the calculation.
    """
    import numpy as np

    method = str(method or "").lower()
    df, stage = _calc_state(bundle)
    results = bundle.get("results") or {}
    FS = _num(results.get("FS"))
    if df is None or not len(df) or FS is None or FS <= 0:
        return None, ""
    if method not in MOMENT_METHODS + FORCE_METHODS:
        return None, ""
    if "n_eff" not in df.columns or not np.all(np.isfinite(df["n_eff"].values)):
        return None, NO_BASE_NORMAL_NOTE
    if method in MOMENT_METHODS and "xo" not in df.columns:
        return None, NO_CENTER_NOTE

    A = _calc_arrays(df)
    right_facing = bool(df.attrs.get(
        "right_facing", df["y_lb"].iat[0] > df["y_rb"].iat[-1]))

    if method == "spencer":
        return _spencer_calculation(df, A, FS, stage)

    if method in MOMENT_METHODS:
        # Passive capacity contributes a resisting MOMENT here, so the moment
        # methods can show it. In the force methods it divides by F on the
        # driving side, which no quotient can print.
        res_terms, terms = _moment_terms(df, A, right_facing)
        res_key, drv_key = "m_res", "m_drv"
    else:
        if any(_any(_column(df, name)) for name in PASSIVE_COLUMNS):
            return None, PASSIVE_NOTE
        res_terms, terms = _force_terms(A)
        res_key, drv_key = "f_res", "f_drv"

    scale = max([float(np.max(np.abs(values)))
                 for _s, _sym, values in terms + res_terms])
    kept = _keep(terms, scale)
    kept_res = _keep(res_terms, scale)
    resisting = sum((sign * values for sign, _s, values in kept_res),
                    np.zeros(A["n"]))
    driving = sum((sign * values for sign, _s, values in kept),
                  np.zeros(A["n"]))
    sum_res = float(np.sum(resisting))
    sum_drv = float(np.sum(driving))
    if not np.isfinite(sum_res) or not np.isfinite(sum_drv) or sum_drv == 0:
        return None, NO_QUOTIENT_NOTE

    quotient = sum_res / sum_drv
    fo = _num(results.get("fo")) if method == "janbu" else None
    computed = quotient * fo if fo else quotient
    if not _closes(computed - FS, FS, method):
        return None, _mismatch_note(computed, FS)

    out = df.copy()
    out[res_key] = resisting
    out[drv_key] = driving

    absent = _absent_features(A, df)
    full, reduced, sentence = _transcription(
        method, A, [sym for _s, sym, _v in kept + kept_res], absent)

    return {
        "method": method, "slice_df": out, "FS": FS, "stage": stage,
        "resisting": sum_res, "driving": sum_drv, "quotient": quotient,
        "fo": fo, "spencer": None,
        "transcribed": full, "reduced": reduced, "reduction": sentence,
        "theta": _num(results.get("theta")) if method in ("corps", "lowe") else None,
        "lambda": _num(results.get("lambda")),
        "residuals": _mp_residuals_for(df, results) if method == "mprice" else None,
        "equation": (f"F = frac{{{_sum_notation(kept_res)}}}"
                     f"{{{_sum_notation(kept)}}}"),
        "kept": kept, "absent": absent,
        "res_key": res_key, "drv_key": drv_key,
        "normal_force": _normal_force_equations(A, method),
        "equilibrium": None, "force_sums": None, "symbols": {},
    }, ""


def _spencer_calculation(df, A, FS, stage):
    """Spencer's calculation, which does not end in a quotient.

    Every other method the report documents closes on a ratio of two sums, and
    that ratio IS the method: reproduce it and you have reproduced the factor of
    safety. Spencer's has no such form. F and θ are the pair at which the two
    equilibrium equations (27) and (28) both vanish, found by Newton's method,
    and the section says so rather than dividing two sums that happen to be
    lying around — which is a general horizontal balance, true of any solution,
    and no more Spencer's method than it is anybody else's.

    So this carries no ``F_R``/``F_D`` columns: those two exist to be divided.
    What takes the quotient's place, here and in the section's verification, is
    the pair of residuals, each of which has to vanish against the sum of the
    magnitudes it cancels within, or against the tolerance the Newton pair was
    driven to (:func:`_closes`).
    """
    state = _spencer_state(df)
    if state is None:
        return None, NO_INTERSLICE_NOTE
    # Passive support mobilizes with the soil and carries 1/F, which the printed
    # equations cannot show — the same exclusion the force methods make.
    if any(_any(_column(df, name)) for name in PASSIVE_COLUMNS):
        return None, PASSIVE_NOTE
    for residual, scale, what in (("R1", "scale", "force"),
                                  ("R2", "m_scale", "moment")):
        if not _closes(state[residual], state[scale], "spencer"):
            from .columns import format_residual, format_sum
            return None, (
                f"At the reported solution the {what} imbalance of the sliding "
                f"mass is {format_residual(state[residual])}, against the "
                f"{format_sum(state[scale])} the magnitudes it cancels within "
                f"come to. The solution is the pair at which it vanishes, and "
                f"no working is printed for one at which it does not.")

    absent = _absent_features(A, df)
    force_sums, force_symbols = _spencer_force_sums(A, absent)
    out = df.copy()
    out["F_h"] = state["F_h"]
    out["F_v"] = state["F_v"]
    out["q_s"] = state["Q"]
    out["y_q"] = state["y_q"]
    return {
        "method": "spencer", "slice_df": out, "FS": FS, "stage": stage,
        "spencer": state, "absent": absent,
        "equilibrium": list(SPENCER_EQUILIBRIUM),
        "force_sums": force_sums, "symbols": force_symbols,
        "transcribed": [], "reduced": [], "reduction": "",
        # No quotient, and so no columns to divide and no correction factor.
        "resisting": None, "driving": None, "quotient": None, "fo": None,
        "res_key": None, "drv_key": None, "equation": None, "kept": None,
        "theta": None, "lambda": None, "residuals": None, "normal_force": None,
    }, ""


def _normal_force_equations(A, method):
    """The base-normal equations the iterative methods solve alongside the
    quotient — Bishop's equation (8) / Janbu's equation (6) — carrying the
    vertical components this model actually has.

    Janbu's page names the denominator m_α and Bishop's writes it out; each is
    printed the way its own page writes it.
    """
    C = _Calc(None, A)
    kept = [(t.sign, t.symbol, None) for t in _equation_terms("normal", C)
            if t.always or _any(t.values(C))]
    numerator = _signed_notation(kept)
    if method == "janbu":
        return [f"N' = frac{{{numerator}}}{{m_α}}",
                "m_α = cos α + frac{sin α·tan φ}{F}"]
    return [f"N' = frac{{{numerator}}}{{cos α + frac{{sin α·tan φ}}{{F}}}}"]


#: How closely a sum rebuilt from the PRINTED per-slice values has to close, as
#: a fraction of the sum of their magnitudes. The solver's own residual is
#: vanishingly small; a sum re-formed from values rounded to a tenth of a force
#: unit carries that rounding, and this is the size it can carry. Stated in the
#: section itself, in words, because a reader who adds the column up gets this
#: number and not zero.
PRINTED_RESIDUAL_TOLERANCE = 1e-3


def _transcribed_blocks(calc, then=(), always=False):
    """The published equations, the sentence that reduces them to this model, and
    what is left of them — the shape every method's section prints.

    ``then`` is the model's own equation where the transcription is not already
    of that shape: the moment methods transcribe a balance and print that balance
    solved for F. ``always`` prints it whether or not anything was dropped, which
    is what those two need — a model that carries every force still has to be
    shown the quotient its arithmetic follows. Where the two forms are the same
    shape, a model that drops nothing has already been shown its own equation,
    and printing it twice would say the second one was different.
    """
    spec = TRANSCRIPTIONS[calc["method"]]
    url = docs_url(spec.link[1]) if spec.link else ""
    blocks = [Prose(spec.lead, links=[(spec.link[0], url)] if url else [])]
    blocks += [Math(line) for line in calc["transcribed"]]
    below = list(calc["reduced"]) or list(then)
    if calc["reduction"]:
        blocks.append(Prose(calc["reduction"]))
        blocks += [Math(line) for line in below]
    elif always:
        blocks.append(Prose(spec.solved))
        blocks += [Math(line) for line in below]
    return blocks


def _method_preamble(calc, method, figure_number=0):
    """What this method solves, ahead of the equation — one short block, in the
    method's own terms.

    ``figure_number`` is the free-body diagram the section opens on, cited by
    every method where the forces it draws are first spoken of. Zero where the
    section is built without a figure counter, and then not cited: a report that
    names Figure 0 is worse than one that names none.
    """
    from .columns import format_fs, format_residual, format_sum

    where, links = cite("Figure", figure_number)
    on_slice = f" of {where}" if where else ""
    blocks = []
    if method == "oms":
        # What the diagram is for, in the one method whose block would otherwise
        # not speak of it: the summary above states the method's assumptions,
        # and this states where the terms of the sums come from.
        blocks.append(Prose(
            f"Each slice's contribution to the two sums below is formed from "
            f"the forces on the slice{on_slice}.", links=links))
    elif method in ("bishop", "janbu"):
        # The base-normal equation is the derivation's own — Bishop's (8),
        # Janbu's (6) with the m_α of its (1) — carrying the terms this model
        # has, so the number it is named by is the number a reader looks it up
        # under.
        named = ("equation (8) of the derivation" if method == "bishop"
                 else "equation (6) of the derivation, whose denominator m_α is "
                      "its equation (1)")
        blocks.append(Prose(
            f"The base normal force N' comes from vertical equilibrium of the "
            f"slice{on_slice} and depends on the factor of safety itself, so it "
            f"and the quotient are solved together by iteration — {named}:",
            links=links))
        blocks.extend(Math(line) for line in calc["normal_force"])
        blocks.append(Prose(
            "Every N' below is that value at the converged factor of safety, so "
            "the sums formed from it return the F it was evaluated at."))
    elif method in ("corps", "lowe"):
        theta = calc.get("theta")
        stated = (f", averaging {theta:.2f} degrees on this surface"
                  if theta is not None else "")
        blocks.append(Prose(
            f"The interslice forces on the slice{on_slice} are taken at an "
            f"inclination θ that the method's own convention fixes from the "
            f"geometry{stated}. The solver marches the slices at a trial factor "
            f"of safety and adjusts it until the interslice force left over at "
            f"the far end is zero; at that value the interslice forces cancel "
            f"over the whole sliding mass, and the horizontal balance of the "
            f"mass is the quotient below.", links=links))
        blocks.extend(_transcribed_blocks(calc))
    elif method == "mprice":
        lam = calc.get("lambda")
        blocks.append(Prose(
            f"The interslice forces on the slice{on_slice} are inclined at a θ "
            f"that varies along the surface as tan θ = λ·f(x)"
            f"{f', with λ = {lam:.4f} at the solution' if lam is not None else ''}"
            f". λ and F are solved together so that force and moment "
            f"equilibrium of the whole sliding mass are satisfied at once. With "
            f"the interslice forces cancelling in the sum, the horizontal "
            f"balance is the quotient below.", links=links))
        blocks.extend(_transcribed_blocks(calc))
    elif method == "spencer":
        state = calc["spencer"]
        blocks.append(Prose(
            f"Spencer's method lumps the interslice forces on each slice into a "
            f"single resultant Q acting at the constant inclination θ, so that "
            f"force and moment equilibrium of the whole sliding mass are two "
            f"equations in the two unknowns F and θ. F_h and F_v are the sums "
            f"of the forces on the slice{on_slice} other than the base normal, "
            f"the base shear and the interslice forces, as equations (1) and "
            f"(2) of the derivation write them:", links=links))
        blocks.extend(calc["force_sums"])
        blocks.append(Prose(
            "Equations (23) and (24) give Q on each slice from them and from "
            "the strength mobilized on its base:"))
        blocks.append(Math(
            "Q = [−F_v·sin α − F_h·cos α − frac{c·Δl}{F} + "
            "(F_v·cos α − F_h·sin α + u·Δl)·frac{tan φ}{F}]·m_α"))
        blocks.append(Math(
            "m_α = frac{1}{cos (α − θ) + sin (α − θ)·frac{tan φ}{F}}"))
        blocks.append(Prose(
            "There is no equation for F: it appears inside Q, in the mobilized "
            "cohesion c·Δl/F and the friction tan φ/F, and again with θ in m_α. "
            "A trial pair (F, θ) fixes Q on every slice."))
        blocks.append(Prose(
            "F_h, F_v, Q and y_Q are columns F_h, F_v, Q_s and y_Q of the slice "
            "table; each row's Q follows from that row through the two "
            "equations above, at the converged F and θ."))
        if state["right_facing"]:
            # The moment equation is written for a slice of a left-facing
            # slope, so a right-facing section is solved as its mirror image.
            # Which quantities have to be mirrored to reproduce a row, and which
            # are ALREADY mirrored, is the whole content of this paragraph:
            # α, c and tan φ are printed as the slope has them, and the solver
            # negates them for the mirror; F_h and F_v are formed AFTER that
            # negation — the horizontal forces summed into F_h have already
            # changed sign — so they, and the Q and y_Q that follow from them,
            # go into equations (23) and (24) exactly as printed. Naming the
            # horizontal forces among the reversals had them reversed twice.
            blocks.append(Prose(
                "The derivation is written for a slice of a left-facing slope; "
                "this slope faces right and is solved as its mirror image. α, "
                "c and tan φ enter the equations above with their signs "
                "reversed. F_h, F_v, Q_s and y_Q are the mirrored slice's own "
                "values and enter as printed. The factor of safety, a ratio, "
                "is unaffected."))
    return blocks


#: The magnitudes a closure ratio is spoken in. A residual is one part in some
#: number of these; the word is what makes the size of that number readable.
_MAGNITUDES = ((1e12, "trillion"), (1e9, "billion"),
               (1e6, "million"), (1e3, "thousand"))

#: The small counts, spelled. "One part in a billion" and "one part in eight
#: billion" are sentences; "1e9" and "8e9" are not.
_SPELLED = {1: "a", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
            7: "seven", 8: "eight", 9: "nine", 10: "ten"}


def _closure_phrase(residual, scale):
    """How small ``residual`` is against ``scale``, as "one part in eight
    billion".

    A residual is only ever small compared with something. The count is one
    figure — the sentence is an order of magnitude, and the exact ratio is
    available to anyone who divides the two numbers printed beside it.
    """
    size, against = abs(float(residual)), abs(float(scale))
    if not size or not against:
        return ""
    ratio = against / size
    for unit, word in _MAGNITUDES:
        if ratio >= unit:
            count = round(ratio / unit)
            spelled = _SPELLED.get(count, f"{count:,}")
            return f"one part in {spelled} {word}"
    return f"one part in {round(ratio):,}"


def _spencer_close(calc, table_number, bookmark):
    """How Spencer's calculation ends: the pair the iteration reached, and each
    of the two equilibrium equations evaluated there.

    This is what the other methods' arithmetic is — the last line of the working,
    the one a reviewer checks the answer on. There it is a division; here it is
    the two equations balancing, because the solution was found by driving them
    to zero and not by evaluating a formula. So each is printed the way the
    quotient is: the equation, the number it comes out at, and the size of the
    sums that number cancels within — a residual against no scale is a digit a
    reader has no way to judge. Both are re-formable from the printed table, and
    the section says which columns.
    """
    from .columns import format_fs, format_residual, format_sum

    state = calc["spencer"]
    fs = format_fs(state["FS"])
    where = f"Table {table_number}" if table_number else "the slice table"
    link = [(where, f"#{bookmark}")] if table_number else []
    blocks = [Prose(
        f"The iteration converged at F = {fs} with θ = {state['theta']:.2f} "
        f"degrees. The force equation is the Q_s column of {where} added up, "
        f"printed as the total at the foot of that column:",
        bold=[fs], links=link)]
    blocks.append(Math(f"R_1 = sum{{Q}} = {format_residual(state['R1'])}"))
    blocks.append(Prose(
        f"That is {_closure_phrase(state['R1'], state['scale'])} of the "
        f"{format_sum(state['scale'])} the magnitudes of Q come to. The moment "
        f"equation weights the same column by x_c — the x_b of equation (28) — "
        f"and by y_Q, at the θ above:"))
    blocks.append(Math(
        f"R_2 = sum{{Q·(x_b·sin θ − y_Q·cos θ)}} = "
        f"{format_residual(state['R2'])}"))
    if _moment_arms_vanish(state):
        blocks.append(Prose(
            f"Q acts along a line through the coordinate origin on every slice "
            f"of this surface, so its moment about that origin vanishes term by "
            f"term; the terms come to {format_sum(state['m_scale'])} in total, "
            f"and the moment equation is satisfied identically rather than by "
            f"cancellation. Sums re-formed from the printed, rounded columns "
            f"carry that rounding and close to within a thousandth of the force "
            f"total above."))
    else:
        blocks.append(Prose(
            f"That is {_closure_phrase(state['R2'], state['m_scale'])} of the "
            f"{format_sum(state['m_scale'])} its moment terms come to. Sums "
            f"re-formed from the printed, rounded columns carry that rounding "
            f"and close to within a thousandth of these totals."))
    return blocks


def _moment_arms_vanish(state):
    """Does Q act along a line through the coordinate origin on every slice?

    The moment of Q about the origin is ``Q·(x_b·sin θ − y_Q·cos θ)``, which is
    zero exactly when its line of action passes through that origin. Where that
    happens on every slice the moment terms sum to nothing while the forces they
    are formed from do not, and the moment equation is satisfied identically —
    which is a different statement from a residual that cancels, and the reason
    the ratio of the two is not one to make.
    """
    floor = CALC_SAFETY_FACTOR * _solver_tolerance("spencer")
    return abs(state["m_scale"]) <= floor < abs(state["scale"])


def _calculations_section(calc, slope_data, table_number, unit_labels,
                          bookmark=SLICE_TABLE_BOOKMARK, counter=None):
    """The Calculations section: what the method solves, the equation, the sums,
    the arithmetic, the factor of safety.

    ``bookmark`` is the slice table this method's per-slice terms are columns of.
    A report that documents several methods carries a slice table for each, and a
    cross-reference has to land on the one whose numbers the sums came from.

    The blocks are built before any of them is placed, because the nomenclature
    in the middle of the section has to be drawn from ALL of them — including the
    arithmetic below it.
    """
    method = calc["method"]
    label = method_label(method)
    sec = Section("Calculations")

    # --- the slice the equations are written about ---
    #
    # The section heads on the free body the derivation draws, so every symbol
    # below meets its picture before it is used rather than a page or a browser
    # away from it. It is the first block of the section and the section is only
    # built where the working is: a method whose equilibrium cannot be worked
    # through gets no diagram, because it gets no Calculations section.
    diagram = force_diagram(method)
    figure_number = 0
    if diagram:
        figure_number = counter.next_figure() if counter is not None else 0
        sec.blocks.append(Figure(
            diagram, f"Forces on a slice — {label}", figure_number,
            source=f"{method} force diagram",
            width_in=FORCE_DIAGRAM_WIDTH_IN))

    preamble = _method_preamble(calc, method, figure_number)

    url = method_doc_url(method)
    # Three methods print an equation their own derivation does not publish, and
    # say so: what reaches their solution is a slice-by-slice march, and the
    # quotient under it is a balance that holds because the march closed.
    if method in WHOLE_MASS_BALANCE_METHODS:
        intro = (f"The slice-by-slice march that reaches this solution is the "
                 f"derivation published for {label}. The quotient below it is "
                 f"the horizontal force balance of the whole sliding mass at "
                 f"the converged solution, in the symbols of the XSLOPE "
                 f"documentation; the numbers in it are the converged values.")
    else:
        intro = (f"The equations below are in the symbols of the derivation "
                 f"published for {label} in the XSLOPE documentation; the "
                 f"numbers in them are the converged values.")
    # What the model does not carry is stated once, in the sentence that reduces
    # the published equations to it. Saying it here as well would say it twice,
    # the first time about equations the reader has not reached yet.
    if calc["stage"]:
        intro += (f" The governing stage is {calc['stage']}, and every number "
                  f"below is that stage's.")
    sec.blocks.append(Prose(intro, links=[(label, url)] if url else []))

    # --- the equation, or the two equations the solution is the root of ---
    #
    # Every method but one closes on a quotient, and that quotient IS the method.
    # Spencer's does not: F and θ are the pair at which two equilibrium equations
    # both vanish, which is how the derivation presents it and the only honest
    # thing to print.
    #
    # The published equation comes first and this model's below it, except for
    # the three that march: what those three print here is the balance of the
    # whole mass at the converged solution, and the equations their derivation
    # publishes are the per-slice pair the march solves, transcribed in the
    # preamble where the march is described.
    equation = []
    if calc["equation"] and method not in WHOLE_MASS_BALANCE_METHODS:
        equation.extend(_transcribed_blocks(
            calc, then=[calc["equation"]],
            always=TRANSCRIPTIONS[method].build == "balance"))
    elif calc["equation"]:
        equation.append(Math(calc["equation"]))
    else:
        equation.append(Prose(
            "Substituting Q and y_Q into the equilibrium of the whole sliding "
            "mass gives equations (27) and (28) of the derivation, two "
            "equations in the two unknowns. R_1 is the force imbalance and R_2 "
            "the moment imbalance, and the solution is the pair (F, θ) at which "
            "both are zero:"))
        equation.extend(Math(line) for line in calc["equilibrium"])
        equation.append(Prose(
            "F and θ are iterated together: each trial pair recomputes both "
            "sums, and the pair is adjusted until both vanish."))

    # --- what the solution came out at ---
    close = (_spencer_close(calc, table_number, bookmark)
             if calc["equation"] is None
             else _quotient_close(calc, table_number, bookmark, unit_labels))

    # --- what every letter means ---
    #
    # An equation printed without its nomenclature is a wall the reader either
    # already knows how to climb or does not. Only the symbols the section's own
    # equations carry are defined, and the ones that are slice-table columns are
    # marked as such, so a reader who wants a value knows to go to the table.
    #
    # Every equation the section PRINTS, which is why the whole of it is built
    # before this point: the preamble's, the equation itself, and the arithmetic
    # that closes it — Janbu's f_o line and Morgenstern-Price's two residuals are
    # down there, and drawing the nomenclature from the equation alone left their
    # letters standing in the section undefined. Nothing that is not printed goes
    # into the pool: calc["normal_force"] on Spencer's section, which prints no
    # such equation, put in a row for the distributed load under both of the
    # letters the two derivations give it.
    symbols = equation_symbols(
        " ".join(b.notation for b in preamble + equation + close
                 if b.kind == "math"),
        calc["slice_df"], unit_labels, calc.get("symbols"))
    nomenclature = []
    if symbols:
        where = f"Table {table_number}" if table_number else "the slice table"
        nomen_number = counter.next_table() if counter is not None else 0
        nomen_where, links = cite("Table", nomen_number)
        # The slice table's own cross-reference lands on the bookmark that names
        # THIS method's table, not on the number alone: a report of several
        # methods carries one slice table each.
        links = links + ([(where, f"#{bookmark}")] if table_number else [])
        defines = f"{nomen_where} defines" if nomen_where else "Defined below are"
        nomenclature.append(Prose(
            f"{defines} the symbols above, in the order they appear. Those that "
            f"are columns of {where} carry a value for every slice.",
            links=links))
        nomenclature.append(Table(
            ["Symbol", "Meaning"], [[s, m] for s, m in symbols],
            f"Nomenclature — {label}", nomen_number))

    sec.blocks.extend(preamble + equation + nomenclature + close)
    return sec


def _quotient_close(calc, table_number, bookmark, unit_labels):
    """How every method but Spencer's ends: the two sums, the division, and — for
    the methods that carry one — the correction or the second residual."""
    from .columns import BY_KEY, format_fs, format_residual, format_sum, unit_label

    # --- the sums, and where their per-slice terms are ---
    res_col = BY_KEY.get(calc["res_key"])
    drv_col = BY_KEY.get(calc["drv_key"])
    unit = unit_label(res_col, unit_labels) if res_col is not None else ""
    n_slices = len(calc["slice_df"])
    blocks = []
    if table_number and res_col is not None and drv_col is not None:
        where = f"Table {table_number}"
        in_units = f", both in {unit}" if unit else ""
        blocks.append(Prose(
            f"Each slice's contribution to the two sums is a column of {where}: "
            f"{res_col.label} is the resisting term and {drv_col.label} is the "
            f"net driving term{in_units}. Summed over the {n_slices} slices:",
            links=[(where, f"#{bookmark}")]))
    else:
        blocks.append(Prose(
            f"Summing the per-slice terms over the {n_slices} slices:"))

    blocks.append(Math(
        f"F = frac{{{format_sum(calc['resisting'])}}}"
        f"{{{format_sum(calc['driving'])}}}"))
    if calc["fo"]:
        blocks.append(Math(f"F = {format_fs(calc['quotient'])}"))
        blocks.append(Prose(
            "Janbu's correction factor f_o compensates for the neglected "
            "interslice shear. Equation (4) of the derivation reads it from the "
            "method's chart fit for this surface's depth-to-length ratio and "
            "the soil type, and equation (5) applies it to the factor of safety "
            "above:"))
        blocks.append(Math(
            f"F_corr = f_o·F = {format_sum(calc['fo'])}·"
            f"{format_fs(calc['quotient'])} = {format_fs(calc['FS'])}"))
    else:
        blocks.append(Math(f"F = {format_fs(calc['FS'])}"))

    # --- Morgenstern-Price: the moment condition the force balance is solved
    # jointly with, and what each residual came out at ---
    residuals = calc.get("residuals")
    if residuals is not None:
        blocks.append(Prose(
            "The moment of the whole sliding mass about the coordinate origin "
            "closes at the same (F, λ). At the solution the interslice force "
            "left at the far end of the march, and the moment sum, are:"))
        blocks.append(Math(f"Z_n = {format_residual(residuals[0])}"))
        blocks.append(Math(f"sum{{M_o}} = {format_residual(residuals[1])}"))
    return blocks


def _method_section(slope_data, bundle, note, method, opts, counter, figure_dir,
                    progress=None):
    """Everything the report says about ONE method, under that method's heading.

    A report that documents several methods is this section repeated, and the
    heading is the method's name — a reader who opens the document at a slice
    table must never have to page backwards to learn whose numbers are on it. The
    figures and tables inside carry the name for the same reason, and take their
    numbers from the report-wide counter, so the sequence runs unbroken through
    however many methods there are.

    **The search belongs here, not above.** A search finds the critical surface
    FOR A METHOD: run the same model under Spencer and under Bishop and the two
    searches settle on different surfaces. Documenting one search above all the
    method blocks said that every method below shared it, which is true only when
    one method was run. Each block now carries the search that produced ITS
    surface, and a method the report solved on another method's surface says so
    instead.
    """
    label = method_label(method)
    sec = Section(label)
    summary = method_summary(method)
    if summary:
        sec.blocks.append(Prose(summary))
    if bundle is None:
        # The note is a whole sentence and names the method itself; the heading
        # above it does too, so it is printed as it stands.
        sec.blocks.append(Prose(note))
        return sec

    results = bundle.get("results") or {}
    slice_df = bundle.get("slice_df")

    # --- the search that found THIS method's surface ---
    if opts["lem_search"] and bundle.get("search"):
        search = _search_section(slope_data, bundle, opts, counter, figure_dir,
                                 method, progress)
        if search is not None:
            sec.children.append(search)

    res = Section("Results")

    # Whether this method's surface was SEARCHED for or entered decides how the
    # answer may be described. "Critical" means least of a family, which is a
    # word a search earns; a factor of safety on a surface the user drew is the
    # factor of safety of that surface and of nothing else.
    searched = bool(bundle.get("search"))
    named = "critical surface" if searched else "specified surface"

    # The figure is rendered ahead of the sentence that reports the answer, so
    # that sentence can name the figure the surface is drawn on.
    figure = None
    if opts["lem_solution_figure"] and slice_df is not None:
        fpath = os.path.join(figure_dir, f"solution_{method or 'lem'}.png")

        def draw(fig):
            from .plot import plot_solution
            plot_solution(slope_data, slice_df, bundle.get("failure_surface"),
                          results, fig=fig, show_title=False,
                          style=opts.get("style"))

        if progress:
            progress(f"the {named} — {label}")
        if _render(draw, fpath, opts):
            figure = Figure(
                fpath,
                f"{named.capitalize()} and slice forces — {label}",
                # The source names what the figure IS, not what this run's
                # surface is called, so a caller looking for a method's solution
                # plot finds it whether or not that method searched.
                counter.next_figure(), source=f"{method} solution surface")
    where, links = cite("Figure", figure.number if figure is not None else 0)

    fs = _num(results.get("FS"))
    if fs is not None:
        surface = f"the {named} of {where}" if where else f"the {named}"
        text = f"{label} gives a factor of safety of {fs:.3f} on {surface}."
        if not searched:
            text += (" No search for a critical surface was performed: this "
                     "factor of safety is for the surface entered in the input "
                     "and is not a minimum over any family of surfaces.")
        res.blocks.append(Prose(
            text + (f" {note}" if note else ""),
            bold=[f"{fs:.3f}"], links=links))
    elif figure is not None:
        text = (f"{where} draws the {named} with the force on the base of each "
                f"slice.")
        if not searched:
            text += (" The surface is the one entered in the input; no search "
                     "for a critical surface was performed.")
        res.blocks.append(Prose(text + (f" {note}" if note else ""), links=links))

    warns = results.get("warnings") or []
    if warns:
        res.blocks.append(Prose(
            "The solution reported the following admissibility notes, where the "
            "computed stresses depart from what the method assumes:"))
        res.blocks.append(Bullets([str(w) for w in warns]))

    if figure is not None:
        res.blocks.append(figure)
    sec.children.append(res)

    # --- rapid drawdown ---
    if opts["lem_rapid"]:
        rapid = _rapid_section(results, counter)
        if rapid is not None:
            sec.children.append(rapid)

    # --- slice table and calculations ---
    # The calculation is worked out first: it adds the per-slice terms of the
    # factor of safety equation to the table, which is how the section can point
    # at a column instead of walking the reader through fifteen slices.
    calc, calc_note = None, ""
    if opts["lem_calculations"]:
        try:
            calc, calc_note = _calculation(slope_data, bundle, method)
        except Exception:
            import traceback
            traceback.print_exc()
    table_df = calc["slice_df"] if calc is not None else slice_df

    bookmark = f"{SLICE_TABLE_BOOKMARK}_{method}"
    table_number = 0
    if opts["lem_slice_table"] and table_df is not None:
        from .columns import slice_table, slice_totals
        headers, rows, legend = slice_table(table_df, _unit_labels(slope_data))
        totals = slice_totals(table_df)
        sub_tab = Section("Slice Table")

        # The key to the table's first column, immediately above it: the same
        # plot the results section carries, framed on the sliced mass alone and
        # with every slice labeled, so a row of the table can be found on the
        # section it describes. It is drawn — and the table's own number taken —
        # before the paragraph that introduces the two, which cites both.
        key = None
        if opts["lem_slice_key"]:
            kpath = os.path.join(figure_dir, f"slice_key_{method or 'lem'}.png")

            def draw_key(fig):
                from .plot import plot_solution
                plot_solution(slope_data, table_df, bundle.get("failure_surface"),
                              results, fig=fig, show_title=False,
                              slice_numbers=True, frame="slices",
                              style=opts.get("style"))

            if progress:
                progress(f"the slice key — {label}")
            if _render(draw_key, kpath, opts):
                # A portrait figure at text width, like every other plot in the
                # report. The key is a picture of fifteen numbered slices; it
                # needs a page no more than the critical-surface plot does, and
                # taking the landscape page the table needs cost a sheet of its
                # own for a figure that reads at a sixth of it.
                key = Figure(
                    kpath, f"Slice numbering for the table below — {label}",
                    counter.next_figure(), source=f"{method} slice key")

        table_number = counter.next_table()
        table_where, links = cite("Table", table_number)
        key_where, key_links = cite("Figure", key.number if key is not None else 0)
        links = links + key_links
        numbered = f", numbered as in {key_where}" if key_where else ""
        sub_tab.blocks.append(Prose(
            f"{table_where} holds the geometry, forces and strengths of every "
            f"slice on the {named} as solved by {label}{numbered}. "
            f"Forces are per unit thickness of section.", links=links))
        if key is not None:
            sub_tab.blocks.append(key)
        sub_tab.blocks.append(Table(
            headers, rows, f"Slice data — {label}", table_number,
            landscape=True, legend=legend, bookmark=bookmark,
            # Every column holds one short number, so the report's own policy
            # centers them all; nothing here has to say so. The totals are the
            # sums the calculations divide, at the foot of the table the terms
            # came from.
            totals=totals))
        sec.children.append(sub_tab)

    if calc is not None:
        # The force diagram is a figure like any other as far as a caller
        # counting them is concerned, even though it is drawn by hand and shipped
        # rather than rendered here.
        if progress and force_diagram(method):
            progress(f"the force diagram — {label}")
        sec.children.append(_calculations_section(
            calc, slope_data, table_number, _unit_labels(slope_data), bookmark,
            counter))
    elif calc_note:
        # A method whose equilibrium cannot be worked through says what is true
        # of it instead, at the foot of its own block — where the working would
        # have stood. A factor of safety with nothing after it reads as an
        # omission, and a passive model's is not one.
        (sec.children[-1] if sec.children else sec).blocks.append(
            Prose(calc_note))
    return sec


def _lem_section(slope_data, solutions, opts, counter, figure_dir, progress=None):
    bundle = select_bundle(solutions, opts.get("method"))
    if bundle is None:
        return None
    methods = featured_methods(solutions, opts)
    slice_df = bundle.get("slice_df")

    sec = Section("Limit Equilibrium Analysis")
    sec.blocks.append(Prose(
        "The stability of the section was evaluated by the method of slices. The "
        "factor of safety is the factor by which the available shear strength "
        "along the failure surface must be divided to bring the sliding mass into "
        "limiting equilibrium."))

    # --- engine inputs ---
    #
    # A circle a search started from and a circle that IS the analysis are the
    # same row of the input sheet and two different statements about the answer,
    # so the row says which it was.
    any_search = any((select_bundle(solutions, m) or {}).get("search")
                     for m in methods)
    items = [(("Method" if len(methods) == 1 else "Methods") + " reported in detail",
              _join([method_label(m) for m in methods]))]
    if slice_df is not None:
        items.append(("Slices", str(len(slice_df))))
    circular = bool(slope_data.get("circular"))
    items.append(("Surface family", "circular" if circular else "non-circular"))
    if circular and slope_data.get("circles"):
        c = slope_data["circles"][0]
        items.append(("Starting circle" if any_search else "Specified circle",
                      f"center ({_fmt(c.get('Xo'), '{:g}')}, "
                      f"{_fmt(c.get('Yo'), '{:g}')}), R = {_fmt(c.get('R'), '{:g}')}"))
    elif slope_data.get("non_circ"):
        items.append(("Non-circular surface",
                      f"{len(slope_data['non_circ'])} defining points"))
    k = _num(slope_data.get("k_seismic"))
    items.append(("Surface", "located by search" if any_search
                  else "specified in the input; no search was performed"))
    items.append(("Seismic coefficient", f"{k:g}" if k else "none"))
    tc = _num(slope_data.get("tcrack_depth"))
    if tc:
        items.append(("Tension crack depth", f"{tc:g}"))
        items.append(("Tension crack water", f"{_num(slope_data.get('tcrack_water')) or 0:g}"))
    md = _num(slope_data.get("max_depth"))
    if md:
        items.append(("Maximum surface depth (elevation)", f"{md:g}"))
    sub_inputs = Section("Analysis Inputs")
    sub_inputs.blocks.append(KeyValues(items))
    sec.children.append(sub_inputs)

    # --- the properties and the water this engine reads ---
    #
    # A material's shear strength and its pore pressure are limit equilibrium
    # inputs, so they are stated where the limit equilibrium analysis is
    # documented rather than as a general description of the section: the
    # seepage analysis reads neither, and gives its conductivities in its own
    # section.
    feats = water_features(slope_data)
    if opts["lem_materials"]:
        sub = Section("Materials")
        table = _materials_table(slope_data, counter)
        if table is not None:
            where, links = cite("Table", table.number)
            sub.blocks.append(Prose(
                f"Every material the section geometry references is given in "
                f"{where}, with the strength option it is analyzed under, the "
                f"properties that option uses, and how its pore pressure is "
                f"taken.", links=links))
            sub.blocks.append(table)
        else:
            sub.blocks.append(Prose("The model defines no materials."))
        if not feats["any"]:
            sub.blocks.append(Prose(
                "The model defines no groundwater and no external water; the "
                "section is analyzed dry, with zero pore pressure throughout."))
        else:
            if feats["pore"]:
                # In words, not in the sheet's own codes: the table prints the
                # code against each material, and the sentence says what it is.
                said = _join([PORE_SOURCES.get(p, p) for p in feats["pore"]])
                sub.blocks.append(Prose(
                    f"Pore pressure on the base of a slice is taken from "
                    f"{said}, as the table states for each material."))
            rows = _water_items(slope_data, feats)
            if rows:
                sub.blocks.append(KeyValues(rows))
        sec.children.append(sub)

    if opts["lem_loads"]:
        loads = _loads_section(slope_data, feats, counter)
        sec.children.append(loads)
        printed = [b.number for b in loads.blocks if b.kind == "table"]
        if printed:
            # The finite element section applies the same blocks; it cites this
            # table rather than setting the same numbers again.
            opts["_loads_table_number"] = printed[0]

    # --- every method's answer, once, ahead of the detail ---
    #
    # The search does NOT belong here. It finds the critical surface for one
    # METHOD, and two methods searched separately settle on different surfaces;
    # each method's block carries its own search.
    table = _fs_table(slope_data, solutions, opts, counter)
    if table is not None:
        sub_fs = Section("Factors of Safety")
        searched = [m for m in methods
                    if (select_bundle(solutions, m) or {}).get("search")]
        where, links = cite("Table", table.number)
        text = (f"{where} gives the factor of safety each method reported. Every "
                f"method finds its own surface, so each row is that method's "
                f"answer on the surface stated beside it and no row is another "
                f"method's answer on the same surface.")
        if searched and len(searched) < len(methods):
            text += (f" {_join([method_label(m) for m in searched])} searched for "
                     f"{'its' if len(searched) == 1 else 'their'} surface; the "
                     f"rest were solved on the surface specified in the input, "
                     f"which is not a minimum over any family of surfaces.")
        elif not searched:
            text += (" No search was performed: every factor of safety here is "
                     "for the surface specified in the input, and is not a "
                     "minimum over any family of surfaces.")
        text += " Each method is then reported in full below."
        sub_fs.blocks.append(Prose(text, links=links))
        sub_fs.blocks.append(table)
        sec.children.append(sub_fs)

    # --- one full detail block per featured method ---
    for name in methods:
        detail, note = detail_bundle(slope_data, solutions, name)
        sec.children.append(_method_section(
            slope_data, detail, note, name, opts, counter, figure_dir, progress))

    return sec


# ---------------------------------------------------------------------------
# Seepage
# ---------------------------------------------------------------------------

def _bc_counts(seep_data):
    """``(specified-head nodes, exit-face nodes)`` — the boundary of the flow
    problem as the solver reads it, counted off the same array the solve used."""
    bc_type = (seep_data or {}).get("bc_type")
    if bc_type is None:
        return 0, 0
    head = exit_face = 0
    for t in bc_type:
        if int(t) == 1:
            head += 1
        elif int(t) == 2:
            exit_face += 1
    return head, exit_face


def _seep_results_section(slope_data, bundle, title, tag, named, opts, counter,
                          figure_dir, mesh_numbers, progress=None):
    """One solved boundary condition set: what it was, its flow net, its flow.

    ``named`` is what the figure caption calls this set, and is empty for a model
    solved for one: a caption reading "Flow net" on the only flow net there is
    needs no qualifier. ``mesh_numbers`` maps each set's tag to the figure number
    its mesh and boundary conditions were drawn under, among the inputs.
    """
    seep_data = bundle.get("seep_data") or {}
    solution = bundle.get("solution") or {}
    sub = Section(title)

    # The mesh and its boundary conditions are an INPUT to this solve and were
    # drawn with the rest of the inputs; the number is carried here so the
    # paragraph that reports the solve can point back at the boundaries it names.
    mesh_number = mesh_numbers.get(tag, 0)

    figure = None
    if opts["seep_flownet"]:
        path = os.path.join(figure_dir, f"seep_{tag}.png")

        def draw(fig):
            from .plot_seep import plot_seep_solution, flownet_base_material
            # mesh=False: this is a flow net. Element edges chop the head
            # contours and the flow lines into a dashed look and hide the field
            # under a grid — the same call the shipped seepage figures make. The
            # base material is the zone the net is scaled to, chosen by the same
            # rule those figures were built on.
            plot_seep_solution(seep_data, solution, fig=fig, show_title=False,
                               mesh=False,
                               base_mat=flownet_base_material(seep_data, solution),
                               style=opts.get("style"))

        if progress:
            progress("the flow net" + (f" — {named}" if named else ""))
        if _render(draw, path, opts):
            figure = Figure(path, "Flow net" + (f" — {named}" if named else ""),
                            counter.next_figure(), source=f"seepage {tag}")
    where, links = cite("Figure", figure.number if figure is not None else 0)
    mesh_where, mesh_links = cite("Figure", mesh_number)
    links = list(links) + mesh_links

    n_head, n_exit = _bc_counts(seep_data)
    if n_exit:
        text = (f"Flow was solved as an unconfined problem: the phreatic surface "
                f"is located as part of the solution rather than prescribed. "
                f"{n_head:,} nodes carry a specified head and {n_exit:,} lie on "
                f"an exit face, where water leaves the section at atmospheric "
                f"pressure.")
        drawn = "the head contours, the phreatic surface and the flowlines"
    else:
        text = (f"Flow was solved as a confined problem: every node of the mesh "
                f"flows saturated, and {n_head:,} of them carry a specified "
                f"head.")
        drawn = "the head contours and the flowlines"
    if mesh_where:
        text += f" {mesh_where} shows where each of those boundaries falls."
    if where:
        text += f" {where} draws {drawn}."
    sub.blocks.append(Prose(text, links=links))

    q = _num(solution.get("flowrate"))
    if q is not None:
        lbl = _unit_labels(slope_data) or {}
        unit = lbl.get("flowrate") or ""
        amount = f"{q:.4g} {unit}".strip()
        tail = "" if unit else " per unit thickness of section"
        sub.blocks.append(Prose(
            f"The flow through the section is {amount}{tail}.", bold=[amount]))

    if figure is not None:
        sub.blocks.append(figure)
    return sub


def _seep_section(slope_data, solutions, opts, counter, figure_dir, progress=None):
    """What the seepage analysis solved, on what mesh, and what came out of it.

    The section stands ahead of the stability analyses because that is the order
    the numbers flow in: a material whose pore pressure is taken from seepage is
    analyzed on the field computed here.
    """
    bundles = seep_bundles(solutions)
    if not bundles:
        return None

    sec = Section("Seepage Analysis")
    text = ("Flow through the section was solved by the finite element method. "
            "The unknown at every node is total head, and the pore pressure at a "
            "node is the height of head above it times the unit weight of water.")
    if "seep" in water_features(slope_data)["pore"]:
        text += (" Every material whose pore pressure is taken from seepage is "
                 "analyzed on this field.")
    sec.blocks.append(Prose(text))

    # --- engine inputs ---
    sub_inputs = Section("Analysis Inputs")
    seep_data = bundles[0].get("seep_data") or {}

    # The model as the flow solver reads it: the same section the stability
    # analysis runs on, but carrying the water surfaces the head boundaries
    # state instead of the trial surfaces and the loads.
    model = None
    if opts["seep_inputs_figure"]:
        mpath = os.path.join(figure_dir, "seep_inputs.png")

        def draw_model(fig):
            from .plot import plot_inputs
            plot_inputs(slope_data, fig=fig, mode="seep", show_title=False,
                        frame="content", style=opts.get("style"),
                        show_mesh=False)

        if progress:
            progress("the seepage model")
        if _render(draw_model, mpath, opts):
            model = Figure(mpath, "Seepage model", counter.next_figure(),
                           source="seep model")
    if model is not None:
        where, links = cite("Figure", model.number)
        sub_inputs.blocks.append(Prose(
            f"{where} shows the flow domain: its material zones and the water "
            f"surface each specified-head boundary states.", links=links))
        sub_inputs.blocks.append(model)

    items = []
    summary = mesh_summary(seep_data)
    if summary:
        items.append(("Mesh", summary))
    gamma_w = _num(seep_data.get("unit_weight"))
    if gamma_w is not None:
        lbl = _unit_labels(slope_data) or {}
        items.append(("Unit weight of water",
                      f"{gamma_w:g} {lbl.get('unit_weight', '')}".strip()))
    if items:
        sub_inputs.blocks.append(KeyValues(items))
    if opts["seep_materials"]:
        table = _seep_materials_table(slope_data, counter)
        if table is not None:
            where, links = cite("Table", table.number)
            # What the table actually carries, column for column. Naming only
            # the conductivities left the unsaturated columns — the ones that
            # decide where the phreatic surface settles — printed and
            # unaccounted for.
            heads = " ".join(table.headers)
            text = (f"{where} gives the properties of every material the flow "
                    f"domain carries: the major and minor saturated "
                    f"conductivities, and the angle the major axis makes with "
                    f"the horizontal.")
            unsat = []
            if "Unsaturated" in heads:
                unsat.append("the unsaturated model it is assigned")
            for key, name in (("k_r0", "the relative conductivity it falls to "
                               "when dry"),
                              ("h₀", "the pressure head it falls off over"),
                              ("a", "the curve parameters that shape the "
                               "fall-off")):
                if any(h == key or h.startswith(key + " ") for h in table.headers):
                    unsat.append(name)
            if unsat:
                text += (f" Above the phreatic surface the conductivity is "
                         f"reduced, and the table also gives {_join(unsat)}.")
            sub_inputs.blocks.append(Prose(text, links=links))
            sub_inputs.blocks.append(table)

    # The unsaturated models, drawn: the parameters in the table are three
    # different functions, and what they mean is the shape of the curve. All the
    # materials go on one axes, so they are read against each other. A model
    # whose materials are all saturated has no curve, and no figure.
    kr_materials = _kr_materials(slope_data) if opts["seep_kr_figure"] else []
    if kr_materials:
        kpath = os.path.join(figure_dir, "seep_kr.png")

        def draw_kr(fig):
            from .plot import plot_material_kr_set
            plot_material_kr_set(kr_materials, fig=fig, show_title=False,
                                 style=opts.get("style"),
                                 unit_labels=_unit_labels(slope_data))

        if progress:
            progress("the unsaturated conductivity curves")
        if _render(draw_kr, kpath, opts):
            figure = Figure(kpath, "Unsaturated relative conductivity",
                            counter.next_figure(), source="seep kr")
            where, links = cite("Figure", figure.number)
            sub_inputs.blocks.append(Prose(
                f"{where} is the reduction each material's unsaturated model "
                f"applies: the factor its saturated conductivity is multiplied "
                f"by at a given matric suction, evaluated by the same functions "
                f"the flow solver evaluates.", links=links))
            sub_inputs.blocks.append(figure)

    # The mesh and the boundary conditions on it: an input to the flow problem,
    # not an outcome of it, so it stands with the inputs. One per solved set — a
    # rapid drawdown model is two different boundary problems on one mesh, and a
    # figure of the mesh alone would show neither.
    tags = []
    for i, bundle in enumerate(bundles):
        bc = (bundle.get("options") or {}).get("bc")
        number = bc if bc is not None else i + 1
        tags.append((bundle, f"bc{number}",
                     "" if len(bundles) == 1
                     else f"boundary condition set {number}", number))

    mesh_numbers = {}
    if opts["seep_mesh_figure"]:
        for bundle, tag, named, _number in tags:
            data = bundle.get("seep_data") or {}
            mpath = os.path.join(figure_dir, f"seep_mesh_{tag}.png")

            def draw_mesh(fig, data=data):
                from .plot_seep import plot_seep_data
                plot_seep_data(data, fig=fig, show_title=False, show_bc=True,
                               style=opts.get("style"))

            if progress:
                progress("the seepage mesh" + (f" — {named}" if named else ""))
            if _render(draw_mesh, mpath, opts):
                figure = Figure(
                    mpath,
                    "Seepage mesh and boundary conditions"
                    + (f" — {named}" if named else ""),
                    counter.next_figure(), source=f"seepage {tag} mesh")
                mesh_numbers[tag] = figure.number
                where, links = cite("Figure", figure.number)
                on = f" — {summary} —" if summary else ""
                sub_inputs.blocks.append(Prose(
                    f"{where} is the mesh the flow was solved on{on} colored by "
                    f"material, with every specified-head and exit-face node "
                    f"marked"
                    + (f" for {named}." if named else "."), links=links))
                sub_inputs.blocks.append(figure)
    sec.children.append(sub_inputs)

    # --- one block per solved boundary condition set ---
    #
    # The runner emits one bundle per set, and a rapid drawdown model is run for
    # two: the full pool and the drawn-down pool are different flow problems on
    # the same mesh, so each gets its own block rather than a shared one that
    # would have to describe both.
    for bundle, tag, named, number in tags:
        title = ("Results" if len(bundles) == 1
                 else f"Boundary Condition Set {number}")
        sec.children.append(_seep_results_section(
            slope_data, bundle, title, tag, named, opts, counter,
            figure_dir, mesh_numbers, progress))
    return sec


# ---------------------------------------------------------------------------
# Finite element deformation and strength reduction
# ---------------------------------------------------------------------------

#: What each finite element results panel is called and what it draws, in the
#: order the report prints them. The panel names are :func:`plot_fem_results`'s
#: own ``plot_type`` values.
FEM_PANELS = (
    ("deformation", "Deformed mesh",
     "the deformed mesh over the original section"),
    ("shear_strain", "Maximum shear strain",
     "the viscoplastic shear strain, which is where the section is shearing"),
    ("displace_vector", "Displacement vectors",
     "the displacement of every node as an arrow, which is how the section is "
     "moving"),
)

#: The published page each kind of one-dimensional member is modelled after.
#: Both formulations differ from the limit equilibrium treatment of the same
#: member, so the paragraph that describes one links the page that derives it.
FEM_DETAIL_DOC_PAGES = {
    "reinforcement": "fem/reinforcement.md",
    "pile": "fem/piles.md",
}

#: The field the member profiles are read at. A strength reduction run is asked
#: about the mechanism it developed, and where no at-failure snapshot was
#: captured :func:`xslope.fem_details.effective_field_state` falls back to the
#: converged field — the selection, and the fallback, the results view's own
#: Field state control opens on.
DETAIL_FIELD_STATE = "failure"

#: The two kinds of one-dimensional member a finite element run can carry: the
#: option that prints the member's subsection, the option that prints its detail
#: figures, the heading it takes, the stem its figures are written under, what
#: one member is called and what several are.
DETAIL_KINDS = {
    "reinforcement": {
        "option": "fem_reinforcement",
        "figure_option": "fem_reinforcement_figure",
        # Headed apart from the Project Definition section of the same name: one
        # is the reinforcement as entered, this one is what the analysis put in
        # it, and a contents page with two "Reinforcement" entries says neither.
        "title": "Reinforcement Forces",
        "tag": "reinf",
        "one": "line",
        "many": "lines",
    },
    "pile": {
        "option": "fem_piles",
        "figure_option": "fem_piles_figure",
        "title": "Pile Forces",
        "tag": "pile",
        "one": "pile",
        "many": "piles",
    },
}

#: How each kind of member is modelled, in the terms its documentation page
#: uses. The linked phrase is the element formulation itself, which is what
#: separates this treatment from the limit equilibrium one.
DETAIL_MODELLING = {
    "reinforcement": (
        "two-node truss elements",
        "Each reinforcement line is discretized into two-node truss elements on "
        "the mesh's own nodes, carrying axial tension only. What one can hold at "
        "a point along the line is the pullout resistance developed from the "
        "nearer free end over its development length, or the tensile capacity "
        "T_max where enough length has developed; the force the ground hands the "
        "bar per unit of its length is the gradient of the axial force along "
        "it."),
    "pile": (
        "Euler-Bernoulli beam elements",
        "Each pile is discretized into Euler-Bernoulli beam elements on the "
        "mesh's own nodes, with a rotational degree of freedom at each. Its "
        "resistance to the moving ground follows from its bending stiffness "
        "rather than from a force applied to it, and is limited by the shear and "
        "moment capacities the model declares."),
}

#: What one detail figure draws, per kind, for the sentence that cites it.
DETAIL_FIGURE_SHOWS = {
    "reinforcement": ("the mobilized axial force over the declared capacity "
                      "envelope, with the bond transfer rate beneath it"),
    "pile": ("the lateral displacement, the shear, the bending moment and the "
             "mobilized soil reaction against depth"),
}

#: What one detail figure's caption calls it, before the member's own name.
DETAIL_FIGURE_CAPTIONS = {
    "reinforcement": "Axial force and bond transfer along",
    "pile": "Displacement, shear, moment and soil reaction along",
}


def _percent(value):
    """``0.62`` -> ``"62%"``, and ``""`` for an unmeasurable utilization."""
    n = _num(value)
    return "" if n is None else f"{n:.0%}"


def _series(profile, key):
    """One of a profile's along-the-member series, as a plain list. The series
    are numpy arrays, which have no truth value, so an empty one is asked for by
    its length and never by ``or []``."""
    values = (profile or {}).get(key)
    return [] if values is None or len(values) == 0 else list(values)


def _detail_members(slope_data, bundle, kind):
    """Every reinforcement line, or every pile, one finite element run carries.

    Enumerated through :func:`xslope.fem_details.list_lines`, which is what the
    results view's details panel lists, at the same field state — so a member the
    report describes is one the solved model actually owns elements for, and not
    a row of an input sheet the mesh never reached.
    """
    try:
        from .fem_details import list_lines
        found = list_lines(bundle.get("fem_data") or {},
                           bundle.get("solution") or {}, slope_data,
                           field_state=DETAIL_FIELD_STATE,
                           failure_solution=bundle.get("failure_solution"))
    except Exception:
        import traceback
        traceback.print_exc()
        return []
    return [m for m in found if m.get("kind") == kind]


def _detail_profiles(slope_data, bundle, kind):
    """The along-the-member profile of every member of one kind, in list order.

    One profile is everything both the table and the detail figure read, so the
    number in a row and the curve in the figure beside it are the same series.
    """
    from . import fem_details
    out = []
    for member in _detail_members(slope_data, bundle, kind):
        read = (fem_details.pile_profile if kind == "pile"
                else fem_details.reinforcement_profile)
        try:
            out.append(read(bundle.get("fem_data") or {},
                            bundle.get("solution") or {}, member["index"],
                            slope_data, field_state=DETAIL_FIELD_STATE,
                            failure_solution=bundle.get("failure_solution")))
        except Exception:
            import traceback
            traceback.print_exc()
    return out


def _figured_members(profiles):
    """Which members of one kind get a detail figure: every one of them.

    A member the analysis solved is a member the report draws. The detail figure
    is what says where along a bar the force peaks and whether the bond carried
    it there, and that question is asked of each bar separately — a table row
    gives the peak, not the profile that produced it.
    """
    return list(profiles)


def _detail_units(profiles):
    """``(force, length, moment)`` unit suffixes for the member tables, from the
    model's own declared system."""
    u = (profiles[0].get("units") if profiles else None) or {}

    def suffix(key):
        return f" ({u[key]})" if u.get(key) else ""

    return suffix("force"), suffix("length"), suffix("moment")


def _reinforcement_forces_table(slope_data, profiles, counter):
    """What the solution put in every reinforcement line: the capacity the model
    declares, the force at the point of greatest utilization, and where that
    point is."""
    lines = slope_data.get("reinforcement_lines") or []
    fu, lu, _mu = _detail_units(profiles)
    own = [lines[p["index"] - 1] if 0 <= p["index"] - 1 < len(lines) else {}
           for p in profiles]
    # A residual capacity column only where some line declares one: a line that
    # never softens has no residual, and a column of blanks states nothing.
    softens = _populated(own, "t_res")
    # The force and the position are those of the point of GREATEST
    # UTILIZATION, which on a line whose capacity ramps down towards a free end
    # is not the point of greatest force. Headed "Force" and "Position" rather
    # than "Peak force", which would read as the largest force in the bar and
    # is a different number; the sentence that cites the table says which point
    # they belong to, and the detail figure annotates that same point.
    headers = (["Line", f"T_max{fu}"] + ([f"T_res{fu}"] if softens else [])
               + [f"Force{fu}", f"Position{lu}", "Utilization", "State"])
    rows = []
    for profile, line in zip(profiles, own):
        row = [profile["label"], _fmt(line.get("t_max"), "{:,.1f}")]
        if softens:
            row.append(_fmt(line.get("t_res"), "{:,.1f}"))
        row += [_fmt(profile.get("peak_T"), "{:,.1f}"),
                _fmt(profile.get("peak_s"), "{:.2f}"),
                _percent(profile.get("peak_utilization")),
                str(profile.get("status") or "")]
        rows.append(row)
    return Table(headers, rows, "Reinforcement forces", counter.next_table())


def _pile_forces_table(profiles, counter):
    """What the solution put in every pile: the largest shear and moment along
    it, the depth of that moment, and how far its head moved."""
    fu, lu, mu = _detail_units(profiles)
    headers = ["Pile", f"Length{lu}", f"Peak shear{fu}", f"Peak moment{mu}",
               f"At depth{lu}", f"Head movement{lu}", "Utilization", "State"]
    rows = []
    for profile in profiles:
        shear = _series(profile, "shear")
        peak_v = max((abs(_num(v) or 0.0) for v in shear), default=None)
        lateral = _series(profile, "u_lateral")
        head = _num(lateral[0]) if lateral else None
        rows.append([
            profile["label"], _fmt(profile.get("length"), "{:.2f}"),
            _fmt(peak_v, "{:,.1f}"), _fmt(profile.get("max_moment"), "{:,.1f}"),
            _fmt(profile.get("max_moment_depth"), "{:.2f}"),
            _fmt(head, "{:.4g}"),
            _percent(profile.get("peak_utilization")),
            str(profile.get("status") or ""),
        ])
    return Table(headers, rows, "Pile forces", counter.next_table())


def _detail_section(slope_data, bundle, kind, tag, opts, counter, figure_dir,
                    progress=None):
    """One kind of one-dimensional member in one finite element run.

    Returns None where the run carries no member of this kind, or where the
    options switched the subsection off: a model with no reinforcement gets no
    reinforcement heading, on the same rule the Project Definition sections
    follow.
    """
    spec = DETAIL_KINDS[kind]
    if not opts[spec["option"]]:
        return None
    profiles = _detail_profiles(slope_data, bundle, kind)
    if not profiles:
        return None

    sec = Section(spec["title"])
    phrase, modelling = DETAIL_MODELLING[kind]
    url = docs_url(FEM_DETAIL_DOC_PAGES[kind])
    sec.blocks.append(Prose(modelling, links=[(phrase, url)] if url else []))

    # Which field the forces were read from, named the way the results view
    # names it: the developed mechanism where the run captured one, and the
    # converged field where it did not.
    from .fem_details import field_state_label
    state = field_state_label(profiles[0].get("field_state", "converged"))
    read_at = ("The forces are read from the developed mechanism at failure."
               if state == "at failure" else
               "The forces are read from the last converged field.")

    if kind == "pile":
        table = _pile_forces_table(profiles, counter)
        gives = ("its length, the largest shear and bending moment along it, "
                 "the depth of that moment, the lateral displacement of its "
                 "head, and the utilization those reach against the capacity "
                 "each pile is measured by")
    else:
        table = _reinforcement_forces_table(slope_data, profiles, counter)
        gives = ("the capacity the model declares, the axial force at the point "
                 "of greatest utilization, the position of that point measured "
                 "from the first end of the line, and the utilization there")
    where, table_links = cite("Table", table.number)
    sec.blocks.append(Prose(
        f"{where} gives every {spec['one']} the analysis solved: {gives}. "
        f"{read_at}", links=table_links))
    sec.blocks.append(table)

    figures = []
    if opts[spec["figure_option"]]:
        chosen = _figured_members(profiles)
        for profile in chosen:
            path = os.path.join(
                figure_dir, f"fem_{tag}_{spec['tag']}{profile['index']}.png")

            def draw(fig, profile=profile):
                from .plot_fem_details import plot_detail
                plot_detail(profile, fig=fig)

            if progress:
                progress(f"the {spec['one']} detail — {profile['label']}")
            if _render(draw, path, opts):
                figures.append(Figure(
                    path,
                    f"{DETAIL_FIGURE_CAPTIONS[kind]} {profile['label']}",
                    counter.next_figure(),
                    source=f"fem {tag} {kind} {profile['index']}"))

    if figures:
        cites, links = [], list(table_links)
        for figure in figures:
            named, link = cite("Figure", figure.number)
            cites.append(named)
            links += link
        verb = "draws" if len(cites) == 1 else "draw"
        text = (f"Along each {spec['one']}, {_join(cites)} {verb} "
                f"{DETAIL_FIGURE_SHOWS[kind]}.")
        if len(figures) < len(profiles):
            # A profile whose plot could not be produced still has a row: the
            # table is the record, and the sentence says how many are only there.
            rest = len(profiles) - len(figures)
            text += (f" The remaining {rest} of the {len(profiles)} "
                     f"{spec['many']} are given in {where} alone.")
        sec.blocks.append(Prose(text, links=links))
        for figure in figures:
            sec.blocks.append(figure)
    return sec


def _fem_results_section(slope_data, bundle, title, tag, opts, counter,
                         figure_dir, progress=None):
    """One finite element run: its figures and the answer it reached."""
    fem_data = bundle.get("fem_data") or {}
    solution = bundle.get("solution") or {}
    failure = bundle.get("failure_solution")
    ssrm = str(bundle.get("analysis")) == "ssrm"
    label = "strength reduction" if ssrm else "finite element"
    sub = Section(title)

    # The panels are rendered one to a figure rather than stacked into one: each
    # is then the size every other figure in the report is, instead of a third of
    # it. They are drawn ahead of the paragraph that reports the answer, which
    # names them.
    figures = []
    if opts["fem_figure"]:
        for panel, caption, _shows in FEM_PANELS:
            path = os.path.join(figure_dir, f"fem_{tag}_{panel}.png")

            def draw(fig, panel=panel):
                from .plot_fem import plot_fem_results
                plot_fem_results(fem_data, solution, plot_type=[panel],
                                 fig=fig, fs=_num(bundle.get("FS")),
                                 failure_solution=failure)

            if progress:
                progress(f"the {caption.lower()} — {label}")
            if _render(draw, path, opts):
                figures.append(Figure(path, f"{caption} — {label} analysis",
                                      counter.next_figure(),
                                      source=f"fem {tag} {panel}"))

    links = []
    named = []
    for figure, (_panel, _caption, shows) in zip(figures, FEM_PANELS):
        where, link = cite("Figure", figure.number)
        links += link
        named.append(f"{where} draws {shows}")
    drawn = f" {_join(named)}." if named else ""

    fs = _num(bundle.get("FS"))
    if ssrm and fs is not None:
        state = ("The field drawn is the mechanism at failure — the trial the "
                 "section could not reach equilibrium under."
                 if failure is not None else
                 "The field drawn is the last trial that reached equilibrium.")
        sub.blocks.append(Prose(
            f"The strength reduction method gives a factor of safety of "
            f"{fs:.3f}. {state}{drawn}", bold=[f"{fs:.3f}"], links=links))
    else:
        F = _num(solution.get("F"))
        at = f" at a strength reduction factor of {F:.3f}" if F is not None else ""
        u_max = _num(solution.get("max_displacement"))
        lbl = _unit_labels(slope_data) or {}
        length = lbl.get("length") or ""
        moved = (f" The largest computed displacement is {u_max:.4g} "
                 f"{length}.".replace(" .", ".") if u_max is not None else "")
        sub.blocks.append(Prose(
            f"The section was solved for its displacements under gravity{at}. "
            f"No strength reduction was run, so this analysis reports no factor "
            f"of safety.{moved}{drawn}", links=links))

    for figure in figures:
        sub.blocks.append(figure)

    # The members the run carries, after the fields it solved: what a bar or a
    # pile ended up holding is read off the same solution the fields are.
    for kind in DETAIL_KINDS:
        child = _detail_section(slope_data, bundle, kind, tag, opts, counter,
                                figure_dir, progress)
        if child is not None:
            sub.children.append(child)
    return sub


def _fem_section(slope_data, solutions, opts, counter, figure_dir, progress=None):
    """The finite element analysis: how it models the section, and what it found."""
    bundles = fem_bundles(solutions)
    if not bundles:
        return None
    ssrm = any(str(b.get("analysis")) == "ssrm" for b in bundles)

    sec = Section("Deformation and Strength Reduction" if ssrm
                  else "Deformation Analysis")
    text = ("The section was analyzed by the finite element method. Each material "
            "is linearly elastic below its Mohr-Coulomb yield surface and "
            "perfectly plastic on it, and the stresses that satisfy equilibrium "
            "under gravity without violating that surface are found by the "
            "viscoplastic algorithm of Griffiths and Lane. No failure surface is "
            "assumed: where the section shears is an outcome of the solution.")
    if ssrm:
        text += (" In the strength reduction method the cohesion and the tangent "
                 "of the friction angle of every material are divided by a trial "
                 "factor, and the solution is repeated at increasing factors "
                 "until equilibrium can no longer be reached. The factor of "
                 "safety is the factor at that margin.")
    sec.blocks.append(Prose(text))

    # --- engine inputs ---
    sub_inputs = Section("Analysis Inputs")
    fem_data = bundles[0].get("fem_data") or {}

    # The model as the finite element solver reads it: the section with the
    # strength reduction zones and the members it carries, ahead of the mesh it
    # was discretized onto.
    model = None
    if opts["fem_inputs_figure"]:
        mpath = os.path.join(figure_dir, "fem_inputs.png")

        def draw_model(fig):
            from .plot import plot_inputs
            # No mesh underlay: the mesh gets a figure of its own directly
            # below, and drawn twice it is a grid over the zones this figure is
            # there to show.
            plot_inputs(slope_data, fig=fig, mode="fem", show_title=False,
                        frame="content", style=opts.get("style"),
                        show_mesh=False)

        if progress:
            progress("the finite element model")
        if _render(draw_model, mpath, opts):
            model = Figure(mpath, "Finite element model", counter.next_figure(),
                           source="fem model")
    if model is not None:
        where, links = cite("Figure", model.number)
        sub_inputs.blocks.append(Prose(
            f"{where} is the section the analysis was run on: the material "
            f"zones the properties below belong to, and the members the "
            f"solution carries.", links=links))
        sub_inputs.blocks.append(model)

    mesh_figure = None
    if opts["fem_mesh_figure"]:
        gpath = os.path.join(figure_dir, "fem_mesh.png")

        def draw_grid(fig):
            # The mesh as the analysis was set up on it — the same presentation
            # Studio's FEM data view draws, from the same function, so a reader
            # comparing the report to the screen is comparing one figure to
            # itself. The fixities are half of what the figure is for: a mesh
            # without them does not say what was held.
            from .plot_fem import plot_fem_data
            plot_fem_data(fem_data, fig=fig, show_title=False, show_bc=True,
                          show_nodes=False, style=opts.get("style"))

        if progress:
            progress("the finite element mesh")
        if _render(draw_grid, gpath, opts):
            mesh_figure = Figure(gpath,
                                 "Finite element mesh and boundary conditions",
                                 counter.next_figure(), source="fem mesh")

    items = []
    summary = mesh_summary(fem_data)
    if summary:
        items.append(("Mesh", summary))
    k0 = _num(fem_data.get("k0"))
    if k0 is not None:
        items.append(("Initial stress state", f"K₀ = {k0:g}"))
    if items:
        sub_inputs.blocks.append(KeyValues(items))
    if mesh_figure is not None:
        where, links = cite("Figure", mesh_figure.number)
        on = f" — {summary}" if summary else ""
        sub_inputs.blocks.append(Prose(
            f"{where} is the mesh the section was discretized onto{on}, "
            f"colored by the material each element carries, with the fixities "
            f"the solution was found under marked on the nodes that carry "
            f"them.", links=links))
        sub_inputs.blocks.append(mesh_figure)
    if opts["fem_materials"]:
        table = _fem_materials_table(slope_data, counter)
        if table is not None:
            where, links = cite("Table", table.number)
            sub_inputs.blocks.append(Prose(
                f"{where} gives the properties every element is solved with: the "
                f"unit weight and Mohr-Coulomb strength that set when it yields, "
                f"and the Young's modulus and Poisson's ratio that set how it "
                f"deforms before it does.", links=links))
            sub_inputs.blocks.append(table)
    sec.children.append(sub_inputs)

    # The loads this engine applies. A limit equilibrium report has already
    # presented the same blocks — same points, same pressures — so where both
    # engines are documented the finite element section says which loads it
    # carries and leaves the table where it stands rather than printing it twice.
    if opts["fem_loads"]:
        sec.children.append(_loads_section(
            slope_data, water_features(slope_data), counter, seismic=False,
            already=opts.get("_loads_table_number") or 0))

    for i, bundle in enumerate(bundles):
        title = "Results" if len(bundles) == 1 else f"Run {i + 1}"
        sec.children.append(_fem_results_section(
            slope_data, bundle, title, f"run{i + 1}", opts, counter, figure_dir,
            progress))
    return sec


def _model_checks_section(slope_data, solutions, opts, counter):
    """The preflight findings that were live when the analysis ran."""
    report = opts.get("preflight")
    if report is None:
        try:
            from .preflight import preflight
            bundle = select_bundle(solutions, opts.get("method"))
            selection = {"method": bundle_method(bundle) if bundle else None,
                         "search": bool((bundle or {}).get("search"))}
            report = preflight(slope_data, "lem", selection=selection)
        except Exception:
            import traceback
            traceback.print_exc()
            return None

    sec = Section("Model Checks")
    findings = relevant_findings(getattr(report, "findings", []) or [],
                                 report_analyses(solutions, opts))
    if not findings:
        sec.blocks.append(Prose(
            "xslope checks a model against what the selected analysis needs "
            "before it runs. The checks raised no findings for the analyses in "
            "this report."))
        return sec

    order = {"error": 0, "warning": 1, "info": 2}
    findings.sort(key=lambda f: order.get(f.severity, 3))
    words = {"error": "Error", "warning": "Warning", "info": "Note"}
    rows = [[words.get(f.severity, f.severity.title()), f.message, f.rule_id]
            for f in findings]
    number = counter.next_table()
    where, links = cite("Table", number)
    sec.blocks.append(Prose(
        f"xslope checks a model against what the selected analysis needs before "
        f"it runs. The findings in {where} are the ones that concern the "
        f"analyses in this report, in the checker's own words.", links=links))
    sec.blocks.append(Table(["Severity", "Finding", "Check"], rows,
                            "Model check findings", number))
    return sec


# ---------------------------------------------------------------------------
# The builder and the entry point
# ---------------------------------------------------------------------------

def planned_figures(slope_data, solutions, opts):
    """How many figures a build with these options will attempt.

    A report of one method renders three plots and a report of five renders
    eleven, so a caller that puts a wait cursor up has to be told which one it is
    on rather than left with a bar that never moves. Counted from the same option
    flags the sections are built from — and, for the force diagram, from the same
    refusal that decides whether a Calculations section exists at all — and
    checked against what a build actually produced, so the two cannot drift
    apart.
    """
    n = 0
    if opts["project_definition"] and opts["pd_figure"]:
        n += 1
    seep = seep_bundles(solutions) if opts["seep"] else []
    if seep:
        n += 1 if opts["seep_inputs_figure"] else 0
        n += 1 if opts["seep_kr_figure"] and _kr_materials(slope_data) else 0
        n += len(seep) * ((1 if opts["seep_mesh_figure"] else 0)
                          + (1 if opts["seep_flownet"] else 0))
    if opts["fem"] and fem_bundles(solutions):
        n += (1 if opts["fem_inputs_figure"] else 0)
        n += (1 if opts["fem_mesh_figure"] else 0)
    if opts["fem"]:
        for bundle in fem_bundles(solutions):
            if opts["fem_figure"]:
                n += len(FEM_PANELS)
            for kind, spec in DETAIL_KINDS.items():
                if opts[spec["option"]] and opts[spec["figure_option"]]:
                    n += len(_figured_members(
                        _detail_profiles(slope_data, bundle, kind)))
    if opts["lem"] and select_bundle(solutions, opts.get("method")) is not None:
        if (opts["lem_search"] and opts["lem_search_figure"]
                and search_bundle(solutions) is not None):
            n += 1
        per = ((1 if opts["lem_solution_figure"] else 0)
               + (1 if opts["lem_slice_table"] and opts["lem_slice_key"] else 0))
        for name in featured_methods(solutions, opts):
            n += per
            if _diagram_is_printed(slope_data, solutions, name, opts):
                n += 1
    return n


def _diagram_is_printed(slope_data, solutions, method, opts):
    """Will this method's block carry its force diagram?

    The diagram heads the Calculations section, so it prints exactly where the
    working does — which is a question about the model's forces, not about the
    options: a method whose factor of safety is no quotient of two sums gets the
    sentence that says so, and neither the working nor the picture of it. The
    same two calls the builder makes decide it, so the count cannot say one thing
    and the build another.
    """
    if not opts["lem_calculations"] or not force_diagram(method):
        return False
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            bundle, _note = detail_bundle(slope_data, solutions, method)
            if bundle is None:
                return False
            return calculation(slope_data, bundle, method) is not None
    except Exception:
        return False


def _progress_reporter(callback, total):
    """A one-argument ``progress(label)`` that counts its calls off against
    ``total`` and hands ``(done, total, label)`` to ``callback``. None when
    nobody is listening, so the builder can test for it once."""
    if not callable(callback):
        return None
    state = {"done": 0}

    def step(label):
        state["done"] += 1
        try:
            callback(state["done"], total, label)
        except Exception:
            pass                      # a progress line is never worth the report

    return step


def build_report(slope_data, solutions=None, options=None, figure_dir=None):
    """Build the content tree for an Analysis Report.

    Parameters
    ----------
    slope_data : dict
        The model, as :func:`xslope.fileio.load_slope_data` returns it.
    solutions : dict, optional
        What has been solved, keyed by engine. Each value is one bundle or a list
        of them, and each bundle is what that engine's Studio runner emits:

        ``"lem"``
            ``{"slice_df", "failure_surface", "results", "search", "method"}``,
            plus the method's name. One per method.
        ``"seep"``
            ``{"seep_data", "solution", "options"}``, where ``options["bc"]``
            names the boundary condition set. One per set solved.
        ``"fem"``
            ``{"fem_data", "solution", "FS", "analysis", "failure_solution"}``,
            where ``analysis`` is ``"ssrm"`` for a strength reduction run and
            ``"single"`` for one trial at a stated factor.

        An engine that is absent gets no section.
    options : dict, optional
        See :data:`DEFAULT_OPTIONS`. Unlisted keys are ignored.
    figure_dir : str, optional
        Where the rendered PNGs are written. Defaults to the current directory,
        so a caller that wants them kept must say where; :func:`generate_report`
        always names one. The directory is created when the first figure is
        written, and not at all when the figures are switched off.

    Returns
    -------
    Report
        The content tree, with every figure already rendered to disk.
    """
    opts = resolve_options(options)
    figure_dir = figure_dir or os.getcwd()
    counter = _Counter()

    # Long date, built rather than strftime'd: "%-d" is not portable to Windows.
    now = datetime.now()
    date = opts.get("date") or f"{now.strftime('%B')} {now.day}, {now.year}"
    meta = {
        "title": opts.get("title") or DEFAULT_OPTIONS["title"],
        "project_number": opts.get("project_number") or "",
        "organization": opts.get("organization") or "",
        "author": opts.get("author") or "",
        "date": date,
        "signature_lines": bool(opts.get("signature_lines")),
        "document_type": "Slope Stability Analysis Report",
    }
    report = Report(meta=meta)

    progress = _progress_reporter(opts.get("progress"),
                                  planned_figures(slope_data, solutions, opts))

    if opts["traceability"]:
        report.sections.append(_traceability_section(slope_data, solutions, opts))
    if opts["project_definition"]:
        report.sections.append(_project_definition_section(
            slope_data, opts, counter, figure_dir, progress))
    # Seepage stands ahead of the analyses that consume it: a stability analysis
    # whose materials take their pore pressure from seepage is run on the field
    # this section documents.
    if opts["seep"]:
        seep = _seep_section(slope_data, solutions, opts, counter, figure_dir,
                             progress)
        if seep is not None:
            report.sections.append(seep)
    if opts["lem"]:
        lem = _lem_section(slope_data, solutions, opts, counter, figure_dir,
                           progress)
        if lem is not None:
            report.sections.append(lem)
    if opts["fem"]:
        fem = _fem_section(slope_data, solutions, opts, counter, figure_dir,
                           progress)
        if fem is not None:
            report.sections.append(fem)
    if opts["model_checks"]:
        checks = _model_checks_section(slope_data, solutions, opts, counter)
        if checks is not None:
            report.sections.append(checks)
    return report


#: Formats :func:`generate_report` knows about. Only the enabled ones can be
#: rendered; the rest are the seams the later stages fill.
FORMATS = {
    "docx": {"label": "Word document (.docx)", "suffix": ".docx", "enabled": True},
    "pdf": {"label": "PDF (.pdf)", "suffix": ".pdf", "enabled": False},
    "latex": {"label": "LaTeX (.tex)", "suffix": ".tex", "enabled": False},
}


def generate_report(slope_data, solutions=None, options=None, path=None,
                    fmt=None, figure_dir=None):
    """Write an Analysis Report to ``path``.

    Parameters
    ----------
    slope_data, solutions, options
        As for :func:`build_report`.
    path : str
        The document to write. Its suffix selects the format unless ``fmt`` says
        otherwise.
    fmt : str, optional
        One of :data:`FORMATS`.
    figure_dir : str, optional
        Where the figures are kept. The figures are embedded in the document, so
        by default they are rendered into a temporary directory that is removed
        once the document is written — a report leaves one file behind, not a
        file and a folder. Naming a directory here keeps the PNGs in it.

    Returns
    -------
    (bool, dict or str)
        ``(True, {"path", "report", "figures"})`` on success, or ``(False,
        message)`` — the package's error convention. ``figures`` is the caption
        of every figure the document carries: the PNGs themselves are inside the
        document, and their files are gone unless ``figure_dir`` asked for them.
    """
    if not path:
        return False, "No output path was given for the report."
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".") or "docx").lower()
    if fmt == "tex":
        fmt = "latex"
    spec = FORMATS.get(fmt)
    if spec is None:
        return False, (f"Unknown report format {fmt!r}; expected one of "
                       f"{', '.join(sorted(FORMATS))}.")
    if not spec["enabled"]:
        return False, f"{spec['label']} reports are not available yet."

    keep_figures = figure_dir is not None
    if not keep_figures:
        figure_dir = tempfile.mkdtemp(prefix="xslope_report_")
    try:
        try:
            report = build_report(slope_data, solutions, options, figure_dir)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            return False, f"The report content could not be built: {exc}"

        try:
            from .report_docx import render_docx
            render_docx(report, path, template=(options or {}).get("template"))
        except Exception as exc:
            import traceback
            traceback.print_exc()
            return False, f"The report could not be written: {exc}"
    finally:
        if not keep_figures:
            shutil.rmtree(figure_dir, ignore_errors=True)

    return True, {"path": path, "report": report,
                  "figures": [f.caption for f in report.figures()]}
