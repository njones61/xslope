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

    ``notation`` is the equation; ``label`` is an optional tag printed beside it.
    """

    notation: str
    label: str = ""

    def __post_init__(self):
        self.kind = "math"


@dataclass
class Section:
    """One numbered part of the document: a heading, its blocks, its children."""

    title: str
    blocks: list = field(default_factory=list)
    children: list = field(default_factory=list)

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
    "pd_materials": True,
    "pd_water": True,
    "pd_loads": True,
    "pd_reinforcement": True,
    "pd_units": True,
    "lem": True,
    "lem_search": True,
    "lem_search_figure": True,
    "lem_solution_figure": True,
    "lem_slice_table": True,
    "lem_slice_key": True,
    "lem_calculations": True,
    "lem_rapid": True,
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

def lem_bundles(solutions):
    """The LEM bundles in a ``solutions`` mapping, as a list.

    Accepts one bundle or a list of them, so a single-method run and a
    ``solve_all`` sweep are the same argument.
    """
    got = (solutions or {}).get("lem")
    if got is None:
        return []
    if isinstance(got, dict):
        return [got]
    return list(got)


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
        try:
            n_nodes = len(mesh["nodes"])
            n_elems = len(mesh["elements"])
            types = sorted({int(t) for t in mesh["element_types"]})
            names = {3: "tri3", 6: "tri6", 4: "quad4", 8: "quad8", 9: "quad9"}
            kinds = ", ".join(names.get(t, str(t)) for t in types)
            items.append(("Mesh", f"{n_nodes:,} nodes, {n_elems:,} elements ({kinds})"))
        except Exception:
            items.append(("Mesh", "present"))

    items.append(("Report generated", datetime.now().strftime("%Y-%m-%d %H:%M")))

    sec = Section("Traceability")
    sec.blocks.append(Prose(
        "The results in this report were computed by xslope from the input file "
        "below. Its SHA-256 digest identifies that file exactly."))
    sec.blocks.append(KeyValues(items))
    return sec


def _materials_table(slope_data, counter):
    """The materials table: referenced rows, populated columns."""
    rows_src = referenced_materials(slope_data)
    if not rows_src:
        return None
    mats = [m for _i, m in rows_src]
    lbl = _unit_labels(slope_data)

    def unit(key):
        return f" ({lbl[key]})" if lbl and lbl.get(key) else ""

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
    keep = [f for f in fields if f[3] or _populated(mats, f[0])]
    headers = [f[1] for f in keep]
    rows = []
    for idx, m in rows_src:
        row = []
        for key, _h, fmt, _always in keep:
            row.append(str(idx) if key == "__index" else fmt(m))
        rows.append(row)
    return Table(headers, rows, "Material properties", counter.next_table())


def _fmt(value, spec="{:.2f}"):
    n = _num(value)
    return "" if n is None else spec.format(n)


def _water_section(slope_data, feats):
    """What the model says about water, and which sheet says it.

    A model that carries no water says exactly that, in one sentence: the
    key-value rows below describe features, and a row about a feature the model
    does not have would be a statement about the analysis that is not true.
    """
    sub = Section("Water Conditions")
    if not feats["any"]:
        sub.blocks.append(Prose(
            "The model defines no groundwater and no external water; the section "
            "is analyzed dry, with zero pore pressure throughout."))
        return sub

    from .water import water_line_for_stage, water_loads_mode

    items = []
    lbl = _unit_labels(slope_data)
    gw = _num(slope_data.get("gamma_water"))
    if gw is not None:
        suffix = f" {lbl['unit_weight']}" if lbl and lbl.get("unit_weight") else ""
        items.append(("Unit weight of water", f"{gw:g}{suffix}"))

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

    for stage, key, name in ((1, "seepage_bc", "Seepage boundaries"),
                             (2, "seepage_bc2", "Seepage boundaries (set 2)")):
        heads = ((slope_data.get(key) or {}).get("specified_heads") or [])
        if heads:
            items.append((name, f"{len(heads)} specified-head boundary "
                                f"block(s)"))
        if stage in feats["surfaces"]:
            line = water_line_for_stage(slope_data, stage=stage)
            items.append((f"Water surface (stage {stage})",
                          f"from {line['source']}"))

    if feats["pore"]:
        items.append(("Pore pressure by material", ", ".join(feats["pore"])))

    sub.blocks.append(KeyValues(items))
    return sub


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
    from .units import normalize_unit_system
    system = normalize_unit_system(slope_data.get("unit_system"))
    if system is None:
        return Prose(
            "The model declares no unit system. Every quantity is in the units "
            "the inputs were entered in, and results are consistent with those "
            "units throughout.")
    lbl = _unit_labels(slope_data) or {}
    name = "SI" if system == "si" else "US customary"
    return Prose(
        f"All quantities are in {name} units: lengths in {lbl.get('length', '')}, "
        f"stresses and pressures in {lbl.get('stress', '')}, unit weights in "
        f"{lbl.get('unit_weight', '')}, and forces per unit thickness of section "
        f"in {lbl.get('force_per_len', '')}. Angles are in degrees. The analysis "
        f"is two-dimensional: every force is per unit thickness normal to the "
        f"section.")


def _project_definition_section(slope_data, opts, counter, figure_dir,
                                progress=None):
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
        later = ["Trial failure surfaces"]
        if slope_data.get("mesh") is not None:
            later.append("the analysis mesh")
        text += (f" Every analysis in this report is run on the model below: "
                 f"{_join(shows)}. {_join(later)} are shown with the analyses "
                 f"that use them.")
    sec.blocks.append(Prose(text))

    # The units statement leads: a reader meets the numbers knowing what they are
    # in, rather than finding out at the end of the section.
    if opts["pd_units"]:
        sec.blocks.append(_units_prose(slope_data))

    if figure is not None:
        sec.blocks.append(figure)

    if opts["pd_materials"]:
        sub = Section("Materials")
        table = _materials_table(slope_data, counter)
        if table is not None:
            sub.blocks.append(table)
        else:
            sub.blocks.append(Prose("The model defines no materials."))
        sec.children.append(sub)

    if opts["pd_water"]:
        sec.children.append(_water_section(slope_data, feats))

    if opts["pd_loads"]:
        table = _loads_table(slope_data, counter)
        sub = Section("Loads")
        if table is not None:
            sub.blocks.append(table)
        elif feats["surfaces"]:
            sub.blocks.append(Prose(
                "The model carries no distributed loads entered by hand. Any "
                "water standing on the section is measured by the engine from "
                "the water surface and applied as a distributed load."))
        else:
            sub.blocks.append(Prose("The model carries no distributed loads."))
        k = _num(slope_data.get("k_seismic"))
        if k:
            sub.blocks.append(Prose(
                f"A pseudo-static seismic coefficient of k = {k:g} is applied. In "
                f"the limit equilibrium analysis it acts horizontally toward the "
                f"toe on every slice."))
        sec.children.append(sub)

    # Reinforcement and piles are separate features with separate tables, and a
    # section for one is not a heading the other hides under.
    if opts["pd_reinforcement"]:
        reinf = _reinforcement_table(slope_data, counter)
        if reinf is not None:
            sec.children.append(Section("Reinforcement", [reinf]))
        piles = _piles_table(slope_data, counter)
        if piles is not None:
            sec.children.append(Section("Piles", [piles]))

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

    sub.blocks.append(Prose(
        "The critical surface was located by automated search, which refines a "
        "grid of trial surfaces until the factor of safety stops improving. The "
        "reported surface is the lowest the search reached."))
    sub.blocks.append(KeyValues(items))

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
            sub.blocks.append(Figure(
                fpath, "Trial surfaces evaluated by the search, with the critical "
                       "surface highlighted", counter.next_figure(),
                source=f"{method} search"))
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
    """EVERY method xslope offers, and its factor of safety on the critical
    surface — not only the methods the caller happened to run.

    A summary that lists three methods because three were solved invites the
    question the report exists to close: what would the others have said? The
    methods that were run report their own answers; the rest are solved here, on
    the critical surface the report documents, which costs milliseconds on a slice
    table that already exists. A method that cannot apply to this surface family,
    and a method that does not converge on it, each say so in a row of their own
    rather than being dropped — a missing row reads as an answer withheld.

    The methods the report goes on to document in detail are set in BOLD. The
    table is the one place every method stands side by side, and the reader who
    finds the report's own methods there has the comparison and the answer in one
    glance.
    """
    from .preflight import method_surface_reason

    solved = {}
    for b in lem_bundles(solutions):
        name = bundle_method(b)
        if name and name not in solved:
            solved[name] = b.get("results") or {}
    bundle = select_bundle(solutions, opts.get("method")) or {}
    base_df = bundle.get("slice_df")
    rapid = "stage1_FS" in (bundle.get("results") or {})
    family = _surface_family(base_df, slope_data)
    featured = set(featured_methods(solutions, opts))

    rows, bold = [], []
    for name in supported_methods():
        res = solved.get(name)
        row = None
        if res is None:
            if method_surface_reason(name, family):
                row = [method_label(name), "not applicable",
                       "takes moments about a circle center; the surface is "
                       "non-circular"]
            elif base_df is None:
                continue
            else:
                res = _solve_for_summary(name, base_df, rapid)
        if row is None:
            fs = _num((res or {}).get("FS"))
            row = ([method_label(name), "did not converge", ""] if fs is None
                   else [method_label(name), f"{fs:.3f}", _solution_parameters(res)])
        if name in featured:
            bold.append(len(rows))
        rows.append(row)

    if not rows:
        return None
    return Table(["Method", "Factor of safety", "Solution parameters"], rows,
                 "Computed factors of safety", counter.next_table(),
                 bold_rows=bold)


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
             "spencer_h", "spencer_v", "normal")


@dataclass(frozen=True)
class NotApplicable:
    """Why a force contributes nothing to one of the equations.

    Not a blank: a reason, in the mechanics of the equation it is absent from.
    """

    reason: str


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
    """

    name: str
    group: str
    meaning: str
    rank: float = None


#: The nomenclature's groups, in the order it prints them.
SYMBOL_GROUPS = ("angle", "arm", "letter", "spencer")


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
        normal=(Term(+1, "D cos β", lambda C: C.A["D"] * C.cos("beta")),),
    ),
    ForceTerm(
        key="kW",
        columns=("kw",),
        arrays=(("kW", "kw", False), ("y_cg", "y_cg", False)),
        symbols=(Symbol("a_s", "arm",
                        "moment arm of the seismic force, taken at the slice "
                        "center of gravity"),),
        feature="seismic load", passive=False,
        moment_res=NotApplicable("the seismic force drives"),
        moment_drv=(Term(+1, "kW·a_s", lambda C: C.A["kW"] * C.arms["a_s"]),),
        force_res=NotApplicable("the seismic force drives"),
        force_drv=(Term(+1, "kW", lambda C: C.A["kW"]),),
        spencer_h=(Term(-1, "kW", lambda C: C.A["kW"]),),
        spencer_v=_NO_VERTICAL_COMPONENT,
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
                 Symbol("V", "spencer",
                        "resultant force of the water in a tension crack — "
                        "column T_c of the slice table, and written T below",
                        rank=2)),
        feature="tension-crack water force", passive=False,
        moment_res=NotApplicable("the water in the crack drives"),
        moment_drv=(Term(+1, "T·a_t", lambda C: C.A["T"] * C.arms["a_t"]),),
        force_res=NotApplicable("the water in the crack drives"),
        force_drv=(Term(+1, "T", lambda C: C.A["T"]),),
        spencer_h=(Term(-1, "V", lambda C: C.A["T"]),),
        spencer_v=_NO_VERTICAL_COMPONENT,
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
        normal=(Term(-1, "L sin δ", lambda C: C.A["L"] * C.cos("ll_b")),),
    ),
    ForceTerm(
        # In the published equations (1) and (2) and in no xslope model: the
        # solver does not simulate a shear force on the top of the slice, and a
        # force no model can carry is not an omission the report can report.
        key="T_top",
        columns=(), arrays=(), symbols=(),
        feature="", passive=False,
        moment_res=NotApplicable("xslope does not simulate it"),
        moment_drv=NotApplicable("xslope does not simulate it"),
        force_res=NotApplicable("xslope does not simulate it"),
        force_drv=NotApplicable("xslope does not simulate it"),
        spencer_h=NotApplicable(
            "the shear on the top of the slice is in equation (1) as published "
            "and is not simulated"),
        spencer_v=NotApplicable(
            "the shear on the top of the slice is in equation (2) as published "
            "and is not simulated"),
        normal=NotApplicable("xslope does not simulate it"),
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
    "λ": "scaling factor on the interslice force function f(x)",
    "f(x)": "the interslice force function; tan θ = λ·f(x)",
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
    rows = []
    for index, term in enumerate(FORCE_TERMS):
        got = getattr(term, consumer)
        if isinstance(got, NotApplicable):
            continue
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


def _spencer_force_sums(A):
    """Equations (1) and (2) of the Spencer page — the two force sums Q is built
    from — and the symbols they need defining.

    Transcribed from the derivation the section links to, in that page's own
    symbols, carrying only the terms this model has: the same convention the rest
    of the section follows, so that a reader meets no column of zeros and no term
    for a force the model does not apply.

    That page follows UTEXAS, in which P is the distributed-load resultant, R the
    reinforcement force and V the tension-crack water force — three letters the
    other derivations here write D, P and T. The returned symbol definitions are
    what let the two be read side by side; they are handed to
    :func:`equation_symbols`, which prefers them to the column registry.

    Returns ``(lines, symbols)`` — ``lines`` is ``[(notation, label), ...]``.
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

    symbols = {}
    if any(sym.startswith("P ") for sym in printed):
        # P is the one letter the section would otherwise use for two different
        # forces: the load on top of the slice here, the reinforcement below.
        dload, reinf = BY_KEY["dload"].label, BY_KEY["p"].label
        meaning = (f"resultant of the distributed load on the top of the slice, "
                   f"in the symbols of equations (1) and (2) — column {dload} of "
                   f"the slice table")
        if any(sym.startswith("R ") for sym in printed):
            meaning += (f", and written {dload} below, where {reinf} is instead "
                        f"the reinforcement force")
        symbols["P"] = meaning
    return ([(f"F_h = {_signed_notation(kept_h)}", "(1)"),
             (f"F_v = {_signed_notation(kept_v)}", "(2)")], symbols)


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
SPENCER_EQUILIBRIUM = (
    ("R_1 = sum{Q}", "(27)"),
    ("R_2 = sum{Q·(x_b·sin θ − y_Q·cos θ)}", "(28)"),
)


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
        return None, (
            f"Evaluated on this model's converged values the equation gives a "
            f"factor of safety of {computed:.6f}, and the solution reports "
            f"{FS:.6f}. A quotient that does not return the solution is not the "
            f"working behind it, and none is printed.")

    out = df.copy()
    out[res_key] = resisting
    out[drv_key] = driving

    return {
        "method": method, "slice_df": out, "FS": FS, "stage": stage,
        "resisting": sum_res, "driving": sum_drv, "quotient": quotient,
        "fo": fo, "spencer": None,
        "theta": _num(results.get("theta")) if method in ("corps", "lowe") else None,
        "lambda": _num(results.get("lambda")),
        "residuals": _mp_residuals_for(df, results) if method == "mprice" else None,
        "equation": (f"F = frac{{{_sum_notation(kept_res)}}}"
                     f"{{{_sum_notation(kept)}}}"),
        "kept": kept, "absent": _absent_features(A, df),
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

    force_sums, force_symbols = _spencer_force_sums(A)
    out = df.copy()
    out["F_h"] = state["F_h"]
    out["F_v"] = state["F_v"]
    out["q_s"] = state["Q"]
    out["y_q"] = state["y_q"]
    return {
        "method": "spencer", "slice_df": out, "FS": FS, "stage": stage,
        "spencer": state, "absent": _absent_features(A, df),
        "equilibrium": list(SPENCER_EQUILIBRIUM),
        "force_sums": force_sums, "symbols": force_symbols,
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


def _method_preamble(calc, method):
    """What this method solves, ahead of the equation — one short block, in the
    method's own terms."""
    from .columns import format_fs, format_residual, format_sum

    blocks = []
    if method in ("bishop", "janbu"):
        blocks.append(Prose(
            "The base normal force N' comes from vertical equilibrium of the "
            "slice and depends on the factor of safety itself, so it and the "
            "quotient are solved together by iteration:"))
        blocks.extend(Math(line) for line in calc["normal_force"])
        blocks.append(Prose(
            "Every N' below is that value at the converged factor of safety, so "
            "the sums formed from it return the F it was evaluated at."))
    elif method in ("corps", "lowe"):
        theta = calc.get("theta")
        stated = (f", averaging {theta:.2f} degrees on this surface"
                  if theta is not None else "")
        blocks.append(Prose(
            f"The interslice forces are taken at an inclination θ that the "
            f"method's own convention fixes from the geometry{stated}. The "
            f"solver marches the slices at a trial factor of safety and adjusts "
            f"it until the interslice force left over at the far end is zero; "
            f"at that value the interslice forces cancel over the whole sliding "
            f"mass, and the horizontal balance of the mass is the quotient "
            f"below."))
    elif method == "mprice":
        lam = calc.get("lambda")
        blocks.append(Prose(
            f"The interslice inclination varies along the surface as "
            f"tan θ = λ·f(x)"
            f"{f', with λ = {lam:.4f} at the solution' if lam is not None else ''}"
            f". λ and F are solved together so that force and moment "
            f"equilibrium of the whole sliding mass are satisfied at once. With "
            f"the interslice forces cancelling in the sum, the horizontal "
            f"balance is the quotient below."))
    elif method == "spencer":
        state = calc["spencer"]
        blocks.append(Prose(
            "Spencer's method lumps the interslice forces on each slice into a "
            "single resultant Q acting at the constant inclination θ, so that "
            "force and moment equilibrium of the whole sliding mass are two "
            "equations in the two unknowns F and θ. F_h and F_v are the sums of "
            "the forces on the slice other than the base normal, the base shear "
            "and the interslice forces:"))
        blocks.extend(Math(line, label) for line, label in calc["force_sums"])
        blocks.append(Prose(
            "Q on each slice follows from them and from the strength mobilized "
            "on its base:"))
        blocks.append(Math(
            "Q = [−F_v·sin α − F_h·cos α − frac{c·Δl}{F} + "
            "(F_v·cos α − F_h·sin α + u·Δl)·frac{tan φ}{F}]·m_α", "(23)"))
        blocks.append(Math(
            "m_α = frac{1}{cos (α − θ) + sin (α − θ)·frac{tan φ}{F}}", "(24)"))
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
    if diagram:
        sec.blocks.append(Figure(
            diagram, f"Forces on a slice — {label}",
            counter.next_figure() if counter is not None else 0,
            source=f"{method} force diagram",
            width_in=FORCE_DIAGRAM_WIDTH_IN))

    preamble = _method_preamble(calc, method)

    url = method_doc_url(method)
    # Singular where the section closes on one equation, which is every method
    # but Spencer's: its answer is the root of two, and it prints the working
    # that gets there.
    if method in WHOLE_MASS_BALANCE_METHODS:
        intro = (f"The equation below is the horizontal force balance of the "
                 f"whole sliding mass at the converged solution, in the symbols "
                 f"of the XSLOPE documentation; the numbers in it are the "
                 f"converged values. The slice-by-slice march that reaches that "
                 f"solution is the derivation published for {label}.")
    elif calc["equation"]:
        intro = (f"The equation below is the one the solver evaluates, in the "
                 f"symbols of the derivation published for {label} in the "
                 f"XSLOPE documentation; the numbers in it are the converged "
                 f"values.")
    else:
        intro = (f"The equations below are the ones the solver evaluates, in "
                 f"the symbols of the derivation published for {label} in the "
                 f"XSLOPE documentation; the numbers in them are the converged "
                 f"values.")
    if calc["absent"]:
        intro += (f" The model carries no {_join(calc['absent'])}; those terms, "
                  f"and any other that is zero on every slice, are dropped "
                  f"rather than printed as zeros.")
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
    equation = []
    if calc["equation"]:
        equation.append(Math(calc["equation"]))
    else:
        equation.append(Prose(
            "Substituting Q and y_Q into the equilibrium of the whole sliding "
            "mass gives two equations in the two unknowns. R_1 is the force "
            "imbalance and R_2 the moment imbalance, and the solution is the "
            "pair (F, θ) at which both are zero:"))
        equation.extend(Math(line, lab) for line, lab in calc["equilibrium"])
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
        nomenclature.append(Prose(
            f"The symbols above, in the order they appear. Those that are "
            f"columns of {where} carry a value for every slice.",
            links=([(where, f"#{bookmark}")] if table_number else [])))
        nomenclature.append(Table(
            ["Symbol", "Meaning"], [[s, m] for s, m in symbols],
            f"Nomenclature — {label}",
            counter.next_table() if counter is not None else 0))

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
            "interslice shear. It is read from the method's chart fit for this "
            "surface's depth-to-length ratio and the soil type, and multiplies "
            "the factor of safety above:"))
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
    fs = _num(results.get("FS"))
    if fs is not None:
        res.blocks.append(Prose(
            f"{label} gives a factor of safety of {fs:.3f} on the critical "
            f"surface." + (f" {note}" if note else ""),
            bold=[f"{fs:.3f}"]))

    warns = results.get("warnings") or []
    if warns:
        res.blocks.append(Prose(
            "The solution reported the following admissibility notes, where the "
            "computed stresses depart from what the method assumes:"))
        res.blocks.append(Bullets([str(w) for w in warns]))

    if opts["lem_solution_figure"] and slice_df is not None:
        fpath = os.path.join(figure_dir, f"solution_{method or 'lem'}.png")

        def draw(fig):
            from .plot import plot_solution
            plot_solution(slope_data, slice_df, bundle.get("failure_surface"),
                          results, fig=fig, show_title=False,
                          style=opts.get("style"))

        if progress:
            progress(f"the critical surface — {label}")
        if _render(draw, fpath, opts):
            res.blocks.append(Figure(
                fpath, f"Critical surface and slice forces — {label}",
                counter.next_figure(), source=f"{method} critical surface"))
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
        sub_tab.blocks.append(Prose(
            f"Slice geometry, forces and strengths for the critical surface as "
            f"solved by {label}. Forces are per unit thickness of section."))

        # The key to the table's first column, immediately above it: the same
        # plot the results section carries, framed on the sliced mass alone and
        # with every slice labeled, so a row of the table can be found on the
        # section it describes.
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
                sub_tab.blocks.append(Figure(
                    kpath, f"Slice numbering for the table below — {label}",
                    counter.next_figure(), source=f"{method} slice key",
                    width_in=0, landscape=True))

        table_number = counter.next_table()
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
    items = [(("Method" if len(methods) == 1 else "Methods") + " reported in detail",
              _join([method_label(m) for m in methods]))]
    if slice_df is not None:
        items.append(("Slices", str(len(slice_df))))
    circular = bool(slope_data.get("circular"))
    items.append(("Surface family", "circular" if circular else "non-circular"))
    if circular and slope_data.get("circles"):
        c = slope_data["circles"][0]
        items.append(("Starting circle",
                      f"center ({_fmt(c.get('Xo'), '{:g}')}, "
                      f"{_fmt(c.get('Yo'), '{:g}')}), R = {_fmt(c.get('R'), '{:g}')}"))
    elif slope_data.get("non_circ"):
        items.append(("Non-circular surface",
                      f"{len(slope_data['non_circ'])} defining points"))
    k = _num(slope_data.get("k_seismic"))
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

    # --- every method's answer, once, ahead of the detail ---
    #
    # The search does NOT belong here. It finds the critical surface for one
    # METHOD, and two methods searched separately settle on different surfaces;
    # each method's block carries its own search.
    table = _fs_table(slope_data, solutions, opts, counter)
    if table is not None:
        sub_fs = Section("Factors of Safety")
        featured = _join([method_label(m) for m in methods])
        # Which surface the filled-in rows were computed on has to be said, not
        # implied: a method that was RUN reports its own critical surface, and a
        # method that was not is solved on the surface named here.
        base = select_bundle(solutions, opts.get("method"))
        base_label = method_label(bundle_method(base)) if base else ""
        on_surface = (f"the critical surface {base_label} found"
                      if base_label else "the critical surface")
        run = [method_label(m) for m in solved_methods(solutions)]
        sub_fs.blocks.append(Prose(
            f"Every limit equilibrium method xslope offers is listed below. "
            f"{_join(run) or 'The method that was run'} reported "
            f"{'their' if len(run) > 1 else 'its'} own answer on the surface "
            f"{'each' if len(run) > 1 else 'it'} searched; every other method was "
            f"solved on {on_surface}, so those rows compare methods rather than "
            f"surfaces. {featured} "
            f"{'is' if len(methods) == 1 else 'are'} set in bold and reported in "
            f"full below, with the search that found "
            f"{'its' if len(methods) == 1 else 'each'} surface."))
        sub_fs.blocks.append(table)
        sec.children.append(sub_fs)

    # --- one full detail block per featured method ---
    for name in methods:
        detail, note = detail_bundle(slope_data, solutions, name)
        sec.children.append(_method_section(
            slope_data, detail, note, name, opts, counter, figure_dir, progress))

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
    sec.blocks.append(Prose(
        "xslope checks a model against what the selected analysis needs before "
        "it runs. The findings below are the ones that concern the analyses in "
        "this report, in the checker's own words."))

    findings = relevant_findings(getattr(report, "findings", []) or [],
                                 report_analyses(solutions, opts))
    if not findings:
        sec.blocks.append(Prose(
            "The model checks raised no findings for the analyses in this report."))
        return sec

    order = {"error": 0, "warning": 1, "info": 2}
    findings.sort(key=lambda f: order.get(f.severity, 3))
    words = {"error": "Error", "warning": "Warning", "info": "Note"}
    rows = [[words.get(f.severity, f.severity.title()), f.message, f.rule_id]
            for f in findings]
    sec.blocks.append(Table(["Severity", "Finding", "Check"], rows,
                            "Model check findings", counter.next_table()))
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
        ``{"lem": bundle}`` or ``{"lem": [bundle, ...]}``. A bundle is
        ``{"slice_df", "failure_surface", "results", "search", "method"}``, which
        is what Studio's LEM runner emits plus the method's name.
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
    if opts["lem"]:
        lem = _lem_section(slope_data, solutions, opts, counter, figure_dir,
                           progress)
        if lem is not None:
            report.sections.append(lem)
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
