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
unless the caller names a ``figure_dir`` to keep them in.
"""

from __future__ import annotations

import hashlib
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
    """

    text: str
    links: list = field(default_factory=list)

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
    """

    headers: list
    rows: list
    caption: str = ""
    number: int = 0
    landscape: bool = False
    legend: list = field(default_factory=list)
    bookmark: str = ""
    bold_rows: list = field(default_factory=list)

    def __post_init__(self):
        self.kind = "table"


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
        "This report was produced by xslope from the input file identified below. "
        "The file digest identifies the exact inputs the results were computed "
        "from, so the analysis can be reproduced or audited."))
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
            "is analysed dry, with zero pore pressure throughout."))
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
    """The distributed loads as entered, one row per block."""
    blocks = slope_data.get("dloads") or []
    if not blocks:
        return None
    lbl = _unit_labels(slope_data)
    su = f" ({lbl['stress']})" if lbl and lbl.get("stress") else ""
    lu = f" ({lbl['length']})" if lbl and lbl.get("length") else ""
    dirs = slope_data.get("dload_dirs") or []
    headers = ["#", f"x from{lu}", f"x to{lu}", f"Max pressure{su}", "Direction"]
    rows = []
    for i, blk in enumerate(blocks):
        xs = [_num(p.get("X")) for p in blk]
        ps = [_num(p.get("Normal")) or 0.0 for p in blk]
        xs = [x for x in xs if x is not None]
        rows.append([str(i + 1),
                     f"{min(xs):g}" if xs else "",
                     f"{max(xs):g}" if xs else "",
                     f"{max(ps):g}" if ps else "",
                     str(dirs[i] if i < len(dirs) else "normal")])
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
            "The model declares no unit system. Every quantity in this report is "
            "in the units the inputs were entered in; results are consistent with "
            "those units throughout.")
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
    if opts["pd_figure"]:
        path = os.path.join(figure_dir, "model.png")

        def draw(fig):
            from .plot import plot_inputs
            plot_inputs(slope_data, fig=fig, mode="shared", show_title=False,
                        frame="content", style=opts.get("style"))

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
        later = ["Trial failure surfaces"]
        if slope_data.get("mesh") is not None:
            later.append("the analysis mesh")
        text += (f" The figure below shows the model every analysis in this report "
                 f"is run on: {_join(shows)}. {_join(later)} are shown with the "
                 f"analyses that use them.")
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
        "The critical surface was located by automated search. The search refines "
        "a grid of trial surfaces until the factor of safety stops improving, so "
        "the reported surface is the lowest the search reached rather than a "
        "surface chosen by hand."))
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
        "The surface was analysed for rapid drawdown by the three-stage procedure "
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

    Read from the slice table's own circle centre where there is one, because that
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
            "between surfaces.")


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
                       "takes moments about a circle centre; the surface is "
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
#    safety before anything is printed (:data:`CALC_TOLERANCE`). A model whose
#    terms this module does not carry fails that comparison and gets no
#    calculation rather than a wrong one.
#
# The equations themselves are the ones on the method's documentation page, in
# that page's own symbols, so the report and the derivation can be read side by
# side. Where the code and the page differ the code wins and the page is the bug.
# ---------------------------------------------------------------------------

#: The Word bookmark placed on the slice table, and linked to from the
#: calculations — the per-slice terms of every sum are columns of that table.
SLICE_TABLE_BOOKMARK = "xslope_slice_table"

#: How closely the printed quotient must reproduce the solver's factor of safety
#: before the calculation is shown at all. This is the arithmetic identity, in
#: exact floats — not the printed-precision reproduction, which is looser and is
#: what the checks assert. A relative 1e-6 catches a missing term (which moves a
#: factor of safety in the third digit or worse) while leaving room for the
#: root-finder's own residual in the force-equilibrium methods.
CALC_TOLERANCE = 1e-6

#: Methods that take moments about the center of rotation, and methods that
#: balance horizontal forces over the whole sliding mass. Every method xslope
#: offers is in one list or the other.
MOMENT_METHODS = ("oms", "bishop")
FORCE_METHODS = ("janbu", "corps", "lowe", "spencer", "mprice")

#: Support features whose forces mobilize with the soil and so carry 1/F. They
#: put the factor of safety on both sides of every term they touch, which the
#: compact form cannot show; a model that uses one gets no calculation section.
PASSIVE_COLUMNS = ("p_pt", "pp_cx", "pp_cy", "pp_mx", "pp_my", "h_pile_pas")


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
    symbol the documentation gives it."""
    import numpy as np

    n = len(df)
    alpha = np.radians(df["alpha"].values.astype(float))
    c = df["c"].values.astype(float)
    if "c_suction" in df.columns:
        c = c + df["c_suction"].values.astype(float)
    return {
        "n": n,
        "alpha": alpha, "sin_a": np.sin(alpha), "cos_a": np.cos(alpha),
        "tan_phi": np.tan(np.radians(df["phi"].values.astype(float))),
        "c": c, "dl": _column(df, "dl"), "W": _column(df, "w"),
        "u": _column(df, "u"), "N": _column(df, "n_eff"),
        "D": _column(df, "dload"), "beta": np.radians(_column(df, "beta")),
        "kW": _column(df, "kw"), "T": _column(df, "t"),
        "P": _column(df, "p"), "pa_cx": _column(df, "pa_cx"),
        "pa_cy": _column(df, "pa_cy"), "pa_mx": _column(df, "pa_mx"),
        "pa_my": _column(df, "pa_my"),
        "H": _column(df, "h_pile"), "theta_p": _column(df, "theta_p"),
        "x_pile": _column(df, "x_pile"), "y_pile": _column(df, "y_pile"),
        "L": _column(df, "lload"), "ll_b": np.radians(_column(df, "ll_beta")),
        "ll_x": _column(df, "ll_x"), "ll_y": _column(df, "ll_y"),
        "d_x": _column(df, "d_x"), "d_y": _column(df, "d_y"),
        "y_cg": _column(df, "y_cg"), "y_t": _column(df, "y_t"),
        "x_c": _column(df, "x_c"), "y_cb": _column(df, "y_cb"),
    }


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


#: What each per-slice force in the equations is called when a model does not
#: have one, in the order the omission sentence lists them.
FEATURE_NAMES = (
    ("D", "distributed load"),
    ("kW", "seismic load"),
    ("T", "tension-crack water force"),
    ("P", "reinforcement crossing the failure surface"),
    ("H", "pile force"),
    ("L", "line load"),
)


def _absent_features(A):
    """The forces of the general equation this model does not carry at all."""
    out = []
    for key, name in FEATURE_NAMES:
        values = A[key] if key != "P" else (A["P"] + A["pa_cx"] + A["pa_cy"])
        if not _any(values):
            out.append(name)
    return out


def _moment_terms(df, A, right_facing):
    """``(resisting_terms, driving_terms)`` for a moment method, as equation (8a)
    of the Ordinary Method of Slices page writes them.

    Each list is ``[(sign, symbol, values, plain name), ...]`` and each moment is
    the signed sum of its terms' values. Passive support — capacity that
    mobilizes with the soil — contributes a resisting moment of its own, which is
    why the resisting side is a list rather than one term.
    """
    import numpy as np
    from .solve import _moment_arms

    Xo = float(df["xo"].iloc[0])
    Yo = float(df["yo"].iloc[0])
    x_r, a_S, a_N = _moment_arms(df, Xo, Yo, A["alpha"], right_facing)
    a_dx = (A["d_x"] - Xo) * (-1.0 if right_facing else 1.0)
    a_dy = Yo - A["d_y"]
    a_s = Yo - A["y_cg"]
    a_t = Yo - A["y_t"]
    pile_x = (A["x_pile"] - Xo) * (-1.0 if right_facing else 1.0)
    ll_x = (A["ll_x"] - Xo) * (-1.0 if right_facing else 1.0)

    resisting = [
        (+1, "(c·Δl + N'·tan φ)·a_S",
         (A["c"] * A["dl"] + A["N"] * A["tan_phi"]) * a_S, ""),
        (+1, "P_p·a_S", _column(df, "p_pt") * a_S, ""),
        (+1, "(P_p cos ψ·a_ry + P_p sin ψ·a_rx)",
         _column(df, "pp_cx") * Yo - _column(df, "pp_my")
         + _column(df, "pp_mx") - Xo * _column(df, "pp_cy"), ""),
        (+1, "(H_p cos θ_p·a_ey + H_p sin θ_p·a_ex)",
         _column(df, "h_pile_pas") * np.cos(A["theta_p"]) * (Yo - A["y_pile"])
         + _column(df, "h_pile_pas") * np.sin(A["theta_p"]) * pile_x, ""),
    ]
    terms = [
        (+1, "W·x_r", A["W"] * x_r, ""),
        (-1, "(N' + u·Δl)·a_N", (A["N"] + A["u"] * A["dl"]) * a_N, ""),
        (+1, "D·cos β·a_dx", A["D"] * np.cos(A["beta"]) * a_dx, "the distributed load"),
        (-1, "D·sin β·a_dy", A["D"] * np.sin(A["beta"]) * a_dy, ""),
        (+1, "kW·a_s", A["kW"] * a_s, "the seismic force"),
        (+1, "T·a_t", A["T"] * a_t, "the tension-crack water force"),
        (-1, "P·a_S", A["P"] * a_S, "the tangent reinforcement force"),
        (-1, "(P cos ψ·a_ry + P sin ψ·a_rx)",
         A["pa_cx"] * Yo - A["pa_my"] + A["pa_mx"] - Xo * A["pa_cy"],
         "the axial reinforcement force"),
        (-1, "(H cos θ_p·a_ey + H sin θ_p·a_ex)",
         (A["H"] - _column(df, "h_pile_pas")) * np.cos(A["theta_p"])
         * (Yo - A["y_pile"])
         + (A["H"] - _column(df, "h_pile_pas")) * np.sin(A["theta_p"]) * pile_x,
         "the pile force"),
        (+1, "L·a_fx", A["L"] * np.cos(A["ll_b"]) * ll_x, "the line load"),
        (-1, "L·a_fy", A["L"] * np.sin(A["ll_b"]) * (Yo - A["ll_y"]), ""),
    ]
    return resisting, terms


def _force_terms(A):
    """``(resisting, driving_terms)`` for a force-equilibrium method.

    Janbu's page writes the balance directly (its equation 7); the march the
    other four run is the force-equilibrium page's equations (6) and (10), and
    the same march under Spencer's and Morgenstern-Price's converged interslice
    forces. All five carry the horizontal component of the surface loads on the
    resisting side: a load written as (D·sin β, −D·cos β) has a SIGNED horizontal
    component in the sliding-direction frame, and a load normal to a slope face
    (β > 0) pushes into the slope.
    """
    import numpy as np

    load_sign = -1
    resisting = [(+1, "(c·Δl + N'·tan φ)·cos α",
                  (A["c"] * A["dl"] + A["N"] * A["tan_phi"]) * A["cos_a"], "")]
    terms = [
        (+1, "(N' + u·Δl)·sin α", (A["N"] + A["u"] * A["dl"]) * A["sin_a"], ""),
        (+1, "kW", A["kW"], "the seismic force"),
        (+1, "T", A["T"], "the tension-crack water force"),
        (load_sign, "D·sin β", A["D"] * np.sin(A["beta"]), "the distributed load"),
        (-1, "P cos ψ", A["P"] * A["cos_a"] + A["pa_cx"], "the reinforcement force"),
        (-1, "H cos θ_p", A["H"] * np.cos(A["theta_p"]), "the pile force"),
        (load_sign, "L cos δ", A["L"] * np.sin(A["ll_b"]), "the line load"),
    ]
    return resisting, terms


def _sum_notation(kept):
    """One side of the printed equation: one Σ per live term, signed."""
    out = ""
    for sign, symbol, _values, _name in kept:
        joiner = "" if not out else (" + " if sign > 0 else " − ")
        if not out and sign < 0:
            joiner = "−"
        out += f"{joiner}sum{{{symbol}}}"
    return out or "0"


def _spencer_state(df):
    """Spencer's per-slice ``Q`` and ``y_Q`` at the converged solution.

    Taken from the solver rather than rebuilt: its debug level writes ``Q``,
    ``y_q`` and the slice moment ``Mo`` into the table it is given, and the
    residuals the section prints are sums of exactly those numbers. The solve is
    repeated on a copy for that reason alone, and its console output — the whole
    point of a debug level — is swallowed.
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
    if not ok or "Q" not in work.columns:
        return None
    import numpy as np
    Q = work["Q"].values.astype(float)
    y_q = work["y_q"].values.astype(float)
    theta = np.radians(res["theta"])
    x_b = work["x_c"].values.astype(float)
    return {
        "Q": Q, "y_q": y_q, "FS": res["FS"], "theta": res["theta"],
        # Equations (27) and (28): the force and moment imbalances the solver
        # drove to zero, recomputed from the same per-slice values the table
        # prints so that the section's numbers and its columns are one thing.
        "R1": float(np.sum(Q)),
        "R2": float(np.sum(Q * (x_b * np.sin(theta) - y_q * np.cos(theta)))),
        "scale": float(np.sum(np.abs(Q))),
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


def calculation(slope_data, bundle, method):
    """The factor of safety calculation the report prints, or None.

    Returns a dict carrying the slice table with the per-slice contribution
    columns added, the sums, the equation notation, and everything the section
    needs to write itself. None means the calculation could not be shown for this
    model — see :data:`CALC_TOLERANCE` and :data:`PASSIVE_COLUMNS`.
    """
    import numpy as np

    method = str(method or "").lower()
    df, stage = _calc_state(bundle)
    results = bundle.get("results") or {}
    FS = _num(results.get("FS"))
    if df is None or not len(df) or FS is None or FS <= 0:
        return None
    if method not in MOMENT_METHODS + FORCE_METHODS:
        return None
    if "n_eff" not in df.columns or not np.all(np.isfinite(df["n_eff"].values)):
        return None
    if method in MOMENT_METHODS and "xo" not in df.columns:
        return None

    A = _calc_arrays(df)
    right_facing = bool(df.attrs.get(
        "right_facing", df["y_lb"].iat[0] > df["y_rb"].iat[-1]))
    if method in MOMENT_METHODS:
        # Passive capacity contributes a resisting MOMENT here, so the moment
        # methods can show it. In the force methods it divides by F on the
        # driving side, which no quotient can print; those models get no section.
        res_terms, terms = _moment_terms(df, A, right_facing)
        res_key, drv_key = "m_res", "m_drv"
    else:
        if any(_any(_column(df, name)) for name in PASSIVE_COLUMNS):
            return None
        res_terms, terms = _force_terms(A)
        res_key, drv_key = "f_res", "f_drv"

    scale = max([float(np.max(np.abs(values)))
                 for _s, _sym, values, _n in terms + res_terms])
    kept = _keep(terms, scale)
    kept_res = _keep(res_terms, scale)
    resisting = sum((sign * values for sign, _s, values, _n in kept_res),
                    np.zeros(A["n"]))
    driving = sum((sign * values for sign, _s, values, _n in kept),
                  np.zeros(A["n"]))
    sum_res = float(np.sum(resisting))
    sum_drv = float(np.sum(driving))
    if not np.isfinite(sum_res) or not np.isfinite(sum_drv) or sum_drv == 0:
        return None

    quotient = sum_res / sum_drv
    fo = _num(results.get("fo")) if method == "janbu" else None
    computed = quotient * fo if fo else quotient
    if abs(computed - FS) > CALC_TOLERANCE * FS:
        return None

    out = df.copy()
    out[res_key] = resisting
    out[drv_key] = driving

    spencer_state = None
    if method == "spencer":
        spencer_state = _spencer_state(df)
        if spencer_state is None:
            return None
        out["q_s"] = spencer_state["Q"]
        out["y_q"] = spencer_state["y_q"]

    return {
        "method": method, "slice_df": out, "FS": FS, "stage": stage,
        "resisting": sum_res, "driving": sum_drv, "quotient": quotient,
        "fo": fo, "spencer": spencer_state,
        "theta": _num(results.get("theta")) if method in ("corps", "lowe") else None,
        "lambda": _num(results.get("lambda")),
        "residuals": _mp_residuals_for(df, results) if method == "mprice" else None,
        "equation": (f"F = frac{{{_sum_notation(kept_res)}}}"
                     f"{{{_sum_notation(kept)}}}"),
        "kept": kept, "absent": _absent_features(A),
        "res_key": res_key, "drv_key": drv_key,
        "normal_force": _normal_force_equations(A, method),
    }


def _normal_force_equations(A, method):
    """The base-normal equations the iterative methods solve alongside the
    quotient — Bishop's equation (8) / Janbu's equation (6) — carrying the
    vertical components this model actually has.

    Janbu's page names the denominator m_α and Bishop's writes it out; each is
    printed the way its own page writes it.
    """
    import numpy as np

    numerator = "W"
    for values, symbol in (
            (A["D"] * np.cos(A["beta"]), " + D cos β"),
            (A["P"] * A["sin_a"] + A["pa_cy"], " − P sin ψ"),
            (A["H"] * np.sin(A["theta_p"]), " − H sin θ_p"),
            (A["L"] * np.cos(A["ll_b"]), " − L sin δ"),
            (A["u"] * A["dl"] * A["cos_a"], " − u·Δl·cos α")):
        if _any(values):
            numerator += symbol
    numerator += " − frac{c·Δl·sin α}{F}"
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
            "Every N' below is that value at the converged factor of safety, "
            "which is what makes the quotient close on itself: the sums formed "
            "from it return the F it was evaluated at."))
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
        blocks.append(Math(
            "Q = [−F_v·sin α − F_h·cos α − frac{c·Δl}{F} + "
            "(F_v·cos α − F_h·sin α + u·Δl)·frac{tan φ}{F}]·m_α"))
        blocks.append(Math(
            "m_α = frac{1}{cos (α − θ) + sin (α − θ)·frac{tan φ}{F}}"))
        blocks.append(Prose(
            f"The solution converged at F = {format_fs(state['FS'])} with "
            f"θ = {state['theta']:.2f} degrees, and at that pair both "
            f"equilibrium sums close on zero:"))
        blocks.append(Math(f"sum{{Q}} = {format_residual(state['R1'])}"))
        blocks.append(Math(
            f"sum{{Q·(x_b·sin θ − y_Q·cos θ)}} = "
            f"{format_residual(state['R2'])}"))
        blocks.append(Prose(
            f"Those are the residuals the solver converged on, and neither is "
            f"exactly zero — that is what a converged iteration leaves behind. "
            f"The per-slice values of Q and y_Q are columns Q_s and y_Q of the "
            f"slice table; added up as printed, rounded to a tenth of a force "
            f"unit, the force sum closes to within a thousandth of "
            f"{format_sum(state['scale'])}, the total the magnitudes of Q come "
            f"to. Since ΣQ = 0 the interslice forces cancel over the whole mass, "
            f"and the horizontal balance of the mass is the quotient below."))
    return blocks


def _calculations_section(calc, slope_data, table_number, unit_labels,
                          bookmark=SLICE_TABLE_BOOKMARK):
    """The Calculations section: what the method solves, the equation, the sums,
    the arithmetic, the factor of safety.

    ``bookmark`` is the slice table this method's per-slice terms are columns of.
    A report that documents several methods carries a slice table for each, and a
    cross-reference has to land on the one whose numbers the sums came from.
    """
    from .columns import BY_KEY, format_fs, format_residual, format_sum, unit_label

    method = calc["method"]
    label = method_label(method)
    sec = Section("Calculations")

    url = method_doc_url(method)
    intro = (f"The factor of safety reported above is worked through below. The "
             f"equation is the one the solver evaluates — the derivation "
             f"published for {label} in the XSLOPE documentation, in that page's "
             f"own symbols — and the numbers in it are the converged values.")
    if calc["absent"]:
        intro += (f" The model carries no {_join(calc['absent'])}; those terms, "
                  f"and any other that is zero on every slice, are dropped from "
                  f"the equation rather than printed as zeros.")
    if calc["stage"]:
        intro += (f" The governing stage is {calc['stage']}, and every number "
                  f"below is that stage's.")
    sec.blocks.append(Prose(intro, links=[(label, url)] if url else []))

    for block in _method_preamble(calc, method):
        sec.blocks.append(block)

    sec.blocks.append(Math(calc["equation"]))

    # --- the sums, and where their per-slice terms are ---
    res_col = BY_KEY.get(calc["res_key"])
    drv_col = BY_KEY.get(calc["drv_key"])
    unit = unit_label(res_col, unit_labels) if res_col is not None else ""
    n_slices = len(calc["slice_df"])
    if table_number and res_col is not None and drv_col is not None:
        where = f"Table {table_number}"
        in_units = f", both in {unit}" if unit else ""
        sec.blocks.append(Prose(
            f"Each slice's contribution to the two sums is a column of {where}: "
            f"{res_col.label} is the resisting term and {drv_col.label} is the "
            f"net driving term{in_units}. Summed over the {n_slices} slices:",
            links=[(where, f"#{bookmark}")]))
    else:
        sec.blocks.append(Prose(
            f"Summing the per-slice terms over the {n_slices} slices:"))

    sec.blocks.append(Math(
        f"F = frac{{{format_sum(calc['resisting'])}}}"
        f"{{{format_sum(calc['driving'])}}}"))
    if calc["fo"]:
        sec.blocks.append(Math(f"F = {format_fs(calc['quotient'])}"))
        sec.blocks.append(Prose(
            "Janbu's correction factor f_o compensates for the neglected "
            "interslice shear. It is read from the method's chart fit for this "
            "surface's depth-to-length ratio and the soil type, and multiplies "
            "the factor of safety above:"))
        sec.blocks.append(Math(
            f"F_corr = f_o·F = {format_sum(calc['fo'])}·"
            f"{format_fs(calc['quotient'])} = {format_fs(calc['FS'])}"))
    else:
        sec.blocks.append(Math(f"F = {format_fs(calc['FS'])}"))

    # --- Morgenstern-Price: the moment condition the force balance is solved
    # jointly with, and what each residual came out at ---
    residuals = calc.get("residuals")
    if residuals is not None:
        sec.blocks.append(Prose(
            "The moment of the whole sliding mass about the coordinate origin "
            "closes at the same (F, λ). Both residuals at the solution — the "
            "interslice force left at the far end of the march, and the moment "
            "sum — are what a converged iteration leaves behind:"))
        sec.blocks.append(Math(f"Z_n = {format_residual(residuals[0])}"))
        sec.blocks.append(Math(f"sum{{M_o}} = {format_residual(residuals[1])}"))
    return sec


def _method_section(slope_data, bundle, note, method, opts, counter, figure_dir,
                    progress=None):
    """Everything the report says about ONE method, under one heading.

    A report that documents several methods is this section repeated, and the
    heading names the method — a reader who opens the document at a slice table
    must never have to page backwards to learn whose numbers are on it. The
    figures and tables inside carry the method's name for the same reason, and
    take their numbers from the report-wide counter, so the sequence runs
    unbroken through however many methods there are.
    """
    label = method_label(method)
    sec = Section(f"Results — {label}")
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

    fs = _num(results.get("FS"))
    if fs is not None:
        sec.blocks.append(Prose(
            f"{label} gives a factor of safety of {fs:.3f} on the critical "
            f"surface." + (f" {note}" if note else "")))

    warns = results.get("warnings") or []
    if warns:
        sec.blocks.append(Prose(
            "The solution reported the following admissibility notes, which "
            "describe where the computed stresses depart from what the method "
            "assumes:"))
        sec.blocks.append(Bullets([str(w) for w in warns]))

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
            sec.blocks.append(Figure(
                fpath, f"Critical surface and slice forces — {label}",
                counter.next_figure(), source=f"{method} critical surface"))

    # --- rapid drawdown ---
    if opts["lem_rapid"]:
        rapid = _rapid_section(results, counter)
        if rapid is not None:
            sec.children.append(rapid)

    # --- slice table and calculations ---
    # The calculation is worked out first: it adds the per-slice terms of the
    # factor of safety equation to the table, which is how the section can point
    # at a column instead of walking the reader through fifteen slices.
    calc = None
    if opts["lem_calculations"]:
        try:
            calc = calculation(slope_data, bundle, method)
        except Exception:
            import traceback
            traceback.print_exc()
    table_df = calc["slice_df"] if calc is not None else slice_df

    bookmark = f"{SLICE_TABLE_BOOKMARK}_{method}"
    table_number = 0
    if opts["lem_slice_table"] and table_df is not None:
        from .columns import slice_table
        headers, rows, legend = slice_table(table_df, _unit_labels(slope_data))
        sub_tab = Section("Slice Table")
        sub_tab.blocks.append(Prose(
            f"The table below lists the slice geometry, forces and strengths for "
            f"the critical surface as solved by {label}. Forces are per unit "
            f"thickness of section."))

        # The key to the table's first column, immediately above it: the same
        # plot the results section carries, framed on the sliced mass alone and
        # with every slice labelled, so a row of the table can be found on the
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
            landscape=True, legend=legend, bookmark=bookmark))
        sec.children.append(sub_tab)

    if calc is not None:
        sec.children.append(_calculations_section(
            calc, slope_data, table_number, _unit_labels(slope_data), bookmark))
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
                      f"centre ({_fmt(c.get('Xo'), '{:g}')}, "
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

    # --- search: once, whichever method ran it ---
    if opts["lem_search"]:
        found = search_bundle(solutions)
        if found is not None:
            search = _search_section(slope_data, found, opts, counter, figure_dir,
                                     bundle_method(found) or methods[0], progress)
            if search is not None:
                sec.children.append(search)

    # --- every method's answer, once, ahead of the detail ---
    table = _fs_table(slope_data, solutions, opts, counter)
    if table is not None:
        sub_fs = Section("Factors of Safety")
        featured = _join([method_label(m) for m in methods])
        sub_fs.blocks.append(Prose(
            f"The table lists every limit equilibrium method xslope offers on the "
            f"critical surface. The methods that were run report their own answers; "
            f"the rest were solved on the same surface, so the comparison is between "
            f"methods rather than between surfaces. {featured} "
            f"{'is' if len(methods) == 1 else 'are'} set in bold and "
            f"{'is' if len(methods) == 1 else 'are'} reported in full below."))
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
        "xslope checks a model against what the selected analysis needs before it "
        "runs, and reports what it finds. The findings below are the ones that "
        "concern the analyses this report contains, in the checker's own words — a "
        "check that belongs to an engine the report does not document would be "
        "noise here. They are reported so that a reviewer sees them, not only the "
        "engineer who ran the analysis."))

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
    flags the sections are built from, and checked against what a build actually
    produced, so the two cannot drift apart.
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
        n += per * len(featured_methods(solutions, opts))
    return n


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
