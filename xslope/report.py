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
into a directory beside the document (a temporary one unless the caller names it),
so what is in the report is what was on screen.
"""

from __future__ import annotations

import hashlib
import os
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
    """A paragraph of running text."""

    text: str

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
    """

    path: str
    caption: str
    number: int = 0
    source: str = ""
    width_in: float = 6.5

    def __post_init__(self):
        self.kind = "figure"


@dataclass
class Table(Block):
    """Headers, rows of strings, and a caption.

    ``landscape`` asks the renderer for a rotated page — the slice table's whole
    reason for existing. ``legend`` is ``[(term, definition), ...]``, printed
    beneath the table.
    """

    headers: list
    rows: list
    caption: str = ""
    number: int = 0
    landscape: bool = False
    legend: list = field(default_factory=list)

    def __post_init__(self):
        self.kind = "table"


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
    "lem_rapid": True,
    "model_checks": False,            # opt-in (Norm: off by default)

    # --- what the report documents ---
    "method": None,                   # which solved method the detail follows
    "input_path": None,               # the .xlsx, for the traceability stamp
    "solved_at": None,                # datetime of the solve; None = now
    "style": None,                    # Studio's live display style
    "preflight": None,                # a PreflightReport captured at solve time
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


def select_bundle(solutions, method=None):
    """The bundle the report's critical-surface figure and slice table follow.

    The named method wins; otherwise the first bundle. Returns None when there is
    nothing to select from.
    """
    bundles = lem_bundles(solutions)
    if not bundles:
        return None
    want = str(method).lower() if method else None
    if want:
        for b in bundles:
            if bundle_method(b) == want:
                return b
    return bundles[0]


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
    take the whole report down)."""
    fig = _new_figure(opts["figsize"])
    try:
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
            "The model defines no groundwater and no external water: there are no "
            "piezometric lines, no seepage head boundaries, and no material takes "
            "a pore pressure. Pore pressures are zero on every slice and the "
            "section is analysed dry."))
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


def _project_definition_section(slope_data, opts, counter, figure_dir):
    sec = Section("Project Definition")
    feats = water_features(slope_data)

    # The figure is drawn first so the prose that introduces it can be written
    # from what was actually produced — a report never announces a figure that
    # the option switched off or the plot failed to render.
    figure = None
    if opts["pd_figure"]:
        path = os.path.join(figure_dir, "01_model.png")

        def draw(fig):
            from .plot import plot_inputs
            plot_inputs(slope_data, fig=fig, mode="shared", show_title=False,
                        frame="content", style=opts.get("style"))

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


def _search_section(slope_data, bundle, opts, counter, figure_dir, method):
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
        fpath = os.path.join(figure_dir, "02_search.png")

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


def _solve_for_summary(name, slice_df, rapid=False):
    """One method's answer on an already-built slice table, or None.

    ``solve_selected`` is the solver's own entry point, prints its result and
    returns the error string rather than a dict when a method does not converge;
    the printing is swallowed here because a report is not a console session, and
    a non-dict return is what "did not converge" means. The slice table is copied:
    the solvers write their working columns into it, and the reported method's own
    table must not pick up another method's arithmetic.
    """
    import contextlib
    import io

    from .solve import solve_selected
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = solve_selected(name, slice_df.copy(), rapid=rapid)
    except Exception:
        return None
    return out if isinstance(out, dict) else None


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

    rows = []
    for name in supported_methods():
        res = solved.get(name)
        if res is None:
            if method_surface_reason(name, family):
                rows.append([method_label(name), "not applicable",
                             "takes moments about a circle centre; the surface is "
                             "non-circular"])
                continue
            if base_df is None:
                continue
            res = _solve_for_summary(name, base_df, rapid)
        fs = _num((res or {}).get("FS"))
        if fs is None:
            rows.append([method_label(name), "did not converge", ""])
            continue
        rows.append([method_label(name), f"{fs:.3f}", _solution_parameters(res)])

    if not rows:
        return None
    return Table(["Method", "Factor of safety", "Solution parameters"], rows,
                 "Computed factors of safety", counter.next_table())


def _lem_section(slope_data, solutions, opts, counter, figure_dir):
    bundle = select_bundle(solutions, opts.get("method"))
    if bundle is None:
        return None
    method = bundle_method(bundle) or (opts.get("method") or DEFAULT_METHOD)
    results = bundle.get("results") or {}
    slice_df = bundle.get("slice_df")

    sec = Section("Limit Equilibrium Analysis")
    sec.blocks.append(Prose(
        "The stability of the section was evaluated by the method of slices. The "
        "factor of safety is the factor by which the available shear strength "
        "along the failure surface must be divided to bring the sliding mass into "
        "limiting equilibrium."))

    # --- engine inputs ---
    items = [("Method reported in detail", method_label(method))]
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

    # --- search ---
    if opts["lem_search"]:
        search = _search_section(slope_data, bundle, opts, counter, figure_dir, method)
        if search is not None:
            sec.children.append(search)

    # --- results ---
    sub_res = Section("Results")
    fs = _num(results.get("FS"))
    if fs is not None:
        sub_res.blocks.append(Prose(
            f"{method_label(method)} gives a factor of safety of {fs:.3f} on the "
            f"critical surface."))
    table = _fs_table(slope_data, solutions, opts, counter)
    if table is not None:
        sub_res.blocks.append(Prose(
            "The table lists every limit equilibrium method xslope offers. The "
            "methods that were run report their own answers; the rest were solved "
            "on the critical surface documented here, so the comparison is between "
            "methods rather than between surfaces."))
        sub_res.blocks.append(table)

    warns = results.get("warnings") or []
    if warns:
        sub_res.blocks.append(Prose(
            "The solution reported the following admissibility notes, which "
            "describe where the computed stresses depart from what the method "
            "assumes:"))
        sub_res.blocks.append(Bullets([str(w) for w in warns]))

    if opts["lem_solution_figure"]:
        fpath = os.path.join(figure_dir, f"03_solution_{method or 'lem'}.png")

        def draw(fig):
            from .plot import plot_solution
            plot_solution(slope_data, slice_df, bundle.get("failure_surface"),
                          results, fig=fig, show_title=False,
                          style=opts.get("style"))

        if slice_df is not None and _render(draw, fpath, opts):
            sub_res.blocks.append(Figure(
                fpath, f"Critical surface and slice forces — {method_label(method)}",
                counter.next_figure(), source=f"{method} critical surface"))
    sec.children.append(sub_res)

    # --- rapid drawdown ---
    if opts["lem_rapid"]:
        rapid = _rapid_section(results, counter)
        if rapid is not None:
            sec.children.append(rapid)

    # --- slice table ---
    if opts["lem_slice_table"] and slice_df is not None:
        from .columns import slice_table
        headers, rows, legend = slice_table(slice_df, _unit_labels(slope_data))
        sub_tab = Section("Slice Table")
        sub_tab.blocks.append(Prose(
            f"The table below lists the slice geometry, forces and strengths for "
            f"the critical surface as solved by {method_label(method)}. Forces are "
            f"per unit thickness of section."))
        sub_tab.blocks.append(Table(
            headers, rows,
            f"Slice data — {method_label(method)}", counter.next_table(),
            landscape=True, legend=legend))
        sec.children.append(sub_tab)

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
        always names one.

    Returns
    -------
    Report
        The content tree, with every figure already rendered to disk.
    """
    opts = resolve_options(options)
    figure_dir = figure_dir or os.getcwd()
    os.makedirs(figure_dir, exist_ok=True)
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

    if opts["traceability"]:
        report.sections.append(_traceability_section(slope_data, solutions, opts))
    if opts["project_definition"]:
        report.sections.append(
            _project_definition_section(slope_data, opts, counter, figure_dir))
    if opts["lem"]:
        lem = _lem_section(slope_data, solutions, opts, counter, figure_dir)
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
        Where the figures are written. Defaults to a ``<stem>_figures`` directory
        beside the document, so the PNGs survive for reuse.

    Returns
    -------
    (bool, dict or str)
        ``(True, {"path", "report", "figures"})`` on success, or ``(False,
        message)`` — the package's error convention.
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

    if figure_dir is None:
        stem = os.path.splitext(path)[0]
        figure_dir = f"{stem}_figures"
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

    return True, {"path": path, "report": report,
                  "figures": [f.path for f in report.figures()]}
