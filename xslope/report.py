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

#: Display names for the limit-equilibrium methods.
METHOD_NAMES = {
    "oms": "Ordinary Method of Slices",
    "bishop": "Bishop's Simplified Method",
    "janbu": "Janbu's Simplified Method",
    "spencer": "Spencer's Method",
    "corps": "Corps of Engineers Method",
    "lowe": "Lowe and Karafiath Method",
    "mprice": "Morgenstern-Price Method",
}


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
    "model_checks": True,

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


def _water_block(slope_data):
    """What the model says about water, and which sheet says it."""
    from .water import water_line_for_stage, water_loads_mode

    items = []
    lbl = _unit_labels(slope_data)
    gw = _num(slope_data.get("gamma_water"))
    if gw is not None:
        suffix = f" {lbl['unit_weight']}" if lbl and lbl.get("unit_weight") else ""
        items.append(("Unit weight of water", f"{gw:g}{suffix}"))

    mode = water_loads_mode(slope_data)
    items.append(("Water loads", (
        "derived by the engine from the model's own water surface (automatic)"
        if mode == "auto" else "taken from the distributed loads as entered (manual)")))

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
        line = water_line_for_stage(slope_data, stage=stage)
        if line.get("points"):
            items.append((f"Water surface (stage {stage})",
                          f"from {line['source']}"))

    methods = sorted({str((m.get("u") or "none")).lower()
                      for _i, m in referenced_materials(slope_data)})
    if methods:
        items.append(("Pore pressure by material", ", ".join(methods)))

    if not items:
        return None
    return KeyValues(items)


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
               f"L_p1{lu}", f"L_p2{lu}", f"Spacing{lu}", "Direction", "Application"]
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
               f"D{lu}", f"Spacing{lu}", "Application"]
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

    geometry = ("material zone polygons"
                if slope_data.get("polygons") and not slope_data.get("profile_lines")
                else "profile lines")
    mats = referenced_materials(slope_data)
    sec.blocks.append(Prose(
        f"The section is defined by {len(mats)} material "
        f"{'zone' if len(mats) == 1 else 'zones'} described with {geometry}. "
        f"The figure below shows the model every analysis in this report is run "
        f"on: the geometry and materials, the water surfaces, and any loads, "
        f"reinforcement and piles the model carries. Trial failure surfaces and "
        f"analysis meshes are shown with the analyses that use them."))

    if opts["pd_figure"]:
        path = os.path.join(figure_dir, "01_model.png")

        def draw(fig):
            from .plot import plot_inputs
            plot_inputs(slope_data, fig=fig, mode="shared", show_title=False,
                        frame="content", style=opts.get("style"))

        if _render(draw, path, opts):
            sec.blocks.append(Figure(path, "Analysis model", counter.next_figure(),
                                     source="shared model"))

    if opts["pd_materials"]:
        sub = Section("Materials")
        table = _materials_table(slope_data, counter)
        if table is not None:
            sub.blocks.append(table)
        else:
            sub.blocks.append(Prose("The model defines no materials."))
        sec.children.append(sub)

    if opts["pd_water"]:
        sub = Section("Water Conditions")
        block = _water_block(slope_data)
        if block is not None:
            sub.blocks.append(block)
        else:
            sub.blocks.append(Prose("The model defines no water conditions; the "
                                    "section is analysed dry."))
        sec.children.append(sub)

    if opts["pd_loads"]:
        table = _loads_table(slope_data, counter)
        sub = Section("Loads")
        if table is not None:
            sub.blocks.append(table)
        else:
            sub.blocks.append(Prose(
                "The model carries no distributed loads entered by hand. Any "
                "water standing on the section is measured by the engine from "
                "the water surface and applied as a distributed load."))
        k = _num(slope_data.get("k_seismic"))
        if k:
            sub.blocks.append(Prose(
                f"A pseudo-static seismic coefficient of k = {k:g} is applied. In "
                f"the limit equilibrium analysis it acts horizontally toward the "
                f"toe on every slice."))
        sec.children.append(sub)

    if opts["pd_reinforcement"]:
        reinf = _reinforcement_table(slope_data, counter)
        piles = _piles_table(slope_data, counter)
        if reinf is not None or piles is not None:
            sub = Section("Reinforcement and Piles")
            if reinf is not None:
                sub.blocks.append(reinf)
            if piles is not None:
                sub.blocks.append(piles)
            sec.children.append(sub)

    if opts["pd_units"]:
        sub = Section("Units")
        sub.blocks.append(_units_prose(slope_data))
        sec.children.append(sub)

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


def _fs_table(solutions, counter):
    """Every solved method and its factor of safety. Independent of which method
    the report's detail follows, so two reports of the same run always agree on
    the answers."""
    rows = []
    for b in lem_bundles(solutions):
        res = b.get("results") or {}
        fs = _num(res.get("FS"))
        name = bundle_method(b)
        extra = []
        if _num(res.get("theta")) is not None:
            extra.append(f"θ = {_num(res['theta']):.2f} deg")
        if _num(res.get("fo")) is not None:
            extra.append(f"f₀ = {_num(res['fo']):.3f}")
        rows.append([method_label(name), "" if fs is None else f"{fs:.3f}",
                     ", ".join(extra)])
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
    table = _fs_table(solutions, counter)
    if table is not None:
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
        "runs, and reports what it finds. The findings that were active for this "
        "analysis are listed below, in the checker's own words. They are reported "
        "here so that a reviewer sees them, not only the engineer who ran the "
        "analysis."))

    findings = list(getattr(report, "findings", []) or [])
    if not findings:
        sec.blocks.append(Prose("The model checks raised no findings."))
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
