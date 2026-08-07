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

"""Checks for the Analysis Report.

What is being defended:

  A. THE CONTENT TREE — a solved model produces the sections the report is
     specified to have, in order, with the tables holding the model's own
     numbers and every figure rendered to a real file.

  B. THE TOGGLES ARE REAL — turning a section off removes it from the tree, and
     turning a parent off takes its children with it. A checkbox that changed
     nothing would pass a "the report has sections" test.

  C. THE METHOD PICKER — the picked method drives the critical-surface figure
     and the slice table, and NOT the factor-of-safety summary: EVERY method the
     solver offers is reported whichever one the detail follows.

  H. NOTHING IS SAID THAT IS NOT SO — the report describes the features the model
     actually has. A dry section gets one plain statement about water, not a
     key-value block explaining how its water loads are derived; a model with no
     piles gets no Piles section; the model checks are filtered to the analyses
     the report contains; and an empty title-page field prints no row.

  D. THE DOCUMENT — the .docx is a real OOXML package: the title is in it, the
     template's styles are referenced, the slice table sits in a landscape
     section, the table of contents is a TOC field whose cached result is the
     report's own heading list, and the running head and foot carry live
     fields. Its tables are fitted to their CONTENT — fixed columns, measured,
     ending where the content ends rather than ruled across the page, cut down
     only where the content would overrun the text width, and indented so their
     borders line up with the body text — and generating one writes one file.

  E. THE COLUMN REGISTRY — every column the registry marks report-worthy exists
     in a solved slice_df, so the registry cannot drift away from the solver.

  F. THE SHARED PLOT — ``mode="shared"`` suppresses the trial surfaces and draws
     the water line a reservoir boundary states, which is where the derived
     water loads on the same figure come from.

  J. EVERY FIGURE AND TABLE IS CITED — the convention a technical report is read
     under. Every numbered block is named by number in prose that prints
     whenever the block does, from its own section or one above it; nothing is
     cited that the options left out; and each citation is a live
     cross-reference to the bookmark on the block's caption.

  K. THE OTHER TWO ENGINES — a seepage run is documented with what it solved,
     the mesh and conductivities it solved on, and the flow that came out of it;
     a strength reduction run states its factor of safety in bold and a single
     trial states none. Neither section is built without its engine's solution.

  G. THE DIALOG — it constructs offscreen, its toggles move the options it hands
     over, and what it remembers survives a second construction.

One small LEM model is solved once and shared. A second model (a dam with a
reservoir head boundary) is loaded, not solved: the shared-plot check reads its
inputs only. The seepage and finite element models are loaded and their SHIPPED
solutions read back, so no check here spends a seepage iteration or an SSRM
bisection.
"""

import contextlib
import copy
import io
import itertools
import json
import math
import os
import re
import shutil
import sys
import tempfile
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REINF_XLSX = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_reinf.xlsx")
DAM_XLSX = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_dam.xlsx")
RFACE_XLSX = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_rface.xlsx")
SIMPLE1_XLSX = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_simple1.xlsx")

# Models carrying one feature each, for the checks that can only be made where
# the feature is present: a right-facing slope with real horizontal forces on its
# slices (xslope_rface carries none, so its mirror rule reads correct however it
# is written), a slope reinforced entirely by passive geosynthetic, one reinforced
# axially, and one with water in a tension crack.
MIRROR_XLSX = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                           "vp039b.xlsx")
PASSIVE_XLSX = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                            "vp088.xlsx")
AXIAL_XLSX = os.path.join(_REPO, "docs", "inputs", "slope",
                          "xslope_nail_axial.xlsx")
TENSION_XLSX = os.path.join(_REPO, "docs", "lem", "files",
                            "xslope_tension_KEY.xlsx")

#: The width a report figure prints at: the text column of a US Letter portrait
#: page inside the report's own margins. A figure is rendered wider than this and
#: placed at it, so type in a figure reaches the page in that ratio.
_TEXT_WIDTH_IN = 6.5

#: The smallest a printed number may be and still be read. Below this the digits
#: of a slice key close up on a laser printer and the key stops keying anything.
_LEGIBLE_PT = 5.0

_SOLVED = {}
_MODELS = {}


def load_slope_data_cached(xlsx):
    """One model, loaded once however many checks read it."""
    if xlsx not in _MODELS:
        from xslope.fileio import load_slope_data
        with contextlib.redirect_stdout(io.StringIO()):
            _MODELS[xlsx] = load_slope_data(xlsx)
    return _MODELS[xlsx]


def _solved():
    """The sample model with one circular surface solved by two methods, plus a
    small synthetic search bundle. Solved once and cached: nothing here reads a
    converged value, only the report's treatment of one."""
    if "it" in _SOLVED:
        return _SOLVED["it"]
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    slope_data = load_slope_data(REINF_XLSX)
    ok, out = generate_slices(slope_data, circle=slope_data["circles"][0],
                              num_slices=15)
    if not ok:
        raise RuntimeError(f"the sample model produced no slices: {out}")
    slice_df, surface = out

    bundles = []
    for name in ("spencer", "bishop"):
        df = slice_df.copy()
        results = solve_selected(name, df)
        bundles.append({"slice_df": df, "failure_surface": surface,
                        "results": results, "search": None, "method": name})
    # A search bundle built from the surfaces that were solved: the search
    # section reads counts and factors of safety, not the search's internals, so
    # this exercises every line of it without spending a real search.
    bundles[0]["search"] = {
        "kind": "circular",
        "fs_cache": [{"Xo": 10.0, "Yo": 50.0, "Depth": 0.0,
                      "FS": b["results"]["FS"], "slices": b["slice_df"],
                      "failure_surface": surface,
                      "solver_result": b["results"]} for b in bundles],
        "search_path": [{"x": 10.0, "y": 50.0, "FS": None},
                        {"x": 9.0, "y": 49.0, "FS": bundles[0]["results"]["FS"]}],
        "circle_cache": None,
    }
    _SOLVED["it"] = (slope_data, {"lem": bundles})
    return _SOLVED["it"]


def _build(options=None, figure_dir=None, fast=True):
    """A content tree from the sample model.

    ``fast`` builds it with the figures switched off: a tree check reads
    sections, tables and numbers, and rendering three 300-dpi plots for each of
    the twenty trees these checks build would be most of the run time for
    nothing. The checks that are ABOUT the figures pass ``fast=False`` and get
    the shipped defaults, resolution included.
    """
    from xslope.report import build_report
    slope_data, solutions = _solved()
    opts = {"input_path": REINF_XLSX, "title": "Sample Levee"}
    if fast:
        opts.update({"pd_figure": False, "lem_inputs_figure": False,
                     "lem_search_figure": False, "lem_solution_figure": False})
    opts.update(options or {})
    if figure_dir is None:
        figure_dir = tempfile.mkdtemp(prefix="xslope_report_")
    return build_report(slope_data, solutions, opts, figure_dir)


# --------------------------------------------------------------------------
# A. the content tree
# --------------------------------------------------------------------------

#: The sections a report of this solved LEM model has, in order, with the shipped
#: defaults. The sample carries reinforcement lines and no piles, so Reinforcement
#: is here and Piles is not — the two are separate sections, each present only
#: where the model has that feature. Model Checks is opt-in and is absent.
#:
#: Project Definition holds only what is true of the model whatever is run on it
#: — the geometry and the units. Everything an engine READS stands where that
#: engine is documented: a material's shear strength, the pore pressures a
#: stability analysis takes, the loads it applies, and the reinforcement it
#: resolves onto the slice bases are all limit equilibrium inputs and are all in
#: the limit equilibrium section.
#:
#: There is no factor of safety summary: the shipped default documents ONE
#: method, and a one-row comparison table restates the number that method's own
#: section states. Everything after Loads is that method's block, headed by its
#: name. The SEARCH is inside the block: it found the critical surface for that
#: method, and two methods searched separately settle on different surfaces.
EXPECTED_SECTIONS = [
    (1, "Traceability"),
    (1, "Project Definition"),
    (1, "Limit Equilibrium Analysis"),
    (2, "Analysis Inputs"),
    (2, "Materials"),
    (2, "Loads"),
    # The reinforcement stands with the engine that reads it, below the inputs
    # that engine reads and above the answers it produced.
    (2, "Reinforcement"),
    (2, "Spencer's Method"),
    (3, "Search for the Critical Surface"),
    (3, "Results"),
    (3, "Slice Table"),
    (3, "Calculations"),
]


def test_tree():
    """The tree has the specified sections, in order, and the units statement
    leads the Project Definition rather than trailing it."""
    fails = []
    report = _build()
    got = report.section_titles()
    if got != EXPECTED_SECTIONS:
        fails.append(f"section order is {got}, expected {EXPECTED_SECTIONS}")
    if not report.meta.get("title"):
        fails.append("the report carries no title")
    if not report.meta.get("date"):
        fails.append("the report carries no date")

    # The units statement is a lead block of Project Definition: a reader meets
    # the section's numbers already knowing what they are in.
    pd_sec = next((s for s in report.sections if s.title == "Project Definition"),
                  None)
    if pd_sec is None:
        return fails + ["there is no Project Definition section"]
    kinds = [b.kind for b in pd_sec.blocks]
    prose = [b.text for b in pd_sec.blocks if b.kind == "prose"]
    units = [i for i, t in enumerate(prose) if "units" in t or "unit system" in t]
    if not units:
        fails.append(f"no units statement in the Project Definition: {prose}")
    elif units[0] != 1:
        fails.append(f"the units statement is prose {units[0]} of the section; it "
                     f"belongs immediately after the intro")
    with_fig = _build({"pd_figure": True})
    fsec = next(s for s in with_fig.sections if s.title == "Project Definition")
    kinds = [b.kind for b in fsec.blocks]
    if "figure" not in kinds:
        fails.append("the model figure is not a block of Project Definition")
    else:
        texts = [b.text for b in fsec.blocks[:kinds.index("figure")]
                 if b.kind == "prose"]
        if not any("units" in t or "unit system" in t for t in texts):
            fails.append("the units statement comes after the model figure")
    return fails


def test_tables_carry_the_model():
    """The tables hold the model's own numbers, not placeholders."""
    fails = []
    slope_data, solutions = _solved()
    report = _build()
    tables = {t.caption: t for t in report.tables()}

    mats = tables.get("Material properties")
    if mats is None:
        return ["there is no materials table"]
    if len(mats.rows) != len(slope_data["materials"]):
        fails.append(f"the materials table has {len(mats.rows)} rows for "
                     f"{len(slope_data['materials'])} materials")
    names = {r[1] for r in mats.rows}
    for m in slope_data["materials"]:
        if m["name"] not in names:
            fails.append(f"material {m['name']!r} is missing from the table")
    # Only populated columns: this model states no gamma_sat and no r_u.
    heads = " ".join(mats.headers)
    for absent in ("γ_sat", "r_u"):
        if absent in heads:
            fails.append(f"the materials table prints {absent}, which no "
                         f"material populates")
    for present in ("γ", "c", "φ"):
        if present not in heads:
            fails.append(f"the materials table omits {present}")

    # The comparison table exists only where there is a comparison, so it is
    # asked for from a report that documents both methods.
    both = _build({"method": ["spencer", "bishop"]})
    fs = {t.caption: t for t in both.tables()}.get("Computed factors of safety")
    if fs is None:
        fails.append("a two-method report carries no factor-of-safety table")
    else:
        for b in solutions["lem"]:
            want = f"{b['results']['FS']:.3f}"
            if not any(want in row for row in fs.rows):
                fails.append(f"FS {want} ({b['method']}) is not in the summary table")

    reinf = tables.get("Reinforcement lines")
    if reinf is None:
        fails.append("there is no reinforcement table")
    elif len(reinf.rows) != len(slope_data["reinforcement_lines"]):
        fails.append(f"the reinforcement table has {len(reinf.rows)} rows for "
                     f"{len(slope_data['reinforcement_lines'])} lines")

    # No header word wider than its column: Word does not wrap a word that will
    # not fit, it breaks it, and "Applicatio n" is what that looks like in print.
    from xslope.report import max_header_word
    for t in report.tables():
        if t.landscape:
            continue                        # the slice table gets a whole page
        limit = max_header_word(len(t.headers))
        for head in t.headers:
            for word in str(head).split():
                if len(word) > limit:
                    fails.append(f"{t.caption!r}: the header word {word!r} is "
                                 f"{len(word)} characters, over the {limit} a "
                                 f"{len(t.headers)}-column table allows, and will "
                                 f"break mid-word")

    slices = [t for t in report.tables() if t.landscape]
    if len(slices) != 1:
        fails.append(f"{len(slices)} landscape tables; the slice table should be "
                     f"the only one")
    elif len(slices[0].rows) != len(solutions["lem"][0]["slice_df"]):
        fails.append(f"the slice table has {len(slices[0].rows)} rows for "
                     f"{len(solutions['lem'][0]['slice_df'])} slices")
    elif not slices[0].legend:
        fails.append("the slice table carries no column legend")
    return fails


def test_figures_rendered():
    """Every figure block points at a real, non-trivial PNG."""
    fails = []
    from xslope.report import DEFAULT_OPTIONS, FIGURE_DPI

    if DEFAULT_OPTIONS["dpi"] != FIGURE_DPI or FIGURE_DPI < 300:
        fails.append(f"figures default to {DEFAULT_OPTIONS['dpi']} dpi; a report "
                     f"figure is a print figure")
    tmp = tempfile.mkdtemp(prefix="xslope_report_fig_")
    report = _build(figure_dir=tmp, fast=False)
    figures = report.figures()
    if len(figures) < 3:
        fails.append(f"only {len(figures)} figures were built (model, search, "
                     f"solution expected)")
    numbers = [f.number for f in figures]
    if numbers != sorted(set(numbers)) or (numbers and numbers[0] != 1):
        fails.append(f"figures are not numbered 1..n in order: {numbers}")
    for f in figures:
        if not os.path.exists(f.path):
            fails.append(f"figure {f.number} was not written to {f.path}")
        elif os.path.getsize(f.path) < 20000:
            fails.append(f"figure {f.number} is {os.path.getsize(f.path)} bytes — "
                         f"too small to be a rendered plot")
        if not f.caption:
            fails.append(f"figure {f.number} has no caption")
    return fails


def test_traceability():
    """The stamp identifies the inputs the answers came from."""
    from xslope.report import file_digest
    from xslope._version import __version__

    report = _build()
    kv = [b for b in report.blocks("keyvalues")]
    if not kv:
        return ["the traceability block is missing"]
    items = dict(kv[0].items)
    fails = []
    if items.get("xslope version") != __version__:
        fails.append(f"the stamp reports version {items.get('xslope version')!r}")
    if items.get("Input file") != os.path.basename(REINF_XLSX):
        fails.append(f"the stamp names {items.get('Input file')!r} as the input")
    digest = file_digest(REINF_XLSX)
    if items.get("Input file SHA-256") != digest:
        fails.append("the stamp's digest is not the input file's SHA-256")
    if len(digest) != 64:
        fails.append(f"the digest is {len(digest)} characters, not a SHA-256")
    return fails


# --------------------------------------------------------------------------
# B. the toggles
# --------------------------------------------------------------------------

def test_toggles():
    """Every content toggle removes what it names, and a parent takes its
    children with it."""
    fails = []
    full = [t for _lvl, t in _build().section_titles()]

    cases = [
        ({"traceability": False}, "Traceability"),
        ({"lem_search": False}, "Search for the Critical Surface"),
        ({"lem_slice_table": False}, "Slice Table"),
        ({"lem_calculations": False}, "Calculations"),
        ({"lem_materials": False}, "Materials"),
        ({"lem_loads": False}, "Loads"),
        ({"lem_members": False}, "Reinforcement"),
    ]
    for opts, gone in cases:
        titles = [t for _lvl, t in _build(opts).section_titles()]
        if gone in titles:
            fails.append(f"{opts} left {gone!r} in the report")
        if len(titles) != len(full) - 1:
            fails.append(f"{opts} changed the section count by "
                         f"{len(full) - len(titles)}, not 1")

    # Model checks are opt-in: the default report does not carry the section, and
    # asking for it adds one.
    if "Model Checks" in full:
        fails.append("Model Checks is in the default report; it is opt-in")
    on = [t for _lvl, t in _build({"model_checks": True}).section_titles()]
    if "Model Checks" not in on:
        fails.append("model_checks=True produced no Model Checks section")
    if len(on) != len(full) + 1:
        fails.append(f"model_checks=True changed the section count by "
                     f"{len(on) - len(full)}, not 1")

    # The units statement is a block, not a section: its toggle takes the
    # paragraph and leaves the section it leads.
    def units_prose(report):
        return [b.text for b in report.blocks("prose")
                if "units" in b.text or "unit system" in b.text]

    if not units_prose(_build()):
        fails.append("the default report states no units")
    off = _build({"pd_units": False})
    if units_prose(off):
        fails.append("pd_units=False still stated the units")
    if "Project Definition" not in [t for _l, t in off.section_titles()]:
        fails.append("pd_units=False removed the Project Definition section")

    # A parent off takes the whole branch.
    titles = [t for _lvl, t in _build({"project_definition": False}).section_titles()]
    if "Project Definition" in titles:
        fails.append("project_definition=False left 'Project Definition' in the "
                     "report")
    # And takes nothing that is not its own. The reinforcement is read by the
    # stability analysis and stands with it: a report built without the general
    # description of the model still documents what that analysis was run on.
    if "Reinforcement" not in titles:
        fails.append("project_definition=False took the Reinforcement section "
                     "out of the limit equilibrium analysis with it")
    titles = [t for _lvl, t in _build({"lem": False}).section_titles()]
    for gone in ("Limit Equilibrium Analysis", "Slice Table", "Results",
                 "Calculations", "Materials", "Loads", "Reinforcement"):
        if gone in titles:
            fails.append(f"lem=False left {gone!r} in the report")

    fails += _member_switch_stands_under_its_engine()
    return fails


def _member_switch_stands_under_its_engine():
    """The member tables are switched on and off per engine, and each switch is
    a child of the engine that prints them.

    The dialog forces a sub-item false whenever its parent section is off, so
    where a switch sits IS what it can mean. Under Project Definition — where it
    was — a report built without the general description of the model lost the
    properties its stability analysis was run on. Under one engine and shared
    with the other it would be the same defect turned sideways: the finite
    element tables would go when the limit equilibrium section did.

    So there is one switch per engine, each under its own engine, and neither
    reaches across.
    """
    fails = []
    from studio.report_dialog import CONTENT_TREE
    from xslope.report import DEFAULT_OPTIONS

    parents = {}
    for key, _label, _tip, children in CONTENT_TREE:
        for child_key, _l, _t in children:
            parents[child_key] = key
    for key, engine in (("lem_members", "lem"), ("fem_members", "fem")):
        if key not in DEFAULT_OPTIONS:
            fails.append(f"{key} is not an option the builder reads")
            continue
        if DEFAULT_OPTIONS[key] is not True:
            fails.append(f"{key} is off by default; the properties an analysis "
                         f"was run on are part of documenting it")
        if key not in parents:
            fails.append(f"the dialog offers no {key} row")
        elif parents[key] != engine:
            fails.append(f"the {key} row is a child of {parents[key]!r}; it "
                         f"belongs under {engine!r}, because the dialog turns a "
                         f"child off with its parent and this switch may only "
                         f"be turned off by the engine that prints it")
    # The retired key is gone from both the tree and the builder: a row nothing
    # reads is a box that does nothing.
    if "pd_reinforcement" in parents:
        fails.append("the dialog still offers pd_reinforcement, which the "
                     "builder no longer reads")
    if "pd_reinforcement" in DEFAULT_OPTIONS:
        fails.append("pd_reinforcement is still a builder option")

    # And the switches really are independent: one engine's tables survive the
    # other engine's section being switched off entirely.
    slope_data, bundle = _fem_1d_bundle(FEM_REINF_XLSX)
    with tempfile.TemporaryDirectory() as tmp:
        from xslope.report import build_report
        report = build_report(slope_data, {"fem": [bundle]},
                              dict(FAST_FIGURES, pd_figure=False,
                                   lem=False, lem_members=False), tmp)
    if "Reinforcement" not in [t for _l, t in report.section_titles()]:
        fails.append("switching the limit equilibrium members off took the "
                     "finite element section's reinforcement with them")
    with tempfile.TemporaryDirectory() as tmp:
        off = build_report(slope_data, {"fem": [bundle]},
                           dict(FAST_FIGURES, pd_figure=False,
                                fem_members=False), tmp)
    if "Reinforcement" in [t for _l, t in off.section_titles()]:
        fails.append("fem_members=False left the reinforcement section in the "
                     "finite element report")
    return fails

    # A figure toggle removes the figure, not the section that holds it.
    on = _build({"pd_figure": True})
    off = _build({"pd_figure": False})
    if not any(f.source == "shared model" for f in on.figures()):
        fails.append("pd_figure=True drew no model figure, so its 'off' case "
                     "proves nothing")
    if any(f.source == "shared model" for f in off.figures()):
        fails.append("pd_figure=False still drew the model figure")
    if "Project Definition" not in [t for _lvl, t in off.section_titles()]:
        fails.append("pd_figure=False removed the Project Definition section")
    return fails


# --------------------------------------------------------------------------
# C. the method picker
# --------------------------------------------------------------------------

def test_method_picker():
    """The picked method drives the detail and leaves the summary alone."""
    fails = []
    from xslope.report import method_label, solved_methods

    slope_data, solutions = _solved()
    if solved_methods(solutions) != ["spencer", "bishop"]:
        fails.append(f"the sample offers {solved_methods(solutions)}")

    reports = {}
    for name in ("spencer", "bishop"):
        reports[name] = _build({"method": name,
                                "lem_search_figure": False,
                                "lem_solution_figure": True})

    a, b = reports["spencer"], reports["bishop"]

    # The detail follows the pick: figure provenance and slice-table caption.
    for name, report in reports.items():
        sources = [f.source for f in report.figures()]
        if f"{name} solution surface" not in sources:
            fails.append(f"method={name}: no solution-surface figure for it "
                         f"(sources {sources})")
        captions = [t.caption for t in report.tables() if t.landscape]
        if not captions or method_label(name) not in captions[0]:
            fails.append(f"method={name}: the slice table is captioned "
                         f"{captions}")
    if ([f.source for f in a.figures()] == [f.source for f in b.figures()]):
        fails.append("the two methods produced the same figure provenance — the "
                     "picker does not reach the figures")
    if ([f.path for f in a.figures()] == [f.path for f in b.figures()]):
        fails.append("the two methods wrote the same figure files")

    # The summary follows the pick, because the summary IS the pick: a report
    # of one method carries none, and a report of two lists those two.
    for report, name in ((a, "spencer"), (b, "bishop")):
        if _fs_rows(report) is not None:
            fails.append(f"method={name}: a single-method report was given a "
                         f"comparison table")

    # An unknown pick falls back rather than failing.
    report = _build({"method": "no_such_method"})
    if "Limit Equilibrium Analysis" not in [t for _l, t in report.section_titles()]:
        fails.append("an unrecognised method dropped the LEM section")
    return fails


def _fs_rows(report):
    """The factor-of-safety summary as ``(headers, rows)``, or None."""
    for t in report.tables():
        if t.caption == "Computed factors of safety":
            return (t.headers, t.rows)
    return None


def _fs_table_block(report):
    for t in report.tables():
        if t.caption == "Computed factors of safety":
            return t
    return None


def _fs_prose(report):
    """The paragraphs of the Factors of Safety section, as one string."""
    for section in report.sections:
        for _l, sub in section.walk():
            if sub.title == "Factors of Safety":
                return " ".join(b.text for b in sub.blocks if b.kind == "prose")
    return ""


#: The two things a summary paragraph can say about where its numbers came from,
#: and they are alternatives. A report that says both says that every method
#: found its own surface AND that no surface was searched for.
_SEARCHED_CLAIM = "searched for"
_UNSEARCHED_CLAIM = "No search was performed"


def _searched_bundles(bundles, which):
    """The same bundles with a search attached to the named methods and taken off
    the rest — the three provenances, off one solved model."""
    out = []
    for bundle in bundles:
        copy = dict(bundle)
        if copy.get("method") in which:
            copy["search"] = copy.get("search") or {
                "kind": "circular", "fs_cache": [], "search_path": []}
        else:
            copy["search"] = None
        out.append(copy)
    return out


def test_the_fs_summary_says_where_its_numbers_came_from():
    """The summary paragraph states one provenance for its rows, and it is the
    true one.

    Three cases, and exactly one sentence is true of any table. Every method
    searched: each row is that method's own critical factor of safety, on its own
    critical surface. Some searched: those rows are minima and the rest are the
    specified surface. None searched: every row is the specified surface, and no
    row is a minimum over anything.

    The paragraph used to open by saying every method finds its own surface,
    whatever the bundles held, and then append the no-search sentence under it —
    so a report solved entirely on one specified circle told the reader both at
    once. The two claims are alternatives, and no paragraph may carry both.
    """
    fails = []
    from xslope.report import build_report, method_label

    slope_data, solutions = _solved()
    bundles = solutions.get("lem") or []
    methods = [b.get("method") for b in bundles]
    if len(bundles) < 2:
        return ["the sample solves fewer than two methods; the summary is untested"]

    def built(which):
        with tempfile.TemporaryDirectory() as tmp:
            return build_report(
                slope_data, {"lem": _searched_bundles(bundles, which)},
                {"pd_figure": False, "lem_search_figure": False,
                 "lem_solution_figure": False, "lem_slice_key": False,
                 "lem_slice_table": False, "lem_calculations": False,
                 "method": methods}, tmp)

    cases = [("every method searched", set(methods), True, False),
             ("one method searched", {methods[0]}, True, True),
             ("no method searched", set(), False, True)]
    for where, which, wants_searched, wants_specified in cases:
        said = _fs_prose(built(which))
        if not said:
            fails.append(f"{where}: the report carries no summary paragraph")
            continue
        # Never both. This is the defect: the two claims contradict each other.
        if _SEARCHED_CLAIM in said and _UNSEARCHED_CLAIM in said:
            fails.append(f"{where}: the paragraph says a search was made and "
                         f"that none was: {said!r}")
        if ("finds its own surface" in said
                and _UNSEARCHED_CLAIM in said):
            fails.append(f"{where}: the paragraph says every method finds its "
                         f"own surface and that no search was performed: "
                         f"{said!r}")
        if wants_searched and _SEARCHED_CLAIM not in said:
            fails.append(f"{where}: the paragraph does not say a search was "
                         f"made: {said!r}")
        if not wants_searched and _SEARCHED_CLAIM in said:
            fails.append(f"{where}: no bundle carries a search and the "
                         f"paragraph says one was made: {said!r}")
        if wants_specified and "specified in the input" not in said:
            fails.append(f"{where}: rows stand on the specified surface and the "
                         f"paragraph does not say so: {said!r}")
        if not wants_specified and "specified in the input" in said:
            fails.append(f"{where}: every method searched and the paragraph "
                         f"still sends the reader to the specified surface: "
                         f"{said!r}")

    # Where every method searched, the paragraph says what the rows ARE: each
    # method's own critical value. That is the summary the owner asked for, and
    # the wording is what makes the column readable as one.
    said = _fs_prose(built(set(methods)))
    for phrase in ("critical factor of safety", "its own critical surface",
                   "minimum over the family"):
        if phrase not in said:
            fails.append(f"every method searched: the paragraph does not say "
                         f"{phrase!r}: {said!r}")

    # And where only some did, the ones that did are named.
    said = _fs_prose(built({methods[0]}))
    if method_label(methods[0]) not in said:
        fails.append(f"one method searched: the paragraph does not name "
                     f"{method_label(methods[0])}: {said!r}")
    return fails


def test_multi_method_detail():
    """Several methods, several full detail blocks — in the order asked for, each
    headed and captioned with its own method's name.

    The failure this defends against is the one that makes a multi-method report
    worse than a single-method one: two blocks of numbers a reader cannot tell
    apart. Every heading, figure caption and table caption in a block has to name
    its method, and the blocks have to come in a stated order.

    A block is headed by the method's bare name, so the blocks are found by
    matching the level-2 headings against the labels of every method the solver
    offers — which is what makes a block for a method nobody asked for a failure
    as well.
    """
    fails = []
    from xslope.report import detail_bundle, method_label, supported_methods

    slope_data, solutions = _solved()
    labels = {method_label(m) for m in supported_methods()}

    wanted = ["bishop", "spencer"]
    report = _build({"method": wanted, "lem_solution_figure": True})
    titles = [t for lvl, t in report.section_titles() if lvl == 2]

    heads = [t for t in titles if t in labels]
    expected = [method_label(m) for m in wanted]
    if heads != expected:
        fails.append(f"the detail blocks are {heads}, expected {expected}")

    # One of everything per method, each naming its own.
    for m in wanted:
        label = method_label(m)
        sources = [f.source for f in report.figures()]
        if f"{m} solution surface" not in sources:
            fails.append(f"{m}: no solution-surface figure ({sources})")
        captions = [f.caption for f in report.figures()
                    if f.source == f"{m} solution surface"]
        if not captions or label not in captions[0]:
            fails.append(f"{m}: the solution-surface caption is {captions}")
        slices = [t for t in report.tables() if t.landscape and label in t.caption]
        if len(slices) != 1:
            fails.append(f"{m}: {len(slices)} slice tables captioned for it")

    # Each block carries the whole subtree, not just a heading — and where the
    # method searched, its own search stands at the head of that subtree.
    blocks = {}
    for sec in report.sections:
        for _lvl, node in sec.walk():
            if node.title in labels:
                blocks[node.title] = node
    for m in wanted:
        node = blocks.get(method_label(m))
        if node is None:
            continue
        detail, _note = detail_bundle(slope_data, solutions, m)
        want_kids = (["Search for the Critical Surface"]
                     if (detail or {}).get("search") else [])
        want_kids += ["Results", "Slice Table", "Calculations"]
        kids = [c.title for c in node.children]
        if kids != want_kids:
            fails.append(f"{node.title!r} carries {kids}, expected {want_kids}")

    # Figures and tables are numbered straight through, whatever the block.
    numbers = [f.number for f in report.figures()]
    if numbers != list(range(1, len(numbers) + 1)):
        fails.append(f"the figures are not numbered 1..n across methods: {numbers}")
    numbers = [t.number for t in report.tables()]
    if numbers != list(range(1, len(numbers) + 1)):
        fails.append(f"the tables are not numbered 1..n across methods: {numbers}")

    # The summary comes ONCE, and before the first detail block: it lists every
    # method, so it belongs to the analysis. The search does NOT — it belongs to
    # the block whose surface it found, so there is one per method that searched
    # and none above the blocks.
    if sum(1 for t in report.tables()
           if t.caption == "Computed factors of safety") != 1:
        fails.append("the factor-of-safety summary is not printed exactly once")
    searched = [m for m in wanted
                if (detail_bundle(slope_data, solutions, m)[0] or {}).get("search")]
    levels = [lvl for lvl, t in report.section_titles()
              if t == "Search for the Critical Surface"]
    if len(levels) != len(searched):
        fails.append(f"the report documents {len(levels)} searches for the "
                     f"{len(searched)} method(s) that searched: {searched}")
    if any(lvl != 3 for lvl in levels):
        fails.append(f"a search is written at heading level {levels}; it is a "
                     f"child of the method block whose surface it found, not a "
                     f"section standing over all of them")
    if "Factors of Safety" not in titles:
        fails.append(f"there is no factor-of-safety section: {titles}")
    elif heads and titles.index("Factors of Safety") > titles.index(heads[0]):
        fails.append("the summary comes after the detail blocks")

    # Order is the caller's, not the solver's.
    other = _build({"method": ["spencer", "bishop"]})
    heads2 = [t for lvl, t in other.section_titles() if lvl == 2 and t in labels]
    if heads2 != list(reversed(expected)):
        fails.append(f"reversing the request gave {heads2}")

    # A bare string is still one method — every caller written before the list
    # option existed keeps working.
    one = _build({"method": "bishop"})
    heads3 = [t for lvl, t in one.section_titles() if lvl == 2 and t in labels]
    if heads3 != [method_label("bishop")]:
        fails.append(f"method='bishop' produced {heads3}")

    # A method that was never RUN is reported on the critical surface, and says
    # so — the same thing the summary already does for it.
    from xslope.report import solved_methods
    unrun = next(m for m in ("janbu", "lowe", "corps")
                 if m not in solved_methods(solutions))
    extra = _build({"method": [unrun]})
    prose = " ".join(b.text for b in extra.blocks("prose"))
    if method_label(unrun) not in [t for _l, t in extra.section_titles()]:
        fails.append(f"a method that was not run ({unrun}) got no detail block")
    if "It was not run in the analysis" not in prose:
        fails.append(f"the {unrun} block does not say the report solved it")
    if not [t for t in extra.tables() if t.landscape]:
        fails.append(f"the {unrun} block carries no slice table")
    return fails


def test_each_engine_presents_the_loads_it_applies():
    """The loads are an analysis input, and the second engine cites the first
    rather than printing the table twice.

    A distributed load is not a property of the section; it is something an
    analysis puts on it. So it stands with the engine that applies it, not in a
    general description of the model. Where both engines are documented the
    blocks are identical — same points, same pressures — so the finite element
    section points at the limit equilibrium table: two copies of one table are
    two things to keep in step.
    """
    fails = []
    from xslope.report import build_report

    slope_data, solutions = _load_bearing_model()
    if slope_data is None:
        return ["no loaded model with both engines; the citation is untested"]

    with tempfile.TemporaryDirectory() as tmp:
        report = build_report(slope_data, solutions,
                              dict(FAST_FIGURES, pd_figure=False), tmp)

    # Not in Project Definition.
    pd_sec = next((s for s in report.sections
                   if s.title == "Project Definition"), None)
    if pd_sec is not None:
        under = [t for _l, t in pd_sec.walk()]
        if "Loads" in under:
            fails.append(f"Project Definition still carries a Loads section: "
                         f"{under}")

    # One table, in the limit equilibrium section; the other engine cites it.
    loads = [t for t in report.tables() if t.caption == "Distributed loads"]
    if len(loads) != 1:
        fails.append(f"{len(loads)} distributed-load tables were printed; the "
                     f"second engine must cite the first")
    for head in ("Limit Equilibrium Analysis",
                 "Deformation and Strength Reduction"):
        sec = next((s for s in report.sections if s.title == head), None)
        if sec is None:
            fails.append(f"the report has no {head} section")
            continue
        sub = next((s for _l, s in sec.walk() if s.title == "Loads"), None)
        if sub is None:
            fails.append(f"{head} presents no loads")
            continue
        if head.startswith("Deformation"):
            if any(b.kind == "table" for b in sub.blocks):
                fails.append("the finite element section printed its own copy "
                             "of the loads table")
            said = " ".join(b.text for b in sub.blocks if b.kind == "prose")
            if loads and f"Table {loads[0].number}" not in said:
                fails.append(f"the finite element section does not cite the "
                             f"loads table: {said!r}")

    fails += _water_load_mechanism_checks()
    return fails


def _loads_prose(xlsx):
    """The Loads subsection's prose for a model, as one string."""
    from xslope.report import _loads_section, water_features, _Counter
    slope_data = load_slope_data_cached(xlsx)
    section = _loads_section(slope_data, water_features(slope_data), _Counter())
    return " ".join(b.text for b in section.blocks if b.kind == "prose")


#: What the Loads prose has to say about a load the user never entered: how big
#: it is, and which way it points. Anything less leaves an engineer reading a
#: factor of safety that stands partly on a force the report never described.
_WATER_MECHANISM = ("depth of water", "unit weight of water",
                    "normal to the ground surface")


def _water_load_mechanism_checks():
    """A load the engine derived is a load the report states the derivation of.

    Where the water surface stands above the ground surface, the engine applies
    the standing water to the ground: pressure equal to the depth of water above
    each point times the unit weight of water, acting normal to the surface and
    tapering to zero at the shoreline. That force enters the factor of safety and
    the user never typed it, so the analysis is not documented until the report
    says what it is.

    Said only where there is one. A section with no water standing on it carries
    no such load, and describing one would describe an analysis that did not
    happen.
    """
    fails = []
    from xslope.water import with_water_loads, has_derived_loads
    import xslope.report as report_mod

    ponding, dry = DAM_XLSX, SIMPLE1_XLSX
    if not has_derived_loads(with_water_loads(load_slope_data_cached(ponding))):
        fails.append(f"{os.path.basename(ponding)} derives no water load, so "
                     f"the sentence describing one is untested")
    if has_derived_loads(with_water_loads(load_slope_data_cached(dry))):
        fails.append(f"{os.path.basename(dry)} derives a water load, so it "
                     f"cannot show the sentence staying away")

    said = _loads_prose(ponding)
    for named in _WATER_MECHANISM:
        if named not in said:
            fails.append(f"the Loads prose of a ponded section does not say "
                         f"{named!r}: {said!r}")

    # The unit weight is named, not printed again. The units paragraph declares
    # its value; a second copy is a number that can disagree with itself.
    gamma = load_slope_data_cached(ponding).get("gamma_water")
    if gamma and f"{float(gamma):g}" in said:
        fails.append(f"the Loads prose restates the unit weight of water "
                     f"({float(gamma):g}), which the units paragraph declares")

    quiet = _loads_prose(dry)
    for named in _WATER_MECHANISM:
        if named in quiet:
            fails.append(f"a section with no water standing on it is told how "
                         f"standing water is loaded ({named!r}): {quiet!r}")

    # Mutation: the check has to be able to fail. With the mechanism withheld
    # the ponded section's Loads prose goes back to naming a derived load
    # without saying what it is, which is the report the defect was found in.
    real = report_mod._water_load_mechanism
    report_mod._water_load_mechanism = lambda slope_data: None
    try:
        muted = _loads_prose(ponding)
    finally:
        report_mod._water_load_mechanism = real
    if any(named in muted for named in _WATER_MECHANISM):
        fails.append(f"the mechanism sentence is written somewhere other than "
                     f"the one place that decides whether there is one: "
                     f"{muted!r}")
    return fails


def _load_bearing_model():
    """A loaded model solved by BOTH engines, or ``(None, None)``.

    The finite element sample carries distributed loads and a shipped strength
    reduction solution; one limit equilibrium trial on its own starting circle
    puts the other engine beside it, which is the only shape in which the
    citation exists to be checked.
    """
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    slope_data, fem = _fem_bundle()
    if not slope_data.get("dloads") or not slope_data.get("circles"):
        return None, None
    with contextlib.redirect_stdout(io.StringIO()):
        ok, out = generate_slices(slope_data, circle=slope_data["circles"][0])
        if not ok:
            return None, None
        df, surface = out
        results = solve_selected("spencer", df)
    lem = [{"slice_df": df, "failure_surface": surface, "results": results,
            "search": None, "method": "spencer"}]
    return slope_data, {"lem": lem, "fem": fem}


def test_critical_is_a_word_a_search_earns():
    """A surface is called critical only where a search found it.

    "Critical" means least of a family. A factor of safety computed on a surface
    the user entered is the factor of safety of THAT surface, and printing it
    against the word critical invites an engineer to read a specified-surface
    answer as a searched minimum — the one misreading a stability report must not
    permit. Where no search ran the report says so outright rather than softening
    the adjective.

    Both shapes are built from the same sample: its Spencer bundle carries a
    search and its Bishop bundle does not.
    """
    fails = []
    _slope_data, solutions = _solved()

    searched = {b["method"]: bool(b.get("search")) for b in solutions["lem"]}
    if set(searched.values()) != {True, False}:
        return [f"the sample's bundles are all the same shape ({searched}); "
                f"one of the two wordings is untested"]

    report = _build({"method": ["spencer", "bishop"]})
    from xslope.report import method_label
    for name, ran in searched.items():
        block = next((s for _l, s in _sections(report)
                      if s.title == method_label(name)), None)
        if block is None:
            fails.append(f"{name}: no method block")
            continue
        said = " ".join(b.text for _l, node in block.walk()
                        for b in node.blocks if b.kind == "prose")
        captions = " ".join(f.caption for _l, node in block.walk()
                            for f in node.blocks if f.kind == "figure")
        if ran:
            if "critical surface" not in said:
                fails.append(f"{name} searched and its block never says critical "
                             f"surface: {said!r}")
            if "no search" in said.lower():
                fails.append(f"{name} searched and its block says no search ran")
        else:
            if "critical surface" in captions:
                fails.append(f"{name} did not search and its figure is captioned "
                             f"{captions!r}")
            if "specified surface" not in said:
                fails.append(f"{name} did not search and its block never calls "
                             f"the surface specified: {said!r}")
            if "no search" not in said.lower():
                fails.append(f"{name} did not search and the block never says so "
                             f"outright: {said!r}")
            # …and it says what that means for the number.
            if "not a minimum" not in said:
                fails.append(f"{name}: nothing says the reported factor is not a "
                             f"minimum over a family of surfaces: {said!r}")

    # The engine's own inputs say which it was, too.
    inputs = next((s for _l, s in _sections(report)
                   if s.title == "Analysis Inputs"), None)
    rows = dict(item for b in (inputs.blocks if inputs else [])
                if b.kind == "keyvalues" for item in b.items)
    if "Surface" not in rows:
        fails.append(f"the analysis inputs do not say how the surface was "
                     f"arrived at: {list(rows)}")

    # A report with no search anywhere never uses the word.
    nosearch = {"lem": [dict(b, search=None) for b in solutions["lem"]]}
    plain = _build({"method": ["spencer", "bishop"]})
    with tempfile.TemporaryDirectory() as tmp:
        from xslope.report import build_report
        plain = build_report(_slope_data, nosearch,
                             {"pd_figure": False, "lem_search_figure": False,
                              "lem_solution_figure": False,
                              "method": ["spencer", "bishop"]}, tmp)
    for block in plain.blocks("prose"):
        if "critical surface" in block.text and "No search" not in block.text \
                and "no search" not in block.text:
            fails.append(f"a report with no search calls a surface critical: "
                         f"{block.text!r}")
    if "Search for the Critical Surface" in [t for _l, t in plain.section_titles()]:
        fails.append("a report with no search carries a search section")
    return fails


def test_lem_inputs_figure():
    """The limit equilibrium section opens on the model the method of slices
    reads, as the other two engines' sections do.

    Every engine section states its own inputs, and the first of them is the
    section drawn the way that engine sees it. The shared figure of the Project
    Definition is not that figure: it deliberately suppresses the surface family,
    which is exactly what a limit equilibrium analysis is run over.
    """
    fails = []
    from xslope.report import DEFAULT_OPTIONS, build_report, planned_figures, \
        resolve_options

    if DEFAULT_OPTIONS.get("lem_inputs_figure") is not True:
        fails.append("the limit equilibrium model figure is not on by default")

    slope_data, solutions = _solved()

    def built(**over):
        opts = {"input_path": REINF_XLSX, "method": ["spencer", "bishop"],
                "pd_figure": False, "lem_search_figure": False,
                "lem_solution_figure": False}
        opts.update(FAST_FIGURES)
        opts.update(over)
        tmp = tempfile.mkdtemp(prefix="xslope_leminputs_")
        with contextlib.redirect_stdout(io.StringIO()):
            return opts, build_report(slope_data, solutions, opts, tmp)

    opts, report = built()
    figures = [f for f in report.figures() if f.source == "lem model"]
    if not figures:
        fails.append(f"the limit equilibrium section drew no model figure: "
                     f"{[f.source for f in report.figures()]}")
        return fails

    # One per section. Two methods are documented, and both are run on this one
    # model; a copy under each method's heading is the same picture twice.
    if len(figures) != 1:
        fails.append(f"the report drew {len(figures)} limit equilibrium model "
                     f"figures; it is an input of the section, not of a method")
    figure = figures[0]
    if "limit equilibrium" not in figure.caption.lower():
        fails.append(f"the figure is captioned {figure.caption!r}")

    # It stands in the section's own Analysis Inputs, ahead of the rows.
    lem = next((s for s in report.sections
                if s.title == "Limit Equilibrium Analysis"), None)
    inputs = next((c for c in (lem.children if lem else [])
                   if c.title == "Analysis Inputs"), None)
    if inputs is None:
        fails.append("the limit equilibrium section has no Analysis Inputs")
    else:
        kinds = [b.source if b.kind == "figure" else b.kind
                 for b in inputs.blocks if b.kind in ("figure", "keyvalues")]
        if kinds != ["prose", "lem model", "keyvalues"][1:]:
            if kinds[:2] != ["lem model", "keyvalues"]:
                fails.append(f"the limit equilibrium inputs are ordered {kinds}, "
                             f"not the model figure and then the rows")
        cites = [b for b in inputs.blocks
                 if b.kind == "prose" and f"Figure {figure.number}" in b.text]
        if not cites:
            fails.append(f"nothing cites Figure {figure.number}: "
                         f"{[b.text for b in inputs.blocks if b.kind == 'prose']}")
        else:
            said = cites[0].text
            if "method of slices" not in said:
                fails.append(f"the sentence does not say whose view of the model "
                             f"the figure is: {said!r}")
            # What it actually draws, and only that: the sample carries a
            # searched circular family and reinforcement.
            for named in ("the circle the search starts from",
                          "the reinforcement it carries"):
                if named not in said:
                    fails.append(f"the sentence does not name {named!r}: {said!r}")

    # It is the LEM view, not the shared one: the surface family the shared
    # figure suppresses is drawn here.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.figure as mplfig
    from xslope.plot import plot_inputs
    gids = {}
    for mode in ("lem", "shared"):
        f = mplfig.Figure()
        plot_inputs(slope_data, fig=f, mode=mode, show_title=False,
                    frame="content", show_mesh=False)
        gids[mode] = {a.get_gid() for a in f.axes[0].get_children()
                      if a.get_gid()}
    if not (gids["lem"] - gids["shared"]):
        fails.append(f"the limit equilibrium view draws nothing the shared view "
                     f"does not, so it is the same figure twice: {gids['lem']}")

    # …and it is the view the SECTION asks for. Reading the two modes apart says
    # nothing about which one the report drew, and the report ships a PNG the
    # check cannot interrogate, so the request itself is watched.
    import xslope.plot as xp
    real = xp.plot_inputs
    asked = []

    def spy(sd, **kw):
        asked.append(kw.get("mode"))
        return real(sd, **kw)

    xp.plot_inputs = spy
    try:
        built()
    finally:
        xp.plot_inputs = real
    if "lem" not in asked:
        fails.append(f"the limit equilibrium model figure was drawn in mode(s) "
                     f"{asked}, none of them the engine's own view")

    # The sentence claims only what this view draws. A pool a head boundary
    # states reaches the limit equilibrium figure as the derived load on the
    # ground surface, not as a water line, and a sentence promising a water
    # surface sends the reader looking for one that is not there.
    dam, dam_solutions = _restored(SEEP_XLSX)
    from xslope.report import water_features
    dam_feats = water_features(dam)
    if dam_feats["piezo"]:
        fails.append("the seepage sample now carries a piezometric line; the "
                     "no-water-line wording is untested")
    else:
        f = mplfig.Figure()
        plot_inputs(dam, fig=f, mode="lem", show_title=False, frame="content",
                    show_mesh=False)
        drawn = {a.get_gid() for a in f.axes[0].get_children() if a.get_gid()}
        if any("PIEZO" in (g or "") for g in drawn):
            fails.append(f"the seepage sample does draw a water line after all: "
                         f"{sorted(drawn)}")
        else:
            report2 = _cite_report(SEEP_XLSX, ("spencer",), engines=("seep",))
            lem2 = next(s for s in report2.sections
                        if s.title == "Limit Equilibrium Analysis")
            inputs2 = next(c for c in lem2.children
                           if c.title == "Analysis Inputs")
            said2 = " ".join(b.text for b in inputs2.blocks if b.kind == "prose")
            for claimed in ("water surface", "piezometric"):
                if claimed in said2.lower():
                    fails.append(f"the limit equilibrium figure's sentence "
                                 f"promises a {claimed} the view does not draw: "
                                 f"{said2!r}")
            if "the loads on it" not in said2:
                fails.append(f"the derived water load the figure does draw is "
                             f"not named: {said2!r}")

    # The toggle takes the figure and leaves the section.
    off_opts, off = built(lem_inputs_figure=False)
    if any(f.source == "lem model" for f in off.figures()):
        fails.append("lem_inputs_figure=False still drew the model figure")
    titles = [t for _l, t in off.section_titles()]
    if "Limit Equilibrium Analysis" not in titles or \
            "Analysis Inputs" not in titles:
        fails.append("lem_inputs_figure=False removed the section, not only the "
                     "figure")
    if len(off.figures()) != len(report.figures()) - 1:
        fails.append(f"lem_inputs_figure=False changed the figure count by "
                     f"{len(report.figures()) - len(off.figures())}, not 1")

    # And the count a caller is promised matches what was built, both ways.
    for these, made in ((opts, report), (off_opts, off)):
        planned = planned_figures(slope_data, solutions, resolve_options(these))
        if planned != len(made.figures()):
            fails.append(f"lem_inputs_figure="
                         f"{these.get('lem_inputs_figure', True)}: planned "
                         f"{planned} figures and built {len(made.figures())}")
    return fails


def test_search_figure_is_read_for_the_engineer():
    """The search subsection names what the search plot draws.

    Several thousand overlaid trials are unreadable to anyone who has not been
    told which mark is which, and the marks are not self-explanatory: a gray arc,
    a black dot and a green arrow are three different statements about the
    search. The prose names each, and names it as the artists actually draw it.
    """
    fails = []
    report = _build(dict({"method": ["spencer"], "lem_search_figure": True},
                         **FAST_FIGURES))
    search = next((s for _l, s in _sections(report)
                   if s.title == "Search for the Critical Surface"), None)
    if search is None:
        return ["the sample's report carries no search section"]
    figure = next((b for b in search.blocks if b.kind == "figure"), None)
    if figure is None:
        return ["the search section drew no figure, so its prose proves nothing"]
    said = " ".join(b.text for b in search.blocks if b.kind == "prose")
    if f"Figure {figure.number}" not in said:
        fails.append(f"the search figure is not cited: {said!r}")

    # Every element the circular search plot draws is named. The colours are the
    # ones plot.py hardcodes — nothing in the style names them — so the prose
    # can and must say which mark is which.
    for element in ("black dot", "red", "green arrow", "center"):
        if element not in said:
            fails.append(f"the search prose does not name the {element!r} the "
                         f"figure draws: {said!r}")
    if "refine" not in said:
        fails.append(f"the search prose never says the grid is refined: {said!r}")

    # The claim is checked against the artists, not against the wording: a plot
    # that stopped drawing centers would leave the sentence describing a figure
    # that is not there.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.figure as mplfig
    from xslope.plot import plot_circular_search_results
    slope_data, solutions = _solved()
    bundle = next(b for b in solutions["lem"] if b.get("search"))
    run = bundle["search"]
    f = mplfig.Figure()
    plot_circular_search_results(slope_data, run.get("fs_cache") or [],
                                 search_path=run.get("search_path"),
                                 circle_cache=run.get("circle_cache"), fig=f,
                                 show_title=False, show_legend=False)
    ax = f.axes[0]
    drawn = {a.get_gid() for a in ax.get_children() if a.get_gid()}
    if "CIRCLE_CENTERS" not in drawn:
        fails.append(f"the search plot draws no circle centers, and the prose "
                     f"says it does: {sorted(drawn)}")
    if "CRITICAL_SURFACE" not in drawn:
        fails.append(f"the search plot highlights no critical surface: "
                     f"{sorted(drawn)}")
    centers = [a for a in ax.get_children()
               if a.get_gid() == "CIRCLE_CENTERS" and hasattr(a, "get_color")]
    colors = {str(a.get_color()) for a in centers}
    if not {"k", "black"} & colors:
        fails.append(f"the trial centers are drawn {colors}, and the prose calls "
                     f"them black")
    if not {"r", "red"} & colors:
        fails.append(f"no center is drawn red, and the prose says the critical "
                     f"one is: {colors}")
    if run.get("search_path") and len(run["search_path"]) > 1:
        arrows = [a for a in ax.get_children() if a.get_gid() == "SEARCH_PATH"]
        if not arrows:
            fails.append("the search plot draws no path, and the prose describes "
                         "arrows running from stage to stage")
        elif not any("green" in str(a.get_facecolor()) or
                     a.get_facecolor()[:3] == (0.0, 0.5019607843137255, 0.0)
                     for a in arrows):
            fails.append(f"the search path is not green, and the prose calls it "
                         f"green: {[a.get_facecolor() for a in arrows]}")

    # A model with no circle cache says only what it drew: the centers, not the
    # arcs it did not keep.
    from xslope.report import _search_section, _Counter, resolve_options
    thin = dict(bundle, search=dict(run, circle_cache=None))
    tmp = tempfile.mkdtemp(prefix="xslope_searchprose_")
    with contextlib.redirect_stdout(io.StringIO()):
        sub = _search_section(slope_data, thin,
                              resolve_options({"input_path": REINF_XLSX}),
                              _Counter(), tmp, "spencer")
    thin_said = " ".join(b.text for b in sub.blocks if b.kind == "prose")
    if "drawn in gray" in thin_said:
        fails.append(f"a search that kept no trial circles still says they are "
                     f"drawn: {thin_said!r}")
    if "black dot" not in thin_said:
        fails.append(f"a search that kept no trial circles no longer names its "
                     f"centers: {thin_said!r}")
    return fails


def test_fs_summary_reaches_the_document():
    """The summary the tree holds is the summary the document prints: the
    documented methods, and no others.

    The rows used to be set in bold to pick the documented methods out of a table
    that also listed every method that had not been run. There is nothing left to
    pick out — every row is a method the report documents — so what the document
    owes the reader is that the table it prints is that table and not a longer
    one it inherited.
    """
    fails = []
    from docx import Document

    from xslope.report import generate_report, method_label, supported_methods

    slope_data, solutions = _solved()
    featured = ["bishop", "spencer"]
    expect = [method_label(m) for m in featured]

    block = _fs_table_block(_build({"method": featured}))
    if block is None:
        return ["there is no factor-of-safety summary"]
    if [r[0] for r in block.rows] != expect:
        fails.append(f"the tree's summary lists {[r[0] for r in block.rows]}, "
                     f"not {expect}")

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "summary.docx")
        ok, out = generate_report(
            slope_data, solutions,
            {"input_path": REINF_XLSX, "method": featured, "pd_figure": False,
             "lem_search_figure": False, "lem_solution_figure": False,
             "lem_slice_key": False}, path)
        if not ok:
            return fails + [f"generate_report failed: {out}"]

        doc = Document(path)
        table = next((t for t in doc.tables
                      if [c.text for c in t.rows[0].cells][:2]
                      == ["Method", "Factor of safety"]), None)
        if table is None:
            return fails + ["the summary table is not in the document"]
        printed = [row.cells[0].text for row in table.rows[1:]]
        if printed != expect:
            fails.append(f"the document's summary lists {printed}, not {expect}")
        if len(printed) >= len(supported_methods()):
            fails.append(f"the summary still enumerates the solver "
                         f"({len(printed)} rows for {len(featured)} documented "
                         f"methods)")
    return fails


def test_slice_key_figure():
    """Every slice table is preceded by a figure of the slices it lists, numbered.

    The table's first column is a slice number and nothing else in the report
    says where slice 9 is. The key is a figure of the sliced mass with the numbers
    on it, and it has to sit immediately BEFORE the table it keys — a key printed
    after the table is a key the reader has already needed.
    """
    fails = []
    from xslope.report import method_label

    tmp = tempfile.mkdtemp(prefix="xslope_key_")
    methods = ["spencer", "bishop"]
    report = _build({"method": methods, "lem_solution_figure": False},
                    figure_dir=tmp, fast=False)

    for m in methods:
        label = method_label(m)
        keys = [f for f in report.figures() if f.source == f"{m} slice key"]
        if len(keys) != 1:
            fails.append(f"{m}: {len(keys)} slice-key figures")
            continue
        key = keys[0]
        if label not in key.caption:
            fails.append(f"{m}: the slice key is captioned {key.caption!r}")
        if not os.path.exists(key.path):
            fails.append(f"{m}: the slice key was not written to {key.path}")
        elif os.path.getsize(key.path) < 20000:
            fails.append(f"{m}: the slice key is {os.path.getsize(key.path)} "
                         f"bytes — too small to be a rendered plot")
        # A portrait figure at text width, the size every other plot is. The
        # landscape page belongs to the table's twenty columns; a figure that
        # took it too spent a sheet on a picture that reads at a sixth of one.
        if key.landscape or key.width_in <= 0:
            fails.append(f"{m}: the slice key takes a landscape page instead of "
                         f"printing at text width (landscape={key.landscape}, "
                         f"width_in={key.width_in})")
        from xslope.report import Figure as ReportFigure
        text_width = ReportFigure("", "").width_in
        if key.width_in != text_width:
            fails.append(f"{m}: the slice key is {key.width_in} in wide, not the "
                         f"{text_width} in every other figure prints at")

    # Immediately before the table, in the same section, with nothing between.
    for sec in report.sections:
        for _lvl, node in sec.walk():
            if node.title != "Slice Table":
                continue
            kinds = [b.kind for b in node.blocks]
            if kinds[-2:] != ["figure", "table"]:
                fails.append(f"the Slice Table section runs {kinds}; the key must "
                             f"be the block immediately before the table")

    # And it is the slice-key option that puts it there.
    off = _build({"method": methods, "lem_slice_key": False})
    if [f for f in off.figures() if "slice key" in f.source]:
        fails.append("lem_slice_key=False still drew a slice key")
    if not [t for t in off.tables() if t.landscape]:
        fails.append("lem_slice_key=False took the slice table with it")

    # The frame is the SLICED MASS, not the model: the key is tighter than the
    # solution plot of the same surface, which is the whole point of it.
    from xslope.plot import sliced_mass_bounds
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MplFigure
    from xslope.plot import plot_solution

    slope_data, solutions = _solved()
    bundle = solutions["lem"][0]
    spans = {}
    for frame in ("fill", "slices"):
        fig = MplFigure(figsize=(9.0, 5.5))
        FigureCanvasAgg(fig)
        plot_solution(slope_data, bundle["slice_df"], bundle["failure_surface"],
                      bundle["results"], fig=fig, show_title=False,
                      slice_numbers=(frame == "slices"), frame=frame)
        ax = fig.axes[0]
        spans[frame] = (ax.get_xlim(), ax.get_ylim(), sliced_mass_bounds(ax),
                        [t for t in ax.texts if t.get_gid() == "SLICE_NUMBER"])
    (fx0, fx1), _fy, _fb, fill_labels = spans["fill"]
    (sx0, sx1), (sy0, sy1), bounds, labels = spans["slices"]
    if fill_labels:
        fails.append("frame='fill' drew slice numbers without being asked")
    if len(labels) != len(bundle["slice_df"]):
        fails.append(f"{len(labels)} slice labels for "
                     f"{len(bundle['slice_df'])} slices")
    if not (sx1 - sx0) < (fx1 - fx0):
        fails.append(f"the slice frame ({sx0:.1f}, {sx1:.1f}) is no tighter than "
                     f"the model frame ({fx0:.1f}, {fx1:.1f})")
    if bounds is None:
        fails.append("the sliced mass has no measurable bounds")
    else:
        bx0, bx1, by0, by1 = bounds
        # Everything the sliced mass draws is inside the frame — the base-stress
        # bars stand off the base and are the reason a frame on the slice corners
        # alone would clip.
        if not (sx0 <= bx0 and bx1 <= sx1 and sy0 <= by0 and by1 <= sy1):
            fails.append(f"the frame ({sx0:.1f}..{sx1:.1f}, {sy0:.1f}..{sy1:.1f}) "
                         f"clips the sliced mass ({bx0:.1f}..{bx1:.1f}, "
                         f"{by0:.1f}..{by1:.1f})")
        # And the cushion is uniform and modest, not a slab of white.
        pads = (bx0 - sx0, sx1 - bx1, by0 - sy0, sy1 - by1)
        size = max(bx1 - bx0, by1 - by0)
        if max(pads) > 0.10 * size:
            fails.append(f"the cushion is {max(pads) / size:.0%} of the mass; the "
                         f"frame is not tight")
        if max(pads) - min(pads) > 1e-6 * size:
            fails.append(f"the cushion is not uniform: {pads}")

    # Labels never overlap, and they are big enough to read ON THE PAGE.
    #
    # The key is drawn at the report's figure size and printed at the text width
    # of a portrait page, so type in the figure reaches the page in that ratio.
    # Both models are checked: the sample's fifteen slices, and the dam's forty,
    # where the crest carries a pair of slivers a fifth of a slice wide — one
    # tight pair used to take the whole key down to three points printed.
    from xslope.report import FIGURE_SIZE
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    printed = _TEXT_WIDTH_IN / FIGURE_SIZE[0]
    with contextlib.redirect_stdout(io.StringIO()):
        dam = load_slope_data_cached(SEEP_XLSX)
        ok, out = generate_slices(dam, circle=dam["circles"][0], num_slices=40)
        dam_df = out[0] if ok else None
        dam_surface = out[1] if ok else None
        dam_res = solve_selected("spencer", dam_df) if ok else None
    cases = [("the sample", slope_data, bundle["slice_df"],
              bundle["failure_surface"], bundle["results"])]
    if dam_df is not None:
        cases.append(("the dam", dam, dam_df, dam_surface, dam_res))
    else:
        fails.append("the dam produced no slices; the crowded key is untested")

    for name, sd, df, surface, results in cases:
        fig = MplFigure(figsize=FIGURE_SIZE)
        FigureCanvasAgg(fig)
        with contextlib.redirect_stdout(io.StringIO()):
            plot_solution(sd, df, surface, results, fig=fig, show_title=False,
                          slice_numbers=True, frame="slices")
        fig.canvas.draw()
        ax = fig.axes[0]
        renderer = fig.canvas.get_renderer()
        marks = [t for t in ax.texts if t.get_gid() == "SLICE_NUMBER"]
        boxes = [t.get_window_extent(renderer) for t in marks]
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                if boxes[i].overlaps(boxes[j]):
                    fails.append(f"{name}: slice numbers "
                                 f"{marks[i].get_text()} and "
                                 f"{marks[j].get_text()} are printed over "
                                 f"each other")
                    break
            else:
                continue
            break
        smallest = min((t.get_fontsize() for t in marks), default=0.0)
        if smallest * printed < _LEGIBLE_PT:
            fails.append(f"{name}: the slice numbers are set at "
                         f"{smallest:.1f} pt, which is "
                         f"{smallest * printed:.1f} pt on the page — under the "
                         f"{_LEGIBLE_PT} pt a number has to be to be read")
    return fails


def test_figure_progress_counts():
    """The figure count a caller is shown is the number of figures there are.

    N methods mean N critical-surface plots and N slice keys; a progress line
    that said "3" while eleven were rendering would be worse than none.
    """
    fails = []
    from xslope.report import planned_figures, resolve_options

    seen = []
    tmp = tempfile.mkdtemp(prefix="xslope_prog_")
    opts = {"method": ["spencer", "bishop"],
            "progress": lambda done, total, label: seen.append((done, total, label))}
    report = _build(opts, figure_dir=tmp, fast=False)

    slope_data, solutions = _solved()
    planned = planned_figures(slope_data, solutions,
                              resolve_options({"input_path": REINF_XLSX,
                                               "title": "Sample Levee", **opts}))
    drawn = len(report.figures())
    if planned != drawn:
        fails.append(f"{planned} figures were planned and {drawn} were built")
    if [d for d, _t, _l in seen] != list(range(1, len(seen) + 1)):
        fails.append(f"the progress steps are {[d for d, _t, _l in seen]}")
    if seen and seen[-1][0] != seen[-1][1]:
        fails.append(f"the last step is {seen[-1][0]} of {seen[-1][1]}")
    if {t for _d, t, _l in seen} != {planned}:
        fails.append(f"the total reported was {[t for _d, t, _l in seen]}, not "
                     f"{planned}")
    labels = [l for _d, _t, l in seen]
    if len([l for l in labels if "slice key" in l]) != 2:
        fails.append(f"the two slice keys are not both counted: {labels}")
    return fails


def test_fs_table_compares_only_what_was_documented():
    """The summary lists the methods the report documents, each with its own
    answer, and only where there is something to compare.

    It used to list every method xslope offers, filling in the ones that had not
    been run by solving them on the surface the featured method had found. Each
    method has its own critical surface, so those rows were not that method's
    factor of safety — a column of numbers that read as a comparison and was not
    one. What the reader wants compared is what the reader asked the report to
    document, so those are the rows, and a report of ONE method carries no table
    at all: one row is the number that method's own section already states, and a
    second statement of a number is a second number to keep in step with it.
    """
    fails = []
    from xslope.report import build_report, method_label

    slope_data, solutions = _solved()

    # One method: no table, and the factor of safety is still stated.
    single = _build({"method": "spencer"})
    if _fs_rows(single) is not None:
        fails.append("a single-method report still carries a comparison table")
    if "Factors of Safety" in [t for _l, t in single.section_titles()]:
        fails.append("a single-method report still has a Factors of Safety "
                     "section")
    fs = solutions["lem"][0]["results"]["FS"]
    if not any(f"{fs:.3f}" in b.text for b in single.blocks("prose")):
        fails.append(f"a single-method report states its factor of safety "
                     f"{fs:.3f} nowhere")

    # Two methods: a table of exactly those two, each carrying its own answer.
    pair = _build({"method": ["spencer", "bishop"]})
    got = _fs_rows(pair)
    if got is None:
        return fails + ["a two-method report carries no summary table"]
    _headers, rows = got
    want = {method_label(b["method"]): b["results"]["FS"]
            for b in solutions["lem"]}
    if [r[0] for r in rows] != list(want):
        fails.append(f"the summary lists {[r[0] for r in rows]}, not the "
                     f"documented methods {list(want)}")
    for row in rows:
        expect = want.get(row[0])
        if expect is None:
            fails.append(f"{row[0]} is in the summary and is not documented")
        elif row[1] != f"{expect:.3f}":
            fails.append(f"{row[0]}: the summary says {row[1]!r}, its own run "
                         f"gave {expect:.3f}")

    # And every row says which surface its number belongs to. The sample's
    # bundles carry a search on one method and none on the other, so both
    # wordings are exercised on one table.
    surfaces = [r[2].lower() for r in rows]
    if not any("critical" in s for s in surfaces):
        fails.append(f"no row names a searched critical surface: {surfaces}")
    if not any("specified" in s for s in surfaces):
        fails.append(f"no row names a specified surface: {surfaces}")

    # A method that does not converge is stated, not dropped.
    bundle = dict(solutions["lem"][0])
    df = bundle["slice_df"].copy()
    df["w"] = float("nan")             # no method can produce an answer from this
    bad = {"lem": [{"slice_df": df, "failure_surface": bundle["failure_surface"],
                    "results": {}, "search": None, "method": "spencer"},
                   {"slice_df": df, "failure_surface": bundle["failure_surface"],
                    "results": {}, "search": None, "method": "bishop"}]}
    with tempfile.TemporaryDirectory() as tmp:
        broken = build_report(slope_data, bad,
                              {"pd_figure": False, "lem_solution_figure": False,
                               "lem_search_figure": False,
                               "method": ["spencer", "bishop"]}, tmp)
    got = _fs_rows(broken)
    if got is None:
        fails.append("a report whose methods all failed lost its summary table")
    else:
        _h, rows = got
        if not [r for r in rows if "did not converge" in r[1]]:
            fails.append(f"no method reported 'did not converge' on an "
                         f"unsolvable slice table: {rows}")
        if len(rows) != 2:
            fails.append(f"the unsolvable run dropped rows: {len(rows)} of 2")
    return fails


# --------------------------------------------------------------------------
# D. the document
# --------------------------------------------------------------------------

def _docx_parts(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        return names, {n: z.read(n).decode("utf-8", "replace")
                       for n in names if n.endswith(".xml")}


def test_docx():
    """The written .docx is a real OOXML package with the structure the report
    depends on."""
    fails = []
    from xslope.report import generate_report

    slope_data, solutions = _solved()
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "report.docx")
        ok, out = generate_report(
            slope_data, solutions,
            {"input_path": REINF_XLSX, "title": "Sample Levee Report",
             "project_number": "26-001", "organization": "Example Engineering",
             "author": "A. Engineer", "method": "spencer",
             "signature_lines": True, "lem_search_figure": False},
            out_path)
        if not ok:
            return [f"generate_report failed: {out}"]
        if not os.path.exists(out_path):
            return ["generate_report reported success but wrote no file"]

        names, xml = _docx_parts(out_path)
        doc = xml.get("word/document.xml", "")
        if not doc:
            return ["the package has no word/document.xml"]

        if "Sample Levee Report" not in doc:
            fails.append("the document does not contain the project title")
        if "A. Engineer" not in doc:
            fails.append("the document does not contain the author")
        if 'DOCPROPERTY "Title"' not in doc:
            fails.append("the title page does not use a document-property field")
        if "TOC \\o" not in doc:
            fails.append("the document has no table-of-contents field")
        if 'w:dirty="true"' in doc:
            fails.append("a field is marked dirty; Word would prompt about "
                         "updating fields on every open")
        if "Prepared by" not in doc or "Checked by" not in doc:
            fails.append("signature_lines=True produced no signature lines")

        # The template's styles are the ones referenced.
        for style in ("Heading1", "Caption", "TableGrid"):
            if style not in doc:
                fails.append(f"the document references no {style} style")

        # The slice table gets its own landscape section, and the report returns
        # to portrait after it.
        orients = [chunk.split('w:orient="')[1].split('"')[0]
                   if 'w:orient="' in chunk else "portrait"
                   for chunk in doc.split("<w:pgSz ")[1:]]
        if orients.count("landscape") != 1:
            fails.append(f"page orientations are {orients}; exactly one landscape "
                         f"section is expected")
        if not orients or orients[-1] != "portrait":
            fails.append(f"the report ends on a {orients[-1:]} page")

        # Running head and foot, with live fields.
        header = next((v for k, v in xml.items() if k.startswith("word/header")), "")
        footer = next((v for k, v in xml.items() if k.startswith("word/footer")), "")
        if 'DOCPROPERTY "Title"' not in header:
            fails.append("the running head carries no title field")
        if "PAGE" not in footer or "NUMPAGES" not in footer:
            fails.append("the footer has no page-of-pages fields")

        # The figures actually went in.
        media = [n for n in names if n.startswith("word/media/")]
        if len(media) < 2:
            fails.append(f"only {len(media)} images were embedded")

        # And the returned bundle names what was produced.
        if out["path"] != out_path:
            fails.append(f"generate_report returned {out['path']!r}")
        if not out["figures"]:
            fails.append("generate_report reported no figures")

    # An unavailable format is refused, in words, rather than half-written.
    with tempfile.TemporaryDirectory() as tmp:
        ok, msg = generate_report(slope_data, solutions, {},
                                  os.path.join(tmp, "r.pdf"))
        if ok:
            fails.append("a PDF report was reported as written; PDF is S3")
        elif "not available" not in str(msg):
            fails.append(f"the PDF refusal reads {msg!r}")
        if os.listdir(tmp):
            fails.append(f"the refused format left files behind: {os.listdir(tmp)}")
    return fails


def _toc_result(doc_xml):
    """The text of the contents field's cached result, entry by entry.

    Everything between the field's ``separate`` and its ``end`` is what a reader
    sees before Word ever updates the field — so this is exactly what is on the
    contents page of a document nobody has pressed F9 in.
    """
    import re
    instr = doc_xml.index("TOC \\o")
    start = doc_xml.index('w:fldCharType="separate"', instr)
    end = doc_xml.index('w:fldCharType="end"', start)
    return re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc_xml[start:end])


def _sections_usable(doc_xml):
    """Every section's usable text width in twips, with the position in the XML
    where the section ends.

    A table belongs to the first section that ends after it — that is how Word
    reads a document, and it is what says whether the slice table had a landscape
    page's width to fill.
    """
    import re
    out = []
    for m in re.finditer(r"<w:sectPr[ >].*?</w:sectPr>", doc_xml, re.S):
        sect = m.group(0)
        page = re.search(r'<w:pgSz[^>]*\bw:w="(\d+)"', sect)
        margins = re.search(r"<w:pgMar[^>]*/>", sect)
        if page is None or margins is None:
            continue
        left = int(re.search(r'w:left="(\d+)"', margins.group(0)).group(1))
        right = int(re.search(r'w:right="(\d+)"', margins.group(0)).group(1))
        out.append((m.start(), int(page.group(1)) - left - right))
    return out


def _cell_margin(tbl_xml, styles_xml, style_id):
    """The leading cell margin a table carries, in twips.

    Its own declaration where it makes one — the report states its padding on
    every table it writes, because the template's is a body table's and starves a
    twenty-two column slice table of the width its numbers need — and the table
    style's otherwise. Whichever it is, it is the number the indent and the
    column measurement both have to be built on.
    """
    import re
    own = re.search(r"<w:tblCellMar>.*?<w:left [^>]*w:w=\"(\d+)\"",
                    tbl_xml, re.S)
    if own:
        return int(own.group(1))
    style = re.search(rf'<w:style [^>]*w:styleId="{style_id}".*?</w:style>',
                      styles_xml, re.S)
    if style is None:
        return None
    left = re.search(r"<w:tblCellMar>.*?<w:left [^>]*w:w=\"(\d+)\"",
                     style.group(0), re.S)
    return int(left.group(1)) if left else None


def _table_cell_texts(tbl_xml):
    """Everything that prints in each column of a table, one list per column.

    The header row, the body rows and the totals row — the same strings the
    renderer measured its columns from, read back out of the document.
    """
    import re
    columns = []
    for row in re.findall(r"<w:tr[ >].*?</w:tr>", tbl_xml, re.S):
        for j, cell in enumerate(re.findall(r"<w:tc>.*?</w:tc>", row, re.S)):
            while len(columns) <= j:
                columns.append([])
            columns[j].append("".join(
                re.findall(r"<w:t[^>]*>([^<]*)</w:t>", cell)))
    return columns


def _table_point_size(tbl_xml):
    """The point size the table's cells are set in, as the document states it.

    Word stores a run's size in half-points; the size that appears most often in
    the table is the one its cells carry.
    """
    import re
    sizes = [int(v) for v in re.findall(r'<w:sz w:val="(\d+)"/>', tbl_xml)]
    return max(set(sizes), key=sizes.count) / 2.0 if sizes else None


def _content_widths(columns, family, size_pt, usable, margin):
    """What each column's content needs, in twips.

    Its widest line plus the cell margins either side, floored at its longest
    single word — Word breaks a word that will not fit rather than widening the
    column — and that floor capped at an equal share of the page, past which
    refusing to break would starve every other column. An empty column still
    reads as a column, so it floors at one em.
    """
    from xslope.report_docx import TWIPS_PER_PT, _text_width

    pad = 2 * margin
    fair_share = usable / max(1, len(columns))
    out = []
    for texts in columns:
        widest = max((_text_width(t, family, size_pt) for t in texts), default=0.0)
        longest_word = max((_text_width(w, family, size_pt)
                            for t in texts for w in t.split()), default=0.0)
        out.append(max(widest + pad,
                       min(max(longest_word, size_pt * TWIPS_PER_PT) + pad,
                           fair_share)))
    return out


def test_table_geometry():
    """Every table is laid out to its content, on the page it sits on.

    Three defects this defends against, none of which Word fixes on its own: a
    table with no indent hangs its left border out into the margin; an
    autofitting table gives "#" the same width as "Material"; and a table ruled
    across the full text width when it holds three columns of factors of safety
    is a table pretending to be a page. The cure is a fixed layout with measured
    columns, an indent of exactly one cell margin, and a width that is the
    content's — cut down to the text width only where the content would overrun
    it.
    """
    import re
    fails = []
    from xslope.report import generate_report
    from xslope.report_docx import DEFAULT_CELL_MARGIN

    slope_data, solutions = _solved()
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "report.docx")
        # The model checks are on: their Finding column is the long-text column
        # the widths have to wrap rather than let it starve its neighbours.
        ok, out = generate_report(
            slope_data, solutions,
            {"input_path": REINF_XLSX, "title": "Sample Levee Report",
             "project_number": "26-001", "author": "A. Engineer",
             "method": "spencer", "signature_lines": True, "model_checks": True,
             "pd_figure": False, "lem_search_figure": False,
             "lem_solution_figure": False},
            out_path)
        if not ok:
            return [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(out_path)
        doc = xml.get("word/document.xml", "")
        styles = xml.get("word/styles.xml", "")
        # The face the columns were measured in, read off the document's own
        # Normal style — the widths mean nothing without it.
        from docx import Document

        from xslope.report_docx import _table_font
        family = _table_font(Document(out_path))

    sections = _sections_usable(doc)
    if not sections:
        return ["the document declares no page size to fit the tables to"]
    tables = list(re.finditer(r"<w:tbl>.*?</w:tbl>", doc, re.S))
    if len(tables) < 5:
        return [f"only {len(tables)} tables were written; the report has more"]

    for i, m in enumerate(tables):
        tbl = m.group(0)
        tbl_pr = re.search(r"<w:tblPr>.*?</w:tblPr>", tbl, re.S).group(0)
        style = re.search(r'<w:tblStyle w:val="([^"]+)"', tbl_pr)
        where = f"table {i + 1} ({style.group(1) if style else 'unstyled'})"
        usable = next((u for pos, u in sections if pos > m.start()),
                      sections[-1][1])

        if '<w:tblLayout w:type="fixed"/>' not in tbl_pr:
            fails.append(f"{where} is not fixed-layout; Word will autofit it and "
                         f"give every column the same width")

        # The indent belongs to a BORDERED table, whose border is the edge the
        # reader sees and belongs on the text margin. A borderless table's edge is
        # its text, and that belongs where a line of prose begins — so it carries
        # no indent, which is what leaves its text flush with the body.
        margin = _cell_margin(tbl_pr, styles,
                              style.group(1) if style else "TableNormal")
        indent = re.search(r"<w:tblInd [^>]*/>", tbl_pr)
        if style is None:
            if indent is not None and int(
                    re.search(r'w:w="(-?\d+)"', indent.group(0)).group(1)):
                fails.append(f"{where} is indented {indent.group(0)}; a "
                             f"borderless table's text would then start inside "
                             f"the body text edge")
        elif indent is None:
            fails.append(f"{where} carries no indent; its left border will sit in "
                         f"the margin, left of the body text")
        elif margin is None:
            fails.append(f"{where} names a style with no cell margin to match")
        else:
            got = re.search(r'w:w="(-?\d+)"', indent.group(0))
            kind = re.search(r'w:type="(\w+)"', indent.group(0))
            if got is None or int(got.group(1)) != margin or (
                    kind is None or kind.group(1) != "dxa"):
                fails.append(f"{where} is indented {indent.group(0)}, not "
                             f"{margin} dxa — the cell margin its border "
                             f"overhangs by")

        grid = [int(w) for w in re.findall(r'<w:gridCol w:w="(\d+)"/>', tbl)]
        if not grid:
            fails.append(f"{where} declares no column grid")
            continue
        if sum(grid) > usable:
            fails.append(f"{where} spans {sum(grid)} twips of a {usable}-twip "
                         f"text width; it is wider than the page")
        if min(grid) <= 0:
            fails.append(f"{where} has a column of {min(grid)} twips")

        # The table declares the width its own columns come to, not the page's:
        # that is the number Word stops the last column at.
        declared = re.search(r'<w:tblW [^>]*w:w="(-?\d+)"[^>]*w:type="(\w+)"',
                             tbl_pr) or re.search(
            r'<w:tblW [^>]*w:type="(\w+)"[^>]*w:w="(-?\d+)"', tbl_pr)
        if declared is None:
            fails.append(f"{where} declares no width of its own")
        else:
            width, kind = declared.groups()
            if not width.lstrip("-").isdigit():
                width, kind = kind, width
            if int(width) != sum(grid) or kind != "dxa":
                fails.append(f"{where} declares a width of {width} {kind}, not "
                             f"the {sum(grid)} dxa its grid comes to")

        # The width IS the content's. Each column is measured in the document's
        # own face at the size the document states, so a column stretched past
        # what it holds shows up as a column wider than its text, and a table
        # padded out to the margin shows up as every column at once.
        #
        # The signature blocks are the one table that asks for the whole page —
        # two of them belong one at each margin, not side by side in the middle —
        # so that one is required to span it and the rest are required not to.
        columns = _table_cell_texts(tbl)
        size = _table_point_size(tbl)
        if "Prepared by" in tbl:
            if abs(sum(grid) - usable) > 1:
                fails.append(f"the signature blocks span {sum(grid)} twips of a "
                             f"{usable}-twip text width; they belong one at each "
                             f"margin")
        elif len(columns) != len(grid) or size is None:
            fails.append(f"{where}: {len(columns)} columns of text for "
                         f"{len(grid)} grid columns at {size} pt")
        else:
            needs = _content_widths(
                columns, family, size, usable,
                DEFAULT_CELL_MARGIN if margin is None else margin)
            for j, width in enumerate(grid):
                if width > needs[j] + 1:
                    fails.append(f"{where} gives column {j} {width} twips for "
                                 f"content {needs[j]:.0f} wide; the table is "
                                 f"being stretched rather than fitted")
            want = min(int(round(sum(needs))), usable)
            if sum(grid) < want - 1:
                fails.append(f"{where} spans {sum(grid)} twips where its content "
                             f"needs {want}; it is being starved")

            # A NUMBER is never printed on two lines. Word breaks a line after a
            # hyphen-minus, so a column a twip too narrow prints "-" over
            # "1416.5" and a total as "78357." over "0" — which is what the
            # slice table did. Two things have to hold at once: the column is at
            # least as wide as the widest value in it, and every cell holding a
            # number says it may not be broken. The first alone is a measurement
            # that a layout is free to disagree with; the second alone would
            # overflow a column too narrow to hold what it forbids breaking.
            from xslope.report_docx import _text_width as _tw
            from xslope.columns import is_number
            for j, texts in enumerate(columns):
                values = [t for t in texts if is_number(t)]
                if not values or any(t.strip() and not is_number(t)
                                     for t in texts[1:]):
                    continue
                widest = max(_tw(t, family, size) for t in values)
                if grid[j] < widest + 2 * margin - 1:
                    fails.append(f"{where} gives column {j} {grid[j]} twips for "
                                 f"a value {widest + 2 * margin:.0f} wide; the "
                                 f"number will be broken across two lines")
            numbers = sum(1 for texts in columns for t in texts if is_number(t))
            if tbl.count("<w:noWrap/>") != numbers:
                fails.append(f"{where} holds {numbers} numbers but marks "
                             f"{tbl.count('<w:noWrap/>')} cells unbreakable")

        # Word wants the widths in the cells as well as in the grid.
        for row in re.findall(r"<w:tr>.*?</w:tr>", tbl, re.S):
            cells = [int(w) for w in
                     re.findall(r'<w:tcW w:type="dxa" w:w="(\d+)"/>', row)]
            if cells != grid:
                fails.append(f"{where} has a row whose cell widths {cells} do "
                             f"not match its grid {grid}")
                break

        # The header row is a header: bold, and repeated on every page it spans.
        if style is not None:
            head = re.search(r"<w:tr>.*?</w:tr>", tbl, re.S).group(0)
            if "<w:tblHeader" not in head:
                fails.append(f"{where} does not repeat its header row")
            if "<w:b/>" not in head:
                fails.append(f"{where} has no bold in its header row")

    # The columns are measured, not equal: "#" is narrower than "Material", and
    # the long Finding column is wide without taking everything.
    materials = next((t.group(0) for t in tables if "Material" in t.group(0)
                      and "Strength" in t.group(0)), "")
    grid = [int(w) for w in re.findall(r'<w:gridCol w:w="(\d+)"/>', materials)]
    if len(grid) > 2 and not grid[0] < grid[1]:
        fails.append(f"the materials table gives '#' {grid[0]} twips and "
                     f"'Material' {grid[1]}: the columns are not measured")
    findings = next((t.group(0) for t in tables if "Severity" in t.group(0)), "")
    grid = [int(w) for w in re.findall(r'<w:gridCol w:w="(\d+)"/>', findings)]
    if len(grid) == 3:
        from xslope.report_docx import CELL_MARGIN, TABLE_PT, _text_width
        if grid[1] != max(grid):
            fails.append(f"the long Finding column is not the widest: {grid}")
        # It is wide by wrapping, not by starving: "Warning" beside it still
        # prints on one line.
        needed = _text_width("Warning", "Calibri", TABLE_PT) + 2 * CELL_MARGIN[0]
        if grid[0] < needed:
            fails.append(f"the Severity column is {grid[0]} twips, under the "
                         f"{needed:.0f} 'Warning' needs — the Finding column "
                         f"starved it")
    elif not findings:
        fails.append("the model-check findings table was not written")
    return fails


def test_section_breaks_take_no_room():
    """A section break can never manufacture a page of its own.

    Word writes one as a paragraph, and that paragraph is a real empty line at
    the end of the section it closes. It is invisible until the section's content
    ends near the foot of a page — which a table spanning pages does sooner or
    later — and then it has nowhere to go but a fresh page, and the report grows
    a blank landscape sheet between the slice table and the section after it.

    So every paragraph that carries one is made to vanish: the smallest size the
    format expresses, an exact line height to match, no space above or below, and
    a hidden paragraph mark. Checked on all of them, not on the one that showed.
    """
    import re
    fails = []
    from xslope.report import generate_report
    from xslope.report_docx import SECT_BREAK_PT

    slope_data, solutions = _solved()
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "report.docx")
        ok, out = generate_report(
            slope_data, solutions,
            {"input_path": REINF_XLSX, "method": "spencer", "pd_figure": False,
             "lem_search_figure": False, "lem_solution_figure": False},
            out_path)
        if not ok:
            return [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(out_path)
    doc = xml.get("word/document.xml", "")

    carriers = [p for p in re.findall(r"<w:p\b[^>]*>.*?</w:p>", doc, re.S)
                if "<w:sectPr" in p]
    if not carriers:
        return ["the report opens no section of its own; there is no section "
                "break to check and the landscape slice table is missing"]
    half_points = SECT_BREAK_PT * 2
    for i, p in enumerate(carriers):
        where = f"section break {i + 1} of {len(carriers)}"
        p_pr = re.search(r"<w:pPr>.*?</w:pPr>", p, re.S)
        r_pr = re.search(r"<w:rPr>.*?</w:rPr>", p_pr.group(0) if p_pr else "",
                         re.S)
        if r_pr is None or "<w:vanish/>" not in r_pr.group(0):
            fails.append(f"{where} does not hide its paragraph mark; an empty "
                         f"line will print where it sits")
            continue
        size = re.search(r'<w:sz w:val="(\d+)"/>', r_pr.group(0))
        if size is None or int(size.group(1)) > half_points:
            fails.append(f"{where} sets its paragraph mark at "
                         f"{size.group(1) if size else 'the body size'} "
                         f"half-points, over the {half_points} it may be")
        spacing = re.search(r"<w:spacing [^>]*/>", p_pr.group(0))
        if spacing is None:
            fails.append(f"{where} names no spacing; it keeps the style's")
            continue
        for attr in ("before", "after"):
            got = re.search(rf'w:{attr}="(\d+)"', spacing.group(0))
            if got is None or int(got.group(1)) != 0:
                fails.append(f"{where} leaves {attr}-spacing "
                             f"{got.group(1) if got else 'unset'} on a line "
                             f"that is meant to take no room")
        line = re.search(r'w:line="(\d+)"', spacing.group(0))
        rule = re.search(r'w:lineRule="(\w+)"', spacing.group(0))
        if line is None or rule is None or rule.group(1) != "exact" or \
                int(line.group(1)) > half_points * 10:
            fails.append(f"{where} sets no exact line height; a hidden mark in "
                         f"a body-height line is still a body-height line")
    return fails


#: Where LibreOffice is looked for. It renders the written document to PDF, and
#: the text drawn on those pages is what a reader can see — the only evidence
#: that a paragraph of the tree reached a page. Without it that leg is skipped
#: and the structural one still runs.
SOFFICE = (os.environ.get("XSLOPE_SOFFICE"), "soffice",
           "/Applications/LibreOffice.app/Contents/MacOS/soffice")

NONCIRC_XLSX = os.path.join(_REPO, "docs", "lem", "files",
                            "xslope_noncircular.xlsx")


def _soffice():
    """LibreOffice, or None on a machine without it."""
    import shutil
    for cand in SOFFICE:
        if cand and (shutil.which(cand) or os.path.isfile(cand)):
            return shutil.which(cand) or cand
    return None


def _rendered_text(path, title):
    """The text of a written report as it comes out on the page, or None when
    this machine cannot render one.

    The document is converted to PDF and the text read off the pages, so a block
    the renderer dropped is missing from the result. The running head and foot go
    with it, and the wrapped lines are rejoined, which puts a paragraph back
    together across a page break it was broken over.
    """
    import subprocess
    soffice = _soffice()
    if soffice is None:
        return None
    try:
        import pypdf
    except Exception:
        return None
    try:
        subprocess.run([soffice, "--headless", "--convert-to", "pdf",
                        "--outdir", os.path.dirname(path), path],
                       capture_output=True, timeout=300)
    except Exception:
        return None
    pdf = os.path.splitext(path)[0] + ".pdf"
    if not os.path.exists(pdf):
        return None
    lines = []
    for page in pypdf.PdfReader(pdf).pages:
        for line in page.extract_text().splitlines():
            line = line.strip()
            if not line or line == title or re.match(r"^Page \d+ of \d+", line):
                continue
            lines.append(line)
    return re.sub(r"\s+", " ", " ".join(lines))


def _refused_solutions():
    """A model, and a solution whose working the report has to refuse.

    Janbu's method is solved on a copy of the slices that is then discarded, so
    the report is handed the unsolved frame beside the converged factor of
    safety. The equation evaluated on those values does not return the solution,
    the calculation is refused, and the refusal sentence stands as the last block
    of the landscape section the slice table opens — where a section break
    written onto a paragraph of content swallows it.
    """
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    slope_data = load_slope_data(NONCIRC_XLSX)
    ok, out = generate_slices(slope_data, non_circ=slope_data["non_circ"],
                              num_slices=15)
    if not ok:
        raise RuntimeError(f"the non-circular model produced no slices: {out}")
    slice_df, surface = out
    solved = slice_df.copy()
    with contextlib.redirect_stdout(io.StringIO()):
        results = solve_selected("janbu", solved)
    bundle = {"slice_df": slice_df.copy(), "failure_surface": surface,
              "results": results, "search": None, "method": "janbu"}
    return slope_data, {"lem": [bundle]}


def _paragraph_text(p):
    """What a paragraph of the document puts on the page — the text of its runs,
    and a marker for a picture, which a text scan would otherwise miss."""
    text = "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S))
    if "<w:drawing>" in p:
        text += "[picture]"
    return text.strip()


def _collapsed_with_content(doc, where):
    """Failures for every paragraph that is collapsed to a hidden mark and still
    holds something — a section break's carrier, or any other hidden mark."""
    fails = []
    for p in re.findall(r"<w:p\b[^>]*>.*?</w:p>", doc, re.S):
        p_pr = re.search(r"<w:pPr>.*?</w:pPr>", p, re.S)
        if p_pr is None:
            continue
        if "<w:sectPr" in p_pr.group(0):
            why = "carries a section break"
        elif "<w:vanish/>" in p_pr.group(0):
            why = "is marked hidden"
        else:
            continue
        text = _paragraph_text(p)
        if text:
            fails.append(f"{where}: a paragraph that {why} holds {text[:90]!r}. "
                         f"Collapsed to a hidden mark at an exact line height of "
                         f"a point, that never prints")
    return fails


def _missing_from_the_page(report, text, where):
    """Failures for every prose block of the tree that is not in the rendered
    text."""
    fails = []
    for block in report.blocks("prose"):
        want = re.sub(r"\s+", " ", block.text).strip()
        if want and want not in text:
            fails.append(f"{where}: a paragraph of the report is not on a page "
                         f"of it: {block.text!r}")
    return fails


def _written(where, slope_data, solutions, opts, tmp):
    """``(report, document.xml, path)`` for one report written to ``tmp``."""
    from xslope.report import generate_report
    path = os.path.join(tmp, "report.docx")
    with contextlib.redirect_stdout(io.StringIO()):
        ok, out = generate_report(slope_data, solutions, dict(opts), path)
    if not ok:
        raise RuntimeError(f"{where}: generate_report failed: {out}")
    _names, xml = _docx_parts(path)
    return out["report"], xml.get("word/document.xml", ""), path


#: The reports this check writes and reads back: the sample model, and the one
#: whose landscape section ends in a sentence.
_PAGE_REPORTS = (
    ("the sample report", REINF_XLSX,
     {"title": "Sample Levee", "method": "spencer"}),
    ("a refused calculation", NONCIRC_XLSX,
     {"title": "Non-circular Surface", "method": "janbu"}),
)


def test_every_block_reaches_the_page():
    """Every paragraph of the report is on a page of the report.

    A block can be in the tree, be written into the document, and still not reach
    a reader. The paragraph that carries a closing section break is collapsed to
    a hidden mark at an exact line height of a point, and a sentence set in a
    paragraph like that is cropped to nothing — which is what becomes of the
    sentence a refused calculation prints when it is the last block of the
    landscape section the slice table opens.

    Two things are required of every report written here. Nothing is collapsed
    but an empty paragraph: a section break rides its own, never a sentence's.
    And every prose block of the tree is found in the text of the rendered pages,
    so a block dropped or hidden by any path through the renderer is caught by
    what came out rather than by what was intended.
    """
    fails = []
    rendered = 0
    for where, xlsx, extra in _PAGE_REPORTS:
        slope_data, solutions = (_solved() if xlsx is REINF_XLSX
                                 else _refused_solutions())
        opts = {"input_path": xlsx, "pd_figure": False,
                "lem_search_figure": False, "lem_solution_figure": False}
        opts.update(extra)
        with tempfile.TemporaryDirectory() as tmp:
            report, doc, path = _written(where, slope_data, solutions, opts, tmp)
            fails += _collapsed_with_content(doc, where)
            text = _rendered_text(path, opts["title"])
            if text is None:
                continue
            rendered += 1
            fails += _missing_from_the_page(report, text, where)
    if not rendered:
        print("Report: LibreOffice not installed — the rendered pages were not "
              "read back.")

    # Mutation. The defect this check is built on is the section break's
    # treatment reaching a paragraph of content: put it back, on the report whose
    # landscape section ends in a sentence, and the sentence has to be reported
    # missing.
    import xslope.report_docx as report_docx
    saved = report_docx._sect_break_carrier
    report_docx._sect_break_carrier = lambda doc: next(
        (p for p in reversed(doc.paragraphs) if p.text.strip()), None)
    try:
        slope_data, solutions = _refused_solutions()
        opts = {"input_path": NONCIRC_XLSX, "title": "Non-circular Surface",
                "method": "janbu", "pd_figure": False,
                "lem_search_figure": False, "lem_solution_figure": False}
        with tempfile.TemporaryDirectory() as tmp:
            _report, doc, _path = _written("the mutation", slope_data,
                                           solutions, opts, tmp)
        if not _collapsed_with_content(doc, "the mutation"):
            fails.append("a sentence was collapsed to a hidden mark and this "
                         "check passed the document")
    finally:
        report_docx._sect_break_carrier = saved
    return fails


def test_contents_page():
    """The contents page lists the report's own headings, from generation.

    The field is never marked dirty — that flag is what makes Word for Mac ask
    about fields on every open — so nothing updates it before a reader sees it,
    and what is cached in it has to be the contents themselves. They come from
    the same tree the document is written from, which is what makes a section
    that was toggled off absent from both.
    """
    import re
    fails = []
    from xslope.report import generate_report
    from xslope.report_docx import TOC_LEVELS

    slope_data, solutions = _solved()
    opts = {"input_path": REINF_XLSX, "title": "Sample Levee",
            "method": "spencer", "pd_figure": False,
            "lem_search_figure": False, "lem_solution_figure": False}
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "report.docx")
        ok, out = generate_report(slope_data, solutions, dict(opts), path)
        if not ok:
            return [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(path)
        doc = xml.get("word/document.xml", "")
        entries = _toc_result(doc)
        expected = [f"{number} {sec.title}"
                    for number, lvl, sec in out["report"].section_numbers()
                    if lvl <= TOC_LEVELS]

        if entries[:-1] != expected:
            fails.append(f"the contents list {entries[:-1]}, not the report's "
                         f"headings {expected}")
        hint = entries[-1] if entries else ""
        if "Update Field" not in hint or "F9" not in hint:
            fails.append(f"the last thing in the field result is {hint!r}, not "
                         f"the line about updating it")
        if hint in doc.split('w:fldCharType="end"')[-1]:
            fails.append("the update hint sits outside the field result; Word "
                         "would leave it behind when the table is updated")
        if "Page numbers" not in hint:
            fails.append("the hint does not say what an update adds")
        # Levels are levels: a sub-section is indented, in the template's style.
        styles = re.findall(r'<w:pStyle w:val="(TOC\d)"/>', doc)
        levels = [lvl for lvl, _t in out["report"].section_titles()
                  if lvl <= TOC_LEVELS]
        if styles != [f"TOC{lvl}" for lvl in levels]:
            fails.append(f"the contents are styled {styles}, not one TOC style "
                         f"per heading level {levels}")
        if "<w:t>1</w:t>" in doc[doc.index("TOC \\o"):doc.index("Traceability")]:
            fails.append("the cached contents carry a page number; the pagination "
                         "is Word's to compute")

    # A section switched off is in neither the document nor its contents.
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "no_slices.docx")
        ok, out = generate_report(slope_data, solutions,
                                  dict(opts, lem_slice_table=False), path)
        _names, xml = _docx_parts(path)
        entries = _toc_result(xml.get("word/document.xml", ""))
        if "Slice Table" in entries:
            fails.append("the contents list a Slice Table section the report was "
                         "told not to write")
    return fails


def _headings_in(doc_xml):
    """``[(style_id, text), ...]`` for every heading paragraph of a document."""
    import re
    out = []
    for para in re.findall(r"<w:p[ >].*?</w:p>", doc_xml, re.S):
        style = re.search(r'<w:pStyle w:val="(Heading\d)"/>', para)
        if style:
            text = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", para))
            out.append((style.group(1), text))
    return out


def test_heading_numbers():
    """The headings are numbered 1, 1.1, 1.1.1 by Word's own multilevel list.

    Numbering that is typed into the heading strings is numbering that goes wrong
    the first time a section is switched off. The document instead carries the
    numbering definition, bound to the heading styles, that Word writes when a
    user applies a multilevel list from the ribbon — so Word computes every
    number, the contents field picks them up, and a cross-reference to a section
    resolves to whatever its number has become.
    """
    import re
    fails = []
    from xslope.report import HEADING_LEVELS, generate_report
    from xslope.report_docx import HEADING_NUM_ID

    slope_data, solutions = _solved()
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "numbered.docx")
        ok, out = generate_report(
            slope_data, solutions,
            {"input_path": REINF_XLSX, "method": "spencer", "pd_figure": False,
             "lem_search_figure": False, "lem_solution_figure": False}, path)
        if not ok:
            return [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(path)
        doc = xml.get("word/document.xml", "")
        numbering = xml.get("word/numbering.xml", "")
        styles = xml.get("word/styles.xml", "")

        # The definition: one level per heading depth, each naming its style and
        # printing the whole path to it.
        if not numbering:
            return ["the package carries no word/numbering.xml, so nothing "
                    "numbers the headings"]
        block = re.search(rf'<w:abstractNum w:abstractNumId="{HEADING_NUM_ID}".*?'
                          rf'</w:abstractNum>', numbering, re.S)
        if block is None:
            return [f"the numbering part defines no list "
                    f"{HEADING_NUM_ID}: {numbering[:400]!r}"]
        text = block.group(0)
        want = [".".join(f"%{i + 1}" for i in range(lvl + 1))
                for lvl in range(HEADING_LEVELS)]
        got = re.findall(r'<w:lvlText w:val="([^"]*)"/>', text)
        if got != want:
            fails.append(f"the heading list prints {got}, not {want}")
        named = re.findall(r'<w:pStyle w:val="([^"]*)"/>', text)
        if named != [f"Heading{i + 1}" for i in range(HEADING_LEVELS)]:
            fails.append(f"the heading list's levels name {named}, not the "
                         f"heading styles")
        if f'<w:num w:numId="{HEADING_NUM_ID}"' not in numbering:
            fails.append(f"list {HEADING_NUM_ID} is defined but never made a "
                         f"numbering a style can point at")

        # And the binding: every heading style carries it, at its own level.
        for level in range(1, HEADING_LEVELS + 1):
            style = re.search(rf'<w:style [^>]*w:styleId="Heading{level}".*?'
                              rf'</w:style>', styles, re.S)
            if style is None:
                fails.append(f"the template defines no Heading{level} style")
                continue
            num_pr = re.search(r"<w:numPr>.*?</w:numPr>", style.group(0), re.S)
            if num_pr is None:
                fails.append(f"the Heading{level} style carries no numbering, so "
                             f"its headings print unnumbered")
                continue
            ids = re.findall(r'<w:numId w:val="(\d+)"/>', num_pr.group(0))
            ilvl = re.findall(r'<w:ilvl w:val="(\d+)"/>', num_pr.group(0))
            if ids != [str(HEADING_NUM_ID)]:
                fails.append(f"the Heading{level} style is numbered by list "
                             f"{ids}, not the heading list")
            if ilvl != [str(level - 1)]:
                fails.append(f"the Heading{level} style is bound to level {ilvl}, "
                             f"not level {level - 1} of the heading list")

        # Nothing is numbered twice: no heading string starts with a number.
        headings = _headings_in(doc)
        if not headings:
            fails.append("the document carries no heading paragraphs")
        for style_id, text in headings:
            if re.match(r"^\s*\d+(\.\d+)*[.\s)]", text):
                fails.append(f"the {style_id} {text!r} carries its own number; "
                             f"Word puts one in front of it as well")

        # The headings in the document are the report's own titles, and the
        # numbers the report caches for them are what this tree numbers to.
        report = out["report"]
        titles = [t for lvl, t in report.section_titles()
                  if lvl <= HEADING_LEVELS]
        if [t for _s, t in headings] != titles:
            fails.append(f"the document's headings are {[t for _s, t in headings]}, "
                         f"not the report's {titles}")
        numbers = report.section_numbers()
        top = [n for n, lvl, _s in numbers if lvl == 1]
        if top != [str(i + 1) for i in range(len(top))]:
            fails.append(f"the top-level sections number {top}")
        seen = {}
        for number, _lvl, _sec in numbers:
            seen[number] = seen.get(number, 0) + 1
            parent = number.rsplit(".", 1)[0]
            if "." in number and parent not in seen:
                fails.append(f"section {number} numbers under {parent}, which "
                             f"the report does not have")
        repeated = [n for n, count in seen.items() if count > 1]
        if repeated:
            fails.append(f"two sections are numbered the same: {repeated}")

        # The contents page carries them too, so a document nobody has pressed
        # F9 in already reads as a numbered report.
        entries = _toc_result(doc)
        if not entries or not entries[0].startswith("1 "):
            fails.append(f"the contents page opens on {entries[:1]!r}, which "
                         f"carries no heading number")
    return fails


def test_report_writes_one_file():
    """A report is a document, not a document and a folder.

    The figures are embedded in the .docx, so the PNGs are rendered into a
    temporary directory and go with it; a report generated with the figures
    switched off creates no directory at all.
    """
    fails = []
    from xslope.report import generate_report

    slope_data, solutions = _solved()
    opts = {"input_path": REINF_XLSX, "title": "Sample Levee",
            "method": "spencer"}
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "report.docx")
        ok, out = generate_report(slope_data, solutions, dict(opts), path)
        if not ok:
            return [f"generate_report failed: {out}"]
        left = sorted(os.listdir(tmp))
        if left != ["report.docx"]:
            fails.append(f"a default report left {left} behind, not the document "
                         f"alone")
        if not out["figures"]:
            fails.append("the result names none of the figures the document "
                         "carries")
        captions = [f.caption for f in out["report"].figures()]
        if list(out["figures"]) != captions:
            fails.append(f"the result reports {out['figures']}, not the "
                         f"captions {captions}")

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "nofig.docx")
        ok, _out = generate_report(
            slope_data, solutions,
            dict(opts, pd_figure=False, lem_search_figure=False,
                 lem_solution_figure=False), path)
        left = sorted(os.listdir(tmp))
        if not ok or left != ["nofig.docx"]:
            fails.append(f"a report with the figures off left {left} behind")

    # A caller that asks for the PNGs still gets them, in the directory it named.
    with tempfile.TemporaryDirectory() as tmp:
        figures = os.path.join(tmp, "figures")
        ok, _out = generate_report(slope_data, solutions, dict(opts),
                                   os.path.join(tmp, "kept.docx"),
                                   figure_dir=figures)
        kept = sorted(os.listdir(figures)) if os.path.isdir(figures) else []
        if not ok or not [f for f in kept if f.endswith(".png")]:
            fails.append(f"an explicit figure_dir kept {kept}")
    return fails


def test_docx_template():
    """The shipped template exists, and the document is built on it."""
    fails = []
    from xslope.report_docx import DEFAULT_TEMPLATE

    if not os.path.exists(DEFAULT_TEMPLATE):
        return [f"the default template is missing: {DEFAULT_TEMPLATE}"]
    _names, xml = _docx_parts(DEFAULT_TEMPLATE)
    styles = xml.get("word/styles.xml", "")
    for style in ("Heading1", "Title", "Caption", "TableGrid"):
        if style not in styles:
            fails.append(f"the template defines no {style} style")

    # The title page's look is the TEMPLATE's: ranged left, and closed by the rule
    # the Title style carries. A company template put in its place decides both,
    # which is only true while the renderer does not draw them itself.
    import re
    title_style = re.search(r'<w:style [^>]*w:styleId="Title".*?</w:style>',
                            styles, re.S)
    if title_style is None:
        fails.append("the template has no Title style to read")
    else:
        block = title_style.group(0)
        if "w:pBdr" not in block:
            fails.append("the Title style carries no rule under the title")
        if 'w:jc w:val="left"' not in block:
            fails.append("the Title style is not ranged left")
    body = xml.get("word/document.xml", "")
    if "<w:p " in body or "<w:p>" in body:
        fails.append("the template carries body content; it must ship empty")

    # The template is reproducible from its build script.
    import importlib.util
    script = os.path.join(_REPO, "tools", "build_report_template.py")
    if not os.path.exists(script):
        fails.append("tools/build_report_template.py is missing — the template "
                     "would be an unreproducible binary")
    else:
        spec = importlib.util.spec_from_file_location("build_report_template", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as tmp:
            rebuilt = mod.build(os.path.join(tmp, "t.docx"))
            _n, rebuilt_xml = _docx_parts(rebuilt)
            if rebuilt_xml.get("word/styles.xml") != styles:
                fails.append("rebuilding the template produces different styles "
                             "from the committed one")
    return fails


# --------------------------------------------------------------------------
# E. the column registry
# --------------------------------------------------------------------------

def test_column_registry():
    """Every report-worthy column is a column a solved model actually has."""
    fails = []
    from xslope import columns as cols

    _slope_data, solutions = _solved()
    slice_df = solutions["lem"][0]["slice_df"]
    have = set(slice_df.columns)
    for column in cols.solver_columns():
        if column.key not in have:
            fails.append(f"report column {column.key!r} is not in slice_df")
    for column in cols.report_columns():
        if not column.description.strip().endswith("."):
            fails.append(f"column {column.key!r} has no complete description")
        if not column.label.strip():
            fails.append(f"column {column.key!r} has no label")
    keys = [c.key for c in cols.SLICE_COLUMNS]
    if len(keys) != len(set(keys)):
        fails.append("the registry declares a column twice")

    # The units reach the headers, and only where they belong.
    from xslope.units import labels
    lbl = labels("imperial")
    headers, rows, legend = cols.slice_table(slice_df, lbl)
    if "W (lb/ft)" not in headers:
        fails.append(f"the weight column is not unit-labeled: {headers}")
    if "α (deg)" not in headers:
        fails.append(f"the base-angle column is not in degrees: {headers}")
    if "Slice" not in headers:
        fails.append(f"the slice number column is missing: {headers}")
    if any(h.startswith("Slice (") for h in headers):
        fails.append("the slice number carries a unit")
    if len(rows) != len(slice_df):
        fails.append(f"the table has {len(rows)} rows for {len(slice_df)} slices")
    if len(legend) != len(headers):
        fails.append(f"{len(legend)} legend entries for {len(headers)} columns")

    # An undeclared unit system leaves the headers bare, exactly as the plots do.
    bare, _rows, _legend = cols.slice_table(slice_df, None)
    if any("(lb/ft)" in h for h in bare):
        fails.append("an undeclared model got unit labels anyway")
    if "α (deg)" not in bare:
        fails.append("degrees are not a declared unit and must always be shown")

    # drop_empty removes a column of zeros, and never a geometry one.
    kept = {h.split(" (")[0] for h in headers}
    if "kW" in kept:
        fails.append("a seismic column was printed for a model with no seismic load")
    for needed in ("Δx", "α", "W", "c", "φ"):
        if needed not in kept:
            fails.append(f"drop_empty removed {needed}, which is never optional")
    return fails


# --------------------------------------------------------------------------
# I. the calculations
#
# The section's whole claim is that the printed numbers ARE the solution. So the
# checks do what a reviewer would: read the numbers off the page, divide them,
# and see whether the factor of safety comes back. Nothing here reads an
# internal float — every operand is parsed out of the string the document
# carries, and the reproduction is asserted at that printed precision.
# --------------------------------------------------------------------------

#: The methods a calculation is checked on. The three the ruling names, plus the
#: rest of the solver's methods, because each writes a different equation and a
#: silent failure on one of them is a section that quietly disappears.
CALC_METHODS = ("oms", "bishop", "spencer", "janbu", "corps", "lowe", "mprice")

#: The methods whose factor of safety IS a quotient of two sums, and which are
#: therefore checked by dividing the printed operands and getting F back.
#:
#: Spencer's is not one of them. Its F and θ are the pair at which two
#: equilibrium equations both vanish, reached by iterating the two together, and
#: the ratio of two force sums over that mass is a general horizontal balance
#: that any solution satisfies — printing it as "the equation" said something
#: true about the answer and nothing at all about the method. Spencer is verified
#: instead on its two equations balancing
#: (:func:`test_calculation_residuals`) and on reproducing each row's Q
#: (:func:`test_spencer_force_sums`).
QUOTIENT_METHODS = tuple(m for m in CALC_METHODS if m != "spencer")

#: How Spencer's section writes an equilibrium equation, transcribed from the
#: derivation — ``R_1 = sum{Q}``. The number the derivation gives it is in the
#: sentence above it, not on the equation.
SPENCER_EQUATION = r"^R_([12]) = (sum\{.+\})$"

#: And how it writes that same equation EVALUATED at the pair the iteration
#: reached: the equation restated, then the number it comes out at. Restated,
#: because the close is the section's demonstration that the equations balance —
#: a bare ``R_2 = 9.136e-05`` a page below the equation it belongs to is a digit
#: the reader has to go back and pair up for themselves.
SPENCER_EVALUATION = r"^R_([12]) = (sum\{.+\}) = (-?\d[^ ]*)$"

#: And the third thing the section prints: equation (1) or (2) reduced to the
#: terms this model carries. No page publishes it — the transcribed pair above it
#: is the derivation's, this is what is left of them here.
SPENCER_REDUCED_SUM = r"^F_[hv] = (0|[^=]+)$"

_CALC = {}


def _calc_report(method, options=None, xlsx=None):
    """A report of the sample model solved by ``method``, and its solved bundle.

    Cached per method: each of these is a full solve, and several checks read the
    same one. ``xlsx`` runs another model through the same path, for a check that
    needs a section built on a different shape of slope.
    """
    xlsx = xlsx or REINF_XLSX
    key = (method, tuple(sorted((options or {}).items())), xlsx)
    if key in _CALC:
        return _CALC[key]
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.report import build_report
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    slope_data = load_slope_data(xlsx)
    # A model carries circles or a non-circular surface, and the corpus models
    # the tolerance check reads carry both kinds.
    circles = slope_data.get("circles") or []
    surface_kw = ({"circle": circles[0]} if circles
                  else {"non_circ": slope_data.get("non_circ")})
    ok, out = generate_slices(slope_data, num_slices=15, **surface_kw)
    if not ok:
        raise RuntimeError(f"the sample model produced no slices: {out}")
    df, surface = out[0].copy(), out[1]
    with contextlib.redirect_stdout(io.StringIO()):
        results = solve_selected(method, df)
    if not isinstance(results, dict):
        _CALC[key] = (None, None)
        return _CALC[key]
    bundle = {"slice_df": df, "failure_surface": surface, "results": results,
              "search": None, "method": method}
    opts = {"method": method, "pd_figure": False, "lem_search_figure": False,
            "lem_solution_figure": False}
    opts.update(options or {})
    with tempfile.TemporaryDirectory() as tmp:
        report = build_report(slope_data, {"lem": [bundle]}, opts, tmp)
    _CALC[key] = (report, bundle)
    return _CALC[key]


def _calc_section(report):
    """The Calculations section of a report, or None."""
    for section in report.sections:
        for _lvl, node in section.walk():
            if node.title == "Calculations":
                return node
    return None


def _numbers(text):
    """Every number in a string, in order, as they are printed."""
    import re
    return [float(m) for m in re.findall(r"-?\d+\.?\d*(?:e[+-]?\d+)?", text)]


def _operands(section):
    """The printed operands of the factor of safety, read off the equations.

    Returns ``(numerator, denominator, quotient, corrected)`` as STRINGS, exactly
    as the document prints them — the corrected pair is Janbu's f_o line and is
    ``(None, None)`` for every other method.
    """
    import re
    num = den = quotient = corrected = factor = None
    for block in section.blocks:
        if block.kind != "math":
            continue
        text = block.notation
        frac = re.match(r"^F = frac\{([^{}]+)\}\{([^{}]+)\}$", text)
        if frac:
            num, den = frac.group(1), frac.group(2)
            continue
        plain = re.match(r"^F = ([\d.]+)$", text)
        if plain:
            quotient = plain.group(1)
            continue
        corr = re.match(r"^F_corr = f_o·F = ([\d.]+)·([\d.]+) = ([\d.]+)$", text)
        if corr:
            factor, quotient, corrected = corr.group(1), corr.group(2), corr.group(3)
    return num, den, quotient, factor, corrected


def _reproduces(num, den, quotient, factor, corrected, fs):
    """Does the factor of safety come back out of the printed operands?

    The tolerance is one unit in the last digit the factor of safety is printed
    to — the strictest statement that can be made about a rounded number.
    """
    from xslope.columns import format_fs

    tolerance = 10 ** -len(format_fs(1.0).split(".")[-1])
    if num is None or den is None or quotient is None:
        return False, "the equations do not print a quotient"
    computed = float(num) / float(den)
    if abs(computed - float(quotient)) > tolerance:
        return False, (f"{num}/{den} = {computed:.6f}, printed as {quotient}")
    if corrected is not None:
        product = float(factor) * float(quotient)
        if abs(product - float(corrected)) > tolerance:
            return False, (f"{factor}x{quotient} = {product:.6f}, printed as "
                           f"{corrected}")
        if abs(product - fs) > tolerance:
            return False, (f"the printed operands give {product:.6f}, the solver "
                           f"{fs:.6f}")
        return True, ""
    if abs(computed - fs) > tolerance:
        return False, (f"the printed operands give {computed:.6f}, the solver "
                       f"{fs:.6f}")
    return True, ""


def test_force_term_registry():
    """Every force is in every equation, or the registry says why it is not.

    :data:`xslope.report.FORCE_TERMS` is the one declaration of the forces a
    slice can carry, and the equations the section prints — the two moment sums,
    the horizontal balance, Spencer's two force sums and the base-normal
    equation — are assembled from it. A force carried into some of them and
    forgotten in the rest is what prints a term in one equation and denies it in
    the next, so every entry has to hold either a contribution or a stated reason
    for each of them.

    And the contributions are all evaluated on a solved model: a term that names
    an array or a column nothing writes would never appear, and would never be
    missed, because a term with no values is a term that reads as absent. That
    includes the terms a stated reason carries as PUBLISHED — the shear on the
    top of the slice, which Spencer's section transcribes and no model exercises.
    """
    fails = []
    import numpy as np
    from xslope.report import (CONSUMERS, EQUATION_SYMBOLS, FORCE_TERMS,
                               ForceTerm, NotApplicable, PASSIVE_COLUMNS,
                               SECTION_SYMBOL_GROUPS, SYMBOL_GROUPS, Term,
                               _Calc, _calc_arrays)

    _slope_data, solutions = _solved()
    df = solutions["lem"][0]["slice_df"]
    A = _calc_arrays(df)
    C = _Calc(df, A, right_facing=False)
    have = set(df.columns)

    keys = [t.key for t in FORCE_TERMS]
    if len(set(keys)) != len(keys):
        fails.append(f"two entries share a key: {keys}")

    for term in FORCE_TERMS:
        for column in term.columns + tuple(c for _n, c, _d in term.arrays):
            if column not in have:
                fails.append(f"{term.key}: names column {column!r}, which a "
                             f"solved slice table does not carry")
        for consumer in CONSUMERS:
            got = getattr(term, consumer)
            if isinstance(got, NotApplicable):
                if not got.reason.strip():
                    fails.append(f"{term.key}: {consumer} is not applicable for "
                                 f"no stated reason")
                # A published term of a force the solver does not carry is held
                # against the same rules as any other, and has to be zero on
                # every slice: it is printed in the transcription and summed
                # nowhere, and a nonzero one is a term gone missing from the
                # arithmetic.
                for contribution in got.published:
                    values = np.asarray(contribution.values(C), dtype=float)
                    if values.shape != (len(df),):
                        fails.append(f"{term.key}: {consumer} publishes "
                                     f"{contribution.symbol!r}, which gives "
                                     f"{values.shape}, not one value per slice")
                    elif values.any():
                        fails.append(f"{term.key}: {consumer} publishes "
                                     f"{contribution.symbol!r} as unexercised, "
                                     f"and it is not zero on every slice")
                continue
            if not got:
                fails.append(f"{term.key}: {consumer} carries no term and no "
                             f"reason for carrying none")
                continue
            for contribution in got:
                if not isinstance(contribution, Term):
                    fails.append(f"{term.key}: {consumer} holds "
                                 f"{contribution!r}, which is not a term")
                    continue
                if contribution.sign not in (+1, -1):
                    fails.append(f"{term.key}: {consumer} term "
                                 f"{contribution.symbol!r} has sign "
                                 f"{contribution.sign}")
                try:
                    values = np.asarray(contribution.values(C), dtype=float)
                except Exception as exc:
                    fails.append(f"{term.key}: {consumer} term "
                                 f"{contribution.symbol!r} does not evaluate: "
                                 f"{exc!r}")
                    continue
                if values.shape != (len(df),):
                    fails.append(f"{term.key}: {consumer} term "
                                 f"{contribution.symbol!r} gives "
                                 f"{values.shape}, not one value per slice")
        for symbol in term.symbols:
            if symbol.group in SECTION_SYMBOL_GROUPS:
                # A letter one derivation gives a force the others give another.
                # It is defined for its own section and must NOT be the report's
                # definition, which is the whole reason it is held apart: a
                # section symbol that agrees with the nomenclature is a row the
                # nomenclature already had.
                if EQUATION_SYMBOLS.get(symbol.name) == symbol.meaning:
                    fails.append(f"{term.key}: symbol {symbol.name!r} is held "
                                 f"apart for one section and is what the "
                                 f"nomenclature already defines")
                continue
            if symbol.group not in SYMBOL_GROUPS:
                fails.append(f"{term.key}: symbol {symbol.name!r} is in group "
                             f"{symbol.group!r}, which the nomenclature has no "
                             f"place for")
            if EQUATION_SYMBOLS.get(symbol.name) != symbol.meaning:
                fails.append(f"{term.key}: symbol {symbol.name!r} is not the "
                             f"one the nomenclature defines")

    passive = tuple(c for t in FORCE_TERMS if t.passive for c in t.columns)
    if passive != PASSIVE_COLUMNS:
        fails.append(f"the passive gate is {PASSIVE_COLUMNS}, and the registry's "
                     f"passive entries carry {passive}")
    if not passive:
        fails.append("no entry is marked passive, so the gate tests nothing")

    # Mutation: a force added to some equations and not the rest has to be
    # impossible to write. Every consumer is a required field, so the half-added
    # entry does not construct.
    whole = {"key": "x", "columns": (), "arrays": (), "symbols": (),
             "feature": "", "passive": False}
    whole.update({c: NotApplicable("nothing") for c in CONSUMERS})
    for consumer in CONSUMERS:
        half = {k: v for k, v in whole.items() if k != consumer}
        try:
            ForceTerm(**half)
        except TypeError:
            continue
        fails.append(f"a force declared for every equation but {consumer} was "
                     f"accepted; nothing stops a term being added half-way")
    return fails


#: Models whose converged solution the report used to refuse, with the method
#: that solved them: the quotient reproduced the factor of safety to a few parts
#: in a million, or Spencer's two imbalances vanished to the tolerance the Newton
#: pair was driven to, and a relative 1e-6 turned each of them into a method
#: section with no working in it. The last two are the collapsed-scale case,
#: where Q acts through the coordinate origin on every slice and the moment terms
#: sum to two parts in a billion of one force-length unit.
_CONVERGED_BUT_REFUSED = (
    ("geostudio", "gs2_46", "spencer"),
    ("rocscience", "vp037", "spencer"),
    ("rocscience", "vp040", "janbu"),
    ("rocscience", "vp061a", "bishop"),
    ("rocscience", "vp061a", "janbu"),
    ("geostudio", "gs2_26", "spencer"),
    ("rocscience", "vp043", "spencer"),
)


def test_calculation_tolerance_follows_the_solver():
    """A solution that converged as far as it was asked to gets its working.

    The report evaluates its equation and compares the answer with the solver's
    before printing anything, which is what keeps a wrong calculation off the
    page. The comparison has to be made against what the solver delivers:
    Bishop's and Janbu's iterations stop at a step in F of 1e-6 and Spencer's
    Newton pair at imbalances of 1e-4, and a gate of a relative 1e-6 demanded
    more than either and refused seven converged model-method pairs outright.

    Every one of them is required to print here, and the gate is required to
    still refuse an answer that is wrong rather than rounded.
    """
    fails = []
    from xslope.report import (CALC_SAFETY_FACTOR, CALC_TOLERANCE, _closes,
                               _solver_tolerance)

    for vendor, model, method in _CONVERGED_BUT_REFUSED:
        xlsx = os.path.join(_REPO, "docs", "verification", "files", vendor,
                            f"{model}.xlsx")
        report, _bundle = _calc_report(method, xlsx=xlsx)
        if report is None:
            fails.append(f"{model} did not solve with {method}")
            continue
        if _calc_section(report) is None:
            fails.append(f"{model} under {method} converged and still prints no "
                         f"calculation")

    # The tolerances are read off the solvers themselves, so retuning a solver
    # moves the gate that judges its answers with it.
    for method, wanted in (("bishop", 1e-6), ("janbu", 1e-6), ("spencer", 1e-4),
                           ("corps", 1e-6), ("lowe", 1e-6), ("mprice", 1e-6)):
        got = _solver_tolerance(method)
        if got != wanted:
            fails.append(f"{method} converges to {wanted:.0e} and the gate "
                         f"reads {got:.0e}")
    if _solver_tolerance("oms") != 0.0:
        fails.append("the Ordinary Method of Slices is closed form; its gate "
                     "should allow no iteration tolerance at all")

    # Mutation, the quotient: a factor of safety out by one percent is not a
    # converged solution rounding, and no allowance may let it through.
    for method, scale in (("oms", 2.0), ("bishop", 2.0), ("janbu", 2.0),
                          ("mprice", 2.0), ("spencer", 144.0)):
        if _closes(0.01 * scale, scale, method):
            fails.append(f"{method}: a residual one percent of {scale} was "
                         f"accepted as a converged solution")
    # Mutation, the residual: an imbalance a thousand times what the solver
    # itself stops at, against a scale at which the relative test is no looser
    # than the absolute one, so neither statement can excuse it.
    for method in ("spencer", "bishop"):
        allowance = CALC_SAFETY_FACTOR * _solver_tolerance(method)
        if _closes(1000 * allowance, allowance / CALC_TOLERANCE, method):
            fails.append(f"{method}: an imbalance a thousand times the "
                         f"allowance was accepted")
    # And the allowance is a real widening: what the solver delivers has to be
    # more than the relative test alone would take.
    if not _closes(15 * _solver_tolerance("bishop"), 1.94, "bishop"):
        fails.append("a Bishop solution fifteen tolerances from its own fixed "
                     "point is still refused")
    return fails


def _method_block(report, method):
    """The whole of one method's block — its section and every descendant."""
    from xslope.report import method_label

    label = method_label(method)
    for section in report.sections:
        for _lvl, node in section.walk():
            if node.title == label:
                return node
    return None


def _subtree_prose(section):
    """Every paragraph in a section and its children, as one list."""
    return [b.text for _lvl, node in section.walk()
            for b in node.blocks if b.kind == "prose"]


def test_a_method_block_never_goes_quiet():
    """A method whose equilibrium cannot be worked through says what is true of
    it instead.

    The passive-support model is the case: its capacity mobilizes with the soil,
    so it enters divided by the factor of safety and stands on both sides of the
    balance the force methods solve. There is no quotient behind that factor of
    safety, and for five methods the block simply stopped after the slice table
    — a factor of safety with nothing after it, which reads as an omission.

    The moment methods CAN show passive support: it makes a resisting moment of
    its own. So the same model is required to carry the working under Bishop and
    the Ordinary Method of Slices, and the sentence under neither.
    """
    fails = []
    import xslope.report as report

    for method in ("janbu", "corps", "lowe", "mprice", "spencer"):
        built, _bundle = _calc_report(method, xlsx=PASSIVE_XLSX)
        if built is None:
            fails.append(f"the passive model did not solve with {method}")
            continue
        block = _method_block(built, method)
        if block is None:
            fails.append(f"{method}: the report has no block for the method")
            continue
        if any(node.title == "Calculations" for _lvl, node in block.walk()):
            fails.append(f"{method}: a passive model has no quotient to work "
                         f"through and the block carries a Calculations section")
        said = [p for p in _subtree_prose(block) if p == report.PASSIVE_NOTE]
        if len(said) != 1:
            fails.append(f"{method}: the block says why no working is printed "
                         f"{len(said)} times, not once")

    for method in ("bishop", "oms"):
        built, _bundle = _calc_report(method, xlsx=PASSIVE_XLSX)
        if built is None:
            fails.append(f"the passive model did not solve with {method}")
            continue
        block = _method_block(built, method)
        if not any(node.title == "Calculations" for _lvl, node in block.walk()):
            fails.append(f"{method}: passive support makes a resisting moment "
                         f"and the block prints no calculation")
        if report.PASSIVE_NOTE in _subtree_prose(block):
            fails.append(f"{method}: the block prints its working and says it "
                         f"cannot")

    # A tolerance refusal states the two numbers. The gate is widened to the
    # solver's own, so this is exercised on a residual constructed to fail it
    # rather than on a model: what is required is that the sentence carries the
    # quotient and the solution it does not return.
    denied = report._calculation(None, {"slice_df": None, "results": {}}, "bishop")
    if denied != (None, ""):
        fails.append(f"a model with no slices produced {denied!r}")

    # Mutation: with nothing to say, the block goes quiet again, which is what
    # this check exists to catch.
    saved = report.PASSIVE_NOTE
    report.PASSIVE_NOTE = ""
    try:
        quiet, _bundle = _calc_report("janbu", options={"title": "quiet"},
                                      xlsx=PASSIVE_XLSX)
        block = _method_block(quiet, "janbu") if quiet is not None else None
        if block is not None and saved in _subtree_prose(block):
            fails.append("the passive sentence is printed from somewhere other "
                         "than the refusal that withheld the calculation, so "
                         "the two can disagree")
    finally:
        report.PASSIVE_NOTE = saved
    return fails


def _quoted_numbers(sentence, fs, bounds):
    """The numbers a refusal sentence prints that are not factors of safety.

    The reported factor of safety is one whatever it is; every other number in
    the sentence has to lie in the range one can take.
    """
    low, high = bounds
    return [n for n in _numbers(sentence)
            if abs(n - fs) > 5e-7 and not low <= n <= high]


def test_a_refusal_prints_no_number_it_cannot_stand_behind():
    """A refused calculation states the mismatch; it does not typeset a
    degenerate quotient as a measurement.

    The equation is evaluated on the values the bundle carries, and a bundle
    whose slices were never solved carries initial ones: the driving sum is near
    zero and the quotient comes out at 1e14. Printed to six decimals beside the
    factor of safety the solution reports, that reads as a computed result and is
    an artifact of arithmetic on an unsolved frame. A mismatch between two
    factors of safety still states both.
    """
    fails = []
    import xslope.report as report
    from xslope.report import CREDIBLE_FS, _mismatch_note

    fs = 1.759948
    near = _mismatch_note(1.752341, fs)
    for want in ("1.752341", f"{fs:.6f}"):
        if want not in near:
            fails.append(f"a mismatch between two factors of safety does not "
                         f"state {want}: {near!r}")

    for computed in (2.7206093492393162e14, 1e-9, -3.4, CREDIBLE_FS[1] * 1.5):
        said = _mismatch_note(computed, fs)
        loose = _quoted_numbers(said, fs, CREDIBLE_FS)
        if loose:
            fails.append(f"the equation evaluated at {computed:g} is reported "
                         f"as {loose}, which no slope has: {said!r}")
        if f"{fs:.6f}" not in said:
            fails.append(f"the refusal at {computed:g} does not state the "
                         f"factor of safety the solution reports: {said!r}")

    # And on the model it was measured on: a non-circular slope whose slices the
    # report is handed unsolved.
    slope_data, solutions = _refused_solutions()
    from xslope.report import build_report
    with tempfile.TemporaryDirectory() as tmp:
        with contextlib.redirect_stdout(io.StringIO()):
            built = build_report(
                slope_data, solutions,
                {"input_path": NONCIRC_XLSX, "method": "janbu",
                 "pd_figure": False, "lem_search_figure": False,
                 "lem_solution_figure": False}, tmp)
    refusals = [p for p in _prose(built)
                if "does not return the solution" in p]
    if len(refusals) != 1:
        fails.append(f"the unsolved frame produced {len(refusals)} refusals, "
                     f"not one")
    fs_reported = solutions["lem"][0]["results"]["FS"]
    for said in refusals:
        loose = _quoted_numbers(said, fs_reported, CREDIBLE_FS)
        if loose:
            fails.append(f"the refusal prints {loose}: {said!r}")

    # Mutation: with the range opened, the degenerate quotient is typeset again,
    # which is what this check exists to catch.
    saved = report.CREDIBLE_FS
    report.CREDIBLE_FS = (0.0, float("inf"))
    try:
        said = report._mismatch_note(2.7206093492393162e14, fs)
        if not _quoted_numbers(said, fs, saved):
            fails.append("the raw quotient was printed and this check passed "
                         "the sentence")
    finally:
        report.CREDIBLE_FS = saved
    return fails


def test_calculation_reproduces_fs():
    """The factor of safety is re-derived from the operands as PRINTED, for
    every method, and matches the solver to the last digit it prints.

    This is the check that keeps the section from drifting: if the equation ever
    stops being the one the solver evaluates, or a sum is printed to too few
    digits to divide, the number will not come back.
    """
    fails = []
    for method in QUOTIENT_METHODS:
        report, bundle = _calc_report(method)
        if report is None:
            fails.append(f"{method}: the sample model did not solve")
            continue
        section = _calc_section(report)
        if section is None:
            fails.append(f"{method}: the report carries no Calculations section")
            continue
        fs = float(bundle["results"]["FS"])
        num, den, quotient, factor, corrected = _operands(section)
        ok, why = _reproduces(num, den, quotient, factor, corrected, fs)
        if not ok:
            fails.append(f"{method}: {why}")

    # And Spencer prints no quotient at all. A ratio of two force sums is true of
    # its answer and true of everybody else's, which is exactly why printing it
    # as Spencer's equation was wrong: it is not the method.
    report, _bundle = _calc_report("spencer")
    section = _calc_section(report) if report is not None else None
    if section is None:
        fails.append("spencer: the report carries no Calculations section")
    else:
        num, den, quotient, _factor, _corrected = _operands(section)
        if num is not None or den is not None:
            fails.append(f"spencer prints a quotient of two sums, "
                         f"{num}/{den}, as though it were the method")
        if quotient is not None:
            fails.append(f"spencer closes on the arithmetic F = {quotient}; its "
                         f"solution is the root of two equilibrium equations, "
                         f"not the value of an expression")

    # Mutation: the check has to be able to fail. A sum wrong in its third
    # significant digit moves the factor of safety in its third decimal, which is
    # the smallest error the printed precision can be asked to catch.
    report, bundle = _calc_report("bishop")
    section = _calc_section(report)
    num, den, quotient, factor, corrected = _operands(section)
    fs = float(bundle["results"]["FS"])
    for scale in (1.01, 0.999):
        nudged = f"{float(num) * scale:.6g}"
        ok, _why = _reproduces(nudged, den, quotient, factor, corrected, fs)
        if ok:
            fails.append(f"a numerator moved from {num} to {nudged} still "
                         f"reproduced the factor of safety; the check cannot "
                         f"fail")
    return fails


def test_calculation_sums_are_printed_precisely_enough():
    """The sums carry enough digits to divide, and that number is derived rather
    than picked: :data:`xslope.columns.SUM_DIGITS` follows from the format the
    factor of safety is printed in."""
    fails = []
    from xslope.columns import FS_FMT, SUM_DIGITS, format_sum

    decimals = len(FS_FMT.format(1.0).split(".")[-1])
    # The quotient's relative error is at most 10^(1-s); the factor of safety
    # needs it under 10^-decimals / F. Six digits must clear that for F = 100.
    if 10 ** (1 - SUM_DIGITS) > 10 ** -decimals / 100:
        fails.append(f"{SUM_DIGITS} significant digits cannot reproduce a "
                     f"factor of safety printed to {decimals} decimals")
    for value, want in ((1234567.89, 7), (0.000123456789, 6), (1.5, 2)):
        text = format_sum(value)
        got = len(text.replace("-", "").replace(".", "").lstrip("0"))
        if got < min(want, SUM_DIGITS):
            fails.append(f"format_sum({value}) = {text!r} keeps {got} digits")
        if abs(float(text) - value) > abs(value) * 10 ** (1 - SUM_DIGITS):
            fails.append(f"format_sum({value}) = {text!r} loses too much")
    return fails


def test_calculation_terms_follow_the_model():
    """The equation carries the terms the model exercises and no others.

    Two models of the same shape: one with reinforcement crossing the failure
    surface and one without. The reinforcement terms are in the first equation
    and absent from the second — the section describes THIS analysis, not the
    general case.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.report import build_report
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    def equation(xlsx, method="bishop"):
        slope_data = load_slope_data(xlsx)
        ok, out = generate_slices(slope_data, circle=slope_data["circles"][0],
                                  num_slices=15)
        if not ok:
            return None, None
        df, surface = out[0].copy(), out[1]
        with contextlib.redirect_stdout(io.StringIO()):
            results = solve_selected(method, df)
        if not isinstance(results, dict):
            return None, None
        with tempfile.TemporaryDirectory() as tmp:
            report = build_report(
                slope_data,
                {"lem": [{"slice_df": df, "failure_surface": surface,
                          "results": results, "search": None, "method": method}]},
                {"method": method, "pd_figure": False,
                 "lem_search_figure": False, "lem_solution_figure": False}, tmp)
        section = _calc_section(report)
        if section is None:
            return None, report
        maths = [b.notation for b in section.blocks if b.kind == "math"]
        prose = " ".join(b.text for b in section.blocks if b.kind == "prose")
        return (maths, prose), report

    nailed = os.path.join(_REPO, "docs", "inputs", "slope",
                          "xslope_nail_axial.xlsx")
    got, _report = equation(nailed)
    if got is None:
        fails.append("the nailed sample produced no calculation to compare")
    else:
        maths, prose = got
        equation_line = next((m for m in maths if m.startswith("F = frac{sum")), "")
        if "P" not in equation_line:
            fails.append(f"a model with reinforcement crossing the surface "
                         f"prints no P term: {equation_line}")
        # The absent list is the CLAUSE that follows, and a semicolon ends a
        # clause. The same sentence goes on to say what the model does carry —
        # passive support, which mobilizes with the soil and joins the numerator
        # — and reading a fixed 200 characters past "carries no" read that as
        # part of the list, calling a reinforced model unreinforced on the
        # strength of the sentence that says it is reinforced.
        absent = prose.split("The model carries no")[-1].split(";")[0][:200]
        if "reinforcement" in absent:
            fails.append("a reinforced model is described as carrying no "
                         "reinforcement")
        # And the reading still catches the thing it is for.
        named = ("The model carries no reinforcement crossing the failure "
                 "surface; the piles in this model mobilize with the soil")
        if "reinforcement" not in \
                named.split("The model carries no")[-1].split(";")[0][:200]:
            fails.append("the absent-force clause is read in a way that would "
                         "miss a reinforced model called unreinforced")

    got, _report = equation(REINF_XLSX)
    if got is None:
        fails.append("the sample model produced no calculation")
    else:
        maths, prose = got
        equation_line = next((m for m in maths if m.startswith("F = frac{sum")), "")
        for absent in ("P ", "P·", "kW", "H cos", "T·a_t"):
            if absent in equation_line:
                fails.append(f"the equation prints a {absent!r} term for a model "
                             f"with none: {equation_line}")
        if "carries no" not in prose:
            fails.append("the section does not say what was left out")
        for named in ("seismic load", "reinforcement", "pile force"):
            if named not in prose:
                fails.append(f"the omission sentence does not name {named!r}")
    return fails


def test_calculation_columns():
    """The per-slice terms of the sums are columns of the slice table.

    That is what lets the section be four lines instead of a walk through every
    slice, so the columns have to be really there — in the table the report
    prints, under the labels the section names.
    """
    fails = []
    from xslope import columns as cols

    # Spencer's F_h and F_v are in the list for the same reason as Q_s: the
    # preamble prints the equation that turns them into Q, and a reader can only
    # check a row against it if the row carries them. F_R and F_D are NOT: they
    # are the two halves of a quotient, and Spencer's section prints none.
    for method, wanted in (("bishop", ("M_R", "M_D")),
                           ("spencer", ("F_h", "F_v", "Q_s", "y_Q"))):
        report, _bundle = _calc_report(method)
        if report is None:
            fails.append(f"{method}: the sample model did not solve")
            continue
        table = next((t for t in report.tables() if t.landscape), None)
        if table is None:
            fails.append(f"{method}: there is no slice table")
            continue
        labels = [h.split(" (")[0] for h in table.headers]
        for label in wanted:
            if label not in labels:
                fails.append(f"{method}: the slice table has no {label} column "
                             f"({labels})")
        legend = dict(table.legend)
        for head in table.headers:
            if head not in legend:
                fails.append(f"{method}: {head} has no legend entry")
        # And the section points at every one of them by the label the table
        # prints — a column nothing in the prose sends the reader to is a column
        # they have no reason to read.
        section = _calc_section(report)
        prose = " ".join(b.text for b in section.blocks if b.kind == "prose")
        for label in wanted:
            if label not in prose:
                fails.append(f"{method}: the calculation never names column "
                             f"{label}")
        if method == "spencer":
            for unwanted in ("F_R", "F_D"):
                if unwanted in labels:
                    fails.append(f"spencer's slice table carries a {unwanted} "
                                 f"column; it exists to be divided, and there "
                                 f"is no quotient in the section to divide it")

    # The computed columns are declared as the report's own, so a solved
    # slice_df is not expected to carry them.
    _slope_data, solutions = _solved()
    have = set(solutions["lem"][0]["slice_df"].columns)
    for column in cols.computed_columns():
        if column.key in have:
            fails.append(f"{column.key!r} is declared as computed by the report "
                         f"but the solver already writes it")
    if not cols.computed_columns():
        fails.append("no column is declared as the report's own; the split is "
                     "not being tested")
    return fails


def test_calculation_residuals():
    """Spencer's two equilibrium equations, rebuilt from the values as PRINTED.

    This is Spencer's answer to the check the other methods get from dividing
    two printed sums: its solution is the pair (F, θ) at which R_1 and R_2 both
    vanish, so what a reviewer verifies is that they DO vanish — at the printed
    θ, from the printed columns, against the size of the sums they cancel within.
    Both are re-formed here the way the section says a reader can re-form them:
    R_1 is the Q_s column added up, R_2 is that column weighted by x_c and y_Q.

    The solver's own residuals are vanishing. A reader adding up the printed
    column gets something larger — the rounding of every row — and the section
    says how large that is allowed to be. This checks both statements.
    """
    import math
    import re
    fails = []
    from xslope.report import CALC_TOLERANCE, PRINTED_RESIDUAL_TOLERANCE

    report, bundle = _calc_report("spencer")
    if report is None:
        return ["the sample model did not solve with Spencer's method"]
    section = _calc_section(report)
    table = next((t for t in report.tables() if t.landscape), None)
    if section is None or table is None:
        return ["there is no calculation or no slice table to read"]

    labels = [h.split(" (")[0] for h in table.headers]
    needed = ("Q_s", "x_c", "y_Q")
    if not set(needed) <= set(labels):
        return [f"the slice table carries no {sorted(set(needed) - set(labels))} "
                f"column, so neither equilibrium sum can be re-formed: {labels}"]
    at = {name: labels.index(name) for name in needed}
    values = [float(row[at["Q_s"]]) for row in table.rows
              if row[at["Q_s"]].strip()]
    if len(values) != len(table.rows):
        fails.append(f"{len(values)} of {len(table.rows)} Q_s cells hold a number")
    total = sum(values)
    scale = sum(abs(v) for v in values) or 1.0
    if abs(total) > PRINTED_RESIDUAL_TOLERANCE * scale:
        fails.append(f"the printed Q_s column sums to {total:.4g}, which is "
                     f"{abs(total) / scale:.1e} of the {scale:.6g} its "
                     f"magnitudes come to — past the stated "
                     f"{PRINTED_RESIDUAL_TOLERANCE:.0e}")

    prose = " ".join(b.text for b in section.blocks if b.kind == "prose")
    if "thousandth" not in prose:
        fails.append("the section does not state, in words, how closely the "
                     "printed values close")

    # The converged F is the answer a reader is looking for, so the sentence
    # that states it sets it in bold — findable without reading the sentence.
    from xslope.columns import format_fs
    printed_fs = format_fs(bundle["results"]["FS"])
    said = [b for b in section.blocks
            if b.kind == "prose" and f"F = {printed_fs}" in b.text]
    if not said:
        fails.append(f"the section never states the converged F = {printed_fs}")
    elif not any(printed_fs in b.bold for b in said):
        fails.append("the converged factor of safety is not set in bold where "
                     "the section states it")

    # θ comes off the page too: the section states it, and without it neither
    # the moment sum nor a single row of Q can be reproduced by a reader.
    stated = re.search(r"θ = (-?\d+\.\d+) degrees", prose)
    if stated is None:
        fails.append("the section never prints the interslice inclination it "
                     "converged at, which the moment sum cannot be re-formed "
                     "without")
    else:
        theta = math.radians(float(stated.group(1)))
        moments = [float(row[at["Q_s"]]) * (float(row[at["x_c"]]) * math.sin(theta)
                                            - float(row[at["y_Q"]]) * math.cos(theta))
                   for row in table.rows]
        m_total, m_scale = sum(moments), sum(abs(v) for v in moments) or 1.0
        if abs(m_total) > PRINTED_RESIDUAL_TOLERANCE * m_scale:
            fails.append(f"the moment sum re-formed from the printed columns is "
                         f"{m_total:.4g}, which is {abs(m_total) / m_scale:.1e} "
                         f"of the {m_scale:.6g} its terms come to — past the "
                         f"stated {PRINTED_RESIDUAL_TOLERANCE:.0e}")

    # The two equilibrium equations, and each of them evaluated at the pair the
    # iteration reached. Their numbers are in the sentence above them — no
    # equation on the page carries one — and that sentence is what
    # :func:`test_equation_numbers_are_in_the_prose` holds to the derivation.
    maths = [b.notation for b in section.blocks if b.kind == "math"]
    transcribed = {}
    for notation in maths:
        found = re.match(SPENCER_EQUATION, notation)
        if found:
            transcribed[found.group(1)] = found.group(2)
    if sorted(transcribed) != ["1", "2"]:
        fails.append(f"the section does not print the two equilibrium equations "
                     f"of the derivation: {maths}")
    if "equations (27) and (28)" not in " ".join(
            b.text for b in section.blocks if b.kind == "prose"):
        fails.append("the section prints the two equilibrium equations without "
                     "saying which equations of the derivation they are")
    evaluated = [re.match(SPENCER_EVALUATION, n) for n in maths
                 if re.match(SPENCER_EVALUATION, n)]
    if len(evaluated) != 2:
        fails.append(f"the section evaluates {len(evaluated)} of its two "
                     f"equilibrium equations at the converged pair: {maths}")
    printed = {}
    for found in evaluated:
        which, notation, value = found.groups()
        name = f"R_{which}"
        # The evaluated line has to BE the equation it is evaluating. An
        # evaluation whose left side drifted from the transcription is a number
        # a reader cannot tie to anything.
        stated_eq = transcribed.get(which, "")
        if notation != stated_eq:
            fails.append(f"{name} is evaluated as {notation!r} but transcribed "
                         f"as {stated_eq!r}; the two do not read as one equation")
        if "e-" not in value and "e+" not in value:
            fails.append(f"a residual is not in scientific notation: {value}")
        printed[name] = float(value)
    # And they really are zero, to the tolerance the section is gated on — the
    # same 1e-6 the other methods' quotients must reproduce F to.
    for name, against in (("R_1", scale), ("R_2", m_scale if stated else None)):
        if name not in printed or against is None:
            continue
        if abs(printed[name]) > CALC_TOLERANCE * against:
            fails.append(f"{name} = {printed[name]:.4g} against a scale of "
                         f"{against:.6g} is {abs(printed[name]) / against:.1e} — "
                         f"past the {CALC_TOLERANCE:.0e} a converged solution "
                         f"has to close to")

    # Both equations close AGAINST SOMETHING, and the section says what. A
    # residual printed on its own is a digit a reader has no way to judge, so
    # each is required to be stated against the sum of the magnitudes it cancels
    # within — the solver's own totals, in the report's own format — and the
    # closure each is claimed to reach is required to be true and not overstated.
    fails += _spencer_closures_are_true(bundle, prose, printed)

    # Morgenstern-Price prints its two residuals as well.
    report, _bundle = _calc_report("mprice")
    section = _calc_section(report) if report is not None else None
    if section is None:
        fails.append("Morgenstern-Price produced no calculation section")
    else:
        maths = [b.notation for b in section.blocks if b.kind == "math"]
        if not any(m.startswith("Z_n = ") for m in maths):
            fails.append(f"Morgenstern-Price prints no force residual: {maths}")
        if not any(m.startswith("sum{M_O} = ") for m in maths):
            fails.append(f"Morgenstern-Price prints no moment residual: {maths}")
    return fails


#: The counts and magnitudes a closure is spoken in — "one part in eight
#: billion". Read back off the page here and turned into a number, so the phrase
#: is a claim the check can measure and not decoration.
_CLOSURE_COUNTS = {"a": 1, "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
                   "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10}
_CLOSURE_MAGNITUDES = {"thousand": 1e3, "million": 1e6,
                       "billion": 1e9, "trillion": 1e12}


def _spencer_closures_are_true(bundle, prose, printed):
    """Each equilibrium equation is stated against the size of the sums it
    cancels within, and the closure claimed for it is the closure it has.

    Zero is not a number a residual can be judged against — 9.136e-05 is tiny in
    a moment sum of ninety thousand and enormous in one of a tenth. So the
    section states both totals and what fraction of each the residual is, and
    both statements are checked here: the totals against the solver's own, the
    claimed closure against the arithmetic. The phrase carries one figure, so it
    is allowed the rounding of one figure either way; past that — a closure
    overstated by more than that rounding, or understated by a factor of ten —
    the sentence has stopped describing this solution.
    """
    import re
    fails = []
    from xslope.columns import format_sum
    from xslope.report import _spencer_state

    state = _spencer_state(bundle["slice_df"])
    if state is None:
        return ["the Spencer state the section is written from cannot be rebuilt"]
    for name, key, what in (("R_1", "scale", "the magnitudes of Q"),
                            ("R_2", "m_scale", "its moment terms")):
        total = format_sum(state[key])
        if total not in prose:
            fails.append(f"{name} is printed against nothing: the section never "
                         f"states the {total} that {what} come to, so a reader "
                         f"has no size to judge the residual by")

    claims = re.findall(r"one part in ([a-z]+|[\d,]+) "
                        r"(thousand|million|billion|trillion)", prose)
    if len(claims) != 2:
        return fails + [f"the section states {len(claims)} closures for its two "
                        f"equations; each has to say how completely it closes"]
    for (count, word), (name, key) in zip(claims, (("R_1", "scale"),
                                                   ("R_2", "m_scale"))):
        if count in _CLOSURE_COUNTS:
            claimed = _CLOSURE_COUNTS[count] * _CLOSURE_MAGNITUDES[word]
        elif count.replace(",", "").isdigit():
            claimed = int(count.replace(",", "")) * _CLOSURE_MAGNITUDES[word]
        else:
            fails.append(f"{name}'s closure is stated as 'one part in {count} "
                         f"{word}', which is not a number")
            continue
        residual = abs(printed.get(name, state[name.replace("R_", "R")]))
        if not residual:
            continue
        actual = abs(state[key]) / residual
        if claimed > actual * 1.1:
            fails.append(f"{name} is claimed to close to one part in {count} "
                         f"{word}, but {residual:.4g} against {state[key]:.6g} "
                         f"is one part in {actual:.3g} — the section overstates "
                         f"its own closure")
        if claimed < actual / 10:
            fails.append(f"{name} is claimed to close to one part in {count} "
                         f"{word} when it closes to one part in {actual:.3g}; "
                         f"the phrase is more than a factor of ten out")
    return fails


#: Spencer's page in LaTeX, in the notation the report prints. Everything the
#: two published force sums are written with, and nothing else: a macro the map
#: does not carry survives the translation as a backslash and fails, rather than
#: quietly comparing equal to something.
_LATEX_TO_NOTATION = {
    r"\sin": "sin", r"\cos": "cos", r"\beta": "β", r"\psi": "ψ",
    r"\theta_p": "θ_p", r"\delta": "δ", "-": "−",
}


def _page_equation(page, number):
    """One numbered equation of a documentation page, in the report's notation.

    The pin behind the transcription: the report assembles equations (1) and (2)
    from its own registry, and this reads them off the page they came from, so
    that a registry entry renamed, resigned or reordered stops matching the
    published form instead of quietly publishing a new one.
    """
    import re

    found = re.search(r"\$([^$]*?)\\qquad \(%s\)\$" % number, page)
    if not found:
        return None, f"{METHOD_DOC_PAGES['spencer']} has no equation ({number})"
    text = found.group(1)
    for latex, plain in _LATEX_TO_NOTATION.items():
        text = text.replace(latex, plain)
    # A leading unary minus is set against its term. LaTeX spaces it or not as
    # the author typed it — equation (2) is written "$F_v = - W ...$" — and that
    # is typesetting, not notation.
    text = " ".join(text.split()).replace("= − ", "= −")
    if "\\" in text:
        return None, (f"equation ({number}) of "
                      f"{METHOD_DOC_PAGES['spencer']} is written with a macro "
                      f"the notation map does not carry: {text!r}")
    return text, ""


#: The base normal, as the two iterative methods' sections write it. Their
#: sections print it transcribed-then-reduced exactly as they print the equation
#: they solve, so a section carries TWO such pairs and every check that reads one
#: has to say which. The left-hand sides tell them apart: no method's own
#: transcription writes a line beginning ``N' =`` or ``m_α =``, and the march
#: writes ``N'·(tan φ_m …``, which is neither.
_NORMAL_FORCE_PREFIXES = ("N' = ", "m_α = ")


def _transcription_pairs(section):
    """``[(sentence, published, reduced), ...]`` — every published-then-reduced
    pair a section prints, in the order it prints them.

    What marks a pair is what marks it for a reader: the sentence between the two
    forms, which names the forces this model does not carry and the equation that
    loses them. The published form is the run of equations directly above that
    sentence and this model's is the run directly below it.
    """
    blocks = list(section.blocks)
    out = []
    for i, block in enumerate(blocks):
        if block.kind != "prose" or not _REDUCTION_LEAD.search(block.text):
            continue
        published, j = [], i - 1
        while j >= 0 and blocks[j].kind == "math":
            published.insert(0, blocks[j].notation)
            j -= 1
        reduced, j = [], i + 1
        while j < len(blocks) and blocks[j].kind == "math":
            reduced.append(blocks[j].notation)
            j += 1
        out.append((block.text, published, reduced))
    return out


def _is_normal_force(run):
    """Is this run of equations the base normal rather than a method's own?"""
    return bool(run) and all(
        any(n.startswith(p) for p in _NORMAL_FORCE_PREFIXES) for n in run)


def _method_pair(section):
    """``(sentence, published, reduced)`` for the equation the METHOD solves.

    Empty where nothing was dropped: the published form is then this model's and
    stands alone, with no sentence between.
    """
    for pair in _transcription_pairs(section):
        if not _is_normal_force(pair[1]):
            return pair
    return "", [], []


def _normal_force_pair(section):
    """``(sentence, published, reduced)`` for the base normal, or empties."""
    for pair in _transcription_pairs(section):
        if _is_normal_force(pair[1]):
            return pair
    return "", [], []


def _transcription_split(section, prefixes=("",)):
    """``(transcribed, this model's)`` — the equations a section prints as its
    derivation publishes the equation it SOLVES, and every other equation it
    prints.

    Everything else — the base normal reduced to this model, the reduction of the
    method's own equation, the arithmetic — carries only the terms this model
    has, and is what a claim about this model is checked against. Every published
    form is kept out of it, the base normal's included: a transcription is the
    derivation speaking and says nothing about this model.
    """
    _sentence, full, _reduced = _method_pair(section)
    published = {n for _s, pub, _r in _transcription_pairs(section) for n in pub}
    keep = [n for n in full if any(n.startswith(p) for p in prefixes)]
    others = [b.notation for b in section.blocks
              if b.kind == "math" and b.notation not in published
              and any(b.notation.startswith(p) for p in prefixes)]
    return keep, others


#: The sentence between a transcription and this model's reduction of it: it
#: names what the model does not carry and then the equation that loses it.
#: Everything printed below it is the model's own.
_REDUCTION_LEAD = re.compile(r", so equations? \(")


#: The consumers whose equations are printed as a quotient of named sums rather
#: than as one expression. A term of those is a piece of one of several lines,
#: not a whole term of one, so which terms a form carries is read by looking each
#: registry symbol up in the lines rather than by cutting the lines into terms.
_PART_CONSUMERS = ("oms_num", "bishop_num", "page_drv")


def _both_forms(full, reduced, consumers):
    """``(published, reduced)`` — the registry terms each printed form carries.

    Every method but the two moment ones prints one equation per form, and a
    term is what stands between its top-level operators. The moment methods
    print theirs as a quotient of named sums spread over several lines, where
    the same reading would return the part letters; those are read by asking,
    of each symbol the registry declares, whether the lines print it.
    """
    if any(c in _PART_CONSUMERS for c in consumers):
        from xslope.report import FORCE_TERMS, NotApplicable
        declared = []
        for term in FORCE_TERMS:
            for consumer in consumers:
                got = getattr(term, consumer)
                published = (got.published if isinstance(got, NotApplicable)
                             else got)
                declared += [c.symbol for c in published]
        return ([s for s in declared if any(s in n for n in full)],
                [s for s in declared if any(s in n for n in reduced)])
    return ([s for n in full for s in _terms_printed(n)],
            [s for n in reduced for s in _terms_printed(n)])


def _reduction_is_true(section, full, reduced, consumers, where):
    """The sentence between the published equations and this model's reduction of
    them says what actually went, and nothing else.

    Two ways for it to lie, and both are tested. A term the published form
    carries and the reduced form does not, with nothing in the sentence to
    account for it, is a term that vanished off the page. A force the sentence
    says the model does not carry, still standing in the reduced form, is the
    denial that the omission sentence has already been caught making once.

    The forces and the terms they are printed as come from the registry the
    section assembles both forms from, so a force added to it is covered here
    without being named here.
    """
    from xslope.report import FORCE_TERMS, NotApplicable

    fails = []
    owner = {}
    for term in FORCE_TERMS:
        for consumer in consumers:
            got = getattr(term, consumer)
            published = (got.published if isinstance(got, NotApplicable) else got)
            for contribution in published:
                owner[contribution.symbol] = term

    published_terms, reduced_terms = _both_forms(full, reduced, consumers)
    # The sentence that reduces THIS pair. A section that also prints the base
    # normal transcribed-then-reduced carries two, and reading the first one
    # judged the factor-of-safety equation against the base normal's sentence.
    sentence = _method_pair(section)[0]
    dropped = [s for s in published_terms if s not in reduced_terms]

    if not reduced:
        # Nothing was dropped, so the published equations ARE this model's and
        # stand alone. A sentence reducing them would be reducing nothing.
        if sentence:
            fails.append(f"{where}: the section prints no reduced equation and "
                         f"still says {sentence!r}")
        return fails
    if dropped and not sentence:
        fails.append(f"{where}: {dropped} are dropped from the published "
                     f"equations with no sentence saying so")

    for symbol in dropped:
        term = owner.get(symbol)
        if term is None:
            fails.append(f"{where}: {symbol!r} is dropped and belongs to no "
                         f"force in the registry")
        elif term.key == "T_top":
            if "shear on the top of the slice" not in sentence:
                fails.append(f"{where}: {symbol!r} is dropped and the sentence "
                             f"does not say that no shear on the top of the "
                             f"slice is simulated: {sentence!r}")
        elif term.feature not in sentence and symbol not in sentence:
            fails.append(f"{where}: {symbol!r} is dropped from the published "
                         f"equations and the sentence accounts for neither it "
                         f"nor the {term.feature}: {sentence!r}")

    for term in FORCE_TERMS:
        if not term.feature or term.feature not in sentence:
            continue
        carried = [s for s, t in owner.items()
                   if t is term and s in reduced_terms]
        if carried:
            fails.append(f"{where}: the sentence says the model carries no "
                         f"{term.feature}, and the reduced equations print "
                         f"{carried}")
    return fails


def test_spencer_force_sums():
    """Spencer's preamble prints equations (1) and (2) as published, reduces them
    to this model, and a row of the slice table reproduces its own Q.

    The section's claim is that the printed numbers ARE the solution, and Q is
    where a reader would have had to take that on trust: the force sums behind it
    were named in words and shown nowhere. So this does what a reviewer does —
    reads F_h, F_v, α, c, Δl, u and φ off one row, puts them through the Q
    equation as printed at the converged F and θ, and requires the row's own Q_s
    back, to the precision the columns are printed at.

    Above that row-level check sits the transcription. The full forms are pinned
    against the derivation page symbol for symbol, in that page's own order, so
    that the registry the report assembles them from cannot drift from the
    equations it claims to be printing. And the reduction below them is held to
    what actually went: every term the full form carries and the reduced form
    does not has to be accounted for in the sentence between them, and nothing
    the sentence names may still be standing in the reduced form.
    """
    import math
    fails = []

    report, bundle = _calc_report("spencer")
    if report is None:
        return ["the sample model did not solve with Spencer's method"]
    section = _calc_section(report)
    table = next((t for t in report.tables() if t.landscape), None)
    if section is None or table is None:
        return ["there is no calculation or no slice table to read"]

    maths = [b.notation for b in section.blocks if b.kind == "math"]
    full, reduced = _transcription_split(section, ("F_h = ", "F_v = "))
    if len(full) != 2:
        fails.append(f"the preamble does not print equations (1) and (2) of the "
                     f"derivation: {maths}")

    from xslope.report import METHOD_DOC_PAGES
    with open(os.path.join(_REPO, "docs", METHOD_DOC_PAGES["spencer"]),
              encoding="utf-8") as f:
        page = f.read()

    # --- the transcription, against the page it is transcribed from ---
    #
    # Symbol for symbol, in that page's own order: this is what keeps the
    # registry the two forms are assembled from from drifting away from the
    # equations the section says it is printing.
    for notation, number in zip(full, ("1", "2")):
        published, why = _page_equation(page, number)
        if why:
            fails.append(why)
        elif notation != published:
            fails.append(f"the section prints equation ({number}) as "
                         f"{notation!r}; {METHOD_DOC_PAGES['spencer']} publishes "
                         f"{published!r}")

    # --- the reduction, against what it says went ---
    fails += _reduction_is_true(section, full, reduced,
                                ("spencer_h", "spencer_v"), "the sample model")
    if not any("W" in n for n in reduced):
        fails.append(f"neither reduced force sum carries the slice weight: "
                     f"{reduced}")
    # This model has a distributed load and nothing else, so no seismic,
    # reinforcement, pile, line load or tension-crack term survives the
    # reduction.
    for notation in reduced:
        for absent in ("kW", "R cos", "R sin", "H cos", "H sin", "L cos",
                       "L sin", "V"):
            if absent in notation:
                fails.append(f"{notation!r} prints a {absent!r} term for a model "
                             f"with none")

    # And the row-level reproduction, on this model and on a right-facing one —
    # which the section is solved as the mirror image of, and which the preamble
    # therefore has to warn the reader about before they check a row.
    fails += _spencer_rows_reproduce(table, bundle["results"], mirror=1)

    report, bundle = _calc_report("spencer", xlsx=RFACE_XLSX)
    if report is None:
        return fails + ["the right-facing model did not solve with Spencer's "
                        "method"]
    table = next((t for t in report.tables() if t.landscape), None)
    section = _calc_section(report)
    prose = " ".join(b.text for b in section.blocks if b.kind == "prose")
    if "mirror image" not in prose:
        fails.append("a right-facing section does not say that it is solved as "
                     "the mirror image of the slope, which is the only way its "
                     "Q column can be checked against the equations")
    fails += _spencer_rows_reproduce(table, bundle["results"], mirror=-1)

    # The transcription and its reduction on every shape of model there is a
    # Spencer section for: one carrying a distributed load, one carrying water in
    # a tension crack, one right-facing and carrying neither. Each drops a
    # different set of terms, and each sentence has to be true of its own.
    for label, xlsx in (("the right-facing model", RFACE_XLSX),
                        ("the tension-crack model", TENSION_XLSX)):
        report, _bundle = _calc_report("spencer", xlsx=xlsx)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{label} produced no Spencer calculation to reduce")
            continue
        full, reduced = _transcription_split(section, ("F_h = ", "F_v = "))
        for notation, number in zip(full, ("1", "2")):
            published, why = _page_equation(page, number)
            if why:
                fails.append(f"{label}: {why}")
            elif notation != published:
                fails.append(f"{label}: the section prints equation ({number}) "
                             f"as {notation!r}; {METHOD_DOC_PAGES['spencer']} "
                             f"publishes {published!r}")
        fails += _reduction_is_true(section, full, reduced,
                                    ("spencer_h", "spencer_v"), label)
    return fails


def _spencer_rows_reproduce(table, results, mirror):
    """Every row of a Spencer slice table put back through the printed Q
    equation, and the mutation that proves the arithmetic is being done.

    ``mirror`` is -1 on a right-facing surface, where the derivation is applied
    to the mirrored section and α, c and tan φ enter with reversed signs — as
    the preamble states.
    """
    import math
    fails = []

    labels = [h.split(" (")[0] for h in table.headers]
    needed = ("α", "Δl", "u", "c", "φ", "F_h", "F_v", "Q_s")
    if not set(needed) <= set(labels):
        return [f"the slice table is missing {sorted(set(needed) - set(labels))}"]
    at = {name: labels.index(name) for name in needed}
    F = float(results["FS"])
    theta = math.radians(float(results["theta"]))

    def Q_of(cell, F_v):
        a = mirror * math.radians(cell["α"])
        tan_p = mirror * math.tan(math.radians(cell["φ"]))
        c = mirror * cell["c"]
        m_a = 1.0 / (math.cos(a - theta) + math.sin(a - theta) * tan_p / F)
        return (- F_v * math.sin(a) - cell["F_h"] * math.cos(a)
                - c * cell["Δl"] / F
                + (F_v * math.cos(a) - cell["F_h"] * math.sin(a)
                   + cell["u"] * cell["Δl"]) * tan_p / F) * m_a

    # Every operand is printed rounded, and the reproduction carries that
    # rounding: a thousandth of the magnitudes on the row covers the hundredth
    # of a degree α is printed to — the term that dominates — on a slice of any
    # size, with an order of magnitude to spare. The floor is for a sliver at
    # the toe with almost nothing on it.
    def tolerance(cell):
        return 0.2 + 1e-3 * (abs(cell["F_v"]) + abs(cell["F_h"])
                             + (abs(cell["c"]) + abs(cell["u"])) * cell["Δl"])

    caught = False
    for row in table.rows:
        cell = {name: float(row[at[name]]) for name in needed}
        Q = Q_of(cell, cell["F_v"])
        if abs(Q - cell["Q_s"]) > tolerance(cell):
            fails.append(f"slice {row[0]}: the printed F_h, F_v and base "
                         f"properties give Q = {Q:.3f}, the table {cell['Q_s']}")
        # Mutation: an F_v that did not come from the solve would not reproduce.
        # It has to be caught somewhere rather than on a nominated row — a slice
        # whose Q passes through zero is insensitive to F_v on its own, and
        # requiring the failure there would be requiring luck.
        caught = caught or abs(Q_of(cell, cell["F_v"] * 1.01)
                               - cell["Q_s"]) > tolerance(cell)
    if not caught:
        fails.append("an F_v column wrong by a percent still reproduced every "
                     "Q; the check cannot fail")
    return fails


#: The quantities a row's Q is rebuilt from that the mirror sentence has to
#: place, and the words it can name each of them in. Every one of them is either
#: reversed by the reader or already reversed in the column, and the sentence is
#: only true if it says which. "the horizontal forces" is F_h — the phrase the
#: rule was wrong in, which named F_h among the reversals when the column is
#: printed reversed already, so a reader following it reversed F_h twice.
_MIRROR_NAMES = {
    "α": (r"α",),
    "c": (r"c",),
    "tan φ": (r"tan φ",),
    "F_h": (r"F_h", r"horizontal forces?"),
    "F_v": (r"F_v", r"vertical forces?"),
}


def _mirror_rule(prose):
    """The sign rule a right-facing section STATES, as ``{quantity: reversed}``.

    Read out of the sentence itself rather than assumed, so that the arithmetic
    below is the arithmetic the words describe: a quantity named in the sentence
    that reverses signs is reversed, one named in the sentence that says the
    columns enter as printed is not. Returns ``(rule, complaints)``.
    """
    import re

    fails = []
    para = next((s for s in re.split(r"(?<=\.)\s+", prose)
                 if "mirror image" in s), None)
    if para is None:
        return {}, ["no sentence says the section is solved as a mirror image"]
    where = prose[prose.index(para):]
    sentences = [s for s in re.split(r"(?<=\.)\s+", where) if s]
    rule = {}
    for sentence in sentences:
        if "reversed" in sentence:
            reversed_by = True
        elif "as printed" in sentence:
            reversed_by = False
        else:
            continue
        for quantity, spellings in _MIRROR_NAMES.items():
            if quantity in rule:
                continue
            for spelling in spellings:
                if re.search(rf"(?<![A-Za-zα-ωΑ-Ω_]){spelling}"
                             rf"(?![A-Za-zα-ωΑ-Ω_])", sentence):
                    rule[quantity] = reversed_by
                    break
    # A quantity the sentence never places is a quantity a reader would take as
    # it is printed. That is a defect in the sentence and it is reported, and the
    # arithmetic is still run under that reading, so a wording that misplaces one
    # quantity and omits another fails on both counts.
    missing = [q for q in _MIRROR_NAMES if q not in rule]
    if missing:
        fails.append(f"the mirror sentence never says whether {missing} enter "
                     f"the equations reversed or as printed: {para!r}")
    return rule, fails


def test_spencer_mirror_rule():
    """A right-facing section's mirror sentence IS the arithmetic behind its Q_s
    column.

    The rule is parsed out of the sentence and applied to the printed columns,
    and the row's own Q has to come back. The model is one with real horizontal
    forces on its slices: where F_h is zero on every slice — the sample
    right-facing slope — a rule that reverses F_h and a rule that does not give
    the same Q, and the sentence can say either and read correct.
    """
    import math
    fails = []

    report, bundle = _calc_report("spencer", xlsx=MIRROR_XLSX)
    if report is None:
        return ["the right-facing model did not solve with Spencer's method"]
    section = _calc_section(report)
    table = next((t for t in report.tables() if t.landscape), None)
    if section is None or table is None:
        return ["there is no calculation or no slice table to read"]

    prose = " ".join(b.text for b in section.blocks if b.kind == "prose")
    rule, fails = _mirror_rule(prose)
    if not rule:
        return fails

    labels = [h.split(" (")[0] for h in table.headers]
    needed = ("α", "Δl", "u", "c", "φ", "F_h", "F_v", "Q_s")
    if not set(needed) <= set(labels):
        return [f"the slice table is missing {sorted(set(needed) - set(labels))}"]
    at = {name: labels.index(name) for name in needed}
    F = float(bundle["results"]["FS"])
    theta = math.radians(float(bundle["results"]["theta"]))

    def Q_of(cell, rule):
        def signed(quantity, value):
            return -value if rule.get(quantity) else value
        a = signed("α", math.radians(cell["α"]))
        tan_p = signed("tan φ", math.tan(math.radians(cell["φ"])))
        c = signed("c", cell["c"])
        F_h = signed("F_h", cell["F_h"])
        F_v = signed("F_v", cell["F_v"])
        m_a = 1.0 / (math.cos(a - theta) + math.sin(a - theta) * tan_p / F)
        return (- F_v * math.sin(a) - F_h * math.cos(a) - c * cell["Δl"] / F
                + (F_v * math.cos(a) - F_h * math.sin(a)
                   + cell["u"] * cell["Δl"]) * tan_p / F) * m_a

    # The tolerance of the sibling check: the operands are printed rounded and
    # the reproduction carries that rounding.
    def tolerance(cell):
        return 0.2 + 1e-3 * (abs(cell["F_v"]) + abs(cell["F_h"])
                             + (abs(cell["c"]) + abs(cell["u"])) * cell["Δl"])

    cells = [{name: float(row[at[name]]) for name in needed} for row in table.rows]
    if not any(abs(cell["F_h"]) > 1e-6 for cell in cells):
        return ["every printed F_h is zero on this model, so the mirror rule "
                "cannot be told from a rule that reverses the horizontal forces"]
    reproduced = True
    for row, cell in zip(table.rows, cells):
        Q = Q_of(cell, rule)
        if abs(Q - cell["Q_s"]) > tolerance(cell):
            reproduced = False
            fails.append(f"slice {row[0]}: the rule the mirror sentence states "
                         f"gives Q = {Q:.4f} from the printed columns, the "
                         f"table {cell['Q_s']}")
    # Mutation: the same rule with F_h reversed the other way — which is the
    # rule the sentence used to state — must NOT reproduce, or the sentence's
    # F_h clause is a clause nothing here tests.
    other = dict(rule, F_h=not rule.get("F_h"))
    if reproduced and all(abs(Q_of(cell, other) - cell["Q_s"]) <= tolerance(cell)
                          for cell in cells):
        fails.append("the same Q comes back with F_h reversed the other way, so "
                     "what the sentence says about the horizontal forces is "
                     "not tested")
    return fails


#: The letters each force of the general equation is printed as, so that a
#: sentence saying the model carries none of one can be tested against the
#: equations under it. Spencer's page follows UTEXAS and gives three of them
#: other letters, which is why the two are listed apart: P is the distributed
#: load there and the reinforcement everywhere else.
_FEATURE_SYMBOLS = {
    "distributed load": ("D",),
    "seismic load": ("kW",),
    "tension-crack water force": ("T",),
    "reinforcement crossing the failure surface": ("P", "P_p"),
    "pile force": ("H", "H_p"),
    "line load": ("L",),
}
_SPENCER_FEATURE_SYMBOLS = dict(
    _FEATURE_SYMBOLS,
    **{"distributed load": ("P",),
       "tension-crack water force": ("V",),
       "reinforcement crossing the failure surface": ("R",)})


def test_absent_features_are_really_absent():
    """No section says the model carries a force its own equations print.

    The omission sentence is what tells a reader that a missing term is a term
    the model does not have. A model reinforced entirely by passive capacity has
    nothing in the tangent and axial columns and its P_p terms in the quotient
    directly below, and the sentence denied the reinforcement anyway. So every
    force the sentence names as absent is required absent from every equation the
    section prints.
    """
    fails = []
    from xslope.report import EQUATION_SYMBOLS, _present_symbols

    models = [("the passive-reinforcement model", PASSIVE_XLSX,
               ("bishop", "oms")),
              ("the right-facing reinforced model",
               os.path.join(_REPO, "docs", "lem", "files",
                            "xslope_reinforce_rface.xlsx"), ("bishop", "oms")),
              ("the axially reinforced model", AXIAL_XLSX, ("bishop", "oms")),
              ("the tension-crack model", TENSION_XLSX, CALC_METHODS),
              ("the sample model", REINF_XLSX, CALC_METHODS)]

    claimed = set()
    for where, xlsx, methods in models:
        for method in methods:
            report, _bundle = _calc_report(method, xlsx=xlsx)
            if report is None:
                continue
            section = _calc_section(report)
            if section is None:
                continue
            intro = next((b.text for b in section.blocks
                          if b.kind == "prose" and "carries no" in b.text), "")
            names = (_SPENCER_FEATURE_SYMBOLS if method == "spencer"
                     else _FEATURE_SYMBOLS)
            absent = [name for name in names if name in intro]
            claimed.update(absent)
            # Every section transcribes its equations in full before reducing
            # them, and the transcription is the derivation's — it carries every
            # force there is, and says nothing about this model. What the
            # sentence is a claim about is the reduction below it.
            _full, model_own = _transcription_split(section)
            text = " ".join(model_own)
            # Every symbol the report knows goes in as a candidate, not just the
            # ones being looked for: the longest match wins, and R_1 has to claim
            # its characters before the reinforcement force R can be read out of
            # them.
            present = _present_symbols(
                text, [s for syms in names.values() for s in syms]
                + list(EQUATION_SYMBOLS))
            for name in absent:
                carried = [s for s in names[name] if s in present]
                if carried:
                    fails.append(f"{where} under {method}: the section says it "
                                 f"carries no {name}, and its equations print "
                                 f"{carried}")
    # A passive-reinforced model has to be one of the models this ran on, and it
    # has to print its reinforcement: without that the check above passes on a
    # sentence that names nothing.
    report, _bundle = _calc_report("bishop", xlsx=PASSIVE_XLSX)
    section = _calc_section(report) if report is not None else None
    if section is None:
        fails.append("the passive-reinforcement model produced no calculation, "
                     "so nothing here tests a passive force against its prose")
    else:
        text = " ".join(b.notation for b in section.blocks if b.kind == "math")
        if "P_p" not in text:
            fails.append("the passive-reinforcement model prints no P_p term")
    if not claimed:
        fails.append("no section named an absent force, so nothing was tested")
    return fails


#: Every equation number a method's Calculations section names, and the
#: documentation page each one is published on. Written out here rather than read
#: from the report: a section that starts citing a different equation has to
#: fail, not pass on its own new claim.
#:
#: Corps of Engineers, Lowe & Karafiath and Morgenstern-Price transcribe the
#: per-slice equilibrium of the force-equilibrium derivation, which is the page
#: the first two are documented on and the page Morgenstern-Price's own sends the
#: reader to for its march.
#:
#: The two moment methods each print their OWN page's factor-of-safety equation
#: — the Ordinary Method's (8), Bishop's (10) — and then the sums the arithmetic
#: is formed from, which are that equation in the general moment arms. The
#: Ordinary Method numbers that form (8a); Bishop's page gives it under
#: Composite Surfaces without a number, so Bishop's section names no number for
#: it. Neither section may cite the other's page: an equation of Bishop's method
#: assembled out of the Ordinary Method's numbering is what this is here to stop.
_EQUATION_NUMBERS = {
    "oms": (("lem/oms.md", ("4", "8", "8a")),),
    "bishop": (("lem/bishop.md", ("8", "10")),),
    "janbu": (("lem/janbu.md", ("1", "4", "5", "6", "7")),),
    "corps": (("lem/force_eq.md", ("6", "7", "12")),),
    "lowe": (("lem/force_eq.md", ("6", "7", "12")),),
    "mprice": (("lem/force_eq.md", ("6", "7", "12")),
               ("lem/mprice.md", ("2", "4", "8"))),
    "spencer": (("lem/spencer.md", ("1", "2", "23", "24", "27", "28")),),
}


def _terms_printed(notation):
    """Every term one printed equation carries, as the registry writes it.

    A term is what stands between the top-level operators: inside a Σ where the
    equation sums over the slices, bare where it is written for one slice, and
    inside a mobilized fraction where the equation divides it by F. Whichever it
    is, what comes back is the symbol the registry declares, so a term can be
    followed from the published form into the model's own however the two are
    arranged.
    """
    def bare(piece):
        piece = piece.strip().lstrip("+−").strip()
        if piece.startswith("sum{") and piece.endswith("}"):
            piece = piece[4:-1]
        if piece.startswith("frac{") and piece.endswith("}{F}"):
            piece = piece[5:-4]
        return piece

    if "sum{" in notation:
        out, at = [], notation.find("sum{")
        while at >= 0:
            depth, j = 0, at + 3
            while j < len(notation):
                if notation[j] == "{":
                    depth += 1
                elif notation[j] == "}":
                    depth -= 1
                    if not depth:
                        break
                j += 1
            out.append(notation[at:j + 1])
            at = notation.find("sum{", j + 1)
        return [bare(p) for p in out]
    body = notation.split(" = ", 1)[-1].strip()
    if body == "0":
        return []
    return [bare(p) for p in re.split(r" [−+] ", body.lstrip("−"))]


#: How the prose names a documentation equation. One form, so that a number
#: reaching the page any other way is a number this cannot check.
_EQUATION_REFERENCE = re.compile(
    r"equations? \(([0-9a-z]+)\)(?: and \(([0-9a-z]+)\))?", re.I)

#: And how a page numbers one. The pages set the number with \qquad or \quad, so
#: the spacing is not what a reference is looked up by.
_PAGE_NUMBER = r"\\q?quad ?\(%s\)"


def _cited_equations(section, method):
    """``[(page, number), ...]`` — every equation number the section names, each
    against the page the sentence that names it resolves on.

    A section that transcribes from two derivations cites numbers from two
    documents, and both documents number an equation (6). What decides which is
    meant is the sentence: one that links a page cites THAT page's numbering, and
    one that links none continues the page the sentence before it established —
    a reduction reads on from the equations introduced above it — starting from
    the method's own derivation. Reading the numbers without their page lets a
    citation drift to the other document and still validate, because the number
    exists on both.
    """
    from xslope.report import METHOD_DOC_PAGES, docs_url

    known = set(METHOD_DOC_PAGES.values())
    for entries in _EQUATION_NUMBERS.values():
        known |= {path for path, _numbers in entries}
    by_url = {docs_url(path).rstrip("/"): path for path in known}
    own = METHOD_DOC_PAGES.get(method, "")

    out, page = [], own
    for block in section.blocks:
        if block.kind != "prose":
            continue
        for _text, target in (getattr(block, "links", None) or ()):
            linked = by_url.get(str(target).rstrip("/"))
            if linked:
                page = linked
                break
        out += [(page, n) for match in _EQUATION_REFERENCE.finditer(block.text)
                for n in match.groups() if n]
    return out


def _citations_are_declared(method, cited, wanted, page_of):
    """The pages and numbers the section cites are the ones declared for it, and
    each number is published on the page it is cited against."""
    fails = []
    if cited != wanted:
        fails.append(f"{method}: the prose cites {sorted(cited)}; the section "
                     f"transcribes {sorted(wanted)}")
    for path, number in sorted(cited | wanted):
        if not re.search(_PAGE_NUMBER % re.escape(number), page_of(path)):
            fails.append(f"{method}: the prose names equation ({number}) of "
                         f"{path}, which publishes no such number")
    return fails


def test_equation_numbers_are_in_the_prose():
    """No equation is printed with a number beside it, and every number the prose
    names is published where it says.

    A report numbers its own equations consecutively. These come from seven
    derivations that number theirs independently, so printing those numbers on
    the page runs (1), (2), (7), (1) down a section with two different equations
    under the same number. The number belongs in the sentence that introduces the
    equation — "equation (7) of the derivation" — where it reads as the reference
    it is, and the section's opening statement links the page it resolves on.

    So: nothing beside the equations, every number named in the prose really on
    the page it is attributed to, and every number a section is expected to name
    really named.
    """
    fails = []
    from xslope.report import METHOD_DOC_PAGES, Math

    # The mechanism itself: an equation has nothing to print a number with.
    if "label" in getattr(Math, "__dataclass_fields__", {}):
        fails.append("Math still carries a label; an equation number can be set "
                     "beside the equation again")

    pages = {}

    def page_of(path):
        if path not in pages:
            with open(os.path.join(_REPO, "docs", path), encoding="utf-8") as f:
                pages[path] = f.read()
        return pages[path]

    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the numbering of")
            continue
        for block in section.blocks:
            if block.kind == "math" and getattr(block, "label", ""):
                fails.append(f"{method}: {block.notation!r} is printed with "
                             f"{block.label!r} beside it")
        declared = dict(_EQUATION_NUMBERS.get(method, ()))
        wanted = {(path, n) for path, numbers in declared.items()
                  for n in numbers}
        cited = set(_cited_equations(section, method))
        fails += _citations_are_declared(method, cited, wanted, page_of)
        # And an equation the section prints has to be introduced by a sentence
        # that says which one it is: a number named nowhere near the equation is
        # a reference the reader cannot follow.
        if wanted and not any(
                b.kind == "prose" and _EQUATION_REFERENCE.search(b.text)
                and b.text.rstrip().endswith(":")
                for b in section.blocks):
            fails.append(f"{method}: no sentence introduces an equation by the "
                         f"number the derivation gives it")
        if method not in _EQUATION_NUMBERS:
            continue
        if not METHOD_DOC_PAGES.get(method):
            fails.append(f"{method}: names equation numbers with no page mapped")

    # The mutation: a citation that drifted to the OTHER page the section
    # transcribes from. Morgenstern-Price cites (8) of its own derivation and
    # (6) and (7) of the force-equilibrium one, and both pages number an
    # equation (6), (7) and (8) — so a number checked without its page passes
    # wherever it is pointed. This moves (8) onto the force-equilibrium page,
    # which publishes an equation (8) of its own, and requires the drift to be
    # caught anyway.
    declared = dict(_EQUATION_NUMBERS["mprice"])
    wanted = {(path, n) for path, numbers in declared.items() for n in numbers}
    drifted = {("lem/force_eq.md", n) if (path, n) == ("lem/mprice.md", "8")
               else (path, n) for path, n in wanted}
    if drifted == wanted:
        fails.append("the mutation moved no citation, so it tests nothing")
    elif not _citations_are_declared("mprice", drifted, wanted, page_of):
        fails.append("a citation moved to the other page the section "
                     "transcribes from still validated")
    return fails


def test_the_published_equation_comes_first():
    """Every method prints the equation its derivation publishes, in full, before
    it prints this model's.

    A section that opens on the reduced form shows a reader an equation with the
    seismic term missing and nothing to say a seismic term exists. So each of the
    seven transcribes the published equation first — every force a slice can
    take, assembled from the registry the model's own equation is assembled from
    — then states in one sentence what this model does not carry, then prints
    what is left.

    Three things are required of that transcription. It carries every
    contribution the registry declares for the equation, so it is the published
    form and not a longer reduction. The forces it carries beyond the ones the
    page's own equation writes are only those the registry marks passive, which
    the pages describe in prose and put on the mobilized side. And the sentence
    below it accounts for every term that then goes.
    """
    fails = []
    from xslope.report import FORCE_TERMS, NotApplicable, TRANSCRIPTIONS

    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the transcription of")
            continue
        spec = TRANSCRIPTIONS[method]
        full, model_own = _transcription_split(section)
        reduced = bool(full)
        if not reduced:
            # Nothing was dropped: the published equation IS this model's, and
            # the section prints it once with no sentence between.
            full = model_own
        printed, _model = _both_forms(full, [], spec.consumers)
        for term in FORCE_TERMS:
            for consumer in spec.consumers:
                got = getattr(term, consumer)
                published = (got.published if isinstance(got, NotApplicable)
                             else got)
                for contribution in published:
                    if contribution.symbol not in printed:
                        fails.append(
                            f"{method}: the published {consumer} carries "
                            f"{contribution.symbol!r} and the transcription "
                            f"does not: {full}")
        # And nothing in it that the registry does not declare for those
        # equations — a term from another equation would be a claim the page
        # never made.
        declared = {c.symbol for t in FORCE_TERMS for consumer in spec.consumers
                    for c in (getattr(t, consumer).published
                              if isinstance(getattr(t, consumer), NotApplicable)
                              else getattr(t, consumer))}
        # The interslice resultant carried in from the slice before is not a
        # force the model applies but what the march computes, so it is written
        # by the equation rather than declared in the registry. Named here so
        # that it is the only thing the registry does not have to account for.
        declared |= {"Z_i·cos θ_i", "Z_i·sin θ_i"}
        if reduced and not any(c in _PART_CONSUMERS for c in spec.consumers):
            # The two moment methods print theirs as a quotient of named sums,
            # where a line carries part of a term rather than a term; what keeps
            # anything the registry does not declare out of THEIR transcription
            # is the recomposition, which has to give the published equation
            # term for term (test_the_moment_quotient_recomposes).
            for symbol in printed:
                if symbol not in declared:
                    fails.append(f"{method}: the transcription prints "
                                 f"{symbol!r}, which the registry declares for "
                                 f"no equation it transcribes")
        if reduced:
            fails += _reduction_is_true(section, full, model_own,
                                        spec.consumers, method)
    return fails


#: Every LaTeX the derivation pages write their transcribed equations with, in
#: the notation the report prints. Applied longest key first. Anything left with
#: a backslash after the pass has not been translated, and the comparisons below
#: refuse rather than quietly matching something else.
_PAGE_LATEX = {
    r"\Delta \ell": "Δl", r"\Delta\ell": "Δl", r"\ell": "l",
    r"\theta_p": "θ_p", r"\theta_{i+1}": "θ_{i+1}", r"\theta_{i}": "θ_i",
    r"\theta_j": "θ_j",
    r"\theta_i": "θ_i", r"\alpha": "α", r"\beta": "β", r"\phi'": "φ",
    r"\phi_m": "φ_m", r"\phi": "φ", r"\psi": "ψ", r"\delta": "δ",
    r"\lambda": "λ", r"\pi": "π",
    r"\dfrac": "frac", r"\frac": "frac", r"\sum": "sum",
    r"\sin": "sin", r"\cos": "cos", r"\tan": "tan",
    r"\left[": "[", r"\right]": "]", r"\left(": "(", r"\right)": ")",
    r"\,": " ", r"\;": " ", r"\ ": " ",
    "-": "−",
}

#: What a term is compared on: the symbols and operators, with the grouping and
#: the spacing taken out. The report writes ``frac{1}{R}·sum{D cos β·a_dx}`` and
#: the page ``\frac{1}{R}\sum D \cos \beta \, a_{dx}``; they are the same term,
#: and the two ways of setting it are not what this is checking.
_GROUPING = "{}[]() ·\t\n"


def _page_latex(text):
    """One page's equation in the report's notation, or ``(None, why)``."""
    out = text
    for latex in sorted(_PAGE_LATEX, key=len, reverse=True):
        out = out.replace(latex, _PAGE_LATEX[latex])
    out = re.sub(r"_\{([A-Za-z0-9+]+)\}", r"_\1", out)
    if "\\" in out:
        return None, f"a macro the notation map does not carry: {out!r}"
    return " ".join(out.split()), ""


def _signed_terms(text):
    """``[(sign, term), ...]`` — one side of an equation cut at its top-level
    operators, with everything inside a fraction, a sum or a bracket left whole.
    """
    out, depth, start, sign = [], 0, 0, +1
    for i, ch in enumerate(text):
        if ch in "{[(":
            depth += 1
        elif ch in "}])":
            depth -= 1
        elif depth == 0 and ch in "+−" and i > start:
            out.append((sign, text[start:i]))
            sign, start = (+1 if ch == "+" else -1), i + 1
    out.append((sign, text[start:]))
    return [(s, t) for s, t in
            [(s, "".join(c for c in t if c not in _GROUPING)) for s, t in out]
            if t]


def _as_quotient(text, where):
    """``(numerator, denominator)`` of ``F = frac{A}{B}``."""
    body = text.split("=", 1)[-1].strip()
    if not body.startswith("frac{"):
        return None, None, f"{where} is not a quotient: {text!r}"
    pieces, at = [], len("frac")
    for _side in (0, 1):
        if at >= len(body) or body[at] != "{":
            return None, None, f"{where} is not a quotient: {text!r}"
        depth, j = 0, at
        while j < len(body):
            if body[j] == "{":
                depth += 1
            elif body[j] == "}":
                depth -= 1
                if not depth:
                    break
            j += 1
        pieces.append(body[at + 1:j])
        at = j + 1
    return pieces[0], pieces[1], ""


def _recompose(lines):
    """The named parts of one printed quotient substituted back into it.

    ``lines`` are the equations the section prints for the published form: the
    quotient itself, the numerator's sum, the group the base normal is formed
    from, and one part per force. Returns ``(numerator, denominator, why)`` with
    every part letter replaced by what the section defines it as — which is what
    is then held against the page.

    A part standing as a whole term of the quotient is substituted bare, so that
    the terms inside it become terms of the quotient, as they are on the page.
    The group inside the numerator is multiplied by tan φ, so it is bracketed —
    also as the page writes it.
    """
    defined = {}
    top = ""
    for line in lines:
        name, _sep, body = line.partition(" = ")
        if name == "F":
            top = line
        else:
            defined[name] = body
    if not top:
        return None, None, f"no quotient among the printed equations: {lines}"
    num, den, why = _as_quotient(top, "the printed quotient")
    if why:
        return None, None, why

    def substitute(text, bracket):
        for name in sorted(defined, key=len, reverse=True):
            if name in text:
                body = defined[name]
                text = text.replace(name, f"({body})" if bracket else body)
        return text

    # The quotient's own parts first, bare; then the group inside the numerator,
    # bracketed, which is the only substitution that lands inside a product.
    num, den = substitute(num, False), substitute(den, False)
    return substitute(num, True), substitute(den, True), ""


def _page_quotient_sides(page, number, where):
    """One numbered equation of a documentation page as ``(numerator terms,
    denominator terms)``, in the report's notation."""
    found = re.search(r"\$([^$]*?)\\q?quad ?\(%s\)\$" % re.escape(number), page)
    if not found:
        return None, None, f"{where} publishes no equation ({number})"
    text, why = _page_latex(found.group(1))
    if why:
        return None, None, f"equation ({number}) of {where} is written with {why}"
    num, den, why = _as_quotient(text, f"equation ({number}) of {where}")
    if why:
        return None, None, why
    return _signed_terms(num), _signed_terms(den), ""


def _recomposes_to(lines, page, number, where):
    """The printed parts, put back together, against the published equation."""
    fails = []
    num, den, why = _recompose(lines)
    if why:
        return [f"{where}: {why}"]
    want_num, want_den, why = _page_quotient_sides(page, number, where)
    if why:
        return [f"{where}: {why}"]
    for side, got, want in (("numerator", _signed_terms(num), want_num),
                            ("denominator", _signed_terms(den), want_den)):
        if sorted(got) != sorted(want):
            missing = [t for t in want if t not in got]
            extra = [t for t in got if t not in want]
            fails.append(
                f"{where}: the named parts recompose to a {side} of {got}; "
                f"equation ({number}) publishes {want}"
                + (f" — missing {missing}" if missing else "")
                + (f", printed instead {extra}" if extra else ""))
    return fails


#: Which equation each moment method's Calculations section prints in named
#: parts, on which page. The pin: the parts substituted back have to give this
#: equation, term for term and sign for sign, on the page's own side of the bar.
_RECOMPOSED = {"oms": ("lem/oms.md", "8"), "bishop": ("lem/bishop.md", "10")}


def _printed_parts(section):
    """The FIRST quotient-in-named-parts a section prints — the published form.

    Read by shape rather than by the sentence below it: the quotient opens the
    run, and every line after it that defines one of its letters belongs to it.
    A section that drops nothing prints the run once and has no sentence, and
    this finds it either way.
    """
    # The part letters come from the module that prints them, so a part renamed
    # there is a part this still finds.
    from xslope.report import MOMENT_PART_SYMBOLS

    maths = [b.notation for b in section.blocks if b.kind == "math"]
    at = next((i for i, n in enumerate(maths) if n.startswith("F = frac{N_S}")),
              None)
    if at is None:
        return []
    out = [maths[at]]
    for notation in maths[at + 1:]:
        name = notation.partition(" = ")[0]
        if name not in ("N_S", "N_v", "N'") and name not in MOMENT_PART_SYMBOLS:
            break
        out.append(notation)
    return out


def _force_of_symbol():
    """``{symbol: force key}`` over every equation the registry feeds.

    One force is written differently in different equations — the reinforcement
    is ``P sin ψ`` in the base-normal group, ``P cos ψ·a_ry + P sin ψ·a_rx`` in
    the moment sums and ``P·a_S`` where it is tangent — so a term is followed
    from one printed form to another by the FORCE it belongs to and not by the
    letters it is set in.
    """
    from xslope.report import CONSUMERS, FORCE_TERMS, NotApplicable

    out = {}
    for term in FORCE_TERMS:
        for consumer in CONSUMERS:
            got = getattr(term, consumer)
            published = (got.published if isinstance(got, NotApplicable)
                         else got)
            for contribution in published:
                out[contribution.symbol] = term.key
    return out


#: A force whose printed form above is not one of its registry symbols, and what
#: introduces it instead. Both pages write the cohesion on the base OUTSIDE the
#: group the normal force is formed from, so it belongs to the numerator's own
#: frame rather than to a registry contribution — and the strength arrives on the
#: page as c·Δl, under whichever arm the equation it is in gives it.
_INTRODUCED_BY = {"strength": ("c·Δl",)}


def _terms_are_introduced(evaluated, above, where):
    """Every force in the equation the arithmetic is formed from is one the
    equations above it have already put on the page."""
    owner = _force_of_symbol()
    fails = []
    introduced = {key for symbol, key in owner.items()
                  if any(symbol in line for line in above)}
    introduced |= {key for key, letters in _INTRODUCED_BY.items()
                   if any(letter in line for letter in letters
                          for line in above)}
    for symbol in _terms_printed(evaluated):
        key = owner.get(symbol)
        if key is None:
            fails.append(f"{where}: the evaluated equation prints {symbol!r}, "
                         f"which belongs to no force in the registry")
        elif key not in introduced:
            fails.append(f"{where}: the evaluated equation prints {symbol!r} "
                         f"and no equation above it carries the {key} force: "
                         f"{above}")
    return fails


def test_the_evaluated_equation_introduces_its_terms():
    """Nothing is used in the equation the arithmetic follows that the equations
    above it have not already shown.

    The two moment methods print their page's own quotient, this model's
    reduction of it, and then the same equilibrium in the general moment arms —
    and it is that last one the two sums are formed from. A force reaching it
    without appearing in either form above is a number with no path back to an
    equation: the passive reinforcement moment did exactly that, standing for
    about 1700 of a 20810 numerator on the passive benchmark while the sentence
    directly above it said the reinforcement terms were zero on every slice.

    The pages write their factor-of-safety equations for the active case and
    give the passive rule in prose — capacity that mobilizes with the soil joins
    the mobilized side and is divided by F — so the reduced form names that
    moment before the evaluated form uses it.
    """
    fails = []
    models = (("the sample model", REINF_XLSX),
              ("the passive-reinforcement model", PASSIVE_XLSX),
              ("the axially reinforced model", AXIAL_XLSX))
    exercised = False

    for label, xlsx in models:
        for method in _RECOMPOSED:
            report, _bundle = _calc_report(method, xlsx=xlsx)
            section = _calc_section(report) if report is not None else None
            if section is None:
                fails.append(f"{label} under {method}: no calculation")
                continue
            maths = [b.notation for b in section.blocks if b.kind == "math"]
            at = max((i for i, n in enumerate(maths)
                      if n.startswith("F = frac{sum{")), default=None)
            if at is None:
                fails.append(f"{label} under {method}: no equation for the "
                             f"arithmetic to follow: {maths}")
                continue
            if "P_p" in maths[at] or "H_p" in maths[at]:
                exercised = True
            fails += _terms_are_introduced(maths[at], maths[:at],
                                           f"{label} under {method}")

    if not exercised:
        fails.append("no model put passive support in an evaluated equation, "
                     "so the case this is here for was not reached")

    # The mutation: take the part that introduces the passive moment back out,
    # which is the state this was written against.
    report, _bundle = _calc_report("bishop", xlsx=PASSIVE_XLSX)
    section = _calc_section(report)
    maths = [b.notation for b in section.blocks if b.kind == "math"]
    at = max(i for i, n in enumerate(maths) if n.startswith("F = frac{sum{"))
    stripped = [n for n in maths[:at] if not n.startswith("N_P = ")]
    if len(stripped) == at:
        fails.append("the passive-reinforcement model prints no N_P part, so "
                     "the mutation removed nothing")
    elif not _terms_are_introduced(maths[at], stripped, "the mutation"):
        fails.append("an evaluated equation using a force no equation above it "
                     "carries was not caught")
    return fails


def test_the_moment_quotient_recomposes():
    """The two moment methods print their own page's factor-of-safety equation.

    Equation (8) of the Ordinary Method of Slices derivation and equation (10)
    of Bishop's each run several times the width of a page, and what went onto
    the page instead was a reconstruction: Bishop's section printed a moment
    balance rewritten in another derivation's moment arms, which is not Bishop's
    equation and does not carry its number.

    Each now prints its own equation as a quotient of named sums, one part per
    force. This is what holds those parts to the page: they are substituted back
    into the quotient and the result compared with the published equation, term
    for term and sign for sign, numerator against numerator. A term dropped from
    a part, or a part moved from one side of the bar to the other, cannot
    survive it.
    """
    fails = []

    for method, (path, number) in _RECOMPOSED.items():
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to recompose")
            continue
        full = _printed_parts(section)
        if len(full) < 3:
            printed = [b.notation for b in section.blocks if b.kind == "math"]
            fails.append(f"{method}: the section prints no quotient in named "
                         f"parts: {printed}")
            continue
        with open(os.path.join(_REPO, "docs", path), encoding="utf-8") as f:
            page = f.read()
        fails += _recomposes_to(full, page, number, path)

    # The mutations: the comparison has to catch a part that lost a term and a
    # part that changed sides. Run on Bishop's own printed parts, so what is
    # mutated is the thing the check reads.
    report, _bundle = _calc_report("bishop")
    section = _calc_section(report) if report is not None else None
    if section is None:
        return fails + ["there is no Bishop calculation to mutate"]
    full = _printed_parts(section)
    with open(os.path.join(_REPO, "docs", "lem", "bishop.md"),
              encoding="utf-8") as f:
        page = f.read()

    dropped = [re.sub(r" − frac\{1\}\{R\}·sum\{D sin β·a_dy\}", "", line)
               for line in full]
    if dropped == full:
        fails.append("the mutation removed nothing, so it tests nothing")
    elif not _recomposes_to(dropped, page, "10", "the mutation"):
        fails.append("a term dropped from a named part still recomposed to "
                     "equation (10)")

    moved = [line.replace("F = frac{N_S}{D_W", "F = frac{N_S − D_P}{D_W")
             .replace(" − D_P − D_H", " − D_H") for line in full]
    if moved == full:
        fails.append("the mutation moved nothing, so it tests nothing")
    elif not _recomposes_to(moved, page, "10", "the mutation"):
        fails.append("a part moved across the fraction bar still recomposed to "
                     "equation (10)")
    return fails


#: The four methods whose transcription is one printed equation per published
#: one, and where those are published: the page, the numbers in the order the
#: section prints them, and the identities the page itself states between its
#: letters and the report's. Each identity is checked to be ON the page before it
#: is used, so a substitution here cannot excuse a notation the page never
#: introduced.
#:
#: The other three are pinned elsewhere and for the same purpose: the two moment
#: methods by recomposing their named parts (test_the_moment_quotient_recomposes)
#: and Spencer's two force sums against its own page (test_spencer_force_sums).
_FULL_FORMS = {
    "janbu": ("lem/janbu.md", ("7",),
              (("sumNsinα", "sumN'+uΔlsinα", r"N = N' + u\,\Delta\ell"),)),
    "corps": ("lem/force_eq.md", ("6", "7"), ()),
    "lowe": ("lem/force_eq.md", ("6", "7"), ()),
    "mprice": ("lem/force_eq.md", ("6", "7"), ()),
}


def _canonical(text):
    """One equation reduced to what it says: the symbols and the operators, with
    the grouping and the spacing taken out."""
    out, why = _page_latex(text)
    if why:
        return None, why
    return "".join(c for c in out if c not in _GROUPING), ""


def test_the_full_forms_match_their_pages():
    """Every method prints its own page's equation, symbol for symbol.

    A section that transcribes the published equation is making a claim about
    another document, and the claim is checkable. Each of the four methods that
    print one equation per published one has its printed form reduced to its
    symbols and operators and held against the same reduction of the equation on
    the page it names — the whole equation, not a list of the letters in it.

    Janbu's page writes the total base normal N and states, in the sentence under
    equation (7), that N = N' + u·Δl; the report writes it out. That identity is
    the one substitution made here, and the page is required to carry it.
    """
    fails = []

    for method, (path, numbers, identities) in _FULL_FORMS.items():
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the transcription of")
            continue
        full, _model_own = _transcription_split(section)
        if len(full) != len(numbers):
            fails.append(f"{method}: the section transcribes {full}, and "
                         f"{path} is cited for {len(numbers)} equation(s)")
            continue
        with open(os.path.join(_REPO, "docs", path), encoding="utf-8") as f:
            page = f.read()
        for number, printed in zip(numbers, full):
            found = re.search(r"\$([^$]*?)\\q?quad ?\(%s\)\$" % re.escape(number),
                              page)
            if not found:
                fails.append(f"{method}: {path} publishes no equation ({number})")
                continue
            want, why = _canonical(found.group(1))
            if why:
                fails.append(f"{method}: equation ({number}) of {path} is "
                             f"written with {why}")
                continue
            for was, becomes, stated in identities:
                if stated not in page:
                    fails.append(f"{method}: {path} does not state {stated!r}, "
                                 f"so {was!r} may not be read as {becomes!r}")
                    continue
                want = want.replace(was, becomes)
            got, why = _canonical(printed)
            if why:
                fails.append(f"{method}: the section prints {printed!r}, which "
                             f"carries {why}")
            elif got != want:
                fails.append(f"{method}: the section prints equation ({number}) "
                             f"as {got!r}; {path} publishes {want!r}")

    # The mutation: a printed form that has drifted from its page has to be
    # caught. One term dropped from the horizontal balance Janbu publishes.
    report, _bundle = _calc_report("janbu")
    section = _calc_section(report)
    full, _model_own = _transcription_split(section)
    with open(os.path.join(_REPO, "docs", "lem", "janbu.md"),
              encoding="utf-8") as f:
        page = f.read()
    found = re.search(r"\$([^$]*?)\\q?quad ?\(7\)\$", page)
    want, _why = _canonical(found.group(1))
    want = want.replace("sumNsinα", "sumN'+uΔlsinα")
    drifted = full[0].replace(" + sum{kW}", "")
    if drifted == full[0]:
        fails.append("the mutation dropped nothing, so it tests nothing")
    elif _canonical(drifted)[0] == want:
        fails.append("a term dropped from the printed equation still matched "
                     "the page")
    return fails


#: The base normal, by the method that solves it alongside its quotient: the page
#: and the number the page publishes it under. Janbu's page names the denominator
#: m_α in its equation (1) and writes equation (6) over it; the report prints both
#: lines, so the two are joined before the comparison.
#: ``method -> (page, number, denominator)``, where ``denominator`` is the
#: equation that page names the base-normal denominator in, for a page that names
#: it rather than writing it out. Janbu's equation (1) introduces m_α beside the
#: base normal it divides; the report prints that definition on a line of its own
#: and it is held against the equation that introduces it.
_NORMAL_FORMS = {
    "bishop": ("lem/bishop.md", "8", ""),
    "janbu": ("lem/janbu.md", "6", "1"),
}

#: The strength mobilized on the slice base, as the two pages and the report each
#: write it. One product — c·Δl·sin α divided by F — with the fraction bar in
#: three places: Bishop's page divides c·Δl and multiplies by sin α, Janbu's
#: takes 1/F outside the product, and the report divides the whole product. They
#: are the same term, and where the bar is set is not what a transcription is
#: checked on; every other difference still is, because nothing but these three
#: spellings collapses.
_MOBILIZED_COHESION = re.compile(
    r"frac\{1\}\{F\}\s*c\s*Δl\s*sin\s*α"
    r"|frac\{c\s*Δl\}\{F\}\s*sin\s*α"
    r"|frac\{c·Δl·sin α\}\{F\}")


#: How a page spells a symbol the report prints in Unicode, for the one place a
#: symbol has to be FOUND in a page rather than compared with one.
_MATH_NAMES = {"m_α": r"m_\alpha"}


def _normal_canonical(text):
    """One base-normal equation reduced to what it says, with the mobilized
    cohesion written one way."""
    out, why = _page_latex(text)
    if why:
        return None, why
    out = _MOBILIZED_COHESION.sub("MOB", out)
    return "".join(c for c in out if c not in _GROUPING), ""


def test_the_normal_force_is_published_then_reduced():
    """Bishop's and Janbu's sections print the base normal their page publishes,
    then this model's reduction of it.

    Both cite their derivation's own number for it — Bishop's (8), Janbu's (6) —
    and both used to print, under that number, an equation with the pore
    pressure, the reinforcement, the pile force and the line load already taken
    out. A reader following the reference met a shorter equation than the one it
    named. So the base normal is held to the same discipline as the equation the
    method solves: the page's form first, carrying every vertical force a slice
    can take, then the sentence saying what this model does not carry, then what
    is left.

    Three things are required. The published form is the page's, symbol for
    symbol. The reduction really is a reduction — every term it drops is one the
    published form had. And the sentence between them accounts for every one.
    """
    fails = []
    from xslope.report import FORCE_TERMS, NotApplicable

    for method, (path, number, denominator) in sorted(_NORMAL_FORMS.items()):
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the base normal of")
            continue
        sentence, full, reduced = _normal_force_pair(section)
        if not full:
            fails.append(f"{method}: the section prints no base normal above a "
                         f"sentence reducing it: "
                         f"{[b.notation for b in section.blocks if b.kind == 'math']}")
            continue
        if not reduced:
            fails.append(f"{method}: the base normal is transcribed and never "
                         f"reduced to this model")
            continue

        # --- the published form is the page's ---
        with open(os.path.join(_REPO, "docs", path), encoding="utf-8") as f:
            page = f.read()
        # Janbu's m_α is a line of its own here and a numbered equation there;
        # the comparison is against the equation the number names, so the
        # denominator is substituted back before the two are held together.
        def against(source, printed, named):
            """The printed equation reduced to what it says, against the page's."""
            want, why = _normal_canonical(source)
            got, why_got = _normal_canonical(printed)
            if why or why_got:
                return [f"{method}: equation ({named}) is written with "
                        f"{why or why_got}"]
            if got != want:
                return [f"{method}: the section prints equation ({named}) as "
                        f"{got!r}; {path} publishes {want!r}"]
            return []

        found = re.search(r"\$([^$]*?)\\q?quad ?\(%s\)\$" % re.escape(number),
                          page)
        if not found:
            fails.append(f"{method}: {path} publishes no equation ({number})")
        else:
            fails += against(found.group(1), full[0], number)
        # A page that NAMES the denominator rather than writing it out: the
        # report prints the definition on a line of its own, and the equation
        # that introduces it on the page has to carry that definition.
        if denominator:
            if len(full) < 2:
                fails.append(f"{method}: {path} writes the base normal over a "
                             f"named denominator and the section never gives it")
            else:
                # The definition itself, out of the equation that introduces it:
                # that equation carries the base normal as well, and what the
                # report prints on this line is the denominator alone.
                name = full[1].split(" = ", 1)[0]
                found = re.search(
                    r"(%s\s*=\s*.*?)\s*\\q?quad ?\(%s\)"
                    % (re.escape(_MATH_NAMES.get(name, name)),
                       re.escape(denominator)), page)
                if not found:
                    fails.append(f"{method}: equation ({denominator}) of {path} "
                                 f"gives no {name}")
                else:
                    fails += against(found.group(1), full[1], denominator)

        # --- the reduction drops only what the published form had, and the
        # sentence accounts for every one of them ---
        published_terms = _terms_printed(full[0].split("}{", 1)[0] + "}")
        reduced_terms = _terms_printed(reduced[0].split("}{", 1)[0] + "}")
        dropped = [s for s in published_terms if s not in reduced_terms]
        if not dropped:
            fails.append(f"{method}: the reduced base normal drops nothing and "
                         f"is printed under a sentence that says it does: "
                         f"{sentence!r}")
        owner = {}
        for term in FORCE_TERMS:
            got = term.normal
            published = got.published if isinstance(got, NotApplicable) else got
            for contribution in published:
                owner[contribution.symbol] = term
        for symbol in dropped:
            term = owner.get(symbol)
            if term is None:
                fails.append(f"{method}: the base normal drops {symbol!r}, which "
                             f"belongs to no force in the registry")
            elif term.feature not in sentence and symbol not in sentence:
                fails.append(f"{method}: the base normal drops {symbol!r} and "
                             f"the sentence accounts for neither it nor the "
                             f"{term.feature}: {sentence!r}")
        # And nothing the sentence calls absent survives into the reduced form.
        for term in FORCE_TERMS:
            if not term.feature or term.feature not in sentence:
                continue
            carried = [s for s, t in owner.items()
                       if t is term and s in reduced_terms]
            if carried:
                fails.append(f"{method}: the sentence says the model carries no "
                             f"{term.feature}, and the reduced base normal "
                             f"prints {carried}")

    # The mutation: the section that prints the reduced form under the page's
    # number and nothing else — the report the defect was found in. With the
    # published run taken away there is no pair, and the check has to say so.
    report, _bundle = _calc_report("bishop")
    section = _calc_section(report)
    _sentence, full, reduced = _normal_force_pair(section)
    if full == reduced:
        fails.append("Bishop's published and reduced base normals are the same "
                     "equation, so the pair tests nothing")
    if not _is_normal_force(reduced):
        fails.append(f"the reduced base normal is not read as one: {reduced}")
    return fails


#: The equations a method's section prints from its OWN documentation page, and
#: the number each is published under there.
#:
#: Morgenstern-Price is the method this exists for. Its march is the
#: force-equilibrium page's, pinned against that page by :data:`_FULL_FORMS`;
#: what defines the METHOD is published on its own page and printed here: the
#: interslice assumption, the force function the solution was solved with, and
#: the moment about the coordinate origin the second condition sums over.
_OWN_PAGE_FORMS = {
    "mprice": ("lem/mprice.md", (
        ("2", "tan θ_j = λ·f(x_j)"),
        ("4", "f(x_j) = sin(π·frac{x_j − x_L}{x_R − x_L})"),
        ("8", "M_O = x_F·F_y − y_F·F_x"),
    )),
}


def _prints_its_page(method, path, page, number, notation, printed):
    """``[]`` where the section prints that equation and it is the page's, or the
    failures saying which of the two it is not."""
    if notation not in printed:
        return [f"{method}: {path} publishes equation ({number}) and the "
                f"section prints no {notation!r}: {printed}"]
    found = re.search(r"\$([^$]*?)\\q?quad ?\(%s\)\$" % re.escape(number), page)
    if not found:
        return [f"{method}: {path} publishes no equation ({number})"]
    want, why = _canonical(found.group(1))
    if why:
        return [f"{method}: equation ({number}) of {path} is written with {why}"]
    got, why = _canonical(notation)
    if why:
        return [f"{method}: the section prints {notation!r}, which carries {why}"]
    if got != want:
        return [f"{method}: the section prints equation ({number}) as {got!r}; "
                f"{path} publishes {want!r}"]
    return []


def test_the_method_prints_its_own_pages_equations():
    """Every equation a section takes from its own method's page is that page's,
    symbol for symbol.

    A method whose derivation publishes the assumption that DEFINES it has to
    print that assumption, not describe it: Morgenstern-Price's interslice
    inclination tan θ_j = λ·f(x_j), the f it was solved with, and the moment
    about the coordinate origin its second condition sums, each held against the
    equation its page numbers it.
    """
    fails = []
    pages = {}
    for method, (path, wanted) in _OWN_PAGE_FORMS.items():
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the equations of")
            continue
        printed = [b.notation for b in section.blocks if b.kind == "math"]
        with open(os.path.join(_REPO, "docs", path), encoding="utf-8") as f:
            pages[method] = page = f.read()
        for number, notation in wanted:
            fails += _prints_its_page(method, path, page, number, notation,
                                      printed)

    # The mutations. A section that stopped printing the assumption, and a page
    # whose display has drifted from what the section prints, both have to be
    # caught — the first is what the check is for, the second is what makes it a
    # pin on the page and not on the string written here.
    path, wanted = _OWN_PAGE_FORMS["mprice"]
    page = pages.get("mprice")
    if page is None:
        return fails + ["Morgenstern-Price produced no section to mutate"]
    number, notation = wanted[0]
    report, _bundle = _calc_report("mprice")
    printed = [b.notation for b in _calc_section(report).blocks
               if b.kind == "math"]
    dropped = [n for n in printed if n != notation]
    if dropped == printed:
        fails.append("the mutation dropped nothing, so it tests nothing")
    elif not _prints_its_page("mprice", path, page, number, notation, dropped):
        fails.append("a section that prints no interslice assumption still "
                     "passed")
    drifted = page.replace(r"\tan \theta_j = \lambda f(x_j)",
                           r"\tan \theta_j = \lambda")
    if drifted == page:
        fails.append("the mutation changed no equation on the page, so it tests "
                     "nothing")
    elif not _prints_its_page("mprice", path, drifted, number, notation,
                              printed):
        fails.append("a printed equation that drifted from its page still "
                     "matched it")
    return fails


def _compiled_scripts(notation):
    """Every subscript and superscript the compiled equation really carries, as
    ``(mark, base, script)`` — read back out of the Word math, not out of the
    notation it was written in."""
    from docx.oxml.ns import qn
    from xslope.report_docx import omath

    def text(node):
        return "".join(t.text or "" for t in node.iter(qn("m:t")))

    out = []
    math = omath(notation)
    for tag, mark in ((qn("m:sSub"), "_"), (qn("m:sSup"), "^")):
        for node in math.iter(tag):
            script = node.find(qn("m:sub") if mark == "_" else qn("m:sup"))
            out.append((mark, text(node.find(qn("m:e"))),
                        text(script) if script is not None else ""))
    return out


#: A symbol written with a subscript or superscript of more than one character,
#: which is what the compiler used to cut short.
_SCRIPTED = re.compile(r"([A-Za-zΑ-Ωα-ω])([_^])([A-Za-z0-9]{1,})")


def _scripts_are_whole(notation, where):
    """Every scripted symbol in one equation, compiled and read back."""
    compiled = _compiled_scripts(notation)
    fails = []
    for base, mark, script in _SCRIPTED.findall(notation):
        if (mark, base, script) not in compiled:
            fails.append(f"{where}: {base}{mark}{script} in {notation!r} "
                         f"compiles to {compiled}, and not to {script!r} under "
                         f"{base!r}")
    return fails


def test_scripts_are_not_cut_short():
    """A subscript runs to the end of the word it opens.

    ``F_corr`` is the corrected factor of safety, one symbol. It reached the page
    as F subscript c followed by the roman letters orr — "F c orr" — and so did
    every moment arm of the general equations: a_dx, a_ry, a_ey. The compiler
    took one character after the mark and set the rest as text beside it.

    Every equation every method prints is compiled here and read back out of the
    Word math, so the rule is checked where it is used and not only on the case
    that was reported.
    """
    fails = []
    import xslope.report_docx as report_docx

    # The pins: the symbol that was reported, and a moment arm, whose subscripts
    # have to arrive whole and under the right letter.
    for notation, want in (
            ("F_corr = f_o·F = 1.07702·1.634 = 1.760",
             [("_", "F", "corr"), ("_", "f", "o")]),
            ("sum{D cos β·a_dx} − sum{(H cos θ_p·a_ey + H sin θ_p·a_ex)}",
             [("_", "a", "dx"), ("_", "θ", "p"), ("_", "a", "ey"),
              ("_", "a", "ex")])):
        compiled = _compiled_scripts(notation)
        for one in want:
            if one not in compiled:
                fails.append(f"{notation!r} compiles to {compiled}, without "
                             f"{one}")

    # And the sweep: every equation of every method's calculation.
    swept = 0
    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            continue
        for block in section.blocks:
            if block.kind == "math":
                swept += 1
                fails += _scripts_are_whole(block.notation, method)
    if swept < len(CALC_METHODS):
        fails.append(f"only {swept} equations were compiled; the sweep is not "
                     f"reaching the sections")

    # Mutation: put the one-character rule back and the sweep has to go red.
    saved = report_docx._script_span
    report_docx._script_span = lambda src, at: at + 1
    try:
        caught = _scripts_are_whole("F_corr = f_o·F", "the mutation")
    finally:
        report_docx._script_span = saved
    if not caught:
        fails.append("a subscript cut to one character was not caught, so "
                     "nothing here tests the rule")
    return fails


#: What a printed equation carries besides symbols: the operators, the fence
#: characters, and the three functions and two macros the notation is written
#: with. Everything else that is a letter has to be a symbol the section defines.
_MATH_WORDS = ("frac", "sum", "sin", "cos", "tan", "π")
_MATH_PUNCTUATION = "·−-+=/{}[]()., '"


def test_printed_symbols_resolve():
    """Every symbol a printed equation carries is defined where it is printed.

    An equation prints a letter and the nomenclature under it says what the
    letter means; a letter that reaches the page with no row is a term the reader
    cannot look up. Both models here print one that used to: T, the tension-crack
    water force, which only Spencer's letter V was defined for, and P on a model
    whose reinforcement is axial, where the tangent column is zero on every slice
    and is not printed.
    """
    import re
    fails = []

    for where, xlsx in (("the tension-crack model", TENSION_XLSX),
                        ("the axially reinforced model", AXIAL_XLSX)):
        for method in CALC_METHODS:
            report, _bundle = _calc_report(method, xlsx=xlsx)
            if report is None:
                continue
            section = _calc_section(report)
            if section is None:
                continue
            defined = []
            for block in section.blocks:
                if block.kind == "table" and "Nomenclature" in block.caption:
                    defined += [row[0].split(" (")[0] for row in block.rows]
            for block in section.blocks:
                if block.kind != "math":
                    continue
                # The symbols first, longest first, and the numbers after them:
                # a subscript is a digit, and stripping the numbers ahead of the
                # symbols leaves R_ where R_1 was printed.
                text = block.notation
                for symbol in sorted(set(defined) | set(_MATH_WORDS),
                                     key=len, reverse=True):
                    text = text.replace(symbol, " ")
                text = re.sub(r"\d[\d,]*\.?\d*(?:e[+-]?\d+)?", " ", text)
                left = "".join(ch for ch in text
                               if ch not in _MATH_PUNCTUATION
                               and not ch.isspace())
                if left:
                    fails.append(f"{where} under {method}: {block.notation!r} "
                                 f"prints {left!r}, which the nomenclature "
                                 f"under it does not define")
    return fails


#: The symbols the report's equations use, and the LaTeX each is written as on
#: the documentation pages. The report prints Unicode and the pages print LaTeX,
#: so a notation check has to translate before it compares; where a page spells
#: one of them more than one way, every spelling it uses is listed.
DOC_SYMBOLS = {
    "α": (r"\alpha",),
    "β": (r"\beta",),
    "φ": (r"\phi",),
    "ψ": (r"\psi",),
    "θ": (r"\theta",),
    "Δl": (r"\Delta \ell", r"\Delta\ell"),
    "N'": ("N'",),
    "a_S": ("a_S",),
    "a_N": ("a_N",),
    "x_r": ("x_r",),
    "a_dx": ("a_{dx}",),
    "a_dy": ("a_{dy}",),
    "a_s": ("a_s",),
    "a_t": ("a_t",),
    "a_ry": ("a_{ry}",),
    "a_rx": ("a_{rx}",),
    "a_ey": ("a_{ey}",),
    "a_ex": ("a_{ex}",),
    "a_fx": ("a_{fx}",),
    "a_fy": ("a_{fy}",),
    "m_α": (r"m_\alpha", r"m_{\alpha}"),
    "θ_p": (r"\theta_p",),
    "y_Q": ("y_Q",),
    "x_b": ("x_b",),
    "R_1": ("R_1",),
    "R_2": ("R_2",),
    "Z_n": ("Z_n", "Z_{i+1}"),
    "M_O": ("M_O",),
    "F_x": ("F_x",),
    "F_y": ("F_y",),
    "x_F": ("x_F",),
    "y_F": ("y_F",),
    "θ_j": (r"\theta_j",),
    "λ": (r"\lambda",),
    "f(x_j)": ("f(x_j)",),
    "x_j": ("x_j",),
    "x_L": ("x_L",),
    "x_R": ("x_R",),
    "f_o": ("f_o",),
    "F_corr": ("F_{corr}",),
    "c_m": ("c_m",),
    "φ_m": (r"\phi_m",),
    "Z_i": ("Z_{i}", "Z_i"),
    "Z_{i+1}": ("Z_{i+1}",),
    "θ_i": (r"\theta_i",),
    "θ_{i+1}": (r"\theta_{i+1}",),
    "F_v": ("F_v",),
    "F_h": ("F_h",),
    "P_p": (r"P\cos \psi", r"P \cos \psi"),
    "H_p": (r"H \cos \theta_p",),
}


#: Constructions that describe the DOCUMENT instead of the analysis. The report
#: is engineering documentation: it states facts about the slope, the model and
#: the solution, and a sentence about the section, the table, the process of
#: presenting or what the reader should do with it is a defect (Norm). Each
#: phrase here was really in the prose and was cut; the list is what keeps the
#: voice from drifting back.
#:
#: "section" alone is not bannable — it is also the cross-section of the slope —
#: so only the unambiguous forms are listed.
#: The second group is about the PAGE — how much of it an equation takes, and
#: what was therefore done about it. "It runs several times the width of the
#: page, so it is written here as a quotient of its named sums" states a fact
#: about the typesetting and reaches the reader as an apology for the equation;
#: the equation's own decomposition is a fact about the equation and needs no
#: such excuse. Every phrase here was really in the prose and was cut.
_DOCUMENT_VOICE = (
    "a reader", "the reader", "a reviewer", "the engineer who",
    "this section", "the section prints", "the section says",
    "the section arrives", "the table lists", "the table below",
    "the figure below shows", "is worked through below",
    "the evaluations are not numbered",
    "width of the page", "on the page", "fits on", "too wide", "one line",
    "written here", "printed here", "shown here", "set out here",
    "for brevity", "abbreviated here",
)


def test_prose_is_about_the_analysis():
    """No paragraph of the report describes the report.

    The prose states facts about the section, the model and the solution. Self
    narration — what the section is doing, what a reader should re-form, why the
    document is arranged as it is — reads as an explanation of the work rather
    than the work, and every instance of it was cut. This is the guard.
    """
    fails = []

    reports = [("the default report", _build()),
               ("the seepage report", _engine_report("seep")),
               ("the strength reduction report", _engine_report("fem")),
               ("the reinforcement report",
                _engine_report("fem", xlsx=FEM_REINF_XLSX)),
               ("the piles report", _engine_report("fem", xlsx=FEM_PILES_XLSX))]
    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        if report is not None:
            reports.append((f"the {method} report", report))

    def narration(report):
        return [(b.text, phrase) for b in report.blocks("prose")
                for phrase in _DOCUMENT_VOICE if phrase in b.text.lower()]

    for where, report in reports:
        for text, phrase in narration(report):
            fails.append(f"{where} says {phrase!r}, which describes the "
                         f"document and not the analysis: {text!r}")

    # The mutation: the lead that explained how much of the page an equation
    # takes, and what was done about it, put back. It stood in Bishop's and the
    # Ordinary Method's sections and this check did not see it.
    import dataclasses
    from xslope.report import TRANSCRIPTIONS
    spec = TRANSCRIPTIONS["bishop"]
    key = ("bishop", (), REINF_XLSX)
    saved = _CALC.pop(key, None)
    TRANSCRIPTIONS["bishop"] = dataclasses.replace(spec, lead=(
        "Equation (10) of the derivation is the factor of safety this method "
        "solves. It runs several times the width of the page, so it is written "
        "here as a quotient of its named sums, one per force:"))
    try:
        narrated, _bundle = _calc_report("bishop")
        caught = narration(narrated) if narrated is not None else []
    finally:
        TRANSCRIPTIONS["bishop"] = spec
        _CALC.pop(key, None)
        if saved is not None:
            _CALC[key] = saved
    if not caught:
        fails.append("a lead narrating the width of the page and what was done "
                     "about it went through, so nothing here would stop that "
                     "phrasing coming back")
    return fails


def test_the_equation_is_cited_for_what_it_is():
    """No section cites a derivation for an equation that derivation does not
    publish.

    The force-equilibrium page derives a slice-by-slice march — its equations
    (6) and (10) — and Morgenstern-Price's derives no quotient at all. What
    Corps of Engineers, Lowe & Karafiath and Morgenstern-Price print is the
    horizontal force balance of the whole sliding mass at the converged
    solution, which is true because the march closed and is not what the page
    publishes. Those three say so, and still link the derivation for the march
    that reaches the solution.

    Janbu's page writes that balance directly as its equation (7), and Bishop's,
    the Ordinary Method of Slices' and Spencer's pages publish the equations
    their sections print, so all four cite the derivation for the equation.
    """
    fails = []
    from xslope.report import (WHOLE_MASS_BALANCE_METHODS, method_doc_url,
                               method_label)

    # Written out here rather than imported: a check that reads the same list
    # the prose is chosen from moves with it, and would pass on all seven
    # sections citing the derivation for an equation none of them prints.
    march = ("corps", "lowe", "mprice")
    if tuple(WHOLE_MASS_BALANCE_METHODS) != march:
        fails.append(f"{tuple(WHOLE_MASS_BALANCE_METHODS)} print the balance of "
                     f"the whole mass, and {march} is what publishes a march")

    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        section = _calc_section(report) if report is not None else None
        if section is None:
            fails.append(f"{method}: no calculation to read the citation of")
            continue
        intro = next((b for b in section.blocks if b.kind == "prose"), None)
        if intro is None:
            fails.append(f"{method}: the calculation opens with no statement of "
                         f"what its equation is")
            continue
        label = method_label(method)
        if (label, method_doc_url(method)) not in intro.links:
            fails.append(f"{method}: the opening statement does not link the "
                         f"published derivation")
        published = "the derivation published for" in intro.text
        if method in march:
            if "symbols of the derivation published" in intro.text:
                fails.append(f"{method}: cites the derivation for an equation "
                             f"that page does not publish: {intro.text!r}")
            if "march" not in intro.text:
                fails.append(f"{method}: the derivation is linked without "
                             f"saying what it derives: {intro.text!r}")
            fails += _quotient_is_introduced(section, method)
        elif not published:
            fails.append(f"{method}: prints the equation its page publishes and "
                         f"does not cite it: {intro.text!r}")
    return fails


#: What the sentence above the quotient has to establish, in the words it uses.
#: The reader's question is where the equation came from — the owner's, on the
#: report that printed it bare — and each of these is part of the answer: what
#: summing the march does, what is left, which numbered equation of which
#: derivation that is, and that it is not what the factor of safety was computed
#: from.
_QUOTIENT_LEAD = ("Summing the march's equation (6) over the slices",
                  "cancels the interslice forces",
                  "horizontal equilibrium of the whole sliding mass",
                  "equation (12) of the force-equilibrium derivation",
                  "not solved directly for F",
                  "holds at the converged factor of safety")


def _quotient_is_introduced(section, method):
    """The quotient a marching method prints is introduced where it is printed.

    Not at the head of the section: the reader meets the equation pages after the
    sentence that was supposed to have accounted for it, with the march and its
    reduction in between. The sentence has to be the block directly above the
    equation, it has to say what the equation is and which numbered equation of
    which derivation, and it has to link that derivation — the number is a
    reference, and a reference the reader cannot follow is not one.
    """
    from xslope.report import WHOLE_MASS_BALANCE_PAGE, docs_url

    fails = []
    blocks = list(section.blocks)
    at = next((i for i, b in enumerate(blocks)
               if b.kind == "math" and b.notation.startswith("F = frac")), None)
    if at is None:
        return [f"{method}: prints no quotient for the factor of safety"]
    lead = blocks[at - 1] if at and blocks[at - 1].kind == "prose" else None
    if lead is None:
        return [f"{method}: the quotient is printed with no sentence above it "
                f"saying what it is"]
    for phrase in _QUOTIENT_LEAD:
        if phrase not in lead.text:
            fails.append(f"{method}: the sentence above the quotient does not "
                         f"say {phrase!r}: {lead.text!r}")
    url = docs_url(WHOLE_MASS_BALANCE_PAGE)
    if not any(str(target).rstrip("/") == url.rstrip("/")
               for _text, target in (lead.links or ())):
        fails.append(f"{method}: the sentence names an equation of "
                     f"{WHOLE_MASS_BALANCE_PAGE} and does not link it: "
                     f"{lead.links}")
    return fails


def test_calculation_notation_matches_the_docs():
    """Every symbol the printed equations use is a symbol the method's own
    documentation page uses.

    The report and the derivation have to be readable side by side. This is a
    light guard — it compares symbols, not equations — but it is what catches a
    silent drift in notation, which is the way the two would come apart.

    A section is read against its own page AND against any page it transcribes an
    equation from: Morgenstern-Price's march is the force-equilibrium
    derivation's per-slice system, and its own page writes none of that system's
    symbols.
    """
    import re
    fails = []
    from xslope.report import METHOD_DOC_PAGES, MOMENT_PART_SYMBOLS

    for method in CALC_METHODS:
        report, _bundle = _calc_report(method)
        if report is None:
            continue
        section = _calc_section(report)
        if section is None:
            fails.append(f"{method}: no calculation to check the notation of")
            continue
        paths = [METHOD_DOC_PAGES[method]]
        paths += [p for p, _n in _EQUATION_NUMBERS.get(method, ())
                  if p not in paths]
        source = ""
        for path in paths:
            page = os.path.join(_REPO, "docs", path)
            if not os.path.exists(page):
                fails.append(f"{method}: the documentation page {page} is "
                             f"missing")
                continue
            with open(page, encoding="utf-8") as f:
                source += f.read()
        if not source:
            continue
        for block in section.blocks:
            if block.kind != "math":
                continue
            # Longest first: Δl must not be tested as Δ and then l.
            text = block.notation
            # The letters the two moment methods' quotient is named in parts by
            # are the report's own and appear on no page. What pins them is the
            # recomposition — substituted back they have to give the published
            # equation — so they are taken out before the sweep for symbols the
            # documentation does not carry.
            for part in sorted(MOMENT_PART_SYMBOLS, key=len, reverse=True):
                text = text.replace(part, " ")
            for symbol in sorted(DOC_SYMBOLS, key=len, reverse=True):
                if symbol not in text:
                    continue
                text = text.replace(symbol, " ")
                if not any(spelling in source
                           for spelling in DOC_SYMBOLS[symbol]):
                    fails.append(
                        f"{method}: the equations use {symbol!r} but "
                        f"{', '.join(paths)} never writes "
                        f"{DOC_SYMBOLS[symbol][0]!r}")
            # Nothing left over that looks like a symbol we have not declared.
            for stray in re.findall(r"[A-Za-zΑ-Ωα-ω]_\{?[A-Za-zΑ-Ωα-ω0-9]+", text):
                fails.append(f"{method}: {stray!r} in {block.notation!r} is not "
                             f"in the notation map, so nothing checks it "
                             f"against the documentation")
    return fails


#: A phrase from each method's summary that only that method's summary carries —
#: the assumption or the equilibrium condition that distinguishes it. Written out
#: here rather than derived, so that one paragraph pasted over another fails: the
#: check requires each phrase in its own method's block AND out of every other's.
_SUMMARY_MARKS = {
    "oms": ("about the center of the circle, and no other condition",
            "no iteration"),
    "bishop": ("vertical force equilibrium of each slice and moment equilibrium",
               "assumed horizontal"),
    "janbu": ("horizontal force equilibrium of the sliding mass as a whole",
              "correction factor f_o"),
    "corps": ("the ground surface at each slice boundary",),
    "lowe": ("average of the ground-surface slope",),
    "spencer": ("one constant inclination θ",),
    "mprice": ("tan θ = λ·f(x)",),
}


def test_method_summary_opens_each_block():
    """Every method's block opens by saying what the method assumes and which
    equilibrium conditions it satisfies.

    The block has to stand on its own: a reader who arrives at a factor of
    safety or a slice table must be able to tell what produced it without paging
    back to the Calculations section, which may be switched off. So the paragraph
    is required as the FIRST block of the method's own section — above its
    search, its results and its slice table alike — and required again with the
    calculations off.

    It also has to be about the method it opens. Each method carries a phrase
    that appears in its summary and in no other, which is what makes a
    copy-pasted paragraph a failure rather than a pass.
    """
    fails = []
    from xslope.report import METHOD_SUMMARIES, method_label, method_summary

    def method_section(report, method):
        want = method_label(method)
        for section in report.sections:
            for _lvl, node in section.walk():
                if node.title == want:
                    return node
        return None

    for method in CALC_METHODS:
        for label, options in (("calculations on", None),
                               ("calculations off", {"lem_calculations": False})):
            report, _bundle = _calc_report(method, options)
            if report is None:
                continue
            node = method_section(report, method)
            if node is None:
                fails.append(f"{method} ({label}): no {method_label(method)} "
                             f"section to open")
                continue
            if not node.blocks or node.blocks[0].kind != "prose":
                kind = node.blocks[0].kind if node.blocks else "nothing"
                fails.append(f"{method} ({label}): the section opens with "
                             f"{kind}, not the method summary")
                continue
            opening = node.blocks[0].text
            if opening != method_summary(method):
                fails.append(f"{method} ({label}): the opening paragraph is not "
                             f"the method summary: {opening[:80]!r}")
                continue
            if "equilibrium" not in opening:
                fails.append(f"{method} ({label}): the opening paragraph names "
                             f"no equilibrium condition")

    # Every method the report can document has a summary, and each summary is
    # about its own method: its marks are in it and out of all the others.
    for method in CALC_METHODS:
        if not method_summary(method):
            fails.append(f"{method}: no summary paragraph is written for it")
            continue
        for mark in _SUMMARY_MARKS[method]:
            if mark not in METHOD_SUMMARIES[method]:
                fails.append(f"{method}: its summary no longer says {mark!r}, so "
                             f"nothing distinguishes it from the others")
            for other in CALC_METHODS:
                if other != method and mark in METHOD_SUMMARIES.get(other, ""):
                    fails.append(f"{other}: its summary carries {mark!r}, which "
                                 f"belongs to {method} — the two paragraphs say "
                                 f"the same thing")
    return fails


def test_docs_links():
    """The published documentation is linked, at the URL the site really uses."""
    fails = []
    from xslope.report import (DOCS_BASE_URL, METHOD_DOC_PAGES, docs_url,
                               method_doc_url)

    mkdocs = os.path.join(_REPO, "mkdocs.yml")
    with open(mkdocs, encoding="utf-8") as f:
        for line in f:
            if line.startswith("site_url:"):
                declared = line.split(":", 1)[1].strip()
                if declared.rstrip("/") != DOCS_BASE_URL.rstrip("/"):
                    fails.append(f"the report links to {DOCS_BASE_URL}, the site "
                                 f"is published at {declared}")
                break
        else:
            fails.append("mkdocs.yml declares no site_url to check against")

    if docs_url("lem/oms.md") != DOCS_BASE_URL + "lem/oms/":
        fails.append(f"docs_url('lem/oms.md') = {docs_url('lem/oms.md')!r}")
    if docs_url("usage/claude/index.md") != DOCS_BASE_URL + "usage/claude/":
        fails.append(f"an index page maps to {docs_url('usage/claude/index.md')!r}")

    for method, page in METHOD_DOC_PAGES.items():
        if not os.path.exists(os.path.join(_REPO, "docs", page)):
            fails.append(f"{method} cites docs/{page}, which does not exist")
        if not method_doc_url(method).endswith("/"):
            fails.append(f"{method}'s URL is {method_doc_url(method)!r}")

    # The finite element members cite the pages their formulations come from.
    from xslope.report import FEM_DETAIL_DOC_PAGES
    for kind, page in FEM_DETAIL_DOC_PAGES.items():
        if not os.path.exists(os.path.join(_REPO, "docs", page)):
            fails.append(f"the {kind} section cites docs/{page}, which does "
                         f"not exist")
        if not docs_url(page).endswith("/"):
            fails.append(f"the {kind} page's URL is {docs_url(page)!r}")

    from xslope.report import supported_methods
    for method in supported_methods():
        if method not in METHOD_DOC_PAGES:
            fails.append(f"the solver offers {method!r} and no documentation "
                         f"page is mapped for it")
    return fails


def test_calculation_in_the_document():
    """The section reaches the .docx as Word math, with its links live.

    Three things a reader depends on: the equations are text (Word's own math,
    not pictures), the reference to the slice table jumps to the slice table, and
    the citation of the derivation opens the published page.
    """
    import re
    fails = []
    from xslope.report import SLICE_TABLE_BOOKMARK, method_doc_url

    def write(method, options=None):
        import matplotlib
        matplotlib.use("Agg")
        from xslope.fileio import load_slope_data
        from xslope.report import generate_report
        from xslope.slice import generate_slices
        from xslope.solve import solve_selected

        slope_data = load_slope_data(REINF_XLSX)
        ok, out = generate_slices(slope_data, circle=slope_data["circles"][0],
                                  num_slices=15)
        df, surface = out[0].copy(), out[1]
        with contextlib.redirect_stdout(io.StringIO()):
            results = solve_selected(method, df)
        opts = {"method": method, "title": "Calculation", "pd_figure": False,
                "lem_search_figure": False, "lem_solution_figure": False}
        opts.update(options or {})
        tmp = tempfile.mkdtemp(prefix="xslope_calc_")
        path = os.path.join(tmp, f"{method}.docx")
        ok, info = generate_report(
            slope_data,
            {"lem": [{"slice_df": df, "failure_surface": surface,
                      "results": results, "search": None, "method": method}]},
            opts, path)
        if not ok:
            return None, None
        with zipfile.ZipFile(path) as z:
            return (z.read("word/document.xml").decode(),
                    z.read("word/_rels/document.xml.rels").decode())

    doc, rels = write("bishop")
    if doc is None:
        return ["the report could not be written"]

    if "<m:oMath>" not in doc:
        fails.append("the equations are not Word math")
    if doc.count("<m:oMathPara>") < 3:
        fails.append(f"only {doc.count('<m:oMathPara>')} displayed equations "
                     f"reached the document")
    for wanted in ("<m:f>", "<m:nary>", "<m:sSub>"):
        if wanted not in doc:
            fails.append(f"no {wanted} in the document: the notation is not "
                         f"compiling to real math")
    if "m:val=\"∑\"" not in doc:
        fails.append("no summation sign in the math")
    # The numbers are in the math, not in a picture of it.
    if not re.search(r"<m:t[^>]*>[^<]*1\.9", doc):
        fails.append("the factor of safety is not text inside the equation")

    # One bookmark per method's slice table: a report that documents several
    # carries several, and each calculation links to its own.
    mark = f"{SLICE_TABLE_BOOKMARK}_bishop"
    if f'w:name="{mark}"' not in doc:
        fails.append("the slice table carries no bookmark to link to")
    if f'w:anchor="{mark}"' not in doc:
        fails.append("nothing links to the slice table's bookmark")
    else:
        order = (doc.index(f'w:anchor="{mark}"'), doc.index(f'w:name="{mark}"'))
        if order[0] < order[1]:
            fails.append("the cross-reference is written before its bookmark")

    external = re.findall(r'Target="([^"]+)"[^>]*TargetMode="External"', rels)
    if method_doc_url("bishop") not in external:
        fails.append(f"the document links to {external}, not to Bishop's page "
                     f"{method_doc_url('bishop')}")

    # And the citation follows the method the report documents.
    doc2, rels2 = write("spencer")
    external2 = re.findall(r'Target="([^"]+)"[^>]*TargetMode="External"', rels2)
    if method_doc_url("spencer") not in external2:
        fails.append(f"a Spencer report links to {external2}")
    if method_doc_url("bishop") in external2:
        fails.append("a Spencer report still links to Bishop's page")

    # Switched off, the section and its equations are gone.
    doc3, _rels3 = write("bishop", {"lem_calculations": False})
    if "Calculations" in doc3:
        fails.append("lem_calculations=False left the section in the document")
    if "<m:oMath>" in doc3:
        fails.append("lem_calculations=False left equations in the document")
    if f'w:anchor="{mark}"' in doc3:
        fails.append("lem_calculations=False left a link to the slice table")
    return fails


#: The slice force diagram each method's Calculations section opens on, written
#: here from the documentation rather than imported from the report: it is the
#: figure the method's own derivation page displays as the complete set of forces
#: on a slice. Janbu's derivation takes the Ordinary Method's figure, the two
#: force-equilibrium methods share the one on the page they share, and
#: Morgenstern-Price takes Spencer's, which its page displays. A mapping that
#: quietly changed — Bishop's section headed by the force-equilibrium slice —
#: would print an equation over a picture of different forces, and only a
#: second copy of the mapping can catch that.
EXPECTED_DIAGRAMS = {
    "oms": "oms_complete.png",
    "janbu": "oms_complete.png",
    "bishop": "bishop_complete.png",
    "corps": "slice_fe_complete.png",
    "lowe": "slice_fe_complete.png",
    "spencer": "spencer3_forces.png",
    "mprice": "spencer3_forces.png",
}


def _calculations_section(block):
    """The Calculations section of a method's block, or None."""
    if block is None:
        return None
    return next((node for _lvl, node in block.walk()
                 if node.title == "Calculations"), None)


def _smallest_label_px(path):
    """The height, in pixels, of the small lettering on a drawing.

    Every dark blob between four and forty pixels tall and three and forty wide
    is a glyph or part of one; the arrowheads, the slice outline and the
    dimension lines are outside that band. The tenth percentile of what is left
    is the small lettering — the subscripts, which are what binds the printed
    size — rather than a stray mark, which one blob at the bottom of the list
    would be.
    """
    import numpy as np
    from PIL import Image
    from scipy import ndimage

    grey = np.array(Image.open(path).convert("L"))
    labels, _n = ndimage.label(grey < 100)
    heights = []
    for rows, cols in ndimage.find_objects(labels):
        h, w = rows.stop - rows.start, cols.stop - cols.start
        if 4 <= h <= 40 and 3 <= w <= 40:
            heights.append(h)
    if not heights:
        return None, grey.shape[1]
    return float(np.percentile(sorted(heights), 10)), grey.shape[1]


def test_force_diagram_is_as_small_as_it_reads():
    """Every force diagram prints at the width its own lettering earns, and no
    wider.

    One width for all four printed the Ordinary Method's eight large labels half
    again the size of the body text and cost every section the page area
    Spencer's thirty small ones need. So each drawing is scaled until its
    SMALLEST label — its subscripts, which bind — reaches the height at which
    those labels are read on the page, floored at the narrowest a free body can
    be drawn and still be taken from.

    The two rules are stated in the module the report reads them from, and the
    widths are re-derived here from the PNGs themselves: a redrawn diagram cannot
    keep a width its new lettering no longer earns.
    """
    fails = []
    from xslope.report import (FORCE_DIAGRAM_LABEL_IN, FORCE_DIAGRAM_MIN_IN,
                               FORCE_DIAGRAM_WIDTHS, FORCE_DIAGRAMS,
                               force_diagram, force_diagram_width)

    drawn = set(FORCE_DIAGRAMS.values())
    if set(FORCE_DIAGRAM_WIDTHS) != drawn:
        fails.append(f"the widths are stated for {sorted(FORCE_DIAGRAM_WIDTHS)} "
                     f"and the sections print {sorted(drawn)}")
    for method, name in sorted(FORCE_DIAGRAMS.items()):
        path = force_diagram(method)
        if not path:
            fails.append(f"{method}: the diagram {name} is not in the package")
            continue
        smallest, pixels = _smallest_label_px(path)
        if smallest is None:
            fails.append(f"{name}: no lettering could be measured on it")
            continue
        earned = FORCE_DIAGRAM_LABEL_IN * pixels / smallest
        # Up to the twentieth of an inch, which is the grain the widths are set
        # on. Up, not down: a width rounded down prints the lettering under the
        # height the width was derived from.
        want = max(FORCE_DIAGRAM_MIN_IN, math.ceil(earned * 20) / 20)
        got = force_diagram_width(method)
        if abs(got - want) > 1e-9:
            fails.append(f"{name}: prints {got} in wide; its smallest label is "
                         f"{smallest:g} px of {pixels}, which earns {want} in "
                         f"({FORCE_DIAGRAM_LABEL_IN} in per label, floored at "
                         f"{FORCE_DIAGRAM_MIN_IN} in)")
        # And the rule itself holds on the printed page: no label smaller than
        # the floor the widths were derived from.
        printed = smallest * got / pixels
        if printed < FORCE_DIAGRAM_LABEL_IN - 1e-9:
            fails.append(f"{name}: at {got} in its smallest label prints "
                         f"{printed:.4f} in, under the {FORCE_DIAGRAM_LABEL_IN} "
                         f"in it is legible at")

    # The widths are not all the same: one width for four drawings that carry
    # three times the lettering of each other is what this replaced.
    if len(set(FORCE_DIAGRAM_WIDTHS.values())) < 2:
        fails.append(f"every diagram prints at one width again: "
                     f"{FORCE_DIAGRAM_WIDTHS}")

    # And no diagram is wider than the one width they all used to share, which is
    # the complaint this answers: the figures took more of the page than they
    # needed.
    if any(w > 3.25 for w in FORCE_DIAGRAM_WIDTHS.values()):
        fails.append(f"a diagram is wider than the 3.25 in they all shared: "
                     f"{FORCE_DIAGRAM_WIDTHS}")
    return fails


def _diagram_figures(report):
    """Every force-diagram figure in a report, by caption."""
    return [f for f in report.figures() if f.caption.startswith("Forces on a slice")]


def test_force_diagram_heads_the_calculations():
    """Every Calculations section opens on the free body its equations are
    written about.

    The equations print in the symbols of the derivation, and the derivation
    draws those symbols on a slice. Printed without it, the letters meet the
    reader for the first time in a quotient. The figure is therefore the FIRST
    block of the section — before the sentence that introduces the equation —
    and it is the diagram that method's own documentation page displays.

    It prints exactly where the working prints: with the calculations switched
    off there is no section and no diagram, and a model whose equilibrium cannot
    be worked through gets the sentence that says so and no picture of a
    quotient that was never printed.
    """
    fails = []
    from xslope.report import METHOD_DOC_PAGES, method_label

    for method in CALC_METHODS:
        built, _bundle = _calc_report(method)
        if built is None:
            fails.append(f"the sample model did not solve with {method}")
            continue
        label = method_label(method)
        sec = _calculations_section(_method_block(built, method))
        if sec is None:
            fails.append(f"{method}: the block carries no Calculations section")
            continue
        if not sec.blocks or sec.blocks[0].kind != "figure":
            fails.append(f"{method}: the Calculations section opens on "
                         f"{[b.kind for b in sec.blocks[:3]]}, not its force "
                         f"diagram")
            continue
        fig = sec.blocks[0]
        if fig.caption != f"Forces on a slice — {label}":
            fails.append(f"{method}: the diagram is captioned {fig.caption!r}")
        if fig.number < 1:
            fails.append(f"{method}: the diagram carries figure number "
                         f"{fig.number}")
        if fig.landscape:
            fails.append(f"{method}: the diagram asks for a landscape page")
        if not 2.0 <= fig.width_in <= 6.5:
            fails.append(f"{method}: the diagram prints {fig.width_in} in wide")

        # And the section cites it, by the number the counter gave it — a
        # technical report refers to every figure it prints. The citation belongs
        # in the sentence that first sums the forces the diagram draws, which is
        # where a reader is sent to the picture. Spencer's is the section that
        # carries one today; the other six are their own round.
        if method == "spencer":
            import re
            says = next((b.text for b in sec.blocks
                         if b.kind == "prose" and "F_h and F_v" in b.text), "")
            cited = set(re.findall(r"Figure (\d+)", says))
            if str(fig.number) not in cited:
                fails.append(f"{method}: the diagram is Figure {fig.number} and "
                             f"the sentence that sums its forces cites "
                             f"{sorted(cited) or 'no figure'}: {says!r}")
            if cited - {str(fig.number)}:
                fails.append(f"{method}: the sentence that sums the forces of "
                             f"Figure {fig.number} cites "
                             f"{sorted(cited - {str(fig.number)})} as well")

        want = os.path.join(_REPO, "xslope", "resources",
                            EXPECTED_DIAGRAMS[method])
        if not os.path.exists(fig.path):
            fails.append(f"{method}: the diagram file {fig.path} is not there")
        elif open(fig.path, "rb").read() != open(want, "rb").read():
            fails.append(f"{method}: the section prints "
                         f"{os.path.basename(fig.path)}; the derivation draws "
                         f"{EXPECTED_DIAGRAMS[method]}")

        # The mapping's premise: the diagram is what the method's own published
        # derivation displays. A redraw under a new name would leave the report
        # printing a figure the documentation no longer shows.
        page = os.path.join(_REPO, "docs", METHOD_DOC_PAGES[method])
        if os.path.exists(page):
            with open(page, encoding="utf-8") as fh:
                text = fh.read()
            if EXPECTED_DIAGRAMS[method] not in text:
                fails.append(f"{method}: {METHOD_DOC_PAGES[method]} does not "
                             f"display {EXPECTED_DIAGRAMS[method]}")

        # One diagram per method block, at the head of the calculations and
        # nowhere else.
        found = _diagram_figures(built)
        if len(found) != 1:
            fails.append(f"{method}: the report carries {len(found)} force "
                         f"diagrams for one method")

    # --- switched off: no section, so no diagram ---
    off, _bundle = _calc_report("spencer", options={"lem_calculations": False})
    if off is not None and _diagram_figures(off):
        fails.append("with the calculations switched off the report still "
                     "carries a force diagram")

    # --- a model whose equilibrium cannot be worked through ---
    #
    # Passive support divides by the factor of safety on the driving side of the
    # force methods, so those blocks print the sentence instead of a quotient.
    # The moment methods can show it, and their diagrams stand.
    from xslope.fileio import load_slope_data
    from xslope.report import planned_figures, resolve_options

    passive = load_slope_data(PASSIVE_XLSX)
    for method, wanted in (("janbu", 0), ("corps", 0), ("lowe", 0),
                           ("mprice", 0), ("spencer", 0),
                           ("bishop", 1), ("oms", 1)):
        built, bundle = _calc_report(method, xlsx=PASSIVE_XLSX)
        if built is None:
            fails.append(f"the passive model did not solve with {method}")
            continue
        got = len(_diagram_figures(built))
        if got != wanted:
            fails.append(f"{method}: the passive model prints {got} force "
                         f"diagrams and {wanted} calculations")
        # The count a caller is shown follows the same refusal: a diagram that
        # is not printed is not planned either.
        opts = resolve_options({"input_path": PASSIVE_XLSX, "method": method,
                                "pd_figure": False, "lem_search_figure": False,
                                "lem_solution_figure": False})
        with contextlib.redirect_stdout(io.StringIO()):
            planned = planned_figures(passive, {"lem": [bundle]}, opts)
        drawn = len(built.figures())
        if planned != drawn:
            fails.append(f"{method}: {planned} figures were planned for the "
                         f"passive model and {drawn} were built")
    return fails


# --------------------------------------------------------------------------
# F. the shared plot
# --------------------------------------------------------------------------

def test_profile_lines_name_their_materials():
    """Every profile line's legend entry names the material it bounds.

    The model figure draws numbered profile lines and the materials table names
    materials; with nothing tying the two together a reader has no way to tell
    which line is the core and which the shell. The name goes in the legend
    entry — "Profile 2 (Core)" — and it is the same base geometry every mode of
    the plot draws, so every mode carries it.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MplFigure
    from xslope.plot import plot_inputs

    dam = load_slope_data_cached(SEEP_XLSX)
    lines = dam.get("profile_lines") or []
    mats = dam.get("materials") or []
    if len(lines) < 2 or not mats:
        return ["the dam carries no layered profile lines; the legend is untested"]

    # What each line SHOULD say, read off the model rather than spelled here.
    want = []
    for i, line in enumerate(lines):
        idx = line.get("mat_id")
        idx = idx if isinstance(idx, int) and 0 <= idx < len(mats) else i
        name = str((mats[idx] or {}).get("name") or "").strip()
        want.append(f"Profile {i + 1}" + (f" ({name})" if name else ""))
    if len({w for w in want}) < 2:
        fails.append(f"every profile line would carry the same legend: {want}")

    for mode in ("shared", "lem", "seep", "fem"):
        fig = MplFigure(figsize=(9, 5.5))
        FigureCanvasAgg(fig)
        with contextlib.redirect_stdout(io.StringIO()):
            plot_inputs(dam, fig=fig, mode=mode, show_title=False,
                        frame="content")
        labels = fig.axes[0].get_legend_handles_labels()[1]
        for entry in want:
            if entry not in labels:
                fails.append(f"mode={mode!r}: the legend has no {entry!r} "
                             f"({[l for l in labels if l.startswith('Profile')]})")
    return fails


def test_shared_plot():
    """mode="shared" draws the model without the trial surfaces, and draws the
    water line the seepage head boundaries state."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from xslope.fileio import load_slope_data
    from xslope.plot import plot_derived_water_lines, plot_inputs
    from xslope.water import water_line_for_stage

    def draw(slope_data, mode):
        fig = Figure(figsize=(9, 5.5))
        FigureCanvasAgg(fig)
        plot_inputs(slope_data, fig=fig, mode=mode, show_title=False,
                    frame="content")
        return fig.axes[0]

    reinf = load_slope_data(REINF_XLSX)
    lem_ax = draw(reinf, "lem")
    shared_ax = draw(reinf, "shared")
    lem_labels = set(lem_ax.get_legend_handles_labels()[1])
    shared_labels = set(shared_ax.get_legend_handles_labels()[1])
    if not any("Circle" in l or "Surface" in l for l in lem_labels):
        fails.append(f"the LEM plot draws no trial circles: {sorted(lem_labels)}")
    for suppressed in ("Starting Circle", "Non-Circular Surface"):
        if suppressed in shared_labels:
            fails.append(f"the shared plot still draws {suppressed!r}")
    for shared_thing in ("Distributed Load", "Reinforcement Line"):
        if shared_thing not in shared_labels:
            fails.append(f"the shared plot dropped {shared_thing!r}")

    # A model with a reservoir/head boundary: the water line is on the figure,
    # and it is the same line the load derivation measures against.
    dam = load_slope_data(DAM_XLSX)
    derived = water_line_for_stage(dam, stage=1)
    if not derived.get("points"):
        return fails + ["the dam sample no longer states a pool; the water-line "
                        "check has nothing to see"]
    ax = draw(dam, "shared")
    lines = [ln for ln in ax.lines if ln.get_gid() == "WATER_LINE"]
    if not lines:
        fails.append("the shared plot of a model with a reservoir boundary draws "
                     "no water line")
    else:
        # The drawn pool stands at the level the boundary states — read from the
        # boundary itself, not from the drawn line, so this is an oracle and not
        # a restatement. (The derived line lies ON the ground where it is dry, so
        # its own maximum is a hilltop, not a water level.)
        from xslope.water import bc_pool_levels
        levels = bc_pool_levels(dam["seepage_bc"], dam["ground_surface"], dam)
        want = max(level for level, _anchors, _label in levels)
        ys = [y for ln in lines for y in ln.get_ydata()]
        if abs(max(ys) - want) > 1e-6:
            fails.append(f"the drawn water line tops out at {max(ys)}, the "
                         f"boundary's pool at {want}")
        # The pool ends at its shoreline: at the drawn line's downstream end the
        # ground stands at the water level, so the load that follows it is not
        # run up a dry face.
        from xslope.water import _y_on
        x_end = max(x for ln in lines for x in ln.get_xdata())
        y_ground = _y_on(dam["ground_surface"], float(x_end))
        if y_ground is None or abs(y_ground - want) > 1e-6:
            fails.append(f"the water line ends at x = {x_end} where the ground is "
                         f"at {y_ground}, not at the {want} shoreline")
    if "Water Surface" not in set(ax.get_legend_handles_labels()[1]):
        fails.append("the water line has no legend key")

    # The layer itself, called directly: a dry model draws nothing.
    fig = Figure(figsize=(6, 4))
    FigureCanvasAgg(fig)
    if plot_derived_water_lines(fig.add_subplot(111), reinf):
        fails.append("a model with no head boundaries drew a water line")

    # The mesh belongs to the mesh figures, not to the shared model — and an
    # engine's model figure can ask for the same, which is what the report's
    # seepage and finite element inputs do: each gives the mesh a figure of its
    # own immediately below, and drawn twice the underlay is a grid over the
    # zones the model figure is there to show.
    if dam.get("mesh") is not None:
        from matplotlib.collections import LineCollection

        def meshed(mode, **kw):
            fig = Figure(figsize=(9, 5.5))
            FigureCanvasAgg(fig)
            plot_inputs(dam, fig=fig, mode=mode, show_title=False,
                        frame="content", **kw)
            return any(isinstance(c, LineCollection)
                       for c in fig.axes[0].collections)

        if meshed("shared"):
            fails.append("the shared plot drew the analysis mesh")
        if not meshed("fem"):
            fails.append("an engine view drew no mesh even when asked for one; "
                         "the suppression below proves nothing")
        for mode in ("fem", "seep"):
            if meshed(mode, show_mesh=False):
                fails.append(f"mode={mode!r} with show_mesh=False still drew the "
                             f"analysis mesh")
    return fails


def test_model_figure_coordinate_labels():
    """The model figure names its points, and the toggle really reaches the plot.

    Geometry is one of the things the model figure is there to show, so the
    coordinate labels are on by default and the dialog's checkbox is for the
    models they crowd. A checkbox wired to nothing looks exactly like one that
    works, so this follows the option from the report's own builder to the
    keyword ``plot_inputs`` is called with, and then to the artists that keyword
    produces: with it on, a vertex the model actually has is annotated with its
    own coordinates; with it off, not one coordinate label is drawn and the
    figure — which stays in the report either way — is a different PNG.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from xslope import plot as plot_mod
    from xslope.report import DEFAULT_OPTIONS, build_report

    slope_data, solutions = _solved()
    if DEFAULT_OPTIONS.get("pd_coords") is not True:
        fails.append("the point coordinates are not on by default")

    # A vertex the model really has, written the way the annotator writes it.
    lines = slope_data.get("profile_lines") or []
    if not lines:
        return fails + ["the sample model has no profile lines; the coordinate "
                        "labels have nothing to name"]
    vx, vy = lines[0]["coords"][0]
    wanted = f"({float(vx):g}, {float(vy):g})"

    passed, drawn, pngs = [], [], []
    real = plot_mod.plot_inputs

    def spy(sd, *args, **kw):
        passed.append(kw.get("label_coordinates", "<not passed>"))
        out = real(sd, *args, **kw)
        fig = kw.get("fig")
        drawn.append([t.get_text() for ax in (fig.axes if fig else [])
                      for t in ax.texts if t.get_gid() == "COORD_LABEL"])
        return out

    base = {"traceability": False, "lem": False, "model_checks": False}
    plot_mod.plot_inputs = spy
    try:
        for state in (True, False):
            with tempfile.TemporaryDirectory() as tmp:
                report = build_report(slope_data, solutions,
                                      dict(base, pd_coords=state), tmp)
                figs = report.figures()
                if len(figs) != 1:
                    fails.append(f"pd_coords={state} left the project definition "
                                 f"with {len(figs)} figures, not one — the labels "
                                 f"option moved the figure itself")
                    pngs.append(None)
                else:
                    with open(figs[0].path, "rb") as fh:
                        pngs.append(fh.read())
                said = " ".join(b.text for b in report.blocks("prose"))
                names_them = "labeled with its coordinates" in said
                if names_them is not state:
                    fails.append(f"pd_coords={state} and the prose "
                                 f"{'announces' if names_them else 'does not announce'} "
                                 f"the coordinate labels")

        # And with no figure at all the plot is not called, so nothing is added
        # to `passed` and the count below still reads [True, False].
        with tempfile.TemporaryDirectory() as tmp:
            build_report(slope_data, solutions,
                         dict(base, pd_figure=False, pd_coords=True), tmp)
    finally:
        plot_mod.plot_inputs = real

    if passed != [True, False]:
        fails.append(f"the plot was asked for label_coordinates={passed}, not "
                     f"[True, False] — the toggle does not reach the figure")
    if len(drawn) >= 2:
        if wanted not in drawn[0]:
            fails.append(f"the labeled figure does not name the vertex "
                         f"{wanted}; it drew {len(drawn[0])} coordinate labels")
        if drawn[1]:
            fails.append(f"the unlabeled figure still drew {len(drawn[1])} "
                         f"coordinate labels")
    if None not in pngs and pngs[0] == pngs[1]:
        fails.append("the figure is byte-identical with the labels on and off")
    return fails


def _coord_labels(xlsx, figsize=(11.0, 5.0)):
    """``(axes, [(text, box, leader, point), ...])`` for one model's report
    figure.

    The boxes are the TEXT extents, measured off the drawn artists: an
    annotation's own window extent is the union of its text and its leader, and
    a leader is a hairline that is meant to run under things.

    ``leader`` is ``None`` where the label carries none, and otherwise the two
    ends of the line as DRAWN — taken off the arrow patch's own path, so what is
    measured is where the hairline stops rather than where it was aimed, with
    the vertex end second.
    """
    import contextlib as _c
    import io as _io
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MplFigure
    from matplotlib.text import Text
    from xslope.fileio import load_slope_data
    from xslope.plot import plot_inputs

    with _c.redirect_stdout(_io.StringIO()):
        slope_data = load_slope_data(xlsx)
    fig = MplFigure(figsize=figsize)
    FigureCanvasAgg(fig)
    with _c.redirect_stdout(_io.StringIO()):
        plot_inputs(slope_data, fig=fig, mode="shared", show_title=False,
                    frame="content", label_coordinates=True)
    fig.canvas.draw()
    ax = fig.axes[0]
    renderer = fig.canvas.get_renderer()

    def drawn(t, point):
        patch = getattr(t, "arrow_patch", None)
        if patch is None:
            return None
        ends = patch.get_path().transformed(patch.get_transform()).vertices
        a, b = tuple(ends[0]), tuple(ends[-1])
        if math.dist(a, point) < math.dist(b, point):
            a, b = b, a
        return a, b

    out = []
    for t in ax.texts:
        if t.get_gid() != "COORD_LABEL":
            continue
        point = tuple(ax.transData.transform(t.xy))
        out.append((t.get_text(), Text.get_window_extent(t, renderer),
                    drawn(t, point), point))
    return ax.get_window_extent(renderer), out


#: How far a coordinate label may sit from the point it names, in multiples of
#: its own text height, when nothing is near enough to push it out.
#:
#: The placement offers a label a character's width off its vertex and walks
#: outward only while that is objectionable, so an uncrowded vertex is labelled
#: at the tightest offset there is. Over the ninety labels the five sample
#: sections below carry at two figure sizes, that offset measures 0.56 to 0.61
#: text heights — a character is a little over half the height of the line it
#: sits on. The bound is set at 0.85: forty per cent clear of the widest hug
#: measured, so a change of face or of figure size does not trip it, and well
#: under the 1.08 heights an offset ring keyed to the line height rather than to
#: the character put every uncrowded label at, so the bound tells the two apart.
HUG = 0.85

#: How near its vertex a leader may stop, and how far short of it it must, in
#: multiples of the label's own text height.
#:
#: A callout approaches the point it names and stops just short: a hairline run
#: into a corner is ink on the section — at a dam crest the leader arrives along
#: the crest line and the two read as one bent line — while one that halts well
#: out has to be read across the gap to see what it points at. The trim is half
#: a character, which measures 0.31 text heights over the leaders the sample
#: sections below draw at two figure sizes. The bounds are set a third under
#: that and half again over it, so a change of face or of figure size does not
#: trip them while a pull-back of a single point — 0.14 heights, which the width
#: of the line drawn through the corner swallows whole — does not pass. The
#: upper bound is well inside the 0.56 heights the tightest label offset
#: measures, so a leader always stops nearer its corner than any label stands.
LEADER_GAP = (0.35, 0.55)


def test_coordinate_labels_are_placed_clear():
    """Every coordinate label is readable where it landed.

    The labels name the vertices of the section, and on a dam several of those
    vertices are a few feet apart — a crest corner, the core beneath it — so the
    placement is solved rather than offset by a constant. What the solution owes
    a reader is checked on the artists it produced, not on the pixels: no label
    lies over another, none runs off the axes it belongs to, any label the
    solver had to move away from its point carries a leader back to it, since a
    coordinate that could belong to either of two corners belongs to neither,
    and NO TWO LEADERS CROSS — a reader who follows the wrong line out of a
    crossing reads the wrong coordinate, and nothing on the figure says so. A
    leader that is drawn stops short of its vertex rather than running into it,
    by a gap measured off the type and bounded both ways: a line that meets the
    corner is ink on the section, and one that halts well out is a line the
    reader has to carry across the white.

    And a label HUGS the point it names unless something made it travel. A
    coordinate floating in the white is read against the drawing rather than
    against its corner, so the placement offers every label a character's width
    off its vertex and moves it only where that offer is objectionable — and a
    label it did move carries a leader. The two facts are checked from both ends:
    a label with no leader has not travelled, so it stands within :data:`HUG` of
    its own text height of its point; and on a section where every vertex has
    somewhere to put its coordinate, NO label carries a leader at all, so a
    placement that walked past a clear position for a softer preference is caught
    even where the label it moved is honestly labelled as having moved.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from xslope.plot import _segments_cross

    # A zoned dam (the crest and the core top are four points inside one label
    # width), a reinforced slope with two vertices two feet apart, a layered
    # section, and two bare slopes: five placements the same rule has to solve.
    # All but the dam have somewhere to put every coordinate — the reinforced
    # slope's two close pairs are resolved by putting the two labels on opposite
    # sides of their points, not by moving either — so nothing on them should
    # travel. Each is drawn at the report's figure size and again at three
    # quarters of it — the same section with the type forty per cent larger
    # against it, which is what a smaller window gives Studio.
    for (label, xlsx, roomy), size in itertools.product(
            (("the dam", SEEP_XLSX, False),
             ("the reinforced slope", REINF_XLSX, True),
             ("the layered section", NONCIRC_XLSX, True),
             ("the simple slope", SIMPLE1_XLSX, True),
             ("the zoned dam", DAM_XLSX, False)),
            ((11.0, 5.0), (8.0, 3.6))):
        label = f"{label} at {size[0]:g}x{size[1]:g}"
        frame, labels = _coord_labels(xlsx, figsize=size)
        if len(labels) < 4:
            fails.append(f"{label}: {len(labels)} coordinate labels — too few to "
                         f"exercise the placement")
            continue
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if labels[i][1].overlaps(labels[j][1]):
                    fails.append(f"{label}: {labels[i][0]} and {labels[j][0]} "
                                 f"are printed over each other")
        # The leader a reader follows: the vertex, to the middle of the box the
        # coordinate is printed in. Taken for every label, drawn or not — a
        # label sitting beside its point is read along that line just the same.
        runs = [(text, ((ax_px, ay_px),
                        ((box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2)))
                for text, box, _leader, (ax_px, ay_px) in labels]
        for i, (t1, r1) in enumerate(runs):
            for t2, r2 in runs[i + 1:]:
                if _segments_cross(r1, r2):
                    fails.append(f"{label}: the leaders of {t1} and {t2} cross")
        # And a label nearer some other vertex than the one it names carries a
        # leader: along a stack of vertices one label-height apart every label
        # can land beside its neighbour's point without any two leaders ever
        # meeting, and only the line back says which coordinate is whose.
        for text, box, leader, own in labels:
            def gap(p):
                return math.hypot(max(box.x0 - p[0], 0.0, p[0] - box.x1),
                                  max(box.y0 - p[1], 0.0, p[1] - box.y1))
            nearer = [t for t, _b, _l, p in labels if p != own and gap(p) < gap(own)]
            if nearer and not leader:
                fails.append(f"{label}: {text} sits nearer the vertex "
                             f"{nearer[0]} names than its own, with no leader "
                             f"tying it back")
        for text, box, leader, (ax_px, ay_px) in labels:
            if not frame.contains(box.x0, box.y0) or not frame.contains(box.x1, box.y1):
                fails.append(f"{label}: {text} is drawn off the axes "
                             f"({box.x0:.0f}..{box.x1:.0f}, "
                             f"{box.y0:.0f}..{box.y1:.0f} in "
                             f"{frame.x0:.0f}..{frame.x1:.0f}, "
                             f"{frame.y0:.0f}..{frame.y1:.0f})")
            # A label that had to travel needs a leader: the reach it may sit
            # inside without one is its own height, so the rule scales with the
            # type rather than with the model.
            reach = 1.5 * box.height
            gap = min(abs(ax_px - box.x0), abs(ax_px - box.x1)) \
                if not (box.x0 <= ax_px <= box.x1) else 0.0
            gap = max(gap, 0.0 if box.y0 <= ay_px <= box.y1
                      else min(abs(ay_px - box.y0), abs(ay_px - box.y1)))
            if gap > reach and not leader:
                fails.append(f"{label}: {text} sits {gap:.0f} px from the point "
                             f"it names with no leader tying it back")

        # A label with no leader HUGS the point it names. The two are the same
        # statement read from either end: the figure carries a leader for every
        # coordinate that had to move, so a coordinate standing on its own has
        # not moved, and it is beside its corner rather than adrift in the white.
        # Which of the two a given label gets is the placement's business; that
        # it is one or the other is the reader's.
        for text, box, leader, own in labels:
            if leader:
                continue
            gap = math.hypot(
                0.0 if box.x0 <= own[0] <= box.x1
                else min(abs(own[0] - box.x0), abs(own[0] - box.x1)),
                0.0 if box.y0 <= own[1] <= box.y1
                else min(abs(own[1] - box.y0), abs(own[1] - box.y1)))
            if gap > HUG * box.height:
                fails.append(f"{label}: {text} carries no leader and yet sits "
                             f"{gap / box.height:.2f} text heights off the point "
                             f"it names (the bound is {HUG:g})")

        # And a leader that IS drawn approaches its vertex without meeting it,
        # by :data:`LEADER_GAP` of its own text height — measured on the path
        # the hairline was drawn along, so what is checked is where the line
        # stops rather than where it was aimed. A leader shorter than the gap it
        # stops short by is a smudge and not a line: that label is against its
        # point already and carries no leader at all, which the hug bound above
        # then holds to.
        for text, box, leader, own in labels:
            if not leader:
                continue
            stop = math.dist(leader[1], own) / box.height
            run = math.dist(leader[0], leader[1]) / box.height
            if not LEADER_GAP[0] <= stop <= LEADER_GAP[1]:
                fails.append(f"{label}: the leader of {text} stops {stop:.2f} "
                             f"text heights from the point it names, outside "
                             f"{LEADER_GAP[0]:g} to {LEADER_GAP[1]:g}")
            elif run < stop:
                fails.append(f"{label}: the leader of {text} is {run:.2f} text "
                             f"heights long and stops {stop:.2f} short of its "
                             f"point — the trim has eaten the line")

        # And read the other way: on a section with room for every coordinate,
        # a leader is the mark of a label that travelled, and nothing had cause
        # to. One here means the placement passed over a clear position beside
        # the point for something it liked better further out.
        if roomy:
            travelled = [t for t, _b, leader, _p in labels if leader]
            if travelled:
                fails.append(f"{label}: {travelled} carry leaders on a section "
                             f"where every vertex has room for its own label")

    # The dam's crest cluster is the case the placement exists for: four
    # vertices inside one label width, so something has to move.
    _frame, labels = _coord_labels(SEEP_XLSX)
    from xslope.plot import _label_geometry
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(SEEP_XLSX)
    points, _segments = _label_geometry(slope_data)
    named = {f"({x:g}, {y:g})" for x, y in points}
    if {t for t, _b, _l, _p in labels} != named:
        fails.append(f"the labels {sorted(t for t, _b, _l, _p in labels)} are not "
                     f"the model's vertices {sorted(named)}")
    if not any(leader for _t, _b, leader, _p in labels):
        fails.append("no label on the dam carries a leader; the crest cluster is "
                     "the case the placement exists for and nothing moved")
    return fails


# --------------------------------------------------------------------------
# K. the seepage and finite element sections
#
# Both engines are documented from SHIPPED solutions: the sample models carry
# the companion files their solvers write ({base}_seep.csv, {base}_fem_*.csv),
# and reading one back is the same bundle a run emits. Nothing here re-solves —
# a check of what the report SAYS about a solution has no business spending a
# seepage iteration or an SSRM bisection to get one.
#
# The reading is xslope.report.solutions_from_sidecars, the same call a caller
# makes to report an already-solved model, so every check below runs on what a
# caller gets rather than on a recipe kept only here.
# --------------------------------------------------------------------------

#: A dam with a solved unconfined seepage analysis beside it. Its four materials
#: are all on the linear-front conductivity model.
SEEP_XLSX = os.path.join(_REPO, "docs", "lem", "files", "xslope_gsat_seep.xlsx")

#: The same dam with its two zones on van Genuchten instead — the corpus's model
#: whose relative conductivity genuinely falls through decades, and the one that
#: settles what the conductivity figures' ordinate scale is chosen by.
VG_SEEP_XLSX = os.path.join(_REPO, "docs", "seep", "files",
                            "xslope_earth_dam1_vg.xlsx")

#: A loaded slope with a solved SSRM run beside it.
FEM_XLSX = os.path.join(_REPO, "docs", "fem", "files",
                        "xslope_griffiths1_load.xlsx")

#: The two finite element models that carry one-dimensional members: six
#: reinforcement lines, and two piles. Neither ships a solution beside it, so
#: these two are the exception to the rule above — one gravity trial each, a
#: second or two, which is the only way to put a real bar force or a real pile
#: moment in front of the report.
FEM_REINF_XLSX = os.path.join(_REPO, "docs", "fem", "files",
                              "xslope_reinforce_fem.xlsx")
FEM_PILES_XLSX = os.path.join(_REPO, "docs", "fem", "files",
                              "xslope_piles_fem.xlsx")
FEM_1D_MODELS = (FEM_REINF_XLSX, FEM_PILES_XLSX)

_ENGINE = {}


def _restored(xlsx, notes=None):
    """``(slope_data, solutions)`` for a model, with every engine solution that
    ships beside it read back through the public helper.

    One recipe, used here and by anyone assembling a report of an already-solved
    model: whatever the checks below assert about a restored bundle, they assert
    about the bundle a caller gets.
    """
    from xslope.fileio import load_slope_data
    from xslope.report import solutions_from_sidecars

    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(xlsx)
        solutions = solutions_from_sidecars(xlsx, slope_data, notes)
    return slope_data, solutions


def _seep_bundle(xlsx=SEEP_XLSX):
    """``(slope_data, bundle)`` for a seepage analysis, read back from the
    solution shipped beside the model. The first boundary condition set: a model
    that carries two is read by the checks that are about two."""
    key = ("seep", xlsx)
    if key not in _ENGINE:
        slope_data, solutions = _restored(xlsx)
        _ENGINE[key] = (slope_data, solutions["seep"][0])
    return _ENGINE[key]


def _fem_bundle(xlsx=FEM_XLSX):
    """``(slope_data, bundle)`` for a strength reduction run, read back from the
    solution shipped beside the model — the shape Studio's FEM runner emits."""
    key = ("fem", xlsx)
    if key not in _ENGINE:
        slope_data, solutions = _restored(xlsx)
        _ENGINE[key] = (slope_data, solutions["fem"])
    return _ENGINE[key]


def _fem_1d_bundle(xlsx):
    """``(slope_data, bundle)`` for a model carrying reinforcement or piles.

    Neither model ships a solution beside it, so one gravity trial is solved
    here and cached for every check that reads it. A trial that left every bar
    and every pile at zero force would exercise the report's member sections on
    a mechanism that never engaged, so the checks that read them assert the
    forces are real.
    """
    key = ("fem1d", xlsx)
    if key in _ENGINE:
        return _ENGINE[key]
    from xslope.fem import build_fem_data, solve_fem
    from xslope.fileio import load_slope_data

    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(xlsx)
        fem_data = build_fem_data(slope_data, slope_data["mesh"])
        solution = solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=60)
    bundle = {"fem_data": fem_data, "solution": solution, "FS": None,
              "analysis": "single", "failure_solution": None}
    _ENGINE[key] = (slope_data, bundle)
    return _ENGINE[key]


def _fem_any_bundle(xlsx):
    """The finite element bundle for a model: read back where a solution ships
    beside it, solved where none does."""
    return _fem_1d_bundle(xlsx) if xlsx in FEM_1D_MODELS else _fem_bundle(xlsx)


def _built_report(slope_data, solutions, options):
    """A report of a model built here rather than read from a file — how a check
    asks what a model with one property changed would report."""
    from xslope.report import build_report

    opts = dict(FAST_FIGURES)
    opts.update(options or {})
    tmp = tempfile.mkdtemp(prefix="xslope_built_")
    with contextlib.redirect_stdout(io.StringIO()):
        return build_report(slope_data, solutions, opts, tmp)


def _engine_report(engine, options=None, bundle=None, xlsx=None):
    """A report of one engine's model with that engine's section built.

    The limit equilibrium section is switched off: these models are loaded, not
    solved for stability, and a report of a seepage run is a real thing to ask
    for.
    """
    from xslope.report import build_report

    xlsx = xlsx or (SEEP_XLSX if engine == "seep" else FEM_XLSX)
    slope_data, default = (_seep_bundle(xlsx) if engine == "seep"
                           else _fem_any_bundle(xlsx))
    opts = {"input_path": xlsx, "lem": False, "pd_figure": False}
    opts.update(FAST_FIGURES)
    opts.update(options or {})
    tmp = tempfile.mkdtemp(prefix=f"xslope_{engine}_")
    with contextlib.redirect_stdout(io.StringIO()):
        return build_report(slope_data, {engine: bundle or default}, opts, tmp)


def _titles(report):
    return [t for _lvl, t in report.section_titles()]


def _prose(report):
    return [b.text for b in report.blocks("prose")]


def _planned_matches(report, engine, options=None, bundle=None, xlsx=None,
                     slope_data=None):
    """``planned_figures`` against what the build produced, for one engine."""
    from xslope.report import planned_figures, resolve_options

    xlsx = xlsx or (SEEP_XLSX if engine == "seep" else FEM_XLSX)
    read, default = (_seep_bundle(xlsx) if engine == "seep"
                     else _fem_any_bundle(xlsx))
    slope_data = read if slope_data is None else slope_data
    opts = {"input_path": xlsx, "lem": False, "pd_figure": False}
    opts.update(FAST_FIGURES)
    opts.update(options or {})
    planned = planned_figures(slope_data, {engine: bundle or default},
                              resolve_options(opts))
    return planned, len(report.figures())


def _kr_axes(materials, abscissa="suction"):
    """The Axes a set of materials' conductivity curves are drawn on, off the
    real figure rather than off the code that draws it."""
    import matplotlib.figure as mplfig
    from xslope.plot import plot_material_kr_set
    fig = mplfig.Figure()
    plot_material_kr_set(materials, fig=fig, abscissa=abscissa,
                         show_legend=False)
    return fig.axes[0]


def _kr_ordinate_scale_checks():
    """The conductivity figures are read on the scale their models are shaped by.

    The linear front is a straight ramp: kr falls linearly from 1 at the ground
    water surface to kr0 one head interval below it, and the position of that
    front is the whole of what the model says. Drawn against a logarithmic
    ordinate the ramp becomes a decaying curve that falls off a cliff, and the
    front cannot be read off it at all. van Genuchten and Gardner do fall through
    decades, and a linear ordinate lays all of that flat against zero.

    So the scale follows the models assigned — linear where every curve on the
    axes is a linear front, logarithmic as soon as one of them spans decades —
    and the same rule governs the report's paired figures and the material
    editor's single one, so an engineer reads one curve on one axis wherever it
    is met.
    """
    fails = []
    from xslope.plot import (plot_material_kr, _kr_yscale, _kr_model)
    from xslope.report import _kr_materials
    import xslope.plot as plot_mod
    import matplotlib.figure as mplfig
    import numpy as np

    linear_mats = _kr_materials(load_slope_data_cached(SEEP_XLSX))
    vg_mats = _kr_materials(load_slope_data_cached(VG_SEEP_XLSX))
    # A set of both is not in the corpus and does not need to be: the rule is
    # about composition, so the composition is built here.
    mixed_mats = (linear_mats[:1] or []) + (vg_mats[:1] or [])

    if not linear_mats or not all(_kr_model(m) == "lf" for m in linear_mats):
        fails.append(f"the seepage sample is not all linear-front — "
                     f"{[_kr_model(m) for m in linear_mats]} — so the scale "
                     f"the defect was found on is untested")
    if not vg_mats or not any(_kr_model(m) == "vg" for m in vg_mats):
        fails.append(f"{os.path.basename(VG_SEEP_XLSX)} carries no van "
                     f"Genuchten material, so the log ordinate is untested")

    # Each composition, on both abscissae: the pair is one statement of the
    # model, and a pair read on two ordinates does not mirror.
    for what, mats, want in (("all linear front", linear_mats, "linear"),
                             ("all van Genuchten", vg_mats, "log"),
                             ("mixed", mixed_mats, "log")):
        if not mats:
            continue
        for abscissa in ("suction", "head"):
            ax = _kr_axes(mats, abscissa)
            got = ax.get_yscale()
            if got != want:
                fails.append(f"the {what} materials are drawn against "
                             f"{abscissa} on a {got} kr ordinate, not {want}")
            lo, hi = ax.get_ylim()
            if want == "linear" and not (abs(lo) < 1e-12 and 1.0 < hi <= 1.2):
                fails.append(f"the {what} kr ordinate spans {lo:g} to {hi:g}; "
                             f"kr runs 0 to 1 and the axis should show that")
            if want == "log" and not (0 < lo < 1 and hi > 1.0):
                fails.append(f"the {what} kr ordinate spans {lo:g} to {hi:g} on "
                             f"a log axis, which does not frame the curves")

        # The material editor draws the same curve on the same axis. Studio and
        # the report disagreeing about a scale is one model with two shapes.
        one = mats[0]
        fig = mplfig.Figure()
        ax = fig.add_subplot(111)
        plot_material_kr(ax, one)
        alone = _kr_yscale([one])
        if ax.get_yscale() != alone:
            fails.append(f"the material editor draws {_kr_model(one)!r} on a "
                         f"{ax.get_yscale()} ordinate; the report's rule gives "
                         f"{alone}")

    # A linear front on a mixed axes is still its own numbers — a log ordinate
    # is the right axis for that set, and the ramp is correctly plotted on it.
    if mixed_mats and linear_mats:
        alone = {ln.get_label(): ln.get_data()[1]
                 for ln in _kr_axes(linear_mats[:1]).get_lines()}
        together = {ln.get_label(): ln.get_data()[1]
                    for ln in _kr_axes(mixed_mats).get_lines()}
        for label, kr in alone.items():
            other = together.get(label)
            if other is None or not np.allclose(kr, other):
                fails.append(f"the linear-front curve for {label!r} changes when "
                             f"it is drawn beside a model that spans decades")

    # Mutation: the check has to be able to fail. Forcing the log ordinate back
    # onto the all-linear-front set is the figure the defect was reported on,
    # and it has to go red here.
    real = plot_mod._kr_yscale
    plot_mod._kr_yscale = lambda materials: "log"
    try:
        forced = [_kr_axes(linear_mats, a).get_yscale()
                  for a in ("suction", "head")] if linear_mats else []
    finally:
        plot_mod._kr_yscale = real
    if linear_mats and set(forced) != {"log"}:
        fails.append(f"the ordinate scale is not read from the models: forcing "
                     f"log still drew {forced}")
    if linear_mats and _kr_yscale(linear_mats) != "linear":
        fails.append("the rule itself does not call the all-linear-front set "
                     "linear, so the mutation above proves nothing")
    return fails


def test_seep_section():
    """A report of a seepage run says what was solved, on what, and what flowed."""
    fails = []
    _slope_data, bundle = _seep_bundle()
    report = _engine_report("seep")

    expected = [(1, "Traceability"), (1, "Project Definition"),
                (1, "Seepage Analysis"), (2, "Analysis Inputs"), (2, "Results")]
    got = report.section_titles()
    if got != expected:
        fails.append(f"the seepage report's sections are {got}, expected {expected}")

    # What was solved: the sample is unconfined, and both boundary counts are
    # read off the same array the solver used rather than stated as a constant
    # here.
    bc_type = list(bundle["seep_data"]["bc_type"])
    n_head = sum(1 for t in bc_type if int(t) == 1)
    n_exit = sum(1 for t in bc_type if int(t) == 2)
    if not n_exit:
        fails.append("the sample is not an unconfined problem, so the report's "
                     "unconfined wording is never exercised")
    text = " ".join(_prose(report))
    if "unconfined" not in text:
        fails.append("the seepage section does not say the problem was unconfined")
    for count in (f"{n_head:,}", f"{n_exit:,}"):
        if count not in text:
            fails.append(f"the boundary node count {count} is not stated: {text!r}")

    # The flow: the number the solution carries, in the prose, in bold.
    q = bundle["solution"]["flowrate"]
    stated = [b for b in report.blocks("prose") if f"{q:.4g}" in b.text]
    if not stated:
        fails.append(f"the computed flow of {q:.4g} is stated nowhere in the "
                     f"seepage section: {text!r}")
    elif not any(f"{q:.4g}" in " ".join(b.bold) for b in stated):
        fails.append("the flow is stated but not set in bold")

    # The inputs: the mesh, and the conductivities the flow was solved with.
    inputs = next((s for s in report.sections if s.title == "Seepage Analysis"),
                  None)
    inputs = next((c for c in (inputs.children if inputs else [])
                   if c.title == "Analysis Inputs"), None)
    if inputs is None:
        fails.append("the seepage section has no Analysis Inputs")
    else:
        kv = [b for b in inputs.blocks if b.kind == "keyvalues"]
        labels = [l for b in kv for l, _v in b.items]
        if "Mesh" not in labels:
            fails.append(f"the seepage inputs name no mesh: {labels}")
        table = next((b for b in inputs.blocks if b.kind == "table"), None)
        if table is None:
            fails.append("the seepage inputs carry no material properties table")
        elif not any("k" in h for h in table.headers):
            fails.append(f"the seepage material table has no conductivity "
                         f"column: {table.headers}")

    # The figures the build produced are the figures it planned, and each is a
    # different reading of the run: the model the flow was solved on, the
    # conductivity curves it reduces each material by, the mesh with the boundary
    # conditions on it, and the field that came out.
    planned, drawn = _planned_matches(report, "seep")
    if planned != drawn:
        fails.append(f"the seepage report planned {planned} figures and built {drawn}")
    sources = [f.source for f in report.figures()]
    for wanted in ("seep model", "seep kr", "seep kr_head", "seepage bc1 mesh",
                   "seepage bc1"):
        if wanted not in sources:
            fails.append(f"the seepage report has no {wanted!r} figure: {sources}")
    if drawn != 5:
        fails.append(f"the seepage report drew {drawn} figures, expected the "
                     f"model, the unsaturated conductivity curves against "
                     f"suction and against pressure head, the mesh with its "
                     f"boundary conditions, and the flow net")

    # The unsaturated curves stand in the inputs, after the properties table
    # whose parameters they draw and before the mesh they are solved on, and the
    # two conventions stand together: they are one statement of the conductivity
    # model, and a page between them makes them two.
    if inputs is not None:
        order = [b.source if b.kind == "figure" else b.kind
                 for b in inputs.blocks if b.kind in ("figure", "table")]
        want = ["seep model", "table", "seep kr", "seep kr_head",
                "seepage bc1 mesh"]
        if order != want:
            fails.append(f"the seepage inputs are ordered {order}, not {want}")

    # Both conductivity figures are cited, from one sentence that presents them
    # as the same models in the two conventions. A figure drawn and not named is
    # a page the reader is never sent to.
    numbers = {f.source: f.number for f in report.figures()}
    citing = [b for b in report.blocks("prose")
              if all(f"Figure {numbers[s]}" in b.text
                     for s in ("seep kr", "seep kr_head") if s in numbers)
              and "conductivity" in b.text]
    if not citing:
        fails.append("no sentence cites both conductivity figures: "
                     + repr([b.text for b in report.blocks("prose")
                             if "conductivity" in b.text]))
    else:
        said = citing[0].text
        for named in ("matric suction", "pressure head"):
            if named not in said:
                fails.append(f"the sentence citing the conductivity figures "
                             f"does not name {named!r}: {said!r}")
        from xslope.report import cite_anchor
        targets = [t for _text, t in citing[0].links]
        for source in ("seep kr", "seep kr_head"):
            if f"#{cite_anchor('Figure', numbers[source])}" not in targets:
                fails.append(f"the sentence naming the {source!r} figure does "
                             f"not link to it: {targets}")

    # The two are the same curves under two sign conventions, so a material keeps
    # its color and its dash across the pair. Read off the artists, not asserted
    # of the code that draws them: the report ships PNGs, and a color that only
    # matches in the source is a pair the reader cannot follow.
    from xslope.plot import plot_material_kr_set
    from xslope.report import _kr_materials
    import matplotlib.figure as mplfig
    import numpy as np
    mats = _kr_materials(_slope_data)
    if not mats:
        fails.append("the seepage sample carries no unsaturated material, so "
                     "the conductivity figures prove nothing")
    else:
        keyed = []
        for abscissa in ("suction", "head"):
            f = mplfig.Figure()
            plot_material_kr_set(mats, fig=f, abscissa=abscissa,
                                 show_legend=False)
            ax = f.axes[0]
            keyed.append({ln.get_label(): (ln.get_color(), ln.get_linestyle())
                          for ln in ax.get_lines()})
        if keyed[0] != keyed[1]:
            fails.append(f"the conductivity pair keys its materials "
                         f"differently: {keyed[0]} against {keyed[1]}")
        if len(keyed[0]) < 1:
            fails.append("the conductivity figures drew no labelled curve")

    # And they are the same numbers: the head figure is the suction figure
    # mirrored, drawn over the negative abscissa and reaching zero, where the
    # material is saturated.
    if mats:
        f = mplfig.Figure()
        plot_material_kr_set(mats, fig=f, abscissa="head", show_legend=False,
                             unit_labels={"length": "ft"})
        ax = f.axes[0]
        if "pressure head" not in ax.get_xlabel():
            fails.append(f"the head figure's abscissa is labelled "
                         f"{ax.get_xlabel()!r}")
        if "(ft)" not in ax.get_xlabel():
            fails.append(f"the head abscissa carries no length unit: "
                         f"{ax.get_xlabel()!r}")
        lo, hi = ax.get_xlim()
        if not (lo < 0 and abs(hi) < 1e-9):
            fails.append(f"the head figure spans {lo:g} to {hi:g}, not the "
                         f"negative range up to saturation at zero")
        g = mplfig.Figure()
        plot_material_kr_set(mats, fig=g, abscissa="suction", show_legend=False)
        by_label = {ln.get_label(): ln for ln in g.axes[0].get_lines()}
        for line in ax.get_lines():
            twin = by_label.get(line.get_label())
            if twin is None:
                continue
            x, y = line.get_data()
            u, v = twin.get_data()
            if not (np.allclose(x, -np.asarray(u))
                    and np.allclose(y, np.asarray(v))):
                fails.append(f"the head curve for {line.get_label()!r} is not "
                             f"the suction curve mirrored")
            if abs(float(y[0]) - 1.0) > 1e-6:
                fails.append(f"the head curve for {line.get_label()!r} does not "
                             f"reach kr = 1 at saturation: kr = {float(y[0]):g}")

    fails += _kr_ordinate_scale_checks()

    # The flow net is a flow net: no element edges over the field, and the base
    # material the flow lines are scaled to is chosen, not left at one.
    from xslope.plot_seep import flownet_base_material
    passed = {}
    import xslope.plot_seep as ps
    real = ps.plot_seep_solution

    def spy(sd, sol, **kw):
        passed.update(kw)
        return real(sd, sol, **kw)

    ps.plot_seep_solution = spy
    try:
        _engine_report("seep", options={"seep_inputs_figure": False,
                                        "seep_mesh_figure": False})
    finally:
        ps.plot_seep_solution = real
    if passed.get("mesh") is not False:
        fails.append(f"the flow net is drawn with mesh={passed.get('mesh')!r}; "
                     f"element edges chop the contours and hide the field")
    chosen = flownet_base_material(bundle["seep_data"], bundle["solution"])
    if passed.get("base_mat") != chosen:
        fails.append(f"the flow net is scaled to material "
                     f"{passed.get('base_mat')!r}, not the {chosen} its "
                     f"conductivities call for")

    # Each option carries its own figures, and switching one off takes only
    # those. The conductivity option governs the pair — the same models in the
    # two conventions are one statement, and half of it is not a report of the
    # conductivity model — so it takes two and the rest take one.
    for option, gone in (("seep_inputs_figure", ("seep model",)),
                         ("seep_kr_figure", ("seep kr", "seep kr_head")),
                         ("seep_mesh_figure", ("seepage bc1 mesh",)),
                         ("seep_flownet", ("seepage bc1",))):
        off = _engine_report("seep", options={option: False})
        got = [f.source for f in off.figures()]
        for source in gone:
            if source in got:
                fails.append(f"{option}=False still drew the {source!r} figure")
        if len(got) != 5 - len(gone):
            fails.append(f"{option}=False left {len(got)} figures, not the "
                         f"other {5 - len(gone)}: {got}")
        planned, drawn = _planned_matches(off, "seep", options={option: False})
        if planned != drawn:
            fails.append(f"{option}=False planned {planned} figures and built "
                         f"{drawn}")

    # A model that carries no unsaturated parameters — a confined problem, where
    # the flow never leaves the saturated zone — has no curve to draw, and
    # NEITHER figure is drawn: it is absent rather than blank, with the sentence
    # that would cite it.
    saturated = copy.deepcopy(_slope_data)
    for m in saturated.get("materials") or []:
        m.update(kr0=0.0, h0=0.0, vg_a=0.0, vg_n=0.0)
    dry = _built_report(saturated, {"seep": bundle},
                        {"input_path": SEEP_XLSX, "lem": False,
                         "pd_figure": False})
    for source in ("seep kr", "seep kr_head"):
        if source in [f.source for f in dry.figures()]:
            fails.append(f"a model with no unsaturated material still drew the "
                         f"{source!r} conductivity figure")
    dry_text = " ".join(_prose(dry))
    for named in ("matric suction", "pressure head"):
        if named in dry_text:
            fails.append(f"a model with no unsaturated material still describes "
                         f"the {named} conductivity figure it did not draw")
    planned, drawn = _planned_matches(dry, "seep", slope_data=saturated)
    if planned != drawn:
        fails.append(f"a saturated model planned {planned} figures and built "
                     f"{drawn}")
    return fails


def test_fem_section():
    """A report of a strength reduction run states its factor of safety, in bold,
    and a report of a single trial states no factor of safety at all."""
    fails = []
    from xslope.report import FEM_PANELS
    _slope_data, bundle = _fem_bundle()
    report = _engine_report("fem")

    expected = [(1, "Traceability"), (1, "Project Definition"),
                (1, "Deformation and Strength Reduction"),
                (2, "Analysis Inputs"), (2, "Loads"), (2, "Results")]
    got = report.section_titles()
    if got != expected:
        fails.append(f"the SSRM report's sections are {got}, expected {expected}")

    fs = bundle["FS"]
    stated = [b for b in report.blocks("prose") if f"{fs:.3f}" in b.text]
    if not stated:
        fails.append(f"the strength reduction factor of safety {fs:.3f} is stated "
                     f"nowhere: {_prose(report)!r}")
    elif not any(f"{fs:.3f}" in " ".join(b.bold) for b in stated):
        fails.append(f"the factor of safety {fs:.3f} is stated but not set in "
                     f"bold, which is how the report states one")
    if "factor of safety" not in " ".join(_prose(report)):
        fails.append("the SSRM section never uses the words factor of safety")

    # The inputs: the stiffness columns are what make this table the finite
    # element one rather than a second copy of the strength table.
    sec = next((s for s in report.sections
                if s.title == "Deformation and Strength Reduction"), None)
    inputs = next((c for c in (sec.children if sec else [])
                   if c.title == "Analysis Inputs"), None)
    if inputs is None:
        fails.append("the finite element section has no Analysis Inputs")
    else:
        table = next((b for b in inputs.blocks if b.kind == "table"), None)
        if table is None:
            fails.append("the finite element inputs carry no properties table")
        else:
            for header in ("E", "ν"):
                if not any(h == header or h.startswith(header + " ")
                           for h in table.headers):
                    fails.append(f"the finite element material table has no "
                                 f"{header} column: {table.headers}")

    planned, drawn = _planned_matches(report, "fem")
    if planned != drawn:
        fails.append(f"the SSRM report planned {planned} figures and built {drawn}")
    sources = [f.source for f in report.figures()]
    for wanted in ("fem model", "fem mesh", "fem run1 deformation",
                   "fem run1 shear_strain", "fem run1 displace_vector"):
        if wanted not in sources:
            fails.append(f"the SSRM report has no {wanted!r} figure: {sources}")
    if drawn != 2 + len(FEM_PANELS):
        fails.append(f"the SSRM report drew {drawn} figures, expected the model, "
                     f"the mesh and the {len(FEM_PANELS)} result panels")

    # The mesh figure is Studio's FEM data view, from the same function, with the
    # boundary conditions on it: a mesh drawn without them does not say what the
    # section was held by, and a report that draws its own version of a plot
    # Studio already draws is a second answer to one question.
    import xslope.plot_fem as pf
    real_fem_data = pf.plot_fem_data
    seen = {}

    def spy_fem_data(fem_data, **kw):
        seen.update(kw)
        return real_fem_data(fem_data, **kw)

    pf.plot_fem_data = spy_fem_data
    try:
        _engine_report("fem", options={"fem_figure": False,
                                       "fem_inputs_figure": False})
    finally:
        pf.plot_fem_data = real_fem_data
    if not seen:
        fails.append("the finite element mesh figure is not drawn by "
                     "plot_fem_data, the function Studio's FEM data view draws")
    elif seen.get("show_bc") is not True:
        fails.append(f"the mesh figure is drawn with show_bc="
                     f"{seen.get('show_bc')!r}; the boundary conditions are "
                     f"half of what the figure says")
    mesh_caption = next((f.caption for f in report.figures()
                         if f.source == "fem mesh"), "")
    if "boundary condition" not in mesh_caption.lower():
        fails.append(f"the mesh figure is captioned {mesh_caption!r}, which "
                     f"does not name the boundary conditions it carries")

    # Each figure carries its own option, and switching one off takes only it.
    for option, gone in (("fem_inputs_figure", 1), ("fem_mesh_figure", 1),
                         ("fem_figure", len(FEM_PANELS))):
        off = _engine_report("fem", options={option: False})
        planned, off_drawn = _planned_matches(off, "fem", options={option: False})
        if planned != off_drawn:
            fails.append(f"{option}=False planned {planned} figures and built "
                         f"{off_drawn}")
        if off_drawn != drawn - gone:
            fails.append(f"{option}=False left {off_drawn} figures, not the "
                         f"{drawn - gone} the other options draw")

    # A single trial is not a strength reduction run: it reports displacements
    # and claims no factor of safety, and the heading says so.
    single = dict(bundle, analysis="single", FS=None)
    plain = _engine_report("fem", bundle=single)
    titles = _titles(plain)
    if "Deformation and Strength Reduction" in titles:
        fails.append("a single-trial run is headed as a strength reduction run")
    if "Deformation Analysis" not in titles:
        fails.append(f"a single-trial run has no deformation section: {titles}")
    text = " ".join(_prose(plain))
    if f"{fs:.3f}" in text:
        fails.append("a single-trial run reports the strength reduction factor "
                     "of safety anyway")
    if "no factor of safety" not in text:
        fails.append(f"a single-trial run does not say it reports no factor of "
                     f"safety: {text!r}")
    return fails


def test_member_detail_figures_are_readable():
    """Every member's detail figure states its peak where a reader can read it.

    The peak is a point ON the profile, and on a member well inside its capacity
    that profile runs along the envelope's own descent — so a label offset by
    quadrant lands on the dashes it is meant to be read against. Checked on the
    artists: the label's box stays inside its panel, covers no other label and
    not the legend, and crosses none of the curves the panel draws. Six
    reinforcement lines and two piles, because where the peak falls is what
    decides where its label can go.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MplFigure
    from matplotlib.text import Text
    from xslope.plot import _box_crosses_segment
    from xslope.plot_fem_details import plot_detail
    from xslope.report import FIGURE_SIZE

    def utilized(profile, share=0.89):
        """The same line worked up to ``share`` of its capacity.

        Under gravity alone the sample's bars sit at a few percent, and a peak
        near the axis is the easy case: the label has the whole panel over it. A
        bar near capacity puts its peak ON the envelope, which is where a rule
        that picks a quadrant lands the label on the dashes it is meant to be
        read against — the case this exists for, made from the real profile by
        scaling the one series it is a peak of.
        """
        import numpy as np
        peak = profile.get("peak_utilization")
        if not peak or not np.isfinite(peak) or peak <= 0:
            return None
        k = share / peak
        out = dict(profile)
        out["T"] = np.asarray(profile["T"], dtype=float) * k
        out["peak_T"] = float(profile["peak_T"]) * k
        out["peak_utilization"] = share
        return out

    for xlsx, kind in ((FEM_REINF_XLSX, "reinforcement"),
                       (FEM_PILES_XLSX, "pile")):
        slope_data, bundle = _fem_1d_bundle(xlsx)
        profiles = _profiles(slope_data, bundle, kind)
        if not profiles:
            fails.append(f"{os.path.basename(xlsx)} carries no {kind} profile")
            continue
        if kind == "reinforcement":
            worked = [p for p in (utilized(q) for q in profiles) if p]
            if not worked:
                fails.append("no reinforcement line could be worked up to near "
                             "its capacity; the crowded case goes untested")
            profiles = profiles + worked
        for profile in profiles:
            fig = MplFigure(figsize=FIGURE_SIZE)
            FigureCanvasAgg(fig)
            with contextlib.redirect_stdout(io.StringIO()):
                plot_detail(profile, fig=fig)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            where = f"{os.path.basename(xlsx)} {profile['label']}"

            # The peak label is the bold one; it is the only annotation these
            # panels place against the data rather than against the axes.
            peaks = [(ax, t) for ax in fig.axes for t in ax.texts
                     if t.get_fontweight() == "bold"]
            if len(peaks) != 1:
                fails.append(f"{where}: {len(peaks)} peak labels on the figure")
                continue
            ax, label = peaks[0]
            box = Text.get_window_extent(label, renderer)
            frame = ax.get_window_extent(renderer)
            if not (frame.x0 <= box.x0 and box.x1 <= frame.x1
                    and frame.y0 <= box.y0 and box.y1 <= frame.y1):
                fails.append(f"{where}: the peak label {label.get_text()!r} runs "
                             f"outside the panel it belongs to")
            for other in ax.texts:
                if other is label:
                    continue
                if box.overlaps(Text.get_window_extent(other, renderer)):
                    fails.append(f"{where}: the peak label is printed over "
                                 f"{other.get_text()!r}")
            legend = ax.get_legend()
            if legend is not None and box.overlaps(legend.get_window_extent(renderer)):
                fails.append(f"{where}: the peak label is printed over the legend")
            # The two series the peak is read AGAINST: the profile it is a point
            # of, and the capacity it is a fraction of. A label over either
            # hides the comparison it was printed to make. (A hairline the label
            # cannot avoid — a step's vertical riser in a crowded panel — sits
            # under an opaque backing and hides nothing, so it is not asked
            # about here.)
            crossed = []
            for line in ax.lines:
                if line.get_gid() not in ("DETAIL_PROFILE", "DETAIL_CAPACITY"):
                    continue
                name = str(line.get_label())
                pts = line.get_xydata()
                if pts is None or len(pts) < 2:
                    continue
                px = ax.transData.transform(pts)
                if any(_box_crosses_segment(box, a, b)
                       for a, b in zip(px[:-1], px[1:])):
                    crossed.append(name)
            if crossed:
                fails.append(f"{where}: the peak label lies across {crossed}")
            if label.get_bbox_patch() is None:
                fails.append(f"{where}: the peak label carries no backing; where "
                             f"a panel leaves it nowhere clear it dissolves into "
                             f"what is behind it")
    return fails


def _profiles(slope_data, bundle, kind):
    """The member profiles the report builds its member table from — read
    through the builder's own reader, so a check that compares a printed number
    to a profile is comparing it to the series it was printed from."""
    from xslope.report import _detail_profiles
    return _detail_profiles(slope_data, bundle, kind)


def _member_section(report, title):
    """The finite element run's Reinforcement or Piles subsection, or None."""
    for section in report.sections:
        for _lvl, sec in section.walk():
            if sec.title == title:
                return sec
    return None


def _column(table, prefix):
    """The index of the column whose header starts with ``prefix``, or None."""
    for i, head in enumerate(table.headers):
        if str(head) == prefix or str(head).startswith(prefix + " "):
            return i
    return None


def _member_faults(table, profiles, force_prefix, force_key):
    """What a member table has to say about the profiles it was built from.

    One row per member, in the order the solver assigns them, carrying that
    member's own peak force and its own utilization. The last rule is the
    honest one: a row cannot claim a share of capacity that the force printed
    beside it could not have produced.
    """
    faults = []
    if len(table.rows) != len(profiles):
        return [f"{table.caption!r} has {len(table.rows)} rows for "
                f"{len(profiles)} members"]
    force_col = _column(table, force_prefix)
    util_col = _column(table, "Utilization")
    if force_col is None or util_col is None:
        return [f"{table.caption!r} has no {force_prefix} or Utilization "
                f"column: {table.headers}"]
    for row, profile in zip(table.rows, profiles):
        where = f"{table.caption!r}, {profile['label']}"
        if row[0] != profile["label"]:
            faults.append(f"{where}: the row is headed {row[0]!r}")
        util = profile.get("peak_utilization")
        stated = row[util_col]
        want = "" if util is None else f"{util:.0%}"
        if stated != want:
            faults.append(f"{where}: the utilization is printed {stated!r}, "
                          f"and the profile peaks at {want!r}")
        force = profile.get(force_key)
        try:
            printed = float(row[force_col].replace(",", ""))
        except ValueError:
            faults.append(f"{where}: the force column reads {row[force_col]!r}")
            continue
        if force is not None and abs(printed - float(force)) > (
                0.05 + 0.005 * abs(float(force))):
            faults.append(f"{where}: the force is printed {printed:,.1f} and "
                          f"the profile peaks at {float(force):,.1f}")
        if stated not in ("", "0%") and printed == 0.0:
            faults.append(f"{where}: the row claims {stated} of capacity with "
                          f"no force in the member")
    return faults


def test_fem_members_are_reported():
    """A finite element run that carries reinforcement lines or piles reports
    what the solution put in them, and one that carries neither reports nothing.

    The forces have to be REAL: a bar the mesh never loaded would let every
    number in the table be zero and every rule about it hold, so the mechanism
    is asserted to have engaged before anything is asked about how it is
    printed.
    """
    fails = []
    from xslope.fem_details import list_lines

    # --- reinforcement: six lines, each drawn --------------------------------
    slope_data, bundle = _fem_1d_bundle(FEM_REINF_XLSX)
    lines = list_lines(bundle["fem_data"], bundle["solution"], slope_data,
                       field_state="failure")
    n_lines = len(lines)
    if n_lines != len(slope_data.get("reinforcement_lines") or []):
        fails.append(f"the model owns elements for {n_lines} of its "
                     f"{len(slope_data.get('reinforcement_lines') or [])} "
                     f"reinforcement lines")
    if n_lines < 2:
        fails.append(f"the reinforcement model carries {n_lines} line; a report "
                     f"of one member proves nothing about drawing them all")

    report = _engine_report("fem", xlsx=FEM_REINF_XLSX)
    sec = _member_section(report, "Reinforcement Forces")
    if sec is None:
        fails.append(f"a run carrying {n_lines} reinforcement lines has no "
                     f"Reinforcement subsection: {_titles(report)}")
        return fails

    profiles = _profiles(slope_data, bundle, "reinforcement")
    peak = max((abs(p["peak_T"]) for p in profiles if p.get("peak_T")),
               default=0.0)
    if peak <= 0.0:
        fails.append("every reinforcement line came out of the solve at zero "
                     "force; the fixture never engages the bars and proves "
                     "nothing about reporting them")
    table = next((b for b in sec.blocks if b.kind == "table"), None)
    if table is None:
        fails.append("the Reinforcement subsection carries no table")
    else:
        fails += _member_faults(table, profiles, "Force", "peak_T")
        for header in ("T_max",):
            if _column(table, header) is None:
                fails.append(f"the reinforcement table declares no {header} "
                             f"capacity: {table.headers}")

    figures = [b for b in sec.blocks if b.kind == "figure"]
    if len(figures) != len(profiles):
        fails.append(f"{n_lines} lines drew {len(figures)} detail figures; every "
                     f"line the analysis solved is drawn")
    else:
        captioned = " | ".join(f.caption for f in figures)
        for profile in profiles:
            if profile["label"] not in captioned:
                fails.append(f"{profile['label']} has no detail figure: "
                             f"{captioned!r}")
    planned, drawn = _planned_matches(report, "fem", xlsx=FEM_REINF_XLSX)
    if planned != drawn:
        fails.append(f"the reinforcement report planned {planned} figures and "
                     f"built {drawn}")

    # --- piles: two, both drawn ---------------------------------------------
    pile_data, pile_bundle = _fem_1d_bundle(FEM_PILES_XLSX)
    pile_profiles = _profiles(pile_data, pile_bundle, "pile")
    if not pile_profiles:
        fails.append("the pile model carries no pile the solver owns elements "
                     "for")
        return fails
    if max((abs(p.get("max_moment") or 0.0) for p in pile_profiles),
           default=0.0) <= 0.0:
        fails.append("every pile came out of the solve at zero moment; the "
                     "fixture never engages the beams")
    pile_report = _engine_report("fem", xlsx=FEM_PILES_XLSX)
    pile_sec = _member_section(pile_report, "Pile Forces")
    if pile_sec is None:
        fails.append(f"a run carrying {len(pile_profiles)} piles has no Piles "
                     f"subsection: {_titles(pile_report)}")
    else:
        pile_table = next((b for b in pile_sec.blocks if b.kind == "table"), None)
        if pile_table is None:
            fails.append("the Piles subsection carries no table")
        else:
            fails += _member_faults(pile_table, pile_profiles, "Peak moment",
                                    "max_moment")
        pile_figures = [b for b in pile_sec.blocks if b.kind == "figure"]
        if len(pile_figures) != len(pile_profiles):
            fails.append(f"{len(pile_profiles)} piles drew "
                         f"{len(pile_figures)} detail figures; every pile the "
                         f"analysis solved is drawn")
    planned, drawn = _planned_matches(pile_report, "fem", xlsx=FEM_PILES_XLSX)
    if planned != drawn:
        fails.append(f"the piles report planned {planned} figures and built "
                     f"{drawn}")

    # --- neither member, and neither subsection ------------------------------
    plain = _engine_report("fem")
    for title in ("Reinforcement Forces", "Pile Forces"):
        if _member_section(plain, title) is not None:
            fails.append(f"a run carrying no member was given a {title} "
                         f"subsection")
    if _member_section(_build(), "Reinforcement Forces") is not None:
        fails.append("a limit equilibrium report carries a finite element "
                     "Reinforcement subsection")

    # --- the toggles ---------------------------------------------------------
    off = _engine_report("fem", {"fem_reinforcement": False},
                         xlsx=FEM_REINF_XLSX)
    if _member_section(off, "Reinforcement Forces") is not None:
        fails.append("the reinforcement toggle left the subsection standing")
    no_figures = _engine_report("fem", {"fem_reinforcement_figure": False},
                                xlsx=FEM_REINF_XLSX)
    quiet = _member_section(no_figures, "Reinforcement Forces")
    if quiet is None or not [b for b in quiet.blocks if b.kind == "table"]:
        fails.append("switching the detail figures off took the table with them")
    elif [b for b in quiet.blocks if b.kind == "figure"]:
        fails.append("the detail figures were drawn with their toggle off")

    # --- the field the forces are read at ------------------------------------
    #
    # An SSRM run is asked about the mechanism it developed. Hand the same run
    # an at-failure snapshot whose bar forces are twice the converged ones, and
    # the table has to move to it and say which field it read.
    import numpy as np
    snapshot = dict(bundle["solution"])
    snapshot["forces_1d"] = 2.0 * np.asarray(bundle["solution"]["forces_1d"],
                                             dtype=float)
    at_failure = dict(bundle, analysis="ssrm", FS=1.207,
                      failure_solution=snapshot)
    failed = _engine_report("fem", bundle=at_failure, xlsx=FEM_REINF_XLSX)
    sec = _member_section(failed, "Reinforcement Forces")
    if sec is None:
        fails.append("the at-failure run lost its Reinforcement subsection")
    else:
        said = " ".join(b.text for b in sec.blocks if b.kind == "prose")
        if "at failure" not in said:
            fails.append(f"the profiles were read at failure and the section "
                         f"does not say so: {said!r}")
        table = next((b for b in sec.blocks if b.kind == "table"), None)
        force_col = _column(table, "Force") if table else None
        if force_col is not None:
            printed = max(float(r[force_col].replace(",", ""))
                          for r in table.rows)
            if abs(printed - 2.0 * peak) > 0.05 + 0.005 * 2.0 * peak:
                fails.append(f"the at-failure table peaks at {printed:,.1f} "
                             f"where the snapshot carries {2.0 * peak:,.1f}; "
                             f"the converged field was read instead")

    # Mutation: the subsection is absent because the model owns no member, not
    # because an option happened to be off. Hand the builder a member for the
    # model that has none, and the absence rule above has to notice.
    import xslope.report as report_mod
    invented = {"kind": "reinforcement", "index": 1, "label": "invented",
                "field_state": "converged", "peak_T": 0.0, "peak_s": 0.0,
                "peak_utilization": 0.0, "status": "within capacity",
                "units": {}}
    saved = report_mod._detail_profiles
    report_mod._detail_profiles = (
        lambda sd, b, kind: [dict(invented)] if kind == "reinforcement"
        else saved(sd, b, kind))
    try:
        mutated = _engine_report("fem", {"fem_reinforcement_figure": False})
        if _member_section(mutated, "Reinforcement Forces") is None:
            fails.append("a member invented for a model with none printed no "
                         "subsection; the absence rule cannot fail")
    finally:
        report_mod._detail_profiles = saved

    # Mutation: a force column zeroed while the utilization beside it still
    # claims a share of capacity is the defect the table exists to make
    # impossible to hide.
    saved = report_mod._detail_profiles
    report_mod._detail_profiles = (
        lambda sd, b, kind: [dict(p, peak_T=0.0) for p in saved(sd, b, kind)]
        if kind == "reinforcement" else saved(sd, b, kind))
    try:
        zeroed = _engine_report("fem", {"fem_reinforcement_figure": False},
                                xlsx=FEM_REINF_XLSX)
        sec = _member_section(zeroed, "Reinforcement Forces")
        table = next((b for b in (sec.blocks if sec else [])
                      if b.kind == "table"), None)
        if table is None or not _member_faults(table, profiles, "Force",
                                               "peak_T"):
            fails.append("every peak force was zeroed and the table still "
                         "agreed with its profiles; the honesty rule has no "
                         "teeth")
    finally:
        report_mod._detail_profiles = saved
    return fails


def test_engine_sections_follow_their_solutions():
    """Neither section is built without its engine's solution, and each toggle
    removes what it names."""
    fails = []
    from xslope.report import FEM_PANELS, build_report

    # The LEM sample carries neither a seepage nor a finite element solution, and
    # the default report of it has neither section.
    titles = [t for _lvl, t in _build().section_titles()]
    for gone in ("Seepage Analysis", "Deformation and Strength Reduction",
                 "Deformation Analysis"):
        if gone in titles:
            fails.append(f"a report with no seepage or FEM solution carries "
                         f"{gone!r}")

    # And the same model, reported with the engine's own solution absent from
    # the mapping: the option is on, and there is still no section.
    for engine, heading in (("seep", "Seepage Analysis"),
                            ("fem", "Deformation and Strength Reduction")):
        slope_data, _bundle = (_seep_bundle() if engine == "seep"
                               else _fem_bundle())
        tmp = tempfile.mkdtemp(prefix="xslope_absent_")
        with contextlib.redirect_stdout(io.StringIO()):
            report = build_report(slope_data, {}, {"lem": False,
                                                   "pd_figure": False}, tmp)
        if heading in [t for _lvl, t in report.section_titles()]:
            fails.append(f"{heading!r} was built from a solutions mapping that "
                         f"carries no {engine!r}")

    # Each case names what its option takes out: a whole section, every block of
    # one kind, or the figures with these sources — and nothing else goes with it.
    cases = [
        ("seep", {"seep": False}, "Seepage Analysis", None, ()),
        ("seep", {"seep_materials": False}, None, "table", ()),
        ("seep", {"seep_inputs_figure": False}, None, None, ("seep model",)),
        # One option, both conventions: the pair is one statement of the
        # conductivity model and goes or stays together.
        ("seep", {"seep_kr_figure": False}, None, None,
         ("seep kr", "seep kr_head")),
        ("seep", {"seep_mesh_figure": False}, None, None, ("seepage bc1 mesh",)),
        ("seep", {"seep_flownet": False}, None, None, ("seepage bc1",)),
        ("fem", {"fem": False}, "Deformation and Strength Reduction", None, ()),
        # The loads table has to go with the properties table for the section to
        # carry none: this model is loaded, and its loads are an input of the
        # analysis rather than a property of a material.
        ("fem", {"fem_materials": False, "fem_loads": False}, None, "table", ()),
        ("fem", {"fem_loads": False}, "Loads", None, ()),
        ("fem", {"fem_inputs_figure": False}, None, None, ("fem model",)),
        ("fem", {"fem_mesh_figure": False}, None, None, ("fem mesh",)),
        ("fem", {"fem_figure": False}, None, None,
         tuple(f"fem run1 {panel}" for panel, _c, _s in FEM_PANELS)),
    ]
    heads = {"seep": "Seepage Analysis",
             "fem": "Deformation and Strength Reduction"}
    for engine, options, gone_section, gone_kind, gone_figures in cases:
        full = _engine_report(engine)
        off = _engine_report(engine, options)
        if gone_section is not None:
            if gone_section in _titles(off):
                fails.append(f"{options} left {gone_section!r} in the report")
            continue
        if heads[engine] not in _titles(off):
            fails.append(f"{options} removed the {heads[engine]!r} section, not "
                         f"only what it names")
        if gone_kind is not None:
            under = _blocks_under(off, heads[engine], gone_kind)
            was = _blocks_under(full, heads[engine], gone_kind)
            if not was:
                fails.append(f"the {engine} report has no {gone_kind} to remove, "
                             f"so {options} proves nothing")
            if under:
                fails.append(f"{options} left {len(under)} {gone_kind}(s) in the "
                             f"{heads[engine]} section")
        if gone_figures:
            was = {f.source for f in full.figures()}
            now = {f.source for f in off.figures()}
            for source in gone_figures:
                if source not in was:
                    fails.append(f"the {engine} report draws no {source!r}, so "
                                 f"{options} proves nothing")
                if source in now:
                    fails.append(f"{options} left the {source!r} figure standing")
            if now != was - set(gone_figures):
                fails.append(f"{options} changed the figures to {sorted(now)}, "
                             f"not only by dropping {sorted(gone_figures)}")
        planned, drawn = _planned_matches(off, engine, options)
        if planned != drawn:
            fails.append(f"{options} planned {planned} figures and built {drawn}")
    return fails


def _blocks_under(report, heading, kind):
    """Every block of one kind in the section headed ``heading`` and below it."""
    sec = next((s for s in report.sections if s.title == heading), None)
    if sec is None:
        return []
    return [b for _lvl, node in sec.walk() for b in node.blocks if b.kind == kind]


#: A rapid drawdown analysis: two boundary condition sets, two saved solutions.
RAPID_SEEP_XLSX = os.path.join(_REPO, "docs", "lem", "files",
                               "xslope_earth_dam_rapid.xlsx")

#: A model whose saved seepage solution has no mesh beside it to be placed on.
NOMESH_SEEP_XLSX = os.path.join(_REPO, "docs", "seep", "files",
                                "xslope_levee1.xlsx")


def _shipped_flowrate(path):
    """The total flowrate a saved seepage solution records in its own footer."""
    with open(path) as f:
        for line in f:
            if line.startswith("# Total Flowrate:"):
                return float(line.split(":", 1)[1])
    return None


def _sidecar_copy(stem, tmp, meta_edit=None):
    """A copy of a model and its solution sidecars in ``tmp``, for the checks
    that damage one and ask what is made of it. ``meta_edit`` rewrites the FEM
    run metadata where it is given. Returns the copied stem."""
    import glob
    for path in glob.glob(stem + "*"):
        if os.path.isfile(path):
            shutil.copy(path, tmp)
    out = os.path.join(tmp, os.path.basename(stem))
    meta_path = f"{out}_fem_meta.json"
    if meta_edit is not None:
        with open(meta_path) as f:
            meta = json.load(f)
        with open(meta_path, "w") as f:
            json.dump(meta_edit(meta), f)
    return out


def _truncated(path, rows=5):
    """``path`` with its last few rows cut off — a solution that no longer has a
    row per node, which is what a saved solution looks like once its mesh has
    been rebuilt under it."""
    with open(path) as f:
        lines = f.readlines()
    with open(path, "w") as f:
        f.writelines(lines[:-rows])


def test_sidecars_assemble_the_solutions():
    """A model solved once is reported without solving it again: the helper reads
    every engine solution saved beside it back into the mapping the report
    consumes, and nothing it returns was computed here.
    """
    fails = []
    from xslope.report import seep_bundles

    # The dam: one boundary condition set, and the flow the saved solution
    # states is the flow the restored bundle carries.
    _sd, solutions = _restored(SEEP_XLSX)
    if set(solutions) != {"seep"}:
        fails.append(f"the dam restored {sorted(solutions)}, expected seepage "
                     f"alone")
    bundles = seep_bundles(solutions)
    shipped = _shipped_flowrate(f"{os.path.splitext(SEEP_XLSX)[0]}_seep.csv")
    got = bundles[0]["solution"]["flowrate"] if bundles else None
    if shipped is None:
        fails.append("the dam's saved seepage solution records no flowrate, so "
                     "there is nothing to restore it against")
    elif got != shipped:
        fails.append(f"the restored seepage solution flows {got}, and the file "
                     f"beside the model says {shipped}")
    if bundles and set(bundles[0]) != {"seep_data", "solution", "options"}:
        fails.append(f"the restored seepage bundle is {sorted(bundles[0])}, not "
                     f"the shape the report reads")

    # Rapid drawdown: two sets solved, two solutions saved, and they come back in
    # the order the section documents them.
    _sd, solutions = _restored(RAPID_SEEP_XLSX)
    bcs = [b["options"]["bc"] for b in seep_bundles(solutions)]
    if bcs != [1, 2]:
        fails.append(f"the rapid drawdown model restored boundary condition "
                     f"sets {bcs}, expected [1, 2]")

    # The strength reduction run: its factor of safety and what it was.
    _sd, solutions = _restored(FEM_XLSX)
    fem = solutions.get("fem")
    if fem is None:
        fails.append("the strength reduction run saved beside the FEM model was "
                     "not restored")
    else:
        if fem["analysis"] != "ssrm":
            fails.append(f"the restored run calls itself {fem['analysis']!r}, "
                         f"not a strength reduction run")
        if not (1.0 < (fem["FS"] or 0) < 2.0):
            fails.append(f"the restored factor of safety is {fem['FS']}")

    # A model that was never solved has nothing beside it, and asking is not an
    # error — the caller gets an empty mapping and reports what it has.
    notes = []
    _sd, solutions = _restored(REINF_XLSX, notes)
    if solutions != {}:
        fails.append(f"a model with no saved solutions restored "
                     f"{sorted(solutions)}")
    if notes:
        fails.append(f"a model with no saved solutions was faulted for it: "
                     f"{notes}")
    return fails


def test_unusable_sidecars_are_reported_not_raised():
    """A saved solution that no longer fits its model costs its own section and
    nothing else: the engine is left out, the caller is told which file and why,
    and the call returns.
    """
    fails = []

    # Stale: the solution was saved against a mesh that has since been rebuilt,
    # so its node count no longer matches. The message names both counts.
    with tempfile.TemporaryDirectory() as tmp:
        stem = _sidecar_copy(os.path.splitext(FEM_XLSX)[0], tmp)
        _truncated(f"{stem}_fem_nodes.csv")
        notes = []
        _sd, solutions = _restored(f"{stem}.xlsx", notes)
        if "fem" in solutions:
            fails.append("a finite element solution saved against a different "
                         "mesh was restored anyway")
        if not any(os.path.basename(stem) in n and "does not match this mesh" in n
                   for n in notes):
            fails.append(f"a stale finite element solution was passed over "
                         f"without saying why: {notes}")

    # The same for a seepage solution, and the reason is the node count again.
    with tempfile.TemporaryDirectory() as tmp:
        stem = _sidecar_copy(os.path.splitext(SEEP_XLSX)[0], tmp)
        _truncated(f"{stem}_seep.csv")
        notes = []
        _sd, solutions = _restored(f"{stem}.xlsx", notes)
        if "seep" in solutions:
            fails.append("a seepage solution saved against a different mesh was "
                         "restored anyway")
        if not any(f"{os.path.basename(stem)}_seep.csv" in n
                   and "does not match this mesh" in n for n in notes):
            fails.append(f"a stale seepage solution was passed over without "
                         f"saying why: {notes}")

    # A saved solution with no mesh beside it has nothing to be placed on, and
    # that is a different reason, said differently.
    notes = []
    _sd, solutions = _restored(NOMESH_SEEP_XLSX, notes)
    if solutions:
        fails.append(f"a model with no mesh restored {sorted(solutions)}")
    if not any("no mesh" in n for n in notes):
        fails.append(f"a saved solution with no mesh to place it on was passed "
                     f"over without saying why: {notes}")

    # A run saved under the older key still names itself. The finite element
    # companion Studio writes says "analysis"; the ones the benchmark figures
    # were built with say "analysis_type", and a strength reduction run read
    # back as a loaded one is one whose factor of safety the report will not
    # state.
    with tempfile.TemporaryDirectory() as tmp:
        stem = _sidecar_copy(
            os.path.splitext(FEM_XLSX)[0], tmp,
            lambda meta: {("analysis_type" if k == "analysis" else k): v
                          for k, v in meta.items()})
        _sd, solutions = _restored(f"{stem}.xlsx")
        got = (solutions.get("fem") or {}).get("analysis")
        if got != "ssrm":
            fails.append(f"a strength reduction run saved under the older key "
                         f"reads back as {got!r}")
    return fails


# --------------------------------------------------------------------------
# H. nothing is said that is not so
# --------------------------------------------------------------------------

def test_water_prose_is_conditional():
    """A dry model gets one plain statement about water; a model with a pool gets
    the rows that describe it.

    The sample defines no piezometric line, no specified-head boundary and no
    material pore pressure, and the report still explained how its water loads
    were derived. Prose that describes a feature the model does not have is not a
    harmless extra sentence — it is a statement about the analysis that is false.
    """
    fails = []
    from xslope.fileio import load_slope_data
    from xslope.report import build_report, water_features

    slope_data, solutions = _solved()
    feats = water_features(slope_data)
    if feats["any"]:
        return ["the sample model now carries water; the dry-model check has "
                "nothing to see"]

    report = _build()
    # Water is a stability input, so it is stated where the stability analysis
    # is documented, beside the materials whose pore pressures it sets.
    lem = next((s for s in report.sections
                if s.title == "Limit Equilibrium Analysis"), None)
    water = next((s for _l, s in (lem.walk() if lem else [])
                  if s.title == "Materials"), None)
    if water is None:
        return ["the limit equilibrium section has no Materials section"]
    dry = [b.text for b in water.blocks
           if b.kind == "prose" and "dry" in b.text]
    if not dry:
        fails.append(f"a dry model gets no statement that it is analyzed dry: "
                     f"{[b.text for b in water.blocks if b.kind == 'prose']}")
    else:
        text = dry[0]
        for want in ("no groundwater", "dry"):
            if want not in text:
                fails.append(f"the dry-model statement does not say {want!r}: "
                             f"{text!r}")
        if "zero" not in text:
            fails.append("the dry-model statement does not say pore pressures are "
                         "zero")
        # One sentence carries it (Norm: the long inventory of absent features
        # was itself unnecessary).
        if text.count(".") > 1:
            fails.append(f"the dry-model statement is more than one sentence: "
                         f"{text!r}")
    if any(b.kind == "keyvalues" for b in water.blocks):
        fails.append("a dry model was given water rows")

    # No sentence anywhere in a dry report describes water the model lacks.
    banned = ("water surface", "water loads", "piezometric line", "seepage",
              "ponded", "standing water", "water level")
    for block in report.blocks("prose"):
        low = block.text.lower()
        for phrase in banned:
            if phrase in low and "no " + phrase not in low and \
                    "no piezometric lines" not in low and \
                    "no seepage head boundaries" not in low:
                fails.append(f"a dry report says {phrase!r}: {block.text!r}")
    for block in report.blocks("keyvalues"):
        for label, _value in block.items:
            if "water" in label.lower():
                fails.append(f"a dry report carries the key-value row {label!r}")

    # The positive path: a model that really has a pool states it.
    dam = load_slope_data(DAM_XLSX)
    wet = water_features(dam)
    if not wet["any"]:
        fails.append("the dam sample no longer carries water; the positive path "
                     "is untested")
    else:
        with tempfile.TemporaryDirectory() as tmp:
            dam_report = build_report(dam, solutions,
                                      {"pd_figure": False,
                                       "lem_search_figure": False,
                                       "lem_solution_figure": False}, tmp)
        sub = next((s for _l, s in _sections(dam_report)
                    if s.title == "Materials"), None)
        if sub is None:
            fails.append("the dam report has no Materials section")
        else:
            kv = [b for b in sub.blocks if b.kind == "keyvalues"]
            if not kv:
                fails.append("a model with a reservoir got the dry statement")
            else:
                labels = [l for l, _v in kv[0].items]
                if not any("Water surface" in l for l in labels):
                    fails.append(f"the pool is not stated: {labels}")
                if not any(l == "Water loads" for l in labels):
                    fails.append(f"how the water loads are applied is not stated: "
                                 f"{labels}")
    fails += _water_stage_language()
    fails += _water_surface_cites_the_seepage_section()
    return fails


def _water_stage_language():
    """Stage numbering appears only where the analysis has stages.

    "Water surface (stage 1)" on a model solved at one water state is rapid
    drawdown vocabulary applied to something that is not staged: it tells the
    reader there is a stage 2 to compare it against, and there is not.
    """
    fails = []
    from xslope.fileio import load_slope_data
    from xslope.report import water_features, _water_items

    # A single-state model: one pool, held. Its one water surface is not
    # numbered, and no row it carries says "stage" at all.
    for name, xlsx in (("the zoned dam", DAM_XLSX), ("the dam", SEEP_XLSX)):
        sd = load_slope_data(xlsx)
        feats = water_features(sd)
        if len(feats["surfaces"]) != 1:
            fails.append(f"{name} no longer carries exactly one water surface, "
                         f"so the unstaged wording is untested")
            continue
        for seepage_section in (False, True):
            rows = _water_items(sd, feats, seepage_section=seepage_section)
            staged = [(l, v) for l, v in rows if "stage" in f"{l} {v}".lower()]
            if staged:
                fails.append(f"{name} is solved at one water state and its rows "
                             f"still use stage language: {staged}")
            if not any(l == "Water surface" for l, _v in rows):
                fails.append(f"{name} states no plain 'Water surface' row: "
                             f"{[l for l, _v in rows]}")

    # And the staged analyses keep their numbering: the two states of a rapid
    # drawdown are told apart by it, and a report that dropped it would be
    # naming two different water surfaces the same thing.
    for name, xlsx in (("the earth dam rapid drawdown", RAPID_SEEP_XLSX),
                       ("the transient reservoir drawdown",
                        os.path.join(_REPO, "docs", "lem", "files",
                                     "xslope_johnson_res_rapid.xlsx"))):
        sd = load_slope_data(xlsx)
        feats = water_features(sd)
        rows = _water_items(sd, feats, seepage_section=True)
        surfaces = [l for l, _v in rows if l.startswith("Water surface")]
        if not surfaces:
            fails.append(f"{name} states no water surface at all: "
                         f"{[l for l, _v in rows]}")
        for label in surfaces:
            if "stage" not in label.lower():
                fails.append(f"{name} is a staged analysis and its water "
                             f"surface row is unnumbered: {label!r}")
        if len(set(surfaces)) != len(surfaces):
            fails.append(f"{name} labels two water surfaces the same: {surfaces}")
    return fails


def _water_surface_cites_the_seepage_section():
    """Where the report solves the flow, the stability section cites it for the
    water surface instead of restating the head boundaries.

    The elevations of the head boundaries are stated once, in the section that
    imposes them. A second copy under the stability analysis is a number that can
    only ever come to disagree with the first.
    """
    fails = []
    from xslope.report import water_features, SEEPAGE_ANCHOR, section_anchor

    slope_data, _bundle = _seep_bundle()
    feats = water_features(slope_data)
    if not (feats["heads"] and feats["surfaces"]):
        return ["the seepage sample's head boundaries no longer put water above "
                "the ground surface; the citation path is untested"]

    report = _cite_report(SEEP_XLSX, ("spencer",), engines=("seep",))
    titles = [t for _l, t in report.section_titles()]
    if "Seepage Analysis" not in titles or \
            "Limit Equilibrium Analysis" not in titles:
        return [f"the seepage sample no longer reports both engines: {titles}"]

    seep = next(s for s in report.sections if s.title == "Seepage Analysis")
    if seep.anchor != section_anchor(SEEPAGE_ANCHOR):
        fails.append(f"the Seepage Analysis section carries the anchor "
                     f"{seep.anchor!r}, which nothing can cite it by")
    number = next(n for n, _l, s in report.section_numbers() if s is seep)

    lem = next(s for s in report.sections
               if s.title == "Limit Equilibrium Analysis")
    mats = next((s for _l, s in lem.walk() if s.title == "Materials"), None)
    if mats is None:
        return fails + ["the limit equilibrium section has no Materials section"]

    cites = [b for b in mats.blocks
             if b.kind == "prose" and f"Section {number}" in b.text]
    if not cites:
        fails.append(f"the stability analysis does not send the reader to "
                     f"Section {number} for its water surface: "
                     f"{[b.text for b in mats.blocks if b.kind == 'prose']}")
    else:
        targets = [t for _text, t in cites[0].links]
        if f"#{section_anchor(SEEPAGE_ANCHOR)}" not in targets:
            fails.append(f"the reference to the seepage section is not a live "
                         f"cross-reference: {targets}")

    # And the boundary elevations are set down once. The row names the analysis
    # the surface comes from; the numbers stay in the section that imposes them.
    from xslope.water import water_line_for_stage
    stated = water_line_for_stage(slope_data, stage=1)["source"]
    detail = stated[stated.find("(") + 1:stated.rfind(")")] if "(" in stated else ""
    if not detail:
        fails.append(f"the seepage sample's water source names no boundary "
                     f"detail ({stated!r}), so the duplication is untested")
    for block in mats.blocks:
        if block.kind != "keyvalues":
            continue
        for label, value in block.items:
            if detail and detail in value:
                fails.append(f"the stability analysis restates the seepage "
                             f"boundaries the seepage section already gives: "
                             f"{label!r} = {value!r}")
    return fails


def _sections(report):
    """Every section in the report, as ``(level, section)``."""
    out = []
    for s in report.sections:
        out.extend(s.walk())
    return out


def _under(report, heading):
    """The titles of every section under one top-level heading."""
    sec = next((s for s in report.sections if s.title == heading), None)
    return [] if sec is None else [s.title for _l, s in sec.walk()]


def _table_in(report, heading, caption):
    """The table with ``caption`` printed under ``heading``, or None."""
    sec = next((s for s in report.sections if s.title == heading), None)
    if sec is None:
        return None
    for _l, sub in sec.walk():
        for block in sub.blocks:
            if block.kind == "table" and block.caption == caption:
                return block
    return None


#: What each engine's member tables must and must not carry, in the words their
#: headers are built with. Written out here rather than read off the registry the
#: tables are built from: a registry that starts handing the finite element
#: analysis a direction, or the method of slices a Young's modulus, has to fail
#: rather than pass on its own new claim.
#:
#: Each entry is the properties that engine reads and the ones the OTHER engine
#: reads and it does not. The evidence is the code: the limit equilibrium force
#: on a slice is built in ``slice.py`` from the capacity, the pullout envelope,
#: the direction and the application; the finite element element is assembled in
#: ``fem.py`` from the modulus and area, softening to the residual capacity, and
#: reads neither the direction nor whether the force is applied or mobilized.
_MEMBER_COLUMNS = {
    ("reinforcement", "lem"): (("T_max", "L_p1", "L_p2", "Direction",
                                "Applied"), ("T_res", "E", "Area")),
    ("reinforcement", "fem"): (("T_max", "T_res", "L_p1", "L_p2",
                                "E", "Area"), ("Direction", "Applied")),
    ("piles", "lem"): (("H", "θ (deg)", "V_cap", "M_cap", "Applied"),
                       ("E", "I", "Head fixity")),
    ("piles", "fem"): (("V_cap", "M_cap", "E", "Head fixity"),
                       ("H", "θ (deg)", "Applied")),
}


def _member_columns_are_the_engines(table, member, engine, where):
    """One member table carries the properties its own engine reads, and none of
    the ones only the other engine reads."""
    fails = []
    wanted, refused = _MEMBER_COLUMNS[(member, engine)]
    heads = [h.split(" (")[0].strip() for h in table.headers]
    for name in wanted:
        if name.split(" (")[0] not in heads:
            fails.append(f"{where}: the {engine} {member} table does not give "
                         f"{name}, which that engine reads: {table.headers}")
    for name in refused:
        if name.split(" (")[0] in heads:
            fails.append(f"{where}: the {engine} {member} table gives {name}, "
                         f"which only the other engine reads: {table.headers}")
    return fails


def test_members_stand_with_the_engine_that_reads_them():
    """Reinforcement and piles are stated where the analysis that reads them is
    documented, with the properties THAT analysis reads.

    They were in the general description of the model, which said they were a
    property of the section. They are not: they are structure an analysis acts
    on, and the two engines read different things off the same line. The method
    of slices resolves a capacity onto the slice base and needs the pullout
    envelope, the direction and whether the force is applied or mobilized; the
    finite element analysis carries the line as an element and needs the modulus
    and area it is assembled from and the residual it softens to. So each engine
    states its own, and a report of both prints two tables that share a geometry
    and nothing else.

    Each is its own section, present only where the model carries that member:
    one heading for the reinforcement, one for the piles, and neither hiding
    under the other.
    """
    fails = []
    from xslope.report import MEMBER_CAPTIONS, build_report

    caption = MEMBER_CAPTIONS
    slope_data, solutions = _solved()
    report = _build()
    if not slope_data.get("reinforcement_lines"):
        return ["the sample carries no reinforcement; the split is untested"]
    if slope_data.get("pile_lines"):
        return ["the sample now carries piles; the 'absent' half is untested"]
    # Two tables off the same lines need two names in the list of tables.
    if len(set(caption.values())) != len(caption):
        fails.append(f"two member tables share a caption: {caption}")

    # --- not in the general description of the model ---
    for gone in ("Reinforcement", "Piles", "Reinforcement and Piles"):
        if gone in _under(report, "Project Definition"):
            fails.append(f"Project Definition still carries {gone!r}")
    lem_titles = _under(report, "Limit Equilibrium Analysis")
    if "Reinforcement" not in lem_titles:
        fails.append(f"the limit equilibrium section does not state the "
                     f"reinforcement it was run on: {lem_titles}")
    if "Piles" in lem_titles:
        fails.append("a model with no piles was given a Piles section")

    # --- the columns are the engine's own ---
    table = _table_in(report, "Limit Equilibrium Analysis",
                      caption[("reinforcement", "lem")])
    if table is None:
        fails.append("the limit equilibrium section prints no reinforcement table")
    else:
        fails += _member_columns_are_the_engines(
            table, "reinforcement", "lem", "the sample")

    # --- a model with piles and no reinforcement gets the other half ---
    piled = dict(slope_data)
    piled["reinforcement_lines"] = []
    piled["pile_lines"] = load_slope_data_cached(PILES_XLSX)["pile_lines"]
    with tempfile.TemporaryDirectory() as tmp:
        report = build_report(piled, solutions, {"pd_figure": False}, tmp)
    titles = _under(report, "Limit Equilibrium Analysis")
    if "Piles" not in titles:
        fails.append(f"a model with piles got no Piles section: {titles}")
    if "Reinforcement" in titles:
        fails.append("a model with no reinforcement got a Reinforcement section")
    table = _table_in(report, "Limit Equilibrium Analysis",
                      caption[("piles", "lem")])
    if table is None:
        fails.append("the limit equilibrium section prints no piles table")
    else:
        fails += _member_columns_are_the_engines(table, "piles", "lem",
                                                 "a piled model")

    # --- the finite element analysis states its own, off the same lines ---
    for xlsx, member in ((FEM_REINF_XLSX, "reinforcement"),
                         (FEM_PILES_XLSX, "piles")):
        fem_data, bundle = _fem_1d_bundle(xlsx)
        with tempfile.TemporaryDirectory() as tmp:
            report = build_report(fem_data, {"fem": [bundle]},
                                  dict(FAST_FIGURES, pd_figure=False), tmp)
        # One gravity trial is a deformation analysis, not a strength reduction,
        # and the section is headed for what it is.
        head = next((s.title for s in report.sections
                     if s.title.startswith("Deformation")), "")
        titles = _under(report, head)
        wanted = "Reinforcement" if member == "reinforcement" else "Piles"
        if wanted not in titles:
            fails.append(f"{os.path.basename(xlsx)}: the finite element section "
                         f"does not state the {member} it carries: {titles}")
        named = caption[(member, "fem")]
        table = _table_in(report, head, named)
        if table is None:
            fails.append(f"{os.path.basename(xlsx)}: the finite element section "
                         f"prints no {named!r} table")
            continue
        fails += _member_columns_are_the_engines(
            table, member, "fem", os.path.basename(xlsx))
        # A finite-element-only report puts them in the only engine section
        # there is, and nowhere else.
        if _table_in(report, "Project Definition", named) is not None:
            fails.append(f"{os.path.basename(xlsx)}: the {named!r} table is in "
                         f"the Project Definition of a report that has an "
                         f"engine to state it under")

    # --- the two engines' tables really are different ---
    #
    # Both sides of the split assembled from one registry can agree by accident;
    # what makes the split real is that the same line prints different columns
    # under each engine.
    from xslope.report import _reinforcement_table, _Counter
    lem_t = _reinforcement_table(slope_data, _Counter(), "lem")
    fem_t = _reinforcement_table(slope_data, _Counter(), "fem")
    if lem_t is not None and fem_t is not None and lem_t.headers == fem_t.headers:
        fails.append(f"the two engines print the same reinforcement columns, so "
                     f"the split states nothing: {lem_t.headers}")
    return fails


def test_model_checks_default_and_filtering():
    """The model checks are opt-in, and what they carry is scoped to the report.

    Relevance is not guessed from the wording: every preflight rule declares the
    analyses it applies to, so a finding is kept exactly when its own rule applies
    to an analysis the report documents. The rules used here are picked off the
    registry rather than named, so the check follows the registry.
    """
    fails = []
    from xslope.preflight import Finding, PreflightReport, rules
    from xslope.report import DEFAULT_OPTIONS, relevant_findings, report_analyses

    if DEFAULT_OPTIONS["model_checks"] is not False:
        fails.append("the model checks are on by default")

    _slope_data, solutions = _solved()
    if report_analyses(solutions, {"lem": True}) != ["lem"]:
        fails.append(f"a LEM report says it documents "
                     f"{report_analyses(solutions, {'lem': True})}")

    lem_only = next((r for r in rules()
                     if "lem" in r.analyses and "fem" not in r.analyses
                     and "*" not in r.analyses), None)
    fem_only = next((r for r in rules()
                     if "fem" in r.analyses and "lem" not in r.analyses
                     and "*" not in r.analyses), None)
    every = next((r for r in rules() if "*" in r.analyses), None)
    if lem_only is None or fem_only is None or every is None:
        return fails + ["the rule registry no longer offers a LEM-only, a FEM-only "
                        "and an every-analysis rule to test the filter with"]

    findings = [Finding(lem_only.id, "warning", "a limit equilibrium finding"),
                Finding(fem_only.id, "error", "a finite element finding"),
                Finding(every.id, "info", "a finding about the model itself"),
                Finding("no.such.rule", "info", "a finding this build cannot place")]
    kept = [f.rule_id for f in relevant_findings(findings, ["lem"])]
    if fem_only.id in kept:
        fails.append(f"a LEM report kept the FEM-only finding {fem_only.id!r}")
    for want in (lem_only.id, every.id, "no.such.rule"):
        if want not in kept:
            fails.append(f"a LEM report dropped {want!r}")

    # And the filter really reaches the section the report builds.
    report = _build({"model_checks": True,
                     "preflight": PreflightReport(analysis="fem",
                                                  findings=findings)})
    checks = next((s for _l, s in _sections(report) if s.title == "Model Checks"),
                  None)
    if checks is None:
        return fails + ["model_checks=True built no Model Checks section"]
    tables = [b for b in checks.blocks if b.kind == "table"]
    if not tables:
        fails.append("the findings did not reach a table")
    else:
        ids = {r[2] for r in tables[0].rows}
        if fem_only.id in ids:
            fails.append(f"the built section carries the FEM-only finding "
                         f"{fem_only.id!r}")
        if lem_only.id not in ids:
            fails.append(f"the built section dropped {lem_only.id!r}")

    # Nothing relevant left: the section says so in those terms.
    report = _build({"model_checks": True,
                     "preflight": PreflightReport(
                         analysis="fem",
                         findings=[Finding(fem_only.id, "error", "a FEM finding")])})
    checks = next((s for _l, s in _sections(report) if s.title == "Model Checks"),
                  None)
    texts = " ".join(b.text for b in (checks.blocks if checks else [])
                     if b.kind == "prose")
    if "no findings for the analyses in this report" not in texts:
        fails.append(f"a report whose only findings were filtered out does not say "
                     f"so: {texts!r}")
    return fails


def test_title_page_omits_empty_rows():
    """A title-page row with nothing in it is not printed.

    Not every project has a number, an organization or a named author, and a label
    with a blank beside it is a question the title page asks and does not answer.
    """
    fails = []
    from xslope.report import generate_report

    slope_data, solutions = _solved()
    opts = {"input_path": REINF_XLSX, "title": "Bare Title Page",
            "pd_figure": False, "lem_search_figure": False,
            "lem_solution_figure": False, "lem_slice_table": False}
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "bare.docx")
        ok, out = generate_report(slope_data, solutions, opts, out_path)
        if not ok:
            return [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(out_path)
        doc = xml.get("word/document.xml", "")
        for absent in ("Project:", "Author:"):
            if absent in doc:
                fails.append(f"a report with no project number and no author "
                             f"still prints {absent!r}")
        if "Date:" not in doc:
            fails.append("the date row is missing; it always has a value")
        if "Bare Title Page" not in doc:
            fails.append("the title is not on the title page")

    # With the fields filled, the rows are there and read as asked.
    opts = dict(opts, project_number="26-001", author="A. Engineer",
                organization="Example Engineering")
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "full.docx")
        ok, out = generate_report(slope_data, solutions, opts, out_path)
        if not ok:
            return fails + [f"generate_report failed: {out}"]
        _names, xml = _docx_parts(out_path)
        doc = xml.get("word/document.xml", "")
        for want in ("Project:", "Author:", "Date:", "26-001", "A. Engineer"):
            if want not in doc:
                fails.append(f"the filled title page is missing {want!r}")
        if "Project number" in doc:
            fails.append("the title page still labels the row 'Project number'")
    return fails


# --------------------------------------------------------------------------
# J. every figure and every table is cited
#
# For a technical report it is customary to cite every table and every figure
# by number, and a block nothing points at is a block the reader meets without
# being sent to it. Three things are checked, over the option combinations and
# the models that switch the numbered blocks on and off:
#
#   1. EVERY block that prints is cited. The citation has to sit in the block's
#      own section or in a section above it — Bishop's block naming "Table 5"
#      says nothing about Spencer's slice table, and a sentence a page away
#      under another method is not a citation of this one.
#   2. NOTHING is cited that does not print. An option that suppresses a table
#      has to take the sentence that names it with it, or the report sends the
#      reader to a number that is not in the document.
#   3. A citation is a live cross-reference: the phrase carries a link to a
#      bookmark, which report_docx places on the caption line of every numbered
#      block.
#
# The numbers are read off the built tree, never off a literal here: a citation
# is only right if it is the number the counter assigned.
# --------------------------------------------------------------------------

#: A citation, as the prose writes it.
CITATION = re.compile(r"\b(Figure|Table) (\d+)\b")

#: What a section citation writes in place of a number the builder cannot know
#: yet. Nothing that survives the build should still carry it.
_MARK_CHAR = "\ue000"

#: Low-resolution figures: these checks read the tree a build produced, not the
#: pixels, and every combination below draws the full set.
FAST_FIGURES = {"dpi": 60, "figsize": (4.0, 2.5)}

PILES_XLSX = os.path.join(_REPO, "docs", "lem", "files", "xslope_piles.xlsx")

_CITE_REPORTS = {}


def _cite_report(xlsx, methods, options=None, engines=()):
    """A report of ``xlsx`` solved by ``methods``, with the figures drawn.

    ``engines`` names the other engines whose solutions the report documents
    alongside — ``("seep",)`` for a model with a solved flow field. Each is read
    back from the solution shipped beside the model. A case naming no method
    switches the limit equilibrium section off: a report of a seepage run or a
    strength reduction run is a report in its own right.

    Cached on its arguments: several combinations share a model, and each is a
    solve plus a set of plots.
    """
    key = (xlsx, tuple(methods), tuple(sorted((options or {}).items())),
           tuple(engines))
    if key in _CITE_REPORTS:
        return _CITE_REPORTS[key]
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.report import build_report
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(xlsx)
    solutions = {}
    if methods:
        circles = slope_data.get("circles") or []
        surface_kw = ({"circle": circles[0]} if circles
                      else {"non_circ": slope_data.get("non_circ")})
        ok, out = generate_slices(slope_data, num_slices=15, **surface_kw)
        if not ok:
            raise RuntimeError(f"{os.path.basename(xlsx)} produced no slices: {out}")
        df, surface = out[0], out[1]
        bundles = []
        for name in methods:
            work = df.copy()
            with contextlib.redirect_stdout(io.StringIO()):
                results = solve_selected(name, work)
            if not isinstance(results, dict):
                continue
            bundles.append({"slice_df": work, "failure_surface": surface,
                            "results": results, "search": None, "method": name})
        if not bundles:
            raise RuntimeError(f"no method converged on {os.path.basename(xlsx)}")
        # A synthetic search on the first bundle, so the search figure and the
        # sentence that cites it are built by every combination that asks for them.
        bundles[0]["search"] = {
            "kind": "circular" if circles else "noncircular",
            "fs_cache": [{"Xo": 10.0, "Yo": 50.0, "Depth": 0.0,
                          "FS": b["results"]["FS"], "slices": b["slice_df"],
                          "failure_surface": surface,
                          "solver_result": b["results"]} for b in bundles],
            "search_path": [{"x": 10.0, "y": 50.0,
                             "FS": bundles[0]["results"]["FS"]}],
            "circle_cache": None,
        }
        solutions["lem"] = bundles
    for engine in engines:
        _sd, bundle = (_seep_bundle(xlsx) if engine == "seep"
                       else _fem_any_bundle(xlsx))
        solutions[engine] = bundle
    opts = {"input_path": xlsx, "method": list(methods)}
    if not methods:
        opts["lem"] = False
    opts.update(FAST_FIGURES)
    opts.update(options or {})
    tmp = tempfile.mkdtemp(prefix="xslope_cite_")
    with contextlib.redirect_stdout(io.StringIO()):
        report = build_report(slope_data, solutions, opts, tmp)
    _CITE_REPORTS[key] = report
    return report


#: The reports the citation rule is checked on. Each names an option
#: combination or a model feature that puts a numbered block into the tree, or
#: takes one out and has to take its citation with it.
CITATION_CASES = [
    ("the shipped defaults", REINF_XLSX, ("spencer", "bishop"), {}, ()),
    ("the calculations switched off", REINF_XLSX, ("spencer",),
     {"lem_calculations": False}, ()),
    ("the slice table switched off", REINF_XLSX, ("spencer",),
     {"lem_slice_table": False}, ()),
    ("the slice key switched off", REINF_XLSX, ("spencer",),
     {"lem_slice_key": False}, ()),
    ("every figure switched off", REINF_XLSX, ("spencer",),
     {"pd_figure": False, "lem_search_figure": False,
      "lem_solution_figure": False, "lem_slice_key": False}, ()),
    ("the search switched off", REINF_XLSX, ("spencer",),
     {"lem_search": False}, ()),
    ("three methods in detail", REINF_XLSX, ("spencer", "bishop", "oms"), {}, ()),
    ("every method in detail", REINF_XLSX,
     ("oms", "bishop", "janbu", "spencer", "corps", "lowe", "mprice"), {}, ()),
    ("the model checks reported", REINF_XLSX, ("spencer",),
     {"model_checks": True}, ()),
    ("a model carrying piles", PILES_XLSX, ("spencer",), {}, ()),
    ("a model with neither reinforcement nor loads", DAM_XLSX, ("spencer",),
     {}, ()),
    ("a model whose working is refused", PASSIVE_XLSX, ("spencer",), {}, ()),
    # The other two engines: a seepage run reported beside the stability it
    # feeds, a seepage run reported on its own, and a strength reduction run.
    ("a seepage run beside the stability analysis", SEEP_XLSX, ("spencer",),
     {}, ("seep",)),
    ("a seepage run on its own", SEEP_XLSX, (), {}, ("seep",)),
    ("a seepage run with no flow net", SEEP_XLSX, (),
     {"seep_flownet": False}, ("seep",)),
    ("a seepage run with no model figure", SEEP_XLSX, (),
     {"seep_inputs_figure": False}, ("seep",)),
    ("a seepage run with no mesh figure", SEEP_XLSX, (),
     {"seep_mesh_figure": False}, ("seep",)),
    ("a seepage run with no conductivity curves", SEEP_XLSX, (),
     {"seep_kr_figure": False}, ("seep",)),
    ("a strength reduction run", FEM_XLSX, (), {}, ("fem",)),
    ("a strength reduction run with no figures", FEM_XLSX, (),
     {"fem_figure": False}, ("fem",)),
    ("a strength reduction run with no model figure", FEM_XLSX, (),
     {"fem_inputs_figure": False}, ("fem",)),
    ("a strength reduction run with no mesh figure", FEM_XLSX, (),
     {"fem_mesh_figure": False}, ("fem",)),
    # The members a run can carry: each puts a table and its detail figures into
    # the tree, and the reinforcement model has more lines than the figure
    # budget, so its table is cited a second time by the sentence that says so.
    ("a run carrying reinforcement", FEM_REINF_XLSX, (), {}, ("fem",)),
    ("a run carrying piles", FEM_PILES_XLSX, (), {}, ("fem",)),
    ("a run whose member figures are switched off", FEM_REINF_XLSX, (),
     {"fem_reinforcement_figure": False}, ("fem",)),
    ("a strength reduction run with no properties table", FEM_XLSX, (),
     {"fem_materials": False}, ("fem",)),
    # Both engines on one model: the finite element section carries the loads the
    # stability section already presented, and cites that section for them.
    ("both engines on one model", FEM_REINF_XLSX, ("spencer",), {}, ("fem",)),
]


def _numbered_blocks(report):
    """Every figure and table that prints, as ``(kind, number, path, caption)``
    where ``path`` is the tuple of section titles it sits under."""
    out = []

    def walk(node, path):
        here = path + (node.title,)
        for block in node.blocks:
            if block.kind in ("figure", "table"):
                out.append(("Figure" if block.kind == "figure" else "Table",
                            block.number, here, block.caption))
        for child in node.children:
            walk(child, here)

    for node in report.sections:
        walk(node, ())
    return out


#: A citation of a section, by the number Word prints in front of its heading.
SECTION_CITATION = re.compile(r"\bSection (\d+(?:\.\d+)*)\b")


def _section_citations(report):
    """Every citation of a section the prose makes, as ``(number, path, block)``."""
    out = []

    def walk(node, path):
        here = path + (node.title,)
        for block in node.blocks:
            if block.kind == "prose":
                for number in SECTION_CITATION.findall(block.text):
                    out.append((number, here, block))
        for child in node.children:
            walk(child, here)

    for node in report.sections:
        walk(node, ())
    return out


def _citations(report):
    """Every citation the prose makes, as ``(kind, number, path, block)``."""
    out = []

    def walk(node, path):
        here = path + (node.title,)
        for block in node.blocks:
            if block.kind == "prose":
                for kind, number in CITATION.findall(block.text):
                    out.append((kind, int(number), here, block))
        for child in node.children:
            walk(child, here)

    for node in report.sections:
        walk(node, ())
    return out


def _reaches(citing, block_path):
    """Does a citation made under ``citing`` reach a block under
    ``block_path``? Only from the block's own section or one above it."""
    return tuple(citing) == tuple(block_path[:len(citing)])


def test_every_block_is_cited():
    """Every numbered block is cited, from its own section or one above it, and
    nothing is cited that does not print."""
    fails = []

    # What the rule refuses, stated on the shape it is applied to: a sentence in
    # one method's block is not a citation of another method's table, and a
    # sentence in a sibling section is not a citation of what stands here.
    lem = ("Limit Equilibrium Analysis",)
    spencer_table = lem + ("Spencer's Method", "Slice Table")
    for citing, allowed in ((spencer_table, True),
                            (lem + ("Spencer's Method",), True),
                            (lem, True),
                            (lem + ("Spencer's Method", "Calculations"), False),
                            (lem + ("Bishop's Simplified Method", "Slice Table"),
                             False)):
        if _reaches(citing, spencer_table) is not allowed:
            fails.append(f"a citation under {' > '.join(citing)} "
                         f"{'does not reach' if allowed else 'reaches'} a block "
                         f"under {' > '.join(spencer_table)}")

    for label, xlsx, methods, options, engines in CITATION_CASES:
        try:
            report = _cite_report(xlsx, methods, options, engines)
        except Exception as exc:
            fails.append(f"{label}: the report could not be built: {exc!r}")
            continue
        blocks = _numbered_blocks(report)
        cites = _citations(report)
        if not blocks:
            fails.append(f"{label}: the report carries no numbered block")
            continue

        for kind, number, path, caption in blocks:
            if not number:
                fails.append(f"{label}: {kind} {caption!r} printed with no number")
                continue
            reached = [c for c in cites
                       if c[0] == kind and c[1] == number and _reaches(c[2], path)]
            if not reached:
                elsewhere = [" > ".join(c[2]) for c in cites
                             if c[0] == kind and c[1] == number]
                fails.append(
                    f"{label}: {kind} {number} ({caption!r}), under "
                    f"{' > '.join(path)}, is cited by no sentence of its own "
                    f"section or of a section above it"
                    + (f"; it is named only under {elsewhere}" if elsewhere else ""))

        printed = {(k, n) for k, n, _p, _c in blocks}
        for kind, number, path, block in cites:
            if (kind, number) not in printed:
                fails.append(
                    f"{label}: a sentence under {' > '.join(path)} cites {kind} "
                    f"{number}, which this report does not print: {block.text!r}")

        # The numbers run 1..n with no gap, in each series: a citation is only
        # useful if the reader can find the number by counting captions.
        for kind in ("Figure", "Table"):
            got = [n for k, n, _p, _c in blocks if k == kind]
            if got != list(range(1, len(got) + 1)):
                fails.append(f"{label}: the {kind.lower()}s are numbered {got}")

        # A section citation names a section this report carries, and every
        # number a citation had to wait for was filled in.
        carried = {n for n, _lvl, _s in report.section_numbers()}
        for number, path, block in _section_citations(report):
            if number not in carried:
                fails.append(
                    f"{label}: a sentence under {' > '.join(path)} cites Section "
                    f"{number}, which this report does not carry: {block.text!r}")
        for block in report.blocks("prose"):
            if _MARK_CHAR in block.text:
                fails.append(f"{label}: a section citation was never numbered: "
                             f"{block.text!r}")
    return fails


def test_citations_are_cross_references():
    """Every cited number is a live cross-reference in the prose, and the .docx
    carries the bookmark it lands on."""
    fails = []
    from xslope.report import cite, cite_anchor

    phrase, links = cite("Figure", 4)
    if phrase != "Figure 4" or links != [("Figure 4", "#xslope_figure_4")]:
        fails.append(f"cite('Figure', 4) = {(phrase, links)!r}")
    if cite("Table", 0) != ("", []):
        fails.append("a block with no number cites as something")

    report = _cite_report(REINF_XLSX, ("spencer", "bishop"), {})
    for kind, number, path, block in _citations(report):
        targets = [t for text, t in (block.links or [])
                   if text == f"{kind} {number}"]
        if not targets:
            fails.append(f"{kind} {number} is named under {' > '.join(path)} "
                         f"without a link to it: {block.text!r}")
        elif not any(t.startswith("#") for t in targets):
            fails.append(f"{kind} {number} links to {targets!r}, which is not a "
                         f"bookmark in this document")

    # And the bookmarks those links name are really placed, on the caption of
    # the block whose number they carry.
    from xslope.report_docx import render_docx

    tmp = tempfile.mkdtemp(prefix="xslope_cite_docx_")
    path = render_docx(report, os.path.join(tmp, "cited.docx"))
    with zipfile.ZipFile(path) as z:
        doc = z.read("word/document.xml").decode("utf-8")
    for kind, number, _path, _caption in _numbered_blocks(report):
        name = cite_anchor(kind, number)
        if f'w:name="{name}"' not in doc:
            fails.append(f"{kind} {number}'s caption carries no bookmark for a "
                         f"citation to land on")
        elif f'w:anchor="{name}"' not in doc:
            fails.append(f"nothing in the document links to {kind} {number}")

    # A link covers its own phrase and nothing else. A renderer that moves runs
    # into a hyperlink element by hand can take the words beside them too, or
    # nest one link inside another, and both read as a sentence in which half a
    # paragraph is a link to one figure.
    fails += _link_spans_fails(report, doc)
    fails += _section_reference_fails()
    return fails


def _link_spans_fails(report, doc_xml):
    """Every hyperlink in the document is exactly one of the report's own link
    phrases, and no link is inside another."""
    import re
    fails = []
    declared = {text for block in report.blocks("prose")
                for text, _target in (block.links or [])}
    for link in re.findall(r"<w:hyperlink[ >].*?</w:hyperlink>", doc_xml, re.S):
        if link.count("<w:hyperlink") > 1:
            inner = re.findall(r"<w:t[^>]*>([^<]*)</w:t>", link)
            fails.append(f"a hyperlink is nested inside another: {inner!r}")
            continue
        text = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", link))
        if text and text not in declared:
            fails.append(f"a hyperlink covers {text!r}, which is no phrase the "
                         f"report asked to be linked")
    return fails


def _section_reference_fails():
    """A citation of a section is a link AND a field: Word computes the number.

    Checked on a report carrying both engines, which is where one section cites
    another — the finite element analysis carries the loads the limit equilibrium
    section already presented, and points at the section rather than printing the
    table twice.
    """
    import re
    fails = []
    from xslope.report import section_anchor

    report = _cite_report(FEM_REINF_XLSX, ("spencer",), {}, ("fem",))
    cited = _section_citations(report)
    if not cited:
        return ["a report of both engines cites no section, so the section "
                "cross-reference is never exercised"]

    numbers = {sec.anchor: number for number, _lvl, sec in report.section_numbers()}
    for number, path, block in cited:
        targets = [t for text, t in (block.links or [])
                   if text == f"Section {number}"]
        if not targets:
            fails.append(f"Section {number} is named under {' > '.join(path)} "
                         f"without a link to it: {block.text!r}")
        elif not any(t.startswith("#" + section_anchor("")) for t in targets):
            fails.append(f"Section {number} links to {targets!r}, which is no "
                         f"heading of this document")
        elif numbers.get(targets[0][1:]) != number:
            fails.append(f"the citation reads Section {number} but points at "
                         f"{targets[0]}, which is Section "
                         f"{numbers.get(targets[0][1:])!r}")

    from xslope.report_docx import render_docx
    tmp = tempfile.mkdtemp(prefix="xslope_secref_")
    path = render_docx(report, os.path.join(tmp, "sections.docx"))
    with zipfile.ZipFile(path) as z:
        doc = z.read("word/document.xml").decode("utf-8")

    # This report puts three figure links in one sentence, which is where a
    # renderer that moves runs by hand takes the words between them.
    fails += _link_spans_fails(report, doc)

    # Every heading is bookmarked, so a reference written later has a target.
    marked = set(re.findall(r'w:name="(xslope_section_[^"]*)"', doc))
    for _number, _lvl, sec in report.section_numbers():
        if sec.anchor not in marked:
            fails.append(f"the heading {sec.title!r} carries no bookmark")

    # And each citation is a REF field on that bookmark, printing its number,
    # with the number Word will compute already cached in it.
    for number, _path, block in cited:
        anchor = next((t[1:] for text, t in (block.links or [])
                       if text == f"Section {number}"), "")
        if not anchor:
            continue
        instr = f"REF {anchor} \\r"
        if instr not in doc:
            fails.append(f"Section {number} is written as text, not as a "
                         f"cross-reference field: no {instr!r} in the document")
            continue
        if f'w:anchor="{anchor}"' not in doc:
            fails.append(f"the field for Section {number} is in no hyperlink, so "
                         f"clicking the reference goes nowhere")
        start = doc.index(instr)
        after = doc.index('w:fldCharType="separate"', start)
        end = doc.index('w:fldCharType="end"', after)
        cached = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc[after:end]))
        if cached != number:
            fails.append(f"the field for Section {number} carries {cached!r} as "
                         f"its result; a reader sees that before Word updates it")
    return fails


# --------------------------------------------------------------------------
# G. the dialog
# --------------------------------------------------------------------------

def _app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_dialog():
    """The dialog constructs offscreen, and its controls reach the options."""
    fails = []
    _app()
    from PySide6.QtCore import Qt
    from studio.report_dialog import ReportDialog, default_output_path

    _slope_data, solutions = _solved()
    dlg = ReportDialog(slope_data=_slope_data, solutions=solutions,
                       model_path=REINF_XLSX, default_method="bishop")
    try:
        if dlg.output_format() != "docx":
            fails.append(f"the dialog opens on {dlg.output_format()!r}")
        if not dlg.output_path().endswith("_report.docx"):
            fails.append(f"the default output path is {dlg.output_path()!r}")
        # The Methods list is a MULTI-select over every method the solver offers,
        # opening ticked on the one the results view is showing.
        from xslope.report import supported_methods
        offered = [i.data(Qt.UserRole) for i in dlg._method_items()]
        if offered != supported_methods():
            fails.append(f"the methods list offers {offered}, not the solver's "
                         f"own {supported_methods()}")
        if dlg.selected_methods() != ["bishop"]:
            fails.append(f"the methods list opens on {dlg.selected_methods()}, "
                         f"not the results view's method")
        if dlg.options().get("method") != ["bishop"]:
            fails.append(f"the options carry {dlg.options().get('method')!r}")

        # Ticking a second one reports both, in list order.
        items = {i.data(Qt.UserRole): i for i in dlg._method_items()}
        items["spencer"].setCheckState(Qt.Checked)
        got = dlg.selected_methods()
        if got != ["bishop", "spencer"]:
            fails.append(f"two ticked methods came back as {got}")

        # And at least one is always ticked: unticking the last re-ticks it.
        for name in got:
            items[name].setCheckState(Qt.Unchecked)
        if not dlg.selected_methods():
            fails.append("every method could be unticked; a report with no "
                         "method has no results in it")
        items["bishop"].setCheckState(Qt.Checked)
        for name in dlg.selected_methods():
            if name != "bishop":
                items[name].setCheckState(Qt.Unchecked)

        # Only DOCX can be picked in S1; the rest are listed and dimmed.
        for i in range(dlg.format.count()):
            key = dlg.format.itemData(i)
            enabled = dlg.format.model().item(i).isEnabled()
            if (key == "docx") != enabled:
                fails.append(f"format {key!r} is {'enabled' if enabled else 'dimmed'}")

        opts = dlg.options()
        for key in ("traceability", "project_definition", "lem",
                    "lem_slice_table", "lem_materials"):
            if opts.get(key) is not True:
                fails.append(f"{key} is not on by default")
        if opts.get("signature_lines") is not False:
            fails.append("signature lines are on by default")
        # The boxes open on the builder's own defaults — one declaration, not two.
        from xslope.report import DEFAULT_OPTIONS
        if opts.get("model_checks") is not False:
            fails.append("the model checks box opens checked; they are opt-in")
        for key in dlg._items:
            if opts.get(key) is not bool(DEFAULT_OPTIONS.get(key, True)):
                fails.append(f"the {key!r} box opens on {opts.get(key)!r}, not the "
                             f"builder's default "
                             f"{bool(DEFAULT_OPTIONS.get(key, True))!r}")

        # The point-coordinate labels are a row of their own under the model
        # figure, on by default, and the box moves the option it names.
        coords = dlg._items.get("pd_coords")
        if coords is None:
            fails.append("the dialog has no Point coordinates row")
        else:
            if coords.parent() is not dlg._items["project_definition"]:
                fails.append("the Point coordinates row is not under Project "
                             "definition")
            if "Point coordinates" not in coords.text(0):
                fails.append(f"the coordinates row reads {coords.text(0)!r}")
            if opts.get("pd_coords") is not True:
                fails.append("the point coordinates open off; geometry is one of "
                             "the things the model figure is for")
            for state in (Qt.Unchecked, Qt.Checked):
                coords.setCheckState(0, state)
                want = state == Qt.Checked
                if dlg.options().get("pd_coords") is not want:
                    fails.append(f"the coordinates box set to {want} came back "
                                 f"as {dlg.options().get('pd_coords')!r}")
            # And it is the figure's option: it never carries the figure with it.
            coords.setCheckState(0, Qt.Unchecked)
            if dlg.options().get("pd_figure") is not True:
                fails.append("turning the coordinate labels off took the model "
                             "figure with them")
            coords.setCheckState(0, Qt.Checked)

        # A toggle reaches the options, and a parent takes its children.
        dlg._items["lem_slice_table"].setCheckState(0, Qt.Unchecked)
        if dlg.options().get("lem_slice_table") is not False:
            fails.append("unchecking the slice table did not reach the options")
        dlg._items["project_definition"].setCheckState(0, Qt.Unchecked)
        after = dlg.options()
        if after.get("pd_coords") is not False:
            fails.append("a section turned off left its children on")
        if not dlg._items["pd_coords"].isDisabled():
            fails.append("a section turned off left its children live")
        # And took nothing that is not its own: the reinforcement is read by the
        # stability analysis and its box is that analysis's.
        if after.get("lem_members") is not True:
            fails.append("turning Project definition off took the limit "
                         "equilibrium section's member properties with it")

        # And the options it produces really build a report.
        from xslope.report import build_report
        dlg._items["project_definition"].setCheckState(0, Qt.Checked)
        opts = dlg.options()
        opts["lem_search_figure"] = False
        with tempfile.TemporaryDirectory() as tmp:
            report = build_report(_slope_data, solutions, opts, tmp)
        titles = [t for _l, t in report.section_titles()]
        if "Slice Table" in titles:
            fails.append("the dialog's options did not carry the slice table off")
        if "Reinforcement" not in titles:
            fails.append("the dialog's options lost the reinforcement section")
    finally:
        dlg.close()

    # Without a saved project the path still points somewhere writable.
    if not default_output_path(None).endswith(".docx"):
        fails.append("an unsaved project gets no default report path")
    return fails


def test_dialog_settings():
    """What the dialog remembers survives a second construction."""
    fails = []
    _app()
    from PySide6.QtCore import QSettings, Qt
    from studio.report_dialog import SETTINGS_PREFIX, ReportDialog

    _slope_data, solutions = _solved()
    with tempfile.TemporaryDirectory() as tmp:
        ini = os.path.join(tmp, "settings.ini")
        settings = QSettings(ini, QSettings.IniFormat)

        first = ReportDialog(slope_data=_slope_data, solutions=solutions,
                             model_path=REINF_XLSX, settings=settings)
        first.organization.setText("Example Engineering")
        first.author.setText("A. Engineer")
        first.signature_lines.setChecked(True)
        first._items["lem_slice_table"].setCheckState(0, Qt.Unchecked)
        for item in first._method_items():
            if item.data(Qt.UserRole) in ("bishop", "corps"):
                item.setCheckState(Qt.Checked)
            elif item.flags() & Qt.ItemIsUserCheckable:
                item.setCheckState(Qt.Unchecked)
        first.title.setText("A One-Off Title")
        first.remember()
        first.close()

        again = ReportDialog(slope_data=_slope_data, solutions=solutions,
                             model_path=REINF_XLSX, settings=settings)
        try:
            if again.organization.text() != "Example Engineering":
                fails.append(f"the organization was not remembered "
                             f"({again.organization.text()!r})")
            if again.author.text() != "A. Engineer":
                fails.append("the author was not remembered")
            if not again.signature_lines.isChecked():
                fails.append("the signature-lines choice was not remembered")
            if again._items["lem_slice_table"].checkState(0) != Qt.Unchecked:
                fails.append("the content selections were not remembered")
            if again._items["traceability"].checkState(0) != Qt.Checked:
                fails.append("a selection that was left on came back off")
            if again.selected_methods() != ["bishop", "corps"]:
                fails.append(f"the methods were not remembered "
                             f"({again.selected_methods()})")
            # The project's own fields are NOT carried between projects.
            if again.title.text() == "A One-Off Title":
                fails.append("the project title was remembered; it belongs to the "
                             "project, not to the user")
        finally:
            again.close()

        # And it really went to disk, under this dialog's own prefix.
        stored = QSettings(ini, QSettings.IniFormat)
        if stored.value(SETTINGS_PREFIX + "author") != "A. Engineer":
            fails.append("the remembered fields are not under the report prefix")

        # A dialog with no settings store remembers nothing and does not fail.
        plain = ReportDialog(slope_data=_slope_data, solutions=solutions,
                             model_path=REINF_XLSX)
        try:
            if plain.author.text():
                fails.append("a dialog with no settings pre-filled the author")
            plain.remember()
        finally:
            plain.close()
    return fails


def test_open_output():
    """The finished document is opened, and the formats offered are the ones
    that exist.

    A LaTeX row that could never be picked was a promise on the screen where the
    choice is made; the report is a Word document, with PDF still to come.
    """
    fails = []
    _app()
    from studio import report_dialog

    formats = [key for key, *_rest in report_dialog.FORMATS]
    if "latex" in formats:
        fails.append(f"LaTeX is still offered: {formats}")
    if formats != ["docx", "pdf"]:
        fails.append(f"the formats offered are {formats}, not ['docx', 'pdf']")

    opened = []
    real = report_dialog.QDesktopServices.openUrl
    report_dialog.QDesktopServices.openUrl = lambda url: opened.append(url.toLocalFile())
    try:
        with tempfile.TemporaryDirectory() as tmp:
            docx = os.path.join(tmp, "r.docx")
            open(docx, "wb").close()
            if report_dialog.open_output(docx, "docx") != "document":
                fails.append("a .docx did not open as a document")
            if opened[-1] != docx:
                fails.append(f"the document opened was {opened[-1]!r}")
    finally:
        report_dialog.QDesktopServices.openUrl = real
    return fails


#: What a runner check builds: two methods, so the figures the progress bar
#: counts are several and two of them are the same figure of different methods,
#: at a resolution that costs nothing to render.
_RUNNER_OPTIONS = {"input_path": REINF_XLSX, "title": "Sample Levee",
                   "method": ["spencer", "bishop"], "lem_search_figure": False,
                   "dpi": 60, "figsize": (4.0, 2.5)}


def _run_report_runner(path, options=None, cancel_after=None, fmt=None):
    """Drive a :class:`studio.runners.ReportRunner` to completion and collect
    what it emitted.

    ``run()`` is called on this thread rather than started as a thread: the
    signals then arrive as direct calls, so the check reads them in order
    without an event loop. What is being checked is what the worker emits, which
    is the same either way.

    ``cancel_after`` cancels the run from inside the progress handler once that
    many steps have arrived — which is what the Cancel button does, at the point
    in the build where it does it. The Word finish is never asked for.
    """
    from studio.runners import ReportRunner

    slope_data, solutions = _solved()
    opts = dict(_RUNNER_OPTIONS)
    opts.update(options or {})
    runner = ReportRunner(slope_data, solutions, opts, path, fmt=fmt,
                          finalize=False)
    seen = {"progress": [], "ok": [], "failed": [], "cancelled": []}

    def on_progress(done, total, label):
        seen["progress"].append((done, total, label))
        if cancel_after is not None and len(seen["progress"]) >= cancel_after:
            runner.cancel()

    runner.progress.connect(on_progress)
    runner.succeeded.connect(lambda out: seen["ok"].append(out))
    runner.failed.connect(lambda msg: seen["failed"].append(msg))
    runner.cancelled.connect(lambda: seen["cancelled"].append(True))
    runner.run()
    return seen


def test_report_runner_progress():
    """The report builds on a worker thread, counting the figures by name.

    The bar is determinate over the figures plus the document itself, every
    figure the build renders is announced with its own label, and the count the
    bar is ranged on is the count the builder plans — the same agreement
    ``planned_figures`` promises, arriving where the user can see it.
    """
    fails = []
    _app()
    import matplotlib
    matplotlib.use("Agg")
    from studio.runners import (REPORT_FINALIZE_LABEL, REPORT_WRITE_LABEL,
                                REPORT_WRITE_STEPS)
    from xslope.report import planned_figures, resolve_options

    slope_data, solutions = _solved()
    planned = planned_figures(slope_data, solutions,
                              resolve_options(dict(_RUNNER_OPTIONS)))
    if planned < 2:
        fails.append(f"the runner check plans {planned} figures; it is meant to "
                     f"exercise several")

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "runner.docx")
        seen = _run_report_runner(path)

        if seen["failed"] or seen["cancelled"]:
            return fails + [f"the run did not succeed: {seen['failed']} "
                            f"{seen['cancelled']}"]
        if len(seen["ok"]) != 1:
            fails.append(f"succeeded fired {len(seen['ok'])} times")
        out = seen["ok"][0] if seen["ok"] else {}
        if not os.path.isfile(path):
            fails.append("the runner emitted success with no document written")
        if out.get("path") != path:
            fails.append(f"the result names {out.get('path')!r}")
        if out.get("finalized") is not False:
            fails.append("a run that was not asked to finish in Word says it did")

        steps = seen["progress"]
        total = planned + REPORT_WRITE_STEPS
        if {t for _d, t, _l in steps} != {total}:
            fails.append(f"the bar was ranged on {sorted({t for _d, t, _l in steps})}, "
                         f"not {total} (= {planned} figures + the document)")
        # The figure steps: one per figure, numbered 1..N in order.
        figure_steps = [(d, l) for d, _t, l in steps
                        if l not in (REPORT_WRITE_LABEL, "rendering the figures")]
        if [d for d, _l in figure_steps] != list(range(1, planned + 1)):
            fails.append(f"the figure steps are {[d for d, _l in figure_steps]}, "
                         f"not 1..{planned}")
        # And they are the figures that were built: the same count, each named.
        drawn = out.get("figures") or []
        if len(figure_steps) != len(drawn):
            fails.append(f"{len(figure_steps)} steps were announced and "
                         f"{len(drawn)} figures reached the document")
        if any(not l.strip() for _d, l in figure_steps):
            fails.append("a figure was announced with no label")
        labels = [l for _d, l in figure_steps]
        # Two methods mean the per-method figures are announced per method, so
        # the labels name which one is rendering rather than repeating.
        for name in ("Spencer", "Bishop"):
            if not any(name in l for l in labels):
                fails.append(f"no step names {name}: {labels}")
        if len(set(labels)) != len(labels):
            fails.append(f"two steps carry the same label: {labels}")

        # The document is a step of its own, named, and the bar reaches the end
        # on it. Word was not asked for and is not announced.
        write = [(d, t) for d, t, l in steps if l == REPORT_WRITE_LABEL]
        if not write:
            fails.append("the document write is not a phase of its own")
        elif write[-1] != (total, total):
            fails.append(f"the bar ends at {write[-1]}, not ({total}, {total})")
        if any(l == REPORT_FINALIZE_LABEL for _d, _t, l in steps):
            fails.append("a run that was not asked to finish in Word announced it")
    return fails


def test_report_runner_cancel():
    """Cancel stops the build at the next figure, and leaves nothing behind.

    The builder takes no cancel argument; what it has is a progress callback it
    calls before each figure, guarded by ``except Exception`` so that a broken
    progress line can never cost a report. ``ReportCancelled`` is a
    ``BaseException`` for exactly that reason — it is the one thing that guard
    lets through — and the builder's own ``finally`` removes the figures on the
    way out. Widening either guard would turn Cancel into a button that does
    nothing, so the seam itself is checked here, not just the button.
    """
    fails = []
    _app()
    import matplotlib
    matplotlib.use("Agg")
    import glob
    from studio.runners import ReportCancelled

    if isinstance(ReportCancelled(), Exception):
        fails.append("ReportCancelled is an Exception; the builder's progress "
                     "guard would swallow it and the run would not stop")

    # The builder renders into a temporary directory of its own and removes it
    # on the way out; a cancelled build must do the same. Counted as a
    # difference, because other checks here keep their figure directories.
    figures = os.path.join(tempfile.gettempdir(), "xslope_report_*")
    before = set(glob.glob(figures))

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "cancelled.docx")
        # Two steps in: the first figure has been announced, so the next one the
        # builder announces is where the run unwinds.
        seen = _run_report_runner(path, cancel_after=2)
        if not seen["cancelled"]:
            fails.append("cancelling mid-build did not report a cancelled run")
        if seen["ok"]:
            fails.append("a cancelled run reported success")
        if seen["failed"]:
            fails.append(f"a cancelled run reported a failure: {seen['failed']}")
        if os.path.exists(path):
            fails.append("a cancelled run left a half-written document")
        # It stopped where it was cancelled rather than running to the end.
        from xslope.report import planned_figures, resolve_options
        slope_data, solutions = _solved()
        planned = planned_figures(slope_data, solutions,
                                  resolve_options(dict(_RUNNER_OPTIONS)))
        if len(seen["progress"]) > planned:
            fails.append(f"the cancelled run announced {len(seen['progress'])} "
                         f"steps of {planned} figures; it did not stop")

    leaked = set(glob.glob(figures)) - before
    if leaked:
        fails.append(f"a cancelled build left its figures behind: {sorted(leaked)}")
    return fails


def test_report_runner_failure():
    """A build that cannot produce a document says so, in the builder's own
    words, instead of succeeding quietly or taking the window down with it."""
    fails = []
    _app()
    import matplotlib
    matplotlib.use("Agg")

    # A format that does not exist yet: refused before anything is rendered.
    seen = _run_report_runner("nowhere.pdf", fmt="pdf")
    if seen["ok"]:
        fails.append("a report in an unavailable format reported success")
    if not seen["failed"]:
        fails.append("an unavailable format produced no failure message")
    elif "not available" not in seen["failed"][0]:
        fails.append(f"the failure reads {seen['failed'][0]!r}")

    # A document that cannot be written: the figures render, the write fails,
    # and the message is the one the builder produced.
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "no_such_folder", "r.docx")
        with contextlib.redirect_stdout(io.StringIO()):
            seen = _run_report_runner(
                path, options={"lem": False, "project_definition": False,
                               "seep": False, "fem": False})
        if seen["ok"]:
            fails.append("a report that could not be written reported success")
        if not seen["failed"]:
            fails.append("an unwritable path produced no failure message")
        elif "could not be written" not in seen["failed"][0]:
            fails.append(f"the write failure reads {seen['failed'][0]!r}")
        if not seen["progress"]:
            fails.append("a failing run showed no progress at all")
    return fails


def test_report_runs_off_the_gui_thread():
    """Generating a report never blocks the window.

    The window hands the build to the runner and waits for its signals: no
    ``generate_report`` on the GUI thread, and no wait cursor standing in for a
    progress bar. While it runs, the report action is dimmed with the rest of
    the run actions, so a second report cannot be started over the first, and
    Cancel reaches the report runner like any other.
    """
    fails = []
    import inspect
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio import report_dialog
    from studio.main_window import MainWindow
    from studio.runners import ReportRunner

    src = inspect.getsource(MainWindow.generate_report)
    if "ReportRunner(" not in src:
        fails.append("the window does not hand the report to the runner")
    if "generate_report as _generate" in src or "_generate(" in src:
        fails.append("the window still builds the report on the GUI thread")
    if "setOverrideCursor" in src:
        fails.append("the report still puts a wait cursor up; the window stays "
                     "live now")
    # The dialog collects options; it does not build reports either.
    dialog_src = inspect.getsource(report_dialog)
    if "generate_report(" in dialog_src:
        fails.append("the dialog builds the report itself")
    # The Word finish the runner calls has no Qt in it.
    finish_src = inspect.getsource(report_dialog.word_finish)
    for qt in ("QApplication", "setOverrideCursor", "_status_line"):
        if qt in finish_src:
            fails.append(f"word_finish touches {qt} off the GUI thread")

    slope_data, solutions = _solved()
    mw = MainWindow()
    try:
        mw.doc.slope_data = slope_data
        mw.doc.results["lem_solution"] = solutions["lem"][0]
        mw._last_lem_opts = {"method": "spencer"}
        mw._update_run_actions()
        if not mw.act_report.isEnabled():
            fails.append("the report action is dimmed with a solved model")

        # A run in flight dims it, and the Run action with it.
        mw._report_runner = ReportRunner(slope_data, solutions, {}, "x.docx",
                                         finalize=False)
        mw._update_run_actions()
        if mw.act_report.isEnabled():
            fails.append("a second report can be started over the first")
        if mw.act_run.isEnabled():
            fails.append("an analysis can be started while the report builds")
        # Cancel goes to the report runner — the one in flight, which here is
        # this one and not started, so it says so.
        cancelled = []
        mw._report_runner.cancel = lambda: cancelled.append(True)
        mw._report_runner.isRunning = lambda: True
        mw.cancel_btn.setVisible(True)
        mw._cancel_run()
        if not cancelled:
            fails.append("Cancel does not reach the report runner")
        # Word's stretch is indeterminate and takes Cancel out with it.
        mw.cancel_btn.setEnabled(True)
        mw._on_report_progress(0, -1, "finalizing…")
        if mw.cancel_btn.isEnabled():
            fails.append("Cancel is still live while Word holds the document")

        mw._report_runner = None
        mw._update_run_actions()
        if not mw.act_report.isEnabled():
            fails.append("the report action stayed dimmed after the run")
    finally:
        mw._report_runner = None
        mw.close()
    return fails


def _noncircular_solutions():
    """The cached solve, presented as a non-circular surface.

    A slice table states its own surface family through the circle radius it
    carries, and blanking that is what makes the same slices non-circular — so
    this costs no second solve and exercises the path the dialog and the report
    both read (``surface_family``).
    """
    slope_data, solutions = _solved()
    bundle = dict(solutions["lem"][0])
    df = bundle["slice_df"].copy()
    df["r"] = float("nan")
    bundle["slice_df"] = df
    bundle["search"] = None
    return dict(slope_data, circular=False), {"lem": [bundle]}


def test_noncircular_dims_the_moment_methods():
    """On a non-circular surface the circular-only methods are dimmed, and the
    list still opens on something the surface can take.

    A moment method needs a center of rotation. Offering one on a non-circular
    surface invites a report section that can only say no, and the ways it goes
    wrong are quiet: the results view may be SHOWING a circular-only method when
    the dialog opens, and a remembered selection may name nothing this surface
    can run. Both have to land on a method that works, because the list always
    keeps one ticked.
    """
    fails = []
    _app()
    from PySide6.QtCore import QSettings, Qt
    from studio.report_dialog import SETTINGS_PREFIX, ReportDialog
    from xslope.preflight import method_surface_reason
    from xslope.report import supported_methods, surface_family

    slope_data, solutions = _noncircular_solutions()
    if surface_family(slope_data, solutions) != "noncircular":
        return ["the fixture is not read as a non-circular surface"]

    circular_only = [n for n in supported_methods()
                     if method_surface_reason(n, "noncircular")]
    if not circular_only:
        return ["no method is circular-only; the check proves nothing"]

    dlg = ReportDialog(slope_data=slope_data, solutions=solutions,
                       model_path=REINF_XLSX, default_method="spencer")
    try:
        for item in dlg._method_items():
            name = item.data(Qt.UserRole)
            checkable = bool(item.flags() & Qt.ItemIsUserCheckable)
            if (name in circular_only) == checkable:
                fails.append(f"{name} is {'live' if checkable else 'dimmed'} on a "
                             f"non-circular surface")
            if name in circular_only:
                if item.checkState() == Qt.Checked:
                    fails.append(f"{name} opens ticked but cannot run")
                if not item.toolTip().startswith(
                        method_surface_reason(name, "noncircular")[:30]):
                    fails.append(f"{name} is dimmed without saying why: "
                                 f"{item.toolTip()!r}")
        # The same model, circular: the very same methods are live again, so the
        # dimming is the surface's doing and not a permanently disabled row.
        circular = ReportDialog(slope_data=dict(slope_data, circular=True),
                                solutions=_solved()[1], model_path=REINF_XLSX)
        try:
            live = [i.data(Qt.UserRole) for i in circular._method_items()
                    if i.flags() & Qt.ItemIsUserCheckable]
            if live != supported_methods():
                fails.append(f"a circular surface dims {set(supported_methods()) - set(live)}")
        finally:
            circular.close()
    finally:
        dlg.close()

    # The results view showing a circular-only method must not leave the dialog
    # opening on one — nor on an arbitrary method, when one that was RUN works.
    dlg = ReportDialog(slope_data=slope_data, solutions=solutions,
                       model_path=REINF_XLSX, default_method="bishop")
    try:
        if dlg.selected_methods() != ["spencer"]:
            fails.append(f"with the view on Bishop and Spencer run, the list "
                         f"opened on {dlg.selected_methods()}")
    finally:
        dlg.close()

    # And a remembered selection this surface cannot take falls back rather than
    # leaving the dialog with nothing ticked.
    with tempfile.TemporaryDirectory() as tmp:
        settings = QSettings(os.path.join(tmp, "s.ini"), QSettings.IniFormat)
        settings.setValue(SETTINGS_PREFIX + "methods", list(circular_only))
        dlg = ReportDialog(slope_data=slope_data, solutions=solutions,
                           model_path=REINF_XLSX, default_method="spencer",
                           settings=settings)
        try:
            got = dlg.selected_methods()
            if not got:
                fails.append("a remembered circular-only selection left the "
                             "dialog with no method ticked")
            if set(got) & set(circular_only):
                fails.append(f"the dialog restored {got}, which this surface "
                             f"cannot run")
        finally:
            dlg.close()

    # Asked for one anyway, the report says so under its own heading rather than
    # dropping the block — a missing section reads as an answer withheld.
    from xslope.report import build_report, method_label
    with tempfile.TemporaryDirectory() as tmp:
        report = build_report(slope_data, solutions,
                              {"method": [circular_only[0], "spencer"],
                               "pd_figure": False, "lem_search_figure": False,
                               "lem_solution_figure": False,
                               "lem_slice_key": False}, tmp)
    titles = [t for _l, t in report.section_titles()]
    head = method_label(circular_only[0])
    if head not in titles:
        fails.append(f"{circular_only[0]} got no section at all: {titles}")
    else:
        said = " ".join(b.text for b in report.blocks("prose"))
        if "cannot be used with a non-circular surface" not in said:
            fails.append(f"the {circular_only[0]} block does not say why it is "
                         f"not reported")
        if [t for t in report.tables()
                if t.landscape and method_label(circular_only[0]) in t.caption]:
            fails.append(f"{circular_only[0]} got a slice table on a surface it "
                         f"cannot run on")
    return fails


def test_slice_numbers_display_option():
    """Studio's LEM solution view can label the slices, and the toggle really
    reaches the plot.

    The same labels the report's slice key carries. A display checkbox wired to
    nothing looks exactly like one that works, so this follows the value from the
    panel through the display-options dict to the keyword the canvas hands the
    plotting function.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio import canvas as canvas_mod
    from studio.display_panels import SolutionDisplayPanel

    panel = SolutionDisplayPanel()
    try:
        if "slice_numbers" not in panel.options():
            return ["the LEM solution display panel has no Slice numbers option"]
        if panel.options()["slice_numbers"] is not False:
            fails.append("slice numbers are on by default")
        box = panel._boxes["slice_numbers"]
        if box.text() != "Slice numbers":
            fails.append(f"the checkbox reads {box.text()!r}")

        # The canvas defers its raster until it is visible, so the stored draw
        # function is called here with a figure of its own — the same call the
        # canvas makes, without needing a screen.
        from matplotlib.figure import Figure as MplFigure

        seen = []
        real = canvas_mod.plot_solution
        canvas_mod.plot_solution = lambda *a, **kw: seen.append(kw)
        try:
            slope_data, solutions = _solved()
            bundle = solutions["lem"][0]
            view = canvas_mod.MplCanvas()
            for state in (True, False):
                box.setChecked(state)
                view.render_solution(slope_data, bundle["slice_df"],
                                     bundle["failure_surface"], bundle["results"],
                                     opts=panel.options())
                view._draw_fn(MplFigure())
            view.deleteLater()
        finally:
            canvas_mod.plot_solution = real

        got = [kw.get("slice_numbers") for kw in seen]
        if got != [True, False]:
            fails.append(f"the canvas was asked for slice_numbers={got}, not "
                         f"[True, False] — the toggle does not reach the plot")
    finally:
        panel.deleteLater()
    return fails


def test_main_window_action():
    """The menu item and toolbar button exist, and are gated on a solution."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio.main_window import MainWindow

    from PySide6.QtWidgets import QMenu, QToolBar

    _slope_data, solutions = _solved()
    mw = MainWindow()
    try:
        act = getattr(mw, "act_report", None)
        if act is None:
            return ["there is no Generate Report action"]
        if "Report" not in act.text():
            fails.append(f"the action reads {act.text()!r}")
        menus = [m for m in mw.findChildren(QMenu) if act in m.actions()]
        if not any("File" in m.title() for m in menus):
            fails.append(f"the action is not on the File menu "
                         f"(found on {[m.title() for m in menus]})")
        if not any(act in tb.actions() for tb in mw.findChildren(QToolBar)):
            fails.append("the action is not on a toolbar")
        if act.isEnabled():
            fails.append("the action is live with nothing solved")
        if "Run an analysis" not in act.toolTip():
            fails.append(f"the dimmed action gives no reason: {act.toolTip()!r}")

        mw.doc.slope_data = _slope_data
        mw.doc.results["lem_solution"] = solutions["lem"][0]
        mw._last_lem_opts = {"method": "bishop"}
        mw._update_run_actions()
        if not act.isEnabled():
            fails.append("the action stays dimmed with a solved model")
        got = mw.report_solutions()
        if not got.get("lem"):
            fails.append("the window offers no LEM solution to the report")
        elif got["lem"][0].get("method") != "bishop":
            fails.append(f"the solution carries method "
                         f"{got['lem'][0].get('method')!r}, not the run's")
    finally:
        mw.close()
    return fails


def _handed_over(slope_data, results):
    """What a main window with ``results`` in its document hands the report."""
    from studio.main_window import MainWindow
    mw = MainWindow()
    try:
        mw.doc.slope_data = slope_data
        mw.doc.results.update(results)
        return mw.report_solutions()
    finally:
        mw.close()


def test_report_solutions_carry_every_engine():
    """Every engine the session has solved reaches the report.

    Studio runs three: limit equilibrium, seepage, and the finite element
    analysis. The report has a section for each, and each reads its engine's own
    bundle — so a window that hands over the LEM solution alone documents a
    seepage run and a strength reduction run nowhere, whatever was solved in
    them.

    What the sidecar says was run is part of it. The metadata Studio writes
    beside a solution names it ``analysis``, the metadata the benchmark figures
    were built with names it ``analysis_type``, and a restored strength reduction
    run that arrives as neither is documented as a deformation analysis with no
    factor of safety — which is a run that reached 1.345 reported as a run that
    reached nothing.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio.main_window import MainWindow
    from xslope.report import (build_report, fem_bundles, lem_bundles,
                               seep_bundles)

    slope_data, solutions = _solved()
    seep_data, seep = _seep_bundle()
    fem_data, fem = _fem_bundle()

    mw = MainWindow()
    try:
        mw.doc.slope_data = slope_data
        mw.doc.results["lem_solution"] = solutions["lem"][0]
        mw._last_lem_opts = {"method": "bishop"}
        # Two boundary condition sets, put in out of order: the report documents
        # them in the order it is handed them.
        mw.doc.results["seep_solutions"] = {
            2: dict(seep, options={"bc": 2}), 1: seep}
        mw.doc.results["fem_solution"] = fem

        got = mw.report_solutions()
        for engine, bundles in (("lem", lem_bundles(got)),
                                ("seep", seep_bundles(got)),
                                ("fem", fem_bundles(got))):
            if not bundles:
                fails.append(f"a session that solved {engine} hands the report "
                             f"nothing for it: {sorted(got)}")
        if len(seep_bundles(got)) != 2:
            fails.append(f"two boundary condition sets were solved and "
                         f"{len(seep_bundles(got))} reached the report")
        bcs = [b.get("options", {}).get("bc") for b in seep_bundles(got)]
        if bcs != sorted(b for b in bcs if b is not None):
            fails.append(f"the seepage sets reach the report in the order {bcs}")

        # A solution from any engine is a report: the action no longer waits on
        # a limit equilibrium run.
        for key in ("lem_solution", "seep_solutions"):
            mw.doc.results.pop(key, None)
        mw._update_run_actions()
        if not mw.act_report.isEnabled():
            fails.append("a strength reduction run alone leaves Generate Report "
                         "dimmed, and its section is the report")
    finally:
        mw.close()

    # And what is handed over is what the sections are built from.
    for engine, model, xlsx, results, title in (
            ("seep", seep_data, SEEP_XLSX, {"seep_solutions": {1: seep}},
             "Seepage Analysis"),
            ("fem", fem_data, FEM_XLSX, {"fem_solution": fem},
             "Deformation and Strength Reduction")):
        opts = {"input_path": xlsx, "lem": False, "pd_figure": False}
        opts.update(FAST_FIGURES)
        tmp = tempfile.mkdtemp(prefix=f"xslope_studio_{engine}_")
        with contextlib.redirect_stdout(io.StringIO()):
            built = build_report(model, _handed_over(model, results), opts, tmp)
        if title not in _titles(built):
            fails.append(f"the {engine} bundle Studio hands over builds "
                         f"{_titles(built)}, with no {title!r} section")

    # The restored run: the sidecar of the sample SSRM model, rewritten to name
    # what was run the way the benchmark builder names it.
    from xslope.fileio import load_slope_data
    with tempfile.TemporaryDirectory() as tmp:
        stem = _sidecar_copy(
            os.path.splitext(FEM_XLSX)[0], tmp,
            lambda meta: {("analysis_type" if k == "analysis" else k): v
                          for k, v in meta.items()})
        mw = MainWindow()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mw.doc.slope_data = load_slope_data(f"{stem}.xlsx")
                mw._restore_fem_sidecar(mw.doc.slope_data["mesh"], stem)
            restored = mw.doc.results.get("fem_solution")
            if not restored:
                fails.append("the FEM sidecar beside the sample model did not "
                             "restore, so what it says was run is untested")
            elif str(restored.get("analysis")) != "ssrm":
                fails.append(f"a restored strength reduction run says it was a "
                             f"{restored.get('analysis')!r} analysis, and its "
                             f"factor of safety {restored.get('FS')} goes "
                             f"unstated")
        finally:
            mw.close()
    return fails


CHECKS = [
    ("the content tree and its section order", test_tree),
    ("the tables carry the model's numbers", test_tables_carry_the_model),
    ("every figure is rendered to a file", test_figures_rendered),
    ("the traceability stamp", test_traceability),
    ("the content toggles remove what they name", test_toggles),
    ("the method picker drives the detail only", test_method_picker),
    ("the summary compares only what is documented",
     test_fs_table_compares_only_what_was_documented),
    ("the summary says where its numbers came from",
     test_the_fs_summary_says_where_its_numbers_came_from),
    ("one full detail block per method", test_multi_method_detail),
    ("the summary reaches the document", test_fs_summary_reaches_the_document),
    ("each engine presents the loads it applies",
     test_each_engine_presents_the_loads_it_applies),
    ("critical is a word a search earns",
     test_critical_is_a_word_a_search_earns),
    ("the stability section opens on its own model", test_lem_inputs_figure),
    ("the search plot is read for the engineer",
     test_search_figure_is_read_for_the_engineer),
    ("the slice key stands before its table", test_slice_key_figure),
    ("the figures are counted for the caller", test_figure_progress_counts),
    ("the seepage section", test_seep_section),
    ("the strength reduction section", test_fem_section),
    ("reinforcement and piles in the finite element section",
     test_fem_members_are_reported),
    ("every member detail figure is readable",
     test_member_detail_figures_are_readable),
    ("each engine's section follows its solution",
     test_engine_sections_follow_their_solutions),
    ("a solved model's companions are the solutions",
     test_sidecars_assemble_the_solutions),
    ("an unusable companion is reported, not raised",
     test_unusable_sidecars_are_reported_not_raised),
    ("the water prose follows the model", test_water_prose_is_conditional),
    ("members stand with the engine that reads them",
     test_members_stand_with_the_engine_that_reads_them),
    ("the model checks are opt-in and scoped",
     test_model_checks_default_and_filtering),
    ("an empty title-page field prints no row", test_title_page_omits_empty_rows),
    ("the .docx and its structure", test_docx),
    ("the tables are fitted to their content", test_table_geometry),
    ("section breaks take no room", test_section_breaks_take_no_room),
    ("every paragraph is on a page of the report",
     test_every_block_reaches_the_page),
    ("the contents page lists the report", test_contents_page),
    ("Word numbers the headings", test_heading_numbers),
    ("the report writes one file", test_report_writes_one_file),
    ("the shipped template is reproducible", test_docx_template),
    ("the slice-column registry", test_column_registry),
    ("every force is in every equation or says why not",
     test_force_term_registry),
    ("a converged solution gets its working",
     test_calculation_tolerance_follows_the_solver),
    ("a method block never goes quiet", test_a_method_block_never_goes_quiet),
    ("a refusal prints no number it cannot stand behind",
     test_a_refusal_prints_no_number_it_cannot_stand_behind),
    ("the factor of safety from the printed operands",
     test_calculation_reproduces_fs),
    ("the sums carry the digits to divide",
     test_calculation_sums_are_printed_precisely_enough),
    ("the equation follows the model", test_calculation_terms_follow_the_model),
    ("the per-slice terms are table columns", test_calculation_columns),
    ("the equilibrium residuals close", test_calculation_residuals),
    ("Spencer's force sums are printed and check out", test_spencer_force_sums),
    ("the mirror rule is the arithmetic behind Q", test_spencer_mirror_rule),
    ("no force the equations print is called absent",
     test_absent_features_are_really_absent),
    ("the equation numbers are in the prose",
     test_equation_numbers_are_in_the_prose),
    ("the published equation comes before this model's",
     test_the_published_equation_comes_first),
    ("the base normal is published then reduced",
     test_the_normal_force_is_published_then_reduced),
    ("the moment quotient recomposes to its page's equation",
     test_the_moment_quotient_recomposes),
    ("the evaluated equation introduces its terms",
     test_the_evaluated_equation_introduces_its_terms),
    ("every printed full form matches its page",
     test_the_full_forms_match_their_pages),
    ("each method prints its own page's equations",
     test_the_method_prints_its_own_pages_equations),
    ("a subscript is not cut short", test_scripts_are_not_cut_short),
    ("every printed symbol is defined where it is printed",
     test_printed_symbols_resolve),
    ("the prose is about the analysis", test_prose_is_about_the_analysis),
    ("the equation is cited for what it is",
     test_the_equation_is_cited_for_what_it_is),
    ("the notation matches the documentation",
     test_calculation_notation_matches_the_docs),
    ("each method's block opens by saying what it satisfies",
     test_method_summary_opens_each_block),
    ("the documentation links resolve", test_docs_links),
    ("the calculations reach the document", test_calculation_in_the_document),
    ("each calculation opens on its force diagram",
     test_force_diagram_heads_the_calculations),
    ("a force diagram is as small as it reads",
     test_force_diagram_is_as_small_as_it_reads),
    ("every figure and table is cited", test_every_block_is_cited),
    ("a citation is a live cross-reference",
     test_citations_are_cross_references),
    ("the shared-model plot", test_shared_plot),
    ("every profile line names its material",
     test_profile_lines_name_their_materials),
    ("the model figure's point-coordinate toggle",
     test_model_figure_coordinate_labels),
    ("the coordinate labels are placed clear",
     test_coordinate_labels_are_placed_clear),
    ("the dialog and its toggles", test_dialog),
    ("the dialog remembers the right things", test_dialog_settings),
    ("the finished report is opened", test_open_output),
    ("the report builds off the GUI thread", test_report_runs_off_the_gui_thread),
    ("the progress bar counts the figures", test_report_runner_progress),
    ("cancelling stops the build", test_report_runner_cancel),
    ("a build that cannot finish says so", test_report_runner_failure),
    ("a non-circular surface dims the moment methods",
     test_noncircular_dims_the_moment_methods),
    ("the slice-numbers display toggle", test_slice_numbers_display_option),
    ("the menu item and its gate", test_main_window_action),
    ("every engine's solution reaches the report",
     test_report_solutions_carry_every_engine),
]

#: Checks that need the Studio layer; skipped when PySide6 is absent.
_STUDIO_ONLY = {test_dialog, test_dialog_settings, test_open_output,
                test_report_runs_off_the_gui_thread, test_report_runner_progress,
                test_report_runner_cancel, test_report_runner_failure,
                test_noncircular_dims_the_moment_methods,
                test_slice_numbers_display_option, test_main_window_action,
                test_report_solutions_carry_every_engine}


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    checks = CHECKS
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("Report: PySide6 not installed — Studio checks skipped.")
        checks = [c for c in CHECKS if c[1] not in _STUDIO_ONLY]

    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:48s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("Analysis Report:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll report checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app_ = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
