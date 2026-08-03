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
     and the slice table, and NOT the factor-of-safety summary: every solved
     method's answer is reported whichever one the detail follows.

  D. THE DOCUMENT — the .docx is a real OOXML package: the title is in it, the
     template's styles are referenced, the slice table sits in a landscape
     section, the table of contents is a TOC field, and the running head and
     foot carry live fields.

  E. THE COLUMN REGISTRY — every column the registry marks report-worthy exists
     in a solved slice_df, so the registry cannot drift away from the solver.

  F. THE SHARED PLOT — ``mode="shared"`` suppresses the trial surfaces and draws
     the water line a reservoir boundary states, which is where the derived
     water loads on the same figure come from.

  G. THE DIALOG — it constructs offscreen, its toggles move the options it hands
     over, and what it remembers survives a second construction.

One small LEM model is solved once and shared. A second model (a dam with a
reservoir head boundary) is loaded, not solved: the shared-plot check reads its
inputs only.
"""

import os
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

_SOLVED = {}


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
        opts.update({"pd_figure": False, "lem_search_figure": False,
                     "lem_solution_figure": False})
    opts.update(options or {})
    if figure_dir is None:
        figure_dir = tempfile.mkdtemp(prefix="xslope_report_")
    return build_report(slope_data, solutions, opts, figure_dir)


# --------------------------------------------------------------------------
# A. the content tree
# --------------------------------------------------------------------------

#: The sections a fully-selected report of a solved LEM model has, in order.
EXPECTED_SECTIONS = [
    (1, "Traceability"),
    (1, "Project Definition"),
    (2, "Materials"),
    (2, "Water Conditions"),
    (2, "Loads"),
    (2, "Reinforcement and Piles"),
    (2, "Units"),
    (1, "Limit Equilibrium Analysis"),
    (2, "Analysis Inputs"),
    (2, "Search for the Critical Surface"),
    (2, "Results"),
    (2, "Slice Table"),
    (1, "Model Checks"),
]


def test_tree():
    """The tree has the specified sections, in order."""
    fails = []
    report = _build()
    got = report.section_titles()
    if got != EXPECTED_SECTIONS:
        fails.append(f"section order is {got}, expected {EXPECTED_SECTIONS}")
    if not report.meta.get("title"):
        fails.append("the report carries no title")
    if not report.meta.get("date"):
        fails.append("the report carries no date")
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

    fs = tables.get("Computed factors of safety")
    if fs is None:
        fails.append("there is no factor-of-safety table")
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
        ({"model_checks": False}, "Model Checks"),
        ({"traceability": False}, "Traceability"),
        ({"lem_search": False}, "Search for the Critical Surface"),
        ({"lem_slice_table": False}, "Slice Table"),
        ({"pd_materials": False}, "Materials"),
        ({"pd_water": False}, "Water Conditions"),
        ({"pd_units": False}, "Units"),
    ]
    for opts, gone in cases:
        titles = [t for _lvl, t in _build(opts).section_titles()]
        if gone in titles:
            fails.append(f"{opts} left {gone!r} in the report")
        if len(titles) != len(full) - 1:
            fails.append(f"{opts} changed the section count by "
                         f"{len(full) - len(titles)}, not 1")

    # A parent off takes the whole branch.
    titles = [t for _lvl, t in _build({"project_definition": False}).section_titles()]
    for gone in ("Project Definition", "Materials", "Water Conditions", "Units"):
        if gone in titles:
            fails.append(f"project_definition=False left {gone!r} in the report")
    titles = [t for _lvl, t in _build({"lem": False}).section_titles()]
    for gone in ("Limit Equilibrium Analysis", "Slice Table", "Results"):
        if gone in titles:
            fails.append(f"lem=False left {gone!r} in the report")

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
        if f"{name} critical surface" not in sources:
            fails.append(f"method={name}: no critical-surface figure for it "
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

    # The summary does NOT follow the pick.
    def fs_table(report):
        for t in report.tables():
            if t.caption == "Computed factors of safety":
                return (t.headers, t.rows)
        return None

    if fs_table(a) != fs_table(b):
        fails.append("the factor-of-safety summary changed with the picked "
                     "method; it must report every solved method either way")
    if fs_table(a) is None:
        fails.append("there is no factor-of-safety summary")

    # An unknown pick falls back rather than failing.
    report = _build({"method": "no_such_method"})
    if "Limit Equilibrium Analysis" not in [t for _l, t in report.section_titles()]:
        fails.append("an unrecognised method dropped the LEM section")
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
        if 'w:dirty="true"' not in doc:
            fails.append("no field is marked for refresh; the TOC would never build")
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
    for column in cols.report_columns():
        if column.key not in have:
            fails.append(f"report column {column.key!r} is not in slice_df")
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
        fails.append(f"the weight column is not unit-labelled: {headers}")
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
# F. the shared plot
# --------------------------------------------------------------------------

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

    # The mesh belongs to the mesh figures, not to the shared model.
    if dam.get("mesh") is not None:
        meshed = draw(dam, "shared")
        from matplotlib.collections import LineCollection
        if any(isinstance(c, LineCollection) for c in meshed.collections):
            fails.append("the shared plot drew the analysis mesh")
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
        if dlg.method.currentData() != "bishop":
            fails.append(f"the method picker opens on "
                         f"{dlg.method.currentData()!r}, not the results view's")
        offered = [dlg.method.itemData(i) for i in range(dlg.method.count())]
        if offered != ["spencer", "bishop"]:
            fails.append(f"the method picker offers {offered}")

        # Only DOCX can be picked in S1; the rest are listed and dimmed.
        for i in range(dlg.format.count()):
            key = dlg.format.itemData(i)
            enabled = dlg.format.model().item(i).isEnabled()
            if (key == "docx") != enabled:
                fails.append(f"format {key!r} is {'enabled' if enabled else 'dimmed'}")

        opts = dlg.options()
        for key in ("traceability", "project_definition", "lem", "model_checks",
                    "lem_slice_table", "pd_materials"):
            if opts.get(key) is not True:
                fails.append(f"{key} is not on by default")
        if opts.get("signature_lines") is not False:
            fails.append("signature lines are on by default")

        # A toggle reaches the options, and a parent takes its children.
        dlg._items["lem_slice_table"].setCheckState(0, Qt.Unchecked)
        if dlg.options().get("lem_slice_table") is not False:
            fails.append("unchecking the slice table did not reach the options")
        dlg._items["project_definition"].setCheckState(0, Qt.Unchecked)
        after = dlg.options()
        if after.get("pd_materials") is not False:
            fails.append("a section turned off left its children on")
        if not dlg._items["pd_materials"].isDisabled():
            fails.append("a section turned off left its children live")

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
        if "Materials" not in titles:
            fails.append("the dialog's options lost the materials section")
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
    """A document is opened; a LaTeX source reveals its folder instead."""
    fails = []
    _app()
    from studio import report_dialog

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
            tex = os.path.join(tmp, "r.tex")
            if report_dialog.open_output(tex, "latex") != "folder":
                fails.append("a .tex opened as a document; it should reveal the "
                             "folder")
            if opened[-1] != tmp:
                fails.append(f"the folder revealed was {opened[-1]!r}")
    finally:
        report_dialog.QDesktopServices.openUrl = real
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


CHECKS = [
    ("the content tree and its section order", test_tree),
    ("the tables carry the model's numbers", test_tables_carry_the_model),
    ("every figure is rendered to a file", test_figures_rendered),
    ("the traceability stamp", test_traceability),
    ("the content toggles remove what they name", test_toggles),
    ("the method picker drives the detail only", test_method_picker),
    ("the .docx and its structure", test_docx),
    ("the shipped template is reproducible", test_docx_template),
    ("the slice-column registry", test_column_registry),
    ("the shared-model plot", test_shared_plot),
    ("the dialog and its toggles", test_dialog),
    ("the dialog remembers the right things", test_dialog_settings),
    ("the finished report is opened", test_open_output),
    ("the menu item and its gate", test_main_window_action),
]

#: Checks that need the Studio layer; skipped when PySide6 is absent.
_STUDIO_ONLY = {test_dialog, test_dialog_settings, test_open_output,
                test_main_window_action}


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
