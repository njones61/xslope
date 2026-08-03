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
     section, the table of contents is a TOC field, and the running head and
     foot carry live fields. Its tables are fitted to the page they sit on —
     fixed columns, measured, summing to the text width, indented so their
     borders line up with the body text.

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

#: The sections a report of this solved LEM model has, in order, with the shipped
#: defaults. The sample carries reinforcement lines and no piles, so Reinforcement
#: is here and Piles is not — the two are separate sections, each present only
#: where the model has that feature. Model Checks is opt-in and is absent.
EXPECTED_SECTIONS = [
    (1, "Traceability"),
    (1, "Project Definition"),
    (2, "Materials"),
    (2, "Water Conditions"),
    (2, "Loads"),
    (2, "Reinforcement"),
    (1, "Limit Equilibrium Analysis"),
    (2, "Analysis Inputs"),
    (2, "Search for the Critical Surface"),
    (2, "Results"),
    (2, "Slice Table"),
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
        ({"pd_materials": False}, "Materials"),
        ({"pd_water": False}, "Water Conditions"),
        ({"pd_reinforcement": False}, "Reinforcement"),
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
    for gone in ("Project Definition", "Materials", "Water Conditions",
                 "Reinforcement"):
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
    if _fs_rows(a) != _fs_rows(b):
        fails.append("the factor-of-safety summary changed with the picked "
                     "method; it must report every method either way")
    if _fs_rows(a) is None:
        fails.append("there is no factor-of-safety summary")

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


def test_fs_table_lists_every_method():
    """The summary reports EVERY method the solver offers, not only the ones the
    caller happened to run.

    The sample was solved by two methods. A summary listing two rows is the bug
    Norm reported: the reader is left asking what the other five would have said.
    The unsolved ones are solved here, on the report's own critical surface — so
    the check is that every supported method has a row, that the two that WERE run
    carry their own answers unchanged, and that the filled-in rows are real
    numbers rather than blanks.
    """
    fails = []
    from xslope.report import method_label, supported_methods

    _slope_data, solutions = _solved()
    report = _build()
    got = _fs_rows(report)
    if got is None:
        return ["there is no factor-of-safety summary"]
    _headers, rows = got
    by_method = {r[0]: r for r in rows}

    supported = supported_methods()
    if len(supported) < 7:
        fails.append(f"the solver offers only {supported}; the report enumerates "
                     f"the solver, so this list is the solver's own")
    for name in supported:
        label = method_label(name)
        if label not in by_method:
            fails.append(f"{label} has no row in the summary")
    if len(rows) != len(supported):
        fails.append(f"{len(rows)} rows for {len(supported)} supported methods")

    # The methods that were run report their own numbers, not a re-solve.
    for b in solutions["lem"]:
        row = by_method.get(method_label(b["method"]))
        if row is None:
            continue
        if row[1] != f"{b['results']['FS']:.3f}":
            fails.append(f"{b['method']}: the summary says {row[1]!r}, the run "
                         f"gave {b['results']['FS']:.3f}")

    # And the filled-in ones are answers, not blanks.
    solved_here = [r for r in rows
                   if r[0] not in {method_label(b["method"])
                                   for b in solutions["lem"]}]
    if not solved_here:
        fails.append("nothing was solved at report-build time; the check proves "
                     "nothing")
    for row in solved_here:
        if not row[1].strip():
            fails.append(f"{row[0]} has an empty factor of safety cell; a method "
                         f"that did not converge must say so")

    # A method that does not converge is stated, not dropped. Forced by handing
    # the builder a slice table the solvers cannot use: every unsolved method
    # fails on it, and every one of them still gets a row.
    from xslope.report import build_report
    slope_data, _sol = _solved()
    bundle = dict(solutions["lem"][0])
    df = bundle["slice_df"].copy()
    df["w"] = float("nan")             # no method can produce an answer from this
    bad = {"lem": [{"slice_df": df, "failure_surface": bundle["failure_surface"],
                    "results": {"FS": 1.5}, "search": None, "method": "spencer"}]}
    with tempfile.TemporaryDirectory() as tmp:
        broken = build_report(slope_data, bad,
                              {"pd_figure": False, "lem_solution_figure": False,
                               "lem_search_figure": False}, tmp)
    got = _fs_rows(broken)
    if got is None:
        fails.append("a report whose extra methods all failed lost its summary "
                     "table entirely")
    else:
        _h, rows = got
        stated = [r for r in rows if "did not converge" in r[1]]
        if not stated:
            fails.append(f"no method reported 'did not converge' on an "
                         f"unsolvable slice table: {rows}")
        if len(rows) != len(supported):
            fails.append(f"the unsolvable run dropped rows: {len(rows)} of "
                         f"{len(supported)}")
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


def _style_cell_margin(styles_xml, style_id):
    """The leading cell margin a table style declares, in twips."""
    import re
    style = re.search(rf'<w:style [^>]*w:styleId="{style_id}".*?</w:style>',
                      styles_xml, re.S)
    if style is None:
        return None
    left = re.search(r"<w:tblCellMar>.*?<w:left [^>]*w:w=\"(\d+)\"",
                     style.group(0), re.S)
    return int(left.group(1)) if left else None


def test_table_geometry():
    """Every table is laid out to the page it sits on.

    Two defects this defends against, both of which show in Word and neither of
    which Word fixes on its own: a table with no indent hangs its left border out
    into the margin, and an autofitting table gives "#" the same width as
    "Material". The cure is a fixed layout with measured columns that sum to the
    text width, and an indent of exactly one cell margin.
    """
    import re
    fails = []
    from xslope.report import generate_report

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

        margin = _style_cell_margin(styles, style.group(1) if style else
                                    "TableNormal")
        indent = re.search(r"<w:tblInd [^>]*/>", tbl_pr)
        if indent is None:
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
        if abs(sum(grid) - usable) > 1:
            fails.append(f"{where} spans {sum(grid)} twips of a {usable}-twip "
                         f"text width")
        if sum(grid) > usable:
            fails.append(f"{where} is wider than the page's text width")
        if min(grid) <= 0:
            fails.append(f"{where} has a column of {min(grid)} twips")

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
        from xslope.report_docx import DEFAULT_CELL_MARGIN, TABLE_PT, _text_width
        if grid[1] != max(grid):
            fails.append(f"the long Finding column is not the widest: {grid}")
        # It is wide by wrapping, not by starving: "Warning" beside it still
        # prints on one line.
        needed = _text_width("Warning", "Calibri", TABLE_PT) + 2 * DEFAULT_CELL_MARGIN
        if grid[0] < needed:
            fails.append(f"the Severity column is {grid[0]} twips, under the "
                         f"{needed:.0f} 'Warning' needs — the Finding column "
                         f"starved it")
    elif not findings:
        fails.append("the model-check findings table was not written")
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
    water = next((s for _l, s in _sections(report) if s.title == "Water Conditions"),
                 None)
    if water is None:
        return ["there is no Water Conditions section"]
    if [b.kind for b in water.blocks] != ["prose"]:
        fails.append(f"a dry model's Water Conditions holds "
                     f"{[b.kind for b in water.blocks]}; it collapses to one "
                     f"statement")
    else:
        text = water.blocks[0].text
        for want in ("no groundwater", "dry"):
            if want not in text:
                fails.append(f"the dry-model statement does not say {want!r}: "
                             f"{text!r}")
        if "zero" not in text:
            fails.append("the dry-model statement does not say pore pressures are "
                         "zero")

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
                                      {"pd_figure": False, "lem": False}, tmp)
        sub = next((s for _l, s in _sections(dam_report)
                    if s.title == "Water Conditions"), None)
        if sub is None:
            fails.append("the dam report has no Water Conditions section")
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
    return fails


def _sections(report):
    """Every section in the report, as ``(level, section)``."""
    out = []
    for s in report.sections:
        out.extend(s.walk())
    return out


def test_reinforcement_and_piles_split():
    """Reinforcement and piles are separate sections, each present only where the
    model has that feature."""
    fails = []
    from xslope.report import build_report

    slope_data, solutions = _solved()
    titles = [t for _l, t in _build().section_titles()]
    if not slope_data.get("reinforcement_lines"):
        return ["the sample carries no reinforcement; the split is untested"]
    if slope_data.get("pile_lines"):
        return ["the sample now carries piles; the 'absent' half is untested"]
    if "Reinforcement" not in titles:
        fails.append("the reinforcement section is missing")
    if "Piles" in titles:
        fails.append("a model with no piles was given a Piles section")
    if "Reinforcement and Piles" in titles:
        fails.append("the combined section is still built")

    # A model with piles and no reinforcement gets the other half, and only it.
    piled = dict(slope_data)
    piled["reinforcement_lines"] = []
    piled["pile_lines"] = [{"label": "P1", "x1": 20.0, "y1": 40.0,
                            "x2": 20.0, "y2": 20.0, "H": 5000.0,
                            "theta_p": 0.0, "D_pile": 2.0, "S": 8.0,
                            "appl": "active"}]
    with tempfile.TemporaryDirectory() as tmp:
        report = build_report(piled, solutions,
                              {"pd_figure": False, "lem": False}, tmp)
    titles = [t for _l, t in report.section_titles()]
    if "Piles" not in titles:
        fails.append("a model with piles got no Piles section")
    if "Reinforcement" in titles:
        fails.append("a model with no reinforcement got a Reinforcement section")
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
        for key in ("traceability", "project_definition", "lem",
                    "lem_slice_table", "pd_materials"):
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
    ("every method is in the summary table", test_fs_table_lists_every_method),
    ("the water prose follows the model", test_water_prose_is_conditional),
    ("reinforcement and piles are separate", test_reinforcement_and_piles_split),
    ("the model checks are opt-in and scoped",
     test_model_checks_default_and_filtering),
    ("an empty title-page field prints no row", test_title_page_omits_empty_rows),
    ("the .docx and its structure", test_docx),
    ("the tables are fitted to the page", test_table_geometry),
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
