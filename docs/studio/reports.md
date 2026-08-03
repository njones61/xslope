# Analysis Report

Studio writes the analysis up as a formal document. **File → Generate Report…**
(also on the toolbar) turns the open model and its results into a Word document
with a title page, a table of contents, running headers and footers, numbered
figures and tables, and the slice table on its own landscape page.

The report is generated from the model and the solution, not from screenshots:
every figure is rendered at 300 dpi by the same plotting code the canvas draws
with, and every table is built from the same DataFrames the solvers produce.

![The Generate Report dialog](images/reports_dialog.png){width="900"}

---

## Generating a report

The action becomes available once an analysis has been run — a report documents
results, so there is nothing to write until there are some. The dialog collects
four things.

**Output.** The format and where the file goes. Word (`.docx`) is the format
available today; PDF and LaTeX are listed and dimmed.

**Analysis.** Which solved method the report follows in detail. The
critical-surface figure and the slice table are that method's. The summary table
lists **every** method xslope offers, whichever one is picked here.

**Title page.** Project title, project number, organization and author. These
become Word document properties, so the title page, the running header and any
company template all read the same values. The organization and the author are
remembered between sessions; the project title and number belong to the project
and are not.

**Contents.** A checkbox tree of everything the report can contain. Turning a
section off removes it; turning a parent off removes its whole branch. Each box
opens on that section's own default — everything is in except the model checks,
which are opt-in. The selections are remembered, so a report composed once keeps
its shape.

Generating takes a few seconds — most of it rendering figures — and the finished
document opens in whatever the system uses for Word files. It writes one file:
the figures are embedded in the document, so nothing is left beside it.

---

## What is in the report

### Title page and contents

The title page is a left-ranged block: the organization above the project title,
a rule under it, the document type beneath, and the **Project**, **Author** and
**Date** rows below that. A field left empty prints no row — not every project
has a number, an organization or a named author. Ruled **prepared by / checked
by** signature lines are available from a checkbox and are off by default.

The table of contents lists the report's own headings, written when the report is
generated, so the contents are readable the moment the document is opened and in
any viewer. Where the sections fall is Word's to compute, so Word is asked: the
written document is opened in Word, its fields are updated and it is saved again
— a few seconds, reported in the status bar — before it is shown to you. What
you open is a finished contents page with real page numbers.

Where Word is not installed, or declines to be driven, the report is shown as
written: the heading list stands, without page numbers, and one update in Word —
right-click the contents and choose **Update Field**, or select all and press F9
— fills them in. Nothing is lost either way, and the page numbers in the footer
keep themselves current.

The footer carries *page N of M* and the report date; the header carries the
project title.

### Traceability

The reproducibility stamp: the xslope version, the input file name and its
SHA-256 digest, when the analysis was run, a mesh summary where the model carries
a mesh, and when the report was generated. The digest identifies the exact inputs
the numbers came from, which is what makes the report auditable a year later.

### Project definition

The model every analysis in the report shares, in one section:

- **A units statement** — which unit system every number in the report is in. It
  leads the section, so the numbers below it are read in known units.
- **A model figure** — geometry with material colours, and whatever else the
  model carries: water surfaces, distributed loads, reinforcement lines, piles.
  Trial failure surfaces and analysis meshes are deliberately absent; they appear
  with the analyses that use them.
- **A materials table** — only the materials the geometry references, and only
  the columns the model populates. A model with no saturated unit weights and no
  pore-pressure ratios prints neither column.
- **Water conditions** — piezometric lines, seepage head boundaries, the unit
  weight of water, whether water loads are derived or entered by hand, and the
  pore-pressure method each material uses. A model with none of these collapses
  to a single statement that the section is analysed dry.
- **Loads** — distributed loads as entered, and the seismic coefficient.
- **Reinforcement** and **Piles** — geometry and capacities, one row per member,
  in a section each. A model with no piles gets no Piles section.

The report describes what the model has and nothing else: the sentences and rows
about water, loads, reinforcement and piles are written from the model, so a
feature that is not there is not mentioned.

Where a model states its water level with a seepage head or reservoir boundary,
the water line that boundary implies is drawn on the model figure. It is the same
line the engine measures the ponded-water load against, so the load and its
source are both visible.

### Limit equilibrium analysis

- **Analysis inputs** — the method, the number of slices, the surface family and
  its defining geometry, the seismic coefficient, and the tension crack.
- **Search documentation** — how many trial surfaces the search evaluated, over
  how many refinement stages, the range of factors of safety it saw, and the
  search window it worked within. Present only when the surface was found by
  search.
- **A search results plot** — every trial surface, with the critical one
  highlighted.
- **Results** — a statement of the factor of safety, a table of **every** limit
  equilibrium method xslope offers with its solution parameters, and any
  admissibility notes the solver reported. The methods that were run report their
  own answers; the rest are solved on the critical surface the report documents,
  which takes milliseconds on a slice table that already exists. A method that
  cannot apply to the surface family, or that does not converge on it, says so in
  a row of its own rather than being left out.
- **A critical surface plot** — the failure surface with its slices and base
  stresses.
- **Rapid drawdown** — the three stage factors of safety and which one governs,
  when the run was a rapid drawdown analysis.
- **The slice table** — slice geometry, forces and strengths on a landscape page,
  with a legend beneath it defining every column. Columns that carry nothing for
  the model are left out: a section with no seismic load prints no seismic
  column.
- **Calculations** — the factor of safety worked through, in Word's own equation
  format.

### Calculations

The calculations section shows the equation the solver evaluated and the
converged numbers in it. The equation is the one published for the reported
method — the derivation on that method's documentation page, in that page's own
symbols, which the section links to — carrying only the terms the model
exercises: a section with no seismic load, no tension crack and no reinforcement
prints none of those terms, and says so.

It is compact by design. Each slice's contribution to the two sums is a column of
the slice table (`M_R` and `M_D` for the moment methods, `F_R` and `F_D` for the
force methods, plus `Q_s` and `y_Q` for Spencer), and the section references
those columns with a link to the table rather than walking through every slice.
What it prints is the equation, the sums with each Σ replaced by its value, and
the division:

$$F = \dfrac{\sum (c \Delta \ell + N' \tan \phi)\, a_S}{\sum W x_r + \sum D \cos \beta\, a_{dx}} = \dfrac{3074565}{1592951} = 1.930$$

The iterative methods show what they iterate on: Bishop and Janbu print the base
normal force N' with the factor of safety in it, so the quotient visibly closes
on itself; Janbu adds the f₀ correction line; the force-equilibrium methods state
the interslice inclination; Spencer prints its per-slice resultant Q and both
equilibrium sums evaluating to their residuals at the converged (F, θ), as
Morgenstern–Price does at its (F, λ). A rapid drawdown report works through the
governing stage and names it.

Every number in the section is printed at a precision the factor of safety can be
rebuilt from: divide the two sums as printed and the printed factor of safety
comes back, to the last digit it carries. That reproduction is a test, so the
section cannot drift away from the solver. Where a model uses a feature the
compact form cannot show — passive support in a force-equilibrium method, whose
capacity mobilizes with the soil and so carries 1/F on both sides — the section
is left out rather than printed with an equation that does not reproduce.

### Model checks

The model-check findings that were live when the analysis ran, in the checker's
own words. Reporting them is deliberate: a reviewer reads the same warnings the
engineer saw rather than having to guess whether any were raised.

The section is **off by default** — turn it on for a submittal where the checks
belong on the record. What it carries is scoped to the report: every check
declares which analyses it applies to, so a limit equilibrium report never prints
a finding about the finite element engine.

---

## Templates

The document is built on a Word template shipped with xslope. The template owns
the page size and margins, the Title, Heading, Body Text and Caption styles, the
rule under the title, and the header and footer frames; the report supplies only
content. Body text is 10.5 pt, a first-level heading 14 pt and the title 24 pt
ranged left, on one-inch margins — sized for a submittal rather than a
presentation.

Tables are fitted to the template's own text width. Each column is measured
against what it has to print — its header and every cell, in the template's
font — so a `#` column stays narrow and a `Material` column does not, and the
set is scaled to fill the width exactly. A column is never made narrower than its
longest word, except where that word is wider than an equal share of the page (a
file digest, say), which wraps instead. A bordered table is indented by one cell
margin so its border lines up with the body text rather than hanging out into the
margin; a borderless one — the label/value blocks — is not, so its text starts
where a line of prose does.

Pointing the report at a company template is planned; the metadata already maps
onto Word document properties, so a template can place them wherever it likes.

---

## From Python

The same report is available without Studio:

```python
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import solve_selected
from xslope.report import generate_report

slope_data = load_slope_data("levee.xlsx")
ok, (slice_df, surface) = generate_slices(slope_data, circle=slope_data["circles"][0])
results = solve_selected("spencer", slice_df)

bundle = {"slice_df": slice_df, "failure_surface": surface,
          "results": results, "search": None, "method": "spencer"}
ok, out = generate_report(
    slope_data, {"lem": bundle},
    {"title": "North Levee", "author": "A. Engineer",
     "organization": "Example Engineering", "input_path": "levee.xlsx"},
    "north_levee_report.docx")
```

`generate_report` returns the package's usual `(success, result)` pair; on
success the result carries the document path, the content tree, and the caption
of every figure the document embeds. The PNGs are drawn in a temporary directory
that goes away with the call — pass `figure_dir` to keep them. Pass several
bundles as a list to report more than one method, and use the `method` option to
choose which one the detail follows.

The options dictionary is the dialog's checkbox tree — one key per box, all
documented in `xslope.report.DEFAULT_OPTIONS`. Report generation is headless: it
renders through the Agg backend, opens no windows, and never opens the finished
document.

### The slice-table columns

The slice table's columns are declared in `xslope.columns`, which states each
column's label, definition, physical quantity and format, and whether it belongs
in a report. The table headers, the legend printed beneath the table, and the
choice of which columns to print all come from that one declaration. A few
columns are marked as the report's own — the per-slice terms of the factor of
safety equation, which the calculations section computes and adds to the table it
prints.
