---
title: "Tutorial W-3 — Generating the Analysis Report"
description: "A solved model turned into a Word calculation package in XSLOPE Studio — the Generate Report dialog walked field by field, the document it writes, and the same report built a second time on a company letterhead."
---

# Tutorial W-3 — Generating the Analysis Report

An analysis that stays on the canvas is not yet a deliverable. **File → Generate
Report…** writes one: a Word document with a title page, a table of contents,
running heads and footers, numbered figures and tables, and a section for every
engine that was run — down to the equation each stability method evaluated and
the converged numbers in it. Nothing in it is a screenshot. The plots are drawn
at 300 dpi by the same code the canvas draws with, and every table comes from the
same DataFrame the solver produced.

We generate two reports on one model. The first documents Bishop and Spencer on
the page layout XSLOPE ships. The second goes out on a firm's letterhead, built
on a Word template this page ships.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Seepage + limit equilibrium + finite element</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~20 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Generate a report from a solved model, read what each of its
sections states, and put the same report on a company template.
</div>
<p><span class="tg-pill">generate report</span><span class="tg-pill">attached solutions</span><span class="tg-pill">method selection</span><span class="tg-pill">title page</span><span class="tg-pill">contents tree</span><span class="tg-pill">flow net</span><span class="tg-pill">calculations</span><span class="tg-pill">slice table</span><span class="tg-pill">company template</span><span class="tg-pill">Word styles</span><span class="tg-pill">letterhead</span></p>
<div class="tgm-model" markdown>
**Model** — [xslope_johnson_res_solved.xlsx](files/xslope_johnson_res_solved.xlsx),
COMBO-1's dam with its seepage and finite element solutions saved beside it

**Reports** — [w03_report.docx](files/w03_report.docx) and
[w03_report_example_template.docx](files/w03_report_example_template.docx)

**Template** — [report_template_example.docx](files/report_template_example.docx)
</div>
</div>

Everything below happens in Studio. The same document can be written from a
script, which [Analysis Report](../studio/reports.md#from-python) covers.

---

## The model

[xslope_johnson_res_solved.xlsx](files/xslope_johnson_res_solved.xlsx) holds the
Johnson Reservoir dam of [COMBO-1](combo01_seepage_stability.md) — an 80 ft
zoned embankment with a shell, a keyed clay core and a foundation, under 60 ft of
water. Eight companion files sit beside it: `..._mesh.json` carries the mesh all
three analyses share, `..._seep.csv` the steady head and pore pressure at every
node, and six `..._fem_*` files the strength reduction solution and the mechanism
it failed on. Download the workbook and its companions into one folder, then open
the workbook with **File → Open…**.

Studio reads those companions on open, so the seepage and finite element
solutions are attached before anything is run. A report **documents** an
attached seepage or finite element solution; it never re-solves either, which
is why this file ships solved — the strength reduction run behind it is the
kind of solve no report should spend again. Limit equilibrium is the one
engine a report runs itself: the slice computations are fast, so any method
ticked without a result is solved on the fly — which is what lets several
methods go into one report.

No companion file carries a limit equilibrium solution, and none is needed —
the report will run whatever methods we tick. Still, a quick run confirms the
model arrived intact. Switch the mode strip to **LEM** (`Ctrl+1`) and click
**Run → Run LEM…**. Leave **Method** on **Spencer**, **Analysis** on
**Auto search** and **Number of slices** at 40, and click **Run**. The search
returns **FS = 1.248**, COMBO-1's own answer on this model.

---

## Generating the report

With three engines' results in the session, the report has something to document.
Click **File → Generate Report…**.

![The Generate Report dialog on the solved model, with Bishop and Spencer ticked and the title page filled in](images/w03_report_dialog.png){width=1000}

Three groups sit on the left — Output, Analysis and Title page — and a contents
tree on the right. Every field opens on a default that works; this report changes
the methods and the title page, and nothing else.

**Output.** **Format** offers Word (`.docx`); PDF is listed and dimmed. **Save
as** defaults to `<model>_report.docx` beside the model. **Template** reads
*Shipped template*, which gives the layout we want here; the company template
below changes it. Leave all three.

**Analysis.** **Methods** lists every method XSLOPE offers, ticked on the one the
results view is showing — Spencer, which we just ran. Tick **Bishop's Simplified
Method** as well. A method that has not been run is run by the report, exactly as
ticking it in the Run dialog would have run it: its own search for its own
critical surface. So Bishop costs a search here, and the two settle on different
surfaces.

**Title page.** These four fields become Word document properties, which the
title page, the running head and a company template all read. Fill in **Project
title** `Johnson Reservoir Dam — Steady Seepage and Stability`, **Project number**
`2026-231`, **Organization** `Example Engineering` and **Author** with your own
name, and leave **Signature lines (prepared by / checked by)** unticked — a
submittal needing ruled signature lines gets them from that box.

**Contents.** One checkbox per section, opening on that section's own default. A
parent turned off dims its whole branch. Leave the tree alone; its defaults
document every analysis this model carries.

The dialog is composed. Click **Generate**.

The progress bar counts the figures by name, and the window stays live —
generating runs off the GUI thread, and **Cancel** stops the build at the next
figure. A labeled stretch at the end builds the contents page's page numbers:
Word does it where the machine has Word, LibreOffice where it does not, and where
neither is available the report is complete but its contents page lists the
headings unnumbered.

---

## What the document carries

The report runs to 28 pages and 24 figures:

(The linked copies below carry their figures at reduced resolution to keep the
download small; a report you generate yourself carries them at full
resolution.)
[w03_report.docx](files/w03_report.docx). Its sections follow the analyses, with
seepage ahead of the one stability section that consumes its answer. Both methods
sit inside that section, §4.5 Bishop's Simplified Method and §4.6 Spencer's
Method, under a §4.4 *Factors of Safety* table that lists the two side by side and
records that each searched for a surface of its own.

| Analysis | Where the answer came from | Result |
| --- | --- | :---: |
| Seepage | attached to the file | 1.925 ft³/day per ft |
| Bishop's Simplified Method | searched by the report | 1.222 |
| Spencer's method | the search we ran | 1.248 |
| Strength reduction | attached to the file | 1.230 |

<!-- test: file=files/xslope_johnson_res_solved.xlsx, type=circular_search, method=bishop, num_slices=40, expected_fs=1.222, tolerance=0.005 -->

The seepage section states the conductivities and the mesh it was solved on, then
the four fields the solution carries: the flow net, pore pressure, velocity
magnitude with the vectors over it, and hydraulic gradient.

![The seepage results page: the computed flow through the section, the flow net, the pore pressure field and the velocity magnitude](images/w03_page_seep.png){width=620}

Each method then gets a block of its own: its search, its critical surface with
the slices and base stresses on it, a slice table on a landscape page with column
totals and a legend, and the calculations — the published equation for that
method, in the symbols of its own documentation page, with the converged numbers
in it.

![Spencer's calculations, §4.6.4: the equations the solver evaluated, in the symbols of Spencer's own documentation page, over the nomenclature that defines them](images/w03_page_spencer.png){width=620}

Spencer's method solves for a pair — a factor of safety and an interslice
inclination θ — at which both the force sum and the moment sum vanish, so it
prints both sums rather than a single quotient. On the page after this one they
close at F = 1.248 and θ = −16.27°, with residuals of −2.728e-11 and 6.752e-09
against terms coming to 369995 and 83989264. Every number is printed at a
precision the factor of safety can be rebuilt from: Bishop divides its own two
sums, 475756 and 389345, and prints the 1.222 they give.

### Check its work

- **Open the document.** Its contents page lists the sections with a page number
  against each; headings without them mean the finishing step did not run, and
  Word rebuilds them with F9.
- **Section 4 carries both methods** — §4.5 Bishop at 1.222 and §4.6 Spencer at
  1.248, with §4.4 *Factors of Safety* comparing them. A method missing means its
  box was never ticked.
- **The seepage and finite element sections state the attached solutions'
  numbers**, 1.925 ft³/day per ft and 1.230. Numbers that differ mean something
  was solved again.
- **The traceability stamp names the file the numbers came from**, with its
  SHA-256 digest. A digest that does not match the workbook means the report
  documents a different file.

---

## A company template

A report goes out on the firm's letterhead by being built on the firm's own Word
template. Start from the one XSLOPE ships —
[report_template.docx](../studio/files/report_template.docx) — and edit it in
Word, where its page size and margins, its header and footer and the fonts and
colors of its styles are all yours. Three edits make a letterhead:

- **A logo in the header**, sized to the header frame.
- **The firm's name**, beside the logo.
- **A footer line** naming the firm and the document type.

Where those three go decides whether they survive. The report writes its running
head into the first paragraph of the header and *page N of M* into the first
paragraph of the footer, and deletes every other paragraph in both — so a logo
dropped into a header paragraph, or a firm name typed on a second line under it,
is gone from the finished report. What the template needs to keep has to sit
outside any paragraph, and a one-row table above each does it. Content placed that
way prints on every page after the title page, which carries no head or foot of
its own.

Two rules decide whether a report can be built on the result at all. **Keep the
style names** — the report is written in *Title*, *Heading 1*, *Heading 2*,
*Heading 3*, *Body Text* and *Caption*, and asks the template for them by name.
Restyle them freely; renaming one is refused when Generate is pressed, with the
missing style named. And **leave the metadata to the fields** — the four
title-page entries arrive as Word document properties a `DOCPROPERTY` field can
place anywhere, under Word's own property names rather than the dialog's labels:

| Dialog field | Word document property |
| --- | --- |
| Project title | `Title` |
| Project number | `Subject` |
| Organization | `Category` |
| Author | `Author` |

XSLOPE ships that template for a firm invented for this page —
[report_template_example.docx](files/report_template_example.docx), the shipped
file with a logo and the firm's name in a header table and a footer line above the
page numbers, nothing else changed. Download it and open the dialog again with
**File → Generate Report…**.

![The dialog's Output group with the example template chosen](images/w03_report_template.png){width=644}

Under **Output**, click **Browse…** on the **Template** row and choose
`report_template_example.docx`; the field replaces *Shipped template* with the
path. Change **Save as** to `w03_report_example_template.docx` so the first report
survives. Under **Analysis**, untick **Bishop's Simplified Method**, leaving
Spencer's finished search as the only stability run, and under **Title page** set
**Organization** to `ACME Geotechnical`. Click **Generate**.

![The first body page of the report built on the example template: the logo and firm name in the head, the firm's line in the foot](images/w03_page_template.png){width=620}

The letterhead prints on every page after the title page, with the report's own
running head under it. **Template** is remembered between sessions, so the next
report on any model opens on the same letterhead. The 22-page result downloads as
[w03_report_example_template.docx](files/w03_report_example_template.docx).

### Check its work

- **The logo and firm name are in the head of page 2 and every page after it**,
  including the landscape slice table. A letterhead that reaches the contents page
  and stops was put in a paragraph rather than a table.
- **The footer carries the firm's line above *page N of M*.** Both print, because
  the report's own paragraph was left where it was.
- **The title page reads ACME Geotechnical**, from the Organization field rather
  than from the template — which is what lets one template serve every project.

---

## Conclusion

This tutorial covered:

- A report generated from a model whose seepage and finite element solutions were
  attached rather than re-solved, documenting four results from one run.
- The Generate Report dialog's three groups and its contents tree, and what
  ticking an unrun method costs: Bishop searched for its own surface and reported
  1.222 against Spencer's 1.248, the pair §4.4 sets side by side.
- What the document states — the flow net and the fields beside it, and each
  method's equation with its converged numbers in it.
- A company template made from the shipped one: keep the style names, place the
  metadata with the document properties the dialog fills, and put the letterhead
  outside the paragraphs the report deletes.

**Where to go next:** [Analysis Report](../studio/reports.md) documents every
section, every content option and the script route;
[COMBO-1](combo01_seepage_stability.md) builds and runs the three analyses this
page reports on.
