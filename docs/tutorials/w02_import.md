---
title: "Tutorial W-2 — Bringing in Models from CAD and Other Programs"
description: "Three importers walked end to end in XSLOPE Studio — a layered DXF section mapped layer by layer through the import wizard, a GeoStudio SLOPE/W verification model whose analysis and solved circle come across whole, and a Rocscience Slide2 scenario that arrives without a failure surface — with the notes each import reports and the factor of safety each imported model returns."
---

# Tutorial W-2 — Bringing in Models from CAD and Other Programs

A section that already exists somewhere else does not have to be redrawn. XSLOPE
Studio reads three outside formats: a **DXF** drawing from any CAD program, a
**GeoStudio SLOPE/W** project (`.gsz`), and a **Rocscience Slide2** model
(`.sli`, `.slim` or `.slmd`). The three differ in how much they can say. A DXF
carries lines on layers and nothing about what those lines mean, so its import
asks; the two program files already know their own regions, materials and water
conditions, so theirs ask only which analysis or scenario to take. All three end
the same way — a new unsaved project on the canvas, and a list of anything that
did not come across.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~25 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Import a layered CAD drawing and map each layer to an input
feature, import a GeoStudio and a Slide2 model, read the notes each import reports
against what the source model defines, and check every imported model before
trusting the factor of safety it returns.
</div>
<p><span class="tg-pill">DXF import</span><span class="tg-pill">layer mapping</span><span class="tg-pill">material zones</span><span class="tg-pill">piezometric line</span><span class="tg-pill">distributed load</span><span class="tg-pill">failure circles</span><span class="tg-pill">GeoStudio import</span><span class="tg-pill">analysis picker</span><span class="tg-pill">imported trial circle</span><span class="tg-pill">Slide2 import</span><span class="tg-pill">scenario picker</span><span class="tg-pill">import notes</span><span class="tg-pill">Spencer</span><span class="tg-pill">circular search</span></p>
<div class="tgm-model" markdown>
**Part 1 drawing** — [w02_section.dxf](files/w02_section.dxf), the layered CAD
section we import below

**Part 1 completed model** — [w02_section_imported.xlsx](files/w02_section_imported.xlsx),
that drawing after the import and after the properties a drawing cannot carry are
typed in

**Parts 2 and 3** run on the vendors' own published models, which are theirs to
distribute — each part says where to download the file.
</div>
</div>

Everything below happens in Studio. DXF and GeoStudio also have a script route,
written up under [DXF Import/Export](../usage/dxf.md) and
[GeoStudio Import/Export](../usage/geostudio.md); all three importers are covered
as reference at
[Studio → DXF import and export](../studio/analysis.md#dxf-import-and-export).

---

## Part 1 — A CAD drawing

[w02_section.dxf](files/w02_section.dxf) holds a highway embankment: 30 ft of
compacted fill on 2:1 side slopes, built over 15 ft of silty clay on dense sand,
with the ground water table a few feet down and a 600 psf surcharge across a 30 ft
strip of the crest. Dimensions are in feet.

A DXF says where the lines are and which layer each one sits on. It says nothing
about unit weights, strengths, or the pressure under a load — CAD has no place to
put them. So the import brings the geometry across and leaves placeholders for the
rest.

Six layers carry the section:

| Layer | What it holds |
| --- | --- |
| `EMBANKMENT` | the fill zone, as a closed outline |
| `SILTY_CLAY` | the upper foundation zone |
| `DENSE_SAND` | the lower foundation zone |
| `PIEZO` | the water table, as an open polyline |
| `DLOADS` | the strip the surcharge acts over |
| `SEARCH_CIRCLES` | two starting circles, each an arc plus its center point |

Those names are XSLOPE's own reserved layer names, which the wizard recognizes.
A drawing from outside CAD will use whatever names its author chose, and the
wizard handles that too — every row can be overridden.

### Importing it

With the drawing downloaded, we start the import from the File menu. Click
**File → Import DXF…** and choose `w02_section.dxf`. The wizard lists every layer
in the file:

![The DXF import wizard on w02_section.dxf: six layers, each with its contents, its suggested target and its material name](images/w02_dxf_wizard.png){width=690}

The table carries four columns. **Layer** and **Contents** are read out of the
file — Contents counts what the layer holds, which is how a layer of closed zones
is told from one open polyline. **Import as** carries the choice, and **Material**
applies to the two targets that need one.

The **Import as** column offers seven targets:

| Import as | What the layer becomes |
| --- | --- |
| Ignore | nothing — the layer is skipped |
| Material zone | one closed material region per ring |
| Profile line | the top of a material layer |
| Piezometric line | the piezometric surface (a second one becomes Piezo Line 2) |
| Distributed load | a load block, at zero pressure |
| Reinforcement | a reinforcement line, at zero capacity |
| Failure circles | starting circles for the search |

Material zone and Profile line both use the **Material** column, and two layers
given the same material name merge into one material. The other targets leave that
cell grayed out.

Every row here is already right, because the drawing uses the reserved names. We
change one thing anyway, and it is cosmetic: the **Material** column defaults to
the layer name, and DXF layer names are uppercase. Type `embankment`,
`silty clay` and `dense sand` into the three material rows, then click **OK**.

### What arrives

The import replaces whatever project was open with a new, unsaved one:

![Studio after the DXF import, with the properties filled in: three material zones, the piezometric line, the surcharge and two starting circles](images/w02_dxf_window.png){width=1000}

The Inputs tree counts what came across — 3 materials, 3 polygons, 2 circles, the
piezometric line's 3 points and 1 distributed load. Studio also reports anything
it could not carry cleanly, in a box and in the Log pane. This drawing produces
one note: load magnitudes and reinforcement strengths import as zero.

Two things still have to be typed, and neither was ever in the file. We open the
materials editor and fill in the three rows:

| name | γ (pcf) | option | c (psf) | φ (deg) | u |
| --- | :---: | --- | :---: | :---: | --- |
| `embankment` | 125 | `mc` | 250 | 28 | `piezo` |
| `silty clay` | 118 | `mc` | 700 | 0 | `piezo` |
| `dense sand` | 132 | `mc` | 0 | 36 | `piezo` |

Then we open the distributed loads editor and set the surcharge to 600 psf at both
of its points. The circles arrived with their centers and radii, and with Depth 0;
the search moves both on its own, so they need nothing further.

The result ships as
[w02_section_imported.xlsx](files/w02_section_imported.xlsx), which opens directly
for anyone skipping the typing.

### Running it

With the properties in, the model is complete. Click **Run LEM**, choose
**Spencer** and **Auto search**, leave the slice count at 40, and click **OK**.

The search returns **FS = 1.185**, on a circle bottoming out at elevation −15 —
the contact between the silty clay and the dense sand. The critical surface cuts
through the foundation clay rather than staying in the fill, which is what the
strengths point to: 700 psf with no friction, under 30 ft of embankment.

<!-- test: file=files/w02_section_imported.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.185, tolerance=0.005, benchmark=W-2-dxf -->

### Check its work

- **The geometry closed.** Three zones, stacked with no gap and no overlap. A DXF
  ring that failed to close would arrive as a missing zone, and the Inputs tree
  counts them.
- **The water table sits below the ground surface everywhere.** It runs from −4 ft
  at the left edge to +2 ft under the crest — two feet into the base of the fill,
  and 28 ft below the crest.
- **The surcharge reads 600 psf, not 0.** The import warns about this, and a load
  left at zero changes nothing and looks like it worked.
- **The critical surface is deep-seated.** It runs tangent to the top of the dense
  sand, daylights about 8 ft beyond the toe on the natural ground rather than on
  the face, and exits on the crest under the surcharge strip.

---

## Part 2 — A GeoStudio SLOPE/W model

A `.gsz` needs no mapping wizard. Its regions already know which material they
are, its piezometric line already knows it carries water, and — when the file was
saved after solving — it carries SLOPE/W's own answer on every trial surface it
tried. Only one question is left: which analysis to import, because a GeoStudio
file routinely holds several over one geometry and they can differ in materials as
well as in slip surface.

**Where to get the file.** Part 2 runs on problem §2.25 of Seequent's
[GeoStudio slope stability verification manual](https://files.seequent.com/GeoStudio/Manuals/Slope%20Stability%20Verification%20Manual.pdf),
Baker and Leshchinsky's safety-map earth dam. Seequent publishes the manual and
the models behind it, and the file downloads as
[Baker and Leschinsky - Earth Dam.gsz](https://files.seequent.com/GeoStudio/SlopeW/Baker%20and%20Leschinsky%20-%20Earth%20Dam.gsz)
from Seequent's own file server. It belongs to Seequent, not to us, and is not
redistributed here — download it from that link.

The dam has a diamond clay core inside a granular fill on a hard base, a reservoir
standing part way up one face, a phreatic surface dropping through the core to the
opposite toe, and a 5 m cracked layer at the crest modeled as a dry tension crack.

Once the file is downloaded, the import runs from the same menu. Click
**File → Import GeoStudio (SLOPE/W)…** and choose the `.gsz`. Because this file
holds two analyses, the picker appears:

![The GeoStudio analysis picker: the file's two Spencer analyses, with their type and method](images/w02_gsz_analyses.png){width=640}

Take the first, **Spencer**, and click **OK**. (A file holding one analysis skips
this step — there is nothing to choose.)

### What the import reports

The model lands on the canvas and Studio lists what it could not carry across
cleanly:

![The notes reported by the GeoStudio import: four of them, on the imported circle, the reservoir, the tension crack and the unit system](images/w02_gsz_notes.png){width=740}

Read those four against what the problem defines:

- **The failure surface came from one of SLOPE/W's own trial circles**, not from a
  search — and the note gives SLOPE/W's factor of safety on it, 1.934 by Spencer.
  It comes across so the two programs can be compared on one surface.
- **Water stands above the ground surface**, so the reservoir counts as ponded
  water. GeoStudio carries its weight implicitly and so does XSLOPE, from the
  imported piezometric line — nothing goes on the loads sheet, and adding a block
  by hand would count the reservoir twice.
- **The tension crack imported at depth 5.00**, which is the cracked crest layer
  the manual describes.
- **The model is metric.** XSLOPE is unit-agnostic and imports the numbers as they
  stand, so the answers come back in kN/m³, m and kPa.

Here is the imported dam:

![Studio after the GeoStudio import: the clay core inside the granular fill on the hard base, the piezometric line, the derived reservoir load and SLOPE/W's own circle](images/w02_gsz_window.png){width=1000}

Three material zones, the piezometric line through the core, the reservoir drawn
as the load the engine derives from it, and the one imported circle.

### Running it

We want the two programs on the same surface, so this run goes on the imported
circle rather than on a search. Click **Run LEM**, choose **Spencer**, set
**Analysis** to **Single surface** — which analyzes the first circle in the file,
and this file has one — and click **OK**.

XSLOPE returns **FS = 1.939** on that circle, against SLOPE/W's own 1.934 —
a difference of 0.3% on identical geometry, materials and water. That agreement is
already published here: see
[§2.25 of the GeoStudio verification page](../verification/geostudio.md#gs-2-25),
which reaches 1.939 on the same circle from an input file built by hand from the
manual.

### Check its work

- **The materials match the manual's table.** Clay core c′ = 20 kPa, φ′ = 20°,
  γ = 20 kN/m³; granular fill c′ = 0, φ′ = 40°, γ = 21.5; hard base c′ = 200,
  φ′ = 45°, γ = 24.
- **The reservoir is carried once.** The loads sheet is empty and water loads are
  automatic — the note says so, and the canvas shows the derived block.
- **The imported circle is SLOPE/W's, not a search result.** A free search on the
  same model finds its own minimum at 1.925, below the 1.939 on SLOPE/W's circle.

### What does not come across

Material zones, strengths, water conditions, surcharge and line loads, tension
cracks and reinforcement all import. Three things do not, and each is reported
rather than dropped in silence. SLOPE/W's search definition has no XSLOPE
equivalent, so a file saved without a solved surface arrives with no failure
surface at all — and a critical surface SLOPE/W found that is not a circle cannot
be taken as a trial circle either. Strength and pore-pressure options XSLOPE does
not carry come in as whatever fits, named in the notes: an unsupported strength
model as Mohr-Coulomb, GeoStudio's impenetrable Bedrock as an elastic material,
and SLOPE/W's `Ru` pore pressures as zero — that last one under a note saying the
factor of safety will be wrong. And a probabilistic analysis brings its standard
deviations across but not the correlations or the truncated sampling ranges
SLOPE/W applies to them.

---

## Part 3 — A Slide2 model

A Slide2 model behaves like a `.gsz`: it already knows what its geometry means, so
it prompts only for the scenario to take. A `.slmd` bundles several — a base case
plus variants — the way a GeoStudio file bundles analyses.

**Where to get the file.** Rocscience publishes its
[Slide2 verification manuals](https://www.rocscience.com/help/slide2/verification-theory/verification-manuals)
as PDFs, but not the 111 verification models behind them; the download links on
that page cover two hand-calculation examples only. Problem 104 can still be
followed all the way, because the manual builds it on Slide2's own *Tutorial 28 —
Seismic Analysis with the Newmark Method*, and the tutorial models are a public
download: [the Slide2 tutorials page](https://www.rocscience.com/help/slide2/tutorials)
links the older tutorial set as a zip, and `Tutorial 28 Seismic.slmd` sits in its
Model Files folder. That file belongs to Rocscience, not to us, and is not
redistributed here.

It models a 10 m slope at 2:1 in three layers over a horizontal base, dry, with no
external loads — soil 1 at c = 0 and φ = 38°, soil 2 at c = 5.3 kPa and φ = 23°,
soil 3 at c = 7.2 kPa and φ = 20°, all at γ = 19.5 kN/m³.

With the zip unpacked, we start the import. Click **File → Import Slide2…** and
choose `Tutorial 28 Seismic.slmd`. The scenario picker lists what the file holds:

![The Slide2 scenario picker: the master scenario and the four the tutorial builds on it](images/w02_slide2_scenarios.png){width=640}

Take **No Seismic** and click **OK**.

### What the import reports

Three notes come back. The model is metric, as Part 2's was. Slide2 was set to run
Spencer, so we solve with Spencer too and the comparison stays like-for-like. And
the third note carries the one difference that matters between this importer and
the GeoStudio one:

> no failure surface was imported — this scenario defines a SEARCH (grid, block,
> path or metaheuristic), which has no xslope equivalent, and no specified circle
> or surface to take one from.

Slide2 tutorial models define a search rather than a surface, and a search
definition does not translate. So the model arrives complete except for the one
thing a stability run cannot do without, and a starting circle has to be added. We
open the circles editor and add a row: **Xo** = 40, **Yo** = 45, **Option** =
`Intercept`, **Xi** = 30, **Yi** = 25. That places the center over the middle of
the face and two slope heights above the toe, on a circle through the toe itself —
the placement any hand-built model would start from, and the one
[LEM-3](lem03_layered_slope.md) explains.

![Studio after the Slide2 import, with the starting circle added: three soil layers over a horizontal base](images/w02_slide2_window.png){width=1000}

### Running it

The model now has everything a search needs. Click **Run LEM**, choose **Spencer**
and **Auto search**, leave the slice count at 40, and click **OK**.

The search returns **FS = 1.372**. Slide2's own published answer for problem 104 is
1.360, a difference of 0.9%, and the two disagree in a way the verification entry
explains: Slide2's Surface Altering optimization refines the surface away from a
circle and so reaches slightly lower than a circular search can. The same 1.372 is
what
[problem 104 of the Slide2 verification page](../verification/rocscience.md#vp104)
publishes, from an input file built by hand.

### Check its work

- **Three materials, three zones**, in the order the manual lists them, all at the
  same unit weight.
- **The model is dry.** Pore pressure reads `none` on every material, which is
  what the source defines — an import that invented a water table would be wrong
  here.
- **The starting circle is yours, not Slide2's.** The import says so; a run
  launched without adding one has no surface to solve.
- **The critical surface comes out a toe circle.** It enters at the toe, exits
  just behind the crest break, and bottoms about 4 m clear of the horizontal base.

---

## Conclusion

This tutorial covered:

- A DXF drawing imported layer by layer, with each layer mapped to a material
  zone, a piezometric line, a load or a set of starting circles — and the
  properties no CAD file can hold typed in afterward.
- A GeoStudio SLOPE/W model imported whole, including one of SLOPE/W's own solved
  trial circles, which puts both programs on one surface: 1.939 against 1.934.
- A Slide2 scenario imported the same way, arriving without a failure surface
  because Slide2 stores a search rather than a surface, and searching to 1.372
  against Slide2's published 1.360 once a starting circle is added.
- The notes every import reports, read against what the source model defines —
  because those notes name the one thing that would otherwise make an imported
  answer quietly wrong.

**Where to go next:** [DXF Import/Export](../usage/dxf.md) and
[GeoStudio Import/Export](../usage/geostudio.md) carry the full layer and mapping
tables and the script route for both. The
[GeoStudio](../verification/geostudio.md) and
[Slide2](../verification/rocscience.md) verification pages hold the corpora the
two vendor models here come from, problem by problem.
