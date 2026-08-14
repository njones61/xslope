---
title: "Tutorial LEM-10 — Finding the Global Minimum"
description: "Search a cohesionless embankment on soft clay in XSLOPE from two different starting circles and get two different answers — a surficial sliver on the sand face and the deep foundation failure — then use grid seeding and a minimum slip depth to find the one a design is sized for."
---

# Tutorial LEM-10 — Finding the Global Minimum

A 15 ft cohesionless embankment on 20 ft of soft clay. Two mechanisms compete
here — a shallow slide in the sand face and a deep one through the foundation —
and an automated search returns whichever of them its starting circle is near.
**The lower of the two numbers is not the answer.**

![A cohesionless embankment on a soft clay foundation](../lem/sample_images/mult_min_inputs1.png){width=700}

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Assistant</span><p>~5 min</p></div>
<div class="tgt-tile"><span class="tg-label">By hand</span><p>10–15 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Build a section that holds two competing failure mechanisms,
search it from a circle in the embankment and from a circle in the foundation,
and read why the two searches disagree — then find the deep mechanism with grid
seeding and a minimum slip depth.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">profile lines</span><span class="tg-pill">starting circles</span><span class="tg-pill">circular search</span><span class="tg-pill">grid seeding</span><span class="tg-pill">minimum slip depth</span></p>
<div class="tgm-model" markdown>**Completed model** — [xslope_mult_min_KEY.xlsx](../lem/files/xslope_mult_min_KEY.xlsx) — the same file used by [LEM Sample Problem 13](../lem/samples.md#13-multiple-local-minima)</div>
</div>

---

## The problem

**Materials** — a drained cohesionless fill on an undrained clay:

| Mat ID | Name | γ (pcf) | c (psf) | φ (deg) | Pore pressure |
|---|---|:---:|:---:|:---:|---|
| 1 | `embankment` | 120 | 0 | 30 | `none` |
| 2 | `foundation` | 120 | 450 | 0 | `none` |

**c = 0 in the fill.** Nothing in the embankment resists a shallow slide except
friction along its base, so a slab of sand of any thickness — including a
vanishingly thin one — has the same factor of safety against sliding down the
face.

**Geometry** — two profile lines, each the *top* of its material layer, as in
[LEM-3](lem03_layered_slope.md#the-problem). Maximum depth = `-20`, the elevation
of the bottom of the soft clay.

**Profile Line 1 — material 1 (`embankment`):**

| x (ft) | y (ft) |
|:---:|:---:|
| 0 | 0 |
| 33.75 | 15 |
| 83.75 | 15 |

**Profile Line 2 — material 2 (`foundation`):**

| x (ft) | y (ft) |
|:---:|:---:|
| -50 | 0 |
| 83.75 | 0 |

Line 1 is the ground surface: toe, crest break at the top of the 2.25:1 face,
back of the crest. Line 2 is the contact, carried 50 ft past the toe so a deep
surface has ground to come up on.

**Starting circle** — one, centered above the face and reaching the bottom of the
foundation:

| Xo | Yo | Option | Depth |
|:---:|:---:|---|:---:|
| 11.25 | 40 | Depth | -20 |

The tables are the model, and each is laid out exactly as its destination is —
the template's worksheets and Studio's editors, same columns in the same order.
Select a table's block of values, copy, and paste it straight into the sheet or
editor rather than retyping it.

---

## Building the model

Three ways in, all producing the same file; they meet again at
[Running the analysis](#running-the-analysis).

**The AI assistant.** Paste the drawing at the top of this page into the chat box
and type `Build this model`, or describe it:

<div class="prompt-block" markdown>
```text
Build a model for a 15 ft embankment with a 2.25:1 face on a 20 ft soft clay foundation. The embankment is 120 pcf, cohesionless, with phi = 30; the foundation is 120 pcf with c = 450 psf, phi = 0. The crest runs 50 ft back from the top of the face, the ground continues 50 ft beyond the toe, and the bottom of the clay is at elevation -20. Add a starting circle tangent to the bottom of the clay.
```
</div>

Audit what it entered against the tables above, as
[LEM-3](lem03_layered_slope.md#a-building-it-with-the-ai-assistant) does: the fill
first so the Mat IDs match the profile lines, **c = 0 and φ = 30 in the fill**
(a drawing that says "sand" says neither), φ = 0 in the clay, and a maximum depth
of `-20` — an elevation, not a thickness.

**The Excel file.** Start from
[input_template.xlsx](../inputs/input_template.xlsx) and fill `mat`, `profile`
and `circles` from the tables above, as
[LEM-3's Excel path](lem03_layered_slope.md#b-building-the-excel-file) walks them.
`profile!B2` **Max Depth** = `-20` is the one cell worth checking twice.

**Studio.** Work down the **Inputs** tree — **Materials**, **Profile lines**,
**Circles** — entering the same tables in the same order, as in
[LEM-3's Studio path](lem03_layered_slope.md#c-building-it-in-studio).

---

## Running the analysis

However you built it, you now hold the same model:

![The finished model](images/lem10_inputs.png){width=1000}

The two profile lines are drawn in their materials' colors and the dashed red arc
is the starting circle, bottoming out on the hatched maximum depth at elevation
−20. Click **Run LEM…** and choose **Method** = `Spencer` and **Analysis** =
`Auto search`, with the slice count left at 40.

---

## Exploring the results

### Searching from a circle in the embankment

[LEM-3's rule](lem03_layered_slope.md#guarding-against-local-minima) is one
starting circle per layer, and Studio's **Generate starting circles…** builds
that set from the geometry. On this section it reports *"4 on the left-facing
face (toe at x = 0, height 15), one of them skimming its 24 degree cohesionless
face"*. Take the one that is the embankment's own mechanism — **Xo** `16.875`,
**Yo** `30`, **Option** `Depth`, **Depth** `0`, tangent to the top of the
foundation — and search from it alone:

![The search from the embankment circle](images/lem10_search_shallow.png){width=1000}

**FS = 1.299.** The grey arcs are the trial circles, fanning across the fill and
shrinking toward the face as the search walks its center up and left; the red
critical circle is the short mark high on the slope — 1.2 ft long, never more
than 0.01 ft below the ground surface, moving 1 lb/ft of sand.

That number is not an artifact of the search. On a cohesionless face the factor
of safety against an infinitely shallow slide is tan φ′ / tan β, which is
tan 30° / tan 23.96° = **1.299** for this 2.25:1 face. The search is converging on
a real limit — and the limit belongs to a slide with no mass in it. It is the
minimum of this model and not a mechanism anything is designed against.

### Searching from a circle in the foundation

The circle the file carries reaches the bottom of the clay. Searching from that
one instead, with everything else unchanged:

![Spencer on the deep foundation surface](images/lem10_solution_deep.png){width=1000}

**FS = 1.376**, on a circle centered at (16.90, 26.22) with a radius of 46.22 —
tangent to the limiting depth at elevation −20, entering at the crest and exiting
21 ft beyond the toe, 82.9 ft of failure surface carrying 204,041 lb/ft of soil.
Two searches on one model, 6% apart, on mechanisms that share almost nothing: the
sliver rides the sand face; this one cuts the full depth of the clay.

The starting circle is a seed and not an answer — solved as entered it gives
1.426, and the search deepens and lengthens it to 1.376.

### The methods on the seeded circle

All seven solution methods answer on the starting circle as entered. **This is a
single-surface run — no search**, so the geometry is identical in every column:

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.354 | 1.434 | 1.417 | 1.719 | 1.524 | 1.426 | 1.431 |

Spencer, Morgenstern-Price and Bishop land within 0.01 of each other, OMS 5%
below them and the Corps of Engineers method 21% above — interslice-force
assumptions on one surface, not a disagreement about which surface it is.

### Two tools for exactly this section

**Grid search (auto-seed the circular search)**, a checkbox in the **Run LEM…**
dialog, sweeps a grid of circle centers against a range of tangent elevations
before refining, instead of refining only the neighborhood of the circles on the
sheet. It removes the dependence on where the seeds were put, and here it returns
**1.327** on a face slide 11.8 ft long and 0.75 ft deep — a global sweep still
lands in the surficial family, because that family really does hold the model's
lowest surfaces.

**Ignore surficial (skin) failures**, with a **Min slip depth** beside it, is the
filter that removes them: a trial surface whose deepest point is less than that
far below the ground is rejected during the search. With grid seeding on and a
minimum slip depth of `5`, the search returns **1.376** on the deep circle — the
foundation mechanism, found without being seeded. From the embankment circle
alone the same 5 ft filter returns 1.455, on a 33 ft circle inside the fill, and
only at 8 ft reaches 1.376. The filter decides which surfaces are admissible; the
seeding decides which of them the search ever looks at.

### Which answer the model reports

Run every circle in one search and it reports the lowest surface any seed
reached — here, the 1.299 sliver. **Several starting circles, and agreement
between them, is the check**; where they disagree, the surfaces decide, and the
design number is the deep one.

---

## Conclusion

This tutorial demonstrated:

- A section with two genuine minima — a cohesionless embankment on soft clay —
  where the search returns whichever mechanism its starting circle sits nearest.
- The surficial limit a c = 0 face produces: **1.299** on a surface 1.2 ft long
  carrying 1 lb/ft, matching tan φ′ / tan β for the face angle — a real minimum
  and not a design case.
- The deep foundation failure at **1.376**, tangent to the bottom of the clay and
  moving 204,041 lb/ft, found by seeding a circle at that depth.
- **Grid search** as the seeding-independent sweep (**1.327** here, still
  surficial) and **Min slip depth** as the filter that rejects the skin, the two
  together returning **1.376**.
- Reading a search by the surface it converged on, since two answers on one model
  can be 6% apart and describe entirely different failures.

**Where to go next:** the [tutorials index](index.md) lists the series.
[Sample Problem 13](../lem/samples.md#13-multiple-local-minima) catalogues this
model, [LEM-3](lem03_layered_slope.md) is where the per-layer starting-circle rule
is built, and [Automated Search](../lem/search.md) documents the search itself.
