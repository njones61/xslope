---
title: "Tutorial LEM-4 — Water in the Slope"
description: "Build a three-layer slope with a piezometric line in XSLOPE, hold one deep circle fixed, and measure what the water is worth: pore pressure on the slice bases, the γ/γ_sat unit-weight split, the r_u alternative, and a parametric sweep on a saturated unit weight."
---

# Tutorial LEM-4 — Water in the Slope

A 44 ft slope in three soil layers, with a piezometric line falling from
elevation 80 behind the crest to the toe. Every strength here is an **effective
stress** strength — a c′ and a φ′ — so every one of them is read against the
pore pressure on the failure surface, and the water is the input the answer
turns on: the same slope, on the same circle, reads **2.471 dry and 1.579
wet.**

![Slope with three layers and a piezometric line](../lem/sample_images/method_slices_problem.png){width=700}

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Assistant</span><p>~5 min</p></div>
<div class="tgt-tile"><span class="tg-label">By hand</span><p>15–20 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Build a layered section whose materials read their pore
pressure from a piezometric line, hold one specified deep circle fixed while
the water assumptions change around it, and learn what each input means: the
line and the u it produces on every slice base, the two unit-weight columns,
the r<sub>u</sub> alternative, and a fixed-surface parametric sweep on a
saturated unit weight.
</div>
<p><span class="tg-pill">three materials</span><span class="tg-pill">piezometric line</span><span class="tg-pill">effective stress</span><span class="tg-pill">single surface</span><span class="tg-pill">saturated unit weight</span><span class="tg-pill">sensitivity sweep</span></p>
<div class="tgm-model" markdown>**Completed model** — [xslope_method_slices_problem.xlsx](../lem/files/xslope_method_slices_problem.xlsx) — the same file used by [LEM Sample Problem 5](../lem/samples.md#5-slope-with-multiple-materials-and-piezometric-line)</div>
</div>

---

## The problem

**Materials** — three Mohr-Coulomb (`mc`) soils, listed top down, every one an
effective-stress strength reading the piezometric line:

| Mat ID | Name | γ (pcf) | c′ (psf) | φ′ (deg) | Pore pressure |
|---|---|---:|---:|---:|---|
| 1 | `soil 1` | 130 | 200 | 28 | `piezo` |
| 2 | `soil 2` | 120 | 100 | 32 | `piezo` |
| 3 | `soil 3` | 132 | 400 | 27 | `piezo` |

**Geometry** — one profile line per material, and the maximum depth at `0`,
the elevation of the rigid base:

| Line | Material | Points (x, y) |
|---|---|---|
| 1 | 1 `soil 1` | (0, 84), (150, 84), (174.7, 64) |
| 2 | 2 `soil 2` | (0, 64), (174.7, 64), (204.3, 40) |
| 3 | 3 `soil 3` | (0, 40), (320, 40) |

The ground surface is the top edge of the section the three lines make
together: level crest at 84, a face falling through the break at (174.7, 64)
to the toe at (204.3, 40), then flat ground at 40 running out to x = 320. By
the [top-of-a-layer rule](lem03_layered_slope.md#the-problem), lines 2 and 3
are the contacts — soil 1 is 20 ft thick, soil 2 is 24 ft, and soil 3 runs 40
ft down to the base.

**Water** — a piezometric line of eight points, spanning the full section:

| x | 0 | 75 | 112 | 140 | 170 | 189.5 | 204.3 | 320 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| y | 80 | 79 | 76 | 70 | 61 | 52 | 40 | 40 |

It starts 4 ft below the crest, steepens as it approaches the face, and comes
down to the toe — from x = 189.5 on, it lies *on* the lower face and the flat
ground beyond it. That stretch is where the water reaches the surface, and it
decides what the search finds below.

**The failure surface** — one circle, specified rather than searched for:
center (195, 150), **Depth** = `18.1`. This problem is a single-surface
exercise: the circle is the surface of interest, entered as given, and
everything this page measures is measured on it.

---

## Choose how you want to build it

Three ways to get this model into XSLOPE. They produce the same model — pick
one and skip the other two.

1. **[Build it with the AI assistant](#a-building-it-with-the-ai-assistant)** —
   describe the layers and the water line in a sentence and audit what it
   entered.
2. **[Build the Excel input file](#b-building-the-excel-file)** — four
   worksheets, one of them the piezometric line.
3. **[Build it in Studio](#c-building-it-in-studio)** — the editors, with the
   line drawing itself across the section as you type.

Whichever you choose, rejoin at [Running the analysis](#running-the-analysis).

---

## A — Building it with the AI assistant {#a-building-it-with-the-ai-assistant}

The drawing at the top of this page carries the layers, the strengths and the
water line, but not the coordinates. Paste it into the chat box with the
numbers, or describe the whole model:

<div class="prompt-block" markdown>
```text
Build a model with three horizontal soil layers over a rigid base at elevation 0. Layer 1: top at elevation 84, 130 pcf, c = 200 psf, phi = 28. Layer 2: top at 64, 120 pcf, c = 100 psf, phi = 32. Layer 3: top at 40, 132 pcf, c = 400 psf, phi = 27. The ground surface runs level at 84 from x = 0 to 150, drops through (174.7, 64) to a toe at (204.3, 40), then runs level at 40 out to x = 320. All three soils take their pore pressure from a piezometric line through (0, 80), (75, 79), (112, 76), (140, 70), (170, 61), (189.5, 52), (204.3, 40), (320, 40). Add one starting circle centered at (195, 150) with the bottom of the circle at elevation 18.1.
```
</div>

### Check its work

- **Three materials, top down**, and every one of them `u` = `piezo`. The
  strengths are effective-stress pairs — a c′ *and* a φ′ — and an
  effective-stress strength without a pore-pressure source reads u = 0
  everywhere. If any material came back with `none`, say: *"All three soils
  read the piezometric line — set every pore pressure option to piezo."*
- **One piezometric line, eight points, spanning x = 0 to 320.** A line that
  stops short of the surface's reach is an error at run time, not a warning.
- **The line never rises above the ground.** Here it touches the ground at the
  toe and runs along it. A piezometric line *above* the ground surface means
  standing water, and the weight of standing water is a distributed load —
  XSLOPE derives it automatically when `main!D23` **Water loads** is `auto` —
  not something to skip.
- **γ_sat is blank.** One unit weight per material is this exercise's
  assumption; the second column is an edit this page makes
  [further down](#the-two-unit-weight-columns), deliberately.
- **The maximum depth is `0`** — the elevation of the rigid base, not a
  thickness.
- **The circle's Depth is `18.1`**, an *elevation*: the loader forms
  R = Yo − Depth = 131.9. A circle entered with a 18.1 ft radius is a
  different, far smaller surface.
- **Units are declared Imperial.**

Continue at [Running the analysis](#running-the-analysis).

---

## B — Building the Excel file {#b-building-the-excel-file}

Start from [input_template.xlsx](../inputs/input_template.xlsx) and save a copy
under a name of your own. Confirm `main!D8` **Units** is `Imperial`, confirm
`main!D23` **Water loads** is `auto`, and leave the rest of that sheet alone.

### 1. The `mat` worksheet

One row per material, and the row number is the ID:

1. `mat!A11` = `1`, `mat!B11` = `soil 1`, `mat!C11` **g** = `130`,
   `mat!E11` **option** = `mc`, `mat!F11` **c** = `200`, `mat!G11` **f** =
   `28`, `mat!O11` **u** = `piezo`.
2. `mat!A12` = `2`, `mat!B12` = `soil 2`, g = `120`, option = `mc`, c = `100`,
   f = `32`, u = `piezo`.
3. `mat!A13` = `3`, `mat!B13` = `soil 3`, g = `132`, option = `mc`, c = `400`,
   f = `27`, u = `piezo`.

![The finished mat worksheet](images/lem04_sheet_mat.png)

Two of these columns are the subject of this page. **`u` is where each
material's pore pressure comes from** — `piezo` reads the piezometric line,
`none` reads nothing, and the [other two options](#r_u-a-pore-pressure-without-a-line)
appear later. And the pair **C:D** are the unit-weight columns: **g** is the
unit weight this model uses everywhere, and **gsat**, left blank here, is the
saturated unit weight the engine would use below the water table if it were
given one. Blank means *g throughout* — which is this problem's assumption.

### 2. The `profile` worksheet

1. `profile!B2` **Max Depth** = `0`. This is an *elevation*.
2. Profile Line #1 — `profile!B5` **Mat ID** = `1`, then `0, 84` / `150, 84` /
   `174.7, 64` in the `x` / `y` columns.
3. Profile Line #2 — `profile!E5` **Mat ID** = `2`, then `0, 64` /
   `174.7, 64` / `204.3, 40`.
4. Profile Line #3 — `profile!H5` **Mat ID** = `3`, then `0, 40` and
   `320, 40`.

![The finished profile worksheet](images/lem04_sheet_profile.png)

Lines 1 and 2 each carry a stretch of the face, which is why they take three
points; line 3 is the flat top of the bottom layer, and two points draw it.

### 3. The `piezo` worksheet

`piezo!B3` **Type:** = `piezo`, then the eight points in the `x` / `y`
columns, left to right: `0, 80` / `75, 79` / `112, 76` / `140, 70` /
`170, 61` / `189.5, 52` / `204.3, 40` / `320, 40`.

![The finished piezo worksheet](images/lem04_sheet_piezo.png)

A piezometric line is the model's water table, given as its own polyline: the
pore pressure at a point on the failure surface is γ<sub>w</sub> times the
vertical distance below the line, and zero above it. The sheet's own legend
(columns G) lists the two **Type** settings: `piezo` is that static-head
reading, and `phreatic` reduces the head under an inclined stretch of line by
cos² of its local slope — the correction for steady seepage flowing parallel
to the line. This problem's source states heads, so `piezo` is the right
reading. **Piezometric Line 2 stays empty**; it is the second stage of a rapid
drawdown analysis, which this model is not.

### 4. The `circles` worksheet

One circle, as the problem specifies it: `circles!B3` **Xo** = `195`,
`circles!C3` **Yo** = `150`, `circles!D3` **Option** = `Depth`,
`circles!E3` **Depth** = `18.1` — the elevation of the circle's lowest point,
18.1 ft above the rigid base and 21.9 ft below the toe.

![The finished circles worksheet](images/lem04_sheet_circles.png)

Save the file and continue at [Running the analysis](#running-the-analysis).

---

## C — Building it in Studio {#c-building-it-in-studio}

Start with **File → New** and work down the **Inputs** tree.

### 1. Materials and profile lines

Click **Materials**, and in **Table view** press **Add row** three times for
the three soils in the order the table above lists them — the row order fixes
the Mat IDs. Set every row's **u** to `piezo`. Then **Profile lines**: set
**Max depth (bottom boundary elevation):** to `0`, and press **Add line**
three times, each line taking its **Material:** and its vertices from the
geometry table above. The mechanics are
[LEM-3's](lem03_layered_slope.md#c-building-it-in-studio), one line more.

### 2. The piezometric line

Click **Piezometric lines**. The editor holds the two lines on tabs — the
water table on **Line 1**, and **Line 2 (rapid drawdown)**, which stays empty
here. Press **Add row** eight times and enter the points from the table above,
left to right.

![The piezometric-lines editor holding the water table](images/lem04_studio_piezo.png)

The preview draws the line on the section as you type, which is where a
mistyped elevation shows up: the points should fall smoothly from 80 at the
left edge to 40 at the toe, not jump. The help line above the table states the
one rule the table has — points ordered left to right. Click **OK**.

### 3. The circle

Click **Circles**, press **Add row**, and enter **Xo** `195`, **Yo** `150`,
**Option** `Depth`, **Depth** `18.1`. This page does not press the generator:
the circle is *specified* — the problem names it — and the surface of a
single-surface study is an input, not a guess to be refined. The preview draws
it bottoming out 21.9 ft below the toe elevation. Click **OK**.

Continue below.

## Running the analysis

However you built it, you now hold the same model:

![The finished model](images/lem04_inputs.png){width=1000}

The three profile lines are drawn in their materials' colors, the heavy blue
line is the piezometric line — note it merging with the ground surface from
the lower face outward — and the red dashed arc is the circle, its center at
(195, 150) above the frame.

Click **Run LEM…** and choose **Method** = `Spencer` and **Analysis** =
`Auto search`, with the slice count left at 40. **Surface** reads `Circular`
as a fixed label — the model defines circles and no non-circular surface.
From Python:

```python
from xslope.fileio import load_slope_data
from xslope.search import circular_search
from xslope.plot import plot_solution

sd = load_slope_data("my_water_slope.xlsx")
fs_cache, _, path, circles = circular_search(sd, "spencer", num_slices=40)
crit = fs_cache[0]
plot_solution(sd, crit["slices"], crit["failure_surface"], crit["solver_result"])
```

---

## Exploring the results

### What the search finds on a wet slope

![The circular search](images/lem04_search.png){width=1000}

**FS = 0.764** — and it is not the circle you entered. The search started at
the deep seed, walked its center down and to the right, and settled on a
sliver of the lower face: 30 ft of section from x = 174, just under the slope
break, to the toe at 204.3, nowhere more than 8.1 ft deep. Its base is 41.7 ft
long and 41.1 ft of that is in soil 2; the mass above it weighs 20,143 lb/ft —
3% of what the deep circle carries.

The face is where this model is weakest, and the water is why. Soil 2 has the
least cohesion of the three soils, and from x = 189.5 to the toe the
piezometric line lies *on* the ground surface: the lower face is a seepage
face, with pore pressures up to 472 psf on a surface skimming a few feet
below it. Effective stress near that face is close to zero, so what holds it
up is c′ = 100 psf — and 100 psf cannot hold it. A search ranks every circle
it tries by factor of safety and nothing else, so a shallow failure in a weak,
wet face outranks every deep mechanism the section has. The result is honest:
this model really does say the lower face sloughs. It is also not the
question this problem asks.

The question is the specified deep circle, and the way to ask it is to stop
searching: run again with **Analysis** = `Single surface`, which solves the
first circle exactly as entered. Holding the surface fixed is what makes the
rest of this page possible — change the water and re-search, and every answer
comes back about whatever surface the search picked that time; change the
water on a pinned surface, and the difference *is* the water. From Python:

```python
from xslope.slice import generate_slices
from xslope.solve import solve_selected

ok, (slice_df, surface) = generate_slices(sd, circle=sd["circles"][0],
                                          num_slices=40)
result = solve_selected("spencer", slice_df)
plot_solution(sd, slice_df, surface, result)
```

### The pinned circle, dry

Start by asking what the circle would be worth with no water at all. Set every
material's pore pressure option to `none` — leave the piezometric line where
it is; a line no material reads contributes nothing:

- **Studio** — open **Materials** and set each row's **u** to `none`.
- **Excel** — `mat!O11`, `mat!O12`, `mat!O13` = `none`.
- **Assistant** — say: *"Set every material's pore pressure option to none."*

Run the single-surface analysis:

![Spencer on the pinned circle, dry](images/lem04_solution_dry.png){width=1000}

**FS = 2.471.** The circle enters the crest at x = 80.8, bottoms out at
elevation 18.1 in soil 3, and exits on the flat ground at x = 267.8 — 215 ft
of base, 154 ft of it in soil 3, carrying 717,261 lb/ft. The green bars drawn
on the base are the legend's **Eff Normal Stress (σ′)**: with u = 0
everywhere, the effective normal stress *is* the total normal stress, and
every pound of it earns frictional strength at tan φ′.

### The piezometric line, and what it costs

Now put the water back — the same three cells, returned to `piezo` (or undo)
— and run the same surface again:

![Spencer on the pinned circle, wet](images/lem04_solution_wet.png){width=1000}

**FS = 1.579.** Same circle, same 215 ft of base, same 717,261 lb/ft of soil —
the water table takes 36% off the factor of safety. The figure shows where it
went: the blue bars are the legend's **Pore Pressure (u)**, and on 39 of the
40 slice bases they are not zero. Each one is the piezometric reading —
u = γ<sub>w</sub> × the vertical distance from the base to the line — and the
largest, 2571 psf at x = 157, is a base sitting 41.2 ft below the line
(62.4 × 41.2). Only the first slice, whose base hugs the crest above the
line's left end, reads zero.

The mechanics are one line long: τ = c′ + (σ − u) tan φ′. The water carries
part of the normal load, the part it carries produces no friction, and on
this circle the pore pressure bars visibly hollow out the effective-stress
bars that were there in the dry run. Nothing about the soil changed — c′ and
φ′ are what they were — and the driving weight did not change either. The
strength did.

Run the other methods on both states and the pattern holds everywhere:

| | OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|---|---:|---:|---:|---:|---:|---:|---:|
| Dry (u = `none`) | 2.212 | 2.477 | 2.411 | 2.689 | 2.561 | 2.471 | 2.471 |
| Piezometric line | 1.303 | 1.576 | 1.533 | 1.766 | 1.641 | 1.579 | 1.579 |

Two readings worth taking. **Bishop, Spencer and Morgenstern-Price agree to
0.2% wet or dry** — on a circular surface with φ > 0 they are the reliable
cluster, and Spencer is the one to report. And **the water widens the
spread**: OMS sits 10% under Spencer dry but 17% under it wet. OMS computes
each slice's normal force from the weight alone and then subtracts the full
u·dl from it, and on a deep circle under high pore pressure that
approximation sheds too much normal force — a known conservatism of the
method, at its worst exactly here.

### The two unit-weight columns {#the-two-unit-weight-columns}

The `mat` sheet carries two unit weights per material and this model filled
one. **g** is the total unit weight of the soil as it sits above the water
table — solids plus whatever moisture it holds. **gsat** is the total unit
weight when the voids are full, which is what the soil below the water table
actually weighs; a saturated soil is typically a few pcf heavier. With
**gsat** blank, XSLOPE uses **g** for the whole column of every slice —
that is this exercise's assumption, one weight per soil.

Fill the column in and the engine splits the weight: each slice's column uses
γ above the water table and γ_sat below it, with the piezometric line serving
as the water table. Enter `mat!D11` = `135`, `mat!D12` = `125`,
`mat!D13` = `140` (in Studio, the **gsat** column of the materials table) and
re-run the pinned circle:

| | ΣW (lb/ft) | FS (Spencer) |
|---|---:|---:|
| γ only (gsat blank) | 717,261 | 1.579 |
| γ above the line, γ_sat below | 747,092 | 1.621 |

The sliding mass gains 4.2% of its weight, all of it below the water table —
and the factor of safety goes *up* 2.6%. That direction surprises most people
once: added weight adds driving moment, but it also adds normal stress on a
base whose pore pressure is fixed by the line, and on this deep, flat-bottomed
circle the frictional gain outruns the extra drive. Neither effect is large.

What γ_sat is **not** is a buoyant unit weight. Both columns are *total* unit
weights, and the water's effects enter the analysis explicitly — the weight of
the pore water through γ_sat, and its pressure through u — never by
pre-subtracting them into an effective weight. A γ_sat entered as
γ − 62.4 builds a model that is wrong twice: the slices weigh too little,
and u is still subtracted from a normal stress that no longer contains the
water it removes. If the two effects are to cancel, the analysis is where
they cancel, not the input file.

The split is worth having when the answer is worth refining — the two rows
above bound what it buys on this section. Leave **gsat** blank and the model
is the exercise as published; this page leaves it filled for one more section,
because a filled column is something a sweep can hold.

### r_u — a pore pressure without a line {#r_u-a-pore-pressure-without-a-line}

The `u` column has two more settings. `seep` reads pore pressures from a
finite element seepage solution — the subject of the
[seepage documentation](../seep/overview.md), and the step past this page
when the flow field matters. `ru` needs no companion at all: it computes
u = r<sub>u</sub> × σ<sub>v</sub> on each slice base, where σ<sub>v</sub> is
the vertical stress of the soil column above it and r<sub>u</sub> is a
per-material ratio entered in the column beside `u` (`mat!P`). It is the
classic chart parameter — Bishop and Morgenstern's stability coefficients are
tabulated against it — and it survives as a quick global estimate
("pore pressures run about a third of overburden") where no line has been
drawn. Note what it couples: an r<sub>u</sub> pore pressure *scales with the
soil weight*, where a piezometric line's u is pure geometry — the two models
answer weight questions differently, which is one reason the line, when you
have one, is the better input.

### Does the γ_sat guess matter?

The 135 / 125 / 140 entered above are estimates, and an input you had to
guess is an input to sweep. The question is narrow — *on this circle, how much
does the answer move as one saturated unit weight moves?* — and narrow is
exactly what a fixed-surface sweep answers.

Click **Run → Parametric…** and set it up:

![The Parametric dialog set up for the γ_sat sweep](images/lem04_studio_parametric.png)

1. **Mode** = `Sensitivity (tornado + plots)`, **Method** = `Spencer`.
2. Under **Parameter**, **Material** = `soil 3` and **Property** =
   `gamma_sat` — the deep layer, the one 154 ft of this base runs through.
3. **Default ±%** = `5`, **Points** = `7` — a ±7 pcf band about 140, about
   the width of a defensible guess. Press **Add parameter**, and the
   parameter lands in the table.
4. **Untick "Re-search the critical surface at each step."** This is the
   opposite of the advice in [LEM-2](lem02_loads_on_the_crest.md), and it is
   not a shortcut — it is the question. Re-solving the entered surface is
   what a prescribed-surface study means.
5. **Run**, then double-click the tornado's one bar — the dialog's own note
   says so: *"Double-click a bar for that parameter's curve."*

From Python, the same sweep:

```python
from xslope.sensitivity import sensitivity
from xslope.plot import plot_sensitivity

success, result = sensitivity(sd, param="mat:soil 3:gamma_sat",
                              rel_range=0.05, n=7, search=False,
                              methods=("spencer",))
plot_sensitivity(result["df"])
```

![The γ_sat sweep on the fixed circle](images/lem04_sweep.png){width=1000}

The curve is a straight line from **1.576 at 133 pcf to 1.666 at 147 pcf** —
0.006 of factor of safety per pcf, the base case marked at (140, 1.62). Read
it against the section above it: the entire ±5% band of saturated-unit-weight
guesses is worth about ±3% of the answer, on the same slope where the
existence of the water table was worth 36%. The sweep is the license to stop
agonizing over the third digit of γ_sat — and its direction is not a law of
nature but a property of *this* circle: sweep soil 2's γ_sat instead and the
answer creeps the other way, 1.625 down to 1.614 across 120–135 pcf, because
that soil's weight sits over the toe end of the arc. Unit weights move a
wet-slope answer by little, in whichever direction the geometry decides.

One honest wrinkle closes the loop. Leave **Re-search the critical surface at
each step** ticked here and the sweep returns a flat line at **0.776** for
every value of γ_sat: each step re-finds the face sliver, whose soil 1 and
soil 2 mass never contains soil 3's saturated weight at all. The curve is
correct, and it is a curve about the wrong surface — the checkbox is the
difference between *"how sensitive is my circle"* and *"how sensitive is
whatever the search finds."* (With the γ_sat column filled, the sliver reads
0.776 rather than the 0.764 of the file as shipped — the splash of extra
saturated weight reaches even it.)

Before saving, put the model back the way the exercise defines it: clear the
three **gsat** cells (or undo), and confirm the three `u` cells read `piezo`.

---

## Conclusion

This tutorial demonstrated:

- A piezometric line as the model's water table — its own polyline, spanning
  the section, entered once and read by every material whose `u` option is
  `piezo` — and pore pressure as γ<sub>w</sub> times depth below it, drawn on
  every slice base of the solution plot.
- Effective stress doing the work: the same circle reads 2.471 dry and 1.579
  wet, with c′, φ′ and the driving weight identical — the water table costs
  36% of the factor of safety through u alone.
- A search on a wet slope honestly reporting a shallow seepage-face failure
  at 0.764, and the single-surface analysis as the way to hold a specified
  deep circle fixed while the water assumptions change.
- The two unit-weight columns: **g** above the water table, **gsat** below it
  (blank = **g** throughout), both *total* unit weights — never buoyant —
  worth +2.6% on this circle when filled with realistic values.
- r<sub>u</sub> as the chart-era alternative — u = r<sub>u</sub> σ<sub>v</sub>
  per material, no line required, scaling with weight where a line does not.
- A fixed-surface parametric sweep: re-search unticked because the surface is
  prescribed, a ±5% γ_sat band moving Spencer by about ±3% — and the same
  sweep with re-search left on answering a different question entirely, flat
  at the sliver's 0.776.

**Where to go next:** the [tutorials index](index.md) lists the rest of the
series. [Sample Problem 5](../lem/samples.md#5-slope-with-multiple-materials-and-piezometric-line)
carries this same model into a multi-seed search of the sloughing failure the
single seed found here, and the
[seepage documentation](../seep/overview.md) replaces the hand-drawn line
with a computed flow field — the `seep` setting this page's `u` column left
unused.
