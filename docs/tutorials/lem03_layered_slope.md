---
title: "Tutorial LEM-3 — A Layered Slope"
description: "Build a two-layer slope in XSLOPE — a weak embankment on a stronger foundation over rigid rock — with a profile line at the material boundary and a starting circle at the base of each layer, then read why the critical surface settles on the contact between them."
---

# Tutorial LEM-3 — A Layered Slope

A 20 ft embankment with a 2:1 face, built on a 10 ft foundation layer over rigid
rock. Two soils, both undrained, and the fill is the weaker of them: c = 400 psf
against the foundation's 800. Nothing in the geometry says how deep the failure
surface will go — **the contact between the two layers decides that**, and the
model has to be built so the search can ask about both sides of it.

![Simple slope with multiple layers](../lem/sample_images/simple_mult_layers.png){width=700}

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Assistant</span><p>~5 min</p></div>
<div class="tgt-tile"><span class="tg-label">By hand</span><p>15–20 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Build a section with two material zones, put a profile line on
the boundary between them and a starting circle at the base of each, and read the
search result at depth: which layer the critical surface runs in, why it stops
where it does, and what would have to change for the deeper circle to be the
answer.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">profile lines</span><span class="tg-pill">per-layer starting circles</span><span class="tg-pill">generated circles</span><span class="tg-pill">circular search</span></p>
<div class="tgm-model" markdown>**Completed model** — [xslope_simple_mult_layers.xlsx](../lem/files/xslope_simple_mult_layers.xlsx) — the same file used by [LEM Sample Problem 3](../lem/samples.md#3-simple-slope-with-multiple-layers)</div>
</div>

---

## The problem

**Materials** — two Mohr-Coulomb (`mc`) soils, both undrained:

| Mat ID | Name | γ (pcf) | c (psf) | φ (deg) | Pore pressure |
|---|---|---:|---:|---:|---|
| 1 | `embankment` | 130 | 400 | 0 | `none` |
| 2 | `foundation` | 135 | 800 | 0 | `none` |

**Geometry** — one profile line per material, listed top down:

| Line | Material | Points (x, y) |
|---|---|---|
| 1 | 1 `embankment` | (0, 0), (40, 20), (90, 20) |
| 2 | 2 `foundation` | (−30, 0), (90, 0) |

Line 1 is the ground surface: the toe, the crest break at the top of the 2:1
face, and the back of the crest. Line 2 is the contact — the top of the
foundation, running the full width of the section. Maximum depth = `-10`, the
elevation of the rigid rock.

A profile line is the *top* of a material layer, so a two-layer section takes two
of them and the second one is the boundary between the layers. **Lines go where
the material changes, and nowhere else.** The toe and the crest break are corners
of the ground surface, not layer boundaries: they are vertices of line 1.

The ground runs 30 ft past the toe at elevation 0 because the foundation
continues there. A surface that dips into the foundation has to come back up
somewhere, and it comes up on that flat ground — a section that stopped at the
toe would have nowhere to put the deep mechanism this problem is about.

---

## Choose how you want to build it

Three ways to get this model into XSLOPE. They produce the same model — pick one
and skip the other two.

1. **[Build it with the AI assistant](#a-building-it-with-the-ai-assistant)** —
   describe both layers in a sentence and audit what it entered.
2. **[Build the Excel input file](#b-building-the-excel-file)** — three
   worksheets, a second row in each.
3. **[Build it in Studio](#c-building-it-in-studio)** — the editors, with the
   section redrawing beside them.

Whichever you choose, rejoin at [Running the analysis](#running-the-analysis).

---

## A — Building it with the AI assistant {#a-building-it-with-the-ai-assistant}

The drawing at the top of this page carries what the assistant needs — layer
thicknesses, the face inclination, both unit weights, both strengths. Paste it
into the chat box and type `Build this model`, or describe it:

<div class="prompt-block" markdown>
```text
Build a model for a 20 ft embankment with a 2:1 face on a 10 ft foundation layer over rigid rock. The embankment is 130 pcf with c = 400 psf, phi = 0; the foundation is 135 pcf with c = 800 psf, phi = 0. The crest runs 50 ft back from the top of the face, and the ground continues 30 ft beyond the toe. Add starting circles for a critical-surface search.
```
</div>

### Check its work

- **Two materials, the fill first.** The order fixes the Mat IDs, and the profile
  lines reference materials by ID — so an embankment entered second is a section
  built upside down.
- **φ = 0 in both.** The drawing gives a cohesion and no friction angle, which is
  an undrained strength; it never says so. If a friction angle appeared, say:
  *"Both strengths are undrained. Set phi to 0 in both materials and leave the
  pore pressure option at none."*
- **Two profile lines, one per material.** The ground surface on material 1, and
  the contact at y = 0 on material 2. A model with a single line has one layer,
  whatever the material table says.
- **The foundation's line spans the whole section**, from x = −30 to x = 90 —
  including the 30 ft beyond the toe, where the foundation is the ground surface.
- **The maximum depth is the elevation of the rock**, `-10`, not a thickness and
  not the toe elevation. If it came back at 0, say: *"The rigid base is 10 ft
  below the toe — set the maximum depth to -10."*
- **One starting circle at the base of each layer**, both centered above the
  middle of the face. If only one arrived, say: *"Add a second starting circle
  tangent to the base of the foundation, at the same center."*

Continue at [Running the analysis](#running-the-analysis).

---

## B — Building the Excel file {#b-building-the-excel-file}

Start from [input_template.xlsx](../inputs/input_template.xlsx) and save a copy
under a name of your own. Confirm `main!D8` **Units** is `Imperial` and leave the
rest of that sheet alone — no tension crack, no seismic coefficient, and the run
options blank.

### 1. The `mat` worksheet

One row per material, and the row number is the ID:

1. `mat!A11` = `1`, `mat!B11` = `embankment`, `mat!C11` **γ** = `130`,
   `mat!E11` **option** = `mc`, `mat!F11` **c** = `400`, `mat!G11` **φ** = `0`,
   `mat!O11` **u** = `none`.
2. `mat!A12` = `2`, `mat!B12` = `foundation`, `mat!C12` **γ** = `135`,
   `mat!E12` **option** = `mc`, `mat!F12` **c** = `800`, `mat!G12` **φ** = `0`,
   `mat!O12` **u** = `none`.

![The finished mat worksheet](images/lem03_sheet_mat.png)

**Rows are ordered top down.** Row 11 is the material at the surface and row 12
the one beneath it, which is the order the profile lines below assume.

### 2. The `profile` worksheet

1. `profile!B2` **Max Depth** = `-10`. This is an *elevation*: it puts the rigid
   base 10 ft below the toe.
2. Profile Line #1 — `profile!B5` **Mat ID** = `1`, then the three
   ground-surface points in the `x` / `y` columns: `0, 0` then `40, 20` then
   `90, 20`.
3. Profile Line #2 — `profile!E5` **Mat ID** = `2`, then two points: `-30, 0`
   and `90, 0`.

![The finished profile worksheet](images/lem03_sheet_profile.png)

The second line needs only its two endpoints because it is straight, and it is
the one input that turns a single-material section into a layered one:
everything between it and the maximum depth is now the foundation.

### 3. The `circles` worksheet

Two starting circles, sharing a center above the middle of the face at twice the
slope height — x = 20, y = 40 — and differing only in how deep they reach:

1. `circles!B3` **Xo** = `20`, `circles!C3` **Yo** = `40`,
   `circles!D3` **Option** = `Depth`, `circles!E3` **Depth** = `0`.
2. `circles!B4` **Xo** = `20`, `circles!C4` **Yo** = `40`,
   `circles!D4` **Option** = `Depth`, `circles!E4` **Depth** = `-10`.

![The finished circles worksheet](images/lem03_sheet_circles.png)

**Depth is the elevation of the circle's lowest point**, so `0` is a circle
tangent to the top of the foundation and `-10` one tangent to the rock. One per
layer base is the rule on layered ground: each circle asks about the mechanism
that rides that layer's bottom, and the search reports whichever turns out to be
worse.

Save the file and continue at [Running the analysis](#running-the-analysis).

---

## C — Building it in Studio {#c-building-it-in-studio}

Start with **File → New** and work down the **Inputs** tree.

### 1. Materials

Click **Materials**. The editor opens on **Table view**, which mirrors the `mat`
worksheet one material per row — the right shape here, because the two soils are
read against each other and the row order is what fixes the Mat IDs.

Press **Add row** twice and fill them: `embankment`, γ `130`, option `mc`,
c `400`, f `0`, u `none`; then `foundation`, γ `135`, option `mc`, c `800`,
f `0`, u `none`.

![The materials editor on this problem's two materials](images/lem03_studio_materials.png)

The swatch at the left of each row is the color that material's zone takes on
every plot. Click **OK**.

### 2. Profile lines

Click **Profile lines**, and set **Max depth (bottom boundary elevation):** to
`-10` — an elevation, which puts the rigid base 10 ft below the toe.

1. Press **Add line**, set **Material:** to `1: embankment`, and **Add row**
   three times for `0, 0` / `40, 20` / `90, 20`.
2. Press **Add line** again, set **Material:** to `2: foundation`, and **Add
   row** twice for `-30, 0` and `90, 0`.

A profile line is the *top* of a material layer, so the second line is the
contact between the two soils, and everything below it down to the maximum depth
is foundation. The preview redraws as you type: two lines in the two materials'
colors, the second running the full width of the section, and the hatched
**Max depth** line beneath both. Click **OK**.

### 3. Starting circles

Click **Circles**, then press **Generate starting circles…**. The generator
reads the geometry and proposes a set from it — one circle through the toe, one
at the base of each layer, all centered above the middle of the face at twice the
slope height.

![The circles editor on the generated circles](images/lem03_studio_circles.png)

On this section that is **three circles**, and the summary line under the button
says so: *"Generated 3 circles: 3 on the left-facing face (toe at x = 0, height
20)."* Read the **Depth** column, which is the elevation each one reaches:

| # | Depth | What it is |
|---|---:|---|
| 1 | −4.72136 | the toe circle — from this center, the arc through the toe bottoms out 4.7 ft inside the foundation |
| 2 | −10 | tangent to the rigid rock, the base of the foundation |
| 3 | 0 | tangent to the top of the foundation, the base of the embankment |

**Nothing was dropped here.** The toe circle is the candidate that usually goes —
on [LEM-1](lem01_simple_embankment.md)'s embankment, whose rigid base is at the
toe, the arc through the toe bottoms out below the model and the same summary
line reports *"1 candidate dropped for not daylighting inside the model"*. This
section has 10 ft of foundation under the toe and 30 ft of ground beyond it, so
all three candidates survive. Click **OK**.

Continue below.

## Running the analysis

However you built it, you now hold the same model:

![The finished model](images/lem03_inputs.png){width=1000}

The two profile lines are drawn in their materials' colors, and the two starting
circles share a center at (20, 40): the shallow one bottoming out on the contact,
the deep one on the rock.

Click **Run LEM…** and choose **Method** = `Spencer` and **Analysis** =
`Auto search`, with the slice count left at 40. From Python:

```python
from xslope.fileio import load_slope_data
from xslope.search import circular_search
from xslope.plot import plot_solution

sd = load_slope_data("my_layered_slope.xlsx")
fs_cache, _, path, circles = circular_search(sd, "spencer", num_slices=40)
crit = fs_cache[0]
plot_solution(sd, crit["slices"], crit["failure_surface"], crit["solver_result"])
```

---

## Exploring the results

### What the starting circles say before anything is searched

Each starting circle is a complete failure surface in its own right, and solving
the three of them — the two the file carries, plus the generator's toe circle —
is the fastest read on which mechanism this section has:

| Circle | Depth | Factor of safety | Weight of the sliding mass (lb/ft) |
|---|---:|---:|---:|
| Base of the embankment | 0 | **1.247** | 61,224 |
| Through the toe | −4.72 | 1.646 | 100,872 |
| Base of the foundation | −10 | 1.656 | 157,486 |

The deep circles carry up to two and a half times the weight of the shallow one
and are still a third safer, because everything they add is in the foundation:
twice the cohesion of the fill along every extra foot of base. On this section
the shallow mechanism controls before the search starts.

### The search result

![The circular search](images/lem03_search.png){width=1000}

The grey circles show both families being tried — the shallow set clustered on
the contact, the deep set running down to the rock. The green arrows trace the
center's walk from the two circles' shared start at (20, 40) to (18.5, 43.75):
both seeds are refined together, and both end on the same surface.

![Spencer's critical surface](images/lem03_solution.png){width=1000}

**FS = 1.244**, on a circle centered at (18.5, 43.75) with a radius of 43.75. Its
lowest point is at elevation 0.000: the critical surface is tangent to the
contact, running along the top of the foundation without entering it. Every one
of the 40 slice bases carries c = 400 psf, the fill's strength.

The solution carries two admissibility warnings — interslice tension (a most
tensile −1436 lb/ft against a largest compression of 7382) and a line of thrust
outside the slice on 15% of the boundaries. This is the crest tension of a φ = 0
slope that [LEM-1](lem01_simple_embankment.md) diagnoses and fixes with a tension
crack; adding one here at z<sub>c</sub> = 2c/γ = 6.15 ft gives a clean solution at
**FS = 1.175**, on a circle still tangent to the contact. It changes nothing
about which surface is critical, and unlike LEM-1's uncracked embankment it never
prevented one from being solved: Spencer answered on all 529 trial surfaces the
search tried.

### The audit at depth

The tangent depth is the dimension the layering acts on, so it is the one worth
auditing. Hold the critical circle's center and vary only how deep it reaches:

| Tangent depth (elevation) | +2 | +1 | +0.5 | **0** | −0.25 | −0.5 | −1 | −2 | −10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Factor of safety (Spencer) | 1.389 | 1.311 | 1.276 | **1.244** | 1.426 | 1.487 | 1.556 | 1.619 | 1.670 |

The minimum sits exactly on the contact, and the curve is not symmetric about it:
approaching from above the factor of safety falls smoothly, and a quarter of a
foot past it the answer jumps 15%. **A circle is flat near its lowest point**, so
dipping 3 in. below the contact does not put 3 in. of base into the foundation —
it puts 9.4 ft of base there, 16% of the surface, at twice the cohesion. The
strong layer acts as a floor, and the search settles on the deepest surface that
stays off it.

That single basin is why this model is forgiving about where the search starts.
Seeded with either of the file's circles alone, with the generator's toe circle
alone, or with any of seven deliberately bad guesses — a shallow circle on the
face, a small one in the crest, centers pushed out to x = 70 and x = −10, a 130
ft arc centered 120 ft up — Spencer's search
returns 1.244 on the same surface every time. The depth refinement is free to
walk across the contact, and from every one of those starts the downhill
direction points to the same place. **The per-layer circles are not what finds
this answer; they are what makes it checkable**, and the table above is the
check.

A section whose depth profile has two basins is not so forgiving, and layering is
the usual way to get one — [Sample Problem 13](../lem/samples.md#13-multiple-local-minima)
is a cohesionless embankment on soft clay where a free search collapses onto a
shallow sliver on the face and the deep foundation mechanism has to be seeded to
be found at all. Which is why the rule is per layer rather than per model: the
circle you can argue was unnecessary costs one row on a worksheet, and the one
you left out costs the answer.

### When the deep circle is the answer

Which layer is weak is a property of the soils, not of the geometry, and it is
the whole question. Take the same section and give the foundation c = 300 psf —
below the fill's 400 — with everything else unchanged:

![Spencer on the same section with a weak foundation](images/lem03_solution_weak.png){width=1000}

**FS = 0.792**, and it is a different mechanism entirely. The critical circle is
now tangent to the rock at elevation −10; it exits 9 ft beyond the toe, on the
flat ground the section provides for it, and reaches back to x = 64 in the crest.
It weighs 152,758 lb/ft against the 61,393 of the answer above, and its base runs
through both soils. Run the same depth audit and the curve has turned over: at
the contact 1.250, a foot below it 1.105, five feet below 0.875, and at the rock
0.792.

Two things follow. **A layered model has one candidate mechanism per layer**, and
a set of starting circles that names them all is how the model states that — the
circle that looked redundant here is the whole answer under a different soil.
And **a search cannot be audited by its own number**: 1.244 and 0.792 come from
the same geometry, the same three build paths and the same run settings, and the
only thing that tells them apart is reading how deep the reported surface went
and what it ran through.

### What the other methods say

The search finds each method its own critical circle:

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|---:|---:|---:|---:|---:|---:|---:|
| 1.244 | 1.244 | 1.313 | 1.326 | 1.285 | 1.244 | 1.244 |

The four methods that satisfy moment equilibrium — OMS, Bishop, Spencer and
Morgenstern-Price — land on the same circle and the same 1.244. With φ = 0 they
cannot disagree about a circle, and here nothing stopped any of them solving the
critical one. The three force-equilibrium procedures each settle on a flatter,
larger-radius circle of their own (Janbu R = 45.6, Lowe 47.7, Corps 52.9, against
Spencer's 43.75) and come out 3–7% higher. Spencer satisfies both force and
moment equilibrium and is the one to report.

---

## Conclusion

This tutorial demonstrated:

- A section with two material zones: one profile line per material, placed at the
  boundary where the material changes and nowhere else.
- Starting circles at the base of each layer, sharing a center above the middle
  of the face — and the generator that derives that set, plus a toe circle, from
  the geometry.
- Reading a search result at depth: the critical surface is tangent to the
  contact, its base entirely in the weaker fill, and the factor of safety jumps
  15% a quarter of a foot below it.
- Why a strong layer under a weak one puts a floor under the answer, and why
  softening it to 300 psf moves the critical surface to the base of the
  foundation and the factor of safety from 1.244 to 0.792.
- Moment-equilibrium methods agreeing exactly on a φ = 0 circle, with the
  force-equilibrium procedures 3–7% higher on circles of their own.

**Where to go next:** the [tutorials index](index.md) lists the rest of the
series. The sample problems carry the layering further —
[three layers with a piezometric line](../lem/samples.md#5-slope-with-multiple-materials-and-piezometric-line)
through them, a section whose
[layers are polygons](../lem/samples.md#11-polygon-input-with-a-sloping-bottom)
rather than profile lines because its base dips, and
[the two-basin slope](../lem/samples.md#13-multiple-local-minima) where the
starting circles decide the answer.
