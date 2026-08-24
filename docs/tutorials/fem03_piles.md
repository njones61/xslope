---
title: "Tutorial FEM-3 — Piles: LEM vs FEM"
description: "Two rows of drilled shafts and a sheet pile wall, each put through both of XSLOPE's engines — why a plane-strain finite element model turns a discrete pile row into a continuous wall, what pile spacing does to each engine's answer and what it cannot do to the finite element one, and the moment, shear and deflection a continuous wall carries, which no limit equilibrium analysis can report."
---

# Tutorial FEM-3 — Piles: LEM vs FEM

This tutorial shows which of XSLOPE's two engines models which kind of pile,
and why that is a question about the member rather than a matter of preference.
Both engines read the same pile rows out of the same file, and they build
different things from them.

A **limit equilibrium** method adds one force where a trial surface crosses a
pile. With the **Ito & Matsui (1975)** method that force is not entered: it is
computed from the pile diameter `D` and the center-to-center spacing `S`, from a
plasticity solution for soil squeezing *between* adjacent piles, and recomputed
for every trial surface. Spacing enters that calculation as the mechanism it
governs.

A **finite element** method meshes the pile into a chain of beam elements
sharing nodes with the soil around it, and divides the beam's axial stiffness
*EA* and its bending stiffness *EI* by `S` so that both are stated per unit
width of section. Nothing is prescribed. The beam carries whatever the deforming
soil pushes onto it, over its whole length rather than at one point, and the run
reports the moment, shear, deflection and soil reaction down the member.

What separates the two is that a two-dimensional analysis is **plane strain**:
every member in it is continuous out of plane. There is no gap between piles for
soil to move through, so a discrete row put into the finite element engine is a
wall. The rule that follows is measured on two models below. A discrete row at
spacing belongs in the limit equilibrium engine with Ito & Matsui, which models
the flow between piles directly. A member that is genuinely continuous — a sheet
pile, diaphragm or secant wall — or a row close enough to act as one belongs in
the finite element engine, which is also the only route that reports what the
member carries.

Work through [LEM-12](lem12_piles.md) before this page: it builds the pile model
used in the first half and covers the Ito & Matsui force, the structural
capacities that cap it, and what spacing does to a limit equilibrium answer.
[FEM-1](fem01_strength_reduction.md) covers strength reduction, meshing for a
stability run and the controls on the Run FEM dialog. Neither is repeated here.
[LEM-9](lem09_tieback_wall.md) is useful background for the second half, where a
wall appears: it enters its soldier pile the other way, with the force stated
directly per foot of wall.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Finite element (strength reduction)</p></div>
<div class="tgt-tile"><span class="tg-label">Build &amp; explore</span><p>~75 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn which engine models which kind of pile: why a plane-strain
model turns a discrete row into a wall, what spacing is worth to each engine,
and how to model a continuous wall and read the moment, shear and deflection it
carries.
</div>
<p><span class="tg-pill">piles</span><span class="tg-pill">Ito &amp; Matsui</span><span class="tg-pill">pile spacing</span><span class="tg-pill">beam elements</span><span class="tg-pill">smeared stiffness</span><span class="tg-pill">plane strain</span><span class="tg-pill">sheet pile wall</span><span class="tg-pill">continuous member</span><span class="tg-pill">strength reduction</span><span class="tg-pill">quadratic triangles</span><span class="tg-pill">steady seepage</span><span class="tg-pill">shear strain</span><span class="tg-pill">1D details</span><span class="tg-pill">moment and shear profiles</span></p>
<div class="tgm-model" markdown>
**Pile row model** — [xslope_piles.xlsx](../lem/files/xslope_piles.xlsx), the
completed model from [LEM-12](lem12_piles.md), which is also
[FEM sample problem 2](../fem/samples.md). It already carries both engines'
inputs, so nothing is entered on the first half of this page

**Wall starter file** — [xslope_pile_wall_start.xlsx](files/xslope_pile_wall_start.xlsx),
the sheet pile wall model with the wall itself left out; this is the file the
second half starts from

**Wall completed model** — [xslope_pile_wall.xlsx](files/xslope_pile_wall.xlsx),
the same model with the wall row entered; open it to skip to
[the wall's run](#re-meshing-and-running-again). Neither wall file carries a mesh
or a seepage solution, so both are built on the page
</div>
</div>

---

## The discrete pile row

![A 20 ft clay slope stabilized by two rows of 2 ft drilled shafts at 6 ft centers](images/fem03_piles_problem_sketch.png){width=1000}

The slope is a single medium-stiff clay — γ = 120 pcf, c = 200 psf, φ = 20° —
standing 20 ft at 1:1 over a rigid base 10 ft below the toe, with no water in
it. Two rows of 2 ft drilled shafts at 6 ft centers run down through the face to
the base, the lower at 5 ft from the toe and the upper at 10 ft. Each shaft is a
reinforced concrete section with 46,000 lb of shear capacity and 60,000 lb·ft of
moment capacity.

The sketch also carries what the finite element run reads and the limit
equilibrium run does not: the clay's elastic properties, E = 2.0 × 10<sup>6</sup> psf
and ν = 0.3, and the shafts' own modulus, E = 5.184 × 10<sup>8</sup> psf. Unlike
the starter files in [FEM-1](fem01_strength_reduction.md) and
[FEM-2](fem02_reinforcement.md), this model already has all of them, so both
engines run on it as it stands.

### Opening the model and reading the pile rows

Download [xslope_piles.xlsx](../lem/files/xslope_piles.xlsx) and open it with
**File → Open…**. The mode strip opens on **LEM**, which is where the first run
happens. Units are `imperial`, so lengths read in feet and stresses and
stiffnesses in psf.

![The pile model as it opens: two rows through the face and the starting circle](images/fem03_inputs_piles.png){width=1000}

The Inputs plot draws the section, the two pile rows as green bars running from
the face down to the hatched maximum-depth line at elevation −10, and the dashed
red starting circle the file carries for the search.

Open **Piles** in the **Inputs** dock and press **Table view**. Leave both
**Show parameters for:** toggles ticked, because this model is the case where
both bands matter at once:

![The two pile rows with both usage bands shown](images/fem03_studio_piles_table.png)

The columns are colored by which engine reads them: the red ones are read only by
the limit equilibrium engine, the blue ones only by the finite element engine.
`D` and `S` are black because both engines read them, and on this file that
matters. Every stiffness the finite element model uses comes from those two
cells: `I` and `Area` are blank, so the engine derives the solid circular section
from the diameter, I = πD⁴/64 and Area = πD²/4, and then divides *EA* and *EI* by
the spacing. Hovering `S` says what each engine does with it:

> Center-to-center pile spacing. In the LEM, spacing is physics: it sets the
> arching between piles (Ito & Matsui) and makes Vcap/Mcap per-pile. In the FEM,
> S converts per-pile properties to per-unit-width (EA/S, EI/S) — the correct
> plane-strain stiffness of the smeared row; what it cannot model is the soil
> arching between piles.

`H` is blank on both rows, which is what tells the limit equilibrium engine to
compute the force by Ito & Matsui rather than use a stated one, and `Appl` is
`Active`, meaning that force enters the equilibrium equations as it stands.
[LEM-12](lem12_piles.md#the-pile-rows) reads the whole row column by column.
Click **Cancel**.

### The limit equilibrium answer

The limit equilibrium run comes first because it is the reference the finite
element run is read against, and because the file is already complete for it.

Click **Run → Run LEM…**, set **Method** to `Spencer` and **Analysis** to
`Auto search`, and leave the slice count at 40. Spencer is the method to compare
against because it satisfies both force and moment equilibrium, which makes it
the closest limit equilibrium statement of what a strength reduction run solves.
Click **Run**.

![Spencer's critical circle, with the two pile crossings marked](images/fem03_lem_solution_piles.png){width=1000}

**FS = 1.842**, on a deep circle that leaves the ground 14 ft beyond the toe and
exits on the crest at x = 35. The two red dots are where it crosses the pile
rows. The rows deliver **4,367.6 lb/ft** to the slices there, both of them
limited by the shafts' moment capacity rather than by what the soil could
supply; [LEM-12](lem12_piles.md#how-the-force-is-computed) works through that
calculation.

### The finite element answer

The same file, the same rows, the other engine. Switch the mode strip to
**FEM**. **Run → Run FEM…** stays disabled until a mesh exists, so build one
first with **Run → Build Mesh…**

Set **Element type** to **Quadratic triangles (tri6)** — a stability mesh has to
be quadratic, and [FEM-1](fem01_strength_reduction.md#building-the-mesh)
measures what linear elements cost. Untick **Auto-size from geometry** and set
**Target element size** to `2`, which is what the
[FEM sample problem](../fem/samples.md) for this model uses. Leave the rest of
the dialog alone and click **Build**.

The mesh comes out at **3,180 nodes and 1,521 triangles**, with the two pile
rows carried in as constraints and discretized into **18 beam elements** — 8 on
the lower row and 10 on the upper one.

Click **Run → Run FEM…**

![Run FEM on the meshed pile model](images/fem03_studio_run_fem_piles.png)

**Model checks — 1 warning**, and **Run** is enabled. The warning is the blank
tensile cutoff `t_cut` on the clay, which lets it carry tension up to its
Mohr-Coulomb apex, c/tan(φ) = 549.5 psf;
[FEM-1](fem01_strength_reduction.md#running-the-strength-reduction) covers what
that means and when to enter a cap. The note below it says that the two pile
rows carry a diameter with no `Area` or `I`, so the engine is deriving the solid
circular section — the derivation described above, reported rather than assumed.

Set **F min (SSRM)** to `1.00` and **F max (SSRM)** to `1.60`, which brackets
this model's answer, and leave everything else as it opens: tolerance 0.0100,
**Max iterations per trial** 12,000, **Iteration ceiling** 50,000, **Rollers**
on the sides, and **Non-convergence** as the failure criterion. Click **Run**.
The run takes about 40 seconds.

**FS = 1.361**, from a final bracket of [1.3563, 1.3656] in eight bisection
steps. Spencer's method gave 1.842 on the same file.

![The mechanism at failure, with the pile rows colored by shear force](images/fem03_fem_shear_piles.png){width=1000}

The contours are viscoplastic shear strain — the shearing left after the elastic
response is subtracted. The band runs up from the flat ground beyond the toe,
past the lower row, and concentrates on the uphill side of the upper row before
breaking out of the face just above it at about x = 12, y = 8. Nowhere does it
pass *through* a row, because in plane strain there is nothing to pass through:
each row is a continuous obstruction over the full length of the slope, and the
mechanism has to climb over it. That is the statement at the top of this page,
seen in a field rather than argued.

### What the row is worth in each engine

Neither number above says what the piles contributed, and the two engines have
to be compared on that rather than on the piled answer alone. Measuring it takes
the same slope without them, run both ways.

Open **Piles**, and on **Table view** select both rows and click **Remove
selected**, then **OK**. Run Spencer again exactly as before — **Run → Run LEM…**, `Spencer`,
`Auto search`, 40 slices:

**FS = 1.149**, on a shallower circle from the toe.
[LEM-12](lem12_piles.md#what-the-two-rows-are-worth) shows that surface.

Now the same model in the finite element engine. The pile rows were mesh
constraints, so the mesh has to be rebuilt without them: switch to **FEM**,
**Run → Build Mesh…**, the same tri6 at 2 ft, **Build**. It comes back at
**3,154 nodes and 1,509 triangles** and no beam elements. **Run → Run FEM…**,
the same 1.00 to 1.60 bracket, **Run**:

**FS = 1.164**.

| | Without the rows | With them | Credit for the rows |
|---|:---:|:---:|:---:|
| **LEM** — Spencer, Ito & Matsui | 1.149 | 1.842 | ×1.60 |
| **FEM** — strength reduction | 1.164 | 1.361 | ×1.17 |

Read the first column first. On the bare slope the two engines land 1.3% apart,
so nothing structural separates them on this section, and whatever the second
column adds is the pile row and only the pile row. There they credit it by a
factor of 1.60 and 1.17 — a disagreement on the quantity being designed, not a
rounding.

Reopen [xslope_piles.xlsx](../lem/files/xslope_piles.xlsx) with **File → Open…**
to put the rows back before the next section.

### How the pile force is applied

Part of that gap is a convention rather than an idealization, and it is worth
separating out before the rest is attributed to plane strain. The `Appl` column
on a pile row has two settings. `Active` hands the computed force to the
equilibrium equations as it stands, which is what this file does. `Passive`
treats it as a resistance and divides it by the factor of safety, so the support
carries the same margin as the soil strength does. Slide2 applies a pile force
the passive way, and XSLOPE's own [VP106](../verification/rocscience.md#vp106)
files are built that way to match it.

Set `Appl` to `Passive` on both rows and search again. Spencer returns
**1.574** instead of 1.842. The choice is worth 0.268 of factor of safety on
this model, against 0.481 between the limit equilibrium and finite element
answers — over half of the difference between the engines. It also moves the
critical surface: under `Passive` the deep circle loses to a shallow one that
skims the toe. Set the column back to `Active` before continuing.

### What spacing does to each engine

A designer adjusts the spacing once the diameter is set, and the two engines part
company over it more sharply than over anything else. The sweep below runs both
of them across the applicable band — 3, 6 and 12 ft, S/D from 1.5 to 6 — with
nothing changing but the `S` cell on the two rows. Each limit equilibrium point
is its own Spencer search, because spacing changes which surface governs.

![Factor of safety against pile spacing, both engines on one axis](images/fem03_spacing_sweep.png){width=800}

The limit equilibrium answer falls from 2.193 to 1.409 over that range, 36%. The
strength reduction answer is 1.361 at all three spacings, to six decimals, from
the same bisection interval and with the same verdict at every trial.

Run one point of it to see both. Open **Piles**, set `S` to `12` on both rows
with `H` still blank, and **OK**. Search with Spencer again: **FS = 1.409**, and
not on the deep circle the 6 ft spacing found — the deepest point of the
critical surface rises from elevation −5.1 to −0.8. Widening the gap did not
only lower the answer, it changed which surface governs, which
[LEM-12](lem12_piles.md#what-the-spacing-is-worth) covers in full. The pile
lines have not moved, so the mesh still fits: switch to **FEM** and
**Run → Run FEM…** with the same bracket. **FS = 1.361**, the same answer as at
6 ft.

That flat line is worth explaining carefully, because the obvious reading of it
is wrong. Dividing *EA* and *EI* by `S` is not bookkeeping the finite element
engine does to dispose of the spacing column — it is the correct plane-strain
stiffness of a row smeared into a wall, and it responds to `S` exactly as it
should. Over this sweep the smeared stiffness falls by a factor of four, and the
run's internals move with it:

| S (ft) | *EI*/S (lb·ft²) | Peak moment per unit width (lb·ft/ft) | Peak moment per pile (lb·ft) | Fraction of M<sub>cap</sub> | Elements yielded in bending |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 3 | 1.357 × 10<sup>8</sup> | 6,518 | 19,553 | 33% | 0 |
| 6 | 6.786 × 10<sup>7</sup> | 6,319 | 37,912 | 63% | 0 |
| 12 | 3.393 × 10<sup>7</sup> | 5,000 | 60,000 | 100% | 8 |

Read across the two moment columns. The moment *per unit width of slope* barely
moves, which is the quantity the smear holds nearly constant. The moment *per
shaft* is that number times the spacing, and it triples: at 12 ft each shaft
reaches its full 60,000 lb·ft of capacity and eight beam elements yield in
bending. The model is responding to spacing in every quantity a designer would
check. What does not move is the factor of safety.

Two things explain that. A strength reduction analysis finds the factor by which
the soil's strength has to be divided before equilibrium becomes impossible, so
the answer is decided by strength, and stiffness reaches it only through the
ratio of the member's stiffness to the soil's. Across this sweep that ratio
stays high — a steel-stiff shaft against a soft clay — so the wall stays
effectively rigid and the mechanism has to climb over it at every spacing. Spread
the row far enough and the smeared wall would go flexible, the mechanism would
begin to bend it, and the factor of safety would fall back toward the unpiled
1.164. What the two-dimensional model cannot represent at any spacing is the
mechanism that spacing actually governs — soil arching onto the shafts and
squeezing through the gaps between them, which is exactly what Ito & Matsui
computes. Plane strain has no gaps.

Set `S` back to `6` on both rows before the next section, or reopen the file.

### The one problem with a three-dimensional answer

Neither number in the comparison table is a three-dimensional answer, and the
direction of the plane-strain error is only known where a three-dimensional
reference exists. One benchmark in XSLOPE's corpus has one. Cai & Ugai (2000)
analyzed a pile-stabilized slope with a strength reduction finite element model
that meshes the individual piles, the soil between them, and the slip interface
on each pile's surface. XSLOPE runs the same slope through both of its engines
in [VP106](../verification/rocscience.md#vp106) and
[the VP106 finite-element diagnostic](../verification/rocscience.md#vp106-fem):

| Case | XSLOPE SSRM (2D beam) | Cai & Ugai 3D FE |
|---|:---:|:---:|
| No pile | 1.136 | 1.14 (−0.4%) |
| Pile at D<sub>1</sub>/D = 3, free head | 1.453 | 1.36 (+6.8%) |
| Pile, head rotation restrained | 1.570 | 1.45 (+8.3%) |

The unpiled row agrees to 0.4%, which is what makes the other two readable. With
the row in place the plane-strain model credits it with multiplying the
unreinforced factor of safety by 1.279, where the three-dimensional model
credits 1.193. On the same slope a Bishop search with the Ito & Matsui force
reads 1.451, a credit of 1.269. Both two-dimensional routes sit above the
three-dimensional answer by a comparable amount, and neither recovers it. What
this benchmark settles is the direction of the plane-strain error, not a ranking
of the two routes against each other.

So the reason to take a discrete row's factor of safety from the limit
equilibrium engine is not that its number is larger or smaller. It is that Ito &
Matsui is a theory *of* the three-dimensional mechanism, while the beam model is
a two-dimensional substitute for it. The finite element run on a discrete row
still tells you what the shafts are carrying, which is why the moment table
above is worth reading; it is the factor of safety that carries the
idealization.

---

## The continuous wall

![A 20 m slope over a weak clay band, held by a sheet pile wall driven from the bench](images/fem03_wall_problem_sketch.png){width=1000}

The second half is the case the finite element engine is right for, and it is a
different model with a different unit system. Everything below is **SI** —
meters, kN/m³ and kPa — because this model is a transcription of GeoStudio's
SIGMA/W *slope stabilization with piles* example and converting it would cut it
loose from the published answer it is checked against.

The section is 64 m wide and 24 m tall, from a rigid base at elevation −4 to a
crest at elevation 20. The face runs from the crest at 1.6:1 down to a bench at
elevation 10, then at 2:1 down to level ground at elevation 8. Sandy clay fills
everything above elevation 6. Under it a **1 m band of weak clay** runs the full
width between elevations 5 and 6, with no cohesion and φ′ = 10°, and under that
a stiff overconsolidated soil, modeled as linear elastic, reaches the base. The
weak band is what makes the slope marginal. A **sheet pile wall** is driven from
the middle of the bench at x = 38, from elevation 10 down to elevation 1 — 9 m,
through the band and 4 m into the stiff soil under it.

Pore pressures come from a **steady-state seepage solution** rather than from a
water table drawn on the section: total head is held at 14 m on the upstream
face and 8 m downstream, and the seepage run computes the field between them.
[SEEP-2](seep02_johnson_dam.md) covers how a steady seepage solve is set up and
read; this page runs it and uses the field.

### Opening the starter file

Download [xslope_pile_wall_start.xlsx](files/xslope_pile_wall_start.xlsx) and
open it with **File → Open…**. It is the whole model except the wall: the
geometry, the three materials with their elastic properties, the seepage
boundary conditions, and the element type and target size the mesh needs.

![The wall starter: three profile lines and no wall](images/fem03_inputs_wall.png){width=1000}

The Inputs plot draws the three material zones as profile lines — the sandy clay
down to elevation 6, the weak clay band between 6 and 5, and the overconsolidated
soil below it — over the hatched maximum-depth line at elevation −4. The band is
thin enough on this scale that its two lines nearly touch, which is why the file
carries a local element size override of 0.3 m on it: a controlling band needs
several element rows across it whatever the rest of the mesh is doing.

### Meshing and solving the seepage

Switch the mode strip to **Seepage** and click **Run → Build Mesh…** The file
already declares **Quadratic triangles (tri6)** with **Auto-size from geometry**
off at `1.5`, so both are set when the dialog opens. Click **Build**.

The mesh comes out at **13,127 nodes and 6,484 triangles**, refined through the
weak band by the size override the file carries.

Click **Run → Run Seep…** and **Run**, with nothing on the dialog changed. It
converges in about half a second, at a flow rate of
1.9432 × 10<sup>−5</sup> m³/s per meter of section. The solved pore pressures go
into the model, where the strength reduction run reads them: all three materials
have their `u` set to `seep`.

The order of those two steps matters for the rest of the page. Building a mesh
discards any seepage solution computed on the old one, so a re-mesh is always
followed by a re-solve.

### The slope without the wall

Switch the mode strip to **FEM** and click **Run → Run FEM…** The **Model
checks** panel reads *No problems found for this run*. Set **F min (SSRM)** to
`0.95` and **F max (SSRM)** to `1.25`, leave everything else as it opens, and
click **Run**.

This run takes about **20 minutes**. It is a 6,484-element mesh and the trials
near the critical factor iterate tens of thousands of times each before they
settle, so it can be left going.

**FS = 1.020**, from a final bracket of [1.0156, 1.0250] in seven trials. The
slope stands, barely. SIGMA/W reads about 1.025 on the same model, and two other
analyses in the same project put it at 1.033 and 1.035 —
[the SIGMA/W wall benchmark](../verification/geostudio.md#sigmaw-wall) states
how each of those published figures is read.

### Adding the wall

The wall goes in as one row of the piles table, and which cells it fills is the
content of this step. Open **Piles** in the **Inputs** dock, press **List
view**, and click **Add**.

![The wall row: S = 1, no diameter, no capacities, section constants entered directly](images/fem03_studio_wall_row.png)

Enter the geometry — `x1` 38, `y1` 10, `x2` 38, `y2` 1 — and then the three
section values, all of them per meter of wall:

| E | I | Area |
|---:|---:|---:|
| 200000000 | 0.0005 | 0.02 |

That is 2.0 × 10<sup>8</sup> kPa of modulus, 5.0 × 10<sup>−4</sup> m⁴/m of
moment of inertia and 0.02 m²/m of area, giving *EA* = 4.0 × 10<sup>6</sup> kN/m
and *EI* = 1.0 × 10<sup>5</sup> kN·m²/m.

Set **S** to `1`, and leave `H`, `D`, `Vcap` and `Mcap` empty.

`S` = 1 is what makes this member continuous. A sheet pile wall has no
out-of-plane gap, so its section constants already are per meter of wall, and a
spacing of 1 passes them into the beam elements unchanged instead of smearing
one member's stiffness over a spacing. The blank cells follow from the same
fact. `D` is blank because there is no circular shaft to derive a section from —
`I` and `Area` are given directly, which is the other way the engine accepts a
section. `Vcap` and `Mcap` are blank because this run is asked to report what the
wall carries rather than to check it against a declared capacity, and, as the
last section shows, the detail panel draws a capacity line only where one is
entered. `H` is blank because a finite element run never uses a stated pile
force.

**Appl** reads `passive` on this row, and **Fixity** `free`. `Appl` is a limit
equilibrium input and does nothing in a strength reduction run. `Fixity` `free`
leaves both ends of the wall able to rotate, which is what a driven sheet pile
with no capping beam does, and the moment profile below is read against it.

Click **OK**.

### Re-meshing and running again

The wall line is a mesh constraint — the beam elements have to lie on element
edges and share their nodes with the soil — so the mesh has to be built again
with it in place, and the seepage solution has to be recomputed on the new mesh.

Switch to **Seepage**, **Run → Build Mesh…** and **Build**. The mesh comes out
at **13,223 nodes and 6,532 triangles**, 48 more than before, with the wall
discretized into **18 beam elements** over its 9 m. Then **Run → Run Seep…** and
**Run**: it converges at 1.9429 × 10<sup>−5</sup> m³/s per meter, within 0.02%
of the flow rate without the wall. A pile row is a structural member in XSLOPE
and not a flow barrier, so the seepage solve sees only the slightly different
mesh; the wall changes the stability run, not the water.

Switch to **FEM** and click **Run → Run FEM…**

![Run FEM on the meshed wall model](images/fem03_studio_run_fem_wall.png)

Set **F min (SSRM)** to `1.15` and **F max (SSRM)** to `1.65` and click **Run**.
This run takes about **30 minutes**.

**FS = 1.420**, from a final bracket of [1.4156, 1.4234] in eight trials. The
wall takes the slope from 1.020 to 1.420, a credit of ×1.39, against SIGMA/W's
own reading of about 1.4 for the same model.

One note on that number before the figures. The
[verification page](../verification/geostudio.md#sigmaw-wall) records 1.451 for
this model rather than 1.420. The two runs differ in the **Failure criterion**
setting: the Run FEM dialog opens on **Non-convergence**, the plain reading that
a trial which cannot reach equilibrium has failed, and the verification run uses
the engine's hybrid criterion, which also weighs how the displacements are
behaving when a trial spends its iteration ceiling. On this model five of the
eight trials reach that ceiling without settling, so the two criteria adjudicate
them differently and the two searches close on different intervals. Both answers
sit inside the published "about 1.4". [FEM-1](fem01_strength_reduction.md)
compares the criteria the list offers. The two figures below are drawn from the
verification run, so their titles read 1.45.

![The mechanism with the wall in place, running along the weak clay band](images/fem03_wall_shear.png){width=1000}

The shear strain band is nothing like the pile model's. It runs almost flat
along the weak clay band between elevations 5 and 6, from under the crest out to
beyond the toe, because that is the only weak material in the section. The wall
crosses it at x = 38, colored by shear force against the second color bar, with
its darkest segment just below the band, where the shear peaks.

### What the wall carries

This run produces two things. The factor of safety is the first, and the second
is what the wall itself is carrying, which no limit equilibrium analysis can
report. Click **1D Details…** on the results toolbar and select the wall in the
list on the left.

![The 1D details panel on the wall](images/fem03_studio_wall_1d_details.png)

![Lateral displacement, shear, moment and soil reaction down the wall](images/fem03_wall_profiles.png){width=1000}

Four profiles, all plotted against depth below the pile head, which is at
elevation 10 — so a depth of 5.00 m on these axes is elevation 5.00, the base of
the weak clay band. The orange stretch shaded from 4.00 to 4.25 m — **Shear
band crossing** in the legend — is where the mechanism crosses the wall,
elevation 6.00 down to 5.75, inside the clay. The panel opens on
**Field state = At failure**, the state the run re-solves once past the factor
of safety so the collapse develops far enough to read; the selector at the
bottom switches it to the last converged state, where every quantity below is
about a quarter smaller.

The **bending moment** is zero at the head and zero at the toe, which is the
check that the profile is being read the right way round — both ends are free,
so neither can carry a moment — and it peaks at **1,118 kN·m/m at a depth of
5.00 m**, the base of the band.

The **shear** holds one sign above that elevation, peaking at 505 kN/m, and
reverses below it to a peak of −852 kN/m at 5.49 m. The **soil reaction** panel
shows the same reversal as a spike, 3,216 kN/m per meter of wall at the same
5.00 m: the soil pushes hard on the wall from upslope through the band, and the
stiff overconsolidated material below pushes back. Between them the wall is
being driven by the weak band and reacting against what it is embedded in, which
is what a sheet pile driven through a sliding layer is for. The **lateral
displacement** panel is the shape that produces: 0.118 m at the head, decaying
smoothly with depth and down to 1 mm at the toe, without ever changing
direction — the wall is bending like a cantilever fixed in the stiff soil, not
sliding with the block.

The published example's own moment and shear diagrams have the same form and the
same turning point at the base of the band, with peaks about half again as
large; the benchmark page
[states how far each sits from the published value and why](../verification/geostudio.md#sigmaw-wall).

Two things the panel does *not* draw are worth noticing.
There is no Ito & Matsui limiting-resistance envelope on the moment axis and no
capacity lines on the shear or moment panels. Both are correct here. The
envelope is a limit pressure for soil flowing between piles, and this member has
no gaps for soil to flow through; the capacity lines would need a `Vcap` and
`Mcap` the row does not declare, and drawing a capacity that was not entered
would be inventing one. The status beside the field selector reads *no capacity
declared* for that reason.

---

## Which engine for which member

| Member | Out of plane | Engine | What it gives |
|---|---|---|---|
| Sheet pile, diaphragm or secant wall | continuous | FEM, spacing `S` = 1 | factor of safety, plus moment, shear, deflection and soil reaction down the member |
| Contiguous or very closely spaced row | nearly continuous | FEM, with the smear stated | the same, with the gap unrepresented |
| Discrete row at spacing | discrete | LEM with Ito & Matsui | factor of safety per spacing, force per row, capacity checks |

The middle row is a judgment rather than a threshold, and the bottom row has a
band of its own: the Ito & Matsui theory applies for S/D between about 2 and 8,
and XSLOPE says so in the **Model checks** panel when a spacing falls outside it
— [LEM-12](lem12_piles.md#what-the-spacing-is-worth) runs the sweep to both
edges.

Two things follow for a discrete row, and both are stated in
[LEM vs. FEM pile modeling](../lem/piles.md#lem-vs-fem-pile-modeling). Take its
factor of safety from the limit equilibrium analysis, and read the finite
element counterpart as a stiffness-and-force study — what the shafts carry and
where they yield — rather than as a competing factor of safety. And do not try
to repair the plane-strain smear by adjusting the pile stiffness or by imposing
the Ito & Matsui limit pressure on the beam: that limit pressure is a theory of
the very mechanism the two-dimensional model does not contain, so applying it
there counts the same resistance twice.

---

## Conclusion

This tutorial covered:

- What a pile row is to each engine — one computed force at the crossing in the
  limit equilibrium method, a chain of beam elements carrying whatever the soil
  pushes onto them in the finite element method.
- Why a two-dimensional model turns a discrete row into a continuous wall, seen
  in a shear strain field where the mechanism climbs over the rows rather than
  passing through them.
- What spacing does to each engine: a 36% swing in the limit equilibrium answer
  and a change in which surface governs, against a strength reduction answer
  that does not move, even though the smeared stiffness falls fourfold and each
  shaft reaches its moment capacity.
- What the one benchmark with a published three-dimensional answer establishes,
  and what it does not.
- How a continuous wall is entered — spacing 1, section constants per meter, no
  diameter and no capacities — and the re-mesh and re-solve a new structural
  line forces on a model whose pore pressures come from seepage.
- The moment, shear, deflection and soil reaction down a wall, which a limit
  equilibrium analysis cannot produce, and why the detail panel draws no arching
  envelope for a member with no gaps.

**Where to go next:** [Piles and concrete piers in FEM](../fem/piles.md) carries
the beam formulation, the assembly and the applicability rule in full, and
[stabilizing piles in LEM](../lem/piles.md) carries the Ito & Matsui theory and
the same rule from the other side;
[the SIGMA/W wall benchmark](../verification/geostudio.md#sigmaw-wall) is the
second half of this page checked against its published source;
[LEM-12](lem12_piles.md) is the pile row on its own. COMBO-1, Seepage →
Stability, is the next tutorial on the roster and takes up as its subject the
seepage-then-stability sequence that appeared here as a step.
