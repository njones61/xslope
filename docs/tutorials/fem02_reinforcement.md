---
title: "Tutorial FEM-2 — Reinforcement: LEM vs FEM"
description: "Six layers of geogrid in a 24 ft sand fill, solved by Spencer's method and then by finite element strength reduction — the two limits every reinforcement model carries, the axial stiffness the finite element engine refuses to start without, the elastic-perfectly-plastic run against the peak-residual one, what the residual capacity is worth on this slope, and what changes when the bond is read from the overburden instead of a stated length."
---

# Tutorial FEM-2 — Reinforcement: LEM vs FEM

Here we look at how a finite element analysis models soil reinforcement, and how
its answer differs from the limit equilibrium answer for the same slope.
Both engines read the same reinforcement lines out of the same file, and they do
different things with them. A limit equilibrium method has no strain in it:
wherever a trial surface crosses a line it adds the tensile force that line can
develop at the crossing point, and the rest of the line does not enter the
calculation. A finite element method meshes each line into a row of bar elements
sharing nodes with the soil around them, and the force in every one of those
elements comes out of how far the soil beside it has moved.

The example is the geogrid-reinforced sand fill from
[LEM-8](lem08_reinforced_slope.md) — six layers of geogrid in a 24 ft fill under
a crest surcharge — and we solve it three times. Spencer's method gives the
reference answer. Then we mesh the same model and run it by **strength
reduction** twice: once with the layers holding their full strength after they
yield, which is what every mainstream finite element code assumes unless told
otherwise, and once with a residual capacity that takes effect after the layer
ruptures.

Work through [LEM-8](lem08_reinforced_slope.md) first if you have not: that is
where we built this model, and it covers reinforcement lines, the capacity
envelope and pullout lengths, all of which we lean on here. Strength reduction,
meshing for a stability run and the controls that decide whether a trial is
allowed to finish are covered in [FEM-1](fem01_strength_reduction.md). Neither is
repeated here. We start from a **starter file** that carries the whole
LEM-8 model and nothing else, so the inputs we add are the ones the finite
element engine needs and the limit equilibrium engine does not: the soils' two
elastic properties, and three columns on every reinforcement line.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Finite element (strength reduction)</p></div>
<div class="tgt-tile"><span class="tg-label">Build &amp; explore</span><p>~30 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn how reinforcement enters a finite element stability
analysis: the two limits a reinforcement model carries, the axial stiffness the
run needs, how to read what each layer is doing at the factor of safety, and how
much the post-peak assumption costs.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">reinforcement lines</span><span class="tg-pill">capacity envelope</span><span class="tg-pill">pullout length</span><span class="tg-pill">bar elements</span><span class="tg-pill">axial stiffness</span><span class="tg-pill">elastic-perfectly-plastic</span><span class="tg-pill">residual capacity</span><span class="tg-pill">strength reduction</span><span class="tg-pill">quadratic triangles</span><span class="tg-pill">thin zones</span><span class="tg-pill">iteration budget</span><span class="tg-pill">overburden pullout</span><span class="tg-pill">1D details</span><span class="tg-pill">shear strain</span></p>
<div class="tgm-model" markdown>
**Starter file** — [xslope_reinforced_slope_start.xlsx](files/xslope_reinforced_slope_start.xlsx),
the LEM-8 model with every finite element input blank: Young's modulus and
Poisson's ratio on both soils, and `Tres`, `E (psf)` and `Area` on all six
reinforcement lines; this is the file the page starts from

**Completed model** — [xslope_reinforced_slope.xlsx](files/xslope_reinforced_slope.xlsx),
the same model with the soils' elastic pair and the three reinforcement columns
filled in and the element type and target size declared; open it to skip to
[the runs](#the-iteration-budget). Neither file carries a mesh, so the meshing
step below is done on either one
</div>
</div>

---

## The problem

![A 24 ft reinforced sand fill: six geogrid layers developing their tension over 4 ft at each end, a cohesive face band, and a crest surcharge](images/fem02_problem_sketch.png){width=1000}

The fill stands **24 ft** high at **1.25:1** on a foundation that runs 10 ft
below the toe, with a **240 psf** surcharge over the 70 ft of crest behind it.
Two Mohr-Coulomb soils at the same unit weight: the fill itself is clean sand
with no cohesion at all, and the face is wrapped in a band of cohesive fill —
the `shell` — 2 ft measured horizontally, which on a 1.25:1 face is 1.25 ft
perpendicular to it; its 300 psf of cohesion keeps the search off the face. Six geogrid layers, each **20 ft** long and 4 ft apart vertically, start at
the face. Each develops **800 lb/ft** of tension over a pullout length of
**4 ft at both ends**, which is the published example's own value. It treats the
two ends as alike, and they are not — the buried end carries 4 to 16 ft of fill
and the face end almost none — so in
[pullout from the overburden](#pullout-from-the-overburden) below we run the same
model with the bond read from the depth of burial instead.

The sketch shows the soils' elastic properties, E = 1.0 × 10<sup>6</sup> psf and
ν = 0.3, the geogrid's axial stiffness, *EA* = 80,000 lb/ft, and its residual
capacity, 600 lb/ft, because the finished model carries them; none of the four is
in the starter file, and we enter all of them as the finite element run calls for
them. The geometry, the soils and the reinforcement are Example 5 from the
UTEXASED user manual, built from nothing in
[LEM-8](lem08_reinforced_slope.md).

---

## Opening the starter file

Download [xslope_reinforced_slope_start.xlsx](files/xslope_reinforced_slope_start.xlsx)
and open it with **File → Open…**. The toolbar's mode strip opens on **LEM**,
which is where our first run happens. Units are `imperial`, so lengths read in
feet, unit weights in pcf, and stresses, strengths and stiffnesses in psf.

![The starter file: two profile lines, six reinforcement lines, the crest surcharge and two starting circles](images/fem02_inputs.png){width=1000}

The Inputs plot shows the face band and the fill as two profile lines, the
hatched maximum-depth line at elevation −10, the crest surcharge as purple
arrows, and the six reinforcement lines as gray layers stepping up the face. The
two dashed red arcs are the starting circles the search begins from.

---

## The limit equilibrium answer

We run the limit equilibrium analysis first because it gives the reference the
two finite element runs are read against, and because the starter file is already
complete for it.

Click **Run → Run LEM…**, choose **Method** = `Spencer` and **Analysis** =
`Auto search`, and leave the slice count at 40. **Model checks** reads
*No problems found for this run*.
We compare against Spencer because it satisfies both force and moment
equilibrium, which makes it the closest limit equilibrium statement of what a
finite element run solves. Click **Run**.

![Spencer's critical circle through the reinforced block](images/fem02_lem_solution.png){width=1000}

**Spencer's method gives FS = 1.587**, on a circle centered at (−5.13, 46.98)
with a radius of 47.26 ft. The surface daylights at the toe, cuts up through the
reinforced block and comes out on the crest at **x = 36.2**, behind the last
geogrid. It crosses five of the six layers and takes 800 lb/ft from each of
them, for 4,000 lb/ft of tension in total.
That crossing-by-crossing is worked through in
[LEM-8](lem08_reinforced_slope.md#where-the-surface-meets-the-lines).

<!-- test: file=files/xslope_reinforced_slope_start.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.587, tolerance=0.005 -->

---

## What a reinforcement line is to each engine

A reinforcement element is an elastic layer with **two separate limits** on the
force it can carry. The **bond limit** is the force the soil can transfer into
the layer through friction along its embedded length; when the layer force reaches
it, the layer slips. The **rupture limit** is the strength of the reinforcement
itself; when the force reaches it, the material yields. Both limits are
**perfectly plastic**: once the force reaches the limit it stays there while
the layer keeps stretching, neither rising nor falling. Below either limit the layer
is **elastic** — its force grows in proportion to how far it is stretched — so
the whole behavior is called **elastic-perfectly-plastic**: force proportional
to stretch up to the limit, then constant at it. A material can also be given a
**residual**: a lower force it drops to after rupture, instead of holding at the
limit. That two-limit structure is the standard idealization the mainstream
finite element codes share; PLAXIS, RS2 and FLAC all model reinforcement this
way. XSLOPE handles rupture the same way: `Tmax` is the rupture limit. With
`Tres` left blank, a layer that reaches it holds there — the elastic-perfectly-
plastic geogrid those codes use by default. With a `Tres` entered, a layer that
reaches `Tmax` drops to that residual and carries no more than it from then on.
The bond limit can be stated either way. Those codes compute bond from the
stress on the layer-soil interface, so the length a layer needs to develop its full
capacity depends on how deep it is buried; XSLOPE offers that form through the
`Adhesion` and `Delta` columns, and reads the length directly from the input —
the pullout lengths `Lp1` and `Lp2` — when they are left blank, which is what
this model does. The distinction between the two limits matters when the
results are read: the 1D Details panel reports whether a line reached its bond
limit (**pullout**) or its rupture limit (**yielded**), and the two mean
different things for the design.

**Bond** is the first limit. A geogrid is held by friction against the soil it
is buried in, and near a free end there is not much soil to be held by, so the
force the line can carry there is limited by its embedment rather than by its
own strength. That is what the pullout lengths `Lp1` and `Lp2` describe: the available force
ramps from zero at an end up to the full 800 lb/ft over that length, and past it
the embedment can develop the whole capacity. Here the ramp is 4 ft at each end.
(Stating `Adhesion` and `Delta` instead makes that ramp follow the effective
overburden along the line rather than rise at a fixed rate; this model uses the
lengths, and in [the section below](#pullout-from-the-overburden) we run it the
other way.) The ramp
and the plateau together are the **capacity envelope**. Bond is **perfectly plastic**:
when the force reaches what the embedment can develop, the line slips there and
goes on carrying that force. Interface friction does not disappear once it has
been overcome, so nothing drops to zero.

**Rupture** is the second limit, and it is a property of the material rather
than of the burial. `Tmax`, 800 lb/ft here, is the tensile capacity of the
geogrid itself. What happens after it is reached is the `Tres` column: leave it
blank and the line holds its capacity indefinitely while the soil around it goes
on straining, which is the **elastic-perfectly-plastic** model and the industry
default; enter a number and a ruptured element drops to that residual; enter
zero and it ruptures brittly and carries nothing afterwards. We run the blank
case first, because it is what a published reinforced-slope analysis assumes
unless it says otherwise.

![The same layer, seen by each engine](images/fem02_bar_two_engines.png){width=1000}

Both engines apply the same envelope. What differs is where they apply it.

The **limit equilibrium** engine evaluates the envelope at **one point** — where the
trial surface crosses the line — and applies that force to the sliding mass,
tangent to the surface for a flexible support like a geogrid, or along the
line's own axis for a rigid member like a soil nail (the **Dir** column). The
force is prescribed, not computed: nothing about the soil's stiffness, the layer's
stiffness or how far anything moved enters it. That is also why `Tres` has no
meaning in a limit equilibrium run: a post-peak drop depends on how far the
reinforcement has strained, and the method never computes a strain.

The **finite element** engine ties each line into the mesh as a row of tension-only
bar elements sharing nodes with the soil. The force in an element is *EA* times the
strain — the axial stiffness times how much it stretched per unit length —
capped by the envelope value at that element's own position along the line. The
layer therefore develops a force *profile*, low at the ends, high where the
mechanism crosses it, and the load that one element cannot carry is shed into
its neighbors through the soil. That profile depends on two inputs the limit equilibrium run does not use:
the reinforcement's elastic modulus `E` and its cross-sectional area `Area`.
Their product *EA* is the axial stiffness.

[Soil reinforcement in LEM](../lem/reinforcement.md) and
[soil reinforcement in FEM](../fem/reinforcement.md) carry the formulations,
the four end conditions of the envelope, and typical values by reinforcement
type.

---

## Building the mesh

We have to mesh the model before the finite element engine can run. On a
reinforced slope the mesh also sets where each bar element falls along its
line, and with it the capacity that element gets from the envelope. Two of the
Build Mesh defaults change here because of that.

Switch the mode strip to **FEM**. **Run → Run FEM…** stays disabled until a mesh
exists, so build one first. Click **Run → Build Mesh…**

![Build Mesh, set for this reinforced run](images/fem02_studio_build_mesh.png)

Set **Element type** to **Quadratic triangles (tri6)**. A finite element slope
stability analysis needs quadratic elements — linear ones lock and report a
factor of safety that is too high; the
[element type table](../fem/overview.md#element-type-selection-and-volumetric-locking)
in the FEM overview shows how much.

Untick **Auto-size from geometry** and set **Target element size** to `2`. That
is the size the [FEM sample problem](../fem/samples.md) for this model uses, and
it puts 10 elements on each 20 ft geogrid, which is enough to draw the force
profile along a layer. Auto-sizing would divide the 130 ft section width by 100
divisions for a 1.3 ft element instead, which is finer than this mechanism needs
and, as the next two paragraphs show, is not a neutral choice on a reinforced
model.

Untick **Refine thin zones** as well. Checked — and it is checked by default —
the mesher guarantees about four element rows across every thin material zone;
the thin zone here is the `shell`, the cohesive band down the face. We leave the
box off on this model, for a reason particular to this slope. The fill is a
37° sand on a 1.25:1 face — a face slightly steeper than its friction angle —
and only the cohesive band holds the face up. Resolve that band finely enough,
by ticking the box or by shrinking the whole mesh, and the model begins to find
a second, real mechanism: shallow sloughing of the sand face, the same one a
tiny starting circle on the face would find in the limit equilibrium search.
The reported factor of safety then drifts down toward the sloughing answer
(from 1.51 at this page's mesh to 1.43 at 1 ft elements, on the run with a
residual) — a different failure, not a better estimate of this one. We are after
the general mechanism through the reinforced block, the one Spencer's circle
finds, so we stay on the published sample's mesh, where that mechanism
governs. To study a fine mesh and still hold the search on the general
mechanism, the Run FEM dialog has the tools: **Ignore surficial (skin)
failures** with a **Min slip depth**, or an **SSR exclusion** over the face
band, which holds that zone at full strength while the rest of the model is
reduced. Leave the rest of the dialog alone and click **Build**.

The mesh comes out at **4,364 nodes and 2,101 triangles**, with the six lines
discretized into **60 bar elements**, 10 per line:

![The mesh, with the reinforcement lines carried onto element edges](images/fem02_mesh.png){width=1000}

The reinforcement lines are drawn in red because they are part of the mesh, not
an overlay on it: each line was carried into the mesh generator as a constraint,
so every bar element lies on an edge shared with the triangles on either side of
it and its end nodes are soil nodes. The shared node is what couples the two: soil
movement stretches the layer, and the layer's tension pulls back on the soil.

The `shell` shows as a one-element-wide blue veneer down the face, which is the
thin zone the dialog offered to refine. Red triangles along the base are fixed
nodes, blue circles up both ends are rollers, and the green arrows on the crest
are the 240 psf surcharge.

---

## The stiffnesses the run needs

The mesh exists, but the run cannot start on it yet. The starter carries the
limit equilibrium model and nothing beyond it, so every input the continuum
needs is still blank: the two elastic properties on each soil, and two of the
three finite element columns on each reinforcement line. Opening **Run → Run
FEM…** now would show the same refusal we walked through in
[FEM-1](fem01_strength_reduction.md#youngs-modulus-and-poissons-ratio) — errors
naming each blank input, with **Run** disabled — so we go straight to entering
them, soils first, because what a layer's stiffness does depends on the stiffness
of the soil it is buried in.

Click **Materials** in the Inputs dock. On **Table view**, set
the **Show parameters for:** toggles to **FEM** alone, and give both soils an
elastic modulus of 1.0 × 10<sup>6</sup> psf and a Poisson's ratio of 0.3, in the
columns headed `E (psf)` and `n` (the Poisson's ratio column is labeled with a
plain `n`);
[FEM-1](fem01_strength_reduction.md#where-e-comes-from-and-what-it-changes)
covers what each of the two is and where a nominal modulus comes from when the
problem gives you nothing better.

| E (psf) | n |
| :---: | :---: |
| 1000000 | 0.3 |

![The materials table with both soils' elastic properties entered](images/fem02_studio_materials.png)

We use a single round modulus for both layers here rather than values from a
soil-type table, and we measure
[what a stiffness does and does not change](#what-the-two-stiffnesses-change)
further down, because on a reinforced model the answer differs from the one in
FEM-1.

The third finite element column, `t_cut`, already reads 0 on both rows and needs
nothing entered. It caps the tension a soil may carry, and 0 caps it at none,
which is what a fill of this kind holds.
[Tensile strength in the SSRM](../fem/overview.md#tensile-strength-in-ssrm)
covers what the cap does to a strength reduction run and when a material takes a
value above zero. Click **OK**.

The reinforcement's two columns are next. Open **Reinforcement** in the Inputs
dock, and set the **Show parameters for:** toggles to **FEM** alone, as on the
materials table — the limit equilibrium inputs are complete and can stay out of
view. The three columns at the right end, in blue, are the ones only the finite
element engine reads.

![The reinforcement table on the starter: Tres, E and Area empty](images/fem02_studio_reinforce_blank.png)

The six lines are complete on everything both engines share: 800 lb/ft of
capacity, a 4 ft pullout length at each end, no end anchorage, and a `Spacing` of
1 because geogrid properties are already per unit width. The three blue
columns are empty.

Fill `E (psf)` and `Area` on all six rows, and leave `Tres` empty:

| E (psf) | Area |
| :---: | :---: |
| 800000 | 0.1 |

An elastic modulus of 800,000 psf and an area of 0.1 ft² per foot of wall give
an axial stiffness *EA* = **80,000 lb/ft**, which is a typical uniaxial geogrid.
*EA* is what the run actually uses; the two columns are separate so that a
discrete support — a nail or a tieback — can carry a real modulus and a real
cross-section with its spacing dividing the area.
[Axial stiffness (EA)](../fem/reinforcement.md#axial-stiffness-ea) tabulates
values by reinforcement type.

Leaving `Tres` blank makes every line elastic-perfectly-plastic, which is the
first of the two runs. Click **OK** and save the model with **File → Save As…**
under a name of your own.

---

## The iteration budget

Two controls on the Run FEM dialog decide how long each trial may iterate:
**Max iterations per trial**, which opens at 12,000, and the **Iteration
ceiling**, which opens at 50,000. Leave both as they open. The budget extends
itself while a trial is still converging, the Newton corrector finishes most
slow trials at its 300, 1,000 and 3,000-iteration checkpoints, and a trial that
neither settles by the ceiling is reported **inconclusive** and treated as the
undecided edge of the bracket. On this model the answer does not depend on the
budget: a budget of 3,000 returns the same 1.566 and 1.535 as 12,000, from the
same bracket. What a smaller budget does change is the captured failure
picture, because the capture is one more solve that stops at its budget: from
3,000 the elastic-perfectly-plastic mechanism has moved 1.92 ft, from 12,000
7.05 ft, and the smaller one shows a collapse that has barely started. Each run
takes a few minutes on an ordinary laptop.

---

## The elastic-perfectly-plastic run

The first of the two runs is the one that leaves `Tres` blank, which makes every
layer elastic-perfectly-plastic. It gives the finite element answer to compare
against Spencer's before any post-peak assumption is added on top of it. Click
**Run → Run FEM…** again.

![Run FEM on your saved model, Tres still blank: one warning, and Run enabled](images/fem02_studio_run_fem.png)

**Model checks — 1 warning**, and **Run** is enabled. The warning names the thin
`shell` zone, which is the refinement declined during meshing, and it is
expected on this model.

The rest of the dialog opens on the defaults this run wants: **SSRM (find FS)**,
a bracket of 1.00 to 2.00, a tolerance of 0.0100, **Rollers** on the sides, and
**Non-convergence** as the failure criterion — the plain reading, that a trial
which cannot reach equilibrium has failed. In FEM-1 we compare it against the
three other criteria the list offers. Click **Run**.

**FS = 1.566**, from a final bracket of [1.5625, 1.5703], in seven bisection
steps. This is **1.3% below** Spencer's 1.587. The two engines solve the
problem in completely different ways, so agreement within a couple of percent
is a good check on both.

### Viscoplastic shear strain

![The mechanism, and the force in every bar element it crosses](images/fem02_shear_strain_epp.png){width=1000}

The **FEM · Results** tab opens on **At failure**, which is the state the run
captures by re-solving once at 1.15 times the factor of safety so the collapse
develops far enough to draw. The contours are viscoplastic shear strain — the
shearing left after the elastic response is subtracted — and the band they draw
runs from the toe, up behind the reinforced block, and out onto the crest near
**x = 48**. Spencer's circle came out at x = 36.2. The two mechanisms start in
the same place and end about 12 ft apart: the finite element band passes
*behind* the reinforced block, while Spencer's circle cuts through the back of
it. The two engines treat the layers differently — the limit equilibrium run
applies each layer's full capacity at the crossing, the finite element run lets
the force develop from strain — and the critical mechanism each finds follows
from that.

The layers are drawn on the same figure, colored by their axial force against the
second color layer. Every one of them reads deep red — 800 lb/ft, the full
capacity — through its middle, with pale ends of matching length where the two
4 ft pullout ramps allow less. At the captured failure state every line is fully
mobilized: the geogrid is carrying everything it has and the slope still cannot
find equilibrium.

### Reading what each layer is doing

The color layers on the shear strain plot give a first reading of each
geogrid layer's state. The **1D Details** panel gives a much more detailed
one: for every layer, its utilization, the force it carries along its whole
length against its capacity envelope, and a verdict naming how it is working.
Click **1D Details…** on the results toolbar. It opens on the **At failure**
field, where all six layers stand at 100% and read alike; set **Field state**
to **Last converged** — here the trial at *F* = 1.5625 — because that is the
last state that reached equilibrium, and it is where the layers separate.

Each row of the list names a line, gives its utilization, and says what state
the line is in. At that state four of the six — lines 2, 3, 4 and 5, the middle of
the stack — read **yielded**: somewhere away from the ends an element is at the
full 800 lb/ft and holding it, which is the whole tensile strength of the geogrid
mobilized. The other two read **pullout**. Line 1, the bottom layer at the toe,
is the clearest of them: its hardest-worked element sits 1 ft in from the face,
carrying the 200 lb/ft its embedment develops there and slipping rather than
carrying more. That point is working against its bond limit, not its rupture
limit.

Line 5 is one of the four that yield, and it carries both limits at once, so we
open that one:

![Line 5 at the last converged trial, elastic-perfectly-plastic](images/fem02_bar_profile_epp.png){width=1000}

Position along the line runs from the face at s = 0 to the buried end at
s = 20 ft. The dashed line is the capacity envelope: it ramps from zero over the
first 4 ft, sits at 800 lb/ft along the middle of the line, and ramps back to
zero over the last 4 ft. The blue line is the force actually mobilized, with a
point at the center of each of the ten 2 ft bar elements the mesh laid along the
layer — that is where the solver reports each element's force. It climbs from almost nothing at the face, reaches the envelope at s = 9 and holds
it at 11 — those two elements at the full 800 lb/ft are what the panel
titles **yielded** — then falls away down the buried ramp and rides the ramp
itself at s = 17 and 19, where 3 ft and 1 ft of embedment allow only 600 and
200 lb/ft. The rings there are the solver's record that those elements reached
their limit during the run — the force is capped the moment it tries to exceed
the envelope. One line, both limits: rupture in the middle, bond at the tail.

The lower strip is the bond transfer, dT/ds: how much the layer's force
changes per foot of length. The only way force gets into or out of a layer is
through the shear between its surface and the soil, so this strip is that shear,
per foot, and its sign is the direction the soil pulls. From the face out to the
force peak the strip is positive: the soil is moving toward the face and drags
the layer with it, so force builds up with every foot. Beyond the peak it is
negative: the stable soil around the buried part holds the layer, and the force
drains back into the ground foot by foot until the tip. The crossing at zero is
where the pull changes direction, which is why the force is greatest there.

---

## The peak-residual run

The second run changes exactly one thing: what a layer does after it ruptures.
Open **Reinforcement** again and enter `600` in the `Tres` column on all six
rows, leaving everything else as it stands.

![The same table with the residual capacity entered](images/fem02_studio_reinforce_filled.png)

A residual of 600 lb/ft against a peak of 800 is a ratio of 0.75, a little above
the 0.3 to 0.7 usually quoted for geosynthetics. Click **OK**, then
**Run → Run FEM…** and **Run**, with nothing else on the dialog touched.

**FS = 1.535**, from a bracket of [1.5312, 1.5391]. Post-peak behavior cost
**0.031** of factor of safety, 2.0% — half again the 0.021 between
Spencer's 1.587 and the elastic-perfectly-plastic 1.566.

<!-- test: file=files/xslope_reinforced_slope.xlsx, type=fem_ssrm, expected_fs=1.535, element_type=tri6, target_size=2, tolerance=0.01, f_min=1.0, f_max=2.0, benchmark=FEM-2-ssrm -->

Where the cost comes from is plain in the two solutions. Without a residual, a
layer that reaches its 800 lb/ft keeps carrying it, and at the
elastic-perfectly-plastic answer lines 2, 3, 4 and 5 are doing exactly that. With a
residual, a layer that reaches 800 lb/ft drops to 600 and its lost load has to
be picked up by the soil and the other layers; the solver redistributes it and
checks whether the slope still stands. At the peak-residual answer no layer has
ruptured yet: the largest force on any line is 787 lb/ft, on line 2, and every
layer is held by its bond. Reduce the soil strength a little further, to the
next trial in the search, and the first layer reaches 800 and drops — and on
this slope the redistribution cannot be carried:
lines 3, 4 and 5 follow it down and the trial fails. That cascade is the whole
0.031. The residual does not weaken the layers; it means the first rupture is
the last one this slope can afford.

### What changed in the results, and what did not

![The mechanism with the residual in place](images/fem02_shear_strain_pr.png){width=1000}

The band is in the same place, at a lower factor of safety, and the layers are
still red through the middle. What has changed is inside the layers, and the
converged field is where to look for it.

![Line 5 with a residual capacity entered](images/fem02_bar_profile_pr.png){width=1000}

What is new on this profile is the dotted purple line, the **residual capacity**,
which the panel draws because there is one. It is not flat at 600: it runs
along the middle at 600 lb/ft and follows the capacity envelope down each ramp,
because a residual is capped by the bond the embedment can develop at that
point. An element near the end of a layer cannot hold 600 lb/ft after rupturing
when its embedment could only ever develop 200.

Nothing on this line has dropped to it, and nothing should have. A layer drops
only when a converged trial asks an element for more than its 800 lb/ft, and
this trial did not: line 5's interior peak is 763 lb/ft, its tip at s = 19 is
on the 200 lb/ft bond limit — which is why the panel calls the line pullout
rather than yielded — and there is no purple Softened square, the panel's mark
for an element that has shed. The drops happened in the trials above this one,
the ones that failed. Each trial is a separate solve from the same starting
model with every element back at its full 800 lb/ft, so what dropped in a
failed trial does not carry into this one. This state looks like the
elastic-perfectly-plastic one because no layer has ruptured in either; the
residual lowers the answer by deciding the trial above, not this one.

That shows in the other five lines. All six now read **pullout**: each has an
element at its bond limit, slipping there rather than carrying more, and none of
the six reaches the 800 lb/ft that four of them held in the run before. In the
elastic-perfectly-plastic field four lines were at their rupture limit and two
at their bond limit; here the whole stack is at its bond limit. **Pullout** and
**yielded** are both "the line has reached a limit", and only the word says which
limit, which is why the panel prints it. The same panel, left on its default **At failure** state,
shows what the converged state cannot:

![The 1D Details panel, line 5 selected](images/fem02_studio_1d_details.png)

The list on the left carries all six lines with their utilization, a badge
colored by it, and the state each line is in; the map below it shows which line
is selected; and the profile is drawn on the right for the state the **Field
state** selector at the bottom names. It opens on **At failure**, the field of
the capture solve the run makes beyond critical, and here that field carries the
drops. Line 5 reads **softened**: the interior elements at s = 11 and s = 13 sit
on the 600 lb/ft residual line, each marked by a purple *Softened* square, while
their neighbors hold 800. Lines 3 and 4 carry drops too, one element and two, at
the same stations. A slope beyond critical never passes
through equilibrium, so the capture solve cannot decide softening on its own; it
starts from the set the bracket's failed-edge trial shed to at *F* = 1.539. The bond transfer strip shows the consequence: it drops
to −100 lb/ft per ft at 10 ft, where the force falls from 800 to 600, and climbs
back to +100 at 14 ft, where it recovers — the 200 lb/ft the softened elements
shed goes to
the elements on either side. The state where the comparison between
layers lives is still **Last converged**; the at-failure field is where the
residual is seen doing its work.

### The failed state

![The deformed mesh at failure](images/fem02_deformed_pr.png){width=1000}

The deformed panel draws the mesh at its displaced position over a dashed
outline of where it started. The reinforced block has slid out over the toe and
settled at the crest, and the six layers — red where they started gray — are
stretched and rotated with it, hinging where the shear band crosses them. The
title's **Scale = 1.0x** is the exaggeration, which the panel picks so the
collapse reads at this figure size. **Scale ×** and **Auto size** on
the Display panel control that multiplier.

![Displacement vectors at failure](images/fem02_displacement_vectors_pr.png){width=1000}

The vectors give the direction the same field moved in. They point down and
outward at the crest, swing through the body of the slope, and come out nearly
horizontal at the toe, and they die out a few feet behind the buried tips — the
soil beyond that is not part of the mechanism at this reduction factor.

---

## Pullout from the overburden

Everything above states the bond as a **development length**: 4 ft at each end,
with the capacity climbing at a constant rate over it. That number is a
judgement about how deep the reinforcement is buried, and it is the same
judgement for a tail under 16 ft of fill and for a tip at the face. The
`Adhesion` and `Delta` columns, left blank until now, let the model make that
judgement itself.

With them filled, the model no longer needs to be told how far the bond takes
to develop. It reads the resistance at every point along the line from the
soil standing above that point. The interface between soil and reinforcement
is given a strength the same way a soil is: an adhesion *a* and a friction
angle δ, so the shear it can carry per unit area is *a* + σ′<sub>v</sub> tan δ,
where σ′<sub>v</sub> is the effective vertical stress at that depth. A sheet
has soil bearing on both faces, so per foot of length the resistance rate is

*r*(s) = 2 (*a* + σ′<sub>v</sub>(s) tan δ)

and the capacity available at any point is that rate integrated in from the
nearer end, capped at `Tmax`. Deep under the fill the rate is high and the
capacity develops within a foot or two; near the face, under little cover, it
develops slowly.

Where do *a* and δ come from? Without site pullout tests, [FHWA-NHI-10-024](https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi10024/nhi10024.pdf)
gives a default for geosynthetics: the pullout friction factor *F*\* = (2/3) tan φ′,
reduced by a scale-effect correction α = 0.8 for geogrids, with no adhesion.
Its resistance rate, 2 *F*\* α σ′<sub>v</sub>, is the law above with
**Adhesion = 0** and **Delta = arctan(*F*\*α)**. For this sand, φ′ = 37° gives
*F*\* = 0.502, *F*\*α = 0.402, and δ = **22°**.

Open **Reinforcement** once more. Clear the `Tres` column, so the layers are
elastic-perfectly-plastic again and the only change from the first run is the
pullout law, and enter `0` under `Adhesion` and `22` under `Delta` on all six
rows. Once those two are filled, `Lp1` and `Lp2` are no longer read. Click
**OK**, then **Run → Run FEM…** and **Run**.

**FS = 1.535**, against 1.566 with the stated 4 ft lengths — 2.0% lower, four
bisection cells. To see why, open **Results → 1D Details** and select line 2,
the second from the bottom, which is the line the new law changes most:

![Line 2 under the overburden law: a curved capacity envelope](images/fem02_bar_profile_law.png){width=1000}

The envelope stops being a pair of straight ramps. σ′<sub>v</sub> grows with
depth of cover, so on a line running back from the face the resistance per foot
grows with it and the integral of a growing rate is a curve: on line 2 the law
allows **42, 168 and 378 lb/ft** at 1, 2 and 3 ft in from the face, where the
stated 4 ft ramp allowed 200, 400 and 600. At the buried end it goes the other
way — under 16 ft of fill the law reaches the full 800 lb/ft within half a foot
of the tip, and holds it at 17, 18 and 19 ft where the stated ramp is ramping
back down through 600, 400 and 200. One foot in from a tip, the two ends of the
same line are allowed 42 lb/ft and 800 lb/ft — a factor of nearly twenty, where
the stated length allows 200 at both.

The run puts line 2 at capacity in one place only — the element 3 ft in from
the face, riding the law's 378 lb/ft — while its interior peaks at 797 lb/ft
9 ft in, just under the 800 the envelope allows there. Under the stated
lengths the same line held the full 800 lb/ft over five interior elements,
from 7 to 15 ft. That is the engineering content of the
comparison: a stated development length is depth-blind, so it under-rates a
deeply buried tail and over-rates a shallow one, and reinforced slopes are
critical at the face.

The factor of safety barely moves, and the field shows why that number
understates what changed:

![The mechanism at failure under the overburden law](images/fem02_shear_strain_law.png){width=1000}

The band is where it was — from the toe up behind the buried tips to the
crest, its center within half a foot of the stated-length run's — and the
strain in it is more concentrated, 52 elements above half the peak against 117.
The story is in the colors on the layers. Under the stated lengths every layer
faded from red to white over its last 4 ft, because the ramp allowed less and
less toward the tip, and the band cut through those weakening tails. Under the
law the layers are red to the tip: the buried ends develop the full 800 lb/ft
within a foot, and the band passes just outside them. At the face the reverse
happens — the outer few feet of every layer are blue and white, because under
little cover the law allows almost nothing there, where the stated ramp allowed
200, 400 and 600. The factor of safety reads one number for both fields; the
layers say the weak end moved from the tail to the face.

Both answers move, and by comparable amounts: the finite element one by 2.0%,
and Spencer, searched again on the same file, from 1.587 to **1.559**, 1.8%.
Spencer's move is not a loss of capacity. On the critical circle from the
stated-length search, every one of the five crossings sits 8.4 ft or more from
the nearer end, out where both laws allow the full 800 lb/ft, so that surface
cannot tell them apart.
What happens is that the critical circle shifts slightly — it now daylights
8.4 ft beyond the toe rather than at it, and clips line 1 two feet from the
face, where the law offers 168 lb/ft against the stated ramp's 400. The red dots on the
solution plot are where each line's available tension first reaches the full
800 lb/ft from either end — the ends of its development lengths. Under the
stated lengths they sat 4 ft in from both ends of every line; under the law they
sit where the curved envelope reaches 800, about 4.4 ft in at the face and
within a foot of the buried tip, so the stretch a line carries its full
capacity over runs almost to its tail. One sentence of caution on
σ′<sub>v</sub>: it is the weight of the soil column standing above the point and
nothing else, so the 240 psf crest surcharge does not count toward pullout
resistance — which is what FHWA directs for a live load.

![Spencer's critical circle under the overburden law: daylighting beyond the toe and clipping line 1 near the face](images/fem02_lem_solution_law.png){width=1000}


---

## The three answers together

Our three runs are finished, and all three land within 4% of one another —
close, for two engines and two assumptions about the layers. What separates
them is still worth knowing, and it has two sources: the post-peak assumption,
which the section above already measured, and the different ways the two
engines decide what force a reinforcement line carries.

Three readings of the same slope, from the same file:

| Reading | FS |
| --- | :---: |
| Spencer's method, searched on this page | 1.587 |
| Strength reduction, layers elastic-perfectly-plastic | 1.566 |
| Strength reduction, layers with a 600 lb/ft residual | 1.535 |

The gap between Spencer and the peak-residual run is 0.052, and the section
above measured what post-peak behavior costs: **0.031 of it**. The other
0.021 separates Spencer from the elastic-perfectly-plastic run, which is the
same physical assumption about the layers — hold capacity once yielded — applied
two different ways.

The difference is where the force is decided. Spencer takes the envelope value
at one crossing point and hands the sliding mass 800 lb/ft, five times over, as
a known force on a surface it chose. The finite element run lets the force
emerge along the whole of each layer from displacement compatibility, so a layer
contributes what the soil around it actually mobilized: line 1 carries 521 lb/ft
at its most-worked point against the 800 the envelope allows, all six lines are
held at their bond limits once the residual is in play, and the
force in every layer tapers away from the shear band instead of standing at its
envelope value everywhere. Prescribing the maximum available force at one point
is the more generous of the two, and 1.8% is what that generosity amounts to on
this slope.

The mechanisms differ in the same direction. Spencer's circle had to be a
circle, and the most critical circle available cuts out through the crest at
x = 36.2. The finite element band, free to take any shape, goes around the back
of the reinforced block and reaches the crest near x = 49 — a path no circular
search had on offer.

---

## What the two stiffnesses change

In [FEM-1](fem01_strength_reduction.md#where-e-comes-from-and-what-it-changes)
we measured that a hundredfold sweep in Young's modulus left the factor of safety
identical to every printed digit. A reinforced slope has a second stiffness in
it — the layers' — and the two do not have to move together, so that result is
checked here rather than assumed.

The same run, at the same mesh and bracket, with the soils' modulus and the
layers' modulus scaled separately:

| Soil E | Layer E | FS |
| :---: | :---: | :---: |
| ×0.1 | ×1 | 1.5273 |
| ×1 (the file) | ×1 | 1.5352 |
| ×10 | ×1 | 1.2305 |
| ×10 | ×10 | 1.5352 |

The result does not carry over. Ten times the soil's modulus, with the geogrid
left as it is, costs **0.305 of factor of safety — 20%**. It moves in the
direction the mechanism says it should: a stiffer soil strains less under the
same load, the layers are stretched less, and the reinforcement is never asked for
the tension it is capable of. One tenth of the modulus costs one bisection
cell, next to nothing, because at that stiffness the layers are already carrying
close to everything the envelope allows and there is little further to mobilize.

Scale the layers by the same ten as the soil and the answer returns to 1.5352, to
the printed digit and from the same bracket. What a reinforced finite element
model responds to is the **ratio** of soil stiffness to layer stiffness, not the
absolute value of either. The invariance measured in FEM-1 holds only where there
is one stiffness to sweep; put a second one in the model and the pair matters. A round
soil modulus against a catalog number for the geogrid is a modeling decision with
a real consequence for the factor of safety, and the two values should be chosen
against each other rather than one at a time.

---

## How much the residual strength changes the answer

The residual we entered above, 600 lb/ft, was one choice out of a range. Running
the same model across the column — and running it blank, and running it at zero,
which is brittle rupture — measures how much the answer depends on that choice:

| `Tres` | FS |
| :---: | :---: |
| blank | 1.5664 |
| 800 | 1.5664 |
| 600 | 1.5352 |
| 400 | 1.5195 |
| 0 | 1.4570 |

The answer steps down once for each drop in the residual. `Tres` = 800
reproduces the blank run exactly — the same factor of safety, the same bracket,
the same trials — which is the model saying what it should: a residual equal to
the peak means nothing ever drops below what it was already carrying, so the
post-peak branch is never entered. Below that the answer keeps falling, in uneven steps: 600 costs
0.031, 400 another 0.016, and a brittle zero a further 0.063.

Two things follow for a real design. The first is that the size of the residual
matters as much as its presence: the step from holding capacity to shedding to
three quarters of it is 2.0% of factor of safety, and the whole way down to a
brittle zero costs three times that again. The second is that the
whole range, from a blank cell to a brittle zero, spans 0.109, about 7%. Leaving `Tres` blank claims the geogrid
holds its capacity once it yields, which is what the mainstream codes assume and
what most published capacities describe; entering zero claims it snaps. Neither
claim is one a catalog value settles, and the meshing step showed that the answer
under either moves with the discretization, so a residual is better treated as a
back-analysis parameter than as a design default.

---

## Conclusion

This tutorial covered:

- The two limits every reinforcement model carries — bond, which limits what the
  embedment can develop and slips perfectly plastically once reached, and
  rupture of the reinforcement itself, whose aftermath is the `Tres` column.
- What each engine does with those limits: a limit equilibrium method applies
  the envelope as a prescribed force at one crossing, and a finite element method
  meshes the line into bar elements whose force emerges from the movement of the
  soil around them.
- The three finite element reinforcement inputs, the axial stiffness *EA* two of
  them amount to, and the model checks that will not start a run without them.
- The per-trial iteration budget, the ceiling it extends to, and why neither
  decides the answer on this model — only how far the captured failure state
  develops.
- Reading a reinforced result: the mechanism out of shear strain, and each
  layer's force profile, bond transfer and verdict out of the 1D details panel at
  the last converged trial.
- How the two engines' answers stand against each other, a residual capacity that
  lands on steps rather than a curve, and what changes when the bond is read from
  the depth of burial instead of a stated length.

**Where to go next:** [Soil reinforcement in FEM](../fem/reinforcement.md) and
[soil reinforcement in LEM](../lem/reinforcement.md) carry both formulations; in
[LEM-8](lem08_reinforced_slope.md) we build this model from nothing and measure
what the geogrid adds to the bare section, and in
[FEM-1](fem01_strength_reduction.md) the method is run on an unreinforced
embankment. In [FEM-3](fem03_piles.md) we put stabilizing piles through the same
two-engine comparison, and measure which engine a discrete row belongs in and
which a continuous wall does.
