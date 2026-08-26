---
title: "Tutorial COMBO-2 — Rapid Drawdown"
description: "The Johnson Reservoir dam with its pool taken from elevation 160 to elevation 100 — the Duncan, Wright and Wong three-stage procedure run on one circle from three different statements of where the water is: a sketched piezometric pair, two steady seepage solutions, and two frames of a transient march."
---

# Tutorial COMBO-2 — Rapid Drawdown

[COMBO-1](combo01_seepage_stability.md) ran the Johnson Reservoir dam under a
standing reservoir. This page takes that reservoir away.

Lowering a pool removes the water pressing on the upstream face, which was
holding the slope up. If the soil behind the face could shed its own pore water
just as fast, nothing much would happen: the driving weight and the resisting
strength would fall together. A compacted clay core cannot. It keeps the pore
pressure the full reservoir put into it while the load that balanced that
pressure disappears, and the upstream slope is left carrying its own weight
against a strength it no longer has. Rapid drawdown is the analysis of that
window.

The procedure needs two states of the water — before the drawdown and after —
and there are three ways to say where the water was in each. This page runs all
three on one dam, one slip circle and one set of strengths, so the only thing
that changes between the answers is the statement about the water.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium + seepage</p></div>
<div class="tgt-tile"><span class="tg-label">Build &amp; run</span><p>~35 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn what the three-stage rapid-drawdown procedure computes and
what it needs on the inputs; how to supply its two water states from a
piezometric pair, from two steady seepage solutions and from two frames of a
transient march; how to read which stage governed; and what changes the answer
between the three sources.
</div>
<p><span class="tg-pill">three materials</span><span class="tg-pill">rapid drawdown</span><span class="tg-pill">three-stage procedure</span><span class="tg-pill">Kc = 1 envelope</span><span class="tg-pill">d and ψ</span><span class="tg-pill">piezometric lines</span><span class="tg-pill">second boundary set</span><span class="tg-pill">transient seepage</span><span class="tg-pill">stage times</span><span class="tg-pill">automatic water loads</span><span class="tg-pill">Spencer</span><span class="tg-pill">parametric sweep</span></p>
<div class="tgm-model" markdown>
**Starter file** — [xslope_johnson_rapid_start.xlsx](files/xslope_johnson_rapid_start.xlsx),
the dam with its strengths, its undrained core envelope and a sketched
piezometric pair; no seepage boundary sets and no schedule, which is what the
build below adds

**Completed model** — [xslope_johnson_rapid.xlsx](files/xslope_johnson_rapid.xlsx),
the same dam with both seepage boundary sets and the pool schedule on it; open it
to skip the construction and start at [Two steady solutions](#two-steady-solutions)
</div>
</div>

---

## What rapid drawdown is

A reservoir standing against an upstream slope does two things at once. Its
weight presses on the face, and that pressure is a stabilizing load. And it
holds a head inside the embankment, which raises the pore pressure and lowers
the effective stress the strength is computed from. Drawing the pool down
removes the first immediately and the second only as fast as the soil can drain.

How fast that is comes from the dimensionless time factor
$T = c_v t / D^2$, where $c_v$ is the coefficient of consolidation, $t$ the time
over which the pool falls and $D$ the longest path pore water has to travel to
get out. [Rapid drawdown analysis](../lem/rapid.md#when-does-rapid-drawdown-apply)
carries the rubric and the $c_v$ table: above $T = 3$ the material drains as fast
as the pool falls and the drawdown is not rapid for it; below 3 it does not, and
the material has to be carried through the drawdown at an undrained strength.

The test is applied per material, not per dam, which is what makes a zoned
embankment the interesting case. This dam's sand shell and silty-sand foundation
drain freely over a 45-day drawdown; its compacted-clay core does not. So the
core is the one zone that arrives at the end of the drawdown still holding the
reservoir's head, and it is the zone that governs.

### The three-stage procedure

XSLOPE uses the three-stage procedure of Duncan, Wright and Wong (1990),
documented in full on the [rapid drawdown](../lem/rapid.md) page. In outline:

**Stage 1 — full pool, drained.** An ordinary effective-stress analysis of the
slope before the drawdown, with the reservoir load on the face and the
pre-drawdown pore pressures inside. Its factor of safety is not the answer. What
it produces is the effective normal stress $\sigma'_{fc}$ and the mobilized shear
stress $\tau_{fc}$ on the base of every slice — the consolidation stress state
the soil arrives at the drawdown in. The slope has to be stable here: at
$FS < 1$ the mobilized stress lies above the failure envelope, there is no
equilibrium consolidation state to read, and the run is
[refused](../lem/rapid.md#the-stage-1-slope-must-be-stable) rather than answered.

**Stage 2 — drawn down, undrained.** The same surface with the drawn-down water
— the reservoir load gone and the post-drawdown pore pressures inside — and, in
every material carrying a $K_c = 1$ envelope, an undrained shear strength
$\tau_{ff}$ computed from that material's stage-1 stresses. The envelope is
entered as two numbers, an intercept $d$ and a slope $\psi$, and the strength
actually used is interpolated between it and the drained envelope according to
the stress ratio $K_1$ the slice consolidated at. Materials with no $d$ / $\psi$
keep their drained $c'$ and $\phi'$ through this stage.

**Stage 3 — drawn down, drained.** The same drawn-down section again, with the
drained strength substituted on every slice where the drained strength turns out
to be *lower* than the undrained one stage 2 used. Because stage 3 only ever
lowers a strength, its factor of safety can never exceed stage 2's, and the
reported answer is the lower of the two. If no slice qualifies, stage 3 is not
run and stage 2 is the answer.

Which of the two governs is a result rather than bookkeeping. Stage 2 governing
says the undrained strengths control: on every slice crossing the core the
undrained strength is the lower of the two, and the drawdown is a strength
problem. Stage 3 governing says the opposite — on at least some of those slices
the undrained strength assigned in stage 2 came out *higher* than the drained
strength at the same stresses, so what limits the slope is the ordinary drained
strength at the drawn-down pore pressures. The
[flip section](#the-governing-stage-flip) below moves this dam from one to the
other and reports what it takes.

### The procedure by hand

Duncan, Wright and Brandon (2014) work the whole procedure by hand on a
submerged infinite slope, and the intermediate quantities are printed in their
Table 9.2. The slope is 3:1 ($\beta = 18.4°$) under 100 ft of water,
$\gamma = 125$ pcf, $c' = 0$, $\phi' = 40°$, with a $K_c = 1$ envelope of
$d = 2000$ psf and $\psi = 20°$. Two slip planes are checked, one 5 ft deep and
one 30 ft deep, and the drawdown is total — the water is fully removed and the
pore pressures afterward are zero.

| | | $z$ = 5 ft | $z$ = 30 ft |
|---|---|---:|---:|
| **Stage 1** | Total normal stress, $\sigma_{fc}$ | 6834 psf | 9803 psf |
| | Pore water pressure, $u$ | 6552 psf | 8112 psf |
| | Effective normal stress, $\sigma'_{fc}$ | 282 psf | 1691 psf |
| | Shear stress after consolidation, $\tau_{fc}$ | 94 psf | 562 psf |
| **Stage 2** | Shear strength, $\tau_{ff(K_c=1)}$ | 2103 psf | 2615 psf |
| | Shear strength, $\tau_{ff(K_c=K_f)}$ | 237 psf | 1419 psf |
| | Consolidation stress ratio, $K_1$ | 2.0 | 2.0 |
| | Undrained strength interpolated, $\tau_{ff}$ | 1585 psf | 2283 psf |
| | Shear stress after drawdown, $\tau$ | 187 psf | 1123 psf |
| | Factor of safety | **8.48** | **2.03** |
| **Stage 3** | Effective normal stress after drawdown | 563 psf | 3376 psf |
| | Drained shear strength | 472 psf | 2833 psf |
| **Final** | Governing strength | 472 psf (drained) | 2283 psf (undrained) |
| | Factor of safety | **2.52** (stage 3) | **2.03** (stage 2) |

*Duncan, J. M., Wright, S. G., and Brandon, T. L. (2014),* Soil Strength and
Slope Stability, *2nd ed., Wiley, Table 9.2 (p. 175).*

Two planes in the same soil, and the two stages split them. At 5 ft the
consolidation stress is small, so the $K_c = 1$ envelope's 2000 psf intercept
dominates and the interpolated undrained strength of 1585 psf is more than three
times the 472 psf the drained envelope allows at the post-drawdown stress: the
drained strength is lower, stage 3 runs, and the factor of safety falls from 8.48
to 2.52. At 30 ft the consolidation stress is six times larger, the drained
strength has grown to 2833 psf while the undrained one has only reached 2283, and
stage 3 has nothing to substitute. Nothing about the soil differs between the two
— only the depth, and with it the stress the soil consolidated at.

### What the procedure needs on the inputs

Three things, beyond an ordinary limit equilibrium model:

| Input | Where | What it does |
|---|---|---|
| $d$ and $\psi$ | Materials table, `d` and `psi` | The $K_c = 1$ envelope, on every material that does not drain during the drawdown. Both or neither: a material carrying only one of them reverts silently to drained. |
| A second water state | Piezometric Line 2, or a second seepage solution, or a transient frame | The pore pressures stages 2 and 3 read |
| The drawn-down water load | `main!D23` **Water loads** = `auto` | The reservoir load on the face at stage 1 and whatever remains of it at stage 2, both derived from the model's own water definition |

The second and third are the same statement seen twice, which is why leaving
**Water loads** on `auto` matters here more than anywhere else: the engine reads
the pool from the same place at both stages, so the water pressing on the face
and the pore pressure inside the dam can never disagree about where it stood.

---

## The dam

![The Johnson Reservoir dam](images/seep02_problem_sketch.png){width=1000}

The section is 750 ft long: a 100 ft foundation on rock at elevation 0, and an
80 ft embankment on it with a crest at elevation 180. A sand shell rises 2:1
upstream and falls at about 2.1:1 downstream, with a compacted-clay core through
the middle carried 40 ft into the foundation as a cutoff key. The reservoir
stands at elevation 160 and the tailwater at 100.
[SEEP-2](seep02_johnson_dam.md) builds the model and covers the geometry in full.

Download [xslope_johnson_rapid_start.xlsx](files/xslope_johnson_rapid_start.xlsx)
and open it with **File → Open…**. The toolbar's mode strip reads
LEM | Seepage | FEM; leave it on **LEM** for now.

### The materials, and the one undrained zone

Click **Materials** in the Inputs tree and set the **Show parameters for:**
toggles to **LEM** alone:

![The three zones, with d and psi on the core alone](images/combo02_studio_materials.png)

| mat | name | γ (pcf) | c′ (psf) | φ′ (deg) | d (psf) | ψ (deg) | u |
|:---:|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | `shell` | 130 | 100 | 35 | — | — | `piezo` |
| 2 | `core` | 125 | 400 | 18 | 250 | 14 | `piezo` |
| 3 | `foundation` | 127 | 100 | 27 | — | — | `piezo` |

The **d** and **psi** columns are what make this a drawdown model, and only the
core carries them. Its $K_c = 1$ envelope is $d = 250$ psf and $\psi = 14°$, well
under its drained envelope of $c' = 400$ psf and $\phi' = 18°$: at the stresses
this circle produces it makes the core weaker undrained than drained, which is
why stage 2 is the stage that governs every run below.
[The sweep](#the-governing-stage-flip) at the end of the page finds the value of
$d$ that reverses it. The shell and the foundation are left blank on both
columns, which declares them free-draining: they keep 100 psf and 35°, and
100 psf and 27°, through every stage.

That declaration is reported back rather than assumed. Every run on this page
carries one warning — *"2 of 3 material(s) a failure surface can cross carry no
d / psi and are treated as free-draining through the drawdown, keeping their
drained strength"* — which is the model saying what it understood. The neighboring
case is an error rather than a warning: a material carrying **d** without
**psi**, or the reverse, would be treated as free-draining with nothing said, so
the checks refuse the run rather than let a half-entered envelope pass.

### The two piezometric lines

Click **Piezometric lines**. The editor has two tabs, and both are filled:

![Line 2, the drawn-down surface, with Line 1 dimmed behind it](images/combo02_studio_piezo.png)

**Line 1 — the full-pool surface.** Five points, running from the reservoir
surface at elevation 160, across the upstream shell at that level, down through
the core, and out to the tailwater at 100:

| x | y |
|---:|---:|
| 0 | 160 |
| 360 | 160 |
| 410 | 120 |
| 550 | 100 |
| 750 | 100 |

**Line 2 — the surface at the end of the drawdown.** Six points. The pool is now
at elevation 100, so the line starts there and stays there across the upstream
foreshore; it climbs to 150 over the core, which is still holding most of the
head the full reservoir put into it; and from x = 410 out the two lines are
identical, because the drawdown is an upstream event and the downstream shell has
not felt it yet:

| x | y |
|---:|---:|
| 0 | 100 |
| 200 | 100 |
| 360 | 150 |
| 410 | 120 |
| 550 | 100 |
| 750 | 100 |

Five and six points is the whole of it, and the coarseness is deliberate. **A
piezometric line is a statement about the water, not a trace of a solved
field.** In practice Line 1 comes
from piezometers in a dam that has been standing at full pool, and Line 2 comes
from judgment about how far the core will have drained by the end of the
drawdown — a straightedge sketch between a few readings and a few assumptions.
The two above were read off the solved seepage fields further down this page and
then rounded to something a designer would draw, which is why they follow those
fields in shape without matching them anywhere in particular. Where a seepage
solution is available it replaces them, and the rest of this page measures what
that replacement changes.

Two rules govern the pair. Line 2 may not stand above Line 1 anywhere — a
drawdown that raised the pool is not a drawdown, and the checks say so. And
Line 1 has to describe the same pool the rest of the model does: on this file it
leaves the upstream shell at exactly elevation 160, which is where the seepage
boundary will put the reservoir when it is built later. A line that came off the
shell a foot low would put two different pools on one model, and the checks flag
that too.

Both lines run the full width of the section, and the tail from x = 550 out to
750 sits exactly on the ground surface. That matters for the water load: with
**Water loads** on `auto` the engine measures each line against the ground and
turns whatever stands above it into a distributed load, so a tail a few feet high
over the downstream foreshore would invent a pond that is not there.

Click **OK**, and the Inputs plot draws the whole model:

![The model the drawdown runs on](images/combo02_inputs.png){width=1000}

Three profile lines carve the section into shell, core and foundation. The heavy
dark blue line is Piezometric Line 1 and the pale one Line 2, each carrying the
inverted-triangle water-table symbol at its midpoint: 60 ft apart over the
upstream foreshore, still 22 ft apart where the upstream slope meets full pool at
x = 320, and the same line from x = 410 out.

The blue arrows on the upstream face and across the foreshore are the **derived
water load** — the reservoir standing against the slope, 3744 psf at the toe
under 60 ft of water and tapering to zero at the elevation-160 waterline, derived
from Line 1 rather than entered anywhere. There is no second set of arrows,
because Line 2 stands nowhere above the ground surface: at the end of this
drawdown no water is left on the slope at all, and the removal of that 3744 psf
is the load the whole analysis is about. The red dashed arc is the circle every
run below uses.

---

## The drawdown from the piezometric pair

Switch to **LEM** (`Ctrl+1`) if the mode strip is elsewhere, and click
**Run → Run LEM…**

![Run LEM, with Rapid drawdown ticked](images/combo02_studio_run_lem.png)

Three fields change from the dialog's defaults:

**Method → Spencer.** Spencer satisfies both force and moment equilibrium, and
every number on this page is Spencer's.

**Analysis → Single surface.** The dialog opens on **Auto search**, which finds
each run its own critical circle. That is the right setting for a design check
and the wrong one here: this page compares three statements about the water, and
if each run is allowed to move to its own surface the comparison carries a change
of geometry as well. **Single surface** runs the circle the file carries — center
(275, 235), radius 160 ft, tangent to elevation 160 — which was located by a
rapid-drawdown search on this dam and sits on the upstream face, toeing near the
upstream base and daylighting just past the crest. All three runs below use it.

**Rapid drawdown → ticked.** This is the control that runs the three-stage
procedure instead of a single solve, and puts the run through the drawdown checks
as well as the ordinary limit equilibrium ones. On this file the checks report
one warning, the free-draining declaration read back off the materials table.

There is no **stage time** field on the form yet. Those belong to a model that
carries a transient schedule, and appear when
[one is built](#the-pool-schedule-and-the-stage-times).

Leave **Number of slices** at 40 and click **Run**. The run returns immediately,
and the Log pane carries the three-stage report:

```text
Running SPENCER — single circular surface (Xo=275.0, Yo=235.0, R=160), 40 slices (rapid drawdown)…
Generated 42 slices; solving…
=== RAPID DRAWDOWN ANALYSIS ===
Stage 1: Pre-drawdown conditions...
Stage 1 FS = 1.9354
Stage 2: Post-drawdown conditions with undrained strengths...
Stage 2 FS = 1.3010
Stage 3: Checking drained strengths...
Stage 3: No drained strength adjustments needed
Final rapid drawdown FS = 1.3010
=== END RAPID DRAWDOWN ANALYSIS ===
Spencer: FS=1.301, theta=12.89
FS = 1.301
Sliding mass = 1,553,767.2 lb/ft over 266.45 ft of failure surface
```

**42 slices, not the 40 asked for.** The slicer puts a boundary wherever a
feature crosses the surface, and the piezometric lines' own vertices are
features, so two extra boundaries appear. The same 42 slices come out of every
run on this page, 14 of them cutting the core.

**Stage 3 was not run.** *"No drained strength adjustments needed"* means that on
all 14 core slices the drained strength at the stage-2 stresses came out higher
than the undrained strength stage 2 used, so there was nothing lower to
substitute. Stage 2 governs, and **the rapid-drawdown factor of safety is 1.301**
against a full-pool 1.935.

![The three-stage answer from the piezometric pair](images/combo02_solution_piezo.png){width=1000}

The figure is the governing stage's own state, not stage 1's: the pale blue band
under the failure surface is the stage-2 pore pressure on each slice base,
reaching 3499 psf, and the green hatched band above it is the effective normal
stress the stage-2 strengths were computed from. The stage-1 reservoir load is
still drawn on the upstream face, because that is where the consolidation
stresses came from.

The analysis report writes the same three stages under **Rapid Drawdown**,
together with which of the two drawn-down stages the reported factor came from.

---

## Two steady solutions

A piezometric line says where the water table is. A seepage solution says what
the head is everywhere, which is a stronger statement: it accounts for the flow
through the core, the anisotropy of the zones, and the seepage face on the
downstream slope, none of which a line drawn between piezometers knows about.
Rapid drawdown can be given two of them — one at full pool and one at the
drawn-down pool — and that is what this section builds.

### The second boundary set

Switch the mode strip to **Seepage** (`Ctrl+2`) and click **Seep BC** in the
Inputs dock. The editor has two tabs.

**Set 1** is the ordinary boundary set: it describes the dam under the full
reservoir. Press **Add head** twice, then select the **Exit face** entry already
in the list.

*Head 1 — the reservoir*, at **Head value (ft)** `160`, along the upstream
foreshore and up the submerged part of the face:

| x | y |
|---:|---:|
| 0 | 100 |
| 200 | 100 |
| 320 | 160 |

*Head 2 — the tailwater*, at `100`, along the downstream foreshore:

| x | y |
|---:|---:|
| 550 | 100 |
| 750 | 100 |

*Exit face*, over the whole downstream slope, from the crest to the tailwater:

| x | y |
|---:|---:|
| 380 | 180 |
| 550 | 100 |

[SEEP-2](seep02_johnson_dam.md#4-boundary-conditions) is where the three boundary
types and the no-flow default are worked through, and it solves this same set.

**Set 2 (rapid drawdown)** is a second, independent boundary set on the same
section. Its job is to describe the dam under the drawn-down pool, so that a
second steady solve can be made from it. Click the tab and enter two heads:

![Set 2, the drawn-down boundary set](images/combo02_studio_seep_bc2.png)

*Head 1 — the drawn-down pool*, at `100`, along the upstream foreshore:

| x | y |
|---:|---:|
| 0 | 100 |
| 200 | 100 |

*Head 2 — the tailwater*, at `100`, unchanged:

| x | y |
|---:|---:|
| 550 | 100 |
| 750 | 100 |

Three things about Set 2 decide what it computes.

**It carries plain heads only.** The `reservoir` type — the submerged-only
boundary built for a pool that moves — and any value bound to a time series both
belong on Set 1, because there is one timeline in a model and it is Set 1's. Set
2 is a constant, steady state by construction, and the editor refuses the other
kinds.

**It has no exit face.** With the pool drawn to the tailwater datum, both heads
are at elevation 100, there is no head difference across the section, and no
water flows. An exit face is the boundary where water *may* leave; with nothing
to leave, drawing one asks the solver to decide which of its nodes are
discharging when none are. Leaving it off states the problem as what it is.

**The pool is drawn all the way down.** Elevation 100 is the tailwater datum and
the upstream toe, so at the end of this drawdown no water remains against the
upstream face at all. That is the most severe version of the load removal, and it
is also what makes Set 2's answer the fully re-equilibrated limiting case: the
head everywhere settles to 100, the embankment above that elevation is dry, and
the core has given up everything the reservoir put into it.

Click **OK**.

### Meshing, and both solves

Click **Run → Build Mesh…** and set **Element type** to **Linear triangles
(tri3)**. Head is a scalar field, so there is nothing for a linear element to
lock up on, and no strength reduction runs on this page — the
[quadratic requirement COMBO-1 meets](combo01_seepage_stability.md#one-mesh-for-all-three-analyses)
comes from the finite element stability engine, not from seepage. Element order
also sets the cost of the transient march below, which takes thousands of steps
on whatever mesh it is given.

Leave **Auto-size from geometry** ticked at **100** size divisions, which makes
the target element size the width of the section over that number,
750/100 = 7.5 ft. Click **Build**. The mesh comes out at **2,080 nodes and 3,923
triangles**:

![The mesh, and Set 1's boundary conditions on it](images/combo02_mesh.png){width=1000}

The blue squares are Set 1's specified-head nodes — along the upstream foreshore
and up the face to elevation 160, and along the downstream foreshore at 100 —
and the red circles down the downstream slope are the exit face.

Click **Run → Run Seep…**, leave **Convergence tol** at `0.0001` and **Max
iterations** at `400`, and click **Run**. Because the file defines two boundary
sets, one steady run solves **both** of them, each keeping its own results tab,
and each writing a companion file beside the workbook — `_seep.csv` for Set 1 and
`_seep2.csv` for Set 2 — so that reopening the file picks them up by name.

(Opened from the completed model rather than built up from the starter, this
dialog also carries a **Run type** selector, because that file already has the
schedule the next section adds. Leave it on **Steady** for this run.)

![The two steady solutions, on one color scale](images/combo02_steady_pair.png){width=1000}

**Set 1** settles in **22 unconfined sweeps** to a total discharge of
**1.9566 ft³/day per ft** of dam. The head runs from 100.000 to 160.000 ft, the
two boundary values and nothing outside them, and the heavy black line is the
phreatic surface: nearly flat across the upstream shell at reservoir level,
falling almost the whole 60 ft inside the core, and running low across the
downstream shell to the tailwater. Almost the entire head loss is in the core,
which is what a cutoff key is for.
([COMBO-1](combo01_seepage_stability.md#solving-it) reports 1.925 for the same
boundary set on a quadratic mesh at the same 7.5 ft target — 8,082 nodes rather
than 2,080, because quadratic elements add a node on every edge. The 1.6%
between the two answers is that discretization, not the physics.)

**Set 2** needs no iteration at all. With both heads at 100 and no flow, the
answer is a head of **100.000 ft everywhere** and a discharge of zero. With no
exit face the solve is linear rather than unconfined, and no phreatic surface is
drawn. What the stability run reads off it is a pore pressure that is hydrostatic
below elevation 100 — and above it, once the negative pressures over the water
table are
[clamped to zero](../seep/seep_slope.md#negative-pore-pressures), nothing at all.

### Pointing the materials at the fields

The seepage solutions exist, and nothing reads them yet. Switch back to **LEM**
(`Ctrl+1`), open **Materials**, and set the **u** column to `seep` on all three
rows. That column is what sends a solved field to every slice base, and it is set
per material rather than once for the model;
[COMBO-1](combo01_seepage_stability.md#the-column-that-connects-the-modes)
covers its four values and what each one costs.

A drawdown run reads **two** fields through that one column: `seep_u` from Set 1
for stage 1, and `seep_u2` from Set 2 for stages 2 and 3. The checks hold the run
to both. A material on `seep` with no solved field at all is an error; a material
on `seep` with only Set 1 solved is a second error, naming the missing
drawn-down field by the file name it would arrive under — because a stage-2 pore
pressure that reads zero everywhere is a drawdown that removed the water *and*
its pressure, and it returns a factor of safety that is too high.

Run **Run → Run LEM…** again with the same three settings — Spencer, Single
surface, Rapid drawdown ticked:

| Stage | | FS |
|---|---|---:|
| 1 | full pool, drained | 2.0542 |
| 2 | drawn down, undrained core | **1.5640** |
| 3 | drawn down, drained | not required |

![The three-stage answer from two steady solutions](images/combo02_solution_steady.png){width=1000}

**The rapid-drawdown factor of safety is 1.564**, against 1.301 from the
piezometric pair. Stage 2 governs again. The stage-1 answer has moved too — 2.054
against the pair's 1.935 — because the solved full-pool field puts a little less
pore pressure on the deep part of the surface than a static column measured down
from Line 1 does. The much larger change is at stage 2: the peak pore pressure on
the 42 slice bases falls from 3499 psf to **1559 psf**, all of it in the core,
which in this field has fully drained.

---

## A transient march

Two independent steady solutions describe two equilibrium states. What they do
not describe is the passage between them: Set 2 is the dam *long after* the
drawdown, when the core has finished draining. On a core with a conductivity of
0.001 ft/day that is a long time, and the drawdown itself takes 45 days. A
transient seepage run resolves that directly — it is a history of pore-pressure
fields as the reservoir is lowered on a schedule, so the two stage fields are two
frames read out of one march rather than two independent solves assembled by
hand.

[SEEP-3](seep03_reservoir_drawdown.md) builds a transient model from nothing and
covers the storage properties, the time series, the initial condition and the
saved-frame schedule. The build below is the same shape on this dam, and the
storage properties it needs are already on the materials table:
*S<sub>s</sub>* and *S<sub>y</sub>* of 1 × 10<sup>−4</sup> /ft and 0.22 on the
shell, 1 × 10<sup>−3</sup> /ft and 0.03 on the core, 2 × 10<sup>−4</sup> /ft and
0.15 on the foundation.

### The pool schedule and the stage times

Switch to **Seepage** and click **Transient** in the Inputs dock.

![The schedule, with both stage times filled](images/combo02_studio_run_transient.png)

**Series names.** Type `pool` over the first column's default name `t1`. A
boundary becomes time-varying by having this exact name typed into its value
cell, which is the last step of this section.

**The breakpoints.** Three rows: the reservoir holds at full pool for five days,
then falls to the tailwater datum over the following 45.

| time | pool |
|---:|---:|
| 0 | 160 |
| 5 | 160 |
| 50 | 100 |

**Duration (day)** `1000`, **Save interval (day)** `200`, and **Extra save
times** `5`, `35`, `50`, `80`, `150`, `300`. The duration carries the run well
past the drawdown so the dam's recovery toward its new steady state is in the
answer too; the extra times put frames where the regular grid would miss them.

**Stage 1 time (day)** `0` and **Stage 2 time (day)** `50`. These two fields are
what a drawdown run reads, and SEEP-3 leaves them deliberately blank because a
seepage page has no use for them. Stage 1 at t = 0 is the full-pool state the
march starts from; stage 2 at t = 50 is the instant the pool reaches the
tailwater datum. Both are forced into the saved-frame schedule, so each is a
*computed* frame rather than a blend of two neighbors. Set both or neither, and
stage 1 must come before stage 2.

The preview draws them as red dashed reference lines across the schedule, which
is the quickest check that they land where they were meant to. Click **OK**.

### Binding the reservoir boundary

Reopen **Seep BC**, select **Head 1** on **Set 1**, and change two fields. Set
**Type:** to `reservoir` — a boundary that holds a node only while that node is
at or below the water level, and releases any node the falling pool has left
standing above it. Then clear **Head value (ft)** and type `pool` in its place. A
value cell holding the name of a series is driven by that series instead of by a
constant, and that is the whole of what makes a boundary time-varying.

Set 2 is untouched by any of this, and cannot follow the series even if asked:
it is the constant-steady set.

Nothing about the numbers already computed changes. A steady solve reads a
series-bound value at t = 0, so Set 1 with `pool` bound to it solves at 160 and
returns the same 1.9566 ft³/day per ft as before, and says so in its log —
*"Series 'pool' read at t = 0 for the steady solve: reservoir value 160"*.

### Marching it

Click **Run → Run Seep…**. The dialog has grown a **Run type** selector, which
appears only on a file carrying a schedule; set it to **Transient
(time-dependent)** and click **Run**. **Convergence tol** and **Max iterations**
gray out — they belong to the steady solve, and the march sets its own step size
from how fast the field is moving.

The march solves its initial condition first — the same unconfined iteration as
the steady run — and then prints a line per saved frame:

```text
  t=5: frame saved (steps so far=6, mass-balance closure=1.39e-05)
  t=35: frame saved (steps so far=374, mass-balance closure=4.34e-02)
  t=50: frame saved (steps so far=594, mass-balance closure=4.62e-02)
  t=80: frame saved (steps so far=958, mass-balance closure=6.11e-02)
  t=150: frame saved (steps so far=1376, mass-balance closure=2.25e-02)
  t=200: frame saved (steps so far=1536, mass-balance closure=1.19e-02)
  t=300: frame saved (steps so far=1720, mass-balance closure=6.19e-02)
  t=400: frame saved (steps so far=1815, mass-balance closure=3.49e-02)
  t=600: frame saved (steps so far=1910, mass-balance closure=9.33e-03)
  t=800: frame saved (steps so far=2019, mass-balance closure=2.52e-02)
  t=1000: frame saved (steps so far=2062, mass-balance closure=3.12e-02)
```

This is the long run of the page — a few minutes, where every other run here has
been seconds. **Twelve frames** come out of it, at t = 0, 5, 35, 50, 80, 150,
200, 300, 400, 600, 800 and 1000: the union of the save-interval grid, the extra
save times, the schedule's own breakpoints, and the two stage times. **2,062
steps** sit behind them, and where the solver spent them is itself a reading:
**594 of them, 29%, are inside the first 50 days**, which are 5% of the run's
duration. The field moves while the pool is falling and barely at all afterward.

The two frames the drawdown will read are the first and the fourth:

![The stage 1 and stage 2 frames, on one color scale](images/combo02_frames.png){width=1000}

At **stage 1** the reservoir is full and the field is the steady one solved
above. At **stage 2** the pool has reached the tailwater datum — the upstream
water level marker has dropped 60 ft — the upstream shell has drained with it,
and the foundation under the upstream slope has begun to unload. The core has
not. It still holds a pocket of head near 150 ft, in the middle of the section,
exactly where the slip circle crosses it. That pocket is the whole difference
between this route and the two steady solves, where the same region reads 100.

### The third answer

Switch to **LEM** and run **Run → Run LEM…** once more — Spencer, Single surface,
Rapid drawdown ticked. The form has grown the group the schedule adds:

![Run LEM on the model with a schedule, showing the stage times](images/combo02_studio_run_lem_staged.png)

**Stage 1 time** and **Stage 2 time** hold 0 and 50, read from the tseep sheet
where the Transient editor stored them, and editing them here changes the model
the same way editing them there does. They replace the single seepage-time
selector an ordinary run on a transient model gets, because a drawdown reads two
instants rather than one; a stage time the solved march never saved re-marches
with it added to the save schedule.

Click **Run**. The two named frames come out of the solved march and their pore
pressures go to the three stages in memory — no `_seep.csv` files are written or
read, and staged frames take precedence over any that are sitting beside the
workbook.

| Stage | | FS |
|---|---|---:|
| 1 | full pool, drained | 2.0542 |
| 2 | drawn down, undrained core | **1.3137** |
| 3 | drawn down, drained | not required |

![The three-stage answer from the transient frames](images/combo02_solution_transient.png){width=1000}

**The rapid-drawdown factor of safety is 1.314.** Stage 1 is identical to the
two-steady run's 2.0542, because the frame at t = 0 *is* that steady solution.
Everything that separates the two answers is at stage 2, where the peak pore
pressure on the slice bases is **3394 psf** against the two-steady route's 1559 —
the core's undissipated head, read straight off the frame.

---

## The three answers, and what brackets them

All three runs used the same dam, the same circle, the same 42 slices, the same
strengths and the same Spencer solver. Only the statement about the water
changed:

| Where the two water states came from | Stage 1 | Stage 2 | Stage 3 | **Rapid drawdown FS** |
|---|---:|---:|---:|---:|
| Piezometric Lines 1 and 2 | 1.9354 | 1.3010 | not required | **1.301** |
| Two steady seepage solutions | 2.0542 | 1.5640 | not required | **1.564** |
| Two frames of a transient march | 2.0542 | 1.3137 | not required | **1.314** |

The spread is 20% of the lowest answer, and all of it is stage 2. The same
circle solved as an ordinary drained problem at each end of the drawdown puts
that spread in context:

| Run | Core strength | Pore pressures | FS |
|---|---|---|---:|
| Full pool | drained | full-pool field | 2.054 |
| Drawn down, long term | drained | re-equilibrated field | 1.706 |
| Drawn down, long term | **undrained** | re-equilibrated field | **1.564** |
| Drawn down, day 50 | drained | transient frame | 1.342 |
| Drawn down, day 50 | **undrained** | transient frame | **1.314** |

The two bold rows are the two-steady and transient drawdown answers; the other
three are ordinary drained or ordinary undrained runs on the same circle, and
between them they separate the drawdown into its parts.

**The load removal alone costs 0.35.** Taking the reservoir off a dam whose
pore pressures have fully caught up moves the slope from 2.054 to 1.706, at
drained strength throughout. That is the long-term condition the dam settles into
months after the drawdown.

**The core's retained water costs a further 0.36.** Holding the day-50 pore
pressures instead of the re-equilibrated ones, still at drained strength, moves
1.706 to 1.342. Nothing about the strengths differs between those two runs.

**The undrained strength costs 0.14 — or 0.03.** On the re-equilibrated field
it moves 1.706 to 1.564; on the day-50 field it moves 1.342 to only 1.314. The
same envelope, and five times the effect on one field as on the other, because
the retained pore pressure has already cut the effective stress the *drained*
strength is computed from, bringing the two strengths close together. Where the
water is decides how much the strength choice changes, which is why the three
sources cannot be compared without saying which field each was read on.

---

## The governing-stage flip

Stage 2 governed all three runs above, and stage 3 never had a lower strength to
substitute. That is a property of this core's envelopes, not of the method, and
it moves.

Hold the transient route fixed — same fields, same circle, same everything — and
sweep the core's undrained cohesion intercept $d$ upward from its 250 psf, with
$\psi$ held at 14°. Raising $d$ raises the undrained strength stage 2 assigns,
which raises the stage-2 factor of safety. It does nothing to the drained
strength. So each step up makes it likelier that on some slice the drained
strength is the lower of the two, which is exactly the test stage 3 applies.

| d (psf) | Stage 2 | Stage 3 | Core slices on drained strength | Reported FS | Governs |
|---:|---:|---:|:---:|---:|:---:|
| 250 | 1.3137 | 1.3137 | 0 of 14 | 1.3137 | stage 2 |
| 300 | 1.3194 | 1.3194 | 0 of 14 | 1.3194 | stage 2 |
| 325 | 1.3223 | 1.3223 | 0 of 14 | 1.3223 | stage 2 |
| 350 | 1.3252 | 1.3251 | 3 of 14 | 1.3251 | stage 3 |
| 375 | 1.3280 | 1.3273 | 5 of 14 | 1.3273 | stage 3 |
| 400 | 1.3309 | 1.3291 | 6 of 14 | 1.3291 | stage 3 |
| 450 | 1.3367 | 1.3320 | 8 of 14 | 1.3320 | stage 3 |
| 500 | 1.3425 | 1.3341 | 10 of 14 | 1.3341 | stage 3 |
| 600 | 1.3540 | 1.3368 | 12 of 14 | 1.3368 | stage 3 |
| 700 | 1.3656 | 1.3385 | 12 of 14 | 1.3385 | stage 3 |
| 800 | 1.3771 | 1.3396 | 13 of 14 | 1.3396 | stage 3 |
| 1000 | 1.4004 | 1.3413 | 13 of 14 | 1.3413 | stage 3 |
| 1200 | 1.4236 | 1.3422 | 14 of 14 | 1.3422 | stage 3 |

![Stage 2, stage 3 and the reported answer against the core's Kc = 1 intercept](images/combo02_sweep.png){width=1000}

**The handover is at d = 341 psf**, 91 psf above the value on the file. Below it
every one of the 14 core slices is stronger undrained than drained, stage 3 is
not run at all, and the two curves are the same line. Above it the crossover
spreads slice by slice — 3 slices at 350 psf, 6 at 400, 10 at 500, all 14 by
1200 — and the two curves separate, with the reported answer following the lower
one.

The shapes of the two curves are the reason the handover does not reverse.
Stage 2 climbs in a straight line, because $d$ enters the undrained strength
directly. Stage 3 flattens, because it can only ever *substitute* the drained
strength, and once a slice has taken the drained value more $d$ changes nothing
on it. Its plateau at **1.3422** is the answer with the whole core drained — and
it is exactly the ordinary drained factor of safety on the same circle at the
same day-50 pore pressures — the fourth row of the bracket table. Stage 3 cannot
buy the slope anything past that, however strong the undrained envelope is
declared to be.

The engineering reading is that the $K_c = 1$ envelope is not a safety margin.
Overstating it does not make the drawdown answer better; past 341 psf it makes no
difference at all beyond a hundredth or two, because the drained strength has
become the limit. What the envelope decides is *which* mechanism is being
checked, and a drawdown run that reports stage 3 as its governing stage is
reporting that the undrained strengths never mattered.

---

## Which source to use

Three answers on one dam, and they are not equally good.

**The transient frames are the most faithful statement of a drawdown.** They
carry the actual lowering rate and the actual elapsed time, so the stage-2 field
is how far dissipation has genuinely progressed by day 50 rather than an
assumption about it. Everything the time factor $T$ estimates at the top of this
page is resolved directly by the solve. On this dam that moves the factor of
safety by 0.25 against the two-steady route, and all of it is the core's
retained head.

**Two steady solutions state a limiting case, and on this model it is the
optimistic one.** Set 2 is the dam after the pore pressures have fully caught up
with the pool — the long-term condition, arrived at instantaneously. That is a
legitimate thing to compute and it is the classic hand-calculation route, but on
a section whose whole difficulty is a core that does not drain, assuming the core
has drained answers a different question. Here it reads **1.564** where the
transient reads **1.314**, 19% high. Where the two-steady route is right is on a
section with no low-conductivity zone to hold water — where the re-equilibrated
field is very nearly what exists a day after the drawdown.

**The piezometric pair is a judgment, and its answer is only as good as
Line 2.** On this model it landed at **1.301**, within 1% of the transient
answer, because Line 2 was drawn from the transient field and the core's held
head was put into it deliberately. Drawn optimistically — flat at the pool level
across the whole section, say — it would have reproduced the two-steady answer or
worse. The pair's real use is where there is no seepage model to solve: an
existing dam with piezometers in it, a preliminary check, or a section that
cannot be meshed. It should not be preferred to a solution when one can be had.

None of which makes the drawdown check itself a matter of taste. The procedure
being run is the same in all three cases, and it is verified against published
answers on four dams:
[VP96](../verification/rocscience.md#vp96) is the USACE EM 1110-2-1902 Appendix G
example (1.434 against the manual's 1.44),
[VP97](../verification/rocscience.md#vp97) is Pilarcitos Dam, which failed in
drawdown (1.044 against 1.05),
[VP98](../verification/rocscience.md#vp98) is the Walter Bouldin Dam failure
(1.046 against 1.04), and
[VP99](../verification/rocscience.md#vp99) is the pumped-storage dam of Duncan,
Wright and Wong's own paper (1.527 against 1.56). All four supply their two water
states as piezometric pairs, which is how their sources published them.

---

## Conclusion

This tutorial covered:

- The three-stage rapid-drawdown procedure and what each stage is for —
  consolidation stresses at full pool, undrained strengths under the drawn-down
  water, and a drained re-check that can only lower the answer — with the
  textbook's own hand calculation, in which two planes in one soil are governed
  by different stages.
- The three inputs the procedure needs: a $d$ / $\psi$ pair on every material
  that does not drain during the drawdown, a second statement of where the water
  is, and **Water loads** on `auto` so the load on the face and the pressure
  inside the dam come from the same place at both stages.
- Three ways to supply the two water states on one dam and one circle — a
  sketched piezometric pair, two steady seepage solutions from two boundary sets,
  and two frames of a transient march — reading 1.301, 1.564 and 1.314.
- The bracket those sit in, and the parts it separates the drawdown into: 0.35
  for the load removal alone (2.054 to 1.706), a further 0.36 for the water the
  core has not released by day 50 (1.706 to 1.342), and 0.14 or 0.03 for the
  undrained strength depending on which of the two fields it is applied to.
- The governing stage as a property of the material, not the method: raising the
  core's undrained intercept past 341 psf hands control from stage 2 to stage 3,
  and past that point a stronger undrained envelope buys nothing, because the
  drained strength has become the limit at 1.342.

**Where to go next:** [Rapid drawdown analysis](../lem/rapid.md) carries the
equations, the time-factor rubric and the interpolation between the two envelopes
in full. [SEEP-3](seep03_reservoir_drawdown.md) builds a transient drawdown model
from nothing and reads it as frames, histories and a water budget;
[SEEP-2](seep02_johnson_dam.md) builds this dam;
[COMBO-1](combo01_seepage_stability.md) runs it under a standing pool. The
[tutorials index](index.md) lists the series.
