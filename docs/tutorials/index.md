---
title: "XSLOPE Tutorials"
description: "Guided, self-contained slope stability problems built from scratch in XSLOPE — each one three ways: the Excel template, Studio's editors, and the AI assistant."
---

# Tutorials

Each tutorial here is one problem, built from scratch and solved. Each is
self-contained: it states its own geometry and materials, and
it ends with the completed model as a file you can open beside your own.

Every tutorial is written three times over — build it with the **AI assistant**,
build the **Excel input file**, or build it in **Studio's editors** — and the three
paths join at the run. Pick whichever suits you and skip the other two; all three
produce the same file.

If you have not used XSLOPE before, read
[Building Models Three Ways](building_models.md) first. It covers what every
tutorial below assumes — the template, the editors, the assistant, what **Save As**
writes, and how to send a whole project to somebody else.

<div class="tut-band" markdown>

## Getting started

<span class="tut-band-sub">1 guide</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/t0_studio_window.png)](building_models.md)
<span class="tut-num">0</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">Start here · All analyses · 10 min read</p>

### [Building Models Three Ways](building_models.md)

The three build paths, the template's worksheets, what **Save As** writes, and `.xslz` packages for sending a project.

<p class="tut-tags"><span class="tut-tag">Excel template</span><span class="tut-tag">Studio editors</span><span class="tut-tag">AI assistant</span><span class="tut-tag">Save As</span><span class="tut-tag">project packages</span></p>

</div>
</div>

<div class="tut-band" markdown>

## Limit equilibrium

<span class="tut-band-sub">13 tutorials</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem01_solution_cracked.png)](lem01_simple_embankment.md)
<span class="tut-num">1</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-1 · Limit equilibrium · 10–30 min</p>

### [Simple Embankment](lem01_simple_embankment.md)

The smallest complete model — one material, one profile line, a rigid base — searched for its critical circle with Spencer's method, then read past the factor of safety to the crest tension its warnings expose.

<p class="tut-tags"><span class="tut-tag">profile lines</span><span class="tut-tag">one material</span><span class="tut-tag">Mohr-Coulomb</span><span class="tut-tag">starting circles</span><span class="tut-tag">circular search</span><span class="tut-tag">tension crack</span><span class="tut-tag">water in crack</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem02_solution_load.png)](lem02_loads_on_the_crest.md)
<span class="tut-num">2</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-2 · Limit equilibrium · 10–25 min</p>

### [Loads on the Crest](lem02_loads_on_the_crest.md)

A surcharge added to LEM-1's embankment — spread over the crest, gathered onto a point, pushed normal to the ground or straight down, and shaken — then a sweep for the cohesion that would carry it at a target factor of safety.

<p class="tut-tags"><span class="tut-tag">distributed loads</span><span class="tut-tag">line loads</span><span class="tut-tag">load direction</span><span class="tut-tag">seismic coefficient</span><span class="tut-tag">design mode</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem03_solution.png)](lem03_layered_slope.md)
<span class="tut-num">3</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-3 · Limit equilibrium · 10–20 min</p>

### [A Layered Slope](lem03_layered_slope.md)

An embankment on a foundation layer over rigid rock — a profile line at the material boundary, a starting circle at the base of each layer, and a search result read at depth: why the critical surface settles on the contact, and what moves it to the rock.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">profile lines</span><span class="tut-tag">per-layer starting circles</span><span class="tut-tag">generated circles</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem04_solution_wet.png)](lem04_water_in_the_slope.md)
<span class="tut-num">4</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-4 · Limit equilibrium · 10–25 min</p>

### [Water in the Slope](lem04_water_in_the_slope.md)

A three-layer slope on a soft foundation clay, searched for its critical circle and then measured on it — the same surface with the pore pressures on and off, u on every slice base, the γ/γ_sat unit-weight split, and r<sub>u</sub>.

<p class="tut-tags"><span class="tut-tag">three materials</span><span class="tut-tag">piezometric line</span><span class="tut-tag">effective stress</span><span class="tut-tag">saturated unit weight</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem05_solution.png)](lem05_weak_layer_noncircular.md)
<span class="tut-num">5</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-5 · Limit equilibrium · 30–60 min</p>

### [A Weak Layer, Non-Circular](lem05_weak_layer_noncircular.md)

A sand slope over a 2 ft seam of soft clay, with the failure surface entered as a table of vertices instead of a circle and then searched — what the Movement column does to each vertex, the methods a polyline rules out, and how steep an end ramp a search will start from.

<p class="tut-tags"><span class="tut-tag">four materials</span><span class="tut-tag">piezometric line</span><span class="tut-tag">non-circular surface</span><span class="tut-tag">non-circular search</span><span class="tut-tag">Movement options</span><span class="tut-tag">weak-zone generator</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem06_solution.png)](lem06_polygon_geometry.md)
<span class="tut-num">6</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-6 · Limit equilibrium · 10–20 min</p>

### [Polygon Geometry](lem06_polygon_geometry.md)

An embankment on bedrock that dips across the section, entered as closed material-zone polygons instead of profile lines and with no maximum depth — then searched, and read for what the base does to it: the circles it refuses, the option that truncates one against it, and the soil change that puts the critical surface on it.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">polygons</span><span class="tut-tag">dipping base</span><span class="tut-tag">composite surfaces</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem07_baker_pow.png)](lem07_strength_envelopes.md)
<span class="tut-num">7</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-7 · Limit equilibrium · 30–40 min</p>

### [Strength Options Beyond Mohr-Coulomb](lem07_strength_envelopes.md)

Two slopes whose strength is not a pair of numbers — a compacted clay fitted once as a curved power envelope and once as the straight Mohr-Coulomb line through the same triaxial data, and a layered clay whose undrained strength grows with depth, flattened to a constant to see how the search responds when strength no longer increases with depth.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">three materials</span><span class="tut-tag">power-curve envelope</span><span class="tut-tag">Mohr-Coulomb</span><span class="tut-tag">strength with depth</span><span class="tut-tag">undrained strength</span><span class="tut-tag">starting circles</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem08_solution.png)](lem08_reinforced_slope.md)
<span class="tut-num">8</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-8 · Limit equilibrium · 15–35 min</p>

### [A Reinforced Slope (Geogrids)](lem08_reinforced_slope.md)

A 24 ft sand fill held at 1.25:1 by six geogrid layers — lines carrying a tensile capacity, a pullout length at each end and a support-type preset — searched, then read against the same slope with the lines removed, against the force each crossing mobilizes, and against the one layer the critical surface never touches.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">distributed load</span><span class="tut-tag">reinforcement lines</span><span class="tut-tag">capacity envelope</span><span class="tut-tag">pullout length</span><span class="tut-tag">support types</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem09_solution.png)](lem09_tieback_wall.md)
<span class="tut-num">9</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-9 · Limit equilibrium · 10–25 min</p>

### [A Tieback Wall](lem09_tieback_wall.md)

A 30 ft soldier-pile wall held by two rows of grouted tiebacks — discrete anchors entered as axial reinforcement with a bearing plate at one end and a grout bond length at the other, plus the pile that carries their heads — searched for its critical wedge, then read against the wedge a reference manual specifies, against the same wall with the anchors removed, and against the two ways an anchor force can be applied.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">reinforcement lines</span><span class="tut-tag">support types</span><span class="tut-tag">anchors</span><span class="tut-tag">piles</span><span class="tut-tag">non-circular search</span><span class="tut-tag">active vs passive</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem10_solution_deep.png)](lem10_global_minimum.md)
<span class="tut-num">10</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-10 · Limit equilibrium · 20–25 min</p>

### [Finding the Global Minimum](lem10_global_minimum.md)

Two slopes with competing failure mechanisms: an embankment on soft clay where moving the starting circle moves the answer, and the James Bay dyke, where a single credible seed reads 23% high — resolved with a minimum slip depth and grid seeding.

<p class="tut-tags"><span class="tut-tag">starting circles</span><span class="tut-tag">circular search</span><span class="tut-tag">grid seeding</span><span class="tut-tag">minimum slip depth</span><span class="tut-tag">James Bay dyke</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem11_mc.png)](lem11_reliability.md)
<span class="tut-num">11</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-11 · Limit equilibrium · ~30 min</p>

### [Reliability](lem11_reliability.md)

A submerged clay slope whose strength and unit weight carry standard deviations — the factor of safety of 1.354 turned into a reliability index and a probability of failure two ways, by the Taylor series over 1 + 2N searches and by a 10,000-realization Monte Carlo campaign, then the dominant uncertainty measured and halved.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">undrained strength</span><span class="tut-tag">distributed load</span><span class="tut-tag">standard deviations</span><span class="tut-tag">Taylor series (TSPM)</span><span class="tut-tag">Monte Carlo</span><span class="tut-tag">variance Pareto</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem12_solution.png)](lem12_piles.md)
<span class="tut-num">12</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-12 · Limit equilibrium · ~25 min</p>

### [Piles](lem12_piles.md)

A 20 ft clay slope stabilized by two rows of drilled shafts whose force column is left blank, so Ito & Matsui computes it from the diameter and the spacing at every trial surface — read in the report the run writes, capped by the shaft's own bending capacity, swept across the spacing band, compared against the same force stated outright, and checked against the shallow surface that slides over the pile row.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">piles</span><span class="tut-tag">Ito &amp; Matsui</span><span class="tut-tag">pile spacing</span><span class="tut-tag">structural capacity</span><span class="tut-tag">specified pile force</span><span class="tut-tag">circular search</span><span class="tut-tag">grid seeding</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/lem13_shear_strain.png)](lem13_rock_slope.md)
<span class="tut-num">13</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">LEM-13 · Limit equilibrium, finite element (strength reduction) · ~60 min</p>

### [A Rock Slope (Hoek-Brown)](lem13_rock_slope.md)

A 10 m rock cut whose strength is described by four field observations rather than a cohesion and a friction angle — the intact strength, the Geological Strength Index, the intact constant and the blast disturbance — solved by Spencer's method and then by strength reduction on the same file, read once more with its intact strength typed in the megapascals the literature quotes, and then swept across the Geological Strength Index and the blast-disturbance factor to find which of the four field inputs controls the answer.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">Hoek-Brown</span><span class="tut-tag">GSI</span><span class="tut-tag">disturbance factor</span><span class="tut-tag">instantaneous tangent</span><span class="tut-tag">circular search</span><span class="tut-tag">strength reduction</span><span class="tut-tag">quadratic triangles</span><span class="tut-tag">parametric study</span><span class="tut-tag">design sweep</span><span class="tut-tag">model checks</span></p>

</div>
</div>

<div class="tut-band" markdown>

## Seepage

<span class="tut-band-sub">4 tutorials</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/seep01_solution_filled.png)](seep01_sheetpile.md)
<span class="tut-num">1</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">SEEP-1 · Seepage · 10–30 min</p>

### [Confined Seepage Under a Sheetpile](seep01_sheetpile.md)

A sheetpile driven 3 m into a 10 m silt foundation behind an impervious clay blanket, built from scratch and solved for the discharge under it — what makes a problem confined, the three boundary condition types and the no-flow default that models the blanket, meshing and element order, the flow net set up so its channels and head drops come out whole, a mesh study on a problem where the discharge keeps changing with refinement, and a conductivity sweep that is exactly proportional.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">confined flow</span><span class="tut-tag">specified head</span><span class="tut-tag">no-flow boundaries</span><span class="tut-tag">mesh generation</span><span class="tut-tag">element types</span><span class="tut-tag">flow net</span><span class="tut-tag">feature refinement</span><span class="tut-tag">parametric study</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/seep02_solution.png)](seep02_johnson_dam.md)
<span class="tut-num">2</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">SEEP-2 · Seepage · ~40 min</p>

### [Unconfined Seepage Through a Zoned Dam](seep02_johnson_dam.md)

An 80 ft earth dam with a shell, a keyed clay core and 60 ft of water behind it, opened and explored — the seepage face whose wet extent is an output rather than an input, where the head drops and how much of the flow goes under the core and above the phreatic surface, the base material a zoned flow net has to be scaled to, the three unsaturated conductivity models run against each other and the one parameter that accounts for the difference between them, and a run that stops without converging.

<p class="tut-tags"><span class="tut-tag">three materials</span><span class="tut-tag">unconfined flow</span><span class="tut-tag">seepage face</span><span class="tut-tag">phreatic surface</span><span class="tut-tag">unsaturated models</span><span class="tut-tag">relative conductivity</span><span class="tut-tag">flow net base material</span><span class="tut-tag">convergence</span><span class="tut-tag">underseepage</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/seep03_frames.png)](seep03_reservoir_drawdown.md)
<span class="tut-num">3</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">SEEP-3 · Seepage · ~30 min</p>

### [Transient Seepage: Reservoir Drawdown](seep03_reservoir_drawdown.md)

A 22 m dam with a granular shell and a clay core, whose reservoir is lowered to the tailwater elevation over 45 days — the storage properties that decide how much water a soil gives up as its head falls, the full-pool steady solution the transient run starts from, a pool schedule entered as a named time series and bound to the upstream boundary, and the three ways a transient answer is read: frames on a play bar, a head history where the core's head passes above the shell's, and an outflow ten times the steady discharge drawn entirely from storage.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">transient seepage</span><span class="tut-tag">specific storage</span><span class="tut-tag">specific yield</span><span class="tut-tag">time series</span><span class="tut-tag">reservoir boundary</span><span class="tut-tag">initial condition</span><span class="tut-tag">saved frames</span><span class="tut-tag">head history</span><span class="tut-tag">mass balance</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/seep04_wet.png)](seep04_dam_infiltration.md)
<span class="tut-num">4</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">SEEP-4 · Seepage · ~35 min</p>

### [Infiltration and Flux Boundaries](seep04_dam_infiltration.md)

A 12 m earth dam with a reservoir, a toe drain and steady rain falling on its exposed surface, solved twice with the rain as the only difference — the specified flux that carries rainfall, the projection that turns a vertical rain rate into the normal flux a boundary takes, why a flux boundary is drawn on the section's corners rather than on a mesh's nodes, what happens at the nodes where it runs into the reservoir and the drain, and a water budget where the discharge rises 75% while the reservoir's share falls 36% — then a sweep of six rain rates that finds the ceiling past which the soil can no longer take what the boundary prescribes, and finally a transient run in which the same three flux boundaries are driven from a time series so the rain starts, holds and stops.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">steady seepage</span><span class="tut-tag">specified flux</span><span class="tut-tag">infiltration</span><span class="tut-tag">normal flux</span><span class="tut-tag">specified head</span><span class="tut-tag">exit face</span><span class="tut-tag">toe drain</span><span class="tut-tag">van Genuchten</span><span class="tut-tag">phreatic surface</span><span class="tut-tag">water budget</span><span class="tut-tag">transient seepage</span><span class="tut-tag">time series</span><span class="tut-tag">storage</span></p>

</div>
</div>

<div class="tut-band" markdown>

## Finite element

<span class="tut-band-sub">5 tutorials</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/fem01_shear_strain.png)](fem01_strength_reduction.md)
<span class="tut-num">1</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">FEM-1 · Finite element · ~25 min</p>

### [Strength Reduction Basics](fem01_strength_reduction.md)

A 50 ft embankment on a rigid base solved twice — once with Spencer's method for the reference answer, once by meshing the slope and weakening its soil until it can no longer stand — the two elastic properties the run requires and the measurement showing E scales every displacement while changing the factor of safety not at all, why a stability mesh has to be quadratic, the bracket the search validates and bisects with the iteration count climbing toward failure, the per-trial budget that decides which near-critical trials are allowed to finish, and the mechanism read out of shear strain, a deformed mesh and displacement vectors against the converged state on one color scale.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">Mohr-Coulomb</span><span class="tut-tag">strength reduction</span><span class="tut-tag">Young's modulus</span><span class="tut-tag">Poisson's ratio</span><span class="tut-tag">quadratic triangles</span><span class="tut-tag">side boundary conditions</span><span class="tut-tag">tensile cutoff</span><span class="tut-tag">bisection bracket</span><span class="tut-tag">iteration budget</span><span class="tut-tag">shear strain</span><span class="tut-tag">deformed mesh</span><span class="tut-tag">displacement vectors</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/fem02_shear_strain_pr.png)](fem02_reinforcement.md)
<span class="tut-num">2</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">FEM-2 · Finite element · ~30 min</p>

### [Reinforcement: LEM vs FEM](fem02_reinforcement.md)

Six layers of geogrid in a 24 ft sand fill solved by Spencer's method and then by strength reduction — the two limits every reinforcement model has, bond and rupture, and what each engine does with them; the axial stiffness the finite element run requires; the per-trial iteration budget without which both runs return the same number; the elastic-perfectly-plastic run against the peak-residual one, read out of shear strain and each layer's own force profile; and two measurements of how much the modeling choices change the answer — the soil-to-bar stiffness ratio, and a residual capacity that lands on three values rather than a curve.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">reinforcement lines</span><span class="tut-tag">capacity envelope</span><span class="tut-tag">pullout length</span><span class="tut-tag">bar elements</span><span class="tut-tag">axial stiffness</span><span class="tut-tag">elastic-perfectly-plastic</span><span class="tut-tag">residual capacity</span><span class="tut-tag">strength reduction</span><span class="tut-tag">quadratic triangles</span><span class="tut-tag">thin zones</span><span class="tut-tag">iteration budget</span><span class="tut-tag">1D details</span><span class="tut-tag">shear strain</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/fem03_fem_blocks_failure.png)](fem03_block_wall_joints.md)
<span class="tut-num">3</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">FEM-3 · Finite element · ~60 min (2 h with the long run)</p>

### [A Block Wall on Slip Joints](fem03_block_wall_joints.md)

A 3.6 m segmental block wall whose every contact is a slip joint — under the base of the block column, on its back face against the fill, and between each pair of courses — so the wall can slide and tip the way a stack of blocks does instead of standing as one solid; the two settings a jointed run needs and what leaving them alone costs; three geogrid layers tied into the facing and what they add; and then two small embankments, each run twice with one cell changed, to show when a geosynthetic should be a bonded bar and when a slip surface — a base geotextile the failure surface cuts across, where a bonded bar and a slip surface agree, and a smooth liner the fill slides out along, where they do not.

<p class="tut-tags"><span class="tut-tag">four materials</span><span class="tut-tag">slip joints</span><span class="tut-tag">joints worksheet</span><span class="tut-tag">interface element</span><span class="tut-tag">block wall</span><span class="tut-tag">elastic blocks</span><span class="tut-tag">geogrid</span><span class="tut-tag">end ties</span><span class="tut-tag">bonded bar</span><span class="tut-tag">jointed sheet</span><span class="tut-tag">hybrid criterion</span><span class="tut-tag">sweep budget</span><span class="tut-tag">local element size</span><span class="tut-tag">joint slip</span><span class="tut-tag">deformed blocks</span><span class="tut-tag">1D details</span><span class="tut-tag">model checks</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/fem04_fem_shear_piles.png)](fem04_piles.md)
<span class="tut-num">4</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">FEM-4 · Finite element · ~60 min</p>

### [Piles: LEM vs FEM](fem04_piles.md)

Two rows of drilled shafts and a sheet pile wall on the same slope, each put through both engines — why a plane-strain finite element model has no gap between piles and turns a discrete row into a continuous wall, seen in the mechanism it forces; a spacing sweep in which the limit equilibrium answer swings 36% and changes which surface governs while the strength reduction answer does not move; the effect of fixing the pile tip, 0.047; the one benchmark with a published three-dimensional answer; and then a case suited to the finite element engine — a wall entered as a continuous member at spacing 1, run with its toe free and then fixed, and read for the moment, shear, deflection and soil reaction no limit equilibrium analysis can report.

<p class="tut-tags"><span class="tut-tag">piles</span><span class="tut-tag">Ito &amp; Matsui</span><span class="tut-tag">pile spacing</span><span class="tut-tag">head and tip fixity</span><span class="tut-tag">beam elements</span><span class="tut-tag">smeared stiffness</span><span class="tut-tag">plane strain</span><span class="tut-tag">sheet pile wall</span><span class="tut-tag">continuous member</span><span class="tut-tag">strength reduction</span><span class="tut-tag">quadratic triangles</span><span class="tut-tag">1D element size</span><span class="tut-tag">shear strain</span><span class="tut-tag">moment and shear profiles</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/fem05_fem_blocks_topple_colored.png)](fem05_rock_slope_joints.md)
<span class="tut-num">5</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">FEM-5 · Finite element · ~30 min</p>

### [A Rock Slope on Its Joints](fem05_rock_slope_joints.md)

A 12 m rock face cut by a bedding plane that daylights at the toe and a release joint behind the crest, with the rock declared elastic so the strength reduction has nothing but the two surfaces to weaken — a factor of safety of 1.199 against the closed-form tan φ / tan β, and a deformed-block figure showing the slab sliding out with the release joint carrying nothing; then the same face cut by a set of joints generated from a dip and a spacing over a joint region instead of typed line by line, where the rock stands as columns that tip rather than slide, and a measurement of what halving the spacing costs.

<p class="tut-tags"><span class="tut-tag">one material</span><span class="tut-tag">elastic rock</span><span class="tut-tag">joints worksheet</span><span class="tut-tag">bedding plane</span><span class="tut-tag">release joint</span><span class="tut-tag">plane failure</span><span class="tut-tag">closed-form check</span><span class="tut-tag">joint region</span><span class="tut-tag">Build network</span><span class="tut-tag">parallel set</span><span class="tut-tag">Voronoi</span><span class="tut-tag">block toppling</span><span class="tut-tag">strength reduction</span><span class="tut-tag">hybrid criterion</span><span class="tut-tag">sweep budget</span><span class="tut-tag">joint slip</span><span class="tut-tag">deformed blocks</span></p>

</div>
</div>

<div class="tut-band" markdown>

## Combined

<span class="tut-band-sub">3 tutorials</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/combo01_seepage.png)](combo01_seepage_stability.md)
<span class="tut-num">1</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">COMBO-1 · Seepage + limit equilibrium + finite element · ~20 min</p>

### [Seepage &rarr; LEM &rarr; FEM](combo01_seepage_stability.md)

The Johnson Reservoir dam run three ways from one file — a seepage solution on a quadratic mesh, then a Spencer search and a strength reduction on that same mesh and the same pore pressures — the element type a shared mesh has to be chosen for, the materials table's `u` column that decides whether the stability engines read the seepage answer at all and what leaving it on `none` costs, and the staleness rule that stops the three results belonging to three different models.

<p class="tut-tags"><span class="tut-tag">three materials</span><span class="tut-tag">unconfined seepage</span><span class="tut-tag">shared mesh</span><span class="tut-tag">quadratic triangles</span><span class="tut-tag">pore-pressure option</span><span class="tut-tag">u = seep</span><span class="tut-tag">Spencer</span><span class="tut-tag">circular search</span><span class="tut-tag">strength reduction</span><span class="tut-tag">shear strain</span><span class="tut-tag">stale results</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/combo02_solution_transient.png)](combo02_rapid_drawdown.md)
<span class="tut-num">2</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">COMBO-2 · Limit equilibrium + seepage · ~35 min</p>

### [Rapid Drawdown](combo02_rapid_drawdown.md)

The same dam with its reservoir lowered 50 ft, from elevation 160 to a residual pool at 110 — the Duncan, Wright and Wong three-stage procedure and its `d` / `ψ` undrained envelope, then its two water states supplied three ways, each searched from the same starting circle: a sketched piezometric pair, two steady seepage solutions from a second boundary set, and two frames of a transient run, reading 1.196, 1.206, and 1.034 — with the drained runs that bracket them, and a sweep of the core's undrained intercept that moves the governing stage from stage 2 to stage 3.

<p class="tut-tags"><span class="tut-tag">three materials</span><span class="tut-tag">rapid drawdown</span><span class="tut-tag">three-stage procedure</span><span class="tut-tag">Kc = 1 envelope</span><span class="tut-tag">d and ψ</span><span class="tut-tag">piezometric lines</span><span class="tut-tag">second boundary set</span><span class="tut-tag">transient seepage</span><span class="tut-tag">stage times</span><span class="tut-tag">automatic water loads</span><span class="tut-tag">Spencer</span><span class="tut-tag">circular search</span><span class="tut-tag">parametric sweep</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/combo03_curve.png)](combo03_fs_vs_time.md)
<span class="tut-num">3</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">COMBO-3 · Seepage + limit equilibrium (drained, and rapid drawdown) · ~35 min</p>

### [Factor of Safety vs Time](combo03_fs_vs_time.md)

Two dams run for stability against every saved frame of a transient drawdown. Part 1 gives SEEP-3's cored earth dam drained strengths — how `u = seep` and Run LEM's **Seepage time** group send one instant of a transient run to a limit equilibrium search, how to place one starting circle per face so a curve spanning a full reservoir and an empty one can report either, and a minimum of 1.331 on day 35 on the upstream face against 1.531 on the downstream face at full pool. Part 2 turns the same sweep into a three-stage rapid drawdown at every instant of COMBO-2's Johnson Reservoir transient run: 1.034 on day 50, which is COMBO-2's own answer on the same circle, against 1.054 for the drained analysis of that instant — with both curves on one figure and the reason the gap is narrowest where the slope is weakest.

<p class="tut-tags"><span class="tut-tag">two materials</span><span class="tut-tag">three materials</span><span class="tut-tag">transient seepage</span><span class="tut-tag">u = seep</span><span class="tut-tag">saved frames</span><span class="tut-tag">seepage time</span><span class="tut-tag">FS vs time</span><span class="tut-tag">rapid drawdown</span><span class="tut-tag">three-stage procedure</span><span class="tut-tag">Kc = 1 envelope</span><span class="tut-tag">d and ψ</span><span class="tut-tag">stage times</span><span class="tut-tag">Spencer</span><span class="tut-tag">circular search</span><span class="tut-tag">starting circles</span><span class="tut-tag">minimum slip depth</span><span class="tut-tag">automatic water loads</span></p>

</div>
</div>

<div class="tut-band" markdown>

## Workflow

<span class="tut-band-sub">3 tutorials</span>

</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/w1_build_from_image_window.png)](w01_ai_assistant.md)
<span class="tut-num">1</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">W-1 · Any · ~30 min</p>

### [The AI Assistant](w01_ai_assistant.md)

Studio's built-in assistant put through eight requests on the layered slope of LEM-3 — setting it up with a provider, a model and a key, and what a turn costs; then building the model from nothing but its drawing, editing its face, its water and its strength, a sweep the Parametric dialog offers and a sweep it does not, material stiffnesses suggested and then run through strength reduction, a question answered out of this documentation and the verification corpus, a broken file diagnosed one input at a time, and the analysis report — with what to check after each one.

<p class="tut-tags"><span class="tut-tag">AI assistant</span><span class="tut-tag">provider setup</span><span class="tut-tag">images</span><span class="tut-tag">model edits</span><span class="tut-tag">parameter studies</span><span class="tut-tag">strength reduction</span><span class="tut-tag">documentation</span><span class="tut-tag">analysis report</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/w02_dxf_wizard.png)](w02_import.md)
<span class="tut-num">2</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">W-2 · Limit equilibrium · ~40 min</p>

### [Importing Models from CAD and Other Programs](w02_import.md)

The three limit-equilibrium importers of Studio's four walked end to end — a layered DXF section of a highway embankment mapped layer by layer through the import wizard, with the seven targets a layer can take and the properties no CAD file can hold typed in afterward (1.184 on a surface tangent to the dense sand); a GeoStudio SLOPE/W verification dam that arrives whole, including one of SLOPE/W's own solved trial circles, so the two programs can be compared on one surface (1.939 against 1.934); and a Slide2 scenario that arrives with no failure surface at all, because Slide2 stores a search rather than a surface (1.372 against 1.360 once a starting circle is added) — with the note list every import reports, read against what each source model defines.

<p class="tut-tags"><span class="tut-tag">DXF import</span><span class="tut-tag">layer mapping</span><span class="tut-tag">material zones</span><span class="tut-tag">piezometric line</span><span class="tut-tag">distributed load</span><span class="tut-tag">failure circles</span><span class="tut-tag">GeoStudio import</span><span class="tut-tag">analysis picker</span><span class="tut-tag">imported trial circle</span><span class="tut-tag">Slide2 import</span><span class="tut-tag">scenario picker</span><span class="tut-tag">import notes</span><span class="tut-tag">Spencer</span><span class="tut-tag">circular search</span></p>

</div>
</div>

<div class="tut-row" markdown>
<div class="tut-thumb" markdown>
[![](images/w03_page_spencer.png)](w03_report.md)
<span class="tut-num">3</span>
</div>
<div class="tut-body" markdown>
<p class="tut-kicker">W-3 · Seepage + limit equilibrium + finite element · ~20 min</p>

### [Generating the Analysis Report](w03_report.md)

A solved model turned into a Word calculation package — COMBO-1's Johnson Reservoir dam opened with its seepage and strength reduction solutions already attached, so the report documents them instead of solving them again; the Generate Report dialog walked group by group, with both stability methods documented in full (Spencer 1.248 from the search we run, Bishop 1.222 from the search the report runs for it); what the document states, from the flow net to each method's own equation with its converged numbers in it; and the same report built a second time on a company letterhead made from the template XSLOPE ships — the style names that have to survive the edit, the document properties that carry the title page into the header, and where a logo has to sit to print on every page.

<p class="tut-tags"><span class="tut-tag">generate report</span><span class="tut-tag">attached solutions</span><span class="tut-tag">method selection</span><span class="tut-tag">title page</span><span class="tut-tag">contents tree</span><span class="tut-tag">flow net</span><span class="tut-tag">calculations</span><span class="tut-tag">slice table</span><span class="tut-tag">company template</span><span class="tut-tag">Word styles</span><span class="tut-tag">letterhead</span></p>

</div>
</div>
