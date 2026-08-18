---
title: "XSLOPE Tutorials"
description: "Guided, self-contained slope stability problems built from scratch in XSLOPE — each one three ways: the Excel template, Studio's editors, and the AI assistant."
---

# Tutorials

A tutorial here is one problem, built from nothing and carried through to a result
you can read. Each is self-contained: it states its own geometry and materials, and
it ends with the completed model as a file you can open beside your own.

Every tutorial is written three times over — build it with the **AI assistant**,
build the **Excel input file**, or build it in **Studio's editors** — and the three
paths meet again at the run. Pick whichever suits you and skip the other two; you are
not choosing a workflow, only an entry point, since all three produce the same file.

If you have not used XSLOPE before, read
[Building Models Three Ways](building_models.md) first. It covers the machinery every
tutorial below assumes — the template, the editors, the assistant, what **Save As**
writes, and how to send a whole project to somebody else.

---

## Getting started

<div class="tut-summary" markdown>

| # | Tutorial | What it covers | Analysis | Features | Time |
|---|---|---|---|---|---|
| 0 | [Building Models Three Ways](building_models.md) | The three build paths, the template's worksheets, what **Save As** writes, and `.xslz` packages for sending a project | All | Excel template, Studio editors, AI assistant, Save As, project packages | 10 min read |

</div>

## Limit equilibrium

<div class="tut-summary" markdown>

| # | Tutorial | What it covers | Analysis | Features | Time |
|---|---|---|---|---|---|
| 1 | [Simple Embankment](lem01_simple_embankment.md) | The smallest complete model — one material, one profile line, a rigid base — searched for its critical circle with Spencer's method, then read past the factor of safety to the crest tension its warnings expose | Limit equilibrium | profile lines, one material, Mohr-Coulomb, starting circles, circular search, tension crack | 10–30 min |
| 2 | [Loads on the Crest](lem02_loads_on_the_crest.md) | A surcharge added to LEM-1's embankment — spread over the crest, gathered onto a point, pushed normal to the ground or straight down, and shaken — then a sweep for the cohesion that would carry it at a target factor of safety | Limit equilibrium | distributed loads, line loads, load direction, seismic coefficient, design mode | 5–20 min |
| 3 | [A Layered Slope](lem03_layered_slope.md) | An embankment on a foundation layer over rigid rock — a profile line at the material boundary, a starting circle at the base of each layer, and a search result read at depth: why the critical surface settles on the contact, and what moves it to the rock | Limit equilibrium | two materials, profile lines, per-layer starting circles, generated circles, circular search | 5–20 min |
| 4 | [Water in the Slope](lem04_water_in_the_slope.md) | A three-layer slope on a soft foundation clay, searched for its critical circle and then measured on it — the same surface with the pore pressures on and off, u on every slice base, the γ/γ_sat unit-weight split, and r<sub>u</sub> | Limit equilibrium | three materials, piezometric line, effective stress, saturated unit weight, circular search | 5–20 min |
| 5 | [A Weak Layer, Non-Circular](lem05_weak_layer_noncircular.md) | A sand slope over a 2 ft seam of soft clay, with the failure surface entered as a table of vertices instead of a circle and then searched — what the Movement column does to each vertex, the methods a polyline rules out, and how steep an end ramp a search will start from | Limit equilibrium | four materials, piezometric line, non-circular surface, non-circular search, Movement options, weak-zone generator | 5–20 min |
| 6 | [Polygon Geometry](lem06_polygon_geometry.md) | An embankment on bedrock that dips across the section, entered as closed material-zone polygons instead of profile lines and with no maximum depth — then searched, and read for what the base does to it: the circles it refuses, the option that truncates one against it, and the soil change that puts the critical surface on it | Limit equilibrium | two materials, polygons, dipping base, composite surfaces, circular search | 5–20 min |
| 7 | [Strength Options Beyond Mohr-Coulomb](lem07_strength_envelopes.md) | Two slopes whose strength is not a pair of numbers — a compacted clay fitted once as a curved power envelope and once as the straight Mohr-Coulomb line through the same triaxial data, and a layered clay whose undrained strength grows with depth, flattened to a constant to see what the search does when depth stops buying strength | Limit equilibrium | one material, three materials, power-curve envelope, Mohr-Coulomb, strength with depth, undrained strength, starting circles, circular search | ~15 min |
| 8 | [A Reinforced Slope (Geogrids)](lem08_reinforced_slope.md) | A 24 ft sand fill held at 1.25:1 by six geogrid layers — lines carrying a tensile capacity, a pullout length at each end and a support-type preset — searched, then read against the same slope with the lines removed, against the force each crossing mobilizes, and against the one layer the critical surface never touches | Limit equilibrium | two materials, distributed load, reinforcement lines, capacity envelope, pullout length, support types, circular search | 5–20 min |
| 9 | [A Tieback Wall](lem09_tieback_wall.md) | A 30 ft soldier-pile wall held by two rows of grouted tiebacks — discrete anchors entered as axial reinforcement with a bearing plate at one end and a grout bond length at the other, plus the pile that carries their heads — searched for its critical wedge, then read against the wedge a reference manual specifies, against the same wall with the anchors removed, and against the two ways an anchor force can be applied | Limit equilibrium | two materials, reinforcement lines, support types, anchors, piles, non-circular search, active vs passive | 5–20 min |
| 10 | [Finding the Global Minimum](lem10_global_minimum.md) | Two slopes with competing failure mechanisms: an embankment on soft clay where moving the starting circle moves the answer, and the James Bay dyke, where a single credible seed reads 23% high — resolved with a minimum slip depth and grid seeding | Limit equilibrium | starting circles, circular search, grid seeding, minimum slip depth, James Bay dyke | ~15 min |
| 11 | [Reliability](lem11_reliability.md) | A submerged clay slope whose strength and unit weight carry standard deviations — the factor of safety of 1.354 turned into a reliability index and a probability of failure two ways, by the Taylor series over 1 + 2N searches and by a 10,000-realization Monte Carlo campaign, then the dominant uncertainty measured and halved | Limit equilibrium | one material, undrained strength, distributed load, standard deviations, Taylor series (TSPM), Monte Carlo, variance Pareto, circular search | ~15 min |
| 12 | [Piles](lem12_piles.md) | A 20 ft clay slope stabilized by two rows of drilled shafts whose force column is left blank, so Ito & Matsui computes it from the diameter and the spacing at every trial surface — read in the report the run writes, capped by the shaft's own bending capacity, swept across the spacing band, compared against the same force stated outright, and checked against the shallow surface that slides over the pile row | Limit equilibrium | one material, piles, Ito &amp; Matsui, pile spacing, structural capacity, specified pile force, circular search, grid seeding | ~15 min |

</div>

## Seepage

<div class="tut-summary" markdown>

| # | Tutorial | What it covers | Analysis | Features | Time |
|---|---|---|---|---|---|
| 1 | [Seepage Under a Sheetpile](seep01_sheetpile.md) | A sheetpile driven 3 m into a 10 m sand foundation behind an impervious clay blanket, built from scratch and solved for the discharge under it — what makes a problem confined, the three boundary condition types and the no-flow default that models the blanket, meshing and element order, the flow net set up so its channels and head drops come out whole, how much of the answer is mesh on a problem where the sequence never settles, and a conductivity sweep that is exactly proportional | Seepage | one material, confined flow, specified head, no-flow boundaries, mesh generation, element types, flow net, feature refinement, parametric study | 25–30 min |
| 2 | [Unconfined Seepage Through a Zoned Dam](seep02_johnson_dam.md) | An 80 ft earth dam with a shell, a keyed clay core and 60 ft of water behind it, opened and explored — the seepage face whose wet extent is an output rather than an input, where the head drops and how much of the flow goes under the core and above the phreatic surface, the base material a zoned flow net has to be scaled to, the three unsaturated conductivity models run against each other and the one parameter the difference between them traces to, and a run that stops without converging | Seepage | three materials, unconfined flow, seepage face, phreatic surface, unsaturated models, relative conductivity, flow net base material, convergence, underseepage | ~25 min |

</div>
