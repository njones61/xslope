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

</div>
