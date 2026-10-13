---
title: "XSLOPE Videos"
description: "Short video walkthroughs of XSLOPE, each with the model files it starts from and ends with."
---

# Videos

The videos are published on the [xslope YouTube channel](https://www.youtube.com/@xslope-org).
Each entry below gives what the video covers, a link to it, and the completed
model it builds. When a video starts from an existing model, that starting file is
listed too, so you can download it and follow along.

## Getting started

### Video 1. A Simple Embankment from Scratch

A 20 ft embankment of undrained clay, built in Studio's editors starting from
**File → New**: the global parameters, one material, one profile line, and a
starting circle from the generator. A Spencer search finds the critical circle,
its warnings show the crest in tension, and an 8 ft tension crack removes the
tension and lowers the factor of safety.

- **Watch** — coming soon
- **Related tutorial** — [LEM-1. Simple Embankment](../tutorials/lem01_simple_embankment.md)

### Video 2. Building a Model Three Ways

The Video 1 embankment built two more ways: in the Excel input template and by
Studio's AI assistant. All three routes produce the same file and the same factor
of safety.

- **Watch** — coming soon
- **Related tutorial** — [Building Models Three Ways](../tutorials/building_models.md)

## Limit equilibrium

### Video 3. Adding Layers

An embankment on a foundation layer over rigid rock, with a profile line at each
material boundary and a starting circle at the base of each layer. The search
result shows why the critical circle settles on the contact between layers.

- **Watch** — coming soon
- **Related tutorial** — [LEM-3. A Layered Slope](../tutorials/lem03_layered_slope.md)

### Video 4. Water in the Slope

A piezometric line added to a layered slope, with saturated unit weights below
it. The same critical surface is solved with the pore pressures on and off to
show what the water costs.

- **Watch** — coming soon
- **Related tutorial** — [LEM-4. Water in the Slope](../tutorials/lem04_water_in_the_slope.md)

### Video 5. Loads and Earthquake

A surcharge on the Video 1 embankment, entered as a distributed load and as a line
load, followed by a seismic coefficient.

- **Watch** — coming soon
- **Related tutorial** — [LEM-2. Loads on the Crest](../tutorials/lem02_loads_on_the_crest.md)

### Video 6. A Non-Circular Surface Through a Weak Layer

A sand slope over a thin seam of soft clay, with the failure surface entered as
a table of vertices and then searched. The Movement column controls how each
vertex is allowed to move during the search.

- **Watch** — coming soon
- **Related tutorial** — [LEM-5. A Weak Layer, Non-Circular](../tutorials/lem05_weak_layer_noncircular.md)

## Seepage

### Video 7. Confined Seepage Under a Sheetpile

A sheetpile in a silt foundation behind an impervious clay blanket, built from a
blank model: the geometry, the three boundary condition types and the no-flow
default, the mesh, the solution, the discharge, and the flow net.

- **Watch** — coming soon
- **Related tutorial** — [SEEP-1. Confined Seepage Under a Sheetpile](../tutorials/seep01_sheetpile.md)

### Video 8. Unconfined Seepage Through a Zoned Dam

An earth dam with a shell and a clay core, built from a blank model, with the
reservoir head on the upstream face and a seepage face downstream. The solution
locates the phreatic surface.

- **Watch** — coming soon
- **Related tutorial** — [SEEP-2. Unconfined Seepage Through a Zoned Dam](../tutorials/seep02_johnson_dam.md)

### Video 9. Seepage Pore Pressures in a Stability Run

The Video 8 dam searched with Spencer's method on the pore pressures from its
seepage solution, by setting `u` to `seep` in the materials table.

- **Watch** — coming soon
- **Related tutorial** — [COMBO-1. Seepage → LEM → FEM](../tutorials/combo01_seepage_stability.md)

## Finite element

### Video 10. Strength Reduction on an Embankment

An embankment solved by the shear strength reduction method: the elastic
properties the run needs, a quadratic mesh, the run, and the deformed mesh and
shear strain at failure, compared with Spencer's method on the same slope.

- **Watch** — coming soon
- **Related tutorial** — [FEM-1. Strength Reduction Basics](../tutorials/fem01_strength_reduction.md)

### Video 11. Seepage, LEM, and FEM on One Mesh

The Video 8 dam run through all three engines on one quadratic mesh: the seepage
solution, a Spencer search on its pore pressures, and a strength reduction run.

- **Watch** — coming soon
- **Related tutorial** — [COMBO-1. Seepage → LEM → FEM](../tutorials/combo01_seepage_stability.md)


<!--
Entry format (one per video, grouped under a series heading):

### Video 1. A Simple Embankment from Scratch

Two or three sentences on what the video builds and what it shows.

- **Watch** — [YouTube](https://youtu.be/VIDEO_ID) · 12 min
- **Starting model** — [name_start.xlsx](files/name_start.xlsx) (only when the video does not start from a blank model)
- **Completed model** — [name.xlsx](files/name.xlsx)
- **Related tutorial** — [LEM-1 Simple Embankment](../tutorials/lem01_simple_embankment.md)
-->
