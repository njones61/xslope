---
title: "XSLOPE Videos"
description: "Short video walkthroughs of XSLOPE, each with the model files it starts from and ends with."
---

# Videos

<div class="vid-hero" markdown>

Short walkthroughs of XSLOPE, each 15 minutes or less, published on the
**[xslope YouTube channel](https://www.youtube.com/@xslope-org)**. Each video
below lists the completed model it builds, and the starting model when it doesn't
begin from a blank one, so you can download the file and follow along.

</div>

<div class="vid-band" markdown>

## Getting started

<span class="vid-band-sub">Videos 1–2</span>

</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
[![](../tutorials/images/lem01_solution_cracked.png)](https://youtu.be/bNy9WkYUvLc)
<span class="vid-num">1</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 1 · Limit equilibrium</p>

### A Simple Embankment from Scratch {#video-1}

A 20 ft embankment of undrained clay, built in Studio's editors starting from
**File → New**: the global parameters, one material, one profile line, and a
starting circle from the generator. A Spencer search finds the critical circle,
its warnings show the crest in tension, and an 8 ft tension crack removes the
tension and lowers the factor of safety.

<p class="vid-meta"><a class="vid-pill vid-watch" href="https://youtu.be/bNy9WkYUvLc">▶ Watch · 18 min</a></p>

**Completed model** — [xslope_video01_simple_embankment.xlsx](files/xslope_video01_simple_embankment.xlsx)

<!-- test: file=files/xslope_video01_simple_embankment.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.084, tolerance=0.005 -->

**Related tutorial** — [LEM-1. Simple Embankment](../tutorials/lem01_simple_embankment.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
[![](../lem/sample_images/simple_foundation_results.png)](https://youtu.be/mjhU-9LGG4k)
<span class="vid-num">2</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 2 · Limit equilibrium</p>

### Building a Model Three Ways {#video-2}

A 20 ft slope of undrained clay, with the soil continuing 10 ft below the toe to
a rigid base, built three ways: in the Excel input template, in Studio's editors,
and by Studio's AI assistant from the problem sketch. All three builds give the
same model and the same factor of safety, 0.964: the slope as drawn is unstable.

<p class="vid-meta"><a class="vid-pill vid-watch" href="https://youtu.be/mjhU-9LGG4k">▶ Watch · 17 min</a></p>

**Completed model** — [xslope_simple_foundation.xlsx](../lem/files/xslope_simple_foundation.xlsx)

**Related pages** — [Building Models Three Ways](../tutorials/building_models.md) ·
[LEM sample problem 1](../lem/samples.md#1-simple-slope-with-foundation)

</div>
</div>

<div class="vid-band" markdown>

## Limit equilibrium

<span class="vid-band-sub">Videos 3–8</span>

</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
[![](../tutorials/images/lem03_solution.png)](https://youtu.be/mNE2UAJeb4Y)
<span class="vid-num">3</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 3 · Limit equilibrium</p>

### Adding Layers {#video-3}

An embankment on a foundation layer over rigid rock, with a profile line at each
material boundary and a starting circle at the base of each layer. The search
result shows why the critical circle settles on the contact between layers.

<p class="vid-meta"><a class="vid-pill vid-watch" href="https://youtu.be/mNE2UAJeb4Y">▶ Watch · 8 min</a></p>

**Completed model** — [xslope_simple_mult_layers.xlsx](../lem/files/xslope_simple_mult_layers.xlsx)

**Related tutorial** — [LEM-3. A Layered Slope](../tutorials/lem03_layered_slope.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/lem04_solution_wet.png)
<span class="vid-num">4</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 4 · Limit equilibrium</p>

### Water in the Slope {#video-4}

A piezometric line added to a layered slope, with saturated unit weights below
it. The same critical surface is solved with the pore pressures on and off to
show what the water costs.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Starting model** — [xslope_video04_water_start.xlsx](files/xslope_video04_water_start.xlsx)

<!-- test: file=files/xslope_video04_water_start.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.667, tolerance=0.005 -->

**Related tutorial** — [LEM-4. Water in the Slope](../tutorials/lem04_water_in_the_slope.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../lem/sample_images/submerged_slope_results.png)
<span class="vid-num">5</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 5 · Limit equilibrium</p>

### A Submerged Slope {#video-5}

A clay slope under 10 ft of standing water. The water is entered as a
piezometric line above the ground, and Studio converts it into the load the water
applies to the slope surface. Because the clay is analyzed with its undrained
strength, the water adds load without adding pore pressure in the soil.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related page** — [LEM sample problem 2](../lem/samples.md#2-submerged-slope)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/lem02_solution_load.png)
<span class="vid-num">6</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 6 · Limit equilibrium</p>

### Loads and Earthquake {#video-6}

A surcharge on the Video 1 embankment, entered as a distributed load and as a line
load, followed by a seismic coefficient.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [LEM-2. Loads on the Crest](../tutorials/lem02_loads_on_the_crest.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/lem05_solution.png)
<span class="vid-num">7</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 7 · Limit equilibrium</p>

### A Non-Circular Surface Through a Weak Layer {#video-7}

A sand slope over a thin seam of soft clay, with the failure surface entered as
a table of vertices and then searched. The Movement column controls how each
vertex is allowed to move during the search.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [LEM-5. A Weak Layer, Non-Circular](../tutorials/lem05_weak_layer_noncircular.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/lem06_solution_composite.png)
<span class="vid-num">8</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 8 · Limit equilibrium</p>

### Polygon Geometry {#video-8}

An embankment on bedrock that dips across the section, entered as one closed
polygon per material zone instead of profile lines. The dipping base limits the
circles a search can try, and the composite-surface option lets a circle run
along the bedrock instead of being refused.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [LEM-6. Polygon Geometry](../tutorials/lem06_polygon_geometry.md)

</div>
</div>

<div class="vid-band" markdown>

## Seepage

<span class="vid-band-sub">Videos 9–11</span>

</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/seep01_solution_filled.png)
<span class="vid-num">9</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 9 · Seepage</p>

### Confined Seepage Under a Sheetpile {#video-9}

A sheetpile in a silt foundation behind an impervious clay blanket, built from a
blank model: the geometry, the three boundary condition types and the no-flow
default, the mesh, the solution, the discharge, and the flow net.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [SEEP-1. Confined Seepage Under a Sheetpile](../tutorials/seep01_sheetpile.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/seep02_solution.png)
<span class="vid-num">10</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 10 · Seepage</p>

### Unconfined Seepage Through a Zoned Dam {#video-10}

An earth dam with a shell and a clay core, built from a blank model, with the
reservoir head on the upstream face and a seepage face downstream. The solution
locates the phreatic surface.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [SEEP-2. Unconfined Seepage Through a Zoned Dam](../tutorials/seep02_johnson_dam.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/combo01_lem_solution.png)
<span class="vid-num">11</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 11 · Seepage + limit equilibrium</p>

### Seepage Pore Pressures in a Stability Run {#video-11}

The Video 10 dam searched with Spencer's method on the pore pressures from its
seepage solution, by setting `u` to `seep` in the materials table.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [COMBO-1. Seepage → LEM → FEM](../tutorials/combo01_seepage_stability.md)

</div>
</div>

<div class="vid-band" markdown>

## Finite element

<span class="vid-band-sub">Videos 12–13</span>

</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/fem01_shear_strain.png)
<span class="vid-num">12</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 12 · Finite element</p>

### Strength Reduction on an Embankment {#video-12}

An embankment solved by the shear strength reduction method: the elastic
properties the run needs, a quadratic mesh, the run, and the deformed mesh and
shear strain at failure, compared with Spencer's method on the same slope.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [FEM-1. Strength Reduction Basics](../tutorials/fem01_strength_reduction.md)

</div>
</div>

<div class="vid-row" markdown>
<div class="vid-thumb" markdown>
![](../tutorials/images/combo01_fem_shear.png)
<span class="vid-num">13</span>
</div>
<div class="vid-body" markdown>
<p class="vid-kicker">Video 13 · Seepage + LEM + FEM</p>

### Seepage, LEM, and FEM on One Mesh {#video-13}

The Video 10 dam run through all three engines on one quadratic mesh: the seepage
solution, a Spencer search on its pore pressures, and a strength reduction run.

<p class="vid-meta"><span class="vid-pill vid-soon">Coming soon</span></p>

**Related tutorial** — [COMBO-1. Seepage → LEM → FEM](../tutorials/combo01_seepage_stability.md)

</div>
</div>

<!--
Entry format. Once a video is uploaded, replace the "Coming soon" pill with
<a class="vid-pill vid-watch" href="https://youtu.be/VIDEO_ID">▶ Watch · 12 min</a>,
wrap the thumbnail image in the same link, and add the file lines:

**Starting model** — [name_start.xlsx](files/name_start.xlsx)  (only when the video does not start from a blank model)
**Completed model** — [name.xlsx](files/name.xlsx)
-->
