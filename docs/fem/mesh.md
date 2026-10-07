---
title: "Mesh generation — XSLOPE"
description: "How XSLOPE meshes a section for seepage and finite element analysis: element types, quadrilateral styles, element size and local refinement, and embedded reinforcement and pile lines."
---

# Mesh Generation

Seepage and finite-element analyses run on a mesh of triangles or quadrilaterals
covering the section. Limit-equilibrium analysis does not need one. You build the mesh
explicitly — in Studio with **Build Mesh**, or in a script with
`build_mesh_from_polygons()` — and Studio saves it beside the model as `{stem}_mesh.json`,
so one mesh serves every run until the geometry changes.

Meshing happens in two stages. First the model's geometry becomes a set of closed
material-zone polygons whose shared boundaries match exactly. Then
[gmsh](https://gmsh.info) fills each polygon with elements and welds the zones together
along those shared boundaries, so the mesh conforms across every material interface and
no element straddles one. gmsh is an optional dependency:
`pip install "xslope[fem]"`.

## Building a mesh

A mesh is built in Studio or from a script; both run the same mesher.

### In Studio

In **Seepage** or **FEM** mode, **Build Mesh** opens a dialog that sets the element type, the
target element size (entered, or auto-sized from the section width), the 1D element size, the
quadrilateral style, and refinement near features and in thin zones:

![Build Mesh dialog](../studio/images/analysis_build_mesh_dialog.png)

Auto-sizing divides the width of the section by a number of divisions, 100 by default. The
mesh is built on a background thread and appears in a **Mesh** tab; **Run** stays disabled
until a mesh exists, and a geometry edit that invalidates the mesh clears it. The dialog is
described control by control under
[Studio → Building a mesh](../studio/analysis.md#building-a-mesh).

### From a script

The same build in Python:

```python
from xslope.fileio import load_slope_data
from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,
                         extract_constraint_line_geometry, extract_size_regions,
                         get_material_polygons)

slope_data = load_slope_data("xslope_earth_dam1.xlsx")

lines, n_reinf, n_pile = extract_constraint_line_geometry(slope_data)
polygons = get_material_polygons(slope_data, reinf_lines=lines or None)

mesh = build_mesh_from_polygons(
    polygons,
    target_size=1.1,
    element_type="tri6",
    lines=lines or None,
    size_regions=extract_size_regions(slope_data),
)
export_mesh_to_json(mesh, "xslope_earth_dam1_mesh.json")
```

The returned mesh is a dictionary of NumPy arrays:

| Key | Contents |
| --- | --- |
| `nodes` | `(n_nodes, 2)` node coordinates |
| `elements` | `(n_elements, 9)` node indices per element; unused slots are 0 |
| `element_types` | nodes per element — 3, 4, 6, 8 or 9 |
| `element_materials` | material ID per element, 1-based in `mat` sheet order |
| `elements_1d`, `element_types_1d`, `element_materials_1d` | the same three arrays for embedded 1D elements, present only when `lines` is given |

## From input file to material polygons

A model defines its geometry one of two ways, and the two are mutually exclusive — a file
that populates both sheets is rejected:

- **Polygon sheet.** Each material zone is drawn directly as a closed polygon. The
  polygons *are* the geometry, and they go to the mesher as they were entered.
- **Profile sheet.** Each material boundary is drawn as a profile line running left to
  right, entered from the top of the section downwards. `build_polygons()` closes them
  into zones: it projects each line's endpoints onto the boundary below to create the
  connecting vertices, then traces each zone between one line and the next, with the
  bottom zone closed off at `max_depth`.

Either way, `get_material_polygons()` is the single entry point that returns
mesh-ready zones, and it does the same preparation for both paths:

- Reinforcement and pile line vertices are inserted into the boundary edges they touch,
  so the mesh can recover the line as a chain of element edges.
- Distributed-load endpoints and seepage boundary-condition vertices are inserted the
  same way. Without them an element edge straddling a load end is dropped by the edge-load
  integrator and less load is applied than the model asked for.
- A zone pinched out to zero thickness — a dam core whose top line runs along the
  foundation beyond the core, for example — produces a self-touching ring that gmsh
  refuses. Those rings are repaired, and a zone split in two by a trench or key is kept as
  two regions carrying the same material.

Inside the mesher, one more pass makes adjacent zones conforming: where a zone's vertex
lands in the interior of a neighbor's edge, that edge is split there, so the interface
meshes without a slit. Every polygon becomes one gmsh surface tagged with its material,
and points shared between zones are created once and reused, so the resulting mesh is
continuous across material boundaries.

## Element types

Five element types are available. The linear types, `tri3` and `quad4`, carry their
solution variable linearly across the element; the quadratic types, `tri6`, `quad8` and
`quad9`, add midside nodes and carry it quadratically, which resolves a curving field
with far fewer elements.

For the soil types, node ordering is corner nodes first, counterclockwise, then the
midside nodes in edge order (0–1, 1–2, and so on), then — for `quad9` alone — the center
node. Counterclockwise corners give a positive Jacobian everywhere, which the element
integration relies on. Each row of `elements` holds the nodes in that order, and the
matching entry of `element_types` says how many of the nine slots are used.

Use the local indices below to read the node order in the mesh arrays.

![Soil and line elements with local node indices](images/all_element_nodes.png){width=800}

The indices are zero-based, and a line element numbers its two ends before its midpoint.

| Type | Nodes | Variation within the element | Built from |
| --- | --- | --- | --- |
| `tri3` | 3 | linear | — |
| `tri6` | 6 | quadratic | `tri3` |
| `quad4` | 4 | bilinear | — |
| `quad8` | 8 | quadratic (serendipity) | `quad4` |
| `quad9` | 9 | quadratic (Lagrange) | `quad4` |

The last column gives the linear type each quadratic type is built from
([Quadratic elements](#quadratic-elements)).

### Triangles or quadrilaterals

Both fill any section. The difference is what they cost and how they behave:

![The same section meshed with triangles and with quadrilaterals](images/mesh_tri_vs_quad.png)

At one requested element size the two meshes carry a similar number of nodes — which is
what sets the size of the solve — but the quadrilateral mesh reaches it with roughly half
as many elements. Triangles fit an irregular boundary with less distortion and never fail
to fill a corner; quadrilaterals give a more regular element layout wherever the shape
allows.

A quadrilateral mesh is quad-*dominant* rather than pure: a small fraction of elements,
typically one or two percent, stays triangular where no pairing exists. XSLOPE carries
mixed meshes end to end, so those triangles solve like any other element.

### Element choice for FEM analyses {#element-choice-for-fem-analyses}

Use a quadratic type for every stress analysis: linear elements lock volumetrically and
overstate the factor of safety
([volumetric locking](overview.md#element-type-selection-and-volumetric-locking)). The
choice among the three:

- `tri6` conforms to complex geometry where quads would distort and is preferred for
  submerged problems, where `quad8`'s reduced integration allows an hourglass mode —
  a zero-energy deformation pattern that 2×2 integration cannot resist.
- `quad8` with 2×2 reduced integration is the element and integration Griffiths & Lane
  (1999) used for their SSRM benchmarks; it avoids
  locking and gives accurate stress fields. It gives a regular layout in block-like sections.
- `quad9` with full 3×3 integration is also suitable for a block-like section, at the cost
  of extra Gauss points and a center node.

`build_mesh_from_polygons()` defaults to `tri6`, and Studio's *Build mesh* dialog opens
on it. A blank **Element type** cell (`main!D18`) gives `tri6`.
An explicit type on the call or main sheet overrides the default.

For seepage, `tri3` is the lighter choice, with a smaller system that solves faster.

The model checks issue a warning, not an error, before a FEM or SSRM solve starts on a linear mesh.

## Quadratic elements

Quadratic meshes are built in two steps: gmsh generates the linear mesh (`tri3` for
`tri6`, `quad4` for `quad8` and `quad9`), then `convert_linear_to_quadratic_mesh()` adds
the extra nodes.

Each unique edge in the linear mesh receives one midside node at its midpoint, shared by
both elements that use the edge, so the mesh stays conforming. A `quad9` also gets a
center node at the average of its four corners. A leftover triangle in a `quad8` or
`quad9` mesh becomes a `tri6` — the mesh stays mixed rather than forcing a bad
quadrilateral. Embedded 1D elements pick up the midside node of the 2D edge they lie on
and become 3-node: a bar or beam couples to the soil only where they share a node, so an
element on the corners alone would leave the soil free to bow away from it between them.
An element on an edge that has no midside node — a linear mesh — stays 2-node.

## Quadrilateral meshing styles

A quadrilateral element type can be meshed in either of two styles, chosen per run:
`build_mesh_from_polygons(..., quad_style='free' | 'structured')` in a script, or the
**Quadrilateral style** radio group in Studio's *Build mesh* dialog. Triangular element
types ignore the setting. Both styles are driven by the same requested element size, so
switching between them changes how the elements are arranged, not how big they are.

**Free** is the default and the right choice for most sections. gmsh's
Frontal-Delaunay-for-quads algorithm lays down a triangulation whose points are placed so
that pairs of triangles form near-square quadrilaterals, and the Blossom algorithm of
[Remacle et al. (2012)](https://doi.org/10.1002/nme.3279) decides which triangles to pair
by solving a global minimum-cost perfect matching over the whole zone. It works on any
shape and needs nothing of the geometry. On a model with reinforcement, pile or joint
lines the pairing is gmsh's simple recombination instead of Blossom.

**Structured where possible** adds a sweep in front of that. Each material zone is tested
against a conservative mappability check — four logical sides, opposite sides of
compatible length, corners near square, and a predicted element aspect ratio no worse
than 2:1 — and every zone that passes is filled with a regular grid of rows and columns.
The row and column counts come from the requested element size, and a boundary shared by
two zones gets exactly one count, so a swept zone and its free-meshed neighbor meet at the
same node spacing with no hanging nodes. Choose it when the section is built of block-like
zones — a layered foundation, a cutoff or grout curtain, a rectangular core — where rows
of aligned elements are useful.

A zone the check declines is meshed exactly as in the default style; the other zones are
unaffected. The three sections below show the effect. Each is meshed with `quad4`
elements at the same requested size in both panels, the section width divided by 100.
AR95 is the 95th-percentile element aspect ratio — longest corner edge over shortest — so
lower is better, and the last figure in each caption is the median delivered element size
as a multiple of the size requested.

**A levee on a blocked foundation.** The three foundation blocks sweep as exact grids. The
embankment is a trapezoid whose base and crest differ by about 5:1, and a grid swept
through that shape would have to stretch its elements by the same ratio, so it is declined
and stays free-meshed in both panels.

![Free and structured quadrilateral meshes of the levee section](images/quad_styles_levee.png)

**An earth dam with a central clay core.** Only the core is mappable; the shell that wraps
around it is not. The swept core is visibly regular, and the aspect ratio of the mesh as a
whole improves, but the free-meshed shell now has to close against a fixed grid, which
leaves a few more triangles behind than the free mesh did.

![Free and structured quadrilateral meshes of the earth dam](images/quad_styles_earth_dam.png)

**A dredged trench in silt.** No zone passes the check — the trench gives each of them
more corners than the four a swept grid needs — so every zone falls back to the free
mesher and the two panels are the same mesh, element for element.

![Free and structured quadrilateral meshes of the sea trench section](images/quad_styles_sea_trench.png)

The style is not stored in the input file. A model re-opened elsewhere meshes in the
default style unless a mesh built at the other setting travels with it as the companion
`{stem}_mesh.json`.

## Element size

The element size comes from a global target, a Size on a zone, automatic refinement near
features and in thin zones, and the 1D element size along reinforcement and pile lines
([below](#reinforcement-and-pile-lines)).

### The global target size

`target_size` is the requested element size, in model length units, and both triangular
and quadrilateral meshing work at that size directly. The **Mesh target size** cell on the
`main` sheet sets the value the *Build mesh* dialog opens with.

One background size field decides the element size everywhere outside a structured
quadrilateral sweep, for both element families and whether or not anything is being
refined: the requested size in the far field, and a graded band around each refined
feature. Delivered node spacing along the geometry runs 0.75 to 1.00 times the requested
size across the sample and verification sections; the shortfall is gmsh rounding a curve
up to a whole number of divisions, which can only make an edge finer.

### A Size on one zone

A material zone can be given a finer size of its own through the **Size** cell on its
polygon (or on its profile line). Everything outside keeps the global target size, and the
size grades back to it across the zone boundary:

![A local Size on the dam core](images/mesh_zone_size.png)

A Size only ever refines — a value at or above the global target has no effect, and
XSLOPE warns rather than silently ignoring it. A polygon whose **Type** is `refine`
carries nothing *but* a Size: it is not a material, never becomes a mesh region, and only
makes the mesh finer where it is drawn. Refine polygons reach the mesher through
`extract_size_regions()`, passed as `size_regions`. Both are described on the
[Input Template](../usage/input_template.md#refine-regions) page.

A Size resolves a zone it is entered on, however thin the zone is — the size field
reaches the zone's boundary as well as its interior. A Size does not locate thin zones; to
find them, and for a section whose thin bands have not been checked by hand, use the
automatic thin-zone refinement below.

### Refining near features

Some models concentrate their numerically demanding behavior in a few small features —
load transfer along a reinforcement or pile line, the singularity at a crack or notch tip,
the kinematics inside a thin soft band. Refining the whole domain to resolve them is
wasteful. `refine_factor` shrinks elements only where a feature is present and lets them
grow smoothly back to the target size elsewhere:

```python
mesh = build_mesh_from_polygons(
    polygons,
    target_size=2.9,
    element_type="tri3",
    lines=lines,
    refine_factor=3.0,             # local size = target_size / 3; None (default) = off
    refine_features=None,          # None = reinforcement, piles, cracks, thin zones
)
```

Both panels are the same reinforced slope at the same target size.

![A mesh refined near reinforcement lines](images/mesh_refine_features.png)

In the lower one every one of the six reinforcement lines carries a band of finer
elements, and the thin shell along the face — detected as a thin zone, not as a line — is
refined along its whole length.

Features are detected from the model geometry; there are no coordinates to enter by hand.
The classes are:

- **Reinforcement and pile lines** — the embedded 1D lines get a distance-based band,
  finest along the polyline and coarsening away from it.
- **Crack and notch tips** — the deepest vertex of a slit, where the material wraps
  almost all the way around the point, such as the notch a sheet-pile wall is modeled
  with. Tips are refined twice as strongly, to resolve the singularity that governs
  convergence there. A sharp *convex* corner is not a tip and is not refined: a layer
  tapering to a point across the section, or an embankment toe, has the same edge
  directions as a slit and none of the physics.
- **Thin material zones** — a zone whose local width cannot fit three elements at the
  target size, refined as described under [Thin material zones](#thin-material-zones).
- **High-contrast material interfaces** — a boundary between two zones whose major
  hydraulic conductivities differ by 100× or more, where a seepage solve must resolve a
  steep head gradient. This class is **opt-in**: it needs conductivities that only a
  seepage model carries, so select it explicitly and pass them,
  `refine_features=['interfaces'], material_k={0: 1.67e-5, 1: 1.67e-7}`, keyed by 0-based
  material index (the first `mat` row is 0).

`refine_factor=None`, the default, adds no band: the background size field is the
requested size everywhere and the mesh is what it would be without the option. Refinement
composes with the field rather than editing the geometry, and detection is pure geometry,
so a refined mesh is reproducible from run to run. Each band holds its local size for
about two element widths and then grows back to the target at 1.2 per element, so a
refined region joins the far field gradually instead of stepping.

Because the field alone sets the size, a feature is meshed at the size that was
requested: at a sheet-pile tip and along a reinforcement line the delivered size is
within about 15 % of the request on either family.

### Thin material zones {#thin-material-zones}

A thin zone is the one refinement case that is not an efficiency question. A soft seam
one element thick cannot carry a shear band, so the model finds no mechanism through it
and the analysis returns a factor of safety that is too high, with nothing in the factor of
safety itself to show that the mesh was the reason. Because of this, Studio's
**Build mesh** dialog carries a **Refine thin zones** checkbox that is **on by default**;
it sizes every thin zone for about four element rows across its local width. A section
with no thin zone is meshed exactly as it would be with the box clear, and the Log names
each zone that was refined with the local size it received.

What counts as the zone's width is the **material's** thickness, not a polygon's. A
single layer is routinely stored as several polygons — a benched face or a step in the
ground surface cuts it into pieces — and each piece is thinner than the layer. Polygons
sharing a material are measured together, so a layer the global element size already
resolves is left alone even where its individual pieces would not be.

The refinement factor plays no part: a thin zone's size is its own thickness over four,
and the same at any factor. Two limits apply to that size. It is never coarser than the
global target, so a zone the global size already resolves is not refined. And it is never
finer than one sixth of the global target.

The model checks report a zone that this one-sixth cap leaves under-resolved. They measure the mesh that
was actually built and name any zone carrying fewer than three element rows, whatever the
reason, so an under-resolved band shows up in the checks panel rather than in a factor of
safety nobody can explain. The two ways to give such a zone the size it needs are a
**Size** on the zone, which is not capped, and a finer global target. A **Size** declared
on the zone composes with the derived size by taking the smaller of the two.

On the Griffiths soft-band section
([xslope_griffiths3_r0p2_thin.xlsx](files/xslope_griffiths3_r0p2_thin.xlsx)) at a 3-unit
target size, the band carries 1.6 element rows on a triangular mesh with the box clear and
1.2 on a quadrilateral one; with the option on both carry 4.0.

From a script, pass `refine_features=['thin_zones']` together with any `refine_factor`
above 1, which switches refinement on. A `size` on the polygon dict, or a `size_regions`
entry, is the manual equivalent for a zone you have measured yourself.

## Reinforcement and pile lines

Reinforcement lines, pile lines and other constraint lines are passed to
`build_mesh_from_polygons()` as `lines`, a list of polylines.
`extract_constraint_line_geometry(slope_data)` returns them from a loaded model in that
form, along with the count of each kind.

The lines are embedded in the 2D mesh rather than meshed separately: gmsh is required to
place element edges along each polyline, and the 1D elements are then extracted from those
edges. For bonded reinforcement and piles, every 1D node is a node of the surrounding 2D mesh, and load transfers
between the reinforcement and the soil through these shared nodes. A line that runs along the domain
boundary — on the base at `max_depth`, or along the ground surface — or that extends
outside the section cannot be embedded; it is rejected before meshing starts, with a
message naming the line, rather than producing a mesh the line is silently missing from.
A line end within one twentieth of the target size of a polygon vertex is moved onto that
vertex, and a line no longer than that is refused.

Node spacing along a line comes from the size field like everything else.
`element_size_1d` (the main sheet's **1D element size**, edited in Studio's *Build mesh*
dialog) sets a finer size along the lines, applied as a graded band as described under
[Structural elements](overview.md#structural-elements). It composes with the field by
taking the minimum, so a value at or above `target_size` has no effect, a zone's finer
local Size is kept, and the 1D size refines only the stretches of a line that are coarser
than it. `element_materials_1d` numbers the lines 1, 2, 3, … in the order they were
passed, so each line's elements can be given their own properties.

On the pile model (`docs/fem/files/xslope_piles_fem.xlsx`, tri6 at a target size of 2,
two pile lines totaling 35 ft) the two lines carry 18 elements at their default spacing;
`element_size_1d=0.5` brings them to 70, and the mesh grows from 3,180 nodes / 1,521
elements to 6,029 / 2,928.

## gmsh options

`mesh_params` is a pass-through for quadrilateral meshes: every entry is a gmsh option
name set verbatim before the mesh is generated, overriding the corresponding default. It
has no effect on a triangular mesh.

```python
mesh = build_mesh_from_polygons(
    polygons, target_size=1.5, element_type="quad4",
    mesh_params={"Mesh.Smoothing": 20},        # extra smoothing passes
)
```

The defaults are those described under
[Quadrilateral meshing styles](#quadrilateral-meshing-styles), so this option is rarely
needed.

## Saving and reusing a mesh

`export_mesh_to_json(mesh, path)` and `import_mesh_from_json(path)` write and read the
mesh as JSON. Studio writes `{stem}_mesh.json` next to the model automatically, and that
file is what a seepage or FEM run loads, so a mesh is built once and reused. The
[sample models](samples.md) ship with their meshes for exactly this reason.

`verify_mesh_connectivity(mesh)` and `print_mesh_connectivity_report(mesh)` inspect a
mesh after the fact: nodes duplicated at one location, nodes no element uses, and
elements that reference the same node twice.
