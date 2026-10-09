---
title: "Piles and walls — XSLOPE"
description: "How piles, piers and walls enter an XSLOPE model: the piles sheet and Studio's pile editor, what each column means, how the limit equilibrium and finite element analyses treat a pile line, and which analysis suits which member."
---

# Piles and Walls

A pile or a wall resists a slide by shear and bending rather than by tension. A row of stabilizing piles or drilled
piers is installed through the sliding mass into stable ground below it; as the mass moves, it pushes against the
piles, and they push back where the slip surface crosses them. A sheet-pile or soldier-pile wall does the same over
its full height and often holds tiebacks. A segmental block wall or a concrete gravity wall is a mass of its own,
which stands by its weight and its grip on the ground beneath it.

![Three kinds of pile and wall: two rows of drilled shafts in a slope, a sheet-pile wall in a slope, and a segmental block wall with joint lines and geogrid](images/piles_walls_types.png){width=1000}

From left to right: two rows of drilled shafts stabilizing a slope ([Tutorial LEM-12](../tutorials/lem12_piles.md)), a
sheet-pile wall on a bench in a slope ([the SIGMA/W wall benchmark](../verification/geostudio.md#sigmaw-wall)), and a
segmental block wall with its joint lines and geogrid ([Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md)).

XSLOPE's limit equilibrium (LEM) and finite element (FEM) analyses take these in two ways. A pile, a micropile, a pier, or a sheet-pile or soldier-pile wall is a straight line
on the `piles` sheet. A block wall or a gravity wall is drawn as polygons of its own material, with joint lines where
it can slide or open in the finite element analysis. What to enter for each type, with typical values from the
design manuals, and how to connect a tieback or a facing to a wall, is on [Pile and Wall Types](types.md).

## Entering a Pile Line

Each pile line is one row of the `piles` sheet of the
[input template](../usage/input_template.md#worksheet-piles):

![The piles sheet of the input template](../usage/images/sheet_piles.png)

In XSLOPE Studio the same rows are edited in the pile editor, whose list view shows one pile at a time beside a
preview of the section ([Editing Inputs](../studio/editing.md#the-inputs-tree-and-the-editors)):

![The pile editor in list view](../studio/images/editing_piles_editor.png)

A pile runs between (`x1`, `y1`) and (`x2`, `y2`), entered in either order: the higher end is the head and the lower
end the tip. A vertical pile has `x1` = `x2`; a battered pile leans. The template is formatted for 20 piles, and
more rows can be added below them.

The user picks a unit system, SI or Imperial, on the `main` sheet
([Worksheet: main](../usage/input_template.md#worksheet-main)), and enters every value in it: kN and m, or lb and ft.
The pile force `H` is entered per unit width of slope, as every force in a two-dimensional analysis is. The section
and capacity columns, `Vcap`, `Mcap`, `I` and `Area`, are entered for a single pile, and the spacing `S` converts
between the two: the force on one pile is `H` × `S`, and the FEM divides a pile's flexural stiffness EI and axial
stiffness EA by `S`. A continuous wall is entered per unit length of wall, with `S` = 1.

## The Columns

Each name carries the sheet's header color, which shows the analysis that reads the column. In the Units column F
is a force and L a length.

<p class="rc-legend">Used by: <span class="rc rc-geom">geometry (both)</span><span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Name | Units | Meaning |
|---|---|---|---|
| B | <code class="rc rc-geom">Label</code> | — | name used in messages, plots and reports; optional |
| C, D | <code class="rc rc-geom">x1</code>, <code class="rc rc-geom">y1</code> | L | one end of the pile |
| E, F | <code class="rc rc-geom">x2</code>, <code class="rc rc-geom">y2</code> | L | the other end |
| G | <code class="rc rc-lem">H</code> | F/L | pile force per unit width of slope, acting perpendicular to the pile where the slip surface crosses it. Blank, with `D` and `S`, computes it for a vertical pile by Ito & Matsui's closed-form solution for the soil squeezing between piles ([Ito & Matsui (1975) Theory](lem.md#ito-matsui-1975-theory)). |
| H | <code class="rc rc-lem">Appl</code> | — | Active: `H` is an allowable force, not divided by the factor of safety.<br>Passive: `H` is a nominal force, divided by it with the soil's strength.<br>Blank reads as Active. |
| I | <code class="rc rc-both">D</code> | L | pile diameter. Ito & Matsui needs it; the FEM computes `I` and `Area` from it when those are blank. |
| J | <code class="rc rc-both">S</code> | L | center-to-center spacing of the piles in the row. Ito & Matsui and the capacity checks need it; the FEM divides EI and EA by it. 1 for a continuous wall. |
| K | <code class="rc rc-both">Vcap</code> | F per pile | shear capacity; blank for no limit ([Structural Capacity Checks](lem.md#structural-capacity-checks)) |
| L | <code class="rc rc-both">Mcap</code> | F·L per pile | moment capacity; blank for no limit. The LEM caps the force on a pile at `Mcap` divided by the moment arm; the FEM forms a plastic hinge where the moment reaches it, a point at which the section yields and rotates freely. |
| M | <code class="rc rc-fem">E</code> | F/L² | elastic modulus of the pile material |
| N | <code class="rc rc-fem">I</code> | L⁴ per pile | moment of inertia of the section; blank, with `D`, is πD⁴/64 |
| O | <code class="rc rc-fem">Area</code> | L² per pile | cross-sectional area; blank, with `D`, is πD²/4 |
| P | <code class="rc rc-fem">Head</code> | `free`, `pinned`, `unrotated` or `fixed` | restraint at the head; blank is free ([Head and Tip Fixity](fem.md#head-and-tip-fixity)) |
| Q | <code class="rc rc-fem">Tip</code> | the same four | restraint at the tip; blank is free |

## LEM vs FEM

The LEM applies one force, `H`, at the point where a trial slip surface crosses the pile, perpendicular to the pile
and against the movement of the sliding mass. `H` is entered, or computed by Ito & Matsui from the soil squeezing
between neighboring piles, and capped by `Vcap` and `Mcap`. The FEM builds the pile into the mesh as a chain of beam
elements that share the soil's nodes. No force is prescribed: the beam carries whatever the deforming soil pushes
onto it, along its whole length, and the analysis returns the moment, shear and deflection down the member. The
figure compares the two on the [LEM-12](../tutorials/lem12_piles.md) slope and its critical circle:

![The pile force in limit equilibrium, one force H at each pile's crossing of the slip surface, and in the finite element analysis, the soil's pressure along the whole pile](images/pile_lem_fem.png){width=1000}

On the left, the LEM applies one force $H$ at each pile's crossing of the circle. On the right, the FEM's beam takes
the soil's pressure along its whole length: the sliding mass pushes it downslope above the slip surface, and the
stable ground holds it below.

A two-dimensional analysis is plane strain, so every member in it is continuous out of plane. A continuous wall is
exactly that, but a row of separate piles is not:

![A continuous wall in section, with S = 1, beside a row of piles of diameter D at spacing S in plan, whose stiffnesses EA and EI per pile become EA/S and EI/S per unit width](images/pile_row_plane_strain.png){width=880}

On the left, a wall in section, continuous out of plane, entered with `S` = 1. On the right, a row of piles in plan,
of diameter `D` at spacing `S`. The FEM divides each pile's EA and EI by `S`, which makes the row a wall of the same
average stiffness, with no gap for the soil to move through between piles. In the ground, the soil arches onto the
piles and, at wide enough spacing, moves between them, and a two-dimensional analysis cannot represent either. That
decides which analysis suits which member:

| Member | Out of plane | Analysis | What it gives |
|---|---|---|---|
| Sheet-pile, diaphragm or secant wall | continuous | FEM, `S` = 1 | factor of safety, plus moment, shear, deflection and soil reaction down the member |
| Contiguous or very closely spaced row | nearly continuous | FEM, which treats the row as a continuous wall | the same, with the gaps unrepresented |
| Row of separate piles | discrete | LEM with Ito & Matsui | factor of safety for the spacing, force per pile, capacity checks |

Take the factor of safety for a discrete pile row from the LEM with Ito & Matsui, and read the FEM's result for the
row as a study of stiffness and member forces. The two analyses are compared on the same slopes, including the one
pile-stabilized slope with a published three-dimensional answer, under
[LEM vs FEM Pile Modeling](lem.md#lem-vs-fem-pile-modeling). The formulations are on [Piles in LEM](lem.md) and
[Piles in FEM](fem.md).

## Worked Examples

[Tutorial LEM-12](../tutorials/lem12_piles.md) stabilizes a slope with two rows of drilled shafts, with `H` from
Ito & Matsui, and [Tutorial FEM-4](../tutorials/fem04_piles.md) runs the same slope by finite elements.
[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) builds a soldier-pile wall held by tiebacks, and
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) a segmental block wall on slip joints, with and without
geogrid.
