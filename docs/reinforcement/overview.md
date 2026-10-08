---
title: "Reinforcement — XSLOPE"
description: "How reinforcement enters an XSLOPE model: the reinforce sheet and Studio's reinforcement editor, what each column means, the capacity along a line, and how the limit equilibrium and finite element analyses use it."
---

# Reinforcement

Soil carries compression and shear but little or no tension, and reinforcement adds it. A member laid or drilled
across the zone where a slip surface would form, long enough to grip the stable ground behind that surface, pulls
back on the sliding mass. The common kinds are geosynthetic layers (geotextiles and geogrids) built into a fill as
it is placed, soil nails grouted into a cut as it is excavated, and tiebacks holding a wall with a grouted length
deep behind the slip surface.

![Three kinds of reinforcement: geosynthetic layers in a fill, soil nails in a cut, tiebacks behind a wall](../fem/images/reinf_types.png){width=1000}

From left to right: geosynthetic layers in a reinforced fill, soil nails in a cut, and tiebacks behind a wall,
each with a slip surface they cross.

XSLOPE models each of these, and an end-anchored bar held by a plate or deadman at each end, as a straight line
with a tensile capacity. The limit equilibrium (LEM) and finite element (FEM) analyses read the same line but use
it differently ([The LEM and the FEM](#the-lem-and-the-fem)). What to enter for each kind of support, with typical
values from the design manuals, is on [Reinforcement Types](types.md). A micropile, pile or pier resists by shear
and bending rather than tension and is entered on the `piles` sheet instead ([LEM](../lem/piles.md),
[FEM](../fem/piles.md)).

## Entering a Line

Each line is one row of the `reinforce` sheet of the
[input template](../usage/input_template.md#worksheet-reinforce):

![The reinforce sheet of the input template](../usage/images/sheet_reinforce.png)

In XSLOPE Studio the same rows are edited in the reinforcement editor, whose list view shows one line at a time
beside a preview of the section ([Editing Inputs](../studio/editing.md#the-inputs-tree-and-the-editors)):

![The reinforcement editor in list view](../studio/images/editing_reinforcement_editor.png)

A line runs from end 1 at (`x1`, `y1`) to end 2 at (`x2`, `y2`). `Lp1` and `Tend1` describe the anchorage at end 1,
and `Lp2` and `Tend2` the anchorage at end 2. The program accepts the ends in either order. `Type` is a preset:
picking `Geosynthetic`, `Nail`, `Tieback` or `Anchor` fills `Dir` and `Appl` with that support's usual pair, and
either can then be changed. A blank `Type` is a generic tensile line, with `Dir` set to Tangent and `Appl` to
Active ([The Columns](#the-columns)).

XSLOPE does not enforce a unit system: every value is entered in one consistent set, such as kN and m or lb and
ft. A continuous sheet is entered per unit width of slope, with `Spacing` blank. A discrete member (a nail, a
tieback, a bar) is entered per member, with its out-of-plane spacing in `Spacing`; the program divides its
capacities and `Area` by that spacing and reports forces per unit width
([Per-unit-width convention and spacing](../lem/reinforcement.md#per-unit-width-convention-and-spacing)). `E`,
`Adhesion` and `Delta` are entered as they are, and the pullout resistance computed from `Adhesion` and `Delta` is
divided by `Spacing` in the same way.

## The Columns

The sheet colors each header by the analysis that reads it: green for the LEM only, red for both, blue for the FEM
only, and black for the geometry. In the Units column F is a force and L a length; elsewhere on this page F is
the factor of safety.

| Column | Name | Units | Read by | Meaning |
|---|---|---|---|---|
| B | `Label` | — | both | name used in messages, plots and reports; optional |
| C, D | `x1`, `y1` | L | both | end 1 |
| E, F | `x2`, `y2` | L | both | end 2 |
| G | `Type` | — | LEM | support preset that fills `Dir` and `Appl` ([Support Type Presets](../lem/reinforcement.md#support-type-presets)) |
| H | `Dir` | — | LEM | direction of the force at a crossing: Tangent to the slip surface, or Axial along the line ([Force Direction](../lem/reinforcement.md#force-direction-dir)) |
| I | `Appl` | — | LEM | Active: allowable capacities, not divided by the factor of safety; Passive: nominal capacities, divided by it ([Force Application](../lem/reinforcement.md#force-application-appl)) |
| J | `Tmax` | F per member, or F/L | both | tensile capacity |
| K, L | `Lp1`, `Lp2` | L | both | length over which friction develops `Tmax` from end 1 and from end 2; 0 makes the full `Tmax` available at that end |
| M, N | `Adhesion`, `Delta` | F/L², degrees | both | interface adhesion and friction angle; filled together, they replace `Lp1` and `Lp2`, and filling only one is an input error ([Pullout from the effective overburden](../lem/reinforcement.md#pullout-from-the-effective-overburden)) |
| O, P | `Tend1`, `Tend2` | F per member, or F/L | both | capacity of a plate, connection or anchorage at end 1 and at end 2; 0 for none |
| Q | `Spacing` | L | both | out-of-plane spacing of discrete members; blank for a sheet |
| R | `Tres` | F per member, or F/L | FEM | tension the bar keeps after it ruptures: blank for no rupture (the bar holds its capacity), 0 for a brittle break, a value between for a bar that keeps that much ([Force Behavior and Failure Modes](../fem/reinforcement.md#force-behavior-and-failure-modes)) |
| S | `E` | F/L² | FEM | elastic modulus of the member |
| T | `Area` | L² per member, or L²/L | FEM | cross-sectional area; `E` × `Area` is the axial stiffness ([Axial Stiffness (EA)](../fem/reinforcement.md#axial-stiffness-ea)) |
| U | `Joint` | `Yes` or blank | FEM | `Yes` makes the line a slip surface, with `Adhesion` and `Delta` as its interface strength ([Two Ways to Represent a Sheet](../fem/reinforcement.md#two-ways-to-represent-a-sheet)) |
| V, W | `kn`, `ks` | F/L³ | FEM | normal and shear stiffness of a jointed line's interfaces; blank derives them from the adjacent soil ([Stiffness](../fem/joints.md#stiffness)) |
| X | `Jred` | `Yes`, `No` or blank | FEM | blank or `Yes` reduces a jointed line's interface strength with the soil's in a strength reduction; `No` holds it at full strength |

## Capacity Along a Line

The tension a line can carry varies along it. In the middle it is the member's own `Tmax`. Toward each end it is
limited by what that end can develop: the plate or connection capacity `Tend1` or `Tend2`, plus the friction or bond
along the line from that end, which on its own rises linearly to `Tmax` over `Lp1` or `Lp2`. The smallest of the three
at each point is the line's envelope, and both analyses use the same one ([Capacity
Envelope](../lem/reinforcement.md#capacity-envelope)).

`Adhesion` and `Delta` are the alternative to `Lp1` and `Lp2`. Instead of a fixed development length they state the
interface strength, and the pullout resistance then follows the effective overburden along the line,
$2(a + \sigma'_v\tan\delta)$ per unit length from the sheet's two faces, with $a$ the adhesion, $\delta$ the
friction angle and $\sigma'_v$ the vertical effective stress at each point.

## The LEM and the FEM

The LEM applies the envelope's value at the point where a trial slip surface crosses the line, as a force on the
sliding mass, in the direction `Dir` sets: tangent to the slip surface or along the line. The FEM builds the line into
the mesh as a row of bar elements. Its force is a result rather than an input: each bar carries the tension its
stretch produces, up to the envelope, and acts along its own length, so `Dir` and `Appl` have no effect there. The
figure compares the two on the [FEM-2](../tutorials/fem02_reinforcement.md) slope and its critical circle:

![The direction of the reinforcement force in limit equilibrium and in the finite element analysis](../fem/images/reinf_force_direction.png){width=1000}

On the left, the limit equilibrium force at each crossing acts tangent to the circle, and the axial alternative is
drawn at the middle layer. On the right, each bar carries tension along its own length, from how much it stretches.

`Appl` decides what the LEM's factor of safety F divides. With Appl Active the LEM does not divide the line's force
by F, so every capacity entered is an allowable value, the nominal capacity divided by its own factor of safety:
`Tmax`, the pullout entry and `Tend1` and `Tend2`. With Appl Passive every capacity is nominal, and the LEM divides
all of them by F together with the soil's strength. A design that applies a different factor to each part, such
as a [soil nail](types.md#soil-nail) design with one factor on the bar and another on pullout, enters allowable
values with Appl Active.

The FEM's strength reduction divides the soil's strength only. Its bars start at zero force in the finished
geometry, with no construction sequence. A bonded line shares the soil's nodes along its whole length, and its
capacities act as entered, as they do in the LEM with Appl Active; with Appl Passive the LEM divides them by F and
the FEM does not. A line with `Joint` = `Yes`, a jointed line, is split from the soil instead: its interface
strength is reduced with the soil's unless `Jred` is `No`, and its `Tmax`, `Tend1` and `Tend2` act as entered.

The formulations are on [Soil Reinforcement in LEM](../lem/reinforcement.md) and
[Soil Reinforcement in Finite Element Analysis](../fem/reinforcement.md).

## Worked Examples

[Tutorial LEM-8](../tutorials/lem08_reinforced_slope.md) reinforces a slope with six geogrid layers, and
[Tutorial FEM-2](../tutorials/fem02_reinforcement.md) runs the same slope by finite elements and compares the two.
[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) builds a soldier-pile wall held by tiebacks, and
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) ties geogrid layers into a block wall built on slip
joints and shows when a sheet is better modeled as a jointed line than as a bonded one.
