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
it differently ([LEM vs FEM](#lem-vs-fem)). What to enter for each kind of support, with typical
values from the design manuals, is on [Reinforcement Types](types.md). A micropile, pile or pier resists by shear
and bending rather than tension and is entered on the `piles` sheet instead ([LEM](../piles/lem.md),
[FEM](../piles/fem.md)).

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
([Per-unit-width convention and spacing](lem.md#per-unit-width-convention-and-spacing)). `E`,
`Adhesion` and `Delta` are entered as they are, and the pullout resistance computed from `Adhesion` and `Delta` is
divided by `Spacing` in the same way.

## The Columns

Each name carries the sheet's header color, which shows the analysis that reads the column. In the Units column F
is a force and L a length; elsewhere on this page F is the factor of safety.

<p class="rc-legend">Used by: <span class="rc rc-geom">geometry (both)</span><span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Name | Units | Meaning |
|---|---|---|---|
| B | <code class="rc rc-geom">Label</code> | — | name used in messages, plots and reports; optional |
| C, D | <code class="rc rc-geom">x1</code>, <code class="rc rc-geom">y1</code> | L | end 1 |
| E, F | <code class="rc rc-geom">x2</code>, <code class="rc rc-geom">y2</code> | L | end 2 |
| G | <code class="rc rc-lem">Type</code> | — | support preset that fills `Dir` and `Appl` ([Support Type Presets](lem.md#support-type-presets)) |
| H | <code class="rc rc-lem">Dir</code> | — | direction of the force at a crossing: Tangent to the slip surface, or Axial along the line ([Force Direction](lem.md#force-direction-dir)) |
| I | <code class="rc rc-lem">Appl</code> | — | Active: allowable capacities, not divided by the factor of safety.<br>Passive: nominal capacities, divided by it.<br>([Force Application](lem.md#force-application-appl)) |
| J | <code class="rc rc-both">Tmax</code> | F per member, or F/L | tensile capacity |
| K, L | <code class="rc rc-both">Lp1</code>, <code class="rc rc-both">Lp2</code> | L | length over which friction develops `Tmax` from end 1 and from end 2; 0 makes the full `Tmax` available at that end |
| M, N | <code class="rc rc-both">Adhesion</code>, <code class="rc rc-both">Delta</code> | F/L², degrees | interface adhesion and friction angle; filled together, they replace `Lp1` and `Lp2`, and filling only one is an input error ([Pullout from the effective overburden](lem.md#pullout-from-the-effective-overburden)) |
| O, P | <code class="rc rc-both">Tend1</code>, <code class="rc rc-both">Tend2</code> | F per member, or F/L | capacity of a plate, connection or anchorage at end 1 and at end 2; 0 for none |
| Q | <code class="rc rc-both">Spacing</code> | L | out-of-plane spacing of discrete members; blank for a sheet |
| R | <code class="rc rc-fem">Tres</code> | F per member, or F/L | tension the bar keeps after it ruptures: blank for no rupture (the bar holds its capacity), 0 for a brittle break, a value between for a bar that keeps that much ([Force Behavior and Failure Modes](fem.md#force-behavior-and-failure-modes)) |
| S | <code class="rc rc-fem">E</code> | F/L² | elastic modulus of the member |
| T | <code class="rc rc-fem">Area</code> | L² per member, or L²/L | cross-sectional area; `E` × `Area` is the axial stiffness ([Axial Stiffness (EA)](fem.md#axial-stiffness-ea)) |
| U | <code class="rc rc-fem">Joint</code> | `Yes` or blank | `Yes` makes the line a slip surface, with `Adhesion` and `Delta` as its interface strength ([Two Ways to Represent a Sheet](fem.md#two-ways-to-represent-a-sheet)) |
| V, W | <code class="rc rc-fem">kn</code>, <code class="rc rc-fem">ks</code> | F/L³ | normal and shear stiffness of a jointed line's interfaces; blank derives them from the adjacent soil ([Stiffness](../fem/joints.md#stiffness)) |
| X | <code class="rc rc-fem">Jred</code> | `Yes`, `No` or blank | blank or `Yes` reduces a jointed line's interface strength with the soil's in a strength reduction; `No` holds it at full strength |

## Capacity Along a Line

The tension a line can carry varies along it. In the middle it is the member's own `Tmax`. Toward each end it is
limited by what that end can develop: the plate or connection capacity `Tend1` or `Tend2`, plus the friction or bond
along the line from that end, which on its own rises linearly to `Tmax` over `Lp1` or `Lp2`. The smallest of the three
at each point is the line's envelope, and both analyses use the same one ([Capacity
Envelope](lem.md#capacity-envelope)).

`Adhesion` and `Delta` are the alternative to `Lp1` and `Lp2`. Instead of a fixed development length they state the
interface strength, and the pullout resistance then follows the effective overburden along the line,
$2(a + \sigma'_v\tan\delta)$ per unit length from the sheet's two faces, with $a$ the adhesion, $\delta$ the
friction angle and $\sigma'_v$ the vertical effective stress at each point.

## LEM vs FEM

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

The formulations are on [Soil Reinforcement in LEM](lem.md) and
[Soil Reinforcement in Finite Element Analysis](fem.md).

## Worked Examples

[Tutorial LEM-8](../tutorials/lem08_reinforced_slope.md) reinforces a slope with six geogrid layers, and
[Tutorial FEM-2](../tutorials/fem02_reinforcement.md) runs the same slope by finite elements and compares the two.
[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) builds a soldier-pile wall held by tiebacks, and
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) ties geogrid layers into a block wall built on slip
joints and shows when a sheet is better modeled as a jointed line than as a bonded one.
