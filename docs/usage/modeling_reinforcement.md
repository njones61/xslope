---
title: "Modeling reinforcement — XSLOPE"
description: "What to enter on the reinforce sheet for a geosynthetic layer, a soil nail, a tieback or an end-anchored bar, how to connect reinforcement to a wall or facing, and when to make a line a joint, with typical values from the FHWA manuals."
---

# Modeling Reinforcement

A support that works in tension enters a model as one row on the `reinforce` sheet of the
[input template](input_template.md#worksheet-reinforce): a straight line, and the columns that give its strength,
its anchorage at each end and its stiffness. What goes in those columns depends on the support (a geosynthetic
layer, a soil nail, a tieback or an end-anchored bar) and on whether the model is analyzed by limit equilibrium
(LEM) or by finite elements (FEM). The mechanics behind each column are set out on
[Soil Reinforcement in LEM](../lem/reinforcement.md) and
[Soil Reinforcement in Finite Element Analysis](../fem/reinforcement.md). A micropile, pile or pier resists by
shear and bending and goes on the `piles` sheet instead ([LEM](../lem/piles.md), [FEM](../fem/piles.md)).

Each row runs from end 1 at (`x1`, `y1`) to end 2 at (`x2`, `y2`); `Lp1` and `Tend1` describe the anchorage at
end 1, and `Lp2` and `Tend2` the anchorage at end 2. The program accepts the ends in either order; in the recipes
below, end 1 is the end at the face or wall where there is one. The LEM applies the line's force at the point
where a trial slip surface crosses it. The FEM builds the line into the mesh as a row of bar elements that share
the soil's nodes along the whole line and start at zero force in the finished geometry, with no construction
sequence ([Initial state and EA selection](../fem/reinforcement.md#initial-state-and-ea-selection)); a bar carries
force along its own length, so `Dir` and `Appl` have no effect there.

`Appl` decides what the LEM's factor of safety F divides. With Appl Active the LEM does not divide the line's force
by F, so every capacity entered is an allowable value, the nominal capacity divided by its own factor of safety:
`Tmax`, the pullout entry (`Lp1` and `Lp2` from an allowable load transfer, or `Adhesion` and `Delta`) and `Tend1`
and `Tend2`. With Appl Passive every capacity is nominal, and the LEM divides all of them by F together with the
soil's strength. A design that applies different factors to different parts, such as GEC 7's 1.8 on a nail bar and
2.0 on pullout (Table 5.1, p. 108), enters allowable values with Appl Active. The FEM's strength reduction divides
the soil's strength only: a bonded line's capacities act as entered whatever `Appl` says, and a jointed line's
interface strength is reduced with the soil unless `Jred` is `No`.

XSLOPE does not enforce a unit system: every value is entered in one consistent set, such as kN and m or lb and
ft. A continuous sheet is entered per unit width of slope, with `Spacing` blank. A discrete member (a nail, a
tieback, a bar) is entered per member, with its out-of-plane spacing in `Spacing`; the program divides `Tmax`,
`Tres`, `Tend1`, `Tend2` and `Area` by `Spacing` when the file is read, and the forces it reports are per unit
width ([Per-unit-width convention and spacing](../lem/reinforcement.md#per-unit-width-convention-and-spacing)).
`E` is not divided. `Adhesion` and `Delta` are entered as they are, and the pullout resistance they produce is
divided by `Spacing` with the forces. Where a manual cited here gives a value in one unit system only, the value
in parentheses is converted from it.

The sheet colors each header by the analysis that reads it. The columns, in sheet order, with F for force and L
for length:

| Column | Name | Units | Read by | Meaning |
|---|---|---|---|---|
| B | `Label` | — | both | name used in messages, plots and reports |
| C, D | `x1`, `y1` | L | both | end 1 |
| E, F | `x2`, `y2` | L | both | end 2 |
| G | `Type` | — | LEM | support preset that fills `Dir` and `Appl` ([Support Type Presets](../lem/reinforcement.md#support-type-presets)) |
| H | `Dir` | — | LEM | direction of the force at a crossing: Tangent to the slip surface, or Axial along the line ([Force Direction](../lem/reinforcement.md#force-direction-dir)) |
| I | `Appl` | — | LEM | Active, an allowable force not divided by the factor of safety; Passive, an ultimate force that is ([Force Application](../lem/reinforcement.md#force-application-appl)) |
| J | `Tmax` | F per member, or F/L | both | tensile capacity |
| K, L | `Lp1`, `Lp2` | L | both | length over which friction develops `Tmax` from end 1 and from end 2; 0 is a fully anchored end |
| M, N | `Adhesion`, `Delta` | F/L², degrees | both | interface adhesion and friction angle; filled together, they replace `Lp1` and `Lp2` ([Pullout from the effective overburden](../lem/reinforcement.md#pullout-from-the-effective-overburden)) |
| O, P | `Tend1`, `Tend2` | F per member, or F/L | both | capacity of a plate, connection or anchorage at end 1 and at end 2 |
| Q | `Spacing` | L | both | out-of-plane spacing of discrete members; blank for a sheet |
| R | `Tres` | F per member, or F/L | FEM | tension the bar keeps after it ruptures: blank for no rupture (the bar holds its capacity), 0 for a brittle break ([Force Behavior and Failure Modes](../fem/reinforcement.md#force-behavior-and-failure-modes)) |
| S | `E` | F/L² | FEM | elastic modulus of the member |
| T | `Area` | L² per member, or L²/L | FEM | cross-sectional area; `E` × `Area` is the axial stiffness ([Axial Stiffness (EA)](../fem/reinforcement.md#axial-stiffness-ea)) |
| U | `Joint` | `Yes` or blank | FEM | `Yes` makes the line a slip surface ([Two Ways to Represent a Sheet](../fem/reinforcement.md#two-ways-to-represent-a-sheet)) |
| V, W | `kn`, `ks` | F/L³ | FEM | normal and shear stiffness of a jointed line's interfaces |
| X | `Jred` | `Yes`, `No` or blank | FEM | whether strength reduction weakens a jointed line's interfaces |

## Geosynthetic Layer

A geotextile or geogrid is a sheet laid flat on a lift of fill as a slope or wall is built. It carries tension
where a slip surface cuts it, and develops that tension by friction on both faces of the length buried beyond the
surface. Its line runs along the layer from the face of the slope, or from the back of a facing, to the buried
end, so the line is as long as the layer. The figure shows one layer in a reinforced slope.

<!-- figure: reinf-geosynthetic -->

The capacity starts at `Tend1` at the face (zero where the layer ends free) and at zero at the buried end, grows
with the friction on both faces of the layer, and is capped at `Tmax`.

| Column | Entry | Typical values |
|---|---|---|
| `Type` | `Geosynthetic`: Dir Tangent, Appl Active; set `Appl` to `Passive` to enter the nominal values instead | — |
| `Tmax` | Appl Active: the long-term allowable strength divided by the target factor of safety, `Tal` ÷ FS<sub>R</sub>. Appl Passive: `Tal`. | `Tal` = `Tult` ÷ RF, with RF the product of the creep, durability and installation-damage factors. In FHWA Example E1, `Tult` = 3,000, 6,000 and 9,000 lb/ft (43.8, 87.6 and 131.3 kN/m) gives `Tal` = 1,085, 2,169 and 3,525 lb/ft (15.8, 31.7 and 51.4 kN/m), with RF<sub>CR</sub> = 1.85, RF<sub>D</sub> = 1.15 and RF<sub>ID</sub> = 1.3, 1.3 and 1.2 (NHI-10-025 Table E1-7.3, p. E1-15). RF = 7 for preliminary design of routine structures in granular fill (p. 9-5). |
| `Adhesion`, `Delta` | `Adhesion` = 0. Appl Active: `Delta` = arctan(F\*α ÷ FS<sub>PO</sub>), with FS<sub>PO</sub> the factor of safety against pullout (NHI-10-025 Eq. 9-9, p. 9-13). Appl Passive: `Delta` = arctan(F\*α). F\* is the pullout resistance factor and α the scale-effect correction. | F\* = 0.67 tan φ with α = 0.6, the most conservative defaults (NHI-10-025 p. E8-6); α = 0.6 to 0.8 for extensible reinforcement without pullout tests (p. B-2); F\* = 0.45 and α = 0.8 for the geogrids of Example E1 (p. E1-16); FS<sub>PO</sub> = 1.5 in granular soil and 2 in cohesive soil, and minimum embedment beyond the critical surface 3 ft (1 m) (p. 9-5) |
| `Lp1`, `Lp2` | in place of `Adhesion` and `Delta`: the length over which friction develops `Tmax` at each end | — |
| `Tend1` | Appl Active: the long-term connection strength `Talc` divided by the factor of safety the design applies to the connection. Appl Passive: `Talc`. 0 where the layer ends free at the face. | see [A geosynthetic and facing blocks, panels or a wrapped face](#a-geosynthetic-and-facing-blocks-panels-or-a-wrapped-face) |
| `Tend2` | 0 | — |
| `Spacing` | blank | — |
| `E`, `Area` (FEM) | `E` × `Area` = the sheet's tensile stiffness per unit width | see [Initial state and EA selection](../fem/reinforcement.md#initial-state-and-ea-selection) |
| `Tres` (FEM) | blank | — |
| `Joint` (FEM) | blank; `Yes` where the soil can slide along the sheet ([A sheet the soil slides along](#a-sheet-the-soil-slides-along)) | — |

The two `Tmax` entries are FHWA's two conventions (NHI-10-025, pp. 8-7 and 8-8): a program that adds the
reinforcement force to the resisting moment takes `Tal` as it stands, and one that subtracts it from the driving
moment takes `Tal` divided by the target factor of safety.

The `Geosynthetic` preset's tangent direction is the one FHWA gives for continuous sheets (p. 8-6). In the FEM the
fill is not placed in lifts, so the tension that construction puts into a layer is absent, and the soil cannot
slide along a layer unless its `Joint` is `Yes`.

[Tutorial LEM-8](../tutorials/lem08_reinforced_slope.md) builds a slope reinforced with six geogrid layers for
limit equilibrium, and [Tutorial FEM-2](../tutorials/fem02_reinforcement.md) runs it by finite elements with
`E` = 800,000 psf and `Area` = 0.1 ft² (`EA` = 80,000 lb/ft), then with `Adhesion` = 0 and `Delta` = 22°.
[FHWA Example E1](../verification/published.md#fhwa-e1) enters a geogrid wall's pullout through `Adhesion` and
`Delta`.

## Soil Nail

A soil nail is a steel bar grouted into a hole drilled into a cut as the cut is excavated from the top down, with a
plate at the head bearing on a shotcrete facing. It carries tension where a slip surface crosses it, developed by
bond between the grout and the soil along its length and, at the head, by the facing. Its line runs from the nail
head on the face, end 1, along the nail's inclination to its tip, end 2; nails are installed 10 to 20 degrees
below horizontal, most commonly 15 (GEC 7, p. 150). The figure shows one nail in a cut.

<!-- figure: reinf-nail -->

At end 1 the capacity starts at the plate's `Tend1` and at end 2 at zero, and bond adds to each over `Lp1` and
`Lp2` until the bar's `Tmax` caps it.

| Column | Entry | Typical values |
|---|---|---|
| `Type` | `Nail`: Dir Axial, Appl Passive, with the nominal values below; set `Appl` to `Active` to enter GEC 7's allowable values instead | — |
| `Tmax` | nominal tensile resistance of the bar per nail, `At` × `fy` (GEC 7 Eq. 6.5, p. 163); with Appl Active, `At` × `fy` ÷ FS<sub>T</sub> | solid threaded bars #6 to #14: area 0.44 to 2.25 in² (284 to 1,452 mm²); yield load 26 to 135 kip (116 to 601 kN) in Grade 60 and 33 to 168 kip (147 to 747 kN) in Grade 75 (GEC 7 Tables A.1a and A.1b, p. 286); FS<sub>T</sub> = 1.8 for Grades 60 and 75 (Table 5.1, p. 108) |
| `Lp1`, `Lp2` | `Tmax` ÷ `rPO` at both ends, with `rPO` = π × `qu` × `DDH` the nominal pullout resistance per unit length, from the bond strength `qu` and the drill-hole diameter `DDH` (GEC 7 Eq. 6.1, p. 161); with Appl Active, `Tmax` ÷ (`rPO` ÷ FS<sub>PO</sub>) | `rPO` = 2 to 20 kip/ft (30 to 290 kN/m) for small-diameter gravity-grouted holes, by soil type and density (GEC 7 Table 4.6, p. 86; GEC 4 Table 6, p. 71); `qu` = 3 to 70 psi (21 to 483 kPa) in soil, by soil type and drilling method (GEC 7 Tables 4.4a and 4.4b, pp. 84–85); FS<sub>PO</sub> = 2.0 (Table 5.1, p. 108) |
| `Adhesion`, `Delta` | blank: GEC 7 states bond as a constant `qu` along the nail, so the development length applies; the `Adhesion` and `Delta` law is written for a sheet with soil on both faces | — |
| `Tend1` | nail-head capacity: the smallest of the facing's resistances in flexure, in punching shear and in headed-stud tension (GEC 7 pp. 104–105); with Appl Active, each resistance divided by its own factor before the smallest is taken | nominal resistances: flexure 12 to 143 kip (53 to 636 kN) for the initial facing; punching shear 32 to 111 kip (142 to 494 kN) initial and 21 to 55 kip (93 to 245 kN) final; headed studs 28 to 146 kip (125 to 649 kN) (GEC 7 Tables 6.6 to 6.8, pp. 171–178); flexure and stud values for Grade 60 steel, × 1.24 for Grade 75; factors of safety 1.5 for flexure and punching shear, 2.0 for A307 and 1.7 for A325 headed studs (Table 5.1, p. 108) |
| `Tend2` | 0 | — |
| `Spacing` | horizontal nail spacing | 4 to 6 ft (1.22 to 1.83 m), routinely 5 ft (1.52 m) (GEC 7 p. 148) |
| `E` (FEM) | modulus of the steel bar | 29,000 ksi (GEC 7 p. 250): 4.176 × 10⁹ psf, or about 2.0 × 10⁸ kPa (200 GPa) |
| `Area` (FEM) | bar area per nail | as for `Tmax` |
| `Tres`, `Joint` (FEM) | blank | — |

The `Nail` preset enters nominal values. GEC 7's allowable stress design, in which programs such as SNAILZ use an
allowable bond stress (p. 117), is entered with `Appl` set to Active and the allowable values in the table; F then
applies to the soil alone, for which GEC 7 sets a minimum of 1.5 (Table 5.1, p. 108). The FEM has no facing member
(see [A soil nail and a shotcrete facing](#a-soil-nail-and-a-shotcrete-facing)).

[VP47](../verification/rocscience.md#vp47) enters this envelope: `Tmax` = 118 kN per nail, `Tend1` = 86 kN,
`Lp1` = `Lp2` = 118 ÷ 15 = 7.87 m from a bond of 15 kN per meter, and `Spacing` = 1.5 m. The 4.9 m nails are
shorter than `Lp1` and `Lp2`, so pullout governs along their whole length.
[VP48](../verification/rocscience.md#vp48) enters each nail as a constant 15 kN tension, with `Lp1` = `Lp2` = 0
and `Spacing` = 1.15 m, the simplification its source analysis makes, so its envelope is flat.

## Tieback (Grouted Ground Anchor)

A tieback, or grouted ground anchor, is a prestressed bar or strand tendon in a grouted hole. Its anchorage (the
anchor head and bearing plate) bears on the wall. Over the unbonded length the tendon is sleeved so that it does
not bond to the grout, and carries the force between the anchorage and the bond length, which transfers it to the
ground behind the critical slip surface (GEC 4, pp. 4–5). Its line runs from the anchor head at the wall, end 1,
through the unbonded length to the far end of the bond length, end 2. The figure shows one tieback through a
wall.

<!-- figure: reinf-tieback -->

With `Lp1` = 0 the full `Tmax` is available from the anchor head to within `Lp2` of end 2, and over that last
`Lp2` the capacity falls linearly to zero, as in GEC 4's limit equilibrium treatment of an anchor (p. 100). Where
the bond governs, as drawn, `Lp2` is the bond length.

| Column | Entry | Typical values |
|---|---|---|
| `Type` | `Tieback`: Dir Axial, Appl Active | — |
| `Tmax` | allowable anchor load per anchor: the smallest of the tendon's design load, the allowable capacity of the head's connection to the wall, and the bond length times the allowable load transfer per unit length | design load at most 0.6 × the tendon's specified minimum tensile strength (GEC 4 p. 77); design loads of 260 to 1,160 kN (58.5 to 260.8 kip) are typical (p. 70) |
| `Lp1` | 0: the head holds the full `Tmax`, and the sleeved unbonded length adds no friction | — |
| `Lp2` | `Tmax` ÷ the allowable load transfer per unit length, which is the ultimate load transfer ÷ 2.0 in soil or ÷ 3.0 in rock (GEC 4 pp. 71, 74): the bond length where the bond governs, shorter where the tendon or the head governs. A 580 kN anchor in medium dense sand (145 kN/m ultimate, 72.5 kN/m allowable) has `Lp2` = 580 ÷ 72.5 = 8.0 m. | ultimate load transfer of small-diameter gravity-grouted anchors: 30 to 290 kN/m (2 to 20 kip/ft) in soil (GEC 4 Table 6, p. 71; GEC 7 Table 4.6, p. 86) and 150 to 730 kN/m (10.3 to 50.0 kip/ft) in rock (GEC 4 Table 8, p. 74); bond lengths 4.5 to 12 m (14.8 to 39.4 ft) in soil and 3 to 10 m (9.8 to 32.8 ft) in rock (pp. 71, 74) |
| `Tend1`, `Tend2`, `Adhesion`, `Delta` | blank: with `Lp1` = 0 the program does not read `Tend1`; for `Adhesion` and `Delta`, see grouted tiebacks under [Pullout from the effective overburden](../lem/reinforcement.md#pullout-from-the-effective-overburden) | — |
| `Spacing` | horizontal anchor spacing | on a soldier-beam wall the anchors connect to the soldier beams, directly or through wales (GEC 4 pp. 13–14); soldier beams are typically 1.5 to 3 m (4.9 to 9.8 ft) apart when driven and up to 3 m apart when drilled in (p. 76) |
| `E` (FEM) | modulus of the tendon steel | 29,000 ksi for a bar tendon (GEC 7 p. 250); for strand, the manufacturer's value, which PTI allows to be reduced 3 to 5 percent for a long multistrand tendon when checking apparent free length (GEC 4 p. 151) |
| `Area` (FEM) | tendon area per anchor | Grade 150 bars 26 to 64 mm (1 to 2½ in.): 548 to 3,348 mm² (0.85 to 5.19 in²), ultimate strength 568 to 3,461 kN (127.5 to 778.0 kip) (GEC 4 Table 9, p. 77); 15-mm strand: 140 mm² (0.217 in²) and 260.7 kN (58.6 kip) per strand (GEC 4 Table 10, p. 78) |
| `Tres`, `Joint` (FEM) | blank | — |

The unbonded length is at least 3 m (9.8 ft) for a bar tendon and 4.5 m (14.8 ft) for strand (GEC 4 p. 70), and
the bond length starts at least one fifth of the wall height or 1.5 m (4.9 ft) behind the critical slip surface
(p. 65). The critical surface the search reports should therefore cross each tieback on its unbonded length, at
least that distance in front of the bond length.

The LEM has no lock-off load: it applies the envelope value where a trial surface crosses the tendon, and nothing
to a surface beyond end 2. The FEM has none either, where designers lock off at 75 to 100 percent of the design
load (GEC 4 p. 154). In the FEM the unbonded length grips the soil as the bond length does, the load transfer
along the unbonded length that an apparent free length below GEC 4's minimum may indicate in a load test
(p. 151). How the tendon reaches the wall is under
[A tieback and a soldier-pile or sheet-pile wall](#a-tieback-and-a-soldier-pile-or-sheet-pile-wall).

[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) builds a soldier-pile tieback wall, the model of
[VP49](../verification/rocscience.md#vp49), with this envelope. It enters capacities per foot of wall (the
per-anchor value ÷ the 8 ft spacing, with `Spacing` = 1), which gives the same force as the per-anchor `Tmax`
with `Spacing` = 8. It enters no separate unbonded length: the bar governs, and `Lp2` is the length of grout that
develops the bar's capacity, 8.87 ft and 12.1 ft. [VP58](../verification/rocscience.md#vp58) has tiebacks 88 ft
long at 20° with a 40 ft bond length that governs at 40,000 lb per foot of wall, so `Lp2` = 40 ft and the 48 ft
in front of the bond length carry the full `Tmax`. Both files enter `Type` = `Anchor`, whose presets are those of
`Tieback`, and neither carries `E` or `Area`, so both run in the LEM only.

## End-Anchored Bar

An end-anchored bar is a steel bar or tie rod held by a plate or a deadman at each end. It carries tension along
its axis when the ground at its two ends moves apart, and is held by its two anchorages, with little or no bond
along the shaft. Its line runs from the anchorage at one end, end 1, to the anchorage at the other, end 2. The
figure shows one bar held by a plate at each end.

<!-- figure: reinf-end-anchored-bar -->

With `Lp1` = `Lp2` = 0 the bar delivers `Tmax` wherever a slip surface crosses it between the two plates.

| Column | Entry | Typical values |
|---|---|---|
| `Type` | `Anchor`: Dir Axial, Appl Active | — |
| `Tmax` | the smallest of the bar's allowable tension and the allowable capacities of its two anchorages | bar sizes and strengths as for a nail or a tieback bar (GEC 7 Tables A.1a and A.1b, p. 286; GEC 4 Table 9, p. 77). The tabulated strengths are nominal; divide by the factor of safety the design uses for the bar (GEC 7 uses 1.8 for Grade 60 and 75 nail bars, Table 5.1, p. 108). |
| `Lp1`, `Lp2` | 0 | — |
| `Tend1`, `Tend2`, `Adhesion`, `Delta` | blank | — |
| `Spacing` | out-of-plane (horizontal) spacing of the bars | — |
| `E`, `Area` (FEM) | steel modulus; bar area per bar | `E` = 29,000 ksi (GEC 7 p. 250) |
| `Tres`, `Joint` (FEM) | blank | — |

With both development lengths 0 the envelope does not read `Tend1` or `Tend2`: the shaft develops no friction
between its anchorages, and their capacities enter through `Tmax`. A bar whose shaft also grips the soil is
entered with the allowable anchorage capacities in `Tend1` and `Tend2` and the development lengths of the
allowable friction in `Lp1` and `Lp2` ([Capacity Envelope](../lem/reinforcement.md#capacity-envelope)). In the FEM the shaft grips the soil
along its whole length, and no plate or deadman is built into the mesh.

No tutorial or verification model has an end-anchored bar.

## Connecting Reinforcement to a Wall or Facing

A support that bears on a wall or a facing is connected at end 1. The two analyses treat the connection
differently: the LEM applies each member's force on its own, and the FEM connects a bar to another member only at
a node the two share or, for a jointed sheet, through a tie at its end.

### A tieback and a soldier-pile or sheet-pile wall

The wall is entered on the `piles` sheet ([Piles and Concrete Piers in LEM](../lem/piles.md)) and each tieback on
the `reinforce` sheet, with end 1 on the wall face. The wall's resistance is its shear force `H`, which GEC 4
takes as the smaller of the wall's shear capacity and the passive force the soil develops below the surface,
divided by the soldier beam spacing (p. 101). The figure shows one tieback through a soldier-pile wall, with a
trial surface that crosses both.

<!-- figure: connect-tieback-wall-lem -->

The tieback's force acts along the tendon at its crossing and the wall's `H` at the pile line.

A trial surface that passes below the toe of the wall receives no force from the wall. In
[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) each tieback starts on the wall face at x = 0, 0.5 ft in
front of the soldier pile line; the offset has no effect in the LEM.

In the FEM the wall is a row of beam elements
([Piles and Concrete Piers in Finite Element Analysis](../fem/piles.md)) and a tieback is a row of bar elements.
A bar shares a node with the wall only where its end 1 is
placed exactly at an end of the pile line, the pile's head or its tip. A bar that starts partway down the pile,
stops short of it, or starts in front of it and passes through it is not attached to the wall, and its force
reaches the wall only through the soil around both. The figure shows one tieback whose end 1 lies on the pile
line below the pile's head.

<!-- figure: connect-tieback-wall-fem -->

The bar's first node lies on the pile line but belongs to the soil alone, so the bar and the wall share no
node.

### A soil nail and a shotcrete facing

A nail's head plate bears on the shotcrete facing at end 1. The LEM takes the facing into account in two ways:
the head's capacity, in `Tend1`, and the facing's weight, entered as a vertical line load at the top of the face
([Worksheet: lloads](input_template.md#worksheet-lloads)). On a vertical face the facing's weight acts along the
face, so one vertical force at the top of the face gives the same force and moment as the weight spread down it,
provided the slip surface exits at the toe; on a battered face it is approximate. The figure shows one nail head
on a shotcrete facing.

<!-- figure: connect-nail-facing -->

The head capacity acts at end 1 and the facing's weight at the crest.

[VP47](../verification/rocscience.md#vp47) and [VP48](../verification/rocscience.md#vp48) enter their facings
this way, with line loads of 14.6 kN/m and 13.2 kN/m. The FEM has no facing member: a beam cannot be laid along
the face, because the mesher rejects a line that runs along the ground surface
([Reinforcement and pile lines](../fem/mesh.md#reinforcement-and-pile-lines)), so each nail head is a soil node on
the face and `Tend1` only raises the capacity of the bar at that end.

### A geosynthetic and facing blocks, panels or a wrapped face

A layer connected to facing blocks or panels starts at the back of the facing, end 1. FHWA takes the long-term
connection strength `Talc` from connection tests on the facing unit and the geosynthetic (NHI-10-025 Eq. 4-41,
p. B-13), and it rises with the normal pressure on the connection: in Example E1 it runs from 533 lb/ft
(7.8 kN/m) at the top layer to 2,550 lb/ft (37.2 kN/m) at the bottom, against `Tal` = 1,085 and 2,169 lb/ft
(15.8 and 31.7 kN/m) for the two geogrid grades (Table E1-7.3, p. E1-15; connection strengths Table E1-7.6,
p. E1-18), so on eight of the wall's eleven layers the connection limits the force at the face. The figure shows
one layer connected to a block facing.

<!-- figure: connect-geosynthetic-facing -->

The connection acts at end 1; end 2 is buried and free.

In the LEM, and in the FEM for a bonded layer, `Tend1` only raises the capacity at end 1
([Capacity Envelope](../lem/reinforcement.md#capacity-envelope)); a jointed layer's `Tend1` ties its end to the
facing it runs into, up to that capacity ([Ends, ties and the bar](../fem/reinforcement.md#ends-ties-and-the-bar)).
In the block wall of [Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md#part-2-the-same-wall-with-geogrid)
the back face of the facing is a joint line, and a bonded line may not end on a joint line, so all three layers
are jointed, each tied to the blocks at `Tend1` = 40 kN/m. A wrapped face has no facing unit: the sheet is folded
back over the face and buried under the next lift
([When a Model Needs a Joint](../fem/joints.md#when-a-model-needs-a-joint)). Its line starts at the face, and
`Tend1` is the pullout resistance of the return length folded back and buried under the next lift.

## Jointed Reinforcement and Joint Lines

`Joint` = `Yes` makes a reinforcement line a slip surface, and `kn`, `ks` and `Jred` set its interfaces; the
`joints` sheet holds slip surfaces with no reinforcement in them. The LEM reads neither: it treats a jointed line
as an ordinary reinforcement line, from `Tmax`, `Tend1`, `Tend2` and either `Lp1` and `Lp2` or `Adhesion` and
`Delta`, and it ignores the `joints` sheet.

### A sheet the soil slides along

Where a slip surface can run along a sheet (a base geotextile under an embankment on soft clay, a smooth
geomembrane liner, the sheets of a wrapped or block-faced wall), `Joint` = `Yes` splits the mesh along the line
and the soil on each side slides on the sheet at the interface strength
([Choosing a Bonded Bar or a Joint](../fem/reinforcement.md#bonded-bar-or-joint)). The figure shows one base sheet
under an embankment.

<!-- figure: joint-base-sheet -->

The sheet lies between two interfaces, one against the embankment above it and one against the foundation below;
both ends are free, so the sheet can pull out at either end.

| Column | Entry |
|---|---|
| `Joint` | `Yes` |
| `Adhesion`, `Delta` | both required: the cohesion and friction angle of the two interfaces, whose tension cutoff is zero |
| `kn`, `ks` | blank: derived from the softer adjacent material over a notional thickness of one tenth of the element length ([Stiffness](../fem/joints.md#stiffness)) |
| `Jred` | blank, so strength reduction weakens the interface along with the soil; `No` holds it at full strength |
| `Tend1`, `Tend2` | a value above zero ties that end at that capacity; blank or 0 leaves the end free |
| `Lp1`, `Lp2` | not read |

The model checks warn when a bonded sheet looks like a slip surface
([Signs that a joint is needed](../fem/reinforcement.md#what-says-a-joint-was-needed)).
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md#part-3-when-a-sheet-is-a-slip-surface-and-when-it-is-bonded)
runs a base geotextile and a liner both ways, with `kn`, `ks` and `Jred` blank on every jointed line.

### A joint with no reinforcement

A slip surface with no member in it (a rock joint, a contact between facing blocks, the back of a wall against its
fill) goes on the `joints` sheet, with its own `c`, `phi` and `t_cut`
([Joints and Interface Elements](../fem/joints.md)). The figure shows one contact between two facing blocks.

<!-- figure: joint-no-reinforcement -->

The mesh splits along the line into two faces with one interface between them, and no member.

A bonded reinforcement line or a pile may not end on or cross a jointed line
([A reinforcement line as a joint](../fem/joints.md#a-reinforcement-line-as-a-joint)), so a sheet that ends on such
a contact is jointed itself. The FEM-3 block wall enters its base, its back face and the five course joints between
its six blocks on the `joints` sheet.

## References

Berg, R.R., Christopher, B.R., & Samtani, N.C. (2009). *Design of Mechanically Stabilized Earth Walls and
Reinforced Soil Slopes – Volume II*. FHWA-NHI-10-025 (FHWA GEC 011, Vol. II). Federal Highway Administration,
Washington, D.C.

Lazarte, C.A., Robinson, H., Gómez, J.E., Baxter, A., Cadden, A., & Berg, R. (2015). *Geotechnical Engineering
Circular No. 7: Soil Nail Walls Reference Manual*. FHWA-NHI-14-007. Federal Highway Administration, Washington,
D.C.

Sabatini, P.J., Pass, D.G., & Bachus, R.C. (1999). *Geotechnical Engineering Circular No. 4: Ground Anchors and
Anchored Systems*. FHWA-IF-99-015. Federal Highway Administration, Washington, D.C.
