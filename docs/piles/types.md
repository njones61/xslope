---
title: "Pile and wall types — XSLOPE"
description: "What to enter for a row of stabilizing piles or drilled shafts, micropiles, a sheet-pile wall, a soldier-pile wall, a segmental block wall and a concrete gravity wall, with typical values from the FHWA manuals; load-bearing piles near a slope; and how to connect tiebacks, nails and geogrid to a wall or facing."
---

# Pile and Wall Types

A row of stabilizing piles or drilled shafts, a row of micropiles, a sheet-pile wall and a soldier-pile wall each
enter the model as one line on the `piles` sheet, with the columns described on the
[Piles and Walls Overview](overview.md#the-columns), but each fills those columns its own way. A segmental block
wall and a concrete gravity wall are drawn as polygons instead, with joint lines where they can slide. A wall that
holds tiebacks, nails or geogrid connects to them as under
[Connecting Reinforcement to a Wall or Facing](#connecting-reinforcement-to-a-wall-or-facing). F is the factor of
safety, and LEM and FEM are the limit equilibrium and finite element analyses.

Typical values come from the Federal Highway Administration (FHWA) manuals on drilled shafts (Geotechnical
Engineering Circular 10, cited as [GEC 10][gec10]), micropiles (the [micropile manual][mp]), ground anchors and
anchored walls ([GEC 4][gec4]), soil nail walls ([GEC 7][gec7]) and mechanically stabilized earth walls
([GEC 11][gec11v1]), and, for steel
sections, one sheet-pile maker's published tables ([Nucor Skyline][nucor]) ([References](#references)). Where a
source gives a value in one unit system only, the value in parentheses is converted from it, and a value given in
inches or millimeters (ksi, in², in⁴, mm²) is also given in model units, built on feet or meters.

## A Row of Stabilizing Piles or Drilled Shafts

A row of drilled shafts, or piers, is installed through the sliding mass and into stable ground below it, at a
spacing that lets the soil arch between them; the shafts act as shear dowels across the slip surface
([GEC 10][gec10] p. 12-59). A discrete row is analyzed with the LEM and the Ito & Matsui force
([LEM vs FEM](overview.md#lem-vs-fem)). The figure shows one shaft of a row in section, and the row in plan.

![In section, a drilled shaft through a slope, from its head at the ground surface, across the slip surface, to its tip in stable ground, with the force H at the crossing; in plan, a row of shafts of diameter D at spacing S](images/pw_shaft.png){width=842}

In section, the pile line runs from the head to the tip, and where a trial slip surface crosses the shaft, the LEM
applies `H`, pointing into the slope, against the movement of the sliding mass. In plan, the shafts of diameter `D`
stand in a row across the slope at spacing `S`, center to center, and the sliding mass moves past them.

<p class="rc-legend">Used by: <span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Entry | Typical values |
|---|---|---|
| <code class="rc rc-lem">H</code> | Blank: Ito & Matsui computes it from `D` and `S` for each trial surface (vertical piles only).<br>Otherwise the force per unit width from a lateral analysis of the pile, the force on one pile ÷ `S`. [GEC 10][gec10] sizes a row for the smaller of the passive force on the shafts and the force the slope needs to reach its target factor of safety (p. 12-58). | None given in the sources. |
| <code class="rc rc-lem">Appl</code> | `Active`: `H` is not divided by F, as in [Tutorial LEM-12](../tutorials/lem12_piles.md). | — |
| <code class="rc rc-both">D</code> | shaft diameter | Drilled-shaft casing and tools come in 6 in (152 mm), or 0.5 ft (0.152 m), steps ([GEC 10][gec10] p. 4-11); GEC 10's wall example uses 4 ft (1.22 m) shafts (pp. 12-14, 12-15). |
| <code class="rc rc-both">S</code> | center-to-center spacing along the row | 3 diameters is routine practice for groups of shafts, 2.5 sometimes advantageous ([GEC 10][gec10] p. 14-2). |
| <code class="rc rc-both">Vcap</code> | shear capacity of one shaft; blank for no limit | [GEC 10][gec10] gives the design method, not typical values. |
| <code class="rc rc-both">Mcap</code> | nominal moment capacity of one shaft | About 27 D³ kip-ft with 1 percent longitudinal steel and 40 D³ with 1.5 percent, D in ft ([GEC 10][gec10] p. 12-14): 1,290 D³ and 1,915 D³ kN·m with D in m. For D = 4 ft (1.22 m), 1.73 × 10⁶ to 2.56 × 10⁶ lb-ft (2,340 to 3,470 kN·m). Longitudinal steel is typically 1 to 2 percent of the gross area (pp. 16-4, 16-5). |
| <code class="rc rc-fem">E</code> | modulus of the concrete | E<sub>c</sub> = 1,820 √f′<sub>c</sub> ksi, with f′<sub>c</sub> the concrete's compressive strength in ksi, generally 3.5 to 5.0 ([GEC 10][gec10] pp. 16-3, 16-4): 3,400 to 4,070 ksi, or 4.9 × 10⁸ to 5.9 × 10⁸ psf (2.35 × 10⁷ to 2.81 × 10⁷ kPa). |
| <code class="rc rc-fem">I</code>, <code class="rc rc-fem">Area</code> | Blank: computed from `D` for the gross section, πD⁴/64 and πD²/4. | For D = 4 ft: I = 12.6 ft⁴ (0.109 m⁴), Area = 12.6 ft² (1.17 m²). |
| <code class="rc rc-fem">Head</code> | Blank (free), unless a cap beam ties the heads (`unrotated`) or anchors hold them (`pinned` or `fixed`). | — |
| <code class="rc rc-fem">Tip</code> | Blank (free): the embedment below the slip surface holds the shaft. `fixed` only for a tip socketed into rock. | — |

## Micropiles

A micropile is a small drilled and grouted pile, typically less than 300 mm (12 in) in diameter, reinforced with a
steel casing or bar ([micropile manual][mp] p. 1-4). To stabilize a slope, micropiles are installed in rows, often
in pairs battered across the slip surface, with their heads tied together by a concrete cap beam at the ground
surface (pp. 6-44, 6-54). Ito & Matsui applies to vertical piles only, so a battered micropile's `H` is entered,
not computed. The figure shows one battered pair under its cap beam in section, and three pairs in plan.

![In section, a pair of micropiles battered in opposite directions from a cap beam on a bench in a slope, each crossing the slip surface with its own force H perpendicular to it; in plan, three pairs along the cap beam at spacing S](images/pw_micropiles.png){width=801}

In section, each leg is a pile line of its own, from its head in the cap beam to its tip in stable ground, with its
own `H` perpendicular to the leg where the slip surface crosses it. In plan, the pairs stand along the cap beam at
spacing `S`, and each pair's legs run upslope and downslope from it below the ground.

<p class="rc-legend">Used by: <span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Entry | Typical values |
|---|---|---|
| <code class="rc rc-lem">H</code> | the resistance per unit width: the force one micropile develops at the slip surface ÷ `S`. Each leg of a pair is a line of its own, with its own `H`. | In the manual's slope example, 365 kN (82 kip) for the upslope leg and 450 kN (101 kip) for the downslope leg, against 650 kN/m (44.5 kip/ft) required (pp. 6-51, 6-52). |
| <code class="rc rc-lem">Appl</code> | `Active` for an allowable `H`.<br>`Passive` for an ultimate `H`, such as the example's, which the LEM then divides by F. | — |
| <code class="rc rc-both">D</code> | Blank: `H`, `I` and `Area` are entered, so the diameter is not read. | Grouted diameter typically less than 0.3 m (1 ft) (p. 1-4); 0.2 m (0.66 ft) in the manual's group example (p. 5-25). |
| <code class="rc rc-both">S</code> | spacing along the row, between micropiles or between pairs | At least 0.76 m (2.5 ft) or 3 diameters, whichever is greater (p. 5-8); 1.25 m (4.1 ft) between pairs in the slope example (p. 6-52). |
| <code class="rc rc-both">Vcap</code> | shear capacity of one micropile | 365 kN (82 kip) at zero axial load for the slope example's 177.8 mm (7 in) casing (p. 6-51). |
| <code class="rc rc-both">Mcap</code> | moment capacity of one micropile | 161 kN·m (119,000 lb-ft) at zero axial load for the same casing ([micropile manual][mp] Table 6-2, p. 6-47). |
| <code class="rc rc-fem">E</code> | modulus of the steel casing | 29,000 ksi: 4.18 × 10⁹ psf, or about 2.0 × 10⁸ kPa ([micropile manual][mp] p. 5-16) |
| <code class="rc rc-fem">Area</code> | area of the steel casing | 3,760 to 8,760 mm² (5.8 to 13.6 in²) for 139.7 to 244.5 mm (5.5 to 9.625 in) casings (Table 4-5, p. 4-33): 0.00376 to 0.00876 m², or 0.0405 to 0.0943 ft². |
| <code class="rc rc-fem">I</code> | moment of inertia of the casing, π(OD⁴ − ID⁴)/64, with OD and ID its outside and inside diameters | For the same casings, 8.05 × 10⁻⁶ to 5.94 × 10⁻⁵ m⁴, or 9.32 × 10⁻⁴ to 6.88 × 10⁻³ ft⁴. |
| <code class="rc rc-fem">Head</code> | `unrotated` where a cap beam ties the heads, as in the slope example (p. 6-54); blank (free) otherwise | — |
| <code class="rc rc-fem">Tip</code> | Blank (free): the bond length below the slip surface holds it. `fixed` for a tip socketed into rock. | The slope example sockets its tips 4.5 m (15 ft) into bedrock (p. 6-57). |

The grout inside the casing is left out of `I` and `Area` here; the moment capacity above includes it.

## Sheet-Pile Wall

A sheet-pile wall is a line of interlocking steel sheets driven to form a continuous wall, cantilevered or held by
tiebacks. It is continuous out of plane, so the FEM represents it directly, and it is entered per unit length of
wall, with `S` = 1. The figure shows a wall held by a tieback in section, and its sheets in plan.

![A sheet-pile wall in section, retaining soil above an excavation, held by a tieback and embedded below the excavation's base, beside the interlocked Z-shaped sheets in plan](images/pw_sheet_pile.png){width=766}

The pile line runs from the top of the wall to its tip, which lies below the excavation's base by the embedment.
The tieback is a line on the `reinforce` sheet, connected as under
[A tieback and a soldier-pile or sheet-pile wall](#a-tieback-and-a-soldier-pile-or-sheet-pile-wall). In plan the
interlocked sheets form one continuous wall, so `I`, `Area` and `Mcap` are entered per unit length of it.

<p class="rc-legend">Used by: <span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Entry | Typical values |
|---|---|---|
| <code class="rc rc-lem">H</code> | the wall's resistance per unit length where a slip surface crosses it: the smaller of its shear capacity and the passive force the soil develops below the excavation ([GEC 4][gec4] p. 101) | — |
| <code class="rc rc-lem">Appl</code> | `Active`: `H` is not divided by F (an allowable resistance). | — |
| <code class="rc rc-both">D</code> | Blank: not read for a wall, whose `I` and `Area` are entered. | — |
| <code class="rc rc-both">S</code> | 1: the values are per unit length of wall. | — |
| <code class="rc rc-both">Vcap</code> | Blank: no limit. | The sources give no typical value. |
| <code class="rc rc-both">Mcap</code> | plastic moment per unit length of wall, F<sub>y</sub> × Z, with F<sub>y</sub> the steel's yield strength and Z the plastic section modulus | F<sub>y</sub> = 50 ksi (345 MPa) is the most common grade ([Nucor Skyline][nucor] p. 6); values for four sections in the table below. |
| <code class="rc rc-fem">E</code> | modulus of the steel | 29,000 ksi: 4.18 × 10⁹ psf, or about 2.0 × 10⁸ kPa ([micropile manual][mp] p. 5-16) |
| <code class="rc rc-fem">I</code>, <code class="rc rc-fem">Area</code> | per unit length of wall, from the section table | the table below |
| <code class="rc rc-fem">Head</code> | Blank (free); tiebacks are entered on the `reinforce` sheet and connected as under [A tieback and a soldier-pile or sheet-pile wall](#a-tieback-and-a-soldier-pile-or-sheet-pile-wall). | — |
| <code class="rc rc-fem">Tip</code> | Blank (free): the embedment holds the toe. | — |

Four common hot-rolled PZ (Z-shaped) sections, per foot and per meter of wall ([Nucor Skyline][nucor] p. 5; `Mcap` computed with
F<sub>y</sub> = 50 ksi):

| Section | `Area` per ft (per m) | `I` per ft (per m) | `Mcap` per ft (per m) |
|---|---|---|---|
| PZ 22 | 0.0449 ft² (0.0137 m²) | 0.00407 ft⁴ (1.15 × 10⁻⁴ m⁴) | 90,800 lb-ft (404 kN·m) |
| PZ 27 | 0.0551 ft² (0.0168 m²) | 0.00888 ft⁴ (2.52 × 10⁻⁴ m⁴) | 152,000 lb-ft (676 kN·m) |
| PZ 35 | 0.0715 ft² (0.0218 m²) | 0.0174 ft⁴ (4.93 × 10⁻⁴ m⁴) | 238,000 lb-ft (1,060 kN·m) |
| PZ 40 | 0.0817 ft² (0.0249 m²) | 0.0237 ft⁴ (6.70 × 10⁻⁴ m⁴) | 300,000 lb-ft (1,330 kN·m) |

The manufacturer states these as 6.47 to 11.77 in² and 84.4 to 491 in⁴ per foot of wall.

## Soldier-Pile Wall

A soldier-pile wall is a row of steel beams, driven H-piles or pairs of channels or wide-flange beams set in
concrete-filled drilled holes, with timber lagging spanning between them to hold the soil ([GEC 4][gec4] p. 13).
The beams are discrete, but the lagging makes the wall continuous, so the wall is entered per beam with `S` the
beam spacing. [Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) builds one held by tiebacks. The figure shows
the wall in section and in plan.

![A soldier-pile wall in section, the lagging down to the excavation's base and the beam embedded below it, and in plan, steel beams at spacing S with timber lagging between them](images/pw_soldier_pile.png){width=658}

In section, the lagging stops at the excavation's base, and the beam continues below it as the embedment. In plan,
the soldier piles, steel beams, stand at spacing `S`, and the lagging spans between them and holds the soil. The
pile line is one beam, and `I`, `Area` and `Mcap` are entered for one beam.

<p class="rc-legend">Used by: <span class="rc rc-lem">LEM only</span><span class="rc rc-both">LEM and FEM</span><span class="rc rc-fem">FEM only</span></p>

| Column | Entry | Typical values |
|---|---|---|
| <code class="rc rc-lem">H</code> | the wall's resistance per unit length: the smaller of the beam's allowable shear capacity and the passive force the soil develops below the excavation, divided by the beam spacing ([GEC 4][gec4] p. 101) | — |
| <code class="rc rc-lem">Appl</code> | `Active`: `H` is not divided by F (an allowable resistance), as in [Tutorial LEM-9](../tutorials/lem09_tieback_wall.md). | — |
| <code class="rc rc-both">D</code> | Blank: `H`, `I` and `Area` are entered. | A drilled-in beam's hole is 0.61 m (2 ft) in GEC 4's examples (pp. A-10, A-29). |
| <code class="rc rc-both">S</code> | center-to-center spacing of the beams | 1.5 to 3 m (4.9 to 9.8 ft) for driven beams, up to 3 m for drilled-in beams ([GEC 4][gec4] p. 76); 2.5 m (8.2 ft) in GEC 4's examples (p. A-7). |
| <code class="rc rc-both">Vcap</code> | Blank: no limit. | The sources give no typical value. |
| <code class="rc rc-both">Mcap</code> | moment capacity of one beam about its strong axis, F<sub>y</sub> × Z. For allowable stress design GEC 4 takes the allowable bending stress F<sub>b</sub> = 0.55 F<sub>y</sub> (p. A-10). | Grade 50 (F<sub>y</sub> = 50 ksi) steel H-piles (HP sections): HP 12×84, 500,000 lb-ft (678 kN·m); HP 14×117, 808,000 lb-ft (1,096 kN·m) ([Nucor Skyline][nucor] pp. 6, 33). |
| <code class="rc rc-fem">E</code> | modulus of the steel | 29,000 ksi: 4.18 × 10⁹ psf, or about 2.0 × 10⁸ kPa ([micropile manual][mp] p. 5-16) |
| <code class="rc rc-fem">I</code>, <code class="rc rc-fem">Area</code> | of one beam, strong axis | HP 12×84: `I` = 650 in⁴, or 0.0314 ft⁴ (2.71 × 10⁻⁴ m⁴); `Area` = 24.6 in², or 0.171 ft² (0.0159 m²). HP 14×117: `I` = 1,220 in⁴, or 0.0588 ft⁴ (5.08 × 10⁻⁴ m⁴); `Area` = 34.4 in², or 0.239 ft² (0.0222 m²) ([Nucor Skyline][nucor] p. 33). |
| <code class="rc rc-fem">Head</code> | Blank (free); tiebacks are connected as under [A tieback and a soldier-pile or sheet-pile wall](#a-tieback-and-a-soldier-pile-or-sheet-pile-wall). | — |
| <code class="rc rc-fem">Tip</code> | Blank (free): the embedment holds the toe. | — |

## Segmental Block Wall

A segmental, or modular block, wall is a column of dry-stacked concrete units, usually with geogrid layers laid
between the courses and running back into a reinforced fill. Its strength lies in its contacts: each course can
slide on the one below it, the column can slide on its foundation, and its back face can part from the fill. The
blocks themselves do not fail. The figure shows the wall of [Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) in section, and in
plan.

![In section, a segmental block wall with a joint line under its base, up its back face and between each pair of courses, and geogrid layers tied into the blocks; in plan, the units laid over the joints of the course below, and a geogrid layer running back into the reinforced fill](images/pw_block_wall.png){width=653}

In section, each course is a polygon of block material, and a joint line runs under the base, up the back face and
between each pair of courses; each geogrid layer ends in a block, at the dot, where it is tied to it. The plan
looks down on the wall with its face at the bottom: the solid lines are the joints between the units of the top
course, the dashed lines those of the course below, offset by half a unit, and the geogrid runs back from the blocks
into the reinforced fill, here along the whole length of the wall.

The wall is drawn, not entered on the `piles` sheet:

- **Blocks.** One polygon per course on the `polygon` sheet, in a block material whose strength option is
  `elastic` ([Worksheet: mat](../usage/input_template.md#worksheet-mat)). An elastic material cannot fail: in the
  FEM the strength reduction has nothing in the blocks to weaken, and in the LEM a trial slip surface may run along
  the wall's boundary but not through it. Give each block polygon a local `Size` of about half the course height,
  so the mesh resolves every course without refining the whole section.
- **Contacts** (FEM). One row on the `joints` sheet for the base, one for the back face, and one for each course
  joint ([The joints worksheet](../fem/joints.md#the-joints-worksheet)). A battered wall, each course set back from
  the one below, has a step in its back face at every course, and each step is a joint line of its own.
- **Geogrid.** Each layer is a line on the `reinforce` sheet, entered as under
  [A geosynthetic and facing blocks, panels or a wrapped face](#a-geosynthetic-and-facing-blocks-panels-or-a-wrapped-face).

The LEM does not read the `joints` sheet, so sliding on a course joint or at the base is a mechanism only the FEM
can find; the LEM takes the geogrid as reinforcement where its slip surfaces cross it.

Typical values for the blocks and their contacts:

| Quantity | Entry | Typical values |
|---|---|---|
| Block size | the polygons' height and depth | Units 0.33 to 1.25 ft (0.10 to 0.375 m) high and 0.67 to 2 ft (0.20 to 0.60 m) deep ([GEC 11][gec11v1] Vol. I p. 3-43). |
| Batter | a setback at each course | From near vertical up to 15 degrees (Vol. I p. 2-40). |
| Course joints | `c` = 0, `phi` from the manufacturer's inter-unit shear tests, `t_cut` 0 | GEC 11 requires the inter-unit shear capacity from tests on the unit (ASTM D6916) (Vol. I p. 4-58) and gives no typical value. |
| Base and back face | `c` = 0, `phi` the friction of concrete on the soil there, `t_cut` 0 | As for a [gravity wall](#gravity-or-cantilever-wall) below. |
| `kn`, `ks`, `Jred` | blank: the stiffnesses come from the softer adjacent material, and the strength reduction weakens the contacts with the soil | — |

## Gravity or Cantilever Wall

A concrete gravity wall stands by its weight, and a cantilever wall by its weight and the soil on its heel. Either
is drawn as a polygon of concrete in a material whose strength option is `elastic`, so that it cannot fail. In the
LEM a trial slip surface may run around it or along its base but not through it, so a non-circular surface along
the base represents sliding on it. In the FEM the wall needs a joint line under its base and one up its back face, so that it can slide and
part from the soil.

![A concrete gravity wall twice: for the LEM, a trial surface along its base and up through the backfill; for the FEM, a joint line under its base and one up its back face](images/pw_gravity_wall.png){width=758}

On the left, the LEM's non-circular trial surface runs along the wall's base and up through the backfill behind it.
On the right, the FEM's base joint and back-face joint.

| Quantity | Entry | Typical values |
|---|---|---|
| Unit weight | γ of the concrete material | 150 lb/ft³ (23.6 kN/m³) for reinforced concrete ([GEC 4][gec4] p. A-13). |
| `E` (FEM) | modulus of the concrete material | E<sub>c</sub> = 1,820 √f′<sub>c</sub> ksi: 4.9 × 10⁸ to 5.9 × 10⁸ psf (2.35 × 10⁷ to 2.81 × 10⁷ kPa) for f′<sub>c</sub> = 3.5 to 5.0 ksi ([GEC 10][gec10] p. 16-3). |
| Base joint (FEM) | `c` = 0, since the sources give friction only; `phi` = δ, the friction angle of concrete on the foundation; `t_cut` 0 | Mass concrete on clean sound rock 35°; on clean gravel or coarse sand 29 to 31°; on clean fine to medium sand 24 to 29°; on fine sand, silty or clayey 19 to 24°; on very stiff to hard clay 22 to 26°; on medium stiff to stiff clay 17 to 19° ([GEC 10][gec10] p. 12-50). |
| Back-face joint (FEM) | `c` = 0, `phi` = δ of the wall against the backfill, `t_cut` 0 | About two thirds of the backfill's φ′ in GEC 10's example: 24° against sand at 36°, 19° against clay at 28° (pp. 12-52, 12-53). |


## Load-Bearing Piles Near a Slope

Load-bearing piles carry structural loads (vertical forces from foundations) and transfer them to the subsurface through a combination of **skin friction** along the pile shaft and **end bearing** at the pile tip. The key question for slope stability is: does the structural load contribute to the driving forces on the failure surface?

![A footing on a slope's crest on a pile, twice: in Case 1 the pile's tip is above the failure surface, in Case 2 the pile crosses it and its tip is in stable ground; in plan, footings of width B along the crest at spacing s](images/pw_load_bearing.png){width=871}

In Case 1 the pile ends inside the sliding mass, so the pile and its load move with it. In Case 2 the pile reaches
stable ground below the failure surface. In plan, the footings stand along the crest, the line with ticks hanging
down the slope, at spacing *s* center to center, each *B* wide across the crest, with its pile dashed below it.

### Case 1: Pile tip above the failure surface

If the pile tip is entirely within the sliding mass (a friction pile in weak soil, for example), the entire pile and its load are part of the sliding mass. The structural load **does** contribute to driving forces and should be included in the analysis. Standard practice is to apply the structural load as a **distributed surface surcharge** using the distributed loads (`dloads`) sheet in the XSLOPE input template. For a row of footings along the crest, each *B* wide across the crest and carrying a force *P*, at spacing *s* along it, the surcharge is *P* ÷ (*B* × *s*). It is entered in the `Normal` column at two points on the ground at the footing's edges, with Direction `vertical` ([Worksheet: dloads](../usage/input_template.md#worksheet-dloads)). This is slightly conservative because it places all the weight at the surface rather than distributing it with depth through skin friction, but the conservatism is generally small and accepted in practice.

### Case 2: Pile tip below the failure surface

This is the usual design intent for load-bearing piles near slopes — the pile is embedded in stable ground below the failure surface. The pile shaft necessarily passes **through** the sliding mass to reach that stable ground, and skin friction is mobilized along the full shaft length, both above and below the failure surface. The portion of the structural load transferred via skin friction **above** the failure surface loads the sliding mass; the remainder (skin friction below the failure surface plus end bearing) bypasses it.

In principle, determining the split requires a load-transfer analysis (t-z, or load-transfer, curves or similar). In practice, this is rarely done in the context of slope stability because the complexity is not justified. Instead, two bounding assumptions are used:

**Lower bound (omit the load)**: Assume the pile delivers all of its load to stable ground below the failure surface. The structural load is omitted entirely from the slope stability model. This is the approach recommended by FHWA and the American Association of State Highway and Transportation Officials (AASHTO), and used by commercial slope stability software (SLOPE/W, Slide2). It is appropriate when:

- The pile is designed as an end-bearing pile in competent material (rock, dense sand) — most of the load genuinely reaches the tip
- The skin friction above the failure surface is small relative to the total pile capacity (shallow failure surface relative to the pile length, or weak soil in the sliding mass)
- The structural load is modest relative to the soil driving forces

**Upper bound (full surcharge)**: Treat the full structural load as a surface surcharge, as in Case 1. This is conservative — it assumes all of the load enters the sliding mass, ignoring the load that bypasses via end bearing and deep skin friction. This approach is appropriate when:

- A significant portion of the pile shaft is above the failure surface
- The soil above the failure surface has high skin friction capacity (the pile sheds substantial load before reaching the failure surface)
- The structural load is large relative to the soil driving forces, and the lower-bound assumption would meaningfully affect the computed factor of safety

For most practical cases with end-bearing piles through a shallow sliding mass, the lower-bound (omit) approach is standard and the error is small. When in doubt, run both assumptions to bracket the answer.

### Summary

For load-bearing piles near slopes, the recommended approach in XSLOPE is:

1. If the pile tip is above the failure surface, apply the structural load as a distributed surface load on the `dloads` sheet
2. If the pile tip is below the failure surface, omit the structural load from the slope stability model (lower bound). If the structural load is significant, also run with full surcharge (upper bound) to bracket the result.
3. If the pile also provides lateral resistance to sliding, model that separately as a stabilizing pile force $H$

The distributed loads in XSLOPE handle the surcharge case, so load-bearing piles need no additional input.

A finite element run does not remove the need for these bounds: its pile is bonded to the soil with no shaft interface ([FEM piles](fem.md#pile-soil-interface-and-load-transfer)).

## Connecting Reinforcement to a Wall or Facing

A support that bears on a wall or a facing is connected at end 1. The two analyses treat the connection
differently: the LEM applies each member's force on its own, and the FEM connects a bar to another member only at
a node the two share or, for a jointed sheet, through a tie at its end.

### A tieback and a soldier-pile or sheet-pile wall

The wall is entered on the `piles` sheet ([Piles in LEM](lem.md)) and each tieback on
the `reinforce` sheet, with end 1 on the wall face. The wall's resistance is its shear force `H`, entered as under
[Soldier-Pile Wall](#soldier-pile-wall). The `piles` sheet has its own `Appl`, read as on the
`reinforce` sheet; [Tutorial LEM-9](../tutorials/lem09_tieback_wall.md) enters its stated `H` with Appl Active. The
figure shows one tieback through a soldier-pile wall, with a trial surface that crosses both.

![A tieback through a wall, with a trial surface from the excavation corner crossing its unbonded length, and the forces Tmax on the tendon and H on the wall](../usage/images/mr_connect_tieback_wall_lem.png){width=450}

The trial surface starts at the corner of the excavation, where it passes through the pile line and the wall's
`H` acts, and crosses the tieback on its unbonded length, where `Tmax` acts along the tendon.

A trial surface that passes below the toe of the wall receives no force from the wall.

In the FEM the wall is a row of beam elements
([Piles in FEM](fem.md)) and a tieback is a row of bar elements.
A bar that ends on a pile line, or crosses it, shares a node with the pile at that point, so the tieback pulls on
the wall at that node. The figure shows one tieback whose end 1 lies on the pile line below the pile's head.

![A tieback whose first node is a node of the wall's pile line below its head, with the bar's other nodes running back into the soil](../usage/images/mr_connect_tieback_wall_fem.png){width=447}

End 1 of the tieback is a node of the pile, and the bar's other nodes run back into the soil.

A bar that crosses a pile line, rather than ending on it, draws a note in the
[model checks](../studio/analysis.md#model-checks-before-a-run) that the two are joined at the crossing; a bar
meant to pass the pile unconnected ends short of it. A tieback entered as in
[Tutorial LEM-9](../tutorials/lem09_tieback_wall.md), starting on the wall face at x = 0 with the pile line 0.5 ft
behind it, crosses the pile line and is joined to it there in the FEM; in the LEM the offset has no effect.

### A soil nail and a shotcrete facing

A nail's head plate bears on the shotcrete facing at end 1, where `Tend1` sets the head's capacity in both
analyses ([Soil Nail](../reinforcement/types.md#soil-nail)). The LEM and the FEM need the facing itself entered
differently, as below. Both apply every line load and polygon in a workbook, so a wall analyzed both ways needs two
workbooks, one for each analysis.

#### In the LEM

The facing is not drawn. The cut face is the ground surface, each nail's end 1 is on it, and the facing's weight is
a vertical line load on the crest, at the crest's elevation and a few inches (about 0.1 m) behind the top of the
face ([Worksheet: lloads](../usage/input_template.md#worksheet-lloads)). The nails carry the facing
([GEC 7][gec7] p. 200), so its weight bears on the reinforced soil behind it, which is the soil that slides. On a
vertical face, one force at the top gives the same force and moment as the weight spread down the face, provided
the slip surface exits at the toe; on a battered face it is approximate. Verification problems
[VP47](../verification/rocscience.md#vp47) and [VP48](../verification/rocscience.md#vp48) enter their facings this
way, with line loads of 14.6 kN/m and 13.2 kN/m.

![A soil nail on a cut face, with the facing a dashed outline not in the model, and the facing's weight as a line load on the crest just behind the top of the face](../usage/images/mr_connect_nail_facing.png){width=478}

Drawing the facing as a polygon of concrete, as in the FEM, does not work here. A trial surface cannot cut through
a material whose strength option is `elastic`, so to reach the toe it passes under the facing, where it runs nearly
level and the facing's weight adds almost nothing to the force driving the slide. A circle can pass under the
facing only as a long, deep arc, so a circular search misses the shallow surfaces that govern.

#### In the FEM

The facing is a polygon of shotcrete, drawn on the air side of the cut face from the toe to the crest, as thick as
the whole facing, so that its outer face becomes the ground surface. Its material's strength
option is `elastic`, as for a [gravity wall](#gravity-or-cantilever-wall). Each nail is extended to end 1 inside
the polygon, at its head plate, and there is no line load. The FEM applies a line load as a force at one node
([Loads](../fem/overview.md#distributed-loads)), and a facing's whole weight on one node can fail the soil under
that node before the slope fails. The polygon spreads the weight down the face and gives the facing the stiffness
of concrete.

![The same soil nail and facing for the FEM: the facing a polygon of shotcrete with the nail's end 1 inside it, and no line load](../usage/images/mr_connect_nail_facing_fem.png){width=466}

A facing is a foot thick or less, so the mesher's thin-zone refinement sizes it for about four elements across
([Thin material zones](../fem/mesh.md#thin-material-zones)). That refinement is never finer than a sixth of the
global element size, so on a large model the facing gets fewer elements across and the model checks report it.
Give its polygon a Size of a quarter of its thickness then; a Size has no such limit
([A Size on one zone](../fem/mesh.md#a-size-on-one-zone)).

#### Typical values

| Quantity | Entry | Typical values |
|---|---|---|
| Facing thickness | the initial and final facings together: in the line load (LEM), and as the polygon's thickness normal to the face (FEM) | Initial facing of shotcrete: typically 4 in. (0.33 ft, 0.10 m), or 6 in. (0.5 ft, 0.15 m). Final facing of shotcrete: 6 or 8 in. (0.5 or 0.67 ft, 0.15 or 0.20 m), typically 8. Final facing cast in place: 10 in. (0.83 ft, 0.25 m) and thicker ([GEC 7][gec7] p. 164). |
| Unit weight | in the line load (LEM); γ of the shotcrete material (FEM) | 150 lb/ft³ (23.6 kN/m³), as for reinforced concrete ([GEC 4][gec4] p. A-13). |
| Line load (LEM) | `P` on the `lloads` sheet, with `Angle` blank (straight down) | Unit weight × facing thickness × face height: for a 4 in. initial and an 8 in. final facing, 1 ft in all, on a 20 ft face, 150 lb/ft³ × 1 ft × 20 ft = 3,000 lb/ft (43.8 kN/m). |
| `E` (FEM) | modulus of the shotcrete material | E<sub>c</sub> = 1,820 √f′<sub>c</sub> ksi, with f′<sub>c</sub> in ksi ([GEC 10][gec10] p. 16-3): 3,640 ksi, or 5.24 × 10⁸ psf (2.51 × 10⁷ kPa), for f′<sub>c</sub> = 4 ksi (4,000 psi). Shotcrete is typically 3,000 to 4,000 psi, more commonly 4,000 ([GEC 7][gec7] p. 164). |

### A geosynthetic and facing blocks, panels or a wrapped face

A layer connected to facing blocks or panels starts at the back of the facing, end 1. The figure shows one layer
connected to a block facing.

![A geosynthetic layer running from the back of a block facing, end 1, into the fill, end 2](../usage/images/mr_connect_geosynthetic_facing.png){width=455}

End 1 is on the back of the block facing, where `Tend1` is the connection, and end 2 is in the fill.

FHWA takes the long-term strength of the connection, per unit width of the layer, from connection tests on the
facing unit and the geosynthetic:

>$T_{alc} = \dfrac{T_{ult} \times CR_{cr}}{RF_D}$

where T<sub>ult</sub> is the layer's ultimate tensile strength, CR<sub>cr</sub> the fraction of it the connection
keeps over the long term, measured in those tests, and RF<sub>D</sub> the reduction factor for chemical and
biological degradation ([GEC 11][gec11] Eq. 4-41, p. B-13). T<sub>alc</sub> rises with the normal pressure on the
connection: in Example E1 it runs
from 533 lb/ft (7.8 kN/m) at the top layer to 2,550 lb/ft (37.2 kN/m) at the bottom, against the layers' long-term strength T<sub>al</sub> ([Geosynthetic Layer](../reinforcement/types.md#geosynthetic-layer)) = 1,085 and
2,169 lb/ft (15.8 and 31.7 kN/m) for the two grades the wall uses, GG-I and GG-II (Table E1-7.3, p. E1-15;
connection strengths Table E1-7.6, p. E1-18), so on eight of the wall's eleven layers the connection limits the
force at the face.

In the LEM, and in the FEM for a bonded layer, `Tend1` only raises the capacity at end 1
([Capacity Envelope](../reinforcement/lem.md#capacity-envelope)); a jointed layer's `Tend1` ties its end to the
facing at end 1, up to that capacity ([Ends, ties and the bar](../reinforcement/fem.md#ends-ties-and-the-bar)).
In the block wall of [Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md#part-2-the-same-wall-with-geogrid)
the back face of the facing is a `joints`-sheet line, and a bonded line may not end on one, so all three layers
are jointed, each tied to the blocks at `Tend1` = 40 kN/m. A wrapped face has no facing unit: the sheet is folded
back over the face and buried under the next lift
([When a Model Needs a Joint](../fem/joints.md#when-a-model-needs-a-joint)). Its line starts at the face, and
`Tend1` is the pullout resistance of that folded-back return.

## References

**GEC 4:** Sabatini, P.J., Pass, D.G., & Bachus, R.C. (1999). *Geotechnical Engineering Circular No. 4: Ground
Anchors and Anchored Systems*. [FHWA-IF-99-015](https://www.fhwa.dot.gov/engineering/geotech/pubs/if99015.pdf).
Federal Highway Administration, Washington, D.C.

**GEC 7:** Lazarte, C.A., Robinson, H., Gómez, J.E., Baxter, A., Cadden, A., & Berg, R. (2015). *Geotechnical
Engineering Circular No. 7: Soil Nail Walls Reference Manual*.
[FHWA-NHI-14-007](https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi14007.pdf). Federal Highway Administration,
Washington, D.C.

**GEC 10:** Brown, D.A., Turner, J.P., & Castelli, R.J. (2010). *Drilled Shafts: Construction Procedures and LRFD
Design Methods* (Geotechnical Engineering Circular No. 10).
[FHWA-NHI-10-016](https://rosap.ntl.bts.gov/view/dot/40746/dot_40746_DS1.pdf). Federal Highway Administration,
Washington, D.C.

**GEC 11:** Berg, R.R., Christopher, B.R., & Samtani, N.C. (2009). *Design of Mechanically Stabilized Earth Walls
and Reinforced Soil Slopes*, Volume I, [FHWA-NHI-10-024](https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi10024/nhi10024.pdf),
and Volume II, [FHWA-NHI-10-025](https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi10025/nhi10025.pdf)
(Geotechnical Engineering Circular No. 11). Federal Highway Administration, Washington, D.C.

**Micropile manual:** Sabatini, P.J., Tanyu, B., Armour, T., Groneck, P., & Keeley, J. (2005). *Micropile Design and
Construction Reference Manual*. [FHWA-NHI-05-039](https://rosap.ntl.bts.gov/view/dot/50231/dot_50231_DS1.pdf).
Federal Highway Administration, Washington, D.C.

**Nucor Skyline:** Nucor Skyline (2026). *Technical Product Manual*.
[nucorskyline.com](https://www.nucorskyline.com/file%20library/document%20library/english/brochures/product_manual_en.pdf).

[gec4]: https://www.fhwa.dot.gov/engineering/geotech/pubs/if99015.pdf
[gec7]: https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi14007.pdf
[gec10]: https://rosap.ntl.bts.gov/view/dot/40746/dot_40746_DS1.pdf
[gec11v1]: https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi10024/nhi10024.pdf
[gec11]: https://www.fhwa.dot.gov/engineering/geotech/pubs/nhi10025/nhi10025.pdf
[mp]: https://rosap.ntl.bts.gov/view/dot/50231/dot_50231_DS1.pdf
[nucor]: https://www.nucorskyline.com/file%20library/document%20library/english/brochures/product_manual_en.pdf
