# Rocscience RS2 Joint Corpus

The [RS2 Joint Verification Manual](https://www.rocscience.com/help/rs2/verification-theory/verification-manuals)
(Rocscience) publishes 23 problems on jointed rock: block and flexural toppling, plane failure,
Alejano's sliding and ploughing slabs, step-path failure, a Voronoi-tessellated mass, a jointed
tunnel, and two shear-box problems that exercise the joint's own constitutive law rather than a
slope. Every one of them is a mesh split along discontinuities with interface elements carrying
the traction between the faces, which is what XSLOPE's [joint lines](../fem/reinforcement.md#joints-without-reinforcement)
are. The rows below verify XSLOPE's FEM/**SSRM** solver on that corpus.

The wall and embankment rows that use the same element but reach it through the reinforce
sheet's `Joint` column are on the [RS2 corpus page](rs2.md) — [RS2-24](rs2.md#rs2-24) and
[RS2-48–55](rs2.md#rs2-48). Full bibliographic details for the author-year citations here are on
the shared [References](references.md) page.

## Methodology

**Where the inputs come from.** Geometry, materials, joint properties, restraints and loads are
read from the vendor's own `.fez` model rather than from the manual's tables. The manual's tables
carry errata the models do not — problem 7's slope angle is printed as 5° where its figure and its
model are 55°, and problem 19's joint inclination is printed as 59° where both give 56° — and the
model is what RS2 solved. The manual supplies what the tables are for: the referee each problem is
scored against, and RS2's own two reported factors.

**Units.** The vendor models are stated in MPa with unit weight in MN/m³. These files carry the
metric kPa / kN/m³ the rest of the corpus uses, so every stress is multiplied by a thousand.

**The referee.** Each problem is scored against its own reference solution — a closed form
(Goodman & Bray, Alejano et al., Lorig & Varona) where the problem has one, and UDEC where that is
all the manual publishes. Both are strength-reduction or limit-equilibrium answers on the same
mechanism, so the pairing is like for like. RS2's factors are recorded beside every row and score
none of them, for the reason the next paragraph gives.

**RS2's two factors.** The manual reports each problem twice, *without* and *with* the "joint
improvement" option. The vendor's rerun models show what that option is: `improve_joint_convergence`
turns on, the convergence criterion and the SRF convergence type change, `CoupledSSR` turns **off**
— so the joints are no longer reduced with the rock — and the iteration limit and tolerance go back
to their defaults. It is a different analysis, not a better-converged one, and on problem 4 it
raises the factor while on problems 17, 18 and 19 it lowers it. XSLOPE reduces the joint with the
rock, which is the *without* setting, and both of RS2's numbers are recorded so a reader can see
the spread.

**The budget.** A joint reaches equilibrium by growing slip, so a jointed model settles over tens
of thousands of viscoplastic sweeps where a bonded one settles over hundreds; a strength-reduction
trial that runs out of sweeps is recorded undecided, the bracket reads that as not standing, and
the factor comes out low. Every row here states its own sweep budget and its trial record is
checked (`tools/ssrm_trial_audit.py`), the same rule the [geotextile wall family](rs2.md#rs2-48)
runs under.

`benchmarks/rocscience/build_joint_problems.py` writes the input files, which are named `rjNNN`
by the manual's own problem number — `rj018.xlsx` is problem 18 — with a letter suffix where the
manual letters its cases. `make_rs2_joint_figures.py` writes the figures.

## Status

Status terms follow the [shared definitions](index.md#status-terms) and match dots the
[shared scoring](index.md#how-the-match-dots-are-scored). Of the manual's 23 problems, 21 report a
factor of safety or a tilt angle; problems 22 and 23 report neither, being shear-box tests of the
joint model whose output is a stress-displacement curve.

<div class="corpus-summary match" markdown>

| # | Match | Problem | Referee | RS2 without / with improvement | Notes |
|---:|:-:|---|---|---|---|
| 1a | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1a | Goodman & Bray 1.0 · UDEC 0.99 | 0.99 / 0.97 | *blocked* — the stepped base puts each column's basal contact partway along its neighbour's side joint; see [where a joint ends on another](#joint-terminations). |
| 1b | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1b | Goodman & Bray 1.0 · UDEC 0.99 | 0.97 / 0.94 | *blocked* — as 1a, and its 2013 kN toe force is a concentrated line load, which the loads sheet does not carry. |
| 1c | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1c | Goodman & Bray 1.02 · UDEC 1.01 | 1.01 / 0.99 | *blocked* — as 1a. |
| 1d | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1d | Goodman & Bray 1.23 · UDEC 1.22 | 1.19 / 1.16 | *blocked* — as 1b. |
| 2 | <span class="nodata">⊘</span> | Alejano & Alonso block toppling | UDEC 0.87 · Goodman 0.76 | 0.86 / 0.82 | *planned* — a 64° parallel set at 1.6 m over a 30° stepped basal surface, 517 joint elements. |
| 3 | <span class="nodata">⊘</span> | Lorig & Varona forward block toppling | UDEC 1.13 | 1.12 / 1.09 | *planned* — two parallel sets, 789 joint elements. |
| 4 | <span class="nodata">⊘</span> | Lorig & Varona flexural toppling | UDEC 1.3 | 1.19 / 1.27 | *planned* — one 70° set, 489 joint elements. |
| 5 | <span class="nodata">⊘</span> | Lorig & Varona backward block toppling | UDEC 1.7 | 1.65 / 1.86 | *planned* — two crossing sets, 1748 joint elements. |
| 6 | <span class="nodata">⊘</span> | Plane failure, daylighting | UDEC 1.27 | 1.25 / 1.31 | *planned* — one −35° set, 1690 joint elements. |
| 7 | <span class="nodata">⊘</span> | Plane failure, non-daylighting | UDEC 1.5 | 1.57 / 1.59 | *planned* — one −70° set, 594 joint elements. |
| 8 | <span class="nodata">⊘</span> | Flexural toppling, base friction model | UDEC 0.76 | 0.75 / 0.75 | *planned* — a −60° set at 5.08 m with a horizontal basal joint. |
| 9 | <span class="nodata">⊘</span> | Bilinear slab failure, example 1a | UDEC 1.03 | 1.01 / 1.09 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 10 | <span class="nodata">⊘</span> | Bilinear slab failure, example 1b | UDEC 1.03 | 0.92 / 1.08 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 11 | <span class="nodata">⊘</span> | Ploughing sliding slab failure | UDEC 1.21 | 1.22 / 1.3 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 12 | <span class="nodata">⊘</span> | Ploughing toppling slab failure | UDEC 1.78 | 1.39 / 1.75 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 13 | <span class="nodata">⊘</span> | Ploughing sliding slab, example 4 | UDEC 1.0 | 1.0 / 1.05 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 14 | <span class="nodata">⊘</span> | Ploughing sliding slab, example 5 | UDEC 0.9 | 0.89 / 1.09 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 15 | <span class="nodata">⊘</span> | Partially joint-controlled footwall | Slide 1.25 · UDEC 1.6 | 1.28 / 1.42 | *planned* — bedding parallel to the face at 2 m, 2143 joint elements. |
| 16 | <span class="nodata">⊘</span> | Barla et al. tilt-table block toppling | experiment 9° · UDEC 11° | 9° / 7° | *blocked* — the problem scores the TILT ANGLE at which a block grid topples, found by rotating gravity through a staged sweep; XSLOPE's seismic coefficient tilts the load but the row needs the sweep and a toppling criterion, neither of which is a strength reduction. |
| 17 | <span class="nodata">⊘</span> | Step-path, en-echelon joints | UDEC 1.29 | 1.24 / 1.2 | *blocked* — the vendor model's second, elastic material has no boundary of its own. Recovered from the element-material map it is an 83-vertex staircase of element edges, trending vertical near x = 12 and horizontal near y = 3 with excursions of about one element either side, so the boundary is a property of the vendor's mesh rather than of its model. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | 1.01 / 1.0 | |
| [19](#rj-19) | <span class="nodata">⊘</span> | Bi-planar step-path failure | UDEC 1.46 | 1.5 / 1.41 | *reported, no lock* — two trials of the bracket reach the sweep budget without a verdict. |
| 20 | <span class="nodata">⊘</span> | Hammah & Yacoub Voronoi slope | UDEC 2.46 | 2.21 / 2.37 | *blocked* — the manual states no block size and no seed, and the vendor model carries the tessellation as 523 digitized traces rather than as a generated network, so the input is not reproducible from anything published. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | 8.27 / 8.5 | *planned* — a two-stage model whose second stage excavates a 2 m opening. |
| 22 | <span class="nodata">⊘</span> | Joint model: hyperbolic softening | — | — | *not supported* — the problem exercises RS2's hyperbolic displacement- and work-softening joint law, which XSLOPE's interface element does not have; it reports no factor of safety. |
| 23 | <span class="nodata">⊘</span> | Joint model: residual strength and dilation | — | — | *no lock possible* — the problem reports no factor of safety, and its six vendor models all carry `include_dilation: no`. See [The dilation problem](#the-dilation-problem). |

</div>

---

## Where a Joint Ends on Another Joint {#joint-terminations}

Ten of the manual's twenty-one scorable problems are held out of the corpus by one property of
their geometry: a joint that stops **on** another joint rather than crossing it.

A joint is an interface, so the mesh has to carry two elements on every edge of it — one on each
face. Where two joints cross, the shared node sits at the middle of four wedges of rock and the
split copies it once per wedge, which leaves two elements on each of the four edges leaving it.
Where one joint *terminates* on another there are three wedges, not four, and the edge at the
junction comes out carrying four elements. XSLOPE refuses that mesh by name rather than solving
a section it has meshed wrongly.

Problem 1's four cases meet it exactly. Goodman & Bray's sixteen columns stand on a **stepped**
base, so each column's basal contact begins partway along its downslope neighbour's side joint:
fifteen terminations in one section, at every mesh size from 1.0 m to 4.0 m.

Problems 9 to 14 meet it by a rounding. Each of those models runs a short release trace from the
crest down to a bedding plane, and the vendor states the trace's lower tip to six decimal places,
so it lands between 2 × 10⁻⁷ and 2 × 10⁻⁶ from the bedding trace it belongs on — a part in 10⁷ of
the section. What is left between the two lines is a sliver no mesh can resolve: on problems 9 and
10 both elements on an edge come out on the same side of the joint, one of them with its centroid
3 × 10⁻⁸ from the line; on problem 11 the edge again carries four elements; on problem 12 gmsh
does not finish at all, reporting `Impossible to recover edge 113 113` and splitting the offending
curves level after level. Snapping the tip onto the bedding plane would make the near termination
an exact one, which is problem 1's case, so both families are waiting on the same thing.

Every junction in the rows that **are** built is a crossing.

## The Dilation Problem {#the-dilation-problem}

Problem 23 is the manual's statement of the joint constitutive law XSLOPE implements: a
Mohr-Coulomb interface whose limit falls from $S_{max} = c - \sigma_n \tan\phi$ to
$S_{res} = c_{res} - \sigma_n \tan\phi_{res}$ once it has slipped, with a dilation angle that opens
the joint as it slides. It is a two-block shear box — a single horizontal joint, a normal load
raised from 3 to 9 MPa part way through, and a prescribed shear displacement — and it reports a
stress-displacement curve rather than a factor of safety, so there is nothing here to lock.

The six vendor models are named for dilation angles of 0, 10, 20, 20 directional, 20
non-directional and 30 degrees. All six carry `include_dilation: no`, so the angle their names
state never reaches the solver and every one of them runs at zero dilation; the only compute-level
difference among them is the `directional` flag on one. Two of the three "20 degree" cases also
differ in when the normal load steps from 3 to 9 MPa — stage 22 against stage 15 — so they are not
a constant-history comparison either. The problem therefore verifies nothing about dilation that
could be scored, and XSLOPE's dilation is verified instead against the kinematic identity it
states: on a sliding interface the normal opening per unit slip is $\tan(\text{dil})$
(`test/joint_element_check.py`, row 6).

---

## The Rows

### 🟢 RJ-18: Step-path failure, continuous joints (rj018) {#rj-18}

A 45 × 20 m section of one Mohr-Coulomb rock (γ = 19.62 kN/m³, E = 20 GPa, ν = 0.3, c = 25 kPa,
φ = 25°, no tensile capacity) with a slope face rising from (17, 8.2) to (26.9, 20), cut by three
parallel joints at 36.1° that run from the face to the crest at a perpendicular spacing of
0.883 m. The joints carry c = 1 kPa, φ = 35°, k<sub>n</sub> = 10<sup>8</sup> kPa/m,
k<sub>s</sub> = 10<sup>7</sup> kPa/m and are reduced with the rock in the strength reduction.

| XSLOPE SSRM | UDEC referee | RS2 without / with improvement |
|---|---|---|
| **0.998** | 1.01 (−1.2%) | 1.01 / 1.00 |

<!-- test: file=files/rocscience/joints/rj018.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-18, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.7 m, which takes the mesh from 3 486 nodes and 48 interface
elements to 6 707 and 66 — does not move the factor at all, and every trial of both brackets
reaches a verdict. The budget is what the row needs rather than what it has to spare: the longest
trial to decide takes about 226 000 sweeps of the 250 000 allowed.

**Input file:** [rj018.xlsx](files/rocscience/joints/rj018.xlsx).

![RJ-18: step-path failure through three continuous joints (rj018) — FEM inputs, mesh, max shear strain and displacement vectors at the critical SRF. Every joint is slipping along its lower half and open along its upper, and the three slabs between them slide out together down the 36.1° path; the rock itself carries almost no plastic strain](images/RJ-18.png)

### ⊘ RJ-19: Bi-planar step-path failure (rj019) {#rj-19}

A 120 × 70 m section of one Mohr-Coulomb rock (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 10 500 kPa,
φ = 35°, tensile capacity 200 kPa) with a slope face from (30, 20) to (60, 70), cut by two
discontinuous joints with a rock bridge between them: a basal joint at 28.4° from (39.0149,
35.0248) to (63, 48), and an upper joint at 56.3° from (62, 49) to (76, 70). Both carry c = 0,
φ = 40° and the same stiffness pair as RJ-18. Referee: UDEC 1.46.

Two trials of the bracket reach 250 000 sweeps without a verdict, so the factor this variant
brackets is a statement about the budget as much as about the slope, and the row is reported
without a lock.

The manual's table for this problem states one joint inclination as 59°. Its figure dimensions 56°
and 28°, and the vendor model's own endpoints give 56.3° and 28.4°, so the table is the outlier and
the model is what is built.

**Input file:** [rj019.xlsx](files/rocscience/joints/rj019.xlsx).

![RJ-19: bi-planar step-path failure with a rock bridge (rj019) — FEM inputs, mesh, max shear strain and displacement vectors at the critical SRF. Both joints have opened, the basal one is slipping at its lower end, and the only plastic strain in the rock is the patch at the bridge between the two joint tips, where the block above has to break through to move](images/RJ-19.png)

---

## Where the Vendor Models Depart from the Manual

Each of these was found by reading the `.fez` against the manual page it belongs to, and each
changes what a faithful transcription is.

- **Problem 7's slope angle.** The table prints 5°; the figure and the model are the same 55°
  slope its siblings use.
- **Problem 19's joint inclination.** The table prints 59°; the figure and the model give 56°.
- **Problem 23's dilation.** All six "dilation" models run at zero dilation — see
  [The dilation problem](#the-dilation-problem).
- **Problem 1's toe force.** The manual describes a stabilizing force at the toe of the lowest
  block. It is 2013 kN in cases b and d — about twice that block's own weight — and 0.5 kN in cases
  a and c, which is 0.05% of it and does nothing. The vendor's own rerun models move the force from
  the toe at (−0.5, 0.866) to the upper-left block corner at (−2.5, 4.330) in all four cases, so
  the "with joint improvement" factors are not the same load case as the "without" ones.
- **Problem 20's network.** The manual calls it a Voronoi tessellation generated in UDEC and
  imported. The vendor model carries it as 523 digitized traces with no block size, density or seed
  recorded, so nothing published reproduces it.
- **Problem 16's boundary conditions.** The 0° case runs on rollers; the nine tilted cases pin
  every exterior node in both directions, and use a convergence tolerance two orders tighter.
- **Problems 3 and 5 are elastic.** Their manual pages describe Mohr-Coulomb rock and their
  material rows carry c = 0.675 MPa and φ = 43°, but both files set the plasticity specification to
  none, so the rock cannot yield and only the joints can.
- **Problems 9 to 14 use two joint strengths, not one.** The network's 1928–1962 elements carry one
  friction angle and the one or two short explicit crest traces carry another — 40°, 20° or 30°
  against the network's 30°, 25° or 20° — so a single quoted joint friction angle for these
  problems is incomplete.
- **Problem 22's units.** The joint numbers (c = 143, residual c = 76, load 345) read as kPa in a
  file declared in MPa: at MPa a 143 MPa joint under 345 MPa never slips at the prescribed 0.3 m of
  shear, and at kPa it is an ordinary rock-joint shear test.
- **The manual states no rock strength for problems 3 to 7.** Their tables carry the slope
  geometry, the joint friction angle and the rock's tensile strength, and nothing else. The models
  carry c = 0.675 MPa and φ = 43° on problems 4, 6 and 7, and on problems 3 and 5 no plasticity
  specification at all — an elastic rock, which only the joints can fail. The same tables state no
  joint cohesion either, where all five models carry 0.1 MPa.
- **Problem 8's tensile strength never binds.** The model caps the rock at 75 kPa, and the
  Mohr-Coulomb apex its own c and φ imply is c/tan φ = 74.1 kPa, so the cap sits above the strength
  the envelope already has. Problem 15 is the same: a 1000 kPa cap over a 285.6 kPa apex.
- **Problem 8's joint network is not the whole section.** Its thirteen columns are clipped to the
  block above the basal joint; generated over the section they would be sixteen, three of them
  under a plane the model has no columns below.
- **Problems 9 to 15 run on a rock a thousand times stiffer than steel.** E = 2 × 10⁸ MPa is the
  manual's own device, and it says so: the UDEC models these are scored against use rigid blocks,
  and the modulus is how RS2 reproduces one.
- **No problem from 1 to 21 states a joint residual strength or a dilation angle.** Those appear
  only in problems 22 and 23, and material residual values, where they appear, always equal the
  peak.
