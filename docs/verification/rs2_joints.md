# Rocscience RS2 Joint Corpus

The [RS2 Joint Verification Manual](https://www.rocscience.com/help/rs2/verification-theory/verification-manuals)
(Rocscience) publishes 23 problems on jointed rock: block and flexural toppling, plane failure,
Alejano's sliding and ploughing slabs, step-path failure, a Voronoi-tessellated mass, a jointed
tunnel, and two shear-box problems that exercise the joint's own constitutive law rather than a
slope. Every one of them is a mesh split along discontinuities with interface elements carrying the
traction between the faces, which is what XSLOPE's [joint lines](../fem/joints.md) are; how they are
modeled, and the reach of the element that carries them, is documented there. The rows below verify
XSLOPE's FEM/**SSRM** solver on that corpus.

The wall and embankment rows that use the same element but reach it through the reinforce sheet's
`Joint` column are on the [RS2 corpus page](rs2.md) — [RS2-24](rs2.md#rs2-24) and
[RS2-48–55](rs2.md#rs2-48). Full bibliographic details for the author-year citations here are on the
shared [References](references.md) page.

## Methodology

How the manual's problems reach this corpus:

- Geometry, materials, joint properties, restraints and loads are read from the vendor's own `.fez`
  model rather than from the manual's tables, which carry errata the models do not; the differences
  are listed under [where the vendor models depart from the
  manual](#where-the-vendor-models-depart-from-the-manual). The manual is the source of the referee
  each problem is scored against and of RS2's own two reported factors. The vendor models are stated
  in MPa with unit weight in MN/m³; these files carry the metric kPa and kN/m³ the rest of the
  corpus uses.
- **The referee.** Where a closed-form rigid-block limit equilibrium exists for a problem, that is
  what scores it, recomputed from the inputs the model carries rather than quoted from the source
  and validated first on the source's own worked examples: Goodman & Bray's iterative column
  analysis on problem 1's four cases and on problem 2, Alejano's ploughing equation on problems 11
  to 14, and Alejano's footwall equations on problem 15. Where no closed form exists the referee is
  the one the manual names, which is UDEC in every such case. Each row records both the recomputed
  value and the one its source prints, and the Slide2 factor Rocscience reports beside problem 15 is
  recorded and labeled as the vendor's own companion program.
- **The rigid-block bound.** A limit-equilibrium formula prices one mechanism and may be
  conservative, but it may not exceed what rigid-block statics admits: the two blocks each ploughing
  problem cuts out have a highest factor at which they admit any set of contact forces lying inside
  their friction cones. So the referee for problems 11 to 14 is Alejano's Eq. (7) where it sits on
  or under that bound, which is problems 12, 13 and 14, and the bound itself on problem 11, where
  Eq. (7) stands half again above it.
- **Apples to apples.** Every row carries XSLOPE's deviation from its referee beside RS2's deviation
  from the same referee. RS2's "without joint improvement" factor is the vendor's default and the
  same method family as XSLOPE's — a strength reduction on a continuum with interfaces in it — so it
  is the yardstick for how close a finite element code of this kind gets. The manual reports each
  problem twice, and the "with joint improvement" run turns `CoupledSSR` **off** so that the joints
  are no longer reduced with the rock, which makes it a different analysis rather than a
  better-converged one. XSLOPE reduces the joint with the rock, which is the *without* setting, and
  both of RS2's numbers are recorded so a reader can see the spread.
- **The mesh.** Every row is meshed at its own joint spacing — the column width, the bedding
  spacing, the mean block width of a tessellation — so that one element spans the rock between one
  discontinuity and the next, and every row states what a step of refinement does to its factor.
  Problem 1's four cases are the exception: at their block width of 10 m that rule puts about one
  element across a column, so those four are cut at a 2D size of 5 m, where three successively finer
  meshes return the same bracket.
- **The budget.** A joint reaches equilibrium by growing slip, so a jointed model settles over tens
  of thousands of viscoplastic sweeps where a bonded one settles over hundreds; a trial that runs
  out of sweeps is recorded undecided, the bracket reads that as not standing, and the factor comes
  out low. Every row here runs at a budget of 250,000 sweeps and states its own trial record, which
  `tools/ssrm_trial_audit.py` reports and `test/corrector_certified_check.py` holds to the rule
  below.
- **What answers a bracket edge.** {#what-answers-a-bracket-edge}
  A bracket asks of each trial factor whether the model stands there, and a trial that runs out of
  sweeps with nothing to say has not answered. Two readings answer one the force tolerance alone
  does not. In the standing direction, the solver offers the viscoplastic loop's own state to a
  Newton corrector at its checkpoints, and where the corrector reaches equilibrium from that state
  inside its force, yield and displacement gates it records the certification on the trial. In the
  failing direction, an interface still growing its slip at a steady rate while the displacement
  field goes nowhere is read as failing. A trial with neither reading is decided only inside its
  budget, and a corrector refusal is not a verdict of any kind: it is the absence of one.
- `benchmarks/rocscience/build_joint_problems.py` writes the input files, which are named `rjNNN` by
  the manual's own problem number — `rj018.xlsx` is problem 18 — with a letter suffix where the
  manual letters its cases. `make_rs2_joint_figures.py` writes the figures.

## Status

Status terms follow the [shared definitions](index.md#status-terms) and match dots the
[shared scoring](index.md#how-the-match-dots-are-scored). Of the manual's 23 problems, 21 report a
factor of safety or a tilt angle; problems 22 and 23 report neither, being shear-box tests of the
joint model whose output is a stress-displacement curve.

<div class="corpus-summary match match8" markdown>

| # | Match | Problem | Referee | RS2 vs referee | Also published | RS2 without / with improvement | Notes |
|---:|:-:|---|---|---|---|---|---|
| [1a](#rj-1a) | 🟢 | Goodman & Bray block toppling, case 1a | SSRM 1.018 vs Goodman & Bray 1.0000 (+1.8%) | 0.99 vs 1.0000 (−1.0%) | UDEC 0.99 (+2.8%) | 0.99 / 0.97 | |
| [1b](#rj-1b) | 🟢 | Goodman & Bray block toppling, case 1b | SSRM 1.018 vs Goodman & Bray 1.0000 (+1.8%) | 0.97 vs 1.0000 (−3.0%) | UDEC 0.99 (+2.8%) | 0.97 / 0.94 | |
| [1c](#rj-1c) | 🟢 | Goodman & Bray block toppling, case 1c | SSRM 1.037 vs Goodman & Bray 1.0185 (+1.8%) | 1.01 vs 1.0185 (−0.8%) | UDEC 1.01 (+2.7%) | 1.01 / 0.99 | |
| [1d](#rj-1d) | 🟢 | Goodman & Bray block toppling, case 1d | SSRM 1.252 vs Goodman & Bray 1.2308 (+1.7%) | 1.19 vs 1.2308 (−3.3%) | UDEC 1.22 (+2.6%) | 1.19 / 1.16 | |
| [2](#rj-2) | 🟢 | Alejano & Alonso block toppling | SSRM 0.764 vs Goodman & Bray 0.7734 (−1.2%) | 0.86 vs 0.7734 (+11.2%) | UDEC 0.87 (−12.2%) | 0.86 / 0.82 | UDEC and RS2 both stand above the closed form; the manual states no UDEC settings for this model. |
| [3](#rj-3) | <span class="nodata">⊘</span> | Lorig & Varona forward block toppling | UDEC 1.13 | 1.12 vs 1.13 (−0.9%) | — | 1.12 / 1.09 | *reported, no lock* — the failing edge of the bracket reaches the sweep budget without a verdict, and the corrector refuses it. |
| [4](#rj-4) | 🟢 | Lorig & Varona flexural toppling | SSRM 1.311 vs UDEC 1.3 (+0.8%) | 1.19 vs 1.3 (−8.5%) | — | 1.19 / 1.27 | |
| [5](#rj-5) | <span class="nodata">⊘</span> | Lorig & Varona backward block toppling | UDEC 1.7 | 1.65 vs 1.7 (−2.9%) | — | 1.65 / 1.86 | *reported, no lock* — four trials of the bracket, including both its edges, reach the sweep budget without a verdict, and the corrector certifies none of them; a finer mesh decides two more and leaves its own standing edge undecided. |
| [6](#rj-6) | <span class="nodata">⊘</span> | Plane failure, daylighting | UDEC 1.27 | 1.25 vs 1.27 (−1.6%) | — | 1.25 / 1.31 | *reported, no lock* — every trial of both brackets decides, and a step of refinement moves the factor by twice the row's tolerance. |
| [7](#rj-7) | 🟡 | Plane failure, non-daylighting | SSRM 1.564 vs UDEC 1.5 (+4.3%) | 1.57 vs 1.5 (+4.7%) | — | 1.57 / 1.59 | Both finite element codes land above the referee, on the same side and within half a point of each other. |
| [8](#rj-8) | 🟢 | Flexural toppling, base friction model | SSRM 0.764 vs UDEC 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | — | 0.75 / 0.75 | |
| [9](#rj-9) | 🟢 | Bilinear slab failure, example 1a | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | LE (Alejano) 0.40–1.45 | 1.01 / 1.09 | |
| [10](#rj-10) | 🟢 | Bilinear slab failure, example 1b | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | LE (Alejano) 0.43–1.45 | 0.92 / 1.08 | |
| [11](#rj-11) | 🟢 | Ploughing sliding slab failure | SSRM 1.213 vs rigid-block bound 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | UDEC 1.21 (+0.2%) · Alejano Eq. (7) 1.7582 | 1.22 / 1.3 | Alejano's Eq. (7) returns a factor above the bound on this problem, so the bound is what scores it; all three programs sit on the bound. |
| [12](#rj-12) | 🟡 | Ploughing toppling slab failure | SSRM 2.033 vs Alejano Eq. (7) 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | UDEC 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 | |
| [13](#rj-13) | 🟢 | Ploughing sliding slab, example 4 | SSRM 0.998 vs Alejano Eq. (7) 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | UDEC 1.0 (−0.2%) · Alejano prints 1.0 | 1.0 / 1.05 | |
| [14](#rj-14) | 🟢 | Ploughing sliding slab, example 5 | SSRM 1.232 vs Alejano Eq. (7) 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | UDEC 0.9 (+36.9%) · Alejano prints 1.00 | 0.89 / 1.09 | The 1.00 the paper prints for this example does not follow from the inputs it prints; Eq. (7) on them gives 1.2034. |
| [15](#rj-15) | 🔴 | Partially joint-controlled footwall | SSRM 1.271 vs Alejano Eqs. (9)–(10) 1.7985 (−29.3%) | 1.28 vs 1.7985 (−28.8%) | UDEC 1.6 (−20.6%) · Alejano prints 1.72 · Slide2 (vendor) 1.25 | 1.28 / 1.42 | The closed form drives a wedge out through a single 2 m bed and its factor rises with bed thickness; the three programs free to search for a surface agree at 1.25–1.28. |
| [16](#rj-16) | <span class="nodata">⊘</span> | Barla et al. tilt-table block toppling | UDEC 11° | 9° vs 11° (−18.2%) | Experiment 9° | 9° / 7° | *reported, no lock* — the problem scores a tilt angle rather than a factor of safety. Swept as a seismic coefficient at full strength, the stack stands at 8.40° and topples at 8.55°. |
| [17](#rj-17) | 🟡 | Step-path, en-echelon joints | SSRM 1.213 vs UDEC 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | — | 1.24 / 1.2 | No closed form: three rock bridges decide the factor, and both finite element codes read them below the distinct-element run, 2.1 points apart. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | 1.01 vs 1.01 (0.0%) | — | 1.01 / 1.0 | |
| [19](#rj-19) | 🔴 | Bi-planar step-path failure | SSRM 1.623 vs UDEC 1.46 (+11.2%) | 1.5 vs 1.46 (+2.7%) | — | 1.5 / 1.41 | The factor is set by the tensile cap on the rock bridge rather than by its cohesion, and the cap is reduced with the trial factor here as the vendor reduces it. |
| [20](#rj-20) | <span class="nodata">⊘</span> | Hammah & Yacoub Voronoi slope | UDEC 2.46 | 2.21 vs 2.46 (−10.2%) | — | 2.21 / 2.37 | *reported, no lock* — four trials of the bracket reach the sweep budget without a verdict, both of the edges it closes on among them, and the corrector certifies none of them. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | 8.27 vs 8.16 (+1.3%) | — | 8.27 / 8.5 | *blocked* — the second stage of the vendor model excavates a 2 m opening and the strength reduction runs on the excavated state, which carries the stress the first stage left behind; XSLOPE has no staged construction. |
| 22 | <span class="nodata">⊘</span> | Joint model: hyperbolic softening | — | — | — | — | *not supported* — the problem exercises RS2's hyperbolic displacement- and work-softening joint law, which XSLOPE's interface element does not have; it reports no factor of safety. |
| 23 | <span class="nodata">⊘</span> | Joint model: residual strength and dilation | — | — | — | — | *no lock possible* — the problem reports no factor of safety, and its six vendor models all carry `include_dilation: no`. See [The dilation problem](#the-dilation-problem). |

</div>

---

## The Rows

### 🟢 RJ-1a: Goodman & Bray block toppling, case a (rj001a) {#rj-1a}

Goodman & Bray's own toppling example, and the section all four of problem 1's cases share: sixteen
rock columns 10 m wide and 4 to 40 m tall standing on a base that steps up one metre per column at
30°, their sides at 120° — normal to that base — and their tops cut off by a 56.6° face above column
ten. The toe is at (−0.5, 0.866) and the crest at (130.564, 93.856), inside the 261 × 140 m
rectangle the vendor cuts, with everything above the columns deleted. The rock is elastic
(γ = 25 kN/m³, E = 20 GPa, ν = 0.3), the vendor's own material type, so the columns cannot break and
every mechanism the model has is a joint one — which is the idealization the closed form makes.

The joints are not a list. They are the **shared edges** of the seventeen polygons the vendor's
sixteen column outlines and the rock beneath them make: 31 contacts over 460.0 m, merged into 31
straight lines, so a pair of columns that touches over three metres gets a three-metre joint rather
than a full-height one. Fifteen of those contacts **end on** another, because a stepped base begins
each column's basal contact partway along its downslope neighbour's side joint; see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint). They carry no cohesion, φ = 38.15° — the angle
this case is posed at — and the corpus's standard stiffness pair.

The referee is Goodman & Bray's own iterative column analysis, recomputed on this section rather
than quoted from the manual. Block by block from the crest down, the force a column needs from the
one below it is the larger of what toppling about its downslope base corner demands and what sliding
on its base demands. On this case the recursion reproduces the method's published mode pattern —
columns 14 to 16 stable, 13 down to 3 toppling, 2 and 1 sliding — and the horizontal toe force it
requires for limit equilibrium is 0.36 kN/m against the 0.5 kN the manual states. Carried into a
strength reduction it gives the referee each of the four cases is scored against below, and each of
the four is within a rounding of the value the manual prints for it.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.99 vs 1.0000 (−1.0%) | 0.99 (+2.8%) | 0.99 / 0.97 |

<!-- test: file=files/rocscience/joints/rj001a.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1a, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

**The excess over the closed form is one assumption, and it is where the thrust between two columns
acts.** Goodman & Bray hand each column-to-column thrust to the top corner of its contact. A contact
cannot do that: the two columns lean together, the face parts from the block's base upward and stays
closed only over its upper quarter to half, and the traction there is distributed, so the resultant
stands a tenth to a sixth of the face below the corner. Given those heights read off the solved
state, and nothing else changed, the same recursion returns **1.0279** — this row's own standing
edge. The closed form's other three assumptions the solution obeys exactly: every block balances in
force and in moment on the interface tractions alone, every closed side pair is at its friction
limit, and the base reaction of every toppling block sits on its downslope corner.

**The mesh moves this row by one bisection step, and the lock stands past that step.** At the block
width of 10 m — 2,072 nodes and 67 interface elements — the bracket reads 1.037; at 2D sizes of
7.0 m, 5.0 m and 3.5 m it reads 1.018, and the three finer meshes agree with one another. The row
is locked at 5.0 m, where the bracket is decided on all nine trials: the failing edge is read as failed
on a steadily slipping interface at 225,001 sweeps of the 250,000 allowed, and four of the five
standing trials are certified by the Newton corrector 302 to 311 sweeps in. What refinement does to the interface is let each contact open a
little further, so the thrust it carries sits a little higher and implies a little less: the closed
form the measured heights imply is 1.0279 at 10 m and 1.0271 at 5.0 m, and 1.02734375 — the value
the two readings straddle — is one of the bisection's own grid points.

The manual's figure for this case labels the side boundaries as rollers. The model's restraint list
does not: all 191 of its restrained nodes carry both components fixed, 111 of them on the base and
80 on the two sides. The restraints follow the model. Its stabilizing toe force is 0.5 kN against a
toe column weighing about a thousand times that, and the file carries no load at all; see the
departures table.

**Input file:** [rj001a.xlsx](files/rocscience/joints/rj001a.xlsx).

![RJ-1a: Goodman & Bray block toppling, case a (rj001a) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own, so the whole figure is the joints: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope, and the deformed section shows the stack rotating forward over the face column by column while the toe column slides out along its own stretch of base](images/RJ-1a.png)

### 🟢 RJ-1b: Goodman & Bray block toppling, case b (rj001b) {#rj-1b}

Case a's section and rock with the joints at φ = 33.0239°, the lowest of the four, held up by the
2013 kN horizontal force the case is posed with — about twice the lowest column's own weight, and
what lets a stack on 33° joints stand where case a's needs 38°. It enters as a line load on the
`lloads` sheet, at the point the departures table names.

Nothing else about the mechanism changes: the stack still rotates forward over the face on slip up
every column contact and along the stepped base.

On this case the recursion requires a horizontal toe force of 2,012.86 kN/m for limit equilibrium,
against the 2,013 kN the manual states — 0.007% on a stack weighing 83,500 kN/m.

The thrust heights here are not case a's. The toe force pushes the two lowest columns back into
their own step risers and carries those two thrusts far down their faces, and given this case's own
measured heights the recursion returns **1.0078**, which is this row's standing edge.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.97 vs 1.0000 (−3.0%) | 0.99 (+2.8%) | 0.97 / 0.94 |

<!-- test: file=files/rocscience/joints/rj001b.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1b, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

The mesh does not move this row. At the block width of 10 m and at the 5.0 m of its lock it
returns the same factor, and at 5.0 m every trial of the bracket decides inside its own budget: the
standing edge converges under its own steam at 115,020 sweeps and the longest trial to reach a
verdict takes 177,381 of the 250,000 allowed. Case a's one-step move with refinement is not
available here, for the reason the toe force gives: on this case the factor sits where a few hundred
kilonewtons per metre of extra capacity is worth under two points of factor of safety, and the
interface change refinement makes is smaller than that.

Every input class matches the vendor model: the rock's
E, ν and γ, the joints' stiffness pair, cohesion, friction angle and tensile cap, the 31 contacts the
column outlines make, the side restraint the vendor clamps in both directions, and the force's
magnitude and direction.

**Input file:** [rj001b.xlsx](files/rocscience/joints/rj001b.xlsx).

![RJ-1b: Goodman & Bray block toppling, case b (rj001b) — FEM inputs with the 2013 kN toe force drawn at the block corner it acts on, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's — slip up every column contact and along the stepped base, the stack rotating forward over the face — reached on joints four degrees weaker, which is what the force at the toe buys](images/RJ-1b.png)

### 🟢 RJ-1c: Goodman & Bray block toppling, case c (rj001c) {#rj-1c}

Case a's section and rock with the joints half a degree steeper, φ = 38.6598°, and the same 0.5 kN
toe force that does nothing. Half a degree is worth about two points of factor of safety here, which
is what the closed form's own pair of answers for cases a and c says as well — 1.0 against 1.02.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.0185 (+1.8%) | 1.01 vs 1.0185 (−0.8%) | 1.01 (+2.7%) | 1.01 / 0.99 |

<!-- test: file=files/rocscience/joints/rj001c.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1c, tier=gate, f_stand=1.02734375, f_fail=1.046875, check=edges -->

Given the thrust heights read off the solved state at the block width of 10 m, the recursion
returns **1.0469**, that mesh's own standing edge — case a's mechanism and case a's assumption, half
a degree of friction further on.

**The mesh moves this row the same way it moves case a.** At the block width of 10 m the bracket
reads 1.057; at the 5.0 m of its lock it reads 1.037, one bisection step down and inside the
row's own tolerance. Every trial of the 5.0 m bracket decides: five stand, four of them certified by
the Newton corrector 302 to 312 sweeps in, and the failing edge is read as failed on a steadily
slipping interface at 225,001 sweeps of the 250,000 allowed.

**Input file:** [rj001c.xlsx](files/rocscience/joints/rj001c.xlsx).

![RJ-1c: Goodman & Bray block toppling, case c (rj001c) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's at half a degree more friction: the same forward rotation of the column stack over the face, and the same slip on every column contact and along the stepped base](images/RJ-1c.png)

### 🟢 RJ-1d: Goodman & Bray block toppling, case d (rj001d) {#rj-1d}

Case c's joint friction angle with case b's 2013 kN force: the same stack, stabilized. It is the
strongest of the four, and the closed form and UDEC agree that it is.

Given the thrust heights read off this case's solved state, the recursion returns **1.2422**,
which is this row's standing edge. Across all four cases XSLOPE stands above the closed form and
RS2 below it, with UDEC between them, and on each of the four the recursion lands on the row's own
standing edge once the thrust is put where the solution puts it rather than at the corner of the
contact.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.252** | 1.2308 (+1.7%) | 1.19 vs 1.2308 (−3.3%) | 1.22 (+2.6%) | 1.19 / 1.16 |

<!-- test: file=files/rocscience/joints/rj001d.xlsx, type=fem_ssrm, expected_fs=1.252, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1d, tier=gate, f_stand=1.2421875, f_fail=1.26171875, check=edges -->

The mesh does not move this row either: the block width of 10 m and the 5.0 m of its lock
return the same bracket, edge for edge, and at 5.0 m every trial decides inside its own budget — the
standing edge converging at 181,301 sweeps and the failing edge running away at 147,501. Like case
b, this case is posed where the toe-force curve is steep, so the interface change that refinement
makes cannot carry it across a bracket step.

Every input class matches the vendor model, the force reaching the same point case b's does.

**Input file:** [rj001d.xlsx](files/rocscience/joints/rj001d.xlsx).

![RJ-1d: Goodman & Bray block toppling, case d (rj001d) — FEM inputs with the 2013 kN toe force, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The right-hand panels are the state past the critical factor: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope and at the toe, and the deformed section shows the stack rotating forward over the face column by column — case b's mechanism at the strength this case's toe force holds](images/RJ-1d.png)

### 🟢 RJ-2: Alejano & Alonso block toppling (rj002) {#rj-2}

A 30 × 19.85 m section whose 9.85 m face rises at 58.65° from (20, 10) to (14, 19.85). A basal
joint runs from the toe of that face up to the crest at (2.9393, 19.85) at 30° — the stepped
surface the columns stand on — and twenty-two columns at 64°, 1.6 m apart, pass through the same
toe. Both carry φ = 31° and no cohesion. The rock is elastic (γ = 25 kN/m³, E = 20 GPa, ν = 0.3),
the vendor's own `Plasticity Specifications: Non`, so the columns cannot yield and every mechanism
the model has is a joint one.

The referee is Goodman & Bray's iterative column analysis, recomputed on this section. The
twenty-two column lines at 64° and the 30° basal joint cut **thirteen** columns out of the toppling
mass — the rest of the set lies below the basal plane — and what stands on it is the 54.47 m²
triangle above it, weighing 1,361.8 kN/m. The recursion gives **0.7734**, against the 0.76 Alejano
& Alonso publish for the same method on the same problem.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **0.764** | 0.7734 (−1.2%) | 0.86 vs 0.7734 (+11.2%) | 0.87 (−12.2%) | 0.86 / 0.82 |

<!-- test: file=files/rocscience/joints/rj002.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-2, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.35 m, which takes the mesh from 11,537 nodes to 21,110 —
moves the factor by one bracket step, inside the row's own tolerance, and all nine trials of both
brackets reach a verdict. On the corpus mesh every standing trial is certified by the Newton
corrector 313 to 350 sweeps in and the longest trial of the bracket takes 6,281 sweeps of the
250,000 allowed.

**The two programs that are not the closed form land together above it.** Alejano & Alonso publish
their own UDEC run at 0.87 beside their Goodman & Bray 0.76, and RS2 reports 0.86 and 0.82: both of
the other programs stand above the recomputed closed form, where this row sits below it. Neither of
the two can be interrogated on this problem: the manual states nothing about the model's UDEC
settings, no block rounding, no deformability and no stiffness, where the rigid-block note it gives
for problems 9 to 14 is the only statement of the kind it makes anywhere.

The joints reconstruct from the vendor's own 517 joint elements to 311.4 m of trace, 291.7 m at 64° and 19.7 m
at 30°, against this file's 311.4 m at the same angles over the same extent; the rock's E, ν and γ,
the joints' k<sub>n</sub>, k<sub>s</sub>, c, φ and tensile cap, and the clamped side restraint all
match.

The relief RS2 applies to a slipping joint's stiffness is off here as it is on every row of this
corpus — see the departures table. Turned on, it moves this row's factor UP, toward UDEC, where on
the problem 1 stacks it moves the factor down.

**Input file:** [rj002.xlsx](files/rocscience/joints/rj002.xlsx).

![RJ-2: Alejano & Alonso block toppling (rj002) — FEM inputs, mesh, joint slip at the critical SRF and the section deformed 17×. The rock carries no strain of its own because it cannot yield, so every movement in the section is on a joint: slip gathers where the basal joint reaches the toe of the face, the columns standing on it open along their upper halves, and the deformed section shows them rotating out over the face while the rock below the basal joint stays put](images/RJ-2.png)

### ⊘ RJ-3: Lorig & Varona forward block toppling (rj003) {#rj-3}

The 260 m section at 55° that problems 3 to 7 share, cut by two sets: columns at 70° at 20 m
spacing and a cross set at −20° at 30 m, both through the origin. The manual states the pair as
"70 and 160" degrees, which is the same two planes measured the other way round the half circle.
The rock is elastic — the vendor's `Plasticity Specifications: Non` — so only the joints can fail;
γ = 26.0946 kN/m³, E = 9072 MPa, ν = 0.26. The joints carry c = 100 kPa and φ = 40°.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| *no lock* | 1.13 | 1.12 vs 1.13 (−0.9%) | 1.12 / 1.09 |

The row carries no lock, so it prints no factor of its own. The FAILING edge of its bracket,
F = 1.222656, reaches 250,000 sweeps `AMBIGUOUS`, and a bracket edge nothing ruled on cannot
define a factor of safety. Every other trial decides, the longest of them taking 206,973 sweeps.

That edge is one the Newton corrector was offered and refused, as [problem 5](#rj-5)'s four are and
[problem 20](#rj-20)'s four are. A refusal is the absence of a verdict rather than a
verdict of its own — see [what answers a bracket edge](#what-answers-a-bracket-edge) — so the
trial stays undecided and the row stays reported.

Every input class matches the vendor model: the
rock's E, ν and γ, the joints' normal and shear stiffness, cohesion, friction angle and tensile
cap, the two sets' dips and spacings, and the side restraint the vendor clamps in both directions.

**Input file:** [rj003.xlsx](files/rocscience/joints/rj003.xlsx).

![RJ-3: Lorig & Varona forward block toppling (rj003) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own, so the whole mechanism is on the two sets: the steep 70° joints slip behind the crest while the −20° cross joints open along it, and the deformed section shows the blocks between them rotating forward over the face](images/RJ-3.png)

### 🟢 RJ-4: Lorig & Varona flexural toppling (rj004) {#rj-4}

The same section cut by one set of columns at 70° at 20 m spacing — problem 3's first set without
its cross joints, so the columns bend rather than topple as blocks. Here the rock is Mohr-Coulomb
and carries a tensile cutoff of zero, which is what lets a column break in flexure: γ = 26.1 kN/m³,
E = 9072 MPa, ν = 0.26, c = 675 kPa, φ = 43°.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.311** | 1.3 (+0.8%) | 1.19 vs 1.3 (−8.5%) | 1.19 / 1.27 |

<!-- test: file=files/rocscience/joints/rj004.xlsx, type=fem_ssrm, expected_fs=1.311, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-4, f_stand=1.30078125, f_fail=1.3203125, check=edges, tier=gate -->

A step of refinement — a 2D size of 8.4 m, which takes the mesh from 9,993 nodes and 936 interface
elements to 19,130 and 1,325 — does not move the factor at all: the finer mesh returns the same
bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking 186,870
sweeps of the 250,000 allowed.

Its lowest trial, F = 0.5, is the corpus's one `JOINT_SETTLED` verdict: the slip, the displacement
field and the soil residual have all stopped and what is left is a limit cycle on the joint degrees
of freedom, which no budget brings down. That is a decision and not a budget running out, and the
bracket reads it as standing.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj004.xlsx](files/rocscience/joints/rj004.xlsx).

![RJ-4: Lorig & Varona flexural toppling (rj004) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. Here the rock can yield, and it does: a band of shear strain climbs from the toe across the columns, and drawn without exaggeration the columns are bent through that band rather than rotated about it, which is what separates flexural toppling from the block toppling of problem 3](images/RJ-4.png)

### ⊘ RJ-5: Lorig & Varona backward block toppling (rj005) {#rj-5}

The shared 260 m section cut by two sets: one at −55° at 10 m spacing through the toe at
(560, 140), dipping out of the face so the blocks lean back rather than forward, and a horizontal
set at 40 m spacing. The rock is elastic — the vendor's `Plasticity Specifications: Non` —
γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26; the joints carry c = 100 kPa and φ = 40°.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| *no lock* | 1.7 | 1.65 vs 1.7 (−2.9%) | 1.65 / 1.86 |

This is the most budget-limited row in the corpus, and the one row here that nothing in the solver
has moved. **Four** of its nine corpus trials reach 250,000 sweeps without a verdict, two of them
the edges of the final bracket, so the factor the bracket encloses is conditioned on the sweep
limit at both ends. The row prints no factor and carries no lock.

**The corrector was offered every one of those four trials and certified none of them.** That is
what separates this row from the rest of the corpus: where a bracket edge elsewhere is settled
either by reading a steadily slipping interface as failing or by the corrector reaching
equilibrium from the loop's own field, neither happens here. A corrector refusal is the absence of
a verdict, not a verdict of its own, so the trials stay undecided.

A step of refinement — a 2D size of 8.4 m, which takes the mesh from 21,659 nodes and 3,056
interface elements to 29,765 and 3,514 — decides two more of the nine and brackets a factor one
step below the corpus mesh's. Two trials still do not decide there, and one of them is the
standing edge of that bracket, so the finer mesh does not settle the row either and neither
bracket prints a factor. Between them the two take about 16 hours, more than twice any other row
in this corpus.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj005.xlsx](files/rocscience/joints/rj005.xlsx).

![RJ-5: Lorig & Varona backward block toppling (rj005) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own; slip runs along the −55° joints in the wedge behind the face while the horizontal bedding opens, and the deformed section shows the slabs stepping out over one another down the face, each leaning back into the slope as it goes](images/RJ-5.png)

### ⊘ RJ-6: Plane failure with daylighting discontinuities (rj006) {#rj-6}

The same section cut by one set at −35° at 10 m spacing through the origin. The joints dip out of
the 55° face at a shallower angle than the face itself, so every one of them daylights and the
slabs between them are free to slide out. The rock is Mohr-Coulomb (γ = 26.1 kN/m³, E = 9072 MPa,
ν = 0.26, c = 675 kPa, φ = 43°, no tensile capacity).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| *no lock* | 1.27 | 1.25 vs 1.27 (−1.6%) | 1.25 / 1.31 |

The row prints no factor and carries no lock, and what holds it back is the mesh rather than the
sweep budget. Every trial of both brackets reaches a verdict — the corpus bracket's longest takes
225,001 sweeps of the 250,000 allowed, and four of its nine trials are settled by a corrector
certification — but a step of refinement to a 2D size of 8.4 m, which takes the mesh from 13,134
nodes and 1,864 interface elements to 25,232 and 2,654, moves the factor by two bracket steps.
The row's own tolerance is one.

It is the only row in the corpus whose factor moves with the mesh: settled, reproducible, and the
answer of its discretization rather than of its section. Every other row that states a refinement
step either returns the same bracket edge for edge or moves inside its own tolerance.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj006.xlsx](files/rocscience/joints/rj006.xlsx).

![RJ-6: plane failure with daylighting discontinuities (rj006) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Slip runs the full length of every joint that reaches the face, over a wedge bounded below by the joint through the toe, and the rock between them carries only a faint strain: the slabs slide out along the joints rather than breaking through anything](images/RJ-6.png)

### 🟡 RJ-7: Plane failure with non-daylighting discontinuities (rj007) {#rj-7}

The same section and the same rock as problem 6, cut by one set at −70° at 20 m spacing through
the origin. The joints now dip out of the face more steeply than the 55° face itself, so none of
them daylights: a slab cannot slide out along one without shearing rock, and the slope stands
higher than problem 6's.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.564** | 1.5 (+4.3%) | 1.57 vs 1.5 (+4.7%) | 1.57 / 1.59 |

<!-- test: file=files/rocscience/joints/rj007.xlsx, type=fem_ssrm, expected_fs=1.564, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-7, f_stand=1.5546875, f_fail=1.57421875, check=edges, tier=gate -->

A step of refinement — a 2D size of 8.4 m, which takes the mesh from 9,986 nodes and 934 interface
elements to 19,179 and 1,323 — does not move the factor at all: the finer mesh returns the same
bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking 118,821
sweeps of the 250,000 allowed.

**Sweeps are not what this row is short of.** A trial of this bracket re-solved on its own at four
times the budget — a million sweeps against 250,000 — returns `STABLE_STUCK`, the same verdict it
returns at 250,000: an undecided jointed trial of this kind is a statement about the convergence
criterion rather than about the iteration limit. What settles five of the nine trials of each
bracket here is the Newton corrector reaching equilibrium from the loop's own state.

Every transcribed input class matches the vendor model, including the side restraint. The manual's
own table for this problem prints the slope angle as 5°; the figure and the model are the same 55°
slope problem 6 uses.

**Input file:** [rj007.xlsx](files/rocscience/joints/rj007.xlsx).

![RJ-7: plane failure with non-daylighting discontinuities (rj007) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. No joint reaches the face at a shallower angle than the face itself, so the failure cannot slide out along one: a strong band of shear strain cuts across the steep joints from the toe, and the mass above it moves out over rock it has had to break](images/RJ-7.png)

### 🟢 RJ-8: Flexural toppling in a base friction model (rj008) {#rj-8}

Pritchard & Savigny's base-friction table model, scaled up a hundred times: a 72.4 × 36.5 m
section whose 30.5 m face rises at 78° from (15, 6) to (21.48, 36.5). Twelve columns at −60°,
5.08 m apart, stand on a horizontal joint at y = 6 that runs the width of the model, with a
vertical joint at x = 68.4 closing the back of the stack; the three highest columns reach that back
joint rather than the base, which is where the vendor's own twelve chains end. The column set is
bounded by those two joints rather than by the section, which makes it 301.555 m of column trace
and 389.462 m of joint in all — the vendor's own totals to the millimetre. The rock is
Mohr-Coulomb (γ = 25.506
kN/m³, E = 22.771 GPa, ν = 0.139, c = 60 kPa, φ = 39°). The joints carry no cohesion, φ = 39°, and
the softest normal stiffness in the corpus bar one: k<sub>n</sub> = 1.5 × 10<sup>7</sup> kPa/m
against the set's usual 10<sup>8</sup>.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **0.764** | 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | 0.75 / 0.75 |

<!-- test: file=files/rocscience/joints/rj008.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-8, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 5,389 nodes to 10,684 —
moves the factor by one bracket step, inside the row's own tolerance, and all nine trials of both
brackets reach a verdict. The budget is not what settles this row: its longest trial takes 1,221
sweeps of the 250,000 allowed, where the [geotextile wall family](rs2.md#rs2-48) exhausts that
budget on five trials of eight rows.

The two meshes disagree about one trial, and it is the one the bracket closes on. At F = 0.753906
the corpus mesh reaches equilibrium in 312 sweeps and the finer mesh runs away in 1,981, so the
corpus bracket closes a step above the refined one. The row is cut on the corpus mesh, as every row
here is, with the refinement inside its tolerance.

The model's stated tensile strength of 75 kPa is above the Mohr-Coulomb apex its own c and φ imply
(c/tan φ = 74.1 kPa), so the cap never binds and the envelope's own apex governs.

**Input file:** [rj008.xlsx](files/rocscience/joints/rj008.xlsx).

![RJ-8: flexural toppling in a base friction model (rj008) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The strain gathers into one lobe per column along a band that climbs from the toe across the stack; above that band the columns are visibly bent rather than merely tilted, and the rock below it is unstrained. That is the break surface of flexural toppling rather than sliding along any one joint](images/RJ-8.png)

### 🟢 RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) {#rj-9}

A 50 m slope at 50° cut by bedding dipping **out of the face** at −50° at 3 m spacing, φ = 30°, with
a two-segment release trace at the toe at φ = 40° that undercuts the lowest slab. Two joint
strengths, which the manual's own geometry table prints in the reverse order from its RS2 legend.
The rock is elastic at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3 — not a rock modulus but the manual's
device for reproducing UDEC's rigid blocks, which it says outright.

The face and the bedding dip at the same angle, so no slab can slide out along a single plane: the
mechanism is the bilinear one the problem is named for, sliding on a basal plane combined with
sliding along the release that the face undercuts. The release trace's lower tip is stated to six
decimals, so it lands 2 × 10⁻⁷ from the bedding plane it belongs on; see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | 0.40–1.45 | 1.01 / 1.09 |

<!-- test: file=files/rocscience/joints/rj009.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-9, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

A step of refinement — a 2D size of 2.1 m, which takes the mesh from 12 620 nodes and 1 806
interface elements to 26 565 and 2 579 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
17,761 sweeps of the 250,000 allowed.

**The bracket closes on a trial the plain viscoplastic loop cannot settle.** At F = 1.027344 the
loop from a cold start runs away; the corrector, seeded from that same path, reaches equilibrium in
338 sweeps. Both are true of the same discrete model — it has
an admissible static equilibrium at that strength, and the cold-start path does not find it — and
the factor is the strength at which the equilibrium stops existing rather than the strength at
which that path first runs away. Five of the nine trials of each bracket are certified this way.

RS2's own two factors straddle this row — 1.01 without the joint improvement option and 1.09 with
it — which is true of only three other problems in this corpus.

The paper's own limit equilibrium for this problem is a range rather than an answer — 0.40 to 1.45,
as the manual prints it — and it is recorded beside the referee rather than scoring, as every source
limit equilibrium here is. The distinct-element run is what the manual verifies RS2 against.

**Input file:** [rj009.xlsx](files/rocscience/joints/rj009.xlsx).

![RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: the release trace under the crest carries the brightest slip, a bedding plane below the toe carries the rest, and the deformed section shows the slab between them sliding out over the bench](images/RJ-9.png)

### 🟢 RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) {#rj-10}

Example 1a with the release joint moved five metres up the face, which the manual says is the whole
difference between the two problems. The section, the bedding and both friction angles are
[problem 9](#rj-9)'s: a 50 m slope at 50° cut by bedding dipping out of the face at −50° at 3 m
spacing at φ = 30°, release at φ = 40°, and the same elastic rigid-block stand-in for the rock at
E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3.

Moving the release upslope leaves the toe undercut where it was and cuts a second block out of the
face above it: the trace that ran from the face to the toe in example 1a stays, and a new one leaves
the face at (−9.977, 11.890) five metres higher. Both end on the bedding plane they are released by;
see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | 0.43–1.45 | 0.92 / 1.08 |

<!-- test: file=files/rocscience/joints/rj010.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-10, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

A step of refinement — a 2D size of 2.1 m, which takes the mesh from 12,608 nodes and 1,806
interface elements to 26,577 and 2,579 — does not move the factor at all: the finer mesh returns the
same bracket, edge for edge. Every trial of both brackets reaches a verdict. The corpus bracket's
longest trial takes 60,715 sweeps of the 250,000 allowed and the refinement's longest takes 18,471,
and the two brackets get there differently: on the corpus mesh every trial converges or diverges
under its own steam, while on the finer one the five standing trials are all settled by the Newton
corrector at 300 sweeps. The bracket is the same either way, which is what the row is cut on.

RS2's own two factors straddle this row — 0.92 without the joint improvement option and 1.08 with it,
around 1.037 — as they do on [problem 9](#rj-9), [problem 5](#rj-5) and [problem 6](#rj-6) and on no
other problem here. The paper's own limit equilibrium is again a range rather than an answer, 0.43 to
1.45, recorded beside the referee rather than scoring.

**Input file:** [rj010.xlsx](files/rocscience/joints/rj010.xlsx).

![RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) — FEM inputs, mesh, joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and, as in example 1a, almost none of it slips: the release trace high on the face carries the brightest slip, one bedding plane below the toe carries the rest, and the deformed section shows the block between them moving out over the bench](images/RJ-10.png)

### 🟢 RJ-11: Alejano et al. ploughing sliding slab failure (rj011) {#rj-11}

A 25 m slope at 50° with bedding dipping out of the face at −50° at 1.5 m spacing, φ = 30°, and two
release traces at φ = 20°: one at the toe running below the bench and one from the face down onto
the bedding plane that releases the toe block. Ploughing failure is the paper's name for sliding on
a primary discontinuity combining with sliding on a joint sub-parallel to the face, which lifts the
toe block and eventually rotates it out of the slope. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The lower tip of a release trace is stated to six
decimals, so it lands a part in 10⁷ from the bedding plane it belongs on; see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | Rigid-block bound referee | RS2 vs referee | UDEC | Alejano Eq. (7) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.213** | 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | 1.21 (+0.2%) | 1.7582 · the paper prints 1.75 | 1.22 / 1.3 |

<!-- test: file=files/rocscience/joints/rj011.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-11, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 12 326 nodes and 1 768
interface elements to 26 052 and 2 524 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
59,761 sweeps of the 250,000 allowed.

**This is the problem where Alejano's Eq. (7) is not a rigid-block answer.** The two blocks the
mechanism cuts out — the slab and the toe block below it — admit no set of contact forces inside
their friction cones above 1.2148, whatever the distribution; Eq. (7) returns 1.7582, half a factor
above that, so it cannot be what scores this row and the bound is the referee instead. XSLOPE reads
1.213, RS2 1.22 and the paper's own UDEC run 1.21 — three codes inside 0.8% of one another, all
three on the bound. [Problem 15](#rj-15) has the same shape at a problem where no bound is
available.

Every input class matches the vendor model: the rock's
E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile cap, the
bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor clamps in
both directions.

**Input file:** [rj011.xlsx](files/rocscience/joints/rj011.xlsx).

![RJ-11: Alejano et al. ploughing sliding slab failure (rj011) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. One bedding plane from the crest to the toe carries almost all of the slip, and at its foot the two release traces cut out a small wedge: the deformed section shows that wedge lifted and rotated out over the bench while the slab above it slides down the plane, which is the ploughing mechanism the paper names. The wedge is driven out and up by the slab above it, and at the panel's exaggeration a movement of centimetres draws as metres, so the block appears to leave the slope](images/RJ-11.png)

### 🟡 RJ-12: Alejano et al. ploughing toppling slab failure (rj012) {#rj-12}

A 25 m slope at 60° cut by bedding dipping out of the face at −60° at 1.5 m spacing, φ = 30°, with
two short release traces at φ = 40°: one at the toe running below the bench, and one from the face
down onto the bedding plane that releases the toe block. Ploughing failure is the paper's name for
sliding on a primary discontinuity combining with sliding on a joint sub-parallel to the face, which
lifts the toe block and rotates it out; at 60° the rotation rather than the sliding governs, which
is what separates this row from problem 11. The rock is elastic at E = 2 × 10⁸ MPa, γ = 25 kN/m³,
ν = 0.3 — not a rock modulus but the manual's own device for reproducing UDEC's rigid blocks.

Alejano's Eq. (7) is the referee for this problem: a moment balance about the toe for exactly this
mechanism, the slab thrusting on the toe block at one point and the toe block rotating out about
the toe. Recomputed on the inputs the vendor file states, Eq. (7) gives **1.9659** against the 2.00
the paper prints for the same example, and the same implementation reproduces the paper's own table
on its other worked examples.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **2.033** | 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | 2.0324 | 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 |

<!-- test: file=files/rocscience/joints/rj012.xlsx, type=fem_ssrm, expected_fs=2.033, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-12, f_stand=2.0234375, f_fail=2.04296875, check=edges, tier=gate -->

The lower tip of each release trace is stated to six decimals, so it lands a part in 10⁷ from the
bedding plane it belongs on; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint). This row and
problem 14 state their generated bedding network to six decimals as well, which is the precision
the mesher's own crossing arithmetic carries — see the departures table.

The factor is mesh independent over two steps of refinement. The corpus mesh at a 2D size of 1.5 m
gives 2.033; a step to 1.05 m, which takes the mesh from 12,294 nodes and 1,818 interface elements
to 26,407 and 2,572, gives 2.033; and a second step to 0.735 m gives 2.033 again, the same bracket
edge for edge on all three. Every trial of every one of those brackets reaches a verdict, and the
longest takes 16,061 sweeps of the 250,000 allowed, which makes this the least budget-limited row
in the corpus.

**The three independent answers for this problem do not agree with one another, and this row is the
one nearest the closed form.** Alejano's limit equilibrium is 1.9659 recomputed, the paper's UDEC run is 1.78 and
RS2 reads 1.39, on a mechanism whose governing mode — the toe block rotating out rather than
sliding — every source agrees on.

**Input file:** [rj012.xlsx](files/rocscience/joints/rj012.xlsx).

![RJ-12: Alejano et al. ploughing toppling slab failure (rj012) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section, and only a few of its traces carry any slip: one release trace under the crest and the bedding beneath the toe block, which is the ploughing pair. The two right-hand panels are the last standing trial of the bracket, the state below the factor rather than past it — on a rock this stiff the model moves by microns until it does not, so the deformed section is drawn at tens of thousands of times scale. The capture past the factor is not drawn: it was stopped in its first sweep, before the section had moved at all](images/RJ-12.png)

### 🟢 RJ-13: Alejano et al. ploughing sliding slab, example 4 (rj013) {#rj-13}

A 25 m slope at 55° cut by bedding dipping out of the face at −55° at 1.5 m spacing, φ = 25°, with
two release traces at φ = 20°: one at the toe running below the bench and one from the face down
onto the bedding plane that releases the toe block. It is problem 11's mechanism at a steeper
bedding and a weaker one — sliding on the primary discontinuity combining with sliding on a joint
sub-parallel to the face, which lifts the toe block. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The lower tip of each release trace is stated to
six decimals, so it lands a part in 10⁷ from the bedding plane it belongs on; see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

Alejano's Eq. (7) recomputed on this problem's inputs gives **1.0002**, its sliding mode governing,
and the paper prints 1.00 for the same example. It sits on the rigid-block bound for these two
blocks, which is 0.9988, so the closed form and statics agree to three figures here.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **0.998** | 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | 0.9988 | 1.0 (−0.2%) | 1.0 / 1.05 |

<!-- test: file=files/rocscience/joints/rj013.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-13, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 12 141 nodes and 1 788
interface elements to 26,079 and 2,539 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
23,161 sweeps of the 250,000 allowed.

Every input class matches the vendor model: the
rock's E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile
cap, the bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor
clamps in both directions.

**Input file:** [rj013.xlsx](files/rocscience/joints/rj013.xlsx).

![RJ-13: Alejano et al. ploughing sliding slab, example 4 (rj013) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it carries slip: one bedding plane from the crest down to the toe does, and the two release traces at its foot cut out the toe block, which the deformed section shows lifted and pushed out over the bench while the slab above it slides down the plane](images/RJ-13.png)

### 🟢 RJ-14: Alejano et al. ploughing sliding slab, example 5 (rj014) {#rj-14}

Example 4's section at 60° with the two joint strengths the other way round: bedding dipping out of
the face at −60° at 1.5 m spacing at φ = 20° — the weakest bedding of the six — and two release
traces at φ = 30°, one at the toe below the bench and one from the face down onto the bedding plane
it releases. The rock is the family's elastic rigid-block stand-in at E = 2 × 10⁸ MPa,
γ = 25 kN/m³, ν = 0.3. Like problem 12, this row states its generated bedding network to six
decimals; see the departures table.

Alejano's Eq. (7), recomputed on the inputs the vendor file states, gives **1.2034** for this
example, under the rigid-block bound for its two blocks. It is the referee.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.232** | 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | 1.2686 | 0.9 (+36.9%), 0.9994 with the paper's corner rounding corrected · the paper prints 1.00 | 0.89 / 1.09 |

<!-- test: file=files/rocscience/joints/rj014.xlsx, type=fem_ssrm, expected_fs=1.232, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-14, f_stand=1.22265625, f_fail=1.2421875, check=edges, tier=gate -->

**The 1.00 the paper prints for this example does not follow from the inputs it prints.** Eq. (7)
reads the bedding dip and spacing, the release joint's dip, the two friction angles, the length of
the toe block's base and the unit weight. Every one of them matches the vendor's model — the base to
five figures, 3.589 m against the 3.588 m printed — and on them the equation returns 1.2034. For it
to return 1.00 the toe block's base would have to be 6.08 m, which is a different problem. The
manual's other value, UDEC's 0.9, is corner rounding by the paper's own account: with the rounding
radius reduced it reports 0.9994 for the same model.

The factor is mesh independent over two steps of refinement. The corpus mesh at a 2D size of 1.5 m
gives 1.232; a step to 1.05 m, which takes the mesh from 12,316 nodes and 1,819 interface elements
to 26,414 and 2,572, gives 1.232; and a second step to 0.735 m gives 1.232 again, the same bracket
edge for edge on all three. Every trial of every one of those brackets reaches a verdict, the
longest taking 15,511 sweeps of the 250,000 allowed.

Every input class matches the vendor model: the
rock's E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile
cap, the bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor
clamps in both directions.

**Input file:** [rj014.xlsx](files/rocscience/joints/rj014.xlsx).

![RJ-14: Alejano et al. ploughing sliding slab, example 5 (rj014) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: one bedding plane from the crest to the toe carries the slip, with the release trace at the toe opening as the slab above it slides out over the bench](images/RJ-14.png)

### 🔴 RJ-15: Partially joint-controlled footwall slope (rj015) {#rj-15}

A 40 m footwall at 40° whose bedding dips in the same direction at the same angle, 2 m apart, so
the slabs lie parallel to the face and a failure has to break rock at the toe to get out. This is
the one problem in the Alejano family whose rock can yield: Mohr-Coulomb, c = 200 kPa, φ = 35°,
γ = 28 kN/m³, E = 1 GPa, ν = 0.3. The joints are the corpus's softest — k<sub>n</sub> = 5 × 10<sup>6</sup>
kPa/m and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m, twenty times below the set's standard pair — with
no cohesion and φ = 25°.

The referee is Alejano's footwall limit equilibrium, Eqs. (9)–(10) of the source paper, recomputed
on the inputs the vendor file states. It resolves one block along and normal to its own break-out
plane, and at the optimum the paper states — a break-out inclined 14° to the bedding and emerging at
55°, with the failure taken at the face — it gives **1.7985** against the value the paper prints.
Minimized over its own two angles it settles within half a degree of that same optimum.

| XSLOPE SSRM | Alejano Eqs. (9)–(10) referee | RS2 vs referee | UDEC-SSRT (Alejano) | Slide2 LEM (vendor) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.271** | 1.7985 (−29.3%) | 1.28 vs 1.7985 (−28.8%) | 1.6 (−20.6%) | 1.25 | 1.28 / 1.42 |

<!-- test: file=files/rocscience/joints/rj015.xlsx, type=fem_ssrm, expected_fs=1.271, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-15, f_stand=1.26171875, f_fail=1.28125, check=edges, tier=gate -->

**The closed form prices one mechanism, and every program free to search for a surface finds a
weaker one.** Eqs. (9)–(10) admit a single slab of the stated bedding thickness driving a wedge out
through the rock, and their own lengths hard-wire the break-out to cross exactly one 2 m bed; the
factor they return rises with that thickness, so the shallowest case the formula allows is also the
lowest it can report. Three programs of three different kinds are under no such restriction: this
SSRM at 1.271, RS2's SSR at 1.28 and Rocscience's own Slide2 limit-equilibrium search at 1.25, the
three inside 2.4% of one another and all three far below the closed form, with the paper's
distinct-element run at 1.6 between. **Both finite element codes miss the referee by the same
amount and in the same direction**, which makes the gap a property of what the closed form is
allowed to consider rather than of either program.

The upper edge of the corpus bracket is read as failed on a steadily slipping interface at 225,001
sweeps, and every trial of the bracket decides. A step of refinement — a 2D
size of 1.4 m, which takes the mesh from 14,964 nodes and 2,179 interface elements to 32,131 and
3,090 — moves the factor by one bracket step, inside the row's own tolerance, and decides on all
nine of its trials too, so the row is locked on the corpus mesh. It is the corpus's longest row: its
bracket and at-failure capture together take about an hour and a half.

Every input class matches the vendor model: the rock's E, ν, γ, c, φ and
tensile cap; the joints' normal and shear stiffness, cohesion, friction angle and tensile cap; and
the side restraint the vendor clamps in both directions. The stated tensile strength of 1000 kPa is
above the Mohr-Coulomb apex its own c and φ imply (c/tan φ = 285.6 kPa), so it never binds.

**Input file:** [rj015.xlsx](files/rocscience/joints/rj015.xlsx).

![RJ-15: partially joint-controlled footwall slope (rj015) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The bedding slips over a long stretch behind the face, and the rock's only strain is a small patch at the toe where the slab has to break through to get out — the coupled mechanism the source paper describes](images/RJ-15.png)

### <span class="nodata">⊘</span> RJ-16: Barla et al. tilt-table block toppling (rj016) {#rj-16}

A laboratory test rather than a slope. Fourteen columns of 9 cm blocks are stacked into a 63.4°
staircase on a plate, and the plate is tilted until the stack topples; what the problem reports is
the angle at which it goes. The blocks are elastic (E = 350 MPa, ν = 0.2, γ = 28 kN/m³) on a plate
that is elastic and effectively rigid (E = 200 GPa), and the 23 joints — thirteen vertical, eight
horizontal, the plate contact along the base of the stack and the backstop behind it — carry no
cohesion, φ = 38°, and the softest stiffness pair in the corpus: k<sub>n</sub> = 5 × 10<sup>6</sup>
and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m. The 0° model is also restrained with rollers, where
the manual's slope models are clamped on their sides.

The vendor tilts the model itself: ten files, `#016_0deg` through `#016_9deg`, each the whole
section rotated one degree further with its boundary pinned. XSLOPE turns the load instead. The
finite element engine carries a horizontal seismic coefficient k, applied as a body force kγ whose
sign is its direction, so a tilt of θ is k = tan θ toward the face, and the row is measured by
sweeping k at full strength — F = 1, no reduction — until the stack stops standing. Every point is
one solve on the corpus mesh, 0.09 m, which is the block size.

| XSLOPE tilt | UDEC referee | RS2 vs referee | Experiment | RS2 without / with improvement |
|---|---|---|---|---|
| **8.5°** | 11° | 9° vs 11° (−18.2%) | 9° | 9° / 7° |

The stack stands at k = 0.1477, which is 8.40°, and goes at k = 0.1504, which is 8.55°. It stands
at every coefficient the sweep tried below that and fails at every one above, through 0.35, and
each verdict is the solver's own: the standing trials are certified by the Newton corrector at 300
sweeps, and the failing ones diverge inside a hundred.

**A coefficient is not a rotation, and here the difference is measurable and measures zero.**
Tilting the model by θ turns the body force through θ and leaves its magnitude at γ; a coefficient
k = tan θ turns it through the same angle and multiplies its magnitude by 1/cos θ, which is 1.1% at
8.5°. On this model that cannot move the answer, because every joint carries zero cohesion and no
block can yield, so what the state depends on is the direction of the body force and not its size.
Re-solving the two bracketing coefficients with every unit weight scaled by 1.011 and then by 0.5
returns the same two verdicts at the same sweep counts.

**The row is reported rather than locked**, because what it measures is not a factor of safety: the
sweep asks whether the model stands at full strength under a stated load, which is a different
question from the strength reduction every other row here answers, and the corpus's locking rule is
written for the second. Its figure is drawn the same way — two solves at the two coefficients that
bracket the tilt, titled by the angle rather than by a factor.

The plate is weightless in the vendor model (`BodyForceSolid: 0` — it is the apparatus, not rock),
and `build_fem_data` requires a positive unit weight, so it is built at the 27 kN/m³ its own
property row states. Its weight is carried by its own restraints and the contact stress along the
base of the stack is the weight of the blocks above it.

**Input file:** [rj016.xlsx](files/rocscience/joints/rj016.xlsx).

![RJ-16: Barla et al. tilt-table block toppling (rj016) — FEM inputs with the seismic coefficient that stands for the tilt, mesh with the vendor's rollers, joint slip at the first coefficient the stack goes at, and the deformed section. The slip is on the vertical joints between the columns and on the bedding under the crest of the stack, and the deformed section at 156x shows what that adds up to: every column leaning downslope about its own base, the tall ones at the back furthest over, which is toppling rather than the stack sliding along the plate](images/RJ-16.png)

### 🟡 RJ-17: Step-path failure, en-echelon joints (rj017) {#rj-17}

[Problem 18](#rj-18)'s section — a 45 × 20 m block of one Mohr-Coulomb rock (γ = 19.62 kN/m³,
E = 20 GPa, ν = 0.3, c = 25 kPa, φ = 25°, no tensile capacity) whose face rises from (17, 8.2) to
(26.9, 20) — cut by three joints at 36.1° that stop short of one another instead of running
through: (18.33, 9.37) to (22.33, 12.25), (22.33, 13.5) to (26.26, 16), and (26.26, 17.16) to
(29.33, 19.33). The rock bridges between them are what a step-path failure has to break through,
and they are why this slope stands where problem 18's continuous joints let it go. The joints carry
problem 18's own strength: c = 1 kPa, φ = 35°, and the corpus's standard stiffness pair.

**The strength reduction is confined to a rectangle.** The vendor model carries an SSR search area —
a polygon from (12.221, 3.19011) to (40.2504, 20.532), which the manual's geometry figure draws as
a dashed box — and RS2 applies it by holding every element whose centroid falls outside it linear
elastic. The model file shows that done: an auto-generated elastic twin of the rock, the same
modulus and no plasticity of its own, on 1,173 of its 2,579 elements. XSLOPE states the same
constraint as a polygon overlay on the model, classified by the same element-centroid test, so what
is transcribed is the rectangle the vendor states.

No closed form exists for a step path through rock bridges, so the referee is the manual's single
distinct-element run.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.213** | 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | 1.24 / 1.2 |

<!-- test: file=files/rocscience/joints/rj017.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-17, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.7 m, which takes the mesh from 3,284 nodes and 14 interface
elements to 6,484 and 21 — moves the factor by one bracket step, inside the row's own tolerance, to
1.193. Every trial of both brackets reaches a verdict. The row is cut on the corpus mesh, as every
row here is.

Four of the nine trials of each bracket are settled by something other than the force test: two are
certified by the Newton corrector 300 sweeps in, the standing edge among them, and the two just
above that edge run to 225,001 sweeps and are read as failed on a steadily slipping interface.

**What decides this row is the rock bridges.** The rock has no tensile capacity at all, so the
intact ligaments between the joint segments carry nothing across them, and the strength reduction
takes their 25 kPa cohesion down alongside the joint friction: the factor is the strength at which
three short bridges shear through. Both finite element codes read that below the distinct-element
run and on the same side of it, 2.1 points apart.

**The same zone, transcribed the way the vendor's mesh resolved it, is the same model.** The
element map of the vendor file is that rectangle rasterized onto element edges: a staircase
wandering about half an element either side of it, from x = 11.85 to 12.73 where the rectangle says
12.221 and from y = 2.73 to 3.66 where it says 3.19011. Transcribed as an 84-vertex ring and run as
a file of its own, it holds 693 elements of this mesh elastic where the rectangle holds 694 — the
zone is an analysis overlay and never meshed, so both files mesh identically and only the
membership can differ — and the bracket comes back identical trial for trial, sweep count for sweep
count. What a mesh does to the edge of a search area is worth one element here, and nothing at all
to the factor.

**Input files:** [rj017.xlsx](files/rocscience/joints/rj017.xlsx), and the staircase variant
[rj017_staircase.xlsx](files/rocscience/joints/rj017_staircase.xlsx).

![RJ-17: step-path failure with en-echelon joints (rj017) — FEM inputs with the vendor's search area drawn as the held-elastic region around it, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. All three joints slip along their whole length, and the strain runs between them: a band climbs from the toe of the face through the rock bridge below the lowest joint and on past the upper two, which is the step-path the problem is named for — shear along the joints and broken rock between them](images/RJ-17.png)

### 🟢 RJ-18: Step-path failure, continuous joints (rj018) {#rj-18}

A 45 × 20 m section of one Mohr-Coulomb rock (γ = 19.62 kN/m³, E = 20 GPa, ν = 0.3, c = 25 kPa,
φ = 25°, no tensile capacity) with a slope face rising from (17, 8.2) to (26.9, 20), cut by three
parallel joints at 36.1° that run from the face to the crest at a perpendicular spacing of
0.883 m. The joints carry c = 1 kPa, φ = 35°, k<sub>n</sub> = 10<sup>8</sup> kPa/m,
k<sub>s</sub> = 10<sup>7</sup> kPa/m and are reduced with the rock in the strength reduction.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **0.998** | 1.01 (−1.2%) | 1.01 vs 1.01 (0.0%) | 1.01 / 1.00 |

<!-- test: file=files/rocscience/joints/rj018.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-18, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.7 m, which takes the mesh from 3,486 nodes and 48 interface
elements to 6,707 and 66 — does not move the factor at all, and every trial of both brackets
reaches a verdict. The budget is what the row needs rather than what it has to spare: the longest
trial to decide takes 196,201 sweeps of the 250,000 allowed.

**Input file:** [rj018.xlsx](files/rocscience/joints/rj018.xlsx).

![RJ-18: step-path failure through three continuous joints (rj018) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Every joint is slipping along its lower half and open along its upper, and the three slabs between them slide out together down the 36.1° path; the rock itself carries almost no plastic strain](images/RJ-18.png)

### 🔴 RJ-19: Bi-planar step-path failure (rj019) {#rj-19}

A 120 × 70 m section of one Mohr-Coulomb rock (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 10,500 kPa,
φ = 35°, tensile capacity 200 kPa) with a slope face from (30, 20) to (60, 70), cut by two
discontinuous joints with a rock bridge between them: a basal joint at 28.4° from (39.0149,
35.0248) to (63, 48), and an upper joint at 56.3° from (62, 49) to (76, 70). Both carry c = 0,
φ = 40° and the same stiffness pair as RJ-18.

**What decides this row is a tensile cap rather than a cohesion.** The two joints leave a rock
bridge 1.414 m long between their tips. At the rock's 10,500 kPa cohesion that bridge carries
14,849 kN/m in shear, which is 1.2 times the whole 11,929 kN/m sliding mass and 2.6 times its
component down the basal joint, so it cannot be sheared through; what lets the block move is the
same rock's 200 kPa tensile capacity, worth 283 kN/m across the bridge. Every model in this manual
reduces that capacity with the trial factor, and this row is run the same way — see the departures
table. It is the one problem in the corpus where that setting reaches the answer.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.623** | 1.46 (+11.2%) | 1.5 vs 1.46 (+2.7%) | 1.5 / 1.41 |

<!-- test: file=files/rocscience/joints/rj019.xlsx, type=fem_ssrm, expected_fs=1.623, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=true, k0=1, benchmark=RJ-19, f_stand=1.61328125, f_fail=1.6328125, check=edges, tier=gate -->

Every trial of the corpus bracket reaches a verdict. Four of the nine are settled by the Newton
corrector reaching equilibrium from the loop's own state 303 to 385 sweeps in, the standing edge
among them, and the failing edge runs away at 194,241 sweeps of the 250,000 allowed.

A step of refinement — a 2D size of 2.1 m, which takes the mesh from 3,485 nodes to 6,910 — moves
the factor by one bracket step, inside the row's own tolerance, but it does not settle it: the
refined bracket's standing edge is answered only by a corrector certification at the budget exit and
its failing edge reaches 250,000 sweeps with nothing to say. The row is cut on the corpus mesh,
where every trial of the bracket decides inside its budget.

The manual's table for this problem states one joint inclination as 59°. Its figure dimensions 56°
and 28°, and the vendor model's own endpoints give 56.3° and 28.4°, so the table is the outlier and
the model is what is built.

**Input file:** [rj019.xlsx](files/rocscience/joints/rj019.xlsx).

![RJ-19: bi-planar step-path failure with a rock bridge (rj019) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Both joints have opened, the basal one is slipping at its lower end, and the only strain in the rock is the patch at the bridge between the two joint tips, where the block above has to break through to move](images/RJ-19.png)

---

### ⊘ RJ-20: Hammah & Yacoub Voronoi slope (rj020) {#rj-20}

An 80 × 70 m section with a 60 m face at 71.6°, from the toe at (10, 10) to the crest at (30, 70),
tessellated into Voronoi blocks over the whole of it. The rock is Mohr-Coulomb (γ = 27 kN/m³,
E = 20 GPa, ν = 0.3, c = 1,000 kPa, φ = 35°, no tensile capacity) and every block wall is a joint at
c = 500 kPa, φ = 20° with the corpus's standard stiffness pair. The paper this comes from asks how
the failure of a slope in blocky rock changes with the scale of its blocks; the manual's answer at
this scale is UDEC's 2.46, with RS2 reporting 2.21 without joint improvement and 2.37 with it.

**The network is in the model file.** The tessellation was generated in UDEC and imported, so RS2
holds it as 523 joint boundaries — polylines of two to seven points — rather than as anything a
block size and a seed regenerate. They are transcribed verbatim: 1,177 segments, one row of the
joints sheet each, and every block corner a three-way junction of them. Measured on the traces
themselves: 525 blocks over the section, a mean block area of 8.381 m² and a mean block width of
2.895 m. That width is this row's mesh size — the block scale, standing in for the joint spacing
every other row meshes at — and it gives 10,864 nodes, 3,567 elements and 4,582 interface elements
against the vendor's own 3,251 elements.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| *no lock* | 2.46 | 2.21 vs 2.46 (−10.2%) | 2.21 / 2.37 |

The row prints no factor, and what holds it back is the sweep budget. Five of the nine trials
decide: three stand, each settled by the joint verdict on a slip field that has stopped growing,
and two fail by running away. The other four reach 250,000 sweeps with nothing to say — and both of
the trials the bisection closes on are among them, one exiting `STABLE_STUCK` and the other
inconclusive. The corrector certified none of them. A factor closed on two such trials would be a
statement about the budget rather than about the slope. The bracket and the at-failure capture together take 3.7 hours, the longest single run here.

The figure shows how a mass of blocks at this scale fails: the slip picks its way from the toe up through
the block walls on a curved path to the crest, and the only rock strained is a patch at the toe
where the path has to turn. A mass of blocks at this scale fails on a surface, not along any
plane it contains — which is the observation the source paper is about.

**How much of that is the block size and how much is this particular tessellation** is the paper's
own question, and answering it takes a second network of the same blocks drawn another way. The one
`xslope.joints.voronoi` draws at that block size stops in the mesher: the generator keeps a trace
down to a thousandth of the section diagonal, 0.106 m here, and a trace that short is pulled onto
its own junction, where the line it leaves has no length left. Three seeds and two block sizes stop
the same way, so this row is scored on the vendor's own tessellation alone.

**Input file:** [rj020.xlsx](files/rocscience/joints/rj020.xlsx).

![RJ-20: Hammah & Yacoub Voronoi slope (rj020) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The slip runs from the toe up through the block walls on a curved path to the crest, taking whichever wall of each block lies nearest that line, and the rock carries almost no strain except a patch at the toe where the path turns: the mass fails on a surface picked out of the tessellation rather than along any one joint in it](images/RJ-20.png)

## The Dilation Problem {#the-dilation-problem}

Problem 23 is a direct shear test on one joint rather than a slope: two blocks are pressed together,
first at 3 MPa and then at 9 MPa, and one of them is dragged sideways. What it reports is shear
stress against slip, a curve, so there is no factor of safety in it to lock. The manual uses the
test to show the joint law working — the peak strength, the drop to a residual strength once the
joint has slipped, and dilation.

The six vendor models are named for dilation angles of 0, 10, 20, 20, 20 and 30 degrees, and every
one of them carries `include_dilation: no`. The angle in a file's name never reaches the solver, so
all six ran at zero dilation, and two of the three 20-degree cases also step the normal load from 3
to 9 MPa at different points in the test, which makes them different tests. The manual's dilation
comparison therefore never exercised dilation, and XSLOPE's dilation cannot be scored against it.
It is checked against the kinematics instead: on a sliding joint the opening per unit slip is the
tangent of the dilation angle, which `test/joint_element_check.py` measures at row 6, and the peak
and residual strengths are checked there as well, each against its closed form.

## Where the Vendor Models Depart from the Manual

Each of these was found by reading the `.fez` against the manual page it belongs to, and each
changes what a faithful transcription is.

- **Problem 7's slope angle.** The table prints 5°; the figure and the model are the same 55°
  slope its siblings use.
- **Problem 19's joint inclination.** The table prints 59°; the figure dimensions 56° and 28°, and
  the model's own endpoints give 56.3° and 28.4°.
- **Problem 15's slope height.** The table prints 25 m. The vendor's section rises 40 m from the toe
  at (0, 0) to the crest at (−47.6701, 40), which is the height the source paper states for the same
  model, and the section is what is built.
- **Problem 23's dilation.** All six "dilation" models run at zero dilation — see
  [The dilation problem](#the-dilation-problem).
- **Problem 1's toe force.** The manual describes a stabilizing force at the toe of the lowest
  block. It is 2013 kN in cases b and d — about twice that block's own weight — and 0.5 kN in cases
  a and c, which is 0.05% of it and does nothing. The force is horizontal and points into the
  slope in every one of the six models. The vendor's own rerun models move it from the toe at
  (−0.5, 0.866) to the upper-left block corner at (−2.5, 4.330) in all four cases, so the "with
  joint improvement" factors are not the same load case as the "without" ones. Cases b and d here
  carry the force at the **corner**: the toe is the end of the basal joint, where the split gives
  the node one copy per wedge of material around it and a load applied there has no defined side to
  act on. Cases a and c carry no load, the 0.5 kN being what it is.
- **Problem 20's network.** The manual calls it a Voronoi tessellation generated in UDEC and
  imported, and publishes no block size, density or seed for it. The model file carries it as 523
  joint boundaries — 1,177 segments — so those traces are the only statement of the network there
  is, and they are what [RJ-20](#rj-20) transcribes.
- **Problem 16's boundary conditions.** The 0° case runs on rollers; the nine tilted cases pin
  every exterior node in both directions, and use a convergence tolerance two orders tighter. The
  tilt plate carries no body force in any of the ten: it is apparatus, not rock.
- **Problem 17 confines its strength reduction.** The model carries an SSR search area — a
  rectangle from (12.221, 3.19011) to (40.2504, 20.532), stated in the file and drawn as a dashed
  box on the manual's geometry figure — and RS2 applies it by holding every element outside it
  linear elastic, which the file shows as an auto-generated elastic twin of the rock on 1,173 of
  its 2,579 elements. The manual's tables say nothing about it, and a transcription that left it
  out would be reducing the strength of the whole section.
- **Problems 3 and 5 are elastic.** Both files set the plasticity specification to none and carry
  no strength values under it, so the rock cannot yield and only the joints can.
- **Problems 12 and 14 state their generated network to six decimals.** A trace generated from a
  stated dip and spacing is clipped to the section, so its ends are computed points on the
  boundary, and the mesher's own crossing arithmetic carries six decimals: at ten, the outline
  gains the rounded crossing while the joint keeps its own end 2 × 10⁻⁷ away, and gmsh does not
  finish recovering the 1D mesh at these two rows' joint spacing. Six decimals is 10⁻⁶ m on a 90 m
  section. The other four rows of the family mesh as the generator states them and are left that
  way, so that the meshes their own factors were measured on do not move.
- **Problems 9 to 14 use two joint strengths, not one.** The network's 1928–1962 elements carry one
  friction angle and the one or two short explicit crest traces carry another — 40°, 20° or 30°
  against the network's 30°, 25° or 20° — so a single quoted joint friction angle for these
  problems is incomplete.
- **Problem 22's units.** The joint numbers (c = 143, residual c = 76, load 345) read as kPa in a
  file declared in MPa: at MPa a 143 MPa joint under 345 MPa never slips at the prescribed 0.3 m of
  shear, and at kPa it is an ordinary rock-joint shear test.
- **The manual states no rock strength for problems 3 to 7.** Their tables carry the slope
  geometry, the joint friction angle and the rock's tensile strength, and nothing else, where the
  models of problems 4, 6 and 7 carry c = 0.675 MPa and φ = 43°. The same tables state no joint
  cohesion either, where all five models carry 0.1 MPa — on a 260 m slope, not a rounding.
- **Problem 8's tensile strength never binds.** The model caps the rock at 75 kPa, and the
  Mohr-Coulomb apex its own c and φ imply is c/tan φ = 74.1 kPa, so the cap sits above the strength
  the envelope already has. Problem 15 is the same: a 1000 kPa cap over a 285.6 kPa apex.
- **Problem 8's joint network is not the whole section.** Its twelve columns are clipped to the
  block the basal joint and the back joint bound, which is how the vendor file stores them: the
  three highest end **on** the back joint at (68.4, 6.032), (68.4, 16.192) and (68.4, 26.352). The
  same set generated over the whole section would be sixteen chains, three of them below a plane the
  model has no columns under and three more running on past x = 68.4 into the strip the vendor
  leaves uncut.
- **Problems 9 to 14 run on a rock a thousand times stiffer than steel.** E = 2 × 10⁸ MPa is the
  manual's own device, and it says so: the UDEC models these are scored against use rigid blocks,
  and the modulus is how RS2 reproduces one. Problem 15, whose rock can yield, runs at 1 GPa.
- **No problem from 1 to 21 states a joint residual strength or a dilation angle.** Those appear
  only in problems 22 and 23, and material residual values, where they appear, always equal the
  peak.
- **Every model relieves a slipping joint's stiffness.** All 23 vendor files carry
  `joint_stiffness_flag: 1` with `joint_stiffness_factor: 0.01`: the stiffness of a joint past its
  strength criterion is dropped a hundredfold in the assembled matrix. XSLOPE has the same relief
  (`joint_tangent`, at the vendor's own factor) and the corpus runs without it. Turned on, it takes
  [problem 1a](#rj-1a)'s bracket down two bisection steps, onto the factor RS2 and UDEC both report
  for that case, so it is the departure that places the two of them below the closed form on the
  toppling rows. It is off here because of what it does to the
  iteration on a column stack, where almost every closed pair is slipping: taking all but a
  hundredth of the interface shear stiffness out of the matrix leaves the assembly with very little
  lateral stiffness in it, and the loop runs away rather than settling sooner.
- **Every model reduces the rock's tensile strength with the trial factor.** All 23 vendor files
  carry `tensilestrength_SRF: 1`, so RS2 divides the tensile cap by the trial factor as it divides
  cohesion and friction. The corpus runs that setting on [problem 19](#rj-19), the one problem where
  the cap governs the answer. Everywhere else it is inert and the rows are cut without it: the
  elastic rows have no cap to reduce, problems 4, 6, 7, 17, 18 and 20 carry a cap of zero, and
  problems 8 and 15 carry one above the Mohr-Coulomb apex their own c and φ imply — an apex that
  does not move under a strength reduction, because c and tan φ are divided by the same factor.
