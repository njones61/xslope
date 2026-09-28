# Rocscience RS2 Joint Corpus

The [RS2 Joint Verification Manual](https://www.rocscience.com/help/rs2/verification-theory/verification-manuals)
(Rocscience) publishes 23 problems on jointed rock: block and flexural toppling, plane failure,
Alejano's sliding and plowing slabs (spelled "ploughing" in the papers and the manual), step-path failure, a Voronoi-tessellated mass, a jointed
tunnel, and two shear-box problems that exercise the joint's own constitutive law rather than a
slope. Every one of them is a mesh split along discontinuities with interface elements carrying the
normal and shear stress between the faces, which is what XSLOPE's [joint lines](../fem/joints.md)
are; how they are
modeled, and the reach of the element that carries them, is documented there. The rows below verify
XSLOPE's FEM/**SSRM** solver on that corpus.

The wall and embankment rows that use the same element but reach it through the reinforce sheet's
`Joint` column are on the [RS2 corpus page](rs2.md) — [RS2-24](rs2.md#rs2-24) and
[RS2-48–55](rs2.md#rs2-48). Full bibliographic details for the author-year citations here are on the
shared [References](references.md) page.

## Methodology

- **Models.** Geometry, materials, joint properties, restraints and loads come from the vendor's
  own `.fez` files rather than the manual's tables, which carry errata the models do not; each
  row states where its model departs from the manual. The vendor models are in MPa and MN/m³;
  these files use kPa and kN/m³.
- **Transcription.** No problem from 1 to 21 states a joint residual strength or a dilation angle,
  and material residual values, where they appear, equal the peak. Every vendor file drops a
  slipping joint's stiffness a hundredfold (`joint_stiffness_factor: 0.01`); XSLOPE has the same
  relief at the same factor, and the corpus runs without it because on the reinforced walls and
  embankments of the [RS2 corpus](rs2.md) it moves brackets by one to four steps of the search.
  Every vendor file also divides the rock's tensile cap by the trial factor, which the corpus does
  only on [problem 19](#rj-19), the one problem where the cap governs the answer.
- **Referee.** Where a closed-form rigid-block limit equilibrium exists for a problem it is the
  referee, recomputed from the inputs the model carries: Goodman & Bray's column analysis on
  problems 1 and 2, Alejano's plowing equation on problems 11 to 14, the sliding block with its tensile bridge on
  problem 19, and, on problem 15, the vendor's Slide2 limit-equilibrium search, since Alejano's
  footwall equations price a mechanism the slope does not take. Where none exists the referee is the one the manual names, UDEC in every such case.
  Each row shows the recomputed value beside the one its source prints.
- **Rigid-block bound.** Each plowing problem is two rigid blocks. The highest reduction factor at
  which some set of joint forces, each inside its friction cone, can still hold both blocks in
  place is a ceiling on any rigid-block answer: above it the blocks must move. Alejano's Eq. (7)
  assumes one particular mechanism, so it may sit below that ceiling but not above it. On problems
  12, 13 and 14 it sits on or below the ceiling and is the referee; on problem 11 it gives 1.76
  against a ceiling of 1.21, so the ceiling is the referee there.
- **RS2's two factors.** The manual reports each problem with and without the vendor's
  `Improve Joint Convergence` option. The run without it is the vendor's default and the same
  method as XSLOPE's, a strength reduction on a continuum with interfaces, so it is the yardstick.
  Both are recorded in every row; their spread is the width of the vendor's own answer.
- **Mesh.** Each row is meshed at its joint spacing, one element across the rock between one
  discontinuity and the next, and states what a step of refinement does to its factor. Problem 1's
  four cases are meshed at 5 m, half their 10 m block width, where three finer meshes return the
  same bracket.
- **Search.** The factor of safety is found by bisection on the bracket [a, b] between the highest
  factor at which the slope stands and the lowest at which it fails. A jointed model takes tens of
  thousands of iterations to settle, so every row allows 250,000. A search that closes on a trial
  still undecided at that limit reports the bracket without a value. See [running a jointed
  model](../fem/joints.md#running-a-jointed-model).

<!--
For maintainers (kept out of the page text):
- The decision rule for a trial: it stands when the viscoplastic loop converges, when the joint
  verdict reads slip and field both stopped, or when the Newton corrector certifies the loop's
  state; it fails on runaway displacement or on a steadily slipping interface with the field
  static. Every row's trial record is written by the figure producer or tools/lock_edges.py; the
  audit is tools/ssrm_trial_audit.py and test/corrector_certified_check.py holds the records to
  that rule. A corrector refusal is not a verdict: a refused trial is decided only by the sweep's
  own readings inside its budget, never by the refusal itself (r19, r22, r30 in the private
  reports). The page text says none of this; it says only whether a trial settled inside the
  budget.
- benchmarks/rocscience/build_joint_problems.py writes the input files and
  benchmarks/rocscience/make_rs2_joint_figures.py the figures, sidecars and records; builders are
  authoritative, corpus workbooks are never patched by hand.
- The referee values are recomputed by the scripts under the private reports r25 and r27
  (Goodman & Bray, Alejano Eqs. 7 and 9–10, the rigid-block bound), validated on the sources' own
  worked examples before use.
-->

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
| [2](#rj-2) | 🟢 | Alejano & Alonso block toppling | SSRM 0.783 vs Goodman & Bray 0.7734 (+1.2%) | 0.86 vs 0.7734 (+11.2%) | UDEC 0.87 (−10.0%) | 0.86 / 0.82 | All three programs stand above the closed form, this one by a percent and the other two by a tenth; the manual states no UDEC settings for this model. |
| [3](#rj-3) | <span class="nodata">⊘</span> | Lorig & Varona forward block toppling | SSRM 1.213 reported vs UDEC 1.13 (+7.3%) | 1.12 vs 1.13 (−0.9%) | — | 1.12 / 1.09 | *reported, no lock* — the trial at the top of the bracket, F = 1.222656, does not settle within the 250,000-iteration limit, so the search cannot close. The vendor's two settings disagree with each other. |
| [4](#rj-4) | 🟢 | Lorig & Varona flexural toppling | SSRM 1.311 vs UDEC 1.3 (+0.8%) | 1.19 vs 1.3 (−8.5%) | — | 1.19 / 1.27 | |
| [5](#rj-5) | <span class="nodata">⊘</span> | Lorig & Varona backward block toppling | SSRM 1.818 reported vs UDEC 1.7 (+6.9%) | 1.65 vs 1.7 (−2.9%) | — | 1.65 / 1.86 | *reported, no lock* — four of the nine trials do not settle within the 250,000-iteration limit, two of them the ends of the bracket; a finer mesh settles two more and leaves an end of its own bracket unsettled. The reported value sits between the vendor's two settings. |
| [6](#rj-6) | <span class="nodata">⊘</span> | Plane failure, daylighting | SSRM 1.271 reported vs UDEC 1.27 (+0.1%) | 1.25 vs 1.27 (−1.6%) | — | 1.25 / 1.31 | *reported, no lock* — every trial settles, but a step of mesh refinement moves the factor by twice the row's tolerance. The reported value sits between the vendor's two settings. |
| [7](#rj-7) | 🟡 | Plane failure, non-daylighting | SSRM 1.564 vs UDEC 1.5 (+4.3%) | 1.57 vs 1.5 (+4.7%) | — | 1.57 / 1.59 | Both finite element codes land above the referee, on the same side and within half a point of each other. |
| [8](#rj-8) | 🟢 | Flexural toppling, base friction model | SSRM 0.764 vs UDEC 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | — | 0.75 / 0.75 | |
| [9](#rj-9) | 🟢 | Bilinear slab failure, example 1a | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | LE (Alejano) 0.40–1.45 | 1.01 / 1.09 | |
| [10](#rj-10) | 🟢 | Bilinear slab failure, example 1b | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | LE (Alejano) 0.43–1.45 | 0.92 / 1.08 | |
| [11](#rj-11) | 🟢 | Plowing sliding slab failure | SSRM 1.213 vs rigid-block bound 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | UDEC 1.21 (+0.2%) · Alejano Eq. (7) 1.7582 | 1.22 / 1.3 | Alejano's Eq. (7) returns a factor above the bound on this problem, so the bound is what scores it; all three programs sit on the bound. |
| [12](#rj-12) | 🟡 | Plowing toppling slab failure | SSRM 2.033 vs Alejano Eq. (7) 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | UDEC 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 | |
| [13](#rj-13) | 🟢 | Plowing sliding slab, example 4 | SSRM 0.998 vs Alejano Eq. (7) 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | UDEC 1.0 (−0.2%) · Alejano prints 1.0 | 1.0 / 1.05 | |
| [14](#rj-14) | 🟢 | Plowing sliding slab, example 5 | SSRM 1.232 vs Alejano Eq. (7) 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | UDEC 0.9 (+36.9%) · Alejano prints 1.00 | 0.89 / 1.09 | The 1.00 the paper prints for this example does not follow from the inputs it prints; Eq. (7) on them gives 1.2034. |
| [15](#rj-15) | 🟢 | Partially joint-controlled footwall | SSRM 1.271 vs Slide2 LE search 1.25 (+1.7%) | 1.28 vs 1.25 (+2.4%) | Alejano Eqs. (9)–(10) 1.7985 (single-bed formula; the paper prints 1.72) · UDEC 1.6 (−20.6%) | 1.28 / 1.42 | Alejano's closed form drives a wedge out through a single 2 m bed and its factor rises with bed thickness; the three programs free to search for a surface agree at 1.25–1.28, so the limit-equilibrium search is the referee, as the rigid-block bound is on problem 11. |
| [16](#rj-16) | 🔴 | Barla et al. tilt-table block toppling | Tilt 10.24° vs UDEC 11° (−6.9%) | 9° vs 11° (−18.2%) | Experiment 9° · Goodman & Bray, plate tilted, 7.6° | 9° / 7° | A tilt angle rather than a factor of safety: pushed with a seismic coefficient at full strength, the stack stands at 10.20° and topples at 10.28°, above the rigid-column bound and the physical stack and below the distinct-element code. |
| [17](#rj-17) | 🟡 | Step-path, en-echelon joints | SSRM 1.213 vs UDEC 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | — | 1.24 / 1.2 | No closed form: three rock bridges decide the factor, and both finite element codes read them below the distinct-element run, 2.1 points apart. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | 1.01 vs 1.01 (0.0%) | — | 1.01 / 1.0 | |
| [19](#rj-19) | 🟢 | Bi-planar step-path failure | SSRM 1.623 vs rigid-block limit equilibrium 1.5914 (+2.0%) | 1.5 vs 1.5914 (−5.7%) | UDEC 1.46 (+11.2%) | 1.5 / 1.41 | One block sliding on the basal joint, held by friction and by the rock bridge in tension; its statics, with the tensile cap reduced with the trial factor as the vendor reduces it, is the referee, and both vendor numbers sit below it. |
| [20](#rj-20) | <span class="nodata">⊘</span> | Hammah & Yacoub Voronoi slope | UDEC 2.46 | 2.21 vs 2.46 (−10.2%) | — | 2.21 / 2.37 | *reported, no lock* — four of the nine trials do not settle within the 250,000-iteration limit, including both of the trials the search closes on. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | 8.27 vs 8.16 (+1.3%) | — | 8.27 / 8.5 | *not supported* — the second stage of the vendor model excavates a 2 m opening and the strength reduction runs on the excavated state, which carries the stress the first stage left behind; staged excavation is outside a slope program's scope. |
| 22 | <span class="nodata">⊘</span> | Joint model: hyperbolic softening | — | — | — | — | *not supported* — the problem exercises RS2's hyperbolic displacement- and work-softening joint law, which XSLOPE's interface element does not have; it reports no factor of safety. |
| 23 | <span class="nodata">⊘</span> | Joint model: residual strength and dilation | — | — | — | — | *no lock possible* — a shear test on one joint that reports no factor of safety. Its six vendor models all carry `include_dilation: no`, so the manual's dilation comparison never exercised dilation; XSLOPE's dilation and its peak and residual strengths are checked against their closed forms in `test/joint_element_check.py`. |

</div>

---

## The Rows

### 🟢 RJ-1a: Goodman & Bray block toppling, case a (rj001a) {#rj-1a}

Goodman & Bray's own toppling example, the section all four of problem 1's cases share: sixteen
rock columns 10 m wide and 4 to 40 m tall on a base that steps up at 30°, their sides at 120°,
normal to that base, and their tops cut off by a 56.6° face. The rock is elastic (γ = 25 kN/m³,
E = 20 GPa, ν = 0.3), so every mechanism the model has is a joint one, the idealization the closed
form makes. The joints are the shared edges of the column
outlines (see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint))
and carry no cohesion, φ = 38.15°, the angle this case is posed at, and the corpus's standard
stiffness pair.

The referee is Goodman & Bray's iterative column analysis, recomputed on this section. It
reproduces the method's published mode pattern and requires a horizontal toe force of 0.36 kN/m
for limit equilibrium, against the 0.5 kN the manual states. Carried into a strength reduction, it gives each of the four cases a
referee within a rounding of the value the manual prints for it.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.99 vs 1.0000 (−1.0%) | 0.99 (+2.8%) | 0.99 / 0.97 |

<!-- test: file=files/rocscience/joints/rj001a.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1a, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

**The excess over the closed form is one assumption, and it is where the thrust between two columns
acts.** Goodman & Bray hand each column-to-column thrust to the top corner of its contact. A contact
cannot do that: the two columns lean together, the face stays closed only over its upper quarter to
half, so the resultant stands a tenth to a sixth of the face
below the corner. Given those heights read off the solved state, and nothing else changed, the same
recursion returns **1.0279**, the bottom of this row's own bracket. The closed form's other three
assumptions the solution obeys exactly: every block balances in force and in moment on the
interface stresses alone, every closed side pair is at its friction limit, and the base reaction of
every toppling block sits on its downslope corner.

At the block width of 10 m the bracket reads 1.037; at 2D sizes of 7.0 m, 5.0 m and 3.5 m it reads
1.018, one step of the search lower, and at the 5.0 m of the lock all nine trials settle.

The manual's figure labels the side boundaries as rollers, but the model clamps them and the
restraints follow the model. The 0.5 kN toe force, negligible against the lowest block's weight,
is not carried; the vendor's reruns of all four cases move the force from the toe at (−0.5, 0.866)
to the block corner at (−2.5, 4.330), so their "with improvement" factors are a different load case.

**Input file:** [rj001a.xlsx](files/rocscience/joints/rj001a.xlsx).

![RJ-1a: Goodman & Bray block toppling, case a (rj001a) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own, so the whole figure is the joints: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope, and the deformed section shows the stack rotating forward over the face column by column while the toe column slides out along its own stretch of base](images/RJ-1a.png)

### 🟢 RJ-1b: Goodman & Bray block toppling, case b (rj001b) {#rj-1b}

Case a's section and rock with the joints at φ = 33.0239°, the lowest of the four, held up by the
2013 kN horizontal force the case is posed with — about twice the lowest column's own weight, and
what lets a stack on 33° joints stand where case a's needs 38°. It enters as a line load on the
`lloads` sheet, at the upper-left block corner, (−2.5, 4.330): the toe is the end of the basal
joint, where a load has no defined side to act on.

On this case the recursion requires a horizontal toe force of 2,012.86 kN/m for limit equilibrium,
against the 2,013 kN the manual states — 0.007% on a stack weighing 83,500 kN/m.

The thrust heights here are not case a's. The toe force pushes the two lowest columns back into
their own step risers and carries those two thrusts far down their faces, and given this case's own
measured heights the recursion returns **1.0078**, the bottom of this row's bracket.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.97 vs 1.0000 (−3.0%) | 0.99 (+2.8%) | 0.97 / 0.94 |

<!-- test: file=files/rocscience/joints/rj001b.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1b, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

The mesh does not move this row: the block width of 10 m and the 5.0 m of its lock return the same
factor, and at 5.0 m every trial settles within the iteration limit. On this case a few hundred
kilonewtons per meter of extra capacity is worth under two points of factor of safety, and the
interface change refinement makes is smaller than that.

Every input class matches the vendor model, including the side restraint and the force's magnitude
and direction.

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
returns **1.0469**, the bottom of that mesh's bracket — case a's mechanism and case a's
assumption, half a degree of friction further on.

**The mesh moves this row the same way it moves case a.** At the block width of 10 m the bracket
reads 1.057; at the 5.0 m of its lock it reads 1.037, one step of the search down and inside the
row's own tolerance. Every trial at 5.0 m settles within the iteration limit.

**Input file:** [rj001c.xlsx](files/rocscience/joints/rj001c.xlsx).

![RJ-1c: Goodman & Bray block toppling, case c (rj001c) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's at half a degree more friction: the same forward rotation of the column stack over the face, and the same slip on every column contact and along the stepped base](images/RJ-1c.png)

### 🟢 RJ-1d: Goodman & Bray block toppling, case d (rj001d) {#rj-1d}

Case c's joint friction angle with case b's 2013 kN force: the same stack, stabilized. It is the
strongest of the four, and the closed form and UDEC agree that it is.

Given the thrust heights read off this case's solved state, the recursion returns **1.2422**,
the bottom of this row's bracket. Across all four cases XSLOPE stands above the closed form and
RS2 below it, with UDEC between them, and on each of the four the recursion lands on the row's own
bottom of the row's own bracket once the thrust is put where the solution puts it rather than at
the corner of the contact.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.252** | 1.2308 (+1.7%) | 1.19 vs 1.2308 (−3.3%) | 1.22 (+2.6%) | 1.19 / 1.16 |

<!-- test: file=files/rocscience/joints/rj001d.xlsx, type=fem_ssrm, expected_fs=1.252, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1d, tier=gate, f_stand=1.2421875, f_fail=1.26171875, check=edges -->

The mesh does not move this row either: the block width of 10 m and the 5.0 m of its lock
return the same bracket, end for end, and at 5.0 m every trial settles within the iteration limit. Like case
b, this case is posed where the toe-force curve is steep, so the interface change that refinement
makes cannot carry it across a step of the search.

Every input class matches the vendor model, the force reaching the same point case b's does.

**Input file:** [rj001d.xlsx](files/rocscience/joints/rj001d.xlsx).

![RJ-1d: Goodman & Bray block toppling, case d (rj001d) — FEM inputs with the 2013 kN toe force, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The right-hand panels are the state past the critical factor: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope and at the toe, and the deformed section shows the stack rotating forward over the face column by column — case b's mechanism at the strength this case's toe force holds](images/RJ-1d.png)

### 🟢 RJ-2: Alejano & Alonso block toppling (rj002) {#rj-2}

A 9.85 m face at 58.65°, a basal joint rising from its toe at 30° as the stepped surface the columns
stand on, and twenty-two columns at 64°, 1.6 m apart, through the same toe. Both sets carry φ = 31°
and no cohesion. The rock is elastic (γ = 25 kN/m³, E = 20 GPa, ν = 0.3), the vendor's own
`Plasticity Specifications: Non`, so every mechanism the model has is a joint one.

The referee is Goodman & Bray's iterative column analysis, recomputed on the thirteen columns the
two sets cut out of the toppling mass: **0.7734**, against the 0.76 Alejano & Alonso publish for the
same method on the same problem.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **0.783** | 0.7734 (+1.2%) | 0.86 vs 0.7734 (+11.2%) | 0.87 (−10.0%) | 0.86 / 0.82 |

<!-- test: file=files/rocscience/joints/rj002.xlsx, type=fem_ssrm, expected_fs=0.783, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-2, f_stand=0.7734375, f_fail=0.79296875, check=edges, tier=gate -->

A step of refinement to 0.35 m moves the factor by one step of the search, inside the row's own
tolerance, and all nine trials decide on both meshes.

**The two programs that are not the closed form land together above it.** Alejano & Alonso publish
their own UDEC run at 0.87 beside their Goodman & Bray 0.76, and RS2 reports 0.86 and 0.82: all
three programs stand above the recomputed closed form, this row by 1.2% and the other two by a
tenth. Neither of the two can be interrogated on this problem: the manual states nothing about the
model's UDEC settings, no block rounding, no deformability and no stiffness, where the rigid-block
note it gives for problems 9 to 14 is the only statement of the kind it makes anywhere.

The joints total the vendor's own 311.4 m of trace at 64° and 30°, and every other input class
matches, including the clamped side restraint. Turned on, RS2's relief of a slipping joint's
stiffness returns this row's bracket unchanged.

**Input file:** [rj002.xlsx](files/rocscience/joints/rj002.xlsx).

![RJ-2: Alejano & Alonso block toppling (rj002) — FEM inputs, mesh, joint slip at the critical SRF and the section deformed 17×. The rock carries no strain of its own because it cannot yield, so every movement in the section is on a joint: slip gathers where the basal joint reaches the toe of the face, the columns standing on it open along their upper halves, and the deformed section shows them rotating out over the face while the rock below the basal joint stays put](images/RJ-2.png)

### ⊘ RJ-3: Lorig & Varona forward block toppling (rj003) {#rj-3}

Problems 3 to 6 are the toppling and plane-failure examples of Lorig & Varona (2004). The 260 m
section at 55° that problems 3 to 7 share, cut by two sets: columns at 70° at 20 m
spacing and a cross set at −20° at 30 m, both through the origin. The manual states the pair as
"70 and 160" degrees, which is the same two planes measured the other way round the half circle.
The rock is elastic, the vendor file setting its plasticity to none, so only the joints can fail;
γ = 26.0946 kN/m³, E = 9072 MPa, ν = 0.26. The joints carry c = 100 kPa and φ = 40°.

The manual's tables for problems 3 to 7 state only the slope geometry, the joint friction angle and
the rock's tensile strength, where the models of problems 4, 6 and 7 carry a rock strength of
c = 675 kPa and φ = 43° and all five carry a joint cohesion of 100 kPa.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| 1.213 *reported, no lock* | 1.13 (+7.3%) | 1.12 vs 1.13 (−0.9%) | 1.12 / 1.09 |

The row reports a factor but does not lock one. The search brackets 1.213 between 1.203, where
the slope stands, and 1.222656, where the trial does not settle within the 250,000-iteration limit:
the slope neither comes to rest nor runs away there, so the top of the bracket is a statement
about the iteration limit and the search cannot close. Every other trial settles. [Problem 5](#rj-5)
and [problem 20](#rj-20) are reported rather than locked because trials at the ends of their
brackets do not settle either.

The reported value stands above UDEC and above the vendor's default. The vendor's own two
solution schemes give 1.12 and 1.09 for this problem, so the vendor's answer moves with how it is
solved; this is one of the three problems the vendor reruns under its `Improve Joint
Convergence` option, described under [Methodology](#methodology).

Every input class matches the vendor model, including the side restraint the vendor clamps in both
directions.

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
bracket, end for end. Every trial settles on both meshes.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj004.xlsx](files/rocscience/joints/rj004.xlsx).

![RJ-4: Lorig & Varona flexural toppling (rj004) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. Here the rock can yield, and it does: a band of shear strain climbs from the toe across the columns, and drawn without exaggeration the columns are bent through that band rather than rotated about it, which is what separates flexural toppling from the block toppling of problem 3](images/RJ-4.png)

### ⊘ RJ-5: Lorig & Varona backward block toppling (rj005) {#rj-5}

The shared 260 m section cut by two sets: one at −55° at 10 m spacing through the toe, dipping out
of the face so the blocks lean back rather than forward, and a horizontal set at 40 m spacing. The
rock is elastic, its plasticity set to none in the vendor file: γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26; the joints carry
c = 100 kPa and φ = 40°.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| 1.818 *reported, no lock* | 1.7 (+6.9%) | 1.65 vs 1.7 (−2.9%) | 1.65 / 1.86 |

This row is the corpus's hardest to settle. Four of its nine trials do not settle within the
250,000-iteration limit, and two of the four are the ends of the bracket, 1.808594 and 1.828125, so
the search cannot close: the row reports the bracket's midpoint but does not lock it.

The reported value stands above UDEC and above the vendor's default, and between the vendor's
own two numbers: its two solution schemes give 1.65 and 1.86 on the same model, the widest spread
in the manual. The vendor needed its `Improve Joint Convergence` option to rerun
this problem, and it reaches a different answer with it.

A step of refinement to 8.4 m settles two more of the nine and brackets a factor one step lower. Two
trials still do not settle there, one of them an end of that bracket, so neither mesh closes the
search. Together the two runs take about 16 hours, more than twice any other row here.

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
| 1.271 *reported, no lock* | 1.27 (+0.1%) | 1.25 vs 1.27 (−1.6%) | 1.25 / 1.31 |

The row reports a factor but does not lock one, and what holds it back is the mesh rather than
the iteration limit. Every trial settles, but a step of refinement to 8.4 m moves the factor by two
steps of the search, where the row's tolerance is one.

The reported value matches UDEC and sits between the vendor's two numbers, 1.25 and 1.31: the
vendor's own answer moves with its solution scheme on this model by more than the mesh moves
XSLOPE's.

It is the only row in the corpus whose factor moves with the mesh: repeatable, but an answer that
belongs to the mesh rather than to the slope. Every other row that states a refinement step either
returns the same bracket, end for end, or moves inside its own tolerance.

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
bracket, end for end. Every trial settles on both meshes.

Every transcribed input class matches the vendor model, including the side restraint. The manual's
own table for this problem prints the slope angle as 5°; the figure and the model are the same 55°
slope problem 6 uses.

**Input file:** [rj007.xlsx](files/rocscience/joints/rj007.xlsx).

![RJ-7: plane failure with non-daylighting discontinuities (rj007) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. No joint reaches the face at a shallower angle than the face itself, so the failure cannot slide out along one: a strong band of shear strain cuts across the steep joints from the toe, and the mass above it moves out over rock it has had to break](images/RJ-7.png)

### 🟢 RJ-8: Flexural toppling in a base friction model (rj008) {#rj-8}

Pritchard & Savigny's base-friction table model, scaled up a hundred times: a 30.5 m face at 78°
with twelve columns at −60°, 5.08 m apart, standing on a horizontal basal joint and closed at the
back by a vertical joint. The three highest columns end on that back joint, as the vendor's own do,
and the network matches the vendor's totals to the millimeter. The vendor file stores the columns
clipped to the block the basal and back joints bound; generated over the whole section, the same set
would be sixteen chains, three of them below a plane the model has no columns under and three more
running on past x = 68.4 into the strip the vendor leaves uncut. The rock
is Mohr-Coulomb (γ = 25.506 kN/m³, E = 22.771 GPa, ν = 0.139, c = 60 kPa, φ = 39°). The joints
carry no cohesion, φ = 39°, and the softest normal stiffness in the corpus bar one:
k<sub>n</sub> = 1.5 × 10<sup>7</sup> kPa/m against the set's usual 10<sup>8</sup>.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **0.764** | 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | 0.75 / 0.75 |

<!-- test: file=files/rocscience/joints/rj008.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-8, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement to 1.05 m moves the factor by one step of the search, inside the row's own
tolerance, and all nine trials settle on both meshes. The two meshes disagree about one trial, the
one the search closes on: at F = 0.754 the slope stands on the corpus mesh and fails on the finer
one. This row settles fastest of any here, where the [geotextile wall family](rs2.md#rs2-48) runs
out of its iteration limit on five trials across its eight rows.

The model's 75 kPa tensile strength sits above the Mohr-Coulomb apex its own c and φ imply, so the
cap never binds.

**Input file:** [rj008.xlsx](files/rocscience/joints/rj008.xlsx).

![RJ-8: flexural toppling in a base friction model (rj008) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The strain gathers into one lobe per column along a band that climbs from the toe across the stack; above that band the columns are visibly bent rather than merely tilted, and the rock below it is unstrained. That is the break surface of flexural toppling rather than sliding along any one joint](images/RJ-8.png)

### 🟢 RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) {#rj-9}

A 50 m slope at 50° cut by bedding dipping **out of the face** at −50° at 3 m spacing, φ = 30°, with
a two-segment release trace at the toe at φ = 40° that undercuts the lowest slab. Two joint
strengths, which the manual's own geometry table prints in the reverse order from its RS2 legend;
on all six of problems 9 to 14 the bedding network and the short crest traces carry different
friction angles, so a single quoted joint friction angle for them is incomplete. The rock is elastic
at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3, the manual's rigid-block device.

The face and the bedding dip at the same angle, so no slab can slide out along a single plane: the
mechanism is the bilinear one the problem is named for, sliding on a basal plane combined with
sliding along the release that the face undercuts. The release trace ends on the bedding plane it
belongs on; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | 0.40–1.45 | 1.01 / 1.09 |

<!-- test: file=files/rocscience/joints/rj009.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-9, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

A step of refinement to 2.1 m returns the same bracket, end for end, and every trial settles on both
meshes.

RS2's own two factors straddle this row — 1.01 without the joint improvement option and 1.09 with
it — which is true of only three other problems in this corpus.

The paper's own limit equilibrium for this problem is a range rather than an answer — 0.40 to 1.45,
as the manual prints it — and it is recorded beside the referee rather than scoring, as every source
limit equilibrium here is. The distinct-element run is what the manual verifies RS2 against.

**Input file:** [rj009.xlsx](files/rocscience/joints/rj009.xlsx).

![RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: the release trace under the crest carries the brightest slip, a bedding plane below the toe carries the rest, and the deformed section shows the slab between them sliding out over the bench](images/RJ-9.png)

### 🟢 RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) {#rj-10}

Example 1a with the release joint moved five meters up the face, which the manual says is the whole
difference between the two problems. The section, the bedding, both friction angles and the rock
are [problem 9](#rj-9)'s.

Moving the release upslope leaves the toe undercut where it was and cuts a second block out of the
face above it: the trace that ran from the face to the toe in example 1a stays, and a new one leaves
the face five meters higher. Both end on the bedding plane they are released by;
see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | 0.43–1.45 | 0.92 / 1.08 |

<!-- test: file=files/rocscience/joints/rj010.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-10, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

A step of refinement to 2.1 m returns the same bracket, end for end, and every trial settles on both
meshes.

RS2's own two factors straddle this row — 0.92 without the joint improvement option and 1.08 with it,
around 1.037 — as they do on [problem 9](#rj-9), [problem 5](#rj-5) and [problem 6](#rj-6) and on no
other problem here. The paper's own limit equilibrium is again a range rather than an answer, 0.43 to
1.45, recorded beside the referee rather than scoring.

**Input file:** [rj010.xlsx](files/rocscience/joints/rj010.xlsx).

![RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) — FEM inputs, mesh, joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and, as in example 1a, almost none of it slips: the release trace high on the face carries the brightest slip, one bedding plane below the toe carries the rest, and the deformed section shows the block between them moving out over the bench](images/RJ-10.png)

### 🟢 RJ-11: Alejano et al. plowing sliding slab failure (rj011) {#rj-11}

A 25 m slope at 50° with bedding dipping out of the face at −50° at 1.5 m spacing, φ = 30°, and two
release traces at φ = 20°: one at the toe running below the bench and one from the face down onto
the bedding plane that releases the toe block. Plowing failure is the paper's name for sliding on
a primary discontinuity combining with sliding on a joint sub-parallel to the face, which lifts the
toe block and eventually rotates it out of the slope. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The release traces end on the bedding plane
they belong on; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | Rigid-block bound referee | RS2 vs referee | UDEC | Alejano Eq. (7) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.213** | 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | 1.21 (+0.2%) | 1.7582 · the paper prints 1.75 | 1.22 / 1.3 |

<!-- test: file=files/rocscience/joints/rj011.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-11, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

A step of refinement to 1.05 m returns the same bracket, end for end, and every trial settles on both
meshes.

**This is the problem where Alejano's Eq. (7) is not a rigid-block answer.** The two blocks the
mechanism cuts out — the slab and the toe block below it — admit no set of contact forces inside
their friction cones above 1.2148, whatever the distribution; Eq. (7) returns 1.7582, half a factor
above that, so it cannot be what scores this row and the bound is the referee instead. XSLOPE reads
1.213, RS2 1.22 and the paper's own UDEC run 1.21 — three codes inside 0.8% of one another, all
three on the bound. [Problem 15](#rj-15) has the same shape at a problem where no bound is
available.

Every input class matches the vendor model, including the release traces' endpoints and the side
restraint the vendor clamps in both directions.

**Input file:** [rj011.xlsx](files/rocscience/joints/rj011.xlsx).

![RJ-11: Alejano et al. plowing sliding slab failure (rj011) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. One bedding plane from the crest to the toe carries almost all of the slip, and at its foot the two release traces cut out a small wedge: the deformed section shows that wedge lifted and rotated out over the bench while the slab above it slides down the plane, which is the plowing mechanism the paper names. The wedge is driven out and up by the slab above it, and at the panel's exaggeration a movement of centimeters draws as meters, so the block appears to leave the slope](images/RJ-11.png)

### 🟡 RJ-12: Alejano et al. plowing toppling slab failure (rj012) {#rj-12}

A 25 m slope at 60° cut by bedding dipping out of the face at −60° at 1.5 m spacing, φ = 30°, with
two short release traces at φ = 40°: one at the toe running below the bench, and one from the face
down onto the bedding plane that releases the toe block. It is problem 11's plowing mechanism, but
at 60° the rotation of the toe block rather than the sliding governs. The rock is the family's
elastic rigid-block stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3.

Alejano's Eq. (7), a moment balance about the toe for exactly this mechanism, is the referee.
Recomputed on the inputs the vendor file states, it gives **1.9659** against the 2.00 the paper
prints for the same example, and the same implementation reproduces the paper's own table on its
other worked examples.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **2.033** | 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | 2.0324 | 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 |

<!-- test: file=files/rocscience/joints/rj012.xlsx, type=fem_ssrm, expected_fs=2.033, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-12, f_stand=2.0234375, f_fail=2.04296875, check=edges, tier=gate -->

The release traces end on the bedding plane they belong on (see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint)).

Two steps of refinement, to 1.05 m and 0.735 m, return the same bracket, end for end, and every
trial settles on all three meshes.

**The three independent answers for this problem do not agree with one another, and this row is the
one nearest the closed form.** Alejano's limit equilibrium is 1.9659 recomputed, the paper's UDEC run is 1.78 and
RS2 reads 1.39, on a mechanism whose governing mode — the toe block rotating out rather than
sliding — every source agrees on.

**Input file:** [rj012.xlsx](files/rocscience/joints/rj012.xlsx).

![RJ-12: Alejano et al. plowing toppling slab failure (rj012) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section, and only a few of its traces carry any slip: one release trace under the crest and the bedding beneath the toe block, which is the plowing pair. The two right-hand panels are the last standing trial of the bracket, the state below the factor rather than past it — on a rock this stiff the model moves by microns until it does not, so the deformed section is drawn at tens of thousands of times scale. The capture past the factor is not drawn: it was stopped in its first iteration, before the section had moved at all](images/RJ-12.png)

### 🟢 RJ-13: Alejano et al. plowing sliding slab, example 4 (rj013) {#rj-13}

A 25 m slope at 55° cut by bedding dipping out of the face at −55° at 1.5 m spacing, φ = 25°, with
two release traces at φ = 20°: one at the toe running below the bench and one from the face down
onto the bedding plane that releases the toe block. It is problem 11's mechanism at a steeper
bedding and a weaker one — sliding on the primary discontinuity combining with sliding on a joint
sub-parallel to the face, which lifts the toe block. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The release traces end on the bedding plane they
belong on; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

Alejano's Eq. (7) recomputed on this problem's inputs gives **1.0002**, its sliding mode governing,
and the paper prints 1.00 for the same example. It sits on the rigid-block bound for these two
blocks, which is 0.9988, so the closed form and statics agree to three figures here.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **0.998** | 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | 0.9988 | 1.0 (−0.2%) | 1.0 / 1.05 |

<!-- test: file=files/rocscience/joints/rj013.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-13, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

A step of refinement to 1.05 m returns the same bracket, end for end, and every trial settles on both
meshes.

Every input class matches the vendor model, including the release traces' endpoints and the side
restraint the vendor clamps in both directions.

**Input file:** [rj013.xlsx](files/rocscience/joints/rj013.xlsx).

![RJ-13: Alejano et al. plowing sliding slab, example 4 (rj013) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it carries slip: one bedding plane from the crest down to the toe does, and the two release traces at its foot cut out the toe block, which the deformed section shows lifted and pushed out over the bench while the slab above it slides down the plane](images/RJ-13.png)

### 🟢 RJ-14: Alejano et al. plowing sliding slab, example 5 (rj014) {#rj-14}

Example 4's section at 60° with the two joint strengths the other way round: bedding dipping out of
the face at −60° at 1.5 m spacing at φ = 20°, the weakest bedding of the six, and two release traces
at φ = 30°. The rock is the family's elastic rigid-block stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³,
ν = 0.3.

Alejano's Eq. (7), recomputed on the inputs the vendor file states, gives **1.2034** for this
example, under the rigid-block bound for its two blocks. It is the referee.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.232** | 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | 1.2686 | 0.9 (+36.9%), 0.9994 with the paper's corner rounding corrected · the paper prints 1.00 | 0.89 / 1.09 |

<!-- test: file=files/rocscience/joints/rj014.xlsx, type=fem_ssrm, expected_fs=1.232, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-14, f_stand=1.22265625, f_fail=1.2421875, check=edges, tier=gate -->

**The 1.00 the paper prints for this example does not follow from the inputs it prints.** Every
input Eq. (7) reads matches the vendor's model — the toe block's base to five figures, 3.589 m
against the 3.588 m printed — and on them the equation returns 1.2034. For it to return 1.00 the toe
block's base would have to be 6.08 m, which is a different problem. The manual's other value,
UDEC's 0.9, is corner rounding by the paper's own account: with the rounding radius reduced it
reports 0.9994 for the same model.

Two steps of refinement, to 1.05 m and 0.735 m, return the same bracket, end for end, and every
trial settles on all three meshes. Every input class matches the vendor model, including the side
restraint the vendor clamps in both directions.

**Input file:** [rj014.xlsx](files/rocscience/joints/rj014.xlsx).

![RJ-14: Alejano et al. plowing sliding slab, example 5 (rj014) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: one bedding plane from the crest to the toe carries the slip, with the release trace at the toe opening as the slab above it slides out over the bench](images/RJ-14.png)

### 🟢 RJ-15: Partially joint-controlled footwall slope (rj015) {#rj-15}

A 40 m footwall at 40° whose bedding dips in the same direction at the same angle, 2 m apart, so
the slabs lie parallel to the face and a failure has to break rock at the toe to get out. This is
the one problem in the Alejano family whose rock can yield: Mohr-Coulomb, c = 200 kPa, φ = 35°,
γ = 28 kN/m³, E = 1 GPa, ν = 0.3. The joints are the corpus's softest — k<sub>n</sub> = 5 × 10<sup>6</sup>
kPa/m and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m, twenty times below the set's standard pair — with
no cohesion and φ = 25°.

The manual's table prints the slope height as 25 m; the vendor's section rises 40 m from the toe at
(0, 0) to the crest at (−47.6701, 40), the height the source paper states, and the section is what
is built.

Alejano's footwall limit equilibrium, Eqs. (9)–(10), recomputed on the vendor file's inputs at the
optimum the paper states (a break-out inclined 14° to the bedding and
emerging at 55°), gives **1.7985** against the value the paper prints, and minimized over its own
two angles it settles within half a degree of that optimum. The referee is Rocscience's own Slide2
limit-equilibrium search at 1.25.

| XSLOPE SSRM | Slide2 LE search referee | RS2 vs referee | Alejano Eqs. (9)–(10) | UDEC-SSRT (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.271** | 1.25 (+1.7%) | 1.28 vs 1.25 (+2.4%) | 1.7985 (single-bed formula; the paper prints 1.72) | 1.6 (−20.6%) | 1.28 / 1.42 |

<!-- test: file=files/rocscience/joints/rj015.xlsx, type=fem_ssrm, expected_fs=1.271, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-15, f_stand=1.26171875, f_fail=1.28125, check=edges, tier=gate -->

**The closed form prices one mechanism, and every program free to search for a surface finds a
weaker one.** Eqs. (9)–(10) drive a wedge out through a single slab and hard-wire the break-out to
cross exactly one 2 m bed; their factor rises with that thickness, so the shallowest case the
formula allows is also the lowest it can report. This SSRM at 1.271, RS2's SSR at 1.28 and the
Slide2 search at 1.25 are under no such restriction and sit inside 2.4% of one another, far below
the closed form. Both finite element codes miss it by the same amount and in the same direction, so
the gap belongs to what the formula is allowed to consider rather than to either program, and the
searched limit-equilibrium answer is the referee, as the rigid-block bound is on
[problem 11](#rj-11).

A step of refinement to 1.4 m moves the factor by one step of the search, inside the row's own
tolerance, and every trial settles on both meshes.

Every input class matches the vendor model, including the side restraint; the 1000 kPa tensile cap
sits above the Mohr-Coulomb apex and never binds.

**Input file:** [rj015.xlsx](files/rocscience/joints/rj015.xlsx).

![RJ-15: partially joint-controlled footwall slope (rj015) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The bedding slips over a long stretch behind the face, and the rock's only strain is a small patch at the toe where the slab has to break through to get out — the coupled mechanism the source paper describes](images/RJ-15.png)

### 🔴 RJ-16: Barla et al. tilt-table block toppling (rj016) {#rj-16}

A laboratory test rather than a slope: fourteen columns of 9 cm blocks stacked into a 63.4°
staircase on a plate, and the plate tilted until the stack topples; the problem reports the angle
at which it goes. The blocks are elastic (E = 350 MPa, ν = 0.2, γ = 28 kN/m³) on an effectively
rigid plate, and the joints carry no cohesion, φ = 38°, and the corpus's softest stiffness
pair. The 0° model runs on rollers; the nine tilted models pin every exterior node in both
directions and converge to a tolerance two orders tighter, and the plate carries no body force in
any of the ten.

The vendor tilts the model itself, one file per degree. XSLOPE turns the load instead: for a tilt θ, a
horizontal seismic coefficient k = tan θ toward the face, raised at full strength (F = 1) until the stack stops
standing, on the corpus mesh of 0.09 m, the block size.

| XSLOPE tilt | UDEC referee | RS2 vs referee | Experiment | Goodman & Bray, plate tilted | RS2 without / with improvement |
|---|---|---|---|---|---|
| **10.24°** | 11° (−6.9%) | 9° vs 11° (−18.2%) | 9° | 7.6° | 9° / 7° |

<!-- test: file=files/rocscience/joints/rj016.xlsx, type=fem_tilt, k_stand=0.179846, k_fail=0.181396, expected_tilt=10.24, tolerance=0.1, element_type=tri6, target_size=0.09, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-16 -->

The stack stands at k = 0.1798, which is 10.20°, and goes at k = 0.1814, which is 10.28°. It stands
at every coefficient tried below that and fails at every one above, through 0.25, and the row's
check re-solves both coefficients at full strength. Run on with no stopping rule at all, the plain
iteration agrees: the stack settles at every tilt through 9.9° and moves at a steady rate from
10.2°.

**A coefficient is not a rotation, and here the difference is measurable and measures zero.**
Tilting the model by θ turns the body force through θ and leaves its magnitude at γ; a coefficient
k = tan θ turns it through the same angle and multiplies its magnitude by 1/cos θ, which is 1.6% at
10.2°. With zero joint cohesion and no block able to yield, the state depends on the body force's
direction and not its size: re-solving the two bracketing coefficients with every unit weight
scaled by 1.016 and then by 0.5 leaves both unchanged.

**Where the five numbers sit.** Goodman & Bray's column analysis, with the plate tilted, puts the
toppling tilt of these rigid columns at 7.6°. It allows no contact pressure below a column's
corner, so it is a lower bound, and everything else sits above it: the physical
stack and RS2 at 9°, XSLOPE at 10.2°, and the distinct-element code, whose contacts roll and re-form
as the columns lean, at 11°.

The vendor's weightless plate is built at the 27 kN/m³ its property row states, because
`build_fem_data` requires a positive unit weight.

**Input file:** [rj016.xlsx](files/rocscience/joints/rj016.xlsx).

![RJ-16: Barla et al. tilt-table block toppling (rj016) — FEM inputs with the seismic coefficient that stands for the tilt, mesh with the vendor's rollers, joint slip at the first coefficient the stack goes at, and the deformed section. The slip is on the vertical joints between the lower columns and along the bedding one block above the plate, and the deformed section at 81x shows what that adds up to: every column leaning downslope about its own base, the tall ones at the back furthest over, which is toppling rather than the stack sliding along the plate](images/RJ-16.png)

### 🟡 RJ-17: Step-path failure, en-echelon joints (rj017) {#rj-17}

[Problem 18](#rj-18)'s section and Mohr-Coulomb rock (c = 25 kPa, φ = 25°, no tensile capacity),
cut by three joints at 36.1° that stop short of one another instead of running through, at problem
18's joint strength, c = 1 kPa and φ = 35°. The rock bridges between them are what a step-path
failure has to break through, and they are why this slope stands where problem 18's continuous
joints let it go.

The strength reduction is confined to the vendor's SSR search area, a rectangle drawn on the
manual's figure but absent from its tables, outside which RS2 holds every element linear elastic,
1,173 of the model's 2,579; XSLOPE states it as a polygon overlay classified by the same
element-centroid test. No closed form exists for a step path through
rock bridges, so the referee is the manual's single distinct-element run.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.213** | 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | 1.24 / 1.2 |

<!-- test: file=files/rocscience/joints/rj017.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-17, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

A step of refinement to 0.7 m moves the factor by one step of the search, inside the row's
tolerance, to 1.193, and every trial settles on both meshes.

**What decides this row is the rock bridges.** The rock has no tensile capacity at all, so the
intact ligaments between the joint segments carry nothing across them, and the strength reduction
takes their 25 kPa cohesion down alongside the joint friction: the factor is the strength at which
three short bridges shear through. Both finite element codes read that below the distinct-element
run and on the same side of it, 2.1 points apart.

Transcribed instead as the staircase the vendor's mesh rasterized the rectangle into, the zone holds
693 elements of this mesh elastic where the rectangle holds 694, and the bracket comes back
identical trial for trial.

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
elements to 6,707 and 66 — does not move the factor at all, and every trial on both meshes
settles.

**Input file:** [rj018.xlsx](files/rocscience/joints/rj018.xlsx).

![RJ-18: step-path failure through three continuous joints (rj018) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Every joint is slipping along its lower half and open along its upper, and the three slabs between them slide out together down the 36.1° path; the rock itself carries almost no plastic strain](images/RJ-18.png)

### 🟢 RJ-19: Bi-planar step-path failure (rj019) {#rj-19}

One Mohr-Coulomb rock (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 10,500 kPa, φ = 35°, tensile capacity
200 kPa) cut by two discontinuous joints with a rock bridge between them: a basal joint at 28.4° and
an upper joint at 56.3°. Both carry c = 0, φ = 40° and the same stiffness pair as RJ-18.

**What decides this row is a tensile cap rather than a cohesion.** The two joints leave a rock
bridge 1.414 m long between their tips. At the rock's 10,500 kPa cohesion that bridge carries
14,849 kN/m in shear, which is 1.2 times the whole 11,929 kN/m sliding mass and 2.6 times its
component down the basal joint, so it cannot be sheared through; what lets the block move is the
same rock's 200 kPa tensile capacity, worth 283 kN/m across the bridge. This row reduces that
capacity with the trial factor, as every vendor model does.

The referee is the rigid-block limit equilibrium of that mechanism, recomputed from the inputs the
workbook carries: the block slides on the basal joint, opens the upper joint behind it, and pulls
the bridge apart in tension, with the joint's friction and the rock's tensile capacity both reduced
by the trial factor. It holds up to 1.5914; friction alone holds it to 1.551. XSLOPE's 1.623 sits
2.0% above 1.5914, and both vendor numbers sit below what the statics of the block admits, so the
distinct-element run and RS2 fail this block by some route those statics do not contain.

| XSLOPE SSRM | Rigid-block limit equilibrium referee | RS2 vs referee | UDEC (manual) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.623** | 1.5914 (+2.0%) | 1.5 vs 1.5914 (−5.7%) | 1.46 (+11.2%) | 1.5 / 1.41 |

<!-- test: file=files/rocscience/joints/rj019.xlsx, type=fem_ssrm, expected_fs=1.623, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=true, k0=1, benchmark=RJ-19, f_stand=1.61328125, f_fail=1.6328125, check=edges, tier=gate -->

A step of refinement to 2.1 m moves the factor by one step of the search, inside the row's own
tolerance, but on that mesh the trial at the top of the bracket does not settle within the
iteration limit. The row is locked on the corpus mesh, where every trial settles.

The joint inclinations are the model's: the manual's table prints 59°, where its figure dimensions
56° and 28° and the model's own endpoints give 56.3° and 28.4°.

**Input file:** [rj019.xlsx](files/rocscience/joints/rj019.xlsx).

![RJ-19: bi-planar step-path failure with a rock bridge (rj019) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Both joints have opened, the basal one is slipping at its lower end, and the only strain in the rock is the patch at the bridge between the two joint tips, where the block above has to break through to move](images/RJ-19.png)

---

### ⊘ RJ-20: Hammah & Yacoub Voronoi slope (rj020) {#rj-20}

An 80 × 70 m section with a 60 m face at 71.6°, tessellated into Voronoi blocks over the whole of it.
The rock is Mohr-Coulomb (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 1,000 kPa, φ = 35°, no tensile
capacity) and every block wall is a joint at c = 500 kPa, φ = 20° with the corpus's standard
stiffness pair. The source paper asks how the failure of a slope in blocky rock changes with the
scale of its blocks.

The tessellation was generated in UDEC and imported, with no block size, density or seed
published for it; the model file holds it as 523 joint boundaries, 1,177 segments, the only
statement of the network there is, and they are transcribed verbatim. Their mean block width, 2.895 m, is this
row's mesh size, standing in for the joint spacing every other row meshes at.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| *no lock* | 2.46 | 2.21 vs 2.46 (−10.2%) | 2.21 / 2.37 |

The row prints no factor, and what holds it back is the iteration limit. Four of the nine trials do
not settle within the 250,000 iterations allowed, and both of the trials the search closes on are
among them, so a factor read off this bracket would be a statement about the iteration limit rather
than about the slope. The other five settle: three stand and two fail. The search, and the solve
past the critical factor its figure is drawn from, take 3.7 hours between them, the longest run
here. The figure shows the mass failing on a surface picked through the block walls rather than
along any one plane it contains, which is the observation the source paper is about.

Separating the block size from this particular tessellation, the paper's own question, takes a
second network of the same blocks drawn another way. The one `xslope.joints.voronoi` draws at that
block size stops in the mesher: the generator keeps a trace down to a thousandth of the section
diagonal, 0.106 m here, and a trace that short is pulled onto its own junction, where the line it
leaves has no length left. Three seeds and two block sizes stop the same way, so this row is scored
on the vendor's own tessellation alone.

**Input file:** [rj020.xlsx](files/rocscience/joints/rj020.xlsx).

![RJ-20: Hammah & Yacoub Voronoi slope (rj020) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The slip runs from the toe up through the block walls on a curved path to the crest, taking whichever wall of each block lies nearest that line, and the rock carries almost no strain except a patch at the toe where the path turns: the mass fails on a surface picked out of the tessellation rather than along any one joint in it](images/RJ-20.png)
