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
  `.fez` files, not the manual's tables, which carry errata; each row notes where its model departs
  from the manual. Every row's geometry, materials, joint network, restraints and loads match the
  vendor's model unless the row's notes say otherwise. The vendor's MPa and MN/m³ are converted to kPa and kN/m³.
- **Transcription.** In transcribing the vendor's files to XSLOPE models, the following
  decisions were made:
    - *Strength after failure.* No problem states a residual joint strength or dilation angle, and
      listed material residuals equal the peak, so joints and rock keep their full strength after
      failure.
    - *Two vendor settings not carried over.* Every file reduces a slipping joint's stiffness a
      hundredfold (`joint_stiffness_factor: 0.01`) and divides the rock's tensile cap by the trial
      factor, the factor that each run of the strength reduction divides the strengths by. The corpus does neither: the stiffness reduction moves brackets by one to four steps
      on the [RS2 corpus](rs2.md)'s walls and embankments, and the cap governs only on
      [problem 19](#rj-19), where it is applied. Turned on, the stiffness reduction leaves the bracket
      unchanged on every row of this page it has been run on: problems 1, 2, 8 to 15 and 17 to 19.
    - *Joint strength on problems 3 to 7.* Their source, Lorig & Varona's chapter in Wyllie & Mah
      (2004), gives the joints friction only; the vendor's files add 100 kPa of cohesion. Problems
      3, 4, 5 and 7 use the chapter's joints because the UDEC results they are scored against come
      from the chapter. The rock stays as the vendor has it, elastic on 3 and 5, because with the
      chapter's yielding rock problem 5 reproduces neither the chapter's answer nor its mechanism.
      RS2's numbers used the 100 kPa and are shown beside each row for comparison only. Each row
      also shows XSLOPE on the vendor's file as given and on the chapter's full list.
- **Referee.** A closed-form rigid-block solution, recomputed from the model's own inputs, is the
  referee, the one value each row is scored against, where one exists: Goodman & Bray's column analysis on problems 1 and 2, Alejano's plowing
  equation on 11 to 14, the sliding block with its tensile bridge on 19 and, on 15, the vendor's
  Slide2 search, since Alejano's footwall equations describe a mechanism the slope does not take.
  Otherwise the referee is the one the manual names, in every case UDEC, a distinct-element program
  that models the rock as separate blocks in contact. Each row shows the
  recomputed value beside the printed one.
- **Rigid-block bound.** On the plowing problems, the highest reduction factor at which any
  admissible set of joint forces can hold both rigid blocks in place is a ceiling on any
  rigid-block answer. Alejano's Eq. (7) assumes one mechanism and sits on or below that ceiling on
  problems 12 to 14, where it is the referee; on problem 11 it exceeds the ceiling (1.76 against
  1.21), so the ceiling is the referee there.
- **RS2's two factors.** The manual reports each problem with and without the vendor's
  `Improve Joint Convergence` option. The run without it is the vendor's default and the same
  method as XSLOPE's, so it is the yardstick; both are recorded, and their spread is the width of
  the vendor's own answer.
- **Mesh.** Every row is meshed at its joint spacing, one element across the rock between
  discontinuities; problem 1's cases at 5 m, half their 10 m block width. No row on this page
  changes its factor with a finer mesh.
- **Search.** Bisection on the bracket between the highest factor at which the slope stands, coming
  to rest, and the lowest at which it fails, still moving. Each trial, one run at a single factor, is allowed 250,000 iterations, and up to a million if
  it is still slowing. A trial that reaches the million with no movement left except a few
  contacts repeating the same cycle stands; one that stops moving without balancing its forces and
  without such a cycle leaves the bracket unconfirmed. See [running a jointed
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
| [2](#rj-2) | 🟢 | Alejano & Alonso block toppling | SSRM 0.783 vs Goodman & Bray 0.7734 (+1.2%) | 0.86 vs 0.7734 (+11.2%) | UDEC 0.87 (−10.0%) | 0.86 / 0.82 | XSLOPE, RS2 and UDEC all give factors above the closed form, XSLOPE by about a percent and the other two by about a tenth. The manual states no UDEC settings for this model. |
| [3](#rj-3) | 🟢 | Lorig & Varona forward block toppling | SSRM 1.115 vs UDEC 1.13 (−1.3%) | 1.12 vs 1.13 (−0.9%) | — | 1.12 / 1.09 | Joints friction only, as the source chapter lists them. XSLOPE on the vendor's file as given, joints at c = 100 kPa: 1.232. RS2's two factors were computed with the 100 kPa. |
| [4](#rj-4) | 🟡 | Lorig & Varona flexural toppling | SSRM 1.232 vs UDEC 1.3 (−5.2%) | 1.19 vs 1.3 (−8.5%) | — | 1.19 / 1.27 | Joints friction only, as the source chapter lists them. XSLOPE on the vendor's file as given, joints at c = 100 kPa: 1.311. RS2's two factors were computed with the 100 kPa. |
| [5](#rj-5) | 🟡 | Lorig & Varona backward block toppling | SSRM 1.799 vs UDEC 1.7 (+5.8%) | 1.65 vs 1.7 (−2.9%) | — | 1.65 / 1.86 | Joints friction only, as the source chapter lists them. XSLOPE on the vendor's file as given, joints at c = 100 kPa: 1.857. RS2's two factors were computed with the 100 kPa. |
| [6](#rj-6) | 🟢 | Plane failure, daylighting | SSRM 1.271 vs UDEC 1.27 (+0.1%) | 1.25 vs 1.27 (−1.6%) | — | 1.25 / 1.31 | |
| [7](#rj-7) | 🟢 | Plane failure, non-daylighting | SSRM 1.486 vs UDEC 1.5 (−0.9%) | 1.57 vs 1.5 (+4.7%) | — | 1.57 / 1.59 | Joints friction only, as the source chapter lists them. XSLOPE on the vendor's file as given, joints at c = 100 kPa: 1.564. RS2's two factors were computed with the 100 kPa. |
| [8](#rj-8) | 🟢 | Flexural toppling, base friction model | SSRM 0.764 vs UDEC 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | — | 0.75 / 0.75 | |
| [9](#rj-9) | 🟢 | Bilinear slab failure, example 1a | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | LE (Alejano) 0.40–1.45 | 1.01 / 1.09 | |
| [10](#rj-10) | 🟢 | Bilinear slab failure, example 1b | SSRM 1.037 vs UDEC 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | LE (Alejano) 0.43–1.45 | 0.92 / 1.08 | |
| [11](#rj-11) | 🟢 | Plowing sliding slab failure | SSRM 1.213 vs rigid-block bound 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | UDEC 1.21 (+0.2%) · Alejano Eq. (7) 1.7582 | 1.22 / 1.3 | Alejano's Eq. (7) returns a factor above the bound on this problem, so the bound is what scores it; all three programs sit on the bound. |
| [12](#rj-12) | 🟡 | Plowing toppling slab failure | SSRM 2.033 vs Alejano Eq. (7) 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | UDEC 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 | |
| [13](#rj-13) | 🟢 | Plowing sliding slab, example 4 | SSRM 0.998 vs Alejano Eq. (7) 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | UDEC 1.0 (−0.2%) · Alejano prints 1.0 | 1.0 / 1.05 | |
| [14](#rj-14) | 🟢 | Plowing sliding slab, example 5 | SSRM 1.232 vs Alejano Eq. (7) 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | UDEC 0.9 (+36.9%) · Alejano prints 1.00 | 0.89 / 1.09 | The 1.00 the paper prints for this example does not follow from the inputs it prints; Eq. (7) on them gives 1.2034. |
| [15](#rj-15) | 🟢 | Partially joint-controlled footwall | SSRM 1.271 vs Slide2 LE search 1.25 (+1.7%) | 1.28 vs 1.25 (+2.4%) | Alejano Eqs. (9)–(10) 1.7985 (single-bed formula; the paper prints 1.72) · UDEC 1.6 (−20.6%) | 1.28 / 1.42 | Alejano's closed form drives a wedge out through a single 2 m bed and its factor rises with bed thickness; the three programs free to search for a surface agree at 1.25–1.28, so the limit-equilibrium search is the referee, as the rigid-block bound is on problem 11. |
| [16](#rj-16) | 🔴 | Barla et al. tilt-table block toppling | Tilt 10.24° vs UDEC 11° (−6.9%) | 9° vs 11° (−18.2%) | Experiment 9° · Goodman & Bray, plate tilted, 7.6° | 9° / 7° | Scored as a tilt angle rather than a factor of safety. With the tilt represented by a seismic coefficient and the strengths unreduced, the stack stands at 10.20° and topples at 10.28°. That is above Goodman & Bray's rigid-column lower bound and the physical test, and below UDEC. |
| [17](#rj-17) | 🟡 | Step-path, en-echelon joints | SSRM 1.213 vs UDEC 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | — | 1.24 / 1.2 | No closed form exists. Three rock bridges decide the factor, and the two finite element programs, XSLOPE and RS2, both give factors below UDEC's, their differences from it 2.1 percentage points apart. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | 1.01 vs 1.01 (0.0%) | — | 1.01 / 1.0 | |
| [19](#rj-19) | 🟢 | Bi-planar step-path failure | SSRM 1.623 vs rigid-block limit equilibrium 1.5914 (+2.0%) | 1.5 vs 1.5914 (−5.7%) | UDEC 1.46 (+11.2%) | 1.5 / 1.41 | One block slides on the basal joint, held by friction and by the rock bridge in tension. The referee is that block's rigid-block equilibrium, with the tensile cap divided by the trial factor as in the vendor's model, and UDEC's and RS2's factors are both below it. |
| [20](#rj-20) | 🟡 | Hammah & Yacoub Voronoi slope | SSRM 2.541 vs UDEC 2.46 (+3.3%) | 2.21 vs 2.46 (−10.2%) | — | 2.21 / 2.37 | At the lower end of the bracket the section stops moving, while three contacts inside it keep opening and closing in a repeating cycle. With no net movement over the cycle, the trial counts as standing. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | 8.27 vs 8.16 (+1.3%) | — | 8.27 / 8.5 | *not supported* — the vendor's model excavates a 2 m opening in its second stage and runs the strength reduction on the excavated state, which carries the stresses left by the first stage; staged excavation is outside the scope of a slope program. |
| 22 | <span class="nodata">⊘</span> | Joint model: hyperbolic softening | — | — | — | — | *not supported* — the problem exercises RS2's hyperbolic displacement- and work-softening joint law, which XSLOPE's interface element does not have; it reports no factor of safety. |
| 23 | <span class="nodata">⊘</span> | Joint model: residual strength and dilation | — | — | — | — | *no reference value* — a shear test on a single joint. It reports a stress-displacement curve, not a factor of safety, so there is nothing to score. The manual presents it as a comparison of residual strength and dilation, but the vendor's six models for it all have dilation switched off (`include_dilation: no`), so that comparison never exercised dilation. |

</div>

---

## The Rows

### 🟢 RJ-1a: Goodman & Bray block toppling, case a (rj001a) {#rj-1a}

Problem 1 is Goodman & Bray's toppling example, and all four of its cases share one section:
sixteen rock columns, 10 m wide and 4 to 40 m tall, standing on a base that rises at 30°. The
column sides are normal to the base, and the tops are cut off by a 56.6° slope face. The rock is
elastic (γ = 25 kN/m³, E = 20 GPa, ν = 0.3), so the model can only fail on its joints, which is
the assumption the closed-form solution makes. The joints are the column boundaries (see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint)). They have
no cohesion, a friction angle of 38.15°, which is the value this case is posed at, and the
corpus's standard stiffness pair, a normal stiffness k<sub>n</sub> = 10<sup>8</sup> kPa/m and a shear
stiffness k<sub>s</sub> = 10<sup>7</sup> kPa/m.

The referee is Goodman & Bray's iterative column analysis, recomputed on this section. It
reproduces the method's published mode pattern and requires a horizontal toe force of 0.36 kN/m
for limit equilibrium, against the 0.5 kN the manual states. Carried into a strength reduction, it gives each of the four cases a
referee within a rounding of the value the manual prints for it.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.99 vs 1.0000 (−1.0%) | 0.99 (+2.8%) | 0.99 / 0.97 |

<!-- test: file=files/rocscience/joints/rj001a.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-1a, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

**XSLOPE's 1.018 sits 1.8% above the closed form because of one assumption in Goodman & Bray's
method: where the thrust between two columns acts. They place each thrust at the top corner of
the contact. In the solved state the two columns lean together and the contact face stays closed
only over its upper quarter to half, so the resultant acts a tenth to a sixth of the face height
below the corner. With those heights read from the solved state, and nothing else changed,
Goodman & Bray's block-by-block limit equilibrium analysis, in which the force needed to hold
each block is carried down to the block below it, gives **1.0279**, the bottom of this row's bracket. The method's other three
assumptions are met exactly by the solution: every block balances in force and moment on the
interface stresses alone, every closed side contact is at its friction limit, and the base
reaction of every toppling block acts at its downslope corner.

The manual's figure shows the side boundaries as rollers, but in the vendor's model they are
fixed. XSLOPE uses the vendor's model, so its sides are fixed too. The 0.5 kN toe force, negligible against the lowest block's weight,
is not carried; the vendor's reruns of all four cases move the force from the toe at (−0.5, 0.866)
to the block corner at (−2.5, 4.330), so their "with improvement" factors are a different load case.

**Input file:** [rj001a.xlsx](files/rocscience/joints/rj001a.xlsx).

![RJ-1a: Goodman & Bray block toppling, case a (rj001a) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF (strength reduction factor), and the deformed section. The rock is elastic and carries no strain of its own, so the whole figure is the joints: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope, and the deformed section shows the stack rotating forward over the face column by column while the toe column slides out along its own stretch of base](images/RJ-1a.png)

### 🟢 RJ-1b: Goodman & Bray block toppling, case b (rj001b) {#rj-1b}

Case b uses case a's section and rock, with the joints at φ = 33.0239°, the lowest friction of
the four cases. The stack stands only because of the 2013 kN horizontal force the case is posed
with. That force is about twice the weight of the lowest column, and it is what lets a stack on 33°
joints stand where case a's needs 38°. It enters the model as a line load on the lloads sheet at
the upper-left corner of the lowest block, (−2.5, 4.330). It cannot be placed at the toe itself,
which is the end of the basal joint, where a load has no defined side to act on.

For limit equilibrium, the block-by-block analysis on this case needs a horizontal toe force of
2,012.86 kN/m, against the 2,013 kN the manual states: a difference of 0.007% on a stack weighing
83,500 kN/m.

The thrust heights differ from case a's. The toe force pushes the two lowest columns back into
their own step risers, which carries their thrusts far down the contact faces. With this case's
own measured heights, the block-by-block analysis gives **1.0078**, the bottom of this row's bracket.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.018** | 1.0000 (+1.8%) | 0.97 vs 1.0000 (−3.0%) | 0.99 (+2.8%) | 0.97 / 0.94 |

<!-- test: file=files/rocscience/joints/rj001b.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-1b, tier=gate, f_stand=1.0078125, f_fail=1.02734375, check=edges -->

The force's magnitude and direction match the vendor's model.

**Input file:** [rj001b.xlsx](files/rocscience/joints/rj001b.xlsx).

![RJ-1b: Goodman & Bray block toppling, case b (rj001b) — FEM inputs with the 2013 kN toe force drawn at the block corner it acts on, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's — slip up every column contact and along the stepped base, the stack rotating forward over the face — reached on joints four degrees weaker, which is what the force at the toe buys](images/RJ-1b.png)

### 🟢 RJ-1c: Goodman & Bray block toppling, case c (rj001c) {#rj-1c}

Case a's section and rock with the joints half a degree steeper, φ = 38.6598°, and the same 0.5 kN
toe force that does nothing. Half a degree is worth about two points of factor of safety here, which
is what the closed form's own pair of answers for cases a and c says as well — 1.0 against 1.02.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.0185 (+1.8%) | 1.01 vs 1.0185 (−0.8%) | 1.01 (+2.7%) | 1.01 / 0.99 |

<!-- test: file=files/rocscience/joints/rj001c.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-1c, tier=gate, f_stand=1.02734375, f_fail=1.046875, check=edges -->

Re-running Goodman & Bray's analysis with the thrust heights taken from XSLOPE's solved state at
the 10 m block width, instead of the top-corner assumption, gives **1.0469**, the lower edge of
that mesh's bracket. This case has case a's mechanism, half a degree of friction further on, and
the same single assumption accounts for the difference.

**Input file:** [rj001c.xlsx](files/rocscience/joints/rj001c.xlsx).

![RJ-1c: Goodman & Bray block toppling, case c (rj001c) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's at half a degree more friction: the same forward rotation of the column stack over the face, and the same slip on every column contact and along the stepped base](images/RJ-1c.png)

### 🟢 RJ-1d: Goodman & Bray block toppling, case d (rj001d) {#rj-1d}

Case d combines case c's joint friction angle with case b's 2013 kN horizontal force at the toe,
so the same stack is held by both. It is the strongest of the four, based on the closed form and
UDEC solutions.

Re-running Goodman & Bray's analysis with the thrust heights taken from XSLOPE's solved state,
instead of the top-corner assumption, gives **1.2422**, the lower edge of XSLOPE's own bracket.
Across all four cases XSLOPE stands above the closed form and RS2 below it, with UDEC between
them, and on each of the four the difference between XSLOPE and the closed form comes entirely
from where the thrust acts.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **1.252** | 1.2308 (+1.7%) | 1.19 vs 1.2308 (−3.3%) | 1.22 (+2.6%) | 1.19 / 1.16 |

<!-- test: file=files/rocscience/joints/rj001d.xlsx, type=fem_ssrm, expected_fs=1.252, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-1d, tier=gate, f_stand=1.2421875, f_fail=1.26171875, check=edges -->

The force reaches the same point case b's does.

**Input file:** [rj001d.xlsx](files/rocscience/joints/rj001d.xlsx).

![RJ-1d: Goodman & Bray block toppling, case d (rj001d) — FEM inputs with the 2013 kN toe force, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The right-hand panels are the state past the critical factor: slip runs up every column contact and along the stepped base, brightest on the columns at mid-slope and at the toe, and the deformed section shows the stack rotating forward over the face column by column — case b's mechanism at the strength this case's toe force holds](images/RJ-1d.png)

### 🟢 RJ-2: Alejano & Alonso block toppling (rj002) {#rj-2}

Problem 2 is Alejano & Alonso's block-toppling example. The section has a 9.85 m slope face at
58.65° and two joint sets, each a family of parallel joints. A basal joint rises at 30° from the toe
of the face and forms the stepped surface the columns stand on. Twenty-two column joints at 64°,
1.6 m apart, start from the same toe. Both sets have φ = 31° and no cohesion. The rock is elastic
(γ = 25 kN/m³, E = 20 GPa, ν = 0.3), as in the vendor's file, where its plasticity is set to none
(`Plasticity Specifications: Non`), so the model can only fail on its joints.

The referee is Goodman & Bray's block-by-block analysis, recomputed on the thirteen columns that the
two joint sets cut out of the toppling mass. It gives **0.7734**, against the 0.76 Alejano & Alonso publish for the
same method on the same problem.

| XSLOPE SSRM | Goodman & Bray referee | RS2 vs referee | UDEC | RS2 without / with improvement |
|---|---|---|---|---|
| **0.783** | 0.7734 (+1.2%) | 0.86 vs 0.7734 (+11.2%) | 0.87 (−10.0%) | 0.86 / 0.82 |

<!-- test: file=files/rocscience/joints/rj002.xlsx, type=fem_ssrm, expected_fs=0.783, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-2, f_stand=0.7734375, f_fail=0.79296875, check=edges, tier=gate -->

**RS2 and UDEC give nearly the same factor, and both are above the closed form.** Alejano & Alonso publish
their own UDEC run at 0.87 beside their Goodman & Bray 0.76, and RS2's two factors are 0.86 and
0.82. All three programs give factors above the recomputed closed form: XSLOPE by 1.2% and the
other two by about a tenth. RS2's model is available from the vendor, and this row's model
is transcribed from it. UDEC's result cannot be examined further on this problem, because the
manual says nothing about the settings of its UDEC model: no block rounding, no deformability and
no stiffness. The only statement of that kind anywhere in the manual is its
note on rigid blocks for problems 9 to 14.

**Input file:** [rj002.xlsx](files/rocscience/joints/rj002.xlsx).

![RJ-2: Alejano & Alonso block toppling (rj002) — FEM inputs, mesh, joint slip at the critical SRF and the section deformed 17×. The rock carries no strain of its own because it cannot yield, so every movement in the section is on a joint: slip gathers where the basal joint reaches the toe of the face, the columns standing on it open along their upper halves, and the deformed section shows them rotating out over the face while the rock below the basal joint stays put](images/RJ-2.png)

### 🟢 RJ-3: Lorig & Varona forward block toppling (rj003) {#rj-3}

Problems 3 to 7 are the toppling and plane-failure examples of Lorig & Varona (2004), and share
one section: a 260 m high slope at 55°. Problem 3 cuts it with two joint sets, both passing through
the origin: columns at 70°, 20 m apart, and a cross set at −20°, 30 m apart. The manual gives the
pair as "70 and 160" degrees, which describes the same two planes with the second angle measured
from the other end of the half circle.

The model this row is scored on is the vendor's file with one change. The rock is elastic, as the
vendor's file has it (plasticity set to none; γ = 26.0946 kN/m³, E = 9072 MPa, ν = 0.26), so only
the joints can fail. The joints have φ = 40° and no cohesion, as the source chapter lists them; the
vendor's file gives them c = 100 kPa, and that is the one change. The UDEC result the row is scored
against comes from the chapter, whose joints have friction only. RS2's two factors were computed
with the vendor's 100 kPa, so they are shown beside the row but are not like for like.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.115** | 1.13 (−1.3%) | 1.12 vs 1.13 (−0.9%) | 1.12 / 1.09 |

<!-- test: file=files/rocscience/joints/rj003.xlsx, type=fem_ssrm, expected_fs=1.115, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-3, f_stand=1.10546875, f_fail=1.125, check=edges, tier=gate -->

Two other readings of the same problem were run for comparison:

- On the vendor's file as given, with the joints at c = 100 kPa, the factor of safety is 1.232
  in XSLOPE, against UDEC's 1.13 (+9.0%).
- On the chapter's full printed list (joints friction only, and a rock that can yield, c = 675 kPa,
  φ = 43°), XSLOPE gives no settled value: the slope stands at 1.027 and fails at 1.125, and the
  trial between them at 1.047 stopped moving without balancing its forces, so the result is
  *unconfirmed*, at least 1.027 against UDEC's 1.13.

**Input file:** [rj003.xlsx](files/rocscience/joints/rj003.xlsx).

![RJ-3: Lorig & Varona forward block toppling (rj003) — FEM inputs, mesh, joint slip and the deformed section. The rock is elastic and carries no strain of its own, so the whole mechanism is on the two sets: slip runs along one −20° cross joint from the toe back under the crest, and the steep 70° joints above it slip, most strongly near the face, so the columns resting on that cross joint lean forward over the face. The two right-hand panels show the last trial at which the slope stands, F = 1.105, which is the state just below the factor of safety rather than past it, and the deformed section is drawn at 52 times scale. The run past the factor of safety that these panels usually show is not drawn: it was stopped in its first iteration, before the section had moved at all](images/RJ-3.png)

### 🟡 RJ-4: Lorig & Varona flexural toppling (rj004) {#rj-4}

Problem 4 uses the same section, cut by one joint set: columns at 70°, 20 m apart. This is
problem 3's first set without its cross joints, so the columns bend rather than topple as blocks.

The model this row is scored on is the vendor's file with one change. The rock can yield, as the
vendor's file has it: Mohr-Coulomb with a tensile cutoff of zero, which lets a column break in
bending (γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26, c = 675 kPa, φ = 43°). The joints have φ = 40° and
no cohesion, as the source chapter lists them; the vendor's file gives them c = 100 kPa, and that is
the one change. The UDEC result the row is scored against comes from the chapter, whose joints have
friction only. RS2's two factors were computed with the vendor's 100 kPa, so they are shown beside
the row but are not like for like.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.232** | 1.3 (−5.2%) | 1.19 vs 1.3 (−8.5%) | 1.19 / 1.27 |

<!-- test: file=files/rocscience/joints/rj004.xlsx, type=fem_ssrm, expected_fs=1.232, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-4, f_stand=1.22265625, f_fail=1.2421875, check=edges, tier=gate -->

Two other readings of the same problem were run for comparison:

- On the vendor's file as given, with the joints at c = 100 kPa, the factor of safety is 1.311
  in XSLOPE, against UDEC's 1.3 (+0.8%).
- On the chapter's full printed list (joints friction only, and a rock that can yield, c = 675 kPa,
  φ = 43°), the model is the same as this row's, because the vendor's rock already has the chapter's
  properties; the factor of safety is the row's own 1.232 in XSLOPE, against UDEC's 1.3 (−5.2%).

<!-- test: file=files/rocscience/joints/rj004_chapter.xlsx, type=fem_ssrm, expected_fs=1.232, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-4c, f_stand=1.22265625, f_fail=1.2421875, check=edges, tier=gate -->

**Input file:** [rj004.xlsx](files/rocscience/joints/rj004.xlsx).

![RJ-4: Lorig & Varona flexural toppling (rj004) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Here the rock can yield, and it does: a band of shear strain climbs from the toe across the columns to behind the crest, the joints slip through that band, and the columns are bent through it rather than rotated about it, which is what separates flexural toppling from the block toppling of problem 3](images/RJ-4.png)

### 🟡 RJ-5: Lorig & Varona backward block toppling (rj005) {#rj-5}

Problem 5 uses the same section, cut by two joint sets. The first is at −55°, 10 m apart, passing
through the toe and dipping out of the face, so the blocks lean back into the slope rather than
forward. The second is horizontal, 40 m apart.

The model this row is scored on is the vendor's file with one change. The rock is elastic, as the
vendor's file has it (plasticity set to none; γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26), so only the
joints can fail. The joints have φ = 40° and no cohesion, as the source chapter lists them; the
vendor's file gives them c = 100 kPa, and that is the one change. The UDEC result the row is scored
against comes from the chapter, whose joints have friction only. RS2's two factors were computed
with the vendor's 100 kPa, so they are shown beside the row but are not like for like.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.799** | 1.7 (+5.8%) | 1.65 vs 1.7 (−2.9%) | 1.65 / 1.86 |

<!-- test: file=files/rocscience/joints/rj005.xlsx, type=fem_ssrm, expected_fs=1.799, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-5, f_stand=1.7890625, f_fail=1.80859375, check=edges, tier=gate -->

Two other readings of the same problem were run for comparison:

- On the vendor's file as given, with the joints at c = 100 kPa, the factor of safety is 1.857
  in XSLOPE, against UDEC's 1.7 (+9.2%).
- On the chapter's full printed list (joints friction only, and a rock that can yield, c = 675 kPa,
  φ = 43°), the factor of safety is 1.096 in XSLOPE, against UDEC's 1.7 (−35.5%).

<!-- test: file=files/rocscience/joints/rj005_chapter.xlsx, type=fem_ssrm, expected_fs=1.096, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-5c, f_stand=1.0859375, f_fail=1.10546875, check=edges, tier=gate -->

With the chapter's yielding rock the mechanism changes: the slabs slide down the joints parallel to
the face and the rock at the toe yields beneath them, instead of the blocks toppling back into the
slope.

**Input file:** [rj005.xlsx](files/rocscience/joints/rj005.xlsx).

![RJ-5: Lorig & Varona backward block toppling (rj005) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own; slip runs along the −55° joints in the wedge behind the face while the horizontal bedding opens, and the deformed section shows the slabs stepping out over one another down the face, each leaning back into the slope as it goes](images/RJ-5.png)

### 🟢 RJ-6: Plane failure with daylighting discontinuities (rj006) {#rj-6}

Problem 6 uses the same section, cut by one joint set at −35°, 10 m apart, passing through the
origin. The joints dip out of the 55° face at a shallower angle than the face itself, so every one
of them daylights, that is, emerges on the face, and the slabs between them are free to slide out.
The rock is Mohr-Coulomb (γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26, c = 675 kPa, φ = 43°, no tensile
capacity).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.271** | 1.27 (+0.1%) | 1.25 vs 1.27 (−1.6%) | 1.25 / 1.31 |

<!-- test: file=files/rocscience/joints/rj006.xlsx, type=fem_ssrm, expected_fs=1.271, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-6, f_stand=1.26171875, f_fail=1.28125, check=edges, tier=gate -->

XSLOPE's factor matches UDEC's and lies between RS2's two factors, 1.25 and 1.31.

**Input file:** [rj006.xlsx](files/rocscience/joints/rj006.xlsx).

![RJ-6: plane failure with daylighting discontinuities (rj006) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Slip runs the full length of every joint that reaches the face, over a wedge bounded below by the joint through the toe, and the rock between them carries only a faint strain: the slabs slide out along the joints rather than breaking through anything](images/RJ-6.png)

### 🟢 RJ-7: Plane failure with non-daylighting discontinuities (rj007) {#rj-7}

Problem 7 uses the same section as problem 6, cut by one joint set at −70°, 20 m apart, passing
through the origin. The joints now dip out of the face more steeply than the 55° face itself, so
none of them daylights. A slab cannot slide out along one without shearing through rock, and the
factor of safety is higher than problem 6's.

The model this row is scored on is the vendor's file with one change. The rock can yield, as the
vendor's file has it, and is problem 6's: Mohr-Coulomb with no tensile capacity (γ = 26.1 kN/m³,
E = 9072 MPa, ν = 0.26, c = 675 kPa, φ = 43°). The joints have φ = 40° and no cohesion, as the
source chapter lists them, and its text calls these planes "cohesionless"; the vendor's file gives
them c = 100 kPa, and that is the one change. The UDEC result the row is scored against comes from
the chapter, whose joints have friction only. RS2's two factors were computed with the vendor's
100 kPa, so they are shown beside the row but are not like for like.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.486** | 1.5 (−0.9%) | 1.57 vs 1.5 (+4.7%) | 1.57 / 1.59 |

<!-- test: file=files/rocscience/joints/rj007.xlsx, type=fem_ssrm, expected_fs=1.486, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-7, f_stand=1.4765625, f_fail=1.49609375, check=edges, tier=gate -->

Two other readings of the same problem were run for comparison:

- On the vendor's file as given, with the joints at c = 100 kPa, the factor of safety is 1.564
  in XSLOPE, against UDEC's 1.5 (+4.3%).
- On the chapter's full printed list (joints friction only, and a rock that can yield, c = 675 kPa,
  φ = 43°), the model is the same as this row's, because the vendor's rock already has the chapter's
  properties; the factor of safety is the row's own 1.486 in XSLOPE, against UDEC's 1.5 (−0.9%).

<!-- test: file=files/rocscience/joints/rj007_chapter.xlsx, type=fem_ssrm, expected_fs=1.486, element_type=tri6, target_size=12.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-7c, f_stand=1.4765625, f_fail=1.49609375, check=edges, tier=gate -->

The manual's table for this problem prints the slope angle as 5°; its figure and the vendor's
model both have the same 55° slope as problem 6.

**Input file:** [rj007.xlsx](files/rocscience/joints/rj007.xlsx).

![RJ-7: plane failure with non-daylighting discontinuities (rj007) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. No joint reaches the face at a shallower angle than the face itself, so the failure cannot slide out along one: a strong band of shear strain cuts across the steep joints from the toe, the joints above it slip, and the mass above it moves out over rock it has had to break](images/RJ-7.png)

### 🟢 RJ-8: Flexural toppling in a base friction model (rj008) {#rj-8}

Problem 8 is Pritchard & Savigny's base-friction table model, scaled up a hundred times. It has a
30.5 m face at 78° and twelve columns at −60°, 5.08 m apart, standing on a horizontal basal joint
and closed at the back by a vertical joint. The three highest columns end on that back joint, as in
the vendor's model. The vendor's
file holds the column joints only inside the block that the basal and back joints bound. Generated
across the whole section instead, the same set would be sixteen column joints. Three of them
would lie entirely in the rock in front of the toe, below the level of the basal joint, where the
model has no columns, and three more would run on past the back joint at x = 68.4 into the strip
the vendor leaves uncut. The rock is Mohr-Coulomb (γ = 25.506 kN/m³,
E = 22.771 GPa, ν = 0.139, c = 60 kPa, φ = 39°). The joints have no cohesion, φ = 39°, and a
normal stiffness below the standard pair's: k<sub>n</sub> = 1.5 × 10<sup>7</sup> kPa/m, against
the usual 10<sup>8</sup>. Only problems 15 and 16 have softer joints.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **0.764** | 0.76 (+0.5%) | 0.75 vs 0.76 (−1.3%) | 0.75 / 0.75 |

<!-- test: file=files/rocscience/joints/rj008.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-8, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

The trials on this row come to rest faster than on any other row on this page. By comparison, the
[geotextile wall rows](rs2.md#rs2-48) of the RS2 corpus, which use the same joint element, reach
their iteration limit on five trials across their eight rows.

The rock's 75 kPa tensile strength is above the tensile strength that its own c and φ already
allow, at the apex of its Mohr-Coulomb envelope, so the tensile cap never controls.

**Input file:** [rj008.xlsx](files/rocscience/joints/rj008.xlsx).

![RJ-8: flexural toppling in a base friction model (rj008) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The strain gathers into one lobe per column along a band that climbs from the toe across the stack; above that band the columns are visibly bent rather than merely tilted, and the rock below it is unstrained. That is the break surface of flexural toppling rather than sliding along any one joint](images/RJ-8.png)

### 🟢 RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) {#rj-9}

Problem 9 is example 1a of Alejano et al.: a 50 m slope at 50°, cut by bedding joints dipping
**out of the face** at −50°, 3 m apart, with φ = 30°. A two-segment release trace at the toe, a
short joint with φ = 40°, undercuts the lowest slab. The model therefore has two joint strengths,
and the manual's geometry table lists them in the reverse order from its RS2 legend. On all six of
problems 9 to 14 the bedding and the short release traces at the toe have different friction
angles, so a single quoted joint friction angle for these problems is incomplete. The rock is
elastic with E = 2 × 10⁸ MPa, γ = 25 kN/m³ and ν = 0.3; the very high stiffness is the manual's
way of making the rock act as rigid blocks.

The face and the bedding dip at the same angle, so no slab can slide out along a single plane. The
mechanism is the bilinear one the problem is named for: sliding on a basal plane combined with
sliding along the release trace that the face undercuts. The release trace ends exactly on the
bedding plane it is meant to meet; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 1.01 vs 1.03 (−1.9%) | 0.40–1.45 | 1.01 / 1.09 |

<!-- test: file=files/rocscience/joints/rj009.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-9, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

RS2's two factors, 1.01 without the `Improve Joint Convergence` option and 1.09 with it, lie on
either side of XSLOPE's factor, as they do on six other problems in this corpus: 3, 4, 5, 6, 10
and 17.

Alejano's paper gives its limit equilibrium for this problem as a range rather than a single
value: 0.40 to 1.45, as the manual prints it. The range is shown beside the referee and does not
score the row, as with every limit-equilibrium value printed by a source on this page. The manual
verifies RS2 against the UDEC run, which is this row's referee.

**Input file:** [rj009.xlsx](files/rocscience/joints/rj009.xlsx).

![RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: the release trace near the toe carries the brightest slip, the bedding plane that starts at the crest and passes below the toe carries the rest, and the deformed section shows the slab between them sliding out over the bench](images/RJ-9.png)

### 🟢 RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) {#rj-10}

Problem 10 is example 1b: example 1a with the release joint moved five meters up the face, which
the manual says is the only difference between the two problems. The section, the bedding, both
friction angles and the rock are [problem 9](#rj-9)'s.

Moving the release upslope leaves the toe undercut as before and cuts a second block out of the
face above it. The release trace that ran from the face to the toe in example 1a stays, and a new
one leaves the face five meters higher. Both end exactly on the bedding plane that releases them;
see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | UDEC referee | RS2 vs referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.037** | 1.03 (+0.7%) | 0.92 vs 1.03 (−10.7%) | 0.43–1.45 | 0.92 / 1.08 |

<!-- test: file=files/rocscience/joints/rj010.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-10, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

RS2's two factors, 0.92 without the `Improve Joint Convergence` option and 1.08 with it, lie on
either side of XSLOPE's 1.037, as they do on problems [3](#rj-3), [4](#rj-4), [5](#rj-5),
[6](#rj-6), [9](#rj-9) and [17](#rj-17) and on no other problem here. Alejano's limit equilibrium is again a range
rather than a single value, 0.43 to 1.45, shown beside the referee and not scoring the row.

**Input file:** [rj010.xlsx](files/rocscience/joints/rj010.xlsx).

![RJ-10: Alejano et al. bilinear slab failure, example 1b (rj010) — FEM inputs, mesh, joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and, as in example 1a, almost none of it slips: the release trace high on the face carries the brightest slip, one bedding plane below the toe carries the rest, and the deformed section shows the block between them moving out over the bench](images/RJ-10.png)

### 🟢 RJ-11: Alejano et al. plowing sliding slab failure (rj011) {#rj-11}

Problem 11 is Alejano et al.'s plowing sliding slab: a 25 m slope at 50°, with bedding dipping out
of the face at −50°, 1.5 m apart, with φ = 30°, and two release traces with φ = 20°. One release
trace is at the toe and runs below the bench; the other runs from the face down onto the bedding
plane that releases the toe block. Plowing failure is the paper's name for sliding on a primary
discontinuity combined with sliding on a joint nearly parallel to the face, which lifts the toe
block and eventually rotates it out of the slope. As on problems 9 and 10, the rock is elastic and
stiff enough to act as rigid blocks: E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The release traces
end exactly on the bedding plane they are meant to meet; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

| XSLOPE SSRM | Rigid-block bound referee | RS2 vs referee | UDEC | Alejano Eq. (7) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.213** | 1.2148 (−0.1%) | 1.22 vs 1.2148 (+0.4%) | 1.21 (+0.2%) | 1.7582 · the paper prints 1.75 | 1.22 / 1.3 |

<!-- test: file=files/rocscience/joints/rj011.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-11, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

**On this problem Alejano's Eq. (7) gives a factor that no rigid-block mechanism can reach.** The
mechanism cuts out two blocks, the slab and the toe block below it. Above a reduction factor of
1.2148, no set of contact forces that stays within the friction limit on every joint can hold
those two blocks in place, however the forces are distributed. Eq. (7) returns 1.7582, about half
a unit of factor of safety above that, so it cannot score this row, and the rigid-block bound is
the referee instead. XSLOPE gives 1.213, RS2 1.22 and the paper's own UDEC run 1.21: the three
programs are inside 0.8% of one another, and all three are on the bound. [Problem 15](#rj-15) has
the same pattern, a closed form well above the programs, on a problem where no rigid-block bound is
available.

**Input file:** [rj011.xlsx](files/rocscience/joints/rj011.xlsx).

![RJ-11: Alejano et al. plowing sliding slab failure (rj011) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. One bedding plane from the crest to the toe carries almost all of the slip, and at its foot the two release traces cut out a small wedge: the deformed section shows that wedge lifted and rotated out over the bench while the slab above it slides down the plane, which is the plowing mechanism the paper names. The wedge is driven out and up by the slab above it, and at the panel's exaggeration a movement of centimeters draws as meters, so the block appears to leave the slope](images/RJ-11.png)

### 🟡 RJ-12: Alejano et al. plowing toppling slab failure (rj012) {#rj-12}

Problem 12 is the plowing toppling slab: a 25 m slope at 60°, cut by bedding dipping out of the
face at −60°, 1.5 m apart, with φ = 30°, and two short release traces with φ = 40°. One is at the
toe and runs below the bench; the other runs from the face down onto the bedding plane that
releases the toe block. The mechanism is problem 11's, but at 60° the rotation of the toe block,
rather than its sliding, governs. As on problems 9 to 11, the rock is elastic and stiff enough to
act as rigid blocks: E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3.

Alejano's Eq. (7), a moment balance about the toe for exactly this mechanism, is the referee.
Recomputed on the inputs in the vendor's file, it gives **1.9659**, against the 2.00 the paper
prints for the same example. The same calculation reproduces the paper's own table for its other
worked examples.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **2.033** | 1.9659 (+3.4%) | 1.39 vs 1.9659 (−29.3%) | 2.0324 | 1.78 (+14.2%) · Alejano prints 2.00 | 1.39 / 1.75 |

<!-- test: file=files/rocscience/joints/rj012.xlsx, type=fem_ssrm, expected_fs=2.033, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-12, f_stand=2.0234375, f_fail=2.04296875, check=edges, tier=gate -->

The release traces end exactly on the bedding plane they are meant to meet (see
[where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint)).

**The three independent answers for this problem disagree with one another, and XSLOPE's is the
one nearest the closed form.** Alejano's limit equilibrium, recomputed, gives 1.9659, the paper's
UDEC run 1.78 and RS2 1.39. Every source agrees on the mode that governs: the toe block rotates
out rather than sliding.

**Input file:** [rj012.xlsx](files/rocscience/joints/rj012.xlsx).

![RJ-12: Alejano et al. plowing toppling slab failure (rj012) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section, and only a few of its traces carry any slip: one release trace under the crest and the bedding beneath the toe block, which is the plowing pair. The two right-hand panels show the last trial at which the slope stands, the state just below the factor of safety rather than past it. On a rock this stiff the model moves only by microns until it fails, so the deformed section is drawn at tens of thousands of times scale. The run past the factor of safety that these panels usually show is not drawn: it was stopped in its first iteration, before the section had moved at all](images/RJ-12.png)

### 🟢 RJ-13: Alejano et al. plowing sliding slab, example 4 (rj013) {#rj-13}

Problem 13 is the paper's example 4, a plowing sliding slab: a 25 m slope at 55°, cut by bedding
dipping out of the face at −55°, 1.5 m apart, with φ = 25°, and two release traces with φ = 20°.
One is at the toe and runs below the bench; the other runs from the face down onto the bedding
plane that releases the toe block. The mechanism is problem 11's, with steeper and weaker bedding:
sliding on the primary discontinuity combined with sliding on a joint nearly parallel to the face,
which lifts the toe block. The rock is again elastic and stiff enough to act as rigid blocks:
E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The release traces end exactly on the bedding plane they
are meant to meet; see [where a joint ends on another](../fem/joints.md#where-a-joint-ends-on-another-joint).

Alejano's Eq. (7), recomputed on this problem's inputs, gives **1.0002**, with its sliding mode
governing; the paper prints 1.00 for the same example. The rigid-block bound for these two blocks
is 0.9988, so Eq. (7) sits on the bound, and the closed form and the rigid-block bound agree to
three figures here.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **0.998** | 1.0002 (−0.2%) | 1.0 vs 1.0002 (0.0%) | 0.9988 | 1.0 (−0.2%) | 1.0 / 1.05 |

<!-- test: file=files/rocscience/joints/rj013.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-13, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

**Input file:** [rj013.xlsx](files/rocscience/joints/rj013.xlsx).

![RJ-13: Alejano et al. plowing sliding slab, example 4 (rj013) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it carries slip: one bedding plane from the crest down to the toe does, and the two release traces at its foot cut out the toe block, which the deformed section shows lifted and pushed out over the bench while the slab above it slides down the plane](images/RJ-13.png)

### 🟢 RJ-14: Alejano et al. plowing sliding slab, example 5 (rj014) {#rj-14}

Problem 14 is the paper's example 5: example 4's section at 60°, with the two joint strengths
swapped. The bedding dips out of the face at −60°, 1.5 m apart, with φ = 20°, the weakest bedding
of the six problems 9 to 14, and the two release traces have φ = 30°. The rock is again elastic and
stiff enough to act as rigid blocks: E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3.

Alejano's Eq. (7), recomputed on the inputs in the vendor's file, gives **1.2034** for this
example, below the rigid-block bound for its two blocks. It is the referee.

| XSLOPE SSRM | Alejano Eq. (7) referee | RS2 vs referee | Rigid-block bound | UDEC | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.232** | 1.2034 (+2.4%) | 0.89 vs 1.2034 (−26.0%) | 1.2686 | 0.9 (+36.9%), 0.9994 with the paper's corner rounding corrected · the paper prints 1.00 | 0.89 / 1.09 |

<!-- test: file=files/rocscience/joints/rj014.xlsx, type=fem_ssrm, expected_fs=1.232, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-14, f_stand=1.22265625, f_fail=1.2421875, check=edges, tier=gate -->

**The 1.00 the paper prints for this example does not follow from the inputs it prints.** Every
input that Eq. (7) uses matches the vendor's model; the toe block's base, for example, matches to
three figures, 3.589 m against the 3.588 m printed. On those inputs the equation returns 1.2034. For
it to return 1.00, the toe block's base would have to be 6.08 m, which is a different problem. The
manual's other value, UDEC's 0.9, comes from the rounding of block corners in the UDEC model, by
the paper's own account: with the rounding radius reduced, the paper reports 0.9994 for the same
model.

**Input file:** [rj014.xlsx](files/rocscience/joints/rj014.xlsx).

![RJ-14: Alejano et al. plowing sliding slab, example 5 (rj014) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: one bedding plane from the crest to the toe carries the slip, with the release trace at the toe opening as the slab above it slides out over the bench](images/RJ-14.png)

### 🟢 RJ-15: Partially joint-controlled footwall slope (rj015) {#rj-15}

Problem 15 is a 40 m footwall slope at 40° whose bedding dips in the same direction at the same
angle, 2 m apart. The slabs therefore lie parallel to the face, and a failure has to break through
rock at the toe to get out. This is the one problem in the Alejano family whose rock can yield:
Mohr-Coulomb, c = 200 kPa, φ = 35°, γ = 28 kN/m³, E = 1 GPa, ν = 0.3. The joints share the corpus's
softest stiffness pair with problem 16, k<sub>n</sub> = 5 × 10<sup>6</sup> kPa/m and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m,
twenty times below the corpus's standard pair, with no cohesion and φ = 25°.

The manual's table prints the slope height as 25 m. The vendor's section rises 40 m, from the toe
at (0, 0) to the crest at (−47.6701, 40), which is the height the source paper states, and XSLOPE
builds the vendor's section.

Alejano's limit equilibrium for the footwall, Eqs. (9)–(10), recomputed on the vendor file's inputs
at the optimum the paper states (a break-out surface inclined 14° to the bedding and emerging at
55°), gives **1.7985**, against the value the paper prints, which is in the table below. Minimized
over its own two angles, it settles within half a degree of that optimum. The referee is the
limit-equilibrium search by Slide2, Rocscience's own slope program, at 1.25.

| XSLOPE SSRM | Slide2 LE search referee | RS2 vs referee | Alejano Eqs. (9)–(10) | UDEC-SSRT (Alejano) | RS2 without / with improvement |
|---|---|---|---|---|---|
| **1.271** | 1.25 (+1.7%) | 1.28 vs 1.25 (+2.4%) | 1.7985 (single-bed formula; the paper prints 1.72) | 1.6 (−20.6%) | 1.28 / 1.42 |

<!-- test: file=files/rocscience/joints/rj015.xlsx, type=fem_ssrm, expected_fs=1.271, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-15, f_stand=1.26171875, f_fail=1.28125, check=edges, tier=gate -->

**Alejano's equations evaluate one fixed mechanism, and every program free to search for its own
failure surface finds a weaker one.** In Eqs. (9)–(10) a wedge is driven out through a single slab,
and the break-out surface always crosses exactly one 2 m bed. Their factor rises with the bed
thickness, so the thinnest break-out the formula allows also gives the lowest factor it can report.
XSLOPE's strength reduction at 1.271, RS2's at 1.28 and the Slide2 search at 1.25 have no such
restriction, and they sit inside 2.4% of one another, far below the closed form. Both finite
element programs differ from the closed form by the same amount and in the same direction, so the
gap comes from what the formula is able to consider, not from either program. The searched
limit-equilibrium answer is therefore the referee, as the rigid-block bound is on
[problem 11](#rj-11).

The rock's 1000 kPa
tensile cap is above the apex of its Mohr-Coulomb envelope and never controls.

**Input file:** [rj015.xlsx](files/rocscience/joints/rj015.xlsx).

![RJ-15: partially joint-controlled footwall slope (rj015) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The bedding slips over a long stretch behind the face, and the rock's only strain is a small patch at the toe where the slab has to break through to get out — the coupled mechanism the source paper describes](images/RJ-15.png)

### 🔴 RJ-16: Barla et al. tilt-table block toppling (rj016) {#rj-16}

Problem 16 is a laboratory test rather than a slope. Fourteen columns of 9 cm blocks are stacked
into a 63.4° staircase on a plate, and the plate is tilted until the stack topples; the problem
reports the tilt angle at which that happens. The blocks are elastic (E = 350 MPa, ν = 0.2,
γ = 28 kN/m³) on an effectively rigid plate, and the joints have no cohesion, φ = 38°, and the
corpus's softest stiffness pair, the same as problem 15's. In the vendor's files the 0° model sits on rollers, while the
nine tilted models fix every exterior node in both directions and converge to a tolerance a
hundred times tighter. The plate carries no body force in any of the ten.

The vendor tilts the model itself, with one file per degree of tilt. XSLOPE keeps the model level
and turns the load instead. A tilt θ is represented by a horizontal seismic coefficient k = tan θ
acting toward the face, and k is raised, with the strengths unreduced (F = 1), until the stack no
longer stands. The mesh size is the block size, 0.09 m, following the corpus's rule of meshing at
the joint spacing.

| XSLOPE tilt | UDEC referee | RS2 vs referee | Experiment | Goodman & Bray, plate tilted | RS2 without / with improvement |
|---|---|---|---|---|---|
| **10.24°** | 11° (−6.9%) | 9° vs 11° (−18.2%) | 9° | 7.6° | 9° / 7° |

<!-- test: file=files/rocscience/joints/rj016.xlsx, type=fem_tilt, k_stand=0.179846, k_fail=0.181396, expected_tilt=10.24, tolerance=0.1, element_type=tri6, target_size=0.09, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-16 -->

The stack stands at k = 0.1798, which is 10.20°, and fails at k = 0.1814, which is 10.28°. It
stands at every coefficient tried below that and fails at every one above, through 0.25, and the
test for this row re-solves both coefficients at full strength. When the solver is left to iterate with
no rule deciding whether a trial stands or fails, the result is the same: the stack comes to rest
at every tilt up to 9.9° and moves at a steady rate from 10.2° on.

**A seismic coefficient is not exactly a rotation, but on this model the difference has no
measurable effect.** Tilting the model by θ turns the body force through θ and leaves its magnitude
at γ. A coefficient k = tan θ turns it through the same angle and also multiplies its magnitude by
1/cos θ, which is 1.6% more at 10.2°. With zero joint cohesion and no block able to yield, the
state depends on the direction of the body force and not on its size. Re-solving the standing and
failing coefficients with every unit weight scaled by 1.016 and then by 0.5 leaves both results
unchanged.

Goodman & Bray's column analysis, with the plate tilted, puts the toppling tilt of these rigid
columns at 7.6°. It allows no contact pressure below a column's corner, so it is a lower bound, and
every other value is above it: the physical stack and RS2 at 9°, XSLOPE at 10.2°, and UDEC, whose
block contacts roll and re-form as the columns lean, at 11°.

The vendor's plate has no weight, but XSLOPE's finite element model requires a positive unit
weight, so the plate is given the 27 kN/m³ stated in its row of the vendor's material properties.

**Input file:** [rj016.xlsx](files/rocscience/joints/rj016.xlsx).

![RJ-16: Barla et al. tilt-table block toppling (rj016) — FEM inputs with the seismic coefficient that stands for the tilt, mesh with the vendor's rollers, joint slip at the first coefficient at which the stack fails, and the deformed section. The slip is on the vertical joints between the lower columns and along the bedding one block above the plate, and the deformed section at 81x shows what that adds up to: every column leaning downslope about its own base, the tall ones at the back furthest over, which is toppling rather than the stack sliding along the plate](images/RJ-16.png)

### 🟡 RJ-17: Step-path failure, en-echelon joints (rj017) {#rj-17}

Problem 17 uses [problem 18](#rj-18)'s section and Mohr-Coulomb rock (c = 25 kPa, φ = 25°, no
tensile capacity), cut by three joints at 36.1° that stop short of one another instead of running
through, with problem 18's joint strength, c = 1 kPa and φ = 35°. The intact rock between the end
of one joint and the start of the next is a rock bridge. A step-path failure has to break through
the bridges, and they are why this slope stands where problem 18, with continuous joints, fails.

In the vendor's model the strength reduction applies only inside the SSR search area, a rectangle
drawn on the manual's figure but absent from its tables. Outside it RS2 holds every element linear
elastic: 1,173 of the model's 2,579 elements. XSLOPE enters the rectangle as a polygon overlay and
decides which elements fall inside it by the same test, the position of each element's centroid.
No closed form exists for a step path through rock bridges, so the referee is the manual's single
UDEC run.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **1.213** | 1.29 (−6.0%) | 1.24 vs 1.29 (−3.9%) | 1.24 / 1.2 |

<!-- test: file=files/rocscience/joints/rj017.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-17, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

**The rock bridges control this row.** The rock has no tensile capacity at all, so the bridges
between the joint segments carry no tension. The strength reduction lowers their 25 kPa cohesion
together with the joint friction, and the factor of safety is the reduction at which the three
short bridges shear through. The two finite element programs, XSLOPE and RS2, both give factors
below UDEC's, and their percentage differences from it are 2.1 percentage points apart.

In the vendor's mesh the edge of the rectangle follows element edges, so the search area is a
staircase rather than a rectangle. A second model, transcribed with that staircase instead of the
rectangle, holds 693 elements of XSLOPE's mesh elastic, where the rectangle holds 694. Every trial
of its search gives the same result as the rectangle's, so its bracket is identical.

**Input files:** [rj017.xlsx](files/rocscience/joints/rj017.xlsx), and the staircase variant
[rj017_staircase.xlsx](files/rocscience/joints/rj017_staircase.xlsx).

![RJ-17: step-path failure with en-echelon joints (rj017) — FEM inputs with the vendor's search area drawn as the held-elastic region around it, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. All three joints slip along their whole length, and the strain runs between them: a band climbs from the toe of the face through the rock bridge below the lowest joint and on past the upper two, which is the step-path the problem is named for — shear along the joints and broken rock between them](images/RJ-17.png)

### 🟢 RJ-18: Step-path failure, continuous joints (rj018) {#rj-18}

Problem 18 is a 45 × 20 m section of one Mohr-Coulomb rock (γ = 19.62 kN/m³, E = 20 GPa, ν = 0.3,
c = 25 kPa, φ = 25°, no tensile capacity) with a slope face rising from (17, 8.2) to (26.9, 20).
Three parallel joints at 36.1° run from the face to the crest, 0.883 m apart measured perpendicular
to the joints. The joints have c = 1 kPa, φ = 35°, k<sub>n</sub> = 10<sup>8</sup> kPa/m and
k<sub>s</sub> = 10<sup>7</sup> kPa/m, and their strength is reduced together with the rock's in the
strength reduction.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **0.998** | 1.01 (−1.2%) | 1.01 vs 1.01 (0.0%) | 1.01 / 1.00 |

<!-- test: file=files/rocscience/joints/rj018.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-18, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

**Input file:** [rj018.xlsx](files/rocscience/joints/rj018.xlsx).

![RJ-18: step-path failure through three continuous joints (rj018) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Every joint is slipping along its lower half and open along its upper, and the three slabs between them slide out together down the 36.1° path; the rock itself carries almost no plastic strain](images/RJ-18.png)

### 🟢 RJ-19: Bi-planar step-path failure (rj019) {#rj-19}

Problem 19 is one Mohr-Coulomb rock (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 10,500 kPa, φ = 35°,
tensile capacity 200 kPa) cut by two joints that do not meet, with a rock bridge between them: a
basal joint at 28.4° and an upper joint at 56.3°. Both have c = 0, φ = 40° and the same stiffness
pair as RJ-18.

**A tensile cap, rather than a cohesion, controls this row.** The two joints leave a rock
bridge 1.414 m long between their tips. At the rock's 10,500 kPa cohesion that bridge can carry
14,849 kN/m in shear. That is 1.2 times the weight of the whole sliding block, 11,929 kN/m, and
2.6 times the component of that weight down the basal joint, so the bridge cannot be sheared
through. What lets the block move is the same rock's 200 kPa tensile capacity, which gives
283 kN/m across the bridge. On this row that capacity is divided by the trial factor, as in every
vendor model.

The referee is the rigid-block limit equilibrium of that mechanism, recomputed from the inputs in
the workbook: the block slides on the basal joint, opens the upper joint behind it, and pulls the
bridge apart in tension, with the joint's friction and the rock's tensile capacity both divided by
the trial factor. The block holds up to a factor of 1.5914; with friction alone it holds up to
1.551. XSLOPE's 1.623 is 2.0% above that limit. Both values in the vendor's manual, UDEC's and
RS2's, are below the limit that the block's statics allow, so UDEC and RS2 let this block fail by
some mechanism that the rigid-block analysis does not contain.

| XSLOPE SSRM | Rigid-block limit equilibrium referee | RS2 vs referee | UDEC (manual) | RS2 without / with improvement |
|---|---|---|---|---|
| **1.623** | 1.5914 (+2.0%) | 1.5 vs 1.5914 (−5.7%) | 1.46 (+11.2%) | 1.5 / 1.41 |

<!-- test: file=files/rocscience/joints/rj019.xlsx, type=fem_ssrm, expected_fs=1.623, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=true, k0=1, benchmark=RJ-19, f_stand=1.61328125, f_fail=1.6328125, check=edges, tier=gate -->

The joint inclinations are taken from the vendor's model. The manual's table prints 59°, its figure
is dimensioned at 56° and 28°, and the joint endpoints in the model give 56.3° and 28.4°.

**Input file:** [rj019.xlsx](files/rocscience/joints/rj019.xlsx).

![RJ-19: bi-planar step-path failure with a rock bridge (rj019) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Both joints have opened, the basal one is slipping at its lower end, and the only strain in the rock is the patch at the bridge between the two joint tips, where the block above has to break through to move](images/RJ-19.png)

---

### 🟡 RJ-20: Hammah & Yacoub Voronoi slope (rj020) {#rj-20}

An 80 × 70 m section with a 60 m face at 71.6°, tessellated into Voronoi blocks over the whole of it.
The rock is Mohr-Coulomb (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 1,000 kPa, φ = 35°, no tensile
capacity) and every block wall is a joint at c = 500 kPa, φ = 20° with the corpus's standard
stiffness pair. The source is Hammah, Yacoub & Curran (2009), whose third example is a rock mass
of Voronoi blocks. The paper holds the block size fixed and follows how the failure changes as the
slope gets higher, from 10 m to 120 m.

The vendor's model departs from that paper in two ways. Its network is coarser: 525 blocks over
the 4,400 m² section, 0.119 blocks per square meter, where the paper's network has 0.2, so the
vendor's mean block is 8.4 m² against the paper's 5 m². And its rock friction is 35°, where the
paper's is 30°. The tessellation was generated in UDEC and imported; the model file holds it as
523 joint boundaries, 1,177 segments, and they are transcribed verbatim. Their mean block width,
2.895 m, is this row's mesh size, standing in for the joint spacing every other row meshes at.

The referee is UDEC's 2.46. That is the vendor's own UDEC run on the vendor's network, and the only
value computed on this network with these inputs. The paper ran only its own finite element
program, on its finer network with 30° rock. At this 60 m height it gives 1.5 with the joints open
where they reach the slope surfaces and 1.6 with them closed. Those two values belong to a
different problem and do not score this row.

| XSLOPE SSRM | UDEC referee | RS2 vs referee | RS2 without / with improvement |
|---|---|---|---|
| **2.541** | 2.46 (+3.3%) | 2.21 vs 2.46 (−10.2%) | 2.21 / 2.37 |

<!-- test: file=files/rocscience/joints/rj020.xlsx, type=fem_ssrm, expected_fs=2.541, element_type=tri6, target_size=2.895, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, max_iter_ceiling=1000000, tension_srf=false, k0=1, benchmark=RJ-20, f_stand=2.53125, f_fail=2.55078125, check=edges, tier=gate -->

The slope stands at the lower end of the bracket and fails at the upper end. At the lower end the
section stops moving, but three contacts inside the sliding mass keep opening and closing in a
cycle that repeats every eight iterations, so the forces never balance exactly; after a million
iterations with no net movement over the cycle, the trial counts as standing. At the upper end the
section is still moving at nearly the same rate after 225,000 iterations, so that trial fails. The
whole search takes nearly five hours. The figure shows the mass failing on a surface picked through the block walls rather than
along any one plane it contains. The paper reports the same of its own Voronoi example: every
failure mechanism it found had a curved overall shape.

XSLOPE's `xslope.joints.voronoi` generates a Voronoi joint network of this kind at any block size.
Networks it generates on the section as the manual states it, an 80 × 70 m block with the face cut
from it, can be meshed. This file's outline is different: it also has a vertex wherever one of the
vendor's joints meets the boundary. A generated joint that ends within a fraction of a millimeter
of one of those vertices leaves a sliver of boundary that the mesher cannot split along. A generated
network should therefore be placed on the plain section as the manual states it, not on this file's
outline. This row is scored on the vendor's own tessellation.

**Input file:** [rj020.xlsx](files/rocscience/joints/rj020.xlsx).

![RJ-20: Hammah & Yacoub Voronoi slope (rj020) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The slip runs from the toe up through the block walls on a curved path toward the crest, where the walls behind the crest open, taking whichever wall of each block lies nearest that line, and the rock carries almost no strain except a patch at the toe where the path turns: the mass fails on a surface picked out of the tessellation rather than along any one joint in it](images/RJ-20.png)
