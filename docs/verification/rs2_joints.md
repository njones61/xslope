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

**The referee.** Each problem is scored against the referee its own manual page names — the answer
RS2 is verified against there. That is UDEC on every problem but the first, whose four cases are
Goodman & Bray's Example 1 and whose stabilizing toe force is the force that method computes for
limit equilibrium, and which the manual scores against Goodman & Bray. Where a source paper
publishes a limit-equilibrium answer of its own, it is recorded in the column beside the referee,
labeled, and scores nothing; so is the Slide2 factor Rocscience reports beside problem 15, which is
the vendor's own companion program on its own model. The reason for that separation is the method:
XSLOPE's SSRM reduces strength on a continuum with interfaces in it, which is what RS2's SSR and
UDEC's strength-reduction technique do, and the [shared scoring](index.md#how-the-match-dots-are-scored)
derives a dot only from a same-method pairing. Problem 1 is the one place a limit equilibrium
scores here, because it is what the manual scores against there; its UDEC run is recorded beside
it. RS2's factors are recorded beside every row and score none of them, for the reason the next
paragraph gives.

**RS2's two factors.** The manual reports each problem twice, *without* and *with* the "joint
improvement" option. The vendor's rerun models show what that option is: `improve_joint_convergence`
turns on, the convergence criterion and the SRF convergence type change, `CoupledSSR` turns **off**
— so the joints are no longer reduced with the rock — and the iteration limit and tolerance go back
to their defaults. It is a different analysis, not a better-converged one, and on problem 4 it
raises the factor while on problems 17, 18 and 19 it lowers it. XSLOPE reduces the joint with the
rock, which is the *without* setting, and both of RS2's numbers are recorded so a reader can see
the spread.

**What answers a bracket edge.** {#what-answers-a-bracket-edge}
A bracket asks of each trial factor whether the model stands there, and a trial that runs out of
sweeps with nothing to say has not answered: a lock closed on one is a statement about the budget.
There is one exception, and three rows here turn on it. On a jointed model the solver offers the
viscoplastic loop's own state to a Newton corrector at its checkpoints, the last of them the budget
exit; where the corrector reaches equilibrium from that state inside its force, yield and
displacement gates, it records the certification on the trial. That is an independent driver saying
the model has a static equilibrium at that strength, which answers the question however many sweeps
preceded it, so such a trial is read as decided whatever its sweep count. A trial with no
certification is decided only inside its budget, as before, and a corrector REFUSAL is not a
verdict of any kind — it is the absence of one. `tools/ssrm_trial_audit.py` reports the reading per
row and `test/corrector_certified_check.py` holds the rule.

**The budget, and what it does not buy.** A joint reaches equilibrium by growing slip, so a jointed
model settles over tens of thousands of viscoplastic sweeps where a bonded one settles over
hundreds; a strength-reduction trial that runs out of sweeps is recorded undecided, the bracket
reads that as not standing, and the factor comes out low. Every row here states its own sweep
budget and its trial record is checked (`tools/ssrm_trial_audit.py`), the same rule the
[geotextile wall family](rs2.md#rs2-48) runs under.

What a trial costs varies by two orders of magnitude across this corpus. The longest trial to reach
a verdict takes 3 745 sweeps on problem 8 and 196 201 on problem 18, both at a budget of 250 000.
Six of the ten rows measured carry at least one trial that does not decide inside that budget, and
on problem 5 four trials do not, including both edges of its bracket.

That is not, on the evidence, simply a matter of allowing more sweeps. Problem 7's undecided trial
was re-solved on its own at **four times** the budget — a million sweeps against 250 000 — and came
back with the same verdict it had before, `STABLE_STUCK`: the model neither converges nor diverges
there, and the extra sweeps changed nothing about that. An undecided trial of this kind is a
statement about the convergence criterion, not only about the iteration limit, and a row held back
by one is not waiting on machine time alone.

`benchmarks/rocscience/build_joint_problems.py` writes the input files, which are named `rjNNN`
by the manual's own problem number — `rj018.xlsx` is problem 18 — with a letter suffix where the
manual letters its cases. `make_rs2_joint_figures.py` writes the figures.

## Status

Status terms follow the [shared definitions](index.md#status-terms) and match dots the
[shared scoring](index.md#how-the-match-dots-are-scored). Of the manual's 23 problems, 21 report a
factor of safety or a tilt angle; problems 22 and 23 report neither, being shear-box tests of the
joint model whose output is a stress-displacement curve.

<div class="corpus-summary match" markdown>

| # | Match | Problem | Referee | Also published | RS2 without / with improvement | Notes |
|---:|:-:|---|---|---|---|---|
| [1a](#rj-1a) | 🟡 | Goodman & Bray block toppling, case 1a | SSRM 1.037 vs Goodman & Bray 1.0 (+3.7%) | UDEC 0.99 | 0.99 / 0.97 | |
| [1b](#rj-1b) | 🟢 | Goodman & Bray block toppling, case 1b | SSRM 1.018 vs Goodman & Bray 1.0 (+1.8%) | UDEC 0.99 | 0.97 / 0.94 | |
| [1c](#rj-1c) | 🟡 | Goodman & Bray block toppling, case 1c | SSRM 1.057 vs Goodman & Bray 1.02 (+3.6%) | UDEC 1.01 | 1.01 / 0.99 | |
| [1d](#rj-1d) | 🟢 | Goodman & Bray block toppling, case 1d | SSRM 1.252 vs Goodman & Bray 1.23 (+1.8%) | UDEC 1.22 | 1.19 / 1.16 | |
| [2](#rj-2) | 🔴 | Alejano & Alonso block toppling | SSRM 0.764 vs UDEC 0.87 (−12.2%) | LE (Goodman & Bray) 0.76 | 0.86 / 0.82 | |
| [3](#rj-3) | <span class="nodata">⊘</span> | Lorig & Varona forward block toppling | UDEC 1.13 | — | 1.12 / 1.09 | *reported, no lock* — the upper edge of the bracket reaches the sweep budget without a verdict. |
| [4](#rj-4) | <span class="nodata">⊘</span> | Lorig & Varona flexural toppling | UDEC 1.3 | — | 1.19 / 1.27 | *reported, no lock* — one trial of the refinement step reaches the sweep budget without a verdict. |
| [5](#rj-5) | <span class="nodata">⊘</span> | Lorig & Varona backward block toppling | UDEC 1.7 | — | 1.65 / 1.86 | *reported, no lock* — four trials of the bracket, including both its edges, reach the sweep budget without a verdict. |
| [6](#rj-6) | <span class="nodata">⊘</span> | Plane failure, daylighting | UDEC 1.27 | — | 1.25 / 1.31 | *reported, no lock* — both edges of the bracket reach the sweep budget without a verdict. |
| [7](#rj-7) | <span class="nodata">⊘</span> | Plane failure, non-daylighting | UDEC 1.5 | — | 1.57 / 1.59 | *reported, no lock* — one trial of the refinement step reaches the sweep budget without a verdict. |
| [8](#rj-8) | 🟢 | Flexural toppling, base friction model | SSRM 0.764 vs UDEC 0.76 (+0.5%) | — | 0.75 / 0.75 | |
| [9](#rj-9) | 🟢 | Bilinear slab failure, example 1a | SSRM 1.018 vs UDEC 1.03 (−1.2%) | LE (Alejano) 0.40–1.45 | 1.01 / 1.09 | |
| 10 | <span class="nodata">⊘</span> | Bilinear slab failure, example 1b | UDEC 1.03 | LE (Alejano) 0.43–1.45 | 0.92 / 1.08 | *planned* — the corpus bracket is run and decided on every trial; the refinement step a lock needs is not yet measured. |
| [11](#rj-11) | 🟢 | Ploughing sliding slab failure | SSRM 1.213 vs UDEC 1.21 (+0.2%) | LE (Alejano) 1.75 | 1.22 / 1.3 | |
| [12](#rj-12) | <span class="nodata">⊘</span> | Ploughing toppling slab failure | UDEC 1.78 | LE (Alejano) 2.0 | 1.39 / 1.75 | *reported, no lock* — a step of refinement moves the factor, so it is not the section's answer but its mesh's. |
| [13](#rj-13) | 🟢 | Ploughing sliding slab, example 4 | SSRM 0.998 vs UDEC 1.0 (−0.2%) | LE (Alejano) 1.0 | 1.0 / 1.05 | |
| [14](#rj-14) | <span class="nodata">⊘</span> | Ploughing sliding slab, example 5 | UDEC 0.9 | LE (Alejano) 1.0 | 0.89 / 1.09 | *reported, no lock* — a step of refinement moves the factor, as it does on problem 12. |
| [15](#rj-15) | <span class="nodata">⊘</span> | Partially joint-controlled footwall | UDEC 1.6 | LE (Alejano) 1.72 · Slide2 (vendor) 1.25 | 1.28 / 1.42 | *reported, no lock* — the upper edge of the bracket reaches the sweep budget without a verdict. |
| 16 | <span class="nodata">⊘</span> | Barla et al. tilt-table block toppling | UDEC 11° | experiment 9° | 9° / 7° | *blocked* — the problem scores the TILT ANGLE at which a block grid topples, found by rotating gravity through a staged sweep; XSLOPE's seismic coefficient tilts the load but the row needs the sweep and a toppling criterion, neither of which is a strength reduction. |
| 17 | <span class="nodata">⊘</span> | Step-path, en-echelon joints | UDEC 1.29 | — | 1.24 / 1.2 | *blocked* — the vendor model's second, elastic material has no boundary of its own. Recovered from the element-material map it is an 83-vertex staircase of element edges, trending vertical near x = 12 and horizontal near y = 3 with excursions of about one element either side, so the boundary is a property of the vendor's mesh rather than of its model. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | — | 1.01 / 1.0 | |
| [19](#rj-19) | <span class="nodata">⊘</span> | Bi-planar step-path failure | UDEC 1.46 | — | 1.5 / 1.41 | *reported, no lock* — the upper edge of the bracket reaches the sweep budget without a verdict, at the corpus mesh and at a finer one. |
| 20 | <span class="nodata">⊘</span> | Hammah & Yacoub Voronoi slope | UDEC 2.46 | — | 2.21 / 2.37 | *blocked* — the manual states no block size and no seed, and the vendor model carries the tessellation as 523 digitized traces rather than as a generated network, so the input is not reproducible from anything published. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | — | 8.27 / 8.5 | *blocked* — the second stage of the vendor model excavates a 2 m opening and the strength reduction runs on the excavated state, which carries the stress the first stage left behind; XSLOPE has no staged construction. |
| 22 | <span class="nodata">⊘</span> | Joint model: hyperbolic softening | — | — | — | *not supported* — the problem exercises RS2's hyperbolic displacement- and work-softening joint law, which XSLOPE's interface element does not have; it reports no factor of safety. |
| 23 | <span class="nodata">⊘</span> | Joint model: residual strength and dilation | — | — | — | *no lock possible* — the problem reports no factor of safety, and its six vendor models all carry `include_dilation: no`. See [The dilation problem](#the-dilation-problem). |

</div>

---

## Where a Joint Ends on Another Joint {#joint-terminations}

Ten of the manual's twenty-one scorable problems turn on one property of their geometry: a joint
that stops **on** another joint rather than crossing it.

A joint is an interface, so the mesh carries two elements on every edge of it, one on each face,
and the split gives every node on a jointed curve one copy per wedge of material around it — four
at a crossing, where the shared node sits at the middle of four wedges of rock, and three at a
termination. The wedge count is read from the mesh's own topology, so a termination needs nothing
of its own there. What it needs is for the two lines to meet at **one point**: the meeting point
is a vertex of the ending line and of nothing else, and the mesher only places a node where the
geometry carries one.

Whether they meet at one point is a question about arithmetic rather than about the drawing. Two
tips that are meant to lie on a line land beside it instead — Goodman & Bray's fifteen column
contacts by between 7 × 10⁻¹⁶ and 3 × 10⁻¹⁴, because they are computed from the same base angle
the line is, and Alejano's release traces by between 2 × 10⁻⁷ and 2 × 10⁻⁶, because the vendor
states their lower tip to six decimals — and the sliver between the two lines is thinner than any
mesh resolves. So **every jointed line's end within a millionth of the section of another jointed
line is moved onto it**, and the line it stops on is given a vertex at that same point, before the
mesher runs. The through line is not moved: it is the plane the ending line belongs to and it keeps
the geometry it was stated with. An end already on another line is left exactly where it is.

Problem 1's stepped base is fifteen such terminations in one section — each column's basal contact
begins partway along its downslope neighbour's side joint — and problems 9 to 14 are one or two
each, where a release trace runs from the crest down onto the bedding plane that releases it.

Problem 8's thirteen columns stand on their basal joint the same way, and there the tips land on
it exactly; every other junction in the corpus is a crossing.

## The Dilation Problem {#the-dilation-problem}

Problem 23 is not a slope. It is a direct shear test on one joint: two blocks are pressed
together, first at 3 MPa and then at 9 MPa, and one of them is dragged sideways. What the test
reports is shear stress against slip, a curve. There is no factor of safety in it and nothing to
lock. The manual uses the test to show the joint law working — the peak strength, the drop to a
residual strength once the joint has slipped, and dilation, the opening of the joint as it
slides.

The six vendor models are named for dilation angles of 0, 10, 20, 20, 20 and 30 degrees. Every
one of them carries `include_dilation: no`. The angle in a file's name never reaches the solver,
so all six ran at zero dilation. Two of the three 20-degree cases also step the normal load from
3 to 9 MPa at different points in the test, which makes them different tests. The manual's
dilation comparison never exercised dilation, so XSLOPE's dilation cannot be scored against it.

XSLOPE's dilation is checked against the kinematics instead. On a sliding joint the opening per
unit slip is the tangent of the dilation angle, which `test/joint_element_check.py` measures at
row 6. The peak and residual strengths are checked there as well, each against its closed form.

---

## The Rows

### 🟡 RJ-1a: Goodman & Bray block toppling, case a (rj001a) {#rj-1a}

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
[where a joint ends on another](#joint-terminations). They carry no cohesion, φ = 38.15° — the angle
this case is posed at — and the corpus's standard stiffness pair.

| XSLOPE SSRM | Goodman & Bray referee | UDEC | RS2 without / with improvement |
|---|---|---|---|
| **1.037** | 1.0 (+3.7%) | 0.99 | 0.99 / 0.97 |

<!-- test: file=files/rocscience/joints/rj001a.xlsx, type=fem_ssrm, expected_fs=1.037, element_type=tri6, target_size=10.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1a, f_stand=1.02734375, f_fail=1.046875, check=edges, tier=gate -->

The corpus bracket is decided on all nine trials and its longest takes 84 401 sweeps of the
250 000 allowed. A step of refinement — a 2D size of 7.0 m, which takes the mesh from 2 072
nodes and 67 interface elements to 3 573 and 94 — moves the factor by one bracket step, inside
the row's own tolerance, and every trial of that bracket reaches a verdict too.

The trial that used to stop this row was the finer bracket's upper edge. It reached the budget with
nothing to say; it now comes back FAILED at 225 001 sweeps on a steady slip. The four standing
trials beside it cost 303 to 309 sweeps each, where the same four cost 3 801 to 109 666 before.

The manual's figure for this case labels the side boundaries as rollers. The model's restraint list
does not: all 191 of its restrained nodes carry both components fixed, 111 of them on the base and
80 on the two sides. The model is what is transcribed. Its stabilizing toe force is 0.5 kN against a
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

| XSLOPE SSRM | Goodman & Bray referee | UDEC | RS2 without / with improvement |
|---|---|---|---|
| **1.018** | 1.0 (+1.8%) | 0.99 | 0.97 / 0.94 |

<!-- test: file=files/rocscience/joints/rj001b.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=10.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1b, f_stand=1.0078125, f_fail=1.02734375, check=edges, tier=gate -->

The standing edge of this bracket is the one trial in the corpus that reaches its sweep budget and
is answered anyway. At F = 1.007812 the viscoplastic loop runs the full 250 000 sweeps without
meeting its force tolerance, and the Newton corrector, handed that state at the budget exit,
reaches equilibrium in two iterations at an out-of-balance ratio of 1.1 × 10⁻¹¹ with no Gauss point
outside its surface. That is a stronger statement about the trial than the loop's own convergence
would have been, and the edge is read as decided on it — see [what answers a bracket
edge](#what-answers-a-bracket-edge).

A step of refinement — a 2D size of 7.0 m, which takes the mesh from 2 072 nodes and 67 interface
elements to 3 573 and 94 — returns the identical factor with every trial decided inside its budget,
the same standing edge taking 138 840 sweeps there. So the two meshes agree, and the finer one
needs no certification to say so.

Every input class the corpus transcribes was diffed against the vendor model and matches: the rock's
E, ν and γ, the joints' stiffness pair, cohesion, friction angle and tensile cap, the 31 contacts the
column outlines make, the side restraint the vendor clamps in both directions, and the force's
magnitude and direction.

**Input file:** [rj001b.xlsx](files/rocscience/joints/rj001b.xlsx).

![RJ-1b: Goodman & Bray block toppling, case b (rj001b) — FEM inputs with the 2013 kN toe force drawn at the block corner it acts on, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's — slip up every column contact and along the stepped base, the stack rotating forward over the face — reached on joints four degrees weaker, which is what the force at the toe buys](images/RJ-1b.png)

### 🟡 RJ-1c: Goodman & Bray block toppling, case c (rj001c) {#rj-1c}

Case a's section and rock with the joints half a degree steeper, φ = 38.6598°, and the same 0.5 kN
toe force that does nothing. Half a degree is worth about two points of factor of safety here, which
is what the closed form's own pair of answers for cases a and c says as well — 1.0 against 1.02.

| XSLOPE SSRM | Goodman & Bray referee | UDEC | RS2 without / with improvement |
|---|---|---|---|
| **1.057** | 1.02 (+3.6%) | 1.01 | 1.01 / 0.99 |

<!-- test: file=files/rocscience/joints/rj001c.xlsx, type=fem_ssrm, expected_fs=1.057, element_type=tri6, target_size=10.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1c, f_stand=1.046875, f_fail=1.06640625, check=edges, tier=gate -->

A step of refinement — a 2D size of 7.0 m, which takes the mesh from 2 072 nodes and 67
interface elements to 3 573 and 94 — moves the factor by one bracket step, inside the row's own
tolerance, and every trial of both brackets reaches a verdict. The corpus bracket's longest takes
82 261 sweeps of the 250 000 allowed, and the refinement's upper edge — the trial that held this
row back, which used to reach the budget with nothing to say — comes back FAILED at 225 001
sweeps on a steady slip.

The four standing trials of that refinement cost 302 to 311 sweeps each, where on the same mesh
and the same bracket they cost 11 196 to 131 095 before. Shortening a standing trial is what the
interface Newton corrector does, and it is why this row can be cut at all.

**Input file:** [rj001c.xlsx](files/rocscience/joints/rj001c.xlsx).

![RJ-1c: Goodman & Bray block toppling, case c (rj001c) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The mechanism is case a's at half a degree more friction: the same forward rotation of the column stack over the face, and the same slip on every column contact and along the stepped base](images/RJ-1c.png)

### 🟢 RJ-1d: Goodman & Bray block toppling, case d (rj001d) {#rj-1d}

Case c's joint friction angle with case b's 2013 kN force: the same stack, stabilized. It is the
strongest of the four, and the closed form and UDEC agree that it is.

| XSLOPE SSRM | Goodman & Bray referee | UDEC | RS2 without / with improvement |
|---|---|---|---|
| **1.252** | 1.23 (+1.8%) | 1.22 | 1.19 / 1.16 |

<!-- test: file=files/rocscience/joints/rj001d.xlsx, type=fem_ssrm, expected_fs=1.252, element_type=tri6, target_size=10.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-1d, f_stand=1.2421875, f_fail=1.26171875, check=edges, tier=gate -->

The standing edge of this bracket, F = 1.242188, runs the full 250 000 sweeps without meeting its
force tolerance, and the corrector reaches equilibrium from that state at the budget exit, which is
what answers it — see [what answers a bracket edge](#what-answers-a-bracket-edge). The plain
viscoplastic loop read the same trial as failing there, so the factor is one bracket step above
what that loop gave.

A step of refinement — a 2D size of 7.0 m, which takes the mesh from 2 072 nodes and 67 interface
elements to 3 573 and 94 — returns the identical bracket edge for edge, with every trial decided
inside its budget and the same standing edge taking 191 554 sweeps there. Two meshes, one answer,
and the finer one needs no certification to give it.

Every transcribed input class matches the vendor model, the force reaching the same point case b's
does.

**Input file:** [rj001d.xlsx](files/rocscience/joints/rj001d.xlsx).

![RJ-1d: Goodman & Bray block toppling, case d (rj001d) — FEM inputs with the 2013 kN toe force, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. As on case b the right-hand panels are one sweep past the critical factor rather than a developed mechanism, so they show where the stack starts to move rather than where it ends up](images/RJ-1d.png)

### 🔴 RJ-2: Alejano & Alonso block toppling (rj002) {#rj-2}

A 30 × 19.85 m section whose 9.85 m face rises at 58.65° from (20, 10) to (14, 19.85). A basal
joint runs from the toe of that face up to the crest at (2.9393, 19.85) at 30° — the stepped
surface the columns stand on — and twenty-two columns at 64°, 1.6 m apart, pass through the same
toe. Both carry φ = 31° and no cohesion. The rock is elastic (γ = 25 kN/m³, E = 20 GPa, ν = 0.3),
the vendor's own `Plasticity Specifications: Non`, so the columns cannot yield and every mechanism
the model has is a joint one.

| XSLOPE SSRM | UDEC referee | LE (Goodman & Bray) | RS2 without / with improvement |
|---|---|---|---|
| **0.764** | 0.87 (−12.2%) | 0.76 | 0.86 / 0.82 |

<!-- test: file=files/rocscience/joints/rj002.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-2, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.35 m, which takes the mesh from 11 537 nodes to 21 110 —
moves the factor by one bracket step, inside the row's own tolerance, and all nine trials of both
brackets reach a verdict. The longest takes 185 381 sweeps of the 250 000 allowed.

**The gap to the referee is 12%, and the three continuum codes do not bracket this row.** Alejano
& Alonso publish Goodman & Bray's limit equilibrium at 0.76 and their own UDEC run at 0.87; RS2
reports 0.86 and 0.82, so both of its own factors sit above this one rather than straddling it.
XSLOPE lands on the limit equilibrium, and RS2 and UDEC land together above it. That is a
difference between codes on one problem, stated here because it is measured, and it is not
explained: every input class the corpus transcribes was diffed against the vendor model and
matches. The joints reconstruct from the vendor's own 517 joint elements to 311.4 m of trace,
291.7 m at 64° and 19.7 m at 30°, against this file's 311.4 m at the same angles over the same
extent; the rock's E, ν and γ, the joints' k<sub>n</sub>, k<sub>s</sub>, c, φ and tensile cap, and
the clamped side restraint all match.

What is not transcribed is a class of RS2 solver switch that every model in this manual carries and
the corpus has never read: a joint stiffness recalculated when a joint violates its strength
criterion, which the manual's own release notes describe as the "joint improvement" option's
mechanism. Whether it bears on this gap is untested.

**Input file:** [rj002.xlsx](files/rocscience/joints/rj002.xlsx).

![RJ-2: Alejano & Alonso block toppling (rj002) — FEM inputs, mesh, joint slip at the critical SRF and the section deformed 17×. The rock carries no strain of its own because it cannot yield, so every movement in the section is on a joint: slip gathers where the basal joint reaches the toe of the face, the columns standing on it open along their upper halves, and the deformed section shows them rotating out over the face while the rock below the basal joint stays put](images/RJ-2.png)

### ⊘ RJ-3: Lorig & Varona forward block toppling (rj003) {#rj-3}

The 260 m section at 55° that problems 3 to 7 share, cut by two sets: columns at 70° at 20 m
spacing and a cross set at −20° at 30 m, both through the origin. The manual states the pair as
"70 and 160" degrees, which is the same two planes measured the other way round the half circle.
The rock is elastic — the vendor's `Plasticity Specifications: Non` — so only the joints can fail;
γ = 26.0946 kN/m³, E = 9072 MPa, ν = 0.26. The joints carry c = 100 kPa and φ = 40°.

The row carries no lock, so it prints no factor of its own. The upper edge of its final bracket
reaches 250 000 sweeps without a verdict, and a bracket edge nothing ruled on cannot define a
factor of safety.

Every input class the corpus transcribes was diffed against the vendor model and matches: the
rock's E, ν and γ, the joints' normal and shear stiffness, cohesion, friction angle and tensile
cap, the two sets' dips and spacings, and the side restraint the vendor clamps in both directions.

**Input file:** [rj003.xlsx](files/rocscience/joints/rj003.xlsx).

![RJ-3: Lorig & Varona forward block toppling (rj003) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own, so the whole mechanism is on the two sets: the steep 70° joints slip behind the crest while the −20° cross joints open along it, and the deformed section shows the blocks between them rotating forward over the face](images/RJ-3.png)

### ⊘ RJ-4: Lorig & Varona flexural toppling (rj004) {#rj-4}

The same section cut by one set of columns at 70° at 20 m spacing — problem 3's first set without
its cross joints, so the columns bend rather than topple as blocks. Here the rock is Mohr-Coulomb
and carries a tensile cutoff of zero, which is what lets a column break in flexure: γ = 26.1 kN/m³,
E = 9072 MPa, ν = 0.26, c = 675 kPa, φ = 43°.

The corpus bracket is decided on all nine trials. What keeps the row reported is the refinement
step: one of its trials reaches 250 000 sweeps without a verdict, so the claim that a finer mesh
does not move the factor is not established, even though the two meshes bracket values one step
apart. The row prints no factor.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj004.xlsx](files/rocscience/joints/rj004.xlsx).

![RJ-4: Lorig & Varona flexural toppling (rj004) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. Here the rock can yield, and it does: a band of shear strain climbs from the toe across the columns, and drawn without exaggeration the columns are bent through that band rather than rotated about it, which is what separates flexural toppling from the block toppling of problem 3](images/RJ-4.png)

### ⊘ RJ-5: Lorig & Varona backward block toppling (rj005) {#rj-5}

The shared 260 m section cut by two sets: one at −55° at 10 m spacing through the toe at
(560, 140), dipping out of the face so the blocks lean back rather than forward, and a horizontal
set at 40 m spacing. The rock is elastic — the vendor's `Plasticity Specifications: Non` —
γ = 26.1 kN/m³, E = 9072 MPa, ν = 0.26; the joints carry c = 100 kPa and φ = 40°.

This is the most budget-limited row in the corpus. **Four** of its nine trials reach 250 000 sweeps
without a verdict, and two of those four are the edges of the final bracket, so the factor the
bracket encloses is conditioned on the sweep limit at both ends. The row prints no factor and
carries no lock. It is also the corpus's longest run at about 5.7 hours on the corpus mesh.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj005.xlsx](files/rocscience/joints/rj005.xlsx).

![RJ-5: Lorig & Varona backward block toppling (rj005) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The rock is elastic and carries no strain of its own; slip runs along the −55° joints in the wedge behind the face while the horizontal bedding opens, and the deformed section shows the slabs stepping out over one another down the face, each leaning back into the slope as it goes](images/RJ-5.png)

### ⊘ RJ-6: Plane failure with daylighting discontinuities (rj006) {#rj-6}

The same section cut by one set at −35° at 10 m spacing through the origin. The joints dip out of
the 55° face at a shallower angle than the face itself, so every one of them daylights and the
slabs between them are free to slide out. The rock is Mohr-Coulomb (γ = 26.1 kN/m³, E = 9072 MPa,
ν = 0.26, c = 675 kPa, φ = 43°, no tensile capacity).

Both edges of the final bracket reach 250 000 sweeps without a verdict — one of them
`STABLE_STUCK`, the other stopped at the cap — so the factor the bracket encloses is a statement
about the sweep budget as much as about the slope. The row prints no factor and carries no lock.

Every transcribed input class matches the vendor model, including the side restraint.

**Input file:** [rj006.xlsx](files/rocscience/joints/rj006.xlsx).

![RJ-6: plane failure with daylighting discontinuities (rj006) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Slip runs the full length of every joint that reaches the face, over a wedge bounded below by the joint through the toe, and the rock between them carries only a faint strain: the slabs slide out along the joints rather than breaking through anything](images/RJ-6.png)

### ⊘ RJ-7: Plane failure with non-daylighting discontinuities (rj007) {#rj-7}

The same section and the same rock as problem 6, cut by one set at −70° at 20 m spacing through
the origin. The joints now dip out of the face more steeply than the 55° face itself, so none of
them daylights: a slab cannot slide out along one without shearing rock, and the slope stands
higher than problem 6's.

The corpus bracket is decided on all nine trials. What keeps the row reported is the refinement
step, whose lower bracket edge reaches 250 000 sweeps without a verdict — so although the two
meshes enclose the same factor, the claim that refinement does not move it rests on an edge nothing
ruled on. The row prints no factor.

Every transcribed input class matches the vendor model, including the side restraint. The manual's
own table for this problem prints the slope angle as 5°; the figure and the model are the same 55°
slope problem 6 uses.

**Input file:** [rj007.xlsx](files/rocscience/joints/rj007.xlsx).

![RJ-7: plane failure with non-daylighting discontinuities (rj007) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. No joint reaches the face at a shallower angle than the face itself, so the failure cannot slide out along one: a strong band of shear strain cuts across the steep joints from the toe, and the mass above it moves out over rock it has had to break](images/RJ-7.png)

### 🟢 RJ-8: Flexural toppling in a base friction model (rj008) {#rj-8}

Pritchard & Savigny's base-friction table model, scaled up a hundred times: a 72.4 × 36.5 m
section whose 30.5 m face rises at 78° from (15, 6) to (21.48, 36.5). Thirteen columns at −60°,
5.08 m apart, stand on a horizontal joint at y = 6 that runs the width of the model, with a
vertical joint at x = 68.4 closing the back of the stack. The rock is Mohr-Coulomb (γ = 25.506
kN/m³, E = 22.771 GPa, ν = 0.139, c = 60 kPa, φ = 39°). The joints carry no cohesion, φ = 39°, and
the softest normal stiffness in the corpus bar one: k<sub>n</sub> = 1.5 × 10<sup>7</sup> kPa/m
against the set's usual 10<sup>8</sup>.

| XSLOPE SSRM | UDEC referee | RS2 without / with improvement |
|---|---|---|
| **0.764** | 0.76 (+0.5%) | 0.75 / 0.75 |

<!-- test: file=files/rocscience/joints/rj008.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-8, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 5 547 nodes to 10 802 —
moves the factor by one bracket step, inside the row's own tolerance, and all nine trials of both
brackets reach a verdict. The budget is not what settles this row: its longest trial takes 1 981
sweeps of the 250 000 allowed, where the [geotextile wall family](rs2.md#rs2-48) exhausts that
budget on five trials of eight rows.

The two meshes disagree about one trial, and it is the one the bracket closes on. At F = 0.753906
the corpus mesh reaches equilibrium in 311 sweeps and the finer mesh runs away in 1 981, so the
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
[where a joint ends on another](#joint-terminations).

| XSLOPE SSRM | UDEC referee | Alejano's limit equilibrium | RS2 without / with improvement |
|---|---|---|---|
| **1.018** | 1.03 (−1.2%) | 0.40–1.45 | 1.01 / 1.09 |

<!-- test: file=files/rocscience/joints/rj009.xlsx, type=fem_ssrm, expected_fs=1.018, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-9, f_stand=1.0078125, f_fail=1.02734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 2.1 m, which takes the mesh from 12 620 nodes and 1 806
interface elements to 26 565 and 2 579 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
78 928 sweeps of the 250 000 allowed. This is the row that shows the family can be
mesh-independent at the corpus size, where problems 12 and 14 are not.

The paper's own limit equilibrium for this problem is a range rather than an answer — 0.40 to 1.45,
as the manual prints it — and it is recorded beside the referee rather than scoring, as every source
limit equilibrium here is. The distinct-element run is what the manual verifies RS2 against.

**Input file:** [rj009.xlsx](files/rocscience/joints/rj009.xlsx).

![RJ-9: Alejano et al. bilinear slab failure, example 1a (rj009) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: the release trace under the crest carries the brightest slip, a bedding plane below the toe carries the rest, and the deformed section shows the slab between them sliding out over the bench](images/RJ-9.png)

### 🟢 RJ-11: Alejano et al. ploughing sliding slab failure (rj011) {#rj-11}

A 25 m slope at 50° with bedding dipping out of the face at −50° at 1.5 m spacing, φ = 30°, and two
release traces at φ = 20°: one at the toe running below the bench and one from the face down onto
the bedding plane that releases the toe block. Ploughing failure is the paper's name for sliding on
a primary discontinuity combining with sliding on a joint sub-parallel to the face, which lifts the
toe block and eventually rotates it out of the slope. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The lower tip of a release trace is stated to six
decimals, so it lands a part in 10⁷ from the bedding plane it belongs on; see
[where a joint ends on another](#joint-terminations).

| XSLOPE SSRM | UDEC referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|
| **1.213** | 1.21 (+0.2%) | 1.75 | 1.22 / 1.3 |

<!-- test: file=files/rocscience/joints/rj011.xlsx, type=fem_ssrm, expected_fs=1.213, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-11, f_stand=1.203125, f_fail=1.22265625, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 12 326 nodes and 1 768
interface elements to 26 052 and 2 524 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
61 171 sweeps of the 250 000 allowed.

**The manual's two published values are half a factor apart, and this row lands on the one that
scores it.** Alejano et al. publish a limit-equilibrium 1.75 and a UDEC run at 1.21 for the same
problem. UDEC is the referee — it is what the manual verifies RS2 against here, and it is the
same-method pairing — and RS2's own two factors, 1.22 and 1.3, sit beside it rather than beside the
limit equilibrium, with this row between RS2's lower factor and UDEC's. Against the source paper's
limit equilibrium the same factor reads −30.7%, which is a cross-method difference and scores
nothing. [Problem 15](#rj-15) has the same shape, with three published values on one problem.

Every input class the corpus transcribes was diffed against the vendor model and matches: the rock's
E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile cap, the
bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor clamps in
both directions. Nothing here explains the gap to the limit-equilibrium value; it is recorded as
measured.

**Input file:** [rj011.xlsx](files/rocscience/joints/rj011.xlsx).

![RJ-11: Alejano et al. ploughing sliding slab failure (rj011) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. One bedding plane from the crest to the toe carries almost all of the slip, and at its foot the two release traces cut out a small wedge: the deformed section shows that wedge lifted and rotated out over the bench while the slab above it slides down the plane, which is the ploughing mechanism the paper names](images/RJ-11.png)

### ⊘ RJ-12: Alejano et al. ploughing toppling slab failure (rj012) {#rj-12}

A 25 m slope at 60° cut by bedding dipping out of the face at −60° at 1.5 m spacing, φ = 30°, with
two short release traces at φ = 40°: one at the toe running below the bench, and one from the face
down onto the bedding plane that releases the toe block. Ploughing failure is the paper's name for
sliding on a primary discontinuity combining with sliding on a joint sub-parallel to the face, which
lifts the toe block and rotates it out; at 60° the rotation rather than the sliding governs, which
is what separates this row from problem 11. The rock is elastic at E = 2 × 10⁸ MPa, γ = 25 kN/m³,
ν = 0.3 — not a rock modulus but the manual's own device for reproducing UDEC's rigid blocks.

The lower tip of each release trace is stated to six decimals, so it lands a part in 10⁷ from the
bedding plane it belongs on; see [where a joint ends on another](#joint-terminations). This row and
problem 14 state their generated bedding network to six decimals as well, which is the precision
the mesher's own crossing arithmetic carries — see the departures table.

The row carries no lock and prints no factor, and what holds it back is the mesh rather than the
sweep budget. Every trial of both brackets reaches a verdict — the longest takes 16 061 sweeps of
the 250 000 allowed, which is the least budget-limited row in this corpus — but a step of
refinement to a 2D size of 1.05 m, which takes the mesh from 12 294 nodes and 1 818 interface
elements to 26 407 and 2 572, moves the factor by 0.27. That is fourteen bracket steps, where the
row's own tolerance is one. The corpus mesh size for this family is the bedding spacing, which puts
about one element between one bedding plane and the next, and this row says that is not enough to
settle the factor.

The manual publishes two referees for this problem and they do not agree with each other, with
RS2's own two factors falling between them.

**Input file:** [rj012.xlsx](files/rocscience/joints/rj012.xlsx).

![RJ-12: Alejano et al. ploughing toppling slab failure (rj012) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section, and only a few of its traces carry any slip: one release trace under the crest and the bedding beneath the toe block, which is the ploughing pair. The two right-hand panels are the state one sweep past the critical factor rather than a developed mechanism — on a rock this stiff the model moves by microns until it does not, so the deformed section is drawn at tens of thousands of times scale](images/RJ-12.png)

### 🟢 RJ-13: Alejano et al. ploughing sliding slab, example 4 (rj013) {#rj-13}

A 25 m slope at 55° cut by bedding dipping out of the face at −55° at 1.5 m spacing, φ = 25°, with
two release traces at φ = 20°: one at the toe running below the bench and one from the face down
onto the bedding plane that releases the toe block. It is problem 11's mechanism at a steeper
bedding and a weaker one — sliding on the primary discontinuity combining with sliding on a joint
sub-parallel to the face, which lifts the toe block. The rock is the family's elastic rigid-block
stand-in at E = 2 × 10⁸ MPa, γ = 25 kN/m³, ν = 0.3. The lower tip of each release trace is stated to
six decimals, so it lands a part in 10⁷ from the bedding plane it belongs on; see
[where a joint ends on another](#joint-terminations).

| XSLOPE SSRM | UDEC referee | LE (Alejano) | RS2 without / with improvement |
|---|---|---|---|
| **0.998** | 1.0 (−0.2%) | 1.0 | 1.0 / 1.05 |

<!-- test: file=files/rocscience/joints/rj013.xlsx, type=fem_ssrm, expected_fs=0.998, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-13, f_stand=0.98828125, f_fail=1.0078125, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m, which takes the mesh from 12 141 nodes and 1 788
interface elements to 26 079 and 2 539 — does not move the factor at all: the finer mesh returns
the same bracket, edge for edge. Every trial of both brackets reaches a verdict, the longest taking
23 891 sweeps of the 250 000 allowed.

This is the third row of the Alejano family to hold its factor under refinement at the corpus mesh
size, against two that move ([12](#rj-12) and [14](#rj-14)). Its bedding spacing is 1.5 m, the same
as both movers and as [problem 11](#rj-11), which holds.

Every input class the corpus transcribes was diffed against the vendor model and matches: the
rock's E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile
cap, the bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor
clamps in both directions.

**Input file:** [rj013.xlsx](files/rocscience/joints/rj013.xlsx).

![RJ-13: Alejano et al. ploughing sliding slab, example 4 (rj013) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it carries slip: one bedding plane from the crest down to the toe does, and the two release traces at its foot cut out the toe block, which the deformed section shows lifted and pushed out over the bench while the slab above it slides down the plane](images/RJ-13.png)

### ⊘ RJ-14: Alejano et al. ploughing sliding slab, example 5 (rj014) {#rj-14}

Example 4's section at 60° with the two joint strengths the other way round: bedding dipping out of
the face at −60° at 1.5 m spacing at φ = 20° — the weakest bedding of the six — and two release
traces at φ = 30°, one at the toe below the bench and one from the face down onto the bedding plane
it releases. The rock is the family's elastic rigid-block stand-in at E = 2 × 10⁸ MPa,
γ = 25 kN/m³, ν = 0.3. Like problem 12, this row states its generated bedding network to six
decimals; see the departures table.

The row carries no lock and prints no factor, and as on problem 12 it is the mesh rather than the
sweep budget that says so. Every trial of both brackets reaches a verdict — the longest takes
32 926 sweeps of the 250 000 allowed — and a step of refinement to a 2D size of 1.05 m, which
takes the mesh from 12 316 nodes and 1 819 interface elements to 26 414 and 2 572, moves the
factor by 0.059. That is three bracket steps where the row's tolerance is one.

Every input class the corpus transcribes was diffed against the vendor model and matches: the
rock's E, ν and γ, both joint friction angles, the joints' stiffness pair, cohesion and tensile
cap, the bedding dip and spacing, the release traces' endpoints, and the side restraint the vendor
clamps in both directions.

**Input file:** [rj014.xlsx](files/rocscience/joints/rj014.xlsx).

![RJ-14: Alejano et al. ploughing sliding slab, example 5 (rj014) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. The bedding set runs the whole section and almost none of it moves: one bedding plane from the crest to the toe carries the slip, with the release trace at the toe opening as the slab above it slides out over the bench](images/RJ-14.png)

### ⊘ RJ-15: Partially joint-controlled footwall slope (rj015) {#rj-15}

A 25 m footwall at 40° whose bedding dips in the same direction at the same angle, 2 m apart, so
the slabs lie parallel to the face and a failure has to break rock at the toe to get out. This is
the one problem in the Alejano family whose rock can yield: Mohr-Coulomb, c = 200 kPa, φ = 35°,
γ = 28 kN/m³, E = 1 GPa, ν = 0.3. The joints are the corpus's softest — k<sub>n</sub> = 5 × 10<sup>6</sup>
kPa/m and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m, twenty times below the set's standard pair — with
no cohesion and φ = 25°.

The manual publishes five factors for this problem, from one paper and three programs:

| UDEC-SSRT referee (Alejano) | LE (Alejano) | RS2 without / with improvement | Slide2 LEM (vendor) |
|---|---|---|---|
| 1.6 | 1.72 | 1.28 / 1.42 | 1.25 |

The first two are the source paper's own answers, its UDEC run and its limit equilibrium, and the
UDEC run is what scores this row: it is what the manual verifies RS2 against, and it is the
same-method pairing. The last is Rocscience's own companion program on its own model, labeled and
recorded rather than scoring. The five run from 1.25 to 1.72, a wider spread than separates any two
codes anywhere else in this corpus.

The row carries no lock, so it prints no factor of its own. The upper edge of its final bracket
reaches 250 000 sweeps without a verdict, and a bracket edge nothing ruled on cannot define a
factor of safety.

Every input class that was transcribed matches the vendor model: the rock's E, ν, γ, c, φ and
tensile cap; the joints' normal and shear stiffness, cohesion, friction angle and tensile cap; and
the side restraint the vendor clamps in both directions. The stated tensile strength of 1000 kPa is
above the Mohr-Coulomb apex its own c and φ imply (c/tan φ = 285.6 kPa), so it never binds.

**Input file:** [rj015.xlsx](files/rocscience/joints/rj015.xlsx).

![RJ-15: partially joint-controlled footwall slope (rj015) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The bedding slips over a long stretch behind the face, and the rock's only strain is a small patch at the toe where the slab has to break through to get out — the coupled mechanism the source paper describes](images/RJ-15.png)

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

![RJ-18: step-path failure through three continuous joints (rj018) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Every joint is slipping along its lower half and open along its upper, and the three slabs between them slide out together down the 36.1° path; the rock itself carries almost no plastic strain](images/RJ-18.png)

### ⊘ RJ-19: Bi-planar step-path failure (rj019) {#rj-19}

A 120 × 70 m section of one Mohr-Coulomb rock (γ = 27 kN/m³, E = 20 GPa, ν = 0.3, c = 10 500 kPa,
φ = 35°, tensile capacity 200 kPa) with a slope face from (30, 20) to (60, 70), cut by two
discontinuous joints with a rock bridge between them: a basal joint at 28.4° from (39.0149,
35.0248) to (63, 48), and an upper joint at 56.3° from (62, 49) to (76, 70). Both carry c = 0,
φ = 40° and the same stiffness pair as RJ-18. Referee: UDEC 1.46.

A step of refinement — a 2D size of 2.1 m, which takes the mesh from 3 485 nodes to 6 910 —
returns the same bracket edge for edge, so the factor this row brackets is not its mesh's. What it
still is, in part, is its budget's: the upper edge of both brackets reaches 250 000 sweeps without a
verdict, and the row is reported without a lock.

The lower edge is where the two meshes differ, and it is the sharper measurement. On the corpus
mesh it converges at 250 003 sweeps — three past the budget, so the trial record reads it as
undecided — and on the finer mesh the same trial converges in 348. A bracket edge that needs a
quarter of a million sweeps at one mesh and three hundred at a finer one is not a slope that is
close to failing there.

The manual's table for this problem states one joint inclination as 59°. Its figure dimensions 56°
and 28°, and the vendor model's own endpoints give 56.3° and 28.4°, so the table is the outlier and
the model is what is built.

**Input file:** [rj019.xlsx](files/rocscience/joints/rj019.xlsx).

![RJ-19: bi-planar step-path failure with a rock bridge (rj019) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section. Both joints have opened, the basal one is slipping at its lower end, and the only strain in the rock is the patch at the bridge between the two joint tips, where the block above has to break through to move](images/RJ-19.png)

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
  a and c, which is 0.05% of it and does nothing. The force is horizontal and points into the
  slope in every one of the six models. The vendor's own rerun models move it from the toe at
  (−0.5, 0.866) to the upper-left block corner at (−2.5, 4.330) in all four cases, so the "with
  joint improvement" factors are not the same load case as the "without" ones. Cases b and d here
  carry the force at the **corner**: the toe is the end of the basal joint, where the split gives
  the node one copy per wedge of material around it and a load applied there has no defined side to
  act on. Cases a and c carry no load, the 0.5 kN being what it is.
- **Problem 20's network.** The manual calls it a Voronoi tessellation generated in UDEC and
  imported. The vendor model carries it as 523 digitized traces with no block size, density or seed
  recorded, so nothing published reproduces it.
- **Problem 16's boundary conditions.** The 0° case runs on rollers; the nine tilted cases pin
  every exterior node in both directions, and use a convergence tolerance two orders tighter.
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
- **Problem 8's joint network is not the whole section.** Its thirteen columns are clipped to the
  block above the basal joint; generated over the section they would be sixteen, three of them
  under a plane the model has no columns below.
- **Problems 9 to 14 run on a rock a thousand times stiffer than steel.** E = 2 × 10⁸ MPa is the
  manual's own device, and it says so: the UDEC models these are scored against use rigid blocks,
  and the modulus is how RS2 reproduces one. Problem 15, whose rock can yield, runs at 1 GPa.
- **No problem from 1 to 21 states a joint residual strength or a dilation angle.** Those appear
  only in problems 22 and 23, and material residual values, where they appear, always equal the
  peak.
