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

| # | Match | Problem | Referee | RS2 without / with improvement | Notes |
|---:|:-:|---|---|---|---|
| 1a | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1a | Goodman & Bray 1.0 · UDEC 0.99 | 0.99 / 0.97 | *blocked* — the stepped base puts each column's basal contact partway along its neighbour's side joint; see [where a joint ends on another](#joint-terminations). |
| 1b | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1b | Goodman & Bray 1.0 · UDEC 0.99 | 0.97 / 0.94 | *blocked* — as 1a; its 2013 kN stabilizing toe force enters as a line load on the `lloads` sheet. |
| 1c | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1c | Goodman & Bray 1.02 · UDEC 1.01 | 1.01 / 0.99 | *blocked* — as 1a. |
| 1d | <span class="nodata">⊘</span> | Goodman & Bray block toppling, case 1d | Goodman & Bray 1.23 · UDEC 1.22 | 1.19 / 1.16 | *blocked* — as 1b, with the same 2013 kN toe force. |
| [2](#rj-2) | 🟢 | Alejano & Alonso block toppling | SSRM 0.764 vs Goodman 0.76 (+0.5%) | 0.86 / 0.82 | |
| [3](#rj-3) | <span class="nodata">⊘</span> | Lorig & Varona forward block toppling | UDEC 1.13 | 1.12 / 1.09 | *reported, no lock* — the upper edge of the bracket reaches the sweep budget without a verdict. |
| [4](#rj-4) | <span class="nodata">⊘</span> | Lorig & Varona flexural toppling | UDEC 1.3 | 1.19 / 1.27 | *reported, no lock* — one trial of the refinement step reaches the sweep budget without a verdict. |
| [5](#rj-5) | <span class="nodata">⊘</span> | Lorig & Varona backward block toppling | UDEC 1.7 | 1.65 / 1.86 | *reported, no lock* — four trials of the bracket, including both its edges, reach the sweep budget without a verdict. |
| [6](#rj-6) | <span class="nodata">⊘</span> | Plane failure, daylighting | UDEC 1.27 | 1.25 / 1.31 | *reported, no lock* — both edges of the bracket reach the sweep budget without a verdict. |
| [7](#rj-7) | <span class="nodata">⊘</span> | Plane failure, non-daylighting | UDEC 1.5 | 1.57 / 1.59 | *reported, no lock* — one trial of the refinement step reaches the sweep budget without a verdict. |
| [8](#rj-8) | 🟢 | Flexural toppling, base friction model | SSRM 0.744 vs UDEC 0.76 (−2.1%) | 0.75 / 0.75 | |
| 9 | <span class="nodata">⊘</span> | Bilinear slab failure, example 1a | UDEC 1.03 | 1.01 / 1.09 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 10 | <span class="nodata">⊘</span> | Bilinear slab failure, example 1b | UDEC 1.03 | 0.92 / 1.08 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 11 | <span class="nodata">⊘</span> | Ploughing sliding slab failure | UDEC 1.21 | 1.22 / 1.3 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 12 | <span class="nodata">⊘</span> | Ploughing toppling slab failure | UDEC 1.78 | 1.39 / 1.75 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 13 | <span class="nodata">⊘</span> | Ploughing sliding slab, example 4 | UDEC 1.0 | 1.0 / 1.05 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| 14 | <span class="nodata">⊘</span> | Ploughing sliding slab, example 5 | UDEC 0.9 | 0.89 / 1.09 | *blocked* — the release trace ends on a bedding plane; see [where a joint ends on another](#joint-terminations). |
| [15](#rj-15) | <span class="nodata">⊘</span> | Partially joint-controlled footwall | LE (Alejano) 1.72 · UDEC 1.6 · Slide2 1.25 | 1.28 / 1.42 | *reported, no lock* — the upper edge of the bracket reaches the sweep budget without a verdict. |
| 16 | <span class="nodata">⊘</span> | Barla et al. tilt-table block toppling | experiment 9° · UDEC 11° | 9° / 7° | *blocked* — the problem scores the TILT ANGLE at which a block grid topples, found by rotating gravity through a staged sweep; XSLOPE's seismic coefficient tilts the load but the row needs the sweep and a toppling criterion, neither of which is a strength reduction. |
| 17 | <span class="nodata">⊘</span> | Step-path, en-echelon joints | UDEC 1.29 | 1.24 / 1.2 | *blocked* — the vendor model's second, elastic material has no boundary of its own. Recovered from the element-material map it is an 83-vertex staircase of element edges, trending vertical near x = 12 and horizontal near y = 3 with excursions of about one element either side, so the boundary is a property of the vendor's mesh rather than of its model. |
| [18](#rj-18) | 🟢 | Step-path, continuous joints | SSRM 0.998 vs UDEC 1.01 (−1.2%) | 1.01 / 1.0 | |
| [19](#rj-19) | <span class="nodata">⊘</span> | Bi-planar step-path failure | UDEC 1.46 | 1.5 / 1.41 | *reported, no lock* — two trials of the bracket reach the sweep budget without a verdict. |
| 20 | <span class="nodata">⊘</span> | Hammah & Yacoub Voronoi slope | UDEC 2.46 | 2.21 / 2.37 | *blocked* — the manual states no block size and no seed, and the vendor model carries the tessellation as 523 digitized traces rather than as a generated network, so the input is not reproducible from anything published. |
| 21 | <span class="nodata">⊘</span> | Shallow excavation, jointed tunnel | UDEC 8.16 | 8.27 / 8.5 | *blocked* — the second stage of the vendor model excavates a 2 m opening and the strength reduction runs on the excavated state, which carries the stress the first stage left behind; XSLOPE has no staged construction. |
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
| **0.744** | 0.76 (−2.1%) | 0.75 / 0.75 |

<!-- test: file=files/rocscience/joints/rj008.xlsx, type=fem_ssrm, expected_fs=0.744, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-8, f_stand=0.734375, f_fail=0.75390625, check=edges, tier=gate -->

A step of refinement — a 2D size of 1.05 m — moves the factor by exactly zero, and all nine trials
of both brackets reach a verdict. The budget is not what settles this row: its longest trial takes
3 639 sweeps of the 250 000 allowed, where the [geotextile wall family](rs2.md#rs2-48) exhausts
that budget on five trials of eight rows.

The model's stated tensile strength of 75 kPa is above the Mohr-Coulomb apex its own c and φ imply
(c/tan φ = 74.1 kPa), so the cap never binds and the envelope's own apex governs.

**Input file:** [rj008.xlsx](files/rocscience/joints/rj008.xlsx).

![RJ-8: flexural toppling in a base friction model (rj008) — FEM inputs, mesh, viscoplastic shear strain with joint slip at the critical SRF, and the deformed section at true scale. The strain gathers into one lobe per column along a band that climbs from the toe across the stack; above that band the columns are visibly bent rather than merely tilted, and the rock below it is unstrained. That is the break surface of flexural toppling rather than sliding along any one joint](images/RJ-8.png)

### 🟢 RJ-2: Alejano & Alonso block toppling (rj002) {#rj-2}

A 30 × 19.85 m section whose 9.85 m face rises at 58.65° from (20, 10) to (14, 19.85). A basal
joint runs from the toe of that face up to the crest at (2.9393, 19.85) at 30° — the stepped
surface the columns stand on — and twenty-two columns at 64°, 1.6 m apart, pass through the same
toe. Both carry φ = 31° and no cohesion. The rock is elastic (γ = 25 kN/m³, E = 20 GPa, ν = 0.3),
the vendor's own `Plasticity Specifications: Non`, so the columns cannot yield and every mechanism
the model has is a joint one.

| XSLOPE SSRM | Goodman & Bray referee | UDEC | RS2 without / with improvement |
|---|---|---|---|
| **0.764** | 0.76 (+0.5%) | 0.87 | 0.86 / 0.82 |

<!-- test: file=files/rocscience/joints/rj002.xlsx, type=fem_ssrm, expected_fs=0.764, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=3.0, max_iter=250000, tension_srf=false, k0=1, benchmark=RJ-2, f_stand=0.75390625, f_fail=0.7734375, check=edges, tier=gate -->

A step of refinement — a 2D size of 0.35 m, which takes the mesh from 11 537 nodes to 21 110 —
moves the factor by one bracket step, inside the row's own tolerance, and all nine trials of both
brackets reach a verdict. The longest takes 185 381 sweeps of the 250 000 allowed.

**The two other codes agree with each other, and not with this row.** Alejano & Alonso publish
Goodman & Bray's limit equilibrium at 0.76 and their own UDEC run at 0.87; RS2 reports 0.86 and
0.82. XSLOPE lands on the limit equilibrium, and RS2 and UDEC land together above it. That is a
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

### ⊘ RJ-15: Partially joint-controlled footwall slope (rj015) {#rj-15}

A 25 m footwall at 40° whose bedding dips in the same direction at the same angle, 2 m apart, so
the slabs lie parallel to the face and a failure has to break rock at the toe to get out. This is
the one problem in the Alejano family whose rock can yield: Mohr-Coulomb, c = 200 kPa, φ = 35°,
γ = 28 kN/m³, E = 1 GPa, ν = 0.3. The joints are the corpus's softest — k<sub>n</sub> = 5 × 10<sup>6</sup>
kPa/m and k<sub>s</sub> = 5 × 10<sup>5</sup> kPa/m, twenty times below the set's standard pair — with
no cohesion and φ = 25°.

The manual publishes five factors for this problem, from one paper and three programs:

| LE (Alejano) referee | UDEC-SSRT (Alejano) | RS2 without / with improvement | Slide2 LEM |
|---|---|---|---|
| 1.72 | 1.6 | 1.28 / 1.42 | 1.25 |

The first two are the source paper's own answers, its limit equilibrium and its UDEC run, and the
limit equilibrium is what scores this row. The last is Rocscience's own companion program run on
its own model, recorded for completeness and not used as a referee. The five run from 1.25 to
1.72, a wider spread than separates any two codes anywhere else in this corpus.

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

Two trials of the bracket reach 250 000 sweeps without a verdict, so the factor this variant
brackets is a statement about the budget as much as about the slope, and the row is reported
without a lock.

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
  a and c, which is 0.05% of it and does nothing. The vendor's own rerun models move the force from
  the toe at (−0.5, 0.866) to the upper-left block corner at (−2.5, 4.330) in all four cases, so
  the "with joint improvement" factors are not the same load case as the "without" ones.
- **Problem 20's network.** The manual calls it a Voronoi tessellation generated in UDEC and
  imported. The vendor model carries it as 523 digitized traces with no block size, density or seed
  recorded, so nothing published reproduces it.
- **Problem 16's boundary conditions.** The 0° case runs on rollers; the nine tilted cases pin
  every exterior node in both directions, and use a convergence tolerance two orders tighter.
- **Problems 3 and 5 are elastic.** Both files set the plasticity specification to none and carry
  no strength values under it, so the rock cannot yield and only the joints can.
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
