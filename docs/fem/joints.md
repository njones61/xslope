# Joints and Interface Elements

A **joint** is a surface two bodies meet on and can slide along: a rock joint, a bedding plane, the
contact between one facing block and the next, the back of a retaining wall against the soil it
holds, a geosynthetic sheet the fill above it can slide on. XSLOPE models one as a **joint line** —
a line the finite element mesh is split along, with an interface element carrying the normal and
shear stress between the two faces (mechanics texts call these the tractions). The faces can then
slide on each other, part, and come back together, which a single bonded mesh cannot do.

![A rock cut in bedded rock: a bedding set dipping out of the face and daylighting in it, a steeper cross-joint set and two release joints, which between them cut the mass into blocks, with a talus of shed blocks at the toe](images/joints_rock_slope.png){width=900}

Joint lines come from two places. A line on the **joints** worksheet is a joint and nothing else —
the rock joint, the bedding plane, the block contact above. A line on the **reinforce** worksheet
whose `Joint` column reads **Yes** is a reinforcement sheet that is *also* a slip surface: the mesh
splits along it, and the sheet keeps its own tension. Everything on this page applies to both
kinds; [Two ways to represent a sheet](reinforcement.md#two-ways-to-represent-a-sheet) covers what
is particular to the reinforced one.

A joint line is finite element geometry. The limit equilibrium engines do not read the joints
worksheet at all.

## When a Model Needs a Joint

A joint belongs where failure **runs along** a surface rather than cutting across it.

**Rock.** A jointed rock mass fails on its discontinuities, not through intact rock: block and
flexural toppling of a columnar face, a plane failure on a bedding plane that daylights, a step-path
surface running from one joint to the next through short rock bridges, a slab plowing into the
block below it. The rock itself is often modeled as elastic, so every way such a slope can fail is
a movement on its joints.

**Blocks and walls.** A segmental block wall is a stack of separate blocks. Modeled as one solid
it cannot move the way a stack does. With a joint under the bottom block, a joint on the back of the
wall against the fill, and a joint between every course of blocks, the blocks can slide on one
another, separate, and tip. The same applies to a gravity or gabion wall, and to the contact between
a structure and the ground it stands on.

**Sheets that are the slip surface.** A base geotextile under an embankment on soft clay, a smooth
geomembrane or liner, a wrap-around geotextile wall (no facing blocks: each lift of fill sits on its
sheet, which is folded back over the face and buried under the next lift): the fill slides *on* the
sheet at the interface friction. Those are reinforce-sheet lines with `Joint = Yes`; which of the two a sheet
needs is set out under
[Choosing a bonded bar or a joint](reinforcement.md#bonded-bar-or-joint).

![A segmental block wall: the blocks stand on a joint under the base, a joint on the back face against the fill and a joint between every course of blocks, with geogrid layers tied into the blocks and running back through the reinforced fill](images/joints_block_wall.png){width=900}

## The Split Mesh

Where a joint line runs, the mesher gives every node on it **one copy for each piece of material
around it**. A node in the middle of a line has material above and below, so it becomes two nodes at
the same point, one on each side. Where two joint lines cross, the node sits at the middle of four
wedges of material and becomes four; where one joint line ends on another, three wedges and three
copies. Each element keeps the copy on its own side of every line through the point, so the pieces
are free to move relative to one another, and the interface elements between the copies hold them
together.

![The mesh around one node on a joint line, at a crossing of two joint lines and at a termination, with the pieces of material drawn pulled apart along the joint traces and one node copy standing in each](images/joint_mesh_split.png){width=760}

The split is invisible on a mesh plot, because the copies stand at one point; a jointed line is
drawn in a style of its own, so it can be told from a line the mesh is not split along.

### Where a joint ends on another joint

Joints may meet: a release joint can end on a bedding plane, and a block's base can start partway
along its neighbor's side. Where two joint lines are meant to meet but the coordinates miss by a
hair, the end is moved onto the other line before the mesh is built, so nothing has to be entered
to more decimals than the drawing gives.

A joint line may not do two things, and the model checks reject both, naming the line: run along
another joint line for part of its length, and run along the outer boundary of the section, where there is
rock on one side only.

## The Interface Element

Each interface element is a zero-thickness joint (Goodman, Taylor & Brekke, 1968) spanning the node
pairs of one mesh edge: three pairs on a quadratic mesh, two on a linear one. It has no area and no
thickness. Its state is the relative displacement of the two faces, resolved along the element's own
chord into a tangential component $\Delta_t$ (sliding) and a normal component $\Delta_n$ (closing,
compression positive).

![The interface element between two elements of the mesh: its three node pairs, the normal and shear stiffness carried at one pair, and the relative movement of its two faces resolved into sliding and closing](images/joint_element.png){width=920}

While the joint is intact the normal and shear stress on it, $t_n$ and $t_s$, grow with that
movement:

$$t_n = k_n \Delta_n, \qquad t_s = k_s \Delta_t,$$

and they are integrated at the element's own nodes rather than at Gauss points, which keeps the
stress along a stiff interface free of the oscillation Gauss quadrature produces there
(Schellekens & de Borst, 1993).

### Stiffness

A joint has two stiffnesses. The normal stiffness `kn` acts across the joint: how much the two
faces compress together under a given normal stress. The shear stiffness `ks` acts along it: how
far the faces shift past each other under a given shear stress before the joint slips. They are not
soil properties. They exist so that a joint that has not slipped or opened behaves as if it were
not there, its two faces moving together like the material on either side, whether that is rock on
rock or fill on a geosynthetic sheet; they need to be just large enough for that and no larger. Enter them in the `kn` and `ks` columns of the joints
sheet, or of the reinforce sheet for a reinforcement line that is a slip surface. Left blank,
they are derived from the softer of the two materials the joint runs between: the normal stiffness
is that material's Young's modulus divided by a notional joint thickness, the shear stiffness its
shear modulus divided by the same thickness, and the thickness is one tenth of the joint element's
length,

$$k_n = \frac{E}{d_v}, \qquad k_s = \frac{G}{d_v}, \qquad d_v = 0.1\,L_{elem}.$$

Where a joint line crosses a material boundary, the softer material can differ from one end of the
line to the other, and so can the derived stiffness.

How much the factor of safety depends on them depends on the mechanism. For a block sliding on a
plane it changes by about 1.5% over a hundredfold range of stiffness. For a stack of columns that
lean on one another, the joints' give is part of how the load passes from column to column, and
the factor moved by about 4% over the same range on a toppling problem of the RS2 joint
verification set. What they always change is run time. A joint ten times stiffer takes about ten times as many iterations to
settle, so a model that states stiffnesses well above the defaults needs its iteration budget
raised to match.

### Strength

The shear stress on the joint is limited by Mohr-Coulomb,

$$|t_s| \le c_j + t_n \tan\phi_j ,$$

and the normal stress by a **tension cutoff** `t_cut`. Past the shear limit the joint **slips**:
the shear stress stays at the limit while the tangential offset grows, which is perfectly plastic
slip. Past the tension cutoff it **opens**: both stresses and both stiffnesses go to zero, and it
carries nothing until the two faces come back into contact, when it closes again and carries stress
as before.

![The joint's strength envelope: the Coulomb limit on the shear stress in either direction, the cohesion intercept, the tension cutoff at which the joint opens, and the residual envelope it drops to once it has slipped](images/joint_envelope.png){width=800}

A joints-sheet line states `c`, `phi` and `t_cut` in columns of its own. A jointed reinforcement
line takes its `Adhesion` and `Delta` as the joint's cohesion and friction angle, and its tension
cutoff is zero: a
soil-geosynthetic contact carries no tension, where a rock joint may hold a little.

### Residual strength

A rough joint surface shears off its bumps and ridges (its asperities, in rock-mechanics terms) the
first time it slips and does not grow them back, so a joint that has slipped may be weaker than one that has not. `c_res` and
`phi_res` state what it keeps. The drop is **immediate** — the limit falls from
$c + t_n \tan\phi$ to $c_{res} + t_n \tan\phi_{res}$ on the sweep after the one that first found the
pair at its limit — and **permanent**: the pair stays on the residual branch for the rest of the
run, even where it later closes or unloads.

![Shear stress against slip: elastic at k_s to the peak limit, then an immediate and permanent drop to the residual limit](images/joint_residual.png){width=820}

Blank residuals mean no residual branch, and the peak carries throughout. Neither residual may
exceed its peak.

### Dilation

A rough joint rides up over the bumps on its surface as it slides. `dil` is that angle. As the joint slips it also
opens, in proportion to the slip: an opening of $\Delta_n = \Delta_t \tan(\text{dil})$ for a slip
of $\Delta_t$, and that opening is permanent. Where the joint is held shut, the
opening it cannot make is taken up as extra compression across it, so the normal stress on the
joint grows while it slides.

![Dilation: the block rides up over the bumps on the joint surface as it slides, opening the joint; held closed it builds normal stress, free to lift it rises instead](images/joint_dilation.png){width=950}

Where the material around the joint holds it closed, that is **dilatant hardening**: sliding builds
normal stress and with it shear strength. Where the sliding block is free to lift, it lifts instead,
and the normal stress stays at whatever equilibrium with the block's weight requires. The joint opens whichever direction it slides, since the opening depends on how far it slips and
not on which way. Blank is zero.

**There is no limit on the opening.** A real joint stops dilating once it has climbed or sheared
off its bumps, so its opening is limited to about the bump height. The law here has no such limit
and keeps opening at the stated angle for as long as the joint slides. That simplification is
harmless at the slips a standing model reaches, which are millimeters, and a joint that slides far
enough for the limit to matter is on a slope that has already failed.

### Strength reduction

A strength reduction divides $c_j$ and $\tan\phi_j$ by the trial factor along with the soil's, on
both the peak and the residual branch, on every jointed line. A line that sets `Jred = No` is
left at full strength while everything else is weakened. Use that for a contact whose strength is
known and not in question, such as a wall's base on a prepared bedding or a liner with a measured
interface friction; a natural joint whose strength is uncertain is part of the margin the search is
looking for, and should be reduced with the rock. The stiffnesses `kn` and `ks` are not
strengths, so they are not reduced, any more than the bar's are.

The element's own law — the peak limit, the residual drop and the opening a dilating joint produces
per unit of slip — is checked against its closed forms by `test/joint_element_check.py`.

## What the Method Can and Cannot Model

The element is a small-strain interface between **fixed node pairs**. A pair carries compression
across the joint and Coulomb shear along it, opens when the normal stress goes into tension, slips
when the shear reaches its limit, and re-closes when the two faces meet again — and it keeps the
partner it started with for the whole solve.

![What a fixed node pair carries — sliding, opening, rocking and re-closing on the same contact — and what it does not: a contact that migrates, and a new contact between faces that were never paired](images/joint_reach.png){width=950}

So blocks slide on their contacts, open at them, tip about them and settle back onto them, and
where the blocks keep the same contacts throughout, the answer is the one rigid-block statics gives.
What the element cannot do follows from the same construction: a corner cannot travel along a face,
so a contact cannot move or shorten as a block moves; no new contact forms between two faces that
were not paired to begin with; rotations have to stay small, because the element is written on the
undeformed geometry; there is no excavation stage; and the joint obeys Mohr-Coulomb with a residual
strength and a dilation angle, rather than a hyperbolic or work-softening law.

That covers the mechanisms a jointed slope usually fails by — block toppling, flexural toppling,
plane failure, step-path failure through rock bridges, plowing slabs, a mass cut into many small
blocks, a block wall sliding and tipping on its courses, an embankment sliding on its base sheet.

The factor of safety a jointed strength reduction reports is
valid: it is the factor by which the joint strengths can be reduced before the mass starts to move,
which is the same quantity the closed forms and a distinct element strength reduction report. On
the problems that have a closed form the method reproduces them — a slab on a daylighting bedding
plane to within a bisection step of tan φ / tan β, the block toppling stacks of Goodman & Bray to
within a few percent ([Tutorial FEM-5](../tutorials/fem05_rock_slope_joints.md) works the first by
hand). What the method does not answer is anything about the movement after that point: how far
the blocks travel, where they come to rest, or what they strike on the way. The deformed section
at failure shows the mechanism that starts, not the run-out. Those questions need a distinct
element program such as UDEC, in which blocks separate, rotate through large angles and make new
contacts as they go, or a rockfall program that follows each block down the slope.

A jointed model needs far more iterations than one without joints. A joint reaches equilibrium
by slipping a little at a time, so a jointed model settles over tens or hundreds of thousands of
iterations where a model without joints settles in hundreds, and a trial cut off before it has
settled is classified as a failure or as undecided when it would have come to rest. The closing
summary in the Log reports when that has happened: it gives the factor of safety as "at least" some
value, or states that the answer depends on the iteration limit. Either way the run needs more iterations;
[Running a Jointed Model](#how-many-iterations-a-jointed-model-needs) says how to set the limit
and how to continue a run that stopped short.
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) shows the effect on a geogrid wall:
at a limit of 100,000 the search reports a factor of safety of at least 1.56, and with the
limit raised to a million the wall gives way at 2.0, after 30 cm of movement.

## Inputs

A joint enters a model as a row on the **joints** worksheet of the
[input template](../usage/input_template.md#worksheet-joints), one line per row with its
endpoints and properties. You can type those rows into the workbook in Excel, edit them in
Studio's joints editor, or have XSLOPE generate a whole network of them from a pattern you
describe. A reinforcement line can also be declared a joint on the reinforce worksheet, so that
the soil slides on the sheet instead of bonding to it.

### The joints worksheet

One row per line: a `Label`, the endpoints `x1, y1, x2, y2`, and the properties above — `c`, `phi`,
`c_res`, `phi_res`, `dil`, `t_cut`, `kn`, `ks` and `Jred`. Only `phi` and the endpoints are
required; every other column is blank for the ordinary case. The columns and their units are
documented with the rest of the template under
[Worksheet: joints](../usage/input_template.md#worksheet-joints). In Studio the same rows are
edited in the joints editor, as a table or one line at a time with the section drawn beside it.

![The joints editor in Studio: the lines of a model, one selected, with its properties and the section beside it](../studio/images/editing_joints_editor.png){width=1240}

### Generating a network of joints

A jointed rock mass may have hundreds of joints. To create a network of joints in Studio, press
[Build network](../studio/editing.md#build-network) on the joints editor, choose the pattern,
enter its numbers and pick the region; the canvas previews the joints before anything is written.

![The Build network dialog: the kind of network, its parameters, and the region it is clipped to](../studio/images/editing_joint_network_dialog.png){width=1180}

A joint network is built from one of three kinds of pattern, chosen in the dialog's **Kind**
selector:

- **Parallel set.** One family of joints at a given dip and spacing, such as bedding planes every
  2 m dipping 35° out of the face.
- **Cross-jointed (two sets).** Two such families crossing, such as bedding plus a steeper joint
  set, which cuts the rock into blocks.
- **Voronoi (blocky mass).** A random pattern of blocks of a given size, for a rock mass with no
  preferred joint direction.

Each kind has its own numbers, which appear below the selector. **Within** says where the set
exists: the whole section, one material, or a region you draw. A drawn region is a polygon on
the polygon sheet with its **Type** set to `joints` (see
[joint regions](../usage/input_template.md#joint-regions) on the template page), and
**Elevation band** limits the set to a range of elevations. No joints are written outside the
region.

Enter a name for the set in the **Name** field, or keep the one the dialog offers (`set1`,
`set2`, …). The generator writes one row per joint on the joints worksheet and labels each row
with that name and a number, so a set called `bed` gives `bed-01`, `bed-02` and so on. To
remove the network later, select any of its rows in the joints editor and press **Remove set**.
XSLOPE keeps the lines, not the pattern that made them, so to change a network you remove the
set and build a new one.

To generate a joint network from a Python script, call `parallel_set`, `cross_jointed` or
`voronoi` from `xslope.joints`; they write the same rows.

### A reinforcement line as a joint

Joints are also used with reinforcement lines. To make a reinforcement line a slip surface,
set its **Joint** column to `Yes` on the reinforce worksheet, or the Joint field in Studio's
reinforcement editor. The mesh then
splits along the line and the sheet sits between two interfaces, one against the soil above and
one against the soil below, whose strength comes from the line's own **Adhesion** and **Delta**.
The sheet's two interfaces act in series where a line on the joints worksheet has a single
contact. How such a sheet is anchored, and what the bar carries once the interfaces carry the
grip, is on the reinforcement page under
[Ends, ties and the bar](reinforcement.md#ends-ties-and-the-bar); when a sheet should be a
joint and when it should be a bonded bar is worked through in
[Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md#part-3-when-a-sheet-is-a-slip-surface-and-when-it-is-bonded).

A reinforcement line that ends on a joint line, or crosses one, must be a joint itself or stop
short of the joint, and a pile, which cannot be a joint, has to stop short. The split gives every
node on the joint one copy for each side, and a bonded bar attaches to the soil's own nodes, so a
bar standing on one of those nodes has no single side to attach to. A geogrid running back from a
block facing ends on the back-face joint, which is why each sheet in such a wall is set to
`Joint = Yes`. The model checks name the line and the joint before the mesh is built, and count
an end that misses the joint by less than the mesher tells apart as being on it.

## Running a Jointed Model

A jointed model is run the way any other finite element model is: build the mesh, open
**Run → Run FEM…**, choose the strength reduction, and press **Run**. The differences are in
how many iterations the run needs, how each trial is judged, and what the results show.

### How many iterations a jointed model needs

A joint reaches equilibrium by slipping a little at a time, so a jointed model settles over tens
or hundreds of thousands of iterations where a model without joints settles over hundreds. A
strength-reduction trial that runs out of iterations before it has settled is recorded undecided.
When that trial is still the top of the final bracket, the run found no failure and reports the
factor of safety as "at least" the highest strength the slope came to rest at, since no trial
above it was shown to fail.

Set **Max iterations per trial** in Studio's Run FEM dialog to at least **100,000** for a jointed
model, in place of the 12,000 it opens with (`max_iterations` on `solve_fem()` and
`solve_ssrm()` in Python). Allowing more costs almost nothing, because a trial that settles stops
as soon as it has.

When a run stops short, the closing summary in the Log reports it: it gives the factor of safety
as "at least" some value, or states that the answer depends on the iteration limit. There are two ways
to give the run more iterations. On the FEM · Results toolbar, press
**Continue with a higher limit…**, enter a new **Max iterations per trial** (the dialog offers
five times the limit the run stopped at), and the search picks up where its trials stopped,
keeping every trial it has already decided. Or open **Run → Run FEM…**, raise **Max iterations
per trial** and the **Iteration ceiling** with it, and run again from the start; the answer is
the same. The geogrid wall of [Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) reports
a factor of safety of at least 1.56 at 100,000 iterations and gives way at 2.0 with the limit at
a million; the tutorial shows both runs.

### How a trial is decided

<span id="jointed-model-solver-policy"></span>

A jointed trial can converge under the same [force and displacement checks](solver.md#convergence-criterion)
as a model without joints. A contact at its slip limit can also keep switching between slipping and holding even
after the slope has settled, leaving a small imbalance that does not disappear with more iterations.
A trial that neither converges nor runs away is therefore judged by what the slope is doing:
whether the joints keep slipping and the ground keeps moving, or both have come to rest.

The imbalance can be almost entirely on the joints. On one jointed rock slope, the force check's
normalized soil reading was $3\times10^{-8}$ while its joint reading was
$3\times10^{-2}$, six orders of magnitude larger on the same mesh in the same solve. The
force test reads the change in the added body load; on a jointed model that change comes
mainly from the shear force the interfaces could not hold on the last iteration. The
switching took tens of iterations, rather than the two-iteration oscillation the ordinary
ten-iteration force average removes. Even after the slope settled, the joint reading varied
between $1.2\times10^{-2}$ and $4.7\times10^{-2}$, with its mean changing only 0.04% over the
last half of the solve. More iterations do not remove an imbalance of that kind.

Use the default `hybrid` [failure criterion](solver.md#ssrm-failure-criteria) for a jointed model:
it accepts a trial shown to be standing even when the force tolerance is not met. The
`non_convergence` criterion requires convergence and counts an otherwise settled jointed trial
as failed, except where the solver [leaves the trial undecided](solver.md#creep-trend); it is
kept for reproducing published results obtained that way.

#### Slipping or standing {#the-joint-verdict}

A jointed trial is judged from the slip along the interfaces as well as the movement of the
ground. Movement is measured in **elastic displacements**: the maximum displacement the same
model would make under the same loads if its response were entirely elastic. Slip is measured
against the total slip already accumulated, so neither reading depends on the units, stiffness
or mesh size.

The displacement field alone can miss a mechanism. A measured rock slope gaining 0.4 elastic
displacements every 25,000 iterations was classified `AMBIGUOUS` because its maximum movement
had not passed 1.5 times an elastic response set by loading stiff rock. Its joint slip showed
that the slope was moving.

For sliding, the [trend reading](solver.md#creep-trend) uses five equal blocks spanning nominally
half the original **Max iterations per trial** allowance. A failing reading needs at least
0.02 elastic displacements gained over the window. The rate ratio must be at least 0.9,
averaged across ground-movement blocks with a positive gain in each, or taken between the two halves of the
slip window; the slip reading also needs a gain of at least 2% of the accumulated slip.
A before-limit decision waits until both 25,000 iterations have run and the trial
has entered the last tenth of its current allowance; the trend is also read at the limit.

That wait matters on a slow trial. One trial reached equilibrium at 185,381 iterations, but at
50,000 had gained 30% of its slip and 1.76 elastic displacements, was moving at an accelerating
rate and stood at 4.3 times its elastic response. The slip rate also separates two trials that
look alike early on: at 25,000 iterations a trial that later converged after 203,000 had gained
16.2% of its slip and 0.375 elastic displacements, while one that never converged had gained
18.7% and 0.488. On the first, the rate ratio fell from 0.64 at 25,000 to 0.44 at 100,000 and
0.21 at 200,000; on a real mechanism it was 1.0002. A trial that is still slowing is left to
the iteration-limit rules, including their corrector attempts and extensions.

For a settled trial, the readings cover the trailing half of its sampled history. The slip
must gain no more than 0.01% of itself, and the ground no more than $10^{-4}$ elastic displacements.
The imbalance at the free nodes carrying no joint must not exceed the solve's force tolerance
(`force_tol`) throughout that window. The joint imbalance must also have stopped falling: its mean in the second half of the
window must be at least 0.85 of the mean in the first half. A trial whose joint imbalance is
still falling can yet converge; in one measured trial it fell 35% per window after the slip,
ground and soil were already at rest, and convergence followed nine thousand iterations later.

The standing result is recorded as `JOINT_SETTLED`, with `exit_reason = 'joint_settled'`;
`converged` remains `False`, because the force tolerance was not met. It moves the hybrid
search's bracket upward, just as a converged trial does. Sliding is recorded as `FAILED`, with
`exit_reason = 'not_slowing'`, and moves the bracket downward. If neither reading decides the
trial, the displacement classifier and the iteration-limit rules apply. These joint readings
are not taken on a model without joints or on a trial already decided by convergence or runaway.

#### Contacts that cycle

A contact can repeat a cycle rather than settle in one state. Where a joint with cohesion but
no tensile strength holds a contact at zero normal stress, the contact switches every few
iterations between closed, carrying shear, and open, carrying nothing, while the rest of the
slope stays still. Held closed, the rest of the slope balances with the contact in tension;
held open, it balances with the faces overlapping. Neither satisfies the joint law.

At the hard iteration ceiling, a trial of this kind counts as standing only if the whole set
of contact states repeats exactly with a period no longer than 64 iterations over the last 256.
Between one and four contacts may change state within a period, and the field must return to
itself after each period to within $5\times10^{-8}$ elastic displacements per iteration. That
movement bound is about ten times the measured standing trial's $5\times10^{-9}$ and about a
tenth of the slowest failing trial with cycling contacts, $5.4\times10^{-7}$.

The trial is recorded as `JOINT_SETTLED`, with `exit_reason = 'joint_settled'` and a
`contact_cycle` reading. The trial line, closing summary and report state that force balance
was not met: "stands: the only movement is 3 contacts cycling (period 8 iterations), no net
movement; force balance not met", with each contact named by its line and location.

#### Finishing a slow jointed trial

While a trial is slowing down, the run periodically takes a shortcut: from the state already
reached, the [Newton corrector](solver.md#finishing-a-trial-with-the-newton-corrector) solves
directly for a state in which the forces balance. Its interface starts with the slip, permanent
dilational opening, residual strength and opening history the ordinary iteration has reached.
A pair that first reaches its limit during that solve moves onto its residual branch; one
that slips further adds the corresponding dilation. The corrector does not start with a
pristine joint and discard the path the trial has followed.

The corrector also carries the change in shear stress caused by normal movement,
$\partial t_s/\partial\Delta_n = \pm k_n\tan\phi_j$, in its stiffness matrix. That term makes
the matrix non-symmetric, and it is handled by the Newton solve's general factorization.
The ordinary iteration keeps the elastic stiffness instead: a slipping pair has zero
tangential stiffness but the matrix still carries its full $k_s$. That mismatch makes a trial
slow. On one rock-toppling trial at the edge of the final bracket, the error shrank by a
factor of only 0.99947 each iteration, taking 13,000 iterations to fall by a factor of a thousand.
A failing block can likewise slide at a few billionths of an elastic displacement per iteration
after both residuals have gone flat within a few hundred iterations, taking a long time to
reach the runaway reading of eight elastic displacements.

A state the corrector finds must pass its force, yield and displacement checks. It then gets
a **hold test**: the ordinary iteration continues from that state, with the same joint history,
for up to 3,000 more iterations, moving no more than a hundredth of an elastic displacement.
The state is accepted only if that continuation counts as standing. An equilibrium on one set
of contacts can be a state the slope leaves as soon as it is allowed to move; the hold test
rejects that state. Models without joints do not run it. If the corrector or hold test finds
nothing acceptable, the ordinary iteration carries on unchanged.

One rock-toppling bracket edge that took 185,381 ordinary iterations was certified from
300 iterations, with a force imbalance of $3\times10^{-11}$ and no integration point outside
its yield surface. The whole bracket closed at the same factor as with the ordinary iteration,
in a thirteenth of the time.

#### Solver settings for joints

The default `fem_solver='auto'` uses the ordinary iteration with the corrector. On
[joint interfaces](reinforcement.md#joints-without-reinforcement), `joint_newton=False` on
`solve_fem()` or `solve_ssrm()` turns the corrector off; `fem.JOINT_NEWTON_ON = False` does the
same for a whole process. `fem_solver='viscoplastic'` turns it off on every model. Trials then
end on the slip and movement readings above and the iteration-limit rules.

The cold-start `fem_solver='newton'` has no accumulated slip behind it. On the rock-joint
benchmarks it diverges on trials the ordinary iteration converges. That divergence is
reported as a failed trial rather than an error, but this driver cannot bracket those models.

The **interface relief**, `joint_tangent='slip'`, is off by default. It lowers the assembled
shear stiffness of slipping pairs and both stiffnesses of open ones to `joint_tangent_factor`
(default 0.01) of their elastic values, and refactorizes the matrix when that set changes.
The traction limit, slip return and equilibrium state are unchanged. The relieved iteration
runs until it settles or spends its own allowance, then turns off and passes its state to
the ordinary iteration, which decides the trial from its own slip and movement history.
The thresholds above are calibrated on that history, not on the relieved iteration.
Relief cut one benchmark's failing trial from 196,201 to 20,991 iterations, but also moved a
bracket the corrector reproduced exactly, so it remains an option. It has no effect without joints.

#### Acceleration

Jointed trials are accelerated by default (`accelerate=None`); `accelerate=False` uses the
ordinary iteration, and `True` also accelerates models without joints. With the corrector on,
acceleration starts after its last checkpoint. Each step is multiplied by a factor between
1 and 50 read from the last two steps (Irons and Tuck's extrapolation), with the plastic
strains and joint slip scaled with it. A step that would change a pair's open or slipping
state, move it onto residual strength or add dilation keeps its ordinary length.

Acceleration changes the number of iterations to the balanced state, not the state itself.
Over the 32 jointed verification rows the answers were the same and the set ran about a
fifth faster. The [K0 in-situ solve](solver.md#in-situ-equilibration) and the hold test always
use the ordinary iteration. The opening lines in the Log show whether acceleration was on.

### What the results show

For jointed models, the results panels include features that show what the joints are doing: on
the shear strain panel every joint is drawn and colored by how far it has slipped, the deformed
section is drawn as the blocks the joints cut it into, and **1D Details…** lists every jointed
line with the stresses along it.

On the shear strain panel every joint is a thin line along its trace, colored on a green scale
titled *Joint slip* by how far its two faces have slid. A closed joint, one that has not slipped,
is a gray hairline. A stretch that has **opened** is drawn as two thin gray lines with a white gap
between them, since opening is a condition rather than an amount; no other line carries white. A
jointed reinforcement sheet is drawn as its two faces, one on each side of the bar. A key names the
three states in a corner of the panel that the section does not reach, or below the panel when the
section reaches every corner, and a model where no joint slipped carries no slip colorbar. On a
model whose soil or rock can yield, the joints are drawn over the strain field and the panel keeps
its strain title and colorbar, with the slip colorbar beside it; on a model whose every material is
elastic there is no strain to draw, and the panel is the joints' own, titled *Joint slip*. The
lines are kept thin because a generated network puts hundreds of them over the field.

![Joint slip on a toppling stack at its critical factor: each joint colored by how far its faces have slid, gray where it has not slipped, drawn as two lines with a gap between them where it has opened](../tutorials/images/fem05_joint_slip_topple.png){width=900}

**Show joints** is set separately for the shear strain panel and the deformed mesh panel, and
is on for both by default on a jointed model. On the **Deformed mesh** panel, with it on, the
deformed section is drawn as the **blocks** the joints cut it into: the section under a faint
tint, one per material, the joints drawn as on the shear strain panel, and the deformed outline as a
dark line against the dashed undeformed outline. Two things differ from the shear strain panel: a
closed joint keeps the full joint-face weight, so the thick gray joints outline the blocks, and the
slipping joints are drawn at that same weight. The panel carries its own slip colorbar and key. A block is a piece of the mesh that moves as one body; the mesh split gives
the two sides of a joint their own nodes, so the pieces between joints are found by following
which elements still share nodes. A jointed failure looks like this, blocks moving as bodies with
all of the movement taken up at the joints, and a joint that has slipped or opened shows as two
lines that no longer lie on each other. The element grid is drawn under the blocks, in light
gray, when the model has eight or fewer jointed lines and left out when it has more, so a
generated network stays readable; **Element edges** in the display panel turns it on or off
either way. **Color by block** sets the fill: off (the default), one tint per material; on, each
block under its own tint so the bodies can be told apart. The toppling stack of Tutorial FEM-5
with the grid on and off:

![The Deformed mesh panel of the toppling stack at its critical factor with Element edges on and Color by block off: the grid under the blocks, one tint for the one material](../tutorials/images/fem05_fem_blocks_topple.png){width=900}

![The same panel with Element edges off and Color by block on: the four columns alternating two tints and the mass a third](../tutorials/images/fem05_fem_blocks_topple_colored.png){width=900}

With Show joints off the panel draws the ordinary deformed mesh. The **Displacement vectors**
panel is the same as on any other model.

**1D Details…** lists every jointed line, from the joints sheet and from the reinforce sheet
alike, under a *Joints* heading, and draws panels along each line: the normal stress on the
interface, the shear stress with its Mohr-Coulomb limit drawn beside it, and the slip. A
reinforcement line that is also a joint gets a fourth panel above those, the bar's tension against
its capacity; a joints-sheet line has no bar and shows the three. Where the shear stress meets its
limit is where the interface is slipping. A
generated report gives the same information as a table: the share of each line's length standing at
its limit, and the largest offset the two faces reached.

![1D Details for a reinforcement line built as a slip surface, a geogrid layer in a block wall. Because this line carries a bar, it has the four panels: the bar's tension against its capacity, then the normal stress on the interface, the shear stress against its Mohr-Coulomb limit, and the slip along the line. A joints-sheet line has no bar and shows the last three only](../tutorials/images/fem03_1d_details.png){width=900}

**Shear traction zigzag on a gripping interface.** Along a stretch of interface
that is not slipping, the shear traction plotted station by station can alternate high and low:
the end stations of each interface element read low and the middle station reads high. The
stations are the element's three node pairs, and each pair's traction is the force passed
through that node divided by the length of interface the node stands for, one sixth of the
element at each end and two thirds in the middle. Where the interface grips, its relative
displacement is tiny, so the traction at a pair is set by the force the neighboring soil
elements pass through that node. Quadratic soil elements pass the forces carried through their
bodies to their mid-side nodes (the weight of a six-node triangle, for instance, is carried
entirely at its three mid-side nodes), and that is what the middle stations show. The zigzag comes
from that sharing rather than from a variation of the stress along the sheet: the element average, one sixth of
each end station plus two thirds of the middle one, is the force the element actually
transfers, it runs smoothly, and a ten times stiffer interface does not change it. It
disappears where the interface slips, because every slipping station is held at its own
Mohr-Coulomb limit, which follows the normal stress. Read a gripping stretch by its element
averages, not station by station.
{ #shear-zigzag }
## References

Goodman, R.E., Taylor, R.L., & Brekke, T.L. (1968). A model for the mechanics of jointed rock.
*Journal of the Soil Mechanics and Foundations Division*, 94(SM3), 637-659.

Schellekens, J.C.J., & de Borst, R. (1993). On the numerical integration of interface elements.
*International Journal for Numerical Methods in Engineering*, 36(1), 43-66.
