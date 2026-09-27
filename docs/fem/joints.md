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
[Bonded bar or joint?](reinforcement.md#bonded-bar-or-joint).

![A segmental block wall: the blocks stand on a joint under the base, a joint on the back face against the fill and a joint between every course of blocks, with geogrid layers tied into the blocks and running back through the reinforced fill](images/joints_block_wall.png){width=900}

## The Split Mesh

Where a joint line runs, the mesher gives every node on it **one copy for each piece of material
around it**. A node in the middle of a line has material above and below, so it becomes two nodes at
the same point, one on each side. Where two joint lines cross, the node sits at the middle of four
wedges of material and becomes four; where one joint line ends on another, three wedges and three
copies. Each element keeps the copy on its own side of every line through the point, so the pieces
are free to move relative to one another, and the interface elements between the copies are what
hold them together.

![The mesh around one node on a joint line, at a crossing of two joint lines and at a termination, with the pieces of material drawn pulled apart along the joint traces and one node copy standing in each](images/joint_mesh_split.png){width=760}

The split is invisible on a mesh plot, because the copies stand at one point; a jointed line is
drawn in a style of its own, so it can be told from a line the mesh is not split along.

### Where a joint ends on another joint

Joints may meet: a release joint can end on a bedding plane, and a block's base can start partway
along its neighbor's side. Where two joint lines are meant to meet but the coordinates miss by a
hair, the end is moved onto the other line before the mesh is built, so nothing has to be entered
to more decimals than the drawing gives.

Two things a joint line may not do, and the model checks refuse both by name: run along another
joint line for part of its length, and run along the outer boundary of the section, where there is
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
settled reads as a failure or as undecided when it would have come to rest. The closing summary
in the Log says when that has happened: it gives the factor of safety as "at least" some value,
or says the answer depends on the iteration limit. Either way the run needs more iterations;
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

When a run stops short, the closing summary in the Log says so: it gives the factor of safety
as "at least" some value, or says the answer depends on the iteration limit. There are two ways
to give the run more iterations. On the FEM · Results toolbar, press
**Continue with a higher limit…**, enter a new **Max iterations per trial** (the dialog offers
five times the limit the run stopped at), and the search picks up where its trials stopped,
keeping every trial it has already decided. Or open **Run → Run FEM…**, raise **Max iterations
per trial** and the **Iteration ceiling** with it, and run again from the start; the answer is
the same. The geogrid wall of [Tutorial FEM-3](../tutorials/fem03_block_wall_joints.md) reports
a factor of safety of at least 1.56 at 100,000 iterations and gives way at 2.0 with the limit at
a million; the tutorial shows both runs.

### How a trial is decided

On a model without joints, a trial converges when the forces come into balance everywhere. On a
jointed model they never quite do: a contact at its slip limit keeps flickering between slipping
and holding as the ground around it settles, so a small imbalance remains however long the run
goes on. A jointed trial that neither converges nor runs away is therefore judged by what the
slope is doing. If the joints are still slipping at a steady rate and the slope keeps moving, the
slope is failing on its joints. If the slip has stopped growing and the movement has stopped,
the slope is standing. [The joint verdict](overview.md#the-joint-verdict) and
[the trend reading](overview.md#creep-trend) on the overview page give the thresholds.

That standing reading counts only under the default `hybrid`
[failure criterion](overview.md#ssrm-failure-criteria). Under `non_convergence`, a trial that
has stopped moving but never converged is still counted as failed.

The run also tries to finish a slowing trial directly, at set checkpoints and whenever the
movement is dying away. Starting from the slip, the opening and the strength each joint has
reached, the [Newton corrector](overview.md#finishing-a-trial-with-the-newton-corrector) solves
for a balanced state at that strength. If it finds one, the iteration resumes from that state to
check that the slope stays there, and if it does the trial stands, usually within a few hundred
more iterations. If no balanced state is found, the trial is judged by its movement as above;
not finding one this way does not mean none exists.

### What the results show

Each interface is drawn as a thin line on the line it runs along, colored by how far its two faces
have slid, on a colorbar titled *Joint slip*: a green ramp, because the strain field under it runs
blue through white to red. A slipping span is backed by a thin white stroke so it reads over a dark
field; a joint that is not slipping is a neutral gray hairline; a stretch that has **opened** is
drawn as its two faces apart — two thin lines with a white gap between them, along the whole
stretch that parted — rather than given a color, because opening is a condition and not a quantity.
A key in the corner of the panel names the three states. A model where no joint slipped carries no
colorbar. The weight is deliberate: a generated network puts hundreds of traces over the field.

![Joint slip on a toppling stack at its critical factor: each joint colored by how far its faces have slid, gray where it has not slipped, drawn as two lines with a gap between them where it has opened](../tutorials/images/fem05_joint_slip_topple.png){width=900}

On a jointed model the displacement panel is the scaled deformed mesh rather than an arrow field,
drawn as the **blocks** the joints cut the section into — each block under a faint tint of its own,
its joint faces in the same green, the outside of the deformed mesh as a dark line against the
dashed undeformed outline. A block is a piece of the mesh that moves as one body, found by following
element adjacency: the split gives the two sides of a joint their own nodes, so they are no longer
neighbors. That is what a jointed failure looks like — blocks moving as bodies, with all of the
movement taken up at the joints, where a slipped or opened joint shows as two lines that no longer
lie on each other. An arrow field samples that at nodes and misses exactly the thing that happened.

![The deformed mesh of the same stack, drawn as blocks: each block moves as one body, the joints between them shown in green, and the dashed outline is the undeformed section](../tutorials/images/fem05_fem_blocks_topple.png){width=900}

**1D Details…** lists every jointed line, from the joints sheet and from the reinforce sheet
alike, under a *Joints* heading, and draws panels along each line: the normal stress on the
interface, the shear stress with its Mohr-Coulomb limit drawn beside it, and the slip. A
reinforcement line that is also a joint gets a fourth panel above those, the bar's tension against
its capacity; a joints-sheet line has no bar and shows the three. Where the shear stress meets its
limit is where the interface is slipping. A
generated report carries the same reading as a table: the share of each line's length standing at
its limit, and the largest offset the two faces reached.

![1D Details for a reinforcement line built as a slip surface, a geogrid layer in a block wall. Because this line carries a bar, it has the four panels: the bar's tension against its capacity, then the normal stress on the interface, the shear stress against its Mohr-Coulomb limit, and the slip along the line. A joints-sheet line has no bar and shows the last three only](../tutorials/images/fem03_1d_details.png){width=900}

**Why the shear traction can zigzag where an interface grips.** Along a stretch of interface
that is not slipping, the shear traction plotted station by station can alternate high and low:
the end stations of each interface element read low and the middle station reads high. The
stations are the element's three node pairs, and each pair's traction is the force passed
through that node divided by the length of interface the node stands for, one sixth of the
element at each end and two thirds in the middle. Where the interface grips, its relative
displacement is tiny, so the traction at a pair is set by the force the neighboring soil
elements pass through that node. Quadratic soil elements hand the forces carried through their
bodies to their mid-side nodes (the weight of a six-node triangle, for instance, is carried
entirely at its three mid-side nodes), and that is what the middle stations show. The zigzag is
that sharing, not a variation of the stress along the sheet: the element average, one sixth of
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
