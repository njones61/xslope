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

**The short version for practice.** The factor of safety a jointed strength reduction reports is
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

The element's own law — the peak limit, the residual drop and the opening a dilating joint produces
per unit of slip — is checked against its closed forms by `test/joint_element_check.py`.

## Inputs

### The joints worksheet

One row per line: a `Label`, the endpoints `x1, y1, x2, y2`, and the properties above — `c`, `phi`,
`c_res`, `phi_res`, `dil`, `t_cut`, `kn`, `ks` and `Jred`. Only `phi` and the endpoints are
required; every other column is blank for the ordinary case. The columns and their units are
documented with the rest of the template under
[Worksheet: joints](../usage/input_template.md#worksheet-joints). In Studio the same rows are
edited in the joints editor, as a table or one line at a time with the section drawn beside it.

![The joints editor in Studio: the lines of a model, one selected, with its properties and the section beside it](../studio/images/editing_joints_editor.png){width=900}

### Generating a network of joints

A jointed rock mass has hundreds of joints, and nobody enters them one line at a time. Instead you
describe the pattern and the program writes the lines for you:

- **Parallel set.** One family of joints at a given dip and spacing, such as bedding planes every
  2 m dipping 35° out of the face.
- **Cross-jointed.** Two such families crossing, such as bedding plus a steeper joint set, which
  cuts the rock into blocks.
- **Voronoi.** A random pattern of blocks of a given size, for a rock mass with no preferred joint
  direction.

You also say where the pattern applies: the whole section, one material, or a region you draw
yourself. A drawn region is a polygon on the polygon sheet with its **Type** set to `joints` (see
[joint regions](../usage/input_template.md#joint-regions) on the template page), and it can be
limited to a band of elevations. Outside the region no joints are written.

In Studio this is the [Build network](../studio/editing.md#build-network) button on the joints
editor: choose the pattern, enter its numbers, pick the region, and the canvas previews the joints
before anything is written. In a script the same three patterns are the functions `parallel_set`,
`cross_jointed` and `voronoi` in `xslope.joints`.

What gets written is ordinary rows on the joints worksheet, labeled `set-01`, `set-02` and so on.
The lines are the model; the pattern that made them is not remembered. To change a network, remove
the set and build another.

![The Build network dialog: the kind of network, its parameters, and the region it is clipped to](../studio/images/editing_joint_network_dialog.png){width=760}

### A reinforcement line as a joint

Setting `Joint = Yes` on a reinforcement line splits the mesh along it and puts the sheet between
two interfaces, one against the soil above and one against the soil below, so the sheet's two
interfaces act in series where a joints-sheet line carries one contact.
[Ends, ties and the bar](reinforcement.md#ends-ties-and-the-bar) covers how such a sheet is
anchored, and what the bar does and does not carry once the interfaces carry the grip.

## Running a Jointed Model

### The sweep budget

A joint reaches equilibrium by **growing slip**, a little per sweep, so a jointed model settles over
tens of thousands of viscoplastic sweeps where a bonded one settles over hundreds. A
strength-reduction trial that runs out of sweeps is recorded undecided, the bracket reads that as
not standing, and the factor of safety comes out low — a reading of the budget rather than of the
slope.

Raise `max_iterations` on `solve_fem()` and `solve_ssrm()`, or the sweep limit in Studio's Run FEM
dialog, from its default of 12,000 to **100,000**. The model checks warn below that number.
Allowing more costs almost nothing, because a trial that settles stops early. A six-course block
wall with three geogrid layers, run at the default, reads one trial of its bracket as standing where
more sweeps show it still moving, so its factor of safety comes out too high; at 36,000 sweeps it
comes out too low; at 100,000 it is settled.

### How a trial is decided

On a jointed model almost all of the out-of-balance force sits on the joints, and a pair of faces
at its slip limit alternates between slipping and sticking as the material around it breathes. What
that leaves is a steady back-and-forth the sweeps never damp out, so its average never falls under
a force tolerance no matter how many sweeps are allowed. A jointed trial that neither converges nor
runs away is therefore read from the **joints** rather than from the displacement field: slip still
growing at an undiminishing rate while the slope keeps moving is a slope failing on its joints, and
slip and movement both stopped with the soil in equilibrium is a slope standing behind that
back-and-forth. [The joint verdict](overview.md#the-joint-verdict) gives both readings and their
thresholds.

The standing reading counts only under the default `hybrid`
[failure criterion](overview.md#ssrm-failure-criteria); under `non_convergence` a trial that has
stopped moving but never converged is still a non-converged trial, and the search reads it as
failed.

A jointed trial is also handed to the [Newton
corrector](overview.md#finishing-a-trial-with-the-newton-corrector), which starts from the state
the sweeps have reached — the slip, the opening and the residual strength each joint has arrived at
— and looks for equilibrium at that strength directly. Where it finds one, the sweeps are
restarted from it to confirm that the slope stays there (the
[hold test](overview.md#finishing-a-trial-with-the-newton-corrector)), and the trial stands,
usually within a few hundred more sweeps. Where it does not, the trial is left as the sweeps read
it: not finding an equilibrium this way is not evidence that none exists.

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

## References

Goodman, R.E., Taylor, R.L., & Brekke, T.L. (1968). A model for the mechanics of jointed rock.
*Journal of the Soil Mechanics and Foundations Division*, 94(SM3), 637-659.

Schellekens, J.C.J., & de Borst, R. (1993). On the numerical integration of interface elements.
*International Journal for Numerical Methods in Engineering*, 36(1), 43-66.
