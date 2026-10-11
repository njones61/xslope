---
title: "Piles in finite element analysis — XSLOPE"
description: "Piles and concrete piers as embedded Euler-Bernoulli beam elements in XSLOPE's finite element analysis: formulation, structural capacities, end restraints, the bonded pile-soil interface and the results panel."
---

# Piles and Concrete Piers in Finite Element Analysis

In the finite element analysis, a pile or concrete pier is a chain of Euler-Bernoulli beam elements embedded in the
soil mesh and sharing its nodes. Unlike the LEM, it applies no pile force. The beam carries axial force in tension or
compression, shear and bending moment as the soil deforms around it, up to the capacities `Vcap` and `Mcap` where they are given. Which members suit the FEM
is on the [Piles and Walls Overview](overview.md#lem-vs-fem), and what to enter for each is on
[Pile and Wall Types](types.md).

## The Beam Element

Each node of a beam element has three degrees of freedom: two translations, $u_x$ and $u_y$, and a rotation,
$\theta$. A reinforcement bar, by contrast, is a truss element that carries tension only and has no rotation
([Soil Reinforcement in Finite Element Analysis](../reinforcement/fem.md)). On a linear soil mesh, the beam element
has two nodes, and its stiffness in local coordinates, with the beam axis along local $x$, is

>$\mathbf{K}_{\text{local}} = \begin{bmatrix} \frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\ 0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\ -\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\ 0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L} \end{bmatrix}$

where $E$ is Young's modulus, $A$ the cross-sectional area, $I$ the moment of inertia and $L$ the element length.
A rotation matrix $\mathbf{R}$, with the direction cosines of the pile axis on the translations and the identity on
the rotations, takes it to global coordinates: $\mathbf{K}_{\text{global}} = \mathbf{R}^T \mathbf{K}_{\text{local}}
\mathbf{R}$.

A pile line is meshed like any constraint line: element edges follow it, and its beam elements are those edges
([Reinforcement and pile lines](../fem/mesh.md#reinforcement-and-pile-lines)). Each pile's $EA$ and $EI$ are divided
by its spacing $S$ before assembly, so the beam carries the row's stiffness per unit width. `Vcap` and `Mcap` are
divided by $S$ the same way.

### Three-Node Elements on a Quadratic Mesh

On a quadratic soil mesh (tri6, quad8 or quad9), the edge a beam element lies on has a midside node as well as its
two corners, and the beam element uses all three. It then has nine degrees of freedom, and its local stiffness is
9×9 ([Quadratic elements](../fem/mesh.md#quadratic-elements)). Its deflection uses quintic Hermite shape functions
$H_i$ that match a value and a slope at each of the three nodes. The bending block is
$EI \int H_i'' H_j'' \, dx$, so the element stays exactly Euler-Bernoulli and reproduces the cubic and quartic deflections of the classical beam
solutions. Axial action is the quadratic bar over the same three nodes, uncoupled from bending.

Unlike a two-node element, the three-node element gives the soil reaction along its own length. A two-node element's cubic deflection has
a zero fourth derivative, so a chain of them can report the reaction only as the step in shear from one element to
the next. The quintic's fourth derivative is linear along the element, and $EI$ times it is the reaction the element
carries. A fourth derivative magnifies small unevenness in the nodal loads, so each element's reaction is reported
as the value, at its own depth, of a least-squares line through its reaction and those of its two neighbors.

### Degrees of Freedom {#mixed-dof-system}

Nodes on a pile have three degrees of freedom and all others two, so XSLOPE gives each node its own count and an
offset into the global displacement vector. The total is $2 \times n_{\text{all nodes}} + n_{\text{pile nodes}}$. The
figure shows a vertical three-node beam element on the edge of a quadratic soil element:

![Translations and rotations of a vertical three-node pile beam](../fem/images/pile_beam_dofs.png){width=640px}

The midside node the beam shares with the soil carries a rotation, as the two ends do.

### Actions in the Beam

At each iteration, each beam element's actions come from its nodal displacements, rotated into local coordinates.
With $u_1, u_2$ the axial displacements of the end nodes, $v_1, v_2$ their transverse displacements and
$\theta_1, \theta_2$ their rotations:

>**Axial force**: $T = \dfrac{EA}{L} (u_2 - u_1)$

>**Shear force**: $V = \dfrac{12EI}{L^3}(v_1 - v_2) + \dfrac{6EI}{L^2}(\theta_1 + \theta_2)$ (two-node element; on a three-node element the shear is the quintic deflection's at the element center)

>**Bending moments** at nodes 1 and 2: $M_1, M_2 = \mathbf{K}_{\text{local}}[\text{rows 3 and 6}] \cdot \mathbf{u}_{\text{local}}$ (the two rotation rows)

## Structural Capacity Checks

Both capacities are imposed as a body-load correction, as soil yield and a bar's tension cap are. The stiffness
matrix $\mathbf{K}$ keeps the member's full elastic stiffness, and the correction $\mathbf{c}$ is the part of the
elastic action above the capacity. The solver solves $\mathbf{K}\mathbf{u} = \mathbf{F} + \mathbf{c}$, so the member
carries $\mathbf{K}\mathbf{u} - \mathbf{c}$, which equals the capacity where the capacity is reached. The pile
actions XSLOPE reports are $\mathbf{K}\mathbf{u} - \mathbf{c}$.

A moment capacity is a plastic hinge. Once the moment at an element end reaches $M_{\text{cap}}/S$, that end rotates
by a plastic rotation $\mathbf{p}$ in addition to the nodal rotation, and the element delivers

>$\mathbf{s} = \mathbf{K}_{\text{local}} (\mathbf{u}_{\text{local}} - \mathbf{p})$

with $\mathbf{K}_{\text{local}}\,\mathbf{p}$ added to the load vector, on the element's translational rows as well
as its rotational ones. $\mathbf{p}$ comes from one small linear solve on the rotational block of
$\mathbf{K}_{\text{local}}$ that puts each hinged end exactly on the capacity. Releasing one end raises the moment at
the other, so the set of hinged ends is grown until it is stable. At an interior node, where two element ends meet,
only one end takes the plastic rotation, and equilibrium puts the other on the capacity with the opposite sign. The
released shaft deflects further, and the factor of safety falls.

Where the elastic shear in an element exceeds $V_{\text{cap}}/S$, the correction

>$\Delta V = V - \text{sign}(V) \cdot V_{\text{cap}}/S$

is applied on the shear's own internal-force pattern, $+1$ and $-1$ on the transverse rows of the element's two end
nodes, rotated into global coordinates. The element then delivers exactly $V_{\text{cap}}/S$. The hinge check runs
first, and the shear is then read off the released element:

![One plastic release at a shared pile node and the ordered capacity checks](../fem/images/pile_capacity_checks.png){width=880px}

The third step caps the shear without changing the capped end moments.

Each element records whether it has yielded in shear, in bending or both, and the plastic rotation at each end node.
Strength reduction leaves the pile's $E$, $I$, $A$, `Vcap` and `Mcap` unchanged, as it does every
[structural property](../fem/overview.md#structural-elements).

## Head and Tip Fixity

The `Head` column sets the restraint at the pile's top node and the `Tip` column at its bottom node. An end can hold
its translation, its rotation, both or neither, so both columns offer the same four settings:

| Setting | Translation | Rotation | At the head | At the tip |
|---|---|---|---|---|
| **free** (default) | free | free | no structural connection at the top | floating in the soil, or resting on the model boundary |
| **pinned** | held | free | tie-rods or anchors holding the head in place | bearing on a hard stratum inside the mesh |
| **unrotated** | free | held | a cap beam tying the heads together | rarely meaningful; offered so the two lists match |
| **fixed** | held | held | cap beam and anchors | socketed into rock |

![Free, pinned, unrotated and fixed pile end restraints](../fem/images/pile_end_restraints.png){width=880px}

Red arrows mark the degrees of freedom an end leaves free, and black stop bars mark those it holds. Cai & Ugai (2000)
define the same four conditions, calling pinned "hinged".

Which tip condition is right depends on where the pile ends. A shaft that continues well below the slip surface is
held by the soil it passes through, so leave the tip free. A shaft whose bottom node lands on a fixed boundary of the
model is already pinned by that boundary: its translations are held, but it can still rotate about its toe. There,
`pinned` changes nothing, and `fixed` is the socketed case. A shaft that ends on a hard stratum inside the mesh needs
`pinned`, since the soil elements below it would otherwise let the tip move.

## Limits of the Beam Model {#applicability-continuous-walls-and-discrete-pile-rows}

The beam is a plane-strain member, continuous out of plane. A sheet-pile, diaphragm or secant wall is continuous, so
the beam represents it directly, with $S = 1$. The element is checked against closed-form beam theory, and a wall
analysis is compared with GeoStudio's in [the SIGMA/W wall benchmark](../verification/geostudio.md#sigmaw-wall). For
a row of separate piles, dividing $EA$ and $EI$ by $S$ gives only the row's average stiffness, so the beam overstates
the row's resistance. It misses the soil arching onto the piles and, at wide enough spacing, moving between them
([LEM vs FEM](overview.md#lem-vs-fem)).

The beam elements stand on the nodes of the soil element edges they lie on, the midside node included on a quadratic
mesh, so pile and soil move together at every node. The interface is perfectly bonded. The shaft cannot slip, and no
interface strength limits the load passed between pile and soil. Only the soil's strength limits it, and, for
sideways load, the pile's `Vcap` and `Mcap`. XSLOPE has no interface element along a pile shaft:
[interface elements](../fem/joints.md) apply only to joint lines and jointed reinforcement lines.

For a stabilizing pile, the soil pushes sideways on the pile, and the soil yielding around it limits the load. The
bond overstates the pile's resistance and makes the factor of safety unconservative. On the one benchmark with a
three-dimensional answer, the FEM's factor of safety is above the three-dimensional value
([LEM vs. FEM Pile Modeling](lem.md#lem-vs-fem-pile-modeling)). The bond and the spreading of the row into a wall
both raise the FEM value. That comparison does not separate the two.

For a load-bearing pile, an axial load at the head reaches the soil in proportion to the stiffness of the beam and
of the soil elements at each node, not through skin friction and end bearing. The pile cannot slip or punch through,
and the depth at which the load is transferred depends on the mesh. So, in the FEM, do not enter a load-bearing pile
near a slope on the `piles` sheet. Apply its load as described under
[Load-Bearing Piles Near a Slope](types.md#load-bearing-piles-near-a-slope).

## Inspecting the Results

The FEM results view colors pile elements by the shear they carry. The **1D Details…** panel, described under
[Inspecting the results](../reinforcement/fem.md#inspecting-the-results) on the reinforcement page, lists each member
with a utilization badge and draws a selected pile's profiles. Piles whose rows share a label are numbered so they
can be told apart. The screenshot is a strength reduction run on the two pile rows of
[LEM-12](../tutorials/lem12_piles.md), used again in [FEM-4](../tutorials/fem04_piles.md), at the mechanism it
developed:

![Pile detail for the lower pile of the piles sample](../fem/images/piles_fem_details.png){width=1000}

Four panels share one depth axis, with the pile head at the top. The first shows the lateral displacement, the
component normal to the pile axis. The second shows the shear $V$ after the `Vcap` limit, with its largest value
and its depth marked. The third shows the bending moment, assembled from the elements' end moments, with its maximum
and its depth marked. The
moment is zero at a free head and a free toe, which is a useful check on the profile.

The fourth panel shows the soil reaction, the sideways resistance the ground mobilizes per unit length of pile. It
comes from the three-node elements as described above, or on a linear mesh from the step in shear between elements.
Beside it, dashed, is the Ito & Matsui limit $p(z) = (c A_1 + \gamma z A_2)/S$, with $c$ and $\gamma$ those of the
soil at depth $z$ below the pile head, and $A_1$ and $A_2$ the coefficients the [LEM](lem.md#ito-matsui-1975-theory)
uses. The panel states the largest ratio of reaction to limit. The limit grows with depth and is often far above the
mobilized reaction. In that case, the panel is scaled to the reaction, and the limit runs off its sides. When the
limit misses the panel entirely, it is left out, and the note gives how far outside the panel the limit lies. The limit does not change with the
[Field state](../reinforcement/fem.md#inspecting-the-results) setting.

Capacity lines appear only where the model gives `Vcap` or `Mcap`. XSLOPE does not compute a capacity from an assumed
section. Where the model gives `D` and `S` but no capacities, the utilization badge compares the mobilized soil
reaction with the Ito & Matsui limit. With neither, the badge stays neutral.

The panels do not mark where the shear band crosses the pile. The soil loads a pile along its whole length, and its
moment peaks where the shear passes through zero, generally some distance from the crossing. The crossing is on the
shear-strain figure. The panel's title says which field the profiles come from: the mechanism a strength reduction
run captured, or the shear strain in a section that is standing.

## References

Cai, F., & Ugai, K. (2000). Numerical analysis of the stability of a slope reinforced with piles. *Soils and Foundations*, 40(1), 73-84. [doi:10.3208/sandf.40.73](https://doi.org/10.3208/sandf.40.73)

Ito, T., & Matsui, T. (1975). Methods to estimate lateral force acting on stabilizing piles. *Soils and Foundations*, 15(4), 43-59.
