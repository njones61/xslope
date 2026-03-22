# Pile Elements in Finite Element Analysis

## Introduction

In the finite element method, piles and concrete piers are modeled as **one-dimensional beam elements** embedded within the two-dimensional soil continuum. Unlike the limit equilibrium approach where the user provides a single force $H$ at the failure surface, the FEM approach models the pile as a structural member with distributed stiffness that naturally interacts with the surrounding soil through shared nodes. The pile's resistance to lateral soil movement emerges from the analysis itself — no pre-specified pile force is needed.

This approach captures several behaviors that are difficult or impossible to model in LEM:

- The pile's resistance is mobilized progressively as the soil deforms, rather than applied as a fixed force
- Both axial (along the pile) and lateral (perpendicular to the pile) stiffness are modeled
- The pile can carry both tension and compression (unlike flexible reinforcement, which is tension-only)
- The bending stiffness of the pile ($EI$) naturally resists lateral soil movement
- During SSRM strength reduction, the pile stiffness remains constant while soil strength is reduced, and the analysis reveals when the soil can no longer transfer load to the pile (failure)

For background on the general finite element slope stability methodology in XSLOPE, see the [FEM Overview](overview.md). For the LEM treatment of piles, including force resolution, typical parameters, and the Ito & Matsui method, see [LEM Piles](../lem/piles.md).


## Comparison with Reinforcement (Truss) Elements

The existing FEM module models flexible reinforcement as **2-node truss elements** (axial only, tension only). Piles reuse much of this infrastructure but differ in key ways:

| Property | Reinforcement (Truss) | Pile (Beam) |
|----------|----------------------|-------------|
| Axial stiffness | $EA/L$ | $EA/L$ |
| Lateral stiffness | None | $3EI/L^3$ (from bending) |
| Compression | Not allowed (zeroed via body-force corrections) | Allowed |
| Orientation | Inclined (along reinforcement) | Vertical or battered |
| Failure mode | Tension rupture / pullout | Shear / bending capacity |
| DOFs per node | 2 (translational only) | 2 (translational only, rotational condensed out) |

For details on the truss element formulation used for reinforcement, see [Soil Reinforcement](reinforcement.md).


## Beam Element Formulation

### Full Beam Stiffness

A standard Euler-Bernoulli beam element has 3 DOFs per node ($u_x$, $u_y$, $\theta$), giving a 6-DOF element with the following stiffness matrix in local coordinates (beam axis along local $x$):

>$\mathbf{K}_{\text{full}} = \begin{bmatrix} \frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\ 0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\ -\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\ 0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L} \end{bmatrix}$

where $E$ is Young's modulus, $A$ is the cross-sectional area, $I$ is the moment of inertia, and $L$ is the element length.

### Static Condensation

The 2D soil mesh in XSLOPE has only 2 DOFs per node ($u_x$, $u_y$). Adding rotational DOFs ($\theta$) at pile nodes would complicate the global system and require special handling. Instead, we use **static condensation** to eliminate the rotational DOFs from the beam stiffness matrix, expressing the element stiffness entirely in terms of translational DOFs.

Partitioning the full stiffness into translational ($tt$) and rotational ($rr$) blocks:

>$\mathbf{K}_{\text{condensed}} = \mathbf{K}_{tt} - \mathbf{K}_{tr} \, \mathbf{K}_{rr}^{-1} \, \mathbf{K}_{rt}$

For a beam element aligned along the local $x$-axis, this yields a $4 \times 4$ condensed stiffness matrix in terms of $[u_1, v_1, u_2, v_2]$:

>$\mathbf{K}_{\text{condensed}}^{\text{local}} = \begin{bmatrix} \frac{EA}{L} & 0 & -\frac{EA}{L} & 0 \\ 0 & \frac{3EI}{L^3} & 0 & -\frac{3EI}{L^3} \\ -\frac{EA}{L} & 0 & \frac{EA}{L} & 0 \\ 0 & -\frac{3EI}{L^3} & 0 & \frac{3EI}{L^3} \end{bmatrix}$

The lateral stiffness reduces from $12EI/L^3$ to $3EI/L^3$ after condensation. This corresponds to a beam with rotations free at both ends (pin-pin for rotation), which is a conservative single-element approximation. In practice, with multiple short elements along the pile connected at shared nodes, the chain of elements recovers the correct bending behavior because intermediate node translations are implicitly constrained by the series of elements.

### Coordinate Transformation

The condensed local stiffness matrix is transformed to global coordinates using the same direction cosine approach used for truss elements:

>$\mathbf{K}_{\text{global}} = \mathbf{T}^T \, \mathbf{K}_{\text{condensed}}^{\text{local}} \, \mathbf{T}$

where $\mathbf{T}$ is the $4 \times 4$ transformation matrix constructed from the element direction cosines $\cos\theta$ and $\sin\theta$.

For a **vertical pile** ($\cos\theta = 0$, $\sin\theta = 1$), the axial stiffness $EA/L$ acts in the $y$-direction (vertical) and the lateral stiffness $3EI/L^3$ acts in the $x$-direction (horizontal) — exactly the behavior needed for resisting lateral slope movement.

### Assembly

The global stiffness matrix combines contributions from all element types:

>$[\mathbf{K}]_{\text{global}} = \sum_{\text{soil}} [\mathbf{K}_e]_{\text{soil}} + \sum_{\text{truss}} [\mathbf{K}_e]_{\text{truss}} + \sum_{\text{beam}} [\mathbf{K}_e]_{\text{beam}}$

This matrix is factored once (via sparse LU decomposition) and reused for all viscoplastic iterations.


## Mesh Generation

Pile elements are integrated into the finite element mesh using the same approach as reinforcement elements. The pile line geometry — defined as two endpoints $(x_1, y_1)$ to $(x_2, y_2)$ in the `piles` sheet of the input template — is passed to the mesh generator as a constraint. The mesh generator ensures that element edges align along the pile line, so that 1D beam elements can be extracted from the edges of adjacent 2D soil elements.

The existing `extract_1d_elements_from_2d_edges()` function identifies mesh edges that are collinear with a given line and extracts them as 1D elements with shared nodes. Since pile lines use the same two-endpoint format as reinforcement lines, this function handles both vertical and battered piles without modification.

For details on the mesh generation process and 1D element extraction, see [Automated Mesh Generation](mesh.md).


## Integration with Viscoplastic Iteration

The pile beam elements participate in the same viscoplastic iteration loop used for soil elements and reinforcement truss elements (see [FEM Overview](overview.md) for details on the viscoplastic algorithm). The key differences from reinforcement elements are:

### Force Computation

At each viscoplastic iteration, after solving for the updated displacement field, the forces in each beam element are computed from the current nodal displacements:

>**Axial force**: $T = \dfrac{EA}{L} \cdot \delta_{\text{axial}}$

>**Lateral (shear) force**: $V = \dfrac{3EI}{L^3} \cdot \delta_{\text{lateral}}$

where $\delta_{\text{axial}}$ and $\delta_{\text{lateral}}$ are the projections of relative nodal displacement along and perpendicular to the element axis, respectively.

### Compression Allowed

Unlike reinforcement truss elements, which zero out compressive forces through body-force corrections, beam elements carry **both tension and compression**. No compression correction is applied.

### Structural Capacity

Beam elements are currently **linearly elastic** — they have no structural strength limit. The SSRM finds the factor of safety at which soil fails around the pile, which is the correct result when the pile is strong enough that soil failure governs.

Users should check the reported beam element forces (max axial force, max lateral force) against the pile's structural shear and moment capacity as a post-processing step. If the beam forces at the SSRM failure state exceed the pile's structural capacity, the pile would fail before the soil and the reported FS is unconservative.

Elastic-perfectly-plastic capacity checks (V_cap, M_cap) are planned for a future release. These would cap beam element forces and redistribute the excess through body-force corrections, following the same viscoplastic pattern used for soil elements.


## SSRM Treatment

During the Shear Strength Reduction Method:

- Pile material properties ($E$, $I$, $A$, yield strength) are **not reduced** — they remain constant throughout the SSRM bisection
- Only soil strength parameters ($c$ and $\tan\varphi$) are reduced by the trial factor $F$
- As soil strength is reduced, the soil deforms more around the pile; the pile's stiffness naturally resists this additional deformation
- Convergence failure (slope collapse) occurs when the reduced soil strength can no longer transfer load to/from the piles effectively
- The resulting factor of safety represents the margin of safety in the **soil** strength, given the pile reinforcement as-designed

This treatment is consistent with how reinforcement elements are handled in the SSRM (see [Soil Reinforcement](reinforcement.md)) and follows standard practice in commercial geotechnical FEM software.


## Input Parameters

Pile properties for FEM analysis are specified in the same `piles` sheet used for LEM analysis. The FEM-specific columns are:

| Column | Field | Description |
|--------|-------|-------------|
| B-E | $(x_1, y_1)$, $(x_2, y_2)$ | Pile geometry (top and tip coordinates) |
| H | $D_{\text{pile}}$ | Pile diameter — used to compute $I$ and $A$ if not provided |
| I | $S$ | Center-to-center spacing — used to scale $EI$ and $EA$ by $1/S$ for per-unit-width |
| J | $E$ | Young's modulus of pile material |
| K | $I$ | Moment of inertia of pile cross-section |
| L | $Area$ | Cross-sectional area of pile cross-section |

If $D_{\text{pile}}$ is provided and $I$/$Area$ are omitted, a solid circular section is assumed:

>$A = \dfrac{\pi D^2}{4}, \qquad I = \dfrac{\pi D^4}{64}$

The LEM-specific columns ($H$, $\theta_p$) are not used for FEM analysis. See [LEM Piles](../lem/piles.md) for typical material property values.

## References

Griffiths, D.V., & Lane, P.A. (1999). Slope stability analysis by finite elements. *Geotechnique*, 49(3), 387-403.

Smith, I.M., & Griffiths, D.V. (2004). *Programming the Finite Element Method* (4th ed.). John Wiley & Sons.

Sun, G., Lin, S., Jiang, W., & Yang, Y. (2021). A simplified solution for determining the factor of safety of a slope reinforced with piles based on the shear strength reduction method. *Bulletin of Engineering Geology and the Environment*, 80, 7719-7730.
