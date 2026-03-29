# Finite Element Method for Slope Stability Analysis

## Introduction

The finite element method (FEM) provides a powerful numerical technique for slope stability analysis that overcomes 
many fundamental limitations of traditional limit equilibrium methods. While limit equilibrium approaches require the engineer to assume a failure surface geometry and then check whether equilibrium conditions are satisfied, FEM allows potential failure mechanisms to emerge naturally through rigorous stress analysis and progressive failure development (Griffiths & Lane, 1999; Duncan, 1996). Rather than imposing kinematic constraints through assumed failure surfaces, FEM solves the complete stress-strain problem throughout the slope domain. The method can capture the complex stress redistribution that occurs as soil elements progressively reach failure, leading to the natural development of failure zones without prior assumptions about their geometry or location. Perhaps most importantly, FEM uses realistic stress-strain constitutive models that can capture the actual behavior of soil materials, including nonlinear elastic behavior, plastic yielding, strain softening, and progressive failure. This provides a much more accurate representation of soil response compared to the rigid-perfectly plastic assumptions typically used in limit equilibrium methods.

## Governing Equations

### Equilibrium Equations

The foundation of finite element slope stability analysis rests on the fundamental equilibrium equations that govern the static behavior of continuum mechanics. In two dimensions, these equilibrium equations express the requirement that forces acting on any infinitesimal element of soil must be in balance:

>>$\dfrac{\partial \sigma_x}{\partial x} + \dfrac{\partial \tau_{xy}}{\partial y} + b_x = 0$

>>$\dfrac{\partial \tau_{xy}}{\partial x} + \dfrac{\partial \sigma_y}{\partial y} + b_y = 0$

Here, $\sigma_x$ and $\sigma_y$ represent the normal stresses acting in the x and y directions respectively, while $\tau_{xy}$ denotes the shear stress. The body force terms $b_x$ and $b_y$ account for forces distributed throughout the volume of the material, with gravity being the most common example where $b_x = 0$ and $b_y = -\gamma$, where $\gamma$ is the unit weight of the soil.

These equations must be satisfied at every point within the slope domain for the system to be in static equilibrium. The challenge in slope stability analysis arises because soil materials exhibit nonlinear, inelastic behavior that violates these equilibrium conditions when failure occurs, leading to the progressive development of failure zones.

### Stress-Strain Relations

The constitutive behavior of soil in finite element slope stability analysis is typically modeled using an elastic-perfectly plastic framework that combines linear elastic behavior with Mohr-Coulomb plasticity. This approach recognizes that soil behaves elastically under small stress changes but exhibits permanent deformation once failure is reached.

During the elastic phase, before any yielding occurs, the relationship between stress and strain follows Hooke's law expressed in matrix form:

>>$\{\sigma\} = [D_e] \{\varepsilon\}$

The elastic constitutive matrix $[D_e]$ for plane strain conditions, which is most appropriate for slope stability problems, takes the form:

>>$[D_e] = \dfrac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & \dfrac{1-2\nu}{2}
\end{bmatrix}$

This formulation requires several fundamental material properties that must be determined through laboratory testing or empirical correlations. Young's modulus $E$ governs the stiffness of the soil under loading, while Poisson's ratio $\nu$ controls the relationship between axial and lateral strains. For slope stability analysis, additional strength parameters are critical: the cohesion $c$ and friction angle $\phi$ define the failure envelope, while the unit weight $\gamma$ determines the gravitational body forces. The coefficient of earth pressure at rest $K_0$ is often needed to establish initial stress conditions, particularly for natural slopes that have developed under gravitational loading over geological time.

#### Typical Elastic Parameters for Finite Element Analysis

The following table provides typical ranges of elastic parameters for common soil types under **drained** conditions. 
These values should be used as initial estimates and refined based on site-specific testing when available.

| Soil Type | Young's Modulus $E$ [MPa] | Young's Modulus $E$ [kPa] | Young's Modulus $E$ [psf] | Poisson's Ratio $\nu$ | Notes |
|-----------|:-------------------------:|:-------------------------:|:-------------------------:|:--------------------:|-----------------|
| **Soft Clay** | 2 - 15 | 2,000 - 15,000 | 41,800 - 313,000 | 0.40 - 0.50 | Use lower E values for very soft clays |
| **Medium Clay** | 15 - 50 | 15,000 - 50,000 | 313,000 - 1,044,000 | 0.35 - 0.45 | Plasticity index affects stiffness |
| **Stiff Clay** | 50 - 200 | 50,000 - 200,000 | 1,044,000 - 4,175,000 | 0.20 - 0.40 | Overconsolidated clays have higher E |
| **Loose Sand** | 10 - 25 | 10,000 - 25,000 | 209,000 - 522,000 | 0.25 - 0.35 | Depends on relative density |
| **Medium Sand** | 25 - 75 | 25,000 - 75,000 | 522,000 - 1,565,000 | 0.30 - 0.40 | Well-graded sands toward upper range |
| **Dense Sand** | 75 - 200 | 75,000 - 200,000 | 1,565,000 - 4,175,000 | 0.25 - 0.35 | Angular particles give higher stiffness |
| **Loose Silt** | 5 - 20 | 5,000 - 20,000 | 104,000 - 418,000 | 0.30 - 0.45 | Non-plastic silts toward lower ν |
| **Dense Silt** | 20 - 100 | 20,000 - 100,000 | 418,000 - 2,088,000 | 0.25 - 0.40 | Cementation increases stiffness |
| **Gravel** | 100 - 500 | 100,000 - 500,000 | 2,088,000 - 10,440,000 | 0.15 - 0.30 | Well-graded, dense materials |
| **Rock Fill** | 50 - 300 | 50,000 - 300,000 | 1,044,000 - 6,260,000 | 0.20 - 0.35 | Depends on gradation and compaction |
| **Soft Rock** | 1,000 - 10,000 | 1,000,000 - 10,000,000 | 20,880,000 - 208,800,000 | 0.15 - 0.30 | Weathered or fractured rock |

The first column for E is for reference only. When working in metric units, E should always be entered in $Kpa$ to be consistent with the unit weights and cohesion values. For English units, E should always be entered in $psf$. 

   For undrained shearing conditions, the undrained modulus can be determined through several approaches:
  
  - **Direct Laboratory Testing**: Unconsolidated undrained (UU) triaxial tests or unconfined compression tests provide direct measurement of $E_u$. The initial tangent modulus from stress-strain curves gives the most representative value.
  
  - **Empirical Correlations**: For clays, $E_u$ can be estimated from undrained shear strength using:

    >>$E_u = (150-1500) \times S_u$

    where $S_u$ is the undrained shear strength. Use lower multipliers (150-400) for soft clays and higher values (400-1500) for stiff clays.
  
  - **Relationship to Drained Modulus**: For saturated clays, the undrained modulus is typically higher than the drained modulus:

    >>$E_u = \dfrac{E'(1-2\nu')}{(1-2\nu_u)(1+\nu')} \times \dfrac{(1+\nu_u)}{(1-\nu')}$

    where $E'$ and $\nu'$ are drained parameters, and $\nu_u \approx 0.5$ for undrained conditions. This often simplifies to $E_u \approx (1.2-2.0) \times E'$ for typical clay parameters.

- **Laboratory vs. Field Values**: Laboratory-derived moduli are often higher than field values due to sample disturbance and scale effects. Field moduli (from pressuremeter, plate load tests) may be more representative.

- **Empirical Correlations**: When direct testing is unavailable, Young's modulus can be estimated from standard penetration test (SPT) or cone penetration test (CPT) data using published correlations.

!!! NOTE
    In finite element slope stability analysis using the Shear Strength Reduction Method (SSRM), the elastic modulus E primarily affects the magnitude of computed deformations but has minimal impact on the calculated factor of safety. The factor of safety is governed by the strength parameters (cohesion c and friction angle φ) rather than the elastic response. Therefore, approximate values of E are often sufficient for slope stability calculations, and extensive effort to precisely determine elastic moduli may not be warranted unless deformation predictions are also required.

### Mohr-Coulomb Failure Criterion

The constitutive behavior described above is valid only while the soil remains in the elastic domain. Once the stress state reaches the failure envelope, plastic yielding occurs and the material behavior changes fundamentally. The Mohr-Coulomb failure criterion forms the theoretical foundation for determining when this failure occurs in finite element slope stability analysis.

This criterion, developed from extensive experimental observations of soil behavior, recognizes that soil failure is fundamentally a shear phenomenon that depends on both the normal stress acting on the failure plane and the inherent strength properties of the material. The following figure illustrates the Mohr-Coulomb failure envelope in stress space, which defines the boundary between stable and unstable states for a given soil material:

![mc_envelope.png](images/mc_envelope.png){width=800px}

The line defined by the Mohr-Coulomb criterion represents the maximum shear stress that can be sustained by the soil at a given effective normal stress. The slope of this line is determined by the angle of internal friction $\phi$, while the intercept on the shear stress axis is defined by the cohesion $c$ of the soil. The basic form of the Mohr-Coulomb criterion expresses the relationship between shear strength and normal stress on any potential failure plane:

>>$\tau_f = c + \sigma \tan \phi$

In this formulation, $\tau_f$ represents the shear strength available to resist failure, $c$ is the cohesion representing the portion of strength that is independent of normal stress, $\sigma$ is the normal stress acting perpendicular to the failure plane, and $\phi$ is the angle of internal friction that governs how strength increases with normal stress. When written in terms of effective stresses, the criterion becomes:

>>$\tau_f = c + \sigma' \tan \phi$

>>$\tau_f = c + (\sigma - u) \tan \phi$

where $\sigma'$ is the effective normal stress, defined as the total normal stress minus pore water pressure (u). This effective stress formulation is critical for understanding how changes in pore water pressure, such as those caused by rainfall or groundwater fluctuations, can significantly affect slope stability.

**Yield Function for Finite Element Implementation:**

For computational implementation in finite element analysis, it is more convenient to express the failure criterion in terms of principal stresses. This transformation yields the yield function:

>>$f(\sigma_1', \sigma_3') = \dfrac{\sigma_1' - \sigma_3'}{2} - \left(\dfrac{\sigma_1' + \sigma_3'}{2} \sin \phi + c \cos \phi\right) = 0$

where $\sigma_1'$ and $\sigma_3'$ are the major and minor principal effective stresses respectively. This yield function defines the boundary between elastic and plastic behavior in the soil:

- When $f < 0$: stress state lies within the elastic domain
- When $f = 0$: stress state lies exactly on the yield surface (incipient yielding)
- When $f > 0$: stress state has exceeded the yield strength (plastic deformation required)

The yield function can be visualized in principal stress space, where the failure envelope is a hyperbolic surface that separates stable and unstable states:

![yield_surface.png](images/yield_surface.png)

This principal stress formulation allows direct evaluation of the yield function using the principal stresses computed at each integration point within the finite element mesh. The evaluation of principal stresses $\sigma_1'$ and $\sigma_3'$ from the general stress tensor requires solution of the eigenvalue problem, which can be computationally intensive but is essential for accurate yield detection.

The implementation of this failure criterion within the finite element framework will be discussed after the basic finite element formulation is presented.

## Finite Element Formulation

### Discretization

The transformation from continuous domain to discrete finite element system begins with dividing the slope domain 
into a collection of simple geometric elements, typically triangles or quadrilaterals.

![sample_mesh.png](images/sample_mesh.png)

This discretization process 
is fundamental to the finite element method because it allows the complex, continuous displacement field throughout 
the slope to be approximated using simple polynomial functions defined over each element. The element types 
supported in XSLOPE include linear triangular elements, quadratic triangular elements, linear quadrilateral elements, and quadratic quadrilateral elements. Each element type has its own set of shape functions that define how displacements vary within the element based on the nodal values.

![element_types.png](images/element_types.png){width=600px}

Within each element, the displacement field is interpolated from the nodal displacement values using shape functions. For a typical two-dimensional element, the horizontal and vertical displacements at any point within the element are expressed as:

>>$u = [N] \{u_e\}$

>>$v = [N] \{v_e\}$

The shape function matrix $[N]$ contains the interpolation functions that define how displacements vary spatially within the element, while $\{u_e\}$ and $\{v_e\}$ are vectors containing the nodal displacement values. For linear triangular elements, the shape functions are simply the area coordinates that ensure displacement compatibility between adjacent elements and provide a linear variation of displacement within each element.

The choice of element type significantly impacts both accuracy and computational efficiency. Triangular elements 
with linear shape functions are particularly well-suited for slope stability problems because they can easily conform to irregular slope geometries and provide adequate accuracy for capturing the stress distributions that govern failure development. The linear displacement variation within each element leads to constant strain and stress fields, which is appropriate for modeling the elastic-perfectly plastic soil behavior typically assumed in slope stability analysis. However, linear triangles are susceptible to **volumetric locking** — an artificial stiffness that can significantly overestimate the factor of safety, particularly for narrow elements or materials approaching incompressibility (see [Element Type Selection and Volumetric Locking](#element-type-selection-and-volumetric-locking) below). The tools for building finite element meshes in XSLOPE are described in the [Automated Mesh Generation](mesh.md) page. 

### Element Stiffness Matrix

The element stiffness matrix represents the fundamental relationship between nodal forces and nodal displacements for each finite element. This matrix is derived through application of the principle of virtual work, which states that for a system in equilibrium, the virtual work done by external forces must equal the virtual work done by internal stresses for any kinematically admissible virtual displacement field. The mathematical expression for the element stiffness matrix emerges from this principle as:

>>$[K_e] = \int_{A_e} [B]^T [D_e] [B] \, dA$

This elegant formulation embodies the essential physics of the problem. The strain-displacement matrix $[B]$ 
transforms nodal displacements into strains throughout the element, while the constitutive matrix $[D_e]$ relates 
these strains to stresses according to the material's stress-strain behavior. The integration over the element area $A_e$ ensures that the stiffness contribution from every point within the element is properly accounted for.

For the commonly used linear triangular elements, the strain-displacement matrix takes the specific form:

>>$[B] = \dfrac{1}{2A} \begin{bmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{bmatrix}$

The coefficients $b_i$ and $c_i$ are geometric constants determined by the nodal coordinates of the triangle, and $A$ represents the triangle area. This matrix remains constant throughout the element because of the linear nature of the shape functions, which simplifies the integration process and leads to computational efficiency.

### Global System Assembly

The transition from individual element stiffness matrices to the global system of equations represents one of the most elegant aspects of the finite element method. Each element contributes to the overall structural behavior according to its connectivity with other elements, creating a sparse but symmetric global stiffness matrix that captures the mechanical interaction throughout the entire slope domain.

The global system of equations takes the familiar form:

>>$[K] \{U\} = \{F\}$

The global stiffness matrix $[K]$ is assembled by systematically adding each element's stiffness contribution to the appropriate locations corresponding to the degrees of freedom associated with that element's nodes. This assembly process ensures displacement compatibility between adjacent elements and force equilibrium at every node in the mesh.

The global displacement vector $\{U\}$ contains the unknown nodal displacements for the entire mesh, while the global force vector $\{F\}$ represents the applied loads including both external forces and body forces due to gravity. The sparsity of the global stiffness matrix, where most entries are zero due to the local connectivity of finite elements, allows efficient solution algorithms to be employed even for large-scale slope stability problems.

The solution of this global system provides the displacement field throughout the slope under the applied loading conditions. From these displacements, the strain and stress fields can be computed at every point in the domain, enabling assessment of the proximity to failure according to the chosen yield criterion.

## Boundary Conditions

The proper specification of boundary conditions is crucial for obtaining physically meaningful solutions in finite element slope stability analysis. Boundary conditions define how the slope interacts with its surroundings and constrain the displacement field to reflect realistic physical constraints. The choice of boundary conditions significantly affects both the stress distribution within the slope and the computed factor of safety.

### Displacement Boundary Conditions

Displacement boundary conditions are applied where the motion of the soil mass is constrained by physical limitations or where the model boundaries must represent the behavior of the extended soil mass beyond the computational domain.

**Fixed supports** represent locations where both horizontal and vertical displacements are completely prevented, 
typically expressed as $u = 0$ and $v = 0$. These conditions are most commonly applied at the base of the finite element model when the analysis extends to sufficient depth that the displacement of deep soil layers has negligible effect on slope stability. The depth required for this assumption depends on the slope geometry and soil properties, but generally the model should extend at least one slope height below the toe and preferably to bedrock or very stiff soil layers.

**Roller supports** constrain displacement in only one direction while allowing free movement in the perpendicular direction. Along vertical side boundaries, horizontal displacement is typically prevented ($u = 0$) while vertical movement is allowed, reflecting the assumption that the slope extends laterally beyond the model boundaries with similar geometry and loading conditions. At the model base, vertical displacement may be prevented ($v = 0$) while allowing horizontal movement, which is appropriate when the analysis does not extend to truly fixed boundary conditions.

**Free boundaries** occur along the ground surface and slope face where no external constraints are applied. These boundaries represent the natural boundary condition of zero traction, meaning that no external forces act normal or tangential to these surfaces except for applied loads such as surcharge loads or foundation pressures.

For slope stability analysis, the most common displacement boundary conditions are fixed supports at the 
base of the 
model, roller supports (free movement in the vertical direction) along the left and right vertical boundaries, and free 
boundaries along the slope face and ground surface. These 
conditions ensure 
that the 
model 
accurately reflects the physical constraints of the slope while allowing for realistic stress distributions and 
potential failure mechanisms to develop. These displacement boundary conditions are automatically applied in XSLOPE.

### Force Boundary Conditions

Force boundary conditions specify the external loads acting on the slope. In XSLOPE, force boundary conditions 
correspond to distributed loads which are defined in the input template as as line loads with a sequence of coordinates and corresponding load values (force per unit length). These can represent:

- **Traffic loads** (vehicles, equipment)
- **Structural loads** (buildings, foundations)  
- **Hydrostatic pressure** (water on slope face)

The distributed loads in the XSLOPE input template can be used either for limit equilibrium analysis or finite element 
analysis. For limit equilibrium analysis, each distributed load is converted to a resultant force applied at the top of each 
slice. The total load on each slice is calculated by integrating the distributed load over the slice width. 

For finite element analysis, the same distributed loads are converted to equivalent nodal forces. Since distributed loads in XSLOPE always vary linearly along the ground surface (creating simple trapezoidal load distributions), the conversion to nodal forces can be simplified using trapezoidal integration rather than general finite element shape functions. For a linear load distribution between two adjacent nodes with load intensities $q_1$ and $q_2$ separated by distance $L$, the equivalent nodal forces perpendicular to the ground surface are:

>>$F_1 = \frac{L}{6}(2q_1 + q_2)$

>>$F_2 = \frac{L}{6}(q_1 + 2q_2)$

Since distributed loads always act perpendicular to the ground surface, these forces must be decomposed into horizontal and vertical components at each node. For a ground surface segment with slope angle $\beta$ (measured from horizontal), the nodal force components are:

>>$F_{1x} = -F_1 \sin \beta$

>>$F_{1y} = -F_1 \cos \beta$

>>$F_{2x} = -F_2 \sin \beta$

>>$F_{2y} = -F_2 \cos \beta$

where the negative signs indicate that the loads act downward and into the slope (typical for traffic or structural loads). The slope angle $\beta$ is calculated from the coordinates of adjacent nodes along the ground surface. For a given node on the ground surface, the total resultant force is the sum of the contributions of the distributed load on both the left and right sides of the node.

#### Body Forces

Body forces act throughout the volume of soil elements, primarily gravitational forces (self-weight of soil). For gravitational loading, the body force components are:

>>$b_x = 0$<br>
$b_y = -\gamma$

where $\gamma$ is the unit weight of the soil. These body forces are incorporated into the equilibrium equations:

>>$\dfrac{\partial \sigma_x}{\partial x} + \dfrac{\partial \tau_{xy}}{\partial y} + b_x = 0$

>>$\dfrac{\partial \tau_{xy}}{\partial x} + \dfrac{\partial \sigma_y}{\partial y} + b_y = 0$

In the finite element formulation, body forces are converted to equivalent nodal forces using:

>>$\{F\}_b = \sum_{e=1}^{N_{elem}} \int_{A_e} [N]^T \{b\} \, dA$

where $[N]$ are the shape functions and the integration is performed over each element area $A_e$, then summed over all elements in the mesh. This integration distributes the self-weight of the soil to the nodes of each element, ensuring that the gravitational loading is properly represented throughout the slope domain.

### Implementation of Boundary Conditions

The boundary conditions described above must be incorporated into the global system of equations $[K]{U} = {F}$ to 
obtain a solvable system. The implementation of boundary conditions fundamentally modifies both the stiffness matrix and force vector.

**Displacement Boundary Conditions:**

For prescribed displacements (such as u = 0 or v = 0), the most common implementation approach is the penalty method or direct modification:

1. **Direct Modification Method:**<br>
- For a node with prescribed displacement $U_i = 0$, replace row i of the stiffness matrix with zeros except for the diagonal term, which is set to a large number<br>
- Set the corresponding force term $F_i = 0$<br>
- This forces the solution to satisfy the constraint $U_i = 0$<br><br>

2. **Example:** For a node at the base with both u = 0 and v = 0:
   >$\begin{bmatrix} K_{11} & K_{12} & \cdots \\ 0 & 1 \times 10^{12} & 0 & \cdots \\ 0 & 0 & 1 \times 10^{12} & \cdots \\ \vdots & \vdots & \vdots & \ddots \end{bmatrix} \begin{bmatrix} U_1 \\ 0 \\ 0 \\ \vdots \end{bmatrix} = \begin{bmatrix} F_1 \\ 0 \\ 0 \\ \vdots \end{bmatrix}$

**Force Boundary Conditions:**

Applied forces and distributed loads are incorporated directly into the global force vector {F} through the integration processes described above. The stiffness matrix [K] remains unchanged for force boundary conditions.

### Automatic Boundary Condition Assignment in XSLOPE

XSLOPE automatically assigns displacement boundary conditions in the `build_fem_data()` function based on the geometry of the mesh. The user does not need to specify boundary conditions manually — they are determined entirely from the node coordinates using the following rules applied in order:

1. **All nodes start as free** (no displacement constraints). The ground surface, slope face, and any internal nodes are unconstrained by default, representing the natural zero-traction boundary condition.

2. **Fixed supports at the base.** All nodes at the global minimum $y$-coordinate are assigned fixed boundary conditions ($u = 0$, $v = 0$). These nodes are identified by comparing each node's $y$-coordinate to the minimum value in the mesh within a small numerical tolerance ($10^{-6}$). This represents rigid bedrock or a sufficiently deep boundary where displacements are negligible.

3. **X-roller supports on the left and right sides.** All nodes at the global minimum $x$-coordinate (left boundary) and maximum $x$-coordinate (right boundary) are assigned x-roller conditions ($u = 0$, $v$ free). This allows vertical settlement along the side boundaries while preventing lateral movement, reflecting the assumption that the slope extends indefinitely in both directions beyond the model domain. Corner nodes where the side boundaries meet the base retain their fixed condition — fixed supports take precedence over rollers.

4. **Force boundary conditions from distributed loads.** If distributed loads are defined in the input template, nodes along the load lines are identified and assigned equivalent nodal forces computed using tributary-area integration (described above). If a force node coincides with a roller or fixed boundary, both the displacement constraint and the applied force are preserved.

The figure below shows the resulting boundary conditions for the reinforced slope example from the [FEM Samples](samples.md) page (Problem 2). Fixed supports (triangles) line the base of the mesh. X-roller supports (circles) line the left and right vertical boundaries, allowing vertical movement but preventing horizontal displacement. The ground surface and slope face are free. Force boundary conditions from a 240 psf surcharge are shown as arrows along the slope crest. Reinforcement elements are shown as red lines within the slope body.

![reinforce_fem_mesh.png](images/reinforce_fem_mesh.png){width=1000}

## Elastic-Plastic Behavior: Viscoplastic Algorithm

The finite element formulation described above assumes purely elastic behavior governed by the elastic constitutive
matrix $[D_e]$. However, when soil elements reach the failure envelope defined by the Mohr-Coulomb criterion,
plastic deformation occurs. XSLOPE handles this using the **viscoplastic algorithm** described by Griffiths & Lane (1999) and Smith & Griffiths (2004), which provides a robust and elegant approach to elastic-perfectly plastic finite element analysis.

### Overview of the Viscoplastic Approach

Unlike traditional elastic-plastic methods that modify the stiffness matrix when elements yield, the viscoplastic algorithm keeps the elastic stiffness matrix $[K]$ constant throughout the analysis. Plastic behavior is instead handled through accumulated viscoplastic strains that generate body load corrections added to the right-hand side of the equilibrium equations. This has two major advantages:

1. The stiffness matrix is assembled and factorized only once, then reused for all iterations via back-substitution
2. The algorithm is unconditionally stable for appropriate choice of the time step parameter

### Viscoplastic Iteration Process

For a given set of strength parameters (c, φ) and applied loads, the solution proceeds as follows:

**1. Setup (performed once):**<br>
>- Assemble global elastic stiffness matrix $[K]$ using $[D_e]$ for all elements<br>
- Build gravity load vector $\{F\}_{gravity}$<br>
- Pre-factorize $[K]$ using sparse LU decomposition (reused for all iterations)<br>
- Initialize viscoplastic strains $\{\varepsilon^{vp}\} = 0$ at all Gauss points

**2. Initial Elastic Solution:**<br>
>- Solve $[K]\{U\} = \{F\}_{gravity}$ using the pre-factored matrix

**3. Viscoplastic Iteration Loop:**<br>

Each iteration builds a corrected load vector and solves the full system:

>a. Start with gravity loads: $\{F\} = \{F\}_{gravity}$<br>
>b. At each Gauss point in every element:<br>
>>- Compute total strains from current displacements: $\{\varepsilon\} = [B]\{u_e\}$<br>
>>- Compute elastic strains: $\{\varepsilon^{el}\} = \{\varepsilon\} - \{\varepsilon^{vp}\}$<br>
>>- Compute stress from elastic strains: $\{\sigma\} = [D_e]\{\varepsilon^{el}\}$<br>
>>- Evaluate Mohr-Coulomb yield function using effective stress: $f = \tau_{max} - (\bar{\sigma} - u)\sin\phi - c\cos\phi$<br>
>>- If $f > 0$ (yielding): compute viscoplastic strain increment $\Delta\varepsilon^{vp} = f \cdot \dfrac{\partial Q}{\partial \sigma} \cdot \Delta t$<br>
>>- Accumulate: $\{\varepsilon^{vp}\} \mathrel{+}= \Delta\varepsilon^{vp}$<br>
>
>c. Add body load corrections from viscoplastic strains:
>>$\{F\} \mathrel{+}= \sum_{elements} \int [B]^T [D_e] \{\varepsilon^{vp}\} \, dA$<br>
>
>d. Solve $[K]\{U_{new}\} = \{F\}$ using the pre-factored matrix<br>
>e. Check convergence (see below)<br>
>f. If not converged, repeat from step (a) with updated displacements

**Pore Pressure at Gauss Points:**

>The pore pressure $u$ in the yield function is evaluated at each Gauss point using the pore pressure option specified for each material in the input template. Three options are available:
>
>- **None** ($u = 0$): No pore pressure. The yield check uses total stress.
>- **Piezometric line** (`"piezo"`): The physical coordinates of the Gauss point are computed via the shape functions ($x_{gp} = \sum N_i x_i$, $y_{gp} = \sum N_i y_i$), and the pore pressure is calculated as $u = \gamma_w (z_{piezo} - y_{gp})$ where $z_{piezo}$ is the elevation of the piezometric surface at the Gauss point's horizontal position. Negative values (above the piezometric surface) are clamped to zero.
>- **Seepage solution** (`"seep"`): Pore pressures from a prior seepage analysis (stored as nodal values) are interpolated to each Gauss point using the element shape functions: $u_{gp} = \sum N_i \cdot u_i$. Negative values are clamped to zero, since the FEM stress analysis uses effective stress and suction (negative pore pressure) is conservatively ignored.

**Key Parameters:**<br>

>- **Time step** $\Delta t$: A numerical parameter (not physical time) that controls stability. Following Smith & Griffiths: $\Delta t = \dfrac{4(1+\nu)}{3E}$<br>
>- **Flow rule** $\dfrac{\partial Q}{\partial \sigma}$: XSLOPE uses non-associated flow with dilation angle $\psi = 0$, producing purely deviatoric plastic strains (no volume change)<br>
>- **Maximum iterations**: 500 (sufficient for convergence near the critical factor of safety, where hundreds of iterations may be needed)

**Key Points:**<br>
>- **Constant stiffness**: $[K]$ never changes — all nonlinearity enters through the body load corrections<br>
>- **Mixed response**: The slope contains both elastic Gauss points ($f < 0$) and yielding Gauss points ($f > 0$) simultaneously<br>
>- **Load redistribution**: As Gauss points yield and accumulate viscoplastic strains, stress redistributes to surrounding regions through the body load correction mechanism<br>
>- **Elastic strains matter**: The yield check uses stress from elastic strains $\{\varepsilon\} - \{\varepsilon^{vp}\}$, not total strains — this correctly accounts for stress relief from already-accumulated viscoplastic deformation

### Convergence Criterion

The viscoplastic iteration loop requires a convergence criterion to determine when the stress redistribution has reached equilibrium. XSLOPE uses an elastic-relative displacement criterion that normalizes the iteration-to-iteration displacement change by the elastic displacement magnitude:

>>$\dfrac{||\{U\}_{i+1} - \{U\}_i||}{||\{U\}_{elastic}||} < \text{tol}$

where $\{U\}_{elastic}$ is the displacement vector from the initial elastic solution (before any viscoplastic iterations). This differs from the conventional relative criterion of Dawson et al. (1999), which normalizes by the current displacement $||\{U\}_{i+1}||$. The elastic-relative formulation has a critical advantage: the denominator is a fixed reference that does not grow as plastic displacements accumulate. With the conventional criterion, when a slope is failing and displacements become very large, the relative change $||\Delta U|| / ||U||$ can appear small even though the absolute change is enormous — producing **false convergence** where the solver reports a converged solution despite unbounded plastic flow. By normalizing against the elastic displacement, the convergence check remains meaningful regardless of the displacement magnitude.

**Implementation in XSLOPE:**

>- Default tolerance: $\text{tol} = 10^{-3}$<br>
>- Maximum iterations: 500<br>
>- Failure to converge within the maximum number of iterations indicates that the current strength reduction factor corresponds to an unstable configuration

Near the critical factor of safety, the viscoplastic algorithm may require hundreds of iterations to converge. For example, Griffiths & Lane (1999) report 792 iterations at a reduction factor just below failure for their Example 1. The maximum iteration count of 500 in XSLOPE provides sufficient margin for most practical problems.

### The `solve_fem()` Function

The viscoplastic algorithm described above is implemented in the `solve_fem()` function in the `fem.py` module. This function takes a FEM data dictionary (built by `build_fem_data()`) and an optional strength reduction factor $F$, assembles and factors the elastic stiffness matrix once, then runs the viscoplastic iteration loop until convergence or the maximum iteration count is reached.

```python
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons
from xslope.fem import build_fem_data, solve_fem

# Load slope data from input template
slope_data = load_slope_data("inputs/slope/input_template.xlsx")

# Build mesh (quad8 elements recommended — see Element Type Selection below)
mesh = build_mesh_from_polygons(polygons, target_size=4, element_type='quad8')

# Build FEM data dictionary
fem_data = build_fem_data(slope_data, mesh)

# Solve with unreduced strength (F=1.0)
solution = solve_fem(fem_data, F=1.0, debug_level=1)

if solution['converged']:
    print(f"Converged in {solution['iterations']} iterations")
    print(f"Max displacement: {solution['max_displacement']:.6f}")
else:
    print(f"Did not converge after {solution['iterations']} iterations")
```

The key parameters of `solve_fem()` are:

>- **`F`** (default 1.0): Strength reduction factor applied as $c_r = c/F$ and $\tan\phi_r = \tan\phi / F$. When called directly, $F = 1.0$ gives the unreduced solution; values of $F > 1.0$ can be used to test stability at specific reduction levels.<br>
>- **`max_iterations`** (default 500): Maximum viscoplastic iterations before declaring non-convergence.<br>
>- **`tolerance`** (default $10^{-3}$): Convergence tolerance for the elastic-relative displacement criterion described above.<br>
>- **`max_disp_factor`** (default 0.1): If the maximum viscoplastic displacement exceeds this fraction of the mesh height, the solve is terminated early and declared non-converged. Set to `None` to disable this check.<br>
>- **`debug_level`** (default 0): Controls output verbosity — 0 for silent, 1 for a summary, 2 for per-iteration details.

The returned solution dictionary contains the convergence status, nodal displacements, element stresses and strains, viscoplastic shear strains, yield function values, and the unbalanced force ratio — all of which can be used for post-processing and visualization with `plot_fem_results()`.

## Shear Strength Reduction Method (SSRM)

The Shear Strength Reduction Method (SSRM) represents the most widely adopted approach for determining factors of 
safety in finite element slope stability analysis, providing a rigorous and theoretically sound alternative to traditional limit equilibrium methods (Matsui & San, 1992; Griffiths & Lane, 1999). This method elegantly bridges the gap between the limit equilibrium concept of factor of safety and the stress-strain framework of finite element analysis. The fundamental principle underlying SSRM is conceptually straightforward yet mathematically sophisticated. Rather than assuming a failure surface and checking equilibrium conditions, SSRM systematically reduces the soil's shear strength parameters until the finite element system can no longer maintain equilibrium under the applied loading conditions. The reduction factor required to bring the slope to the brink of failure represents the factor of safety, defined consistently with traditional limit equilibrium approaches.

### Methodology

The SSRM procedure follows a systematic approach that progressively weakens the soil until failure occurs. The process begins by reducing both cohesion and friction angle by a trial factor $F$ according to the relationships:

>>$c_r = \dfrac{c}{F}$<br>
$\tan \phi_r = \dfrac{\tan \phi}{F}$

This reduction scheme ensures that both components of shear strength are diminished proportionally, maintaining the fundamental character of the Mohr-Coulomb failure criterion while systematically reducing the available resistance to shear failure. The choice to reduce the tangent of the friction angle rather than the friction angle itself ensures mathematical consistency and avoids complications that arise when the friction angle approaches zero.

With the reduced strength parameters, the finite element system is solved using the same equilibrium equations and constitutive relationships employed in conventional stress analysis. However, as the reduction factor increases, the soil's capacity to resist the applied gravitational and external loads diminishes, leading to progressively larger deformations and increasing numbers of elements reaching the yield condition.

The iterative nature of SSRM requires careful monitoring of the solution behavior to identify the onset of failure. Convergence characteristics provide the primary indicator of impending failure, as the finite element system transitions from stable equilibrium solutions to unstable behavior characterized by rapidly increasing displacements and failure to achieve force equilibrium.

The critical factor of safety is determined when the iterative solution process fails to converge within acceptable tolerances, indicating that the reduced strength parameters are insufficient to maintain equilibrium under the applied loading conditions. This point represents the transition from stable to unstable behavior and corresponds to the classical definition of factor of safety as the ratio of available strength to required strength for equilibrium.

### SSRM Failure Criteria

Determining the critical factor of safety in the SSRM requires a criterion to distinguish stable from unstable configurations as the strength reduction factor $F$ increases. The choice of failure criterion can significantly influence the computed factor of safety. XSLOPE implements four failure criteria, selectable via the `failure_criterion` parameter in `solve_ssrm()`:

#### 1. Non-Convergence (`"non_convergence"`)

This is the classical approach of Griffiths & Lane (1999). The viscoplastic iteration loop is run with no displacement limit; failure is identified purely by whether the solver converges within the maximum number of iterations. The SSRM bisection brackets the transition between convergence (stable) and non-convergence (failure).

This criterion is the most theoretically pure — it relies entirely on the mathematical behavior of the equilibrium equations. However, it can be sensitive to the convergence tolerance and maximum iteration count. With the elastic-relative convergence criterion described above, the non-convergence transition is sharper than with the conventional relative criterion, because false convergence at high displacement levels is eliminated.

#### 2. Displacement Limit (`"displacement_limit"`)

This criterion adds a physical failure check on top of the convergence criterion. During the viscoplastic iteration, the viscoplastic (plastic) displacement is computed as the difference between the current total displacement and the initial elastic displacement:

>>$\{U\}_{vp} = \{U\}_{total} - \{U\}_{elastic}$

If the maximum nodal VP displacement magnitude exceeds a specified fraction of the mesh height, the slope is declared as failed regardless of whether the convergence criterion is satisfied:

>>$\max_i ||\{U\}_{vp,i}|| > \alpha \cdot H_{mesh}$

where $\alpha$ is the displacement limit factor (default 0.1, i.e., 10% of mesh height) and $H_{mesh}$ is the vertical extent of the mesh. This prevents the solver from reporting convergence when plastic displacements have grown to physically unreasonable levels.

The displacement limit is controlled by the `max_disp_factor` parameter. The default value of 0.1 (10% of mesh height) is empirical but has physical meaning: a slope with plastic displacements exceeding 10% of its height has undergone significant deformation indicative of failure.

#### 3. Displacement Catastrophe (`"displacement_increase"`)

This criterion implements the displacement mutation method described by Sun et al. (2021). Rather than using a fixed displacement threshold, it detects the **sudden acceleration** in displacement growth as the strength reduction factor increases — a hallmark of catastrophic failure.

The procedure has three phases:

>**Phase 1 — Coarse sweep:** The solver runs at $n$ evenly-spaced values of $F$ from $F_{min}$ to $F_{max}$, recording the maximum VP displacement at each. A generous displacement cap (50% of mesh height) is used solely as an early termination to avoid wasting iterations on obviously failed cases.

>**Phase 2 — Catastrophe detection:** The ratio of VP displacement between consecutive $F$ values is computed. The interval with the largest ratio identifies where the displacement growth accelerates most rapidly — the catastrophe point.

>**Phase 3 — Refinement:** The catastrophe interval is refined by bisection, using displacement ratios to determine which half contains the sharper transition, until the interval width falls below the specified tolerance.

This criterion is self-calibrating in that it does not require an absolute displacement threshold. However, it requires the displacement-versus-$F$ curve to exhibit a sufficiently sharp inflection. For problems where the displacement grows smoothly (gradual transition rather than sudden catastrophe), the method may identify the acceleration zone at a higher $F$ than other criteria.

#### 4. Unbalanced Force Ratio (`"unbalanced_force"`)

This criterion is inspired by the unbalanced force ratio used in FLAC (Itasca, 2019) as the primary indicator of mechanical equilibrium. The unbalanced force ratio (UFR) measures the magnitude of the viscoplastic body load corrections relative to the total applied gravity load:

>>$\text{UFR} = \dfrac{||\{F\}_{vp}||}{||\{F\}_{gravity}||}$

where $\{F\}_{vp} = \{F\}_{loads} - \{F\}_{gravity}$ represents the accumulated body load corrections from all yielding Gauss points. A large UFR indicates that significant force redistribution is still occurring — the soil with reduced strength cannot achieve equilibrium under the applied loads.

The UFR method works by establishing a baseline UFR at $F_{min}$ (the unreduced slope), then bisecting to find the $F$ where the UFR exceeds a specified multiple of the baseline. The `ufr_threshold` parameter (default 2.0) defines this multiplier: failure is declared when $\text{UFR} > \text{ufr\_threshold} \times \text{UFR}_{baseline}$.

Unlike the displacement-based criteria, the UFR criterion directly measures force equilibrium rather than deformation. However, the multiplier threshold is itself a parameter that must be chosen, analogous to the displacement limit factor.

#### Comparison of Failure Criteria

The following table summarizes the characteristics of each criterion:

| Criterion | Basis | Arbitrary Parameter | Strengths | Limitations |
|---|---|---|---|---|
| Non-convergence | Mathematical convergence | Tolerance, max iterations | Theoretically pure | Sensitive to convergence definition |
| Displacement limit | Physical deformation | `max_disp_factor` (default 0.1) | Simple, physically intuitive | Threshold is empirical |
| Displacement catastrophe | Rate of displacement growth | `n_sweep` points | Self-calibrating, no absolute threshold | Requires sharp inflection in displacement curve |
| Unbalanced force ratio | Force equilibrium | `ufr_threshold` (default 2.0) | Measures equilibrium directly | Multiplier must be chosen |

The **non-convergence** criterion is the default in XSLOPE. It is the most theoretically rigorous approach — based directly on the foundational work of Griffiths & Lane (1999) — and requires no arbitrary parameters beyond the convergence tolerance and maximum iteration count. It is important to recognize that FEM-SSRM and limit equilibrium methods are fundamentally different approaches to slope stability, and some difference in computed factors of safety is expected. The alternative criteria are provided for users who wish to explore the sensitivity of results to the failure definition. Users are encouraged to compare multiple criteria to understand this sensitivity for their specific problem.

### The `solve_ssrm()` Function

The SSRM procedure is implemented in the `solve_ssrm()` function, which repeatedly calls `solve_fem()` at different strength reduction factors to bracket and refine the critical factor of safety.

```python
from xslope.fem import solve_ssrm

# Find critical factor of safety using SSRM
result = solve_ssrm(
    fem_data,
    F_min=1.0,
    F_max=2.0,
    tolerance=0.05,
    debug_level=1,
    failure_criterion="non_convergence"
)

if result['converged']:
    print(f"Factor of Safety: {result['FS']:.2f}")
    print(f"SSRM iterations: {result['iterations_ssrm']}")
    print(f"Final interval: {result['final_interval']}")
else:
    print(f"SSRM failed: {result.get('error', 'Unknown error')}")
```

The key parameters of `solve_ssrm()` are:

>- **`F_min`** (default 1.0): Lower bound for the bisection search. The slope must be stable (converge) at this reduction factor.<br>
>- **`F_max`** (default 2.0): Upper bound for the bisection search. The slope should be unstable (not converge) at this reduction factor.<br>
>- **`tolerance`** (default 0.05): Bisection stops when $F_{right} - F_{left} <$ tolerance. The critical FS is reported as $F_{left}$ (the last stable value).<br>
>- **`failure_criterion`** (default `"non_convergence"`): Selects the failure criterion — `"non_convergence"`, `"displacement_limit"`, `"displacement_increase"`, or `"unbalanced_force"` as described above.<br>
>- **`max_disp_factor`** (default 0.1): Displacement limit fraction, used with the `"displacement_limit"` criterion.<br>
>- **`n_sweep`** (default 10): Number of coarse sweep points for the `"displacement_increase"` criterion.<br>
>- **`ufr_threshold`** (default 2.0): UFR multiplier for the `"unbalanced_force"` criterion.<br>
>- **`convergence_tol`** (default $10^{-3}$) and **`max_iterations`** (default 500): Passed through to `solve_fem()` for each trial.

The returned result dictionary contains the critical factor of safety (`FS`), the last converged `solve_fem()` solution (`last_solution`), the final bisection interval, and the number of SSRM iterations. The `last_solution` can be passed directly to `plot_fem_results()` for visualization of the failure mechanism at the critical state.

## Element Type Selection and Volumetric Locking

### The Problem: Volumetric Locking in Low-Order Elements

A critical consideration in finite element slope stability analysis is the choice of element type. Low-order elements — specifically 3-node linear triangles (tri3) and 4-node bilinear quadrilaterals (quad4) — suffer from a well-known numerical pathology called **volumetric locking** that causes them to significantly overestimate the factor of safety.

Volumetric locking occurs because plastic deformation under the Mohr-Coulomb criterion (particularly with a non-associated flow rule where the dilation angle $\psi = 0$) produces nearly incompressible plastic strains — the material deforms in shear without changing volume. Low-order elements have too few degrees of freedom to simultaneously satisfy this incompressibility constraint and represent the displacement field accurately. The result is an artificially stiff response: the elements resist plastic deformation more than they should, requiring a larger strength reduction factor before failure occurs. This manifests as an unconservative (too high) factor of safety.

The severity of locking depends on the element formulation:

- **Constant-strain triangles (tri3)** have a single integration point and only 6 displacement DOFs per element. The constant strain field cannot represent the complex shear deformation patterns that develop in plastic zones. This is the most severely affected element type.

- **Bilinear quadrilaterals (quad4)** with full 2×2 integration have 8 displacement DOFs but still exhibit significant locking because the bilinear interpolation cannot adequately represent incompressible deformation modes.

- **Reduced integration** (using fewer Gauss points than required for exact integration) can alleviate locking in quadrilateral elements. The 8-node quadrilateral (quad8) with 2×2 reduced integration is a particularly effective combination that avoids locking while maintaining accuracy.

### Recommended Elements: Quadratic Formulations

Quadratic elements — 6-node triangles (tri6), 8-node quadrilaterals (quad8), and 9-node quadrilaterals (quad9) — provide sufficient degrees of freedom to represent incompressible plastic deformation without artificial stiffness. These element types should always be used for slope stability analysis with the SSRM.

The following table shows SSRM results for the Griffiths & Lane (1999) Example 1 benchmark problem (homogeneous slope with $c/\gamma H = 0.05$, $\phi = 20°$, slope angle 26.57°) using each element type with a target mesh size of 5. The expected factor of safety is approximately 1.40 (Griffiths & Lane report FS = 1.4 by FEM, Spencer's method gives FS = 1.376).

| Element Type | Nodes per Element | SSRM Factor of Safety | Error vs. Reference | Recommendation |
|:---:|:---:|:---:|:---:|:---|
| tri3 | 3 | 1.70 | +21% | Not recommended — severe locking |
| quad4 | 4 | 1.56 | +11% | Not recommended — significant locking |
| **tri6** | **6** | **1.41** | **< 1%** | **Recommended** |
| **quad8** | **8** | **1.41** | **< 1%** | **Recommended (default)** |
| **quad9** | **9** | **1.41** | **< 1%** | **Recommended** |

The three quadratic element types all converge to the same factor of safety (FS = 1.41), consistent with the published reference solution. The low-order elements overestimate FS by 11–21%, which would be unconservative and potentially dangerous in practice.

### Practical Guidance

For slope stability analysis in XSLOPE:

1. **Always use quadratic elements** (tri6, quad8, or quad9) for SSRM analysis. The default element type in XSLOPE is quad8, which provides the best balance of accuracy and efficiency through reduced integration.

2. **Do not use tri3 or quad4 for computing factors of safety.** While these elements may be adequate for purely elastic problems or for qualitative assessment of stress distributions, they produce unconservative factors of safety in elastic-plastic slope stability analysis.

3. **quad8 with reduced integration** is the preferred choice following Griffiths & Lane (1999). The 2×2 Gauss integration (instead of the full 3×3 required for exact integration of the quadratic shape functions) effectively eliminates volumetric locking while providing accurate stress and strain fields.

4. **tri6 elements** are useful when the mesh must conform to complex geometries where quadrilateral elements would become highly distorted. The 3-point triangle integration rule provides sufficient accuracy for the quadratic interpolation.

5. **quad9 elements** with full 3×3 integration also produce correct results and may be preferred for problems requiring higher stress accuracy, though at a computational cost proportional to the additional Gauss points and the extra center node per element.

## Seismic Forces

For both limit equilibrium and finite element slope stability analysis in XSLOPE, seismic loading is simulated using the pseudo-static method, which is a widely accepted approach for incorporating seismic effects into slope stability assessments. This method simplifies the complex dynamic response of soil during earthquakes by representing seismic loading as equivalent static forces applied to the slope mass. The pseudo-static approach assumes that the earthquake ground acceleration can be represented by a constant horizontal acceleration applied throughout the slope mass. This acceleration generates inertial forces that act on every element of soil, creating additional driving forces that tend to destabilize the slope.

In the finite element formulation, seismic loading is incorporated by modifying the body force vector to include both gravitational and seismic components:

>>$\{b\}_{total} = \{b\}_{gravity} + \{b\}_{seismic}$

For a horizontal seismic coefficient $k$, representing the ratio of horizontal acceleration to gravitational acceleration, the seismic body forces are:

>>$b_{x,seismic} = k \gamma$<br>
$b_{y,seismic} = 0$

The direction of the horizontal seismic force is chosen to maximize the driving forces that promote slope failure. For typical left-facing slopes, this corresponds to a leftward (negative x-direction) seismic acceleration that increases the shear stresses along potential failure surfaces. For a right-facing slope, the force acts in the positive x-direction. For some problems, such as dams or levees, the slope geometry can include both right- and left-facing slopes. With the limit equilibrium method, XSLOPE determines if a slope is right- or left-facing based on the location and geometry of the failure surface. But with the finite element method, both sides of a dam or levee are analyzed at the same time. Therefore, the user should adjust the sign of the seismic coefficient (k) in the Excel input template to indicate the direction of the psuedostatic seismic force. For the left-facing slopes, k should be entered as a negative value indicating a force in the negative x direction, and for right-facing slopes, k should be entered as a positive value. 

The equilibrium equations are modified to include seismic body forces as follows:

>>$\dfrac{\partial \sigma_x}{\partial x} + \dfrac{\partial \tau_{xy}}{\partial y} + b_x + k\gamma = 0$

>>$\dfrac{\partial \tau_{xy}}{\partial x} + \dfrac{\partial \sigma_y}{\partial y} + b_y = 0$

where $k\gamma$ represents the additional horizontal body force due to seismic loading, with $k$ being the seismic coefficient and $\gamma$ the unit weight of the soil. Since $b_x = 0$ and $b_y = \gamma$, the equations simplify to:

>>$\dfrac{\partial \sigma_x}{\partial x} + \dfrac{\partial \tau_{xy}}{\partial y} + k\gamma = 0$

>>$\dfrac{\partial \tau_{xy}}{\partial x} + \dfrac{\partial \sigma_y}{\partial y} + \gamma = 0$

The seismic body forces are incorporated into the global force vector through element-level integration:

>>$\{F\}_{seismic} = \sum_{e=1}^{N_{elem}} \int_{A_e} [N]^T \{b\}_{seismic} \, dA$

For linear triangular elements, this integration yields nodal forces that distribute the seismic loading according to the element shape functions. The seismic forces are added to the gravitational body forces and any external applied loads to form the complete loading vector.

## Structural Elements

XSLOPE supports two types of one-dimensional structural elements embedded within the two-dimensional soil mesh: **truss elements** for flexible reinforcement and **beam elements** for piles and concrete piers. Both element types share nodes with the surrounding 2D soil elements and participate in the viscoplastic iteration loop through body-force corrections.

- **[Soil Reinforcement](reinforcement.md)**: Flexible reinforcement (geotextiles, soil nails, ground anchors) modeled as tension-only truss elements with axial stiffness $EA/L$. Includes failure modes (pullout, peak-residual, complete failure), wished-in-place analysis considerations, and typical material properties.

- **[Piles and Concrete Piers](piles.md)**: Rigid structural elements modeled as beam elements with both axial stiffness ($EA/L$) and lateral bending stiffness ($3EI/L^3$). Unlike reinforcement, piles carry both tension and compression. Rotational DOFs are eliminated through static condensation to maintain compatibility with the 2-DOF-per-node soil mesh.

In both cases, the structural element properties are **not reduced** during SSRM strength reduction — only soil $c$ and $\tan\varphi$ are reduced. The resulting factor of safety represents the margin of safety in the soil strength, given the structural elements as-designed.

## Integration with XSLOPE Framework

The integration of finite element capabilities into the XSLOPE framework represents leverages the established infrastructure for limit equilibrium and seepage analysis while extending the analysis 
capabilities to include 
rigorous 
stress-strain based slope stability assessment. 

### Input Integration

The finite element implementation uses XSLOPE's  Excel-based input system. The slope geometry definitions, material 
property specifications, and boundary condition inputs that currently support limit equilibrium analysis provide an 
excellent foundation for finite element modeling. Material property definitions can leverage the existing cohesion 
$c$, friction angle $\phi$, and unit weight $\gamma$ specifications that are  in the XSLOPE input system. Additional 
parameters required for finite element analysis, such as Young's modulus $E$ and Poisson's ratio $\nu$, are included 
in the material property table in the Material tab.

The existing seepage analysis infrastructure provides a particularly valuable foundation for coupled seepage-stability finite element analysis. The current mesh generation capabilities and pore pressure calculation tools can be directly utilized to provide the groundwater input conditions required for effective stress analysis in the finite element system.

In the reinforcement table, some of the properties for each reinforcement line such as the pullout length $L_p$ and 
maximum tensile force $T_{max}$ are used for both limit equilibrium and finite element analysis. The additional 
parameters required for finite element analysis, such as the modulus of elasticity $E$ and cross-sectional area $A$, are also included in the reinforcement table. This allows the same reinforcement definitions to be used seamlessly across both analysis methods.

### Pore Pressure Options

Pore pressures reduce the effective stress in the soil, which in turn reduces the available shear strength. XSLOPE supports three pore pressure options for each material zone, specified via the `u` column in the material property table of the input template:

**None** (`u = "none"`): No pore pressure is applied. The yield function uses total stress. This is appropriate for dry slopes or total stress analyses with undrained shear strength parameters.

**Piezometric Line** (`u = "piezo"`): A piezometric surface is defined in the input template as a series of coordinate points. The pore pressure at any point below the piezometric surface is computed as:

>>$u = \gamma_w (z_{piezo} - z)$

where $\gamma_w$ is the unit weight of water, $z_{piezo}$ is the elevation of the piezometric surface directly above the point, and $z$ is the elevation of the point. For points above the piezometric surface, $u = 0$ (suction is conservatively ignored).

**Seepage Solution** (`u = "seep"`): Pore pressures are obtained from a prior finite element seepage analysis performed with the `seep.py` module. The seepage analysis solves the groundwater flow equation on a triangular mesh, producing pore pressure values at all nodes:

>>$u = \gamma_w (h - z)$

where $h$ is the hydraulic head from the seepage solution and $z$ is the elevation coordinate. These nodal pore pressures are stored in the slope data and transferred to the structural mesh during `build_fem_data()`. Negative pore pressures (suction above the phreatic surface) are clamped to zero.

For both the piezometric and seepage options, pore pressures are precomputed at each Gauss point during `build_fem_data()` using the element shape functions to interpolate from nodal values (seep) or by computing the physical coordinates of each Gauss point and projecting onto the piezometric surface (piezo). These precomputed values are then used directly in the effective stress yield check during the viscoplastic iteration, avoiding repeated interpolation at each iteration step.

## Visualization of Results

The `plot_fem_results()` function provides a flexible interface for visualizing FEM solutions. It accepts a `plot_type` parameter that specifies one or more plot types to display as vertically stacked subplots. The available plot types are:

| Plot Type | Description |
|-----------|-------------|
| `deformation` | Deformed mesh overlay on the original mesh. The original mesh is shown in light gray and the deformed shape in blue. Viscoplastic displacements (total minus elastic) are used when available, so the plot shows the failure mechanism rather than gravity settlement. The `deform_scale` parameter amplifies displacements for visibility and is auto-calculated by default so the maximum deformation is approximately 10% of the mesh height. |
| `shear_strain` | Viscoplastic maximum shear strain contours. This is the most important plot for identifying the failure mechanism — high shear strain concentrations reveal the failure surface without any prior assumption about its shape or location. Falls back to total shear strain if viscoplastic data is not available. |
| `displace_vector` | Displacement vectors at corner nodes. Like the deformed mesh plot, viscoplastic displacements are used when available to show the failure mechanism. Vectors below a threshold fraction of the maximum displacement are hidden to reduce clutter. |
| `displace_mag` | Displacement magnitude contours using the viridis colormap. Shows total displacement magnitude at each node as filled contours. |
| `stress` | Von Mises stress contours with yielded elements highlighted. Useful for identifying stress concentrations and the extent of the plastic zone. |
| `strain` | Von Mises equivalent strain contours computed from total strains. |
| `yield` | Mohr-Coulomb yield function contours. Positive values indicate yielding/failure, negative values indicate an elastic state. |

The default combination is `['deformation', 'shear_strain', 'displace_vector']`, which provides a comprehensive view of the failure mechanism. The following example shows the SSRM results for the Griffiths & Lane (1999) homogeneous slope from the [FEM Samples](samples.md) page (Problem 1). The top subplot shows the deformed mesh at the last converged solution, and the bottom subplot shows the viscoplastic shear strain concentration revealing the circular failure mechanism:

![griffiths1_results.png](images/griffiths1_results.png){width=1000}

Additional options control the appearance of all plot types:

- `show_mesh` — show or hide mesh lines (default: `True`)
- `show_reinforcement` — show or hide reinforcement and pile elements (default: `True`)
- `label_elements` — show element ID labels at centroids (default: `False`)
- `figsize` — figure size as `(width, height)` in inches (default: `(12, 8)`)
- `save_png` — save the figure to a PNG file (default: `False`)
- `dpi` — resolution for saved PNG (default: `300`)

A typical call for an SSRM analysis uses the default plot types:

```python
plot_fem_results(fem_data, result['last_solution'],
                 plot_type=['deformation', 'shear_strain', 'displace_vector'],
                 save_png=True)
```

A single plot type can also be specified as a string rather than a list:

```python
plot_fem_results(fem_data, solution, plot_type='shear_strain')
```

## References

Dawson, E.M., Roth, W.H., & Drescher, A. (1999). Slope stability analysis by strength reduction. *Géotechnique*, 49(6), 835-840.

Duncan, J.M. (1996). State-of-the-art: Limit equilibrium and finite element analysis of slopes. *Journal of Geotechnical Engineering*, 122(7), 577-596.

Duncan, J.M., & Wright, S.G. (2005). *Soil Strength and Slope Stability*. John Wiley & Sons.

Dyson, A.P., & Tolooiyan, A. (2018). Comparative approaches to probabilistic finite element methods for slope stability analysis. *Innovative Infrastructure Solutions*, 3(1), 1-11.

Griffiths, D.V., & Lane, P.A. (1999). Slope stability analysis by finite elements. *Géotechnique*, 49(3), 387-403.

Itasca Consulting Group. (2019). *FLAC — Fast Lagrangian Analysis of Continua, Version 8.1, User's Guide*. Itasca Consulting Group, Inc., Minneapolis, Minnesota.

Matsui, T., & San, K.C. (1992). Finite element slope stability analysis by shear strength reduction technique. *Soils and Foundations*, 32(1), 59-70.

Smith, I.M., & Griffiths, D.V. (2004). *Programming the Finite Element Method* (4th ed.). John Wiley & Sons.

Sun, G., Lin, S., Jiang, W., & Yang, Y. (2021). A simplified solution for determining the factor of safety of a slope reinforced with piles based on the shear strength reduction method. *Bulletin of Engineering Geology and the Environment*, 80, 7719-7730.

Zheng, H., Liu, D.F., & Li, C.G. (2005). Slope stability analysis based on elasto‐plastic finite element method. *International Journal for Numerical Methods in Engineering*, 64(14), 1871-1888.