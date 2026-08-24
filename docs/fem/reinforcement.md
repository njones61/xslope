# Soil Reinforcement in Finite Element Analysis

The integration of soil reinforcement elements such as geotextiles, soil nails, and ground anchors into finite element slope stability analysis represents a significant advancement in modeling stabilized slopes. These reinforcement systems fundamentally alter the stress distribution and failure mechanisms within slopes, requiring sophisticated modeling approaches to capture their beneficial effects accurately (Duncan & Wright, 2005).

![reinf_layers.png](images/reinf_layers.png)

The modeling of reinforced slopes presents unique challenges because the reinforcement elements typically have dramatically different mechanical properties compared to the surrounding soil. Reinforcement elements are usually much stiffer in tension and often have negligible compressive strength, creating a highly anisotropic composite material that requires specialized finite element formulations.

## Truss Element Approach

While there are numerous ways to simulate soil reinforcement in the finite element method including the equivalent
force method and interface element modeling, the most straightforward method involves representing reinforcement elements as one-dimensional truss elements embedded within the two-dimensional soil continuum. These truss elements are characterized by their axial stiffness $EA/L$, where $E$ is the elastic modulus of the reinforcement material, $A$ is the cross-sectional area, and $L$ is the element length.

This approach is particularly effective for modeling geosynthetic reinforcement, soil nails, and tie-back anchors.
The truss elements can only carry tension loads up to a specified tensile strength limit $T_{max}$, beyond which
they either yield plastically to a residual strength $T_{res}$, or fail completely. The inability to carry compression loads
accurately
reflects the behavior of flexible reinforcement materials like geotextiles and ensures that the reinforcement cannot resist compressive buckling. The truss elements are oriented along the centerline of the physical reinforcement and connected to the surrounding soil elements through shared nodes. This connection ensures that the reinforcement participates in the overall deformation pattern of the slope while contributing its tensile resistance to improve stability.

Truss elements are incorporated into the XSLOPE finite element mesh by passing the geometry of the reinforcement
lines from the input template to the mesh generation process. The reinforcement lines are discretized into multiple
truss elements based on the specified mesh density (target_size), or on the **1D element size** on the main sheet
where a model states one — the element size along the reinforcement and pile lines, which refines the truss elements
and the soil sharing their nodes together. The 1D elements are fully integrated with the 2D
elements - each 1D element corresponds to the edge of two adjacent 2D elements and both the 1D and 2D elements share
the same nodes. The 1D elements have their own set of material properties corresponding to the properties of the
corresponding reinforcement lines input by the user and include $T_{max}$, $T_{res}$, $E$, and cross-sectional area
$A$.

Only the two end nodes of each 1D element are used for computing truss stiffness and axial forces. When the mesh uses
quadratic (8-node) 2D elements, the 1D elements may have a mid-side node shared with the adjacent 2D elements. This
mid-node participates in the 2D element shape functions for displacement interpolation but is ignored for the truss
element formulation. A 2-node linear truss element is exact for a prismatic bar with constant axial stiffness, so
the quadratic interpolation adds no physical fidelity for the 1D element.

The meshing algorithms used in XSLOPE, including the integration of 1D and 2D elements for problems involving soil
reinforcement are documented in the [Mesh Generation](mesh.md) page.

## Mathematical Formulation

**Truss Element Stiffness Matrix:** Each 1D truss element contributes to the global stiffness matrix through its element stiffness matrix. For a truss element with nodes $i$ and $j$, the element stiffness matrix in local coordinates is:

>>$[K_e]_{local} = \dfrac{AE}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$

where $A$ is the cross-sectional area, $E$ is the elastic modulus, and $L$ is the element length.

**Coordinate Transformation:** The local stiffness matrix must be transformed to global coordinates using the transformation matrix $[R]$:

>>$[K_e]_{global} = [R]^T [K_e]_{local} [R]$

The transformation is built from $\psi$, the inclination of the reinforcement line to the horizontal — the same
angle the LEM formulation uses for the direction of an axial reinforcement force:

>>$[R] = \begin{bmatrix} \cos\psi & \sin\psi & 0 & 0 \\ 0 & 0 & \cos\psi & \sin\psi \end{bmatrix}$

**Assembly Process:** The global stiffness matrix combines contributions from both 2D soil elements and 1D truss elements:

>>$[K]_{global} = \sum_{soil} [K_e]_{soil} + \sum_{truss} [K_e]_{truss}$

**Force Vector Assembly:** The global force vector includes both soil body forces and any applied forces on reinforcement:

>>$\{F\}_{global} = \{F\}_{soil} + \{F\}_{reinforcement}$

## Force Behavior and Failure Modes

The forces and failure modes in 1D truss elements are analyzed in an iterative fashion. Each truss element
has a maximum allowable tensile capacity $T_{allow}$, derived from the user-specified reinforcement parameters,
and optionally a residual tensile capacity $T_{res}$. The axial force in each element is calculated as
$T = (AE/L) \cdot \delta$, where $\delta$ is the element elongation (the component of relative nodal displacement
along the element axis).

The global stiffness matrix carries each bar's *full elastic* stiffness, so $K u$ always contains the uncapped
elastic force. The capacity is imposed the same way plasticity is imposed on the soil — through a viscoplastic
body load equal to the part of the elastic force the element **cannot** carry:

>>$f_{body} = (T - T_{true}) \cdot [-\cos\psi,\; -\sin\psi,\; +\cos\psi,\; +\sin\psi]$

where $T_{true}$ is the force the bar can actually deliver: the elastic $T$ clipped into $[0, T_{cap}]$. Because
equilibrium is solved as $K u - f_{body}$, this leaves exactly $T_{true}$ in the bar. (The sign matters. Adding
the *opposite* correction makes an overloaded bar carry $2T - T_{cap}$ — it gets **stiffer** the more it is
overloaded, an "anti-cap" under which a reinforced slope can never be driven to failure and the SSR factor is
insensitive to $T_{allow}$ altogether.)

![reinf_bar_law.png](images/reinf_bar_law.png)

Which of the three post-peak branches a bar follows is decided entirely by the $T_{res}$ column of its
reinforcement line.

**Elastic-Perfectly-Plastic Model (the default):**

If $T_{res}$ is left **blank**, the bar yields at $T_{allow}$ and holds it indefinitely while the surrounding soil
keeps straining. This is the default because it is what the mainstream FEM codes do (PLAXIS geogrids and anchors
are elastoplastic with a maximum axial force), and it is what published reinforced-slope analyses assume.

A blank $T_{res}$ means *no post-peak drop* — it does **not** mean zero.

**Peak-Residual Model:**

Entering a value for $T_{res}$ turns on post-peak behavior: an element that yields drops from $T_{allow}$ to its
residual capacity, which is $T_{res}$ or the capacity its embedment can develop, whichever is smaller. Appropriate
for ductile materials where the published capacity is a peak rather than a plateau; typical residual ratios for
geosynthetics are $T_{res}/T_{allow} = 0.3-0.7$.

The drop is decided **only on a converged equilibrium state**, never inside the viscoplastic iteration. This
matters: the first iterate of a viscoplastic solve is the elastic predictor, whose bar forces overshoot badly
before the soil sheds load into them, so a mid-iteration trigger would condemn bars for a transient that never
physically existed, and the answer would depend on the path the solver happened to take. Instead the solver
converges with the bars capped at $T_{allow}$, then drops any bar whose elastic demand exceeded its capacity to
$T_{res}$ and re-solves. Shedding that load can push neighbors over, so the process repeats until the softened set
stops growing — a genuine progressive-failure fixed point, and one that is independent of the solution path.

**Complete (Brittle) Failure Model:**

Setting $T_{res} = 0$ explicitly is the brittle case: a yielding element ruptures and carries nothing afterwards.
Appropriate for brittle materials (some steel cables, fiber reinforcement).

!!! warning "Post-peak behavior makes the SSR factor mesh-sensitive"
    Once $T_{res} < T_{allow}$ actually engages, the reinforcement is strain-**softening**. A softening system in
    an unregularized continuum has no length scale to arrest localization, so the computed factor of safety can
    drift with mesh refinement instead of converging, and the SSRM bracket becomes less crisp. This is physics, not
    a numerical defect — but it means $T_{res}$ is best treated as a forensics/back-analysis parameter rather than
    a design default. Leave it blank unless you specifically intend to model post-peak strength loss.

**Pullout Failure Model:**

For each reinforcement line, it is assumed that the tension force in the reinforcement is zero at the two ends and increases linearly with distance along the line as frictional resistance between the reinforcement and the surrounding soil develops. Full tension force develops over a pullout distance, $L_p$.

- Pullout failure may occur in elements where the embedment length is less than the pullout length $L_p$<br>
- The available strength is limited by pullout resistance rather than material strength<br>
- For elements at distance $d$ from the reinforcement end where $d < L_p$:<br>
  >>$T_{available} = T_{allow} \times \frac{d}{L_p}$<br>
- Pullout is **perfectly plastic**. An element that reaches its embedment-limited capacity slips at that force and
  goes on carrying it: interface friction does not vanish once it has been overcome, so there is no drop to zero.
  This is the standard cable and geogrid treatment, and it is the same assumption the LEM envelope makes.<br>
- Pullout spreads along a line as elements near the ends reach their capacity one after another and shed the
  balance of the demand into the interior

**Tension-Only Behavior:**

Truss elements are restricted to carry only tension forces. This is implemented through body-force corrections within the viscoplastic iteration loop:

- After each iteration, the axial force in each element is computed from the current displacement field<br>
- If compression develops ($T < 0$), a corrective body force is applied that cancels the compressive force<br>
- The element remains in the stiffness matrix at full elastic stiffness; the correction enters through the load vector<br>
- This approach is consistent with the viscoplastic initial-stiffness method used for soil elements, where the stiffness matrix is factored once and all nonlinearity is driven through load corrections

## Integration with Viscoplastic Iteration

The 1D truss element nonlinearity (tension-only behavior, capacity limits, and failure) is handled through body-force
corrections within the [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387) viscoplastic iteration loop, using the same initial-stiffness approach
that governs the 2D soil elements. The key principle is that the global stiffness matrix $[K]$ is assembled once with
the full elastic stiffness of all elements (both 2D soil and 1D truss) and pre-factored for efficient repeated solves.
All nonlinear behavior is then driven entirely through corrections to the right-hand-side load vector.

**Assembly:** The global stiffness matrix includes contributions from both 2D soil elements and 1D truss elements:

>>$[K]_{global} = \sum_{soil} [K_e]_{soil} + \sum_{truss} [K_e]_{truss}$

This matrix is factored once (via sparse LU decomposition) and reused for all viscoplastic iterations.

**Iteration procedure:** At each viscoplastic iteration, after solving for the updated displacement field $\{u\}$:

1. For each 1D truss element, compute the axial force from the current displacements:
>>$\delta = (u_{x,j} - u_{x,i})\cos\psi + (u_{y,j} - u_{y,i})\sin\psi$
>>$T = \dfrac{AE}{L} \cdot \delta$

2. Determine whether a correction is needed:
>>- If $T < 0$ (compression): set $\Delta T = -T$ (cancel the compressive force entirely)
>>- If $T > T_{allow}$ and the element has not previously failed: mark the element as failed and set $\Delta T = T_{res} - T$
>>- If $T > T_{res}$ and the element has previously failed: set $\Delta T = T_{res} - T$
>>- Otherwise: no correction ($\Delta T = 0$)

3. Convert the axial force correction to equivalent nodal forces and add to the load vector:
>>$\{F\}_{correction} = \Delta T \begin{Bmatrix} -\cos\psi \\ -\sin\psi \\ \cos\psi \\ \sin\psi \end{Bmatrix}$

These corrections are added to the same load vector that receives the soil viscoplastic strain corrections ($[B]^T [D] \{\varepsilon_{vp}\}$). The factored stiffness matrix then solves the corrected system, and the process repeats until convergence.

**Failure irreversibility:** Once an element is marked as failed (having exceeded $T_{allow}$), it remains failed for all subsequent iterations within that analysis. Its effective capacity permanently drops from $T_{allow}$ to the residual assigned to it. This models the irreversible nature of material yielding.

## Strength Reduction and Reinforcement

In the Shear Strength Reduction Method (SSRM), only the soil shear strength parameters $c$ and $\tan\phi$ are reduced
by the strength reduction factor $F$. The reinforcement properties ($T_{max}$, $T_{res}$, $E$, $A$) are held constant
throughout the SSRM bisection. The resulting factor of safety represents the margin of safety in the soil strength,
given the structural reinforcement as-designed. This follows standard practice (Duncan & Wright, 2005) where
reinforcement capacity is treated as a structural property independent of the soil strength reduction.

This also means that the truss element contributions to the global stiffness matrix do not change between SSRM
bisection steps — only the soil yield parameters change — so the same pre-factored stiffness matrix approach remains
efficient.

## Reinforcement Line Input Parameters and Element Properties

In the Excel input template used by XSLOPE, the user can define up to 20 reinforcement lines by entering the reinforcement line geometry and properties into the lines of a table. Each row of the table includes the following:

| Item | Description |
|:----:|-------------|
| x1, y1 | The x and y coordinates of the left end of the line |
| x2, y2 | The x and y coordinates of the right end of the line |
| Tmax | Maximum allowable tensile force |
| Tres | Residual tensile force the reinforcement retains after it ruptures, capped by the capacity its embedment can develop. **Leave blank for no post-peak drop** (elastic-perfectly-plastic — the usual choice, and the default). An explicit `0` means brittle rupture. Used by the FEM only. |
| Lp1  | The pullout length on the left side |
| Lp2  | The pullout length on the right side |
| Adhesion | Soil-reinforcement interface adhesion. Filled together with Delta, it replaces Lp1/Lp2 with the overburden-dependent pullout law. Blank (with Delta blank) uses the pullout lengths. |
| Delta | Soil-reinforcement interface friction angle, degrees. |
| E    | The modulus of elasticity of reinforcement material  |
| Area | The cross-sectional area of the reinforcement material  |

The units for E and Area need to be compatible with each other and with the other weight and length units used. For
metric units, E should be in $kPa$ and Area should be in $m^2$. For English units, E should be in $psf$ and Area
should be in $ft^2$. Alternately, E could be in $psi$ as long as Area is in $in^2$.

### Element Discretization and Capacity Assignment

A separate pullout length (Lp) is used for each end since each end may be embedded in a separate soil with different shear resistance values. A line may instead state its interface strength through Adhesion and Delta, in which case the resistance follows the effective overburden along the line and the pullout lengths are not used.

During mesh generation, each reinforcement line is discretized into multiple truss elements based on the specified mesh density. The discretization process follows these steps:

1. **Material Property Assignment**: Each truss element along the line receives the same material properties:

>>Cross-sectional area: $Area$<br>
Elastic modulus: $E$<br>
Element stiffness: $K_e = AE/L$ (where L varies based on element length)<br>

2. **Tensile Capacity Assignment**: Each truss element is assigned an allowable $T_{allow}$ and residual $T_{res}$ tensile capacity. $T_{allow}$ is the capacity envelope at the element centroid — the **same** envelope the limit-equilibrium engine applies at a slip-surface crossing, evaluated by the same function, so the two engines cannot drift:

>>For an element whose centroid is at distances $d_1$ and $d_2$ from the two ends of a line of length $L$:
>>
>>$T_{allow} = \min\left(T_{max},\;\; T_{end1} + \displaystyle\int_0^{d_1} r,\;\; T_{end2} + \int_{L-d_2}^{L} r\right)$
>>
>>where $r$ is the pullout resistance per unit length. Under the development-length law $r = T_{max}/L_p$ at each end, and the integrals are the linear ramps $T_{max}d_1/L_{p1}$ and $T_{max}d_2/L_{p2}$. Under the overburden-dependent law (Adhesion and Delta both filled) $r = 2(a + \sigma'_v\tan\delta)$ varies along the line with the effective overburden, and $L_{p1}$/$L_{p2}$ are not read. Both laws, and the end anchorage capacities $T_{end}$, are set out in **[Soil Reinforcement in LEM](../lem/reinforcement.md#capacity-envelope)**.

>>
>>$T_{res} = \begin{cases}
\text{unset (no post-peak drop)} & \text{if } T_{res}\ \text{is blank in the input} \\
\min\left(T_{residual},\ T_{allow}\right) & \text{otherwise}
\end{cases}$

Elements near a free end therefore have reduced capacity — zero at the end itself, unless an anchorage capacity $T_{end}$ is entered there — while elements far enough from both ends carry the full design strength. The taper is the gradual development of pullout resistance through interface friction, and the minimum is taken over BOTH ends rather than the nearer one, which is the same answer wherever the two zones do not overlap and the correct one where they do. Each end has its own $L_{p1}$/$L_{p2}$ because each may be embedded in a different soil; under the overburden law the same variation comes from $\sigma'_v$ instead, without the soils having to be named.

The residual capacity is only assigned at all when the user has entered a $T_{res}$ for the line. Where post-peak behavior *is* switched on, two independent mechanisms can limit what an element retains, and the smaller of the two governs. Bond slip is perfectly plastic, so the embedment goes on developing $T_{allow}$ — the ramped envelope, end anchorage included — however far the bar is pulled. $T_{residual}$ is the rupture residual, a property of the reinforcement itself and not of its embedment. Beyond the ramps $T_{allow} = T_{max}$ and the element takes the user's residual strength; inside a ramp it takes whichever of the two is less.

### Axial Stiffness (EA)

The analysis depends only on the product $EA$ (the axial stiffness, sometimes called the tensile stiffness or $J$ in
geosynthetic specifications), not on $E$ and $A$ independently. Any combination of $E$ and $A$ that produces the same
$EA$ will give identical results. The axial stiffness controls how much the reinforcement must elongate before
mobilizing its tensile capacity:

>>$EA = \dfrac{T_{max}}{\varepsilon_{rupture}}$

where $\varepsilon_{rupture}$ is the strain at which the reinforcement reaches its ultimate tensile strength.

### Determining Reinforcement Line Pullout Lengths

The pullout length $L_p$ represents the distance from each end of the reinforcement over which the full tensile strength is mobilized. This variation captures the physical reality that pullout resistance must develop over a finite distance from the reinforcement ends through interface friction between the reinforcement and surrounding soil. This friction cannot be mobilized instantaneously but requires relative displacement to develop, creating the gradual strength buildup characteristic of all reinforcement systems. Pullout length can be estimated as follows:

**For Soil Nails:**
>>$L_p = \dfrac{T_{max}}{\alpha \pi D \sigma_n' \tan \phi_{interface}}$

where:<br>
>>$T_{max}$ = design tensile capacity of the nail<br>
$\alpha$ = surface roughness factor (0.5-1.0 for grouted nails)<br>
$D$ = effective nail diameter <br>
$\sigma_n'$ = average effective normal stress along the nail<br>
$\phi_{interface}$ = interface friction angle (typically 0.8-1.0 times soil friction angle)

**For Geotextiles:**
>>$L_p = \dfrac{T_{max}}{2 \alpha \sigma_n' \tan \phi_{interface}}$

where the factor of 2 accounts for friction on both sides of the geotextile.

These equations are a general guide that can be used to come up with reasonable estimates of Lp. Typical values are as follows:

|Reinforcement Type | Pullout Length $L_p$ (m) | Notes |
|-------------------|--------------------------|-------|
| **Soil Nails** | 1.5 - 3.0 | Depends on soil conditions and nail diameter |
| **Geotextiles** | 0.5 - 1.5 | Depends on normal stress and surface texture |
| **Geogrid** | 1.0 - 2.0 | Depends on aperture size and bearing resistance |

### Wished-in-Place Analysis and EA Selection

XSLOPE currently uses a **wished-in-place** approach: the entire slope (all soil layers and all reinforcement) is
assumed to exist in its final geometry, and gravity is applied in a single step. The reinforcement starts at zero
strain and zero force. Tension develops only through deformation that occurs during the gravity application and
subsequent SSRM strength reduction. This differs from reality, where reinforcement accumulates tension progressively
as each soil lift is placed during construction.

The wished-in-place approach is **conservative** — it underestimates the reinforcement contribution because it misses
the construction-induced pre-tension. However, for the reinforcement to provide meaningful stabilization in a
wished-in-place SSRM analysis, it must be stiff enough to develop significant force from the relatively small
displacements that occur near the incipient failure state. Low-stiffness reinforcement may undergo insufficient
strain during SSRM to mobilize its capacity, leading to an unrealistically low factor of safety.

Parametric studies show that the computed factor of safety increases with $EA$ and then plateaus above a threshold
stiffness, beyond which further increases in $EA$ have negligible effect. For typical reinforced slope geometries,
this plateau is reached at approximately $EA \approx 100$–$200 \times T_{max}$. Below approximately
$EA \approx 50 \times T_{max}$, the reinforcement may not mobilize enough force to significantly improve the factor
of safety in a wished-in-place analysis.

The following table provides recommended $EA$ values for wished-in-place SSRM analysis. These values are deliberately
at the stiffer end of the physical range for each material type, because the wished-in-place approach requires
sufficient stiffness to compensate for the absence of construction-induced pre-tension:

| Reinforcement Type | Recommended $EA/T_{max}$ | $\varepsilon_{rupture}$ | Notes |
|---|---|---|---|
| **Woven geotextiles** | $50$–$100$ | 1–2% | Use stiffer end for SSRM |
| **HDPE geogrids** | $50$–$100$ | 1–2% | Uniaxial, reinforcement grade |
| **PET geogrids** | $100$–$200$ | 0.5–1% | Higher stiffness than HDPE at same strength |
| **Steel strips** | $500$–$2{,}000$ | 0.05–0.2% | Very stiff, minimal elongation |
| **Soil nails (grouted)** | $1{,}000$–$5{,}000$ | 0.02–0.1% | Based on steel bar + grout composite |

### Typical E and Area Values

The following tables provide representative $E$ and $Area$ values for common reinforcement materials. These are
intended as starting points when manufacturer-specific data is not available. Any combination of $E$ and $Area$
producing the target $EA$ will give identical results.

**English Units:**

| Material | $E$ (psf) | $Area$ (ft$^2$/ft) | $EA$ (lb/ft) | $T_{max}$ (lb/ft) |
|---|---|---|---|---|
| Woven geotextile (light) | 500,000 | 0.01 | 5,000 | 100 |
| Woven geotextile (heavy) | 2,000,000 | 0.02 | 40,000 | 500 |
| HDPE geogrid | 2,000,000 | 0.03 | 60,000 | 800 |
| PET geogrid | 5,000,000 | 0.02 | 100,000 | 1,000 |
| Steel strip (galvanized) | 400,000,000 | 0.0003 | 120,000 | 5,000 |
| Soil nail (grouted, #8 bar) | 600,000,000 | 0.0006 | 360,000 | 15,000 |

**Metric Units:**

| Material | $E$ (kPa) | $Area$ (m$^2$/m) | $EA$ (kN/m) | $T_{max}$ (kN/m) |
|---|---|---|---|---|
| Woven geotextile (light) | 25,000 | 0.003 | 75 | 1.5 |
| Woven geotextile (heavy) | 100,000 | 0.006 | 600 | 7 |
| HDPE geogrid | 100,000 | 0.009 | 900 | 12 |
| PET geogrid | 250,000 | 0.006 | 1,500 | 15 |
| Steel strip (galvanized) | 20,000,000 | 0.00009 | 1,800 | 75 |
| Soil nail (grouted, #8 bar) | 30,000,000 | 0.0002 | 6,000 | 220 |

### Staged Construction Alternative

In practice, reinforced slopes and walls are built in lifts. Each soil layer is placed and compacted on top of
previously installed reinforcement, which develops tension in the reinforcement before the next lift is added. By the
time the slope reaches its final geometry, the lower reinforcement layers may have accumulated significant
construction-induced pre-tension.

A **staged construction analysis** models this process by activating soil layers and reinforcement elements
sequentially in the FEM, solving for equilibrium at each stage. The stress state — including locked-in reinforcement
tension — carries forward from each stage to the next. The SSRM analysis then begins from the end-of-construction
stress state rather than from a zero-strain condition.

Staged construction analysis is supported by all major commercial geotechnical FEM packages (PLAXIS, FLAC,
RS2/Phase2, SIGMA/W) and is the recommended approach when:

- The reinforcement has low axial stiffness (extensible geosynthetics)
- The slope is tall with many reinforcement layers
- Accurate prediction of reinforcement forces is important (not just the factor of safety)
- The analysis needs to match instrumented field measurements

For the wished-in-place approach currently used by XSLOPE, selecting $EA$ values at the stiffer end of the
recommended range (see table above) partially compensates for the absence of construction-induced pre-tension and
produces conservative but reasonable factors of safety. Staged construction analysis may be implemented in a future
version of XSLOPE.

## Inspecting the Results

The FEM results view colors each reinforcement element by the force it carries, which shows at a glance which
lines are working hardest. To read one line along its length, use the **1D Details…** button on that view's
toolbar. It opens a panel listing every reinforcement line and pile in the model with its utilization and a badge
colored by it — a reinforcement row also names the state the line is in — and draws the selected member's profiles
beside the list. Under the list is a map of the section with the selected member picked out, so a name in a list is
a place on the slope. The button is dimmed for a model with no reinforcement lines and no piles.

![Reinforcement detail for Line 4 of the reinforcement sample](images/reinforce_fem_details.png){width=1000}

The main plot is the mobilized axial force $T$ against position along the line, drawn over the dashed capacity
envelope of the [pullout section above](#determining-reinforcement-line-pullout-lengths): the friction ramp
developing from each free end over its pullout length $L_p$, the tensile plateau at $T_{max}$ in the middle, and
the step to the end anchorage capacity $T_{end}$ where one is declared. That envelope is the same expression the
solver evaluates at each element centroid to set $T_{allow}$, so the curve and the element capacities cannot
disagree. Where $T_{res}$ is filled in, the residual capacity is drawn as a dotted step beneath it — flat at
$T_{res}$ along the middle of the line, and following the friction ramp wherever the embedment develops less than
that — and elements that have softened onto it are marked. An element left with no residual at all, and therefore
carrying no force, is marked at zero.

The greatest utilization along a line is usually held over a stretch rather than at a point, the force being capped
by a flat envelope. Where it is a point, that point is ringed. Where it is a stretch, every sample on the stretch is
ringed and the run of curve between them is thickened — and a stretch with a sample inside it that stands below the
rest is drawn as the runs it really is, so a break in the thickened curve is where the line comes off capacity. The
legend calls the mark **At capacity**, or **Peak utilization** on a line that never reaches its envelope, and the
title states the fraction of capacity the peak reaches.

A band is shaded behind the profile where the viscoplastic shear strain concentrates along the line, called
**Shear band crossing** in the legend: what a profile along one line can show is where the band crosses it. The
extent is measured by walking the line itself and sampling the soil's shear strain field at every step, so the
shading is where the band crosses and not which bar elements hold the crossing. A line the concentration does not
reach carries no mark. Which field the band was read from, the mechanism an SSRM run captured or the shear strain
in a section that is standing, is what the title says.

Every mark on the panel is named in the legend, and nothing is labeled over the curves: the panel is wide and
shallow, and a label placed in it stands over the profile it describes.

Beneath the force profile is the bond transfer rate $dT/ds$: the force the ground hands the bar per unit of its
length, which is the gradient of the profile above it. There is no companion slip series because the formulation
has no slip degree of freedom — a reinforcement element is a truss bar on the continuum's own nodes, so bar and
soil displacement are the same number at every node. Load transfer is expressed through the capacity envelope, not through a slip law.

A **Field state** control at the foot of the panel selects which field the profiles are read from — the at-failure
mechanism an SSRM run captured, or the last converged solution — and is the same switch, with the same default, as
the one on the results view, so the two views can be set to the same instant of the analysis. It is dimmed for a
run that captured no mechanism, where there is only one field to read, and neither the capacity envelope nor the
shaded band moves with it. On a softening line the two fields can differ in kind: an element drops to its residual
only when an equilibrium state demands more than its capacity, so the last converged field may show no softened
element at all while the at-failure field — which starts from the set the failed-edge trial shed to — shows the
elements that gave way sitting on the residual line, marked *Softened*.

**Export** writes the current view as a PNG and its plotted series as a CSV named from the model, the line and the
field state, with that state also recorded in the CSV's header, so the picture and the numbers behind it stay
together. The panel is non-modal and reads the solution it was opened
with, so it can stay open beside the results view; it works the same on a solution reloaded from its saved
sidecar files as on a fresh solve.

The screenshot above is the reinforcement sample's own strength reduction run, $FS = 1.51$, read at the mechanism
it developed (see [FEM sample problems](samples.md)).

### The state of a line

One line is in one state, named the same way everywhere XSLOPE reports it: on the panel's list rows and under its
plot, in the title of a detail figure, in the table `print_reinforcement_summary()` prints, and in a generated
report.

| State | The line |
|-------|----------|
| within capacity | is below the capacity available to it everywhere along its length |
| near capacity | is below capacity everywhere, but close to it where it is most utilized |
| pullout | is slipping near an end at the capacity its embedment can develop there |
| yielded | is at its full tensile capacity away from the ends and holding it |
| softened | has dropped off its peak capacity onto its residual |
| ruptured | has softened with no residual capacity left and now carries nothing |
| inactive | carries no tension anywhere and is not engaged |

The two middle states are the ones worth separating, and they read alike on a badge: both are a line standing at
100% of what is available to it. **Pullout** is an end element at the reduced capacity its embedment can develop —
the friction ramp doing what a friction ramp does, with the interior of the line still below capacity. **Yielded**
is an element out on the $T_{max}$ plateau, where the whole tensile strength of the geosynthetic is mobilized. A
line in both states at once is reported yielded, the more serious of the two. **Softened** and **ruptured** need a
$T_{res}$: a line that declares none cannot reach them.


## References

Duncan, J.M., & Wright, S.G. (2005). *Soil Strength and Slope Stability*. John Wiley & Sons.

Griffiths, D.V., & Lane, P.A. (1999). Slope stability analysis by finite elements. *Geotechnique*, 49(3), 387-403.

Smith, I.M., & Griffiths, D.V. (2004). *Programming the Finite Element Method* (4th ed.). John Wiley & Sons.
