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

Each truss element stands on every node of the 2D element edge it lies on. On a linear mesh that is the edge's two
end nodes and the element is a 2-node bar. On a quadratic mesh (tri6, quad8, quad9) the edge also carries a midside
node, and the truss element is a 3-node bar carrying that node too.

The midside node is what ties the bar to the soil in the middle of the edge. The soil's displacement along a quadratic
edge is a parabola through all three nodes, so a bar attached at the corners alone leaves the edge free to bow away
from it between them, and leaves the midside node free to slide along it: the bar and the soil around it displace
together at only half the stations the edge has. Carrying the node closes that gap, and it costs no node and moves
none, because the node is already there as part of the 2D element.

The 3-node bar is the standard isoparametric quadratic bar, and its axial force at the element center is
$EA(u_2 - u_1)/L$ — the same chord expression the 2-node bar uses, so an element's reported force means what it always
did.

The meshing algorithms used in XSLOPE, including the integration of 1D and 2D elements for problems involving soil
reinforcement are documented in the [Mesh Generation](mesh.md) page.

## Mathematical Formulation

**Truss Element Stiffness Matrix:** Each 1D truss element contributes to the global stiffness matrix through its element stiffness matrix. On a linear mesh the element has two nodes $i$ and $j$, and its stiffness in local (axial) coordinates is:

>>$[K_e]_{local} = \dfrac{AE}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$

On a quadratic mesh the element also carries the midside node $m$ of its soil edge, and its stiffness is the quadratic
bar's, in the node order $(i, j, m)$:

>>$[K_e]_{local} = \dfrac{AE}{3L} \begin{bmatrix} 7 & 1 & -8 \\ 1 & 7 & -8 \\ -8 & -8 & 16 \end{bmatrix}$

where $A$ is the cross-sectional area, $E$ is the elastic modulus, and $L$ is the element length.

**Coordinate Transformation:** The local stiffness matrix must be transformed to global coordinates using the transformation matrix $[R]$:

>>$[K_e]_{global} = [R]^T [K_e]_{local} [R]$

The transformation is built from $\psi$, the inclination of the reinforcement line to the horizontal — the same
angle the LEM formulation uses for the direction of an axial reinforcement force:

>>$[R] = \begin{bmatrix} \cos\psi & \sin\psi & 0 & 0 \\ 0 & 0 & \cos\psi & \sin\psi \end{bmatrix}$

with one more row, $\begin{bmatrix} 0 & 0 & 0 & 0 & \cos\psi & \sin\psi \end{bmatrix}$, for the midside node of a
three-node bar.

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

The profile draws no mark for where the shear band crosses the line. The field figure is where a crossing is
read — the band is drawn there as strain contours — and the profile is where what the line carries is read; a
shaded stretch derived from a threshold on the sampled strain added a rule between the two that the legend
could not explain.

Every mark on the panel is named in the legend, and nothing is labeled over the curves: the panel is wide and
shallow, and a label placed in it stands over the profile it describes.

Beneath the force profile is the bond transfer rate $dT/ds$: the force the ground hands the bar per unit of its
length, which is the gradient of the profile above it. There is no companion slip series because the formulation
has no slip degree of freedom — a reinforcement element is a truss bar on the continuum's own nodes, so bar and
soil displacement are the same number at every node. Load transfer is expressed through the capacity envelope, not through a slip law.

A **Field state** control at the foot of the panel selects which field the profiles are read from — the at-failure
mechanism an SSRM run captured, or the last converged solution — and is the same switch, with the same default, as
the one on the results view, so the two views can be set to the same instant of the analysis. It is dimmed for a
run that captured no mechanism, where there is only one field to read, and the capacity envelope does not move with it. On a softening line the two fields can differ in kind: an element drops to its residual
only when an equilibrium state demands more than its capacity, so the last converged field may show no softened
element at all while the at-failure field — which starts from the set the failed-edge trial shed to — shows the
elements that gave way sitting on the residual line, marked *Softened*.

**Export** writes the current view as a PNG and its plotted series as a CSV named from the model, the line and the
field state, with that state also recorded in the CSV's header, so the picture and the numbers behind it stay
together. The panel is non-modal and reads the solution it was opened
with, so it can stay open beside the results view; it works the same on a solution reloaded from its saved
sidecar files as on a fresh solve.

The screenshot above is a strength reduction run on the reinforced slope built in
[FEM-2](../tutorials/fem02_reinforcement.md), read at the mechanism it developed.

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


## Two Ways to Represent a Sheet

Everything above is the **bonded bar**: the truss element shares the soil's own nodes, so the soil above the line
and the soil below it are one body and the sheet cannot slide in the soil except by the bar reaching the capacity
its embedment develops. "Bond" there is a cap on the tension, read from the pullout envelope at each element's
position, and not a sliding surface.

The alternative is the **joint**. Setting `Joint = Yes` on a reinforcement line makes that line a slip surface: the
mesh is split along it, the sheet becomes a bar with its own nodes between an upper and a lower interface, and each
interface carries the line's own `Adhesion` and `Delta` as a Mohr-Coulomb strength. The soil on the two sides can
then slide on the sheet, and on each other, at that interface strength.

```
  bonded (the default)                joint (Joint = Yes)

   soil above                         soil above
   ----o----o----o----   the bar      ----o----o----o----  upper face nodes
       |    |    |       shares       ~~~~~~~~~~~~~~~~~~~  upper interface
   ----o----o----o----   these        ====b====b====b====  bar nodes (their own)
   soil below            nodes        ~~~~~~~~~~~~~~~~~~~  lower interface
                                      ----o'---o'---o'---  lower face nodes
                                      soil below           (all three at one point)
```

Every node on a jointed line exists three times at the same point — a copy for the soil above, a copy for the bar,
a copy for the soil below — and two interface elements connect them at each station: the soil above against the
bar, and the bar against the soil below. On a quadratic mesh the midside nodes are tripled too, so an interface
element has three node pairs and matches the adjacent triangles' edges. The split is invisible on a mesh plot,
because the three copies stand at one point; the jointed line is drawn in its input style with short ticks on both
sides so it can be told from a bonded one.

### The interface element

Each interface element is a zero-thickness joint (Goodman, Taylor & Brekke, 1968). Its state is the relative
displacement of the two faces, resolved into a normal component $\Delta_n$ (closing, compression positive) and a
tangential component $\Delta_t$ (sliding), and its constitutive law is a traction-displacement relation:

$$
t_n = k_n \Delta_n, \qquad t_s = k_s \Delta_t, \qquad |t_s| \le c_j + t_n \tan\phi_j
$$

with $c_j$ and $\phi_j$ the line's `Adhesion` and `Delta`. Slip past the limit is perfectly plastic: the shear
traction stays at the limit while the tangential offset grows. A normal traction below the tension cutoff — zero on
a reinforcement line, because a soil-geosynthetic contact carries no tension — **opens** the joint, which then
carries neither traction until the faces come back into contact.

The tractions are integrated at the element's own nodes rather than at Gauss points. Gauss quadrature on a
zero-thickness interface produces traction oscillations that grow with the penalty stiffness; nodal integration
does not, and the traction spread along an element measured on the direct-shear case *falls* from 15% to 1.8% when
$k_n$ is multiplied by a hundred.

$k_n$ and $k_s$ are penalty stiffnesses: large enough that an intact joint does not visibly deform, small enough
not to ill-condition the system. Left blank they are derived as $E_{adj}/d_v$ and $G_{adj}/d_v$ over a virtual
thickness $d_v = 0.1\,L_{1D}$, with $E_{adj}$ and $G_{adj}$ those of the softer of the two soils the element stands
between. The factor of safety is insensitive to them — it moves by about 1.5% over two orders of magnitude — but
the **cost** is not: the slip a viscoplastic sweep puts into a joint is the excess traction divided by $k_s$, so a
model that states a $k_n$ / $k_s$ an order of magnitude above the derived default needs its iteration limit raised
by the same factor, or the strength reduction reports the iteration budget instead of the slope.

### Ends, ties and the bar

The two soil faces rejoin at each end of the line: one shared soil node, a crack tip. The bar's end node is not
that node — it is the bar's own, connected to the soil only through the interfaces along it — so **a sheet end is
free by default**. It carries no load, it can pull out, and the tension the bar develops at any station is the
interface shear integrated from the free end, which is the pullout envelope produced by the elements instead of by
hand.

An end is **tied** when that end's `Tend1` / `Tend2` is filled in: a spring from the bar's end node to the soil (or
facing) node at the same point, perfectly plastic at the stated capacity, restraining both components with the
capacity read on the resultant. That is what those columns have always meant, and it is how a wall sheet is
connected to its facing.

On a jointed line the bar keeps its own law — tension only, rupture at $T_{max}$, softening to $T_{res}$ where
stated — **minus the bond-slip cap**. The grip on the soil is now the interface traction the joints integrate, so
applying the pullout envelope as well would count it twice; `Lp1` and `Lp2` are not read there. Preflight says so.

### Strength reduction

A strength reduction divides $c_j$ and $\tan\phi_j$ by the trial factor along with the soil's, on every jointed
line unless that line sets `Jred = No`. The stiffnesses $k_n$ and $k_s$, the ties and $T_{max}$ are structural and
are not reduced, exactly as the bar's properties are not.

A jointed model takes the viscoplastic solver's verdict directly: the Newton corrector is not offered one. On a
slipping interface the shear traction is held at $c_j + t_n \tan\phi_j$ and does not depend on the tangential
displacement, so a state the viscoplastic loop reached by growing the slip leaves a corrector — which may only move
displacements — with nothing to move.

What ends a jointed trial is therefore the viscoplastic loop alone, and the force test it has to pass is measured
almost entirely on the joints. A trial that neither converges nor fails is read from the **slip** instead — steady
slip with the field gaining is a slope moving on its joints and is `FAILED`; a slip and a field that have both
stopped, with the soil in equilibrium, is a slope standing behind a limit cycle no budget brings down. Both readings,
and why the ordinary displacement evidence cannot make the distinction, are in
[the joint verdict](overview.md#the-joint-verdict).

### Making a jointed trial cheaper

The slip a sweep puts into a joint is the excess traction over $k_s$, so at $\Delta t = 1$ one sweep returns the
shear traction exactly to its limit at the current displacement field. `joint_slip_stiffness_factor` scales the
stiffness that division uses **on the pairs that are at their limit**, and nothing else: the traction limit, the
assembled elastic block, and a pair that re-sticks are all untouched. A factor below 1 drives the interface past its
own return and is therefore an over-relaxation, not a change of physics.

It is **off by default**, and the reason is a stability limit rather than a preference. Because one sweep already
returns the traction exactly, the iteration is *at* its boundary at a factor of 1 and a factor $f$ is a relaxation of
$1/f$. RS2 sets its equivalent to 0.01 — a relaxation of 100 — and at that value this scheme diverges outright.
Values near 1 are stable and buy nothing. Reach for the iteration budget and the [joint
verdict](overview.md#the-joint-verdict) instead; the factor exists so the setting can be reproduced, not because it
is a knob worth turning.

## Joints Without Reinforcement

Not every slip surface has a sheet in it. A rock joint, a bedding plane, the contact between a concrete facing block
and the one beneath it, the back of a retaining wall against the soil it holds: each is a surface two bodies meet on
and can slide along, with no member between them. Those go on the **joints** worksheet, one row per line:

| column | what it is |
|---|---|
| `Label` | the name the plots, the details view and the report use |
| `x1`, `y1`, `x2`, `y2` | the line's endpoints |
| `c` | the interface's cohesion (blank = 0) |
| `phi` | its friction angle, in degrees — required |
| `c_res`, `phi_res` | what it drops to once it has slipped (blank = the peak, i.e. no residual branch) |
| `dil` | its dilation angle, in degrees: a unit of slip opens it by $\tan$`dil` (blank = 0) |
| `t_cut` | the tension cutoff at which the two faces part (blank = 0) |
| `kn`, `ks` | the penalty stiffnesses (blank = derived, as above) |
| `Jred` | blank or **Yes** reduces the joint with the soil in a strength reduction; **No** holds it |

The mesh split, the interface element, the derived stiffnesses and the strength reduction are all exactly what a
jointed reinforcement line gets, described above. Two things this sheet states that the reinforce sheet does not
are described under [Residual strength and dilation](#residual-strength-and-dilation) below. The other difference
is that there is no bar between the two faces, so a station carries two coincident nodes instead of three and
**one** interface element spans them instead of two:

```
  Joint = Yes on a sheet             a joints-sheet line

   ----o----o----o----  upper        ----o----o----o----  upper face nodes
   ~~~~~~~~~~~~~~~~~~~  interface    ~~~~~~~~~~~~~~~~~~~  the interface
   ====b====b====b====  the bar      ----o'---o'---o'---  lower face nodes
   ~~~~~~~~~~~~~~~~~~~  interface
   ----o'---o'---o'---  lower
```

A sheet's two interfaces act in series, so the pair carries the stated stiffness twice over; a joints-sheet line
carries it once, which is what a single contact is. Neither difference reaches the strength: `c` and `phi` are the
Mohr-Coulomb limit the faces slide at either way. The tension cutoff is a column of its own here, where a
reinforcement line's is fixed at zero, because a rock joint may hold a little tension across it and a
soil-geosynthetic contact does not.

Joint lines may **meet** — at a T, at a crossing, at a corner, end to end. Where they do, the mesh split counts the
wedges of material around the shared node and gives the node one copy per wedge, so each element keeps the material
on its own side of every line through the point. That is what makes a block column with a joint on its back face, a
joint under its base and a course joint at every mortar line into a stack that can slide, part and rock, rather than
a notched solid. A joint that ENDS on another — a column's base beginning partway along its neighbour's side joint,
a release trace running down onto a bedding plane — is that same rule with three wedges instead of four, and it needs
the two lines to meet at one point: an end that lies within a millionth of the section of another joint line is
moved onto it, and that line is given a vertex there, before the mesh is built. So a tip stated to fewer decimals
than the line it belongs on still lands on it. The through line keeps the geometry it was stated with.

What joint lines may not do is lie **on** one another over a stretch, or run along the outside of the section, where
there is material on one side only and nothing for the other face to be; preflight refuses both by name.

A joint line is finite element geometry: the limit equilibrium engines do not read the joints worksheet at all.

### Residual Strength and Dilation

A rock joint is not a smooth plane, and two of its columns say so.

**`c_res` and `phi_res`** are the strength the joint keeps once it has slipped. A rough surface shears through its
asperities the first time it reaches its limit and does not rebuild them, so the drop is instantaneous and
permanent: the pair's limit falls from

$$S_{peak} = c + t_n \tan\phi \qquad\text{to}\qquad S_{res} = c_{res} + t_n \tan\phi_{res}$$

on the sweep after the one that first found it at its limit, and it stays on the residual branch for the rest of the
run even where it later closes or unloads. A blank leaves the joint on its peak strength throughout, which is what a
joint with no residual branch is. Neither residual may exceed its peak; preflight refuses one that does. A strength
reduction divides both branches by the trial factor, on the same `Jred` switch, so a joint that opts out of the
reduction opts both branches out.

**`dil`** is the dilation angle. A rough joint rides up on its asperities as it slides, so a slip increment
$|\Delta u_s|$ opens it by $|\Delta u_s| \tan(\text{dil})$. That opening is a plastic normal offset: the elastic
part of the normal displacement, and with it the normal traction

$$t_n = k_n (\Delta_n + u_{open})$$

grows while the joint slides. Where the material around the joint holds it closed, that is dilatant hardening —
sliding builds normal stress and with it shear strength. Where the sliding block is free to lift, it lifts instead,
and the normal traction stays at whatever equilibrium with the block's weight requires. The dilation is
**non-directional**: the joint opens whichever way it slides, and it does not die out with accumulated
slip — a joint that slides a long way keeps riding up at the stated angle. A blank is zero and the normal
traction is $k_n \Delta_n$ exactly.

Neither column exists on the reinforce sheet. A soil-geosynthetic contact is a frictional interface with one
strength, and the sheet between its two faces is what carries the mechanism.

### Both kinds in one model

A model may carry both, and the multi-tiered geotextile wall is the case that needs both at once: each sheet is a
reinforcement line with `Joint = Yes`, so the fill can slide on it while the sheet still carries tension, and each
facing column stands on joints-sheet lines that let the blocks slide on each other, on the fill behind them and on
the foundation beneath. A sheet whose front end stops on the back face of a column is **tied** there, and the tie
takes the material the line runs into past its end — the column — so the wrap's connection to the facing is what the
tie represents.

## Bonded Bar or Joint?

The choice is about the mechanism, not about the material: does the slip surface **cut** the reinforcement, or run
**along** it?

**A bonded bar is right where the surface crosses the layers.** A circle through a geogrid slope, a nail wall, a
pile row: the soil on both sides of each layer moves together, the layer carries tension across the surface, and
the only interface question is pullout, which is a capacity the bar's cap answers. It is the finite element twin of
the limit equilibrium treatment — a force where the surface crosses the line — so the two engines compare like for
like.

**A joint is right where the surface can run along the layer.** A reinforced embankment on soft clay sliding on
its base geotextile, a wrapped-face or block-faced wall where the fill between the sheets moves relative to them
and each sheet anchors to a facing, a smooth geomembrane or liner whose interface friction is well below the soil's,
any long flat sheet under a sliding mass. The interface shear strength along the sheet governs, the two sides move
differently, and a bonded bar cannot represent it: it reports the bars at their cap while the mesh decides the
answer.

Joints are off by default, because tripled nodes and a penalty stiffness cost something. On a wall or a base sheet,
running both ways settles it: a bonded answer that matches the jointed one says the mechanism is crossing, and one
that sits above it says the mechanism was sliding and the bond was holding it up.

### What says a joint was needed

Preflight cannot know the mechanism, but four of these cases show in the inputs, and it reports each of them as a
warning naming the line:

- a sheet **lying on a material boundary** over most of its length — the base geotextile of a reinforced
  embankment, which slides on its own interface;
- a **long, flat sheet** within a few degrees of horizontal spanning most of the width of the zone above it;
- a **smooth interface**, `Delta` below about 0.6 of the surrounding soil's $\phi$ — a geomembrane or liner rather
  than an ordinary soil-geosynthetic contact;
- the **wall pattern**: four or more near-horizontal sheets at a vertical spacing under a meter, behind a face at
  least 70° steep, with their front ends inside a thin facing column.

Two stronger signals come out of a bonded run itself, and a strength reduction reports them when it finds them: the
shear strain band at the critical factor running **along** a sheet rather than across it, and **every bar element
on one sheet at its capacity** — a sheet held up by a cap the bond set rather than by a grip the mesh resolved,
which refines as the bar elements refine instead of converging.

### What the results show

The results view draws each interface as a thin line on the line it runs along, colored by how far its two faces have
slid, on a colorbar titled *Joint slip*. The ramp is green — bright lime at the smallest slip through to dark green at
the largest — because the strain field under it runs blue through white to red, and a slip color the field can also
produce is a slip color you cannot find. A slipping span is backed by a thin white stroke, so it still reads where the
field goes dark blue or red. A joint that is not slipping is a lighter neutral gray hairline with no backing, and a
stretch that has **opened** is marked with a short tick across the line rather than given a color of its own, because
opening is a condition and not a quantity. A model where no joint slipped carries no colorbar. The weight is
deliberate: a generated network puts hundreds of traces over the field, and anything heavier covers what it is drawn
on. A joint carries no strain of its own, so it appears in the shear-strain field only through what it does to the
soil beside it; this is its own reading.

On a jointed model the displacement panel is the scaled deformed mesh rather than an arrow field, drawn as the
**blocks** the joints cut the section into: each block under a faint tint of its own, its joint faces in the same
green, the outside of the deformed mesh as a dark line against the dashed undeformed outline, and the element edges
in light gray behind them — until there are more than eight jointed lines, when the edges come off entirely and the
block outlines carry the panel. A block is a piece of the mesh that moves as one body, found by following element
adjacency: the split gives the two sides of a joint their own nodes, so they are no longer neighbors, while a joint
that stops inside the mass leaves the material wrapped around its tip as one block.

The exaggeration is bounded: no point of the deformed mesh is drawn more than 4% of the section away from where it
started, whatever multiplier that takes, and where the largest displacement is below a hundred-millionth of the
section there is nothing to draw — the panel shows the undeformed mesh and its title says the deformation is below
drawing resolution, rather than magnifying the solver's own residue into a shape. That is what a
jointed failure looks like — blocks moving as bodies, with all of the movement taken up at the joints, where a
slipped or opened joint shows as two lines that no longer lie on each other. An arrow field samples that at nodes
and misses exactly the thing that happened.

**1D Details…** lists every jointed line under a *Joints* heading beside the reinforcement lines, and draws four
panels along the line: the bar's tension over its capacity, the normal traction the interface carries, the shear
traction with the Mohr-Coulomb limit $c_j + t_n \tan\phi_j$ drawn beside it, and the slip. Where the two shear
curves meet is where the interface is at its limit. A generated report carries the same reading as a table: the
share of each line's length standing at its limit, and the largest offset the two faces reached.

## References

Duncan, J.M., & Wright, S.G. (2005). *Soil Strength and Slope Stability*. John Wiley & Sons.

Griffiths, D.V., & Lane, P.A. (1999). Slope stability analysis by finite elements. *Geotechnique*, 49(3), 387-403.

Goodman, R.E., Taylor, R.L., & Brekke, T.L. (1968). A model for the mechanics of jointed rock. *Journal of the Soil Mechanics and Foundations Division*, 94(SM3), 637-659.

Schellekens, J.C.J., & de Borst, R. (1993). On the numerical integration of interface elements. *International Journal for Numerical Methods in Engineering*, 36(1), 43-66.

Smith, I.M., & Griffiths, D.V. (2004). *Programming the Finite Element Method* (4th ed.). John Wiley & Sons.
