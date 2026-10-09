# Piles and Concrete Piers in LEM Slope Stability

A pile resists a slide by shear and bending where the slip surface crosses it. The kinds of pile and wall, and how a
pile line is entered, are on the [Piles and Walls Overview](overview.md), and what to enter for each, with typical
values, on [Pile and Wall Types](types.md). In a limit equilibrium analysis the pile is one force on the sliding
mass, at the point where a trial slip surface crosses it:

![A row of piles through a sliding mass, each pushing back on it with a force H where the failure surface crosses the pile](../lem/images/pile_diagram.png){width=598}

Each pile pushes back on the sliding mass with a force $H$ at the point where the failure surface crosses it,
against the direction the mass moves.

## Pile Force in Limit Equilibrium Analysis

The pile force enters each slice's equilibrium through its magnitude, its direction and the point where it acts.

### Force Definition

In the limit equilibrium framework, a pile is characterized by:

- A **force magnitude** $H$ (per unit width of slope, i.e., force/length) acting at the point where the pile intersects the failure surface
- A **force angle** $\theta_p$, the inclination of that force from horizontal (positive upward), which XSLOPE computes from the pile's end points so that the force is perpendicular to the pile (see [Force Direction](#force-direction))

The force is decomposed into horizontal and vertical components:

>$H_h = H \cos\theta_p \qquad \text{(horizontal)}$

>$H_v = H \sin\theta_p \qquad \text{(vertical, positive upward)}$

For a vertical pile, the most common case, $\theta_p = 0$ and the force is purely horizontal: $H_h = H$ and $H_v = 0$.

### Force Resolution on the Slice

The pile force acts at point $e$ on the failure surface where the pile intersects the slice base. Relative to the slice base (inclined at angle $\alpha$), the force resolves into:

>**Normal to base**: $H \sin(\alpha - \theta_p)$

>**Tangential to base** (resists sliding): $H \cos(\alpha - \theta_p)$

The tangential component directly opposes the sliding force. The normal component adds to the normal force on the base where $\alpha > \theta_p$, raising the frictional resistance, and subtracts from it where $\alpha < \theta_p$: near the toe of a circle, where the base rises toward the toe, or under a force tilted upward more steeply than the base. The Ordinary Method of Slices takes $N'$ from this resolution; the other methods take it from their own equilibrium equations, listed under [Integration with LEM Methods](#integration-with-lem-methods).

### Moment Contribution

For methods that use moment equilibrium about a circle center $(X_o, Y_o)$, the pile force creates a resisting moment. The horizontal and vertical components of $H$ each contribute through their respective moment arms:

>$M_{\text{pile}} = H \cos\theta_p \cdot (Y_o - y_e) + H \sin\theta_p \cdot (x_e - X_o)$

where $(x_e, y_e)$ is the pile-failure surface intersection point. The equation is written for a slope that descends to the left; on one that descends to the right the arm $x_e - X_o$ changes sign. When $\theta_p = 0$, this reduces to $M_{\text{pile}} = H(Y_o - y_e)$.

How the pile force enters the factor of safety is set by the row's **Appl** entry:

- **Active** (the default, and how a blank cell is read): $H$ is an allowable force, applied at full value and **not** divided by the factor of safety $F$. It is a known force, like the weight of the slice: in the moment equations its moment reduces the driving moment (the denominator), and its components enter the force equations, including the one for $N'$, at full value.
- **Passive**: $H$ is an ultimate capacity that is mobilized with the soil strength, so its components are divided by $F$, as $c$ and $\tan\phi$ are. In the moment equations the passive pile moment joins the resisting side (the numerator). The Ordinary Method of Slices, which computes $F$ directly without iteration, takes that moment only and leaves the passive force out of $N'$.

### Force Direction

A pile is entered by its two end points, $(x_1, y_1)$ and $(x_2, y_2)$, in either order: XSLOPE takes the higher end as the head and the lower end as the tip. The pile force is directed perpendicular to the pile's axis, with its horizontal component against the movement of the sliding soil, so the pile resists the slide whether the slope descends to the left or to the right. The force angle is therefore the pile's inclination from vertical:

>$\theta_p = \tan^{-1}\!\left(\dfrac{d_u}{y_h - y_t}\right)$

where $y_h$ and $y_t$ are the elevations of the head and the tip, and $d_u$ is the horizontal distance by which the tip lies upslope of the head (negative when the tip lies downslope, toward the toe). Upslope is opposite to the movement of the sliding mass, which XSLOPE determines for each failure surface, so a model and its mirror image give a pile the same angle.

- **Vertical pile** ($x_2 = x_1$): $\theta_p = 0$, a horizontal force pushing back against the moving soil.
- **Battered pile, tip upslope of the head** (the head leans out toward the toe): $\theta_p > 0$, and the force tilts **upward**.
- **Battered pile, tip downslope of the head** (the head leans back into the slope): $\theta_p < 0$, and the force tilts **downward**.

The Ito & Matsui computation described below applies to vertical piles only, so a battered pile needs its $H$ entered. If a workbook's `piles` sheet carries a $\theta_p$ column, an angle entered there is used in place of the computed one, measured the same way, and a blank cell leaves the computed angle.

### Per-Unit-Width Convention

All forces in 2D limit equilibrium analysis are expressed per unit width of slope (perpendicular to the cross-section). If a row of piles has individual capacity $H_{\text{single}}$ at center-to-center spacing $S$, the equivalent force per unit width is:

>$H = \dfrac{H_{\text{single}}}{S}$

### Non-Circular Failure Surfaces

The pile force formulation for Janbu, Corps of Engineers, Lowe-Karafiath, and Spencer's method uses only per-slice quantities ($\alpha$, force components) and has no dependence on circle geometry. These methods work with any failure surface shape, and the pile terms carry over without modification. The only circle-dependent pile terms appear in OMS and Bishop (the moment term), but those methods inherently require circular surfaces.

### Integration with LEM Methods

The pile force $H$ at angle $\theta_p$ is incorporated into each limit equilibrium method supported by XSLOPE. The specific modifications for each method are presented in the respective method documentation pages:

- [**OMS**](../lem/oms.md): $H\sin(\alpha-\theta_p)$ added to $N'$; pile moment terms added to the denominator
- [**Bishop**](../lem/bishop.md): $-H\sin\theta_p$ enters vertical equilibrium for $N'$; pile moment terms added to the denominator
- [**Janbu**](../lem/janbu.md): $-H\sin\theta_p$ enters vertical equilibrium for $N'$; $-H\cos\theta_p$ enters the horizontal force balance
- [**Force Equilibrium** (Corps of Engineers, Lowe-Karafiath)](../lem/force_eq.md): $-H\cos\theta_p$ added to horizontal equilibrium ($b_0$); $-H\sin\theta_p$ added to vertical equilibrium ($b_1$)
- [**Spencer**](../lem/spencer.md): $H\cos\theta_p$ added to $F_h$; $H\sin\theta_p$ added to $F_v$; moment terms added to $M_o$

In all methods, for a vertical pile ($\theta_p = 0$) the equations reduce to the simpler horizontal-force-only case.


## Determining the Pile Force $H$

### User-Specified Force

The simplest approach is for the user to specify $H$ directly based on external analysis. The pile force may come from:

- p-y curve analysis software (e.g., LPILE, GROUP, RSPile)
- Structural analysis of the pile section
- Published design charts or empirical correlations
- Full 3D finite element analysis

When using user-specified forces, the user enters $H$ (per unit width) in the `piles` sheet of the input template or in Studio's Piles editor, and XSLOPE computes $\theta_p$ from the pile's end points. This approach gives the user full control and is appropriate when detailed pile analysis has already been performed.

### Ito & Matsui (1975) Theory

The Ito & Matsui method is the most widely used closed-form approach for computing the lateral force that soil exerts on passive stabilizing piles. It models the soil between adjacent piles as being in a state of **plastic equilibrium** — the soil deforms plastically as it squeezes between the piles, like material flowing through a constriction. Using Mohr-Coulomb plasticity theory, Ito & Matsui derived closed-form equations for the lateral pressure on the piles as a function of depth.

#### Setup and Notation

Consider a row of piles embedded in a slope:

![Ito & Matsui plan view](../lem/images/pile_ito_matsui_plan.png)

- $D$ = pile diameter (or width of the pile cross-section)
- $S$ = center-to-center spacing between piles
- $D_1 = S - D$ = clear spacing between adjacent pile faces
- $z$ = depth below the ground surface
- $z_f$ = depth from the ground surface to the failure surface at the pile location
- Soil properties: cohesion $c$, friction angle $\phi$, unit weight $\gamma$
- Passive earth pressure coefficient: $N_\phi = \tan^2\!\left(45° + \dfrac{\phi}{2}\right)$

The theory applies to the portion of the pile **above** the failure surface — this is the zone where soil is actively moving and pushing against the pile.

**Notation**: The original paper uses $d$ for pile diameter, $D_1$ for center-to-center spacing, and $D_2$ for clear spacing. XSLOPE uses $D$, $S$, and $D_1$ respectively. The equations are identical; only the symbols differ.

#### General $c$-$\phi$ Soil

For a soil with both cohesion and friction ($c > 0$, $\phi > 0$), the distributed lateral force $p(z)$ (force per unit length of pile) at depth $z$ is:

>$p(z) = c \cdot A_1 + \gamma z \cdot A_2$

where $A_1$ and $A_2$ are arching coefficients with units of length. Let $S = D_1 + D$ (center-to-center spacing). The coefficients are computed from three intermediate quantities:

>$R = \left(\dfrac{S}{D_1}\right)^{\sqrt{N_\phi}\,\tan\phi + N_\phi - 1}$

>$\mathcal{E} = \exp\!\left(\dfrac{D}{D_1}\,N_\phi\tan\phi\,\tan\!\left(\dfrac{\pi}{8}+\dfrac{\phi}{4}\right)\right)$

>$F = \dfrac{2\tan\phi + 2\sqrt{N_\phi} + N_\phi^{-1/2}}{\sqrt{N_\phi}\,\tan\phi + N_\phi - 1}$

$R$ is the geometric amplification from the ratio of center-to-center to clear spacing, $\mathcal{E}$ captures the exponential plastic flow between piles, and $F$ is a dimensionless grouping of friction and earth pressure terms.

The overburden coefficient (from Eq. 14, the $c = 0$ specialization of Eq. 13) is:

>$A_2 = \dfrac{S \cdot R \cdot \mathcal{E}\, -\, D_1}{N_\phi}$

The cohesion coefficient (from Eq. 13) is:

>$A_1 = \dfrac{S \cdot R\,(\mathcal{E} - 2\sqrt{N_\phi}\,\tan\phi - 1)}{N_\phi\tan\phi} + S \cdot F\,(R - 1) + \dfrac{2\,D_1}{\sqrt{N_\phi}}$

These expressions were verified against the original paper's parametric charts (Figs. 7–9) and field measurements (Table 1, Figs. 13–14).

Key behavior of $p(z)$:

- **Increases with depth** through the $\gamma z$ term — deeper soil mobilizes more pressure against the pile
- **Increases as $D_1/D$ decreases** (closer piles = more arching = more force per pile)
- **Increases with $\phi$** — higher friction angle produces stronger soil arching between piles
- **Increases with $c$** — cohesion contributes a constant (depth-independent) component

#### Cohesionless Soil ($c = 0$)

For a purely frictional soil with $c = 0$, the cohesion term vanishes and the lateral pressure is:

>$p(z) = \gamma z \cdot A_2$

where $A_2$ is the same expression as above. The pressure increases linearly from zero at the ground surface.

#### Undrained Clay ($\phi = 0$)

For a purely cohesive (undrained) soil with $\phi = 0$, the general $c$-$\phi$ expressions become indeterminate because $N_\phi = 1$ and $\tan\phi = 0$. Deriving the solution independently for $\phi = 0$ (Ito & Matsui Eq. 23) yields:

>$p(z) = c_u \left[S\left(3\ln\dfrac{S}{D_1} + \dfrac{D}{D_1}\tan\dfrac{\pi}{8} - 2\right) + 2D_1\right] + \gamma z \cdot D$

where $S = D_1 + D$ is the center-to-center spacing and $c_u$ is the undrained shear strength. The first term represents the cohesion contribution (constant with depth) and the second term represents the overburden contribution (linear with depth), which simplifies to $\gamma z \cdot D$ (the pile diameter).

#### Total Force per Pile

![Ito & Matsui pressure distribution](../lem/images/pile_ito_matsui_pressure.png)

The total lateral force on a single pile is obtained by integrating $p(z)$ from the ground surface down to the failure surface depth $z_f$:

>$F_{\text{pile}} = \int_0^{z_f} p(z) \, dz$

Since $p(z)$ is linear in $z$ within a homogeneous layer, the integration is straightforward:

>$F_{\text{pile}} = c \cdot A_1 \cdot z_f + \gamma \cdot A_2 \cdot \dfrac{z_f^2}{2}$

#### Force per Unit Width

The 2D plane-strain equivalent force used in LEM (per unit width of slope) is:

>$H = \dfrac{F_{\text{pile}}}{S}$

This is the value entered (or computed) for the pile force in the slope stability analysis.

#### Multi-Layer Soils

When the pile passes through multiple material zones above the failure surface (common in practice), the integration is performed piecewise. For each layer $j$ with properties $c_j$, $\phi_j$, $\gamma_j$ between depths $z_{\text{top},j}$ and $z_{\text{bot},j}$:

>$F_j = c_j \cdot A_{1,j} \cdot (z_{\text{bot},j} - z_{\text{top},j}) + \gamma_j \cdot A_{2,j} \cdot \dfrac{z_{\text{bot},j}^2 - z_{\text{top},j}^2}{2}$

The total force per pile is the sum over all layers:

>$F_{\text{pile}} = \sum_j F_j$

Note that $A_{1,j}$ and $A_{2,j}$ must be recomputed for each layer since $\phi$ may differ between layers. The pile geometry ($D$, $D_1$) remains the same for all layers.

#### Computation at Each Trial Surface

An important characteristic of the Ito & Matsui calculation is that **$H$ depends on the failure surface location**. A deeper failure surface means more soil above it pushing on the pile, giving a higher $H$. Therefore, $H$ should be recomputed for each trial failure surface during an automated search. Since the computation involves only closed-form expressions and simple integration, it is essentially instantaneous and adds no meaningful computational cost.

In XSLOPE, when $H$ is left blank in the `piles` sheet but the pile diameter $D$ and spacing $S$ are provided, the Ito & Matsui force is computed automatically at slice generation time for each trial surface. If the user provides an explicit $H$ value, that value is used instead (override mode).

#### Soil Arching Between Piles

A critical aspect of pile-stabilized slopes is the three-dimensional soil arching that develops between adjacent piles. As the sliding soil mass pushes against the pile row, stress concentrations develop around each pile, and the soil "arches" between piles in a manner analogous to arching above a tunnel. This is the mechanism captured by the Ito & Matsui theory.

The effectiveness of soil arching depends on:

- **Pile spacing**: Closer spacing produces stronger arching and higher force per pile. The optimal spacing balances structural efficiency (fewer piles) against arching effectiveness (closer piles).
- **Soil strength**: Stronger soils develop more effective arching. In very weak soils (soft clay), arching may be minimal and the soil may flow between the piles without mobilizing significant resistance.
- **Pile rigidity**: Rigid piles provide fixed points for arch development. Flexible piles may deflect enough to reduce arching effectiveness.

For design purposes, $S/D$ ratios of 3 to 6 are typical for slope stabilization applications.

#### Low Friction Angle Floor

The general $c$-$\phi$ equation (Eq. 13) and the undrained clay equation (Eq. 23) were derived independently using different mathematical approaches. As $\phi \to 0$, the $c$-$\phi$ equation does **not** converge to the $\phi = 0$ result — the overburden coefficient $A_2$ drops to near zero for small $\phi$ before recovering and exceeding the $\phi = 0$ value at approximately $\phi = 12$–$15°$. This creates an unphysical discontinuity where a soil with $\phi = 2°$ would produce less pile force than one with $\phi = 0°$.

Since friction can only strengthen soil arching (and thus increase the lateral force on the pile), XSLOPE enforces the $\phi = 0$ coefficients as a lower bound for all friction angles:

>$A_1 = \max(A_{1,c\text{-}\phi},\; A_{1,\phi=0}) \qquad A_2 = \max(A_{2,c\text{-}\phi},\; A_{2,\phi=0})$

This ensures that the computed pile force increases monotonically with $\phi$. For further discussion of limitations in the Ito & Matsui formulation, see [Ukritchon & Keawsawasvong (2017)](https://doi.org/10.1061/(ASCE)GT.1943-5606.0001753).

#### Applicability and Limitations

The Ito & Matsui method has the following characteristics and limitations:

- **Rigid pile assumption**: The theory assumes piles do not deflect significantly. This is conservative for flexible piles, which mobilize less soil pressure than rigid piles.
- **Plastic flow assumption**: The method gives an **upper bound** on the soil's capacity to push on the pile. The actual mobilized resistance may be lower if the soil has not fully reached the plastic state.
- **Spacing ratio**: The theory is applicable for $S/D$ between approximately **2 and 8**. Below $S/D \approx 2$, the piles act more like a continuous retaining wall. Above $S/D \approx 8$, soil arching between piles becomes negligible and the method overestimates the force. XSLOPE's model checks warn before a run when an auto-computed row sits outside this band.
- **Originally derived for horizontal ground**: The overburden term $\gamma z$ in $p(z)$ assumes the vertical stress at depth $z$ equals $\gamma z$, which is exact only for horizontal ground. On a slope face, the actual vertical stress at the pile location is less than $\gamma z$ because the soil column above is truncated by the slope geometry. This means the method can overestimate the overburden contribution for piles on the slope face, particularly near the toe where the soil column is shallowest relative to a horizontal surface at the same elevation. For piles behind the crest on level ground, the approximation is accurate. In practice, the theory is routinely applied to slopes and the overestimation is generally accepted as conservative (it increases the computed pile resistance, not the driving forces). XSLOPE uses the vertical depth from the ground surface at the pile location to the failure surface, consistent with standard practice in commercial software (Slide2, SLOPE/W).
- **Upper bound on soil force**: The computed $H$ represents the soil's capacity to push on the pile. The actual pile resistance used in the LEM is the **lesser** of the Ito-Matsui soil force and the pile's structural shear/bending capacity. See [Structural Capacity Checks](#structural-capacity-checks) below for how XSLOPE enforces this limit when $V_{\text{cap}}$ and $M_{\text{cap}}$ are provided.


## Structural Capacity Checks

The pile resistance used in LEM should not exceed the structural capacity of the pile. Two structural failure modes are checked when the optional $V_{\text{cap}}$ and $M_{\text{cap}}$ columns are provided in the `piles` sheet:

- **Shear capacity** ($V_{\text{cap}}$): The maximum lateral shear force that the pile cross-section can resist. For concrete piles, this is governed by the concrete and steel reinforcement; for steel piles, by the web and flange dimensions.
- **Moment capacity** ($M_{\text{cap}}$): The maximum bending moment the pile can resist. The limiting lateral force from bending is $M_{\text{cap}} / L_m$, where $L_m$ is the moment arm from the pressure centroid to the failure surface.

Both $V_{\text{cap}}$ and $M_{\text{cap}}$ are properties of a **single pile** (not per unit width). The capacity check compares them against the per-pile force $F_{\text{pile}}$, not the per-unit-width force $H$. The pile forces xslope reports back — and the FEM pile-shear colorbar — are per unit width of slope; multiply by the spacing $S$ to recover the per-pile force for comparison against the single-pile $V_{\text{cap}}$ / $M_{\text{cap}}$.

### Capacity Check Procedure

The capacity check applies regardless of how the pile force was obtained, but the details differ between the two cases.

**Common steps** (both cases):

1. If $V_{\text{cap}}$ is provided: $\;F_{\text{pile}} = \min(F_{\text{pile}},\; V_{\text{cap}})$
2. If $M_{\text{cap}}$ is provided: $\;F_{\text{pile}} = \min(F_{\text{pile}},\; M_{\text{cap}} / L_m)$
3. Convert back to per-unit-width: $\;H = F_{\text{pile}} / S$

The capped $H$ is then used in the slice equilibrium equations.

### Case 1: Ito & Matsui Auto-Computed $H$

When $H$ is left blank and $D$ and $S$ are provided, XSLOPE computes $F_{\text{pile}}$ by integrating the Ito & Matsui pressure distribution $p(z) = c \cdot A_1 + \gamma z \cdot A_2$ from the ground surface to the failure surface. Because the full pressure distribution is known, XSLOPE also computes the **exact moment arm** $L_m$ from the centroid of that distribution:

>$L_m = \dfrac{\displaystyle\int_0^{z_f} (z_f - z)\, p(z)\, dz}{F_{\text{pile}}}$

The integration is performed piecewise over each soil layer (the same segments used for the force calculation). Some limiting cases:

- **Uniform pressure** ($c > 0$, $\gamma = 0$): $L_m = z_f / 2$
- **Triangular pressure** ($c = 0$, $\gamma > 0$): $L_m = z_f / 3$
- **General** $c$-$\phi$ **soil**: $z_f / 3 < L_m < z_f / 2$

The controlling design value is:

>$H = \dfrac{1}{S}\min(F_{\text{Ito-Matsui}},\; V_{\text{cap}},\; M_{\text{cap}} / L_m)$

The summary output reports the soil force, each capacity check with $[\text{GOVERNS}]$ or $[\text{OK}]$ status, and the capped values if the structural capacity controls.

### Case 2: User-Specified $H$

When the user provides $H$ directly, the per-pile force is computed as $F_{\text{pile}} = H \times S$. The $V_{\text{cap}}$ check is straightforward — it is a direct comparison of $F_{\text{pile}}$ against the shear capacity.

For the $M_{\text{cap}}$ check, the pressure distribution behind the pile is unknown, so XSLOPE cannot compute $L_m$ from integration. Instead, it uses a default of:

>$L_m = z_f / 3$

This corresponds to a triangular pressure distribution (linearly increasing with depth), which gives the **largest** $M_{\text{cap}} / L_m$ and therefore the **least restrictive** cap on $F_{\text{pile}}$. If the actual pressure distribution is more uniform (top-heavy), the true $L_m$ would be larger and the $M_{\text{cap}}$ check would be more restrictive. Users who know their pressure distribution can account for this by pre-computing the capped force externally:

>$H = \dfrac{1}{S}\min(F_{\text{soil}},\; V_{\text{cap}},\; M_{\text{cap}} / L_m) \qquad \text{(enter this value directly)}$

### Summary of Differences

| | $F_{\text{pile}}$ source | $L_m$ for $M_{\text{cap}}$ check | Summary detail |
|---|---|---|---|
| **Ito & Matsui** | From integration of $p(z)$ | Exact, from pressure centroid | Full Ito & Matsui summary with capacity check |
| **User-specified** $H$ | $H \times S$ | Default $z_f / 3$ | Capacity check only (no Ito & Matsui summary) |

### Run-Summary Output

A limit equilibrium run on a model with pile rows prints one line per row
with the factor of safety, so the applied forces are visible without
generating a report: the Ito & Matsui soil force per pile, the capacity
that governed (`bending governs (Mcap/Lm, Lm = ...)` or
`shear governs (Vcap)`), and the force applied per unit width, with a total
when more than
one row contributes. A row with a user-specified $H$ prints its stated
value, or — when a capacity binds it — the governing check and the applied
force. The full per-slice accounting remains in the
[Analysis Report](../studio/reports.md).

## LEM vs. FEM Pile Modeling

Which analysis suits which member is set out under [LEM vs FEM](overview.md#lem-vs-fem). The two are compared here
on the slopes where both have been run.

For a **continuous member** the beam formulation is an exact description rather than an idealization, its $EA$ and $EI$ already are per unit width, and it returns internal actions that a limit equilibrium analysis cannot produce at all. It is compared with GeoStudio's SIGMA/W sheet pile wall example: XSLOPE gives 1.048 without the wall and 1.691 with it, against about 1.025 and 1.4 from SIGMA/W, a gap with the wall that remains — see [the SIGMA/W wall benchmark](../verification/geostudio.md#sigmaw-wall) and [Applicability](fem.md#applicability-continuous-walls-and-discrete-pile-rows) in the FEM pile documentation.

For a **discrete row** the limit equilibrium analysis models the actual mechanism, with the soil moving between the piles. The size of the difference is measured on the pile model of [Tutorial LEM-12](../tutorials/lem12_piles.md) and [FEM-4](../tutorials/fem04_piles.md) — a 1:1 slope in c = 200 psf, $\phi$ = 20° soil with two rows of 2 ft drilled shafts at 6 ft spacing — which is solved by both engines on the same section, soil and pile rows:

| | Without piles | With piles | Credit for the row |
|---|---|---|---|
| **LEM** (Spencer) | 1.149 | 1.842 | ×1.60 |
| **FEM** (SSRM) | 1.137 | 1.363 | ×1.20 |

Without the piles the two engines agree to about 1%, so the difference in the second column comes from the pile row alone. The row raises the factor of safety by a factor of 1.60 in the limit equilibrium analysis and 1.20 in the finite element analysis, a difference in the quantity being designed that is far larger than rounding.

Neither of those is a three-dimensional answer, and the direction of the error is only known where a three-dimensional reference exists. Cai & Ugai (2000) analyzed a pile-stabilized slope with a shear-strength-reduction finite element model that meshes the individual piles, the soil between them and the slip interfaces on each pile's surface. XSLOPE solves the same slope with both of its engines:

| Case | XSLOPE SSRM (2D beam) | Cai & Ugai 3D FE |
|---|---|---|
| No pile | 1.136 | 1.14 (−0.4%) |
| Pile at $D_1/D$ = 3, free head | 1.578 | 1.36 (+16.0%) |
| Pile, head rotation restrained | 1.594 | 1.45 (+9.9%) |

The unpiled case agrees to 0.4%, so the differences in the other two rows come from the pile. With the row in place the plane-strain model gives higher values: it credits the row with multiplying the unreinforced factor of safety by 1.389 where the three-dimensional model credits 1.193. On the same slope a Bishop search with the Ito & Matsui force gives 1.451 against the paper's own limit-equilibrium value of 1.37 and Slide2's 1.43, a credit of 1.269. Both two-dimensional credits stand above the three-dimensional one, the beam's by 0.196 and the limit-equilibrium search's by 0.076, so neither recovers it and the limit-equilibrium search lands nearer. This is the only benchmark with a published three-dimensional answer. Both comparisons are quantified in [VP106](../verification/rocscience.md#vp106) and [the VP106 finite-element diagnostic](../verification/rocscience.md#vp106-fem).

## References

Ito, T., & Matsui, T. (1975). Methods to estimate lateral force acting on stabilizing piles. *Soils and Foundations*, 15(4), 43-59.

Poulos, H.G. (1995). Design of reinforcing piles to increase slope stability. *Canadian Geotechnical Journal*, 32(5), 808-818.

Hassiotis, S., Chameau, J.L., & Gunaratne, M. (1997). Design method for stabilization of slopes with piles. *Journal of Geotechnical and Geoenvironmental Engineering*, 123(4), 314-323.
