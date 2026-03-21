# Piles and Concrete Piers in Slope Stability

## Introduction

Piles and concrete piers are rigid structural elements used to stabilize slopes by resisting lateral soil movement through shear and bending at the failure surface. Unlike flexible reinforcement (geotextiles, soil nails) which deforms with the soil and provides tension along its axis, piles act as passive structural inclusions that the sliding soil mass pushes against. The pile develops a lateral resisting force at the point where the failure surface intersects the pile, directly opposing the driving forces that cause slope instability.

Pile stabilization is widely used in practice for:

- Remediation of active landslides
- Stabilization of slopes for new construction (roads, buildings, retaining structures)
- Protection of existing infrastructure from slope movement
- Increasing the factor of safety of marginally stable slopes

The following figure shows the basic concept. A row of piles is installed through the sliding mass and embedded into stable ground below the failure surface. As the soil mass attempts to move, it pushes against the piles, which resist through their lateral stiffness and embedment in the stable zone.

## How Piles Differ from Flexible Reinforcement

In limit equilibrium analysis, flexible reinforcement (geotextiles, soil nails) is modeled as a tensile force $P$ acting **parallel to the slice base** in the direction that resists sliding. This is appropriate because flexible reinforcement deforms with the soil and its resistance develops along the reinforcement axis.

Piles behave fundamentally differently:

| Property | Flexible Reinforcement | Pile / Pier |
|----------|----------------------|-------------|
| **Deformation** | Deforms with the soil | Resists soil deformation through rigidity |
| **Force direction** | Along reinforcement axis (parallel to base) | At the failure surface intersection, typically horizontal or at a user-specified angle |
| **Force mechanism** | Tension along the reinforcement | Lateral shear and bending resistance |
| **Force components** | Tangential to base only ($P$) | Both normal and tangential to base ($H\sin(\alpha - \theta_p)$ and $H\cos(\alpha - \theta_p)$) |
| **Compression** | Cannot resist compression | Resists both tension and compression |
| **Moment contribution** | None (flexible) | Resisting moment about circle center (moment-based methods) |

These differences require a distinct formulation for incorporating pile forces into each limit equilibrium method. The pile force $H$ is applied at the point where the pile intersects the failure surface and resolved into components appropriate for each method's equilibrium equations.


## Pile Force in Limit Equilibrium Analysis

### Force Definition

In the limit equilibrium framework, a pile is characterized by:

- A **force magnitude** $H$ (per unit width of slope, i.e., force/length) acting at the point where the pile intersects the failure surface
- A **force angle** $\theta_p$ measured from horizontal (positive = counterclockwise/upward, default = 0)

The force is decomposed into horizontal and vertical components:

>$H_h = H \cos\theta_p \qquad \text{(horizontal)}$

>$H_v = H \sin\theta_p \qquad \text{(vertical, positive upward)}$

When $\theta_p = 0$ (the default and most common case), the force is purely horizontal: $H_h = H$ and $H_v = 0$.

### Force Resolution on the Slice

The pile force acts at point $e$ on the failure surface where the pile intersects the slice base. Relative to the slice base (inclined at angle $\alpha$), the force resolves into:

>**Normal to base** (increases effective stress): $H \sin(\alpha - \theta_p)$

>**Tangential to base** (resists sliding): $H \cos(\alpha - \theta_p)$

The normal component increases the effective normal stress on the failure surface, which in turn increases frictional resistance. The tangential component directly opposes the sliding force. Both components contribute to stability.

### Moment Contribution

For methods that use moment equilibrium about a circle center $(X_o, Y_o)$, the pile force creates a resisting moment. The horizontal and vertical components of $H$ each contribute through their respective moment arms:

>$M_{\text{pile}} = H \cos\theta_p \cdot (Y_o - y_e) + H \sin\theta_p \cdot (x_e - X_o)$

where $(x_e, y_e)$ is the pile-failure surface intersection point. When $\theta_p = 0$, this reduces to $M_{\text{pile}} = H(Y_o - y_e)$.

The pile force is a **known applied force** and is **not** divided by the factor of safety $F$. This is because $H$ represents the structural resistance of the pile, not a soil strength parameter. It appears in the denominator of the factor of safety equation (reducing driving forces) rather than in the numerator (which contains soil shear strength divided by $F$).

### Per-Unit-Width Convention

All forces in 2D limit equilibrium analysis are expressed per unit width of slope (perpendicular to the cross-section). If a row of piles has individual capacity $H_{\text{single}}$ at center-to-center spacing $S$, the equivalent force per unit width is:

>$H = \dfrac{H_{\text{single}}}{S}$

### Non-Circular Failure Surfaces

The pile force formulation for Janbu, Corps of Engineers, Lowe-Karafiath, and Spencer's method uses only per-slice quantities ($\alpha$, force components) and has no dependence on circle geometry. These methods work with any failure surface shape, and the pile terms carry over without modification. The only circle-dependent pile terms appear in OMS and Bishop (the moment term), but those methods inherently require circular surfaces.

### Integration with LEM Methods

The pile force $H$ at angle $\theta_p$ is incorporated into each limit equilibrium method supported by XSLOPE. The specific modifications for each method are presented in the respective method documentation pages:

- [**OMS**](oms.md): $H\sin(\alpha-\theta_p)$ added to $N'$; pile moment terms added to the denominator
- [**Bishop**](bishop.md): $-H\sin\theta_p$ enters vertical equilibrium for $N'$; pile moment terms added to the denominator
- [**Janbu**](janbu.md): $H\sin(\alpha-\theta_p)$ added to $N'$; $H\cos(\alpha-\theta_p)$ subtracted from the denominator
- [**Force Equilibrium** (Corps of Engineers, Lowe-Karafiath)](force_eq.md): $-H\cos\theta_p$ added to horizontal equilibrium ($b_0$); $-H\sin\theta_p$ added to vertical equilibrium ($b_1$)
- [**Spencer**](spencer.md): $H\cos\theta_p$ added to $F_h$; $H\sin\theta_p$ added to $F_v$; moment terms added to $M_o$

In all methods, when $\theta_p = 0$, the equations reduce to the simpler horizontal-force-only case.


## Stabilizing Piles vs. Load-Bearing Piles

The discussion above focuses on **stabilizing piles** (also called passive piles) — piles installed specifically to resist lateral soil movement and improve slope stability. However, piles near slopes may also serve as **load-bearing piles** that carry structural loads (buildings, bridges, retaining walls) through the soil to deeper bearing strata. The treatment of these two cases in slope stability analysis is fundamentally different.

### Stabilizing (Passive) Piles

Stabilizing piles are the primary focus of the pile implementation in XSLOPE. These piles:

- Are installed through the sliding mass and embedded in stable ground below the failure surface
- Resist lateral soil movement through shear and bending
- Provide a lateral force $H$ at the failure surface that is incorporated into the LEM equations as described above
- Their own self-weight is typically negligible relative to the soil mass and is ignored

### Load-Bearing Piles

Load-bearing piles carry structural loads (vertical forces from foundations) and transfer them to the subsurface through a combination of **skin friction** along the pile shaft and **end bearing** at the pile tip. The key question for slope stability is: does the structural load contribute to the driving forces on the failure surface?

The answer depends on where the pile tip is relative to the failure surface:

**Case 1: Pile tip below the failure surface (most common)**

This is the usual design intent for load-bearing piles near slopes. The pile transfers the structural load through skin friction and end bearing to stable ground **below** the failure surface. In this case:

- The structural load effectively **bypasses** the sliding mass — it is delivered to stable ground below the failure surface and does not contribute to driving forces
- The structural load should be **omitted** from the slope stability model (i.e., do not apply it as a surcharge on the slope surface)
- The pile may still provide lateral resistance to sliding, which can be modeled separately as a stabilizing pile force $H$
- This is the approach recommended by FHWA, AASHTO, and used by commercial slope stability software (SLOPE/W, Slide2)

**Case 2: Pile tip above the failure surface**

If the pile tip is entirely within the sliding mass (a friction pile in weak soil, for example), the entire pile and its load are part of the sliding mass. In this case:

- The structural load **does** contribute to driving forces and should be included in the analysis
- Standard practice is to apply the structural load as a **distributed surface surcharge** using the distributed loads (`dloads`) sheet in the XSLOPE input template
- This is slightly conservative because it places all the weight at the surface rather than distributing it with depth through skin friction, but the conservatism is generally small and accepted in practice

**Case 3: Failure surface intersects the pile mid-shaft**

This is the most complex scenario. Some of the structural load has been transferred above the failure surface through skin friction (and thus loads the sliding mass), while some has been transferred below (and bypasses it). In practice:

- Most engineers apply the **conservative upper bound** — treat the full structural load as a surface surcharge, assuming all load is within the sliding mass
- A refined approach would compute the cumulative skin friction transferred above the failure surface (requiring a load-transfer or t-z analysis), but this is rarely done in practice because the complexity is not justified by the improvement in accuracy
- Modeling the depth-proportional weight distribution within the method of slices is not practical — the slice formulation assigns a single weight to each slice based on its cross-sectional area and unit weight

**Summary**: For load-bearing piles near slopes, the recommended approach in XSLOPE is:

1. If the pile tip is below the failure surface, omit the structural load from the slope stability model
2. If the pile tip is above the failure surface, apply the structural load as a distributed surface load on the `dloads` sheet
3. If the pile also provides lateral resistance to sliding, model that separately as a stabilizing pile force $H$

The existing distributed load capability in XSLOPE handles the surcharge case. No additional code is needed for load-bearing piles — only clear guidance on when and how to apply the structural load.


## Determining the Pile Force $H$

### User-Specified Force

The simplest approach is for the user to specify $H$ directly based on external analysis. The pile force may come from:

- p-y curve analysis software (e.g., LPILE, GROUP, RSPile)
- Structural analysis of the pile section
- Published design charts or empirical correlations
- Full 3D finite element analysis

When using user-specified forces, the user enters $H$ (per unit width) and $\theta_p$ in the `piles` sheet of the input template. This approach gives the user full control and is appropriate when detailed pile analysis has already been performed.

### Ito & Matsui (1975) Theory

The Ito & Matsui method is the most widely used closed-form approach for computing the lateral force that soil exerts on passive stabilizing piles. It models the soil between adjacent piles as being in a state of **plastic equilibrium** — the soil deforms plastically as it squeezes between the piles, like material flowing through a constriction. Using Mohr-Coulomb plasticity theory, Ito & Matsui derived closed-form equations for the lateral pressure on the piles as a function of depth.

#### Setup and Notation

Consider a row of piles embedded in a slope:

- $D$ = pile diameter (or width of the pile cross-section)
- $S$ = center-to-center spacing between piles
- $D_1 = S - D$ = clear spacing between adjacent pile faces
- $z$ = depth below the ground surface
- $z_f$ = depth from the ground surface to the failure surface at the pile location
- Soil properties: cohesion $c$, friction angle $\varphi$, unit weight $\gamma$
- Passive earth pressure coefficient: $N_\varphi = \tan^2\!\left(45° + \dfrac{\varphi}{2}\right)$

The theory applies to the portion of the pile **above** the failure surface — this is the zone where soil is actively moving and pushing against the pile.

#### General $c$-$\varphi$ Soil

For a soil with both cohesion and friction ($c > 0$, $\varphi > 0$), the distributed lateral force $p(z)$ (force per unit length of pile) at depth $z$ is:

>$p(z) = c \cdot A_1 + \gamma z \cdot A_2$

where $A_1$ and $A_2$ are dimensionless arching coefficients. To define them, first compute the geometric amplification factor:

>$G = \left(\dfrac{D_1 + D}{D_1}\right)^{N_\varphi \tan\varphi + N_\varphi^{1/2} - 1}$

This factor captures the plastic flow amplification as soil squeezes between the piles. It increases as the piles get closer together ($D_1/D$ decreases) and as $\varphi$ increases.

The arching coefficients are then:

>$A_1 = D_1 \left[\dfrac{N_\varphi^{1/2} \cdot G}{N_\varphi \tan\varphi - 1}\left(N_\varphi^{1/2} \tan\varphi + N_\varphi^{1/2} - 1\right) - \dfrac{2 N_\varphi^{1/2} \tan\varphi + 1}{N_\varphi \tan\varphi - 1}\right] + \dfrac{D}{2}\left(N_\varphi^{1/2} \cdot G - 1\right)$

>$A_2 = D_1 \left[\dfrac{N_\varphi \cdot G - 1}{N_\varphi \tan\varphi - 1}\right] + \dfrac{D}{2}\left(N_\varphi \cdot G - 1\right)$

The full derivation is given in Ito & Matsui (1975). These expressions emerge from solving the plastic equilibrium equations for soil flowing between two rigid boundaries under Mohr-Coulomb failure.

Key behavior of $p(z)$:

- **Increases with depth** through the $\gamma z$ term — deeper soil mobilizes more pressure against the pile
- **Increases as $D_1/D$ decreases** (closer piles = more arching = more force per pile)
- **Increases with $\varphi$** — higher friction angle produces stronger soil arching between piles
- **Increases with $c$** — cohesion contributes a constant (depth-independent) component

#### Cohesionless Soil ($c = 0$)

For a purely frictional soil with $c = 0$, the cohesion term vanishes and the lateral pressure is:

>$p(z) = \gamma z \cdot A_2$

where $A_2$ is the same expression as above. The pressure increases linearly from zero at the ground surface.

#### Undrained Clay ($\varphi = 0$)

For a purely cohesive (undrained) soil with $\varphi = 0$, the general $c$-$\varphi$ expressions become indeterminate because $N_\varphi = 1$ and $\tan\varphi = 0$, causing the denominator $N_\varphi \tan\varphi - 1 = -1$. However, taking the appropriate limit or deriving the solution independently for $\varphi = 0$ yields simplified expressions:

>$p(z) = c_u \cdot D_1 \left[2\left(\dfrac{D_1 + D}{D_1}\right)\left(\dfrac{D_1 + D}{2 D_1} + \dfrac{\pi}{4}\right) - 1\right] + \gamma z \cdot D \left(\dfrac{D_1 + D}{D_1} - 1\right)$

where $c_u$ is the undrained shear strength. The first term represents the cohesion contribution (constant with depth) and the second term represents the overburden contribution (linear with depth). For undrained clay, the overburden term is typically small relative to the cohesion term.

#### Total Force per Pile

The total lateral force on a single pile is obtained by integrating $p(z)$ from the ground surface down to the failure surface depth $z_f$:

>$F_{\text{pile}} = \int_0^{z_f} p(z) \, dz$

Since $p(z)$ is linear in $z$ within a homogeneous layer, the integration is straightforward:

>$F_{\text{pile}} = c \cdot A_1 \cdot z_f + \gamma \cdot A_2 \cdot \dfrac{z_f^2}{2}$

#### Force per Unit Width

The 2D plane-strain equivalent force used in LEM (per unit width of slope) is:

>$H = \dfrac{F_{\text{pile}}}{S}$

This is the value entered (or computed) for the pile force in the slope stability analysis.

#### Multi-Layer Soils

When the pile passes through multiple material zones above the failure surface (common in practice), the integration is performed piecewise. For each layer $j$ with properties $c_j$, $\varphi_j$, $\gamma_j$ between depths $z_{\text{top},j}$ and $z_{\text{bot},j}$:

>$F_j = c_j \cdot A_{1,j} \cdot (z_{\text{bot},j} - z_{\text{top},j}) + \gamma_j \cdot A_{2,j} \cdot \dfrac{z_{\text{bot},j}^2 - z_{\text{top},j}^2}{2}$

The total force per pile is the sum over all layers:

>$F_{\text{pile}} = \sum_j F_j$

Note that $A_{1,j}$ and $A_{2,j}$ must be recomputed for each layer since $\varphi$ may differ between layers. The pile geometry ($D$, $D_1$) remains the same for all layers.

#### Computation at Each Trial Surface

An important characteristic of the Ito & Matsui calculation is that **$H$ depends on the failure surface location**. A deeper failure surface means more soil above it pushing on the pile, giving a higher $H$. Therefore, $H$ should be recomputed for each trial failure surface during an automated search. Since the computation involves only closed-form expressions and simple integration, it is essentially instantaneous and adds no meaningful computational cost.

In XSLOPE, when $H$ is left blank in the `piles` sheet but the pile diameter $D$ and spacing $S$ are provided, the Ito & Matsui force is computed automatically at slice generation time for each trial surface. If the user provides an explicit $H$ value, that value is used instead (override mode).

#### Applicability and Limitations

The Ito & Matsui method has the following characteristics and limitations:

- **Rigid pile assumption**: The theory assumes piles do not deflect significantly. This is conservative for flexible piles, which mobilize less soil pressure than rigid piles.
- **Plastic flow assumption**: The method gives an **upper bound** on the soil's capacity to push on the pile. The actual mobilized resistance may be lower if the soil has not fully reached the plastic state.
- **Spacing ratio**: The theory is applicable for $S/D$ between approximately **2 and 8**. Below $S/D \approx 2$, the piles act more like a continuous retaining wall. Above $S/D \approx 8$, soil arching between piles becomes negligible and the method overestimates the force.
- **Originally derived for horizontal ground**: Applying the method to sloping ground introduces some approximation, though in practice the theory is routinely used for slopes. The error is small for typical slope angles.
- **No structural check**: The computed $H$ represents the soil's capacity to push on the pile. The actual pile resistance used in the LEM should be the **lesser** of the Ito-Matsui soil force and the pile's structural shear/bending capacity:

>$H_{\text{design}} = \min(H_{\text{Ito-Matsui}},\; H_{\text{structural}})$

### Structural Capacity

The pile resistance used in LEM should not exceed the structural capacity of the pile. Two structural failure modes should be checked:

- **Shear capacity**: The nominal shear strength $V_n$ of the pile cross-section. For concrete piles, this is governed by the concrete and steel reinforcement; for steel piles, by the web and flange dimensions.
- **Moment capacity**: The pile must be able to resist the bending moment that develops from the soil pressure distribution. The limiting lateral force from bending is approximately $M_p / L_m$, where $M_p$ is the plastic moment capacity and $L_m$ is the moment arm from the point of maximum moment to the fixity point below the failure surface.

The controlling design value is:

>$H = \min(H_{\text{soil}},\; V_n,\; M_p / L_m)$

where $H_{\text{soil}}$ is from Ito & Matsui or other soil-pile interaction analysis.

### Soil Arching Between Piles

A critical aspect of pile-stabilized slopes is the three-dimensional soil arching that develops between adjacent piles. As the sliding soil mass pushes against the pile row, stress concentrations develop around each pile, and the soil "arches" between piles in a manner analogous to arching above a tunnel. This is the mechanism captured by the Ito & Matsui theory.

The effectiveness of soil arching depends on:

- **Pile spacing**: Closer spacing produces stronger arching and higher force per pile. The optimal spacing balances structural efficiency (fewer piles) against arching effectiveness (closer piles).
- **Soil strength**: Stronger soils develop more effective arching. In very weak soils (soft clay), arching may be minimal and the soil may flow between the piles without mobilizing significant resistance.
- **Pile rigidity**: Rigid piles provide fixed points for arch development. Flexible piles may deflect enough to reduce arching effectiveness.

For design purposes, $S/D$ ratios of 3 to 6 are typical for slope stabilization applications.


## Typical Parameter Values

The following tables provide typical ranges of pile and pier parameters for use in slope stability analysis. These values are intended as guidance for filling in the `piles` sheet of the input template and for preliminary design estimates.

### Pile Types and Typical Dimensions

| Pile Type | Typical Diameter/Width | Typical Length | Notes |
|-----------|----------------------|----------------|-------|
| Steel H-pile | 200-400 mm (8-14 in) | 6-30 m (20-100 ft) | Wide-flange sections; high strength-to-weight ratio |
| Steel pipe pile | 300-900 mm (12-36 in) | 6-40 m (20-130 ft) | Open-ended or closed-ended; often concrete-filled |
| Precast concrete pile | 250-600 mm (10-24 in) | 6-25 m (20-80 ft) | Square or octagonal cross-section |
| Drilled shaft (caisson) | 450-3000 mm (18-120 in) | 3-60 m (10-200 ft) | Cast-in-place; larger diameters common for slope stabilization |
| Concrete pier | 600-1500 mm (24-60 in) | 3-15 m (10-50 ft) | Often used for slope stabilization; rectangular or circular |
| Micropile | 150-300 mm (6-12 in) | 6-30 m (20-100 ft) | Drilled and grouted; used in tight access or existing structures |
| Timber pile | 200-400 mm (8-16 in) | 6-20 m (20-65 ft) | Tapered; limited to lighter loads |

### Typical Spacing

| Application | Typical $S/D$ Ratio | Typical Spacing $S$ | Notes |
|-------------|---------------------|---------------------|-------|
| Slope stabilization | 3-6 | 1.5-6 m (5-20 ft) | Closer spacing = more arching between piles |
| Retaining structures | 2-4 | 1-4 m (3-12 ft) | Often soldier piles with lagging |
| Ito & Matsui applicability | 2-8 | -- | Theory assumes plastic flow between piles |

### Typical Section Properties

| Section | $A$ (Area) | $I$ (Moment of Inertia) |
|---------|-----------|------------------------|
| Solid circular, $D$ = 0.6 m (24 in) | 0.283 m$^2$ (452 in$^2$) | 6.36 x 10$^{-3}$ m$^4$ (16,286 in$^4$) |
| Solid circular, $D$ = 0.9 m (36 in) | 0.636 m$^2$ (1,018 in$^2$) | 3.22 x 10$^{-2}$ m$^4$ (82,448 in$^4$) |
| Solid circular, $D$ = 1.2 m (48 in) | 1.131 m$^2$ (1,810 in$^2$) | 1.02 x 10$^{-1}$ m$^4$ (260,576 in$^4$) |
| HP 14x117 (steel H-pile) | 0.022 m$^2$ (34.4 in$^2$) | 4.43 x 10$^{-4}$ m$^4$ (1,063 in$^4$) |
| HP 12x84 (steel H-pile) | 0.016 m$^2$ (24.6 in$^2$) | 2.18 x 10$^{-4}$ m$^4$ (524 in$^4$) |
| Pipe pile, $D$ = 0.6 m, $t$ = 12 mm | 0.022 m$^2$ (34.6 in$^2$) | 9.7 x 10$^{-4}$ m$^4$ (2,330 in$^4$) |

For solid circular sections: $A = \pi D^2 / 4$, $\quad I = \pi D^4 / 64$.

### Material Properties

| Material | $E$ (kPa) | $E$ (psf) | $\nu$ |
|----------|-----------|-----------|-------|
| Structural steel | 2.0 x 10$^8$ | 4.18 x 10$^9$ | 0.3 |
| Reinforced concrete ($f'_c$ = 4000 psi) | 2.5 x 10$^7$ | 5.2 x 10$^8$ | 0.2 |
| Reinforced concrete ($f'_c$ = 5000 psi) | 2.8 x 10$^7$ | 5.8 x 10$^8$ | 0.2 |
| Reinforced concrete ($f'_c$ = 6000 psi) | 3.0 x 10$^7$ | 6.3 x 10$^8$ | 0.2 |
| Timber (Douglas Fir) | 1.2 x 10$^7$ | 2.5 x 10$^8$ | 0.3 |

Concrete modulus computed as $E_c = 57{,}000 \sqrt{f'_c}$ (psi) per ACI 318.

### Typical Lateral Resistance $H$

Lateral resistance depends heavily on soil conditions, pile geometry, and embedment depth. The following ranges are approximate for preliminary estimates only.

| Soil Type | Pile Type | Typical $H$ per pile | Notes |
|-----------|-----------|---------------------|-------|
| Stiff clay ($c$ = 50-100 kPa) | Drilled shaft, $D$ = 0.9 m, $S$ = 3 m | 100-400 kN/m (7-27 kip/ft) | Per unit width = $H_{\text{pile}} / S$ |
| Medium dense sand ($\varphi$ = 30-35 deg) | Steel H-pile, $S$ = 2 m | 50-200 kN/m (3-14 kip/ft) | Depends on depth above failure surface |
| Soft clay ($c$ = 15-30 kPa) | Concrete pier, $D$ = 1.2 m, $S$ = 3 m | 30-100 kN/m (2-7 kip/ft) | Lower bound; may govern over structural capacity |
| Weathered rock | Drilled shaft, $D$ = 0.9 m, $S$ = 3 m | 200-800 kN/m (14-55 kip/ft) | High capacity but expensive installation |

These values are for preliminary guidance only. Actual $H$ should be determined from Ito & Matsui theory, p-y analysis, or structural analysis of the pile.


## Input Template

Pile data is entered in the **piles** sheet of the Excel input template. Each row defines one pile with the following columns:

| Column | Field | Required | Description |
|--------|-------|----------|-------------|
| A | Label | No | Name/identifier for the pile |
| B | $x_1$ | Yes | X-coordinate of pile top |
| C | $y_1$ | Yes | Y-coordinate of pile top |
| D | $x_2$ | Yes | X-coordinate of pile tip (bottom) |
| E | $y_2$ | Yes | Y-coordinate of pile tip (bottom) |
| F | $H$ | LEM | Pile force magnitude per unit width of slope (kN/m or lb/ft) |
| G | $\theta_p$ | No | Force angle from horizontal in degrees (default 0; positive = upward) |
| H | $D_{\text{pile}}$ | No | Pile diameter or width |
| I | $S$ | No | Center-to-center spacing |
| J | $E$ | FEM | Young's modulus of pile material |
| K | $I$ | FEM | Moment of inertia (computed from $D_{\text{pile}}$ if omitted) |
| L | $Area$ | FEM | Cross-sectional area (computed from $D_{\text{pile}}$ if omitted) |

The two-endpoint geometry $(x_1, y_1)$ to $(x_2, y_2)$ supports both **vertical** and **battered (inclined) piles**. Vertical piles are simply the case where $x_1 = x_2$.

The force magnitude $H$ is per unit width of slope, consistent with the 2D plane-strain assumption used throughout XSLOPE. If the user has a row of piles at spacing $S$ with individual capacity $H_{\text{single}}$, the input value is $H = H_{\text{single}} / S$.

If $H$ is left blank and $D_{\text{pile}}$ and $S$ are provided, XSLOPE can auto-compute $H$ using the Ito & Matsui method (see above). If $H$ is provided, it is used directly regardless of other columns.

Reading stops when column B ($x_1$) is empty, following the same pattern as the reinforcement sheet.

Columns J-L ($E$, $I$, $Area$) are used only for FEM beam element modeling (see [FEM Piles](../fem/piles.md)). If $D_{\text{pile}}$ is provided and $I$/$Area$ are omitted, a solid circular section is assumed: $A = \pi D^2/4$, $I = \pi D^4/64$.


## References

Ito, T., & Matsui, T. (1975). Methods to estimate lateral force acting on stabilizing piles. *Soils and Foundations*, 15(4), 43-59.

Poulos, H.G. (1995). Design of reinforcing piles to increase slope stability. *Canadian Geotechnical Journal*, 32(5), 808-818.

FHWA. (2009). *Design and Construction of Driven Pile Foundations*. Publication No. FHWA-NHI-05-042/043, Federal Highway Administration.

Hassiotis, S., Chameau, J.L., & Gunaratne, M. (1997). Design method for stabilization of slopes with piles. *Journal of Geotechnical and Geoenvironmental Engineering*, 123(4), 314-323.
