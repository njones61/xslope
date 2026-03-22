# Plan: Pile and Concrete Pier Support in Slope Stability Analysis


## 1. Problem Overview

Flexible reinforcement (geogrids, soil nails) resists slope movement through **tension** — the reinforcement deforms with the soil and the force acts along the reinforcement orientation. In xslope, this is handled by resolving the tension force as a component **P parallel to the slice base**, which directly reduces driving forces.

Piles and concrete piers are fundamentally different. They are **rigid structural elements** that resist soil movement through **lateral shear and bending** at the failure surface intersection. Key differences from flexible reinforcement:

1. **Force direction**: A pile provides a resisting force $H$ at the point where the failure surface crosses the pile, typically horizontal but optionally at a user-specified angle $\theta_p$. Flexible reinforcement provides a force along its own (typically inclined) axis.

2. **Force resolution**: The pile force $H$ at angle $\theta_p$ must be resolved into **both** normal and tangential components relative to the slice base:
   - Tangential (parallel to base): $H \cos(\alpha - \theta_p)$ — directly resists sliding
   - Normal (perpendicular to base): $H \sin(\alpha - \theta_p)$ — increases normal stress, boosting frictional resistance

   This is different from flexible reinforcement where $P$ acts purely parallel to the base.

3. **Bending moment**: A pile with flexural stiffness contributes a **resisting moment** at the failure surface intersection. This matters for methods that satisfy moment equilibrium (OMS, Bishop, Spencer) but not for force-only methods (Janbu, Corps of Engineers, Lowe-Karafiath).

4. **Soil arching**: Soil arches between piles in a row, creating a 3D effect that is typically simplified to a 2D equivalent force per unit width of slope.

5. **Two failure modes**: The available resistance is the **lesser** of:
   - Soil flow-around capacity (soil squeezes between piles)
   - Structural capacity of the pile itself


## 2. Determining Pile Properties

There are several established approaches for computing the lateral force $H$ that a pile provides at the failure surface:

### 2.1 Ito & Matsui (1975) Theory

The most widely used closed-form method for passive piles in slopes. It models the soil between adjacent piles as being in a state of **plastic equilibrium** — the soil deforms plastically as it "squeezes" between the piles, like material flowing through a constriction. Using Mohr-Coulomb plasticity, Ito & Matsui derived closed-form equations for the lateral pressure on the piles as a function of depth.

#### Setup

Consider a row of piles embedded in a slope:
- $D$ = pile diameter (or width)
- $S$ = center-to-center spacing
- $D_1 = S - D$ = clear spacing between piles
- Soil with properties $c$, $\varphi$, $\gamma$

The theory applies to the portion of the pile **above** the failure surface — this is the zone where soil is actively moving and pushing against the pile.

#### Distributed Lateral Force

The key result is the distributed lateral force $p(z)$ (force per unit length of pile) at depth $z$ below the ground surface. For a $c$-$\varphi$ soil:

$$p(z) = c \cdot A_1(D_1, D, \varphi) + \gamma z \cdot A_2(D_1, D, \varphi)$$

where $A_1$ and $A_2$ are dimensionless coefficients that depend on the pile geometry ratio $D_1/D$ and the friction angle $\varphi$. They involve the passive earth pressure coefficient $N_\varphi = \tan^2(45° + \varphi/2)$ and exponential terms from the plastic flow solution. Key behavior:

- $p(z)$ **increases with depth** (through the $\gamma z$ term) — deeper soil mobilizes more arching pressure
- $p(z)$ **increases as $D_1/D$ decreases** (closer piles = more arching = more force per pile)
- $p(z)$ **increases with $\varphi$** — higher friction angle mobilizes stronger arching
- For $\varphi = 0$ (undrained clay), the equations simplify to a function of $c$ and geometry only

#### Total Force per Pile

Integrate $p(z)$ from the ground surface down to the failure surface:

$$F_{\text{pile}} = \int_0^{z_f} p(z) \, dz$$

where $z_f$ is the depth from the ground surface to the failure surface at the pile location. Since $p(z)$ is linear in $z$ within a homogeneous layer, this integration is straightforward:

$$F_{\text{pile}} = c \cdot A_1 \cdot z_f + \gamma \cdot A_2 \cdot \frac{z_f^2}{2}$$

#### Force per Unit Width

The 2D plane-strain equivalent used in LEM:

$$H = \frac{F_{\text{pile}}}{S}$$

#### Multi-Layer Soils

When the pile passes through multiple material zones above the failure surface (common in practice), the integration is performed piecewise. For each layer $j$ with properties $c_j$, $\varphi_j$, $\gamma_j$ between depths $z_{\text{top},j}$ and $z_{\text{bot},j}$:

$$F_j = c_j \cdot A_{1,j} \cdot (z_{\text{bot},j} - z_{\text{top},j}) + \gamma_j \cdot A_{2,j} \cdot \frac{z_{\text{bot},j}^2 - z_{\text{top},j}^2}{2}$$

$$F_{\text{pile}} = \sum_j F_j$$

Note that $A_{1,j}$ and $A_{2,j}$ must be recomputed for each layer since $\varphi$ may differ.

#### Applicability and Limitations

- **Rigid pile assumption** — the theory assumes piles do not deflect. This is conservative for flexible piles (which mobilize less force).
- **Plastic flow assumption** — gives an upper bound on soil resistance. The actual resistance may be lower if soil hasn't fully mobilized.
- **Spacing ratio**: The theory is applicable for $S/D$ between approximately 2 and 8. Below 2, the piles act more like a wall; above 8, arching becomes negligible.
- **Originally derived for horizontal ground** — applying to sloping ground introduces some approximation, though in practice the theory is routinely used for slopes.
- **No pile structural check** — the computed $H$ is the soil's capacity to push on the pile. The actual resistance is $\min(H_{\text{Ito-Matsui}},\; H_{\text{structural}})$.

#### Phase 3 Implementation in xslope

All the data needed for Ito & Matsui is already available or planned in the input template:

**From the `piles` sheet**: $D_{\text{pile}}$, $S$, pile endpoint coordinates

**From existing `slope_data`**: Ground surface geometry, material zone boundaries, material properties ($c$, $\varphi$, $\gamma$)

**From slice generation**: Failure surface elevation at the pile location ($y_{\text{pile}}$)

The implementation would involve:

1. **`intersect_pile_with_materials()`** — For each pile, intersect the pile line (from ground surface down to failure surface) with material zone boundaries using Shapely. Returns a list of segments with their material properties.

2. **`ito_matsui_coefficients(D1, D, phi)`** — Compute the dimensionless coefficients $A_1$, $A_2$ from the closed-form expressions.

3. **`compute_ito_matsui_force(pile, segments)`** — Integrate $p(z)$ piecewise over each material segment and sum. Divide by $S$ for force per unit width.

4. **Auto-compute logic in `fileio.py` or `slice.py`**: If $H$ is left blank in the `piles` sheet but $D_{\text{pile}}$ and $S$ are provided, automatically compute $H$ using Ito & Matsui. If $H$ is specified, use the user's value (override mode).

### 2.2 User-Specified Force

The simplest approach — the user provides $H$ directly based on external analysis (e.g., from LPILE, GROUP, or other pile analysis software). This is the most flexible option and should be supported regardless of whether we also implement Ito & Matsui.

### 2.3 p-y Curve Analysis

Full soil-pile interaction modeling using p-y curves. This is typically done in separate software and the resulting lateral resistance at the failure surface is extracted and input to the slope stability model. For our purposes, this falls under "user-specified force."

### 2.4 Structural Capacity Check

The pile resistance used in the LEM analysis should not exceed the structural capacity of the pile:
- **Shear capacity**: $V_{\text{cap}}$ of the pile cross-section
- **Moment capacity**: $M_{\text{cap}}$ / (moment arm from fixity point to failure surface)

The controlling value is $\min(\text{soil capacity},\; \text{structural capacity})$.

### 2.5 Recommended Approach for xslope

Three methods are supported, corresponding to implementation phases (see Section 6):

- **User-specified $H$** (Phase 2): The user provides the pile force magnitude and angle directly, based on external analysis. This is the simplest and most flexible option.
- **Ito & Matsui auto-computation** (Phase 3): If $H$ is left blank and $D_{\text{pile}}$/$S$ are provided, $H$ is computed automatically using Ito & Matsui theory. The computation is performed at slice generation time for each trial failure surface, since the integration depth depends on where the surface crosses the pile. See Section 2.1 for theory.
- **FEM beam elements** (Phase 4): The FEM/SSRM approach computes pile resistance naturally from the beam element stiffness ($EI$, $EA$) — no user-specified $H$ needed.


## 3. Input: "piles" Sheet

Both LEM and FEM use the same `piles` sheet in the Excel input template. The sheet structure parallels the "reinforce" sheet but with pile-specific parameters.

### 3.1 Sheet Layout

```
Row 1: Sheet title ("Pile/Pier Data")
Row 2: Headers
Row 3+: Data

     A        B      C      D      E       F      G       H        I       J        K       L
  ──────────────────────────────────────────────────────────────────────────────────────────────
2 │ Label  │  x1  │  y1  │  x2  │  y2   │  H   │ theta │ D │  S    │ E │  I    │  Area
  ──────────────────────────────────────────────────────────────────────────────────────────────
3 │ Pile 1 │ 45.0 │ 32.0 │ 45.0 │ 15.0  │ 50.0 │       │ 0.9    │ 3.0   │        │       │
4 │ Pile 2 │ 55.0 │ 35.0 │ 55.0 │ 12.0  │ 75.0 │  10   │ 0.9    │ 3.0   │        │       │
5 │ Batter │ 40.0 │ 30.0 │ 42.0 │ 10.0  │ 60.0 │       │ 0.9    │ 3.0   │        │       │
6 │        │      │      │      │       │      │       │        │       │        │       │
```

Here is a shot of the actual sheet in Excel:

![alt text](pile_sheet.png)

The sheet name is **"piles"**. The user can input any number of piles, with one pile per row. Each pile is defined by its top and bottom coordinates, the force magnitude $H$, and optional parameters for angle, diameter, spacing, and FEM properties.

### 3.2 Column Definitions

| Column | Field | Required | Description |
|--------|-------|----------|-------------|
| A | Label | No | Name/identifier for the pile |
| B | $x_1$ | Yes | X-coordinate of pile top |
| C | $y_1$ | Yes | Y-coordinate of pile top |
| D | $x_2$ | Yes | X-coordinate of pile tip (bottom) |
| E | $y_2$ | Yes | Y-coordinate of pile tip (bottom) |
| F | $H$ | LEM only | Pile force magnitude per unit width of slope (kN/m or lb/ft) |
| G | $\theta_p$ | No | Force angle from horizontal in degrees (default 0°; positive = upward) |
| H | $D$ | No | Pile diameter |
| I | $S$ | No | Center-to-center spacing |
| J | $E$ | FEM only | Young's modulus of pile material |
| K | $I$ | FEM only | Moment of inertia (computed from $D_{\text{pile}}$ if omitted) |
| L | $Area$ | FEM only | Cross-sectional area (computed from $D_{\text{pile}}$ if omitted) |

### 3.3 Column Usage by Analysis Type

| Column | LEM Use | FEM Use |
|--------|---------|---------|
| $(x_1, y_1)$ | Pile top — defines pile line geometry | Top of 1D element chain |
| $(x_2, y_2)$ | Pile tip — defines pile line geometry | Bottom of 1D element chain |
| $H$ | Pile force magnitude (direct input) | Not used (FEM computes resistance naturally) |
| $\theta_p$ | Force direction (0° = horizontal) | Not used |
| $D$ | Not used in Phase 1 | Compute $I = \pi D^4 / 64$ for circular section |
| $S$ | User divides $H_{\text{single}}/S$ | Scale $EI$ and $EA$ by $1/S$ for per-unit-width |
| $E$, $I$, $Area$ | Not used | Beam element stiffness properties |

If $D$ is provided and $I$, $Area$ are omitted, assume a solid circular section:

$$A = \frac{\pi D^2}{4}, \qquad I = \frac{\pi D^4}{64}$$

### 3.4 Design Decisions

- **Two-endpoint geometry** $(x_1, y_1)$ to $(x_2, y_2)$ supports both vertical and **battered (inclined) piles** from the start. Vertical piles are simply the case where $x_1 = x_2$. This avoids a breaking template change later.
- **Force magnitude + angle** ($H$, $\theta_p$) matches the convention used by Slide2 and SLOPE/W. The user specifies the pile force magnitude and its direction. When $\theta_p = 0°$ (default), the force is purely horizontal. This is more general than assuming horizontal-only and avoids a future template change. The angle is independent of the pile geometry — the user can batter a pile one way and apply force in another direction.
- **$H$ is per unit width** (force/length), consistent with the 2D plane-strain assumption used throughout xslope. If the user has a row of piles at spacing $S$ with individual capacity $H_{\text{single}}$, they input $H = H_{\text{single}} / S$.
- **Pile-failure surface intersection** is computed geometrically during slice generation by intersecting the pile line with the failure surface. This works for both vertical and battered piles.
- The pile must extend below the failure surface to be effective. If the failure surface does not intersect the pile line, the pile provides no resistance.
- Reading stops when column B ($x_1$) is empty (same pattern as reinforcement).
- This format mirrors the `reinforcement_lines` data structure (endpoint pairs), simplifying reuse of mesh extraction and geometric intersection code.

### 3.5 Validation

- Pile endpoints must define a line with nonzero length
- $(x_1, y_1)$ should be above $(x_2, y_2)$: $y_1 > y_2$
- $H > 0$ (required for LEM)
- $\theta_p$ defaults to 0° if omitted

### 3.6 Data Structure

Store as `slope_data['pile_lines']` — a list of dictionaries:
```python
{
    "x1": x1, "y1": y1,       # pile top
    "x2": x2, "y2": y2,       # pile tip
    "H": H,                   # pile force magnitude per unit width (LEM)
    "theta_p": theta_p,       # force angle from horizontal, degrees (default 0)
    "D_pile": D_pile,          # diameter (optional)
    "S": S,                    # spacing (optional)
    "E": E_pile,               # Young's modulus (FEM)
    "I": I,                    # moment of inertia (FEM)
    "area": A,                 # cross-sectional area (FEM)
    "label": label
}
```

### 3.7 Visualization (`plot.py`)

Add a `plot_piles()` function (similar to `plot_reinforcement_lines()`):
- Draw each pile as a thick line from $(x_1, y_1)$ to $(x_2, y_2)$
- Mark the failure surface intersection point
- Annotate with $H$ value
- Use a distinct style (e.g., solid brown/gray rectangle or thick line) to differentiate from reinforcement


## 4. LEM Implementation

### 4.1 Force Components and General Approach

The pile applies a force of magnitude $H$ at angle $\theta_p$ from horizontal at the point where it intersects the failure surface. The force decomposes into:

- $H_h = H\cos\theta_p$: horizontal component
- $H_v = H\sin\theta_p$: vertical component (positive upward)

When $\theta_p = 0°$ (default), $H_h = H$ and $H_v = 0$, recovering the purely horizontal case.

These components can be further resolved relative to the slice base (inclination $\alpha$):

- **Normal to base** (increases effective stress): $H\sin(\alpha - \theta_p) = H_h \sin\alpha - H_v \cos\alpha$
- **Tangential to base** (resists sliding): $H\cos(\alpha - \theta_p) = H_h \cos\alpha + H_v \sin\alpha$

For methods with moment equilibrium about a circle center $(X_o, Y_o)$, the pile force creates a **resisting moment**. The horizontal component contributes $H_h \cdot (Y_o - y_{\text{pile}})$ and the vertical (upward) component contributes $H_v \cdot (x_{\text{pile}} - X_o)$:

$$M_{\text{pile}} = H\cos\theta_p \cdot (Y_o - y_{\text{pile}}) + H\sin\theta_p \cdot (x_{\text{pile}} - X_o)$$

When $\theta_p = 0$, this reduces to $H(Y_o - y_{\text{pile}})$.

**Non-circular failure surfaces**: The pile modifications for Janbu, Corps of Engineers, Lowe-Karafiath, and Spencer have no dependence on circle geometry — they use only per-slice quantities ($\alpha$, $y_b$, force components). These methods work with any failure surface shape, and the pile terms carry over without modification. The only circle-dependent pile terms appear in OMS and Bishop (the moment term above), but those methods inherently require circular surfaces. The `slice.py` pile assignment logic is also geometry-agnostic — it only checks whether the pile line intersects the slice base.

### 4.2 Equation Modifications by Method

#### 4.2.1 OMS (Ordinary Method of Slices) — `solve.py:121`

OMS resolves forces normal to the slice base for $N'$ and uses moment equilibrium about the circle center for FS.

**Modified $N'$** — pile force resolved normal to base:

$$N'_i = W_i \cos\alpha_i + D_i \cos(\alpha_i - \beta_i) - k W_i \sin\alpha_i - T_i \sin\alpha_i - u_i \, \Delta\ell_i + H_i \sin(\alpha_i - \theta_{p,i})$$

The $H_i \sin(\alpha_i - \theta_{p,i})$ term is the pile force resolved normal to the base. When $\theta_p = 0$, this reduces to $H_i \sin\alpha_i$.

**Modified denominator** — pile resisting moment about circle center:

$$\cdots - \sum P_i - \frac{1}{R} \sum \left[ H_i \cos\theta_{p,i} (Y_o - y_{\text{pile},i}) + H_i \sin\theta_{p,i} (x_{\text{pile},i} - X_o) \right] - \frac{1}{R} \sum D_i \sin\beta_i (Y_o - d_{y,i})$$

The horizontal component creates a resisting moment through the vertical moment arm $(Y_o - y_{\text{pile}})$, and the vertical (upward) component creates a resisting moment through the horizontal moment arm $(x_{\text{pile}} - X_o)$. When $\theta_p = 0$, the second term vanishes and this reduces to $-\frac{1}{R} \sum H_i (Y_o - y_{\text{pile},i})$.

#### 4.2.2 Bishop's Simplified Method — `solve.py:256`

Bishop derives $N'$ from **vertical equilibrium** on each slice, then substitutes into **moment equilibrium** about the circle center to solve for $F$ iteratively.

**Modified denominator** (moment equilibrium) — same pile moment term as OMS:

$$\cdots - \sum P_i - \frac{1}{R}\sum \left[ H_i \cos\theta_{p,i} (Y_o - y_{\text{pile},i}) + H_i \sin\theta_{p,i} (x_{\text{pile},i} - X_o) \right] - \frac{1}{R}\sum D_i \sin\beta_i (Y_o - d_{y,i})$$

The pile force is a **known external force** (not a strength parameter), so it is NOT factored by $F$.

**Modified $N'$ from vertical equilibrium** (lines 325–331):

Bishop's vertical equilibrium on the slice is:

$$W_i + D_i \cos\beta_i - H_i \sin\theta_{p,i} - S_i \sin\alpha_i - N_i \cos\alpha_i = 0$$

The vertical component $H_i \sin\theta_{p,i}$ (upward) reduces the net downward force on the slice. Substituting $S = (c \, \Delta\ell + N' \tan\varphi) / F$ and solving for $N'$:

$$N'_i = \frac{W_i + D_i \cos\beta_i - H_i \sin\theta_{p,i} - u_i \, \Delta\ell_i \cos\alpha_i - \dfrac{c_i \, \Delta\ell_i \sin\alpha_i}{F}}{m_{\alpha,i}}$$

where $m_{\alpha,i} = \cos\alpha_i + \dfrac{\sin\alpha_i \tan\varphi_i}{F}$.

**When $\theta_p = 0$**: $H \sin(0) = 0$, so $N'$ is unchanged and the pile affects Bishop only through the moment equilibrium denominator. This is the key simplification for purely horizontal piles — they don't enter vertical equilibrium.

**When $\theta_p \neq 0$**: The vertical component enters $N'$ in the same way that reinforcement $P \sin\alpha$ does — reducing the effective normal force. This also modifies the **numerator** (shear resistance), which uses the same $N'$ terms:

$$\text{shear}_i = \frac{c_i \, \Delta\ell_i \cos\alpha_i + (W_i + D_i \cos\beta_i - H_i \sin\theta_{p,i} - u_i \, \Delta\ell_i \cos\alpha_i) \tan\varphi_i}{m_{\alpha,i}}$$

#### 4.2.3 Janbu's Simplified Method — `solve.py:354`

Janbu uses **force equilibrium parallel to the base**, not moment equilibrium.

**Modified denominator** — pile force resolved parallel to base reduces driving forces:

$$\cdots - \sum P_i - \sum H_i \cos(\alpha_i - \theta_{p,i})$$

All forces in the Janbu denominator are resolved parallel to the slice base. When $\theta_p = 0$, this reduces to $-\sum H_i \cos\alpha_i$.

**Modified $N'$** — full pile force resolved normal to the base:

$$N'_i = W_i \cos\alpha_i - k W_i \sin\alpha_i + D_i \cos(\beta_i - \alpha_i) - T_i \sin\alpha_i + H_i \sin(\alpha_i - \theta_{p,i}) - u_i \, \Delta\ell_i$$

The $H_i \sin(\alpha_i - \theta_{p,i})$ term is the pile force resolved normal to the base. It increases the effective normal force, which increases frictional resistance in the numerator. When $\theta_p = 0$, this reduces to $H_i \sin\alpha_i$.

#### 4.2.4 Force Equilibrium Methods (Corps of Engineers, Lowe-Karafiath) — `solve.py:468`

These methods solve $X$ and $Y$ force equilibrium on each slice.

**Modified $b_0$** (horizontal equilibrium) — subtract horizontal component of pile force (resisting, same sign convention as $-P\cos\alpha$):

$$b_0 = \cdots + k W_i + T_i - H_i \cos\theta_{p,i}$$

**Modified $b_1$** (vertical equilibrium) — subtract upward vertical component:

$$b_1 = \cdots + W_i - H_i \sin\theta_{p,i} - Z_i \sin\theta_i + D_i \cos\beta_i$$

When $\theta_p = 0$: $b_0$ gets $+H_i$ and $b_1$ is unchanged (matching the purely horizontal case).

#### 4.2.5 Spencer's Method — `solve.py:683`

Spencer satisfies both force and moment equilibrium simultaneously. The pile force is treated as an additional external force with known magnitude and direction.

Note: In Spencer's implementation, $P$ = distributed load, $R$ = reinforcement, $V$ = tension crack force.

**Modified $F_h$** — add horizontal component of pile force:

$$F_h = \cdots + R \cos\psi + H_i \cos\theta_{p,i}$$

**Modified $F_v$** — add vertical component (upward):

$$F_v = \cdots + R \sin\psi + H_i \sin\theta_{p,i}$$

**Modified $M_o$** — add pile moment about slice base center, following the same sign convention as reinforcement $R$:

$$M_o = \cdots - H_i \cos\theta_{p,i} \, (y_{\text{pile}} - y_b) + H_i \sin\theta_{p,i} \, (x_{\text{pile}} - x_b)$$

Since Spencer takes moments about the base center of each slice, and the pile force acts at the failure surface (which is the base of the slice containing the pile), the moment arms $(y_{\text{pile}} - y_b)$ and $(x_{\text{pile}} - x_b)$ are typically very small. The moment effect is captured primarily through interslice force propagation.

When $\theta_p = 0$: $F_h$ gets $+H_i$, $F_v$ is unchanged, and $M_o$ gets $-H_i(y_{\text{pile}} - y_b) \approx 0$.

#### 4.2.6 Summary Table

| Method | $N'$ Change | Driving Force/Moment Change | Force Equil. Change |
|--------|-------------|---------------------------|-------------------|
| OMS | $+H\sin(\alpha - \theta_p)$ | $-\frac{1}{R}[H\cos\theta_p(Y_o - y_p) + H\sin\theta_p(x_p - X_o)]$ | N/A |
| Bishop | $-H\sin\theta_p$ in numerator | Same moment term as OMS | N/A |
| Janbu | $+H\sin(\alpha - \theta_p)$ | $-H\cos(\alpha - \theta_p)$ from denominator | N/A |
| Corps/L-K | N/A (implicit) | N/A | $-H\cos\theta_p$ to $b_0$; $-H\sin\theta_p$ to $b_1$ |
| Spencer | N/A (implicit) | N/A | $+H\cos\theta_p$ to $F_h$; $+H\sin\theta_p$ to $F_v$; moment to $M_o$ |

All entries reduce to the previously discussed horizontal-only formulation when $\theta_p = 0$.

### 4.3 Code Changes

#### `fileio.py` — New Sheet Parser

Add a `piles` sheet parser following the same pattern as the `reinforce` sheet parser. Reading stops when column B ($x_1$) is empty.

#### `slice.py` — Pile Force Assignment to Slices

During slice generation, for each pile:

1. Construct a Shapely `LineString` from $(x_1, y_1)$ to $(x_2, y_2)$ for the pile
2. For each slice, construct a `LineString` for the slice base from $(x_l, y_{lb})$ to $(x_r, y_{rb})$
3. Compute the geometric intersection of the pile line with each slice base
4. If intersection exists (a `Point`), assign $H$ to that slice's `'h_pile'` column and store the intersection $y$-coordinate in `'y_pile'`

This approach works identically for vertical and battered piles, and mirrors how reinforcement intersections are already computed in `slice.py` (lines 908–932).

**New columns in slice_df**:
- `'h_pile'`: Pile force magnitude per unit width (0 for slices without piles)
- `'theta_p'`: Pile force angle from horizontal in radians (0 for slices without piles)
- `'y_pile'`: $Y$-coordinate of pile-failure surface intersection (0 for slices without piles)
- `'x_pile'`: $X$-coordinate of pile-failure surface intersection (0 for slices without piles)

If multiple piles intersect the same slice base, sum their $H$ values (though this should be unusual in practice).

#### `solve.py` — Method Modifications

Each method needs to extract and use the new `h_pile` column:

```python
H_pile  = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
theta_p = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)
y_pile  = slice_df['y_pile'].values  if 'y_pile'  in slice_df.columns else np.zeros(n)
x_pile  = slice_df['x_pile'].values  if 'x_pile'  in slice_df.columns else np.zeros(n)
```

This ensures backward compatibility — if no pile data exists, the terms are zero and equations reduce to the current form.

**Method-specific changes**:

- **`oms()`**: Add $H \sin(\alpha - \theta_p)$ to $N'$; add pile moment term (both components) to denominator
- **`bishop()`**: Add $-H \sin\theta_p$ to $N'$ numerator; add pile moment term to denominator
- **`janbu()`**: Add $H \sin(\alpha - \theta_p)$ to $N'$; subtract $\sum H \cos(\alpha - \theta_p)$ from denominator
- **`force_equilibrium()`**: Subtract $H \cos\theta_p$ from $b_0$; subtract $H \sin\theta_p$ from $b_1$
- **`spencer()`**: Add $H \cos\theta_p$ to $F_h$; add $H \sin\theta_p$ to $F_v$; add moment terms to $M_o$


## 5. FEM / SSRM Implementation

### 5.1 Comparison with Existing Reinforcement Elements

The current FEM module models flexible reinforcement as **2-node truss elements** (axial only, tension only):

- Extracted from 2D mesh edges that lie along reinforcement lines (`extract_1d_elements_from_2d_edges()`)
- Shared nodes with the 2D soil mesh — natural displacement coupling, no constraint equations
- Element stiffness: $k = EA/L$, assembled as a $4 \times 4$ global matrix using direction cosines
- **Tension only** — compression is zeroed out via body-force corrections during viscoplastic iteration
- **Not reduced during SSRM** — only soil $c$ and $\varphi$ are reduced; reinforcement $E$, $A$, $T_{\max}$ remain constant
- Failure: when $T > T_{\text{allow}}$, element drops to residual capacity $T_{\text{res}}$

Piles reuse much of this infrastructure but differ in key ways:

| Property | Reinforcement (truss) | Pile (beam) |
|----------|----------------------|-------------|
| Axial stiffness | $EA/L$ | $EA/L$ |
| Lateral stiffness | None | $12EI/L^3$ (from bending) |
| Compression | Not allowed (zeroed) | Allowed |
| Orientation | Inclined (along reinforcement) | Vertical or battered (any angle) |
| Failure mode | Tension rupture / pullout | Shear / bending capacity |
| DOFs per node | 2 (translational only) | 2 (translational only, rotational condensed) |

### 5.2 Beam Element Formulation

The current 2D soil mesh has 2 DOFs per node ($u_x$, $u_y$). A full Euler-Bernoulli beam element has 3 DOFs per node ($u_x$, $u_y$, $\theta$). To avoid adding rotational DOFs to the global system, we use **static condensation** to express the beam stiffness in terms of translational DOFs only.

**Full beam element stiffness** (6 DOFs: $u_1, v_1, \theta_1, u_2, v_2, \theta_2$, beam aligned along local $x$-axis):

$$\mathbf{K}_{\text{full}} = \begin{bmatrix} \frac{EA}{L} & 0 & 0 & -\frac{EA}{L} & 0 & 0 \\ 0 & \frac{12EI}{L^3} & \frac{6EI}{L^2} & 0 & -\frac{12EI}{L^3} & \frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{4EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{2EI}{L} \\ -\frac{EA}{L} & 0 & 0 & \frac{EA}{L} & 0 & 0 \\ 0 & -\frac{12EI}{L^3} & -\frac{6EI}{L^2} & 0 & \frac{12EI}{L^3} & -\frac{6EI}{L^2} \\ 0 & \frac{6EI}{L^2} & \frac{2EI}{L} & 0 & -\frac{6EI}{L^2} & \frac{4EI}{L} \end{bmatrix}$$

**Static condensation** partitions into translational ($tt$) and rotational ($rr$) blocks:

$$\mathbf{K}_{\text{condensed}} = \mathbf{K}_{tt} - \mathbf{K}_{tr} \, \mathbf{K}_{rr}^{-1} \, \mathbf{K}_{rt}$$

For a beam element aligned along the local $x$-axis, this yields a $4 \times 4$ condensed stiffness matrix in terms of $[u_1, v_1, u_2, v_2]$:

$$\mathbf{K}_{\text{condensed}}^{\text{local}} = \begin{bmatrix} \frac{EA}{L} & 0 & -\frac{EA}{L} & 0 \\ 0 & \frac{3EI}{L^3} & 0 & -\frac{3EI}{L^3} \\ -\frac{EA}{L} & 0 & \frac{EA}{L} & 0 \\ 0 & -\frac{3EI}{L^3} & 0 & \frac{3EI}{L^3} \end{bmatrix}$$

The lateral stiffness reduces from $12EI/L^3$ to $3EI/L^3$ after condensation — this corresponds to a beam with rotations free at both ends (pin-pin for rotation), which is conservative. In practice, with multiple short elements along the pile connected at shared nodes, the chain of elements recovers the correct bending behavior because intermediate node rotations are implicitly constrained by the series of elements.

**Transformation to global coordinates** uses the same direction cosine approach as the existing truss elements:

$$\mathbf{K}_{\text{global}} = \mathbf{T}^T \, \mathbf{K}_{\text{condensed}}^{\text{local}} \, \mathbf{T}$$

For a vertical pile ($\cos\theta = 0$, $\sin\theta = 1$ if directed upward), the axial stiffness $EA/L$ acts in the $y$-direction and the lateral stiffness $3EI/L^3$ acts in the $x$-direction — which is exactly the behavior we want (lateral resistance to slope movement).

### 5.3 Mesh Generation

The existing `extract_1d_elements_from_2d_edges()` function can be reused directly. It identifies 2D mesh edges that are collinear with a given line and extracts them as 1D elements with shared nodes. Since pile lines are defined as two endpoints $(x_1, y_1)$ to $(x_2, y_2)$ — the same format as reinforcement lines — the function handles both vertical and battered piles without modification.

- The mesh generator must be constrained so that element edges align along the pile line (same as how reinforcement lines are handled)
- The extracted 1D elements become beam elements instead of truss elements
- The distinction between truss vs. beam behavior is made during `build_fem_data()` based on element type

### 5.4 Changes to `build_fem_data()`

Add a parallel path for pile elements, similar to the existing reinforcement loop (`fem.py:258–341`):

1. **Extract 1D elements** along pile lines (reuse `extract_1d_elements_from_2d_edges()`)
2. **Compute element properties**:
   - Axial stiffness: $k_a = EA/L$
   - Lateral stiffness: $k_l = 3EI/L^3$
   - Direction cosines: $\cos\theta$, $\sin\theta$ (typically $0, \pm 1$ for vertical piles)
3. **Build condensed $4 \times 4$ stiffness matrix** per element (Section 5.2)
4. **Store pile-specific data**:
   - `k_axial_by_pile_elem`, `k_lateral_by_pile_elem`
   - `cos_theta_pile`, `sin_theta_pile`
   - `dof_indices_pile`
   - `K_global_pile_elems` (list of $4 \times 4$ matrices)
   - Failure capacities: shear capacity $V_n$, moment capacity $M_p$ (optional)

### 5.5 Changes to Viscoplastic Iteration

The body-force correction loop (`fem.py:893–948`) needs a parallel section for pile elements:

1. **Compute axial and lateral forces** from nodal displacements:
   - Axial: $T = k_a \cdot \delta_{\text{axial}}$ (projection along element axis)
   - Lateral (shear): $V = k_l \cdot \delta_{\text{lateral}}$ (projection perpendicular to element axis)

2. **Key difference from reinforcement**: Piles carry **both tension and compression** — do NOT zero out compression forces.

3. **Failure checks**:
   - **Shear capacity**: If $|V| > V_n$, cap and apply correction forces
   - **Axial capacity**: If $|T| > T_{\max}$ (structural), cap and apply correction
   - Failure criteria could be simplified to just checking if nodal forces exceed the user-specified $H$ value, or left unchecked for Phase 1 (let the SSRM find the natural failure mode)

4. **Correction forces**: Same pattern as reinforcement — compute the overshoot and apply as nodal load corrections, but now in both axial and lateral directions.

### 5.6 SSRM Treatment

- Pile material properties ($E$, $I$, $A$, yield strength) are **NOT reduced** during strength reduction — same as reinforcement
- The global stiffness matrix $\mathbf{K}$ includes pile contributions and is factored once with `splu`, then reused
- As soil strength is reduced, the soil deforms more around the pile; the pile's stiffness naturally resists this deformation
- Convergence failure (i.e., slope collapse) occurs when the soil is too weak to transfer load to/from the piles


## 6. Implementation Phases

### Phase 1: Documentation and Equations

Update the documentation first so the theory is reviewed and agreed upon before writing code. This also produces the free body diagrams and equation derivations that guide the implementation.

1. **LEM overview** (`docs/lem/overview.md`): Add "Piles and Structural Elements" subsection with annotated free body diagram showing $H$ at angle $\theta_p$ on a slice
2. **Individual method pages**: Integrate pile terms into existing equation derivations for each method:
   - `docs/lem/oms.md` — $H\sin(\alpha - \theta_p)$ in $N'$, pile moment in denominator
   - `docs/lem/bishop.md` — $-H\sin\theta_p$ in $N'$ numerator, pile moment in denominator
   - `docs/lem/janbu.md` — $H\sin(\alpha - \theta_p)$ in $N'$, $-H\cos\theta_p$ in denominator
   - `docs/lem/force_eq.md` — $-H\cos\theta_p$ in $b_0$, $-H\sin\theta_p$ in $b_1$
   - `docs/lem/spencer.md` — pile terms in $F_h$, $F_v$, $M_o$
3. **Input template** (`docs/usage/input_template.md`): Add `piles` sheet section with column definitions and examples

### Phase 2: LEM with User-Specified $H$

Core LEM implementation — user provides $H$ and $\theta_p$ directly.

1. **`fileio.py`**: Add `piles` sheet parser. Read $(x_1, y_1)$, $(x_2, y_2)$, $H$, $\theta_p$, and optional columns. Store in `slope_data['pile_lines']`. Validate inputs.
2. **`slice.py`**: Add pile-to-slice assignment:
   - Construct Shapely `LineString` for each pile
   - Intersect with each slice base to find the crossing point
   - Add `h_pile`, `theta_p`, `x_pile`, `y_pile` columns to `slice_df`
3. **`solve.py`**: Add pile terms to each method (start with Janbu as simplest to verify, then OMS/Bishop, force equilibrium, Spencer). Use backward-compatible pattern: `H_pile = slice_df['h_pile'].values if 'h_pile' in slice_df.columns else np.zeros(n)`
4. **`plot.py`**: Add `plot_piles()` — draw pile lines, mark failure surface intersection, annotate with $H$
5. **Input template**: Create `piles` sheet in the Excel template
6. **Testing**:
   - Create test input file with vertical piles at toe of an existing sample problem
   - Verify FS increases with piles vs. without
   - Add to `run_tests.py` regression suite with `<!-- test: -->` tags
   - Test with $\theta_p = 0$ and $\theta_p \neq 0$
7. **Sample problem**: Add pile example to `docs/lem/samples.md`

### Phase 3: Ito & Matsui Auto-Computation (Optional)

Auto-compute $H$ from pile geometry and soil properties when $H$ is left blank.

**Computation characteristics**: Ito & Matsui is a closed-form equation — no iteration, no convergence issues. It evaluates coefficients from $\tan\varphi$ and $N_\varphi = \tan^2(45° + \varphi/2)$ plus an exponential, then integrates a linear function over each soil layer. Computation is essentially instantaneous (microseconds per pile). Well-behaved for all practical friction angles; $\varphi = 0$ (pure clay) has a separate simpler equation.

**When to compute**: At **slice generation time** (in `slice.py`), not at file load time. The reason: Ito & Matsui integrates the lateral pressure from the ground surface down to the failure surface. The failure surface location changes with each trial circle during the search. Computing $H$ per trial surface is physically correct — a deeper failure surface means more soil above it pushing on the pile, giving a higher $H$. Since the computation is instantaneous, there is no performance penalty for recomputing on every trial. This is the approach used by commercial codes that implement Ito & Matsui (e.g., ReSSA by ADAMA Engineering).

**Implementation**:

1. **`ito_matsui.py`** (new module):
   - `ito_matsui_coefficients(D1, D, phi)` — closed-form $A_1$, $A_2$ coefficients
   - `intersect_pile_with_materials(pile, slope_data, y_failure)` — find soil layers along the pile between the ground surface and failure surface intersection
   - `compute_ito_matsui_force(pile, segments, S)` — piecewise integration over layers, divide by spacing $S$
2. **`slice.py`**: During pile-to-slice assignment, if `H` is blank and `D_pile`/`S` are provided:
   - Determine the failure surface $y$-coordinate at the pile intersection (already computed in Phase 2)
   - Call `compute_ito_matsui_force()` with the pile geometry, soil layers, and failure surface depth
   - Assign the computed $H$ to the slice's `h_pile` column
   - If `H` is provided by the user, use that value instead (user override)
3. **Validation**: Compare auto-computed $H$ against hand calculations and published examples
4. **Documentation**: Add Ito & Matsui theory section, worked example showing auto-computed vs. user-specified $H$

### Phase 4: FEM Beam Elements

Add pile support to the FEM/SSRM module.

1. **`mesh.py`**: Extend mesh generation to constrain edges along pile lines (same pattern as reinforcement lines). Extract 1D beam elements via `extract_1d_elements_from_2d_edges()`.
2. **`fem.py` — `build_fem_data()`**: Add beam element path parallel to existing reinforcement:
   - Compute axial stiffness $k_a = EA/L$ and lateral stiffness $k_l = 3EI/L^3$
   - Build condensed $4 \times 4$ stiffness matrix per element
   - Store pile-specific arrays (`k_axial_by_pile_elem`, `k_lateral_by_pile_elem`, etc.)
3. **`fem.py` — viscoplastic iteration**: Add pile body-force corrections:
   - Compute axial and lateral forces from nodal displacements
   - Allow both tension and compression (unlike reinforcement)
   - Optional shear/axial capacity checks
4. **`fem.py` — SSRM**: Ensure pile properties are NOT reduced during strength reduction
5. **`plot_fem.py`**: Visualize pile elements in deformed mesh and force plots
6. **Testing**: Run FEM on same pile problem as LEM, compare FS values
7. **Sample problem**: Add FEM pile example to `docs/fem/samples.md`


## 7. Design Decisions

1. **No separate moment input.** LEM uses force magnitude $H$ and angle $\theta_p$ only — no moment capacity column. This matches Slide2 and SLOPE/W. PLAXIS handles moments naturally through FEM beam elements, which we will also do in Phase 2. Adding a separate $M$ term to LEM introduces ambiguity without meaningful improvement.

2. **Force magnitude + angle from the start.** The $\theta_p$ column (defaulting to $0° =$ horizontal) matches the Slide2 and SLOPE/W convention. All LEM equations are generalized with both $H\cos\theta_p$ and $H\sin\theta_p$ components. This avoids a future template-breaking change if angled forces are needed, and costs almost nothing to implement since the $\theta_p = 0$ case reduces to the simpler horizontal-only formulation.

3. **Multiple failure surfaces** — automatic. Piles that don't intersect a given trial surface contribute $H_{\text{pile}} = 0$ for those slices. No special handling needed in the search algorithm.

4. **Pile below failure surface** — automatic. If the failure surface does not intersect the pile line, Shapely intersection returns empty and the pile contributes nothing. Handled naturally during slice generation.


## 8. Typical Parameter Values

The following tables provide typical values for pile/pier parameters. These are intended as guidance for users filling in the `piles` sheet and will be included in the documentation.

### 8.1 Pile Types and Typical Dimensions

| Pile Type | Typical Diameter/Width | Typical Length | Notes |
|-----------|----------------------|----------------|-------|
| Steel H-pile | 200–400 mm (8–14 in) | 6–30 m (20–100 ft) | Wide-flange sections; high strength-to-weight ratio |
| Steel pipe pile | 300–900 mm (12–36 in) | 6–40 m (20–130 ft) | Can be open-ended or closed-ended; often concrete-filled |
| Precast concrete pile | 250–600 mm (10–24 in) | 6–25 m (20–80 ft) | Square or octagonal cross-section |
| Drilled shaft (caisson) | 450–3000 mm (18–120 in) | 3–60 m (10–200 ft) | Cast-in-place; larger diameters common for slope stabilization |
| Concrete pier | 600–1500 mm (24–60 in) | 3–15 m (10–50 ft) | Often used for slope stabilization; rectangular or circular |
| Micropile | 150–300 mm (6–12 in) | 6–30 m (20–100 ft) | Drilled and grouted; used in tight access or existing structures |
| Timber pile | 200–400 mm (8–16 in) | 6–20 m (20–65 ft) | Tapered; limited to lighter loads |

### 8.2 Typical Spacing

| Application | Typical $S/D$ Ratio | Typical Spacing $S$ | Notes |
|-------------|---------------------|---------------------|-------|
| Slope stabilization | 3–6 | 1.5–6 m (5–20 ft) | Closer spacing = more arching between piles |
| Retaining structures | 2–4 | 1–4 m (3–12 ft) | Often soldier piles with lagging |
| Ito & Matsui applicability | 2–8 | — | Theory assumes plastic flow between piles |

### 8.3 Material Properties (FEM)

| Material | $E$ (kPa) | $E$ (psf) | $\nu$ |
|----------|-----------|-----------|-------|
| Structural steel | 2.0 × 10⁸ | 4.18 × 10⁹ | 0.3 |
| Reinforced concrete ($f'_c$ = 4000 psi) | 2.5 × 10⁷ | 5.2 × 10⁸ | 0.2 |
| Reinforced concrete ($f'_c$ = 5000 psi) | 2.8 × 10⁷ | 5.8 × 10⁸ | 0.2 |
| Reinforced concrete ($f'_c$ = 6000 psi) | 3.0 × 10⁷ | 6.3 × 10⁸ | 0.2 |
| Timber (Douglas Fir) | 1.2 × 10⁷ | 2.5 × 10⁸ | 0.3 |

Concrete modulus computed as $E_c = 57{,}000 \sqrt{f'_c}$ (psi) per ACI 318. Values should be entered in the same unit system used for the rest of the input (kPa or psf).

### 8.4 Typical Section Properties

| Section | $A$ (Area) | $I$ (Moment of Inertia) |
|---------|-----------|------------------------|
| Solid circular, $D$ = 0.6 m (24 in) | 0.283 m² (452 in²) | 6.36 × 10⁻³ m⁴ (16,286 in⁴) |
| Solid circular, $D$ = 0.9 m (36 in) | 0.636 m² (1,018 in²) | 3.22 × 10⁻² m⁴ (82,448 in⁴) |
| Solid circular, $D$ = 1.2 m (48 in) | 1.131 m² (1,810 in²) | 1.02 × 10⁻¹ m⁴ (260,576 in⁴) |
| HP 14×117 (steel H-pile) | 0.022 m² (34.4 in²) | 4.43 × 10⁻⁴ m⁴ (1,063 in⁴) |
| HP 12×84 (steel H-pile) | 0.016 m² (24.6 in²) | 2.18 × 10⁻⁴ m⁴ (524 in⁴) |
| Pipe pile, $D$ = 0.6 m, $t$ = 12 mm | 0.022 m² (34.6 in²) | 9.7 × 10⁻⁴ m⁴ (2,330 in⁴) |

For solid circular sections: $A = \pi D^2 / 4$, $I = \pi D^4 / 64$.

### 8.5 Typical Lateral Resistance $H$

Lateral resistance depends heavily on soil conditions, pile geometry, and embedment. The following ranges are approximate for preliminary estimates only.

| Soil Type | Pile Type | Typical $H$ per pile | Notes |
|-----------|-----------|---------------------|-------|
| Stiff clay ($c$ = 50–100 kPa) | Drilled shaft, $D$ = 0.9 m, $S$ = 3 m | 100–400 kN/m (7–27 kip/ft) | Per unit width = $H_{\text{pile}} / S$ |
| Medium dense sand ($\varphi$ = 30–35°) | Steel H-pile, $S$ = 2 m | 50–200 kN/m (3–14 kip/ft) | Depends on depth above failure surface |
| Soft clay ($c$ = 15–30 kPa) | Concrete pier, $D$ = 1.2 m, $S$ = 3 m | 30–100 kN/m (2–7 kip/ft) | Lower bound; may govern over structural capacity |
| Weathered rock | Drilled shaft, $D$ = 0.9 m, $S$ = 3 m | 200–800 kN/m (14–55 kip/ft) | High capacity but expensive installation |

These values are for guidance only. Actual $H$ should be determined from Ito & Matsui theory, p-y analysis, or structural analysis of the pile.


## 9. Documentation Updates

The following documentation pages need to be created or updated when pile support is implemented. The existing docs follow a consistent pattern: MathJax equations, annotated figures, Excel template screenshots, and worked sample problems.

### 8.1 Input Template (`docs/usage/input_template.md`)

This is the primary reference for all Excel worksheets. Currently documents 10 sheets (main, plot, mat, profile, piezo, circles, non-circ, dloads, reinforce, seep bc).

**Add a new "piles" sheet section** covering:
- Sheet layout diagram (matching the format used for the reinforce sheet)
- Column-by-column descriptions: $(x_1, y_1)$, $(x_2, y_2)$, $H$, $\theta_p$, $D_{\text{pile}}$, $S$, and FEM columns
- Which columns are required vs. optional, and which apply to LEM vs. FEM
- Default values ($\theta_p = 0°$, FEM columns computed from $D_{\text{pile}}$ if omitted)
- Guidance on computing $H$: user-specified from external analysis, or auto-computed via Ito & Matsui (Phase 3)
- Note on per-unit-width convention: $H = H_{\text{single}} / S$
- Example entries for vertical piles, battered piles

### 8.2 LEM Overview (`docs/lem/overview.md`)

The overview page introduces the method of slices and covers all loading conditions (distributed loads, seismic, reinforcement, tension cracks). Currently has a section on reinforcement describing the $P$ force parallel to the slice base.

**Add a "Piles and Structural Elements" subsection** covering:
- How piles differ from flexible reinforcement (shear/bending vs. tension)
- Force $H$ at angle $\theta_p$ applied at the failure surface intersection
- Resolution into horizontal and vertical components ($H\cos\theta_p$, $H\sin\theta_p$)
- Resolution relative to the slice base (normal and tangential)
- Annotated **free body diagram** of a slice with pile force — this is the key figure, showing $H$ at angle $\theta_p$ acting on the slice base alongside existing forces ($W$, $N'$, $S$, $P$, $u\Delta\ell$, interslice forces)
- Note on non-circular surfaces (works for Janbu, force equilibrium, Spencer)

### 8.3 Individual LEM Method Pages

Each method page (`oms.md`, `bishop.md`, `janbu.md`, `force_eq.md`, `spencer.md`) currently derives the FS equation and shows how reinforcement $P$ enters the formulation.

**For each method, add the pile terms to the existing equations:**

- **`docs/lem/oms.md`**: Add $H\sin(\alpha - \theta_p)$ to $N'$ equation; add pile moment term to denominator; show reduction to $\theta_p = 0$ case
- **`docs/lem/bishop.md`**: Add $-H\sin\theta_p$ to $N'$ vertical equilibrium; add pile moment to denominator; explain why $N'$ is unchanged when $\theta_p = 0$
- **`docs/lem/janbu.md`**: Add $H\sin(\alpha - \theta_p)$ to $N'$; subtract $H\cos(\alpha - \theta_p)$ from denominator
- **`docs/lem/force_eq.md`**: Subtract $H\cos\theta_p$ from $b_0$; subtract $H\sin\theta_p$ from $b_1$
- **`docs/lem/spencer.md`**: Add pile terms to $F_h$, $F_v$, and $M_o$

For each page, the updates should be integrated into the existing equation derivations (not a separate section). The pile terms should appear naturally alongside the existing reinforcement, distributed load, and seismic terms.

### 8.4 FEM Documentation

**`docs/fem/overview.md`** — Add a section on pile beam elements:
- How piles are modeled as 1D beam elements with shared nodes
- Static condensation of rotational DOFs (the condensed $4 \times 4$ stiffness matrix)
- Comparison table: truss (reinforcement) vs. beam (piles)
- Key behavioral differences: compression allowed, lateral stiffness from $EI$
- SSRM treatment: pile properties not reduced during strength reduction

**`docs/fem/mesh.md`** — Update mesh generation documentation:
- How pile lines are treated as mesh constraints (same as reinforcement lines)
- Extraction of 1D beam elements from 2D mesh edges along pile lines

### 8.5 Sample Problems

**New LEM sample problem** (`docs/lem/samples.md`):
- Simple slope with a row of vertical piles at the toe
- Show the input template with the `piles` sheet filled in
- Run all methods and compare FS with and without piles
- Verify that FS increases appropriately
- Good candidate: take an existing sample problem and add piles to it, so the reader can see the incremental effect

**New FEM sample problem** (`docs/fem/samples.md`):
- Same slope geometry as the LEM sample, but with FEM beam elements
- Compare FEM FS with LEM FS (they should be in reasonable agreement)
- Show the deformed mesh with pile elements highlighted
- Show pile force distribution along the pile length

**Battered pile example** (could be in either LEM or FEM samples):
- Demonstrate a battered pile with $\theta_p \neq 0$
- Show how the force angle affects FS compared to a vertical pile with the same $H$

### 8.6 API Reference

Update auto-generated API docs for modified modules:
- **`docs/api/fileio.md`** — new `piles` sheet parser
- **`docs/api/slice.md`** — new `h_pile`, `theta_p`, `y_pile`, `x_pile` columns in slice_df
- **`docs/api/solve.md`** — updated method signatures (no changes, but docstrings updated)
- **`docs/api/plot.md`** — new `plot_piles()` function
- **`docs/api/fem.md`** — new beam element data structures and functions

### 8.7 mkdocs.yml Navigation

No new top-level sections needed. The pile documentation fits naturally into existing sections:
- Input template → add piles sheet description
- LEM method pages → integrate pile terms into equations
- FEM pages → add beam element section
- Sample problems → add pile examples to existing sample pages

---

## 9. Future: Pile Structural Capacity Checks

### Motivation

Neither LEM nor FEM currently enforces structural capacity limits on piles. In LEM, the
Ito & Matsui force H can exceed the pile's shear or moment capacity — the documentation
recommends `H_design = min(H_soil, V_cap, M_cap/L_m)` but this is not enforced in code.
In FEM, beam elements are linearly elastic with infinite strength — the SSRM finds the FS
at which soil fails around the pile, which is unconservative if the pile would fail structurally
first. Commercial software (PLAXIS, RS2, FLAC, Slide2) supports structural capacity limits.

### Input Changes

Add two optional columns to the `piles` sheet:

| Column | Field | Description |
|--------|-------|-------------|
| N | V_cap | Shear capacity per pile (force units). Blank = no limit. |
| O | M_cap | Moment capacity per pile (force × length units). Blank = no limit. |

Both are per-pile values. If blank, current behavior is preserved (no structural limit).
These columns are shared by both LEM and FEM analyses.

### LEM Implementation

V_cap and M_cap are properties of a single pile, so the comparison must be against
F_pile (total force per pile), not H (force per unit width = F_pile / S).

The capacity check applies regardless of how the pile force was obtained — whether
auto-computed by Ito & Matsui or specified directly by the user.

When V_cap and/or M_cap are provided, cap F_pile before converting to H:

1. Obtain F_pile:
   - If H is auto-computed by Ito & Matsui: F_pile is already computed
   - If H is user-specified: F_pile = H × S
2. If V_cap is provided: `F_pile = min(F_pile, V_cap)`
3. If M_cap is provided: estimate the moment arm L_m from the pressure distribution
   and cap `F_pile = min(F_pile, M_cap / L_m)`. For the Ito & Matsui linear pressure
   distribution (p = c*A1 + gamma*z*A2), L_m can be computed from the centroid of the
   pressure diagram. For user-specified H where the distribution is unknown, a
   conservative default (e.g., L_m = h/3) could be used.
4. Convert to per-unit-width: `H = F_pile / S`
5. Use the capped H in the slice equilibrium equations

The Ito & Matsui summary should report both the soil force and the capped force when
capacity governs.

### FEM Implementation

The existing viscoplastic correction loop for soil elements provides the pattern:

1. After each SSRM solve iteration, compute beam element forces (already done for summary)
2. For each beam element, check lateral force against V_cap and moment against M_cap
3. If exceeded, compute a correction force (excess beyond capacity) as equivalent nodal loads
4. Add corrections to the RHS vector and re-iterate (same as soil viscoplastic corrections)
5. The beam carries a constant force at the capacity value (elastic-perfectly-plastic)

### Scope Estimate

**LEM**:
- Input handling: ~10 lines (read V_cap, M_cap from piles sheet)
- Force capping in slice.py: ~15 lines (min check after Ito & Matsui computation)
- Summary output: ~5 lines (report when capacity governs)

**FEM**:
- Capacity check in SSRM loop: ~30 lines (compare forces, compute excess)
- Force correction / nodal load computation: ~30 lines (convert excess to equivalent nodal forces)
- Summary output updates: ~10 lines (report yielded elements)
- Testing: need a problem where structural failure governs (small pile in strong soil)

### What This Enables

- Realistic modeling of pile structural failure (plastic hinge formation)
- Proper FS when pile structural capacity governs over soil failure
- Comparison with LEM's `H_design = min(H_ito_matsui, H_structural)` approach
- More meaningful LEM-FEM comparisons since both methods would respect structural limits
