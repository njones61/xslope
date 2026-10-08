# Soil Reinforcement in LEM Slope Stability

Soil and rock carry compression and shear well but little or no tension. Reinforcement supplies the tension:
members placed across the zone where a slip surface would form and anchored in the stable ground beyond it.
Geosynthetic layers (geotextiles and geogrids) are built into a fill as it is placed, so a slope can stand steeper
than the fill alone would; soil nails are drilled and grouted into a cut as it is excavated; tiebacks hold a wall
with a grouted bond length deep behind the slip surface, and end-anchored bars with a plate or deadman at each end.
In each case the part of the
member beyond the slip surface grips the stable ground, and the member pulls back on the sliding mass.

![Three kinds of reinforcement: geosynthetic layers in a fill, soil nails in a cut, tiebacks behind a wall](../fem/images/reinf_types.png){width=1000}

From left to right: geosynthetic layers in a reinforced fill, soil nails in a cut, and tiebacks behind a wall,
each with a slip surface they cross.

In a limit equilibrium analysis each reinforcement line is a straight line defined by its end points, and wherever a
trial slip surface crosses a line, a tensile force is applied to the sliding mass at the crossing point. Three
questions determine that force, and they are independent of one another:

1. **How large is the force?** — governed by the *capacity envelope*: the tensile strength of the element, the
   frictional pullout development from each end, and any end anchorage (plates, connections, anchors).
2. **In what direction does it act?** — governed by the **Dir** setting: tangent to the slip surface (flexible
   reinforcement) or along the reinforcement's own axis (rigid supports).
3. **Is it factored by the safety factor?** — governed by the **Appl** setting: active (a known allowable force,
   not divided by $F$) or passive (an ultimate capacity that mobilizes with the soil, divided by $F$).

This decomposition follows the convention used by Slide2 and other commercial programs, which allows xslope
results to be compared directly against them. The **Type** column in the input template is a *preset* over these
settings — selecting a support type fills Dir and Appl with the appropriate defaults — not a separate mechanism.

What to enter for a particular support (a geosynthetic layer, a soil nail, a tieback or an end-anchored bar), and
how to connect it to a wall or make it a joint, is set out column by column on
[Modeling Reinforcement](../usage/modeling_reinforcement.md).

The finite element treatment, in which each reinforcement line is a row of bar elements in the mesh and carries the
tension its stretch produces, is on the [FEM reinforcement](../fem/reinforcement.md) page.

## Capacity Envelope

### Force magnitude at the crossing point

The tensile force available at any point along a reinforcement line is limited by three mechanisms, and the
available force is the smallest of them:

>$T(x) = \min\left(T_{max},\;\; T_{end1} + T_{max}\dfrac{d_1}{L_{p1}},\;\; T_{end2} + T_{max}\dfrac{d_2}{L_{p2}}\right)$

where:

- $T_{max}$ = tensile capacity of the element (rupture limit)
- $d_1$, $d_2$ = distances from the point to end 1 and end 2 of the line
- $L_{p1}$, $L_{p2}$ = pullout lengths at each end — the distance over which interface friction develops the full
  tensile capacity
- $T_{end1}$, $T_{end2}$ = anchorage capacity at each end: a bearing plate, a facing connection, or an end anchor
  (default 0)

Special cases:

- **$T_{end} = 0$, $L_p > 0$** (the classical friction-only taper): tension is zero at the free end and develops
  linearly over the pullout length. This is the correct model for geosynthetics and for the embedded end of nails.
- **$L_p = 0$**: the end is fully anchored — the full capacity is available immediately at the end.
- **$T_{end} > 0$**: the end starts at the anchorage capacity and frictional development adds to it. This models
  a nail with a bearing plate at the wall face, a geosynthetic connected to facing panels, or a bar anchored at
  both ends. (If $T_{end} \geq T_{max}$, the end is effectively fully anchored — the tendon governs.)
- **Line shorter than $L_{p1} + L_{p2}$** with no anchorage: the envelopes from the two ends intersect below
  $T_{max}$ and only partial tension is mobilized.

![reinf_envelope.png](images/reinf_envelope.png)

The envelope for each of the four end conditions. The force available where a trial surface crosses the line is
the envelope value at the crossing point, so a surface that clips a line near a free end mobilizes only a
fraction of $T_{max}$.

### Pullout from the effective overburden

$L_p$ states the bond as a development length: the capacity grows at a constant rate $T_{max}/L_p$ no matter how
deep the reinforcement is buried. Interface friction does not work that way — it grows with the normal stress
pressing the soil onto the reinforcement — so the **Adhesion** and **Delta** columns state the interface strength
instead and let the resistance follow the depth of burial. Per unit length of a planar reinforcement, with soil
bearing on both faces:

>$r(s) = 2\left(a + \sigma'_v(s)\tan\delta\right)$

where $a$ is the soil–reinforcement adhesion (stress units), $\delta$ the interface friction angle (degrees), and
$\sigma'_v(s)$ the **effective** vertical stress at the point $s$ along the line: the weight of the soil column
standing above that point — every material zone it crosses at that material's unit weight, saturated below the
water table where the material declares a $\gamma_{sat}$ — less the pore pressure the model declares there
(piezometric line, $r_u$, or seepage field, in the same way as for a slice base).

The envelope is then the same three-way minimum, with the ramps integrated rather than assumed linear:

>$T(s) = \min\left(T_{max},\;\; T_{end1} + \displaystyle\int_0^s r,\;\; T_{end2} + \int_s^L r\right)$

A constant $r$ recovers the straight ramps above, so the two laws are one formula. Both columns filled selects
this law and $L_{p1}$/$L_{p2}$ are then not read; both blank is the default and selects the development-length
law. Filling one column and leaving the other blank is an input error. LEM and FEM use the same envelope under
either law.

**FHWA pullout capacity.** The FHWA form $F^{*}\alpha\sigma'_v$ per unit area is this law with $a = 0$ and
$\delta = \arctan(F^{*}\alpha)$. Written out, FHWA's nominal pullout resistance of a layer is
$P_r = F^{*}\alpha\,\sigma'_v L_e C R_c$, where $C = 2$ counts the two bearing faces of a sheet and $R_c$ is the
fraction of the wall the reinforcement covers. For a continuous geosynthetic ($R_c = 1$) that is the integral
above, term for term: the factor of two is already in $r(s)$, and the per-unit-width convention is what $R_c = 1$
means. In FHWA's Example E1 — a 20 ft geogrid-reinforced wall — the geogrids take $F^{*} = 0.45$ and
$\alpha = 0.8$ from the manual's Table 3-6 ($\alpha$ is 0.8 for geogrids, 0.6 for geotextiles, 1.0 for metallic
reinforcement), so the two columns are Adhesion = 0 and
Delta = $\arctan(0.45 \times 0.8) = 19.80°$, and nothing else about the bond is entered. Reading the envelope
where the design failure surface crosses each of that wall's eleven layers reproduces the manual's whole pullout
table; the entry is under
[published problems](../verification/published.md#fhwa-e1).

**Grouted tiebacks with a bonded length.** A tieback develops pullout resistance only over its grouted bond
length; the sleeved unbonded length transfers no load to the ground. Its envelope is the development-length law,
with $L_{p1} = 0$ at the head and the ramp $L_{p2}$ at the bonded end, $T_{max}$ divided by the load transfer per
unit length. The overburden law would accumulate resistance from the head along the unbonded length. The entries
are under [Tieback](../usage/modeling_reinforcement.md#tieback-grouted-ground-anchor).

### Per-unit-width convention and spacing

All LEM forces are per unit width of slope. Geosynthetic properties are already per unit width (kN/m or lb/ft), so
for them the **Spacing** column is left blank (or 1). Discrete supports — nails, tiebacks — have per-element
capacities (kN per nail) installed at a horizontal spacing $S$; enter the per-element values and the spacing, and
xslope divides all capacity terms ($T_{max}$, $T_{res}$, $T_{end1}$, $T_{end2}$, and the FEM stiffness $EA$)
by $S$. The forces xslope reports back — the LEM line forces, the FEM axial forces, and the
reinforcement-force colorbar — are likewise per unit width; multiply by the spacing $S$ to recover the
per-element (per-nail) force for comparison against a per-element capacity.

## Force Direction (Dir)

Where a reinforcement line crosses the base of a slice at point $r = (x_r, y_r)$, a force of magnitude $T(x_r)$ is
applied to the sliding mass at that point. In the slice free-body diagram it is the force $P$, drawn at the angle
$\psi$ measured from the horizontal — the same reference the slice base inclination $\alpha$ is measured from.
(The $T$ in that diagram is the tension-crack water force, a separate quantity.)

![slice_adv.png](images/slice_adv.png)

The **Dir** setting sets $\psi$:

- **Tangent to slip surface** ($\psi = \alpha$) — the **default**. Flexible
  reinforcement cannot resist bending; as the sliding mass moves, the reinforcement deforms with it and the force
  reorients tangent to the slip surface, whatever the line's own inclination. The angle between the force and the
  base, $\alpha - \psi$, is then zero. This is the appropriate (and conservative) assumption for geotextiles and
  geogrids, and is discussed by Duncan & Wright (2005).
- **Axial** ($\psi$ = the inclination of the reinforcement line itself) — rigid supports such as soil nails,
  grouted tiebacks, and anchored bars carry their force along their own axis; the soil cannot reorient them.

The direction affects each solution method the same way the pile force does: the force is resolved into components
normal and tangential to the slice base — $P\sin(\alpha - \psi)$ normal (zero for tangent) and
$P\cos(\alpha - \psi)$ tangential — and for moment-based methods it contributes a moment about the circle center
through its real moment arm at point $r$.

![The force direction for Dir = Tangent and Dir = Axial at a slice base](images/reinf_direction.png){width=874}

Tangent reinforcement acts with its whole magnitude along the base and has no component across it. Axial
reinforcement has a smaller component along the base, and the remainder presses the sliding mass onto the base,
where it adds frictional resistance $P\sin(\alpha - \psi)\tan\phi$. Which of the two gives the larger factor of
safety therefore depends on $\phi$ and on the angle between the line and the surface it crosses. For tangent reinforcement on a circular surface the force is tangent to
the circle and its moment arm is exactly $R$, which is why the classical formulation reduces to a bare $\sum P$ in
the OMS and Bishop denominators. The per-method equations are given on the
[OMS](oms.md), [Bishop](bishop.md), [Janbu](janbu.md), [force equilibrium](force_eq.md), [Spencer](spencer.md),
and [Morgenstern-Price](mprice.md) pages.

## Force Application (Appl)

Two conventions exist for how a support force enters the factor of safety, and published solutions use both,
so the choice is set per line:

- **Active** (Slide2's "Method A", the **default**): the force is a known, *allowable* working load. It is applied
  to the driving side of the equilibrium equations and is **not** divided by $F$ — the factor of safety applies to
  the soil strength only. Appropriate for pre-tensioned supports (tiebacks) and whenever the entered capacity
  already carries its own safety factor.
- **Passive** (Slide2's "Method B"): the force is an *ultimate* capacity that mobilizes together with the soil
  strength. It is added to the resisting side and **is** divided by $F$. Appropriate when the support only develops
  force as the soil deforms (nails, geosynthetics in some formulations) and the entered capacity is unfactored.

The distinction matters numerically: on the classic Duncan & Wright tieback example (their Fig. 6.34), the same
9,000 lb/ft support gives FS = 1.51 active and FS = 1.32 passive. It also changes what you should enter in the
$T_{max}$ column — an **allowable** force for active, an **ultimate** force for passive.

## Support Type Presets

The **Type** column fills Dir and Appl automatically (either can be overridden by typing over the value):

| Type | Dir | Appl | Typical use |
|---|---|---|---|
| Geosynthetic | Tangent | Active | geotextile / geogrid layers |
| Nail | Axial | Passive | drilled and grouted soil nails |
| Tieback | Axial | Active | pre-tensioned grouted anchors |
| Anchor | Axial | Active | end-anchored bars |

Leave Type blank for a generic tensile line with the defaults (Tangent, Active). A micropile, pile or pier resists by
shear and bending rather than tension and is entered on the [piles](piles.md) sheet.

## LEM vs. FEM

Both engines use the same reinforcement lines, but the mechanics differ:

- **LEM** applies the capacity envelope value as a *prescribed* force at the crossing point, in the Dir direction,
  factored per Appl. The residual strength $T_{res}$ is not used — LEM has no strain compatibility, so there is no
  notion of an element loading past peak.
- **FEM** models each line as tension-only truss elements whose force *emerges* from displacement compatibility;
  the same capacity envelope caps each element's allowable force. An element that reaches it yields and holds that
  force (elastic-perfectly-plastic) — unless $T_{res}$ has been filled in, in which case it drops to that residual
  where the residual is the lower of the two, and holds the envelope value where the envelope is. Both engines
  therefore treat bond slip the same way; what $T_{res}$ adds in the FEM is rupture of the reinforcement itself.
  Dir and Appl have no meaning in the FEM.
  See [Soil Reinforcement in FEM](../fem/reinforcement.md).

## Typical Anchorage Capacities

Approximate ranges for the $T_{end}$ columns, for preliminary estimates only:

| End condition | Typical capacity | Notes |
|---|---|---|
| Soil nail bearing plate | 50-150 kN (10-35 kip) per nail | plate punching or facing flexure governs |
| Geosynthetic facing connection | 30-80% of $T_{max}$ | per connection test data (wrap-around, bodkin, panel) |
| Tieback anchor head / connection | — | enters through $T_{max}$ with $L_{p1} = 0$ (see grouted tiebacks under [Pullout from the effective overburden](#pullout-from-the-effective-overburden)) |
| Free (no plate) | 0 | the friction-only default |

Capacities are per element; with a Spacing entry they are converted to per-unit-width automatically.

## References

Duncan, J.M., & Wright, S.G. (2005). *Soil Strength and Slope Stability*. John Wiley & Sons.

Rocscience Inc. *Slide2 Documentation — Support: Active/Passive Force Application; Define Support Properties.*

Wright, S.G. (1999). *UTEXAS4 — A Computer Program for Slope Stability Calculations.* Shinoak Software, Austin.
