# Finite-Element Slope Stability (SSRM) Benchmarks

The SSRM implementation uses the Smith & Griffiths 4-component plane-strain
Mohr-Coulomb viscoplastic formulation. The factor of safety is found by
bisection on the **equilibrium (non-convergence) criterion** — the default —
which brackets the strength reduction at which the viscoplastic iteration can
no longer reach true equilibrium. The displacement-versus-$F$ catastrophe sweep
is available as a secondary diagnostic and is shown below for Example 1. Pore
pressures (where present) enter through the effective-stress formulation, and
reservoir loads are applied as consistent boundary tractions, so submerged
problems converge without any special criterion (see
[FEM Overview](../fem/overview.md)).

Full bibliographic details for the author-year citations on this page are on the
shared [References](references.md) page.

### Griffiths & Lane (1999) Example 1 — Homogeneous Slope {#verification-griffiths1}

This is the benchmark problem from Griffiths & Lane (1999), "Slope stability analysis by finite elements,"
*Geotechnique*, 49(3), 387-403. It features a homogeneous slope with the following properties:

| Property | Value |
|----------|-------|
| Cohesion, $c$ | 312.5 psf |
| Friction angle, $\phi$ | 20 degrees |
| Unit weight, $\gamma$ | 125 pcf |
| Young's modulus, $E$ | 700,000 psf |
| Poisson's ratio, $\nu$ | 0.3 |

The dimensionless parameter $c/\gamma H = 0.05$ with $\phi = 20°$ gives an expected factor of safety of
approximately 1.4 (Griffiths & Lane, 1999, Table 1).

Excel input file: [xslope_griffiths1.xlsx](../fem/files/xslope_griffiths1.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![griffiths1_inputs.png](../fem/images/griffiths1_inputs.png){width=1000}

FEM mesh with boundary conditions. Fixed supports (triangles) at the base, x-rollers (circles) on the sides:

![griffiths1_mesh.png](../fem/images/griffiths1_mesh.png){width=1000}

SSRM results. The computed factor of safety is **FS = 1.35** under XSLOPE's strict
true-equilibrium convergence criterion, with the displacement-vs-F upturn at **F ≈ 1.40** —
bracketing the published values: [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387)
report FE FOS = 1.4 (their tolerant convergence check accepts slow residual creep that
XSLOPE's equilibrium criterion rejects; their Table 2 converges at F = 1.35 and fails at
1.40), and the [Bishop & Morgenstern (1960)](https://doi.org/10.1680/geot.1960.10.4.129)
stability chart gives 1.380. All three readings agree within ±4%. The plots below show the
solution at the computed factor of safety (F = 1.35). The top plot shows the deformed mesh;
the middle plot shows the viscoplastic shear strain concentration, which reveals the circular
failure mechanism without any prior assumption about its shape or location; the bottom plot
shows the displacement vectors.

![griffiths1_results.png](../fem/images/griffiths1_results.png){width=1000}

The displacement-versus-F sweep — the failure evidence Griffiths & Lane themselves present
(their Fig. 2) — shows the upturn exactly at F ≈ 1.40: the maximum displacement is
essentially flat through F = 1.35, then grows by 3× between F = 1.40 and 1.45 and an
order of magnitude by F = 1.6.

![griffiths1_sweep.png](../fem/images/griffiths1_sweep.png){width=700}

This benchmark also appears on the
[Verification](../verification/ssrm.md) page.

<!-- test: file=../fem/files/xslope_griffiths1.xlsx, type=fem_ssrm, expected_fs=1.35, element_type=quad8, target_size=3.5, tolerance=0.01, f_min=1.0, f_max=1.8, max_iter=16000, benchmark=SSRM-1 -->
<!-- Element-type coverage: SSRM on each quadratic type (tri6, quad8, quad9). Slower (SSRM x3), so benchmark-gated. -->
<!-- test: file=../fem/files/xslope_griffiths1.xlsx, type=fem_elements, expected_fs=1.36, tolerance=0.04, target_size=3.5, f_min=1.0, f_max=1.8, max_iter=4000, benchmark=SSRM-elements -->
<!-- SSRM auto-bracketing: a deliberately-wrong [F_min,F_max] must still find the FS. Coarse tri6 mesh (~1.39, fast) so these run un-gated. -->
<!-- test: file=../fem/files/xslope_griffiths1.xlsx, type=fem_ssrm, expected_fs=1.39, element_type=tri6, target_size=6, tolerance=0.05, f_min=1.5, f_max=1.9, max_iter=4000 -->
<!-- test: file=../fem/files/xslope_griffiths1.xlsx, type=fem_ssrm, expected_fs=1.39, element_type=tri6, target_size=6, tolerance=0.05, f_min=0.5, f_max=0.9, max_iter=4000 -->

### Griffiths & Lane (1999) Example 2 — Homogeneous Slope with a Foundation Layer {#verification-griffiths2}

This is Example 2 of [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387)
(their Fig. 5): the Example 1 slope with a foundation layer of the **same soil** added
beneath it, thickness $H/2$, so the firm base now sits at depth $D = 1.5\,H$ below the
crest. Every material property and the slope geometry are unchanged from Example 1 — only
the foundation is added.

| Property | Value |
|----------|-------|
| Cohesion, $c'$ | 312.5 psf |
| Friction angle, $\phi'$ | 20 degrees |
| Unit weight, $\gamma$ | 125 pcf |
| Young's modulus, $E$ | 700,000 psf |
| Poisson's ratio, $\nu$ | 0.3 |
| Foundation | $H/2$ = 25 ft of the same soil below the toe ($D = 1.5$) |

The dimensionless strength is again $c'/\gamma H = 0.05$ with $\phi' = 20°$. The domain uses
the printed dimensions of Fig. 1 — a $1.2\,H$ (60 ft) crest platform, the 2:1 face, and a
$2\,H$ (100 ft) runout past the toe — with the $H/2$ foundation from the Example 2 text
($D = 1.5$). (Example 1 as modelled here omits the runout because with the base at toe
level a mechanism cannot reach it; with a foundation beneath the toe the runout gives any
deeper surface room to form.)

Excel input file: [xslope_griffiths2.xlsx](../fem/files/xslope_griffiths2.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![griffiths2_inputs.png](../fem/images/griffiths2_inputs.png){width=1000}

FEM mesh with boundary conditions. Fixed supports (triangles) at the base, x-rollers
(circles) on the sides:

![griffiths2_mesh.png](../fem/images/griffiths2_mesh.png){width=1000}

SSRM results. The computed factor of safety is **FS = 1.33** under XSLOPE's strict
true-equilibrium convergence criterion, with the displacement upturn at **F ≈ 1.4** — the
same value as Example 1 (there FS = 1.35), confirming Griffiths & Lane's central point:
adding the foundation layer does **not** change the factor of safety, because the critical
mechanism stays a toe failure. Griffiths & Lane report FE FOS ≈ 1.4, unchanged from
Example 1. The plots below show the solution at the developed failure state. The middle
panel — the viscoplastic shear-strain concentration — is the key image: the shear band
runs from the crest and **exits at the toe** (x = 160, y = 0), staying well above the
foundation base at y = −25. The mechanism finds the toe on its own, with no assumption
about the slip surface.

![griffiths2_results.png](../fem/images/griffiths2_results.png){width=1000}

**The false base circle.** This example is Griffiths & Lane's signature demonstration that
finite elements find the true mechanism where a limit-equilibrium search can be misled. For
this slope-plus-foundation, the Bishop & Morgenstern (1960) base-circle charts give
FOS = 1.752, and a proprietary slip-circle program returned FOS = 1.7 when its failure
circle was forced tangent to the base of the foundation — both well above the true value.
Cousins' (1978) toe-circle charts, on the other hand, agree with the FE result at 1.4. The
lesson is that a limit-equilibrium method requires the user to steer the search onto the
correct family of surfaces; assume a base circle and the answer is 25% unconservative.

XSLOPE's own limit-equilibrium search reproduces both sides of this exactly. An
unconstrained global circular search (Spencer, grid seed) settles on a **toe circle** —
critical centre (143.3, 114.1), lowest point y = −1.2, just below the toe and far above the
foundation base — at **FS = 1.37**, agreeing with the SSRM result and with Cousins' toe
charts. When the same search is confined to circles tangent to the foundation base (the
assumption behind the chart value), it returns **FS = 1.70** — the paper's false base-circle
result. XSLOPE's default search is not misled: left free, it finds the toe.

<!-- test: file=../fem/files/xslope_griffiths2.xlsx, type=fem_ssrm, expected_fs=1.334, element_type=quad8, target_size=3.5, tolerance=0.01, f_min=1.0, f_max=1.8, max_iter=16000, benchmark=SSRM-G2 -->
<!-- Coarse tri6 quick SSRM (ungated): confirms the foundation layer leaves the toe-failure FS unchanged. -->
<!-- test: file=../fem/files/xslope_griffiths2.xlsx, type=fem_ssrm, expected_fs=1.36, element_type=tri6, target_size=6, tolerance=0.05, f_min=1.0, f_max=1.8, max_iter=4000 -->
<!-- LEM teaching point: the unconstrained global search finds the TOE circle (true mechanism, ~1.37); forcing tangency to the foundation base reproduces the paper's false base circle (~1.70). -->
<!-- test: file=../fem/files/xslope_griffiths2.xlsx, type=circular_search, method=spencer, seed=grid, num_slices=40, expected_fs=1.366, tolerance=0.02 -->
<!-- test: file=../fem/files/xslope_griffiths2.xlsx, type=circular_search, method=spencer, seed=grid, num_slices=40, tangent_depth=-25;-23, expected_fs=1.702, tolerance=0.02 -->

### Griffiths & Lane (1999) Example 4 — Undrained Clay Slope over a Weak Foundation {#verification-griffiths4}

This is Example 4 of [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387)
(their Fig. 9): an **undrained** ($\phi_u = 0$) clay slope resting on a foundation layer,
with the two layers assigned different undrained strengths. The slope body carries a
constant $c_{u1}/\gamma H = 0.25$; the foundation strength $c_{u2}$ is varied, and the
behaviour is governed by the strength ratio $c_{u2}/c_{u1}$. The firm base sits at depth
$D = 2$ (a full $2H$ below the crest); the foundation layer is $H$ thick, and the material
boundary is the toe elevation, which is also the ground surface of the runout.

| Property | Value |
|----------|-------|
| Slope undrained strength, $c_{u1}$ | 1562.5 psf ($c_{u1}/\gamma H = 0.25$) |
| Foundation undrained strength, $c_{u2}$ | $c_{u2}/c_{u1} \times c_{u1}$ (ratio varied) |
| Friction angle, $\phi_u$ | 0 degrees |
| Unit weight, $\gamma$ | 125 pcf |
| Slope height, $H$ | 50 ft (crest at $y=100$, toe at $y=50$, firm base at $y=0$) |
| Geometry | $2H$ crest platform, 2:1 face, $2H$ runout, $D = 2$ |

The domain is tiled by two explicit material polygons — the slope body $c_{u1}$ and the
foundation $c_{u2}$ — because the layer boundary coincides with the runout surface, so a
profile-line stack would leave a zero-thickness $c_{u1}$ sliver over the runout. Elastic
constants are assigned by soil type from the [FEM Overview](../fem/overview.md) table:
$c_{u1} = 1562.5$ psf classifies Medium Clay ($E = 668{,}300$ psf, $\nu = 0.40$), and the
$c_{u2} = 3125$ psf foundation at ratio 2 classifies Stiff Clay ($E = 2{,}610{,}700$ psf,
$\nu = 0.30$). The SSRM factor of safety is independent of the elastic constants; they only
set the displacement scale.

Two bracket cases are built: **$c_{u2}/c_{u1} = 1$** (the foundation is as weak as the
slope) and **$c_{u2}/c_{u1} = 2$** (the foundation is firmly stronger). Griffiths & Lane's
Fig. 10 shows the computed factor of safety flattening onto a toe-circle plateau for
$c_{u2}/c_{u1} \gtrsim 1.5$; the ratio-2 case sits firmly on that plateau (figure-confirmed),
and the ratio-1 case is the homogeneous base-circle limit.

Excel input files:
[xslope_griffiths4_r1.xlsx](../fem/files/xslope_griffiths4_r1.xlsx) ($c_{u2}/c_{u1} = 1$),
[xslope_griffiths4_r2.xlsx](../fem/files/xslope_griffiths4_r2.xlsx) ($c_{u2}/c_{u1} = 2$)

Inputs plotted with the XSLOPE plot_inputs() function (geometry is identical for both
cases; only the foundation strength differs):

![griffiths4_inputs.png](../fem/images/griffiths4_inputs.png){width=1000}

FEM mesh with boundary conditions. Fixed supports (triangles) at the base, x-rollers
(circles) on the sides:

![griffiths4_mesh.png](../fem/images/griffiths4_mesh.png){width=1000}

Results. The two cases straddle a change of failure mechanism, which is the point of the
example. The published anchors are Taylor's (1937) classical solutions, quoted by Griffiths
& Lane on Fig. 10: FOS = 1.47 for the deep **base circle** at $c_{u2} = c_{u1}$, and
FOS = 2.10 for the shallow **toe circle** at $c_{u2} \gg c_{u1}$.

| Case | $c_{u2}/c_{u1}$ | XSLOPE SSRM (quad8) | Taylor (1937) | Mechanism |
|---|---|---|---|---|
| r1 | 1 | 1.44 | 1.47 (base circle) | deep base |
| r2 | 2 | 2.00 | 2.10 (toe circle) | shallow toe |

Both XSLOPE values sit a few percent below the Taylor anchors — the same offset seen in
Examples 1 and 2, where XSLOPE's strict true-equilibrium convergence criterion rejects the
slow residual creep that a tolerant check accepts. The relative jump between the two cases
(1.44 → 2.00, a factor of 1.39) closely tracks the published jump (1.47 → 2.10, 1.43).

The critical mechanism flips between the two cases, exactly as in Griffiths & Lane's Fig. 11.
For $c_{u2}/c_{u1} = 1$ (base case), the shear-strain concentration dips to the firm base at
$y = 0$ and passes along it — a deep-seated base mechanism (their Fig. 11a):

![griffiths4_r1_results.png](../fem/images/griffiths4_r1_results.png){width=1000}

For $c_{u2}/c_{u1} = 2$ (toe case), the shear band runs from the crest and exits at the toe
($x = 200$, $y = 50$), staying entirely within the weaker upper clay and never entering the
stronger foundation — a shallow toe mechanism (their Fig. 11c). The stronger foundation has
lifted the factor of safety and forced the failure surface up out of the foundation layer:

![griffiths4_r2_results.png](../fem/images/griffiths4_r2_results.png){width=1000}

XSLOPE's own limit-equilibrium search reproduces the same flip. An unconstrained global
circular search (Spencer, grid seed) settles on a **base circle** for the ratio-1 case —
critical surface bottoming at $y \approx 0$, tangent to the firm base — at **FS = 1.47**,
matching Taylor's base-circle chart. For the ratio-2 case the same search settles on a
**toe circle** bottoming at $y \approx 50$ (the toe elevation), confined to the upper clay,
at **FS = 2.02**. The limit-equilibrium method finds the correct mechanism family on its own
here, and the SSRM and Spencer results agree on both the factor of safety and the base→toe
transition.

<!-- test: file=../fem/files/xslope_griffiths4_r1.xlsx, type=fem_ssrm, expected_fs=1.441, element_type=quad8, target_size=3.5, tolerance=0.01, f_min=1.0, f_max=1.8, max_iter=16000, benchmark=SSRM-G4 -->
<!-- test: file=../fem/files/xslope_griffiths4_r2.xlsx, type=fem_ssrm, expected_fs=2.002, element_type=quad8, target_size=3.5, tolerance=0.01, f_min=1.8, f_max=2.4, max_iter=16000, benchmark=SSRM-G4 -->
<!-- Coarse tri6 quick SSRM (ungated): base case (cu2=cu1) and toe case (cu2=2cu1); confirms the mechanism flip lifts the FS from ~1.44 to ~2.0. -->
<!-- test: file=../fem/files/xslope_griffiths4_r1.xlsx, type=fem_ssrm, expected_fs=1.44, element_type=tri6, target_size=6, tolerance=0.05, f_min=1.0, f_max=1.8, max_iter=4000 -->
<!-- test: file=../fem/files/xslope_griffiths4_r2.xlsx, type=fem_ssrm, expected_fs=2.02, element_type=tri6, target_size=6, tolerance=0.05, f_min=1.6, f_max=2.4, max_iter=4000 -->
<!-- LEM companions: the unconstrained global search finds the BASE circle (~1.47, tangent to the firm base) at cu2=cu1 and the TOE circle (~2.02, confined to the upper clay) at cu2=2cu1 — the same base->toe flip as the SSRM and Taylor's charts. -->
<!-- test: file=../fem/files/xslope_griffiths4_r1.xlsx, type=circular_search, method=spencer, seed=grid, num_slices=40, expected_fs=1.468, tolerance=0.02 -->
<!-- test: file=../fem/files/xslope_griffiths4_r2.xlsx, type=circular_search, method=spencer, seed=grid, num_slices=40, expected_fs=2.022, tolerance=0.02 -->

### Griffiths & Lane (1999) Example 6 — Two-Sided Earth Dam {#verification-griffiths6}

The second SSRM verification benchmark, from [Griffiths, D.V. & Lane, P.A. (1999)](https://doi.org/10.1680/geot.1999.49.3.387), *Géotechnique* 49(3),
Example 6: an actual earth dam cross-section (Torres & Coffman, 1997) with
homogenized properties, analyzed both with the reservoir full (free surface
sloping from the upstream face to the downstream toe) and before filling (no
free surface).

| Property | Value |
|---|---|
| Cohesion, $c'$ | 13.8 kPa |
| Friction angle, $\phi'$ | 37° |
| Unit weight, $\gamma$ | 18.2 kN/m³ (above and below the water table) |
| Foundation layer | 7.3 m thick |
| Dam height | 21.3 m above foundation, crest 7.3 m wide |
| Faces | upstream ≈ 18°, downstream ≈ 23° |
| Reservoir | 17.1 m above foundation level |

Pore pressures are taken as $\gamma_w$ × vertical depth below the free surface
(a piezometric line, per the paper), and the reservoir water load is applied as
a normal pressure on the submerged upstream boundary — both exactly as
described by Griffiths & Lane.

Excel input files:
[xslope_griffiths6_full.xlsx](../fem/files/xslope_griffiths6_full.xlsx) (reservoir full),
[xslope_griffiths6_dry.xlsx](../fem/files/xslope_griffiths6_dry.xlsx) (before filling)

![griffiths6_full_inputs.png](../fem/images/griffiths6_full_inputs.png){width=1000}

Results:

| Case | XSLOPE FOS | G&L FOS | Diff |
|---|---|---|---|
| Full reservoir (free surface) | 1.86 | ~1.9 | −2% |
| Before filling (no free surface) | 2.40 | ~2.4 | 0% |

Solution for the before-filling (dry) case at the computed factor of safety (F = 2.40). The
shear strain concentration and displacement vectors show the critical mechanism passing
beneath the crest and exiting on the downstream face:

![griffiths6_dry_results.png](../fem/images/griffiths6_dry_results.png){width=1000}

Solution for the full-reservoir case at the computed factor of safety (F = 1.86). With the
free surface in place, the downstream slope is the weaker side: the shear strain band runs
from the crest to the downstream toe, and the displacement vectors show the rotational
sliding mass — the same surface found by Griffiths & Lane and by XSLOPE's own Spencer
analysis. The wet case uses quadratic triangles (tri6): the submerged upstream skin
carries small persistent stresses near the yield surface, and the quad8 element's
reduced-integration hourglass mode is susceptible to such forcing (see the
[FEM Overview](../fem/overview.md) discussion of submerged boundaries).

![griffiths6_full_results.png](../fem/images/griffiths6_full_results.png){width=1000}

The wet case is a strong test of the pore-pressure treatment. Under the
effective-stress formulation with consistently integrated boundary loads, the
submerged soil simply carries its buoyant weight: a solve at F = 1 converges in
a handful of iterations with an essentially elastic strain field (flooded
ground at working strength sits quietly — a sanity check worth running on any
submerged model), and the failure boundary emerges sharply at F = 1.86 under
the default non-convergence criterion. The agreement with limit equilibrium is
striking: XSLOPE's own Spencer analysis of the same section gives 1.915 (vs the
paper's limit-equilibrium 1.90), with the same downstream critical surface, and
the relative reservoir effect matches the paper (wet/dry = 0.77 vs 0.79). See
the [Verification](../verification/ssrm.md) page.

<!-- test: file=../fem/files/xslope_griffiths6_dry.xlsx, type=fem_ssrm, expected_fs=2.40, element_type=quad8, target_size=2, tolerance=0.01, f_min=2.0, f_max=2.8, max_iter=16000, benchmark=SSRM-2 -->
<!-- test: file=../fem/files/xslope_griffiths6_full.xlsx, type=fem_ssrm, expected_fs=1.858, element_type=tri6, target_size=2, tolerance=0.01, f_min=1.6, f_max=2.2, max_iter=16000, benchmark=SSRM-2 -->

