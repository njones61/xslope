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

