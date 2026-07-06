# Sample Problems - Finite Element Method

The following examples illustrate how to use XSLOPE's finite element capabilities for slope stability analysis using
the Shear Strength Reduction Method (SSRM). Each of the Excel input files below can be uploaded and used with the following Google Colab notebook which has been set up specifically for running FEM slope stability analyses:

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_fem.ipynb" target="_"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

The FEM implementation is described in the [FEM Overview](overview.md) page.

### 1. Reinforced Slope with Geogrid Reinforcement

This problem features an engineered slope with six layers of geogrid reinforcement. It is the FEM counterpart of
the LEM reinforced slope example described in the [LEM Samples](../lem/samples.md) page (Problem 9). The slope
geometry and soil properties are the same, with the addition of elastic modulus and Poisson's ratio for the FEM
analysis:

| Property | Shell | Base |
|----------|-------|------|
| Cohesion, $c$ (psf) | 300 | 0 |
| Friction angle, $\phi$ (degrees) | 37 | 37 |
| Unit weight, $\gamma$ (pcf) | 130 | 130 |
| Young's modulus, $E$ (psf) | 1,000,000 | 1,000,000 |
| Poisson's ratio, $\nu$ | 0.3 | 0.3 |

A 240 psf surcharge is applied along the slope crest from $x = 30$ ft to $x = 100$ ft.

Six reinforcement lines are defined with the following properties:

| Property | Value |
|----------|-------|
| $T_{max}$ | 800 lb/ft |
| $T_{res}$ | 600 lb/ft |
| $L_{p1}$, $L_{p2}$ | 4 ft |
| $E$ | 800,000 psf |
| $Area$ | 0.1 ft$^2$/ft |
| $EA$ | 80,000 lb/ft |

Excel input file: [xslope_reinforce_fem.xlsx](files/xslope_reinforce_fem.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![reinforce_fem_inputs.png](images/reinforce_fem_inputs.png){width=1000}

FEM mesh with boundary conditions and reinforcement elements (red lines):

![reinforce_fem_mesh.png](images/reinforce_fem_mesh.png){width=1000}

SSRM results. The computed factor of safety is **FS = 1.67**, in good agreement with the
companion LEM analysis (fixed reinforcement forces), which gives **FS = 1.59** by Spencer's
method (see [LEM sample problem 9](../lem/samples.md)) — the FEM reads ~5% above the LEM,
typical of reinforced systems where the FE solution mobilizes reinforcement through
deformation rather than assuming fixed forces. Long-run equilibrium verification confirms a
genuine failure boundary: at F = 1.6 the viscoplastic field settles to true equilibrium
(displacement constant to six digits), while at F = 1.75 it creeps indefinitely.

The plots below show the solution at the computed factor of safety (**F = 1.67**). The
top plot shows the deformed mesh with original and deformed reinforcement positions. The
middle plot shows the viscoplastic shear strain concentration with reinforcement elements
colored by axial force (blue = low, red = high); green elements are inactive (no tension)
and black elements at the ends have pulled out. The bottom plot shows the displacement
vectors. The reinforcement summary table is shown below.

![reinforce_fem_results.png](images/reinforce_fem_results.png){width=1000}

Reinforcement summary:

```bash
=== Reinforcement Summary ===
Line  Elems     Max T     Avg T  Tension  In Lp  At Tres  Broken  Status
--------------------------------------------------------------------------------
   1      9     722.8     497.9        7      4        3       2  YIELDED
   2      9     663.1     569.0        6      4        4       2  YIELDED
   3      9     644.0     489.4        7      4        4       2  YIELDED
   4      9     774.7     524.2        7      4        4       2  YIELDED
   5      9     696.1     523.8        7      4        4       2  YIELDED
   6      9     725.5     596.4        7      4        0       1  PULLOUT
--------------------------------------------------------------------------------

  PULLOUT: Elements near the reinforcement ends (within Lp) have failed due to insufficient embedment length. Interior elements are intact.
  YIELDED: One or more elements have exceeded Tallow and dropped to residual capacity Tres. The line is still carrying load at reduced strength.
```

<!-- test: file=files/xslope_reinforce_fem.xlsx, type=fem_ssrm, expected_fs=1.67, element_type=tri6, target_size=2, tolerance=0.01, f_min=1.2, f_max=1.9, max_iter=4000 -->

At the factor of safety, the reinforcement is heavily mobilized: lines 1-5 have yielded
(interior elements at residual capacity $T_{res}$ = 600 lb/ft) and line 6 shows end pullout,
with peak forces of 644-775 lb/ft approaching $T_{max}$ = 800 lb/ft. This is the expected
state at incipient failure — the system fails when the soil's reduced strength and the
reinforcement's residual capacity can no longer balance the driving forces together.

### 2. Slope Stabilized with Drilled Shaft Piles

This problem features a 1:1 slope in a medium-stiff clay stabilized by two rows of drilled shafts.

![piles_fem_inputs.png](images/piles_fem_inputs.png){width=1000}

Excel input file: [xslope_piles_fem.xlsx](files/xslope_piles_fem.xlsx)

The soil properties are:

| Property | Value |
|----------|-------|
| Cohesion, $c$ | 200 psf |
| Friction angle, $\phi$ | 20 degrees |
| Unit weight, $\gamma$ | 120 pcf |
| Young's modulus, $E$ | 2,000,000 psf |
| Poisson's ratio, $\nu$ | 0.3 |

Two rows of vertical drilled shafts are placed at $x = 5$ ft and $x = 10$ ft along the slope face, both extending from the ground surface to $y = -10$ ft:

| Property | Value |
|----------|-------|
| Diameter, $D$ | 2.0 ft |
| Spacing, $S$ | 6.0 ft |
| Young's modulus, $E_{\text{pile}}$ | 518,400,000 psf (concrete, $f'_c$ = 4000 psi) |
| Moment of inertia, $I$ | $\pi D^4 / 64$ = 0.785 ft$^4$ (auto-computed from $D$) |
| Cross-sectional area, $A$ | $\pi D^2 / 4$ = 3.14 ft$^2$ (auto-computed from $D$) |
| Shear capacity, $V_{\text{cap}}$ | 46,000 lb |
| Moment capacity, $M_{\text{cap}}$ | 60,000 ft·lb |
| Fixity | free |

Each pile is modeled as a chain of 6-DOF Euler-Bernoulli beam elements with rotational DOFs at each node (see [FEM Piles](piles.md) for the formulation). The pile stiffness ($EI$ and $EA$) is scaled by $1/S$ to convert from per-pile to per-unit-width quantities. Unlike the LEM approach where the user provides a single force $H$, the FEM beam elements naturally develop resistance as the soil deforms around the pile. Bending moments are computed directly at each node, and structural capacity limits ($V_{\text{cap}}$, $M_{\text{cap}}$) are enforced through the viscoplastic correction loop.

FEM mesh with boundary conditions. The piles are shown as green line elements along the pile axes:

![piles_fem_mesh.png](images/piles_fem_mesh.png){width=1000}

SSRM results without piles (**FS = 1.18**). The shear strain concentration shows a failure mechanism passing through the toe:

![piles_fem_results_no_pile.png](images/piles_fem_results_no_pile.png){width=1000}

SSRM results with two rows of piles (**FS = 1.38**). The pile elements are colored by lateral (shear) force in the shear strain plot. The piles resist the sliding mass and the failure mechanism is modified by their presence:

![piles_fem_results.png](images/piles_fem_results.png){width=1000}

Pile summary:

```text
=== Pile Summary ===
Pile  Elems   Max |T|   Max |V|   Max |M|     V_cap     M_cap  Yielded  Status
--------------------------------------------------------------------------------
   1      7     482.1    2116.7    6352.4    7666.7   10000.0    0/7  OK
   2      9    1280.6    2323.6    7227.6    7666.7   10000.0    0/9  OK
--------------------------------------------------------------------------------
```

The two rows of piles increase the factor of safety from 1.18 to 1.38 — a 17% improvement. The maximum bending moment (7228 per unit width in Pile 2) reaches about 72% of the moment capacity ($M_{\text{cap}}/S$ = 10,000), and the maximum shear about 30% of $V_{\text{cap}}/S$, so the structural capacity does not govern for this problem. The soil's ability to transfer lateral load to the piles is the limiting factor, not the pile strength.

This is typical behavior for piles in relatively weak soil — the pile is much stiffer than the surrounding soil, and increasing the pile diameter or stiffness beyond a certain point produces diminishing returns. The 2D plane-strain model also does not capture the three-dimensional soil arching between piles that the Ito & Matsui theory accounts for in LEM, which can make the FEM result more conservative than the LEM result.

<!-- test: file=files/xslope_piles_fem.xlsx, type=fem_ssrm, expected_fs=1.38, element_type=tri6, target_size=2, tolerance=0.01, f_min=1.0, f_max=1.6, max_iter=4000 -->

### 3. Non-Circular Failure Surface with Thin Weak Layer

This is the FEM counterpart of the LEM non-circular failure surface example described in the [LEM Samples](../lem/samples.md)
page (Problem 7). The problem features a thin weak clay layer in the foundation of a slope, which controls the
failure mechanism. This problem was also featured in the user manual for the UTEXASED slope stability analysis
software developed by Stephen G. Wright at the University of Texas at Austin.

![noncircular.png](../lem/sample_images/noncircular.png){width=900}

The slope geometry and strength properties are the same as the LEM problem. Young's modulus ($E$) and Poisson's
ratio ($\nu$) are estimated from typical correlations for each soil type:

| Soil | $c'$ (psf) | $\phi'$ (deg) | $\gamma$ (pcf) | $E$ (psf) | $\nu$ |
|------|:----------:|:--------------:|:---------------:|:---------:|:-----:|
| Sand Fill | 0 | 37 | 120 | 1,000,000 | 0.30 |
| Sand | 0 | 33 | 123 | 700,000 | 0.30 |
| Soft Clay ($S_u$ = 200) | 0 ($\phi = 0$) | 0 | 118 | 60,000 | 0.40 |
| Dense Sand | 0 | 37 | 131 | 1,500,000 | 0.28 |

The soft clay is modeled as an undrained material ($\phi = 0$) with $E/S_u \approx 300$. A Poisson's ratio of 0.40
is used rather than the theoretical undrained value of 0.5 to avoid numerical issues with near-incompressibility.

Excel input file: [xslope_noncircular_fem.xlsx](files/xslope_noncircular_fem.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![non_circ_inputs.png](images/non_circ_inputs.png){width=1000}

FEM mesh with boundary conditions and material zones. **Mesh resolution matters for this
problem**: the soft clay layer is only 2 ft thick, and the mesh must place at least two
elements through its thickness to resolve the shear band that controls the failure
mechanism — a target element size of 1.0 ft (or finer) is required. A coarser mesh
stiffens the thin layer artificially and distorts the strain field within it:

![non_circ_mesh.png](images/non_circ_mesh.png){width=1000}

SSRM results. The computed factor of safety is **FS = 1.68**. The plots show the solution
at the computed factor of safety. The middle plot shows the viscoplastic shear strain
concentration, which clearly reveals the non-circular failure mechanism passing through the
thin weak clay layer — matching the expected behavior without any prior assumption about
the failure surface shape. The bottom plot shows the displacement vectors, confirming
lateral sliding of the slope mass along the clay layer.

![non_circ_results.png](images/non_circ_results.png){width=1000}

The FEM result of FS = 1.68 is about 3% below the LEM result of FS = 1.74 obtained using
Spencer's method — both analyses use the same piezometric surface in the foundation sand.
Differences of this order between SSRM and LEM are typical: the FEM develops the failure
mechanism freely through the global stress field, while the LEM evaluates rigid-block
equilibrium on a prescribed surface, and the two methods answer subtly different questions.
The FS shows a mild residual mesh sensitivity characteristic of thin-shear-band
localization (1.70 / 1.68 / 1.67 at target sizes 2.0 / 1.0 / 0.75): the finer the mesh, the
more sharply the band through the 2-ft layer is resolved.

<!-- mesh resolution: the 2-ft soft clay layer needs >=2 elements through its thickness;
     target_size=1.0 or finer (ts=2.0 gives 1.697, ts=1.0 gives 1.684, ts=0.75 gives
     1.672 — mild thin-band localization sensitivity) -->
<!-- test: file=files/xslope_noncircular_fem.xlsx, type=fem_ssrm, expected_fs=1.68, element_type=tri6, target_size=1, tolerance=0.01, f_min=1.4, f_max=2.2, max_iter=4000 -->


---

The remaining problems are **verification benchmarks**: published or
analytically-anchored cases used to validate the FEM-SSRM implementation.
Each is locked into the automated regression suite. See also the
[Verification](../verification.md) page.

### 4. Verification: Griffiths & Lane (1999) Example 1 — Homogeneous Slope {#verification-griffiths1}

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

Excel input file: [xslope_griffiths1.xlsx](files/xslope_griffiths1.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![griffiths1_inputs.png](images/griffiths1_inputs.png){width=1000}

FEM mesh with boundary conditions. Fixed supports (triangles) at the base, x-rollers (circles) on the sides:

![griffiths1_mesh.png](images/griffiths1_mesh.png){width=1000}

SSRM results. The computed factor of safety is **FS = 1.36** under XSLOPE's strict
true-equilibrium convergence criterion, with the displacement-vs-F upturn at **F ≈ 1.40** —
bracketing the published values: [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387)
report FE FOS = 1.4 (their tolerant convergence check accepts slow residual creep that
XSLOPE's equilibrium criterion rejects; their Table 2 converges at F = 1.35 and fails at
1.40), and the [Bishop & Morgenstern (1960)](https://doi.org/10.1680/geot.1960.10.4.129)
stability chart gives 1.380. All three readings agree within ±3%. The plots below show the
solution at the computed factor of safety (F = 1.36). The top plot shows the deformed mesh;
the middle plot shows the viscoplastic shear strain concentration, which reveals the circular
failure mechanism without any prior assumption about its shape or location; the bottom plot
shows the displacement vectors.

![griffiths1_results.png](images/griffiths1_results.png){width=1000}

The displacement-versus-F sweep — the failure evidence Griffiths & Lane themselves present
(their Fig. 2) — shows the upturn exactly at F ≈ 1.40: the maximum displacement is
essentially flat through F = 1.35, then grows by 3× between F = 1.40 and 1.45 and an
order of magnitude by F = 1.6.

![griffiths1_sweep.png](images/griffiths1_sweep.png){width=700}

This benchmark also appears on the
[Verification](../verification.md#finite-element-slope-stability-ssrm) page.

<!-- test: file=files/xslope_griffiths1.xlsx, type=fem_ssrm, expected_fs=1.36, element_type=quad8, target_size=3.5, tolerance=0.01, f_min=1.0, f_max=1.8, max_iter=4000, benchmark=SSRM-1 -->
<!-- Element-type coverage: SSRM on each quadratic type (tri6, quad8, quad9). Slower (SSRM x3), so benchmark-gated. -->
<!-- test: file=files/xslope_griffiths1.xlsx, type=fem_elements, expected_fs=1.36, tolerance=0.04, target_size=3.5, f_min=1.0, f_max=1.8, max_iter=4000, benchmark=SSRM-elements -->

### 5. Verification: Griffiths & Lane (1999) Example 6 — Two-Sided Earth Dam {#verification-griffiths6}

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
[xslope_griffiths6_full.xlsx](files/xslope_griffiths6_full.xlsx) (reservoir full),
[xslope_griffiths6_dry.xlsx](files/xslope_griffiths6_dry.xlsx) (before filling)

![griffiths6_full_inputs.png](images/griffiths6_full_inputs.png){width=1000}

Results:

| Case | XSLOPE FOS | G&L FOS | Diff |
|---|---|---|---|
| Full reservoir (free surface) | 1.91 | ~1.9 | +1% |
| Before filling (no free surface) | 2.44 | ~2.4 | +2% |

Solution for the before-filling (dry) case at the computed factor of safety (F = 2.44). The
shear strain concentration and displacement vectors show the critical mechanism passing
beneath the crest and exiting on the downstream face:

![griffiths6_dry_results.png](images/griffiths6_dry_results.png){width=1000}

Solution for the full-reservoir case at the computed factor of safety (F = 1.91). With the
free surface in place, the downstream slope is the weaker side: the shear strain band runs
from the crest to the downstream toe, and the displacement vectors show the rotational
sliding mass — the same surface found by Griffiths & Lane and by XSLOPE's own Spencer
analysis. The wet case uses quadratic triangles (tri6): the submerged upstream skin
carries small persistent stresses near the yield surface, and the quad8 element's
reduced-integration hourglass mode is susceptible to such forcing (see the
[FEM Overview](overview.md) discussion of submerged boundaries).

![griffiths6_full_results.png](images/griffiths6_full_results.png){width=1000}

The wet case is a strong test of the pore-pressure treatment. Under the
effective-stress formulation with consistently integrated boundary loads, the
submerged soil simply carries its buoyant weight: a solve at F = 1 converges in
a handful of iterations with an essentially elastic strain field (flooded
ground at working strength sits quietly — a sanity check worth running on any
submerged model), and the failure boundary emerges sharply at F = 1.91 under
the default non-convergence criterion. The agreement with limit equilibrium is
striking: XSLOPE's own Spencer analysis of the same section gives 1.915 (vs the
paper's limit-equilibrium 1.90), with the same downstream critical surface, and
the relative reservoir effect matches the paper (wet/dry = 0.78 vs 0.79). See
the [Verification](../verification.md#finite-element-slope-stability-ssrm) page.

<!-- test: file=files/xslope_griffiths6_dry.xlsx, type=fem_ssrm, expected_fs=2.45, element_type=quad8, target_size=2, tolerance=0.01, f_min=2.0, f_max=2.8, max_iter=4000, benchmark=SSRM-2 -->
<!-- test: file=files/xslope_griffiths6_full.xlsx, type=fem_ssrm, expected_fs=1.91, element_type=tri6, target_size=2, tolerance=0.01, f_min=1.6, f_max=2.2, max_iter=4000, benchmark=SSRM-2 -->
