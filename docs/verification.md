# Verification and Validation

XSLOPE's three analysis modes — limit equilibrium, finite-element seepage, and
finite-element slope stability (SSRM) — are verified against a two-tier
benchmark suite:

1. **Analytical anchors** — problems with exact closed-form solutions, where
   agreement is limited only by discretization.
2. **Established-code and published cross-checks** — problems from the
   verification literature and from codes practitioners already trust, run on
   matched geometry and properties.

All benchmark models, build scripts, and runners are in the repository
(`benchmarks/` and the sample files under `docs/*/files/`), so every number
below can be regenerated. Many are also locked into the automated regression
suite (`run_tests.py`).

---

## Limit Equilibrium

### ACADS simple homogeneous slope (circular search)

Source: [GeoStudio SLOPE/W Verification Manual (Oct 2022)](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf),
ACADS suite (Donald & Giam, 1989; Giam & Donald, 1992). Full problem
description, input file, and figures:
[LEM sample problem 12](lem/samples.md#verification-acads-simple). Single soil, c′ = 3.0 kPa, φ′ = 19.6°,
γ = 20.0 kN/m³; 2:1 slope, 10 m high; automated critical-circle search with 50
slices per method. The ACADS consensus answer is FOS ≈ 1.00.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 0.942 | 1.00 | −5.8% |
| Bishop's Simplified | 0.985 | 1.00 | −1.5% |
| Simplified Janbu | 0.986 | 1.00 | −1.4% |
| Corps of Engineers | 0.990 | 1.00 | −1.0% |
| Lowe & Karafiath | 0.987 | 1.00 | −1.3% |
| Spencer | 0.984 | 1.00 | −1.6% |
| Morgenstern-Price | 0.984 | 1.00 | −1.6% |

All rigorous methods fall within the ACADS accepted band; OMS reads low, as
expected for a legacy method (its known conservative bias on effective-stress
problems is why it is reported for completeness only).

### ACADS weak-layer slope (non-circular search)

Source: [SLOPE/W Verification Manual](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf),
ACADS weak-layer case (sec. 2.7). Full details:
[LEM sample problem 13](lem/samples.md#verification-acads-weak-layer). A 2:1
slope with a thin c′ = 0, φ′ = 10° interlayer and piezometric line; the
critical surface runs along the weak layer. ACADS accepted band ≈ 1.26.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Spencer | 1.258 | ~1.26 | −0.2% |
| Morgenstern-Price | 1.248 | ~1.26 | −1.0% |
| Corps of Engineers | 1.336 | ~1.26 | +6.0% |
| Lowe & Karafiath | 1.249 | ~1.26 | −0.9% |
| Simplified Janbu | 1.278 | ~1.26 | +1.4% |

The non-circular slip is seeded just above the base of the weak layer (the
standard placement for a weak interlayer — see the sample-problem write-up).
This is exactly where complete equilibrium earns its keep: XSLOPE's Spencer and
Morgenstern-Price (half-sine) both land within ~1% of SLOPE/W's named
Morgenstern-Price value (1.261) — Spencer at 1.258 (−0.2%) and M-P at 1.248
(−1.0%).
Corps of Engineers reads modestly high here, consistent with ground-parallel
side-force inclinations on this geometry (XSLOPE uses the standard "Corps #2"
convention — see [Force Equilibrium Methods](lem/force_eq.md)).

### Arai & Tagyo homogeneous slope (circular search)

Source: [Arai & Tagyo (1985)](https://doi.org/10.3208/sandf1972.25.43), Soils
and Foundations 25(1); republished in Greco (1996), Malkawi et al. (2001), Kim
et al. (2002). Full details:
[LEM sample problem 14](lem/samples.md#verification-arai-tagyo). Homogeneous 1.5:1 slope,
20 m high, c = 41.65 kPa, φ = 15.0°, γ = 18.82 kN/m³. Published FOS ≈ 1.451.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 1.344 | 1.451 | −7.4% |
| Bishop's Simplified | 1.404 | 1.451 | −3.2% |
| Simplified Janbu | 1.411 | 1.451 | −2.8% |
| Corps of Engineers | 1.476 | 1.451 | +1.7% |
| Lowe & Karafiath | 1.438 | 1.451 | −0.9% |
| Spencer | 1.401 | 1.451 | −3.4% |
| Morgenstern-Price | 1.400 | 1.451 | −3.5% |

---

## Finite-Element Seepage

### Confined radial flow (analytical anchor — exact)

Full details: [seepage sample problem 7](seep/samples.md#verification-confined-radial).

A quarter-annulus confined flow problem in **plan view** (one quadrant of the
Thiem radial-flow-to-a-well geometry): inner arc (r = 10) at head 30, outer
arc (r = 30) at head 10, straight radial edges as no-flow streamlines. Confined
saturated flow obeys Laplace's equation in head with no gravity term, so the
model-plane orientation is irrelevant to the mathematics. Steady
flow has the exact solution q = k·(π/2)·Δh / ln(r₂/r₁) = 28.596
(k = 1) with a logarithmic head profile.

| Quantity | XSLOPE (tri6) | Exact | Diff |
|---|---|---|---|
| Discharge q | 28.5961 | 28.5960 | <0.01% |
| Max nodal head error | 0.004 | 0 | 0.02% of total drop |

Mesh-converged (identical at 2k and 6k nodes); quad8 gives the same result.
The only error source is faceting of the curved arcs by the polygon boundary.

### Partially penetrating sheetpile (analytical anchor — exact)

Full details: [seepage sample problem 8](seep/samples.md#verification-sheetpile).

Pavlovsky's conformal-mapping solution for a cutoff wall of depth s in a
confined stratum of thickness T (Harr, 1962; Polubarinova-Kochina, 1962).
Boundary heads are 30 upstream / 20 downstream (downstream head at the stratum
top, so pressures are non-negative throughout the vertical section):
q = k·H·K(λ′)/(2·K(λ)) with λ = sin(πs/2T). At s/T = ½ the modulus is
self-dual and q = k·H/2 **exactly**. The closed form was additionally verified
with an independent finite-difference solution of the same boundary-value
problem (agreement ~0.4–0.5% at three penetration ratios).

| Case | XSLOPE q | Exact q | Diff | Head below wall tip |
|---|---|---|---|---|
| s/T = 0.50 | 5.010 | 5.000 | +0.20% | 25.0000 (exact: 25) |
| s/T = 0.75 | 3.412 | 3.403 | +0.27% | 25.0000 (exact: 25) |

The error halves with mesh refinement (set by the r^−½ singularity at the wall
tip) and converges to the exact value from above. The head on the wall plane
below the tip equals (h₁+h₂)/2 exactly, an antisymmetry property of the exact
solution that the FE solution reproduces to four decimals.

### SEEP2D cross-check — Johnson Reservoir (established code)

Full details: [seepage sample problem 4](seep/samples.md#johnson-reservoir).

The Johnson Reservoir zoned earth dam (permeable shell, low-permeability core,
foundation; reservoir at el. 160 ft, tailwater at el. 100 ft) was exported to a
SEEP2D input file — the **exact same tri3 mesh topology, boundary conditions,
and material parameters** — and solved with the original USACE/WES SEEP2D
Fortran program. Identical-mesh comparison over all 2,604 nodes:

| Quantity | XSLOPE | SEEP2D | Diff |
|---|---|---|---|
| Total discharge q (ft³/day per ft) | 1.9575 | 1.9603 | −0.14% |
| Nodal heads | RMS Δh = 0.105 ft | (60-ft head range) | 0.18% |

The largest local head difference (~2 ft) occurs adjacent to the free surface,
where the two codes' unsaturated relative-permeability treatments differ in
detail; the bulk flow field agrees to about 0.1 ft.

---

## Finite-Element Slope Stability (SSRM)

The SSRM implementation uses the Smith & Griffiths 4-component plane-strain
Mohr-Coulomb viscoplastic formulation. The factor of safety is found by
bisection on the **equilibrium (non-convergence) criterion** — the default —
which brackets the strength reduction at which the viscoplastic iteration can
no longer reach true equilibrium. The displacement-versus-$F$ catastrophe sweep
is available as a secondary diagnostic and is shown below for Example 1. Pore
pressures (where present) enter through the effective-stress formulation, and
reservoir loads are applied as consistent boundary tractions, so submerged
problems converge without any special criterion (see
[FEM Overview](fem/overview.md)).

### Griffiths & Lane (1999) Example 1 — homogeneous slope

Full details: [FEM sample problem 4](fem/samples.md#verification-griffiths1).

The canonical FE slope stability benchmark: 2:1 slope, φ′ = 20°,
c′/γH = 0.05, no foundation layer.
[Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387) report
FE FOS = 1.4; the
[Bishop & Morgenstern (1960)](https://doi.org/10.1680/geot.1960.10.4.129)
chart gives 1.380.

| | XSLOPE | Reference | Diff |
|---|---|---|---|
| FOS (equilibrium criterion) | 1.36 | 1.4 (G&L FE) | −2.9% |
| FOS (displacement upturn) | ~1.40 | 1.4 (G&L FE) | ~0% |

The displacement-versus-F sweep shows the failure upturn exactly at F ≈ 1.40:
the maximum displacement is essentially flat through F = 1.35, then grows by
3× between F = 1.40 and 1.45 and an order of magnitude by F = 1.6 — the same
diagnostic Griffiths & Lane present (their Fig. 2).

### Griffiths & Lane (1999) Example 6 — two-sided earth dam

Full details: [FEM sample problem 5](fem/samples.md#verification-griffiths6).

An actual dam cross-section (Torres & Coffman, 1997), homogenized: c′ = 13.8
kPa, φ′ = 37°, γ = 18.2 kN/m³, 7.3-m foundation layer, reservoir 17.1 m above
foundation level with a sloping free surface to the downstream toe (pore
pressures from vertical depth below the free surface; reservoir water applied
as a boundary pressure — both per the paper).

| Case | XSLOPE FOS | G&L FOS | Diff |
|---|---|---|---|
| Full reservoir (free surface) | 1.91 | ~1.9 | +1% |
| Before filling (no free surface) | 2.45 | ~2.4 | +2% |

The wet case runs under the default non-convergence criterion with the
effective-stress pore-pressure formulation and consistent boundary-load
integration: the submerged section converges to true equilibrium in a handful
of iterations at F = 1 (the flooded soil carries its buoyant weight) and fails
sharply above F = 1.91. The input model is independently validated by XSLOPE's
Spencer analysis — 1.915, in essentially exact agreement with the SSRM result —
vs the paper's limit-equilibrium 1.90, and the relative reservoir effect
matches the paper (wet/dry = 0.78 vs 0.79).

---

## References

- Arai, K. & Tagyo, K. (1985). Determination of noncircular slip surface giving
  the minimum factor of safety in slope stability analysis. *Soils and
  Foundations* 25(1).
- Bishop, A.W. & Morgenstern, N. (1960). Stability coefficients for earth
  slopes. *Géotechnique* 10(4).
- Donald, I.B. & Giam, P. (1989). *Soil slope stability programs review*. ACADS.
- GEO-SLOPE / Seequent (2022). *Stability Modeling with GeoStudio — SLOPE/W
  Verification Manual*.
- Griffiths, D.V. & Lane, P.A. (1999). Slope stability analysis by finite
  elements. *Géotechnique* 49(3).
- Harr, M.E. (1962). *Groundwater and Seepage*. McGraw-Hill.
- Polubarinova-Kochina, P.Ya. (1962). *Theory of Ground Water Movement*.
  Princeton University Press.
- Smith, I.M. & Griffiths, D.V. (1998). *Programming the Finite Element
  Method*, 3rd ed. Wiley.
- Sun, et al. (2021). Displacement-catastrophe criterion for strength-reduction
  finite-element slope stability.
- Tracy, F.T. *SEEP2D* (USACE Waterways Experiment Station).


---

## Rocscience Slide2 Verification Corpus

The [Rocscience Slide2 verification manual](https://www.rocscience.com/help/slide2/verification-theory/verification-manuals)
contains 111 slope stability problems drawn from the published literature, each with Slide2's computed
factors of safety and (in most cases) independent reference values from the original authors. XSLOPE is
being verified against this corpus problem by problem: each **built** entry links an XSLOPE input file
reproducing the problem, reports the comparison, and is locked into the automated regression suite via a
test tag. Status values: **built** (input file + verified results below), *covered* (equivalent to an
existing XSLOPE sample), *partial* (some data extracted, some still needed), *blocked* (source data not
yet available), *planned*.

Problems are built from the manual's tabulated data and coordinate-labeled figures; where a problem's
geometry exists only as an unlabeled figure, the original source publication is consulted before the
problem is marked *built* — no digitized guesses are used for benchmark inputs.

<!-- test: file=files/rocscience/vp002.xlsx, type=circular_search, num_slices=40, fs_bishop=1.589, fs_spencer=1.585, fs_janbu=1.495, fs_mprice=1.586, benchmark=VP2 -->
<!-- test: file=files/rocscience/vp003.xlsx, type=circular_search, num_slices=40, fs_bishop=1.403, fs_spencer=1.372, fs_janbu=1.354, fs_mprice=1.371, benchmark=VP3 -->
<!-- test: file=files/rocscience/vp004.xlsx, type=circular_search, num_slices=40, fs_bishop=1.013, fs_spencer=0.989, fs_janbu=0.963, fs_mprice=0.987, benchmark=VP4 -->
<!-- test: file=files/rocscience/vp005.xlsx, type=single_circle, num_slices=60, fs_bishop=1.955, fs_spencer=1.955, fs_janbu=1.965, fs_mprice=1.955, benchmark=VP5 -->
<!-- test: file=files/rocscience/vp006.xlsx, type=single_circle, num_slices=60, fs_bishop=2.206, fs_spencer=2.290, fs_janbu=2.073, fs_mprice=2.299, benchmark=VP6 -->
<!-- test: file=files/rocscience/vp008.xlsx, type=single_noncirc, num_slices=50, fs_spencer=1.276, fs_janbu=1.294, fs_mprice=1.260, benchmark=VP8 -->
<!-- test: file=files/rocscience/vp009.xlsx, type=noncircular_search, num_slices=50, fs_spencer=0.724, fs_janbu=0.718, benchmark=VP9 -->
<!-- test: file=files/rocscience/vp015.xlsx, type=circular_search, num_slices=40, fs_bishop=0.419, fs_spencer=0.422, fs_janbu=0.436, fs_mprice=0.420, benchmark=VP15 -->
<!-- test: file=files/rocscience/vp016.xlsx, type=circular_search, num_slices=40, fs_bishop=1.112, fs_spencer=1.113, fs_janbu=1.122, fs_mprice=1.111, benchmark=VP16 -->
<!-- test: file=files/rocscience/vp017.xlsx, type=circular_search, num_slices=50, fs_oms=1.274, fs_bishop=1.342, fs_spencer=1.340, benchmark=VP17 -->
<!-- test: file=files/rocscience/vp018.xlsx, type=noncircular_search, num_slices=50, fs_spencer=1.033, fs_mprice=1.024, benchmark=VP18 -->
<!-- test: file=files/rocscience/vp019.xlsx, type=circular_search, num_slices=50, fs_bishop=1.448, fs_spencer=1.429, benchmark=VP19 -->
<!-- test: file=files/rocscience/vp020.xlsx, type=circular_search, num_slices=50, fs_bishop=1.086, fs_spencer=1.091, benchmark=VP20-circ -->
<!-- test: file=files/rocscience/vp020.xlsx, type=noncircular_search, num_slices=50, fs_spencer=1.082, benchmark=VP20-noncirc -->
<!-- test: file=files/rocscience/vp023.xlsx, type=circular_search, num_slices=50, fs_oms=1.357, fs_bishop=1.130, benchmark=VP23 -->
<!-- test: file=files/rocscience/vp024.xlsx, type=circular_search, num_slices=50, fs_oms=1.433, fs_bishop=1.433, benchmark=VP24 -->
<!-- test: file=files/rocscience/vp036.xlsx, type=circular_search, num_slices=50, fs_bishop=1.333, benchmark=VP36-fs -->
<!-- test: file=files/rocscience/vp036.xlsx, type=reliability, method=bishop, expected_beta=2.263, tolerance=0.03, benchmark=VP36-beta -->
<!-- test: file=files/rocscience/vp041.xlsx, type=circular_search, num_slices=50, fs_bishop=1.668, fs_spencer=1.670, fs_janbu=1.660, benchmark=VP41 -->
<!-- test: file=files/rocscience/vp045a.xlsx, type=circular_search, num_slices=50, fs_spencer=2.801, benchmark=VP45-mc -->
<!-- test: file=files/rocscience/vp045b.xlsx, type=circular_search, num_slices=50, fs_spencer=2.649, benchmark=VP45-pow -->
<!-- test: file=files/rocscience/vp050.xlsx, type=single_noncirc, num_slices=60, fs_janbu=1.448, fs_spencer=1.576, benchmark=VP50 -->
<!-- test: file=files/rocscience/vp054a.xlsx, type=single_circle, num_slices=50, fs_bishop=1.100, benchmark=VP54-nopile -->
<!-- test: file=files/rocscience/vp054b.xlsx, type=single_circle, num_slices=50, fs_bishop=1.185, benchmark=VP54-pile -->
<!-- test: file=files/rocscience/vp086.xlsx, type=circular_search, num_slices=50, fs_bishop=1.617, fs_spencer=1.611, benchmark=VP86 -->
<!-- test: file=files/rocscience/vp062a.xlsx, type=circular_search, num_slices=50, fs_spencer=1.001, fs_bishop=0.991, benchmark=VP62-dry -->
<!-- test: file=files/rocscience/vp062b.xlsx, type=circular_search, num_slices=50, fs_spencer=1.001, fs_bishop=0.986, benchmark=VP62-ru -->
<!-- test: file=files/rocscience/vp096.xlsx, type=single_circle, rapid=true, num_slices=60, fs_spencer=1.434, fs_bishop=1.432, benchmark=VP96 -->
<!-- test: file=files/rocscience/vp085a.xlsx, type=single_circle, num_slices=60, fs_bishop=1.567, fs_spencer=1.567, benchmark=VP85-active -->
<!-- test: file=files/rocscience/vp085b.xlsx, type=single_circle, num_slices=60, fs_oms=1.319, fs_bishop=1.319, benchmark=VP85-passive -->
<!-- test: file=files/rocscience/vp021a.xlsx, type=single_circle, num_slices=60, fs_oms=1.927, fs_bishop=2.075, fs_spencer=2.071, fs_mprice=2.071, benchmark=VP21-dry -->
<!-- test: file=files/rocscience/vp021b.xlsx, type=single_circle, num_slices=60, fs_oms=1.606, fs_bishop=1.759, fs_spencer=1.757, fs_mprice=1.756, benchmark=VP21-ru -->

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| 1 | Slope, homogenous | covered | [LEM sample 12](lem/samples.md#verification-acads-simple) (`xslope_acads_simple.xlsx`); ACADS table above. Bishop 0.985 vs Slide 0.987. |
| 2 | Slope, homogenous, tension crack | **built** | [vp002.xlsx](files/rocscience/vp002.xlsx). Bishop 1.589 / Spencer 1.585 / Janbu(corr) 1.495 / M-P 1.586 vs Slide 1.596 / 1.592 / 1.489 / 1.592 (±0.4%); Giam reference 1.65. |
| 3 | Slope, (3) materials | **built** | [vp003.xlsx](files/rocscience/vp003.xlsx). Bishop 1.403 / Spencer 1.372 / Janbu(corr) 1.354 / M-P 1.371 vs Slide 1.405 / 1.375 / 1.357 / 1.374 (±0.3%); ACADS reference 1.39. Interface coordinates read from the labeled GeoStudio verification-manual figure of the same ACADS 1(c) problem. |
| 4 | Slope, (3) materials, seismic | **built** | [vp004.xlsx](files/rocscience/vp004.xlsx). Problem #3 + k=0.15g. Bishop 1.013 / Spencer 0.989 / Janbu(corr) 0.963 / M-P 0.987 vs Slide 1.016 / 0.991 / 0.965 / 0.989 (±0.3%); ACADS reference 1.00. |
| 5 | Dam, (4) materials | **built** | [vp005.xlsx](files/rocscience/vp005.xlsx). Talbingo Dam, end of construction (polygon-zone geometry). Critical mechanism is the infinite-slope limit on the upstream face: stored shallow circle gives 1.955 (all methods) vs Slide 1.948-1.949 and the tan φ/tan β limit 1.9475. |
| 6 | Dam, (4) materials, predefined slip surface | **built** | [vp006.xlsx](files/rocscience/vp006.xlsx). Specified circle (100.3, 291, R=278.8) through the inclined core. Bishop 2.206 / Spencer 2.290 / Janbu(corr) 2.073 / M-P 2.299 vs Slide 2.208 / 2.292 / 2.073 / 2.301 (±0.1%); ACADS reference 2.29. This problem exposed (and now guards) the folded-zone slice-weight bug fixed in `slice.py`. |
| 7 | Slope, (2) materials, weak layer | covered | [LEM sample 13](lem/samples.md#verification-acads-weak-layer) (`xslope_acads_weak_layer.xlsx`) is this exact problem (ACADS 3(a)). Non-circular search: Spencer 1.258 / M-P 1.248 vs Slide 1.246 / 1.275; Giam reference 1.24-1.27. |
| 8 | Slope, (2) materials, weak layer, predefined slip surface | **built** | [vp008.xlsx](files/rocscience/vp008.xlsx). Specified 4-point surface (Table 8.2). Spencer 1.276 / Janbu(corr) 1.294 / M-P 1.260 vs Slide 1.277 / 1.294 / 1.262 (exact to ±0.002); SLOPE/W M-P 1.261; Giam reference 1.34. |
| 9 | Slope, (2) materials, weak layer, water table, distributed load | **built** | [vp009.xlsx](files/rocscience/vp009.xlsx). ACADS 4: inclined 0.6 m seam (geometry from the labeled GeoStudio figure), 8-point piezometric line, two surcharge strips. Non-circular search: Spencer 0.724 / Janbu(corr) 0.718 (M-P reaches 0.707 from a wider seed but does not solve the stored seed, so it is not tagged) vs Slide 0.760/0.720/0.734 (block search) and 0.707/0.683/0.699 (optimized); SLOPE/W 0.699-0.689; ACADS references 0.78 [Giam], 0.6878 [Slope 2000], 20-program mean 0.808. Published spread is wide; XSLOPE sits mid-band. |
| 10 | Slope, homogenous, pore pressure grid, ponded water | planned | Tier A via u='seep' per plan (ACADS consensus 1.53), not the pore-pressure grid |
| 11 | Embankment, (2) materials, pore pressure grid | planned |  |
| 12 | Embankment, (4) materials, tension crack, pore pressure grid | planned |  |
| 13 | Embankment, (3) materials, pore pressure grid | planned |  |
| 14 | Slope, homogenous | covered | Arai & Tagyo (1985) ex. 1 — [LEM sample 14](lem/samples.md#verification-arai-tagyo) (`xslope_arai_tagyo.xlsx`); Bishop ref 1.451, Janbu 1.265 |
| 15 | Slope, (3) materials, weak layer | **built** | [vp015.xlsx](files/rocscience/vp015.xlsx). Arai & Tagyo (1985) ex. 2, weak middle band. Circular search: Bishop 0.419 / Spencer 0.422 / Janbu(corr) 0.436 / M-P 0.420 vs Slide 0.420 / 0.409 / 0.423 / (GLE) 0.437; A&T Bishop 0.417; Kim et al. 0.43. |
| 16 | Slope, homogenous, water table | **built** | [vp016.xlsx](files/rocscience/vp016.xlsx). Arai & Tagyo (1985) ex. 3, piezometric line. Circular search: Bishop 1.112 / Spencer 1.113 / Janbu(corr) 1.122 / M-P 1.111 vs Slide 1.118 / 1.118 / 1.131; A&T Bishop 1.138. |
| 17 | Slope, homogenous | **built** | [vp017.xlsx](files/rocscience/vp017.xlsx). Yamagami & Ueta (1988). Circular search: Ordinary 1.274 / Bishop 1.342 / Spencer 1.340 vs Slide 1.278 / 1.344 and Y&U 1.282 / 1.348; published non-circular Spencer 1.325-1.339 (our local non-circular search reaches 1.394 — same search-power note as #19/#20). |
| 18 | Slope, homogenous slope, ru pore pressure | **built** | [vp018.xlsx](files/rocscience/vp018.xlsx). Spencer (1969)/Baker (1980) slope, ru=0.5, non-circular search (right-facing). Spencer 1.033 / M-P 1.024 vs Slide 1.010 (random search + Monte-Carlo optimization), Baker 1.02, Spencer (1969) 1.08. |
| 19 | Slope, (4) materials | **built** | [vp019.xlsx](files/rocscience/vp019.xlsx). Greco (1996) ex. 4 / Yamagami & Ueta (1988) four-layer slope. Circular search: Spencer 1.429 / Bishop 1.448 vs published Spencer 1.40-1.42. Non-circular: XSLOPE's local search plateaus at ~1.45 from the stored seed while Slide's Monte-Carlo optimization reaches 1.398 — a search-power gap (noted for future search work), not a model difference. |
| 20 | Slope, (4) materials, weak layer, water table | **built** | [vp020.xlsx](files/rocscience/vp020.xlsx). Greco (1996) ex. 5 / Chen & Shao (1988): 0.5 m weak seam along the inclined base (polygon zones), water table. Circular: Bishop 1.086 / Spencer 1.091 vs Slide 1.087 / 1.093 (exact). Non-circular seam block: local search 1.082 vs Slide Monte-Carlo 1.010, Chen & Shao 1.01-1.03, Greco 0.973-1.1 — same search-power gap as #19. |
| 21 | Slope, homogenous, ru pore pressure | partial | [vp021a.xlsx](files/rocscience/vp021a.xlsx) (dry), [vp021b.xlsx](files/rocscience/vp021b.xlsx) (ru=0.25) — Fredlund & Krahn (1977) classic, fixed circle (120,90,R=80), imperial units. Dry: OMS 1.927 / Bishop 2.075 / Spencer 2.071 / M-P 2.071 vs F&K 1.928 / 2.080 / 2.073 / 2.076. ru: OMS 1.606 / Bishop 1.759 / Spencer 1.757 / M-P 1.756 vs F&K 1.607 / 1.766 / 1.761 / 1.764 (XSLOPE matches the F&K OMS-ru value exactly; Slide reports 1.687 there). Case 3 (water table) pending the phreatic-line coordinates. |
| 22 | Slope, (2) materials, weak layer, ru pore pressure | planned |  |
| 23 | Slope, (3) materials | **built** | [vp023.xlsx](files/rocscience/vp023.xlsx). Low (1989): undrained layers, lower cu grows 15→30 kPa with depth (`cp` linear-strength option). Circular search: Ordinary 1.357 / Bishop 1.130 vs Low 1.36 / 1.14 (Slide 1.370 / 1.192; Kim 1.17 — the published Bishop values themselves spread 1.14-1.19). |
| 24 | Slope, (3) materials | **built** | [vp024.xlsx](files/rocscience/vp024.xlsx). Low (1989) three-layer undrained slope (φ=0). Circular search: Ordinary 1.433 / Bishop 1.433 vs Slide 1.439 / 1.439; Low reference 1.44. |
| 25 | Bearing capacity test slope, homogenous, distributed load, predefined slip surface | blocked | Prandtl bearing-capacity mechanism (weightless φ=0 soil, surface load, theoretical FS=1.0; Slide Spencer 1.051). XSLOPE's surface-validity checks reject failure surfaces whose two ends sit at equal elevation (the 'flat arc' guard) — flat-ground bearing mechanisms cannot currently be evaluated. Feature gap noted for a relaxed guard when driving comes from loads. |
| 26 | Bearing capacity test prism, homogenous, distributed load, predefined slip surface | blocked | Prandtl bearing-capacity mechanism (weightless φ=0 soil, surface load, theoretical FS=1.0; Slide Spencer 0.940). XSLOPE's surface-validity checks reject failure surfaces whose two ends sit at equal elevation (the 'flat arc' guard) — flat-ground bearing mechanisms cannot currently be evaluated. Feature gap noted for a relaxed guard when driving comes from loads. |
| 27 | Slope, (2) materials, tension crack, water table (auto Hu) | planned |  |
| 28 | Excavated slope and embankment, (3) materials and (5) materials, probabilistic analysis | planned |  |
| 29 | Submerged slope, homogenous, probabilistic analysis, water table | planned |  |
| 30 | Reinforced embankment, (4) materials, tension crack, geosynthetic | planned |  |
| 31 | Reinforced embankment, (5) materials, geosynthetic | planned |  |
| 32 | Reinforced embankment, (7) materials, geosynthetic | planned |  |
| 33 | Dike, (5) materials, probabilistic analysis, water table | planned |  |
| 34 | Dam, (3) materials, probabilistic analysis, water table | planned |  |
| 35 | Dam, (5) materials, probabilistic analysis, reliability index | planned |  |
| 36 | Slope, homogenous, probabilistic analysis, ru pore pressure, reliability index | **built** | [vp036.xlsx](files/rocscience/vp036.xlsx). Li & Lumb (1987) / Hassan & Wolff (1999) reliability benchmark (c′=18±3.6, φ′=30±3, γ=18±0.9, ru=0.2). Deterministic Bishop 1.333 vs H&W 1.334 (Slide 1.340). Taylor-series β_ln on the critical surface 2.263 vs H&W (FOSM) 2.336 and Slide (Monte-Carlo) 2.482 — β estimates legitimately spread by method; xslope does not yet perturb ru (σ=0.02, minor). |
| 37 | Slope, homogenous, distributed load, back analysis of required support force and length | planned |  |
| 38 | Excavated slope, homogenous, finite element groundwater seepage analysis, matric suction | planned |  |
| 39 | Reinforced embankment, (2) materials, tension crack, geosynthetic | planned |  |
| 40 | Slope, homogenous, sensitivity analysis | planned |  |
| 41 | Slope, homogenous, ru pore pressure | **built** | [vp041.xlsx](files/rocscience/vp041.xlsx). Jiang, Baker & Yamagami (2003): power-curve strength τ=1.4·σ′^0.8 with ru=0.3 — exercises the v12 `pow` and `ru` options together. Circular search: Bishop 1.668 / Spencer 1.670 / Janbu(corr) 1.660 vs Slide Bishop 1.656 (non-linear path search), Charles & Soares 1.66, published range 1.56-1.67. |
| 42 | Dam, (3) materials, water table, ponded water, tension crack | partial | Baker & Leshchinsky (2001) safety-map dam. Geometry and core polygon fully labeled (extracted); reservoir level 30 on the right face; Slide Spencer 1.925 circular / 1.877 non-circular, Baker 1.91. Blocked on the internal phreatic-line coordinates through the core (drawn but unlabeled) — needs the B&L paper. |
| 43 | Slope, homogenous, planar surface, RocPlane comparison | planned |  |
| 44 | Slope, homogenous | partial | Baker (2003) ex. 1, linear vs power-curve envelopes. Geometry extracted ((0,0)-(6.43,6)-(20,6)); power-curve case buildable (a=1.107, b=0.86, Spencer 0.960/Baker 0.97), but the Mohr-Coulomb rows (Spencer 1.536) do not reconcile with the extracted property table (c'=0, φ=38 gives ~0.84 on a 43° face) — needs the Baker paper to resolve the intended MC properties before building. |
| 45 | Slope, homogenous | **built** | [vp045a.xlsx](files/rocscience/vp045a.xlsx) (Mohr-Coulomb), [vp045b.xlsx](files/rocscience/vp045b.xlsx) (power curve). Baker (2003) ex. 2: linear vs non-linear envelope on the same 4:1 slope. Spencer: MC 2.801 vs Slide 2.794; power curve 2.649 vs Slide 2.662. (Slide's Janbu values are simplified/uncorrected; ours carry the fo correction and agree once scaled.) |
| 46 | Dam, (2) materials, rapid drawdown, finite element groundwater seepage analysis, ponded water | partial | Baker (1993) three-stage dam (dry / steady-state seep / rapid drawdown). The manual itself calls this a validation problem: permeabilities were estimated by Rocscience and the stage-3 undrained strengths live in discrete .fn6 functions. Stage 1 (dry, Spencer 2.534 / Baker 2.41 / theory 2.5) buildable once Figure 46.1 coordinates are read; stages 2-3 map onto xslope's seep + rapid-drawdown pipeline but need those source functions. |
| 47 | Retaining wall, homogenous, planar failure, line load, shotcrete, soil nails | planned |  |
| 48 | Retaining wall, homogenous, planar failure, line load , soil nails, shotcrete | planned |  |
| 49 | Retaining wall, (2) materials, grouted tiebacks, soldier piles | planned |  |
| 50 | Reinforced slope, (2) materials, predefined slip surface, geosynthetic | **built** | [vp050.xlsx](files/rocscience/vp050.xlsx). SNAILZ reference-manual nail wall: 14 rows with per-row length/tensile/bond values, evaluated on the printed deep wedge (-15.8,0)-(0,-5)-(41.7,25). With Slide's nail defaults (tangent orientation, force factored by FS): Janbu(corr) 1.448 vs SNAILZ 1.46 and Slide 1.417. The capacity envelope reproduces the hand-computed available tension at every crossing (Σ 10.6 kip). The shallow (0,0) surface's kink is not printed — only the deep case is tagged. |
| 51 | Slope, (4) materials, water table, tension crack, seismic | partial | geometry + circle (18.058, 66.744, R=86) read from figures; blocked on Layer 4 properties and water-table coordinates (Zhu 2003 paper) |
| 52 | Slope, (4) materials, water table, tension crack | planned |  |
| 53 | Slope, homogenous, water table, tension crack, planar failure, RocPlane comparison | planned |  |
| 54 | Slope, homogenous, micro piles | **built** | [vp054a.xlsx](files/rocscience/vp054a.xlsx) (no pile), [vp054b.xlsx](files/rocscience/vp054b.xlsx) (with pile). Yamagami (2000): micro-pile row at the crest, 10.7 kN shear per pile at 1 m spacing. On the printed critical circle: no-pile Bishop 1.100 vs Slide 1.102 / Yamagami 1.10; with-pile 1.185 vs Slide 1.193 / Yamagami 1.20 (Slide adds the pile shear un-factored, i.e. active application). A free search with the pile finds 1.113 on a circle exiting upslope of the pile — the published comparison is per-circle. |
| 55 | Slope, homogenous, water table | partial | Pockoski & Duncan (2000) test slope 1 — geometry and properties extracted (ground (-75,100)-(0,100)-(100,150)-(170,150), c=300 psf, φ=30, γ=120 pcf; multi-program table: Spencer 1.30, Bishop 1.29, Janbu 1.15, Lowe 1.32); blocked on the paper's water-table coordinates (figure trace only). |
| 56 | Slope, homogenous, water table, tension crack | planned |  |
| 57 | Slope, (2) materials, water table, tension crack, composite surfaces | planned |  |
| 58 | Retaining wall, (8) materials, water table, grouted tieback | planned |  |
| 59 | Retaining wall, homogenous, water table, grouted tieback | planned |  |
| 60 | Retaining wall, (2) materials, tension crack, distributed load, soil nails | planned |  |
| 61 | Slope, homogenous, composite surfaces | planned |  |
| 62 | Slope, homogenous, ru pore pressure, seismic | **built** | [vp062a.xlsx](files/rocscience/vp062a.xlsx) (dry, kc=0.432), [vp062b.xlsx](files/rocscience/vp062b.xlsx) (ru=0.5, kc=0.132). Loukidis et al. (2003) critical-seismic-coefficient benchmark: FS should be 1.0 at kc. Circular search: Spencer 1.001 / 1.001 and Bishop 0.991 / 0.986 vs Slide 1.001 / 1.001 and 0.991 / 0.987 — exact. |
| 63 | Slope, (3) materials, seismic | planned |  |
| 64 | Embankment, (4) materials, water table, tension crack | planned |  |
| 65 | Embankment, (4) materials, water table, ponded water | planned |  |
| 66 | Embankment, (4) materials, water table, ponded water | planned |  |
| 67 | Embankment, (2) materials | planned |  |
| 68 | Embankment, (3) materials, ponded water | planned |  |
| 69 | Embankment, (2) materials, water table, ponded water | planned |  |
| 70 | Submerged slope, homogenous, water table, ponded water | planned |  |
| 71 | Slope, homogenous, finite element groundwater seepage analysis, water table | planned |  |
| 72 | Embankment dam, (4) materials, finite element groundwater seepage analysis, ponded water | planned |  |
| 73 | Excavated slope, (4) materials, tension crack | planned |  |
| 74 | Embankment, (2) materials | planned |  |
| 75 | Dyke, (4) materials | planned |  |
| 76 | Embankment dam, homogenous, finite element groundwater seepage analysis, ponded water | planned |  |
| 77 | Dam, (2) materials, finite element groundwater seepage analysis, ponded water | planned |  |
| 78 | Slope, homogenous | planned |  |
| 79 | Slope, (2) materials, infinite slope failure | planned |  |
| 80 | Embankment, (6) materials | planned |  |
| 81 | Embankment, (2) materials, infinite slope failure | planned |  |
| 82 | Embankment, (2) materials, water table | planned |  |
| 83 | Embankment, (2) materials | planned |  |
| 84 | Embankment, (2) materials | planned |  |
| 85 | Reinforced slope, homogenous, grouted tieback | **built** | [vp085a.xlsx](files/rocscience/vp085a.xlsx) (active), [vp085b.xlsx](files/rocscience/vp085b.xlsx) (passive). Duncan & Wright (2005) Fig. 6.34: one 9,000 lb/ft horizontal tieback at mid-height of an undrained clay slope — the reinforcement acceptance benchmark from the design plan. On Slide's printed critical circles: active 1.567 vs Slide GLE 1.575 (D&W 1.51); passive Bishop 1.319 vs Slide 1.324 (D&W 1.32). Slide's own method table scatters 1.42-2.05 here (concentrated force strains interslice assumptions), so per-circle comparison is the meaningful one. |
| 86 | Reinforced slope, homogenous, grouted tieback | **built** | [vp086.xlsx](files/rocscience/vp086.xlsx). Duncan & Wright (2005) Fig. 7.28 / STABGM reinforced fill on rock: five 800 lb/ft geogrids. Circular search: Bishop 1.617 / Spencer 1.611 vs Slide 1.629 / 1.620; D&W reference 1.61. |
| 87 | Retaining wall, (3) materials, geotextile | planned |  |
| 88 | Retaining wall, (3) materials, geotextile | planned |  |
| 89 | Retaining wall, (3) materials, geotextile | planned |  |
| 90 | Retaining wall, (3) materials, geotextile | planned |  |
| 91 | Retaining wall, (3) materials, geotextile | planned |  |
| 92 | Retaining wall, (3) materials, geotextile | planned |  |
| 93 | Retaining wall, (3) materials, distributed load, geotextile | planned |  |
| 94 | Retaining wall, (3) materials, geotextile | planned |  |
| 95 | Embankment dam, homogenous, rapid drawdown, water table | partial | USACE EM 1110-2-1902 App. G example with the **Corps 1970 2-stage** method (Slide 1.347, USACE 1.35). XSLOPE implements the Duncan-Wright-Wong 3-stage procedure (see #96, same model file) — the 2-stage variant would be a new option. |
| 96 | Embankment dam, homogenous, rapid drawdown, water table | **built** | [vp096.xlsx](files/rocscience/vp096.xlsx). USACE EM 1110-2-1902 (2003) Appendix G example (Figure G-5: 3:1 / 2.5:1 face, pool 103→24, Kc=1 envelope d=1379 psf ψ=18.2°), specified circle (169.5, 210, R=210), Duncan-Wright-Wong 3-stage: Spencer 1.434 / Bishop 1.432 vs Slide 1.443 and USACE 1.44 (Modified Swedish). First corpus problem through the rapid-drawdown pipeline. |
| 97 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 98 | Embankment dam, (5) materials, rapid drawdown, water table | planned |  |
| 99 | Embankment dam, (3) materials, rapid drawdown, water table | planned |  |
| 100 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 101 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 102 | Embankment dam, homogenous, rapid drawdown | planned |  |
| 103 | Undrained slope, multi-model optimization (MMO) | planned |  |
| 104 | Newmark analysis, seismic analysis, multi-modal optimization (MMO) | planned |  |
| 105 | Anisotropic surface, multi-modal optimization (MMO) | planned |  |
| 106 | Support, Ito & Matsui pile | partial | Cai & Ugai (2000): Bishop with Ito & Matsui limit-pressure pile forces (FS 1.13-1.54 vs spacing). XSLOPE models pile resistance via Vcap/Mcap, not the Ito & Matsui pressure equation, so only the no-pile case (1.13/1.14) is directly reproducible; an Ito & Matsui capacity option would be a new feature. |
| 107 | Retaining walls, gabion walls, supports | planned |  |
| 108 | Retaining walls, gabion walls, supports | planned |  |
| 109 | Retaining walls, gabion walls, weak layers | planned |  |
| 110 | Retaining walls, equivalent fluid pressure | planned |  |
| 111 | Helical anchor | planned |  |
