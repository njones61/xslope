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

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| 1 | Slope, homogenous | covered | [LEM sample 12](lem/samples.md#verification-acads-simple) (`xslope_acads_simple.xlsx`); ACADS table above. Bishop 0.985 vs Slide 0.987. |
| 2 | Slope, homogenous, tension crack | **built** | [vp002.xlsx](files/rocscience/vp002.xlsx). Bishop 1.589 / Spencer 1.585 / Janbu(corr) 1.495 / M-P 1.586 vs Slide 1.596 / 1.592 / 1.489 / 1.592 (±0.4%); Giam reference 1.65. |
| 3 | Slope, (3) materials | blocked | layer-interface coordinates are unlabeled in Fig 3.1 — needs ACADS 1(c) source data |
| 4 | Slope, (3) materials, seismic | planned |  |
| 5 | Dam, (4) materials | planned |  |
| 6 | Dam, (4) materials, predefined slip surface | planned |  |
| 7 | Slope, (2) materials, weak layer | planned |  |
| 8 | Slope, (2) materials, weak layer, predefined slip surface | planned |  |
| 9 | Slope, (2) materials, weak layer, water table, distributed load | planned |  |
| 10 | Slope, homogenous, pore pressure grid, ponded water | planned | Tier A via u='seep' per plan (ACADS consensus 1.53), not the pore-pressure grid |
| 11 | Embankment, (2) materials, pore pressure grid | planned |  |
| 12 | Embankment, (4) materials, tension crack, pore pressure grid | planned |  |
| 13 | Embankment, (3) materials, pore pressure grid | planned |  |
| 14 | Slope, homogenous | covered | Arai & Tagyo (1985) ex. 1 — [LEM sample 14](lem/samples.md#verification-arai-tagyo) (`xslope_arai_tagyo.xlsx`); Bishop ref 1.451, Janbu 1.265 |
| 15 | Slope, (3) materials, weak layer | planned |  |
| 16 | Slope, homogenous, water table | planned |  |
| 17 | Slope, homogenous | planned |  |
| 18 | Slope, homogenous slope, ru pore pressure | planned |  |
| 19 | Slope, (4) materials | planned |  |
| 20 | Slope, (4) materials, weak layer, water table | planned |  |
| 21 | Slope, homogenous, ru pore pressure | planned |  |
| 22 | Slope, (2) materials, weak layer, ru pore pressure | planned |  |
| 23 | Slope, (3) materials | planned |  |
| 24 | Slope, (3) materials | planned |  |
| 25 | Bearing capacity test slope, homogenous, distributed load, predefined slip surface | planned |  |
| 26 | Bearing capacity test prism, homogenous, distributed load, predefined slip surface | planned |  |
| 27 | Slope, (2) materials, tension crack, water table (auto Hu) | planned |  |
| 28 | Excavated slope and embankment, (3) materials and (5) materials, probabilistic analysis | planned |  |
| 29 | Submerged slope, homogenous, probabilistic analysis, water table | planned |  |
| 30 | Reinforced embankment, (4) materials, tension crack, geosynthetic | planned |  |
| 31 | Reinforced embankment, (5) materials, geosynthetic | planned |  |
| 32 | Reinforced embankment, (7) materials, geosynthetic | planned |  |
| 33 | Dike, (5) materials, probabilistic analysis, water table | planned |  |
| 34 | Dam, (3) materials, probabilistic analysis, water table | planned |  |
| 35 | Dam, (5) materials, probabilistic analysis, reliability index | planned |  |
| 36 | Slope, homogenous, probabilistic analysis, ru pore pressure, reliability index | planned |  |
| 37 | Slope, homogenous, distributed load, back analysis of required support force and length | planned |  |
| 38 | Excavated slope, homogenous, finite element groundwater seepage analysis, matric suction | planned |  |
| 39 | Reinforced embankment, (2) materials, tension crack, geosynthetic | planned |  |
| 40 | Slope, homogenous, sensitivity analysis | planned |  |
| 41 | Slope, homogenous, ru pore pressure | planned |  |
| 42 | Dam, (3) materials, water table, ponded water, tension crack | planned |  |
| 43 | Slope, homogenous, planar surface, RocPlane comparison | planned |  |
| 44 | Slope, homogenous | planned |  |
| 45 | Slope, homogenous | planned |  |
| 46 | Dam, (2) materials, rapid drawdown, finite element groundwater seepage analysis, ponded | planned |  |
| 47 | Retaining wall, homogenous, planar failure, line load, shotcrete, soil nails | planned |  |
| 48 | Retaining wall, homogenous, planar failure, line load , soil nails, shotcrete | planned |  |
| 49 | Retaining wall, (2) materials, grouted tiebacks, soldier piles | planned |  |
| 50 | Reinforced slope, (2) materials, predefined slip surface, geosynthetic | planned |  |
| 51 | Slope, (4) materials, water table, tension crack, seismic | partial | geometry + circle (18.058, 66.744, R=86) read from figures; blocked on Layer 4 properties and water-table coordinates (Zhu 2003 paper) |
| 52 | Slope, (4) materials, water table, tension crack | planned |  |
| 53 | Slope, homogenous, water table, tension crack, planar failure, RocPlane comparison | planned |  |
| 54 | Slope, homogenous, micro piles | planned | runnable today via the piles sheet (H = 10.7 kN/m, passive; Slide 1.193, Yamagami 1.20) |
| 55 | Slope, homogenous, water table | planned |  |
| 56 | Slope, homogenous, water table, tension crack | planned |  |
| 57 | Slope, (2) materials, water table, tension crack, composite surfaces | planned |  |
| 58 | Retaining wall, (8) materials, water table, grouted tieback | planned |  |
| 59 | Retaining wall, homogenous, water table, grouted tieback | planned |  |
| 60 | Retaining wall, (2) materials, tension crack, distributed load, soil nails | planned |  |
| 61 | Slope, homogenous, composite surfaces | planned |  |
| 62 | Slope, homogenous, ru pore pressure, seismic | planned |  |
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
| 85 | Reinforced slope, homogenous, grouted tieback | planned | active/passive tieback acceptance pair (D&W 1.51/1.32); geometry from Fig 85.1 |
| 86 | Reinforced slope, homogenous, grouted tieback | planned |  |
| 87 | Retaining wall, (3) materials, geotextile | planned |  |
| 88 | Retaining wall, (3) materials, geotextile | planned |  |
| 89 | Retaining wall, (3) materials, geotextile | planned |  |
| 90 | Retaining wall, (3) materials, geotextile | planned |  |
| 91 | Retaining wall, (3) materials, geotextile | planned |  |
| 92 | Retaining wall, (3) materials, geotextile | planned |  |
| 93 | Retaining wall, (3) materials, distributed load, geotextile | planned |  |
| 94 | Retaining wall, (3) materials, geotextile | planned |  |
| 95 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 96 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 97 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 98 | Embankment dam, (5) materials, rapid drawdown, water table | planned |  |
| 99 | Embankment dam, (3) materials, rapid drawdown, water table | planned |  |
| 100 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 101 | Embankment dam, homogenous, rapid drawdown, water table | planned |  |
| 102 | Embankment dam, homogenous, rapid drawdown | planned |  |
| 103 | Undrained slope, multi-model optimization (MMO) | planned |  |
| 104 | Newmark analysis, seismic analysis, multi-modal optimization (MMO) | planned |  |
| 105 | Anisotropic surface, multi-modal optimization (MMO) | planned |  |
| 106 | Support, Ito & Matsui pile | planned | Ito & Matsui pile — direct comparison, no new code |
| 107 | Retaining walls, gabion walls, supports | planned |  |
| 108 | Retaining walls, gabion walls, supports | planned |  |
| 109 | Retaining walls, gabion walls, weak layers | planned |  |
| 110 | Retaining walls, equivalent fluid pressure | planned |  |
| 111 | Helical anchor | planned |  |
