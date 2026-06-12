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
| Bishop's Simplified | 0.987 | 1.00 | −1.3% |
| Simplified Janbu | 0.992 | 1.00 | −0.8% |
| Corps of Engineers | 0.991 | 1.00 | −0.9% |
| Lowe & Karafiath | 0.987 | 1.00 | −1.3% |
| Spencer | 0.986 | 1.00 | −1.4% |

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
| Spencer | 1.279 | ~1.26 | +1.5% |
| Corps of Engineers | 1.355 | ~1.26 | +7.6% |
| Lowe & Karafiath | 1.268 | ~1.26 | +0.6% |
| Simplified Janbu | 1.278 | ~1.26 | +1.4% |

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
| Simplified Janbu | 1.441 | 1.451 | −0.7% |
| Corps of Engineers | 1.477 | 1.451 | +1.8% |
| Lowe & Karafiath | 1.439 | 1.451 | −0.8% |
| Spencer | 1.402 | 1.451 | −3.4% |

---

## Finite-Element Seepage

### Confined radial flow (analytical anchor — exact)

Full details: [seepage sample problem 8](seep/samples.md#verification-confined-radial).

A quarter-annulus confined flow problem: inner arc (r = 10) at head 30, outer
arc (r = 30) at head 10, straight radial edges as no-flow streamlines. Steady
saturated flow has the exact solution q = k·(π/2)·Δh / ln(r₂/r₁) = 28.596
(k = 1) with a logarithmic head profile.

| Quantity | XSLOPE (tri6) | Exact | Diff |
|---|---|---|---|
| Discharge q | 28.5961 | 28.5960 | <0.01% |
| Max nodal head error | 0.004 | 0 | 0.02% of total drop |

Mesh-converged (identical at 2k and 6k nodes); quad8 gives the same result.
The only error source is faceting of the curved arcs by the polygon boundary.

### Partially penetrating sheetpile (analytical anchor — exact)

Full details: [seepage sample problem 9](seep/samples.md#verification-sheetpile).

Pavlovsky's conformal-mapping solution for a cutoff wall of depth s in a
confined stratum of thickness T (Harr, 1962; Polubarinova-Kochina, 1962):
q = k·H·K(λ′)/(2·K(λ)) with λ = sin(πs/2T). At s/T = ½ the modulus is
self-dual and q = k·H/2 **exactly**. The closed form was additionally verified
with an independent finite-difference solution of the same boundary-value
problem (agreement ~0.4–0.5% at three penetration ratios).

| Case | XSLOPE q | Exact q | Diff | Head below wall tip |
|---|---|---|---|---|
| s/T = 0.50 | 5.010 | 5.000 | +0.20% | 20.0000 (exact: 20) |
| s/T = 0.75 | 3.412 | 3.403 | +0.27% | 20.0000 (exact: 20) |

The error halves with mesh refinement (set by the r^−½ singularity at the wall
tip) and converges to the exact value from above. The head on the wall plane
below the tip equals (h₁+h₂)/2 exactly, an antisymmetry property of the exact
solution that the FE solution reproduces to four decimals.

### Kozeny dam with toe drain (free-surface check)

A homogeneous dam on an impervious base with a horizontal toe drain, built with
the upstream face as the exact confocal-parabola equipotential of the Kozeny
flow net. The computed free surface tracks the analytical parabola, and the
discharge brackets the closed form q = k·y₀ = 4.0 within about ±4% — the
spread reflecting the sensitivity of the discharge to exactly where the
drain-start boundary lands relative to the free-surface/base tangent point, a
known feature of this idealization rather than a solver characteristic. Full
details: [seepage sample problem 7](seep/samples.md#verification-kozeny).

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
Mohr-Coulomb viscoplastic formulation with the displacement-catastrophe failure
criterion (see [FEM Overview](fem/overview.md)).

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
| Full reservoir (free surface) | 1.94 | ~1.9 | +2% |
| Before filling (no free surface) | 2.45 | ~2.4 | +2% |

The wet case is measured at an automatically-selected characteristic point on
the failure mechanism (the criterion picks the node with the fastest plastic-
displacement growth — here the downstream toe) rather than by the global
displacement maximum:
submerged-boundary problems carry a benign boundary-corner creep artifact that
dominates the global maximum and, read naively, biases the apparent failure
point about 10% high. The characteristic-point displacement is flat to
F = 1.8 and jumps tenfold between 1.9 and 2.0. The input model is
independently validated by XSLOPE's Spencer analysis (1.915 vs the paper's
limit-equilibrium 1.90), and the relative reservoir effect matches the paper
(wet/dry = 0.80 vs 0.79).

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
