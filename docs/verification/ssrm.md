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

### Griffiths & Lane (1999) Example 1 — homogeneous slope

Full details: [FEM sample problem 4](../fem/samples.md#verification-griffiths1).

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

Full details: [FEM sample problem 5](../fem/samples.md#verification-griffiths6).

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
