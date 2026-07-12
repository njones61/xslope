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

## Benchmark classes

| Page | Scope |
|---|---|
| [Limit Equilibrium](lem.md) | Analytical anchors and published LEM cross-checks |
| [FE Seepage](seep.md) | Confined/unconfined flow benchmarks and analytical anchors |
| [FE Slope Stability (SSRM)](ssrm.md) | Strength-reduction benchmarks |
| [Rocscience Slide2 Corpus](rocscience.md) | The 111-problem Slide2 verification manual, problem by problem |
| [Rocscience Groundwater Corpus](rocscience_groundwater.md) | The 21-problem Slide2 groundwater (FE seepage) verification manual — *stub, planned* |
| [GeoStudio (SLOPE/W) Corpus](geostudio.md) | The 47-problem SLOPE/W verification manual, cross-referenced |
| [PLAXIS LE (SVSLOPE) Corpus](plaxis_le.md) | The ~105-problem SVSLOPE verification manual (third-source cross-bearing) — *stub, planned* |

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
