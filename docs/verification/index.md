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
| [Rocscience RS2 (SSR) Corpus](rs2.md) | The 68-problem RS2 shear-strength-reduction manual (Parts I–III), many sharing geometry with the Slide2 corpus — *stub, planned* |
| [GeoStudio (SLOPE/W) Corpus](geostudio.md) | The 47-problem SLOPE/W verification manual, cross-referenced |
| [PLAXIS LE (SVSLOPE) Corpus](plaxis_le.md) | The ~105-problem SVSLOPE verification manual (third-source cross-bearing) — *stub, planned* |

---

## References

Full bibliographic details for every author-year citation across the corpus
pages — the published benchmark papers, the vendor verification manuals, and the
analytical-anchor sources — are collected on the shared
[**References**](references.md) page, alphabetical by author.


---
