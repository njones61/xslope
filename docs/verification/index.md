---
title: "Verification — 800+ published benchmark cases — XSLOPE"
description: "XSLOPE verification: more than 800 benchmark cases checked problem by problem against the published verification manuals of Rocscience Slide2 and RS2 and GeoStudio SLOPE/W, plus analytical solutions — every comparison published."
---

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
below can be regenerated, and the benchmarks are re-verified automatically
whenever XSLOPE changes.

Every comparison the suite runs is declared by a test tag on the page that
publishes it, so the corpus counts itself: `python3 tools/count_verification_cases.py`
prints how many tags, locked comparisons, locked values and distinct models
these pages currently hold, broken down by source.

---

## Benchmark classes

| Page | Scope |
|---|---|
| [FE Seepage](seep.md) | Confined/unconfined flow benchmarks and analytical anchors |
| [FE Slope Stability (SSRM)](ssrm.md) | Strength-reduction benchmarks |
| [Rocscience Slide2 Corpus](rocscience.md) | The 111-problem Slide2 verification manual, problem by problem |
| [Rocscience Groundwater Corpus](rocscience_groundwater.md) | The 21-problem Slide2 groundwater (FE seepage) verification manual, problem by problem |
| [Rocscience RS2 (SSRM) Corpus](rs2.md) | The RS2 shear-strength-reduction manual, Parts I–IV — the corpus's FEM/SSRM backbone |
| [GeoStudio (SLOPE/W) Corpus](geostudio.md) | The 47-problem SLOPE/W verification manual, cross-referenced |
| [Published Problems](published.md) | Worked hand calculations from design manuals and the literature, table by table |

---

## References

Full bibliographic details for every author-year citation across the corpus
pages — the published benchmark papers, the vendor verification manuals, and the
analytical-anchor sources — are collected on the shared
[**References**](references.md) page, alphabetical by author.


---
