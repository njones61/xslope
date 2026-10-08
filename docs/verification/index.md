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
prints how many tags, comparisons, values and distinct models these pages
currently hold, broken down by source.

A strength-reduction factor of safety is **proved** once, by the bisection
that cut it, and **checked** afterwards. The two are not the same work: a
bisection is nine solves rediscovering a number nobody has moved, while the pair
of trials it closed on — the highest strength the slope stood at and the lowest
it failed at — is what the factor rests on, and re-solving those two says whether
it still reproduces. Rows whose tag carries that pair are checked in two
solves; the rest are re-bisected. The heaviest rows of all, where one model is
hours of iterations, are marked for the release gate and run there.
`python run_tests.py --gate` is that gate: every row, every factor re-proved.

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

## How the match dots are scored {#how-the-match-dots-are-scored}

Every corpus page's summary table carries a match dot per problem.

| Symbol | Meaning |
|---|---|
| 🟢 | within 3% of the vendor and/or reference figure |
| 🟡 | 3–6% |
| 🔴 | more than 6% |
| <span class="nodata">⊘</span> | one of the [status terms](#status-terms) in place of a dot |

The dot scores the **match quality of what is built**, not how much of a problem is built: a
partly built problem is scored on the stages that are built, and the partial or blocked detail is
in the row text. **Only same-method pairings derive a dot.** Most rows here are strength-reduction
rows, and their pairing is XSLOPE's SSRM against RS2's own SSR column — the same method under two
names; another program's strength-reduction factor (PLAXIS, Z-Soil, GEO FEM, a published FEM/FDM
referee) pairs the same way. On the rows verified with limit equilibrium instead
([#51](rs2.md#p4-vp51), [#60](rs2.md#rs2-60), [#61](rs2.md#rs2-61) cases 1 and 3, [#68](rs2.md#rs2-68)) the method the
source itself names governs, and where the source names no method the fallback is XSLOPE's Spencer
or Morgenstern-Price against the published headline value. A pairing whose two sides are different
methods is reported as information only and never governs a dot; neither does an unconstrained
XSLOPE search against an unconstrained search of the vendor's, since two programs' searches may
settle on different mechanisms, nor a band stitched together from several programs' answers. A
comparison is scored at the source's own precision, so a difference smaller than the source's
printed or figure resolution counts as a match. The problem's *published* answer — the
referee/consensus value, or the source author's own factor — is a reference authority in its own
right whatever engine produced it, as is a closed form, which governs only where XSLOPE is itself
within band of it. Where a row has more than one valid pairing the dot takes the **best of them**;
where a row scores several cases, the worst scored case sets it. These conventions apply to every
summary table on this page.

**How the tables show it.** Every valid pairing carries its difference inline, in parentheses,
computed source-relative, (XSLOPE − source) / source, to one decimal. Where a table gives each
authority a column of its own the difference sits beside the value it is measured against —
`RS2 SSRM 1.33 (−2.0%)`; where a table gives the authority one column and the readings several, it
sits beside each reading instead. So a column carries a percentage exactly when it is a pairing the
dot could rest on, and a column that is **cross-method** for the row shows bare values, because it
is context rather than a pairing. Against a published *range* the entry reads `(inside)` where
XSLOPE falls within it and otherwise carries the difference to the nearer bound. A source author's
numbers fall on either side of that line depending on how the source published them: a single
headline factor for the problem is the published answer and takes a percentage whatever engine
produced it (Low's factor at [#19](rs2.md#rs2-19), Perry's at [#30](rs2.md#rs2-30)), while a per-method table
from the same author is a set of method-specific values, so each entry stays bare (Yamagami &
Ueta's Bishop, Fredlund & Krahn's four methods).

---

## Status terms {#status-terms}

Every row on a corpus page's summary table carries exactly one of these:

- **a match dot** — built, checked, and how close it lands on its referee: 🟢 within 3%, 🟡 3–6%, 🔴 more than 6%;
- *unconfirmed* — a number was found but the run did not settle it; the row says why in one clause;
- *planned* — can be built with the program and the source data as they stand, and is not built yet;
- *blocked* — cannot be built yet; the row names what is missing;
- *not supported* — left out on purpose;
- *no reference value* — the source publishes nothing comparable to the run.

A row with a status term shows <span class="nodata">⊘</span> in place of a dot. A row whose
problem is built on another corpus page reads *covered* and links there.

---

## References

Full bibliographic details for every author-year citation across the corpus
pages — the published benchmark papers, the vendor verification manuals, and the
analytical-anchor sources — are collected on the shared
[**References**](references.md) page, alphabetical by author.


---
