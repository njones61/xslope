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

Every row on a corpus page's summary table carries a match dot or one of the
[status terms](#status-terms) below.

**Each row has one referee.** Where a closed-form solution exists for the stated inputs,
it is the referee; where that closed form prices a mechanism the slope does not take, the
searched limit-equilibrium answer is the referee instead. Where no closed form exists, the
referee is the reference the source names: the vendor's own result, the published referee
or consensus value, or the source author's headline factor of safety. The dot is scored
against that referee alone. An experiment (a field failure, a tilt table, a laboratory
column) is a validation datum shown beside the referee, never the referee itself.

The dot scores the match on what is built, not how much of the problem is built: a partly built problem
is scored on the cases that are built, and the row text names what remains. Where a row
carries several cases, the worst case sets the dot. A comparison is scored at the source's
own precision, so a difference smaller than the source's printed or figure resolution
counts as a match.

**How the tables show it.** A comparison carries its difference inline, in parentheses,
computed source-relative, (XSLOPE − source) / source, to one decimal. Where a table gives
each authority a column of its own, the difference sits beside the value it is measured
against — `RS2 SSRM 1.33 (−2.0%)`; where a table gives the authority one column and the
readings several, it sits beside each reading instead. Apart from the referee, a value
from a method other than XSLOPE's on that row is context and stays bare. Against a
published *range* the entry reads `(inside)` where XSLOPE falls within it and otherwise
carries the difference to the nearer bound. A source author's single headline factor for the problem takes a percentage
whatever engine produced it (Low's factor at [#19](rs2.md#rs2-19), Perry's at
[#30](rs2.md#rs2-30)), while a per-method table from the same author is a set of
method-specific values, so each entry stays bare (Yamagami & Ueta's Bishop, Fredlund &
Krahn's four methods).

---

## Status terms {#status-terms}

Every row on a corpus page's summary table carries exactly one of these:

- **a match dot** — built, checked, and how close it lands on its referee: 🟢 within 3%, 🟡 3–6%, 🔴 more than 6%;
- *unconfirmed* — a number was found but the run did not settle it; the row says why in one clause;
- *planned* — can be built with the program and the source data as they stand, and is not built yet;
- *blocked* — cannot be built yet; the row names what is missing;
- *not supported* — left out on purpose;
- *no reference value* — the source publishes nothing to compare against.

A row with a status term shows <span class="nodata">⊘</span> in place of a dot. A row whose
problem is built on another corpus page reads *covered* and links there.

---

## References

Full bibliographic details for every author-year citation across the corpus
pages — the published benchmark papers, the vendor verification manuals, and the
analytical-anchor sources — are collected on the shared
[**References**](references.md) page, alphabetical by author.


---
