# Export LEM Calcs

**File → Export LEM Calcs…** writes the factor of safety out as a spreadsheet
whose sums and quotient are **live formulas**. Every slice is a row, every sum is
a real `=SUM()` over those rows, and the factor of safety is a real quotient of
two cells — so the calculation can be checked, and re-checked, cell by cell.

The action becomes available once a limit equilibrium analysis has been run. The
method exported is the one the results view is showing, and the default filename
is `<model>_calcs_<method>.xlsx` beside the model file.

---

## The four sheets

**README.** What the workbook is, which method, which model file, the date, the
xslope version, and a link to that method's published derivation. It also
explains how to read the other three sheets, and which values are *stated* — the
converged answers the solver reported — as opposed to computed by formula.

**Slice Plot.** The critical surface with its slices numbered, so a row in the
slice table can be found on the section.

**Slice Table.** One row per slice, left to right along the failure surface, in
the same curated columns and under the same symbols the
[Analysis Report](reports.md) prints, with the model's own units in each header.
The last columns are the method's per-slice contributions to the two sides of the
factor of safety equation, and under each of those is a `=SUM()`.

**FS Calculation.** The governing equation in the notation the method's
documentation page uses, then the arithmetic: the two sums pulled across from the
Slice Table, their quotient, the factor of safety, and a check cell reading the
workbook's own answer against the solver's.

---

## Methods that iterate

An iterative method cannot be written as one pass of arithmetic — the base normal
force on each slice depends on the factor of safety, which depends on the base
normal forces. The converged values are therefore stated, and the sheet then
demonstrates the closure by recomputing, with formulas, the per-slice quantities
that depend on them:

- **Bishop's Simplified Method** rebuilds m<sub>α</sub> and N′ on every slice from
  the stated factor of safety, reforms both moments from that N′, and divides.
  The last cell is the recomputed factor of safety minus the stated one.
- **Spencer's Method** rebuilds Q on every slice from the stated pair (F, θ) and
  sums it twice — once as a force, once as a moment about the origin. Both sums
  must read zero.

Every check cell is labelled in words and says what it should read. None of them
is exactly zero: what a converged iteration leaves behind is a small residual,
and the sheet shows it rather than rounding it away.

The closure block is written where the model allows it. A model carrying passive
support — capacity that mobilises with the soil, and so divides by the factor of
safety on both sides — gets the equation, the sums and the quotient, with the
converged unknowns stated.

---

## Notes

- The workbook holds formulas and no cached results, and asks to be recomputed
  when it is opened, so what you see is always what the formulas say.
- Nothing is written as a cell comment. Every explanation is a cell you can read,
  print and search.
- The two sums and the factor of safety carry workbook names —
  `Sum_resisting`, `Sum_driving` and `Factor_of_safety` — so a formula of your
  own can refer to them by name.
- Headless, the same export is `xslope.calc_export.export_lem_calcs(slope_data,
  bundle, path)`, where `bundle` is a solved LEM solution.
