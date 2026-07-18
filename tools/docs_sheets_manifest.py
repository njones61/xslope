"""Manifest of the input-template docs images: which workbook / sheet / window
each ``docs/usage/images/sheet_*.png`` is rendered from.

Consumed by ``tools/render_sheet.py`` (``python tools/render_sheet.py`` rebuilds
every entry). Sources are either committed sample files (``docs/lem/files`` etc.,
whose sheet layouts are identical in template v12 and v15) or the current-version
files (re)built by ``tools/build_docs_sheet_samples.py`` under
``docs/usage/sample_sheets/`` (mat / rapid / seep-bc, the sheets that changed
shape or need data no committed file carries).

The renderer auto-frames each image to its formatted-content bounding box plus one
spare row/column of margin, and auto-extends the right edge so Excel-style overflow
(the option/help legends) renders complete. So ``rows`` is only needed where the row
window must be pinned (the three ``mat`` views, kept row-aligned), and ``cols`` only
where a horizontal *selection* is required — the multi-table sheets (show N tables,
not all 15/6 slots), the ``seep bc`` split (skip the empty M:S slots), and the sheets
whose real content stops before a block of hidden dropdown-list data (reinforce, piles).

Each entry:
  out            output PNG name in docs/usage/images/
  src            source workbook, repo-relative
  sheet          worksheet name
  rows           (first, last) 1-based inclusive row window, pinned verbatim (optional)
  cols           column selection as letters, e.g. "A:U" or "A:L,T:Z"      (optional)
  identity_cols  columns re-shown at the left of a split view              (optional)
  renderer       "grid" (default), "libreoffice", or "manual" (keep-manual)
  note           free text (used with renderer="manual")
"""

TEMPLATE = "xslope/resources/input_template.xlsx"
MAT = "docs/usage/sample_sheets/sheets_mat.xlsx"
RAPID = "docs/usage/sample_sheets/sheets_rapid.xlsx"
SEEPBC = "docs/usage/sample_sheets/sheets_seepbc.xlsx"

SHEETS = [
    # main comes from the blank template (what that section describes); auto-framed.
    {"out": "sheet_main.png", "src": TEMPLATE, "sheet": "main"},
    # plot is a live Excel chart (a drawing object, not cell data). LibreOffice
    # rasterises it with a clipped axis and a monster legend, so the hand-taken
    # capture is kept and NOT regenerated.
    {"out": "sheet_plot.png", "src": RAPID, "sheet": "plot",
     "renderer": "manual", "note": "live Excel chart; manual capture stays"},

    # mat (v16, 38 cols A:AL) — one wide sheet shown as three views, each matching
    # one of the sheet's own row-9 band headers exactly, so the split needs no
    # hand-picked column break: "Shear Strength/Stiffness" (C:X, which now also
    # carries t_cut/E/nu and so all four option legends: strength options, the
    # color legend, pore-pressure options, and the elastic row), "Standard
    # Deviations" (Y:AD), and "Seepage" (AE:AL, which carries the unsat-model
    # legend). Each view re-shows the mat/name identity columns on the left; rows
    # auto-frame to that view's own content (the row-number gutter keeps the
    # material rows aligned across views by absolute row number, not by shared
    # framing top).
    {"out": "sheet_mat1.png", "src": MAT, "sheet": "mat",
     "cols": "C:X", "identity_cols": "A:B"},
    {"out": "sheet_mat2.png", "src": MAT, "sheet": "mat",
     "cols": "Y:AD", "identity_cols": "A:B"},
    {"out": "sheet_mat3.png", "src": MAT, "sheet": "mat",
     "cols": "AE:AL", "identity_cols": "A:B"},

    # profile / dloads carry many blank table slots; select the filled tables and
    # let the renderer keep their full (bordered) height.
    {"out": "sheet_profile.png", "src": RAPID, "sheet": "profile", "cols": "A:N"},
    {"out": "sheet_piezo.png", "src": RAPID, "sheet": "piezo"},
    {"out": "sheet_circles.png", "src": RAPID, "sheet": "circles"},
    {"out": "sheet_noncirc.png", "src": "docs/lem/files/xslope_noncircular.xlsx",
     "sheet": "non-circ"},
    {"out": "sheet_dloads.png", "src": RAPID, "sheet": "dloads", "cols": "A:P"},
    {"out": "sheet_dloads2.png", "src": RAPID, "sheet": "dloads (2)", "cols": "A:P"},
    # reinforce/piles: stop at the real table edge — the columns beyond hold hidden
    # dropdown-list source data (Z:AB). piles keeps its LEM/FEM colour key (to T).
    {"out": "sheet_reinforce.png", "src": "docs/lem/files/xslope_reinforce.xlsx",
     "sheet": "reinforce", "cols": "A:R"},
    {"out": "sheet_piles.png", "src": "docs/lem/files/xslope_piles.xlsx",
     "sheet": "piles", "cols": "A:T"},
    {"out": "sheet_lloads.png", "src": "docs/inputs/slope/xslope_nail_axial.xlsx",
     "sheet": "lloads"},
    # exit face + three BC blocks (one flux, two head) + the BC-option legend;
    # the two empty BC slots (M:S) are skipped so the legend is not pushed off-image.
    {"out": "sheet_seepbc.png", "src": SEEPBC, "sheet": "seep bc",
     "cols": "A:L,T:Z"},
    {"out": "sheet_seepbc2.png", "src": SEEPBC, "sheet": "seep bc (2)",
     "cols": "A:L,T:Z"},
]
