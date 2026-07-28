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

# The editable master the docs page links for download (docs/inputs/input_template.xlsx),
# NOT the bundled package copy (xslope/resources/): the two can lag between a master edit
# and the next resource sync, and the main-sheet image must match the template a reader
# actually downloads and the version the docs prose describes. Same sourcing rationale as
# build_mat() in tools/build_docs_sheet_samples.py.
TEMPLATE = "docs/inputs/input_template.xlsx"
MAT = "docs/usage/sample_sheets/sheets_mat.xlsx"
RAPID = "docs/usage/sample_sheets/sheets_rapid.xlsx"
SEEPBC = "docs/usage/sample_sheets/sheets_seepbc.xlsx"
TSEEP = "docs/usage/sample_sheets/sheets_tseep.xlsx"
POLY = "docs/usage/sample_sheets/sheets_polygon.xlsx"

SHEETS = [
    # main comes from the blank template (what that section describes); auto-framed.
    {"out": "sheet_main.png", "src": TEMPLATE, "sheet": "main"},
    # plot is a live Excel chart (a drawing object, not cell data). LibreOffice
    # rasterises it with a clipped axis and a monster legend, so the hand-taken
    # capture is kept and NOT regenerated.
    {"out": "sheet_plot.png", "src": RAPID, "sheet": "plot",
     "renderer": "manual", "note": "live Excel chart; manual capture stays"},

    # mat (v20, 42 cols A:AP) — one wide sheet shown as three views, each matching
    # one of the sheet's own row-9 band headers exactly, so the split needs no
    # hand-picked column break: "Shear Strength/Stiffness" (C:Z, which carries the
    # v17 matric-suction pair phi_b/s_cap alongside t_cut/E/nu, and so all four
    # option legends: strength options, the color legend, pore-pressure options, and
    # the elastic row), "Standard Deviations" (AA:AF), and "Seepage" (AG:AP, which
    # carries the unsat-model legend plus the v18 transient-storage pair Ss/Sy at
    # AO/AP). Each view re-shows the mat/name identity columns on the left; rows
    # auto-frame to that view's own content (the row-number gutter keeps the material
    # rows aligned across views by absolute row number, not by shared framing top).
    #
    # These windows MUST track the row-9 merges: v19's ssr_zone insert at column O
    # pushed the last two bands one column right, and v20's DELETION of that column
    # pulled them back, so both moves had to be tracked here — a stale window makes
    # the renderer look up a merge anchor outside its own frame (KeyError) rather
    # than quietly cropping.
    {"out": "sheet_mat1.png", "src": MAT, "sheet": "mat",
     "cols": "C:Z", "identity_cols": "A:B"},
    {"out": "sheet_mat2.png", "src": MAT, "sheet": "mat",
     "cols": "AA:AF", "identity_cols": "A:B"},
    {"out": "sheet_mat3.png", "src": MAT, "sheet": "mat",
     "cols": "AG:AP", "identity_cols": "A:B"},

    # profile / polygon / dloads carry many blank table slots; select the filled
    # tables and let the renderer keep their full (bordered) height.
    {"out": "sheet_profile.png", "src": RAPID, "sheet": "profile", "cols": "A:N"},
    # polygon: the levee's four material zones followed by the three v20 SSR zone
    # rows (Mat ID -1/-2/-3), so the image shows the sentinels and the display codes
    # the row-6 formula echoes for them next to ordinary material rows.
    {"out": "sheet_polygon.png", "src": POLY, "sheet": "polygon", "cols": "A:U"},
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
    # tseep (v18) — the transient-seepage time-series table (time axis in column B,
    # named series across C..) plus the run-controls block in columns I/J. Rendered
    # from a filled sample so the layout and an example drawdown series both read. The
    # template pre-borders an empty fill-in grid down to row 100, so the rows are pinned
    # to the populated region (like the mat views) instead of auto-framing to that grid.
    {"out": "sheet_tseep.png", "src": TSEEP, "sheet": "tseep",
     "rows": (1, 13), "cols": "A:K"},
]
