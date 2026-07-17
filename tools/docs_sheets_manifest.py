"""Manifest of the input-template docs images: which workbook / sheet / window
each ``docs/usage/images/sheet_*.png`` is rendered from.

Consumed by ``tools/render_sheet.py`` (``python tools/render_sheet.py`` rebuilds
every entry). Sources are either committed sample files (``docs/lem/files`` etc.,
whose sheet layouts are identical in template v12 and v15) or the current-version
files (re)built by ``tools/build_docs_sheet_samples.py`` under
``docs/usage/sample_sheets/`` (mat / rapid / seep-bc, the sheets that changed
shape or need data no committed file carries).

Each entry:
  out            output PNG name in docs/usage/images/
  src            source workbook, repo-relative
  sheet          worksheet name
  rows           (first, last) 1-based inclusive row window        (optional)
  cols           column window as letters, e.g. "A:U" or "A:R,Z:AB" (optional)
  identity_cols  columns re-shown at the left of a split view       (optional)
  renderer       "grid" (default) or "libreoffice" (the live-chart plot sheet)
"""

TEMPLATE = "xslope/resources/input_template.xlsx"
MAT = "docs/usage/sample_sheets/sheets_mat.xlsx"
RAPID = "docs/usage/sample_sheets/sheets_rapid.xlsx"
POLYGON = "docs/usage/sample_sheets/sheets_polygon.xlsx"
SEEPBC = "docs/usage/sample_sheets/sheets_seepbc.xlsx"

SHEETS = [
    # main + plot come from the blank template (what those sections describe);
    # plot is a live Excel chart, rendered via LibreOffice.
    {"out": "sheet_main.png", "src": TEMPLATE, "sheet": "main",
     "rows": (1, 24), "cols": "A:P"},
    {"out": "sheet_plot.png", "src": RAPID, "sheet": "plot",
     "renderer": "libreoffice"},

    # mat — one wide sheet shown as three views, each re-showing the mat/name
    # identity columns: strength, variability, seepage+stiffness.
    {"out": "sheet_mat1.png", "src": MAT, "sheet": "mat",
     "rows": (1, 16), "cols": "C:U", "identity_cols": "A:B"},
    {"out": "sheet_mat2.png", "src": MAT, "sheet": "mat",
     "rows": (1, 16), "cols": "V:AA", "identity_cols": "A:B"},
    {"out": "sheet_mat3.png", "src": MAT, "sheet": "mat",
     "rows": (1, 16), "cols": "AB:AK", "identity_cols": "A:B"},

    # profile / polygon: crop to the four filled tables (the template carries 15
    # blank table slots that would otherwise make the image mostly empty).
    {"out": "sheet_profile.png", "src": RAPID, "sheet": "profile",
     "rows": (1, 12), "cols": "A:K"},
    {"out": "sheet_polygon.png", "src": POLYGON,
     "sheet": "polygon", "rows": (1, 15), "cols": "A:K"},
    {"out": "sheet_piezo.png", "src": RAPID, "sheet": "piezo",
     "rows": (1, 14)},
    {"out": "sheet_circles.png", "src": RAPID, "sheet": "circles",
     "rows": (1, 14), "cols": "A:K"},
    {"out": "sheet_noncirc.png", "src": "docs/lem/files/xslope_noncircular.xlsx",
     "sheet": "non-circ", "rows": (1, 14)},
    {"out": "sheet_dloads.png", "src": RAPID, "sheet": "dloads",
     "rows": (1, 14), "cols": "A:H"},
    {"out": "sheet_dloads2.png", "src": RAPID, "sheet": "dloads (2)",
     "rows": (1, 14), "cols": "A:H"},
    {"out": "sheet_reinforce.png", "src": "docs/lem/files/xslope_reinforce.xlsx",
     "sheet": "reinforce", "rows": (1, 12), "cols": "A:R"},
    {"out": "sheet_piles.png", "src": "docs/lem/files/xslope_piles.xlsx",
     "sheet": "piles", "rows": (1, 9), "cols": "A:Q"},
    {"out": "sheet_lloads.png", "src": "docs/inputs/slope/xslope_nail_axial.xlsx",
     "sheet": "lloads", "rows": (1, 9)},
    # exit face + three BC blocks (one flux, two head) + the BC-option legend;
    # the two empty BC slots (M:S) are skipped so the legend is not pushed off-image.
    {"out": "sheet_seepbc.png", "src": SEEPBC, "sheet": "seep bc",
     "rows": (1, 12), "cols": "A:L,T:Z"},
    {"out": "sheet_seepbc2.png", "src": SEEPBC, "sheet": "seep bc (2)",
     "rows": (1, 12), "cols": "A:L,T:Z"},
]
