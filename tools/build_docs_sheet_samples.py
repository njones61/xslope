"""Build the per-sheet SAMPLE workbooks that the input-template docs images are
rendered from (see ``tools/render_sheet.py`` and ``tools/docs_sheets_manifest.py``).

Most template sheets (profile, polygon, circles, reinforce, piles, lloads,
non-circ, dloads) have an identical layout in template v12 and the current v15,
so the manifest points the renderer straight at existing committed sample files
in ``docs/lem/files`` / ``docs/seep/files``. Three sheets changed shape (or need
data no committed file carries) and so are (re)built here as current-version
files under ``docs/usage/sample_sheets/``:

  * ``sheets_mat.xlsx``   — a materials table exercising every strength option
    (mc / cp / pow / hb), plus gamma_sat, the ru pore option, both unsaturated
    laws (lf and van Genuchten), a couple of standard deviations, and stiffness.
    Written cell-by-cell into a fresh template via :func:`mat_header_cols` +
    :func:`write_cells_to_xlsx`, so N/A cells stay *blank* (matching how the
    conditional formatting greys them) rather than the ``0`` that a full
    ``save_slope_data_to_xlsx`` round-trip would leave.
  * ``sheets_rapid.xlsx`` — a rapid-drawdown model (two piezometric lines, two
    distributed-load stages) round-tripped to v15 so the piezo sheet carries the
    v13+ Type row; line 2 is flipped to ``phreatic`` to show both line types.
  * ``sheets_seepbc.xlsx`` — a two-stage seepage model round-tripped to v15, with
    a specified-flux (infiltration) boundary added to stage 1 so the seep-bc image
    shows all three BC kinds (head, flux, exit face).

Run:  python tools/build_docs_sheet_samples.py     # regenerate the three files
"""

from __future__ import annotations

import os

from xslope.fileio import (
    load_slope_data, save_slope_data_to_xlsx, default_template_path,
    mat_header_cols, write_cells_to_xlsx, cell_ref,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(REPO_ROOT, "docs", "usage", "sample_sheets")


# --------------------------------------------------------------------------- #
# 1. mat — write only the applicable cells (blank N/A cells, like a hand-filled
#    template), so the conditional-format greying reads clean in the docs image.
# --------------------------------------------------------------------------- #
# Each material is a {header-key: value} dict; header keys match mat_header_cols()
# (underscore-stripped: 'powa', 'hbsci', 's(g)', ...). Only keys present are
# written; everything else is left blank.
MAT_SHOWCASE = [
    # Embankment fill — Mohr-Coulomb, moist+saturated unit weights, a Kc=1
    # rapid-drawdown envelope (d, psi), and two reliability std deviations.
    {"name": "Embankment", "g": 125, "gsat": 128, "option": "mc",
     "c": 50, "f": 32, "d": 300, "psi": 20, "u": "piezo",
     "s(c)": 15, "s(f)": 3,
     "k1": 5e-5, "k2": 5e-5, "alpha": 0, "unsat": "lf", "kr0": 0.001, "h0": -1,
     "E": 30000, "nu": 0.33},
    # Soft foundation clay — c/p undrained strength increasing below a reference
    # elevation.
    {"name": "Soft Clay", "g": 118, "gsat": 122, "option": "cp",
     "c": 600, "c/p": 12, "r-elev": 95, "u": "piezo",
     "s(c)": 60, "s(c/p)": 2,
     "k1": 1e-7, "k2": 1e-7, "alpha": 0, "unsat": "lf", "kr0": 0.001, "h0": -1,
     "E": 8000, "nu": 0.40},
    # Rockfill shell — nonlinear power-curve envelope.
    {"name": "Rockfill", "g": 140, "option": "pow",
     "powa": 5.5, "powb": 0.82, "powc": 0, "powd": 10, "u": "none",
     "k1": 1e-2, "k2": 1e-2, "alpha": 0, "unsat": "lf", "kr0": 0.001, "h0": -1,
     "E": 50000, "nu": 0.30},
    # Weathered bedrock — generalized Hoek-Brown.
    {"name": "Weathered Rock", "g": 155, "option": "hb",
     "hbsci": 25000, "hbgsi": 45, "hbmi": 8, "hbd": 0, "u": "none",
     "k1": 1e-6, "k2": 1e-6, "alpha": 0, "unsat": "lf", "kr0": 0.001, "h0": -1,
     "E": 500000, "nu": 0.25},
    # Drainage sand — Mohr-Coulomb with the ru pore-pressure option and a van
    # Genuchten unsaturated curve.
    {"name": "Drain Sand", "g": 120, "option": "mc",
     "c": 0, "f": 34, "u": "ru", "ru": 0.15,
     "k1": 5e-2, "k2": 5e-2, "alpha": 0, "unsat": "vg", "a": 0.05, "n": 1.8,
     "E": 40000, "nu": 0.30},
]


def build_mat(out_path):
    template = default_template_path()
    import shutil
    shutil.copy(template, out_path)
    header_row, cols = mat_header_cols(out_path)

    def col_for(key):
        return cols.get(str(key).replace("_", ""))

    updates = {}
    for i, mat in enumerate(MAT_SHOWCASE):
        row = header_row + 1 + i
        updates[cell_ref(row, cols["mat"])] = i + 1
        for key, val in mat.items():
            c = col_for(key)
            if c is None:
                raise KeyError("mat header %r not found in template" % key)
            updates[cell_ref(row, c)] = val
    write_cells_to_xlsx(out_path, {"mat": updates})
    return out_path


# --------------------------------------------------------------------------- #
# 2 & 3. round-trip existing rich files to v15, augmenting where a feature that
#        no committed file carries is needed for the docs.
# --------------------------------------------------------------------------- #
def _fill_material_names(path, sheet, materials):
    """Write the profile/polygon row-6 material NAMES as literal text.

    Those cells hold an XLOOKUP against the mat sheet; a file written by xslope
    (never opened in Excel) has no cached formula result, so the name row renders
    blank. Overwriting with the resolved name is what Excel would display and keeps
    the docs image complete. Mat-ID values sit in row 5 every third column (B, E,
    H, ...); the name is a merged cell in row 6 anchored one column to the LEFT of
    its Mat-ID value (A6:B6 looks up B5, D6:E6 looks up E5, ...), written at col-1."""
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet]
    id_to_name = {i + 1: str(m.get("name", "")) for i, m in enumerate(materials)}
    updates = {}
    for c in range(2, ws.max_column + 1, 3):
        mid = ws.cell(row=5, column=c).value
        if isinstance(mid, (int, float)):
            name = id_to_name.get(int(mid))
            if name:
                updates[cell_ref(6, c - 1)] = name
    wb.close()
    if updates:
        write_cells_to_xlsx(path, {sheet: updates})


def build_rapid(out_path):
    sd = load_slope_data(os.path.join(REPO_ROOT, "docs/lem/files/xslope_gsat_rapid.xlsx"))
    # Show both piezometric-line types: keep line 1 as a piezometric line, mark
    # line 2 (the drawdown line) as a phreatic surface.
    sd["piezo_phreatic2"] = True
    save_slope_data_to_xlsx(sd, out_path)
    _fill_material_names(out_path, "profile", sd.get("materials", []))
    return out_path


def build_polygon(out_path):
    sd = load_slope_data(os.path.join(REPO_ROOT, "docs/seep/files/xslope_levee_poly.xlsx"))
    save_slope_data_to_xlsx(sd, out_path)
    _fill_material_names(out_path, "polygon", sd.get("materials", []))
    return out_path


def build_seepbc(out_path):
    sd = load_slope_data(os.path.join(REPO_ROOT, "docs/lem/files/xslope_johnson_rapid_KEY.xlsx"))
    # Add a specified-flux (infiltration) boundary on the crest of stage 1, so the
    # seep-bc image shows all three BC kinds alongside the two heads + exit face.
    bc = sd.setdefault("seepage_bc", {})
    bc.setdefault("specified_fluxes", []).append(
        {"flux": 1e-6, "coords": [(320.0, 160.0), (360.0, 180.0), (380.0, 180.0)]}
    )
    save_slope_data_to_xlsx(sd, out_path)
    return out_path


BUILDERS = {
    "sheets_mat.xlsx": build_mat,
    "sheets_rapid.xlsx": build_rapid,
    "sheets_polygon.xlsx": build_polygon,
    "sheets_seepbc.xlsx": build_seepbc,
}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, fn in BUILDERS.items():
        out = os.path.join(OUT_DIR, name)
        fn(out)
        print("built", os.path.relpath(out, REPO_ROOT))


if __name__ == "__main__":
    main()
