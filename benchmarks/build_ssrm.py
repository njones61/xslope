"""Build the SSRM benchmark input files (SSRM-2: Griffiths & Lane Ex. 6).

Griffiths & Lane (1999), Geotechnique 49(3):387-403, Example 6 (Fig. 16):
two-sided earth embankment on a foundation layer, original cross-section from
Torres & Coffman (1997), homogenized: c' = 13.8 kPa, phi' = 37 deg,
gamma = 18.2 kN/m3 above and below the water table.

Geometry (metres; the figure dimensions are ft conversions and slightly
over-determined — base width 124.4 with 18/23 deg faces doesn't close by 1.4 m.
We anchor on H = 21.3, crest = 7.3 and base = 124.4, splitting the residual by
the angle ratio, giving face runs 66.33 / 50.77 (17.8 / 22.8 deg, within 0.25
deg of the stated angles):

  foundation: y = 0 (fixed base) to 7.3; aprons 33.5 each side
  upstream toe (33.5, 7.3) -> crest (99.83, 28.6)-(107.13, 28.6)
  -> downstream toe (157.9, 7.3); right edge x = 191.4

Full-reservoir case: reservoir level 17.1 above foundation (y = 24.4);
free surface from the upstream face (waterline x = 86.75) to the downstream
toe, used as a piezometric line (G&L take u = gamma_w x vertical depth below
the free surface — identical to xslope's piezo option). Reservoir water load
applied as normal stress on the submerged boundary (G&L Fig. 13) via dloads.
Dry case (before filling): no piezo line, no dload.

Reference FOS (G&L Figs 18-19): full ~1.9, dry ~2.4 (LE companions 1.90/2.42).

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_ssrm.py
"""
from benchmarks._xlsx_writer import (
    new_file, write_cells_to_xlsx,
    main_cells, material_cells, profile_line_cells, piezo_cells, dload_cells,
    circle_cells,
)

TEMPLATE = "docs/inputs/input_template.xlsx"
OUTDIR = "docs/fem/files"

GAMMA_W = 9.81
Y_FND = 7.3                      # foundation surface
Y_CREST = 28.6                   # 7.3 + 21.3
Y_RES = 24.4                     # 7.3 + 17.1 reservoir level
X_TOE_US, X_CREST_L = 33.5, 99.83
X_CREST_R, X_TOE_DS = 107.13, 157.9
X_R = 191.4
# waterline on the upstream face
X_WL = X_TOE_US + (Y_RES - Y_FND) / (Y_CREST - Y_FND) * (X_CREST_L - X_TOE_US)

GROUND = [(0, Y_FND), (X_TOE_US, Y_FND), (X_CREST_L, Y_CREST),
          (X_CREST_R, Y_CREST), (X_TOE_DS, Y_FND), (X_R, Y_FND)]


def _build(dst, wet):
    new_file(dst, TEMPLATE)
    u = {}
    u['main'] = main_cells(gamma_w=GAMMA_W)
    u['mat'] = material_cells(9, 1, "Embankment", 18.2, "mc", c=13.8, phi=37.0,
                              u=("piezo" if wet else "none"),
                              E=1.0e5, nu=0.3)
    prof = {'B2': 0}              # base of foundation at y = 0
    prof.update(profile_line_cells(1, 1, GROUND))
    u['profile'] = prof
    # placeholder circle (required by the loader; FEM/SSRM ignores it) — set
    # plausibly over the downstream face, tangent to the foundation surface
    u['circles'] = circle_cells(1, 130, 60, option="Depth", depth=Y_FND)
    if wet:
        u['piezo'] = piezo_cells([(0, Y_RES), (round(X_WL, 3), Y_RES),
                                  (X_TOE_DS, Y_FND), (X_R, Y_FND)])
        p = (Y_RES - Y_FND) * GAMMA_W     # 167.75 kPa at foundation level
        u['dloads'] = dload_cells(1, [
            (0, Y_FND, round(p, 3)),
            (X_TOE_US, Y_FND, round(p, 3)),
            (round(X_WL, 3), Y_RES, 0.0),
        ])
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths6_seep():
    """Full-reservoir case with pore pressures from an actual FE seepage
    solution instead of the piezometric-line shortcut. G&L's simplified
    u = gamma_w x (vertical depth below the free surface) is statically
    inconsistent with the elastic field at the upstream boundary and sustains
    a permanent corner creep; the coupled seepage field is consistent.
    Seepage BCs: reservoir head on the submerged upstream boundary, exit face
    on the downstream face and apron."""
    from benchmarks._xlsx_writer import seep_bc_cells
    dst = f"{OUTDIR}/xslope_griffiths6_seep.xlsx"
    new_file(dst, TEMPLATE)
    u = {}
    u['main'] = main_cells(gamma_w=GAMMA_W)
    u['mat'] = material_cells(9, 1, "Embankment", 18.2, "mc", c=13.8, phi=37.0,
                              u="seep", k1=1.0, k2=1.0, alpha=0.0,
                              kr0=0.01, h0=-1.0, E=1.0e5, nu=0.3)
    prof = {'B2': 0}
    prof.update(profile_line_cells(1, 1, GROUND))
    u['profile'] = prof
    u['circles'] = circle_cells(1, 130, 60, option="Depth", depth=Y_FND)
    u['seep bc'] = seep_bc_cells(
        head1=Y_RES,
        head1_pts=[(0, Y_FND), (X_TOE_US, Y_FND), (round(X_WL, 3), Y_RES)],
        exit_face=[(X_CREST_R, Y_CREST), (X_TOE_DS, Y_FND), (X_R, Y_FND)],
    )
    p_res = (Y_RES - Y_FND) * GAMMA_W
    u['dloads'] = dload_cells(1, [
        (0, Y_FND, round(p_res, 3)),
        (X_TOE_US, Y_FND, round(p_res, 3)),
        (round(X_WL, 3), Y_RES, 0.0),
    ])
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths6_full():
    return _build(f"{OUTDIR}/xslope_griffiths6_full.xlsx", wet=True)


def build_griffiths6_dry():
    return _build(f"{OUTDIR}/xslope_griffiths6_dry.xlsx", wet=False)


if __name__ == "__main__":
    build_griffiths6_full()
    build_griffiths6_dry()
    build_griffiths6_seep()
    print(f"\nwaterline on upstream face: x = {X_WL:.3f}")
    print("SSRM benchmark input files built.")
