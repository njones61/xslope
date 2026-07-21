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
    circle_cells, polygon_cells,
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
    u['mat'] = material_cells(1, "Embankment", 18.2, "mc", c=13.8, phi=37.0,
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
    u['mat'] = material_cells(1, "Embankment", 18.2, "mc", c=13.8, phi=37.0,
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


def build_griffiths2():
    """Griffiths & Lane (1999) Example 2 (Fig. 5): the Example-1 homogeneous
    2:1 slope WITH a foundation layer of the SAME material, thickness H/2 (so
    the firm base sits at depth D = 1.5 H below the crest). English units, all
    Example-1 conventions carried over: H = 50 ft, gamma = 125 pcf, c' = 312.5
    psf (c'/gamma-H = 0.05), phi' = 20 deg, E = 700000 psf, nu = 0.3, gamma_w =
    62.4 pcf. quad8 SSRM at target_size 3.5 (SSRM-1 settings).

    Domain uses the PRINTED Fig. 1 dimensions — 1.2 H crest platform (60 ft),
    2:1 face (toe at x = 160), 2 H runout past the toe (to x = 260) — with the
    H/2 = 25 ft foundation from the Example-2 text (D = 1.5). Example 1 as built
    here omitted the runout because with the base at toe level (D = 1) it is a
    y = 0 sliver that cannot affect a toe mechanism; with a foundation below the
    toe the runout gives any deeper mechanism lateral room, so it is restored.

    Published (paper text, p. 393-394): FE FOS ~= 1.4, unchanged from Example 1
    because the mechanism stays a toe failure; Bishop & Morgenstern (1960)
    base-circle charts give 1.752 and a proprietary slip-circle program forced
    tangent to the base of the foundation returned a false 1.7, while Cousins
    (1978) toe-circle charts agree with 1.4."""
    dst = f"{OUTDIR}/xslope_griffiths2.xlsx"
    new_file(dst, TEMPLATE)
    H = 50.0
    ground = [(0.0, H), (60.0, H), (160.0, 0.0), (260.0, 0.0)]
    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    u['mat'] = material_cells(1, "soil", 125.0, "mc", c=312.5, phi=20.0,
                              u="none", E=7.0e5, nu=0.3)
    prof = {'B2': -H / 2.0}          # firm base at y = -25 (H/2 foundation)
    prof.update(profile_line_cells(1, 1, ground))
    u['profile'] = prof
    # placeholder circle (loader requires one; FEM/SSRM ignores it) — a toe
    # circle over the slope face, tangent into the foundation
    u['circles'] = circle_cells(1, 150.0, 70.0, option="Depth", depth=-H / 2.0)
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths4(ratio, tag):
    """Griffiths & Lane (1999) Example 4 (Fig. 9): an undrained (phi_u = 0) clay
    slope over a WEAK foundation layer. Two clean, uniform-Su layers — the slope
    body (cu1) and the foundation (cu2) — with the strength ratio cu2/cu1 varied.

    Geometry is Fig. 9's, English units carried from Examples 1-2 (H = 50 ft,
    gamma = 125 pcf): D = 2 (firm base 2H = 100 ft below the crest), a 2H = 100 ft
    crest platform, a 2:1 face dropping H = 50 ft to the toe, and a 2H = 100 ft
    runout past the toe. The foundation layer is H = 50 ft thick; its top (y = 50)
    is the toe elevation and the cu1/cu2 material boundary, which is also the ground
    surface of the runout. The section is tiled by two explicit polygons (the
    boundary coincides with the runout surface, so a profile-line stack would leave
    a zero-thickness cu1 sliver over the runout):

      cu1 slope body:  (0,100)-(100,100)-(200,50)-(0,50)
      cu2 foundation:  (0,50)-(300,50)-(300,0)-(0,0)

    Undrained strength: cu1/gamma-H = 0.25 -> cu1 = 0.25 x 125 x 50 = 1562.5 psf,
    phi = 0 (option 'mc', u = 'none' — the established uniform-Su pattern). The
    foundation strength is cu2 = ratio x cu1. E and nu are the elastic_props.py
    soil-type values (imperial): cu1 = 1562.5 psf classifies Medium Clay
    (E = 668300 psf, nu = 0.40); cu2 at ratio 2 = 3125 psf classifies Stiff Clay
    (E = 2610700 psf, nu = 0.30). SSRM FS is E-invariant; the elastic constants
    only set the displacement scale.

    Published anchors (Taylor 1937, quoted by G&L Fig. 10): FOS = 1.47 at
    cu2/cu1 = 1 (a deep BASE circle) and FOS = 2.10 at cu2 >> cu1 (a shallow TOE
    circle). G&L's FE curve flattens to the toe-circle plateau at cu2/cu1 ~ 1.5;
    the mechanism flips from base (Fig. 11a) to toe (Fig. 11c) across that
    transition. The two bracket ratios (1 and 2) sit firmly on the base-circle and
    toe-circle plateaus respectively (Fig. 10 confirmed)."""
    dst = f"{OUTDIR}/xslope_griffiths4_{tag}.xlsx"
    new_file(dst, TEMPLATE)
    cu1 = 1562.5
    cu2 = ratio * cu1
    # elastic_props.py soil-type classification (imperial), per material strength
    from benchmarks.rocscience.elastic_props import classify
    _, E1, nu1 = classify({'option': 'mc', 'c': cu1, 'phi': 0.0}, imperial=True)
    _, E2, nu2 = classify({'option': 'mc', 'c': cu2, 'phi': 0.0}, imperial=True)
    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    mat = material_cells(1, "Slope clay", 125.0, "mc", c=cu1, phi=0.0,
                         u="none", E=E1, nu=nu1)
    mat.update(material_cells(2, "Foundation clay", 125.0, "mc", c=cu2, phi=0.0,
                              u="none", E=E2, nu=nu2))
    u['mat'] = mat
    poly = polygon_cells(1, 1, [(0.0, 100.0), (100.0, 100.0), (200.0, 50.0),
                                (0.0, 50.0), (0.0, 100.0)])
    poly.update(polygon_cells(2, 2, [(0.0, 50.0), (300.0, 50.0), (300.0, 0.0),
                                     (0.0, 0.0), (0.0, 50.0)]))
    u['polygon'] = poly
    # placeholder circle (loader requires one; FEM/SSRM ignores it, and the LEM
    # companion runs seed='grid'). A base circle over the slope, tangent to y = 0.
    u['circles'] = circle_cells(1, 150.0, 90.0, option="Depth", depth=0.0)
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths4_r1():
    return build_griffiths4(1.0, "r1")


def build_griffiths4_r2():
    return build_griffiths4(2.0, "r2")


def build_griffiths6_full():
    return _build(f"{OUTDIR}/xslope_griffiths6_full.xlsx", wet=True)


def build_griffiths6_dry():
    return _build(f"{OUTDIR}/xslope_griffiths6_dry.xlsx", wet=False)


if __name__ == "__main__":
    build_griffiths2()
    build_griffiths4_r1()
    build_griffiths4_r2()
    build_griffiths6_full()
    build_griffiths6_dry()
    build_griffiths6_seep()
    print(f"\nwaterline on upstream face: x = {X_WL:.3f}")
    print("SSRM benchmark input files built.")
