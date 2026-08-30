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
    circle_cells, polygon_cells, noncirc_cells,
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
        # The reservoir is stated ONCE, as the piezometric line, and main!D23
        # (written 'auto' by main_cells) hands the engine the job of turning it
        # into the surface load. The file used to carry both -- the line and a
        # hand-computed hydrostatic block from (0, Y_FND) at 167.75 kPa to the
        # waterline at (X_WL, Y_RES) -- which is the same reservoir written twice
        # and had to be kept in step by hand.
        u['piezo'] = piezo_cells([(0, Y_RES), (round(X_WL, 3), Y_RES),
                                  (X_TOE_DS, Y_FND), (X_R, Y_FND)])
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
    # The reservoir load comes from the same head boundary that supplies the
    # seepage field (main!D23 = auto), so the two can never describe different
    # pools. It used to be a hand-computed hydrostatic block alongside them.
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths1():
    """Griffiths & Lane (1999) Example 1 (Fig. 2): the homogeneous 2:1 slope, the
    root case the later examples build on. English units carried through the whole
    G&L set: H = 50 ft, gamma = 125 pcf, c' = 312.5 psf (c'/gamma-H = 0.05),
    phi' = 20 deg, gamma_w = 62.4 pcf. The firm base sits at the toe elevation
    (D = 1), so no runout is modelled (a y = 0 sliver cannot affect a toe
    mechanism). quad8 SSRM at target_size 3.5 (SSRM-1 settings).

    Geometry (Fig. 1 dimensions): crest platform (0,50)-(60,50), 2:1 face to the
    toe (160,0), firm base at y = 0. The section is a single material described by
    a profile line; get_material_polygons derives the domain from it and the max
    depth. Placeholder starting circle over the slope face (the loader requires one;
    FEM/SSRM ignores it).

    Elastic constants are the paper's printed nominal values (Griffiths & Lane 1999,
    p.390 -- "in the absence of meaningful data for E' and nu', they can be given
    nominal values (e.g. E' = 10^5 kN/m^2 and nu' = 0.3)"), which the paper applies
    throughout. Converted to English units: E' = 1e5 kPa = 2.0885e6 psf, nu' = 0.3.
    SSRM FS is E-invariant, so this does not move the locked factor of safety; it
    sets the displacement scale to the paper's basis.

    Published (Griffiths & Lane, Table 1): FE FOS ~= 1.4 at c'/gamma-H = 0.05,
    phi' = 20 deg; Bishop & Morgenstern (1960) chart gives 1.380. XSLOPE's strict
    true-equilibrium quad8 lock reads FS = 1.35 with the displacement-vs-F upturn at
    F ~= 1.40, bracketing both.
    """
    dst = f"{OUTDIR}/xslope_griffiths1.xlsx"
    new_file(dst, TEMPLATE)
    ground = [(0.0, 50.0), (60.0, 50.0), (160.0, 0.0)]
    # Elastic constants: the paper's printed nominal values (Griffiths & Lane 1999,
    # p.390 -- E' = 10^5 kN/m^2, nu' = 0.3), converted to English units:
    # E' = 1e5 kPa = 2.0885e6 psf. SSRM FS is E-invariant; this only sets the
    # displacement scale to the paper's basis.
    E = 2.0885e6
    nu = 0.3
    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    u['mat'] = material_cells(1, "soil", 125.0, "mc", c=312.5, phi=20.0,
                              u="none", E=E, nu=nu)
    prof = {'B2': 0.0}              # firm base at y = 0 (D = 1, base at toe level)
    prof.update(profile_line_cells(1, 1, ground))
    u['profile'] = prof
    # placeholder circle (loader requires one; FEM/SSRM ignores it)
    u['circles'] = circle_cells(1, 120.0, 80.0, option="Depth", depth=0.0)
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths2():
    """Griffiths & Lane (1999) Example 2 (Fig. 5): the Example-1 homogeneous
    2:1 slope WITH a foundation layer of the SAME material, thickness H/2 (so
    the firm base sits at depth D = 1.5 H below the crest). English units, all
    Example-1 conventions carried over: H = 50 ft, gamma = 125 pcf, c' = 312.5
    psf (c'/gamma-H = 0.05), phi' = 20 deg, gamma_w = 62.4 pcf. Elastic constants
    are the paper's printed nominal values (Griffiths & Lane 1999, p.390: E' = 10^5
    kN/m^2 = 2.0885e6 psf, nu' = 0.3), same as Example 1. quad8 SSRM at target_size
    3.5 (SSRM-1 settings).

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
    # Elastic constants: the paper's printed nominal values (Griffiths & Lane 1999,
    # p.390 — E' = 10^5 kN/m^2, nu' = 0.3), converted to English units:
    # E' = 1e5 kPa = 2.0885e6 psf. SSRM FS is E-invariant; this only sets the
    # displacement scale to the paper's basis.
    u['mat'] = material_cells(1, "soil", 125.0, "mc", c=312.5, phi=20.0,
                              u="none", E=2.0885e6, nu=0.3)
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
    foundation strength is cu2 = ratio x cu1. Elastic constants are the paper's
    printed nominal values (Griffiths & Lane 1999, p.390: "in the absence of
    meaningful data for E' and nu', they can be given nominal values (e.g.
    E' = 10^5 kN/m^2 and nu' = 0.3)"), which the paper applies throughout; carried
    to both materials here and converted to English units: E' = 1e5 kPa =
    2.0885e6 psf, nu' = 0.3. SSRM FS is E-invariant; the elastic constants only set
    the displacement scale.

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
    # Elastic constants: the paper's printed nominal values (Griffiths & Lane 1999,
    # p.390 — "in the absence of meaningful data for E' and nu', they can be given
    # nominal values (e.g. E' = 10^5 kN/m^2 and nu' = 0.3)"), which the paper applies
    # throughout. Carried to BOTH materials and converted to English units:
    # E' = 1e5 kPa = 2.0885e6 psf, nu' = 0.3. SSRM FS is E-invariant, so this does not
    # move the locked factor of safety; it sets the displacement scale to the paper's basis.
    E = 2.0885e6
    nu = 0.3
    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    mat = material_cells(1, "Slope clay", 125.0, "mc", c=cu1, phi=0.0,
                         u="none", E=E, nu=nu)
    mat.update(material_cells(2, "Foundation clay", 125.0, "mc", c=cu2, phi=0.0,
                              u="none", E=E, nu=nu))
    u['mat'] = mat
    poly = polygon_cells(1, 1, [(0.0, 100.0), (100.0, 100.0), (200.0, 50.0),
                                (0.0, 50.0), (0.0, 100.0)])
    poly.update(polygon_cells(2, 2, [(0.0, 50.0), (300.0, 50.0), (300.0, 0.0),
                                     (0.0, 0.0), (0.0, 50.0)]))
    u['polygon'] = poly
    # Starting circle (FEM/SSRM ignores it; the LEM companion runs seed='grid',
    # which still carries it in as an extra seed). A deep base circle over the
    # slope: center (180, 165), tangent at y = 20, so it daylights at (50.4, 100)
    # on the crest and (268.3, 50) on the runout — BOTH crossings below the center,
    # which is what makes it a circle the slicer can build. The earlier (150, 90,
    # tangent y = 0) met the crest at (60.6, 100), ten feet ABOVE its own center;
    # that arc is longer than a semicircle, so it produced no failure surface at
    # all and a search seeded with it started from its launch grid instead. See
    # preflight rule surface.circle_daylights_above_center.
    u['circles'] = circle_cells(1, 180.0, 165.0, option="Depth", depth=20.0)
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


def build_griffiths4_r1():
    return build_griffiths4(1.0, "r1")


def build_griffiths4_r2():
    return build_griffiths4(2.0, "r2")


# ---------------------------------------------------------------------------
# Example 3: undrained clay slope with a THIN WEAK LAYER (Fig. 6)
# ---------------------------------------------------------------------------
# The outer profile is Example 4's (H = 50 ft, D = 2): crest platform (0,100)-
# (100,100), 2:1 face to the toe (200,50), 2H runout to (300,50), firm base y=0.
# A thin weaker band (cu2) runs parallel to the slope face, flattens horizontally
# through the foundation, and outcrops at 45 deg beyond the toe.
#
# THE GEOMETRY IS MEASURED OFF FIG. 6 (a schematic — the band THICKNESS is never
# dimensioned in the paper). It was read at 600 dpi by scaling from the labelled
# H-fractions, which close exactly:
#   * Crest (top surface y=100): the band's lower boundary daylights at x=0.6H=30,
#     its upper boundary at x=0.8H=40  (0.6H + 0.2H + 1.2H = 2H, the crest platform).
#   * Sloped reach: both boundaries run parallel to the 2:1 face.
#   * Horizontal foundation reach: the band top sits 0.4H below the runout surface
#     (y=30) and the band bottom 0.4H above the firm base (y=20), a 0.2H=10 ft gap
#     (0.4H + 0.2H + 0.4H = H, runout-surface-to-base).
#   * Outcrop: both boundaries kick UP at 45 deg (1:1) beyond the toe; the upper
#     daylights at x=260 (toe + 1.2H), the lower at x=270 (right edge - 0.6H)
#     (1.2H + 0.2H + 0.6H = 2H, the runout).
# DECLARED GOVERNING THICKNESS: 0.2H = 10 ft in the horizontal foundation reach
# (the clearly-dimensioned 0.4H/0.2H/0.4H stack), tapering to ~0.1H perpendicular
# where the band parallels the face and fanning to meet the ground at the two
# daylights. The band is bounded by two measured polylines; a THICK multiplier
# scales both boundaries about the band centerline for the thickness-sensitivity
# variant (thick=0.5 halves the band everywhere).
#
# Measured boundary polylines (thick = 1.0), y increasing upward, H = 50:
#   upper (nearer the face):  (40,100) (180,30) (240,30) (260,50)
#   lower (deeper):           (30,100) (190,20) (240,20) (270,50)
# centerline = vertex-wise midpoint; each boundary vertex is scaled toward it by
# `thick`.
_G3_U = [(40.0, 100.0), (180.0, 30.0), (240.0, 30.0), (260.0, 50.0)]
_G3_L = [(30.0, 100.0), (190.0, 20.0), (240.0, 20.0), (270.0, 50.0)]
_G3_C = [((u[0] + l[0]) / 2.0, (u[1] + l[1]) / 2.0) for u, l in zip(_G3_U, _G3_L)]
# Local element size for the weak band (polygon sheet Size row). See build_griffiths3.
G3_BAND_SIZE = 2.5


def _g3_boundaries(thick):
    """Return (upper, lower) weak-layer boundary polylines scaled by `thick`
    about the centerline (thick=1.0 is the measured band, 0.5 is half-thickness)."""
    up = [(c[0] + thick * (u[0] - c[0]), c[1] + thick * (u[1] - c[1]))
          for u, c in zip(_G3_U, _G3_C)]
    lo = [(c[0] + thick * (l[0] - c[0]), c[1] + thick * (l[1] - c[1]))
          for l, c in zip(_G3_L, _G3_C)]
    return up, lo


def build_griffiths3(ratio, tag, thick=1.0):
    """Griffiths & Lane (1999) Example 3 (Fig. 6): an undrained (phi_u = 0) clay
    slope on a foundation layer (D = 2) containing a THIN WEAK LAYER. The
    surrounding clay is held at cu1/gamma-H = 0.25; the thin layer strength cu2 is
    varied, and the behaviour is governed by the ratio cu2/cu1.

    Same outer geometry, units and elastic constants as Example 4 (H = 50 ft,
    gamma = 125 pcf, cu1 = 1562.5 psf, phi_u = 0; E'/nu' the paper's printed nominal
    values, Griffiths & Lane 1999 p.390, carried to both materials). The domain is
    tiled by THREE explicit polygons so the thin band is its own material region:
      1. cu1 outer wedge  — between the slope face and the band's upper boundary
      2. cu2 thin layer   — the band itself
      3. cu1 inner body   — everything below the band's lower boundary
    Regions 1 and 3 share material cu1; region 2 is cu2. The weak-layer boundaries
    are traced off Fig. 6 (see module header); `thick` scales the band thickness for
    the sensitivity variant.

    Published (Griffiths & Lane, Fig. 7, FE points, all GRAPHICAL): FOS = 1.47 at
    cu2/cu1 = 1 (Taylor 1937, circular base failure). The curve holds a base-circle
    plateau (~1.4-1.47) for cu2/cu1 > 0.6, then falls roughly linearly once
    cu2/cu1 < 0.6 as a non-circular mechanism follows the weak layer: 0.8->~1.45,
    0.6->~1.40 (transition), 0.5->~1.25, 0.4->~1.05, 0.2->~0.60. The mechanism flips
    from a base circle (Fig. 8a, cu2/cu1=1) to a layer-following two-wedge slide
    (Fig. 8c, cu2/cu1=0.2). At cu2/cu1 = 1 the band is the same material as its
    surroundings, so the model is materially identical to Example 4 r1 (homogeneous
    D = 2) and the anchor should land near its 1.441 (mesh-only difference)."""
    dst = f"{OUTDIR}/xslope_griffiths3_{tag}.xlsx"
    new_file(dst, TEMPLATE)
    cu1 = 1562.5
    cu2 = ratio * cu1
    # Elastic constants: the paper's printed nominal values (Griffiths & Lane 1999,
    # p.390 — E' = 10^5 kN/m^2, nu' = 0.3), converted to English units:
    # E' = 1e5 kPa = 2.0885e6 psf. Carried to BOTH materials (same as Example 4).
    # SSRM FS is E-invariant; this only sets the displacement scale to the paper's basis.
    E = 2.0885e6
    nu = 0.3

    up, lo = _g3_boundaries(thick)

    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    mat = material_cells(1, "Surrounding clay", 125.0, "mc", c=cu1, phi=0.0,
                         u="none", E=E, nu=nu)
    mat.update(material_cells(2, "Thin weak layer", 125.0, "mc", c=cu2, phi=0.0,
                              u="none", E=E, nu=nu))
    u['mat'] = mat

    # 1. cu1 outer wedge: crest surface -> face -> runout -> back up the upper band edge
    outer = [up[0], (100.0, 100.0), (200.0, 50.0), up[3], up[2], up[1], up[0]]
    # 2. cu2 thin weak layer band (upper edge forward, lower edge back)
    band = [up[0], up[1], up[2], up[3], lo[3], lo[2], lo[1], lo[0], up[0]]
    # 3. cu1 inner body: outer domain -> runout -> up the lower band edge -> crest
    inner = [(0.0, 100.0), (0.0, 0.0), (300.0, 0.0), (300.0, 50.0),
             lo[3], lo[2], lo[1], lo[0], (0.0, 100.0)]
    # Local mesh size on the band. The band is 0.2H = 10 ft thick in its horizontal
    # foundation reach and it is the feature that governs the answer at low cu2/cu1:
    # the failure follows it end to end, so the factor of safety is set by how many
    # element rows resolve the shear through it. At the global target sizes these
    # models run (quad8 3.5, tri6 6) the band would otherwise carry 2-3 rows. 2.5 ft
    # puts about four rows across the band on both element types, leaving the rest of
    # the domain at the requested global size.
    #
    # Scaled by `thick` so the half-thickness variant resolves ITS band to the same
    # number of rows. Holding the size fixed would halve the resolution exactly where
    # the thickness is halved, and the thickness-sensitivity comparison would then be
    # reading a mesh difference rather than a thickness difference.
    poly = polygon_cells(1, 1, outer)
    poly.update(polygon_cells(2, 2, band, size=G3_BAND_SIZE * thick))
    poly.update(polygon_cells(3, 1, inner))
    u['polygon'] = poly

    # Starting circle (FEM/SSRM ignores it) — the same deep base circle Example 4
    # carries, on the same outer profile: center (180, 165), tangent at y = 20,
    # daylighting at (50.4, 100) and (268.3, 50), both below the center. See
    # build_griffiths4 for what the earlier (150, 90, tangent y = 0) did.
    u['circles'] = circle_cells(1, 180.0, 165.0, option="Depth", depth=20.0)

    # Limit-equilibrium companion for the weak-ratio station: the paper's own
    # three-line wedge (Griffiths & Lane's Janbu comparison, ~0.47) as a starting
    # surface for the non-circular search — down the band parallel to the face,
    # along the horizontal foundation reach, and up the 45 deg outcrop. Seeded on
    # the band CENTERLINE (_G3_C) so every slice base sits strictly inside the cu2
    # material rather than on a polygon boundary. Only the weak-ratio station gets
    # one: at the other ratios the critical mechanism is the base circle, which the
    # circular search already covers. FEM/SSRM ignores non_circ, so the SSRM locks
    # on this file are unaffected.
    if ratio <= 0.2 and thick == 1.0:
        cl = [(c[0], c[1]) for c in _G3_C]
        u['non-circ'] = noncirc_cells([
            (cl[0][0], cl[0][1], "Free"),   # crest daylight, band centerline
            (cl[1][0], cl[1][1], "Free"),   # foot of the slope reach
            (cl[2][0], cl[2][1], "Free"),   # end of the horizontal reach
            (cl[3][0], cl[3][1], "Free"),   # outcrop daylight beyond the toe
        ])

    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


# Example-3 sweep stations (cu2/cu1, tag), read off Fig. 7: the anchor at 1.0, the
# base-circle plateau (0.8), the transition (0.6), and the falling weak-layer
# regime (0.5, 0.4, 0.2). Gated quad8 runs at the anchor (1.0) and the clearest
# weak-layer mechanism (0.2, Fig. 8c); coarse tri6 at the rest.
GRIFFITHS3_STATIONS = [
    (1.0, "r1"), (0.8, "r0p8"), (0.6, "r0p6"),
    (0.5, "r0p5"), (0.4, "r0p4"), (0.2, "r0p2"),
]


def build_griffiths3_all():
    files = [build_griffiths3(r, tag) for r, tag in GRIFFITHS3_STATIONS]
    # thickness-sensitivity variant: half-thickness band at the representative
    # weak ratio (cu2/cu1 = 0.2), same everything else.
    files.append(build_griffiths3(0.2, "r0p2_thin", thick=0.5))
    return files


def build_griffiths5(l_over_h, tag):
    """Griffiths & Lane (1999) Example 5 (Figs 12-15): the Example-1 homogeneous
    2:1 slope with a HORIZONTAL free surface at depth L below the crest, swept over
    the drawdown ratio L/H ("slow" drawdown). One input file per L/H station; the
    sweep reproduces the Fig. 15 factor-of-safety curve.

    Geometry and material are Example 1's, English units (H = 50 ft, crest at y = 50,
    toe at (160, 0), 2:1 face, firm base at y = 0): ground (0,50)-(60,50)-(160,0),
    gamma = 125 pcf, c' = 312.5 psf (c'/gamma-H = 0.05), phi' = 20 deg, gamma_w = 62.4.
    A CONSTANT TOTAL unit weight is used above and below the water level (paper text):
    gravity loads use total gamma and pore pressures are subtracted at the Gauss points
    (xslope's effective-stress formulation). Elastic constants are the paper's printed
    nominal values (Griffiths & Lane 1999, p.390: E' = 10^5 kN/m^2 = 2.0885e6 psf,
    nu' = 0.3), same as Example 1. SSRM FS is E-invariant; the elastic constants only
    set the displacement scale.

    Water treatment (Figs 12-13), stated once as the free surface exactly as the
    Example 6 full-reservoir build:
      - Free surface: a HORIZONTAL piezometric line at y_fs = 50 - L across the domain.
        u = gamma_w x (vertical depth below the free surface); xslope clamps u = 0 above
        it, so a free surface at the toe (L/H = 1) is the dry slope.
      - Reservoir load: derived by the engine from that same line (main!D23 = auto).
        It is a normal pressure on the submerged outer face, zero at the waterline
        and increasing linearly to gamma_w x y_fs at the toe (Fig. 13), applied as
        consistent boundary tractions. For L/H < 0 (water above the crest, slope
        fully submerged) the crest platform also carries the constant overburden
        pressure gamma_w x (y_fs - 50); for L/H >= 1 there is no submerged face and
        no reservoir load. All three cases follow from the geometry, so the station
        no longer states which one it is.

    Published anchors (Fig. 15): L/H = 0 -> F = 1.85 (Morgenstern, 1963); L/H = 1 ->
    FOS = 1.4 (Bishop & Morgenstern, 1960); the FE curve reaches a minimum FOS ~= 1.3
    at L/H ~= 0.7 (paper text). The curve shape itself is figure-read from Fig. 15."""
    dst = f"{OUTDIR}/xslope_griffiths5_{tag}.xlsx"
    new_file(dst, TEMPLATE)
    H = 50.0
    CREST_Y, TOE_X, TOE_Y = 50.0, 160.0, 0.0
    CREST_R_X = 60.0
    ground = [(0.0, CREST_Y), (CREST_R_X, CREST_Y), (TOE_X, TOE_Y)]
    # Elastic constants: the paper's printed nominal values (Griffiths & Lane 1999,
    # p.390 — E' = 10^5 kN/m^2, nu' = 0.3), converted to English units:
    # E' = 1e5 kPa = 2.0885e6 psf. SSRM FS is E-invariant; this only sets the
    # displacement scale to the paper's basis.
    E = 2.0885e6
    nu = 0.3
    y_fs = CREST_Y - l_over_h * H            # free-surface elevation

    u = {}
    u['main'] = main_cells(gamma_w=62.4)
    u['mat'] = material_cells(1, "soil", 125.0, "mc", c=312.5, phi=20.0,
                              u="piezo", E=E, nu=nu)
    prof = {'B2': TOE_Y}                      # firm base at y = 0 (Example 1, D = 1)
    prof.update(profile_line_cells(1, 1, ground))
    u['profile'] = prof
    # horizontal free surface across the domain (u clamped to 0 above it)
    u['piezo'] = piezo_cells([(0.0, round(y_fs, 3)), (TOE_X, round(y_fs, 3))])
    # placeholder circle (loader requires one; FEM/SSRM ignores it) — a toe circle
    u['circles'] = circle_cells(1, 100.0, 60.0, option="Depth", depth=TOE_Y)

    # The water above the slope is the free surface written above, and nothing
    # else: main!D23 = auto, so the engine measures the piezometric line against
    # the ground and applies the hydrostatic load itself, at whatever L/H this
    # station sits at. Each station used to carry a hand-built block for its own
    # case -- a waterline on the face for 0 <= L/H < 1, a fully submerged crest
    # below that -- which is the branch the derivation now takes from the geometry.
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


# Example-5 drawdown stations (L/H, tag): flat submerged plateau, descent, the
# ~0.7 minimum, and the rise back toward the dry slope — spanning Fig. 15.
GRIFFITHS5_STATIONS = [
    (-0.2, "m0p2"), (0.0, "0"), (0.2, "0p2"), (0.4, "0p4"),
    (0.5, "0p5"), (0.7, "0p7"), (0.9, "0p9"), (1.0, "1"),
]


def build_griffiths5_all():
    return [build_griffiths5(lh, tag) for lh, tag in GRIFFITHS5_STATIONS]


def build_griffiths6_full():
    return _build(f"{OUTDIR}/xslope_griffiths6_full.xlsx", wet=True)


def build_griffiths6_dry():
    return _build(f"{OUTDIR}/xslope_griffiths6_dry.xlsx", wet=False)


if __name__ == "__main__":
    build_griffiths1()
    build_griffiths2()
    build_griffiths3_all()
    build_griffiths4_r1()
    build_griffiths4_r2()
    build_griffiths5_all()
    build_griffiths6_full()
    build_griffiths6_dry()
    build_griffiths6_seep()
    print(f"\nwaterline on upstream face: x = {X_WL:.3f}")
    print("SSRM benchmark input files built.")
