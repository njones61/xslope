"""Builders for RS2-corpus problems that have NO Slide2 counterpart file
(docs/verification/rs2.md). Most RS2 problems run on the Slide2 corpus files;
these are the exceptions - currently the Pruska (2003) multi-program
cross-bearing slopes (RS2 #56-58).

Run from the repo root:  python benchmarks/rocscience/build_rs2.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'files', 'rocscience')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')


def _pruska_slope_data(H, gamma, c, phi, nu):
    """RS2 #56-58 / Pruska (2003) software-comparison slopes: 15 m toe flat,
    slope rising H over a 10 m run, 15 m crest flat, foundation 8 m below the
    toe (fully dimensioned in each section's Figure 1). PUBLISHED elastic
    constants (E = 5000 kPa, nu per case, psi = 0) - not the Griffiths
    convention - since the comparison paper specifies them. Each section's
    lock pair brackets its material family (weakest and strongest case); the
    full 17-case tables live in the rs2.md section."""
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name='soil', c=float(c), phi=float(phi), gamma=float(gamma),
             gamma_sat=float(gamma), option='mc', u='none',
             E=5000.0, nu=float(nu), psi=0.0)
    sd['materials'] = [m]
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 8.0), (15.0, 8.0),
                                                    (25.0, 8.0 + H),
                                                    (40.0, 8.0 + H)]}]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [{'Xo': 20.0, 'Yo': 8.0 + H + 5.0, 'Depth': 6.0,
                      'R': H + 7.0}]
    return sd


def rs2_56a():
    """RS2 #56 case 2 (H=7, gamma=18, c=5, phi=10): SSRM 0.667 vs RS2 0.67,
    Z-Soil 0.71, PLAXIS 0.68, GEO FEM 0.73."""
    sd = _pruska_slope_data(7.0, 18, 5, 10, 0.35)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_56a.xlsx'))
    return 'rs2_56a.xlsx'


def rs2_56b():
    """RS2 #56 case 5 (H=7, gamma=24, c=20, phi=30): SSRM 2.131 vs RS2 2.14,
    Z-Soil 1.98, PLAXIS 2.09, GEO FEM 2.19."""
    sd = _pruska_slope_data(7.0, 24, 20, 30, 0.3)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_56b.xlsx'))
    return 'rs2_56b.xlsx'


def rs2_57a():
    """RS2 #57 case 1 (H=10.5, gamma=18, c=5, phi=10): SSRM 0.449 vs RS2
    0.44, Z-Soil 0.46, PLAXIS 0.44, GEO FEM 0.48."""
    sd = _pruska_slope_data(10.5, 18, 5, 10, 0.35)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_57a.xlsx'))
    return 'rs2_57a.xlsx'


def rs2_57b():
    """RS2 #57 case 6 (H=10.5, gamma=24, c=20, phi=30): SSRM 1.411 vs RS2
    1.42, Z-Soil 1.52, PLAXIS 1.45, GEO FEM 1.54."""
    sd = _pruska_slope_data(10.5, 24, 20, 30, 0.3)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_57b.xlsx'))
    return 'rs2_57b.xlsx'


def rs2_58a():
    """RS2 #58 case 1 (H=14, gamma=18, c=5, phi=10): SSRM 0.342 vs RS2 0.33,
    Z-Soil 0.34, PLAXIS 0.35, GEO FEM 0.35."""
    sd = _pruska_slope_data(14.0, 18, 5, 10, 0.35)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_58a.xlsx'))
    return 'rs2_58a.xlsx'


def rs2_58b():
    """RS2 #58 case 6 (H=14, gamma=24, c=20, phi=30): SSRM 1.057 vs RS2 1.06,
    Z-Soil 1.07, PLAXIS 1.06, GEO FEM 1.10."""
    sd = _pruska_slope_data(14.0, 24, 20, 30, 0.3)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_58b.xlsx'))
    return 'rs2_58b.xlsx'


def hammah_hb1():
    """Hammah, Yacoub, Corkum & Curran (2005), ARMA/USRMS 05-810, Example 1 — the
    verification case for the Hoek-Brown ('hb') strength option, in BOTH the LEM and
    the SSRM. A 10 m high homogeneous 45 deg rock slope in a very weak rock mass.

    Published (their Tables 1-2):
        E = 5000 MPa, nu = 0.3, gamma = 0.025 MN/m3, sigma_ci = 30 MPa,
        GSI = 5, mi = 2, D = 0   ->   mb = 0.067, s = 2.5e-5, a = 0.619
        FE SSR (generalized Hoek-Brown) 1.15 | FE SSR (equivalent Mohr-Coulomb) 1.15
        Bishop simplified 1.153             | Spencer 1.152

    Units are kPa / kN/m3 here (sigma_ci = 30 MPa = 30000 kPa, gamma = 25 kN/m3).

    The paper dimensions only the slope itself (10 m high, 10 m run); the foundation
    depth below the toe and the lateral extents are not given. The answer turns out not
    to care -- the mechanism is toe-exiting, and foundation depths of 2, 4, 6 and 10 m
    all return Bishop 1.150 / Spencer 1.152 -- so 5 m is used and the LEM lock is
    unconditional.
    """
    H, y0, L1, L2 = 10.0, 5.0, 20.0, 20.0
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name='rock', option='hb', hb_sci=30000.0, hb_gsi=5.0, hb_mi=2.0,
             hb_d=0.0, gamma=25.0, gamma_sat=25.0, u='none',
             E=5.0e6, nu=0.3, psi=0.0)
    sd['materials'] = [m]
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [
        (0.0, y0), (L1, y0), (L1 + H, y0 + H), (L1 + H + L2, y0 + H)]}]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    # 'Depth' is the ELEVATION of the circle's lowest point; R must satisfy R = Yo - Depth.
    Xo, Yo, Depth = L1 + H / 2, y0 + H + 6.0, y0 - 1.0
    sd['circles'] = [{'Xo': Xo, 'Yo': Yo, 'Depth': Depth, 'R': Yo - Depth}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'hammah_hb1.xlsx'))
    return 'hammah_hb1.xlsx'


def _rs2_60_slope_data(beta_deg, sci, H=1.0, dfound=2.0, flat=4.0):
    """RS2 #60 -- generalized Hoek-Brown, homogeneous slope, after Li, Merifield &
    Lyamin (2008), "Stability charts for rock slopes based on the Hoek-Brown failure
    criterion", IJRMMS 45, 689-700.

    THE MANUAL DOES NOT STATE sigma_ci. RS2's Table 1 gives only H = 1 m,
    gamma = 23 kN/m3, nu = 0.3, GSI = 70, mi = 15 (and D = 0 comes from Li's text).
    sigma_ci is back-computed from Li's Table 1, which tabulates the CRITICAL ratio
    (sigma_ci / gamma H)_crit -- the value at which collapse has just occurred, i.e.
    F = 1. So sigma_ci is *chosen* to put the slope at limiting equilibrium, and

        THE VERIFICATION TARGET IS FS = 1.0 BY CONSTRUCTION.

    It is not an independently computed reference factor of safety.

        beta = 15 deg -> (sci/gH)crit = 0.026 -> sci = 0.598 kPa
        beta = 30 deg -> (sci/gH)crit = 0.075 -> sci = 1.725 kPa
        beta = 45 deg -> (sci/gH)crit = 0.176 -> sci = 4.048 kPa

    UNITS. These look absurd for rock, and that is the trap. H = 1 m makes
    gamma*H = 23 kPa, and the critical ratio is LESS THAN ONE, so sigma_ci is a
    fraction of gamma*H -- sub-kPa to a few kPa. The problem is normalized (only the
    ratio sigma_ci/(gamma H) matters); a 1 m slope in 0.6 kPa rock is the same problem
    as a 100 m slope in 60 kPa rock. Entering sigma_ci in MPa, as Hoek-Brown convention
    would invite, overstates the strength by 1000x and the slope becomes trivially
    stable.

    Li's Table 1 prints its last block as "beta = 10", but the body text, the charts
    (Fig. 5 is beta = 15, there is no beta = 10 chart) and the base-failure discussion
    all say 15 -- it is a typo in the paper. RS2 evidently read it as 15 too: its Slide2
    value for case 1 (1.011) reproduces Li's own F1 for that row (1.010).

    Published: RS2 SSRM 1.02 / 1.02 / 1.10 (its figure reads 1.11 for the 45 deg case);
    Slide2 Spencer 1.011 / 0.992 / 1.035; Li's reference F = 1 in all three.

    Domain: Li reports the depth factor d/H is insignificant and RS2 says the figure's
    overall extents "were shown to be insignificant", so neither is dimensioned. The
    foundation still has to be deep enough to let the mechanism form: beta = 15 deg is
    the one case Li reports as a BASE failure rather than a toe failure.
    """
    import math
    run = H / math.tan(math.radians(beta_deg))
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name='rock', option='hb', hb_sci=float(sci), hb_gsi=70.0, hb_mi=15.0,
             hb_d=0.0, gamma=23.0, gamma_sat=23.0, u='none',
             E=1.0e5, nu=0.3, psi=0.0)
    sd['materials'] = [m]
    y0 = dfound
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [
        (0.0, y0), (flat, y0), (flat + run, y0 + H), (flat + run + flat, y0 + H)]}]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    Xo, Yo = flat + run / 2.0, y0 + H + 0.8 * H
    Depth = y0 - 0.5 * H
    sd['circles'] = [{'Xo': Xo, 'Yo': Yo, 'Depth': Depth, 'R': Yo - Depth}]
    return sd


def rs2_60a():
    """RS2 #60 case 1: beta = 15 deg, sigma_ci = 0.598 kPa. Base failure (Li)."""
    sd = _rs2_60_slope_data(15.0, 0.026 * 23.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_60a.xlsx'))
    return 'rs2_60a.xlsx'


def rs2_60b():
    """RS2 #60 case 2: beta = 30 deg, sigma_ci = 1.725 kPa."""
    sd = _rs2_60_slope_data(30.0, 0.075 * 23.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_60b.xlsx'))
    return 'rs2_60b.xlsx'


def rs2_60c():
    """RS2 #60 case 3: beta = 45 deg, sigma_ci = 4.048 kPa."""
    sd = _rs2_60_slope_data(45.0, 0.176 * 23.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_60c.xlsx'))
    return 'rs2_60c.xlsx'


def _poly_slope_data(polygons, materials, circle, max_depth,
                     piezo_line=None, piezo_phreatic=False):
    """Assemble a polygon-input slope_data for RS2 Part-III problems whose
    geometry is a set of material zones (not stacked profile lines).

    ``polygons`` is a list of (mat_id, [(x, y), ...]) exterior rings; ``materials``
    a list of dicts merged onto the ACADS-1a template material (so every physics
    field the loader expects is present); ``circle`` a single starting circle for
    the LEM search; ``max_depth`` the search floor. Geometry coordinates are
    transcribed from the vendor RS2 .fez models (dev-side cross-check only) and
    hard-coded here so the corpus rebuilds without them.

    ``piezo_line`` optionally supplies a single piezometric surface [(x, y), ...];
    when given, each material spec should carry ``u='piezo'`` and pore pressure is
    computed from the surface. ``piezo_phreatic`` sets the phreatic-inclination
    (cos^2) correction flag, matching the RS2 groundwater 'phreatic' surface type."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    mats = []
    for spec in materials:
        m = dict(base)
        m.update(option='mc', u='none', psi=0.0)
        m.update(spec)
        mats.append(m)
    sd['materials'] = mats
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': mid, 'polygon': Polygon(coords)}
                      for mid, coords in polygons]
    sd['max_depth'] = float(max_depth)
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['piezo_line'] = list(piezo_line) if piezo_line else []
    sd['piezo_phreatic'] = bool(piezo_phreatic)
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [circle]
    return sd


def rs2_59():
    """RS2 #59 (Part III) -- Stability of a Three-Layered Soil Slope, the Budapest
    (Rozsadomb) landslide after

        Gorog, P. & Torok, A. (2007), "Slope stability assessment of weathered clay
        by using field data and computer modelling: a case study from Budapest,"
        Natural Hazards 45.

    A ~415 m wide x ~75 m tall real-coordinate hillslope in three layers: a
    YellowClay/Debris cover (c = 50, phi = 15, gamma = 19), a thin weak WASTE lens
    (c = 1, phi = 5, gamma = 14) and a strong GreyClay base (c = 250, phi = 30,
    gamma = 22). The waste lens is a wedge that exists only in the left ~136 m,
    ~12.6 m thick at the toe face and tapering to zero at x = 136; it daylights on
    the toe and the critical mechanism is a NON-CIRCULAR translational slip riding
    the top of the lens. A circular search misfinds the deeper competing surface
    (FS ~ 1.9), so this is an SSRM problem -- the SSRM localizes the shear band
    through the c = 1 / phi = 5 lens on its own.

    Published (manual results table, Case 1 constant E = 50000 kPa, nu = 0.4):
    Slide2 1.567 | RS2 SSR 1.57 | PLAXIS 1.6. Case 2 (varying E: GreyClay 20000,
    YellowClay/Debris 18000, Waste 2000) reports RS2 SSR 1.56 / Slide2 1.567 /
    PLAXIS 1.6 -- an E-only change that SSRM FS is insensitive to, so it is not a
    separate XSLOPE case. Elastic constants are the published Case-1 values (the
    .fez reader does not parse E/nu); psi = 0 (Griffiths convention). Dry -- no
    water table, tension crack, seismic or loads in the model.

    Geometry transcribed from the RS2 vendor model 'slope stability #059_01.fez'
    via the .fez zone polygonizer; the circle is an inert placeholder the loader
    requires (SSRM only -- no LEM surface is defined on an RS2 SSR model)."""
    zones = [
        (0, [(0.858, 120.076), (0.858, 122.016), (115.713, 126.355),
             (132.493, 136.770), (147.248, 136.192), (267.021, 157.022),
             (359.600, 164.833), (414.986, 175.537), (414.986, 158.868),
             (350.052, 148.921), (223.914, 135.613), (135.965, 116.229)]),
        (1, [(0.858, 107.469), (0.858, 120.076), (135.965, 116.229),
             (73.474, 106.682)]),
        (2, [(414.986, 158.868), (414.986, 100.317), (0.858, 100.317),
             (0.858, 107.469), (73.474, 106.682), (135.965, 116.229),
             (223.914, 135.613), (350.052, 148.921)]),
    ]
    sd = _poly_slope_data(
        polygons=zones,
        materials=[
            dict(name='YellowClay/Debris', c=50.0, phi=15.0, gamma=19.0,
                 gamma_sat=19.0, E=50000.0, nu=0.4),
            dict(name='Waste', c=1.0, phi=5.0, gamma=14.0, gamma_sat=14.0,
                 E=50000.0, nu=0.4),
            dict(name='GreyClay', c=250.0, phi=30.0, gamma=22.0, gamma_sat=22.0,
                 E=50000.0, nu=0.4),
        ],
        circle={'Xo': 100.0, 'Yo': 160.0, 'Depth': 100.0, 'R': 60.0},
        max_depth=0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_59.xlsx'))
    return 'rs2_59.xlsx'


def rs2_63():
    """RS2 #63 -- Slope Stability Assessment of a Homogeneous Slope, after Cheng,
    Lansivaara & Wei (2007). An 11 m homogeneous slope (c = 10 kPa, phi = 30 deg,
    gamma = 20 kN/m3). Published: Slide2 1.380, RS2 SSRM 1.38, Cheng LEM 1.3830.
    Geometry from RS2 vendor model 'slope stability #063.fez'."""
    poly = [(46.0, 0.0), (0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (16.0, 16.0),
            (26.0, 21.0), (46.0, 21.0), (46.0, 0.0)]
    sd = _poly_slope_data(
        polygons=[(0, poly)],
        materials=[dict(name='soil', c=10.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                        E=14000.0, nu=0.3)],
        circle={'Xo': 8.0, 'Yo': 24.0, 'Depth': 10.0, 'R': 14.0},
        max_depth=0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_63.xlsx'))
    return 'rs2_63.xlsx'


def rs2_61a():
    """RS2 #61 case 1 -- Local and Global Minima for a Homogeneous Slope, after
    Cheng, Lansivaara & Wei (2007). Case 1 is the UNCONSTRAINED global minimum
    (the other three cases fence a Polygon Search Area onto successive local
    minima, which xslope has no direct analog for -- see the rs2.md row).
    Homogeneous benched slope, c = 5 kPa, phi = 30 deg, gamma = 20 kN/m3.
    Published for case 1: Slide2 1.336, Cheng LEM 1.327, RS2 SSRM 1.35.

    A grid seed traps on a steeper local circle (FS 1.44) -- exactly the trap the
    paper is about -- so the search is seeded with a toe-to-crest circle that
    refines onto the global minimum. Geometry from 'slope stability #061_01.fez'."""
    poly = [(54.324, 0.073), (-0.057, 0.073), (-0.057, 9.965), (10.118, 9.965),
            (20.011, 14.883), (26.003, 19.858), (31.939, 19.858), (34.878, 22.797),
            (42.905, 26.811), (54.324, 26.811), (54.324, 0.073)]
    sd = _poly_slope_data(
        polygons=[(0, poly)],
        materials=[dict(name='soil', c=5.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                        E=14000.0, nu=0.3)],
        circle={'Xo': 25.0, 'Yo': 40.0, 'Depth': 0.5, 'R': 39.5},
        max_depth=0.073)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_61a.xlsx'))
    return 'rs2_61a.xlsx'


def _rs2_66_slope_data(h1):
    """RS2 #66 (Part III) -- Embankment Basal Stability, after

        Nakamura, A., Cai, F. & Ugai, K. (2008), "Embankment basal stability
        analysis using shear strength reduction finite element method."

    A 10 m high embankment (1.5:1 side slopes, 20 m crest) is built on a two-layer
    foundation: a SOFT upper stratum (phi = 0, c = 35 kPa) of thickness h1, over a
    firm 10 m bearing stratum (phi = 0, c = 100 kPa). The soft-layer thickness h1 is
    the varied parameter (2, 4, 6, 8, 10 m); the embankment and the bearing stratum
    are held constant. Fill is a c = 0, phi = 35 deg granular material. gamma = 18.82
    kN/m3 throughout. The mechanism is a basal squeeze through the soft clay, so the
    factor of safety is governed by the phi = 0 band, not by the fill.

    The paper (and RS2) set the dilation angle psi = phi for all soils; XSLOPE's
    SSRM runs non-associated (psi = 0, the Griffiths convention the corpus uses).
    For a basal mechanism inside the phi = 0 clay the flow rule is dilationless there
    regardless, so the difference is confined to the granular fill and is small.

    Geometry transcribed from the RS2 vendor models 'slope stability #066_0N (h=..).fez':
    external domain 150 m wide, bearing stratum 0..10 m, soft stratum 10..(10+h1),
    embankment base at y = 10 + h1 from x = 50..100, crest y = 20 + h1 from x = 65..85.
    """
    top = 10.0 + h1
    crest = 20.0 + h1
    bearing = [(0.0, 0.0), (150.0, 0.0), (150.0, 10.0), (0.0, 10.0)]
    soft = [(0.0, 10.0), (150.0, 10.0), (150.0, top), (0.0, top)]
    embank = [(50.0, top), (100.0, top), (85.0, crest), (65.0, crest)]
    sd = _poly_slope_data(
        polygons=[(0, embank), (1, soft), (2, bearing)],
        materials=[
            dict(name='embankment', c=0.0, phi=35.0, gamma=18.82, gamma_sat=18.82,
                 E=1.0e5, nu=0.3),
            dict(name='soft ground', c=35.0, phi=0.0, gamma=18.82, gamma_sat=18.82,
                 E=1.0e5, nu=0.3),
            dict(name='bearing stratum', c=100.0, phi=0.0, gamma=18.82,
                 gamma_sat=18.82, E=1.0e5, nu=0.3),
        ],
        circle={'Xo': 75.0, 'Yo': crest + 20.0, 'Depth': 9.0,
                'R': crest + 20.0 - 9.0},
        max_depth=0.0)
    return sd


def rs2_66a():
    """RS2 #66 case 1, h1 = 2 m. Published: Slide2 Spencer 1.05, RS2 SSR 1.13,
    Nakamura LEM 1.21, Nakamura FEM 1.24."""
    sd = _rs2_66_slope_data(2.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_66a.xlsx'))
    return 'rs2_66a.xlsx'


def rs2_66b():
    """RS2 #66 case 2, h1 = 4 m. Published: Slide2 Spencer 1.16, RS2 SSR 1.19,
    Nakamura LEM 1.22, Nakamura FEM 1.16."""
    sd = _rs2_66_slope_data(4.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_66b.xlsx'))
    return 'rs2_66b.xlsx'


def rs2_66c():
    """RS2 #66 case 3, h1 = 6 m. Published: Slide2 Spencer 1.10, RS2 SSR 1.13,
    Nakamura LEM 1.22, Nakamura FEM 1.16."""
    sd = _rs2_66_slope_data(6.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_66c.xlsx'))
    return 'rs2_66c.xlsx'


def rs2_66d():
    """RS2 #66 case 4, h1 = 8 m. Published: Slide2 Spencer 1.13, RS2 SSR 1.08,
    Nakamura LEM 1.10, Nakamura FEM 1.10."""
    sd = _rs2_66_slope_data(8.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_66d.xlsx'))
    return 'rs2_66d.xlsx'


def rs2_66e():
    """RS2 #66 case 5, h1 = 10 m. Published: Slide2 Spencer 1.05, RS2 SSR 1.05,
    Nakamura LEM 1.08, Nakamura FEM 1.08."""
    sd = _rs2_66_slope_data(10.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_66e.xlsx'))
    return 'rs2_66e.xlsx'


_RS2_62_SOILS = [
    dict(name='Soil 1', c=20.0, phi=35.0, gamma=19.0, gamma_sat=19.0, E=14000.0, nu=0.3),
    dict(name='Soil 2 (soft band)', c=0.0, phi=25.0, gamma=19.0, gamma_sat=19.0,
         E=14000.0, nu=0.3),
    dict(name='Soil 3', c=10.0, phi=35.0, gamma=19.0, gamma_sat=19.0, E=14000.0, nu=0.3),
]


def _rs2_62_slope_data(polys):
    """RS2 #62 (Part III) -- Stability of a Three-Layered Slope With a Soft Band, after

        Cheng, Y.M., Lansivaara, T. & Wei, W.B. (2007), "Two-dimensional slope stability
        analysis by limit equilibrium and strength reduction methods." Comput. Geotech. 34.

    A 10 m slope carries a thin SOFT band (Soil 2: c = 0, phi = 25 deg) dipping through it,
    between a stronger cap (Soil 1: c = 20, phi = 35) and base (Soil 3: c = 10, phi = 35);
    gamma = 19, E = 14 MPa, nu = 0.3 throughout. Three geometries vary the band's DAYLIGHT
    width (Analysis I / II / III = 28 / 20 / 12 m domains), and each was run at two dilation
    angles (Case 1 psi = 0, Case 2 psi = phi). XSLOPE's SSRM is non-associated only (psi = 0),
    so ONLY the Case-1 (psi = 0) column is reproducible here; the Case-2 associated-flow column
    and Flac3D's much higher associated-flow values are recorded in rs2.md for context.

    ``polys`` are the three zone polygons (mat_id 0 = Soil 1, 1 = soft band, 2 = Soil 3),
    transcribed from the RS2 vendor models 'slope stability #062_0N.fez' via the .fez zone
    polygonizer and hard-coded so the corpus rebuilds without them. SSRM only; the circle is
    an inert placeholder the loader requires."""
    sd = _poly_slope_data(
        polygons=polys, materials=[dict(m) for m in _RS2_62_SOILS],
        circle={'Xo': 10.0, 'Yo': 20.0, 'Depth': 3.0, 'R': 17.0}, max_depth=0.0)
    return sd


def rs2_62a():
    """RS2 #62 Analysis I (28 m domain), psi = 0. Published (psi = 0): RS2 SSR 0.88,
    Plaxis 0.86 | Flac3D 1.64 (associated). Case 2 psi = phi: RS2 0.98 (not reproducible)."""
    polys = [
        (0, [(5.0, 5.0), (8.0, 8.0), (20.0, 15.0), (28.0, 15.0), (28.0, 10.0), (8.0, 7.5)]),
        (1, [(28.0, 10.0), (28.0, 9.5), (8.0, 7.1), (5.0, 4.5), (0.0, 5.0), (5.0, 5.0),
             (8.0, 7.5)]),
        (2, [(28.0, 9.5), (28.0, 0.0), (0.0, 0.0), (0.0, 5.0), (5.0, 4.5), (8.0, 7.1)]),
    ]
    save_slope_data_to_xlsx(_rs2_62_slope_data(polys), os.path.join(OUT, 'rs2_62a.xlsx'))
    return 'rs2_62a.xlsx'


def rs2_62b():
    """RS2 #62 Analysis II (20 m domain), psi = 0. Published (psi = 0): RS2 SSR 0.89,
    Plaxis 0.85 | Flac3D 1.30 (associated). Case 2 psi = phi: RS2 0.98 (not reproducible)."""
    polys = [
        (0, [(5.0, 5.0), (8.0, 8.0), (20.0, 15.0), (20.0, 8.911), (8.0, 7.5)]),
        (1, [(20.0, 8.911), (20.0, 8.287), (8.0, 7.1), (5.0, 4.5), (0.0, 5.0), (5.0, 5.0),
             (8.0, 7.5)]),
        (2, [(20.0, 8.287), (20.0, 0.0), (0.0, 0.0), (0.0, 5.0), (5.0, 4.5), (8.0, 7.1)]),
    ]
    save_slope_data_to_xlsx(_rs2_62_slope_data(polys), os.path.join(OUT, 'rs2_62b.xlsx'))
    return 'rs2_62b.xlsx'


def rs2_62c():
    """RS2 #62 Analysis III (12 m domain), psi = 0. Published (psi = 0): RS2 SSR 0.81,
    Plaxis 0.82 | Flac3D 1.03 (associated). Case 2 psi = phi: RS2 0.93 (not reproducible)."""
    polys = [
        (0, [(5.0, 5.0), (8.0, 8.0), (12.0, 10.333), (12.0, 7.97), (8.0, 7.5)]),
        (1, [(12.0, 7.97), (12.0, 7.496), (8.0, 7.1), (5.0, 4.5), (0.0, 5.0), (5.0, 5.0),
             (8.0, 7.5)]),
        (2, [(12.0, 7.496), (12.0, 0.0), (0.0, 0.0), (0.0, 5.0), (5.0, 4.5), (8.0, 7.1)]),
    ]
    save_slope_data_to_xlsx(_rs2_62_slope_data(polys), os.path.join(OUT, 'rs2_62c.xlsx'))
    return 'rs2_62c.xlsx'


# RS2 #65 -- eight materials of the Padina tailings dam (Table 1 / the .fez).
# c, phi, gamma are Table 1 = the .fez densities (g/cm3 x 10); E, nu are Table 1
# (absent from the .fez, which imports E = 0). psi = 0 (Griffiths convention).
# For a static-groundwater effective-stress run the pore pressure comes from the
# piezometric surface, so gamma_sat = gamma (one unit weight + explicit water).
_RS2_65_SOILS = [
    dict(name='Rockfill Lyulyaka', c=20.0, phi=38.0, gamma=18.6, gamma_sat=18.6,
         E=75000.0, nu=0.30, u='piezo'),
    dict(name='Fill', c=22.5, phi=33.7, gamma=18.9, gamma_sat=18.9,
         E=70000.0, nu=0.31, u='piezo'),
    dict(name='Rockfill G.Sakar', c=20.0, phi=38.0, gamma=18.6, gamma_sat=18.6,
         E=75000.0, nu=0.30, u='piezo'),
    dict(name='Counterfill', c=22.5, phi=33.7, gamma=18.9, gamma_sat=18.9,
         E=70000.0, nu=0.31, u='piezo'),
    dict(name='Tailings', c=0.0, phi=34.8, gamma=13.3, gamma_sat=13.3,
         E=16100.0, nu=0.35, u='piezo'),
    dict(name='Alluvial Clay', c=0.0, phi=24.65, gamma=19.8, gamma_sat=19.8,
         E=16300.0, nu=0.34, u='piezo'),
    dict(name='Marly Clay', c=0.0, phi=19.5, gamma=22.2, gamma_sat=22.2,
         E=38000.0, nu=0.33, u='piezo'),
    dict(name='Marl', c=30.0, phi=24.5, gamma=24.0, gamma_sat=24.0,
         E=75000.0, nu=0.30, u='piezo'),
]

# The 14-point phreatic surface (RS2 groundwater, 'phreatic' type -> cos^2 flag).
_RS2_65_PIEZO = [
    (-70.0, 9.0), (-17.5, 9.0), (-1.769, 12.739), (15.821, 13.65),
    (80.77, 13.65), (79.941, 16.754), (82.344, 21.909), (87.613, 28.198),
    (94.921, 34.997), (103.165, 40.181), (106.171, 42.117), (107.511, 42.117),
    (116.167, 46.0), (155.5, 46.0),
]


def rs2_65():
    """RS2 #65 (Part III) -- Slope Stability Assessment of a Tailings Dam, the
    Padina tailings dam after Tzenkov (2008), RS2 Slope Stability Verification
    Manual Part III Problem 65 (Table 1, pp. 230-231).

    A 225 m wide x 77 m tall cross-section of an eight-material tailings dam with a
    phreatic surface. Twelve zones tile the domain (union area = domain area =
    13262 m2, no gaps/overlaps): a Marl base (mat 7), Marly-Clay and Alluvial-Clay
    bands (mat 6, 5), a Counterfill body (mat 3), the Tailings core (mat 4) and the
    Rockfill/Fill embankment shells (mats 0-2). Pore pressure is applied from a
    single 14-point phreatic surface connected to every material (static
    groundwater; the .fez carries no FE seepage solution to read).

    Published (manual results table): RS2 SSRM 1.29 (the FE target) | Slide2
    circular 1.41 | Slide2 non-circular 1.33 | reference LEM 1.39 | reference FEM
    1.41. Geometry + phreatic line transcribed from 'slope stability #065.fez' via
    the .fez zone polygonizer; the circle is an inert placeholder the loader
    requires (SSRM only -- no LEM surface is defined on an RS2 SSR model)."""
    zones = [
        (7, [(155.5, 6.808), (155.5, -30.0), (-70.0, -30.0), (-70.0, -11.429),
             (-37.925, -8.257), (31.51, -1.914), (86.672, 3.183), (111.139, 4.769),
             (155.5, 6.808)]),
        (6, [(155.5, 8.766), (155.5, 6.808), (111.139, 4.769), (86.672, 3.183),
             (31.51, -1.914), (-37.925, -8.257), (-70.0, -11.429), (-70.0, -2.367),
             (-55.142, -1.008), (-3.377, 3.07), (49.746, 5.335), (115.103, 7.488),
             (155.5, 8.766)]),
        (5, [(155.5, 10.933), (155.5, 8.766), (115.103, 7.488), (49.746, 5.335),
             (-3.377, 3.07), (-55.142, -1.008), (-70.0, -2.367), (-70.0, 9.0),
             (-17.5, 9.0), (10.263, 9.896), (36.87, 10.265), (84.954, 10.933),
             (155.5, 10.933)]),
        (4, [(101.0, 46.0), (116.167, 46.0), (155.5, 46.0), (155.5, 10.933),
             (84.954, 10.933), (66.734, 22.767), (66.734, 23.868), (74.159, 23.868),
             (64.72, 30.287), (79.823, 30.287), (69.754, 36.957), (92.786, 36.957),
             (84.102, 42.117), (103.165, 42.117), (107.511, 42.117), (101.0, 46.0)]),
        (1, [(96.01, 47.0), (100.0, 47.0), (101.0, 46.0), (107.511, 42.117),
             (103.165, 42.117), (96.01, 47.0)]),
        (2, [(82.5, 43.0), (89.0, 47.0), (96.01, 47.0), (103.165, 42.117),
             (84.102, 42.117), (82.5, 43.0)]),
        (2, [(68.0, 37.5), (76.0, 43.0), (82.5, 43.0), (84.102, 42.117),
             (92.786, 36.957), (69.754, 36.957), (68.0, 37.5)]),
        (0, [(44.0, 30.0), (56.0, 37.5), (68.0, 37.5), (69.754, 36.957),
             (79.823, 30.287), (64.72, 30.287), (44.0, 30.0)]),
        (0, [(23.0, 18.5), (40.5, 30.0), (44.0, 30.0), (55.037, 22.767),
             (36.87, 10.265), (10.263, 9.896), (23.0, 18.5)]),
        (3, [(-17.5, 9.0), (-2.0, 18.5), (23.0, 18.5), (10.263, 9.896),
             (-17.5, 9.0)]),
        (0, [(84.954, 10.933), (36.87, 10.265), (55.037, 22.767), (66.734, 22.767),
             (84.954, 10.933)]),
        (0, [(64.72, 30.287), (74.159, 23.868), (66.734, 23.868), (66.734, 22.767),
             (55.037, 22.767), (44.0, 30.0), (64.72, 30.287)]),
    ]
    sd = _poly_slope_data(
        polygons=zones,
        materials=[dict(m) for m in _RS2_65_SOILS],
        circle={'Xo': 75.0, 'Yo': 70.0, 'Depth': 0.0, 'R': 70.0},
        max_depth=0.0,
        piezo_line=_RS2_65_PIEZO, piezo_phreatic=True)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_65.xlsx'))
    return 'rs2_65.xlsx'


def rs2_51():
    """RS2 Part-4 Verification Problem #51 -- a four-material slope with a water
    table, a tension crack in the top layer and a horizontal seismic coefficient
    k = 0.1, evaluated on a GIVEN circular surface with twelve LEM methods, after

        Zhu, D.Y., Lee, C.F. & Jiang, H.D. (2003), "Generalised framework of
        limit equilibrium methods for slope stability analysis," Geotechnique
        53(4), 377-395 (the RS2/Slide2 Slope Stability Verification Manual's
        Problem #51).

    Four Mohr-Coulomb layers dip parallel to the 1V:2H slope face (toe (0,0) ->
    crest (60,30)): Layer 1 (top, c = 20, phi = 32, gamma = 18.2), Layer 2
    (c = 25, phi = 30, gamma = 18.0), Layer 3 (the weak band, c = 40, phi = 18,
    gamma = 18.5) and Layer 4 (bottom, c = 40, phi = 28, gamma = 18.8). Pore
    pressure is taken from a 9-point piezometric surface connected to every
    material; the horizontal seismic coefficient is 0.1 and a dry tension crack
    of the Rankine active depth h_c = 2c/(gamma*sqrt(Ka)) ~ 3.97 m sits in the top
    layer. gamma_water = 9.81; psi = 0 / E, nu placeholders (LEM only -- the
    elastic constants are not used).

    Published (manual Table 51.2) Slide2 | Zhu: Ordinary 1.145 | 1.066; Bishop
    1.278 | 1.278; Janbu Simplified 1.112 | 1.112; Corps of Engineers 1.422 |
    1.377; Lowe & Karafiath 1.288 | 1.290; Spencer 1.293 | 1.293; GLE/M-P 1.304 |
    1.303. RS2 SSR (a different, FE mechanism) = 1.22.

    PARTIAL reconstruction: the vendor '.fez' is an RS2 SSR model that carries NO
    LEM slip surface, and the given circle + tension-crack depth are figure-only
    (Figs 51.1-51.3, not in the '.fez' or the manual text). The circle here is a
    slope-face circle recovered by inversion against the rigorous methods: centre
    (32, 36), tangent at y = 1.0 (Depth option -> R = Yo - Depth = 35), which
    daylights from the lower slope face (x ~ 13) to the back plateau (x ~ 66). On
    it Spencer reproduces to +0.5% and Janbu-simplified (xslope reports Janbu-
    corrected, so divide out fo ~ 1.08) to within 0.5%. Geometry, materials, piezo
    line and k = 0.1 are transcribed from 'slope stability #051.fez' (k = 0.1 is
    the '.fea' body force bx = -0.1, which the reader does not auto-import)."""
    import math
    zones = [
        (2, [(-40.0, -15.0), (-40.0, 0.0), (0.0, 0.0), (10.0, 5.0), (100.0, 5.0),
             (100.0, -5.0), (100.0, -6.25)]),                        # Layer 3 (weak)
        (3, [(100.0, -6.25), (100.0, -30.0), (-40.0, -30.0), (-40.0, -15.0)]),  # Layer 4
        (1, [(100.0, 20.0), (100.0, 5.0), (10.0, 5.0), (40.0, 20.0)]),          # Layer 2
        (0, [(100.0, 25.0), (100.0, 20.0), (40.0, 20.0), (50.0, 25.0)]),        # Layer 1
        (0, [(50.0, 25.0), (60.0, 30.0), (100.0, 30.0), (100.0, 25.0)]),        # Layer 1 (top split)
    ]
    materials = [
        dict(name='Layer 1 (top)', c=20.0, phi=32.0, gamma=18.2, gamma_sat=18.2,
             E=1.0e5, nu=0.3, u='piezo'),
        dict(name='Layer 2', c=25.0, phi=30.0, gamma=18.0, gamma_sat=18.0,
             E=1.0e5, nu=0.3, u='piezo'),
        dict(name='Layer 3 (weak)', c=40.0, phi=18.0, gamma=18.5, gamma_sat=18.5,
             E=1.0e5, nu=0.3, u='piezo'),
        dict(name='Layer 4 (bottom)', c=40.0, phi=28.0, gamma=18.8, gamma_sat=18.8,
             E=1.0e5, nu=0.3, u='piezo'),
    ]
    piezo = [(-40.0, 0.0), (0.0, 0.0), (8.0, 3.0), (20.0, 7.0), (30.0, 10.0),
             (43.0, 12.0), (55.0, 14.0), (70.0, 15.0), (100.0, 15.0)]
    sd = _poly_slope_data(
        polygons=zones, materials=materials,
        circle={'Xo': 32.0, 'Yo': 36.0, 'R': 35.0, 'Depth': 1.0},
        max_depth=0.0, piezo_line=piezo)
    sd['k_seismic'] = 0.1                      # .fea body force bx = -0.1 (not auto-imported)
    # dry Rankine tension crack in the top layer (c = 20, phi = 32, gamma = 18.2)
    ka = (1.0 - math.sin(math.radians(32.0))) / (1.0 + math.sin(math.radians(32.0)))
    sd['tcrack_depth'] = round(2.0 * 20.0 / (18.2 * math.sqrt(ka)), 3)   # ~3.965 m
    sd['tcrack_water'] = 0.0
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_51.xlsx'))
    return 'rs2_51.xlsx'


if __name__ == '__main__':
    for fn in (rs2_56a, rs2_56b, rs2_57a, rs2_57b, rs2_58a, rs2_58b, hammah_hb1,
               rs2_60a, rs2_60b, rs2_60c, rs2_61a, rs2_59, rs2_63,
               rs2_66a, rs2_66b, rs2_66c, rs2_66d, rs2_66e,
               rs2_62a, rs2_62b, rs2_62c, rs2_65, rs2_51):
        print(fn())
