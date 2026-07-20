"""Builders for RS2-corpus problems that have NO Slide2 counterpart file
(docs/verification/rs2.md). Most RS2 problems run on the Slide2 corpus files;
these are the exceptions - currently the Pruska (2003) multi-program
cross-bearing slopes (RS2 #56-58).

Run from the repo root:  python benchmarks/rocscience/build_rs2.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.dirname(__file__))

from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx as _write_xlsx  # noqa: E402
from elastic_props import assign_elastic_props  # noqa: E402


def save_slope_data_to_xlsx(slope_data, path):
    """Write an input file, assigning physical FEM elastic constants (E, nu) by soil
    type to any material lacking deliberately-set values, then delegating to fileio.

    Most RS2 builders set E/nu explicitly (published vendor .fez values, the Pruska
    Poisson-ratio study, HB rock), and those non-default values are preserved by
    assign_elastic_props(); only a material still carrying the inherited unit-blind
    default (E = 100,000, nu = 0.3) or no E at all is (re)assigned a soil-type value
    in the file's own unit system. See build_problems.save_slope_data_to_xlsx and
    elastic_props.py for the full rationale.
    """
    assign_elastic_props(slope_data.get('materials', []))
    return _write_xlsx(slope_data, path)


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
    sigma_ci is taken directly from the RS2 vendor .fez model (the material's
    Generalized Hoek-Brown `ucs` field), which is the value RS2/Slide2 actually
    solved. Li's Table 1 tabulates a CRITICAL ratio (sigma_ci / gamma H)_crit -- the
    value at which collapse has just occurred, i.e. F = 1 -- and only case a's ratio
    reproduces the vendor sigma_ci; cases b and c do not (Li's 0.075 / 0.176 give
    1.725 / 4.048 kPa, but the vendor files carry 1.61 / 4.37 kPa, i.e. ratios of
    0.070 / 0.190). Using the vendor values makes the reconstruction faithful to the
    RS2 model and reproduces Slide2's own Spencer results for cases b and c.

        beta = 15 deg -> vendor ucs = 0.598 kPa  (== Li ratio 0.026 * 23)
        beta = 30 deg -> vendor ucs = 1.61  kPa  (Li ratio 0.075 -> 1.725 does NOT match)
        beta = 45 deg -> vendor ucs = 4.37  kPa  (Li ratio 0.176 -> 4.048 does NOT match)

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
    """RS2 #60 case 1: beta = 15 deg, vendor sigma_ci = 0.598 kPa. Base failure (Li)."""
    sd = _rs2_60_slope_data(15.0, 0.598)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_60a.xlsx'))
    return 'rs2_60a.xlsx'


def rs2_60b():
    """RS2 #60 case 2: beta = 30 deg, vendor sigma_ci = 1.61 kPa (.fez #060_02)."""
    sd = _rs2_60_slope_data(30.0, 1.61)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_60b.xlsx'))
    return 'rs2_60b.xlsx'


def rs2_60c():
    """RS2 #60 case 3: beta = 45 deg, vendor sigma_ci = 4.37 kPa (.fez #060_03)."""
    sd = _rs2_60_slope_data(45.0, 4.37)
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


# ======================================================================================
# RS2 #67 (Part III) — Stability of an earth dam under steady & transient unsaturated
# seepage, after Huang, M. & Jia, C-Q. (2009), "Strength reduction FEM in stability
# analysis of soil slopes subjected to transient unsaturated seepage."
#
# A HOMOGENEOUS earth dam (single Mohr-Coulomb material, Table 1: c = 13.8 kPa,
# phi = 37 deg, gamma = 18.2 kN/m3, E = 1e5 kPa, nu = 0.3), ~28 m tall on a ~191 m base,
# ~1V:3H upstream / ~1V:2.4H downstream with toe berms at el 6.66 / 6.86. The manual runs
# six SSR stages: Case 1 dry, Case 2 steady downstream, Case 3 (90 h after rapid drawdown)
# downstream + upstream, Case 4 (1500 h) downstream + upstream.
#
# XSLOPE has NO transient-seepage solver, so the drawdown pore-pressure fields are not
# solved — they are IMPORTED from RS2's own computed '.fea' nodal pore-pressure blocks
# through the u='seep' seep-solution path. benchmarks/rocscience/rs2_transient_seep.py
# reads a vendor computed '.fea' and writes the '*_mesh.json' / '*_seep.csv' sidecars
# (RS2's own mesh + RS2's own snapshot u), which are committed corpus artifacts like the
# vp077a_* seep sidecars. This verifies the SSR-under-transient-u MECHANICS; the transient
# seepage SOLUTION itself is RS2's, imported, not XSLOPE's.
#
# Only the stages whose vendor computed '.fea' carries a recoverable pore-pressure field
# are built: dry (01, u = 0 everywhere -> u='none'), 90 h downstream (03) and 90 h upstream
# (04) (each carries a full nodal block). The steady (02) and 1500 h (05/06) computed files
# carry NEITHER a nodal block NOR a phreatic-line geometry (empty gwgrid, zero material
# piezos), so those stages cannot be reconstructed from the snapshot and stay blocked
# (see rs2.md). The RS2 nodal field for each snapshot is EXACTLY hydrostatic below a
# 14-point phreatic surface (residual < 1e-3 kPa at every node) — RS2 represents each
# transient snapshot as a water-table position, not a spatially complex field.
#
# The 8-vertex dam outline is transcribed from the vendor mesh free-edge boundary and
# hard-coded so the '.xlsx' rebuilds without the vendor files; the circle is an inert
# placeholder (SSRM only). Case 3 downstream (03) runs an UNCONSTRAINED SSR (downstream is
# the weaker face); Case 3 upstream (04) confines strength reduction to RS2's upstream SSR
# Search Area (rectangle x in [-6.96, 102.32]) via the tag's ssr_zone, so the mechanism is
# the upstream face — the two stages share the same snapshot field and geometry.
# ======================================================================================

_RS2_67_DAM = [(0.256, 6.663), (0.266, 0.091), (191.582, 0.091), (191.582, 6.86),
               (157.459, 6.86), (107.161, 28.162), (99.272, 28.162), (33.787, 6.663)]


def _rs2_67_slope_data(u='none'):
    """One RS2 #67 stage: the homogeneous dam with pore-pressure option ``u`` ('none'
    for the dry case, 'seep' for a snapshot whose nodal field is imported through the
    committed '*_mesh.json' / '*_seep.csv' sidecars)."""
    return _poly_slope_data(
        polygons=[(0, _RS2_67_DAM)],
        materials=[dict(name='rock1', c=13.8, phi=37.0, gamma=18.2, gamma_sat=18.2,
                        E=1.0e5, nu=0.3, u=u)],
        circle={'Xo': 130.0, 'Yo': 34.0, 'Depth': 4.0, 'R': 30.0},
        max_depth=0.0)


def rs2_67a():
    """RS2 #67 Case 1 — dry dam (no free water surface), u='none'. Published: RS2 SSR
    2.48 | Slide2 Bishop 2.45 / Janbu 2.32 / Spencer 2.44 / GLE-MP 2.42 | ref LEM 2.43 /
    FEM 2.50."""
    sd = _rs2_67_slope_data(u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_67a.xlsx'))
    return 'rs2_67a.xlsx'


def rs2_67c():
    """RS2 #67 Case 3 downstream — free surface 90 h after rapid drawdown. Pore pressure
    IMPORTED from the vendor '#067_03 (90h).fea' nodal field via u='seep'
    (rs2_67c_mesh.json / rs2_67c_seep.csv). Unconstrained SSR — downstream is the weaker
    face. Published: RS2 SSR 1.83 | Slide2 Bishop 1.77 / Janbu 1.68 / Spencer 1.88 /
    GLE-MP 1.85 | ref LEM 1.92 / FEM 2.08."""
    sd = _rs2_67_slope_data(u='seep')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_67c.xlsx'))
    return 'rs2_67c.xlsx'


def rs2_67d():
    """RS2 #67 Case 3 upstream — same 90 h snapshot field as rs2_67c; the tag confines
    the SSR to RS2's upstream Search Area (rectangle x in [-6.96, 102.32]) so the critical
    mechanism is the upstream face. Pore pressure IMPORTED from '#067_04 (90h).fea'
    (rs2_67d_mesh.json / rs2_67d_seep.csv). Published: RS2 SSR 2.04 | Slide2 Bishop 1.99 /
    Janbu 1.89 / Spencer 2.07 / GLE-MP 2.06 | ref LEM 2.03."""
    sd = _rs2_67_slope_data(u='seep')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_67d.xlsx'))
    return 'rs2_67d.xlsx'


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


def rs2_68a():
    """RS2 #68 Case 1 (Ru = 0.5, homogeneous) -- Stability of Seismically Loaded
    Slopes, after Loukidis, Bandini & Salgado (2003), Geotechnique 53(5), 463-479.

    The verification target is a CRITICAL seismic coefficient k_c: the horizontal
    pseudo-static coefficient at which the slope is just stable (searched minimum
    FS = 1). It is not an FS -- the `critical_kc` test type bisects k until the
    circular search returns FS = 1 and locks the resulting k_c. The builder must
    therefore NOT set k_seismic; the harness sweeps it.

    A 25 m high homogeneous slope, 1V:3H face (toe (35, 25) -> crest (110, 50)),
    c = 25 kPa, phi = 30 deg, gamma = 20 kN/m3. Case 1 carries a pore-pressure
    ratio r_u = 0.5 (u = r_u * sigma_v) -- the manual's Case-1 condition, which the
    RS2 SSR .fez does not carry through the reader (it imports r_u = 0), so it is
    set here from the manual. Left-facing (ground descends leftward); the seismic
    auto-direction drives the failure.

    Published k_c (Slide2 LEM columns are the XSLOPE targets): Slide2 Bishop 0.118,
    Spencer 0.132; RS2 SSR 0.125; reference FEM 0.132 / Bishop 0.127 / Spencer 0.131
    / log-spiral 0.132; limit-analysis UB 0.145 / LB 0.126. Geometry transcribed
    from the RS2 vendor model 'slope stability #068_01.fez'."""
    poly = [(0.0, 0.0), (0.0, 25.0), (35.0, 25.0), (110.0, 50.0), (300.0, 50.0),
            (300.0, 0.0)]
    sd = _poly_slope_data(
        polygons=[(0, poly)],
        materials=[dict(name='rock1', c=25.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                        E=1.0e5, nu=0.3, u='ru', ru=0.5)],
        circle={'Xo': 72.0, 'Yo': 85.0, 'Depth': 20.0, 'R': 65.0},
        max_depth=0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_68a.xlsx'))
    return 'rs2_68a.xlsx'


def rs2_68b():
    """RS2 #68 Case 2 (Ru = 0, homogeneous) -- same 25 m, 1V:3H, c = 25 / phi = 30
    / gamma = 20 slope as Case 1 (rs2_68a) but dry (u = none). Critical seismic
    coefficient k_c is the verification target (see rs2_68a).

    Published k_c: Slide2 Bishop 0.425, Spencer 0.431; RS2 SSR 0.413; reference FEM
    0.433 / Bishop 0.426 / Spencer 0.431 / log-spiral 0.432; UB 0.454 / LB 0.423.
    Geometry from 'slope stability #068_02.fez' (200 m wide domain vs Case 1's
    300 m; the critical mechanism is near the slope, so the extra width is inert)."""
    poly = [(0.0, 0.0), (0.0, 25.0), (35.0, 25.0), (110.0, 50.0), (200.0, 50.0),
            (200.0, 0.0)]
    sd = _poly_slope_data(
        polygons=[(0, poly)],
        materials=[dict(name='rock1', c=25.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                        E=1.0e5, nu=0.3, u='none')],
        circle={'Xo': 72.0, 'Yo': 85.0, 'Depth': 20.0, 'R': 65.0},
        max_depth=0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_68b.xlsx'))
    return 'rs2_68b.xlsx'


def rs2_68c():
    """RS2 #68 Case 3 (Ru = 0, non-homogeneous 3-layer) -- a benched slope in three
    dipping rock bands, dry. Critical seismic coefficient k_c is the target (see
    rs2_68a). Materials (mat_id 0/1/2): rock1 c = 4 / phi = 30 / gamma = 17 (upper
    wedge), rock2 c = 25 / phi = 15 / gamma = 19 (the weak-friction middle band that
    the mechanism rides), rock3 c = 15 / phi = 45 / gamma = 19 (strong base). Ground
    (-30, 20) -> (150, 55); left-facing.

    Published k_c: Slide2 Bishop 0.155, Spencer 0.151; RS2 SSR 0.161; reference FEM
    0.161 / Bishop 0.155 / UB 0.172 / LB 0.148. Zone polygons transcribed from the
    RS2 vendor model 'slope stability #068_03.fez' via the .fez zone polygonizer."""
    zones = [
        (0, [(74.693, 41.077), (109.5, 55.0), (150.0, 55.0), (150.0, 44.5),
             (74.693, 41.077)]),                                        # rock1 (upper)
        (1, [(150.0, 44.5), (150.0, 33.0), (35.6, 27.8), (60.0, 40.0), (72.0, 40.0),
             (74.693, 41.077)]),                                        # rock2 (weak band)
        (2, [(150.0, 33.0), (150.0, -20.0), (-30.0, -20.0), (-30.0, 20.0),
             (20.0, 20.0), (35.6, 27.8)]),                              # rock3 (base)
    ]
    sd = _poly_slope_data(
        polygons=zones,
        materials=[
            dict(name='rock1', c=4.0, phi=30.0, gamma=17.0, gamma_sat=17.0,
                 E=1.0e5, nu=0.3, u='none'),
            dict(name='rock2', c=25.0, phi=15.0, gamma=19.0, gamma_sat=19.0,
                 E=1.0e5, nu=0.3, u='none'),
            dict(name='rock3', c=15.0, phi=45.0, gamma=19.0, gamma_sat=19.0,
                 E=1.0e5, nu=0.3, u='none'),
        ],
        circle={'Xo': 78.0, 'Yo': 82.0, 'Depth': 27.0, 'R': 55.0},
        max_depth=0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_68c.xlsx'))
    return 'rs2_68c.xlsx'


# ======================================================================================
# RS2 #64 -- Slope Stability Assessment of Three Homogeneous Landslides (Ankara E90),
# after Teoman, M.B., Topal, T. & Isik, N.S. (2004), and the RS2 Slope Stability
# Verification Manual Part III Problem 64 (pp. 219-227).
#
# Twelve cases: three Ankara-clay road-cut landslides (Slopes 1/2/3), each in its
# ORIGINAL and its FAILED (post-slide, with scarp) geometry, evaluated SHORT-TERM
# (total-stress, dry, no seismic; cases 1-6) and LONG-TERM (effective-stress with a
# piezometric line + a 0.03 g horizontal pseudo-static coefficient; cases 7-12).
# Each model is a SINGLE homogeneous Mohr-Coulomb zone; strengths are Tables 1-2 and
# are read straight from the vendor '.fez'. The manual reports three FS columns: RS2
# SSR, the Teoman reference (Bishop on a digitized PROPOSED slip surface), and Slide2
# Bishop on that same surface. The proposed surface exists in NO vendor file -- the
# '.fez' carry only an SSR Search-Area region -- so XSLOPE can reproduce ONLY the RS2
# SSR column (an SSRM run); the Ref/Slide2 Bishop columns are recorded as context.
#
# Elastic constants are the corpus convention (E = 14000 kPa, nu = 0.3, psi = 0 --
# Griffiths non-associated): the '.fez' reader imports E/nu = 0. The 0.03 g seismic is
# stored in the '.fez' as a per-element 'seismic forces: bx 0.03' body force that the
# reader does not extract (it returns k_seismic = 0), so it is set here by hand. For
# these LEFT-high slopes (crest at x = 0, toe descending to +x) RS2's bx = +0.03
# points downslope (+x = destabilizing); XSLOPE's FEM applies +k_seismic in +x, so the
# sign carries over directly (k_seismic = +0.03). Geometry (external boundary) and the
# piezo line are transcribed from 'slope stability #064_NN.fez' via read_fez (dev-side
# cross-check, hard-coded so the corpus rebuilds without the vendor files); the circle
# is an inert placeholder the loader requires (SSRM only -- no LEM surface is defined).
# ======================================================================================

def _rs2_64_slope_data(poly, c, phi, gamma, piezo=None, k=0.0):
    """One RS2 #64 case: a single homogeneous Mohr-Coulomb zone. ``poly`` is the
    external boundary; ``piezo`` (long-term cases) a piezometric line -> u='piezo';
    ``k`` the horizontal seismic coefficient (long-term = 0.03, set by hand)."""
    u = 'piezo' if piezo else 'none'
    sd = _poly_slope_data(
        polygons=[(0, poly)],
        materials=[dict(name='soil', c=float(c), phi=float(phi), gamma=float(gamma),
                        gamma_sat=float(gamma), E=14000.0, nu=0.3, u=u)],
        circle={'Xo': 8.0, 'Yo': 16.0, 'Depth': -4.0, 'R': 20.0},
        max_depth=0.0,
        piezo_line=list(piezo) if piezo else None)
    if k:
        sd['k_seismic'] = float(k)
    return sd


# Short-term (dry, total-stress) external boundaries -- cases 1-6.
_RS2_64A = [(25.5, -4.0), (25.5, 0.16), (20.915, 0.16), (4.111, 7.0), (0.0, 7.0),
            (0.0, -4.0), (4.111, -4.0), (20.915, -4.0)]
_RS2_64B = [(29.5, -4.0), (29.52, 0.641), (24.488, 0.641), (22.045, 1.663),
            (19.42, 2.812), (8.23, 5.016), (4.115, 6.995), (0.001, 6.995),
            (0.001, -4.047)]
_RS2_64C = [(23.0, -4.0), (23.0, 0.041), (18.062, 0.041), (6.607, 5.67), (4.161, 6.0),
            (0.0, 6.0), (0.0, -4.0), (4.161, -4.0), (6.607, -4.0), (18.062, -4.0)]
_RS2_64D = [(23.0, -4.0), (23.0, 0.0), (18.0, 0.0), (16.245, 1.909), (12.963, 2.833),
            (11.475, 2.833), (9.526, 3.961), (3.678, 4.885), (2.755, 6.0), (0.0, 6.0),
            (0.0, -4.0), (2.755, -4.0), (3.678, -4.0), (9.526, -4.0), (11.475, -4.0),
            (12.963, -4.0), (16.245, -4.0), (18.0, -4.0)]
_RS2_64E = [(18.0, -4.0), (18.0, 0.0), (13.0, 0.0), (3.02, 4.585), (0.0, 4.585),
            (0.0, -4.0), (3.02, -4.0), (13.0, -4.0)]
_RS2_64F = [(20.0, -4.0), (20.0, 0.352), (14.698, 0.352), (12.844, 1.647),
            (8.653, 2.83), (6.554, 2.495), (3.579, 2.786), (3.105, 4.491),
            (0.0, 4.491), (0.0, -4.0)]

# Long-term (saturated + 0.03 g) external boundaries and piezo lines -- cases 7-12.
_RS2_64G = [(26.111, -4.0), (26.111, -0.06), (20.975, -0.06), (9.294, 4.529),
            (3.669, 6.739), (0.0, 6.739), (0.0, -4.0)]
_RS2_64G_PZ = [(0.0, 4.529), (9.294, 4.529), (20.975, -0.06), (26.111, -2.093)]
_RS2_64H = [(28.445, -4.0), (28.445, 0.65), (23.359, 0.65), (20.877, 0.93),
            (18.627, 2.68), (7.814, 5.055), (6.927, 5.5), (3.938, 7.0), (0.0, 7.0),
            (0.0, -4.0), (3.938, -4.0), (6.927, -4.0), (7.814, -4.0), (18.627, -4.0),
            (20.877, -4.0), (23.359, -4.0)]
_RS2_64H_PZ = [(0.0, 5.5), (6.927, 5.5), (7.814, 5.055), (18.627, 2.68),
               (20.877, 0.93), (28.445, -2.162)]
_RS2_64I = [(23.07, -4.0), (23.07, 0.08), (18.055, 0.08), (9.27, 4.532), (6.042, 6.167),
            (3.29, 6.45), (0.0, 6.45), (0.0, -4.0)]
_RS2_64I_PZ = [(0.0, 4.532), (9.27, 4.532), (18.055, 0.08), (23.07, -1.973)]
_RS2_64J = [(22.8, -4.0), (22.8, 0.055), (18.017, 0.055), (16.214, 1.602), (12.8, 2.79),
            (11.437, 2.79), (8.968, 3.892), (3.822, 4.922), (3.551, 4.922),
            (2.859, 6.0), (0.0, 6.0), (0.0, -4.0)]
_RS2_64J_PZ = [(0.0, 4.922), (3.822, 4.922), (8.968, 3.892), (11.437, 2.79),
               (12.8, 2.79), (16.214, 1.602), (18.017, 0.055), (22.8, -2.053)]
_RS2_64K = [(18.0, -4.0), (18.0, 0.0), (12.823, 0.0), (6.07, 3.357), (3.57, 4.6),
            (0.0, 4.6), (0.0, -4.0)]
_RS2_64K_PZ = [(0.0, 3.357), (6.07, 3.357), (12.823, 0.0), (18.0, -2.086)]
_RS2_64L = [(20.0, -4.0), (20.0, 0.566), (14.424, 0.566), (12.074, 1.634),
            (8.995, 2.457), (7.457, 2.868), (5.531, 2.457), (3.668, 2.457),
            (3.603, 2.745), (3.144, 4.8), (0.0, 4.8), (0.0, -4.0)]
_RS2_64L_PZ = [(0.0, 2.745), (3.603, 2.745), (3.668, 2.457), (5.531, 2.457),
               (8.995, 2.457), (12.074, 1.634), (14.424, 0.566), (20.0, -0.783)]

# --- Vendor MATERIAL PARTITION corridors (elastic-materials run option) -----------
# The RS2 #64 long-term-Failed sub-unity cases C8/C12 constrain their SSR two ways at
# once: an SSR_polygonal_zones search area AND a material partition — the Mohr-Coulomb
# material ('rock1') is placed only in a corridor hugging the proposed slip surface,
# while the rest is 'rock2' with "Plasticity Specifications: None" (linear-elastic,
# cannot yield). ssr_zone can only hold the outside at FULL STRENGTH, so for these
# sub-unity slopes the failing skin outside the corridor could not be suppressed and
# no critical SRF bracketed (see rs2.md). solve_ssrm's elastic_materials option now
# makes the outside genuinely elastic, replicating rock2. These two corridors are the
# per-element rock1 footprints of 'slope stability #064_08.fez' / '#064_12.fez',
# recovered via rs2_ssr_zones.read_mc_footprint (the vendor 'elements:' material tags),
# simplified (Douglas-Peucker 0.2) and snapped to the domain boundary, then hard-coded
# so the corpus rebuilds without the vendor files. The elastic outer zone is the domain
# MINUS the corridor (shapely difference), an identical-strength material named for
# elastic_materials at solve time. Do NOT edit the single-material rs2_64h/rs2_64l.
_RS2_64H_MC = [(8.532, 3.477), (8.192, 4.244), (7.805, 4.196), (7.385, 4.659),
               (7.814, 5.055), (10.217, 4.527), (10.219, 4.058), (12.086, 3.872),
               (12.206, 3.631), (12.62, 3.999), (12.752, 3.411), (13.832, 3.421),
               (13.881, 3.045), (14.622, 3.559), (15.269, 2.726), (15.371, 3.072),
               (17.458, 2.613), (17.619, 2.256), (17.961, 2.504), (18.627, 2.68),
               (19.094, 1.749), (19.27, 2.18), (20.877, 0.93), (22.532, 0.743),
               (22.398, 0.287), (21.917, 0.417), (20.937, -0.246), (20.808, 0.184),
               (19.806, -0.003), (19.743, 0.62), (18.17, 0.286), (18.145, 0.814),
               (17.165, 0.704), (16.961, 1.179), (15.96, 0.976), (14.865, 1.695),
               (14.473, 1.238), (14.233, 1.648), (11.87, 1.959), (11.695, 2.411),
               (11.228, 2.239), (11.026, 2.822), (10.763, 2.493), (9.797, 2.884),
               (9.857, 3.25)]
_RS2_64L_MC = [(7.716, 1.02), (7.608, 1.692), (7.17, 1.672), (7.457, 2.868),
               (8.995, 2.457), (9.611, 2.292), (10.155, 1.505), (12.323, 0.635),
               (13.249, 1.1), (14.424, 0.566), (15.663, 0.566), (15.746, 0.205),
               (14.68, -0.023), (14.799, -0.355), (13.802, -0.582), (12.351, -0.392),
               (12.338, -0.739), (11.976, -0.701), (10.905, -0.511), (10.953, -0.115),
               (9.299, -0.0), (9.125, 0.485), (8.786, 0.138), (8.123, 0.715),
               (8.228, 0.991)]


def _rs2_64_split_slope_data(ext, mc_poly, c, phi, g_mc, g_el, piezo, k=0.03,
                             min_frac=0.02):
    """One RS2 #64 long-term case as a two-material PARTITION: the Mohr-Coulomb
    material ('rock1') fills only the corridor ``mc_poly``; the elastic outer
    zone(s) ('rock2', 'rock2b', ...) are the domain MINUS the corridor, carrying
    IDENTICAL strength but designated pure linear elastic at solve time via
    ``solve_ssrm(elastic_materials=['rock2', ...])``. The outer region is split
    into simple tiling polygons by shapely difference so gmsh meshes one conforming
    domain (the corridor is snapped to the boundary, so no interior hole/sliver is
    produced). Any outer fragment below ``min_frac`` of the domain (a sub-mesh-scale
    skin) is dropped into the corridor. Unit weights are the vendor's per-zone values
    (rock1 = g_mc, rock2 = g_el)."""
    dom = Polygon(ext)
    corr = Polygon(mc_poly)
    diff = dom.difference(corr)
    pieces = list(diff.geoms) if diff.geom_type == 'MultiPolygon' else [diff]
    pieces = [p for p in pieces if p.area > min_frac * dom.area]

    def _ring(p):
        return [(round(x, 3), round(y, 3)) for x, y in list(p.exterior.coords)[:-1]]

    polygons = [(0, _ring(corr))]
    mats = [dict(name='rock1', c=float(c), phi=float(phi), gamma=float(g_mc),
                 gamma_sat=float(g_mc), E=14000.0, nu=0.3, u='piezo')]
    elastic_names = []
    for i, p in enumerate(pieces):
        nm = 'rock2' if len(pieces) == 1 else f'rock2{chr(ord("a") + i)}'
        polygons.append((i + 1, _ring(p)))
        mats.append(dict(name=nm, c=float(c), phi=float(phi), gamma=float(g_el),
                         gamma_sat=float(g_el), E=14000.0, nu=0.3, u='piezo'))
        elastic_names.append(nm)
    sd = _poly_slope_data(
        polygons=polygons, materials=mats,
        circle={'Xo': 8.0, 'Yo': 16.0, 'Depth': -4.0, 'R': 20.0}, max_depth=0.0,
        piezo_line=list(piezo))
    if k:
        sd['k_seismic'] = float(k)
    return sd


def rs2_64a():
    """RS2 #64 Case 1 -- Slope 1 ORIGINAL, short-term (dry). c = 40.9, phi = 40.2,
    gamma = 20.5. Published: RS2 SSR 5.14 | Ref 5.25 | Slide2 5.24."""
    sd = _rs2_64_slope_data(_RS2_64A, 40.9, 40.2, 20.5)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64a.xlsx'))
    return 'rs2_64a.xlsx'


def rs2_64b():
    """RS2 #64 Case 2 -- Slope 1 FAILED, short-term (dry). c = 27.8, phi = 34,
    gamma = 20.5. Published: RS2 SSR 6.10 | Ref 6.67 | Slide2 6.64."""
    sd = _rs2_64_slope_data(_RS2_64B, 27.8, 34.0, 20.5)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64b.xlsx'))
    return 'rs2_64b.xlsx'


def rs2_64c():
    """RS2 #64 Case 3 -- Slope 2 ORIGINAL, short-term (dry). c = 33.6, phi = 41.4,
    gamma = 20. Published: RS2 SSR 4.69 | Ref 4.87 | Slide2 4.89."""
    sd = _rs2_64_slope_data(_RS2_64C, 33.6, 41.4, 20.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64c.xlsx'))
    return 'rs2_64c.xlsx'


def rs2_64d():
    """RS2 #64 Case 4 -- Slope 2 FAILED, short-term (dry). c = 28.4, phi = 33,
    gamma = 20. Published: RS2 SSR 4.95 | Ref 5.32 | Slide2 5.32."""
    sd = _rs2_64_slope_data(_RS2_64D, 28.4, 33.0, 20.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64d.xlsx'))
    return 'rs2_64d.xlsx'


def rs2_64e():
    """RS2 #64 Case 5 -- Slope 3 ORIGINAL, short-term (dry). c = 33.6, phi = 41.4,
    gamma = 20. Published: RS2 SSR 5.47 | Ref 5.44 | Slide2 5.45."""
    sd = _rs2_64_slope_data(_RS2_64E, 33.6, 41.4, 20.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64e.xlsx'))
    return 'rs2_64e.xlsx'


def rs2_64f():
    """RS2 #64 Case 6 -- Slope 3 FAILED, short-term (dry). c = 28.4, phi = 33,
    gamma = 20. Published: RS2 SSR 6.97 | Ref 7.02 | Slide2 6.96."""
    sd = _rs2_64_slope_data(_RS2_64F, 28.4, 33.0, 20.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64f.xlsx'))
    return 'rs2_64f.xlsx'


def rs2_64g():
    """RS2 #64 Case 7 -- Slope 1 ORIGINAL, long-term (saturated + 0.03 g). c = 12,
    phi = 26, gamma = 20.5. Published: RS2 SSR 1.7 | Ref 1.79 | Slide2 1.68."""
    sd = _rs2_64_slope_data(_RS2_64G, 12.0, 26.0, 20.5, piezo=_RS2_64G_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64g.xlsx'))
    return 'rs2_64g.xlsx'


def rs2_64h():
    """RS2 #64 Case 8 -- Slope 1 FAILED, long-term (saturated + 0.03 g). c = 3,
    phi = 19, gamma = 20.5. Published: RS2 SSR 0.99 | Ref 1.13 | Slide2 1.09."""
    sd = _rs2_64_slope_data(_RS2_64H, 3.0, 19.0, 20.5, piezo=_RS2_64H_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64h.xlsx'))
    return 'rs2_64h.xlsx'


def rs2_64i():
    """RS2 #64 Case 9 -- Slope 2 ORIGINAL, long-term (saturated + 0.03 g). c = 7,
    phi = 32, gamma = 20. Published: RS2 SSR 1.30 | Ref 1.30 | Slide2 1.30."""
    sd = _rs2_64_slope_data(_RS2_64I, 7.0, 32.0, 20.0, piezo=_RS2_64I_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64i.xlsx'))
    return 'rs2_64i.xlsx'


def rs2_64j():
    """RS2 #64 Case 10 -- Slope 2 FAILED, long-term (saturated + 0.03 g). c = 4,
    phi = 25, gamma = 20. Published: RS2 SSR 1.09 | Ref 1.08 | Slide2 1.07."""
    sd = _rs2_64_slope_data(_RS2_64J, 4.0, 25.0, 20.0, piezo=_RS2_64J_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64j.xlsx'))
    return 'rs2_64j.xlsx'


def rs2_64k():
    """RS2 #64 Case 11 -- Slope 3 ORIGINAL, long-term (saturated + 0.03 g). c = 7,
    phi = 32, gamma = 20. Published: RS2 SSR 1.46 | Ref 1.51 | Slide2 1.51."""
    sd = _rs2_64_slope_data(_RS2_64K, 7.0, 32.0, 20.0, piezo=_RS2_64K_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64k.xlsx'))
    return 'rs2_64k.xlsx'


def rs2_64l():
    """RS2 #64 Case 12 -- Slope 3 FAILED, long-term (saturated + 0.03 g). c = 2,
    phi = 25, gamma = 20. Published: RS2 SSR 1.22 | Ref 1.13 | Slide2 1.15."""
    sd = _rs2_64_slope_data(_RS2_64L, 2.0, 25.0, 20.0, piezo=_RS2_64L_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64l.xlsx'))
    return 'rs2_64l.xlsx'


def rs2_64h_split():
    """RS2 #64 Case 8 (Slope 1 FAILED, long-term) as the VENDOR MATERIAL PARTITION:
    rock1 Mohr-Coulomb corridor (c=3, phi=19, gamma=20.5) + rock2 elastic-None outer
    (identical strength, gamma=20.0). Run constrained via
    solve_ssrm(elastic_materials=['rock2']). Published RS2 SSR 0.99. Separate from the
    single-material rs2_64h.xlsx (unchanged)."""
    sd = _rs2_64_split_slope_data(_RS2_64H, _RS2_64H_MC, 3.0, 19.0, 20.5, 20.0,
                                  _RS2_64H_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64h_split.xlsx'))
    return 'rs2_64h_split.xlsx'


def rs2_64l_split():
    """RS2 #64 Case 12 (Slope 3 FAILED, long-term) as the VENDOR MATERIAL PARTITION:
    rock1 Mohr-Coulomb corridor (c=2, phi=25, gamma=20.0) + rock2 elastic-None outer
    (identical strength, gamma=20.0). Run constrained via
    solve_ssrm(elastic_materials=['rock2']). Published RS2 SSR 1.22. Separate from the
    single-material rs2_64l.xlsx (unchanged)."""
    sd = _rs2_64_split_slope_data(_RS2_64L, _RS2_64L_MC, 2.0, 25.0, 20.0, 20.0,
                                  _RS2_64L_PZ, k=0.03)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'rs2_64l_split.xlsx'))
    return 'rs2_64l_split.xlsx'


if __name__ == '__main__':
    for fn in (rs2_56a, rs2_56b, rs2_57a, rs2_57b, rs2_58a, rs2_58b, hammah_hb1,
               rs2_60a, rs2_60b, rs2_60c, rs2_61a, rs2_59, rs2_63,
               rs2_66a, rs2_66b, rs2_66c, rs2_66d, rs2_66e,
               rs2_62a, rs2_62b, rs2_62c, rs2_65, rs2_51,
               rs2_67a, rs2_67c, rs2_67d,
               rs2_68a, rs2_68b, rs2_68c,
               rs2_64a, rs2_64b, rs2_64c, rs2_64d, rs2_64e, rs2_64f,
               rs2_64g, rs2_64h, rs2_64i, rs2_64j, rs2_64k, rs2_64l,
               rs2_64h_split, rs2_64l_split):
        print(fn())
