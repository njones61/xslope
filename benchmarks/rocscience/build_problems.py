"""Builders for the Rocscience Slide2 verification corpus (docs/verification.md).

Each vpNNN() function reproduces one manual problem as an XSLOPE input file in
docs/files/rocscience/, from the manual's tabulated data and coordinate-labeled
figures (see the corpus section of docs/verification.md for the methodology and
status table). Run this script to regenerate every built problem.

Extraction pipeline for new problems:
  1. text: pdftotext -layout Slide_SlopeStabilityVerification.pdf (properties,
     search settings, results tables)
  2. figures: pdftoppm -png -r 150 -f <page> -l <page> ... then read the
     coordinate labels off the geometry figure (printed page == pdf page)
  3. build via xslope.fileio.save_slope_data_to_xlsx, verify against the
     manual's Slide2/reference values, add a test tag in docs/verification.md
"""

import math
import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'files', 'rocscience')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')


def vp002():
    """ACADS 1(b) (Giam & Donald 1989): the 1(a) slope with c'=32, phi'=10,
    gamma=20 and a water-filled tension crack of depth 2c/(gamma*sqrt(ka))
    [Craig 1997]. Slide2: Bishop 1.596, Spencer 1.592, Janbu corrected 1.489,
    GLE 1.592; Giam reference 1.65."""
    ka = (1 - math.sin(math.radians(10.0))) / (1 + math.sin(math.radians(10.0)))
    depth = 2 * 32.0 / (20.0 * math.sqrt(ka))     # 3.8136 m
    sd = load_slope_data(ACADS_1A)                 # ACADS 1(a) geometry
    m = sd['materials'][0]
    m.update(c=32.0, phi=10.0, gamma=20.0, name='ACADS 1(b) soil')
    sd['tcrack_depth'] = round(depth, 3)
    sd['tcrack_water'] = round(depth, 3)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp002.xlsx'))
    return 'vp002.xlsx'




LEVEE_POLY = os.path.join(os.path.dirname(__file__), '..', '..',
                          'docs', 'seep', 'files', 'xslope_levee_poly.xlsx')

# Talbingo Dam (ACADS 2(a)/2(b)) — Table 5.2 point coordinates, connectivity
# read from Figures 5.1/5.2 (upstream shell | transition | core | filter seam |
# downstream transition | downstream shell).
_T = {1: (0, 0), 2: (315.5, 162), 3: (319.5, 162), 4: (321.6, 162), 5: (327.6, 162),
      6: (386.9, 130.6), 7: (394.1, 130.6), 8: (453.4, 97.9), 9: (460.6, 97.9),
      10: (515, 65.3), 11: (521.1, 65.3), 12: (577.9, 31.4), 13: (585.1, 31.4),
      14: (648, 0), 15: (168.1, 0), 16: (302.2, 130.6), 17: (200.7, 0),
      18: (311.9, 130.6), 19: (307.1, 0), 20: (331.3, 130.6), 21: (328.8, 146.1),
      22: (310.7, 0), 23: (333.7, 130.6), 24: (331.3, 146.1), 25: (372.4, 0),
      26: (347, 130.6)}

_TALBINGO_ZONES = [  # (mat_id, point numbers, ccw or cw — shapely doesn't care)
    (0, [1, 2, 16, 15]),                          # upstream rockfill shell
    (1, [2, 3, 18, 17, 15, 16]),                  # upstream transition
    (3, [3, 4, 21, 20, 19, 17, 18]),              # core (inclined)
    (2, [4, 5, 24, 23, 22, 19, 20, 21]),          # filter (very thin seam)
    (1, [5, 26, 25, 22, 23, 24]),                 # downstream transition
    (0, [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25, 26]),  # downstream rockfill shell
]


def _talbingo_slope_data():
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    sd = load_slope_data(LEVEE_POLY)   # polygon-mode base
    base_mat = dict(sd['materials'][0])
    props = [('Rockfill', 0.0, 45.0, 20.4), ('Transitions', 0.0, 45.0, 20.4),
             ('Filter', 0.0, 45.0, 20.4), ('Core', 85.0, 23.0, 18.1)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base_mat)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='none')
        sd['materials'].append(m)
    sd['polygons'] = [{'polygon': Polygon([_T[i] for i in pts]), 'mat_id': mid}
                      for mid, pts in _TALBINGO_ZONES]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['circular'] = True
    return sd


def vp005():
    """ACADS 2(a) (Giam & Donald 1989): Talbingo Dam at end of construction,
    4 zones, critical circular surface. Slide2: Bishop 1.948, Spencer 1.948,
    GLE 1.948, Janbu corrected 1.949; Giam reference 1.95. The minimum is a
    shallow slide parallel to the (steeper) upstream face."""
    sd = _talbingo_slope_data()
    sd['circles'] = [{'Xo': 100.0, 'Yo': 290.0, 'Depth': 12.0, 'R': 278.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp005.xlsx'))
    return 'vp005.xlsx'


def vp006():
    """ACADS 2(b): Talbingo Dam, single specified circle Xc=100.3, Yc=291.0,
    R=278.8 (Table 6.1). Slide2: Bishop 2.208, Spencer 2.292, GLE 2.301,
    Janbu corrected 2.073; Giam reference 2.29."""
    sd = _talbingo_slope_data()
    sd['circles'] = [{'Xo': 100.3, 'Yo': 291.0, 'Depth': 291.0 - 278.8, 'R': 278.8}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp006.xlsx'))
    return 'vp006.xlsx'


def _acads_nonhom_slope_data(k_seismic=0.0):
    """ACADS 1(c)/1(d) three-layer slope. Geometry from the GeoStudio
    verification manual Figure 7 (fully coordinate-labeled); identical problem
    in the Slide2 manual (#3/#4) whose figure has unlabeled interface nodes.
    Domain extended to x=0..90 (interfaces flat beyond their last point) to
    give the circular search room, matching the ACADS 1(a) sample extents."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Soil #1', 0.0, 38.0), ('Soil #2', 5.3, 23.0), ('Soil #3', 7.2, 20.0)]
    sd['materials'] = []
    for name, c, phi in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=19.5, option='mc', u='none')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 25.0), (30.0, 25.0), (50.0, 35.0), (90.0, 35.0)]},
        {'mat_id': 1, 'coords': [(0.0, 25.0), (30.0, 25.0), (50.0, 29.0), (54.0, 31.0), (90.0, 31.0)]},
        # Soil 3's surface coincides with soil 1's base from x=30 to the (40,27)
        # hump (soil 2 pinches out there), then drops to (52,24) — see the Slide
        # manual Figure 3.1 axes and the GeoStudio Figure 7 labels.
        {'mat_id': 2, 'coords': [(0.0, 25.0), (30.0, 25.0), (40.0, 27.0), (52.0, 24.0), (90.0, 24.0)]},
    ]
    sd['max_depth'] = 20.0
    sd['k_seismic'] = k_seismic
    return sd


def vp003():
    """ACADS 1(c): non-homogeneous three-layer slope, critical circle.
    Slide2: Bishop 1.405, Spencer 1.375, GLE 1.374, Janbu corrected 1.357;
    SLOPE/W: Bishop 1.414, M-P 1.382; ACADS reference 1.39."""
    sd = _acads_nonhom_slope_data()
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp003.xlsx'))
    return 'vp003.xlsx'


def vp004():
    """ACADS 1(d): problem #3 plus horizontal seismic coefficient 0.15.
    Slide2: Bishop 1.016, Spencer 0.991, GLE 0.989, Janbu corrected 0.965;
    SLOPE/W: Bishop 1.02, M-P 0.989; ACADS reference 1.00."""
    sd = _acads_nonhom_slope_data(k_seismic=0.15)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp004.xlsx'))
    return 'vp004.xlsx'


WEAK_LAYER = os.path.join(os.path.dirname(__file__), '..', '..',
                          'docs', 'lem', 'files', 'xslope_acads_weak_layer.xlsx')


def vp008():
    """ACADS 3(b): the weak-layer slope (= LEM sample 13 / Slide #7) with the
    fully specified non-circular surface of Table 8.2. Slide2: Spencer 1.277,
    GLE 1.262, Janbu corrected 1.294; SLOPE/W: Bishop 1.259, M-P 1.261;
    Giam reference 1.34."""
    sd = load_slope_data(WEAK_LAYER)
    pts = [(41.85, 27.75, 'Free'), (44.0, 26.5, 'Horiz'),
           (63.5, 27.0, 'Horiz'), (73.31, 40.0, 'Free')]
    sd['non_circ'] = [{'X': x, 'Y': y, 'Movement': m} for x, y, m in pts]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp008.xlsx'))
    return 'vp008.xlsx'


def vp009():
    """ACADS 4 (Slide #9): weak-layer slope + piezometric surface (Table 9.3)
    + two surcharge strips (Table 9.2: 20 kPa on the lower bench x=23-43,
    20->40 kPa ramp on the crest x=70-80). Non-circular search. Slide2 (block
    search, no optimization): Spencer 0.760, GLE 0.720, Janbu corrected 0.734;
    with optimization 0.683-0.707; SLOPE/W: Bishop 0.699, M-P 0.689; Giam
    reference 0.78; Slope 2000 GLE reference 0.6878. The published spread is
    wide - this is a search-difficulty benchmark."""
    sd = load_slope_data(WEAK_LAYER)
    # ACADS 4's weak layer is NOT problem 3's horizontal seam: it is a 0.6 m
    # band inclined roughly parallel to the face, (20,18.28-18.88) up to
    # (84,36.2-36.8) - read from the labeled GeoStudio Figure 25 (the Slide
    # figure shows the same diagonal seam). The bedrock floor drops to y=15.
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(20.0, 27.75), (43.0, 27.75), (67.5, 40.0), (84.0, 40.0)]},
        {'mat_id': 1, 'coords': [(20.0, 18.88), (84.0, 36.8)]},
        {'mat_id': 0, 'coords': [(20.0, 18.28), (84.0, 36.2)]},
    ]
    sd['max_depth'] = 15.0
    sd['piezo_line'] = [(20.0, 27.75), (43.0, 27.75), (49.0, 29.8), (60.0, 34.0),
                        (66.0, 35.8), (74.0, 37.6), (80.0, 38.4), (84.0, 38.4)]
    sd['dloads'] = [
        [{'X': 23.0, 'Y': 27.75, 'Normal': 20.0}, {'X': 43.0, 'Y': 27.75, 'Normal': 20.0}],
        [{'X': 70.0, 'Y': 40.0, 'Normal': 20.0}, {'X': 80.0, 'Y': 40.0, 'Normal': 40.0}],
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp009.xlsx'))
    return 'vp009.xlsx'


def _fk_slope_data():
    """Fredlund & Krahn (1977) homogeneous slope (Slide #21), imperial units:
    ground (0,60)-(60,60)-(140,20)-(180,20), bedrock at y=0, c'=600 psf,
    phi'=20, gamma=120 pcf, specified circle xc=120, yc=90, R=80."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='F&K soil', c=600.0, phi=20.0, gamma=120.0, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 60.0), (60.0, 60.0), (140.0, 20.0), (180.0, 20.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circles'] = [{'Xo': 120.0, 'Yo': 90.0, 'Depth': 10.0, 'R': 80.0}]
    return sd


def vp021a():
    """Slide #21 case 1 (dry). F&K: OMS 1.928, Bishop 2.080, Spencer 2.073,
    M-P 2.076; Slide: 1.931 / 2.079 / 2.075 / 2.075."""
    sd = _fk_slope_data()
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp021a.xlsx'))
    return 'vp021a.xlsx'


def vp021b():
    """Slide #21 case 2 (ru = 0.25). F&K: OMS 1.607, Bishop 1.766, Spencer
    1.761, M-P 1.764; Slide: 1.687 / 1.763 / 1.760 / 1.760. Exercises the
    template v12 ru pore-pressure option."""
    sd = _fk_slope_data()
    sd['materials'][0].update(u='ru', ru=0.25)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp021b.xlsx'))
    return 'vp021b.xlsx'


def vp018():
    """Slide #18: Spencer (1969) / Baker (1980) homogeneous slope with ru=0.5,
    non-circular critical surface. The slope descends left-to-right (a
    right-facing case). Slide2 (random search + Monte Carlo optimization):
    Spencer 1.010; Baker reference 1.02; Spencer (1969) 1.08."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Spencer 1969 soil', c=10.8, phi=40.0, gamma=18.0,
             option='mc', u='ru', ru=0.5)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 40.0), (10.0, 40.0), (70.0, 10.0), (80.0, 10.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = []
    sd['circular'] = False
    sd['non_circ'] = [
        {'X': 12.0, 'Y': 39.0, 'Movement': 'Free'},
        {'X': 25.0, 'Y': 24.0, 'Movement': 'Horiz'},
        {'X': 45.0, 'Y': 13.0, 'Movement': 'Horiz'},
        {'X': 60.0, 'Y': 10.5, 'Movement': 'Horiz'},
        {'X': 70.0, 'Y': 10.0, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp018.xlsx'))
    return 'vp018.xlsx'


def vp015():
    """Slide #15: Arai & Tagyo (1985) example 2 - three layers with a weak
    (c=9.8, phi=5) middle band, no water. Circular search. Slide2 (auto
    refine): Bishop 0.420, Spencer 0.409, GLE 0.437, Janbu corrected 0.423;
    A&T Bishop 0.417; Kim et al. 0.43."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Upper Layer', 29.4, 12.0), ('Middle Layer', 9.8, 5.0),
             ('Lower Layer', 294.0, 40.0)]
    sd['materials'] = []
    for name, c, phi in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=18.82, option='mc', u='none')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 15.0), (18.0, 15.0), (48.0, 35.0), (96.0, 35.0)]},
        {'mat_id': 1, 'coords': [(0.0, 15.0), (18.0, 15.0), (24.0, 19.0), (72.0, 35.0), (96.0, 35.0)]},
        {'mat_id': 2, 'coords': [(0.0, 3.0), (96.0, 35.0)]},
    ]
    sd['max_depth'] = 3.0
    sd['circles'] = [{'Xo': 30.0, 'Yo': 40.0, 'Depth': 12.0, 'R': 28.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp015.xlsx'))
    return 'vp015.xlsx'


def vp016():
    """Slide #16: Arai & Tagyo (1985) example 3 - homogeneous slope with a
    water table. Circular search. Slide2 (auto refine): Bishop 1.118,
    Janbu simplified 1.046, Janbu corrected 1.131, Spencer 1.118;
    A&T Bishop 1.138."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='A&T ex3 soil', c=41.65, phi=15.0, gamma=18.82,
             option='mc', u='piezo')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 15.0), (18.0, 15.0), (48.0, 35.0), (66.0, 35.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['piezo_line'] = [(0.0, 15.0), (18.0, 15.0), (30.0, 23.0), (48.0, 29.0), (66.0, 32.0)]
    sd['circles'] = [{'Xo': 30.0, 'Yo': 45.0, 'Depth': 10.0, 'R': 35.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp016.xlsx'))
    return 'vp016.xlsx'


def vp019():
    """Slide #19: Greco (1996) ex. 4 / Yamagami & Ueta (1988) four-layer
    slope, no water, non-circular critical surface. Slide2 (random search +
    Monte-Carlo optimization, convex): Spencer 1.398; Greco and Yamagami &
    Ueta references 1.40-1.42."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Upper Layer', 49.0, 29.0, 20.38), ('Layer 2', 0.0, 30.0, 17.64),
             ('Layer 3', 7.84, 20.0, 20.38), ('Bottom Layer', 0.0, 30.0, 17.64)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='none')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 40.0), (60.0, 40.0), (180.0, 100.0), (200.0, 100.0), (260.0, 70.0)]},
        {'mat_id': 1, 'coords': [(0.0, 40.0), (260.0, 40.0)]},
        {'mat_id': 2, 'coords': [(0.0, 35.0), (260.0, 35.0)]},
        {'mat_id': 3, 'coords': [(0.0, 26.0), (260.0, 26.0)]},
    ]
    sd['max_depth'] = 0.0
    # Both search families are stored: the circular minimum (Spencer 1.429)
    # lands within the published 1.40-1.42 range; the local non-circular
    # search plateaus ~1.45 from this seed, above Slide's Monte-Carlo 1.398 -
    # a search-power difference, not a model difference.
    sd['circular'] = True
    sd['circles'] = [{'Xo': 86.0, 'Yo': 134.0, 'Depth': 30.0, 'R': 104.0}]
    sd['non_circ'] = [
        {'X': 40.0, 'Y': 40.0, 'Movement': 'Free'},
        {'X': 70.0, 'Y': 28.0, 'Movement': 'Horiz'},
        {'X': 100.0, 'Y': 28.0, 'Movement': 'Horiz'},
        {'X': 130.0, 'Y': 28.0, 'Movement': 'Horiz'},
        {'X': 160.0, 'Y': 45.0, 'Movement': 'Horiz'},
        {'X': 180.0, 'Y': 75.0, 'Movement': 'Horiz'},
        {'X': 196.0, 'Y': 100.0, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp019.xlsx'))
    return 'vp019.xlsx'


def vp020():
    """Slide #20: Greco (1996) ex. 5 / Chen & Shao (1988): four layers with a
    0.5 m weak seam along the inclined model base, water table. Polygon-zone
    geometry (the base is not horizontal). Slide2: circular (toe focus grid)
    Bishop 1.087, Spencer 1.093; non-circular (block search in seam + MC
    optimization) Spencer 1.010; Chen & Shao 1.01-1.03; Greco 0.973-1.1."""
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    sd = load_slope_data(LEVEE_POLY)   # polygon-mode base
    base_mat = dict(sd['materials'][0])
    props = [('Layer 1', 9.8, 35.0, 20.0), ('Layer 2', 58.8, 25.0, 19.0),
             ('Layer 3', 19.8, 30.0, 21.5), ('Layer 4 (seam)', 9.8, 16.0, 21.5)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base_mat)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='piezo')
        sd['materials'].append(m)
    zones = [
        (0, [(100,40),(145,70),(240,70),(240,55),(190,55)]),
        (1, [(75,30),(95,40),(100,40),(190,55),(240,55),(240,30)]),
        (2, [(0,20),(55,20),(75,30),(240,30),(240,16.5),(150,15.5),(0,0.5)]),
        (3, [(0,0.5),(150,15.5),(240,16.5),(240,16),(150,15),(0,0)]),
    ]
    sd['polygons'] = [{'polygon': Polygon(pts), 'mat_id': mid} for mid, pts in zones]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['piezo_line'] = [(0.0,20.0),(55.0,20.0),(75.0,30.0),(95.0,40.0),(100.0,40.0),(190.0,55.0),(240.0,55.0)]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 90.0, 'Yo': 60.0, 'Depth': 15.0, 'R': 45.0}]
    # Non-circular seed for the seam-block mechanism: the local search reaches
    # ~1.08 from here (below the circular minimum); Slide's Monte-Carlo block
    # optimization reaches 1.010 - same search-power gap noted on #19.
    sd['non_circ'] = [
        {'X': 60.0, 'Y': 22.5, 'Movement': 'Free'},
        {'X': 75.0, 'Y': 8.1, 'Movement': 'Horiz'},
        {'X': 140.0, 'Y': 14.7, 'Movement': 'Horiz'},
        {'X': 190.0, 'Y': 55.0, 'Movement': 'Horiz'},
        {'X': 210.0, 'Y': 70.0, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp020.xlsx'))
    return 'vp020.xlsx'


def vp017():
    """Slide #17: Yamagami & Ueta (1988) homogeneous slope, dry. Circular:
    Slide Bishop 1.344, Ordinary 1.278 (Y&U 1.348 / 1.282). Non-circular:
    Slide Spencer 1.325 (Y&U 1.339, Greco 1.33)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Y&U soil', c=9.8, phi=10.0, gamma=17.64, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 5.0), (5.0, 5.0), (15.0, 10.0), (25.0, 10.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 10.0, 'Yo': 14.0, 'Depth': 3.0, 'R': 11.0}]
    sd['non_circ'] = [
        {'X': 4.5, 'Y': 5.0, 'Movement': 'Free'},
        {'X': 8.0, 'Y': 4.3, 'Movement': 'Horiz'},
        {'X': 12.0, 'Y': 5.0, 'Movement': 'Horiz'},
        {'X': 16.5, 'Y': 10.0, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp017.xlsx'))
    return 'vp017.xlsx'


def vp024():
    """Slide #24: Low (1989) three-layer undrained slope (phi=0). Circular
    search. Slide2: Ordinary 1.439, Bishop 1.439; Low reference 1.44 both."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Upper Layer', 30.0), ('Middle Layer', 20.0), ('Bottom Layer', 150.0)]
    sd['materials'] = []
    for name, c in props:
        m = dict(base)
        m.update(name=name, c=c, phi=0.0, gamma=18.0, option='mc', u='none')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 14.0), (20.0, 14.0), (34.0, 9.0), (38.0, 8.0), (60.0, 8.0)]},
        {'mat_id': 1, 'coords': [(0.0, 9.0), (34.0, 9.0), (38.0, 8.0), (60.0, 8.0)]},
        {'mat_id': 2, 'coords': [(0.0, 5.0), (60.0, 5.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 28.0, 'Yo': 20.0, 'Depth': 5.0, 'R': 15.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp024.xlsx'))
    return 'vp024.xlsx'


def vp023():
    """Slide #23: Low (1989) slope over two undrained layers; the lower
    layer's cu grows linearly 15->30 kPa from y=4 to y=0 (xslope 'cp'
    option: Su = c + cp*(r_elev - y)). Circular search. Slide2:
    Ordinary 1.370, Bishop 1.192; Low 1.36 / 1.14; Kim (2002) 1.17."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    sd['materials'] = []
    m = dict(base); m.update(name='Upper Soil', c=95.0, phi=15.0, gamma=20.0, option='mc', u='none')
    sd['materials'].append(m)
    m = dict(base); m.update(name='Middle Soil', c=15.0, phi=0.0, gamma=20.0, option='mc', u='none')
    sd['materials'].append(m)
    m = dict(base); m.update(name='Lower Soil', c=15.0, phi=0.0, gamma=20.0, option='cp',
                             cp=3.75, r_elev=4.0, u='none')
    sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 8.0), (10.0, 8.0), (26.0, 16.0), (40.0, 16.0)]},
        {'mat_id': 1, 'coords': [(0.0, 8.0), (10.0, 8.0), (40.0, 8.0)]},
        {'mat_id': 2, 'coords': [(0.0, 4.0), (40.0, 4.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 18.0, 'Yo': 20.0, 'Depth': 2.0, 'R': 18.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp023.xlsx'))
    return 'vp023.xlsx'


def _dw634_slope_data(appl):
    """Duncan & Wright (2005) Fig. 6.34 (Slide #85): 20-ft undrained clay
    slope (c=350 psf, phi=0, gamma=98 pcf), one horizontal 9,000 lb/ft
    tieback at mid-height, fully anchored (Lp=0 -> constant capacity)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='saturated clay', c=350.0, phi=0.0, gamma=98.0, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(15.0, 10.0), (25.0, 30.0), (57.0, 30.0)]},
    ]
    sd['max_depth'] = 10.0
    sd['gamma_water'] = 62.4
    # Slide's printed critical circles (GLE fig 85.2 for active, Bishop
    # fig 85.3 for passive) are stored by the vp085* builders below so the
    # regression tags evaluate a deterministic surface.
    sd['circles'] = [{'Xo': 30.0, 'Yo': 40.0, 'Depth': 12.0, 'R': 28.0}]
    sd['reinforcement_lines'] = [{
        'x1': 20.0, 'y1': 20.0, 'x2': 57.0, 'y2': 20.0,
        't_max': 9000.0, 't_res': 0.0, 'lp1': 0.0, 'lp2': 0.0,
        'E': float('nan'), 'area': float('nan'), 'label': 'Tieback',
        'type': 'anchor', 'dir': 'axial', 'appl': appl,
        'tend1': 0.0, 'tend2': 0.0, 'spacing': 1.0,
    }]
    sd['reinforce_lines'] = sd['reinforcement_lines']
    return sd


def vp085a():
    """Slide #85 case 1 (active). D&W reference 1.51; Slide circular Bishop
    1.531."""
    sd = _dw634_slope_data('active')
    sd['circles'] = [{'Xo': 15.446, 'Yo': 37.624, 'Depth': 37.624 - 27.594, 'R': 27.594}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp085a.xlsx'))
    return 'vp085a.xlsx'


def vp085b():
    """Slide #85 case 2 (passive). D&W reference 1.32; Slide circular Bishop
    1.324."""
    sd = _dw634_slope_data('passive')
    sd['circles'] = [{'Xo': 17.169, 'Yo': 34.480, 'Depth': 34.480 - 24.465, 'R': 24.465}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp085b.xlsx'))
    return 'vp085b.xlsx'


def _yamagami_pile_slope(with_pile):
    """Yamagami (2000) micro-pile slope (Slide #54): homogeneous c=4.9,
    phi=10, gamma=15.68; single vertical micro-pile row at x=9 (crest),
    10.7 kN shear per pile at 1 m spacing."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Material 1', c=4.9, phi=10.0, gamma=15.68, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-6.0, 0.0), (0.0, 0.0), (8.0, 4.0), (12.0, 4.0)]},
    ]
    sd['max_depth'] = -5.0
    sd['circles'] = [{'Xo': 2.674, 'Yo': 7.573, 'Depth': 7.573 - 8.031, 'R': 8.031}]
    if with_pile:
        sd['pile_lines'] = [{
            'x1': 9.0, 'y1': 4.0, 'x2': 9.0, 'y2': -2.0,
            'H': None, 'theta_p': 0.0, 'D_pile': 0.3, 'S': 1.0,
            'E': None, 'I': None, 'area': None,
            'V_cap': 10.7, 'M_cap': 1.0e6, 'fixity': 'free',   # shear governs; Mcap>0 required
            # Slide applies the micro-pile shear in the ACTIVE sense (added
            # to the resisting sum un-factored): active reproduces its 1.193
            # on the printed circle (1.185); passive-(/F) gives 1.172.
            'appl': 'active', 'label': 'micro-pile',
        }]
    return sd


def vp054a():
    """Slide #54, unreinforced case on the printed critical circle
    (2.674, 7.573, R=8.031). Slide Bishop 1.102; Yamagami 1.10."""
    sd = _yamagami_pile_slope(False)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp054a.xlsx'))
    return 'vp054a.xlsx'


def vp054b():
    """Slide #54 with the micro-pile row. Slide 1.193; Yamagami 1.20."""
    sd = _yamagami_pile_slope(True)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp054b.xlsx'))
    return 'vp054b.xlsx'


def vp041():
    """Slide #41: Jiang, Baker & Yamagami (2003) homogeneous clay slope with
    power-curve strength tau = 1.4*(sigma')^0.8 and ru = 0.3 - exercises the
    v12 pow option and ru together. Slide2: Bishop 1.656, Janbu simplified
    1.563; Charles & Soares 1.66; Baker 1.56-1.60; Perry 1.67."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='clay (power curve)', c=0.0, phi=0.0, gamma=20.0,
             option='pow', pow_a=1.4, pow_b=0.8, pow_c=0.0, pow_d=0.0,
             u='ru', ru=0.3)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 30.0), (5.0, 30.0), (85.0, 10.0), (93.0, 10.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 45.0, 'Yo': 45.0, 'Depth': 8.0, 'R': 37.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp041.xlsx'))
    return 'vp041.xlsx'


def _baker2_slope_data():
    """Baker (2003) example 2 (Slide #45): gentle 4:1 clay slope, 12 m high,
    ground (0,0)-(48,12)-(100,12), bedrock at y=0."""
    sd = load_slope_data(ACADS_1A)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (48.0, 12.0), (100.0, 12.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 25.0, 'Yo': 25.0, 'Depth': 2.0, 'R': 23.0}]
    return sd


def vp045a():
    """Slide #45, Mohr-Coulomb case: c'=11.64, phi'=24.7, gamma=18.
    Slide2: Janbu simplified 2.662, Spencer 2.794."""
    sd = _baker2_slope_data()
    sd['materials'][0].update(name='clay (MC)', c=11.64, phi=24.7, gamma=18.0,
                              option='mc', u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp045a.xlsx'))
    return 'vp045a.xlsx'


def vp045b():
    """Slide #45, power-curve case: tau = 1.107*(sigma')^0.86 (Baker's
    A=0.58, n=0.86, T=0). Slide2: Janbu simplified 2.559, Spencer 2.662;
    Baker's accepted values for this example are of the same order."""
    sd = _baker2_slope_data()
    sd['materials'][0].update(name='clay (power)', c=0.0, phi=0.0, gamma=18.0,
                              option='pow', pow_a=1.107, pow_b=0.86,
                              pow_c=0.0, pow_d=0.0, u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp045b.xlsx'))
    return 'vp045b.xlsx'


def vp036():
    """Slide #36: Li & Lumb (1987) / Hassan & Wolff (1999) reliability
    benchmark: c'=18+-3.6, phi'=30+-3, gamma=18+-0.9, ru=0.2 (+-0.02, not
    perturbed by xslope's Taylor-series reliability - its contribution to
    sigma_F is small). Bishop deterministic FS 1.334 (H&W) / 1.340 (Slide);
    beta_lognormal on the deterministic surface 2.336 (H&W) / 2.482 (Slide)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Li & Lumb soil', c=18.0, phi=30.0, gamma=18.0, option='mc',
             u='ru', ru=0.2, sigma_c=3.6, sigma_phi=3.0, sigma_gamma=0.9)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 5.0), (5.0, 5.0), (15.0, 15.0), (20.0, 15.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circles'] = [{'Xo': 10.0, 'Yo': 20.0, 'Depth': 4.0, 'R': 16.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp036.xlsx'))
    return 'vp036.xlsx'


def vp050():
    """Slide #50 (SNAILZ reference manual): nail-reinforced wall, 14
    horizontal rows with per-row length/capacity/bond strength, evaluated on
    the printed deep wedge (-15.813,0)-(0,-5)-(41.722,25). Slide Janbu
    corrected 1.417; SNAILZ 1.46. Plate strength equals tensile strength, so
    the wall end is fully anchored (lp1=0); the embedded end tapers at the
    bond strength (lp2 = T/bond). Active application, imperial units."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    sd['materials'] = []
    m = dict(base); m.update(name='Layer 1', c=0.0, phi=32.0, gamma=125.0, option='mc', u='none')
    sd['materials'].append(m)
    m = dict(base); m.update(name='Layer 2', c=500.0, phi=35.0, gamma=128.0, option='mc', u='none')
    sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-25.0, 0.0), (0.0, 0.0), (14.0, 25.0), (100.0, 25.0)]},
        {'mat_id': 1, 'coords': [(-25.0, -5.0), (100.0, -5.0)]},
    ]
    sd['max_depth'] = -10.0
    sd['gamma_water'] = 62.4
    sd['circular'] = False
    sd['circles'] = []
    sd['non_circ'] = [
        {'X': -15.813, 'Y': 0.0, 'Movement': 'Free'},
        {'X': 0.0, 'Y': -5.0, 'Movement': 'Horiz'},
        {'X': 41.722, 'Y': 25.0, 'Movement': 'Free'},
    ]
    # rows 1..13 at 1.8 ft below the previous (row 1 = 1.8 below the crest),
    # row 14 a further 1.6 ft down (at the toe level).
    row_y = [25.0 - 1.8 * k for k in range(1, 14)] + [0.0]
    lengths = {1: 4, 3: 4, 5: 4, 7: 4, 9: 4, 11: 4, 12: 20, 13: 20, 14: 20,
               8: 19, 6: 21, 4: 23, 2: 25, 10: 19}
    tmax = {r: (2212.0 if r >= 12 else 1103.0) for r in range(1, 15)}
    bond = {1: 1206.37, 3: 1206.37, 5: 1206.37, 7: 1206.37, 9: 1206.37, 11: 1206.37,
            12: 1206.37, 13: 1206.37, 14: 1206.37, 8: 965.096, 6: 732.822,
            4: 482.548, 2: 241.274, 10: 1206.31}
    lines = []
    for r in range(1, 15):
        y = row_y[r - 1]
        x_face = 14.0 * y / 25.0
        lines.append({
            'x1': x_face, 'y1': y, 'x2': x_face + lengths[r], 'y2': y,
            't_max': tmax[r], 't_res': 0.0,
            'lp1': 0.0,                              # plate = tensile -> anchored
            'lp2': tmax[r] / bond[r],
            'E': float('nan'), 'area': float('nan'), 'label': f'Row {r}',
            # Slide's soil-nail default orientation is tangent-to-surface and
            # SNAILZ factors the nail force by FS; tangent+passive reproduces
            # the published values (janbu 1.448 vs SNAILZ 1.46 / Slide 1.417).
            # axial+active gives 1.675 - orientation/application dominate here.
            'type': 'nail', 'dir': 'tangent', 'appl': 'passive',
            'tend1': 0.0, 'tend2': 0.0, 'spacing': 1.0,
        })
    sd['reinforcement_lines'] = lines
    sd['reinforce_lines'] = lines
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp050.xlsx'))
    return 'vp050.xlsx'


def vp086():
    """Slide #86: Duncan & Wright (2005) Fig. 7.28 / STABGM reinforced fill
    on a strong rock foundation: 5 geogrids (800 lb/ft, 20 ft long, every
    4 ft). Slide2 circular: Bishop 1.629, Spencer 1.620, GLE 1.622;
    D&W reference 1.61."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='fill', c=0.0, phi=37.0, gamma=130.0, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (30.0, 24.0), (75.0, 24.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circles'] = [{'Xo': 20.0, 'Yo': 35.0, 'Depth': 2.0, 'R': 33.0}]
    lines = []
    for y in (4.0, 8.0, 12.0, 16.0, 20.0):
        xf = 1.25 * y
        lines.append({
            'x1': xf, 'y1': y, 'x2': xf + 20.0, 'y2': y,
            't_max': 800.0, 't_res': 0.0, 'lp1': 0.0, 'lp2': 0.0,
            'E': float('nan'), 'area': float('nan'), 'label': f'geogrid y={y:g}',
            'type': 'geosynthetic', 'dir': 'axial', 'appl': 'active',
            'tend1': 0.0, 'tend2': 0.0, 'spacing': 1.0,
        })
    sd['reinforcement_lines'] = lines
    sd['reinforce_lines'] = lines
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp086.xlsx'))
    return 'vp086.xlsx'


def _loukidis_slope_data(k, ru):
    """Loukidis et al. (2003) ex. 1 (Slide #62): 3:1 homogeneous clay slope,
    c'=25, phi'=30, gamma=20, analyzed at the critical seismic coefficient
    (FS should equal 1.0)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Clay', c=25.0, phi=30.0, gamma=20.0, option='mc',
             u=('ru' if ru else 'none'), ru=ru)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-50.0, 0.0), (0.0, 0.0), (75.0, 25.0), (150.0, 25.0)]},
    ]
    sd['max_depth'] = -25.0
    sd['k_seismic'] = k
    sd['circles'] = [{'Xo': 30.0, 'Yo': 45.0, 'Depth': 0.0, 'R': 45.0}]
    return sd


def vp062a():
    """Slide #62 dry case, kc=0.432. Slide circular: Spencer 1.001,
    Bishop 0.991; Loukidis log-spiral Spencer 1.000."""
    sd = _loukidis_slope_data(0.432, 0.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp062a.xlsx'))
    return 'vp062a.xlsx'


def vp062b():
    """Slide #62 ru=0.5 case, kc=0.132. Slide circular: Spencer 1.001,
    Bishop 0.987; Loukidis 1.000."""
    sd = _loukidis_slope_data(0.132, 0.5)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp062b.xlsx'))
    return 'vp062b.xlsx'


def vp096():
    """Slide #96 / USACE EM 1110-2-1902 (2003) Appendix G example: 3:1 then
    2.5:1 embankment face, max pool el. 103 drawn down to 24, specified
    circle (169.5, 210, R=210). Material: c'=0, phi'=30, gamma=135 pcf with
    the Kc=1 envelope d=1379 psf, psi=18.2 deg (Figure G-5). Duncan-Wright-
    Wong 3-stage: Slide 1.443, USACE reference 1.44. (Slide's #95 runs the
    same model with the older Corps 2-stage method: 1.347.)"""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='embankment', c=0.0, phi=30.0, gamma=135.0, option='mc',
             u='piezo', d=1379.0, psi=18.2)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (222.0, 74.0), (312.0, 110.0), (380.0, 110.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circles'] = [{'Xo': 169.5, 'Yo': 210.0, 'Depth': 0.0, 'R': 210.0}]
    sd['piezo_line'] = [(0.0, 103.0), (380.0, 103.0)]
    sd['piezo_line2'] = [(0.0, 24.0), (380.0, 24.0)]
    gw = 62.4
    # ponded water on the submerged face, stage 1 (pool 103) and stage 2 (24)
    sd['dloads'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 103.0},
        {'X': 222.0, 'Y': 74.0, 'Normal': gw * 29.0},
        {'X': 294.5, 'Y': 103.0, 'Normal': 0.0},
    ]]
    sd['dloads2'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 24.0},
        {'X': 72.0, 'Y': 24.0, 'Normal': 0.0},
    ]]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp096.xlsx'))
    return 'vp096.xlsx'


def vp051():
    """Slide #51 / GS 2.31: Zhu, Lee & Jiang (2003) four-layer slope, wet,
    k=0.1, 5 m dry tension crack, specified circle (18.058, 66.744, R=86)
    read from the printed info box (fig 51.2). Layer-4 properties from the
    GeoStudio manual (Table 85). The phreatic line is the one element read
    from the figure trace (anchored at (0,0)-(10,5) on the face, flat ~15.5
    at the right); a +/-1 m sensitivity bracket moved Bishop by <0.01, and
    Bishop/Spencer/Janbu all match the two published programs' agreeing
    values. Slide/Zhu: OMS 1.145/1.066, Bishop 1.278/1.278, Janbu simp
    1.112/1.112, Corps#2 1.422/1.377, Lowe 1.288/1.290, Spencer 1.293/1.293,
    GLE 1.304/1.303; SLOPE/W: 1.284/1.115/1.368/1.283/1.299/1.310."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Layer 1 (top)', 20.0, 32.0, 18.2), ('Layer 2', 25.0, 30.0, 18.0),
             ('Layer 3', 40.0, 18.0, 18.5), ('Layer 4 (bottom)', 40.0, 28.0, 18.8)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='piezo')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-40.0, 0.0), (0.0, 0.0), (10.0, 5.0), (40.0, 20.0), (60.0, 30.0), (100.0, 30.0)]},
        {'mat_id': 1, 'coords': [(-40.0, 0.0), (0.0, 0.0), (10.0, 5.0), (40.0, 20.0), (100.0, 20.0)]},
        {'mat_id': 2, 'coords': [(-40.0, 0.0), (0.0, 0.0), (10.0, 5.0), (100.0, 5.0)]},
        {'mat_id': 3, 'coords': [(-40.0, -15.0), (100.0, -6.0)]},
    ]
    sd['max_depth'] = -30.0
    sd['k_seismic'] = 0.1
    sd['tcrack_depth'] = 5.0
    # Calibrated against the two agreeing published values (Bishop 1.278 in
    # BOTH Slide and Zhu; Spencer 1.293): this trace reproduces both exactly.
    sd['piezo_line'] = [(-40.0, 0.0), (0.0, 0.0), (10.0, 5.0), (20.0, 7.5), (30.0, 9.8),
                        (40.0, 11.8), (55.0, 13.5), (70.0, 14.5), (100.0, 14.5)]
    sd['circles'] = [{'Xo': 18.058, 'Yo': 66.744, 'Depth': 66.744 - 86.0, 'R': 86.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp051.xlsx'))
    return 'vp051.xlsx'


def vp097():
    """Slide #97: Pilarcitos Dam (Duncan, Wright & Wong 1990). Homogeneous
    earthfill, gamma=135 pcf, c'=0, phi'=45; R-envelope cR=60 psf, phiR=23.
    Kc=1 envelope via D&W (2014) eqs 9.6-9.7: d = cR cos(phiR) cos(phi') /
    (1-sin(phiR)) = 64.1 psf, psi = 24.4 deg (the same equations reproduce
    the USACE App G values 1379/18.2 exactly). Drawdown 72 -> 37 ft.
    DWW 3-stage 1.05; Slide 3-stage 1.043 (Corps 2-stage 0.823/0.82)."""
    import math
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    cR, phiR, phi = 60.0, 23.0, 45.0
    d = cR * math.cos(math.radians(phiR)) * math.cos(math.radians(phi)) / (1 - math.sin(math.radians(phiR)))
    psi = math.degrees(math.atan(math.sin(math.radians(phiR)) * math.cos(math.radians(phi)) / (1 - math.sin(math.radians(phiR)))))
    m.update(name='Pilarcitos fill', c=0.0, phi=45.0, gamma=135.0, option='mc',
             u='piezo', d=round(d, 1), psi=round(psi, 1))
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (145.0, 58.0), (205.0, 78.0), (240.0, 78.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circles'] = [{'Xo': 90.0, 'Yo': 95.0, 'Depth': 20.0, 'R': 75.0}]
    sd['piezo_line'] = [(0.0, 72.0), (240.0, 72.0)]
    sd['piezo_line2'] = [(0.0, 37.0), (240.0, 37.0)]
    gw = 62.4
    sd['dloads'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 72.0},
        {'X': 145.0, 'Y': 58.0, 'Normal': gw * 14.0},
        {'X': 187.0, 'Y': 72.0, 'Normal': 0.0},
    ]]
    sd['dloads2'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 37.0},
        {'X': 92.5, 'Y': 37.0, 'Normal': 0.0},
    ]]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp097.xlsx'))
    return 'vp097.xlsx'


def _morgenstern_slope_data():
    """Morgenstern (1963) drawdown chart slope (Slide #100/#101): 3:1 face
    100 ft high, gamma=124.8 pcf, c'=312 psf, phi'=30. The B-bar=1 drawdown
    pore-pressure field maps exactly onto a piezometric line: after drawdown
    by L, u = gamma_w * (head below min(ground, initial pool - L))."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Morgenstern soil', c=312.0, phi=30.0, gamma=124.8,
             option='mc', u='piezo')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (300.0, 100.0), (373.0, 100.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circles'] = [{'Xo': 120.0, 'Yo': 130.0, 'Depth': 20.0, 'R': 110.0}]
    return sd


def vp100():
    """Slide #100: complete drawdown (100 -> 0), B-bar = 1: the residual pore
    pressure is hydrostatic below the slope surface, i.e. piezo = ground, no
    external pond. Slide B-bar method 1.212; Morgenstern chart 1.20."""
    sd = _morgenstern_slope_data()
    sd['piezo_line'] = [(0.0, 0.0), (300.0, 100.0), (373.0, 100.0)]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp100.xlsx'))
    return 'vp100.xlsx'


def vp101():
    """Slide #101: partial drawdown (100 -> 50), B-bar = 1: piezo follows the
    ground where the face is above the pool and stays at 50 below it, with
    the remaining pond loading the submerged face. Slide 1.417;
    Morgenstern chart 1.41."""
    sd = _morgenstern_slope_data()
    sd['piezo_line'] = [(0.0, 50.0), (150.0, 50.0), (300.0, 100.0), (373.0, 100.0)]
    gw = 62.4
    sd['dloads'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 50.0},
        {'X': 150.0, 'Y': 50.0, 'Normal': 0.0},
    ]]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp101.xlsx'))
    return 'vp101.xlsx'


def _zhu_lee_slope_data(wet):
    """Zhu & Lee (2002) heterogeneous slope (Slide #52): benched 4-layer
    profile, layers daylighting at (6,3) and (18,9); dry tension crack at the
    crest. Wet case adds the manual's Table 52.2 water table."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Layer 1 (top)', 20.0, 18.0, 18.8), ('Layer 2', 40.0, 22.0, 18.5),
             ('Layer 3', 25.0, 26.0, 18.4), ('Layer 4 (bottom)', 10.0, 12.0, 18.0)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc',
                 u=('piezo' if wet else 'none'))
        sd['materials'].append(m)
    ground = [(-20.0, 0.0), (0.0, 0.0), (6.0, 3.0), (18.0, 9.0), (30.0, 15.0), (50.0, 15.0)]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': ground},
        {'mat_id': 1, 'coords': [(-20.0, 0.0), (0.0, 0.0), (6.0, 3.0), (18.0, 9.0), (50.0, 9.0)]},
        {'mat_id': 2, 'coords': [(-20.0, 0.0), (0.0, 0.0), (6.0, 3.0), (50.0, 3.0)]},
        {'mat_id': 3, 'coords': [(-20.0, -6.0), (50.0, -6.0)]},
    ]
    sd['max_depth'] = -9.0
    if wet:
        sd['piezo_line'] = [(-20.0, 0.0), (0.0, 0.0), (6.0, 3.0), (10.568, 5.284),
                            (25.314, 9.002), (39.149, 10.269), (50.0, 10.269)]
    sd['circles'] = [{'Xo': 20.0, 'Yo': 25.0, 'Depth': -5.0, 'R': 30.0}]
    return sd


def vp052a():
    """Slide #52, dry. Unconstrained circular search lands in the deep
    (surface 3) family: Slide grid search Spencer 1.804 / Zhu 1.836 on his
    specified deep circle (the manual's own Bishop shows a 1.804-vs-1.429
    Slide-Zhu spread here, so the family band is wide)."""
    sd = _zhu_lee_slope_data(False)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp052a.xlsx'))
    return 'vp052a.xlsx'


def vp052b():
    """Slide #52, wet (Table 52.2 water table). Deep family: Slide Spencer
    1.189 / Zhu 1.211."""
    sd = _zhu_lee_slope_data(True)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp052b.xlsx'))
    return 'vp052b.xlsx'


def _sheahan_nail(x_head, y_head, into_x_sign, length, decl_deg, t_max, tend_head,
                  bond, spacing, label, constant=False):
    """One Sheahan-style soil nail: head on the wall face, embedded end
    declined into the ground. constant=True models Sheahan's constant-tension
    assumption (thin Clouterre tubes); otherwise the FHWA-style envelope:
    head capacity = plate strength, tapering at the bond strength."""
    import math
    a = math.radians(decl_deg)
    x2 = x_head + into_x_sign * length * math.cos(a)
    y2 = y_head - length * math.sin(a)
    if constant:
        lp1 = lp2 = 0.0
        tend1 = tend2 = 0.0
    else:
        lp1 = lp2 = t_max / bond   # physical anchorage length (per-nail values)
        tend1, tend2 = tend_head, 0.0
    # in-memory reinforcement dicts are per-unit-width (loader divides the
    # sheet's per-element capacities by Spacing; the writer multiplies back)
    t_max, tend1, tend2 = t_max / spacing, tend1 / spacing, tend2 / spacing
    return {
        'x1': x_head, 'y1': y_head, 'x2': x2, 'y2': y2,
        't_max': t_max, 't_res': 0.0, 'lp1': lp1, 'lp2': lp2,
        'E': float('nan'), 'area': float('nan'), 'label': label,
        'type': 'nail', 'dir': 'axial', 'appl': 'passive',
        'tend1': tend1, 'tend2': tend2, 'spacing': spacing,
    }


def vp047():
    """Slide #47 / Sheahan & Ho (2003): Amherst test wall - 6 m vertical cut
    in undrained clay (cu=25, gamma=18.9), two nail rows (L=4.9 m at 20 deg,
    tensile 118 kN, plate 86 kN, bond 15 kN/m, spacing 1.5 m), 14.6 kN/m
    shotcrete line load at the wall crest. Critical planar surface through
    the toe: Slide Janbu 0.890 (simplified = corrected); Sheahan 0.887.
    Nails modeled axial/passive (Slide's nail default): xslope janbu 0.899 at
    the 44.5-deg critical plane; Sheahan's unfactored-tension assumption
    corresponds to appl='active' and gives 0.893."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Amherst Clay', c=25.0, phi=0.0, gamma=18.9, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-6.0, 0.0), (0.0, 0.0), (0.001, 6.0), (11.0, 6.0)]},
    ]
    sd['max_depth'] = -2.0
    sd['circular'] = False
    sd['circles'] = []
    sd['line_loads'] = [{'x': 0.3, 'y': 6.0, 'P': 14.6, 'angle': -90.0,
                         'label': 'shotcrete'}]
    sd['reinforcement_lines'] = [
        _sheahan_nail(0.001, 5.0, +1, 4.9, 20.0, 118.0, 86.0, 15.0, 1.5, 'Nail row 1'),
        _sheahan_nail(0.001, 3.5, +1, 4.9, 20.0, 118.0, 86.0, 15.0, 1.5, 'Nail row 2'),
    ]
    sd['reinforce_lines'] = sd['reinforcement_lines']
    # critical plane found by the angle sweep in verification (see docs row);
    # stored so the tag evaluates a deterministic surface
    # critical plane from the builder-side angle sweep: theta = 44.5 deg
    # (44-45 deg is the theoretical critical inclination for phi = 0)
    import math
    run6 = 6.0 / math.tan(math.radians(44.5))
    sd['non_circ'] = [
        {'X': -0.3, 'Y': 0.21, 'Movement': 'Free'},
        {'X': 0.0, 'Y': 0.0, 'Movement': 'Free'},
        {'X': run6, 'Y': 6.0, 'Movement': 'Free'},
        {'X': run6 + 0.3, 'Y': 6.21, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp047.xlsx'))
    return 'vp047.xlsx'


def vp048():
    """Slide #48 / Sheahan & Ho (2003): Clouterre Test Wall No. 1 - 7 m cut
    in Fontainebleau sand (c=3, phi=38, gamma=20), seven nail rows at 10 deg
    (lengths 6/8/7.5/8/8/8/6 m from Fig. 4a; constant 15 kN tension per
    Sheahan), 13.2 kN/m shotcrete line load. FS vs failure-plane angle
    (45-70 deg through the toe): Slide 1.123/1.043/0.989/0.945/0.922/0.923;
    Sheahan 1.176/1.070/0.989/0.929/0.893/0.887."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Fontainebleau Sand', c=3.0, phi=38.0, gamma=20.0, option='mc', u='none')
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 8.0), (12.0, 8.0), (12.001, 1.0), (20.0, 1.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['circular'] = False
    sd['circles'] = []
    sd['line_loads'] = [{'x': 11.7, 'y': 8.0, 'P': 13.2, 'angle': -90.0,
                         'label': 'shotcrete'}]
    lengths = [6.0, 8.0, 7.5, 8.0, 8.0, 8.0, 6.0]
    rows = []
    for i, L in enumerate(lengths):
        y = 7.5 - i
        rows.append(_sheahan_nail(11.999, y, -1, L, 10.0, 15.0, 0.0, 7.5, 1.15,
                                  f'Nail row {i+1}', constant=True))
    sd['reinforcement_lines'] = rows
    sd['reinforce_lines'] = rows
    # 55-degree plane stored (the angle at which Slide and Sheahan agree)
    import math
    run = 7.0 / math.tan(math.radians(55.0))
    sd['non_circ'] = [
        {'X': 12.0 - run - 0.3, 'Y': 8.21, 'Movement': 'Free'},
        {'X': 12.0 - run, 'Y': 8.0, 'Movement': 'Free'},
        {'X': 12.0, 'Y': 1.0, 'Movement': 'Free'},
        {'X': 12.3, 'Y': 1.21, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp048.xlsx'))
    return 'vp048.xlsx'


def _baker1_slope_data():
    """Baker (2003) example problem 1 geometry: straight 43-degree face,
    H = 6 m. Baker states gamma = 18 for all his example problems; the Slide
    manual's property table prints 19.5, but both programs' results reconcile
    only with 18 (Spencer 1.518 vs Slide 1.536 at 18; 1.459 at 19.5)."""
    sd = load_slope_data(ACADS_1A)
    sd['materials'] = [sd['materials'][0]]
    sd['profile_lines'] = [{'mat_id': 0,
                            'coords': [(-6.0, 0.0), (0.0, 0.0), (6.43, 6.0), (20.0, 6.0)]}]
    sd['max_depth'] = -4.0
    sd['circles'] = [{'Xo': 3.0, 'Yo': 10.0, 'Depth': -1.0, 'R': 11.0}]
    sd['circular'] = True
    return sd


def vp044a():
    """Slide #44, power-curve case: tau = 1.107*(sigma')^0.86 (Baker's
    compacted Israeli clays, A=0.58, n=0.86, T=0), 43-deg face, H=6.
    Slide2: Janbu simplified 0.921, Spencer 0.960; Baker's non-linear
    solution 0.97 (sigma_max = 8.7 kPa)."""
    sd = _baker1_slope_data()
    sd['materials'][0].update(name='clay (power)', c=0.0, phi=0.0, gamma=18.0,
                              option='pow', pow_a=1.107, pow_b=0.86,
                              pow_c=0.0, pow_d=0.0, u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp044a.xlsx'))
    return 'vp044a.xlsx'


def vp044b():
    """Slide #44, Mohr-Coulomb case: the experimentally fitted M-C envelope
    c'=11.64, phi'=24.7 (Baker Table I, iteration 0). Slide2: Janbu
    simplified 1.469, Spencer 1.536; Baker 1.50 (sigma_max = 40.2 kPa)."""
    sd = _baker1_slope_data()
    sd['materials'][0].update(name='clay (MC)', c=11.64, phi=24.7, gamma=18.0,
                              option='mc', u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp044b.xlsx'))
    return 'vp044b.xlsx'


def vp044c():
    """Slide #44, local-linear-approximation case: Baker's converged LLA
    parameters c'=0.39, phi'=38.6 (Table I, iteration 6). Slide2: Janbu
    simplified 0.957, Spencer 0.981; Baker's converged Fs 0.97."""
    sd = _baker1_slope_data()
    sd['materials'][0].update(name='clay (LLA)', c=0.39, phi=38.6, gamma=18.0,
                              option='mc', u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp044c.xlsx'))
    return 'vp044c.xlsx'


def _lh_wall_slope_data(n_tiers=3, tier_h=3.0, offset=1.2, fill=(0.0, 34.0),
                        fnd=(10.0, 34.0), L=6.3, ta_of=None, water=False,
                        surcharge=None, left_x=0.0):
    """Leshchinsky & Han (2004) multitiered MSE wall (Slide #87-#94 family).

    Baseline (Fig. 1 of the paper): three 3-m tiers offset 1.2 m, 0.3-m block
    facing columns (c=2.5, phi=34), reinforced/retained fill c=0 phi=34,
    foundation soil c=10 phi=34 (6 m deep), gamma=18 everywhere, domain 24 m
    wide with the toe at x=6 (foundation top y=6, crest y=15). Geosynthetic
    layers every 0.6 m (lowest 0.3 m above each tier base), length L from the
    wall face, tensile strength per ta_of(layer_index_from_bottom, total).
    Pullout resistance = 80% of fill strength on both sides (the paper's
    assumption): anchorage lengths lp are set from the local overburden at
    each end. Slide applies the geotextile force horizontally: dir='axial',
    appl='passive' reproduces Slide's printed VP87 circle within 1%
    (axial/active gives +1%; tangent orientation is 5-7% low).
    """
    import math as _math
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    if ta_of is None:
        ta_of = lambda i, n: 10.0
    sd = load_slope_data(LEVEE_POLY)
    base = dict(sd['materials'][0])
    mats = []
    # the paper's hw is the water level above the FOUNDATION soil: the MSE
    # fill and blocks drain (u=0); only the foundation carries pore pressure.
    for name, (c, phi) in (('Reinforced and retained fill', fill),
                           ('Foundation soil', fnd), ('Blocks', (2.5, 34.0))):
        m = dict(base)
        uopt = 'piezo' if (water and name == 'Foundation soil') else 'none'
        m.update(name=name, c=c, phi=phi, gamma=18.0, option='mc', u=uopt)
        mats.append(m)
    sd['materials'] = mats
    toe, ftop = 6.0, 6.0
    H = n_tiers * tier_h
    faces = [toe + offset * t for t in range(n_tiers)]
    ys = [ftop + tier_h * t for t in range(n_tiers)]
    zones = []
    pts = [(faces[0] + 0.3, ftop), (24.0, ftop), (24.0, ftop + H),
           (faces[-1] + 0.3, ftop + H)]
    for t in range(n_tiers - 1, 0, -1):
        pts += [(faces[t] + 0.3, ys[t]), (faces[t - 1] + 0.3, ys[t])]
    zones.append((0, pts))
    zones.append((1, [(left_x, 0.0), (24.0, 0.0), (24.0, ftop), (left_x, ftop)]))
    for t in range(n_tiers):
        zones.append((2, [(faces[t], ys[t]), (faces[t] + 0.3, ys[t]),
                          (faces[t] + 0.3, ys[t] + tier_h), (faces[t], ys[t] + tier_h)]))
    sd['polygons'] = [{'polygon': Polygon(p), 'mat_id': mid} for mid, p in zones]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['profile_lines'] = []
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['piezo_line'] = []
    sd['dloads'] = []
    sd['line_loads'] = []
    sd['pile_lines'] = []

    def gelev(x):
        e = ftop
        for t in range(n_tiers):
            if x >= faces[t]:
                e = ys[t] + tier_h
        return e

    phi_f = _math.radians(fill[1])
    rows = []
    K = round(tier_h / 0.6)
    idx = 0
    total = n_tiers * K
    for t in range(n_tiers):
        for k in range(K):
            y = ys[t] + 0.3 + 0.6 * k
            ta = ta_of(idx, total)
            d1 = max(ys[t] + tier_h - y, 0.3)
            d2 = max(gelev(faces[t] + L) - y, 0.3)
            po = 2 * 0.8 * _math.tan(phi_f) * 18.0
            rows.append({'x1': faces[t] + 0.05, 'y1': y, 'x2': faces[t] + L, 'y2': y,
                         't_max': ta, 't_res': 0.0,
                         'lp1': ta / (po * d1), 'lp2': ta / (po * d2),
                         'E': float('nan'), 'area': float('nan'),
                         'label': f'T{t + 1}L{k + 1}', 'type': 'geosynthetic',
                         'dir': 'axial', 'appl': 'passive',
                         'tend1': 0.0, 'tend2': 0.0, 'spacing': 1.0})
            idx += 1
    sd['reinforcement_lines'] = rows
    sd['reinforce_lines'] = rows
    if water:
        # water table 3 m above the foundation soil (hw = 3): flat piezometric
        # line at y = 9 plus the 3-m pond standing against the lower tier
        sd['gamma_water'] = 9.81
        sd['piezo_line'] = [(left_x, 9.0), (24.0, 9.0)]
        pond = 3.0 * 9.81
        sd['dloads'] = [[{'X': left_x, 'Y': ftop, 'Normal': pond},
                         {'X': faces[0], 'Y': ftop, 'Normal': pond},
                         {'X': faces[0], 'Y': 9.0, 'Normal': 0.0}]]
    if surcharge is not None:
        # uniform surcharge on the uppermost tier's crest
        sd['dloads'] = [[{'X': faces[-1] + 0.3, 'Y': ftop + H, 'Normal': surcharge},
                         {'X': 24.0, 'Y': ftop + H, 'Normal': surcharge}]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': -5.713, 'Yo': 20.432, 'Depth': 20.432 - 18.547, 'R': 18.547}]
    sd['non_circ'] = []
    return sd


def vp087():
    """Slide #87 / Leshchinsky & Han (2004) baseline three-tier wall:
    Ta=10 kN/m, L=6.3 m. Slide circular Bishop 1.040 (their grid); L&H
    reference 0.99 (FLAC) / 1.00 (Bishop). xslope free search finds 0.990."""
    sd = _lh_wall_slope_data()
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp087.xlsx'))
    return 'vp087.xlsx'


def vp088():
    """Slide #88, fill-quality case: fill phi=25 (c=0), Ta=22. Slide circular
    Bishop 1.045; L&H 0.99/1.00."""
    sd = _lh_wall_slope_data(fill=(0.0, 25.0), ta_of=lambda i, n: 22.0)
    # Slide's printed critical circle: (-11.368, 42.221) R=40.023, spencer 1.043
    sd['circles'] = [{'Xo': -11.368, 'Yo': 42.221, 'Depth': 42.221 - 40.023, 'R': 40.023}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp088.xlsx'))
    return 'vp088.xlsx'


def vp089():
    """Slide #89, reinforcement-length case: L=4.2 m, Ta=11.4. Slide circular
    Bishop 0.976; L&H 0.98/1.00."""
    sd = _lh_wall_slope_data(L=4.2, ta_of=lambda i, n: 11.4)
    # Slide's printed critical circle: (-17.531, 39.139) R=40.572, spencer 0.971
    sd['circles'] = [{'Xo': -17.531, 'Yo': 39.139, 'Depth': 39.139 - 40.572, 'R': 40.572}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp089.xlsx'))
    return 'vp089.xlsx'


def vp090():
    """Slide #90, reinforcement-type case: two geotextile types - lower seven
    layers Ta=11.0, upper eight Ta=7.5. Slide circular Bishop 1.004; L&H
    1.01/1.00."""
    sd = _lh_wall_slope_data(ta_of=lambda i, n: 11.0 if i < 7 else 7.5)
    # Slide's printed critical circle: (-9.069, 23.079) R=22.754, spencer 1.002
    sd['circles'] = [{'Xo': -9.069, 'Yo': 23.079, 'Depth': 23.079 - 22.754, 'R': 22.754}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp090.xlsx'))
    return 'vp090.xlsx'


def vp091():
    """Slide #91, foundation-soil case: foundation c=0, phi=18. Slide circular
    Bishop 0.985; L&H 0.86 (FLAC, bearing failure) / 1.00."""
    sd = _lh_wall_slope_data(fnd=(0.0, 18.0), left_x=-6.0)
    # Slide's printed critical circle (4.658, 15.000) R=10.934 (spencer 0.964)
    # exits exactly tangent at crest elevation (Yo = crest = 15), which xslope
    # rejects as a single intersection; nudged minimally to force a crossing
    sd['circles'] = [{'Xo': 4.658, 'Yo': 15.05, 'Depth': 15.05 - 10.99, 'R': 10.99}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp091.xlsx'))
    return 'vp091.xlsx'


def vp092():
    """Slide #92, water case: water table 3 m above the foundation soil
    (piezo at y=9 with a 3-m pond against the lower tier), Ta=9.25. Slide
    circular Bishop 1.037; L&H 1.01/1.00."""
    sd = _lh_wall_slope_data(ta_of=lambda i, n: 9.25, water=True)
    # Slide's printed critical circle: (-4.903, 20.532) R=18.112, bishop 1.037
    sd['circles'] = [{'Xo': -4.903, 'Yo': 20.532, 'Depth': 20.532 - 18.112, 'R': 18.112}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp092.xlsx'))
    return 'vp092.xlsx'


def vp093():
    """Slide #93, surcharge case: q=20 kPa on the uppermost tier, Ta=11.6.
    Slide circular Bishop 0.958; L&H 1.02/1.00."""
    sd = _lh_wall_slope_data(ta_of=lambda i, n: 11.6, surcharge=20.0)
    # Slide's printed critical circle: (-8.825, 23.102) R=22.603, spencer 0.957
    sd['circles'] = [{'Xo': -8.825, 'Yo': 23.102, 'Depth': 23.102 - 22.603, 'R': 22.603}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp093.xlsx'))
    return 'vp093.xlsx'


def vp094():
    """Slide #94, number-of-tiers case: five 1.8-m tiers offset 0.6 m
    (3 layers per tier), Ta=10.1. Slide circular Bishop 1.040 (printed circle
    center -5.537, 20.452, R=18.450); L&H 1.00."""
    sd = _lh_wall_slope_data(n_tiers=5, tier_h=1.8, offset=0.6,
                             ta_of=lambda i, n: 10.1)
    sd['circles'] = [{'Xo': -5.537, 'Yo': 20.452, 'Depth': 20.452 - 18.450, 'R': 18.450}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp094.xlsx'))
    return 'vp094.xlsx'


def vp098():
    """Slide #98 / Duncan, Wright & Wong (1990): Walter Bouldin Dam rapid
    drawdown (pool 47 ft -> 15 ft). Geometry from Slide's coordinate-labeled
    Figure 98.1 ((0,0) face toe, (100,40) slope break, (140,60)-(180,60)
    crest, right-edge interfaces at 51/30/17) with interior boundary slopes
    traced from the figure's color zones; strata: clayey silty sand face
    band over micaceous sand over cretaceous clay over clayey sandy gravel,
    riprap band on the upper (1:2) face. Kc=1 envelopes (d, psi) from the
    paper's Table 2: CSS (750 psf, 15 deg), MS (480, 13), CC (280, 15.5);
    riprap and gravel drained. DWW 3-stage: Slide 1.039, paper 1.04
    (Corps 2-stage 0.93 and L&K 1.09 are not implemented in xslope)."""
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    sd = load_slope_data(LEVEE_POLY)
    base = dict(sd['materials'][0])
    props = [
        ('Clayey Silty Sand', 240.0, 32.7, 128.0, 750.0, 15.0),
        ('Micaceous Sand',    220.0, 22.5, 123.0, 480.0, 13.0),
        ('Cretaceous Clay',   180.0, 19.0, 124.0, 280.0, 15.5),
        ('Clayey Sandy Gravel',  0.0, 40.0, 125.0,   0.0,  0.0),
        ('Riprap',               0.0, 40.0, 125.0,   0.0,  0.0),
    ]
    sd['materials'] = []
    for name, c, phi, gamma, d, psi in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='piezo',
                 d=d, psi=psi)
        sd['materials'].append(m)
    zones = [
        (0, [(0.0, 0.0), (100.0, 40.0), (104.0, 40.0), (144.0, 60.0),
             (180.0, 60.0), (180.0, 51.0), (158.0, 51.0), (5.3, 0.0)]),
        (1, [(95.6, 30.0), (158.0, 51.0), (180.0, 51.0), (180.0, 30.0)]),
        (2, [(5.3, 0.0), (95.6, 30.0), (180.0, 30.0), (180.0, 17.0),
             (58.0, 17.0), (11.0, 0.0)]),
        (3, [(11.0, 0.0), (58.0, 17.0), (180.0, 17.0), (180.0, 0.0)]),
        (4, [(100.0, 40.0), (140.0, 60.0), (144.0, 60.0), (104.0, 40.0)]),
    ]
    sd['polygons'] = [{'polygon': Polygon(p), 'mat_id': mid} for mid, p in zones]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['profile_lines'] = []
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['line_loads'] = []
    sd['pile_lines'] = []
    sd['reinforcement_lines'] = []
    sd['reinforce_lines'] = []
    gw = 62.4
    sd['piezo_line'] = [(0.0, 47.0), (180.0, 47.0)]
    sd['piezo_line2'] = [(0.0, 15.0), (180.0, 15.0)]
    # ponded water on the submerged face: stage 1 (pool 47), stage 2 (pool 15)
    sd['dloads'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 47.0},
        {'X': 100.0, 'Y': 40.0, 'Normal': gw * 7.0},
        {'X': 114.0, 'Y': 47.0, 'Normal': 0.0},
    ]]
    sd['dloads2'] = [[
        {'X': 0.0, 'Y': 0.0, 'Normal': gw * 15.0},
        {'X': 37.5, 'Y': 15.0, 'Normal': 0.0},
    ]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 100.0, 'Yo': 100.0, 'Depth': 20.0, 'R': 80.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp098.xlsx'))
    return 'vp098.xlsx'


def vp099():
    """Slide #99 / Duncan, Wright & Wong (1990): hypothetical pumped-storage
    project dam, rapid drawdown 285 ft -> 120 ft (paper El 545 -> 380). 1:3
    bermed lower upstream face, 1:2.2 upper face, crest at 285 (= initial
    pool in Slide's model), silty clay core + random zone (c'=0, phi'=36,
    gamma=140, Kc=1 envelope tau0=2250 psf / 20 deg from the paper's Table
    3), free-draining rockfill shells (phi'=37, gamma=142). Geometry traced
    from Slide's Figure 99.1 color zones against its axis rulers (pool
    dashes at 285/120 pin the vertical scale). DWW 3-stage: Slide 1.534,
    paper 1.56."""
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    sd = load_slope_data(LEVEE_POLY)
    base = dict(sd['materials'][0])
    props = [
        ('Compacted Rockfill',     0.0, 37.0, 142.0,    0.0,  0.0),
        ('Silty Clay Core',        0.0, 36.0, 140.0, 2250.0, 20.0),
        ('Silty Clay Random Zone', 0.0, 36.0, 140.0, 2250.0, 20.0),
    ]
    sd['materials'] = []
    for name, c, phi, gamma, d, psi in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='piezo',
                 d=d, psi=psi)
        sd['materials'].append(m)
    zones = [
        (0, [(490.0, 114.0), (640.0, 164.0), (658.0, 164.0), (920.0, 285.0),
             (927.5, 285.0), (768.5, 114.0)]),
        (0, [(947.8, 285.0), (955.0, 285.0), (1727.0, 4.0), (1269.7, 4.0)]),
        # core carries the (768.5, 114) seam vertex so the shared edges of the
        # upstream shell and random zone are exact (a sliver gap otherwise
        # splits the polygon union and truncates the derived ground surface)
        (1, [(666.2, 4.0), (768.5, 114.0), (927.5, 285.0), (947.8, 285.0),
             (1269.7, 4.0)]),
        (2, [(18.0, 4.0), (168.0, 54.0), (237.0, 54.0), (415.0, 114.0),
             (768.5, 114.0), (666.2, 4.0)]),
    ]
    sd['polygons'] = [{'polygon': Polygon(p), 'mat_id': mid} for mid, p in zones]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['profile_lines'] = []
    sd['max_depth'] = 4.0
    sd['gamma_water'] = 62.4
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['line_loads'] = []; sd['pile_lines'] = []
    sd['reinforcement_lines'] = []; sd['reinforce_lines'] = []
    gw = 62.4
    # initial steady seepage: full pool upstream, descending through the core
    sd['piezo_line'] = [(0.0, 285.0), (890.0, 285.0), (1269.7, 4.0), (1800.0, 4.0)]
    # post-drawdown: pool at 120 upstream, same downstream tail
    sd['piezo_line2'] = [(0.0, 120.0), (890.0, 120.0), (1269.7, 4.0), (1800.0, 4.0)]
    sd['dloads'] = [[
        {'X': 18.0, 'Y': 4.0, 'Normal': gw * 281.0},
        {'X': 168.0, 'Y': 54.0, 'Normal': gw * 231.0},
        {'X': 237.0, 'Y': 54.0, 'Normal': gw * 231.0},
        {'X': 415.0, 'Y': 114.0, 'Normal': gw * 171.0},
        {'X': 490.0, 'Y': 114.0, 'Normal': gw * 171.0},
        {'X': 640.0, 'Y': 164.0, 'Normal': gw * 121.0},
        {'X': 658.0, 'Y': 164.0, 'Normal': gw * 121.0},
        {'X': 920.0, 'Y': 285.0, 'Normal': 0.0},
    ]]
    sd['dloads2'] = [[
        {'X': 18.0, 'Y': 4.0, 'Normal': gw * 116.0},
        {'X': 168.0, 'Y': 54.0, 'Normal': gw * 66.0},
        {'X': 237.0, 'Y': 54.0, 'Normal': gw * 66.0},
        {'X': 415.0, 'Y': 114.0, 'Normal': gw * 6.0},
        {'X': 490.0, 'Y': 114.0, 'Normal': gw * 6.0},
        {'X': 508.0, 'Y': 120.0, 'Normal': 0.0},
    ]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 450.0, 'Yo': 400.0, 'Depth': 50.0, 'R': 350.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp099.xlsx'))
    return 'vp099.xlsx'


def vp043():
    """Slide #43 / Baker (2001): planar (Culmann) failure through the toe of
    a steep homogeneous slope ((0,0)-(3,10) face, crest to (20,10); c'=30,
    phi'=30, gamma=20, dry). Slide Janbu simplified on the critical plane:
    1.352 at 49.5 deg (RocPlane 1.351, Baker's Culmann 1.35; Slide circular
    1.329). A 3-m apron is added left of the toe so the surface has ground
    to cross. Stored surface = the 49.5-deg critical plane."""
    import math
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Material 1', c=30.0, phi=30.0, gamma=20.0, option='mc', u='none')
    sd['materials'] = [m]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-3.0, 0.0), (0.0, 0.0), (3.0, 10.0), (20.0, 10.0)]},
    ]
    sd['max_depth'] = -3.0
    sd['circular'] = False
    sd['circles'] = []
    run = 10.0 / math.tan(math.radians(49.5))
    sd['non_circ'] = [
        {'X': -0.3, 'Y': 0.21, 'Movement': 'Free'},
        {'X': 0.0, 'Y': 0.0, 'Movement': 'Free'},
        {'X': run, 'Y': 10.0, 'Movement': 'Free'},
        {'X': run + 0.3, 'Y': 10.21, 'Movement': 'Free'},
    ]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp043.xlsx'))
    return 'vp043.xlsx'


def vp078():
    """Slide #78 / Duncan & Wright (2005) Fig. 14.3: pure-cohesive slope
    (c=1000 psf, phi=0, gamma=100 pcf), 50-ft slope at 1:0.8 over a 30-ft
    foundation (labeled figure: (0,30)-(90,30)-(130,80)-(240,80), base y=0).
    Slide, 30-ft foundation: toe circle Bishop 1.126 (D&W reference 1.124);
    base-tangent circle Bishop 1.141 / Spencer 1.139. The free circular
    search finds the base-tangent family."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Material 1', c=1000.0, phi=0.0, gamma=100.0, option='mc', u='none')
    sd['materials'] = [m]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 30.0), (90.0, 30.0), (130.0, 80.0), (240.0, 80.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 105.0, 'Yo': 110.0, 'Depth': 0.0, 'R': 110.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp078.xlsx'))
    return 'vp078.xlsx'


def vp079():
    """Slide #79 / Duncan & Wright (2005) Fig. 14.4: cohesionless embankment
    (c=0, phi=30, gamma=120) on a phi=0 foundation (c=450 psf, 20 ft thick).
    Labeled figure: (0,20)-(40,20)-(78,35)-(130,35), base y=0. Deep circle
    tangent to the base: Slide Bishop 1.412 / Spencer 1.400; D&W 1.40. The
    shallow infinite-slope mechanism is tan(30)/tan(21.5) = 1.46."""
    sd = load_slope_data(ACADS_1A)
    m0 = dict(sd['materials'][0]); m1 = dict(sd['materials'][0])
    m0.update(name='Embankment', c=0.0, phi=30.0, gamma=120.0, option='mc', u='none')
    m1.update(name='Foundation', c=450.0, phi=0.0, gamma=120.0, option='mc', u='none')
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 20.0), (40.0, 20.0), (78.0, 35.0), (130.0, 35.0)]},
        {'mat_id': 1, 'coords': [(0.0, 20.0), (130.0, 20.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 60.0, 'Yo': 60.0, 'Depth': 0.0, 'R': 60.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp079.xlsx'))
    return 'vp079.xlsx'


def vp081():
    """Slide #81 / Duncan & Wright (2005) Fig. 14.7: embankment (c=0,
    phi=30, gamma=124) on a phi=0 foundation (c=500 psf, gamma=98, 15 ft).
    Labeled figure: (0,15)-(35,15)-(73,34)-(128,34), base y=0. Deep circle
    tangent to the base: Slide Bishop 1.230 / Spencer 1.209; D&W 1.21."""
    sd = load_slope_data(ACADS_1A)
    m0 = dict(sd['materials'][0]); m1 = dict(sd['materials'][0])
    m0.update(name='Embankment', c=0.0, phi=30.0, gamma=124.0, option='mc', u='none')
    m1.update(name='Foundation', c=500.0, phi=0.0, gamma=98.0, option='mc', u='none')
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 15.0), (35.0, 15.0), (73.0, 34.0), (128.0, 34.0)]},
        {'mat_id': 1, 'coords': [(0.0, 15.0), (128.0, 15.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 55.0, 'Yo': 55.0, 'Depth': 0.0, 'R': 55.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp081.xlsx'))
    return 'vp081.xlsx'


def vp061a():
    """Slide #61, power-curve case / Baker (2003) ex. 3: London clay
    tau = 3.39344*(sigma'+0.152)^0.6 (Baker A=0.535, n=0.6, T=0.0015) on
    the same 43-deg, H=6 slope as #44. Slide: Janbu simplified 1.348,
    Spencer 1.468; Baker's non-linear solution 1.48 (sigma_max 21.4)."""
    sd = _baker1_slope_data()
    sd['materials'][0].update(name='London clay (power)', c=0.0, phi=0.0,
                              gamma=18.0, option='pow', pow_a=3.39344,
                              pow_b=0.6, pow_c=0.0, pow_d=0.152, u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp061a.xlsx'))
    return 'vp061a.xlsx'


def vp061b():
    """Slide #61, Mohr-Coulomb case: the fitted M-C envelope c'=6.0,
    phi'=32 (Baker's London-clay fit). Slide: Janbu simplified 1.291,
    Spencer 1.366; Baker 1.35 (sigma_max 27.5)."""
    sd = _baker1_slope_data()
    sd['materials'][0].update(name='London clay (MC)', c=6.0, phi=32.0,
                              gamma=18.0, option='mc', u='none')
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp061b.xlsx'))
    return 'vp061b.xlsx'


def vp025():
    """Slide #25 / Chen & Shao (1988): Prandtl's weightless bearing-capacity
    mechanism on a 60-deg, 10-m slope (c=49 kPa, phi=0, gamma=1e-6), strip
    load q=149.31 kPa over 10 m on the crest. Surface = active 45-deg wedge
    from the load's right edge, circular fan centered on the load's left
    edge (R=7.071, tangent to both straight segments), exit through the
    face at Slide's printed (0.773, 1.340). Theoretical FS=1.0; Slide
    Spencer 1.051 (Chen & Shao 1.05); their M-P 1.009 uses a custom f(x)."""
    import numpy as np
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='weightless clay', c=49.0, phi=0.0, gamma=1e-6, option='mc', u='none')
    sd['materials'] = [m]
    xc = 10.0 / np.tan(np.radians(60.0))   # crest edge x = 5.7735
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-3.0, 0.0), (0.0, 0.0), (xc, 10.0), (20.0, 10.0)]},
    ]
    sd['max_depth'] = -2.0
    sd['gamma_water'] = 9.81
    sd['dloads'] = [[
        {'X': xc, 'Y': 10.0, 'Normal': 149.31},
        {'X': xc + 10.0, 'Y': 10.0, 'Normal': 149.31},
    ]]
    sd['circular'] = False
    sd['circles'] = []
    # Prandtl surface: exit point (printed), tangent line to the fan arc,
    # arc to the apex (10.774, 5), 45-deg line to the load's right edge
    C = np.array([xc, 10.0]); R = 10.0 / np.sqrt(2.0)
    P = np.array([0.773, 1.340])
    d = C - P; L2 = d @ d; t2 = L2 - R * R
    # lower tangent point from external point P
    a = t2 / L2; b = R * np.sqrt(t2) / L2
    T1 = P + a * d + b * np.array([d[1], -d[0]])
    ang1 = np.arctan2(*(T1 - C)[::-1])
    ang2 = np.arctan2(-5.0, 5.0)   # apex (xc+5, 5)
    pts = [(-0.05, 2.05), tuple(P)]
    for t in np.linspace(ang1, ang2, 13):
        pts.append((C[0] + R * np.cos(t), C[1] + R * np.sin(t)))
    pts += [(xc + 10.0, 10.0), (xc + 10.3, 10.3)]
    sd['non_circ'] = [{'X': round(x, 4), 'Y': round(y, 4), 'Movement': 'Free'} for x, y in pts]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp025.xlsx'))
    return 'vp025.xlsx'


def vp027():
    """Slide #27 / XSTABL v5 manual (Sharma 1996) via Malkawi et al. (2001):
    two-material slope over undulating bedrock (polygon-mode bottom), zero-
    strength cap layer, water table, soil 1 with distinct moist/saturated
    unit weights (116.4/124.2 pcf). Specified circle (59.52, 219.21,
    R=157.68): Slide Bishop 1.396 / Spencer 1.402; XSTABL 1.397 / 1.403.
    NOTE: Slide and XSTABL apply the phreatic-inclination (Hu) correction
    (u reduced by cos^2 of the phreatic slope); xslope uses the static
    vertical head, so its pore pressures are slightly higher. The water
    table was pixel-traced from the labeled figure (ground-coincident to
    x=63, then departing below the ground line)."""
    from shapely.geometry import Polygon
    from xslope.fileio import build_ground_surface_from_polygons
    sd = load_slope_data(LEVEE_POLY)
    base = dict(sd['materials'][0])
    m1 = dict(base); m2 = dict(base)
    m1.update(name='Soil 1', c=500.0, phi=14.0, gamma=116.4, gamma_sat=124.2,
              option='mc', u='piezo')
    m2.update(name='Soil 2 (zero strength)', c=0.0, phi=0.0, gamma=116.4,
              gamma_sat=116.4, option='mc', u='piezo')
    sd['materials'] = [m1, m2]
    bedrock = [(0.0, 15.0), (29.0, 24.0), (51.0, 26.0), (78.0, 56.0),
               (94.0, 65.0), (113.0, 64.0), (133.0, 56.0), (161.0, 58.0),
               (200.0, 76.0)]
    ground_lo = [(0.0, 68.0), (22.0, 67.0), (38.0, 63.0), (63.0, 73.0), (101.0, 88.0)]
    zones = [
        (0, bedrock + [(200.0, 99.0)] + ground_lo[::-1]),
        (1, [(101.0, 88.0), (200.0, 99.0), (200.0, 110.0), (138.0, 103.0)]),
    ]
    sd['polygons'] = [{'polygon': Polygon(p), 'mat_id': mid} for mid, p in zones]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    sd['profile_lines'] = []
    sd['max_depth'] = 15.0
    sd['gamma_water'] = 62.4
    sd['seepage_bc'] = {'specified_heads': [], 'exit_face': []}
    sd['dloads'] = []; sd['line_loads'] = []; sd['pile_lines'] = []
    sd['reinforcement_lines'] = []; sd['reinforce_lines'] = []
    sd['piezo_phreatic'] = True   # XSTABL/Slide phreatic (Hu) correction
    sd['piezo_line'] = [(0.0, 68.0), (22.0, 67.0), (38.0, 63.0), (63.0, 73.0),
                        (80.0, 78.6), (95.0, 81.8), (110.0, 84.5), (125.0, 86.9),
                        (155.0, 90.2), (185.0, 93.1), (200.0, 94.5)]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 59.52, 'Yo': 219.21, 'Depth': 219.21 - 157.68, 'R': 157.68}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp027.xlsx'))
    return 'vp027.xlsx'


def _dw145_slope_data():
    """Duncan & Wright (2005) Fig. 14.5 (Slide #80): embankment (c=1,
    phi=35, gamma=120) over five flat foundation layers alternating
    phi=0 clays and c~0 sands (950/0/110, 1/32/122, 500/0/98, 1/37/131,
    600/0/103), all labeled in Slide's figure. Imperial units."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Embankment', 1.0, 35.0, 120.0), ('Foundation I', 950.0, 0.0, 110.0),
             ('Foundation II', 1.0, 32.0, 122.0), ('Foundation III', 500.0, 0.0, 98.0),
             ('Foundation IV', 1.0, 37.0, 131.0), ('Foundation V', 600.0, 0.0, 103.0)]
    sd['materials'] = []
    for name, c, phi, gamma in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gamma, option='mc', u='none')
        sd['materials'].append(m)
    tops = [
        [(0.0, 60.0), (80.0, 60.0), (205.0, 110.0), (330.0, 110.0)],
        [(0.0, 60.0), (330.0, 60.0)],
        [(0.0, 45.0), (330.0, 45.0)],
        [(0.0, 35.0), (330.0, 35.0)],
        [(0.0, 25.0), (330.0, 25.0)],
        [(0.0, 20.0), (330.0, 20.0)],
    ]
    sd['profile_lines'] = [{'mat_id': i, 'coords': c} for i, c in enumerate(tops)]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['non_circ'] = []
    return sd


def vp080a():
    """Slide #80 case 1: circle from (142,147) tangent to the foundation
    top (R=87). Slide Bishop 2.549 / Spencer 2.545; D&W 2.56."""
    sd = _dw145_slope_data()
    sd['circles'] = [{'Xo': 142.0, 'Yo': 147.0, 'Depth': 60.0, 'R': 87.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp080a.xlsx'))
    return 'vp080a.xlsx'


def vp080b():
    """Slide #80 case 2: circle tangent to the 15-ft-depth line (R=102).
    Slide Bishop 1.398 / Spencer 1.359; D&W 1.35."""
    sd = _dw145_slope_data()
    sd['circles'] = [{'Xo': 142.0, 'Yo': 147.0, 'Depth': 45.0, 'R': 102.0}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp080b.xlsx'))
    return 'vp080b.xlsx'


def vp074():
    """Slide #74 / Duncan & Wright (2005) Fig. 7.12: 100-ft sand embankment
    (c=0, phi=40, gamma=140) on a 50-ft saturated clay foundation (c=2500,
    phi=0, gamma=140). Labeled figure: (100,50)-(300,150) crest to (400,150),
    down to (600,50); base y=0. Circular search: D&W Bishop 1.22 / Spencer
    1.19; Slide 1.228 / 1.201 (Janbu simplified 1.079, D&W 1.07)."""
    sd = load_slope_data(ACADS_1A)
    m0 = dict(sd['materials'][0]); m1 = dict(sd['materials'][0])
    m0.update(name='Sand', c=0.0, phi=40.0, gamma=140.0, option='mc', u='none')
    m1.update(name='Saturated clay', c=2500.0, phi=0.0, gamma=140.0, option='mc', u='none')
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 50.0), (100.0, 50.0), (300.0, 150.0),
                                 (400.0, 150.0), (600.0, 50.0), (700.0, 50.0)]},
        {'mat_id': 1, 'coords': [(0.0, 50.0), (700.0, 50.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 250.0, 'Yo': 260.0, 'Depth': 10.0, 'R': 250.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp074.xlsx'))
    return 'vp074.xlsx'


def vp064():
    """Slide #64 / USACE EM 1110-2-1902 Figure 4-1: end-of-construction
    Spencer hand-check dam. Symmetric 50-ft embankment at 4H:1V (crest
    x=-5..5, toes at +-205), 10-ft sand blanket over foundation clay
    (-10..-37) and rock (-37..-40), embankment core trench through the sand
    at x=-17..17, groundwater at the sand top (el 0, dipping through the
    trench), 7-ft crest tension crack. Undrained strengths: embankment
    1000/5 (115/120 pcf), sand 0/35 (125/130), clay 3000/0 (110/115), rock
    0/45 (160/165). Specified circle: center (102, 163), tangent to el 0
    (R=163). USACE Spencer 2.44; Slide Bishop 2.447 / Spencer 2.445 /
    Janbu corrected 2.430. xslope Spencer 2.488 (+1.8%; Slide's figure is
    vertex-unlabeled, so the crest/toe placement carries a small
    uncertainty)."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Embankment', 1000.0, 5.0, 115.0, 120.0),
             ('Sand', 0.0, 35.0, 125.0, 130.0),
             ('Foundation Clay', 3000.0, 0.0, 110.0, 115.0),
             ('Rock', 0.0, 45.0, 160.0, 165.0)]
    sd['materials'] = []
    for name, c, phi, gm, gs in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gm, gamma_sat=gs,
                 option='mc', u='piezo')
        sd['materials'].append(m)
    # crest half-width 17 / toes at 217 (4H:1V): pinned by reconciling
    # USACE's Fig 4-1 slice table (slice 1 width 23 ft, avg height 16 ft;
    # total slice span 173 ft matches the crack-truncated arc)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-225.0, 0.0), (-217.0, 0.0), (-17.0, 50.0),
                                 (17.0, 50.0), (217.0, 0.0), (225.0, 0.0)]},
        {'mat_id': 1, 'coords': [(-225.0, 0.0), (-17.0, 0.0), (-8.0, -10.0),
                                 (8.0, -10.0), (17.0, 0.0), (225.0, 0.0)]},
        {'mat_id': 2, 'coords': [(-225.0, -10.0), (225.0, -10.0)]},
        {'mat_id': 3, 'coords': [(-225.0, -37.0), (225.0, -37.0)]},
    ]
    sd['max_depth'] = -40.0
    sd['gamma_water'] = 62.4
    sd['tcrack_depth'] = 7.0
    sd['piezo_line'] = [(-225.0, 0.0), (-17.0, 0.0), (-8.0, -10.0),
                        (8.0, -10.0), (17.0, 0.0), (225.0, 0.0)]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 102.0, 'Yo': 163.0, 'Depth': 0.0, 'R': 163.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp064.xlsx'))
    return 'vp064.xlsx'


def vp067():
    """Slide #67 / USACE EM 1110-2-1902 example F-5: end-of-construction
    embankment on an undrained fine-grained foundation. Labeled figure:
    ground (-100,100)-(0,100)-(174,158)-(257,191)-(301,191)-(400,150),
    foundation top el 100, base el 0. Embankment c=1780 psf, phi=5,
    gamma=135; foundation c=1600, phi=2, gamma=127. Circle centered 259 ft
    above / 101 ft right of the toe (0,100), through the toe (R=278.0).
    USACE 1.33; Slide Bishop 1.332 / Spencer 1.328 / Janbu corr 1.345."""
    import math
    sd = load_slope_data(ACADS_1A)
    m0 = dict(sd['materials'][0]); m1 = dict(sd['materials'][0])
    m0.update(name='Embankment', c=1780.0, phi=5.0, gamma=135.0, option='mc', u='none')
    m1.update(name='Foundation', c=1600.0, phi=2.0, gamma=127.0, option='mc', u='none')
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-100.0, 100.0), (0.0, 100.0), (174.0, 158.0),
                                 (257.0, 191.0), (301.0, 191.0), (400.0, 150.0)]},
        {'mat_id': 1, 'coords': [(-100.0, 100.0), (400.0, 100.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    R = math.hypot(101.0, 259.0)
    sd['circles'] = [{'Xo': 101.0, 'Yo': 359.0, 'Depth': 359.0 - R, 'R': R}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp067.xlsx'))
    return 'vp067.xlsx'


def vp065():
    """Slide #65 / USACE EM 1110-2-1902 Figure 4-2: the #64 dam under
    upstream low-pool steady conditions (pool el 20, no tension crack),
    drained strengths: embankment c=100 psf phi=25 (115/120 pcf), sand
    0/35 (125/130), clay 0/28 (110/115), rock 0/45 (160/165). Printed
    circle center (-102, 163), R=173 (tangent to the clay top at el -10).
    USACE Bishop 2.71; Slide Bishop 2.716 / Spencer 2.736 / GLE 2.744."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Embankment', 100.0, 25.0, 115.0, 120.0),
             ('Sand', 0.0, 35.0, 125.0, 130.0),
             ('Foundation Clay', 0.0, 28.0, 110.0, 115.0),
             ('Rock', 0.0, 45.0, 160.0, 165.0)]
    sd['materials'] = []
    for name, c, phi, gm, gs in props:
        m = dict(base)
        m.update(name=name, c=c, phi=phi, gamma=gm, gamma_sat=gs,
                 option='mc', u='piezo')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-225.0, 0.0), (-217.0, 0.0), (-17.0, 50.0),
                                 (17.0, 50.0), (217.0, 0.0), (225.0, 0.0)]},
        {'mat_id': 1, 'coords': [(-225.0, 0.0), (-17.0, 0.0), (-8.0, -10.0),
                                 (8.0, -10.0), (17.0, 0.0), (225.0, 0.0)]},
        {'mat_id': 2, 'coords': [(-225.0, -10.0), (225.0, -10.0)]},
        {'mat_id': 3, 'coords': [(-225.0, -37.0), (225.0, -37.0)]},
    ]
    sd['max_depth'] = -40.0
    sd['gamma_water'] = 62.4
    sd['piezo_line'] = [(-225.0, 20.0), (225.0, 20.0)]
    gw = 62.4
    # ponded water on the submerged upstream face (pool el 20)
    sd['dloads'] = [[
        {'X': -225.0, 'Y': 0.0, 'Normal': gw * 20.0},
        {'X': -217.0, 'Y': 0.0, 'Normal': gw * 20.0},
        {'X': -137.0, 'Y': 20.0, 'Normal': 0.0},
    ]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': -102.0, 'Yo': 163.0, 'Depth': -10.0, 'R': 173.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp065.xlsx'))
    return 'vp065.xlsx'


def _dw627_slope_data(pool):
    """Duncan & Wright (2005) Fig. 6.27 (Slide #70): fully submerged
    homogeneous slope (c=100 psf, phi=20, gamma=128), ground
    (0,15)-(30,15)-(105,45)-(140,45), base el 0, pool at `pool` (30 or 60 ft
    above the crest). Pond load on the whole submerged surface; FS should be
    independent of the submergence depth (D&W 1.60)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='Material 1', c=100.0, phi=20.0, gamma=128.0, option='mc', u='piezo')
    sd['materials'] = [m]
    ground = [(0.0, 15.0), (30.0, 15.0), (105.0, 45.0), (140.0, 45.0)]
    sd['profile_lines'] = [{'mat_id': 0, 'coords': ground}]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['piezo_line'] = [(0.0, pool), (140.0, pool)]
    gw = 62.4
    sd['dloads'] = [[{'X': x, 'Y': y, 'Normal': gw * (pool - y)} for x, y in ground]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 60.0, 'Yo': 80.0, 'Depth': 5.0, 'R': 70.0}]
    sd['non_circ'] = []
    return sd


def vp070a():
    """Slide #70 case 1 (pool 30 ft above the crest, el 75). Slide circular
    Bishop 1.603 / Spencer 1.599; D&W 1.60."""
    sd = _dw627_slope_data(75.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp070a.xlsx'))
    return 'vp070a.xlsx'


def vp070b():
    """Slide #70 case 2 (pool 60 ft above the crest, el 105): identical FS
    to case 1 — the submergence-independence demonstration. Slide Bishop
    1.603 / Spencer 1.599; D&W 1.60."""
    sd = _dw627_slope_data(105.0)
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp070b.xlsx'))
    return 'vp070b.xlsx'


def vp066():
    """Slide #66 / USACE EM 1110-2-1902 Figure 4-3: the #64/#65 dam with the
    chart-check property set (single unit weights): embankment c=200 psf
    phi=25 gamma=115; sand 0/35/130; clay 0/27/115; rock 0/45/160. Pool at
    el 20, printed circle center (-135, 169), R=169 (tangent to el 0).
    USACE 2.30; Slide Bishop 2.307 / Spencer 2.307 / Janbu corr 2.290."""
    sd = vp065.__wrapped__() if hasattr(vp065, '__wrapped__') else None
    sd = load_slope_data(os.path.join(OUT, 'vp065.xlsx'))
    props = [(200.0, 25.0, 115.0), (0.0, 35.0, 130.0), (0.0, 27.0, 115.0), (0.0, 45.0, 160.0)]
    for m, (c, phi, gamma) in zip(sd['materials'], props):
        m.update(c=c, phi=phi, gamma=gamma, gamma_sat=gamma)
    # Slide's own VP66 geometry (its printed slip endpoints prove a
    # -222 toe / -15 crest-edge face, slope 1:4.14 - different from the
    # #64/#65 models); R nudged +0.1 past the exact crest-corner tangency
    sd['profile_lines'][0]['coords'] = [(-225.0, 0.0), (-222.0, 0.0), (-15.0, 50.0),
                                        (15.0, 50.0), (222.0, 0.0), (225.0, 0.0)]
    sd['dloads'] = [[{'X': -225.0, 'Y': 0.0, 'Normal': 62.4 * 20.0},
                     {'X': -222.0, 'Y': 0.0, 'Normal': 62.4 * 20.0},
                     {'X': -139.2, 'Y': 20.0, 'Normal': 0.0}]]
    sd['circles'] = [{'Xo': -135.0, 'Yo': 169.0, 'Depth': 169.0 - 169.1, 'R': 169.1}]
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp066.xlsx'))
    return 'vp066.xlsx'


def vp068():
    """Slide #68 / USACE EM 1110-2-1902 example E-10: undrained (phi=0)
    three-layer slope with 8 ft of ponded water (pool el 0). Labeled figure:
    ground (-10,-20)-(0,-8)-(40,-8)-(50,4)-(60,16)-(100,16); soil 1
    (c=600, gamma=120) above el 4, soil 2 (400/100) to el -8, soil 3
    (500/105) to the base at el -20. Circle center 8.4 right / 36 above the
    toe (40,-8) -> (48.4, 28), tangent to the base (R=48).
    Slide Bishop/Spencer 1.241; USACE E-10 chart solution 1.33."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    props = [('Soil 1', 600.0, 120.0), ('Soil 2', 400.0, 100.0), ('Soil 3', 500.0, 105.0)]
    sd['materials'] = []
    for name, c, gamma in props:
        m = dict(base)
        m.update(name=name, c=c, phi=0.0, gamma=gamma, option='mc', u='piezo')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-10.0, -20.0), (0.0, -8.0), (40.0, -8.0),
                                 (50.0, 4.0), (60.0, 16.0), (100.0, 16.0)]},
        {'mat_id': 1, 'coords': [(-10.0, -20.0), (0.0, -8.0), (40.0, -8.0),
                                 (50.0, 4.0), (100.0, 4.0)]},
        {'mat_id': 2, 'coords': [(-10.0, -20.0), (0.0, -8.0), (100.0, -8.0)]},
    ]
    sd['max_depth'] = -20.0
    sd['gamma_water'] = 62.4
    sd['piezo_line'] = [(-10.0, 0.0), (100.0, 0.0)]
    gw = 62.4
    sd['dloads'] = [[
        {'X': -10.0, 'Y': -20.0, 'Normal': gw * 20.0},
        {'X': 0.0, 'Y': -8.0, 'Normal': gw * 8.0},
        {'X': 40.0, 'Y': -8.0, 'Normal': gw * 8.0},
        {'X': 46.667, 'Y': 0.0, 'Normal': 0.0},
    ]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 48.4, 'Yo': 28.0, 'Depth': -20.0, 'R': 48.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp068.xlsx'))
    return 'vp068.xlsx'


def vp069():
    """Slide #69 / USACE EM 1110-2-1902 example F-6: steady-seepage embankment
    on a granular foundation. Geometry read off Slide's Figure 69.1 (the
    figure has axis ticks and vertex markers, so the vertices are exact):
    ground (0,100)-(38.4,100)-(60.8,112)-(81,112)-(194.9,73.7)-(400,0)-(450,0),
    foundation el 0 down to the model base el -75. Embankment c'=0, phi'=34,
    gamma=130; foundation c'=0, phi'=35, gamma=125 (total unit weights). The
    piezometric line runs from the pool surface (el 100) down to the chimney
    drain and out to the tailwater at el 22.5; USACE takes u = gamma_w times
    the vertical distance to that line, so the line is a plain piezo line (no
    cos^2 correction). Tailwater ponds the toe from x=337.4 to the right edge.
    Circle: center (269, 248), R = 280 (131 ft left of / 248 ft above the toe;
    it bottoms out exactly on USACE's el -32 stratum line).
    USACE 2.01; Slide Bishop 2.011 / Spencer 2.026 / GLE 2.027 / Janbu corr 1.830.

    The small zones in the real section (rip-rap, chimney drain, drainage
    blanket) are all given the embankment properties, as USACE does -- the
    circle misses them anyway. Slide's own section flattens the upstream face
    to a bench at the pool elevation; the circle enters on that bench, so
    everything further upstream is outside the sliding mass."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    sd['materials'] = []
    for name, phi, gamma in [('Embankment', 34.0, 130.0), ('Foundation', 35.0, 125.0)]:
        m = dict(base)
        m.update(name=name, c=0.0, phi=phi, gamma=gamma, gamma_sat=gamma,
                 option='mc', u='piezo')
        sd['materials'].append(m)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 100.0), (38.4, 100.0), (60.8, 112.0),
                                 (81.0, 112.0), (194.9, 73.7), (400.0, 0.0),
                                 (450.0, 0.0)]},
        {'mat_id': 1, 'coords': [(0.0, 0.0), (450.0, 0.0)]},
    ]
    sd['max_depth'] = -75.0
    sd['gamma_water'] = 62.4
    sd['piezo_line'] = [(0.0, 100.0), (38.4, 100.0), (121.1, 69.0),
                        (179.4, 22.5), (450.0, 22.5)]
    sd['piezo_phreatic'] = False
    gw = 62.4
    sd['dloads'] = [[{'X': 337.4, 'Y': 22.5, 'Normal': 0.0},
                     {'X': 400.0, 'Y': 0.0, 'Normal': gw * 22.5},
                     {'X': 450.0, 'Y': 0.0, 'Normal': gw * 22.5}]]
    sd['circular'] = True
    sd['circles'] = [{'Xo': 269.0, 'Yo': 248.0, 'Depth': -32.0, 'R': 280.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp069.xlsx'))
    return 'vp069.xlsx'


def _vp071_slope_data():
    """Common geometry for Slide #71 / D&W Figs. 6.37-6.38: a homogeneous
    2H:1V slope, ground (0,40)-(120,40)-(200,80)-(440,80) over a base at el 0,
    c'=200 psf, phi'=20 deg, gamma=125 pcf. Water stands at el 40 on the left
    (toe) boundary and el 75 on the right."""
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name='Material 1', c=200.0, phi=20.0, gamma=125.0, gamma_sat=125.0,
             option='mc', u='piezo', k1=1.0, k2=1.0, alpha=0.0, kr0=0.001, h0=-1.0)
    sd['materials'] = [m]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 40.0), (120.0, 40.0), (200.0, 80.0),
                                 (440.0, 80.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 190.0, 'Yo': 140.0, 'Depth': 20.0, 'R': 120.0}]
    sd['non_circ'] = []
    return sd


def vp071a():
    """Slide #71 case 1 / D&W Fig. 6.37: pore pressures from a finite-element
    seepage analysis (u='seep'). Heads: 40 ft on the left boundary (x=0, the
    toe-side water level), 75 ft on the right boundary (x=440); the ground
    surface from the left edge up over the slope face is an exit face.
    Slide Bishop/Spencer/GLE 1.141 (circular); D&W 1.138."""
    sd = _vp071_slope_data()
    sd['materials'][0]['u'] = 'seep'
    sd['piezo_line'] = []
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 40.0, 'coords': [(0.0, 0.0), (0.0, 40.0)]},
            {'head': 75.0, 'coords': [(440.0, 0.0), (440.0, 75.0)]},
        ],
        'exit_face': [(0.0, 40.0), (120.0, 40.0), (200.0, 80.0), (440.0, 80.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp071a.xlsx'))
    return 'vp071a.xlsx'


def vp071b():
    """Slide #71 case 2 / D&W Fig. 6.38: the same slope with the piezometric
    line approximation. The line was read off Slide's Figure 71.2 vertex
    markers: (0,40)-(120,40)-(127.4,43.2)-(147.4,46.2)-(181.7,50.7)-
    (232.3,56.7)-(326.8,65.7)-(440,75).
    Slide Bishop/Spencer 1.142, GLE 1.141 (circular); D&W 1.141."""
    sd = _vp071_slope_data()
    sd['piezo_line'] = [(0.0, 40.0), (120.0, 40.0), (127.4, 43.2), (147.4, 46.2),
                        (181.7, 50.7), (232.3, 56.7), (326.8, 65.7), (440.0, 75.0)]
    sd['piezo_phreatic'] = False
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp071b.xlsx'))
    return 'vp071b.xlsx'


def _vp083_slope_data(fnd):
    """Slide #83 / D&W (2005) Fig. 14.20-b: an embankment wall, ground
    (0,40)-(55,40)-(75,30)-(140,30), foundation el 30 down to a base at el 0.
    Embankment c'=0, phi'=36, gamma=123 pcf. ``fnd`` supplies the foundation's
    undrained strength profile."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    m0 = dict(base); m0.update(name='Embankment', c=0.0, phi=36.0, gamma=123.0,
                               gamma_sat=123.0, option='mc', u='none')
    m1 = dict(base); m1.update(name='Foundation', phi=0.0, gamma=97.0,
                               gamma_sat=97.0, u='none', **fnd)
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 40.0), (55.0, 40.0), (75.0, 30.0), (140.0, 30.0)]},
        {'mat_id': 1, 'coords': [(0.0, 30.0), (140.0, 30.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 65.0, 'Yo': 80.0, 'Depth': 5.0, 'R': 75.0}]
    sd['non_circ'] = []
    return sd


def vp083a():
    """Slide #83 case 1: foundation strength profile I, cu = 200 + 15*depth psf
    below the foundation surface (el 30) -- XSLOPE's 'cp' option with c=200,
    cp=15, r_elev=30. Circular search.
    Slide Bishop 1.313 / Spencer 1.285 / GLE 1.294; D&W average 1.300."""
    sd = _vp083_slope_data(dict(option='cp', c=200.0, cp=15.0, r_elev=30.0))
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp083a.xlsx'))
    return 'vp083a.xlsx'


def vp083b():
    """Slide #83 case 2: foundation strength profile II, cu = 300 psf constant.
    The critical circle runs down to the base of the foundation.
    Slide Bishop 1.335 / Spencer 1.330 / GLE 1.331; D&W average 1.312."""
    sd = _vp083_slope_data(dict(option='mc', c=300.0))
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'vp083b.xlsx'))
    return 'vp083b.xlsx'


def _vp084(cz, tag):
    """Slide #84 / D&W (2005) Fig. 15.9: an earth embankment, ground
    (0,20)-(40,20)-(90,40)-(140,40), foundation el 20 down to a base at el 0.
    Embankment c'=0, phi'=35, gamma=125 pcf; foundation phi=0, gamma=100 pcf,
    cu = 300 + cz*z with z the depth below the foundation surface -- XSLOPE's
    'cp' option with c=300, cp=cz, r_elev=20. Four strength profiles,
    cz = 0/5/10/15 psf/ft. Circular search."""
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    m0 = dict(base); m0.update(name='Embankment', c=0.0, phi=35.0, gamma=125.0,
                               gamma_sat=125.0, option='mc', u='none')
    m1 = dict(base); m1.update(name='Foundation', c=300.0, cp=cz, r_elev=20.0,
                               phi=0.0, gamma=100.0, gamma_sat=100.0,
                               option='cp', u='none')
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 20.0), (40.0, 20.0), (90.0, 40.0), (140.0, 40.0)]},
        {'mat_id': 1, 'coords': [(0.0, 20.0), (140.0, 20.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 62.4
    sd['circular'] = True
    sd['circles'] = [{'Xo': 65.0, 'Yo': 70.0, 'Depth': 5.0, 'R': 65.0}]
    sd['non_circ'] = []
    save_slope_data_to_xlsx(sd, os.path.join(OUT, f'vp084{tag}.xlsx'))
    return f'vp084{tag}.xlsx'


def vp084a():
    """Profile I (cz=0): Slide Bishop 0.761 / Spencer 0.756 / GLE 0.762; D&W 0.75."""
    return _vp084(0.0, 'a')


def vp084b():
    """Profile II (cz=5): Slide Bishop 0.909 / Spencer 0.898 / GLE 0.908; D&W 0.90."""
    return _vp084(5.0, 'b')


def vp084c():
    """Profile III (cz=10): Slide Bishop 1.045 / Spencer 1.032 / GLE 1.034; D&W 1.03."""
    return _vp084(10.0, 'c')


def vp084d():
    """Profile IV (cz=15): Slide Bishop 1.154 / Spencer 1.134 / GLE 1.138; D&W 1.13."""
    return _vp084(15.0, 'd')


BUILDERS = [vp002, vp003, vp004, vp005, vp006, vp008, vp009, vp015, vp016, vp017, vp018, vp019, vp020, vp021a, vp021b, vp023, vp024, vp025, vp027, vp036, vp041, vp043, vp044a, vp044b, vp044c, vp061a, vp061b, vp064, vp065, vp066, vp067, vp068, vp069, vp070a, vp070b, vp071a, vp071b, vp083a, vp083b, vp084a, vp084b, vp084c, vp084d, vp045a, vp045b, vp047, vp048, vp050, vp051, vp052a, vp052b, vp054a, vp054b, vp062a, vp062b, vp074, vp078, vp079, vp080a, vp080b, vp081, vp085a, vp085b, vp086, vp087, vp088, vp089, vp090, vp091, vp092, vp093, vp094, vp096, vp098, vp099, vp097, vp100, vp101]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
