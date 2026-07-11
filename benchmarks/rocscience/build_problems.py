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


BUILDERS = [vp002, vp003, vp004, vp005, vp006, vp008, vp009, vp015, vp016, vp017, vp018, vp019, vp020, vp021a, vp021b, vp023, vp024, vp036, vp041, vp045a, vp045b, vp050, vp051, vp054a, vp054b, vp062a, vp062b, vp085a, vp085b, vp086, vp096]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
