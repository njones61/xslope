"""Build the FHWA Appendix E steel-reinforcement MSE wall models, E3 to E7.

Berg, R.R., Christopher, B.R., and Samtani, N.C. (2009), "Design and
Construction of Mechanically Stabilized Earth Walls and Reinforced Soil
Slopes", FHWA-NHI-10-025 (Volume II), Appendix E, Examples E3 to E7.

Example E1, already in the corpus, is a geogrid wall: a continuous sheet whose
pullout law is stated per unit width of wall.  The five examples built here are
the manual's steel-reinforcement designs, and every one of them reinforces the
wall with DISCRETE elements at a horizontal spacing:

  E3  segmental precast panel wall, 2H:1V sloping backfill, ribbed steel strips
  E4  segmental precast panel wall, level backfill and a live load surcharge,
      steel bar mats
  E5  bridge abutment on a spread footing on top of a panel wall, ribbed strips
  E6  the E4 wall re-checked for a traffic barrier impact (Extreme Event II)
  E7  the E4 wall re-checked for earthquake loading (Extreme Event I)

HOW A DISCRETE ELEMENT IS ENTERED
---------------------------------
FHWA states the nominal pullout resistance of one element as

    Pr = F* alpha (2b) Le sigma'v

where b is the element's bearing width -- the strip width for E3 and E5, the
bar mat width for E4, E6 and E7 -- and 2b carries the two bearing faces.
XSLOPE states the same resistance as a rate per unit length of line,

    r(s) = 2 (Adhesion + sigma'v(s) tan(Delta)) / Spacing

and integrates it along the embedment.  The two are the same statement with
Adhesion = 0,

    Delta = arctan(F* alpha b)     (b in feet)

and Spacing set to the element's own horizontal spacing, so that everything
entered on the reinforce sheet is the manual's own per-element number and the
loader divides it once, on the way in, to reach the per-unit-width-of-wall
convention every engine downstream works in.  XSLOPE's factor of two carries
FHWA's C = 2, and the bearing width rides in Delta because that is where a
coefficient on the overburden term belongs -- it is the same transformation a
coverage ratio takes, Delta = arctan(F* alpha Rc), read on the element rather
than on the wall.

F* varies with depth for steel reinforcement.  A reinforcement layer in an MSE
wall is horizontal, so its depth is a single number and its F* is a single
number: one Delta per layer, taken from the manual's own F* column.

Tend1 = Tmax at the face end.  FHWA counts only the resisting-zone side of each
layer; the front side is carried by the facing connection and is checked
separately.  Setting the face-end anchorage to the layer's own tensile strength
makes the face branch of the envelope non-governing, so what the resisting side
develops is what the envelope reports.  XSLOPE does not model the
panel-to-reinforcement connection, so no connection capacity is entered.

The reinforcement Type vocabulary on the reinforce sheet has no steel-strip or
bar-mat entry, so these lines are entered as generic tensile lines.  Dir is set
to axial rather than the generic tangent: inextensible steel reinforcement
carries its tension along the strip, which is the direction FHWA computes Tmax
in.

WATER
-----
None of the five examples has a water table -- each states level ground in front
of the wall and a foundation soil with no water table -- so no piezometric line
is entered and every surcharge is a distributed load, never a unit weight.

Run from the repo root:
    PYTHONPATH=. python3 benchmarks/published/build_fhwa_e3_e7.py
"""
import math
import os

from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'docs', 'verification', 'files', 'published')
SEED = os.path.join(ROOT, 'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

GAMMA_WATER = 62.4

#: Ribbed strip width, ft (1.969 in), Examples E3 Step 7.5 and E5 Step 8.5.
B_STRIP = 0.164
#: Scale correction factor for inextensible reinforcement, AASHTO Table
#: 11.10.6.3.2-1, quoted by every one of these examples.
ALPHA = 1.0


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------

def _material(name, gamma, phi):
    """One drained Mohr-Coulomb soil, every column stated.

    These examples publish a unit weight and a friction angle and nothing else,
    so the cohesion is zero, the uncertainties are zero, and the stiffness,
    tensile cutoff and unsaturated columns are left empty rather than carrying
    invented values.
    """
    return {
        'name': name, 'gamma': gamma, 'gamma_sat': None, 'option': 'mc',
        'c': 0.0, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,
        't_cut': None, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': 'none', 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': 0.0, 'k2': 0.0, 'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0, 'h0': 0.0, 'vg_a': 0.0, 'vg_n': 0.0,
        'E': None, 'nu': None,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _zone(mat_id, coords):
    """One material zone, in the shape the loader hands the rest of XSLOPE."""
    return {'mat_id': mat_id, 'polygon': Polygon(coords), 'size': None}


def _delta(f_star, width):
    """Delta for an element of bearing width ``width`` at pullout factor F*.

    ``arctan(F* alpha b)`` in degrees: the manual's whole per-element pullout
    coefficient, expressed as the interface friction angle XSLOPE's law reads.
    """
    return math.degrees(math.atan(f_star * ALPHA * width))


def _layer(label, y, length, t_element, f_star, width, spacing):
    """One horizontal reinforcement layer, face to the back of the zone.

    ``t_element`` is the manual's own per-element tensile resistance -- per
    strip, or per bar mat.  ``slope_data`` carries capacities per unit width of
    wall (the loader divides by Spacing on the way in, and the writer multiplies
    back), so it is divided here and the reinforce sheet ends up holding the
    manual's number.
    """
    t_max = t_element / spacing
    return {
        'label': label,
        'x1': 0.0, 'y1': round(y, 4),
        'x2': length, 'y2': round(y, 4),
        # The Type vocabulary has no steel entry; a blank Type is the generic
        # tensile line, and Dir is set to axial because inextensible steel
        # carries its tension along the strip.
        'type': '', 'dir': 'axial', 'appl': 'active',
        't_max': t_max, 't_res': None,
        # Lp1/Lp2 are not read once Adhesion and Delta are filled: the
        # overburden law carries the whole development, so a development length
        # would be a second, contradictory statement of the same bond.
        'lp1': 0.0, 'lp2': 0.0,
        'adhesion': 0.0, 'delta': _delta(f_star, width),
        'tend1': t_max, 'tend2': 0.0,
        'spacing': spacing,
        'E': 0.0, 'area': 0.0,   # no FEM run on these models
    }


def _base(x_left, x_right, y_base, materials):
    """A seeded slope_data with everything these models do not use cleared."""
    sd = load_slope_data(SEED)
    sd['unit_system'] = 'imperial'
    sd['gamma_water'] = GAMMA_WATER
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = 0.0
    sd['k0'] = None          # limit equilibrium only: no in-situ stress state
    sd['materials'] = materials
    sd['profile_lines'] = []
    sd['max_depth'] = y_base
    sd['piezo_line'] = []
    sd['piezo_line2'] = []
    sd['dloads'] = []
    sd['dload_dirs'] = []
    sd['dloads2'] = []
    sd['dload2_dirs'] = []
    sd['reinforce_lines'] = []
    sd['pile_lines'] = []
    sd['line_loads'] = []
    sd['non_circ'] = []
    sd['circular'] = True
    sd['water_loads'] = 'manual'
    sd['element_type'] = None
    sd['target_size'] = None
    sd['search_window'] = {}
    sd['_extents'] = (x_left, x_right, y_base)
    return sd


def _circles(sd, x_o, y_top, height, y_grade):
    """The two starting circles these walls open the search with.

    The center sits over the middle of the reinforced zone at the crest
    elevation plus a quarter of the wall height, one circle tangent to the
    finished grade in front of the wall and one tangent to the top of the
    leveling pad, so the search opens with a surface through the reinforced zone
    and one that undercuts it.  That center is higher than the toe-plus-2H a
    slope would use because the face here is vertical, and on a vertical face
    the arc a high center draws through the toe is far flatter than any
    mechanism the wall has.
    """
    y_o = y_top + height / 4.0
    sd['circles'] = [
        {'Xo': x_o, 'Yo': y_o, 'Depth': y_grade, 'R': y_o - y_grade},
        {'Xo': x_o, 'Yo': y_o, 'Depth': 0.0, 'R': y_o},
    ]


def _write(sd, name):
    del sd['_extents']
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, dest)
    return dest


# ---------------------------------------------------------------------------
# Example E3 -- sloping backfill, ribbed steel strips
# ---------------------------------------------------------------------------
# Steps 1-3 (printed pages E3-1 to E3-2): exposed height He = 28 ft, embedment
# d = 2 ft, so H = 30 ft; reinforcement length L = 0.8H = 24 ft at every level;
# 2H:1V backslope.  Step 2: reinforced fill phi = 34 deg, retained backfill and
# foundation phi = 30 deg, all at 125 pcf, no water table.
E3_H = 30.0
E3_D = 2.0
E3_L = 24.0
E3_TAN_BETA = 0.5
#: The backslope rises L tan(beta) = 12 ft across the reinforced zone, which is
#: the manual's own h = H + L tan(beta) = 42 ft at the back of that zone, and
#: its V2 is the triangle of retained fill standing on the reinforced zone -- a
#: triangle that ends where the reinforcement does.  The model levels the ground
#: at that station.  Nothing in the pullout table reads past it: the deepest
#: overburden any layer averages is Zp at x = L.
E3_Y_CREST = E3_H + E3_L * E3_TAN_BETA      # 42 ft

#: Table E3-7.3 (printed page E3-19), one row per level: depth Z below the top
#: of the wall, the pullout resistance factor F*, the resisting length Le, and
#: the horizontal spacing Sh of the strips the design selects.
E3_LEVELS = [
    # Z,     F*,     Le,     Sh
    (1.25,  1.917, 13.41,  2.50),
    (3.75,  1.751, 13.41,  2.50),
    (6.25,  1.586, 13.41,  2.50),
    (8.75,  1.420, 13.41,  2.50),
    (11.25, 1.254, 13.41,  2.50),
    (13.75, 1.089, 14.25,  2.50),
    (16.25, 0.923, 15.75,  2.50),
    (18.75, 0.757, 17.25,  2.50),
    (21.25, 0.675, 18.75,  1.67),
    (23.75, 0.675, 20.25,  1.67),
    (26.25, 0.675, 21.75,  1.67),
    (28.75, 0.675, 23.25,  1.67),
]
#: Step 7.4 (printed page E3-16): 1.969 in x 0.157 in Grade 65 strip, 75-year
#: design life, Tn = 65 ksi (0.200 in2) = 13.00 k/strip.
E3_TN = 13000.0


def build_fhwa_e3():
    sd = _base(-55.0, 85.0, -25.0,
               [_material('Reinforced fill', 125.0, 34.0),
                _material('Retained backfill', 125.0, 30.0),
                _material('Foundation', 125.0, 30.0)])
    x_left, x_right, y_base = sd['_extents']
    # Polygon geometry: the vertical wall face and the vertical zone boundary at
    # the back of the reinforcement are edges of the zones themselves, which
    # profile lines (one elevation per station) cannot express.
    sd['polygons'] = [
        _zone(0, [(0.0, 0.0), (E3_L, 0.0), (E3_L, E3_Y_CREST), (0.0, E3_H)]),
        _zone(1, [(E3_L, 0.0), (x_right, 0.0), (x_right, E3_Y_CREST),
                  (E3_L, E3_Y_CREST)]),
        _zone(2, [(x_left, E3_D), (0.0, E3_D), (0.0, 0.0), (x_right, 0.0),
                  (x_right, y_base), (x_left, y_base)]),
    ]
    sd['reinforcement_lines'] = [
        _layer('Level %d' % (i + 1), E3_H - z, E3_L, E3_TN, f_star,
               B_STRIP, sh)
        for i, (z, f_star, _le, sh) in enumerate(E3_LEVELS)]
    _circles(sd, round(E3_L / 2.0, 4), E3_Y_CREST, E3_H, E3_D)
    return _write(sd, 'fhwa_e3.xlsx')


# ---------------------------------------------------------------------------
# Example E4 -- level backfill with a live load surcharge, steel bar mats
# ---------------------------------------------------------------------------
# Steps 1-3 (printed pages E4-2 to E4-3): He = 23.64 ft, d = 2 ft, H = 25.64 ft,
# L = 0.7H = 18 ft, level backfill, live load surcharge heq = 2 ft of soil.
# Step 2: reinforced fill phi = 34 deg, retained backfill and foundation
# phi = 30 deg, all at 125 pcf, no water table.
E4_H = 25.64
E4_D = 2.0
E4_L = 18.0
#: Step 2: heq = 2 ft of soil at 125 pcf.
E4_Q = 250.0
#: Step 7.6 (printed page E4-17): longitudinal wires at Sl = 6 in.
E4_SL = 0.5
#: Step 1: 5 ft wide precast panels.
E4_WP = 5.0

#: Table E4-7.4 (printed page E4-18), one row per level: depth Z below the top
#: of the wall, the pullout resistance factor F*, the resisting length Le, the
#: number of longitudinal wires Ng the design selects, and the nominal tensile
#: resistance Tn of one of those wires (7.42 k for W15, 5.17 k for W11,
#: Step 7.4, printed pages E4-14 and E4-15).
E4_LEVELS = [
    # Z,     F*,    Le,     Ng, Tn (lb/wire)
    (1.87,  1.188, 10.31,   4, 5170.0),
    (4.37,  1.110, 10.31,   3, 5170.0),
    (6.87,  1.033, 10.31,   4, 5170.0),
    (9.37,  0.955, 10.31,   4, 5170.0),
    (11.87, 0.438, 10.31,   4, 7420.0),
    (14.37, 0.399, 11.24,   4, 7420.0),
    (16.87, 0.360, 12.74,   4, 7420.0),
    (19.37, 0.214, 14.24,   4, 7420.0),
    (21.87, 0.208, 15.74,   4, 7420.0),
    (24.37, 0.208, 17.24,   4, 7420.0),
]


def _e4_mat_width(n_wires):
    """Bar mat width, ft, for a mat of ``n_wires`` longitudinal wires.

    Step 7.6 sizes the mat from pullout as Np = 1 + (Tmax/Prr)/Sl, which is the
    manual's own statement that a mat of N longitudinal wires spans (N-1)Sl.
    """
    return (n_wires - 1) * E4_SL


def _e4_zones(sd, y_top):
    """The three zones of the E4 wall, up to ``y_top``."""
    x_left, x_right, y_base = sd['_extents']
    return [
        _zone(0, [(0.0, 0.0), (E4_L, 0.0), (E4_L, y_top), (0.0, y_top)]),
        _zone(1, [(E4_L, 0.0), (x_right, 0.0), (x_right, y_top), (E4_L, y_top)]),
        _zone(2, [(x_left, E4_D), (0.0, E4_D), (0.0, 0.0), (x_right, 0.0),
                  (x_right, y_base), (x_left, y_base)]),
    ]


def _e4_base():
    return _base(-55.0, 75.0, -22.0,
                 [_material('Reinforced fill', 125.0, 34.0),
                  _material('Retained backfill', 125.0, 30.0),
                  _material('Foundation', 125.0, 30.0)])


def build_fhwa_e4():
    sd = _e4_base()
    x_left, x_right, y_base = sd['_extents']
    sd['polygons'] = _e4_zones(sd, E4_H)
    # The traffic surcharge is a distributed load over the whole finished top of
    # the wall.  Step 7.5 uses the unfactored SOIL stress for pullout and
    # excludes the live load from it, and XSLOPE's overburden law reads material
    # zones and pore pressure only, so a distributed load never enters sigma'v.
    sd['dloads'] = [[{'X': 0.0, 'Y': E4_H, 'Normal': E4_Q},
                     {'X': x_right, 'Y': E4_H, 'Normal': E4_Q}]]
    sd['dload_dirs'] = ['normal']
    sd['reinforcement_lines'] = [
        _layer('Level %d (%dW)' % (i + 1, ng), E4_H - z, E4_L, ng * tn,
               f_star, _e4_mat_width(ng), E4_WP)
        for i, (z, f_star, _le, ng, tn) in enumerate(E4_LEVELS)]
    _circles(sd, round(E4_L / 2.0, 4), E4_H, E4_H, E4_D)
    return _write(sd, 'fhwa_e4.xlsx')


# ---------------------------------------------------------------------------
# Example E5 -- bridge abutment on a spread footing, ribbed steel strips
# ---------------------------------------------------------------------------
# Steps 1-3 (printed pages E5-4 to E5-5): abutment height Ha = 23 ft measured
# from finished grade to the bottom of the spread footing, embedment d = 2.5 ft,
# so the design height H = 25.5 ft; reinforcement length L = 26 ft at every
# level; footing height h = 10.35 ft standing on top of the wall, base width
# bf = 10.75 ft set back cf = 0.5 ft from the back of the panels.
# Step 2: reinforced fill and the backfill around the footing phi = 34 deg at
# 125 pcf, retained backfill phi = 30 deg at 125 pcf, foundation phi = 30 deg at
# 120 pcf, clayey sand with no water table.
E5_H = 25.5
E5_D = 2.5
E5_L = 26.0
E5_HF = 10.35              # footing height standing on the wall
E5_CF = 0.5                # footing set-back from the back of the panels
E5_BF = 10.75              # footing base width
E5_Y_TOP = E5_H + E5_HF    # 35.85 ft: the roadway surface
#: Step 2: nominal dead and live load reactions from the bridge, k/ft of wall.
E5_DL = 10600.0
E5_LL = 5700.0
#: Step 2: live load on the bridge approach, heqM = 2 ft of soil at 125 pcf.
E5_Q_APPROACH = 250.0
#: Step 1: 10 ft wide precast panels.
E5_WP = 10.0

#: Table E5-8.3 (printed page E5-29), one row per level: depth Z below the
#: bottom of the footing, the pullout resistance factor F*, the resisting length
#: Le, and the horizontal spacing Sh of the strips the design selects.
E5_LEVELS = [
    # Z,     F*,    Le,     Sh
    (1.12,  1.240, 14.75,  1.7),
    (2.35,  1.158, 14.75,  1.7),
    (4.81,  0.995, 14.75,  1.3),
    (7.27,  0.832, 15.06,  1.4),
    (9.73,  0.675, 16.54,  1.4),
    (12.19, 0.675, 18.01,  1.4),
    (14.65, 0.675, 19.49,  1.4),
    (17.11, 0.675, 20.97,  1.4),
    (19.57, 0.675, 22.44,  1.3),
    (22.03, 0.675, 23.92,  1.3),
    (24.49, 0.675, 25.39,  1.3),
]
#: Step 8.4 (printed page E5-27): the same strip section as E3 over a 100-year
#: design life, Tn = 65 ksi (0.154 in2) = 10.00 k/strip.
E5_TN = 10000.0


def build_fhwa_e5():
    sd = _base(-60.0, 85.0, -25.0,
               [_material('Reinforced fill', 125.0, 34.0),
                _material('Retained backfill', 125.0, 30.0),
                _material('Footing backfill', 125.0, 34.0),
                _material('Foundation', 120.0, 30.0)])
    x_left, x_right, y_base = sd['_extents']
    # Step 8.5 states the pullout overburden as sigma'v = gamma (Z + h): the
    # full height above the reinforcement, footing included, read as backfill at
    # 125 pcf, with the footing's own net pressure kept out of the pullout check
    # and carried separately as a spread load.  The model is that statement --
    # the block standing on the wall is backfill, and the bridge reactions are a
    # distributed load, which XSLOPE's overburden law does not read.
    sd['polygons'] = [
        _zone(0, [(0.0, 0.0), (E5_L, 0.0), (E5_L, E5_H), (0.0, E5_H)]),
        _zone(1, [(E5_L, 0.0), (x_right, 0.0), (x_right, E5_H), (E5_L, E5_H)]),
        _zone(2, [(0.0, E5_H), (x_right, E5_H), (x_right, E5_Y_TOP),
                  (0.0, E5_Y_TOP)]),
        _zone(3, [(x_left, E5_D), (0.0, E5_D), (0.0, 0.0), (x_right, 0.0),
                  (x_right, y_base), (x_left, y_base)]),
    ]
    # The bridge reactions over the footing's own base width, and the approach
    # live load beyond it.
    q_bridge = (E5_DL + E5_LL) / E5_BF
    sd['dloads'] = [
        [{'X': E5_CF, 'Y': E5_Y_TOP, 'Normal': q_bridge},
         {'X': E5_CF + E5_BF, 'Y': E5_Y_TOP, 'Normal': q_bridge}],
        [{'X': E5_CF + E5_BF, 'Y': E5_Y_TOP, 'Normal': E5_Q_APPROACH},
         {'X': x_right, 'Y': E5_Y_TOP, 'Normal': E5_Q_APPROACH}],
    ]
    sd['dload_dirs'] = ['normal', 'normal']
    sd['reinforcement_lines'] = [
        _layer('Level %d' % (i + 1), E5_H - z, E5_L, E5_TN, f_star,
               B_STRIP, sh)
        for i, (z, f_star, _le, sh) in enumerate(E5_LEVELS)]
    _circles(sd, round(E5_L / 2.0, 4), E5_Y_TOP, E5_H, E5_D)
    return _write(sd, 'fhwa_e5.xlsx')


# ---------------------------------------------------------------------------
# Example E6 -- the E4 wall under a traffic barrier impact
# ---------------------------------------------------------------------------
# Step 7.5 (printed page E6-6) states the pullout overburden of the impact check
# as sigma'v = gamma (Z + heq): the live load surcharge is inside the vertical
# stress here, where Example E4's static check leaves it out.  The model carries
# it the way heq itself is defined -- as two feet of equivalent soil on the
# reinforced and retained fill alike -- so the wall stands 2 ft taller than the
# E4 model and no distributed load is entered.
#
# Step 7.5 also reads the impact pullout over the FULL length of each layer
# rather than beyond the internal failure surface, so the station these layers
# are read at is the wall face.
E6_HEQ = 2.0
#: F* for the two layers Example E6 checks, printed on page E6-7; the rest of
#: the wall keeps Example E4's ladder, which E6 does not restate.
E6_F_STAR = {0: 1.189, 1: 1.111}


def build_fhwa_e6():
    sd = _e4_base()
    y_top = E4_H + E6_HEQ
    sd['polygons'] = _e4_zones(sd, y_top)
    sd['reinforcement_lines'] = [
        _layer('Level %d (%dW)' % (i + 1, ng), E4_H - z, E4_L, ng * tn,
               E6_F_STAR.get(i, f_star), _e4_mat_width(ng), E4_WP)
        for i, (z, f_star, _le, ng, tn) in enumerate(E4_LEVELS)]
    _circles(sd, round(E4_L / 2.0, 4), y_top, E4_H, E4_D)
    return _write(sd, 'fhwa_e6.xlsx')


# ---------------------------------------------------------------------------
# Example E7 -- the E4 wall under earthquake loading
# ---------------------------------------------------------------------------
# External Step 5 (printed page E7-4): kav = alpha kmax = 1.024 (0.206 g)
# = 0.211 g, the average peak ground acceleration within the reinforced zone.
# Internal Step 4 (printed page E7-14): "For seismic loading conditions, the
# value of F*, the pullout resistance factor, is reduced to 80 percent of the
# value used for static design", so every layer's Delta is taken on 0.8 F*.
# Everything else is the E4 wall, live load surcharge included, because the
# static pullout resistances E7 reduces are E4's own.
E7_KAV = 0.211
E7_F_STAR_FACTOR = 0.8


def build_fhwa_e7():
    sd = _e4_base()
    x_left, x_right, y_base = sd['_extents']
    sd['k_seismic'] = E7_KAV
    sd['polygons'] = _e4_zones(sd, E4_H)
    sd['dloads'] = [[{'X': 0.0, 'Y': E4_H, 'Normal': E4_Q},
                     {'X': x_right, 'Y': E4_H, 'Normal': E4_Q}]]
    sd['dload_dirs'] = ['normal']
    sd['reinforcement_lines'] = [
        _layer('Level %d (%dW)' % (i + 1, ng), E4_H - z, E4_L, ng * tn,
               E7_F_STAR_FACTOR * f_star, _e4_mat_width(ng), E4_WP)
        for i, (z, f_star, _le, ng, tn) in enumerate(E4_LEVELS)]
    _circles(sd, round(E4_L / 2.0, 4), E4_H, E4_H, E4_D)
    return _write(sd, 'fhwa_e7.xlsx')


BUILDERS = [build_fhwa_e3, build_fhwa_e4, build_fhwa_e5, build_fhwa_e6,
            build_fhwa_e7]


if __name__ == '__main__':
    for fn in BUILDERS:
        print('built', fn())
