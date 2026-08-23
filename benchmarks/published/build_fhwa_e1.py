"""Build the FHWA Example E1 MSE wall model.

Berg, R.R., Christopher, B.R., and Samtani, N.C. (2009), "Design and
Construction of Mechanically Stabilized Earth Walls and Reinforced Soil
Slopes", FHWA-NHI-10-025 (Volume II), Appendix E, Example E1.

Example E1 is a 20 ft design-height MSE wall faced with modular block units and
reinforced with eleven geogrid layers, carrying a 2H:1V broken backslope 9 ft
high and a 250 psf traffic surcharge on the level retained backfill behind the
crest.  Step 7.8 of the example checks every layer for pullout, and the model
built here is that check expressed as an XSLOPE input file.

Geometry (all of it derived from the text of Steps 1-4, none scaled off a
figure).  x = 0 at the wall face, y = 0 at the top of the leveling pad:

  - Exposed height He = 18 ft, embedment 2 ft, so the design height H = 20 ft
    and the finished grade in front of the wall stands at y = 2.
  - The 3 degree facing batter is treated as vertical, which is what the
    example itself instructs in Step 4 ("assume a vertical face").
  - The backslope rises 2H:1V from the top of the wall at (0, 20).  Step 4
    places the slope crest at the back end of the reinforcement, which is what
    fixes the total height h = 20 + 9 = 29 ft; the crest is therefore at
    (18, 29) and the retained backfill is level beyond it.
  - Reinforcement length L = 18 ft (0.9H) at every level, and the eleven layers
    sit at the depths Z of Table E1-7.5 below the top of the wall.

Soils are the three of Step 2, all at 125 pcf and drained: reinforced wall fill
phi = 34 degrees, retained backfill phi = 30 degrees, foundation phi = 30
degrees.  There is no groundwater.  Because all three unit weights are equal,
the vertical stress the pullout law reads along a reinforcement layer is
unaffected by where the reinforced/retained zone boundary is drawn.

The pullout law.  FHWA states pullout resistance as

    Pr = F* alpha sigma'v Le C Rc

with F* = 0.45 and alpha = 0.8 for these geogrids, C = 2 (two bearing faces),
and Rc = 1.0 (continuous geogrid, full coverage).  XSLOPE's overburden-dependent
law states the same resistance per unit length of line as

    r(s) = 2 (Adhesion + sigma'v(s) tan(Delta))

so the two are the same statement with Adhesion = 0 and

    Delta = arctan(F* alpha) = arctan(0.45 x 0.8) = arctan(0.36) = 19.80 degrees

and XSLOPE's factor of two carrying FHWA's C.  Rc = 1 needs no representation:
the law is already per unit width of a continuous sheet.

Tend1 = Tmax at the face end.  FHWA's Step 7.8 counts only the resisting-zone
side of each layer; the front side is carried by the facing connection and is
checked separately in Step 7.9.  Setting the face-end anchorage to the layer's
own tensile strength makes the face branch of the envelope non-governing, so
what the resisting side develops is what the envelope reports.  XSLOPE does not
model the geogrid-to-block connection, so no connection capacity is entered.

Tmax is the nominal long-term strength Tal of Table E1-7.3: 1,085 lb/ft for the
GG-I layers (1-4) and 2,169 lb/ft for the GG-II layers (5-11).  These are the
example's own values, and they are far below the pullout resistance at every
level, which is the example's finding: pullout capacity-demand ratios run 11 to
31 while rupture runs 1.00 to 2.82.

The traffic surcharge is entered as a distributed load of 250 psf over the level
retained backfill from the crest at x = 18 out to the right edge of the domain.
AASHTO excludes live load from the vertical stress used for pullout, and this
placement cannot reach the reinforcement two ways over: the surcharge begins at
x = 18, where the reinforcement ends, and XSLOPE's overburden law reads material
zones and pore pressure only, so a distributed load never enters sigma'v.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/published/build_fhwa_e1.py
"""
import math
import os

from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'docs', 'verification', 'files', 'published')
SEED = os.path.join(ROOT, 'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')
DEST = os.path.join(OUT, 'fhwa_e1.xlsx')

# --- Step 1-3: geometry -----------------------------------------------------
H = 20.0                 # design height, ft (18 ft exposed + 2 ft embedment)
D_EMBED = 2.0            # embedment below finished grade in front of the wall
L_REINF = 18.0           # reinforcement length, ft (0.9H)
H_SLOPE = 9.0            # backslope height, ft
SLOPE_RUN = 2.0          # 2H:1V broken backslope
Y_CREST = H + H_SLOPE    # 29 ft: the crest sits at the back of the reinforcement
# Step 4 states three things about the backslope that have to agree: it rises
# 2H:1V, it is 9 ft high, and its crest is at the back end of the reinforcement.
# Any two of them fix the third, so the agreement is asserted rather than assumed.
assert abs(L_REINF / SLOPE_RUN - H_SLOPE) < 1e-9

#: Model extents.  The example dimensions the wall and nothing around it, so the
#: domain is opened until the searched mechanism sits clear of every edge: 35 ft
#: of level ground in front of the toe, 42 ft of retained backfill behind the
#: crest, and 15 ft of foundation below the leveling pad.
X_LEFT, X_RIGHT, Y_BASE = -35.0, 60.0, -15.0

# --- Step 2: soils ----------------------------------------------------------
GAMMA = 125.0
PHI_R, PHI_B, PHI_F = 34.0, 30.0, 30.0

# --- Step 7: pullout law ----------------------------------------------------
F_STAR = 0.45            # pullout resistance factor for these geogrids
ALPHA = 0.8              # scale correction factor for these geogrids
#: arctan(F* alpha) = arctan(0.36), to the two decimals the reinforce sheet takes
DELTA = round(math.degrees(math.atan(F_STAR * ALPHA)), 2)     # 19.80 degrees

# --- Step 7.4 / 7.6: layout and nominal long-term strength ------------------
#: Depth Z of each layer below the top of the wall, ft (Table E1-7.5).  The top
#: layer sits 8 in below the top of the wall and the bottom layer 8 in above the
#: leveling pad, with the block courses in between.
LAYER_Z = [0.67, 2.67, 4.67, 6.67, 8.67, 10.67, 12.67, 14.67, 16.67, 18.67,
           19.33]
#: Grade at each level (Table E1-7.4) and its nominal long-term strength Tal
#: (Table E1-7.3), lb/ft.
TAL = {'GG-I': 1085.0, 'GG-II': 2169.0}
LAYER_GRADE = ['GG-I'] * 4 + ['GG-II'] * 7

#: The traffic surcharge of Step 4: 2.0 ft of equivalent soil height at 125 pcf.
Q_TRAFFIC = 250.0


def _material(name, phi):
    """One drained Mohr-Coulomb soil, every column stated.

    The example publishes a unit weight and a friction angle and nothing else,
    so the cohesion is zero, the uncertainties are zero, and the stiffness,
    tensile cutoff and unsaturated columns are left empty rather than carrying
    invented values.
    """
    return {
        'name': name, 'gamma': GAMMA, 'gamma_sat': None, 'option': 'mc',
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


def _layer(i):
    """Reinforcement layer i (0-based), horizontal, face to back of zone."""
    y = H - LAYER_Z[i]
    t_max = TAL[LAYER_GRADE[i]]
    return {
        'label': 'Layer %d (%s)' % (i + 1, LAYER_GRADE[i]),
        'x1': 0.0, 'y1': round(y, 4),
        'x2': L_REINF, 'y2': round(y, 4),
        'type': 'geosynthetic', 'dir': 'tangent', 'appl': 'active',
        't_max': t_max, 't_res': None,
        # Lp1/Lp2 are not read once Adhesion and Delta are filled: the
        # overburden law carries the whole development, so a development length
        # would be a second, contradictory statement of the same bond.
        'lp1': 0.0, 'lp2': 0.0,
        'adhesion': 0.0, 'delta': DELTA,
        # The face end is anchored to the layer's own strength so the face
        # branch of the envelope never governs; the connection itself is
        # checked in Step 7.9 of the example and is not modeled here.
        'tend1': t_max, 'tend2': 0.0,
        'spacing': 1.0,          # continuous sheet, already per unit width
        'E': 0.0, 'area': 0.0,   # no FEM run on this model
    }


def _slope_data():
    sd = load_slope_data(SEED)
    sd['unit_system'] = 'imperial'
    sd['gamma_water'] = 62.4
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = 0.0
    sd['materials'] = [_material('Reinforced fill', PHI_R),
                       _material('Retained backfill', PHI_B),
                       _material('Foundation', PHI_F)]
    # Polygon geometry: the vertical wall face and the two vertical zone
    # boundaries are edges of the zones themselves, which profile lines (one
    # elevation per station) cannot express.
    sd['polygons'] = [
        _zone(0, [(0.0, 0.0), (L_REINF, 0.0), (L_REINF, Y_CREST), (0.0, H)]),
        _zone(1, [(L_REINF, 0.0), (X_RIGHT, 0.0), (X_RIGHT, Y_CREST),
                  (L_REINF, Y_CREST)]),
        _zone(2, [(X_LEFT, D_EMBED), (0.0, D_EMBED), (0.0, 0.0),
                  (X_RIGHT, 0.0), (X_RIGHT, Y_BASE), (X_LEFT, Y_BASE)]),
    ]
    sd['profile_lines'] = []
    sd['max_depth'] = Y_BASE
    sd['piezo_line'] = []
    sd['piezo_line2'] = []
    # Traffic surcharge on the level retained backfill behind the crest.
    sd['dloads'] = [[{'X': L_REINF, 'Y': Y_CREST, 'Normal': Q_TRAFFIC},
                     {'X': X_RIGHT, 'Y': Y_CREST, 'Normal': Q_TRAFFIC}]]
    sd['dload_dirs'] = ['normal']
    sd['dloads2'] = []
    sd['dload2_dirs'] = []
    sd['reinforcement_lines'] = [_layer(i) for i in range(len(LAYER_Z))]
    sd['reinforce_lines'] = []
    sd['pile_lines'] = []
    sd['line_loads'] = []
    sd['non_circ'] = []
    sd['circular'] = True
    sd['water_loads'] = 'manual'
    sd['element_type'] = None
    sd['target_size'] = None
    sd['search_window'] = {}
    # Starting circles: center over the middle of the reinforced zone, one
    # tangent to the finished grade in front of the wall and one tangent to the
    # top of the leveling pad, so the search opens with a surface through the
    # reinforced zone and one that undercuts it.  The center is placed at the
    # crest elevation plus a fifth of the height rather than at the toe plus 2H
    # that a slope would use: the face here is vertical, and on a vertical face
    # the arc through the toe that a high center draws is far flatter than any
    # mechanism the wall has.  The search converges on the same surface from
    # either seed.
    x_o = round(L_REINF / 2.0, 4)
    y_o = Y_CREST + H / 4.0
    sd['circles'] = [
        {'Xo': x_o, 'Yo': y_o, 'Depth': D_EMBED, 'R': y_o - D_EMBED},
        {'Xo': x_o, 'Yo': y_o, 'Depth': 0.0, 'R': y_o},
    ]
    return sd


def build_fhwa_e1():
    sd = _slope_data()
    os.makedirs(OUT, exist_ok=True)
    save_slope_data_to_xlsx(sd, DEST)
    return DEST


if __name__ == '__main__':
    print('built', build_fhwa_e1())
