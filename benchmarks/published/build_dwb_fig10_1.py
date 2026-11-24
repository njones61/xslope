"""Build the closed-form yield coefficient check: a rigid wedge on a planar surface.

Duncan, J.M., Wright, S.G., and Brandon, T.L. (2014), *Soil Strength and Slope
Stability*, 2nd ed., Wiley, Section 10.1.2, Figure 10.1.

Figure 10.1 resolves the forces on a block sliding on a plane inclined at beta,
with the seismic force kW through its center of gravity:

    N = W cos(beta) - kW sin(beta)          (perpendicular to the plane)
    T = W sin(beta) + kW cos(beta)          (along the plane)
    F = (cL + N tan(phi)) / T

The book applies it to an infinite slope. The same two equations hold exactly for
any mass sliding on one plane, because every slice base shares the plane's
direction and the interslice forces cancel in the sum: every method that
satisfies force equilibrium (Spencer, Morgenstern-Price, Corps, Lowe and
Karafiath) returns this F on a planar surface, at any k, with any slicing. Setting
F = 1 and solving for k gives the yield coefficient in closed form:

    k_y = (cL/W + cos(beta) tan(phi) - sin(beta)) / (cos(beta) + sin(beta) tan(phi))

which is tan(phi - beta) when c = 0, the yield acceleration of Newmark's (1965)
block on a plane.

The model: a 1H:1V slope 10 m high, gamma = 20 kN/m3, dry. The surface is one
plane at beta = 30 degrees from the toe (0, 0) to the crest at (10/tan 30, 10),
so the wedge is the triangle (0, 0), (10, 10), (17.3205, 10):

    W = 20 x 36.6025 = 732.05 kN/m,  L = 20 m.

Two materials, each picked by a tag's ``material=<name>``:

  - "c-phi": c = 10 kPa, phi = 30 deg. k_y = 0.273205 / 1.154701 = 0.23660.
  - "frictional": c = 0, phi = 35 deg. k_y = tan 5 deg = 0.08749.
"""
import math
import os

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'docs', 'verification', 'files', 'published')
SEED = os.path.join(ROOT, 'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')
NAME = 'dwb_fig10_1_wedge.xlsx'

H = 10.0                    # slope height, m
GAMMA = 20.0                # kN/m3
BETA = 30.0                 # the plane, deg
X_CREST = H                 # 1H:1V face from (0, 0) to (10, 10)
X_PLANE = H / math.tan(math.radians(BETA))   # where the plane meets the crest

#: name: (c in kPa, phi in deg)
SOILS = {'c-phi': (10.0, 30.0), 'frictional': (0.0, 35.0)}


def closed_form_ky(c, phi):
    """k_y of the wedge from Figure 10.1's two equations with F = 1."""
    area = 0.5 * (X_PLANE - X_CREST) * H
    w = GAMMA * area
    length = math.hypot(X_PLANE, H)
    b, p = math.radians(BETA), math.radians(phi)
    return ((c * length / w + math.cos(b) * math.tan(p) - math.sin(b))
            / (math.cos(b) + math.sin(b) * math.tan(p)))


def _material(name, c, phi):
    return {
        'name': name, 'gamma': GAMMA, 'gamma_sat': None, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,
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


def _slope_data():
    sd = load_slope_data(SEED)
    sd['unit_system'] = 'metric'
    sd['gamma_water'] = 9.81
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = 0.0
    sd['k0'] = None
    sd['materials'] = [_material(n, c, phi) for n, (c, phi) in SOILS.items()]
    sd['profile_lines'] = [{'coords': [(-2 * H, 0.0), (0.0, 0.0), (X_CREST, H),
                                       (4 * H, H)],
                            'mat_id': 0, 'size': None}]
    sd['polygons'] = []
    sd['max_depth'] = -H
    sd['piezo_line'] = []
    sd['piezo_line2'] = []
    sd['water_loads'] = 'auto'
    for key in ('dloads', 'dload_dirs', 'dloads2', 'dload2_dirs', 'reinforcement_lines',
                'reinforce_lines', 'pile_lines', 'line_loads', 'circles'):
        sd[key] = []
    sd['circular'] = False
    sd['non_circ'] = [{'X': 0.0, 'Y': 0.0, 'Movement': 'Fixed'},
                      {'X': X_PLANE, 'Y': H, 'Movement': 'Fixed'}]
    sd['element_type'] = None
    sd['target_size'] = None
    sd['search_window'] = {}
    return sd


def build_dwb_fig10_1():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, NAME)
    save_slope_data_to_xlsx(_slope_data(), dest)
    return dest


if __name__ == '__main__':
    print('built', build_dwb_fig10_1())
    for n, (c, phi) in SOILS.items():
        print(f'  {n}: k_y = {closed_form_ky(c, phi):.5f}')
