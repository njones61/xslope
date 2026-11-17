"""Build the two models of Duncan, Wright & Brandon's Table 10.4.

Duncan, J.M., Wright, S.G., and Brandon, T.L. (2014), *Soil Strength and Slope
Stability*, 2nd ed., Wiley, Section 10.4.2, Figure 10.7 and Tables 10.3 and 10.4.

The table compares the simplified single-stage pseudostatic analysis (the R
envelope) with the two-stage analysis (consolidation stresses from a static first
stage, undrained strengths in the second) on one slope, with k = 0.15, for four
soils and two water conditions: dry, and submerged before the earthquake.

Geometry (Figure 10.7): a 2.5H:1V slope 30 ft high, gamma = 135 pcf, one soil
throughout. x = 0 at the toe, y = 0 at the toe's elevation. The figure gives no
extents and no foundation depth, so:

  - Level ground runs 90 ft (3H) in front of the toe and 90 ft behind the crest,
    far enough that no searched circle reaches either end.
  - The base of the model is 30 ft (one slope height) below the toe. The figure
    draws a foundation and the text gives its depth nowhere. On the dry slope the
    critical circles stay well above any base that deep, and the answers do not
    depend on it. On the submerged slope the stronger soils' circles run deep, and
    a base 60 ft down gave values 2 to 3% below the table's, where one slope height
    brings them within about 0.5%.

Submerged (Case II): the water stands 15 ft above the crest (a piezometric line at
y = 45 across the model), with Water loads = auto, so the ponded water is the
derived distributed load. The 15 ft is an assumption; the figure shows the slope
under water without dimensioning the depth.

Materials (Table 10.3): each soil appears twice, so one workbook per water case
carries all eight strength sets, and a verification tag picks one with
``material=<name>`` (run_tests makes every row a copy of the named one):

  - "Soil N two-stage": c', phi', and the K_c = 1 envelope's d and psi, for the
    two-stage analysis.
  - "Soil N R envelope": c_R and phi_R entered as c and phi, with the pore
    pressures before the earthquake (the book's footnote 1), for the single-stage
    analysis.

The profile line points at the first row; Studio users pick a soil by pointing
it at another.
"""
import os
import sys

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'docs', 'verification', 'files', 'published')
SEED = os.path.join(ROOT, 'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')
NAMES = {'dry': 'dwb_t10_4_dry.xlsx', 'submerged': 'dwb_t10_4_submerged.xlsx'}

H = 30.0                  # slope height, ft
RUN = 2.5                 # 2.5H:1V
GAMMA = 135.0             # pcf, every soil
K = 0.15                  # seismic coefficient
X_LEFT, X_RIGHT = -3 * H, RUN * H + 3 * H
Y_BASE = -H               # one slope height below the toe
Y_POOL = H + 15.0         # submerged case: water 15 ft above the crest

#: Table 10.3: c' (psf), phi' (deg), c_R (psf), phi_R (deg), d (psf), psi (deg).
SOILS = {
    1: (0.0, 45.0, 60.0, 23.0, 64.0, 24.4),      # Pilarcitos Dam fill, 0-10 psi
    2: (200.0, 31.0, 700.0, 15.0, 782.0, 16.7),  # Rio Blanco sandy clay
    3: (0.0, 34.0, 300.0, 15.5, 327.0, 16.8),    # Pilarcitos fill, 0-100 psi
    4: (0.0, 35.0, 1400.0, 22.5, 1716.0, 26.9),  # Hirfanli Dam fill
}


def _material(name, c, phi, d, psi, u):
    """One Mohr-Coulomb soil, every column stated; the table publishes strengths
    and one unit weight, so the rest is zero or empty."""
    return {
        'name': name, 'gamma': GAMMA, 'gamma_sat': None, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': d, 'psi': psi,
        't_cut': None, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': u, 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': 0.0, 'k2': 0.0, 'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0, 'h0': 0.0, 'vg_a': 0.0, 'vg_n': 0.0,
        'E': None, 'nu': None,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _slope_data(case):
    u = 'none' if case == 'dry' else 'piezo'
    sd = load_slope_data(SEED)
    sd['unit_system'] = 'imperial'
    sd['gamma_water'] = 62.4
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = K
    sd['k0'] = None
    mats = []
    for n, (c, phi, c_r, phi_r, d, psi) in SOILS.items():
        mats.append(_material(f'Soil {n} two-stage', c, phi, d, psi, u))
        mats.append(_material(f'Soil {n} R envelope', c_r, phi_r, 0.0, 0.0, u))
    sd['materials'] = mats
    sd['profile_lines'] = [{'coords': [(X_LEFT, 0.0), (0.0, 0.0), (RUN * H, H), (X_RIGHT, H)],
                            'mat_id': 0, 'size': None}]
    sd['polygons'] = []
    sd['max_depth'] = Y_BASE
    sd['piezo_line'] = [] if case == 'dry' else [(X_LEFT, Y_POOL), (X_RIGHT, Y_POOL)]
    sd['piezo_line2'] = []
    sd['water_loads'] = 'auto'
    for key in ('dloads', 'dload_dirs', 'dloads2', 'dload2_dirs', 'reinforcement_lines',
                'reinforce_lines', 'pile_lines', 'line_loads', 'non_circ'):
        sd[key] = []
    sd['circular'] = True
    sd['element_type'] = None
    sd['target_size'] = None
    sd['search_window'] = {}
    # One starting circle, centered above the middle of the face and tangent to the
    # toe's elevation; the search moves it.
    sd['circles'] = [{'Xo': RUN * H / 2.0, 'Yo': 70.0, 'Depth': 0.0, 'R': 70.0}]
    return sd


def build_dwb_t10_4():
    os.makedirs(OUT, exist_ok=True)
    out = []
    for case, name in NAMES.items():
        dest = os.path.join(OUT, name)
        save_slope_data_to_xlsx(_slope_data(case), dest)
        out.append(dest)
    return out


if __name__ == '__main__':
    for p in build_dwb_t10_4():
        print('built', p)
