"""Builders for RS2-corpus problems that have NO Slide2 counterpart file
(docs/verification/rs2.md). Most RS2 problems run on the Slide2 corpus files;
these are the exceptions - currently the Pruska (2003) multi-program
cross-bearing slopes (RS2 #56-58).

Run from the repo root:  python benchmarks/rocscience/build_rs2.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

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


if __name__ == '__main__':
    for fn in (rs2_56a, rs2_56b, rs2_57a, rs2_57b, rs2_58a, rs2_58b, hammah_hb1,
               rs2_60a, rs2_60b, rs2_60c):
        print(fn())
