"""Finite-element variants of Slide #106 / Cai & Ugai (2000) -- the 2D-vs-3D
pile diagnostic on docs/verification/rocscience.md (VP106).

Cai & Ugai analysed this slope BOTH ways: with a three-dimensional
shear-strength-reduction finite-element model that meshes the individual piles
and the soil between them (with interface elements that can slip), and with
limit equilibrium using their own Ito & Matsui limit pressure. The
limit-equilibrium side is already verified -- vp106a..e reproduce Slide2 and the
paper across four spacings.

These three files drive the OTHER side of the paper through XSLOPE's SSRM, at the
one spacing where the published 3D answer spreads widest from limit equilibrium
(D1/D = 3):

  vp106a_fem.xlsx        no pile                       Cai & Ugai 3D FE  1.14
  vp106c_fem.xlsx        pile, free head               Cai & Ugai 3D FE  1.36
  vp106c_fem_fix.xlsx    pile, head rotation restrained Cai & Ugai 3D FE 1.45

The comparison is a DIAGNOSTIC, not a verification: their model is
three-dimensional, so the soil arches between piles and slips against them, while
XSLOPE's pile is a plane-strain beam whose stiffness is smeared over the spacing.
The deltas are the point -- they measure what the 2D idealization costs on a
discrete pile row. The docs entry says so and the row carries no match dot.

Geometry and soil are vp106's, unchanged (see _vp106_slope_data in
build_problems.py): 10 m toe flat, 1V:1.5H slope 10 m high, 10 m crest flat, 10 m
of soil over bedrock; gamma = 20 kN/m3, c' = 10 kPa, phi' = 20 deg, dry.

What the FE variants add:
  * soil stiffness E = 2e5 kPa, nu = 0.25 (the paper's elastic properties; the
    limit-equilibrium files carry the template's donor values, which no LEM
    method reads).
  * pile stiffness E = 6e7 kPa on the 0.8 m section. I and A are left blank so
    build_fem_data derives them from the diameter -- pi D^4/64 = 0.0201062 m4 and
    pi D^2/4 = 0.5026548 m2, the paper's own section constants. The spacing
    S = 3D = 2.4 m divides both, which is the plane-strain smear the diagnostic
    is about.
  * the third file sets the pile head fixity to 'fixed', which restrains the head
    ROTATION only -- the paper's "unrotated" pile.

The Ito & Matsui force (the LEM path's H) is NOT set on these files: in a finite
element analysis the pile carries what the soil pushes onto it, and setting a
limit pressure as well would count the pile twice.

Run:  PYTHONPATH=. python3 benchmarks/rocscience/build_vp106_fem.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
OUT = os.path.abspath(os.path.join(REPO, 'docs', 'verification', 'files',
                                   'rocscience'))
DONOR = os.path.abspath(os.path.join(REPO, 'docs', 'lem', 'files',
                                     'xslope_acads_simple.xlsx'))

E_SOIL = 2.0e5      # kPa
NU_SOIL = 0.25
E_PILE = 6.0e7      # kPa
D_PILE = 0.8        # m
SD_RATIO = 3        # D1/D -- the spacing the diagnostic is stated at


def _slope_data(with_pile, head_fixity='free'):
    sd = load_slope_data(DONOR)
    # The donor is RS2-1 and declares an isotropic at-rest initial stress
    # (main!D16, K0 = 1) that RS2 solved its SSR under. It belongs to that
    # problem alone, so it is cleared here rather than riding into files whose
    # own analysis is authored without one (benchmarks/tag_k0.py).
    sd['k0'] = None
    m = sd['materials'][0]
    m.update(name='Soil', c=10.0, phi=20.0, gamma=20.0, gamma_sat=20.0,
             option='mc', u='none', ru=0.0, E=E_SOIL, nu=NU_SOIL, t_cut=None)
    sd['materials'] = [m]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (10.0, 0.0), (25.0, 10.0),
                                 (35.0, 10.0)], 'size': None},
    ]
    sd['max_depth'] = -10.0
    sd['circles'] = [{'Xo': 17.0, 'Yo': 16.0, 'Depth': -2.0, 'R': 18.0}]
    sd['non_circ'] = []
    sd['circular'] = True
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['pile_lines'] = []
    if with_pile:
        sd['pile_lines'] = [{
            'x1': 17.5, 'y1': 5.0, 'x2': 17.5, 'y2': -10.0,
            # H blank: the FE pile carries what the soil pushes onto it. I and
            # area blank: derived from D as pi D^4/64 and pi D^2/4.
            'H': None, 'theta_p': 0.0, 'D_pile': D_PILE, 'S': SD_RATIO * D_PILE,
            'E': E_PILE, 'I': None, 'area': None,
            'V_cap': None, 'M_cap': None, 'head_fixity': head_fixity,
            'appl': 'passive', 'label': f'pile row D1={SD_RATIO}D ({head_fixity} head)',
        }]
    return sd


def vp106a_fem():
    """No-pile baseline. Cai & Ugai 3D FE: 1.14."""
    save_slope_data_to_xlsx(_slope_data(False),
                            os.path.join(OUT, 'vp106a_fem.xlsx'))
    return 'vp106a_fem.xlsx'


def vp106c_fem():
    """D1/D = 3, free pile head. Cai & Ugai 3D FE: 1.36."""
    save_slope_data_to_xlsx(_slope_data(True, 'free'),
                            os.path.join(OUT, 'vp106c_fem.xlsx'))
    return 'vp106c_fem.xlsx'


def vp106c_fem_fix():
    """D1/D = 3, head rotation restrained. Cai & Ugai 3D FE: 1.45."""
    save_slope_data_to_xlsx(_slope_data(True, 'fixed'),
                            os.path.join(OUT, 'vp106c_fem_fix.xlsx'))
    return 'vp106c_fem_fix.xlsx'


BUILDERS = [vp106a_fem, vp106c_fem, vp106c_fem_fix]


if __name__ == '__main__':
    for b in BUILDERS:
        print('wrote', b())
