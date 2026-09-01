"""Builder for GeoStudio SLOPE/W Verification Manual problem 2.45 —
Eurocode 7: Cutting in Clay.

Source: GeoStudio Slope Stability Verification Manual (Oct 2022) §2.45
(printed pp. 90-91), after the example on p.202 of "Designers' Guide to
Eurocode 7: Geotechnical Design". A 1:2 cutting in a homogeneous clay with a
water table and a permanent 35 kPa surcharge on the crest, analysed with
Eurocode 7 Design Approach 3 (DA3).

  Geometry (Fig 133): ground surface (0,5)-(30,5)-(50,15)-(80,15); model base
    at y=0 (domain 0..80 wide). Units metric (kPa, kN/m3, m).
  Material (Table 115, characteristic): Clay c'=10 kPa, phi'=28 deg,
    gamma=20 kN/m3.
  Water table: the book gave no coordinates, so SLOPE/W used an approximate
    line. We take that line from the .gsz: (0,5)-(30,5)-(45.5,9)-(70,12)-(82,12).
  Surcharge: 35 kPa permanent, on the crest strip (52,15)-(62,15) — the load's
    footprint read from the .gsz surcharge datapoints (52,16)-(62,16), matching
    the hatched symbol in Fig 133 (it covers ~10 m of crest, not the whole top).

FACTORED-PARAMETER EMULATION
  xslope has no native EC7 partial-factor feature, so we bake the DA3 (set M2)
  partial factors into the authored material and load, then run an ordinary
  analysis. The resulting factor of safety of the *factored* model equals
  GeoStudio's "Overdesign Factor" (ODF) — the ratio it reports for an EC7 run.
  Partial factors read from the .gsz DesignFactors block ("Eurocode 7 - DA3"),
  which matches EN 1997 DA3 set M2:
      EffectiveCohesion  gamma_c'    = 1.25   -> c*   = 10 / 1.25      = 8.0 kPa
      EffectiveTanPhi    gamma_phi'  = 1.25   -> phi* = atan(tan28/1.25) = 23.043 deg
      SoilWeight         gamma_gamma = 1.0    -> gamma unchanged (20)
      PermanentLoads     gamma_G     = 1.0    -> surcharge unchanged (35 kPa)
  (In this model the crest surcharge is a *permanent* action factored by 1.0, so
  the emulation is a pure material-strength reduction — no action scaling needed.
  Variable loads would carry 1.3, but there are none here.)

Targets (Table, printed p.91): ODF Spencer 1.174 / Bishop 1.173 / Janbu 1.033;
book (Bishop) 1.193. The .gsz solves Spencer ODF = 1.1745 at the critical
circle center (35.4, 25.467), R = 22.099.
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx as _write_xlsx  # noqa: E402
from _gs2_donor import donor_material, load_donor  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- DA3 (set M2) partial factors, from the .gsz DesignFactors block ---
GAMMA_C = 1.25
GAMMA_PHI = 1.25

C_CHAR, PHI_CHAR, GAMMA = 10.0, 28.0, 20.0
C_STAR = C_CHAR / GAMMA_C
PHI_STAR = math.degrees(math.atan(math.tan(math.radians(PHI_CHAR)) / GAMMA_PHI))


def build():
    sd = load_donor(ACADS_1A)   # valid metric single-material template

    m = donor_material(sd)
    m.update(name='Clay (DA3 factored)', c=round(C_STAR, 4), phi=round(PHI_STAR, 4),
             gamma=GAMMA, gamma_sat=GAMMA, option='mc', u='piezo')
    sd['materials'] = [m]

    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 5.0), (30.0, 5.0), (50.0, 15.0), (80.0, 15.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 9.807

    # water table (approximate line SLOPE/W used; book gave no coords)
    sd['piezo_line'] = [(0.0, 5.0), (30.0, 5.0), (45.5, 9.0), (70.0, 12.0), (82.0, 12.0)]

    # 35 kPa permanent surcharge on the crest strip x=52..62 (footprint from the
    # .gsz surcharge datapoints; factored by 1.0 under DA3)
    sd['dloads'] = [[
        {'X': 52.0, 'Y': 15.0, 'Normal': 35.0},
        {'X': 62.0, 'Y': 15.0, 'Normal': 35.0},
    ]]

    # circular search seeded at SLOPE/W's own critical center (35.4, 25.467),
    # R = 22.099 -> tangent elevation Depth = Yo - R = 3.368
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [
        {'Xo': 35.4, 'Yo': 25.467, 'Depth': 3.368, 'R': 22.099},
        {'Xo': 33.0, 'Yo': 22.0, 'Depth': 3.0, 'R': 19.0},
        {'Xo': 38.0, 'Yo': 28.0, 'Depth': 3.5, 'R': 24.5},
    ]

    path = os.path.join(OUT, 'gs2_45.xlsx')
    _write_xlsx(sd, path)
    return path


if __name__ == '__main__':
    p = build()
    print('wrote', p)
    print('factored: c*=%.4f kPa, phi*=%.4f deg, gamma=%.1f' % (C_STAR, PHI_STAR, GAMMA))
