"""Builder for GeoStudio SLOPE/W Verification Manual sec. 2.26 -- Baker - Planar
Homogeneous (docs/verification/geostudio.md, "2.26 -- Baker - Planar
Homogeneous").

Baker (2001) evaluates the factor of safety of every PLANAR slip surface
through the toe of a homogeneous dry slope, as a function of where the plane
daylights on the crest/top-of-slope ground surface (Baker's X = x/H). The
manual (p.58-60, Figure 89/Table 73/Figure 91-92) reproduces this with
SLOPE/W's Morgenstern-Price solver: a single Mohr-Coulomb material
(c'=30 kPa, phi'=30 deg, gamma=20 kN/m3, dry), and reports the critical plane
at X = x/H = 0.85 (FS = 1.352), matching Baker & Leshchinsky (2001)'s own
curve (min. FS approx 1.35 at the same X).

Geometry is pinned to the vendor .gsz (dev-side cross-check only, never
committed/redistributed -- see xslope.geostudio.read_gsz), which stores the
model's four boundary points exactly:
    1 = (0, 0)     toe (origin of the trial-plane fan)
    2 = (20, 0)    bottom-right corner of the analysis domain
    3 = (20, 10)   top-right corner (flat ground extends well past the crest)
    4 = (2.5, 10)  crest (top of the H=10 m slope face)
The ground surface is therefore (0,0) -> (2.5,10) [slope face] -> (20,10)
[flat backslope], with a flat base at y=0 (max_depth=0.0) -- exactly the gsz's
own rectangle/toe-triangle boundary. SLOPE/W's own trial-plane fan in the gsz
(the 'SlipSurface' data points under the stability analysis) daylights at
x = 6.0 .. 20.0 in 0.5 m steps -- i.e. X = x/H = 0.6 .. 2.0 -- and its saved
per-surface FS confirms the critical plane at x = 8.5 (X = 0.85), FS =
1.352158, weight = 600 kN/m (a 30 m^2 wedge x 20 kN/m3), matching the printed
Figure 91(c)/92 values to the last digit.

This is the same problem as Rocscience corpus #43 (docs/verification/
rocscience.md, vp043.xlsx -- also RocPlane), which is 'partial' there because
the manual's printed c'/phi'/gamma table did not reproduce Slide/RocPlane/
Baker's approx.1.35 with xslope's own Culmann/Janbu hand check (Slide's face
was 49.5 deg from a *different*, unlabeled H/crest-offset geometry). Reading
the exact geometry out of the GeoStudio .gsz (H=10, crest offset=2.5 m, so
face angle = atan(10/2.5) = 75.96 deg -- steeper than vp043's inferred 73.3
deg) is what resolves it: with this geometry, SLOPE/W's own M-P solve
reproduces Baker's curve, and xslope's M-P/Spencer on the identical planar
surface should as well. See docs/verification/geostudio.md sec. 2.26 for the
result.

Run: PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_26.py
"""

import math
import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

H = 10.0            # slope height (gsz points 1 -> 4)
CREST_X = 2.5        # crest x-offset from the toe (gsz point 4)
DOMAIN_X = 20.0       # right edge of the analysis domain (gsz points 2, 3)
CRIT_X_OVER_H = 0.85  # Baker (2001) / manual Fig. 91(c)/92 critical plane


def gs2_26():
    """Baker (2001) planar-homogeneous slope, evaluated on the manual's
    critical plane (X = x/H = 0.85 through the toe). c'=30, phi'=30, gamma=20,
    dry, Mohr-Coulomb (manual Table 73)."""
    sd = load_slope_data(ACADS_1A)
    m = sd['materials'][0]
    m.update(name='MC Material', c=30.0, phi=30.0, gamma=20.0, option='mc', u='none')
    sd['materials'] = [m]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (CREST_X, H), (DOMAIN_X, H)]},
    ]
    sd['max_depth'] = 0.0
    sd['circular'] = False
    sd['circles'] = []
    # Critical planar surface: toe (0,0) to the daylight point at x/H=0.85 on
    # the flat backslope (8.5, 10) -- both points lie exactly on the ground
    # surface / domain corner, matching the gsz's own trial-plane endpoints.
    x_daylight = CRIT_X_OVER_H * H
    sd['non_circ'] = [
        {'X': 0.0, 'Y': 0.0, 'Movement': 'Free'},
        {'X': x_daylight, 'Y': H, 'Movement': 'Free'},
    ]
    path = os.path.join(OUT, 'gs2_26.xlsx')
    save_slope_data_to_xlsx(sd, path)
    return 'gs2_26.xlsx'


BUILDERS = [gs2_26]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
