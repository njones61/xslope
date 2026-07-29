"""Builder for GeoStudio SLOPE/W Verification Manual sec. 2.33 -- Priest -
Rigid Blocks (docs/verification/geostudio.md, "2.33 -- Priest - Rigid
Blocks").

Priest (1993), *Discontinuity Analysis for Rock Engineering*: a classical
planar-failure rock slope with a tension crack -- a single rigid block
translating on ONE straight failure plane, not a multi-block/wedge assembly
despite the section title. The manual (printed p.71-72, Figure 106/Table
91/Table 92) reproduces this with SLOPE/W's specified-surface solver and
reports Janbu Simplified 1.049 (Priest's own hand value AND SLOPE/W, exact
match) and Morgenstern-Price 1.049 (SLOPE/W only -- Priest's method has no M-P
analogue). A single Mohr-Coulomb material: c'=20 kPa, phi'=30 deg, gamma=25
kN/m3 (Table 91, printed).

Geometry -- domain and material corners are printed as explicit coordinate
labels on Figure 106, transcribed directly:
    (-25,0) -> (0,0) [toe] -> (17.321,30) [crest, 60 deg face] -> (60,30)
    domain floor at y=-15, right edge at x=60 -- see Figure 106.
The 60 deg face angle is exact: atan(30/17.321) = 60.000 deg.

The FAILURE PLANE and the tension-crack / water-table geometry are NOT
printed as text/table coordinates anywhere in sec. 2.33 -- only drawn, at low
resolution, in Figure 106 (a solid red line and a dashed blue line). Reading
pixel coordinates off the embedded figure (1100x624 source raster) gives a
plane dipping at roughly 30 deg from the toe and a water table kinking near
(20-26, 15-19) -- consistent with, but not precise enough to build from
_alone_. The precise numbers below were confirmed against the vendor .gsz's
own SOLVED SLICE TABLE (dev-side cross-check only, per repo convention --
never committed/redistributed; see xslope.geostudio.read_gsz), which resolves
every pixel-reading ambiguity exactly:

    SLOPE&W Analysis/001/column_1.csv (30 slices, SLOPE/W's own solve):
      slice 1 (uphill end):  XRight=25.981, YBackBotRight=15.0, GroundSurfaceY=30
      slice 30 (toe end):    XLeft=0.0,     YBackBotLeft=0.0

  i.e. the solved failure-plane BASE runs from the toe (0,0) to (25.981, 15.0)
  -- a dead-straight line at slope 15/25.981 = tan(30 deg) EXACTLY. That
  truncation point is not an independent input: it is precisely where a plane
  dipping at 30 deg from the toe crosses "ground surface minus tcrack_depth"
  on the flat crest (30 - 15 = 15) -- i.e. this is xslope's own tcrack_depth
  auto-truncation mechanism (xslope.slice.generate_failure_surface, the
  "left_facing" branch, since the crest is upslope/right here), applied to a
  plane that would otherwise daylight naturally at (30/tan(30deg), 30) =
  (51.9615, 30) on the flat crest. So the surface below is built UNTRUNCATED,
  exactly like every other tcrack problem in the corpus (e.g. Rocscience
  vp002's ACADS 1(b) crack depth), and tcrack_depth=15 does the truncation --
  which reproduces SLOPE/W's own solved (25.981, 15.0) endpoint to 4
  significant figures with no truncation point entered by hand.

  Water table (piezo_line), also confirmed against the same .gsz:
      [(-25,0), (0,0), (25.981,18.75), (60,18.75)]
  18.75 = 15.0 + 3.75 = the tcrack floor elevation (15) plus tcrack_water
  (3.75 = 25% of the 15 m crack depth, matching the manual's "water fills the
  tension crack 25% at the line of failure" and the .gsz's own imported
  tcrack_water=3.75) -- i.e. the water table is flat at the water level
  standing in the crack out to x=25.981 (directly above the crack floor),
  then a straight line down to the toe (0,0), then flat at y=0 out to the
  left domain edge. This is exactly Figure 106's dashed blue line -- the
  short diagonal segment from the toe up to the y=15 boundary that is
  visible in the printed figure is the dip segment described here (the
  horizontal-at-elevation-15ish portion is not separately visible because it
  runs behind the black region boundary at that elevation).

Because the failure surface is a single straight line (not curved), the
Janbu correction factor fo = 1.0 exactly, which is why the manual shows
Priest's uncorrected hand value and SLOPE/W's fo-corrected value agreeing to
3 decimals (1.049 = 1.049) -- and it is why xslope's own (fo-corrected)
`janbu` and Priest's Janbu Simplified should also agree.

Run: PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_33.py
"""

import math
import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx as _write_xlsx  # noqa: E402
from _gs2_donor import donor_material  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# Figure 106 corner coordinates (printed labels)
TOE = (0.0, 0.0)
CREST = (17.321, 30.0)          # top of the 60 deg face
DOMAIN_LEFT = (-25.0, 0.0)
DOMAIN_TOP_RIGHT = (60.0, 30.0)
DOMAIN_BOTTOM = -15.0

TCRACK_DEPTH = 15.0
TCRACK_FRACTION_FULL = 0.25
TCRACK_WATER = TCRACK_DEPTH * TCRACK_FRACTION_FULL     # 3.75

PLANE_DIP_DEG = 30.0             # confirmed against the vendor .gsz's solved slice table
PLANE_SLOPE = math.tan(math.radians(PLANE_DIP_DEG))

# Natural (untruncated) daylight point of the failure plane on the flat crest
# top (y=30): xslope's tcrack_depth mechanism truncates this automatically.
X_DAYLIGHT = DOMAIN_TOP_RIGHT[1] / PLANE_SLOPE          # 51.9615

# Water table kink point: directly above the (auto-truncated) crack floor,
# at the water level standing in the crack.
X_WT_KINK = TCRACK_DEPTH / PLANE_SLOPE                  # 25.981
Y_WT_KINK = TCRACK_DEPTH + TCRACK_WATER                 # 18.75


def gs2_33():
    """Priest (1993) rigid-block planar failure with a 15 m tension crack,
    25% water-filled. Single MC material c'=20, phi'=30, gamma=25 (Table 91).
    Manual: Janbu Simplified 1.049 (Priest / SLOPE/W), M-P 1.049 (SLOPE/W)."""
    sd = load_slope_data(ACADS_1A)
    m = donor_material(sd)
    m.update(name='Material 1', c=20.0, phi=30.0, gamma=25.0, option='mc', u='piezo')
    sd['materials'] = [m]

    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [DOMAIN_LEFT, TOE, CREST, DOMAIN_TOP_RIGHT]},
    ]
    sd['max_depth'] = DOMAIN_BOTTOM

    sd['gamma_water'] = 9.807
    sd['piezo_line'] = [DOMAIN_LEFT, TOE, (X_WT_KINK, Y_WT_KINK), (DOMAIN_TOP_RIGHT[0], Y_WT_KINK)]
    # NOTE: the vendor .gsz tags this piezometric line Type=phreatic, which the
    # xslope importer maps to its own cos^2(local piezo-line slope) inclination
    # correction (xslope.slice, the "Lines declared Type='phreatic'" branch).
    # Cross-checked against SLOPE/W's OWN solved per-slice PoreWaterPressure in
    # column_1.csv: SLOPE/W's values match a PLAIN vertical-head lookup on this
    # line exactly (e.g. slice at X=25.548: SLOPE/W u=36.163 kPa vs an
    # uncorrected hw*gamma_w=36.19; xslope's cos^2-corrected value is 23.78 --
    # 34% low). So SLOPE/W did not apply an inclination correction to this
    # particular (steep, kinked) piezometric line, and piezo_phreatic is left
    # unset (plain static-head 'piezo') to match SLOPE/W's own solve, not the
    # importer's default mapping.

    sd['tcrack_depth'] = TCRACK_DEPTH
    sd['tcrack_water'] = TCRACK_WATER

    # Untruncated planar failure surface, toe -> natural daylight point on the
    # flat crest; tcrack_depth=15 truncates it to SLOPE/W's own (25.981, 15.0)
    # (see module docstring).
    sd['circular'] = False
    sd['circles'] = []
    sd['non_circ'] = [
        {'X': TOE[0], 'Y': TOE[1], 'Movement': 'Free'},
        {'X': X_DAYLIGHT, 'Y': DOMAIN_TOP_RIGHT[1], 'Movement': 'Free'},
    ]

    path = os.path.join(OUT, 'gs2_33.xlsx')
    _write_xlsx(sd, path)
    return path


BUILDERS = [gs2_33]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
    print('X_DAYLIGHT=%.4f  X_WT_KINK=%.4f  Y_WT_KINK=%.4f' % (X_DAYLIGHT, X_WT_KINK, Y_WT_KINK))
