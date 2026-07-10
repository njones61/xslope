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


BUILDERS = [vp002]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
