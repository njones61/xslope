"""Builder for GeoStudio SLOPE/W Verification Manual secs. 2.40-2.43 -- the four staged
rapid drawdown examples (docs/verification/geostudio.md, "2.40 -- Rapid Drawdown --
Walter Bouldin Dam" to "2.43 -- Rapid Drawdown -- Pilarcitos Dam").

SLOPE/W solves each of them by the three stages of Duncan, Wright and Wong (1990),
the same method as XSLOPE's rapid drawdown, and states the problem the same way, so
each workbook is SLOPE/W's own model: its geometry and materials, its two
piezometric surfaces, its drawdown strengths, and its critical circle. The page's
rows solve that circle as a rapid drawdown with SLOPE/W's own method.

Everything below is transcribed from the vendor's model files (``Rapid drawdown with
multi-stage.gsz`` and ``Staged rapid drawdown - CE Example / Pumped Storage Project
Dam / Pilarcitos Dam.gsz``; dev-side cross-checks only, never committed or
redistributed -- see xslope.geostudio.read_gsz), element by element:

  * geometry: each <Region>'s <PointIDs>, resolved through the geometry's <Points>,
    with the material the analysis's <Contexts> gives it, in region order;
  * materials: <StressStrain> <UnitWeight>, <CohesionPrime>, <PhiPrime>, and the
    drawdown strengths <DrawdownCohesion> and <DrawdownPhi> (c_R, phi_R) on each soil
    that does not drain. A free-draining soil has neither, and keeps d = psi = 0;
  * water: the analysis's two <PiezometricSurfaces>, resolved through its own
    <DataPoints>. SLOPE/W reads the higher as the water before drawdown and the
    lower as the water after it (Slope Stability Modeling, sec. 9.3), so the higher
    is piezo line 1 and the lower piezo line 2. <MaterialUsesPiezs> connects every
    material to a surface, so every material takes u = piezo. Unit weight of water
    62.4 (<UnitWaterWeight>), imperial units (ft, pcf, psf);
  * the circle: SLOPE/W's critical circle, the lowest-FS row of the analysis's solved
    slip_surface.csv (SlipCenterX, SlipCenterY, SlipRadiusX).

SLOPE/W's drawdown strengths are the R envelope drawn tangent to the Mohr circles.
SLOPE/W converts them to the K_c = 1 envelope by its Eqs. 77-78 (the same as Duncan,
Wright and Brandon 2014, Eqs. 9.6-9.7), and that envelope is XSLOPE's d and psi:

    d   = c_R cos(phi_R) cos(phi') / (1 - sin(phi_R))
    psi = atan( sin(phi_R) cos(phi') / (1 - sin(phi_R)) )

Three details of the vendor files that the transcription settles:

  * the pumped-storage core and random zone carry <Cohesion> 2000 beside their
    Mohr-Coulomb strength. A Mohr-Coulomb material reads <CohesionPrime>, which they
    do not set, so c' = 0; the 2000 is their c_R;
  * the pumped-storage rockfill carries <DryWeight> 128, a lighter unit weight above
    the water table. SLOPE/W's own solved column weights for this analysis use its
    <UnitWeight> 142 throughout (to 0.0003% over the critical circle), so the
    rockfill weighs 142 here too and gamma_sat is left blank;
  * the Walter Bouldin gravel is SLOPE/W's impenetrable Bedrock, which is XSLOPE's
    elastic option.

SLOPE/W's rapid drawdown factor of safety on its critical circle, and XSLOPE's on the
same circle (rapid drawdown, SLOPE/W's method, 60 slices):

    sec.  dam                  method   SLOPE/W  XSLOPE
    2.40  Walter Bouldin       Spencer  1.0260   1.0227
    2.41  USACE benchmark      Spencer  1.4614   1.4338
    2.42  pumped storage       M-P      1.5496   1.5378
    2.43  Pilarcitos           M-P      1.0499   1.0498

Run: PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_rapid.py
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

from shapely.geometry import Polygon                       # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx           # noqa: E402
from _gs2_donor import donor_material, load_donor           # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

GAMMA_W = 62.4                     # <UnitWaterWeight>, pcf


def kc1_envelope(c_r, phi_r, phi):
    """(d, psi) from the R envelope (c_R, phi_R): SLOPE/W's Eqs. 77-78."""
    sr = math.sin(math.radians(phi_r))
    cos_p = math.cos(math.radians(phi))
    return (c_r * math.cos(math.radians(phi_r)) * cos_p / (1.0 - sr),
            math.degrees(math.atan(sr * cos_p / (1.0 - sr))))


# Each model: materials (name, option, <UnitWeight>, <CohesionPrime>, <PhiPrime>,
# c_R, phi_R), the regions as (material index, vertices) in region order, the two
# surfaces, and SLOPE/W's critical circle (Xo, Yo, R).
MODELS = {
    # Walter Bouldin Dam, 'Walter Bouldin Dam' analysis (Spencer). Critical circle:
    # slip 83 of 'Walter Bouldin Dam/001/slip_surface.csv', SLOPE/W FS 1.02596.
    'gs2_40.xlsx': {
        'materials': [
            ('Riprap',              'mc',      140.0,   0.0, 40.0,   0.0,  0.0),
            ('Clayey Silty Sand',   'mc',      128.0, 240.0, 32.7, 650.0, 13.0),
            ('Micoceous Silt',      'mc',      123.0, 220.0, 22.5, 450.0, 11.0),
            ('Cretoceous Clay',     'mc',      124.0, 180.0, 19.0, 230.0, 13.0),
            ('Clayey Sandy Gravel', 'elastic', 128.0,   0.0,  0.0,   0.0,  0.0),
        ],
        'regions': [
            (0, [(104.0, 246.0), (141.0, 265.0), (146.0, 265.0), (109.0, 246.0)]),
            (1, [(0.0, 205.0), (104.0, 246.0), (109.0, 246.0), (146.0, 265.0),
                 (180.0, 265.0), (203.0, 255.0), (146.0, 255.0), (117.0, 242.0),
                 (87.5, 235.0), (55.0, 222.0), (12.5, 205.0)]),
            (2, [(203.0, 255.0), (216.0, 250.0), (216.0, 235.0), (87.5, 235.0),
                 (117.0, 242.0), (146.0, 255.0)]),
            (3, [(87.5, 235.0), (55.0, 222.0), (216.0, 222.0), (216.0, 235.0)]),
            (4, [(216.0, 222.0), (216.0, 205.0), (12.5, 205.0), (55.0, 222.0)]),
        ],
        'before': [(0.0, 252.0), (211.0, 252.0)],                  # surface 1
        'after': [(0.0, 220.0), (216.0, 220.0)],                   # surface 2
        'circle': (77.5, 316.5, 94.49652000000002),
    },
    # USACE benchmark (EM 1110-2-1902 App. G), 'Benchmark Example' (Spencer). The
    # analysis saves one trial, USACE's circle as SLOPE/W placed it: slip 1 of
    # 'Benchmark Example/001/slip_surface.csv', SLOPE/W FS 1.46138.
    'gs2_41.xlsx': {
        'materials': [
            ('New Material', 'mc', 135.0, 0.0, 30.0, 1200.0, 16.0),
        ],
        'regions': [
            (0, [(400.0, 110.0), (312.0, 110.0), (294.5, 103.0), (222.0, 74.0),
                 (72.0, 24.0), (0.0, 0.0), (400.0, 0.0)]),
        ],
        'before': [(0.0, 103.0), (400.0, 103.0)],                  # surface 1
        'after': [(0.0, 24.0), (72.0, 24.0), (400.0, 24.0)],       # surface 2
        'circle': (169.95, 219.47, 218.6383392728732),
    },
    # Pumped storage project dam, 'Staged rapid drawdown' (Morgenstern-Price).
    # Critical circle: slip 50 of 'Staged rapid drawdown/001/slip_surface.csv',
    # SLOPE/W FS 1.54962. The rockfill is two regions, one in each shell.
    'gs2_42.xlsx': {
        'materials': [
            ('Compacted Rockfill',     'mc', 142.0, 0.0, 37.0,    0.0,  0.0),
            ('Silty Clay Core',        'mc', 140.0, 0.0, 36.0, 2000.0, 18.0),
            ('Silty Clay Random Zone', 'mc', 140.0, 0.0, 36.0, 2000.0, 18.0),
        ],
        'regions': [
            (0, [(1005.0, 550.0), (870.0, 380.0), (550.0, 380.0), (700.0, 450.0),
                 (770.0, 450.0)]),
            (1, [(1005.0, 550.0), (1030.0, 550.0), (1300.0, 250.0), (770.0, 250.0),
                 (870.0, 380.0)]),
            (2, [(870.0, 380.0), (770.0, 250.0), (0.0, 250.0), (265.0, 320.0),
                 (320.0, 320.0), (490.0, 380.0), (550.0, 380.0)]),
            (0, [(1030.0, 550.0), (1800.0, 250.0), (1300.0, 250.0)]),
        ],
        'before': [(0.0, 545.0), (980.0, 545.0), (1035.0, 545.0), (1300.0, 260.0),
                   (1805.0, 250.0)],                               # surface 1
        'after': [(0.0, 380.0), (870.0, 380.0), (1300.0, 260.0),
                  (1805.0, 250.0)],                                # surface 2
        'circle': (700.0, 835.0, 485.0),
    },
    # Pilarcitos Dam, 'Pilarcitos Dam' (Morgenstern-Price). Critical circle: slip 62
    # of 'Pilarcitos Dam/001/slip_surface.csv', SLOPE/W FS 1.04990. The surface
    # before drawdown starts at x = 1, as the vendor drew it.
    'gs2_43.xlsx': {
        'materials': [
            ('New Material', 'mc', 135.0, 0.0, 45.0, 60.0, 23.0),
        ],
        'regions': [
            (0, [(0.0, 622.0), (100.0, 657.0), (157.06564, 678.82839), (193.0, 692.0),
                 (209.0, 698.0), (250.0, 698.0), (290.0, 680.0), (290.0, 622.0)]),
        ],
        'before': [(1.0, 692.0), (193.0, 692.0), (290.0, 680.0)],  # surface 1
        'after': [(0.0, 657.0), (290.0, 657.0)],                   # surface 2
        'circle': (97.5, 798.75, 146.77129317410814),
    },
}


def build_one(name, spec):
    sd = load_donor(ACADS_1A)
    mats = []
    for mname, opt, gamma, c, phi, c_r, phi_r in spec['materials']:
        d, psi = kc1_envelope(c_r, phi_r, phi) if (c_r or phi_r) else (0.0, 0.0)
        mats.append(donor_material(
            sd, name=mname, option=opt, gamma=gamma, gamma_sat=None, c=c, phi=phi,
            cp=0.0, r_elev=0.0, d=d, psi=psi, u='piezo', ru=0.0))
    sd['materials'] = mats
    sd['polygons'] = [{'mat_id': mid, 'polygon': Polygon(pts)}
                      for mid, pts in spec['regions']]
    sd['profile_lines'] = []
    sd['max_depth'] = None
    sd['unit_system'] = 'imperial'
    sd['gamma_water'] = GAMMA_W
    sd['water_loads'] = 'auto'
    sd['piezo_line'] = list(spec['before'])
    sd['piezo_line2'] = list(spec['after'])
    sd['piezo_phreatic'] = sd['piezo_phreatic2'] = False
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['k_seismic'] = 0.0
    sd['k_vertical'] = 0.0
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['circular'] = True
    sd['non_circ'] = []
    xo, yo, r = spec['circle']
    # Depth is the tangent ELEVATION (Yo - R), the form the circles sheet stores.
    sd['circles'] = [{'Xo': xo, 'Yo': yo, 'R': r, 'Depth': yo - r}]
    path = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, path)
    return path


def build():
    os.makedirs(OUT, exist_ok=True)
    return [build_one(name, spec) for name, spec in MODELS.items()]


if __name__ == '__main__':
    for p in build():
        print('built', os.path.basename(p))
