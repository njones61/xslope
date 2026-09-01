"""Builder for GeoStudio SLOPE/W Verification Manual sec. 2.22 -- Cannon Dam #2
(docs/verification/geostudio.md, "2.22 -- Cannon Dam #2").

Hassan and Wolff (1999)'s Cannon Dam reliability study. The paper searches for the
minimum-reliability-index surface and, along the way, prints factors of safety on
NINE fixed circular surfaces -- its Figure 7 surfaces A-E and its Figure 8 surfaces
B, F, G and H. The SLOPE/W verification model stores all nine as nine saved
analyses of one geometry, each solved with Morgenstern-Price, so every surface
comes with the vendor's own factor of safety on the vendor's own geometry.

This file is the geometry plus the nine circles. It is the like-for-like companion
to the Slide2 corpus' vp035.xlsx, which digitizes the same dam from the paper's
printed section: that frame is ~7% anisotropic (a least-squares fit of its ground
surface against this one gives x scale 0.958 against y scale 1.011, 1.6 m rms shape
residual over nine shared vertices), so the paper's fixed circles cannot be carried
onto it -- slice weight lands 33-372% high. The circles belong here, on the exact
frame, and vp035.xlsx keeps its own free-search row.

Geometry, materials and circles are transcribed from the vendor model
``Cannon Dam #2.gsz`` (dev-side cross-check only, never committed or
redistributed -- see xslope.geostudio.read_gsz). Seven material zones over six
Mohr-Coulomb materials; the model is a TOTAL-stress analysis -- the vendor's own
solved slice tables carry pore-water pressure identically zero on all nine
surfaces, and there is no piezometric line in the file.

Vendor per-surface answers, read from each analysis' own ``slip_surface.csv``
(Morgenstern-Price, moment = force at the tabulated lambda):

    surface   SLOPE/W FS   center (x, y)          R         slices   sum W
    Fig7 A       2.5598    (320.175, 313.336)   167.989       34     93921.5
    Fig7 B       2.8062    (295.235, 277.589)   121.291       30     54914.2
    Fig7 C       2.7713    (303.683, 310.206)   154.278       33     59347.6
    Fig7 D       2.5888    (322.152, 310.206)   164.846       32     91402.0
    Fig7 E       2.7031    (327.693, 324.018)   176.855       32     84630.7
    Fig8 B       2.6726    (304.239, 307.858)   152.547       33     60703.5
    Fig8 F       3.5861    (313.114, 304.784)   144.416       30     35874.4
    Fig8 G       6.0694    (314.856, 305.604)   139.848       31     18837.4
    Fig8 H      11.5825    (316.198, 305.604)   135.067       29      7200.2

Run: PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_22.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from shapely.geometry import Polygon                       # noqa: E402
from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from _gs2_donor import donor_material, load_donor  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# (name, c kPa, phi deg, gamma kN/m3) -- the vendor's material list, in its order
MATS = [
    ('Phase I Fill',  117.79,  8.5, 22.0),
    ('Sand Filter',     0.00, 35.0, 22.0),
    ('Phase II Fill', 143.64, 15.0, 22.0),
    ('Material 4',      5.00, 35.0, 25.0),
    ('Material 5',      5.00, 18.0, 20.0),
    ('Material 6',    100.00, 35.0, 20.0),
]

# (mat_id, vertices) -- the vendor's seven region polygons, bottom-up
POLYS = [
    (5, [(0.937, 121.205), (423.373, 120.824), (423.373, 139.853),
         (225.856, 139.853), (204.163, 139.853), (0.937, 139.853)]),
    (4, [(0.937, 156.598), (40.517, 156.598), (174.098, 156.598),
         (204.163, 139.853), (0.937, 139.853)]),
    (4, [(255.921, 155.837), (379.607, 155.837), (423.373, 155.837),
         (423.373, 139.853), (225.856, 139.853)]),
    (3, [(423.373, 155.837), (423.373, 166.112), (349.162, 166.112),
         (357.596, 163.266), (379.607, 155.837)]),
    (0, [(357.596, 163.266), (379.607, 155.837), (255.921, 155.837),
         (225.856, 139.853), (204.163, 139.853), (174.098, 156.598),
         (40.517, 156.598), (66.125, 165.456), (215.2, 165.456),
         (218.625, 165.456)]),
    (1, [(215.2, 193.894), (218.244, 193.894), (218.625, 165.456),
         (215.2, 165.456)]),
    (2, [(320.999, 174.865), (349.162, 166.112), (357.596, 163.266),
         (218.625, 165.456), (218.244, 193.894), (215.2, 193.894),
         (215.2, 165.456), (66.125, 165.456), (101.028, 177.529),
         (158.875, 182.857), (210.633, 198.842), (219.767, 198.842),
         (258.205, 185.141), (306.537, 175.627)]),
]

# Hassan & Wolff's nine fixed surfaces, as SLOPE/W stores them. Order matters:
# the docs tags address a circle by its index on the 'circles' sheet.
CIRCLES = [
    ('Fig7 A', 320.175, 313.336, 167.989),
    ('Fig7 B', 295.235, 277.589, 121.291),
    ('Fig7 C', 303.683, 310.206, 154.278),
    ('Fig7 D', 322.152, 310.206, 164.846),
    ('Fig7 E', 327.693, 324.018, 176.855),
    ('Fig8 B', 304.239, 307.858, 152.547),
    ('Fig8 F', 313.114, 304.784, 144.416),
    ('Fig8 G', 314.856, 305.604, 139.848),
    ('Fig8 H', 316.198, 305.604, 135.067),
]


def build():
    sd = load_donor(ACADS_1A)
    base = donor_material(sd)
    sd['materials'] = [
        dict(base, name=name, option='mc', c=c, phi=phi, gamma=g, gamma_sat=g,
             cp=0.0, r_elev=0.0, u='none', ru=0.0)
        for (name, c, phi, g) in MATS
    ]
    sd['gamma_water'] = 9.807
    sd['piezo_line'] = []
    sd['dloads'] = []
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': mid, 'polygon': Polygon(pts)} for mid, pts in POLYS]
    sd['max_depth'] = None
    sd['circular'] = True
    sd['non_circ'] = []
    # Depth is the tangent ELEVATION (Yo - R), the form the circles sheet stores.
    sd['circles'] = [{'Xo': x, 'Yo': y, 'R': r, 'Depth': round(y - r, 4)}
                     for (_n, x, y, r) in CIRCLES]

    path = os.path.join(OUT, 'gs2_22.xlsx')
    save_slope_data_to_xlsx(sd, path)
    return path


if __name__ == '__main__':
    print('wrote', build())
