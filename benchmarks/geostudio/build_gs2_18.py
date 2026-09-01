"""Builder for GeoStudio SLOPE/W Verification Manual sec. 2.18 -- Borges and
Cardoso - Geosynthetic Embankment #2 (docs/verification/geostudio.md,
"2.18 -- Borges And Cardoso - Geosynthetic Embankment #2").

Borges and Cardoso (2002), their Case 2: a geosynthetic-reinforced embankment
placed over soft clay whose undrained shear strength varies with depth. The
manual (p.42-43, Figure 62 / Table 55 / Table 56) reproduces this with
SLOPE/W's Morgenstern-Price solver and reports FS = 1.171 (M-P, moment = force),
Bishop 1.170, Janbu 1.233; Borges and Cardoso's own published value is 1.15.

Geometry and strengths are taken from the manual's printed Table 55 and Figure
62, and pinned to the vendor .gsz (dev-side cross-check only, never
committed/redistributed -- see xslope.geostudio.read_gsz), whose region
polygons store the model's boundaries exactly:
    Embankment  (0,0)->(6.15,5)->(36.41,5)->(42.56,0), gamma=20, c'=0, phi'=35
    Clay 1  top y= 0.0 (to -1.0)   Su = 33            constant     gamma=17
    Clay 2  top y=-1.0 (to -4.5)   Su = 16            constant     gamma=17
    Clay 3  top y=-4.5 (to -7.0)   Su = 16 + 0.96*z              gamma=17
    Clay 4  top y=-7.0 (to -21.0)  Su = 18.4 + 1.314*z           gamma=17
with z = depth below the layer top. The ground surface is flat at y=0 from
x=-30 to the embankment toe (0,0), up the 5 m embankment, and flat again to
x=45; the domain base is at y=-21.

The two rate-of-increase clays are GeoStudio "Su-varying-with-depth" (SFnDepth)
materials: CTopOfLayer and CRateOfIncrease map exactly onto xslope's 'cp'
option, Su = c + cp*max(0, r_elev - y), with c = CTopOfLayer, cp =
CRateOfIncrease and r_elev = the layer top. (Table 55's "cu bottom" of 55.1 for
Clay 4 is the value extrapolated to the deposit's nominal base well below this
model's floor; within the -7..-21 model the rate 1.314 is what governs, and it
is reproduced exactly.)

The geosynthetic (Table 55 / sec. 2.18 text) is Tmax = 200 kN/m, interface
friction 33.7 deg, NOT anchored and NO adhesion -- i.e. it develops force by
frictional pullout only. It is laid at y=1.0 from x=1.23 to x=41.33 (gsz). It
is modelled as an 'axial'/'active' geosynthetic (the convention
xslope.geostudio uses for GeoStudio's Geosynthetic type, which reproduces
SLOPE/W's own factor of safety). The frictional pullout is carried as a finite
bond length lp; on the critical surface the bar is fully embedded on both sides
of the crossing, so the full 200 kN/m acts and FS is insensitive to the exact
bond length (verified over lp = 0..3 m -- all give FS = 1.153). lp = 1.9 m is a
representative value for delta = 33.7 deg under the crest overburden.

Verified with Morgenstern-Price (half-sine) on SLOPE/W's own critical circle
(Xo=3.47, Yo=5.66, R=13.28), which a grid search confirms is also xslope's
critical: xslope M-P = 1.153 vs SLOPE/W 1.171 (-1.5%) and B&C 1.15 (+0.3%);
xslope Bishop 1.154 vs SLOPE/W 1.170. See docs/verification/geostudio.md
sec. 2.18.

Run: PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_18.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from _gs2_donor import donor_material, load_donor  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# SLOPE/W's own critical circle (from the solved .gsz), also xslope's critical.
CRIT_CIRCLE = {'Xo': 3.465515637162291, 'Yo': 5.664150600821959,
               'R': 13.276107294827842}


def gs2_18():
    """Borges & Cardoso (2002) Case 2 geosynthetic-reinforced embankment on
    depth-varying soft clay (manual Table 55 / Figure 62 / Table 56)."""
    sd = load_donor(ACADS_1A)
    base = donor_material(sd)

    def mat(name, opt, c, phi, g, cp=0.0, r_elev=0.0):
        m = dict(base)
        m.update(name=name, option=opt, c=c, phi=phi, cp=cp, r_elev=r_elev,
                 gamma=g, gamma_sat=g, u='none')
        return m

    sd['materials'] = [
        mat('Embankment', 'mc', 0.0, 35.0, 20.0),
        mat('Clay 1', 'cp', 33.0, 0.0, 17.0, cp=0.0,   r_elev=0.0),
        mat('Clay 2', 'cp', 16.0, 0.0, 17.0, cp=0.0,   r_elev=-1.0),
        mat('Clay 3', 'cp', 16.0, 0.0, 17.0, cp=0.96,  r_elev=-4.5),
        mat('Clay 4', 'cp', 18.4, 0.0, 17.0, cp=1.314, r_elev=-7.0),
    ]

    # Profile lines top-to-bottom (each is the top boundary of its material).
    # The embankment band has zero thickness left/right of the toe, so the flat
    # foundation surface at y=0 is Clay 1.
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(-30.0, 0.0), (0.0, 0.0), (6.15, 5.0),
                                 (36.41, 5.0), (42.56, 0.0), (45.0, 0.0)]},
        {'mat_id': 1, 'coords': [(-30.0, 0.0), (45.0, 0.0)]},
        {'mat_id': 2, 'coords': [(-30.0, -1.0), (45.0, -1.0)]},
        {'mat_id': 3, 'coords': [(-30.0, -4.5), (45.0, -4.5)]},
        {'mat_id': 4, 'coords': [(-30.0, -7.0), (45.0, -7.0)]},
    ]
    sd['max_depth'] = -21.0
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    c = CRIT_CIRCLE
    sd['circles'] = [{'Xo': c['Xo'], 'Yo': c['Yo'], 'R': c['R'],
                      'Depth': c['Yo'] - c['R']}]

    # Geosynthetic: Tmax=200, frictional pullout (delta=33.7 deg, no adhesion,
    # not anchored). Finite bond length; FS is insensitive to it here.
    lines = [{'x1': 1.23, 'y1': 1.0, 'x2': 41.33, 'y2': 1.0, 't_max': 200.0,
              't_res': float('nan'), 'lp1': 1.9, 'lp2': 1.9,
              'tend1': 0.0, 'tend2': 0.0, 'E': 2e4, 'area': 0.1,
              'label': 'geosynthetic', 'type': 'geosynthetic',
              'dir': 'axial', 'appl': 'active', 'spacing': 1.0}]
    sd['reinforcement_lines'] = lines
    sd['reinforce_lines'] = lines

    path = os.path.join(OUT, 'gs2_18.xlsx')
    save_slope_data_to_xlsx(sd, path)
    return 'gs2_18.xlsx'


BUILDERS = [gs2_18]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
