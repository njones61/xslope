"""Builder for the SEEP/W transient verification example "Verification --
Infiltration into dry soil" (SEEPW-T02) -> docs/verification/geostudio.md,
"Transient seepage (SEEP/W)" section.

The canonical hard-infiltration benchmark: a tall, initially dry column is
ponded at the surface and a sharp wetting front advances downward. It exercises
XSLOPE's UNSATURATED transient path end to end -- the van Genuchten moisture
capacity C(psi) storage term and the vG-Mualem relative permeability kr(psi) --
where the consolidation port (SEEPW-T01) exercised only the saturated Ss term.

Model (from the vendor .gsz, read-only oracle -- never committed):
    1-D column 1 m tall, 101 nodes (SatUnsat).  van Genuchten retention
    theta_s = 0.363, theta_r = 0.186  ->  Sy = theta_s - theta_r = 0.177;
    GeoStudio A = 9.81 kPa, N = 1.53.  Converting the air-entry parameter from
    matric suction (kPa) to pressure head (m): A = 9.81 kPa = 1.0 m, so XSLOPE's
    vg_a (units 1/length) = 1/1.0 = 1.0 /m and vg_n = 1.53.  Ksat = 1e-6 m/s;
    Beta = 1e-5 /kPa  ->  Ss = Beta*gamma_w = 9.81e-5 /m (the small saturated
    storage, active only behind the front).  Duration 46 800 s, 100 steps.

Initial and boundary conditions:
    IC = uniform pressure head psi = -8 m (GeoStudio InitPWP = -78.48 kPa).
    Top (y = 1): ponded, psi = 0.  Bottom (y = 0): held dry at psi = -8 m
    (far field).

The uniform-psi dry IC is reached WITHOUT an explicit h_init: a uniform pressure
head with a unit downward hydraulic gradient IS a steady state (gravity
drainage, div(kr K grad h) = 0), so holding the top at h = -8 + 1 = -7 and the
bottom at h = -8 + 0 = -8 at t = 0 makes the t = 0 steady solve return the
linear h = y - 8, i.e. psi = -8 everywhere.  For t > 0 the top steps to the
ponded head h = 0 + 1 = 1 (psi = 0) while the bottom stays at h = -8, and the
front advances.  This is the same repeated-time step-series idiom the GW15-21
saturated ports use, extended to a second (top) series.

Oracle: SEEP/W's own solved node.csv (pore-water pressure at every node, every
saved step) is the primary transient answer; the published external reference is
the Warrick, Lomen & Yates (1985) semi-analytical profile.  The locked values
are SEEP/W's final-time (t = 46 800 s) total heads at points chosen off the
steepest part of the front (the census flags that XSLOPE's lumped HRZ mass vs
SEEP/W's consistent mass can differ slightly right at a sharp front); the docs
report the achieved deltas honestly.

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_infil.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_infil.py --locks
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from _gs2_donor import donor_material  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- model constants --------------------------------------------------------
_GW = 9.81                       # gamma_w [kN/m3]
_K = 1e-6                        # saturated conductivity [m/s]
_THETA_S = 0.363
_THETA_R = 0.186
_SY = _THETA_S - _THETA_R        # specific yield (drainable porosity) [-]
_VG_A = 1.0                      # van Genuchten alpha [1/m]  (GeoStudio A=9.81 kPa)
_VG_N = 1.53                     # van Genuchten n [-]
_SS = 1e-5 * _GW                 # saturated specific storage [1/m]
_HEIGHT = 1.0                    # column height [m]
_WIDTH = 0.05                    # column width (arbitrary; 1-D) [m]
_PSI0 = -8.0                     # uniform initial pressure head [m]
_DURATION = 46800.0             # total time [s]
_SAVE = 46800.0                  # locked save time [s]


def _base_sd():
    sd = load_slope_data(ACADS_1A)
    sd['gamma_water'] = _GW
    sd['time_unit'] = 'sec'
    sd['unit_system'] = 'si'
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['profile_lines'] = []
    sd['max_depth'] = None
    m = donor_material(sd)
    m.update(name='Sandy Clay Loam', c=1.0, phi=30.0, gamma=20.0,
             gamma_sat=20.0, option='mc', u='seep', k1=_K, k2=_K, alpha=0.0,
             unsat='vg', vg_a=_VG_A, vg_n=_VG_N, kr0=1e-3, h0=-0.4,
             Ss=_SS, Sy=_SY)
    sd['materials'] = [m]
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(
        [(0.0, 0.0), (_WIDTH, 0.0), (_WIDTH, _HEIGHT), (0.0, _HEIGHT)])}]
    # Starting circle: the 'circles' sheet is unused by this deck's locked
    # analysis (type=tseep_head reads only the transient seepage field), but
    # xslope requires circles[0] to actually slice. This is a 1-D column
    # (ground is the two points (0,H)-(W,H), no slope): a circle here can only
    # daylight twice on that short top segment if its half-span stays well
    # inside the 0.05 m column width, which in turn bounds how deep it can
    # reach (past roughly half the width the arc bulges out past the vertical
    # side walls and the domain-containment check rejects it). Half-span
    # 0.02 m, nadir 0.016 m below the surface -- comfortably inside that
    # bound, with margin from the ~0.025 m limit.
    sd['circles'] = [{'Xo': 0.025, 'Yo': 1.0045, 'Depth': 0.984, 'R': 0.0205}]
    return sd


def gs2_infil():
    """SEEPW-T02: 1-D infiltration into a dry van Genuchten column."""
    sd = _base_sd()
    top_ic = _PSI0 + _HEIGHT          # h so that psi = -8 at y = 1 (IC)
    top_pond = 0.0 + _HEIGHT          # h so that psi = 0 at y = 1 (t > 0, ponded)
    bot = _PSI0 + 0.0                 # h so that psi = -8 at y = 0 (held dry)
    sd['seepage_bc'] = {'specified_heads': [
        {'head': 'top', 'coords': [(0.0, _HEIGHT), (_WIDTH, _HEIGHT)]},
        {'head': bot, 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}], 'exit_face': []}
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'top': [top_ic, top_pond]},
                   'duration': _DURATION * 1.001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': [_SAVE]}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_infil.xlsx'))
    return 'gs2_infil.xlsx'


# SEEP/W final-time (t=46 800 s) node.csv oracle, distilled to TOTAL HEAD
# h = psi + y at four points in the WETTED zone behind the front (values embedded
# as locked literals; the .gsz itself is never committed).  These points are
# chosen off the steepest part of the front: XSLOPE reproduces them to within
# 0.05 m of SEEP/W (the mid-front position differs by ~0.04 m, the lumped-HRZ vs
# consistent-mass front diffusion the census predicts).  The tolerance (0.08 m)
# is deliberately wider than the saturated GW page norm (0.05) to carry that
# front-diffusion offset; it is flagged in the docs.
_LOCK_POINTS = [(0.6, 0.4344), (0.7, 0.6358), (0.8, 0.7803), (0.9, 0.8967)]
_LOCK_TOL = 0.08


def _print_locks():
    pts = ';'.join(f'{_WIDTH / 2.0:g}:{y:g}:{h:.4f}' for y, h in _LOCK_POINTS)
    print(f'<!-- test: file=files/geostudio/gs2_infil.xlsx, type=tseep_head, '
          f'target_size=0.01, time={_SAVE:g}, max_head_change_frac=0.01, '
          f'points={pts}, tolerance={_LOCK_TOL:g}, benchmark=SEEPW-INF -->')


BUILDERS = [gs2_infil]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--locks' in sys.argv:
        _print_locks()
    else:
        for b in BUILDERS:
            print('built', b())
