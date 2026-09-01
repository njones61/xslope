"""Builder for the SEEP/W transient verification example "Simulating
consolidation with SEEP/W" (SEEPW-T01) -> docs/verification/geostudio.md,
"Transient seepage (SEEP/W)" section.

This is the SEEP/W engineering example that reproduces Terzaghi one-dimensional
consolidation with a pure seepage (uncoupled) transient run: a saturated column
with a uniform initial excess pore-water pressure that dissipates by drainage.
It is the ideal end-to-end test of XSLOPE's saturated storage term (S = Ss),
because the problem is linear (no unsaturated nonlinearity) and has a
closed-form target.

Model (from the vendor .gsz, read-only oracle -- never committed):
    1-D clay column 0.05 m tall, 51 nodes (SatOnly), K = 1e-8 m/s,
    mv = 0.005 /kPa  ->  Ss = mv*gamma_w = 0.005*9.81 = 0.04905 /m,
    cv = K/Ss = 2.039e-7 m2/s.  Uniform initial excess PWP u0 = 10 kPa (the
    GeoStudio InitPWP attribute).  Both ends drained (PWP = 0): DOUBLE drainage,
    drainage path H = 0.025 m.

Target -- closed form (Terzaghi 1943, Eq. 17.3):
    ue(Z,Tv)/u0 = sum_m (2/M) sin(M Z) exp(-M^2 Tv),  M = (pi/2)(2m+1),
    Z = z/H (z from a drained face), Tv = cv t / H^2.
At the published times t = 150 / 604 / 1460 s the average degree of
consolidation is U = 25 / 50 / 75 %.  These analytical excess pressures are the
locked values.  SEEP/W's own solved node.csv is a second, independent oracle;
the docs report agreement with both.  SEEP/W lags Terzaghi at late time (its own
U = 24 / 48 / 69 %) because it runs only 10 exponential steps -- XSLOPE, stepping
finely, sits on the closed form.

Modelled saturated on a datum offset H_ref = 100 m (the same device the GW15-21
consolidation ports use): the total head is held uniform at H_ref + u0/gamma_w
at t = 0 (the t = 0 steady solve sets the uniform IC) and both faces step to
H_ref for t > 0, so the excess head h - H_ref = (u0/gamma_w)*ue/u0 follows
Terzaghi and the excess PWP gamma_w*(h - H_ref) = u0*ue/u0 is directly
comparable to SEEP/W's node.csv in kPa.

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_cons.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_cons.py --locks
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from _gs2_donor import donor_material, load_donor  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- model constants --------------------------------------------------------
_GW = 9.81                       # gamma_w [kN/m3]
_H_REF = 100.0                   # inert datum offset (keeps the column saturated)
_MV = 0.005                      # coeff. of volume compressibility [1/kPa]
_K = 1e-8                        # saturated conductivity [m/s]
_U0 = 10.0                       # uniform initial excess PWP [kPa]
_HEIGHT = 0.05                   # column height [m]
_WIDTH = 0.01                    # column width (arbitrary; 1-D) [m]
_H = _HEIGHT / 2.0               # double-drainage path [m]
_SAVES = [150.0, 604.0, 1460.0]  # published U = 25/50/75 % times [s]


def _ss():
    return _MV * _GW


def _cv():
    return _K / _ss()


def terzaghi_ue(Z, Tv, nterms=800):
    """Terzaghi (1943) Eq. 17.3 excess-PWP ratio ue/u0 at normalized depth Z
    (z/H from a drained face, Z in [0,2] double drainage) and time factor Tv."""
    s = 0.0
    for m in range(nterms):
        M = (np.pi / 2.0) * (2 * m + 1)
        s += (2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv)
    return s


def _base_sd():
    sd = load_donor(ACADS_1A)
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
    m.update(name='Clay', c=1.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
             option='mc', u='seep', k1=_K, k2=_K, alpha=0.0, kr0=1e-3,
             h0=-0.4, Ss=_ss(), Sy=0.1)
    sd['materials'] = [m]
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(
        [(0.0, 0.0), (_WIDTH, 0.0), (_WIDTH, _HEIGHT), (0.0, _HEIGHT)])}]
    # Starting circle: the 'circles' sheet is unused by this deck's locked
    # analysis (type=tseep_head reads only the transient seepage field), but
    # xslope requires circles[0] to actually slice. This is a 1-D column
    # (ground is the two points (0,H)-(W,H), no slope): the deepest circle
    # whose chord still fits within the column width and daylights twice on
    # the top -- a shallow sliver near the crest (chord half-span 0.0045 m
    # against a 0.01 m column, depth 0.0025 m below the surface).
    sd['circles'] = [{'Xo': 0.005, 'Yo': 0.05195, 'Depth': 0.0475, 'R': 0.00445}]
    return sd


def gs2_cons():
    """SEEPW-T01: 1-D Terzaghi consolidation, double drainage."""
    sd = _base_sd()
    base = _H_REF + _U0 / _GW
    # top (y=H) and bottom (y=0) both drained: series holds base at t=0 (uniform
    # IC via the t=0 steady solve) then steps to H_ref for t>0.
    sd['seepage_bc'] = {'specified_heads': [
        {'head': 'res', 'coords': [(0.0, _HEIGHT), (_WIDTH, _HEIGHT)]},
        {'head': 'res', 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}], 'exit_face': []}
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'res': [base, _H_REF]},
                   'duration': max(_SAVES) * 1.001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_cons.xlsx'))
    return 'gs2_cons.xlsx'


def _print_locks():
    """Print the docs test tags locking the Terzaghi closed-form excess heads."""
    cv = _cv()
    u0h = _U0 / _GW
    ys = (0.0125, 0.025, 0.0375)
    lines = []
    for t in _SAVES:
        Tv = cv * t / (_H * _H)
        pts = []
        for y in ys:
            Z = y / _H                      # z from the bottom drained face
            h = _H_REF + u0h * terzaghi_ue(Z, Tv)
            pts.append(f'{_WIDTH / 2.0:g}:{y:g}:{h:.4f}')
        lines.append(
            f'<!-- test: file=files/geostudio/gs2_cons.xlsx, type=tseep_head, '
            f'target_size=0.001, time={t:g}, max_head_change_frac=0.02, '
            f'points={";".join(pts)}, tolerance=0.02, benchmark=SEEPW-CONS-t{t:g} -->')
    print('\n'.join(lines))


BUILDERS = [gs2_cons]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--locks' in sys.argv:
        _print_locks()
    else:
        for b in BUILDERS:
            print('built', b())
