"""Builder for the SEEP/W transient verification example "Mineral heap leaching"
(SEEPW-T05) -> docs/verification/geostudio.md, "Transient seepage (SEEP/W)".

A one-dimensional column of leach ore under an applied surface irrigation flux --
the transient-seepage test of XSLOPE's SPECIFIED-FLUX (Neumann) boundary path and
the van Genuchten moisture-capacity storage in a gravity-drained unsaturated
column. The vendor example carries several analyses; the transient "rate" cases
(2a "Low to high rate with time", 3a "High to low rate with time") are the ones
that march in time. This port builds the low-to-high case (2a): a uniform silty-
sand ore column, initially at a low-rate steady infiltration, is stepped to a
higher irrigation rate at t = 0 and the extra water works its way down.

(The example's OTHER, non-transient, analyses -- "4 Coarse layers" and "5 Fine
layers" -- are the ones that layer the coarse/fine ore; the transient 2a/3a rate
cases both run on the single uniform "Mineral ore" silty sand, as read from the
vendor Contexts. So this is a one-material column, not a layered one.)

Model (from the vendor .gsz, read-only oracle -- never committed):
    1-D column 8 m tall (SatUnsat), Ksat = 5e-6 m/s.  The vendor VWC "Ore W/C -
    silty sand" (theta_s = 0.5, theta_r = 0.0251 -> Sy = 0.475) is a 20-point
    spline; a van Genuchten fit (suction kPa -> pressure head m) gives vg_a =
    1.393 /m, vg_n = 1.897 (RMS 0.005 in effective saturation).  No Beta (mv) is
    set on the material, so the small saturated storage does not enter (the column
    never saturates); Ss = 1e-4 /m is used as a nominal floor.

Boundary + initial conditions:
    Base (y = 0): specified head = 0 (the free-draining "zero pressure" outlet).
    Top (y = 8): specified irrigation flux q, a Darcy velocity applied downward
    into the column -- bound to the tseep series "leach".  IC = the low-rate
    steady state (q = 3e-7 m/s), reached by the t = 0 steady solve; the series
    then steps to the high rate (q = 3e-6 m/s) for t > 0.

    The low-rate IC is a gravity-drained UNIT-GRADIENT profile: away from the
    base the steady pressure head is the constant psi with kr(psi)*Ksat = q, i.e.
    kr = 0.06 -> psi ~ -0.78 m, with a thin boundary layer relaxing to psi = 0 at
    the drained base.  That nonlinear unsaturated steady state is only found by
    XSLOPE's unsaturated (exit-face) IC solver, not the linear confined one, so
    the column's impermeable right side is declared a potential seepage-exit face:
    it never activates (psi < 0 over the whole run, so it is a no-flow boundary in
    fact) but it routes the t = 0 initial-condition solve through the unsaturated
    Picard solver, which correctly returns the unit-gradient profile.  Without it
    the linear confined IC would return a dry hydrostatic column and the front
    would never advance.

Oracle: SEEP/W's own solved node.csv (pore-water pressure at every node, every
saved step).  The published external answer is a graphical VWC/flow-rate response
(no closed form), so the seepage lock is the PWP/head field, not a headline scalar.
The locked values are XSLOPE's own solved total heads at interior elevations at
the IC and the near-steady end state; the docs tabulate the SEEP/W comparison and
report the achieved deltas honestly (IC and early frames agree within ~0.04 m of
head; the high-rate near-steady end state within ~0.12 m -- the recurring SWCC-
mapping timing caveat, the vG fit vs SEEP/W's tabulated retention/kr spline).

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_heap.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_heap.py --locks
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- model constants --------------------------------------------------------
_GW = 9.81                       # gamma_w [kN/m3]
_K = 5e-6                        # saturated conductivity [m/s]
_THETA_S = 0.5
_THETA_R = 0.0251
_SY = _THETA_S - _THETA_R        # 0.475
_VG_A = 1.3927                   # van Genuchten alpha [1/m]  (fit to Ore silty sand)
_VG_N = 1.8965                   # van Genuchten n [-]
_SS = 1e-4                       # nominal saturated storage [1/m] (column never saturates)
_HEIGHT = 8.0                    # column height [m]
_WIDTH = 0.2                     # column width (arbitrary; 1-D) [m]
_Q_LOW = 3e-7                    # low irrigation flux (IC) [m/s]
_Q_HIGH = 3e-6                   # high irrigation flux (t > 0) [m/s]
_DURATION = 345600.0            # 4 days [s]
_SAVES = [13527.0, 67383.0, 345600.0]   # vendor saved steps 003 / 008 / 015 [s]
_TARGET = 0.2                    # mesh target size [m]


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
    m = dict(sd['materials'][0])
    m.update(name='Mineral ore (silty sand)', c=1.0, phi=30.0, gamma=20.0,
             gamma_sat=20.0, option='mc', u='seep', k1=_K, k2=_K, alpha=0.0,
             unsat='vg', vg_a=_VG_A, vg_n=_VG_N, kr0=1e-3, h0=-0.4,
             Ss=_SS, Sy=_SY)
    sd['materials'] = [m]
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(
        [(0.0, 0.0), (_WIDTH, 0.0), (_WIDTH, _HEIGHT), (0.0, _HEIGHT)])}]
    sd['circles'] = [{'Xo': _WIDTH / 2.0, 'Yo': 10.0, 'Depth': 0.0, 'R': 10.0}]
    return sd


def gs2_heap():
    """SEEPW-T05: 1-D leach-ore column, irrigation flux stepped low -> high."""
    sd = _base_sd()
    # base drained (head 0); top irrigation flux bound to the "leach" series; the
    # impermeable right side is a (never-activating) exit face so the IC solve
    # uses the unsaturated Picard solver -> unit-gradient low-rate steady state.
    sd['seepage_bc'] = {
        'specified_heads': [{'head': 0.0, 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}],
        'specified_fluxes': [{'flux': 'leach',
                              'coords': [(0.0, _HEIGHT), (_WIDTH, _HEIGHT)]}],
        'exit_face': [(_WIDTH, 1.0), (_WIDTH, _HEIGHT)]}
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'leach': [_Q_LOW, _Q_HIGH]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_heap.xlsx'))
    return 'gs2_heap.xlsx'


def _solve_own(save_times):
    """Solve the built model and return (nodes, {t: head_frame})."""
    import contextlib
    import io as _io
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import (build_seep_data, build_tseep_data,
                             run_transient_seepage, transient_frame_index)
    sd = _base_sd()
    sd['seepage_bc'] = {
        'specified_heads': [{'head': 0.0, 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}],
        'specified_fluxes': [{'flux': 'leach',
                              'coords': [(0.0, _HEIGHT), (_WIDTH, _HEIGHT)]}],
        'exit_face': [(_WIDTH, 1.0), (_WIDTH, _HEIGHT)]}
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'leach': [_Q_LOW, _Q_HIGH]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    polys = get_material_polygons(sd)
    with contextlib.redirect_stdout(_io.StringIO()):
        mesh = build_mesh_from_polygons(polys, _TARGET, 'tri3')
    seep = build_seep_data(mesh, sd)
    td = build_tseep_data(sd)
    sol = run_transient_seepage(seep, td, verbose=False)
    nodes = seep['nodes']
    out = {}
    for t in save_times:
        fi = transient_frame_index(sol, t)
        out[t] = np.asarray(sol['frames'][fi]['head'], dtype=float)
    return nodes, out


def _head_at(nodes, h, xq, yq):
    d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
    idx = np.argsort(d2)[:4]
    w = 1.0 / np.maximum(d2[idx], 1e-12)
    return float(np.sum(w * h[idx]) / np.sum(w))


def _print_locks():
    ys = (2.0, 4.0, 6.0)
    xc = _WIDTH / 2.0
    lock_times = [0.0, 345600.0]
    nodes, frames = _solve_own(lock_times)
    for t in lock_times:
        h = frames[t]
        pts = ';'.join(f'{xc:g}:{y:g}:{_head_at(nodes, h, xc, y):.4f}' for y in ys)
        tag = 'ic' if t == 0.0 else 'end'
        print(f'<!-- test: file=files/geostudio/gs2_heap.xlsx, type=tseep_head, '
              f'target_size={_TARGET:g}, time={t:g}, max_head_change_frac=0.05, '
              f'points={pts}, tolerance=0.03, benchmark=SEEPW-HEAP-{tag} -->')


BUILDERS = [gs2_heap]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--locks' in sys.argv:
        _print_locks()
    else:
        for b in BUILDERS:
            print('built', b())
