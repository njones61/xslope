"""Builder for the SEEP/W transient verification example "GeoStudio-PEST --
Multistep Outflow" (SEEPW-T07).

T07's base boundary is a STEPPED SUCTION -- a specified PRESSURE head that is
negative at every stage (IC -0.0728 m, stepping to -0.0932 ... -0.1749 m). It is
carried by a time-varying "head" (plain-Dirichlet) series: the boundary is held
at h(t) at every node of the polyline at all times, so a negative-pressure
(suction) head is applied faithfully. (The base polyline sits at y = 0, so its
total head equals its pressure head -- the series values are the base total
heads.) This is the plain "head" companion to the submerged-only "reservoir"
series; the reservoir path could not hold a suction head (an unsubmerged node
flips to a pressure-head-0 exit face), so T07 requires the "head" type.

--- model as read from the vendor .gsz (read-only oracle, never committed) ---

A classic multistep-outflow laboratory experiment: a small coarse-sand sample
sits on a saturated porous ceramic plate, and the suction applied at the base of
the plate is stepped progressively more negative in five stages.  After each step
the sample drains a little more water toward the base; the ceramic plate (two
orders of magnitude less permeable than the sand) meters the outflow, so the
whole test plays out over ~61 hours.  It exercises XSLOPE's unsaturated storage
term C(psi) under a STEPPED specified-pressure-head boundary -- the transient
stepped-drawdown forward model.  (The example ships a PEST inverse loop that
calibrates the van Genuchten parameters from the measured outflow; that external
optimizer is out of scope -- only the forward seepage solve is ported.)

Model (from the vendor .gsz, read-only oracle -- never committed):
    1-D column, base at y = 0.  Porous plate  y = 0..0.007 m (7 mm); sample
    y = 0.007..0.1234 m.
    Sample = coarse sand: Ksat = 9.009e-4 m/s, theta_s = 0.348, theta_r = 0.029
    -> Sy = 0.319.  The vendor VWC is an 80-point spline; a van Genuchten fit
    (suction kPa -> pressure head m) gives vg_a = 8.913 /m, vg_n = 10.19
    (RMS 1e-4 in effective saturation -- a near-exact fit).  The sample
    desaturates sharply across the operating band (suction 0.9->1.7 kPa gives
    Se 0.89->0.02), so the moisture-capacity storage term drives the outflow.
    Porous plate: Ksat = 5.833e-6 m/s, high air-entry -> stays saturated (kr = 1)
    over the whole suction range; it is the rate-limiting element.

Boundary + initial conditions:
    Base (y = 0): specified PRESSURE head stepped in five stages --
        -0.0932 (0 -> 46700 s), -0.1136 (-> 81000), -0.1341 (-> 133000),
        -0.1545 (-> 168670), -0.1749 m (-> 219600 s end).
    IC = uniform TOTAL head H = -0.0728 m (the vendor's linear-PWP initial state,
    hydrostatic equilibrium with the pre-test base suction).  Reached with the
    repeated-time step-series idiom: the base series holds H = -0.0728 at t = 0
    (the t = 0 steady solve makes the column uniform, since only the base is
    Dirichlet and the rest is no-flow) and steps to -0.0932 for t > 0.

Oracle: SEEP/W's own solved node.csv (pore-water pressure at every node, every
saved step).  The published external answer is the lab outflow curve (a SLOPE/W-
free scalar the PEST loop fits) -- not a seepage headline number -- so the lock
is the PWP field.  The locked heads are XSLOPE's own solved values at three
elevations at representative save times; the docs tabulate the SEEP/W comparison
and report the achieved deltas honestly (the vG fit reproduces the retention
curve to RMS 1e-4, so the residual is storage-discretization, not SWCC mapping).

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_mso.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_mso.py --locks
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
_GW = 9.81
_WIDTH = 0.006                   # column width (arbitrary; 1-D) [m]
_Y_PLATE = 0.007                 # top of porous plate [m]
_Y_TOP = 0.1234                  # top of sample [m]
_K_SAMPLE = 9.009e-4            # sample saturated conductivity [m/s]
_K_PLATE = 5.8333e-6            # porous-plate conductivity [m/s]
_THETA_S = 0.348
_THETA_R = 0.029
_SY = _THETA_S - _THETA_R        # 0.319
_VG_A = 8.9125                   # van Genuchten alpha [1/m]  (fit to sample VWC)
_VG_N = 10.1875                  # van Genuchten n [-]
_SS = 1e-5                       # small saturated storage [1/m]
_H_IC = -0.072801                # uniform initial total head [m]
_DURATION = 219600.0            # total time [s]
# five-stage base pressure-head schedule (m) and the step-change times (s)
_STEP_T = [46700.0, 81000.0, 133000.0, 168670.0]
_STEP_H = [-0.093215, -0.113625, -0.13408, -0.154472, -0.174902]
_SAVES = [46000.0, 132000.0, 219600.0]
_TARGET_SIZE = 0.004             # 1-D column mesh size [m]


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
    sample = donor_material(sd)
    sample.update(name='Sample (coarse sand)', c=1.0, phi=30.0, gamma=20.0,
                  gamma_sat=20.0, option='mc', u='seep', k1=_K_SAMPLE,
                  k2=_K_SAMPLE, alpha=0.0, unsat='vg', vg_a=_VG_A, vg_n=_VG_N,
                  kr0=1e-3, h0=-0.4, Ss=_SS, Sy=_SY)
    plate = donor_material(sd)
    plate.update(name='Porous plate', c=1.0, phi=30.0, gamma=20.0,
                 gamma_sat=20.0, option='mc', u='seep', k1=_K_PLATE, k2=_K_PLATE,
                 alpha=0.0, kr0=1e-3, h0=-0.4, Ss=_SS, Sy=0.01)
    sd['materials'] = [sample, plate]
    sd['polygons'] = [
        {'mat_id': 0, 'polygon': Polygon(  # sample
            [(0.0, _Y_PLATE), (_WIDTH, _Y_PLATE),
             (_WIDTH, _Y_TOP), (0.0, _Y_TOP)])},
        {'mat_id': 1, 'polygon': Polygon(  # porous plate
            [(0.0, 0.0), (_WIDTH, 0.0), (_WIDTH, _Y_PLATE), (0.0, _Y_PLATE)])}]
    # Starting circle: the 'circles' sheet is unused by this deck's locked
    # analysis (type=tseep_head reads only the transient seepage field), but
    # xslope requires circles[0] to actually slice. This is a 1-D column
    # (ground is the two points (0,H)-(W,H), no slope): a circle here can only
    # daylight twice on that short top segment if its half-span stays well
    # inside the 0.006 m column width, which in turn bounds how deep it can
    # reach (past roughly half the width the arc bulges out past the vertical
    # side walls and the domain-containment check rejects it). Half-span
    # 0.0024 m, nadir 0.00192 m below the surface -- comfortably inside that
    # bound, with margin from the ~0.003 m limit.
    sd['circles'] = [{'Xo': 0.003, 'Yo': 0.12394, 'Depth': 0.12148, 'R': 0.00246}]
    return sd


def _step_series():
    """Base pressure-head series as repeated-time step anchors; the leading
    t=0 pair carries the uniform IC (H_IC) then jumps to the first stage."""
    times = [0.0, 0.0]
    vals = [_H_IC, _STEP_H[0]]
    edges = _STEP_T + [_DURATION]
    for i, hv in enumerate(_STEP_H):
        t_end = edges[i]
        times.append(t_end)
        vals.append(hv)
        if i + 1 < len(_STEP_H):
            times.append(t_end)          # jump to next stage at the same instant
            vals.append(_STEP_H[i + 1])
    return times, vals


def _mso_sd():
    """The T07 model: stepped-suction base as a plain 'head' (Dirichlet) series."""
    sd = _base_sd()
    times, vals = _step_series()
    sd['seepage_bc'] = {'specified_heads': [
        # base (y=0): plain time-varying head series (a suction Dirichlet held at
        # every base node at all times -- NOT submerged-only reservoir).
        {'head': 'base', 'kind': 'head', 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}],
        'exit_face': []}
    sd['tseep'] = {'times': times, 'series': {'base': vals},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    return sd


def gs2_mso():
    """SEEPW-T07: stepped-suction multistep-outflow forward drainage model."""
    save_slope_data_to_xlsx(_mso_sd(), os.path.join(OUT, 'gs2_mso.xlsx'))
    return 'gs2_mso.xlsx'


BUILDERS = [gs2_mso]


def _solve():
    """Mesh + march the T07 model; return (seep_data, solution)."""
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import build_seep_data, build_tseep_data, run_transient_seepage
    sd = _mso_sd()
    polygons = get_material_polygons(sd)
    mesh = build_mesh_from_polygons(polygons, _TARGET_SIZE, 'tri3')
    seep_data = build_seep_data(mesh, sd)
    solution = run_transient_seepage(seep_data, build_tseep_data(sd), verbose=False)
    return seep_data, solution


def _print_locks():
    """Solve and print the tseep_head test tags: XSLOPE's OWN solved total head at
    three elevations up the column at each save time (the honest-fallback lock the
    parked notes specify -- the vendor node.csv is the read-only oracle the docs
    compare to; the lock is a regression guard on XSLOPE's own field)."""
    import numpy as np
    from xslope.seep import transient_frame_index
    seep_data, solution = _solve()
    nodes = seep_data['nodes']
    xq = _WIDTH / 2.0
    ys = [0.02, 0.06, 0.10]                 # three elevations up the column

    def head_at(h, yq):
        d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
        idx = np.argsort(d2)[:4]
        w = 1.0 / np.maximum(d2[idx], 1e-12)
        return float(np.sum(w * h[idx]) / np.sum(w))

    lines = []
    for t in _SAVES:
        fi = transient_frame_index(solution, t)
        h = np.asarray(solution['frames'][fi]['head'], dtype=float)
        pts = ';'.join(f'{xq:g}:{y:g}:{head_at(h, y):.4f}' for y in ys)
        lines.append(
            f'<!-- test: file=files/geostudio/gs2_mso.xlsx, type=tseep_head, '
            f'target_size={_TARGET_SIZE:g}, time={t:g}, points={pts}, '
            f'tolerance=0.01, benchmark=SEEPW-T07-t{int(round(t))} -->')
    print('\n'.join(lines))


if __name__ == '__main__':
    if '--locks' in sys.argv:
        _print_locks()
    else:
        os.makedirs(OUT, exist_ok=True)
        print('built', gs2_mso())
