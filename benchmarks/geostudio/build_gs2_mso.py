"""PARKED / BLOCKED builder for the SEEP/W transient verification example
"GeoStudio-PEST -- Multistep Outflow" (SEEPW-T07).

    *** This model is NOT built into the corpus. It is blocked by the transient
    solver envelope and writes NO .xlsx by default. ***

BLOCK (discovered at build time, 2026-07-24): T07's base boundary is a STEPPED
SUCTION -- a specified PRESSURE head that is negative at every stage (IC
-0.0728 m, stepping to -0.0932 ... -0.1749 m). XSLOPE's transient solver has
exactly one time-varying head BC: the SUBMERGED-ONLY reservoir series
(xslope/seep.py `resolve_bc`), which applies the Dirichlet head only where the
node is submerged (elevation <= h(t), i.e. pressure head >= 0) and flips every
node above the water line to a pressure-head-0 exit face. A base under suction
is never "submerged", so a series head there is silently turned into an exit
face at pressure head 0 -- the applied suction is dropped and the solve returns
a uniform total head H = 0 field (verified: XSLOPE psi = -y at all times,
missing the vendor's -0.07 .. -0.24 m field by 0.07-0.17 m, ~50-100% relative).
A NUMERIC constant head is a true Dirichlet regardless of pressure-head sign, so
a single constant suction would work -- but T07 is a five-stage TRANSIENT, so
the time-varying path is required.

Unblocking needs a solver + schema change (gated on Norm, per the template-edit
rule): an opt-in NON-submerged time-varying head series -- e.g. a
seepage_bc specified-head entry carrying {'submerged': False} that `resolve_bc`
applies as a plain Dirichlet bv[mask] = h(t) for all masked nodes. Default stays
submerged=True so every existing transient tag is bit-identical. The persistence
of that per-BC flag through the .xlsx seepage_bc schema is the piece that needs
approval. Until then T07 stays PLANNED (the SEEPW-T01/T02/T03 unsaturated ports
already cover the vG storage + kr path against clean oracles).

The faithful model is retained below as documented reconnaissance so the port is
one flag away once the capability lands.

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
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from shapely.geometry import Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

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
    sample = dict(sd['materials'][0])
    sample.update(name='Sample (coarse sand)', c=1.0, phi=30.0, gamma=20.0,
                  gamma_sat=20.0, option='mc', u='seep', k1=_K_SAMPLE,
                  k2=_K_SAMPLE, alpha=0.0, unsat='vg', vg_a=_VG_A, vg_n=_VG_N,
                  kr0=1e-3, h0=-0.4, Ss=_SS, Sy=_SY)
    plate = dict(sd['materials'][0])
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
    sd['circles'] = [{'Xo': _WIDTH / 2.0, 'Yo': 1.0, 'Depth': 0.0, 'R': 1.0}]
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


def gs2_mso():
    """SEEPW-T07: stepped-suction multistep-outflow forward drainage model."""
    sd = _base_sd()
    times, vals = _step_series()
    sd['seepage_bc'] = {'specified_heads': [
        {'head': 'base', 'coords': [(0.0, 0.0), (_WIDTH, 0.0)]}],
        'exit_face': []}
    sd['tseep'] = {'times': times, 'series': {'base': vals},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_mso.xlsx'))
    return 'gs2_mso.xlsx'


# NOT registered as a corpus builder: T07 is blocked by the solver envelope
# (stepped-suction base head; see the module docstring). gs2_mso() is retained
# only so the port is one opt-in flag away once a non-submerged time-varying head
# BC lands. It writes a corpus file ONLY if explicitly forced with --force-build.
BUILDERS = []  # intentionally empty -- do not emit a misleading gs2_mso.xlsx

if __name__ == '__main__':
    if '--force-build' in sys.argv:
        os.makedirs(OUT, exist_ok=True)
        print('WARNING: gs2_mso.xlsx solves to a WRONG (uniform H=0) field; the '
              'stepped suction is dropped by the submerged-only head path.')
        print('built', gs2_mso())
    else:
        print('SEEPW-T07 is PARKED/BLOCKED by the transient solver envelope '
              '(stepped-suction base head). No corpus file written. See the '
              'module docstring for the unblock (non-submerged time-varying head '
              'BC, gated on Norm).')
