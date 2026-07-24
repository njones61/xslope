"""Builder for the SEEP/W transient verification example "Rapid drawdown"
(SEEPW-T03) -> docs/verification/geostudio.md, "Transient seepage (SEEP/W)".

This is the census flagship for XSLOPE's reservoir feature set: a silty-clay
embankment with a toe drain, initially at steady state under a full reservoir
(head 8 m), is drained two ways --

  * INSTANTANEOUS: the reservoir is removed at t = 0 (head steps 8 -> 0). The
    submerged upstream face becomes a potential seepage (exit) face and the
    trapped pore pressures dissipate over 30 days.
  * SLOW: the reservoir falls 8 -> 0 m linearly over 5 days, then stays empty
    to 30 days.

Both exercise exactly the pair the "submerged-only Dirichlet series head" path
was built for: a boundary node on the upstream face is held at the reservoir
head h(t) only while submerged (y <= h(t)) and flips to an exit-face node once
the falling water line drops below it (xslope/seep.py resolve_bc).

Model (from the vendor .gsz, read-only oracle -- never committed):
    Embankment on a y = 0 base from x = 3 to x = 47 m; crest 10 m tall between
    x = 23 and x = 27; upstream slope (3,0)->(19,8)->(23,10) with the reservoir
    surface at el. 8 (point (19,8)); downstream slope (27,10)->(47,0).  A toe
    drain (region 4) is a wedge below the downstream toe, (39.97,0)-(40,-1)-
    (50,-1)-(50,0)-(47,0), held at total head 0 (free draining).

    Dam fill = "silty clay":  Ksat = 1e-6 m/s (isotropic), theta_s = 0.4,
    theta_r = 0.05 -> Sy = 0.35; Beta (mv) = 1e-4 /kPa -> Ss = Beta*gamma_w =
    9.81e-4 /m.  The vendor VWC is a SampleMatls "Silty Clay" spline; a van
    Genuchten fit to its 20 points (suction kPa -> pressure head m) gives
    vg_a = 0.338 /m, vg_n = 1.846 (RMS 0.007 in effective saturation).  This is
    the recurring SWCC-mapping caveat: the vG fit reproduces the retention curve
    closely but shifts the unsaturated-zone drainage TIMING slightly vs SEEP/W's
    tabulated spline.
    Toe drain: Ksat = 1.157e-5 m/s (~11x the fill), saturated (kr = 1), Ss =
    9.81e-4 /m; its boundary is pinned at total head 0.

Initial condition = the t = 0 steady solve with the reservoir series at 8 m
(the standard repeated-time step-series idiom the GW15-21 / T01 / T02 ports
use): the first two series anchors sit at t = 0, value [8, 0], so the t = 0
steady state is the full-reservoir field and the reservoir then drops.

Oracle: SEEP/W's own solved node.csv (pore-water pressure at every node, every
saved step) is the transient answer (the example's PUBLISHED numbers are FS-vs-
time from a downstream SLOPE/W coupling -- out of scope here, so the seepage
lock is the PWP field, not a headline scalar).  The locked heads are XSLOPE's
own solved values at robust interior stations at the IC and near the drained
end state (both near-steady, so SWCC timing barely enters); the docs tabulate
the SEEP/W comparison and report the achieved deltas honestly.

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_rdd.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_rdd.py --locks
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
_GW = 9.81                       # gamma_w [kN/m3]
_K_FILL = 1e-6                    # dam-fill saturated conductivity [m/s]
_K_DRAIN = 1.157e-5              # toe-drain saturated conductivity [m/s]
_THETA_S = 0.4
_THETA_R = 0.05
_SY = _THETA_S - _THETA_R        # 0.35
_VG_A = 0.3378                   # van Genuchten alpha [1/m]  (fit to Silty Clay)
_VG_N = 1.8462                   # van Genuchten n [-]
_SS = 1e-4 * _GW                 # 9.81e-4 /m  (Beta = mv = 1e-4 /kPa)
_RES_FULL = 8.0                  # full reservoir total head [m]
_DURATION = 2.592e6             # 30 days [s]
_RAMP = 4.32e5                   # slow-drawdown ramp 8->0 [s]  (5 days)

# geometry vertices
_TOE_U = (3.0, 0.0)              # upstream toe
_WL = (19.0, 8.0)                # reservoir surface on the upstream face
_CREST_L = (23.0, 10.0)
_CREST_R = (27.0, 10.0)
_TOE_D = (47.0, 0.0)             # downstream toe
_DRAIN_L = (39.973763, 0.0)      # drain/fill interface, left
_DRAIN = [(39.973763, 0.0), (40.0, -1.0), (50.0, -1.0), (50.0, 0.0), (47.0, 0.0)]


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
    fill = dict(sd['materials'][0])
    fill.update(name='Silty clay (dam fill)', c=5.0, phi=30.0, gamma=20.0,
                gamma_sat=20.0, option='mc', u='seep', k1=_K_FILL, k2=_K_FILL,
                alpha=0.0, unsat='vg', vg_a=_VG_A, vg_n=_VG_N, kr0=1e-3,
                h0=-0.4, Ss=_SS, Sy=_SY)
    drain = dict(sd['materials'][0])
    drain.update(name='Toe drain', c=0.0, phi=35.0, gamma=20.0, gamma_sat=20.0,
                 option='mc', u='seep', k1=_K_DRAIN, k2=_K_DRAIN, alpha=0.0,
                 kr0=1e-3, h0=-0.4, Ss=_SS, Sy=0.3)
    sd['materials'] = [fill, drain]
    sd['polygons'] = [
        {'mat_id': 0, 'polygon': Polygon(
            [_TOE_U, _DRAIN_L, _TOE_D, _CREST_R, _CREST_L, _WL])},
        {'mat_id': 1, 'polygon': Polygon(_DRAIN)}]
    # inert failure surface (transient seepage only; no LEM here)
    sd['circles'] = [{'Xo': 25.0, 'Yo': 30.0, 'Depth': 0.0, 'R': 30.0}]
    # reservoir face + toe-drain-held-at-0 BCs (shared by both cases)
    sd['seepage_bc'] = {'specified_heads': [
        {'head': 'res', 'coords': [_TOE_U, _WL]},
        {'head': 0.0, 'coords': _DRAIN + [_DRAIN_L]}], 'exit_face': []}
    return sd


def gs2_rdd_inst():
    """SEEPW-T03 instantaneous drawdown: reservoir head steps 8 -> 0 at t = 0."""
    sd = _base_sd()
    sd['tseep'] = {'times': [0.0, 0.0],
                   'series': {'res': [_RES_FULL, 0.0]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None,
                   'save_times': [21600.0, 291181.0, 2592000.0]}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_rdd_inst.xlsx'))
    return 'gs2_rdd_inst.xlsx'


def gs2_rdd_slow():
    """SEEPW-T03 slow drawdown: reservoir head falls 8 -> 0 linearly over 5 d."""
    sd = _base_sd()
    sd['tseep'] = {'times': [0.0, _RAMP, _DURATION],
                   'series': {'res': [_RES_FULL, 0.0, 0.0]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None,
                   'save_times': [103638.0, 461858.0, 2592000.0]}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_rdd_slow.xlsx'))
    return 'gs2_rdd_slow.xlsx'


BUILDERS = [gs2_rdd_inst, gs2_rdd_slow]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--locks' in sys.argv:
        pass  # locks are printed by the solve-probe (own-value regression)
    else:
        for b in BUILDERS:
            print('built', b())
