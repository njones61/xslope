"""Builder for the SEEP/W transient verification example "Leakage from pond with
clay liner" (SEEPW-T04) -> docs/verification/geostudio.md, "Transient seepage
(SEEP/W)".

The unconfined two-dimensional transient with an exit face -- the census flagship
for XSLOPE's storage-driven water-table RISE on a LINEAR mesh (the solver's
quadratic-exit-face caveat asks for tri3/quad4 here).  A clay-lined pond on a
hillside is filled to a constant level; water leaks down through the low-
permeability liner into the embankment and the phreatic surface rises over ~240
days from its initial far-field position toward a new near-steady leaking state,
draining out the downstream seepage face.

Model (from the vendor .gsz, read-only oracle -- never committed; a symmetric
half-model, x = 0 the pond centre-line):
    Embankment body (SatUnsat): Ksat = 1.157e-6 m/s, theta_s = 0.35, theta_r =
    0.032 -> Sy = 0.318.
    Clay liner (SatUnsat): Ksat = 9.259e-8 m/s (~12x less than the fill),
    theta_s = 0.45, theta_r = 0.131 -> Sy = 0.319.  The liner is the thin wedge
    under the pond floor.
    Both VolWCFns carry Beta (mv) = 0.001 /kPa, so Ss = gamma_w*mv = 9.81e-3 /m.

    THE VAN GENUCHTEN PAIRS ARE FITTED TO THE VENDOR'S CONDUCTIVITY TABLES, not to
    its retention splines.  SEEP/W ships theta(psi) and k(psi) as two independent
    20-point tables; XSLOPE has one (a, n) that drives both kr and the moisture
    capacity, so the fit has to choose which table to honour, and kr is what this
    problem turns on -- the liner is the throttle and its k(psi) at the ~7 kPa the
    vendor's own write-up says prevails beneath the pond sets the leakage rate.
    Fitted to k(psi) (fill vg_a = 0.6611, vg_n = 1.9883; liner vg_a = 0.1675,
    vg_n = 1.6029) the relative-conductivity residual is 0.004 and 0.011 decades
    rms against 0.059 and 0.151 for a retention fit, at a cost of 0.016 and 0.079
    RMS in effective saturation against 0.008 and 0.022.  Measured on the answer:
    the 240-day interior heads move from +0.35 m to +0.12 m on this change alone.

    MESH: the vendor sets a 0.5 m global edge length with a RelativeLength 0.25
    constraint on Regions-1 (the liner), i.e. 0.125 m through the 0.25 m liner,
    and its write-up says why ("in order to simulate the influence of the clay
    liner on the movement of water accurately").  A uniform 0.5 m tri3 mesh puts
    ELEVEN triangles in the entire liner with zero nodes interior to it -- the
    whole head drop the problem is about carried by one layer of constant-gradient
    elements.  The tags therefore run refine_factor=2, whose thin-zone field lands
    the liner at a 0.125 m mean element edge (145 triangles, 42 interior nodes),
    the vendor's own resolution.

Boundary + initial conditions:
    Far-field water table (downstream toe, y = 4): specified head = 4 (constant).
    Downstream slope face: potential seepage (exit) face.
    Pond floor (top of the liner, y = 10 -> 10.5): the reservoir series "pond".
    IC = the pre-fill steady state (pond series at head 4, BELOW the pond floor
    at y = 10, so the floor nodes are unsubmerged -> inactive exit faces and only
    the far water table at el. 4 sets the field -> uniform total head 4).  For
    t > 0 the pond series steps to head 10.5 (the floor submerges -> Dirichlet
    h = 10.5) and the pond leaks.  This is the same submerged-only reservoir
    series idiom the rapid-drawdown port (T03) uses, run in the filling direction.

Oracle: SEEP/W's own solved node.csv (pore-water pressure at every node, every
saved step).  The published external answer is a water-table-rise-vs-time graph
(no closed form), so the seepage lock is the PWP/head field.  The locked values
are XSLOPE's own solved total heads at interior stations at the IC and the near-
steady leaking end state; the docs tabulate the SEEP/W comparison and report the
achieved deltas honestly (the recurring SWCC-mapping caveat perturbs the transient
timing, bounded by the two near-steady end members that anchor the locks).

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_pond.py
      PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_pond.py --locks
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
from _gs2_donor import donor_material  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- model constants --------------------------------------------------------
_GW = 9.81                       # gamma_w [kN/m3]
_K_FILL = 1.157e-6              # embankment saturated conductivity [m/s]
_K_LINER = 9.259e-8            # clay-liner saturated conductivity [m/s]
_SS = 9.81e-3                    # gamma_w * Beta(mv=0.001 /kPa) [1/m]
_RES_FULL = 10.5                # pond total head [m]
_WT = 4.0                        # far-field water-table elevation [m]
_DURATION = 20736000.0         # 240 days [s]
_SAVES = [2081483.0, 20736000.0]   # vendor saved steps 012 / 030 [s]
_TARGET = 0.5                    # mesh target size [m] (linear tri3)
_REFINE = 2.0                    # thin-zone refinement -> 0.125 m in the liner

# embankment fill: vG fitted to the vendor "Embankment K" table
_F_TS, _F_TR = 0.35, 0.032
_F_A, _F_N = 0.6611, 1.9883
# clay liner: vG fitted to the vendor "Clay Liner K" table
_L_TS, _L_TR = 0.45, 0.131
_L_A, _L_N = 0.1675, 1.6029

# geometry (vendor points; half-model, x = 0 is the pond centre-line)
_LINER = [(0.0, 10.0), (1.75, 10.0), (2.15, 10.5), (2.95, 11.5),
          (3.25, 11.5), (1.9, 9.75), (0.0, 9.75)]
_FILL = [(0.0, 9.75), (1.9, 9.75), (3.25, 11.5), (7.0, 11.5),
         (21.0, 4.0), (30.0, 4.0), (30.0, 0.0), (0.0, 0.0)]
_POND_FLOOR = [(0.0, 10.0), (1.75, 10.0), (2.15, 10.5)]   # Lines 1 + 11
_WT_FACE = [(21.0, 4.0), (30.0, 4.0)]                     # Line 7
_SEEP_FACE = [(21.0, 4.0), (7.0, 11.5)]                   # Line 13


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
    fill = donor_material(sd)
    fill.update(name='Embankment', c=5.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                option='mc', u='seep', k1=_K_FILL, k2=_K_FILL, alpha=0.0,
                unsat='vg', vg_a=_F_A, vg_n=_F_N, kr0=1e-3, h0=-0.4,
                Ss=_SS, Sy=_F_TS - _F_TR)
    liner = donor_material(sd)
    liner.update(name='Clay liner', c=5.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
                 option='mc', u='seep', k1=_K_LINER, k2=_K_LINER, alpha=0.0,
                 unsat='vg', vg_a=_L_A, vg_n=_L_N, kr0=1e-3, h0=-0.4,
                 Ss=_SS, Sy=_L_TS - _L_TR)
    sd['materials'] = [fill, liner]
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(_FILL)},
                      {'mat_id': 1, 'polygon': Polygon(_LINER)}]
    sd['circles'] = [{'Xo': 15.0, 'Yo': 40.0, 'Depth': 0.0, 'R': 40.0}]
    sd['seepage_bc'] = {'specified_heads': [
        {'head': 'pond', 'kind': 'reservoir', 'coords': _POND_FLOOR},
        {'head': _WT, 'coords': _WT_FACE}], 'exit_face': _SEEP_FACE}
    return sd


def gs2_pond():
    """SEEPW-T04: clay-lined pond fills; water table rises through the liner."""
    sd = _base_sd()
    # pond series: head 4 at t=0 (below the y=10 floor -> unsubmerged -> IC has
    # no pond, only the far water table) then step to 10.5 for t>0 (pond fills).
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'pond': [_WT, _RES_FULL]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(_SAVES)}
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gs2_pond.xlsx'))
    return 'gs2_pond.xlsx'


def _solve_own(save_times):
    import contextlib
    import io as _io
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import (build_seep_data, build_tseep_data,
                             run_transient_seepage, transient_frame_index)
    sd = _base_sd()
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'pond': [_WT, _RES_FULL]},
                   'duration': _DURATION * 1.0001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': list(save_times)}
    polys = get_material_polygons(sd)
    with contextlib.redirect_stdout(_io.StringIO()):
        mesh = build_mesh_from_polygons(polys, _TARGET, 'tri3',
                                        refine_factor=_REFINE)
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
    stations = [(5.0, 2.0), (10.0, 3.0), (15.0, 4.0)]
    lock_times = [0.0, 20736000.0]
    nodes, frames = _solve_own(lock_times)
    for t in lock_times:
        h = frames[t]
        pts = ';'.join(f'{x:g}:{y:g}:{_head_at(nodes, h, x, y):.4f}'
                       for x, y in stations)
        tag = 'ic' if t == 0.0 else 'end'
        print(f'<!-- test: file=files/geostudio/gs2_pond.xlsx, type=tseep_head, '
              f'target_size={_TARGET:g}, refine_factor={_REFINE:g}, time={t:g}, '
              f'max_head_change_frac=0.05, '
              f'points={pts}, tolerance=0.03, benchmark=SEEPW-POND-{tag} -->')


BUILDERS = [gs2_pond]

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    if '--locks' in sys.argv:
        _print_locks()
    else:
        for b in BUILDERS:
            print('built', b())
