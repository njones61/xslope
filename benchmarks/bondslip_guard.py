"""Guard for the BOND-SLIP load-transfer model in solve_fem / solve_ssrm.

The default reinforcement idealization develops the bar's tensile capacity by a fixed
linear end-ramp over lp1/lp2 at each end (fileio.reinforce_available_tension). The
optional bond-slip model (bond_slip={line: (bond_c, bond_phi_deg, perimeter)}) instead
caps the force GRADIENT along the embedded length by a stress-dependent Coulomb bond:

    df/ds <= q(s) = perimeter * max(0, bond_c + sigma_n(s) * tan(bond_phi))

with sigma_n = local vertical overburden. The tension a bar can carry at a point is the
smaller of the two one-sided integrals of q from each free end, still capped by the
material axial capacity t_max (see fem._bond_slip_caps). This guard asserts, partly on
tiny synthetic fem_data (exact envelope math) and partly on a small reinforced slope that
CONVERGES under gravity at F = 1:

  1. INVARIANCE. bond_slip=None is bit-identical to omitting the argument entirely.

  2. ENVELOPE MATH. Under uniform overburden the per-element cap reduces EXACTLY to the
     classical double-ended pull-out ramp q*min(d1, d2), and is capped at t_max where
     the two ramps would exceed it.

  3. AXIAL CAP. A very strong bond (the envelope everywhere above t_max) caps every
     element at t_max.

  4. ENGAGE + CONVERGE. On a stable slope whose reinforcement is not load-critical, a
     tight bond starves the bar: the converged max bar force drops well below the
     end-ramp value and never exceeds the bond cap, while the slope still converges.

  5. GRADIENT MONOTONE. Progressively tighter bonds give monotonically non-increasing
     bar force (finer bond -> less developed capacity).

  6. UNKNOWN LINE REJECTED. solve_ssrm(bond_slip={bad: ...}) raises ValueError rather
     than silently no-opping.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/bondslip_guard.py
Exits non-zero on any failure.
"""
import contextlib
import hashlib
import io
import os
import sys
import tempfile

import numpy as np
from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints)
from xslope.fem import build_fem_data, solve_fem, solve_ssrm, _bond_slip_caps

# A small embankment in a comfortably-stable c-phi soil (converges at F=1 without leaning
# on the reinforcement), with one horizontal geosynthetic layer at mid-height. The bar
# carries a modest tension the bond model can then starve without failing the slope.
_C, _PHI, _E, _NU = 8.0, 24.0, 30_000.0, 0.3
_GEOM = [(0.0, 0.0), (20.0, 0.0), (20.0, 8.0), (8.0, 8.0), (0.0, 0.0)]
_TMAX = 400.0
_MAX_ITER = 5000


def _fixture():
    template = os.path.join(os.path.dirname(__file__), '..', 'docs', 'lem',
                            'files', 'xslope_acads_simple.xlsx')
    sd = load_slope_data(template)
    base = dict(sd['materials'][0])
    mat = dict(base, c=_C, phi=_PHI, gamma=20.0, gamma_sat=20.0, E=_E, nu=_NU,
               u='none', option='mc', psi=0.0, name='soil')
    sd['materials'] = [mat]
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(_GEOM)}]
    sd['max_depth'] = 0.0
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [{'Xo': 8.0, 'Yo': 14.0, 'Depth': 0.0, 'R': 14.0}]
    sd['reinforcement_lines'] = [{
        'x1': 4.0, 'y1': 4.0, 'x2': 18.0, 'y2': 4.0,
        't_max': _TMAX, 't_res': float('nan'), 'lp1': 1.0, 'lp2': 1.0,
        'tend1': 0.0, 'tend2': 0.0, 'E': 2.0e5, 'area': 0.1,
        'label': 'geo', 'spacing': 1.0}]
    sd['reinforce_lines'] = sd['reinforcement_lines']
    tmp = os.path.join(tempfile.gettempdir(), 'bondslip_guard.xlsx')
    save_slope_data_to_xlsx(sd, tmp)
    sd = load_slope_data(tmp)
    cl, _, _ = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=cl)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(polys, target_size=1.0, element_type='tri6',
                                        lines=cl,
                                        point_constraints=extract_point_constraints(sd))
        fem = build_fem_data(sd, mesh)
    return fem


def _digest(sol):
    def h(a):
        return hashlib.sha256(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()[:16]
    return '|'.join([str(sol['converged']), str(sol['iterations']),
                     h(sol['displacements']), h(sol['stresses']),
                     h(sol['forces_1d'])])


def _synthetic():
    """Six unit-length elements along a 6 m line under uniform overburden — the exact
    envelope has a closed form (q*min(d1,d2), capped at t_max)."""
    n = 6
    return {
        'element_materials_1d': np.array([1, 1, 1, 1, 1, 1]),
        't_allow_by_1d_elem': np.full(n, _TMAX),
        'reinforce_line_labels': ['geo'],
        'elem_length_1d': np.full(n, 1.0),
        'dist_end1_1d': np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5]),
        't_max_1d': np.full(n, _TMAX),
        'sigma_v_1d': np.full(n, 100.0),
    }


def main():
    failures = []

    # (2) + (3) ENVELOPE MATH / AXIAL CAP on synthetic uniform overburden.
    fd = _synthetic()
    q = 2.0 * 100.0 * np.tan(np.radians(30.0))   # perimeter*sigma*tan(phi)
    d1 = fd['dist_end1_1d']
    expect = np.minimum(_TMAX, q * np.minimum(d1, 6.0 - d1))
    cap = _bond_slip_caps(fd, {'geo': (0.0, 30.0, 2.0)})
    if not np.allclose(cap, expect):
        failures.append(f"envelope math wrong: {np.round(cap,2)} vs {np.round(expect,2)}")
    else:
        print(f"[envelope] double-ended ramp reproduced, capped at t_max "
              f"({np.round(cap,1).tolist()})")
    star = _bond_slip_caps(fd, {'*': (0.0, 30.0, 2.0)})
    if not np.allclose(star, cap):
        failures.append("'*' key differs from the explicit label")
    huge = _bond_slip_caps(fd, {'geo': (0.0, 60.0, 10.0)})   # envelope >> t_max
    if not np.allclose(huge, _TMAX):
        failures.append(f"strong bond did not cap at t_max: {np.round(huge,1)}")
    else:
        print(f"[axial cap] a strong bond caps every element at t_max = {_TMAX}")

    fem = _fixture()

    def solve(**kw):
        with contextlib.redirect_stdout(io.StringIO()):
            return solve_fem(fem, F=1.0, max_iterations=_MAX_ITER, tolerance=1e-3,
                             debug_level=0, **kw)

    base = solve()
    base_t = float(np.max(base['forces_1d']))
    if not base['converged']:
        failures.append(f"fixture did not converge at F=1 ({base['iterations']} iters)")
    else:
        print(f"[fixture] converged in {base['iterations']} iters, "
              f"end-ramp max bar force = {base_t:.2f} kN")

    # (1) INVARIANCE.
    none_sol = solve(bond_slip=None)
    if _digest(none_sol) != _digest(base):
        failures.append("bond_slip=None differs from the default (end-ramp) path")
    else:
        print("[invariance] bond_slip=None is bit-identical to the default path")

    # (4) ENGAGE + CONVERGE: a tight bond starves the bar but the slope still converges.
    tight = {'geo': (0.0, 3.0, 0.2)}
    cap_t = _bond_slip_caps(fem, tight)
    cap_max = float(cap_t[fem['elem_length_1d'] > 0].max())
    eng = solve(bond_slip=tight)
    eng_t = float(np.max(eng['forces_1d']))
    if not eng['converged']:
        failures.append(f"tight-bond solve did not converge ({eng['iterations']} iters)")
    if _digest(eng) == _digest(base):
        failures.append("tight bond left the solution unchanged (bond never fired)")
    if eng_t > cap_max + 1e-6:
        failures.append(f"reported bar force {eng_t:.2f} exceeds the bond cap {cap_max:.2f}")
    if eng_t > 0.6 * base_t:
        failures.append(f"tight bond barely moved the bar force ({base_t:.2f} -> {eng_t:.2f})")
    if not [f for f in failures if 'tight' in f or 'bond cap' in f or 'bond never' in f]:
        print(f"[engage] tight bond (cap <= {cap_max:.2f}): max bar force "
              f"{base_t:.2f} -> {eng_t:.2f} (converged in {eng['iterations']} iters)")

    # (5) GRADIENT MONOTONE: progressively finer bond -> non-increasing bar force.
    forces = []
    for phi, perim in [(10.0, 0.4), (6.0, 0.3), (3.0, 0.2)]:
        s = solve(bond_slip={'geo': (0.0, phi, perim)})
        forces.append(float(np.max(s['forces_1d'])))
    if any(forces[i + 1] > forces[i] + 1e-6 for i in range(len(forces) - 1)):
        failures.append(f"bar force not monotone in bond tightness: {np.round(forces,2)}")
    else:
        print(f"[gradient] tighter bond -> less force: {np.round(forces,2).tolist()} kN")

    # (6) UNKNOWN LINE REJECTED.
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            solve_ssrm(fem, F_min=1.0, F_max=1.0, max_iterations=5, debug_level=0,
                       bond_slip={'no_such_line': (0.0, 30.0, 2.0)})
        failures.append("unknown bond_slip line name did not raise")
    except ValueError:
        print("[guard] unknown bond_slip line name rejected")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: bond-slip caps the reinforcement force gradient by the Coulomb bond "
          "envelope, reproduces the classical double-ended ramp under uniform overburden, "
          "respects the axial cap, is bit-identical to the end-ramp path when off, and "
          "rejects unknown line names.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
