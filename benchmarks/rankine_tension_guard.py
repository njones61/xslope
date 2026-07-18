"""Guard for the RANKINE tension cutoff in solve_fem / solve_ssrm.

solve_fem caps the major (most-tensile) principal stress at a per-element strength
T through a second viscoplastic yield surface F_t = sigma_1 - T (associated flow,
summed with the Mohr-Coulomb shear flow at the corner by Koiter's rule). It is
driven either globally (tension_cutoff=True -> T = 0 everywhere) or per element
(tension_cap_by_elem, which solve_ssrm builds from a material-name -> T dict via
tension_cutoff_by_material). This guard asserts, on a small slope that develops
crest tension and CONVERGES at F = 1:

  1. INVARIANCE. tension_cap_by_elem=None with tension_cutoff=False is
     bit-identical to omitting the arguments (the historical no-cutoff path).

  2. ENGAGE. A finite cap T below the uncut major principal tension actually fires:
     the converged max major-principal tension drops from ~13.6 to <= ~T (the
     surface is returned to), a large reduction the mean-stress primitive could
     not produce on this compression-dominated mean field.

  3. GLOBAL T = 0. tension_cutoff=True drives the max principal tension to ~0 (no
     principal tensile stress permitted anywhere).

  4. SRF NO-OP AT F = 1. tension_srf=True divides T by the trial F; at F = 1 that
     is T/1 = T, so it is bit-identical to the fixed cap. (The SRF reduction of T
     with F>1 is exercised in the measurement campaign, not here.)

  5. UNKNOWN NAME REJECTED. solve_ssrm(tension_cutoff_by_material={bad: T}) raises
     ValueError rather than silently no-opping.

The cutoff is DAMPED (relaxes ~10% of the tension overshoot per iteration), so the
capped solves take more iterations than the uncut solve but still converge well
within budget; the iteration counts are printed so a regression in convergence is
visible.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/rankine_tension_guard.py
Exits non-zero on any failure.
"""
import hashlib
import os
import sys
import tempfile

import numpy as np
from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints)
from xslope.fem import build_fem_data, solve_fem, solve_ssrm

# A small single-material slope with a steep crest that develops principal tension
# near the crest yet still CONVERGES under gravity at F = 1 (c/phi chosen in the
# stable-but-near-failure regime: too strong -> no tension, too weak -> fails).
_C, _PHI, _E, _NU = 12.0, 18.0, 5000.0, 0.3
_GEOM = [(0.0, 0.0), (10.0, 0.0), (8.0, 6.0), (0.0, 6.0)]
_CAP = 6.0            # cap below the uncut max principal tension (~13.6)
_MAX_ITER = 4000


def _fixture():
    """Build the fem_data through the standard save/load + mesh pipeline so the
    base BCs are exactly what solve_fem consumes (deterministic mesh)."""
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
    sd['circles'] = [{'Xo': 5.0, 'Yo': 12.0, 'Depth': 0.0, 'R': 12.0}]
    tmp = os.path.join(tempfile.gettempdir(), 'rankine_tension_guard.xlsx')
    save_slope_data_to_xlsx(sd, tmp)
    sd = load_slope_data(tmp)
    cl, _, _ = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=cl)
    mesh = build_mesh_from_polygons(polys, target_size=1.0, element_type='tri3',
                                    lines=cl, point_constraints=extract_point_constraints(sd))
    return build_fem_data(sd, mesh)


def _max_principal_tension(sol):
    """Max major-principal tensile stress (tension-positive) over elements, from the
    reported (compression-positive) [sx, sy, txy, svm] element stresses. u='none'
    here, so total == effective (what the Rankine surface caps)."""
    st = np.asarray(sol['stresses'])
    sx, sy, txy = st[:, 0], st[:, 1], st[:, 2]
    R = np.sqrt(((sx - sy) / 2.0) ** 2 + txy ** 2)
    return float((R - (sx + sy) / 2.0).max())


def _digest(sol):
    def h(a):
        return hashlib.sha256(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()[:16]
    return '|'.join([str(sol['converged']), str(sol['iterations']),
                     h(sol['displacements']), h(sol['stresses']),
                     h(sol['plastic_elements'].astype(float))])


def main():
    fem = _fixture()
    n_elem = len(fem['element_materials'])
    cap = np.full(n_elem, _CAP)
    failures = []

    def solve(**kw):
        return solve_fem(fem, F=1.0, max_iterations=_MAX_ITER, tolerance=1e-3,
                         debug_level=0, **kw)

    base = solve()
    base_t = _max_principal_tension(base)
    if not base['converged']:
        failures.append(f"fixture did not converge at F=1 ({base['iterations']} iters)")
    if base_t < 2.0 * _CAP:
        failures.append(f"fixture max principal tension {base_t:.2f} too low to guard a cap of {_CAP}")
    else:
        print(f"[fixture] converged in {base['iterations']} iters, "
              f"max principal tension = {base_t:.3f} (cap under test = {_CAP})")

    # (1) INVARIANCE: None + tension_cutoff=False is bit-identical to the default.
    none_sol = solve(tension_cap_by_elem=None, tension_cutoff=False)
    if _digest(none_sol) != _digest(base):
        failures.append("tension_cap_by_elem=None / tension_cutoff=False differs from the default path")
    else:
        print("[invariance] tension_cap_by_elem=None is bit-identical to the default path")

    # (2) ENGAGE: a finite cap returns the major principal stress to ~T.
    capped = solve(tension_cap_by_elem=cap)
    capped_t = _max_principal_tension(capped)
    if not capped['converged']:
        failures.append(f"capped solve did not converge ({capped['iterations']} iters)")
    if _digest(capped) == _digest(base):
        failures.append("finite cap left the solution unchanged (cutoff never fired)")
    if capped_t > 1.1 * _CAP:
        failures.append(f"capped major principal tension {capped_t:.3f} exceeds ~T={_CAP} (not returned to surface)")
    if capped_t > 0.6 * base_t:
        failures.append(f"cap barely moved the tension ({base_t:.2f} -> {capped_t:.2f})")
    if not [f for f in failures if 'capped' in f or 'cap ' in f or 'cutoff' in f]:
        print(f"[engage] cap T={_CAP}: max principal tension {base_t:.3f} -> {capped_t:.3f} "
              f"(converged in {capped['iterations']} iters)")

    # (3) GLOBAL T = 0: no principal tension survives.
    zero = solve(tension_cutoff=True)
    zero_t = _max_principal_tension(zero)
    if not zero['converged']:
        failures.append(f"global T=0 cutoff did not converge ({zero['iterations']} iters)")
    if zero_t > 0.05 * base_t:
        failures.append(f"global T=0 cutoff left principal tension {zero_t:.3f} (should be ~0)")
    else:
        print(f"[global T=0] max principal tension -> {zero_t:.3f} "
              f"(converged in {zero['iterations']} iters)")

    # (4) SRF NO-OP AT F=1: tension_srf divides T by F; at F=1, T/1 == T exactly.
    srf1 = solve(tension_cap_by_elem=cap, tension_srf=True)
    if _digest(srf1) != _digest(capped):
        failures.append("tension_srf=True at F=1 differs from the fixed cap (T/1 should equal T)")
    else:
        print("[srf] tension_srf=True is bit-identical to the fixed cap at F=1 (T/1 = T)")

    # (5) UNKNOWN NAME REJECTED: solve_ssrm must raise, not silently no-op.
    try:
        solve_ssrm(fem, F_min=1.0, F_max=1.0, max_iterations=5, debug_level=0,
                   tension_cutoff_by_material={'NoSuchMat': 5.0})
        failures.append("unknown tension_cutoff_by_material name did not raise")
    except ValueError:
        print("[guard] unknown tension_cutoff_by_material name rejected")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: the Rankine tension cutoff returns the major principal stress to the "
          "per-element cap T (and to ~0 under the global flag), is bit-identical to the "
          "default path when off, and rejects unknown material names.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
