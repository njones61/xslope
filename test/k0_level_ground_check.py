"""K0 initial stress must be an EXACT equilibrium on level ground.

On flat ground with the base fixed and the sides on x-rollers, the field

    sigma'_v = -gamma' z      sigma'_h = sigma'_z = K0 sigma'_v      tau_xy = 0

satisfies equilibrium under self weight identically, for ANY K0 — the vertical
equilibrium equation contains only sigma_v, and the horizontal one only the
x-derivative of sigma_h, which is zero when the ground is flat. So a correct
initial-stress implementation has nothing to redistribute here: it must converge
on the first iteration, leave the mesh undisplaced to machine precision,
reproduce the imposed stress field element by element, and yield nowhere.

That makes level ground the one configuration where the K0 answer is known in
closed form, and it is therefore the standing check on the whole path — the
overburden integration, the initial-stress load term int B^T sigma_0 dV, the
addend at the yield check, and the pore-pressure convention. Any of them wrong
by a term shows up here as a displacement that should not exist. Four
configurations are pinned (K0 = 0.5, 1.0, 2.0 dry, and K0 = 1.0 with the water
table at the surface, where the K0 relation holds between EFFECTIVE stresses and
the recovered total horizontal stress is K0 sigma'_v + u).

A fifth leg pins the SSRM sequencing: a trial that starts from the equilibrated
in-situ state must, on level ground, be handed back exactly the state it would
have built for itself. There is nothing to equilibrate here, so carrying the
state forward has to be an exact no-op — same one-iteration convergence, same
zero displacement, same stresses.

Run directly:  PYTHONPATH=. python3 test/k0_level_ground_check.py
"""

import copy
import warnings
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Polygon

warnings.filterwarnings('ignore')

from xslope.fem import build_fem_data, solve_fem, _prepare_fem_model
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons, get_material_polygons

W, H = 20.0, 10.0        # block width and height
GAMMA = 20.0             # kN/m3
GAMMA_W = 9.81
TARGET = 1.0             # mesh target size (tri6)

# Tolerances. These are machine-precision assertions, not engineering ones: the
# field is an exact equilibrium, so the only thing between the solver and zero is
# floating-point round-off on stresses of order 200 kPa.
MAX_U = 1e-12            # absolute displacement (metres) — the mesh must not move
REL_STRESS = 1e-9        # relative to the peak vertical stress
MAX_ITER_EXPECTED = 1    # equilibrium on the first iteration, nothing to redistribute

# The base file supplies the boilerplate (units, gamma_water, solver options); the
# geometry and the single material are replaced wholesale below.
BASE_FILE = Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files' / 'xslope_griffiths1.xlsx'


def _slope_data(wet):
    """A W x H rectangle of one soil, flat top, optionally with the water table
    at the ground surface."""
    d = copy.deepcopy(load_slope_data(str(BASE_FILE)))
    d['unit_system'] = 'metric'
    d['gamma_water'] = GAMMA_W
    poly = Polygon([(0.0, H), (W, H), (W, 0.0), (0.0, 0.0)])
    d['profile_lines'] = []
    d['polygons'] = [{'polygon': poly, 'mat_id': 0}]
    d['domain_polygon'] = poly
    d['ground_surface'] = LineString([(0.0, H), (W, H)])
    d['max_depth'] = 0.0
    d['circles'] = []
    d['non_circ'] = []
    d['piezo_line'] = []
    d['piezo_phreatic'] = False
    m = d['materials'][0]
    m.update(dict(name='soil', gamma=GAMMA, gamma_sat=None, option='mc',
                  c=25.0, phi=30.0, E=30000.0, nu=0.3, t_cut=None,
                  u='none', ru=0.0))
    if wet:
        d['piezo_line'] = [(0.0, H), (W, H)]
        m['u'] = 'piezo'
        m['gamma_sat'] = GAMMA
    return d


def _build(wet):
    d = _slope_data(wet)
    mesh = build_mesh_from_polygons(get_material_polygons(d), target_size=TARGET,
                                    element_type='tri6')
    return build_fem_data(d, mesh)


def _gauss_point_xy(fem_data, prep):
    """(element, Gauss point) -> (x, y), from the shape functions the solver
    integrates with — so the expected stresses are compared at the same points."""
    nodes, elements = fem_data['nodes'], fem_data['elements']
    et = fem_data['element_types']
    out = []
    for e, gps in enumerate(prep['elem_gp_data']):
        xy = nodes[elements[e][:et[e]]]
        out.append([(float(g['N'] @ xy[:, 0]), float(g['N'] @ xy[:, 1])) for g in gps])
    return out


def _expected(fem_data, prep, k0, wet):
    """Element-average TOTAL stresses (compression-positive) for the imposed field."""
    coords = _gauss_point_xy(fem_data, prep)
    y_top = float(np.max(fem_data['nodes'][:, 1]))
    n_el = len(fem_data['elements'])
    sv = np.zeros(n_el)
    sh = np.zeros(n_el)
    for e in range(n_el):
        z = np.array([y_top - y for _, y in coords[e]])
        u_w = GAMMA_W * z if wet else np.zeros_like(z)
        sv_t = GAMMA * z
        sv[e] = sv_t.mean()
        sh[e] = (k0 * (sv_t - u_w) + u_w).mean()
    return sv, sh


def _check(tag, sol, sv_exp, sh_exp, failures, expect_iterations=MAX_ITER_EXPECTED):
    st = sol['stresses']              # [sig_x, sig_y, tau_xy, sig_vm], compression +
    ref = max(float(sv_exp.max()), 1e-30)
    err_sy = float(np.abs(st[:, 1] - sv_exp).max()) / ref
    err_sx = float(np.abs(st[:, 0] - sh_exp).max()) / ref
    err_txy = float(np.abs(st[:, 2]).max()) / ref
    max_u = float(sol['max_displacement'])
    n_plastic = int(np.count_nonzero(sol['plastic_elements']))
    iters = int(sol['iterations'])

    print(f"  {tag:<28} iters={iters:<3d} max|u|={max_u:.2e}  "
          f"err(sy,sx,txy)={err_sy:.1e},{err_sx:.1e},{err_txy:.1e}  "
          f"yielded={n_plastic}")

    if not sol['converged']:
        failures.append(f"{tag}: did not converge")
    if iters > expect_iterations:
        failures.append(f"{tag}: converged in {iters} iterations, expected "
                        f"{expect_iterations} — the K0 field is being redistributed, "
                        "so it is not in equilibrium as built")
    if max_u > MAX_U:
        failures.append(f"{tag}: max|u| = {max_u:.3e} > {MAX_U:.0e} — the mesh moved "
                        "under a field that is an exact equilibrium")
    for name, err in (('sigma_y', err_sy), ('sigma_x', err_sx), ('tau_xy', err_txy)):
        if err > REL_STRESS:
            failures.append(f"{tag}: {name} off by {err:.3e} (relative) — the "
                            "recovered stress is not the imposed K0 field")
    if n_plastic:
        failures.append(f"{tag}: {n_plastic} yielded elements under an equilibrium "
                        "state well inside the Mohr-Coulomb envelope")
    return dict(tag=tag, iterations=iters, max_u=max_u, err_sy=err_sy, err_sx=err_sx,
                err_txy=err_txy, n_plastic=n_plastic)


def run():
    """Returns a list of failure strings (empty = pass)."""
    failures = []
    print(f"K0 level-ground check ({W:g} x {H:g} m block, tri6 target {TARGET:g}, "
          f"gamma = {GAMMA:g}):")

    cases = [(0.5, False), (1.0, False), (2.0, False), (1.0, True)]
    fem_dry = _build(False)
    fem_wet = _build(True)

    for k0, wet in cases:
        fem_data = fem_wet if wet else fem_dry
        prep = _prepare_fem_model(fem_data, k0=k0)
        sv_exp, sh_exp = _expected(fem_data, prep, k0, wet)
        sol = solve_fem(fem_data, F=1.0, k0=k0, max_iterations=2000,
                        fast_kernel=False, _prepared=prep)
        tag = f"K0={k0:g} {'wet' if wet else 'dry'}"
        _check(tag, sol, sv_exp, sh_exp, failures)

        if k0 == 1.0 and not wet:
            # SSRM sequencing: hand the equilibrated state back in as the initial
            # state (what solve_ssrm does for every trial). On level ground the
            # equilibration had nothing to do, so this must reproduce the same
            # solve exactly.
            sol2 = solve_fem(fem_data, F=1.0, k0=k0, max_iterations=2000,
                             fast_kernel=False, _prepared=prep,
                             _init_state=sol['_k0_state'])
            _check(f"K0={k0:g} dry, pre-equilibrated", sol2, sv_exp, sh_exp, failures)
            d_stress = float(np.abs(sol2['stresses'] - sol['stresses']).max())
            if d_stress > REL_STRESS * float(sv_exp.max()):
                failures.append(
                    f"carrying the equilibrated state changed the stresses by "
                    f"{d_stress:.3e} — pre-equilibration is not a no-op on level "
                    "ground")

    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe K0 field is carried as an exact equilibrium on level ground, "
          "before and after in-situ equilibration.")


if __name__ == '__main__':
    main()
