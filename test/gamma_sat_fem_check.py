"""A saturated unit weight must weigh the same as the zoning that spells it out.

Soil below the water table weighs gamma_sat, above it gamma. A model can say that
two ways:

* **sidecar** — one material carrying both weights, and one water table for the
  engine to split the soil at;
* **zoned** — the same slope cut into two material zones at the water table, the
  upper one weighing gamma and the lower one gamma_sat, with no gamma_sat anywhere.

They are the same soil described twice, so the finite element engine has to give
one answer for them: the same gravity load, the same overburden, the same
factor of safety, the same displaced mesh. That is what this checks, and it is
the one comparison that can tell a per-Gauss-point split from an element-constant
one — the sidecar model reaches its answer by weighing Gauss points, the zoned
model by weighing whole elements, and only a split taken at the water table
itself makes the two agree.

The two models are meshed ONCE and share the mesh element for element, so nothing
in the comparison is discretization: the zone boundary and the water table are
the same line, and the only difference between the runs is which of the two ways
the weights are declared. Agreement is therefore exact, not approximate, and the
checks below are written that way.

Four legs:

1. **The gravity load vector.** The load each formulation builds, node for node.
   This is the split itself, before any solving.
2. **The strength reduction run.** Factor of safety and displacement field.
3. **K0 initial stress.** The overburden integral is a second, independent path
   through the same two weights: a column standing in the sidecar model is split
   at the water table, in the zoned model it crosses two zones. Both the initial
   stress field and the solved displacements must match.
4. **The pore-pressure ratio.** u = ru * sigma_v weighs the soil column too, so
   the ru option is checked on the same pair.

A fifth leg proves the comparison is not vacuous: weighing the sidecar model
moist throughout — gamma_sat dropped, everything else identical — must change the
answer. Without it, an engine that ignored gamma_sat entirely would pass legs 1-4
by weighing both models 18 kN/m3 and never be caught.

Run directly:  PYTHONPATH=. python3 test/gamma_sat_fem_check.py
"""

import copy
import warnings
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Polygon

warnings.filterwarnings('ignore')

from xslope.fem import build_fem_data, solve_fem, solve_ssrm, _prepare_fem_model
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons

#: Slope geometry. A 10 m high 1:1 slope on a 5 m foundation, and a horizontal
#: water table at Y_W, which crosses the slope face at x = 16.
Y_W = 4.0
GROUND = [(0.0, 10.0), (10.0, 10.0), (20.0, 0.0), (40.0, 0.0)]
UPPER = [(0.0, 10.0), (10.0, 10.0), (16.0, Y_W), (0.0, Y_W)]
LOWER = [(0.0, Y_W), (16.0, Y_W), (20.0, 0.0), (40.0, 0.0),
         (40.0, -5.0), (0.0, -5.0)]

GAMMA = 18.0             # kN/m3, above the water table
GAMMA_SAT = 21.0         # kN/m3, below it
GAMMA_W = 9.81
C, PHI = 15.0, 25.0
E, NU = 30000.0, 0.3
TARGET = 3.0             # mesh target size (tri6)
RU = 0.15                # for the pore-pressure-ratio leg

#: The two formulations agree exactly: same mesh, same weights, same arithmetic.
#: These are the round-off allowances on quantities that are expected to be
#: identical, not engineering tolerances.
F_TOL = 1e-12            # factor of safety
U_TOL = 1e-12            # displacement, metres
LOAD_TOL = 1e-9          # gravity load / initial stress, kN (relative to the peak)

BASE_FILE = Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files' / 'xslope_griffiths1.xlsx'


def _material(name, gamma, gamma_sat, u, ru):
    return dict(name=name, gamma=gamma, gamma_sat=gamma_sat, option='mc',
                c=C, phi=PHI, cp=0.0, r_elev=0.0, d=0.0, psi=0.0,
                E=E, nu=NU, t_cut=None, u=u, ru=ru,
                pow_a=0.0, pow_b=0.0, pow_c=0.0, pow_d=0.0,
                hb_sci=0.0, hb_gsi=0.0, hb_mi=0.0, hb_d=0.0)


def _slope_data(zoned, u='none', ru=0.0, moist_only=False):
    """The slope, described either way.

    ``zoned`` builds two materials split at the water table; otherwise one
    material carries gamma and gamma_sat. ``moist_only`` drops gamma_sat from the
    sidecar model, which is the control for leg 5.
    """
    d = copy.deepcopy(load_slope_data(str(BASE_FILE)))
    d['unit_system'] = 'metric'
    d['gamma_water'] = GAMMA_W
    d['profile_lines'] = []
    d['ground_surface'] = LineString(GROUND)
    d['max_depth'] = 0.0
    d['circles'] = []
    d['non_circ'] = []
    d['dloads'] = []
    d['dloads2'] = []
    # Both models carry the same piezometric line, so the water table is the same
    # sheet in both; in the zoned model nothing reads it, which is the point.
    d['piezo_line'] = [(0.0, Y_W), (40.0, Y_W)]
    d['piezo_phreatic'] = False
    d['water_loads'] = 'manual'
    if zoned:
        d['materials'] = [_material('upper', GAMMA, None, u, ru),
                          _material('lower', GAMMA_SAT, None, u, ru)]
        mat_ids = (0, 1)
    else:
        gsat = None if moist_only else GAMMA_SAT
        d['materials'] = [_material('soil', GAMMA, gsat, u, ru)]
        mat_ids = (0, 0)
    d['polygons'] = [{'polygon': Polygon(UPPER), 'mat_id': mat_ids[0]},
                     {'polygon': Polygon(LOWER), 'mat_id': mat_ids[1]}]
    d['domain_polygon'] = Polygon(UPPER).union(Polygon(LOWER))
    return d


def _mesh():
    """One mesh, built from the zone geometry, shared by both formulations.

    Element materials come back 1 and 2 for the upper and lower zone; the sidecar
    model reads the same elements as one material.
    """
    polys = [{'coords': list(Polygon(UPPER).exterior.coords), 'mat_id': 0},
             {'coords': list(Polygon(LOWER).exterior.coords), 'mat_id': 1}]
    return build_mesh_from_polygons(polys, target_size=TARGET, element_type='tri6')


def _build(mesh, zoned, u='none', ru=0.0, moist_only=False):
    m = dict(mesh)
    if not zoned:
        m['element_materials'] = np.ones_like(np.asarray(mesh['element_materials']))
    return build_fem_data(_slope_data(zoned, u=u, ru=ru, moist_only=moist_only), m)


def _rel(a, b):
    """Largest absolute difference, and the same relative to the peak magnitude."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    d = float(np.abs(a - b).max()) if a.size else 0.0
    scale = max(float(np.abs(a).max()) if a.size else 0.0, 1e-30)
    return d, d / scale


def _check_gravity(prep_s, prep_z, failures, tag):
    d, rel = _rel(prep_s['F_gravity'], prep_z['F_gravity'])
    print(f"  {tag}: gravity load max |diff| = {d:.3e} ({rel:.2e} of peak)")
    if rel > LOAD_TOL:
        failures.append(f"{tag}: the sidecar and zoned gravity loads differ by "
                        f"{d:.3e} ({rel:.2e} of peak) — the weight split is not "
                        f"taken at the water table")


def _check_solution(sol_s, sol_z, failures, tag):
    if not (sol_s.get('converged', True) and sol_z.get('converged', True)):
        failures.append(f"{tag}: a solve did not converge, so there is nothing to "
                        f"compare")
        return
    du, _ = _rel(sol_s['displacements'], sol_z['displacements'])
    it_s, it_z = sol_s.get('iterations'), sol_z.get('iterations')
    print(f"  {tag}: displacement max |diff| = {du:.3e} m, "
          f"iterations {it_s} vs {it_z}")
    if du > U_TOL:
        failures.append(f"{tag}: the two formulations settle to different "
                        f"displacement fields (max |diff| = {du:.3e} m)")
    if it_s != it_z:
        failures.append(f"{tag}: the two formulations took different iteration "
                        f"counts ({it_s} vs {it_z}) — same problem, same path")


def run():
    """Returns a list of failure strings (empty = pass)."""
    failures = []
    mesh = _mesh()
    n_el = len(mesh['elements'])
    print(f"gamma_sat FEM pair (1:1 slope, water table at y = {Y_W:g}, "
          f"gamma {GAMMA:g} / gamma_sat {GAMMA_SAT:g}, tri6 x {n_el}):")

    fem_s = _build(mesh, zoned=False)
    fem_z = _build(mesh, zoned=True)

    # --- leg 1: the gravity load vector ---------------------------------------
    prep_s = _prepare_fem_model(fem_s)
    prep_z = _prepare_fem_model(fem_z)
    _check_gravity(prep_s, prep_z, failures, "self weight")

    # --- leg 2: the strength reduction run ------------------------------------
    ssrm_s = solve_ssrm(fem_s, F_min=0.8, F_max=2.5, tolerance=0.005, debug_level=0)
    ssrm_z = solve_ssrm(fem_z, F_min=0.8, F_max=2.5, tolerance=0.005, debug_level=0)
    fs_s = ssrm_s.get('FS') if ssrm_s.get('converged') else None
    fs_z = ssrm_z.get('FS') if ssrm_z.get('converged') else None
    print(f"  SSRM: sidecar FS = {fs_s}, zoned FS = {fs_z}")
    if fs_s is None or fs_z is None:
        failures.append("SSRM: one of the two formulations did not bracket a "
                        "factor of safety")
    elif abs(fs_s - fs_z) > F_TOL:
        failures.append(f"SSRM: sidecar FS = {fs_s:.9f} but zoned FS = {fs_z:.9f} "
                        f"(differ by {abs(fs_s - fs_z):.3e})")

    sol_s = solve_fem(fem_s, F=1.0, fast_kernel=False)
    sol_z = solve_fem(fem_z, F=1.0, fast_kernel=False)
    _check_solution(sol_s, sol_z, failures, "solve at F = 1")

    # --- leg 3: K0 initial stress ---------------------------------------------
    k0 = 1.0
    prep_s0 = _prepare_fem_model(fem_s, k0=k0)
    prep_z0 = _prepare_fem_model(fem_z, k0=k0)
    sv_s = np.array([v for gl in prep_s0['sv0_gp'] for v in gl])
    sv_z = np.array([v for gl in prep_z0['sv0_gp'] for v in gl])
    d, rel = _rel(sv_s, sv_z)
    print(f"  K0 = {k0:g}: overburden max |diff| = {d:.3e} kPa ({rel:.2e} of peak)")
    if rel > LOAD_TOL:
        failures.append(f"K0: the overburden integral differs by {d:.3e} kPa "
                        f"({rel:.2e} of peak) — the column is not split at the "
                        f"water table")
    sol_s0 = solve_fem(fem_s, F=1.0, k0=k0, fast_kernel=False, _prepared=prep_s0)
    sol_z0 = solve_fem(fem_z, F=1.0, k0=k0, fast_kernel=False, _prepared=prep_z0)
    _check_solution(sol_s0, sol_z0, failures, f"solve at F = 1, K0 = {k0:g}")

    # --- leg 4: the pore-pressure ratio ---------------------------------------
    fem_s_ru = _build(mesh, zoned=False, u='ru', ru=RU)
    fem_z_ru = _build(mesh, zoned=True, u='ru', ru=RU)
    d, rel = _rel(fem_s_ru['sigma_v'], fem_z_ru['sigma_v'])
    print(f"  ru = {RU:g}: column stress max |diff| = {d:.3e} kPa "
          f"({rel:.2e} of peak)")
    if rel > LOAD_TOL:
        failures.append(f"ru: the soil column stress differs by {d:.3e} kPa "
                        f"({rel:.2e} of peak) — u = ru * sigma_v is weighed from "
                        f"a column the two formulations do not agree on")
    _check_solution(solve_fem(fem_s_ru, F=1.0, fast_kernel=False),
                    solve_fem(fem_z_ru, F=1.0, fast_kernel=False),
                    failures, f"solve at F = 1, u = ru")

    # --- leg 5: the split is not a no-op --------------------------------------
    fem_dry = _build(mesh, zoned=False, moist_only=True)
    prep_dry = _prepare_fem_model(fem_dry)
    d, rel = _rel(prep_s['F_gravity'], prep_dry['F_gravity'])
    print(f"  control: dropping gamma_sat moves the gravity load by {d:.3e} kN "
          f"({rel:.2%} of peak)")
    if rel < 1e-3:
        failures.append("control: weighing the model moist throughout gives the "
                        "same gravity load as the split — the comparison above "
                        "proves nothing")
    ssrm_dry = solve_ssrm(fem_dry, F_min=0.8, F_max=2.5, tolerance=0.005,
                          debug_level=0)
    fs_dry = ssrm_dry.get('FS') if ssrm_dry.get('converged') else None
    print(f"  control: moist-throughout FS = {fs_dry}")
    if fs_s is not None and fs_dry is not None and abs(fs_s - fs_dry) <= F_TOL:
        failures.append(f"control: the moist-throughout model gives the same "
                        f"factor of safety ({fs_dry}) as the split one")

    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe sidecar and zoned formulations are one model to the FEM: same "
          "gravity load, same overburden, same factor of safety, same mesh "
          "displacement.")


if __name__ == '__main__':
    main()
