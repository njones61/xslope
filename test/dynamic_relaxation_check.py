"""The explicit dynamic-relaxation driver, and that it changes nothing switched off.

A second per-trial engine sits beside the viscoplastic sweep, reached by
``fem_solver='dynamic'``. Every node is given a mass, Newton's second law is
integrated forward by central differences and the motion is damped; a slope that
can stand comes to rest and a slope that cannot keeps moving. Nothing selects it
by default, and every locked factor of safety in the repository stays defined by
the sweep.

Four legs:

  1. **The static answer.** An unjointed elastic block on an elastic foundation,
     gravity turn-on, settled from zero displacement. The state the integrator
     reaches must BE the direct elastic solution ``K u = f``, to 1e-6 relative on
     the displacement field and on the stress field, under all three damping
     options. This is what says the mass scaling, the step and the assembled
     internal force are right before any nonlinearity is involved.
  2. **A single block's toppling threshold.** Goodman & Bray's block on a
     Coulomb joint at ``b/h = 0.400`` under a horizontal load coefficient: it
     must stand at ``k = 0.40`` and topple at ``k = 0.42``, which is where the
     sweep puts it at this mesh. The engine's failing test carries no absolute
     displacement level at all, so a threshold this sharp — 0.15 % of the weight
     separates the two sides of the sliding form — is the guard against damping
     that makes a failing slope look at rest.
  3. **The verdict rule, with the loop left to read it alone.** Every case in
     this leg runs with the Newton corrector switched OFF, because with it on a
     standing single block is certified at the first step rung -- 300 steps --
     before the verdict's own 500-step window has closed even once, and a rule
     that reads a standing slope as failing is invisible. The elastic block at
     the shipped start must read standing; the toppling block at ``k = 0.38`` and
     ``k = 0.40``, both standing by the closed form, must not be read
     ``diverging``; the sliding block must stand at ``k = 0.55`` and ``k = 0.57``
     and be read ``diverging`` at ``k = 0.58``.
  4. **Goodman direct shear.** ``test/joint_element_check.py``'s row 1 run on the
     dynamic driver at the tolerances that check already holds the sweep to: the
     elastic branch, every slipping pair exactly on the Mohr-Coulomb limit, the
     slip load against its closed form, and no traction oscillation at a hundred
     times the normal stiffness.

     The two STRENGTH-REDUCTION rows — the block on an inclined plane
     (``FS = tan phi_j / tan beta``) and the infinite-slope form — pass on the
     dynamic driver and are still not here, for cost alone: they return 1.5874
     against 1.5863 (+0.07 %) and 1.5581 against 1.5837 (-1.62 %), inside the 3 %
     the check holds, and the two bisections take 656 s and 405 s because a
     standing trial on the explicit path needs 8,565 to 20,181 steps to satisfy
     the standing test. Run them with
     ``python3 test/joint_element_check.py --driver dynamic``.

     They read -14.70 % and -7.17 % until round A2, and the cause was not the
     verdict: re-run on the settled rule they were unmoved to four figures. It
     was the bisection's budget, which is the sweep's 3,000 iterations scaled
     with ``k_s``, and on this driver 3,000 steps is less work than a STANDING
     trial needs — so every trial above 1.3530 ended its budget, reported FAILED,
     and walked the bracket down. ``joint_element_check._ssrm`` now gives the
     dynamic driver a flat budget.
  5. **Switched off, nothing moved.** With the driver unselected, a jointed
     sweep solve is byte-identical to the same solve on a pristine package built
     from ``git show joints-fix:``. The comparison is on the raw bytes of the
     displacement, stress and interface fields, not on a tolerance.

Run directly:  PYTHONPATH=. python3 test/dynamic_relaxation_check.py
"""

import contextlib
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np

from xslope.fem import build_fem_data, solve_fem, resolve_fem_solver
from xslope.mesh import (build_mesh_from_polygons, extract_constraint_line_geometry,
                         extract_joint_options, extract_point_constraints,
                         extract_size_regions, get_material_polygons)

#: The branch the switched-off comparison is made against.
BASELINE_REF = os.environ.get('XSLOPE_DR_BASELINE') or 'joints-fix'


# --------------------------------------------------------------------------
# Leg 1 — the static answer
# --------------------------------------------------------------------------

def _elastic_block(ts=1.0):
    """A 4 x 3 elastic block welded onto a 12 x 4 elastic foundation. No joints,
    no plasticity: the only equilibrium is the elastic one."""
    sys.path.insert(0, os.path.join(_ROOT, 'benchmarks', 'rocscience'))
    import build_joint_problems as B
    block = [(4.0, 0.0), (8.0, 0.0), (8.0, 3.0), (4.0, 3.0)]
    found = [(0.0, -4.0), (12.0, -4.0), (12.0, 0.0), (8.0, 0.0), (4.0, 0.0), (0.0, 0.0)]
    sd = B._base()
    sd['side_bc'] = 'rollers'
    mats = [B._rock('Block', 25.0, 1.0e6, 0.25, 0.0, 0.0, option='elastic'),
            B._rock('Foundation', 22.0, 5.0e5, 0.30, 0.0, 0.0, option='elastic')]
    B._finish(sd, [(block, 0), (found, 1)], mats)
    sd['non_circ'] = B._surface([(4.0, 0.0), (8.0, 0.0)])
    lines, _a, _b = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=ts, element_type='tri6', lines=lines,
            element_size_1d=None, point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd),
            joint_lines=extract_joint_options(sd))
        fd = build_fem_data(sd, mesh)
    assert fd.get('joint_data') is None
    return fd


def _leg_static(failures, results):
    """The integrator's fixed point IS the direct elastic solution."""
    fd = _elastic_block(1.0)
    with contextlib.redirect_stdout(io.StringIO()):
        ref = solve_fem(fd, F=1.0, max_iterations=20000, fast_kernel=False,
                        fem_solver='viscoplastic')
    s_ref = np.asarray(ref['stresses'], dtype=float)
    for damping in ('local', 'kinetic', 'viscous'):
        with contextlib.redirect_stdout(io.StringIO()):
            dyn = solve_fem(fd, F=1.0, max_iterations=60000,
                            max_iterations_ceiling=120000, force_tol=1e-7,
                            fast_kernel=False, fem_solver='dynamic',
                            _corrector=False, dr_damping=damping, dr_start='zero')
        u_dyn = np.asarray(dyn['displacements'], dtype=float)
        u_dir = np.asarray(dyn['displacements_elastic'], dtype=float)
        den = float(np.max(np.abs(u_dir)))
        err_u = float(np.max(np.abs(u_dyn - u_dir))) / den if den > 0 else np.inf
        s_dyn = np.asarray(dyn['stresses'], dtype=float)
        sden = float(np.max(np.abs(s_ref)))
        err_s = float(np.max(np.abs(s_dyn - s_ref))) / sden if sden > 0 else np.inf
        if not dyn['converged']:
            failures.append(f"leg 1 ({damping}): the elastic block does not "
                            f"report standing (exit {dyn.get('exit_reason')})")
        if not (err_u <= 1e-6):
            failures.append(f"leg 1 ({damping}): the settled displacement field "
                            f"is {err_u:.2e} from the direct elastic solution, "
                            f"past 1e-6")
        if not (err_s <= 1e-6):
            failures.append(f"leg 1 ({damping}): the settled stress field is "
                            f"{err_s:.2e} from the direct one, past 1e-6")
        results.append(f"static   {damping:<8} {dyn.get('dr_steps'):>6} steps, "
                       f"displacement {err_u:.1e}, stress {err_s:.1e} from the "
                       f"direct elastic answer")


# --------------------------------------------------------------------------
# Leg 2 — the single block's toppling threshold
# --------------------------------------------------------------------------

def _block_on_joint(k, b=1.0, h=2.5, phi_j=45.0, ts=0.4):
    """Goodman & Bray's block: an elastic rectangle on a Coulomb joint over a
    stiff plate, with gravity turned through atan(k) by the seismic coefficient.
    The block's faces are free (the ground surface steps over it) and the joint's
    ends sit on the external boundary, so no corner is a rigid pin."""
    sys.path.insert(0, os.path.join(_ROOT, 'benchmarks', 'rocscience'))
    import build_joint_problems as B
    L, D = 2.0, 1.0
    A, Bp = (0.0, 0.0), (b, 0.0)
    xR = b + L
    block = [A, Bp, (b, h), (0.0, h)]
    plate = [(-L, 0.0), (-L, -D), (xR, -D), (xR, 0.0), Bp, A]
    sd = B._base()
    sd['side_bc'] = 'rollers'
    sd['k_seismic'] = float(k)
    mats = [B._rock('Block', 25.0, 1.0e6, 0.25, 0.0, 0.0, option='elastic'),
            B._rock('Plate', 25.0, 1.0e7, 0.25, 0.0, 0.0, option='elastic')]
    B._finish(sd, [(block, 0), (plate, 1)], mats)
    sd['joint_lines'] = [B._joint('base', A, Bp, 0.0, phi_j, kn=1.0e8, ks=1.0e7,
                                  t_cut=0.0)]
    sd['non_circ'] = B._surface([A, Bp])
    lines, _a, _b = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=ts, element_type='tri6', lines=lines,
            element_size_1d=None, point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd),
            joint_lines=extract_joint_options(sd))
        fd = build_fem_data(sd, mesh)
    return fd


def _leg_threshold(failures, results):
    """b/h = 0.400: standing at k = 0.40, toppling at k = 0.42."""
    for k, want in ((0.40, True), (0.42, False)):
        fd = _block_on_joint(k)
        with contextlib.redirect_stdout(io.StringIO()):
            sol = solve_fem(fd, F=1.0, max_iterations=40000,
                            max_iterations_ceiling=80000, tension_srf=False,
                            k0=None, fast_kernel=False, fem_solver='dynamic')
        stands = bool(sol.get('stable', sol['converged']))
        if stands != want:
            failures.append(
                f"leg 2: the block at k = {k:.2f} (b/h = 0.400) reads "
                f"{'standing' if stands else 'failing'} on the dynamic driver "
                f"and the closed form says {'standing' if want else 'toppling'} "
                f"(exit {sol.get('exit_reason')}, {sol.get('iterations')} steps)")
        results.append(f"block    k = {k:.2f}: "
                       f"{'stands' if stands else 'topples'} "
                       f"({sol.get('exit_reason')}, {sol.get('iterations')} steps, "
                       f"max|u| {sol['max_displacement']:.3e})")


# --------------------------------------------------------------------------
# Leg 3 — the verdict rule, with the loop left to read it alone
# --------------------------------------------------------------------------

def _leg_verdict(failures, results):
    """The cases round A measured the verdict wrong on, corrector OFF.

    The corrector is switched off on every case here and that is the whole point:
    with it on, every standing single block is certified at the first step rung,
    300 steps, before the verdict's own 500-step window has closed even once, so
    a rule that reads a standing slope as failing is invisible. Round A's
    measurement of that is section 7.1 of its report; these are the cases it
    named, run against the rule that replaced the one it refuted.

    Four readings, and each is a regression guard on a defect that was measured:

      * the elastic block at the SHIPPED start reads standing. Started at its own
        elastic answer it never moves, so its peak kinetic energy is round-off and
        the kinetic floor -- a ratio to that peak -- could never be met: the block
        spent a 40,000-step budget reading `dr_undecided` with a residual of
        1e-12;
      * the toppling block at k = 0.38 and k = 0.40, both standing by the closed
        form and by the sweep, must not be read `diverging`. The design's rule
        read them failing at 950 and 1,150 steps;
      * the sliding block at k = 0.55 and k = 0.57 stands, read by the loop alone
        rather than by the corrector;
      * the sliding block at k = 0.58 -- a real runaway, its kinetic energy and
        its displacement both growing quadratically -- is still read `diverging`,
        and the step it costs is recorded, because this is the case the sweep
        never decided at all.
    """
    # the elastic block at the shipped start: standing, by the loop alone
    fd = _elastic_block(1.0)
    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_fem(fd, F=1.0, max_iterations=5000,
                        max_iterations_ceiling=5000, fast_kernel=False,
                        fem_solver='dynamic', _corrector=False)
    if not sol['converged']:
        failures.append(
            f"leg 3: the elastic block started at its own elastic answer does "
            f"not read standing with the corrector off (exit "
            f"{sol.get('exit_reason')}, {sol.get('dr_steps')} steps, "
            f"out-of-balance {sol.get('unbalanced_force_ratio'):.1e})")
    results.append(f"verdict  elastic block, shipped start: "
                   f"{sol.get('exit_reason')} at {sol.get('dr_steps')} steps")

    # the two standing toppling blocks round A read as failing
    for k in (0.38, 0.40):
        fd = _block_on_joint(k, ts=0.2)
        with contextlib.redirect_stdout(io.StringIO()):
            sol = solve_fem(fd, F=1.0, max_iterations=5000,
                            max_iterations_ceiling=5000, tension_srf=False,
                            k0=None, fast_kernel=False, fem_solver='dynamic',
                            _corrector=False)
        if sol.get('exit_reason') == 'diverging':
            failures.append(
                f"leg 3: the toppling block at k = {k:.2f} (b/h = 0.400) stands "
                f"by the closed form and the loop read it `diverging` at "
                f"{sol.get('dr_steps')} steps with the corrector off")
        results.append(f"verdict  topple k = {k:.2f}, corrector off: "
                       f"{sol.get('exit_reason')} at {sol.get('dr_steps')} steps")

    # the sliding threshold, both sides, read by the loop alone
    for k, want in ((0.55, True), (0.57, True), (0.58, False)):
        fd = _block_on_joint(k, b=1.5, h=1.0, phi_j=30.0, ts=0.2)
        with contextlib.redirect_stdout(io.StringIO()):
            sol = solve_fem(fd, F=1.0, max_iterations=20000,
                            max_iterations_ceiling=20000, tension_srf=False,
                            k0=None, fast_kernel=False, fem_solver='dynamic',
                            _corrector=False)
        exit_reason = sol.get('exit_reason')
        stands = bool(sol['converged'])
        if want and not stands:
            failures.append(
                f"leg 3: the sliding block at k = {k:.2f} (tan phi = 0.5774) "
                f"stands and the loop read it {exit_reason} at "
                f"{sol.get('dr_steps')} steps with the corrector off")
        if not want and exit_reason != 'diverging':
            failures.append(
                f"leg 3: the sliding block at k = {k:.2f} runs away -- its "
                f"kinetic energy and its displacement both grow quadratically -- "
                f"and the loop read it {exit_reason} at {sol.get('dr_steps')} "
                f"steps rather than `diverging`")
        results.append(f"verdict  slide  k = {k:.2f}, corrector off: "
                       f"{exit_reason} at {sol.get('dr_steps')} steps")


# --------------------------------------------------------------------------
# Leg 4 — two closed forms, on the element check's own rows
# --------------------------------------------------------------------------

def _leg_closed_forms(failures, results):
    """Goodman direct shear on the dynamic driver, at the element check's own
    tolerances. See this module's docstring for why the two strength-reduction
    rows are not here yet."""
    import joint_element_check as J
    saved = J.DRIVER
    J.DRIVER = 'dynamic'
    try:
        f1, r1 = [], []
        J._leg_goodman(f1, r1)
    finally:
        J.DRIVER = saved
    for f in f1:
        failures.append("leg 4 (dynamic): " + f)
    for r in r1:
        results.append("closed   " + r)


# --------------------------------------------------------------------------
# Leg 5 — switched off, nothing moved
# --------------------------------------------------------------------------

_PROBE = r'''
import contextlib, hashlib, io, json, os, sys
import numpy as np
root = sys.argv[1]
sys.path.insert(0, root)
sys.path.insert(0, os.path.join(root, 'test'))
import xslope
assert xslope.__file__.startswith(root), xslope.__file__
from xslope.fem import solve_fem
import joint_element_check as J

def h(a):
    a = np.ascontiguousarray(np.asarray(a, dtype=float))
    return hashlib.sha256(a.tobytes()).hexdigest()[:32]

out = {'xslope': xslope.__file__}
_d, fd, _g = J.slab_model(**dict(J.ROW2, k_seismic=0.0))
for driver in ('viscoplastic', None):
    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_fem(fd, F=1.0, max_iterations=4000,
                        max_iterations_ceiling=8000, fast_kernel=False,
                        fem_solver=driver)
    key = driver or 'auto'
    out[key] = {
        'converged': bool(sol['converged']),
        'verdict': str(sol.get('verdict')),
        'exit_reason': str(sol.get('exit_reason')),
        'iterations': int(sol['iterations']),
        'oob': repr(float(sol.get('unbalanced_force_ratio', 0.0))),
        'u': h(sol['displacements']),
        'u_elastic': h(sol['displacements_elastic']),
        'stresses': h(sol['stresses']),
        'strains': h(sol['strains']),
        'vp_shear_strain': h(sol['vp_shear_strain']),
        'joint_tn': h(sol['joint_tn']),
        'joint_ts': h(sol['joint_ts']),
        'joint_slip': h(sol['joint_slip']),
        'joint_open': h(np.asarray(sol['joint_open'], dtype=float)),
    }
print(json.dumps(out))
'''


def _run_probe(root):
    with tempfile.NamedTemporaryFile('w', suffix='.py', delete=False) as fh:
        fh.write(_PROBE)
        path = fh.name
    try:
        env = dict(os.environ)
        for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env.setdefault(var, '2')
        env.pop('PYTHONPATH', None)
        env.pop('XSLOPE_FEM_SOLVER', None)
        out = subprocess.run([sys.executable, path, root], capture_output=True,
                             text=True, env=env, cwd=root, timeout=1800)
        if out.returncode != 0:
            return None, (out.stderr or out.stdout)[-800:]
        return json.loads(out.stdout.strip().splitlines()[-1]), None
    finally:
        os.unlink(path)


def _leg_switched_off(failures, results):
    """The sweep path, byte for byte, against a pristine `joints-fix` package."""
    try:
        subprocess.run(['git', '-C', _ROOT, 'rev-parse', '--verify', BASELINE_REF],
                       capture_output=True, check=True)
    except Exception:
        results.append(f"identity  SKIPPED: {BASELINE_REF} is not a revision in "
                       f"this repository")
        return
    base = tempfile.mkdtemp(prefix='xslope-dr-baseline-')
    try:
        tar = subprocess.run(['git', '-C', _ROOT, 'archive', BASELINE_REF],
                             capture_output=True, check=True)
        subprocess.run(['tar', '-x', '-C', base], input=tar.stdout, check=True)
        mine, err_m = _run_probe(_ROOT)
        theirs, err_t = _run_probe(base)
    finally:
        shutil.rmtree(base, ignore_errors=True)
    if mine is None or theirs is None:
        failures.append(f"leg 4: the identity probe did not run "
                        f"({err_m or ''} {err_t or ''})".strip())
        return
    for key in ('viscoplastic', 'auto'):
        a, b = mine.get(key), theirs.get(key)
        if a != b:
            diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
            failures.append(
                f"leg 4: the '{key}' path is NOT byte-identical to "
                f"{BASELINE_REF} — {', '.join(diff)} differ "
                f"({ {k: (a.get(k), b.get(k)) for k in diff} })")
        else:
            results.append(f"identity  the '{key}' path is byte-identical to "
                           f"{BASELINE_REF} on all fourteen fields")


# --------------------------------------------------------------------------

def run():
    """Returns a list of failure strings (empty = pass)."""
    failures, results = [], []
    t0 = time.time()
    if resolve_fem_solver('dynamic') != 'dynamic':
        failures.append("resolve_fem_solver does not accept 'dynamic'")
    _leg_static(failures, results)
    _leg_threshold(failures, results)
    _leg_verdict(failures, results)
    _leg_closed_forms(failures, results)
    _leg_switched_off(failures, results)
    print(f"Dynamic relaxation check ({time.time() - t0:.0f} s):")
    for line in results:
        print("  " + line)
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe explicit dynamic driver reaches the static elastic answer, puts "
          "Goodman & Bray's block on the right side of its toppling threshold, "
          "reproduces Goodman direct shear, and changes nothing on the sweep "
          "path when it is not selected.")


if __name__ == '__main__':
    main()
