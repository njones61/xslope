"""The two SSRM drivers must agree on one locked case.

The strength-reduction bisection can be driven two ways (see SPIKE.md): the
viscoplastic initial-stiffness iteration, which is the default and the definition
of every locked and published factor of safety, and a Newton-Raphson iteration on
a consistent Mohr-Coulomb tangent, selected by ``fem_solver='newton'``. They are
different numerical routes to the same physical question — at what strength
reduction does the slope stop standing — so they have to reach the same answer.

This pins that on Griffiths & Lane Example 1 (docs/fem/files/xslope_griffiths1.xlsx)
at a deliberately coarse tri6 mesh, which the verification page locks at FS = 1.39.
Three things are asserted:

  * both drivers land within the bisection tolerance of the locked value;
  * they land within the bisection tolerance of EACH OTHER;
  * the Newton run leaves no trial undecided — every trial the bisection visits
    comes back converged or failed, never 'inconclusive'.

The last one is the point of the whole exercise. The viscoplastic path can stop a
trial at its iteration ceiling with the residual still falling and no verdict to
give; the Newton path decides a trial by whether the full gravity load is reached
in equilibrium, so it cannot.

The mesh is coarse on purpose — this check exists to catch a regression in either
driver, not to re-measure the benchmark, and it has to stay cheap enough to keep
in the suite.

Run directly:  PYTHONPATH=. python3 test/nr_ssrm_check.py
"""

import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

from xslope.fem import build_fem_data, solve_ssrm
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons, get_material_polygons

MODEL = (Path(__file__).resolve().parents[1] / 'docs' / 'fem' / 'files'
         / 'xslope_griffiths1.xlsx')
ELEMENT_TYPE = 'tri6'
TARGET_SIZE = 6.0
LOCKED_FS = 1.39         # docs/verification/ssrm.md, the coarse tri6 row
F_MIN, F_MAX = 1.30, 1.50
TOLERANCE = 0.05         # the bisection cell width, and the agreement tolerance
MAX_ITER = 4000


def _fem_data():
    slope_data = load_slope_data(str(MODEL))
    mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                    target_size=TARGET_SIZE,
                                    element_type=ELEMENT_TYPE)
    return build_fem_data(slope_data, mesh)


def run():
    failures = []
    fem_data = _fem_data()
    results = {}
    for solver in ('viscoplastic', 'newton'):
        results[solver] = solve_ssrm(
            fem_data, F_min=F_MIN, F_max=F_MAX, tolerance=TOLERANCE,
            max_iterations=MAX_ITER, fem_solver=solver,
            capture_failure_state=False)

    for solver, res in results.items():
        fs = res.get('FS')
        if fs is None:
            failures.append(f"{solver}: no factor of safety returned "
                            f"({res.get('message', 'no message')})")
            continue
        if abs(fs - LOCKED_FS) > TOLERANCE:
            failures.append(
                f"{solver}: FS = {fs:.4f} is {abs(fs - LOCKED_FS):.4f} from the "
                f"locked {LOCKED_FS:.2f}, outside the bisection tolerance "
                f"{TOLERANCE:g}")

    fs_vp = results['viscoplastic'].get('FS')
    fs_nr = results['newton'].get('FS')
    if fs_vp is not None and fs_nr is not None and abs(fs_vp - fs_nr) > TOLERANCE:
        failures.append(
            f"the two drivers disagree: viscoplastic {fs_vp:.4f} vs Newton "
            f"{fs_nr:.4f}, a gap of {abs(fs_vp - fs_nr):.4f} against the "
            f"bisection tolerance {TOLERANCE:g}")

    undecided = [t['F'] for t in results['newton'].get('trials', [])
                 if t.get('exit_reason') == 'inconclusive']
    if undecided:
        failures.append(
            "the Newton driver left trial(s) undecided at F = "
            + ", ".join(f"{f:.4f}" for f in undecided)
            + " — its verdict is meant to be binary")

    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nBoth SSRM drivers reproduce the locked factor of safety, agree with "
          "each other, and the Newton run decided every trial.")


if __name__ == '__main__':
    main()
