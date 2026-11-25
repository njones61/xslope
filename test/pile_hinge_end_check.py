# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A pile's moment at failure is read from the element end that carries the hinge.

Two beam elements meet at every interior node of a pile, and each states its own
end moment there. The beam law lets exactly one of the two release when the
moment reaches the capacity; that end's moment is the section's. The other end
reaches the capacity only through equilibrium at the node. In a converged field
the two agree. The at-failure field is not converged -- it is the iterate the
capture solve was stopped at -- and on Tutorial FEM-4's tip-fixed piles the 1D
Details panel, which drew both ends, showed moments up to 2,684 times the
capacity (its list rows: ``lower row 142282%``, ``upper row 268367%``) while
every hinged end sat on the capacity.

Shear and axial force are read once per element (at its middle), so they carry
no such pair and are not affected.

  A. **at failure, no section passes its capacity.** FEM-4's pile model with both
     tips fixed, run as the page runs it (tri6 at 2 ft, F 1.00-2.00): on both
     rows, every value of the at-failure moment profile and the list's
     utilization stand at or under the capacity, and moment hinges have formed
     (so the case is the one this check is for);
  B. **converged readings do not move.** On that run, and on the two pile runs
     shipped as sidecars (the FEM piles sample and the GeoStudio sheet-pile
     wall), every value of the converged moment profile equals the reading of
     both element ends to a part in 10^9.

Part A solves one strength-reduction run, about a minute.

Run directly:  PYTHONPATH=. python3 test/pile_hinge_end_check.py
"""

import contextlib
import io
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent
PILES = ROOT / 'docs' / 'lem' / 'files' / 'xslope_piles.xlsx'
SHIPPED = (ROOT / 'docs' / 'fem' / 'files' / 'xslope_piles_fem.xlsx',
           ROOT / 'docs' / 'verification' / 'files' / 'geostudio' / 'gs2_wall.xlsx')
#: The capacity a moment may stand above, relatively, and still be "on" it.
ON_CAP = 1e-9


def _both_ends(fem_details, fem_data, solution, i, sd, state, fsol):
    """The profile as it reads with BOTH element ends at every node (the reading
    before the hinge-end rule), for comparison."""
    real = getattr(fem_details, '_hinge_end_moments', None)
    if real is not None:
        fem_details._hinge_end_moments = lambda _fd, _sel, _n, m: m
    try:
        return fem_details.pile_profile(fem_data, solution, i, sd,
                                        field_state=state, failure_solution=fsol)
    finally:
        if real is not None:
            fem_details._hinge_end_moments = real


def _tips_fixed_run():
    """FEM-4's piles with both tips fixed, meshed and solved as the page does."""
    from xslope.fem import build_fem_data, solve_ssrm
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_point_constraints, extract_size_regions,
                             get_material_polygons)

    sd = load_slope_data(str(PILES))
    for row in sd['pile_lines']:
        row['tip_fixity'] = 'fixed'
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            get_material_polygons(sd, reinf_lines=lines), 2.0, 'tri6',
            lines=lines or None,
            point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd))
        fem_data = build_fem_data(sd, mesh)
        result = solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.01,
                            debug_level=0, capture_failure_state=True,
                            failure_criterion='non_convergence',
                            max_iterations=12000)
    solution = dict(result['last_solution'])
    solution['failure_solution'] = result.get('failure_solution')
    return sd, fem_data, solution


def _converged_unchanged(fails, tag, fem_details, fem_data, solution, sd, fsol):
    for i in range(len(sd.get('pile_lines') or [])):
        new = fem_details.pile_profile(fem_data, solution, i, sd,
                                       field_state='converged',
                                       failure_solution=fsol)
        old = _both_ends(fem_details, fem_data, solution, i, sd, 'converged', fsol)
        a, b = np.asarray(old['moment'], float), np.asarray(new['moment'], float)
        scale = max(1.0, float(np.max(np.abs(a))) if a.size else 1.0)
        if a.shape != b.shape or (a.size and float(np.max(np.abs(a - b))) > 1e-9 * scale):
            fails.append(f"B: {tag}, {new['label']}: the converged moment profile "
                         f"moved by {float(np.max(np.abs(a - b))) if a.shape == b.shape else 'shape'}"
                         f" of {scale:.6g}")


def check_tips_fixed(fails):
    """A, and B on the same run."""
    from xslope import fem_details

    sd, fem_data, solution = _tips_fixed_run()
    fsol = solution.get('failure_solution')
    if not fsol:
        fails.append("A: the tip-fixed run captured no at-failure state")
        return
    if not np.asarray(fsol.get('yielded_pile_M', []), bool).any():
        fails.append("A: no moment hinge formed at failure; the case does not "
                     "exercise the hinge-end reading")
    caps = np.asarray(fem_data['M_cap_by_pile_elem'], float)
    for i, row in enumerate(sd['pile_lines']):
        prof = fem_details.pile_profile(fem_data, solution, i, sd,
                                        field_state='failure', failure_solution=fsol)
        cap = float(prof['M_cap']) if prof.get('M_cap') else float(np.min(caps))
        peak = float(np.max(np.abs(prof['moment'])))
        if peak > cap * (1.0 + ON_CAP):
            fails.append(f"A: {prof['label']} at failure reads a moment of "
                         f"{peak:,.0f} against a capacity of {cap:,.0f} "
                         f"({peak / cap:.0%})")
    for entry in fem_details.list_lines(fem_data, solution, sd,
                                        field_state='failure', failure_solution=fsol):
        util = entry.get('utilization')
        if entry['kind'] == 'pile' and util is not None and util > 1.0 + ON_CAP:
            fails.append(f"A: the 1D Details list reads {entry['label']} at "
                         f"{util:.0%} at failure")
    _converged_unchanged(fails, "FEM-4 tips fixed", fem_details, fem_data,
                         solution, sd, fsol)


def check_shipped_runs(fails):
    """B on the pile runs shipped as sidecars."""
    from xslope import fem_details
    from xslope.fileio import load_slope_data
    from xslope.report import solutions_from_sidecars

    for path in SHIPPED:
        if not path.exists():
            fails.append(f"B: missing {path}")
            continue
        sd = load_slope_data(str(path))
        bundle = solutions_from_sidecars(str(path), sd, None).get('fem')
        if not bundle:
            fails.append(f"B: no solved run beside {path.name}")
            continue
        _converged_unchanged(fails, path.name, fem_details, bundle['fem_data'],
                             bundle['solution'], sd, bundle.get('failure_solution'))


def run():
    fails = []
    for check in (check_shipped_runs, check_tips_fixed):
        try:
            check(fails)
        except Exception as exc:                             # pragma: no cover
            fails.append(f"{check.__name__} raised {type(exc).__name__}: {exc}")
    return fails


if __name__ == '__main__':
    failures = run()
    for f in failures:
        print("FAIL", f)
    print("pile_hinge_end:", "FAIL" if failures else "PASS")
    sys.exit(1 if failures else 0)
