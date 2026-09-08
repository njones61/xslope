"""A reinforcement line's mesh geometry is its two endpoints, whatever the law.

A reinforcement line carries two different descriptions of itself. One is the
line: two endpoints, a straight member the mesh has to conform to. The other is
the LEM tension distribution — the list of points at which the available-tension
envelope is stored, so the limit-equilibrium engine can read the capacity at the
crossing by interpolating between them. That list depends on the pullout law:
under the constant-rate law it holds the envelope's breakpoints (the ramp knees),
and under the overburden-dependent law the envelope is a curve, so the list is a
dense sampling of it — forty-one points on a line whose geometry is two.

``extract_reinforcement_line_geometry`` used to hand gmsh the second list. Every
stored tension point became a mesh vertex, so filling in the Adhesion and Delta
columns of the reinforce sheet — a statement about bond strength, not about
discretization — silently re-meshed the model: on the FEM-2 tutorial slope, 2,101
triangles and 60 bar elements became 6,569 and 240, and the factor of safety fell
8%. A capacity law cannot be allowed to choose a discretization.

The rule this check pins: whatever the law, a reinforcement line reaches the
mesher as exactly two points, and they are the line's own endpoints. The 1D
element size decides the subdivision, as it does for a pile.

Three things are pinned, on one real reinforced model:

  1. under the constant-rate law the mesh geometry is 2 points per line, at the
     line's stated endpoints;
  2. with Adhesion and Delta filled in — the overburden law — the mesh geometry
     is still 2 points per line, at the same endpoints, and byte-identical to
     the constant-rate geometry; and
  3. the law really is live on the fixture, so (2) is not passing because
     nothing changed: the LEM tension-distribution list itself must GROW when
     the law is entered. Without this leg a loader that ignored Adhesion/Delta
     entirely would pass the check.

Run directly:  PYTHONPATH=. python3 test/reinforce_mesh_geometry_check.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from xslope.fileio import (load_slope_data, ensure_reinforce_pullout,
                           build_reinforce_lines)
from xslope.mesh import extract_reinforcement_line_geometry

#: A shipped model with reinforcement lines on the constant-rate law. The docs/fem
#: sample rather than a synthetic one, so the check is scored on geometry a user
#: really has: eight geotextile layers, all of them straight two-endpoint members
#: whose stored point lists carry ramp knees.
MODEL = os.path.join(_ROOT, 'docs/fem/files/xslope_reinforce_fem.xlsx')

#: Adhesion and interface friction angle written onto every line to put the
#: fixture on the overburden law. Any finite pair does it; these are the values
#: the FEM-2 tutorial uses (adhesion 0, delta 22 degrees), which is the case that
#: exposed the defect.
ADHESION = 0.0
DELTA = 22.0

TOL = 1e-9


def _endpoints(slope_data):
    """Each raw reinforcement row's stated endpoints, in file order."""
    return [((r['x1'], r['y1']), (r['x2'], r['y2']))
            for r in slope_data['reinforcement_lines']]


def _check_geometry(slope_data, law, failures):
    """Legs 1 and 2: two points per line, at the line's own endpoints."""
    geom = extract_reinforcement_line_geometry(slope_data)
    ends = _endpoints(slope_data)
    if len(geom) != len(ends):
        failures.append(f"{law}: {len(geom)} mesh lines for "
                        f"{len(ends)} reinforcement rows")
        return geom
    for i, (pts, (p1, p2)) in enumerate(zip(geom, ends), start=1):
        stored = len(slope_data['reinforce_lines'][i - 1])
        if len(pts) != 2:
            failures.append(
                f"{law}: line {i} reaches the mesher as {len(pts)} points, not 2 "
                f"(its stored tension distribution has {stored})")
            continue
        for got, want, end in ((pts[0], p1, 'start'), (pts[1], p2, 'end')):
            if (abs(got[0] - want[0]) > TOL) or (abs(got[1] - want[1]) > TOL):
                failures.append(
                    f"{law}: line {i} {end} point is {got}, not the stated {want}")
    return geom


def main():
    failures = []

    if not os.path.exists(MODEL):
        return [f"fixture model missing: {MODEL}"]

    # --- constant-rate law, as the file states it -------------------------
    rate = load_slope_data(MODEL)
    n_lines = len(rate.get('reinforcement_lines') or [])
    if n_lines < 1:
        return [f"fixture model has no reinforcement lines: {MODEL}"]
    rate_stored = [len(pts) for pts in rate['reinforce_lines']]
    rate_geom = _check_geometry(rate, 'constant-rate law', failures)

    # --- the same model with Adhesion/Delta filled ------------------------
    law = load_slope_data(MODEL)
    for r in law['reinforcement_lines']:
        r['adhesion'] = ADHESION
        r['delta'] = DELTA
    ensure_reinforce_pullout(law)
    law['reinforce_lines'] = build_reinforce_lines(law['reinforcement_lines'])
    law_stored = [len(pts) for pts in law['reinforce_lines']]
    law_geom = _check_geometry(law, 'overburden law', failures)

    # --- leg 3: the law is live on this fixture ---------------------------
    if not any(b > a for a, b in zip(rate_stored, law_stored)):
        failures.append(
            "the overburden law did not change any line's stored tension "
            f"distribution ({rate_stored} -> {law_stored}), so the check would "
            "pass on a loader that ignored Adhesion/Delta — the fixture no "
            "longer exercises the rule")

    # --- the two geometries are the same line -----------------------------
    if len(rate_geom) == len(law_geom):
        for i, (a, b) in enumerate(zip(rate_geom, law_geom), start=1):
            if a != b:
                failures.append(f"line {i} meshes differently under the two "
                                f"laws: {a} vs {b}")
    else:
        failures.append(f"{len(rate_geom)} mesh lines under the constant-rate "
                        f"law, {len(law_geom)} under the overburden law")

    if not failures:
        print(f"{n_lines} reinforcement lines, stored tension points "
              f"{rate_stored} -> {law_stored}; mesh geometry 2 points per line "
              f"under both laws.")
    return failures


def run():
    """Failures as a list, for run_tests.py."""
    return main()


def _cli():
    failures = main()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nA reinforcement line meshes by its endpoints, whatever the "
          "pullout law.")


if __name__ == "__main__":
    _cli()
