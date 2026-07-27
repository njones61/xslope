"""Distributed-load direction must come from geometry, not from point order.

A pressure on a boundary is a physical thing: which way it pushes cannot depend
on which end of the line the points happen to be typed from. It did. The FEM
built the inward normal by rotating the load line's tangent 90 degrees clockwise
— correct only for a line written left-to-right — so RS2-67 b/e/f, whose
downstream pool is traced from the far right boundary back toward the dam toe,
had the reservoir pressure applied as UPLIFT on the submerged bench, leaving a
skin of soil in effective tension that no deformation could relieve.

The main assembly path (Pass 2a, consistent edge-load integration) was fixed by
taking the inside from the element that owns each loaded boundary edge. This
check covers the OTHER path: Pass 2b, the tributary-lumping fallback used when
no complete element edge is found on the load line. No corpus file exercises it
— every real mesh puts complete edges on its load lines — so it needs a
synthetic mesh, which is what this file builds:

  two disjoint triangles hanging below a horizontal load line, arranged so that
  consecutive load-line nodes never share an element. Pass 2a then finds no
  loaded edge and the load falls through to Pass 2b.

Three things are pinned:

  1. Pass 2b is really the path under test (the run must not report the
     "consistent" edge assembly), otherwise the check silently tests Pass 2a;
  2. writing the load line in both directions gives byte-identical assembled
     nodal forces; and
  3. the force points INTO the material both ways — checked with the material
     below the line AND with the mesh mirrored so the material is above it, so
     the rule is geometric rather than a hardcoded "loads push down".

Run directly:  PYTHONPATH=. python3 test/dload_direction_check.py
"""

import contextlib
import io
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from xslope.fem import build_fem_data

Y_LINE = 10.0       # the load line elevation
P = 100.0           # uniform normal pressure on the line
LENGTH = 10.0       # load line length -> total force P * LENGTH


def _mesh(material_above):
    """Two disjoint triangles on one side of the horizontal load line y = Y_LINE.

    The two nodes ON the line belong to different triangles, so no element owns
    an edge lying in the line and the distributed load falls through to Pass 2b.
    ``material_above`` puts the triangles above the line instead of below.
    """
    dy = 10.0 if material_above else -10.0
    nodes = np.array([
        [0.0, Y_LINE],                 # 0  on the load line
        [LENGTH, Y_LINE],              # 1  on the load line
        [0.0, Y_LINE + dy],            # 2
        [0.4 * LENGTH, Y_LINE + dy],   # 3
        [0.6 * LENGTH, Y_LINE + dy],   # 4
        [LENGTH, Y_LINE + dy],         # 5
    ])
    elements = np.zeros((2, 9), dtype=int)
    elements[0, :3] = [0, 2, 3]
    elements[1, :3] = [1, 4, 5]
    return {
        'nodes': nodes,
        'elements': elements,
        'element_types': np.array([3, 3]),
        'element_materials': np.array([1, 1]),
        'elements_1d': np.zeros((0, 3), dtype=int),
        'element_types_1d': np.zeros(0, dtype=int),
        'element_materials_1d': np.zeros(0, dtype=int),
    }


def _slope_data(reversed_order):
    line = [{'X': 0.0, 'Y': Y_LINE, 'Normal': P},
            {'X': LENGTH, 'Y': Y_LINE, 'Normal': P}]
    if reversed_order:
        line = list(reversed(line))
    return {
        'materials': [{'name': 'soil', 'c': 10.0, 'phi': 30.0, 'gamma': 20.0,
                       'E': 1.0e5, 'nu': 0.3, 'u': 'none'}],
        'dloads': [line],
        'dloads2': [],
        'gamma_water': 9.81,
        'k_seismic': 0.0,
        'max_depth': 0.0,
        'reinforce_lines': [],
        'pile_lines': [],
        'line_loads': [],
        'seepage_solution': None,
    }


def _assemble(material_above, reversed_order):
    """Return (nodal force array on the two load nodes, verbose report)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        fem_data = build_fem_data(_slope_data(reversed_order),
                                  mesh=_mesh(material_above), verbose=True)
    return fem_data['bc_values'][[0, 1]].copy(), buf.getvalue()


def main():
    print("Distributed-load direction check (fem.py Pass 2b, tributary lumping):")
    failures = []

    for material_above in (False, True):
        side = "above" if material_above else "below"
        # into the material is +y when the soil is above the line, -y when below
        inward_y = 1.0 if material_above else -1.0
        f_fwd, report = _assemble(material_above, False)
        f_rev, _ = _assemble(material_above, True)

        if 'consistent' in report:
            failures.append(
                f"material {side}: the load was assembled by Pass 2a (consistent "
                "edge integration), so Pass 2b was never exercised — the synthetic "
                "mesh no longer forces the fallback")

        print(f"  material {side} the load line:")
        print(f"    left-to-right  fy = {f_fwd[:, 1]}")
        print(f"    right-to-left  fy = {f_rev[:, 1]}")

        if not np.allclose(f_fwd, f_rev, rtol=0, atol=1e-9):
            failures.append(
                f"material {side}: reversing the point order changed the assembled "
                f"nodal forces ({f_fwd.tolist()} vs {f_rev.tolist()})")

        expected_total = P * LENGTH
        for label, f in (("left-to-right", f_fwd), ("right-to-left", f_rev)):
            fy = float(f[:, 1].sum())
            fx = float(f[:, 0].sum())
            if fy * inward_y <= 0.0:
                failures.append(
                    f"material {side}, {label}: total vertical force {fy:+.1f} points "
                    f"AWAY from the material (expected sign {inward_y:+.0f}) — the "
                    "pressure is being applied as uplift")
            if abs(abs(fy) - expected_total) > 1e-6 or abs(fx) > 1e-9:
                failures.append(
                    f"material {side}, {label}: assembled force ({fx:+.3f}, {fy:+.3f}) "
                    f"does not match the applied pressure resultant "
                    f"(0, {inward_y * expected_total:+.3f})")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nPass 2b applies the load into the material regardless of point order.")


if __name__ == "__main__":
    main()
