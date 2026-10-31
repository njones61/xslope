"""The water table read from a seepage field must be the field's u = 0 contour.

A seepage solution carries negative pore pressure above the phreatic surface, and
the slicer and the finite element engine both weigh soil by which side of that
surface it lies on. ``seep_water_table_profile`` therefore has to read the signed
field. Read with suction clamped to zero, the pressure at the top of a column looks
like "not negative", and every column that reaches the top elevation of the mesh
(a flat crest, say) is taken as saturated to the ground surface.

The field here is linear, u = gamma_w (h(x) - y) on a flat-topped block, so the
contour is known exactly and tri3 elements reproduce it exactly:

1. A water table falling across the block, below the top everywhere.
2. Control: a head above the top reads the top of the mesh.
3. Control: a head below the base reads the base of the mesh.

Run directly:  PYTHONPATH=. python3 test/water_table_profile_check.py
"""

import numpy as np

from xslope.water import seep_water_table_profile

GAMMA_W = 62.4
LENGTH, HEIGHT = 100.0, 20.0


def _block(nx=20, ny=4):
    xs, ys = np.meshgrid(np.linspace(0.0, LENGTH, nx + 1), np.linspace(0.0, HEIGHT, ny + 1))
    nodes = np.column_stack([xs.ravel(), ys.ravel()])
    tris = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            b, c, d = a + 1, a + nx + 1, a + nx + 2
            tris += [[a, b, d], [a, d, c]]
    elements = np.zeros((len(tris), 8), dtype=int)
    elements[:, :3] = tris
    return {'nodes': nodes, 'elements': elements,
            'element_types': np.full(len(tris), 3, dtype=int)}


def _profile(head):
    mesh = _block()
    x, y = mesh['nodes'][:, 0], mesh['nodes'][:, 1]
    return seep_water_table_profile({'mesh': mesh, 'seep_u': GAMMA_W * (head(x) - y)})


def run():
    failures = []

    falling = lambda x: 15.0 - 0.05 * np.asarray(x, dtype=float)
    xs, wt = _profile(falling)
    worst = float(np.max(np.abs(wt - falling(xs))))
    if worst > 1e-6:
        failures.append(f"falling water table misread by {worst:.3g}: every column reaches "
                        f"the mesh top, so a clamped read returns {HEIGHT:g}")

    xs, wt = _profile(lambda x: np.full_like(x, HEIGHT + 5.0))
    if not np.allclose(wt, HEIGHT):
        failures.append(f"a head above the top should read the top, got {wt.min():g} to {wt.max():g}")

    xs, wt = _profile(lambda x: np.full_like(x, -5.0))
    if not np.allclose(wt, 0.0):
        failures.append(f"a head below the base should read the base, got {wt.min():g} to {wt.max():g}")

    return failures


if __name__ == '__main__':
    found = run()
    print("\n".join(found) if found else "water table profile: all legs pass")
    raise SystemExit(1 if found else 0)
