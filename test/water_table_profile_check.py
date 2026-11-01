"""The water table read from a seepage field must be the highest point in each
column where the pressure is positive.

A seepage solution carries negative pore pressure above the phreatic surface, and
the slicer and the finite element engine both weigh soil by which side of that
surface it lies on. Each column has to be read between its own top and bottom of
mesh, on the signed field, with zero counted as dry.

The fields here are piecewise linear with kinks on node rows, so tri3 elements
reproduce them exactly and the expected surface is known:

1.  A water table falling across a flat block, below the top everywhere.
2.  Control: a head above the top reads the top of the mesh.
3.  Control: a head below the base reads the base of the mesh.
4.  A rising base under a flat top, where the base climbs past the water table.
5.  A rising base and a rising top.
6.  A field with no suction in it (an imported saturated-only solution).
7.  A wet band over suction, as under a dry foundation.
8.  A perched zone above a second wet zone: the highest wet point governs.
9.  The seepage field replaced after a first read.
10. A slope face on a flat base, wet at the toe.

Run directly:  PYTHONPATH=. python3 test/water_table_profile_check.py
"""

import numpy as np

from xslope.water import seep_water_table_profile

GAMMA_W = 62.4
LENGTH, HEIGHT = 100.0, 20.0
TOL = 1e-4  # the end columns are read a nudge inside the mesh


def _block(nx=20, ny=4, top=None, bot=None):
    top = top or (lambda x: HEIGHT + 0.0 * x)
    bot = bot or (lambda x: 0.0 * x)
    xs = np.linspace(0.0, LENGTH, nx + 1)
    nodes = [(x, b + j / ny * (t - b)) for j in range(ny + 1)
             for x, t, b in zip(xs, top(xs), bot(xs))]
    tris = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            b, c, d = a + 1, a + nx + 1, a + nx + 2
            tris += [[a, b, d], [a, d, c]]
    elements = np.zeros((len(tris), 8), dtype=int)
    elements[:, :3] = tris
    return {'nodes': np.array(nodes), 'elements': elements,
            'element_types': np.full(len(tris), 3, dtype=int)}


def _profile(head, mesh=None):
    mesh = mesh or _block()
    x, y = mesh['nodes'][:, 0], mesh['nodes'][:, 1]
    return seep_water_table_profile({'mesh': mesh, 'seep_u': GAMMA_W * (head(x) - y)})


def _column_field(f):
    """A field of y alone, on a block with a node row every foot."""
    mesh = _block(ny=int(HEIGHT))
    return seep_water_table_profile({'mesh': mesh, 'seep_u': GAMMA_W * f(mesh['nodes'][:, 1])})


def _worst(wt, expected):
    return float(np.max(np.abs(wt - expected)))


def run():
    failures = []
    falling = lambda x: 15.0 - 0.05 * np.asarray(x, dtype=float)

    def check(name, prof, expected, tol=TOL):
        xs, wt = prof
        if not np.isfinite(wt).all():
            failures.append(f"{name}: {int(np.isnan(wt).sum())} columns came back blank")
            return
        worst = _worst(wt, expected(xs))
        if worst > tol:
            failures.append(f"{name}: misread by {worst:.3g}")

    check("falling water table", _profile(falling), falling)
    check("head above the top", _profile(lambda x: np.full_like(x, HEIGHT + 5.0)),
          lambda x: np.full_like(x, HEIGHT))
    check("head below the base", _profile(lambda x: np.full_like(x, -5.0)),
          lambda x: np.zeros_like(x))

    rising = lambda x: 0.1 * x
    check("rising base, flat top",
          _profile(lambda x: np.full_like(x, 8.0), _block(40, 4, bot=rising)),
          lambda x: np.maximum(8.0, rising(x)))
    check("rising base and top",
          _profile(falling, _block(40, 4, top=lambda x: HEIGHT + 0.1 * x, bot=rising)), falling)

    mesh = _block(ny=int(HEIGHT))
    x, y = mesh['nodes'][:, 0], mesh['nodes'][:, 1]
    xs, wt = seep_water_table_profile({'mesh': mesh,
                                       'seep_u': np.maximum(0.0, GAMMA_W * (falling(x) - y))})
    check("field with no suction", (xs, wt), falling, tol=1.0 + TOL)  # one row of nodal data

    band = lambda y: np.minimum(y - 12.0, 16.0 - y)
    check("wet band over suction", _column_field(band), lambda x: np.full_like(x, 16.0))
    perched = lambda y: np.where(y <= 9.0, 6.0 - y, band(y))
    check("perched zone", _column_field(perched), lambda x: np.full_like(x, 16.0))

    mesh = _block()
    x, y = mesh['nodes'][:, 0], mesh['nodes'][:, 1]
    model = {'mesh': mesh, 'seep_u': GAMMA_W * (falling(x) - y)}
    seep_water_table_profile(model)
    model['seep_u'] = GAMMA_W * (falling(x) - 5.0 - y)
    check("replaced field", seep_water_table_profile(model), lambda x: falling(x) - 5.0)

    face = lambda x: np.where(x < 50.0, HEIGHT, HEIGHT - 0.3 * (x - 50.0))
    check("slope face", _profile(falling, _block(40, 4, top=face)),
          lambda x: np.minimum(falling(x), face(x)))

    return failures


if __name__ == '__main__':
    found = run()
    print("\n".join(found) if found else "water table profile: all legs pass")
    raise SystemExit(1 if found else 0)
