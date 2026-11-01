"""Reading a nodal field at a point must honor where the point sits in its element.

Seepage pore pressures reach the slice bases and the water table through
``interpolate_at_point``. An isoparametric element reproduces a linear field
exactly, whatever the shape of the element, so a linear field on one straight-sided
element of each type must read exactly at every point inside it. Quadratic elements
on a rectangle also reproduce a quadratic field exactly.

1. A linear field on a rectangle, for every element type.
2. A linear field on a distorted quadrilateral, or a general triangle.
3. A quadratic field on a rectangle, for the quadratic elements.
4. A hydrostatic field on a two-element column of quad8s, read across both.

Run directly:  PYTHONPATH=. python3 test/mesh_interpolation_check.py
"""

import numpy as np

from xslope.mesh import interpolate_at_point

TOL = 1e-8
RECTANGLE = [(0.0, 0.0), (8.0, 0.0), (8.0, 5.0), (0.0, 5.0)]
DISTORTED = [(0.0, 0.0), (8.0, 1.0), (9.0, 6.0), (1.0, 4.0)]
TRIANGLE = [(0.0, 0.0), (8.0, 1.0), (3.0, 6.0)]


def _element(corners, n):
    """Nodes of one straight-sided element of type n, in XSLOPE's node order."""
    c = np.array(corners, dtype=float)
    if n in (3, 4):
        return c
    mids = 0.5 * (c + np.roll(c, -1, axis=0))
    return np.vstack([c, mids] + ([c.mean(axis=0, keepdims=True)] if n == 9 else []))


def _points(corners, count=200):
    """Random points strictly inside a convex polygon."""
    c = np.array(corners, dtype=float)
    rng = np.random.default_rng(7)
    w = rng.dirichlet(np.ones(len(c)), count)
    return 0.05 * c.mean(axis=0) + 0.95 * w @ c


def _worst(nodes, elements, types, field):
    values = field(nodes[:, 0], nodes[:, 1])
    corners = nodes[elements[0][:4 if types[0] in (4, 8, 9) else 3]]
    worst = 0.0
    for x, y in _points(corners):
        got, found = interpolate_at_point(nodes, elements, types, values, (x, y),
                                          return_found=True, signed=True)
        worst = max(worst, np.inf if not found else abs(got - field(x, y)))
    return worst


def _one(corners, n):
    nodes = _element(corners, n)
    elements = np.zeros((1, 9), dtype=int)
    elements[0, :n] = np.arange(n)
    return nodes, elements, np.array([n])


def run():
    failures = []
    linear = lambda x, y: 3.0 + 2.0 * x - 5.0 * y
    quadratic = lambda x, y: 1.0 + x + 2.0 * y + 0.5 * x * x - x * y + 0.25 * y * y

    def check(name, worst):
        if worst > TOL:
            failures.append(f"{name}: misread by {worst:.3g}")

    for n in (3, 6):
        check(f"type {n}, linear field", _worst(*_one(TRIANGLE, n), linear))
    for n in (4, 8, 9):
        check(f"type {n}, linear field on a rectangle", _worst(*_one(RECTANGLE, n), linear))
        check(f"type {n}, linear field on a distorted quad", _worst(*_one(DISTORTED, n), linear))
    check("type 6, quadratic field", _worst(*_one(TRIANGLE, 6), quadratic))
    for n in (8, 9):
        check(f"type {n}, quadratic field", _worst(*_one(RECTANGLE, n), quadratic))

    # two quad8s stacked: 0-1-2-3 below, 3-2-4-5 above, midsides after the corners
    nodes = np.array([(0, 0), (5, 0), (5, 5), (0, 5), (5, 10), (0, 10), (2.5, 0), (5, 2.5),
                      (2.5, 5), (0, 2.5), (5, 7.5), (2.5, 10), (0, 7.5)], dtype=float)
    elements = np.array([[0, 1, 2, 3, 6, 7, 8, 9, 0], [3, 2, 4, 5, 8, 10, 11, 12, 0]])
    types = np.array([8, 8])
    head = lambda x, y: 62.4 * (7.0 - y)
    worst = 0.0
    values = head(nodes[:, 0], nodes[:, 1])
    for y in (0.2, 2.4, 4.8, 5.2, 7.6, 9.8):
        got = interpolate_at_point(nodes, elements, types, values, (1.7, y), signed=True)
        worst = max(worst, abs(got - head(1.7, y)))
    check("hydrostatic column of quad8s", worst)
    return failures


if __name__ == '__main__':
    found = run()
    print("\n".join(found) if found else "mesh interpolation: all legs pass")
    raise SystemExit(1 if found else 0)
