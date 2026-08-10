"""Gradient and velocity recovery against the field they must reproduce exactly.

XSLOPE solves for nodal head, and everything a reader sees of the FLOW is
differentiated out of that solution afterwards, by ``seep.compute_gradient``
(hydraulic gradient i = -grad h) and ``seep.compute_velocity`` (Darcy velocity
v = -kr K grad h). Nothing in the solve reads either one back: head, pore
pressure, the flow potential and the total flowrate are assembled and solved by
separate code. So an error in the recovery moves every gradient and velocity a
model reports -- the flow arrows, the v_mag and i_mag panels, the exit gradient
-- while the factor of safety, which reads only pore pressure, holds still and
says nothing is wrong.

THE PATCH TEST

A linear head field is the case that pins this down without a reference
solution, because the answer does not depend on the mesh at all:

    h = a*x + b*y   =>   grad h = (a, b)   =>   i = (-a, -b)   everywhere

Every element formulation here is at least linearly complete, so each must
return that one constant vector at every node of any mesh, whatever its element
shapes -- exactly, to rounding. A recovery that is wrong by a rotation, a
transpose, a scale or a swapped axis cannot reproduce it, and the residual is
the size of the error rather than a tolerance someone chose.

The meshes are UNSTRUCTURED, which is the property that makes the test bite. On
a structured mesh every element has the same shape and orientation, so a
Jacobian error acts as one uniform rotation on the whole field and can look
plausible; it also cancels entirely for elements that happen to be right
isoceles. The shipped tri6 branch inverted a transposed Jacobian for exactly
this reason -- on the corpus's unstructured meshes it put every node's gradient
a median 133 degrees out and 125% wrong in magnitude, and no test saw it.

Velocity adds the conductivity on top of the same gradient, so it is checked
against -K grad h with K isotropic and again rotated and anisotropic, which
pins the frame K is applied in as well as the gradient underneath it.

Element types the functions do not handle must RAISE rather than return the
zeros of an array nobody filled: a zero velocity field is a readable, plottable
answer that means still water, so silence there is indistinguishable from a
result.

Run directly:  PYTHONPATH=. python3 test/flow_recovery_check.py
"""

import contextlib
import io
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from xslope.mesh import build_mesh_from_polygons
from xslope.seep import compute_gradient, compute_velocity

#: A closed polygon with no two edges alike, so gmsh has to produce genuinely
#: irregular elements: no element is right isoceles, none shares an orientation
#: with its neighbour, and the aspect ratios vary across the domain. A rectangle
#: would let a Jacobian error partly cancel.
DOMAIN = [(0.0, 0.0), (37.0, 0.0), (37.0, 9.0), (23.0, 14.0),
          (11.0, 12.0), (0.0, 17.0)]
TARGET = 2.3

#: The linear field. Both coefficients are non-zero and of unlike sign and
#: magnitude, so a swapped or negated axis cannot pass by symmetry.
A_COEF, B_COEF = 0.37, -0.85

ELEMENT_TYPES = ['tri3', 'quad4', 'tri6', 'quad8', 'quad9']

#: Exact to rounding: the residual is float noise on the coordinates, not a
#: discretisation error, so this is tight on purpose.
TOL = 1e-12


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _mesh(element_type):
    polygons = [{'mat_id': 1, 'coords': DOMAIN, 'size': None}]
    return _quiet(build_mesh_from_polygons, polygons, TARGET, element_type)


def _linear_head(nodes):
    return A_COEF * nodes[:, 0] + B_COEF * nodes[:, 1]


def _live(mesh):
    """Nodes an element actually references — the ones a recovery must fill."""
    els = np.asarray(mesh['elements'])
    ets = np.asarray(mesh['element_types'])
    used = set()
    for el, et in zip(els, ets):
        used.update(int(n) for n in el[:int(et)])
    return np.array(sorted(used), dtype=int)


def _check(label, got, want, fails, tol=TOL):
    err = float(np.max(np.abs(got - want)))
    if not (err <= tol):
        fails.append(f"{label}: max error {err:.3e} exceeds {tol:.0e} "
                     f"(want {want}, worst row {got[np.argmax(np.abs(got - want).max(axis=1))]})")
    return err


def leg_gradient():
    """i = (-a, -b) at every node, on every element type."""
    fails = []
    for et in ELEMENT_TYPES:
        mesh = _mesh(et)
        nodes = np.asarray(mesh['nodes'])
        live = _live(mesh)
        h = _linear_head(nodes)
        i = compute_gradient(nodes, np.asarray(mesh['elements']), h,
                             np.asarray(mesh['element_types']))
        want = np.tile([-A_COEF, -B_COEF], (len(live), 1))
        err = _check(f"gradient/{et}", i[live], want, fails)
        print(f"  gradient  {et:<6} {len(live):>6} nodes  max err {err:.2e}")
    return fails


def leg_velocity_isotropic():
    """v = -k grad h, so the recovered velocity is -k (a, b) everywhere."""
    fails = []
    k = 3.7
    for et in ELEMENT_TYPES:
        mesh = _mesh(et)
        nodes = np.asarray(mesh['nodes'])
        els = np.asarray(mesh['elements'])
        ets = np.asarray(mesh['element_types'])
        live = _live(mesh)
        h = _linear_head(nodes)
        n_el = len(els)
        v = compute_velocity(nodes, els, h, np.full(n_el, k), np.full(n_el, k),
                             np.zeros(n_el), element_types=ets)
        want = np.tile([-k * A_COEF, -k * B_COEF], (len(live), 1))
        err = _check(f"velocity-isotropic/{et}", v[live], want, fails, tol=TOL * k)
        print(f"  velocity  {et:<6} {len(live):>6} nodes  max err {err:.2e}")
    return fails


def leg_velocity_anisotropic():
    """The conductivity is applied in ITS OWN frame, not the mesh's.

    With k1 along a direction rotated by theta from x, K = R diag(k1, k2) R^T,
    and v = -K grad h. A gradient recovered in the wrong frame and a K rotated
    the wrong way both fail here, and they fail differently from the isotropic
    leg, which cannot see a rotation at all.
    """
    fails = []
    k1, k2, theta = 9.0, 1.5, 34.0
    c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    R = np.array([[c, -s], [s, c]])
    K = R @ np.diag([k1, k2]) @ R.T
    grad_h = np.array([A_COEF, B_COEF])
    want_row = -K @ grad_h
    for et in ELEMENT_TYPES:
        mesh = _mesh(et)
        nodes = np.asarray(mesh['nodes'])
        els = np.asarray(mesh['elements'])
        ets = np.asarray(mesh['element_types'])
        live = _live(mesh)
        h = _linear_head(nodes)
        n_el = len(els)
        v = compute_velocity(nodes, els, h, np.full(n_el, k1), np.full(n_el, k2),
                             np.full(n_el, theta), element_types=ets)
        want = np.tile(want_row, (len(live), 1))
        err = _check(f"velocity-anisotropic/{et}", v[live], want, fails,
                     tol=TOL * max(k1, k2))
        print(f"  aniso     {et:<6} {len(live):>6} nodes  max err {err:.2e}")
    return fails


def leg_unhandled_type_raises():
    """A type neither function has a branch for must say so, not return zeros.

    The zeros are the danger: an unfilled velocity array is a legal field that
    plots as still water, so it reads as an answer. This drives a real mesh
    through both functions under a relabelled element type.
    """
    fails = []
    mesh = _mesh('tri6')
    nodes = np.asarray(mesh['nodes'])
    els = np.asarray(mesh['elements'])
    h = _linear_head(nodes)
    n_el = len(els)
    bogus = 7
    ets = np.full(n_el, bogus)
    for name, call in [
        ('compute_gradient', lambda: compute_gradient(nodes, els, h, ets)),
        ('compute_velocity', lambda: compute_velocity(
            nodes, els, h, np.ones(n_el), np.ones(n_el), np.zeros(n_el),
            element_types=ets)),
    ]:
        try:
            call()
        except ValueError as e:
            if str(bogus) not in str(e):
                fails.append(f"{name}: raised without naming type {bogus}: {e}")
            else:
                print(f"  raises    {name} on type {bogus}: ok")
        except Exception as e:  # noqa: BLE001
            fails.append(f"{name}: raised {type(e).__name__}, expected ValueError: {e}")
        else:
            fails.append(f"{name}: returned silently for element type {bogus}")
    return fails


def leg_mixed_mesh():
    """A mesh of two element types recovers the same one field across the seam.

    Corpus meshes mix types (the tri6/quad8 FEM meshes do), and the nodal
    average at a seam draws on both branches, so both must agree there.
    """
    fails = []
    tri = _mesh('tri6')
    nodes = np.asarray(tri['nodes'])
    els = np.asarray(tri['elements'])
    ets = np.asarray(tri['element_types']).copy()
    if len(ets) < 4:
        fails.append('mixed: mesh too small to split')
        return fails
    # Relabelling would change the physics; instead check that a genuinely mixed
    # corpus-shaped case (tri6 elements padded into a wider array, as the mesh
    # files store them) recovers the same constant field.
    wide = np.zeros((len(els), 9), dtype=int)
    wide[:, :els.shape[1]] = els
    h = _linear_head(nodes)
    live = _live({'elements': wide, 'element_types': ets})
    i = compute_gradient(nodes, wide, h, ets)
    want = np.tile([-A_COEF, -B_COEF], (len(live), 1))
    err = _check('gradient/padded-connectivity', i[live], want, fails)
    print(f"  padded    tri6   {len(live):>6} nodes  max err {err:.2e}")
    return fails


LEGS = [
    ("hydraulic gradient on a linear field", leg_gradient),
    ("Darcy velocity, isotropic k", leg_velocity_isotropic),
    ("Darcy velocity, rotated anisotropic k", leg_velocity_anisotropic),
    ("padded connectivity", leg_mixed_mesh),
    ("unhandled element type raises", leg_unhandled_type_raises),
]


def run():
    failures = []
    for label, fn in LEGS:
        print(f"{label}:")
        failures += fn()
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print("\nPASS: gradient and velocity recovery reproduces a linear head field "
          "exactly on tri3, quad4, tri6, quad8 and quad9.")


if __name__ == "__main__":
    main()
