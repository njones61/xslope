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

"""Reinforcement and pile elements against the soil edge they lie on.

An embedded one-dimensional element shares its nodes with the two-dimensional
soil elements around it, and that shared node is the whole of the coupling: the
member and the soil displace together wherever they meet, and nowhere else. On a
quadratic soil mesh each edge carries a midside node as well as its two corners,
so an element that stands only on the corners is tied to the soil at half the
stations the edge has. The soil edge can then bow away from the member between
the corners, and the member is free to slide past the node in the middle of it.

What is checked here:

  A. THE MESH ATTACHES THE MIDSIDE NODE -- on a quadratic mesh every embedded 1D
     element records three nodes, the third is the midpoint of the first two, and
     it is a midside node of an adjacent 2D element. Checked on a freshly meshed
     model and on a mesh written to JSON and read back, since a mesh file stored
     before the elements carried the node must come back repaired rather than
     silently keeping the old element.

  B. THE PATCH TEST -- a bar on a soil edge under a uniform axial strain carries
     exactly EA*strain, the midside node takes no share of the nodal forces (a
     constant axial force is [-N, +N, 0] on the three-node bar), and the midside
     node's degrees of freedom really are in the assembled stiffness. A bar
     assembled against the corners alone passes the first of those and fails the
     last, which is the defect this element replaced.

  C. THE LINEAR MESH DOES NOT MOVE -- on a tri3 mesh there is no midside node to
     attach, the elements stay two-node, and the assembled stiffness and the
     recovered member forces are bit-for-bit what they are with the
     XSLOPE_LINEAR_1D escape hatch set. Nothing about a linear mesh changed.

Run directly:  PYTHONPATH=. python3 test/one_d_compatibility_check.py
"""

import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np                                              # noqa: E402


# --------------------------------------------------------------------------- #
# A small meshed model: one rectangular soil block with a line through it
# --------------------------------------------------------------------------- #

#: The block, and the line the member runs along. The line's ends sit on the
#: block's own boundary vertices so the mesher has nothing to recover.
_BLOCK = [(0.0, 0.0), (12.0, 0.0), (12.0, 8.0), (0.0, 8.0)]
_LINE = [(0.0, 4.0), (12.0, 4.0)]


def _polygon_with_line():
    """The block with the line's ends inserted into its boundary, which is what
    ``get_material_polygons`` does for a real model."""
    return [[(0.0, 0.0), (12.0, 0.0), (12.0, 4.0), (12.0, 8.0),
             (0.0, 8.0), (0.0, 4.0)]]


def build_block(element_type='tri6', size=2.0):
    """One meshed block with a constraint line across it."""
    from xslope.mesh import build_mesh_from_polygons
    return build_mesh_from_polygons(_polygon_with_line(), target_size=size,
                                    element_type=element_type, lines=[_LINE])


def reinforced_slope_data(E=2.0e7, area=0.002, t_max=1.0e6):
    """A slope_data carrying one reinforcement line along ``_LINE``.

    The capacity is set far above anything the patch test mobilizes, so the bar
    stays elastic and the check reads the element rather than the cap.
    """
    return {
        'materials': [{'name': 'soil', 'gamma': 20.0, 'E': 30000.0, 'nu': 0.3,
                       'option': 'mc', 'c': 25.0, 'phi': 30.0, 'u': 'none'}],
        'gamma_water': 9.81,
        'reinforcement_lines': [{
            'x1': _LINE[0][0], 'y1': _LINE[0][1],
            'x2': _LINE[1][0], 'y2': _LINE[1][1],
            'E': E, 'area': area, 't_max': t_max, 'lp1': 0.0, 'lp2': 0.0,
            'tend1': t_max, 'tend2': t_max,
        }],
    }


# --------------------------------------------------------------------------- #
# A. the mesh attaches the midside node
# --------------------------------------------------------------------------- #

def _edge_midside_map(mesh):
    """(corner, corner) -> midside node, over every quadratic 2D element."""
    out = {}
    spec = {6: (((0, 1), 3), ((1, 2), 4), ((2, 0), 5)),
            8: (((0, 1), 4), ((1, 2), 5), ((2, 3), 6), ((3, 0), 7)),
            9: (((0, 1), 4), ((1, 2), 5), ((2, 3), 6), ((3, 0), 7))}
    for element, elem_type in zip(np.asarray(mesh['elements']),
                                  np.asarray(mesh['element_types'])):
        for (a, b), m in spec.get(int(elem_type), ()):
            n_a, n_b = int(element[a]), int(element[b])
            out[(min(n_a, n_b), max(n_a, n_b))] = int(element[m])
    return out


def _check_attached(mesh, label, fails):
    """Every 1D element of a quadratic mesh stands on three nodes, and the third
    is the midside node of the 2D edge its two corners span."""
    nodes = np.asarray(mesh['nodes'], dtype=float)
    e1d = np.asarray(mesh['elements_1d'], dtype=int)
    types = np.asarray(mesh['element_types_1d'], dtype=int)
    midside = _edge_midside_map(mesh)

    n_two = int((types < 3).sum())
    if n_two:
        fails.append(f"{label}: {n_two} of {len(e1d)} 1D elements are still "
                     f"two-node on a quadratic mesh")
    for i in range(len(e1d)):
        if int(types[i]) < 3:
            continue
        n_0, n_1, n_m = int(e1d[i, 0]), int(e1d[i, 1]), int(e1d[i, 2])
        want = midside.get((min(n_0, n_1), max(n_0, n_1)))
        if want is None:
            fails.append(f"{label}: 1D element {i} does not lie on a 2D edge")
            continue
        if want != n_m:
            fails.append(f"{label}: 1D element {i} records node {n_m}, but the "
                         f"2D edge's midside node is {want}")
            continue
        offset = np.linalg.norm(nodes[n_m] - 0.5 * (nodes[n_0] + nodes[n_1]))
        if offset > 1e-9:
            fails.append(f"{label}: 1D element {i}'s third node is {offset:.2e} "
                         f"from the midpoint of its ends")


def test_mesh_attaches_the_midside_node():
    """A quadratic mesh delivers three-node 1D elements, from the mesher and
    from a mesh file read back."""
    fails = []
    from xslope.mesh import export_mesh_to_json, import_mesh_from_json

    mesh = build_block('tri6')
    if 'elements_1d' not in mesh or len(mesh['elements_1d']) == 0:
        return ["the meshed block carries no 1D element, so nothing is proven"]
    _check_attached(mesh, "tri6 mesh", fails)

    # A stored mesh: written with the elements stripped back to two nodes, the
    # way a file saved before the midside node was attached carries them.
    stale = {k: (v.copy() if isinstance(v, np.ndarray) else v)
             for k, v in mesh.items()}
    stale['elements_1d'][:, 2] = 0
    stale['element_types_1d'][:] = 2
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'stale_mesh.json')
        import io, contextlib
        with contextlib.redirect_stdout(io.StringIO()):
            export_mesh_to_json(stale, path)
            reloaded = import_mesh_from_json(path)
    _check_attached(reloaded, "stored mesh read back", fails)

    # ... and a linear mesh has no midside node to attach, so its elements stay
    # two-node rather than being given a node that is not there.
    linear = build_block('tri3')
    types = np.asarray(linear.get('element_types_1d', []), dtype=int)
    if len(types) and int(types.max()) != 2:
        fails.append(f"tri3 mesh: 1D elements are recorded with "
                     f"{int(types.max())} nodes, but a linear mesh has no "
                     f"midside node")
    return fails


def test_the_escape_hatch_holds_the_old_element():
    """XSLOPE_LINEAR_1D=1 keeps the elements two-node on a quadratic mesh.

    The hatch exists so a model can be solved both ways while results measured
    against the two-node element are re-measured. It is temporary, and this check
    is what says it still does what it says.
    """
    fails = []
    old = os.environ.get('XSLOPE_LINEAR_1D')
    os.environ['XSLOPE_LINEAR_1D'] = '1'
    try:
        mesh = build_block('tri6')
    finally:
        if old is None:
            os.environ.pop('XSLOPE_LINEAR_1D', None)
        else:
            os.environ['XSLOPE_LINEAR_1D'] = old
    types = np.asarray(mesh.get('element_types_1d', []), dtype=int)
    if not len(types):
        return ["the meshed block carries no 1D element under the hatch"]
    if int(types.max()) != 2:
        fails.append(f"XSLOPE_LINEAR_1D=1 still produced "
                     f"{int(types.max())}-node 1D elements")
    return fails


# --------------------------------------------------------------------------- #
# B. the patch test
# --------------------------------------------------------------------------- #

def test_uniform_axial_strain():
    """A bar on a soil edge under a uniform axial strain.

    Every node of the mesh is displaced by u = (eps*x, 0), a state of uniform
    axial strain along the horizontal member line. Each bar element must then
    carry exactly EA*eps at its center, and the internal nodal forces of the
    three-node bar under a constant axial force are [-N, +N, 0] in its own node
    order -- the midside node takes no share of them, which is why the capped
    bar's body-load correction is applied at the ends alone.
    """
    fails = []
    from xslope.fem import build_fem_data

    mesh = build_block('tri6')
    slope_data = reinforced_slope_data()
    fem_data = build_fem_data(slope_data, mesh)

    n_1d = len(fem_data['elements_1d'])
    if not n_1d:
        return ["the meshed block carries no 1D element, so nothing is proven"]

    eps = 1.0e-4
    nodes = np.asarray(fem_data['nodes'], dtype=float)
    dof_offset = np.asarray(fem_data['dof_offset'], dtype=int)
    u = np.zeros(int(dof_offset[-1]))
    for nid in range(len(nodes)):
        u[dof_offset[nid]] = eps * nodes[nid, 0]

    n_dof_1d = np.asarray(fem_data['n_dof_by_1d_elem'], dtype=int)
    dof_1d = np.asarray(fem_data['dof_indices_1d'], dtype=int)
    k_1d = np.asarray(fem_data['k_by_1d_elem'], dtype=float)
    lengths = np.asarray(fem_data['elem_length_1d'], dtype=float)
    cos_t = np.asarray(fem_data['cos_theta_1d'], dtype=float)
    sin_t = np.asarray(fem_data['sin_theta_1d'], dtype=float)

    n_three = 0
    worst_force, worst_mid = 0.0, 0.0
    for i in range(n_1d):
        n_dof = int(n_dof_1d[i])
        idx = dof_1d[i][:n_dof]
        # Axial force at the element center, the same chord expression the solver
        # uses. EA = k*L, so the reference is EA*eps.
        u_e = u[idx]
        delta = (u_e[2] - u_e[0]) * cos_t[i] + (u_e[3] - u_e[1]) * sin_t[i]
        force = k_1d[i] * delta
        want = k_1d[i] * lengths[i] * eps
        worst_force = max(worst_force, abs(force - want) / max(abs(want), 1e-30))

        # Internal nodal forces of the element in this state.
        f_int = np.asarray(fem_data['K_global_1d_elems'][i], dtype=float) @ u_e
        if n_dof == 6:
            n_three += 1
            mid = np.hypot(f_int[4], f_int[5])
            worst_mid = max(worst_mid, mid / max(abs(want), 1e-30))

    if n_three == 0:
        fails.append("no three-node bar in the patch, so nothing is proven")
    if worst_force > 1e-9:
        fails.append(f"a bar under uniform axial strain carries EA*eps to only "
                     f"{worst_force:.2e} relative")
    if worst_mid > 1e-9:
        fails.append(f"the midside node of a bar under a constant axial force "
                     f"carries {worst_mid:.2e} of it, and should carry none")

    # ... and those midside DOFs are really in the assembled stiffness. A bar
    # tied to the corners alone leaves them untouched by the member, which is
    # exactly what a diagonal contributed by the soil alone cannot show -- so the
    # bar's OWN blocks are assembled here and the midside rows read off them.
    n_dof_total = int(dof_offset[-1])
    reach = np.zeros(n_dof_total)
    for i in range(n_1d):
        n_dof = int(n_dof_1d[i])
        if n_dof < 6:
            continue
        idx = dof_1d[i][:n_dof]
        reach[idx] += np.abs(np.asarray(fem_data['K_global_1d_elems'][i],
                                        dtype=float)).sum(axis=1)
    mid_nodes = sorted({int(np.asarray(fem_data['elements_1d'])[i, 2])
                        for i in range(n_1d) if int(n_dof_1d[i]) == 6})
    silent = [nid for nid in mid_nodes if reach[dof_offset[nid]] <= 0.0]
    if silent:
        fails.append(f"{len(silent)} midside node(s) of a three-node bar carry "
                     f"no member stiffness at all")
    return fails


# --------------------------------------------------------------------------- #
# C. the linear mesh does not move
# --------------------------------------------------------------------------- #

def _tri3_state():
    """(K, forces_1d) for the linear-mesh block, through the ordinary path."""
    import contextlib
    import io
    from xslope.fem import build_fem_data, build_global_stiffness, solve_fem

    mesh = build_block('tri3')
    slope_data = reinforced_slope_data()
    fem_data = build_fem_data(slope_data, mesh)
    K = build_global_stiffness(fem_data['nodes'], fem_data['elements'],
                               fem_data['element_types'],
                               fem_data['element_materials'],
                               fem_data['E_by_mat'], fem_data['nu_by_mat'],
                               fem_data=fem_data)
    with contextlib.redirect_stdout(io.StringIO()):
        solution = solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=25)
    types = np.asarray(fem_data['element_types_1d'], dtype=int)
    return K, np.asarray(solution['forces_1d'], dtype=float), types


def test_a_linear_mesh_is_untouched():
    """On tri3 the elements stay two-node, and the stiffness and member forces
    are bit-for-bit what the two-node escape hatch gives."""
    fails = []
    K_now, forces_now, types_now = _tri3_state()

    old = os.environ.get('XSLOPE_LINEAR_1D')
    os.environ['XSLOPE_LINEAR_1D'] = '1'
    try:
        K_hatch, forces_hatch, _types = _tri3_state()
    finally:
        if old is None:
            os.environ.pop('XSLOPE_LINEAR_1D', None)
        else:
            os.environ['XSLOPE_LINEAR_1D'] = old

    if len(types_now) and int(types_now.max()) != 2:
        fails.append(f"a tri3 mesh produced {int(types_now.max())}-node 1D "
                     f"elements")
    if K_now.shape != K_hatch.shape:
        fails.append(f"the stiffness changed shape on tri3: {K_now.shape} vs "
                     f"{K_hatch.shape}")
    else:
        diff = np.abs((K_now - K_hatch).data) if hasattr(K_now - K_hatch, "data") \
            else np.abs(np.asarray(K_now - K_hatch))
        if diff.size and float(diff.max()) != 0.0:
            fails.append(f"the tri3 stiffness moved by {float(diff.max()):.3e}, "
                         f"and should be bit-identical")
    if forces_now.shape != forces_hatch.shape:
        fails.append("the tri3 member force array changed length")
    elif forces_now.size and not np.array_equal(forces_now, forces_hatch):
        worst = float(np.max(np.abs(forces_now - forces_hatch)))
        fails.append(f"the tri3 member forces moved by {worst:.3e}, and should "
                     f"be bit-identical")
    return fails


# --------------------------------------------------------------------------- #

CHECKS = [
    ("the mesh attaches the midside node", test_mesh_attaches_the_midside_node),
    ("the escape hatch holds the old element",
     test_the_escape_hatch_holds_the_old_element),
    ("a bar under uniform axial strain", test_uniform_axial_strain),
    ("a linear mesh is untouched", test_a_linear_mesh_is_untouched),
]


def run():
    failures = []
    for label, fn in CHECKS:
        try:
            found = fn() or []
        except Exception as exc:                                # noqa: BLE001
            found = [f"{label}: raised {type(exc).__name__}: {exc}"]
        print(f"  {label:48s} {'FAIL' if found else 'ok'}")
        for message in found:
            print(f"      {message}")
        failures += found
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print("\nPASS: one-dimensional elements stand on every node of the soil "
          "edge they lie on.")


if __name__ == "__main__":
    main()
