"""The pile beam element against closed-form beam theory.

XSLOPE models a pile as a two-node Euler-Bernoulli beam element: a 6x6 local
stiffness in [u1, v1, theta1, u2, v2, theta2] built at ``fem.py`` in
``build_fem_data``, rotated into global coordinates by a 6x6 transformation, and
assembled into the mixed-DOF global system (pile nodes carry three DOFs, soil
nodes two). The element's internal actions are recovered afterwards from the
local end displacements.

That element is the whole of XSLOPE's structural-member physics, so it is worth
pinning against solutions that do not depend on any of XSLOPE's own machinery.
The reference values are the closed forms quoted in GeoStudio's SIGMA/W
verification example "Beams and Bars in a Frame":

    simply supported, L = 8 m, E = 1e8 kPa, I = 0.01 m^4 (EI = 1e6 kN m^2)
        central point load P = 1000 kN   ->  M = PL/4      = 2000 kN m
                                             d = PL^3/48EI = 0.0106667 m
        uniform load w = 100 kN/m        ->  M = wL^2/8    = 800 kN m
                                             d = 5wL^4/384EI = 0.00533333 m
    cantilever, L = 3 m, same EI
        uniform load w = 100 kN/m        ->  M = wL^2/2    = 450 kN m
                                             V = wL        = 300 kN
                                             d = wL^4/8EI  = 0.0010125 m

HOW THE ELEMENT IS EXERCISED FREE-STANDING

``build_fem_data`` returns the assembled per-element beam stiffnesses
(``K_global_pile_elems``), their global DOF maps (``dof_indices_pile``) and the
section constants it derived (``EI_by_pile_elem``, ``EA_by_pile_elem``). So the
shipped element is obtained through the ordinary public path -- including the
1/S per-unit-width scaling and the rotation -- and is then assembled here into a
beam-only system. No soil stiffness enters, nothing is monkey-patched, and no
part of the element or its force recovery is reimplemented: the recovery step
below calls the same expressions the solver uses, quoted from ``fem.py``.

The mesh is hand-built rather than meshed by gmsh, so the nodes sit at exactly
the stations the closed forms are written for (a node at midspan, equal element
lengths) and every comparison is exact rather than mesh-dependent.

DISTRIBUTED LOAD

XSLOPE never applies an element-level distributed load to a pile: the soil loads
it through shared nodes. The two uniformly loaded cases are therefore applied as
the consistent (work-equivalent) nodal load vector, which is what a distributed
load becomes in any displacement-based beam formulation. Nodal displacements and
support reactions are then exact for any subdivision, and they are compared with
no correction of any kind. The element's recovered END ACTIONS are the equivalent
nodal-load problem's actions, so the true internal actions are recovered the
standard way -- by subtracting the element's own consistent load vector,
S = K u - f, a step that involves only the element length and the load
intensity. Both routes are checked, so the uniform-load benchmarks are reached
once without any correction (reactions) and once with it (end actions).

Run directly:  PYTHONPATH=. python3 test/beam_element_check.py
"""

import warnings

import numpy as np

warnings.filterwarnings('ignore')

from xslope.fem import build_fem_data

# --- the SIGMA/W "Beams and Bars in a Frame" section ------------------------ #
E_BEAM = 1.0e8          # kPa
I_BEAM = 0.01           # m^4
A_BEAM = 0.02           # m^2
EI = E_BEAM * I_BEAM    # 1e6 kN m^2

# Machine-precision tolerances: every comparison below is an exact identity of
# the discretization, not an engineering approximation, so the only thing between
# the element and the closed form is round-off on quantities of order 1e3.
RTOL = 1e-9


# --------------------------------------------------------------------------- #
# Building a beam-only model through the ordinary build_fem_data path
# --------------------------------------------------------------------------- #

def build_beam(length, n_elem, angle_deg=0.0, spacing=1.0, E=E_BEAM, I=I_BEAM,
               area=A_BEAM):
    """A straight chain of ``n_elem`` equal pile beam elements of total
    ``length``, laid at ``angle_deg`` to the x axis.

    Returns ``(fem_data, node_ids)`` where ``node_ids`` runs from the start of
    the line to its end, so station k sits at arc length k*length/n_elem.
    """
    s = length / n_elem
    t = np.radians(angle_deg)
    ax, ay = np.cos(t), np.sin(t)
    stations = np.arange(n_elem + 1) * s
    beam_nodes = np.column_stack([stations * ax, stations * ay])

    # One dummy 2D element, parked well away from the beam. build_fem_data wants
    # a 2D mesh to exist; its stiffness is never assembled here.
    soil = np.array([[1.0e3, 0.0], [1.0e3 + 1.0, 0.0], [1.0e3, 1.0]])
    nodes = np.vstack([soil, beam_nodes])
    elements = np.zeros((1, 6), dtype=int)
    elements[0, :3] = [0, 1, 2]

    n0 = len(soil)
    mesh = {
        'nodes': nodes,
        'elements': elements,
        'element_types': np.array([3]),
        'element_materials': np.array([1]),
        'elements_1d': np.array([[n0 + i, n0 + i + 1, 0] for i in range(n_elem)],
                                dtype=int),
        'element_types_1d': np.full(n_elem, 2, dtype=int),
        'element_materials_1d': np.ones(n_elem, dtype=int),
    }
    slope_data = {
        'materials': [{'name': 'soil', 'gamma': 20.0, 'E': 30000.0, 'nu': 0.3,
                       'option': 'mc', 'c': 25.0, 'phi': 30.0, 'u': 'none'}],
        'gamma_water': 9.81,
        'pile_lines': [{'x1': float(beam_nodes[0, 0]), 'y1': float(beam_nodes[0, 1]),
                        'x2': float(beam_nodes[-1, 0]), 'y2': float(beam_nodes[-1, 1]),
                        'E': E, 'I': I, 'area': area, 'S': spacing,
                        'head_fixity': 'free'}],
    }
    fem_data = build_fem_data(slope_data, mesh)
    return fem_data, list(range(n0, n0 + n_elem + 1))


def assemble(fem_data):
    """The beam-only global stiffness: the shipped per-element matrices
    scattered through the shipped DOF map. Nothing else is added.

    Returns ``(K, beam_dofs)`` -- the DOFs the beam chain actually reaches. The
    global vector is sized for the whole (mixed-DOF) model, so the dummy 2D
    element's DOFs are in it and carry no stiffness at all; they are held out of
    the solve rather than given a fictitious stiffness.
    """
    n_dof = int(fem_data['dof_offset'][-1])
    K = np.zeros((n_dof, n_dof))
    beam_dofs = set()
    for Ke, idx in zip(fem_data['K_global_pile_elems'], fem_data['dof_indices_pile']):
        K[np.ix_(idx, idx)] += Ke
        beam_dofs.update(int(i) for i in idx)
    return K, np.array(sorted(beam_dofs), dtype=int)


def solve(K, beam_dofs, f, fixed):
    """Solve K u = f over ``beam_dofs`` with the DOFs in ``fixed`` held at zero;
    return ``(u, reactions)`` where the reaction at a held DOF is (K u - f)."""
    free = np.setdiff1d(beam_dofs, np.asarray(fixed, dtype=int))
    u = np.zeros(K.shape[0])
    u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free])
    return u, K @ u - f


def dofs(fem_data, node, comp):
    """Global DOF index of component ``comp`` (0=ux, 1=uy, 2=theta) at ``node``."""
    return int(fem_data['dof_offset'][node]) + comp


# --------------------------------------------------------------------------- #
# Internal actions, recovered exactly as fem.py recovers them
# --------------------------------------------------------------------------- #

def end_actions(fem_data, u, p_idx):
    """(V, M1, M2, axial) for pile element ``p_idx`` -- the expressions in
    ``fem.py``'s "Compute final pile beam element forces" step, applied to the
    same local displacement vector."""
    idx = fem_data['dof_indices_pile'][p_idx]
    c = fem_data['cos_theta_pile'][p_idx]
    s = fem_data['sin_theta_pile'][p_idx]
    L = fem_data['elem_length_by_pile_elem'][p_idx]
    EI_val = fem_data['EI_by_pile_elem'][p_idx]
    EA_val = fem_data['EA_by_pile_elem'][p_idx]

    # The T matrix and the axial/V/M1/M2 expressions below are copies of
    # fem.py's Step 10c, "Compute final pile beam element forces"
    # (xslope/fem.py lines ~4845-4866), verified line-identical at this
    # commit. A change there must be mirrored here.
    T = np.array([
        [c, s, 0, 0, 0, 0],
        [-s, c, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0],
        [0, 0, 0, -s, c, 0],
        [0, 0, 0, 0, 0, 1],
    ])
    ul = T @ u[idx]

    L2, L3 = L * L, L * L * L
    axial = EA_val / L * (ul[3] - ul[0])
    V = 12 * EI_val / L3 * (ul[1] - ul[4]) + 6 * EI_val / L2 * (ul[2] + ul[5])
    M1 = EI_val * (6.0 / L2 * ul[1] + 4.0 / L * ul[2]
                   - 6.0 / L2 * ul[4] + 2.0 / L * ul[5])
    M2 = EI_val * (6.0 / L2 * ul[1] + 2.0 / L * ul[2]
                   - 6.0 / L2 * ul[4] + 4.0 / L * ul[5])
    return V, M1, M2, axial


def udl_consistent(fem_data, node_ids, w):
    """Global load vector for a uniform transverse load of intensity ``w``
    (acting along -local y) applied to every element of the chain.

    Per element of length s the work-equivalent nodal load is
    [-ws/2, -ws^2/12] at end 1 and [-ws/2, +ws^2/12] at end 2, rotated into
    global coordinates.
    """
    n_dof = int(fem_data['dof_offset'][-1])
    f = np.zeros(n_dof)
    for p_idx, idx in enumerate(fem_data['dof_indices_pile']):
        c = fem_data['cos_theta_pile'][p_idx]
        s_ = fem_data['sin_theta_pile'][p_idx]
        L = fem_data['elem_length_by_pile_elem'][p_idx]
        f_local = np.array([0.0, -w * L / 2.0, -w * L * L / 12.0,
                            0.0, -w * L / 2.0, +w * L * L / 12.0])
        T = np.array([
            [c, s_, 0, 0, 0, 0],
            [-s_, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s_, 0],
            [0, 0, 0, -s_, c, 0],
            [0, 0, 0, 0, 0, 1],
        ])
        f[idx] += T.T @ f_local
    return f


def udl_element_load(fem_data, p_idx, w):
    """The element's own consistent load vector in LOCAL coordinates, the term
    subtracted to turn recovered end actions into true internal actions."""
    L = fem_data['elem_length_by_pile_elem'][p_idx]
    return np.array([0.0, -w * L / 2.0, -w * L * L / 12.0,
                     0.0, -w * L / 2.0, +w * L * L / 12.0])


# --------------------------------------------------------------------------- #
# Comparison helper
# --------------------------------------------------------------------------- #

def cmp(failures, label, got, want, rtol=RTOL, unit='', atol=None):
    """Compare ``got`` to ``want``. A target of zero has no scale to be relative
    to, so those comparisons state an absolute tolerance instead (the quantity is
    normalised by its own reference before it gets here)."""
    if atol is not None:
        err, ok, how = abs(got - want), abs(got - want) <= atol, f"abs<={atol:g}"
    else:
        err = abs(got - want) / max(abs(want), 1e-30)
        ok, how = err <= rtol, f"rel<={rtol:g}"
    print(f"    {label:50s} {got:14.8f} vs {want:12.8f} {unit:5s} "
          f"err={err:.2e} {'ok' if ok else 'FAIL'}")
    if not ok:
        failures.append(f"{label}: {got:.10g} vs {want:.10g} ({err:.3e}, need {how})")


# --------------------------------------------------------------------------- #
# Leg 1 -- simply supported, central point load
# --------------------------------------------------------------------------- #

def leg_point_load():
    """L = 8 m, P = 1000 kN at midspan. Point loads only, so the cubic element
    is exact everywhere and the recovered end actions need no correction."""
    print("  simply supported, L=8 m, central P=1000 kN:")
    fails = []
    L, P, n = 8.0, 1000.0, 16
    fd, nid = build_beam(L, n)

    mid = nid[n // 2]
    K, bd = assemble(fd)
    f = np.zeros(int(fd['dof_offset'][-1]))
    f[dofs(fd, mid, 1)] = -P
    # pinned both ends (uy held; ux held at one end to remove the axial mode)
    fixed = [dofs(fd, nid[0], 1), dofs(fd, nid[-1], 1), dofs(fd, nid[0], 0)]
    u, R = solve(K, bd, f, fixed)

    cmp(fails, "midspan deflection  PL^3/48EI", -u[dofs(fd, mid, 1)],
        P * L ** 3 / (48 * EI), unit='m')
    cmp(fails, "support reaction    P/2", R[dofs(fd, nid[0], 1)], P / 2.0, unit='kN')

    # The element to the LEFT of midspan: its end-2 moment is the midspan moment,
    # its shear is the constant P/2 over that half-span.
    p_left = n // 2 - 1
    V, M1, M2, _ = end_actions(fd, u, p_left)
    cmp(fails, "midspan moment      PL/4", abs(M2), P * L / 4.0, unit='kN m')
    cmp(fails, "shear, left half    P/2", abs(V), P / 2.0, unit='kN')

    # ... and the element to the right carries the mirrored shear.
    V_r, _, _, _ = end_actions(fd, u, n // 2)
    cmp(fails, "shear, right half   -P/2", V_r, -abs(V), unit='kN')
    return fails


# --------------------------------------------------------------------------- #
# Leg 2 -- simply supported, uniform load
# --------------------------------------------------------------------------- #

def leg_udl_simply_supported():
    """L = 8 m, w = 100 kN/m."""
    print("  simply supported, L=8 m, w=100 kN/m:")
    fails = []
    L, w, n = 8.0, 100.0, 16
    fd, nid = build_beam(L, n)
    K, bd = assemble(fd)
    f = udl_consistent(fd, nid, w)
    fixed = [dofs(fd, nid[0], 1), dofs(fd, nid[-1], 1), dofs(fd, nid[0], 0)]
    u, R = solve(K, bd, f, fixed)

    mid = nid[n // 2]
    # No correction on either of these.
    cmp(fails, "midspan deflection  5wL^4/384EI", -u[dofs(fd, mid, 1)],
        5 * w * L ** 4 / (384 * EI), unit='m')
    cmp(fails, "support reaction    wL/2", R[dofs(fd, nid[0], 1)], w * L / 2.0, unit='kN')

    # True end actions: S = K u - f_element, the element's own consistent load
    # vector removed.
    p_left = n // 2 - 1
    V, M1, M2, _ = end_actions(fd, u, p_left)
    fe = udl_element_load(fd, p_left, w)
    cmp(fails, "midspan moment      wL^2/8", abs(M2 - fe[5]), w * L ** 2 / 8.0,
        unit='kN m')

    # End element: the moment at the support must vanish (a pinned end).
    V0, M1_0, _, _ = end_actions(fd, u, 0)
    fe0 = udl_element_load(fd, 0, w)
    cmp(fails, "moment at pinned end / (wL^2/8)",
        abs(M1_0 - fe0[2]) / (w * L ** 2 / 8.0), 0.0, atol=1e-12, unit='-')
    cmp(fails, "shear at support    wL/2", abs(V0 - fe0[1]), w * L / 2.0, unit='kN')
    return fails


# --------------------------------------------------------------------------- #
# Leg 3 -- cantilever, uniform load
# --------------------------------------------------------------------------- #

def leg_udl_cantilever():
    """L = 3 m, w = 100 kN/m, built in at the first node."""
    print("  cantilever, L=3 m, w=100 kN/m:")
    fails = []
    L, w, n = 3.0, 100.0, 12
    fd, nid = build_beam(L, n)
    K, bd = assemble(fd)
    f = udl_consistent(fd, nid, w)
    root = nid[0]
    fixed = [dofs(fd, root, 0), dofs(fd, root, 1), dofs(fd, root, 2)]
    u, R = solve(K, bd, f, fixed)

    cmp(fails, "tip deflection      wL^4/8EI", -u[dofs(fd, nid[-1], 1)],
        w * L ** 4 / (8 * EI), unit='m')
    # Support reactions -- exact, and reached with no correction at all.
    cmp(fails, "base reaction shear wL", R[dofs(fd, root, 1)], w * L, unit='kN')
    cmp(fails, "base reaction moment wL^2/2", abs(R[dofs(fd, root, 2)]),
        w * L ** 2 / 2.0, unit='kN m')

    # The same two numbers as element end actions.
    V, M1, M2, _ = end_actions(fd, u, 0)
    fe = udl_element_load(fd, 0, w)
    cmp(fails, "fixed-end moment    wL^2/2", abs(M1 - fe[2]), w * L ** 2 / 2.0,
        unit='kN m')
    cmp(fails, "fixed-end shear     wL", abs(V - fe[1]), w * L, unit='kN')

    # The tip is free: zero shear and zero moment there.
    V_t, _, M2_t, _ = end_actions(fd, u, n - 1)
    fe_t = udl_element_load(fd, n - 1, w)
    cmp(fails, "moment at free tip / (wL^2/2)",
        abs(M2_t - fe_t[5]) / (w * L ** 2 / 2.0), 0.0, atol=1e-12, unit='-')
    return fails


# --------------------------------------------------------------------------- #
# Leg 4 -- orientation invariance (the rotation matrix)
# --------------------------------------------------------------------------- #

def leg_orientation():
    """The same simply supported beam laid at several angles must give the same
    answer in its own local frame -- the only thing that changes is T."""
    print("  orientation invariance (cantilever, w=100 kN/m normal to the axis):")
    fails = []
    L, w, n = 3.0, 100.0, 12
    for angle in (0.0, 37.0, 90.0, 143.0, -65.0):
        fd, nid = build_beam(L, n, angle_deg=angle)
        K, bd = assemble(fd)
        f = udl_consistent(fd, nid, w)          # already applied along -local y
        root = nid[0]
        u, R = solve(K, bd, f,
                     [dofs(fd, root, 0), dofs(fd, root, 1), dofs(fd, root, 2)])

        # tip deflection measured in the local transverse direction
        c = fd['cos_theta_pile'][0]
        s = fd['sin_theta_pile'][0]
        tip = nid[-1]
        v_local = -s * u[dofs(fd, tip, 0)] + c * u[dofs(fd, tip, 1)]
        V, M1, _, _ = end_actions(fd, u, 0)
        fe = udl_element_load(fd, 0, w)
        tag = f"[{angle:+6.1f} deg]"
        cmp(fails, f"{tag} tip deflection  wL^4/8EI", -v_local,
            w * L ** 4 / (8 * EI), unit='m')
        cmp(fails, f"{tag} fixed-end moment wL^2/2", abs(M1 - fe[2]),
            w * L ** 2 / 2.0, unit='kN m')
        cmp(fails, f"{tag} fixed-end shear  wL", abs(V - fe[1]), w * L, unit='kN')
    return fails


# --------------------------------------------------------------------------- #
# Leg 5 -- axial stiffness and the 1/S per-unit-width scaling
# --------------------------------------------------------------------------- #

def leg_axial_and_spacing():
    """EA/L bar action, and the rule that a pile spacing S divides both EA and
    EI so the 2D model carries per-unit-width stiffness."""
    print("  axial stiffness and 1/S per-unit-width scaling:")
    fails = []
    L, N, n = 8.0, 1000.0, 16
    fd, nid = build_beam(L, n)
    K, bd = assemble(fd)
    f = np.zeros(int(fd['dof_offset'][-1]))
    f[dofs(fd, nid[-1], 0)] = N
    u, _ = solve(K, bd, f, [dofs(fd, nid[0], 0), dofs(fd, nid[0], 1),
                            dofs(fd, nid[-1], 1)])
    cmp(fails, "axial elongation    NL/EA", u[dofs(fd, nid[-1], 0)],
        N * L / (E_BEAM * A_BEAM), unit='m')
    _, _, _, axial = end_actions(fd, u, 0)
    cmp(fails, "axial force         N", axial, N, unit='kN')

    # Spacing: S = 2.5 must divide EI and EA exactly.
    S = 2.5
    fd2, nid2 = build_beam(L, n, spacing=S)
    cmp(fails, f"EI at S={S}          EI/S", fd2['EI_by_pile_elem'][0], EI / S,
        unit='kN m2')
    cmp(fails, f"EA at S={S}          EA/S", fd2['EA_by_pile_elem'][0],
        E_BEAM * A_BEAM / S, unit='kN')

    # ... so the same load produces S times the deflection.
    K2, bd2 = assemble(fd2)
    f2 = udl_consistent(fd2, nid2, 100.0)
    u2, _ = solve(K2, bd2, f2, [dofs(fd2, nid2[0], 1), dofs(fd2, nid2[-1], 1),
                                dofs(fd2, nid2[0], 0)])
    fd1, nid1 = build_beam(L, n, spacing=1.0)
    K1, bd1 = assemble(fd1)
    f1 = udl_consistent(fd1, nid1, 100.0)
    u1, _ = solve(K1, bd1, f1, [dofs(fd1, nid1[0], 1), dofs(fd1, nid1[-1], 1),
                                dofs(fd1, nid1[0], 0)])
    mid2 = u2[dofs(fd2, nid2[n // 2], 1)]
    mid1 = u1[dofs(fd1, nid1[n // 2], 1)]
    cmp(fails, f"deflection ratio at S={S}", mid2 / mid1, S, unit='-')
    return fails


# --------------------------------------------------------------------------- #
# Leg 6 -- section constants derived from a diameter
# --------------------------------------------------------------------------- #

def leg_circular_section():
    """A pile stated by diameter alone: I = pi D^4/64, A = pi D^2/4."""
    print("  circular section derived from the diameter:")
    fails = []
    D = 0.8
    s = 8.0 / 16
    t = np.radians(0.0)
    stations = np.arange(17) * s
    beam_nodes = np.column_stack([stations, np.zeros(17)])
    soil = np.array([[1.0e3, 0.0], [1.0e3 + 1.0, 0.0], [1.0e3, 1.0]])
    nodes = np.vstack([soil, beam_nodes])
    elements = np.zeros((1, 6), dtype=int)
    elements[0, :3] = [0, 1, 2]
    mesh = {
        'nodes': nodes, 'elements': elements,
        'element_types': np.array([3]), 'element_materials': np.array([1]),
        'elements_1d': np.array([[3 + i, 4 + i, 0] for i in range(16)], dtype=int),
        'element_types_1d': np.full(16, 2, dtype=int),
        'element_materials_1d': np.ones(16, dtype=int),
    }
    slope_data = {
        'materials': [{'name': 'soil', 'gamma': 20.0, 'E': 30000.0, 'nu': 0.3,
                       'option': 'mc', 'c': 25.0, 'phi': 30.0, 'u': 'none'}],
        'gamma_water': 9.81,
        'pile_lines': [{'x1': 0.0, 'y1': 0.0, 'x2': 8.0, 'y2': 0.0,
                        'E': E_BEAM, 'D_pile': D, 'S': 1.0, 'head_fixity': 'free'}],
    }
    fd = build_fem_data(slope_data, mesh)
    cmp(fails, "I from D            pi D^4/64", fd['EI_by_pile_elem'][0] / E_BEAM,
        np.pi * D ** 4 / 64.0, unit='m4')
    cmp(fails, "A from D            pi D^2/4", fd['EA_by_pile_elem'][0] / E_BEAM,
        np.pi * D ** 2 / 4.0, unit='m2')
    return fails


# --------------------------------------------------------------------------- #

LEGS = [
    ("simply supported, central point load", leg_point_load),
    ("simply supported, uniform load", leg_udl_simply_supported),
    ("cantilever, uniform load", leg_udl_cantilever),
    ("orientation invariance", leg_orientation),
    ("axial stiffness and spacing", leg_axial_and_spacing),
    ("circular section from diameter", leg_circular_section),
]


def run():
    failures = []
    for label, fn in LEGS:
        failures += fn()
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print("\nPASS: the pile beam element reproduces closed-form beam theory.")


if __name__ == "__main__":
    main()
