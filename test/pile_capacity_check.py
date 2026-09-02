"""A pile's structural capacities, ENFORCED by the equilibrium rather than
reported on top of it.

Two things a capacity has to do, and neither can be checked by reading the
solver's own reported arrays, because those are what a broken cap gets right:

**A moment capacity is a release.** Two beam elements meet at every interior node
of a pile, and at equilibrium their end moments there are equal and opposite. A
correction applied to the rotational degree of freedom they SHARE is therefore
equal and opposite too and cancels exactly -- it moves nothing, at any capacity,
while the reported moments are clipped to the capacity and read as enforced. A
plastic hinge instead lets the element end rotate freely once its moment reaches
``M_cap``: the correction is the full element vector ``K_local p``, it has a
TRANSLATIONAL part that no nodal moment has, and the moment diagram is bounded by
the equilibrium.

**A capped member delivers the cap, not more.** The viscoplastic scheme solves
``K u = base_loads + corrections``, so the internal force the state is in
equilibrium with is ``K u - corrections``. The correction has to be
``action - action_true``, the bar's sign; the opposite sign is an anti-cap and the
member ends up delivering ``2 s - s_true``, growing STIFFER the tighter the cap is
drawn.

Both are read here off ``K_local u_local - correction_local`` -- the internal force,
never the reported array -- and compared against the capacity and against what the
solver reports. Cost: one coarse mesh of a shipped wall model and six fixed-F
solves, a few seconds each.
"""

import os
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

MODEL = _REPO / 'docs' / 'tutorials' / 'files' / 'xslope_pile_wall.xlsx'
TARGET_SIZE = 4.0        # coarse: this check is about the capacity law, not the mesh
F_TRIAL = 1.4            # below the model's own factor of safety, so uncapped converges


# ---------------------------------------------------------------- element level

def _element_law_failures():
    """The capacity law on one element, with no soil and no solve."""
    from xslope.fem import (_beam_local_stiffness, _pile_element_actions,
                            _pile_element_capacity)
    fails = []
    rng = np.random.default_rng(20260901)
    branch = {'elastic': 0, 'M': 0, 'V': 0, 'MV': 0}

    for _ in range(400):
        n_node = int(rng.integers(2, 4))
        L = float(rng.uniform(0.4, 6.0))
        EA = float(rng.uniform(1e4, 1e8))
        EI = float(rng.uniform(1e4, 1e8))
        K = _beam_local_stiffness(EA, EI, L, n_node == 3)
        u = rng.normal(size=3 * n_node) * 1e-4

        _, V_e, M1_e, M2_e, _ = _pile_element_actions(
            None, 1.0, 0.0, L, EA, EI, K, n_node, u_local=u)
        # Draw the caps from the demand this displacement actually produces, so
        # every branch is reached: no cap, moment only, shear only, both.
        m_ref = max(abs(M1_e), abs(M2_e))
        v_ref = abs(V_e)
        M_cap = float(rng.choice([np.inf, 10.0 * m_ref, 0.8 * m_ref, 0.3 * m_ref]))
        V_cap = float(rng.choice([np.inf, 10.0 * v_ref, 0.8 * v_ref, 0.3 * v_ref]))

        ax, V, M1, M2, corr, p_rot, yV, yM = _pile_element_capacity(
            u, 1.0, 0.0, L, EA, EI, K, n_node, V_cap, M_cap)
        branch['MV' if (yV and yM) else 'M' if yM else 'V' if yV else 'elastic'] += 1

        f = K @ u - corr           # the internal force this correction leaves

        # What the element reports IS what it delivers. The end moments are rows
        # of the internal force, so they compare on any displacement; the shear
        # the element reports is the midspan value of its own deflected shape,
        # which coincides with the pattern reading (f[1] - f[4]) / 2 only on the
        # two-node element, whose shear is constant along it. On the three-node
        # element the two agree once the displacement satisfies the assembled
        # equilibrium, which is where the model-level check reads them.
        pairs = [('M1', M1, f[2]), ('M2', M2, f[5])]
        if n_node == 2:
            pairs.append(('V', V, 0.5 * (f[1] - f[4])))
        f_scale = float(np.max(np.abs(f)))
        for name, got, want in pairs:
            scale = max(abs(want), f_scale, 1e-12)
            if abs(got - want) / scale > 1e-9:
                fails.append(f"element {name} reported {got:.6e} but delivers "
                             f"{want:.6e} (n_node={n_node})")

        # Nothing over its capacity.
        if max(abs(M1), abs(M2)) > M_cap * (1 + 1e-9):
            fails.append(f"end moment {max(abs(M1), abs(M2)):.6e} over M_cap "
                         f"{M_cap:.6e}")
        if abs(V) > V_cap * (1 + 1e-9):
            fails.append(f"shear {abs(V):.6e} over V_cap {V_cap:.6e}")

        if yV:
            if abs(abs(V) / V_cap - 1) > 1e-9:
                fails.append(f"shear-capped element carries {V:.6e}, not V_cap "
                             f"{V_cap:.6e}")
            # The correction REMOVES shear: on the shear's own pattern it acts in
            # the direction of the shear it is taking off, which is the bar's
            # sign. The opposite sign is the anti-cap, and it ADDS shear. Read
            # clear of the hinge's own contribution to the same rows.
            p_local = np.zeros(3 * n_node)
            p_local[2], p_local[5] = p_rot
            shear_leg = corr - K @ p_local
            sym = abs(shear_leg[1] + shear_leg[4]) > 1e-9 * max(
                abs(shear_leg[1]), float(np.max(np.abs(corr))))
            if shear_leg[1] * V <= 0 or sym:
                fails.append(f"the shear correction {shear_leg[1]:.6e} opposes "
                             f"the {V:.6e} it caps: the anti-cap sign, which "
                             f"makes the member deliver 2 V - V_cap")

        if yM:
            # A hinged end sits ON the capacity, and it got there by rotating.
            on_cap = [abs(abs(m) / M_cap - 1) < 1e-9 for m in (M1, M2)]
            if not any(on_cap):
                fails.append(f"hinged element carries {M1:.6e}/{M2:.6e}, neither "
                             f"on M_cap {M_cap:.6e}")
            if not np.any(np.abs(p_rot) > 0):
                fails.append("hinged element shows no plastic rotation")
            # THE structural point: a hinge is a release, so its correction acts on
            # the element's translational rows too. A moment applied at the shared
            # rotational degree of freedom has nothing there -- and cancels.
            trans = np.abs(corr[[i for i in range(3 * n_node) if i % 3 != 2]])
            rot = np.abs(corr[[2, 5]])
            if trans.max() <= 1e-12 * max(rot.max(), 1.0):
                fails.append("the moment correction is a bare nodal moment (no "
                             "translational part): it cancels between the two "
                             "elements that share the node and enforces nothing")
        else:
            if np.any(p_rot != 0.0):
                fails.append("plastic rotation on an element that did not hinge")

    missing = [k for k, v in branch.items() if v == 0]
    if missing:
        fails.append(f"capacity branches never exercised: {missing} ({branch})")
    return fails, branch


def _hinge_ownership_failures():
    """One release per pile NODE, not one per element end.

    A hinge is a single release at a single section, and every interior node of a
    pile carries two element ends. Release both and the node's rotation sees the
    two capacities, equal and opposite, whatever it does: no longer determined by
    the beam, and the initial-stiffness iteration drifts along that direction
    without ever settling — which reads as a failed trial and pulls the factor of
    safety below the unreinforced slope's.
    """
    from xslope.fem import _pile_hinge_ends, _pile_moment_hinge, _beam_local_stiffness
    fails = []

    # A chain of 4 elements on 5 nodes: 3 interior nodes, 2 ends.
    chain = np.array([[0, 1, -1], [1, 2, -1], [2, 3, -1], [3, 4, -1]])
    ends = _pile_hinge_ends(chain, 4)
    per_node = {}
    for p in range(4):
        for e in (0, 1):
            per_node.setdefault(int(chain[p][e]), []).append(bool(ends[p, e]))
    for node, flags in sorted(per_node.items()):
        if sum(flags) != 1:
            fails.append(f"node {node} carries {sum(flags)} releases, not 1 "
                         f"({flags})")

    # A disallowed end does not hinge, whatever its moment.
    K = _beam_local_stiffness(1e7, 1e7, 2.0, False)
    p_rot, m = _pile_moment_hinge(5e5, -5e5, 1e4, K, allowed=[True, False])
    if p_rot[1] != 0.0:
        fails.append("a disallowed end took a plastic rotation")
    if abs(abs(m[0]) / 1e4 - 1) > 1e-9:
        fails.append(f"the allowed end is at {m[0]:.6e}, not the capacity 1e4")
    return fails


# ------------------------------------------------------------------ model level

def _build():
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_point_constraints, extract_size_regions)
    slope_data = load_slope_data(str(MODEL))
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=lines)
    mesh = build_mesh_from_polygons(
        polygons, target_size=TARGET_SIZE, element_type='tri6', lines=lines,
        element_size_1d=slope_data.get('element_size_1d'),
        point_constraints=extract_point_constraints(slope_data),
        size_regions=extract_size_regions(slope_data))
    return build_fem_data(slope_data, mesh)


def _solve(fem_data, V_cap, M_cap):
    """Solve at F_TRIAL under the given caps and read the pile back off the
    internal force ``K_local u_local - correction_local``."""
    from xslope.fem import (solve_fem, _beam_rotation, _pile_element_capacity)
    n = fem_data['n_pile_elements']
    fem_data['V_cap_by_pile_elem'] = np.full(n, V_cap)
    fem_data['M_cap_by_pile_elem'] = np.full(n, M_cap)
    sol = solve_fem(fem_data, F=F_TRIAL, debug_level=0, fast_kernel=False,
                    max_iterations=16000)
    u = sol['displacements']

    rows = []
    for p in range(n):
        n_dof = int(fem_data['n_dof_by_pile_elem'][p])
        n_node = n_dof // 3
        dof = fem_data['dof_indices_pile'][p][:n_dof]
        cos_t = fem_data['cos_theta_pile'][p]
        sin_t = fem_data['sin_theta_pile'][p]
        L = fem_data['elem_length_by_pile_elem'][p]
        EA = fem_data['EA_by_pile_elem'][p]
        EI = fem_data['EI_by_pile_elem'][p]
        K = fem_data['K_local_by_pile_elem'][p]
        u_local = _beam_rotation(cos_t, sin_t, n_node) @ u[dof]
        _ax, V, M1, M2, corr, p_rot, yV, yM = _pile_element_capacity(
            u_local, cos_t, sin_t, L, EA, EI, K, n_node, V_cap, M_cap)
        f = K @ u_local - corr
        rows.append(dict(V=V, M1=M1, M2=M2,
                         V_del=0.5 * (f[1] - f[4]), M1_del=f[2], M2_del=f[5],
                         p_rot=p_rot, yV=yV, yM=yM))
    return sol, rows


def _maxes(rows):
    return (max(abs(r['V_del']) for r in rows),
            max(max(abs(r['M1_del']), abs(r['M2_del'])) for r in rows),
            max(np.max(np.abs(r['p_rot'])) for r in rows),
            sum(r['yM'] for r in rows), sum(r['yV'] for r in rows))


def _delivered_matches_reported(rows, tag):
    """On a solved state the pile's reported actions ARE the ones the equilibrium
    carries. Under the anti-cap sign a capped member reports its cap and delivers
    ``2 V - V_cap``, so this is what separates the two."""
    out = []
    for i, r in enumerate(rows):
        pairs = [('M1', r['M1'], r['M1_del']), ('M2', r['M2'], r['M2_del'])]
        if not r['yM']:
            # A hinged element's reported shear is the midspan value of its
            # RELEASED shape, which no longer satisfies the assembled equilibrium
            # the pattern reading is exact on; everywhere else the two coincide.
            pairs.append(('V', r['V'], r['V_del']))
        for name, got, want in pairs:
            scale = max(abs(want), 1e-6)
            if abs(got - want) / scale > 1e-6:
                out.append(f"{tag}: element {i} reports {name}={got:.6e} but the "
                           f"equilibrium carries {want:.6e}")
    return out


def _model_failures():
    fem_data = _build()
    if fem_data.get('n_pile_elements', 0) == 0:
        return [f"{MODEL.name} built no pile elements"], {}
    notes = {}
    fails = []

    # --- uncapped reference -------------------------------------------------
    sol0, rows0 = _solve(fem_data, np.inf, np.inf)
    V0, M0, prot0, nM0, nV0 = _maxes(rows0)
    u0 = float(np.max(np.abs(sol0['displacements'])))
    notes['uncapped'] = dict(V=V0, M=M0, maxu=u0, converged=bool(sol0['converged']))
    if not sol0['converged']:
        fails.append("the uncapped reference solve does not converge; the check's "
                     "trial F is no longer below this model's factor of safety")
    if nM0 or nV0 or prot0:
        fails.append("the uncapped solve reports a yielded pile element")
    fails += _delivered_matches_reported(rows0, 'uncapped')

    # --- a capacity above the demand changes nothing at all -----------------
    sol_hi, rows_hi = _solve(fem_data, np.inf, 10.0 * M0)
    u_hi = float(np.max(np.abs(sol_hi['displacements'])))
    if u_hi != u0:
        fails.append(f"M_cap above the elastic demand moved the answer: "
                     f"max|u| {u_hi!r} vs {u0!r}")

    # --- the binding case ---------------------------------------------------
    M_bind = 0.9 * M0
    sol_b, rows_b = _solve(fem_data, np.inf, M_bind)
    Vb, Mb, protb, nMb, nVb = _maxes(rows_b)
    u_b = float(np.max(np.abs(sol_b['displacements'])))
    notes['binding'] = dict(M_cap=M_bind, M_delivered=Mb, plastic_rotation=protb,
                            hinges=nMb, maxu=u_b, converged=bool(sol_b['converged']))
    if abs(Mb / M_bind - 1) > 1e-3:
        fails.append(f"the equilibrium's largest end moment is {Mb:.6f} against a "
                     f"capacity of {M_bind:.6f} ({abs(Mb / M_bind - 1):.2e} "
                     f"relative): the moment capacity is not being enforced")
    if nMb == 0 or protb <= 0.0:
        fails.append(f"no plastic hinge formed at M_cap = {M_bind:.6f} "
                     f"(hinges={nMb}, max plastic rotation={protb:.3e})")
    if not (u_b > u0):
        fails.append(f"releasing a hinge did not soften the wall: max|u| "
                     f"{u_b:.9f} against {u0:.9f} uncapped")
    fails += _delivered_matches_reported(rows_b, 'binding')

    # --- the shear capacity is a cap, not an anti-cap -----------------------
    ladder = []
    for frac in (0.9, 0.6, 0.4):
        cap = frac * V0
        sol_v, rows_v = _solve(fem_data, cap, np.inf)
        Vd, _Md, _pr, _nM, nV = _maxes(rows_v)
        ladder.append((cap, Vd, nV))
        fails += _delivered_matches_reported(rows_v, f"V_cap={cap:.1f}")
        if nV and abs(Vd / cap - 1) > 1e-6:
            fails.append(f"V_cap = {cap:.4f} binds but the pile delivers "
                         f"{Vd:.4f}")
    notes['shear_ladder'] = [(V0, None)] + [(c, v) for c, v, _ in ladder]
    seq = [V0] + [v for _c, v, _n in ladder]
    for a, b in zip(seq, seq[1:]):
        if b > a * (1 + 1e-9):
            fails.append(f"delivered shear ROSE as the capacity tightened: "
                         f"{a:.4f} -> {b:.4f} (the anti-cap sign)")
    return fails, notes


def run():
    os.environ.setdefault('MPLBACKEND', 'Agg')
    failures, branch = _element_law_failures()
    failures += _hinge_ownership_failures()
    model_failures, notes = _model_failures()
    return failures + model_failures


if __name__ == '__main__':
    el_f, br = _element_law_failures()
    el_f += _hinge_ownership_failures()
    print("element branches:", br)
    for f in el_f:
        print("  FAIL", f)
    mo_f, notes = _model_failures()
    for k, v in notes.items():
        print(k, "=", v)
    for f in mo_f:
        print("  FAIL", f)
    print("TOTAL FAILURES:", len(el_f) + len(mo_f))
