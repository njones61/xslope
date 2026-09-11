"""Zero-thickness interface (joint) elements and the end ties that anchor a bar.

A jointed constraint line is a slip surface: the mesher splits the mesh along it
(``mesh.split_mesh_along_joints``) so the two faces can move relative to one
another, and joint elements carry the traction between them. This module turns
that mesh into the arrays the finite element solve reads, and holds the
interface's kinematics, constitutive law and element stiffness.

A line comes in two kinds. A REINFORCEMENT line whose ``Joint`` column reads yes
has a sheet between its faces, so every station carries three coincident nodes —
the soil above, the bar, the soil below — and a pair of joint elements spans
them, an upper one between the soil above and the bar and a lower one between the
bar and the soil below. A line off the ``joints`` sheet — a rock joint, a bedding
plane, a block-on-block contact — has nothing between its faces, so a station
carries two nodes and ONE joint element spans them. The two kinds share
everything below; the pair carries its stiffness as two springs in series, the
single element carries it once.

**Kinematics.** A joint element has no thickness and no area: its state is the
RELATIVE displacement of the two faces it connects. Side ``a`` of
``elements_joint`` is always the upper of the pair, so with

    t = the element chord, from its first node to its second,
    n = t rotated a quarter turn counter-clockwise (from side b to side a),

the relative displacement ``d = u_a - u_b`` resolves into

    delta_t = t . d          (sliding)
    delta_n = -n . d         (closing; COMPRESSION POSITIVE)

**Constitutive law.** Elastic tractions ``t_n = k_n delta_n`` and
``t_s = k_s delta_t``, with a Mohr-Coulomb limit on the shear

    |t_s| <= c_j + t_n tan phi_j     (t_n >= 0)

and a tension cutoff on the normal: at ``t_n < -t_cut`` the joint OPENS — both
tractions and both stiffnesses go to zero — and it closes again when the faces
return to contact. Slip past the limit is perfectly plastic: the shear traction
stays at the limit while the tangential offset grows.

**Integration.** The tractions are integrated at the element's own NODES
(Newton-Cotes / Lobatto), not at Gauss points: L/6, L/6, 2L/3 on the three-pair
element that matches a tri6 edge and L/2, L/2 on the two-pair element that
matches a tri3 edge. Nodal integration keeps the node pairs uncoupled, which is
what keeps the traction profile along a stiff interface free of the oscillation
Gauss quadrature produces there (Schellekens & de Borst 1993).

**Stiffness.** ``k_n`` and ``k_s`` are penalty-like: large enough that an intact
joint does not visibly deform, small enough not to ill-condition the system.
Blank on the line, they are derived as ``E_adj / d_v`` and ``G_adj / d_v`` over a
virtual thickness ``d_v = 0.1 x`` the 1D element length, with ``E_adj`` and
``G_adj`` those of the SOFTER of the two soils the joint element stands between —
which differ where the line crosses a material boundary.

**Strength.** ``c_j`` and ``tan phi_j`` are the line's ``Adhesion`` and ``Delta``,
the same two numbers the overburden-dependent pullout law reads
(``fileio.reinforce_pullout_profile``), in the same per-unit-width convention:
that law divides its rate by the line's ``Spacing``, and so do these, so a
discrete support's interface strength and its pullout envelope stay in one
convention. A continuous sheet — the only thing a mesh split represents — has a
blank Spacing and the division is inert. Its tension cutoff is zero: a
soil-geotextile interface carries no tension. A ``joints`` sheet line states
``c``, ``phi`` and ``t_cut`` in columns of its own, per unit width, with no
spacing to divide by.

**Tips.** The two soil faces rejoin at each end of the line, so the end station
carries one soil node, not two, and the two joint elements there stand between
the SAME pair of nodes with opposite orientation. Their shear contributions add
(the sheet still has two faces at its tip) but their normal tractions are equal
and opposite, so neither measures the confining stress: a relative normal
displacement between the two soil faces does not exist where they are one node.
A tip pair therefore never opens, and the normal traction its yield limit is
taken at is read from the rest of its own element.
"""

import warnings

import numpy as np

#: Virtual thickness for the derived joint stiffnesses, as a fraction of the 1D
#: element length (PLAXIS's idiom): k_n = E_adj / d_v, k_s = G_adj / d_v.
JOINT_VIRTUAL_THICKNESS_FRAC = 0.1

#: The iteration budget a jointed model needs before a strength-reduction trial
#: can be trusted to have decided.
#:
#: A joint reaches equilibrium by growing slip, and the slip a viscoplastic sweep
#: puts into an interface is the excess traction divided by k_s — so a jointed
#: model converges by tens of thousands of sweeps where a bonded one converges by
#: hundreds. A trial that runs out of budget is recorded undecided, the bracket
#: reads it as not standing, and the factor of safety comes out LOW: it is then a
#: statement about the budget rather than about the slope.
#:
#: 100 000 is measured, not assumed. On the geotextile wall family the longest
#: trial that reached a verdict was 96 738 sweeps (the weak-foundation variant's
#: equilibrium just below its critical factor), and at the 50 000 the solver
#: reaches by default four of that row's nine trials never decided and the factor
#: came out ten percent low. It costs almost nothing to allow, because a trial
#: that decides stops.
JOINT_DECIDED_BUDGET = 100000

#: The viscoplastic pseudo-time step of the joint's slip increment. At 1.0 one
#: sweep returns the shear traction exactly to its limit at the current
#: displacement field, which is the interface twin of the bar's tension cap.
JOINT_VP_DT = 1.0

#: Newton-Cotes (Lobatto) weights on [0, 1] in the node order the mesh writes,
#: (start, end, midside): Simpson's rule for the three-pair element and the
#: trapezoidal rule for the two-pair one.
_LOBATTO_3 = (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0)
_LOBATTO_2 = (0.5, 0.5, 0.0)


def mesh_has_joints(mesh):
    """True when this mesh was split along at least one jointed line."""
    conn = None if mesh is None else mesh.get("elements_joint")
    return conn is not None and len(conn) > 0


def solution_has_joint_state(solution, n):
    """True when a solution carries the interface state of ``n`` joint elements.

    A field read back from a saved sidecar carries the soil, the bars and the
    piles and nothing about the joints, so its joint arrays are absent. Reading
    an absent array as zeros would report every interface intact — a state
    nothing measured — so every consumer asks this first and draws, tabulates
    and lists nothing where it is false.
    """
    ts = (solution or {}).get("joint_ts")
    if ts is None:
        return False
    ts = np.asarray(ts)
    return ts.ndim == 2 and ts.shape[0] == int(n) and ts.shape[1] == 3


#: ``element_side_joint`` on a joint element that is the whole interface — one
#: element between the two faces, on a line with no bar. The mesher's own name
#: for it is :data:`xslope.mesh.JOINT_SIDE_WHOLE`; kept here as a plain number so
#: the element module does not import the mesher to read a mesh key.
_SIDE_WHOLE = 2


def _node_element_map(elements, element_types, n_nodes):
    """node id -> list of 2D element indices standing on it."""
    out = [[] for _ in range(n_nodes)]
    for ei in range(len(elements)):
        for k in range(int(element_types[ei])):
            out[int(elements[ei, k])].append(ei)
    return out


def _dof_pair(node, dof_offset):
    base = int(dof_offset[node]) if dof_offset is not None else 2 * int(node)
    return base, base + 1


def build_joint_data(slope_data, mesh, nodes, E_by_mat, nu_by_mat,
                     elements, element_types, element_materials,
                     dof_offset=None):
    """Per-joint-element arrays for the solve, or ``None`` when the mesh has none.

    Reads the split the mesher wrote (``elements_joint``, ``element_types_joint``,
    ``element_materials_joint``, ``element_side_joint``, ``joints``, ``ties``) and
    the jointed lines' own properties — from ``slope_data['reinforcement_lines']``
    for a line whose ``Joint`` column reads yes, whose interface strength is its
    Adhesion and Delta, and from ``slope_data['joint_lines']`` for a line off the
    ``joints`` sheet, which states c, phi and a tension cutoff of its own — and
    returns the dictionary ``fem_data['joint_data']`` carries:

    ``n``            number of joint elements
    ``conn``         (n, 6) node ids, [a0, a1, a2, b0, b1, b2]; side a is the upper
    ``n_pairs``      (n,) 3 on a quadratic mesh, 2 on a linear one
    ``line_id``      (n,) 1-based constraint-line index
    ``side``         (n,) 1 on the upper joint of a bar's pair, 0 on the lower,
                     2 on a bar-less line's single element, which IS the interface
    ``dof``          (n, 12) global degrees of freedom in the ``conn`` order
    ``w``            (n, 3) Lobatto weights, LENGTHS (zero on a padded pair)
    ``tx, ty``       (n,) the unit chord
    ``nx, ny``       (n,) its counter-clockwise normal, running from b to a
    ``L``            (n,) element length
    ``kn, ks``       (n,) normal and shear stiffness
    ``cj, tanphi``   (n,) interface cohesion and friction coefficient
    ``tcut``         (n,) tension cutoff — the joints sheet's own column; zero on
                     a reinforcement line, whose interface carries no tension
    ``jred``         (n,) reduce this joint's strength in the SSR
    ``tip``          (n, 3) pairs standing at an end station, where the two soil
                     faces are one node
    ``K``            (n, 12, 12) elastic element stiffness
    ``jointed_lines`` the 1-based line indices that are jointed

    and, when the model states any tied end, a ``ties`` sub-dictionary with
    ``dof`` (m, 4) in the order (bar node, soil node), ``k`` (m,), ``cap`` (m,)
    and ``K`` (m, 4, 4).
    """
    if not mesh_has_joints(mesh):
        return None
    conn = np.asarray(mesh["elements_joint"], dtype=int).reshape(-1, 6)
    n = len(conn)
    if n == 0:
        return None
    n_pairs = np.asarray(mesh["element_types_joint"], dtype=int)
    line_id = np.asarray(mesh["element_materials_joint"], dtype=int)
    side = np.asarray(mesh["element_side_joint"], dtype=int)
    nodes = np.asarray(nodes, dtype=float)

    lines = slope_data.get("reinforcement_lines") or []
    # The constraint-line numbering the mesh keys off: reinforcement lines, then
    # piles, then the joints sheet's own lines. A joint element's line id lands in
    # the first block or the third; the middle one is the piles, which are never
    # jointed.
    n_reinf = len(lines)
    n_pile = len(slope_data.get("pile_lines") or [])
    sheet_joints = slope_data.get("joint_lines") or []

    # ---- geometry: the chord, its normal, and the nodal integration weights ----
    p0 = nodes[conn[:, 0], :2]
    p1 = nodes[conn[:, 1], :2]
    dvec = p1 - p0
    L = np.hypot(dvec[:, 0], dvec[:, 1])
    if np.any(L <= 0.0):
        bad = int(np.flatnonzero(L <= 0.0)[0])
        raise ValueError(
            f"Joint element {bad} on constraint line {int(line_id[bad])} has zero "
            "length; the mesh split produced a degenerate interface element.")
    tx, ty = dvec[:, 0] / L, dvec[:, 1] / L
    nx, ny = -ty, tx

    w = np.zeros((n, 3))
    w[n_pairs >= 3] = _LOBATTO_3
    w[n_pairs < 3] = _LOBATTO_2
    w *= L[:, None]

    dof = np.zeros((n, 12), dtype=int)
    for i in range(n):
        row = []
        for k in range(6):
            row.extend(_dof_pair(conn[i, k], dof_offset))
        dof[i] = row

    # ---- which pairs stand at a tip, where the two soil faces are one node ----
    tip_stations = set()
    for rec in mesh.get("joints") or []:
        for st in rec.get("stations", ()):
            if int(st[0]) == int(st[2]):
                tip_stations.add(int(st[0]))
    tip = np.zeros((n, 3), dtype=bool)
    if tip_stations:
        for p in range(3):
            a_ok = np.array([int(v) in tip_stations for v in conn[:, p]])
            b_ok = np.array([int(v) in tip_stations for v in conn[:, 3 + p]])
            tip[:, p] = a_ok | b_ok
        tip &= (w > 0.0)

    # ---- the adjacent soils, per joint element: the two sides of the line ----
    node_elems = _node_element_map(elements, element_types, len(nodes))
    E_by_mat = np.asarray(E_by_mat, dtype=float)
    nu_by_mat = np.asarray(nu_by_mat, dtype=float)
    G_by_mat = E_by_mat / (2.0 * (1.0 + nu_by_mat))

    # The upper and lower joint of a station span share the bar's nodes, so the
    # two are paired on that node set and each reads BOTH soils. A bar-less joint
    # element (side WHOLE) is the whole interface on its own and has no partner:
    # both of its sides are soil, and it reads them directly.
    bar_key = {}
    for i in range(n):
        if side[i] == _SIDE_WHOLE:
            continue
        cols = (3, 4, 5) if side[i] == 1 else (0, 1, 2)
        bar_key.setdefault(tuple(int(conn[i, c]) for c in cols), []).append(i)
    partner = np.full(n, -1, dtype=int)
    for members in bar_key.values():
        if len(members) == 2:
            partner[members[0]], partner[members[1]] = members[1], members[0]

    E_adj = np.zeros(n)
    G_adj = np.zeros(n)
    for i in range(n):
        if side[i] == _SIDE_WHOLE:
            soil_nodes = ([int(conn[i, c]) for c in range(int(n_pairs[i]))]
                          + [int(conn[i, 3 + c]) for c in range(int(n_pairs[i]))])
            j = -1
        else:
            soil_cols = (0, 1, 2) if side[i] == 1 else (3, 4, 5)
            soil_nodes = [int(conn[i, c]) for c in soil_cols[:int(n_pairs[i])]]
            j = partner[i]
        if j >= 0:
            cols_j = (0, 1, 2) if side[j] == 1 else (3, 4, 5)
            soil_nodes += [int(conn[j, c]) for c in cols_j[:int(n_pairs[j])]]
        mats = set()
        for nd in soil_nodes:
            for ei in node_elems[nd]:
                mats.add(int(element_materials[ei]) - 1)
        if not mats:
            raise ValueError(
                f"Joint element {i} on constraint line {int(line_id[i])} stands on "
                "no 2D element, so no adjacent soil stiffness can be read.")
        mats = sorted(mats)
        E_adj[i] = float(np.min(E_by_mat[mats]))
        G_adj[i] = float(np.min(G_by_mat[mats]))

    d_v = JOINT_VIRTUAL_THICKNESS_FRAC * L
    kn = E_adj / d_v
    ks = G_adj / d_v

    # ---- the line's own properties ----
    cj = np.zeros(n)
    tanphi = np.zeros(n)
    tcut = np.zeros(n)
    jred = np.ones(n, dtype=bool)
    for li in sorted(set(int(v) for v in line_id)):
        sel = line_id == li
        jsheet = li - n_reinf - n_pile - 1     # index into the joints sheet
        if 0 <= jsheet < len(sheet_joints):
            # A line off the joints sheet states its own strength: c, phi and the
            # tension cutoff are columns of its own, and there is no out-of-plane
            # spacing to divide by — a joint is a surface, not a discrete member.
            jl = sheet_joints[jsheet]
            _c = jl.get("c")
            _phi = jl.get("phi")
            _tc = jl.get("t_cut")
            cj[sel] = 0.0 if _c is None or not np.isfinite(float(_c)) else float(_c)
            tanphi[sel] = (0.0 if _phi is None or not np.isfinite(float(_phi))
                           else np.tan(np.radians(float(_phi))))
            tcut[sel] = (0.0 if _tc is None or not np.isfinite(float(_tc))
                         else float(_tc))
            for key, arr in (("kn", kn), ("ks", ks)):
                v = jl.get(key)
                if v is not None and np.isfinite(float(v)) and float(v) > 0.0:
                    arr[sel] = float(v)
            jred[sel] = str(jl.get("jred", "yes") or "yes").strip().lower() != "no"
            continue
        line = lines[li - 1] if 0 < li <= len(lines) else {}
        label = line.get("label") or f"line {li}"
        # Adhesion and Delta ARE the interface strength. A line that states
        # neither describes a frictionless, cohesionless interface, which is a
        # model the element can carry and preflight is the place that refuses;
        # the warning is here because a silent zero is the reading that would
        # otherwise go unnoticed.
        adhesion = line.get("adhesion")
        delta = line.get("delta")
        adhesion = 0.0 if adhesion is None or not np.isfinite(float(adhesion)) \
            else float(adhesion)
        delta = 0.0 if delta is None or not np.isfinite(float(delta)) \
            else float(delta)
        if adhesion == 0.0 and delta == 0.0:
            warnings.warn(
                f"Reinforcement {label!r} is flagged as a joint, so the mesh is "
                f"split along it and the two faces slide on the interface "
                f"strength stated in the Adhesion and Delta columns. Both are "
                f"blank, so this interface has no strength at all.")
        spacing = float(line.get("spacing") or 1.0)
        if not np.isfinite(spacing) or spacing <= 0.0:
            spacing = 1.0
        cj[sel] = adhesion / spacing
        tanphi[sel] = np.tan(np.radians(delta)) / spacing
        _kn = line.get("kn")
        _ks = line.get("ks")
        if _kn is not None and np.isfinite(float(_kn)) and float(_kn) > 0.0:
            kn[sel] = float(_kn)
        if _ks is not None and np.isfinite(float(_ks)) and float(_ks) > 0.0:
            ks[sel] = float(_ks)
        jred[sel] = str(line.get("jred", "yes") or "yes").strip().lower() != "no"

    K = _joint_element_stiffness(w, tx, ty, nx, ny, kn, ks)

    jd = {
        "n": int(n), "conn": conn, "n_pairs": n_pairs, "line_id": line_id,
        "side": side, "dof": dof, "w": w,
        "tx": tx, "ty": ty, "nx": nx, "ny": ny, "L": L,
        "kn": kn, "ks": ks, "cj": cj, "tanphi": tanphi, "tcut": tcut,
        "jred": jred, "tip": tip, "K": K,
        "jointed_lines": sorted(set(int(v) for v in line_id)),
    }
    ties = _build_tie_data(mesh, lines, nodes, dof_offset)
    if ties is not None:
        jd["ties"] = ties
    return jd


def _joint_element_stiffness(w, tx, ty, nx, ny, kn, ks):
    """The elastic (12, 12) element matrices, block-diagonal over the node pairs.

    Nodal integration leaves the pairs uncoupled, so pair ``p`` contributes only
    the 4x4 block ``[[wD, -wD], [-wD, wD]]`` on its own two nodes, with
    ``D = k_s t (x) t + k_n n (x) n``.
    """
    n = len(w)
    K = np.zeros((n, 12, 12))
    for p in range(3):
        wp = w[:, p]
        D = np.empty((n, 2, 2))
        D[:, 0, 0] = wp * (ks * tx * tx + kn * nx * nx)
        D[:, 0, 1] = wp * (ks * tx * ty + kn * nx * ny)
        D[:, 1, 0] = D[:, 0, 1]
        D[:, 1, 1] = wp * (ks * ty * ty + kn * ny * ny)
        ia, ib = 2 * p, 6 + 2 * p
        K[:, ia:ia + 2, ia:ia + 2] += D
        K[:, ia:ia + 2, ib:ib + 2] -= D
        K[:, ib:ib + 2, ia:ia + 2] -= D
        K[:, ib:ib + 2, ib:ib + 2] += D
    return K


def _build_tie_data(mesh, lines, nodes, dof_offset):
    """The tied ends: a spring from the bar's end node to the soil node there.

    A jointed line's end is FREE unless the line states an end anchorage, in
    which case the bar's end node is tied to the soil (or facing) node at the
    same point by a spring of stiffness ``EA`` over one 1D element length —
    rigid next to the bar — that is perfectly plastic at the stated ``Tend``.
    """
    recs = mesh.get("ties") or []
    if not recs:
        return None
    e1d = np.asarray(mesh.get("elements_1d", np.zeros((0, 3))), dtype=int)
    mats_1d = np.asarray(mesh.get("element_materials_1d", np.zeros(len(e1d))), dtype=int)
    dofs, kk, cap = [], [], []
    for rec in recs:
        li = int(rec["line"])
        line = lines[li - 1] if 0 < li <= len(lines) else {}
        E = float(line.get("E") or 0.0)
        A = float(line.get("area") or 0.0)
        sel = np.flatnonzero(mats_1d == li)
        if sel.size == 0 or E <= 0.0 or A <= 0.0:
            continue
        p0 = nodes[e1d[sel, 0], :2]
        p1 = nodes[e1d[sel, 1], :2]
        le = float(np.mean(np.hypot(p1[:, 0] - p0[:, 0], p1[:, 1] - p0[:, 1])))
        if le <= 0.0:
            continue
        bar_dof = _dof_pair(int(rec["bar_node"]), dof_offset)
        soil_dof = _dof_pair(int(rec["soil_node"]), dof_offset)
        dofs.append([bar_dof[0], bar_dof[1], soil_dof[0], soil_dof[1]])
        kk.append(E * A / le)
        cap.append(float(rec["capacity"]))
    if not dofs:
        return None
    m = len(dofs)
    k = np.asarray(kk, dtype=float)
    K = np.zeros((m, 4, 4))
    for d in range(2):
        K[:, d, d] += k
        K[:, d, 2 + d] -= k
        K[:, 2 + d, d] -= k
        K[:, 2 + d, 2 + d] += k
    return {"n": m, "dof": np.asarray(dofs, dtype=int), "k": k,
            "cap": np.asarray(cap, dtype=float), "K": K}


def joint_reduced_strength(jd, F):
    """``(c_j, tan phi_j)`` after the strength reduction, per joint element.

    Both are divided by the trial factor on every joint whose line does not say
    otherwise; ``k_n``, ``k_s``, the ties and the bar's own capacity are
    structural and are not reduced, the same rule the reinforcement follows.
    """
    F = float(F)
    if F == 1.0:
        return jd["cj"].copy(), jd["tanphi"].copy()
    div = np.where(jd["jred"], F, 1.0)
    return jd["cj"] / div, jd["tanphi"] / div


def joint_kinematics(jd, u):
    """The relative displacement of the two faces at every node pair.

    Returns ``(delta_t, delta_n)``, each (n, 3), with ``delta_n`` COMPRESSION
    POSITIVE. A padded pair on a two-pair element holds one node against itself
    and reads zero.
    """
    dof = jd["dof"]
    ue = u[dof]                                     # (n, 12)
    dx = ue[:, 0:6:2] - ue[:, 6:12:2]               # (n, 3)
    dy = ue[:, 1:6:2] - ue[:, 7:12:2]
    dt = dx * jd["tx"][:, None] + dy * jd["ty"][:, None]
    dn = -(dx * jd["nx"][:, None] + dy * jd["ny"][:, None])
    return dt, dn


def _limit_normal(jd, tn):
    """The normal traction the yield limit is taken at, tips borrowing from their
    own element (see the module docstring)."""
    tip = jd["tip"]
    if not tip.any():
        return tn
    out = tn.copy()
    body = (jd["w"] > 0.0) & ~tip
    cnt = body.sum(axis=1)
    have = cnt > 0
    mean = np.zeros(len(tn))
    mean[have] = (tn * body).sum(axis=1)[have] / cnt[have]
    rows = np.flatnonzero(tip.any(axis=1) & have)
    for i in rows:
        out[i, tip[i]] = mean[i]
    return out


def joint_state(jd, u, cj_r, tanphi_r, slip_p=None, open_prev=None):
    """The joint tractions at displacement ``u``.

    ``slip_p`` is the accumulated plastic tangential offset (the viscoplastic
    driver's state); ``None`` means the stateless form the Newton path uses, in
    which the shear traction is returned straight onto the limit surface.

    Returns a dict with ``dt``, ``dn``, ``tn``, ``ts``, ``tlim``, ``open`` and
    ``slipping``, each (n, 3).
    """
    dt, dn = joint_kinematics(jd, u)
    kn, ks = jd["kn"][:, None], jd["ks"][:, None]
    tn = kn * dn
    ts_el = ks * (dt if slip_p is None else dt - slip_p)

    # Opening: a joint carrying tension past its cutoff parts, and stays parted
    # until the faces come back into contact. A tip pair, where the two soil
    # faces are one node, never opens.
    opened = tn < -jd["tcut"][:, None]
    if open_prev is not None:
        opened = opened | (open_prev & (dn < 0.0))
    opened &= ~jd["tip"]
    opened &= jd["w"] > 0.0

    tn_lim = _limit_normal(jd, tn)
    tlim = cj_r[:, None] + np.maximum(tn_lim, 0.0) * tanphi_r[:, None]
    slipping = (~opened) & (np.abs(ts_el) > tlim)

    tn_true = np.where(opened, 0.0, tn)
    ts_true = np.where(opened, 0.0,
                       np.where(slipping, np.sign(ts_el) * tlim, ts_el))
    return {"dt": dt, "dn": dn, "tn": tn_true, "ts": ts_true, "tlim": tlim,
            "open": opened, "slipping": slipping, "ts_trial": ts_el}


def joint_vp_sweep(jd, u, loads, cj_r, tanphi_r, slip_p, open_state,
                   dt_vp=JOINT_VP_DT):
    """One viscoplastic sweep over the joints: update the slip, load the residual.

    The global stiffness carries every joint's FULL elastic block, so ``K u``
    contains the unrestrained tractions ``k_s delta_t`` and ``k_n delta_n``. The
    part the interface cannot deliver is subtracted as a body load, exactly as
    the soil's plastic strain and the bar's tension cap are: the plastic
    tangential offset ``slip_p`` grows by ``dt (|t_s| - t_lim) sign / k_s`` on a
    slipping pair — at ``dt = 1`` that returns the shear traction to its limit at
    the current displacement field — and an open pair sheds the whole traction
    vector.

    Mutates ``slip_p`` and ``open_state`` in place, adds into ``loads``, and
    returns the number of pairs that are slipping or open.
    """
    st = joint_state(jd, u, cj_r, tanphi_r, slip_p=slip_p, open_prev=open_state)
    np.copyto(open_state, st["open"])

    ks = jd["ks"][:, None]
    excess = np.abs(st["ts_trial"]) - st["tlim"]
    grow = st["slipping"] & (excess > 0.0)
    if np.any(grow):
        slip_p += np.where(grow,
                           dt_vp * excess * np.sign(st["ts_trial"]) / ks, 0.0)

    # The correction, in the same sense as K u's own contribution: the elastic
    # traction minus the one the interface can carry. Open pairs shed both
    # components; closed ones shed only the plastic tangential offset.
    kn = jd["kn"][:, None]
    corr_t = np.where(st["open"], ks * st["dt"], ks * slip_p)
    corr_n = np.where(st["open"], kn * st["dn"], 0.0)
    # Back to global components. The internal force on side a of pair p is
    # w (t_s t - t_n n); the correction carries the same pattern.
    fx = jd["w"] * (corr_t * jd["tx"][:, None] - corr_n * jd["nx"][:, None])
    fy = jd["w"] * (corr_t * jd["ty"][:, None] - corr_n * jd["ny"][:, None])
    dof = jd["dof"]
    np.add.at(loads, dof[:, 0:6:2].ravel(), fx.ravel())
    np.add.at(loads, dof[:, 1:6:2].ravel(), fy.ravel())
    np.add.at(loads, dof[:, 6:12:2].ravel(), (-fx).ravel())
    np.add.at(loads, dof[:, 7:12:2].ravel(), (-fy).ravel())
    return int(np.count_nonzero(st["slipping"] | st["open"])), st


def tie_vp_sweep(td, u, loads):
    """The tied ends' body-load correction: the force past the tie's capacity.

    The tie is perfectly plastic at its capacity, so what it delivers is the
    elastic spring force scaled back onto the capacity circle. Returns the tie
    forces (m, 2) and the number of ties at capacity.
    """
    dof = td["dof"]
    ue = u[dof]
    d = ue[:, 0:2] - ue[:, 2:4]                     # bar minus soil
    f = td["k"][:, None] * d                        # force on the bar node
    mag = np.hypot(f[:, 0], f[:, 1])
    at_cap = mag > td["cap"]
    scale = np.ones(len(mag))
    with np.errstate(divide="ignore", invalid="ignore"):
        scale[at_cap] = td["cap"][at_cap] / mag[at_cap]
    f_true = f * scale[:, None]
    corr = f - f_true
    np.add.at(loads, dof[:, 0], corr[:, 0])
    np.add.at(loads, dof[:, 1], corr[:, 1])
    np.add.at(loads, dof[:, 2], -corr[:, 0])
    np.add.at(loads, dof[:, 3], -corr[:, 1])
    return f_true, int(np.count_nonzero(at_cap))


def joint_internal_force(jd, u, cj_r, tanphi_r, want_tangent=False):
    """The joints' internal force vector contribution, and optionally the tangent.

    The stateless form the Newton driver needs: the shear traction is returned
    onto the Mohr-Coulomb limit and the normal traction onto the tension cutoff,
    both as functions of the current displacement alone, so nothing is committed
    at the end of a step. The tangent drops ``k_s`` on a slipping pair and both
    stiffnesses on an open one.

    Returns ``(f, Ke, state)`` with ``f`` of shape (n, 12).
    """
    st = joint_state(jd, u, cj_r, tanphi_r, slip_p=None, open_prev=None)
    w = jd["w"]
    fx = w * (st["ts"] * jd["tx"][:, None] - st["tn"] * jd["nx"][:, None])
    fy = w * (st["ts"] * jd["ty"][:, None] - st["tn"] * jd["ny"][:, None])
    n = jd["n"]
    f = np.zeros((n, 12))
    f[:, 0:6:2] = fx
    f[:, 1:6:2] = fy
    f[:, 6:12:2] = -fx
    f[:, 7:12:2] = -fy
    Ke = None
    if want_tangent:
        ks_eff = np.where(st["slipping"] | st["open"], 0.0, jd["ks"][:, None])
        kn_eff = np.where(st["open"], 0.0, jd["kn"][:, None])
        Ke = np.zeros((n, 12, 12))
        tx, ty, nx, ny = jd["tx"], jd["ty"], jd["nx"], jd["ny"]
        for p in range(3):
            wp = w[:, p]
            kse, kne = ks_eff[:, p], kn_eff[:, p]
            D = np.empty((n, 2, 2))
            D[:, 0, 0] = wp * (kse * tx * tx + kne * nx * nx)
            D[:, 0, 1] = wp * (kse * tx * ty + kne * nx * ny)
            D[:, 1, 0] = D[:, 0, 1]
            D[:, 1, 1] = wp * (kse * ty * ty + kne * ny * ny)
            ia, ib = 2 * p, 6 + 2 * p
            Ke[:, ia:ia + 2, ia:ia + 2] += D
            Ke[:, ia:ia + 2, ib:ib + 2] -= D
            Ke[:, ib:ib + 2, ia:ia + 2] -= D
            Ke[:, ib:ib + 2, ib:ib + 2] += D
    return f, Ke, st


def tie_internal_force(td, u, want_tangent=False):
    """The ties' internal force (m, 4) and tangent, perfectly plastic at capacity."""
    dof = td["dof"]
    ue = u[dof]
    d = ue[:, 0:2] - ue[:, 2:4]
    f2 = td["k"][:, None] * d
    mag = np.hypot(f2[:, 0], f2[:, 1])
    at_cap = mag > td["cap"]
    scale = np.ones(len(mag))
    with np.errstate(divide="ignore", invalid="ignore"):
        scale[at_cap] = td["cap"][at_cap] / mag[at_cap]
    f2 = f2 * scale[:, None]
    f = np.concatenate([f2, -f2], axis=1)
    Ke = None
    if want_tangent:
        Ke = td["K"] * np.where(at_cap, 0.0, 1.0)[:, None, None]
    return f, Ke, at_cap


def joint_yield_violation(st, cj_r, tanphi_r, floor=0.0):
    """How far past its limit the worst shear traction sits, on its own strength.

    A pair AT its limit reads zero — a slipping joint is admissible — and one
    above it reads the fraction of its own strength it exceeds by.
    """
    den = np.maximum(st["tlim"], floor)
    ok = den > 0.0
    if not np.any(ok):
        return 0.0
    viol = (np.abs(st["ts"]) - st["tlim"])[ok] / den[ok]
    return float(max(0.0, np.max(viol)))


# ---------------------------------------------------------------------------
# What a BONDED run says about the joint it did not have
#
# The selection rule (docs/fem/reinforcement.md, "Bonded bar or joint?") is
# about the mechanism, and a bonded run has already found the mechanism. Two of
# its readings say the surface wanted to run ALONG a sheet rather than across
# it, and both are one pass over a solution the engine already holds.
# ---------------------------------------------------------------------------

#: How much of a sheet's length has to lie in the strain band before the band is
#: read as running ALONG the sheet. A surface that merely CUTS a sheet touches a
#: short stretch of it; one that follows it covers most of it.
JOINT_BAND_COVERAGE = 0.6

#: An element is IN the band when its shear strain is at least this fraction of
#: the largest in the model.
JOINT_BAND_STRAIN_FRAC = 0.5

#: A bar is at its cap within this fraction of it.
JOINT_CAP_TOL = 0.99


def _element_centroids(nodes, elements, element_types):
    """The centroid of every 2D element, from its corner nodes."""
    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=int)
    types = np.asarray(element_types, dtype=int)
    corners = np.where(types >= 6, types // 2, types)
    out = np.zeros((len(elements), 2))
    for k in np.unique(corners):
        rows = corners == k
        out[rows] = nodes[elements[rows][:, :int(k)], :2].mean(axis=1)
    return out


def _line_band_coverage(p1, p2, points, half_width):
    """How much of the segment p1-p2 has one of ``points`` within ``half_width``.

    Measured as the fraction of the segment's length covered by the union of the
    intervals each near point projects onto, so a cluster at one end reads as the
    short stretch it is.
    """
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    d = p2 - p1
    L = float(np.hypot(*d))
    if L <= 0 or len(points) == 0:
        return 0.0
    t = d / L
    rel = np.asarray(points, dtype=float)[:, :2] - p1
    s = rel @ t
    off = np.abs(rel[:, 0] * (-t[1]) + rel[:, 1] * t[0])
    near = (off <= half_width) & (s >= -half_width) & (s <= L + half_width)
    if not np.any(near):
        return 0.0
    lo = np.clip(s[near] - half_width, 0.0, L)
    hi = np.clip(s[near] + half_width, 0.0, L)
    order = np.argsort(lo)
    lo, hi = lo[order], hi[order]
    covered, end = 0.0, -np.inf
    for a, b in zip(lo, hi):
        if a > end:
            covered += b - a
            end = b
        elif b > end:
            covered += b - end
            end = b
    return float(covered / L)


def joint_advisories(fem_data, solution):
    """What a BONDED solution says about a sheet that wanted to be a joint.

    Two readings, each naming the line and suggesting ``Joint``:

    * the shear-strain band at the critical factor lies ALONG a sheet rather
      than across it, and
    * every bar element on one sheet sits at its capacity, so the sheet is
      holding the mass up through a grip the bond cap set rather than through
      an interface the mesh resolved.

    Reads ``fem_data['reinforcement_lines']`` — the per-line label, endpoints
    and ``Joint`` flag ``build_fem_data`` carries for exactly this — so a caller
    that has the solve has everything the reading needs.

    Returns a list of message strings, empty on a model with no reinforcement,
    on a line already flagged ``Joint``, and where neither reading fires.
    """
    from .mesh import line_is_jointed
    lines = (fem_data or {}).get("reinforcement_lines") or []
    if not lines or fem_data is None or not solution:
        return []
    elements_1d = fem_data.get("elements_1d")
    if elements_1d is None or len(elements_1d) == 0:
        return []
    mat_1d = np.asarray(fem_data.get("element_materials_1d", ()), dtype=int)
    pile_mask = np.asarray(fem_data.get("pile_elem_mask",
                                        np.zeros(len(elements_1d), dtype=bool)),
                           dtype=bool)
    forces = np.asarray(solution.get("forces_1d", ()), dtype=float)
    t_allow = np.asarray(fem_data.get("t_allow_by_1d_elem", ()), dtype=float)

    # The strain band, where the solution carries one.
    hot = np.zeros((0, 2))
    strains = solution.get("strains")
    half_width = 0.0
    if strains is not None:
        strains = np.asarray(strains, dtype=float)
        if strains.ndim == 2 and strains.shape[1] >= 4:
            shear = np.abs(strains[:, 3])
            peak = float(np.max(shear)) if len(shear) else 0.0
            if peak > 0:
                cen = _element_centroids(fem_data["nodes"], fem_data["elements"],
                                         fem_data["element_types"])
                hot = cen[shear >= JOINT_BAND_STRAIN_FRAC * peak]
                # The band's own width: one element. Read from the mesh rather
                # than chosen, so a fine mesh judges on a narrow band and a
                # coarse one on a wide one, which is what each can resolve.
                areas = float(np.ptp(cen[:, 0]) * np.ptp(cen[:, 1]))
                half_width = 1.5 * np.sqrt(max(areas, 0.0)
                                           / max(len(cen), 1)) if len(cen) else 0.0

    out = []
    for i, line in enumerate(lines):
        if line_is_jointed(line):
            continue
        label = line.get("label") or f"Reinforcement line {i + 1}"
        try:
            p1 = (float(line["x1"]), float(line["y1"]))
            p2 = (float(line["x2"]), float(line["y2"]))
        except (KeyError, TypeError, ValueError):
            continue

        if len(hot) and half_width > 0:
            cover = _line_band_coverage(p1, p2, hot, half_width)
            if cover >= JOINT_BAND_COVERAGE:
                out.append(
                    f"The shear strain band at the critical factor runs along "
                    f"{cover * 100:.0f}% of '{label}' rather than across it. A "
                    f"surface that follows a sheet is a surface ON the sheet, "
                    f"and a bonded bar cannot carry one: the soil above and "
                    f"below the line share its nodes. Set Joint = Yes on that "
                    f"line to split the mesh along it and give the interface "
                    f"its own strength (reinforce sheet, Adhesion/Delta).")

        rows = np.flatnonzero((mat_1d == i + 1) & ~pile_mask) if len(mat_1d) else []
        if (len(rows) >= 3 and len(forces) == len(elements_1d)
                and len(t_allow) == len(elements_1d)
                and np.all(t_allow[rows] > 0)
                and np.all(forces[rows] >= JOINT_CAP_TOL * t_allow[rows])):
            out.append(
                f"Every one of the {len(rows)} bar elements on '{label}' sits "
                f"at its capacity at the critical factor. A sheet whose whole "
                f"length is at its cap is not carrying a computed grip on the "
                f"soil: the bond cap is what is holding the mass, and refining "
                f"the bar elements refines it instead of converging. Set "
                f"Joint = Yes on that line to make the grip the traction the "
                f"interface elements integrate (reinforce sheet, "
                f"Adhesion/Delta).")
    return out
