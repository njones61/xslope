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

"""Per-line detail profiles for the FEM's 1D elements.

The FEM stores reinforcement and pile results one row per 1D element, indexed by
mesh element id. That is the right storage order for the solver and for the
force-coloured overlays on the result plots, but it is not the order an engineer
reads a structural member in: along the member, from one end to the other.

This module turns the per-element arrays into **per-line profiles** — position
along the member plus every series measured at that position — and is the single
source those profiles come from. The Studio detail dialog plots them, the
dialog's CSV export writes them, and the checks assert against them, so all
three necessarily agree.

Two properties matter and are deliberate:

*Capacity is never re-derived here.* The reinforcement capacity envelope is
evaluated by calling :func:`xslope.fileio.reinforce_available_tension` — the
same function ``build_fem_data`` calls to fill ``t_allow_by_1d_elem``. A plotted
envelope and the solver's element capacities are therefore the same curve
sampled at different resolutions, and cannot drift.

*The field is the caller's choice.* A profile is read from the converged
solution or from the at-failure snapshot SSRM captured, selected by
``field_state`` — the same switch, the same two names, and the same automatic
fallback as ``plot_fem_results``, so a detail view and a results panel set to
the same state show the same instant of the analysis.

*Everything works from a reloaded solution.* Every quantity below comes either
from ``fem_data`` (rebuilt from the model and mesh on load, so identical) or
from solution fields that ``import_fem_solution`` restores from the sidecars.
Nothing reads a value that exists only in the middle of a live solve.
"""

import numpy as np

from .fileio import reinforce_available_tension

# Utilization at or above which a member is reported as at capacity, and the
# threshold between the "comfortable" and "watch" bands. These set the badge
# colours in the Studio list and the status word in the profile dicts.
UTIL_AT_CAPACITY = 0.995
UTIL_WATCH = 0.70

# A 1D element sits inside the failure band when the mechanism field at its
# location reaches this fraction of the field's maximum along that member. Read
# from the at-failure snapshot when one was captured, else from the converged
# field.
BAND_FRACTION = 0.5

# The two fields a profile can be read from, spelled exactly as
# ``plot_fem_results``' ``field_state`` switch spells them, so the detail views
# and the result plots name the same states.
FIELD_STATES = ("converged", "failure")


# --------------------------------------------------------------------------
# which field a profile is read from
# --------------------------------------------------------------------------

def failure_snapshot(solution, failure_solution=None):
    """The at-failure (unconverged) field SSRM captured, or None.

    There is one source for it and this is the only place it is looked up.
    :func:`xslope.fem.solve_ssrm` returns it as ``result['failure_solution']``
    and :func:`xslope.fem.import_fem_solution` restores it — nodal field,
    element field and the 1D result sidecars alike — nested under
    ``solution['failure_solution']``. A caller holding it separately (Studio
    lifts it onto its results bundle, which is the shape ``plot_fem_results``
    takes as its ``failure_solution`` argument) passes it in; an explicit
    argument wins over the nested one.
    """
    snap = failure_solution
    if snap is None:
        snap = (solution or {}).get("failure_solution")
    return snap or None


def has_failure_state(solution, failure_solution=None):
    """True when there is an at-failure snapshot to read profiles from."""
    return failure_snapshot(solution, failure_solution) is not None


def field_solution(solution, field_state="converged", failure_solution=None):
    """The solution dict a profile's member forces and displacements come from.

    ``'failure'`` reads the captured at-failure snapshot, ``'converged'`` the
    converged solution. With no snapshot the two selections are the same field —
    there is only one to read — which is the automatic fallback
    ``plot_fem_results`` makes for the same switch.
    """
    if field_state not in FIELD_STATES:
        raise ValueError(f"Unknown field_state: '{field_state}'. "
                         f"Valid: {', '.join(FIELD_STATES)}.")
    if field_state == "failure":
        snap = failure_snapshot(solution, failure_solution)
        if snap is not None:
            return snap
    return solution


def effective_field_state(solution, field_state="converged", failure_solution=None):
    """The state a profile was actually read at: ``'failure'`` only when a
    snapshot existed to read, otherwise ``'converged'``."""
    if field_state == "failure" and has_failure_state(solution, failure_solution):
        return "failure"
    return "converged"


def field_state_label(field_state):
    """The display name for a field state — the wording the results view's own
    Field state control uses."""
    return "at failure" if field_state == "failure" else "last converged"


# --------------------------------------------------------------------------
# presence / discovery
# --------------------------------------------------------------------------

def has_1d_details(fem_data):
    """True when the model carries any 1D element — reinforcement bar or pile
    beam — that the detail profiles can describe."""
    if not fem_data:
        return False
    elements_1d = fem_data.get("elements_1d", None)
    return elements_1d is not None and len(elements_1d) > 0


def _n_reinforcement_lines(fem_data, slope_data=None):
    labels = fem_data.get("reinforce_line_labels", None)
    if labels is not None:
        return len(labels)
    if slope_data:
        return len(slope_data.get("reinforcement_lines", []) or [])
    return 0


def _reinforcement_line_ids(fem_data):
    """1-based line ids that actually own reinforcement (non-pile) elements."""
    if not has_1d_details(fem_data):
        return []
    n_1d = len(fem_data["elements_1d"])
    mats = np.asarray(fem_data.get("element_materials_1d", np.zeros(n_1d)), dtype=int)
    mask = np.asarray(fem_data.get("pile_elem_mask", np.zeros(n_1d, dtype=bool)), dtype=bool)
    if len(mats) != n_1d:
        return []
    if len(mask) != n_1d:
        mask = np.zeros(n_1d, dtype=bool)
    return sorted(int(v) for v in np.unique(mats[~mask]) if v > 0)


def _pile_line_indices(fem_data):
    """0-based indices into ``slope_data['pile_lines']`` that own pile elements."""
    n_pile = int(fem_data.get("n_pile_elements", 0) or 0)
    if n_pile == 0:
        return []
    by_elem = np.asarray(fem_data.get("pile_line_idx_by_pile_elem", []), dtype=int)
    if len(by_elem) != n_pile:
        return []
    return sorted(int(v) for v in np.unique(by_elem))


def _line_label(fem_data, slope_data, kind, index):
    """A display name for one member.

    The model's own label is used when it has one. Labels need not be unique —
    a piles sheet where every row is called "pile" is perfectly legal — so a
    label shared with another member of the same kind is numbered, otherwise the
    list would show two rows a user cannot tell apart.
    """
    if kind == "reinforcement":
        labels = [str(x) if x else "" for x in
                  (fem_data.get("reinforce_line_labels", None) or [])]
        own = labels[index - 1] if 0 <= index - 1 < len(labels) else ""
        if not own:
            return f"Line {index}"
        return own if labels.count(own) == 1 else f"{own} {index}"
    lines = (slope_data or {}).get("pile_lines", []) or []
    labels = [str(ln.get("label")) if ln.get("label") else "" for ln in lines]
    own = labels[index] if 0 <= index < len(labels) else ""
    if not own:
        return f"Pile {index + 1}"
    return own if labels.count(own) == 1 else f"{own} {index + 1}"


def _badge(util):
    """Badge colour name for a utilization ratio (None when unmeasurable)."""
    if util is None or not np.isfinite(util):
        return "none"
    if util >= UTIL_AT_CAPACITY:
        return "red"
    if util >= UTIL_WATCH:
        return "amber"
    return "green"


def list_lines(fem_data, solution, slope_data=None, field_state="converged",
               failure_solution=None):
    """Every reinforcement line and pile in the model, in list order.

    Returns a list of dicts with ``kind`` (``"reinforcement"`` / ``"pile"``),
    ``index`` (the 1-based reinforcement line id, or the 0-based pile line
    index), ``label``, ``n_elements``, ``utilization`` (peak, or None when the
    model supplies no capacity to measure against), ``badge`` and ``status``.

    ``field_state`` selects the field the utilizations are measured on, exactly
    as it does for the profiles themselves, so the badges and the plotted
    profile always report the same state.

    Grouped output is the caller's job; the two kinds are returned in order,
    reinforcement first, matching the order the solver assigns line ids.
    """
    out = []
    if not has_1d_details(fem_data):
        return out
    for line_id in _reinforcement_line_ids(fem_data):
        prof = reinforcement_profile(fem_data, solution, line_id, slope_data,
                                     field_state=field_state,
                                     failure_solution=failure_solution)
        out.append({
            "kind": "reinforcement",
            "index": line_id,
            "label": prof["label"],
            "n_elements": len(prof["s"]),
            "utilization": prof["peak_utilization"],
            "badge": prof["badge"],
            "status": prof["status"],
        })
    for pidx in _pile_line_indices(fem_data):
        prof = pile_profile(fem_data, solution, pidx, slope_data,
                            field_state=field_state,
                            failure_solution=failure_solution)
        out.append({
            "kind": "pile",
            "index": pidx,
            "label": prof["label"],
            "n_elements": prof["n_elements"],
            "utilization": prof["peak_utilization"],
            "badge": prof["badge"],
            "status": prof["status"],
        })
    return out


# --------------------------------------------------------------------------
# mechanism field sampling (failure-band crossing)
# --------------------------------------------------------------------------

def _element_centroids(fem_data):
    nodes = np.asarray(fem_data["nodes"], dtype=float)
    elements = np.asarray(fem_data["elements"], dtype=int)
    types = np.asarray(fem_data.get("element_types", np.full(len(elements), 3)), dtype=int)
    cents = np.zeros((len(elements), 2))
    for i in range(len(elements)):
        ncor = 3 if int(types[i]) in (3, 6) else 4
        cents[i] = nodes[elements[i][:ncor]].mean(axis=0)
    return cents


def _mechanism_field(fem_data, solution, failure_solution=None):
    """The captured mechanism's viscoplastic shear strain per 2D element, with
    the at-failure snapshot preferred over the converged field when one was
    captured. The failure band is a property of the mechanism, so it is read
    from the snapshot in both field states — the band marks stay put when the
    profiles switch between the converged and the at-failure field. Returns None
    when the field is unavailable."""
    if not solution:
        return None
    src = failure_snapshot(solution, failure_solution) or solution
    field = src.get("vp_shear_strain", None)
    if field is None:
        return None
    field = np.asarray(field, dtype=float)
    n_elem = len(fem_data.get("elements", []))
    if field.shape[0] != n_elem or n_elem == 0:
        return None
    return field


def _sample_mechanism(fem_data, solution, points, failure_solution=None):
    """Mechanism field sampled at ``points`` (n,2) by nearest 2D element, or
    None when there is no mechanism field to sample."""
    field = _mechanism_field(fem_data, solution, failure_solution)
    if field is None or len(points) == 0:
        return None
    cents = _element_centroids(fem_data)
    if len(cents) == 0:
        return None
    try:
        from scipy.spatial import cKDTree
        _, idx = cKDTree(cents).query(np.asarray(points, dtype=float))
    except Exception:
        pts = np.asarray(points, dtype=float)
        idx = np.array([int(np.argmin(((cents - p) ** 2).sum(axis=1))) for p in pts])
    return field[np.asarray(idx, dtype=int)]


def _band_span(positions, mech):
    """Contiguous run of ``positions`` around the mechanism peak whose field
    reaches ``BAND_FRACTION`` of that peak — the failure band as this member
    sees it. Returns ``(lo, hi, peak_position)`` or ``(None, None, None)``."""
    if mech is None or len(mech) == 0:
        return None, None, None
    mech = np.asarray(mech, dtype=float)
    peak = float(np.nanmax(mech))
    if not np.isfinite(peak) or peak <= 0:
        return None, None, None
    k = int(np.nanargmax(mech))
    thresh = BAND_FRACTION * peak
    lo = k
    while lo - 1 >= 0 and mech[lo - 1] >= thresh:
        lo -= 1
    hi = k
    while hi + 1 < len(mech) and mech[hi + 1] >= thresh:
        hi += 1
    return float(positions[lo]), float(positions[hi]), float(positions[k])


# --------------------------------------------------------------------------
# reinforcement
# --------------------------------------------------------------------------

def _reinforcement_line_inputs(fem_data, slope_data, line_id):
    """The capacity-envelope inputs for one reinforcement line, or None.

    Read from ``slope_data['reinforcement_lines']`` — already reduced to per
    unit width of slope by ``load_slope_data``, the same values
    ``build_fem_data`` passed to the envelope function — so the analytic
    envelope and ``t_allow_by_1d_elem`` are two samplings of one curve.
    """
    lines = (slope_data or {}).get("reinforcement_lines", []) or []
    i = line_id - 1
    if not (0 <= i < len(lines)):
        return None
    ln = lines[i]
    x1, y1 = float(ln.get("x1", 0.0)), float(ln.get("y1", 0.0))
    x2, y2 = float(ln.get("x2", 0.0)), float(ln.get("y2", 0.0))
    length = float(np.hypot(x2 - x1, y2 - y1))
    if length <= 0:
        return None
    return {
        "x1": x1, "y1": y1, "x2": x2, "y2": y2, "length": length,
        "t_max": float(ln.get("t_max", 0.0) or 0.0),
        "t_res": float(ln.get("t_res", float("nan"))),
        "lp1": float(ln.get("lp1", 0.0) or 0.0),
        "lp2": float(ln.get("lp2", 0.0) or 0.0),
        "tend1": float(ln.get("tend1", 0.0) or 0.0),
        "tend2": float(ln.get("tend2", 0.0) or 0.0),
    }


def capacity_envelope(inputs, n=241):
    """Sample the reinforcement capacity envelope over a whole line.

    ``inputs`` is the dict from :func:`_reinforcement_line_inputs`. Returns
    ``(s, T)`` where ``s`` runs from end 1 to end 2 and ``T`` is
    :func:`xslope.fileio.reinforce_available_tension` — the pullout ramp from
    each free end over its development length, the flat tensile plateau at
    ``Tmax`` in the middle, and the step to the connection capacity ``Tend`` at
    an anchored end. The sample grid always includes the envelope's own kinks,
    so the polyline drawn through it is exact rather than merely dense.
    """
    L = inputs["length"]
    t_max, lp1, lp2 = inputs["t_max"], inputs["lp1"], inputs["lp2"]
    tend1, tend2 = inputs["tend1"], inputs["tend2"]

    kinks = {0.0, L}
    if lp1 > 0 and tend1 < t_max and t_max > 0:
        kinks.add(min(L, (t_max - tend1) * lp1 / t_max))
    if lp2 > 0 and tend2 < t_max and t_max > 0:
        kinks.add(max(0.0, L - (t_max - tend2) * lp2 / t_max))
    if lp1 > 0 and lp2 > 0 and t_max > 0:
        m1, m2 = t_max / lp1, t_max / lp2
        s_x = (tend2 + m2 * L - tend1) / (m1 + m2)
        if 0.0 < s_x < L and (tend1 + m1 * s_x) < t_max:
            kinks.add(s_x)
    s = np.unique(np.concatenate([np.linspace(0.0, L, n), np.array(sorted(kinks))]))
    T = np.array([reinforce_available_tension(float(v), L - float(v), t_max,
                                              lp1, lp2, tend1, tend2) for v in s])
    return s, T


def reinforcement_profile(fem_data, solution, line_id, slope_data=None,
                          field_state="converged", failure_solution=None):
    """Everything the detail view draws for one reinforcement line.

    ``field_state`` selects the field the mobilized force and the element state
    flags are read from — ``'converged'`` (default) the converged solution,
    ``'failure'`` the at-failure snapshot when one was captured (see
    :func:`field_solution`). The capacity envelope is a property of the model
    and does not move with it, and neither do the failure-band marks, which
    always follow the captured mechanism.

    Keys
    ----
    ``kind``, ``index``, ``label``, ``length``
    ``s`` : arc length from end 1 to each element centroid, ascending
    ``x``, ``y`` : centroid coordinates
    ``T`` : mobilized axial force at each centroid
    ``t_allow`` : the solver's capacity at each centroid (``t_allow_by_1d_elem``)
    ``t_cap`` : the capacity the solve actually enforced — identical to
        ``t_allow`` unless the optional bond-slip model re-capped the line
    ``t_res`` : residual (post-peak) capacity, NaN where the line never softens
    ``failed``, ``softened`` : per-element state flags
    ``field_state`` : the state the series above were read at
    ``env_s``, ``env_T`` : the analytic capacity envelope (None without model
        inputs, in which case the per-element ``t_cap`` is the capacity to draw)
    ``bond_s``, ``bond_q`` : mobilized bond force per unit length, dT/ds, at the
        boundaries between consecutive elements
    ``slip_modelled`` : False — see note below
    ``band_lo``, ``band_hi``, ``band_peak`` : failure-band extents in ``s``
    ``peak_s``, ``peak_T``, ``peak_utilization``, ``badge``, ``status``
    ``pullout_s``, ``softened_s`` : the ``s`` of elements in each state
    ``units`` : the model's display unit strings

    Relative slip has no entry because the model has no slip degree of freedom:
    reinforcement bars are truss elements on the continuum's own nodes, so bar
    and soil displacement are the same number at every node and their difference
    is identically zero. Load transfer is expressed through the capacity
    envelope, not through a slip law, and the bond series below is the mobilized
    transfer rate implied by the axial force gradient.
    """
    state = effective_field_state(solution, field_state, failure_solution)
    field = field_solution(solution, field_state, failure_solution)

    n_1d = len(fem_data.get("elements_1d", []))
    mats = np.asarray(fem_data.get("element_materials_1d", np.zeros(n_1d)), dtype=int)
    mask = np.asarray(fem_data.get("pile_elem_mask", np.zeros(n_1d, dtype=bool)), dtype=bool)
    if len(mask) != n_1d:
        mask = np.zeros(n_1d, dtype=bool)
    idx = np.where((mats == line_id) & (~mask))[0] if len(mats) == n_1d else np.array([], int)

    d1 = _elem_array(fem_data, "dist_end1_1d", n_1d)
    d2 = _elem_array(fem_data, "dist_end2_1d", n_1d)
    elen = _elem_array(fem_data, "elem_length_1d", n_1d)
    t_allow = _elem_array(fem_data, "t_allow_by_1d_elem", n_1d)
    t_res = _elem_array(fem_data, "t_res_by_1d_elem", n_1d, fill=np.nan)
    forces = _sol_array(field, "forces_1d", n_1d)
    failed = _sol_array(field, "failed_1d_elements", n_1d, dtype=bool)
    softened = _sol_array(field, "softened_1d_elements", n_1d, dtype=bool)
    t_cap_all = _sol_array(field, "t_cap_1d", n_1d, fill=np.nan)
    if not np.any(np.isfinite(t_cap_all)):
        t_cap_all = t_allow

    idx = idx[elen[idx] > 0] if len(idx) else idx
    order = np.argsort(d1[idx], kind="stable") if len(idx) else np.array([], int)
    idx = idx[order]

    nodes = np.asarray(fem_data.get("nodes", np.zeros((0, 2))), dtype=float)
    elements_1d = np.asarray(fem_data.get("elements_1d", np.zeros((0, 3))), dtype=int)
    if len(idx):
        pts = np.array([0.5 * (nodes[elements_1d[i][0]] + nodes[elements_1d[i][1]])
                        for i in idx])
    else:
        pts = np.zeros((0, 2))

    s = d1[idx]
    inputs = _reinforcement_line_inputs(fem_data, slope_data, line_id)
    if inputs is not None:
        length = inputs["length"]
        env_s, env_T = capacity_envelope(inputs)
    else:
        length = float(np.max(d1[idx] + d2[idx])) if len(idx) else 0.0
        env_s = env_T = None

    T = forces[idx]
    cap = t_cap_all[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        util = np.where(cap > 1e-12, T / cap, np.nan)

    # Mobilized bond transfer rate dT/ds, at the boundary between consecutive
    # elements — the force the ground hands the bar per unit of its length.
    if len(idx) >= 2:
        bond_s = 0.5 * (s[:-1] + s[1:])
        ds = np.diff(s)
        bond_q = np.where(ds > 1e-12, np.diff(T) / np.where(ds > 1e-12, ds, 1.0), 0.0)
    else:
        bond_s = np.zeros(0)
        bond_q = np.zeros(0)

    mech = _sample_mechanism(fem_data, solution, pts, failure_solution)
    band_lo, band_hi, band_peak = _band_span(s, mech) if len(idx) else (None, None, None)

    peak_i = int(np.nanargmax(util)) if len(idx) and np.any(np.isfinite(util)) else None
    peak_util = float(util[peak_i]) if peak_i is not None else None

    # Pulled out, on the renderer's own definition: the element yielded, its
    # residual capacity is (finitely) zero, and it now carries no force. A NaN
    # residual means the line never softens at all, which is not pullout.
    tr = t_res[idx]
    pulled = (failed[idx] & np.isfinite(tr) & (tr < 1e-6) & (T < 1e-6)) if len(idx) \
        else np.zeros(0, dtype=bool)
    soft = softened[idx] if len(idx) else np.zeros(0, dtype=bool)
    pulled_out = bool(np.any(pulled))

    status = _reinforcement_status(peak_util, pulled_out, bool(np.any(soft)))
    badge = "red" if pulled_out else _badge(peak_util)

    return {
        "kind": "reinforcement",
        "index": int(line_id),
        "label": _line_label(fem_data, slope_data, "reinforcement", line_id),
        "field_state": state,
        "length": length,
        "element_ids": idx,
        "s": s, "x": pts[:, 0] if len(pts) else np.zeros(0),
        "y": pts[:, 1] if len(pts) else np.zeros(0),
        "T": T,
        "t_allow": t_allow[idx],
        "t_cap": cap,
        "t_res": t_res[idx],
        "utilization": util,
        "failed": failed[idx],
        "softened": softened[idx],
        "env_s": env_s, "env_T": env_T,
        "bond_s": bond_s, "bond_q": bond_q,
        "slip_modelled": False,
        "mechanism": mech,
        "band_lo": band_lo, "band_hi": band_hi, "band_peak": band_peak,
        "peak_s": float(s[peak_i]) if peak_i is not None else None,
        "peak_T": float(T[peak_i]) if peak_i is not None else None,
        "peak_utilization": peak_util,
        "pullout_s": s[pulled] if len(idx) else np.zeros(0),
        "softened_s": s[soft] if len(idx) else np.zeros(0),
        "badge": badge,
        "status": status,
        "units": unit_labels(fem_data),
    }


def _reinforcement_status(util, pulled_out, softened):
    if pulled_out:
        return "pulled out"
    if softened:
        return "softened to residual"
    if util is None or not np.isfinite(util):
        return "no capacity declared"
    if util >= UTIL_AT_CAPACITY:
        return "at capacity"
    if util >= UTIL_WATCH:
        return "near capacity"
    return "within capacity"


# --------------------------------------------------------------------------
# piles
# --------------------------------------------------------------------------

def _pile_line_props(slope_data, index):
    lines = (slope_data or {}).get("pile_lines", []) or []
    if 0 <= index < len(lines):
        return lines[index]
    return {}


def _ito_matsui_limit(slope_data, props, depths, y_head, S):
    """Ito & Matsui limiting lateral soil resistance per unit width of slope,
    sampled at ``depths`` below the pile head.

    ``p(z) = (c*A1 + gamma*z*A2) / S`` with the ``A1``/``A2`` coefficients from
    :func:`xslope.ito_matsui.ito_matsui_coefficients` — the same theory and the
    same code the LEM uses for its passive-pile force, evaluated here as a
    per-depth envelope instead of integrated to a resultant. Returns None when
    the inputs the theory needs are absent (no pile diameter, no spacing, or a
    spacing at or below the diameter, where soil arching between piles does not
    apply).
    """
    D = props.get("D_pile", None)
    if D is None or not np.isfinite(D) or D <= 0:
        return None
    if S is None or not np.isfinite(S) or S <= D:
        return None
    polygons = (slope_data or {}).get("polygons", None)
    materials = (slope_data or {}).get("materials", None)
    if not polygons or not materials or len(depths) == 0:
        return None
    try:
        from .ito_matsui import intersect_pile_with_materials, ito_matsui_coefficients
        x_pile = float(props.get("x1", 0.0))
        y_toe = y_head - float(np.max(depths))
        segments = intersect_pile_with_materials(x_pile, y_head, y_toe, polygons, materials)
    except Exception:
        return None
    if not segments:
        return None
    D1 = S - D
    p = np.full(len(depths), np.nan)
    for seg in segments:
        A1, A2 = ito_matsui_coefficients(D1, D, seg["phi"])
        in_seg = (depths >= seg["z_top"] - 1e-9) & (depths <= seg["z_bot"] + 1e-9)
        p[in_seg] = (seg["c"] * A1 + seg["gamma"] * depths[in_seg] * A2) / S
    return p


def pile_profile(fem_data, solution, pile_index, slope_data=None,
                 field_state="converged", failure_solution=None):
    """Everything the detail view draws for one pile, ordered head to toe.

    ``field_state`` selects the field the displacements, shears and moments are
    read from — ``'converged'`` (default) the converged solution, ``'failure'``
    the at-failure snapshot when one was captured (see :func:`field_solution`).
    The Ito & Matsui limiting-resistance envelope is a property of the model and
    does not move with it, and neither does the failure-band depth, which always
    follows the captured mechanism.

    Keys
    ----
    ``kind``, ``index``, ``label``, ``length``, ``n_elements``
    ``node_depth`` : depth below the pile head at each pile node
    ``u_lateral`` : displacement normal to the pile axis at each node
    ``elem_depth`` : depth at each beam element's midpoint
    ``shear`` : element shear force, ``V_cap`` : its capacity where declared
    ``moment_depth``, ``moment`` : the continuous bending-moment profile,
        assembled from the beam elements' end moments; ``M_cap`` where declared
    ``reaction_depth``, ``reaction`` : lateral soil reaction per unit length,
        from the shear discontinuity the soil imposes at each interior node
    ``limit_depth``, ``limit_p`` : the Ito & Matsui limiting resistance
        envelope (None when the model does not supply pile diameter and spacing)
    ``max_moment``, ``max_moment_depth``
    ``band_depth`` : depth at which the failure band crosses the pile
    ``peak_utilization``, ``badge``, ``status``, ``utilization_basis``
    ``field_state`` : the state the series above were read at
    ``units``

    Bending moment is assembled so that it is continuous: the beam formulation
    reports each element's two end moments in the element's own sense, in which
    the second end's value is the negative of the first end's value of the next
    element. No moment capacity is drawn unless the model declares ``Mcap`` —
    the inputs carry force capacities, not section properties, so a capacity
    line without one would be invented.
    """
    state = effective_field_state(solution, field_state, failure_solution)
    field = field_solution(solution, field_state, failure_solution)

    n_pile = int(fem_data.get("n_pile_elements", 0) or 0)
    by_line = np.asarray(fem_data.get("pile_line_idx_by_pile_elem", []), dtype=int)
    sel = np.where(by_line == pile_index)[0] if len(by_line) == n_pile else np.array([], int)

    props = _pile_line_props(slope_data, pile_index)
    units = unit_labels(fem_data)
    empty = {
        "kind": "pile", "index": int(pile_index),
        "label": _line_label(fem_data, slope_data, "pile", pile_index),
        "field_state": state,
        "length": 0.0, "n_elements": 0,
        "node_depth": np.zeros(0), "u_lateral": np.zeros(0),
        "elem_depth": np.zeros(0), "shear": np.zeros(0),
        "moment_depth": np.zeros(0), "moment": np.zeros(0),
        "reaction_depth": np.zeros(0), "reaction": np.zeros(0),
        "limit_depth": None, "limit_p": None, "reaction_ratio": None,
        "V_cap": None, "M_cap": None,
        "max_moment": None, "max_moment_depth": None,
        "band_depth": None, "peak_utilization": None,
        "badge": "none", "status": "no results", "utilization_basis": None,
        "units": units,
    }
    if len(sel) == 0:
        return empty

    nodes = np.asarray(fem_data["nodes"], dtype=float)
    pairs = list(fem_data.get("pile_node_pairs", []))
    cos_t = np.asarray(fem_data.get("cos_theta_pile", np.zeros(n_pile)), dtype=float)
    sin_t = np.asarray(fem_data.get("sin_theta_pile", np.zeros(n_pile)), dtype=float)

    # Order elements head (highest node) to toe.
    def _elem_top(p):
        n0, n1 = pairs[p]
        return max(nodes[n0][1], nodes[n1][1])
    sel = sel[np.argsort([-_elem_top(p) for p in sel], kind="stable")]

    # Chain the node list head to toe.
    node_ids = []
    for p in sel:
        n0, n1 = pairs[p]
        if nodes[n0][1] < nodes[n1][1]:
            n0, n1 = n1, n0
        if not node_ids:
            node_ids.append(int(n0))
        node_ids.append(int(n1))
    node_ids = np.array(node_ids, dtype=int)
    xy = nodes[node_ids]
    y_head = float(xy[0][1])
    node_depth = y_head - xy[:, 1]
    length = float(node_depth[-1])

    # Lateral displacement: the component of nodal displacement normal to the
    # pile axis. The axis has two normals; take the one pointing in +x, so a
    # vertical pile reports plain downslope movement whichever way its two
    # endpoints were entered.
    disp = np.asarray((field or {}).get("displacements", np.zeros(0)), dtype=float)
    dof_offset = fem_data.get("dof_offset", None)
    c, s_ = float(cos_t[sel[0]]), float(sin_t[sel[0]])
    nx, ny = -s_, c
    if nx < 0:
        nx, ny = -nx, -ny
    u_lat = np.zeros(len(node_ids))
    for k, nid in enumerate(node_ids):
        base = int(dof_offset[nid]) if dof_offset is not None else 2 * nid
        if base + 1 < len(disp):
            u_lat[k] = disp[base] * nx + disp[base + 1] * ny

    shear = _sol_array(field, "forces_pile_lateral", n_pile)[sel]
    fm = np.asarray((field or {}).get("forces_pile_moment", np.zeros((n_pile, 2))),
                    dtype=float)
    if fm.shape != (n_pile, 2):
        fm = np.zeros((n_pile, 2))
    fm = fm[sel]
    elem_len = _pile_array(fem_data, "elem_length_by_pile_elem", n_pile)[sel]
    elem_depth = 0.5 * (node_depth[:-1] + node_depth[1:])

    # Continuous bending-moment profile: M1 of each element at its upper node,
    # and -M2 at its lower node (which equals the next element's M1).
    moment_depth = np.empty(2 * len(sel))
    moment = np.empty(2 * len(sel))
    for k in range(len(sel)):
        moment_depth[2 * k] = node_depth[k]
        moment_depth[2 * k + 1] = node_depth[k + 1]
        moment[2 * k] = fm[k, 0]
        moment[2 * k + 1] = -fm[k, 1]

    # Lateral soil reaction per unit length: the shear discontinuity the soil
    # imposes at each interior node, spread over that node's tributary length.
    if len(sel) >= 2:
        trib = 0.5 * (elem_len[:-1] + elem_len[1:])
        reaction_depth = node_depth[1:-1]
        reaction = np.where(trib > 1e-12, (shear[:-1] - shear[1:]) / np.where(trib > 1e-12, trib, 1.0), 0.0)
    else:
        reaction_depth = np.zeros(0)
        reaction = np.zeros(0)

    S = float(props.get("S")) if props.get("S") else (
        float(_pile_array(fem_data, "S_by_pile_elem", n_pile)[sel[0]]) or None)
    limit_depth = np.linspace(0.0, length, 121) if length > 0 else np.zeros(0)
    limit_p = _ito_matsui_limit(slope_data, props, limit_depth, y_head, S)
    if limit_p is None:
        limit_depth = None

    v_caps = _pile_array(fem_data, "V_cap_by_pile_elem", n_pile, fill=np.inf)[sel]
    m_caps = _pile_array(fem_data, "M_cap_by_pile_elem", n_pile, fill=np.inf)[sel]
    V_cap = float(np.min(v_caps)) if np.all(np.isfinite(v_caps)) and len(v_caps) else None
    M_cap = float(np.min(m_caps)) if np.all(np.isfinite(m_caps)) and len(m_caps) else None

    k_max = int(np.argmax(np.abs(moment))) if len(moment) else None
    max_moment = float(moment[k_max]) if k_max is not None else None
    max_moment_depth = float(moment_depth[k_max]) if k_max is not None else None

    mech_pts = np.column_stack([np.full(len(elem_depth), xy[0][0]),
                                y_head - elem_depth]) if len(elem_depth) else np.zeros((0, 2))
    if len(elem_depth):
        mech_pts[:, 0] = 0.5 * (xy[:-1, 0] + xy[1:, 0])
    mech = _sample_mechanism(fem_data, solution, mech_pts, failure_solution)
    _, _, band_depth = _band_span(elem_depth, mech) if len(elem_depth) else (None, None, None)

    reaction_ratio = _reaction_ratio(reaction, reaction_depth, limit_depth, limit_p)
    util, basis = _pile_utilization(shear, moment, V_cap, M_cap, reaction_ratio)

    return {
        "kind": "pile",
        "index": int(pile_index),
        "label": _line_label(fem_data, slope_data, "pile", pile_index),
        "field_state": state,
        "length": length,
        "n_elements": int(len(sel)),
        "pile_indices": sel,
        "node_ids": node_ids,
        "x": xy[:, 0],
        "node_depth": node_depth,
        "u_lateral": u_lat,
        "elem_depth": elem_depth,
        "shear": shear,
        "moment_depth": moment_depth,
        "moment": moment,
        "reaction_depth": reaction_depth,
        "reaction": reaction,
        "limit_depth": limit_depth,
        "limit_p": limit_p,
        "reaction_ratio": reaction_ratio,
        "V_cap": V_cap,
        "M_cap": M_cap,
        "max_moment": max_moment,
        "max_moment_depth": max_moment_depth,
        "mechanism": mech,
        "band_depth": band_depth,
        "peak_utilization": util,
        "utilization_basis": basis,
        "badge": _badge(util),
        "status": _pile_status(util, basis),
        "units": units,
    }


def _reaction_ratio(reaction, reaction_depth, limit_depth, limit_p):
    """Peak mobilized soil reaction as a fraction of the Ito & Matsui limiting
    resistance at the same depth, or None when the limit is unavailable."""
    if limit_p is None or len(reaction) == 0:
        return None
    p_at = np.interp(reaction_depth, limit_depth, limit_p)
    ok = np.isfinite(p_at) & (p_at > 1e-12)
    if not np.any(ok):
        return None
    return float(np.max(np.abs(reaction[ok]) / p_at[ok]))


def _pile_utilization(shear, moment, V_cap, M_cap, reaction_ratio):
    """Peak utilization for a pile and the quantity it was measured against.

    Structural capacities win when the model declares them, because they are the
    stated limits. Otherwise the mobilized soil reaction against the Ito &
    Matsui limiting resistance is reported. When neither exists there is nothing
    to measure and the badge stays neutral rather than guessing at one.
    """
    best, basis = None, None
    if V_cap and np.isfinite(V_cap) and V_cap > 0 and len(shear):
        r = float(np.max(np.abs(shear)) / V_cap)
        best, basis = r, "shear vs Vcap"
    if M_cap and np.isfinite(M_cap) and M_cap > 0 and len(moment):
        r = float(np.max(np.abs(moment)) / M_cap)
        if best is None or r > best:
            best, basis = r, "moment vs Mcap"
    if best is None and reaction_ratio is not None:
        best, basis = reaction_ratio, "soil reaction vs Ito-Matsui limit"
    return best, basis


def _pile_status(util, basis):
    if util is None:
        return "no capacity declared"
    if util >= UTIL_AT_CAPACITY:
        return f"at capacity ({basis})"
    if util >= UTIL_WATCH:
        return f"near capacity ({basis})"
    return f"within capacity ({basis})"


# --------------------------------------------------------------------------
# shared array helpers + units
# --------------------------------------------------------------------------

def _elem_array(fem_data, key, n, fill=0.0):
    val = (fem_data or {}).get(key, None)
    if val is None:
        return np.full(n, fill, dtype=float)
    arr = np.asarray(val, dtype=float)
    if arr.shape[:1] != (n,):
        return np.full(n, fill, dtype=float)
    return arr


def _pile_array(fem_data, key, n, fill=0.0):
    return _elem_array(fem_data, key, n, fill)


def _sol_array(solution, key, n, dtype=float, fill=0.0):
    val = (solution or {}).get(key, None)
    if val is None:
        return np.full(n, fill, dtype=dtype)
    arr = np.asarray(val)
    if arr.shape[:1] != (n,):
        return np.full(n, fill, dtype=dtype)
    return arr.astype(dtype)


def unit_labels(fem_data):
    """Display unit strings for the detail plots and the exported CSV.

    ``force`` and ``moment`` follow the FEM's per-unit-width-of-slope
    convention: reinforcement and pile results are already divided by the
    member spacing, so a force reads as force per unit length of slope and a
    moment as force (moment per unit length of slope).
    """
    from .units import labels, normalize_unit_system
    system = normalize_unit_system((fem_data or {}).get("unit_system"))
    if system is None:
        return {"length": "", "force": "", "moment": "", "line_load": "", "stress": ""}
    lb = labels(system, (fem_data or {}).get("time_unit"))
    fpl = lb.get("force_per_len", "")
    force = fpl.split("/")[0] if fpl else ""
    length = lb.get("length", "")
    return {
        "length": length,
        # Force per unit width of slope, e.g. "kN/m".
        "force": fpl,
        # Moment per unit width of slope: force x length / length. Spelled out
        # rather than cancelled to bare force, which reads as the wrong quantity.
        "moment": f"{force}·{length} per {length}" if force and length else "",
        # Force per unit length of member, per unit width of slope — the units
        # of a bond transfer rate or a soil reaction.
        "line_load": f"{force}/{length} per {length}" if force and length else "",
        "stress": lb.get("stress", ""),
    }


# --------------------------------------------------------------------------
# tabular form (the dialog's CSV export)
# --------------------------------------------------------------------------

def profile_table(profile):
    """The profile's plotted series as ``(columns, rows)`` for CSV export.

    Series sampled at different positions (element centroids, nodes, element
    boundaries) cannot share one row without inventing values, so each is
    written as its own position/value column pair. Column headers carry the
    model's units when it declares a system.
    """
    u = profile.get("units", {}) or {}
    ln = f" ({u['length']})" if u.get("length") else ""
    fo = f" ({u['force']})" if u.get("force") else ""
    mo = f" ({u['moment']})" if u.get("moment") else ""
    ll = f" ({u['line_load']})" if u.get("line_load") else ""

    if profile["kind"] == "reinforcement":
        series = [
            (f"position{ln}", profile["s"]),
            (f"axial_force{fo}", profile["T"]),
            (f"capacity{fo}", profile["t_cap"]),
            (f"residual_capacity{fo}", profile["t_res"]),
            ("utilization", profile["utilization"]),
            ("failed", profile["failed"].astype(int)),
            ("softened", profile["softened"].astype(int)),
        ]
        extra = []
        if profile.get("env_s") is not None:
            extra += [(f"envelope_position{ln}", profile["env_s"]),
                      (f"envelope_capacity{fo}", profile["env_T"])]
        extra += [(f"bond_position{ln}", profile["bond_s"]),
                  (f"bond_transfer{ll}", profile["bond_q"])]
    else:
        series = [
            (f"node_depth{ln}", profile["node_depth"]),
            (f"lateral_displacement{ln}", profile["u_lateral"]),
        ]
        extra = [
            (f"shear_depth{ln}", profile["elem_depth"]),
            (f"shear{fo}", profile["shear"]),
            (f"moment_depth{ln}", profile["moment_depth"]),
            (f"moment{mo}", profile["moment"]),
            (f"reaction_depth{ln}", profile["reaction_depth"]),
            (f"soil_reaction{ll}", profile["reaction"]),
        ]
        if profile.get("limit_p") is not None:
            extra += [(f"limit_depth{ln}", profile["limit_depth"]),
                      (f"ito_matsui_limit{ll}", profile["limit_p"])]

    cols = series + extra
    n = max((len(v) for _, v in cols), default=0)
    columns = [c for c, _ in cols]
    rows = []
    for i in range(n):
        rows.append(["" if i >= len(v) else _fmt(v[i]) for _, v in cols])
    return columns, rows


def _fmt(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(f):
        return ""
    return f"{f:.6g}"


def profile_comment_lines(profile):
    """Provenance lines for the exported CSV, in the ``#`` comment form the FEM
    result sidecars use. The field state comes first: a reader has to be able to
    tell an at-failure profile from a converged one from the file alone."""
    return [f"# field state: {field_state_label(profile.get('field_state'))}\n"]


def write_profile_csv(profile, path):
    """Write :func:`profile_table` to ``path`` as CSV, under the provenance
    comment lines that record which field the profile was read from."""
    import csv
    columns, rows = profile_table(profile)
    with open(path, "w", newline="") as f:
        for line in profile_comment_lines(profile):
            f.write(line)
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(rows)
    return path
