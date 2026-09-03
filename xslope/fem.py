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

import os
import time
import warnings
from math import degrees, sin, cos, sqrt, asin, tan


import numpy as np
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, issparse
from scipy.sparse.linalg import splu
import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from .hoekbrown import hb_constants, hb_tangent_const
from .units import require_gamma_water


# ===================== Phase 0 profiling (SPIKE, "THE FACTORIZATION") ============
# Timers only. Nothing here reads or writes solver state, so a profiled run and an
# unprofiled one execute the same arithmetic in the same order and return the same
# numbers — verified on the 13-row representative set, where a profiled Newton
# bisection reproduces the Sweep 1 column value for value in factor of safety,
# Newton iterations, force evaluations and predictor iterations. The only
# difference is a perf_counter pair around six call sites, and it is off unless
# XSLOPE_NR_PROFILE is set, so the shipped path pays one boolean per call.
#
# Keys are (phase -> seconds) plus (n_<phase> -> count):
#   nr_assemble     the consistent tangent's bincount into the ready-made pattern
#   nr_factorize    splu on it
#   nr_trisolve     lu.solve for the Newton correction
#   nr_const_tan    the constitutive pass that also differences the moduli
#   nr_const_res    the residual-only constitutive pass at the top of an iteration
#                   (an iteration that re-uses its factorization takes this instead)
#   nr_linesearch   the line search's residual-only constitutive passes
#   vp_const / vp_trisolve / vp_factorize   the same for the viscoplastic driver,
#                   whose assembly and factorization are once per trial rather than
#                   once per iteration
_PROF_ON = bool(os.environ.get("XSLOPE_NR_PROFILE", "").strip())
_PROF = {}


def _prof_reset():
    """Start a fresh accumulation. Read it back with :func:`_prof_read`."""
    _PROF.clear()


def _prof_read():
    """A copy of the accumulated timers, safe to ship across a process boundary."""
    return dict(_PROF)


def _prof_add(key, t0, n=1):
    _PROF[key] = _PROF.get(key, 0.0) + (time.perf_counter() - t0)
    _PROF["n_" + key] = _PROF.get("n_" + key, 0) + n



def _extract_nodal_uv(disp, fem_data):
    """Extract per-node translational displacements from a mixed-DOF vector."""
    dof_offset = fem_data.get("dof_offset", None)
    if dof_offset is not None:
        n_nodes = len(fem_data["nodes"])
        u = np.array([disp[dof_offset[i]] for i in range(n_nodes)])
        v = np.array([disp[dof_offset[i] + 1] for i in range(n_nodes)])
    else:
        u = disp[0::2]
        v = disp[1::2]
    return u, v


def _extract_nodal_theta(disp, fem_data):
    """Per-node rotation from a mixed-DOF vector, or None where no node has one.

    Only pile nodes carry a rotational degree of freedom; every other node is
    reported at zero, so the column reads across the whole mesh.
    """
    dof_offset = fem_data.get("dof_offset", None)
    is_pile_node = fem_data.get("is_pile_node", None)
    if dof_offset is None or is_pile_node is None:
        return None
    is_pile_node = np.asarray(is_pile_node, dtype=bool)
    if not is_pile_node.any():
        return None
    n_nodes = len(fem_data["nodes"])
    theta = np.zeros(n_nodes)
    for i in range(n_nodes):
        if i < len(is_pile_node) and is_pile_node[i]:
            base = int(dof_offset[i]) + 2
            if base < len(disp):
                theta[i] = disp[base]
    return theta


# Scalar fields of a solve_fem field carried across a save/reload for the
# at-failure snapshot (the node/element CSVs carry the arrays). The trial F is
# the only one the result-panel titles read; the rest round-trip the snapshot's
# own metadata faithfully. Render-time markers (_ssrm_fs, _at_failure) are NOT
# stored here — plot_fem_results derives them from the SSRM FS and the presence
# of the failure field, so persisting them would only diverge from the original.
_FEM_FAILURE_META_KEYS = (
    "F", "converged", "iterations", "max_displacement", "algorithm",
    "residual", "unbalanced_force_ratio", "plastic_fraction",
)

# What a solve knows about itself that the node/element CSVs do not carry: did it
# close, in how many iterations, on what residual, and how far the section moved.
# :func:`export_fem_solution` writes them into the meta sidecar off the solution
# being exported, and :func:`import_fem_solution` restores them onto the solution
# it rebuilds — WHEN THE FILE RECORDS THEM. A file written before these lines
# existed records none, and the keys stay ABSENT: import used to set
# converged = True on every file it read, which asserted a solve had closed on
# the strength of the file existing. Unknown is not converged.
_FEM_SOLVE_META_KEYS = ("converged", "iterations", "residual",
                        "max_displacement")

# What a strength reduction RUN knows about itself, as against what one solved
# FIELD knows (_FEM_SOLVE_META_KEYS above). These are the choices and the record
# that produced the factor of safety: which criterion decided a trial, how narrow
# the search was driven, the interval it ended on, how many steps that took, and
# the per-trial verdicts themselves. solve_ssrm returns them on its RESULT dict —
# not on result['last_solution'], which is the field the CSVs are written from —
# so a writer that persisted only the field dropped every one of them, and a
# reloaded run could say nothing about how its answer was reached.
_SSRM_RUN_RESULT_KEYS = ("failure_criterion", "method", "final_interval",
                         "interval_width", "iterations_ssrm", "trials")

# The same, for the options the CALLER chose: the tolerance the search was driven
# to (solve_ssrm takes it as a kwarg and does not return it), the search range,
# and the materials held at full strength. Read off the options dict a runner
# already has.
_SSRM_RUN_OPTION_KEYS = ("tolerance", "F_min", "F_max", "ssr_exclude")

# Per-trial keys that stay OUT of a meta sidecar: the Newton cost attribution added
# for SPIKE.md's "THE COST OF THE RESCUE". They are how the answer was reached and
# how long it took on one machine, not what the answer is.
_SSRM_TRIAL_COST_KEYS = frozenset((
    "wall", "bracket", "nr_cold_wall", "nr_cold_iterations", "nr_cold_force_evals",
    "nr_rungs", "nr_cold_skipped"))


def _jsonable(value):
    """``value`` as something :mod:`json` will write: numpy scalars become plain
    numbers, numpy arrays and tuples become lists, and containers are converted
    through. Anything else is returned unchanged, for json to accept or reject on
    its own."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_jsonable(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def ssrm_run_record(result, fem_data=None, options=None):
    """What a strength reduction run chose and what its trials found, as meta.

    ``result`` is :func:`solve_ssrm`'s own return value, ``options`` the dict of
    run options the caller solved with, and ``fem_data`` the model — the strength
    reduction zone overlays are carried on it rather than on either of the other
    two. Every key is recorded only where its source carries it, so a run that
    made no choice is not credited with one and the returned dict layers onto a
    meta sidecar additively: a reader of an older file sees the same keys it
    always saw, and one of a newer file can also say how the answer was reached.

    The zones are recorded as a COUNT PER KIND rather than as polygons. What a
    reader needs from them is whether the factor of safety is a whole-section
    answer or a confined one; the polygons themselves are in the model file that
    the run was made from, and a second copy of them here would be a second
    truth about the same geometry.

    Returns a plain, JSON-writable dict — empty for a run that recorded nothing.
    """
    record = {}
    for key in _SSRM_RUN_RESULT_KEYS:
        value = (result or {}).get(key)
        if value is None:
            continue
        if key == "trials":
            # The Newton driver's cost attribution (SPIKE.md, "THE COST OF THE
            # RESCUE") is measurement, not a fact about the slope: wall times are
            # machine-dependent and the rung breakdown is about how the answer was
            # reached rather than what it is. A meta sidecar is a committed artifact,
            # so it carries neither.
            value = [{k: v for k, v in trial.items() if k not in _SSRM_TRIAL_COST_KEYS}
                     for trial in value]
        record[key] = _jsonable(value)
    for key in _SSRM_RUN_OPTION_KEYS:
        value = (options or {}).get(key)
        if value is not None and value != []:
            record[key] = _jsonable(value)

    # The zones, from wherever the run got them: an explicit search-area polygon
    # passed as ssr_zone IS a reduce zone, and is counted as one, so the record
    # says the same thing however the confinement was expressed.
    counts = {}
    for zone in ((fem_data or {}).get("ssr_zones") or []):
        kind = str(zone.get("kind", "")).strip()
        if kind:
            counts[kind] = counts.get(kind, 0) + 1
    if (options or {}).get("ssr_zone"):
        counts["reduce"] = counts.get("reduce", 0) + 1
    if counts:
        record["ssr_zones"] = counts
    return record


def _fem_solution_dataframes(fem_data, solution):
    """Build the (node_df, element_df) pair persisted for one solve_fem field.

    Shared by the converged export and the at-failure snapshot export so both
    round-trip through the identical CSV schema.
    """
    import pandas as pd

    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    element_materials = fem_data["element_materials"]

    disp_total = solution["displacements"]
    ux_total, uy_total = _extract_nodal_uv(disp_total, fem_data)
    u_mag_total = np.sqrt(ux_total**2 + uy_total**2)

    disp_elastic = solution.get("displacements_elastic")
    if disp_elastic is not None:
        disp_vp = disp_total - disp_elastic
        ux_vp, uy_vp = _extract_nodal_uv(disp_vp, fem_data)
    else:
        ux_vp = np.zeros(len(nodes))
        uy_vp = np.zeros(len(nodes))
    u_mag_vp = np.sqrt(ux_vp**2 + uy_vp**2)

    node_df = pd.DataFrame({
        "node_id": np.arange(1, len(nodes) + 1),
        "x": nodes[:, 0],
        "y": nodes[:, 1],
        "u_x": ux_total,
        "u_y": uy_total,
        "u_mag": u_mag_total,
        "u_x_vp": ux_vp,
        "u_y_vp": uy_vp,
        "u_mag_vp": u_mag_vp,
    })

    # Pile nodes carry a rotation as well as a displacement, and it is part of
    # the field: the moment and the soil reaction a pile reports are read from it.
    # Written only where the model has a pile, so a model without one exports
    # exactly the columns it always did.
    theta = _extract_nodal_theta(disp_total, fem_data)
    if theta is not None:
        node_df["theta"] = theta

    centroids = np.zeros((len(elements), 2))
    for i, elem_nodes in enumerate(elements):
        active_nodes = elem_nodes[:element_types[i]]
        coords = nodes[active_nodes]
        centroids[i] = np.mean(coords, axis=0)

    element_df = pd.DataFrame({
        "element_id": np.arange(1, len(elements) + 1),
        "material_id": element_materials,
        "x_centroid": centroids[:, 0],
        "y_centroid": centroids[:, 1],
        "sigma_x": solution["stresses"][:, 0],
        "sigma_y": solution["stresses"][:, 1],
        "tau_xy": solution["stresses"][:, 2],
        "sigma_vm": solution["stresses"][:, 3],
        "eps_x": solution["strains"][:, 0],
        "eps_y": solution["strains"][:, 1],
        "gamma_xy": solution["strains"][:, 2],
        "max_shear_strain": solution["strains"][:, 3],
        "vp_shear_strain": solution["vp_shear_strain"],
        "plastic": solution["plastic_elements"],
        "yield_function": solution["yield_function"],
    })

    return node_df, element_df


def _as_len(arr, n, dtype=float):
    """Coerce a per-element array to length ``n`` (zeros if it does not match).

    A reloaded/absent field can hand back an empty or short array; the solution
    fields are only trustworthy when their length equals the mesh's element count.
    """
    arr = np.asarray(arr)
    if arr.shape[0] == n:
        return arr.astype(dtype)
    return np.zeros(n, dtype=dtype)


def _fem_reinforcement_dataframe(fem_data, solution):
    """Build the per-reinforcement-element results DataFrame for one solve_fem
    field, or ``None`` when the model has no (non-pile) reinforcement 1D elements.

    One row per reinforcement 1D element with engineer-readable columns: the global
    element id and 1-based line id, the two endpoint coordinates, the axial force,
    the allowable and residual tensile capacities, the mobilization (force /
    T_allow), and the failed / softened flags. This doubles as a human-readable
    results file AND carries everything needed to rebuild
    ``solution['forces_1d']`` / ``['failed_1d_elements']`` / ``['softened_1d_elements']``
    for the reinforcement-force colorbar and the per-line detail profiles on reload.

    The capacity the solve enforces is ``fem_data['t_allow_by_1d_elem']``, which is
    rebuilt from the model on reload, so it is written here for the reader but never
    read back.
    """
    import pandas as pd

    elements_1d = fem_data.get("elements_1d", None)
    if elements_1d is None or len(elements_1d) == 0:
        return None
    n_1d = len(elements_1d)
    pile_elem_mask = _as_len(fem_data.get("pile_elem_mask", np.zeros(n_1d)), n_1d, bool)
    reinf_idx = np.where(~pile_elem_mask)[0]
    if len(reinf_idx) == 0:
        return None

    nodes = fem_data["nodes"]
    element_materials_1d = _as_len(
        fem_data.get("element_materials_1d", np.zeros(n_1d)), n_1d, int)
    t_allow = _as_len(fem_data.get("t_allow_by_1d_elem", np.zeros(n_1d)), n_1d)
    t_res = _as_len(fem_data.get("t_res_by_1d_elem", np.zeros(n_1d)), n_1d)
    forces = _as_len(solution.get("forces_1d", np.zeros(n_1d)), n_1d)
    failed = _as_len(solution.get("failed_1d_elements", np.zeros(n_1d)), n_1d, bool)
    softened = _as_len(solution.get("softened_1d_elements", np.zeros(n_1d)), n_1d, bool)

    n0 = np.array([nodes[elements_1d[i][0]] for i in reinf_idx])
    n1 = np.array([nodes[elements_1d[i][1]] for i in reinf_idx])
    ta = t_allow[reinf_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        mobilization = np.where(ta > 1e-12, forces[reinf_idx] / ta, 0.0)

    return pd.DataFrame({
        "element_id": reinf_idx.astype(int),
        "line_id": element_materials_1d[reinf_idx].astype(int),
        "x_start": n0[:, 0],
        "y_start": n0[:, 1],
        "x_end": n1[:, 0],
        "y_end": n1[:, 1],
        "axial_force": forces[reinf_idx],
        "t_allow": ta,
        "t_res": t_res[reinf_idx],
        "mobilization": mobilization,
        "failed": failed[reinf_idx],
        "softened": softened[reinf_idx],
    })


def _fem_pile_dataframe(fem_data, solution):
    """Build the per-pile-element results DataFrame for one solve_fem field, or
    ``None`` when the model has no pile elements.

    One row per pile beam element with engineer-readable columns: the 0-based pile
    element index (the order the pile-force arrays and the force colorbar use), the
    global 1D element id and 1-based pile-line id, the two endpoint coordinates, the
    axial / lateral(shear) forces, the two end moments, the structural shear/moment
    capacities (V_cap / M_cap, blank-as-``inf`` when uncapped), and the yielded flags.
    Doubles as a human-readable results file AND carries everything needed to rebuild
    ``solution['forces_pile_axial'/'forces_pile_lateral'/'forces_pile_moment']`` and
    the ``yielded_pile*`` masks for the pile-shear colorbar on reload.
    """
    import pandas as pd

    n_pile = fem_data.get("n_pile_elements", 0)
    if n_pile == 0:
        return None

    elements_1d = fem_data["elements_1d"]
    nodes = fem_data["nodes"]
    element_materials_1d = fem_data.get("element_materials_1d", np.array([], dtype=int))
    pile_elem_indices = _as_len(
        fem_data.get("pile_elem_indices", np.arange(n_pile)), n_pile, int)

    forces_axial = _as_len(solution.get("forces_pile_axial", np.zeros(n_pile)), n_pile)
    forces_shear = _as_len(solution.get("forces_pile_lateral", np.zeros(n_pile)), n_pile)
    fm = np.asarray(solution.get("forces_pile_moment", np.zeros((n_pile, 2))))
    forces_moment = fm if fm.shape == (n_pile, 2) else np.zeros((n_pile, 2))
    yielded_V = _as_len(solution.get("yielded_pile_V", np.zeros(n_pile)), n_pile, bool)
    yielded_M = _as_len(solution.get("yielded_pile_M", np.zeros(n_pile)), n_pile, bool)
    V_cap = _as_len(fem_data.get("V_cap_by_pile_elem", np.full(n_pile, np.inf)), n_pile)
    M_cap = _as_len(fem_data.get("M_cap_by_pile_elem", np.full(n_pile, np.inf)), n_pile)
    pr = np.asarray(solution.get("pile_plastic_rotation", np.zeros((n_pile, 2))))
    plastic_rot = pr if pr.shape == (n_pile, 2) else np.zeros((n_pile, 2))

    line_ids = np.array([element_materials_1d[pile_elem_indices[p]]
                         if len(element_materials_1d) else 0 for p in range(n_pile)],
                        dtype=int)
    n0 = np.array([nodes[elements_1d[pile_elem_indices[p]][0]] for p in range(n_pile)])
    n1 = np.array([nodes[elements_1d[pile_elem_indices[p]][1]] for p in range(n_pile)])

    return pd.DataFrame({
        "pile_index": np.arange(n_pile, dtype=int),
        "element_id": pile_elem_indices.astype(int),
        "line_id": line_ids,
        "x_start": n0[:, 0],
        "y_start": n0[:, 1],
        "x_end": n1[:, 0],
        "y_end": n1[:, 1],
        "axial_force": forces_axial,
        "shear_force": forces_shear,
        "moment_1": forces_moment[:, 0],
        "moment_2": forces_moment[:, 1],
        "v_cap": V_cap,
        "m_cap": M_cap,
        "yielded_shear": yielded_V,
        "yielded_moment": yielded_M,
        "yielded": yielded_V | yielded_M,
        # The rotation the moment capacity released at each end node. Zero
        # wherever no hinge formed, which is every element of an uncapped model.
        "plastic_rotation_1": plastic_rot[:, 0],
        "plastic_rotation_2": plastic_rot[:, 1],
    })


def _reconstruct_reinforcement(fem_data, reinf_df, solution):
    """Restore the reinforcement result arrays onto ``solution`` from the sidecar.

    ``forces_1d`` / ``failed_1d_elements`` / ``softened_1d_elements`` are rebuilt at
    full 1D-element length, each row placed at its global ``element_id``; pile slots
    stay zero — the renderer and the print summaries skip them via ``pile_elem_mask``.

    The capacity is not restored from the file: it is ``t_allow_by_1d_elem``, which
    ``build_fem_data`` rebuilds from the model.
    """
    elements_1d = fem_data.get("elements_1d", None)
    n_1d = 0 if elements_1d is None else len(elements_1d)
    forces = np.zeros(n_1d)
    failed = np.zeros(n_1d, dtype=bool)
    softened = np.zeros(n_1d, dtype=bool)

    eid = reinf_df["element_id"].to_numpy()
    f = reinf_df["axial_force"].to_numpy()
    fl = reinf_df["failed"].to_numpy().astype(bool)
    sf = reinf_df["softened"].to_numpy().astype(bool)
    for k in range(len(reinf_df)):
        j = int(eid[k])
        if 0 <= j < n_1d:
            forces[j] = f[k]
            failed[j] = fl[k]
            softened[j] = sf[k]

    solution["forces_1d"] = forces
    solution["failed_1d_elements"] = failed
    solution["softened_1d_elements"] = softened


def _reconstruct_piles(fem_data, pile_df, solution):
    """Restore the pile result arrays onto ``solution`` from the sidecar.

    The pile-force arrays are indexed by the 0-based pile element index (the same
    order the renderer walks pile elements), so each row is placed at its
    ``pile_index``.
    """
    n_pile = fem_data.get("n_pile_elements", 0)
    axial = np.zeros(n_pile)
    shear = np.zeros(n_pile)
    moment = np.zeros((n_pile, 2))
    plastic_rot = np.zeros((n_pile, 2))
    yV = np.zeros(n_pile, dtype=bool)
    yM = np.zeros(n_pile, dtype=bool)
    # Written since the moment capacity became a hinge; a sidecar from before
    # that carries no hinge and reads as zero, which is what it was.
    has_pr = ("plastic_rotation_1" in pile_df.columns
              and "plastic_rotation_2" in pile_df.columns)

    pidx = pile_df["pile_index"].to_numpy()
    for k in range(len(pile_df)):
        p = int(pidx[k])
        if 0 <= p < n_pile:
            axial[p] = pile_df["axial_force"].iloc[k]
            shear[p] = pile_df["shear_force"].iloc[k]
            moment[p, 0] = pile_df["moment_1"].iloc[k]
            moment[p, 1] = pile_df["moment_2"].iloc[k]
            yV[p] = bool(pile_df["yielded_shear"].iloc[k])
            yM[p] = bool(pile_df["yielded_moment"].iloc[k])
            if has_pr:
                plastic_rot[p, 0] = pile_df["plastic_rotation_1"].iloc[k]
                plastic_rot[p, 1] = pile_df["plastic_rotation_2"].iloc[k]

    solution["forces_pile_axial"] = axial
    solution["forces_pile_lateral"] = shear
    solution["forces_pile_moment"] = moment
    solution["pile_plastic_rotation"] = plastic_rot
    solution["yielded_pile_V"] = yV
    solution["yielded_pile_M"] = yM
    solution["yielded_pile"] = yV | yM


def _write_units_header(f, fem_data):
    """Write a ``# units:`` provenance line to an open FEM result CSV -- only when the
    model declares a unit system; a no-op otherwise, so undeclared exports stay
    byte-identical. Readers skip it as an ordinary ``#`` comment."""
    from .units import units_comment_line
    header = units_comment_line(fem_data.get("unit_system"), fem_data.get("time_unit"))
    if header:
        f.write(header)


def _write_per_width_note(f, fem_data):
    """Note the per-unit-width force/capacity convention at the top of a reinforcement
    or pile result CSV -- only when the model declares a unit system, so undeclared
    exports stay byte-identical. An ordinary ``#`` comment readers skip (the sidecar
    reloaders read with ``comment='#'``)."""
    from .units import normalize_unit_system
    if normalize_unit_system(fem_data.get("unit_system")) is None:
        return
    f.write("# forces and capacities are per unit width of slope "
            "(per-element value divided by Spacing / S)\n")


def _write_1d_result_sidecars(fem_data, solution, output_stem, tag):
    """Write the reinforcement / pile per-element result CSVs for one solve_fem
    field. ``tag`` is ``"fem"`` for the converged field or ``"fem_failure"`` for the
    at-failure twin. Each CSV is written only when the model actually has that
    element type. Returns a list of ``(kind, path)`` for the caller to announce.
    """
    written = []
    reinf_df = _fem_reinforcement_dataframe(fem_data, solution)
    if reinf_df is not None:
        path = output_stem.parent / f"{output_stem.name}_{tag}_reinf.csv"
        with open(path, "w") as f:
            _write_units_header(f, fem_data)
            _write_per_width_note(f, fem_data)
            reinf_df.to_csv(f, index=False)
        written.append(("reinforcement", path))
    pile_df = _fem_pile_dataframe(fem_data, solution)
    if pile_df is not None:
        path = output_stem.parent / f"{output_stem.name}_{tag}_piles.csv"
        with open(path, "w") as f:
            _write_units_header(f, fem_data)
            _write_per_width_note(f, fem_data)
            pile_df.to_csv(f, index=False)
        written.append(("pile", path))
    return written


def _1d_sidecar_mismatch(name, df, id_column, n_rows, n_slots, kind,
                         expected=None, also=()):
    """Why a 1D result sidecar does not belong to this model, or ``None``.

    The file is this model's only if it holds one row per element of ``kind`` and
    every row addresses an element OF THAT KIND that the model has. Three things
    are checked, and each catches what the others let through:

    * the row count — a file with too few rows covers a fraction of the members;
    * the id RANGE — a file whose ids have all shifted addresses elements the
      model does not have;
    * the id SET, which is the real test. A model carrying both reinforcement and
      piles splits one 1D element list between the two kinds, and the ids of one
      kind are exactly the slots that kind occupies. A range check accepts a
      reinforcement file whose rows land on the pile slots — same count, every id
      inside the list, every force grafted onto the wrong element — because both
      kinds index the same array. Set membership does not.

    ``n_rows`` is how many rows the model's own export writes for ``kind``;
    ``n_slots`` is the size of the array the ids index (the two differ for
    reinforcement, whose ids are positions in the full 1D element list).
    ``expected`` is the exact set of ids the model's own export would write under
    ``id_column``; ``also`` is a sequence of ``(column, ids)`` pairs holding a
    second id column of the same file to its own set, checked only where the file
    carries that column (older sidecars do not).
    """
    if len(df) != n_rows:
        return (f"{name} records {len(df)} {kind} "
                f"{'result' if len(df) == 1 else 'results'}; the model has "
                f"{n_rows}, so the saved member forces are not this model's and "
                f"were not restored.")
    if len(df) == 0:
        return None
    ids = df[id_column].to_numpy()
    out = ids[(ids < 0) | (ids >= n_slots)]
    if len(out):
        return (f"{name} places a {kind} result at element {int(out[0])}, which "
                f"is outside the {n_slots} the model carries, so the saved member "
                f"forces are not this model's and were not restored.")
    for column, want in ((id_column, expected),) + tuple(also):
        if want is None or column not in df.columns:
            continue
        have = np.asarray(df[column].to_numpy(), dtype=int)
        want = np.asarray(want, dtype=int)
        wrong = np.setdiff1d(have, want)
        missing = np.setdiff1d(want, have)
        if len(wrong) or len(missing):
            stray = int(wrong[0]) if len(wrong) else int(missing[0])
            return (f"{name} addresses element {stray} under {column}, which is "
                    f"not one of the {len(want)} the model's {kind} elements "
                    f"occupy, so the saved member forces are not this model's "
                    f"and were not restored.")
    return None


def _import_1d_result_sidecars(fem_data, solution, output_stem, tag):
    """Restore reinforcement / pile results onto ``solution`` from the ``tag``
    sidecars (``"fem"`` converged / ``"fem_failure"`` twin) when they are present.
    A no-op — leaving the solution unchanged — when the files are absent, so
    solutions saved before these sidecars existed import cleanly.

    A file that is not this model's is refused rather than grafted. The reloader
    reads whatever ``{stem}_{tag}_reinf.csv`` and ``{stem}_{tag}_piles.csv`` sit
    beside the field, and those names carry no model identity: a reinforced model
    reopened next to another model's members would have taken them, and rows
    addressing elements this model does not have were dropped one at a time,
    leaving a partial set of forces that reads as a solved result. Each file is
    checked against the model's own element count before anything is placed, and a
    file that fails is left out whole. The refusal is recorded on the solution
    under ``sidecar_notes``, which is where the report drains its notes from — so a
    report of the model says the member forces were not restored instead of
    printing someone else's.
    """
    import pandas as pd

    elements_1d = fem_data.get("elements_1d", None)
    n_1d = 0 if elements_1d is None else len(elements_1d)
    pile_mask = _as_len(fem_data.get("pile_elem_mask", np.zeros(n_1d)), n_1d, bool)
    n_reinf = int(np.count_nonzero(~pile_mask)) if n_1d else 0
    n_pile = int(fem_data.get("n_pile_elements", 0))

    # The slots each kind occupies in the model's own 1D element list. These are
    # exactly the ids the model's own export writes, so a file whose ids are any
    # other set is not this model's — including one whose rows land on the OTHER
    # kind's slots, which a count-and-range check cannot tell apart.
    reinf_ids = np.flatnonzero(~pile_mask) if n_1d else np.array([], dtype=int)
    pile_ids = np.asarray(
        fem_data.get("pile_elem_indices", np.arange(n_pile)), dtype=int)[:n_pile]

    notes = []
    reinf_path = output_stem.parent / f"{output_stem.name}_{tag}_reinf.csv"
    if reinf_path.exists():
        df = pd.read_csv(reinf_path, comment="#")
        bad = _1d_sidecar_mismatch(reinf_path.name, df, "element_id",
                                   n_reinf, n_1d, "reinforcement",
                                   expected=reinf_ids)
        if bad:
            notes.append(bad)
        else:
            _reconstruct_reinforcement(fem_data, df, solution)
    pile_path = output_stem.parent / f"{output_stem.name}_{tag}_piles.csv"
    if pile_path.exists():
        df = pd.read_csv(pile_path, comment="#")
        # pile_index is the position in the pile-force arrays; element_id is the
        # slot in the 1D element list. The first is 0..n-1 by construction, so it
        # is the second that identifies the model — held here where the file
        # carries it.
        bad = _1d_sidecar_mismatch(pile_path.name, df, "pile_index",
                                   n_pile, n_pile, "pile",
                                   expected=np.arange(n_pile),
                                   also=(("element_id", pile_ids),))
        if bad:
            notes.append(bad)
        else:
            _reconstruct_piles(fem_data, df, solution)
    if notes:
        solution.setdefault("sidecar_notes", []).extend(notes)


def export_fem_solution(fem_data, solution, output_stem, meta=None,
                        failure_solution=None):
    """Export FEM nodal and element results to CSV files using a common stem.

    If ``meta`` (a dict) is given, it is also written to ``{stem}_fem_meta.json`` —
    use it for run metadata that is not in the node/element CSVs, e.g. the SSRM
    factor of safety and the analysis type, so they survive a reload. The unit
    declaration and the solve facts ``solution`` carries (``converged``,
    ``iterations``, ``residual``, ``max_displacement``) are added to it here, so
    every writer records them without naming them; a key the caller already set
    wins, and one the solution does not carry is left out.

    If ``failure_solution`` (a solve_fem field, i.e. ``result['failure_solution']``
    captured by :func:`solve_ssrm`) is given, the at-failure mechanism is persisted
    alongside the converged solution as a second CSV pair,
    ``{stem}_fem_failure_nodes.csv`` / ``{stem}_fem_failure_elements.csv`` (same
    schema), plus its scalar metadata in ``{stem}_fem_failure_meta.json`` (the trial
    F and the snapshot's own diagnostics). These let a reloaded solution re-render
    the deformation / displacement-vector / failure-state contour panels from the
    mechanism instead of the sub-critical last-converged field. When it is ``None``
    (a single solve, or an SSRM run with no capture) nothing extra is written.

    When the model carries reinforcement and/or pile 1D elements, the per-element
    structural results are ALSO written as engineer-readable CSVs —
    ``{stem}_fem_reinf.csv`` (per reinforcement bar: line/element ids, endpoints,
    axial force, capacities, mobilization, failed/softened flags) and
    ``{stem}_fem_piles.csv`` (per pile beam element: ids, endpoints, axial/shear
    forces, end moments, V/M capacities, yielded flags) — plus their at-failure twins
    ``{stem}_fem_failure_reinf.csv`` / ``{stem}_fem_failure_piles.csv`` when a failure
    snapshot is given. These double as results files for reading AND let a reloaded
    solution re-render the reinforcement-force / pile-shear colorbars solve-free.
    They are written only when the corresponding element type is present.
    """
    from pathlib import Path

    output_stem = Path(output_stem)
    nodes_file = output_stem.parent / f"{output_stem.name}_fem_nodes.csv"
    elements_file = output_stem.parent / f"{output_stem.name}_fem_elements.csv"

    node_df, element_df = _fem_solution_dataframes(fem_data, solution)

    with open(nodes_file, "w") as f:
        _write_units_header(f, fem_data)
        node_df.to_csv(f, index=False)
    with open(elements_file, "w") as f:
        _write_units_header(f, fem_data)
        element_df.to_csv(f, index=False)

    print(f"Exported FEM nodal results to {nodes_file}")
    print(f"Exported FEM element results to {elements_file}")

    for kind, path in _write_1d_result_sidecars(fem_data, solution, output_stem, "fem"):
        print(f"Exported FEM {kind} results to {path}")

    if failure_solution is not None:
        import json
        f_nodes_file = output_stem.parent / f"{output_stem.name}_fem_failure_nodes.csv"
        f_elements_file = output_stem.parent / f"{output_stem.name}_fem_failure_elements.csv"
        f_meta_file = output_stem.parent / f"{output_stem.name}_fem_failure_meta.json"

        f_node_df, f_element_df = _fem_solution_dataframes(fem_data, failure_solution)
        with open(f_nodes_file, "w") as f:
            _write_units_header(f, fem_data)
            f_node_df.to_csv(f, index=False)
        with open(f_elements_file, "w") as f:
            _write_units_header(f, fem_data)
            f_element_df.to_csv(f, index=False)

        f_meta = {}
        for key in _FEM_FAILURE_META_KEYS:
            if key in failure_solution:
                val = failure_solution[key]
                # np scalars (max_displacement, residual, plastic_fraction, …) are
                # not JSON-serializable — coerce to plain Python.
                if isinstance(val, np.generic):
                    val = val.item()
                f_meta[key] = val
        with open(f_meta_file, "w") as f:
            json.dump(f_meta, f, indent=2)

        print(f"Exported FEM at-failure nodal results to {f_nodes_file}")
        print(f"Exported FEM at-failure element results to {f_elements_file}")

        for kind, path in _write_1d_result_sidecars(
                fem_data, failure_solution, output_stem, "fem_failure"):
            print(f"Exported FEM at-failure {kind} results to {path}")

    if meta is not None:
        import json
        # Record the unit declaration as provenance (units plan phase 2, sec 2a). Pulled
        # from fem_data, which carries it when the model declared a system; omitted when
        # None (matching the existing meta style), and never overriding a value the
        # caller already put in meta.
        meta = dict(meta)
        for key in ("unit_system", "time_unit"):
            val = fem_data.get(key)
            if val is not None and key not in meta:
                meta[key] = val
        # The solve facts (_FEM_SOLVE_META_KEYS), off the field being written.
        # Same rule as the unit declaration: recorded where the solve knows them,
        # never overriding a value the caller put in meta, and omitted where the
        # solution carries none — so import can tell "did not converge" from
        # "does not say". np scalars are not JSON-serializable.
        for key in _FEM_SOLVE_META_KEYS:
            val = solution.get(key)
            if val is None or key in meta:
                continue
            if isinstance(val, np.generic):
                val = val.item()
            meta[key] = val
        meta_file = output_stem.parent / f"{output_stem.name}_fem_meta.json"
        with open(meta_file, "w") as f:
            json.dump(meta, f, indent=2)
        print(f"Exported FEM run metadata to {meta_file}")


def import_fem_meta(output_stem):
    """Read the ``{stem}_fem_meta.json`` sidecar written by
    :func:`export_fem_solution`, or ``None`` if it is absent or unreadable. Carries
    run metadata not stored in the node/element CSVs (FS, analysis type)."""
    import json
    from pathlib import Path

    output_stem = Path(output_stem)
    meta_file = output_stem.parent / f"{output_stem.name}_fem_meta.json"
    if not meta_file.exists():
        return None
    try:
        with open(meta_file) as f:
            return json.load(f)
    except Exception:
        return None


def _reconstruct_fem_solution(fem_data, node_df, element_df):
    """Rebuild a solve_fem field dict from one persisted node/element CSV pair.

    Only the quantities the result plots need are restored: total and elastic
    displacements (rebuilt into the mixed-DOF vector via ``dof_offset``, including
    the pile nodes' rotations where the CSV carries a ``theta`` column — a file
    written without one leaves them zero, which is what the plots, reading only
    translations, always saw), element stresses/strains, vp shear strain, the
    plastic mask, and the yield function.

    Raises:
        ValueError: if the CSV node/element counts do not match ``fem_data``.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    if len(node_df) != len(nodes):
        raise ValueError(
            f"FEM solution has {len(node_df)} nodes but the mesh has {len(nodes)} — "
            "the saved solution does not match this mesh.")
    if len(element_df) != len(elements):
        raise ValueError(
            f"FEM solution has {len(element_df)} elements but the mesh has "
            f"{len(elements)} — the saved solution does not match this mesh.")

    dof_offset = fem_data.get("dof_offset", None)
    n_dof = int(dof_offset[len(nodes)]) if dof_offset is not None else 2 * len(nodes)
    disp_total = np.zeros(n_dof)
    disp_vp = np.zeros(n_dof)
    ux, uy = node_df["u_x"].to_numpy(), node_df["u_y"].to_numpy()
    uxv, uyv = node_df["u_x_vp"].to_numpy(), node_df["u_y_vp"].to_numpy()
    theta = (node_df["theta"].to_numpy() if "theta" in node_df.columns
             else None)
    is_pile_node = fem_data.get("is_pile_node", None)
    for i in range(len(nodes)):
        base = int(dof_offset[i]) if dof_offset is not None else 2 * i
        disp_total[base], disp_total[base + 1] = ux[i], uy[i]
        disp_vp[base], disp_vp[base + 1] = uxv[i], uyv[i]
        if (theta is not None and is_pile_node is not None
                and i < len(is_pile_node) and is_pile_node[i]):
            disp_total[base + 2] = theta[i]

    return {
        "displacements": disp_total,
        "displacements_elastic": disp_total - disp_vp,
        "stresses": element_df[["sigma_x", "sigma_y", "tau_xy", "sigma_vm"]].to_numpy(),
        "strains": element_df[["eps_x", "eps_y", "gamma_xy", "max_shear_strain"]].to_numpy(),
        "vp_shear_strain": element_df["vp_shear_strain"].to_numpy(),
        "plastic_elements": element_df["plastic"].to_numpy(),
        "yield_function": element_df["yield_function"].to_numpy(),
    }


def import_fem_solution(fem_data, output_stem):
    """Reconstruct an FEM ``solution`` dict from the CSV pair written by
    :func:`export_fem_solution` — the inverse operation. Lets a previously saved
    solution be re-plotted (``plot_fem_results``) without re-running solve_fem.

    When the at-failure sidecars (``{stem}_fem_failure_nodes.csv`` /
    ``{stem}_fem_failure_elements.csv`` / ``{stem}_fem_failure_meta.json``) written
    by :func:`export_fem_solution` are present, the captured mechanism is
    reconstructed the same way and attached to the returned dict under
    ``"failure_solution"`` — the shape ``plot_fem_results`` consumes to render the
    at-failure deformation / vector / contour panels, and the metadata (trial F,
    convergence diagnostics) is restored from the failure meta sidecar. Files
    written before the failure snapshot existed simply lack these sidecars, so the
    returned dict has no ``"failure_solution"`` key and the plots fall back to the
    converged field — backward compatible in both directions.

    When the reinforcement / pile result sidecars (``{stem}_fem_reinf.csv`` /
    ``{stem}_fem_piles.csv``, and their ``_fem_failure_`` twins) are present, the
    per-element structural results are restored onto the matching solution dict where
    the renderers expect them (``forces_1d`` / ``failed_1d_elements`` /
    ``softened_1d_elements`` for the reinforcement-force colorbar; the
    ``forces_pile_*`` / ``yielded_pile*`` arrays for the pile-shear colorbar), so a
    reloaded reinforced/piled solution re-renders those overlays solve-free. Absent
    sidecars are a no-op — backward compatible in both directions. A member sidecar
    whose element count is not this model's is refused whole and noted under
    ``"sidecar_notes"`` rather than grafted (see
    :func:`_import_1d_result_sidecars`).

    The solve facts the meta sidecar records — ``converged``, ``iterations``,
    ``residual``, ``max_displacement`` — are restored onto the returned dict.
    Keys the file does not record are left ABSENT rather than guessed, so a
    consumer can tell a solve that did not close from one that never said.

    Raises:
        ValueError: if the file node/element counts do not match ``fem_data``.
    """
    import json
    import pandas as pd
    from pathlib import Path

    output_stem = Path(output_stem)
    nodes_file = output_stem.parent / f"{output_stem.name}_fem_nodes.csv"
    elements_file = output_stem.parent / f"{output_stem.name}_fem_elements.csv"
    # comment="#" skips the "# units:" provenance header when present (older sidecars
    # have none, so the read is unchanged for them).
    node_df = pd.read_csv(nodes_file, comment="#")
    element_df = pd.read_csv(elements_file, comment="#")

    solution = _reconstruct_fem_solution(fem_data, node_df, element_df)
    # The solve facts, from the meta sidecar that recorded them. A file that
    # records none leaves them absent — see _FEM_SOLVE_META_KEYS.
    meta = import_fem_meta(output_stem) or {}
    for key in _FEM_SOLVE_META_KEYS:
        if meta.get(key) is not None:
            solution[key] = meta[key]
    _import_1d_result_sidecars(fem_data, solution, output_stem, "fem")

    f_nodes_file = output_stem.parent / f"{output_stem.name}_fem_failure_nodes.csv"
    f_elements_file = output_stem.parent / f"{output_stem.name}_fem_failure_elements.csv"
    if f_nodes_file.exists() and f_elements_file.exists():
        f_node_df = pd.read_csv(f_nodes_file, comment="#")
        f_element_df = pd.read_csv(f_elements_file, comment="#")
        failure_solution = _reconstruct_fem_solution(fem_data, f_node_df, f_element_df)
        f_meta_file = output_stem.parent / f"{output_stem.name}_fem_failure_meta.json"
        if f_meta_file.exists():
            try:
                with open(f_meta_file) as f:
                    failure_solution.update(json.load(f))
            except Exception:
                pass
        # The snapshot is by definition the UNCONVERGED at-failure field; keep that
        # honest even if a stale/absent meta sidecar leaves it unset.
        failure_solution.setdefault("converged", False)
        _import_1d_result_sidecars(fem_data, failure_solution, output_stem, "fem_failure")
        solution["failure_solution"] = failure_solution

    return solution


# A domain-ring edge counts as a left/right SIDE face (x-rollers, build_fem_data
# step 3) when its two ends sit within this fraction of the domain span of the
# extreme x. It bounds the digitizing drift of a boundary meant to be vertical —
# rs2_28's truncation face drifts 1.3e-3 of the span, rs2_64b's 7e-4 — not the
# inclination of a real cut face, which runs far further than that in x.
_SIDE_EDGE_DRIFT_FRAC = 5e-3


def _midside_1d_node(elements_1d, element_types_1d, elem_idx):
    """The midside node of 1D element ``elem_idx``, or None where it has none.

    A 1D element on a quadratic soil edge stands on that edge's midside node as
    well as its two corners; ``element_types_1d`` records 3 for those and 2 for
    an element on a linear edge, whose third column is a padding zero. Node 0 is
    a real node, so the node count decides, never the stored value.
    """
    try:
        if int(element_types_1d[elem_idx]) < 3:
            return None
    except (IndexError, TypeError, ValueError):
        return None
    element = elements_1d[elem_idx]
    if len(element) < 3:
        return None
    return int(element[2])


def _quintic_hermite_bending_matrix():
    """The 6x6 bending stiffness of a three-node C1 beam element on the unit
    interval, in the DOF order [v1, dv1, v2, dv2, v3, dv3].

    The deflection is the quintic that matches a value and a slope at both ends
    AND at the midside node -- six conditions, six coefficients, so the shape
    functions are unique. Euler-Bernoulli is kept exactly: the entries are the
    integrals of EI H_i'' H_j'' with EI and the length divided out, and the
    quintic space contains both the cubic and the quartic exact solutions the
    beam benchmarks are written for, so a three-node element reproduces them to
    round-off just as the two-node cubic element does.

    The element of length L with rotations measured in x follows by scaling:
    K_bending = EI / L**3 * S @ M @ S with S = diag(1, L, 1, L, 1, L).
    Integrated with four-point Gauss, which is exact for the degree-six
    integrand.
    """
    # p(xi) = sum_k a_k xi**k, k = 0..5. Rows of C are the six conditions in DOF
    # order, so a = C^-1 q and the shape functions are the columns of C^-1.
    def row_value(x):
        return np.array([1.0, x, x**2, x**3, x**4, x**5])

    def row_slope(x):
        return np.array([0.0, 1.0, 2*x, 3*x**2, 4*x**3, 5*x**4])

    C = np.array([row_value(0.0), row_slope(0.0),
                  row_value(1.0), row_slope(1.0),
                  row_value(0.5), row_slope(0.5)])
    C_inv = np.linalg.inv(C)

    # Four-point Gauss-Legendre on [0, 1].
    g = np.array([-0.8611363115940526, -0.3399810435848563,
                  0.3399810435848563, 0.8611363115940526])
    w = np.array([0.3478548451374538, 0.6521451548625461,
                  0.6521451548625461, 0.3478548451374538])
    xi = 0.5 * (g + 1.0)
    wt = 0.5 * w

    M = np.zeros((6, 6))
    for x, weight in zip(xi, wt):
        b = np.array([0.0, 0.0, 2.0, 6*x, 12*x**2, 20*x**3])   # p''(xi)
        h2 = C_inv.T @ b                                        # H_i''(xi)
        M += weight * np.outer(h2, h2)
    return M


def _quintic_hermite_third_derivative(xi):
    """d3H/dxi3 of the three-node C1 beam shape functions at ``xi``, in the DOF
    order [v1, dv1, v2, dv2, v3, dv3].

    The shear at a station is EI v'''(x) = EI / L**3 * (d3H/dxi3) . q_hat, with
    q_hat the DOFs scaled to the unit interval. On a two-node cubic element that
    quantity is the element's constant shear, so reading it at the element center
    is the same measurement carried onto the three-node element, where the shear
    varies along the element.
    """
    def row_value(x):
        return np.array([1.0, x, x**2, x**3, x**4, x**5])

    def row_slope(x):
        return np.array([0.0, 1.0, 2*x, 3*x**2, 4*x**3, 5*x**4])

    C = np.array([row_value(0.0), row_slope(0.0),
                  row_value(1.0), row_slope(1.0),
                  row_value(0.5), row_slope(0.5)])
    C_inv = np.linalg.inv(C)
    b3 = np.array([0.0, 0.0, 0.0, 6.0, 24.0 * xi, 60.0 * xi**2])
    return C_inv.T @ b3


def _quintic_hermite_fourth_derivative(xi):
    """d4H/dxi4 of the three-node C1 beam shape functions at ``xi``, in the DOF
    order [v1, dv1, v2, dv2, v3, dv3].

    EI times the fourth derivative of the deflection is the distributed lateral
    load the element's own deflected shape carries -- what the soil is passing to
    the pile along it, per unit length. For the quintic that is linear along the
    element and generally not zero; for a two-node cubic element it is identically
    zero, which is why a chain of them can only report the reaction by
    differencing one element's shear against the next.
    """
    def row_value(x):
        return np.array([1.0, x, x**2, x**3, x**4, x**5])

    def row_slope(x):
        return np.array([0.0, 1.0, 2*x, 3*x**2, 4*x**3, 5*x**4])

    C = np.array([row_value(0.0), row_slope(0.0),
                  row_value(1.0), row_slope(1.0),
                  row_value(0.5), row_slope(0.5)])
    C_inv = np.linalg.inv(C)
    b4 = np.array([0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * xi])
    return C_inv.T @ b4


#: The constants the three-node beam element is built from, computed once.
_QUINTIC_BENDING_M = _quintic_hermite_bending_matrix()
_QUINTIC_SHEAR_MID = _quintic_hermite_third_derivative(0.5)
_QUINTIC_LOAD_MID = _quintic_hermite_fourth_derivative(0.5)


def pile_element_reaction(u_elem, cos_t, sin_t, L, EI, n_node, p_rot=None):
    """The lateral soil reaction per unit length at a pile element's center, or
    None on an element that cannot report one.

    A three-node beam element carries its own distributed load: EI times the
    fourth derivative of its deflected shape, which is linear along the element
    and read here at its midpoint. The sign follows the element's local
    transverse axis, the same sense its shear is reported in.

    A two-node element's deflection is a cubic, whose fourth derivative is zero
    everywhere, so it can say nothing about the load along itself -- a chain of
    them reports the reaction at the nodes BETWEEN elements instead, by
    differencing their shears.

    ``p_rot`` is the element's plastic end rotations where a moment capacity has
    hinged it. Only the ELASTIC part of the shape carries EI times a curvature, so
    the hinge rotation is taken out first; read on the raw nodal displacement, a
    hinge reads as a kink and the fourth derivative of a kink is a spike that no
    soil is applying.
    """
    if n_node != 3:
        return None
    u_local = _beam_rotation(cos_t, sin_t, 3) @ u_elem
    if p_rot is not None:
        u_local = u_local.copy()
        u_local[2] -= float(p_rot[0])
        u_local[5] -= float(p_rot[1])
    q_hat = np.array([u_local[1], L * u_local[2],
                      u_local[4], L * u_local[5],
                      u_local[7], L * u_local[8]])
    return -float(EI / (L ** 4) * (_QUINTIC_LOAD_MID @ q_hat))

#: Quadratic bar axial stiffness on the unit interval, node order
#: [end 1, end 2, mid]. K_axial = EA / (3 L) * this.
_QUADRATIC_BAR_AXIAL = np.array([[ 7.0,  1.0, -8.0],
                                 [ 1.0,  7.0, -8.0],
                                 [-8.0, -8.0, 16.0]])


def _beam_local_stiffness(EA, EI, L, three_node):
    """The pile beam element's local stiffness matrix.

    Two-node: the 6x6 Euler-Bernoulli matrix in [u1, v1, th1, u2, v2, th2].

    Three-node: 9x9 in [u1, v1, th1, u2, v2, th2, u3, v3, th3] with node 3 the
    midside node of the soil edge the element lies on. Bending is the quintic
    Hermite block above, axial is the quadratic bar, and the two are uncoupled
    exactly as they are on the two-node element.
    """
    L2 = L * L
    L3 = L2 * L
    if not three_node:
        return np.array([
            [ EA/L,   0.0,          0.0,        -EA/L,  0.0,          0.0       ],
            [ 0.0,    12*EI/L3,     6*EI/L2,     0.0,  -12*EI/L3,    6*EI/L2   ],
            [ 0.0,    6*EI/L2,      4*EI/L,      0.0,  -6*EI/L2,     2*EI/L    ],
            [-EA/L,   0.0,          0.0,         EA/L,   0.0,         0.0       ],
            [ 0.0,   -12*EI/L3,    -6*EI/L2,     0.0,   12*EI/L3,   -6*EI/L2   ],
            [ 0.0,    6*EI/L2,      2*EI/L,      0.0,  -6*EI/L2,     4*EI/L    ],
        ])

    K = np.zeros((9, 9))
    axial_dofs = [0, 3, 6]            # u at end 1, end 2, mid
    K[np.ix_(axial_dofs, axial_dofs)] = (EA / (3.0 * L)) * _QUADRATIC_BAR_AXIAL

    scale = np.diag([1.0, L, 1.0, L, 1.0, L])
    K_bend = (EI / L3) * (scale @ _QUINTIC_BENDING_M @ scale)
    bend_dofs = [1, 2, 4, 5, 7, 8]    # (v, theta) at end 1, end 2, mid
    K[np.ix_(bend_dofs, bend_dofs)] = K_bend
    return K


def _pile_element_actions(u_elem, cos_t, sin_t, L, EA, EI, K_local, n_node,
                          u_local=None):
    """``(axial, V, M1, M2, u_local)`` for one pile beam element.

    ``axial`` is the axial force at the element center and ``V`` the shear there;
    ``M1`` and ``M2`` are the bending moments at the two END nodes, which is what
    the moment profile is drawn from and what the plastic-hinge check compares
    against M_cap.

    On the two-node element the shear is constant along the element and these are
    the closed forms the element has always used. On the three-node element the
    end moments are rows 2 and 5 of K_local u_local -- the same quantity, read off
    the assembled element rather than off a transcribed row -- and the shear is
    EI v'''(L/2), which on a two-node element is exactly that constant shear.

    ``u_local`` may be passed in place of ``u_elem`` when the caller already holds
    the local displacement -- which is what the capacity check does, since it has
    to read the actions a second time off the ELASTIC part of the displacement
    once a plastic hinge has taken its share.
    """
    if u_local is None:
        T = _beam_rotation(cos_t, sin_t, n_node)
        u_local = T @ u_elem

    # Axial force at the element center: EA times the chord strain, on both the
    # linear and the quadratic bar.
    axial = EA / L * (u_local[3] - u_local[0])

    if n_node == 2:
        L2 = L * L
        L3 = L2 * L
        V = (12 * EI / L3 * (u_local[1] - u_local[4])
             + 6 * EI / L2 * (u_local[2] + u_local[5]))
        M1 = EI * (6.0 / L2 * u_local[1] + 4.0 / L * u_local[2]
                   - 6.0 / L2 * u_local[4] + 2.0 / L * u_local[5])
        M2 = EI * (6.0 / L2 * u_local[1] + 2.0 / L * u_local[2]
                   - 6.0 / L2 * u_local[4] + 4.0 / L * u_local[5])
        return axial, V, M1, M2, u_local

    S = K_local @ u_local
    M1, M2 = float(S[2]), float(S[5])
    q_hat = np.array([u_local[1], L * u_local[2],
                      u_local[4], L * u_local[5],
                      u_local[7], L * u_local[8]])
    V = float(EI / (L * L * L) * (_QUINTIC_SHEAR_MID @ q_hat))
    return axial, V, M1, M2, u_local


def _pile_moment_hinge(M1, M2, M_cap, K_local, allowed=None):
    """``(p, m)`` for one beam element at its moment capacity: the plastic end
    rotations and the end moments that survive them.

    A moment capacity is a RELEASE of rotational continuity, not a moment applied
    at a node. Two beam elements meet at every interior node of a pile, and at
    equilibrium their end moments there are equal and opposite; a correction
    applied to the rotational degree of freedom they SHARE is therefore equal and
    opposite too, cancels exactly, and enforces nothing at all. What a plastic
    hinge does instead is let the element end rotate freely once its moment
    reaches ``M_cap``: the element's elastic end rotation is the nodal rotation
    less a plastic rotation ``p``, so the moment it delivers is
    ``K_local (u_local - p)`` and the pile's moment diagram is bounded by
    equilibrium rather than by clipping what is reported.

    Each end moment is linear in ``p``, so ``p`` is solved for directly rather
    than iterated to. With ``A`` the rotational block of ``K_local`` (its rows and
    columns 2 and 5, the two END nodes) and ``m_e`` the elastic end moments, the
    delivered moments are ``m_e - A p``, and requiring the hinged ends to sit on
    the capacity is one 1x1 or 2x2 solve. Releasing one end raises the moment at
    the other, so the hinge set is grown until it is stable -- at most two passes,
    there being two ends. ``p`` is zero at an end that has not hinged.

    ``allowed`` masks the two ends. A hinge is ONE release at one node, and every
    interior node of a pile has two element ends on it: releasing both would leave
    the node's rotation carrying equal and opposite capacities whatever it does,
    unconstrained by the beam, and the iteration drifts along it without ever
    settling. Exactly one end per node is therefore released (see
    ``_pile_hinge_ends``); the other end is left elastic, and equilibrium at the
    node it shares puts it on the capacity too, with the opposite sign -- which is
    what a hinge means and what keeps the node's rotation determined.
    """
    m_e = np.array([float(M1), float(M2)])
    p = np.zeros(2)
    can = np.ones(2, dtype=bool) if allowed is None else np.asarray(allowed, dtype=bool)
    active = (np.abs(m_e) > M_cap) & can
    if not active.any():
        return p, m_e

    A = K_local[np.ix_([2, 5], [2, 5])]
    m = m_e
    for _ in range(2):
        idx = np.flatnonzero(active)
        target = np.sign(m_e[idx]) * M_cap
        p = np.zeros(2)
        p[idx] = np.linalg.solve(A[np.ix_(idx, idx)], m_e[idx] - target)
        m = m_e - A @ p
        newly_over = (np.abs(m) > M_cap * (1.0 + 1e-12)) & ~active & can
        if not newly_over.any():
            break
        active |= newly_over
    return p, m


def _pile_hinge_ends(pile_elem_nodes, n_pile_elements):
    """Which element END may carry a plastic hinge, one per pile node.

    A plastic hinge is a single release at a single section. Two beam elements
    meet at every interior node of a pile, so releasing each element's own end
    there would release the node twice: its rotation would then see the two
    capacities, equal and opposite, whatever it did -- no longer determined by the
    beam, and the initial-stiffness iteration drifts along that direction and
    never settles, which reads as a failed trial and drops the factor of safety.

    The first element to reach a node owns the release there; the other's end is
    left elastic, and node equilibrium puts it on the capacity with the opposite
    sign anyway. Every node the pile line ends on has one element and owns itself.
    """
    ends = np.ones((n_pile_elements, 2), dtype=bool)
    if pile_elem_nodes is None or len(pile_elem_nodes) != n_pile_elements:
        return ends
    owner = set()
    for p in range(n_pile_elements):
        for e in (0, 1):
            node = int(pile_elem_nodes[p][e])
            if node in owner:
                ends[p, e] = False
            else:
                owner.add(node)
    return ends


def _pile_element_capacity(u_local, cos_t, sin_t, L, EA, EI, K_local, n_node,
                           V_cap, M_cap, hinge_ends=None):
    """One pile beam element's actions under its structural capacities, together
    with the body-load correction that enforces them.

    Returns ``(axial, V, M1, M2, corr_local, p_rot, yielded_V, yielded_M)``, with
    ``corr_local`` in the element's LOCAL frame. The viscoplastic scheme solves
    ``K u = base_loads + corrections``, so the internal force the converged state
    is in equilibrium with is ``K u - corrections``: a correction of ``K_local p``
    leaves ``K_local (u_local - p)`` in the member, which is what makes the
    returned actions the ones the equilibrium actually carries rather than a
    clipped report of an uncapped state.

    **Moment.** A plastic hinge (``_pile_moment_hinge``), whose correction is the
    full element vector ``K_local p`` and not a moment at the shared rotational
    degree of freedom. ``hinge_ends`` says which of the element's two ends may
    carry the release, so that each pile NODE is released once (see
    ``_pile_hinge_ends``). The axial force and the shear are re-read off the
    released displacement ``u_local - p``, because a hinge changes them.

    **Shear.** The BAR's convention: the correction is ``V - V_true`` on the
    shear's own internal-force pattern (``+1`` and ``-1`` on the two ends'
    transverse rows), so ``K u - corrections`` leaves ``V_true`` in the member.
    The opposite sign -- which this code carried -- is an anti-cap, exactly as the
    comment beside the bar's correction says: the member ends up delivering
    ``2 V - V_true``, so it gets STIFFER the tighter the cap is drawn and the
    delivered shear RISES as the capacity falls. The pattern has no component on
    the rotational rows, so it leaves the capped end moments untouched.
    """
    if K_local is None:                 # a fem_data built before it was cached
        K_local = _beam_local_stiffness(EA, EI, L, n_node == 3)

    axial, V, M1, M2, _ = _pile_element_actions(
        None, cos_t, sin_t, L, EA, EI, K_local, n_node, u_local=u_local)

    n_dof = 3 * n_node
    corr_local = np.zeros(n_dof)
    p_rot = np.zeros(2)
    yielded_M = False
    yielded_V = False

    if M_cap < float('inf') and (abs(M1) > M_cap or abs(M2) > M_cap):
        p_rot, m = _pile_moment_hinge(M1, M2, M_cap, K_local, allowed=hinge_ends)
        if p_rot.any():
            p_local = np.zeros(n_dof)
            p_local[2] = p_rot[0]
            p_local[5] = p_rot[1]
            corr_local += K_local @ p_local
            M1, M2 = float(m[0]), float(m[1])
            axial, V, _, _, _ = _pile_element_actions(
                None, cos_t, sin_t, L, EA, EI, K_local, n_node,
                u_local=u_local - p_local)
        yielded_M = True

    if abs(V) > V_cap:
        V_true = np.sign(V) * V_cap
        corr_local[1] += V - V_true
        corr_local[4] -= V - V_true
        V = float(V_true)
        yielded_V = True

    return axial, V, M1, M2, corr_local, p_rot, yielded_V, yielded_M


def _beam_rotation(cos_t, sin_t, n_node):
    """Local-from-global rotation for a beam element, block diagonal with one
    [[c, s, 0], [-s, c, 0], [0, 0, 1]] per node."""
    n = 3 * n_node
    T = np.zeros((n, n))
    for k in range(n_node):
        b = 3 * k
        T[b, b], T[b, b + 1] = cos_t, sin_t
        T[b + 1, b], T[b + 1, b + 1] = -sin_t, cos_t
        T[b + 2, b + 2] = 1.0
    return T


def build_fem_data(slope_data, mesh=None, verbose=False):
    """
    Build a fem_data dictionary from slope_data and optional mesh.
    
    This function takes a slope_data dictionary (from load_slope_data) and optionally a mesh
    dictionary and constructs a fem_data dictionary suitable for finite element slope stability
    analysis using the Shear Strength Reduction Method (SSRM).
    
    The function:
    1. Extracts or loads mesh information (nodes, elements, element types, element materials)
    2. Builds material property arrays (c, phi, E, nu, gamma) from the materials table
    3. Computes pore pressure field if needed (piezo or seep options)
    4. Processes reinforcement lines into 1D truss elements with material properties
    5. Constructs boundary conditions (fixed, roller, force) based on mesh geometry
    6. Converts distributed loads to equivalent nodal forces
    
    Parameters:
        slope_data (dict): Data dictionary from load_slope_data containing:
            - materials: list of material dictionaries with c, phi, gamma, E, nu, pp_option, etc.
            - mesh: optional mesh data if mesh argument is None
            - gamma_water: unit weight of water
            - k_seismic: seismic coefficient
            - reinforcement_lines: list of reinforcement line definitions
            - distributed_loads: list of distributed load definitions
            - seepage_solution: pore pressure data if pp_option is 'seep'
            - max_depth: maximum depth for fixed boundary conditions
        mesh (dict, optional): Mesh dictionary from build_mesh_from_polygons containing:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 9) of element node indices  
            - element_types: np.ndarray (n_elements,) indicating 3, 4, 6, 8, or 9 nodes per element
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
            - elements_1d: np.ndarray (n_1d_elements, 3) of 1D element node indices
            - element_types_1d: np.ndarray (n_1d_elements,) indicating 2 or 3 nodes per 1D element  
            - element_materials_1d: np.ndarray (n_1d_elements,) of reinforcement line IDs (1-based)
        verbose (bool): if True, print a per-distributed-load assembly report
            (edges/nodes used and total force vs. expected) — a sanity check for
            debugging load application. Off by default so it doesn't spam callers
            (e.g. Studio) that build fem_data routinely.

    Returns:
        dict: fem_data dictionary with the following structure:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 9) of element node indices
            - element_types: np.ndarray (n_elements,) indicating 3 for tri3 elements, 4 for quad4 elements, etc
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
            - bc_type: np.ndarray (n_nodes,) of boundary condition flags (0=free, 1=fixed, 2=x roller, 3=y roller, 4=force)
            - bc_values: np.ndarray (n_nodes, 2) of boundary condition values (f_x, f_y for type 4)
            - c_by_mat: np.ndarray (n_materials,) of cohesion values
            - phi_by_mat: np.ndarray (n_materials,) of friction angle values (degrees)
            - E_by_mat: np.ndarray (n_materials,) of Young's modulus values
            - nu_by_mat: np.ndarray (n_materials,) of Poisson's ratio values
            - gamma_by_mat: np.ndarray (n_materials,) of unit weight values
            - u: np.ndarray (n_nodes,) of pore pressures (if applicable)
            - elements_1d: np.ndarray (n_1d_elements, 3) of 1D element node indices
            - element_types_1d: np.ndarray (n_1d_elements,) indicating 2 for linear elements and 3 for quadratic elements
            - element_materials_1d: np.ndarray (n_1d_elements,) of material IDs (1-based) corresponding to reinforcement lines
            - t_allow_by_1d_elem: np.ndarray (n_1d_elements,) of maximum tensile forces for reinforcement lines
            - t_res_by_1d_elem: np.ndarray (n_1d_elements,) of residual tensile forces for reinforcement lines
            - k_by_1d_elem: np.ndarray (n_1d_elements,) of axial stiffness values for reinforcement lines
            - cos_theta_1d: np.ndarray (n_1d_elements,) of direction cosines (x) for each 1D element
            - sin_theta_1d: np.ndarray (n_1d_elements,) of direction cosines (y) for each 1D element
            - dof_indices_1d: np.ndarray (n_1d_elements, 6) of global DOF indices using dof_offset,
              ordered [end 1, end 2, mid]; only the first n_dof_by_1d_elem entries of a row are used
            - n_dof_by_1d_elem: np.ndarray (n_1d_elements,) of 4 (two-node bar) or 6 (three-node bar)
            - K_global_1d_elems: list of np.ndarray (4, 4) or (6, 6) global stiffness matrices for each 1D element
            - dof_offset: np.ndarray (n_nodes+1,) cumulative DOF count; pile nodes get 3 DOFs, others get 2
            - is_pile_node: np.ndarray (n_nodes,) boolean, True for nodes belonging to pile elements
            - n_dof_total: int, total number of DOFs (dof_offset[n_nodes])
            - dof_indices_pile: np.ndarray (n_pile_elements, 9) of global DOF indices for beam elements,
              ordered [u, v, theta] at end 1, end 2 and mid; only the first n_dof_by_pile_elem entries are used
            - n_dof_by_pile_elem: np.ndarray (n_pile_elements,) of 6 (two-node beam) or 9 (three-node beam)
            - pile_elem_nodes: np.ndarray (n_pile_elements, 3) of [end 1, end 2, mid] node indices, mid = -1 without one
            - K_local_by_pile_elem: list of local-frame beam stiffness matrices, one per pile element
            - K_global_pile_elems: list of np.ndarray (6, 6) or (9, 9) global stiffness matrices for pile beam elements
            - EI_by_pile_elem: np.ndarray (n_pile_elements,) of flexural rigidity per unit width
            - EA_by_pile_elem: np.ndarray (n_pile_elements,) of axial rigidity per unit width
            - pile_head_nodes: np.ndarray of node indices for each pile line's top node
            - pile_head_fixed: np.ndarray of booleans, True where the pile head's
              rotation is held ('Head' = unrotated or fixed in the piles sheet)
            - pile_head_pinned: np.ndarray of booleans, True where the pile head's
              translations are held ('Head' = pinned or fixed)
            - pile_tip_nodes: np.ndarray of node indices for each pile line's bottom node
            - pile_tip_fixed: np.ndarray of booleans, True where the pile tip's
              rotation is held ('Tip' = unrotated or fixed)
            - pile_tip_pinned: np.ndarray of booleans, True where the pile tip's
              translations are held ('Tip' = pinned or fixed)
            - unit_weight: float, unit weight of water
            - k_seismic: float, seismic coefficient (horizontal acceleration / gravity)
    """
    
    # === WATER LOADS (v22 main!D23) ===
    # In automatic mode the ponded-water surface load is derived from the model's
    # own water definition and applied as tractions below, from the SAME derivation
    # the LEM slice forces use (xslope.water.with_water_loads). That is the point of
    # routing both engines through one function: distributed-load direction was
    # fixed twice and differently once already, because each engine carried its own
    # understanding of the same concept, and water is now structurally immune to
    # that. In manual mode the model comes back by identity and nothing changes.
    from .water import with_water_loads
    slope_data = with_water_loads(slope_data)

    # Get mesh data - either provided or from slope_data
    if mesh is None:
        if 'mesh' not in slope_data or slope_data['mesh'] is None:
            raise ValueError("No mesh provided and no mesh found in slope_data")
        mesh = slope_data['mesh']
    
    # Extract mesh data
    nodes = mesh["nodes"]
    elements = mesh["elements"] 
    element_types = mesh["element_types"]
    element_materials = mesh["element_materials"]
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    
    # Initialize boundary condition arrays
    bc_type = np.zeros(n_nodes, dtype=int)  # 0=free, 1=fixed, 2=x roller, 3=y roller, 4=force
    bc_values = np.zeros((n_nodes, 2))  # f_x, f_y values for type 4
    
    # Build material property arrays
    materials = slope_data["materials"]
    n_materials = len(materials)
    
    c_by_mat = np.zeros(n_materials)
    phi_by_mat = np.zeros(n_materials) 
    E_by_mat = np.zeros(n_materials)
    nu_by_mat = np.zeros(n_materials)
    gamma_by_mat = np.zeros(n_materials)
    material_names = []
    
    # Check for consistent pore pressure options
    # Material "u" key may contain "none", "piezo", "seep", or NaN-like values from empty Excel cells
    def _normalize_pp_option(val):
        if val is None or (isinstance(val, str) and val.lower() == "nan") or (isinstance(val, float) and np.isnan(val)):
            return "none"
        return str(val).lower().strip()

    pp_options = [_normalize_pp_option(mat.get("u", "none")) for mat in materials]
    unique_pp_options = set([opt for opt in pp_options if opt != "none"])
    
    if len(unique_pp_options) > 1:
        raise ValueError(f"Mixed pore pressure options not allowed: {unique_pp_options}")
    
    pp_option = list(unique_pp_options)[0] if unique_pp_options else "none"
    
    for i, material in enumerate(materials):
        strength_option = material.get("option", "mc")
        
        if strength_option == "mc":
            # Mohr-Coulomb: use c and phi directly
            c_by_mat[i] = material.get("c", 0.0)
            phi_by_mat[i] = material.get("phi", 0.0)
        elif strength_option == "cp":
            # 'cp' option: undrained strength c at the reference elevation, increasing
            # by the rate cp per unit elevation below it. Assigned per element (by
            # centroid elevation) in the loop below; store the cp rate temporarily.
            cp_rate = material.get("cp", 0.0)
            c_by_mat[i] = cp_rate  # Store cp rate temporarily (used per-element below)
            phi_by_mat[i] = 0.0     # Undrained analysis
        elif strength_option in ("pow", "hb"):
            # Curved envelopes: the power curve tau = a*(sigma'_n + d)^b + c_p
            # and generalized Hoek-Brown. Both are handled by per-Gauss-point
            # tangent linearization of the F-reduced envelope inside the
            # viscoplastic loop; the c/phi arrays only carry a seed tangent,
            # assigned per element below from the overburden estimate.
            c_by_mat[i] = 0.0
            phi_by_mat[i] = 0.0
        elif strength_option == "elastic":
            # Pure linear elastic / infinite strength (v16): the element is held
            # out of plasticity ENTIRELY via elastic_mask (below), so its c/phi
            # never enter the stress update. Carry zeros -- they are inert.
            c_by_mat[i] = 0.0
            phi_by_mat[i] = 0.0
        else:
            # A blank option is legal on rows that never carry strength; any
            # OTHER option is not implemented in the FEM, and the material's
            # c/phi columns would be zeros - silently running it as
            # zero-strength soil is the failure mode this refuses.
            if strength_option:
                raise ValueError(
                    f"Material {i+1} ({material.get('name', f'Material {i+1}')}): "
                    f"strength option '{strength_option}' is not supported by the "
                    f"FEM (supported: mc, cp, pow, hb, elastic).")
            c_by_mat[i] = material.get("c", 0.0)
            phi_by_mat[i] = material.get("phi", 0.0)

        # Require critical material properties to be explicitly specified
        if "E" not in material:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Young's modulus (E) is required but not specified")
        if "nu" not in material:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Poisson's ratio (nu) is required but not specified")
        if "gamma" not in material:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Unit weight (gamma) is required but not specified")
            
        E_by_mat[i] = material["E"]
        nu_by_mat[i] = material["nu"]
        gamma_by_mat[i] = material["gamma"]
        
        # Validate material property ranges
        if E_by_mat[i] <= 0:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Young's modulus (E) must be positive, got {E_by_mat[i]}")
        if nu_by_mat[i] < 0 or nu_by_mat[i] >= 0.5:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Poisson's ratio (nu) must be in range [0, 0.5), got {nu_by_mat[i]}")
        if gamma_by_mat[i] <= 0:
            raise ValueError(f"Material {i+1} ({material.get('name', f'Material {i+1}')}): Unit weight (gamma) must be positive, got {gamma_by_mat[i]}")
        material_names.append(material.get("name", f"Material {i+1}"))
    
    # Handle c/p strength option - compute actual cohesion per element
    c_by_elem = np.zeros(n_elements)
    phi_by_elem = np.zeros(n_elements)

    # Power-curve (option 'pow') per-element parameters. The envelope is
    # linearized to an instantaneous tangent (c_i, phi_i) at the current
    # effective normal stress inside the viscoplastic loop; here we store the
    # parameters and seed c/phi at the vertical-overburden estimate so any
    # pre-loop consumer of the arrays sees sane values.
    pow_flag_by_elem = np.zeros(n_elements, dtype=bool)
    pow_a_by_elem = np.zeros(n_elements)
    pow_b_by_elem = np.zeros(n_elements)
    pow_cp_by_elem = np.zeros(n_elements)
    pow_d_by_elem = np.zeros(n_elements)

    # Hoek-Brown (option 'hb'), same treatment. mb/s/a are derived once here
    # from GSI/mi/D and carried per element; the VP loop only inverts Balmer's
    # curve for the tangent.
    hb_flag_by_elem = np.zeros(n_elements, dtype=bool)
    hb_sci_by_elem = np.zeros(n_elements)
    hb_mb_by_elem = np.zeros(n_elements)
    hb_s_by_elem = np.zeros(n_elements)
    hb_a_by_elem = np.zeros(n_elements)
    _gs = slope_data.get('ground_surface')
    if _gs is not None and not _gs.is_empty:
        _gxy = np.asarray(_gs.coords)
        _gorder = np.argsort(_gxy[:, 0])
        _gx, _gy = _gxy[_gorder, 0], _gxy[_gorder, 1]
    else:
        _gx = _gy = None

    # Depth of every node below the ground surface, for the optional min_slip_depth
    # surficial-failure filter (see solve_fem). Positive = below ground. The ground
    # profile is single-valued in x, so a vertical interpolation is exact. With no
    # ground surface the depth is zero everywhere; solve_fem then raises if the filter
    # is switched on (rather than silently masking the whole mesh), so a ground surface
    # is required to use min_slip_depth.
    if _gx is not None:
        node_depth = np.interp(nodes[:, 0], _gx, _gy) - nodes[:, 1]
    else:
        node_depth = np.zeros(len(nodes))

    for elem_idx in range(n_elements):
        mat_id = element_materials[elem_idx] - 1  # Convert to 0-based
        material = materials[mat_id]
        strength_option = material.get("option", "mc")
        
        if strength_option == "cp":
            cp_rate = c_by_mat[mat_id]  # cp rate stored above
            c_base = material.get("c", 0.0)
            r_elev = material.get("r_elev", 0.0)

            # Compute element centroid
            elem_nodes = elements[elem_idx]
            elem_type = element_types[elem_idx]
            active_nodes = elem_nodes[:elem_type]  # Only use active nodes
            elem_coords = nodes[active_nodes]
            centroid_y = np.mean(elem_coords[:, 1])

            # Su = c + cp * max(0, r_elev - y): base strength c at r_elev, increasing
            # by the rate cp per unit elevation below it (clamped to c at/above r_elev).
            depth = max(0.0, r_elev - centroid_y)
            c_by_elem[elem_idx] = c_base + cp_rate * depth
            phi_by_elem[elem_idx] = 0.0
        elif strength_option == "pow":
            a = material.get("pow_a", 0.0)
            b = material.get("pow_b", 1.0)
            cp_ = material.get("pow_c", 0.0)
            d_ = material.get("pow_d", 0.0)
            pow_flag_by_elem[elem_idx] = True
            pow_a_by_elem[elem_idx] = a
            pow_b_by_elem[elem_idx] = b
            pow_cp_by_elem[elem_idx] = cp_
            pow_d_by_elem[elem_idx] = d_
            # seed tangent at gamma * depth-below-ground of the centroid
            elem_nodes = elements[elem_idx]
            elem_type = element_types[elem_idx]
            elem_coords = nodes[elem_nodes[:elem_type]]
            cx = float(np.mean(elem_coords[:, 0]))
            cy = float(np.mean(elem_coords[:, 1]))
            y_top = (float(np.interp(cx, _gx, _gy)) if _gx is not None
                     else float(np.max(nodes[:, 1])))
            s_n = max(gamma_by_mat[mat_id] * max(y_top - cy, 0.0), 0.0)
            s_eff = max(s_n + d_, 1e-4 * max(1.0, s_n))
            slope = a * b * s_eff ** (b - 1.0)
            tau = a * s_eff ** b + cp_
            c_by_elem[elem_idx] = tau - s_n * slope
            phi_by_elem[elem_idx] = np.degrees(np.arctan(slope))
        elif strength_option == "hb":
            sci_ = material.get("hb_sci", 0.0)
            mb_, s_, a_ = hb_constants(material.get("hb_gsi", 0.0),
                                       material.get("hb_mi", 0.0),
                                       material.get("hb_d", 0.0))
            hb_flag_by_elem[elem_idx] = True
            hb_sci_by_elem[elem_idx] = sci_
            hb_mb_by_elem[elem_idx] = mb_
            hb_s_by_elem[elem_idx] = s_
            hb_a_by_elem[elem_idx] = a_
            # seed tangent at gamma * depth-below-ground of the centroid
            elem_nodes = elements[elem_idx]
            elem_type = element_types[elem_idx]
            elem_coords = nodes[elem_nodes[:elem_type]]
            cx = float(np.mean(elem_coords[:, 0]))
            cy = float(np.mean(elem_coords[:, 1]))
            y_top = (float(np.interp(cx, _gx, _gy)) if _gx is not None
                     else float(np.max(nodes[:, 1])))
            s_n = max(gamma_by_mat[mat_id] * max(y_top - cy, 0.0), 0.0)
            c_seed, phi_seed = hb_tangent_const(s_n, sci_, mb_, s_, a_)
            c_by_elem[elem_idx] = float(c_seed)
            phi_by_elem[elem_idx] = float(phi_seed)
        else:
            c_by_elem[elem_idx] = c_by_mat[mat_id]
            phi_by_elem[elem_idx] = phi_by_mat[mat_id]
    
    # Process pore pressures
    u = np.zeros(n_nodes)
    # Signed nodal pore pressure for the opt-in matric-suction option: identical to
    # u below the water table, but NOT clamped above it, so the negative (suction)
    # part of an unsaturated seepage field survives (the effective-normal u keeps its
    # own clamp). None unless a seep field carries suction — mirrors the LEM fix in
    # fd7344b, where interpolate_at_point gained signed=True. Consumed only by the
    # suction machinery in solve_fem (gated on suction_phi_b), so it is inert by
    # default.
    u_signed = None
    sigma_v = None          # nodal vertical soil stress (ru option only)
    piezo_line_coords = None

    if pp_option == "piezo":
        # Find nodes and compute pore pressure from piezometric line
        # Look for piezometric line in various possible locations
        if "piezo_line" in slope_data:
            piezo_line_coords = slope_data["piezo_line"]
        elif "profile_lines" in slope_data:
            # Check if one of the profile lines is designated as piezo
            for line in slope_data["profile_lines"]:
                if line.get('type') == 'piezo':
                    piezo_line_coords = line['coords']
                    break
        
        if piezo_line_coords and len(piezo_line_coords) >= 2:
            gamma_water = require_gamma_water(slope_data, "FEM piezometric pore pressure")
            # u = gamma_w * VERTICAL distance below the piezometric line at
            # the node's x - the same convention the LEM slicer uses
            # (slice.get_piezometric_y_coordinates) and the hand/RS2
            # convention. Closest-point projection reads lower u under a
            # sloping water table and was inconsistent with the LEM.
            px = np.array([p[0] for p in piezo_line_coords], dtype=float)
            py = np.array([p[1] for p in piezo_line_coords], dtype=float)
            order = np.argsort(px)
            px, py = px[order], py[order]

            # Lines declared Type='phreatic' on the piezo sheet get the
            # phreatic-inclination correction (same flag and formula as the
            # LEM slicer): u = gamma_w * h_vertical * cos^2(local slope).
            _phreatic = bool(slope_data.get('piezo_phreatic', False))
            # Every node reads this line, so every node must lie inside it.
            # Extrapolating the end elevation would invent a water body the file
            # does not declare; reading zero would delete one it does. Both are
            # silent, so an out-of-extent node is an error (see
            # _check_piezo_extent).
            _check_piezo_extent(nodes[:, 0], px, py, "mesh node")
            for i, node in enumerate(nodes):
                # every node is inside the line (guarded above), so np.interp's
                # end-clamping only ever bites within floating-point noise
                piezo_elevation = float(np.interp(node[0], px, py))
                if node[1] < piezo_elevation:
                    u[i] = gamma_water * (piezo_elevation - node[1])
                    if _phreatic:
                        u[i] *= float(_piezo_cos2(node[0], px, py))
                else:
                    u[i] = 0.0
        else:
            # u='piezo' and no line to read: refused, not silently zeroed.
            raise _no_piezo_line_error(
                [m.get('name') for m in materials
                 if _normalize_pp_option(m.get("u", "none")) == "piezo"])

    elif pp_option == "seep":
        # Use existing seep solution
        if "seep_u" in slope_data:
            seep_u = slope_data["seep_u"]
            if isinstance(seep_u, np.ndarray) and len(seep_u) == n_nodes:
                u = np.maximum(0.0, seep_u)
                # Keep the raw signed field for the suction option (see u_signed
                # above). max(0, signed) == u, so the effective normal is unchanged.
                u_signed = np.asarray(seep_u, dtype=float)
            else:
                n_seep = len(seep_u) if hasattr(seep_u, "__len__") else "?"
                warnings.warn(
                    f"Seepage pore pressures NOT applied: the stored seep solution has {n_seep} "
                    f"values but the FEM mesh has {n_nodes} nodes. Pore pressures default to 0, "
                    "which over-predicts the factor of safety — run the FEM on the same mesh as "
                    "the seepage solution.")

    elif pp_option == "ru":
        # Pore-pressure ratio (template v12): u = ru * sigma_v, where sigma_v
        # is the vertical total stress of the SOIL column above the point —
        # the same definition the LEM slicer uses (slice.py mat_u == 'ru':
        # u = ru * sum(gamma_i * h_i); distributed loads and crack water are
        # excluded by definition, Bishop & Morgenstern). The overburden is
        # integrated by intersecting a vertical ray from each node with the
        # material polygons, which handles multi-band zones exactly; moist
        # gamma throughout, matching the LEM's no-water-table path (ru models
        # carry no piezometric surface). ru itself is PER MATERIAL, so the
        # shared nodal u array (ambiguous on material boundaries) stays zero
        # here; the Gauss-point precompute applies the element material's ru
        # to sigma_v interpolated from the nodes.
        from shapely.geometry import LineString as _RayLS, Polygon as _RayPoly
        from .mesh import get_material_polygons as _gmp
        _ray_polys = []
        for _pi, _pd in enumerate(_gmp(slope_data)):
            _mid = _pd.get('mat_id')
            _midx = _mid if (_mid is not None and 0 <= _mid < len(materials)) else _pi
            _ray_polys.append((_RayPoly(_pd['coords']),
                               float(materials[_midx].get('gamma', 0.0))))
        _y_top = float(np.max(nodes[:, 1])) + 1.0
        sigma_v = np.zeros(n_nodes)
        for i, node in enumerate(nodes):
            _x0, _y0 = float(node[0]), float(node[1])
            _ray = _RayLS([(_x0, _y0), (_x0, _y_top)])
            _sv = 0.0
            for _poly, _g in _ray_polys:
                _minx, _miny, _maxx, _maxy = _poly.bounds
                if _x0 < _minx or _x0 > _maxx or _y0 >= _maxy:
                    continue
                _inter = _ray.intersection(_poly)
                if not _inter.is_empty:
                    _sv += _g * _inter.length
            sigma_v[i] = _sv
    
    # Process 1D reinforcement elements
    elements_1d = np.array([]).reshape(0, 3) if 'elements_1d' not in mesh else mesh['elements_1d']
    element_types_1d = np.array([]) if 'element_types_1d' not in mesh else mesh['element_types_1d'] 
    element_materials_1d = np.array([]) if 'element_materials_1d' not in mesh else mesh['element_materials_1d']
    
    n_1d_elements = len(elements_1d)
    
    t_allow_by_1d_elem = np.zeros(n_1d_elements)
    # NaN = "no post-peak drop" (see the t_res handling below). Zero would mean
    # brittle rupture, which must not be the default for an unset field.
    t_res_by_1d_elem = np.full(n_1d_elements, np.nan)
    k_by_1d_elem = np.zeros(n_1d_elements)
    cos_theta_1d = np.zeros(n_1d_elements)
    sin_theta_1d = np.zeros(n_1d_elements)
    # Six DOF slots per bar: (x, y) at end 1, end 2 and -- on a quadratic mesh --
    # the midside node of the soil edge the bar lies on. n_dof_by_1d_elem says how
    # many of them an element actually uses, 4 or 6, and every reader slices by it.
    dof_indices_1d = np.zeros((n_1d_elements, 6), dtype=int)
    n_dof_by_1d_elem = np.full(n_1d_elements, 4, dtype=int)
    K_global_1d_elems = []
    # Per-1D-element geometry along the reinforcement line, used by the 1D details
    # panel to order a line's elements and place them on its arc length. Zero for
    # pile elements and for elements this loop skips.
    elem_length_1d = np.zeros(n_1d_elements)
    dist_end1_1d = np.zeros(n_1d_elements)   # centroid -> line end 1
    dist_end2_1d = np.zeros(n_1d_elements)   # centroid -> line end 2

    if n_1d_elements > 0 and "reinforcement_lines" in slope_data:
        from .fileio import ensure_reinforce_pullout
        ensure_reinforce_pullout(slope_data)
        reinforcement_lines = slope_data["reinforcement_lines"]

        for elem_idx in range(n_1d_elements):
            line_id = element_materials_1d[elem_idx] - 1  # Convert to 0-based

            if line_id < len(reinforcement_lines):
                line_data = reinforcement_lines[line_id]

                # Element geometry. The chord runs between end nodes [0] and [1];
                # node [2], where the mesh recorded one, is the midside node of
                # the soil edge the bar lies on and is a node of the bar too.
                elem_nodes_1d = elements_1d[elem_idx]
                node_0 = elem_nodes_1d[0]
                node_1 = elem_nodes_1d[1]
                node_m = _midside_1d_node(elements_1d, element_types_1d, elem_idx)
                coord_0 = nodes[node_0]
                coord_1 = nodes[node_1]

                # Compute element length from end nodes
                dx = coord_1[0] - coord_0[0]
                dy = coord_1[1] - coord_0[1]
                elem_length = np.sqrt(dx * dx + dy * dy)
                elem_centroid = 0.5 * (coord_0 + coord_1)

                if elem_length > 1e-12:
                    # Direction cosines
                    cos_t = dx / elem_length
                    sin_t = dy / elem_length
                    cos_theta_1d[elem_idx] = cos_t
                    sin_theta_1d[elem_idx] = sin_t

                    # DOF indices, in the element's own node order
                    # [end 1, end 2, mid]. Rebuilt against dof_offset below.
                    if node_m is None:
                        dof_indices_1d[elem_idx, :4] = [2*node_0, 2*node_0+1,
                                                        2*node_1, 2*node_1+1]
                        n_dof_by_1d_elem[elem_idx] = 4
                    else:
                        dof_indices_1d[elem_idx] = [2*node_0, 2*node_0+1,
                                                    2*node_1, 2*node_1+1,
                                                    2*node_m, 2*node_m+1]
                        n_dof_by_1d_elem[elem_idx] = 6

                    # Compute distance from element centroid to line ends
                    x1, y1 = line_data.get("x1", 0), line_data.get("y1", 0)
                    x2, y2 = line_data.get("x2", 0), line_data.get("y2", 0)

                    dist_to_left = np.linalg.norm(elem_centroid - [x1, y1])
                    dist_to_right = np.linalg.norm(elem_centroid - [x2, y2])

                    # Retain the line geometry for the 1D details panel.
                    elem_length_1d[elem_idx] = elem_length
                    dist_end1_1d[elem_idx] = dist_to_left
                    dist_end2_1d[elem_idx] = dist_to_right

                    # Get reinforcement properties
                    t_max = line_data.get("t_max", 0.0)
                    # NaN = unset: no post-peak drop (elastic-perfectly-plastic).
                    # An explicit 0.0 is different — it means brittle rupture.
                    t_res = line_data.get("t_res", float('nan'))
                    if t_res is None:
                        t_res = float('nan')
                    lp1 = line_data.get("lp1", 0.0)   # Pullout length left end
                    lp2 = line_data.get("lp2", 0.0)   # Pullout length right end
                    tend1 = line_data.get("tend1", 0.0)  # End anchorage capacities
                    tend2 = line_data.get("tend2", 0.0)

                    # Allowable tension from the capacity envelope shared with the
                    # LEM point list (fileio.reinforce_available_tension): tensile
                    # strength, frictional development from BOTH ends, and end
                    # anchorage. With tend = 0 this reproduces the historical
                    # nearest-end taper for any element whose centroid lies in at
                    # most one pullout zone, and is the correct min() when zones
                    # overlap. A line on the overburden-dependent law carries a
                    # pullout profile instead of a development length, and the
                    # same call reads that: both engines take their capacity
                    # from one function under either law.
                    from .fileio import reinforce_available_tension
                    t_allow = reinforce_available_tension(
                        dist_to_left, dist_to_right, t_max, lp1, lp2, tend1, tend2,
                        pullout=line_data.get("_pullout_profile"))
                    t_allow_by_1d_elem[elem_idx] = t_allow

                    if t_res != t_res:
                        # Unset: this line never softens, anywhere along its length.
                        t_res_by_1d_elem[elem_idx] = float('nan')
                    else:
                        # Two independent post-peak mechanisms, and the element is
                        # governed by whichever leaves it the less capacity.
                        #
                        # Bond slip is PERFECTLY PLASTIC: an element inside a
                        # pullout ramp that reaches the capacity its embedment can
                        # develop keeps slipping AT that capacity. The soil-bar
                        # interface is frictional, and friction does not vanish
                        # once it is overcome. t_allow — the ramped envelope,
                        # including end anchorage — is therefore both the yield
                        # force and the post-slip force.
                        #
                        # Tres is the RUPTURE residual: what the bar itself retains
                        # once the material tears. It is a property of the
                        # reinforcement, not of the embedment.
                        #
                        # So the element can never carry more than the smaller of
                        # the two. This is the same envelope the limit-equilibrium
                        # engine applies (fileio.reinforce_available_tension) and
                        # the standard cable/geogrid treatment.
                        #
                        # (Elements inside a ramp were dropped all the way to zero
                        # in earlier versions — "sudden complete pullout" — so
                        # results predating this rule read lower.)
                        t_res_by_1d_elem[elem_idx] = min(t_res, t_allow)

                    # Compute axial stiffness. E and Area are OPTIONAL in the
                    # input (the LEM needs neither — it applies the capacity
                    # envelope directly), so a perfectly valid LEM file can
                    # arrive here with them blank, which load_slope_data turns
                    # into NaN. Left alone, that poisons k, then K_global, and
                    # the solve dies with an opaque "Factor is exactly singular".
                    # Fail here instead, naming the line and what to supply.
                    E = line_data.get("E")
                    A = line_data.get("area")
                    if (E is None or A is None
                            or not np.isfinite(E) or not np.isfinite(A)
                            or E <= 0 or A <= 0):
                        label = line_data.get("label") or f"line {line_id + 1}"
                        raise ValueError(
                            f"Reinforcement '{label}' has no usable axial stiffness "
                            f"(E={E}, Area={A}). The FEM models reinforcement as a bar "
                            f"element, so it needs E and Area (the axial rigidity "
                            f"EA = E*Area per unit width) on the 'reinforce' sheet. "
                            f"The LEM does not — it applies the tensile capacity "
                            f"envelope (Tmax/Lp) directly — so this file can run in the "
                            f"LEM but not the FEM until E and Area are filled in.")
                    # The CHORD stiffness EA/L. It is the axial stiffness of the
                    # two-node bar, and it is also what force recovery multiplies
                    # the chord elongation by on the three-node bar, where
                    # EA(u2-u1)/L is the axial force at the element center.
                    k_val = E * A / elem_length
                    k_by_1d_elem[elem_idx] = k_val

                    if node_m is None:
                        # Two-node bar, DOFs [x1, y1, x2, y2].
                        # T = [[cos, sin, 0, 0], [0, 0, cos, sin]]  (2x4)
                        # K_axial = k * [[1, -1], [-1, 1]]          (2x2)
                        # K_global_elem = T^T @ K_axial @ T         (4x4)
                        c2 = cos_t * cos_t
                        cs = cos_t * sin_t
                        s2 = sin_t * sin_t
                        K_global_elem = k_val * np.array([
                            [ c2,  cs, -c2, -cs],
                            [ cs,  s2, -cs, -s2],
                            [-c2, -cs,  c2,  cs],
                            [-cs, -s2,  cs,  s2]
                        ])
                    else:
                        # Three-node isoparametric quadratic bar, axial DOFs in
                        # the element's node order [end 1, end 2, mid]:
                        #
                        #   K_axial = EA/(3L) * [[ 7,  1, -8],
                        #                        [ 1,  7, -8],
                        #                        [-8, -8, 16]]
                        #
                        # (the exact two-point Gauss integral of EA (dN/dx)^2).
                        # Rotated into global coordinates by the 3x6
                        # T = [[c, s, 0, 0, 0, 0],
                        #      [0, 0, c, s, 0, 0],
                        #      [0, 0, 0, 0, c, s]].
                        K_axial = (E * A / (3.0 * elem_length)) * np.array([
                            [ 7.0,  1.0, -8.0],
                            [ 1.0,  7.0, -8.0],
                            [-8.0, -8.0, 16.0],
                        ])
                        T_bar = np.zeros((3, 6))
                        T_bar[0, 0], T_bar[0, 1] = cos_t, sin_t
                        T_bar[1, 2], T_bar[1, 3] = cos_t, sin_t
                        T_bar[2, 4], T_bar[2, 5] = cos_t, sin_t
                        K_global_elem = T_bar.T @ K_axial @ T_bar
                    K_global_1d_elems.append(K_global_elem)
                else:
                    K_global_1d_elems.append(np.zeros((4, 4)))
            else:
                K_global_1d_elems.append(np.zeros((4, 4)))

    # Reinforcement line labels, in line order — how the 1D details panel names a
    # line when the model gives it one.
    reinforce_line_labels = [ln.get("label") for ln in
                             slope_data.get("reinforcement_lines", [])]

    # === PILE BEAM ELEMENTS ===
    # Pile 1D elements are identified by element_materials_1d values that exceed
    # the number of reinforcement lines (since constraint lines are ordered:
    # reinforcement first, then piles).
    n_reinf_lines = len(slope_data.get("reinforcement_lines", []))
    pile_lines = slope_data.get("pile_lines", [])
    n_pile_lines = len(pile_lines)

    # Identify which 1D elements are pile elements
    pile_elem_mask = np.zeros(n_1d_elements, dtype=bool)
    n_pile_elements = 0
    cos_theta_pile = []
    sin_theta_pile = []
    K_global_pile_elems = []
    K_local_by_pile_elem = []   # local-frame stiffness, for end-action recovery
    pile_elem_indices = []  # maps pile element index to global 1D element index
    V_cap_by_pile_elem = []
    M_cap_by_pile_elem = []
    elem_length_by_pile_elem = []
    S_by_pile_elem = []
    EI_by_pile_elem = []
    EA_by_pile_elem = []
    pile_node_pairs = []  # (node_0, node_1) -- the END nodes of each pile element
    # Every node of each pile element, [end 1, end 2, mid], with mid = -1 where
    # the element has none. pile_node_pairs stays the end pair, so readers that
    # chain elements end to end are unaffected.
    pile_elem_nodes = []
    n_dof_by_pile_elem = []  # 6 on a two-node element, 9 on a three-node one
    pile_line_idx_by_pile_elem = []  # which pile_line each pile element belongs to

    if n_1d_elements > 0 and n_pile_lines > 0:
        for elem_idx in range(n_1d_elements):
            line_id = element_materials_1d[elem_idx] - 1  # 0-based
            pile_line_idx = line_id - n_reinf_lines  # index into pile_lines

            if pile_line_idx < 0 or pile_line_idx >= n_pile_lines:
                continue

            pile_data = pile_lines[pile_line_idx]
            pile_elem_mask[elem_idx] = True
            pile_elem_indices.append(elem_idx)

            # Get element geometry
            elem_nodes = elements_1d[elem_idx]
            node_0 = elem_nodes[0]
            node_1 = elem_nodes[1]
            node_m = _midside_1d_node(elements_1d, element_types_1d, elem_idx)
            coord_0 = nodes[node_0]
            coord_1 = nodes[node_1]

            dx = coord_1[0] - coord_0[0]
            dy = coord_1[1] - coord_0[1]
            elem_length = np.sqrt(dx * dx + dy * dy)

            if elem_length > 1e-12:
                cos_t = dx / elem_length
                sin_t = dy / elem_length

                # Get pile properties
                E_pile = pile_data.get("E", 0.0)
                D_pile = pile_data.get("D_pile")
                S_pile = pile_data.get("S", 1.0)  # default 1.0 = no spacing reduction
                I_pile = pile_data.get("I")
                A_pile = pile_data.get("area")

                # Auto-compute I and A from D if not provided
                if D_pile is not None:
                    if A_pile is None:
                        A_pile = np.pi * D_pile**2 / 4.0
                    if I_pile is None:
                        I_pile = np.pi * D_pile**4 / 64.0
                else:
                    if A_pile is None:
                        A_pile = 0.0
                    if I_pile is None:
                        I_pile = 0.0

                if E_pile is None or E_pile == 0:
                    continue

                # Scale by 1/S for per-unit-width (2D plane strain)
                if S_pile and S_pile > 0:
                    EA = E_pile * A_pile / S_pile
                    EI = E_pile * I_pile / S_pile
                else:
                    EA = E_pile * A_pile
                    EI = E_pile * I_pile

                L = elem_length

                # Local beam stiffness. Two-node: 6x6 Euler-Bernoulli in
                # [u1, v1, theta1, u2, v2, theta2]. Three-node -- on a quadratic
                # mesh, where the element also stands on the midside node of the
                # soil edge -- 9x9 in [u1, v1, th1, u2, v2, th2, u3, v3, th3],
                # quintic Hermite bending with quadratic axial. Rotated into
                # global coordinates by the matching block-diagonal T.
                three_node = node_m is not None
                n_node = 3 if three_node else 2
                K_local = _beam_local_stiffness(EA, EI, L, three_node)
                T = _beam_rotation(cos_t, sin_t, n_node)
                K_beam = T.T @ K_local @ T

                cos_theta_pile.append(cos_t)
                sin_theta_pile.append(sin_t)
                pile_node_pairs.append((node_0, node_1))
                pile_elem_nodes.append((int(node_0), int(node_1),
                                        int(node_m) if three_node else -1))
                n_dof_by_pile_elem.append(3 * n_node)
                K_global_pile_elems.append(K_beam)
                K_local_by_pile_elem.append(K_local)
                elem_length_by_pile_elem.append(elem_length)
                EI_by_pile_elem.append(EI)
                EA_by_pile_elem.append(EA)
                pile_line_idx_by_pile_elem.append(pile_line_idx)

                # Structural capacity (per-unit-width = per-pile / S)
                V_cap_pile = pile_data.get("V_cap")
                M_cap_pile = pile_data.get("M_cap")
                V_cap_by_pile_elem.append(V_cap_pile / S_pile if V_cap_pile is not None else float('inf'))
                M_cap_by_pile_elem.append(M_cap_pile / S_pile if M_cap_pile is not None else float('inf'))
                S_by_pile_elem.append(S_pile)

                n_pile_elements += 1

    cos_theta_pile = np.array(cos_theta_pile)
    sin_theta_pile = np.array(sin_theta_pile)
    pile_elem_nodes = (np.array(pile_elem_nodes, dtype=int) if n_pile_elements > 0
                       else np.zeros((0, 3), dtype=int))
    n_dof_by_pile_elem = (np.array(n_dof_by_pile_elem, dtype=int)
                          if n_pile_elements > 0 else np.zeros(0, dtype=int))
    pile_elem_indices = np.array(pile_elem_indices, dtype=int)
    V_cap_by_pile_elem = np.array(V_cap_by_pile_elem)
    M_cap_by_pile_elem = np.array(M_cap_by_pile_elem)
    elem_length_by_pile_elem = np.array(elem_length_by_pile_elem)
    S_by_pile_elem = np.array(S_by_pile_elem)
    EI_by_pile_elem = np.array(EI_by_pile_elem)
    EA_by_pile_elem = np.array(EA_by_pile_elem)
    pile_line_idx_by_pile_elem = np.array(pile_line_idx_by_pile_elem, dtype=int) if n_pile_elements > 0 else np.array([], dtype=int)

    # === BUILD DOF OFFSET MAP ===
    # Pile nodes get 3 DOFs (ux, uy, theta), all other nodes get 2 DOFs (ux, uy).
    # Every node a pile element stands on, midside nodes included: the midside
    # node carries the beam's deflection and slope like any other node on it.
    is_pile_node = np.zeros(n_nodes, dtype=bool)
    for p_idx in range(n_pile_elements):
        for nd in pile_elem_nodes[p_idx]:
            if nd >= 0:
                is_pile_node[nd] = True

    dof_offset = np.zeros(n_nodes + 1, dtype=int)
    for i in range(n_nodes):
        dof_offset[i + 1] = dof_offset[i] + (3 if is_pile_node[i] else 2)
    n_dof_total = int(dof_offset[n_nodes])

    # Pile element DOF indices against dof_offset: three per node, and nine slots
    # per element so a three-node element fits. n_dof_by_pile_elem says how many
    # of them each element uses, and every reader slices by it.
    dof_indices_pile = np.zeros((n_pile_elements, 9), dtype=int) if n_pile_elements > 0 else np.zeros((0, 9), dtype=int)
    for p_idx in range(n_pile_elements):
        row = []
        for nd in pile_elem_nodes[p_idx]:
            if nd < 0:
                continue
            row += [dof_offset[nd], dof_offset[nd] + 1, dof_offset[nd] + 2]
        dof_indices_pile[p_idx, :len(row)] = row

    # Rebuild 1D truss DOF indices using dof_offset (in case any share nodes with piles)
    for elem_idx in range(n_1d_elements):
        if pile_elem_mask[elem_idx]:
            continue
        elem_nodes_1d = elements_1d[elem_idx]
        node_0 = elem_nodes_1d[0]
        node_1 = elem_nodes_1d[1]
        node_m = _midside_1d_node(elements_1d, element_types_1d, elem_idx)
        row = [dof_offset[node_0], dof_offset[node_0] + 1,
               dof_offset[node_1], dof_offset[node_1] + 1]
        if node_m is not None and n_dof_by_1d_elem[elem_idx] == 6:
            row += [dof_offset[node_m], dof_offset[node_m] + 1]
        dof_indices_1d[elem_idx, :len(row)] = row

    # Identify the pile end nodes and their rotation restraints for boundary
    # conditions. The head is the top node (highest y) of each pile line and the
    # tip is the bottom node (lowest y). Each end carries its own restraint --
    # 'Head' and 'Tip' in the piles sheet -- and 'fixed' constrains that node's
    # ROTATION degree of freedom only; the translations stay with the boundary
    # conditions and the surrounding soil. ``fixity`` is read as the head for a
    # slope_data dict built before the two ends were separated.
    pile_head_nodes = []
    pile_head_fixed = []     # rotation held ('Head' = unrotated or fixed)
    pile_head_pinned = []    # translations held ('Head' = pinned or fixed)
    pile_tip_nodes = []
    pile_tip_fixed = []      # rotation restrained ('Tip' = fixed)
    pile_tip_pinned = []     # translations restrained ('Tip' = pinned or fixed)
    for pl_idx in range(n_pile_lines):
        pile_data = pile_lines[pl_idx]
        head_fixity = pile_data.get("head_fixity", pile_data.get("fixity", "free"))
        tip_fixity = pile_data.get("tip_fixity", "free")

        # Collect all nodes belonging to this pile line, midside nodes included,
        # and the end nodes separately.
        pile_nodes_for_line = set()
        end_nodes_for_line = set()
        for p_idx in range(n_pile_elements):
            if pile_line_idx_by_pile_elem[p_idx] == pl_idx:
                n0, n1 = pile_node_pairs[p_idx]
                end_nodes_for_line.update((int(n0), int(n1)))
                for nd in pile_elem_nodes[p_idx]:
                    if nd >= 0:
                        pile_nodes_for_line.add(int(nd))

        if pile_nodes_for_line:
            # Top node = highest y coordinate; bottom node = lowest y
            top_node = max(pile_nodes_for_line, key=lambda nd: nodes[nd, 1])
            bottom_node = min(pile_nodes_for_line, key=lambda nd: nodes[nd, 1])
            # A midside node sits at the midpoint of its element, strictly
            # between the two ends in y, so the extremes are end nodes. The head
            # and tip restraints are stated for the ends of the pile, and a
            # midside node picked up here would put them mid-element.
            if top_node not in end_nodes_for_line or bottom_node not in end_nodes_for_line:
                raise ValueError(
                    f"Pile line {pl_idx + 1}: the highest or lowest node of the "
                    f"pile is not an end node of one of its elements, so the head "
                    f"and tip restraints cannot be placed. Check that the pile "
                    f"line is not horizontal.")
            pile_head_nodes.append(top_node)
            pile_head_fixed.append(head_fixity in ("unrotated", "fixed"))
            pile_head_pinned.append(head_fixity in ("pinned", "fixed"))
            pile_tip_nodes.append(bottom_node)
            pile_tip_fixed.append(tip_fixity in ("unrotated", "fixed"))
            pile_tip_pinned.append(tip_fixity in ("pinned", "fixed"))

    pile_head_nodes = np.array(pile_head_nodes, dtype=int)
    pile_head_fixed = np.array(pile_head_fixed, dtype=bool)
    pile_head_pinned = np.array(pile_head_pinned, dtype=bool)
    pile_tip_nodes = np.array(pile_tip_nodes, dtype=int)
    pile_tip_fixed = np.array(pile_tip_fixed, dtype=bool)
    pile_tip_pinned = np.array(pile_tip_pinned, dtype=bool)

    # Set up boundary conditions
    
    # Step 1: Default to free (type 0)
    # Already initialized to zeros
    
    # Step 2: Fixed supports along the BOTTOM boundary (type 1) - standard
    # practice. The bottom is the domain polygon's lower boundary POLYLINE,
    # not simply y == min(y): an undulating bedrock base (max_depth absent,
    # lowest profile line forms the bottom - e.g. vp027) must be fixed along
    # its whole length, or the body is restrained at one low corner and
    # never reaches equilibrium at any strength-reduction factor. The bottom
    # polyline is every exterior segment of the domain polygon that is
    # neither on the ground surface nor a vertical side edge at the domain's
    # x-extremes; for a flat-bottomed domain this reproduces the old
    # y == y_min rule node-for-node.
    tolerance = 1e-6
    y_min = float(np.min(nodes[:, 1])) if len(nodes) > 0 else 0.0
    bottom_nodes = np.abs(nodes[:, 1] - y_min) < tolerance   # fallback rule
    _domain = slope_data.get('domain_polygon')
    _ground = slope_data.get('ground_surface')
    if (_domain is not None and _ground is not None and not _ground.is_empty
            and len(nodes) > 0):
        _ring = list(_domain.exterior.coords)
        _rx = [c[0] for c in _ring]
        _dx_min, _dx_max = min(_rx), max(_rx)
        _span = max(_dx_max - _dx_min, float(np.max(nodes[:, 1]) - y_min), 1.0)
        _geom_tol = 1e-6 * _span
        _bottom_segs = []
        for _a, _b in zip(_ring[:-1], _ring[1:]):
            _mid = Point((_a[0] + _b[0]) / 2.0, (_a[1] + _b[1]) / 2.0)
            if _ground.distance(_mid) < _geom_tol:
                continue                       # ground-surface segment
            if (abs(_a[0] - _b[0]) < _geom_tol and
                    (abs(_a[0] - _dx_min) < _geom_tol or
                     abs(_a[0] - _dx_max) < _geom_tol)):
                continue                       # vertical side edge
            _bottom_segs.append(LineString([_a, _b]))
        if _bottom_segs:
            from shapely.ops import unary_union
            _bottom_geom = unary_union(_bottom_segs)
            bottom_nodes = np.array(
                [_bottom_geom.distance(Point(nd[0], nd[1])) < _geom_tol
                 for nd in nodes], dtype=bool)
    bc_type[bottom_nodes] = 1  # Fixed (u=0, v=0)
    
    # Step 3: X-roller supports at left and right sides (type 2) - standard practice
    # The side is the domain polygon's extreme-x boundary EDGE, not simply
    # x == x_min: a truncation boundary digitized slightly off vertical (rs2_28's
    # left face drifts 0.13 over 30.7 of height) otherwise gets a roller at its
    # single extreme node and is left traction-free over the rest, which cannot
    # equilibrate at any strength-reduction factor once a pore-pressure field
    # prescribes uplift on it. A ring edge counts as a side when both its ends
    # sit within _SIDE_EDGE_DRIFT_FRAC of the domain span of that extreme and it
    # runs more vertically than horizontally; an exactly vertical side reproduces
    # the x == x_extreme rule node-for-node, and the rule is a union with it, so
    # a side node is never lost.
    if len(nodes) > 0:
        x_min = float(np.min(nodes[:, 0]))
        x_max = float(np.max(nodes[:, 0]))
        left_nodes = np.abs(nodes[:, 0] - x_min) < tolerance
        right_nodes = np.abs(nodes[:, 0] - x_max) < tolerance
        if _domain is not None:
            _ring = list(_domain.exterior.coords)
            _span = max(x_max - x_min, float(np.max(nodes[:, 1])) - y_min, 1.0)
            _geom_tol = 1e-6 * _span
            _drift = _SIDE_EDGE_DRIFT_FRAC * _span
            for _x_ext, _side in ((x_min, 'left'), (x_max, 'right')):
                _segs = [LineString([_a, _b]) for _a, _b in zip(_ring[:-1], _ring[1:])
                         if abs(_a[0] - _x_ext) < _drift and abs(_b[0] - _x_ext) < _drift
                         and abs(_b[1] - _a[1]) > abs(_b[0] - _a[0])]
                if not _segs:
                    continue
                _side_geom = unary_union(_segs)
                _on_side = np.array(
                    [_side_geom.distance(Point(nd[0], nd[1])) < _geom_tol
                     for nd in nodes], dtype=bool)
                if _side == 'left':
                    left_nodes = left_nodes | _on_side
                else:
                    right_nodes = right_nodes | _on_side

        # v21 main!D22 chooses what the side restraint IS. 'rollers' (the default, and
        # every file that does not declare it) is the historical hardwired behavior:
        # u = 0, v free, so the truncated ground can still settle under its own weight.
        # 'fixed' clamps both components, which is what RS2 does on its side
        # boundaries — a vendor-parity option, not a better model: fixing the sides
        # adds shear restraint the real ground does not have, and stiffens a domain
        # truncated close to the slope.
        _side_bc = str(slope_data.get('side_bc') or 'rollers').strip().lower()
        if _side_bc not in ('rollers', 'fixed'):
            raise ValueError(
                f"Unknown side boundary condition {_side_bc!r}; expected 'rollers' "
                "or 'fixed'.")
        _side_code = 1 if _side_bc == 'fixed' else 2

        # Apply the side restraint but preserve existing boundary conditions (the
        # fully-fixed bottom takes precedence at corners)
        left_not_fixed = left_nodes & (bc_type != 1)
        right_not_fixed = right_nodes & (bc_type != 1)

        bc_type[left_not_fixed] = _side_code
        bc_type[right_not_fixed] = _side_code
    
    # Save displacement constraints before force BCs can overwrite them.
    # Nodes on boundary faces that also receive distributed loads need both
    # their displacement constraint (roller/fixed) AND the applied force.
    fixed_nodes = set(np.where(bc_type == 1)[0])
    roller_x_nodes = set(np.where(bc_type == 2)[0])
    roller_y_nodes = set(np.where(bc_type == 3)[0])

    # Step 4: Convert distributed loads to nodal forces (type 4)
    # Check for distributed loads (could be 'dloads', 'dloads2', or 'distributed_loads')
    # Each entry is (load_line, direction). The direction is the v21 per-block
    # 'normal' | 'vertical' (blank/pre-v21 = 'normal'), carried alongside its own
    # list so a model may mix the two — see fileio's dload_dirs / dload2_dirs.
    distributed_loads = []

    def _with_dirs(key, dirs_key):
        lines = slope_data.get(key) or []
        dirs = slope_data.get(dirs_key) or []
        return [(ln, str(dirs[i] if i < len(dirs) else 'normal').lower())
                for i, ln in enumerate(lines)]

    distributed_loads += _with_dirs("dloads", "dload_dirs")
    distributed_loads += _with_dirs("dloads2", "dload2_dirs")
    distributed_loads += _with_dirs("distributed_loads", "distributed_load_dirs")
    # The derived water loads (empty in manual mode), each stage alongside the user
    # sheet it accompanies and in the same order, so an automatic model assembles
    # exactly the tractions the equivalent manual model does. They carry no dirs
    # list: water acts normal to the surface, which is what a missing entry means.
    distributed_loads += _with_dirs("dloads_derived", "")
    distributed_loads += _with_dirs("dloads2_derived", "")

    if distributed_loads:
        tolerance = 1e-1  # Tolerance for finding nodes on load lines

        for load_idx, (load_line, load_dir) in enumerate(distributed_loads):
            # Handle different possible data structures
            if isinstance(load_line, dict) and "coords" in load_line:
                load_coords = load_line["coords"]
                load_values = load_line["loads"]
            elif isinstance(load_line, list):
                load_coords = [(pt["X"], pt["Y"]) for pt in load_line]
                load_values = [pt["Normal"] for pt in load_line]
            else:
                continue

            if len(load_coords) < 2 or len(load_values) < 2:
                continue

            load_linestring = LineString(load_coords)
            load_total_length = load_linestring.length

            # Pass 1: Collect all nodes on the load line with their projected distances
            load_nodes = []  # list of (node_index, projected_distance)
            for i, node in enumerate(nodes):
                node_point = Point(node)
                if load_linestring.distance(node_point) <= tolerance:
                    proj_dist = load_linestring.project(node_point)
                    load_nodes.append((i, proj_dist))

            if not load_nodes:
                continue

            # Sort by projected distance along the load line
            load_nodes.sort(key=lambda x: x[1])

            # Interpolate the traction magnitude at a projected distance
            def _p_at(proj_dist):
                cumulative_length = 0
                val = load_values[-1]
                for j in range(len(load_coords) - 1):
                    seg_length = np.linalg.norm(np.array(load_coords[j+1]) - np.array(load_coords[j]))
                    cumulative_length += seg_length
                    if proj_dist <= cumulative_length:
                        local_distance = proj_dist - (cumulative_length - seg_length)
                        ratio = local_distance / seg_length if seg_length > 0 else 0
                        val = load_values[j] * (1 - ratio) + load_values[j+1] * ratio
                        break
                return val

            # Pass 2a: CONSISTENT edge-load integration. Tributary lumping is
            # wrong for quadratic edges (consistent distribution is 1/6-2/3-1/6
            # corner-mid-corner, not 1/4-1/2-1/4): the difference is a chain of
            # self-equilibrated nodal force couples of magnitude ~p*L/6 along
            # the loaded boundary, which shows up as spurious near-surface
            # stress oscillations of order p/6. For submerged faces where the
            # traction p is large compared to the soil strength (reservoir
            # loading), those oscillations are strong enough to yield the skin
            # elements and prevent convergence. Integrating N_i * p over each
            # boundary edge eliminates the couples exactly.
            node_proj = dict(load_nodes)
            on_line = set(node_proj)
            _seen_edges = {}
            _load_edges = []
            for _eidx, _elem in enumerate(elements):
                _et = element_types[_eidx]
                if _et == 3:
                    _ledges = [(0, 1), (1, 2), (2, 0)]
                elif _et == 6:
                    _ledges = [(0, 3, 1), (1, 4, 2), (2, 5, 0)]
                elif _et == 4:
                    _ledges = [(0, 1), (1, 2), (2, 3), (3, 0)]
                elif _et in (8, 9):
                    _ledges = [(0, 4, 1), (1, 5, 2), (2, 6, 3), (3, 7, 0)]
                else:
                    continue
                for _le in _ledges:
                    _gn = [int(_elem[i]) for i in _le]
                    if all(n in on_line for n in _gn):
                        _key = tuple(sorted((_gn[0], _gn[-1])))
                        if _key in _seen_edges:
                            # A SECOND element owns this edge, so it is interior
                            # to the mesh and "into the material" has no meaning.
                            # Drop the centroid recorded for the first owner and
                            # let the tangent rule below decide, unchanged.
                            _first = _seen_edges[_key]
                            _load_edges[_first] = (_load_edges[_first][0], None)
                            continue
                        _ncorner = 4 if _et in (4, 8, 9) else 3
                        _cen = nodes[[int(_elem[i])
                                      for i in range(_ncorner)]].mean(axis=0)
                        _seen_edges[_key] = len(_load_edges)
                        _load_edges.append((_gn, _cen))

            if _load_edges:
                nodal_fx = {}
                nodal_fy = {}
                _g = 1.0 / np.sqrt(3.0)
                for _gn, _cen in _load_edges:
                    c1, c2 = _gn[0], _gn[-1]
                    # orient the edge along increasing projection so the
                    # inward normal (tangent rotated 90 degrees CW) points
                    # into the slope for a left-to-right ground surface
                    if node_proj[c1] > node_proj[c2]:
                        _gn = list(reversed(_gn))
                        c1, c2 = _gn[0], _gn[-1]
                    x1, y1 = nodes[c1]
                    x2, y2 = nodes[c2]
                    L = float(np.hypot(x2 - x1, y2 - y1))
                    if L < 1e-12:
                        continue
                    tx, ty = (x2 - x1) / L, (y2 - y1) / L
                    nx, ny = ty, -tx
                    # Rotating the tangent names the INWARD side only when the
                    # load line happens to be written left-to-right. A segment
                    # traced the other way (a downstream pool followed from the
                    # far boundary back toward the toe, say) gets the opposite
                    # normal, and the pressure is applied as UPLIFT — which on a
                    # submerged bench leaves a skin of soil in effective tension
                    # that no deformation can relieve. A boundary edge belongs
                    # to exactly ONE element, so the inside is available
                    # geometrically: point the traction at that element's
                    # centroid and the assembly no longer depends on the order
                    # the load points were typed. Interior edges (two owners)
                    # carry _cen = None and keep the tangent rule.
                    if _cen is not None and (
                            nx * (_cen[0] - 0.5 * (x1 + x2))
                            + ny * (_cen[1] - 0.5 * (y1 + y2))) < 0.0:
                        nx, ny = -nx, -ny
                    if load_dir == 'vertical':
                        # Gravity surcharge (the vendor's 'type: vertical'): the same
                        # traction magnitude, resolved straight DOWN instead of into
                        # the surface. On an inclined crest the surface-normal form
                        # carries a horizontal component of tan(inclination) times the
                        # surcharge — real thrust the load does not have. Integrating
                        # p over the edge LENGTH is unchanged, so the total applied
                        # force has the same magnitude as the LEM slicer's resultant.
                        nx, ny = 0.0, -1.0
                    for tg in ((1.0 - _g) / 2.0, (1.0 + _g) / 2.0):
                        wg = 0.5
                        d = node_proj[c1] * (1.0 - tg) + node_proj[c2] * tg
                        pmag = _p_at(d)
                        if len(_gn) == 3:
                            Nvals = ((2*tg - 1)*(tg - 1), 4*tg*(1 - tg), tg*(2*tg - 1))
                        else:
                            Nvals = (1.0 - tg, tg)
                        for node, Nv in zip(_gn, Nvals):
                            f = pmag * L * wg * Nv
                            nodal_fx[node] = nodal_fx.get(node, 0.0) + f * nx
                            nodal_fy[node] = nodal_fy.get(node, 0.0) + f * ny

                for node, fx in nodal_fx.items():
                    bc_type[node] = 4
                    bc_values[node, 0] += fx
                    bc_values[node, 1] += nodal_fy[node]

                expected_force = np.mean(load_values) * load_total_length
                total_force = np.sqrt(
                    sum(nodal_fx.values())**2 + sum(nodal_fy.values())**2)
                if verbose:
                    print(f"  Distributed load {load_idx}: {len(_load_edges)} edges / "
                          f"{len(load_nodes)} nodes (consistent), "
                          f"total force = {total_force:.1f}, expected ~{expected_force:.1f}")
                continue

            # Pass 2b (fallback when no boundary edges found on the line):
            # tributary lumping per node.
            #
            # The inward direction here needs the SAME geometric rule Pass 2a
            # uses: rotating the tangent 90 degrees CW names the inward side
            # only when the load line happens to be written left-to-right, so a
            # right-to-left line (a downstream pool traced from the far boundary
            # back toward the toe, say) would be applied as UPLIFT. Pass 2a asks
            # the single element that owns the loaded boundary edge; no edge is
            # available here, so ask the same geometry one level coarser: the
            # mean centroid of the elements touching a boundary node lies on the
            # material side, so the normal must point toward it. A node that is
            # interior to the mesh has elements on both sides, its mean centroid
            # sits essentially on top of it, and the test goes inert — the
            # tangent rule then stands, exactly as for Pass 2a's interior edges.
            _load_set = {ni for ni, _ in load_nodes}
            _inside_acc = {}
            for _eidx, _elem in enumerate(elements):
                _et = int(element_types[_eidx])
                if _et not in (3, 4, 6, 8, 9):
                    continue
                _gns = [int(_elem[i]) for i in range(_et)]
                if not any(n in _load_set for n in _gns):
                    continue
                _ncorner = 4 if _et in (4, 8, 9) else 3
                _cen = nodes[[int(_elem[i]) for i in range(_ncorner)]].mean(axis=0)
                for n in _gns:
                    if n in _load_set:
                        _a = _inside_acc.setdefault(n, [np.zeros(2), 0])
                        _a[0] = _a[0] + _cen
                        _a[1] += 1

            # Pass 2: Compute tributary length and load for each node
            # Endpoint nodes extend to the actual line start/end so the
            # full load line length is covered.
            n_load_nodes = len(load_nodes)
            total_trib = 0.0
            for k, (node_idx, proj_dist) in enumerate(load_nodes):
                if n_load_nodes == 1:
                    trib_length = load_total_length
                else:
                    if k == 0:
                        # First node: from line start (0) to midpoint with next node
                        trib_length = (load_nodes[k+1][1] + proj_dist) / 2.0
                    elif k == n_load_nodes - 1:
                        # Last node: from midpoint with prev node to line end
                        trib_length = load_total_length - (proj_dist + load_nodes[k-1][1]) / 2.0
                    else:
                        # Interior node: half-distance to each neighbor
                        trib_length = (load_nodes[k+1][1] - load_nodes[k-1][1]) / 2.0
                total_trib += trib_length

                # Interpolate load value at this position along the load line
                cumulative_length = 0
                load_at_node = load_values[-1]  # default to last value
                for j in range(len(load_coords) - 1):
                    seg_length = np.linalg.norm(np.array(load_coords[j+1]) - np.array(load_coords[j]))
                    cumulative_length += seg_length
                    if proj_dist <= cumulative_length:
                        local_distance = proj_dist - (cumulative_length - seg_length)
                        ratio = local_distance / seg_length if seg_length > 0 else 0
                        load_at_node = load_values[j] * (1 - ratio) + load_values[j+1] * ratio
                        break

                nodal_force_magnitude = load_at_node * trib_length

                # Compute inward normal direction at this point on the load line
                # The tangent is along the load line; rotate 90° CW for inward normal
                # (assumes load line runs left-to-right with slope body below)
                closest_pt = load_linestring.interpolate(proj_dist)
                eps = min(1e-3, load_total_length * 0.01)
                d_back = max(0.0, proj_dist - eps)
                d_fwd = min(load_total_length, proj_dist + eps)
                pt_back = load_linestring.interpolate(d_back)
                pt_fwd = load_linestring.interpolate(d_fwd)
                tx = pt_fwd.x - pt_back.x
                ty = pt_fwd.y - pt_back.y
                t_len = np.sqrt(tx**2 + ty**2)
                if t_len > 1e-15:
                    # Inward normal: rotate tangent 90° clockwise → (ty, -tx)
                    nx = ty / t_len
                    ny = -tx / t_len
                else:
                    nx, ny = 0.0, -1.0  # fallback to vertical

                # Point the normal at the material (see the note above Pass 2b),
                # so the assembled force does not depend on the order the load
                # points were typed.
                _a = _inside_acc.get(node_idx)
                if _a is not None and _a[1]:
                    _to_in = _a[0] / _a[1] - nodes[node_idx]
                    if nx * _to_in[0] + ny * _to_in[1] < 0.0:
                        nx, ny = -nx, -ny

                if load_dir == 'vertical':
                    nx, ny = 0.0, -1.0      # gravity surcharge; see Pass 2a

                # Apply force in inward normal direction (into the slope)
                bc_type[node_idx] = 4  # Applied force
                bc_values[node_idx, 0] = nodal_force_magnitude * nx
                bc_values[node_idx, 1] = nodal_force_magnitude * ny

            # Sanity check: tributary lengths must sum to the full line length
            expected_force = np.mean(load_values) * load_total_length
            total_force = np.sqrt(
                sum(bc_values[ni, 0] for ni, _ in load_nodes)**2 +
                sum(bc_values[ni, 1] for ni, _ in load_nodes)**2
            )
            trib_error = abs(total_trib - load_total_length) / load_total_length
            if trib_error > 0.01:
                warnings.warn(
                    f"Distributed load {load_idx}: tributary lengths sum to {total_trib:.3f} "
                    f"but load line length is {load_total_length:.3f} "
                    f"(error {trib_error:.1%}). Check mesh resolution along load line."
                )
            if n_load_nodes < 2:
                warnings.warn(
                    f"Distributed load {load_idx}: only {n_load_nodes} mesh node(s) found "
                    f"on load line (tolerance={tolerance}). Increase mesh density along load line."
                )
            if verbose:
                print(f"  Distributed load {load_idx}: {n_load_nodes} nodes, "
                      f"total force = {total_force:.1f}, expected ~{expected_force:.1f}, "
                      f"sum(trib) = {total_trib:.2f}")


    # Step 4c: Line loads (v12 'lloads') -> concentrated nodal forces (type 4).
    # The mesh should carry a node exactly at each load point — pass
    # mesh.extract_point_constraints(slope_data) to build_mesh_from_polygons
    # (point_constraints=...) so one is guaranteed; otherwise the force snaps to
    # the nearest node with a warning.
    line_loads_fem = slope_data.get('line_loads') or []
    if line_loads_fem:
        _span = float(np.max(nodes[:, 0]) - np.min(nodes[:, 0])) if len(nodes) else 1.0
        _tol_ll = 1e-3 * max(1.0, _span)
        for _ll_idx, _ll in enumerate(line_loads_fem):
            _d = np.linalg.norm(nodes - np.array([_ll['x'], _ll['y']]), axis=1)
            _i = int(np.argmin(_d))
            if _d[_i] > _tol_ll:
                warnings.warn(
                    f"Line load '{_ll.get('label', _ll_idx + 1)}' at ({_ll['x']}, {_ll['y']}): "
                    f"nearest mesh node is {_d[_i]:.3g} away; the force was applied there. "
                    "Pass mesh.extract_point_constraints(slope_data) as point_constraints "
                    "to build_mesh_from_polygons so a node lands exactly at the load point.")
            _ang = np.radians(_ll.get('angle', -90.0))
            bc_type[_i] = 4
            bc_values[_i, 0] += _ll['P'] * np.cos(_ang)
            bc_values[_i, 1] += _ll['P'] * np.sin(_ang)
            if verbose:
                print(f"  Line load '{_ll.get('label', _ll_idx + 1)}': node {_i} at "
                      f"({nodes[_i,0]:.3f}, {nodes[_i,1]:.3f}), "
                      f"F = ({_ll['P'] * np.cos(_ang):.1f}, {_ll['P'] * np.sin(_ang):.1f})")

    # Material-zone columns for the vertical-overburden integral (K0 option). Built
    # from the same zone polygons and the same moist unit weight the 'ru' option
    # integrates, so the two definitions of "the soil column above this point" can
    # never drift apart.
    try:
        from .mesh import get_material_polygons as _gmp3
        _overburden_columns = []
        for _pi, _pd in enumerate(_gmp3(slope_data)):
            _mid = _pd.get('mat_id')
            _midx = _mid if (_mid is not None and 0 <= _mid < len(materials)) else _pi
            _overburden_columns.append(
                ([(float(x), float(y)) for x, y in _pd['coords']],
                 float(materials[_midx].get('gamma', 0.0))))
    except Exception as _e:
        # A model with no usable zone geometry simply cannot offer K0 initialization;
        # solve_fem raises there rather than silently initializing to zero stress.
        _overburden_columns = []

    # Get other parameters
    unit_weight = require_gamma_water(slope_data, "FEM analysis")
    # SIGN CONVENTION (see the FEM overview page): the FEM analyzes both
    # faces of a dam/levee at once, so unlike the LEM (which takes abs(k) and
    # gets the direction from the failure surface) the USER directs the
    # pseudo-static force with the sign of k in the template - positive k
    # acts in +x (right-facing slopes), negative in -x (left-facing). The
    # value is applied exactly as entered.
    k_seismic = slope_data.get("k_seismic", 0.0)
    
    # Construct fem_data dictionary
    fem_data = {
        "nodes": nodes,
        "node_depth": node_depth,  # depth below ground per node (min_slip_depth filter)
        "elements": elements,
        "element_types": element_types,
        "element_materials": element_materials,
        "bc_type": bc_type,
        "bc_values": bc_values,
        "fixed_nodes": fixed_nodes,
        "roller_x_nodes": roller_x_nodes,
        "roller_y_nodes": roller_y_nodes,
        "c_by_mat": c_by_mat,
        "phi_by_mat": phi_by_mat,
        "E_by_mat": E_by_mat,
        "nu_by_mat": nu_by_mat,
        "gamma_by_mat": gamma_by_mat,
        "material_names": material_names,
        # v16 template-carried run-option defaults. solve_ssrm / solve_fem fall
        # back to these when their tension_cutoff_by_material / elastic_materials
        # (resp. tension_cap_by_elem / elastic_mask) kwargs are left None, so a file
        # that specifies t_cut / option=elastic is honored automatically; an
        # explicit kwarg overrides (pass {} / [] to disable). Blank t_cut is omitted
        # (None -> no cutoff), matching the loader.
        "tension_cutoff_by_material": {
            m["name"]: float(m["t_cut"]) for m in materials
            if m.get("t_cut") is not None},
        "elastic_materials": [
            m["name"] for m in materials
            if str(m.get("option", "")).strip().lower() == "elastic"],
        # v17 template-carried matric-suction strength (Fredlund extended MC),
        # off by default. phi_b (deg) is the unsaturated friction angle that turns
        # matric suction into apparent cohesion in the FEM yield; s_cap (stress
        # units) bounds the credited suction. Blank -> None on the material, mapped
        # here to phi_b = 0 (no suction credit, bit-identical to pre-v17) and
        # s_cap = inf (uncapped). solve_fem / solve_ssrm read these when their
        # suction_phi_b / suction_cap kwargs are left None (an explicit kwarg
        # overrides; pass {} to force suction off regardless of the file). See the
        # suction machinery in solve_fem for how they enter the MC strength.
        "phi_b_by_mat": np.array([
            float(m["phi_b"]) if m.get("phi_b") is not None else 0.0
            for m in materials]),
        "s_cap_by_mat": np.array([
            float(m["s_cap"]) if m.get("s_cap") is not None else np.inf
            for m in materials]),
        # v19 template-carried FEM run options. All three are None/empty on a file
        # that does not declare them (every pre-v19 file), and solve_fem /
        # solve_ssrm read them only when the matching kwarg is left unset — an
        # explicit kwarg always wins, so nothing already scripted changes.
        #   k0            -- at-rest coefficient for the initial stress state
        #   tension_srf   -- reduce the tensile cap with the trial F (None = engine
        #                    default, which is True)
        "k0": slope_data.get("k0"),
        "tension_srf": slope_data.get("tension_srf"),
        # v20 SSR zone overlays, read from the polygon sheet's sentinel Mat ID rows
        # (fileio.SSR_ZONE_SENTINELS). Each entry is {'kind': 'reduce' | 'hold' |
        # 'hold_elastic', 'polygon': [(x, y), ...]}. They are ANALYSIS OVERLAYS, not
        # geometry — the loader keeps them out of slope_data['polygons'], so nothing
        # in the mesh, the material assignment or the slices has seen them. solve_ssrm
        # composes them into the element masks (see _compose_ssr_zone_masks); an
        # explicit ssr_zone kwarg wins. An EMPTY list means "reduce everything".
        "ssr_zones": list(slope_data.get("ssr_zones") or []),
        "c_by_elem": c_by_elem,  # Element-wise cohesion (for c/p option)
        "phi_by_elem": phi_by_elem,  # Element-wise friction angle
        "pow_flag_by_elem": pow_flag_by_elem,  # power-curve elements (tangent-linearized in the VP loop)
        "pow_a_by_elem": pow_a_by_elem,
        "pow_b_by_elem": pow_b_by_elem,
        "pow_cp_by_elem": pow_cp_by_elem,
        "pow_d_by_elem": pow_d_by_elem,
        "hb_flag_by_elem": hb_flag_by_elem,  # Hoek-Brown elements (tangent-linearized in the VP loop)
        "hb_sci_by_elem": hb_sci_by_elem,
        "hb_mb_by_elem": hb_mb_by_elem,      # derived from GSI/mi/D at build time
        "hb_s_by_elem": hb_s_by_elem,
        "hb_a_by_elem": hb_a_by_elem,
        # Material-zone columns for the vertical-overburden integral: (exterior
        # coords, moist gamma) per zone. Carried as plain coordinate lists (cheap to
        # build and to pickle) so the K0 initial-stress option can integrate the soil
        # column above any point WITHOUT needing slope_data back — see
        # _gauss_point_overburden. Same zones and same moist-gamma convention as the
        # 'ru' nodal sigma_v above.
        "overburden_columns": _overburden_columns,
        "u": u,
        "u_signed": u_signed,  # raw (un-clamped) nodal seep field for the suction option; None if no seep suction
        "sigma_v": sigma_v,  # nodal vertical soil stress (ru option; else None)
        "ru_by_mat": np.array([float(m.get("ru", 0.0) or 0.0) for m in materials]),
        "elements_1d": elements_1d,
        "element_types_1d": element_types_1d,
        "element_materials_1d": element_materials_1d,
        "t_allow_by_1d_elem": t_allow_by_1d_elem,
        "t_res_by_1d_elem": t_res_by_1d_elem,
        "k_by_1d_elem": k_by_1d_elem,
        # Per-element geometry along each reinforcement line, and the line labels:
        # what the 1D details panel needs to place and name a line's elements.
        "elem_length_1d": elem_length_1d,
        "dist_end1_1d": dist_end1_1d,
        "dist_end2_1d": dist_end2_1d,
        "reinforce_line_labels": reinforce_line_labels,
        "cos_theta_1d": cos_theta_1d,
        "sin_theta_1d": sin_theta_1d,
        "dof_indices_1d": dof_indices_1d,
        "n_dof_by_1d_elem": n_dof_by_1d_elem,
        "K_global_1d_elems": K_global_1d_elems,
        "unit_weight": unit_weight,
        "k_seismic": k_seismic,
        "pp_option": pp_option,
        "piezo_line_coords": piezo_line_coords,
        "piezo_phreatic": bool(slope_data.get('piezo_phreatic', False)),
        "gamma_water": require_gamma_water(slope_data, "FEM analysis"),
        # DOF offset map (pile nodes get 3 DOFs, others get 2)
        "dof_offset": dof_offset,
        "is_pile_node": is_pile_node,
        "n_dof_total": n_dof_total,
        # Pile beam elements (Euler-Bernoulli)
        "n_pile_elements": n_pile_elements,
        "pile_elem_mask": pile_elem_mask,
        "pile_elem_indices": pile_elem_indices,
        "cos_theta_pile": cos_theta_pile,
        "sin_theta_pile": sin_theta_pile,
        "dof_indices_pile": dof_indices_pile,
        "n_dof_by_pile_elem": n_dof_by_pile_elem,
        "K_global_pile_elems": K_global_pile_elems,
        "K_local_by_pile_elem": K_local_by_pile_elem,
        "V_cap_by_pile_elem": V_cap_by_pile_elem,
        "M_cap_by_pile_elem": M_cap_by_pile_elem,
        "elem_length_by_pile_elem": elem_length_by_pile_elem,
        "S_by_pile_elem": S_by_pile_elem,
        "EI_by_pile_elem": EI_by_pile_elem,
        "EA_by_pile_elem": EA_by_pile_elem,
        "pile_node_pairs": pile_node_pairs,
        "pile_elem_nodes": pile_elem_nodes,
        "pile_line_idx_by_pile_elem": pile_line_idx_by_pile_elem,
        "pile_head_nodes": pile_head_nodes,
        "pile_head_fixed": pile_head_fixed,
        "pile_head_pinned": pile_head_pinned,
        "pile_tip_nodes": pile_tip_nodes,
        "pile_tip_fixed": pile_tip_fixed,
        "pile_tip_pinned": pile_tip_pinned,
    }


    # Carry the unit declaration through as provenance so export_fem_solution can record
    # it in the run meta.json (units plan phase 2, sec 2a). Recorded only when set, so an
    # undeclared model's fem_data is unchanged.
    if slope_data.get("unit_system"):
        fem_data["unit_system"] = slope_data["unit_system"]
    if slope_data.get("time_unit"):
        fem_data["time_unit"] = slope_data["time_unit"]

    return fem_data


# Implementation of Perzyna Visco-Plastic Algorithm for Slope Stability
#
# Based on:
# - Griffiths & Lane (1999) "Slope stability analysis by finite elements"
# - Perzyna (1966) "Fundamental problems in viscoplasticity"
# - Zienkiewicz & Cormeau (1974) visco-plastic algorithm
#
# Key features:
# - Pure non-convergence failure criterion
# - Perzyna stress redistribution algorithm
# - 8-node quadrilateral elements with reduced integration
# - No plastic stiffness reduction

def _resolve_suction_by_elem(fem_data, suction_phi_b, suction_cap, element_materials):
    """Resolve the opt-in matric-suction strength inputs to per-element arrays,
    mirroring the LEM (slice.generate_slices) auto-wire + kwarg-override semantics.

    Returns ``(tanphib_by_elem, scap_by_elem, active)`` where ``tanphib_by_elem``
    is ``tan(phi_b)`` per element (0.0 where a material has no suction angle),
    ``scap_by_elem`` is the per-element suction cap (``inf`` = uncapped), and
    ``active`` is True iff any element carries a positive phi_b (the fast-path
    gate: when False the caller skips ALL suction machinery, so a default run is
    bit-identical to the pre-suction solver).

    ``suction_phi_b`` (name -> deg dict) and ``suction_cap`` (scalar or name ->
    cap dict) default to None = read the v17 template values carried on fem_data
    (``phi_b_by_mat`` / ``s_cap_by_mat``); an explicit kwarg wins over the file,
    and an empty dict forces suction off regardless of the file. A phi_b keyed to
    an unknown material name warns (a typo yields zero suction, not an error) —
    the same policy as the LEM.
    """
    names = list(fem_data.get("material_names", []))
    n_mat = len(names)
    # phi_b per material (degrees): explicit dict, else the file's phi_b_by_mat.
    if suction_phi_b is None:
        phib_by_mat = np.asarray(
            fem_data.get("phi_b_by_mat", np.zeros(n_mat)), dtype=float)
    else:
        phib_by_mat = np.zeros(n_mat)
        for nm, deg in suction_phi_b.items():
            if nm in names:
                phib_by_mat[names.index(nm)] = float(deg) if deg else 0.0
            else:
                warnings.warn(
                    f"suction_phi_b names material '{nm}', which is not in the "
                    f"model (materials: {names}). No suction strength applied for it.")
    # suction cap per material (stress units): scalar (one cap for all), name-dict,
    # else the file's s_cap_by_mat. inf = uncapped.
    if suction_cap is None:
        scap_by_mat = np.asarray(
            fem_data.get("s_cap_by_mat", np.full(n_mat, np.inf)), dtype=float)
    elif isinstance(suction_cap, dict):
        scap_by_mat = np.full(n_mat, np.inf)
        for nm, cap in suction_cap.items():
            if nm in names and cap is not None:
                scap_by_mat[names.index(nm)] = float(cap)
    else:
        scap_by_mat = np.full(n_mat, float(suction_cap))
    tanphib_by_mat = np.tan(np.radians(phib_by_mat))
    tanphib_by_elem = tanphib_by_mat[element_materials - 1]
    scap_by_elem = scap_by_mat[element_materials - 1]
    active = bool(np.any(tanphib_by_elem > 0.0))
    return tanphib_by_elem, scap_by_elem, active


def _suction_capped(u_signed_gp, scap):
    """The suction credited at a Gauss point, before phi_b and the trial F.

    Fredlund's ``s = max(0, -u_w)`` on the SIGNED pore pressure -- not the
    clamped ``u >= 0`` the effective-normal term reads, since the whole point is
    the negative pressure above the water table the clamp throws away -- bounded
    by the material's ``s_cap`` (``inf`` = uncapped).

    F-INDEPENDENT, which is why it is separated from the product below: the
    Newton path caches it per group and re-forms only the ``tan(phi_b) * s / F``
    factor when a ramp step changes the strength.
    """
    return np.minimum(np.maximum(-u_signed_gp, 0.0), scap)


def _suction_apparent_cohesion(tanphib, s_capped, Finv):
    """Fredlund's apparent cohesion at a Gauss point, reduced by the trial F.

        c_suction,r = tan(phi_b) * min(max(0, -u_w), s_cap) / F

    added to c' in the Mohr-Coulomb envelope. It IS divided by the trial strength,
    alongside c' and tan(phi'), and by the same per-element factor -- so an
    ssr_exclude element, whose F is 1.0, keeps its whole credit. That is what
    distinguishes it from the tensile cap, which is reduced only when
    ``tension_srf`` says so.

    BOTH drivers call this, so the only thing a driver-against-driver comparison
    can find is a plumbing difference: which field, which cap, which F.
    """
    return tanphib * s_capped * Finv


# UNUSED (see solve_fem). Cormeau's limits are material-point bounds on the local
# viscoplastic rate recurrence — they contain only E, nu and phi, and are element
# INDEPENDENT, which is why Smith & Griffiths tabulate them per yield criterion
# rather than per element. An earlier comment here claimed the safety factor
# corrected for Cormeau assuming constant-strain triangles; that was fabricated and
# is withdrawn. The real unstated caveat is that the Mohr-Coulomb limit assumes
# ASSOCIATED flow, while XSLOPE runs psi = 0.
# The activation threshold is also nu-dependent, not universal: the clamp bites when
# sin^2(phi_r) > 1.7(1-2nu), i.e. phi_r > 67 deg at nu=0.25 but > 24 deg at nu=0.45.
_DT_CORMEAU_SAFETY = 0.9


def _dt_stability_clamp(dt_base, E_by_elem, nu_by_elem, phi_reduced_rad,
                        safety=_DT_CORMEAU_SAFETY):
    """Clamp the viscoplastic pseudo-timestep to the Zienkiewicz & Cormeau (1974)
    stability limit evaluated at the CURRENT (F-reduced) friction angle.

        dt <= 4(1+nu)(1-2nu) / (E * (1 - 2nu + sin^2(phi)))

    Why this cannot live in the cached prepared model: the limit tightens as
    sin^2(phi) grows, and SSRM raises the reduced angle atan(tan(phi)/F) without
    bound as F falls -- a 35 deg soil is 60 deg at F = 0.4 and 74 deg at F = 0.2.
    The phi-free bound used elsewhere (4(1+nu)/3E) is safe only up to ~63 deg; past
    that the viscoplastic iteration oscillates and never converges, which a
    non-convergence failure criterion then scores as COLLAPSE. That is how a stable
    slope earns an absurd factor of safety: RS2-62c returned 0.208 on one mesh
    because every trial below F ~ 0.45 stalled for this reason while max|u| never
    moved off its elastic value.

    Returns min(dt_base, safety * limit), so dt only ever shrinks and only where the
    caller's value is genuinely unsafe -- at ordinary reduced angles the base value
    is already well inside the limit and is returned untouched, keeping prior
    results bit-identical.
    """
    s2 = np.sin(phi_reduced_rad) ** 2
    denom = E_by_elem * (1.0 - 2.0 * nu_by_elem + s2)
    limit = np.min(4.0 * (1.0 + nu_by_elem) * (1.0 - 2.0 * nu_by_elem) / denom)
    return min(dt_base, safety * float(limit))


def _gauss_point_overburden(fem_data, elem_gp_data):
    """Total vertical soil stress from the overburden at every Gauss point.

    Returns ``sv0[elem_idx][gp_idx]`` = the MAGNITUDE (>= 0) of the total vertical
    stress carried by the soil column directly above that Gauss point, integrated by
    intersecting a vertical ray with the material-zone polygons carried on
    ``fem_data["overburden_columns"]``. Multi-band zones are handled exactly, and
    moist gamma is used throughout — the same definition (and the same code shape)
    as the 'ru' option's nodal ``sigma_v``, so the two can never disagree about what
    "the column above this point" means.

    DELIBERATELY EXCLUDED, matching that definition (Bishop & Morgenstern): applied
    surface tractions (a reservoir load, a distributed load, a footing) and water
    above ground. They are not soil overburden; they enter the solution as boundary
    forces during the equilibrium iteration, which is where a load applied AFTER the
    in-situ state belongs.
    """
    cols = fem_data.get("overburden_columns") or []
    if not cols:
        raise ValueError(
            "K0 initial stress requires material-zone geometry to integrate the "
            "overburden, and this model carries none. Rebuild fem_data from a "
            "slope_data with profile lines or polygons.")
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    polys = [(Polygon(coords), gamma) for coords, gamma in cols]
    y_top = float(np.max(nodes[:, 1])) + 1.0

    sv0 = []
    for elem_idx in range(len(elem_gp_data)):
        elem_nodes_idx = elements[elem_idx][:element_types[elem_idx]]
        elem_coords = nodes[elem_nodes_idx]
        gp_list = []
        for gp_data in elem_gp_data[elem_idx]:
            N = gp_data['N']
            x_gp = float(N @ elem_coords[:, 0])
            y_gp = float(N @ elem_coords[:, 1])
            ray = LineString([(x_gp, y_gp), (x_gp, y_top)])
            sv = 0.0
            for poly, gamma in polys:
                minx, _miny, maxx, maxy = poly.bounds
                if x_gp < minx or x_gp > maxx or y_gp >= maxy:
                    continue
                inter = ray.intersection(poly)
                if not inter.is_empty:
                    sv += gamma * inter.length
            gp_list.append(max(0.0, sv))
        sv0.append(gp_list)
    return sv0


class _CholmodFactorSolver:
    """``.solve(b)`` over a CHOLMOD Cholesky factor, matching SuperLU's interface."""

    kind = "CHOLMOD Cholesky"

    def __init__(self, factor):
        self._factor = factor

    def solve(self, b):
        return self._factor(b)


def _factorize_free_stiffness(K_free):
    """Factorize the free-free stiffness once, for reuse by every iteration.

    After the constrained DOFs are removed, K is symmetric positive definite, so
    the general unsymmetric LU that ``splu`` performs by default stores and
    applies two triangular factors where one would do -- and the back-
    substitution is roughly 60 % of every viscoplastic pass. Two symmetric paths
    are used instead, in order of preference:

    * CHOLMOD (``scikit-sparse``), a true supernodal Cholesky, if it is
      importable. It is an optional accelerator, never a requirement.
    * ``splu``'s documented symmetric mode -- ``permc_spec='MMD_AT_PLUS_A'``
      with the diagonal pivot threshold at zero and ``SymmetricMode`` on -- which
      orders for the symmetric pattern and pivots down the diagonal.

    Either way the matrix itself is unchanged, and a factorization that fails
    (an indefinite or singular K, which means a badly posed model rather than a
    bad solver choice) falls back to the general LU so the failure is reported
    where it always was. ``XSLOPE_FEM_FACTOR`` (``auto``, ``cholmod``,
    ``symmetric``, ``unsymmetric``) pins the path for measurement.

    Returns (factor, kind) where factor exposes ``.solve(b)``.
    """
    # Default is the symmetric path (the full FEM suite ran on it 2026-08-23:
    # 186/189, the two genuine movers being knife-edge models whose near-critical
    # answer is round-off sensitive under any configuration — their locks carry
    # that note). XSLOPE_FEM_FACTOR pins a specific path when reproducing an
    # older configuration.
    mode = os.environ.get("XSLOPE_FEM_FACTOR", "auto").strip().lower()

    if mode in ("auto", "cholmod"):
        try:
            from sksparse.cholmod import cholesky as _cholmod_cholesky
        except Exception:
            if mode == "cholmod":
                raise
        else:
            try:
                return (_CholmodFactorSolver(_cholmod_cholesky(K_free.tocsc())),
                        _CholmodFactorSolver.kind)
            except Exception:
                if mode == "cholmod":
                    raise

    if mode in ("auto", "symmetric"):
        try:
            return (splu(K_free.tocsc(), permc_spec="MMD_AT_PLUS_A",
                         diag_pivot_thresh=0.0,
                         options=dict(SymmetricMode=True)),
                    "SuperLU symmetric mode")
        except Exception:
            if mode == "symmetric":
                raise

    return splu(K_free.tocsc()), "SuperLU general LU"


def _prepare_fem_model(fem_data, *, dt_scale=1.0, suction_phi_b=None,
                       suction_cap=None, elastic_mask=None,
                       tension_cap_by_elem=None, tension_cutoff=False,
                       min_slip_depth=None, k0=None, debug_level=0):
    """Build the strength-reduction-factor-INDEPENDENT setup for solve_fem ONCE.

    Everything a solve_fem call assembles that does NOT depend on the trial
    strength-reduction factor F lives here: the global stiffness assembly and its
    LU factorization, the gravity load vector (+ boundary forces), the free-DOF
    partition, the element/Gauss-point geometry (elem_gp_data), the pore-pressure
    fields (u_gp / u_gp_signed), the F-independent halves of the vectorized
    Gauss-point GROUPS (B, D4, weights, dof maps, the viscoplastic dt_r, and the
    per-GP suction / elastic / power-curve / Hoek-Brown parameter arrays), and the
    Dawson–Roth–Drescher per-node body-force normalization (g_node). solve_ssrm
    builds this ONE prepared-model dict and threads it into every trial's solve_fem
    via the internal ``_prepared=`` kwarg, so the ~10 bisection trials share a
    single K factorization and one geometry precompute instead of rebuilding them.

    Crucially the returned dict carries NO per-solve state and NO F-dependent
    array: the working Gauss-point groups (with their F-reduced c_r/phi_r, the
    accumulated viscoplastic strains evp, and the per-stage pore-pressure loads)
    are rebuilt fresh inside every solve_fem call from these static pieces, so a
    reused prepared model can never serve a stale strength or a stale plastic
    state. A standalone solve_fem builds its own prepared model each call, so its
    behavior is bit-identical to before this cache existed.

    The kwargs here mirror the same-named solve_fem parameters and must match the
    values the trials will use (solve_ssrm holds them fixed across the bisection).
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    element_materials = fem_data["element_materials"]
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]
    fixed_nodes = fem_data.get("fixed_nodes", set())
    roller_x_nodes = fem_data.get("roller_x_nodes", set())
    roller_y_nodes = fem_data.get("roller_y_nodes", set())

    # Template-carried defaults (v16), resolved exactly as solve_fem's direct-call
    # fallback does (an explicit per-element array wins). solve_ssrm always passes
    # explicit arrays, so this only fires on a standalone solve_fem's prepared model.
    if tension_cap_by_elem is None:
        _tc_by_mat = fem_data.get("tension_cutoff_by_material")
        if _tc_by_mat:
            _names = list(fem_data.get("material_names", []))
            tension_cap_by_elem = np.full(len(element_materials), np.inf)
            for _nm, _T in _tc_by_mat.items():
                if _nm in _names:
                    tension_cap_by_elem[element_materials == _names.index(_nm) + 1] = float(_T)
    if elastic_mask is None:
        _el_names = fem_data.get("elastic_materials")
        if _el_names:
            _names = list(fem_data.get("material_names", []))
            _ids = [_names.index(_nm) + 1 for _nm in _el_names if _nm in _names]
            if _ids:
                elastic_mask = np.isin(element_materials, _ids)
        # v20 "SSR elastic" (-3) overlay rows join the elastic set by union — the
        # polygon-addressed twin of elastic_materials. Same direct-call-only rule:
        # solve_ssrm resolves and passes an explicit mask, and strips 'ssr_zones'
        # from the trial fem_data, so this cannot double-apply.
        _, _zone_elastic = _compose_ssr_zone_masks(fem_data, fem_data.get("ssr_zones") or [])
        if _zone_elastic is not None:
            elastic_mask = (_zone_elastic if elastic_mask is None
                            else (np.asarray(elastic_mask, dtype=bool) | _zone_elastic))

    pow_flag_by_elem = fem_data.get("pow_flag_by_elem")
    if pow_flag_by_elem is None:
        pow_flag_by_elem = np.zeros(len(elements), dtype=bool)
    has_pow = bool(np.any(pow_flag_by_elem))
    hb_flag_by_elem = fem_data.get("hb_flag_by_elem")
    if hb_flag_by_elem is None:
        hb_flag_by_elem = np.zeros(len(elements), dtype=bool)
    has_hb = bool(np.any(hb_flag_by_elem))
    E_by_mat = fem_data["E_by_mat"]
    nu_by_mat = fem_data["nu_by_mat"]
    gamma_by_mat = fem_data["gamma_by_mat"]
    k_seismic = fem_data.get("k_seismic", 0.0)

    n_nodes = len(nodes)
    n_elements = len(elements)

    dof_offset = fem_data.get("dof_offset", None)
    if dof_offset is not None:
        n_dof = int(dof_offset[n_nodes])
    else:
        n_dof = 2 * n_nodes

    # Base Rankine tensile cap per element (tension_cutoff -> 0, then per-element
    # overrides). The tension_srf F-scaling is applied per solve in solve_fem, so
    # this base is F-independent.
    t_cap_base = np.full(n_elements, np.inf)
    if tension_cutoff:
        t_cap_base[:] = 0.0
    if tension_cap_by_elem is not None:
        _tc = np.asarray(tension_cap_by_elem, dtype=float)
        _have = np.isfinite(_tc)
        t_cap_base[_have] = _tc[_have]

    elastic_by_elem = (np.asarray(elastic_mask, dtype=bool)
                       if elastic_mask is not None else None)

    suction_tanphib_by_elem, suction_scap_by_elem, suction_active = \
        _resolve_suction_by_elem(fem_data, suction_phi_b, suction_cap, element_materials)

    # ---- K_global (elastic, constant) ----
    K_global = build_global_stiffness(nodes, elements, element_types,
                                      element_materials, E_by_mat, nu_by_mat,
                                      fem_data=fem_data)

    # ---- Gravity load vector ----
    F_gravity = build_gravity_loads(nodes, elements, element_types,
                                    element_materials, gamma_by_mat, k_seismic,
                                    fem_data=fem_data)
    for i in range(n_nodes):
        if bc_type[i] == 4:  # Force boundary condition
            dof_x = dof_offset[i] if dof_offset is not None else 2 * i
            F_gravity[dof_x] += bc_values[i, 0]
            F_gravity[dof_x + 1] += bc_values[i, 1]

    # ---- Free/constrained DOFs ----
    n_pile_elements = fem_data.get("n_pile_elements", 0)
    has_pile_elements = n_pile_elements > 0
    pile_head_nodes = fem_data.get("pile_head_nodes", np.array([], dtype=int))
    pile_head_fixed = fem_data.get("pile_head_fixed", np.array([], dtype=bool))
    pile_head_pinned = fem_data.get("pile_head_pinned", np.array([], dtype=bool))
    pile_tip_nodes = fem_data.get("pile_tip_nodes", np.array([], dtype=int))
    pile_tip_fixed = fem_data.get("pile_tip_fixed", np.array([], dtype=bool))
    pile_tip_pinned = fem_data.get("pile_tip_pinned", np.array([], dtype=bool))

    constraint_dofs = []
    for i in range(n_nodes):
        dof_x = dof_offset[i] if dof_offset is not None else 2 * i
        dof_y = dof_x + 1
        if bc_type[i] == 1 or i in fixed_nodes:  # Fixed
            constraint_dofs.extend([dof_x, dof_y])
        elif bc_type[i] == 2 or i in roller_x_nodes:  # X-roller
            constraint_dofs.append(dof_x)
        elif bc_type[i] == 3 or i in roller_y_nodes:  # Y-roller
            constraint_dofs.append(dof_y)
    if has_pile_elements:
        # Each pile end has four states, the same at the head and the tip:
        # 'free' (nothing held — the soil, or the boundary the node sits on,
        # decides), 'pinned' (translations held, rotation free), 'unrotated'
        # (rotation held, translations free), 'fixed' (all three held). A held
        # rotation constrains the third DOF a pile node carries; held
        # translations constrain its two displacement DOFs.
        for _end_nodes, _rot_held, _trans_held in (
                (pile_head_nodes, pile_head_fixed, pile_head_pinned),
                (pile_tip_nodes, pile_tip_fixed, pile_tip_pinned)):
            for pe_idx in range(len(_end_nodes)):
                pe_node = _end_nodes[pe_idx]
                base = dof_offset[pe_node] if dof_offset is not None else 2 * pe_node
                if _rot_held[pe_idx] and dof_offset is not None:
                    constraint_dofs.append(base + 2)
                if len(_trans_held) > pe_idx and _trans_held[pe_idx]:
                    constraint_dofs.append(base)
                    constraint_dofs.append(base + 1)

    constraint_set = set(constraint_dofs)
    free_dofs = np.array(sorted(set(range(n_dof)) - constraint_set))
    n_free = len(free_dofs)

    # ---- Extract K_free and PRE-FACTORIZE ----
    # Constrained DOFs are eliminated by taking the free-free submatrix, exactly
    # as before -- the selection is done sparsely, so K never exists as a dense
    # n_dof x n_dof array (1.4 GB at the 2 ft mesh, 21.8 GB at 1 ft).
    K_csr = K_global if issparse(K_global) else csr_matrix(K_global)
    K_csr = K_csr.tocsr()
    K_free = K_csr[free_dofs][:, free_dofs]
    K_free.eliminate_zeros()
    _tp = time.perf_counter() if _PROF_ON else None
    K_factor, K_factor_kind = _factorize_free_stiffness(K_free)
    if _PROF_ON:
        _prof_add("vp_factorize", _tp)

    if debug_level >= 1:
        print(f"  DOFs: {n_dof} total, {n_free} free, {len(constraint_dofs)} constrained")
        print(f"  K factorized ({K_factor_kind}, reused for all iterations)")

    # ---- dt from material properties (Smith & Griffiths) + Rankine dt_r ----
    # NOTE: this is the phi-INDEPENDENT part only. It is a safe bound while the
    # friction angle stays modest, but the true viscoplastic stability limit
    # (Zienkiewicz & Cormeau 1974) tightens as sin^2(phi) grows, and SSRM drives the
    # REDUCED angle atan(tan(phi)/F) very high at low F -- 35 deg becomes 74 deg at
    # F = 0.2. Above ~63 deg this bound is no longer safe, the iteration oscillates
    # forever, and a non-convergence failure criterion scores that as collapse. So
    # solve_fem clamps this value per trial F (see _dt_stability_clamp); the element
    # E/nu needed for that clamp are cached here alongside dt.
    dt = 1.0e15
    dt_r = np.zeros(n_elements)
    E_by_elem_dt = np.zeros(n_elements)
    nu_by_elem_dt = np.zeros(n_elements)
    for elem_idx in range(n_elements):
        mat_id = element_materials[elem_idx] - 1
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        E_by_elem_dt[elem_idx] = E
        nu_by_elem_dt[elem_idx] = nu
        ddt = 4.0 * (1.0 + nu) / (3.0 * E)
        if ddt < dt:
            dt = ddt
        dt_r[elem_idx] = 0.1 * (1.0 + nu) * (1.0 - 2.0 * nu) / (E * (1.0 - nu))
    dt *= dt_scale
    if debug_level >= 2:
        print(f"  dt = {dt:.3e}")

    # ---- Pre-compute element B / D matrices at Gauss points ----
    gauss_points_2x2, gauss_weights_2x2 = get_gauss_points_2x2()
    elem_gp_data = []
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        D = build_constitutive_matrix(E, nu)
        D4 = build_constitutive_matrix_4(E, nu)

        elem_nodes_idx = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes_idx]

        gp_list = []
        dof_indices = _elem_dof_indices(elem_nodes_idx, dof_offset=dof_offset)
        if elem_type == 3:
            B, area = compute_B_matrix_triangle(elem_coords)
            N = np.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
            gp_list.append({'B': B, 'weight': area, 'D': D, 'D4': D4, 'dof_indices': dof_indices, 'N': N})
        elif elem_type == 6:
            tri_gp, tri_wt = get_gauss_points_tri3()
            for gp_idx in range(3):
                L1, L2, L3 = tri_gp[gp_idx]
                B, det_J = _compute_B_and_detJ_tri6(elem_coords, L1, L2, L3)
                weight = 0.5 * abs(det_J) * tri_wt[gp_idx]
                N = compute_tri6_shape_functions(L1, L2, L3)
                gp_list.append({'B': B, 'weight': weight, 'D': D, 'D4': D4, 'dof_indices': dof_indices, 'N': N})
        elif elem_type == 4:
            for gp_idx in range(4):
                xi, eta = gauss_points_2x2[gp_idx]
                B, det_J = _compute_B_and_detJ_quad4(elem_coords, xi, eta)
                weight = gauss_weights_2x2[gp_idx] * abs(det_J)
                N = compute_quad4_shape_functions(xi, eta)
                gp_list.append({'B': B, 'weight': weight, 'D': D, 'D4': D4, 'dof_indices': dof_indices, 'N': N})
        elif elem_type == 8:
            for gp_idx in range(4):
                xi, eta = gauss_points_2x2[gp_idx]
                B, det_J = _compute_B_and_detJ_quad8(elem_coords, xi, eta)
                weight = gauss_weights_2x2[gp_idx] * abs(det_J)
                N = compute_quad8_shape_functions(xi, eta)
                gp_list.append({'B': B, 'weight': weight, 'D': D, 'D4': D4, 'dof_indices': dof_indices, 'N': N})
        elif elem_type == 9:
            q9_gp, q9_wt = get_gauss_points_3x3()
            for gp_idx in range(9):
                xi, eta = q9_gp[gp_idx]
                B, det_J = _compute_B_and_detJ_quad9(elem_coords, xi, eta)
                weight = q9_wt[gp_idx] * abs(det_J)
                N = compute_quad9_shape_functions(xi, eta)
                gp_list.append({'B': B, 'weight': weight, 'D': D, 'D4': D4, 'dof_indices': dof_indices, 'N': N})

        elem_gp_data.append(gp_list)

    # ---- Pore pressure at each Gauss point ----
    u_nodes = fem_data.get("u", np.zeros(n_nodes))
    pp_option = fem_data.get("pp_option", "none")

    u_gp = []
    if pp_option == "none":
        for elem_idx in range(n_elements):
            u_gp.append([0.0] * len(elem_gp_data[elem_idx]))
    elif pp_option == "piezo":
        piezo_line_coords = fem_data.get("piezo_line_coords", None)
        gamma_water = require_gamma_water(fem_data, "FEM pore-pressure assembly")
        if piezo_line_coords and len(piezo_line_coords) >= 2:
            px = np.array([p[0] for p in piezo_line_coords], dtype=float)
            py = np.array([p[1] for p in piezo_line_coords], dtype=float)
            order = np.argsort(px)
            px, py = px[order], py[order]
            _phreatic = bool(fem_data.get('piezo_phreatic', False))
            # No piezometric surface outside the line's x-extent (see the nodal path
            # in build_fem_data), and reading zero there would silently delete pore
            # pressure -> raise. Guarded here too, independently of the nodal check,
            # so a fem_data assembled elsewhere cannot slip a short line past it.
            _pz_tol = _piezo_extent_tol(px)
            for elem_idx in range(n_elements):
                elem_type = element_types[elem_idx]
                elem_nodes_idx = elements[elem_idx][:elem_type]
                elem_coords = nodes[elem_nodes_idx]
                gp_u_list = []
                for gp_data in elem_gp_data[elem_idx]:
                    N = gp_data['N']
                    x_gp = N @ elem_coords[:, 0]
                    y_gp = N @ elem_coords[:, 1]
                    if x_gp < px[0] - _pz_tol or x_gp > px[-1] + _pz_tol:
                        _check_piezo_extent(x_gp, px, py,
                                            f"Gauss point of element {elem_idx}")
                    piezo_elev = float(np.interp(x_gp, px, py))
                    u_val = max(0.0, gamma_water * (piezo_elev - y_gp))
                    if _phreatic and u_val > 0.0:
                        u_val *= float(_piezo_cos2(x_gp, px, py))
                    gp_u_list.append(u_val)
                u_gp.append(gp_u_list)
        else:
            # u='piezo' and no line to read (see build_fem_data): refused, not
            # silently zeroed. Guarded here too, so a fem_data assembled elsewhere
            # cannot reach the solver with no water at all.
            raise _no_piezo_line_error(None)
    elif pp_option == "seep":
        for elem_idx in range(n_elements):
            elem_type = element_types[elem_idx]
            elem_nodes_idx = elements[elem_idx][:elem_type]
            u_elem_nodes = u_nodes[elem_nodes_idx]
            gp_u_list = []
            for gp_data in elem_gp_data[elem_idx]:
                N = gp_data['N']
                gp_u_list.append(max(0.0, float(N @ u_elem_nodes)))
            u_gp.append(gp_u_list)
    elif pp_option == "ru":
        sigma_v_nodes = fem_data.get("sigma_v")
        ru_by_mat = fem_data.get("ru_by_mat")
        if sigma_v_nodes is None or ru_by_mat is None:
            for elem_idx in range(n_elements):
                u_gp.append([0.0] * len(elem_gp_data[elem_idx]))
        else:
            ru_by_elem_pp = ru_by_mat[element_materials - 1]
            for elem_idx in range(n_elements):
                ru_e = float(ru_by_elem_pp[elem_idx])
                elem_type = element_types[elem_idx]
                elem_nodes_idx = elements[elem_idx][:elem_type]
                sv_elem = sigma_v_nodes[elem_nodes_idx]
                gp_u_list = []
                for gp_data in elem_gp_data[elem_idx]:
                    N = gp_data['N']
                    gp_u_list.append(max(0.0, ru_e * float(N @ sv_elem)))
                u_gp.append(gp_u_list)
    else:
        for elem_idx in range(n_elements):
            u_gp.append([0.0] * len(elem_gp_data[elem_idx]))

    # Signed pore-pressure field for the opt-in matric-suction option.
    u_gp_signed = None
    if suction_active and pp_option in ("piezo", "seep"):
        u_gp_signed = []
        if pp_option == "piezo":
            piezo_line_coords = fem_data.get("piezo_line_coords", None)
            gamma_water = require_gamma_water(fem_data, "FEM suction pore-pressure assembly")
            if piezo_line_coords and len(piezo_line_coords) >= 2:
                px = np.array([p[0] for p in piezo_line_coords], dtype=float)
                py = np.array([p[1] for p in piezo_line_coords], dtype=float)
                order = np.argsort(px)
                px, py = px[order], py[order]
                _phreatic = bool(fem_data.get('piezo_phreatic', False))
                _pz_tol = _piezo_extent_tol(px)
                for elem_idx in range(n_elements):
                    elem_type = element_types[elem_idx]
                    elem_nodes_idx = elements[elem_idx][:elem_type]
                    elem_coords = nodes[elem_nodes_idx]
                    gp_list = []
                    for gp_data in elem_gp_data[elem_idx]:
                        N = gp_data['N']
                        x_gp = N @ elem_coords[:, 0]
                        y_gp = N @ elem_coords[:, 1]
                        if x_gp < px[0] - _pz_tol or x_gp > px[-1] + _pz_tol:
                            _check_piezo_extent(x_gp, px, py,
                                                f"Gauss point of element {elem_idx}")
                        piezo_elev = float(np.interp(x_gp, px, py))
                        u_val = gamma_water * (piezo_elev - y_gp)
                        if _phreatic and u_val > 0.0:
                            u_val *= float(_piezo_cos2(x_gp, px, py))
                        gp_list.append(u_val)
                    u_gp_signed.append(gp_list)
            else:
                raise _no_piezo_line_error(None)
        else:  # seep
            u_nodes_signed = fem_data.get("u_signed")
            if u_nodes_signed is None:
                u_nodes_signed = u_nodes
            u_nodes_signed = np.asarray(u_nodes_signed, dtype=float)
            for elem_idx in range(n_elements):
                elem_type = element_types[elem_idx]
                elem_nodes_idx = elements[elem_idx][:elem_type]
                u_elem_nodes = u_nodes_signed[elem_nodes_idx]
                gp_list = [float(gp_data['N'] @ u_elem_nodes)
                           for gp_data in elem_gp_data[elem_idx]]
                u_gp_signed.append(gp_list)

    if debug_level >= 1 and pp_option != "none":
        all_u = [u_gp[e][g] for e in range(n_elements) for g in range(len(u_gp[e]))]
        max_u = max(all_u) if all_u else 0.0
        n_nonzero = sum(1 for v in all_u if v > 0.0)
        print(f"  Pore pressure ({pp_option}): max u_gp = {max_u:.3f}, {n_nonzero}/{len(all_u)} GPs with u > 0")
        if u_gp_signed is not None:
            all_s = [max(0.0, -u_gp_signed[e][g]) for e in range(n_elements)
                     for g in range(len(u_gp_signed[e]))]
            max_s = max(all_s) if all_s else 0.0
            n_suc = sum(1 for v in all_s if v > 0.0)
            print(f"  Matric suction (signed field): max s = {max_s:.3f}, "
                  f"{n_suc}/{len(all_s)} GPs with suction > 0")

    # ---- F-INDEPENDENT halves of the vectorized Gauss-point groups ----
    # The F-reduced strengths (c_r, phi_r, F, snph, csph, t_cap, Finv) and the
    # accumulated viscoplastic strains evp are rebuilt per solve in solve_fem; here
    # we cache only the geometry (B, D4, w, dof), the per-GP dt_r, and the
    # F-independent material parameter arrays. e_idx (the per-GP element index) lets
    # solve_fem gather the F-dependent arrays by fancy-indexing, bit-identically to
    # the original per-GP list comprehensions.
    gp_groups_static = []
    _by_ndof = {}
    for _e in range(n_elements):
        for _g, _gpd in enumerate(elem_gp_data[_e]):
            _nd = len(_gpd['dof_indices'])
            _by_ndof.setdefault(_nd, []).append((_e, _g))
    for _nd, _pairs in _by_ndof.items():
        G = len(_pairs)
        grp = {
            'pairs': _pairs,
            'e_idx': np.array([e for e, g in _pairs], dtype=int),
            'n': G,
            'B': np.array([elem_gp_data[e][g]['B'] for e, g in _pairs]),
            'D4': np.array([elem_gp_data[e][g]['D4'] for e, g in _pairs]),
            'w': np.array([elem_gp_data[e][g]['weight'] for e, g in _pairs]),
            'dof': np.array([elem_gp_data[e][g]['dof_indices'] for e, g in _pairs], dtype=int),
            'dt_r': np.array([dt_r[e] for e, g in _pairs]),
        }
        if suction_active:
            grp['tanphib'] = np.array([suction_tanphib_by_elem[e] for e, g in _pairs])
            grp['scap'] = np.array([suction_scap_by_elem[e] for e, g in _pairs])
        if elastic_by_elem is not None:
            _em = np.array([elastic_by_elem[e] for e, g in _pairs])
            if _em.any():
                grp['elastic'] = _em
                grp['has_elastic'] = True
        if has_pow:
            _pm = np.array([pow_flag_by_elem[e] for e, g in _pairs])
            if _pm.any():
                grp['pow_m'] = _pm
                for _k, _key in (('pow_a', 'pow_a_by_elem'),
                                 ('pow_b', 'pow_b_by_elem'),
                                 ('pow_cp', 'pow_cp_by_elem'),
                                 ('pow_d', 'pow_d_by_elem')):
                    grp[_k] = np.array([fem_data[_key][e] for e, g in _pairs])
        if has_hb:
            _hm = np.array([hb_flag_by_elem[e] for e, g in _pairs])
            if _hm.any():
                grp['hb_m'] = _hm
                for _k, _key in (('hb_sci', 'hb_sci_by_elem'),
                                 ('hb_mb', 'hb_mb_by_elem'),
                                 ('hb_s', 'hb_s_by_elem'),
                                 ('hb_a', 'hb_a_by_elem')):
                    grp[_k] = np.array([fem_data[_key][e] for e, g in _pairs])
        gp_groups_static.append(grp)
    n_total_gp = sum(len(g['pairs']) for g in gp_groups_static)

    # ---- K0 initial stress: per-GP overburden (F-independent) ----
    sv0_gp = None
    if k0 is not None:
        sv0_gp = _gauss_point_overburden(fem_data, elem_gp_data)
        if debug_level >= 1:
            _all = [v for gl in sv0_gp for v in gl]
            print(f"  K0 initial stress: K0 = {float(k0):g}, max overburden "
                  f"sigma_v = {max(_all) if _all else 0.0:.3f}")

    # ---- Displacement-limit yardstick (F-independent mesh height only) ----
    # vp_disp_limit itself is deliberately NOT cached: it is max_disp_factor *
    # mesh_height, and max_disp_factor DIFFERS across the calls that share a prepared
    # model (the SSRM trials vs the at-failure capture solve). solve_fem recomputes it
    # per call from this mesh_height, so a reused prepared model cannot serve a stale
    # displacement limit.
    mesh_height = float(np.max(nodes[:, 1]) - np.min(nodes[:, 1]))

    # ---- Per-node out-of-balance normalization (Dawson, Roth & Drescher 1999) ----
    if dof_offset is not None:
        node_dof_x = np.array([dof_offset[i] for i in range(n_nodes)], dtype=int)
    else:
        node_dof_x = 2 * np.arange(n_nodes, dtype=int)
    node_dof_y = node_dof_x + 1

    free_dof_mask = np.zeros(n_dof, dtype=bool)
    free_dof_mask[free_dofs] = True
    node_has_free = free_dof_mask[node_dof_x] | free_dof_mask[node_dof_y]

    _elem_w = np.zeros(n_nodes)
    _k_fac = float(np.sqrt(1.0 + k_seismic ** 2))
    for _e in range(n_elements):
        _et = int(element_types[_e])
        _en = [int(v) for v in elements[_e][:_et]]
        if not _en:
            continue
        _corners = _en[:3] if _et in (3, 6) else _en[:4]
        _xy = nodes[_corners]
        _area = 0.5 * abs(np.dot(_xy[:, 0], np.roll(_xy[:, 1], -1))
                          - np.dot(_xy[:, 1], np.roll(_xy[:, 0], -1)))
        _gam = float(gamma_by_mat[int(element_materials[_e]) - 1])
        _w = _gam * _area * _k_fac / len(_en)
        for _nd in _en:
            _elem_w[_nd] += _w
    g_node = _elem_w
    _g_typ = float(np.median(g_node[g_node > 0])) if np.any(g_node > 0) else 0.0
    _f_ext_max = float(np.max(np.sqrt(F_gravity[node_dof_x] ** 2
                                      + F_gravity[node_dof_y] ** 2))) if n_nodes else 0.0
    _floor = max(1e-3 * _g_typ, 1e-3 * _f_ext_max, 1e-30)
    g_node_den = np.maximum(g_node, _floor)

    _deep_free_mask = None
    _deep_dof_mask = None
    _skin_elem_mask = None
    if min_slip_depth is not None and float(min_slip_depth) > 0:
        _nd = fem_data.get("node_depth")
        if _nd is not None:
            _nd_free = np.asarray(_nd, dtype=float)[node_has_free]
            _deep_free_mask = _nd_free >= float(min_slip_depth)
            if not _deep_free_mask.any():
                raise ValueError(
                    f"min_slip_depth={float(min_slip_depth):g} excludes every node — the "
                    f"deepest lies {float(_nd_free.max()) if _nd_free.size else 0.0:.3g} "
                    f"below the ground surface. Reduce min_slip_depth, or check that a "
                    f"ground surface is defined.")
            # The same filter, in the two other shapes a driver needs it in.
            #
            # `deep_dof_mask` carries it to the DEGREES OF FREEDOM: the x and y
            # translations of every deep node that has at least one free one. The
            # Newton path reads its residual norms and its displacement bound
            # through this, so the skin the filter excludes cannot decide a verdict
            # there either — see the note above _nr_equilibrate's merit.
            _dn = np.flatnonzero(node_has_free)[_deep_free_mask]
            _deep_dof_mask = np.zeros(n_dof, dtype=bool)
            _deep_dof_mask[node_dof_x[_dn]] = True
            _deep_dof_mask[node_dof_y[_dn]] = True
            _deep_dof_mask &= free_dof_mask.astype(bool)
            # `skin_elem_mask` carries it to the ELEMENTS, and an element counts as
            # skin only when EVERY one of its nodes is shallower than the filter, so
            # anything straddling the depth keeps the consistent tangent.
            _nd_all = np.asarray(_nd, dtype=float)
            _en = np.asarray(elements)[:, :]
            _skin_elem_mask = np.zeros(n_elements, dtype=bool)
            for _e in range(n_elements):
                _ids = [int(_n) - 1 for _n in _en[_e] if int(_n) > 0]
                _skin_elem_mask[_e] = bool(np.all(_nd_all[_ids] < float(min_slip_depth)))

    return {
        "K_factor": K_factor,
        "F_gravity": F_gravity,
        "free_dofs": free_dofs,
        "n_free": n_free,
        "n_dof": n_dof,
        "n_constrained": len(constraint_dofs),
        "dt": dt,
        "E_by_elem_dt": E_by_elem_dt,
        "nu_by_elem_dt": nu_by_elem_dt,
        "dt_r": dt_r,
        "elem_gp_data": elem_gp_data,
        "u_gp": u_gp,
        "u_gp_signed": u_gp_signed,
        "pp_option": pp_option,
        "gp_groups_static": gp_groups_static,
        "n_total_gp": n_total_gp,
        "mesh_height": mesh_height,
        "node_dof_x": node_dof_x,
        "node_dof_y": node_dof_y,
        "free_dof_mask": free_dof_mask,
        "node_has_free": node_has_free,
        "g_node_den": g_node_den,
        "deep_free_mask": _deep_free_mask,
        "deep_dof_mask": _deep_dof_mask,
        "skin_elem_mask": _skin_elem_mask,
        "suction_active": suction_active,
        "suction_tanphib_by_elem": suction_tanphib_by_elem,
        "suction_scap_by_elem": suction_scap_by_elem,
        "t_cap_base": t_cap_base,
        "elastic_by_elem": elastic_by_elem,
        # K0 initial stress: the per-Gauss-point overburden magnitude, computed
        # once here (it is F-independent, so the SSRM bisection shares it) and None
        # when K0 initialization is off — which is the default, and why an ordinary
        # run is bit-identical to before.
        "sv0_gp": sv0_gp,
    }


# ===================== Hybrid failure criterion (opt-in) =====================
# Non-convergence on its own is a statement about the SOLVER, not about the slope.
# The hybrid criterion keeps non-convergence as the trigger but requires
# DISPLACEMENT EVIDENCE before a non-converged trial is called a failure, using the
# only yardstick available inside a single solve: the trial's own elastic response,
# which solve_fem already computes (`displacements_elastic`).
#
# CALIBRATION (measured, 2026-07-26; see docs/fem/overview.md "Hybrid" and the
# RS2-62 verification section for the underlying evidence):
#
#   * Stable-but-numerically-stuck trials sit at max|u| = 1.0-1.1 x the elastic
#     displacement and are FROZEN there — identical to three decimals whether the
#     iteration budget is 10,000 or 80,000. Nothing is moving; the solve simply
#     cannot drive the out-of-balance residual under an absolute tolerance.
#   * Genuinely failing trials reach 4-21 x elastic and are STILL GROWING when the
#     budget runs out (they run into the displacement cap when one is set).
#   * Griffiths & Lane Example 1, the reference case for the displacement-vs-F
#     upturn (docs/fem/images/griffiths1_sweep.png), sits at ~1.5-3 x elastic and
#     growing through the bisection-relevant band around F = 1.40.
#
# The thresholds below place the STUCK ceiling above the measured frozen band with
# headroom (1.25) and the FAILED floor at the bottom of the G&L upturn band (1.5),
# leaving 1.25-1.5 deliberately undecided. GROWTH is measured as the gain over the
# trailing window EXPRESSED IN ELASTIC DISPLACEMENTS, so the two signals share one
# yardstick; 2% of the elastic displacement over the final quarter of a solve is far
# above the noise of a frozen state and far below what any runaway produces.
_HYBRID_SAMPLE_EVERY = 10      # iterations between max|u| samples (cost ~0: the
                               # value is already computed by the CHECON test)
_HYBRID_WINDOW_FRAC = 0.25     # trailing fraction of the history used for growth
_HYBRID_MIN_SAMPLES = 8        # below this there is no trend to read
_HYBRID_U_STUCK_MAX = 1.25     # max|u| / max|u|_elastic at or under this = elastic scale
_HYBRID_U_FAIL_MIN = 1.5       # ... at or over this = beyond elastic scale
_HYBRID_GROWTH_MIN = 0.02      # elastic displacements gained over the trailing window
# Both signals are RATIOS to the elastic displacement scale, so that scale must be a
# real length before either can be read. A yardstick far below the model's own size
# is not a small measurement, it is no measurement: dividing by it turns rounding
# noise into arbitrarily large u_ratio and growth, and every non-converged trial
# scores FAILED on evidence that means nothing. The floor below (a fraction of the
# mesh height) is the size at which a displacement stops being a length and becomes
# noise; under it the classifier declines to rule.
_HYBRID_U_SCALE_FLOOR_FRAC = 1e-6

# Iterations without a >1% improvement on the best out-of-balance value seen after
# which the residual is called PLATEAUED. This is a reporting threshold only: a
# plateau is recorded in the result and the solve keeps running (see the no-progress
# watch inside solve_fem for why it may not decide a verdict).
_NO_PROGRESS_WINDOW = 1500

# === Budget extension =========================================================
# Reaching `max_iterations` is not a verdict either. A trial whose out-of-balance
# is still TRENDING DOWN when the budget runs out has not failed; it has run out of
# budget, and the number of iterations a viscoplastic solve needs grows with mesh
# refinement and with proximity to the critical F. So the budget is EXTENDED, one
# chunk at a time, for as long as the residual keeps falling, up to a hard ceiling
# (`max_iterations_ceiling`).
#
# The trend is read from the mean of the last window against the mean of the window
# before it, rather than from a single iterate: the residual oscillates on the yield
# surface (see oob_window) and creeps non-monotonically, so consecutive values say
# nothing. Requiring the same 1% the no-progress watch uses keeps one definition of
# "meaningful improvement" in the file.
_OOB_TREND_WINDOW = 500        # iterations averaged in each of the two windows
_OOB_TREND_MIN = 0.01          # the later window must be this much lower (1%)
# A steady decay is not the only way a solve is worth more iterations. Measured on
# the reinforced slope at F = 1.25 (tri6, 1 ft): the residual falls to 2e-3 by
# iteration 9,000, sits there, then RISES eighty-fold through a burst of plastic
# redistribution around iteration 14,000 before coming back down and reaching
# equilibrium at 16,242 — while max|u| moves from 0.1101 to 0.1139 ft. The slope is
# standing still the whole time; only the residual is thrashing. At iteration 12,000
# that solve is mid-excursion, so a trend test alone reads RISING and stops it 4,000
# iterations short of its answer.
#
# The second signal is therefore the DISPLACEMENT field, read through the same
# classifier the failure criterion uses: a trial whose displacement is not growing
# and whose evidence the classifier cannot rule on (AMBIGUOUS) has not shown itself
# to be failing, and gets more budget. A trial that IS growing is failing and gets
# none — which is also where the iterations are saved. A STABLE_STUCK trial is
# excluded deliberately: the classifier can already rule on it, and under the hybrid
# criterion it counts as standing, so spending the ceiling on it would buy nothing.


# === Early failure ============================================================
# The mirror image of the budget extension. A trial that is running away does not
# need its whole budget to prove it: the evidence is already in the two series the
# loop samples every `_HYBRID_SAMPLE_EVERY` iterations, and once it is unambiguous
# there is nothing left to learn from the remaining iterations.
#
# The thresholds were calibrated pass-by-pass on 45 strength-reduction trials — 24
# that reached equilibrium and 21 that did not — drawn from five complete bisection
# walks on two models (an embankment at two mesh sizes, a reinforced slope under two
# bar rules, and one of those at a deliberately small budget). NEAR the critical
# strength the two sides of the bracket are not separable: a trial that reaches
# equilibrium can grow past five times its elastic response, sit with a flat
# residual for thousands of iterations, and then converge. The intent of the levels
# below is therefore to sit OUTSIDE that region, catching only the gross runaways —
# about half of the failing trials — and leaving the near-critical ones to spend
# their budget as before.
#
# CORRECTED 2026-08-31. This comment used to state the margin as a fact about the
# rule: that the levels sit where no trial that reached equilibrium was observed to
# come near, "the closest came within a factor of 3.0 of the first and 1.59 of the
# second". That is a property of the calibration set, not of the rule, and it does
# not survive contact with a model outside it.
#
# Measured this session on Griffiths & Lane 1 at quad9, mesh size 3.5, force
# tolerance 1e-3. The trial at F = 1.3875 is closed HERE as 'runaway', at 23,861
# iterations, with max|u| at 8.0002 times its elastic response — the level, to four
# figures. Re-run with early_failure=False and one 400,000-iteration budget, the
# same trial on the same mesh reaches equilibrium: 72,580 iterations, out-of-balance
# 9.998e-4 against the 1e-3 tolerance, max|u| = 0.5689 m on a 50 m model, which is
# 1.14% of its height and 10.03 times the elastic response.
#
# So `u_max` = 8.0 does not sit outside the equilibrating trials' range on this
# model; it sits INSIDE the displacement path of one of them, at 80% of where that
# trial finally settles. The margin is not a factor of 3, it is negative, and the
# rule has a demonstrated false-positive mode: it can close a trial that would have
# stood. The two neighbours are closed the same way (F = 1.390625 at 15,231
# iterations, F = 1.396875 at 9,241, both at 8.00 times elastic).
#
# The rule is left EXACTLY as calibrated. Every locked and published factor of
# safety in this repository is defined by it, and moving it would move them; that is
# the owner's call, not a comment's. What changes here is only the claim: this is a
# threshold with a known false-positive mode near the critical strength, not a
# threshold with a proven margin.
#
# Both are measured in the trial's OWN elastic displacement, so the rule carries to
# any model and any mesh without retuning.
_EARLY_FAIL_WINDOW = 2000      # iterations in the trailing window read by (b)
_EARLY_FAIL_WARMUP = 500       # (b) may not fire before this iteration
_EARLY_FAIL_GAIN = 1.0         # (b): elastic displacements gained over that window
_EARLY_FAIL_U_MAX = 8.0        # (f): elastic displacements, the absolute runaway level


def _early_failure(disp_hist, oob_hist, u_elastic_scale,
                   sample_every=_HYBRID_SAMPLE_EVERY,
                   window=_EARLY_FAIL_WINDOW, gain=_EARLY_FAIL_GAIN,
                   u_max=_EARLY_FAIL_U_MAX, grow=_HYBRID_GROWTH_MIN,
                   min_fall=_OOB_TREND_MIN, warmup=_EARLY_FAIL_WARMUP):
    """Has this trial already proved that it is failing?

    Read from the same two sampled series the budget extension uses, at the latest
    sample. Either signal is sufficient:

      * ``'stalled_residual'`` — the residual has NOT fallen by ``min_fall`` over the
        last ``window`` iterations (the mean of the window's second half against the
        mean of its first, i.e. ``_oob_still_falling`` negated) AND the field gained
        at least ``gain`` elastic displacements over that same window. A residual
        that has stopped coming down while the slope moves by a whole elastic
        response is not a solve in difficulty; it is a slope in motion.
      * ``'runaway'`` — max|u| has passed ``u_max`` elastic displacements and is
        still gaining (at least ``grow`` of one over the last doubling of the
        iteration count). No trend shape is asked for, only a level, which is the one
        thing a trend test cannot supply.

    The first is gated by ``warmup``; the second is not, because it is an absolute
    level — eight times the elastic response — that a solve still warming up cannot
    reach, and gating it would only delay the cheapest catches (the grossest
    runaways declare themselves within a few hundred iterations).

    Returns the signal's name, or None while neither is satisfied.
    """
    if not u_elastic_scale or u_elastic_scale <= 0.0:
        return None
    i = len(disp_hist) - 1
    if i < 4:
        return None                      # no doubling to read yet
    d_now = float(disp_hist[i])
    # Still gaining: the last doubling of the iteration count added displacement.
    gaining = (d_now - float(disp_hist[i // 2])) >= grow * u_elastic_scale

    # (f) absolute runaway
    if d_now >= u_max * u_elastic_scale and gaining:
        return 'runaway'

    # (b) residual stalled while the field keeps moving
    k = max(2, int(window // max(1, sample_every)))     # samples in the window
    if i < k or i < (warmup // max(1, sample_every)):
        return None
    h = k // 2
    prev = sum(oob_hist[i - k + 1:i - h + 1]) / float(k - h)
    last = sum(oob_hist[i - h + 1:i + 1]) / float(h)
    if prev > 0.0 and last < (1.0 - min_fall) * prev:
        return None                      # still falling -> not a failure
    if (d_now - float(disp_hist[i - k])) < gain * u_elastic_scale:
        return None                      # not moving fast enough to call it
    return 'stalled_residual'


def _still_progressing(oob_hist, disp_hist, u_elastic_scale, mesh_height,
                       window=_OOB_TREND_WINDOW):
    """Is this solve worth more iterations than its budget allowed?

    Two signals, either of which counts:

      * the residual is TRENDING DOWN — the mean over the last ``window`` iterations
        is at least ``_OOB_TREND_MIN`` below the mean over the ``window`` before it.
        This is the ordinary case: a solve part-way down its decay.
      * the displacement field is STANDING STILL on evidence the failure classifier
        cannot rule on — verdict AMBIGUOUS with no trailing growth. This is the
        excursion case: the residual is off on a rise while the slope does not move.

    Neither signal is present in a trial whose displacement is GROWING while its
    residual sits flat: that is a slope failing, and it stops at its budget exactly as
    it always has. A STABLE_STUCK trial is left alone too — the classifier can
    already rule on it, and under the hybrid criterion it counts as standing, so
    spending the ceiling on it would buy nothing.
    """
    if _oob_still_falling(oob_hist, window=window):
        return True
    verdict, _, growth = classify_nonconvergence(
        disp_hist, u_elastic_scale, 'iteration_cap', model_height=mesh_height)
    return (verdict == 'AMBIGUOUS' and growth is not None
            and growth <= _HYBRID_GROWTH_MIN)


def _oob_still_falling(oob_hist, sample_every=_HYBRID_SAMPLE_EVERY,
                       window=_OOB_TREND_WINDOW, min_fall=_OOB_TREND_MIN):
    """Is the out-of-balance residual still trending down?

    Compares the mean of the last ``window`` iterations against the mean of the
    ``window`` before it, both read from ``oob_hist`` (sampled every
    ``sample_every`` iterations). Returns False when there is not yet enough
    history to judge, which is the conservative answer: no history, no extension.
    """
    k = max(2, int(window // max(1, sample_every)))
    if len(oob_hist) < 2 * k:
        return False
    last = sum(oob_hist[-k:]) / k
    prev = sum(oob_hist[-2 * k:-k]) / k
    if not (prev > 0.0):
        return False
    return last < (1.0 - min_fall) * prev


def classify_nonconvergence(disp_hist, u_elastic_scale, exit_reason,
                            sample_every=_HYBRID_SAMPLE_EVERY, model_height=None):
    """Classify a NON-CONVERGED viscoplastic trial from its displacement history.

    This is the hybrid failure criterion's discriminator. It answers "did this trial
    fail, or did the solver merely fail to settle it?" from two signals measured
    against the trial's own elastic response:

      * ``u_ratio``  = max|u| / max|u|_elastic at the end of the solve — is the
        displacement field beyond elastic scale at all?
      * ``growth``   = (max|u|_end - max|u|_window_start) / max|u|_elastic over the
        trailing ``_HYBRID_WINDOW_FRAC`` of the history — is it still moving?

    Verdicts:

      ``'FAILED'``        both signals present (or the displacement cap was hit) —
                          the slope is failing; the bisection treats it as failed,
                          exactly as the legacy criterion does.
      ``'STABLE_STUCK'``  both signals absent — frozen at elastic scale. The slope
                          is standing still; only the residual test is unsatisfied.
                          The bisection treats it as NOT failed.
      ``'AMBIGUOUS'``     one signal without the other, or too little history. The
                          legacy verdict stands (failed), flagged so it is never
                          silent.

    Requiring BOTH signals in each direction is what keeps the hybrid conservative:
    it overrides the legacy verdict only where the evidence is unambiguous, and every
    override is recorded in the result dict.

    Parameters:
        disp_hist (list of float): max|u| sampled every ``sample_every`` iterations.
        u_elastic_scale (float): max|u| of the purely elastic solution for this trial.
        exit_reason (str): 'iteration_cap', 'inconclusive', 'disp_limit' or
            'diverging' ('no_progress' is still accepted and read as
            'iteration_cap'; solve_fem no longer ends a solve on a no-progress
            plateau). 'disp_limit' and 'diverging' are evidence in their own right
            and return FAILED.
        sample_every (int): sampling stride (documentation only; the window is a
            fraction of the sample count, so the stride does not enter the maths).
        model_height (float or None): mesh height, used only to floor the elastic
            scale at ``_HYBRID_U_SCALE_FLOOR_FRAC`` x height. An elastic scale below
            that floor is noise, not a length, and both signals are ratios to it, so
            the verdict is AMBIGUOUS rather than a FAILED read off a meaningless
            denominator. None disables the floor (the ``<= 0`` guard still applies).

    Returns:
        tuple: ``(verdict, u_ratio, growth)`` — floats are ``None`` when there is no
        elastic scale to divide by.
    """
    n = len(disp_hist)
    scale_floor = (_HYBRID_U_SCALE_FLOOR_FRAC * float(model_height)
                   if model_height else 0.0)
    if (not u_elastic_scale or u_elastic_scale <= 0.0 or n == 0
            or float(u_elastic_scale) < scale_floor):
        # 'disp_limit' is the one verdict that survives a missing yardstick: it is an
        # ABSOLUTE displacement budget (a fraction of mesh height), so it is evidence
        # in its own right and does not divide by the elastic scale. ('diverging' is
        # listed for completeness only — the early-failure rule reads ratios to the
        # elastic scale and cannot fire without one.)
        return ('FAILED' if exit_reason in ('disp_limit', 'diverging')
                else 'AMBIGUOUS'), None, None

    u_ratio = float(disp_hist[-1]) / float(u_elastic_scale)

    # The displacement cap is itself corroborating evidence: the trial did not just
    # fail to settle, it ran away far enough to trip a physical displacement budget.
    if exit_reason == 'disp_limit':
        return 'FAILED', u_ratio, None

    if n < _HYBRID_MIN_SAMPLES:
        return ('FAILED' if exit_reason == 'diverging' else 'AMBIGUOUS'), u_ratio, None

    k = max(2, int(round(n * _HYBRID_WINDOW_FRAC)))
    growth = (float(disp_hist[-1]) - float(disp_hist[n - k])) / float(u_elastic_scale)
    growing = growth > _HYBRID_GROWTH_MIN

    # The early-failure rule (see _early_failure) is a failure verdict already made,
    # on thresholds no trial that reaches equilibrium has been observed to reach.
    # It is not re-litigated here on a truncated history — reading the classifier at
    # a fraction of a budget is exactly the false-FAILED mechanism the no-progress
    # exit used to produce.
    if exit_reason == 'diverging':
        return 'FAILED', u_ratio, growth

    if u_ratio >= _HYBRID_U_FAIL_MIN and growing:
        return 'FAILED', u_ratio, growth
    if u_ratio <= _HYBRID_U_STUCK_MAX and not growing:
        return 'STABLE_STUCK', u_ratio, growth
    return 'AMBIGUOUS', u_ratio, growth


def solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=12000, tolerance=1e-3,
              max_disp_factor=0.1, tension_cutoff=False, dt_scale=1.0,
              force_tol=1e-3, oob_window=10,
              early_exit=True, progress_callback=None, min_slip_depth=None,
              ssr_exclude_mask=None, tension_cap_by_elem=None, tension_srf=None,
              elastic_mask=None,
              suction_phi_b=None, suction_cap=None, _prepared=None,
              fast_kernel='auto', failure_criterion="hybrid", k0=None,
              max_iterations_ceiling=50000, early_failure=True, _init_state=None,
              _softened_seed=None, fem_solver=None, _nr_export=None,
              _nr_rescue_rungs=None, _nr_seed_first=False):
    """
    Solve FEM using the Griffiths & Lane (1999) viscoplastic algorithm.

    Implements the algorithm from the 1999 Geotechnique paper:
    - 8-node quadrilateral elements with reduced integration (4 Gauss points)
    - Viscoplastic stress redistribution with accumulated plastic strains
    - Pre-factored elastic stiffness matrix for efficiency
    - No damping (stability from dt parameter)
    - Direct solve each iteration (not residual-based)

    A trial is CONVERGED (the slope is stable at this F) when both hold:

      1. Displacements have stopped changing — Smith & Griffiths' CHECON test,
         max|du| / max|u| < `tolerance`.
      2. Force equilibrium has been reached — the maximum over all nodes of the
         nodal out-of-balance force, each normalized by that node's own
         gravitational body force, is below `force_tol`. This is the criterion of
         Dawson, Roth & Drescher (Geotechnique 49(6), 1999). Its locality is what
         makes it trustworthy: inert material added to the mesh sits in equilibrium
         and contributes ~0 to the maximum, so padding a model with extra foundation
         or runout cannot dilute the measure. A GLOBAL norm ratio can be diluted
         exactly that way, which is what an earlier version of this code did, and it
         inflated the factor of safety on deeply-founded models.

    Failure to satisfy both within `max_iterations` is failure of the slope, which
    is Griffiths & Lane's non-convergence criterion.

    With `failure_criterion='hybrid'` (the default since 2026-07-26) that last
    sentence is qualified: a non-converged trial is additionally required to show
    displacement evidence of failure before it counts as failed. See
    `classify_nonconvergence` and the `verdict` / `stable` keys of the returned
    dictionary. Pass `failure_criterion='non_convergence'` for the legacy verdict.

    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        F (float): Shear strength reduction factor (c/F, tan(phi)/F)
        debug_level (int): 0=silent, 1=summary, 2=per-iteration
        max_iterations (int): Viscoplastic iteration budget per trial (default
            12000). It is a budget, not a ceiling: a trial that reaches it with the
            out-of-balance residual still trending down is EXTENDED by another
            budget's worth, repeatedly, while the trend holds.
        max_iterations_ceiling (int): Hard stop on that extension (default 50000).
            A trial that reaches the ceiling while still improving stops with
            exit_reason 'inconclusive' - neither converged nor failed - and
            solve_ssrm reports it as the bracket's upper uncertainty rather than
            counting it as a failure.
        early_failure (bool): End a trial as soon as its movement is unambiguously
            running away, rather than spending the rest of its budget proving it
            (default True). exit_reason 'diverging', which the bisection reads as
            failed. The thresholds sit far outside the range any trial that reaches
            equilibrium has been observed to occupy, so no verdict moves; see
            `_early_failure`. Off for the at-failure capture solve, whose whole
            purpose is to let the mechanism develop.
        tolerance (float): Convergence tolerance ||du|| / ||u|| (default 1e-3).
            Normalized by the current displacement (Smith & Griffiths CHECON-style),
            so steady benign viscoplastic creep is accepted as converged; false
            convergence at true failure is guarded by max_disp_factor below.
        max_disp_factor (float): Displacement limit as fraction of mesh height (default 0.1).
            If max VP displacement (total - elastic) exceeds this fraction of the mesh height,
            the slope is declared as failed regardless of convergence. This prevents false
            convergence when large displacements make the relative change appear small.
            Set to None to disable. NOTE this yardstick is the height of the MESH, not of
            the SLOPE, so it grows when a model is given a deeper foundation; prefer
            force_tol, which has no such dependence. The default SSRM path disables it.
        force_tol (float): Force-equilibrium tolerance (default 1e-3, the value used by
            Dawson, Roth & Drescher 1999). The state is in equilibrium when the maximum
            over all nodes of |nodal out-of-balance force| / |nodal body force| falls below
            this. The denominator is a LUMPED tributary weight (sum over the adjacent
            elements of gamma_e * A_e / n_e), NOT the consistent nodal gravity load — the
            consistent load is exactly zero at a tri6 corner, which makes the ratio there
            meaningless. Being a force over a force at the same node, the ratio is
            dimensionless and independent of the unit system and — the property this test
            exists for — of the size of the domain and of the yielding zone. It is NOT
            independent of element size: it goes roughly as 1/h, so a coarser mesh narrows
            the margin to this tolerance.
        oob_window (int): Number of iterations the plastic-flow increment is averaged over
            before the per-node maximum is taken (default 10). Must be >= 2. A one-iteration
            increment does NOT decay on a settled slope — Gauss points resting on the yield
            surface flip flow direction every iteration, giving a period-2 limit cycle whose
            amplitude scales with dt but never vanishes — so with oob_window=1 a stable slope
            can never converge. Averaging cancels that mode exactly and leaves genuine
            plastic drift untouched. The verdict is insensitive to the width (10, 50 and 200
            agree), so this is not a tuning parameter.
        min_slip_depth (float or None): Optional surficial-failure filter (default None
            = off). When set, nodes shallower than this depth below the ground surface are
            excluded from the out-of-balance maximum, so a shallow cohesionless "skin"
            (FS = tan phi / tan beta, depth-independent for c=0) cannot on its own declare
            the slope failing. A genuine deep-seated mechanism still trips the criterion
            through its deep nodes. This is the SSRM analogue of an LEM minimum-slip-depth
            search filter (Slide2's "minimum depth"); the search-side twin lives in
            search.py. Off by default, so the reported FS is the true global minimum
            (surficial skin included) unless the caller opts in.
        tension_cutoff (bool): If True (default False), applies a RANKINE tension
            cutoff at T = 0 everywhere: the major principal (most-tensile) stress is
            capped at zero via a second viscoplastic yield surface F_t = sigma_1 - 0
            (damped, ~10% per iteration; see tension_cap_by_elem for the mechanism).
            No principal tensile stress is permitted anywhere. The psi=0
            Mohr-Coulomb flow is purely deviatoric and cannot return near-apex
            tensile states to the envelope; this surface handles them. Equivalent in
            spirit to the tension cutoff used by commercial SSRM codes; Griffiths &
            Lane (1999) include no tension treatment. This is the T = 0 special case
            of the per-element cap below.
        dt_scale (float): Multiplier on the viscoplastic pseudo-timestep
            (research/diagnostic knob; default 1.0).
        progress_callback (callable or None): If given, called (throttled, ~every
            10 iterations) as ``progress_callback(frac, label)`` where ``frac`` is
            ``(iteration + 1) / max_iterations`` in [0, 1] and ``label`` is a short
            status string. Lets a caller (e.g. SSRM, a GUI) show progress *within*
            a single viscoplastic solve rather than only between solves. ``frac`` is
            a pessimistic estimate — a trial that converges or exits early snaps to
            its step boundary before reaching 1.0. Never allowed to break the solve.
        ssr_exclude_mask (array of bool or None): Per-element mask (length n_elements)
            of elements to hold at FULL strength during reduction — their c and
            tan(phi) are divided by 1.0 rather than F. This is the SSR-exclusion
            mechanism; solve_ssrm builds it from material names (``ssr_exclude``), an
            explicit search-area polygon (``ssr_zone``) and the file's own v20 SSR
            zone overlay rows, unioned. None (default) = the file's ``ssr_zones``
            overlays if it carries any, else every element reduced by F
            (bit-identical to the un-excluded path).
        tension_cap_by_elem (array of float or None): Per-element tensile-strength
            cap T (length n_elements, problem stress units) for the RANKINE tension
            cutoff. A finite entry T caps the major (most-tensile) principal stress
            at T through a second viscoplastic yield surface F_t = sigma_1 - T
            (associated flow, damped by dt_r; summed with the Mohr-Coulomb shear
            flow at the corner by Koiter's rule). An inf (or NaN) entry turns the
            cutoff off for that element. This is the SAME mechanism the global
            ``tension_cutoff`` flag uses — that flag is simply T = 0 for every
            element (no principal tension allowed) — so there is ONE mechanism, not
            two. A per-element cap overrides the global value element by element. It
            is the principal-stress (Rankine) tensile cutoff of RS2/PLAXIS/FLAC
            (solve_ssrm builds this array from a material-name -> T dict). None
            (default) = cutoff governed solely by ``tension_cutoff`` (bit-identical
            to the pre-existing path when that is also False).
        tension_srf (bool): If True (DEFAULT), the tensile cap T is divided by the
            trial strength-reduction factor each solve, exactly as c and tan(phi)
            are, so the factor of safety is the factor by which the WHOLE strength
            envelope — shear and tensile — is reduced. This is RS2's
            ``tensilestrength_SRF=1``, the setting its entire published
            verification set uses (Plaxis reduces the cap the same way).
            NO-OP WITHOUT A CAP: with no per-element cap and no global
            ``tension_cutoff`` there is no T to reduce, so on a model that sets
            neither this flag cannot change a single number — every cap-less run
            (all Griffiths & Lane anchors) is bit-identical either way.
            If False the cap is held at its authored value regardless of F (RS2's
            flag = 0), reproducing the pre-2026-07 static-cap behavior when this
            argument defaulted off. No effect on the T = 0 global cutoff
            (0/F = 0).
        k0 (float or None): At-rest lateral earth pressure coefficient for the
            INITIAL STRESS STATE. None (default) = the historical gravity turn-on:
            the model starts from zero stress and gravity is switched on in one
            step, so the initial lateral stress is whatever plane-strain elasticity
            produces, sigma_h = nu/(1-nu) * sigma_v (K0 ~ 0.43 at nu = 0.3, and it
            is the STIFFNESS, not the soil, that sets it). Real soil is not that
            lightly confined -- normally consolidated sand sits near 1 - sin(phi)
            ~ 0.43, but a compacted fill or an overconsolidated clay runs to 1.0
            and beyond, and every vendor model in the verification corpus is
            authored with Kx = Kz = 1.
            When set, the initial stress at each Gauss point is built from the
            overburden instead:
                sigma'_v = -(soil column weight above the point) + u
                sigma'_h = sigma'_z = k0 * sigma'_v   (in-plane AND out-of-plane)
                tau_xy   = 0
            (tension-positive; u is the model's pore pressure). The state is then
            carried through the classical initial-stress method --
            sigma = sigma_0 + D (B u - evp), with int B^T sigma_0 dV moved to the
            right-hand side -- so the solver still ITERATES TO EQUILIBRIUM under the
            body forces from that starting point. Where sigma_0 already balances
            gravity (level ground with k0 = nu/(1-nu)) the displacement solution is
            ~0; on a slope it is not, and the viscoplastic loop redistributes.
            The overburden is soil only: surface tractions (a reservoir load, a
            distributed load) are NOT part of the in-situ state and enter as
            boundary forces during the iteration. Requires material-zone geometry
            (it integrates a vertical ray through the zones, exactly as the ``ru``
            pore-pressure option does) and raises if the model carries none.
            The compiled Mohr-Coulomb kernel has no slot for an initial stress, so a
            k0 run always takes the NumPy reference path -- the oracle, but slower.
            Off by default, so every locked factor of safety is untouched.
            A SINGLE solve is its own in-situ equilibration: the redistribution of
            the K0 field against the geometry and the response to the applied loads
            are one and the same solve, and its displacements are measured from the
            un-redistributed K0 state. At F = 1 that is exactly what is wanted. At a
            REDUCED strength it is not, because the redistribution then happens
            against weakened soil and its plastic strain is charged to the trial;
            that separation is solve_ssrm's job (see ``_init_state`` and the
            equilibration solve it runs before the bisection), and a direct
            reduced-strength solve_fem call does not perform it.
        _init_state (dict or None): INTERNAL. The equilibrated state of an earlier
            full-strength k0 solve (its ``_k0_state``: the displacement field and the
            accumulated viscoplastic strain). This solve then STARTS from that state
            — same stress field, no in-situ redistribution to repeat — and measures
            displacement from it: reported displacements, the CHECON convergence
            ratio and the hybrid criterion's history are all relative to the
            equilibrated configuration, while stresses and structural forces remain
            functions of the absolute displacement. Requires k0 and the same prepared
            model the state was produced with.
        elastic_mask (array of bool or None): Per-element mask (length n_elements)
            marking elements that are PURE LINEAR ELASTIC — held out of the
            plastic-correction loop entirely. A True element accumulates no
            viscoplastic strain: it never yields, is never tension-relaxed, and is
            never flagged in plastic_elements. Its converged stress is the linear
            elastic D*B*u, even where that state lies outside the Mohr-Coulomb
            envelope. This mirrors RS2's "Plasticity Specifications: None" placed
            materials. It is DISTINCT from ssr_exclude_mask: an SSR-excluded element
            keeps FULL (un-reduced) strength but STILL yields once its stress
            reaches that full envelope, whereas an elastic_mask element has no
            strength surface at all and cannot yield under any stress. The two masks
            are independent and compose freely. None (default) = the file's own
            elastic set — ``option = elastic`` materials plus any v20 "SSR elastic"
            (-3) overlay row — else no elastic elements (bit-identical to the
            pre-existing path).
        suction_phi_b (dict or None): OPT-IN matric-suction strength, {material
            name: phi_b degrees}. Turns the signed (un-clamped) pore pressure's
            negative part s = max(0, -u) into an apparent cohesion s*tan(phi_b),
            reduced by the trial F, added to c' in the MC yield; the effective-
            normal u stays clamped exactly as today. None (default) auto-wires from
            the v17 template (fem_data['phi_b_by_mat']); an explicit dict overrides,
            {} forces off. Off => bit-identical to the pre-suction solver. See
            _resolve_suction_by_elem and the suction note after the reduction block.
        suction_cap (float, dict, or None): Cap on the credited suction s before it
            becomes apparent cohesion (scalar or {name: cap}); None auto-wires from
            fem_data['s_cap_by_mat'] (inf = uncapped). Ignored when phi_b is 0.
        _prepared (dict or None): INTERNAL. A prepared-model dict from
            _prepare_fem_model carrying the strength-reduction-factor-INDEPENDENT setup
            (K factorization, geometry precompute, pore-pressure fields, Dawson g_node
            normalization). solve_ssrm builds it ONCE and threads it through every
            bisection trial so they share one factorization instead of rebuilding it
            ~10 times. None (default) = build it here, so a standalone solve_fem is
            bit-identical to the pre-cache path. It holds no F-dependent or per-solve
            state, so a reused prepared model cannot serve a stale strength or geometry.
        fast_kernel ('auto' | bool): Compiled (Cython) constitutive kernel for the
            Mohr-Coulomb Step-6 Gauss-point update.
              * 'auto' (DEFAULT since 2026-07-26) — use the compiled
                xslope._fem_kernel when it imports, the NumPy reference when it
                does not. Rationale: the installer builds compile the kernel, so
                installed users get the fast path without knowing it exists, while
                pip users (the wheel is pure-Python) silently get the reference.
                Both give the SAME answers — certified corpus-wide on 2026-07-26
                over all 103 FEM benchmarks (fast == reference on every row that
                solves). 'auto' never warns and never fails for a missing module.
              * True — REQUIRE the compiled kernel; if it is not built, warn and
                transparently fall back to NumPy (behavior unchanged).
              * False — never use it; the pure-NumPy reference path.
            ORACLE DOCTRINE (unchanged by the 'auto' default): every locked and
            published factor of safety in XSLOPE is DEFINED by the NumPy reference
            path, permanently; the compiled kernel is an optimization that must
            reproduce it. The suite therefore pins its kernels explicitly rather
            than inheriting this default — fast-first with a reference fallback on
            a lock miss, and --reference-only forced to fast_kernel=False — and
            benchmarks/kernel_xcheck.py remains the required divergence fence.
            USAGE DOCTRINE: the compiled kernel is never the definition of a
            result — it is an optimization that must reproduce the reference. Do
            not use it to define or re-record a lock; the NumPy reference alone
            does that, forever. Its former carve-out ("never in Studio, never on
            by default") was predicated on _fem_kernel.pyx's KNOWN LIMIT: RS2-62c,
            a thin-soft-band case where floating-point re-association flipped a
            bisection verdict. That case diverged by 0.58 while its model was
            MISSING the vendor's tensile-strength caps; with the caps restored the
            verdicts sit on real displacement runaway and the two kernels agree
            (resolved 2026-07-26, then certified over the whole 103-benchmark FEM
            corpus). Divergence is a smell of a knife-edge MODEL, not merely of
            kernel drift. On that evidence the default is 'auto' everywhere,
            Studio included. The standing guard is benchmarks/kernel_xcheck.py,
            which solves cases both ways and fails on any divergence; it MUST NOT
            be removed while 'auto' is the default.
        early_exit (bool): Watch the out-of-balance residual for a no-progress
            plateau and report it (`plateau_iteration`, `plateau_ratio`). The plateau
            is an observation only - it does not end the solve, and a trial always
            runs to convergence, its iteration cap, or the displacement cap. It used
            to end the solve and report the trial as failed, which cut off trials that
            were still converging and biased the factor of safety low; the residual is
            not monotone and the iterations a trial needs grow with mesh refinement,
            so a fixed window cannot decide the outcome (see the watch in the
            iteration loop).
        failure_criterion (str): 'hybrid' (DEFAULT since 2026-07-26) or
            'non_convergence' (the legacy verdict, still fully supported). Under
            'hybrid', a trial that does not reach equilibrium is put through
            classify_nonconvergence() before it is called failed: a state frozen at
            elastic displacement scale is reported STABLE_STUCK and the 'stable'
            key comes back True, so a caller (solve_ssrm's bisection) does not
            treat a numerically stuck solve as a failed slope. The verdict
            metadata is returned on BOTH settings; only 'stable' differs, and only
            for a STABLE_STUCK trial.
            WHY IT IS THE DEFAULT: a corpus-wide A/B over all 103 FEM benchmarks
            (2026-07-26) moved 4 rows and NO row downward. Three were small upward
            shifts where the legacy criterion had truncated the bracket early on a
            stuck-but-standing trial (a truncation bias, not a physical failure),
            and one was an outright rescue — RS2-48, a T=0 reinforced fill the
            legacy machinery cannot bracket at all because its trials are
            stationary rather than collapsing. The remaining 99 rows are
            bit-identical between the two criteria.

    Returns:
        dict: Solution dictionary with keys:
            - converged (bool): Whether iterations converged (true equilibrium)
            - stable (bool): Whether the CALLER should treat the slope as standing at
              this F. Equals `converged` on every criterion except 'hybrid', where it
              is also True for a STABLE_STUCK verdict.
            - verdict (str): 'CONVERGED' | 'FAILED' | 'STABLE_STUCK' | 'AMBIGUOUS'
            - u_ratio (float or None): max|u| / max|u|_elastic at the end of the solve
            - u_growth (float or None): elastic displacements gained over the trailing
              window of the iteration history (the growth signal)
            - exit_reason (str): 'converged' | 'iteration_cap' | 'inconclusive' |
              'disp_limit' | 'diverging' - why this solve stopped
            - diverging_iteration (int or None): iteration at which the early-failure
              rule fired, and diverging_signal (str or None) which of its two tests
              fired ('runaway' or 'stalled_residual'); None on any other exit
            - plateau_iteration (int or None): iteration at which the out-of-balance
              stopped improving by >1% for `_NO_PROGRESS_WINDOW` iterations, or None
              if it never stalled. A plateau does NOT stop the solve; it is reported
              so a trial that ran to its cap can say the residual had gone flat.
            - plateau_ratio (float or None): the out-of-balance ratio it stalled at
            - iterations (int): Number of iterations used
            - displacements (ndarray): Nodal displacement vector
            - stresses (ndarray): Element stress array (n_elements, 4) [sig_x, sig_y, tau_xy, sig_vm] compression-positive
            - strains (ndarray): Element strain array (n_elements, 4) [eps_x, eps_y, gamma_xy, max_shear]
            - plastic_elements (ndarray): Boolean array of yielded elements
            - F (float): Applied strength reduction factor
            - forces_1d (ndarray): Final axial forces in 1D truss elements (empty if none)
            - failed_1d_elements (ndarray): Boolean array of failed 1D elements (empty if none)
    """

    if debug_level >= 1:
        print(f"=== Griffiths & Lane Viscoplastic FEM (F={F:.3f}) ===")

    # Extract data
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    element_materials = fem_data["element_materials"]
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]
    fixed_nodes = fem_data.get("fixed_nodes", set())
    roller_x_nodes = fem_data.get("roller_x_nodes", set())
    roller_y_nodes = fem_data.get("roller_y_nodes", set())

    # Template-carried defaults (v16, tension_cap_by_elem / elastic_mask), the suction
    # resolution, the Rankine cap base, and every other strength-reduction-FACTOR-
    # INDEPENDENT quantity are resolved in _prepare_fem_model below (which honors the
    # same explicit-kwarg-wins semantics). solve_ssrm passes a prepared model so the
    # ~10 bisection trials share one K factorization / geometry precompute; a
    # standalone solve_fem builds its own here, bit-identical to before.

    # Material properties
    c_by_elem = fem_data.get("c_by_elem", fem_data["c_by_mat"][element_materials - 1])
    phi_by_elem = fem_data.get("phi_by_elem", fem_data["phi_by_mat"][element_materials - 1])
    pow_flag_by_elem = fem_data.get("pow_flag_by_elem")
    if pow_flag_by_elem is None:
        pow_flag_by_elem = np.zeros(len(elements), dtype=bool)
    has_pow = bool(np.any(pow_flag_by_elem))
    hb_flag_by_elem = fem_data.get("hb_flag_by_elem")
    if hb_flag_by_elem is None:
        hb_flag_by_elem = np.zeros(len(elements), dtype=bool)
    has_hb = bool(np.any(hb_flag_by_elem))
    E_by_mat = fem_data["E_by_mat"]
    nu_by_mat = fem_data["nu_by_mat"]
    gamma_by_mat = fem_data["gamma_by_mat"]
    k_seismic = fem_data.get("k_seismic", 0.0)

    n_nodes = len(nodes)
    n_elements = len(elements)

    # DOF offset map: pile nodes get 3 DOFs, others get 2
    dof_offset = fem_data.get("dof_offset", None)
    if dof_offset is not None:
        n_dof = int(dof_offset[n_nodes])
    else:
        n_dof = 2 * n_nodes

    # ---- v19 file-carried run options: kwarg > file > engine default ----
    # solve_ssrm resolves these itself and always passes concrete values down, so
    # this only fires on a standalone solve_fem. None means UNSPECIFIED at every
    # level, which is why a pre-v19 file lands on exactly the historical defaults.
    if k0 is None:
        k0 = fem_data.get("k0")
    if tension_srf is None:
        tension_srf = fem_data.get("tension_srf")
    if tension_srf is None:
        tension_srf = True                      # engine default

    # ---- Strength-reduction-factor-INDEPENDENT setup (built once, reused) ----
    # solve_ssrm passes a prepared model in via _prepared so all trials share the K
    # factorization and geometry precompute; a standalone call builds its own here.
    if _prepared is not None:
        prep = _prepared
    else:
        prep = _prepare_fem_model(
            fem_data, dt_scale=dt_scale, suction_phi_b=suction_phi_b,
            suction_cap=suction_cap, elastic_mask=elastic_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_cutoff=tension_cutoff,
            min_slip_depth=min_slip_depth, k0=k0, debug_level=debug_level)

    K_factor = prep["K_factor"]
    F_gravity = prep["F_gravity"]
    free_dofs = prep["free_dofs"]
    n_free = prep["n_free"]
    dt = prep["dt"]
    elem_gp_data = prep["elem_gp_data"]
    u_gp = prep["u_gp"]
    u_gp_signed = prep["u_gp_signed"]
    pp_option = prep["pp_option"]
    n_total_gp = prep["n_total_gp"]
    mesh_height = prep["mesh_height"]
    # K0 initial stress. The prepared model decides whether the overburden exists
    # (solve_ssrm passes its own k0 through to _prepare_fem_model), so a trial that
    # inherits a prepared model built WITHOUT k0 cannot silently turn it on.
    sv0_gp = prep.get("sv0_gp") if k0 is not None else None
    if k0 is not None and sv0_gp is None:
        raise ValueError(
            "k0 was given to solve_fem but the prepared model carries no overburden "
            "field. Pass k0 to _prepare_fem_model / solve_ssrm as well.")
    # ---- Carried-in equilibrated state (internal; see _init_state) ----
    # The state is the displacement field and the accumulated viscoplastic strain at
    # the end of an earlier full-strength solve. Together with the K0 initial stress
    # (rebuilt below exactly as always) they define the equilibrated stress field
    #     sigma = sigma_0 + D (B u_eq - evp_eq),
    # so starting this solve from them starts it from that field. The solve itself
    # stays in ABSOLUTE displacements — every internal force, including the axial
    # force in a reinforcement bar and the end forces in a pile, is a function of the
    # absolute displacement and must keep seeing it. What the datum reset changes is
    # the MEASUREMENT: displacement is reported, and the convergence and failure
    # criteria are read, relative to the equilibrated state, because the in-situ
    # displacement is an artifact of imposing a stress field that the geometry does
    # not hold in equilibrium — the soil did not travel there, it was always there.
    # The evp is a list of per-Gauss-point-group (G, 4) arrays in the group order the
    # prepared model defines, hence the requirement that the state was produced on
    # this same prepared model.
    _init_evp = None
    _init_u = None
    if _init_state is not None:
        if sv0_gp is None:
            raise ValueError("_init_state was given without k0; an equilibrated "
                             "initial state has no meaning without the K0 "
                             "formulation.")
        _init_evp = _init_state["evp"]
        _init_u = np.asarray(_init_state["u"], dtype=float)
        _sizes = [len(_sg['pairs']) for _sg in prep["gp_groups_static"]]
        if [len(a) for a in _init_evp] != _sizes or _init_u.shape != (n_dof,):
            raise ValueError(
                "_init_state does not match this prepared model's Gauss-point "
                "groups / degrees of freedom; the state must be produced on the "
                "same prepared model.")
    # Displacement datum: the state displacements are measured FROM. Zero (the
    # default) leaves every measurement exactly as it was.
    u_datum = np.zeros(n_dof) if _init_u is None else _init_u
    u_datum_free = u_datum[free_dofs]
    # Per-call (max_disp_factor varies across the SSRM trials vs the capture solve
    # that share a prepared model), so this is recomputed here, never cached.
    if max_disp_factor is not None and mesh_height > 0:
        vp_disp_limit = max_disp_factor * mesh_height
    else:
        vp_disp_limit = None
    if debug_level >= 1 and vp_disp_limit is not None:
        print(f"  VP displacement limit: {vp_disp_limit:.2f} ({max_disp_factor:.0%} of mesh height {mesh_height:.1f})")
    node_dof_x = prep["node_dof_x"]
    node_dof_y = prep["node_dof_y"]
    free_dof_mask = prep["free_dof_mask"]
    node_has_free = prep["node_has_free"]
    g_node_den = prep["g_node_den"]
    _deep_free_mask = prep["deep_free_mask"]
    suction_active = prep["suction_active"]
    suction_tanphib_by_elem = prep["suction_tanphib_by_elem"]
    suction_scap_by_elem = prep["suction_scap_by_elem"]
    elastic_by_elem = prep["elastic_by_elem"]

    # Apply strength reduction (Griffiths & Lane 1999): c_r = c/F, phi_r = atan(tan(phi)/F)
    # Note: Only soil strength (c, phi) is reduced by F. Reinforcement properties
    # (T_allow, T_res, EA/L) are NOT reduced — they are structural capacities.
    #
    # Per-element reduction factor. Elements flagged by ssr_exclude_mask keep FULL
    # strength (F = 1) while the rest are reduced by the trial F — the SSR-exclusion
    # semantics (RS2's per-material Apply_SSR / "SSR Exclusion Area"): a stiff
    # foundation kept at full strength forces the mechanism up into the reducible
    # zones. With no mask every entry equals F, so F_by_elem is exactly the scalar F
    # and the result is bit-identical to the un-excluded path.
    # Direct-call fallback for the file-carried v20 SSR zone overlays, mirroring the
    # tension_cap_by_elem / elastic_mask fallbacks in _prepare_fem_model: a bare
    # solve_fem at F != 1 honors the file's zones. solve_ssrm has already resolved
    # them into an explicit mask and strips 'ssr_zones' from the trial fem_data, so
    # this never fires twice.
    if ssr_exclude_mask is None and fem_data.get("ssr_zones"):
        ssr_exclude_mask, _ = _compose_ssr_zone_masks(fem_data, fem_data["ssr_zones"])

    F_by_elem = np.full(n_elements, float(F))
    if ssr_exclude_mask is not None:
        F_by_elem[np.asarray(ssr_exclude_mask, dtype=bool)] = 1.0
    c_reduced = c_by_elem / F_by_elem
    tan_phi_reduced = np.tan(np.radians(phi_by_elem)) / F_by_elem
    phi_reduced = np.arctan(tan_phi_reduced)  # radians

    # The prepared model's dt ignores phi (it has to -- it is shared across trials),
    # but the viscoplastic stability limit depends on the REDUCED angle, which this
    # trial's F has just fixed. Clamp here, per trial. See _dt_stability_clamp.
    # REVERTED 2026-07-26 after adversarial review. The clamp is NOT applied.
    # Its stated mechanism did not survive testing: at F = 0.40 (where it decided
    # RS2-62c's lock) dt_base/limit = 0.96 — the Cormeau limit is NOT violated, so
    # only the arbitrary 0.9 safety factor made it bite. The worst-case overshoot
    # for nu = 0.3 is bounded at ~1.17, not a hard instability, and the largest
    # overshoot in the sweep (F = 0.15) converges in 512 iterations. Turning it on
    # flipped exactly one verdict on the 0.6 mesh — in the WRONG direction (F = 0.20
    # converged before, failed after). It is also non-conservative where it does
    # bite: the residual is linear in dt against an ABSOLUTE force_tol, so damping
    # dt inflates FS (see docs/fem/overview.md on dt_scale), and at nu = 0.45 (soft
    # clay in elastic_props.py) it cuts dt 37% at F = 1.0. _dt_stability_clamp is
    # kept below, unused, because the FORMULA is correct and the real open question
    # — whether Hoek-Brown per-Gauss-point tangent angles reach the unstable regime
    # — needs it. Do not wire it in without evidence of an actual instability.

    # Per-element tensile-strength cap T for the Rankine tension cutoff (caps the
    # major principal stress; see the tension_cap_by_elem docstring). inf = off.
    # The base (global cutoff -> 0 plus per-material caps) is F-independent and lives
    # in the prepared model; tension_srf then divides it by the trial F here alongside
    # c/F and tan(phi)/F (RS2's tensilestrength_SRF=1). When tension_srf is off the
    # shared base is used as-is (read-only — copied only when it is about to be scaled,
    # so the prepared model is never mutated). elastic_by_elem and the suction arrays
    # likewise come from the prepared model (see the unpack above).
    t_cap_by_elem = prep["t_cap_base"]
    if tension_srf:
        t_cap_by_elem = t_cap_by_elem.copy()
        _red = np.isfinite(t_cap_by_elem) & (t_cap_by_elem > 0.0)
        t_cap_by_elem[_red] = t_cap_by_elem[_red] / F_by_elem[_red]

    # ---- Solver switch (INTERNAL, spike; see SPIKE.md) ----------------------
    # The only place fem_solver is read. 'viscoplastic' — the default, and the
    # definition of every locked and published factor of safety — falls straight
    # through and the rest of this function is untouched. 'newton' hands the trial
    # to the Newton-Raphson driver, which returns the same result dictionary.
    if resolve_fem_solver(fem_solver) == 'newton':
        if _nr_export is not None:
            # The ramp driver (see _ssrm_ramp_newton) continues ONE solve across
            # strength steps, so it needs the working state this call is about to
            # build and a way to re-reduce the strengths in place. Only the ramp
            # passes this; every other caller leaves it None and nothing here runs.
            def _restrength(groups, F_new):
                Fb = np.full(n_elements, float(F_new))
                if ssr_exclude_mask is not None:
                    Fb[np.asarray(ssr_exclude_mask, dtype=bool)] = 1.0
                c_r_new = c_by_elem / Fb
                phi_r_new = np.arctan(np.tan(np.radians(phi_by_elem)) / Fb)
                # The tensile cap is reduced with the shear strength when
                # `tension_srf` is on, exactly as it is in the per-trial code
                # above. Without this the ramp would carry the cap the FOOT was
                # solved at all the way up, so the strength being reduced would
                # not be the whole envelope the setting says it is.
                tc_new = prep["t_cap_base"]
                if tension_srf:
                    tc_new = tc_new.copy()
                    _rr = np.isfinite(tc_new) & (tc_new > 0.0)
                    tc_new[_rr] = tc_new[_rr] / Fb[_rr]
                for grp in groups:
                    e = grp['e_idx']
                    cc = c_r_new[e]
                    if 'elastic' in grp:
                        cc = cc.copy()
                        cc[grp['elastic']] = np.inf
                    grp['c_r'] = cc
                    grp['snph'] = np.sin(phi_r_new[e])
                    grp['csph'] = np.cos(phi_r_new[e])
                    # The matric-suction apparent cohesion is re-derived at the
                    # new strength for the same reason the tensile cap is: it
                    # scales as 1/F, so a ramp carrying the FOOT's credit up the
                    # whole ramp would be reducing only part of the envelope. The
                    # capped suction itself is F-independent and stays cached.
                    _nr_group_apparent_cohesion(grp, 1.0 / Fb)
                    _nr_group_tension_cap(grp, tc_new)
                    # A curved envelope is reduced on its TANGENT, which is
                    # re-derived at every evaluation, so what has to be carried
                    # forward is the per-element F the linearization divides by.
                    # Without this the ramp would linearize the whole way up at
                    # the strength its foot was solved at.
                    _nr_group_restrength_envelope(grp, Fb)
                    grp.pop('_env_seed', None)
            _nr_export['restrength'] = _restrength
        _nr_kw = dict(
            c_reduced=c_reduced, phi_reduced=phi_reduced,
            elastic_by_elem=elastic_by_elem, t_cap_by_elem=t_cap_by_elem,
            finv_by_elem=1.0 / F_by_elem, _nr_env_F=F_by_elem,
            force_tol=force_tol,
            min_slip_depth=min_slip_depth, max_iterations=max_iterations,
            max_disp_factor=max_disp_factor, _nr_export=_nr_export,
            k0=k0, _nr_init_state=_init_state,
            debug_level=debug_level, progress_callback=progress_callback)
        # P3: on a model already known to be carried only from a seed, the cold
        # attempt is skipped and the chain starts at the first rung (see
        # _NR_SEED_MEMORY). The seed-first path is never taken by the ramp, which
        # manages its own history.
        _seed_first = bool(_nr_seed_first) and _nr_export is None
        # P4: the cold attempt's own step control, coarsened only when it is going to
        # be handed off anyway (see _NR_COLD_CHEAP). The correctors below take the
        # shipped control from _nr_kw.
        # ... and only where there IS a hand-off. A ramp foot (_nr_export) has no
        # predictor behind it — the ramp manages its own plastic history — so its
        # cold solve keeps the shipped step control.
        _cold_kw = dict(_nr_kw)
        if _NR_COLD_CHEAP and _nr_export is None:
            _cold_kw.update(nr_min_step=_NR_COLD_MIN_STEP,
                            nr_max_iter=_NR_COLD_MAX_ITER)
        _t_cold = time.perf_counter()
        _sol = (None if _seed_first
                else _solve_fem_newton(fem_data, F, prep, **_cold_kw))
        if _sol is None:
            _sol = {"converged": False, "exit_reason": 'diverging', "iterations": 0,
                    "nr_force_evals": 0, "nr_predictor_iterations": 0,
                    "softened_1d_elements": None, "_nr_cold_skipped": True}
        # Cost attribution, bookkeeping only (see SPIKE.md, "THE COST OF THE RESCUE").
        # Nothing below is ever read for a verdict; it exists so a trial's wall time
        # can be split between the cold attempt, the predictor rungs and the seeded
        # correctors without re-deriving it from a log.
        _sol["nr_cold_wall"] = time.perf_counter() - _t_cold
        _sol["nr_cold_iterations"] = int(_sol.get("iterations", 0) or 0)
        _sol["nr_cold_force_evals"] = int(_sol.get("nr_force_evals", 0) or 0)
        _sol["nr_rungs"] = []
        # The viscoplastic predictor (see _NR_VP_PREDICTOR_ITERS). Only a trial that
        # died at the LOAD-STEP FLOOR is retried — a trial that reached full gravity
        # and was refused on force or on displacement has an answer already, and
        # re-deriving it from a different starting state would be re-litigating a
        # verdict rather than reaching one.
        if (_nr_export is None and not _sol["converged"]
                and _sol.get("exit_reason") == 'diverging'):
            # The post-peak set the failed attempt had reached, so both the seed and
            # the corrector that reads it are built on the same constitutive law. A
            # seed grown with every bar at its peak would be a different model.
            _soft_seed = _sol.get("softened_1d_elements", None)
            if _soft_seed is not None and not np.any(_soft_seed):
                _soft_seed = None
            _cold_rec = dict(wall=_sol["nr_cold_wall"],
                             iterations=_sol["nr_cold_iterations"],
                             force_evals=_sol["nr_cold_force_evals"],
                             rungs=[])
            _rungs = _nr_predictor_rungs(
                prep, n_elements, max_iterations, max_iterations_ceiling,
                early_failure)
            if _nr_rescue_rungs is not None:
                # P2: a trial far above a standing bound gets a truncated chain. A
                # seed-first trial (P3) has no cold attempt to fall back on, so it
                # always keeps at least one rung — a trial with no attempt at all
                # would be a verdict nothing measured.
                _keep = int(_nr_rescue_rungs)
                if _seed_first:
                    _keep = max(_keep, 1)
                _rungs = _rungs[:max(0, _keep)]
            for (_chunk, _ceiling, _prep_seed,
                 _pred_early_failure) in _rungs:
                _t_pred = time.perf_counter()
                _vp = solve_fem(
                    fem_data, F=F, debug_level=0, max_iterations=_chunk,
                    _softened_seed=_soft_seed,
                    tolerance=tolerance, max_disp_factor=None,
                    tension_cutoff=tension_cutoff,
                    dt_scale=dt_scale,
                    force_tol=force_tol, oob_window=oob_window,
                    early_exit=early_exit, min_slip_depth=min_slip_depth,
                    ssr_exclude_mask=ssr_exclude_mask,
                    tension_cap_by_elem=tension_cap_by_elem,
                    tension_srf=tension_srf, elastic_mask=elastic_mask,
                    suction_phi_b=suction_phi_b, suction_cap=suction_cap,
                    _prepared=_prep_seed, fast_kernel=fast_kernel,
                    failure_criterion=failure_criterion,
                    max_iterations_ceiling=_ceiling,
                    early_failure=_pred_early_failure,
                    k0=k0, _init_state=_init_state,
                    fem_solver='viscoplastic')
                _pred_wall = time.perf_counter() - _t_pred
                # On a K0 model the predictor's reported displacement is measured
                # from the in-situ datum, so the seed is read off its `_k0_state`,
                # which carries the ABSOLUTE field and the plastic strain in this
                # prepared model's own group order. Without K0 the reported
                # displacement IS the absolute one and the legacy per-element form
                # is used unchanged, which is what keeps that path bit-identical.
                _seed_kw = (dict(_nr_seed_state=_vp["_k0_state"]) if k0 is not None
                            else dict(_nr_seed=(_vp["displacements"],
                                                _vp["plastic_strains"])))
                # The predictor is allowed to have latched further than the cold
                # attempt did; the corrector takes whichever set the seed was
                # actually grown on, and its own latch continues from there.
                _vp_soft = _vp.get("softened_1d_elements", None)
                if _vp_soft is not None and np.any(_vp_soft):
                    _seed_kw['_softened_seed'] = _vp_soft
                elif _soft_seed is not None:
                    _seed_kw['_softened_seed'] = _soft_seed
                _t_seed = time.perf_counter()
                _retry = _solve_fem_newton(fem_data, F, prep, **_seed_kw, **_nr_kw)
                _seed_wall = time.perf_counter() - _t_seed
                _cold_rec['rungs'].append(dict(
                    chunk=int(_chunk),
                    pred_iterations=int(_vp.get("iterations", 0) or 0),
                    pred_wall=_pred_wall,
                    seed_iterations=int(_retry.get("iterations", 0) or 0),
                    seed_force_evals=int(_retry.get("nr_force_evals", 0) or 0),
                    seed_wall=_seed_wall,
                    converged=bool(_retry.get("converged", False))))
                # Work is CUMULATIVE: the failed cold attempt, every predictor run
                # and every corrector are all charged to this trial, so no cost is
                # hidden by the retry succeeding.
                _retry["iterations"] += _sol["iterations"]
                _retry["nr_force_evals"] += _sol["nr_force_evals"]
                _retry["nr_predictor_iterations"] = (
                    _sol.get("nr_predictor_iterations", 0)
                    + int(_vp.get("iterations", 0)))
                _sol = _retry
                if _sol["converged"]:
                    break
            _sol["nr_cold_wall"] = _cold_rec['wall']
            _sol["nr_cold_iterations"] = _cold_rec['iterations']
            _sol["nr_cold_force_evals"] = _cold_rec['force_evals']
            _sol["nr_rungs"] = _cold_rec['rungs']
        # L4: the coarse attempt was a PRE-FILTER, so the trial does not get to fail
        # on it. The full cold attempt the shipped driver would have made runs here,
        # from zero, on the shipped step control, and its work is charged to this
        # trial like every other attempt.
        if (_NR_COLD_FALLBACK and _NR_COLD_CHEAP and _nr_export is None
                and not _sol["converged"] and not _seed_first):
            _t_full = time.perf_counter()
            _full = _solve_fem_newton(fem_data, F, prep, **_nr_kw)
            _full["iterations"] += _sol["iterations"]
            _full["nr_force_evals"] += _sol["nr_force_evals"]
            _full["nr_predictor_iterations"] = _sol.get("nr_predictor_iterations", 0)
            _full["nr_cold_wall"] = (_sol.get("nr_cold_wall", 0.0)
                                     + time.perf_counter() - _t_full)
            _full["nr_rungs"] = _sol.get("nr_rungs", [])
            _sol = _full
        return _sol

    if debug_level >= 1:
        print(f"  c: {c_by_elem[0]:.1f} -> {c_reduced[0]:.1f}")
        print(f"  phi: {phi_by_elem[0]:.1f} -> {np.degrees(phi_reduced[0]):.1f}")
        if suction_active:
            _nsuc = int(np.count_nonzero(suction_tanphib_by_elem > 0.0))
            print(f"  matric suction: phi_b>0 on {_nsuc}/{n_elements} elements "
                  f"(apparent cohesion reduced by F)")

    # Extract 1D truss element data
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    n_1d_elements = len(elements_1d)
    has_1d_elements = n_1d_elements > 0

    if has_1d_elements:
        k_by_1d_elem = fem_data["k_by_1d_elem"]
        t_allow_by_1d_elem = fem_data["t_allow_by_1d_elem"]
        t_res_by_1d_elem = fem_data["t_res_by_1d_elem"]
        cos_theta_1d = fem_data["cos_theta_1d"]
        sin_theta_1d = fem_data["sin_theta_1d"]
        dof_indices_1d = fem_data["dof_indices_1d"]
        # How many of each row's DOF slots the element uses: 4 for a two-node
        # bar, 6 for a three-node bar that also stands on the midside node of
        # its soil edge.
        n_dof_1d = fem_data.get("n_dof_by_1d_elem",
                                np.full(n_1d_elements, 4, dtype=int))

        # Tracking arrays for 1D element status
        forces_1d = np.zeros(n_1d_elements)
        failed_1d = np.zeros(n_1d_elements, dtype=bool)
        # Post-peak set. A bar drops from t_allow to t_res ONLY by entering this
        # set, and it is only ever updated on a CONVERGED state (see the
        # softening fixed point at the convergence break). Never mid-iteration:
        # the first viscoplastic iterate is the elastic predictor, whose bar
        # forces overshoot wildly before the soil sheds load into them, so a
        # force-triggered latch inside the loop condemns bars for a transient
        # that never physically existed. Softening is possible only where t_res
        # is FINITE — an unset (NaN) t_res means the bar is
        # elastic-perfectly-plastic and holds t_allow forever.
        softened_1d = np.zeros(n_1d_elements, dtype=bool)
        # A carried-in post-peak set (INTERNAL; solve_ssrm's at-failure capture).
        # The capture solve runs beyond critical, where the slope never passes
        # through equilibrium, so the fixed point below can never fire and every
        # softening bar would sit at its peak capacity in the at-failure field.
        # Seeding it with the set the bracket's failed-edge trial actually shed
        # to shows those bars at their residual, which is the state that failed.
        if _softened_seed is not None and len(_softened_seed) == n_1d_elements:
            softened_1d |= np.asarray(_softened_seed, dtype=bool)
        can_soften_1d = (np.isfinite(t_res_by_1d_elem)
                         & (t_res_by_1d_elem < t_allow_by_1d_elem - 1e-12))
        n_soften_rounds = 0

        if debug_level >= 1:
            print(f"  1D truss elements: {n_1d_elements}")

    # Extract pile beam element data (Euler-Bernoulli)
    n_pile_elements = fem_data.get("n_pile_elements", 0)
    has_pile_elements = n_pile_elements > 0
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(n_1d_elements, dtype=bool))

    if has_pile_elements:
        cos_theta_pile = fem_data["cos_theta_pile"]
        sin_theta_pile = fem_data["sin_theta_pile"]
        dof_indices_pile = fem_data["dof_indices_pile"]
        # 6 on a two-node beam element, 9 on a three-node one.
        n_dof_pile = fem_data.get("n_dof_by_pile_elem",
                                  np.full(n_pile_elements, 6, dtype=int))
        K_local_pile = fem_data.get("K_local_by_pile_elem", None)
        V_cap_pile = fem_data["V_cap_by_pile_elem"]
        M_cap_pile = fem_data["M_cap_by_pile_elem"]
        L_pile_elem = fem_data["elem_length_by_pile_elem"]
        EI_pile = fem_data["EI_by_pile_elem"]
        EA_pile = fem_data["EA_by_pile_elem"]
        pile_head_nodes = fem_data.get("pile_head_nodes", np.array([], dtype=int))
        pile_head_fixed = fem_data.get("pile_head_fixed", np.array([], dtype=bool))
        forces_pile_axial = np.zeros(n_pile_elements)
        forces_pile_lateral = np.zeros(n_pile_elements)
        forces_pile_moment = np.zeros((n_pile_elements, 2))  # [M1, M2] at each node
        yielded_pile_V = np.zeros(n_pile_elements, dtype=bool)
        yielded_pile_M = np.zeros(n_pile_elements, dtype=bool)
        # Plastic hinge rotation at the two END nodes of each element -- the
        # rotation the moment capacity released, zero wherever no hinge formed.
        plastic_rot_pile = np.zeros((n_pile_elements, 2))
        # Which element end may carry that release: one per pile NODE, since a
        # hinge is one release at one section and every interior node of a pile
        # carries two element ends. F-independent, so it is built once here.
        pile_hinge_ends = _pile_hinge_ends(
            fem_data.get("pile_elem_nodes", None), n_pile_elements)

        if debug_level >= 1:
            _n_node_pile = sorted({int(n) // 3 for n in n_dof_pile}) or [2]
            _pile_kind = "/".join(f"{n}-node" for n in _n_node_pile)
            print(f"  Pile beam elements: {n_pile_elements} "
                  f"({_pile_kind} Euler-Bernoulli)")

    # ---- Working Gauss-point groups: the F-DEPENDENT half, rebuilt each solve ----
    # The prepared model carries the F-INDEPENDENT halves of each group (geometry
    # B/D4/w/dof, the per-GP dt_r, and the suction / elastic / power-curve /
    # Hoek-Brown parameter arrays). Here we attach THIS trial's F-reduced strengths
    # and a FRESH zeroed viscoplastic-strain buffer, so a reused prepared model can
    # never carry a stale strength or a stale plastic state between SSRM trials. The
    # per-GP F-dependent arrays are gathered by fancy-indexing the per-element
    # quantities with the group's cached element index e_idx — bit-identical to the
    # original per-Gauss-point list comprehensions. Static arrays (B, D4, w, dof,
    # dt_r, and the suction/elastic/pow/hb parameters) are shared by reference; the
    # viscoplastic loop only ever reads them and mutates the fresh per-solve arrays
    # (evp, and — for power-curve / Hoek-Brown Gauss points — c_r, snph, csph).
    gp_groups = []
    for _gi, _sg in enumerate(prep["gp_groups_static"]):
        _e_idx = _sg['e_idx']
        _G = _sg['n']
        grp = {
            'pairs': _sg['pairs'],
            'B': _sg['B'], 'D4': _sg['D4'], 'w': _sg['w'], 'dof': _sg['dof'],
            'dt_r': _sg['dt_r'],
            'c_r': c_reduced[_e_idx],
            'phi_r': phi_reduced[_e_idx],
            'F': F_by_elem[_e_idx],
            't_cap': t_cap_by_elem[_e_idx],
            # Fresh per solve — or the carried equilibrated strain, which together
            # with u_datum reproduces the equilibrated stress field at iteration 0.
            'evp': (np.zeros((_G, 4)) if _init_evp is None
                    else np.array(_init_evp[_gi], dtype=float, copy=True)),
        }
        grp['snph'] = np.sin(grp['phi_r'])
        grp['csph'] = np.cos(grp['phi_r'])
        grp['has_cap'] = bool(np.isfinite(grp['t_cap']).any())
        if sv0_gp is not None:
            grp['sv0'] = np.array([sv0_gp[e][g] for e, g in _sg['pairs']])
        if suction_active:
            grp['tanphib'] = _sg['tanphib']
            grp['scap'] = _sg['scap']
            grp['Finv'] = 1.0 / grp['F']
        if _sg.get('has_elastic'):
            grp['elastic'] = _sg['elastic']
            grp['has_elastic'] = True
        if 'pow_m' in _sg:
            grp['pow_m'] = _sg['pow_m']
            for _k in ('pow_a', 'pow_b', 'pow_cp', 'pow_d'):
                grp[_k] = _sg[_k]
        if 'hb_m' in _sg:
            grp['hb_m'] = _sg['hb_m']
            for _k in ('hb_sci', 'hb_mb', 'hb_s', 'hb_a'):
                grp[_k] = _sg[_k]
        gp_groups.append(grp)

    # ---- Compiled Mohr-Coulomb kernel ('auto' by default; NumPy path is the oracle) ----
    # When fast_kernel is on, MC-only groups (no power-curve / Hoek-Brown Gauss
    # points) run their Step-6 constitutive update in the compiled kernel; every
    # other group, and all 1D/pile work, stays on the NumPy reference below. The
    # kernel needs contiguous intp dof indices and uint8 elastic flags plus a
    # shared zero buffer; these are static per solve, so they are built once here.
    # The compiled module is NOT shipped by pip: it is built locally with
    # `python setup_kernel.py build_ext --inplace` (needs Cython), which the
    # installer builds do. Hence the 'auto' default: use it when it is there, use
    # the NumPy reference when it is not, and never make a run's success depend on
    # which of the two a machine happens to have. The two paths were certified to
    # give identical answers over all 103 FEM benchmarks on 2026-07-26.
    #   'auto'  -> compiled if importable, reference otherwise, silently
    #   True    -> require it; if not built, warn and fall back (unchanged)
    #   False   -> never; pure NumPy reference
    #
    # ORACLE DOCTRINE (see the fast_kernel parameter doc above for the full version):
    # the NumPy reference alone DEFINES every locked and published factor of safety;
    # the compiled kernel must reproduce it and is never itself the definition. The
    # verification suite therefore pins its kernel explicitly instead of inheriting
    # this default (fast-first with reference fallback on a lock miss;
    # --reference-only pinned to False), and benchmarks/kernel_xcheck.py is the
    # required divergence fence that keeps 'auto' safe.
    _mc_kernel = None
    _kernel_required = (fast_kernel is True)
    # The compiled Step-6 kernel computes sigma = D(Bu - evp) + u*m internally and
    # has no slot for a per-Gauss-point INITIAL stress, so a K0 run stays entirely on
    # the NumPy reference — which is the oracle anyway. Runs without k0 are untouched.
    if sv0_gp is not None:
        fast_kernel = False
    if fast_kernel:
        try:
            from xslope import _fem_kernel as _mc_kernel
        except ImportError:
            if _kernel_required:
                import warnings
                warnings.warn(
                    "fast_kernel=True but the compiled xslope._fem_kernel is not built; "
                    "falling back to the NumPy reference path. Build it with "
                    "`python setup_kernel.py build_ext --inplace` (requires Cython).",
                    RuntimeWarning, stacklevel=2)
            _mc_kernel = None
    if _mc_kernel is not None:
        for grp in gp_groups:
            if 'pow_m' in grp or 'hb_m' in grp:
                grp['_fast'] = False
                continue
            grp['_fast'] = True
            _G = grp['dof'].shape[0]
            grp['_dof_intp'] = np.ascontiguousarray(grp['dof'], dtype=np.intp)
            grp['_zeroG'] = np.zeros(_G, dtype=np.float64)
            grp['_B_c'] = np.ascontiguousarray(grp['B'], dtype=np.float64)
            grp['_D4_c'] = np.ascontiguousarray(grp['D4'], dtype=np.float64)
            grp['_w_c'] = np.ascontiguousarray(grp['w'], dtype=np.float64)
            grp['_dtr_c'] = np.ascontiguousarray(grp['dt_r'], dtype=np.float64)
            grp['_tcap_c'] = np.ascontiguousarray(grp['t_cap'], dtype=np.float64)
            if grp.get('has_elastic'):
                grp['_elastic_u8'] = np.ascontiguousarray(grp['elastic'], dtype=np.uint8)
            else:
                grp['_elastic_u8'] = np.zeros(_G, dtype=np.uint8)

    # A one-iteration increment never decays on a settled slope (period-2 yield-surface
    # flicker), so a window of 1 would make convergence unreachable. Refuse it loudly rather
    # than silently returning a factor of safety driven by a numerical artifact.
    if int(oob_window) < 2:
        raise ValueError(
            f"oob_window must be >= 2 (got {oob_window}). A single-iteration increment does "
            "not decay on a settled slope: Gauss points on the yield surface flip flow "
            "direction every iteration, and that period-2 mode never vanishes at any dt or "
            "iteration count, so a stable slope would be reported as failing.")
    oob_window = int(oob_window)

    # ---- Step 8: Initialize viscoplastic strains (zero) ----
    # evp[elem_idx][gp_idx] = array of shape (4,): [ex, ey, gxy, ez]
    # (4-component plane strain after Smith & Griffiths: plastic eps_z is
    # tracked so sigma_z can relax; total eps_z = 0 so elastic eps_z = -evp_z)
    evp = []
    for elem_idx in range(n_elements):
        n_gp = len(elem_gp_data[elem_idx])
        evp.append([np.zeros(4) for _ in range(n_gp)])


    # ---- Step 9: Viscoplastic iteration loop ----
    # Self weight, the applied boundary forces (e.g. reservoir pressure) and the
    # pore-pressure field are all applied together, in one load stage.
    stage_list = [(F_gravity, u_gp, u_gp_signed, None)]

    total_iterations = 0
    u = np.zeros(n_dof)
    converged = False
    iteration = 0
    unbalanced_force_ratio = 0.0
    # Hybrid-criterion instrumentation. `disp_hist` samples max|u| every
    # _HYBRID_SAMPLE_EVERY iterations; the value sampled is `norm_u_new`, which the
    # CHECON test already computes, so the only cost is a list append every tenth
    # iteration. Both are reset per STAGE (u_elastic is per-stage, so a history
    # carried across a stage boundary would be measured against the wrong yardstick);
    # the classifier therefore reads the history of whichever stage did not settle.
    disp_hist = []
    u_elastic_scale = 0.0
    exit_reason = 'iteration_cap'
    plateau_iter = None            # iteration at which the residual plateaued
    plateau_ratio = None           # the out-of-balance ratio it plateaued at
    diverging_iter = None          # iteration at which the early-failure rule fired
    diverging_signal = None        # which of its two tests fired
    n_extensions = 0               # budget extensions granted (all stages)
    budget = int(max_iterations)
    sq3 = np.sqrt(3.0)   # loop-invariant constant (hoisted out of the VP iteration)

    for stage_idx, (base_loads, u_gp_active, u_gp_signed_active, stage_label) in enumerate(stage_list):
        if debug_level >= 1 and stage_label is not None:
            print(f"  {stage_label}")

        # flatten this stage's per-GP pore pressures into the groups
        for grp in gp_groups:
            grp['u_gp'] = np.array([u_gp_active[e][g] for e, g in grp['pairs']])
            if suction_active:
                # Matric-suction apparent cohesion for this stage: s = max(0,
                # -u_signed) capped at scap, times tan(phi_b), reduced by the trial
                # F (Finv = 1/F). Independent of the effective-normal u_gp (which is
                # still clamped and — under 'effective' — moved to the load vector
                # below); this reads the SIGNED field so the suction above the water
                # table is not lost to the clamp. Rebuilt each stage so the dry
                # stage credits no suction.
                if u_gp_signed_active is not None:
                    u_sgn = np.array([u_gp_signed_active[e][g] for e, g in grp['pairs']])
                    s_suc = _suction_capped(u_sgn, grp['scap'])
                    grp['c_suc_r'] = _suction_apparent_cohesion(
                        grp['tanphib'], s_suc, grp['Finv'])
                else:
                    grp['c_suc_r'] = np.zeros(len(grp['pairs']))

        # ---- K0 initial stress (per stage) ----
        # Build this stage's initial stress state sigma_0 at every Gauss point and
        # move it to the right-hand side. With
        #     sigma = sigma_0 + D (B u - evp)                       (tension-positive)
        # equilibrium int B^T sigma dV = F_ext becomes
        #     K u = F_ext - int B^T sigma_0 dV + int B^T D evp dV,
        # so the ONLY changes are one extra load term here and one extra addend at
        # the yield check. This is the classical initial-stress method: the solver
        # still iterates to equilibrium under the body forces, it just starts from a
        # K0 state instead of from zero stress. When sigma_0 happens to be in
        # equilibrium with gravity (level ground, K0 = nu/(1-nu)) the displacement
        # solution is ~0; on a slope it is not, and the iteration redistributes.
        #
        # sigma_0 is EFFECTIVE:
        #     sigma'_v = -(soil overburden) + u        (u >= 0, tension-positive)
        #     sigma'_h = sigma'_z = K0 * sigma'_v      (in-plane AND out-of-plane)
        #     tau_xy   = 0
        # matching the load vector the stage actually applies.
        F_sig0 = None
        if sv0_gp is not None:
            k0f = float(k0)
            F_sig0 = np.zeros(n_dof)
            for grp in gp_groups:
                u_st = grp['u_gp']
                sv_eff = -grp['sv0'] + u_st
                sh_eff = k0f * sv_eff
                z = np.zeros_like(sv_eff)
                sig0 = np.stack([sh_eff, sv_eff, z, sh_eff], axis=1)
                grp['sig0'] = sig0
                contrib = np.einsum('gij,gi->gj', grp['B'], sig0[:, :3]) * grp['w'][:, None]
                np.add.at(F_sig0, grp['dof'].ravel(), contrib.ravel())

        # Effective-stress formulation: equilibrium of sigma_total =
        # sigma_eff - u*m (tension-positive) gives
        #   int B^T sigma_eff dV = F_ext + int B^T m u dV,
        # so the pore-pressure term joins the load vector and D*B*u is the
        # EFFECTIVE stress directly (no subtraction at the yield check).
        F_u = np.zeros(n_dof)
        for grp in gp_groups:
            contrib = (grp['w'] * grp['u_gp'])[:, None] * (
                grp['B'][:, 0, :] + grp['B'][:, 1, :])
            np.add.at(F_u, grp['dof'], contrib)
            grp['u_gp'] = np.zeros_like(grp['u_gp'])
        base_loads = base_loads + F_u

        # ---- The hybrid criterion's yardstick ----
        # `loads_grav` is this stage's APPLIED load — self weight, water, surcharge —
        # before the initial-stress term is moved to the right-hand side. Its elastic
        # response is the displacement scale the hybrid thresholds (1.25 stuck ceiling,
        # 1.5 failed floor, 0.02 growth) were calibrated against, and it is the scale a
        # run without K0 measures itself by, so K0-on and K0-off verdicts are read off
        # the same ruler.
        #
        # The load the solver actually applies, `base_loads`, is that load MINUS the
        # initial stress's internal forces. With K0 that difference is only the part of
        # the weight the K0 field does not already carry (near zero once the state is
        # equilibrated, and never the whole weight), so its elastic response is a
        # residual, not a displacement scale. Using it would shrink the denominator by
        # whatever fraction of the load the initial stress happens to balance and move
        # the thresholds by the same factor.
        loads_grav = base_loads
        if F_sig0 is not None:
            base_loads = base_loads - F_sig0   # rebinds; loads_grav keeps the applied load

        # per-stage elastic reference (pure elastic response to this stage's loads)
        u_e_free = K_factor.solve(base_loads[free_dofs])
        u_elastic = np.zeros(n_dof)
        u_elastic[free_dofs] = u_e_free
        # One extra back-substitution on the SAME factorization, and only when K0 is
        # active; without it the two vectors are identical.
        u_e_grav = (u_e_free if F_sig0 is None
                    else K_factor.solve(loads_grav[free_dofs]))
        if stage_idx == 0:
            # Start from the elastic solution — or, with an equilibrated state
            # carried in, from that state, which is already this solve's answer at
            # full strength and its zero of displacement.
            u = u_elastic.copy() if _init_u is None else u_datum.copy()
        if debug_level >= 1 and stage_idx == 0:
            print(f"  Initial elastic: max|u| = {np.max(np.abs(u)):.6f}")

        converged = False
        unbalanced_force_ratio = 0.0
        # Reset per STAGE: base_loads changes at a stage boundary, so a history carried
        # across it would measure a load step, not a residual.
        loads_hist = [base_loads.copy()]
        ufr_best = float('inf')        # lowest out-of-balance seen this stage
        last_progress_iter = 0         # iteration of last meaningful improvement
        disp_hist = []                 # max|u| samples (hybrid criterion)
        u_elastic_scale = float(np.max(np.abs(u_e_grav))) if u_e_grav.size else 0.0
        exit_reason = 'iteration_cap'
        plateau_iter = None            # no-progress watch, reset per stage
        plateau_ratio = None
        diverging_iter = None          # early-failure watch, reset per stage
        diverging_signal = None
        oob_hist = []                  # out-of-balance samples (budget-extension trend)
        budget = int(max_iterations)   # this stage's CURRENT budget; may be extended
        ceiling = max(int(max_iterations_ceiling or 0), budget)
        chunk = int(max_iterations)    # one extension is worth one original budget
        # Two trend windows have to FIT inside the budget, or the trend can never be
        # read and a small budget would never be extended. They are the nominal
        # width, capped at a quarter of the budget.
        trend_window = min(_OOB_TREND_WINDOW,
                           max(2 * _HYBRID_SAMPLE_EVERY, budget // 4))

        iteration = -1
        while True:
            iteration += 1
            if iteration >= budget:
                # Budget reached. Extend it while the solve is still progressing.
                if budget < ceiling and _still_progressing(
                        oob_hist, disp_hist, u_elastic_scale, mesh_height,
                        trend_window):
                    budget = min(ceiling, budget + chunk)
                    n_extensions += 1
                    if debug_level >= 1:
                        print(f"  Budget extended at iteration {iteration}: "
                              f"the solve is still making progress "
                              f"({unbalanced_force_ratio:.2e} against tolerance "
                              f"{force_tol:.1e}); budget now {budget} "
                              f"(ceiling {ceiling})")
                else:
                    if (budget >= ceiling
                            and _still_progressing(oob_hist, disp_hist,
                                                   u_elastic_scale, mesh_height,
                                                   trend_window)
                            and _oob_still_falling(oob_hist, window=trend_window)):
                        # Out of ceiling, not out of progress: this trial has not
                        # failed and has not converged, and nothing here can say
                        # which. The residual must be measurably STILL FALLING for
                        # that claim — an undecidable displacement field on its own
                        # keeps the legacy verdict, so the bisection is only ever
                        # halted on positive evidence of an unfinished convergence.
                        exit_reason = 'inconclusive'
                        if debug_level >= 1:
                            print(f"  Iteration ceiling {ceiling} reached with the "
                                  f"solve still making progress "
                                  f"({unbalanced_force_ratio:.2e} against tolerance "
                                  f"{force_tol:.1e}) - INCONCLUSIVE, neither "
                                  f"converged nor failed")
                    converged = False
                    iteration -= 1          # the last iteration actually performed
                    break
            # Build body load correction from accumulated viscoplastic strains
            _tp = time.perf_counter() if _PROF_ON else None
            loads = base_loads.copy()

            n_yielding = 0

            for grp in gp_groups:
                if _mc_kernel is not None and grp['_fast']:
                    # Compiled Mohr-Coulomb Step-6 (opt-in). Mutates grp['evp'] and
                    # scatters the body-load correction into loads; returns this
                    # group's MC-yielding count.
                    _c_suc = grp.get('c_suc_r')
                    n_yielding += _mc_kernel.mc_step6(
                        u, loads,
                        grp['_B_c'], grp['_D4_c'], grp['_w_c'], grp['_dof_intp'],
                        grp['_dtr_c'], grp['c_r'],
                        grp['_zeroG'] if _c_suc is None else _c_suc,
                        grp['snph'], grp['csph'], grp['_tcap_c'],
                        grp['u_gp'], grp['_elastic_u8'], grp['evp'],
                        dt, 1 if grp['has_cap'] else 0,
                        1 if grp.get('has_elastic') else 0)
                    continue
                Bg, D4g, wg = grp['B'], grp['D4'], grp['w']
                dofg, evpg = grp['dof'], grp['evp']
                u_e = u[dofg]                                   # (G, ndof)
                eps = np.einsum('gij,gj->gi', Bg, u_e)          # (G, 3)
                eps4 = np.empty((len(wg), 4))
                eps4[:, :3] = eps - evpg[:, :3]
                eps4[:, 3] = -evpg[:, 3]
                sig4 = np.einsum('gij,gj->gi', D4g, eps4)       # (G, 4) tension-positive
                _sig0 = grp.get('sig0')
                if _sig0 is not None:
                    sig4 = sig4 + _sig0        # K0 initial stress (see the stage loop)
                sig_eff = sig4.copy()
                sig_eff[:, [0, 1, 3]] += grp['u_gp'][:, None]

                sx, sy, txy, sz = sig_eff.T

                # Power-curve elements: re-linearize the F-reduced envelope
                # tau_F = [a*(s+d)^b + c_p]/F at the CURRENT effective normal
                # stress every iteration. Linearization point: the in-plane
                # Mohr-circle center s' = -(sx+sy)/2 (compression-positive) -
                # the failure-plane normal is implicit in phi_t, and over one
                # circle radius the envelope curvature is small, so the
                # center is a stable, vectorizable abscissa (the LEM uses the
                # slice-base normal; both converge on the same reduced
                # envelope). Guards mirror solve._pow_update_strength.
                pm = grp.get('pow_m')
                if pm is not None:
                    Fpm = grp['F'][pm]           # per-GP F (1.0 where SSR-excluded)
                    s_n = np.maximum(-(sx[pm] + sy[pm]) * 0.5, 0.0)
                    ref = max(1.0, float(s_n.mean()) if s_n.size else 1.0)
                    s_ef = np.maximum(s_n + grp['pow_d'][pm], 1e-4 * ref)
                    bb = grp['pow_b'][pm]
                    slope_t = grp['pow_a'][pm] * bb * s_ef ** (bb - 1.0) / Fpm
                    tau_F = (grp['pow_a'][pm] * s_ef ** bb + grp['pow_cp'][pm]) / Fpm
                    phi_t = np.arctan(slope_t)
                    grp['c_r'][pm] = tau_F - s_n * slope_t
                    grp['snph'][pm] = np.sin(phi_t)
                    grp['csph'][pm] = np.cos(phi_t)

                # Hoek-Brown elements: linearize the envelope at the normal stress on
                # the FAILURE PLANE, sigma_n = s'cos^2(phi) - c sin(phi)cos(phi), taken
                # from the previous iterate's REDUCED tangent. That expression is the
                # exact point at which a Mohr circle touches its tangent line, so it
                # closes as a fixed point inside the VP loop: at convergence the circle,
                # the tangent line and the reduced envelope all meet at the same
                # sigma_n. It is also the abscissa the LEM uses (the slice-base normal
                # stress), so the two solvers linearize the same curve at the same place.
                #
                # Do NOT linearize at the minor principal stress sigma3 instead. Balmer's
                # sigma3 -> tangency mapping is derived for the UNREDUCED envelope, so it
                # names the right point only at F = 1; under strength reduction the
                # F-times-smaller Mohr circle contacts the reduced envelope at a much
                # lower normal stress, and because the HB envelope is CONCAVE the tangent
                # taken at the old abscissa lies strictly above it -- a one-sided,
                # over-strong yield surface. Measured on Hammah et al. (2005) Example 1
                # this cost +6% in FS; the failure-plane abscissa reproduces their
                # published SSR 1.15 to +0.7%.
                # The circle is the IN-PLANE one (sx, sy, txy); the out-of-plane sz is not
                # folded in. That was checked, not assumed: on the Hammah benchmark sz is
                # below the in-plane minor principal at 0.0% of points (and 0.0% of
                # YIELDING points), so the in-plane circle is the critical one wherever
                # the linearization actually matters. It cannot be otherwise in the
                # yielding band -- plane strain gives sz ~ nu(sx+sy), which can only drop
                # below the in-plane minor where the deviatoric radius is small, i.e. in
                # deep near-hydrostatic material that is not yielding. Building the circle
                # from the true 3D major/minor instead moved the SSRM factor of safety by
                # less than the bracket tolerance (1.161 either way).
                hm = grp.get('hb_m')
                if hm is not None:
                    ctr_t = (sx[hm] + sy[hm]) * 0.5
                    s_prime = -ctr_t             # circle center, compression-positive
                    sn_p, cs_p = grp['snph'][hm], grp['csph'][hm]
                    # Clamping sigma_n >= 0 bounds phi_i at its zero-normal-stress value
                    # (~60 deg for a typical rock mass) rather than letting a Gauss point
                    # in tension run to the tensile apex, where dsigma1/dsigma3 diverges
                    # and phi_i -> 90 deg (i.e. effectively infinite strength).
                    s_n = np.maximum(s_prime * cs_p ** 2
                                     - grp['c_r'][hm] * sn_p * cs_p, 0.0)
                    c_i, phi_i = hb_tangent_const(
                        s_n, grp['hb_sci'][hm], grp['hb_mb'][hm],
                        grp['hb_s'][hm], grp['hb_a'][hm], iters=40)
                    # Strength reduction divides the SHEAR strength by F, which divides
                    # both the instantaneous cohesion and tan(phi_i) by F -- dividing a
                    # function by F divides its tangent's slope and intercept by F. The HB
                    # constants themselves are NOT reduced -- sigma_ci/F is a different
                    # envelope, because of the exponent a.
                    Fhm = grp['F'][hm]           # per-GP F (1.0 where SSR-excluded)
                    slope_t = np.tan(np.radians(phi_i)) / Fhm
                    phi_t = np.arctan(slope_t)
                    grp['c_r'][hm] = c_i / Fhm
                    grp['snph'][hm] = np.sin(phi_t)
                    grp['csph'][hm] = np.cos(phi_t)

                sigm = (sx + sy + sz) / 3.0
                dsbar = np.sqrt(((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2
                                 + 6.0 * txy**2) / 2.0)
                dxv, dyv, dzv = sx - sigm, sy - sigm, sz - sigm
                ds3 = np.maximum(dsbar, 1e-10)**3
                sine = np.clip(np.where(dsbar > 1e-10,
                                        -13.5 * (dxv * dyv * dzv - dzv * txy**2) / ds3,
                                        0.0), -1.0, 1.0)
                theta = np.arcsin(sine) / 3.0
                snth, csth = np.sin(theta), np.cos(theta)
                # Cohesion in the MC envelope: the F-reduced c', plus the opt-in
                # matric-suction apparent cohesion c_suc_r (already reduced by F).
                # When suction is inactive the key is absent and c_env is exactly
                # grp['c_r'] (bit-identical to the pre-suction yield).
                _c_suc = grp.get('c_suc_r')
                c_env = grp['c_r'] if _c_suc is None else grp['c_r'] + _c_suc
                f = (sigm * grp['snph']
                     + dsbar * (csth / sq3 - snth * grp['snph'] / 3.0)
                     - c_env * grp['csph'])

                m = (f > 0) & (dsbar > 1e-20)
                if grp.get('has_elastic'):
                    # Pure-elastic elements are held out of plasticity entirely:
                    # drop their Gauss points from the yielding set so they never
                    # accrue viscoplastic strain (evpg stays 0 -> elastic stress).
                    m = m & ~grp['elastic']
                n_yielding += int(np.count_nonzero(m))
                if np.any(m):
                    # vectorized MOCOUQ flow, psi = 0 (dq1 = 0); corner freeze
                    dsb, th = dsbar[m], theta[m]
                    snt, cst = snth[m], csth[m]
                    dxm, dym, dzm, txm = dxv[m], dyv[m], dzv[m], txy[m]
                    xj2 = dsb**2 / 3.0
                    a2 = (3.0 / (2.0 * dsb))[:, None] * np.stack(
                        [dxm, dym, 2.0 * txm, dzm], axis=1)
                    a3 = np.stack([dxm**2 + txm**2 - (2.0/3.0) * xj2,
                                   dym**2 + txm**2 - (2.0/3.0) * xj2,
                                   -2.0 * dzm * txm,
                                   dzm**2 - (2.0/3.0) * xj2], axis=1)
                    corner = np.abs(snt) > 0.49
                    cs3, tn3 = np.cos(3.0 * th), np.tan(np.where(corner, 0.0, 3.0 * th))
                    K = cst / sq3
                    Kp = -snt / sq3
                    C2 = np.where(corner, 0.5, K - Kp * tn3)
                    C3 = np.where(corner, 0.0,
                                  -4.5 * Kp / (np.where(corner, 1.0, cs3) * dsb**2))
                    flow = C2[:, None] * a2 + C3[:, None] * a3
                    evpg[m] += (f[m] * dt)[:, None] * flow

                if grp['has_cap']:
                    # ---- Rankine tension cutoff (second yield surface) ----
                    # Cap the MAJOR (most-tensile) in-plane principal stress at the
                    # per-element cap T (tension-positive): F_t = sigma_1 - T. Where
                    # F_t > 0, accrue viscoplastic strain along the ASSOCIATED flow
                    # normal n = d(sigma_1)/d(sigma), damped by dt_r. This is a
                    # principal-stress (Rankine) cap, the form RS2/PLAXIS/FLAC use,
                    # and a SEPARATE surface from the Mohr-Coulomb shear surface
                    # above: where both are active the two contributions simply SUM
                    # into evpg (Koiter's rule for the MC/tension corner). The psi=0
                    # MC flow is purely deviatoric and cannot relax a near-apex
                    # tensile state, which is exactly what this surface handles.
                    # Capping the in-plane major also bounds sigma_z
                    # (sigma_z <= sigma_1 in plane strain, since sigma_1 >= sx,sy and
                    # sigma_z = nu(sx+sy) here), so EVERY principal stress is held
                    # <= T. inf caps never fire, so capless elements are untouched;
                    # the global tension_cutoff flag is the T = 0 special case
                    # (no principal tension permitted anywhere).
                    # Refs: Smith & Griffiths (viscoplastic Mohr-Coulomb, Ch.6,
                    # Progs 6.11-6.13) for the initial-strain framework; Koiter
                    # (1960) / Owen & Hinton (1980) for multi-surface (corner)
                    # summation of viscoplastic flows.
                    cap = grp['t_cap']
                    ctr = 0.5 * (sx + sy)                        # circle center
                    Rc = np.sqrt((0.5 * (sx - sy))**2 + txy**2)  # circle radius
                    s1 = ctr + Rc                                # major principal (tension +)
                    tm = s1 > cap
                    if grp.get('has_elastic'):
                        # Pure-elastic elements do not tension-relax either.
                        tm = tm & ~grp['elastic']
                    if np.any(tm):
                        # n = d(sigma_1)/d(sx,sy,txy); txy is the direct derivative
                        # (engineering-shear conjugacy, matching a2/a3 above), and
                        # d(sigma_1)/d(sz) = 0. The radius is floored so n stays
                        # bounded at the biaxial apex Rc -> 0 (there sx-sy -> 0 too,
                        # so n -> [1/2, 1/2, 0, 0], an isotropic in-plane relaxation).
                        Rf = np.maximum(Rc[tm], 1e-10)
                        half = 0.5 * (sx[tm] - sy[tm])
                        Ft = (s1[tm] - cap[tm]) * grp['dt_r'][tm]
                        evpg[tm, 0] += Ft * (0.5 + 0.5 * half / Rf)
                        evpg[tm, 1] += Ft * (0.5 - 0.5 * half / Rf)
                        evpg[tm, 2] += Ft * (txy[tm] / Rf)

                # body-load correction: B^T (D4 evp)[:3] * w, scattered to dofs
                s4 = np.einsum('gij,gj->gi', D4g, evpg)
                contrib = np.einsum('gij,gi->gj', Bg, s4[:, :3]) * wg[:, None]
                np.add.at(loads, dofg.ravel(), contrib.ravel())

            if _PROF_ON:
                _prof_add("vp_const", _tp)

            # ---- 1D Truss element body-force corrections ----
            if has_1d_elements:
                n_1d_compression = 0
                n_1d_exceeded = 0

                for elem_idx_1d in range(n_1d_elements):
                    if pile_elem_mask[elem_idx_1d]:
                        continue  # pile elements handled separately below
                    dof_idx = dof_indices_1d[elem_idx_1d][:n_dof_1d[elem_idx_1d]]
                    k = k_by_1d_elem[elem_idx_1d]
                    cos_t = cos_theta_1d[elem_idx_1d]
                    sin_t = sin_theta_1d[elem_idx_1d]

                    # Relative displacement of the two ends, projected along the
                    # element axis. On the three-node bar that chord elongation
                    # gives the axial force AT THE ELEMENT CENTER exactly, which
                    # is the station every reader of forces_1d already places it
                    # at, so the formula and its meaning are unchanged.
                    u_elem = u[dof_idx]  # [u_x0, u_y0, u_x1, u_y1, (u_xm, u_ym)]
                    du_x = u_elem[2] - u_elem[0]
                    du_y = u_elem[3] - u_elem[1]
                    delta = du_x * cos_t + du_y * sin_t

                    # Axial force: T = k * delta (positive = tension)
                    T = k * delta
                    forces_1d[elem_idx_1d] = T

                    # Bar constitutive law: tension-only, elastic-PERFECTLY-
                    # PLASTIC. The yield force is t_allow — the SAME capacity
                    # envelope the LEM applies at a reinforcement crossing
                    # (fileio.reinforce_available_tension: tensile strength
                    # tapered by the pullout ramps at both ends). The force the
                    # bar can actually deliver is the elastic k*delta clipped
                    # into [0, t_allow]; beyond that the bar yields and holds
                    # t_allow while the soil around it keeps straining.
                    #
                    # NOTE 1 (strength reduction): t_allow is NOT divided by F.
                    # Only the SOIL strength is reduced (see the note at the top
                    # of solve_fem). The reinforcement keeps its full structural
                    # capacity, so the reported FS is the factor on soil strength
                    # at which the *supported* slope fails. This is the RS2/Slide
                    # convention and is what the LEM does (it applies the full
                    # available tension as a resisting force, independent of FS).
                    #
                    # NOTE 2 (post-peak): a bar that has entered the softened set
                    # yields at t_res instead of t_allow. Membership in that set is
                    # decided ONLY on a converged state, by the fixed point at the
                    # convergence break below — never here, mid-iteration. See the
                    # comment at softened_1d for why.
                    T_cap = (t_res_by_1d_elem[elem_idx_1d]
                             if softened_1d[elem_idx_1d]
                             else t_allow_by_1d_elem[elem_idx_1d])
                    T_true = min(max(T, 0.0), T_cap)

                    if T < 0:
                        n_1d_compression += 1
                    elif T > T_cap:
                        failed_1d[elem_idx_1d] = True   # "has yielded", for reporting
                        n_1d_exceeded += 1

                    # Viscoplastic body-load correction. The global stiffness
                    # carries the bar's FULL elastic stiffness, so K*u contains
                    # the (uncapped) elastic bar force T. Exactly as for the 2D
                    # soil, where loads += B^T*D*evp so that the true internal
                    # force is K*u - loads_body, the bar's body load must be the
                    # part of the elastic force the bar cannot actually carry:
                    #
                    #     f_body = (T - T_true) * [-cos, -sin, +cos, +sin]
                    #
                    # Then K*u - f_body leaves T_true in the bar. Adding the
                    # OPPOSITE sign (T_true - T) — as this code did — turns the
                    # cap into an anti-cap: the bar ends up carrying 2T - T_true,
                    # i.e. it gets *stiffer* the more it is overloaded, so the
                    # reinforced slope can never be driven to failure and the SSR
                    # factor is insensitive to T_allow.
                    correction_T = T - T_true

                    if abs(correction_T) > 1e-30:
                        # Internal force pattern for a constant tension T:
                        # [-cos, -sin, +cos, +sin] at the two ends. On the
                        # three-node bar the consistent nodal forces for a
                        # constant axial force are [-T, +T, 0] in node order, so
                        # the midside node takes no share of the correction.
                        loads[dof_idx[0]] += correction_T * (-cos_t)
                        loads[dof_idx[1]] += correction_T * (-sin_t)
                        loads[dof_idx[2]] += correction_T * cos_t
                        loads[dof_idx[3]] += correction_T * sin_t

                if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                    print(f"    1D elements: {n_1d_compression} in compression, "
                          f"{n_1d_exceeded} exceeded capacity, "
                          f"{np.sum(failed_1d)} total failed")

            # ---- Pile beam element force computation and capacity checks ----
            if has_pile_elements:
                n_pile_yielded_V = 0
                n_pile_yielded_M = 0
                for p_idx in range(n_pile_elements):
                    n_dof_e = int(n_dof_pile[p_idx])
                    n_node_e = n_dof_e // 3
                    dof_idx = dof_indices_pile[p_idx][:n_dof_e]
                    cos_t = cos_theta_pile[p_idx]
                    sin_t = sin_theta_pile[p_idx]
                    L = L_pile_elem[p_idx]
                    EI_val = EI_pile[p_idx]
                    EA_val = EA_pile[p_idx]
                    Kl = K_local_pile[p_idx] if K_local_pile is not None else None

                    # [ux, uy, theta] at end 1, end 2 and -- on a three-node
                    # element -- the midside node of the soil edge.
                    u_elem = u[dof_idx]
                    T_beam = _beam_rotation(cos_t, sin_t, n_node_e)
                    u_local = T_beam @ u_elem

                    # The capacities enter as a viscoplastic body-load correction,
                    # exactly as the soil's plastic strain and the bar's tension
                    # cap do: the global stiffness carries the member's FULL
                    # elastic stiffness, so `K u - corrections` is the internal
                    # force this state is in equilibrium with. See
                    # _pile_element_capacity for the moment hinge, and for why the
                    # shear correction carries the bar's sign.
                    (T_force, V, M1, M2, corr_local, _p_rot,
                     yielded_V_e, yielded_M_e) = _pile_element_capacity(
                        u_local, cos_t, sin_t, L, EA_val, EI_val, Kl, n_node_e,
                        V_cap_pile[p_idx], M_cap_pile[p_idx],
                        hinge_ends=pile_hinge_ends[p_idx])

                    forces_pile_axial[p_idx] = T_force
                    forces_pile_lateral[p_idx] = V
                    forces_pile_moment[p_idx] = [M1, M2]

                    if yielded_V_e:
                        yielded_pile_V[p_idx] = True
                        n_pile_yielded_V += 1
                    if yielded_M_e:
                        yielded_pile_M[p_idx] = True
                        n_pile_yielded_M += 1

                    if corr_local.any():
                        f_global = T_beam.T @ corr_local
                        for k in range(n_dof_e):
                            loads[dof_idx[k]] += f_global[k]

                if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                    if n_pile_yielded_V > 0 or n_pile_yielded_M > 0:
                        print(f"    Pile elements: {n_pile_yielded_V} V-yielded, {n_pile_yielded_M} M-yielded")

            # ---- Out-of-balance force, per node (Dawson, Roth & Drescher 1999) ----
            #
            # This is an initial-stress viscoplastic scheme, so the solve below
            # enforces
            #       int B^T D (B u - evp) dV  =  F_ext
            # EXACTLY, using the evp that built this iteration's `loads`. The state
            # is therefore always in equilibrium with the stresses it was built from,
            # and what is still "out of balance" is the amount by which the
            # viscoplastic body load is STILL CHANGING: the increment of
            # (loads - base_loads) from one iteration to the next.
            #
            # When plastic flow genuinely ceases, that increment decays to zero and
            # the stress field is both admissible and in equilibrium — a stable slope.
            # When the slope is failing, plastic flow never ceases: the increment
            # plateaus at a non-zero value and keeps feeding displacement forever.
            #
            # The increment is STRICTLY LOCAL — it is non-zero only at nodes adjacent
            # to Gauss points whose plastic strain is still flowing. Elastic padding
            # contributes exactly zero. Taking the MAXIMUM of the per-node value,
            # each normalized by that node's own weight, therefore measures the
            # failure mechanism against itself, and inert material added to the mesh
            # cannot dilute it.
            # The increment is AVERAGED OVER A WINDOW of `oob_window` iterations rather than
            # taken between consecutive ones. A one-iteration increment does not decay on a
            # settled slope: Gauss points sitting exactly on the yield surface flip their
            # flow direction every iteration, and the resulting body-load flicker is a clean
            # PERIOD-2 limit cycle — measured cos(dL_n, dL_n-1) = -1.0000 with |dL| pinned at
            # a constant, on a model whose displacements were frozen to four decimals and
            # whose accumulated body load was frozen to seven significant figures. Damping dt
            # only scales its amplitude (it is proportional to dt); it never removes it, so no
            # iteration budget can clear it. Averaging over the window cancels it exactly,
            # while genuine plastic drift (cos = +1) passes through untouched. The result is
            # insensitive to the width — 10, 50 and 200 converge on the same F at the same
            # iteration — so this rejects a specific numerical mode, it is not a tuning knob.
            # Locality, and hence padding immunity, is unaffected: elastic material contributes
            # exactly zero over any window.
            loads_hist.append(loads.copy())
            if len(loads_hist) > oob_window + 1:
                loads_hist.pop(0)
            d_load = ((loads - loads_hist[0])
                      / min(oob_window, len(loads_hist) - 1)) * free_dof_mask
            r_node = np.sqrt(d_load[node_dof_x] ** 2 + d_load[node_dof_y] ** 2)
            oob_node = (r_node / g_node_den)[node_has_free]
            # min_slip_depth filter: take the maximum only over nodes deep enough to
            # count. With no filter (_deep_free_mask is None) this is the full set.
            _oob_for_max = oob_node if _deep_free_mask is None else oob_node[_deep_free_mask]
            unbalanced_force_ratio = float(np.max(_oob_for_max)) if _oob_for_max.size else 0.0
            if debug_level >= 3:
                n_hot = int(np.count_nonzero(oob_node > force_tol))
                print(f"    OOB dist: max={unbalanced_force_ratio:.2e} "
                      f"p999={np.quantile(oob_node, 0.999):.2e} "
                      f"p99={np.quantile(oob_node, 0.99):.2e} "
                      f"p90={np.quantile(oob_node, 0.90):.2e} "
                      f"n>tol={n_hot}/{oob_node.size} ({100*n_hot/oob_node.size:.2f}%)")

            # Solve K * u_new = loads
            loads_free = loads[free_dofs]
            _tp = time.perf_counter() if _PROF_ON else None
            u_free_new = K_factor.solve(loads_free)
            if _PROF_ON:
                _prof_add("vp_trisolve", _tp)

            u_new = np.zeros(n_dof)
            u_new[free_dofs] = u_free_new

            # Convergence check: max|du| / max|u| < tolerance — Smith & Griffiths'
            # CHECON test (infinity norms), exactly as in p62.f90. The max-norm
            # matters: a small cluster of Gauss points sustaining benign localized
            # creep produces a bounded per-iteration du measured against the GLOBAL
            # maximum displacement, so the ratio decays and the stable state is
            # accepted within the iteration ceiling. A true failure mechanism
            # keeps feeding the global du and stays above tolerance. False
            # convergence from large failure
            # displacements is guarded by the max_disp_factor limit below.
            # Infinity norms over the FREE dofs only. Both u_new and u are exactly
            # zero on constrained dofs (u_new = np.zeros(n_dof) then
            # u_new[free_dofs] = solve; u inherits the same structure), so |u_new - u|
            # and |u_new| vanish there and the max over the free dofs equals the max
            # over all dofs — bit-identical, taken on the already-computed free
            # solution vector without materializing the full-length differences.
            # Measured FROM THE DATUM: with an equilibrated state carried in, the
            # in-situ displacement is not part of this solve's answer, and leaving it
            # in the norm would let a large fixed offset dilute both the CHECON ratio
            # and the hybrid criterion's displacement scale. u_datum_free is zero
            # without a carried state, so this is the same norm as before.
            _u_free_prev = u[free_dofs]
            norm_diff = np.max(np.abs(u_free_new - _u_free_prev))
            norm_u_new = np.max(np.abs(u_free_new - u_datum_free))

            if norm_u_new > 1e-30:
                relative_change = norm_diff / norm_u_new
            else:
                relative_change = norm_diff

            # Displacement history for the hybrid criterion. norm_u_new is max|u|,
            # already computed just above for the CHECON test, so this costs one
            # list append per _HYBRID_SAMPLE_EVERY iterations and nothing else. It
            # is recorded unconditionally: the metadata is reported on every
            # criterion, only the VERDICT's effect on the caller differs.
            if iteration % _HYBRID_SAMPLE_EVERY == 0:
                disp_hist.append(float(norm_u_new))
                oob_hist.append(float(unbalanced_force_ratio))

            # Force-equilibrium condition. The threshold is ABSOLUTE, which is what
            # makes the test immune to the size of the domain and to the size of the
            # yielding zone: the quantity is already dimensionless (a force over a
            # force, both at the same node), so it needs no reference drawn from the
            # run itself.
            #
            # This replaces an earlier PEAK-RELATIVE test — "has the rate of change
            # fallen to 1% of the largest rate this solve has seen?" — whose reference,
            # the first elastic-to-plastic burst, is an EXTENSIVE quantity: it grows
            # with the number of Gauss points that yield at once, and so with the size
            # of the mesh. The creep it was compared against is INTENSIVE, set by the
            # slope. Padding the domain inflated the reference, loosened the threshold,
            # and let a still-creeping slope be called settled, which read out as a
            # non-conservatively HIGH factor of safety.
            plastic_settled = unbalanced_force_ratio < force_tol

            # No-progress WATCH. When `_NO_PROGRESS_WINDOW` iterations pass with no
            # meaningful improvement (>1%) on the best out-of-balance value seen, the
            # residual is called PLATEAUED: the fact is recorded and reported, and the
            # solve carries on to its iteration cap. A plateau is an observation about
            # the residual, never a verdict on the slope.
            #
            # It used to END the solve and report the trial as failed, and that was
            # wrong in the direction that matters. Two measurements say so:
            #
            #   * THE RESIDUAL IS NOT MONOTONE. On the reinforced slope at F = 1.25
            #     (tri6, 1 ft elements) it sits near 2e-3 around iteration 9,500 —
            #     twice the tolerance it is chasing — climbs back above 1e-2, and only
            #     then falls through 1e-3 at iteration 16,242. "No improvement lately"
            #     is not evidence that the solve is finished.
            #   * THE WORK A TRIAL NEEDS GROWS WITH MESH REFINEMENT, A FIXED WINDOW
            #     DOES NOT. Iterations to equilibrium at that same F are 5,054 /
            #     10,555 / 16,242 at 2.5 / 1.5 / 1.0 ft element size, so a fixed
            #     1,500-iteration window is 30% of the required work on the coarse
            #     mesh and 9% on the fine one: refining the mesh tightened the
            #     guillotine without anyone asking it to. The bisection closed on the
            #     false failure and reported FS = 1.238 for a slope that reaches
            #     equilibrium at F = 1.25, and again at F = 1.424, on that same mesh.
            #     With the plateau demoted to an observation it reports 1.434, and the
            #     spread over a 6x range in element count falls from 18% to 5%.
            #
            # A near-tolerance guard ("do not stop while the residual is within 10 x
            # force_tol") was measured and does NOT repair it — the stop simply fires
            # later, on the rebound, at iteration 12,920. The window had to stop
            # deciding, not move.
            #
            # The cost falls on trials that genuinely fail: they now spend the whole
            # budget instead of leaving early. Converged trials are untouched, because
            # a converging solve leaves on convergence, which comes before the cap.
            if unbalanced_force_ratio < 0.99 * ufr_best:
                ufr_best = unbalanced_force_ratio
                last_progress_iter = iteration
            if (early_exit and plateau_iter is None and not plastic_settled
                    and iteration - last_progress_iter > _NO_PROGRESS_WINDOW):
                plateau_iter = iteration + 1
                plateau_ratio = float(unbalanced_force_ratio)
                if debug_level >= 1:
                    print(f"  Out-of-balance plateaued at iteration {iteration+1} "
                          f"({unbalanced_force_ratio:.2e} against tolerance "
                          f"{force_tol:.1e}; no >1% improvement for "
                          f"{_NO_PROGRESS_WINDOW} iterations) - recorded; the solve "
                          f"continues to its iteration cap")

            if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                print(f"  Iter {iteration+1:4d}: max|du|/max|u| = {relative_change:.3e}, "
                      f"max nodal OOB = {unbalanced_force_ratio:.3e}, "
                      f"yielding = {n_yielding}/{n_total_gp}, max|u| = {np.max(np.abs(u_new)):.6f}")

            # Report intra-solve progress (throttled) so a caller can advance a
            # progress bar within this viscoplastic solve, not just between solves.
            if progress_callback is not None and iteration % 10 == 0:
                try:
                    progress_callback((iteration + 1) / budget,
                                      f"vp iter {iteration + 1}/{budget}, "
                                      f"oob={unbalanced_force_ratio:.1e}")
                except Exception:
                    pass

            # Displacement limit check: detect false convergence from unbounded plastic flow
            # When VP displacements exceed a fraction of mesh height, the slope has physically
            # failed even if the relative convergence criterion is satisfied (see FLAC manual;
            # Griffiths & Lane 1999 displacement-vs-F plots).
            if vp_disp_limit is not None:
                # Plastic displacement = total, less the datum it is measured from and
                # less the elastic response to this stage's loads.
                u_vp = u_new - u_datum - u_elastic
                # Extract translational DOFs only for VP displacement check
                if dof_offset is not None:
                    vp_x = np.array([u_vp[dof_offset[nd]] for nd in range(n_nodes)])
                    vp_y = np.array([u_vp[dof_offset[nd] + 1] for nd in range(n_nodes)])
                else:
                    vp_x = u_vp[0::2]
                    vp_y = u_vp[1::2]
                max_vp_disp = float(np.max(np.sqrt(vp_x**2 + vp_y**2)))
                if max_vp_disp > vp_disp_limit:
                    converged = False
                    exit_reason = 'disp_limit'
                    u = u_new
                    if debug_level >= 1:
                        print(f"  Displacement limit exceeded at iteration {iteration+1}: "
                              f"max VP disp = {max_vp_disp:.2f} > limit {vp_disp_limit:.2f}")
                    break

            if relative_change < tolerance and plastic_settled:
                # --- post-peak softening fixed point -------------------------
                # The state is in equilibrium. NOW, and only now, ask which bars
                # have actually yielded: their elastic demand k*delta (forces_1d,
                # which is the UNCAPPED force) exceeds the capacity they were
                # allowed to carry. Those with a finite t_res drop to it, and we
                # keep iterating; shedding their load can push neighbours over,
                # so the process repeats until the softened set stops growing —
                # a genuine progressive-failure fixed point. It terminates: the
                # set only ever grows, and it is bounded by the element count.
                #
                # Deciding this on a converged state (rather than on whichever
                # iterate first overshot) is what makes the result independent of
                # the path the solver took to get here.
                if has_1d_elements and can_soften_1d.any():
                    demand = forces_1d
                    newly = (~softened_1d & can_soften_1d
                             & (demand > t_allow_by_1d_elem + 1e-9))
                    if newly.any() and n_soften_rounds < n_1d_elements:
                        softened_1d |= newly
                        n_soften_rounds += 1
                        if debug_level >= 1:
                            print(f"  Softening round {n_soften_rounds}: "
                                  f"{int(newly.sum())} bar element(s) dropped to "
                                  f"t_res ({int(softened_1d.sum())} total); "
                                  f"re-solving")
                        # reopen the iteration: reset the no-progress tracker so the
                        # new equilibrium is judged on its own decay history
                        ufr_best = np.inf
                        last_progress_iter = iteration
                        plateau_iter = None
                        plateau_ratio = None
                        u = u_new
                        continue
                # -------------------------------------------------------------
                converged = True
                exit_reason = 'converged'
                u = u_new
                if debug_level >= 1:
                    print(f"  Converged after {iteration+1} iterations "
                          f"(max|du|/max|u| = {relative_change:.3e}, "
                          f"max nodal OOB = {unbalanced_force_ratio:.2e} < {force_tol:.1e})")
                break

            # Early failure. Read on the sampled iterations only, from the two series
            # just appended above, and only AFTER the convergence test, so a solve
            # that settles on this very iteration always settles. A trial that fires
            # never reaches the budget check at the top of the loop, so it is never
            # considered for an extension — the two rules meet, but do not interact.
            if early_failure and iteration % _HYBRID_SAMPLE_EVERY == 0:
                _signal = _early_failure(disp_hist, oob_hist, u_elastic_scale)
                if _signal is not None:
                    converged = False
                    exit_reason = 'diverging'
                    diverging_iter = iteration
                    diverging_signal = _signal
                    u = u_new
                    if debug_level >= 1:
                        print(f"  Failing at iteration {iteration}: max|u| = "
                              f"{norm_u_new:.4g} "
                              f"({norm_u_new / u_elastic_scale:.1f}x elastic), "
                              f"out-of-balance {unbalanced_force_ratio:.2e} "
                              f"({_signal}) - the slope is running away; the trial "
                              f"is closed as FAILED without spending the rest of "
                              f"its budget")
                    break

            u = u_new


        total_iterations += iteration + 1
        if not converged:
            break   # stage failed -> overall failure

    if not converged and debug_level >= 1:
        print(f"  Did NOT converge after {total_iterations} iterations "
              f"(max|du|/max|u| = {relative_change:.3e}, exit {exit_reason})")

    # === Hybrid failure criterion: classify a non-converged trial ===
    # Computed on EVERY criterion so the metadata is always available for reporting
    # and for the A/B harness; only `stable` changes behavior, and only under
    # 'hybrid'. `converged` itself is never rewritten — a STABLE_STUCK trial did not
    # reach equilibrium and must not claim it.
    if converged:
        verdict, u_ratio, u_growth = 'CONVERGED', None, None
    else:
        verdict, u_ratio, u_growth = classify_nonconvergence(
            disp_hist, u_elastic_scale, exit_reason, model_height=mesh_height)
    stable = bool(converged or (failure_criterion == 'hybrid'
                                and verdict == 'STABLE_STUCK'))
    if not converged and debug_level >= 1:
        _ur = 'n/a' if u_ratio is None else f"{u_ratio:.2f}x elastic"
        _gr = 'n/a' if u_growth is None else f"{u_growth:+.3f}"
        print(f"  Displacement evidence: {_ur}, trailing growth {_gr} "
              f"-> verdict {verdict}"
              + ("  [HYBRID: treated as STABLE, not a failure]"
                 if stable else ""))

    # Copy grouped viscoplastic strains back into the per-element list used by
    # the post-processing blocks below.
    for grp in gp_groups:
        for k, (e, g) in enumerate(grp['pairs']):
            evp[e][g][:] = grp['evp'][k]

    # K0 initial stress, regrouped per (element, Gauss point) for the reporting pass
    # below. It carries the LAST stage's value, which is the state the reported
    # stresses belong to. None (the default) leaves the reporting arithmetic exactly
    # as it was.
    sig0_by_gp = None
    if sv0_gp is not None:
        sig0_by_gp = [[np.zeros(4) for _ in elem_gp_data[e]]
                      for e in range(n_elements)]
        for grp in gp_groups:
            for k, (e, g) in enumerate(grp['pairs']):
                sig0_by_gp[e][g] = grp['sig0'][k]

    # ---- Equilibrated state, for a later solve to start from (see _init_state) ----
    # The displacement field and the accumulated viscoplastic strain are the whole
    # state: with the K0 initial stress they reconstruct the stress field
    # sigma = sigma_0 + D (B u - evp) exactly, and a solve handed them starts from
    # this one's answer with its displacement measured from here.
    k0_state = None
    if sv0_gp is not None:
        k0_state = {"u": u.copy(),
                    "evp": [grp['evp'].copy() for grp in gp_groups],
                    "F": float(F), "converged": bool(converged)}

    # ---- Step 10: Compute final stresses, strains, plastic elements ----
    final_stresses = np.zeros((n_elements, 4))  # [sig_x, sig_y, tau_xy, sig_vm] compression-positive
    plastic_elements = np.zeros(n_elements, dtype=bool)
    yield_function_out = np.zeros(n_elements)

    for elem_idx in range(n_elements):
        gp_data_list = elem_gp_data[elem_idx]
        n_gp = len(gp_data_list)
        stress_avg_tp4 = np.zeros(4)

        for gp_idx, gp_data in enumerate(gp_data_list):
            B = gp_data['B']
            D4 = gp_data['D4']
            dof_idx = gp_data['dof_indices']

            u_elem = u[dof_idx]
            eps_total = B @ u_elem
            evp_gp = evp[elem_idx][gp_idx]
            eps_elastic4 = np.array([
                eps_total[0] - evp_gp[0],
                eps_total[1] - evp_gp[1],
                eps_total[2] - evp_gp[2],
                -evp_gp[3],
            ])
            stress_avg_tp4 += D4 @ eps_elastic4
            if sig0_by_gp is not None:
                stress_avg_tp4 += sig0_by_gp[elem_idx][gp_idx]

        stress_avg_tp4 /= n_gp
        u_elem_avg = sum(u_gp[elem_idx]) / len(u_gp[elem_idx]) if u_gp[elem_idx] else 0.0
        # D*B*u is effective; report total stresses for output (legacy
        # convention) and use the stresses as-is for the yield check.
        stress_total_tp4 = stress_avg_tp4 - np.array(
            [u_elem_avg, u_elem_avg, 0.0, u_elem_avg])
        sig_eff4 = stress_avg_tp4
        # compression-positive in-plane components for output; sig_vm = dsbar
        sig_x, sig_y = -stress_total_tp4[0], -stress_total_tp4[1]
        tau_xy = stress_total_tp4[2]
        _, sig_vm, _ = stress_invariants(stress_total_tp4)
        final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]

        sigm, dsbar, theta = stress_invariants(sig_eff4)
        _c_rep, _phi_rep = c_reduced[elem_idx], phi_reduced[elem_idx]
        if pow_flag_by_elem[elem_idx]:
            # reporting tangent from the final effective stress state
            _sn = max(-(sig_eff4[0] + sig_eff4[1]) * 0.5, 0.0)
            _sef = max(_sn + fem_data["pow_d_by_elem"][elem_idx],
                       1e-4 * max(1.0, _sn))
            _a = fem_data["pow_a_by_elem"][elem_idx]
            _b = fem_data["pow_b_by_elem"][elem_idx]
            _Fe = F_by_elem[elem_idx]           # 1.0 where SSR-excluded
            _sl = _a * _b * _sef ** (_b - 1.0) / _Fe
            _c_rep = (_a * _sef ** _b + fem_data["pow_cp_by_elem"][elem_idx]) / _Fe \
                     - _sn * _sl
            _phi_rep = np.arctan(_sl)
        elif hb_flag_by_elem[elem_idx]:
            # Reporting tangent from the final effective stress state. This MUST use the
            # same failure-plane abscissa the viscoplastic loop linearized on, or the
            # reported yield function describes a different envelope than the one that
            # was actually solved. The VP loop's converged (c_r, phi) are not carried out
            # here, so recover the abscissa by iterating the same fixed point to closure.
            _s_prime = max(-(sig_eff4[0] + sig_eff4[1]) * 0.5, 0.0)
            _sci = fem_data["hb_sci_by_elem"][elem_idx]
            _mb = fem_data["hb_mb_by_elem"][elem_idx]
            _s = fem_data["hb_s_by_elem"][elem_idx]
            _a_hb = fem_data["hb_a_by_elem"][elem_idx]
            _c_rep, _phi_rep = 0.0, 0.0
            _sn = _s_prime                      # seed at the circle center
            _Fe = F_by_elem[elem_idx]           # 1.0 where SSR-excluded
            for _ in range(40):
                _ci, _phii = hb_tangent_const(_sn, _sci, _mb, _s, _a_hb, iters=40)
                _c_rep = float(_ci) / _Fe
                _phi_rep = np.arctan(np.tan(np.radians(float(_phii))) / _Fe)
                _sn_new = max(_s_prime * np.cos(_phi_rep) ** 2
                              - _c_rep * np.sin(_phi_rep) * np.cos(_phi_rep), 0.0)
                if abs(_sn_new - _sn) <= 1e-9 * max(1.0, abs(_sn)):
                    _sn = _sn_new
                    break
                _sn = _sn_new
        # Add the matric-suction apparent cohesion (reduced by F) to the reported
        # tangent cohesion so the reported yield function matches the envelope the
        # VP loop actually solved. Gated on suction_active -> default runs untouched.
        if (suction_active and u_gp_signed is not None
                and suction_tanphib_by_elem[elem_idx] > 0.0):
            _s_suc = np.minimum(np.maximum(-np.asarray(u_gp_signed[elem_idx]), 0.0),
                                suction_scap_by_elem[elem_idx])
            _c_rep = _c_rep + (suction_tanphib_by_elem[elem_idx]
                               * float(_s_suc.mean()) / F_by_elem[elem_idx])
        f_yield = mc_yield_invariants(sigm, dsbar, theta, _c_rep, _phi_rep)
        yield_function_out[elem_idx] = f_yield
        plastic_elements[elem_idx] = f_yield > 1e-8

    # Pure-elastic elements never entered the plastic loop, so they carry no
    # viscoplastic strain; but their linear elastic stress can still lie outside
    # the Mohr-Coulomb envelope, which the per-element yield check above would read
    # as "plastic". Force the failure flag off — RS2's "Plasticity: None" materials
    # never show as yielded. No-op when elastic_mask is None (bit-identical).
    if elastic_by_elem is not None:
        plastic_elements[elastic_by_elem] = False

    strains = compute_strains(nodes, elements, element_types, u, dof_offset=dof_offset)

    # Compute viscoplastic max shear strain per element (averaged over Gauss points)
    # This is what Griffiths plots — zero in elastic regions, large in failure zone
    vp_shear_strain = np.zeros(n_elements)
    for elem_idx in range(n_elements):
        n_gp = len(evp[elem_idx])
        gp_shear = 0.0
        for gp_idx in range(n_gp):
            evp_gp = evp[elem_idx][gp_idx]
            eps_x, eps_y, gamma_xy = evp_gp[0], evp_gp[1], evp_gp[2]
            gp_shear += sqrt(((eps_x - eps_y) / 2)**2 + (gamma_xy / 2)**2)
        vp_shear_strain[elem_idx] = gp_shear / n_gp

    # ---- Step 10b: Compute final 1D truss element forces ----
    if has_1d_elements:
        for elem_idx_1d in range(n_1d_elements):
            dof_idx = dof_indices_1d[elem_idx_1d][:n_dof_1d[elem_idx_1d]]
            k = k_by_1d_elem[elem_idx_1d]
            cos_t = cos_theta_1d[elem_idx_1d]
            sin_t = sin_theta_1d[elem_idx_1d]

            u_elem = u[dof_idx]
            du_x = u_elem[2] - u_elem[0]
            du_y = u_elem[3] - u_elem[1]
            delta = du_x * cos_t + du_y * sin_t
            T = k * delta

            # Report the force the bar actually delivers, under the same law that
            # was enforced in the viscoplastic loop: tension-only, yielding at the
            # end-ramp envelope t_allow — or at t_res if the bar ended up in
            # the softened set.
            cap = (t_res_by_1d_elem[elem_idx_1d] if softened_1d[elem_idx_1d]
                   else t_allow_by_1d_elem[elem_idx_1d])
            forces_1d[elem_idx_1d] = min(max(T, 0.0), cap)

    # ---- Step 10c: Final pile beam element actions, under the same capacity law
    # the loop enforced ----
    #
    # Read through _pile_element_capacity, the one place the law is written, so the
    # reported actions are the ones the converged equilibrium carries. Reading them
    # any other way -- computing the elastic actions and then CLIPPING them to the
    # capacity -- reports a cap the equilibrium never enforced, which is how a moment
    # capacity that did nothing at all could read as an enforced one.
    if has_pile_elements:
        for p_idx in range(n_pile_elements):
            n_dof_e = int(n_dof_pile[p_idx])
            n_node_e = n_dof_e // 3
            dof_idx = dof_indices_pile[p_idx][:n_dof_e]
            cos_t = cos_theta_pile[p_idx]
            sin_t = sin_theta_pile[p_idx]
            L = L_pile_elem[p_idx]
            EI_val = EI_pile[p_idx]
            EA_val = EA_pile[p_idx]
            Kl = K_local_pile[p_idx] if K_local_pile is not None else None

            u_local = _beam_rotation(cos_t, sin_t, n_node_e) @ u[dof_idx]
            (T_force, V, M1, M2, _corr, p_rot,
             yielded_V_e, yielded_M_e) = _pile_element_capacity(
                u_local, cos_t, sin_t, L, EA_val, EI_val, Kl, n_node_e,
                V_cap_pile[p_idx], M_cap_pile[p_idx],
                hinge_ends=pile_hinge_ends[p_idx])

            forces_pile_axial[p_idx] = T_force
            forces_pile_lateral[p_idx] = V
            forces_pile_moment[p_idx] = [M1, M2]
            plastic_rot_pile[p_idx] = p_rot
            if yielded_V_e:
                yielded_pile_V[p_idx] = True
            if yielded_M_e:
                yielded_pile_M[p_idx] = True

    # Reported displacement is measured from the datum (zero without a carried
    # state). Stresses, strains and every structural force above are functions of the
    # ABSOLUTE displacement and are computed from it; only the displacement field
    # itself is re-referenced, so that a K0 run reports the deformation caused by
    # this solve rather than the fictitious travel of setting up the in-situ state.
    u_reported = u if _init_u is None else u - u_datum

    n_plastic = np.sum(plastic_elements)
    if debug_level >= 1:
        print(f"  Plastic elements: {n_plastic}/{n_elements}")
        print(f"  Max displacement: {np.max(np.abs(u_reported)):.6f}")
        print(f"  Max VP shear strain: {np.max(vp_shear_strain):.6e}")
        print(f"  Unbalanced force ratio: {unbalanced_force_ratio:.3e}")
        if has_1d_elements:
            max_force = np.max(forces_1d) if n_1d_elements > 0 else 0.0
            n_failed = np.sum(failed_1d)
            n_active = np.sum(forces_1d > 0)
            print(f"  1D elements: {n_active} active, {n_failed} failed, "
                  f"max force = {max_force:.2f}")
        if has_pile_elements:
            max_axial = np.max(np.abs(forces_pile_axial))
            max_lateral = np.max(np.abs(forces_pile_lateral))
            max_moment = np.max(np.abs(forces_pile_moment))
            print(f"  Pile elements: {n_pile_elements}, "
                  f"max axial = {max_axial:.2f}, max shear = {max_lateral:.2f}, "
                  f"max moment = {max_moment:.2f}")

    return {
        "converged": converged,
        # Hybrid-criterion metadata (present on every criterion; see
        # classify_nonconvergence). `stable` is what a bisection should read.
        "stable": stable,
        "verdict": verdict,
        "u_ratio": u_ratio,
        "u_growth": u_growth,
        "u_elastic_scale": u_elastic_scale,
        "exit_reason": exit_reason,
        # Equilibrated initial-stress state (K0 runs only; None otherwise). Internal:
        # solve_ssrm's equilibration solve hands this to every trial as _init_state.
        "_k0_state": k0_state,
        # The no-progress watch: the iteration at which the out-of-balance stopped
        # improving, and the value it stalled at, or None if it never stalled. A
        # plateau does not end the solve — a trial that plateaus and then runs out of
        # budget reports exit_reason 'iteration_cap' with these fields set, which is
        # what tells a reader the residual had stopped moving before the cap.
        "plateau_iteration": plateau_iter,
        "plateau_ratio": plateau_ratio,
        # The early-failure rule: the iteration at which it fired and which of its
        # two tests fired, or None on any other exit (see _early_failure).
        "diverging_iteration": diverging_iter,
        "diverging_signal": diverging_signal,
        # Retained name: True when a plateau was seen and the solve carried on anyway.
        "early_exit_suppressed": plateau_iter is not None,
        # Budget extension: how many extra chunks this solve was granted because the
        # solve was still progressing at the budget, and the budget it ended on.
        "budget_extensions": n_extensions,
        "iteration_budget": budget,
        "failure_criterion": failure_criterion,
        "iterations": total_iterations,
        "displacements": u_reported,
        "displacements_elastic": u_elastic,
        "stresses": final_stresses,
        "strains": strains,
        "vp_shear_strain": vp_shear_strain,
        "plastic_elements": plastic_elements,
        "yield_function": yield_function_out,
        "max_displacement": np.max(np.abs(u_reported)),
        "plastic_strains": {i: np.array(evp[i]) for i in range(n_elements)},
        "algorithm": "Griffiths & Lane (1999) Viscoplastic",
        "F": F,
        "residual": relative_change if 'relative_change' in locals() else 0.0,
        "unbalanced_force_ratio": unbalanced_force_ratio,
        "plastic_fraction": n_plastic / n_elements if n_elements > 0 else 0.0,
        "forces_1d": forces_1d if has_1d_elements else np.array([]),
        "failed_1d_elements": failed_1d if has_1d_elements else np.array([], dtype=bool),
        # bars that dropped to their residual capacity (converged-state fixed point)
        "softened_1d_elements": softened_1d if has_1d_elements else np.array([], dtype=bool),
        "forces_pile_axial": forces_pile_axial if has_pile_elements else np.array([]),
        "forces_pile_lateral": forces_pile_lateral if has_pile_elements else np.array([]),
        "forces_pile_moment": forces_pile_moment if has_pile_elements else np.zeros((0, 2)),
        "pile_plastic_rotation": (plastic_rot_pile if has_pile_elements
                                  else np.zeros((0, 2))),
        "yielded_pile_V": yielded_pile_V if has_pile_elements else np.array([], dtype=bool),
        "yielded_pile_M": yielded_pile_M if has_pile_elements else np.array([], dtype=bool),
        "yielded_pile": (yielded_pile_V | yielded_pile_M) if has_pile_elements else np.array([], dtype=bool),
    }


# ===================== Newton-Raphson Mohr-Coulomb path (SPIKE) ==============
#
# An ALTERNATIVE driver for a single strength-reduction trial, selected by the
# internal `fem_solver` switch ('viscoplastic', the default and the definition of
# every locked result, or 'newton'). Nothing below runs unless the switch is
# thrown, and the switch is read in exactly one place — the hook at the top of
# solve_fem — so the viscoplastic path is untouched.
#
# The scheme is the textbook one: an implicit backward-Euler closest-point return
# to the Mohr-Coulomb cone in principal-stress space, its consistent (algorithmic)
# tangent, and a global Newton-Raphson iteration on the equilibrium residual under
# incremental (adaptively sub-stepped) gravity loading. Non-associated flow with
# psi = 0 makes the tangent non-symmetric, so the sparse solve is a general LU.
#
# Reference for the return mapping and its corner/apex cases: de Souza Neto,
# Peric & Owen, "Computational Methods for Plasticity" (2008), Section 8.2 —
# a two-vector return at the sextant corners and a return to the tensile apex.

_FEM_SOLVER_ENV = "XSLOPE_FEM_SOLVER"
# The environment override announces itself once per process. A stale shell
# variable would otherwise silently recompute every factor of safety in a session
# on the non-default driver, and the run would look exactly like a default run.
_FEM_SOLVER_ENV_ANNOUNCED = False


def resolve_fem_solver(fem_solver=None):
    """Which per-trial driver runs: 'viscoplastic' (default) or 'newton'.

    ``None`` means UNSPECIFIED and falls through to the environment variable
    ``XSLOPE_FEM_SOLVER``, then to 'viscoplastic'. The environment hook exists so
    a whole benchmark run can be flipped without touching a call site; an
    explicit argument always wins over it.

    When the ENVIRONMENT — and not an explicit argument — is what selects a
    non-default driver, one warning line is printed per process. The explicit
    argument path stays silent: a caller that named the solver knows which one it
    asked for, while a shell variable left over from an earlier session does not
    announce itself any other way.
    """
    global _FEM_SOLVER_ENV_ANNOUNCED
    from_env = False
    if fem_solver is None:
        env_value = os.environ.get(_FEM_SOLVER_ENV)
        from_env = bool(env_value)
        fem_solver = env_value or "viscoplastic"
    key = str(fem_solver).strip().lower()
    if key in ("vp", "viscoplastic", "griffiths", "default"):
        return "viscoplastic"
    if key in ("nr", "newton", "newton-raphson", "newton_raphson"):
        if from_env and not _FEM_SOLVER_ENV_ANNOUNCED:
            _FEM_SOLVER_ENV_ANNOUNCED = True
            print(f"\n*** {_FEM_SOLVER_ENV}={fem_solver!r} is set in the environment: "
                  f"every FEM/SSRM solve in this process runs on the NON-DEFAULT "
                  f"'newton' driver, and its factors of safety are NOT the locked "
                  f"viscoplastic ones. Unset {_FEM_SOLVER_ENV} to restore the "
                  f"default. ***\n")
        return "newton"
    raise ValueError(
        f"Unknown fem_solver {fem_solver!r}. Supported: 'viscoplastic' (default) "
        "and 'newton'.")


# Branch codes recorded per Gauss point by the return map, for reporting only.
_NR_ELASTIC, _NR_MAIN, _NR_RIGHT, _NR_LEFT, _NR_APEX = 0, 1, 2, 3, 4
# ... and the Rankine tension-cutoff branches, reached only where a Gauss point
# carries a FINITE tensile cap. The naming says which surfaces are active:
# TENS = the cap alone on sigma_1; TENS2 = on sigma_1 and sigma_2; TENS3 = the
# hydrostatic-tension return to (T, T, T), which is what replaces the Mohr-Coulomb
# apex whenever T sits below it; MC_TENS = the Mohr-Coulomb main plane AND the cap,
# which is the intersection edge; CORNER_TENS = a sextant corner and the cap;
# MC_TENS2 = the main plane with the cap on two principal stresses.
(_NR_TENS, _NR_TENS2, _NR_TENS3, _NR_MC_TENS, _NR_CORNER_TENS,
 _NR_MC_TENS2, _NR_TENS_FALLBACK) = 5, 6, 7, 8, 9, 10, 11

_NR_BRANCH_NAMES = {
    _NR_ELASTIC: "elastic", _NR_MAIN: "main", _NR_RIGHT: "right_corner",
    _NR_LEFT: "left_corner", _NR_APEX: "apex",
    _NR_TENS: "tension", _NR_TENS2: "tension_2", _NR_TENS3: "tension_apex",
    _NR_MC_TENS: "mc_tension_edge", _NR_CORNER_TENS: "corner_tension",
    _NR_MC_TENS2: "mc_tension_2", _NR_TENS_FALLBACK: "tension_fallback",
}

# The active sets the combined Mohr-Coulomb / Rankine return tries, in order, with
# the branch code each one records. A surface index is
#   0,1,2 -> the Mohr-Coulomb planes on the ordered pairs (1,3), (1,2), (2,3)
#   3,4,5 -> the Rankine planes sigma_1 <= T, sigma_2 <= T, sigma_3 <= T
#
# Two facts about the ordered sextant sigma_1 >= sigma_2 >= sigma_3 cut the
# combinatorics down to these seven, and both are exact rather than heuristic:
#
#   * f(1,3) >= f(1,2) and f(1,3) >= f(2,3) identically there, so plane (1,3) is in
#     the Mohr-Coulomb active set whenever any Mohr-Coulomb plane is. The only
#     Mohr-Coulomb sets are {}, {13}, {13,12} and {13,23} — which are exactly the
#     main plane and the two corners the plain map already carries.
#   * a RETURNED state must be ordered, so if sigma_2 is at the cap then sigma_1 is
#     too. The only Rankine sets are the prefixes {1}, {1,2} and {1,2,3}.
#
# Three principal stresses admit at most three independent active constraints, which
# rules out the corner-plus-two-cap and main-plus-three-cap combinations.
_NR_RANKINE_SETS = (
    ((3,), _NR_TENS),
    ((0, 3), _NR_MC_TENS),
    ((0, 2, 3), _NR_CORNER_TENS),
    ((0, 1, 3), _NR_CORNER_TENS),
    ((3, 4), _NR_TENS2),
    ((0, 3, 4), _NR_MC_TENS2),
    ((3, 4, 5), _NR_TENS3),
)


def _nr_surface(k, A, Bc, ccp, mu, lam, T):
    """One yield surface's normal, its elastic return vector, and its intercept.

    Every surface here is LINEAR in the ordered principal stress — ``f_k = n_k . s
    - r_k`` — with a CONSTANT flow direction, which is what makes a return for a
    given active set one small linear solve and exact rather than iterated.

    ``v_k = D m_k`` is the stress change per unit multiplier. For a Mohr-Coulomb
    plane the psi = 0 flow ``m = (e_i - e_j)/2`` is deviatoric and the isotropic
    operator maps it to ``mu (e_i - e_j)``, so that return cannot move the mean
    stress. For a Rankine plane the flow is ASSOCIATED, ``m = e_k``, and
    ``D m = 2 mu e_k + lambda 1`` DOES move it — which is the entire reason the
    second surface exists, since a psi = 0 Mohr-Coulomb flow cannot relieve a
    tensile mean stress at all.
    """
    z = np.zeros_like(A)
    o = np.ones_like(A)
    if k == 0:                                   # Mohr-Coulomb on (sigma_1, sigma_3)
        return (np.stack([A, z, -Bc], 1),
                mu[:, None] * np.stack([o, z, -o], 1), ccp)
    if k == 1:                                   # ... on (sigma_1, sigma_2)
        return (np.stack([A, -Bc, z], 1),
                mu[:, None] * np.stack([o, -o, z], 1), ccp)
    if k == 2:                                   # ... on (sigma_2, sigma_3)
        return (np.stack([z, A, -Bc], 1),
                mu[:, None] * np.stack([z, o, -o], 1), ccp)
    e = np.zeros((len(A), 3))
    e[:, k - 3] = 1.0                            # Rankine on sigma_(k-3)
    return e, 2.0 * mu[:, None] * e + lam[:, None], T


def _nr_solve_small(M, rhs):
    """Solve a stack of small dense systems, flagging the singular ones.

    ``M`` is (N, m, m) with m at most 3. A singular member is a degenerate active
    set — two surfaces whose return vectors are parallel — and it is reported
    rather than solved, so the caller moves on to the next candidate instead of
    taking a nonsense multiplier.
    """
    m = M.shape[1]
    det = np.linalg.det(M)
    scale = np.max(np.abs(M), axis=(1, 2)) ** m
    ok = np.abs(det) > 1e-12 * np.maximum(scale, 1e-300)
    Ms = np.where(ok[:, None, None], M, np.eye(m))
    # rhs is stacked explicitly: NumPy 2 reads a (N, m) right-hand side as one
    # m-by-N matrix rather than as N vectors.
    return np.linalg.solve(Ms, rhs[:, :, None])[:, :, 0], ok


def _nr_rankine_return(s, A, Bc, ccp, mu, lam, T):
    """Return ordered principal stresses to Mohr-Coulomb INTERSECTED with the cap.

    ``s`` is (N, 3), already ordered sigma_1 >= sigma_2 >= sigma_3, tension-positive,
    and every other argument is per-point. Returns ``(sigma, branch, unresolved)``.

    This is Koiter's rule done as an ACTIVE-SET SEARCH rather than a case tree: each
    candidate set in ``_NR_RANKINE_SETS`` is returned exactly by one linear solve,
    and a candidate is accepted only if it is CONSISTENT — every multiplier
    non-negative, the returned state still ordered, and BOTH surfaces satisfied. The
    first consistent candidate wins. Nothing here decides a region by inspecting the
    trial state, so a region the design got wrong shows up as an inconsistent return
    rather than as a wrong answer that looks right.

    Two passes over the candidates: the first at a tolerance three orders above
    round-off, the second an order below the branch structure, which is what a point
    sitting exactly ON a region boundary needs (there the two candidates agree to
    round-off, so accepting either is the same answer). Anything still unresolved is
    reported through ``unresolved`` and the caller decides; it is expected to be
    empty and the fuzz asserts that it is.
    """
    n = len(s)
    out = s.copy()
    branch = np.zeros(n, dtype=np.int8)
    todo = np.ones(n, dtype=bool)
    sc = np.maximum.reduce([np.abs(s[:, 0]), np.abs(s[:, 2]), np.abs(T),
                            np.abs(ccp), np.ones(n)])
    for rel_tol in (1e-13, 1e-9):
        for surf, code in _NR_RANKINE_SETS:
            if not todo.any():
                break
            j = np.flatnonzero(todo)
            parts = [_nr_surface(k, A[j], Bc[j], ccp[j], mu[j], lam[j], T[j])
                     for k in surf]
            nmat = np.stack([p[0] for p in parts], axis=1)      # (n, m, 3)
            vmat = np.stack([p[1] for p in parts], axis=1)      # (n, m, 3)
            rvec = np.stack([p[2] for p in parts], axis=1)      # (n, m)
            M = np.einsum('nia,nja->nij', nmat, vmat)
            rhs = np.einsum('nia,na->ni', nmat, s[j]) - rvec
            gam, ok = _nr_solve_small(M, rhs)
            sig = s[j] - np.einsum('ni,nia->na', gam, vmat)
            tol = rel_tol * sc[j]
            # A multiplier is checked through the stress it moves, so the test is in
            # one unit whatever the surface's normalization.
            ok &= np.all(gam * np.linalg.norm(vmat, axis=2)
                         >= -tol[:, None], axis=1)
            a1, a2, a3 = sig[:, 0], sig[:, 1], sig[:, 2]
            ok &= (a1 >= a2 - tol) & (a2 >= a3 - tol)
            ok &= (A[j] * a1 - Bc[j] * a3 - ccp[j]) <= tol
            ok &= a1 <= T[j] + tol
            if ok.any():
                sel = j[ok]
                out[sel] = sig[ok]
                branch[sel] = code
                todo[sel] = False
        if not todo.any():
            break
    return out, branch, todo


def mc_return_map(sig_tr4, c, snph, csph, mu, t_cap=None, lam=None):
    """Closest-point return to the Mohr-Coulomb cone, psi = 0, plane strain.

    ``sig_tr4`` is the (G, 4) elastic trial stress [sx, sy, txy, sz],
    TENSION-POSITIVE, in the component order build_constitutive_matrix_4 uses.
    ``c``, ``snph``, ``csph``, ``mu`` are per-Gauss-point arrays of length G
    (cohesion, sin/cos of the friction angle, shear modulus).

    ``t_cap`` is the per-point Rankine tensile cap T (``inf`` = none) and ``lam``
    the Lame constant the associated Rankine flow needs; both ``None`` — the
    default, and every call on a model that sets no ``t_cut`` — skips that code
    entirely, so the plain Mohr-Coulomb path is untouched arithmetic for
    untouched arithmetic. See :func:`_nr_rankine_return`.

    Returns ``(sig4, branch)`` — the returned stress and a per-point branch code
    (0 elastic, 1 main plane, 2 right corner, 3 left corner, 4 apex, and 5-11 for
    the branches where the tensile cap is active; see ``_NR_BRANCH_NAMES``).

    On the ordered principal stresses sigma_1 >= sigma_2 >= sigma_3
    (tension-positive) the yield function is

        f = A*sigma_1 - Bc*sigma_3 - c*cos(phi),
        A = (1 + sin phi)/2,   Bc = (1 - sin phi)/2,

    and the psi = 0 flow potential is g = (sigma_1 - sigma_3)/2. The isotropic
    elastic operator maps that flow direction to 2*mu times itself, so the
    main-plane multiplier is exactly f/mu and each corner is a 2x2 linear solve:
    the return is LINEAR in the trial stress on every branch, which is what makes
    the algorithmic tangent below exact on that branch.
    """
    sx, sy, txy, sz = sig_tr4[:, 0], sig_tr4[:, 1], sig_tr4[:, 2], sig_tr4[:, 3]
    ctr = 0.5 * (sx + sy)
    half = 0.5 * (sx - sy)
    R = np.sqrt(half * half + txy * txy)

    P = np.stack([ctr + R, ctr - R, sz], axis=1)
    order = np.argsort(-P, axis=1, kind="stable")
    Ps = np.take_along_axis(P, order, axis=1)
    s1, s2, s3 = Ps[:, 0], Ps[:, 1], Ps[:, 2]

    A = 0.5 * (1.0 + snph)
    Bc = 0.5 * (1.0 - snph)
    ccp = c * csph
    f = A * s1 - Bc * s3 - ccp

    branch = np.zeros(len(f), dtype=np.int8)
    yid = f > 0.0
    n1, n2, n3 = s1.copy(), s2.copy(), s3.copy()

    if np.any(yid):
        # --- main plane -------------------------------------------------------
        gm = f                                   # = mu * dgamma
        m1 = s1 - gm
        m3 = s3 + gm
        # An ordering violation is what sends a point to a corner. The slack is
        # relative to the local stress scale, so a point sitting exactly on a
        # sextant boundary is not thrown to a corner by rounding.
        scale = np.maximum(np.abs(s1) + np.abs(s3), 1.0) * 1e-12
        ok = (m1 >= s2 - scale) & (s2 >= m3 - scale)
        take_main = yid & ok
        n1 = np.where(take_main, m1, n1)
        n3 = np.where(take_main, m3, n3)
        branch[take_main] = _NR_MAIN

        # --- corners ----------------------------------------------------------
        bad = yid & ~ok
        right = np.zeros_like(bad)
        if np.any(bad):
            fa = f
            # right corner (sigma_1 collides with sigma_2): planes (1,3) and (2,3)
            right = bad & (m1 < s2 - scale)
            if np.any(right):
                fb = A * s2 - Bc * s3 - ccp
                det = 1.0 - Bc * Bc
                ga = (fa - Bc * fb) / det          # = mu * dgamma_a
                gb = (fb - Bc * fa) / det
                n1 = np.where(right, s1 - ga, n1)
                n2 = np.where(right, s2 - gb, n2)
                n3 = np.where(right, s3 + ga + gb, n3)
                branch[right] = _NR_RIGHT
            # left corner (sigma_2 collides with sigma_3): planes (1,3) and (1,2)
            left = bad & ~right
            if np.any(left):
                fb = A * s1 - Bc * s2 - ccp
                det = 1.0 - A * A
                ga = (fa - A * fb) / det
                gb = (fb - A * fa) / det
                n1 = np.where(left, s1 - ga - gb, n1)
                n2 = np.where(left, s2 + gb, n2)
                n3 = np.where(left, s3 + ga, n3)
                branch[left] = _NR_LEFT

        # --- apex -------------------------------------------------------------
        # With psi = 0 every return direction above is purely deviatoric, so none
        # of them changes the mean stress. A trial state whose mean stress already
        # lies beyond the cone's tensile apex therefore cannot be returned to the
        # cone at all by a deviatoric correction, and the only admissible point is
        # the apex itself, sigma_1 = sigma_2 = sigma_3 = c*cot(phi). This is the
        # case the viscoplastic scheme cannot resolve: its psi = 0 flow chases a
        # yield function that never closes, which is why that path needs a
        # separate Rankine surface to hold tensile states.
        snph_safe = np.where(snph > 1e-12, snph, 1.0)
        s_apex = np.where(snph > 1e-12, ccp / snph_safe, np.inf)
        p_tr = (s1 + s2 + s3) / 3.0
        apex = yid & (p_tr >= s_apex)
        # A corner solve whose result is still inadmissible lands on the apex too.
        apex = apex | (yid & ((n1 < n2 - scale) | (n2 < n3 - scale)))
        if np.any(apex):
            n1 = np.where(apex, s_apex, n1)
            n2 = np.where(apex, s_apex, n2)
            n3 = np.where(apex, s_apex, n3)
            branch[apex] = _NR_APEX

    # --- the Rankine tension cutoff, the SECOND surface -----------------------
    # The cap is a separate yield surface from the Mohr-Coulomb one and the two
    # combine by Koiter's rule, which is what the viscoplastic path does when it
    # SUMS the two flows into evpg. Here the same physics is a multi-surface
    # return: the Mohr-Coulomb-only state above is a valid answer for any point
    # whose returned sigma_1 is already under its cap (the Rankine multiplier is
    # zero there), so only the rest are reworked. That is what keeps a model with
    # no cap on exactly the arithmetic it was on before this surface existed, and
    # it is also why an inert cap — one at or above the Mohr-Coulomb apex, which no
    # admissible state can reach — costs a comparison and nothing else.
    if t_cap is not None:
        T_ = np.broadcast_to(np.asarray(t_cap, dtype=float), f.shape)
        over = np.isfinite(T_) & (n1 > T_)
        if np.any(over):
            if lam is None:
                raise ValueError("mc_return_map needs `lam` when `t_cap` is set: "
                                 "the Rankine flow is associated, so its return "
                                 "vector is 2 mu e_k + lambda 1.")
            k = np.flatnonzero(over)
            A_ = np.broadcast_to(A, f.shape)
            Bc_ = np.broadcast_to(Bc, f.shape)
            ccp_ = np.broadcast_to(ccp, f.shape)
            mu_ = np.broadcast_to(np.asarray(mu, dtype=float), f.shape)
            lam_ = np.broadcast_to(np.asarray(lam, dtype=float), f.shape)
            # From the TRIAL state, not from the Mohr-Coulomb-returned one. Koiter's
            # rule is SIMULTANEOUS — sigma = sigma_trial - sum_k gamma_k D m_k over
            # the whole active set — and returning to one surface and then to the
            # other is a different, wrong operation. It is wrong in a way that
            # hides: a state already returned to the Mohr-Coulomb surface and then
            # capped comes back INSIDE that surface (the Rankine return lowers f by
            # gamma*(2*mu*A + lambda*sin(phi)), which is non-negative at every
            # friction angle), so every point looks admissible and the intersection
            # edge never executes. The branch histogram is what says so.
            sig_r, br_r, todo = _nr_rankine_return(
                np.stack([s1[k], s2[k], s3[k]], axis=1),
                A_[k], Bc_[k], ccp_[k], mu_[k], lam_[k], T_[k])
            if np.any(todo):
                # No candidate active set was consistent. The isotropic point at
                # min(T, apex) satisfies BOTH surfaces, so the state stays
                # admissible and the verdict's yield reading stays honest; the
                # branch code says it happened, and the fuzz asserts it does not.
                snph_s = np.where(snph > 1e-12, snph, 1.0)
                apex_ = np.where(snph > 1e-12, ccp / snph_s, np.inf)
                fb = np.minimum(np.broadcast_to(apex_, f.shape)[k][todo],
                                T_[k][todo])
                sig_r[todo] = fb[:, None]
                br_r[todo] = _NR_TENS_FALLBACK
            n1[k], n2[k], n3[k] = sig_r[:, 0], sig_r[:, 1], sig_r[:, 2]
            branch[k] = br_r

    # --- back to component form ----------------------------------------------
    Pn = np.empty_like(P)
    np.put_along_axis(Pn, order, np.stack([n1, n2, n3], axis=1), axis=1)
    sa_n, sb_n, sz_n = Pn[:, 0], Pn[:, 1], Pn[:, 2]
    ctr_n = 0.5 * (sa_n + sb_n)
    rad_n = 0.5 * (sa_n - sb_n)
    # Mohr-Coulomb flow is coaxial: the principal DIRECTIONS are the trial state's
    # and only the principal VALUES move. Below a degenerate radius the in-plane
    # state is isotropic and the direction is arbitrary; cos(2t) = 1, sin(2t) = 0
    # is then the correct — and continuous — choice.
    big = R > 1e-14 * np.maximum(np.abs(ctr), 1.0)
    Rs = np.where(big, R, 1.0)
    c2 = np.where(big, half / Rs, 1.0)
    s2t = np.where(big, txy / Rs, 0.0)

    out = np.empty_like(sig_tr4)
    out[:, 0] = ctr_n + rad_n * c2
    out[:, 1] = ctr_n - rad_n * c2
    out[:, 2] = rad_n * s2t
    out[:, 3] = sz_n
    # An untouched point is restored EXACTLY rather than reconstructed, so a stress
    # the return map did not act on carries no round-off from the round trip
    # through the principal frame. Read off the branch code, which is zero on
    # exactly the points no surface returned — including a point that yielded
    # neither the shear surface nor the cap.
    el = branch == _NR_ELASTIC
    if np.any(el):
        out[el] = sig_tr4[el]
    return out, branch


# ---- curved strength envelopes: Hoek-Brown and the power curve ---------------
# Both are carried here the way the viscoplastic driver carries them: as a
# per-Gauss-point MOHR-COULOMB TANGENT to the F-reduced envelope, re-derived from
# the current stress, with the ordinary psi = 0 return map running on that tangent.
# Neither driver has a curved yield surface anywhere — solve_fem's Step 6
# overwrites grp['c_r'] / ['snph'] / ['csph'] in place for the flagged Gauss points
# and then runs the same MC yield function and the same MOCOUQ flow on all points
# alike. Reproducing THAT is what makes the two drivers answer the same question;
# an exact curved-surface return would converge somewhere else. See SPIKE.md,
# "CURVED ENVELOPES".
#
# The one thing this path must do that the viscoplastic path need not: the
# linearization has to be a pure function of the displacement, or the residual is
# not one and neither the line search nor the convergence test means anything. So
# the tangent is closed as a SELF-CONSISTENT fixed point inside every residual
# evaluation — linearize, return, re-read the abscissa off the RETURNED stress,
# re-linearize — which is the same fixed point the viscoplastic loop reaches over a
# whole solve, reached per evaluation instead.
_NR_ENV_MC, _NR_ENV_POW, _NR_ENV_HB = 0, 1, 2
_NR_ENV_MAX_ITER = 60      # self-consistency passes allowed per evaluation
_NR_ENV_RELAX1 = 14        # ... after this many, relax by 1/2
_NR_ENV_RELAX2 = 34        # ... and after this many, by 1/4
_NR_ENV_TOL = 1e-10        # ... and the relative change that ends them.
# 1e-10 and not tighter because the Balmer inversion is a 40-step BISECTION: its
# bracket closes at 2^-40 of the initial width, which puts a floor of about 6e-12
# under the tangent it can return, and a fixed point asked for less than that
# never exits. Measured on the Hoek-Brown benchmark, the residual falls by a
# decade a pass — 1.4e-1, 4.8e-3, 4.7e-4, ... — and then sits on that floor. The
# level itself is eight orders below `force_tol` and is not a physical setting.
_NR_HB_BISECT = 40         # Balmer inversion steps — the viscoplastic path's own


def _nr_pow_tangent(s_prime, a, b, cp, d, F, pow_ref=None):
    """The power curve's F-reduced tangent line at the in-plane Mohr-circle centre.

    ``tau = a (sigma_n + d)^b + c_p``, so ``d`` is the tension intercept: it shifts
    the origin, and the ``1e-4 ref`` floor is what stops a tensile point running
    past it. Every line here mirrors solve_fem's Step-6 block, including the
    abscissa: the CENTRE ``s' = -(sx + sy)/2`` clamped at zero, not the
    failure-plane normal, because that is the abscissa the viscoplastic driver
    linearizes at and the two have to be solving the same model.

    ``pow_ref`` is the group-wide scale of the floor. It is computed here when the
    caller has the whole group and passed back in for the tangent's difference
    quotient, so the perturbed evaluation floors against the same number the base
    one did rather than against a mean over the yielded subset.
    """
    s_n = np.maximum(s_prime, 0.0)
    if pow_ref is None:
        pow_ref = max(1.0, float(s_n.mean()) if s_n.size else 1.0)
    s_ef = np.maximum(s_n + d, 1e-4 * pow_ref)
    slope = a * b * s_ef ** (b - 1.0) / F
    tau_F = (a * s_ef ** b + cp) / F
    phi = np.arctan(slope)
    return tau_F - s_n * slope, np.sin(phi), np.cos(phi), pow_ref


def _nr_hb_tangent(s_prime, c_cur, snph_cur, csph_cur, sci, mb, s_, a_, F):
    """Hoek-Brown's F-reduced tangent at the FAILURE-PLANE normal stress.

    ``sigma_n = s' cos^2(phi) - c sin(phi) cos(phi)`` is where a Mohr circle of
    centre ``s'`` touches its own tangent line, which is why the linearization
    closes as a fixed point: at convergence the circle, the tangent line and the
    reduced envelope all meet at one ``sigma_n``. It is also the abscissa the
    limit-equilibrium engine uses (the slice-base normal), so the two solvers
    linearize the same curve at the same place. Balmer's parametric curve is then
    inverted by the same bisection the viscoplastic path runs.

    The ``sigma_n >= 0`` clamp inside ``hb_tangent_const`` is load-bearing and is
    documented there: at the Hoek-Brown tensile apex ``dsigma_1/dsigma_3`` diverges
    and ``phi_i -> 90 deg``, i.e. effectively infinite strength.

    Strength reduction divides the TANGENT — ``c_i/F`` and ``tan(phi_i)/F`` — and
    NOT the Hoek-Brown constants: ``sigma_ci/F`` is a different envelope, because of
    the exponent ``a``.
    """
    s_n = np.maximum(s_prime * csph_cur ** 2 - c_cur * snph_cur * csph_cur, 0.0)
    c_i, phi_i = hb_tangent_const(s_n, sci, mb, s_, a_, iters=_NR_HB_BISECT)
    slope = np.tan(np.radians(phi_i)) / F
    phi = np.arctan(slope)
    return c_i / F, np.sin(phi), np.cos(phi)


def _nr_envelope_step(code, c_cur, sn_cur, cs_cur, s_prime, env, pow_ref):
    """One linearization pass: every curved Gauss point gets its own new tangent.

    This is the per-element dispatch, done per GAUSS POINT: ``code`` carries
    _NR_ENV_MC / _NR_ENV_POW / _NR_ENV_HB from the model's own
    ``pow_flag_by_elem`` / ``hb_flag_by_elem``, so one group — one mesh — can hold
    Mohr-Coulomb, power-curve and Hoek-Brown points at once and each takes its own
    branch in the same pass. A Mohr-Coulomb point is not touched at all, which is
    what keeps a mixed model's MC half on exactly the arithmetic it would be on
    alone.
    """
    c, sn, cs = c_cur.copy(), sn_cur.copy(), cs_cur.copy()
    m = code == _NR_ENV_POW
    if m.any():
        cc, ss, cz, pow_ref = _nr_pow_tangent(
            s_prime[m], env['pow_a'][m], env['pow_b'][m], env['pow_cp'][m],
            env['pow_d'][m], env['F'][m], pow_ref)
        c[m], sn[m], cs[m] = cc, ss, cz
    m = code == _NR_ENV_HB
    if m.any():
        cc, ss, cz = _nr_hb_tangent(
            s_prime[m], c_cur[m], sn_cur[m], cs_cur[m], env['hb_sci'][m],
            env['hb_mb'][m], env['hb_s'][m], env['hb_a'][m], env['F'][m])
        c[m], sn[m], cs[m] = cc, ss, cz
    return c, sn, cs, pow_ref


def _nr_envelope_return(grp, sig_tr, sel=None, pow_ref=None):
    """The return map with any curved envelope linearized self-consistently.

    Returns ``(sig, branch, c, snph, csph, pow_ref)`` — the returned stress, the
    branch codes, and the tangent the return was actually taken on, which is what
    the verdict's yield reading has to be read against.

    The seed is the group's stored ``c_r`` / ``snph`` / ``csph``, which for a
    curved-envelope element is the overburden-estimate tangent build_fem_data
    lays down, already divided by F. It is a fixed, deterministic starting point,
    so the fixed point below is a function of the displacement and of nothing
    else — not of the iteration that came before it.

    ``sel`` restricts the whole thing to a subset of the group, which is what the
    tangent's difference quotient needs; ``pow_ref`` carries the base evaluation's
    floor scale into it.

    ``grp['_env_seed']`` replaces that cold seed once a solve is under way, and it
    is updated ONLY at an accepted iterate (in :func:`_nr_group_state`). Within one
    Newton iteration the seed is therefore frozen, so the residual the line search
    evaluates is a pure function of the displacement; and because the fixed point
    is driven to ``_NR_ENV_TOL`` the answer does not depend on where it started,
    only the number of passes to reach it does — 2 or 3 warm against 6 to 10 cold,
    which is most of this feature's cost.
    """
    env = grp['env']
    seed = grp.get('_env_seed')
    if sel is None:
        code = env['code']
        c, sn, cs = (grp['c_r'], grp['snph'], grp['csph']) if seed is None else seed
        mu, sub = grp['mu'], env
        t_cap, lam = grp.get('t_cap'), grp.get('lam')
    else:
        code = env['code'][sel]
        src = (grp['c_r'], grp['snph'], grp['csph']) if seed is None else seed
        c, sn, cs = src[0][sel], src[1][sel], src[2][sel]
        mu = grp['mu'][sel]
        sub = {k: v[sel] for k, v in env.items() if k != 'code'}
        t_cap = None if grp.get('t_cap') is None else grp['t_cap'][sel]
        lam = None if grp.get('lam') is None else grp['lam'][sel]
    c = np.array(c, dtype=float, copy=True)
    phi = np.arctan2(np.asarray(sn, float), np.asarray(cs, float))
    sn, cs = np.sin(phi), np.cos(phi)
    m_env = code != _NR_ENV_MC
    sig, branch = mc_return_map(sig_tr, c, sn, cs, mu, t_cap=t_cap, lam=lam)
    # Under-relaxation, applied only when the plain iteration stops making
    # progress. It is not a tuning dial: the fixed point it converges to is the
    # same one, and theta only decides how it gets there. It exists because a
    # HANDFUL of trial states put the iteration in a period-2 limit cycle — the
    # tangent moves the returned state across a branch boundary and the new
    # abscissa moves it back — and a cycle at one Gauss point stalls the whole
    # group. Measured on the fuzz below, damping is what takes the worst residual
    # from a stalled 2.2e-1 to convergence.
    theta = 1.0
    err = 0.0
    for _it in range(_NR_ENV_MAX_ITER):
        # The abscissa is read off the RETURNED stress, which is the stress the
        # viscoplastic loop reads it off too (its sig4 is the current state, not an
        # elastic predictor). At the fixed point the tangent is therefore taken
        # where the returned Mohr circle actually touches it.
        s_prime = -0.5 * (sig[:, 0] + sig[:, 1])
        c_new, sn_new, cs_new, pow_ref = _nr_envelope_step(
            code, c, sn, cs, s_prime, sub, pow_ref)
        fin = np.isfinite(c_new) & m_env
        cref = max(1.0, float(np.abs(c_new[fin] * cs_new[fin]).mean())
                   if fin.any() else 1.0)
        err = 0.0
        if fin.any():
            # The residual of the fixed point itself, g(x) - x, read BEFORE any
            # relaxation, so the convergence test means the same thing at every
            # theta.
            err = max(float(np.abs(c_new[fin] * cs_new[fin]
                                   - c[fin] * cs[fin]).max()) / cref,
                      float(np.abs(sn_new[fin] - sn[fin]).max()))
        if err < _NR_ENV_TOL:
            break
        # The schedule is on the PASS COUNT and not on the error trend, because a
        # trend test cannot tell a slow contraction from a cycle and this one does
        # not have to: an iteration that has not closed in _NR_ENV_RELAX1 passes is
        # either cycling, in which case halving breaks it, or converging at a rate
        # that halving costs a factor of two on. On the corpus the loop exits in one
        # to ten passes and never reaches either level.
        theta = (1.0 if _it < _NR_ENV_RELAX1
                 else 0.5 if _it < _NR_ENV_RELAX2 else 0.25)
        if theta == 1.0:
            c, sn, cs = c_new, sn_new, cs_new
        else:
            # Relax on (c, phi) and rebuild the sine and cosine from the angle, so
            # the pair can never drift off the unit circle. Mohr-Coulomb points are
            # excluded from the blend: their c is +inf where a material is held
            # elastic, and inf + theta*(inf - inf) is a NaN.
            phi_new = np.arctan2(sn_new, cs_new)
            c = np.where(m_env, c + theta * (c_new - c), c)
            phi = np.where(m_env, phi + theta * (phi_new - phi), phi)
            sn, cs = np.sin(phi), np.cos(phi)
        sig, branch = mc_return_map(sig_tr, c, sn, cs, mu, t_cap=t_cap, lam=lam)
    return sig, branch, c, sn, cs, pow_ref


def _nr_envelope_by_elem(fem_data, F_by_elem):
    """Per-element curved-envelope parameters, or ``None`` on a plain model.

    ``None`` is the whole of the no-envelope guarantee: a model with no ``pow`` or
    ``hb`` material never gets an ``env`` key on any group, so ``grp.get('env')``
    is None at the one point the linearization could enter and the return map is
    called with the group's own arrays exactly as it always was.
    """
    n = len(fem_data["elements"])
    pm = np.asarray(fem_data.get("pow_flag_by_elem", np.zeros(n, bool)), dtype=bool)
    hm = np.asarray(fem_data.get("hb_flag_by_elem", np.zeros(n, bool)), dtype=bool)
    if not (pm.any() or hm.any()):
        return None
    code = np.zeros(n, dtype=np.int8)
    code[pm] = _NR_ENV_POW
    code[hm] = _NR_ENV_HB
    z = np.zeros(n)
    out = {'code': code, 'F': np.asarray(F_by_elem, dtype=float).copy()}
    for key, src in (('pow_a', 'pow_a_by_elem'), ('pow_b', 'pow_b_by_elem'),
                     ('pow_cp', 'pow_cp_by_elem'), ('pow_d', 'pow_d_by_elem'),
                     ('hb_sci', 'hb_sci_by_elem'), ('hb_mb', 'hb_mb_by_elem'),
                     ('hb_s', 'hb_s_by_elem'), ('hb_a', 'hb_a_by_elem')):
        out[key] = np.asarray(fem_data.get(src, z), dtype=float)
    # pow_b defaults to 1 on a non-power element so the exponent is never 0^-1 in
    # a masked-off lane; the mask means it is never read, and this keeps it finite.
    out['pow_b'] = np.where(code == _NR_ENV_POW, out['pow_b'], 1.0)
    return out


def _nr_group_envelope(grp, env_by_elem):
    """Attach (or refresh) one group's per-Gauss-point envelope arrays.

    A Gauss point in a material held LINEAR ELASTIC is dispatched back to the
    plain Mohr-Coulomb code even where its material declares a curved envelope.
    On this path "elastic" is carried as an unreachable cohesion, ``c_r = inf``,
    and the linearization writes a finite tangent into ``c_r`` — so without this
    line a curved material inside an SSR elastic zone would start yielding, which
    is the opposite of what the zone says. The viscoplastic path never meets the
    problem because it drops elastic points from the yielding mask instead
    (``m & ~grp['elastic']``) and does not care what its ``c_r`` says. Measured on
    `vp040`, where 1,284 of 2,539 elements are held elastic: without it the first
    Newton step leaves 47% of the load uncarried at every load factor down to
    1/64, because half the mesh is yielding a material that cannot yield.
    """
    if env_by_elem is None:
        return
    e = grp['e_idx']
    env = {k: v[e] for k, v in env_by_elem.items()}
    if 'elastic' in grp:
        env['code'] = env['code'].copy()
        env['code'][grp['elastic']] = _NR_ENV_MC
    grp['env'] = env
    if not np.any(env['code'] != _NR_ENV_MC):
        # Nothing curved is left in this group, so it takes the plain path.
        del grp['env']


def _nr_group_restrength_envelope(grp, F_by_elem):
    """Re-reduce one group's envelope for a new trial strength.

    For a curved-envelope point the reduction is 1/F on the TANGENT, and the
    tangent is re-derived at every evaluation anyway, so all the ramp has to carry
    forward is the per-element F the linearization divides by. Without this a ramp
    would linearize the whole way up at the strength its foot was solved at.
    """
    env = grp.get('env')
    if env is not None:
        env['F'] = np.asarray(F_by_elem, dtype=float)[grp['e_idx']]


def _nr_group_state(grp, u, h_eps):
    """One group's stress, branch codes and algorithmic tangent at displacement u.

    Returns ``(sig4, branch, Dep)`` where ``Dep`` is the (G, 3, 3) tangent
    d[sx, sy, txy]/d[ex, ey, gxy] at fixed total eps_z = 0 — the block the global
    tangent stiffness needs. Elastic points carry the elastic D exactly; yielded
    points are differentiated through the return map by a one-sided difference in
    each of the three free strain components. The map is linear in the trial
    stress on each branch, so the quotient is exact there and the only error is
    the (smooth) rotation of the principal frame.
    """
    B, D4, dof = grp['B'], grp['D4'], grp['dof']
    ep = grp['ep']
    eps = (B @ u[dof][:, :, None])[:, :, 0]           # (G, 3) total strain
    eps4 = np.empty((len(eps), 4))
    eps4[:, :3] = eps - ep[:, :3]
    eps4[:, 3] = -ep[:, 3]
    sig_tr = (D4 @ eps4[:, :, None])[:, :, 0]
    # K0 initial stress, at the live load factor (see _nr_group_sig0). The return
    # map is handed the TOTAL trial stress, so sigma_0 sits inside the yield
    # surface evaluation exactly as it does on the viscoplastic path — and, being
    # constant in the strain, it drops out of the difference quotient below.
    _s0 = grp.get('sig0_s')
    if _s0 is not None:
        sig_tr = sig_tr + _s0
    t_cap, lam = grp.get('t_cap'), grp.get('lam')
    if grp.get('env') is None:
        sig, branch = mc_return_map(sig_tr, grp['c_r'], grp['snph'], grp['csph'],
                                    grp['mu'], t_cap=t_cap, lam=lam)
        pref = None
        grp['_c_eff'] = grp['c_r']
        grp['_snph_eff'] = grp['snph']
        grp['_csph_eff'] = grp['csph']
    else:
        sig, branch, _ce, _se, _cse, pref = _nr_envelope_return(grp, sig_tr)
        grp['_c_eff'], grp['_snph_eff'], grp['_csph_eff'] = _ce, _se, _cse
        # The accepted iterate's converged linearization becomes the next
        # evaluation's starting point. Only here, so the seed is frozen for the
        # whole of a line search.
        grp['_env_seed'] = (_ce, _se, _cse)
    yid = branch != _NR_ELASTIC
    # The min_slip_depth skin keeps the ELASTIC tangent, which is the viscoplastic
    # driver's own operator: that path solves K_elastic u = f_ext + f_body(eps^vp)
    # every iteration and never forms an elastoplastic tangent at all. A tangent is
    # a search direction and nothing else — the residual, the return map and the
    # committed plastic strain here are untouched — so this cannot move a converged
    # answer; what it removes is the near-null collapse mode a fully yielded
    # cohesionless skin puts into the consistent tangent, which the linear solve
    # otherwise answers with a correction 1e15 times the elastic response and no
    # backtrack can rescue. Measured on rs2_66b at F = 1.1: with the consistent
    # tangent the second iteration of every load increment reaches max|u| = 2e13
    # and the increment is dead; with this it stays at 0.08.
    _skin = grp.get('skin')
    if _skin is not None:
        yid = yid & ~_skin
    Dep = np.array(grp['D3'], copy=True)
    if np.any(yid):
        idx = np.flatnonzero(yid)
        base = sig[idx, :3]
        s_tr = sig_tr[idx]
        D4y = D4[idx]
        tc_y = None if t_cap is None else t_cap[idx]
        lam_y = None if lam is None else lam[idx]
        for j in range(3):
            if grp.get('env') is None:
                sp, _ = mc_return_map(s_tr + h_eps * D4y[:, :, j],
                                      grp['c_r'][idx], grp['snph'][idx],
                                      grp['csph'][idx], grp['mu'][idx],
                                      t_cap=tc_y, lam=lam_y)
            else:
                # Through the WHOLE map, linearization included, so the quotient
                # carries d(c, phi)/d(sigma) as well. The tangent is then the
                # derivative of the map the residual actually uses, with no
                # separate algebra to keep in step — at the price of the
                # linearization's own first-order truncation, alongside the
                # principal frame's, which is measured rather than assumed.
                sp = _nr_envelope_return(grp, s_tr + h_eps * D4y[:, :, j],
                                         sel=idx, pow_ref=pref)[0]
            Dep[idx, :, j] = (sp[:, :3] - base) / h_eps
    return sig, branch, Dep


def _nr_build_bars(fem_data):
    """The model's reinforcement bar elements, grouped by DOF count, or ``None``.

    One entry per distinct element size — four DOFs for a two-node bar, six for a
    three-node bar standing on the midside node of its soil edge — so each group is
    a dense array pass. Pile elements are excluded; the Newton driver refuses any
    model that has them.

    Every array here is read from ``fem_data`` and never written, so the two solvers
    carry the SAME bar: ``K_global_1d_elems`` is the element matrix the global
    elastic stiffness is assembled from, ``k_by_1d_elem`` the chord stiffness EA/L
    the force recovery uses, and ``t_allow_by_1d_elem`` the capacity envelope shared
    with the limit-equilibrium engine.

    ``p`` is the chord pattern [-cos, -sin, +cos, +sin, 0, 0]: ``p . u_e`` is the
    chord elongation, which on a three-node bar is the axial force station at the
    element center — the same quantity every reader of ``forces_1d`` expects.
    """
    elements_1d = fem_data.get("elements_1d", None)
    if elements_1d is None or len(elements_1d) == 0:
        return None
    n_1d = len(elements_1d)
    pile_mask = np.asarray(fem_data.get("pile_elem_mask", np.zeros(n_1d)), dtype=bool)
    reinf = np.flatnonzero(~pile_mask)
    if reinf.size == 0:
        return None
    K_list = fem_data.get("K_global_1d_elems", [])
    dof_all = np.asarray(fem_data["dof_indices_1d"], dtype=np.int64)
    nd_all = np.asarray(fem_data.get("n_dof_by_1d_elem",
                                     np.full(n_1d, 4, dtype=int)), dtype=int)
    k_all = np.asarray(fem_data["k_by_1d_elem"], dtype=float)
    cos_all = np.asarray(fem_data["cos_theta_1d"], dtype=float)
    sin_all = np.asarray(fem_data["sin_theta_1d"], dtype=float)
    cap_all = np.asarray(fem_data["t_allow_by_1d_elem"], dtype=float)
    # The post-peak residual, and which bars can reach it. `build_fem_data` has
    # already taken min(t_res, t_allow), so a bar inside a pullout ramp where
    # t_allow < t_res cannot soften at all: bond slip stays perfectly plastic and
    # only RUPTURE softens. Blank t_res is NaN, which is "no post-peak drop".
    res_all = np.asarray(fem_data.get("t_res_by_1d_elem",
                                      np.full(n_1d, np.nan)), dtype=float)
    can_all = np.isfinite(res_all) & (res_all < cap_all - 1e-12)

    bars = []
    for nd in (4, 6):
        sel = reinf[nd_all[reinf] == nd]
        if sel.size == 0:
            continue
        p = np.zeros((sel.size, nd))
        p[:, 0], p[:, 1] = -cos_all[sel], -sin_all[sel]
        p[:, 2], p[:, 3] = cos_all[sel], sin_all[sel]
        K = np.array([np.asarray(K_list[i], dtype=float)[:nd, :nd] for i in sel])
        bars.append({
            'idx': sel, 'nd': nd, 'dof': dof_all[sel, :nd], 'K': K, 'p': p,
            'k': k_all[sel],
            # `t_cap` is the WORKING cap, which is the peak until a bar softens.
            't_cap': cap_all[sel].copy(),
            't_peak': cap_all[sel], 't_res': res_all[sel],
            'can_soften': can_all[sel],
            'softened': np.zeros(sel.size, dtype=bool),
            # k p (x) p, the rank-one block the tangent LOSES when a bar goes
            # slack or reaches its capacity. On a two-node bar it is K itself.
            'kpp': k_all[sel][:, None, None] * (p[:, :, None] * p[:, None, :]),
        })
    return bars or None


def _nr_translational_dofs(fem_data, n_dof):
    """Indices of the degrees of freedom that are LENGTHS, or None when all are.

    A pile node carries a third degree of freedom, its rotation, and a rotation is
    not comparable with a displacement: put one in a `max|u|` and the bound is
    comparing a length against a radian. The viscoplastic driver's own
    displacement-limit check already reads translational degrees of freedom only
    ("Extract translational DOFs only for VP displacement check"); this is the same
    reading for the Newton path.

    ``None`` on a model with no pile — every entry is already a length there, so the
    caller takes the plain ``np.max(np.abs(u))`` it always took and the arithmetic is
    bit-identical.
    """
    dof_offset = fem_data.get("dof_offset", None)
    is_pile = fem_data.get("is_pile_node", None)
    if dof_offset is None or is_pile is None:
        return None
    is_pile = np.asarray(is_pile, dtype=bool)
    if not is_pile.any():
        return None
    rot = np.asarray(dof_offset[:len(is_pile)], dtype=np.int64)[is_pile] + 2
    keep = np.ones(n_dof, dtype=bool)
    rot = rot[(rot >= 0) & (rot < n_dof)]
    keep[rot] = False
    return np.flatnonzero(keep)


def _nr_umax(u, trans_dofs):
    """max|u| over the translational degrees of freedom (all of them without a pile)."""
    v = u if trans_dofs is None else u[trans_dofs]
    return float(np.max(np.abs(v))) if v.size else 0.0


def _nr_umax_split(u, trans_dofs, deep_dof_mask):
    """``(deep, skin)`` maxima of |u|, split by the ``min_slip_depth`` filter.

    With no filter (``deep_dof_mask is None``) the deep value is the ordinary
    :func:`_nr_umax` and the skin value is ``None`` — the filtered path adds a
    reading, it does not change the unfiltered one.
    """
    if deep_dof_mask is None:
        return _nr_umax(u, trans_dofs), None
    if trans_dofs is None:
        trans = np.ones(len(u), dtype=bool)
    else:
        trans = np.zeros(len(u), dtype=bool)
        trans[trans_dofs] = True
    deep = trans & deep_dof_mask
    skin = trans & ~deep_dof_mask
    d = float(np.max(np.abs(u[deep]))) if deep.any() else 0.0
    k = float(np.max(np.abs(u[skin]))) if skin.any() else 0.0
    return d, k


def _nr_oob_split(r_full, free_dof_mask, node_dof_x, node_dof_y, node_has_free,
                  g_node_den, deep_free_mask):
    """``(whole model, skin)`` Dawson out-of-balance from one residual vector.

    The skin reading is ``None`` without a ``min_slip_depth`` filter, where there is
    no skin to separate and the whole-model value is the reported one.
    """
    d = r_full * free_dof_mask
    rn = np.sqrt(d[node_dof_x] ** 2 + d[node_dof_y] ** 2)
    v = (rn / g_node_den)[node_has_free]
    if v.size == 0:
        return 0.0, (None if deep_free_mask is None else 0.0)
    whole = float(np.max(v))
    if deep_free_mask is None:
        return whole, None
    sk = v[~deep_free_mask]
    return whole, (float(np.max(sk)) if sk.size else 0.0)


def _nr_build_piles(fem_data):
    """The model's pile beam elements, grouped by DOF count, or ``None``.

    One entry per distinct element size — six degrees of freedom for a two-node beam,
    nine for a three-node one standing on the midside node of its soil edge as well —
    so each group is a dense array pass.

    Every array is read from ``fem_data`` and never written, so the two solvers carry
    the SAME beam: ``K_global_pile_elems`` is the element matrix the global elastic
    stiffness is assembled from, and ``V_cap_by_pile_elem`` / ``M_cap_by_pile_elem``
    are the per-unit-width capacities (the file's single-pile values divided by the
    spacing) the viscoplastic path checks against.

    ``act`` holds the three rows that read the ACTIONS out of the element's global
    displacement — the shear at the element center and the bending moment at each
    end node — so ``act @ u_e`` is ``(V, M1, M2)``. ``patV`` holds the internal
    force pattern the SHEAR capacity is delivered on: the transverse pair at the two
    end nodes, rotated into global coordinates. ``ax`` reads the axial force at the
    element center, which carries no capacity and is reported only.

    The MOMENT capacity is not delivered on a pattern at all — it is a plastic hinge
    (see :func:`_nr_pile_force`), so what the group carries for it is ``S``, the
    rotational block of the element stiffness that the hinge is solved on, and
    ``hinge``, the mask of which of the element's two ends may take the release.
    That mask is built ONCE over the whole pile element list, before the grouping,
    because a hinge is one release at one pile NODE and the two element ends meeting
    at an interior node can land in different DOF-count groups.
    """
    n_pile = int(fem_data.get("n_pile_elements", 0))
    if n_pile == 0:
        return None
    K_list = fem_data["K_global_pile_elems"]
    K_local_list = fem_data.get("K_local_by_pile_elem", None)
    dof_all = np.asarray(fem_data["dof_indices_pile"], dtype=np.int64)
    nd_all = np.asarray(fem_data.get("n_dof_by_pile_elem",
                                     np.full(n_pile, 6, dtype=int)), dtype=int)
    cos_all = np.asarray(fem_data["cos_theta_pile"], dtype=float)
    sin_all = np.asarray(fem_data["sin_theta_pile"], dtype=float)
    L_all = np.asarray(fem_data["elem_length_by_pile_elem"], dtype=float)
    EI_all = np.asarray(fem_data["EI_by_pile_elem"], dtype=float)
    EA_all = np.asarray(fem_data["EA_by_pile_elem"], dtype=float)
    V_all = np.asarray(fem_data["V_cap_by_pile_elem"], dtype=float)
    M_all = np.asarray(fem_data["M_cap_by_pile_elem"], dtype=float)
    # Which element END may carry a plastic hinge, one per pile node. The same
    # function the viscoplastic path uses, called on the WHOLE element list before
    # the DOF-count grouping so the ownership is global: the two element ends that
    # meet at an interior pile node can sit in different groups.
    hinge_all = _pile_hinge_ends(fem_data.get("pile_elem_nodes", None), n_pile)

    piles = []
    for nd in (6, 9):
        sel = np.flatnonzero(nd_all[:n_pile] == nd)
        if sel.size == 0:
            continue
        n_node = nd // 3
        act = np.zeros((sel.size, 3, nd))
        patV = np.zeros((sel.size, nd))
        ax = np.zeros((sel.size, nd))
        K = np.zeros((sel.size, nd, nd))
        Kl_g = np.zeros((sel.size, nd, nd))
        S = np.zeros((sel.size, 2, 2))
        for j, i in enumerate(sel):
            L, EI, EA = L_all[i], EI_all[i], EA_all[i]
            T = _beam_rotation(cos_all[i], sin_all[i], n_node)
            Kl = (np.asarray(K_local_list[i], dtype=float)[:nd, :nd]
                  if K_local_list is not None
                  else _beam_local_stiffness(EA, EI, L, n_node == 3))
            K[j] = np.asarray(K_list[i], dtype=float)[:nd, :nd]
            Kl_g[j] = Kl
            # The rotational block the hinge is solved on, rows and columns 2 and 5
            # — the two END nodes. It is the same block in either frame: the beam
            # rotation is the identity on a rotational row, so
            # (T^T Kl T)[2, 5] == Kl[2, 5].
            S[j] = Kl[np.ix_([2, 5], [2, 5])]
            # --- the three actions, as local rows ---
            rV = np.zeros(nd)
            if n_node == 2:
                rV[1], rV[2] = 12.0 * EI / L ** 3, 6.0 * EI / L ** 2
                rV[4], rV[5] = -12.0 * EI / L ** 3, 6.0 * EI / L ** 2
            else:
                # V = EI/L**3 * (d3H/dxi3 at the center) . q_hat, with q_hat the
                # bending dofs scaled to the unit interval: [v, L*theta] per node.
                scale = np.array([1.0, L, 1.0, L, 1.0, L])
                coef = (EI / L ** 3) * (_QUINTIC_SHEAR_MID * scale)
                for k, slot in enumerate((1, 2, 4, 5, 7, 8)):
                    rV[slot] = coef[k]
            act[j, 0] = rV @ T
            act[j, 1] = Kl[2] @ T          # M1, the end-1 moment
            act[j, 2] = Kl[5] @ T          # M2, the end-2 moment
            rT = np.zeros(nd)
            rT[0], rT[3] = -EA / L, EA / L
            ax[j] = rT @ T
            # --- the pattern the SHEAR capacity is delivered on ---
            pV = np.zeros(nd)
            pV[1], pV[4] = 1.0, -1.0
            patV[j] = T.T @ pV
        piles.append({'idx': sel, 'nd': nd, 'dof': dof_all[sel, :nd], 'K': K,
                      'Kl': Kl_g, 'S': S, 'act': act, 'patV': patV, 'ax': ax,
                      'Vcap': V_all[sel], 'Mcap': M_all[sel],
                      'hinge': hinge_all[sel]})
    return piles or None


def _nr_pile_force(pg, u, want_tangent):
    """One pile group's internal force (and tangent) at displacement ``u``.

    The beam is linear elastic; the only nonlinearities are its two structural
    capacities, and they are two different KINDS of thing.

    **The moment capacity is a plastic hinge**, the same one the viscoplastic path
    carries (:func:`_pile_moment_hinge`, :func:`_pile_hinge_ends`). The element's
    elastic end rotation is the nodal rotation less a plastic rotation ``p``, so the
    element delivers ``K (u_e - p)``. A correction applied instead to the rotational
    degree of freedom the two adjacent beam elements SHARE is equal and opposite
    there and cancels exactly, which is a capacity that enforces nothing; the
    element vector ``K p`` has a translational part and does not. Exactly one of the
    two element ends at a pile node may take the release, and the other is left
    elastic — node equilibrium puts it on the capacity with the opposite sign
    anyway, and releasing both would leave the node's rotation undetermined by the
    beam.

    Written in the element's GLOBAL frame the correction is ``K_global p`` outright,
    because ``p`` carries only ROTATIONAL components and a rotation is invariant
    under the beam's frame change (``T^T K_local p_local = K_global p``). So the
    residual is the element matrix applied to the RELEASED displacement.

    **The shear capacity is the bar's convention** — the part of the elastic action
    the member cannot deliver, subtracted from its own internal force — read on the
    released displacement, because a hinge changes the shear:

        f_int = K (u_e - p) - (V - clip(V, -V_cap, +V_cap)) q_V

    A member at its cap therefore delivers the cap, which is what a cap is.

    **The tangent is exact, with the rotational block condensed.** On a fixed hinge
    set ``A``, ``p_A = S_AA^-1 (m_A - sign(m_A) M_cap)`` is affine in ``u_e`` — every
    end moment is linear in ``u_e`` and ``S`` is the rotational block of the element
    stiffness — so ``dp/du_e`` is the constant matrix ``R`` carrying ``S_AA^-1 G_A``
    in the released rotational rows and zero elsewhere, with ``G`` the two moment
    rows. Then ``dw/du_e = I - R`` and

        dK/du_e = K (I - R) - [|V| > V_cap] q_V (x) (g_V (I - R))

    which is the derivative of the map the residual actually uses. A hinged end's
    delivered moment has derivative exactly zero, which is what a hinge means.

    There is no state. ``p`` is a function of the current displacement alone and the
    ``yielded_pile_*`` masks are read on the reported state, so the pile law is
    nonlinear-ELASTIC and nothing is committed at the end of a step.
    """
    ue = u[pg['dof']]
    n, nd = ue.shape
    gM = pg['act'][:, 1:3]                              # the two end-moment rows
    m_e = np.einsum('ikj,ij->ik', gM, ue)               # elastic end moments
    Mcap = pg['Mcap']
    p = np.zeros((n, nd))
    R = None
    # yielded_M is set wherever the capacity is EXCEEDED, whether or not this
    # element's own end took the release: the other end of a shared node is on the
    # capacity too. Same convention as _pile_element_capacity.
    yM = np.isfinite(Mcap) & (np.abs(m_e) > Mcap[:, None]).any(axis=1)
    if yM.any():
        R = np.zeros((n, nd, nd))
        slots = np.array([2, 5])
        for i in np.flatnonzero(yM):
            p_rot, _m = _pile_moment_hinge(m_e[i, 0], m_e[i, 1], Mcap[i],
                                           pg['Kl'][i], allowed=pg['hinge'][i])
            act_i = np.flatnonzero(p_rot != 0.0)
            if act_i.size == 0:                 # over the cap at a node it does
                continue                        # not own; the other end releases
            p[i, 2], p[i, 5] = p_rot[0], p_rot[1]
            R[i, slots[act_i], :] = np.linalg.solve(
                pg['S'][i][np.ix_(act_i, act_i)], gM[i][act_i])
        if not R.any():
            R = None
    w = ue - p
    V = np.einsum('ij,ij->i', pg['act'][:, 0], w)
    Vcap = pg['Vcap']
    V_true = np.clip(V, -Vcap, Vcap)
    f = np.einsum('ijk,ik->ij', pg['K'], w) - (V - V_true)[:, None] * pg['patV']
    pg['_V'], pg['_V_true'], pg['_m_e'] = V, V_true, m_e
    pg['_M'] = np.einsum('ikj,ij->ik', gM, w)           # the delivered end moments
    pg['_p_rot'] = p[:, [2, 5]]
    pg['_yV'], pg['_yM'] = np.abs(V) > Vcap, yM
    pg['_axial'] = np.einsum('ij,ij->i', pg['ax'], w)
    if want_tangent:
        aV = (np.abs(V) > Vcap)[:, None, None]
        gV = pg['act'][:, 0]
        if R is None:
            pg['_Ke'] = pg['K'] - aV * (pg['patV'][:, :, None] * gV[:, None, :])
        else:
            A = np.eye(nd)[None, :, :] - R
            gVA = np.einsum('ij,ijk->ik', gV, A)
            pg['_Ke'] = (np.einsum('ijk,ikl->ijl', pg['K'], A)
                         - aV * (pg['patV'][:, :, None] * gVA[:, None, :]))
    return f


def _nr_bar_force(bg, u, want_tangent):
    """One bar group's internal force (and tangent) at displacement ``u``.

    The physics is the viscoplastic driver's, restated as a residual. There the
    global stiffness carries the bar's FULL elastic matrix and the part of the
    elastic force the bar cannot deliver is subtracted as a body load, so the bar's
    internal force is

        K_bar u_e - (T - T_true) p,   T = k (p . u_e),   T_true = clip(T, 0, t_cap)

    — tension only, capped at the pullout/anchorage envelope, and on a two-node bar
    that expression is exactly ``T_true p``. On a three-node bar it is not, and this
    reproduces the three-node expression rather than a cleaner one, because the two
    drivers have to be solving the same model.

    The map is piecewise LINEAR in u_e, so its tangent is exact with nothing
    differenced: dT/du = k p, and dT_true/du is k p while the bar is active
    (0 < T < t_cap) and zero otherwise. The element tangent is therefore K_bar for
    an active bar and K_bar - k p (x) p for a slack or capacity-limited one. On a
    two-node bar the second is exactly zero — a yielded bar adds no stiffness. On a
    three-node bar it is a rank-one positive-semidefinite block that resists only
    the midside node leaving the chord.

    There is no state: with softening refused, T_true is a function of the current
    displacement alone, so the bar law is nonlinear-ELASTIC and nothing is committed
    at the end of a step.
    """
    ue = u[bg['dof']]
    T = bg['k'] * np.einsum('ij,ij->i', bg['p'], ue)
    T_true = np.clip(T, 0.0, bg['t_cap'])
    f = np.einsum('ijk,ik->ij', bg['K'], ue) - (T - T_true)[:, None] * bg['p']
    bg['_T'], bg['_T_true'] = T, T_true
    if want_tangent:
        active = (T > 0.0) & (T < bg['t_cap'])
        bg['_Ke'] = np.where(active[:, None, None], bg['K'], bg['K'] - bg['kpp'])
    return f


def _nr_internal_force(gp_groups, u, n_dof, h_eps=None, want_tangent=False,
                       bars=None, piles=None):
    """Assemble the internal force vector (and optionally the per-group tangents).

    ``h_eps=None`` with ``want_tangent=False`` is the cheap residual-only pass the
    line search uses: one return map per group and no differentiation.

    ``bars`` adds the reinforcement bar elements' contribution and ``piles`` the
    pile beam elements'. Their tangents are stashed on their own groups (``_Ke``)
    rather than returned, because they need no differencing and the caller consumes
    them through the same assembly pattern.
    """
    fint = np.zeros(n_dof)
    tangents = [] if want_tangent else None
    n_plastic_gp = 0
    for grp in gp_groups:
        if want_tangent:
            sig, branch, Dep = _nr_group_state(grp, u, h_eps)
            tangents.append(Dep)
        else:
            B, D4, dof = grp['B'], grp['D4'], grp['dof']
            ep = grp['ep']
            eps = (B @ u[dof][:, :, None])[:, :, 0]
            eps4 = np.empty((len(eps), 4))
            eps4[:, :3] = eps - ep[:, :3]
            eps4[:, 3] = -ep[:, 3]
            sig_tr = (D4 @ eps4[:, :, None])[:, :, 0]
            _s0 = grp.get('sig0_s')
            if _s0 is not None:
                sig_tr = sig_tr + _s0
            # The SAME constitutive law the tangent path takes, tensile cap
            # and K0 initial stress included. This branch exists only to skip the differencing, so any
            # divergence between the two is a residual and a tangent for two
            # different materials.
            if grp.get('env') is None:
                sig, branch = mc_return_map(sig_tr, grp['c_r'], grp['snph'],
                                            grp['csph'], grp['mu'],
                                            t_cap=grp.get('t_cap'),
                                            lam=grp.get('lam'))
                grp['_c_eff'] = grp['c_r']
                grp['_snph_eff'] = grp['snph']
                grp['_csph_eff'] = grp['csph']
            else:
                sig, branch, _ce, _se, _cse, _ = _nr_envelope_return(grp, sig_tr)
                grp['_c_eff'], grp['_snph_eff'], grp['_csph_eff'] = _ce, _se, _cse
        grp['_sig'] = sig
        grp['_branch'] = branch
        n_plastic_gp += int(np.count_nonzero(branch))
        contrib = (sig[:, None, :3] @ grp['B'])[:, 0, :] * grp['w'][:, None]
        np.add.at(fint, grp['dof'].ravel(), contrib.ravel())
    if bars is not None:
        for bg in bars:
            fb = _nr_bar_force(bg, u, want_tangent)
            np.add.at(fint, bg['dof'].ravel(), fb.ravel())
    if piles is not None:
        for pg in piles:
            fp = _nr_pile_force(pg, u, want_tangent)
            np.add.at(fint, pg['dof'].ravel(), fp.ravel())
    return fint, tangents, n_plastic_gp


def _nr_commit_plastic_strain(gp_groups):
    """Freeze the converged step: eps^p := eps - D4^-1 sigma, per Gauss point.

    The elastic strain is recovered from the returned stress rather than from the
    increment, so the committed plastic strain is exactly consistent with the
    stress the residual was built from — no drift accumulates across steps.
    """
    for grp in gp_groups:
        eps = (grp['B'] @ grp['_u'][grp['dof']][:, :, None])[:, :, 0]
        # The ELASTIC part of the stress is what D4inv inverts: with a K0 initial
        # stress the constitutive relation is sigma = sigma_0 + D (eps - eps^p), so
        # sigma_0 comes off before the inversion or the whole in-situ field would be
        # booked as elastic strain.
        _sig = grp['_sig']
        _s0 = grp.get('sig0_s')
        if _s0 is not None:
            _sig = _sig - _s0
        eps_e4 = (grp['D4inv'] @ _sig[:, :, None])[:, :, 0]
        ep = np.empty_like(grp['ep'])
        ep[:, :3] = eps - eps_e4[:, :3]
        ep[:, 3] = -eps_e4[:, 3]
        grp['ep'] = ep


def _nr_group_tension_cap(grp, t_cap_by_elem):
    """Attach (or refresh) one group's Rankine tensile cap.

    The keys are set ONLY where the group actually carries a finite cap, so
    ``grp.get('t_cap')`` is None on every model that sets no ``t_cut`` and the
    return map's cap code is never entered there. A material declared ``elastic``
    is held out of the cap exactly as the viscoplastic path holds it out
    (``tm & ~grp['elastic']``): it cannot fail, so it cannot crack either.
    """
    if t_cap_by_elem is None:
        return
    tc = np.asarray(t_cap_by_elem, dtype=float)[grp['e_idx']]
    if 'elastic' in grp:
        tc = tc.copy()
        tc[grp['elastic']] = np.inf
    if not np.isfinite(tc).any():
        grp.pop('t_cap', None)
        grp.pop('lam', None)
        return
    grp['t_cap'] = tc
    # The associated Rankine flow returns 2*mu*e_k + lambda*1, so it needs the
    # second Lame constant. D4[0, 1] is exactly lambda (see
    # build_constitutive_matrix_4), which is why no second material lookup happens.
    grp['lam'] = grp['D4'][:, 0, 1].copy()


def _nr_group_suction(grp, sg, prep):
    """Attach one group's capped matric suction, the F-independent half.

    ``s_suc`` is ``min(max(0, -u_w), s_cap)`` per Gauss point, read off the same
    SIGNED pore-pressure field (``prep['u_gp_signed']``) and the same per-element
    cap the viscoplastic path reads, through the same helper. The keys are set
    only on a model that actually carries suction, so ``grp.get('s_suc')`` is None
    everywhere else and none of the apparent-cohesion arithmetic is entered.

    The credit itself is F-dependent and is folded into the cohesion by
    :func:`_nr_group_apparent_cohesion`, which the ramp re-runs at every step.
    """
    if not prep.get("suction_active"):
        return
    u_sgn_all = prep.get("u_gp_signed")
    if u_sgn_all is None:
        # Suction is switched on but the pore-pressure source carries no signed
        # field (u = none or ru). The viscoplastic path credits zero there; so
        # does this, rather than skipping the keys, so the two agree on a model
        # that asks for suction it cannot have.
        grp['tanphib'] = np.zeros(grp['n'])
        grp['s_suc'] = np.zeros(grp['n'])
        return
    u_sgn = np.array([u_sgn_all[e][g] for e, g in sg['pairs']], dtype=float)
    grp['tanphib'] = sg['tanphib']
    grp['s_suc'] = _suction_capped(u_sgn, sg['scap'])


def _nr_group_apparent_cohesion(grp, finv_by_elem):
    """Fold the F-reduced apparent cohesion into one group's cohesion.

    This is the whole of matric suction on the Newton path. The return map
    already takes ``c`` per Gauss point, so Fredlund's credit is a per-point
    cohesion and nothing else: no surface, no flow direction and no tangent
    changes, and everything downstream that reads ``c_r`` -- the return map, the
    algorithmic tangent, the strength scale ``nr_max_yield_violation`` divides by
    -- reads the envelope the trial is actually solved on.

    A no-op on a model without suction. On an ``elastic`` material ``c_r`` is
    already ``inf``, so the credit is inert there exactly as it is on the
    viscoplastic path, which drops elastic points from the yielding mask.
    """
    s_suc = grp.get('s_suc')
    if s_suc is None:
        return
    c_suc = _suction_apparent_cohesion(grp['tanphib'], s_suc,
                                       finv_by_elem[grp['e_idx']])
    grp['c_suc'] = c_suc
    grp['c_r'] = grp['c_r'] + c_suc


def _nr_group_sig0(grp, sg, prep, k0):
    """Attach one group's K0 initial stress, or leave the group without one.

    The field is the viscoplastic driver's, read from the same cached overburden
    (``prep['sv0_gp']``, F-independent) and built by the same formula, per Gauss
    point and tension-positive:

        sigma'_v = -(soil overburden) + u        (u >= 0)
        sigma'_h = sigma'_z = K0 sigma'_v        (in-plane AND out-of-plane)
        tau_xy   = 0

    The field is EFFECTIVE, as it is on the viscoplastic path: the pore-pressure
    term is already in the load vector and D (B u - eps^p) is the effective stress
    directly, so the initial stress is effective too and needs no correction.

    ``sig0`` is the full in-situ field and ``sig0_s`` is that field at the live
    LOAD FACTOR: the driver walks gravity in increments and the overburden IS
    gravity, so the two scale together and lam = 0 is a genuinely stress-free
    state. The key is set only on a K0 run, so ``grp.get('sig0_s')`` is None on
    every other model and none of the initial-stress arithmetic is entered there.
    """
    if k0 is None or prep.get("sv0_gp") is None:
        return
    sv0_gp, u_gp_all = prep["sv0_gp"], prep["u_gp"]
    pairs = sg['pairs']
    sv0 = np.array([sv0_gp[e][g] for e, g in pairs], dtype=float)
    u_st = np.array([u_gp_all[e][g] for e, g in pairs], dtype=float)
    sv_eff = -sv0 + u_st
    sh_eff = float(k0) * sv_eff
    z = np.zeros_like(sv_eff)
    grp['sig0'] = np.stack([sh_eff, sv_eff, z, sh_eff], axis=1)
    grp['sig0_s'] = grp['sig0']


def _nr_set_load_factor(groups, lam):
    """Scale every group's initial stress to the load factor being attempted.

    A no-op on a model without K0. Within one increment ``lam`` is fixed, so the
    residual, the line search and the tangent all see the same material.
    """
    for grp in groups:
        s0 = grp.get('sig0')
        if s0 is not None:
            grp['sig0_s'] = s0 if lam == 1.0 else lam * s0


def _nr_build_groups(prep, c_reduced, phi_reduced, elastic_by_elem,
                     t_cap_by_elem=None, k0=None, finv_by_elem=None,
                     env_by_elem=None):
    """The Newton path's working Gauss-point groups.

    Geometry (B, D4, w, dof) is shared BY REFERENCE with the prepared model, which
    the viscoplastic path also reads and never writes. Everything the Newton solve
    mutates — the committed plastic strain, the stress, the branch codes — is
    fresh per solve, so a shared prepared model cannot carry state between trials
    or between the two solvers.
    """
    groups = []
    for sg in prep["gp_groups_static"]:
        e_idx = sg['e_idx']
        G = sg['n']
        D4 = sg['D4']
        grp = {
            'pairs': sg['pairs'], 'e_idx': e_idx, 'n': G,
            'B': sg['B'], 'D4': D4, 'w': sg['w'], 'dof': sg['dof'],
            'D3': D4[:, :3, :3],
            'D4inv': np.linalg.inv(D4),
            # mu is the (3,3) entry of the plane-strain operator, E/(2(1+nu)),
            # so no second lookup of the material table is needed.
            'mu': D4[:, 2, 2].copy(),
            'c_r': c_reduced[e_idx],
            'snph': np.sin(phi_reduced[e_idx]),
            'csph': np.cos(phi_reduced[e_idx]),
            'ep': np.zeros((G, 4)),
        }
        # The surficial skin the min_slip_depth filter excludes (None when the
        # filter is off). It changes ONE thing: the tangent at these points stays
        # elastic — see _nr_group_state. Nothing about the internal force, the
        # return map or the committed plastic strain differs, so the equilibrium
        # being solved is the same one.
        _skin = prep.get("skin_elem_mask")
        if _skin is not None:
            _sk = _skin[e_idx]
            if _sk.any():
                grp['skin'] = _sk
        if elastic_by_elem is not None:
            em = elastic_by_elem[e_idx]
            if em.any():
                # A material held linear elastic never yields: give it an
                # unreachable cohesion rather than branching in the return map.
                grp['c_r'] = grp['c_r'].copy()
                grp['c_r'][em] = np.inf
                grp['elastic'] = em
        # Matric suction, after the elastic override: an elastic point's c_r is
        # already inf and adding a finite credit to it leaves it inf, which is the
        # same inertness the viscoplastic path gets by dropping the point from the
        # yielding mask.
        _nr_group_suction(grp, sg, prep)
        if finv_by_elem is not None:
            _nr_group_apparent_cohesion(grp, finv_by_elem)
        _nr_group_tension_cap(grp, t_cap_by_elem)
        _nr_group_sig0(grp, sg, prep, k0)
        # The curved envelopes last: they carry no cohesion of their own, only the
        # parameters the tangent is re-derived from. c_r / snph / csph stay the
        # SEED tangent build_fem_data laid down, which is what the fixed point
        # starts from.
        _nr_group_envelope(grp, env_by_elem)
        groups.append(grp)
    return groups


def _nr_prepare_assembly(groups, free_dofs, n_dof, bars=None, piles=None):
    """Cache the assembly pattern every tangent re-form reuses.

    The sparsity pattern is fixed for the whole solve — only the VALUES change
    from one Newton iteration to the next — so the (row, column) bookkeeping is
    done once here and each assembly becomes a single ``np.bincount`` into a
    ready-made CSC structure. Rebuilding a COO matrix and re-sorting it into CSC
    every iteration, which is the obvious way to write this, costs as much as the
    factorization it feeds.

    Reinforcement bar blocks join the pattern EXPLICITLY, after the soil groups and
    in the order :func:`_nr_assemble_tangent` walks them. On these meshes a bar lies
    along a soil element edge, so its DOF pairs happen to be in the soil pattern
    already — but relying on that would make the tangent silently lose a bar whose
    end nodes did not share an element.

    Pile beam blocks join it last, and for them the explicit join is not a precaution
    but a requirement: a pile node's ROTATIONAL degree of freedom appears in no soil
    element and in no bar, so every rotation would otherwise be missing from the
    tangent entirely.
    """
    n_free = len(free_dofs)
    fidx = np.full(n_dof, -1, dtype=np.int64)
    fidx[free_dofs] = np.arange(n_free)
    lin_parts = []
    for grp in (list(groups) + (list(bars) if bars else [])
                + (list(piles) if piles else [])):
        fd = fidx[grp['dof']]
        G, nd = fd.shape
        rows = np.broadcast_to(fd[:, :, None], (G, nd, nd))
        cols = np.broadcast_to(fd[:, None, :], (G, nd, nd))
        keep = (rows >= 0) & (cols >= 0)
        grp['_keep'] = keep
        # Column-major key, so the unique keys come out in CSC order directly.
        lin_parts.append((np.asarray(cols)[keep] * np.int64(n_free)
                          + np.asarray(rows)[keep]))
    lin = np.concatenate(lin_parts) if lin_parts else np.zeros(0, dtype=np.int64)
    uniq, inv = np.unique(lin, return_inverse=True)
    indices = (uniq % n_free).astype(np.int32)
    col_of = (uniq // n_free).astype(np.int64)
    indptr = np.searchsorted(col_of, np.arange(n_free + 1)).astype(np.int32)
    return {"inv": inv, "nnz": len(uniq), "indices": indices,
            "indptr": indptr, "n_free": n_free}


def _nr_assemble_tangent(groups, tangents, pattern, bars=None, piles=None):
    """Global non-symmetric tangent stiffness on the free degrees of freedom.

    Soil groups first, then the reinforcement bar groups, then the pile beam groups
    — the same order :func:`_nr_prepare_assembly` built the pattern in.
    """
    vals = []
    for grp, Dep in zip(groups, tangents):
        Bt = np.ascontiguousarray(grp['B'].transpose(0, 2, 1))
        ke = (Bt @ Dep) @ grp['B']
        ke *= grp['w'][:, None, None]
        vals.append(ke[grp['_keep']])
    if bars is not None:
        for bg in bars:
            vals.append(bg['_Ke'][bg['_keep']])
    if piles is not None:
        for pg in piles:
            vals.append(pg['_Ke'][pg['_keep']])
    data = np.bincount(pattern["inv"], weights=np.concatenate(vals),
                       minlength=pattern["nnz"])
    n_free = pattern["n_free"]
    return csc_matrix((data, pattern["indices"], pattern["indptr"]),
                      shape=(n_free, n_free))


class _NrPermutedLU:
    """``.solve(b)`` over a factorization taken on the column-permuted tangent."""

    __slots__ = ("_lu", "_pc")

    def __init__(self, lu, pc):
        self._lu = lu
        self._pc = pc

    def solve(self, b):
        return self._lu.solve(b)[self._pc]


def _nr_factorize_tangent(K, cache):
    """``splu`` on a Newton tangent, re-using the column ordering the pattern fixes.

    ``splu``'s default ``permc_spec='COLAMD'`` orders for fill, and that ordering
    is a function of the SPARSITY PATTERN ALONE — it never looks at a value. A
    Newton trial re-forms its tangent hundreds of times into one pattern built
    once by :func:`_nr_prepare_assembly`, so COLAMD is re-derived from the same
    structure on every one of them and returns the same permutation every time.

    The first factorization keeps that permutation. Every later one gathers the
    tangent's values into the already-permuted CSC — the pattern's own gather
    index, built once — and factorizes with ``permc_spec='NATURAL'``, which is
    SuperLU doing exactly what it would have done after its own COLAMD call and
    nothing before it. The factors and the solution are BIT-IDENTICAL to the
    unpermuted call: same column order, same partial pivoting, same arithmetic in
    the same sequence, and the fill measured on Griffiths & Lane 1 (435,014
    nonzeros in L+U) and on RS2-65 (6,193,455) is equal to the last entry. What is
    saved is the analysis, which is 21% of the factorization at 2,788 unknowns and
    12% at 27,230 — the numeric work grows faster than the ordering does.
    """
    if _NR_FACTOR_SYMMETRIC:
        return splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                    options=dict(SymmetricMode=True))
    ipc = cache.get("ipc")
    if ipc is not None and cache["sig"] != (K.shape, K.nnz):
        # The pattern this cache was built on is not the pattern in hand. Nothing
        # in the driver rebuilds a pattern under a live solve, so this cannot
        # happen today; it is here because a wrong gather index would be silent.
        ipc = None
    if ipc is None:
        lu = splu(K)
        cache["sig"] = (K.shape, K.nnz)
        pc = lu.perm_c
        cache["pc"] = pc
        cache["ipc"] = ipc = np.argsort(pc)
        indptr = K.indptr
        cache["take"] = np.concatenate(
            [np.arange(indptr[j], indptr[j + 1]) for j in ipc]
        ).astype(np.intp) if K.nnz else np.empty(0, dtype=np.intp)
        newptr = np.zeros(K.shape[1] + 1, dtype=indptr.dtype)
        np.cumsum(np.diff(indptr)[ipc], out=newptr[1:])
        cache["newptr"] = newptr
        cache["newind"] = K.indices[cache["take"]]
        return lu
    take = cache["take"]
    Kp = csc_matrix((K.data[take], cache["newind"], cache["newptr"]),
                    shape=K.shape)
    return _NrPermutedLU(splu(Kp, permc_spec="NATURAL"), cache["pc"])


def _nr_tangent_factorable(K):
    """Is this tangent something SuperLU can be handed at all?

    A consistent Mohr-Coulomb tangent can lose a degree of freedom ENTIRELY: the
    apex branch returns a zero tangent, so a node every one of whose elements has
    yielded to the apex carries an all-zero row and column. That matrix is
    structurally singular, and ``splu`` does not raise on it — it builds a
    degenerate supernode and the process DIES inside OpenBLAS, printing
    "lda must be >= MAX(N,1)" and taking the whole run with it. Measured on
    Griffiths & Lane 1 (quad8, 3.5) at F = 1.8, well past its limit of about 1.37:
    the second tangent re-form of a seeded trial has 98 zero rows and 98 zero
    columns out of 2,788, every entry finite, and SuperLU aborts on it.

    A degree of freedom with no stiffness at all carries exactly the information the
    ``RuntimeError`` path carries — this load is unreachable — so the caller takes
    the same exit. The scan is a bincount and a reduceat over the stored values,
    which costs a small fraction of the factorization it precedes.
    """
    if K.data.size == 0 or np.any(np.diff(K.indptr) == 0):
        return False
    a = np.abs(K.data)
    if not np.all(np.add.reduceat(a, K.indptr[:-1]) > 0.0):
        return False                    # an all-zero column
    return bool(np.all(np.bincount(K.indices, weights=a,
                                   minlength=K.shape[0]) > 0.0))


# Newton-Raphson controls. These are the spike's defaults; every one of them is a
# solver control, not a modeling choice, and none of them changes what the trial
# means — only how hard the solver works before it declares the load unreachable.
_NR_MAX_ITER = 150         # Newton iterations allowed inside ONE load increment
_NR_MIN_STEP = 1.0 / 64    # smallest load increment; below it the load is unreachable
# A load increment is abandoned when it makes NO PROGRESS — no improvement on the
# best residual it has seen for this many iterations — or when the residual climbs
# far above that best. It is NOT abandoned for being slow. An elastoplastic Newton
# solve characteristically sits on a long plateau while the active set churns and
# then drops quadratically: on the dry-dam benchmark at F = 2.4 the relative
# residual holds between 1.4e-1 and 1.3e-1 for ten iterations and then falls
# 4e-3 -> 6e-5 -> 9e-8 -> 1e-12. An earlier rule that required the residual to
# HALVE over six iterations cut exactly that solve off mid-plateau and reported a
# slope that demonstrably stands (the viscoplastic field there is in equilibrium
# and violates the yield surface by 0.2% of the cohesion, which is a lower-bound
# proof of stability) as failed, moving the factor of safety from 2.42 to 1.75.
_NR_NO_PROGRESS = 30       # iterations without a >1% improvement on the best residual
_NR_DIVERGED = 1.0e3       # ... or the residual standing this far above that best
_NR_INIT_STEP = 1.0        # first load increment attempted (1.0 = full gravity in one step)
_NR_GROW = 1.6             # increment growth after a comfortable step
_NR_COMFORT = 8            # ... which is a step converged in this many iterations
_NR_LS_MAX = 9             # line-search backtracks (down to a step of 1/256)
_NR_REL_TOL = 1e-8         # ||r|| / ||f_ext||, the step-level equilibrium test
_NR_TANGENT_H = 1e-7       # strain perturbation for the algorithmic tangent

# ---- the tangent refresh rule, and what it costs (SPIKE, "THE FACTORIZATION") ----
# _nr_equilibrate holds its factorization while the residual keeps falling by at
# least 1/_NR_REFORM_RATIO per iteration, and re-forms when it does not. Measured on
# the representative set, that rule fires on 96% of iterations: the line search
# regularly accepts a fraction of the Newton step, a fractional step does not cut the
# residual by four, and so the tangent is re-formed and re-factorized on nearly every
# iteration. _NR_REFORM_EVERY overrides it with a fixed hold — reform on the first
# iteration of an increment and every Nth after — for measuring what a genuine
# modified-Newton trade costs. None keeps the shipped rule.
_NR_REFORM_RATIO = 0.25
_NR_REFORM_EVERY = None
# The Newton tangent is NON-symmetric at psi = 0, so SuperLU's symmetric mode
# (MMD_AT_PLUS_A ordering, diagonal pivoting) is an approximation on it rather than
# the exact structural statement it is on the viscoplastic driver's elastic
# stiffness. It is 1.8x faster to factorize and 1.5x faster to back-substitute and
# it is NOT bit-identical, so it is a measurement knob and ships off.
_NR_FACTOR_SYMMETRIC = False

# Admissibility of the ANSWER, not a control on the solve. Reaching equilibrium is
# only half of what "the slope stands at F" means; the other half is that the state
# it stands in is one the model can represent. Everything here is small-strain
# kinematics on an undeformed mesh, so a field that has moved a tenth of the model
# height is outside the theory that produced it whatever its residual says. 0.1 is
# solve_fem's own `max_disp_factor` default — the standard a SINGLE solve applies.
# (It is not one the viscoplastic SSRM path enforces: solve_ssrm runs its VP trials
# with `max_disp_factor=None` and leans on the runaway rule instead.)
#
# This is deliberately NOT the gate that was removed (see the long note in the
# increment loop). That one abandoned a load INCREMENT mid-solve on a displacement
# threshold, so it decided verdicts by stopping the solver early. This one is read
# ONCE, on the final fully-loaded equilibrated state, and it changes nothing about
# how hard the solver works — the trial is driven to convergence exactly as before
# and only then asked whether the converged state is admissible.
#
# Measured need, on Griffiths & Lane 1 (quad9, 3.5): F = 1.400 converges to an
# out-of-balance of 3.05e-8 — machine-precision statics, and a worst yield violation
# of 1.5e-8 — at max|u| = 7.693 m on a 50 m model, 15.4% of its height. The
# neighbouring trials read 1.66% (F = 1.3875) and 6.04% (F = 1.396875) of the height:
# the classic displacement upturn, with the limit between the last two. Griffiths &
# Lane read their own Example 1 at about 1.40 the same way, off the upturn.
#
# An explicit non-None `max_disp_factor` from the caller wins; the constant supplies
# the standard on the hybrid/non-convergence path, where solve_ssrm passes None
# because the viscoplastic driver substitutes its own runaway classifier there and
# the Newton driver has none. Setting this to None turns the bound off entirely,
# which is for experiments, not for runs whose answers are quoted.
_NR_DISP_FACTOR = 0.1

# The viscoplastic predictor. A trial that dies at the load-step floor is retried
# once per budget here: a bounded viscoplastic run at the SAME strength on the SAME
# prepared model supplies a displacement field and, more to the point, a plastic
# strain field, and the Newton iteration corrects that instead of walking the load
# path from zero.
#
# It exists because the Newton driver's difficulty on a cohesionless soil is the
# plastic HISTORY and not the equilibrium. Measured on the reinforced specimen
# (docs/fem/files/xslope_reinforce_fem.xlsx, softening unset, tri6/2.0), F = 1.3:
# from zero the driver fails at every load increment down to the 1/64 floor, and
# does so after 18,914 iterations when the per-increment cap and the no-progress
# watch are both lifted, carrying 3% of gravity. Handed a 200-iteration
# viscoplastic state it reaches full gravity in 27 Newton iterations at an
# out-of-balance of 3.7e-7. The root was always there and always Newton-attracting;
# load control from a zero start could not reach it.
#
# What this is NOT is a fallback to the viscoplastic verdict. The predictor's own
# convergence is never read, and the answer is decided entirely by whether the Newton
# corrector reaches full gravity in equilibrium and passes the same force gate, the
# same displacement bound and the same yield reading every other trial passes. A
# trial that is genuinely past failure gets a predictor state that is running away,
# and the corrector refuses it.
#
# Two SHORT rungs, because a state near the limit needs a longer walk than one well
# below it. The measured need on the specimen is 100 iterations at F = 1.3 and 50 at
# F = 1.5; the rungs are set well above that and the second exists for the trials
# that sit hard against the limit.
_NR_VP_PREDICTOR_ITERS = (250, 1000)

# AND ONE ADAPTIVE RUNG, appended. A fixed ladder cannot express "run until this walk
# has finished", which is what a seed on a cohesionless model near its limit needs:
# on the three-layer reinforced variant six FIXED budgets up to 32,000 iterations all
# failed at F = 1.225 while a converged-state seed succeeded in a handful of
# iterations. The last rung is therefore budgeted exactly as a viscoplastic trial of
# the same model at the same strength would be — the caller's own `max_iterations` as
# the chunk, the caller's own `max_iterations_ceiling` as the hard stop, and the
# caller's own `early_failure` — which puts the stopping decision on
# `_still_progressing`, the viscoplastic driver's already-calibrated budget-extension
# rule: continue while the residual's trailing-window mean is falling by at least 1%,
# or while the displacement field is standing still on evidence the failure
# classifier cannot rule on. There is no new dial, and the hard ceiling is the same
# one every viscoplastic trial on this model already respects.
#
# It is APPENDED rather than substituted, and that is deliberate. The short rungs are
# not merely cheaper, they are sometimes better: on the locked tri6/2.0 reinforced
# mesh at F = 1.55625 a 250-iteration seed corrects in 51 Newton iterations while the
# adaptive rung stops one chunk short of its own convergence and the corrector
# refuses it. Keeping the short rungs first and unchanged means every trial the
# driver already rescued is rescued at the same rung from the same state, so the
# extra rung can only convert a FAILED trial, never move a passing one.
#
# What justifies the progress test is a measurement, not its provenance. On the
# three-layer model at the four strengths that decide its bisection, and at three
# chunk sizes, the corrector converged on EVERY seed the adaptive rung had carried to
# its own convergence and on NO seed it had not.
_NR_VP_PREDICTOR_ADAPTIVE = True

# THE ADAPTIVE RUNG RUNS UNDER A RANKINE CAP. With psi = 0 the Mohr-Coulomb flow is
# purely deviatoric, so it cannot relieve a point's MEAN stress; on a c = 0 material
# the yield surface is a cone through the origin, so a tensile point sheds its
# deviatoric stress, stops flowing, and freezes hundreds of psf outside the surface.
# The viscoplastic force gate has no term that could notice (see SPIKE.md, "THE
# THREE-LAYER DISAGREEMENT"), so an uncapped long predictor hands the corrector a
# state that is in force equilibrium but not admissible — and there is no root of the
# Newton residual near it. Measured on the three-layer model at F = 1.20625: the
# UNCAPPED rung converges in 19,716 iterations and the corrector refuses its state;
# the capped rung converges in 22,964 and the corrector takes it in FOUR iterations,
# to an out-of-balance of 1.1e-5 and a yield violation of 1.2e-15.
#
# The cap is on the SEED only. The corrector runs on the caller's own tensile caps
# and reports the caller's own yield reading, so nothing about the verdict is capped.
_NR_VP_PREDICTOR_TENSION_CAP = True

# ===================== The rescue cost policy (SPIKE, "THE COST OF THE RESCUE") ==
#
# The rescue chain above buys the Newton driver its robustness — across the 191-row
# corpus it left ZERO inconclusive trials against the viscoplastic driver's 37 — and
# it buys it at two thirds of that driver's total constitutive work, because it runs
# on every trial that dies at the load-step floor and half of every bisection is a
# trial the search visits only to prove it fails. The three knobs below spend the
# chain where it can change an answer and not where it cannot.
#
# ALL THREE DEFAULT TO OFF, and with them off every path below is the code that was
# there, evaluation for evaluation.
#
# P2 — budget the chain by how far the trial sits above the highest strength this
# search has SHOWN TO STAND, in units of the bracket it started from. A trial that
# far above a standing bound is one the bisection visits to prove a failure; the
# chain there can confirm the cold attempt's verdict but is not measured to overturn
# one, and the adaptive rung costs a whole viscoplastic trial to do it.
#
# Two thresholds, both read off the Phase 0 slice rather than chosen: across 383
# trials the adaptive rung converted a verdict only at or below a QUARTER of the
# initial bracket above the standing bound, and no rung of any length converted one
# above HALF. The values below sit clear of both measured edges. None = off, and the
# whole chain runs on every trial as it did.
_NR_RESCUE_ADAPTIVE_FRAC = None   # above this distance, no adaptive rung
_NR_RESCUE_NONE_FRAC = None       # above this distance, no rescue at all

# P3 — per-model memory. On a model whose trials are carried only from a predictor
# seed, the cold attempt is a walk down to the load-step floor that has already been
# measured to fail at every increment; once one trial has been rescued that way, the
# later trials skip it and start at the first rung.
_NR_SEED_MEMORY = False

# P4 — a cheaper cold attempt: the load-path walk that is going to be handed off
# anyway stops at a coarser increment floor and a shorter per-increment budget. The
# SEEDED correctors keep the shipped control, so nothing that decides a verdict is
# coarsened — only the attempt whose failure is the trigger for the hand-off.
_NR_COLD_CHEAP = False
_NR_COLD_MIN_STEP = 1.0 / 8
_NR_COLD_MAX_ITER = 60

# L4 — the cheap cold attempt as a PRE-FILTER rather than a replacement, which is
# what the cost round's verdict named as the only verdict-safe shape for it. With
# both this and _NR_COLD_CHEAP on, a trial runs the coarse attempt, then the rescue
# chain, and then — before it may be called FAILED — the full cold attempt the coarse
# one stood in for. Nothing is ever refused on less evidence than the driver refuses
# it on today, so the FAILED verdict is safe by construction rather than by a
# threshold fitted to a sample. What it is not is field-identical: a trial the coarse
# attempt or the chain carries reports the state THEY reached, and on a trial the
# full attempt would have carried that is a different converged field.
_NR_COLD_FALLBACK = False


def _nr_predictor_rungs(prep, n_elements, max_iterations, max_iterations_ceiling,
                        early_failure):
    """The predictor's rungs, in the order they are tried.

    Each rung is ``(chunk, ceiling, prepared_model, early_failure)``. The short fixed
    rungs come first, on the caller's own prepared model and with the predictor's own
    convergence never able to end them early; the adaptive rung comes last, on a
    prepared model that is the caller's with its tensile caps replaced by a Rankine
    T = 0. That capped model shares everything else with ``prep`` — the
    factorization, the geometry precompute, the Gauss-point groups are all
    independent of the caps — so building it costs a dict copy.

    See ``_NR_VP_PREDICTOR_ITERS`` and ``_NR_VP_PREDICTOR_ADAPTIVE``.
    """
    rungs = [(int(b), int(b), prep, False) for b in (_NR_VP_PREDICTOR_ITERS or ())]
    if _NR_VP_PREDICTOR_ADAPTIVE:
        prep_seed = prep
        if _NR_VP_PREDICTOR_TENSION_CAP:
            prep_seed = dict(prep)
            # COMPOSED with the caller's own caps rather than replacing them, so a
            # model that already sets t_cut keeps the tighter of the two. A tensile
            # cap is never negative, so this is the same zeros array on every model
            # that sets none — the seed is the caller's model with no tension
            # permitted anywhere, which is what the seed is for.
            prep_seed["t_cap_base"] = np.minimum(
                np.asarray(prep["t_cap_base"], dtype=float), 0.0)
        rungs.append((int(max_iterations),
                      max(int(max_iterations_ceiling or 0), int(max_iterations)),
                      prep_seed, bool(early_failure)))
    return tuple(rungs)


# ===================== Post-peak softening (SPIKE, "POST-PEAK SOFTENING") ====
#
# The step schedule the DROP is walked down on. A bar that has just been latched has
# its cap taken from t_allow to t_res over eta in [0, 1], re-equilibrating at full
# gravity at each accepted step, with the same halve-on-failure control the gravity
# walk uses. A step the corrector cannot carry at the floor is a LIMIT POINT in the
# drop -- the structure cannot shed that force -- which is the same argument this
# driver already makes for gravity, and it is what keeps the verdict a statement
# about the slope rather than about the solver.
#
# The first attempt is the whole drop in one step, because on a structure that can
# carry it that is one equilibrate and the walk costs nothing.
_NR_SOFTEN_INIT_STEP = 1.0
_NR_SOFTEN_MIN_STEP = 1.0 / 64.0


def _nr_soften_newly(bars):
    """Which bars the latch drops now, per group.

    The viscoplastic expression, unchanged: a bar softens when its UNCAPPED elastic
    demand ``k*delta`` -- the quantity both drivers publish as ``forces_1d`` -- exceeds
    the capacity it was allowed to carry, read on a state that has already reached
    full gravity and passed the force gate. Never mid-iteration, never on a failed
    trial, and the set only ever grows.
    """
    return [(~bg['softened'] & bg['can_soften']
             & (bg['_T'] > bg['t_peak'] + 1e-9)) for bg in bars]


def _nr_soften_set_eta(bars, newly, eta):
    """Put the newly-dropping bars' caps at ``eta`` of the way from peak to residual."""
    for bg, nw in zip(bars, newly):
        if nw.any():
            bg['t_cap'] = np.where(
                nw, bg['t_peak'] - eta * (bg['t_peak'] - bg['t_res']), bg['t_cap'])


def _nr_soften_latch(bars, groups, pattern, u, f_ext, free_dofs, n_dof, h_eps,
                     force_tol, oob_fn, nr_max_iter, u_el, piles=None,
                     trans_dofs=None, debug_level=0, deep_free=None):
    """Run the post-peak latch to its fixed point on a converged full-load state.

    Returns ``(ok, u, iterations, force_evaluations, out_of_balance, rounds, cuts)``.
    ``ok`` is False when a drop reaches the step floor, which is a LIMIT POINT in the
    drop: the structure cannot shed that force. Nothing is committed on failure
    beyond the plastic strains the accepted sub-steps already committed, which the
    caller restores if it rejects the state.

    The law is the viscoplastic driver's, read on the same uncapped demand against
    the same threshold, at the same moment — a state in equilibrium with full
    gravity. The PATH is this driver's: the cap walks down rather than jumping, so a
    drop the structure survives but a single step from the pre-drop state cannot
    reach is carried instead of being reported as a failure.
    """
    n_bars = sum(len(bg['idx']) for bg in bars)
    it_total = fe_total = rounds = cuts = 0
    oob_last = 0.0
    while rounds < n_bars:
        newly = _nr_soften_newly(bars)
        if not any(nw.any() for nw in newly):
            break
        rounds += 1
        eta, deta = 0.0, float(_NR_SOFTEN_INIT_STEP)
        while eta < 1.0 - 1e-12:
            eta_try = min(eta + deta, 1.0)
            _nr_soften_set_eta(bars, newly, eta_try)
            ok, u_try, it, fe, oob_here, _rel = _nr_equilibrate(
                groups, pattern, u, f_ext, free_dofs, n_dof, h_eps, force_tol,
                oob_fn, nr_max_iter, u_el, debug_level=debug_level,
                label=f"soften eta={eta_try:.4f}", bars=bars, piles=piles,
                trans_dofs=trans_dofs, deep_free=deep_free)
            it_total += it
            fe_total += fe
            if ok:
                u = u_try
                for grp in groups:
                    grp['_u'] = u
                _nr_commit_plastic_strain(groups)
                eta = eta_try
                oob_last = oob_here
                if it <= _NR_COMFORT:
                    deta = min(1.0, deta * _NR_GROW)
            else:
                cuts += 1
                deta = (eta_try - eta) * 0.5
                if deta < _NR_SOFTEN_MIN_STEP:
                    return False, u_try, it_total, fe_total, oob_here, rounds, cuts
        # The walk finished: those bars now hold the residual, permanently for this
        # solve. The latch is one-way; the FORCE is not, since a softened bar that
        # unloads carries less than its residual.
        for bg, nw in zip(bars, newly):
            if nw.any():
                bg['softened'] |= nw
                bg['t_cap'] = np.where(nw, bg['t_res'], bg['t_cap'])
        if debug_level >= 1:
            _n = sum(int(nw.sum()) for nw in newly)
            _t = sum(int(bg['softened'].sum()) for bg in bars)
            print(f"  Softening round {rounds}: {_n} bar element(s) dropped to "
                  f"t_res ({_t} total); re-solved")
    return True, u, it_total, fe_total, oob_last, rounds, cuts


def _nr_equilibrate(groups, pattern, u_start, f_ext, free_dofs, n_dof, h_eps,
                    force_tol, oob_fn, nr_max_iter, u_elastic_scale,
                    debug_level=0, label="", bars=None, piles=None,
                    trans_dofs=None, deep_free=None):
    """Drive the equilibrium residual to zero at a FIXED external load.

    The inner Newton iteration, lifted out of :func:`_solve_fem_newton` verbatim so
    the two drivers that need it share one copy. The bisection driver calls it once
    per gravity increment with ``f_ext`` a fraction of the load; the ramp driver
    calls it once per strength step with ``f_ext`` the whole load and ``u_start``
    the previous step's converged field.

    Nothing here knows what is being continued. It starts from ``u_start``, uses
    whatever strengths the groups currently carry, and returns

        (ok, u, iterations, force_evaluations, out_of_balance, relative_du)

    without committing plastic strain — the caller decides whether the state is
    worth keeping.

    ``deep_free`` is the ``min_slip_depth`` filter as a boolean over the free
    degrees of freedom (``None`` when the filter is off, which is every unfiltered
    run and leaves every quantity below bit-identical). When it is given, the MERIT
    the iteration is driven and judged on — the residual norm the line search
    minimizes, the no-progress watch, the divergence test, and the relative
    equilibrium test — is taken over the deep degrees of freedom only. That is the
    same convention the out-of-balance reading already carries, extended to the
    place the Newton driver actually decides: this driver's dominant failure is
    "the load increment cannot be carried", and read on the whole model that
    verdict belongs to a cohesionless surface skin the filter exists to exclude.
    The viscoplastic driver never faces the question — its elastic operator is
    non-singular, its convergence test is a RATIO that a creeping skin's own
    displacement flatters, and the one gate that decides is already filtered — so
    a skin runaway is not a verdict there, and it is not one here either.

    The skin is still free to move: nothing bounds it, and its state is carried
    along untouched. It is only barred from ENDING the increment. The one
    concession to arithmetic is the eligibility test in the line search, which
    refuses a candidate whose whole-model residual has grown by more than
    ``_NR_DIVERGED`` in a single step — the driver's own divergence factor, reused
    rather than a new threshold — so a step that is good for the deep system does
    not carry the skin to 1e26 and take the reported state with it.
    """
    def _merit(v_free):
        return float(np.linalg.norm(v_free if deep_free is None
                                    else v_free[deep_free]))

    f_ext_free = f_ext[free_dofs]
    f_norm = max(_merit(f_ext_free), 1e-30)
    n_fe = 0
    rel_du = 0.0
    u_try = u_start.copy()
    ok = False
    it = 0
    r0_norm = None
    oob_here = np.inf
    r_hist = []
    r_best = float('inf')
    last_progress = 0
    lu = None                  # the live tangent factorization (see below)
    # The column ordering COLAMD derives from this pattern, kept for the whole
    # trial (see _nr_factorize_tangent). The pattern is built once per solve, so
    # this is the same structure every re-form and every load increment sees.
    _order_cache = pattern.setdefault("_order", {})
    prev_r = None
    for it in range(1, nr_max_iter + 1):
        # Re-form and re-factorize the tangent only when the previous step did
        # not earn its keep. The consistent tangent changes only where the
        # active set changes, so once the plastic zone has settled the SAME
        # factorization drives the remaining iterations at essentially the same
        # rate for a fraction of the cost — the assembly and the LU are the
        # whole per-iteration expense, the return map is not. A step whose
        # residual stops falling by 4x per iteration gets a fresh tangent.
        if _NR_REFORM_EVERY is None:
            reform = (lu is None or prev_r is None
                      or r_hist[-1] > _NR_REFORM_RATIO * prev_r)
        else:
            reform = lu is None or (it - 1) % int(_NR_REFORM_EVERY) == 0
        _tp = time.perf_counter() if _PROF_ON else None
        fint, tangents, _ = _nr_internal_force(groups, u_try, n_dof,
                                               h_eps=h_eps, want_tangent=reform,
                                               bars=bars, piles=piles)
        if _PROF_ON:
            _prof_add("nr_const_tan" if reform else "nr_const_res", _tp)
        n_fe += 1
        r = f_ext - fint
        r_free = r[free_dofs]
        r_norm = _merit(r_free)
        if not np.isfinite(r_norm):
            break
        r_full = np.zeros(n_dof)
        r_full[free_dofs] = r_free
        oob_here = oob_fn(r_full)
        if r0_norm is None:
            r0_norm = max(r_norm, 1e-30)
        # Equilibrium for this increment. Either test alone is enough: the
        # relative norm is the Newton-natural one, and the Dawson measure at a
        # tenth of the trial tolerance guarantees the final state passes the
        # same gate the viscoplastic verdict is read on.
        if r_norm / f_norm < _NR_REL_TOL or oob_here < 0.1 * force_tol:
            ok = True
            break
        prev_r = r_hist[-1] if r_hist else None
        r_hist.append(r_norm)
        if r_norm < 0.99 * r_best:
            r_best = r_norm
            last_progress = it
        if it - last_progress > _NR_NO_PROGRESS or r_norm > _NR_DIVERGED * r_best:
            break                       # no progress: this load is unreachable

        if reform:
            _tp = time.perf_counter() if _PROF_ON else None
            K = _nr_assemble_tangent(groups, tangents, pattern, bars=bars,
                                     piles=piles)
            if _PROF_ON:
                _prof_add("nr_assemble", _tp)
            if not _nr_tangent_factorable(K):
                break                   # structurally singular = the limit load
            _tp = time.perf_counter() if _PROF_ON else None
            try:
                lu = _nr_factorize_tangent(K, _order_cache)
            except RuntimeError:
                break                   # singular tangent = the limit load
            finally:
                if _PROF_ON:
                    _prof_add("nr_factorize", _tp)
        _tp = time.perf_counter() if _PROF_ON else None
        du_free = lu.solve(r_free)
        if _PROF_ON:
            _prof_add("nr_trisolve", _tp)
        if not np.all(np.isfinite(du_free)):
            break
        du = np.zeros(n_dof)
        du[free_dofs] = du_free

        # Line search: accept the full Newton step unless it makes the
        # residual worse, then backtrack. Near the limit load the full step
        # regularly overshoots, and halving is what keeps the iteration in the
        # basin instead of throwing it into a spurious runaway.
        #
        # The step actually taken is always one that was MEASURED, and the best
        # of those measured. Taking an unmeasured fallback step — what this did
        # when every backtrack failed — is how the dry-dam solve blew up: with a
        # near-singular tangent the correction was enormous, no tested step
        # reduced the residual, the untested one was applied anyway, and the
        # displacement reached 4.6e13 times the elastic response in a single
        # iteration. A step that cannot be shown to help is a step that has to
        # stay small.
        alpha = 1.0
        best_alpha, best_rc = None, np.inf
        g_now = (None if deep_free is None
                 else max(float(np.linalg.norm(r_free)), 1e-30))
        for _ls in range(_NR_LS_MAX):
            cand = u_try + alpha * du
            _tp = time.perf_counter() if _PROF_ON else None
            f_c, _, _ = _nr_internal_force(groups, cand, n_dof, bars=bars,
                                           piles=piles)
            if _PROF_ON:
                _prof_add("nr_linesearch", _tp)
            n_fe += 1
            rc_free = (f_ext - f_c)[free_dofs]
            rc = _merit(rc_free)
            eligible = bool(np.isfinite(rc))
            if eligible and g_now is not None:
                # See the docstring: the skin may run away, but not by a factor of
                # a thousand in one accepted step.
                eligible = float(np.linalg.norm(rc_free)) <= _NR_DIVERGED * g_now
            if eligible and rc < best_rc:
                best_rc, best_alpha = rc, alpha
            if eligible and rc < r_norm:
                break
            alpha *= 0.5
        alpha = best_alpha if best_alpha is not None else 0.0
        if _PROF_ON:
            # How much of the Newton step the search actually takes, and how often
            # it runs out of backtracks. Both are bookkeeping.
            _PROF["alpha_sum"] = _PROF.get("alpha_sum", 0.0) + alpha
            _PROF["n_alpha"] = _PROF.get("n_alpha", 0) + 1
            _PROF["n_ls_full"] = _PROF.get("n_ls_full", 0) + (1 if alpha == 1.0 else 0)
            _PROF["n_ls_exhausted"] = (_PROF.get("n_ls_exhausted", 0)
                                       + (1 if _ls == _NR_LS_MAX - 1 else 0))
        u_try = u_try + alpha * du
        # Translational degrees of freedom only: a rotation is not a length, and
        # this ratio is reported as `residual` on both drivers. Without a pile every
        # degree of freedom is translational and this is the same number as before.
        _umax = _nr_umax(u_try, trans_dofs)
        rel_du = ((_nr_umax(alpha * du, trans_dofs) / _umax)
                  if _umax > 0.0 else 0.0)
        if debug_level >= 3:
            print(f"      NR {label} it={it:3d} "
                  f"||r||/||f||={r_norm / f_norm:.3e} oob={oob_here:.2e} "
                  f"alpha={alpha:.3g} max|u|={np.max(np.abs(u_try)):.4g} "
                  f"({np.max(np.abs(u_try)) / max(u_elastic_scale, 1e-30):.1f}x "
                  f"elastic){' [tangent re-formed]' if reform else ''}")

        # There is deliberately NO displacement gate here. An earlier version
        # abandoned an increment once max|u| passed 50x the elastic response, on
        # the reasoning that a slope moving that far is running away. It is not a
        # statement about the slope: at the limit load the load-displacement curve
        # flattens, so the last standing increments are exactly the ones with the
        # largest displacements, and the gate cut them off. Measured on Griffiths
        # & Lane 1 (quad9, 3.5) at F = 1.400: with the gate the trial is reported
        # FAILED after 440 iterations; without it the same driver reaches
        # equilibrium in 48, at an out-of-balance of 3e-8 and with no Gauss point
        # more than 1.5e-8 of its strength outside the yield surface -- a
        # statically admissible field, which is the definition of the trial
        # standing. One trial's verdict, and the factor of safety moved 1.3969 ->
        # 1.4031. A load that cannot be carried already proves itself by driving
        # the increment below its floor, which is a statement about the slope; a
        # displacement threshold is a statement about a number nobody chose for a
        # reason. The runaway is still caught, just not by a threshold: a
        # non-finite residual, a singular tangent and the no-progress watch all
        # end the increment, and the failing trials in the benchmark table all
        # terminate at the load-step floor without it.

    return ok, u_try, it, n_fe, oob_here, rel_du


def _solve_fem_newton(fem_data, F, prep, *, c_reduced, phi_reduced,
                      elastic_by_elem, t_cap_by_elem, finv_by_elem,
                      force_tol, min_slip_depth, max_iterations,
                      max_disp_factor=None, _nr_export=None,
                      debug_level=0, progress_callback=None,
                      nr_max_iter=_NR_MAX_ITER, nr_min_step=_NR_MIN_STEP,
                      _nr_seed=None, k0=None, _nr_init_state=None,
                      _nr_seed_state=None, _nr_env_F=None, _softened_seed=None):
    """One strength-reduction trial by Newton-Raphson (SPIKE; see SPIKE.md).

    Returns the same result dictionary solve_fem returns, so the SSRM bisection
    consumes it unchanged. The verdict is BINARY by construction — the full
    gravity load is either reached in equilibrium (`converged`, `stable`) or the
    load increment is driven below its floor without reaching it (`FAILED`,
    exit_reason 'diverging') — so this path can never leave a trial
    'inconclusive'.

    A trial that reaches full gravity in equilibrium is checked once more, against
    the displacement bound (`_NR_DISP_FACTOR`): a converged state that has moved
    more than a tenth of the model height is reported FAILED with exit_reason
    'displacement_limit', which is distinguishable from the force-side failures
    ('diverging' at the load-step floor, 'iteration_cap' at the final force gate).

    Deliberately narrow: Mohr-Coulomb, psi = 0, no strain softening. The Rankine
    tension cutoff IS carried, as a second yield surface combined with the shear
    one by Koiter's rule (see mc_return_map and _nr_rankine_return), which is the
    same physics the viscoplastic path gets by summing the two flows.
    Reinforcement bar elements ARE carried (see _nr_build_bars and _nr_bar_force),
    with the soil strength reduced and the bars keeping their capacity, which is the
    vendor convention and what the LEM does. Pile beam elements ARE carried (see
    _nr_build_piles and _nr_pile_force), rotational degrees of freedom and all, with
    their rigidity and their structural capacity likewise held while only the soil
    strength falls. The K0 initial stress IS carried (see
    _nr_group_sig0): the in-situ field rides inside the internal force rather than
    being moved to the right-hand side, because Newton has an internal force to
    write, and `_nr_init_state` starts a trial from the equilibrated in-situ state
    with its displacement measured from there — the same sequencing solve_ssrm
    gives the viscoplastic driver. Everything else the spike does not cover raises
    rather than being silently ignored, because a switch that quietly drops a
    tensile cap or a softening bar would return a wrong factor of safety that looks
    right.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    n_elements = len(elements)
    n_dof = prep["n_dof"]
    free_dofs = prep["free_dofs"]
    n_free = prep["n_free"]
    dof_offset = fem_data.get("dof_offset", None)

    # ---- what this path does NOT cover -------------------------------------
    unsupported = []
    # Pile beam elements ARE carried (see _nr_build_piles / _nr_pile_force): the same
    # element matrix the global elastic stiffness is assembled from, the same
    # per-unit-width rigidities and capacities, and the same head/tip fixity, which
    # is a constraint on `prep["free_dofs"]` and therefore already in force on this
    # path. Nothing about a pile is reduced by F, as nothing about a bar is.
    # Reinforcement bars ARE carried (see _nr_build_bars), and so is their POST-PEAK
    # drop (see _nr_soften). A bar with a finite t_res below the capacity its
    # embedment develops drops to that residual through the same converged-state
    # latch the viscoplastic path applies — read on the same uncapped demand, against
    # the same threshold — with the drop walked down as a continuation on the cap so
    # that a drop the structure cannot carry reads as a limit point rather than as a
    # solver that gave up. CONTINUUM strain-softening is a different animal and stays
    # refused on both drivers; there is no input for it, so there is nothing to guard.
    # The two CURVED envelopes ARE carried (see _nr_envelope_return and
    # _nr_envelope_by_elem), per GAUSS POINT, so one mesh can hold a Mohr-Coulomb
    # soil and a Hoek-Brown rock and solve them together. Both arrive here the way
    # the viscoplastic driver carries them — as a Mohr-Coulomb tangent to the
    # F-reduced envelope, re-derived from the current stress — because neither
    # driver has a curved yield surface and the two have to solve the same model.
    # Matric suction IS carried (see _nr_group_suction / _nr_group_apparent_
    # cohesion). Fredlund's credit is an apparent cohesion on the same linear
    # envelope, so it arrives here folded into the per-Gauss-point c_r that the
    # return map already takes, reduced by the trial F exactly as c' is. Nothing
    # about its semantics is decided on this path.
    # The Rankine tension cutoff IS carried (see mc_return_map's second surface and
    # _nr_rankine_return). It arrives here as t_cap_by_elem, already reduced by the
    # trial F where `tension_srf` says so, so nothing about the cap's semantics is
    # decided on this path — it is the same per-element array the viscoplastic
    # driver caps sigma_1 against.
    # The K0 initial stress IS carried (see _nr_group_sig0): the same overburden
    # integral the viscoplastic path reads, the same in-plane AND out-of-plane
    # sigma_h = K0 sigma'_v, and the same in-situ pre-equilibration sequencing
    # through `_nr_init_state`. It arrives here as prep['sv0_gp'] plus the k0
    # argument, so nothing about the field's semantics is decided on this path.
    if unsupported:
        raise NotImplementedError(
            "fem_solver='newton' is a plain Mohr-Coulomb spike and does not "
            "handle " + ", ".join(unsupported) + ". Run this model on the "
            "default viscoplastic solver.")

    # The displacement bound (_NR_DISP_FACTOR), the elastic displacement scale, the
    # tangent probe's step and the max|du|/max|u| both drivers report as `residual`
    # are all LENGTHS, so they are read over the translational degrees of freedom
    # only — which is what the viscoplastic driver's own displacement-limit check
    # reads. A pile node carries a rotation as its third degree of freedom; putting
    # one into any of those maxima would compare a length against a radian. Without a
    # pile this index is None and every one of them is the plain max over the whole
    # vector, bit for bit.
    trans_dofs = _nr_translational_dofs(fem_data, n_dof)

    env_by_elem = _nr_envelope_by_elem(fem_data, _nr_env_F)
    groups = _nr_build_groups(prep, c_reduced, phi_reduced, elastic_by_elem,
                              t_cap_by_elem, k0=k0, finv_by_elem=finv_by_elem,
                              env_by_elem=env_by_elem)
    bars = _nr_build_bars(fem_data)
    # A carried-in post-peak set. The bisection never passes one — every trial is a
    # cold solve with an empty set, which is the viscoplastic convention — so this is
    # the retry path (a seeded corrector picks up the set its cold attempt had
    # reached) and the ramp (a bar that ruptured at one strength is ruptured at the
    # next). Softening only ever grows, so the seed is a lower bound and never a
    # verdict.
    if bars is not None and _softened_seed is not None:
        _seed = np.asarray(_softened_seed, dtype=bool)
        if _seed.size == len(fem_data.get("elements_1d", ())):
            for bg in bars:
                take = _seed[bg['idx']] & bg['can_soften']
                if take.any():
                    bg['softened'] |= take
                    bg['t_cap'] = np.where(take, bg['t_res'], bg['t_cap'])
    piles = _nr_build_piles(fem_data)
    pattern = _nr_prepare_assembly(groups, free_dofs, n_dof, bars=bars, piles=piles)

    # ---- external load -------------------------------------------------------
    # Effective-stress formulation, exactly as the viscoplastic path builds it:
    # the pore-pressure term joins the load vector and D*(B u - eps^p) is then the
    # EFFECTIVE stress the yield check reads directly.
    base_loads = prep["F_gravity"].copy()
    F_u = np.zeros(n_dof)
    any_u = False
    for grp, sg in zip(groups, prep["gp_groups_static"]):
        u_gp = np.array([prep["u_gp"][e][g] for e, g in sg['pairs']])
        if np.any(u_gp):
            any_u = True
        contrib = (grp['w'] * u_gp)[:, None] * (grp['B'][:, 0, :] + grp['B'][:, 1, :])
        np.add.at(F_u, grp['dof'], contrib)
    if any_u:
        base_loads = base_loads + F_u

    node_dof_x = prep["node_dof_x"]
    node_dof_y = prep["node_dof_y"]
    node_has_free = prep["node_has_free"]
    g_node_den = prep["g_node_den"]
    deep_free_mask = prep["deep_free_mask"]
    free_dof_mask = prep["free_dof_mask"]
    # The same filter on the degrees of freedom, restricted to the free ones in
    # the order the residual vector carries them. None without a filter, and every
    # quantity it touches is then exactly what it was.
    deep_dof_mask = prep.get("deep_dof_mask")
    deep_free = None if deep_dof_mask is None else deep_dof_mask[free_dofs]

    def _oob(r_full):
        """The Dawson, Roth & Drescher measure, on the TRUE residual.

        Same quantity the viscoplastic path reports — the largest nodal
        out-of-balance force divided by that node's own tributary weight — but
        read off the genuine equilibrium residual rather than off the increment of
        a body-load correction, because Newton-Raphson has one to read.
        """
        d = r_full * free_dof_mask
        rn = np.sqrt(d[node_dof_x] ** 2 + d[node_dof_y] ** 2)
        v = (rn / g_node_den)[node_has_free]
        if deep_free_mask is not None:
            v = v[deep_free_mask]
        return float(np.max(v)) if v.size else 0.0

    # Elastic reference response (the displacement scale the runaway test and the
    # reported u_ratio are measured in — the same yardstick the viscoplastic path
    # uses, so the two solvers' metadata are comparable).
    u_elastic = np.zeros(n_dof)
    u_elastic[free_dofs] = prep["K_factor"].solve(base_loads[free_dofs])
    u_elastic_scale = _nr_umax(u_elastic, trans_dofs) if n_free else 0.0
    h_eps = _NR_TANGENT_H * max(
        1e-12, _nr_umax(u_elastic, trans_dofs) / max(1e-12, prep["mesh_height"]))

    ext_norm = float(np.linalg.norm(base_loads[free_dofs]))
    u = np.zeros(n_dof)
    # A seeded start (see _NR_VP_PREDICTOR_ITERS): the plastic strain field and the
    # displacements come from a bounded viscoplastic run at this same strength, on
    # this same prepared model, so the load path is already walked and this call has
    # only to correct it. `plastic_strains` carries the same quantity on both
    # drivers — eps^p per Gauss point, [ex, ey, gxy, ez], with
    # sigma = D (B u - eps^p) — which is what makes the hand-over meaningful rather
    # than a shared array name.
    if _nr_seed is not None:
        u_seed, ep_seed = _nr_seed
        u = np.asarray(u_seed, dtype=float)[:n_dof].copy()
        for grp in groups:
            ep = np.empty_like(grp['ep'])
            for k, (e, gp) in enumerate(grp['pairs']):
                ep[k] = ep_seed[e][gp]
            grp['ep'] = ep
    # The equilibrated in-situ state (K0 runs only; see solve_ssrm's equilibration
    # solve). It is the displacement field and the plastic strain at the end of a
    # full-strength solve, in this prepared model's own group order, and starting
    # here starts this trial from that solve's answer. Displacement is then measured
    # FROM it — the reported field and the displacement bound both — because the
    # in-situ displacement is the artifact of imposing a stress field the geometry
    # does not hold in equilibrium, not travel the soil made. Stresses stay
    # functions of the ABSOLUTE displacement, which is what the solve iterates on.
    u_datum = np.zeros(n_dof)
    if _nr_init_state is not None:
        if not groups or groups[0].get('sig0') is None:
            raise ValueError("_nr_init_state was given without k0; an equilibrated "
                             "initial state has no meaning without the K0 "
                             "formulation.")
        _sizes = [len(g['pairs']) for g in groups]
        _iu = np.asarray(_nr_init_state["u"], dtype=float)
        _iev = _nr_init_state["evp"]
        if [len(a) for a in _iev] != _sizes or _iu.shape != (n_dof,):
            raise ValueError(
                "_nr_init_state does not match this prepared model's Gauss-point "
                "groups / degrees of freedom; the state must be produced on the "
                "same prepared model.")
        u_datum = _iu.copy()
        u = u_datum.copy()
        for grp, ev in zip(groups, _iev):
            grp['ep'] = np.array(ev, dtype=float, copy=True)
    # A predictor's state on a K0 model, in the same (u, eps^p) form the in-situ
    # state carries. It moves where the solve STARTS and not what displacement is
    # measured from: the datum stays the in-situ state, so a predictor's own travel
    # is not quietly credited to the trial. (The `_nr_seed` form above is the
    # no-K0 route, where the viscoplastic driver's reported displacement IS the
    # absolute one and its plastic strain comes back per element.)
    if _nr_seed_state is not None:
        u = np.asarray(_nr_seed_state["u"], dtype=float)[:n_dof].copy()
        for grp, ev in zip(groups, _nr_seed_state["evp"]):
            grp['ep'] = np.array(ev, dtype=float, copy=True)
    # A trial that starts from a state which already stands at full gravity does not
    # walk the load path — the in-situ solve (or the predictor) already walked it,
    # and cutting the load below full gravity from such a state would be solving a
    # different problem.
    _one_shot = (_nr_seed is not None or _nr_init_state is not None
                 or _nr_seed_state is not None)
    lam = 0.0
    dlam = float(_NR_INIT_STEP)
    total_iterations = 0
    n_steps = 0
    n_cuts = 0
    step_iters = []
    exit_reason = 'converged'
    converged = True
    last_oob = 0.0
    # Work, honestly counted. An "iteration" here buys one residual evaluation AND
    # up to `_NR_LS_MAX` more inside the line search, so an iteration count compared
    # against the viscoplastic driver's — where an iteration is one constitutive
    # pass — understates Newton's constitutive work by up to a factor of ten. Every
    # call to _nr_internal_force in this solve is counted, including the two
    # reporting passes at the end.
    n_force_evals = 0
    # The viscoplastic driver publishes max|du| / max|u| as `residual`. This is the
    # same quantity on this path, recorded at the last accepted load increment, so
    # the key means one thing on both drivers.
    last_rel_du = 0.0
    rel_du = 0.0

    while lam < 1.0 - 1e-12:
        # A seeded trial does not walk the load path — the predictor already did,
        # and cutting the load below full gravity from a state that stands at full
        # gravity would be solving a different problem. It is one attempt at the
        # whole load, and it either corrects the seed or the trial has failed.
        step = 1.0 - lam if _one_shot else min(dlam, 1.0 - lam)
        lam_try = lam + step
        f_ext = lam_try * base_loads
        # The K0 initial stress is gravity's own field, so it rides the load factor
        # with the load: at lam the material carries lam*sigma_0 and lam of the
        # weight, and lam = 1 is the authored in-situ state. A no-op without K0.
        _nr_set_load_factor(groups, lam_try)

        ok, u_try, it, _fe, oob_here, rel_du = _nr_equilibrate(
            groups, pattern, u, f_ext, free_dofs, n_dof, h_eps, force_tol,
            _oob, nr_max_iter, u_elastic_scale, debug_level=debug_level,
            label=f"lam={lam_try:.4f}", bars=bars, piles=piles,
            trans_dofs=trans_dofs, deep_free=deep_free)
        n_force_evals += _fe
        total_iterations += it
        if ok:
            u = u_try
            for grp in groups:
                grp['_u'] = u
            _nr_commit_plastic_strain(groups)
            lam = lam_try
            last_oob = oob_here
            last_rel_du = rel_du
            n_steps += 1
            step_iters.append(it)
            if it <= _NR_COMFORT:
                dlam = min(1.0, dlam * _NR_GROW)
            if debug_level >= 2:
                print(f"    NR load factor {lam:.4f} reached in {it} iterations "
                      f"(oob {oob_here:.2e})")
            if progress_callback is not None:
                try:
                    progress_callback(lam, f"newton lambda={lam:.3f}, {it} iters")
                except Exception:
                    pass
        else:
            n_cuts += 1
            dlam = step * 0.5
            if debug_level >= 2:
                print(f"    NR increment to {lam_try:.4f} failed after {it} "
                      f"iterations; cutting to {dlam:.5f}")
            if _one_shot:
                converged = False
                exit_reason = 'diverging'
                u = u_try
                break
            if dlam < nr_min_step:
                converged = False
                exit_reason = 'diverging'
                # Report the best state actually reached, which is the last
                # equilibrated one plus whatever the failing attempt reached.
                u = u_try
                break

    # A trial that reached full load must still pass the SAME force-equilibrium
    # gate the viscoplastic verdict is read on, or it is not converged.
    if converged:
        fint, _, _ = _nr_internal_force(groups, u, n_dof, bars=bars, piles=piles)
        n_force_evals += 1
        r_full = np.zeros(n_dof)
        r_full[free_dofs] = (base_loads - fint)[free_dofs]
        last_oob = _oob(r_full)
        if last_oob >= force_tol:
            converged = False
            exit_reason = 'iteration_cap'

    # ---- post-peak softening: the viscoplastic latch, on a stepped continuation --
    #
    # The state is in equilibrium at full gravity. NOW, and only now, ask which bars
    # have actually yielded -- their uncapped elastic demand exceeds the capacity they
    # were allowed to carry -- and drop those with a finite residual to it. Shedding
    # their load can push neighbours over, so the process repeats until the set stops
    # growing; it terminates because the set only grows and is bounded by the element
    # count. That is the viscoplastic driver's own fixed point, read on the same
    # quantity against the same threshold.
    #
    # What is different here is the PATH, not the law. The viscoplastic driver applies
    # the whole drop at once and relaxes; Newton is given a residual perturbed by a
    # finite amount with no small parameter to continue in, and a drop the structure
    # survives can sit outside the basin of a single step from the pre-drop state. So
    # the cap is walked down from peak to residual and the step is halved when the
    # corrector cannot carry it. Reaching the floor IS the verdict: the structure
    # cannot shed that force. The end state at eta = 1 is the same law either way.
    n_soften_rounds = 0
    if converged and bars is not None and any(bg['can_soften'].any() for bg in bars):
        _ok, u, _it, _fe, _oobs, n_soften_rounds, _cuts = _nr_soften_latch(
            bars, groups, pattern, u, base_loads, free_dofs, n_dof, h_eps,
            force_tol, _oob, nr_max_iter, u_elastic_scale, piles=piles,
            trans_dofs=trans_dofs, debug_level=debug_level,
            deep_free=deep_free)
        total_iterations += _it
        n_force_evals += _fe
        n_cuts += _cuts
        if not _ok:
            converged = False
            exit_reason = 'diverging'
        # The state that came out of the last drop must pass the SAME force gate.
        if converged:
            fint, _, _ = _nr_internal_force(groups, u, n_dof, bars=bars, piles=piles)
            n_force_evals += 1
            r_full = np.zeros(n_dof)
            r_full[free_dofs] = (base_loads - fint)[free_dofs]
            last_oob = _oob(r_full)
            if last_oob >= force_tol:
                converged = False
                exit_reason = 'iteration_cap'

    # ... and the same displacement bound. Force equilibrium alone is not the
    # question; equilibrium in a state the small-strain model can represent is.
    # See _NR_DISP_FACTOR. Read once, here, on the final state.
    disp_factor = (max_disp_factor if max_disp_factor is not None
                   else _NR_DISP_FACTOR)
    nr_disp_limit = (disp_factor * prep["mesh_height"]
                     if disp_factor is not None and prep["mesh_height"] > 0
                     else None)
    # With a min_slip_depth filter the bound is read on the DEEP degrees of freedom
    # only, for the reason the out-of-balance reading already is: the skin the filter
    # excludes is allowed to fail, and a state whose deep field is small-strain
    # admissible is not disqualified by how far that skin has travelled. The skin's
    # own travel is reported beside it rather than dropped.
    nr_disp_deep, nr_disp_skin = _nr_umax_split(u - u_datum, trans_dofs,
                                               deep_dof_mask)
    if converged and nr_disp_limit is not None:
        if nr_disp_deep > nr_disp_limit:
            converged = False
            exit_reason = 'displacement_limit'

    # ---- reporting ----------------------------------------------------------
    for grp in groups:
        grp['_u'] = u
    # refresh _sig / _branch (and the bar and pile forces) at the reported state
    _fint_rep, _, _ = _nr_internal_force(groups, u, n_dof, bars=bars, piles=piles)
    n_force_evals += 1
    # The out-of-balance at the reported state, over the WHOLE model and over the
    # skin the min_slip_depth filter excludes. Reported, never read for a verdict:
    # `unbalanced_force_ratio` stays the filtered quantity both drivers publish, and
    # this says how much of the model that quantity is not looking at. Without a
    # filter the skin reading is None and the global one is the same number.
    _r_rep = np.zeros(n_dof)
    _r_rep[free_dofs] = (base_loads - _fint_rep)[free_dofs]
    nr_oob_global, nr_oob_skin = _nr_oob_split(
        _r_rep, free_dof_mask, node_dof_x, node_dof_y, node_has_free, g_node_den,
        deep_free_mask)

    # ---- reinforcement diagnostics ------------------------------------------
    # The same three arrays the viscoplastic path returns, at full 1D-element
    # length so every reader indexes them the same way: the force each bar actually
    # delivers (the elastic k*delta clipped into [0, t_allow], which is what
    # forces_1d means on both drivers), which bars are AT their capacity, and the
    # post-peak set, which is what the latch above wrote.
    #
    # One convention difference is deliberate and is not a defect on either side.
    # The viscoplastic driver's failed mask LATCHES: it records every bar that
    # exceeded its capacity at any point in the iteration history, including the
    # elastic predictor's overshoot before the soil sheds load into the bars. This
    # one is read on the reported state, so it says which bars are at capacity in
    # the field being exported. The Newton mask is therefore a subset.
    n_1d_total = len(fem_data.get("elements_1d", np.array([]).reshape(0, 3)))
    forces_1d_out = np.zeros(n_1d_total)
    failed_1d_out = np.zeros(n_1d_total, dtype=bool)
    softened_1d_out = np.zeros(n_1d_total, dtype=bool)
    if bars is not None:
        for bg in bars:
            forces_1d_out[bg['idx']] = bg['_T_true']
            failed_1d_out[bg['idx']] = bg['_T'] > bg['t_cap'] + 1e-9
            softened_1d_out[bg['idx']] = bg['softened']

    # ---- pile diagnostics ---------------------------------------------------
    # The same five arrays the viscoplastic path returns, indexed by pile element
    # in the same order, so every reader — the summary printer, the result CSVs,
    # the pile-shear colorbar — consumes them unchanged. The forces reported are
    # the ones the element actually delivers, which is the capped action, exactly
    # as `forces_1d` reports the capped bar force. The yielded masks are read on
    # the REPORTED state, where the viscoplastic driver's latch every element that
    # was ever over its capacity at any point in the iteration history, so this
    # mask is a subset by construction — the same convention difference the bar
    # masks carry, and for the same reason.
    n_pile_out = int(fem_data.get("n_pile_elements", 0))
    pile_axial = np.zeros(n_pile_out)
    pile_shear = np.zeros(n_pile_out)
    pile_moment = np.zeros((n_pile_out, 2))
    pile_prot = np.zeros((n_pile_out, 2))
    pile_yV = np.zeros(n_pile_out, dtype=bool)
    pile_yM = np.zeros(n_pile_out, dtype=bool)
    if piles is not None:
        for pg in piles:
            i = pg['idx']
            pile_axial[i] = pg['_axial']
            pile_shear[i] = pg['_V_true']
            # The DELIVERED end moments, read on the released displacement, so a
            # hinged end reports the capacity because the equilibrium carries it
            # and not because the report was clipped.
            pile_moment[i] = pg['_M']
            pile_prot[i] = pg['_p_rot']
            pile_yV[i] = pg['_yV']
            pile_yM[i] = pg['_yM']

    # ---- the verdict's own evidence -----------------------------------------
    # A converged Newton trial asserts two things about the slope: that full
    # gravity is carried in equilibrium, and that no Gauss point is outside the
    # yield surface. `unbalanced_force_ratio` already carries the first as the
    # Dawson out-of-balance. This carries the second — the largest yield-function value
    # over every Gauss point, divided by that point's own strength scale, so it
    # reads as a fraction of the strength available there.
    #
    # It is computed from the INVARIANT form of the Mohr-Coulomb function that the
    # viscoplastic path uses, not from the ordered-principal-stress form the return
    # map is written on. The two are the same surface algebraically, so a defect in
    # one cannot hide behind the other, and a converged trial that reports a
    # violation near machine precision is a statically admissible stress field
    # rather than a solver's word for one.
    sq3_ = np.sqrt(3.0)
    max_yield_violation = 0.0
    max_tension_violation = None
    for grp in groups:
        sg = grp['_sig']
        sx_, sy_, txy_, sz_ = sg[:, 0], sg[:, 1], sg[:, 2], sg[:, 3]
        sigm_ = (sx_ + sy_ + sz_) / 3.0
        dsbar_ = np.sqrt(((sx_ - sy_) ** 2 + (sy_ - sz_) ** 2 + (sz_ - sx_) ** 2
                          + 6.0 * txy_ ** 2) / 2.0)
        dx_, dy_, dz_ = sx_ - sigm_, sy_ - sigm_, sz_ - sigm_
        ds3_ = np.maximum(dsbar_, 1e-10) ** 3
        sine_ = np.clip(np.where(dsbar_ > 1e-10,
                                 -13.5 * (dx_ * dy_ * dz_ - dz_ * txy_ ** 2) / ds3_,
                                 0.0), -1.0, 1.0)
        th_ = np.arcsin(sine_) / 3.0
        # The envelope the trial was actually solved on. On a Mohr-Coulomb group
        # these three ARE grp['c_r'] / ['snph'] / ['csph'], the same objects, so
        # this reading is unchanged there; on a curved-envelope group they are the
        # converged linearization the return was taken on, which is the only
        # envelope against which "how far outside" means anything.
        _ce_ = grp.get('_c_eff', grp['c_r'])
        _se_ = grp.get('_snph_eff', grp['snph'])
        _cse_ = grp.get('_csph_eff', grp['csph'])
        fv_ = (sigm_ * _se_
               + dsbar_ * (np.cos(th_) / sq3_ - np.sin(th_) * _se_ / 3.0)
               - _ce_ * _cse_)
        # Strength scale at the point: the two terms the deviatoric radius is held
        # against. A material held linear elastic carries c = inf and is skipped,
        # as is a point with no strength scale to divide by.
        den_ = _ce_ * _cse_ + np.abs(sigm_) * _se_
        ok_ = np.isfinite(den_) & (den_ > 0.0)
        if np.any(ok_):
            max_yield_violation = max(max_yield_violation,
                                      float(np.max(fv_[ok_] / den_[ok_])))
        # The tensile half of the same reading: how far the major principal stress
        # sits above the cap, on the same scale. Computed from the components
        # rather than from the return map's ordered principals, for the same reason
        # the shear reading is — an independent form of the same statement.
        _tc = grp.get('t_cap')
        if _tc is not None:
            _m = np.isfinite(_tc) & ok_
            if np.any(_m):
                _ctr = 0.5 * (sx_ + sy_)
                _r = np.sqrt((0.5 * (sx_ - sy_)) ** 2 + txy_ ** 2)
                _s1 = np.maximum(_ctr + _r, sz_)
                _tv = float(np.max((_s1[_m] - _tc[_m]) / den_[_m]))
                max_tension_violation = (
                    _tv if max_tension_violation is None
                    else max(max_tension_violation, _tv))
    if max_tension_violation is not None:
        max_yield_violation = max(max_yield_violation, max_tension_violation)

    sig_by_gp = [[None] * len(prep["elem_gp_data"][e]) for e in range(n_elements)]
    branch_by_gp = [[0] * len(prep["elem_gp_data"][e]) for e in range(n_elements)]
    ep_by_gp = [[np.zeros(4)] * len(prep["elem_gp_data"][e]) for e in range(n_elements)]
    # The matric-suction apparent cohesion per Gauss point, or None on a model
    # without it. The REPORTED element yield function has to be read on the
    # envelope the trial was solved on, which is c' + c_suction, and the
    # viscoplastic path adds the element mean of exactly this quantity.
    csuc_by_gp = ([[0.0] * len(prep["elem_gp_data"][e]) for e in range(n_elements)]
                  if any(g.get('c_suc') is not None for g in groups) else None)
    for grp in groups:
        _cs = grp.get('c_suc')
        for k, (e, g) in enumerate(grp['pairs']):
            sig_by_gp[e][g] = grp['_sig'][k]
            branch_by_gp[e][g] = int(grp['_branch'][k])
            ep_by_gp[e][g] = grp['ep'][k]
            if csuc_by_gp is not None and _cs is not None:
                csuc_by_gp[e][g] = float(_cs[k])

    final_stresses = np.zeros((n_elements, 4))
    plastic_elements = np.zeros(n_elements, dtype=bool)
    yield_function_out = np.zeros(n_elements)
    vp_shear_strain = np.zeros(n_elements)
    for e in range(n_elements):
        n_gp = len(sig_by_gp[e])
        sig_avg = sum(sig_by_gp[e]) / n_gp
        u_avg = (sum(prep["u_gp"][e]) / len(prep["u_gp"][e])) if prep["u_gp"][e] else 0.0
        stress_total = sig_avg - np.array([u_avg, u_avg, 0.0, u_avg])
        _, sig_vm, _ = stress_invariants(stress_total)
        final_stresses[e] = [-stress_total[0], -stress_total[1], stress_total[2], sig_vm]
        sigm, dsbar, theta = stress_invariants(sig_avg)
        _c_rep = c_reduced[e]
        if csuc_by_gp is not None:
            _c_rep = _c_rep + sum(csuc_by_gp[e]) / n_gp
        yield_function_out[e] = mc_yield_invariants(sigm, dsbar, theta,
                                                    _c_rep, phi_reduced[e])
        plastic_elements[e] = any(b != _NR_ELASTIC for b in branch_by_gp[e])
        ep_avg = sum(ep_by_gp[e]) / n_gp
        vp_shear_strain[e] = float(np.sqrt((ep_avg[0] - ep_avg[1]) ** 2
                                           + ep_avg[2] ** 2))
    if elastic_by_elem is not None:
        plastic_elements[elastic_by_elem] = False

    if _nr_export is not None:
        # Hand the live solve state to the ramp driver, which continues it.
        _nr_export.update({
            'groups': groups, 'bars': bars, 'piles': piles,
            'trans_dofs': trans_dofs,
            'pattern': pattern, 'base_loads': base_loads,
            'oob_fn': _oob, 'h_eps': h_eps, 'u_elastic_scale': u_elastic_scale,
            'free_dofs': free_dofs, 'n_dof': n_dof, 'u': u,
            'disp_limit': nr_disp_limit, 'prep': prep,
            'u_datum': u_datum, 'deep_free': deep_free,
            'deep_dof_mask': deep_dof_mask,
        })

    strains = compute_strains(nodes, elements, element_types, u,
                              dof_offset=dof_offset)
    # Measured from the datum (zero without a carried in-situ state), exactly as the
    # viscoplastic path measures it.
    u_reported = u if _nr_init_state is None else u - u_datum
    max_disp = _nr_umax(u_reported, trans_dofs) if u.size else 0.0
    # The equilibrated state, for a later solve to start from. Same shape and same
    # meaning as the viscoplastic driver's `_k0_state`, over the same Gauss-point
    # groups in the same order, so either driver's state can seed the other.
    nr_k0_state = None
    if groups and groups[0].get('sig0') is not None:
        nr_k0_state = {"u": u.copy(),
                       "evp": [g['ep'].copy() for g in groups],
                       "F": float(F), "converged": bool(converged)}
    u_ratio = (max_disp / u_elastic_scale) if u_elastic_scale > 0 else None
    verdict = 'CONVERGED' if converged else 'FAILED'

    if debug_level >= 1:
        print(f"  Newton-Raphson (F={F:.3f}): {'CONVERGED' if converged else 'FAILED'} "
              f"— {n_steps} load increment(s), {n_cuts} cut(s), "
              f"{total_iterations} Newton iterations, load factor {lam:.4f}, "
              f"out-of-balance {last_oob:.2e}"
              + ('' if converged else f", exit {exit_reason}")
              + (f", max|u| {max_disp:.4g}"
                 + (f" (deep {nr_disp_deep:.4g}, skin {nr_disp_skin:.4g})"
                    if nr_disp_skin is not None else "")
                 + (f" against limit {nr_disp_limit:.4g}"
                    if nr_disp_limit is not None else "")))

    return {
        "converged": bool(converged),
        "stable": bool(converged),
        "verdict": verdict,
        "u_ratio": u_ratio,
        "u_growth": None,
        "u_elastic_scale": u_elastic_scale,
        "exit_reason": exit_reason,
        "_k0_state": nr_k0_state,
        "plateau_iteration": None,
        "plateau_ratio": None,
        "diverging_iteration": (total_iterations if not converged else None),
        # WHY the trial failed, distinguishable at a glance: the load could not be
        # carried at all ('load_step_floor'), it was carried but not to the force
        # tolerance ('force_tolerance'), or it was carried in force but into a state
        # the small-strain model cannot represent ('displacement_limit').
        "diverging_signal": (None if converged else
                             {'diverging': 'load_step_floor',
                              'iteration_cap': 'force_tolerance',
                              'displacement_limit': 'displacement_limit'}
                             .get(exit_reason, 'load_step_floor')),
        # The bound the final state was measured against (None = bound off).
        "nr_disp_limit": nr_disp_limit,
        # The displacement the bound was READ on, and the skin's own. Without a
        # min_slip_depth filter the first is `max_displacement` and the second is
        # None; with one, the deep field is what the verdict rests on and the skin's
        # travel is reported rather than hidden or used.
        "nr_max_disp_deep": nr_disp_deep,
        "nr_max_disp_skin": nr_disp_skin,
        # The same reading for the out-of-balance: the whole model's, and the
        # excluded skin's. Diagnostics, not a verdict.
        "nr_oob_global": nr_oob_global,
        "nr_oob_skin": nr_oob_skin,
        "early_exit_suppressed": False,
        "budget_extensions": 0,
        "iteration_budget": nr_max_iter,
        "failure_criterion": "newton",
        # The largest Mohr-Coulomb violation over all Gauss points at the reported
        # state, as a fraction of the local strength scale (see above). On a
        # converged trial this is the yield half of the verdict's evidence; the
        # force half is `residual`.
        "nr_max_yield_violation": float(max_yield_violation),
        "iterations": int(total_iterations),
        "displacements": u_reported,
        "displacements_elastic": u_elastic,
        "stresses": final_stresses,
        "strains": strains,
        "vp_shear_strain": vp_shear_strain,
        "plastic_elements": plastic_elements,
        "yield_function": yield_function_out,
        "max_displacement": max_disp,
        "plastic_strains": {e: np.array(ep_by_gp[e]) for e in range(n_elements)},
        "algorithm": "Newton-Raphson, consistent Mohr-Coulomb tangent (psi = 0)",
        "F": F,
        # `residual` means the SAME thing on both drivers: the relative change in
        # the displacement field over the last iteration, max|du| / max|u|. It used
        # to carry the out-of-balance here, which reads like the viscoplastic key of
        # the same name but is a different physical quantity by three or more orders
        # of magnitude, so any code comparing the two drivers on it was comparing
        # nothing. The out-of-balance is `unbalanced_force_ratio` — like-for-like on
        # both drivers — and `nr_last_oob` for the Newton-specific reading.
        "residual": float(last_rel_du),
        "unbalanced_force_ratio": last_oob,
        "nr_last_oob": last_oob,
        # Constitutive work actually done: residual evaluations plus every line
        # search backtrack. `iterations` counts only the outer Newton iterations, so
        # it understates the work by up to _NR_LS_MAX per iteration.
        "nr_force_evals": int(n_force_evals),
        "plastic_fraction": (int(np.count_nonzero(plastic_elements)) / n_elements
                             if n_elements else 0.0),
        # Newton-specific accounting, read by the spike's measurement harness.
        # The branch histogram is how a run is read for WHY it behaved as it did:
        # apex points carry a zero tangent, so a field full of them has no
        # stiffness left to iterate on.
        "nr_branch_counts": {
            name: int(sum(int(np.count_nonzero(g['_branch'] == code))
                          for g in groups))
            for code, name in sorted(_NR_BRANCH_NAMES.items())
        },
        # The largest violation of the Rankine cap at the reported state, as a
        # fraction of the same local strength scale, or None on a model with no
        # cap. It is the tensile half of the same evidence `nr_max_yield_violation`
        # carries for the shear surface, and it is folded into that number too, so
        # a reader who checks only one still sees a capped model's whole
        # admissibility.
        "nr_max_tension_violation": max_tension_violation,
        # Whether this trial was solved from a viscoplastic predictor's state
        # rather than from zero, and how many predictor iterations it cost (see
        # _NR_VP_PREDICTOR_ITERS). Zero on every trial the driver solves on its own,
        # which is every trial that does not fail at the load-step floor.
        "nr_predictor_iterations": 0,
        "nr_load_steps": n_steps,
        "nr_step_cuts": n_cuts,
        "nr_step_iterations": step_iters,
        "nr_load_factor": float(lam),
        # Empty arrays on a model with no bars, exactly as the viscoplastic path
        # returns them there (`forces_1d if has_1d_elements else np.array([])`).
        "forces_1d": forces_1d_out,
        "failed_1d_elements": failed_1d_out,
        "softened_1d_elements": softened_1d_out,
        # Empty arrays on a model with no pile, exactly as the viscoplastic path
        # returns them there.
        "forces_pile_axial": pile_axial if n_pile_out else np.array([]),
        "forces_pile_lateral": pile_shear if n_pile_out else np.array([]),
        "forces_pile_moment": pile_moment if n_pile_out else np.zeros((0, 2)),
        "pile_plastic_rotation": pile_prot if n_pile_out else np.zeros((0, 2)),
        "yielded_pile_V": pile_yV if n_pile_out else np.array([], dtype=bool),
        "yielded_pile_M": pile_yM if n_pile_out else np.array([], dtype=bool),
        "yielded_pile": ((pile_yV | pile_yM) if n_pile_out
                         else np.array([], dtype=bool)),
    }


# ===================== The monotonic strength-reduction ramp (SPIKE) =========
#
# An alternative to the bisection, reachable only on the Newton driver:
# solve_ssrm(..., fem_solver='newton', ssrm_driver='ramp'). See SPIKE.md, "RAMP".
#
# The bisection asks an independent question at every trial — it drops the model
# back to zero displacement, re-applies gravity from nothing and rediscovers the
# whole elastoplastic history at the new strength. The ramp carries ONE history:
# it reaches equilibrium once at the starting strength and then reduces strength
# monotonically, each step warm-started from the state before it. The gravity load
# never changes after the first solve; only c and tan(phi) do. So a step is a
# continuation in F rather than in load, and the mechanism develops instead of
# being rebuilt.
#
# Two properties follow by construction rather than by tuning:
#   * the verdict sequence cannot be non-monotone, because there is one history;
#   * no strength more than one step past the highest one carried is ever
#     evaluated, which is where the Newton driver is at its worst (17x to 47x the
#     viscoplastic driver's constitutive work on trials well past failure).
_RAMP_DF_INIT = 0.05       # first strength increment
_RAMP_DF_MIN = 0.005       # below this a failed step means the limit is reached
_RAMP_DF_GROW = 1.6        # increment growth after a comfortable step
# The ramp reports the MIDPOINT of its final interval, (F_stands + F_refused)/2,
# because that is the convention every locked and published factor of safety in
# this repository is defined on: the bisection returns the midpoint of the last
# bracket it could not split further. The ramp used to report the last strength it
# had CARRIED, floored to 0.01, which is the interval's lower edge rounded down —
# two separate downward biases, measured at 0.0031 to 0.0119 below the bisection's
# answer across the eight spike benchmarks, none of it a difference in what the two
# routes found. The raw facts are still reported, unrounded, as `ramp_last_carried`
# and `ramp_first_refused`, so the interval is always readable.


def _ssrm_ramp_newton(fem_data, F_min, F_max, *, prep, force_tol, convergence_tol,
                      max_iterations, oob_window, dt_scale,
                      tension_cutoff, min_slip_depth, ssr_exclude_mask,
                      tension_cap_by_elem, tension_srf, elastic_mask,
                      suction_phi_b, suction_cap, k0, f_adjust, f_min_floor,
                      max_expand, max_iterations_ceiling=50000, early_failure=True,
                      init_state=None, debug_level=0, progress_callback=None,
                      cancel_check=None):
    """Walk F up from a converged state, warm-starting every step.

    Returns the same result shape solve_ssrm's bisection drivers return, so the
    caller consumes it unchanged. ``trials`` records every strength actually
    evaluated, in the order evaluated, including the steps that were rejected and
    retried at half the increment.
    """
    trials = []
    ctx = {}

    def _cold(F, export_only=False):
        # `_nr_export` is what makes a solve a RAMP FOOT: it hands back the working
        # state the ramp then continues, and it is also what holds the viscoplastic
        # predictor off, because a ramp manages its own plastic history. The final
        # export solve continues nothing, so it takes the predictor like any other
        # standalone trial — without it the state written to figures and CSVs would
        # be whatever the cold path could reach, which on a cohesionless model is
        # not the state at the limit the ramp just carried.
        return solve_fem(fem_data, F=F, debug_level=max(0, debug_level - 1),
                         force_tol=force_tol, oob_window=oob_window,
                         dt_scale=dt_scale,
                         max_iterations=max_iterations, tolerance=convergence_tol,
                         max_disp_factor=None,
                         tension_cutoff=tension_cutoff,
                         min_slip_depth=min_slip_depth,
                         ssr_exclude_mask=ssr_exclude_mask,
                         tension_cap_by_elem=tension_cap_by_elem,
                         tension_srf=tension_srf, elastic_mask=elastic_mask,
                         suction_phi_b=suction_phi_b, suction_cap=suction_cap,
                         k0=k0, fem_solver='newton', _prepared=prep,
                         max_iterations_ceiling=max_iterations_ceiling,
                         early_failure=early_failure, _init_state=init_state,
                         _nr_export=(None if export_only else ctx))

    def _record(F, verdict, iters, fevals, oob, maxu, kind, why=None):
        trials.append({"F": float(F), "verdict": verdict, "iterations": int(iters),
                       "nr_force_evals": int(fevals),
                       "unbalanced_force_ratio": float(oob),
                       "max_displacement": float(maxu), "ramp_step": kind,
                       "exit_reason": why})

    # --- the foot of the ramp: one cold solve, walked down if it does not stand --
    F0 = float(F_min)
    n_expand = 0
    sol0 = _cold(F0)
    _record(F0, sol0['verdict'], sol0['iterations'], sol0.get('nr_force_evals', 0),
            sol0['unbalanced_force_ratio'], sol0['max_displacement'], 'cold',
            sol0.get('exit_reason'))
    while not sol0['converged']:
        if F0 <= f_min_floor + 1e-9 or n_expand >= max_expand:
            msg = (f"SSRM (ramp): the slope does not reach equilibrium even at "
                   f"F = {F0:.2f} — it is unstable at or below this strength-"
                   f"reduction factor (FS < {F0:.2f}).")
            print(f"\n{msg}")
            return {"converged": False, "error": msg, "FS": None, "trials": trials}
        F0 = max(f_min_floor, F0 - f_adjust)
        n_expand += 1
        ctx = {}
        sol0 = _cold(F0)
        _record(F0, sol0['verdict'], sol0['iterations'],
                sol0.get('nr_force_evals', 0), sol0['unbalanced_force_ratio'],
                sol0['max_displacement'], 'cold', sol0.get('exit_reason'))

    groups, pattern, bars = ctx['groups'], ctx['pattern'], ctx.get('bars')
    piles, trans_dofs = ctx.get('piles'), ctx.get('trans_dofs')
    base_loads, oob_fn = ctx['base_loads'], ctx['oob_fn']
    free_dofs, n_dof = ctx['free_dofs'], ctx['n_dof']
    h_eps, u_el = ctx['h_eps'], ctx['u_elastic_scale']
    disp_limit = ctx['disp_limit']
    deep_free = ctx.get('deep_free')
    deep_dof_mask = ctx.get('deep_dof_mask')
    restrength = ctx['restrength']
    u = ctx['u']
    # Displacement is measured from the in-situ state on a K0 run and from zero
    # otherwise, the same datum the foot's own verdict was read on.
    u_datum = ctx.get('u_datum')
    if u_datum is None:
        u_datum = np.zeros(n_dof)

    F_stands = F0
    last_solution = sol0
    dF = _RAMP_DF_INIT
    n_steps, n_retries = 0, 0
    # Every solve so far — the foot AND any failed walk-down solves — is already
    # in `trials`; the totals count them all, not just the solve that stood.
    total_iters = int(sum(t['iterations'] for t in trials))
    total_fevals = int(sum(t['nr_force_evals'] for t in trials))
    warm_iters = []
    F_refused = None
    cancelled = False
    pred_iters = 0             # viscoplastic predictor iterations, charged on top

    if debug_level >= 1:
        print("=== SSRM Analysis (Newton, monotonic ramp) ===")
        print(f"  Foot of the ramp: F = {F0:.4f} stands "
              f"({sol0['iterations']} iterations, cold)")

    while True:
        if cancel_check is not None and cancel_check():
            cancelled = True
            break
        F_try = F_stands + dF
        if F_try > F_max + 1e-12:
            # The ramp has carried the whole requested range. There is no failed
            # step above it, so the answer is a lower bound, not a limit.
            msg = (f"SSRM (ramp): the slope still stands at F = {F_stands:.4f}, the "
                   f"top of the requested range — FS > {F_stands:.2f}. Raise F_max.")
            print(f"\n{msg}")
            return {"converged": False, "error": msg, "FS": None, "trials": trials,
                    "last_solution": last_solution}

        restrength(groups, F_try)
        seeded = False
        # The step's own history, kept so a REJECTED step leaves nothing behind —
        # the plastic strains a predictor seed or a softening sub-step committed,
        # and the post-peak set the latch may have grown. Softening only ever grows
        # WITHIN a solve, and the ramp is one solve, so a bar that ruptured at F_k is
        # ruptured at F_{k+1}; but a step that is refused never happened.
        ep_save = [grp['ep'] for grp in groups]
        soft_save = ([(bg['softened'].copy(), bg['t_cap'].copy()) for bg in bars]
                     if bars is not None else None)
        ok, u_try, it, fe, oob, _rel = _nr_equilibrate(
            groups, pattern, u, base_loads, free_dofs, n_dof, h_eps, force_tol,
            oob_fn, _NR_MAX_ITER, u_el, debug_level=debug_level,
            label=f"ramp F={F_try:.4f}", bars=bars, piles=piles,
            trans_dofs=trans_dofs, deep_free=deep_free)
        total_iters += it
        total_fevals += fe
        if not (ok and oob < force_tol):
            # The step is refused, and on a cohesionless model that is the same
            # failure the bisection retries (see _NR_VP_PREDICTOR_ITERS): the warm
            # plastic history this ramp carries is not the field this strength
            # needs, and the corrector cannot grow one from it. So grow one the same
            # way — a viscoplastic predictor at THIS strength, under the same
            # Rankine cap — and give the corrector one attempt from that state.
            #
            # The verdict stays the corrector's: the same force gate and the same
            # displacement bound the cold step just failed, read on the corrected
            # state and not on the seed. The ramp used to walk every step cold,
            # which is why it read 1.0531 on the three-layer model against the
            # bisection's 1.2109.
            for _chunk, _ceiling, _seed_prep, _seed_ef in _nr_predictor_rungs(
                    prep, len(fem_data['elements']), max_iterations,
                    max_iterations_ceiling, early_failure):
                _vp = solve_fem(
                    fem_data, F=F_try, debug_level=0, max_iterations=_chunk,
                    max_iterations_ceiling=_ceiling, tolerance=convergence_tol,
                    max_disp_factor=None, tension_cutoff=tension_cutoff,
                    dt_scale=dt_scale, force_tol=force_tol,
                    oob_window=oob_window, min_slip_depth=min_slip_depth,
                    ssr_exclude_mask=ssr_exclude_mask,
                    tension_cap_by_elem=tension_cap_by_elem,
                    tension_srf=tension_srf, elastic_mask=elastic_mask,
                    suction_phi_b=suction_phi_b, suction_cap=suction_cap,
                    k0=k0, _init_state=init_state,
                    _prepared=_seed_prep, early_failure=_seed_ef,
                    fem_solver='viscoplastic')
                pred_iters += int(_vp.get('iterations', 0) or 0)
                if k0 is not None:
                    # The predictor's reported displacement is measured from the
                    # in-situ datum; its `_k0_state` carries the absolute field and
                    # the plastic strain in this model's own group order.
                    _ks = _vp['_k0_state']
                    _u_seed = np.asarray(_ks['u'], dtype=float)[:n_dof].copy()
                    for grp, _ev in zip(groups, _ks['evp']):
                        grp['ep'] = np.array(_ev, dtype=float, copy=True)
                else:
                    _u_seed = np.asarray(_vp['displacements'],
                                         dtype=float)[:n_dof].copy()
                    _ep_vp = _vp['plastic_strains']
                    for grp in groups:
                        _ep = np.empty_like(grp['ep'])
                        for _k, (_e, _gp) in enumerate(grp['pairs']):
                            _ep[_k] = _ep_vp[_e][_gp]
                        grp['ep'] = _ep
                ok, u_try, it2, fe2, oob, _rel = _nr_equilibrate(
                    groups, pattern, _u_seed, base_loads, free_dofs, n_dof, h_eps,
                    force_tol, oob_fn, _NR_MAX_ITER, u_el, debug_level=debug_level,
                    label=f"ramp F={F_try:.4f} (seeded)", bars=bars, piles=piles,
                    trans_dofs=trans_dofs, deep_free=deep_free)
                # Work is cumulative: the refused cold step, every predictor run and
                # every corrector are all charged to this step.
                it += it2
                fe += fe2
                total_iters += it2
                total_fevals += fe2
                if ok:
                    seeded = True
                    break
                # Nothing was carried, so the ramp's own history goes back exactly
                # as it was before the seed overwrote it.
                for grp, _ep0 in zip(groups, ep_save):
                    grp['ep'] = _ep0
        # The post-peak latch, on the same converged full-load state the bisection
        # reads it on. The set is carried FORWARD across steps by construction — the
        # bar groups are built once for the whole ramp — so a bar that ruptured at a
        # lower strength is ruptured here, which is what one continuous history
        # means. A drop that reaches the step floor refuses the step.
        if ok and oob < force_tol and bars is not None and any(
                bg['can_soften'].any() for bg in bars):
            ok, u_try, it_s, fe_s, oob_s, _r, _c = _nr_soften_latch(
                bars, groups, pattern, u_try, base_loads, free_dofs, n_dof, h_eps,
                force_tol, oob_fn, _NR_MAX_ITER, u_el, piles=piles,
                trans_dofs=trans_dofs, debug_level=debug_level,
                deep_free=deep_free)
            it += it_s
            fe += fe_s
            total_iters += it_s
            total_fevals += fe_s
            if _r:
                oob = oob_s
        # The bound reads the DEEP field when a min_slip_depth filter is on — the
        # same standard the bisection's trials are held to, for the same reason.
        maxu, maxu_skin = (_nr_umax_split(u_try - u_datum, trans_dofs,
                                          deep_dof_mask)
                           if u_try.size else (0.0, None))
        # The SAME admissibility standard a converged bisection trial passes: force
        # equilibrium under the tolerance, and a state the small-strain model can
        # represent.
        admissible = bool(ok) and oob < force_tol and (
            disp_limit is None or maxu <= disp_limit)
        why = ('converged' if admissible else
               'displacement_limit' if (ok and oob < force_tol) else
               'diverging')
        _record(F_try, 'CONVERGED' if admissible else 'FAILED', it, fe, oob, maxu,
                f'step dF={dF:.5f}', why)

        if admissible:
            u = u_try
            for grp in groups:
                grp['_u'] = u
            _nr_commit_plastic_strain(groups)
            F_stands = F_try
            n_steps += 1
            warm_iters.append(it)
            if it <= _NR_COMFORT:
                dF = min(_RAMP_DF_INIT, dF * _RAMP_DF_GROW)
            if debug_level >= 1:
                print(f"    F = {F_stands:.4f} carried in {it} iterations "
                      f"(oob {oob:.2e}, max|u| {maxu:.4g})")
            if progress_callback is not None:
                try:
                    progress_callback(
                        min(1.0, (F_stands - F0) / max(1e-9, F_max - F0)),
                        f"ramp F={F_stands:.3f}")
                except Exception:
                    pass
        else:
            # Reject the step and retry from the SAME converged state at half the
            # increment. Nothing was committed, so there is nothing to undo but the
            # strengths, which the next attempt overwrites — and a predictor seed or
            # a softening sub-step, if one committed a plastic strain or dropped a
            # bar here only for the step to be refused.
            for grp, _ep0 in zip(groups, ep_save):
                grp['ep'] = _ep0
            if soft_save is not None:
                for bg, (_s0, _c0) in zip(bars, soft_save):
                    bg['softened'], bg['t_cap'] = _s0, _c0
            F_refused = F_try
            n_retries += 1
            if debug_level >= 1:
                print(f"    F = {F_try:.4f} refused ({why}, {it} iterations, "
                      f"oob {oob:.2e}, max|u| {maxu:.4g}); halving dF to "
                      f"{dF / 2:.5f}")
            dF *= 0.5
            if dF < _RAMP_DF_MIN:
                break

    if cancelled:
        msg = (f"SSRM (ramp): cancelled at F = {F_stands:.4f} after {n_steps} "
               f"steps; no limit was reached.")
        return {"converged": False, "error": msg, "FS": None, "trials": trials,
                "last_solution": last_solution}

    # Export the state at the limit, not the foot: figures and CSVs are written
    # from `last_solution`, and the ramp's warm state carries no post-processing.
    # One cold solve at the last carried strength matches what the bisection
    # driver exports for its final converged trial.
    sol_end = _cold(F_stands, export_only=True)
    _record(F_stands, sol_end['verdict'], sol_end['iterations'],
            sol_end.get('nr_force_evals', 0), sol_end['unbalanced_force_ratio'],
            sol_end['max_displacement'], 'final_export', sol_end.get('exit_reason'))
    total_iters += int(sol_end['iterations'])
    total_fevals += int(sol_end.get('nr_force_evals', 0))
    pred_iters += int(sol_end.get('nr_predictor_iterations', 0) or 0)
    if sol_end.get('converged'):
        last_solution = sol_end

    # The limit lies between the last strength carried and the first one refused,
    # and is reported as that interval's MIDPOINT — the bisection's convention, and
    # therefore the convention of every locked and published factor of safety here.
    # See the note at _RAMP_DF_GROW.
    FS = 0.5 * (F_stands + F_refused)
    restrength(groups, F_stands)
    if debug_level >= 1:
        print(f"  Ramp limit: {F_stands:.4f} carried, {F_refused:.4f} refused "
              f"-> FS = {FS:.4f} ({n_steps} steps, {n_retries} retries, "
              f"{total_iters} iterations, {total_fevals} force evaluations)")
    return {
        "converged": True,
        "FS": float(FS),
        "last_solution": last_solution,
        "iterations_ssrm": n_steps,
        "final_interval": (float(F_stands), float(F_refused)),
        "interval_width": float(F_refused - F_stands),
        "failed_edge_softened": None,
        "trials": trials,
        "inconclusive": [],
        "note": None,
        "failure_criterion": "newton_ramp",
        "method": "SSRM — monotonic strength-reduction ramp (Newton)",
        # Ramp accounting, read by the spike's measurement harness.
        "ramp_steps": n_steps,
        "ramp_retries": n_retries,
        "ramp_iterations": total_iters,
        "ramp_force_evals": total_fevals,
        # Viscoplastic predictor iterations spent rescuing refused steps, charged
        # on top of the Newton work above (see _NR_VP_PREDICTOR_ITERS).
        "ramp_predictor_iterations": pred_iters,
        "ramp_warm_iterations": warm_iters,
        "ramp_foot": float(F0),
        "ramp_last_carried": float(F_stands),
        "ramp_first_refused": float(F_refused),
    }


def print_reinforcement_summary(fem_data, solution):
    """
    Print a summary table of reinforcement line results.

    Groups 1D elements by reinforcement line and reports per-line statistics
    including element counts, force ranges, and failure modes. Each line's
    verdict is :func:`xslope.fem_details.reinforcement_status`, which is where
    the Studio 1D Details panel and the report read it from as well: the word
    printed here is the word on the screen.
    """
    from .fem_details import (REINFORCEMENT_STATE_ORDER, REINFORCEMENT_STATES,
                              reinforcement_status)

    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    n_1d = len(elements_1d)
    if n_1d == 0:
        return

    element_materials_1d = fem_data["element_materials_1d"]
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(n_1d, dtype=bool))
    t_allow_by_elem = fem_data["t_allow_by_1d_elem"]
    t_res_by_elem = fem_data["t_res_by_1d_elem"]
    forces = solution.get("forces_1d", np.zeros(n_1d))
    failed = solution.get("failed_1d_elements", np.zeros(n_1d, dtype=bool))
    softened = solution.get("softened_1d_elements", np.zeros(n_1d, dtype=bool))
    if len(softened) != n_1d:
        softened = np.zeros(n_1d, dtype=bool)

    # Filter out pile elements — reinforcement reported separately
    reinf_mask = ~pile_elem_mask
    if not np.any(reinf_mask):
        return

    # Get per-line Tmax and Tres from slope_data stored in fem_data
    # element_materials_1d is 1-based line ID
    line_ids = np.unique(element_materials_1d[reinf_mask])

    print("\n=== Reinforcement Summary ===")
    print(f"{'Line':>4}  {'Elems':>5}  {'Max T':>8}  {'Avg T':>8}  "
          f"{'Tension':>7}  {'In Lp':>5}  {'Yielded':>7}  {'Pullout':>7}  {'Status'}")
    print("-" * 80)

    statuses_seen = set()
    for line_id in sorted(line_ids):
        mask = element_materials_1d == line_id
        n_elem = int(mask.sum())
        line_forces = forces[mask]
        line_failed = failed[mask]
        line_softened = softened[mask]
        line_t_allow = t_allow_by_elem[mask]
        line_t_res = t_res_by_elem[mask]

        # Max Tallow for this line (the full-capacity elements)
        t_max_line = line_t_allow.max() if n_elem > 0 else 0.0

        # Active: elements carrying tension
        n_active = int((line_forces > 0).sum())

        # Pullout zone: elements where T_allow < T_max (reduced by proximity to end)
        n_pullout = int(((line_t_allow < t_max_line - 1e-6) & (line_t_allow > 1e-6)).sum())

        # Yielded: elements that reached their allowable capacity and are now
        # holding it (elastic-perfectly-plastic bar — see solve_fem; there is no
        # rupture, so a yielded element still carries T_allow).
        yielded_mask = line_failed

        # Elements inside a pullout ramp carry less than the line's full Tmax
        in_lp_mask = (line_t_allow < t_max_line - 1e-6) & (line_t_allow > 1e-6)
        at_end_mask = line_t_allow < 1e-6   # the very ends develop no tension
        lp_zone_mask = in_lp_mask | at_end_mask

        # Yielded at the pullout (embedment-limited) capacity vs. at full Tmax
        n_yield_in_lp = int((yielded_mask & lp_zone_mask).sum())
        n_yield_outside_lp = int((yielded_mask & ~lp_zone_mask).sum())

        max_t = line_forces.max() if n_elem > 0 else 0.0
        active_forces = line_forces[line_forces > 0]
        avg_t = active_forces.mean() if len(active_forces) > 0 else 0.0

        # The line's verdict, from the one function that decides it — the same
        # call, and so the same word, the Studio 1D Details panel puts on the
        # line and the report writes into its prose.
        status_key, phrase = reinforcement_status(
            line_forces, line_t_allow, t_res=line_t_res, failed=line_failed,
            softened=line_softened)
        status = phrase.upper()

        statuses_seen.add(status_key)

        print(f"{line_id:>4}  {n_elem:>5}  {max_t:>8.1f}  {avg_t:>8.1f}  "
              f"{n_active:>7}  {n_pullout:>5}  {n_yield_outside_lp:>7}  "
              f"{n_yield_in_lp:>7}  {status}")

    print("-" * 80)

    # What each state that appeared means, in the words every other surface
    # defines it in (:data:`xslope.fem_details.REINFORCEMENT_STATES`).
    notes = [f"{REINFORCEMENT_STATES[key][0].upper()}: the line "
             f"{REINFORCEMENT_STATES[key][1]}."
             for key in REINFORCEMENT_STATE_ORDER if key in statuses_seen]
    if notes:
        print()
        for note in notes:
            print(f"  {note}")


def print_pile_summary(fem_data, solution):
    """
    Print a summary table of pile results, grouped by pile line.
    """
    n_pile = fem_data.get("n_pile_elements", 0)
    if n_pile == 0:
        return

    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    nodes = fem_data["nodes"]
    element_materials_1d = fem_data.get("element_materials_1d", np.array([]))
    pile_elem_indices = fem_data.get("pile_elem_indices", np.array([], dtype=int))
    forces_axial = solution.get("forces_pile_axial", np.zeros(n_pile))
    forces_shear = solution.get("forces_pile_lateral", np.zeros(n_pile))
    forces_moment = solution.get("forces_pile_moment", np.zeros((n_pile, 2)))
    yielded = solution.get("yielded_pile", np.zeros(n_pile, dtype=bool))
    yielded_V = solution.get("yielded_pile_V", np.zeros(n_pile, dtype=bool))
    yielded_M = solution.get("yielded_pile_M", np.zeros(n_pile, dtype=bool))
    V_cap_arr = fem_data.get("V_cap_by_pile_elem", np.full(n_pile, float('inf')))
    M_cap_arr = fem_data.get("M_cap_by_pile_elem", np.full(n_pile, float('inf')))

    # Check if any pile has capacity limits
    has_V_cap = np.any(V_cap_arr < float('inf'))
    has_M_cap = np.any(M_cap_arr < float('inf'))
    has_capacity = has_V_cap or has_M_cap

    # Group pile elements by pile line (material ID)
    pile_line_ids = {}
    for p_idx in range(n_pile):
        global_idx = pile_elem_indices[p_idx]
        line_id = element_materials_1d[global_idx]
        if line_id not in pile_line_ids:
            pile_line_ids[line_id] = []
        pile_line_ids[line_id].append(p_idx)

    L_elems = fem_data.get("elem_length_by_pile_elem", np.zeros(n_pile))
    S_arr = fem_data.get("S_by_pile_elem", np.ones(n_pile))

    print(f"\n=== Pile Summary ===")
    header = (f"{'Pile':>4}  {'Elems':>5}  {'Max |T|':>8}  {'Max |V|':>8}  "
              f"{'Max |M|':>8}")
    if has_V_cap:
        header += f"  {'V_cap':>8}"
    if has_M_cap:
        header += f"  {'M_cap':>8}"
    header += f"  {'Yielded':>7}  {'Status'}"
    print(header)
    print("-" * (len(header) + 2))

    statuses_seen = set()
    pile_num = 0
    for line_id in sorted(pile_line_ids.keys()):
        pile_num += 1
        indices = pile_line_ids[line_id]
        n_elem = len(indices)
        n_yielded = np.sum(yielded[indices])

        max_axial = np.max(np.abs(forces_axial[indices]))
        max_shear = np.max(np.abs(forces_shear[indices]))
        max_moment = np.max(np.abs(forces_moment[indices]))

        v_cap = V_cap_arr[indices[0]]
        m_cap = M_cap_arr[indices[0]]

        if n_yielded > 0:
            status = "YIELDED"
        elif (v_cap < float('inf') and max_shear > 0.95 * v_cap) or \
             (m_cap < float('inf') and max_moment > 0.95 * m_cap):
            status = "NEAR CAP"
        else:
            status = "OK"
        statuses_seen.add(status)

        row = (f"{pile_num:>4}  {n_elem:>5}  {max_axial:>8.1f}  {max_shear:>8.1f}  "
               f"{max_moment:>8.1f}")
        if has_V_cap:
            vcap_str = f"{v_cap:.1f}" if v_cap < float('inf') else "-"
            row += f"  {vcap_str:>8}"
        if has_M_cap:
            mcap_str = f"{m_cap:.1f}" if m_cap < float('inf') else "-"
            row += f"  {mcap_str:>8}"
        row += f"  {n_yielded:>3}/{n_elem}  {status}"
        print(row)

    print("-" * (len(header) + 2))

    # Capacity calculation notes
    if has_capacity:
        print()
        pile_num = 0
        for line_id in sorted(pile_line_ids.keys()):
            pile_num += 1
            indices = pile_line_ids[line_id]
            v_cap = V_cap_arr[indices[0]]
            m_cap = M_cap_arr[indices[0]]
            S_pile = S_arr[indices[0]]

            if v_cap < float('inf') or m_cap < float('inf'):
                print(f"  Pile {pile_num} capacity (per unit width = per pile / S):")
                if v_cap < float('inf'):
                    print(f"    V_cap/S = {v_cap:.1f}")
                if m_cap < float('inf'):
                    print(f"    M_cap/S = {m_cap:.1f}")
                n_yV = int(np.sum(yielded_V[indices]))
                n_yM = int(np.sum(yielded_M[indices]))
                if n_yV > 0:
                    print(f"    {n_yV} element(s) reached V_cap (shear hinge)")
                if n_yM > 0:
                    print(f"    {n_yM} element(s) reached M_cap (plastic hinge)")

    # Status notes
    status_notes = {
        "OK": "OK: All elements within structural capacity.",
        "NEAR CAP": "NEAR CAP: Max force/moment exceeds 95% of structural capacity.",
        "YIELDED": "YIELDED: Elements reached structural capacity; forces capped at limit (elastic-perfectly-plastic).",
    }
    notes = [status_notes[s] for s in ["OK", "NEAR CAP", "YIELDED"] if s in statuses_seen]
    if notes:
        print()
        for note in notes:
            print(f"  {note}")


def print_detailed_element_summary(fem_data, solution):
    """
    Print a detailed per-element summary table for reinforcement and pile elements.

    For reinforcement: lists each element with its line, centroid coordinates,
    force, allowable force, and status.
    For piles: lists each element with centroid coordinates, axial force,
    lateral force.
    """
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    n_1d = len(elements_1d)
    if n_1d == 0 and fem_data.get("n_pile_elements", 0) == 0:
        return

    nodes = fem_data["nodes"]
    element_materials_1d = fem_data.get("element_materials_1d", np.array([]))
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(n_1d, dtype=bool))
    t_allow_by_elem = fem_data.get("t_allow_by_1d_elem", np.zeros(n_1d))
    t_res_by_elem = fem_data.get("t_res_by_1d_elem", np.zeros(n_1d))
    forces = solution.get("forces_1d", np.zeros(n_1d))
    failed = solution.get("failed_1d_elements", np.zeros(n_1d, dtype=bool))

    # --- Reinforcement element table ---
    reinf_indices = [i for i in range(n_1d) if not pile_elem_mask[i]]
    if reinf_indices:
        print("\n=== Detailed Reinforcement Element Summary ===")
        print(f"{'Elem':>4}  {'Line':>4}  {'X1':>8}  {'Y1':>8}  {'X2':>8}  {'Y2':>8}  "
              f"{'Force':>8}  {'T_allow':>8}  {'T_res':>8}  {'Status'}")
        print("-" * 96)

        for i in reinf_indices:
            elem = elements_1d[i]
            n0 = nodes[elem[0]]
            n1 = nodes[elem[1]]
            line_id = element_materials_1d[i]
            force = forces[i]
            t_allow = t_allow_by_elem[i]
            t_res = t_res_by_elem[i]
            is_failed = failed[i]

            if is_failed:
                # The bar is at its allowable capacity and holding it. Whether
                # that capacity is the full tensile strength or an
                # embedment-limited (pullout) value depends on where the element
                # sits along the line.
                line_mask = element_materials_1d == line_id
                t_max_line = t_allow_by_elem[line_mask].max()
                status = "PULLOUT" if t_allow < t_max_line - 1e-6 else "YIELDED"
            elif force < -1e-6:
                status = "COMPRESS"
            elif force < 1e-6:
                status = "INACTIVE"
            elif force > 0.95 * t_allow and t_allow > 0:
                status = "NEAR CAP"
            else:
                status = "OK"

            print(f"{i:>4}  {line_id:>4}  {n0[0]:>8.2f}  {n0[1]:>8.2f}  {n1[0]:>8.2f}  {n1[1]:>8.2f}  "
                  f"{force:>8.1f}  {t_allow:>8.1f}  {t_res:>8.1f}  {status}")

        print("-" * 96)

    # --- Pile element table ---
    n_pile = fem_data.get("n_pile_elements", 0)
    if n_pile > 0:
        pile_elem_indices = fem_data.get("pile_elem_indices", np.array([], dtype=int))
        forces_axial = solution.get("forces_pile_axial", np.zeros(n_pile))
        forces_lateral = solution.get("forces_pile_lateral", np.zeros(n_pile))
        forces_moment = solution.get("forces_pile_moment", np.zeros((n_pile, 2)))
        yielded = solution.get("yielded_pile", np.zeros(n_pile, dtype=bool))
        yielded_V = solution.get("yielded_pile_V", np.zeros(n_pile, dtype=bool))
        yielded_M = solution.get("yielded_pile_M", np.zeros(n_pile, dtype=bool))
        V_cap_arr = fem_data.get("V_cap_by_pile_elem", np.full(n_pile, float('inf')))
        M_cap_arr = fem_data.get("M_cap_by_pile_elem", np.full(n_pile, float('inf')))
        L_elems = fem_data.get("elem_length_by_pile_elem", np.zeros(n_pile))
        has_V_cap = np.any(V_cap_arr < float('inf'))
        has_M_cap = np.any(M_cap_arr < float('inf'))
        has_cap = has_V_cap or has_M_cap

        print("\n=== Detailed Pile Element Summary ===")
        header = (f"{'Elem':>4}  {'X1':>8}  {'Y1':>8}  {'X2':>8}  {'Y2':>8}  "
                  f"{'Axial':>10}  {'Shear':>10}  {'M1':>10}  {'M2':>10}")
        if has_V_cap:
            header += f"  {'V_cap':>8}"
        if has_M_cap:
            header += f"  {'M_cap':>8}"
        if has_cap:
            header += f"  {'Status'}"
        print(header)
        sep_len = len(header) + 2
        print("-" * sep_len)

        for p_idx in range(n_pile):
            global_idx = pile_elem_indices[p_idx]
            elem = elements_1d[global_idx]
            n0 = nodes[elem[0]]
            n1 = nodes[elem[1]]
            T = forces_axial[p_idx]
            V = forces_lateral[p_idx]
            M1 = forces_moment[p_idx, 0]
            M2 = forces_moment[p_idx, 1]

            row = (f"{p_idx:>4}  {n0[0]:>8.2f}  {n0[1]:>8.2f}  {n1[0]:>8.2f}  {n1[1]:>8.2f}  "
                   f"{T:>10.1f}  {V:>10.1f}  {M1:>10.1f}  {M2:>10.1f}")

            if has_V_cap:
                vcap_str = f"{V_cap_arr[p_idx]:.1f}" if V_cap_arr[p_idx] < float('inf') else "-"
                row += f"  {vcap_str:>8}"
            if has_M_cap:
                mcap_str = f"{M_cap_arr[p_idx]:.1f}" if M_cap_arr[p_idx] < float('inf') else "-"
                row += f"  {mcap_str:>8}"

            if has_cap:
                if yielded[p_idx]:
                    parts = []
                    if yielded_V[p_idx]:
                        parts.append("V")
                    if yielded_M[p_idx]:
                        parts.append("M")
                    status = "YIELDED(" + "+".join(parts) + ")"
                elif (has_V_cap and V_cap_arr[p_idx] < float('inf') and abs(V) > 0.95 * V_cap_arr[p_idx]) or \
                     (has_M_cap and M_cap_arr[p_idx] < float('inf') and max(abs(M1), abs(M2)) > 0.95 * M_cap_arr[p_idx]):
                    status = "NEAR CAP"
                else:
                    status = "OK"
                row += f"  {status}"

            print(row)

        print("-" * sep_len)
        max_M = np.max(np.abs(forces_moment)) if n_pile > 0 else 0.0
        print(f"{'':>4}  {'':>8}  {'':>8}  {'':>8}  {'max abs:':>8}  "
              f"{np.max(np.abs(forces_axial)):>10.1f}  {np.max(np.abs(forces_lateral)):>10.1f}  "
              f"{max_M:>10.1f}")


def _compute_B_and_detJ_quad8(coords, xi, eta):
    """Compute B matrix and det(J) for 8-node quad at (xi, eta)."""
    dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)

    J = np.zeros((2, 2))
    for a in range(8):
        J[0, 0] += dN_dxi[a] * coords[a, 0]
        J[0, 1] += dN_dxi[a] * coords[a, 1]
        J[1, 0] += dN_deta[a] * coords[a, 0]
        J[1, 1] += dN_deta[a] * coords[a, 1]

    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]

    if abs(det_J) < 1e-14:
        return np.zeros((3, 16)), 0.0

    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J

    dN_dx = np.zeros(8)
    dN_dy = np.zeros(8)
    for a in range(8):
        dN_dx[a] = J_inv[0, 0] * dN_dxi[a] + J_inv[0, 1] * dN_deta[a]
        dN_dy[a] = J_inv[1, 0] * dN_dxi[a] + J_inv[1, 1] * dN_deta[a]

    B = np.zeros((3, 16))
    for a in range(8):
        B[0, 2*a] = dN_dx[a]
        B[1, 2*a+1] = dN_dy[a]
        B[2, 2*a] = dN_dy[a]
        B[2, 2*a+1] = dN_dx[a]

    return B, det_J


def _compute_B_and_detJ_tri6(coords, L1, L2, L3):
    """Compute B matrix (3,12) and det(J) for 6-node triangle at area coordinates."""
    dN_dL1, dN_dL2, dN_dL3 = compute_tri6_shape_derivatives(L1, L2, L3)

    x0, y0 = coords[0]
    x1, y1 = coords[1]
    x2, y2 = coords[2]

    # Jacobian from area coordinates to physical coordinates
    det_J = (x0 - x2) * (y1 - y2) - (x1 - x2) * (y0 - y2)
    area = 0.5 * abs(det_J)

    if area < 1e-14:
        return np.zeros((3, 12)), 0.0

    # Area coordinate derivatives w.r.t. physical coordinates
    dL1_dx = (y1 - y2) / (2 * area)
    dL1_dy = (x2 - x1) / (2 * area)
    dL2_dx = (y2 - y0) / (2 * area)
    dL2_dy = (x0 - x2) / (2 * area)
    dL3_dx = (y0 - y1) / (2 * area)
    dL3_dy = (x1 - x0) / (2 * area)

    # Shape function derivatives in physical coordinates via chain rule
    dN_dx = dN_dL1 * dL1_dx + dN_dL2 * dL2_dx + dN_dL3 * dL3_dx
    dN_dy = dN_dL1 * dL1_dy + dN_dL2 * dL2_dy + dN_dL3 * dL3_dy

    B = np.zeros((3, 12))
    for a in range(6):
        B[0, 2*a] = dN_dx[a]
        B[1, 2*a+1] = dN_dy[a]
        B[2, 2*a] = dN_dy[a]
        B[2, 2*a+1] = dN_dx[a]

    return B, det_J


def _compute_B_and_detJ_quad4(coords, xi, eta):
    """Compute B matrix (3,8) and det(J) for 4-node quad at (xi, eta)."""
    dN_dxi, dN_deta = compute_quad4_shape_derivatives(xi, eta)

    J = np.zeros((2, 2))
    for a in range(4):
        J[0, 0] += dN_dxi[a] * coords[a, 0]
        J[0, 1] += dN_dxi[a] * coords[a, 1]
        J[1, 0] += dN_deta[a] * coords[a, 0]
        J[1, 1] += dN_deta[a] * coords[a, 1]

    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]

    if abs(det_J) < 1e-14:
        return np.zeros((3, 8)), 0.0

    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J

    dN_dx = np.zeros(4)
    dN_dy = np.zeros(4)
    for a in range(4):
        dN_dx[a] = J_inv[0, 0] * dN_dxi[a] + J_inv[0, 1] * dN_deta[a]
        dN_dy[a] = J_inv[1, 0] * dN_dxi[a] + J_inv[1, 1] * dN_deta[a]

    B = np.zeros((3, 8))
    for a in range(4):
        B[0, 2*a] = dN_dx[a]
        B[1, 2*a+1] = dN_dy[a]
        B[2, 2*a] = dN_dy[a]
        B[2, 2*a+1] = dN_dx[a]

    return B, det_J


def _compute_B_and_detJ_quad9(coords, xi, eta):
    """Compute B matrix (3,18) and det(J) for 9-node quad at (xi, eta)."""
    dN_dxi, dN_deta = compute_quad9_shape_derivatives(xi, eta)

    J = np.zeros((2, 2))
    for a in range(9):
        J[0, 0] += dN_dxi[a] * coords[a, 0]
        J[0, 1] += dN_dxi[a] * coords[a, 1]
        J[1, 0] += dN_deta[a] * coords[a, 0]
        J[1, 1] += dN_deta[a] * coords[a, 1]

    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]

    if abs(det_J) < 1e-14:
        return np.zeros((3, 18)), 0.0

    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J

    dN_dx = np.zeros(9)
    dN_dy = np.zeros(9)
    for a in range(9):
        dN_dx[a] = J_inv[0, 0] * dN_dxi[a] + J_inv[0, 1] * dN_deta[a]
        dN_dy[a] = J_inv[1, 0] * dN_dxi[a] + J_inv[1, 1] * dN_deta[a]

    B = np.zeros((3, 18))
    for a in range(9):
        B[0, 2*a] = dN_dx[a]
        B[1, 2*a+1] = dN_dy[a]
        B[2, 2*a] = dN_dy[a]
        B[2, 2*a+1] = dN_dx[a]

    return B, det_J


def _elem_dof_indices(elem_nodes, dof_offset=None):
    """Get global DOF indices for element nodes (translational DOFs only).

    If dof_offset is provided, uses it to map node -> global DOF start.
    Otherwise falls back to the 2*node convention (no pile nodes in mesh).
    """
    dof_idx = np.zeros(2 * len(elem_nodes), dtype=int)
    if dof_offset is not None:
        for i, node in enumerate(elem_nodes):
            dof_idx[2*i] = dof_offset[node]
            dof_idx[2*i+1] = dof_offset[node] + 1
    else:
        for i, node in enumerate(elem_nodes):
            dof_idx[2*i] = 2 * node
            dof_idx[2*i+1] = 2 * node + 1
    return dof_idx




def _ssrm_progress(callback, done, total, label):
    """Report SSRM progress; never let a misbehaving callback break the solve."""
    if callback is not None:
        try:
            callback(int(done), int(total), str(label))
        except Exception:
            pass


def _verdict_note(sol):
    """One-line trial outcome for the SSRM log, naming the hybrid verdict when the
    trial did not converge. Reads the same on the default criterion as it always
    has ("Converged" / "Did NOT converge") with the displacement evidence appended."""
    if sol.get("converged"):
        return "Converged"
    if sol.get("exit_reason") == 'inconclusive':
        return "INCONCLUSIVE at the iteration ceiling (still improving)"
    v = sol.get("verdict") or "FAILED"
    ur = sol.get("u_ratio")
    ur_txt = "" if ur is None else f", max|u| = {ur:.2f}x elastic"
    if v == 'STABLE_STUCK':
        return f"Did NOT converge but STABLE_STUCK -> counted STABLE{ur_txt}"
    if v == 'AMBIGUOUS':
        return f"Did NOT converge (AMBIGUOUS evidence -> failed{ur_txt})"
    return f"Did NOT converge ({v}{ur_txt})"


def _failed_edge_softened(solution):
    """The softened 1D-element set of a failed trial, or None when it has none."""
    if not solution:
        return None
    soft = np.asarray(solution.get("softened_1d_elements", []), dtype=bool)
    return soft if soft.size and soft.any() else None


def _ssrm_bisect_steps(width, tolerance):
    """Number of halvings needed to narrow ``width`` below ``tolerance``."""
    if width <= tolerance:
        return 0
    return max(1, int(np.ceil(np.log2(width / tolerance))))


def _element_centroids(fem_data):
    """(n_elements, 2) array of element centroids: the mean of each element's active
    nodes, matching the convention used when writing element results (save_fem_results)."""
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    centroids = np.zeros((len(elements), 2))
    for i, elem_nodes in enumerate(elements):
        active = elem_nodes[:element_types[i]]
        centroids[i] = np.mean(nodes[active], axis=0)
    return centroids


def _ssr_zone_exclusion_mask(fem_data, zone):
    """Per-element boolean mask flagging elements to EXCLUDE from strength reduction
    because their centroid lies OUTSIDE the SSR search-area polygon ``zone``.

    ``zone`` is an (x, y) vertex list (RS2's "SSR Search Area"; strength reduction
    applies only INSIDE it). Elements whose centroid is inside get mask = False
    (reduced); those outside get mask = True (held at full strength) — the same
    True = excluded convention as ssr_exclude_mask, so the two compose by union.
    """
    poly = Polygon([(float(x), float(y)) for x, y in zone])
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty or poly.area <= 0:
        raise ValueError("ssr_zone polygon is degenerate (zero area); "
                         "supply at least three non-collinear vertices.")
    cen = _element_centroids(fem_data)
    inside = shapely.contains_xy(poly, cen[:, 0], cen[:, 1])
    return ~np.asarray(inside, dtype=bool)


def _zone_union(zones, kinds):
    """Shapely union of every zone in ``zones`` whose 'kind' is in ``kinds``; None
    when no zone matches. Each zone's 'polygon' is a plain (x, y) vertex list."""
    polys = []
    for z in zones:
        if str(z.get("kind", "")).strip() not in kinds:
            continue
        poly = Polygon([(float(x), float(y)) for x, y in z["polygon"]])
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty or poly.area <= 0:
            raise ValueError(
                f"SSR zone ({z.get('kind')}) is degenerate (zero area); supply at "
                "least three non-collinear vertices.")
        polys.append(poly)
    if not polys:
        return None
    return unary_union(polys)


def _compose_ssr_zone_masks(fem_data, zones):
    """Compose the file-carried SSR zone overlays into per-element masks.

    The three overlay kinds (fileio.SSR_ZONE_SENTINELS, displayed as "SSR reduce" /
    "SSR hold" / "SSR elastic") compose by ONE rule:

        reduce set = union(reduce zones) MINUS union(hold zones + elastic zones)

    with the reduce set defaulting to the WHOLE DOMAIN when no reduce zone is drawn.
    So exclusions always carve out — of a search area they sit inside, and of the
    model as a whole when no search area was drawn — and an interior hole is
    expressed simply by drawing the exclusion on top of the search area.

    Membership is by ELEMENT CENTROID, exactly as the ``ssr_zone`` kwarg does: the
    zones are overlays, so the mesh is never made to conform to them, and a zone
    added to a file cannot move a node. That is what makes a file-carried zone and
    the same polygon passed as a kwarg produce bit-identical answers.

    Returns:
        tuple(ndarray or None, ndarray or None): (exclusion mask, elastic mask).
        The exclusion mask uses the same True = held-at-full-strength convention as
        ``ssr_exclude_mask``, so it composes with the other exclusions by union. The
        elastic mask is True inside the "SSR elastic" zones and composes with
        ``elastic_materials`` the same way. Either is None when nothing constrains it.
    """
    if not zones:
        return None, None

    cen = _element_centroids(fem_data)
    x, y = cen[:, 0], cen[:, 1]

    reduce_u = _zone_union(zones, {"reduce"})
    elastic_u = _zone_union(zones, {"hold_elastic"})
    exclude_u = _zone_union(zones, {"hold", "hold_elastic"})

    inside_reduce = (np.ones(len(cen), dtype=bool) if reduce_u is None
                     else np.asarray(shapely.contains_xy(reduce_u, x, y), dtype=bool))
    inside_exclude = (np.zeros(len(cen), dtype=bool) if exclude_u is None
                      else np.asarray(shapely.contains_xy(exclude_u, x, y), dtype=bool))

    # Held at full strength = NOT in the reduce set = outside every search area, or
    # inside an exclusion.
    exclusion_mask = ~(inside_reduce & ~inside_exclude)
    if not exclusion_mask.any():
        exclusion_mask = None          # inert: the overlay constrains nothing

    elastic_mask = None
    if elastic_u is not None:
        elastic_mask = np.asarray(shapely.contains_xy(elastic_u, x, y), dtype=bool)
        if not elastic_mask.any():
            elastic_mask = None

    return exclusion_mask, elastic_mask


def solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.01, debug_level=0, force_tol=1e-3,
               oob_window=10,
               max_iterations=12000, max_iterations_ceiling=50000,
               convergence_tol=1e-3, max_disp_factor=0.1,
               failure_criterion="hybrid", n_sweep=10,
               tension_cutoff=False, char_point=None,
               dt_scale=1.0, cancel_check=None,
               progress_callback=None,
               f_adjust=0.25, f_min_floor=0.1, f_max_ceiling=10.0, max_expand=20,
               grid=None, min_slip_depth=None, ssr_exclude=None, ssr_zone=None,
               tension_cutoff_by_material=None, tension_srf=None, k0=None,
               elastic_materials=None,
               suction_phi_b=None, suction_cap=None,
               capture_failure_state=True, capture_max_iterations=None,
               capture_margin=0.15, early_failure=True, fem_solver=None,
               ssrm_driver='bisection'):
    """
    Shear Strength Reduction Method using bisection on solve_fem convergence.

    Finds the critical strength reduction factor F where the viscoplastic
    algorithm transitions from converging (stable) to non-converging (failure).

    Parameters:
        fem_data (dict): FEM data from build_fem_data
        F_min (float): Lower bound for F (must converge). Default 1.0. If it does
            NOT converge, the bracket auto-expands downward (see f_adjust).
        F_max (float): Upper bound for F (should not converge). Default 2.0. If it
            DOES converge, the bracket auto-expands upward (see f_adjust).
        f_adjust (float): Step by which the bracket is widened when the guess is
            off — F_min lowered / F_max raised by f_adjust and re-checked, so a
            wrong [F_min, F_max] still finds the FS instead of aborting. A good
            guess brackets on the first try and skips this. Default 0.25.
        f_min_floor (float): F_min is never lowered below this (F stays positive).
            Default 0.1; failing to converge even here means FS < f_min_floor.
        f_max_ceiling (float): F_max is never raised above this. Default 10.0;
            still converging here means FS exceeds the ceiling (or the slope
            deforms ductilely without a displacement catastrophe).
        max_expand (int): Cap on expansion steps in each direction. Default 20.
        tolerance (float): Bisection stops when F_right - F_left < tolerance. Default 0.01.
            The reported FS is the midpoint of the final bracket (+/- tolerance/2);
            the bracket itself is returned in 'final_interval'.
        grid (float or None): If set, bisect over a FIXED global grid of step ``grid``
            (F = i*grid) instead of halving the supplied bracket. Every starting
            bracket then converges to the same global cell straddling the failure
            threshold, so the reported FS is INDEPENDENT of F_min/F_max (identical to
            every decimal, not just +/- tolerance/2). ``grid`` becomes the precision
            (cell width). Default None = continuous bisection (bracket-dependent to
            +/- tolerance/2). Used by ``reliability_fem`` for reproducible results.
        min_slip_depth (float or None): Optional surficial-failure filter, threaded to
            solve_fem (default None = off). Excludes failures shallower than this depth
            below the ground surface, so a shallow cohesionless skin does not govern the
            SSRM factor of safety. Off by default; see solve_fem for the full description.
            It steers the mechanism by DEPTH; ssr_exclude / ssr_zone steer it by REGION.
            The filter changes only the failure verdict (which nodes may declare the
            slope failing) — no element is masked and no strength is held back — so a
            shallow zone still yields, it just cannot decide the bisection alone. Set the
            same value on the LEM search (search.py) to compare like-for-like surfaces.
        debug_level (int): Verbosity (0=silent, 1=summary, 2=detailed)
        max_iterations (int): Viscoplastic iteration BUDGET per trial, passed to
            solve_fem (default 12000). A trial that reaches it with the
            out-of-balance residual still trending down is extended by another
            budget's worth, repeatedly, up to max_iterations_ceiling.
        max_iterations_ceiling (int): Hard stop on that extension (default 50000).
            A trial that reaches the ceiling while still improving is INCONCLUSIVE:
            it is not counted as a failure. The bisection carries on below it, the
            factor of safety is reported as the final bracket's midpoint exactly as
            on any other run, and the result carries 'inconclusive' and a 'note'
            recording that the bracket's upper edge is undecided rather than a
            measured failure.

        convergence_tol (float): Convergence tolerance passed to solve_fem
        max_disp_factor (float): Displacement limit (fraction of mesh height) used as a
            backstop/early-termination cap inside solve_fem trials (default 0.1).
        failure_criterion (str): How to determine failure. Default "hybrid"
            (since 2026-07-26).
            "non_convergence" - Bisection on TRUE viscoplastic
                equilibrium: a trial converges only if both the CHECON
                displacement test and the force-equilibrium test are satisfied
                — the latter being the maximum over nodes of |out-of-balance
                force| / |nodal body force|, below `force_tol` (Dawson, Roth &
                Drescher 1999). Being per-node, it cannot be diluted by padding
                the mesh with inert foundation or runout. It is NOT independent
                of element size (it goes as ~1/h), and it needs roughly 3x the
                iterations the old rate-based test did, because it demands real
                equilibrium rather than a decayed rate. Use for problems without
                reservoir loading.
            "hybrid" (DEFAULT) - Same bisection as "non_convergence", but a trial
                that fails to reach equilibrium must also show DISPLACEMENT
                EVIDENCE of failure before the bisection counts it as failed:
                max|u| beyond the trial's own elastic scale AND still growing.
                A trial frozen at elastic scale is reported STABLE_STUCK and
                treated as standing, which stops a numerically stuck solve from
                being read as a failed slope. Everything in between keeps the
                "non_convergence" verdict, flagged AMBIGUOUS. Per-trial verdicts
                are returned in result['trials']. See classify_nonconvergence.
                It is the default on corpus A/B evidence: over all 103 FEM
                benchmarks (2026-07-26) it moved 4 rows and NO row downward —
                three small upward shifts where the legacy criterion truncated
                the bracket early on a stuck-but-standing trial, plus one
                outright rescue (RS2-48, a T=0 reinforced fill whose trials are
                stationary rather than collapsing, which the legacy machinery
                cannot bracket at all). The other 99 rows are unchanged.
            "displacement_limit" - Bisection on whether the max VP displacement
                exceeds max_disp_factor x mesh height within the iteration
                budget. A simple physical backstop; verdict is coupled to the
                iteration budget for slowly creeping states.
            "displacement_increase" - Displacement-catastrophe sweep (cf. Sun,
                Wang & Zhang 2021): locate the upturn of displacement vs F,
                measured at a characteristic point on the failure mechanism
                (auto-selected as the node with the largest plastic-displacement
                growth across the sweep, or supplied via char_point). Produces
                the displacement-vs-F evidence curve; also robust against any
                localized background creep contaminating global measures.
        char_point (tuple or None): (x, y) of a characteristic point for the
            "displacement_increase" measure. Default None = automatic
            selection after the coarse sweep.
        n_sweep (int): Number of points in coarse sweep for "displacement_increase". Default 10.
        ssr_exclude (list of str or None): Names of material zones to EXCLUDE from
            strength reduction. Excluded zones keep their full c and tan(phi) at
            every trial F while the rest are reduced — RS2's per-material Apply_SSR
            flag / "SSR Exclusion Area". Used to force the mechanism up out of a
            zone (e.g. a stiff foundation) so the SSR is evaluated on the surface
            of interest instead of a deeper true-global minimum. Names must match
            the fem_data material names exactly; an unknown name raises ValueError.
            Default None = every zone reduced (today's behavior, bit-identical).
        ssr_zone (list of (x, y) or None): An "SSR Search Area" polygon — strength
            reduction is applied ONLY to elements whose centroid lies INSIDE this
            polygon; everything outside keeps full strength (F = 1). This is RS2's
            SSR-Search-Area constraint, whose native models store it as an exact
            vertex polygon; it confines the strength-reduction mechanism to a chosen
            region (e.g. a band around a proposed slip surface, or one local-minimum
            face) instead of searching the whole domain. Given as a vertex list
            [(x1, y1), (x2, y2), ...] in the model's coordinate system (a closing
            repeat of the first vertex is allowed; the ring is closed automatically).
            Composes with ssr_exclude by UNION of exclusions: an element is held at
            full strength if it is named-excluded OR outside the zone. Default None
            = no search-area constraint from this argument.
            ONE SENSE ONLY -- reduce INSIDE. A vendor polygon flagged as an SSR
            EXCLUSION area (RS2 records the kind per polygon; the downstream shell of
            the RS2-4 Talbingo model is one) means reduce everywhere BUT inside, and
            must be passed as its COMPLEMENT within the model outline. Handing an
            exclusion ring over as-drawn reduces exactly the wrong region.
            FILE-CARRIED TWIN: the input template says the same thing on the polygon
            sheet, as rows with a sentinel Mat ID -- -1 "SSR reduce", -2 "SSR hold",
            -3 "SSR elastic" (fileio.SSR_ZONE_SENTINELS). Those arrive on fem_data as
            ``ssr_zones`` and apply whenever this argument is None -- so "default
            None" means "no polygon FROM THIS ARGUMENT", not necessarily "no zone".
            They are ANALYSIS OVERLAYS: never meshed, never material regions, never
            slice-generating, and classified by the same element-centroid test used
            here, so the same polygon carried by the file and passed as this kwarg
            give bit-identical answers. Composition among the file's own rows:
            reduce set = union(-1 rows) minus union(-2 and -3 rows), defaulting to
            the whole domain when no -1 row is drawn, so an exclusion always carves
            out (of a search area or of the model). -3 rows ADDITIONALLY join the
            elastic set (see elastic_materials). When BOTH this kwarg and the file's
            zones are present, the kwarg wins and the file's zones are ignored, with
            a warning: an explicit polygon has named the search area precisely, and
            silently intersecting it with the file's zones would quietly shrink it.
        tension_cutoff_by_material (dict or None): Per-material tensile-strength
            cutoff as {material name -> T} (T in the model's stress units). Each
            named material's elements get a RANKINE cap on their major principal
            stress at T (second viscoplastic yield surface; see solve_fem's
            tension_cap_by_elem). Names must match a material 'name' exactly.
            None (default) = no per-material cutoff (bit-identical to the path
            without it). This is a RUN OPTION only — it reads nothing from the
            material template.
        tension_srf (bool): Whether the tensile cutoff T is reduced with the trial
            SRF (RS2's ``tensilestrength_SRF``). True (DEFAULT): T -> T/F each
            trial, like c and tan(phi), so FS is the factor by which the whole
            strength envelope — shear and tensile — is reduced. This matches RS2
            (``tensilestrength_SRF=1``, the setting behind its entire published
            verification set) and Plaxis. False: T held at its authored value
            through the bisection, reproducing the pre-2026-07 static-cap behavior
            when this argument defaulted off.
            NO-OP WITHOUT A CAP: it only acts where a cap exists — materials named
            in tension_cutoff_by_material, the file's own ``mat`` t_cut column, or
            the global tension_cutoff. A model with no cap anywhere has nothing to
            reduce and solves bit-identically either way.
        k0 (float or None): At-rest lateral earth pressure coefficient for the FEM
            initial stress state, passed straight through to every trial's
            solve_fem (and to the shared prepared model, which computes the
            overburden once). None (default) = gravity turn-on, i.e. the initial
            lateral stress is the elastic nu/(1-nu) * sigma_v -- an under-confined
            state for thin structural columns such as a reinforced-soil block.
            See solve_fem's ``k0`` for the full formulation.
            When set, the analysis runs in two steps. ONE full-strength solve first
            equilibrates the K0 field against the geometry — on a slope the field is
            not in equilibrium as built, and a sizeable share of the weight has to
            redistribute — and every bisection trial then starts from that
            equilibrated stress state with a zero displacement datum and reduces
            strength from there. In-situ stress and strength reduction are thus
            separate steps, as they are in the field: without the separation each
            trial performs the redistribution against soil already weakened by F and
            charges the resulting plastic strain and displacement to the trial.
            The extra solve is run once and shared by all trials. It must come back
            STABLE on the run's own failure criterion; if it does not, the slope does
            not stand at full strength with that initial stress (FS < 1), a warning
            is issued and the bisection proceeds without a carried in-situ state.
            The result dict carries ``k0_equilibration`` with the outcome either way.
        elastic_materials (list of str or None): Material names whose elements are
            treated as PURE LINEAR ELASTIC — they skip the plastic-correction loop
            entirely and can never yield, mirroring RS2's "Plasticity
            Specifications: None" placed materials. Their stress is the linear
            elastic D*B*u throughout the reduction; they carry no strength surface,
            are never reduced, and never appear as failed. This is DISTINCT from
            ssr_exclude: an ssr_exclude material keeps its FULL strength but STILL
            yields once its (un-reduced) envelope is reached, so it can shed load
            and localize a mechanism; an elastic_materials material has no envelope
            at all and cannot fail under any stress. Composes with ssr_exclude and
            ssr_zone (independent masks). Names must match a material 'name'
            exactly (unknown names raise ValueError). None (default) = the file's own
            ``option = elastic`` materials, if any (bit-identical to the path without
            it when there are none); pass [] to force the elastic set empty.
            POLYGON-ADDRESSED TWIN: a v20 "SSR elastic" (-3) polygon row joins this
            same mask by union, naming the region by outline instead of by material.
        suction_phi_b (dict or None): OPT-IN matric-suction strength (Fredlund
            extended Mohr-Coulomb), {material name: phi_b degrees}, threaded
            unchanged to every solve_fem trial. Above the water table the pore
            pressure is negative (matric suction); for a material named here the
            suction s = max(0, -u) becomes an apparent cohesion s*tan(phi_b) added
            to c' in the MC yield, REDUCED by the trial F alongside c'/tan(phi').
            The effective-normal pore pressure stays clamped at 0 exactly as today.
            None (default) auto-wires from the v17 template phi_b column; an explicit
            dict overrides the file; {} forces suction off. Off by default =>
            bit-identical to the pre-suction solver. Both RS2 (Verification #28,
            Ng & Shi 1998) and SIGMA/W credit suction through the SRF-reduced
            frictional strength, hence the /F reduction here.
        suction_cap (float, dict, or None): Upper bound on the credited suction s
            (stress units) before it becomes apparent cohesion — one scalar for
            every material or a {name: cap} dict. None (default) auto-wires from the
            v17 template s_cap column (uncapped where blank). Ignored when
            suction_phi_b resolves empty.
        capture_failure_state (bool): After the bracket resolves, re-solve ONCE just
            beyond critical (see capture_margin) with the displacement cap OFF and the
            early divergence-exit OFF, letting the unconverged viscoplastic field run
            to a generous iteration ceiling so the failure MECHANISM accumulates and
            dominates the elastic baseline. The resulting field is stored as
            result['failure_solution']. This is the at-failure (unconverged) deformed
            state Griffiths & Lane plot — a rotational mechanism, not the sub-critical
            settlement of the last CONVERGED trial. Purely a rendering extra: FS, the
            bracket, and last_solution are UNAFFECTED, so with this off (or absent)
            results are bit-identical to before. Default True (the figure path); pass
            False on bulk paths that never render the field (reliability, sensitivity)
            to skip the extra solve. Costs one additional non-converging solve per SSRM.
        capture_margin (float): Proportional strength margin above the critical factor
            of safety at which the failure-state field is captured: F = FS x (1 +
            capture_margin), floored at the bracket's failed edge. A margin is needed
            because the bisection resolves the failed edge to within `tolerance` of
            critical, where the viscoplastic runaway is far too slow to develop the
            mechanism in any practical iteration budget — right at the edge the field
            still reads as diffuse settlement. A modest PROPORTIONAL margin (scale-free
            in FS) puts the solve into the developed-mechanism regime, reconstructing
            Griffiths & Lane's convention of plotting the unconverged state a discrete
            F-step beyond critical. Default 0.15 (reproduces the paper's rotational
            deformed mesh); lower toward ~0.05 to stay nearer critical, higher for a
            bolder slip band. Ignored when capture_failure_state is False.
        capture_max_iterations (int or None): Iteration ceiling for the failure-state
            capture solve. None (default) uses max(max_iterations, 3000) — a generous
            budget so the mechanism develops fully. Ignored when
            capture_failure_state is False.
        early_failure (bool): Close a trial as failed as soon as its movement is
            unambiguously running away, instead of spending the rest of its budget
            (default True; see solve_fem). It applies to every bisection trial and
            never to the at-failure capture solve. The 'displacement_increase'
            criterion does not use it at all: that search reads how far each trial
            moves rather than whether it failed, and cutting a trial short shortens
            the very displacement its curve is made of.

    Returns:
        dict: Result with keys FS, converged, last_solution, final_interval, and —
            when capture_failure_state is on — failure_solution (the at-failure
            unconverged field for the deformation/vector figures).
    """

    t_start = time.perf_counter()

    # Template-carried defaults (v16): a t_cut column / an option=elastic material
    # read from the input file populates fem_data['tension_cutoff_by_material'] /
    # ['elastic_materials'] at build time. Honor them automatically when the caller
    # left the run option None; an explicit kwarg wins (pass {} / [] to disable).
    # The substituted values flow through the SAME resolution below as an explicit
    # kwarg, so file-carried and explicitly-passed are bit-identical.
    if tension_cutoff_by_material is None:
        tension_cutoff_by_material = fem_data.get("tension_cutoff_by_material") or None
    if elastic_materials is None:
        elastic_materials = fem_data.get("elastic_materials") or None

    # Template-carried defaults (v19), same rule: kwarg > file > engine default.
    # k0 is threaded on to every trial and to the prepared model; tension_srf is
    # resolved to a concrete bool HERE so the trials never re-resolve it (they are
    # handed fem_data_trials, and a second resolution is how a value drifts).
    if k0 is None:
        k0 = fem_data.get("k0")
    if tension_srf is None:
        tension_srf = fem_data.get("tension_srf")
    if tension_srf is None:
        tension_srf = True                      # engine default

    # Resolve the SSR-exclusion material names to a per-element boolean mask once,
    # up front, so every trial solve shares it. Excluded elements keep full strength
    # (F = 1) inside solve_fem; see the F_by_elem note there.
    ssr_exclude_mask = None
    if ssr_exclude:
        material_names = list(fem_data.get("material_names", []))
        element_materials = fem_data["element_materials"]
        wanted = [str(n).strip() for n in ssr_exclude]
        unknown = [n for n in wanted if n not in material_names]
        if unknown:
            raise ValueError(
                f"ssr_exclude names not found in the model materials {material_names}: "
                f"{unknown}. Names must match a material's 'name' field exactly.")
        excluded_ids = {material_names.index(n) + 1 for n in wanted}  # 1-based mat IDs
        ssr_exclude_mask = np.isin(element_materials, list(excluded_ids))
        if debug_level >= 1:
            print(f"  SSR exclusion: {wanted} "
                  f"({int(ssr_exclude_mask.sum())}/{len(element_materials)} elements "
                  f"held at full strength)")

    # An "SSR Search Area" polygon confines strength reduction to its INTERIOR:
    # elements whose centroid lies OUTSIDE the zone keep full strength (mask = True),
    # matching RS2's SSR-Search-Area semantics. Same True = excluded convention as
    # ssr_exclude_mask, so the two compose by union (an element is held at full
    # strength if it is named-excluded OR outside the zone).
    # File-carried SSR zone overlays (v20 polygon-sheet sentinel rows, arriving on
    # fem_data as 'ssr_zones'). They compose among themselves by the search-area-
    # minus-exclusions rule in _compose_ssr_zone_masks, then join ssr_exclude_mask by
    # union like every other exclusion. An explicit ssr_zone polygon kwarg wins,
    # loudly: a caller who names the search area precisely should not have it
    # silently intersected with the file's own zones.
    _file_zones = list(fem_data.get("ssr_zones") or [])
    _zone_elastic_mask = None
    if _file_zones and ssr_zone is not None:
        warnings.warn(
            f"An explicit ssr_zone polygon was given AND the file carries "
            f"{len(_file_zones)} SSR zone polygon(s) "
            f"({[z.get('kind') for z in _file_zones]}). The kwarg polygon wins; the "
            "file's zones are ignored for this run.")
    elif _file_zones:
        _zone_excl, _zone_elastic_mask = _compose_ssr_zone_masks(fem_data, _file_zones)
        if _zone_excl is not None:
            ssr_exclude_mask = (_zone_excl if ssr_exclude_mask is None
                                else (ssr_exclude_mask | _zone_excl))
        if debug_level >= 1:
            from .fileio import SSR_ZONE_LABELS
            _kinds = [SSR_ZONE_LABELS.get(z.get("kind"), z.get("kind"))
                      for z in _file_zones]
            _n_out = 0 if _zone_excl is None else int(_zone_excl.sum())
            _n_el = 0 if _zone_elastic_mask is None else int(_zone_elastic_mask.sum())
            print(f"  SSR zones (file): {_kinds} "
                  f"({_n_out}/{len(fem_data['element_materials'])} elements held at "
                  f"full strength, {_n_el} held linear elastic)")

    if ssr_zone is not None:
        zone_mask = _ssr_zone_exclusion_mask(fem_data, ssr_zone)
        ssr_exclude_mask = (zone_mask if ssr_exclude_mask is None
                            else (ssr_exclude_mask | zone_mask))
        if debug_level >= 1:
            print(f"  SSR search area: {len(list(ssr_zone))}-vertex polygon "
                  f"({int(zone_mask.sum())}/{len(zone_mask)} elements outside, held "
                  f"at full strength)")

    # Resolve the per-material tensile-strength cutoff dict {name -> T} to a
    # per-element cap array once, up front (like ssr_exclude_mask). inf = no cap.
    # solve_fem reduces it with the trial F when tension_srf is set, so it is built
    # here at full (un-reduced) value.
    tension_cap_by_elem = None
    if tension_cutoff_by_material:
        material_names = list(fem_data.get("material_names", []))
        element_materials = fem_data["element_materials"]
        wanted = {str(n).strip(): float(T) for n, T in tension_cutoff_by_material.items()}
        unknown = [n for n in wanted if n not in material_names]
        if unknown:
            raise ValueError(
                f"tension_cutoff_by_material names not found in the model materials "
                f"{material_names}: {unknown}. Names must match a material's 'name' "
                f"field exactly.")
        tension_cap_by_elem = np.full(len(element_materials), np.inf)
        for name, T in wanted.items():
            mid = material_names.index(name) + 1  # 1-based material ID
            tension_cap_by_elem[element_materials == mid] = T
        if debug_level >= 1:
            print(f"  Tensile cutoff (SRF={'on' if tension_srf else 'off'}): "
                  f"{wanted} "
                  f"({int(np.isfinite(tension_cap_by_elem).sum())}/"
                  f"{len(element_materials)} elements capped)")

    # Resolve the elastic-materials names to a per-element boolean mask once, up
    # front (like ssr_exclude_mask). These materials are held out of plasticity
    # ENTIRELY inside solve_fem — pure linear elastic, cannot yield. DISTINCT from
    # ssr_exclude (full strength but still yields); composes independently with it
    # and with ssr_zone. inf/True = elastic.
    elastic_mask = None
    if elastic_materials:
        material_names = list(fem_data.get("material_names", []))
        element_materials = fem_data["element_materials"]
        wanted = [str(n).strip() for n in elastic_materials]
        unknown = [n for n in wanted if n not in material_names]
        if unknown:
            raise ValueError(
                f"elastic_materials names not found in the model materials "
                f"{material_names}: {unknown}. Names must match a material's 'name' "
                f"field exactly.")
        elastic_ids = {material_names.index(n) + 1 for n in wanted}  # 1-based IDs
        elastic_mask = np.isin(element_materials, list(elastic_ids))
        if debug_level >= 1:
            print(f"  Elastic (no plasticity): {wanted} "
                  f"({int(elastic_mask.sum())}/{len(element_materials)} elements "
                  f"held pure linear elastic)")

    # An "SSR elastic" (-3) overlay row joins the elastic set by union — the same
    # machinery as elastic_materials, addressed by polygon instead of by name. It is
    # ALSO in ssr_exclude_mask (it is an exclusion), which is consistent: an element
    # that cannot yield is not reduced either.
    if _zone_elastic_mask is not None:
        elastic_mask = (_zone_elastic_mask if elastic_mask is None
                        else (np.asarray(elastic_mask, dtype=bool) | _zone_elastic_mask))

    # Warn about volumetric locking with low-order elements
    element_types = fem_data['element_types']
    has_linear = any(t in (3, 4) for t in element_types)
    if has_linear:
        print("\n" + "!" * 72)
        print("!  WARNING: VOLUMETRIC LOCKING — RESULTS MAY BE UNCONSERVATIVE")
        print("!" * 72)
        print("!  This mesh contains low-order elements (tri3 and/or quad4).")
        print("!  These elements have too few DOFs to represent the nearly")
        print("!  incompressible plastic strains produced by Mohr-Coulomb")
        print("!  yielding, causing an artificially stiff response that")
        print("!  overestimates the factor of safety by 10-20% or more.")
        print("!")
        print("!  Use quadratic elements (tri6, quad8, or quad9) for reliable SSRM")
        print("!  results. The default element type quad8 is recommended.")
        print("!" * 72 + "\n")

    # solve_ssrm is the SOLE auto-wiring point for the SSRM path: it resolved the
    # (auto-wired-or-explicit) tension_cutoff_by_material / elastic_materials into
    # the per-element arrays above and passes them to every trial explicitly. Hand
    # the trials a fem_data with the by-material template DEFAULTS stripped, so
    # solve_fem's own direct-call fallback cannot RE-apply them — which would both
    # double-count and, worse, defeat an explicit disable ({} / []), where the
    # resolved arrays are None but the defaults still sit in fem_data.
    # The v19 file-carried options are stripped for the same reason: solve_ssrm has
    # already resolved them (kwarg > file > default) and passes concrete values to
    # every trial, so leaving them on fem_data would give solve_fem a second chance
    # to re-resolve — the path by which an explicitly-disabled option comes back.
    fem_data_trials = {k: v for k, v in fem_data.items()
                       if k not in ("tension_cutoff_by_material", "elastic_materials",
                                    "k0", "tension_srf", "ssr_zones")}

    # Build the strength-reduction-factor-INDEPENDENT prepared model ONCE and share it
    # across every trial (and the capture solve). The trials differ only in F, which
    # this setup does not touch — the K factorization, the geometry precompute, the
    # pore-pressure fields and the Dawson g_node normalization are all reused instead
    # of being rebuilt ~10 times. Built on fem_data_trials with the same F-independent
    # options the trials pass, so it can never serve a stale strength or geometry.
    # (max_disp_factor and tension_srf are per-call scalings applied inside solve_fem,
    # so they are intentionally NOT part of the prepared model.)
    prep = _prepare_fem_model(
        fem_data_trials, dt_scale=dt_scale, suction_phi_b=suction_phi_b,
        suction_cap=suction_cap, elastic_mask=elastic_mask,
        tension_cap_by_elem=tension_cap_by_elem, tension_cutoff=tension_cutoff,
        min_slip_depth=min_slip_depth, k0=k0, debug_level=max(0, debug_level - 1))

    # === K0 in-situ equilibration (once, before the bisection) ===
    # Establishing the in-situ state and reducing the strength are two different
    # steps of the analysis and must not be run as one. The K0 field is built from
    # the overburden alone: it is an exact equilibrium under level ground, but on a
    # slope a substantial fraction of the weight (on the Griffiths & Lane geometry,
    # about a quarter) is left unbalanced and has to redistribute. Left inside the
    # trials, that redistribution happens against SOIL ALREADY WEAKENED BY F, and
    # every trial is charged with the displacement and plastic strain it produces —
    # displacement measured from a configuration the slope was never in (on the G&L
    # geometry at F = 1.2 the reported movement is ~3x the movement the strength
    # reduction actually causes), and a couple of crest elements taking their
    # yielding at reduced strength instead of at the full strength where the K0 field
    # puts them.
    #
    # So: solve ONCE at full strength, let the K0 field settle against the real
    # geometry, and hand that equilibrated stress field to every trial as its initial
    # state. Each trial then starts from the in-situ state with a zero displacement
    # datum and reduces strength from there, which is what strength reduction means.
    # Cost is one extra solve per SSRM, shared by all ~10 trials and by the capture.
    #
    # If the equilibration itself will not converge, the slope does not stand at full
    # strength (FS < 1) and there is no in-situ state to carry. The bisection then
    # runs on the legacy sequencing and finds the sub-unity factor of safety.
    init_state = None
    equilibration = None
    if k0 is not None:
        if debug_level >= 1:
            print(f"  K0 = {float(k0):g}: equilibrating the in-situ stress state at "
                  "full strength before the bisection…")
        # The bisection's progress bar has not started yet, so say what the wait is.
        _ssrm_progress(progress_callback, 0, 1,
                       f"Equilibrating the K0 = {float(k0):g} initial stress state")
        eq = solve_fem(
            fem_data_trials, F=1.0, debug_level=max(0, debug_level - 1),
            force_tol=force_tol, oob_window=oob_window, dt_scale=dt_scale,
            max_iterations=max_iterations,
            max_iterations_ceiling=max_iterations_ceiling,
            tolerance=convergence_tol, max_disp_factor=None,
            tension_cutoff=tension_cutoff, min_slip_depth=min_slip_depth,
            early_exit=True, ssr_exclude_mask=ssr_exclude_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
            elastic_mask=elastic_mask,
            suction_phi_b=suction_phi_b, suction_cap=suction_cap, k0=k0,
            failure_criterion=failure_criterion, early_failure=early_failure,
            # The in-situ state is established on the SAME driver that will solve
            # the trials. It is the state every trial starts from, so equilibrating
            # it with one solver and then continuing it with another would hand the
            # bisection a starting point its own corrector never produced.
            fem_solver=fem_solver,
            _prepared=prep)
        equilibration = {
            "converged": bool(eq["converged"]),
            "stable": bool(eq["stable"]),
            "verdict": eq["verdict"],
            "iterations": int(eq["iterations"]),
            "max_displacement": float(eq["max_displacement"]),
            "n_plastic": int(np.count_nonzero(eq["plastic_elements"])),
            "unbalanced_force_ratio": float(eq["unbalanced_force_ratio"]),
        }
        # The state counts as established on the SAME standard the bisection uses to
        # call a slope standing — `stable`, which under the default criterion also
        # accepts a state frozen at elastic scale. Demanding force equilibrium here
        # would refuse an in-situ state on models that never reach force_tol at ANY
        # strength but sit perfectly still, which is a solver property, not a
        # statement about the slope.
        if eq["stable"]:
            init_state = eq["_k0_state"]
            if debug_level >= 1:
                print(f"    in-situ state established in {eq['iterations']} "
                      f"iterations (max|u| = {eq['max_displacement']:.4g}, "
                      f"{equilibration['n_plastic']} yielded elements); trials start "
                      "from it with a zero displacement datum")
        else:
            warnings.warn(
                "The K0 in-situ stress state could not be equilibrated at full "
                f"strength ({eq['iterations']} iterations, max|u| = "
                f"{eq['max_displacement']:.4g}, verdict {eq['verdict']}): the slope "
                "does not stand under its own weight with this initial stress, so "
                "FS < 1. The bisection runs without a carried in-situ state.")

    # ---- SSRM driver switch (INTERNAL, spike; see SPIKE.md, "RAMP") ---------
    # 'bisection' is the default and the definition of every locked factor of
    # safety; it falls straight through and nothing below it changes. 'ramp' is
    # the monotonic strength-reduction continuation, and it exists only on the
    # Newton per-trial driver — there is no viscoplastic warm start to continue.
    _driver = str(ssrm_driver or 'bisection').strip().lower()
    if _driver not in ('bisection', 'ramp'):
        raise ValueError(
            f"Unknown ssrm_driver {ssrm_driver!r}. Supported: 'bisection' "
            "(default) and 'ramp'.")
    if _driver == 'ramp':
        if resolve_fem_solver(fem_solver) != 'newton':
            raise ValueError(
                "ssrm_driver='ramp' requires fem_solver='newton'. The ramp warm-"
                "starts each strength step from the previous step's converged "
                "state, which the viscoplastic driver does not carry.")
        if failure_criterion not in ("non_convergence", "hybrid"):
            raise ValueError(
                f"ssrm_driver='ramp' does not run failure_criterion="
                f"'{failure_criterion}'. Its verdict is force equilibrium plus the "
                "Newton displacement bound, which is the hybrid/non-convergence "
                "question.")
        result = _ssrm_ramp_newton(
            fem_data_trials, F_min, F_max, prep=prep, force_tol=force_tol,
            convergence_tol=convergence_tol, max_iterations=max_iterations,
            oob_window=oob_window, dt_scale=dt_scale,
            tension_cutoff=tension_cutoff, min_slip_depth=min_slip_depth,
            ssr_exclude_mask=ssr_exclude_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
            elastic_mask=elastic_mask, suction_phi_b=suction_phi_b,
            suction_cap=suction_cap, k0=k0, f_adjust=f_adjust,
            f_min_floor=f_min_floor, max_expand=max_expand,
            max_iterations_ceiling=max_iterations_ceiling,
            early_failure=early_failure, init_state=init_state,
            debug_level=debug_level, progress_callback=progress_callback,
            cancel_check=cancel_check)
    elif failure_criterion in ("non_convergence", "hybrid"):
        # Same driver, same trials — 'hybrid' only changes how a NON-CONVERGED
        # trial's verdict is read (see classify_nonconvergence).
        result = _ssrm_displacement_limit(
            fem_data_trials, F_min=F_min, F_max=F_max, tolerance=tolerance, force_tol=force_tol,
            oob_window=oob_window, hybrid=(failure_criterion == "hybrid"),
            debug_level=debug_level, max_iterations=max_iterations,
            max_iterations_ceiling=max_iterations_ceiling,
            convergence_tol=convergence_tol, max_disp_factor=None,
            tension_cutoff=tension_cutoff, dt_scale=dt_scale,
            cancel_check=cancel_check, progress_callback=progress_callback,
            f_adjust=f_adjust, f_min_floor=f_min_floor, f_max_ceiling=f_max_ceiling,
            max_expand=max_expand, grid=grid, min_slip_depth=min_slip_depth,
            ssr_exclude_mask=ssr_exclude_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
            elastic_mask=elastic_mask,
            suction_phi_b=suction_phi_b, suction_cap=suction_cap, k0=k0,
            early_failure=early_failure, fem_solver=fem_solver,
            _prepared=prep, _init_state=init_state)
    elif failure_criterion == "displacement_limit":
        result = _ssrm_displacement_limit(
            fem_data_trials, F_min=F_min, F_max=F_max, tolerance=tolerance, force_tol=force_tol,
            oob_window=oob_window,
            debug_level=debug_level, max_iterations=max_iterations,
            max_iterations_ceiling=max_iterations_ceiling,
            convergence_tol=convergence_tol, max_disp_factor=max_disp_factor,
            tension_cutoff=tension_cutoff, dt_scale=dt_scale,
            cancel_check=cancel_check, progress_callback=progress_callback,
            f_adjust=f_adjust, f_min_floor=f_min_floor, f_max_ceiling=f_max_ceiling,
            max_expand=max_expand, grid=grid, min_slip_depth=min_slip_depth,
            ssr_exclude_mask=ssr_exclude_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
            elastic_mask=elastic_mask,
            suction_phi_b=suction_phi_b, suction_cap=suction_cap, k0=k0,
            early_failure=early_failure, fem_solver=fem_solver,
            _prepared=prep, _init_state=init_state)
    elif failure_criterion == "displacement_increase":
        result = _ssrm_displacement_increase(
            fem_data_trials, F_min=F_min, F_max=F_max, tolerance=tolerance, force_tol=force_tol,
            oob_window=oob_window,
            debug_level=debug_level, max_iterations=max_iterations,
            convergence_tol=convergence_tol, n_sweep=n_sweep,
            tension_cutoff=tension_cutoff, char_point=char_point, dt_scale=dt_scale,
            cancel_check=cancel_check, progress_callback=progress_callback,
            min_slip_depth=min_slip_depth, ssr_exclude_mask=ssr_exclude_mask,
            tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
            elastic_mask=elastic_mask,
            suction_phi_b=suction_phi_b, suction_cap=suction_cap, k0=k0,
            _prepared=prep, _init_state=init_state)
    else:
        raise ValueError(
            f"Unknown failure_criterion '{failure_criterion}'. Supported: "
            "'non_convergence' (default; bisection on true viscoplastic "
            "equilibrium), 'hybrid' (non-convergence corroborated by "
            "displacement evidence), 'displacement_limit' (displacement-budget "
            "backstop), and "
            "'displacement_increase' (displacement-catastrophe sweep) - see "
            "docs/fem/overview.md, 'Choosing a Failure Criterion'.")

    if equilibration is not None:
        # What the in-situ equilibration cost and found — reported so a run can be
        # read for whether the initial state was established at all.
        result["k0_equilibration"] = equilibration

    # === Post-bracket capture of the at-failure (unconverged) mechanism ===
    # The bisection keeps only the last CONVERGED field, which is sub-critical and
    # reads as diffuse settlement. The deformed-mesh figures Griffiths & Lane plot are
    # the UNCONVERGED runaway at the failure strength — a rotational mechanism. Re-solve
    # ONCE just beyond critical with the displacement cap OFF and the early
    # divergence-exit OFF, letting the mechanism iterate to a generous ceiling so it
    # dominates the elastic baseline. The capture F is FS x (1 + capture_margin),
    # floored at the bracket's failed edge: bisection narrows the failed edge to within
    # `tolerance` of critical, where the runaway is too slow to develop the mechanism
    # in a finite budget, so a modest proportional margin is required to reach the
    # developed regime (see capture_margin). This changes nothing about FS / the
    # bracket / last_solution — it only ADDS 'failure_solution', so results are
    # bit-identical when this is off.
    if (capture_failure_state and result.get("converged")
            and result.get("final_interval") is not None
            and result.get("FS") is not None):
        F_edge = result["final_interval"][1]
        F_fail = max(F_edge, result["FS"] * (1.0 + capture_margin))
        cap_iters = capture_max_iterations or max(max_iterations, 3000)
        if debug_level >= 1:
            print(f"  Capturing at-failure mechanism: solve at F={F_fail:.3f} "
                  f"(FS x {1.0 + capture_margin:.2f}, failed edge {F_edge:.3f}; cap off, "
                  f"early-exit off, up to {cap_iters} iters)…")
        try:
            failure_solution = solve_fem(
                fem_data_trials, F=F_fail, debug_level=max(0, debug_level - 1),
                force_tol=force_tol, oob_window=oob_window, dt_scale=dt_scale,
                max_iterations=cap_iters,
                max_iterations_ceiling=cap_iters, tolerance=convergence_tol,
                max_disp_factor=None,
                tension_cutoff=tension_cutoff, min_slip_depth=min_slip_depth,
                early_exit=False, early_failure=False,
                ssr_exclude_mask=ssr_exclude_mask,
                tension_cap_by_elem=tension_cap_by_elem, tension_srf=tension_srf,
                elastic_mask=elastic_mask,
                suction_phi_b=suction_phi_b, suction_cap=suction_cap, k0=k0,
                _prepared=prep, _init_state=init_state,
                _softened_seed=result.get("failed_edge_softened"))
            result["failure_solution"] = failure_solution
            if debug_level >= 1:
                print(f"    at-failure field: converged={failure_solution['converged']} "
                      f"iters={failure_solution['iterations']} "
                      f"max_disp={failure_solution['max_displacement']:.3g}")
                # The at-failure title no longer names the trial F (see _fs_title in
                # plot_fem.py — "at Failure" already discloses the unconverged state);
                # this line carries that detail into the log instead.
                print(f"    failure snapshot: trial F={F_fail:.3f} "
                      f"(margin {capture_margin:.2f}), "
                      f"{failure_solution['iterations']} iterations, "
                      f"max disp {failure_solution['max_displacement']:.3g}")
        except Exception:
            # A failed capture must never sink a good FS result; the figure path
            # simply falls back to last_solution when failure_solution is absent.
            # (KeyboardInterrupt is a BaseException and still propagates.)
            if debug_level >= 1:
                print("    at-failure capture failed; continuing without it.")

    elapsed = time.perf_counter() - t_start
    result["elapsed_time"] = elapsed
    if debug_level >= 1:
        print(f"  SSRM completed in {elapsed:.1f} seconds")

    return result


def _ssrm_displacement_limit(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05, force_tol=1e-3,
                              oob_window=10, k0=None,
                              debug_level=0, max_iterations=500,
                              max_iterations_ceiling=50000,
                              convergence_tol=1e-3, max_disp_factor=0.1,
                              tension_cutoff=False,
                 dt_scale=1.0, cancel_check=None, progress_callback=None,
                 f_adjust=0.25, f_min_floor=0.1, f_max_ceiling=10.0, max_expand=20,
                 grid=None, min_slip_depth=None, ssr_exclude_mask=None,
                 tension_cap_by_elem=None, tension_srf=False, elastic_mask=None,
                 suction_phi_b=None, suction_cap=None, early_failure=True,
                 fem_solver=None, _prepared=None, _init_state=None, hybrid=False):
    """SSRM using fixed VP displacement limit as failure criterion.

    The [F_min, F_max] bracket auto-expands when the user's guess is off: if F_min
    doesn't converge it is lowered by f_adjust (down to f_min_floor, keeping F
    positive); if F_max converges it is raised by f_adjust (up to f_max_ceiling).
    So a wrong bracket still finds the FS instead of aborting, while a good bracket
    skips the expansion and bisects immediately.

    ``hybrid=True`` (only meaningful with max_disp_factor=None, i.e. the
    non-convergence path) reads each trial's ``stable`` flag instead of its
    ``converged`` flag, so a trial classified STABLE_STUCK — non-converged but
    frozen at elastic displacement scale — moves the bracket UP rather than down.
    Every trial's verdict is recorded in the returned ``trials`` list on both
    settings, so an A/B comparison needs no extra solves."""

    trials = []                        # per-trial verdict metadata (both settings)
    inconclusive = []                  # trials that hit the iteration ceiling, still improving

    # The rescue cost policy's state (SPIKE.md, "THE COST OF THE RESCUE"). With the
    # policy off — the default — none of it is read and every trial is solved the
    # way it was. `_w0` is the bracket the search was ASKED for, which is the natural
    # unit for "far above a standing bound"; `_carried` is the highest strength shown
    # to stand so far, and is None until the lower bound is established, so the
    # walk-down that looks for one is never budgeted against a bound that does not
    # exist yet.
    _w0 = max(float(F_max) - float(F_min), 1e-9)
    _carried = [None]
    _seed_carried = [False]            # has any trial on this model been seed-carried?

    def _rescue_policy(F):
        """(rung budget, seed-first) for a trial at F, from the standing bracket."""
        rungs = None
        if _carried[0] is not None:
            d = (float(F) - _carried[0]) / _w0
            if _NR_RESCUE_NONE_FRAC is not None and d > _NR_RESCUE_NONE_FRAC:
                rungs = 0
            elif (_NR_RESCUE_ADAPTIVE_FRAC is not None
                    and d > _NR_RESCUE_ADAPTIVE_FRAC):
                rungs = len(_NR_VP_PREDICTOR_ITERS or ())
        return rungs, bool(_NR_SEED_MEMORY and _seed_carried[0])

    def _stable(sol):
        """Does the bisection treat this trial as standing at its F?"""
        return bool(sol.get("stable", sol["converged"])) if hybrid else bool(sol["converged"])

    def _inconclusive(sol):
        """Did this trial run out of CEILING rather than out of progress?

        Such a trial is neither converged nor failed: the residual was still coming
        down when the hard ceiling stopped it. Counting it as a failure is what the
        no-progress exit used to do, and it biases the factor of safety low, so the
        bisection refuses to rule on it. A criterion that CAN rule on it — the
        hybrid's STABLE_STUCK verdict — is left to do so, so this asks only about
        trials the bisection would otherwise have counted as failures."""
        return sol.get("exit_reason") == 'inconclusive' and not _stable(sol)

    def _note_inconclusive(F, sol):
        msg = (f"SSRM: trial F = {F:.4f} is inconclusive at the iteration ceiling "
               f"({sol.get('iterations', 0)} iterations, out-of-balance still "
               f"falling) - raise max_iterations_ceiling to decide it. It is NOT "
               f"counted as a failure: the bracket's upper edge carries this trial "
               f"as an uncertainty rather than a measured failure, and the factor of "
               f"safety is reported as the bracket midpoint, as on any other run.")
        inconclusive.append({"F": float(F), "iterations": int(sol.get("iterations", 0)),
                             "message": msg})
        print(f"\n{msg}")
        return msg

    def _record(F, sol, role):
        trials.append({
            "F": float(F),
            "role": role,
            "converged": bool(sol.get("converged", False)),
            "stable": _stable(sol),
            "verdict": sol.get("verdict"),
            "u_ratio": sol.get("u_ratio"),
            "growth": sol.get("u_growth"),
            "exit_reason": sol.get("exit_reason"),
            # The residual went flat before this trial stopped (None if it never did).
            "plateau_iteration": sol.get("plateau_iteration"),
            # How many extra budgets this trial was granted for still improving.
            "budget_extensions": int(sol.get("budget_extensions", 0) or 0),
            # The iteration the early-failure rule fired at (None if it never did).
            "diverging_iteration": sol.get("diverging_iteration"),
            "ee_suppressed": bool(sol.get("early_exit_suppressed", False)),
            "iterations": int(sol.get("iterations", 0)),
            # Newton-only, and zero on the viscoplastic driver: the constitutive
            # work actually done (see the note at `nr_force_evals`), and the
            # viscoplastic predictor iterations the trial was charged (see
            # `_NR_VP_PREDICTOR_ITERS`). Both are here so a run's total work can be
            # read off the trial record rather than re-derived.
            "nr_force_evals": int(sol.get("nr_force_evals", 0) or 0),
            "nr_predictor_iterations": int(sol.get("nr_predictor_iterations", 0) or 0),
            # Cost attribution, bookkeeping only and never read for a verdict (see
            # SPIKE.md, "THE COST OF THE RESCUE"): this trial's wall time, the
            # standing bracket it was asked from, and — on the Newton driver — the
            # split between the cold attempt, the predictor rungs and the seeded
            # correctors. Zero/absent on the viscoplastic driver.
            "wall": float(sol.get("_trial_wall", 0.0) or 0.0),
            "bracket": (float(F_left), float(F_right)),
            "nr_cold_wall": float(sol.get("nr_cold_wall", 0.0) or 0.0),
            "nr_cold_iterations": int(sol.get("nr_cold_iterations", 0) or 0),
            "nr_cold_force_evals": int(sol.get("nr_cold_force_evals", 0) or 0),
            "nr_rungs": list(sol.get("nr_rungs", []) or []),
            "nr_cold_skipped": bool(sol.get("_nr_cold_skipped", False)),
        })
        if _stable(sol):
            _carried[0] = (float(F) if _carried[0] is None
                           else max(_carried[0], float(F)))
            if any(r.get('converged') for r in (sol.get("nr_rungs") or ())):
                _seed_carried[0] = True
        return sol

    if debug_level >= 1:
        label = ("Hybrid" if hybrid else
                 "Non-Convergence" if max_disp_factor is None else "Displacement Limit")
        print(f"=== SSRM Analysis ({label} Method) ===")
        print(f"  Bisection range: [{F_min:.2f}, {F_max:.2f}], tolerance: {tolerance}")
        if max_disp_factor is not None:
            print(f"  Displacement limit: {max_disp_factor:.0%} of mesh height")

    # Progress reported as: the bracket-establishment solves (>=2) + the
    # (deterministic) bisection steps. Each solve_fem is subdivided (SUBDIV) and its
    # intra-solve progress fills its slice — the bar advances *within* each trial.
    # n_bracket / n_steps live in ``prog`` because auto-expansion (below) can add
    # bracket solves and widen the range, so the total is recomputed on the fly.
    SUBDIV = 100
    prog = {"n_bracket": 2, "n_steps": _ssrm_bisect_steps(F_max - F_min, tolerance)}

    def _total():
        return prog["n_bracket"] + prog["n_steps"]

    def _fem_progress(step, prefix):
        """A solve_fem progress_callback mapping its inner [0, 1] fraction into
        ``step``'s slice of the overall SSRM bar (total recomputed live)."""
        if progress_callback is None:
            return None

        def _cb(frac, info):
            total_fine = _total() * SUBDIV
            s = min(step, _total() - 1)
            pos = (s + max(0.0, min(1.0, frac))) * SUBDIV
            _ssrm_progress(progress_callback, min(int(pos), total_fine - 1),
                           total_fine, f"{prefix} · {info}")
        return _cb

    def _solve_at(F, step, prefix):
        _t0 = time.perf_counter()
        _sol_at = _solve_at_inner(F, step, prefix)
        _sol_at["_trial_wall"] = time.perf_counter() - _t0
        return _sol_at

    def _solve_at_inner(F, step, prefix):
        return solve_fem(fem_data, F=F, debug_level=max(0, debug_level - 1), force_tol=force_tol,
                         oob_window=oob_window,
                         dt_scale=dt_scale,
                         max_iterations=max_iterations, tolerance=convergence_tol,
                         max_iterations_ceiling=max_iterations_ceiling,
                         max_disp_factor=max_disp_factor,
                         tension_cutoff=tension_cutoff, min_slip_depth=min_slip_depth,
                         early_exit=(max_disp_factor is None),
                         ssr_exclude_mask=ssr_exclude_mask,
                         tension_cap_by_elem=tension_cap_by_elem,
                         tension_srf=tension_srf, elastic_mask=elastic_mask,
                         suction_phi_b=suction_phi_b, suction_cap=suction_cap,
                         progress_callback=_fem_progress(step, prefix),
                         failure_criterion=("hybrid" if hybrid else "non_convergence"),
                         k0=k0, early_failure=early_failure,
                         fem_solver=fem_solver,
                         _nr_rescue_rungs=_rescue_policy(F)[0],
                         _nr_seed_first=_rescue_policy(F)[1],
                         _prepared=_prepared, _init_state=_init_state)

    F_left = F_min
    F_right = F_max
    bracket_step = 0                       # running index into the progress bar

    def _bump_bracket():
        nonlocal bracket_step
        bracket_step += 1
        prog["n_bracket"] = max(prog["n_bracket"], bracket_step + 1)

    # === Establish a valid bracket, auto-expanding a wrong guess ===
    # The lower bound must converge and the upper must not. If the guess is off,
    # step F_left DOWN / F_right UP by f_adjust until the bracket is valid, bounded
    # by a positive floor (F stays > 0) and a ceiling (F can't grow forever). A
    # good guess brackets on the first try and skips the expansion entirely.

    # -- Lower bound: lower F_left until it converges --
    _ssrm_progress(progress_callback, 0, _total() * SUBDIV, f"Checking lower bound F={F_left:.3f}")
    if debug_level >= 1:
        print(f"  Verifying lower bound F={F_left:.2f} converges...")
    solution_min = _record(F_left, _solve_at(F_left, bracket_step,
                                             f"Lower bound F={F_left:.3f}"), "lower")
    if _inconclusive(solution_min):
        # An inconclusive lower bound has not been shown to stand, so the bracket
        # walks down exactly as it would for a failure — but the trial is recorded
        # as the uncertainty it is rather than silently counted as a failure.
        _note_inconclusive(F_left, solution_min)
    n_expand = 0
    while not _stable(solution_min):
        if F_left <= f_min_floor + 1e-9 or n_expand >= max_expand:
            msg = (f"SSRM: the slope does not reach equilibrium even at F = {F_left:.2f} "
                   f"(lowered to the floor while auto-bracketing) — it is unstable at or "
                   f"below this strength-reduction factor (FS < {F_left:.2f}).")
            print(f"\n{msg}")
            return {"converged": False, "error": msg, "FS": None, "trials": trials}
        F_new = max(f_min_floor, F_left - f_adjust)
        if debug_level >= 1:
            print(f"    -> F={F_left:.2f} did not converge; lowering F_min to {F_new:.2f}")
        F_left = F_new
        n_expand += 1
        _bump_bracket()
        solution_min = _record(F_left, _solve_at(F_left, bracket_step,
                                                f"Lower bound F={F_left:.3f}"), "lower")
        if _inconclusive(solution_min):
            _note_inconclusive(F_left, solution_min)
    F_min = F_left
    if debug_level >= 1:
        print(f"    -> Converged in {solution_min['iterations']} iters (F_min={F_min:.2f})")

    # -- Upper bound: raise F_right until it does NOT converge --
    _bump_bracket()
    _ssrm_progress(progress_callback, bracket_step * SUBDIV, _total() * SUBDIV,
                   f"Checking upper bound F={F_right:.3f}")
    if debug_level >= 1:
        print(f"  Verifying upper bound F={F_right:.2f} does not converge...")
    solution_max = _record(F_right, _solve_at(F_right, bracket_step,
                                              f"Upper bound F={F_right:.3f}"), "upper")
    if _inconclusive(solution_max):
        # The upper bound only has to NOT stand, and an inconclusive trial does not:
        # it is accepted as the bracket's upper edge, carrying its uncertainty, and
        # the bisection proceeds below it.
        _note_inconclusive(F_right, solution_max)
    n_expand = 0
    while _stable(solution_max):
        if F_right >= f_max_ceiling - 1e-9 or n_expand >= max_expand:
            msg = (f"SSRM: the slope still reaches equilibrium at F = {F_right:.2f} "
                   f"(raised to the ceiling while auto-bracketing), so the factor of safety "
                   f"exceeds this. Raise F_max / f_max_ceiling, or the slope may deform "
                   f"ductilely without a displacement catastrophe.")
            print(f"\n{msg}")
            return {"converged": False, "FS": None, "last_solution": solution_max,
                    "error": msg, "trials": trials}
        F_new = min(f_max_ceiling, F_right + f_adjust)
        if debug_level >= 1:
            print(f"    -> F={F_right:.2f} converged; raising F_max to {F_new:.2f}")
        F_right = F_new
        n_expand += 1
        _bump_bracket()
        solution_max = _record(F_right, _solve_at(F_right, bracket_step,
                                                 f"Upper bound F={F_right:.3f}"), "upper")
        if _inconclusive(solution_max):
            _note_inconclusive(F_right, solution_max)
    F_max = F_right
    if debug_level >= 1:
        print(f"    -> Did NOT converge ({solution_max['iterations']} iters, F_max={F_max:.2f})")

    # Recompute the bisection step budget for the (possibly widened) bracket.
    prog["n_steps"] = _ssrm_bisect_steps(F_right - F_left, tolerance)

    last_converged_solution = solution_min
    # The bracket's failed edge, kept for the at-failure capture: the post-peak
    # set the trial at F_right shed to before it gave way (empty when nothing
    # can soften). Updated whenever a trial moves F_right down.
    failed_edge_solution = solution_max
    iteration = 0
    from .search import _check_cancel

    if grid is not None and grid > 0:
        # === Grid bisection (bracket-independent) ===
        # Bisect over integer indices on a FIXED global grid (F = i*grid) instead of
        # halving the supplied bracket. The failure threshold sits between two global
        # grid points (k*·grid converges, (k*+1)·grid does not) — a fact of the
        # slope+mesh, not of the bracket — so EVERY starting bracket narrows to that
        # same cell and the reported FS is identical to every decimal. Snap the
        # bracket outward (floor lo / ceil hi) so the grid endpoints keep the
        # "lo converges, hi doesn't" invariant by monotonicity (no extra solves).
        import math
        i_lo = int(math.floor(F_left / grid + 1e-9))
        i_hi = int(math.ceil(F_right / grid - 1e-9))
        if i_hi <= i_lo:
            i_hi = i_lo + 1
        prog["n_steps"] = max(1, _ssrm_bisect_steps((i_hi - i_lo) * grid, grid))
        while (i_hi - i_lo) > 1 and iteration < 60:
            _check_cancel(cancel_check)
            i_mid = (i_lo + i_hi) // 2
            F_mid = i_mid * grid
            step = min(prog["n_bracket"] + iteration, _total() - 1)
            lo_f, hi_f = i_lo * grid, i_hi * grid
            _ssrm_progress(progress_callback, step * SUBDIV, _total() * SUBDIV,
                           f"F={F_mid:.3f} in [{lo_f:.3f}, {hi_f:.3f}]")
            if debug_level >= 1:
                print(f"\n  SSRM step {iteration+1} (grid {grid:g}): "
                      f"F = {F_mid:.4f}  [{lo_f:.4f}, {hi_f:.4f}]")
            solution = _record(
                F_mid, _solve_at(F_mid, step, f"F={F_mid:.3f} [{lo_f:.3f}, {hi_f:.3f}]"),
                "bisect")
            if _inconclusive(solution):
                # A trial nothing can rule on becomes the bracket's UPPER
                # UNCERTAINTY: the search carries on below it, so the run still
                # reports the highest F that actually reached equilibrium, but that
                # edge is never counted as a measured failure and never sets the
                # answer by being averaged into it.
                _note_inconclusive(F_mid, solution)
                i_hi = i_mid
                iteration += 1
                continue
            if _stable(solution):
                i_lo = i_mid
                last_converged_solution = solution
                if debug_level >= 1:
                    print(f"    -> {_verdict_note(solution)} ({solution['iterations']} iters)")
            else:
                i_hi = i_mid
                failed_edge_solution = solution
                if debug_level >= 1:
                    print(f"    -> {_verdict_note(solution)} ({solution['iterations']} iters)")
            iteration += 1
        F_left, F_right = i_lo * grid, i_hi * grid
    else:
        while (F_right - F_left) > tolerance and iteration < 50:
            _check_cancel(cancel_check)
            F_mid = (F_left + F_right) / 2.0

            step = min(prog["n_bracket"] + iteration, _total() - 1)
            _ssrm_progress(progress_callback, step * SUBDIV, _total() * SUBDIV,
                           f"F={F_mid:.3f} in [{F_left:.3f}, {F_right:.3f}]")
            if debug_level >= 1:
                print(f"\n  SSRM step {iteration+1}: F = {F_mid:.4f}  [{F_left:.4f}, {F_right:.4f}]")

            solution = _record(
                F_mid,
                _solve_at(F_mid, step, f"F={F_mid:.3f} [{F_left:.3f}, {F_right:.3f}]"),
                "bisect")

            if _inconclusive(solution):
                # A trial nothing can rule on becomes the bracket's UPPER
                # UNCERTAINTY: the search carries on below it, so the run still
                # reports the highest F that actually reached equilibrium, but that
                # edge is never counted as a measured failure and never sets the
                # answer by being averaged into it.
                _note_inconclusive(F_mid, solution)
                F_right = F_mid
                iteration += 1
                continue

            if _stable(solution):
                F_left = F_mid
                last_converged_solution = solution
                if debug_level >= 1:
                    print(f"    -> {_verdict_note(solution)} ({solution['iterations']} iters)")
            else:
                F_right = F_mid
                failed_edge_solution = solution
                if debug_level >= 1:
                    print(f"    -> {_verdict_note(solution)} ({solution['iterations']} iters)")

            iteration += 1

    # Report the midpoint of the final bracket (unbiased, +/- tolerance/2, or exactly
    # the grid-cell center when grid bisection is used);
    # the full bracket is returned in 'final_interval'.
    #
    # An INCONCLUSIVE trial does not change how the number is reported: the midpoint
    # is still the estimate, +/- half the bracket. What changes is what the bracket's
    # upper edge means - an undecided trial rather than a measured failure - and that
    # is disclosed in 'inconclusive', in 'note', and on the log.
    critical_FS = 0.5 * (F_left + F_right)

    _ssrm_progress(progress_callback, _total() * SUBDIV, _total() * SUBDIV, f"FS = {critical_FS:.3f}")
    if debug_level >= 1:
        print(f"\n  SSRM result: FS = {critical_FS:.4f}")
        print(f"  Final interval: [{F_left:.4f}, {F_right:.4f}]")

    return {
        "converged": True,
        "FS": critical_FS,
        "last_solution": last_converged_solution,
        "iterations_ssrm": iteration,
        "final_interval": (F_left, F_right),
        "interval_width": F_right - F_left,
        # The post-peak set the failed-edge trial shed to (None when nothing can
        # soften): what the at-failure capture starts from, so a softening bar
        # is shown at its residual in the field that failed.
        "failed_edge_softened": _failed_edge_softened(failed_edge_solution),
        # Per-trial record: F, role, converged/stable, hybrid verdict, u_ratio,
        # growth, exit_reason, iterations. Populated on every criterion so an A/B
        # between criteria costs no extra solves.
        "trials": trials,
        # Trials the iteration ceiling could not decide (empty on a clean run). When
        # this is non-empty the reported FS is still the bracket midpoint, but the
        # bracket's upper edge is an UNCERTAINTY rather than a measured failure, and
        # `note` says so in words.
        "inconclusive": inconclusive,
        "note": (inconclusive[-1]["message"] if inconclusive else None),
        "failure_criterion": ("hybrid" if hybrid else
                              "non_convergence" if max_disp_factor is None
                              else "displacement_limit"),
        # The bisection-on-non-convergence framing is Griffiths & Lane's; the
        # equilibrium test that decides "converged" is Dawson, Roth & Drescher's
        # per-node out-of-balance force. Credit both — G&L's own criterion is the
        # displacement test plus an iteration ceiling, and the displacement test on
        # its own does not discriminate here.
        "method": ("SSRM — Hybrid (non-convergence corroborated by displacement evidence)"
                   if hybrid else
                   "SSRM — Non-Convergence (Griffiths & Lane 1999; equilibrium test "
                   "after Dawson, Roth & Drescher 1999)"
                   if max_disp_factor is None else "SSRM — Displacement Limit")
    }


def _ssrm_displacement_increase(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05, force_tol=1e-3,
                                 oob_window=10, k0=None,
                                 debug_level=0, max_iterations=500,
                                 convergence_tol=1e-3, n_sweep=10,
                                 tension_cutoff=False, char_point=None,
                 dt_scale=1.0, cancel_check=None, progress_callback=None,
                 min_slip_depth=None, ssr_exclude_mask=None,
                 tension_cap_by_elem=None, tension_srf=False, elastic_mask=None,
                 suction_phi_b=None, suction_cap=None, early_failure=False,
                 _prepared=None, _init_state=None):
    # char_point (x, y): when given, the displacement measure is the
    # CHARACTERISTIC-POINT displacement (nearest node) instead of the global
    # maximum — robust when localized background creep away from the
    # mechanism contaminates the global maximum.
    """
    SSRM using Sun et al. (2021) displacement catastrophe method.

    Phase 1: Coarse sweep of n_sweep F values from F_min to F_max, recording
             max displacement at each. A generous displacement cap (50% mesh
             height) is used purely for early termination of hopeless cases.
    Phase 2: Find the interval with the largest displacement ratio (catastrophe).
    Phase 3: Refine within that interval by bisection, tracking displacement
             at each midpoint to determine which half contains the jump.
    """

    # Compute mesh height for early-termination cap
    nodes = fem_data['nodes']
    mesh_height = float(np.max(nodes[:, 1]) - np.min(nodes[:, 1]))
    # Generous cap just to avoid wasting iterations on hopeless cases
    early_term_factor = 0.5

    if debug_level >= 1:
        print("=== SSRM Analysis (Displacement Catastrophe — Sun et al. 2021) ===")
        print(f"  Sweep range: [{F_min:.2f}, {F_max:.2f}], tolerance: {tolerance}")
        print(f"  Coarse sweep points: {n_sweep}")
        print(f"  Early termination cap: {early_term_factor:.0%} of mesh height ({early_term_factor * mesh_height:.1f})")

    dof_offset_local = fem_data.get("dof_offset", None)

    # Characteristic-point holder. If char_point is given, resolve it now;
    # otherwise it is selected AUTOMATICALLY after the coarse sweep (the node
    # whose plastic displacement grows fastest across the sweep — i.e. a node
    # on the failure mechanism). Background zones that creep steadily at all F
    # have a low growth ratio, so the selection lands on the mechanism.
    cp_holder = [None]
    if char_point is not None:
        cp_holder[0] = int(np.argmin(np.hypot(nodes[:, 0] - char_point[0],
                                              nodes[:, 1] - char_point[1])))
        if debug_level >= 1:
            print(f"  Characteristic point (user): node {cp_holder[0]} at "
                  f"({nodes[cp_holder[0], 0]:.2f}, {nodes[cp_holder[0], 1]:.2f})")

    def _node_vp_disp(sol):
        """Per-node VP displacement magnitudes for a solution."""
        u_vp = sol["displacements"] - sol["displacements_elastic"]
        if dof_offset_local is not None:
            _n = len(nodes)
            vp_x = np.array([u_vp[dof_offset_local[nd]] for nd in range(_n)])
            vp_y = np.array([u_vp[dof_offset_local[nd] + 1] for nd in range(_n)])
        else:
            vp_x = u_vp[0::2]
            vp_y = u_vp[1::2]
        return np.sqrt(vp_x**2 + vp_y**2)

    def _measure(sol):
        d = _node_vp_disp(sol)
        return float(d[cp_holder[0]]) if cp_holder[0] is not None else float(d.max())

    def _get_max_vp_disp(F_val, progress_cb=None):
        """Run solve_fem; return the displacement measure (characteristic-point
        VP displacement, or global max before the point is selected)."""
        sol = solve_fem(fem_data, F=F_val, debug_level=max(0, debug_level-1), dt_scale=dt_scale, force_tol=force_tol,
                        oob_window=oob_window, min_slip_depth=min_slip_depth,
                        max_iterations=max_iterations, tolerance=convergence_tol,
                        max_disp_factor=early_term_factor,
                        ssr_exclude_mask=ssr_exclude_mask,
                        tension_cap_by_elem=tension_cap_by_elem,
                        tension_srf=tension_srf, elastic_mask=elastic_mask,
                        suction_phi_b=suction_phi_b, suction_cap=suction_cap,
                        tension_cutoff=tension_cutoff, progress_callback=progress_cb,
                        k0=k0, early_failure=early_failure,
                        _prepared=_prepared, _init_state=_init_state)
        # NOTE on early_failure: this criterion does not read a VERDICT off each
        # trial, it reads a MEASUREMENT — how far the characteristic point moves at
        # this F — and takes the factor of safety from where that measurement jumps.
        # Stopping a running-away trial early shortens the very displacement the
        # curve is made of, so the default here is OFF (set by the caller below);
        # the early-failure thresholds were measured against bisection verdicts,
        # not against this curve.
        # Use VP displacement (total - elastic) to isolate plastic deformation.
        # The elastic component is roughly constant regardless of F and masks
        # the catastrophic growth in plastic displacement at failure.
        return _measure(sol), sol

    # Progress reported as: the coarse sweep + the refining bisection (estimated
    # from one sweep interval's width).
    n_refine = _ssrm_bisect_steps((F_max - F_min) / max(1, n_sweep - 1), tolerance)
    total = n_sweep + n_refine
    SUBDIV = 100
    total_fine = total * SUBDIV

    def _fem_progress(step, prefix):
        """Map a solve_fem inner [0, 1] fraction into ``step``'s slice of the bar."""
        if progress_callback is None:
            return None
        step = min(step, total - 1)

        def _cb(frac, info):
            pos = (step + max(0.0, min(1.0, frac))) * SUBDIV
            _ssrm_progress(progress_callback, min(int(pos), total_fine - 1),
                           total_fine, f"{prefix} · {info}")
        return _cb

    # Phase 1: Coarse sweep
    F_values = np.linspace(F_min, F_max, n_sweep)
    displacements = []
    solutions = []

    if debug_level >= 1:
        print(f"\n  Phase 1: Coarse sweep ({n_sweep} points)")

    for i, F_val in enumerate(F_values):
        from .search import _check_cancel
        _check_cancel(cancel_check)
        _ssrm_progress(progress_callback, i * SUBDIV, total_fine,
                       f"Sweep F={F_val:.3f} ({i + 1}/{n_sweep})")
        max_vp, sol = _get_max_vp_disp(
            F_val, _fem_progress(i, f"Sweep F={F_val:.3f} ({i + 1}/{n_sweep})"))
        displacements.append(max_vp)
        solutions.append(sol)
        if debug_level >= 1:
            print(f"    F = {F_val:.4f}  max_vp_disp = {max_vp:.2f}  "
                  f"(iters={sol['iterations']}, converged={sol['converged']})")

    # Automatic characteristic-point selection: pick the node whose plastic
    # displacement grew the most (ratio) from the first to the last sweep
    # point — a mechanism node — then re-read the sweep curve at that node.
    # Nodes essentially elastic at the high end are excluded via a floor.
    if cp_holder[0] is None and len(solutions) >= 2:
        d_lo = _node_vp_disp(solutions[0])
        d_hi = _node_vp_disp(solutions[-1])
        floor = 1e-4 * mesh_height
        ratio = d_hi / np.maximum(d_lo, floor)
        ratio[d_hi < 10 * floor] = 0.0
        cp_holder[0] = int(np.argmax(ratio))
        if debug_level >= 1:
            print(f"  Characteristic point (auto): node {cp_holder[0]} at "
                  f"({nodes[cp_holder[0], 0]:.2f}, {nodes[cp_holder[0], 1]:.2f}), "
                  f"growth ratio {ratio[cp_holder[0]]:.1f}x")
        displacements = [_measure(sol) for sol in solutions]

    # Phase 2: Find catastrophe interval (largest displacement ratio).
    # Noise floor: ratios computed from a near-zero (essentially elastic)
    # baseline are meaningless — e.g. VP displacement growing from 1e-3 to
    # 2e-2 is a 20x "ratio" but both values are noise relative to the mesh
    # scale. Only intervals whose left endpoint exceeds the floor count.
    disp_floor = 1e-4 * mesh_height
    max_ratio = 0.0
    catastrophe_idx = None
    for i in range(1, len(displacements)):
        if displacements[i-1] <= disp_floor:
            continue
        ratio = displacements[i] / displacements[i-1]
        if ratio > max_ratio:
            max_ratio = ratio
            catastrophe_idx = i
    if catastrophe_idx is None:
        # all baselines below the floor: fall back to the largest absolute jump
        jumps = [displacements[i] - displacements[i-1]
                 for i in range(1, len(displacements))]
        catastrophe_idx = int(np.argmax(jumps)) + 1
        max_ratio = float('nan')

    F_left = F_values[catastrophe_idx - 1]
    F_right = F_values[catastrophe_idx]
    d_left = displacements[catastrophe_idx - 1]
    d_right = displacements[catastrophe_idx]
    last_converged_solution = solutions[catastrophe_idx - 1]

    if debug_level >= 1:
        print(f"\n  Phase 2: Catastrophe detected between F = {F_left:.4f} and {F_right:.4f}")
        print(f"    VP displacement ratio: {max_ratio:.1f}x  ({d_left:.2f} -> {d_right:.2f})")

    # Phase 3: Refine within catastrophe interval by bisection
    iteration = 0
    while (F_right - F_left) > tolerance and iteration < 50:
        from .search import _check_cancel
        _check_cancel(cancel_check)
        F_mid = (F_left + F_right) / 2.0
        step = min(n_sweep + iteration, total - 1)
        _ssrm_progress(progress_callback, step * SUBDIV, total_fine,
                       f"Refine F={F_mid:.3f} [{F_left:.3f}, {F_right:.3f}]")
        d_mid, sol_mid = _get_max_vp_disp(
            F_mid, _fem_progress(step, f"Refine F={F_mid:.3f} [{F_left:.3f}, {F_right:.3f}]"))

        # Compute ratios for each half
        ratio_left_half = d_mid / d_left if d_left > 1e-12 else d_mid
        ratio_right_half = d_right / d_mid if d_mid > 1e-12 else d_right

        if debug_level >= 1:
            print(f"\n  Refine step {iteration+1}: F = {F_mid:.4f}  [{F_left:.4f}, {F_right:.4f}]")
            print(f"    vp_left={d_left:.2f}, vp_mid={d_mid:.2f}, vp_right={d_right:.2f}")
            print(f"    ratio_left_half={ratio_left_half:.2f}x, ratio_right_half={ratio_right_half:.2f}x")

        if ratio_left_half >= ratio_right_half:
            # Catastrophe is in left half
            F_right = F_mid
            d_right = d_mid
        else:
            # Catastrophe is in right half
            F_left = F_mid
            d_left = d_mid
            last_converged_solution = sol_mid

        iteration += 1

    critical_FS = (F_left + F_right) / 2.0

    _ssrm_progress(progress_callback, total_fine, total_fine, f"FS = {critical_FS:.3f}")
    if debug_level >= 1:
        print(f"\n  SSRM result: FS = {critical_FS:.4f}")
        print(f"  Final interval: [{F_left:.4f}, {F_right:.4f}]")

    return {
        "converged": True,
        "FS": critical_FS,
        "last_solution": last_converged_solution,
        "iterations_ssrm": iteration,
        "final_interval": (F_left, F_right),
        "interval_width": F_right - F_left,
        "method": "SSRM — Displacement Catastrophe (Sun et al. 2021)",
        "sweep_F": F_values.tolist(),
        "sweep_vp_displacements": displacements,
    }


def _coo_to_csr_ordered(rows, cols, vals, shape):
    """Sum stacked COO triplets into a CSR matrix, in the order they were staged.

    ``scipy``'s own duplicate summation is free to add the entries of a repeated
    (row, col) in any order, and a floating-point sum is not associative: the
    last bits of the assembled stiffness would depend on how the library chose
    to sort. The stiffness feeds a factorization that feeds a viscoplastic
    fixed-point iteration with a convergence test, and near-critical strength
    reduction trials are decided in exactly those last bits, so the summation
    order is held fixed here instead.

    Entries are grouped by (row, col) with a stable sort -- so within a group
    they keep the element order they were staged in -- and then added one rank
    at a time, which is the same left-to-right accumulation the entry-by-entry
    ``K[i, j] += ...`` fill performed. Exact zeros are dropped, as the dense
    round-trip this replaced also dropped them.
    """
    n_dof = shape[0]
    if not rows:
        return csr_matrix(shape)

    rows = np.concatenate(rows).astype(np.int64, copy=False)
    cols = np.concatenate(cols).astype(np.int64, copy=False)
    vals = np.concatenate(vals)
    if len(rows) == 0:
        return csr_matrix(shape)

    # Stable sort on the flattened (row, col) key: duplicates stay in staging order.
    key = rows * np.int64(shape[1]) + cols
    order = np.argsort(key, kind='stable')
    key = key[order]
    vals = vals[order]

    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    group_of = np.zeros(len(key), dtype=np.int64)
    group_of[starts[1:]] = 1
    np.cumsum(group_of, out=group_of)
    rank = np.arange(len(key), dtype=np.int64) - starts[group_of]

    data = vals[starts].copy()                     # rank 0 of every group
    for r in range(1, int(rank.max()) + 1 if len(rank) else 1):
        sel = rank == r
        if not sel.any():
            break
        data[group_of[sel]] += vals[sel]           # one entry per group per rank

    g_rows = key[starts] // shape[1]
    g_cols = key[starts] - g_rows * shape[1]

    keep = data != 0.0
    if not keep.all():
        data = data[keep]
        g_rows = g_rows[keep]
        g_cols = g_cols[keep]

    indptr = np.zeros(n_dof + 1, dtype=np.int64)
    np.cumsum(np.bincount(g_rows, minlength=n_dof), out=indptr[1:])
    K = csr_matrix((data, g_cols, indptr), shape=shape)
    K.has_sorted_indices = True
    return K


def build_global_stiffness(nodes, elements, element_types, element_materials, E_by_mat, nu_by_mat, fem_data=None):
    """
    Build global stiffness matrix from 2D soil elements and (optionally) 1D truss + pile beam elements.

    Uses dof_offset from fem_data to support mixed DOF systems (pile nodes have 3 DOFs).

    Assembly is COO: each element's stiffness block and its global row/column
    indices are stacked, and the whole matrix is built in one pass with
    ``_coo_to_csr_ordered``. That helper sums each repeated (row, col) entry in
    ELEMENT ORDER, which is the order the earlier entry-by-entry ``lil_matrix``
    fill used, so the assembled matrix is bit-for-bit what it was -- the change
    is the cost of building it, never its value.
    """
    n_nodes = len(nodes)

    # Get DOF offset map from fem_data if available
    dof_offset = fem_data.get("dof_offset", None) if fem_data is not None else None
    if dof_offset is not None:
        n_dof = int(dof_offset[n_nodes])
    else:
        n_dof = 2 * n_nodes

    coo_rows = []
    coo_cols = []
    coo_vals = []

    def _stage(K_block, dof_idx):
        """Queue one element block (n x n) against its n global DOF indices."""
        n = len(dof_idx)
        coo_rows.append(np.repeat(dof_idx, n))
        coo_cols.append(np.tile(dof_idx, n))
        coo_vals.append(np.asarray(K_block, dtype=float).ravel())

    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1

        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]

        # Get element coordinates
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]

        # Build element stiffness matrix using corrected implementation
        try:
            if elem_type == 3:
                K_elem = build_triangle_stiffness_corrected(elem_coords, E, nu)
            elif elem_type == 6:
                K_elem = build_tri6_stiffness(elem_coords, E, nu)
            elif elem_type == 4:
                K_elem = build_quad4_stiffness(elem_coords, E, nu)
            elif elem_type == 8:
                K_elem = build_quad8_stiffness_reduced_integration_corrected(elem_coords, E, nu)
            elif elem_type == 9:
                K_elem = build_quad9_stiffness(elem_coords, E, nu)
            else:
                print(f"Warning: Element type {elem_type} not supported")
                continue
        except Exception as e:
            print(f"Error building stiffness for element {elem_idx}, type {elem_type}: {e}")
            continue

        # Assemble into global matrix using dof_offset for global DOF indices.
        # An element matrix smaller than its 2 x n_node DOF count contributes
        # only its leading block (the guard the entry-by-entry fill applied).
        dof_idx = _elem_dof_indices(elem_nodes, dof_offset=dof_offset)
        n_use = min(len(dof_idx), K_elem.shape[0], K_elem.shape[1])
        _stage(K_elem[:n_use, :n_use], dof_idx[:n_use])

    # Assemble 1D truss element stiffness matrices (reinforcement only — skip pile elements)
    if fem_data is not None:
        K_global_1d_elems = fem_data.get("K_global_1d_elems", [])
        dof_indices_1d = fem_data.get("dof_indices_1d", np.zeros((0, 6), dtype=int))
        pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(K_global_1d_elems), dtype=bool))

        for elem_idx_1d in range(len(K_global_1d_elems)):
            if pile_elem_mask[elem_idx_1d]:
                continue  # pile elements use beam stiffness, assembled below
            # A two-node bar reaches four DOFs and a three-node bar six, so the
            # element's own block size selects the DOFs it is staged against.
            n_use = np.asarray(K_global_1d_elems[elem_idx_1d]).shape[0]
            _stage(K_global_1d_elems[elem_idx_1d],
                   np.asarray(dof_indices_1d[elem_idx_1d])[:n_use])

        # Assemble pile beam element stiffness matrices (Euler-Bernoulli, 6x6 on
        # a two-node element and 9x9 on a three-node one)
        K_global_pile_elems = fem_data.get("K_global_pile_elems", [])
        dof_indices_pile = fem_data.get("dof_indices_pile", np.zeros((0, 9), dtype=int))

        for p_idx in range(len(K_global_pile_elems)):
            n_use = np.asarray(K_global_pile_elems[p_idx]).shape[0]
            _stage(K_global_pile_elems[p_idx],
                   np.asarray(dof_indices_pile[p_idx])[:n_use])

    return _coo_to_csr_ordered(coo_rows, coo_cols, coo_vals, (n_dof, n_dof))



def _no_piezo_line_error(names):
    """u='piezo' with no piezometric line in the file at all.

    The same silent-zero class as sampling one from outside its extent: the model
    asks for pore pressure from a line, the file defines none, and every point
    reads zero. There is no reading of "piezometric line, but there isn't one"
    that means a dry slope -- a dry slope is u='none' -- so this is refused rather
    than defaulted.
    """
    who = ("material(s) " + ", ".join(f"'{n}'" for n in names) if names
           else "the materials of this model")
    return ValueError(
        f"Pore pressure: {who} take pore pressure from a piezometric "
        f"line (u='piezo'), but the file defines no piezometric line -- the piezo "
        f"sheet carries fewer than two points. Every point would read zero pore "
        f"pressure -- an unsafe, non-conservative result. Draw the piezometric "
        f"line, or set u='none' for a dry model (or 'seep' / 'ru' for the other "
        f"pore-pressure sources).")


def _piezo_extent_tol(px):
    """Floating-point slack on a piezometric line's x-extent. A point sitting ON
    the line's own end is inside it, and the corpus is full of lines drawn exactly
    to the domain edge, so the comparison cannot be a bare `<`. A real violation is
    orders of magnitude larger than this."""
    return 1e-9 * max(1.0, float(px[-1] - px[0]))


def _check_piezo_extent(xq, px, py, what):
    """Every point that samples a piezometric line must lie inside it.

    ``xq`` may be a scalar or an array of x-coordinates that are about to read a
    pore pressure off the piezometric line whose ascending x-coordinates are
    ``px``. Outside that line's own x-extent there is no piezometric surface, and
    the interpolation returns NaN -> zero pore pressure. Silently reading zero
    below a water table over-predicts the factor of safety, so it is raised here
    instead, naming the offending point and the extent it fell out of.

    A line that legitimately stops short of the domain (an upstream pool with no
    downstream pond) is fine as long as nothing samples it beyond its end. Where a
    model really is dry past the line, say so explicitly -- carry the line on at an
    elevation below the mesh -- rather than leaning on a silent default.

    ``what`` names the sampling layer for the message (e.g. "mesh node").
    """
    xq = np.atleast_1d(np.asarray(xq, dtype=float))
    tol = _piezo_extent_tol(px)
    outside = (xq < px[0] - tol) | (xq > px[-1] + tol)
    if not np.any(outside):
        return
    bad = np.flatnonzero(outside)
    first = int(bad[0])
    where = f"{what} {first}" if xq.size > 1 else what
    tally = (f" ({bad.size} of {xq.size} sampled points are outside, x from "
             f"{xq[outside].min():.6g} to {xq[outside].max():.6g})"
             if xq.size > 1 else "")
    raise ValueError(
        f"Pore pressure: {where} is at x = {xq[first]:.6g}, outside the "
        f"piezometric line's x-extent [{px[0]:.6g}, {px[-1]:.6g}]"
        f"{tally}. The materials take "
        f"pore pressure from that line (u='piezo'), so there is no piezometric "
        f"surface to read and the pore pressure would silently be zero -- an "
        f"unsafe, non-conservative result. Extend the piezometric line across the "
        f"mesh (drop it below the mesh where the model is dry), or trim the mesh "
        f"to the line.")


def _piezo_cos2(xq, px, py):
    """Phreatic-inclination (Hu) factor cos^2(theta) of the LOCAL piezo-line
    segment at x = xq (XSTABL / Slide "Hu: auto"), mirroring the LEM's _cos2
    in slice.py: 1/(1 + m^2) inside a segment, 1.0 outside the line's x-range.
    px must be ascending; xq may be scalar or array."""
    xq = np.asarray(xq, dtype=float)
    k = np.clip(np.searchsorted(px, xq, side='right') - 1, 0, len(px) - 2)
    dx = px[k + 1] - px[k]
    m = np.where(dx > 0, (py[k + 1] - py[k]) / np.where(dx > 0, dx, 1.0), 0.0)
    cos2 = 1.0 / (1.0 + m * m)
    outside = (xq < px[0]) | (xq > px[-1])
    return np.where(outside, 1.0, cos2)


def build_gravity_loads(nodes, elements, element_types, element_materials, gamma_by_mat, k_seismic, fem_data=None):
    """
    Build gravity load vector using Griffiths & Lane (1999) approach.

    Uses equation 3 from the paper: p(e) = gamma * integral[Ve] N^T d(vol)
    This integrates shape functions over each element to properly distribute gravity loads.

    Uses dof_offset from fem_data when available to support mixed DOF systems.
    """
    n_nodes = len(nodes)

    # Get DOF offset map from fem_data if available
    dof_offset = fem_data.get("dof_offset", None) if fem_data is not None else None
    if dof_offset is not None:
        n_dof = int(dof_offset[n_nodes])
    else:
        n_dof = 2 * n_nodes

    F_gravity = np.zeros(n_dof)

    def _node_dof_x(node):
        return dof_offset[node] if dof_offset is not None else 2 * node

    def _node_dof_y(node):
        return (dof_offset[node] + 1) if dof_offset is not None else 2 * node + 1

    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        gamma = gamma_by_mat[mat_id]

        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]

        if elem_type == 3:  # 3-node triangle
            # For linear triangles, shape function integration gives equal distribution (1/3 each)
            x1, y1 = elem_coords[0]
            x2, y2 = elem_coords[1]
            x3, y3 = elem_coords[2]
            area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

            # Each node gets 1/3 of the element weight (exact for linear shape functions)
            for i, node in enumerate(elem_nodes):
                load = gamma * area / 3.0
                F_gravity[_node_dof_y(node)] -= load  # Vertical (negative = downward)
                F_gravity[_node_dof_x(node)] += k_seismic * load  # Horizontal seismic

        elif elem_type == 8:  # 8-node quad
            # For 8-node quads, use 2x2 Gauss integration as in Griffiths
            # This properly weights corner vs midside nodes

            # Gauss points for 2x2 integration
            gauss_coord = 1.0 / np.sqrt(3.0)
            xi_points = np.array([-gauss_coord, gauss_coord])
            eta_points = np.array([-gauss_coord, gauss_coord])
            weights = np.array([1.0, 1.0])

            # Initialize element load vector (local, 2 DOFs per node)
            elem_loads = np.zeros(2 * elem_type)

            # Numerical integration over Gauss points
            for i in range(2):
                for j in range(2):
                    xi = xi_points[i]
                    eta = eta_points[j]
                    w = weights[i] * weights[j]

                    # Shape functions for 8-node quad at (xi, eta)
                    N = compute_quad8_shape_functions(xi, eta)

                    # Jacobian for coordinate transformation
                    J = compute_quad8_jacobian(elem_coords, xi, eta)
                    det_J = np.linalg.det(J)

                    # Accumulate load contribution: w * det(J) * gamma * N
                    for k in range(8):
                        elem_loads[2*k + 1] -= w * det_J * gamma * N[k]  # Vertical
                        elem_loads[2*k] += w * det_J * gamma * k_seismic * N[k]  # Horizontal

            # Add element loads to global vector
            for i, node in enumerate(elem_nodes):
                F_gravity[_node_dof_x(node)] += elem_loads[2*i]
                F_gravity[_node_dof_y(node)] += elem_loads[2*i + 1]

        elif elem_type == 6:  # 6-node triangle
            gauss_pts_tri, gauss_wts_tri = get_gauss_points_tri3()
            elem_loads = np.zeros(2 * 6)
            for gp_idx in range(3):
                L1, L2, L3 = gauss_pts_tri[gp_idx]
                w = gauss_wts_tri[gp_idx]
                N = compute_tri6_shape_functions(L1, L2, L3)
                x0, y0 = elem_coords[0]
                x1, y1 = elem_coords[1]
                x2, y2 = elem_coords[2]
                det_J = (x0 - x2) * (y1 - y2) - (x1 - x2) * (y0 - y2)
                integration_weight = 0.5 * abs(det_J) * w
                for k in range(6):
                    elem_loads[2*k + 1] -= integration_weight * gamma * N[k]
                    elem_loads[2*k] += integration_weight * gamma * k_seismic * N[k]
            for i, node in enumerate(elem_nodes):
                F_gravity[_node_dof_x(node)] += elem_loads[2*i]
                F_gravity[_node_dof_y(node)] += elem_loads[2*i + 1]

        elif elem_type == 4:  # 4-node quad
            gauss_pts, gauss_wts = get_gauss_points_2x2()
            elem_loads = np.zeros(2 * 4)
            for gp_idx in range(4):
                xi, eta = gauss_pts[gp_idx]
                w = gauss_wts[gp_idx]
                N = compute_quad4_shape_functions(xi, eta)
                _, det_J = _compute_B_and_detJ_quad4(elem_coords, xi, eta)
                for k in range(4):
                    elem_loads[2*k + 1] -= w * abs(det_J) * gamma * N[k]
                    elem_loads[2*k] += w * abs(det_J) * gamma * k_seismic * N[k]
            for i, node in enumerate(elem_nodes):
                F_gravity[_node_dof_x(node)] += elem_loads[2*i]
                F_gravity[_node_dof_y(node)] += elem_loads[2*i + 1]

        elif elem_type == 9:  # 9-node quad
            gauss_pts, gauss_wts = get_gauss_points_3x3()
            elem_loads = np.zeros(2 * 9)
            for gp_idx in range(9):
                xi, eta = gauss_pts[gp_idx]
                w = gauss_wts[gp_idx]
                N = compute_quad9_shape_functions(xi, eta)
                _, det_J = _compute_B_and_detJ_quad9(elem_coords, xi, eta)
                for k in range(9):
                    elem_loads[2*k + 1] -= w * abs(det_J) * gamma * N[k]
                    elem_loads[2*k] += w * abs(det_J) * gamma * k_seismic * N[k]
            for i, node in enumerate(elem_nodes):
                F_gravity[_node_dof_x(node)] += elem_loads[2*i]
                F_gravity[_node_dof_y(node)] += elem_loads[2*i + 1]

    return F_gravity


def compute_quad8_shape_functions(xi, eta):
    """
    Compute shape functions for 8-node serendipity quadrilateral at (xi, eta).
    
    Node numbering:
    3---6---2
    |       |
    7       5
    |       |
    0---4---1
    """
    N = np.zeros(8)
    
    # Corner nodes
    N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1)
    N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1)
    N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1)
    N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1)
    
    # Midside nodes
    N[4] = 0.5 * (1 - xi**2) * (1 - eta)
    N[5] = 0.5 * (1 + xi) * (1 - eta**2)
    N[6] = 0.5 * (1 - xi**2) * (1 + eta)
    N[7] = 0.5 * (1 - xi) * (1 - eta**2)
    
    return N


def compute_quad8_jacobian(coords, xi, eta):
    """
    Compute Jacobian matrix for 8-node quad at (xi, eta).
    """
    # Shape function derivatives
    dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)
    
    # Jacobian matrix
    J = np.zeros((2, 2))
    for i in range(8):
        J[0, 0] += dN_dxi[i] * coords[i, 0]   # dx/dxi
        J[0, 1] += dN_dxi[i] * coords[i, 1]   # dy/dxi
        J[1, 0] += dN_deta[i] * coords[i, 0]  # dx/deta
        J[1, 1] += dN_deta[i] * coords[i, 1]  # dy/deta
    
    return J




def compute_B_matrix_triangle(coords):
    """Compute B matrix and area for triangle element."""
    
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    if area < 1e-12:
        return np.zeros((3, 6)), 0.0
    
    # Shape function derivatives
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    # B matrix (standard linear triangle)
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],  # εx = ∂u/∂x
        [0,  c1, 0,  c2, 0,  c3],  # εy = ∂v/∂y
        [c1, b1, c2, b2, c3, b3]   # γxy = ∂u/∂y + ∂v/∂x
    ]) / (2 * area)
    
    return B, area


def get_gauss_points_2x2():
    """Get 2x2 Gauss quadrature points and weights for reduced integration."""
    # 2x2 Gauss points in natural coordinates
    gp = 1.0 / np.sqrt(3.0)
    gauss_points = [
        (-gp, -gp),  # Gauss point 0
        ( gp, -gp),  # Gauss point 1
        ( gp,  gp),  # Gauss point 2
        (-gp,  gp),  # Gauss point 3
    ]
    weights = [1.0, 1.0, 1.0, 1.0]  # Equal weights for 2x2
    return gauss_points, weights


def get_gauss_points_tri3():
    """Get 3-point Gauss quadrature for triangles (area coordinates).
    Integration: integral = 0.5 * |detJ| * sum(w_i * f_i)."""
    gauss_points = [
        (1.0/6.0, 1.0/6.0, 2.0/3.0),
        (1.0/6.0, 2.0/3.0, 1.0/6.0),
        (2.0/3.0, 1.0/6.0, 1.0/6.0),
    ]
    weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
    return gauss_points, weights


def get_gauss_points_3x3():
    """Get 3x3 Gauss quadrature points and weights for full integration (quad9)."""
    pts_1d = [-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)]
    wts_1d = [5.0/9.0, 8.0/9.0, 5.0/9.0]
    gauss_points = []
    weights = []
    for i in range(3):
        for j in range(3):
            gauss_points.append((pts_1d[i], pts_1d[j]))
            weights.append(wts_1d[i] * wts_1d[j])
    return gauss_points, weights


def compute_tri6_shape_functions(L1, L2, L3):
    """Shape functions for 6-node quadratic triangle at area coordinates (L1, L2, L3).
    Node ordering: corners (0,1,2), edge midpoints (3=edge 0-1, 4=edge 1-2, 5=edge 2-0)."""
    return np.array([
        L1 * (2*L1 - 1),
        L2 * (2*L2 - 1),
        L3 * (2*L3 - 1),
        4*L1*L2,
        4*L2*L3,
        4*L3*L1,
    ])


def compute_tri6_shape_derivatives(L1, L2, L3):
    """Shape function derivatives for 6-node triangle w.r.t. area coordinates.
    Returns dN_dL1(6,), dN_dL2(6,), dN_dL3(6,)."""
    dN_dL1 = np.array([4*L1 - 1, 0, 0, 4*L2, 0, 4*L3])
    dN_dL2 = np.array([0, 4*L2 - 1, 0, 4*L1, 4*L3, 0])
    dN_dL3 = np.array([0, 0, 4*L3 - 1, 0, 4*L2, 4*L1])
    return dN_dL1, dN_dL2, dN_dL3


def compute_quad4_shape_functions(xi, eta):
    """Shape functions for 4-node bilinear quadrilateral at (xi, eta)."""
    return 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])


def compute_quad4_shape_derivatives(xi, eta):
    """Shape function derivatives for 4-node quad. Returns dN_dxi(4,), dN_deta(4,)."""
    dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return dN_dxi, dN_deta


def compute_quad9_shape_functions(xi, eta):
    """Shape functions for 9-node Lagrange quadrilateral at (xi, eta).
    Node ordering: corners (0-3), edge midpoints (4-7), center (8)."""
    return np.array([
        0.25 * xi * (xi - 1) * eta * (eta - 1),
        0.25 * xi * (xi + 1) * eta * (eta - 1),
        0.25 * xi * (xi + 1) * eta * (eta + 1),
        0.25 * xi * (xi - 1) * eta * (eta + 1),
        0.5 * (1 - xi*xi) * eta * (eta - 1),
        0.5 * xi * (xi + 1) * (1 - eta*eta),
        0.5 * (1 - xi*xi) * eta * (eta + 1),
        0.5 * xi * (xi - 1) * (1 - eta*eta),
        (1 - xi*xi) * (1 - eta*eta),
    ])


def compute_quad9_shape_derivatives(xi, eta):
    """Shape function derivatives for 9-node quad. Returns dN_dxi(9,), dN_deta(9,)."""
    dN_dxi = np.array([
        0.25 * (2*xi - 1) * eta * (eta - 1),
        0.25 * (2*xi + 1) * eta * (eta - 1),
        0.25 * (2*xi + 1) * eta * (eta + 1),
        0.25 * (2*xi - 1) * eta * (eta + 1),
        -xi * eta * (eta - 1),
        0.5 * (2*xi + 1) * (1 - eta*eta),
        -xi * eta * (eta + 1),
        0.5 * (2*xi - 1) * (1 - eta*eta),
        -2*xi * (1 - eta*eta),
    ])
    dN_deta = np.array([
        0.25 * xi * (xi - 1) * (2*eta - 1),
        0.25 * xi * (xi + 1) * (2*eta - 1),
        0.25 * xi * (xi + 1) * (2*eta + 1),
        0.25 * xi * (xi - 1) * (2*eta + 1),
        0.5 * (1 - xi*xi) * (2*eta - 1),
        -xi * (xi + 1) * eta,
        0.5 * (1 - xi*xi) * (2*eta + 1),
        -xi * (xi - 1) * eta,
        -2*eta * (1 - xi*xi),
    ])
    return dN_dxi, dN_deta


def compute_triangle_strains_manual(coords, displacements):
    """Manually compute triangle strains from displacements."""
    
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    if area < 1e-12:
        return np.array([0.0, 0.0, 0.0])
    
    # Shape function derivatives
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    # B matrix (standard linear triangle)
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],  # εx = ∂u/∂x
        [0,  c1, 0,  c2, 0,  c3],  # εy = ∂v/∂y
        [c1, b1, c2, b2, c3, b3]   # γxy = ∂u/∂y + ∂v/∂x
    ]) / (2 * area)
    
    # Strains
    strains = B @ displacements
    return strains


def build_constitutive_matrix(E, nu):
    """Build constitutive matrix for plane strain - standard tension-positive convention."""
    # A near-incompressible nu (> 0.45) is reported by preflight
    # (mat.nu_implausible) before any run; nothing is printed per element here.

    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    # Standard tension-positive convention (σ > 0 in tension, σ < 0 in compression)
    return D


def build_constitutive_matrix_4(E, nu):
    """4-component plane-strain D matrix (tension-positive), stress order
    [sig_x, sig_y, tau_xy, sig_z] and strain order [eps_x, eps_y, gamma_xy, eps_z].

    Matches Smith & Griffiths' nst=4 plane-strain formulation (p62.f90): sigma_z
    is carried explicitly so the viscoplastic algorithm can relax it through
    plastic eps_z (total eps_z = 0, elastic eps_z = -eps_z^vp).
    """
    factor = E / ((1 + nu) * (1 - 2 * nu))
    return factor * np.array([
        [1 - nu, nu,     0,              nu    ],
        [nu,     1 - nu, 0,              nu    ],
        [0,      0,      (1 - 2 * nu)/2, 0     ],
        [nu,     nu,     0,              1 - nu],
    ])


def stress_invariants(stress4_tp):
    """Invariants of a 4-component plane-strain stress [sx, sy, txy, sz]
    (tension-positive), after Smith & Griffiths' INVAR.

    Returns:
        sigm  : mean stress (sx+sy+sz)/3
        dsbar : deviatoric stress sqrt(3*J2) (von Mises equivalent)
        theta : Lode angle in radians, in [-pi/6, pi/6];
                sin(3*theta) = -13.5*J3/dsbar^3 (S&G convention)
    """
    sx, sy, txy, sz = stress4_tp
    sigm = (sx + sy + sz) / 3.0
    dsbar = sqrt((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2 + 6.0 * txy**2) / sqrt(2.0)
    if dsbar < 1e-10:
        theta = 0.0
    else:
        dx = (2.0 * sx - sy - sz) / 3.0
        dy = (2.0 * sy - sz - sx) / 3.0
        dz = (2.0 * sz - sx - sy) / 3.0
        xj3 = dx * dy * dz - dz * txy**2
        sine = -13.5 * xj3 / dsbar**3
        sine = min(1.0, max(-1.0, sine))
        theta = asin(sine) / 3.0
    return sigm, dsbar, theta


def mc_yield_invariants(sigm, dsbar, theta, c, phi):
    """Mohr-Coulomb yield function in invariant form (S&G MOCOUF), for
    tension-positive stresses. F > 0 means yielding.

    F = sigm*sin(phi) + dsbar*(cos(theta)/sqrt(3) - sin(theta)*sin(phi)/3)
        - c*cos(phi)
    """
    snph = sin(phi)
    return (sigm * snph
            + dsbar * (cos(theta) / sqrt(3.0) - sin(theta) * snph / 3.0)
            - c * cos(phi))




def compute_strains(nodes, elements, element_types, displacements, dof_offset=None):
    """
    Compute element strains for visualization.

    If dof_offset is provided, uses it for DOF indexing (mixed DOF system with pile nodes).
    """
    n_elements = len(elements)
    strains = np.zeros((n_elements, 4))  # [eps_x, eps_y, gamma_xy, max_shear_strain]

    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]

        # Get element displacements (translational DOFs only)
        elem_disp = np.zeros(2 * elem_type)
        for i, node in enumerate(elem_nodes):
            if dof_offset is not None:
                base = dof_offset[node]
            else:
                base = 2 * node
            elem_disp[2*i] = displacements[base]
            elem_disp[2*i+1] = displacements[base + 1]
        
        # Compute strains at element centroid
        if elem_type == 3:
            element_strains = compute_triangle_strains_manual(elem_coords, elem_disp)
        elif elem_type == 6:
            B, det_J = _compute_B_and_detJ_tri6(elem_coords, 1.0/3.0, 1.0/3.0, 1.0/3.0)
            element_strains = B @ elem_disp
        elif elem_type == 4:
            B, det_J = _compute_B_and_detJ_quad4(elem_coords, 0.0, 0.0)
            element_strains = B @ elem_disp
        elif elem_type == 8:
            xi, eta = 0.0, 0.0
            element_strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta)
        elif elem_type == 9:
            B, det_J = _compute_B_and_detJ_quad9(elem_coords, 0.0, 0.0)
            element_strains = B @ elem_disp
        else:
            element_strains = np.array([0.0, 0.0, 0.0])
        
        eps_x = element_strains[0]
        eps_y = element_strains[1] 
        gamma_xy = element_strains[2]
        
        # Maximum shear strain
        max_shear_strain = sqrt(((eps_x - eps_y) / 2)**2 + (gamma_xy / 2)**2)
        
        strains[elem_idx] = [eps_x, eps_y, gamma_xy, max_shear_strain]
    
    return strains


def compute_quad8_strains_at_xi_eta(coords, displacements, xi, eta):
    """
    Compute strains for 8-node quadrilateral at specific (xi, eta) coordinates.
    """
    # 8-node quad shape function derivatives at (xi, eta)
    dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)
    
    # Jacobian matrix and its inverse
    J = np.zeros((2, 2))
    for i in range(8):
        x, y = coords[i]
        J[0, 0] += dN_dxi[i] * x   # dx/dxi
        J[0, 1] += dN_dxi[i] * y   # dy/dxi
        J[1, 0] += dN_deta[i] * x  # dx/deta
        J[1, 1] += dN_deta[i] * y  # dy/deta
    
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    
    if abs(det_J) < 1e-12:
        return np.array([0.0, 0.0, 0.0])
    
    # Inverse Jacobian
    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J
    
    # Shape function derivatives in physical coordinates
    dN_dx = np.zeros(8)
    dN_dy = np.zeros(8)
    for i in range(8):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    # B matrix for strain calculation (standard tension positive)
    B = np.zeros((3, 16))  # 3 strains x 16 DOFs (8 nodes x 2 DOFs)
    for i in range(8):
        B[0, 2*i] = dN_dx[i]      # εx = ∂u/∂x
        B[1, 2*i+1] = dN_dy[i]    # εy = ∂v/∂y
        B[2, 2*i] = dN_dy[i]      # γxy = ∂u/∂y + ∂v/∂x
        B[2, 2*i+1] = dN_dx[i]    # γxy = ∂u/∂y + ∂v/∂x
    
    # Compute strains
    strains = B @ displacements
    return strains


def compute_quad8_shape_derivatives(xi, eta):
    """
    Compute shape function derivatives for 8-node quadrilateral at (xi, eta).
    
    Uses correct serendipity formulation with CCW node ordering:
    3 --- 6 --- 2
    |           |
    7     +     5
    |           |
    0 --- 4 --- 1
    
    Corner nodes: 0(-1,-1), 1(1,-1), 2(1,1), 3(-1,1) 
    Edge nodes: 4(0,-1), 5(1,0), 6(0,1), 7(-1,0)
    """
    
    # Serendipity shape function derivatives for CCW node ordering
    # (From working implementation in seep.py)
    dN_dxi = np.array([
        -0.25*(1-eta)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta), # Node 0: corner (-1,-1)
        0.25*(1-eta)*(xi-eta-1) + 0.25*(1+xi)*(1-eta),   # Node 1: corner (1,-1)
        0.25*(1+eta)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),   # Node 2: corner (1,1)
        -0.25*(1+eta)*(-xi+eta-1) - 0.25*(1-xi)*(1+eta), # Node 3: corner (-1,1)
        -xi*(1-eta),                                      # Node 4: edge (0,-1)
        0.5*(1-eta*eta),                                  # Node 5: edge (1,0)
        -xi*(1+eta),                                      # Node 6: edge (0,1)
        -0.5*(1-eta*eta)                                  # Node 7: edge (-1,0)
    ])
    
    dN_deta = np.array([
        -0.25*(1-xi)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta),  # Node 0: corner (-1,-1)
        -0.25*(1+xi)*(xi-eta-1) - 0.25*(1+xi)*(1-eta),   # Node 1: corner (1,-1)
        0.25*(1+xi)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),    # Node 2: corner (1,1)
        0.25*(1-xi)*(-xi+eta-1) + 0.25*(1-xi)*(1+eta),   # Node 3: corner (-1,1)
        -0.5*(1-xi*xi),                                   # Node 4: edge (0,-1)
        -eta*(1+xi),                                      # Node 5: edge (1,0)
        0.5*(1-xi*xi),                                    # Node 6: edge (0,1)
        -eta*(1-xi)                                       # Node 7: edge (-1,0)
    ])
    
    return dN_dxi, dN_deta


def build_quad8_stiffness_reduced_integration_corrected(coords, E, nu):
    """
    Build stiffness matrix for 8-node quadrilateral with 2x2 reduced integration.
    
    This follows the Griffiths & Lane (1999) implementation exactly:
    - 8-node serendipity quadrilateral elements
    - 2x2 reduced integration (4 Gauss points) 
    - Prevents volumetric locking in nearly incompressible materials
    """
    # Constitutive matrix for plane strain
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    # 2x2 Gauss points for reduced integration (exactly as in Griffiths paper)
    gauss_coord = 1.0 / np.sqrt(3.0)  # = 0.5773502692
    xi_points = np.array([-gauss_coord, gauss_coord])
    eta_points = np.array([-gauss_coord, gauss_coord])
    weights = np.array([1.0, 1.0, 1.0, 1.0])  # 2D weights = 1 * 1
    
    K = np.zeros((16, 16))  # 8 nodes x 2 DOF = 16x16 matrix
    
    gp_idx = 0
    for i in range(2):
        for j in range(2):
            xi, eta = xi_points[i], eta_points[j] 
            w = weights[gp_idx]
            gp_idx += 1
            
            # Use the existing correct shape function derivatives
            dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)
            
            # Jacobian matrix
            J = np.zeros((2, 2))
            for a in range(8):
                J[0,0] += dN_dxi[a] * coords[a,0]   # dx/dxi
                J[0,1] += dN_dxi[a] * coords[a,1]   # dy/dxi
                J[1,0] += dN_deta[a] * coords[a,0]  # dx/deta
                J[1,1] += dN_deta[a] * coords[a,1]  # dy/deta
            
            det_J = J[0,0] * J[1,1] - J[0,1] * J[1,0]
            
            if abs(det_J) < 1e-12:
                print(f"Warning: Nearly singular Jacobian in quad8 element: det(J) = {det_J}")
                continue
            
            # Inverse Jacobian
            J_inv = np.array([[J[1,1], -J[0,1]], [-J[1,0], J[0,0]]]) / det_J
            
            # Shape function derivatives in physical coordinates
            dN_dx = np.zeros(8)
            dN_dy = np.zeros(8)
            for a in range(8):
                dN_dx[a] = J_inv[0,0]*dN_dxi[a] + J_inv[0,1]*dN_deta[a]
                dN_dy[a] = J_inv[1,0]*dN_dxi[a] + J_inv[1,1]*dN_deta[a]
            
            # B matrix (strain-displacement, standard tension positive)
            B = np.zeros((3, 16))  # 3 strains x 16 DOF
            for a in range(8):
                B[0, 2*a] = dN_dx[a]      # εx = ∂u/∂x
                B[1, 2*a+1] = dN_dy[a]    # εy = ∂v/∂y
                B[2, 2*a] = dN_dy[a]      # γxy = ∂u/∂y + ∂v/∂x
                B[2, 2*a+1] = dN_dx[a]    # γxy = ∂u/∂y + ∂v/∂x
            
            # Element stiffness matrix contribution
            K += w * det_J * (B.T @ D @ B)
    
    return K


def build_triangle_stiffness_corrected(coords, E, nu):
    """
    Build corrected stiffness matrix for triangular element (plane strain).
    """
    x1, y1 = coords[0]
    x2, y2 = coords[1] 
    x3, y3 = coords[2]
    
    # Area
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    if area < 1e-12:
        print(f"Warning: Very small triangle area: {area}")
        return np.zeros((6, 6))
    
    # Shape function derivatives
    b1 = y2 - y3
    b2 = y3 - y1  
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    # B matrix (standard linear triangle)
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],  # εx = ∂u/∂x
        [0,  c1, 0,  c2, 0,  c3],  # εy = ∂v/∂y
        [c1, b1, c2, b2, c3, b3]   # γxy = ∂u/∂y + ∂v/∂x
    ]) / (2 * area)
    
    # Constitutive matrix (plane strain)
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    
    # Element stiffness matrix
    K_elem = area * B.T @ D @ B

    return K_elem


def build_tri6_stiffness(coords, E, nu):
    """Build stiffness matrix (12x12) for 6-node triangle with 3-point integration."""
    D = build_constitutive_matrix(E, nu)
    gauss_points, weights = get_gauss_points_tri3()
    K = np.zeros((12, 12))
    for gp_idx in range(3):
        L1, L2, L3 = gauss_points[gp_idx]
        B, det_J = _compute_B_and_detJ_tri6(coords, L1, L2, L3)
        w = weights[gp_idx] * 0.5 * abs(det_J)  # 0.5 for area coordinate mapping
        K += w * (B.T @ D @ B)
    return K


def build_quad4_stiffness(coords, E, nu):
    """Build stiffness matrix (8x8) for 4-node quad with 2x2 integration."""
    D = build_constitutive_matrix(E, nu)
    gauss_points, weights = get_gauss_points_2x2()
    K = np.zeros((8, 8))
    for gp_idx in range(4):
        xi, eta = gauss_points[gp_idx]
        B, det_J = _compute_B_and_detJ_quad4(coords, xi, eta)
        w = weights[gp_idx] * abs(det_J)
        K += w * (B.T @ D @ B)
    return K


def build_quad9_stiffness(coords, E, nu):
    """Build stiffness matrix (18x18) for 9-node quad with 3x3 integration."""
    D = build_constitutive_matrix(E, nu)
    gauss_points, weights = get_gauss_points_3x3()
    K = np.zeros((18, 18))
    for gp_idx in range(9):
        xi, eta = gauss_points[gp_idx]
        B, det_J = _compute_B_and_detJ_quad9(coords, xi, eta)
        w = weights[gp_idx] * abs(det_J)
        K += w * (B.T @ D @ B)
    return K
