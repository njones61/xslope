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

import logging
import time
import warnings
from collections import defaultdict, deque

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, splu
from shapely.geometry import LineString, Point

from .units import require_gamma_water, infer_system_from_gamma_water

logger = logging.getLogger(__name__)


class SeepInputError(ValueError):
    """Raised when a seepage analysis cannot run because of invalid/missing input
    (e.g. an exit-face BC with materials that lack unsaturated parameters). Carries
    a user-facing message so callers (CLI, Studio) can surface it directly rather
    than crashing downstream on a None solution."""
    pass


def _min_distance_to_polyline(points, polyline):
    """Vectorized min distance from each point to a polyline (array of vertices).

    For each segment of the polyline, projects all points onto the segment
    and computes the clamped distance. Returns the minimum over all segments.
    """
    dists = np.full(len(points), np.inf)
    for i in range(len(polyline) - 1):
        a = polyline[i]
        b = polyline[i + 1]
        ab = b - a
        ab_sq = np.dot(ab, ab)
        if ab_sq < 1e-30:
            continue
        # Project each point onto segment ab, clamped to [0, 1]
        ap = points - a
        t = np.clip(ap @ ab / ab_sq, 0.0, 1.0)
        proj = a + t[:, None] * ab
        d = np.linalg.norm(points - proj, axis=1)
        np.minimum(dists, d, out=dists)
    return dists


def _boundary_edge_map(elements, element_types):
    """Boundary edges of the mesh, keyed by the sorted CORNER-node pair.

    Returns (boundary_edges, midside_map): boundary_edges maps each boundary
    edge (a corner pair owned by exactly one element) to that element index;
    midside_map maps a corner pair to its midside node for quadratic elements
    (absent for linear ones).
    """
    edge_counts = defaultdict(list)
    midside_map = {}
    for idx, en in enumerate(elements):
        et = element_types[idx]
        if et == 3:
            corners = [(en[0], en[1]), (en[1], en[2]), (en[2], en[0])]
        elif et == 6:
            corners = [(en[0], en[1]), (en[1], en[2]), (en[2], en[0])]
            midside_map[tuple(sorted((en[0], en[1])))] = en[3]
            midside_map[tuple(sorted((en[1], en[2])))] = en[4]
            midside_map[tuple(sorted((en[2], en[0])))] = en[5]
        elif et == 4:
            corners = [(en[0], en[1]), (en[1], en[2]), (en[2], en[3]), (en[3], en[0])]
        elif et in (8, 9):
            corners = [(en[0], en[1]), (en[1], en[2]), (en[2], en[3]), (en[3], en[0])]
            midside_map[tuple(sorted((en[0], en[1])))] = en[4]
            midside_map[tuple(sorted((en[1], en[2])))] = en[5]
            midside_map[tuple(sorted((en[2], en[3])))] = en[6]
            midside_map[tuple(sorted((en[3], en[0])))] = en[7]
        else:
            continue
        for a, b in corners:
            edge_counts[tuple(sorted((a, b)))].append(idx)

    boundary_edges = {e: elems[0] for e, elems in edge_counts.items() if len(elems) == 1}
    return boundary_edges, midside_map


def assemble_flux_nodal(nodes, elements, element_types, specified_fluxes, tolerance):
    """Consistent nodal loads for the Neumann (specified-flux) boundary condition.

    Each flux BC is a polyline carrying a uniform normal Darcy velocity q
    (positive = inflow). The load lives on the boundary EDGES that lie on that
    polyline, not on the nodes: a boundary edge belongs to the BC when both of
    its corner nodes are within `tolerance` of the polyline. The consistent
    load vector for the integral of N_i*q along a straight edge of length L is
    q*L/2 at each end for a linear edge, and (q*L/6, q*L/6, 2*q*L/3) for a
    quadratic (corner, corner, midside) edge — each set sums to q*L.

    Returns a (n_nodes,) vector to be ADDED to the RHS.
    """
    flux_nodal = np.zeros(len(nodes))
    if not specified_fluxes:
        return flux_nodal

    boundary_edges, midside_map = _boundary_edge_map(elements, element_types)
    if not boundary_edges:
        return flux_nodal

    edges = np.array(list(boundary_edges.keys()), dtype=int)   # (m, 2) corner pairs
    mids = np.array([midside_map.get(tuple(e), -1) for e in map(tuple, edges)], dtype=int)
    lengths = np.linalg.norm(nodes[edges[:, 1]] - nodes[edges[:, 0]], axis=1)

    for i, bc in enumerate(specified_fluxes):
        coords = np.asarray(bc["coords"], dtype=float)
        if len(coords) < 2:
            raise SeepInputError(
                f"Flux BC #{i + 1} has {len(coords)} coordinate(s); a flux BC is "
                "applied over the edges of a polyline and needs at least 2 points."
            )
        q = float(bc["flux"])

        d1 = _min_distance_to_polyline(nodes[edges[:, 0]], coords)
        d2 = _min_distance_to_polyline(nodes[edges[:, 1]], coords)
        on = (d1 <= tolerance) & (d2 <= tolerance)
        if not np.any(on):
            logger.warning("Flux BC #%d matched no boundary edges (flux=%g)", i + 1, q)
            continue

        # A PARTIAL match silently delivers less water than the user asked for: an
        # edge is loaded only when BOTH its corners lie on the polyline, so
        # endpoints landing mid-edge drop whole edges and q*L quietly shrinks. Warn
        # whenever the matched length misses the specified length by more than one
        # element, which is the most a correct match can lose to endpoint rounding.
        matched_len = float(np.sum(lengths[on]))
        spec_len = float(np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1)))
        h_typ = matched_len / max(int(np.count_nonzero(on)), 1)
        if spec_len - matched_len > h_typ:
            warnings.warn(
                f"Flux BC #{i + 1} (flux={q:g}) matched {matched_len:.4g} of its "
                f"{spec_len:.4g} specified length, so only "
                f"{matched_len / spec_len:.1%} of the intended inflow is applied. "
                "An edge is loaded only when both of its corner nodes lie on the "
                "polyline; check that the polyline endpoints fall on mesh nodes and "
                "that it follows the mesh boundary.",
                stacklevel=2,
            )

        qL = q * lengths[on]
        e_on = edges[on]
        m_on = mids[on]
        lin = m_on < 0
        quad = ~lin
        np.add.at(flux_nodal, e_on[lin, 0], 0.5 * qL[lin])
        np.add.at(flux_nodal, e_on[lin, 1], 0.5 * qL[lin])
        np.add.at(flux_nodal, e_on[quad, 0], qL[quad] / 6.0)
        np.add.at(flux_nodal, e_on[quad, 1], qL[quad] / 6.0)
        np.add.at(flux_nodal, m_on[quad], qL[quad] * 2.0 / 3.0)

    return flux_nodal


def _flux_inflow(flux_nodal, free_mask):
    """Inflow that specified-flux loads actually deliver into the mesh.

    `_dirichlet_system` seeds the RHS with the loads and then OVERWRITES the
    Dirichlet rows, so a load sitting on a Dirichlet node never enters the solve
    and must not be counted. The mask therefore has to be the RUNTIME Dirichlet
    set — for an unconfined problem that is the specified heads plus the ACTIVE
    exit-face nodes, which changes from iteration to iteration. Passing
    `bc_type == 0` instead would miss the loads on INACTIVE exit-face nodes,
    whose rows stay free and whose loads are real inflow (rain landing on the
    unsaturated part of a seepage face does infiltrate).
    """
    if flux_nodal is None:
        return 0.0
    f = np.asarray(flux_nodal, dtype=float)
    return float(np.sum(f[np.asarray(free_mask) & (f > 0)]))


def bc_node_tolerance(nodes):
    """How close a mesh node must lie to a boundary polyline to belong to it.

    A tenth of a typical element size, taken as the smaller of the mesh's x and y
    ranges over the square root of its node count. :func:`build_seep_data` binds
    every boundary condition with it, and preflight reads boundary lines against
    the ground surface with the same number."""
    nodes = np.asarray(nodes, dtype=float)
    x_range = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
    y_range = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    typical_element_size = min(x_range, y_range) / np.sqrt(len(nodes))
    return typical_element_size * 0.1


def build_seep_data(mesh, slope_data, seep_bc=1, check_inputs=True):
    """
    Build a seep_data dictionary from a mesh and data dictionary.
    
    This function takes a mesh dictionary (from build_mesh_from_polygons) and a data dictionary
    (from load_slope_data) and constructs a seep_data dictionary suitable for seep analysis.
    
    The function:
    1. Extracts mesh information (nodes, elements, element types, element materials)
    2. Builds material property arrays (k1, k2, alpha, kr0, h0) from the materials table
    3. Constructs boundary conditions by finding nodes that intersect with specified head
       and seep face lines from the data dictionary
    
    Parameters:
        mesh (dict): Mesh dictionary from build_mesh_from_polygons containing:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 3 or 4) of element node indices
            - element_types: np.ndarray (n_elements,) indicating 3 for triangles, 4 for quads
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
        data (dict): Data dictionary from load_slope_data containing:
            - materials: list of material dictionaries with k1, k2, alpha, kr0, h0 properties
            - seepage_bc: dictionary with "specified_heads" and "exit_face" boundary conditions
            - gamma_water: unit weight of water
        seep_bc (int): 1 for the main 'seep bc' set, 2 for the 'seep bc (2)' set.
        check_inputs (bool): Run :func:`xslope.preflight.preflight` on the model
            first and refuse with a :class:`~xslope.preflight.PreflightError` when it
            finds an error -- a conductivity of zero, a boundary set that drives no
            flow, a mesh built against a different material table. Default True.
            Pass False to bypass the checks.

    Returns:
        dict: seep_data dictionary with the following structure:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 3 or 4) of element node indices
            - element_types: np.ndarray (n_elements,) indicating 3 for triangles, 4 for quads
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
            - bc_type: np.ndarray (n_nodes,) of boundary condition flags (0=free, 1=fixed head, 2=exit face)
            - bc_values: np.ndarray (n_nodes,) of boundary condition values
            - flux_nodal: np.ndarray (n_nodes,) of consistent nodal loads from the
              specified-flux (Neumann) BCs, zero elsewhere (+ = inflow)
            - k1_by_mat: np.ndarray (n_materials,) of major conductivity values
            - k2_by_mat: np.ndarray (n_materials,) of minor conductivity values
            - angle_by_mat: np.ndarray (n_materials,) of angle values (degrees)
            - kr0_by_mat: np.ndarray (n_materials,) of relative conductivity values
            - h0_by_mat: np.ndarray (n_materials,) of suction head values
            - unit_weight: float, unit weight of water
    """
    
    # === PREFLIGHT ===
    # Checked here, but NOT refused here. Building a seep_data is not always a run:
    # `import_seep_solution` takes one purely as the shape of the mesh it needs to
    # read a stored field back for re-plotting, and a model carrying a vendor's
    # solved field legitimately has no conductivities and no boundary conditions of
    # its own. Refusing at the container build would block that. So the findings
    # ride along and `run_seepage_analysis` raises on them -- the gate lands on the
    # solve, which is the thing the preconditions are preconditions FOR.
    #
    # The mesh travels in the selection because a run is handed its mesh as an
    # argument: it is frequently built in memory and never stored on the model, so
    # the mesh rules must check the mesh this run will actually use.
    #
    # A model carrying a tseep sheet is checked as a TRANSIENT run, because that is
    # what it is: the storage, time-control and stage-time rules exist for the march
    # and have to be evaluated before it starts, not after it has finished and the
    # staging is asked for.
    preflight_report = None
    if check_inputs:
        from .preflight import preflight as _preflight
        _analysis = 'tseep' if slope_data.get('tseep') else 'seep'
        preflight_report = _preflight(slope_data, _analysis, {'mesh': mesh})

    # Extract mesh data
    nodes = mesh["nodes"]
    elements = mesh["elements"]
    element_types = mesh["element_types"]
    element_materials = mesh["element_materials"]
    
    # Initialize boundary condition arrays
    n_nodes = len(nodes)
    bc_type = np.zeros(n_nodes, dtype=int)  # 0 = free, 1 = fixed head, 2 = exit face
    bc_values = np.zeros(n_nodes)
    
    # Build material property arrays
    materials = slope_data["materials"]
    n_materials = len(materials)
    
    k1_by_mat = np.zeros(n_materials)
    k2_by_mat = np.zeros(n_materials)
    angle_by_mat = np.zeros(n_materials)
    kr0_by_mat = np.zeros(n_materials)
    h0_by_mat = np.zeros(n_materials)
    unsat_by_mat = np.zeros(n_materials, dtype=int)   # KR_LF / KR_VG per material
    vg_a_by_mat = np.zeros(n_materials)               # van Genuchten alpha
    vg_n_by_mat = np.zeros(n_materials)               # van Genuchten n
    # v26 Mualem pore-connectivity exponent. Blank on the sheet, and every file
    # written before the column existed, reads as Mualem's own 0.5.
    vg_l_by_mat = np.full(n_materials, 0.5)
    ss_by_mat = np.zeros(n_materials)                 # v18 specific storage [1/len]
    sy_by_mat = np.zeros(n_materials)                 # v18 specific yield [-]
    material_names = []

    for i, material in enumerate(materials):
        k1_by_mat[i] = material.get("k1", 1.0)
        k2_by_mat[i] = material.get("k2", 1.0)
        angle_by_mat[i] = material.get("alpha", 0.0)
        kr0_by_mat[i] = material.get("kr0", 0.001)
        h0_by_mat[i] = material.get("h0", -1.0)
        _u = str(material.get("unsat", "lf")).strip().lower()
        unsat_by_mat[i] = {"vg": KR_VG, "gard": KR_GARD}.get(_u, KR_LF)
        vg_a_by_mat[i] = material.get("vg_a", 0.0)
        vg_n_by_mat[i] = material.get("vg_n", 0.0)
        _l = material.get("vg_l")
        vg_l_by_mat[i] = 0.5 if _l is None else float(_l)
        # v18 transient storage: None (blank) -> 0 here; a transient run is refused
        # by preflight's tseep.storage_nonpositive before it gets here, so a 0 only
        # survives on the steady path where it is never read.
        _ss = material.get("Ss")
        _sy = material.get("Sy")
        ss_by_mat[i] = 0.0 if _ss is None else float(_ss)
        sy_by_mat[i] = 0.0 if _sy is None else float(_sy)
        material_names.append(material.get("name", f"Material {i+1}"))
    
    # Process boundary conditions
    if seep_bc == 2:
        seepage_bc = slope_data.get("seepage_bc2", {})
    else:
        seepage_bc = slope_data.get("seepage_bc", {})
    
    # How close a node must lie to a boundary polyline to belong to it.
    tolerance = bc_node_tolerance(nodes)
    
    logger.debug("Mesh tolerance for boundary conditions: %.6f", tolerance)
    
    # Series bindings for time-varying (transient) BCs. A head/flux VALUE that is a
    # string names a `tseep` series (fileio validated it). These stay OUT of the
    # baked scalar arrays: build_seep_data leaves the nodes free, and
    # run_transient_seepage resolves them per time step (by the head binding's kind:
    # the reservoir rule, or a plain time-varying head; scaled unit flux loads for
    # fluxes). Empty on every steady file, so the steady solve is
    # bit-identical — and run_seepage_analysis refuses to solve while bindings
    # exist, so a series-bound node can never be silently left free in a steady
    # solution.
    head_series_bindings = []   # [{"mask": bool(n_nodes), "series": name, "kind": str}]
    flux_series_bindings = []   # [{"unit_flux_nodal": (n_nodes,), "series": name}]

    # Each binding also carries the series' value AT t = 0 ("value_t0"), computed
    # here where the tseep sheet is in reach. A STEADY solve on a series-bound
    # model reads every series at t = 0 — the same values the transient march's
    # own initial condition uses — so the steady run is the initial-condition
    # snapshot, cheap, without the march behind it. It is read by the march's own
    # evaluator (_eval_series on the anchors build_tseep_data makes), so a step
    # anchored at t = 0 reads its value BEFORE the step, exactly as the march's
    # initial condition does.
    _series_t0 = {}
    _ts = (slope_data or {}).get("tseep") or {}
    for _name, _vals in (_ts.get("series") or {}).items():
        _pairs = [(float(t), float(v)) for t, v in zip(_ts.get("times") or [], _vals)
                  if v is not None]
        if _pairs:
            _ta, _va = zip(*_pairs)
            _series_t0[_name] = float(_eval_series(
                (np.asarray(_ta, dtype=float), np.asarray(_va, dtype=float)), 0.0))

    # Process specified head boundary conditions
    # Vectorized: compute distance from all nodes to each BC line at once
    specified_heads = seepage_bc.get("specified_heads", [])
    for bc in specified_heads:
        head_value = bc["head"]
        coords = bc["coords"]
        # v18: two Dirichlet head types. "head" = a plain Dirichlet held at the
        # value/series at every node of the polyline, at all times. "reservoir" =
        # a node is held at the level while submerged (elevation <= level); a node
        # the water has reached and then uncovered is a seepage-exit face; a node
        # the water has never reached is no-flow. A CONSTANT level has no history,
        # so its nodes above the level are no-flow and it never makes an exit face.
        # With the polyline drawn at or below the level (the usual case) the two
        # types bake the same plain Dirichlet.
        kind = str(bc.get("kind", "head")).strip().lower()

        if len(coords) < 2:
            continue

        # Compute min distance from each node to any segment of the BC line
        seg_coords = np.array(coords)
        dists = _min_distance_to_polyline(nodes, seg_coords)
        mask = dists <= tolerance
        if isinstance(head_value, str):
            # Time-varying head: resolved per step by run_transient_seepage (see
            # resolve_transient_bc), which applies the reservoir rule or the plain
            # head rule by the binding's kind. It stays OUT of the baked scalar
            # arrays here.
            head_series_bindings.append({"mask": mask, "series": head_value,
                                         "kind": kind,
                                         "value_t0": _series_t0.get(head_value)})
            continue
        if kind == "reservoir":
            # A constant level: hold the submerged nodes at it. The nodes above it
            # have never been under water, so they stay no-flow (or whatever
            # another boundary makes them).
            _apply_reservoir_level(bc_type, bc_values, mask,
                                   np.asarray(nodes[:, 1], dtype=float),
                                   float(head_value))
        else:
            bc_type[mask] = 1
            bc_values[mask] = head_value

    # Process seep face (exit face) boundary conditions.
    # A node that already carries a specified head (bc_type == 1) keeps it:
    # specified head is a true Dirichlet condition and must win over the
    # free-boundary exit face at shared corner nodes (e.g. the downstream toe,
    # which lies on both the exit-face line and the tailwater head line). If
    # the exit face is allowed to claim that corner, the active-set logic can
    # deactivate it and the fixed tailwater head floats free, raising the whole
    # downstream field and depressing the through-flow.
    exit_face_coords = seepage_bc.get("exit_face", [])
    if len(exit_face_coords) >= 2:
        seg_coords = np.array(exit_face_coords)
        dists = _min_distance_to_polyline(nodes, seg_coords)
        mask = (dists <= tolerance) & (bc_type != 1)
        bc_type[mask] = 2
        bc_values[mask] = nodes[mask, 1]  # Use node's y-coordinate as elevation

    # Specified-flux (Neumann) BCs. These are NOT Dirichlet: they stay bc_type 0
    # and enter the system only through the RHS, as consistent nodal loads
    # assembled over the boundary edges lying on each flux polyline.
    #
    # A string flux VALUE names a tseep series (a time-varying flux). Flux load is
    # LINEAR in the value, so we assemble a UNIT (flux=1) load vector for that
    # polyline once and record it; run_transient_seepage scales it by series(t) each
    # step. Numeric fluxes are baked into flux_nodal exactly as before.
    specified_fluxes = seepage_bc.get("specified_fluxes", [])
    _numeric_fluxes = [b for b in specified_fluxes if not isinstance(b["flux"], str)]
    for _b in specified_fluxes:
        if isinstance(_b["flux"], str):
            unit_fn = assemble_flux_nodal(
                nodes, elements, element_types,
                [{"coords": _b["coords"], "flux": 1.0}], tolerance)
            flux_series_bindings.append({"unit_flux_nodal": unit_fn,
                                         "series": _b["flux"],
                                         "value_t0": _series_t0.get(_b["flux"])})
    flux_nodal = assemble_flux_nodal(nodes, elements, element_types,
                                     _numeric_fluxes, tolerance)
    # A series-bound head resolves to a Dirichlet (or exit-face) node at runtime, so
    # it anchors the head even though its build-time bc_type is 0.
    _has_dirichlet = bool(np.any(bc_type > 0) or head_series_bindings)
    if len(specified_fluxes) > 0 and not _has_dirichlet:
        raise SeepInputError(
            "Seepage problem has specified-flux boundary conditions but no "
            "Dirichlet boundary (no specified head and no exit face). With only "
            "Neumann conditions the head is defined solely up to an additive "
            "constant and the system is singular. Add a specified-head or "
            "exit-face boundary condition."
        )

    # Check for missing unsaturated parameters when exit face BCs are present
    has_exit_face = np.any(bc_type == 2)
    missing_unsat_params = False
    if has_exit_face:
        bad_mats = []
        def _bad(v):
            return v is None or (isinstance(v, (int, float)) and (v == 0 or np.isnan(v)))
        for i, material in enumerate(materials):
            name = material.get("name", f"Material {i+1}")
            if unsat_by_mat[i] == KR_VG:
                # van Genuchten needs alpha > 0 and n > 1.
                a, nn = material.get("vg_a", 0), material.get("vg_n", 0)
                if _bad(a) or nn is None or (isinstance(nn, (int, float)) and nn <= 1):
                    bad_mats.append(f"  Material {i+1} ({name}, unsat=vg): a={a}, n={nn} (need a>0, n>1)")
            elif unsat_by_mat[i] == KR_GARD:
                # Gardner kr = 1/(1 + a*psi^n) needs a > 0 and n > 0. (Unlike van
                # Genuchten, n is not required to exceed 1 — there is no m = 1-1/n.)
                a, nn = material.get("vg_a", 0), material.get("vg_n", 0)
                if _bad(a) or _bad(nn) or (isinstance(nn, (int, float)) and nn <= 0):
                    bad_mats.append(f"  Material {i+1} ({name}, unsat=gard): a={a}, n={nn} (need a>0, n>0)")
            else:
                # Linear front needs kr0 (>0) and h0 (<0).
                mat_kr0, mat_h0 = material.get("kr0", 0), material.get("h0", 0)
                if _bad(mat_kr0) or _bad(mat_h0):
                    bad_mats.append(f"  Material {i+1} ({name}, unsat=lf): kr0={mat_kr0}, h0={mat_h0} (need kr0>0, h0<0)")
        if bad_mats:
            missing_unsat_params = True
            print("\n" + "="*70)
            print("WARNING: Exit face boundary condition detected but the following")
            print("materials are missing valid unsaturated parameters:")
            for line in bad_mats:
                print(line)
            print("\nThe unsaturated seepage solver requires, per material, valid")
            print("kr0 (>0) and h0 (<0) for the linear-front model, a (>0) and n (>1)")
            print("for van Genuchten, or a (>0) and n (>0) for Gardner. Set these in")
            print("the input file.")
            print("="*70 + "\n")

    # Get unit weight of water
    unit_weight = require_gamma_water(slope_data, "seepage analysis")

    # Construct seep_data dictionary
    seep_data = {
        "nodes": nodes,
        "elements": elements,
        "element_types": element_types,
        "element_materials": element_materials,
        "bc_type": bc_type,
        "bc_values": bc_values,
        "flux_nodal": flux_nodal,
        "k1_by_mat": k1_by_mat,
        "k2_by_mat": k2_by_mat,
        "angle_by_mat": angle_by_mat,
        "kr0_by_mat": kr0_by_mat,
        "h0_by_mat": h0_by_mat,
        "unsat_by_mat": unsat_by_mat,
        "vg_a_by_mat": vg_a_by_mat,
        "vg_n_by_mat": vg_n_by_mat,
        "vg_l_by_mat": vg_l_by_mat,
        # v18 transient storage (per material). Unused by the steady solver.
        "ss_by_mat": ss_by_mat,
        "sy_by_mat": sy_by_mat,
        # v18 transient BC series bindings (empty on every steady file).
        "head_series_bindings": head_series_bindings,
        "flux_series_bindings": flux_series_bindings,
        # BC geometry + which set was baked, carried so plot_seep_solution can draw
        # each head/reservoir boundary's instantaneous WATER LEVEL over a solution
        # frame (opt-in show_bc_levels). These are display metadata only — the solver
        # reads the baked bc_type/bc_values arrays, never these — so a seep_data
        # missing them (e.g. reconstructed from JSON) simply draws no levels.
        "seepage_bc": seepage_bc,
        "seep_bc": seep_bc,
        "tseep": slope_data.get("tseep"),
        "material_names": material_names,
        "unit_weight": unit_weight,
        "missing_unsat_params": missing_unsat_params,
        # Per-material detail lines for the missing params (empty when all valid),
        # so a caller can build a specific error message without re-deriving it.
        "missing_unsat_detail": bad_mats if missing_unsat_params else [],
    }

    # Carry the unit declaration through as provenance so it survives into
    # save_seep_data_to_json (units plan phase 2, sec 2a). Recorded only when set, so an
    # undeclared model's seep_data is byte-identical to before.
    if slope_data.get("unit_system"):
        seep_data["unit_system"] = slope_data["unit_system"]
    if slope_data.get("time_unit"):
        seep_data["time_unit"] = slope_data["time_unit"]

    # The preflight findings ride with the data they describe, so the solve can
    # refuse on them while a re-plot or a stored-field import never sees them.
    # Absent entirely when check_inputs=False, which is what makes that the escape
    # hatch for the solve as well as for the build.
    if preflight_report is not None:
        seep_data["_preflight"] = preflight_report

    return seep_data


def import_seep2d(filepath):
    """
    Reads SEEP2D .s2d input file and returns mesh, materials, and BC data.
    Supports both triangular and quadrilateral elements.
    Uses implicit numbering (0-based array indices) instead of explicit node IDs.
    
    Note: All node indices in elements are converted to 0-based indexing during import.
    Material IDs remain 1-based as they appear in the SEEP2D file.

    Returns:
        {
            "nodes": np.ndarray (n_nodes, 2),
            "bc_type": np.ndarray (n_nodes,),   # boundary condition flags
            "bc_values": np.ndarray (n_nodes,),    # boundary condition values (head or elevation)
            "elements": np.ndarray (n_elements, 3 or 4),  # triangle or quad node indices (0-based)
            "element_types": np.ndarray (n_elements,),    # 3 for triangles, 4 for quads
            "element_materials": np.ndarray (n_elements,) # material IDs (1-based)
        }
    """
    import re

    with open(filepath, "r", encoding="latin-1") as f:
        lines = [line.rstrip() for line in f if line.strip()]

    title = lines[0]                  # First line is the title (any text)
    parts = lines[1].split()          # Second line contains analysis parameters

    num_nodes = int(parts[0])         # Number of nodes
    num_elements = int(parts[1])      # Number of elements
    num_materials = int(parts[2])     # Number of materials
    datum = float(parts[3])           # Datum elevation (not used, assume 0.0)

    problem_type = parts[4]           # "PLNE" = planar, otherwise axisymmetric (we only support "PLNE")
    analysis_flag = parts[5]          # Unknown integer (ignore)
    flow_flag = parts[6]              # "F" or "T" = compute flowlines (ignore)
    unit_weight = float(parts[7])     # Unit weight of water (e.g. 62.4 lb/ft³ or 9.81 kN/m³)
    model_type = int(parts[8])        # 1 = linear front, 2 = van Genuchten (we only support 0)

    assert problem_type == "PLNE", "Only planar problems are supported"
    assert model_type == 1, "Only linear front models are supported"

    unit_weight = float(parts[7])   # the unit weight
    mat_props = []
    line_offset = 2
    while len(mat_props) < num_materials:
        nums = [float(n) if '.' in n or 'e' in n.lower() else int(n)
                for n in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', lines[line_offset])]
        if len(nums) >= 6:
            mat_props.append(nums[:6])
        line_offset += 1
    mat_props = np.array(mat_props)
    k1_array = mat_props[:, 1]
    k2_array = mat_props[:, 2]
    angle_array = mat_props[:, 3]
    kr0_array = mat_props[:, 4]
    h0_array = mat_props[:, 5]
    node_lines = lines[line_offset:line_offset + num_nodes]
    element_lines = lines[line_offset + num_nodes:]

    coords = []
    bc_type = []
    bc_values = []

    for line in node_lines:
        try:
            node_id = int(line[0:5])
            bc_type_val = int(line[7:10])
            x = float(line[10:25])
            y = float(line[25:40])

            if bc_type_val == 1 and len(line) >= 41:
                bc_value = float(line[40:55])
            elif bc_type_val == 2:
                bc_value = y
            else:
                bc_value = 0.0

            bc_type.append(bc_type_val)
            bc_values.append(bc_value)
            coords.append((x, y))

        except Exception as e:
            print(f"Warning: skipping node due to error: {e}")

    elements = []
    element_mats = []
    element_types = []

    for line in element_lines:
        nums = [int(n) for n in re.findall(r'\d+', line)]
        if len(nums) >= 6:
            _, n1, n2, n3, n4, mat = nums[:6]
            
            # Convert to 0-based indexing during reading
            n1, n2, n3, n4 = n1 - 1, n2 - 1, n3 - 1, n4 - 1
            
            # Check if this is a triangle (n3 == n4) or quad (n3 != n4)
            if n3 == n4:
                # Triangle: repeat the last node to create 4-node format
                elements.append([n1, n2, n3, n3])
                element_types.append(3)
            else:
                # Quadrilateral: use all 4 nodes
                elements.append([n1, n2, n3, n4])
                element_types.append(4)
            
            element_mats.append(mat)

    nodes_arr = np.array(coords)
    elements_arr = np.array(elements, dtype=int)  # Already 0-based
    element_types_arr = np.array(element_types, dtype=int)

    # Defensive: normalize any clockwise elements to CCW (the tri3 assembly
    # skips CW elements). SEEP2D/GMS meshes are normally CCW already.
    from xslope.mesh import ensure_ccw_elements
    ensure_ccw_elements(nodes_arr, elements_arr, element_types_arr)
    # restore the repeat-last-node padding convention for triangles
    tri_rows = element_types_arr == 3
    elements_arr[tri_rows, 3] = elements_arr[tri_rows, 2]

    return {
        "nodes": nodes_arr,
        "bc_type": np.array(bc_type, dtype=int),
        "bc_values": np.array(bc_values),
        "elements": elements_arr,
        "element_types": element_types_arr,
        "element_materials": np.array(element_mats),
        "k1_by_mat": k1_array,
        "k2_by_mat": k2_array,
        "angle_by_mat": angle_array,
        "kr0_by_mat": kr0_array,
        "h0_by_mat": h0_array,
        "unit_weight": unit_weight,
        # A .s2d file declares no unit system -- only the bare unit weight of water --
        # so infer SI/Imperial from its magnitude (None if off-band). No time_unit: the
        # format carries no time-base declaration for its conductivities (units plan
        # phase 2, sec 2d). LABEL only; the numbers are used exactly as read.
        "unit_system": infer_system_from_gamma_water(unit_weight),
    }


def solve_confined(nodes, elements, bc_type, dirichlet_bcs, k1_vals, k2_vals, angles=None,
                   element_types=None, flux_nodal=None):
    """
    FEM solver for confined seep with anisotropic conductivity.
    Supports triangular and quadrilateral elements with both linear and quadratic shape functions.
    
    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 9) element node indices (padded with zeros for unused nodes)
        bc_type : (n_nodes,) array of boundary condition flags
        dirichlet_bcs : list of (node_id, head_value)
        k1_vals : (n_elements,) or scalar, major axis conductivity
        k2_vals : (n_elements,) or scalar, minor axis conductivity
        angles : (n_elements,) or scalar, angle in degrees (from x-axis)
        element_types : (n_elements,) array indicating:
            3 = 3-node triangle (linear)
            4 = 4-node quadrilateral (bilinear)  
            6 = 6-node triangle (quadratic)
            8 = 8-node quadrilateral (serendipity)
            9 = 9-node quadrilateral (Lagrange)
    Returns:
        head : (n_nodes,) array of nodal heads
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals, angles)
    data = _assembly_data(asm)
    A_full = _csr_from_data(asm, data)

    dir_mask = np.zeros(n_nodes, dtype=bool)
    dir_values = np.zeros(n_nodes)
    for node, value in dirichlet_bcs:
        dir_mask[node] = True
        dir_values[node] = value
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values, neumann=flux_nodal)

    head = spsolve(A, b)
    q = A_full @ head

    # Inflow = what the Dirichlet boundaries supply, plus what the flux loads
    # deliver at the free nodes.
    #
    # Only the loads on FREE rows are ever applied: `_dirichlet_system` seeds b with
    # f_ext and then overwrites the Dirichlet rows, so a load on a Dirichlet node is
    # dropped from the system entirely. `f_eff` is therefore the load the solve
    # actually saw, and A.h - f_eff is the boundary reaction — zero on free rows (the
    # equation just solved) and the full external supply on Dirichlet rows. Summing
    # the positive reactions and the free-node inflow then closes exactly, because A
    # has zero row sums and so the two sides of sum(A.h) = 0 must balance.
    f_ext = (np.zeros(n_nodes) if flux_nodal is None
             else np.asarray(flux_nodal, dtype=float))
    f_eff = np.where(dir_mask, 0.0, f_ext)
    reaction = q - f_eff
    total_flow = float(np.sum(reaction[dir_mask & (reaction > 0)]))
    total_flow += _flux_inflow(f_ext, ~dir_mask)

    return head, A, q, total_flow


def _require_runtime_dirichlet(dir_mask):
    """Guard the unconfined solve against a RUNTIME-singular system.

    The static input guards (build_seep_data / run_seepage_analysis) only see
    the DECLARED boundary conditions: they accept a model as long as it carries
    a specified head OR an exit face. But the exit face is a free boundary whose
    active (Dirichlet) set is recomputed every iteration, and it can empty
    mid-solve — every exit-face node goes unsaturated and drops out. When it does
    and there is no specified head anywhere, the effective Dirichlet set is empty:
    the stiffness matrix is pure Neumann (singular), head is defined only up to
    an additive constant, and spsolve returns garbage (~1e15) instead of failing.

    Testing the effective set each iteration turns that silent divergence into a
    clear, user-actionable SeepInputError — routed exactly like the static input
    guards (studio's SeepRunner prints SeepInputError cleanly, no traceback)."""
    if not np.any(dir_mask):
        raise SeepInputError(
            "Unconfined seepage solve became singular: the exit face fully "
            "deactivated and no specified-head boundary remains, so the head is "
            "defined only up to an additive constant. Add a specified head, or "
            "check the exit-face definition (an exit face entirely above the "
            "water table drains nothing and cannot anchor the head)."
        )


#: Limit-cycle escape, and the one relaxation floor below the ladder's 0.01.
#:
#: A handful of interior nodes on a high-conductivity-contrast interface can settle
#: into an exact periodic orbit: the relaxed iterate returns to a sweep 2 to 25 back,
#: to within :data:`_CYCLE_TOL` x domain height, with the exit-face set unchanged, and
#: repeats that indefinitely. The ladder's terminal relax = 0.01 is what hides the
#: orbit — the iterate barely moves, while the un-relaxed step the convergence gate
#: tests lands far from it, so the solve neither converges nor visibly diverges.
#:
#: Detected only as a RETURN — a repeat at a period of 2 or more that the iterate had
#: travelled further than the tolerance away from at every shorter lag — held for
#: :data:`_CYCLE_PERSIST` consecutive sweeps. A period of 1 is not a cycle: it is
#: ordinary stagnation near a fixed point, which is what vp046b and vp077a do for
#: hundreds of sweeps before converging normally, and dropping the floor under them
#: costs them their convergence. Nor is a repeat the iterate never left: a fine enough
#: creep sits inside the tolerance at every lag and reads as a period of 2. And the
#: persistence window is what separates a true orbit from the short-lived repeats a
#: converging solve passes through (johnson_res holds one for 7 sweeps).
#:
#: THE FLOOR DROPS ONCE, TO 1e-3, AND NO FURTHER — it is a window, not a ladder.
#: Measured over the two fields that cycle (earth_dam_rapid, whose section, mesh and
#: first boundary set gsat_seep and earth_dam_down_seep share, and johnson_rapid_KEY),
#: each given 2000 sweeps: at 3e-3 they close in 348 and 459 sweeps, at 1e-3 in 525 and
#: 472, and at 3e-4 in 1142 and 525, every field within 6e-4 ft of head (0.04 psf) of
#: the 1e-3 answer. Under a floor applied from the start rather than on detection,
#: johnson_rapid_KEY does not converge at 1e-3, 3e-4 or 1e-4, so the floor has to wait
#: for the detection; once it has, 1e-3 sits inside the window on both fields.
_CYCLE_LOOKBACK = 25      # how many sweeps back a repeat is looked for
#: The repeat tolerance sets how soon an orbit is SEEN, not whether it is one: the
#: orbit's return distance shrinks by about half every 65 sweeps, and detection waits
#: for it to fall under the tolerance. At 1e-6 the earth dam's orbit (earth_dam_rapid
#: on its 3.5 ft mesh) is seen at sweep 491, past the 400 sweeps Studio allows a
#: steady solve; at 1e-5 it is seen at 261 and the solve closes at 525, its field
#: within 0.003 ft of head of the 1e-6 answer. Every other steady companion in the
#: corpus (tools/make_seep_sidecars.py, 29 more solves) is bit for bit the same at
#: either value.
_CYCLE_TOL = 1e-5         # repeat tolerance, x domain height
_CYCLE_PERSIST = 30       # consecutive sweeps the repeat must hold before acting
_CYCLE_RELAX_FLOOR = 1e-3
#: Sweeps granted beyond ``max_iter``, once, on the sweep the floor drops and never
#: otherwise. A solve stepping at a thousandth needs the room, and the sweep it
#: detects its orbit on is not something the caller can budget for: earth_dam_rapid
#: detects at 261 and closes at 525, so at the library's own 400-sweep default it
#: would otherwise stop with its orbit diagnosed and not yet damped. The extension is
#: granted from that sweep and is not renewed; a solve that has not closed by then stops.
_CYCLE_EXTRA_SWEEPS = 600

#: Set-revisit escape for a cycling quadratic exit-face edge.
#:
#: The all-or-nothing edge rule below keeps a quadratic seepage side fully active or
#: fully inactive. On earth_dam2 that forbids the one self-consistent set the solve
#: is reaching for (midside active, corner not), so the exit-face active set orbits
#: between two sets it has already visited, for as long as it is allowed to run.
#:
#: The gate is that orbit and nothing else: the active set RETURNS TO A SET IT HAS
#: ALREADY VISITED, after sweep :data:`_SET_REVISIT_SWEEP`. Every converging model in
#: the corpus revisits sets too — while its seepage face is still finding its extent —
#: but is done doing so early: the latest last-revisit measured over the corpus is
#: sweep 41 (vp077a), and earth_dam2 is still revisiting at sweep 992.
_SET_REVISIT_SWEEP = 100
_WET_FACE_FORCE_MAX = 3      # times an inactive exit node under pressure is joined to the face at set stability

#: Pinned-element escape, for a solve the budget leaves unconverged.
#:
#: One element on the wet/dry edge can hold a solve in a loop that no relaxation the
#: solver uses will close. One of the element's kr sampling points sits at the bend of
#: the kr curve: a drier element backs the water up behind it and wets it, a wetter
#: one drains it, and the element's kr comes back from each sweep hundreds to
#: hundreds of thousands of times further from its fixed point than it went in.
#: Relaxed Picard is stable on that loop only below 2 / (1 - gain). Measured on the
#: earth dam of test/transient_seep_check.py (tri6, 9 to 32 divisions across the
#: section), the gain runs from -210 to -2.7e5, so the factor needed runs from 1e-2
#: down to 7e-6: under the ladder's last rung on every one of those meshes, and under
#: the limit-cycle floor on four of them.
#:
#: The escape takes that one element out of the Picard iteration. Every sweep solves
#: for the element's kr EXACTLY against the current kr of every other element: one
#: factorization of the sweep's matrix and a rank-nn update of it for the element
#: (Sherman-Morrison-Woodbury) make the element's nodal heads a closed form in its kr,
#: and the self-consistent kr is then a one-dimensional root, always bracketed by the
#: element's kr floor and 1 because kr itself lies between them. Every other element
#: iterates as before. The fixed point reached is the full problem's own, and it is
#: accepted only by the solver's ordinary gate, on one ordinary sweep at relax 1.
#:
#: The candidates, tried in turn: the element whose kr swung most across the last
#: :data:`_CYCLE_LOOKBACK` sweeps, which is the orbit's own signature, and then the
#: elements whose kr :data:`_PIN_STEPS` un-relaxed steps from the last iterate move
#: most, for a solve creeping at the limit-cycle floor rather than orbiting, where
#: nothing swings. The un-relaxed steps are a coarse reading on their own: each moves
#: the field far enough to flip the culprit's neighbors too, and on five of the ten
#: orbiting meshes below the culprit is not their first pick. A wrong element shows
#: itself: the loop it does not hold keeps cycling, and the sweeps reach the ladder's
#: last rung and set no new low there, which abandons it.
#:
#: Measured on the 11 meshes of that dam the budget leaves unconverged (9 to 17, 20
#: and 32 divisions), at a budget of 400 sweeps and of 2000: on the ten that orbit the
#: swing names the element and it closes in 33 to 132 sweeps; on the one that creeps
#: (14 divisions) nothing swings, the first step candidate is abandoned at sweep 62
#: and the second closes in 257. No right element needed a rung finer than 0.1, and
#: every gate sweep came in under 2e-6 of the tolerance. The two budgets land on the
#: same flow rate to eight digits.
#:
#: Read only once the sweep budget is spent unconverged, so a solve that converges
#: within its budget never reaches it, and a rescue that does not pass the gate leaves
#: the unconverged answer as the budget left it. It spends at most the budget again.
_PIN_MIN_SWING = 0.1     # decades of kr: an element that moved less is not the loop
_PIN_STEPS = 3           # un-relaxed steps read for candidates, and candidates taken from them
#: The pinned sweeps' head tolerance, as a fraction of the gate's. The gate sweep
#: carries what is left of the element's error at up to the loop's gain (2.7e5
#: measured), so the inner tolerance sits that far under the gate.
_PIN_TOL = 1e-6


def _element_kr_avg(asm, p_nodes, kr_kw):
    """Each element's kr averaged at its sampling points, in element order: the factor
    :func:`_assembly_data` scales the element's saturated matrix by in head mode."""
    out = np.empty(sum(len(g['idx']) for g in asm['groups']))
    for g in asm['groups']:
        idx = g['idx']
        p_gp = p_nodes[g['conn']] @ g['N_kr'].T
        kr = kr_relative_vec(p_gp, **{k: None if v is None else v[idx][:, None]
                                      for k, v in kr_kw.items()})
        out[idx] = (kr @ g['w_kr']) / g['wsum']
    return out


def _pinned_element_solve(asm, kr_kw, y, h, e, dir_mask, dir_values, f_ext, tol,
                          max_sweeps):
    """Steady sweeps from ``h`` with element ``e``'s kr solved for exactly at every
    sweep (see :data:`_PIN_MIN_SWING`), on a fixed exit-face set. Returns
    ``(head, converged, sweeps)``; not converged when the ladder stalls on its last
    rung, which is how a wrong element shows itself, or the sweeps run out.

    The element's matrix is its kr times its saturated matrix, so changing its kr by
    ``dk`` adds ``dk`` times the saturated matrix, less the rows the Dirichlet system
    replaces by identity. With ``Z`` the solved unit columns of the element's nodes,
    its nodal heads at any kr are a closed form needing no further factorization."""
    g = next(g for g in asm['groups'] if np.any(g['idx'] == e))
    j = int(np.flatnonzero(g['idx'] == e)[0])
    conn = g['conn'][j]
    nn = len(conn)
    Kt = np.where(dir_mask[conn][:, None], 0.0, g['ke'][j])
    one = np.array([e])
    kw_e = {k: None if v is None else v[one][:, None] for k, v in kr_kw.items()}

    def kr_e(h_e):
        kr = kr_relative_vec(((h_e - y[conn]) @ g['N_kr'].T)[None, :], **kw_e)[0]
        return float(kr @ g['w_kr']) / g['wsum']

    log_floor = np.log10(max(kr_e(y[conn] - 1e12), 1e-300))
    E = np.zeros((len(h), nn))
    E[conn, np.arange(nn)] = 1.0
    eye = np.eye(nn)
    rung, best, stall = 0, float('inf'), 0
    for sweep in range(1, max_sweeps + 1):
        data = _assembly_data(asm, p_nodes=h - y, mode='head', **kr_kw)
        k_now = kr_e(h[conn])
        A, b = _dirichlet_system(asm, data, dir_mask, dir_values, neumann=f_ext)
        lu = splu(A.tocsc())
        h_base = lu.solve(b)
        Z = lu.solve(E)
        Ze, he_base = Z[conn], h_base[conn]

        def update(kr):
            dk = kr - k_now
            return dk, np.linalg.solve(eye + dk * Kt @ Ze, Kt @ he_base)

        # kr out minus kr in falls across [floor, 1]; bisect it in log kr to the last
        # bit, each trial a small dense solve.
        lo, hi = log_floor, 0.0
        while True:
            mid = 0.5 * (lo + hi)
            if not lo < mid < hi:
                break
            dk, w = update(10.0 ** mid)
            if kr_e(he_base - dk * Ze @ w) > 10.0 ** mid:
                lo = mid
            else:
                hi = mid
        dk, w = update(10.0 ** (0.5 * (lo + hi)))
        h_solved = h_base - dk * Z @ w
        step = float(np.max(np.abs(h_solved - h)))
        if step < tol:
            return h_solved, True, sweep
        # The ladder, stepping down only: these sweeps start from a field already
        # near its answer, and a stall on the last rung is the wrong element.
        if step < best:
            best, stall = step, 0
        else:
            stall += 1
            if stall >= _LADDER_PATIENCE:
                if rung + 1 == len(_RELAX_LADDER):
                    return h, False, sweep
                rung, best, stall = rung + 1, float('inf'), 0
        h = h + _RELAX_LADDER[rung] * (h_solved - h)
    return h, False, max_sweeps

#: PROGRESS-FED RELAXATION LADDER. The under-relaxation factor walks this ladder, and
#: it moves on what the head change is DOING, never on the sweep number:
#:
#:   * STEP DOWN one rung after :data:`_LADDER_PATIENCE` consecutive sweeps on which
#:     the head change set no new low for the rung, or the iterate returned to a recent
#:     one -- the solve is stalling or orbiting at this relaxation and will not close
#:     at it;
#:   * HOLD while it keeps setting new lows;
#:   * STEP BACK UP one rung after a run of consecutive new lows -- :data:`_LADDER_RECOVER`
#:     of them the first time a rung is climbed into, doubled at each later attempt on
#:     that same rung.
#:
#: _LADDER_RECOVER is twice _LADDER_PATIENCE: it takes twice the evidence to decide a
#: rung is finer than the solve needs as it takes to decide it is too coarse, which is
#: the asymmetry that keeps a coarse rung from being re-entered on a lucky run.
#:
#: NEW LOW, NOT "FELL SINCE LAST SWEEP". An oscillating solve falls on half its sweeps
#: and would hold a rung forever under the simpler reading: measured on the heap column
#: of SEEPW-T05, whose head change bounces between 0.9 and 5.7 ft at relax = 1, the
#: ladder sat at 1.0 for 600 sweeps before a run of five rises finally happened to
#: occur. A best-so-far test reads that same trajectory as the stall it is and steps
#: down within a dozen sweeps.
#:
#: The doubling on the way up is what bounds the ladder. A rung too coarse for the
#: field as it stands now may be right for the field fifty sweeps later, so it is not
#: shut out for good; but each failed attempt costs twice the evidence of the last, so
#: a solve can only spend a logarithmic number of sweeps re-testing a rung instead of
#: thrashing between two of them.
#:
#: Sweep count is not a measure of how hard a problem is. The fixed schedule this
#: replaces stepped to 0.01 by sweep 121 on every model that got that far and then held
#: it, and the heap column then spent about seven times the sweeps it needed, arriving
#: at the field its own coarsest workable rung reaches in a fraction of them.
#:
#: The limit-cycle escape's floor still caps the factor from below (see
#: :data:`_CYCLE_RELAX_FLOOR`): a solve in an orbit is not progressing by any reading,
#: and that floor is what damps it.
_RELAX_LADDER = (1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)
_LADDER_PATIENCE = 5      # sweeps without a new low before stepping down
_LADDER_RECOVER = 10      # new lows in a row before the first attempt at the rung above

#: TEST-ONLY switch, and the only thing that reads it is run_tseep_exit_cycle_test.
#: Set True to run the transient stepper WITHOUT the free-surface rule -- the set is
#: then allowed to be called stable with the head standing above an inactive exit
#: node, which is how the guard measures the limit cycle the rule prevents. Library
#: code never sets it: the rule is the behavior, not an option.
_TSEEP_NO_WET_FACE_RULE = False


def solve_unsaturated(nodes, elements, bc_type, bc_values, kr0=0.001, h0=-1.0,
                      k1_vals=1.0, k2_vals=1.0, angles=0.0,
                      max_iter=400, tol=1e-6, element_types=None,
                      closure_tol=1e-3, vg_a=None, vg_n=None, vg_l=None, model=None,
                      flux_nodal=None):
    """
    Iterative FEM solver for unconfined flow using linear kr frontal function.
    Supports triangular and quadrilateral elements with both linear and quadratic shape functions.

    Convergence is a HYBRID test of three things, all of which must hold on the same
    sweep: the head change (max-norm, an absolute head), the nonlinear residual
    ``rel_closure`` below closure_tol, and an exit-face active set that did not
    change on this sweep.

    THE HEAD-CHANGE GATE IS DATUM-INDEPENDENT: ``max|h_new - h| < tol * head_scale``,
    where ``head_scale`` is the larger of the mesh height and the range of the
    specified heads (floored at 1, as the transient stepper's ``char_h`` is). Both are
    differences, so the gate is the same length wherever the elevation datum is put,
    and one physical problem is held to one head tolerance at every datum.

    WHAT ``rel_closure`` IS: the L1 nodal residual of the UN-RELAXED step over the
    free rows, divided by the inflow — how far the kr matrix lags the head field it
    was assembled from. It is NOT the boundary flow imbalance. The two are different
    quantities and differ by orders of magnitude on the same solve: the closure error
    this function RETURNS is |net inflow - net outflow| at the final consistency
    solve, an absolute flow rate in the model's own flow units. Passing the gate says
    the nonlinear iteration has stopped moving, not that the reported flowrate
    balances to closure_tol. A model with no flow through it has no such ratio to
    form — both sides of it are round-off — and its closure is 0 (see the no-flow
    rule at the gate).

    THE RELAXATION LADDER MOVES ON PROGRESS, NOT ON THE SWEEP NUMBER: it steps down a
    rung after five consecutive sweeps that made none -- a head change setting no new
    low for the rung, or an iterate returning to a recent one -- holds while progress
    continues, and steps back up a rung after a run of ten consecutive progressing
    sweeps, doubled at each later attempt on that rung (see :data:`_RELAX_LADDER`).

    Two escapes act on a solve that is caught in a cycle rather than converging; both
    are inert on a solve that is not (see :data:`_CYCLE_PERSIST` and
    :data:`_SET_REVISIT_SWEEP` for the gates and what was measured on the corpus). A
    third is read only once the sweep budget is spent unconverged: it solves for the
    kr of the one element holding the solve in a loop at every sweep, and the field it
    settles on is then held to the ordinary gate (see :data:`_PIN_MIN_SWING`).

    Parameters:
        element_types : (n_elements,) array indicating:
            3 = 3-node triangle (linear)
            4 = 4-node quadrilateral (bilinear)  
            6 = 6-node triangle (quadratic)
            8 = 8-node quadrilateral (serendipity)
            9 = 9-node quadrilateral (Lagrange)
    Note: Quadratic elements currently use linear/bilinear approximation pending full implementation.
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    y = nodes[:, 1]
    f_ext = (np.zeros(n_nodes) if flux_nodal is None
             else np.asarray(flux_nodal, dtype=float))
    # Water OFFERED at each node by its boundary condition: the threshold the
    # exit-face switch measures an apparent inflow against (see the switch below).
    # Zero wherever no inflow is prescribed, which is every node of a model with no
    # flux BC, so the switch is then the strict-sign rule bit-for-bit.
    q_offered = np.maximum(f_ext, 0.0)

    # Initialize heads
    fixed_heads = bc_values[bc_type == 1]
    has_fixed_head = bool(np.any(bc_type == 1))  # loop-invariant: anchors the solve
    h_free = np.mean(fixed_heads) if len(fixed_heads) > 0 else np.mean(y)
    h = np.where(bc_type == 1, bc_values, np.where(bc_type == 2, y, h_free))

    # Track which exit face nodes are active (saturated)
    exit_face_active = np.ones(n_nodes, dtype=bool)
    exit_face_active[bc_type != 2] = False


    # Build exit-face topology for the active-set update.
    # Corner nodes keep the existing per-node h/q test. For quadratic boundary
    # edges, however, the active set is tracked per edge instead of letting the
    # midside inherit `corner1 OR corner2`. This keeps a quadratic seepage-face
    # side either fully active or fully inactive, so the transition can occur
    # at a corner but not in the middle of a tri6 boundary edge.
    _exit_is_corner, _exit_linear_corners, _exit_quadratic_edges = \
        _exit_face_topology(elements, element_types, bc_type)

    # Store previous iteration values
    h_last = h.copy()

    # Get material properties per element
    if np.isscalar(kr0):
        kr0 = np.full(len(elements), kr0)
    if np.isscalar(h0):
        h0 = np.full(len(elements), h0)
    if vg_a is not None and np.isscalar(vg_a):
        vg_a = np.full(len(elements), vg_a)
    if vg_n is not None and np.isscalar(vg_n):
        vg_n = np.full(len(elements), vg_n)
    if vg_l is not None and np.isscalar(vg_l):
        vg_l = np.full(len(elements), vg_l)
    if model is not None and np.isscalar(model):
        model = np.full(len(elements), model)

    # HEAD-CHANGE GATE: AN ABSOLUTE HEAD, AGAINST A HEAD SCALE BUILT ONLY FROM
    # DIFFERENCES. The sweep-to-sweep head change is compared as a length,
    # max|h_new - h| < tol * head_scale, and the scale is the LARGER of the mesh
    # height and the range of the specified heads. Both are differences, so neither
    # moves when the elevation datum moves, and one physical problem is therefore held
    # to one head tolerance wherever its ground line is drawn. Measured on the heap
    # model of SEEPW-T05, run to the head gate alone: translated up 1000 ft it used to
    # stop 140x further from its answer than it did at the datum it was drawn on, and
    # the two fields differed by 0.73 ft; they now differ by 5.5e-04 ft.
    #
    # This is the scale the TRANSIENT stepper in this file already builds for the same
    # purpose -- ``char_h``, the larger of the head range its series drive and the mesh
    # height, floored at 1, against which its inner Picard's absolute head change is
    # tested. The steady gate now reads the same way, so a model does not change
    # convergence criterion when it is marched instead of solved once.
    #
    # The scale takes the head range as well as the height because the head test is
    # a statement about the head field, and a model can drive far more head across
    # itself than its own mesh is tall (a pressurized confined layer) or far less (a
    # tall column at unit gradient). Taking the larger of the two makes the gate the
    # looser of the two readings rather than the stricter, which is the reading the
    # rest of the solver already uses: `hyst` and the limit-cycle tolerance are both
    # scaled by the mesh height and never by an absolute head.
    ymin, ymax = np.min(y), np.max(y)
    head_span = (float(np.max(fixed_heads) - np.min(fixed_heads))
                 if len(fixed_heads) > 0 else 0.0)
    head_scale = max(float(ymax - ymin), head_span, 1.0)
    eps = head_scale * tol

    print("Starting unsaturated flow iteration...")
    print(f"Convergence tolerance: {eps:.6e}")

    # Track convergence history
    residuals = []
    converged = False
    relax = _RELAX_LADDER[0]   # un-relaxed until the head change argues otherwise
    prev_residual = float('inf')

    # Relaxation ladder state. _rung indexes _RELAX_LADDER; _rung_best is the lowest
    # head change seen SINCE THIS RUNG WAS ENTERED, which is what "progress" is
    # measured against; _up_tries counts the attempts already made at each rung from
    # below, and sets the price of the next one.
    _rung = 0
    _rung_best = float('inf')
    _no_low = 0                # consecutive sweeps with no new low
    _new_lows = 0              # consecutive sweeps with one
    _up_tries = [0] * len(_RELAX_LADDER)

    # Precompute the saturated element matrices, kr sampling operators, and
    # COO index arrays once; each iteration below reduces to a vectorized kr
    # average, a per-element scaling, and one coo_matrix construction.
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals, angles)
    data = _assembly_data(asm, p_nodes=h - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
    _prev_active = exit_face_active.copy()
    _n_stable = 0   # consecutive iterations with an unchanged exit-face set

    # Limit-cycle escape state (Class A). The two ring buffers hold the last
    # _CYCLE_LOOKBACK relaxed iterates and their exit-face sets; _cyc_run counts the
    # consecutive sweeps a period >= 2 repeat has held.
    _cyc_h, _cyc_a = [], []
    _cyc_run = 0
    _cycle_floor = None      # the relaxation floor once the orbit is detected
    _cycle_fired_at = None

    # Set-revisit escape state (Class B). _seen maps each exit-face set visited (bit
    # packed) to the sweep it was last visited on, and _flip_sweep records when each
    # node last changed state, so the nodes that moved during a revisited orbit are
    # exactly those that flipped after the set was last here.
    _seen = {np.packbits(exit_face_active).tobytes(): 0}
    _flip_sweep = np.zeros(n_nodes, dtype=np.int64)   # last sweep each node flipped
    _free_edges = set()      # quadratic exit-face edges the veto no longer governs
    _wet_forced = np.zeros(n_nodes, dtype=int)   # times a wet inactive node was joined
    _revisit_fired_at = None

    # Pinned-element escape state (see _PIN_MIN_SWING): the last _CYCLE_LOOKBACK
    # relaxed iterates, whose kr swing names the element, and what the budget left,
    # which stands if the rescue's gate sweep fails.
    _tail_h = deque(maxlen=_CYCLE_LOOKBACK)
    _pin_tried = False
    _pin_saved = None

    budget = max_iter        # raised once, and only where the floor drops
    iteration = 0
    while iteration < budget:
        iteration += 1

        # Apply boundary conditions: fixed heads plus the active exit-face set
        dir_mask = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
        _require_runtime_dirichlet(dir_mask)
        dir_values = np.where(bc_type == 1, bc_values, y)
        A, b = _dirichlet_system(asm, data, dir_mask, dir_values, neumann=f_ext)

        h_solved = spsolve(A, b)
        h_new = h_solved

        # Relaxation: the rung the head change has argued for (see _RELAX_LADDER),
        # never above the limit-cycle escape's floor.
        relax = _RELAX_LADDER[_rung]
        if _cycle_floor is not None:
            relax = min(relax, _cycle_floor)

        # Apply relaxation
        h_new = relax * h_solved + (1 - relax) * h_last

        # Nodal REACTION, for the exit-face switch below (not used for closure).
        #
        # The switch asks whether the BOUNDARY would have to push water into the
        # domain — unphysical, since a free-draining face can only let water out —
        # so it must test the reaction, not A.h. Only the loads on free rows were
        # ever applied (the Dirichlet rows of b get overwritten), so the reaction is
        # A.h - f_eff:
        #   inactive exit node (free row): the solve enforces (A.h)_i = f_ext_i, so
        #     the reaction is ~0 and `turn_on` falls back to the pressure test, which
        #     is the pre-flux behaviour. Testing raw A.h here would compare the
        #     applied load against ITSELF, always read "inflow", and lock the node
        #     out of the seepage face however positive its pressure grew.
        #   active exit node (Dirichlet row): its load was dropped, so the whole of
        #     A.h is the drain reaction and nothing may be subtracted.
        # Without a flux BC f_eff is zero and this is bit-for-bit the old test.
        f_eff = np.where(dir_mask, 0.0, f_ext)
        q = _coo_matvec(asm, data, h_new) - f_eff

        # Update exit face boundary conditions with hysteresis
        n_active_before = np.sum(exit_face_active)
        hyst = 0.001 * (ymax - ymin)  # Hysteresis threshold

        # Check exit face activation at corner nodes using per-node q
        # (matching SEEP2D). Quadratic exit-face sides are then updated from
        # these corner candidates plus the midside candidate so the whole side
        # remains edge-consistent.
        # SEEP2D (seep2d.f:2137-2147) activates an inactive node on pressure alone and
        # deactivates an active one on inflow. We keep an extra `q <= 0` term on
        # activation: on a free row q is the previous iterate's residual scaled by
        # (1 - relax), and requiring it to be non-positive damps the activation of
        # nodes the last sweep was still pushing water into. It is load-bearing —
        # dropping it costs convergence on the quadratic earth-dam meshes.
        #
        # What that term must NOT see is the applied flux. An inactive node's row is
        # free, so the solve enforces (A.h)_i = f_ext_i there; testing raw A.h would
        # compare the applied load against ITSELF, read "inflow" for any q > 0, and
        # pin the node out of the seepage face forever, however high its pressure
        # climbed. Subtracting f_eff restores the residual the term is meant to test.
        #
        # THE THRESHOLD IS THE WATER OFFERED AT THE NODE, NOT ZERO. Subtracting f_eff
        # is necessary but not sufficient: on a converged free row the residual it
        # leaves is identically zero, so at a flux-loaded exit candidate the SIGN of q
        # is round-off. Measured on the dam-infiltration model at the rate where its
        # surface first ponds: q = +2.4e-22 against a nodal rain load of +2.2e-08 (a
        # ratio of 1e-14) at a node standing at psi = +0.055 m, which therefore never
        # joined the seepage face; a rate higher still and that noise changes sign from
        # sweep to sweep and the active set cycles instead of converging.
        #
        # So the tests compare q against `q_offered` = max(f_ext, 0), which states one
        # physical thing in both directions of the switch: THE BOUNDARY MAY NOT SUPPLY
        # MORE WATER THAN IS PRESCRIBED THERE. A node held at psi = 0 that takes in no
        # more than the rain falling on it is infiltrating that rain and shedding the
        # rest as runoff — precisely what an active exit node over a flux boundary
        # means. Only a reaction ABOVE the offered flux is water the boundary invented,
        # and only that may veto activation or shed an active node. The threshold is
        # thus scaled by the node's own assembled load, in the model's own flow units,
        # as `hyst` is scaled by the domain height and `eps` by height x tol; nothing
        # here is a fixed number.
        #
        # Both sides use the SAME threshold, and that symmetry is what keeps the set
        # from chattering: a node that activates on a near-zero residual cannot then be
        # shed by a reaction the rain it already carries accounts for. And it is f_ext,
        # not f_eff: the offered water is a property of the boundary condition, so the
        # threshold must not move when the node's own row turns Dirichlet.
        corner_candidate = np.zeros(n_nodes, dtype=bool)
        is_corner = (bc_type == 2) & _exit_is_corner
        stay = is_corner & exit_face_active & ~((h_new < y - hyst) | (q > q_offered))
        turn_on = is_corner & ~exit_face_active & (h_new >= y + hyst) & (q <= q_offered)
        corner_candidate[stay | turn_on] = True

        new_exit_face_active = np.zeros(n_nodes, dtype=bool)
        new_exit_face_active[_exit_linear_corners] = corner_candidate[_exit_linear_corners]

        # A specified-head corner on an exit edge (the toe; see _exit_face_topology)
        # is a Dirichlet row already: it counts as saturated in the edge test and is
        # never switched by it.
        _corner_ok = corner_candidate | (bc_type == 1)
        for c1, mid, c2 in _exit_quadratic_edges:
            if exit_face_active[mid]:
                mid_candidate = not (h_new[mid] < y[mid] - hyst
                                     or q[mid] > q_offered[mid])
            else:
                mid_candidate = (h_new[mid] >= y[mid] + hyst
                                 and q[mid] <= q_offered[mid])

            if (c1, mid, c2) in _free_edges:
                # This edge, and only this edge, is updated per node: it was carrying
                # a cycling active set, which is what the all-or-nothing veto below
                # forbids from settling (see _SET_REVISIT_SWEEP). Each node is set
                # from its own test, still by setting True only, so a corner shared
                # with a governed edge keeps whatever that edge gives it.
                if corner_candidate[c1]:
                    new_exit_face_active[c1] = True
                if mid_candidate:
                    new_exit_face_active[mid] = True
                if corner_candidate[c2]:
                    new_exit_face_active[c2] = True
                continue

            edge_active = bool(_corner_ok[c1] and mid_candidate and _corner_ok[c2])
            if edge_active:
                for _nd in (c1, mid, c2):
                    if bc_type[_nd] == 2:
                        new_exit_face_active[_nd] = True

        # Rescue against a SINGULAR full collapse of a quadratic exit face. The
        # edge-based rule above keeps a tri6/quad8 seepage side all-or-nothing,
        # which is what makes the seepage transition land cleanly on a corner on a
        # well-resolved face. But on a COARSE exit face the phreatic exit point
        # lands on a corner shared by the last wet edge and the first dry edge:
        # that corner sits at pressure ~ 0 and reads a borderline (near-zero,
        # often slightly positive) reaction, so it fails its own h/q test AND
        # vetoes the wet edge below it — dropping EVERY exit-face node to inactive.
        # Once empty the face cannot recover, because the strict turn-on test needs
        # h >= y + hyst but the transition corner's equilibrium head is exactly its
        # elevation. With no specified head anywhere the solve is then left with no
        # Dirichlet row at all and spsolve diverges to ~1e15 (issue #51, the tri6
        # thin-domain divergence).
        #
        # Rescue ONLY in that singular configuration — no fixed head AND every exit
        # node inactive — by falling back to per-corner activation (the rule linear
        # tri3 elements always use) for any corner that PASSES its own SEEP2D h/q
        # test. The transition still occurs at a corner (the midside above stays
        # governed by the edge rule and inactive). Gating on `has_fixed_head`
        # is what makes this INERT on every healthy model: any model carrying a
        # specified head (every earth-dam / reservoir case, where the exit face is
        # free to empty transiently and re-fill as the phreatic front settles)
        # keeps a Dirichlet row regardless, so the clause is unreachable there and
        # the edge rule's tuned all-or-nothing toe behaviour is untouched. When
        # the set empties with nothing to rescue, `_require_runtime_dirichlet`
        # (issue #53) reports the true singularity on the next solve.
        face_active = np.any((bc_type == 2) & new_exit_face_active)
        stranded = is_corner & corner_candidate & ~new_exit_face_active
        if not has_fixed_head and not face_active and np.any(stranded):
            new_exit_face_active[stranded] = True

        # A SEEPAGE FACE WITH THE HEAD ABOVE THE GROUND IS NOT A SOLUTION. The set
        # may only be called stable while every inactive exit node stands at or below
        # its elevation. Where it does not, the active-set tests above have finished
        # damping and still leave pressure on a free face, and the face is resolved
        # by hand: the node joins the face, and any quadratic edge it lies on is
        # updated per node from here (`_free_edges`), which is the representation a
        # partly wet edge needs — the all-or-nothing rule sheds the whole edge when
        # its upper corner dries, and the corner below then reads inflow from the
        # midside it just released.
        #
        # Measured on the Johnson Reservoir dam (tri6, 100 divisions): the toe edges
        # cycled on and off four times, the relaxation floor then damped the cycle
        # with the face DRY, and the solve reported converged with 0.66 to 1.77 ft of
        # pressure head on the last 10 ft of the downstream face — 41 to 110 psf of
        # pore pressure that a stress analysis must then hold as tension at a free
        # surface. Read only when the set would otherwise be judged stable, so the
        # ordinary damping is untouched while the set is still moving; capped per
        # node (`_WET_FACE_FORCE_MAX`) so a node the tests keep shedding cannot hold
        # the solve open forever — after that the old behaviour stands and the
        # pressure is reported by the closure check.
        if np.array_equal(new_exit_face_active, _prev_active):
            _wet_free = ((bc_type == 2) & ~new_exit_face_active
                         & (h_new >= y + hyst) & (_wet_forced < _WET_FACE_FORCE_MAX))
            if np.any(_wet_free):
                _wet_forced[_wet_free] += 1
                new_exit_face_active[_wet_free] = True
                _freed_now = {e for e in _exit_quadratic_edges
                              if _wet_free[e[0]] or _wet_free[e[1]] or _wet_free[e[2]]}
                _free_edges = _free_edges | _freed_now
                _wet_idx = np.flatnonzero(_wet_free)
                print(f"Iteration {iteration}: {len(_wet_idx)} inactive exit-face "
                      f"node(s) stand above the ground (max {float(np.max(h_new[_wet_idx] - y[_wet_idx])):.3f} "
                      f"of head) — joined to the seepage face; "
                      f"{len(_freed_now)} quadratic edge(s) updated per node from here")

        newly_active = new_exit_face_active & ~exit_face_active
        h_new[newly_active] = y[newly_active]
        exit_face_active = new_exit_face_active

        n_active_after = np.sum(exit_face_active)

        # The sweep-to-sweep head change, as a HEAD — the quantity `eps` above is a
        # threshold for. Dividing it by max|h| would make it a ratio measured against
        # a number that carries the elevation datum, so the same physical problem
        # would be held to a different real strictness at every datum.
        residual = float(np.max(np.abs(h_new - h)))
        residuals.append(residual)
        _tail_h.append(h_new.copy())


        # Flow-closure probe, measured on the UNRELAXED iterate: q_chk =
        # A(kr(h_solved)) . h_solved. The free rows of A(kr_prev) . h_solved
        # are zero by construction, so this residual isolates the pure
        # nonlinear (kr) lag — and it is immune to the relaxation factor
        # (probing the relaxed blend instead overstates the error by ~1/relax).
        # Closure is the UNSIGNED nodal residual at free nodes (interior +
        # inactive exit face) relative to the inflow; signed in/out sums
        # cancel to ~0 at any head field (zero row sums) and are useless.
        # The free-node residual is A.h - f_ext (zero without flux BCs), so the
        # applied flux must be subtracted or a flux node reads as pure imbalance.
        data_chk = _assembly_data(asm, p_nodes=h_solved - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
        q_chk = _coo_matvec(asm, data_chk, h_solved)
        # One runtime Dirichlet mask for BOTH sides of the ratio. The numerator sums
        # the residual over the free rows, so the denominator must count the flux
        # delivered on those same free rows — masking it by bc_type instead would
        # drop the loads on inactive exit-face nodes and inflate the ratio. Note the
        # mask is rebuilt from the CURRENT active set, which the switch above may
        # have just changed, so it is not necessarily the one used for this solve.
        free_mask = ~((bc_type == 1) | ((bc_type == 2) & exit_face_active))
        inflow_pos = (float(np.sum(q_chk[(bc_type == 1) & (q_chk > 0)]))
                      + _flux_inflow(f_ext, free_mask))
        # NO-FLOW RULE. A model with no flow through it — every specified head at one
        # elevation and no applied flux — makes both sides of this ratio round-off
        # rather than flow rates: the conduction operator has zero row sums, so A.h
        # for a uniform h is identically zero in exact arithmetic, and the ratio is
        # then 0/0. It takes any value at all and no sweep budget can bring it under
        # closure_tol, so the threshold that decides "there is flow here" has to be
        # the model's own flow scale and not an absolute constant. Every term summed
        # into q_chk is one matrix entry times one nodal head, so
        # eps x (number of terms) x max|entry| x max|h| bounds the round-off the sum
        # can carry; an inflow below that is not a flow, and the closure of a model
        # with no flow is 0. That is the reading run_seepage_analysis's reporting
        # layer already takes, where a solution with no flow through the section
        # records closure_fraction = 0.0 rather than a ratio.
        flow_eps = (np.finfo(float).eps * data_chk.size
                    * float(np.max(np.abs(data_chk)))
                    * float(np.max(np.abs(h_solved))))
        rel_closure = (float(np.sum(np.abs(q_chk[free_mask] - f_ext[free_mask]))) / inflow_pos
                       if inflow_pos > max(flow_eps, 1e-30) else 0.0)
        set_stable = bool(np.array_equal(exit_face_active, _prev_active))
        _flip_sweep[exit_face_active != _prev_active] = iteration
        _prev_active = exit_face_active.copy()
        _n_stable = _n_stable + 1 if set_stable else 0

        # ESCAPE 1 (Class A): the head field is in a limit cycle. Compare the relaxed
        # iterate against the last _CYCLE_LOOKBACK sweeps; a RETURN to a sweep 2 or
        # more back — within _CYCLE_TOL x domain height of it, having been further
        # than that away at every shorter lag, with the exit-face set unchanged —
        # held for _CYCLE_PERSIST consecutive sweeps, is an orbit and not a solve
        # still converging. Drop the relaxation floor once, to _CYCLE_RELAX_FLOOR,
        # and grant the sweeps that stepping that finely then takes.
        #
        # Both halves of "returned" are load-bearing. A solve creeping monotonically
        # toward its answer also sits within the tolerance of its recent sweeps once
        # the creep is fine enough, and lands on the SAME reading as an orbit unless
        # the iterate is required to have gone somewhere in between: earth_dam2 after
        # its exit face settles creeps at ~1e-8 of the head range per sweep, which
        # reads as a period-2 repeat and is nothing of the kind. Requiring a shorter
        # lag to be FURTHER away is what separates the two — distance grows with lag
        # under a creep, and collapses at the period under an orbit.
        repeats = False
        if _cycle_floor is None:
            _cyc_eps = _CYCLE_TOL * (ymax - ymin)
            _cyc_away = float('inf')     # the CLOSEST any shorter lag came
            for k in range(1, min(_CYCLE_LOOKBACK, len(_cyc_h)) + 1):
                _cyc_d = float(np.max(np.abs(h_new - _cyc_h[-k])))
                if (k >= 2 and _cyc_d < _cyc_eps and _cyc_away > _cyc_eps
                        and np.array_equal(exit_face_active, _cyc_a[-k])):
                    repeats = True
                    break
                _cyc_away = min(_cyc_away, _cyc_d)
            _cyc_run = _cyc_run + 1 if repeats else 0
            if _cyc_run >= _CYCLE_PERSIST:
                _cycle_floor = _CYCLE_RELAX_FLOOR
                _cycle_fired_at = iteration
                budget = max(budget, iteration + _CYCLE_EXTRA_SWEEPS)
                print(f"Iteration {iteration}: head field is cycling (period <= "
                      f"{_CYCLE_LOOKBACK}, held {_cyc_run} sweeps) — relaxation floor "
                      f"lowered to {_CYCLE_RELAX_FLOOR:g}, sweeps allowed to {budget}")
                _cyc_h, _cyc_a = [], []
            else:
                _cyc_h.append(h_new.copy())
                _cyc_a.append(exit_face_active.copy())
                if len(_cyc_h) > _CYCLE_LOOKBACK:
                    _cyc_h.pop(0)
                    _cyc_a.pop(0)

        # LADDER FEEDBACK. Two things count as progress at this rung, and both have to
        # hold: the head change set a new low for the rung, AND the iterate did not
        # RETURN to a recent one. The second term is the same repeat the class A escape
        # above is built on, and it is here because a shrinking orbit satisfies the
        # first term while going nowhere: johnson_res at relax = 1 swings 2.3 ft a sweep
        # and comes back to where it was, with each swing a little smaller, so its head
        # change sets a new low nearly every sweep. Reading that as progress holds the
        # ladder at 1.0 until the escape fires at sweep 64 and floors the relaxation at
        # 1e-3, which is far finer than the model needs and does not close in the budget.
        # Reading a repeat as the stall it is steps the ladder down while the orbit is
        # still forming, which damps it — the escape then never fires, and the solve
        # closes in tens of sweeps.
        if residual < _rung_best and not repeats:
            _rung_best = residual
            _new_lows += 1
            _no_low = 0
        else:
            _no_low += 1
            _new_lows = 0
        _step = 0
        if _no_low >= _LADDER_PATIENCE and _rung + 1 < len(_RELAX_LADDER):
            _step = 1
        elif (_rung > 0 and _new_lows
                >= _LADDER_RECOVER * (2 ** _up_tries[_rung - 1])):
            _up_tries[_rung - 1] += 1
            _step = -1
        if _step:
            _rung += _step
            _rung_best = float('inf')
            _no_low = _new_lows = 0

        # ESCAPE 2 (Class B): the exit-face active set is in a limit cycle — it has
        # LEFT a set and come back to it, past the sweep by which every converging
        # model in the corpus has finished revisiting. The set it cannot settle on is
        # one the all-or-nothing edge rule forbids, so the veto is lifted from the
        # edges carrying the nodes that moved during the orbit, and only those. The
        # rule still governs every other edge: it fixes a real divergence (issue #51).
        # The transient stepper meets the same cycle within a step and resolves it the
        # same way, per edge (see _resolve_exit_cycle); it pins the resolved nodes for
        # the rest of its step, where here the freed edge keeps testing itself every
        # sweep and the ordinary set-stability gate decides when it has settled.
        if not set_stable and _exit_quadratic_edges and not _free_edges:
            key = np.packbits(exit_face_active).tobytes()
            visited = _seen.get(key)
            if visited is not None and iteration > _SET_REVISIT_SWEEP:
                moved = _flip_sweep > visited
                _free_edges = {e for e in _exit_quadratic_edges
                               if moved[e[0]] or moved[e[1]] or moved[e[2]]}
                if _free_edges:
                    _revisit_fired_at = iteration
                    print(f"Iteration {iteration}: exit-face set returned to the set "
                          f"of sweep {visited} — {len(_free_edges)} cycling edge(s) "
                          f"updated per node from here")
            _seen[key] = iteration

        # Matrix for the next iteration, from the relaxed head field
        data = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)

        # Print detailed iteration info
        if iteration <= 3 or iteration % 5 == 0 or n_active_before != n_active_after:
            print(f"Iteration {iteration}: residual = {residual:.6e}, closure = {rel_closure:.3e}, relax = {relax:.3f}, {n_active_after}/{np.sum(bc_type == 2)} exit face active")

        # Hybrid convergence: head change AND flow closure AND a settled
        # exit-face active set (the flowrate is not meaningful while seepage
        # face nodes are still switching).
        if residual < eps and rel_closure < closure_tol and set_stable:
            print(f"Converged in {iteration} iterations "
                  f"(residual = {residual:.3e}, closure = {rel_closure:.3e}, "
                  f"exit face stable)")
            converged = True
            break

        # Update for next iteration
        h = h_new.copy()
        h_last = h_new.copy()

        # ESCAPE 3: one element holds the solve in a loop (see _PIN_MIN_SWING). Read
        # once, when the budget is spent unconverged. The element's kr is solved for
        # at every sweep until the field settles, and the settled field then gets one
        # ordinary sweep at relax 1: the gate decides, as for any other solve.
        if iteration >= budget and not _pin_tried:
            _pin_tried = True
            kr_kw = dict(kr0=kr0, h0=h0, vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
            pin_mask = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
            dir_values = np.where(bc_type == 1, bc_values, y)
            log_kr = np.log10(np.maximum(
                [_element_kr_avg(asm, t - y, kr_kw) for t in _tail_h], 1e-300))
            swing = log_kr.max(axis=0) - log_kr.min(axis=0)
            h_try, log_try = h_new, [log_kr[-1]]
            for _ in range(_PIN_STEPS):
                d_try = _assembly_data(asm, p_nodes=h_try - y, mode='head', **kr_kw)
                A_try, b_try = _dirichlet_system(asm, d_try, pin_mask, dir_values,
                                                 neumann=f_ext)
                h_try = spsolve(A_try, b_try)
                log_try.append(np.log10(np.maximum(
                    _element_kr_avg(asm, h_try - y, kr_kw), 1e-300)))
            steps = np.ptp(np.array(log_try), axis=0)
            pinned = [int(e) for e in np.argsort(-swing)[:1] if swing[e] > _PIN_MIN_SWING]
            pinned += [int(e) for e in np.argsort(-steps)[:_PIN_STEPS]
                       if steps[e] > _PIN_MIN_SWING and int(e) not in pinned]
            pin_left = max_iter     # the rescue spends at most the budget it follows
            for e in pinned:
                if pin_left <= 0:
                    break
                try:
                    h_pin, pin_ok, n_pin = _pinned_element_solve(
                        asm, kr_kw, y, h_new, e, pin_mask, dir_values, f_ext,
                        eps * _PIN_TOL, pin_left)
                except (np.linalg.LinAlgError, RuntimeError):   # a singular update
                    pin_ok, n_pin = False, 0
                pin_left -= n_pin
                print(f"Iteration {iteration}: element {e + 1}'s kr solved for directly "
                      f"at every sweep — "
                      + (f"settled in {n_pin} sweeps; one ordinary sweep tests it"
                         if pin_ok else f"did not settle ({n_pin} sweeps)"))
                if pin_ok:
                    _pin_saved = (A, data, exit_face_active.copy())
                    h = h_pin.copy()
                    h_last = h_pin.copy()
                    data = _assembly_data(asm, p_nodes=h_pin - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
                    _rung, _cycle_floor = 0, None
                    budget += 1
                    break

    else:
        # The sweeps actually run, which is max_iter unless an escape extended it.
        print(f"Warning: Did not converge in {budget} iterations")
        print("\nConvergence history:")
        for i, r in enumerate(residuals):
            if i % 20 == 0 or i == len(residuals) - 1:
                print(f"  Iteration {i+1}: residual = {r:.6e}")


    if not converged and _pin_saved is not None:
        # The pinned field failed its gate sweep: the answer is the one the budget left.
        A, data, exit_face_active = _pin_saved

    # Final consistency solve: q was computed before the last exit face update,
    # so inactive nodes retain stale reaction forces. Re-solve with the final
    # exit face status and the final (kr-consistent) matrix to get clean q
    # (inactive nodes become free → q ≈ 0), and report the flowrate from it.
    dir_mask = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
    _require_runtime_dirichlet(dir_mask)
    dir_values = np.where(bc_type == 1, bc_values, y)
    A_final, b_final = _dirichlet_system(asm, data, dir_mask, dir_values, neumann=f_ext)
    h_new = spsolve(A_final, b_final)
    q_final = _coo_matvec(asm, data, h_new)

    # Flowrate: the water crossing INTO the domain, which conservation makes equal to
    # the water crossing out of it (the conduction operator has zero row sums, so the
    # reactions and the delivered loads sum to zero). Two entries, and the enumeration
    # has to be complete or the number is not the throughput.
    #
    #   1. THE REACTION AT EVERY DIRICHLET-HELD ROW THAT IS POSITIVE — the specified
    #      heads AND the active exit-face nodes, not the specified heads alone. An
    #      active exit node over a flux boundary is held at p = 0 with its own load
    #      discarded, and what it then draws is the share of that load the soil can
    #      still accept; the rest runs off. That is water entering the domain at the
    #      boundary and it is counted nowhere else. Without a flux load on the exit
    #      face the term is empty: an exit node is active because the switch measured
    #      it draining. (This is the enumeration `_frame_flows` already uses per frame
    #      on the transient side.)
    #   2. THE FLUX DELIVERED ON THE FREE ROWS. "Free" is the runtime Dirichlet
    #      complement, so it includes the INACTIVE exit-face nodes — their rows stay
    #      free, so their loads do enter the solve and are real inflow (rain landing on
    #      the unsaturated part of a seepage face infiltrates).
    #
    # Nothing is counted twice. `_dirichlet_system` seeds the RHS with the loads and
    # then OVERWRITES the held rows, so a load sitting on a held node never entered the
    # solve: `f_eff` zeroes it here and `_flux_inflow` masks it out there. The only
    # water crossing the boundary at a held node is its reaction, and the only water
    # crossing at a free node is its delivered load.
    free_final = ~dir_mask
    f_eff = np.where(dir_mask, 0.0, f_ext)
    react_final = q_final - f_eff
    total_inflow = float(np.sum(react_final[dir_mask & (react_final > 0)]))
    total_inflow += _flux_inflow(f_ext, free_final)

    # Closure report: kr-consistent imbalance at the converged state
    data_chk = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
    q_chk = _coo_matvec(asm, data_chk, h_new)
    react_chk = q_chk - f_eff
    net_inflow = (float(np.sum(react_chk[bc_type == 1]))
                  + _flux_inflow(f_ext, free_final))
    net_outflow = (-float(np.sum(react_chk[dir_mask & (bc_type == 2)]))
                   - float(np.sum(f_ext[free_final & (f_ext < 0)])))
    closure_error = abs(net_inflow - net_outflow)
    print(f"Flow closure check: inflow = {net_inflow:.6e}, outflow = {net_outflow:.6e}, error = {closure_error:.6e}")

    return h_new, A, q_final, total_inflow, exit_face_active, converged, closure_error

def compute_tri6_centroid_pressure(p_nodes, element_nodes):
    """
    Compute pressure at the centroid of a tri6 element using quadratic shape functions.
    
    For GMSH tri6 ordering at centroid (L1=L2=L3=1/3):
    - Corner nodes (0,1,2): N = L*(2*L-1) = 1/3*(2/3-1) = -1/9
    - Edge midpoint nodes (3,4,5): N = 4*L1*L2 = 4*(1/3)*(1/3) = 4/9
    """
    p_elem_nodes = p_nodes[element_nodes[:6]]
    # Shape function values at centroid for GMSH tri6 ordering
    N_corner = -1.0/9.0  # For nodes 0, 1, 2
    N_edge = 4.0/9.0     # For nodes 3, 4, 5
    
    p_centroid = (N_corner * (p_elem_nodes[0] + p_elem_nodes[1] + p_elem_nodes[2]) + 
                  N_edge * (p_elem_nodes[3] + p_elem_nodes[4] + p_elem_nodes[5]))
    return p_centroid



def kr_frontal(p, kr0, h0):
    """
    Fortran-compatible relative permeability function (front model).
    This matches the fkrelf function in the Fortran code exactly.
    """
    if p >= 0.0:
        return 1.0
    elif p > h0:  # when h0 < p < 0
        return kr0 + (1.0 - kr0) * (p - h0) / (-h0)
    else:
        return kr0


# =============================================================================
# Vectorized assembly core
#
# All element stiffness matrices in this module factor as
#     ke = factor(kr_avg) * ke_saturated
# where kr_avg is a weighted average of the frontal kr function sampled at
# element-type-specific points, and ke_saturated depends only on geometry and
# the (anisotropic) conductivity matrix. The helpers below exploit this:
# ke_saturated, the kr sampling matrices, and the COO index arrays are built
# ONCE per solve (batched over all elements of each type with einsum), and
# each assembly (or nonlinear iteration) reduces to a vectorized kr average,
# a per-element scaling, and a single coo_matrix() construction. This replaces
# the per-element Python loops with scipy.lil insertion that previously
# dominated the runtime.
# =============================================================================

from scipy.sparse import coo_matrix


def kr_frontal_vec(p, kr0, h0):
    """Vectorized frontal kr function (same formula as kr_frontal).

    p, kr0, h0 broadcast together; kr0/h0 are per-element, p may be
    (n_elements, n_sample_points)."""
    # h0 = 0 occurs on materials using a different kr model (e.g. vg) whose
    # rows still pass through here before the model mask applies - guard the
    # division so those rows don't emit a spurious divide-by-zero warning
    safe_h0 = np.where(h0 == 0.0, -1.0, h0)
    lin = kr0 + (1.0 - kr0) * (p - safe_h0) / (-safe_h0)
    return np.where(p >= 0.0, 1.0, np.where(p > h0, lin, kr0))


# Per-material unsaturated relative-permeability model codes.
KR_LF = 0    # linear front  (parameters kr0, h0)
KR_VG = 1    # van Genuchten (parameters a = alpha, n)
KR_GARD = 2  # Gardner       (parameters a, n — the SAME two template columns)


def kr_gardner_vec(p, a, n, kr_min=1e-8):
    """Vectorized Gardner (1958) relative permeability.

        kr = 1 / (1 + a * psi^n),   psi = -p (suction, positive when unsaturated)

    This is the POWER form of Gardner, the one carried as a legacy option by
    SEEP/W and Slide — not Gardner's exponential form kr = exp(alpha*psi), which
    is a different function used mainly to linearize Richards' equation for
    analytical work. Saturated (p >= 0) gives kr = 1.

    ``kr_min`` is a numerical guard, not a material property: it keeps the
    conduction matrix nonsingular where a dry element would otherwise contribute
    nothing. At 1e-8 it sits below the conductivity any element in the
    verification corpus carries, and the fields are converged with respect to it
    — a floor of 1e-10 returns the same answer.

    Shares the a/n parameter columns with van Genuchten: the two laws never apply
    to the same material, and the template selects between them with `unsat`.
    """
    a = np.asarray(a, dtype=float)
    n = np.asarray(n, dtype=float)
    psi = np.abs(np.minimum(p, 0.0))          # suction; 0 in the saturated zone
    kr = 1.0 / (1.0 + a * psi ** n)
    kr = np.where(p >= 0.0, 1.0, kr)
    return np.clip(kr, kr_min, 1.0)


def kr_vg_vec(p, vg_a, vg_n, vg_l=None, kr_min=1e-8):
    """Vectorized relative permeability: van Genuchten retention with Mualem
    conductivity (steady-state form).

    Retention is van Genuchten in alpha (``vg_a``) and n (``vg_n``); conductivity
    is Mualem's capillary-bundle integral of it, whose pore-connectivity exponent
    ``vg_l`` scales the tortuosity factor Se^l. Mualem's own value is 0.5, which
    is what ``vg_l`` None means; a fitted curve (HYDRUS, and the PEST calibrations
    SEEP/W ships) carries its own, and negative values are admissible. The
    residual/saturated water contents scale storage, not kr, and a steady-state
    solve carries no storage term. ``p`` is the pressure head (negative in the
    unsaturated zone) and broadcasts with vg_a/vg_n/vg_l. Saturation at p>=0 caps
    the function at 1 and the ``kr_min`` floor holds it away from 0 in the dry
    limit — a numerical guard that keeps the conduction matrix nonsingular, not a
    material property. At 1e-8 the floor lies below the conductivity any element
    in the verification corpus carries, and the fields are converged with respect
    to it (a floor of 1e-10 returns the same answer).

        Se = [1 + (alpha|psi|)^n]^(-m),  m = 1 - 1/n
        kr = Se^l [1 - (1 - Se^(1/m))^m]^2
    """
    n = np.maximum(vg_n, 1.0 + 1e-6)             # guard m = 1 - 1/n away from 0
    m = 1.0 - 1.0 / n
    ah = vg_a * np.abs(np.minimum(p, 0.0))       # |alpha*psi|, 0 in the saturated zone
    Se = (1.0 + ah ** n) ** (-m)
    l = 0.5 if vg_l is None else vg_l
    # The half power is taken as a square root: it is the exponent every material
    # but a fitted one carries, and numpy's sqrt is correctly rounded where a
    # general power is not, so the default curve is exact to the last bit.
    if np.all(l == 0.5):
        se_l = np.sqrt(Se)
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            se_l = np.where(l == 0.5, np.sqrt(Se), Se ** l)
    kr = se_l * (1.0 - (1.0 - Se ** (1.0 / m)) ** m) ** 2
    # Se underflows to 0 in the far dry tail, where the bracket is 0 as well: at
    # l >= 0 that product is 0 and the floor below takes it, but a negative l
    # sends Se^l to infinity and the product to NaN. Both are the same dry limit.
    kr = np.where(Se > 0.0, kr, kr_min)
    kr = np.where(p >= 0.0, 1.0, kr)
    return np.clip(kr, kr_min, 1.0)


def kr_relative_vec(p, kr0, h0, vg_a=None, vg_n=None, vg_l=None, model=None,
                    kr_min=1e-8):
    """Per-element relative permeability dispatching on the unsaturated model.

    ``model`` is a per-element code array (``KR_LF``/``KR_VG``/``KR_GARD``)
    broadcasting with ``p``. With ``model`` None or all linear-front this returns
    exactly ``kr_frontal_vec`` — the linear-front law carries its own floor in
    ``kr0`` and ``kr_min`` never reaches it. ``vg_l`` is the Mualem exponent of
    the van Genuchten branch; None means Mualem's 0.5 for every element."""
    lf = kr_frontal_vec(p, kr0, h0)
    if model is None or not np.any(model):
        return lf
    out = lf
    model = np.asarray(model)
    if np.any(model == KR_VG):
        out = np.where(model == KR_VG, kr_vg_vec(p, vg_a, vg_n, vg_l, kr_min), out)
    if np.any(model == KR_GARD):
        out = np.where(model == KR_GARD, kr_gardner_vec(p, vg_a, vg_n, kr_min), out)
    return out


def kr_relative(p, kr0, h0, vg_a=None, vg_n=None, vg_l=None, model=KR_LF,
                kr_min=1e-8):
    """Scalar relative permeability with model dispatch (linear front, van
    Genuchten–Mualem or Gardner), for the per-edge flow-potential integration."""
    if model == KR_VG:
        return float(kr_vg_vec(float(p), vg_a, vg_n, vg_l, kr_min))
    if model == KR_GARD:
        return float(kr_gardner_vec(float(p), vg_a, vg_n, kr_min))
    return kr_frontal(p, kr0, h0)


# =============================================================================
# Transient storage (specific-storage / moisture-capacity) functions
#
# The transient governing equation is
#     div(kr(psi) K grad h) + Q = S(psi) dh/dt
# with the storage coefficient S split by zone (plan_transient_seep.md §1):
#   saturated   (psi >= 0):  S = Ss              [1/len] elastic specific storage
#   unsaturated (psi < 0):   S = C(psi) = dtheta/dpsi, the specific moisture capacity.
# The split is exclusive — the elastic Ss does NOT carry above the phreatic surface,
# where the retention curve alone supplies the storage.  This is the RS2 / SEEP/W
# convention (both apply m_v only in the saturated zone) and it is what sets the
# transient timing; summing the two would add a compressibility term to every
# unsaturated element that can swamp a flat retention curve.
# These mirror the kr functions: a per-element average of S at the SAME sampling
# points scales a precomputed unit mass matrix, exactly as kr scales the unit
# stiffness matrix — which also smears the psi = 0 switch over one element.
# =============================================================================

# Residual unsaturated-storage floor, as a fraction of the material's own capacity
# scale (Sy/|h0| for the linear front, Sy*alpha for van Genuchten).  The storage
# analog of kr_min: it keeps the lumped mass strictly positive so a shrinking dt
# still pins the head at a node whose retention capacity has decayed to zero.  Four
# orders below the capacity scale, it stores no meaningful water.
S_FLOOR_FRAC = 1e-4


def storage_capacity_vec(p, Ss, Sy, h0, vg_a=None, vg_n=None, model=None):
    """Vectorized specific-storage coefficient S(psi) [1/len], model-dispatched.

    ``p`` is pressure head (negative in the unsaturated zone) and broadcasts with
    the per-element ``Ss``/``Sy``/``h0``/``vg_a``/``vg_n``/``model`` (same broadcast
    convention as ``kr_relative_vec``).

    The elastic ``Ss`` applies in the SATURATED zone ONLY (``psi >= 0``); above the
    phreatic surface the storage is the moisture capacity of the retention curve
    alone, NOT the sum of the two.  That is the convention RS2 and SEEP/W use (both
    apply m_v only where the soil is saturated), and it is what sets the transient
    timing: carrying Ss into the unsaturated zone as well can dominate a flat
    retention curve by more than an order of magnitude and slow the march
    accordingly.

    Linear-front / Gardner: the drainage capacity is ``Sy/|h0|`` inside the band
    ``h0 < psi < 0`` (the same pressure band the kr front spans) and zero below it.
    van Genuchten: ``Sy * dSe/dpsi``, the analytic vG moisture capacity, with
    ``Sy := theta_s - theta_r``.  Gardner reuses the linear band (it is a kr curve,
    not a retention curve).

    The zone switch is sampled per integration point and averaged over the element,
    exactly as ``kr`` is, so an element straddling the phreatic surface carries a
    blended coefficient and the transition is smeared over one element rather than
    snapping node by node.

    A residual floor ``S_FLOOR_FRAC`` x (the material's own capacity scale) applies in
    the unsaturated zone — the storage analog of ``kr_min``, and needed for the same
    reason.  Both retention curves go to zero far from saturation, and a node with
    NO storage is elliptic: its mass term M/dt vanishes, so halving dt no longer
    pins its head, and the exit-face active set can cycle forever without the
    stepper being able to recover by shrinking the step.  Four orders below the
    material's own capacity, the floor never carries physical water; it only keeps
    dt-reduction effective.
    """
    Ss = np.asarray(Ss, dtype=float)
    Sy = np.asarray(Sy, dtype=float)
    h0 = np.asarray(h0, dtype=float)
    safe_h0 = np.where(h0 == 0.0, -1.0, h0)
    in_band = (p < 0.0) & (p > h0)
    band_C = np.where(in_band, Sy / np.abs(safe_h0), 0.0)
    floor_lf = S_FLOOR_FRAC * np.maximum(Ss, Sy / np.abs(safe_h0))
    lf = np.where(p >= 0.0, Ss, np.maximum(band_C, floor_lf))   # LF and Gardner
    if model is None or not np.any(model):
        return lf
    out = lf
    model = np.asarray(model)
    if np.any(model == KR_VG):
        n = np.maximum(vg_n, 1.0 + 1e-6)
        m = 1.0 - 1.0 / n
        ah = vg_a * np.abs(np.minimum(p, 0.0))          # |alpha*psi|, 0 when saturated
        # dSe/dpsi >= 0 (Se rises as psi -> 0); zero in the saturated zone.
        with np.errstate(divide='ignore', invalid='ignore'):
            dSe = m * n * vg_a * ah ** (n - 1.0) * (1.0 + ah ** n) ** (-m - 1.0)
        dSe = np.where(p < 0.0, np.nan_to_num(dSe, nan=0.0, posinf=0.0), 0.0)
        # vG capacity scale = Sy*alpha (dSe/dpsi peaks at O(alpha) near psi ~ -1/alpha)
        floor_vg = S_FLOOR_FRAC * np.maximum(Ss, Sy * np.abs(vg_a))
        vg = np.where(p >= 0.0, Ss, np.maximum(Sy * dSe, floor_vg))
        out = np.where(model == KR_VG, vg, out)
    return out


def _storage_law(Ss, Sy, h0, vg_a=None, vg_n=None):
    """Per-material constants of the storage law ``storage_capacity_vec`` evaluates,
    in the form its integral needs. All arguments broadcast together.

    Linear front / Gardner: the capacity is piecewise constant, ``Ss`` above
    ``psi = 0``, ``s_band`` in the band ``h0 < psi < 0`` and the residual floor
    ``floor_lf`` below it. van Genuchten: ``Ss`` above ``psi = 0`` and
    ``max(Sy dSe/dpsi, floor_vg)`` below it. ``Sy dSe/dpsi`` rises from zero at
    saturation to a peak and decays again, so it crosses the floor twice, at
    ``p1`` (near saturation) and ``p2`` (deep in the dry range); between them the
    retention curve holds the water, outside them the floor. The crossings are
    found once per material by bisection on ``log(alpha|psi|)``. A curve whose
    peak never clears the floor (or a non-positive alpha) carries the floor alone,
    flagged by ``vg_has_band`` False."""
    Ss = np.asarray(Ss, dtype=float)
    Sy = np.asarray(Sy, dtype=float)
    h0 = np.asarray(h0, dtype=float)
    safe_h0 = np.where(h0 == 0.0, -1.0, h0)
    band_C = Sy / np.abs(safe_h0)
    floor_lf = S_FLOOR_FRAC * np.maximum(Ss, band_C)
    law = {'Ss': Ss, 'Sy': Sy, 'h0': h0, 'floor_lf': floor_lf,
           's_band': np.maximum(band_C, floor_lf)}
    if vg_a is None or vg_n is None:
        return law
    a = np.asarray(vg_a, dtype=float)
    n = np.maximum(np.asarray(vg_n, dtype=float), 1.0 + 1e-6)
    m = 1.0 - 1.0 / n
    floor_vg = S_FLOOR_FRAC * np.maximum(Ss, Sy * np.abs(a))
    Sy_b, a_b, n_b, m_b, f_b = np.broadcast_arrays(Sy, a, n, m, floor_vg)
    with np.errstate(divide='ignore', invalid='ignore'):
        log_c = np.log(Sy_b * m_b * n_b * a_b)

        def log_g(ls):                     # log of Sy dSe/dpsi at s = exp(ls)
            return log_c + (n_b - 1.0) * ls - (m_b + 1.0) * np.logaddexp(0.0, n_b * ls)
        log_f = np.log(f_b)
        ls_pk = np.log(m_b) / n_b          # the capacity peaks at s^n = m
        has = (a_b > 0) & (Sy_b > 0) & np.isfinite(log_c) & (log_g(ls_pk) > log_f)
        lo1 = np.full(a_b.shape, -700.0)
        hi1 = np.where(has, ls_pk, 0.0)
        lo2 = np.where(has, ls_pk, 0.0)
        hi2 = np.where(has, ls_pk + 700.0 / n_b, 1.0)
        low_edge = has & (log_g(lo1) >= log_f)   # crossing below s = e^-700: s1 = 0
        for _ in range(100):
            mid = 0.5 * (lo1 + hi1)
            up = log_g(mid) >= log_f
            hi1 = np.where(up, mid, hi1)
            lo1 = np.where(up, lo1, mid)
            mid = 0.5 * (lo2 + hi2)
            up = log_g(mid) >= log_f
            lo2 = np.where(up, mid, lo2)
            hi2 = np.where(up, hi2, mid)
        s1 = np.where(low_edge, 0.0, np.exp(0.5 * (lo1 + hi1)))
        s2 = np.exp(0.5 * (lo2 + hi2))
        a_safe = np.where(a_b > 0, a_b, 1.0)
        p1 = np.where(has, -s1 / a_safe, 0.0)
        p2 = np.where(has, -s2 / a_safe, 0.0)
        Se1 = (1.0 + s1 ** n_b) ** (-m_b)
        Se2 = np.exp(-m_b * np.logaddexp(0.0, n_b * np.log(np.maximum(s2, 1e-300))))
    law.update({'vg_a': a, 'vg_n': n, 'vg_m': m, 'floor_vg': floor_vg,
                'vg_has_band': has, 'vg_p1': p1, 'vg_p2': p2,
                'vg_Se1': np.where(has, Se1, 1.0), 'vg_Se2': np.where(has, Se2, 1.0)})
    return law


def _vg_point(p, law, want_cap=True):
    """Phi and (``want_cap``) the capacity of the van Genuchten law at pressure head
    p, both with the residual floor, from one evaluation of (alpha |psi|)^n."""
    a, n, m, Sy = law['vg_a'], law['vg_n'], law['vg_m'], law['Sy']
    s = a * np.maximum(-p, 0.0)                    # alpha |psi| (0 where saturated)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        sn = s ** n
        base = 1.0 + sn
        Se = base ** (-m)
    fl, p1, p2 = law['floor_vg'], law['vg_p1'], law['vg_p2']
    off = fl * p1
    band = off + Sy * (Se - law['vg_Se1'])
    dry = off + Sy * (law['vg_Se2'] - law['vg_Se1']) + fl * (p - p2)
    uns = np.where(p >= p1, fl * p, np.where(p >= p2, band, dry))
    uns = np.where(law['vg_has_band'], uns, fl * p)
    phi = np.where(p >= 0.0, law['Ss'] * p, uns)
    if not want_cap:
        return phi, None
    # Sy dSe/dpsi = Sy m n alpha s^(n-1) (1 + s^n)^(-m-1) = Sy m n alpha (s^n/s) Se/base
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        dSe = m * n * a * (sn / np.where(s > 0.0, s, 1.0)) * Se / base
    dSe = np.where(s > 0.0, np.nan_to_num(dSe, nan=0.0, posinf=0.0), 0.0)
    cap = np.where(p >= 0.0, law['Ss'], np.maximum(Sy * dSe, fl))
    return phi, cap


def _phi_lf(p, law):
    """Integral from 0 to p of the linear-front / Gardner capacity, floor included."""
    h0, s_band = law['h0'], law['s_band']
    return np.where(p >= 0.0, law['Ss'] * p,
                    np.where(p > h0, s_band * p,
                             s_band * h0 + law['floor_lf'] * (p - h0)))


def _phi_vg(p, law):
    """Integral from 0 to p of the van Genuchten capacity, floor included."""
    return _vg_point(p, law, want_cap=False)[0]


def _cap_vg(p, law):
    """The van Genuchten branch of ``storage_capacity_vec`` (floor included)."""
    return _vg_point(p, law)[1]


def storage_potential_vec(p, Ss, Sy, h0, vg_a=None, vg_n=None, model=None):
    """Storage potential Phi(psi) = integral_0^psi S(psi') dpsi' [-], vectorized.

    ``Phi`` is the stored water per unit volume, measured from saturation:
    ``Phi(psi_new) - Phi(psi_old)`` is the water a unit volume of soil takes up
    between the two pressure heads. It is the exact integral of
    ``storage_capacity_vec`` -- the same zone-split law (elastic Ss below the
    phreatic surface, the retention capacity above it) with the same residual
    floor -- so the stored water the transient solver conserves and the capacity
    it reports are one law. Arguments broadcast as in ``storage_capacity_vec``.
    """
    law = _storage_law(Ss, Sy, h0)
    lf = _phi_lf(p, law)
    if model is None or not np.any(model):
        return lf
    model = np.asarray(model)
    out = lf
    if np.any(model == KR_VG):
        vlaw = _storage_law(Ss, Sy, h0, vg_a, vg_n)
        out = np.where(model == KR_VG, _phi_vg(p, vlaw), out)
    return out


def _idx_or_none(arr, idx):
    """Index an optional per-element array, or pass through None (no vG model)."""
    return None if arr is None else arr[idx]


def _element_kmats(n_el, k1_vals, k2_vals, angles, flow=False):
    """Batched (n_el, 2, 2) conductivity matrices K = R^T diag(k1,k2) R.

    flow=True returns K/det(K) (the stream-function coefficient matrix)."""
    k1 = np.broadcast_to(np.asarray(k1_vals, dtype=float), (n_el,))
    k2 = np.broadcast_to(np.asarray(k2_vals, dtype=float), (n_el,))
    th = np.radians(np.broadcast_to(np.asarray(angles, dtype=float), (n_el,)))
    c, sn = np.cos(th), np.sin(th)
    K = np.empty((n_el, 2, 2))
    K[:, 0, 0] = k1 * c * c + k2 * sn * sn
    K[:, 1, 1] = k1 * sn * sn + k2 * c * c
    K[:, 0, 1] = K[:, 1, 0] = (k1 - k2) * c * sn
    if flow:
        det = K[:, 0, 0] * K[:, 1, 1] - K[:, 0, 1] * K[:, 1, 0]
        K = K / det[:, None, None]
    return K


def _kr_sampling(et):
    """(N, w) for the kr average of element type et: N is (n_pts, nn) shape
    functions at the sampling points, w the quadrature weights. Identical
    points/weights to the per-element *_stiffness_matrix_kr functions."""
    if et == 3:
        # 7-point symmetric triangle rule (degree 5), linear N = (L1, L2, L3)
        a1, b1 = 0.059715871789770, 0.470142064105115
        a2, b2 = 0.797426985353087, 0.101286507323456
        w0, w1, w2 = 0.1125, 0.066197076394253, 0.062969590272414
        pts = [(1/3, 1/3, 1/3, w0),
               (a1, b1, b1, w1), (b1, a1, b1, w1), (b1, b1, a1, w1),
               (a2, b2, b2, w2), (b2, a2, b2, w2), (b2, b2, a2, w2)]
        N = np.array([[L1, L2, L3] for L1, L2, L3, _ in pts])
        w = np.array([p[3] for p in pts])
    elif et == 6:
        pts = [(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)]
        N = np.array([[L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                       4*L1*L2, 4*L2*L3, 4*L3*L1] for L1, L2, L3 in pts])
        w = np.full(3, 1/3)
    elif et == 4:
        # 4x4 rule (matching SEEP2D's qdflow kr sampling)
        p1 = [-0.86113631, -0.33998104, 0.33998104, 0.86113631]
        w1 = [0.34785485, 0.65214516, 0.65214516, 0.34785485]
        N, w = [], []
        for i, xi in enumerate(p1):
            for j, eta in enumerate(p1):
                N.append([0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta),
                          0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta)])
                w.append(w1[i] * w1[j])
        N, w = np.array(N), np.array(w)
    elif et in (8, 9):
        p1 = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        w1 = [5/9, 8/9, 5/9]
        N, w = [], []
        for i, xi in enumerate(p1):
            for j, eta in enumerate(p1):
                if et == 8:
                    N.append([0.25*(1-xi)*(1-eta)*(-xi-eta-1),
                              0.25*(1+xi)*(1-eta)*(xi-eta-1),
                              0.25*(1+xi)*(1+eta)*(xi+eta-1),
                              0.25*(1-xi)*(1+eta)*(-xi+eta-1),
                              0.5*(1-xi*xi)*(1-eta), 0.5*(1+xi)*(1-eta*eta),
                              0.5*(1-xi*xi)*(1+eta), 0.5*(1-xi)*(1-eta*eta)])
                else:
                    N.append([0.25*xi*(xi-1)*eta*(eta-1), 0.25*xi*(xi+1)*eta*(eta-1),
                              0.25*xi*(xi+1)*eta*(eta+1), 0.25*xi*(xi-1)*eta*(eta+1),
                              0.5*(1-xi*xi)*eta*(eta-1), 0.5*xi*(xi+1)*(1-eta*eta),
                              0.5*(1-xi*xi)*eta*(eta+1), 0.5*xi*(xi-1)*(1-eta*eta),
                              (1-xi*xi)*(1-eta*eta)])
                w.append(w1[i] * w1[j])
        N, w = np.array(N), np.array(w)
    else:
        raise ValueError(f"Unknown element type {et}")
    return N, w


def _quad_gauss(et):
    """[(xi, eta, w)] for quad element type et: 2x2 for quad4, 3x3 for quad8/quad9.

    One definition, shared by the assembly and by the gradient/velocity recovery,
    so an element is differentiated at the points it was integrated at."""
    if et == 4:
        g = 1/np.sqrt(3)
        return [(-g, -g, 1.0), (g, -g, 1.0), (g, g, 1.0), (-g, g, 1.0)]
    p1 = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    w1 = [5/9, 8/9, 5/9]
    return [(xi, eta, w1[i]*w1[j]) for i, xi in enumerate(p1)
            for j, eta in enumerate(p1)]


def _quad_dshape(et, xi, eta):
    """(dN_dxi, dN_deta) for quad element type et at natural point (xi, eta)."""
    if et == 4:
        dxi = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
        deta = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)])
    elif et == 8:
        dxi = np.array([
            -0.25*(1-eta)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta),
            0.25*(1-eta)*(xi-eta-1) + 0.25*(1+xi)*(1-eta),
            0.25*(1+eta)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),
            -0.25*(1+eta)*(-xi+eta-1) - 0.25*(1-xi)*(1+eta),
            -xi*(1-eta), 0.5*(1-eta*eta), -xi*(1+eta), -0.5*(1-eta*eta)])
        deta = np.array([
            -0.25*(1-xi)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta),
            -0.25*(1+xi)*(xi-eta-1) - 0.25*(1+xi)*(1-eta),
            0.25*(1+xi)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),
            0.25*(1-xi)*(-xi+eta-1) + 0.25*(1-xi)*(1+eta),
            -0.5*(1-xi*xi), -eta*(1+xi), 0.5*(1-xi*xi), -eta*(1-xi)])
    else:  # 9
        dxi = np.array([
            0.25*(2*xi-1)*eta*(eta-1), 0.25*(2*xi+1)*eta*(eta-1),
            0.25*(2*xi+1)*eta*(eta+1), 0.25*(2*xi-1)*eta*(eta+1),
            -xi*eta*(eta-1), 0.5*(2*xi+1)*(1-eta*eta),
            -xi*eta*(eta+1), 0.5*(2*xi-1)*(1-eta*eta), -2*xi*(1-eta*eta)])
        deta = np.array([
            0.25*xi*(xi-1)*(2*eta-1), 0.25*xi*(xi+1)*(2*eta-1),
            0.25*xi*(xi+1)*(2*eta+1), 0.25*xi*(xi-1)*(2*eta+1),
            0.5*(1-xi*xi)*(2*eta-1), -eta*xi*(xi+1),
            0.5*(1-xi*xi)*(2*eta+1), -eta*xi*(xi-1), -2*eta*(1-xi*xi)])
    return dxi, deta


def _batched_ke_sat(et, coords, Kmats):
    """Batched saturated element stiffness matrices for one element type.

    coords: (n_e, nn, 2); Kmats: (n_e, 2, 2). Returns (n_e, nn, nn).
    Reproduces the per-element *_stiffness_matrix functions, including the
    degenerate-element guards (zero contribution)."""
    n_e = coords.shape[0]
    if et == 3:
        x, y = coords[:, :, 0], coords[:, :, 1]
        area = 0.5 * np.abs((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0])
                            - (x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
        ok = area > 1e-14
        a_safe = np.where(ok, area, 1.0)
        beta = np.stack([y[:, 1]-y[:, 2], y[:, 2]-y[:, 0], y[:, 0]-y[:, 1]], axis=1)
        gamma = np.stack([x[:, 2]-x[:, 1], x[:, 0]-x[:, 2], x[:, 1]-x[:, 0]], axis=1)
        grad = np.stack([beta, gamma], axis=1) / (2 * a_safe)[:, None, None]
        ke = area[:, None, None] * np.einsum('eai,eab,ebj->eij', grad, Kmats, grad)
        ke[~ok] = 0.0
        return ke
    if et == 6:
        x, y = coords[:, :, 0], coords[:, :, 1]
        detJ = (x[:, 0]-x[:, 2])*(y[:, 1]-y[:, 2]) - (x[:, 1]-x[:, 2])*(y[:, 0]-y[:, 2])
        ok = np.abs(detJ) > 1e-10
        area = 0.5 * np.abs(detJ)
        a_safe = np.where(ok, area, 1.0)
        # dL_i/d(x,y): (n_e, 2, 3)
        dL = np.empty((n_e, 2, 3))
        dL[:, 0, 0] = y[:, 1]-y[:, 2]
        dL[:, 1, 0] = x[:, 2]-x[:, 1]
        dL[:, 0, 1] = y[:, 2]-y[:, 0]
        dL[:, 1, 1] = x[:, 0]-x[:, 2]
        dL[:, 0, 2] = y[:, 0]-y[:, 1]
        dL[:, 1, 2] = x[:, 1]-x[:, 0]
        dL /= (2 * a_safe)[:, None, None]
        pts = [(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)]
        ke = np.zeros((n_e, 6, 6))
        for (L1, L2, L3), w in zip(pts, [1/3]*3):
            dN_dL = np.array([[4*L1-1, 0, 0, 4*L2, 0, 4*L3],
                              [0, 4*L2-1, 0, 4*L1, 4*L3, 0],
                              [0, 0, 4*L3-1, 0, 4*L2, 4*L1]])  # (3, 6)
            gradN = np.einsum('exl,ln->exn', dL, dN_dL)         # (n_e, 2, 6)
            ke += (w * area)[:, None, None] * np.einsum(
                'eai,eab,ebj->eij', gradN, Kmats, gradN)
        ke[~ok] = 0.0
        return ke
    # quads: 2x2 rule for quad4, 3x3 for quad8/quad9 (same as per-element code)
    gps = _quad_gauss(et)
    nn = coords.shape[1]
    ke = np.zeros((n_e, nn, nn))
    for xi, eta, w in gps:
        dxi, deta = _quad_dshape(et, xi, eta)
        J00 = coords[:, :, 0] @ dxi
        J01 = coords[:, :, 1] @ dxi
        J10 = coords[:, :, 0] @ deta
        J11 = coords[:, :, 1] @ deta
        detJ = J00*J11 - J01*J10
        ok = detJ > 0
        d_safe = np.where(ok, detJ, 1.0)
        # Jinv rows applied to (dxi, deta)
        dN_dx = (J11[:, None]*dxi[None, :] - J01[:, None]*deta[None, :]) / d_safe[:, None]
        dN_dy = (-J10[:, None]*dxi[None, :] + J00[:, None]*deta[None, :]) / d_safe[:, None]
        gradN = np.stack([dN_dx, dN_dy], axis=1)               # (n_e, 2, nn)
        contrib = (w * detJ)[:, None, None] * np.einsum(
            'eai,eab,ebj->eij', gradN, Kmats, gradN)
        contrib[~ok] = 0.0
        ke += contrib
    return ke


def _build_assembly(nodes, elements, element_types, k1_vals, k2_vals, angles,
                    flow=False):
    """Precompute everything reusable across assemblies of one mesh.

    Returns a dict with per-element-type groups (saturated ke stacks, kr
    sampling matrices, connectivity) plus the concatenated COO row/col index
    arrays. Build once per solve; each assembly is then a kr scaling plus a
    single coo_matrix construction."""
    elements = np.asarray(elements)
    element_types = np.asarray(element_types)
    n_el = len(elements)
    Kmats = _element_kmats(n_el, k1_vals, k2_vals, angles, flow=flow)
    groups, rows_all, cols_all = [], [], []
    for et in np.unique(element_types):
        idx = np.where(element_types == et)[0]
        nn = int(et)
        conn = elements[idx][:, :nn]
        ke = _batched_ke_sat(int(et), nodes[conn], Kmats[idx])
        N_kr, w_kr = _kr_sampling(int(et))
        groups.append({'et': int(et), 'idx': idx, 'conn': conn, 'ke': ke,
                       'N_kr': N_kr, 'w_kr': w_kr, 'wsum': w_kr.sum()})
        rows_all.append(np.repeat(conn, nn, axis=1).ravel())
        cols_all.append(np.tile(conn, (1, nn)).ravel())
    return {'groups': groups, 'n_nodes': len(nodes),
            'rows': np.concatenate(rows_all), 'cols': np.concatenate(cols_all)}


def _assembly_data(asm, p_nodes=None, kr0=None, h0=None, mode='head',
                   vg_a=None, vg_n=None, vg_l=None, model=None):
    """Concatenated COO data for the global matrix.

    p_nodes=None -> saturated assembly. Otherwise each element's ke is scaled
    by factor(kr_avg) with kr averaged at the type-specific sampling points
    (mode='head': kr_avg; mode='stream': 1/kr_avg, guarded). ``model``/``vg_a``/
    ``vg_n``/``vg_l`` (per-element) select the unsaturated model; with ``model``
    None the kr is the linear-front function exactly as before."""
    parts = []
    for g in asm['groups']:
        if p_nodes is None:
            parts.append(g['ke'].ravel())
            continue
        idx = g['idx']
        p_gp = p_nodes[g['conn']] @ g['N_kr'].T              # (n_e, n_pts)
        kr = kr_relative_vec(
            p_gp, kr0[idx][:, None], h0[idx][:, None],
            None if vg_a is None else vg_a[idx][:, None],
            None if vg_n is None else vg_n[idx][:, None],
            None if vg_l is None else vg_l[idx][:, None],
            None if model is None else model[idx][:, None])
        kr_avg = (kr @ g['w_kr']) / g['wsum']
        if mode == 'head':
            factor = kr_avg
        else:
            factor = np.where(kr_avg > 1e-12, 1.0 / np.maximum(kr_avg, 1e-300), 1e10)
        parts.append((factor[:, None, None] * g['ke']).ravel())
    return np.concatenate(parts)


def _csr_from_data(asm, data):
    return coo_matrix((data, (asm['rows'], asm['cols'])),
                      shape=(asm['n_nodes'], asm['n_nodes'])).tocsr()


def _dirichlet_system(asm, data, dir_mask, dir_values, neumann=None):
    """BC-applied system (A, b): Dirichlet rows replaced by identity.

    Equivalent to the previous LIL row-zeroing, built directly from the COO
    arrays (entries in Dirichlet rows are dropped, then unit diagonals are
    appended). `neumann` is the (n_nodes,) vector of consistent nodal loads
    from specified-flux BCs; it seeds b, and Dirichlet rows then overwrite it —
    a node cannot be both."""
    n = asm['n_nodes']
    keep = ~dir_mask[asm['rows']]
    di = np.where(dir_mask)[0]
    rows = np.concatenate([asm['rows'][keep], di])
    cols = np.concatenate([asm['cols'][keep], di])
    vals = np.concatenate([data[keep], np.ones(len(di))])
    A = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    b = np.zeros(n) if neumann is None else np.asarray(neumann, dtype=float).copy()
    b[di] = dir_values[di]
    return A, b


def _coo_matvec(asm, data, x):
    """y = A @ x directly from the COO arrays (no tocsr / duplicate handling)."""
    return np.bincount(asm['rows'], weights=data * x[asm['cols']],
                       minlength=asm['n_nodes'])


# =============================================================================
# Transient storage (batched, build-once)
#
# The storage term of a time step is the change in each element's stored water,
# the integral of Phi(psi) (storage_potential_vec) over the element with psi carried
# linearly across sub-triangles of it, between the head at the start of the step
# and the head at its end -- the mass-conservative mixed form of Celia, Bouloutas
# and Zarba (1990). The integral is exact for the linear-front / Gardner law on
# triangles and a sum over storage points otherwise. The change splits exactly into
# one term per node, each driven by that node's own head change, so the storage
# matrix stays diagonal; wherever an element stays on one linear branch of the law
# it is exactly the HRZ-lumped m_i * S * dpsi_i (see _storage_residual). Everything
# geometric -- the HRZ lumped weights, the sub-triangles, the storage points and
# their weights -- is built once per solve. HRZ lumping is used rather than row-sum
# lumping because row-sum gives a zero corner weight on quadratic elements; the two
# coincide on tri3/quad4.
# =============================================================================


def _shape_N(et, xi, eta):
    """Shape-function values N(xi, eta) for element type et (natural coords).

    Triangles use area coordinates (xi, eta, 1-xi-eta); quads use (xi, eta) in
    [-1, 1].  Matches the shape functions the stiffness/kr code uses."""
    if et == 3:
        L1, L2, L3 = xi, eta, 1.0 - xi - eta
        return np.array([L1, L2, L3])
    if et == 6:
        L1, L2, L3 = xi, eta, 1.0 - xi - eta
        return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                         4*L1*L2, 4*L2*L3, 4*L3*L1])
    if et == 4:
        return 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                                (1+xi)*(1+eta), (1-xi)*(1+eta)])
    if et == 8:
        return np.array([0.25*(1-xi)*(1-eta)*(-xi-eta-1),
                         0.25*(1+xi)*(1-eta)*(xi-eta-1),
                         0.25*(1+xi)*(1+eta)*(xi+eta-1),
                         0.25*(1-xi)*(1+eta)*(-xi+eta-1),
                         0.5*(1-xi*xi)*(1-eta), 0.5*(1+xi)*(1-eta*eta),
                         0.5*(1-xi*xi)*(1+eta), 0.5*(1-xi)*(1-eta*eta)])
    if et == 9:
        return np.array([0.25*xi*(xi-1)*eta*(eta-1), 0.25*xi*(xi+1)*eta*(eta-1),
                         0.25*xi*(xi+1)*eta*(eta+1), 0.25*xi*(xi-1)*eta*(eta+1),
                         0.5*(1-xi*xi)*eta*(eta-1), 0.5*xi*(xi+1)*(1-eta*eta),
                         0.5*(1-xi*xi)*eta*(eta+1), 0.5*xi*(xi-1)*(1-eta*eta),
                         (1-xi*xi)*(1-eta*eta)])
    raise ValueError(f"Unknown element type {et}")


def _mass_quadrature(et):
    """(points, weights) for integrating the consistent mass matrix of type et.

    Each point is (xi, eta); weights are reference-element weights.  Triangles
    use a degree-5 7-point rule (exact for the degree-4 N_i N_j of tri6); quads
    use 2x2 (quad4) or 3x3 (quad8/quad9) Gauss."""
    if et in (3, 6):
        a1, b1 = 0.059715871789770, 0.470142064105115
        a2, b2 = 0.797426985353087, 0.101286507323456
        w0, w1, w2 = 0.1125, 0.066197076394253, 0.062969590272414
        pts = [(1/3, 1/3, w0),
               (a1, b1, w1), (b1, a1, w1), (b1, b1, w1),
               (a2, b2, w2), (b2, a2, w2), (b2, b2, w2)]
        return [(L1, L2) for L1, L2, _ in pts], [p[2] for p in pts]
    if et == 4:
        g = 1.0 / np.sqrt(3.0)
        return [(-g, -g), (g, -g), (g, g), (-g, g)], [1.0, 1.0, 1.0, 1.0]
    if et in (8, 9):
        p1 = [-np.sqrt(3/5), 0.0, np.sqrt(3/5)]
        w1 = [5/9, 8/9, 5/9]
        pts, wts = [], []
        for i, x in enumerate(p1):
            for j, y in enumerate(p1):
                pts.append((x, y))
                wts.append(w1[i] * w1[j])
        return pts, wts
    raise ValueError(f"Unknown element type {et}")


def _batched_me(et, coords):
    """Batched unit (density-1) consistent element mass matrices for one type.

    coords: (n_e, nn, 2).  Returns (n_e, nn, nn) = integral(N_i N_j dOmega),
    with the degenerate-element guard (zero contribution) matching the stiffness
    builder."""
    n_e, nn = coords.shape[0], coords.shape[1]
    pts, wts = _mass_quadrature(et)
    me = np.zeros((n_e, nn, nn))
    for (xi, eta), w in zip(pts, wts):
        N = _shape_N(et, xi, eta)                       # (nn,)
        if et in (3, 6):
            x, y = coords[:, :, 0], coords[:, :, 1]
            # constant Jacobian (straight-sided): detJ = 2*Area, from the 3 corners
            detJ = np.abs((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0])
                          - (x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
        else:
            dxi, deta = _quad_dshape(et, xi, eta)
            J00 = coords[:, :, 0] @ dxi
            J01 = coords[:, :, 1] @ dxi
            J10 = coords[:, :, 0] @ deta
            J11 = coords[:, :, 1] @ deta
            detJ = np.abs(J00*J11 - J01*J10)
        me += (w * detJ)[:, None, None] * (N[:, None] * N[None, :])[None, :, :]
    ok = np.abs(me).sum(axis=(1, 2)) > 1e-300
    me[~ok] = 0.0
    return me


def _build_mass_assembly(nodes, elements, element_types):
    """Precompute the HRZ-lumped unit mass of every element (parallels
    ``_build_assembly``, same group order and connectivity).

    Returns per-element-type groups carrying the unit consistent mass ``me`` and
    its HRZ-lumped diagonal ``me_lump`` (n_e, nn): the weight with which each node
    of an element carries that element's storage."""
    elements = np.asarray(elements)
    element_types = np.asarray(element_types)
    groups = []
    for et in np.unique(element_types):
        idx = np.where(element_types == et)[0]
        nn = int(et)
        conn = elements[idx][:, :nn]
        me = _batched_me(int(et), nodes[conn])
        diag = np.einsum('eii->ei', me)                 # (n_e, nn)
        total = me.sum(axis=(1, 2))                     # (n_e,)
        dsum = diag.sum(axis=1)
        # HRZ diagonal scaling: lumped_i = diag_i * total / sum(diag). Identical to
        # row-sum on linear elements, strictly positive on quadratics.
        scale = np.where(dsum > 1e-300, total / np.where(dsum > 1e-300, dsum, 1.0), 0.0)
        me_lump = diag * scale[:, None]                 # (n_e, nn)
        groups.append({'et': int(et), 'idx': idx, 'conn': conn,
                       'me': me, 'me_lump': me_lump})
    return {'groups': groups, 'n_nodes': len(nodes)}


#: Storage points per sub-triangle where the drainage water is a point sum: each
#: sub-triangle is cut into k x k equal parts, each sampled by the 7-point rule. The
#: van Genuchten law is smooth, and one 7-point rule per sub-triangle holds its
#: water to well inside the iteration tolerance; the linear-front law on a
#: quadrilateral has kinks, so its sub-triangles are cut finer. (The linear-front
#: law on triangles is integrated exactly; see ``_lf_tri_water``.)
STORAGE_SUBDIV_VG = 1
STORAGE_SUBDIV_LF = 2

# 7-point (degree-5) triangle rule in barycentric coordinates.
_TRI7_L = np.array([[1/3, 1/3, 1/3],
                    [0.059715871789770, 0.470142064105115, 0.470142064105115],
                    [0.470142064105115, 0.059715871789770, 0.470142064105115],
                    [0.470142064105115, 0.470142064105115, 0.059715871789770],
                    [0.797426985353087, 0.101286507323456, 0.101286507323456],
                    [0.101286507323456, 0.797426985353087, 0.101286507323456],
                    [0.101286507323456, 0.101286507323456, 0.797426985353087]])
_TRI7_W = np.array([0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                    0.125939180544827, 0.125939180544827, 0.125939180544827])
_TRI7_W = _TRI7_W / _TRI7_W.sum()


def _subtri_rule(k):
    """Barycentric points (n, 3) and weights (n,) summing to 1: the 7-point rule on
    each of the k*k equal sub-triangles of a triangle."""
    L, W = [], []
    for j in range(k):
        for i in range(k - j):
            cells = [((i, j), (i + 1, j), (i, j + 1))]
            if i + j < k - 1:
                cells.append(((i + 1, j), (i + 1, j + 1), (i, j + 1)))
            for v in cells:
                V = np.array([[a / k, b / k, 1.0 - (a + b) / k] for a, b in v])
                for lam, w in zip(_TRI7_L, _TRI7_W):
                    L.append(lam @ V)
                    W.append(w / (k * k))
    return np.array(L), np.array(W)


def _storage_cells(et):
    """Storage sub-triangles of the parent element, across each of which the
    pressure head is carried linearly.

    Returns ``(P, tris)``: ``P`` (n_pts, nn) the non-negative weights that give the
    head at each sub-point from the nodal heads, and ``tris`` (n_tri, 3) the
    sub-points of each sub-triangle. tri3 is its own single sub-triangle; tri6 is cut
    at its nodes into three corner triangles and a center one; quadrilaterals are
    cut on their nodes (quad9) or on the nodes and a center point (quad4: the mean
    of the corners; quad8: the mean of the midsides)."""
    if et == 3:
        return np.eye(3), np.array([[0, 1, 2]])
    if et == 6:
        # nodes 0-2 corners, 3 = (0,1), 4 = (1,2), 5 = (2,0) midsides
        return np.eye(6), np.array([[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]])
    if et == 4:
        P = np.vstack([np.eye(4), np.full((1, 4), 0.25)])
        return P, np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    if et == 8:
        c = np.zeros(8)
        c[4:] = 0.25
        P = np.vstack([np.eye(8), c[None, :]])
    else:
        P = np.eye(9)
    # corners 0-3, midsides 4 (0-1), 5 (1-2), 6 (2-3), 7 (3-0), center 8
    tris = np.array([[0, 4, 8], [4, 1, 8], [1, 5, 8], [5, 2, 8],
                     [2, 6, 8], [6, 3, 8], [3, 7, 8], [7, 0, 8]])
    return P, tris


# Order-averaged (Shapley) split of a sub-triangle's water change among its three
# vertices: rows are the 8 states with each vertex at its old (bit 0) or new (bit 1)
# head, columns the vertices; share = W(states) @ _SHAPLEY3.
_SHAPLEY3 = np.zeros((8, 3))
for _v in range(3):
    for _s in range(8):
        if _s & (1 << _v):
            continue
        _k = bin(_s).count('1')
        _wt = {0: 1.0 / 3.0, 1: 1.0 / 6.0, 2: 1.0 / 3.0}[_k]
        _SHAPLEY3[_s | (1 << _v), _v] += _wt
        _SHAPLEY3[_s, _v] -= _wt
_SUBSET_MASK = np.array([[(s >> v) & 1 for v in range(3)] for s in range(8)], dtype=bool)


def _ramp_integral(U, area):
    """Exact integral of max(u, 0) over triangles on which u is linear.

    ``U`` (..., 3) holds the vertex values and ``area`` (...) the areas. Returns
    ``(I, G)``: the integral and its derivatives with respect to the three vertex
    values (continuous in U, so the integral is continuously differentiable)."""
    order = np.argsort(U, axis=-1)
    Us = np.take_along_axis(U, order, axis=-1)
    a, b, c = Us[..., 0], Us[..., 1], Us[..., 2]
    full = a >= 0.0
    I = np.where(full, area * (a + b + c) / 3.0, 0.0)
    Gs = np.where(full[..., None], (area / 3.0)[..., None], 0.0)
    # one vertex above zero (b <= 0 < c): a small triangle at c
    up = (c > 0.0) & (b <= 0.0)
    ca = np.where(up, c - a, 1.0)
    cb = np.where(up, c - b, 1.0)
    Ga = area * c ** 3 / (3.0 * ca * ca * cb)
    Gb = area * c ** 3 / (3.0 * ca * cb * cb)
    I = np.where(up, area * c ** 3 / (3.0 * ca * cb), I)
    Gs = np.where(up[..., None],
                  np.stack([Ga, Gb, area * c * c / (ca * cb) - Ga - Gb], axis=-1), Gs)
    # one vertex below zero (a < 0 < b): the whole triangle less a small one at a
    dn = (a < 0.0) & (b > 0.0)
    ba = np.where(dn, b - a, 1.0)
    ca2 = np.where(dn, c - a, 1.0)
    Gb2 = area / 3.0 + area * a ** 3 / (3.0 * ba * ba * ca2)
    Gc2 = area / 3.0 + area * a ** 3 / (3.0 * ba * ca2 * ca2)
    I = np.where(dn, area * (a + b + c) / 3.0 - area * a ** 3 / (3.0 * ba * ca2), I)
    Gs = np.where(dn[..., None],
                  np.stack([area * (1.0 - a * a / (ba * ca2)) - Gb2 - Gc2, Gb2, Gc2],
                           axis=-1), Gs)
    G = np.empty_like(Gs)
    np.put_along_axis(G, order, Gs, axis=-1)
    return I, G


def _lf_tri_water(U, area, Ss, h0, s_band, s_dry, grad=False):
    """Exact water of the linear-front / Gardner law over triangles on which psi is
    linear. Phi is piecewise linear in psi, with kinks at psi = 0 and psi = h0,

        Phi = Ss max(psi, 0) + s_band (psi - max(psi, 0))
              + (s_dry - s_band) [(psi - h0) - max(psi - h0, 0)],

    so its integral is three ramp integrals; ``Ss`` = 0 gives the drainage water
    Phi(min(psi, 0)) alone. ``U`` (..., 3) vertex heads, ``area`` and the law constants
    broadcast against U[..., 0]. Returns the water and, with ``grad``, its
    derivatives with respect to the vertex heads (..., 3)."""
    S_int = area * U.mean(axis=-1)
    I0, G0 = _ramp_integral(U, area)
    Ih, Gh = _ramp_integral(U - h0[..., None], area)
    W = (Ss * I0 + s_band * (S_int - I0)
         + (s_dry - s_band) * (S_int - h0 * area - Ih))
    if not grad:
        return W, None
    a3 = (area / 3.0)[..., None]
    g = (Ss[..., None] * G0 + s_band[..., None] * (a3 - G0)
         + (s_dry - s_band)[..., None] * (a3 - Gh))
    return W, g


def _build_storage_assembly(nodes, masm, Ss, Sy, h0, vg_a, vg_n, model):
    """Precompute the stored-water evaluation of every element: one storage group
    per element type and law.

    The water an element holds is the integral of Phi(psi) over its storage
    sub-triangles at their true areas: exact for the linear-front / Gardner law on
    triangles (the group carries its sub-triangles and their areas), and otherwise a
    sum over storage points (the group carries their non-negative interpolation
    weights ``N`` (n_q, nn) and their physical weights ``w`` (n_e, n_q)). On tri3 the
    true-area thirds are the HRZ weights, so a saturated tri3 stores exactly its
    lumped mass. On tri6 and quadrilaterals (``split``) they are not, and there the
    water splits in two: the elastic water of the saturated soil, Ss max(psi, 0), is
    lumped at the nodes with their HRZ weights ``m`` -- the lumped mass of a
    saturated element, exactly -- and only the drainage water, Phi(min(psi, 0)), is
    integrated. ``omega`` (n_e, nn) is the weight with which the integrated part of a
    field lying on one linear branch falls to each node."""
    groups = []
    is_vg_all = (np.zeros(len(Ss), dtype=bool) if model is None
                 else np.asarray(model) == KR_VG)
    for g in masm['groups']:
        et, gidx, gconn = g['et'], g['idx'], g['conn']
        P, tris = _storage_cells(et)
        xy = np.einsum('pj,ejd->epd', P, nodes[gconn])
        v = xy[:, tris]
        area_all = 0.5 * np.abs(
            (v[..., 1, 0] - v[..., 0, 0]) * (v[..., 2, 1] - v[..., 0, 1])
            - (v[..., 2, 0] - v[..., 0, 0]) * (v[..., 1, 1] - v[..., 0, 1]))
        omega_all = np.repeat(area_all / 3.0, 3, axis=1) @ P[tris.ravel()]  # (n_e, nn)
        vg_rows = is_vg_all[gidx]
        for vg in (False, True):
            rows = np.where(vg_rows == vg)[0]
            if len(rows) == 0:
                continue
            idx = gidx[rows]
            law = _storage_law(Ss[idx], Sy[idx], h0[idx],
                               None if (vg_a is None or not vg) else vg_a[idx],
                               None if (vg_n is None or not vg) else vg_n[idx])
            sg = {'conn': gconn[rows], 'm': g['me_lump'][rows], 'omega': omega_all[rows],
                  'A_w': area_all[rows].sum(axis=1), 'is_vg': vg, 'split': et != 3,
                  'lf': {k: np.broadcast_to(law[k], (len(idx),)).astype(float)[:, None]
                         for k in ('Ss', 'h0', 'floor_lf', 's_band')},
                  'vg': ({k: (np.broadcast_to(val, (len(idx),))[:, None]
                              if np.ndim(val) else val) for k, val in law.items()}
                         if vg else None)}
            if not vg and et in (3, 6):
                sg['exact'] = True
                sg['tris'] = tris
                sg['area'] = area_all[rows]
            else:
                Lq, Wq = _subtri_rule(STORAGE_SUBDIV_VG if vg else STORAGE_SUBDIV_LF)
                sg['exact'] = False
                sg['N'] = np.concatenate([Lq @ P[t] for t in tris])          # (n_q, nn)
                sg['w'] = (np.repeat(area_all[rows], len(Wq), axis=1)
                           * np.tile(Wq, len(tris)))                         # (n_e, n_q)
            groups.append(sg)
    return {'groups': groups, 'n_nodes': masm['n_nodes']}


def _point_drainage(sg, psi_q, rows=None, want_cap=True):
    """The integrated part of the water at the storage points of a point group,
    and (``want_cap``) its slope: Phi(psi) and the law's capacity, or, on a ``split``
    group, the drainage water Phi(min(psi, 0)) and a slope that is zero above
    psi = 0. ``rows`` restricts the evaluation to those elements."""
    pu = np.minimum(psi_q, 0.0) if sg['split'] else psi_q
    if sg['is_vg']:
        vg = sg['vg'] if rows is None else {k: (v[rows] if np.ndim(v) else v)
                                            for k, v in sg['vg'].items()}
        phi, cap = _vg_point(pu, vg, want_cap)
    else:
        lf = sg['lf'] if rows is None else {k: v[rows] for k, v in sg['lf'].items()}
        phi = _phi_lf(pu, lf)
        cap = (np.where(pu >= 0.0, lf['Ss'],
                        np.where(pu > lf['h0'], lf['s_band'], lf['floor_lf']))
               if want_cap else None)
    if want_cap and sg['split']:
        cap = np.where(psi_q >= 0.0, 0.0, cap)
    return phi, cap


def _branch(sg, psi_n):
    """The linear branch of Phi every node of an element lies on -- 2 saturated,
    1 inside the linear-front band, 0 below it -- or -1 where the element's nodes
    span more than one (or the law has no linear branch there: an unsaturated van
    Genuchten element). The storage points interpolate the nodes with non-negative
    weights, so they lie on the same branch."""
    lo, hi = psi_n.min(axis=1), psi_n.max(axis=1)
    if sg['is_vg']:
        return np.where(lo >= 0.0, 2, -1)
    h0 = sg['lf']['h0'][:, 0]
    return np.where(lo >= 0.0, 2, np.where((hi <= 0.0) & (lo >= h0), 1,
                                           np.where(hi <= h0, 0, -1)))


def _storage_state(sasm, p_nodes):
    """Per group: the nodal pressure heads and the branch of every element, plus,
    for a point group, the storage-point heads and drainage water of its elements on
    no single branch. The step's starting state, evaluated once per attempt."""
    out = []
    for sg in sasm['groups']:
        psi_n = p_nodes[sg['conn']]
        br = _branch(sg, psi_n)
        pos = psi_q = phi = None
        if not sg['exact']:
            ix = np.where(br < 0)[0]
            pos = np.full(len(br), -1)
            pos[ix] = np.arange(len(ix))
            psi_q = psi_n[ix] @ sg['N'].T
            phi = (_point_drainage(sg, psi_q, ix, want_cap=False)[0] if len(ix)
                   else psi_q)
        out.append((psi_n, br, pos, psi_q, phi))
    return out


def _drain_slope(sg, br):
    """Slope of the integrated water on a branch: Ss saturated (0 on a ``split``
    group, whose elastic water is lumped instead), s_band in the band, the floor
    below it."""
    lf = sg['lf']
    sat = 0.0 if sg['split'] else lf['Ss'][:, 0]
    return np.where(br == 2, sat, np.where(br == 1, lf['s_band'][:, 0], lf['floor_lf'][:, 0]))


def _drain_icpt(sg, br):
    """Intercept of the drainage water on a branch: only the drained branch of the
    linear front has one, (s_band - floor) h0."""
    lf = sg['lf']
    return np.where(br == 0, (lf['s_band'][:, 0] - lf['floor_lf'][:, 0]) * lf['h0'][:, 0],
                    0.0)


def _lf_rows(sg, rows, extra_dims):
    """(Ss, h0, s_band, floor) of the given rows, shaped to broadcast; Ss is 0 on a
    ``split`` group."""
    lf = sg['lf']
    shape = (len(rows),) + (1,) * extra_dims
    Ss = lf['Ss'][rows, 0].reshape(shape) * (0.0 if sg['split'] else 1.0)
    return (Ss,) + tuple(lf[k][rows, 0].reshape(shape) for k in ('h0', 's_band', 'floor_lf'))


def _stored_water(sasm, p_nodes):
    """Water held in each element (list per group), measured from saturation."""
    out = []
    for sg, (psi_n, br, _, _, _) in zip(sasm['groups'], _storage_state(sasm, p_nodes)):
        elastic = (sg['lf']['Ss'][:, 0] * (sg['m'] * np.maximum(psi_n, 0.0)).sum(axis=1)
                   if sg['split'] else np.zeros(len(br)))
        W = (elastic + _drain_slope(sg, br) * (sg['omega'] * psi_n).sum(axis=1)
             + _drain_icpt(sg, br) * sg['A_w'])
        ix = np.where(br < 0)[0]
        if len(ix):
            W[ix] = elastic[ix]
            if sg['exact']:
                Ss, h0, sb, sd = _lf_rows(sg, ix, 1)
                W[ix] += _lf_tri_water(psi_n[ix][:, sg['tris']], sg['area'][ix],
                                       Ss, h0, sb, sd)[0].sum(axis=1)
            else:
                phi, _ = _point_drainage(sg, psi_n[ix] @ sg['N'].T, ix, want_cap=False)
                W[ix] += (sg['w'][ix] * phi).sum(axis=1)
        out.append(W)
    return out


def _storage_residual(sasm, p_nodes, p_old, old_state):
    """Storage term of one time step, from pressure head ``p_old`` to ``p_nodes``.

    Each element's water change is split exactly into one term per node, each term
    driven by that node's own head change, so a node whose head does not move
    stores nothing and the term stays lumped (see ``_build_storage_assembly`` for
    what the water of an element is):

    * the lumped elastic water of a tri6 or quadrilateral element: node j takes up
      Ss m_j times the change of max(psi_j, 0), m_j its HRZ weight;
    * the integrated water of an element whose nodes stay on one linear branch of
      Phi through the step: node j takes up omega_j S dpsi_j, S the branch's slope
      and omega_j the node's share of the element's area (its HRZ weight on tri3);
    * the linear-front / Gardner law on a triangle otherwise: each sub-triangle's
      water is its exact integral, and its change is split among its three vertices
      by the average, over the six orders in which the vertices could move from
      their old head to their new one, of each vertex's own contribution (8 exact
      evaluations; the shares sum to the change and each has the sign of its own
      head change);
    * a sum over storage points otherwise, W = sum_q w_q Phi(psi_q),
      psi_q = sum_j N_qj psi_j with N_qj >= 0: with each point's own chord C*_q the
      change splits as sum_j (sum_q w_q N_qj C*_q) dpsi_j.

    A saturated element therefore stores m_j Ss dpsi_j at each node: the lumped
    mass, exactly, on every element type.

    Returns ``(R, D, W)``: the nodal storage (volume per unit width), the diagonal
    the iteration linearizes it with (at each node, the larger of its chord
    capacity over the step and its tangent capacity at p), and the per-group
    element water at p."""
    n = sasm['n_nodes']
    R = np.zeros(n)
    D = np.zeros(n)
    W_new = []
    for sg, (psi_o, br_o, pos_o, psi_qo_c, phi_o_c) in zip(sasm['groups'], old_state):
        conn = sg['conn']
        psi_n = p_nodes[conn]
        dpsi = psi_n - psi_o
        # elastic water of the saturated soil, lumped at the nodes (split groups)
        Ssm = sg['lf']['Ss'] * sg['m'] * (1.0 if sg['split'] else 0.0)
        pos_n, pos_old = np.maximum(psi_n, 0.0), np.maximum(psi_o, 0.0)
        Re = Ssm * (pos_n - pos_old)
        small_n = np.abs(dpsi) <= 1e-12 * (1.0 + np.abs(psi_n))
        frac = np.where(small_n, (psi_n > 0.0).astype(float),
                        (pos_n - pos_old) / np.where(small_n, 1.0, dpsi))
        De = Ssm * np.maximum(frac, (psi_n > 0.0).astype(float))
        W = (Ssm * pos_n).sum(axis=1)
        # drainage water
        br = _branch(sg, psi_n)
        lin = (br >= 0) & (br == br_o)
        slope = _drain_slope(sg, br)
        Dl = slope[:, None] * sg['omega']
        Re = Re + np.where(lin[:, None], Dl * dpsi, 0.0)
        De = De + np.where(lin[:, None], Dl, 0.0)
        W = W + np.where(lin, slope * (sg['omega'] * psi_n).sum(axis=1)
                         + _drain_icpt(sg, br) * sg['A_w'], 0.0)
        ix = np.where(~lin)[0]
        if len(ix) and sg['exact']:
            tris = sg['tris']
            uo = psi_o[ix][:, tris]                               # (n_k, n_tri, 3)
            un = psi_n[ix][:, tris]
            U = np.where(_SUBSET_MASK[None, None], un[:, :, None, :], uo[:, :, None, :])
            Ss, h0, sb, sd = _lf_rows(sg, ix, 2)
            area = sg['area'][ix][:, :, None]
            Ws, _ = _lf_tri_water(U, area, Ss, h0, sb, sd)        # (n_k, n_tri, 8)
            share = Ws @ _SHAPLEY3                                 # (n_k, n_tri, 3)
            _, tan = _lf_tri_water(un, area[..., 0], Ss[..., 0], h0[..., 0], sb[..., 0],
                                   sd[..., 0], grad=True)
            du = un - uo
            small = np.abs(du) <= 1e-12 * (1.0 + np.abs(un))
            chord = np.where(small, tan, share / np.where(small, 1.0, du))
            dv = np.maximum(chord, tan)
            Rk = np.zeros((len(ix), conn.shape[1]))
            Dk = np.zeros((len(ix), conn.shape[1]))
            for t in range(len(tris)):
                Rk[:, tris[t]] += share[:, t, :]
                Dk[:, tris[t]] += dv[:, t, :]
            Re[ix] += Rk
            De[ix] += Dk
            W[ix] += Ws[..., 7].sum(axis=1)
        elif len(ix):
            N = sg['N']
            psi_q = psi_n[ix] @ N.T
            phi, cap = _point_drainage(sg, psi_q, ix)
            # starting state: cached where the element was on no single branch, the
            # branch's own linear drainage water where it was
            psi_qo = psi_o[ix] @ N.T
            phi_o = np.empty_like(psi_qo)
            c = pos_o[ix] >= 0
            if np.any(c):
                phi_o[c] = phi_o_c[pos_o[ix][c]]
                psi_qo[c] = psi_qo_c[pos_o[ix][c]]
            if np.any(~c):
                r = ix[~c]
                phi_o[~c] = (_drain_slope(sg, br_o)[r][:, None] * psi_qo[~c]
                             + _drain_icpt(sg, br_o)[r][:, None])
            dq = psi_q - psi_qo
            small = np.abs(dq) <= 1e-12 * (1.0 + np.abs(psi_q))
            if np.any(small):
                cap_mid = _point_drainage(sg, 0.5 * (psi_q + psi_qo), ix)[1]
            else:
                cap_mid = cap
            chord = np.where(small, cap_mid, (phi - phi_o) / np.where(small, 1.0, dq))
            w = sg['w'][ix]
            Re[ix] += ((w * chord) @ N) * dpsi[ix]
            De[ix] += (w * np.maximum(chord, cap)) @ N
            W[ix] += (w * phi).sum(axis=1)
        R += np.bincount(conn.ravel(), weights=Re.ravel(), minlength=n)
        D += np.bincount(conn.ravel(), weights=De.ravel(), minlength=n)
        W_new.append(W)
    return R, D, W_new



def create_flow_potential_bc_from_elements(nodes, elements, element_types, head,
                                           k1_vals, k2_vals, angles,
                                           kr0=None, h0=None, total_flow=None,
                                           bc_type=None, exit_face_active=None,
                                           vg_a=None, vg_n=None, vg_l=None,
                                           model=None):
    """
    Generates Dirichlet BCs for flow potential φ by integrating the Darcy
    velocity flux along each boundary edge from the owning element's shape
    functions. This avoids the reaction-force artifacts that plague the
    q-based approach at sharp transitions (phreatic exit point).

    For each boundary edge, the stream function change is:
        Δψ = ∫(-vy·dx + vx·dy) along the edge
    where v = -kr·K·grad(h) is evaluated from the element's shape functions.

    Parameters:
        nodes, elements, element_types: mesh data
        head: (n_nodes,) nodal head solution
        k1_vals, k2_vals, angles: per-element conductivity properties
        kr0, h0: per-element unsaturated parameters (optional)

    For quadratic boundary edges, the active seepage or fixed-head interval can
    occupy only half the side when the boundary-condition transition falls at a
    midside node. In that case the accumulated stream-function change must be
    restricted to the active half-edge rather than spread over the full side.

    Returns:
        List of (node_id, phi_value) tuples
    """
    from collections import defaultdict

    if element_types is None:
        element_types = np.full(len(elements), 3)

    y = nodes[:, 1]
    p_nodes = head - y

    # Step 1: Build corner-only boundary edges, midside map, and edge→element map
    boundary_edges, midside_map = _boundary_edge_map(elements, element_types)

    # Step 2: Walk boundary (corner nodes only)
    neighbor_map = defaultdict(list)
    for a, b in boundary_edges:
        neighbor_map[a].append(b)
        neighbor_map[b].append(a)
    start_node = list(boundary_edges.keys())[0][0]
    ordered_nodes = [start_node]
    visited = {start_node}
    current = start_node
    while True:
        nbrs = [nn for nn in neighbor_map[current] if nn not in visited]
        if not nbrs:
            break
        nxt = nbrs[0]
        ordered_nodes.append(nxt)
        visited.add(nxt)
        current = nxt
    n = len(ordered_nodes)

    def get_quadratic_edge_interval(c1, mid, c2):
        """Return the active interval on a quadratic boundary edge.

        The interval is inferred from boundary-condition flags on the
        corner-midside-corner triplet. For the stream-function BCs, only a
        fully active quadratic side is treated as a flux-carrying segment.
        Mixed edge patterns are left constant in phi; allowing half-edge
        stream-function jumps on tri6 exit-face transition edges distorts the
        flownet even though the head field itself remains acceptable.
        """
        if bc_type is None or mid is None:
            return None

        edge_nodes = np.array([c1, mid, c2], dtype=int)
        fixed_mask = bc_type[edge_nodes] == 1
        exit_mask = np.zeros(3, dtype=bool)
        if exit_face_active is not None:
            exit_mask = (bc_type[edge_nodes] == 2) & exit_face_active[edge_nodes]

        for mask in (fixed_mask, exit_mask):
            if np.array_equal(mask, [True, True, True]):
                return (0.0, 1.0)
        return None

    # Step 3: For each boundary segment, compute Δψ from element-level velocity
    segment_flux = {}
    midside_flux_from_start = {}
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        elem_idx = boundary_edges.get(ekey)
        mid = midside_map.get(ekey)
        if elem_idx is None:
            segment_flux[c1] = 0.0
            continue

        en = elements[elem_idx]
        et = element_types[elem_idx]
        k1 = k1_vals[elem_idx] if hasattr(k1_vals, '__len__') else k1_vals
        k2 = k2_vals[elem_idx] if hasattr(k2_vals, '__len__') else k2_vals
        theta = angles[elem_idx] if hasattr(angles, '__len__') else angles
        theta_rad = np.radians(theta)
        cs, sn = np.cos(theta_rad), np.sin(theta_rad)
        R = np.array([[cs, sn], [-sn, cs]])
        Kmat = R.T @ np.diag([k1, k2]) @ R

        # kr factor
        kr_elem = 1.0
        if kr0 is not None and h0 is not None:
            if et == 6:
                p_elem = compute_tri6_centroid_pressure(p_nodes, en)
            elif et in (3, 4):
                p_elem = np.mean(p_nodes[en[:3 if et == 3 else 4]])
            else:
                p_elem = np.mean(p_nodes[en[:4]])
            kr_e0 = kr0[elem_idx] if hasattr(kr0, '__len__') else kr0
            h0_e0 = h0[elem_idx] if hasattr(h0, '__len__') else h0
            vga_e0 = vg_a[elem_idx] if hasattr(vg_a, '__len__') else vg_a
            vgn_e0 = vg_n[elem_idx] if hasattr(vg_n, '__len__') else vg_n
            vgl_e0 = vg_l[elem_idx] if hasattr(vg_l, '__len__') else vg_l
            mdl_e0 = model[elem_idx] if hasattr(model, '__len__') else (model or KR_LF)
            kr_elem = kr_relative(p_elem, kr_e0, h0_e0, vga_e0, vgn_e0, vgl_e0, mdl_e0)

        # Edge tangent
        dx = nodes[c2, 0] - nodes[c1, 0]
        dy = nodes[c2, 1] - nodes[c1, 1]

        if et == 3:
            # Constant gradient — single evaluation
            i0, j0, k0 = en[:3]
            xi, yi = nodes[i0]; xj, yj = nodes[j0]; xk, yk = nodes[k0]
            area = 0.5 * abs((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi))
            if area < 1e-30:
                segment_flux[c1] = 0.0
                continue
            beta = np.array([yj-yk, yk-yi, yi-yj])
            gamma = np.array([xk-xj, xi-xk, xj-xi])
            grad = np.array([beta, gamma]) / (2*area)
            grad_h = grad @ head[en[:3]]
            v = -kr_elem * Kmat @ grad_h
            segment_flux[c1] = -v[1]*dx + v[0]*dy

        elif et == 6:
            # Linear gradient — 2-point Gauss on edge
            nodes_elem = nodes[en[:6]]
            h_elem = head[en[:6]]
            x0, y0 = nodes_elem[0]; x1, y1 = nodes_elem[1]; x2, y2 = nodes_elem[2]
            J = np.array([[x0-x2, x1-x2], [y0-y2, y1-y2]])
            detJ = np.linalg.det(J)
            if abs(detJ) < 1e-10:
                segment_flux[c1] = 0.0
                continue
            total_area = 0.5 * abs(detJ)
            dL1_dx = (y1-y2)/(2*total_area); dL1_dy = (x2-x1)/(2*total_area)
            dL2_dx = (y2-y0)/(2*total_area); dL2_dy = (x0-x2)/(2*total_area)
            dL3_dx = (y0-y1)/(2*total_area); dL3_dy = (x1-x0)/(2*total_area)

            # Determine which edge of the element this boundary edge is on
            c1_local = list(en[:3]).index(c1) if c1 in en[:3] else -1
            c2_local = list(en[:3]).index(c2) if c2 in en[:3] else -1
            if c1_local < 0 or c2_local < 0:
                segment_flux[c1] = 0.0
                continue

            # Parameterize edge: t=0 at c1, t=1 at c2
            # Map to area coordinates based on which edge
            def area_coords(t):
                L = [0.0, 0.0, 0.0]
                L[c1_local] = 1.0 - t
                L[c2_local] = t
                # Third coordinate = 0 (on the edge)
                return L[0], L[1], L[2]

            # 2-point Gauss quadrature on an edge interval [t_start, t_end]
            # with per-Gauss-point kr. For quadratic edges, this lets us set
            # midside phi from the integrated half-edge flux instead of
            # forcing it to the average of the corner values.
            p_elem_nodes = p_nodes[en[:6]]
            kr_e0 = kr0[elem_idx] if hasattr(kr0, '__len__') else kr0
            h0_e0 = h0[elem_idx] if hasattr(h0, '__len__') else h0
            vga_e0 = vg_a[elem_idx] if hasattr(vg_a, '__len__') else vg_a
            vgn_e0 = vg_n[elem_idx] if hasattr(vg_n, '__len__') else vg_n
            vgl_e0 = vg_l[elem_idx] if hasattr(vg_l, '__len__') else vg_l
            mdl_e0 = model[elem_idx] if hasattr(model, '__len__') else (model or KR_LF)
            gp = [0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)]

            def integrate_edge_interval(t_start, t_end):
                total = 0.0
                interval = t_end - t_start
                for xi in gp:
                    t = t_start + interval * xi
                    L1, L2, L3 = area_coords(t)
                    N = np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                                  4*L1*L2, 4*L2*L3, 4*L3*L1])
                    p_gp = N @ p_elem_nodes
                    kr_gp = kr_relative(p_gp, kr_e0, h0_e0, vga_e0, vgn_e0, vgl_e0, mdl_e0) if kr0 is not None else 1.0
                    dN_dL1 = np.array([4*L1-1, 0, 0, 4*L2, 0, 4*L3])
                    dN_dL2 = np.array([0, 4*L2-1, 0, 4*L1, 4*L3, 0])
                    dN_dL3 = np.array([0, 0, 4*L3-1, 0, 4*L2, 4*L1])
                    gradN = np.zeros((2, 6))
                    for ii in range(6):
                        gradN[0, ii] = dN_dL1[ii]*dL1_dx + dN_dL2[ii]*dL2_dx + dN_dL3[ii]*dL3_dx
                        gradN[1, ii] = dN_dL1[ii]*dL1_dy + dN_dL2[ii]*dL2_dy + dN_dL3[ii]*dL3_dy
                    grad_h = gradN @ h_elem
                    v = -kr_gp * Kmat @ grad_h
                    total += -v[1]*dx + v[0]*dy
                return total * 0.5 * interval

            edge_interval = get_quadratic_edge_interval(c1, mid, c2)
            segment_flux[c1] = integrate_edge_interval(0.0, 1.0)
            if edge_interval == (0.0, 1.0):
                midside_flux_from_start[c1] = integrate_edge_interval(0.0, 0.5)
            # Transition edges: segment_flux is kept (correct corner phi walk)
            # but midside_flux_from_start is omitted so the midside node gets
            # linear interpolation between corners, avoiding half-edge phi
            # jumps that distort the flownet on mixed BC edges.

        elif et == 4:
            # Bilinear quad — evaluate at edge midpoint
            i0, j0, k0, l0 = en[:4]
            nodes_elem = nodes[en[:4]]
            h_elem = head[en[:4]]
            mid_x = 0.5*(nodes[c1,0]+nodes[c2,0])
            mid_y = 0.5*(nodes[c1,1]+nodes[c2,1])
            # Use centroid gradient as approximation
            xc = np.mean(nodes_elem[:,0]); yc = np.mean(nodes_elem[:,1])
            # Simple: use the tri3-style constant gradient from two triangles
            xi,yi = nodes_elem[0]; xj,yj = nodes_elem[1]; xk,yk = nodes_elem[2]; xl,yl = nodes_elem[3]
            area1 = 0.5*abs((xj-xi)*(yk-yi)-(xk-xi)*(yj-yi))
            area2 = 0.5*abs((xk-xi)*(yl-yi)-(xl-xi)*(yk-yi))
            if area1+area2 < 1e-30:
                segment_flux[c1] = 0.0
                continue
            beta1 = np.array([yj-yk, yk-yi, yi-yj])
            gamma1 = np.array([xk-xj, xi-xk, xj-xi])
            grad1 = np.array([beta1, gamma1])/(2*area1)
            grad_h1 = grad1 @ h_elem[:3]
            beta2 = np.array([yk-yl, yl-yi, yi-yk])
            gamma2 = np.array([xl-xk, xi-xl, xk-xi])
            grad2 = np.array([beta2, gamma2])/(2*area2)
            grad_h2 = grad2 @ h_elem[[0,2,3]]
            grad_h = (grad_h1*area1 + grad_h2*area2)/(area1+area2)
            v = -kr_elem * Kmat @ grad_h
            segment_flux[c1] = -v[1]*dx + v[0]*dy

        elif et in (8, 9):
            # quad8/quad9: use corner nodes same as quad4
            nodes_elem = nodes[en[:4]]
            h_elem = head[en[:4]]
            xi,yi = nodes_elem[0]; xj,yj = nodes_elem[1]; xk,yk = nodes_elem[2]; xl,yl = nodes_elem[3]
            area1 = 0.5*abs((xj-xi)*(yk-yi)-(xk-xi)*(yj-yi))
            area2 = 0.5*abs((xk-xi)*(yl-yi)-(xl-xi)*(yk-yi))
            if area1+area2 < 1e-30:
                segment_flux[c1] = 0.0
                continue
            beta1 = np.array([yj-yk, yk-yi, yi-yj])
            gamma1 = np.array([xk-xj, xi-xk, xj-xi])
            grad1 = np.array([beta1, gamma1])/(2*area1)
            grad_h1 = grad1 @ h_elem[:3]
            beta2 = np.array([yk-yl, yl-yi, yi-yk])
            gamma2 = np.array([xl-xk, xi-xl, xk-xi])
            grad2 = np.array([beta2, gamma2])/(2*area2)
            grad_h2 = grad2 @ h_elem[[0,2,3]]
            grad_h = (grad_h1*area1 + grad_h2*area2)/(area1+area2)
            v = -kr_elem * Kmat @ grad_h
            segment_flux[c1] = -v[1]*dx + v[0]*dy

    # Step 4: Determine walk orientation from signed area
    signed_area = 0.0
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        signed_area += nodes[c1, 0] * nodes[c2, 1] - nodes[c2, 0] * nodes[c1, 1]
    if signed_area < 0:
        segment_flux = {k: -v for k, v in segment_flux.items()}
        midside_flux_from_start = {k: -v for k, v in midside_flux_from_start.items()}

    # Step 4b: Normalize inflow and outflow separately so the accumulated
    # boundary phi remains single-valued around the closed loop.
    # Using one global scale factor matches the inflow sum but can leave a
    # non-zero closure error when the raw boundary flux extraction slightly
    # over/under-estimates the outflow, which is especially noticeable for
    # quadratic triangles near the seepage transition.
    if total_flow is not None:
        pos_sum = sum(f for f in segment_flux.values() if f > 0)
        neg_sum = -sum(f for f in segment_flux.values() if f < 0)
        pos_scale = total_flow / pos_sum if pos_sum > 1e-30 else 1.0
        neg_scale = total_flow / neg_sum if neg_sum > 1e-30 else 1.0

        scaled_segment_flux = {}
        scaled_midside_flux = {}
        for k, v in segment_flux.items():
            if v > 0:
                scale = pos_scale
            elif v < 0:
                scale = neg_scale
            else:
                scale = 1.0
            scaled_segment_flux[k] = v * scale
            if k in midside_flux_from_start:
                scaled_midside_flux[k] = midside_flux_from_start[k] * scale

        segment_flux = scaled_segment_flux
        midside_flux_from_start = scaled_midside_flux

    # Step 5: Find starting point (max inflow)
    start_idx = max(range(n), key=lambda i: segment_flux.get(ordered_nodes[i], 0))

    # Step 6: Accumulate phi at corners
    total_q = total_flow if total_flow is not None else sum(f for f in segment_flux.values() if f > 0)
    phi_val = total_q
    phi = {}
    for i in range(n):
        idx = (start_idx + i) % n
        node = ordered_nodes[idx]
        phi[node] = phi_val
        phi_val -= segment_flux[node]

    # Assign quadratic-edge midside values from the integrated half-edge flux.
    # This keeps the boundary phi field compatible with the quadratic head field.
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        mid = midside_map.get(ekey)
        if mid is not None:
            if c1 in midside_flux_from_start:
                phi[mid] = phi[c1] - midside_flux_from_start[c1]
            else:
                phi[mid] = 0.5 * (phi[c1] + phi[c2])

    # With a known total_flow, separate inflow/outflow normalization keeps the
    # accumulated boundary phi single-valued around the closed loop. Any
    # remaining local inconsistencies are resolved by the stream-function PDE.

    return list(phi.items())


def solve_flow_function_confined(nodes, elements, k1_vals, k2_vals, angles, dirichlet_nodes, element_types=None):
    """
    Solves the stream function (flow function) Phi on the same mesh,
    assigning Dirichlet values along no-flow boundaries.

    For anisotropic permeability, the stream function equation uses K/det(K)
    (not K^(-1)) in the stiffness matrix assembly. This is because the stream
    function PDE has swapped diagonal coefficients compared to the head equation.

    Supports both triangular and quadrilateral elements.
    
    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 3 or 4) triangle or quad node indices
        k1_vals : (n_elements,) or scalar, major axis conductivity
        k2_vals : (n_elements,) or scalar, minor axis conductivity
        angles : (n_elements,) or scalar, angle in degrees (from x-axis)
        dirichlet_nodes : list of (node_id, phi_value)
        element_types : (n_elements,) array indicating 3 for triangles, 4 for quads
    Returns:
        phi : (n_nodes,) stream function (flow function) values
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    # Stream-function coefficient matrix is K/det(K) (flow=True)
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals,
                          angles, flow=True)
    data = _assembly_data(asm)

    dir_mask = np.zeros(n_nodes, dtype=bool)
    dir_values = np.zeros(n_nodes)
    for node, phi_value in dirichlet_nodes:
        dir_mask[node] = True
        dir_values[node] = phi_value
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

    phi = spsolve(A, b)
    return phi

def solve_flow_function_unsaturated(nodes, elements, head, k1_vals, k2_vals, angles, kr0, h0, dirichlet_nodes, element_types=None, vg_a=None, vg_n=None, vg_l=None, model=None):
    """
    Solves the stream function (flow function) Phi for unsaturated flow.

    For anisotropic permeability, the stream function equation uses K/det(K)
    (not K^(-1)) in the stiffness matrix assembly. The relative permeability
    kr is also included in the formulation (linear-front or van Genuchten per
    the per-element ``model``).

    Supports both triangular and quadrilateral elements.
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    y = nodes[:, 1]
    p_nodes = head - y

    kr0 = np.asarray(kr0, dtype=float)
    h0 = np.asarray(h0, dtype=float)
    if kr0.ndim == 0:
        kr0 = np.full(len(elements), float(kr0))
    if h0.ndim == 0:
        h0 = np.full(len(elements), float(h0))
    if vg_a is not None:
        vg_a = np.broadcast_to(np.asarray(vg_a, dtype=float), (len(elements),))
    if vg_n is not None:
        vg_n = np.broadcast_to(np.asarray(vg_n, dtype=float), (len(elements),))
    if vg_l is not None:
        vg_l = np.broadcast_to(np.asarray(vg_l, dtype=float), (len(elements),))
    if model is not None:
        model = np.broadcast_to(np.asarray(model), (len(elements),))

    # Stream-function coefficient matrix is K/det(K); each element is scaled
    # by 1/kr_avg (mode='stream'), matching the per-element *_kr functions.
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals,
                          angles, flow=True)
    data = _assembly_data(asm, p_nodes=p_nodes, kr0=kr0, h0=h0, mode='stream',
                          vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)

    dir_mask = np.zeros(n_nodes, dtype=bool)
    dir_values = np.zeros(n_nodes)
    for node, phi_value in dirichlet_nodes:
        dir_mask[node] = True
        dir_values[node] = phi_value
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

    phi = spsolve(A, b)
    return phi


def compute_velocity(nodes, elements, head, k1_vals, k2_vals, angles, kr0=None, h0=None, element_types=None, vg_a=None, vg_n=None, vg_l=None, model=None):
    """
    Compute nodal velocities by averaging element-wise Darcy velocities.
    If kr0 and h0 are provided, compute kr_elem using kr_frontal; otherwise, kr_elem = 1.0.
    Supports both triangular and quadrilateral elements.
    For quads, velocity is computed at Gauss points and averaged to nodes.
    
    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 3 or 4) triangle or quad node indices
        head : (n_nodes,) nodal head solution
        k1_vals, k2_vals, angles : per-element anisotropic properties (or scalar)
        kr0 : (n_elements,) or scalar, relative permeability parameter (optional)
        h0 : (n_elements,) or scalar, pressure head parameter (optional)
        element_types : (n_elements,) node count per element: 3 (tri3), 4 (quad4),
            6 (tri6), 8 (quad8) or 9 (quad9). Any other value raises.
    
    Returns:
        velocity : (n_nodes, 2) array of nodal velocity vectors [vx, vy]
    """
    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    elements = np.asarray(elements)
    element_types = np.asarray(element_types)
    n_nodes = nodes.shape[0]
    n_el = len(elements)
    velocity = np.zeros((n_nodes, 2))
    count = np.zeros(n_nodes)

    Kmats = _element_kmats(n_el, k1_vals, k2_vals, angles)
    use_kr = kr0 is not None and h0 is not None
    if use_kr:
        kr0 = np.broadcast_to(np.asarray(kr0, dtype=float), (n_el,))
        h0 = np.broadcast_to(np.asarray(h0, dtype=float), (n_el,))
        if vg_a is not None:
            vg_a = np.broadcast_to(np.asarray(vg_a, dtype=float), (n_el,))
        if vg_n is not None:
            vg_n = np.broadcast_to(np.asarray(vg_n, dtype=float), (n_el,))
        if vg_l is not None:
            vg_l = np.broadcast_to(np.asarray(vg_l, dtype=float), (n_el,))
        if model is not None:
            model = np.broadcast_to(np.asarray(model), (n_el,))
    p_all = head - nodes[:, 1]

    for et in np.unique(element_types):
        idx = np.where(element_types == et)[0]
        nn = int(et)
        conn = elements[idx][:, :nn]
        coords = nodes[conn]              # (n_e, nn, 2)
        h_el = head[conn]                 # (n_e, nn)
        K = Kmats[idx]

        if et == 3:
            x, y = coords[:, :, 0], coords[:, :, 1]
            area = 0.5 * ((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0])
                          - (x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
            ok = area > 0
            a_safe = np.where(ok, area, 1.0)
            beta = np.stack([y[:, 1]-y[:, 2], y[:, 2]-y[:, 0], y[:, 0]-y[:, 1]], axis=1)
            gamma = np.stack([x[:, 2]-x[:, 1], x[:, 0]-x[:, 2], x[:, 1]-x[:, 0]], axis=1)
            grad = np.stack([beta, gamma], axis=1) / (2 * a_safe)[:, None, None]
            grad_h = np.einsum('exn,en->ex', grad, h_el)
            if use_kr:
                kr_e = kr_relative_vec(p_all[conn].mean(axis=1), kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(vg_l, idx), _idx_or_none(model, idx))
            else:
                kr_e = np.ones(len(idx))
            v_e = -kr_e[:, None] * np.einsum('exy,ey->ex', K, grad_h)
            v_e[~ok] = 0.0
            np.add.at(velocity, conn.ravel(), np.repeat(v_e, 3, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(ok.astype(float), 3))

        elif et in (4, 8, 9):
            if use_kr:
                # One kr per element, from the mean pressure at its four CORNER
                # nodes. For quad4 that is every node it has. For quad8 and quad9
                # it is an APPROXIMATION and is exact for nothing unsaturated: a
                # serendipity quad8 centroid weights corners -1/4 and midsides
                # +1/2, and a Lagrange quad9 centroid is its centre node alone, so
                # neither equals the corner average, and neither matches what the
                # assembly does either -- _assembly_data averages kr over the same
                # 3x3 Gauss points it integrates on. Recovering kr at those Gauss
                # points, so the scaling here is the scaling the stiffness was
                # built with, is the consistent fix and a known follow-up. No corpus
                # model is on a quad mesh, and the guard runs saturated, where kr is
                # 1 and the choice cannot bite -- but quad meshing is a user option,
                # so the size of the error was measured rather than assumed. On the
                # two seepage sample models, one confined and one with a free
                # surface, quad8 and quad9 put the section flowrate within 0.3% of
                # the tri6 answer on the same geometry (sheetpile 40.100 / 39.886
                # against 40.058; earth dam 38.673 / 38.672 against 38.695), which
                # is inside the 5% tolerance those rows lock at.
                kr_e = kr_relative_vec(p_all[conn[:, :4]].mean(axis=1), kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(vg_l, idx), _idx_or_none(model, idx))
            else:
                kr_e = np.ones(len(idx))
            v_sum = np.zeros((len(idx), 2))
            n_ok = np.zeros(len(idx))
            for xi, eta, _w in _quad_gauss(et):
                dxi, deta = _quad_dshape(et, xi, eta)
                J00 = coords[:, :, 0] @ dxi
                J01 = coords[:, :, 1] @ dxi
                J10 = coords[:, :, 0] @ deta
                J11 = coords[:, :, 1] @ deta
                detJ = J00*J11 - J01*J10
                ok = detJ > 0
                d_safe = np.where(ok, detJ, 1.0)
                dN_dx = (J11[:, None]*dxi[None, :] - J01[:, None]*deta[None, :]) / d_safe[:, None]
                dN_dy = (-J10[:, None]*dxi[None, :] + J00[:, None]*deta[None, :]) / d_safe[:, None]
                grad_h = np.stack([(dN_dx*h_el).sum(axis=1), (dN_dy*h_el).sum(axis=1)], axis=1)
                v_gp = -kr_e[:, None] * np.einsum('exy,ey->ex', K, grad_h)
                v_gp[~ok] = 0.0
                v_sum += v_gp
                n_ok += ok
            np.add.at(velocity, conn.ravel(), np.repeat(v_sum, nn, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(n_ok, nn))

        elif et == 6:
            x, y = coords[:, :, 0], coords[:, :, 1]
            detJ = (x[:, 0]-x[:, 2])*(y[:, 1]-y[:, 2]) - (x[:, 1]-x[:, 2])*(y[:, 0]-y[:, 2])
            ok = np.abs(detJ) > 1e-10
            d_safe = np.where(ok, detJ, 1.0)
            # Jinv of the constant corner-node Jacobian
            Ji00 = (y[:, 1]-y[:, 2]) / d_safe
            Ji01 = -(y[:, 0]-y[:, 2]) / d_safe
            Ji10 = -(x[:, 1]-x[:, 2]) / d_safe
            Ji11 = (x[:, 0]-x[:, 2]) / d_safe
            if use_kr:
                # kr at the element centroid via quadratic shape functions
                Nc = np.array([1/3*(2/3-1)]*3 + [4/9]*3)
                p_c = (p_all[conn] * Nc[None, :]).sum(axis=1)
                kr_e = kr_relative_vec(p_c, kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(vg_l, idx), _idx_or_none(model, idx))
            else:
                kr_e = np.ones(len(idx))
            w_total = np.zeros(len(idx))
            v_sum = np.zeros((len(idx), 2))
            for (L1, L2, L3), w in zip([(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)],
                                       [1/3, 1/3, 1/3]):
                dN_dL1 = np.array([4*L1-1, 0, 0, 4*L2, 0, 4*L3])
                dN_dL2 = np.array([0, 4*L2-1, 0, 4*L1, 4*L3, 0])
                dN_dL3 = np.array([0, 0, 4*L3-1, 0, 4*L2, 4*L1])
                dxi = dN_dL1 - dN_dL3
                deta = dN_dL2 - dN_dL3
                dN_dx = Ji00[:, None]*dxi[None, :] + Ji01[:, None]*deta[None, :]
                dN_dy = Ji10[:, None]*dxi[None, :] + Ji11[:, None]*deta[None, :]
                grad_h = np.stack([(dN_dx*h_el).sum(axis=1), (dN_dy*h_el).sum(axis=1)], axis=1)
                v_gp = -kr_e[:, None] * np.einsum('exy,ey->ex', K, grad_h)
                v_gp[~ok] = 0.0
                v_sum += w * v_gp
                w_total += np.where(ok, w, 0.0)
            np.add.at(velocity, conn.ravel(), np.repeat(v_sum, 6, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(w_total, 6))

        else:
            raise ValueError(
                f"compute_velocity has no branch for element type {int(et)} "
                f"({len(idx)} elements); it recovers velocity on element types "
                f"3 (tri3), 4 (quad4), 6 (tri6), 8 (quad8) and 9 (quad9). "
                f"Falling through would leave those nodes at zero velocity, "
                f"which reads as still water rather than as a missing branch.")

    count[count == 0] = 1  # Avoid division by zero
    velocity /= count[:, None]
    return velocity

def compute_gradient(nodes, elements, head, element_types=None):
    """
    Compute nodal hydraulic gradient by averaging element-wise head gradients.
    The hydraulic gradient i = -grad(h), where grad(h) is the gradient of head.
    Supports both triangular and quadrilateral elements.
    
    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 3 or 4) triangle or quad node indices
        head : (n_nodes,) nodal head solution
        element_types : (n_elements,) node count per element: 3 (tri3), 4 (quad4),
            6 (tri6), 8 (quad8) or 9 (quad9). Any other value raises.
    
    Returns:
        gradient : (n_nodes, 2) array of nodal hydraulic gradient vectors [ix, iy]
    """
    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    elements = np.asarray(elements)
    element_types = np.asarray(element_types)
    n_nodes = nodes.shape[0]
    gradient = np.zeros((n_nodes, 2))
    count = np.zeros(n_nodes)

    for et in np.unique(element_types):
        idx = np.where(element_types == et)[0]
        nn = int(et)
        conn = elements[idx][:, :nn]
        coords = nodes[conn]
        h_el = head[conn]

        if et == 3:
            x, y = coords[:, :, 0], coords[:, :, 1]
            area = 0.5 * ((x[:, 1]-x[:, 0])*(y[:, 2]-y[:, 0])
                          - (x[:, 2]-x[:, 0])*(y[:, 1]-y[:, 0]))
            ok = area > 0
            a_safe = np.where(ok, area, 1.0)
            beta = np.stack([y[:, 1]-y[:, 2], y[:, 2]-y[:, 0], y[:, 0]-y[:, 1]], axis=1)
            gamma = np.stack([x[:, 2]-x[:, 1], x[:, 0]-x[:, 2], x[:, 1]-x[:, 0]], axis=1)
            grad = np.stack([beta, gamma], axis=1) / (2 * a_safe)[:, None, None]
            i_e = -np.einsum('exn,en->ex', grad, h_el)
            i_e[~ok] = 0.0
            np.add.at(gradient, conn.ravel(), np.repeat(i_e, 3, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(ok.astype(float), 3))

        elif et in (4, 8, 9):
            i_sum = np.zeros((len(idx), 2))
            n_ok = np.zeros(len(idx))
            for xi, eta, _w in _quad_gauss(et):
                dxi, deta = _quad_dshape(et, xi, eta)
                J00 = coords[:, :, 0] @ dxi
                J01 = coords[:, :, 1] @ dxi
                J10 = coords[:, :, 0] @ deta
                J11 = coords[:, :, 1] @ deta
                detJ = J00*J11 - J01*J10
                ok = detJ > 0
                d_safe = np.where(ok, detJ, 1.0)
                dN_dx = (J11[:, None]*dxi[None, :] - J01[:, None]*deta[None, :]) / d_safe[:, None]
                dN_dy = (-J10[:, None]*dxi[None, :] + J00[:, None]*deta[None, :]) / d_safe[:, None]
                i_gp = -np.stack([(dN_dx*h_el).sum(axis=1), (dN_dy*h_el).sum(axis=1)], axis=1)
                i_gp[~ok] = 0.0
                i_sum += i_gp
                n_ok += ok
            np.add.at(gradient, conn.ravel(), np.repeat(i_sum, nn, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(n_ok, nn))

        elif et == 6:
            x, y = coords[:, :, 0], coords[:, :, 1]
            detJ = (x[:, 0]-x[:, 2])*(y[:, 1]-y[:, 2]) - (x[:, 1]-x[:, 2])*(y[:, 0]-y[:, 2])
            ok = np.abs(detJ) > 1e-10
            d_safe = np.where(ok, detJ, 1.0)
            Ji00 = (y[:, 1]-y[:, 2]) / d_safe
            Ji01 = -(y[:, 0]-y[:, 2]) / d_safe
            Ji10 = -(x[:, 1]-x[:, 2]) / d_safe
            Ji11 = (x[:, 0]-x[:, 2]) / d_safe
            w_total = np.zeros(len(idx))
            i_sum = np.zeros((len(idx), 2))
            for (L1, L2, L3), w in zip([(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)],
                                       [1/3, 1/3, 1/3]):
                dN_dL1 = np.array([4*L1-1, 0, 0, 4*L2, 0, 4*L3])
                dN_dL2 = np.array([0, 4*L2-1, 0, 4*L1, 4*L3, 0])
                dN_dL3 = np.array([0, 0, 4*L3-1, 0, 4*L2, 4*L1])
                dxi = dN_dL1 - dN_dL3
                deta = dN_dL2 - dN_dL3
                dN_dx = Ji00[:, None]*dxi[None, :] + Ji01[:, None]*deta[None, :]
                dN_dy = Ji10[:, None]*dxi[None, :] + Ji11[:, None]*deta[None, :]
                i_gp = -np.stack([(dN_dx*h_el).sum(axis=1), (dN_dy*h_el).sum(axis=1)], axis=1)
                i_gp[~ok] = 0.0
                i_sum += w * i_gp
                w_total += np.where(ok, w, 0.0)
            np.add.at(gradient, conn.ravel(), np.repeat(i_sum, 6, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(w_total, 6))

        else:
            raise ValueError(
                f"compute_gradient has no branch for element type {int(et)} "
                f"({len(idx)} elements); it recovers the gradient on element "
                f"types 3 (tri3), 4 (quad4), 6 (tri6), 8 (quad8) and 9 "
                f"(quad9). Falling through would leave those nodes at zero "
                f"gradient, which reads as a flat head field rather than as a "
                f"missing branch.")

    count[count == 0] = 1  # Avoid division by zero
    gradient /= count[:, None]
    return gradient


def run_seepage_analysis(seep_data, tol=1e-6, closure_tol=1e-3, max_iter=400):
    """
    Standalone function to run seep analysis.

    Args:
        seep_data: Dictionary containing all the seep data. When it was built by
            :func:`build_seep_data` with ``check_inputs=True`` (the default) it
            carries that model's preflight findings, and this function refuses to
            solve on an error -- a conductivity of zero, a boundary set that drives
            no flow, a mesh built against a different material table. The gate lands
            here rather than at the build because building a seep_data is not always
            a run: :func:`import_seep_solution` takes one purely as the shape of the
            mesh it needs to read a stored field back.
        tol: head-change tolerance, as a fraction of the model's head scale (the
            larger of the mesh height and the range of the specified heads), so the
            gate it sets does not move with the elevation datum
        closure_tol: nonlinear-residual tolerance for unconfined problems —
            iteration continues until the L1 nodal residual of the un-relaxed step
            over the free rows, relative to the inflow, falls below this. It is a
            measure of how far the kr matrix lags the head field, not of boundary
            mass balance; the balance the solution actually closed to is reported
            separately, as ``closure_error`` and ``closure_fraction``
        max_iter: iteration cap for the unconfined (exit-face) solver; raise it
            when a hard problem reports non-convergence near the cap
    
    Returns:
        Dictionary containing solution results with the following keys:
        - 'head': numpy array of hydraulic head values at each node
        - 'u': numpy array of pore pressure values at each node
        - 'velocity': numpy array of shape (n_nodes, 2) containing velocity vectors [vx, vy] at each node
        - 'gradient': numpy array of shape (n_nodes, 2) containing hydraulic gradient vectors [ix, iy] at each node
        - 'v_mag': numpy array of velocity magnitude at each node
        - 'i_mag': numpy array of hydraulic gradient magnitude at each node
        - 'q': numpy array of nodal flow vector
        - 'phi': numpy array of stream function/flow potential values at each node
        - 'flowrate': scalar total flow rate
        - 'flux_nodal': numpy array of consistent nodal loads applied by the
          specified-flux (Neumann) BCs (+ = inflow), zero without flux BCs
        - 'converged': whether the unconfined iteration closed (a confined solve is
          direct, and always True)
        - 'closure_error': |net inflow - net outflow| at the solution, an ABSOLUTE
          flow rate in the model's own flow units
        - 'closure_fraction': the same imbalance divided by the flow through the
          section — the dimensionless reading of it
    """
    report = seep_data.get("_preflight")
    if report is not None:
        report.raise_for_errors()

    # A boundary value bound to a tseep time series has no meaning in a steady
    # solve — there is no time axis to read the series at, and silently leaving
    # those nodes free would solve a different problem than the file describes.
    # Refuse loudly instead. The steady solution at the initial pool level is the
    # transient run's own first frame (t = 0), so nothing is lost by refusing.
    # A boundary value bound to a tseep series is read AT ITS t = 0 VALUE for a
    # steady solve — the same values the transient march's own initial condition
    # uses, so the steady run is the initial-condition snapshot without the march
    # behind it. Never silent: each reading is announced on the log. The overlay
    # works on copies, so the caller's seep_data stays untouched (its arrays keep
    # the binding nodes free, as the transient path expects).
    _hb = seep_data.get("head_series_bindings") or []
    _fb = seep_data.get("flux_series_bindings") or []
    _bc_overlay = None
    if _hb or _fb:
        _unresolved = [b["series"] for b in (_hb + _fb)
                       if b.get("value_t0") is None]
        if _unresolved:
            _names = ", ".join(f"'{s}'" for s in dict.fromkeys(_unresolved))
            raise SeepInputError(
                f"Cannot run a STEADY seepage analysis: boundary values are bound "
                f"to time series ({_names}) whose t = 0 values are not carried on "
                f"this seep_data. Rebuild it with build_seep_data from the loaded "
                f"model (so the tseep sheet is in reach), or replace the series "
                f"name with a number in the seep bc sheet.")
        _nodes_y = np.asarray(seep_data["nodes"], dtype=float)[:, 1]
        _bt = np.asarray(seep_data["bc_type"]).copy()
        _bv = np.asarray(seep_data["bc_values"], dtype=float).copy()
        _fn = seep_data.get("flux_nodal")
        _fn = (np.zeros(len(_bt)) if _fn is None
               else np.asarray(_fn, dtype=float).copy())
        for _b in _hb:
            _v = float(_b["value_t0"])
            _mask = _b["mask"]
            print(f"Series '{_b['series']}' read at t = 0 for the steady solve: "
                  f"{_b.get('kind', 'head')} value {_v:g}")
            if str(_b.get("kind", "head")).strip().lower() == "reservoir":
                # At the t = 0 level, exactly as a constant reservoir bakes: hold
                # the submerged nodes. A series read at t = 0 has no history, so
                # the nodes above the level stay no-flow (or whatever another
                # boundary makes them).
                _apply_reservoir_level(_bt, _bv, _mask, _nodes_y, _v)
            else:
                _bt[_mask] = 1
                _bv[_mask] = _v
        for _b in _fb:
            _v = float(_b["value_t0"])
            print(f"Series '{_b['series']}' read at t = 0 for the steady solve: "
                  f"flux value {_v:g}")
            _fn = _fn + np.asarray(_b["unit_flux_nodal"], dtype=float) * _v
        _bc_overlay = (_bt, _bv, _fn)

    start_time = time.time()

    # Missing unsaturated parameters: raise rather than return None. Returning None
    # only pushed the failure downstream (export/plot index into the solution and
    # crash with a cryptic TypeError); a raised SeepInputError carries a clear
    # message that the CLI and Studio can surface directly.
    if seep_data.get("missing_unsat_params", False):
        detail = seep_data.get("missing_unsat_detail") or []
        msg = ("Cannot run seepage analysis: one or more materials are missing the "
               "unsaturated parameters required for unconfined seepage with an "
               "exit-face boundary condition.")
        if detail:
            msg += "\n\n" + "\n".join(detail)
        msg += ("\n\nSet, per material, either kr0 (>0) and h0 (<0) for the "
                "linear-front model, or vg_a (>0) and vg_n (>1) for the van "
                "Genuchten model.")
        raise SeepInputError(msg)

    # Extract data from seep_data
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    if _bc_overlay is not None:
        # Series-bound model: the t = 0 overlay computed above.
        bc_type, bc_values, flux_nodal = _bc_overlay
    else:
        bc_type = seep_data["bc_type"]
        bc_values = seep_data["bc_values"]
        # Consistent nodal loads from specified-flux BCs (absent on seep_data
        # built by import_seep2d or older callers -> no flux).
        flux_nodal = seep_data.get("flux_nodal")
        if flux_nodal is None:
            flux_nodal = np.zeros(len(bc_type))
        else:
            flux_nodal = np.asarray(flux_nodal, dtype=float)
    has_flux = bool(np.any(flux_nodal != 0.0))
    if has_flux and not np.any(bc_type > 0):
        raise SeepInputError(
            "Seepage problem has specified-flux boundary conditions but no "
            "Dirichlet boundary (no specified head and no exit face); the head is "
            "then defined only up to an additive constant and the system is "
            "singular. Add a specified-head or exit-face boundary condition."
        )
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    k1_by_mat = seep_data["k1_by_mat"]
    k2_by_mat = seep_data["k2_by_mat"]
    angle_by_mat = seep_data["angle_by_mat"]
    kr0_by_mat = seep_data["kr0_by_mat"]
    h0_by_mat = seep_data["h0_by_mat"]
    # Per-material unsaturated model + van Genuchten params (absent on seep_data
    # built before vG support / by import_seep2d -> default to linear front).
    unsat_by_mat = seep_data.get("unsat_by_mat")
    vg_a_by_mat = seep_data.get("vg_a_by_mat")
    vg_n_by_mat = seep_data.get("vg_n_by_mat")
    vg_l_by_mat = seep_data.get("vg_l_by_mat")
    unit_weight = seep_data["unit_weight"]
    
    # Determine if unconfined flow
    is_unconfined = np.any(bc_type == 2)
    flow_type = "unconfined" if is_unconfined else "confined"
    print(f"Solving {flow_type.upper()} seep problem...")
    print("Number of fixed-head nodes:", np.sum(bc_type == 1))
    print("Number of exit face nodes:", np.sum(bc_type == 2))

    # Dirichlet BCs: fixed head (bc_type == 1) and possibly exit face (bc_type == 2)
    bcs = [(i, bc_values[i]) for i in range(len(bc_type)) if bc_type[i] in (1, 2)]

    # Material properties (per element)
    mat_ids = element_materials - 1
    k1 = k1_by_mat[mat_ids]
    k2 = k2_by_mat[mat_ids]
    angle = angle_by_mat[mat_ids]

    # Solve for head, stiffness matrix A, and nodal flow vector q
    if is_unconfined:
        # Get kr0/h0 and the unsaturated-model arrays per element based on material
        kr0_per_element = kr0_by_mat[mat_ids]
        h0_per_element = h0_by_mat[mat_ids]
        model_per_element = None if unsat_by_mat is None else unsat_by_mat[mat_ids]
        vg_a_per_element = None if vg_a_by_mat is None else vg_a_by_mat[mat_ids]
        vg_n_per_element = None if vg_n_by_mat is None else vg_n_by_mat[mat_ids]
        vg_l_per_element = None if vg_l_by_mat is None else vg_l_by_mat[mat_ids]

        head, A, q, total_flow, exit_face_active, converged, closure_error = solve_unsaturated(
            nodes=nodes,
            elements=elements,
            bc_type=bc_type,
            bc_values=bc_values,
            kr0=kr0_per_element,
            h0=h0_per_element,
            k1_vals=k1,
            k2_vals=k2,
            angles=angle,
            element_types=element_types,
            tol=tol,
            max_iter=max_iter,
            closure_tol=closure_tol,
            vg_a=vg_a_per_element,
            vg_n=vg_n_per_element,
            vg_l=vg_l_per_element,
            model=model_per_element,
            flux_nodal=flux_nodal,
        )
        # Compute phi BCs from element-level boundary flux
        dirichlet_phi_bcs = create_flow_potential_bc_from_elements(
            nodes, elements, element_types, head, k1, k2, angle,
            kr0=kr0_per_element, h0=h0_per_element, total_flow=total_flow,
            bc_type=bc_type, exit_face_active=exit_face_active,
            vg_a=vg_a_per_element, vg_n=vg_n_per_element, vg_l=vg_l_per_element,
            model=model_per_element)
        phi = solve_flow_function_unsaturated(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, dirichlet_phi_bcs, element_types,
                                              vg_a=vg_a_per_element, vg_n=vg_n_per_element, vg_l=vg_l_per_element, model=model_per_element)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, element_types,
                                    vg_a=vg_a_per_element, vg_n=vg_n_per_element, vg_l=vg_l_per_element, model=model_per_element)
    else:
        # Confined analysis is a single direct linear solve — always "converged".
        converged, closure_error = True, 0.0
        head, A, q, total_flow = solve_confined(nodes, elements, bc_type, bcs, k1, k2, angle,
                                                element_types, flux_nodal=flux_nodal)
        dirichlet_phi_bcs = create_flow_potential_bc_from_elements(
            nodes, elements, element_types, head, k1, k2, angle,
            total_flow=total_flow, bc_type=bc_type)
        phi = solve_flow_function_confined(nodes, elements, k1, k2, angle, dirichlet_phi_bcs, element_types)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, element_types=element_types)

    # Compute hydraulic gradient i = -grad(h)
    gradient = compute_gradient(nodes, elements, head, element_types)

    # Compute velocity and gradient magnitudes
    v_mag = np.linalg.norm(velocity, axis=1)
    i_mag = np.linalg.norm(gradient, axis=1)

    gamma_w = unit_weight
    u = gamma_w * (head - nodes[:, 1])

    # Ponding check. A specified inflow the soil cannot accept drives the surface
    # pressure positive; the physical response is ponding (the BC would switch to a
    # specified head), which is not modelled here — warn that the result is suspect.
    # Only meaningful on an unconfined problem: a confined domain is saturated by
    # construction and positive pressure there says nothing about ponding.
    if has_flux and is_unconfined:
        # Test the free nodes of the RUNTIME Dirichlet set, not bc_type == 0. An
        # INACTIVE exit-face node is free, still carries its load, and is exactly
        # where an over-specified inflow shows up — filtering on bc_type would skip
        # the nodes the check exists to catch. (Loads on ACTIVE exit nodes are
        # discarded by the solve, and a prescribed p = 0 there cannot pond anyway.)
        _dir_now = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
        ponded = np.where((flux_nodal > 0) & ~_dir_now & (u > 0))[0]
        if len(ponded) > 0:
            warnings.warn(
                f"{len(ponded)} specified-flux node(s) finished with positive pore "
                "pressure (max u = "
                f"{float(np.max(u[ponded])):.4g}): the specified inflow exceeds what "
                "the soil can accept there, so in reality the surface would pond and "
                "the boundary would become a specified head. This solution is suspect.",
                stacklevel=2,
            )

    solution = {
        "head": head,
        "u": u,
        "velocity": velocity,
        "gradient": gradient,
        "v_mag": v_mag,
        "i_mag": i_mag,
        "q": q,
        "phi": phi,
        "flowrate": total_flow,
        # Consistent nodal loads applied by the specified-flux BCs (+ = inflow).
        "flux_nodal": flux_nodal,
        "converged": converged,
        # Two readings of the same imbalance, because one of them alone is not
        # readable. ``closure_error`` is |net inflow - net outflow| — an absolute
        # flow rate, in whatever flow units the model is in, which says nothing on
        # its own about whether it is large. ``closure_fraction`` is that imbalance
        # over the flow through the section, which is the number a reader can judge
        # (and the one comparable with ``closure_tol``). A solve with no flow at all
        # has no fraction to report and records 0.0.
        "closure_error": closure_error,
        "closure_fraction": (abs(closure_error) / abs(total_flow)
                             if total_flow else 0.0),
        # Which branch produced this solution. Consumers need it: a confined solve is
        # fully saturated with kr never evaluated, so its negative pore pressures carry
        # no phreatic surface (see plot_seep_solution).
        "unconfined": bool(is_unconfined),
    }

    if not converged:
        print("WARNING: seepage solution did not converge — flowrate is unreliable "
              "(solution['converged'] is False).")

    elapsed = time.time() - start_time
    print(f"Seepage analysis completed in {elapsed:.2f} seconds.")

    return solution


# =============================================================================
# Transient (time-dependent) seepage
#
# plan_transient_seep.md §1 (formulation), §3 (architecture), §5 (time stepping).
# The transient solver slots into the SAME batched build-once assembly the steady
# solver uses: the saturated stiffness stacks and kr-sampling operators are built
# once (_build_assembly), the unit mass matrices are built once
# (_build_mass_assembly), and each theta-scheme Picard iteration reduces to a kr
# scaling of K and a storage scaling of M, exactly parallel. Steady behavior is
# untouched (this is an additive entry point; run_seepage_analysis is unchanged).
# =============================================================================


def _eval_series(anchors, t):
    """Evaluate one time series at time t. ``anchors`` is (t_arr, v_arr) of the
    series' own non-blank breakpoints (ascending t). Linear between anchors; holds
    the nearest value before the first and after the last anchor; right-continuous
    at a repeated time (a step function)."""
    ta, va = anchors
    n = len(ta)
    if n == 0:
        return 0.0
    if t <= ta[0]:
        return float(va[0])
    if t >= ta[-1]:
        return float(va[-1])
    # index of the last anchor <= t (for a repeated time, the LATER one -> the new
    # value applies from t onward, right-continuous)
    i = int(np.searchsorted(ta, t, side='right')) - 1
    i = min(max(i, 0), n - 2)
    t0, t1 = ta[i], ta[i + 1]
    if t1 == t0:
        return float(va[i + 1])
    return float(va[i] + (va[i + 1] - va[i]) * (t - t0) / (t1 - t0))


def _series_high_water(anchors, t):
    """The highest value a series has reached over the run so far, [0, t].

    The series is piecewise linear on its anchors (:func:`_eval_series`), so its
    running maximum is the largest of the value at t = 0, the value at t, and every
    anchor value at an anchor time in [0, t]. At t = 0 the run has no history, and
    the high water is the t = 0 value itself, so a step anchored at t = 0 (two
    anchors at the same time) counts only once the run has moved past it."""
    ta, va = anchors
    high = max(_eval_series(anchors, 0.0), _eval_series(anchors, t))
    if t > 0.0:
        for tt, vv in zip(ta, va):
            if 0.0 <= tt <= t and vv > high:
                high = float(vv)
    return float(high)


def _apply_reservoir_level(bc_type, bc_values, mask, y, level, high_water=None):
    """Resolve one reservoir boundary in place, at a water level and the highest
    level it has reached so far in the run.

    A node of the boundary (``mask``) at or below ``level`` is submerged and held at
    the level (``bc_type`` 1). A node above the level that the water has reached
    (``y <= high_water``) is a seepage-exit face (``bc_type`` 2, value = its own
    elevation), unless another boundary already holds it at a head. A node above
    ``high_water`` has never been under water and is left as it stands: no-flow,
    or whatever another boundary makes it. ``high_water`` defaults to ``level``,
    which is the steady case (a constant level, or a series read at t = 0, has no
    history), and there no exit face arises."""
    if high_water is None:
        high_water = level
    submerged = mask & (y <= level)
    bc_type[submerged] = 1
    bc_values[submerged] = level
    uncovered = mask & (y > level) & (y <= high_water) & (bc_type != 1)
    bc_type[uncovered] = 2
    bc_values[uncovered] = y[uncovered]


def resolve_transient_bc(seep_data, tseep_data, t):
    """The boundary configuration a transient run applies at time ``t``.

    Returns ``(bc_type, bc_values, flux_nodal)`` as new arrays; ``seep_data`` is not
    modified. Constant (numeric) heads, exit faces and fluxes carry through from the
    baked arrays. Each time-varying head series is resolved by its kind:

    * ``"head"``: a plain Dirichlet at the series value h(t) at every node of the
      polyline, at all times.
    * ``"reservoir"``: a node at or below h(t) is held at h(t). A node above h(t)
      that the water has reached at some time in [0, t] (its elevation is at or
      below the running maximum of the series) is a seepage-exit face. A node above
      that maximum has never been under water and is no-flow, unless another
      boundary governs it.

    Each time-varying flux adds its unit load scaled by the series value at t.
    This is the resolution :func:`run_transient_seepage` applies at its initial
    condition and at every time step."""
    nodes = np.asarray(seep_data["nodes"], dtype=float)
    base_flux = seep_data.get("flux_nodal")
    base_flux = (np.zeros(len(nodes)) if base_flux is None
                 else np.asarray(base_flux, dtype=float))
    return _resolve_series_bc(
        nodes[:, 1], np.asarray(seep_data["bc_type"]),
        np.asarray(seep_data["bc_values"], dtype=float), base_flux,
        seep_data.get("head_series_bindings", []) or [],
        seep_data.get("flux_series_bindings", []) or [],
        tseep_data["series"], t)


def _resolve_series_bc(y, base_bc_type, base_bc_values, base_flux, head_bindings,
                       flux_bindings, series, t):
    """:func:`resolve_transient_bc` on the unpacked arrays (the form the march
    holds them in)."""
    def _anchors(name):
        if name not in series:
            raise SeepInputError(
                f"BC references tseep series {name!r}, not found in tseep_data. "
                f"Available: {', '.join(sorted(series)) or '(none)'}.")
        return series[name]

    bt = np.asarray(base_bc_type).copy()
    bv = np.asarray(base_bc_values, dtype=float).copy()
    for b in head_bindings:
        anchors = _anchors(b["series"])
        ht = _eval_series(anchors, t)
        mask = b["mask"]
        if str(b.get("kind", "reservoir")).strip().lower() == "head":
            # plain Dirichlet at h(t) everywhere on the polyline
            bt[mask] = 1
            bv[mask] = ht
            continue
        _apply_reservoir_level(bt, bv, mask, y, ht, _series_high_water(anchors, t))
    flux = np.asarray(base_flux, dtype=float).copy()
    for b in flux_bindings:
        flux = flux + _eval_series(_anchors(b["series"]), t) * b["unit_flux_nodal"]
    return bt, bv, flux


def build_tseep_data(slope_data):
    """Build the transient-seepage control/series bundle from ``slope_data['tseep']``.

    Returns ``None`` when the file has no ``tseep`` sheet (a steady model). Otherwise
    a dict carrying the time-series evaluator inputs (each series reduced to its own
    non-blank anchors), the run controls (duration, save_interval, save_times, stage
    times) and the derived series breakpoints. The per-BC series bindings live in
    ``seep_data`` (recorded by ``build_seep_data`` against the mesh node masks);
    ``run_transient_seepage`` consumes both.
    """
    ts = slope_data.get('tseep')
    if ts is None:
        return None
    times = np.asarray(ts['times'], dtype=float)
    series = {}
    for name, vals in ts.get('series', {}).items():
        ta, va = [], []
        for t, v in zip(times, vals):
            if v is not None:
                ta.append(float(t))
                va.append(float(v))
        series[name] = (np.asarray(ta, dtype=float), np.asarray(va, dtype=float))
    breakpoints = sorted({float(t) for ta, _ in series.values() for t in ta})
    return {
        'times': times,
        'series': series,
        'duration': ts.get('duration'),
        'save_interval': ts.get('save_interval'),
        'save_times': list(ts.get('save_times') or []),
        'stage_1': ts.get('stage_1'),
        'stage_2': ts.get('stage_2'),
        # The instant a stability run reads out of this march. It is not a control on
        # the march itself, but the stepper still has to LAND on it, so it joins the
        # saved-frame schedule alongside the stage times.
        'stability_time': ts.get('stability_time'),
        'breakpoints': breakpoints,
    }


def _transient_saved_schedule(tseep_data, duration):
    """Saved-frame schedule = UNION(save_interval grid, save_times, stage times,
    the stability time, series breakpoints), clamped to (0, duration], sorted and
    de-duplicated. t=0 is handled separately (the initial condition). The adaptive
    stepper lands exactly on every entry, so saved frames are computed states, never
    interpolated."""
    save_interval = tseep_data.get('save_interval')
    if not save_interval or save_interval <= 0:
        save_interval = duration / 50.0
    targets = set()
    k = 1
    while k * save_interval < duration - 1e-12:
        targets.add(k * save_interval)
        k += 1
    for t in tseep_data.get('save_times', []):
        if 0 < t <= duration:
            targets.add(float(t))
    for key in ('stage_1', 'stage_2', 'stability_time'):
        st = tseep_data.get(key)
        if st is not None and 0 < st <= duration:
            targets.add(float(st))
    for bp in tseep_data.get('breakpoints', []):
        if 0 < bp <= duration:
            targets.add(float(bp))
    targets.add(float(duration))
    return sorted(targets), save_interval


def _transient_system(asm, k_data, m_diag_over_dt, rhs, dir_mask, dir_values):
    """Assemble the BC-applied theta-scheme system (A, b).

    A = theta*K + diag(M/dt) with Dirichlet rows replaced by identity; b = rhs with
    Dirichlet rows overwritten by the prescribed head. Built directly from the COO
    index arrays plus n diagonal mass entries, mirroring ``_dirichlet_system``."""
    n = asm['n_nodes']
    diag_idx = np.arange(n)
    rows = np.concatenate([asm['rows'], diag_idx])
    cols = np.concatenate([asm['cols'], diag_idx])
    vals = np.concatenate([k_data, m_diag_over_dt])
    keep = ~dir_mask[rows]
    di = np.where(dir_mask)[0]
    rows2 = np.concatenate([rows[keep], di])
    cols2 = np.concatenate([cols[keep], di])
    vals2 = np.concatenate([vals[keep], np.ones(len(di))])
    A = coo_matrix((vals2, (rows2, cols2)), shape=(n, n)).tocsr()
    b = np.asarray(rhs, dtype=float).copy()
    b[di] = dir_values[di]
    return A, b


def _exit_face_topology(elements, element_types, bc_type):
    """Exit-face (``bc_type==2``) boundary topology for the active-set update.

    Returns ``(is_corner, linear_corners, quadratic_edges)``:

    * ``is_corner``       boolean mask of exit-face CORNER nodes;
    * ``linear_corners``  boolean mask of corners on linear (tri3/quad4) exit
                          edges — these use the per-corner SEEP2D h/q test;
    * ``quadratic_edges`` list of ``(corner1, midside, corner2)`` triplets for
                          quadratic (tri6/quad8/quad9) exit edges, tracked
                          all-or-nothing so a seepage-face transition lands on a
                          corner, never in the middle of a curved edge.

    An exit edge is one whose BOTH corners are exit-face nodes, or a BOUNDARY edge
    with one exit-face corner and one specified-head corner. The second kind is the
    last edge of a seepage face where it meets a head line at a shared corner — the
    downstream toe, which the specified head claims (see build_seep_data). Without
    it that edge belongs to no exit edge at all: its midside node has no state to be
    given, its upper corner is tied to the dry edge above, and the face ends one
    edge short of the toe with the head standing above the ground there. Measured on
    the Johnson Reservoir dam (tri6): 0.66 to 1.77 ft of pressure head ON the face
    over the last 10 ft above the toe, 41 to 110 psf of pore pressure that the stress
    analysis then has to hold as effective tension at a free surface. The fixed
    corner is Dirichlet already, so in the active-set test it counts as saturated
    and is never itself switched (see the edge loops in the two solvers).

    This is the single source of the exit-face edge topology used by BOTH the steady
    unconfined solver (``solve_unsaturated``) and the transient stepper
    (``run_transient_seepage``); the transient rebuilds it per step because a
    series-head boundary moves the exit-face set as the water level changes."""
    n = len(bc_type)
    is_corner = np.zeros(n, dtype=bool)
    linear_corners = np.zeros(n, dtype=bool)
    quadratic_edges = []
    seen = set()
    # Element-edge multiplicity, so an exit/head edge is admitted only where it is
    # a boundary edge (one element on it). Two exit-face corners are on the face by
    # definition; an exit corner and a head corner can also be joined by an INTERIOR
    # edge of a sliver element spanning the toe, and that edge must not become a
    # seepage face.
    edge_count = {}
    for idx, en in enumerate(elements):
        et = element_types[idx]
        if et == 3:
            _le = [(0, 1), (1, 2), (2, 0)]
        elif et == 6:
            _le = [(0, 1), (1, 2), (2, 0)]
        elif et == 4:
            _le = [(0, 1), (1, 2), (2, 3), (3, 0)]
        elif et in (8, 9):
            _le = [(0, 1), (1, 2), (2, 3), (3, 0)]
        else:
            continue
        for a, b in _le:
            key = tuple(sorted((int(en[a]), int(en[b]))))
            edge_count[key] = edge_count.get(key, 0) + 1
    for idx, en in enumerate(elements):
        et = element_types[idx]
        if et == 3:
            ledges = [(0, 1), (1, 2), (2, 0)]
        elif et == 6:
            ledges = [(0, 3, 1), (1, 4, 2), (2, 5, 0)]
        elif et == 4:
            ledges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        elif et in (8, 9):
            ledges = [(0, 4, 1), (1, 5, 2), (2, 6, 3), (3, 7, 0)]
        else:
            continue
        for le in ledges:
            gnodes = [en[i] for i in le]
            c1, c2 = gnodes[0], gnodes[-1]
            key = tuple(sorted((int(c1), int(c2))))
            if key in seen:
                continue
            both_exit = bc_type[c1] == 2 and bc_type[c2] == 2
            exit_to_head = ({int(bc_type[c1]), int(bc_type[c2])} == {1, 2}
                            and edge_count.get(key, 0) == 1)
            if both_exit or exit_to_head:
                seen.add(key)
                for c in (c1, c2):
                    if bc_type[c] == 2:
                        is_corner[c] = True
                if len(gnodes) == 3:  # quadratic edge: has a midside node
                    quadratic_edges.append((int(c1), int(gnodes[1]), int(c2)))
                else:
                    for c in (c1, c2):
                        if bc_type[c] == 2:
                            linear_corners[c] = True
    return is_corner, linear_corners, quadratic_edges


def _resolve_exit_cycle(act, new_act, corner_candidate, mid_cands,
                        quadratic_edges, linear_corners, pinned, pinned_state):
    """Resolve a cycling exit-face active set, in place into ``pinned`` /
    ``pinned_state``.

    A quadratic (tri6/quad8/quad9) exit edge is tracked ALL-OR-NOTHING, so a partly
    wet edge -- the state the face is genuinely in while the phreatic exit point
    descends through that edge -- has no representation at all. Both available
    states are then inadmissible: held wet the edge over-drains and the rule sheds
    it, held dry the pressure climbs above the face and the rule takes it back. The
    set alternates forever and the Picard set-stability test can never close.

    The resolution is to give the offending edge the missing state. Every node on an
    edge that is oscillating is set to its OWN SEEP2D h/q candidate -- the per-node
    rule linear elements always use, which CAN carry a wet lower half and a dry
    upper half -- and pinned there for the rest of the step, so the seepage-face
    transition lands at the midside node instead of sweeping the whole edge.
    Untouched edges keep the all-or-nothing rule, which is what makes the transition
    land cleanly on a corner everywhere else.

    Parameters
    ----------
    act, new_act : the two alternating active sets (boolean, per node).
    corner_candidate : each exit corner's own h/q test this sweep.
    mid_cands : the midside node's own h/q test, one per entry of ``quadratic_edges``.
    quadratic_edges : ``(corner1, midside, corner2)`` triplets.
    linear_corners : boolean mask of corners on linear exit edges.
    pinned, pinned_state : boolean arrays updated in place (nodes already pinned by
        an earlier resolution in the same step keep their state).
    """
    osc = new_act != act
    for (c1, mid, c2), mc in zip(quadratic_edges, mid_cands):
        if not (osc[c1] or osc[mid] or osc[c2]):
            continue
        for nd, cand in ((c1, bool(corner_candidate[c1])),
                         (mid, bool(mc)),
                         (c2, bool(corner_candidate[c2]))):
            pinned[nd] = True
            pinned_state[nd] = cand
    lin_osc = osc & linear_corners
    pinned[lin_osc] = True
    pinned_state[lin_osc] = corner_candidate[lin_osc]
    return pinned, pinned_state


#: Time-error tolerance of the transient stepper, as a fraction of the
#: characteristic head span: the estimated local error of a backward-Euler step at
#: any free node is held below it (see ``run_transient_seepage``).
TSEEP_TIME_ERROR_FRAC = 1e-3


def run_transient_seepage(seep_data, tseep_data, theta=1.0, lumped=True,
                          dt_max=None, growth=1.3,
                          max_head_change_frac=0.05, h_init=None, verbose=True,
                          progress_callback=None, max_iter=2000,
                          time_error_frac=None):
    """Transient variably-saturated seepage solver (plan §1/§3/§5).

    Solves ``div(kr K grad h) + Q = S dh/dt`` with the theta-method in time
    (default backward Euler, ``theta=1``), reusing the steady batched assembly.
    Storage S is the elastic Ss where saturated and the retention curve's specific
    moisture capacity C(psi) where unsaturated — one or the other, never their sum
    (``storage_capacity_vec``), the same saturated-only convention RS2 and SEEP/W
    use. The storage term of a step is the change in each element's stored water
    between the two heads (``storage_potential_vec`` summed over the element's
    storage points; see ``_storage_residual``), so the march conserves water
    whatever the step size. Each step is solved by iteration: kr is taken from the
    latest iterate and the storage term is linearized on a diagonal, the larger of
    the chord and the tangent capacity at each storage point.
    Time-varying BCs come from the ``tseep`` series bound
    in ``seep_data``, resolved per step by the head binding's kind
    (:func:`resolve_transient_bc`): a "reservoir" series holds a boundary node at
    h(t) while it is submerged, makes it an exit-face node once the level falls
    back below it, and leaves it no-flow until the water first reaches it, while a
    "head" series is a plain Dirichlet at h(t) at every node of the polyline at all
    times (it can hold a suction/negative-pressure head and never converts to an
    exit face). Adaptive dt clamps to the saved-frame schedule so every saved frame is a
    computed state.

    Parameters
    ----------
    seep_data : dict from ``build_seep_data`` (carries ss_by_mat/sy_by_mat and the
        head/flux series bindings).
    tseep_data : dict from ``build_tseep_data``.
    theta : time-integration weight (1 = backward Euler, 0.5 = Crank-Nicolson).
    lumped : lumped (HRZ) mass by default; ``False`` is not implemented in v1.
    h_init : optional explicit initial head field; default is a steady solve at the
        t=0 BC configuration (plan IC rule).
    time_error_frac : time-error tolerance as a fraction of the characteristic head
        span (default :data:`TSEEP_TIME_ERROR_FRAC`). After each step the local
        error of backward Euler at every free node is estimated from the change in
        that node's rate of head change since the step before; a step estimated at
        more than twice the tolerance is redone shorter, and the next step is sized
        to land on it.
    max_iter : sweep budget for that initial-condition solve, the same budget the
        steady runners carry and the only convergence parameter of it a caller sets.
        It is spent once, before the march, and has no bearing on the time steps.
        An initial condition that does not close within it leaves ``converged``
        False on the returned dict, exactly as a time step that has to be
        force-accepted does.
    progress_callback : optional ``callback(t, duration) -> bool``. When supplied it
        is called once after each ACCEPTED time step with the new simulated time and
        the run duration (and once at completion with ``(duration, duration)``) so a
        caller can drive a determinate progress bar (fraction = ``t / duration``). It
        doubles as the cancellation channel: a FALSY return requests an abort, on
        which the march stops cleanly after the current step and the returned dict has
        ``cancelled=True`` (its ``frames`` hold only the states saved before the
        abort). ``None`` (the default) is zero-overhead and preserves today's
        behavior exactly (never cancelled).

    Returns
    -------
    dict with ``times`` (saved times), ``frames`` (list of {time, head, u}),
    ``dt_history``, ``mass_balance`` (cumulative + per-frame ledger), ``converged``
    (the initial-condition solve AND every time step; a force-accepted step or an
    unclosed IC clears it), ``cancelled`` (True only when a ``progress_callback``
    requested an abort), and provenance (theta, lumped, duration).
    """
    if not lumped:
        raise NotImplementedError(
            "run_transient_seepage: consistent-mass (lumped=False) is not implemented "
            "in v1; the lumped (HRZ) mass matrix is the default and only option.")
    if tseep_data is None:
        raise SeepInputError("run_transient_seepage requires tseep_data (this model "
                             "has no 'tseep' sheet).")
    duration = tseep_data.get('duration')
    if duration is None or duration <= 0:
        raise SeepInputError("Transient run needs a positive 'duration' in the tseep "
                             "controls.")

    nodes = seep_data["nodes"]
    elements = np.asarray(seep_data["elements"])
    element_types = seep_data.get("element_types")
    if element_types is None:
        element_types = np.full(len(elements), 3)
    element_types = np.asarray(element_types)
    element_materials = np.asarray(seep_data["element_materials"])
    y = nodes[:, 1]
    n_nodes = len(nodes)

    mat_ids = element_materials - 1
    k1 = seep_data["k1_by_mat"][mat_ids]
    k2 = seep_data["k2_by_mat"][mat_ids]
    angle = seep_data["angle_by_mat"][mat_ids]
    kr0 = seep_data["kr0_by_mat"][mat_ids]
    h0 = seep_data["h0_by_mat"][mat_ids]
    ss = seep_data["ss_by_mat"][mat_ids]
    sy = seep_data["sy_by_mat"][mat_ids]
    _unsat = seep_data.get("unsat_by_mat")
    _vga = seep_data.get("vg_a_by_mat")
    _vgn = seep_data.get("vg_n_by_mat")
    _vgl = seep_data.get("vg_l_by_mat")
    model = None if _unsat is None else _unsat[mat_ids]
    vg_a = None if _vga is None else _vga[mat_ids]
    vg_n = None if _vgn is None else _vgn[mat_ids]
    vg_l = None if _vgl is None else _vgl[mat_ids]
    gamma_w = seep_data["unit_weight"]

    base_bc_type = np.asarray(seep_data["bc_type"]).copy()
    base_bc_values = np.asarray(seep_data["bc_values"], dtype=float).copy()
    base_flux = seep_data.get("flux_nodal")
    base_flux = (np.zeros(n_nodes) if base_flux is None
                 else np.asarray(base_flux, dtype=float))
    head_bindings = seep_data.get("head_series_bindings", []) or []
    flux_bindings = seep_data.get("flux_series_bindings", []) or []
    series = tseep_data["series"]

    # An exit face can be present statically (a downstream seepage face) or appear
    # dynamically as a series-bound reservoir level falls below a node it has
    # submerged. The per-step active-set update (in ``attempt``) resolves the exit
    # face at any element order:
    # corners use the per-corner SEEP2D h/q test and quadratic (tri6/quad8/quad9)
    # exit edges are tracked all-or-nothing via ``_exit_face_topology`` — the same
    # edge logic the steady solver uses — so the seepage-face transition lands on a
    # corner rather than mid-edge.
    # Only a "reservoir" head series can dynamically raise an exit face, as its
    # level falls back below nodes it has submerged; a plain "head" series is a
    # Dirichlet at every node at all times and never converts. Static exit faces
    # (bc_type 2) always count.
    _reservoir_binding = any(
        str(b.get("kind", "reservoir")).strip().lower() != "head"
        for b in head_bindings)
    _may_have_exit = bool(np.any(base_bc_type == 2) or _reservoir_binding)

    def _series(name):
        if name not in series:
            raise SeepInputError(
                f"BC references tseep series {name!r}, not found in tseep_data. "
                f"Available: {', '.join(sorted(series)) or '(none)'}.")
        return series[name]

    def resolve_bc(t):
        """BC configuration at time t: (bc_type, bc_values, flux_nodal), as
        :func:`resolve_transient_bc` states it. A "reservoir" series holds its
        submerged nodes (elevation <= h(t)) at h(t), makes an exit face of each node
        the water has reached and since uncovered (h(t) < elevation <= the series'
        running maximum over [0, t]), and leaves a node the water has never reached
        no-flow. A "head" series is a plain Dirichlet at h(t) at every node of the
        polyline, at all times (it can hold a negative-pressure/suction head and
        never converts to an exit face)."""
        return _resolve_series_bc(y, base_bc_type, base_bc_values, base_flux,
                                  head_bindings, flux_bindings, series, t)

    # Build-once assembly (stiffness + mass).
    asm = _build_assembly(nodes, elements, element_types, k1, k2, angle)
    masm = _build_mass_assembly(nodes, elements, element_types)
    sasm = _build_storage_assembly(nodes, masm, ss, sy, h0, vg_a, vg_n, model)

    def kdata(p):
        return _assembly_data(asm, p_nodes=p, kr0=kr0, h0=h0, mode='head',
                              vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)

    def _stored(hh):
        """Water held in the whole mesh at head hh, measured from saturation."""
        return float(sum(W.sum() for W in _stored_water(sasm, hh - y)))

    def _frame_flows(rate, flux, dir_mask):
        """Boundary inflow and outflow (rates) for a saved frame.

        ``rate`` is the boundary flow at each fixed-head and active exit-face node
        over the step that ended at the frame: the node's flow into the mesh
        (``theta K h + (1 - theta) K_old h_old`` at that row) plus the water its own
        storage took up over the step, divided by the step length. For the initial
        condition, a steady state, it is ``K h``. Positive is inflow. The specified
        fluxes on the free nodes are counted with the boundary flows, inflow and
        outflow by sign, so inflow minus outflow is exactly the rate at which the
        water stored in the mesh changed over that step. A cheap by-product of the
        step's own assembly, not a dual solve. No stream function is computed — a
        transient (storage-release) state is not divergence-free, so it has no flow
        net; the frame's flow direction is read from velocity vectors instead."""
        r = rate[dir_mask]
        f_free = np.asarray(flux, dtype=float)[~dir_mask]
        inflow = float(np.sum(r[r > 0]) + np.sum(f_free[f_free > 0]))
        outflow = float(-np.sum(r[r < 0]) - np.sum(f_free[f_free < 0]))
        return inflow, outflow

    # ---- initial condition ----
    bt0, bv0, flux0 = resolve_bc(0.0)
    if h_init is not None:
        h = np.asarray(h_init, dtype=float).copy()
        ic_converged = True
    else:
        h, ic_converged = _transient_steady_ic(
            nodes, elements, element_types, bt0, bv0,
            k1, k2, angle, kr0, h0, vg_a, vg_n, vg_l, model, flux0, verbose,
            max_iter=max_iter)
        if not ic_converged:
            # Reported the way the steady path reports the same failure: the march
            # is starting from a field that is not a steady state, so every frame
            # downstream of it inherits that, and the flag says so.
            print("WARNING: the transient initial condition did not converge — the "
                  "transient run starts from a field that is not a steady state "
                  "(solution['converged'] is False).")

    # characteristic head scale for the Picard tolerance and the dh limiter
    hmin, hmax = float(np.min(h)), float(np.max(h))
    for b in head_bindings:
        for v in _series(b["series"])[1]:
            hmin, hmax = min(hmin, v), max(hmax, v)
    char_h = max(hmax - hmin, float(np.max(y) - np.min(y)), 1.0)
    picard_max = 25                # Picard iterations allowed within a step
    picard_eps = 1e-5 * char_h     # inner-loop tolerance, scaled to the head range
    dh_limit = max_head_change_frac * char_h
    # Time-error control: the local error of a backward-Euler step at a free node is
    # dt^2/2 h'', estimated from the change in the node's rate of head change between
    # the step and the one before it. Held below err_tol; a step estimated at more than
    # twice that is redone at a shorter step.
    err_tol = (TSEEP_TIME_ERROR_FRAC if time_error_frac is None
               else float(time_error_frac)) * char_h

    targets, save_interval = _transient_saved_schedule(tseep_data, duration)
    if dt_max is None:
        dt_max = save_interval
    first = targets[0] if targets else duration
    dt0 = max(min(first, save_interval) / 10.0, duration * 1e-9)
    dt_min = dt0 * 1e-6

    # Mass-balance ledger. The water held in the mesh is measured with the same
    # integral the steps conserve (_stored), and the boundary flow is accumulated
    # step by step from the accepted states themselves (see ``attempt``), so the
    # two are independent: their difference is the sum, over the accepted steps, of
    # whatever each state leaves unsolved in its own equations.
    stored0 = _stored(h)
    # Physical storage scale (total area x max Ss x head span). The relative
    # closure is normalized against the LARGER of the actual net cumulative inflow
    # and a small fraction of this scale, so a symmetric relaxation with ~zero NET
    # storage change does not produce a 0/0 ratio.
    stored_scale = (float(sum(sg['A_w'].sum() for sg in sasm['groups']))
                    * max(float(np.max(ss)), 1e-30) * char_h)

    # The initial condition is a steady state: its boundary flow is the nodal
    # reaction K h at the fixed-head nodes and at the exit-face nodes the initial
    # head field saturates (h >= y), which are the active exit set at t = 0. A
    # reservoir series contributes none of the latter: at t = 0 it has no history,
    # so each of its nodes is held (submerged) or no-flow. Confined problems have
    # no exit nodes.
    _act0 = (bt0 == 2) & (h >= y)
    _dir0 = (bt0 == 1) | _act0
    _kd0 = kdata(h - y) if bool(np.any(bt0 == 2)) else _assembly_data(asm)
    _in0, _out0 = _frame_flows(_coo_matvec(asm, _kd0, h), flux0, _dir0)
    frames = [{"time": 0.0, "head": h.copy(), "u": gamma_w * (h - y),
               "inflow": _in0, "outflow": _out0}]
    dt_history = []
    mb_frames = [{"time": 0.0, "stored_change": 0.0, "cumulative_inflow": 0.0,
                  "closure": 0.0}]
    cum_inflow = 0.0            # net boundary inflow, accumulated over accepted steps
    last_flow = [None]          # (rate, flux, dir_mask) of the last accepted step
    active = np.zeros(n_nodes, dtype=bool)   # exit-face active set, warm-started
    # Exit-face nodes of the last accepted step. A node that was already on the
    # exit face keeps its active/inactive state from that step as the warm start;
    # only a node that has just joined the face (a reservoir falling past it) is
    # seeded saturated. Empty before the first step, so the first step seeds the
    # whole face, as the steady solver starts it.
    prev_exit = [np.zeros(n_nodes, dtype=bool)]
    all_converged = ic_converged
    #: [cycles detected, edges resolved]. A cycle DETECTED is a sweep in which the
    #: exit-face active set revisited a state it had left; RESOLVED counts the times
    #: the partly-wet fallback actually had to be applied, which happens only on the
    #: retry of a step that already failed. detected > 0 with resolved == 0 is the
    #: ordinary case: the set wobbled and settled on its own. Lists so ``attempt``
    #: can increment them as closures.
    n_cycles = [0]
    n_damped = [0]
    #: Nodes the free-surface rule has joined to the seepage face over the march, and
    #: the pressure head it is there to prevent: the largest ``h - y`` left standing on
    #: an INACTIVE exit node at the end of an accepted step. The second is a property
    #: of the delivered field, not of the solver's path -- a free surface drains at
    #: atmospheric pressure, so anything above the hysteresis band is pore pressure the
    #: stress analysis would have to hold as tension at a traction-free boundary (the
    #: water.seep_pressure_on_free_surface preflight rule reads the same quantity off
    #: the staged field).
    n_wet = [0]
    wet_probe = [0.0]      # the last attempt's reading, recorded when a step is kept
    max_free_p = [0.0]

    def attempt(h_old, dt, t_new, warm_active, force=False, damp=False):
        """One theta-scheme step with the iteration + exit-face active set. Returns
        (ok, h_new, active, picard_count, step_flow, cycling), ``step_flow`` the
        accepted step's ledger entry ``(boundary_inflow, rate, flux, dir_mask,
        exit_mask)``.
        ``h_new`` is always the last computed iterate (not h_old), so the caller may
        accept a best-effort step. ``force=True`` accepts that iterate even if Picard
        did not converge or the head-change limiter would reject it — the last resort
        at dt_min so the stepper always makes progress. A genuinely singular system
        (empty Dirichlet set) still returns ok=False with h_old unchanged.

        ``cycling`` reports that the exit-face active set revisited a state it had
        already left — a limit cycle the set-stability test can never close.
        ``damp=True`` additionally RESOLVES it (see the anti-cycling block below);
        the caller turns that on only for a retry of a step that has already failed,
        so a march whose steps converge normally follows exactly the path it always
        did."""
        bt, bv, flux = resolve_bc(t_new)
        has_exit = bool(np.any(bt == 2))
        # Exit-face topology, rebuilt per step (a series-head boundary moves the
        # exit set as the water level changes). Corner nodes use the per-corner
        # SEEP2D h/q test; quadratic (tri6/quad8/quad9) exit edges are tracked
        # all-or-nothing via _exit_face_topology — the SAME edge machinery the
        # steady solver uses — so a seepage-face transition lands on a corner, not
        # mid-edge. Confined steps skip the per-element Python loop entirely.
        if has_exit:
            exit_is_corner, exit_linear_corners, exit_quadratic_edges = \
                _exit_face_topology(elements, element_types, bt)
            # nodes carrying an exit-edge state (corners + quadratic midsides)
            exit_on_edge = exit_is_corner.copy()
            for _c1, _mid, _c2 in exit_quadratic_edges:
                exit_on_edge[_mid] = True
        else:
            exit_is_corner = np.zeros(n_nodes, dtype=bool)
            exit_linear_corners = np.zeros(n_nodes, dtype=bool)
            exit_quadratic_edges = []
            exit_on_edge = exit_is_corner
        act = warm_active & (bt == 2)
        # seed newly-exit topology nodes (corners AND quadratic midsides) as active
        # (saturated) so the set can shed them — mirrors the steady solver starting
        # the whole exit face saturated. A node already on the exit face at the last
        # accepted step keeps that step's state instead.
        act[(bt == 2) & exit_on_edge & ~warm_active & ~prev_exit[0]] = True
        # Anti-cycling state for the all-or-nothing quadratic exit edge. A tri6/quad8
        # exit edge is either wholly in the active set or wholly out, so a PARTLY WET
        # edge -- exactly what the face carries while the drainage front's exit point
        # descends THROUGH that edge -- has no representation. Neither available state
        # is admissible there, so the two alternate: wet, the edge over-drains and the
        # rule sheds it; dry, pressure climbs above the face and the rule takes it back.
        # The set then runs a clean period-2 limit cycle, the Picard set-stability test
        # can never close, dt collapses to its floor and the march force-accepts steps
        # microseconds wide. ``seen_sets`` counts repeats so that cycle can be told from
        # ordinary settling, and ``pinned`` holds the resolution (below).
        seen_sets = {}
        shed = np.zeros(n_nodes, dtype=bool)    # exit nodes this step has shed
        pinned = np.zeros(n_nodes, dtype=bool)
        pinned_state = np.zeros(n_nodes, dtype=bool)
        cycling = False
        kold_h = None
        if theta < 1.0:
            kold_h = _coo_matvec(asm, kdata(h_old - y), h_old)
        p_old = h_old - y
        old_state = _storage_state(sasm, p_old)
        h_it = h_old.copy()
        pic = 0
        ok = False
        r_prev = None
        relax = None
        for pic in range(1, picard_max + 1):
            p_it = h_it - y
            kd = kdata(p_it)
            # The step's equations on a free row:
            #   R_i(h)/dt + theta (K h)_i + (1 - theta) (K_old h_old)_i = f_i,
            # R the storage term (_storage_residual). Each solve takes K at this
            # iterate's kr and linearizes R about this iterate on the diagonal D:
            #   (theta K + D/dt) h_new = (D h_it - R(h_it))/dt + ...
            R_it, D_it, _ = _storage_residual(sasm, p_it, p_old, old_state)
            dir_mask = (bt == 1) | ((bt == 2) & act)
            if not np.any(dir_mask):
                # singular (all exit nodes shed, no fixed head): cannot force.
                return False, h_old, warm_active, pic, None, cycling
            dir_values = np.where(bt == 1, bv, y)
            md = D_it / dt
            rhs = md * h_it - R_it / dt + theta * flux
            if theta < 1.0:
                rhs = rhs + (1.0 - theta) * flux - (1.0 - theta) * kold_h
            A, b = _transient_system(asm, theta * kd, md, rhs, dir_mask, dir_values)
            h_new = spsolve(A, b)
            # exit-face active-set update (SEEP2D per-corner rule on the reaction,
            # with quadratic exit edges kept all-or-nothing — ported from the
            # steady solve_unsaturated).
            if has_exit:
                f_eff = np.where(dir_mask, 0.0, flux)
                q = _coo_matvec(asm, kd, h_new) - f_eff
                hyst = 1e-3 * char_h
                # An apparent inflow is measured against the water OFFERED at the
                # node, max(flux, 0), not against zero: on a settled free row the
                # residual left after f_eff is round-off, so a strict sign test hands
                # a flux-loaded exit candidate to noise. Same threshold on both tests,
                # same reasoning as the steady switch (see solve_unsaturated); exactly
                # zero, hence bit-for-bit the strict rule, wherever no inflow is
                # prescribed.
                q_offered = np.maximum(flux, 0.0)
                is_c = (bt == 2) & exit_is_corner
                stay = is_c & act & ~((h_new < y - hyst) | (q > q_offered))
                turn_on = is_c & ~act & (h_new >= y + hyst) & (q <= q_offered)
                corner_candidate = np.zeros(n_nodes, dtype=bool)
                corner_candidate[stay | turn_on] = True

                new_act = np.zeros(n_nodes, dtype=bool)
                # linear exit edges: per-corner activation (tri3/quad4 rule)
                new_act[exit_linear_corners] = corner_candidate[exit_linear_corners]
                # quadratic exit edges: all-or-nothing on (corner, midside, corner)
                mid_cands = []
                _corner_ok = corner_candidate | (bt == 1)
                for _c1, _mid, _c2 in exit_quadratic_edges:
                    if act[_mid]:
                        mid_candidate = not (h_new[_mid] < y[_mid] - hyst
                                             or q[_mid] > q_offered[_mid])
                    else:
                        mid_candidate = (h_new[_mid] >= y[_mid] + hyst
                                         and q[_mid] <= q_offered[_mid])
                    mid_cands.append(bool(mid_candidate))
                    if _corner_ok[_c1] and mid_candidate and _corner_ok[_c2]:
                        for _nd in (_c1, _mid, _c2):
                            if bt[_nd] == 2:
                                new_act[_nd] = True
                # Rescue a singular full collapse of a quadratic exit face when NO
                # fixed head anchors the step (steady issue #51): fall back to
                # per-corner activation for any corner that passes its own h/q test.
                # Inert whenever a Dirichlet head exists this step (every reservoir /
                # earth-dam case keeps a specified-head row), so the linear path and
                # every anchored quadratic step are untouched.
                if not np.any(bt == 1):
                    face_active = np.any((bt == 2) & new_act)
                    stranded = is_c & corner_candidate & ~new_act
                    if not face_active and np.any(stranded):
                        new_act[stranded] = True
                # Anything already resolved by the anti-cycling rule below keeps the
                # state it was pinned to, so the cycle cannot restart.
                if pinned.any():
                    new_act[pinned] = pinned_state[pinned]
                set_stable = np.array_equal(new_act, act)
                # A seepage face with the head above the ground is not a solution
                # (same rule as the steady solver's; see solve_unsaturated). A set
                # about to be called stable with an inactive exit node under
                # pressure is resolved by hand: the node joins the face and is
                # pinned there for the rest of the step, the representation a partly
                # wet quadratic edge needs and the one the cycle breaker below hands
                # out when it catches the cycle.
                # The same rule for a node this step has already shed: one that comes
                # back standing at or above its own elevation joins the face and is
                # pinned at once, without waiting for a stable set. A step starts from
                # the face the last step ended with, and a partly wet quadratic edge
                # carried in wet can show a small inflow at its upper corner: shed for
                # it, the edge is under pressure; taken back, the inflow sheds it again,
                # and the set alternates without ever reaching the stable set the first
                # rule waits for. Measured on the exit-cycle guard's GW18 march (tri6,
                # target 3.0): 13 cycles detected, 4 damped, without this; 0 with it.
                if not _TSEEP_NO_WET_FACE_RULE:
                    _rejoin = (bt == 2) & shed & ~pinned & (h_new >= y + hyst)
                    if np.any(_rejoin):
                        new_act[_rejoin] = True
                        pinned[_rejoin] = True
                        pinned_state[_rejoin] = True
                        set_stable = np.array_equal(new_act, act)
                        n_wet[0] += int(np.count_nonzero(_rejoin))
                if set_stable and not _TSEEP_NO_WET_FACE_RULE:
                    _wet_free = (bt == 2) & ~new_act & (h_new >= y + hyst) & ~pinned
                    if np.any(_wet_free):
                        new_act[_wet_free] = True
                        pinned[_wet_free] = True
                        pinned_state[_wet_free] = True
                        set_stable = False
                        n_wet[0] += int(np.count_nonzero(_wet_free))

                # Break a limit cycle in the active set by giving the offending edge
                # the PARTIAL state the all-or-nothing rule cannot express. A set that
                # is unstable and has already been visited three times is cycling, not
                # settling (a stable set is excluded here outright). Every node on a
                # cycling edge is then resolved from its OWN SEEP2D h/q test -- the
                # per-node rule linear elements always use, which CAN represent a wet
                # lower half and a dry upper half -- and pinned there for the rest of
                # this step, so the seepage-face transition lands at the midside node
                # instead of oscillating across the whole edge.
                #
                # Detection is always on; the RESOLUTION only runs under ``damp``,
                # which the caller sets only when a step has already failed with the
                # set cycling. A step that converges the ordinary way is therefore
                # bit-for-bit what it was before this rule existed -- including the
                # steps that cycle for a sweep or two and then settle on their own.
                if not set_stable:
                    key = new_act.tobytes()
                    seen_sets[key] = seen_sets.get(key, 0) + 1
                    if seen_sets[key] >= 3:
                        if not cycling:
                            n_cycles[0] += 1
                        cycling = True
                    if damp and seen_sets[key] >= 3:
                        _resolve_exit_cycle(
                            act, new_act, corner_candidate, mid_cands,
                            exit_quadratic_edges, exit_linear_corners,
                            pinned, pinned_state)
                        new_act[pinned] = pinned_state[pinned]
                        set_stable = np.array_equal(new_act, act)
                        seen_sets.clear()
                        n_damped[0] += 1
                shed |= act & ~new_act
                act = new_act
            else:
                set_stable = True
            change = float(np.max(np.abs(h_new - h_it)))
            if change < picard_eps and set_stable and pic > 1:
                # The step keeps the iterate whose own full update is below the
                # tolerance, so it solves the step's equations to that tolerance;
                # the update itself need not, where the iteration alternates. On
                # lock 2c's tri6 leg, keeping the update left 48 of 816 steps
                # whose next update was up to seven times the tolerance. The first
                # iterate is the starting head, never kept: a step whose first
                # update is below the tolerance takes one more.
                ok = True
                break
            if set_stable:
                # Damp a 2-cycle. A node whose update reverses sign without shrinking
                # by half has its next update halved, again each time it persists,
                # and restored step by step once it stops. kr changes by orders of
                # magnitude across a sharp front, and an element straddling the front
                # can take a different kr at each iterate: on lock 2c's tri6 leg a
                # dry node (pressure head -16.6 ft; the band ends at -1 ft), its
                # storage the same at every iterate, saw its conductance swing
                # fourfold between iterates and its update alternate and grow. Measured on that leg (test/transient_seep_check.py,
                # tri6 at 22 divisions): 781 steps and 18 s with the damping, 2430
                # steps and 63 s without. Only the path to the solution changes.
                r = h_new - h_it
                if r_prev is not None:
                    flip = (r * r_prev < 0.0) & (np.abs(r) > 0.5 * np.abs(r_prev))
                    if np.any(flip):
                        if relax is None:
                            relax = np.ones(n_nodes)
                        relax[flip] *= 0.5
                    if relax is not None:
                        relax[~flip] = np.minimum(2.0 * relax[~flip], 1.0)
                r_prev = r
                h_it = h_it + r if relax is None else h_it + relax * r
            else:
                r_prev = None
                relax = None
                h_it = h_new

        if not ok and not force:
            return False, h_it, act, pic, None, cycling
        # head-change limiter (resolve sharp fronts even when Picard is happy).
        # Measured on FREE nodes only: a Dirichlet node's change is imposed by the
        # BC (a stepped reservoir head jumps by the full step regardless of dt), so
        # counting it would force dt -> dt_min forever with no benefit.
        free = ~dir_mask
        if (not force and np.any(free)
                and float(np.max(np.abs((h_it - h_old)[free]))) > dh_limit):
            return False, h_it, act, pic, None, cycling
        # The pressure this step would deliver on the free part of its exit face.
        # Every accepted step returns through here, force-accepted ones included.
        _free_exit = (bt == 2) & ~act if has_exit else None
        wet_probe[0] = (float(np.max(h_it[_free_exit] - y[_free_exit]))
                        if _free_exit is not None and np.any(_free_exit) else 0.0)
        # Ledger entry, from the accepted state itself (kr and storage both at the
        # accepted head): the boundary flow at each fixed row is that row's own
        # equation -- the water its storage took up plus what it passes into the
        # mesh -- and the specified flux enters at the free rows. The accepted state
        # solves the free rows to the iteration tolerance, so this inflow and the
        # change in stored water differ by exactly what it leaves unsolved there.
        R_new, _, _ = _storage_residual(sasm, h_it - y, p_old, old_state)
        rate = R_new / dt + theta * _coo_matvec(asm, kdata(h_it - y), h_it)
        if theta < 1.0:
            rate = rate + (1.0 - theta) * kold_h
        inflow = dt * float(np.sum(rate[dir_mask]) + np.sum(flux[~dir_mask]))
        return True, h_it, act, pic, (inflow, rate, flux, dir_mask, bt == 2), cycling

    t = 0.0
    dt = dt0
    cancelled = False
    prev_rate = None            # last accepted step's (h_new - h_old)/dt, free nodes
    prev_free = None
    dt_prev = None
    for t_target in targets:
        if cancelled:
            break
        guard = 0
        while t < t_target - 1e-12:
            guard += 1
            if guard > 500000:
                raise RuntimeError("transient stepper exceeded 5e5 sub-steps; "
                                   "check dt_min / convergence")
            dt_try = min(dt, dt_max, t_target - t)
            ok, h_new, new_active, pic, flow, cyc = attempt(h, dt_try, t + dt_try, active)
            # A step that failed WITH the exit-face set cycling is not a step that
            # needs a smaller dt: no dt closes a set-stability test the all-or-nothing
            # edge rule cannot satisfy, which is why this used to collapse to dt_min
            # and force-accept. Retry it once at the SAME dt with the cycle resolved,
            # and only fall back to halving if that still fails.
            if not ok and cyc:
                ok, h_new, new_active, pic, flow, cyc = attempt(
                    h, dt_try, t + dt_try, active, damp=True)
            while not ok and dt_try > dt_min:
                dt_try *= 0.5
                ok, h_new, new_active, pic, flow, cyc = attempt(
                    h, dt_try, t + dt_try, active)
                if not ok and cyc:
                    ok, h_new, new_active, pic, flow, cyc = attempt(
                        h, dt_try, t + dt_try, active, damp=True)
            if not ok:
                all_converged = False
                dt_try = max(dt_try, dt_min)
                if verbose:
                    print(f"  transient step failed to converge at t={t:.6g} "
                          f"(dt down to {dt_try:.3g}); force-accepting last iterate")
                # force-accept the smallest-dt iterate to make progress (flagged
                # non-converged). A singular system still can't be forced -> break.
                ok, h_new, new_active, pic, flow, cyc = attempt(
                    h, dt_try, t + dt_try, active, force=True, damp=True)
                if not ok:
                    raise SeepInputError(
                        f"Transient solve became singular at t={t:.6g} (the exit "
                        f"face fully deactivated with no specified head). Add a "
                        f"specified head or check the exit-face definition.")
                forced = True
            else:
                forced = False
            # Time error of the step just solved (none for the first step, which has
            # no rate before it).
            rate = (h_new - h) / dt_try
            free_now = ~flow[3]
            err = 0.0
            if prev_rate is not None:
                both = free_now & prev_free
                if np.any(both):
                    err = (dt_try * dt_try / (dt_try + dt_prev)
                           * float(np.max(np.abs(rate[both] - prev_rate[both]))))
            if (not forced and err > 2.0 * err_tol and dt_try > 1e3 * dt_min):
                dt = dt_try * max(0.25, 0.9 * np.sqrt(err_tol / err))
                continue
            prev_rate, prev_free, dt_prev = rate, free_now, dt_try
            t += dt_try
            h = h_new
            active = new_active
            max_free_p[0] = max(max_free_p[0], wet_probe[0])
            cum_inflow += flow[0]
            last_flow[0] = flow[1:4]
            prev_exit[0] = flow[4]
            dt_history.append(dt_try)
            # Progress + cancellation: report the simulated-time fraction and honor an
            # abort request (a falsy return) after this accepted step. Zero cost when
            # no callback is supplied.
            if progress_callback is not None and not progress_callback(t, duration):
                cancelled = True
                break
            # adapt dt for the next sub-step: by the iteration count, and no larger
            # than the time error allows
            if pic <= 5:
                dt = min(dt_try * growth, dt_max)
            elif pic > 0.8 * picard_max:
                dt = max(dt_try * 0.5, dt_min)
            else:
                dt = dt_try
            if err > 0.0:
                dt = max(min(dt, dt_try * 0.9 * np.sqrt(err_tol / err)), dt_min)
        if cancelled:
            break
        # landed exactly on t_target -> save frame with the boundary inflow/outflow
        # of the step that ended here (no dual flow-function solve).
        sc = _stored(h) - stored0
        closure = abs(sc - cum_inflow) / max(abs(cum_inflow), 1e-3 * stored_scale, 1e-30)
        in_s, out_s = _frame_flows(*last_flow[0])
        frames.append({"time": t_target, "head": h.copy(), "u": gamma_w * (h - y),
                       "inflow": in_s, "outflow": out_s})
        mb_frames.append({"time": t_target, "stored_change": sc,
                          "cumulative_inflow": cum_inflow, "closure": closure})
        if verbose:
            print(f"  t={t_target:.6g}: frame saved (steps so far={len(dt_history)}, "
                  f"mass-balance closure={closure:.2e})")

    # Final progress tick at 100% (skipped when the run was cancelled — the last
    # accepted-step tick already reported where the abort landed).
    if progress_callback is not None and not cancelled:
        progress_callback(duration, duration)

    if verbose and n_damped[0]:
        print(f"  exit face: {n_damped[0]} partly-wet quadratic edge(s) resolved "
              f"per-node (the all-or-nothing edge rule was cycling)")

    return {
        "times": [f["time"] for f in frames],
        "frames": frames,
        "dt_history": dt_history,
        "mass_balance": {
            "cumulative_inflow": cum_inflow,
            "final_stored_change": _stored(h) - stored0,
            "final_closure": mb_frames[-1]["closure"] if mb_frames else 0.0,
            "per_frame": mb_frames,
        },
        "converged": all_converged,
        "cancelled": cancelled,
        # Exit-face limit cycles: how many steps SAW the active set revisit a state
        # it had left, and how many times the partly-wet fallback had to be applied
        # to break one. The fallback is gated behind a step that has already failed,
        # so a march that converges normally reports damped = 0 however often the
        # set wobbles.
        "exit_face_cycles": n_cycles[0],
        "exit_face_damped": n_damped[0],
        # How many times the free-surface rule joined an inactive exit node standing
        # under pressure to the seepage face, and the worst pressure head left on the
        # free part of an exit face at the end of an accepted step. The second should
        # not exceed the hysteresis band: above it, the field carries pore pressure on
        # a surface that drains at atmospheric.
        "exit_face_wet_joined": n_wet[0],
        "exit_face_free_pressure": max_free_p[0],
        "exit_face_hyst": 1e-3 * char_h,
        "theta": theta,
        "lumped": lumped,
        "duration": duration,
        "unit_weight": gamma_w,
        # An exit face (static seepage face, or a series-head reservoir uncovering
        # a face it has submerged) can appear; kr scales velocity on load exactly
        # as in the steady unconfined branch. Recorded so the frame loader derives
        # velocities the same way.
        "unconfined": bool(_may_have_exit),
        # Rapid-drawdown stage times carried through for the meta ledger and the
        # in-memory drawdown staging (§6); None when the model sets no stages.
        "stage_1": tseep_data.get("stage_1"),
        "stage_2": tseep_data.get("stage_2"),
        # The single-run extraction time, carried the same way and for the same
        # reason: an exported march records which instant its stability run reads.
        "stability_time": tseep_data.get("stability_time"),
    }


def _transient_steady_ic(nodes, elements, element_types, bc_type, bc_values,
                         k1, k2, angle, kr0, h0, vg_a, vg_n, vg_l, model,
                         flux_nodal, verbose, max_iter=2000):
    """Initial condition = a steady solve at the t=0 BC configuration (plan IC rule).
    Confined (no exit face) -> one linear solve; unconfined -> the steady
    exit-face Picard solver. Reuses the existing steady solvers unchanged.

    Returns ``(head, converged)``. The flag is the unconfined solver's own, and it
    belongs to the caller: an initial condition that did not close is a march
    starting from a field that is not a steady state, which is a fact about the
    solution and not a line of console output. A confined IC is one direct linear
    solve and always converged.

    ``max_iter`` is the sweep budget for the unconfined solve. It is the caller's to
    set because the steady runners already let a model that needs more sweeps ask
    for them, and the initial condition is the same solve on the same models."""
    if np.any(bc_type == 2):
        head, _A, _q, _tf, _efa, conv, _ce = solve_unsaturated(
            nodes, elements, bc_type, bc_values, kr0=kr0, h0=h0,
            k1_vals=k1, k2_vals=k2, angles=angle, element_types=element_types,
            vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model, flux_nodal=flux_nodal,
            max_iter=max_iter)
        return head, bool(conv)
    bcs = [(i, bc_values[i]) for i in range(len(bc_type)) if bc_type[i] == 1]
    head, _A, _q, _tf = solve_confined(nodes, elements, bc_type, bcs, k1, k2, angle,
                                       element_types, flux_nodal=flux_nodal)
    return head, True


def export_seep_solution(seep_data, solution, filename):
    """Exports nodal results to a CSV file.
    
    The exported CSV file contains the following columns:
    - node_id: Node identifier (1-based)
    - head: Hydraulic head at each node
    - u: Pore pressure at each node
    - v_x, v_y: Velocity vector components
    - v_mag: Velocity magnitude
    - i_x, i_y: Hydraulic gradient vector components
    - i_mag: Hydraulic gradient magnitude
    - q: Nodal flow vector
    - phi: Stream function/flow potential

    Below the table the file records what the solve was, as ``#`` comment lines the
    CSV readers skip: the total flowrate, and — for a solution that carries them —
    whether the problem was solved unconfined, whether the solve converged, and the
    flow imbalance it closed to, written twice: ``Closure Error`` is the absolute
    imbalance, a flow rate in the model's own units, and ``Closure Fraction`` is that
    imbalance over the flow through the section. Those are what
    :func:`run_seepage_analysis` returns beside the nodal fields and the only part of
    its answer the columns cannot hold; without them a solution read back could not
    say whether its negative pore pressures mean a phreatic surface (unconfined) or
    are the ordinary suction of a saturated potential field (confined). Written only
    when the solution carries them, so a re-exported imported field gains no footer
    it has no answer for.

    Args:
        filename: Path to the output CSV file
        seep_data: Dictionary containing seep data
        solution: Dictionary containing solution results from run_seepage_analysis
    """
    import pandas as pd
    n_nodes = len(seep_data["nodes"])
    df = pd.DataFrame({
        "node_id": np.arange(1, n_nodes + 1),  # Generate 1-based node IDs for output
        "head": solution["head"],
        "u": solution["u"],
        "v_x": solution["velocity"][:, 0],
        "v_y": solution["velocity"][:, 1],
        "v_mag": solution["v_mag"],
        "i_x": solution["gradient"][:, 0],
        "i_y": solution["gradient"][:, 1],
        "i_mag": solution["i_mag"],
        "q": solution["q"],
        "phi": solution["phi"]
    })
    # Write a units-provenance header (only when the model declares a system; blank
    # otherwise, so undeclared exports stay byte-identical), then the table, then the
    # solve footer. Every # line is an ordinary CSV comment the readers skip.
    from .units import units_comment_line
    header = units_comment_line(seep_data.get("unit_system"), seep_data.get("time_unit"))
    with open(filename, "w") as f:
        if header:
            f.write(header)
        df.to_csv(f, index=False)
        f.write(f"# Total Flowrate: {solution['flowrate']:.6f}\n")
        # The three solve facts the columns cannot hold. Appended AFTER the flowrate
        # line so a file written today is a file written before it plus these lines:
        # every reader that skips "#" comments is unaffected, and import_seep_solution
        # reads them where they are and falls back where they are not.
        if solution.get("unconfined") is not None:
            f.write(f"# Unconfined: {bool(solution['unconfined'])}\n")
        if solution.get("converged") is not None:
            f.write(f"# Converged: {bool(solution['converged'])}\n")
        if solution.get("closure_error") is not None:
            f.write(f"# Closure Error: {float(solution['closure_error']):.6e}\n")
            # The same imbalance as a fraction of the flow through the section. The
            # line above it is an absolute flow rate in the model's own units, which
            # a bare number in a footer reads as a proportion when it is nothing of
            # the kind — "0.0266" against a flow of 1.27 is 2% of it, not 2.66%.
            # Written beside it rather than instead of it, so nothing that reads the
            # absolute value changes meaning, and appended AFTER it so a file written
            # today is still a file written before it plus lines.
            if solution.get("closure_fraction") is not None:
                f.write(f"# Closure Fraction: "
                        f"{float(solution['closure_fraction']):.6e}\n")

    print(f"Exported solution to {filename}")


def export_seep_u(nodes, u, filename, gamma_water):
    """Write a bare nodal pore-pressure field in the same CSV format, for a solution
    xslope did not compute.

    ``gamma_water`` is required (no default): it sets the head implied by each
    pressure, and a hardcoded fallback here could silently write heads in the wrong
    unit system. Callers pass the model's declared value (see
    :func:`xslope.units.require_gamma_water`).

    A pore-pressure field can arrive from outside — GeoStudio's importer lifts one
    straight out of a solved SEEP/W analysis. It carries pressures and nothing else: no
    velocities, no gradients, no stream function. Those columns are therefore OMITTED
    rather than filled with zeros, which would render as a flow net that is flat
    everywhere and belongs to nobody. ``load_slope_data`` reads only ``u``, so the file
    is complete for stability analysis; a flow-net plot will fail loudly on the missing
    columns instead of drawing a lie.

    ``head`` is the total head implied by the pressure: elevation + u / gamma_w.
    """
    import pandas as pd
    df = pd.DataFrame({
        "node_id": np.arange(1, len(nodes) + 1),
        "head": np.asarray(nodes)[:, 1] + np.asarray(u, dtype=float) / gamma_water,
        "u": np.asarray(u, dtype=float),
    })
    with open(filename, "w") as f:
        df.to_csv(f, index=False)
        # A total-flowrate footer, matching export_seep_solution's format so the same
        # readers (which skip "#" comment lines) treat both files identically. There is
        # no solved flowrate here -- the field was imported, not computed.
        f.write("# Total Flowrate: not computed (field imported)\n")


def import_seep_solution(seep_data, filename):
    """Reconstruct a seepage ``solution`` dict from a CSV written by
    :func:`export_seep_solution` — the inverse operation. Lets a previously
    saved solution be re-plotted (``plot_seep_solution``) without re-running the
    analysis.

    Args:
        seep_data: Seep data for the mesh the solution was computed on (used to
            validate the node count).
        filename: Path to a ``*_seep.csv`` / ``*_seep2.csv`` file.

    Returns:
        dict: solution with the keys ``plot_seep_solution`` expects — head, u,
        velocity, v_mag, gradient, i_mag, q, phi, flowrate — plus ``unconfined``,
        ``converged``, ``closure_error`` and ``closure_fraction`` for a file whose
        footer records them (``closure_fraction`` is the newest of the four, so a
        companion written before it reads back with the other three and without it).
        Those keys are ABSENT for a file written before the footer existed, or
        by :func:`export_seep_u`, which has no answer for them: absent means unknown,
        and a consumer falls back to what it can see (``plot_seep_solution`` treats an
        unknown solve as unconfined, its long-standing default). For a pressure-only
        file (one written by :func:`export_seep_u`, e.g. a field imported from
        SEEP/W), the flow-net columns are filled with NaN: head and u are real, and a
        flow-net plot fails visibly on the NaNs rather than drawing a flat field that
        is not the solution.

    Raises:
        ValueError: if the file's node count does not match the mesh.
    """
    import pandas as pd
    # The trailing '# …' footer lines are comments, skipped on read.
    df = pd.read_csv(filename, comment="#")
    n_nodes = len(seep_data["nodes"])
    if len(df) != n_nodes:
        raise ValueError(
            f"Seepage solution has {len(df)} nodes but the mesh has {n_nodes} — "
            "the saved solution does not match this mesh.")

    flowrate = None
    # The solve facts the footer may record. Left absent from the returned dict where
    # the file does not record them — every solution written before these lines
    # existed reads back exactly as it did, and "unknown" stays distinguishable from
    # "recorded as False".
    meta = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("# Total Flowrate:"):
                try:
                    flowrate = float(line.split(":", 1)[1])
                except ValueError:
                    pass
            elif line.startswith("# Unconfined:"):
                flag = line.split(":", 1)[1].strip().lower()
                if flag in ("true", "false"):
                    meta["unconfined"] = flag == "true"
            elif line.startswith("# Converged:"):
                flag = line.split(":", 1)[1].strip().lower()
                if flag in ("true", "false"):
                    meta["converged"] = flag == "true"
            elif line.startswith("# Closure Error:"):
                try:
                    meta["closure_error"] = float(line.split(":", 1)[1])
                except ValueError:
                    pass
            elif line.startswith("# Closure Fraction:"):
                try:
                    meta["closure_fraction"] = float(line.split(":", 1)[1])
                except ValueError:
                    pass

    # A field imported from outside (export_seep_u) carries head and u only. Fill the
    # flow-net quantities with NaN rather than inventing zeros -- a plot then fails on
    # the NaN instead of drawing a flat, wrong flow net, while stability, which uses only
    # u, is unaffected.
    nan = np.full(n_nodes, np.nan)
    def col(name):
        return df[name].to_numpy() if name in df.columns else nan.copy()

    return {
        "head": df["head"].to_numpy(),
        "u": df["u"].to_numpy(),
        "velocity": np.column_stack([col("v_x"), col("v_y")]),
        "v_mag": col("v_mag"),
        "gradient": np.column_stack([col("i_x"), col("i_y")]),
        "i_mag": col("i_mag"),
        "q": col("q"),
        "phi": col("phi"),
        "flowrate": flowrate,
        **meta,
    }


# =============================================================================
# Transient output files + frame loader + drawdown staging (plan §4 / §6)
# =============================================================================


def _transient_base(path):
    """Normalize a caller-supplied path to the transient file base (no suffix).
    Accepts a bare base, ``{base}_tseep.csv`` or ``{base}_tseep_meta.json``."""
    for suf in ("_tseep.csv", "_tseep_meta.json"):
        if path.endswith(suf):
            return path[: -len(suf)]
    if path.endswith(".csv"):
        return path[:-4]
    return path


def _elem_props(seep_data):
    """Per-element (k1, k2, angle, kr0, h0, vg_a, vg_n, vg_l, model) arrays from the
    per-material tables in ``seep_data`` — the inputs compute_velocity /
    solve_flow_function_* expect. Shared by the solve-time flow net and the
    derive-on-load frame reconstruction."""
    elements = np.asarray(seep_data["elements"])
    element_types = seep_data.get("element_types")
    if element_types is None:
        element_types = np.full(len(elements), 3)
    element_types = np.asarray(element_types)
    mat_ids = np.asarray(seep_data["element_materials"]) - 1
    k1 = seep_data["k1_by_mat"][mat_ids]
    k2 = seep_data["k2_by_mat"][mat_ids]
    angle = seep_data["angle_by_mat"][mat_ids]
    kr0 = seep_data["kr0_by_mat"][mat_ids]
    h0 = seep_data["h0_by_mat"][mat_ids]
    _u = seep_data.get("unsat_by_mat")
    _a = seep_data.get("vg_a_by_mat")
    _n = seep_data.get("vg_n_by_mat")
    _l = seep_data.get("vg_l_by_mat")
    model = None if _u is None else _u[mat_ids]
    vg_a = None if _a is None else _a[mat_ids]
    vg_n = None if _n is None else _n[mat_ids]
    vg_l = None if _l is None else _l[mat_ids]
    return elements, element_types, k1, k2, angle, kr0, h0, vg_a, vg_n, vg_l, model


def _transient_frame_solution(seep_data, head, u, phi, inflow, outflow,
                              unconfined, time=None):
    """Reconstruct the ``solution`` dict ``plot_seep_solution`` consumes from one
    stored transient frame.

    head and u are STORED; velocity and gradient are DERIVED here from head + mesh
    + k(kr) (the plan's derive-on-load contract). This is valid only while kr is a
    pure function of psi; a future hysteresis / storage-dependent-kr model would
    have to store velocity columns instead. No stream function is stored — a
    transient (storage-release) state is not divergence-free, so it has no flow net
    (flow lines are read from velocity vectors instead); ``phi`` is therefore
    ``None`` for a current run and is accepted only so a legacy 5-column file still
    loads. ``flowrate`` mirrors the boundary inflow; ``inflow``/``outflow`` (which
    differ under transient storage exchange) and ``time`` drive the title variant."""
    nodes = seep_data["nodes"]
    (elements, element_types, k1, k2, angle, kr0, h0,
     vg_a, vg_n, vg_l, model) = _elem_props(seep_data)
    if unconfined:
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle,
                                    kr0=kr0, h0=h0, element_types=element_types,
                                    vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
    else:
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle,
                                    element_types=element_types)
    gradient = compute_gradient(nodes, elements, head, element_types)
    return {
        "head": np.asarray(head, dtype=float),
        "u": np.asarray(u, dtype=float),
        "velocity": velocity,
        "v_mag": np.linalg.norm(velocity, axis=1),
        "gradient": gradient,
        "i_mag": np.linalg.norm(gradient, axis=1),
        "phi": None if phi is None else np.asarray(phi, dtype=float),
        # Mirrors the boundary inflow for API parity with a steady solution; no flow
        # lines are drawn for a transient frame (no stream function), so it is not a
        # channel count. inflow/outflow carry the honest (unequal) boundary rates.
        "flowrate": inflow,
        "inflow": inflow,
        "outflow": outflow,
        "unconfined": bool(unconfined),
        "time": time,
    }


def export_transient_solution(seep_data, transient_solution, base_path,
                              input_file=None, mesh_file=None):
    """Write a transient run to ``{base}_tseep.csv`` (+ ``{base}_tseep_meta.json``).

    The CSV is long format — ``time, node_id, head, u`` — one block of n_nodes
    rows per saved frame (plan §4). Velocities are NOT stored: they are derived on
    load from head + mesh + k(kr) (:func:`_transient_frame_solution`). No stream
    function is stored either — a transient (storage-release) state has no flow net
    (see the note in ``docs/seep/transient.md``), so a per-frame flow-function
    solve would be both wrong to draw and a needless cost at export.

    The meta JSON carries the frame ledger, per-frame inflow/outflow, dt history,
    the mass-balance ledger, series/stage times and provenance (input/mesh files,
    unit system + time unit). ``base_path`` may be a bare base or any of the two
    file names; the suffixes are normalized off.
    """
    import json
    import pandas as pd
    base = _transient_base(base_path)
    nodes = seep_data["nodes"]
    n_nodes = len(nodes)
    frames = transient_solution["frames"]

    blocks = []
    for fr in frames:
        head = np.asarray(fr["head"], dtype=float)
        u = np.asarray(fr.get("u"), dtype=float)
        blocks.append(pd.DataFrame({
            "time": np.full(n_nodes, fr["time"]),
            "node_id": np.arange(1, n_nodes + 1),
            "head": head,
            "u": u,
        }))
    df = pd.concat(blocks, ignore_index=True)

    csv_path = f"{base}_tseep.csv"
    from .units import units_comment_line
    header = units_comment_line(seep_data.get("unit_system"),
                                seep_data.get("time_unit"))
    with open(csv_path, "w") as f:
        if header:
            f.write(header)
        df.to_csv(f, index=False)

    meta = {
        "format": "xslope_tseep_v1",
        "times": [float(fr["time"]) for fr in frames],
        "n_nodes": int(n_nodes),
        "unconfined": bool(transient_solution.get("unconfined", False)),
        "duration": transient_solution.get("duration"),
        "theta": transient_solution.get("theta"),
        "lumped": transient_solution.get("lumped"),
        "converged": bool(transient_solution.get("converged", True)),
        "unit_system": seep_data.get("unit_system"),
        "time_unit": seep_data.get("time_unit"),
        "unit_weight": transient_solution.get("unit_weight"),
        "dt_history": [float(x) for x in transient_solution.get("dt_history", [])],
        "mass_balance": transient_solution.get("mass_balance"),
        "frame_flows": [
            {"time": float(fr["time"]),
             "inflow": fr.get("inflow"),
             "outflow": fr.get("outflow")}
            for fr in frames
        ],
        "provenance": {
            "input_file": input_file,
            "mesh_file": mesh_file,
        },
    }
    # Stage and stability times live on the tseep controls; pass them through.
    for key in ("stage_1", "stage_2", "stability_time"):
        if key in transient_solution:
            meta[key] = transient_solution[key]
    def _json_default(o):
        # Convert numpy scalars/arrays faithfully; raise (never silently null) on
        # anything truly unserializable, so a leaked type surfaces instead of
        # writing a null over a physical value (the u-option lesson).
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(f"transient meta: value of type {type(o).__name__} is not "
                        f"JSON serializable")
    meta_path = f"{base}_tseep_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2, default=_json_default)

    print(f"Exported {len(frames)} transient frame(s) to {csv_path} and {meta_path}")
    return csv_path, meta_path


def import_transient_solution(seep_data, base_path):
    """Reconstruct a transient solution from ``{base}_tseep.csv`` (+ meta).

    Returns a dict shaped like :func:`run_transient_seepage`'s output that the
    plotter and drawdown staging consume: ``times`` (list), ``frames`` (list of
    the solution dicts :func:`plot_seep_solution` takes — velocities/gradients
    derived on load, head/u read from the file), and ``meta``. The lean current
    format is ``time, node_id, head, u``; a legacy 5-column file that still carries
    a ``phi`` column loads fine — the column is ignored (a transient state has no
    flow net). Raises if the file's node count does not match the mesh.
    """
    import json
    import os
    import pandas as pd
    base = _transient_base(base_path)
    csv_path = f"{base}_tseep.csv"
    meta_path = f"{base}_tseep_meta.json"
    meta = {}
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
    n_nodes = len(seep_data["nodes"])
    df = pd.read_csv(csv_path, comment="#")
    # Per-frame inflow/outflow from the meta (keyed by time), else None.
    flows = {round(float(ff["time"]), 12): ff
             for ff in meta.get("frame_flows", [])}
    unconfined = bool(meta.get("unconfined", False))

    # A legacy 5-column file still carries a phi column; tolerate it by ignoring
    # that column (a transient state has no flow net, so phi was dropped from the
    # format). Noted once at debug level — days-old files are graceful, not an error.
    if "phi" in df.columns:
        logger.debug("Transient frames file %s carries a legacy 'phi' column; "
                     "ignoring it (a transient state has no flow net).", csv_path)

    frames = []
    times = []
    for t, g in df.groupby("time", sort=True):
        g = g.sort_values("node_id")
        if len(g) != n_nodes:
            raise ValueError(
                f"Transient frame t={t} has {len(g)} nodes but the mesh has "
                f"{n_nodes} — the saved solution does not match this mesh.")
        ff = flows.get(round(float(t), 12), {})
        frames.append(_transient_frame_solution(
            seep_data, g["head"].to_numpy(), g["u"].to_numpy(), None,
            ff.get("inflow"), ff.get("outflow"), unconfined, time=float(t)))
        times.append(float(t))
    return {"times": times, "frames": frames, "meta": meta,
            "unconfined": unconfined}


def transient_frame_index(transient_solution, t, tol=None):
    """Index of the saved frame nearest time ``t``. Saved frames are computed
    states clamped to the schedule, so a requested stage/analysis time should
    match one exactly; ``tol`` (default: a small fraction of the total duration)
    guards against silently pulling a far-away frame."""
    times = np.asarray(transient_solution["times"], dtype=float)
    if len(times) == 0:
        raise ValueError("transient solution carries no saved frames")
    i = int(np.argmin(np.abs(times - t)))
    if tol is None:
        span = float(times[-1] - times[0]) or 1.0
        tol = 1e-6 * span + 1e-9
    if abs(times[i] - t) > tol:
        raise ValueError(
            f"No saved transient frame at time {t}; nearest is {times[i]} "
            f"(available: {', '.join(f'{x:g}' for x in times)}). Add it to "
            f"save_times / stage_1 / stage_2 / stability_time so the stepper lands "
            f"on it, or re-run the analysis with remarch_for_times().")
    return i


def has_transient_frame(transient_solution, t, tol=None):
    """``True`` when ``t`` names a saved frame. The question the free-entry path asks
    before deciding whether a re-march is needed."""
    try:
        transient_frame_index(transient_solution, t, tol=tol)
    except ValueError:
        return False
    return True


def resolve_stability_time(slope_data, transient_solution, time=None):
    """Which instant a stability run reads, and where that answer came from.

    Returns ``(time, source)`` with ``source`` one of:

    ``"argument"``
        a time was passed in — a Run-dialog choice or an explicit ``time=``;
    ``"file"``
        the model's ``tseep`` ``stability_time``, the reproducibility default that
        makes a headless re-run repeat the choice;
    ``"default"``
        neither, so the engine default stands: the LAST saved frame, usually the
        drained end state. That is what a blank ``stability_time`` means, and callers
        should say so rather than let it pass silently.

    Precedence is argument > file > engine default, per the template spec.
    """
    if time is not None:
        return float(time), "argument"
    ts = slope_data.get("tseep") or {}
    st = ts.get("stability_time")
    if st is not None:
        return float(st), "file"
    return float(transient_solution["times"][-1]), "default"


def _frame_u(transient_solution, i):
    """Pore-pressure array of frame ``i`` (works for both the raw run output —
    frames carry ``u`` — and a loaded solution)."""
    return np.asarray(transient_solution["frames"][i]["u"], dtype=float)


def select_transient_frame_u(slope_data, transient_solution, time=None):
    """Single-time LEM/FEM (plan §6): place the pore pressures of one frame into
    ``slope_data['seep_u']``. The downstream ``u = seep`` machinery is untouched.

    With no ``time`` the instant is resolved by :func:`resolve_stability_time` —
    the model's ``tseep`` ``stability_time`` when it declares one, otherwise the
    last saved frame. That is what makes a scripted re-run reproduce a Studio run's
    choice without repeating it at the call site.
    """
    time, _source = resolve_stability_time(slope_data, transient_solution, time)
    i = transient_frame_index(transient_solution, time)
    slope_data["seep_u"] = _frame_u(transient_solution, i)
    # Record WHICH moment these pore pressures belong to. Under automatic water
    # loads the surface load is derived from the pool at the moment being solved,
    # and this is how the derivation learns it: without it a transient run would
    # carry the field from time t under the reservoir as it stood at t = 0.
    slope_data["seep_u_time"] = float(transient_solution["times"][i])
    return slope_data


def stage_transient_for_drawdown(slope_data, transient_solution):
    """In-memory rapid-drawdown staging (plan §6, no file export).

    Pulls the two frames at the tseep ``stage_1``/``stage_2`` times and builds
    ``slope_data['seep_u']`` / ``['seep_u2']`` — exactly the structures the
    classic two-steady-file path produces — so the existing 3-stage LEM machinery
    runs unchanged. Resolution order (per the plan): a transient solution + stage
    times present takes precedence over the classic ``{base}_seep.csv`` /
    ``_seep2.csv`` files; without stage times the classic path is used.
    """
    ts = slope_data.get("tseep") or {}
    s1, s2 = ts.get("stage_1"), ts.get("stage_2")
    if s1 is None and s2 is None:
        raise ValueError(
            "stage_transient_for_drawdown needs stage_1 and stage_2 on the tseep "
            "controls; neither is set. Set both (stage_1 < stage_2) to couple the "
            "transient solution to rapid drawdown, or use the classic two-file path.")
    if s1 is None or s2 is None:
        raise ValueError("Rapid-drawdown staging needs BOTH stage_1 and stage_2 "
                         "(one is blank).")
    if s1 >= s2:
        raise ValueError(f"stage_1 ({s1}) must be < stage_2 ({s2}).")
    i1 = transient_frame_index(transient_solution, s1)
    i2 = transient_frame_index(transient_solution, s2)
    slope_data["seep_u"] = _frame_u(transient_solution, i1)
    slope_data["seep_u2"] = _frame_u(transient_solution, i2)
    # The moments the two staged fields belong to, for the same reason as in
    # select_transient_frame_u: an automatic water load is derived from the pool as
    # it stood at the instant the field was read, not at t = 0. Stage 2's instant
    # does a second job -- it is the only mark that says this drawdown is running
    # the transient-staged route, and the water derivation reads it to take the
    # drawn-down pool from boundary set 1 at that instant instead of from boundary
    # set 2, which belongs to the other route (xslope.water.stage_water_source).
    slope_data["seep_u_time"] = float(transient_solution["times"][i1])
    slope_data["seep_u2_time"] = float(transient_solution["times"][i2])
    return slope_data


def remarch_for_times(seep_data, slope_data, times, progress_callback=None, **march_kw):
    """Re-run the transient march with ``times`` injected into ``save_times``, and
    return the new solution.

    This is how a time that is not already a saved frame is served — a free-entry
    stability time, or stage times the existing solution never saved. **A field is
    never interpolated between frames:** a blend of two nodal fields is not a
    solution of anything, so the requested instant is made into a computed state
    instead. The stepper already clamps to the saved schedule, so every injected
    time comes back as a frame the solver actually landed on.

    The cost is a full re-solve, which is not uniformly cheap — a long march can take
    minutes. Callers that can show progress should pass ``progress_callback`` (the
    same callable :func:`run_transient_seepage` takes, which doubles as the
    cancellation channel).

    Times at or below zero, or beyond the run duration, are dropped: the march cannot
    land on them.
    """
    tseep_data = build_tseep_data(slope_data)
    if tseep_data is None:
        raise SeepInputError("remarch_for_times: this model has no 'tseep' sheet, so "
                             "there is no transient seepage analysis to re-run.")
    duration = tseep_data.get("duration")
    if duration is None or duration <= 0:
        raise SeepInputError("remarch_for_times: the tseep controls set no positive "
                             "duration.")
    wanted = sorted({float(t) for t in times if t is not None
                     and 0 < float(t) <= float(duration)})
    tseep_data = dict(tseep_data)
    tseep_data["save_times"] = sorted(
        set(tseep_data.get("save_times") or []) | set(wanted))
    return run_transient_seepage(seep_data, tseep_data,
                                 progress_callback=progress_callback, **march_kw)


def apply_steady_stability_field(slope_data, solution, bc=1, verbose=True):
    """Put a solved STEADY seepage field into the model as its pore-pressure field.

    The steady counterpart of :func:`apply_transient_stability_frame`. A stability run
    with ``u = seep`` reads ``slope_data['seep_u']`` (``['seep_u2']`` for the second
    BC set, which supplies rapid drawdown's stage 2). Three things can put a field
    there and all three are the same field: :func:`load_slope_data` reading a saved
    ``{base}_seep.csv``, the transient stagers reading one instant of a march, and
    this — a solution that has just been computed and not yet been written anywhere.

    A solve that filled only a results tab left the model itself carrying no field:
    a run gate then reported, correctly for what it could see, that the model has no
    seepage solution, and a run that got past a gate would have read pore pressures
    from the PREVIOUS solve's saved file. Neither is a display problem. A solved
    field belongs to the model the moment it exists, not the moment it is reloaded
    from disk.

    Parameters
    ----------
    solution : dict
        A solution as returned by :func:`run_seepage_analysis` (or reconstructed by
        :func:`import_seep_solution`) — only its ``u`` array is read.
    bc : int
        Which BC set this solution solved: 1 -> ``seep_u``, 2 -> ``seep_u2``.

    Raises
    ------
    ValueError
        If the field carries a different number of node values than the model's
        mesh has nodes — it was computed on a different mesh, and silently keeping
        it would price every slice base against the wrong geometry.
    """
    if solution is None:
        raise ValueError("apply_steady_stability_field needs a solved seepage "
                         "solution; run or import the analysis first.")
    if bc not in (1, 2):
        raise ValueError(f"seep BC set must be 1 or 2, got {bc!r}.")
    u = np.asarray(solution["u"], dtype=float)
    mesh = slope_data.get("mesh")
    nodes = None if mesh is None else mesh.get("nodes")
    if nodes is not None and len(u) != len(nodes):
        raise ValueError(
            f"The seepage solution carries {len(u)} node value(s) but the model's "
            f"mesh has {len(nodes)} node(s), so it was computed on a different mesh. "
            f"Re-run the seepage analysis on this mesh.")
    key = "seep_u" if bc == 1 else "seep_u2"
    slope_data[key] = u
    # A steady field belongs to no instant. Leaving a stale time behind would date
    # these pressures to a moment of a march they did not come from -- and under
    # automatic water loads that moment is what the surface load is derived from, so
    # the pool would be read at the wrong time. For set 2 it would do worse: the
    # stage-2 instant is also the mark of the transient-staged route, so a stale one
    # would send the drawn-down water load back to boundary set 1 on a run that has
    # just been given set 2's own solved field.
    slope_data.pop("seep_u_time" if bc == 1 else "seep_u2_time", None)
    if verbose:
        print(f"Analysis uses the steady seepage solution (BC set {bc}), "
              f"{len(u)} nodal pore pressures.")
    return slope_data


def apply_transient_stability_frame(slope_data, transient_solution, time=None,
                                    rapid=False, seep_data=None, remarch=False,
                                    progress_callback=None, verbose=True):
    """Point a stability run at one instant of a transient seepage solution.

    This is the entry surface for an LEM or FEM analysis with ``u = seep`` against a
    transient march: it resolves WHICH instant, puts that frame's pore pressures into
    ``slope_data['seep_u']`` (and ``['seep_u2']`` for rapid drawdown), and reports
    what it did.

    Parameters
    ----------
    time : float, optional
        The instant to read. Omitted, the model's ``tseep`` ``stability_time`` is
        used; with that blank too, the LAST saved frame — usually the drained end
        state, which is often but not always the critical one. Ignored for a rapid
        drawdown, whose two instants are the ``stage_1`` / ``stage_2`` times.
    rapid : bool
        Stage the run for rapid drawdown (``stage_1`` and ``stage_2`` frames into
        ``seep_u`` / ``seep_u2``) instead of reading a single instant.
    seep_data, remarch : optional
        With ``remarch=True`` and ``seep_data`` supplied, a requested time with no
        saved frame is served by re-marching with that time injected into
        ``save_times`` (:func:`remarch_for_times`) rather than raising. The march is
        never interpolated.
    verbose : bool
        Print the one-line statement of which instant was used and where the choice
        came from. A run against a transient field should never be silent about the
        instant it read.

    Returns
    -------
    dict
        ``{"times": [...], "source": ..., "remarched": bool, "solution": ...}`` —
        ``solution`` is the (possibly re-marched) solution the frames came from, so a
        caller can keep the re-marched one instead of paying for it twice.
    """
    if transient_solution is None:
        raise ValueError("apply_transient_stability_frame needs a transient solution; "
                         "run or import the transient seepage analysis first.")
    remarched = False
    if rapid:
        ts = slope_data.get("tseep") or {}
        wanted = [ts.get("stage_1"), ts.get("stage_2")]
        if (remarch and seep_data is not None
                and all(t is not None for t in wanted)
                and not all(has_transient_frame(transient_solution, t) for t in wanted)):
            transient_solution = remarch_for_times(
                seep_data, slope_data, wanted, progress_callback=progress_callback)
            remarched = True
        stage_transient_for_drawdown(slope_data, transient_solution)
        used = [float(t) for t in wanted]
        source = "stages"
        if verbose:
            print(f"Rapid drawdown uses transient stages {used[0]:g} / {used[1]:g}"
                  + (" (re-run)." if remarched else "."))
    else:
        t, source = resolve_stability_time(slope_data, transient_solution, time)
        if remarch and seep_data is not None and not has_transient_frame(
                transient_solution, t):
            transient_solution = remarch_for_times(
                seep_data, slope_data, [t], progress_callback=progress_callback)
            remarched = True
        select_transient_frame_u(slope_data, transient_solution, time=t)
        used = [float(slope_data["seep_u_time"])]
        if verbose:
            where = {"argument": "as requested",
                     "file": "from the file's stability_time",
                     "default": "the last saved frame (no stability_time set)"}[source]
            print(f"Analysis uses transient seepage frame t = {used[0]:g} — {where}"
                  + (", re-run." if remarched else "."))
    return {"times": used, "source": source, "remarched": remarched,
            "solution": transient_solution}


def print_seep_data_diagnostics(seep_data):
    """
    Diagnostic function to print out the contents of seep_data after loading.
    
    Args:
        seep_data: Dictionary containing seep data
    """
    print("\n" + "="*60)
    print("SEEP DATA DIAGNOSTICS")
    print("="*60)
    
    # Basic problem information
    print(f"Number of nodes: {len(seep_data['nodes'])}")
    print(f"Number of elements: {len(seep_data['elements'])}")
    print(f"Number of materials: {len(seep_data['k1_by_mat'])}")
    print(f"Unit weight of water: {seep_data['unit_weight']}")
    
    # Element type information
    element_types = seep_data.get('element_types', None)
    if element_types is not None:
        num_triangles = np.sum(element_types == 3)
        num_quads = np.sum(element_types == 4)
        print(f"Element types: {num_triangles} triangles, {num_quads} quadrilaterals")
    else:
        print("Element types: All triangles (legacy format)")
    
    # Coordinate ranges
    coords = seep_data['nodes']
    print(f"\nCoordinate ranges:")
    print(f"  X: {coords[:, 0].min():.3f} to {coords[:, 0].max():.3f}")
    print(f"  Y: {coords[:, 1].min():.3f} to {coords[:, 1].max():.3f}")
    
    # Boundary conditions
    bc_type = seep_data['bc_type']
    bc_values = seep_data['bc_values']
    print(f"\nBoundary conditions:")
    print(f"  Fixed head nodes (bc_type=1): {np.sum(bc_type == 1)}")
    print(f"  Exit face nodes (bc_type=2): {np.sum(bc_type == 2)}")
    print(f"  Free nodes (bc_type=0): {np.sum(bc_type == 0)}")
    
    if np.sum(bc_type == 1) > 0:
        fixed_head_nodes = np.where(bc_type == 1)[0]
        print(f"  Fixed head values: {bc_values[fixed_head_nodes]}")
    
    if np.sum(bc_type == 2) > 0:
        exit_face_nodes = np.where(bc_type == 2)[0]
        print(f"  Exit face elevations: {bc_values[exit_face_nodes]}")
    
    # Material properties
    print(f"\nMaterial properties:")
    for i in range(len(seep_data['k1_by_mat'])):
        print(f"  Material {i+1}:")
        print(f"    k1 (major conductivity): {seep_data['k1_by_mat'][i]:.6f}")
        print(f"    k2 (minor conductivity): {seep_data['k2_by_mat'][i]:.6f}")
        print(f"    angle (degrees): {seep_data['angle_by_mat'][i]:.1f}")
        print(f"    kr0 (relative conductivity): {seep_data['kr0_by_mat'][i]:.6f}")
        print(f"    h0 (suction head): {seep_data['h0_by_mat'][i]:.3f}")
    
    # Element material distribution
    element_materials = seep_data['element_materials']
    unique_materials, counts = np.unique(element_materials, return_counts=True)
    print(f"\nElement material distribution:")
    for mat_id, count in zip(unique_materials, counts):
        print(f"  Material {mat_id}: {count} elements")
    
    # Check for potential issues
    print(f"\nData validation:")
    if np.any(seep_data['k1_by_mat'] <= 0):
        print("  WARNING: Some k1 values are <= 0")
    if np.any(seep_data['k2_by_mat'] <= 0):
        print("  WARNING: Some k2 values are <= 0")
    if np.any(seep_data['k1_by_mat'] < seep_data['k2_by_mat']):
        print("  WARNING: Some k1 values are less than k2 (should be major >= minor)")
    
    # Flow type determination
    is_unconfined = np.any(bc_type == 2)
    flow_type = "unconfined" if is_unconfined else "confined"
    print(f"  Flow type: {flow_type}")
    
    print("="*60 + "\n")

def save_seep_data_to_json(seep_data, filename):
    """Save seep_data dictionary to JSON file."""
    import json
    import numpy as np
    
    # Convert numpy arrays to lists for JSON serialization. The transient BC series
    # bindings are internal runtime structures (lists of dicts holding numpy node
    # masks / unit-flux vectors) that are not JSON-serializable and are rebuilt by
    # build_seep_data — skip them rather than crash. (Transient-solution persistence
    # is a separate format; see plan_transient_seep.md §4.)
    _skip = {"head_series_bindings", "flux_series_bindings", "_preflight"}
    seep_data_json = {}
    for key, value in seep_data.items():
        if key in _skip:
            continue
        if isinstance(value, np.ndarray):
            seep_data_json[key] = value.tolist()
        else:
            seep_data_json[key] = value
    
    with open(filename, 'w') as f:
        json.dump(seep_data_json, f, indent=2)
    
    print(f"Seepage data saved to {filename}")

def load_seep_data_from_json(filename):
    """Load seep_data dictionary from JSON file."""
    import json
    import numpy as np
    
    with open(filename, 'r') as f:
        seep_data_json = json.load(f)
    
    # Convert lists back to numpy arrays
    seep_data = {}
    for key, value in seep_data_json.items():
        if isinstance(value, list):
            seep_data[key] = np.array(value)
        else:
            seep_data[key] = value
    
    return seep_data
