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
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
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
        # v18 transient storage: None (blank) -> 0 here; the transient path
        # validates presence at load time (fileio), so a 0 only survives on the
        # steady path where it is never read. Never affects steady behavior.
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
    
    # Calculate appropriate tolerance based on mesh size
    # Use a fraction of the typical element size
    x_range = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
    y_range = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    typical_element_size = min(x_range, y_range) / np.sqrt(len(nodes))  # Approximate element size
    tolerance = typical_element_size * 0.1  # 10% of typical element size
    
    logger.debug("Mesh tolerance for boundary conditions: %.6f", tolerance)
    
    # Series bindings for time-varying (transient) BCs. A head/flux VALUE that is a
    # string names a `tseep` series (fileio validated it). These stay OUT of the
    # baked scalar arrays: build_seep_data leaves the nodes free, and
    # run_transient_seepage resolves them per time step (by the head binding's kind:
    # submerged-only reservoir Dirichlet, or plain time-varying head; scaled unit
    # flux loads for fluxes). Empty on every steady file, so the steady solve is
    # bit-identical.
    head_series_bindings = []   # [{"mask": bool(n_nodes), "series": name, "kind": str}]
    flux_series_bindings = []   # [{"unit_flux_nodal": (n_nodes,), "series": name}]

    # Process specified head boundary conditions
    # Vectorized: compute distance from all nodes to each BC line at once
    specified_heads = seepage_bc.get("specified_heads", [])
    for bc in specified_heads:
        head_value = bc["head"]
        coords = bc["coords"]
        # v18: two Dirichlet head types. "head" = a plain Dirichlet held at the
        # value/series at every node of the polyline, at all times. "reservoir" =
        # the submerged-only face: a node is held at the level only while submerged
        # (elevation <= level) and becomes a seepage-exit face above the water line.
        # For a CONSTANT numeric value with a polyline drawn at or below the level
        # (the usual case) the two are identical — every node is submerged, so no
        # exit face is created and both bake the same plain Dirichlet.
        kind = str(bc.get("kind", "head")).strip().lower()

        if len(coords) < 2:
            continue

        # Compute min distance from each node to any segment of the BC line
        seg_coords = np.array(coords)
        dists = _min_distance_to_polyline(nodes, seg_coords)
        mask = dists <= tolerance
        if isinstance(head_value, str):
            # Time-varying head: resolved per step by run_transient_seepage, which
            # applies the submerged-only reservoir rule or the plain head rule by
            # the binding's kind. It stays OUT of the baked scalar arrays here.
            head_series_bindings.append({"mask": mask, "series": head_value,
                                         "kind": kind})
            continue
        if kind == "reservoir":
            # Submerged-only against the constant level: hold submerged nodes,
            # convert above-level nodes to seepage-exit faces (bc_type 2).
            submerged = mask & (nodes[:, 1] <= head_value)
            above = mask & (nodes[:, 1] > head_value)
            bc_type[submerged] = 1
            bc_values[submerged] = head_value
            new_exit = above & (bc_type != 1)
            bc_type[new_exit] = 2
            bc_values[new_exit] = nodes[new_exit, 1]
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
                                         "series": _b["flux"]})
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


def solve_unsaturated(nodes, elements, bc_type, bc_values, kr0=0.001, h0=-1.0,
                      k1_vals=1.0, k2_vals=1.0, angles=0.0,
                      max_iter=400, tol=1e-6, element_types=None,
                      closure_tol=1e-3, vg_a=None, vg_n=None, model=None,
                      flux_nodal=None):
    """
    Iterative FEM solver for unconfined flow using linear kr frontal function.
    Supports triangular and quadrilateral elements with both linear and quadratic shape functions.

    Convergence is a HYBRID test: both the relative head change (max-norm,
    scaled by domain height x tol) and the relative flow-closure error
    (|net inflow - net outflow| / inflow < closure_tol) must be satisfied.
    The head test alone is a numerical stationarity check whose relation to
    mass balance varies from problem to problem; requiring closure directly
    guarantees the reported flowrate balances to closure_tol on every problem.
    
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
    if model is not None and np.isscalar(model):
        model = np.full(len(elements), model)

    # Set convergence tolerance based on domain height
    ymin, ymax = np.min(y), np.max(y)
    eps = (ymax - ymin) * tol

    print("Starting unsaturated flow iteration...")
    print(f"Convergence tolerance: {eps:.6e}")

    # Track convergence history
    residuals = []
    converged = False
    relax = 1.0  # Initial relaxation factor
    prev_residual = float('inf')

    # Precompute the saturated element matrices, kr sampling operators, and
    # COO index arrays once; each iteration below reduces to a vectorized kr
    # average, a per-element scaling, and one coo_matrix construction.
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals, angles)
    data = _assembly_data(asm, p_nodes=h - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, model=model)
    _prev_active = exit_face_active.copy()
    _n_stable = 0   # consecutive iterations with an unchanged exit-face set

    for iteration in range(1, max_iter + 1):

        # Apply boundary conditions: fixed heads plus the active exit-face set
        dir_mask = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
        _require_runtime_dirichlet(dir_mask)
        dir_values = np.where(bc_type == 1, bc_values, y)
        A, b = _dirichlet_system(asm, data, dir_mask, dir_values, neumann=f_ext)

        h_solved = spsolve(A, b)
        h_new = h_solved

        # FORTRAN-style relaxation strategy
        if iteration > 20:
            relax = 0.5
        if iteration > 40:
            relax = 0.2
        if iteration > 60:
            relax = 0.1
        if iteration > 80:
            relax = 0.05
        if iteration > 100:
            relax = 0.02
        if iteration > 120:
            relax = 0.01

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
        corner_candidate = np.zeros(n_nodes, dtype=bool)
        is_corner = (bc_type == 2) & _exit_is_corner
        stay = is_corner & exit_face_active & ~((h_new < y - hyst) | (q > 0))
        turn_on = is_corner & ~exit_face_active & (h_new >= y + hyst) & (q <= 0)
        corner_candidate[stay | turn_on] = True

        new_exit_face_active = np.zeros(n_nodes, dtype=bool)
        new_exit_face_active[_exit_linear_corners] = corner_candidate[_exit_linear_corners]

        for c1, mid, c2 in _exit_quadratic_edges:
            if exit_face_active[mid]:
                mid_candidate = not (h_new[mid] < y[mid] - hyst or q[mid] > 0)
            else:
                mid_candidate = (h_new[mid] >= y[mid] + hyst and q[mid] <= 0)

            edge_active = bool(corner_candidate[c1] and mid_candidate and corner_candidate[c2])
            if edge_active:
                new_exit_face_active[c1] = True
                new_exit_face_active[mid] = True
                new_exit_face_active[c2] = True

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

        newly_active = new_exit_face_active & ~exit_face_active
        h_new[newly_active] = y[newly_active]
        exit_face_active = new_exit_face_active

        n_active_after = np.sum(exit_face_active)

        # Compute relative residual
        residual = np.max(np.abs(h_new - h)) / (np.max(np.abs(h)) + 1e-10)
        residuals.append(residual)

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
        data_chk = _assembly_data(asm, p_nodes=h_solved - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, model=model)
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
        rel_closure = (float(np.sum(np.abs(q_chk[free_mask] - f_ext[free_mask]))) / inflow_pos
                       if inflow_pos > 1e-30 else 0.0)
        set_stable = bool(np.array_equal(exit_face_active, _prev_active))
        _prev_active = exit_face_active.copy()
        _n_stable = _n_stable + 1 if set_stable else 0

        # Matrix for the next iteration, from the relaxed head field
        data = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, model=model)

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

    else:
        print(f"Warning: Did not converge in {max_iter} iterations")
        print("\nConvergence history:")
        for i, r in enumerate(residuals):
            if i % 20 == 0 or i == len(residuals) - 1:
                print(f"  Iteration {i+1}: residual = {r:.6e}")


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

    # Flowrate: the reaction at the specified-head boundaries, plus the flux
    # delivered on the free rows. "Free" is the runtime Dirichlet complement, so it
    # includes the INACTIVE exit-face nodes — their rows stay free, so their loads
    # do enter the solve and are real inflow (rain landing on the unsaturated part
    # of a seepage face infiltrates). Loads on ACTIVE exit nodes were dropped by the
    # solve, so f_eff zeroes them and they are counted nowhere: rain falling on a
    # saturated, free-draining face simply runs off.
    free_final = ~dir_mask
    f_eff = np.where(dir_mask, 0.0, f_ext)
    react_final = q_final - f_eff
    total_inflow = float(np.sum(react_final[(bc_type == 1) & (react_final > 0)]))
    total_inflow += _flux_inflow(f_ext, free_final)

    # Closure report: kr-consistent imbalance at the converged state
    data_chk = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head', vg_a=vg_a, vg_n=vg_n, model=model)
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

def compute_quad8_centroid_pressure(p_nodes, element_nodes):
    """
    Compute pressure at the centroid of a quad8 element using serendipity shape functions.
    At centroid (xi=0, eta=0), only corner nodes contribute equally.
    """
    valid_nodes = element_nodes[:8][element_nodes[:8] != 0]
    p_elem_nodes = p_nodes[valid_nodes]
    # For serendipity quad8 at center, corner nodes have N=1/4, edge nodes have N=0
    if len(valid_nodes) == 8:
        # Corner nodes (0,1,2,3) contribute 1/4 each, edge nodes (4,5,6,7) contribute 0
        return 0.25 * (p_elem_nodes[0] + p_elem_nodes[1] + p_elem_nodes[2] + p_elem_nodes[3])
    else:
        return np.mean(p_elem_nodes)  # Fallback for incomplete elements

def compute_quad9_centroid_pressure(p_nodes, element_nodes):
    """
    Compute pressure at the centroid of a quad9 element using biquadratic shape functions.
    At centroid (xi=0, eta=0), only the center node contributes.
    """
    p_elem_nodes = p_nodes[element_nodes[:9]]
    # For biquadratic quad9 at center, only center node (node 8) has N=1, all others have N=0
    return p_elem_nodes[8]

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


def kr_gardner_vec(p, a, n, kr_min=1e-4):
    """Vectorized Gardner (1958) relative permeability.

        kr = 1 / (1 + a * psi^n),   psi = -p (suction, positive when unsaturated)

    This is the POWER form of Gardner, the one carried as a legacy option by
    SEEP/W and Slide — not Gardner's exponential form kr = exp(alpha*psi), which
    is a different function used mainly to linearize Richards' equation for
    analytical work. Saturated (p >= 0) gives kr = 1.

    Shares the a/n parameter columns with van Genuchten: the two laws never apply
    to the same material, and the template selects between them with `unsat`.
    """
    a = np.asarray(a, dtype=float)
    n = np.asarray(n, dtype=float)
    psi = np.abs(np.minimum(p, 0.0))          # suction; 0 in the saturated zone
    kr = 1.0 / (1.0 + a * psi ** n)
    kr = np.where(p >= 0.0, 1.0, kr)
    return np.clip(kr, kr_min, 1.0)


def kr_vg_vec(p, vg_a, vg_n, kr_min=1e-4):
    """Vectorized van Genuchten–Mualem relative permeability (steady-state form).

    Depends only on alpha (``vg_a``) and n (``vg_n``): the residual/saturated water
    contents scale storage, not kr, and a steady-state solve carries no storage
    term. ``p`` is the pressure head (negative in the unsaturated zone) and
    broadcasts with vg_a/vg_n. A ``kr_min`` floor plus saturation at p>=0 keep the
    function numerically tame near the wet end; because suction is conservatively
    neglected in stability, the floor does not affect stability results.

        Se = [1 + (alpha|psi|)^n]^(-m),  m = 1 - 1/n
        kr = Se^(1/2) [1 - (1 - Se^(1/m))^m]^2
    """
    n = np.maximum(vg_n, 1.0 + 1e-6)             # guard m = 1 - 1/n away from 0
    m = 1.0 - 1.0 / n
    ah = vg_a * np.abs(np.minimum(p, 0.0))       # |alpha*psi|, 0 in the saturated zone
    Se = (1.0 + ah ** n) ** (-m)
    kr = np.sqrt(Se) * (1.0 - (1.0 - Se ** (1.0 / m)) ** m) ** 2
    kr = np.where(p >= 0.0, 1.0, kr)
    return np.clip(kr, kr_min, 1.0)


def kr_relative_vec(p, kr0, h0, vg_a=None, vg_n=None, model=None, kr_min=1e-4):
    """Per-element relative permeability dispatching on the unsaturated model.

    ``model`` is a per-element code array (``KR_LF``/``KR_VG``/``KR_GARD``)
    broadcasting with ``p``. With ``model`` None or all linear-front this returns
    exactly ``kr_frontal_vec`` — so the linear-front path is bit-identical to
    before."""
    lf = kr_frontal_vec(p, kr0, h0)
    if model is None or not np.any(model):
        return lf
    out = lf
    model = np.asarray(model)
    if np.any(model == KR_VG):
        out = np.where(model == KR_VG, kr_vg_vec(p, vg_a, vg_n, kr_min), out)
    if np.any(model == KR_GARD):
        out = np.where(model == KR_GARD, kr_gardner_vec(p, vg_a, vg_n, kr_min), out)
    return out


def kr_relative(p, kr0, h0, vg_a=None, vg_n=None, model=KR_LF, kr_min=1e-4):
    """Scalar relative permeability with model dispatch (linear front, van
    Genuchten or Gardner), for the per-edge flow-potential integration."""
    if model == KR_VG:
        return float(kr_vg_vec(float(p), vg_a, vg_n, kr_min))
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


def storage_potential_vec(p, Ss, Sy, h0, vg_a=None, vg_n=None, model=None):
    """Storage potential Phi(psi) = integral_0^psi S(psi') dpsi' [-], vectorized.

    Phi is the (reference-shifted) water content: ``Phi(psi_new) - Phi(psi_old)``
    is the change in stored water per unit volume, used for the SECANT (Celia
    mixed-form) mass-balance ledger rather than the tangent capacity that scales
    the mass matrix.  For a saturated linear problem ``Phi = Ss*psi`` so the secant
    ledger and the tangent solve agree to machine precision.  It integrates the
    SAME zone-split S as ``storage_capacity_vec`` (elastic Ss below the phreatic
    surface, retention capacity above it), so the ledger tracks the solve.
    """
    Ss = np.asarray(Ss, dtype=float)
    Sy = np.asarray(Sy, dtype=float)
    h0 = np.asarray(h0, dtype=float)
    safe_h0 = np.where(h0 == 0.0, -1.0, h0)
    ah0 = np.abs(safe_h0)
    # Linear-front / Gardner piecewise integral of the zone-split S:
    #   psi >= 0      : Ss*psi
    #   h0 < psi < 0  : (Sy/|h0|)*psi
    #   psi <= h0     : -Sy            (the band is fully drained)
    lf = np.where(
        p >= 0.0, Ss * p,
        np.where(p > h0, (Sy / ah0) * p, -Sy))
    if model is None or not np.any(model):
        return lf
    out = lf
    model = np.asarray(model)
    if np.any(model == KR_VG):
        n = np.maximum(vg_n, 1.0 + 1e-6)
        m = 1.0 - 1.0 / n
        ah = vg_a * np.abs(np.minimum(p, 0.0))
        Se = (1.0 + ah ** n) ** (-m)
        Se = np.where(p >= 0.0, 1.0, Se)
        # integral of the zone-split S: Ss*psi saturated, Sy*(Se-1) unsaturated
        vg = np.where(p >= 0.0, Ss * p, Sy * (Se - 1.0))
        out = np.where(model == KR_VG, vg, out)
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
    if et == 4:
        g = 1/np.sqrt(3)
        gps = [(-g, -g, 1.0), (g, -g, 1.0), (g, g, 1.0), (-g, g, 1.0)]
    else:
        p1 = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        w1 = [5/9, 8/9, 5/9]
        gps = [(xi, eta, w1[i]*w1[j]) for i, xi in enumerate(p1)
               for j, eta in enumerate(p1)]
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
                   vg_a=None, vg_n=None, model=None):
    """Concatenated COO data for the global matrix.

    p_nodes=None -> saturated assembly. Otherwise each element's ke is scaled
    by factor(kr_avg) with kr averaged at the type-specific sampling points
    (mode='head': kr_avg; mode='stream': 1/kr_avg, guarded). ``model``/``vg_a``/
    ``vg_n`` (per-element) select the unsaturated model; with ``model`` None the kr
    is the linear-front function exactly as before."""
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
# Transient mass-matrix assembly (batched, build-once, storage-scaled per step)
#
# Parallels the stiffness path: the unit (density-1) consistent element mass
# matrices M_e = integral(N_i N_j dOmega) are built ONCE per solve, batched over
# each element type; each time step scales them by the element-average storage
# coefficient S (from storage_capacity_vec, sampled at the SAME kr points) exactly
# as kr scales the stiffness.  The default LUMPED mass uses HRZ diagonal scaling
# (row-sum lumping is identical to HRZ on linear tri3/quad4, but row-sum gives a
# zero corner mass on quadratic elements, so HRZ is used to keep the lumped
# diagonal strictly positive on every element type).
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
    """Precompute the reusable transient mass data (parallels ``_build_assembly``).

    Returns per-element-type groups carrying the unit consistent mass ``me``, the
    HRZ-lumped unit diagonal ``me_lump`` (n_e, nn), and the kr-sampling operators
    (so the element-average storage coefficient is averaged at the same points as
    kr), plus the geometric lumped mass ``node_lump`` (n_nodes,) = integral(N_i
    dOmega) used by the mass-balance ledger."""
    elements = np.asarray(elements)
    element_types = np.asarray(element_types)
    groups = []
    node_lump = np.zeros(len(nodes))
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
        # geometric row-sum lump (integral of N_i) for the mass-balance ledger
        geo = me.sum(axis=2)                             # (n_e, nn)
        np.add.at(node_lump, conn, geo)
        N_kr, w_kr = _kr_sampling(int(et))
        groups.append({'et': int(et), 'idx': idx, 'conn': conn,
                       'me': me, 'me_lump': me_lump,
                       'N_kr': N_kr, 'w_kr': w_kr, 'wsum': w_kr.sum()})
    return {'groups': groups, 'n_nodes': len(nodes), 'node_lump': node_lump}


def _lumped_mass_diag(masm, S_elem_by_group):
    """Assemble the global LUMPED mass diagonal (n_nodes,) from per-group element
    storage coefficients.  ``S_elem_by_group`` is a list aligned to
    ``masm['groups']`` giving the element-average storage S for each element."""
    n = masm['n_nodes']
    diag = np.zeros(n)
    for g, S in zip(masm['groups'], S_elem_by_group):
        contrib = S[:, None] * g['me_lump']             # (n_e, nn)
        np.add.at(diag, g['conn'], contrib)
    return diag


def _mass_storage_by_group(masm, p_nodes, Ss, Sy, h0, vg_a, vg_n, model):
    """Element-average storage coefficient S per group, sampled at the kr points
    and quadrature-averaged (mirrors ``_assembly_data``'s kr averaging).  Returns a
    list aligned to ``masm['groups']``."""
    out = []
    for g in masm['groups']:
        idx = g['idx']
        p_gp = p_nodes[g['conn']] @ g['N_kr'].T          # (n_e, n_pts)
        S_gp = storage_capacity_vec(
            p_gp, Ss[idx][:, None], Sy[idx][:, None], h0[idx][:, None],
            None if vg_a is None else vg_a[idx][:, None],
            None if vg_n is None else vg_n[idx][:, None],
            None if model is None else model[idx][:, None])
        out.append((S_gp @ g['w_kr']) / g['wsum'])
    return out


def diagnose_exit_face(nodes, bc_type, h, q, bc_values):
    """
    Diagnostic function to understand exit face behavior
    """

    print("\n=== Exit Face Diagnostics ===")
    exit_nodes = np.where(bc_type == 2)[0]
    y = nodes[:, 1]

    print(f"Total exit face nodes: {len(exit_nodes)}")
    print("\nNode | x      | y      | h      | h-y    | q        | Status")
    print("-" * 65)

    for node in exit_nodes:
        x_coord = nodes[node, 0]
        y_coord = y[node]
        head = h[node]
        pressure = head - y_coord
        flow = q[node]

        if head >= y_coord:
            status = "SATURATED"
        else:
            status = "UNSATURATED"

        print(f"{node:4d} | {x_coord:6.2f} | {y_coord:6.2f} | {head:6.3f} | {pressure:6.3f} | {flow:8.3e} | {status}")

    # Summary statistics
    saturated = np.sum(h[exit_nodes] >= y[exit_nodes])
    print(f"\nSaturated nodes: {saturated}/{len(exit_nodes)}")

    # Check phreatic surface
    print("\n=== Phreatic Surface Location ===")
    # Find where the phreatic surface intersects the exit face
    for i in range(len(exit_nodes) - 1):
        n1, n2 = exit_nodes[i], exit_nodes[i + 1]
        if (h[n1] >= y[n1]) and (h[n2] < y[n2]):
            # Interpolate intersection point
            y1, y2 = y[n1], y[n2]
            h1, h2 = h[n1], h[n2]
            y_intersect = y1 + (y2 - y1) * (h1 - y1) / (h1 - y1 - h2 + y2)
            print(f"Phreatic surface exits between nodes {n1} and {n2}")
            print(f"Approximate exit elevation: {y_intersect:.3f}")
            break

def create_flow_potential_bc(nodes, elements, q, debug=False, element_types=None, total_flow=None):
    """
    Generates Dirichlet BCs for flow potential φ by marching around the boundary
    and accumulating q to assign φ, ensuring closed-loop conservation.

    Improved version that handles numerical noise and different boundary types.
    Supports both triangular and quadrilateral elements.

    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 3 or 4) triangle or quad node indices
        q : (n_nodes,) nodal flow vector
        debug : bool, if True prints detailed diagnostic information
        element_types : (n_elements,) array indicating 3 for triangles, 4 for quads

    Returns:
        List of (node_id, phi_value) tuples
    """

    from collections import defaultdict

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    if debug:
        print("=== FLOW POTENTIAL BC DEBUG ===")

    # Step 1: Build corner-only edge dictionary and a midside node map.
    # Corner-only edges give a smooth boundary walk. For quadratic elements,
    # the edge-net q (sum of corner + midside q) is used for phi accumulation
    # to avoid oscillations from consistent nodal forces. Midside nodes get
    # phi interpolated from their adjacent corners afterwards.
    edge_counts = defaultdict(list)
    midside_map = {}  # edge_key (sorted corner pair) -> midside node index
    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]

        if element_type == 3:
            i, j, k = element_nodes[:3]
            edges = [(i, j), (j, k), (k, i)]
        elif element_type == 6:
            i, j, k = element_nodes[:3]
            edges = [(i, j), (j, k), (k, i)]
            midside_map[tuple(sorted((i, j)))] = element_nodes[3]
            midside_map[tuple(sorted((j, k)))] = element_nodes[4]
            midside_map[tuple(sorted((k, i)))] = element_nodes[5]
        elif element_type == 4:
            i, j, k, l = element_nodes[:4]
            edges = [(i, j), (j, k), (k, l), (l, i)]
        elif element_type == 8:
            i, j, k, l = element_nodes[:4]
            edges = [(i, j), (j, k), (k, l), (l, i)]
            midside_map[tuple(sorted((i, j)))] = element_nodes[4]
            midside_map[tuple(sorted((j, k)))] = element_nodes[5]
            midside_map[tuple(sorted((k, l)))] = element_nodes[6]
            midside_map[tuple(sorted((l, i)))] = element_nodes[7]
        elif element_type == 9:
            i, j, k, l = element_nodes[:4]
            edges = [(i, j), (j, k), (k, l), (l, i)]
            midside_map[tuple(sorted((i, j)))] = element_nodes[4]
            midside_map[tuple(sorted((j, k)))] = element_nodes[5]
            midside_map[tuple(sorted((k, l)))] = element_nodes[6]
            midside_map[tuple(sorted((l, i)))] = element_nodes[7]
        else:
            continue

        for a, b in edges:
            edge = tuple(sorted((a, b)))
            edge_counts[edge].append(idx)

    # Step 2: Extract boundary edges (appear only once)
    boundary_edges = [edge for edge, elems in edge_counts.items() if len(elems) == 1]

    if debug:
        print(f"Found {len(boundary_edges)} boundary edges")

    # Step 3: Build connectivity for the boundary edges
    neighbor_map = defaultdict(list)
    for a, b in boundary_edges:
        neighbor_map[a].append(b)
        neighbor_map[b].append(a)

    # Step 4: Walk the boundary in order (clockwise or counterclockwise)
    start_node = boundary_edges[0][0]
    ordered_nodes = [start_node]
    visited = {start_node}
    current = start_node

    while True:
        neighbors = [n for n in neighbor_map[current] if n not in visited]
        if not neighbors:
            break
        next_node = neighbors[0]
        ordered_nodes.append(next_node)
        visited.add(next_node)
        current = next_node
        if next_node == start_node:
            break  # closed loop

    # Debug boundary flow statistics
    if debug:
        boundary_nodes = sorted(set(ordered_nodes))
        print(f"Boundary nodes: {len(boundary_nodes)}")
        print(f"Flow statistics on boundary:")
        q_boundary = [q[node] for node in boundary_nodes]
        print(f"  Min q: {min(q_boundary):.6e}")
        print(f"  Max q: {max(q_boundary):.6e}")
        print(f"  Mean |q|: {np.mean([abs(qval) for qval in q_boundary]):.6e}")
        print(f"  Std |q|: {np.std([abs(qval) for qval in q_boundary]):.6e}")

        # Count "small" flows
        thresholds = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4]
        for thresh in thresholds:
            count = sum(1 for qval in q_boundary if abs(qval) < thresh)
            print(f"  Nodes with |q| < {thresh:.0e}: {count}/{len(boundary_nodes)}")

    # Step 5: Find starting point - improved algorithm
    start_idx = None
    n = len(ordered_nodes)

    # Define threshold for "effectively zero" flow based on the magnitude of flows
    q_boundary = [abs(q[node]) for node in ordered_nodes]
    q_max = max(q_boundary) if q_boundary else 1.0
    q_threshold = max(1e-10, q_max * 1e-6)  # Adaptive threshold

    if debug:
        print(f"Flow analysis: max |q| = {q_max:.3e}, threshold = {q_threshold:.3e}")

    # Find the boundary node with maximum positive flow (main inlet)
    max_positive_q = -float('inf')
    max_positive_idx = None
    
    for i in range(n):
        node = ordered_nodes[i]
        if q[node] > max_positive_q:
            max_positive_q = q[node]
            max_positive_idx = i
    
    if max_positive_idx is not None and max_positive_q > q_threshold:
        start_idx = max_positive_idx
        if debug:
            print(f"Starting at maximum inflow node {ordered_nodes[start_idx]} (q = {max_positive_q:.6f})")
    else:
        # Fallback: start at first node
        start_idx = 0
        if debug:
            print(f"No significant positive flow found, starting at first boundary node {ordered_nodes[start_idx]}")

    # Step 6: Assign flow potential values by walking from inlet to exit.
    # Use edge-net q (corner + midside) for accumulation so that quadratic
    # element oscillations cancel out instead of corrupting phi.
    phi = {}

    # Build edge-net q: for each boundary edge, sum q at corner + midside nodes
    edge_net_q = {}  # corner node -> net q for the edge starting at that corner
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        mid = midside_map.get(ekey)
        net = q[c1] + (q[mid] if mid is not None else 0)
        edge_net_q[c1] = net

    # Use known flowrate if provided, otherwise sum positive edge-net flows
    total_q = total_flow if total_flow is not None else sum(v for v in edge_net_q.values() if v > 0)
    phi_val = total_q  # Start with total flow at inlet

    if debug:
        print(f"Starting flow potential calculation at node {ordered_nodes[start_idx]}")
        print(f"Total positive flow: {total_q:.6f}, starting phi: {phi_val:.6f}")

    for i in range(n):
        idx = (start_idx + i) % n
        node = ordered_nodes[idx]
        phi[node] = phi_val
        phi_val -= edge_net_q[node]  # Subtract edge-net flow

        if debug and (i < 5 or i >= n - 5):
            print(f"  Node {node}: φ = {phi[node]:.6f}, edge_net_q = {edge_net_q[node]:.6f}")

    # Interpolate phi at midside boundary nodes from adjacent corners
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        mid = midside_map.get(ekey)
        if mid is not None:
            phi[mid] = 0.5 * (phi[c1] + phi[c2])

    # Check closure
    starting_phi = phi[ordered_nodes[start_idx]]
    closure_error = phi_val - starting_phi

    # Use a relative threshold based on total positive boundary flow
    rel_tol = 1e-2  # 1%
    scale = max(total_q, 1e-12)
    
    if debug or abs(closure_error) > rel_tol * scale:
        print(f"Flow potential closure check: error = {closure_error:.6e}")

        if abs(closure_error) > rel_tol * scale:
            print(f"Warning: Large flow potential closure error = {closure_error:.6e}")
            print("This may indicate:")
            print("  - Non-conservative flow field")
            print("  - Incorrect boundary identification")
            print("  - Numerical issues in the flow solution")

    if debug:
        print("✓ Flow potential BC creation succeeded")

    return list(phi.items())

def create_flow_potential_bc_from_elements(nodes, elements, element_types, head,
                                           k1_vals, k2_vals, angles,
                                           kr0=None, h0=None, total_flow=None,
                                           bc_type=None, exit_face_active=None,
                                           vg_a=None, vg_n=None, model=None):
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
            mdl_e0 = model[elem_idx] if hasattr(model, '__len__') else (model or KR_LF)
            kr_elem = kr_relative(p_elem, kr_e0, h0_e0, vga_e0, vgn_e0, mdl_e0)

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
                    kr_gp = kr_relative(p_gp, kr_e0, h0_e0, vga_e0, vgn_e0, mdl_e0) if kr0 is not None else 1.0
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


def create_flow_potential_bc_from_velocity(nodes, elements, velocity, element_types=None):
    """
    Generates Dirichlet BCs for flow potential φ by integrating the velocity
    flux (v·n) along boundary edges instead of using FEM reaction forces.

    This avoids consistent-force artifacts (large spikes at sharp transitions)
    that corrupt phi when using q = K @ h at boundary nodes.

    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, n) array of element node indices
        velocity : (n_nodes, 2) array of velocity vectors at each node
        element_types : (n_elements,) array of element type codes

    Returns:
        List of (node_id, phi_value) tuples
    """
    from collections import defaultdict

    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Step 1: Build corner-only boundary edges and midside map
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

    boundary_edges = [e for e, elems in edge_counts.items() if len(elems) == 1]

    # Step 2: Walk boundary (corner nodes only)
    neighbor_map = defaultdict(list)
    for a, b in boundary_edges:
        neighbor_map[a].append(b)
        neighbor_map[b].append(a)

    start_node = boundary_edges[0][0]
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

    # Step 3: Compute stream function change (Δψ) along each boundary segment.
    # Uses the identity: Δψ = ∫(-vy·dx + vx·dy) along the segment.
    # This avoids computing an outward normal — the sign is determined
    # entirely by the velocity field and the walk direction.
    # For CCW walks: Δψ > 0 at inflow, < 0 at outflow.
    # For CW walks: signs are flipped; corrected in step 4.
    segment_flux = {}
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        mid = midside_map.get(ekey)

        dx = nodes[c2, 0] - nodes[c1, 0]
        dy = nodes[c2, 1] - nodes[c1, 1]

        # f(node) = -vy * dx + vx * dy  (integrand for Δψ along segment)
        f1 = -velocity[c1, 1] * dx + velocity[c1, 0] * dy
        f2 = -velocity[c2, 1] * dx + velocity[c2, 0] * dy

        if mid is not None:
            f_mid = -velocity[mid, 1] * dx + velocity[mid, 0] * dy
            segment_flux[c1] = (f1 + 4.0 * f_mid + f2) / 6.0
        else:
            segment_flux[c1] = (f1 + f2) / 2.0

    # Step 4: Ensure sign convention: segment_flux > 0 at inflow.
    # For CCW walk, Δψ > 0 at inflow. For CW, it's flipped.
    # Determine walk orientation from the signed area of the boundary polygon.
    signed_area = 0.0
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        signed_area += nodes[c1, 0] * nodes[c2, 1] - nodes[c2, 0] * nodes[c1, 1]
    if signed_area < 0:
        # CW walk — flip signs so segment_flux > 0 at inflow
        segment_flux = {k: -v for k, v in segment_flux.items()}

    # Step 5: Find starting point (maximum inflow segment)
    start_idx = max(range(n), key=lambda i: segment_flux.get(ordered_nodes[i], 0))

    # Step 6: Accumulate phi
    total_q = sum(f for f in segment_flux.values() if f > 0)
    phi_val = total_q
    phi = {}

    for i in range(n):
        idx = (start_idx + i) % n
        node = ordered_nodes[idx]
        phi[node] = phi_val
        phi_val -= segment_flux[node]

    # Interpolate phi at midside boundary nodes
    for i in range(n):
        c1 = ordered_nodes[i]
        c2 = ordered_nodes[(i + 1) % n]
        ekey = tuple(sorted((c1, c2)))
        mid = midside_map.get(ekey)
        if mid is not None:
            phi[mid] = 0.5 * (phi[c1] + phi[c2])

    closure_error = phi_val - phi[ordered_nodes[start_idx]]
    rel_tol = 1e-2
    scale = max(total_q, 1e-12)
    if abs(closure_error) > rel_tol * scale:
        print(f"Flow potential closure check (velocity-based): error = {closure_error:.6e}")

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

def solve_flow_function_unsaturated(nodes, elements, head, k1_vals, k2_vals, angles, kr0, h0, dirichlet_nodes, element_types=None, vg_a=None, vg_n=None, model=None):
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
    if model is not None:
        model = np.broadcast_to(np.asarray(model), (len(elements),))

    # Stream-function coefficient matrix is K/det(K); each element is scaled
    # by 1/kr_avg (mode='stream'), matching the per-element *_kr functions.
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals,
                          angles, flow=True)
    data = _assembly_data(asm, p_nodes=p_nodes, kr0=kr0, h0=h0, mode='stream',
                          vg_a=vg_a, vg_n=vg_n, model=model)

    dir_mask = np.zeros(n_nodes, dtype=bool)
    dir_values = np.zeros(n_nodes)
    for node, phi_value in dirichlet_nodes:
        dir_mask[node] = True
        dir_values[node] = phi_value
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

    phi = spsolve(A, b)
    return phi


def compute_velocity(nodes, elements, head, k1_vals, k2_vals, angles, kr0=None, h0=None, element_types=None, vg_a=None, vg_n=None, model=None):
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
        element_types : (n_elements,) array indicating 3 for triangles, 4 for quads
    
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
                kr_e = kr_relative_vec(p_all[conn].mean(axis=1), kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(model, idx))
            else:
                kr_e = np.ones(len(idx))
            v_e = -kr_e[:, None] * np.einsum('exy,ey->ex', K, grad_h)
            v_e[~ok] = 0.0
            np.add.at(velocity, conn.ravel(), np.repeat(v_e, 3, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(ok.astype(float), 3))

        elif et == 4:
            if use_kr:
                kr_e = kr_relative_vec(p_all[conn].mean(axis=1), kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(model, idx))
            else:
                kr_e = np.ones(len(idx))
            g = 1/np.sqrt(3)
            v_sum = np.zeros((len(idx), 2))
            n_ok = np.zeros(len(idx))
            for xi, eta in [(-g, -g), (g, -g), (g, g), (-g, g)]:
                dxi, deta = _quad_dshape(4, xi, eta)
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
            np.add.at(velocity, conn.ravel(), np.repeat(v_sum, 4, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(n_ok, 4))

        elif et == 6:
            x, y = coords[:, :, 0], coords[:, :, 1]
            detJ = (x[:, 0]-x[:, 2])*(y[:, 1]-y[:, 2]) - (x[:, 1]-x[:, 2])*(y[:, 0]-y[:, 2])
            ok = np.abs(detJ) > 1e-10
            d_safe = np.where(ok, detJ, 1.0)
            # Jinv of the constant corner-node Jacobian
            Ji00 = (y[:, 1]-y[:, 2]) / d_safe
            Ji01 = -(x[:, 1]-x[:, 2]) / d_safe
            Ji10 = -(y[:, 0]-y[:, 2]) / d_safe
            Ji11 = (x[:, 0]-x[:, 2]) / d_safe
            if use_kr:
                # kr at the element centroid via quadratic shape functions
                Nc = np.array([1/3*(2/3-1)]*3 + [4/9]*3)
                p_c = (p_all[conn] * Nc[None, :]).sum(axis=1)
                kr_e = kr_relative_vec(p_c, kr0[idx], h0[idx], _idx_or_none(vg_a, idx), _idx_or_none(vg_n, idx), _idx_or_none(model, idx))
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
        element_types : (n_elements,) array indicating 3 for triangles, 4 for quads
    
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

        elif et == 4:
            g = 1/np.sqrt(3)
            i_sum = np.zeros((len(idx), 2))
            n_ok = np.zeros(len(idx))
            for xi, eta in [(-g, -g), (g, -g), (g, g), (-g, g)]:
                dxi, deta = _quad_dshape(4, xi, eta)
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
            np.add.at(gradient, conn.ravel(), np.repeat(i_sum, 4, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(n_ok, 4))

        elif et == 6:
            x, y = coords[:, :, 0], coords[:, :, 1]
            detJ = (x[:, 0]-x[:, 2])*(y[:, 1]-y[:, 2]) - (x[:, 1]-x[:, 2])*(y[:, 0]-y[:, 2])
            ok = np.abs(detJ) > 1e-10
            d_safe = np.where(ok, detJ, 1.0)
            Ji00 = (y[:, 1]-y[:, 2]) / d_safe
            Ji01 = -(x[:, 1]-x[:, 2]) / d_safe
            Ji10 = -(y[:, 0]-y[:, 2]) / d_safe
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

    count[count == 0] = 1  # Avoid division by zero
    gradient /= count[:, None]
    return gradient

def tri3_stiffness_matrix(nodes_elem, Kmat):
    """
    Compute the 3x3 local stiffness matrix for a 3-node triangular element.
    
    Args:
        nodes_elem: (3,2) array of nodal coordinates
        Kmat: (2,2) conductivity matrix for the element
    Returns:
        ke: (3,3) element stiffness matrix
    """
    xi, yi = nodes_elem[0]
    xj, yj = nodes_elem[1]
    xk, yk = nodes_elem[2]
    
    area = 0.5 * np.linalg.det([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
    if area <= 0:
        return np.zeros((3, 3))

    beta = np.array([yj - yk, yk - yi, yi - yj])
    gamma = np.array([xk - xj, xi - xk, xj - xi])
    grad = np.array([beta, gamma]) / (2 * area)

    ke = area * grad.T @ Kmat @ grad
    return ke


def tri6_stiffness_matrix(nodes_elem, Kmat):
    """
    Compute the 6x6 local stiffness matrix for a 6-node quadratic triangular element.
    Uses 3-point Gaussian quadrature and quadratic shape functions.
    
    GMSH tri6 node ordering:
    0,1,2: corner vertices
    3: midpoint of edge 0-1
    4: midpoint of edge 1-2  
    5: midpoint of edge 2-0
    
    Args:
        nodes_elem: (6,2) array of nodal coordinates
        Kmat: (2,2) conductivity matrix for the element
    Returns:
        ke: (6,6) element stiffness matrix
    """
    # 3-point Gauss quadrature for triangles (exact for degree 2 polynomials)
    gauss_pts = [(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)]
    weights = [1/3, 1/3, 1/3]  # Standard weights for unit triangle
    
    ke = np.zeros((6, 6))
    
    for (L1, L2, L3), w in zip(gauss_pts, weights):
        # Quadratic shape functions in area coordinates for standard GMSH tri6 ordering
        # N0 = L1*(2*L1-1), N1 = L2*(2*L2-1), N2 = L3*(2*L3-1)
        # N3 = 4*L1*L2 (edge 0-1), N4 = 4*L2*L3 (edge 1-2), N5 = 4*L3*L1 (edge 2-0)
        
        # Shape function derivatives w.r.t. area coordinates (standard GMSH ordering)
        dN_dL1 = np.array([4*L1-1, 0, 0, 4*L2, 0, 4*L3])  # dN/dL1
        dN_dL2 = np.array([0, 4*L2-1, 0, 4*L1, 4*L3, 0])  # dN/dL2
        dN_dL3 = np.array([0, 0, 4*L3-1, 0, 4*L2, 4*L1])  # dN/dL3
        
        # Transform from area coordinates to Cartesian coordinates
        # We need the Jacobian: J = [dx/dL1, dx/dL2; dy/dL1, dy/dL2]
        # where L3 = 1 - L1 - L2 is eliminated
        
        # Calculate coordinate derivatives directly from nodal coordinates (now properly oriented)
        x0, y0 = nodes_elem[0]  # Vertex L1=1
        x1, y1 = nodes_elem[1]  # Vertex L2=1  
        x2, y2 = nodes_elem[2]  # Vertex L3=1
        
        # Jacobian matrix (from area to global coordinates)
        # Since x = L1*x0 + L2*x1 + L3*x2 and L3 = 1-L1-L2:
        # dx/dL1 = x0-x2, dx/dL2 = x1-x2, dy/dL1 = y0-y2, dy/dL2 = y1-y2  
        J = np.array([[x0 - x2, x1 - x2],
                      [y0 - y2, y1 - y2]])
        
        detJ = np.linalg.det(J)
        if abs(detJ) < 1e-10:
            continue
        
        # Handle clockwise node ordering by using signed determinant
        # If detJ < 0, the nodes are ordered clockwise, but we still need proper transformation
        Jinv = np.linalg.inv(J)
        
        # Transform shape function derivatives from area coordinates to global coordinates
        # Use direct method based on area coordinate derivatives
        
        # Total triangle area
        total_area = 0.5 * abs(detJ)
        
        # Direct computation of area coordinate derivatives (exact formulas)
        dL1_dx = (y1 - y2) / (2 * total_area)
        dL1_dy = (x2 - x1) / (2 * total_area)
        dL2_dx = (y2 - y0) / (2 * total_area)
        dL2_dy = (x0 - x2) / (2 * total_area)
        dL3_dx = (y0 - y1) / (2 * total_area)
        dL3_dy = (x1 - x0) / (2 * total_area)
        
        # Transform to global coordinates using chain rule:
        # dNi/dx = (dNi/dL1)*(dL1/dx) + (dNi/dL2)*(dL2/dx) + (dNi/dL3)*(dL3/dx)
        # dNi/dy = (dNi/dL1)*(dL1/dy) + (dNi/dL2)*(dL2/dy) + (dNi/dL3)*(dL3/dy)
        
        gradN = np.zeros((2, 6))  # [dN/dx; dN/dy] for 6 shape functions
        
        for i in range(6):
            gradN[0, i] = dN_dL1[i]*dL1_dx + dN_dL2[i]*dL2_dx + dN_dL3[i]*dL3_dx  # dNi/dx
            gradN[1, i] = dN_dL1[i]*dL1_dy + dN_dL2[i]*dL2_dy + dN_dL3[i]*dL3_dy  # dNi/dy
        
        # Element stiffness contribution at this Gauss point
        # Scale by triangle area (detJ = 2 * area for area coordinate mapping)
        triangle_area = 0.5 * abs(detJ)
        ke += (gradN.T @ Kmat @ gradN) * triangle_area * w
    
    return ke


def quad8_stiffness_matrix(nodes_elem, Kmat):
    """
    Compute the 8x8 local stiffness matrix for an 8-node serendipity quadrilateral element.
    Uses 3x3 Gaussian quadrature and serendipity shape functions.
    
    Args:
        nodes_elem: (8,2) array of nodal coordinates
        Kmat: (2,2) conductivity matrix for the element
    Returns:
        ke: (8,8) element stiffness matrix
    """
    # 3x3 Gauss quadrature points and weights
    pts_1d = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    wts_1d = [5/9, 8/9, 5/9]
    
    ke = np.zeros((8, 8))
    
    for i, xi in enumerate(pts_1d):
        for j, eta in enumerate(pts_1d):
            w = wts_1d[i] * wts_1d[j]
            
            # Serendipity shape function derivatives for CCW node ordering
            # Corner nodes: 0(-1,-1), 1(1,-1), 2(1,1), 3(-1,1) 
            # Edge nodes: 4(0,-1), 5(1,0), 6(0,1), 7(-1,0)
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
            
            # Jacobian
            J = np.zeros((2, 2))
            for a in range(8):
                J[0,0] += dN_dxi[a] * nodes_elem[a,0]
                J[0,1] += dN_dxi[a] * nodes_elem[a,1]
                J[1,0] += dN_deta[a] * nodes_elem[a,0]
                J[1,1] += dN_deta[a] * nodes_elem[a,1]
            
            detJ = np.linalg.det(J)
            if detJ <= 0:
                continue
                
            Jinv = np.linalg.inv(J)
            
            # Shape function derivatives w.r.t. x,y
            dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
            dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
            gradN = np.vstack((dN_dx, dN_dy))  # shape (2,8)
            
            # Element stiffness contribution at this Gauss point
            ke += (gradN.T @ Kmat @ gradN) * detJ * w
    
    return ke


def tri6_stiffness_matrix_inverse_k(nodes_elem, Kmat_inv):
    """
    Compute the 6x6 local stiffness matrix for a 6-node quadratic triangular element
    using the inverse conductivity matrix (for flow function computation).
    """
    return tri6_stiffness_matrix(nodes_elem, Kmat_inv)


def quad8_stiffness_matrix_inverse_k(nodes_elem, Kmat_inv):
    """
    Compute the 8x8 local stiffness matrix for an 8-node serendipity quadrilateral element
    using the inverse conductivity matrix (for flow function computation).
    """
    return quad8_stiffness_matrix(nodes_elem, Kmat_inv)


def quad9_stiffness_matrix_inverse_k(nodes_elem, Kmat_inv):
    """
    Compute the 9x9 local stiffness matrix for a 9-node Lagrange quadrilateral element
    using the inverse conductivity matrix (for flow function computation).
    """
    return quad9_stiffness_matrix(nodes_elem, Kmat_inv)


def quad9_stiffness_matrix(nodes_elem, Kmat):
    """
    Compute the 9x9 local stiffness matrix for a 9-node Lagrange quadrilateral element.
    Uses 3x3 Gaussian quadrature and biquadratic Lagrange shape functions.
    
    Args:
        nodes_elem: (9,2) array of nodal coordinates
        Kmat: (2,2) conductivity matrix for the element
    Returns:
        ke: (9,9) element stiffness matrix
    """
    # 3x3 Gauss quadrature points and weights
    pts_1d = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    wts_1d = [5/9, 8/9, 5/9]
    
    ke = np.zeros((9, 9))
    
    for i, xi in enumerate(pts_1d):
        for j, eta in enumerate(pts_1d):
            w = wts_1d[i] * wts_1d[j]
            
            # Lagrange shape function derivatives (biquadratic) for CCW node ordering
            # Corner nodes: 0(-1,-1), 1(1,-1), 2(1,1), 3(-1,1)
            # Edge nodes: 4(0,-1), 5(1,0), 6(0,1), 7(-1,0)
            # Center node: 8(0,0)
            dN_dxi = np.array([
                0.25*(2*xi-1)*eta*(eta-1),                      # Node 0: corner (-1,-1)
                0.25*(2*xi+1)*eta*(eta-1),                      # Node 1: corner (1,-1)
                0.25*(2*xi+1)*eta*(eta+1),                      # Node 2: corner (1,1)
                0.25*(2*xi-1)*eta*(eta+1),                      # Node 3: corner (-1,1)
                -xi*eta*(eta-1),                                # Node 4: edge (0,-1)
                0.5*(2*xi+1)*(1-eta*eta),                       # Node 5: edge (1,0)
                -xi*eta*(eta+1),                                # Node 6: edge (0,1)
                0.5*(2*xi-1)*(1-eta*eta),                       # Node 7: edge (-1,0)
                -2*xi*(1-eta*eta)                               # Node 8: center (0,0)
            ])
            
            dN_deta = np.array([
                0.25*xi*(xi-1)*(2*eta-1),                       # Node 0: corner (-1,-1)
                0.25*xi*(xi+1)*(2*eta-1),                       # Node 1: corner (1,-1)
                0.25*xi*(xi+1)*(2*eta+1),                       # Node 2: corner (1,1)
                0.25*xi*(xi-1)*(2*eta+1),                       # Node 3: corner (-1,1)
                0.5*(1-xi*xi)*(2*eta-1),                        # Node 4: edge (0,-1)
                -eta*xi*(xi+1),                                 # Node 5: edge (1,0)
                0.5*(1-xi*xi)*(2*eta+1),                        # Node 6: edge (0,1)
                -eta*xi*(xi-1),                                 # Node 7: edge (-1,0)
                -2*eta*(1-xi*xi)                                # Node 8: center (0,0)
            ])
            
            # Jacobian
            J = np.zeros((2, 2))
            for a in range(9):
                J[0,0] += dN_dxi[a] * nodes_elem[a,0]
                J[0,1] += dN_dxi[a] * nodes_elem[a,1]
                J[1,0] += dN_deta[a] * nodes_elem[a,0]
                J[1,1] += dN_deta[a] * nodes_elem[a,1]
            
            detJ = np.linalg.det(J)
            if detJ <= 0:
                continue
                
            Jinv = np.linalg.inv(J)
            
            # Shape function derivatives w.r.t. x,y
            dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
            dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
            gradN = np.vstack((dN_dx, dN_dy))  # shape (2,9)
            
            # Element stiffness contribution at this Gauss point
            ke += (gradN.T @ Kmat @ gradN) * detJ * w
    
    return ke


def quad4_stiffness_matrix(nodes_elem, Kmat):
    """
    Compute the 4x4 local stiffness matrix for a 4-node quadrilateral element
    using 2x2 Gauss quadrature and bilinear shape functions.
    nodes_elem: (4,2) array of nodal coordinates (in order: [i,j,k,l])
    Kmat: (2,2) conductivity matrix for the element
    Returns:
        ke: (4,4) element stiffness matrix
    """
    # 2x2 Gauss points and weights
    gauss_pts = [(-1/np.sqrt(3), -1/np.sqrt(3)),
                 (1/np.sqrt(3), -1/np.sqrt(3)),
                 (1/np.sqrt(3), 1/np.sqrt(3)),
                 (-1/np.sqrt(3), 1/np.sqrt(3))]
    weights = [1, 1, 1, 1]
    ke = np.zeros((4, 4))
    
    for gp_idx, ((xi, eta), w) in enumerate(zip(gauss_pts, weights)):
        # Shape function derivatives w.r.t. natural coords
        dN_dxi = np.array([
            [-(1-eta),  (1-eta),  (1+eta), -(1+eta)]
        ]) * 0.25
        dN_deta = np.array([
            [-(1-xi), -(1+xi),  (1+xi),  (1-xi)]
        ]) * 0.25
        dN_dxi = dN_dxi.flatten()
        dN_deta = dN_deta.flatten()
        
        # Jacobian
        J = np.zeros((2,2))
        for a in range(4):
            J[0,0] += dN_dxi[a] * nodes_elem[a,0]
            J[0,1] += dN_dxi[a] * nodes_elem[a,1]
            J[1,0] += dN_deta[a] * nodes_elem[a,0]
            J[1,1] += dN_deta[a] * nodes_elem[a,1]
        
        detJ = np.linalg.det(J)
        if detJ <= 0:
            continue
        Jinv = np.linalg.inv(J)
        # Shape function derivatives w.r.t. x,y
        dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
        dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
        gradN = np.vstack((dN_dx, dN_dy))  # shape (2,4)
        # Element stiffness contribution at this Gauss point
        ke += (gradN.T @ Kmat @ gradN) * detJ * w
    
    return ke


def _kr_factor(p, kr0, h0, mode):
    """Compute kr weighting factor at a point for head or stream mode."""
    kr = kr_frontal(p, kr0, h0)
    if mode == 'head':
        return kr
    else:  # stream
        return 1.0 / kr if kr > 1e-12 else 1e10


def tri3_stiffness_matrix_kr(nodes_elem, Kmat, p_elem_nodes, kr0, h0, mode='head'):
    """
    Tri3 element stiffness with high-order kr quadrature.

    Gradient is constant for tri3, so ke = factor * area * grad^T @ K @ grad.
    We use 13-point triangle quadrature to integrate the nonlinear kr function
    over the element area, matching SEEP2D's approach of over-integrating kr.
    kr_avg and 1/kr_avg are used for head/stream to maintain consistency.

    Args:
        nodes_elem: (3,2) nodal coordinates
        Kmat: (2,2) conductivity matrix (Kmat for head, Kmat_flow for stream)
        p_elem_nodes: (3,) nodal pressure values
        kr0, h0: unsaturated parameters
        mode: 'head' (multiply by kr_avg) or 'stream' (multiply by 1/kr_avg)
    """
    xi, yi = nodes_elem[0]
    xj, yj = nodes_elem[1]
    xk, yk = nodes_elem[2]

    area = 0.5 * abs((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi))
    if area <= 0:
        return np.zeros((3, 3))

    beta = np.array([yj - yk, yk - yi, yi - yj])
    gamma = np.array([xk - xj, xi - xk, xj - xi])
    grad = np.array([beta, gamma]) / (2 * area)

    # 7-point symmetric triangle quadrature (degree 5)
    # Over-integrates kr for better resolution of the phreatic transition.
    a1 = 0.059715871789770
    b1 = 0.470142064105115
    a2 = 0.797426985353087
    b2 = 0.101286507323456
    w0 = 0.1125
    w1 = 0.066197076394253
    w2 = 0.062969590272414
    gauss_pts = [
        (1/3, 1/3, 1/3, w0),
        (a1, b1, b1, w1), (b1, a1, b1, w1), (b1, b1, a1, w1),
        (a2, b2, b2, w2), (b2, a2, b2, w2), (b2, b2, a2, w2),
    ]

    # Weighted average of kr (weights sum to 0.5 for unit triangle)
    kr_wsum = 0.0
    wsum = 0.0
    for L1, L2, L3, w in gauss_pts:
        p_gp = L1 * p_elem_nodes[0] + L2 * p_elem_nodes[1] + L3 * p_elem_nodes[2]
        kr_wsum += w * kr_frontal(p_gp, kr0, h0)
        wsum += w
    kr_avg = kr_wsum / wsum

    if mode == 'head':
        factor = kr_avg
    else:  # stream
        factor = 1.0 / kr_avg if kr_avg > 1e-12 else 1e10

    return factor * area * grad.T @ Kmat @ grad


def tri6_stiffness_matrix_kr(nodes_elem, Kmat, p_elem_nodes, kr0, h0, mode='head'):
    """
    Tri6 element stiffness with averaged kr from Gauss points.

    Averages kr at 3 Gauss points (using quadratic shape function interpolation
    of pressure), then uses kr_avg for head and 1/kr_avg for stream. This avoids
    the 1/kr blowup that occurs with per-GP evaluation when individual GPs fall
    deep in the unsaturated zone.

    Args:
        nodes_elem: (6,2) nodal coordinates
        Kmat: (2,2) conductivity matrix (Kmat for head, Kmat_flow for stream)
        p_elem_nodes: (6,) nodal pressure values
        kr0, h0: unsaturated parameters
        mode: 'head' or 'stream'
    """
    # 3-point Gauss quadrature for triangles (exact for degree 2 polynomials)
    gauss_pts = [(1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6)]
    weights = [1/3, 1/3, 1/3]

    ke = np.zeros((6, 6))

    x0, y0 = nodes_elem[0]
    x1, y1 = nodes_elem[1]
    x2, y2 = nodes_elem[2]

    J = np.array([[x0 - x2, x1 - x2],
                  [y0 - y2, y1 - y2]])
    detJ = np.linalg.det(J)
    if abs(detJ) < 1e-10:
        return ke

    total_area = 0.5 * abs(detJ)
    dL1_dx = (y1 - y2) / (2 * total_area)
    dL1_dy = (x2 - x1) / (2 * total_area)
    dL2_dx = (y2 - y0) / (2 * total_area)
    dL2_dy = (x0 - x2) / (2 * total_area)
    dL3_dx = (y0 - y1) / (2 * total_area)
    dL3_dy = (x1 - x0) / (2 * total_area)

    # Average kr across Gauss points first, then apply kr_avg or 1/kr_avg.
    # This matches the tri3 approach and avoids 1/kr blowup when individual
    # GPs fall in the unsaturated zone (where kr → 0 makes 1/kr → ∞).
    kr_wsum = 0.0
    wsum = 0.0
    for (L1, L2, L3), w in zip(gauss_pts, weights):
        N = np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                      4*L1*L2, 4*L2*L3, 4*L3*L1])
        p_gp = N @ p_elem_nodes
        kr_wsum += w * kr_frontal(p_gp, kr0, h0)
        wsum += w
    kr_avg = kr_wsum / wsum

    if mode == 'head':
        factor = kr_avg
    else:  # stream
        factor = 1.0 / kr_avg if kr_avg > 1e-12 else 1e10

    for (L1, L2, L3), w in zip(gauss_pts, weights):
        dN_dL1 = np.array([4*L1-1, 0, 0, 4*L2, 0, 4*L3])
        dN_dL2 = np.array([0, 4*L2-1, 0, 4*L1, 4*L3, 0])
        dN_dL3 = np.array([0, 0, 4*L3-1, 0, 4*L2, 4*L1])

        gradN = np.zeros((2, 6))
        for i in range(6):
            gradN[0, i] = dN_dL1[i]*dL1_dx + dN_dL2[i]*dL2_dx + dN_dL3[i]*dL3_dx
            gradN[1, i] = dN_dL1[i]*dL1_dy + dN_dL2[i]*dL2_dy + dN_dL3[i]*dL3_dy

        ke += (gradN.T @ Kmat @ gradN) * total_area * w

    ke *= factor
    return ke


def quad4_stiffness_matrix_kr(nodes_elem, Kmat, p_elem_nodes, kr0, h0, mode='head'):
    """
    Quad4 element stiffness with averaged kr from Gauss points.

    Uses 4x4 Gauss quadrature (matching SEEP2D's qdflow subroutine) to sample
    kr at 16 interior points, then averages to get kr_avg. Uses kr_avg for head
    and 1/kr_avg for stream to maintain head/stream consistency and avoid
    1/kr blowup at individual Gauss points.

    Args:
        nodes_elem: (4,2) nodal coordinates
        Kmat: (2,2) conductivity matrix (Kmat for head, Kmat_flow for stream)
        p_elem_nodes: (4,) nodal pressure values
        kr0, h0: unsaturated parameters
        mode: 'head' or 'stream'
    """
    # 4-point Gauss rule (matching SEEP2D) for kr sampling
    pts_1d = [-0.86113631, -0.33998104, 0.33998104, 0.86113631]
    wts_1d = [0.34785485, 0.65214516, 0.65214516, 0.34785485]

    # First pass: weighted average of kr at 4x4 Gauss points
    kr_wsum = 0.0
    wsum = 0.0
    for i_gp, xi in enumerate(pts_1d):
        for j_gp, eta in enumerate(pts_1d):
            w = wts_1d[i_gp] * wts_1d[j_gp]
            N = np.array([0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta),
                          0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta)])
            p_gp = N @ p_elem_nodes
            kr_wsum += w * kr_frontal(p_gp, kr0, h0)
            wsum += w
    kr_avg = kr_wsum / wsum

    if mode == 'head':
        factor = kr_avg
    else:  # stream
        factor = 1.0 / kr_avg if kr_avg > 1e-12 else 1e10

    # Second pass: assemble stiffness with standard 2x2 quadrature
    gauss_pts = [(-1/np.sqrt(3), -1/np.sqrt(3)),
                 (1/np.sqrt(3), -1/np.sqrt(3)),
                 (1/np.sqrt(3), 1/np.sqrt(3)),
                 (-1/np.sqrt(3), 1/np.sqrt(3))]
    weights = [1, 1, 1, 1]
    ke = np.zeros((4, 4))

    for (xi, eta), w in zip(gauss_pts, weights):
        dN_dxi = np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)]) * 0.25
        dN_deta = np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)]) * 0.25

        J = np.zeros((2, 2))
        for a in range(4):
            J[0,0] += dN_dxi[a] * nodes_elem[a,0]
            J[0,1] += dN_dxi[a] * nodes_elem[a,1]
            J[1,0] += dN_deta[a] * nodes_elem[a,0]
            J[1,1] += dN_deta[a] * nodes_elem[a,1]

        detJ = np.linalg.det(J)
        if detJ <= 0:
            continue
        Jinv = np.linalg.inv(J)

        dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
        dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
        gradN = np.vstack((dN_dx, dN_dy))

        ke += (gradN.T @ Kmat @ gradN) * detJ * w

    return factor * ke


def quad8_stiffness_matrix_kr(nodes_elem, Kmat, p_elem_nodes, kr0, h0, mode='head'):
    """
    Quad8 (serendipity) element stiffness with averaged kr from Gauss points.

    Args:
        nodes_elem: (8,2) nodal coordinates
        Kmat: (2,2) conductivity matrix (Kmat for head, Kmat_flow for stream)
        p_elem_nodes: (8,) nodal pressure values
        kr0, h0: unsaturated parameters
        mode: 'head' or 'stream'
    """
    pts_1d = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    wts_1d = [5/9, 8/9, 5/9]

    # First pass: weighted average of kr at 3x3 Gauss points
    kr_wsum = 0.0
    wsum = 0.0
    for i_gp, xi in enumerate(pts_1d):
        for j_gp, eta in enumerate(pts_1d):
            w = wts_1d[i_gp] * wts_1d[j_gp]
            N = np.array([
                0.25*(1-xi)*(1-eta)*(-xi-eta-1),
                0.25*(1+xi)*(1-eta)*(xi-eta-1),
                0.25*(1+xi)*(1+eta)*(xi+eta-1),
                0.25*(1-xi)*(1+eta)*(-xi+eta-1),
                0.5*(1-xi*xi)*(1-eta),
                0.5*(1+xi)*(1-eta*eta),
                0.5*(1-xi*xi)*(1+eta),
                0.5*(1-xi)*(1-eta*eta)
            ])
            p_gp = N @ p_elem_nodes
            kr_wsum += w * kr_frontal(p_gp, kr0, h0)
            wsum += w
    kr_avg = kr_wsum / wsum

    if mode == 'head':
        factor = kr_avg
    else:  # stream
        factor = 1.0 / kr_avg if kr_avg > 1e-12 else 1e10

    # Second pass: assemble stiffness
    ke = np.zeros((8, 8))

    for i_gp, xi in enumerate(pts_1d):
        for j_gp, eta in enumerate(pts_1d):
            w = wts_1d[i_gp] * wts_1d[j_gp]

            dN_dxi = np.array([
                -0.25*(1-eta)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta),
                0.25*(1-eta)*(xi-eta-1) + 0.25*(1+xi)*(1-eta),
                0.25*(1+eta)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),
                -0.25*(1+eta)*(-xi+eta-1) - 0.25*(1-xi)*(1+eta),
                -xi*(1-eta),
                0.5*(1-eta*eta),
                -xi*(1+eta),
                -0.5*(1-eta*eta)
            ])
            dN_deta = np.array([
                -0.25*(1-xi)*(-xi-eta-1) - 0.25*(1-xi)*(1-eta),
                -0.25*(1+xi)*(xi-eta-1) - 0.25*(1+xi)*(1-eta),
                0.25*(1+xi)*(xi+eta-1) + 0.25*(1+xi)*(1+eta),
                0.25*(1-xi)*(-xi+eta-1) + 0.25*(1-xi)*(1+eta),
                -0.5*(1-xi*xi),
                -eta*(1+xi),
                0.5*(1-xi*xi),
                -eta*(1-xi)
            ])

            J = np.zeros((2, 2))
            for a in range(8):
                J[0,0] += dN_dxi[a] * nodes_elem[a,0]
                J[0,1] += dN_dxi[a] * nodes_elem[a,1]
                J[1,0] += dN_deta[a] * nodes_elem[a,0]
                J[1,1] += dN_deta[a] * nodes_elem[a,1]

            detJ = np.linalg.det(J)
            if detJ <= 0:
                continue
            Jinv = np.linalg.inv(J)

            dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
            dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
            gradN = np.vstack((dN_dx, dN_dy))

            ke += (gradN.T @ Kmat @ gradN) * detJ * w

    return factor * ke


def quad9_stiffness_matrix_kr(nodes_elem, Kmat, p_elem_nodes, kr0, h0, mode='head'):
    """
    Quad9 (Lagrange) element stiffness with averaged kr from Gauss points.

    Args:
        nodes_elem: (9,2) nodal coordinates
        Kmat: (2,2) conductivity matrix (Kmat for head, Kmat_flow for stream)
        p_elem_nodes: (9,) nodal pressure values
        kr0, h0: unsaturated parameters
        mode: 'head' or 'stream'
    """
    pts_1d = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    wts_1d = [5/9, 8/9, 5/9]

    # First pass: weighted average of kr at 3x3 Gauss points
    kr_wsum = 0.0
    wsum = 0.0
    for i_gp, xi in enumerate(pts_1d):
        for j_gp, eta in enumerate(pts_1d):
            w = wts_1d[i_gp] * wts_1d[j_gp]
            N = np.array([
                0.25*xi*(xi-1)*eta*(eta-1),
                0.25*xi*(xi+1)*eta*(eta-1),
                0.25*xi*(xi+1)*eta*(eta+1),
                0.25*xi*(xi-1)*eta*(eta+1),
                0.5*(1-xi*xi)*eta*(eta-1),
                0.5*xi*(xi+1)*(1-eta*eta),
                0.5*(1-xi*xi)*eta*(eta+1),
                0.5*xi*(xi-1)*(1-eta*eta),
                (1-xi*xi)*(1-eta*eta)
            ])
            p_gp = N @ p_elem_nodes
            kr_wsum += w * kr_frontal(p_gp, kr0, h0)
            wsum += w
    kr_avg = kr_wsum / wsum

    if mode == 'head':
        factor = kr_avg
    else:  # stream
        factor = 1.0 / kr_avg if kr_avg > 1e-12 else 1e10

    # Second pass: assemble stiffness
    ke = np.zeros((9, 9))

    for i_gp, xi in enumerate(pts_1d):
        for j_gp, eta in enumerate(pts_1d):
            w = wts_1d[i_gp] * wts_1d[j_gp]

            dN_dxi = np.array([
                0.25*(2*xi-1)*eta*(eta-1),
                0.25*(2*xi+1)*eta*(eta-1),
                0.25*(2*xi+1)*eta*(eta+1),
                0.25*(2*xi-1)*eta*(eta+1),
                -xi*eta*(eta-1),
                0.5*(2*xi+1)*(1-eta*eta),
                -xi*eta*(eta+1),
                0.5*(2*xi-1)*(1-eta*eta),
                -2*xi*(1-eta*eta)
            ])
            dN_deta = np.array([
                0.25*xi*(xi-1)*(2*eta-1),
                0.25*xi*(xi+1)*(2*eta-1),
                0.25*xi*(xi+1)*(2*eta+1),
                0.25*xi*(xi-1)*(2*eta+1),
                0.5*(1-xi*xi)*(2*eta-1),
                -eta*xi*(xi+1),
                0.5*(1-xi*xi)*(2*eta+1),
                -eta*xi*(xi-1),
                -2*eta*(1-xi*xi)
            ])

            J = np.zeros((2, 2))
            for a in range(9):
                J[0,0] += dN_dxi[a] * nodes_elem[a,0]
                J[0,1] += dN_dxi[a] * nodes_elem[a,1]
                J[1,0] += dN_deta[a] * nodes_elem[a,0]
                J[1,1] += dN_deta[a] * nodes_elem[a,1]

            detJ = np.linalg.det(J)
            if detJ <= 0:
                continue
            Jinv = np.linalg.inv(J)

            dN_dx = Jinv[0,0]*dN_dxi + Jinv[0,1]*dN_deta
            dN_dy = Jinv[1,0]*dN_dxi + Jinv[1,1]*dN_deta
            gradN = np.vstack((dN_dx, dN_dy))

            ke += (gradN.T @ Kmat @ gradN) * detJ * w

    return factor * ke


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
        tol: relative head-change tolerance (scaled by domain height)
        closure_tol: relative flow-closure tolerance for unconfined problems —
            iteration continues until |net inflow - net outflow| / inflow is
            below this, so the reported flowrate balances regardless of how
            the head tolerance maps to mass balance on a given problem
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
    """
    report = seep_data.get("_preflight")
    if report is not None:
        report.raise_for_errors()

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
    bc_type = seep_data["bc_type"]
    bc_values = seep_data["bc_values"]
    # Consistent nodal loads from specified-flux BCs (absent on seep_data built by
    # import_seep2d or older callers -> no flux).
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
            model=model_per_element,
            flux_nodal=flux_nodal,
        )
        # Compute phi BCs from element-level boundary flux
        dirichlet_phi_bcs = create_flow_potential_bc_from_elements(
            nodes, elements, element_types, head, k1, k2, angle,
            kr0=kr0_per_element, h0=h0_per_element, total_flow=total_flow,
            bc_type=bc_type, exit_face_active=exit_face_active,
            vg_a=vg_a_per_element, vg_n=vg_n_per_element, model=model_per_element)
        phi = solve_flow_function_unsaturated(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, dirichlet_phi_bcs, element_types,
                                              vg_a=vg_a_per_element, vg_n=vg_n_per_element, model=model_per_element)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, element_types,
                                    vg_a=vg_a_per_element, vg_n=vg_n_per_element, model=model_per_element)
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
        "closure_error": closure_error,
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

    An exit edge is one whose BOTH corners are exit-face nodes. This is the single
    source of the exit-face edge topology used by BOTH the steady unconfined solver
    (``solve_unsaturated``) and the transient stepper (``run_transient_seepage``);
    the transient rebuilds it per step because a series-head boundary moves the
    exit-face set as the water level changes."""
    n = len(bc_type)
    is_corner = np.zeros(n, dtype=bool)
    linear_corners = np.zeros(n, dtype=bool)
    quadratic_edges = []
    seen = set()
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
            if bc_type[c1] == 2 and bc_type[c2] == 2:
                seen.add(key)
                is_corner[c1] = True
                is_corner[c2] = True
                if len(gnodes) == 3:  # quadratic edge: has a midside node
                    quadratic_edges.append((int(c1), int(gnodes[1]), int(c2)))
                else:
                    linear_corners[c1] = True
                    linear_corners[c2] = True
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


def run_transient_seepage(seep_data, tseep_data, theta=1.0, lumped=True,
                          dt0=None, dt_max=None, dt_min=None, growth=1.3,
                          max_head_change_frac=0.05, picard_tol=1e-5,
                          picard_max=25, h_init=None, verbose=True,
                          progress_callback=None):
    """Transient variably-saturated seepage solver (plan §1/§3/§5).

    Solves ``div(kr K grad h) + Q = S dh/dt`` with the theta-method in time
    (default backward Euler, ``theta=1``) and Picard iteration within each step,
    reusing the steady batched assembly. Storage S is the elastic Ss where saturated
    and the retention curve's specific moisture capacity C(psi) where unsaturated —
    one or the other, never their sum (``storage_capacity_vec``), the same
    saturated-only convention RS2 and SEEP/W use. The lumped mass matrix is scaled by
    the element-average S each iteration, exactly as kr scales the stiffness.
    Time-varying BCs come from the ``tseep`` series bound
    in ``seep_data``, resolved per step by the head binding's kind: a "reservoir"
    series is a submerged-only Dirichlet (a boundary node is held at h(t) only while
    submerged, and becomes an exit-face node above the water line), while a "head"
    series is a plain Dirichlet at h(t) at every node of the polyline at all times
    (it can hold a suction/negative-pressure head and never converts to an exit
    face). Adaptive dt clamps to the saved-frame schedule so every saved frame is a
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
    ``dt_history``, ``mass_balance`` (cumulative + per-frame ledger), ``converged``,
    ``cancelled`` (True only when a ``progress_callback`` requested an abort), and
    provenance (theta, lumped, duration).
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
    model = None if _unsat is None else _unsat[mat_ids]
    vg_a = None if _vga is None else _vga[mat_ids]
    vg_n = None if _vgn is None else _vgn[mat_ids]
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
    # dynamically as a series-bound reservoir head falls below a node. The per-step
    # active-set update (in ``attempt``) resolves the exit face at any element order:
    # corners use the per-corner SEEP2D h/q test and quadratic (tri6/quad8/quad9)
    # exit edges are tracked all-or-nothing via ``_exit_face_topology`` — the same
    # edge logic the steady solver uses — so the seepage-face transition lands on a
    # corner rather than mid-edge.
    # Only a "reservoir" head series can dynamically raise an exit face as the
    # water level falls; a plain "head" series is a Dirichlet at every node at all
    # times and never converts. Static exit faces (bc_type 2) always count.
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
        """BC configuration at time t: (bc_type, bc_values, flux_nodal). Constant
        (numeric) heads/exit faces carry through. A time-varying head series is
        resolved by its kind: a "reservoir" series is Dirichlet at h(t) on submerged
        nodes (elevation <= h(t)) and an exit-face node above the water line; a
        "head" series is a plain Dirichlet at h(t) at every node of the polyline,
        at all times (it can hold a negative-pressure/suction head and never
        converts to an exit face)."""
        bt = base_bc_type.copy()
        bv = base_bc_values.copy()
        for b in head_bindings:
            ht = _eval_series(_series(b["series"]), t)
            mask = b["mask"]
            if str(b.get("kind", "reservoir")).strip().lower() == "head":
                # plain Dirichlet at h(t) everywhere on the polyline
                bt[mask] = 1
                bv[mask] = ht
                continue
            submerged = mask & (y <= ht)
            above = mask & (y > ht)
            bt[submerged] = 1
            bv[submerged] = ht
            # nodes above the water line become potential seepage-exit faces
            new_exit = above & (bt != 1)
            bt[new_exit] = 2
            bv[new_exit] = y[new_exit]
        flux = base_flux.copy()
        for b in flux_bindings:
            flux = flux + _eval_series(_series(b["series"]), t) * b["unit_flux_nodal"]
        return bt, bv, flux

    # Build-once assembly (stiffness + mass).
    asm = _build_assembly(nodes, elements, element_types, k1, k2, angle)
    masm = _build_mass_assembly(nodes, elements, element_types)
    node_lump = masm["node_lump"]

    def kdata(p):
        return _assembly_data(asm, p_nodes=p, kr0=kr0, h0=h0, mode='head',
                              vg_a=vg_a, vg_n=vg_n, model=model)

    def mdiag(p, dt):
        S = _mass_storage_by_group(masm, p, ss, sy, h0, vg_a, vg_n, model)
        return _lumped_mass_diag(masm, S) / dt

    def _frame_flows(h_cur, bt, bv, flux, act):
        """Boundary inflow/outflow for a saved frame, from the FE nodal reaction
        (``K·h`` less any applied flux) at this instant. A cheap by-product of the
        assembly already built for the time step — NOT a dual solve.

        Under storage exchange the boundary inflow and outflow DIFFER (a steady
        solve's single "Total Flowrate" no longer applies): the difference is the
        water going into, or coming out of, storage. Both are carried per frame.
        No stream function is computed — a transient (storage-release) state is not
        divergence-free, so it has no flow net; the frame's flow direction is read
        from velocity vectors instead."""
        p = h_cur - y
        unconf = bool(np.any(bt == 2))
        kd = kdata(p) if unconf else _assembly_data(asm)
        dir_mask = (bt == 1) | ((bt == 2) & act)
        f_eff = np.where(dir_mask, 0.0, flux)
        reaction = _coo_matvec(asm, kd, h_cur) - f_eff
        pos = reaction > 0
        inflow = float(np.sum(reaction[dir_mask & pos])) + _flux_inflow(flux, ~dir_mask)
        outflow = float(-np.sum(reaction[dir_mask & ~pos]))
        return inflow, outflow

    # ---- initial condition ----
    bt0, bv0, flux0 = resolve_bc(0.0)
    if h_init is not None:
        h = np.asarray(h_init, dtype=float).copy()
    else:
        h = _transient_steady_ic(nodes, elements, element_types, bt0, bv0,
                                  k1, k2, angle, kr0, h0, vg_a, vg_n, model,
                                  flux0, verbose)

    # characteristic head scale for the Picard tolerance and the dh limiter
    hmin, hmax = float(np.min(h)), float(np.max(h))
    for b in head_bindings:
        for v in _series(b["series"])[1]:
            hmin, hmax = min(hmin, v), max(hmax, v)
    char_h = max(hmax - hmin, float(np.max(y) - np.min(y)), 1.0)
    picard_eps = picard_tol * char_h
    dh_limit = max_head_change_frac * char_h

    targets, save_interval = _transient_saved_schedule(tseep_data, duration)
    if dt_max is None:
        dt_max = save_interval
    if dt0 is None:
        first = targets[0] if targets else duration
        dt0 = max(min(first, save_interval) / 10.0, duration * 1e-9)
    if dt_min is None:
        dt_min = dt0 * 1e-6

    # Nodal storage parameters for the secant mass-balance ledger (smoothing a
    # per-element scalar to nodes is exact for a single material and adequate
    # otherwise — the ledger is a diagnostic, not the solve).
    ss_nodal = _nodal_from_elem(masm, ss, n_nodes)
    sy_nodal = _nodal_from_elem(masm, sy, n_nodes)
    h0_nodal = _nodal_from_elem(masm, h0, n_nodes)
    vga_nodal = None if vg_a is None else _nodal_from_elem(masm, vg_a, n_nodes)
    vgn_nodal = None if vg_n is None else _nodal_from_elem(masm, vg_n, n_nodes)
    model_nodal = None if model is None else _nodal_from_elem(masm, model, n_nodes,
                                                              mode='max')

    def _stored(hh):
        return float(np.sum(node_lump * storage_potential_vec(
            hh - y, ss_nodal, sy_nodal, h0_nodal, vga_nodal, vgn_nodal, model_nodal)))

    stored0 = _stored(h)   # secant mass-balance baseline
    # Physical storage scale (total lumped volume x max Ss x head span). The
    # relative mass-balance closure is normalized against the LARGER of the actual
    # net cumulative inflow and a small fraction of this scale, so a symmetric
    # relaxation with ~zero NET storage change does not produce a 0/0 ratio.
    stored_scale = float(np.sum(node_lump)) * max(float(np.max(ss)), 1e-30) * char_h

    # IC flow net: submerged exit-face nodes (elevation <= head) are the active
    # (saturated) exit set at t=0; confined problems have no exit nodes so this is moot.
    _act0 = (bt0 == 2) & (h >= y)
    _in0, _out0 = _frame_flows(h, bt0, bv0, flux0, _act0)
    frames = [{"time": 0.0, "head": h.copy(), "u": gamma_w * (h - y),
               "inflow": _in0, "outflow": _out0}]
    dt_history = []
    mb_frames = [{"time": 0.0, "stored_change": 0.0, "cumulative_inflow": 0.0,
                  "closure": 0.0}]
    cum_inflow = 0.0            # tangent storage change accumulated (= net boundary inflow)
    active = np.zeros(n_nodes, dtype=bool)   # exit-face active set, warm-started
    all_converged = True
    #: [cycles detected, edges resolved]. A cycle DETECTED is a sweep in which the
    #: exit-face active set revisited a state it had left; RESOLVED counts the times
    #: the partly-wet fallback actually had to be applied, which happens only on the
    #: retry of a step that already failed. detected > 0 with resolved == 0 is the
    #: ordinary case: the set wobbled and settled on its own. Lists so ``attempt``
    #: can increment them as closures.
    n_cycles = [0]
    n_damped = [0]

    def attempt(h_old, dt, t_new, warm_active, force=False, damp=False):
        """One theta-scheme step with Picard + exit-face active set. Returns
        (ok, h_new, active, picard_count, tangent_storage_change, cycling).
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
        # the whole exit face saturated.
        act[(bt == 2) & exit_on_edge & ~warm_active] = True
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
        pinned = np.zeros(n_nodes, dtype=bool)
        pinned_state = np.zeros(n_nodes, dtype=bool)
        cycling = False
        kold_h = None
        if theta < 1.0:
            kold_h = _coo_matvec(asm, kdata(h_old - y), h_old)
        h_it = h_old.copy()
        pic = 0
        ok = False
        for pic in range(1, picard_max + 1):
            p_it = h_it - y
            kd = kdata(p_it)
            md = mdiag(p_it, dt)
            dir_mask = (bt == 1) | ((bt == 2) & act)
            if not np.any(dir_mask):
                # singular (all exit nodes shed, no fixed head): cannot force.
                return False, h_old, warm_active, pic, 0.0, cycling
            dir_values = np.where(bt == 1, bv, y)
            rhs = md * h_old + theta * flux
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
                is_c = (bt == 2) & exit_is_corner
                stay = is_c & act & ~((h_new < y - hyst) | (q > 0))
                turn_on = is_c & ~act & (h_new >= y + hyst) & (q <= 0)
                corner_candidate = np.zeros(n_nodes, dtype=bool)
                corner_candidate[stay | turn_on] = True

                new_act = np.zeros(n_nodes, dtype=bool)
                # linear exit edges: per-corner activation (tri3/quad4 rule)
                new_act[exit_linear_corners] = corner_candidate[exit_linear_corners]
                # quadratic exit edges: all-or-nothing on (corner, midside, corner)
                mid_cands = []
                for _c1, _mid, _c2 in exit_quadratic_edges:
                    if act[_mid]:
                        mid_candidate = not (h_new[_mid] < y[_mid] - hyst
                                             or q[_mid] > 0)
                    else:
                        mid_candidate = (h_new[_mid] >= y[_mid] + hyst
                                         and q[_mid] <= 0)
                    mid_cands.append(bool(mid_candidate))
                    if corner_candidate[_c1] and mid_candidate and corner_candidate[_c2]:
                        new_act[_c1] = True
                        new_act[_mid] = True
                        new_act[_c2] = True
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
                act = new_act
            else:
                set_stable = True
            change = float(np.max(np.abs(h_new - h_it)))
            h_it = h_new
            if change < picard_eps and set_stable:
                ok = True
                break

        def _tsc():
            # tangent storage change over the step (= net boundary inflow, discretely)
            md_final = mdiag(h_it - y, 1.0)
            return float(np.sum(md_final * (h_it - h_old)))

        if not ok and not force:
            return False, h_it, act, pic, 0.0, cycling
        # head-change limiter (resolve sharp fronts even when Picard is happy).
        # Measured on FREE nodes only: a Dirichlet node's change is imposed by the
        # BC (a stepped reservoir head jumps by the full step regardless of dt), so
        # counting it would force dt -> dt_min forever with no benefit.
        free = ~dir_mask
        if (not force and np.any(free)
                and float(np.max(np.abs((h_it - h_old)[free]))) > dh_limit):
            return False, h_it, act, pic, 0.0, cycling
        return True, h_it, act, pic, _tsc(), cycling

    t = 0.0
    dt = dt0
    cancelled = False
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
            ok, h_new, new_active, pic, tsc, cyc = attempt(h, dt_try, t + dt_try, active)
            # A step that failed WITH the exit-face set cycling is not a step that
            # needs a smaller dt: no dt closes a set-stability test the all-or-nothing
            # edge rule cannot satisfy, which is why this used to collapse to dt_min
            # and force-accept. Retry it once at the SAME dt with the cycle resolved,
            # and only fall back to halving if that still fails.
            if not ok and cyc:
                ok, h_new, new_active, pic, tsc, cyc = attempt(
                    h, dt_try, t + dt_try, active, damp=True)
            while not ok and dt_try > dt_min:
                dt_try *= 0.5
                ok, h_new, new_active, pic, tsc, cyc = attempt(
                    h, dt_try, t + dt_try, active)
                if not ok and cyc:
                    ok, h_new, new_active, pic, tsc, cyc = attempt(
                        h, dt_try, t + dt_try, active, damp=True)
            if not ok:
                all_converged = False
                dt_try = max(dt_try, dt_min)
                if verbose:
                    print(f"  transient step failed to converge at t={t:.6g} "
                          f"(dt down to {dt_try:.3g}); force-accepting last iterate")
                # force-accept the smallest-dt iterate to make progress (flagged
                # non-converged). A singular system still can't be forced -> break.
                ok, h_new, new_active, pic, tsc, cyc = attempt(
                    h, dt_try, t + dt_try, active, force=True, damp=True)
                if not ok:
                    raise SeepInputError(
                        f"Transient solve became singular at t={t:.6g} (the exit "
                        f"face fully deactivated with no specified head). Add a "
                        f"specified head or check the exit-face definition.")
            t += dt_try
            h = h_new
            active = new_active
            cum_inflow += tsc
            dt_history.append(dt_try)
            # Progress + cancellation: report the simulated-time fraction and honor an
            # abort request (a falsy return) after this accepted step. Zero cost when
            # no callback is supplied.
            if progress_callback is not None and not progress_callback(t, duration):
                cancelled = True
                break
            # adapt dt for the next sub-step
            if pic <= 5:
                dt = min(dt_try * growth, dt_max)
            elif pic > 0.8 * picard_max:
                dt = max(dt_try * 0.5, dt_min)
            else:
                dt = dt_try
        if cancelled:
            break
        # landed exactly on t_target -> save frame with its boundary inflow/outflow
        # at this instant (from the FE reaction — no dual flow-function solve).
        sc = _stored(h) - stored0
        closure = abs(sc - cum_inflow) / max(abs(cum_inflow), 1e-3 * stored_scale, 1e-30)
        bt_s, bv_s, flux_s = resolve_bc(t_target)
        in_s, out_s = _frame_flows(h, bt_s, bv_s, flux_s, active)
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
        "theta": theta,
        "lumped": lumped,
        "duration": duration,
        "unit_weight": gamma_w,
        # An exit face (static seepage face or a falling series-head reservoir) can
        # appear; kr scales velocity on load exactly as in the steady unconfined
        # branch. Recorded so the frame loader derives velocities the same way.
        "unconfined": bool(_may_have_exit),
        # Rapid-drawdown stage times carried through for the meta ledger and the
        # in-memory drawdown staging (§6); None when the model sets no stages.
        "stage_1": tseep_data.get("stage_1"),
        "stage_2": tseep_data.get("stage_2"),
        # The single-run extraction time, carried the same way and for the same
        # reason: an exported march records which instant its stability run reads.
        "stability_time": tseep_data.get("stability_time"),
    }


def _nodal_from_elem(masm, elem_vals, n_nodes, mode='mean'):
    """Scatter a per-element scalar to nodes (mean over incident elements, or 'max'
    for a model code). Used only for the nodal storage-potential ledger, where a
    smooth nodal parameter field is adequate."""
    if elem_vals is None:
        return None
    acc = np.zeros(n_nodes)
    cnt = np.zeros(n_nodes)
    for g in masm['groups']:
        v = elem_vals[g['idx']]
        if mode == 'max':
            np.maximum.at(acc, g['conn'], v[:, None])
        else:
            np.add.at(acc, g['conn'], np.broadcast_to(v[:, None], g['conn'].shape))
        np.add.at(cnt, g['conn'], 1.0)
    if mode == 'max':
        return acc
    return np.where(cnt > 0, acc / np.maximum(cnt, 1.0), 0.0)


def _transient_steady_ic(nodes, elements, element_types, bc_type, bc_values,
                         k1, k2, angle, kr0, h0, vg_a, vg_n, model, flux_nodal,
                         verbose):
    """Initial condition = a steady solve at the t=0 BC configuration (plan IC rule).
    Confined (no exit face) -> one linear solve; unconfined -> the steady
    exit-face Picard solver. Reuses the existing steady solvers unchanged."""
    if np.any(bc_type == 2):
        head, _A, _q, _tf, _efa, _conv, _ce = solve_unsaturated(
            nodes, elements, bc_type, bc_values, kr0=kr0, h0=h0,
            k1_vals=k1, k2_vals=k2, angles=angle, element_types=element_types,
            vg_a=vg_a, vg_n=vg_n, model=model, flux_nodal=flux_nodal)
        return head
    bcs = [(i, bc_values[i]) for i in range(len(bc_type)) if bc_type[i] == 1]
    head, _A, _q, _tf = solve_confined(nodes, elements, bc_type, bcs, k1, k2, angle,
                                       element_types, flux_nodal=flux_nodal)
    return head


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
    # flowrate footer. Both # lines are ordinary CSV comments the readers skip.
    from .units import units_comment_line
    header = units_comment_line(seep_data.get("unit_system"), seep_data.get("time_unit"))
    with open(filename, "w") as f:
        if header:
            f.write(header)
        df.to_csv(f, index=False)
        f.write(f"# Total Flowrate: {solution['flowrate']:.6f}\n")

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
        f.write("# Total Flowrate: not computed (field imported, not solved)\n")


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
        velocity, v_mag, gradient, i_mag, q, phi, flowrate. For a pressure-only file
        (one written by :func:`export_seep_u`, e.g. a field imported from SEEP/W), the
        flow-net columns are filled with NaN: head and u are real, and a flow-net plot
        fails visibly on the NaNs rather than drawing a flat field that is not the
        solution.

    Raises:
        ValueError: if the file's node count does not match the mesh.
    """
    import pandas as pd
    # The trailing '# Total Flowrate: …' line is a comment, skipped on read.
    df = pd.read_csv(filename, comment="#")
    n_nodes = len(seep_data["nodes"])
    if len(df) != n_nodes:
        raise ValueError(
            f"Seepage solution has {len(df)} nodes but the mesh has {n_nodes} — "
            "the saved solution does not match this mesh.")

    flowrate = None
    with open(filename) as f:
        for line in f:
            if line.startswith("# Total Flowrate:"):
                try:
                    flowrate = float(line.split(":", 1)[1])
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
    """Per-element (k1, k2, angle, kr0, h0, vg_a, vg_n, model) arrays from the
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
    model = None if _u is None else _u[mat_ids]
    vg_a = None if _a is None else _a[mat_ids]
    vg_n = None if _n is None else _n[mat_ids]
    return elements, element_types, k1, k2, angle, kr0, h0, vg_a, vg_n, model


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
     vg_a, vg_n, model) = _elem_props(seep_data)
    if unconfined:
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle,
                                    kr0=kr0, h0=h0, element_types=element_types,
                                    vg_a=vg_a, vg_n=vg_n, model=model)
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
            f"on it, or re-march with remarch_for_times().")
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
    # select_transient_frame_u: an automatic stage-1 water load is derived from the
    # pool at the stage-1 time, not at t = 0. (Stage 2 already reads the tseep
    # sheet's own stage_2 time, which is where i2 came from.)
    slope_data["seep_u_time"] = float(transient_solution["times"][i1])
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
                             "there is no transient march to re-run.")
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
                         "run or import the march first.")
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
                  + (" (re-marched)." if remarched else "."))
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
                  + (", re-marched." if remarched else "."))
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
