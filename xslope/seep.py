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

import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from shapely.geometry import LineString, Point


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


def build_seep_data(mesh, slope_data, seep_bc=1):
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
    
    Returns:
        dict: seep_data dictionary with the following structure:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 3 or 4) of element node indices
            - element_types: np.ndarray (n_elements,) indicating 3 for triangles, 4 for quads
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
            - bc_type: np.ndarray (n_nodes,) of boundary condition flags (0=free, 1=fixed head, 2=exit face)
            - bc_values: np.ndarray (n_nodes,) of boundary condition values
            - k1_by_mat: np.ndarray (n_materials,) of major conductivity values
            - k2_by_mat: np.ndarray (n_materials,) of minor conductivity values
            - angle_by_mat: np.ndarray (n_materials,) of angle values (degrees)
            - kr0_by_mat: np.ndarray (n_materials,) of relative conductivity values
            - h0_by_mat: np.ndarray (n_materials,) of suction head values
            - unit_weight: float, unit weight of water
    """
    
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
    material_names = []
    
    for i, material in enumerate(materials):
        k1_by_mat[i] = material.get("k1", 1.0)
        k2_by_mat[i] = material.get("k2", 1.0)
        angle_by_mat[i] = material.get("alpha", 0.0)
        kr0_by_mat[i] = material.get("kr0", 0.001)
        h0_by_mat[i] = material.get("h0", -1.0)
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
    
    print(f"Mesh tolerance for boundary conditions: {tolerance:.6f}")
    
    # Process specified head boundary conditions
    # Vectorized: compute distance from all nodes to each BC line at once
    specified_heads = seepage_bc.get("specified_heads", [])
    for bc in specified_heads:
        head_value = bc["head"]
        coords = bc["coords"]

        if len(coords) < 2:
            continue

        # Compute min distance from each node to any segment of the BC line
        seg_coords = np.array(coords)
        dists = _min_distance_to_polyline(nodes, seg_coords)
        mask = dists <= tolerance
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
    
    # Check for missing unsaturated parameters when exit face BCs are present
    has_exit_face = np.any(bc_type == 2)
    missing_unsat_params = False
    if has_exit_face:
        bad_mats = []
        for i, material in enumerate(materials):
            mat_kr0 = material.get("kr0", 0)
            mat_h0 = material.get("h0", 0)
            kr0_missing = mat_kr0 is None or (isinstance(mat_kr0, (int, float)) and (mat_kr0 == 0 or np.isnan(mat_kr0)))
            h0_missing = mat_h0 is None or (isinstance(mat_h0, (int, float)) and (mat_h0 == 0 or np.isnan(mat_h0)))
            if kr0_missing or h0_missing:
                name = material.get("name", f"Material {i+1}")
                bad_mats.append(f"  Material {i+1} ({name}): kr0={mat_kr0}, h0={mat_h0}")
        if bad_mats:
            missing_unsat_params = True
            print("\n" + "="*70)
            print("WARNING: Exit face boundary condition detected but the following")
            print("materials are missing unsaturated parameters (kr0 and/or h0):")
            for line in bad_mats:
                print(line)
            print("\nThe unsaturated seepage solver requires valid kr0 (>0) and h0 (<0)")
            print("values for all materials. Please set these in the input file.")
            print("="*70 + "\n")

    # Get unit weight of water
    unit_weight = slope_data.get("gamma_water", 9.81)

    # Construct seep_data dictionary
    seep_data = {
        "nodes": nodes,
        "elements": elements,
        "element_types": element_types,
        "element_materials": element_materials,
        "bc_type": bc_type,
        "bc_values": bc_values,
        "k1_by_mat": k1_by_mat,
        "k2_by_mat": k2_by_mat,
        "angle_by_mat": angle_by_mat,
        "kr0_by_mat": kr0_by_mat,
        "h0_by_mat": h0_by_mat,
        "material_names": material_names,
        "unit_weight": unit_weight,
        "missing_unsat_params": missing_unsat_params
    }

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
        "unit_weight": unit_weight
    }


def solve_confined(nodes, elements, bc_type, dirichlet_bcs, k1_vals, k2_vals, angles=None, element_types=None):
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
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

    head = spsolve(A, b)
    q = A_full @ head

    # Sum only positive q at specified-head nodes (inflow)
    bc_type = np.asarray(bc_type)
    total_flow = float(np.sum(q[(bc_type == 1) & (q > 0)]))

    return head, A, q, total_flow


def solve_unsaturated(nodes, elements, bc_type, bc_values, kr0=0.001, h0=-1.0,
                      k1_vals=1.0, k2_vals=1.0, angles=0.0,
                      max_iter=400, tol=1e-6, element_types=None,
                      closure_tol=1e-3):
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

    # Initialize heads
    fixed_heads = bc_values[bc_type == 1]
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
    _exit_is_corner = np.zeros(n_nodes, dtype=bool)
    _exit_linear_corners = np.zeros(n_nodes, dtype=bool)
    _exit_quadratic_edges = []  # [(corner1, midside, corner2), ...]
    _seen_edges = set()
    for _idx, _en in enumerate(elements):
        _et = element_types[_idx]
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
            _gnodes = [_en[i] for i in _le]
            _c1, _c2 = _gnodes[0], _gnodes[-1]
            _ekey = tuple(sorted((_c1, _c2)))
            if _ekey in _seen_edges:
                continue
            if bc_type[_c1] == 2 and bc_type[_c2] == 2:
                _seen_edges.add(_ekey)
                _exit_is_corner[_c1] = True
                _exit_is_corner[_c2] = True
                if len(_gnodes) == 3:  # has midside node
                    _exit_quadratic_edges.append((_c1, _gnodes[1], _c2))
                else:
                    _exit_linear_corners[_c1] = True
                    _exit_linear_corners[_c2] = True

    # Store previous iteration values
    h_last = h.copy()

    # Get material properties per element
    if np.isscalar(kr0):
        kr0 = np.full(len(elements), kr0)
    if np.isscalar(h0):
        h0 = np.full(len(elements), h0)

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
    data = _assembly_data(asm, p_nodes=h - y, kr0=kr0, h0=h0, mode='head')
    _prev_active = exit_face_active.copy()
    _n_stable = 0   # consecutive iterations with an unchanged exit-face set

    for iteration in range(1, max_iter + 1):

        # Apply boundary conditions: fixed heads plus the active exit-face set
        dir_mask = (bc_type == 1) | ((bc_type == 2) & exit_face_active)
        dir_values = np.where(bc_type == 1, bc_values, y)
        A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

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

        # Compute flows at all nodes (not used for closure, but for exit face logic)
        q = _coo_matvec(asm, data, h_new)

        # Update exit face boundary conditions with hysteresis
        n_active_before = np.sum(exit_face_active)
        hyst = 0.001 * (ymax - ymin)  # Hysteresis threshold

        # Check exit face activation at corner nodes using per-node q
        # (matching SEEP2D). Quadratic exit-face sides are then updated from
        # these corner candidates plus the midside candidate so the whole side
        # remains edge-consistent.
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
        data_chk = _assembly_data(asm, p_nodes=h_solved - y, kr0=kr0, h0=h0, mode='head')
        q_chk = _coo_matvec(asm, data_chk, h_solved)
        free_mask = ~((bc_type == 1) | ((bc_type == 2) & exit_face_active))
        inflow_pos = float(np.sum(q_chk[(bc_type == 1) & (q_chk > 0)]))
        rel_closure = (float(np.sum(np.abs(q_chk[free_mask]))) / inflow_pos
                       if inflow_pos > 1e-30 else 0.0)
        set_stable = bool(np.array_equal(exit_face_active, _prev_active))
        _prev_active = exit_face_active.copy()
        _n_stable = _n_stable + 1 if set_stable else 0

        # Matrix for the next iteration, from the relaxed head field
        data = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head')

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
    dir_values = np.where(bc_type == 1, bc_values, y)
    A_final, b_final = _dirichlet_system(asm, data, dir_mask, dir_values)
    h_new = spsolve(A_final, b_final)
    q_final = _coo_matvec(asm, data, h_new)

    # Flowrate: sum of positive q at specified-head nodes
    total_inflow = float(np.sum(q_final[(bc_type == 1) & (q_final > 0)]))

    # Closure report: kr-consistent imbalance at the converged state
    data_chk = _assembly_data(asm, p_nodes=h_new - y, kr0=kr0, h0=h0, mode='head')
    q_chk = _coo_matvec(asm, data_chk, h_new)
    net_inflow = float(np.sum(q_chk[bc_type == 1]))
    net_outflow = -float(np.sum(q_chk[bc_type == 2]))
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
    lin = kr0 + (1.0 - kr0) * (p - h0) / (-h0)
    return np.where(p >= 0.0, 1.0, np.where(p > h0, lin, kr0))


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


def _assembly_data(asm, p_nodes=None, kr0=None, h0=None, mode='head'):
    """Concatenated COO data for the global matrix.

    p_nodes=None -> saturated assembly. Otherwise each element's ke is scaled
    by factor(kr_avg) with kr averaged at the type-specific sampling points
    (mode='head': kr_avg; mode='stream': 1/kr_avg, guarded)."""
    parts = []
    for g in asm['groups']:
        if p_nodes is None:
            parts.append(g['ke'].ravel())
            continue
        p_gp = p_nodes[g['conn']] @ g['N_kr'].T              # (n_e, n_pts)
        kr = kr_frontal_vec(p_gp, kr0[g['idx']][:, None], h0[g['idx']][:, None])
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


def _dirichlet_system(asm, data, dir_mask, dir_values):
    """BC-applied system (A, b): Dirichlet rows replaced by identity.

    Equivalent to the previous LIL row-zeroing, built directly from the COO
    arrays (entries in Dirichlet rows are dropped, then unit diagonals are
    appended)."""
    n = asm['n_nodes']
    keep = ~dir_mask[asm['rows']]
    di = np.where(dir_mask)[0]
    rows = np.concatenate([asm['rows'][keep], di])
    cols = np.concatenate([asm['cols'][keep], di])
    vals = np.concatenate([data[keep], np.ones(len(di))])
    A = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    b = np.zeros(n)
    b[di] = dir_values[di]
    return A, b


def _coo_matvec(asm, data, x):
    """y = A @ x directly from the COO arrays (no tocsr / duplicate handling)."""
    return np.bincount(asm['rows'], weights=data * x[asm['cols']],
                       minlength=asm['n_nodes'])


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
                                           bc_type=None, exit_face_active=None):
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
            kr_elem = kr_frontal(p_elem, kr_e0, h0_e0)

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
                    kr_gp = kr_frontal(p_gp, kr_e0, h0_e0) if kr0 is not None else 1.0
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

def solve_flow_function_unsaturated(nodes, elements, head, k1_vals, k2_vals, angles, kr0, h0, dirichlet_nodes, element_types=None):
    """
    Solves the stream function (flow function) Phi for unsaturated flow.

    For anisotropic permeability, the stream function equation uses K/det(K)
    (not K^(-1)) in the stiffness matrix assembly. The relative permeability
    kr is also included in the formulation.

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

    # Stream-function coefficient matrix is K/det(K); each element is scaled
    # by 1/kr_avg (mode='stream'), matching the per-element *_kr functions.
    asm = _build_assembly(nodes, elements, element_types, k1_vals, k2_vals,
                          angles, flow=True)
    data = _assembly_data(asm, p_nodes=p_nodes, kr0=kr0, h0=h0, mode='stream')

    dir_mask = np.zeros(n_nodes, dtype=bool)
    dir_values = np.zeros(n_nodes)
    for node, phi_value in dirichlet_nodes:
        dir_mask[node] = True
        dir_values[node] = phi_value
    A, b = _dirichlet_system(asm, data, dir_mask, dir_values)

    phi = spsolve(A, b)
    return phi


def compute_velocity(nodes, elements, head, k1_vals, k2_vals, angles, kr0=None, h0=None, element_types=None):
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
                kr_e = kr_frontal_vec(p_all[conn].mean(axis=1), kr0[idx], h0[idx])
            else:
                kr_e = np.ones(len(idx))
            v_e = -kr_e[:, None] * np.einsum('exy,ey->ex', K, grad_h)
            v_e[~ok] = 0.0
            np.add.at(velocity, conn.ravel(), np.repeat(v_e, 3, axis=0))
            np.add.at(count, conn.ravel(), np.repeat(ok.astype(float), 3))

        elif et == 4:
            if use_kr:
                kr_e = kr_frontal_vec(p_all[conn].mean(axis=1), kr0[idx], h0[idx])
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
                kr_e = kr_frontal_vec(p_c, kr0[idx], h0[idx])
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


def run_seepage_analysis(seep_data, tol=1e-6, closure_tol=1e-3):
    """
    Standalone function to run seep analysis.
    
    Args:
        seep_data: Dictionary containing all the seep data
        tol: relative head-change tolerance (scaled by domain height)
        closure_tol: relative flow-closure tolerance for unconfined problems —
            iteration continues until |net inflow - net outflow| / inflow is
            below this, so the reported flowrate balances regardless of how
            the head tolerance maps to mass balance on a given problem
    
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
    """
    start_time = time.time()

    # Check for missing unsaturated parameters
    if seep_data.get("missing_unsat_params", False):
        print("\n" + "="*70)
        print("ERROR: Cannot run seepage analysis.")
        print("One or more materials are missing unsaturated parameters (kr0, h0)")
        print("which are required for unconfined seepage with an exit face BC.")
        print("Please set valid kr0 (>0) and h0 (<0) values in the input file.")
        print("="*70 + "\n")
        return None

    # Extract data from seep_data
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    bc_type = seep_data["bc_type"]
    bc_values = seep_data["bc_values"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    k1_by_mat = seep_data["k1_by_mat"]
    k2_by_mat = seep_data["k2_by_mat"]
    angle_by_mat = seep_data["angle_by_mat"]
    kr0_by_mat = seep_data["kr0_by_mat"]
    h0_by_mat = seep_data["h0_by_mat"]
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
        # Get kr0 and h0 values per element based on material
        kr0_per_element = kr0_by_mat[mat_ids]
        h0_per_element = h0_by_mat[mat_ids]

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
            closure_tol=closure_tol
        )
        # Compute phi BCs from element-level boundary flux
        dirichlet_phi_bcs = create_flow_potential_bc_from_elements(
            nodes, elements, element_types, head, k1, k2, angle,
            kr0=kr0_per_element, h0=h0_per_element, total_flow=total_flow,
            bc_type=bc_type, exit_face_active=exit_face_active)
        phi = solve_flow_function_unsaturated(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, dirichlet_phi_bcs, element_types)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, element_types)
    else:
        # Confined analysis is a single direct linear solve — always "converged".
        converged, closure_error = True, 0.0
        head, A, q, total_flow = solve_confined(nodes, elements, bc_type, bcs, k1, k2, angle, element_types)
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
        "converged": converged,
        "closure_error": closure_error,
    }

    if not converged:
        print("WARNING: seepage solution did not converge — flowrate is unreliable "
              "(solution['converged'] is False).")

    elapsed = time.time() - start_time
    print(f"Seepage analysis completed in {elapsed:.2f} seconds.")

    return solution

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
    # Write to file, then append flowrate as comment
    with open(filename, "w") as f:
        df.to_csv(f, index=False)
        f.write(f"# Total Flowrate: {solution['flowrate']:.6f}\n")

    print(f"Exported solution to {filename}")

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
    
    # Convert numpy arrays to lists for JSON serialization
    seep_data_json = {}
    for key, value in seep_data.items():
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
