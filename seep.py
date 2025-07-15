import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from shapely.geometry import LineString, Point


def build_seep_data(mesh, data):
    """
    Build a seep_data dictionary from a mesh and data dictionary.
    
    This function takes a mesh dictionary (from build_mesh_from_polygons) and a data dictionary
    (from load_globals) and constructs a seep_data dictionary suitable for seepage analysis.
    
    The function:
    1. Extracts mesh information (nodes, elements, element types, element materials)
    2. Builds material property arrays (k1, k2, alpha, kr0, h0) from the materials table
    3. Constructs boundary conditions by finding nodes that intersect with specified head
       and seepage face lines from the data dictionary
    
    Parameters:
        mesh (dict): Mesh dictionary from build_mesh_from_polygons containing:
            - nodes: np.ndarray (n_nodes, 2) of node coordinates
            - elements: np.ndarray (n_elements, 3 or 4) of element node indices
            - element_types: np.ndarray (n_elements,) indicating 3 for triangles, 4 for quads
            - element_materials: np.ndarray (n_elements,) of material IDs (1-based)
        data (dict): Data dictionary from load_globals containing:
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
    materials = data["materials"]
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
    seepage_bc = data.get("seepage_bc", {})
    
    # Calculate appropriate tolerance based on mesh size
    # Use a fraction of the typical element size
    x_range = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
    y_range = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    typical_element_size = min(x_range, y_range) / np.sqrt(len(nodes))  # Approximate element size
    tolerance = typical_element_size * 0.1  # 10% of typical element size
    
    print(f"Mesh tolerance for boundary conditions: {tolerance:.6f}")
    
    # Process specified head boundary conditions
    specified_heads = seepage_bc.get("specified_heads", [])
    for bc in specified_heads:
        head_value = bc["head"]
        coords = bc["coords"]
        
        if len(coords) < 2:
            continue
            
        # Create LineString from boundary condition coordinates
        bc_line = LineString(coords)
        
        # Find nodes that are close to this line (within tolerance)
        for i, node_coord in enumerate(nodes):
            node_point = Point(node_coord)
            
            # Check if node is on or very close to the boundary condition line
            if bc_line.distance(node_point) <= tolerance:
                bc_type[i] = 1  # Fixed head
                bc_values[i] = head_value
    
    # Process seepage face (exit face) boundary conditions
    exit_face_coords = seepage_bc.get("exit_face", [])
    if len(exit_face_coords) >= 2:
        # Create LineString from exit face coordinates
        exit_face_line = LineString(exit_face_coords)
        
        # Find nodes that are close to this line
        for i, node_coord in enumerate(nodes):
            node_point = Point(node_coord)
            
            # Check if node is on or very close to the exit face line
            if exit_face_line.distance(node_point) <= tolerance:
                bc_type[i] = 2  # Exit face
                bc_values[i] = node_coord[1]  # Use node's y-coordinate as elevation
    
    # Get unit weight of water
    unit_weight = data.get("gamma_water", 9.81)
    
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
        "unit_weight": unit_weight
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

    return {
        "nodes": np.array(coords),
        "bc_type": np.array(bc_type, dtype=int),
        "bc_values": np.array(bc_values),
        "elements": np.array(elements, dtype=int),  # Already 0-based
        "element_types": np.array(element_types, dtype=int),
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
    FEM solver for confined seepage with anisotropic conductivity.
    Supports both triangular and quadrilateral elements.
    Parameters:
        nodes : (n_nodes, 2) array of node coordinates
        elements : (n_elements, 3 or 4) triangle or quad node indices
        bc_type : (n_nodes,) array of boundary condition flags
        dirichlet_bcs : list of (node_id, head_value)
        k1_vals : (n_elements,) or scalar, major axis conductivity
        k2_vals : (n_elements,) or scalar, minor axis conductivity
        angles : (n_elements,) or scalar, angle in degrees (from x-axis)
        element_types : (n_elements,) array indicating 3 for triangles, 4 for quads
    Returns:
        head : (n_nodes,) array of nodal heads
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    A = lil_matrix((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: use first 3 nodes (4th node is repeated)
            i, j, k = element_nodes[:3]
            xi, yi = nodes[i]
            xj, yj = nodes[j]
            xk, yk = nodes[k]

            area = 0.5 * np.linalg.det([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
            if area <= 0:
                continue

            beta = np.array([yj - yk, yk - yi, yi - yj])
            gamma = np.array([xk - xj, xi - xk, xj - xi])
            grad = np.array([beta, gamma]) / (2 * area)

            # Get anisotropic conductivity
            k1 = k1_vals[idx]
            k2 = k2_vals[idx]
            theta = angles[idx]

            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R

            ke = area * grad.T @ Kmat @ grad

            for a in range(3):
                for b_ in range(3):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]
        elif element_type == 4:
            # Quadrilateral: use all 4 nodes
            i, j, k, l = element_nodes
            nodes_elem = nodes[[i, j, k, l], :]
            k1 = k1_vals[idx]
            k2 = k2_vals[idx]
            theta = angles[idx]
            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R
            ke = quad4_stiffness_matrix(nodes_elem, Kmat)
            for a in range(4):
                for b_ in range(4):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]

    A_full = A.copy()  # Keep original matrix for computing q

    for node, value in dirichlet_bcs:
        A[node, :] = 0
        A[node, node] = 1
        b[node] = value

    head = spsolve(A.tocsr(), b)
    q = A_full.tocsr() @ head

    total_flow = 0.0

    for node_idx in range(len(bc_type)):
        if q[node_idx] > 0:  # Positive flow
            total_flow += q[node_idx]

    return head, A, q, total_flow


def solve_unsaturated(nodes, elements, bc_type, bc_values, kr0=0.001, h0=-1.0,
                      k1_vals=1.0, k2_vals=1.0, angles=0.0,
                      max_iter=200, tol=1e-4, element_types=None):
    """
    Iterative FEM solver for unconfined flow using linear kr frontal function.
    Supports both triangular and quadrilateral elements.
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    y = nodes[:, 1]

    # Initialize heads
    h = np.zeros(n_nodes)
    for node_idx in range(n_nodes):
        if bc_type[node_idx] == 1:
            h[node_idx] = bc_values[node_idx]
        elif bc_type[node_idx] == 2:
            h[node_idx] = y[node_idx]
        else:
            fixed_heads = bc_values[bc_type == 1]
            h[node_idx] = np.mean(fixed_heads) if len(fixed_heads) > 0 else np.mean(y)

    # Track which exit face nodes are active (saturated)
    exit_face_active = np.ones(n_nodes, dtype=bool)
    exit_face_active[bc_type != 2] = False

    # Store previous iteration values
    h_last = h.copy()

    # Get material properties per element
    if np.isscalar(kr0):
        kr0 = np.full(len(elements), kr0)
    if np.isscalar(h0):
        h0 = np.full(len(elements), h0)

    # Set convergence tolerance based on domain height
    ymin, ymax = np.min(y), np.max(y)
    eps = (ymax - ymin) * 0.0001

    print("Starting unsaturated flow iteration...")
    print(f"Convergence tolerance: {eps:.6e}")

    # Track convergence history
    residuals = []
    relax = 1.0  # Initial relaxation factor
    prev_residual = float('inf')

    for iteration in range(1, max_iter + 1):
        # Reset diagnostics for this iteration
        kr_diagnostics = []

        # Build global stiffness matrix
        A = lil_matrix((n_nodes, n_nodes))
        b = np.zeros(n_nodes)

        # Compute pressure head at nodes
        p_nodes = h - y

        # Element assembly with element-wise kr computation
        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            
            if element_type == 3:
                # Triangle: use first 3 nodes (4th node is repeated)
                i, j, k = element_nodes[:3]
                xi, yi = nodes[i]
                xj, yj = nodes[j]
                xk, yk = nodes[k]

                # Element area
                area = 0.5 * abs((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi))
                if area <= 0:
                    continue

                # Shape function derivatives
                beta = np.array([yj - yk, yk - yi, yi - yj])
                gamma = np.array([xk - xj, xi - xk, xj - xi])
                grad = np.array([beta, gamma]) / (2 * area)

                # Get material properties for this element
                k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
                k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
                theta = angles[idx] if hasattr(angles, '__len__') else angles

                # Anisotropic conductivity matrix
                theta_rad = np.radians(theta)
                c, s = np.cos(theta_rad), np.sin(theta_rad)
                R = np.array([[c, s], [-s, c]])
                Kmat = R.T @ np.diag([k1, k2]) @ R

                # Compute element pressure (centroid)
                p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k]) / 3.0

                # Get kr for this element based on its material properties
                kr_elem = kr_frontal(p_elem, kr0[idx], h0[idx])

                # Element stiffness matrix with kr
                ke = kr_elem * area * grad.T @ Kmat @ grad

                # Inside the element assembly loop, after kr_elem is computed:
                kr_diagnostics.append({
                    'element': idx,
                    'p_elem': p_elem,
                    'kr_elem': kr_elem,
                    'y_centroid': (yi + yj + yk) / 3.0,
                    'h_centroid': (h[i] + h[j] + h[k]) / 3.0
                })

                # Assembly
                for row in range(3):
                    for col in range(3):
                        A[element_nodes[row], element_nodes[col]] += ke[row, col]
            elif element_type == 4:
                # Quadrilateral: use all 4 nodes
                i, j, k, l = element_nodes
                nodes_elem = nodes[[i, j, k, l], :]
                k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
                k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
                theta = angles[idx] if hasattr(angles, '__len__') else angles
                theta_rad = np.radians(theta)
                c, s = np.cos(theta_rad), np.sin(theta_rad)
                R = np.array([[c, s], [-s, c]])
                Kmat = R.T @ np.diag([k1, k2]) @ R
                # Compute element pressure (centroid)
                p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k] + p_nodes[l]) / 4.0
                kr_elem = kr_frontal(p_elem, kr0[idx], h0[idx])
                ke = kr_elem * quad4_stiffness_matrix(nodes_elem, Kmat)
                for a in range(4):
                    for b_ in range(4):
                        A[element_nodes[a], element_nodes[b_]] += ke[a, b_]

        # Store unmodified matrix for flow computation
        A_full = A.tocsr()

        # Apply boundary conditions
        for node_idx in range(n_nodes):
            if bc_type[node_idx] == 1:
                A[node_idx, :] = 0
                A[node_idx, node_idx] = 1
                b[node_idx] = bc_values[node_idx]
            elif bc_type[node_idx] == 2 and exit_face_active[node_idx]:
                A[node_idx, :] = 0
                A[node_idx, node_idx] = 1
                b[node_idx] = y[node_idx]

        # Convert to CSR and solve
        A_csr = A.tocsr()
        h_new = spsolve(A_csr, b)

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
        h_new = relax * h_new + (1 - relax) * h_last

        # Compute flows at all nodes (not used for closure, but for exit face logic)
        q = A_full @ h_new

        # Update exit face boundary conditions with hysteresis
        n_active_before = np.sum(exit_face_active)
        hyst = 0.001 * (ymax - ymin)  # Hysteresis threshold

        for node_idx in range(n_nodes):
            if bc_type[node_idx] == 2:
                if exit_face_active[node_idx]:
                    # Check if node should become inactive
                    if h_new[node_idx] < y[node_idx] - hyst or q[node_idx] > 0:
                        exit_face_active[node_idx] = False
                else:
                    # Check if node should become active again
                    if h_new[node_idx] >= y[node_idx] + hyst and q[node_idx] <= 0:
                        exit_face_active[node_idx] = True
                        h_new[node_idx] = y[node_idx]  # Reset to elevation

        n_active_after = np.sum(exit_face_active)

        # Compute relative residual
        residual = np.max(np.abs(h_new - h)) / (np.max(np.abs(h)) + 1e-10)
        residuals.append(residual)

        # Print detailed iteration info
        if iteration <= 3 or iteration % 5 == 0 or n_active_before != n_active_after:
            print(f"Iteration {iteration}: residual = {residual:.6e}, relax = {relax:.3f}, {n_active_after}/{np.sum(bc_type == 2)} exit face active")
            #print(f"  BCs: {np.sum(bc_type == 1)} fixed head, {n_active_after}/{np.sum(bc_type == 2)} exit face active")

        # Check convergence
        if residual < eps:
            print(f"Converged in {iteration} iterations")
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


    q_final = q

    # Flow potential closure check - FORTRAN-style
    total_inflow = 0.0
    total_outflow = 0.0
    
    for node_idx in range(n_nodes):
        if bc_type[node_idx] == 1:  # Fixed head boundary
            if q_final[node_idx] > 0:
                total_inflow += q_final[node_idx]
            elif q_final[node_idx] < 0:
                total_outflow -= q_final[node_idx]
        elif bc_type[node_idx] == 2 and exit_face_active[node_idx]:  # Active exit face
            if q_final[node_idx] < 0:
                total_outflow -= q_final[node_idx]

    closure_error = abs(total_inflow - total_outflow)
    print(f"Flow potential closure check: error = {closure_error:.6e}")
    print(f"Total inflow: {total_inflow:.6e}")
    print(f"Total outflow: {total_outflow:.6e}")

    if closure_error > 0.01 * max(abs(total_inflow), abs(total_outflow)):
        print(f"Warning: Large flow potential closure error = {closure_error:.6e}")
        print("This may indicate:")
        print("  - Non-conservative flow field")
        print("  - Incorrect boundary identification")
        print("  - Numerical issues in the flow solution")


    return h, A, q_final, total_inflow

def kr_frontal(p, kr0, h0):
    """
    Fortran-compatible relative permeability function (front model).
    This matches the fkrelf function in the Fortran code exactly.
    """
    if p >= 0.0:
        return 1.0
    elif p > h0:
        return kr0 + (1.0 - kr0) * p / h0
    else:
        return kr0


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

def create_flow_potential_bc(nodes, elements, q, debug=False, element_types=None):
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

    # Step 1: Build edge dictionary and count how many times each edge appears
    edge_counts = defaultdict(list)
    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: 3 edges (use first 3 nodes, 4th is repeated)
            i, j, k = element_nodes[:3]
            edges = [(i, j), (j, k), (k, i)]
        elif element_type == 4:
            # Quadrilateral: 4 edges (use all 4 nodes)
            i, j, k, l = element_nodes
            edges = [(i, j), (j, k), (k, l), (l, i)]
        else:
            continue  # Skip unknown element types
            
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

    # Strategy 1: Look for transition from significant flow to near-zero flow
    for i in range(n):
        current_q = abs(q[ordered_nodes[i]])
        next_q = abs(q[ordered_nodes[(i + 1) % n]])

        if current_q > q_threshold and next_q <= q_threshold:
            start_idx = (i + 1) % n
            if debug:
                print(f"Found transition at node {ordered_nodes[i]} -> {ordered_nodes[start_idx]}")
            break

    # Strategy 2: If no clear transition, find the node with minimum |q|
    if start_idx is None:
        min_q_idx = min(range(n), key=lambda i: abs(q[ordered_nodes[i]]))
        start_idx = min_q_idx
        if debug:
            print(
                f"No clear transition found, starting at minimum |q| node {ordered_nodes[start_idx]} (|q|={abs(q[ordered_nodes[start_idx]]):.3e})")

    # Strategy 3: If all flows are significant, look for flow direction changes
    if start_idx is None or abs(q[ordered_nodes[start_idx]]) > q_threshold:
        # Look for where flow changes sign (inflow vs outflow)
        for i in range(n):
            current_q = q[ordered_nodes[i]]
            next_q = q[ordered_nodes[(i + 1) % n]]

            if current_q * next_q < 0:  # Sign change
                start_idx = i if abs(current_q) < abs(next_q) else (i + 1) % n
                if debug:
                    print(f"Found sign change, starting at node {ordered_nodes[start_idx]}")
                break

    # Strategy 4: Default fallback - start at node 0
    if start_idx is None:
        start_idx = 0
        if debug:
            print(f"Using fallback: starting at first boundary node {ordered_nodes[start_idx]}")

    # Step 6: Assign φ = 0 to the starting node, then accumulate q
    phi = {}
    phi_val = 0.0

    if debug:
        print(f"Starting flow potential calculation at node {ordered_nodes[start_idx]}")

    for i in range(n):
        idx = (start_idx + i) % n
        node = ordered_nodes[idx]
        phi[node] = phi_val
        phi_val += q[node]

        if debug and (i < 5 or i >= n - 5):  # Print first and last few for debugging
            print(f"  Node {node}: φ = {phi[node]:.6f}, q = {q[node]:.6f}")

    # Check closure - should be close to zero for a proper flow field
    closure_error = phi_val - q[ordered_nodes[start_idx]]

    if debug or abs(closure_error) > 1e-3:
        print(f"Flow potential closure check: error = {closure_error:.6e}")

        if abs(closure_error) > 1e-3:
            print(f"Warning: Large flow potential closure error = {closure_error:.6e}")
            print("This may indicate:")
            print("  - Non-conservative flow field")
            print("  - Incorrect boundary identification")
            print("  - Numerical issues in the flow solution")

    if debug:
        print("✓ Flow potential BC creation succeeded")

    return list(phi.items())

def solve_flow_function_confined(nodes, elements, k1_vals, k2_vals, angles, dirichlet_nodes, element_types=None):
    """
    Solves Laplace equation for flow function Phi on the same mesh,
    assigning Dirichlet values along no-flow boundaries.
    Assembles the element matrix using the inverse of Kmat for each element.
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
    A = lil_matrix((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: use first 3 nodes (4th node is repeated)
            i, j, k = element_nodes[:3]
            xi, yi = nodes[i]
            xj, yj = nodes[j]
            xk, yk = nodes[k]

            area = 0.5 * np.linalg.det([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
            if area <= 0:
                continue

            beta = np.array([yj - yk, yk - yi, yi - yj])
            gamma = np.array([xk - xj, xi - xk, xj - xi])
            grad = np.array([beta, gamma]) / (2 * area)

            # Get anisotropic conductivity for this element
            k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
            k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
            theta = angles[idx] if hasattr(angles, '__len__') else angles

            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R

            # Assemble using the inverse of Kmat
            ke = area * grad.T @ np.linalg.inv(Kmat) @ grad

            for a in range(3):
                for b_ in range(3):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]
        elif element_type == 4:
            # Quadrilateral: use all 4 nodes
            i, j, k, l = element_nodes
            nodes_elem = nodes[[i, j, k, l], :]
            k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
            k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
            theta = angles[idx] if hasattr(angles, '__len__') else angles
            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R
            ke = quad4_stiffness_matrix(nodes_elem, np.linalg.inv(Kmat))
            for a in range(4):
                for b_ in range(4):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]

    for node, phi_value in dirichlet_nodes:
        A[node, :] = 0
        A[node, node] = 1
        b[node] = phi_value

    phi = spsolve(A.tocsr(), b)
    return phi

def solve_flow_function_unsaturated(nodes, elements, head, k1_vals, k2_vals, angles, kr0, h0, dirichlet_nodes, element_types=None):
    """
    Solves the flow function Phi using the correct ke for unsaturated flow.
    For flowlines, assemble the element matrix using the inverse of kr_elem and Kmat, matching the FORTRAN logic.
    Supports both triangular and quadrilateral elements.
    """

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    n_nodes = nodes.shape[0]
    A = lil_matrix((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    y = nodes[:, 1]
    p_nodes = head - y

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: use first 3 nodes (4th node is repeated)
            i, j, k = element_nodes[:3]
            xi, yi = nodes[i]
            xj, yj = nodes[j]
            xk, yk = nodes[k]

            area = 0.5 * abs((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi))
            if area <= 0:
                continue

            beta = np.array([yj - yk, yk - yi, yi - yj])
            gamma = np.array([xk - xj, xi - xk, xj - xi])
            grad = np.array([beta, gamma]) / (2 * area)  # grad is (2,3)

            # Get material properties for this element
            k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
            k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
            theta = angles[idx] if hasattr(angles, '__len__') else angles

            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R  # Kmat is (2,2)

            # Compute element pressure (centroid)
            p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k]) / 3.0
            kr_elem = kr_frontal(p_elem, kr0[idx], h0[idx])

            # Assemble using the inverse of kr_elem and Kmat
            # If kr_elem is very small, avoid division by zero
            if kr_elem > 1e-12:
                ke = (1.0 / kr_elem) * area * grad.T @ np.linalg.inv(Kmat) @ grad
            else:
                ke = 1e12 * area * grad.T @ np.linalg.inv(Kmat) @ grad  # Large value for near-zero kr

            for a in range(3):
                for b_ in range(3):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]
        elif element_type == 4:
            # Quadrilateral: use all 4 nodes
            i, j, k, l = element_nodes
            nodes_elem = nodes[[i, j, k, l], :]
            k1 = k1_vals[idx] if hasattr(k1_vals, '__len__') else k1_vals
            k2 = k2_vals[idx] if hasattr(k2_vals, '__len__') else k2_vals
            theta = angles[idx] if hasattr(angles, '__len__') else angles
            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            Kmat = R.T @ np.diag([k1, k2]) @ R  # Kmat is (2,2)
            # Get kr for this element based on its material properties (use centroid)
            p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k] + p_nodes[l]) / 4.0
            kr_elem = kr_frontal(p_elem, kr0[idx], h0[idx])
            # Assemble using the inverse of kr_elem and Kmat
            if kr_elem > 1e-12:
                ke = (1.0 / kr_elem) * quad4_stiffness_matrix(nodes_elem, np.linalg.inv(Kmat))
            else:
                ke = 1e12 * quad4_stiffness_matrix(nodes_elem, np.linalg.inv(Kmat))
            for a in range(4):
                for b_ in range(4):
                    A[element_nodes[a], element_nodes[b_]] += ke[a, b_]

    for node, phi_value in dirichlet_nodes:
        A[node, :] = 0
        A[node, node] = 1
        b[node] = phi_value

    phi = spsolve(A.tocsr(), b)
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

    n_nodes = nodes.shape[0]
    velocity = np.zeros((n_nodes, 2))
    count = np.zeros(n_nodes)

    scalar_k = np.isscalar(k1_vals)
    scalar_kr = kr0 is not None and np.isscalar(kr0)

    y = nodes[:, 1]
    p_nodes = head - y

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: use first 3 nodes (4th node is repeated)
            i, j, k = element_nodes[:3]
            xi, yi = nodes[i]
            xj, yj = nodes[j]
            xk, yk = nodes[k]

            area = 0.5 * np.linalg.det([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
            if area <= 0:
                continue

            beta = np.array([yj - yk, yk - yi, yi - yj])
            gamma = np.array([xk - xj, xi - xk, xj - xi])
            grad = np.array([beta, gamma]) / (2 * area)

            h_vals = head[[i, j, k]]
            grad_h = grad @ h_vals

            if scalar_k:
                k1 = k1_vals
                k2 = k2_vals
                theta = angles
            else:
                k1 = k1_vals[idx]
                k2 = k2_vals[idx]
                theta = angles[idx]

            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            K = R.T @ np.diag([k1, k2]) @ R

            # Compute kr_elem if kr0 and h0 are provided
            if kr0 is not None and h0 is not None:
                p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k]) / 3.0
                kr_elem = kr_frontal(p_elem, kr0[idx] if not scalar_kr else kr0, h0[idx] if not scalar_kr else h0)
            else:
                kr_elem = 1.0

            v_elem = -kr_elem * K @ grad_h

            for node in element_nodes[:3]:  # Only use first 3 nodes for triangles
                velocity[node] += v_elem
                count[node] += 1
        elif element_type == 4:
            # Quadrilateral: use all 4 nodes
            i, j, k, l = element_nodes
            nodes_elem = nodes[[i, j, k, l], :]
            h_elem = head[[i, j, k, l]]
            if scalar_k:
                k1 = k1_vals
                k2 = k2_vals
                theta = angles
            else:
                k1 = k1_vals[idx]
                k2 = k2_vals[idx]
                theta = angles[idx]
            theta_rad = np.radians(theta)
            c, s = np.cos(theta_rad), np.sin(theta_rad)
            R = np.array([[c, s], [-s, c]])
            K = R.T @ np.diag([k1, k2]) @ R
            if kr0 is not None and h0 is not None:
                p_elem = np.mean(p_nodes[[i, j, k, l]])
                kr_elem = kr_frontal(p_elem, kr0[idx] if not scalar_kr else kr0, h0[idx] if not scalar_kr else h0)
            else:
                kr_elem = 1.0
            # 2x2 Gauss points and weights
            gauss_pts = [(-1/np.sqrt(3), -1/np.sqrt(3)),
                         (1/np.sqrt(3), -1/np.sqrt(3)),
                         (1/np.sqrt(3), 1/np.sqrt(3)),
                         (-1/np.sqrt(3), 1/np.sqrt(3))]
            Nvals = [
                lambda xi, eta: np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]) * 0.25
                for _ in range(4)
            ]
            for (xi, eta) in gauss_pts:
                # Shape function derivatives w.r.t. natural coords
                dN_dxi = np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)]) * 0.25
                dN_deta = np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)]) * 0.25
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
                # Compute grad(h) at this Gauss point
                grad_h = gradN @ h_elem
                v_gp = -kr_elem * K @ grad_h  # Darcy velocity at Gauss point
                # Distribute/average to nodes (simple: add to all 4 nodes)
                for node in element_nodes:
                    velocity[node] += v_gp
                    count[node] += 1

    count[count == 0] = 1  # Avoid division by zero
    velocity /= count[:, None]
    return velocity

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
    for (xi, eta), w in zip(gauss_pts, weights):
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

def run_seepage_analysis(seep_data):
    """
    Standalone function to run seepage analysis.
    
    Args:
        seep_data: Dictionary containing all the seepage data 
    
    Returns:
        Dictionary containing solution results
    """
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
    print(f"Solving {flow_type.upper()} seepage problem...")
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

        head, A, q, total_flow = solve_unsaturated(
            nodes=nodes,
            elements=elements,
            bc_type=bc_type,
            bc_values=bc_values,
            kr0=kr0_per_element,
            h0=h0_per_element,
            k1_vals=k1,
            k2_vals=k2,
            angles=angle,
            max_iter=200,
            tol=1e-4,
            element_types=element_types
        )
        # Solve for potential function φ for flow lines
        dirichlet_phi_bcs = create_flow_potential_bc(nodes, elements, q, element_types=element_types)
        phi = solve_flow_function_unsaturated(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, dirichlet_phi_bcs, element_types)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        # Compute velocity, pass element-level kr0 and h0
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, kr0_per_element, h0_per_element, element_types)
    else:
        head, A, q, total_flow = solve_confined(nodes, elements, bc_type, bcs, k1, k2, angle, element_types)
        # Solve for potential function φ for flow lines
        dirichlet_phi_bcs = create_flow_potential_bc(nodes, elements, q, element_types=element_types)
        phi = solve_flow_function_confined(nodes, elements, k1, k2, angle, dirichlet_phi_bcs, element_types)
        print(f"phi min: {np.min(phi):.3f}, max: {np.max(phi):.3f}")
        # Compute velocity, don't pass kr0 and h0
        velocity = compute_velocity(nodes, elements, head, k1, k2, angle, element_types=element_types)

    gamma_w = unit_weight
    u = gamma_w * (head - nodes[:, 1])

    solution = {
        "head": head,
        "u": u,
        "velocity": velocity,
        "q": q,
        "phi": phi,
        "flowrate": total_flow
    }

    return solution

def export_solution_csv(filename, seep_data, solution):
    """Exports nodal results to a CSV file.
    
    Args:
        filename: Path to the output CSV file
        seep_data: Dictionary containing seepage data 
        solution: Dictionary containing solution results from run_analysis
    """
    import pandas as pd
    n_nodes = len(seep_data["nodes"])
    df = pd.DataFrame({
        "node_id": np.arange(1, n_nodes + 1),  # Generate 1-based node IDs for output
        "head": solution["head"],
        "u": solution["u"],
        "v_x": solution["velocity"][:, 0],
        "v_y": solution["velocity"][:, 1],
        "v_mag": np.linalg.norm(solution["velocity"], axis=1),
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
        seep_data: Dictionary containing seepage data 
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

