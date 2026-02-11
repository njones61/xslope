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

import warnings
from math import radians, degrees, sin, cos, tan, sqrt, atan2


import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import splu
from shapely.geometry import LineString, Point


def build_fem_data(slope_data, mesh=None):
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
            - unit_weight: float, unit weight of water
            - k_seismic: float, seismic coefficient (horizontal acceleration / gravity)
    """
    
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
    pp_options = [mat.get("pp_option", "none") for mat in materials]
    unique_pp_options = set([opt for opt in pp_options if opt != "none"])
    
    if len(unique_pp_options) > 1:
        raise ValueError(f"Mixed pore pressure options not allowed: {unique_pp_options}")
    
    pp_option = list(unique_pp_options)[0] if unique_pp_options else "none"
    
    for i, material in enumerate(materials):
        strength_option = material.get("strength_option", "mc")
        
        if strength_option == "mc":
            # Mohr-Coulomb: use c and phi directly
            c_by_mat[i] = material.get("c", 0.0)
            phi_by_mat[i] = material.get("phi", 0.0)
        elif strength_option == "cp":
            # c/p ratio: compute undrained strength based on depth
            cp_ratio = material.get("cp_ratio", 0.0)
            r_elev = material.get("r_elev", 0.0)
            
            # For c/p option, we need to assign strength per element based on element centroid
            # This will be handled when processing elements
            c_by_mat[i] = cp_ratio  # Store cp_ratio temporarily
            phi_by_mat[i] = 0.0     # Undrained analysis
        else:
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
    
    for elem_idx in range(n_elements):
        mat_id = element_materials[elem_idx] - 1  # Convert to 0-based
        material = materials[mat_id]
        strength_option = material.get("strength_option", "mc")
        
        if strength_option == "cp":
            cp_ratio = c_by_mat[mat_id]  # This is actually cp_ratio
            r_elev = material.get("r_elev", 0.0)
            
            # Compute element centroid
            elem_nodes = elements[elem_idx]
            elem_type = element_types[elem_idx]
            active_nodes = elem_nodes[:elem_type]  # Only use active nodes
            elem_coords = nodes[active_nodes]
            centroid_y = np.mean(elem_coords[:, 1])
            
            # Depth below reference elevation
            depth = max(0.0, r_elev - centroid_y)
            c_by_elem[elem_idx] = cp_ratio * depth
            phi_by_elem[elem_idx] = 0.0
        else:
            c_by_elem[elem_idx] = c_by_mat[mat_id]
            phi_by_elem[elem_idx] = phi_by_mat[mat_id]
    
    # Process pore pressures
    u = np.zeros(n_nodes)
    
    if pp_option == "piezo":
        # Find nodes and compute pore pressure from piezometric line
        # Assuming the piezometric line is stored in slope_data
        piezo_line_coords = None
        
        # Look for piezometric line in various possible locations
        if "piezo_line" in slope_data:
            piezo_line_coords = slope_data["piezo_line"]
        elif "profile_lines" in slope_data:
            # Check if one of the profile lines is designated as piezo
            for line in slope_data["profile_lines"]:
                if line.get('type') == 'piezo':
                    piezo_line_coords = line['coords']
                    break
        
        if piezo_line_coords:
            piezo_line = LineString(piezo_line_coords)
            gamma_water = slope_data.get("gamma_water", 9.81)
            
            for i, node in enumerate(nodes):
                node_point = Point(node)
                
                # Find closest point on piezometric line
                closest_point = piezo_line.interpolate(piezo_line.project(node_point))
                piezo_elevation = closest_point.y
                
                # Compute pore pressure (only positive values)
                if node[1] < piezo_elevation:
                    u[i] = gamma_water * (piezo_elevation - node[1])
                else:
                    u[i] = 0.0
    
    elif pp_option == "seep":
        # Use existing seep solution
        if "seepage_solution" in slope_data:
            seepage_solution = slope_data["seepage_solution"]
            if isinstance(seepage_solution, np.ndarray) and len(seepage_solution) == n_nodes:
                u = np.maximum(0.0, seepage_solution)  # Ensure non-negative
            else:
                print("Warning: Seepage solution dimensions don't match mesh nodes")
    
    # Process 1D reinforcement elements
    elements_1d = np.array([]).reshape(0, 3) if 'elements_1d' not in mesh else mesh['elements_1d']
    element_types_1d = np.array([]) if 'element_types_1d' not in mesh else mesh['element_types_1d'] 
    element_materials_1d = np.array([]) if 'element_materials_1d' not in mesh else mesh['element_materials_1d']
    
    n_1d_elements = len(elements_1d)
    
    t_allow_by_1d_elem = np.zeros(n_1d_elements)
    t_res_by_1d_elem = np.zeros(n_1d_elements)
    k_by_1d_elem = np.zeros(n_1d_elements)
    
    if n_1d_elements > 0 and "reinforcement_lines" in slope_data:
        reinforcement_lines = slope_data["reinforcement_lines"]
        
        for elem_idx in range(n_1d_elements):
            line_id = element_materials_1d[elem_idx] - 1  # Convert to 0-based
            
            if line_id < len(reinforcement_lines):
                line_data = reinforcement_lines[line_id]
                
                # Get element geometry
                elem_nodes = elements_1d[elem_idx]
                elem_type = element_types_1d[elem_idx]
                active_nodes = elem_nodes[:elem_type]
                elem_coords = nodes[active_nodes]
                
                # Compute element length and centroid
                if len(elem_coords) >= 2:
                    elem_length = np.linalg.norm(elem_coords[1] - elem_coords[0])
                    elem_centroid = np.mean(elem_coords, axis=0)
                    
                    # Compute distance from element centroid to line ends
                    x1, y1 = line_data.get("x1", 0), line_data.get("y1", 0)
                    x2, y2 = line_data.get("x2", 0), line_data.get("y2", 0)
                    
                    dist_to_left = np.linalg.norm(elem_centroid - [x1, y1])
                    dist_to_right = np.linalg.norm(elem_centroid - [x2, y2])
                    dist_to_nearest_end = min(dist_to_left, dist_to_right)
                    
                    # Get reinforcement properties
                    t_max = line_data.get("t_max", 0.0)
                    t_res = line_data.get("t_res", 0.0)
                    lp1 = line_data.get("lp1", 0.0)  # Pullout length left end
                    lp2 = line_data.get("lp2", 0.0)  # Pullout length right end
                    
                    # Use appropriate pullout length based on which end is closer
                    lp = lp1 if dist_to_left < dist_to_right else lp2
                    
                    # Compute allowable and residual tensile forces
                    if dist_to_nearest_end < lp:
                        # Within pullout zone - linear variation
                        t_allow_by_1d_elem[elem_idx] = t_max * (dist_to_nearest_end / lp)
                        t_res_by_1d_elem[elem_idx] = 0.0  # Sudden pullout failure
                    else:
                        # Beyond pullout zone - full capacity
                        t_allow_by_1d_elem[elem_idx] = t_max
                        t_res_by_1d_elem[elem_idx] = t_res
                    
                    # Compute axial stiffness
                    E = line_data.get("E", 2e11)  # Steel default
                    A = line_data.get("area", 1e-4)  # Default area
                    k_by_1d_elem[elem_idx] = E * A / elem_length
    
    # Set up boundary conditions
    
    # Step 1: Default to free (type 0)
    # Already initialized to zeros
    
    # Step 2: Fixed supports at bottom (type 1) - standard practice
    # Use global minimum y as bottom
    tolerance = 1e-6
    y_min = float(np.min(nodes[:, 1])) if len(nodes) > 0 else 0.0
    bottom_nodes = np.abs(nodes[:, 1] - y_min) < tolerance
    bc_type[bottom_nodes] = 1  # Fixed (u=0, v=0)
    
    # Step 3: X-roller supports at left and right sides (type 2) - standard practice
    # Use global min/max x to identify left/right boundaries
    if len(nodes) > 0:
        x_min = float(np.min(nodes[:, 0]))
        x_max = float(np.max(nodes[:, 0]))
        left_nodes = np.abs(nodes[:, 0] - x_min) < tolerance
        right_nodes = np.abs(nodes[:, 0] - x_max) < tolerance
        
        # Apply X-roller but preserve existing boundary conditions (fixed takes precedence at corners)
        left_not_fixed = left_nodes & (bc_type != 1)
        right_not_fixed = right_nodes & (bc_type != 1)
        
        bc_type[left_not_fixed] = 2   # X-roller (u=0, v=free)
        bc_type[right_not_fixed] = 2  # X-roller (u=0, v=free)
    
    # Step 4: Convert distributed loads to nodal forces (type 4)
    # Check for distributed loads (could be 'dloads', 'dloads2', or 'distributed_loads')
    distributed_loads = []
    if "dloads" in slope_data and slope_data["dloads"]:
        distributed_loads.extend(slope_data["dloads"])
    if "dloads2" in slope_data and slope_data["dloads2"]:
        distributed_loads.extend(slope_data["dloads2"])
    if "distributed_loads" in slope_data and slope_data["distributed_loads"]:
        distributed_loads.extend(slope_data["distributed_loads"])
    
    if distributed_loads:
        tolerance = 1e-1  # Tolerance for finding nodes on load lines (increased for better matching)
        
        for load_idx, load_line in enumerate(distributed_loads):
            # Handle different possible data structures
            if isinstance(load_line, dict) and "coords" in load_line:
                # Expected format: {"coords": [...], "loads": [...]}
                load_coords = load_line["coords"]
                load_values = load_line["loads"]
            elif isinstance(load_line, list):
                # Format from fileio: list of dicts with X, Y, Normal keys
                load_coords = [(pt["X"], pt["Y"]) for pt in load_line]
                load_values = [pt["Normal"] for pt in load_line]
            else:
                continue
            
            if len(load_coords) < 2 or len(load_values) < 2:
                continue
                
            load_linestring = LineString(load_coords)
            nodes_found = 0
            
            # Find nodes that lie on or near the load line
            for i, node in enumerate(nodes):
                node_point = Point(node)
                distance_to_line = load_linestring.distance(node_point)
                
                if distance_to_line <= tolerance:
                    # This node is on the load line
                    nodes_found += 1
                    # Find position along line and interpolate load
                    projected_distance = load_linestring.project(node_point)
                    
                    # Get segments and interpolate load value
                    segment_lengths = []
                    cumulative_length = 0
                    
                    for j in range(len(load_coords) - 1):
                        seg_length = np.linalg.norm(np.array(load_coords[j+1]) - np.array(load_coords[j]))
                        segment_lengths.append(seg_length)
                        cumulative_length += seg_length
                        
                        if projected_distance <= cumulative_length:
                            # Interpolate within this segment
                            local_distance = projected_distance - (cumulative_length - seg_length)
                            ratio = local_distance / seg_length if seg_length > 0 else 0
                            
                            load_at_node = load_values[j] * (1 - ratio) + load_values[j+1] * ratio
                            break
                    else:
                        # Use last load value if beyond end
                        load_at_node = load_values[-1]
                    
                    # Convert to nodal force using tributary length
                    # For simplicity, use average of adjacent segment lengths
                    tributary_length = np.mean(segment_lengths) if segment_lengths else 1.0
                    nodal_force_magnitude = load_at_node * tributary_length
                    
                    # Determine direction (perpendicular to ground surface)
                    # For now, assume vertical loading
                    bc_type[i] = 4  # Applied force
                    bc_values[i, 0] = 0.0  # No horizontal component
                    bc_values[i, 1] = -nodal_force_magnitude  # Downward
            
            pass
    
    # Print boundary condition summary
    bc_summary = np.bincount(bc_type, minlength=5)
    print(f"\nBoundary condition summary:")
    print(f"  Type 0 (free): {bc_summary[0]} nodes")
    print(f"  Type 1 (fixed): {bc_summary[1]} nodes") 
    print(f"  Type 2 (x-roller): {bc_summary[2]} nodes")
    print(f"  Type 3 (y-roller): {bc_summary[3]} nodes")
    print(f"  Type 4 (force): {bc_summary[4]} nodes")
    
    # Count non-zero forces
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        max_force = np.max(np.abs(bc_values[force_nodes]))
        print(f"  Maximum force magnitude: {max_force:.3f}")
    
    # Get other parameters
    unit_weight = slope_data.get("gamma_water", 9.81)
    k_seismic = slope_data.get("k_seismic", 0.0)
    
    # Construct fem_data dictionary
    fem_data = {
        "nodes": nodes,
        "elements": elements,
        "element_types": element_types,
        "element_materials": element_materials,
        "bc_type": bc_type,
        "bc_values": bc_values,
        "c_by_mat": c_by_mat,
        "phi_by_mat": phi_by_mat,
        "E_by_mat": E_by_mat,
        "nu_by_mat": nu_by_mat,
        "gamma_by_mat": gamma_by_mat,
        "material_names": material_names,
        "c_by_elem": c_by_elem,  # Element-wise cohesion (for c/p option)
        "phi_by_elem": phi_by_elem,  # Element-wise friction angle
        "u": u,
        "elements_1d": elements_1d,
        "element_types_1d": element_types_1d,
        "element_materials_1d": element_materials_1d,
        "t_allow_by_1d_elem": t_allow_by_1d_elem,
        "t_res_by_1d_elem": t_res_by_1d_elem,
        "k_by_1d_elem": k_by_1d_elem,
        "unit_weight": unit_weight,
        "k_seismic": k_seismic
    }
    
   
    return fem_data


def apply_boundary_conditions(K_global, F_global, bc_type, nodes):
    """
    Apply boundary conditions to global system using constraint elimination.
    
    This function applies boundary conditions by eliminating constrained degrees
    of freedom from the global stiffness matrix and load vector.
    
    Parameters:
        K_global: Global stiffness matrix (sparse or dense)
        F_global: Global load vector
        bc_type: Array of boundary condition types for each node:
                 0 = free (both u and v free)
                 1 = fixed (both u=0 and v=0) 
                 2 = x-roller (u=0, v free)
                 3 = y-roller (u free, v=0)
                 4 = force (both u and v free, external forces applied)
        nodes: Array of node coordinates (for reference)
    
    Returns:
        K_constrained: Constrained stiffness matrix (only free DOFs)
        F_constrained: Constrained load vector (only free DOFs)
        constraint_dofs: List of constrained DOF indices
    """
    
    n_nodes = len(nodes)
    n_dof = 2 * n_nodes
    
    # Identify constrained DOFs
    constraint_dofs = []
    
    for i in range(n_nodes):
        if bc_type[i] == 1:  # Fixed: both u and v constrained
            constraint_dofs.extend([2*i, 2*i+1])
        elif bc_type[i] == 2:  # X-roller: u constrained, v free
            constraint_dofs.append(2*i)
        elif bc_type[i] == 3:  # Y-roller: u free, v constrained
            constraint_dofs.append(2*i+1)
        # bc_type 0 and 4 are free DOFs - no constraints
    
    # Get free DOFs
    all_dofs = set(range(n_dof))
    constraint_dofs_set = set(constraint_dofs)
    free_dofs = sorted(all_dofs - constraint_dofs_set)
    
    # Extract free DOF submatrices
    if hasattr(K_global, 'toarray'):
        # Sparse matrix
        K_global_dense = K_global.toarray()
    else:
        K_global_dense = K_global
    
    # Extract submatrix for free DOFs only
    K_constrained = K_global_dense[np.ix_(free_dofs, free_dofs)]
    F_constrained = F_global[free_dofs]
    
    # Convert back to sparse if original was sparse and matrix is large
    if hasattr(K_global, 'toarray') and len(free_dofs) > 100:
        K_constrained = csr_matrix(K_constrained)
    
    return K_constrained, F_constrained, constraint_dofs


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

def solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=500, tolerance=1e-3):
    """
    Solve FEM using Griffiths & Lane (1999) viscoplastic algorithm.

    Implements the exact algorithm from the 1999 Geotechnique paper:
    - 8-node quadrilateral elements with reduced integration (4 Gauss points)
    - Viscoplastic stress redistribution with accumulated plastic strains
    - Non-convergence failure criterion
    - Pre-factored elastic stiffness matrix for efficiency
    - No damping (stability from dt parameter)
    - Direct solve each iteration (not residual-based)

    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        F (float): Shear strength reduction factor (c/F, tan(phi)/F)
        debug_level (int): 0=silent, 1=summary, 2=per-iteration
        max_iterations (int): Maximum viscoplastic iterations (default 500)
        tolerance (float): Convergence tolerance ||u_new - u_old|| / ||u_new|| (default 1e-3)

    Returns:
        dict: Solution dictionary with keys:
            - converged (bool): Whether iterations converged
            - iterations (int): Number of iterations used
            - displacements (ndarray): Nodal displacement vector
            - stresses (ndarray): Element stress array (n_elements, 4) [sig_x, sig_y, tau_xy, sig_vm] compression-positive
            - strains (ndarray): Element strain array (n_elements, 4) [eps_x, eps_y, gamma_xy, max_shear]
            - plastic_elements (ndarray): Boolean array of yielded elements
            - F (float): Applied strength reduction factor
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

    # Material properties
    c_by_elem = fem_data.get("c_by_elem", fem_data["c_by_mat"][element_materials - 1])
    phi_by_elem = fem_data.get("phi_by_elem", fem_data["phi_by_mat"][element_materials - 1])
    E_by_mat = fem_data["E_by_mat"]
    nu_by_mat = fem_data["nu_by_mat"]
    gamma_by_mat = fem_data["gamma_by_mat"]
    k_seismic = fem_data.get("k_seismic", 0.0)

    n_nodes = len(nodes)
    n_elements = len(elements)
    n_dof = 2 * n_nodes

    # Apply strength reduction (Griffiths & Lane 1999): c_r = c/F, phi_r = atan(tan(phi)/F)
    c_reduced = c_by_elem / F
    tan_phi_reduced = np.tan(np.radians(phi_by_elem)) / F
    phi_reduced = np.arctan(tan_phi_reduced)  # radians

    if debug_level >= 1:
        print(f"  c: {c_by_elem[0]:.1f} -> {c_reduced[0]:.1f}")
        print(f"  phi: {phi_by_elem[0]:.1f} -> {np.degrees(phi_reduced[0]):.1f}")

    # ---- Step 1: Build K_global (elastic, constant) ----
    K_global = build_global_stiffness(nodes, elements, element_types,
                                      element_materials, E_by_mat, nu_by_mat)

    # ---- Step 2: Build gravity load vector ----
    F_gravity = build_gravity_loads(nodes, elements, element_types,
                                    element_materials, gamma_by_mat, k_seismic)

    # Add boundary condition forces
    for i in range(n_nodes):
        if bc_type[i] == 4:  # Force boundary condition
            F_gravity[2*i] += bc_values[i, 0]
            F_gravity[2*i+1] += bc_values[i, 1]

    # ---- Step 3: Identify free/constrained DOFs ONCE ----
    constraint_dofs = []
    for i in range(n_nodes):
        if bc_type[i] == 1:  # Fixed
            constraint_dofs.extend([2*i, 2*i+1])
        elif bc_type[i] == 2:  # X-roller
            constraint_dofs.append(2*i)
        elif bc_type[i] == 3:  # Y-roller
            constraint_dofs.append(2*i+1)

    constraint_set = set(constraint_dofs)
    free_dofs = np.array(sorted(set(range(n_dof)) - constraint_set))
    n_free = len(free_dofs)

    # ---- Step 4: Extract K_free and PRE-FACTORIZE ----
    if hasattr(K_global, 'toarray'):
        K_dense = K_global.toarray()
    else:
        K_dense = K_global

    K_free = csr_matrix(K_dense[np.ix_(free_dofs, free_dofs)])
    K_factor = splu(K_free.tocsc())

    F_grav_free = F_gravity[free_dofs]

    if debug_level >= 1:
        print(f"  DOFs: {n_dof} total, {n_free} free, {len(constraint_dofs)} constrained")
        print(f"  K factorized (reused for all iterations)")

    # ---- Step 5: Compute dt from material properties (Griffiths p62.f90) ----
    dt = 1.0e15
    for mat_idx in range(len(E_by_mat)):
        E = E_by_mat[mat_idx]
        nu = nu_by_mat[mat_idx]
        ddt = 4.0 * (1.0 + nu) / (3.0 * E)
        if ddt < dt:
            dt = ddt

    if debug_level >= 2:
        print(f"  dt = {dt:.3e}")

    # ---- Step 6: Pre-compute element B matrices and D matrices at Gauss points ----
    gauss_points_2x2, gauss_weights_2x2 = get_gauss_points_2x2()

    # Store per-element, per-GP: B matrix, weight (det_J * gauss_weight), D matrix
    elem_gp_data = []
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        D = build_constitutive_matrix(E, nu)

        elem_nodes_idx = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes_idx]

        gp_list = []
        if elem_type == 8:
            for gp_idx in range(4):
                xi, eta = gauss_points_2x2[gp_idx]
                B, det_J = _compute_B_and_detJ_quad8(elem_coords, xi, eta)
                weight = gauss_weights_2x2[gp_idx] * abs(det_J)
                gp_list.append({'B': B, 'weight': weight, 'D': D, 'dof_indices': _elem_dof_indices(elem_nodes_idx)})
        elif elem_type == 3:
            B, area = compute_B_matrix_triangle(elem_coords)
            gp_list.append({'B': B, 'weight': area, 'D': D, 'dof_indices': _elem_dof_indices(elem_nodes_idx)})
        else:
            B, area = compute_B_matrix_triangle(elem_coords)
            gp_list.append({'B': B, 'weight': area, 'D': D, 'dof_indices': _elem_dof_indices(elem_nodes_idx)})

        elem_gp_data.append(gp_list)

    # ---- Step 7: Initial elastic solution ----
    u_free = K_factor.solve(F_grav_free)

    # Full displacement vector
    u = np.zeros(n_dof)
    u[free_dofs] = u_free

    if debug_level >= 1:
        print(f"  Initial elastic: max|u| = {np.max(np.abs(u)):.6f}")

    # ---- Step 8: Initialize viscoplastic strains (zero) ----
    # evp[elem_idx][gp_idx] = array of shape (3,)
    evp = []
    for elem_idx in range(n_elements):
        n_gp = len(elem_gp_data[elem_idx])
        evp.append([np.zeros(3) for _ in range(n_gp)])

    # ---- Step 9: Viscoplastic iteration loop ----
    converged = False
    iteration = 0

    for iteration in range(max_iterations):
        # Build body load correction from accumulated viscoplastic strains
        loads = F_gravity.copy()

        n_yielding = 0
        n_total_gp = 0

        for elem_idx in range(n_elements):
            gp_data_list = elem_gp_data[elem_idx]
            elem_type = element_types[elem_idx]
            elem_nodes_idx = elements[elem_idx][:elem_type]

            for gp_idx, gp_data in enumerate(gp_data_list):
                B = gp_data['B']
                D = gp_data['D']
                weight = gp_data['weight']
                dof_idx = gp_data['dof_indices']
                n_total_gp += 1

                # a. Total strains from current displacements: eps = B * u_elem
                u_elem = u[dof_idx]
                eps_total = B @ u_elem

                # b. Elastic strains = total - viscoplastic (Smith & Griffiths p62.f90)
                eps_elastic = eps_total - evp[elem_idx][gp_idx]

                # c. Actual stress from elastic strains (tension-positive)
                sigma_tp = D @ eps_elastic

                # d. Convert to compression-positive for yield check
                sigma_cp = np.array([-sigma_tp[0], -sigma_tp[1], sigma_tp[2]])

                # e. Check Mohr-Coulomb yield
                c_r = c_reduced[elem_idx]
                phi_r = phi_reduced[elem_idx]
                f_yield = check_mohr_coulomb_cp(sigma_cp, c_r, phi_r)

                # f. If yielding, compute viscoplastic strain increment
                if f_yield > 0:
                    n_yielding += 1
                    # Flow direction from actual stress (tension-positive)
                    flow_tp = compute_flow_vector_tp(sigma_tp, psi=0.0)

                    # Viscoplastic strain increment: delta_evp = f_yield * flow * dt
                    delta_evp = f_yield * flow_tp * dt
                    evp[elem_idx][gp_idx] += delta_evp

                # f. Body load correction: loads += B^T * D * evp * weight
                evp_gp = evp[elem_idx][gp_idx]
                if np.any(np.abs(evp_gp) > 1e-30):
                    correction = B.T @ (D @ evp_gp) * weight
                    loads[dof_idx] += correction

        # Solve K * u_new = loads
        loads_free = loads[free_dofs]
        u_free_new = K_factor.solve(loads_free)

        u_new = np.zeros(n_dof)
        u_new[free_dofs] = u_free_new

        # Convergence check: ||u_new - u|| / ||u_new|| < tolerance
        norm_u_new = np.linalg.norm(u_new)
        norm_diff = np.linalg.norm(u_new - u)

        if norm_u_new > 1e-30:
            relative_change = norm_diff / norm_u_new
        else:
            relative_change = norm_diff

        if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
            print(f"  Iter {iteration+1:4d}: ||du||/||u|| = {relative_change:.3e}, "
                  f"yielding = {n_yielding}/{n_total_gp}, max|u| = {np.max(np.abs(u_new)):.6f}")

        if relative_change < tolerance:
            converged = True
            u = u_new
            if debug_level >= 1:
                print(f"  Converged after {iteration+1} iterations (||du||/||u|| = {relative_change:.3e})")
            break

        u = u_new

    if not converged and debug_level >= 1:
        print(f"  Did NOT converge after {max_iterations} iterations (||du||/||u|| = {relative_change:.3e})")

    # ---- Step 10: Compute final stresses, strains, plastic elements ----
    final_stresses = np.zeros((n_elements, 4))  # [sig_x, sig_y, tau_xy, sig_vm] compression-positive
    plastic_elements = np.zeros(n_elements, dtype=bool)

    for elem_idx in range(n_elements):
        gp_data_list = elem_gp_data[elem_idx]
        n_gp = len(gp_data_list)
        stress_avg_cp = np.zeros(3)

        for gp_idx, gp_data in enumerate(gp_data_list):
            B = gp_data['B']
            D = gp_data['D']
            dof_idx = gp_data['dof_indices']

            u_elem = u[dof_idx]
            eps_total = B @ u_elem
            eps_elastic = eps_total - evp[elem_idx][gp_idx]
            sigma_tp = D @ eps_elastic
            sigma_cp = np.array([-sigma_tp[0], -sigma_tp[1], sigma_tp[2]])
            stress_avg_cp += sigma_cp

        stress_avg_cp /= n_gp
        sig_x, sig_y, tau_xy = stress_avg_cp
        sig_vm = sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
        final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]

        f_yield = check_mohr_coulomb_cp(stress_avg_cp, c_reduced[elem_idx], phi_reduced[elem_idx])
        plastic_elements[elem_idx] = f_yield > 1e-8

    strains = compute_strains(nodes, elements, element_types, u)

    n_plastic = np.sum(plastic_elements)
    if debug_level >= 1:
        print(f"  Plastic elements: {n_plastic}/{n_elements}")
        print(f"  Max displacement: {np.max(np.abs(u)):.6f}")

    return {
        "converged": converged,
        "iterations": iteration + 1,
        "displacements": u,
        "stresses": final_stresses,
        "strains": strains,
        "plastic_elements": plastic_elements,
        "yield_function": np.array([check_mohr_coulomb_cp(final_stresses[i, :3], c_reduced[i], phi_reduced[i]) for i in range(n_elements)]),
        "max_displacement": np.max(np.abs(u)),
        "plastic_strains": {i: np.array(evp[i]) for i in range(n_elements)},
        "algorithm": "Griffiths & Lane (1999) Viscoplastic",
        "F": F,
        "residual": relative_change if 'relative_change' in locals() else 0.0,
        "plastic_fraction": n_plastic / n_elements if n_elements > 0 else 0.0,
    }


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


def _elem_dof_indices(elem_nodes):
    """Get global DOF indices for element nodes."""
    dof_idx = np.zeros(2 * len(elem_nodes), dtype=int)
    for i, node in enumerate(elem_nodes):
        dof_idx[2*i] = 2 * node
        dof_idx[2*i+1] = 2 * node + 1
    return dof_idx


def compute_flow_vector_tp(stress_tp, psi=0.0):
    """
    Compute non-associated flow direction in tension-positive convention.

    For psi=0 (non-associated, no dilation): purely deviatoric flow.

    The plastic potential is g = tau_max - sigma_mean * sin(psi) (compression-positive).
    We compute dg/d(sigma) in compression-positive, then convert to tension-positive.

    Args:
        stress_tp: [sig_x, sig_y, tau_xy] in tension-positive convention
        psi: Dilation angle in radians (0 for non-associated)

    Returns:
        flow_tp: [dg/dsig_x, dg/dsig_y, dg/dtau_xy] in tension-positive convention
    """
    # Convert to compression-positive
    sig_x_cp = -stress_tp[0]
    sig_y_cp = -stress_tp[1]
    tau_xy = stress_tp[2]

    tau_max = sqrt(((sig_x_cp - sig_y_cp) / 2.0)**2 + tau_xy**2)

    if tau_max < 1e-20:
        return np.zeros(3)

    sin_psi = sin(psi)

    # Derivatives of g w.r.t. compression-positive stresses:
    # dg/dsig_x_cp = (sig_x_cp - sig_y_cp) / (4*tau_max) - sin(psi)/2
    # dg/dsig_y_cp = -(sig_x_cp - sig_y_cp) / (4*tau_max) - sin(psi)/2
    # dg/dtau_xy   = tau_xy / tau_max
    flow_x_cp = (sig_x_cp - sig_y_cp) / (4.0 * tau_max) - sin_psi * 0.5
    flow_y_cp = -(sig_x_cp - sig_y_cp) / (4.0 * tau_max) - sin_psi * 0.5
    flow_xy_cp = tau_xy / tau_max

    # Convert to tension-positive: d/dsig_tp = -d/dsig_cp for normals; shear unchanged
    flow_x_tp = -flow_x_cp
    flow_y_tp = -flow_y_cp
    flow_xy_tp = flow_xy_cp

    return np.array([flow_x_tp, flow_y_tp, flow_xy_tp])


def solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05, debug_level=0,
               max_iterations=500, convergence_tol=1e-3):
    """
    Shear Strength Reduction Method using bisection on solve_fem convergence.

    Finds the critical strength reduction factor F where the viscoplastic
    algorithm transitions from converging (stable) to non-converging (failure).

    Parameters:
        fem_data (dict): FEM data from build_fem_data
        F_min (float): Lower bound for F (must converge). Default 1.0.
        F_max (float): Upper bound for F (should not converge). Default 2.0.
        tolerance (float): Bisection stops when F_right - F_left < tolerance. Default 0.05.
        debug_level (int): Verbosity (0=silent, 1=summary, 2=detailed)
        max_iterations (int): Max viscoplastic iterations passed to solve_fem
        convergence_tol (float): Convergence tolerance passed to solve_fem

    Returns:
        dict: Result with keys FS, converged, last_solution, final_interval, etc.
    """

    if debug_level >= 1:
        print("=== SSRM Analysis (Griffiths & Lane 1999) ===")
        print(f"  Bisection range: [{F_min:.2f}, {F_max:.2f}], tolerance: {tolerance}")

    F_left = F_min
    F_right = F_max

    # Verify lower bound converges
    solution_min = solve_fem(fem_data, F=F_min, debug_level=max(0, debug_level-1),
                             max_iterations=max_iterations, tolerance=convergence_tol)
    if not solution_min["converged"]:
        return {
            "converged": False,
            "error": f"F_min = {F_min} does not converge - slope unstable at F=1",
            "FS": None
        }

    # Verify upper bound does not converge
    solution_max = solve_fem(fem_data, F=F_max, debug_level=max(0, debug_level-1),
                             max_iterations=max_iterations, tolerance=convergence_tol)
    if solution_max["converged"]:
        if debug_level >= 1:
            print(f"  Warning: F_max = {F_max} still converges - increase F_max")
        return {
            "converged": True,
            "FS": F_max,
            "last_solution": solution_max,
            "note": f"Slope stable up to F = {F_max}"
        }

    last_converged_solution = solution_min
    iteration = 0

    while (F_right - F_left) > tolerance and iteration < 50:
        F_mid = (F_left + F_right) / 2.0

        if debug_level >= 1:
            print(f"\n  SSRM step {iteration+1}: F = {F_mid:.4f}  [{F_left:.4f}, {F_right:.4f}]")

        solution = solve_fem(fem_data, F=F_mid, debug_level=max(0, debug_level-1),
                             max_iterations=max_iterations, tolerance=convergence_tol)

        if solution["converged"]:
            F_left = F_mid
            last_converged_solution = solution
            if debug_level >= 1:
                print(f"    -> Converged in {solution['iterations']} iters")
        else:
            F_right = F_mid
            if debug_level >= 1:
                print(f"    -> Did NOT converge ({solution['iterations']} iters)")

        iteration += 1

    critical_FS = F_left

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
        "method": "Griffiths & Lane (1999) Viscoplastic SSRM"
    }


def build_global_stiffness(nodes, elements, element_types, element_materials, E_by_mat, nu_by_mat):
    """
    Build global stiffness matrix using existing FE implementation for proper 8-node quad support.
    """
    # Use existing stiffness functions (now they are in this same file after consolidation)
    
    n_nodes = len(nodes)
    n_dof = 2 * n_nodes
    
    K_global = lil_matrix((n_dof, n_dof))
    
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
            if elem_type == 3:  # Triangular elements
                K_elem = build_triangle_stiffness_corrected(elem_coords, E, nu)
            elif elem_type == 8:  # 8-node quadrilateral elements - use corrected Griffiths version
                K_elem = build_quad8_stiffness_reduced_integration_corrected(elem_coords, E, nu)
            elif elem_type in [4, 6, 9]:  # Other elements - use simple triangle implementation
                K_elem = build_triangle_stiffness_corrected(elem_coords, E, nu)
            else:
                print(f"Warning: Element type {elem_type} not supported")
                continue
        except Exception as e:
            print(f"Error building stiffness for element {elem_idx}, type {elem_type}: {e}")
            continue
        
        # Assemble into global matrix
        for i in range(elem_type):
            for j in range(elem_type):
                node_i = elem_nodes[i]
                node_j = elem_nodes[j]
                
                for di in range(2):
                    for dj in range(2):
                        global_i = 2 * node_i + di
                        global_j = 2 * node_j + dj
                        local_i = 2 * i + di
                        local_j = 2 * j + dj
                        
                        if local_i < K_elem.shape[0] and local_j < K_elem.shape[1]:
                            K_global[global_i, global_j] += K_elem[local_i, local_j]
    
    return K_global.tocsr()



def build_triangle_stiffness(coords, E, nu):
    """
    Build stiffness matrix for triangular element (plane strain).
    """
    x1, y1 = coords[0]
    x2, y2 = coords[1] 
    x3, y3 = coords[2]
    
    # Area
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    if area < 1e-12:
        print(f"Warning: Very small element area: {area}")
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


def build_gravity_loads(nodes, elements, element_types, element_materials, gamma_by_mat, k_seismic):
    """
    Build gravity load vector using Griffiths & Lane (1999) approach.
    
    Uses equation 3 from the paper: p(e) = γ ∫[Ve] N^T d(vol)
    This integrates shape functions over each element to properly distribute gravity loads.
    """
    n_nodes = len(nodes)
    F_gravity = np.zeros(2 * n_nodes)
    
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
                F_gravity[2*node + 1] -= load  # Vertical (negative = downward)
                F_gravity[2*node] += k_seismic * load  # Horizontal seismic
                
        elif elem_type == 8:  # 8-node quad
            # For 8-node quads, use 2x2 Gauss integration as in Griffiths
            # This properly weights corner vs midside nodes
            
            # Gauss points for 2x2 integration
            gauss_coord = 1.0 / np.sqrt(3.0)
            xi_points = np.array([-gauss_coord, gauss_coord])
            eta_points = np.array([-gauss_coord, gauss_coord])
            weights = np.array([1.0, 1.0])
            
            # Initialize element load vector
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
                    
                    # Accumulate load contribution: w * det(J) * γ * N
                    for k in range(8):
                        elem_loads[2*k + 1] -= w * det_J * gamma * N[k]  # Vertical
                        elem_loads[2*k] += w * det_J * gamma * k_seismic * N[k]  # Horizontal
            
            # Add element loads to global vector
            for i, node in enumerate(elem_nodes):
                F_gravity[2*node] += elem_loads[2*i]
                F_gravity[2*node + 1] += elem_loads[2*i + 1]
                
        elif elem_type == 4:  # 4-node quad (if used)
            # For 4-node quads, use 2x2 Gauss integration
            area = compute_quad_area(elem_coords)
            # Simple equal distribution for now (can be refined)
            load_per_node = gamma * area / 4.0
            
            for i, node in enumerate(elem_nodes):
                F_gravity[2*node + 1] -= load_per_node
                F_gravity[2*node] += k_seismic * load_per_node
        else:
            # Fallback for other element types
            if elem_type >= 3:
                # Triangle area calculation
                x1, y1 = elem_coords[0]
                x2, y2 = elem_coords[1]
                x3, y3 = elem_coords[2]
                area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
                
                load_per_node = gamma * area / elem_type
                
                for i, node in enumerate(elem_nodes):
                    F_gravity[2*node + 1] -= load_per_node
                    F_gravity[2*node] += k_seismic * load_per_node
    
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


def compute_quad_area(coords):
    """
    Compute area of quadrilateral (approximate).
    """
    if len(coords) >= 4:
        # Use shoelace formula for polygon area
        x = coords[:4, 0]
        y = coords[:4, 1]
        return 0.5 * abs(sum(x[i]*y[i+1] - x[i+1]*y[i] for i in range(-1, 3)))
    else:
        return 0.0



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


def compute_gauss_point_coordinates_quad8(elem_coords, xi, eta):
    """
    Compute physical coordinates of a Gauss point in an 8-node quadrilateral.
    
    Args:
        elem_coords: Array of element node coordinates (8x2)
        xi, eta: Natural coordinates of Gauss point
        
    Returns:
        Physical coordinates [x, y] of the Gauss point
    """
    # 8-node quadrilateral shape functions
    N = np.zeros(8)
    N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1)
    N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1)
    N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1)
    N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1)
    N[4] = 0.5 * (1 - xi*xi) * (1 - eta)
    N[5] = 0.5 * (1 + xi) * (1 - eta*eta)
    N[6] = 0.5 * (1 - xi*xi) * (1 + eta)
    N[7] = 0.5 * (1 - xi) * (1 - eta*eta)
    
    # Compute physical coordinates
    gauss_coords = np.zeros(2)
    for i in range(8):
        gauss_coords += N[i] * elem_coords[i]
    
    return gauss_coords



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
    # Add numerical stability check for near-incompressible materials
    if nu >= 0.45:
        print(f"Warning: Poisson's ratio {nu:.3f} is close to incompressible limit (0.5)")
        print("Consider using nu <= 0.4 for better numerical stability")
    
    # Optional: Add small regularization to prevent singularity
    # nu_effective = min(nu, 0.495)  # Cap at safe value if needed
    
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    # Standard tension-positive convention (σ > 0 in tension, σ < 0 in compression)
    return D



def check_mohr_coulomb_cp(stress_cp, c, phi):
    """Mohr-Coulomb yield function for compression-positive stresses.
    
    For compression-positive convention (compression > 0, tension < 0):
    F = tau_max - sigma_mean * sin(phi) - c * cos(phi)
    
    Where:
    - tau_max = maximum shear stress = sqrt((sig_x - sig_y)^2/4 + tau_xy^2)
    - sigma_mean = mean normal stress = (sig_x + sig_y)/2
    - Positive F indicates yielding
    
    Args:
        stress_cp: Array [sig_x, sig_y, tau_xy] in compression-positive convention
        c: Cohesion
        phi: Friction angle in radians
        
    Returns:
        F: Yield function value (F > 0 means yielding)
    """
    sig_x, sig_y, tau_xy = stress_cp
    sig_mean = (sig_x + sig_y) / 2.0
    tau_max = sqrt(((sig_x - sig_y) / 2.0)**2 + tau_xy**2)
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    F = tau_max - sig_mean * sin_phi - c * cos_phi
    return F



def compute_strains(nodes, elements, element_types, displacements):
    """
    Compute element strains for visualization.
    """
    n_elements = len(elements)
    strains = np.zeros((n_elements, 4))  # [eps_x, eps_y, gamma_xy, max_shear_strain]
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Get element displacements
        elem_disp = np.zeros(2 * elem_type)
        for i, node in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node]
            elem_disp[2*i+1] = displacements[2*node+1]
        
        # Compute strains
        if elem_type == 3:
            element_strains = compute_triangle_strains_manual(elem_coords, elem_disp)
        elif elem_type == 8:
            # For 8-node quad, compute strain at centroid
            xi, eta = 0.0, 0.0  # Centroid
            element_strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta)
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