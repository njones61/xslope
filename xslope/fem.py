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

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
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
                if hasattr(line, 'type') and line.type == 'piezo':
                    piezo_line_coords = line
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
        # Use existing seepage solution
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

def solve_fem(fem_data, F=1.0, debug_level=0, abort_after=-1, iteration_print_frequency=5, dt_max=1e-5, max_iterations=60, tolerance=5e-4, damping_factor=0.9, plastic_strain_cap=0.005):
    """
    Solve FEM using Perzyna visco-plastic algorithm exactly as in Griffiths & Lane (1999).
    
    This implements the exact algorithm from the 1999 Geotechnique paper:
    - 8-node quadrilateral elements with reduced integration (4 Gauss points)
    - Perzyna visco-plastic stress redistribution
    - Non-convergence failure criterion (1000 iteration limit)
    - No plastic stiffness reduction
    
    Parameters:
        fem_data (dict): FEM data dictionary
        F (float): Shear strength reduction factor
        debug_level (int): Verbosity level
        abort_after (int): Abort after this many iterations. -1 = no abort (default)
                          0 = abort after gravity loading (before plasticity check)
                          1 = abort after first plasticity iteration
                          etc.
        iteration_print_frequency (int): Print iteration info every N iterations (default=1)
        dt_max (float): Maximum pseudo-time step for Perzyna updates; dt = min(dt_base, dt_max).
                        Defaults to 1e-6.
        max_iterations (int): Maximum Perzyna iterations. Defaults to 60.
        tolerance (float): Convergence tolerance on relative displacement change. Defaults to 5e-4.
        damping_factor (float): Under-relaxation factor for displacement update (0<d<=1).
                               Lower for more damping (e.g., 0.8–0.9). Defaults to 0.95.
        plastic_strain_cap (float|None): Max per-Gauss plastic strain increment magnitude per
                               iteration (cap on |erate*dt|). None disables capping.
        
    Returns:
        dict: Solution dictionary with convergence status
    """
    
    if debug_level >= 1:
        print(f"=== Perzyna Visco-Plastic FEM Analysis (F={F:.3f}) ===")
    
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
    u_nodal = fem_data["u"]
    k_seismic = fem_data.get("k_seismic", 0.0)
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    n_dof = 2 * n_nodes
    
    # Apply strength reduction (Griffiths & Lane 1999 approach)
    c_reduced = c_by_elem / F
    tan_phi_reduced = np.tan(np.radians(phi_by_elem)) / F
    phi_reduced = np.arctan(tan_phi_reduced)  # Keep in radians for yield functions
    
    if debug_level >= 2:
        print(f"Original c range: [{np.min(c_by_elem):.1f}, {np.max(c_by_elem):.1f}]")
        print(f"Reduced c range: [{np.min(c_reduced):.1f}, {np.max(c_reduced):.1f}]")
        print(f"Original phi: {phi_by_elem[0]:.1f}°")
        print(f"Reduced phi: {np.degrees(phi_reduced[0]):.1f}°")
        print(f"Original φ range: [{np.min(phi_by_elem):.1f}°, {np.max(phi_by_elem):.1f}°]")
        print(f"Reduced φ range: [{np.min(np.degrees(phi_reduced)):.1f}°, {np.max(np.degrees(phi_reduced)):.1f}°]")
    
    # Build global stiffness matrix (elastic, constant throughout)
    K_global = build_global_stiffness(nodes, elements, element_types, 
                                             element_materials, E_by_mat, nu_by_mat)
    
    # Build gravity load vector
    F_gravity = build_gravity_loads(nodes, elements, element_types, 
                                           element_materials, gamma_by_mat, k_seismic)
    
    # Boundary conditions will be applied in each iteration using apply_boundary_conditions
    

    
    # Pseudo-time step for numerical stability (not real time - this is steady-state)
    # Griffiths & Lane approach: start with large value, then calculate based on material properties
    dt_base = 1.0e15  # Large initial value as in p62.f90 line 19
    
    # Calculate time step based on material properties (p62.f90 lines 72-73)
    # Use the first material's properties for time step calculation
    E = E_by_mat[0]  # Young's modulus of first material
    nu = nu_by_mat[0]  # Poisson's ratio of first material
    ddt = 4.0 * (1.0 + nu) / (3.0 * E)  # d4*(one+prop(3))/(d3*prop(2))
    if ddt < dt_base:
        dt_base = ddt
    
    # Final time step control: dt is min of material-based dt_base and user dt_max
    dt = min(dt_base, dt_max)
    
    # Debug prints for dt calculation
    if debug_level >= 2:
        print(f"Time step estimate: E={E:.3g}, nu={nu:.3g}, ddt={ddt:.3e}")
        print(f"dt_base={dt_base:.3e}, dt_max={dt_max:.3e}, dt={dt:.3e}")
    
    # Phase 1: Establish K0 (gravity) stress state from elastic solution
    if debug_level >= 1:
        print("Phase 1: Establishing K0 stress state by elastic gravity loading...")
    initial_displacements, stress_state = establish_k0_stress_state(
        K_global, F_gravity, bc_type, nodes, elements, element_types,
        element_materials, E_by_mat, nu_by_mat, gamma_by_mat, u_nodal, debug_level=max(0, debug_level-1)
    )
    
    # Debug gravity loads and element areas
    if debug_level >= 1:
        print(f"  Debug: Checking gravity load calculation")
        sample_elem = 129  # should be close to shear zone for tri3 mesh
        if sample_elem < len(elements):
            elem_type = element_types[sample_elem]
            elem_nodes = elements[sample_elem][:elem_type]
            elem_coords = nodes[elem_nodes]
            if elem_type == 8:
                area = compute_quad_area(elem_coords)
            else:
                # Triangle area
                x1, y1 = elem_coords[0]
                x2, y2 = elem_coords[1] 
                x3, y3 = elem_coords[2]
                area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            gamma = gamma_by_mat[element_materials[sample_elem] - 1]
            load_per_node = gamma * area / elem_type
            print(f"    Element {sample_elem}: area={area:.2f}, gamma={gamma}, load_per_node={load_per_node:.1f}")
    
    # Calculate yield function values for all elements after gravity loading
    yield_function_values = np.zeros(n_elements)
    # Collect all σyy values for diagnostics
    all_sigma_yy = []
    total_gauss_points = 0
    yielded_gauss_points = 0
    
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        # Use first Gauss point stress for yield function (or average for quads)
        if elem_type == 8:  # 8-node quad - average over Gauss points
            elem_stress_avg = np.mean(stress_state['element_stresses'][elem_idx, :4, :], axis=0)
            # Collect all σyy from Gauss points
            for gp in range(4):
                all_sigma_yy.append(stress_state['element_stresses'][elem_idx, gp, 1])  # σyy is index 1
                total_gauss_points += 1
                # Check if this Gauss point yields
                gp_yield = check_mohr_coulomb_cp(
                    stress_state['element_stresses'][elem_idx, gp, :], 
                    c_reduced[elem_idx], phi_reduced[elem_idx])
                if gp_yield > 0:
                    yielded_gauss_points += 1
        else:  # Triangle or other - use first Gauss point
            elem_stress_avg = stress_state['element_stresses'][elem_idx, 0, :]
            all_sigma_yy.append(stress_state['element_stresses'][elem_idx, 0, 1])  # σyy
            total_gauss_points += 1
            # Check if yields
            gp_yield = check_mohr_coulomb_cp(
                stress_state['element_stresses'][elem_idx, 0, :],
                c_reduced[elem_idx], phi_reduced[elem_idx])
            if gp_yield > 0:
                yielded_gauss_points += 1
        
        # Calculate yield function with reduced strength parameters
        yield_function_values[elem_idx] = check_mohr_coulomb_cp(
            elem_stress_avg, c_reduced[elem_idx], phi_reduced[elem_idx])
    
    # Diagnostic 1: Min/max σyy after gravity
    all_sigma_yy = np.array(all_sigma_yy)
    min_sigma_yy = np.min(all_sigma_yy)
    max_sigma_yy = np.max(all_sigma_yy)
    
    if debug_level >= 1:
        print(f"\n=== After Gravity Loading ====")
        print(f"  Min σyy: {min_sigma_yy:.3f} kPa")
        print(f"  Max σyy: {max_sigma_yy:.3f} kPa")
    
    # Diagnostic 2: Fraction of yielded Gauss points
    yielding_fraction_initial = yielded_gauss_points / total_gauss_points if total_gauss_points > 0 else 0
    if debug_level >= 1:
        print(f"\n=== After Yield Check (F > 0) ===")
        print(f"  Gauss points meeting yield criterion (F>0): {yielded_gauss_points}/{total_gauss_points} ({yielding_fraction_initial*100:.1f}%)")
        print(f"  Note: These points satisfy F>0 but haven't developed plastic strains yet")
        
    
    # Check for early abort after gravity loading
    if abort_after == 0:
        if debug_level >= 1:
            print("Aborting after gravity loading (abort_after=0)")
        
        # Compute stresses and strains for output
        final_stresses, plastic_elements = compute_final_state_perzyna(
            nodes, elements, element_types, element_materials,
            initial_displacements, {}, c_reduced, phi_reduced,
            E_by_mat, nu_by_mat, u_nodal, stress_state)
        
        strains = compute_strains(nodes, elements, element_types, initial_displacements)
        
        return {
            "converged": True,
            "iterations": 0,
            "displacements": initial_displacements,
            "stresses": final_stresses,
            "strains": strains,
            "plastic_elements": plastic_elements,
            "yield_function": yield_function_values,
            "max_displacement": np.max(np.abs(initial_displacements)),
            "plastic_strains": {},
            "algorithm": "Perzyna Visco-Plastic (aborted after gravity)",
            "aborted": True,
            "abort_after": abort_after
        }
    
    if debug_level >= 1:
        print("Phase 2: Starting Perzyna strength reduction analysis...")
        
    displacements = initial_displacements.copy()
    displacements_prev = initial_displacements.copy()  # Track previous displacements
    plastic_strains = {}  # Store plastic strains at each Gauss point
    
    # Initialize total stress state from K₀ (this will be updated incrementally)
    current_stress_state = {
        'element_stresses': stress_state['element_stresses'].copy(),
        'plastic_state': np.zeros((n_elements, 4), dtype=bool)
    }
    
    # Initialize plastic strain storage
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        if elem_type == 8:  # 8-node quad
            n_gauss = 4  # Reduced integration
        elif elem_type == 3:  # 3-node triangle
            n_gauss = 1
        else:
            n_gauss = 1
        
        plastic_strains[elem_idx] = np.zeros((n_gauss, 3))  # [eps_x, eps_y, gamma_xy] plastic
    
    converged = False
    
    # Track data for CSV output
    csv_data = []
    
    for iteration in range(max_iterations):
        if debug_level >= 3:
            print(f"\n--- Iteration {iteration + 1} ---")
        
        # Build load vector: maintain gravity (constant external load) + plastic corrections
        F_total = F_gravity.copy()
        
        # Add plastic strain corrections (Perzyna stress redistribution)
        F_plastic_correction = compute_plastic_load_correction_perzyna(
            nodes, elements, element_types, element_materials, 
            plastic_strains, E_by_mat, nu_by_mat, dt)
        
        if debug_level >= 2:
            print(f"Plastic correction norm: {np.linalg.norm(F_plastic_correction):.2e}")
        
        F_total += F_plastic_correction
        
        # DEBUG: Detailed analysis for first iteration only
        if iteration == 0 and debug_level >= 1:
            print(f"\n=== DEBUGGING FIRST ITERATION ===")
            print(f"Gravity load norm: {np.linalg.norm(F_gravity):.2e}")
            print(f"Total load norm: {np.linalg.norm(F_total):.2e}")
            

        # Add displacement and stress debugging after load application
        if iteration == 0 and debug_level >= 1:
            # Apply boundary conditions and solve to see what displacements result
            K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
                K_global, F_total, bc_type, nodes)
            
            try:
                if hasattr(K_constrained, 'toarray'):
                    K_constrained = K_constrained.tocsr()
                displacements_free = spsolve(K_constrained, F_constrained)
                
                # Reconstruct full displacement vector
                n_dof = 2 * n_nodes
                displacements_new = np.zeros(n_dof)
                free_dofs = [i for i in range(n_dof) if i not in constraint_dofs]
                displacements_new[free_dofs] = displacements_free
                
                # Check displacement magnitudes
                max_disp = np.max(np.abs(displacements_new))
                print(f"\nFirst iteration displacement analysis:")
                print(f"  Max displacement magnitude: {max_disp:.4f}")
                print(f"  Max vertical displacement: {np.max(displacements_new[1::2]):.4f}")
                print(f"  Min vertical displacement: {np.min(displacements_new[1::2]):.4f}")
                print(f"  Max horizontal displacement: {np.max(displacements_new[0::2]):.4f}")
                print(f"  Min horizontal displacement: {np.min(displacements_new[0::2]):.4f}")
                

                    
            except Exception as e:
                print(f"  Error in displacement analysis: {e}")
        
        # Add boundary condition loads
        for i in range(n_nodes):
            if bc_type[i] == 4:  # Force boundary condition
                F_total[2*i] += bc_values[i, 0]
                F_total[2*i+1] += bc_values[i, 1]
        
        # Apply constraints using proper free DOF extraction
        K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
            K_global, F_total, bc_type, nodes)
        
        # Solve for displacement increment at equilibrium (incremental residual form)
        try:
            if hasattr(K_constrained, 'toarray'):
                K_constrained = K_constrained.tocsr()
            
            # Current free DOF vector
            n_dof = 2 * n_nodes
            free_dofs = [i for i in range(n_dof) if i not in constraint_dofs]
            u_free_curr = displacements[free_dofs]
            
            # Residual: r = F_free - K_free * u_free
            r_free = F_constrained - K_constrained @ u_free_curr
            
            # Solve K * delta_u = r
            delta_u_free = spsolve(K_constrained, r_free)
            
            # Classic update (no extra damping)
            u_free_new = u_free_curr + delta_u_free
            
            # Debug displacement increments
            if debug_level >= 2:
                max_delta_u = np.max(np.abs(delta_u_free))
                max_u_curr = np.max(np.abs(u_free_curr))
                max_u_new = np.max(np.abs(u_free_new))
                print(f"  Displacement increment: max_delta_u={max_delta_u:.2e}, max_u_curr={max_u_curr:.2e}, max_u_new={max_u_new:.2e}")
            
            # Reconstruct full displacement vector
            displacements_new = np.zeros(n_dof)
            displacements_new[free_dofs] = u_free_new
            # Constrained DOFs remain zero
        except Exception as e:
            if debug_level >= 1:
                print(f"Matrix solution failed: {e}")
            return {
                "converged": False,
                "error": f"Matrix solution failed: {e}",
                "iterations": iteration + 1,
                "displacements": displacements,
                "algorithm": "Perzyna"
            }
        
        # Update plastic strains using Perzyna algorithm with incremental approach
        plastic_strains_new, total_plastic_increment, current_stress_state = update_plastic_strains_perzyna_incremental(
            nodes, elements, element_types, element_materials,
            displacements_new, displacements_prev, plastic_strains, current_stress_state,
            c_reduced, phi_reduced, E_by_mat, nu_by_mat, dt, plastic_strain_cap, debug_level)
        
        # Debug plastic strain accumulation
        if debug_level >= 2:
            max_plastic_strain = 0.0
            total_plastic_strain = 0.0
            n_plastic_points = 0
            for elem_idx in plastic_strains_new:
                for gp in range(len(plastic_strains_new[elem_idx])):
                    plastic_magnitude = np.linalg.norm(plastic_strains_new[elem_idx][gp])
                    if plastic_magnitude > 1e-12:
                        max_plastic_strain = max(max_plastic_strain, plastic_magnitude)
                        total_plastic_strain += plastic_magnitude
                        n_plastic_points += 1
            print(f"  Plastic strain stats: max={max_plastic_strain:.2e}, total={total_plastic_strain:.2e}, n_points={n_plastic_points}")
        
        # Check convergence
        disp_change = np.linalg.norm(displacements_new - displacements)
        plastic_change = total_plastic_increment
        residual_norm = disp_change  # Using displacement change as residual
        
        # Calculate current plastic and yielding fractions
        n_plastic_gauss = 0
        n_yielding_gauss = 0
        total_gauss = 0
        for elem_idx in range(n_elements):
            elem_type = element_types[elem_idx]
            n_gauss = 4 if elem_type == 8 else 1
            total_gauss += n_gauss
            
            # Count plastic strains (use tolerance to avoid numerical precision issues)
            if elem_idx in plastic_strains:
                for gp in range(n_gauss):
                    strain_magnitude = np.linalg.norm(plastic_strains[elem_idx][gp])
                    if strain_magnitude > 1e-12:
                        n_plastic_gauss += 1
            
            # Count yielding points (F > 0)
            for gp in range(n_gauss):
                stress_gp = stress_state['element_stresses'][elem_idx, gp, :]
                f_yield = check_mohr_coulomb_cp(stress_gp, c_reduced[elem_idx], phi_reduced[elem_idx])
                if f_yield > 0:
                    n_yielding_gauss += 1
        
        plastic_fraction = n_plastic_gauss / total_gauss if total_gauss > 0 else 0
        yielding_fraction = n_yielding_gauss / total_gauss if total_gauss > 0 else 0
        
        # Diagnostic 3: Per-iteration output (controlled by frequency)
        if debug_level >= 2 and (iteration + 1) % iteration_print_frequency == 0:
            print(f"Iteration {iteration + 1}: F={F:.3f}, Residual={residual_norm:.3e}, Yielding={yielding_fraction*100:.1f}%, Plastic strains={plastic_fraction*100:.1f}%")
        elif debug_level >= 3 and (iteration + 1) % iteration_print_frequency == 0:
            print(f"Displacement change norm: {disp_change:.2e}")
            print(f"Plastic strain increment: {plastic_change:.2e}")
            print(f"Max displacement: {np.max(np.abs(displacements_new)):.6f}")
        
        # Griffiths convergence criterion - check for equilibrium
        # Converge if displacement change is small relative to current displacement magnitude
        max_current_disp = np.max(np.abs(displacements_new))
        relative_disp_change = disp_change / max(max_current_disp, 1e-6)
        
        # Don't converge too early - ensure at least some iterations for plastic development
        if iteration > 10 and relative_disp_change < tolerance and plastic_change < 0.01:
            converged = True
            if debug_level >= 2:
                print(f"Converged after {iteration + 1} iterations")
            break
        
        # Additional check: if displacements become very large, this indicates failure
        max_disp = np.max(np.abs(displacements_new))
        if max_disp > 100.0:  # Much higher threshold to allow expansive plastic zone development
            if debug_level >= 1:
                print(f"Large displacements detected ({max_disp:.3f}) - slope failure")
            converged = False
            break
        
        # Apply light numerical damping to control global instability
        displacements = damping_factor * displacements_new + (1 - damping_factor) * displacements
        displacements_prev = displacements.copy()  # Store previous displacements
        plastic_strains = plastic_strains_new
        
        # Check for early abort
        if abort_after > 0 and iteration + 1 >= abort_after:
            if debug_level >= 1:
                print(f"Aborting after iteration {iteration + 1} (abort_after={abort_after})")
            converged = True  # Mark as converged for early abort
            break
        
        # Check for excessive displacements (numerical instability indicator)
        max_disp = np.max(np.abs(displacements))
        if max_disp > 1e6:  # Much higher threshold to allow expansive plastic zone development
            if debug_level >= 1:
                print(f"Numerical instability detected: max displacement = {max_disp:.2e}")
            break
    
    # Compute final state
    final_stresses, plastic_elements = compute_final_state_perzyna(
        nodes, elements, element_types, element_materials,
        displacements, plastic_strains, c_reduced, phi_reduced,
        E_by_mat, nu_by_mat, u_nodal, current_stress_state)
    
    # Compute strains
    strains = compute_strains(nodes, elements, element_types, displacements)
    
    # Calculate final statistics
    final_sigma_yy = []
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        if elem_type == 8:
            for gp in range(4):
                final_sigma_yy.append(final_stresses[elem_idx, 1])  # Using element average stress
        else:
            final_sigma_yy.append(final_stresses[elem_idx, 1])
    
    final_min_sigma_yy = np.min(final_sigma_yy)
    final_max_sigma_yy = np.max(final_sigma_yy)
    
    # Final plastic fraction
    final_n_plastic_gauss = 0
    final_total_gauss = 0
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        n_gauss = 4 if elem_type == 8 else 1
        final_total_gauss += n_gauss
        if elem_idx in plastic_strains:
            for gp in range(n_gauss):
                strain_magnitude = np.linalg.norm(plastic_strains[elem_idx][gp])
                if strain_magnitude > 1e-12:
                    final_n_plastic_gauss += 1
    
    final_plastic_fraction = final_n_plastic_gauss / final_total_gauss if final_total_gauss > 0 else 0
    
    # Calculate final yield function values
    final_yield_function_values = np.zeros(n_elements)
    for elem_idx in range(n_elements):
        # Use the stress from final_stresses (which includes von Mises as 4th column)
        elem_stress = final_stresses[elem_idx, :3]  # [sig_x, sig_y, tau_xy]
        final_yield_function_values[elem_idx] = check_mohr_coulomb_cp(
            elem_stress, c_reduced[elem_idx], phi_reduced[elem_idx])
    
    # Diagnostic 4: Final summary
    if debug_level >= 1:
        n_yielding = np.sum(final_yield_function_values > 0)
        n_plastic = np.sum(plastic_elements)
        print(f"\n=== Final Summary ===")
        print(f"  F={F:.3f}, Iterations={iteration + 1}, Converged={'Yes' if converged else 'No'}")
        print(f"  Final residual: {residual_norm:.3e}")
        print(f"  Elements with F>0 (yielding): {n_yielding}/{n_elements}")
        print(f"  Elements with plastic strains: {n_plastic}/{n_elements} (based on final F>1e-8)")
        print(f"  Gauss points with accumulated plastic strains: {final_n_plastic_gauss}/{final_total_gauss} ({final_plastic_fraction*100:.1f}%)")
        print(f"  Min σyy: {final_min_sigma_yy:.3f} kPa, Max σyy: {final_max_sigma_yy:.3f} kPa")
    elif debug_level >= 2:
        n_plastic = np.sum(plastic_elements)
        print(f"Final: {n_plastic}/{n_elements} plastic elements")
    
    result = {
        "converged": converged,
        "iterations": iteration + 1,
        "displacements": displacements,
        "stresses": final_stresses,
        "strains": strains,
        "plastic_elements": plastic_elements,
        "yield_function": final_yield_function_values,
        "max_displacement": np.max(np.abs(displacements)),
        "plastic_strains": plastic_strains,
        "algorithm": "Perzyna Visco-Plastic",
        "F": F,
        "residual": residual_norm if 'residual_norm' in locals() else 0.0,
        "plastic_fraction": final_plastic_fraction,
        "min_sigma_yy": final_min_sigma_yy,
        "max_sigma_yy": final_max_sigma_yy
    }
    
    # Add abort information if applicable
    if abort_after > 0 and iteration + 1 >= abort_after:
        result["aborted"] = True
        result["abort_after"] = abort_after
    
    return result


def solve_ssrm(fem_data, F_min=1.0, F_max=3.0, tolerance=0.01, debug_level=0):
    """
    SSRM using Perzyna algorithm with pure non-convergence failure criterion.
    """
    
    if debug_level >= 1:
        print("=== Perzyna SSRM Analysis ===")
        print("Failure criterion: Pure non-convergence (Griffiths & Lane 1999)")
    
    F_left = F_min
    F_right = F_max
    
    # Verify bounds
    solution_min = solve_fem(fem_data, F=F_min, debug_level=max(0, debug_level-1))
    if not solution_min["converged"]:
        return {
            "converged": False,
            "error": f"F_min = {F_min} does not converge - slope unstable",
            "FS": None
        }
    
    solution_max = solve_fem(fem_data, F=F_max, debug_level=max(0, debug_level-1))
    if solution_max["converged"]:
        if debug_level >= 1:
            print(f"Warning: F_max = {F_max} still converges - very stable slope")
        return {
            "converged": True,
            "FS": F_max,
            "last_solution": solution_max,
            "note": f"Slope stable up to F = {F_max}"
        }
    
    iteration = 0
    max_iterations = 50
    last_converged_solution = solution_min  # Initialize with F_min solution
    
    # Bisection search for critical F
    while (F_right - F_left) > tolerance and iteration < max_iterations:
        F_mid = (F_left + F_right) / 2.0
        
        if debug_level >= 1:
            print(f"\nSSRM Iteration {iteration + 1}: Testing F = {F_mid:.4f}")
            print(f"Current interval: [{F_left:.4f}, {F_right:.4f}]")
        
        solution = solve_fem(fem_data, F=F_mid, debug_level=max(0, debug_level-1))
        
        if solution["converged"]:
            # F_mid is stable, critical F is higher
            F_left = F_mid
            last_converged_solution = solution
            if debug_level >= 2:
                print(f"F = {F_mid:.4f} converged (stable)")
        else:
            # F_mid failed, critical F is lower
            F_right = F_mid
            if debug_level >= 2:
                print(f"F = {F_mid:.4f} failed to converge (unstable)")
        
        iteration += 1
    
    critical_FS = F_left
    
    if debug_level >= 1:
        print(f"\nPerszyna SSRM completed: Critical FS = {critical_FS:.4f}")
        print(f"Final interval: [{F_left:.4f}, {F_right:.4f}]")
        print(f"Iterations: {iteration}")
    
    return {
        "converged": True,
        "FS": critical_FS,
        "last_solution": last_converged_solution,
        "iterations_ssrm": iteration,
        "final_interval": (F_left, F_right),
        "interval_width": F_right - F_left,
        "method": "Perzyna Visco-Plastic (Griffiths & Lane 1999)"
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


def build_quad8_stiffness_reduced_integration(coords, E, nu):
    """
    Build stiffness matrix for 8-node quadrilateral with reduced integration (4 Gauss points).
    This is the key element type used in Griffiths & Lane (1999).
    """
    # For now, use a simplified approach - proper implementation would use isoparametric mapping
    # This is a placeholder that should be replaced with full 8-node quad implementation
    
    # Simplified: use average coordinates to create an equivalent triangle
    if len(coords) >= 4:
        # Use first 4 corners for a quad approximation
        quad_coords = coords[:4]
        # Convert to equivalent triangle for now
        tri_coords = np.array([
            quad_coords[0],
            quad_coords[1], 
            quad_coords[2]
        ])
        return build_triangle_stiffness(tri_coords, E, nu)
    else:
        return build_triangle_stiffness(coords, E, nu)


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


def compute_plastic_load_correction_perzyna(nodes, elements, element_types, element_materials, 
                                           plastic_strains, E_by_mat, nu_by_mat, dt):
    """
    Compute plastic load correction vector using Perzyna algorithm.
    
    This computes the internal force vector due to plastic strains:
    F_plastic = ∫ B^T D ε_plastic dV
    """
    n_nodes = len(nodes)
    F_plastic = np.zeros(2 * n_nodes)
    
    for elem_idx in range(len(elements)):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        
        # Get element data
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        if elem_type == 8:
            # 8-node quad with 2x2 Gauss integration
            gauss_points, weights = get_gauss_points_2x2()
            n_gauss = 4
        else:
            n_gauss = 1  # Triangle - single Gauss point
            gauss_points = [(0.0, 0.0)]
            weights = [1.0]
        
        # Element plastic force vector
        elem_f_plastic = np.zeros(2 * elem_type)
        
        # For each Gauss point
        for gp in range(n_gauss):
            # Get plastic strains at this Gauss point
            plastic_strain = plastic_strains[elem_idx][gp, :]
            
            # Skip if no plastic strain
            if np.linalg.norm(plastic_strain) < 1e-20:
                continue
            
            # Constitutive matrix
            D = build_constitutive_matrix(E, nu)
            plastic_stress = D @ plastic_strain  # Tension-positive convention
            
            # Compute B and weight at this Gauss point
            if elem_type == 3:  # Triangle
                B, area = compute_B_matrix_triangle(elem_coords)
                weight = area
            elif elem_type == 8:  # 8-node quad
                xi, eta_local = gauss_points[gp]
                # Build B and detJ at Gauss point
                dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta_local)
                # Jacobian
                J = np.zeros((2, 2))
                for a in range(8):
                    J[0,0] += dN_dxi[a] * elem_coords[a,0]
                    J[0,1] += dN_dxi[a] * elem_coords[a,1]
                    J[1,0] += dN_deta[a] * elem_coords[a,0]
                    J[1,1] += dN_deta[a] * elem_coords[a,1]
                det_J = J[0,0]*J[1,1] - J[0,1]*J[1,0]
                if abs(det_J) < 1e-14:
                    continue
                J_inv = np.array([[J[1,1], -J[0,1]], [-J[1,0], J[0,0]]]) / det_J
                dN_dx = np.zeros(8)
                dN_dy = np.zeros(8)
                for a in range(8):
                    dN_dx[a] = J_inv[0,0]*dN_dxi[a] + J_inv[0,1]*dN_deta[a]
                    dN_dy[a] = J_inv[1,0]*dN_dxi[a] + J_inv[1,1]*dN_deta[a]
                B = np.zeros((3, 16))
                for a in range(8):
                    B[0, 2*a]   = dN_dx[a]
                    B[1, 2*a+1] = dN_dy[a]
                    B[2, 2*a]   = dN_dy[a]
                    B[2, 2*a+1] = dN_dx[a]
                weight = weights[gp] * abs(det_J)
            else:
                # Simplified for other elements
                B = np.zeros((3, 2 * elem_type))
                weight = 1.0
            
            # Add contribution to element force vector
            if B.size > 0:
                elem_f_plastic += B.T @ plastic_stress * weight
        
        # Assemble into global force vector
        for i in range(elem_type):
            node = elem_nodes[i]
            F_plastic[2*node] += elem_f_plastic[2*i]
            F_plastic[2*node + 1] += elem_f_plastic[2*i + 1]
    
    return F_plastic


def compute_B_matrix_quad8_centroid(coords):
    """Compute B matrix and determinant of Jacobian for 8-node quad at centroid."""
    # Evaluate at centroid (xi=0, eta=0)
    xi, eta = 0.0, 0.0
    
    # Shape function derivatives
    dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)
    
    # Jacobian matrix
    J = np.zeros((2, 2))
    for i in range(8):
        x, y = coords[i]
        J[0, 0] += dN_dxi[i] * x   # dx/dxi
        J[0, 1] += dN_dxi[i] * y   # dy/dxi  
        J[1, 0] += dN_deta[i] * x  # dx/deta
        J[1, 1] += dN_deta[i] * y  # dy/deta
    
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    
    if abs(det_J) < 1e-12:
        return np.zeros((3, 16)), 0.0
    
    # Inverse Jacobian
    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J
    
    # Shape function derivatives in physical coordinates
    dN_dx = np.zeros(8)
    dN_dy = np.zeros(8)
    for i in range(8):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    # B matrix (standard tension positive)
    B = np.zeros((3, 16))  # 3 strains x 16 DOFs (8 nodes x 2 DOFs)
    for i in range(8):
        B[0, 2*i] = dN_dx[i]      # εx = ∂u/∂x
        B[1, 2*i+1] = dN_dy[i]    # εy = ∂v/∂y
        B[2, 2*i] = dN_dy[i]      # γxy = ∂u/∂y + ∂v/∂x
        B[2, 2*i+1] = dN_dx[i]    # γxy = ∂u/∂y + ∂v/∂x
    
    return B, abs(det_J)


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


def check_initial_yield_state(stress_state, c_values, phi_values):
    """
    Check how many elements are yielding at the initial stress state.
    
    Args:
        stress_state: Dictionary with 'element_stresses' array
        c_values: Cohesion values for each element
        phi_values: Friction angle values for each element (in radians)
        
    Returns:
        Number of elements that are yielding
    """
    element_stresses = stress_state['element_stresses']
    n_elements = element_stresses.shape[0]
    yield_count = 0
    
    for elem_idx in range(n_elements):
        # Check first Gauss point of each element
        sig_x = element_stresses[elem_idx, 0, 0]
        sig_y = element_stresses[elem_idx, 0, 1]
        tau_xy = element_stresses[elem_idx, 0, 2]
        
        c = c_values[elem_idx]
        phi = phi_values[elem_idx]
        
        # Check Mohr-Coulomb yield criterion
        stress = np.array([sig_x, sig_y, tau_xy])
        F_yield = check_mohr_coulomb_cp_from_tp(stress, c, phi)
        if F_yield > 0:
            yield_count += 1
    
    return yield_count


def update_plastic_strains_perzyna(nodes, elements, element_types, element_materials,
                                  displacements, plastic_strains, c_reduced, phi_reduced,
                                  E_by_mat, nu_by_mat, u_nodal, eta, initial_stresses=None):
    """
    Update plastic strains using Perzyna visco-plastic algorithm with proper Gauss integration.
    """
    plastic_strains_new = {}
    total_increment = 0.0
    
    # Get Gauss points for 8-node quads
    gauss_points_2x2, weights_2x2 = get_gauss_points_2x2()
    
    for elem_idx in range(len(elements)):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        c = c_reduced[elem_idx]
        phi = phi_reduced[elem_idx]
        
        if elem_type == 8:
            n_gauss = 4  # 8-node quad with reduced integration
        else:
            n_gauss = 1  # Triangle - single Gauss point
        
        plastic_strains_new[elem_idx] = plastic_strains[elem_idx].copy()
        
        # Get element data
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Get element displacements
        elem_disp = np.zeros(2 * elem_type)
        for i, node in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node]
            elem_disp[2*i+1] = displacements[2*node+1]
        
        # For each Gauss point
        for gp in range(n_gauss):
            # Compute total strains at this Gauss point
            if elem_type == 3:  # Triangle
                total_strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            elif elem_type == 8:  # 8-node quad with proper Gauss points
                xi, eta_local = gauss_points_2x2[gp]
                total_strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta_local)
            else:
                # Simplified for other element types
                total_strains = np.array([0.0, 0.0, 0.0])
            
            # Elastic trial strains
            plastic_strain_old = plastic_strains[elem_idx][gp, :]
            elastic_strains = total_strains - plastic_strain_old
            
            # Elastic trial stress = initial stress + incremental stress
            D = build_constitutive_matrix(E, nu)
            incremental_stress = D @ elastic_strains  # Tension-positive convention
            
            # Add initial stress if provided
            if initial_stresses is not None:
                initial_stress = initial_stresses['element_stresses'][elem_idx, gp, :]
                trial_stress = initial_stress + incremental_stress
            else:
                trial_stress = incremental_stress
            
            # Check yield criterion with total stress
            f_yield = check_mohr_coulomb_cp_from_tp(trial_stress, c, phi)
            
            if f_yield > 1e-6:  # Plastic loading - higher threshold to reduce initial yielding
                # Perzyna visco-plastic flow as per Griffiths & Lane (1999)
                # Δλ = η * <f> where <f> = max(0, f)
                # Use appropriate viscosity parameter from paper
                delta_lambda = eta * f_yield  # Remove artificial cap
                
                # Flow vector (non-associated: ψ = 0)
                flow_vector = compute_plastic_flow_vector_cp_return_tp(trial_stress, 0.0)  # ψ = 0
                
                # Plastic strain increment - controlled flow to prevent instability
                plastic_increment = delta_lambda * flow_vector
                
                # Apply reasonable limit to prevent numerical explosion
                increment_norm = np.linalg.norm(plastic_increment)
                if increment_norm > 1e-5:  # Much smaller limit for very controlled plastic development
                    plastic_increment *= 1e-5 / increment_norm
                
                # Update plastic strains
                plastic_strains_new[elem_idx][gp, :] += plastic_increment
                
                # Track total plastic increment
                total_increment += np.linalg.norm(plastic_increment)
    
    return plastic_strains_new, total_increment


def update_plastic_strains_perzyna_incremental(nodes, elements, element_types, element_materials,
                                              displacements_new, displacements_prev, plastic_strains, 
                                              current_stress_state, c_reduced, phi_reduced,
                                              E_by_mat, nu_by_mat, dt, plastic_strain_cap, debug_level=0):
    """
    Update plastic strains using incremental Perzyna algorithm with compression-positive stress storage.
    """
    plastic_strains_new = {}
    total_increment = 0.0
    
    # New stress state (compression-positive)
    new_stress_state = {
        'element_stresses': current_stress_state['element_stresses'].copy(),
        'plastic_state': current_stress_state['plastic_state'].copy()
    }
    
    # Displacement increment
    displacement_increment = displacements_new - displacements_prev
    
    for elem_idx in range(len(elements)):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        c = c_reduced[elem_idx]
        phi = phi_reduced[elem_idx]
        
        # Element nodes/coords
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Incremental displacements for element
        elem_disp_increment = np.zeros(2 * elem_type)
        for i, node_idx in enumerate(elem_nodes):
            elem_disp_increment[2*i] = displacement_increment[2*node_idx]
            elem_disp_increment[2*i+1] = displacement_increment[2*node_idx+1]
        
        plastic_strains_new[elem_idx] = plastic_strains[elem_idx].copy()
        
        # Gauss points
        if elem_type == 8:
            gauss_points_2x2, _ = get_gauss_points_2x2()
            n_gauss = 4
        else:
            n_gauss = 1
        
        for gp in range(n_gauss):
            # Incremental strains
            if elem_type == 3:
                incremental_strains = compute_triangle_strains_manual(elem_coords, elem_disp_increment)
            elif elem_type == 8:
                xi, eta_local = gauss_points_2x2[gp]
                incremental_strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp_increment, xi, eta_local)
            else:
                incremental_strains = np.array([0.0, 0.0, 0.0])
            
            D = build_constitutive_matrix(E, nu)
            incremental_stress_tp = D @ incremental_strains
            incremental_stress_cp = -incremental_stress_tp
            
            # Current total stress at GP (compression-positive)
            current_stress_cp = current_stress_state['element_stresses'][elem_idx, gp, :]
            
            # Trial stress (compression-positive)
            trial_stress_cp = current_stress_cp + incremental_stress_cp
            
            # Yield check (compression-positive)
            f_yield = check_mohr_coulomb_cp(trial_stress_cp, c, phi)
            
            if f_yield > 1e-12:
                # Flow direction in tp from cp stress
                n_flow_tp = compute_plastic_flow_vector_cp_return_tp(-trial_stress_cp, 0.0)
                
                # Simple stress return: use Perzyna approach with controlled plastic multiplier
                D = build_constitutive_matrix(E, nu)
                
                # Griffiths & Lane viscoplastic strain method with strain softening
                # Calculate plastic strain rate using their approach
                if f_yield > 0:
                    # Non-associated flow with ψ=0; Perzyna rate: erate ∝ f * n
                    flow_vector = f_yield * n_flow_tp
                    erate = flow_vector
                    # Plastic strain increment: Δεp = erate * dt
                    plastic_increment = erate * dt
                else:
                    plastic_increment = np.zeros(3)
                inc_norm = float(np.linalg.norm(plastic_increment))
                
                # Debug stress return
                # if debug_level >= 3 and inc_norm > 1e-3:
                #     print(f"    GP {gp}: inc_norm={inc_norm:.2e}, f_yield={f_yield:.2e}")
                
                # Enforce plastic_strain_cap per Gauss point
                if plastic_strain_cap is not None and inc_norm > plastic_strain_cap:
                    plastic_increment *= plastic_strain_cap / max(inc_norm, 1e-20)
                    inc_norm = plastic_strain_cap
                
                # Update plastic strains
                plastic_strains_new[elem_idx][gp, :] += plastic_increment
                
                # Update stress state (remove plastic stress contribution)
                plastic_stress_tp = D @ plastic_increment
                plastic_stress_cp = -plastic_stress_tp
                new_stress_state['element_stresses'][elem_idx, gp, :] = trial_stress_cp + plastic_stress_cp
                total_increment += inc_norm
            else:
                new_stress_state['element_stresses'][elem_idx, gp, :] = trial_stress_cp
    
    return plastic_strains_new, total_increment, new_stress_state


def compute_final_state_perzyna(nodes, elements, element_types, element_materials,
                               displacements, plastic_strains, c_reduced, phi_reduced,
                               E_by_mat, nu_by_mat, u_nodal, stress_state):
    """
    Compute final stress state (compression-positive) and identify plastic elements.
    """
    n_elements = len(elements)
    final_stresses = np.zeros((n_elements, 4))
    plastic_elements = np.zeros(n_elements, dtype=bool)
    
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        c = c_reduced[elem_idx]
        phi = phi_reduced[elem_idx]
        
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        elem_disp = np.zeros(2 * elem_type)
        for i, node_idx in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node_idx]
            elem_disp[2*i+1] = displacements[2*node_idx+1]
        
        if elem_type == 3:
            total_strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            plastic_strain = plastic_strains[elem_idx][0, :] if elem_idx in plastic_strains else np.zeros(3)
            elastic_strains = total_strains - plastic_strain
            D = build_constitutive_matrix(E, nu)
            stress_tp = D @ elastic_strains
            if stress_state is not None:
                initial_cp = stress_state['element_stresses'][elem_idx, 0, :]
                stress_cp = initial_cp - stress_tp  # cp = initial_cp + (-tp)
            else:
                stress_cp = -stress_tp
            sig_x, sig_y, tau_xy = stress_cp
            sig_vm = np.sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
            final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]
            f_yield = check_mohr_coulomb_cp(stress_cp, c, phi)
            plastic_elements[elem_idx] = f_yield > 1e-8
        else:
            # 8-node quad: average cp stress over Gauss points from stress_state
            elem_stress_avg_cp = np.zeros(3)
            n_gauss = 4
            for gp in range(n_gauss):
                elem_stress_avg_cp += stress_state['element_stresses'][elem_idx, gp, :]
            elem_stress_avg_cp /= n_gauss
            sig_x, sig_y, tau_xy = elem_stress_avg_cp
            sig_vm = np.sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
            final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]
            f_yield = check_mohr_coulomb_cp(elem_stress_avg_cp, c, phi)
            plastic_elements[elem_idx] = f_yield > 1e-8
    
    return final_stresses, plastic_elements


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


def check_shear_convention_consistency(strain_vec, D_matrix, E, nu, element_id=None, debug=True):
    """
    Diagnostic check for engineering vs tensor shear strain convention consistency.
    
    This catches the common bug where B-matrix computes εxy (tensor shear) but 
    D-matrix expects γxy (engineering shear), causing τxy to be off by factor of 2.
    
    Args:
        strain_vec: [εx, εy, γxy] strain vector from B-matrix
        D_matrix: Constitutive matrix 
        E, nu: Material properties
        element_id: For debugging output
        debug: Whether to print warnings
    """
    if len(strain_vec) < 3:
        return True  # Skip if not enough components
        
    ex, ey, gxy = strain_vec[:3]  # what your code uses
    
    # Expected shear modulus for engineering shear
    mu = E / (2 * (1 + nu))
    
    # What D-matrix actually computes for shear stress
    tau_from_D = (D_matrix @ strain_vec)[2]
    
    # What we expect for engineering shear: τxy = μ * γxy  
    tau_expected = mu * gxy
    
    # Check consistency
    tolerance = 1e-6 * max(1.0, abs(tau_expected))
    is_consistent = abs(tau_from_D - tau_expected) <= tolerance
    
    if not is_consistent and debug:
        elem_str = f" (element {element_id})" if element_id is not None else ""
        print(f"WARNING: Shear convention mismatch{elem_str}!")
        print(f"  γxy from B-matrix: {gxy:.6f}")
        print(f"  τxy from D-matrix: {tau_from_D:.2f}")  
        print(f"  τxy expected (μ*γxy): {tau_expected:.2f}")
        print(f"  Ratio (actual/expected): {tau_from_D/tau_expected if abs(tau_expected) > 1e-12 else 'inf':.3f}")
        print(f"  D[2,2] = {D_matrix[2,2]:.1f}, expected μ = {mu:.1f}")
        print("  Either B builds [εx, εy, εxy] but labels it γxy, or")
        print("  D[2,2] isn't correct for engineering shear convention")
    
    return is_consistent


def validate_bd_matrices_simple_test():
    """
    Simple validation test for B and D matrix consistency using pure shear.
    Apply known displacement field and check if resulting stresses are correct.
    """
    print("\n=== B/D Matrix Validation Test ===")
    
    # Simple triangle for testing
    coords = np.array([[0, 0], [1, 0], [0, 1]])  # Right triangle
    E, nu = 30000, 0.3
    
    # Test case 1: Pure shear displacement field
    # u = 0.001*y, v = 0.001*x (creates γxy = 0.002)
    displacements = np.array([
        0.0,     0.0,      # Node 1: u1=0, v1=0
        0.0,     0.001,    # Node 2: u2=0, v2=0.001  
        0.001,   0.0       # Node 3: u3=0.001, v3=0
    ])
    
    # Compute strains using triangle B-matrix
    strains = compute_triangle_strains_manual(coords, displacements)
    print(f"Computed strains: εx={strains[0]:.6f}, εy={strains[1]:.6f}, γxy={strains[2]:.6f}")
    
    # Expected: εx=0, εy=0, γxy=0.002 (engineering shear)
    expected_gxy = 0.002
    print(f"Expected γxy: {expected_gxy:.6f}")
    
    # Compute stresses  
    D = build_constitutive_matrix(E, nu)
    stresses = D @ strains
    print(f"Computed stresses: σx={stresses[0]:.2f}, σy={stresses[1]:.2f}, τxy={stresses[2]:.2f}")
    
    # Expected: σx=0, σy=0, τxy = G*γxy = E/(2*(1+ν))*γxy
    G = E / (2 * (1 + nu))
    expected_tau_xy = G * expected_gxy
    print(f"Expected τxy: G*γxy = {G:.1f} * {expected_gxy:.6f} = {expected_tau_xy:.2f}")
    
    # Check consistency
    is_consistent = check_shear_convention_consistency(strains, D, E, nu, debug=False)
    
    if is_consistent:
        print("✓ B/D matrices are consistent (engineering shear convention)")
    else:
        print("✗ B/D matrices have shear convention mismatch!")
    
    print("=================================\n")
    return is_consistent

# Uncomment this line to run the validation test:
# validate_bd_matrices_simple_test()


def test_constitutive_matrix_sanity():
    """
    Basic sanity checks for constitutive matrix behavior.
    Tests fundamental stress-strain relationships.
    """
    print("\n=== Constitutive Matrix Sanity Tests ===")
    
    E, nu = 1e5, 0.3
    D = build_constitutive_matrix(E, nu)
    
    # Test (A) Pure vertical shortening → σyy negative
    print("Test A: Pure vertical shortening")
    eps = np.array([0.0, -1e-4, 0.0])   # [εx, εy, γxy]
    sig = D @ eps
    print(f"  Strain: εx={eps[0]:.1e}, εy={eps[1]:.1e}, γxy={eps[2]:.1e}")
    print(f"  Stress: σx={sig[0]:.1f}, σy={sig[1]:.1f}, τxy={sig[2]:.1f}")
    
    try:
        assert sig[1] < 0, "σyy must be negative under vertical shortening (tension-positive)."
        print("  ✓ PASS: σyy is negative (compression) as expected")
    except AssertionError as e:
        print(f"  ✗ FAIL: {e}")
        print(f"    σyy = {sig[1]:.3f}, expected < 0")
    
    # Test (B) Pure shear → τxy = μ * γxy
    print("\nTest B: Pure shear stress")
    mu = E / (2*(1+nu))
    eps = np.array([0.0, 0.0, 1e-3])    # γxy = 1e-3
    sig = D @ eps
    expected_tau = mu * eps[2]
    
    print(f"  Strain: εx={eps[0]:.1e}, εy={eps[1]:.1e}, γxy={eps[2]:.1e}")
    print(f"  Stress: σx={sig[0]:.1f}, σy={sig[1]:.1f}, τxy={sig[2]:.1f}")
    print(f"  Expected τxy = μ*γxy = {mu:.1f} * {eps[2]:.1e} = {expected_tau:.1f}")
    
    try:
        assert abs(sig[2] - mu*eps[2]) < 1e-8*max(1.0, abs(mu*eps[2])), "τxy must equal μ*γxy."
        print("  ✓ PASS: τxy = μ*γxy exactly")
    except AssertionError as e:
        print(f"  ✗ FAIL: {e}")
        print(f"    τxy = {sig[2]:.6f}, expected = {expected_tau:.6f}")
        print(f"    Error = {abs(sig[2] - expected_tau):.2e}")
    
    # Test (C) Check Poisson effect - verify D-matrix formulation
    print("\nTest C: Poisson effect under uniaxial tension")
    eps = np.array([1e-4, 0.0, 0.0])   # Pure εx
    sig = D @ eps
    
    print(f"  Strain: εx={eps[0]:.1e}, εy={eps[1]:.1e}, γxy={eps[2]:.1e}")
    print(f"  Stress: σx={sig[0]:.1f}, σy={sig[1]:.1f}, τxy={sig[2]:.1f}")
    
    # Debug: Check what D-matrix actually looks like
    print(f"  D-matrix:")
    print(f"    [{D[0,0]:.1f}  {D[0,1]:.1f}  {D[0,2]:.1f}]")
    print(f"    [{D[1,0]:.1f}  {D[1,1]:.1f}  {D[1,2]:.1f}]")
    print(f"    [{D[2,0]:.1f}  {D[2,1]:.1f}  {D[2,2]:.1f}]")
    
    # What we expect from your D-matrix formulation
    factor = E / ((1 + nu) * (1 - 2*nu))
    expected_d11 = factor * (1 - nu)
    expected_d12 = factor * nu  
    expected_d22 = factor * (1 - nu)
    expected_d33 = factor * (1 - 2*nu) / 2
    
    print(f"  Expected D-matrix (your formulation):")
    print(f"    [{expected_d11:.1f}  {expected_d12:.1f}  0.0]")
    print(f"    [{expected_d12:.1f}  {expected_d22:.1f}  0.0]")
    print(f"    [0.0  0.0  {expected_d33:.1f}]")
    
    # For your D-matrix, under pure εx strain:
    expected_sigx_your = expected_d11 * eps[0]  # D[0,0] * εx
    expected_sigy_your = expected_d12 * eps[0]  # D[1,0] * εx
    
    print(f"  Your D gives: σx = {expected_sigx_your:.1f}, σy = {expected_sigy_your:.1f}")
    
    # Verify this matches what D actually computed
    if abs(sig[0] - expected_sigx_your) < 1e-10 and abs(sig[1] - expected_sigy_your) < 1e-10:
        print("  ✓ PASS: D-matrix working as designed")
    else:
        print("  ✗ FAIL: D-matrix computation error")
    
    print("==========================================\n")


# Run the tests
# test_constitutive_matrix_sanity()


def check_mohr_coulomb_cp_from_tp(stress_tp, c, phi):
    """Mohr-Coulomb yield using compression-positive convention, input stress in tension-positive.
    F = tau_max - sigma_mean_cp * sin(phi) - c * cos(phi)
    Positive F => yield.
    """
    sig_x_cp, sig_y_cp, tau_xy = stress_tp_to_cp(stress_tp)
    sig_mean_cp = (sig_x_cp + sig_y_cp) / 2.0
    tau_max = sqrt(((sig_x_cp - sig_y_cp) / 2.0)**2 + tau_xy**2)
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    F = tau_max - sig_mean_cp * sin_phi - c * cos_phi
    return F


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


def compute_mohr_coulomb_potential_derivatives(psi_deg, dsbar, theta):
    """
    Compute derivatives of Mohr-Coulomb potential function with respect to invariants.
    Based on mocouq.f90 from Griffiths & Lane implementation.
    
    Args:
        psi_deg: Dilation angle in degrees
        dsbar: Second deviatoric stress invariant
        theta: Lode angle in radians
        
    Returns:
        dq1, dq2, dq3: Derivatives with respect to I1, sqrt(J2), theta
    """
    psi_rad = np.radians(psi_deg)
    sin_theta = np.sin(theta)
    sin_psi = np.sin(psi_rad)
    
    dq1 = sin_psi
    
    if abs(sin_theta) > 0.49:
        c1 = 1.0 if sin_theta >= 0 else -1.0
        dq2 = (np.sqrt(3) * 0.5 - c1 * sin_psi * 0.5 / np.sqrt(3)) * np.sqrt(3) * 0.5 / dsbar
    else:
        cos_theta = np.cos(theta)
        cos_3theta = np.cos(3 * theta)
        tan_3theta = np.tan(3 * theta)
        dq2 = (np.sqrt(3) * 0.5 * cos_theta + sin_psi * (sin_theta + sin_theta * cos_3theta / (3 * cos_3theta))) / (dsbar * cos_theta)
    
    dq3 = 0.0  # Simplified - could include more complex terms
    
    return dq1, dq2, dq3

def compute_plastic_flow_vector_cp_return_tp(stress_tp, psi):
    """Compute flow direction for compression-positive potential and return vector in tension-positive axes.
    Uses g = tau_max - sigma_mean_cp * sin(psi).
    Maps derivatives back to tension-positive components: d/dsig_tp = - d/dsig_cp for normal components.
    """
    sig_x_cp, sig_y_cp, tau_xy = stress_tp_to_cp(stress_tp)
    sig_mean_cp = (sig_x_cp + sig_y_cp) / 2.0
    tau_max = sqrt(((sig_x_cp - sig_y_cp) / 2.0)**2 + tau_xy**2)
    if tau_max < 1e-20:
        return np.array([0.0, 0.0, 0.0])
    # Derivatives in cp convention
    dsig_mean_dsigx_cp = 0.5
    dsig_mean_dsigy_cp = 0.5
    dsig_mean_dtau_cp = 0.0
    dtau_dsigx_cp = (sig_x_cp - sig_y_cp) / (4.0 * tau_max)
    dtau_dsigy_cp = -(sig_x_cp - sig_y_cp) / (4.0 * tau_max)
    dtau_dtau_cp   = tau_xy / tau_max
    sin_psi = sin(psi)
    # ∂g/∂σ_cp = ∂τ * 1 + ∂σ_mean * ( - sin ψ)
    flow_x_cp = dtau_dsigx_cp - sin_psi * dsig_mean_dsigx_cp
    flow_y_cp = dtau_dsigy_cp - sin_psi * dsig_mean_dsigy_cp
    flow_xy_cp = dtau_dtau_cp - sin_psi * dsig_mean_dtau_cp
    # Map back to tension-positive: d/dsig_tp = - d/dsig_cp for normals; shear unchanged
    flow_x_tp = -flow_x_cp
    flow_y_tp = -flow_y_cp
    flow_xy_tp = flow_xy_cp
    return np.array([flow_x_tp, flow_y_tp, flow_xy_tp])


def stress_tp_to_cp(stress_tp):
    """Convert tension-positive stress [sigx, sigy, tau_xy] to compression-positive."""
    sig_x, sig_y, tau_xy = stress_tp
    return np.array([-sig_x, -sig_y, tau_xy])


def establish_k0_stress_state(K_global, F_gravity, bc_type, nodes, elements, element_types, 
                             element_materials, E_by_mat, nu_by_mat, gamma_by_mat, u_nodal, debug_level=0):
    """
    Establish K₀ initial stress state through elastic gravity loading.
    
    This creates the geostatic stress field that exists before applying strength reduction.
    Critical for developing proper rotational failure modes in slopes.
    """
    
    # Apply boundary conditions to gravity loading system
# apply_boundary_conditions is now defined in this same file
    K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
        K_global, F_gravity, bc_type, nodes)
    
    # Solve elastic system under gravity
    try:
        if hasattr(K_constrained, 'toarray'):
            K_constrained = K_constrained.tocsr()
        displacements_free = spsolve(K_constrained, F_constrained)
        
        # Reconstruct full displacement vector
        n_dof = 2 * len(nodes)
        displacements = np.zeros(n_dof)
        free_dofs = [i for i in range(n_dof) if i not in constraint_dofs]
        displacements[free_dofs] = displacements_free
        
    except Exception as e:
        print(f"K₀ stress establishment failed: {e}")
        # Fall back to zero displacement
        displacements = np.zeros(2 * len(nodes))
    
    # Compute stress state from elastic solution
    stress_state = compute_k0_stress_state(
        nodes, elements, element_types, element_materials, displacements,
        E_by_mat, nu_by_mat, gamma_by_mat, u_nodal)
    
    if debug_level >= 2:
        max_disp = np.max(np.abs(displacements))
        print(f"  K₀ solution: max displacement = {max_disp:.6f}")
        
        # Debug: Check actual displacement at a specific node
        node_near_top = nodes[:, 1].argmax()  # Node with highest y coordinate
        disp_x = displacements[2*node_near_top]
        disp_y = displacements[2*node_near_top+1]
        print(f"  Top node {node_near_top} at y={nodes[node_near_top, 1]:.1f}: disp_x={disp_x:.6f}, disp_y={disp_y:.6f}")
        
        n_stress_elements = len(stress_state.get('element_stresses', []))
        print(f"  Stress state established for {n_stress_elements} elements")
    
    return displacements, stress_state


def compute_k0_stress_state(nodes, elements, element_types, element_materials, displacements,
                           E_by_mat, nu_by_mat, gamma_by_mat, u_nodal):
    """
    Compute initial stress state from elastic FEM gravity solution.
    
    Following Griffiths & Lane (1999): "The present work applies gravity in a single 
    increment to an initially stress-free slope" - this means:
    1. Start with zero stress everywhere
    2. Apply gravity loads via FEM
    3. Compute strains from resulting displacements  
    4. Compute stresses from elastic strains: σ = D·ε
    
    Store stresses as compression-positive throughout the codebase.
    """
    n_elements = len(elements)
    max_gauss_points = 4
    element_stresses = np.zeros((n_elements, max_gauss_points, 3))  # [sig_x, sig_y, tau_xy] (compression+)
    
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        
        # Get element data
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Get Gauss points for proper integration
        if elem_type == 8:
            gauss_points_2x2, _ = get_gauss_points_2x2()
            n_gauss = 4
        else:
            n_gauss = 1
        
        # Get element displacements from FEM gravity solution
        elem_disp = np.zeros(2 * elem_type)
        for i, node in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node]
            elem_disp[2*i+1] = displacements[2*node+1]
        
        # Compute stresses at each Gauss point from FEM strains
        for gp in range(n_gauss):
            if elem_type == 3:  # Triangle
                strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            elif elem_type == 8:  # 8-node quad
                xi, eta_local = gauss_points_2x2[gp]
                strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta_local)
            else:
                strains = np.array([0.0, 0.0, 0.0])
            
            # Compute stresses from strains using elastic constitutive matrix
            D = build_constitutive_matrix(E, nu)
            stresses_tp = D @ strains  # tension-positive
            stresses_cp = -stresses_tp  # store compression-positive
            
            element_stresses[elem_idx, gp, :] = stresses_cp
    
    # Debug: Check stress state statistics (compression-positive)
    stress_stats = {
        'sigma_x': {'min': np.min(element_stresses[:, :, 0]), 'max': np.max(element_stresses[:, :, 0]), 'mean': np.mean(element_stresses[:, :, 0])},
        'sigma_y': {'min': np.min(element_stresses[:, :, 1]), 'max': np.max(element_stresses[:, :, 1]), 'mean': np.mean(element_stresses[:, :, 1])},
        'tau_xy': {'min': np.min(element_stresses[:, :, 2]), 'max': np.max(element_stresses[:, :, 2]), 'mean': np.mean(element_stresses[:, :, 2])}
    }
    
    print(f"  Initial stress state statistics (compression+):")
    print(f"    σ_x: min={stress_stats['sigma_x']['min']:.1f}, max={stress_stats['sigma_x']['max']:.1f}, mean={stress_stats['sigma_x']['mean']:.1f}")
    print(f"    σ_y: min={stress_stats['sigma_y']['min']:.1f}, max={stress_stats['sigma_y']['max']:.1f}, mean={stress_stats['sigma_y']['mean']:.1f}")
    print(f"    τ_xy: min={stress_stats['tau_xy']['min']:.1f}, max={stress_stats['tau_xy']['max']:.1f}, mean={stress_stats['tau_xy']['mean']:.1f}")
    
    return {
        'element_stresses': element_stresses,
        'plastic_state': np.zeros((n_elements, max_gauss_points), dtype=bool)
    }


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

def compute_simple_quad4_strains(coords, displacements):
    """
    Simple strain calculation for 4-node quad using bilinear interpolation.
    This is a test to see if the issue is in the isoparametric formulation.
    """
    # Use center point (xi=0, eta=0) for simplicity
    xi, eta = 0.0, 0.0
    
    # 4-node bilinear shape function derivatives
    dN_dxi = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
    dN_deta = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)])
    
    # Jacobian matrix
    J = np.zeros((2, 2))
    for i in range(4):
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
    dN_dx = np.zeros(4)
    dN_dy = np.zeros(4)
    for i in range(4):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    # B matrix (standard tension positive, 3 strains x 8 DOFs for 4 nodes)
    B = np.zeros((3, 8))
    for i in range(4):
        B[0, 2*i] = dN_dx[i]      # εx = ∂u/∂x
        B[1, 2*i+1] = dN_dy[i]    # εy = ∂v/∂y
        B[2, 2*i] = dN_dy[i]      # γxy = ∂u/∂y + ∂v/∂x
        B[2, 2*i+1] = dN_dx[i]    # γxy = ∂u/∂y + ∂v/∂x
    
    # Compute strains
    strains = B @ displacements
    return strains

def compute_quad8_strains_at_gauss_point(coords, displacements, gauss_point):
    """
    Compute strains for 8-node quadrilateral element at specific Gauss point.
    
    This implements the exact formulation used in Griffiths & Lane (1999).
    Uses reduced integration with 4 Gauss points (2x2 rule).
    """
    # 2x2 Gauss points for reduced integration (as per Griffiths paper)
    gauss_coords = [
        (-0.5773502692, -0.5773502692),  # Point 0
        ( 0.5773502692, -0.5773502692),  # Point 1
        ( 0.5773502692,  0.5773502692),  # Point 2
        (-0.5773502692,  0.5773502692)   # Point 3
    ]
    
    if gauss_point >= len(gauss_coords):
        return np.array([0.0, 0.0, 0.0])
    
    xi, eta = gauss_coords[gauss_point]
    
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
    J_inv = np.array([[J[1, 1], -J[0, 1]], 
                      [-J[1, 0], J[0, 0]]]) / det_J
    
    # Shape function derivatives in physical coordinates
    dN_dx = np.zeros(8)
    dN_dy = np.zeros(8)
    for i in range(8):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    # B matrix for strain calculation
    B = np.zeros((3, 16))  # 3 strains x 16 DOFs (8 nodes x 2 DOFs)
    for i in range(8):
        B[0, 2*i]     = dN_dx[i]    # ∂u/∂x
        B[1, 2*i+1]   = dN_dy[i]    # ∂v/∂y  
        B[2, 2*i]     = dN_dy[i]    # ∂u/∂y
        B[2, 2*i+1]   = dN_dx[i]    # ∂v/∂x
    
    # Compute strains: ε = B * u
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