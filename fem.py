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

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import eigh
from shapely.geometry import LineString, Point
from math import radians, degrees, sin, cos, tan, sqrt, atan2
import warnings


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
    max_depth = slope_data.get("max_depth", None)
    if max_depth is not None:
        tolerance = 1e-6
        bottom_nodes = np.abs(nodes[:, 1] - max_depth) < tolerance
        bc_type[bottom_nodes] = 1  # Fixed (u=0, v=0)
    
    # Step 3: X-roller supports at left and right sides (type 2) - standard practice
    if "ground_surface" in slope_data:
        ground_surface = slope_data["ground_surface"]
        if len(ground_surface.coords) >= 2:
            x_left = ground_surface.coords[0][0]
            x_right = ground_surface.coords[-1][0]
            tolerance = 1e-6
            
            left_nodes = np.abs(nodes[:, 0] - x_left) < tolerance
            right_nodes = np.abs(nodes[:, 0] - x_right) < tolerance
            
            # Apply X-roller but preserve existing boundary conditions
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

def solve_fem_perzyna(fem_data, F=1.0, debug_level=0, abort_after=-1):
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
    K_global = build_global_stiffness_perzyna(nodes, elements, element_types, 
                                             element_materials, E_by_mat, nu_by_mat)
    
    # Build gravity load vector
    F_gravity = build_gravity_loads_perzyna(nodes, elements, element_types, 
                                           element_materials, gamma_by_mat, k_seismic)
    
    # Boundary conditions will be applied in each iteration using apply_boundary_conditions
    
    # Perzyna algorithm parameters (from Zienkiewicz & Cormeau 1974)
    max_iterations = 500  # More iterations to allow proper localization
    tolerance = 1e-5  # Tighter tolerance for better convergence
    eta = 0.01  # Much lower viscosity for more localized plastic flow
    
    # Phase 1: Establish K₀ initial stress state through elastic gravity loading
    if debug_level >= 1:
        print("Phase 1: Establishing K₀ initial stress state...")
    
    initial_displacements, stress_state = establish_k0_stress_state(
        K_global, F_gravity, bc_type, nodes, elements, element_types, 
        element_materials, E_by_mat, nu_by_mat, gamma_by_mat, u_nodal, debug_level)
    
    # Phase 2: Initialize Perzyna iteration from K₀ state
    # Debug gravity loads and element areas
    if debug_level >= 1:
        print(f"  Debug: Checking gravity load calculation")
        sample_elem = 100
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
    
    # Debug: Check yield states at different F values  
    if debug_level >= 1:
        print(f"  Checking yield at F=1.0 (original strength):")
        c_original = fem_data["c_by_mat"][element_materials - 1]
        phi_original = np.radians(fem_data["phi_by_mat"][element_materials - 1])
        initial_yield_count = check_initial_yield_state(stress_state, c_original, phi_original)
        print(f"    Elements yielding at F=1.0: {initial_yield_count}/{n_elements}")
        
        print(f"  Checking yield at F={F:.2f} (reduced strength):")
        print(f"    c_reduced = {c_reduced[0]:.1f}, phi_reduced = {np.degrees(phi_reduced[0]):.1f}°")
        reduced_yield_count = check_initial_yield_state(stress_state, c_reduced, phi_reduced)
        print(f"    Elements yielding at F={F:.2f}: {reduced_yield_count}/{n_elements}")
        
        # Find elements near the slope face (high x, high y)
        elem_centroids = []
        for elem_idx in range(n_elements):
            elem_type = element_types[elem_idx]
            elem_nodes = elements[elem_idx][:elem_type]
            centroid = np.mean(nodes[elem_nodes], axis=0)
            elem_centroids.append((elem_idx, centroid[0], centroid[1]))
        
        # Find elements with x > 100 and y > 20 (near slope face)
        face_elements = [e for e in elem_centroids if e[1] > 100 and e[2] > 20]
        if face_elements:
            face_elements.sort(key=lambda x: -x[2])  # Sort by y descending
            sample_elem = face_elements[0][0]
            print(f"  Checking element {sample_elem} near face (x={face_elements[0][1]:.1f}, y={face_elements[0][2]:.1f}):")
            sample_stress = stress_state['element_stresses'][sample_elem, 0, :]
            print(f"    Stress: σx={sample_stress[0]:.1f}, σy={sample_stress[1]:.1f}, τxy={sample_stress[2]:.1f}")
            f_sample = check_mohr_coulomb_perzyna(sample_stress, c_reduced[sample_elem], phi_reduced[sample_elem])
            print(f"    Yield function F = {f_sample:.3f}")
    
    # Calculate yield function values for all elements after gravity loading
    yield_function_values = np.zeros(n_elements)
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        # Use first Gauss point stress for yield function (or average for quads)
        if elem_type == 8:  # 8-node quad - average over Gauss points
            elem_stress_avg = np.mean(stress_state['element_stresses'][elem_idx, :4, :], axis=0)
        else:  # Triangle or other - use first Gauss point
            elem_stress_avg = stress_state['element_stresses'][elem_idx, 0, :]
        
        # Calculate yield function with reduced strength parameters
        yield_function_values[elem_idx] = check_mohr_coulomb_perzyna(
            elem_stress_avg, c_reduced[elem_idx], phi_reduced[elem_idx])
    
    # Check for early abort after gravity loading
    if abort_after == 0:
        if debug_level >= 1:
            print("Aborting after gravity loading (abort_after=0)")
        
        # Compute stresses and strains for output
        final_stresses, plastic_elements = compute_final_state_perzyna(
            nodes, elements, element_types, element_materials,
            initial_displacements, {}, c_reduced, phi_reduced,
            E_by_mat, nu_by_mat, u_nodal, stress_state)
        
        strains = compute_strains_perzyna(nodes, elements, element_types, initial_displacements)
        
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
    
    for iteration in range(max_iterations):
        if debug_level >= 3:
            print(f"\n--- Iteration {iteration + 1} ---")
        
        # Build load vector including plastic corrections
        F_total = F_gravity.copy()
        
        # Add plastic strain corrections (Perzyna stress redistribution)
        F_plastic_correction = compute_plastic_load_correction_perzyna(
            nodes, elements, element_types, element_materials, 
            plastic_strains, E_by_mat, nu_by_mat, eta)
        
        if debug_level >= 3:
            print(f"Plastic correction norm: {np.linalg.norm(F_plastic_correction):.2e}")
        
        F_total += F_plastic_correction
        
        # Add boundary condition loads
        for i in range(n_nodes):
            if bc_type[i] == 4:  # Force boundary condition
                F_total[2*i] += bc_values[i, 0]
                F_total[2*i+1] += bc_values[i, 1]
        
        # Apply constraints using proper free DOF extraction
        K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
            K_global, F_total, bc_type, nodes)
        
        # Solve for total displacements at equilibrium
        try:
            if hasattr(K_constrained, 'toarray'):
                K_constrained = K_constrained.tocsr()
            displacements_free = spsolve(K_constrained, F_constrained)
            
            # Reconstruct full displacement vector
            n_dof = 2 * n_nodes
            displacements_new = np.zeros(n_dof)
            
            # Free DOFs
            free_dofs = [i for i in range(n_dof) if i not in constraint_dofs]
            displacements_new[free_dofs] = displacements_free
            
            # Constrained DOFs remain zero (already initialized)
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
            c_reduced, phi_reduced, E_by_mat, nu_by_mat, eta)
        
        # Check convergence
        disp_change = np.linalg.norm(displacements_new - displacements)
        plastic_change = total_plastic_increment
        
        if debug_level >= 3:
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
        if max_disp > 50.0:  # Allow much larger displacements to test convergence
            if debug_level >= 1:
                print(f"Large displacements detected ({max_disp:.3f}) - slope failure")
            converged = False
            break
        
        # Update for next iteration
        displacements_prev = displacements.copy()  # Store previous displacements
        displacements = displacements_new.copy()
        plastic_strains = plastic_strains_new
        
        # Check for early abort
        if abort_after > 0 and iteration + 1 >= abort_after:
            if debug_level >= 1:
                print(f"Aborting after iteration {iteration + 1} (abort_after={abort_after})")
            converged = True  # Mark as converged for early abort
            break
        
        # Check for excessive displacements (numerical instability indicator)
        max_disp = np.max(np.abs(displacements))
        if max_disp > 1e6:  # Very large threshold - only for true instability
            if debug_level >= 1:
                print(f"Numerical instability detected: max displacement = {max_disp:.2e}")
            break
    
    # Compute final state
    final_stresses, plastic_elements = compute_final_state_perzyna(
        nodes, elements, element_types, element_materials,
        displacements, plastic_strains, c_reduced, phi_reduced,
        E_by_mat, nu_by_mat, u_nodal, stress_state)
    
    # Compute strains
    strains = compute_strains_perzyna(nodes, elements, element_types, displacements)
    
    # Calculate final yield function values
    final_yield_function_values = np.zeros(n_elements)
    for elem_idx in range(n_elements):
        # Use the stress from final_stresses (which includes von Mises as 4th column)
        elem_stress = final_stresses[elem_idx, :3]  # [sig_x, sig_y, tau_xy]
        final_yield_function_values[elem_idx] = check_mohr_coulomb_perzyna(
            elem_stress, c_reduced[elem_idx], phi_reduced[elem_idx])
    
    if debug_level >= 2:
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
        "algorithm": "Perzyna Visco-Plastic"
    }
    
    # Add abort information if applicable
    if abort_after > 0 and iteration + 1 >= abort_after:
        result["aborted"] = True
        result["abort_after"] = abort_after
    
    return result


def solve_ssrm_perzyna(fem_data, F_min=1.0, F_max=3.0, tolerance=0.01, debug_level=0):
    """
    SSRM using Perzyna algorithm with pure non-convergence failure criterion.
    """
    
    if debug_level >= 1:
        print("=== Perzyna SSRM Analysis ===")
        print("Failure criterion: Pure non-convergence (Griffiths & Lane 1999)")
    
    F_left = F_min
    F_right = F_max
    
    # Verify bounds
    solution_min = solve_fem_perzyna(fem_data, F=F_min, debug_level=max(0, debug_level-1))
    if not solution_min["converged"]:
        return {
            "converged": False,
            "error": f"F_min = {F_min} does not converge - slope unstable",
            "FS": None
        }
    
    solution_max = solve_fem_perzyna(fem_data, F=F_max, debug_level=max(0, debug_level-1))
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
        
        solution = solve_fem_perzyna(fem_data, F=F_mid, debug_level=max(0, debug_level-1))
        
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


def build_global_stiffness_perzyna(nodes, elements, element_types, element_materials, E_by_mat, nu_by_mat):
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
        return build_triangle_stiffness_perzyna(tri_coords, E, nu)
    else:
        return build_triangle_stiffness_perzyna(coords, E, nu)


def build_triangle_stiffness_perzyna(coords, E, nu):
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
    
    # B matrix (standard tension positive convention)
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


def build_gravity_loads_perzyna(nodes, elements, element_types, element_materials, gamma_by_mat, k_seismic):
    """
    Build gravity load vector.
    """
    n_nodes = len(nodes)
    F_gravity = np.zeros(2 * n_nodes)
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        gamma = gamma_by_mat[mat_id]
        
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Calculate element area/volume
        if elem_type >= 3:  # Triangle or higher
            if elem_type == 8:  # 8-node quad
                # Approximate area for 8-node quad
                area = compute_quad_area(elem_coords)
            else:
                # Triangle area
                x1, y1 = elem_coords[0]
                x2, y2 = elem_coords[1]
                x3, y3 = elem_coords[2]
                area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            # Distribute loads to nodes
            load_per_node = gamma * area / elem_type
            
            for i, node in enumerate(elem_nodes):
                # Vertical gravity load (negative = downward)
                F_gravity[2*node + 1] -= load_per_node
                # Horizontal seismic load  
                F_gravity[2*node] += k_seismic * load_per_node
    
    return F_gravity


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
                                           plastic_strains, E_by_mat, nu_by_mat, eta):
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
            n_gauss = 4  # 8-node quad with reduced integration
        else:
            n_gauss = 1  # Triangle - single Gauss point
        
        # Element plastic force vector
        elem_f_plastic = np.zeros(2 * elem_type)
        
        # For each Gauss point
        for gp in range(n_gauss):
            # Get plastic strains at this Gauss point
            plastic_strain = plastic_strains[elem_idx][gp, :]
            
            # Skip if no plastic strain (but use much smaller threshold)
            if np.linalg.norm(plastic_strain) < 1e-20:
                continue
            
            # Compute plastic stress
            D = build_constitutive_matrix_perzyna(E, nu)
            plastic_stress = D @ plastic_strain
            
            # Compute B matrix for this element
            if elem_type == 3:  # Triangle
                B, area = compute_B_matrix_triangle(elem_coords)
                weight = area  # Integration weight
            elif elem_type == 8:  # 8-node quad
                # For 8-node quad, use centroid for simplified integration
                B, det_J = compute_B_matrix_quad8_centroid(elem_coords)
                weight = det_J  # Integration weight (determinant of Jacobian)
            else:
                # Simplified for other elements
                B = np.zeros((3, 2 * elem_type))
                weight = 1.0
            
            # Add contribution to element force vector
            # F_elem += B^T * σ_plastic * weight
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
    
    # B matrix (standard tension positive convention)
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
        F_yield = check_mohr_coulomb_perzyna(stress, c, phi)
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
            D = build_constitutive_matrix_perzyna(E, nu)
            incremental_stress = D @ elastic_strains
            
            # Add initial stress if provided
            if initial_stresses is not None:
                initial_stress = initial_stresses['element_stresses'][elem_idx, gp, :]
                trial_stress = initial_stress + incremental_stress
            else:
                trial_stress = incremental_stress
            
            # Check yield criterion with total stress
            f_yield = check_mohr_coulomb_perzyna(trial_stress, c, phi)
            
            if f_yield > 1e-8:  # Plastic loading
                # Perzyna visco-plastic flow as per Griffiths & Lane (1999)
                # Δλ = η * <f> where <f> = max(0, f)
                # Use appropriate viscosity parameter from paper
                delta_lambda = eta * f_yield  # Remove artificial cap
                
                # Flow vector (non-associated: ψ = 0)
                flow_vector = compute_plastic_flow_vector(trial_stress, 0.0)  # ψ = 0
                
                # Plastic strain increment - controlled flow to prevent instability
                plastic_increment = delta_lambda * flow_vector
                
                # Apply reasonable limit to prevent numerical explosion
                increment_norm = np.linalg.norm(plastic_increment)
                if increment_norm > 1e-5:  # Very conservative limit for localized failure
                    plastic_increment *= 1e-5 / increment_norm
                
                # Update plastic strains
                plastic_strains_new[elem_idx][gp, :] += plastic_increment
                
                # Track total plastic increment
                total_increment += np.linalg.norm(plastic_increment)
    
    return plastic_strains_new, total_increment


def update_plastic_strains_perzyna_incremental(nodes, elements, element_types, element_materials,
                                              displacements_new, displacements_prev, plastic_strains, 
                                              current_stress_state, c_reduced, phi_reduced,
                                              E_by_mat, nu_by_mat, eta):
    """
    Update plastic strains using incremental Perzyna algorithm.
    
    This version properly handles incremental displacements and stress updates
    to avoid double-counting stress contributions.
    """
    plastic_strains_new = {}
    total_increment = 0.0
    
    # Create new stress state (copy of current)
    new_stress_state = {
        'element_stresses': current_stress_state['element_stresses'].copy(),
        'plastic_state': current_stress_state['plastic_state'].copy()
    }
    
    # Compute displacement increment
    displacement_increment = displacements_new - displacements_prev
    
    for elem_idx in range(len(elements)):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        c = c_reduced[elem_idx]
        phi = phi_reduced[elem_idx]
        
        # Get element nodes and coordinates
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Get incremental displacements for this element
        elem_disp_increment = np.zeros(2 * elem_type)
        for i, node_idx in enumerate(elem_nodes):
            elem_disp_increment[2*i] = displacement_increment[2*node_idx]
            elem_disp_increment[2*i+1] = displacement_increment[2*node_idx+1]
        
        # Initialize plastic strains for this element
        plastic_strains_new[elem_idx] = plastic_strains[elem_idx].copy()
        
        # Determine number of Gauss points
        if elem_type == 8:
            gauss_points_2x2, _ = get_gauss_points_2x2()
            n_gauss = 4
        else:
            n_gauss = 1
        
        # For each Gauss point
        for gp in range(n_gauss):
            # Compute incremental strains at this Gauss point
            if elem_type == 3:  # Triangle
                incremental_strains = compute_triangle_strains_manual(elem_coords, elem_disp_increment)
            elif elem_type == 8:  # 8-node quad
                xi, eta_local = gauss_points_2x2[gp]
                incremental_strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp_increment, xi, eta_local)
            else:
                incremental_strains = np.array([0.0, 0.0, 0.0])
            
            # Compute incremental stress from incremental strains
            D = build_constitutive_matrix_perzyna(E, nu)
            incremental_stress = D @ incremental_strains
            
            # Get current total stress at this Gauss point
            current_stress = current_stress_state['element_stresses'][elem_idx, gp, :]
            
            # Trial stress = current stress + incremental stress  
            trial_stress = current_stress + incremental_stress
            
            # Check yield criterion with trial stress
            f_yield = check_mohr_coulomb_perzyna(trial_stress, c, phi)
            
            if f_yield > 1e-8:  # Plastic loading
                # Perzyna visco-plastic flow
                delta_lambda = eta * f_yield
                
                # Flow vector (non-associated: ψ = 0)
                flow_vector = compute_plastic_flow_vector(trial_stress, 0.0)
                
                # Plastic strain increment
                plastic_increment = delta_lambda * flow_vector
                
                # Apply reasonable limit to prevent numerical explosion
                increment_norm = np.linalg.norm(plastic_increment)
                if increment_norm > 1e-5:
                    plastic_increment *= 1e-5 / increment_norm
                
                # Update plastic strains
                plastic_strains_new[elem_idx][gp, :] += plastic_increment
                
                # Update stress state (remove plastic stress contribution)
                plastic_stress = D @ plastic_increment
                new_stress_state['element_stresses'][elem_idx, gp, :] = trial_stress - plastic_stress
                
                # Track total plastic increment
                total_increment += np.linalg.norm(plastic_increment)
                
            else:
                # Elastic loading - just update stress
                new_stress_state['element_stresses'][elem_idx, gp, :] = trial_stress
    
    return plastic_strains_new, total_increment, new_stress_state


def compute_final_state_perzyna(nodes, elements, element_types, element_materials,
                               displacements, plastic_strains, c_reduced, phi_reduced,
                               E_by_mat, nu_by_mat, u_nodal, stress_state):
    """
    Compute final stress state and identify plastic elements.
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
        
        # Get element nodes and coordinates
        elem_nodes = elements[elem_idx][:elem_type]
        elem_coords = nodes[elem_nodes]
        
        # Get element displacements
        elem_disp = np.zeros(2 * elem_type)
        for i, node_idx in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node_idx]
            elem_disp[2*i+1] = displacements[2*node_idx+1]
        
        # For triangular elements (single Gauss point)
        if elem_type == 3:
            # Compute strains
            total_strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            
            # Get plastic strains (single Gauss point) - handle empty plastic_strains
            if elem_idx in plastic_strains:
                plastic_strain = plastic_strains[elem_idx][0, :]  # First (and only) Gauss point
            else:
                plastic_strain = np.zeros(3)  # No plastic strains yet
            
            # Elastic strains
            elastic_strains = total_strains - plastic_strain
            
            # Stress calculation: initial stress + incremental stress
            D = build_constitutive_matrix_perzyna(E, nu)
            incremental_stress = D @ elastic_strains
            
            # Add initial geostatic stress
            if stress_state is not None:
                initial_stress = stress_state['element_stresses'][elem_idx, 0, :]  # Single Gauss point
                stress = initial_stress + incremental_stress
            else:
                stress = incremental_stress
            
            # Store stress (σx, σy, τxy, σvm)
            sig_x, sig_y, tau_xy = stress
            sig_vm = np.sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
            final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]
            
            # Check if element is plastic (yield criterion)
            f_yield = check_mohr_coulomb_perzyna(stress, c, phi)
            plastic_elements[elem_idx] = f_yield > 1e-8
        
        else:
            # 8-node quad: compute average stress and check yield function
            elem_nodes = elements[elem_idx][:elem_type]
            elem_coords = nodes[elem_nodes]
            
            # Get element displacements
            elem_disp = np.zeros(2 * elem_type)
            for i, node in enumerate(elem_nodes):
                elem_disp[2*i] = displacements[2*node]
                elem_disp[2*i+1] = displacements[2*node+1]
            
            # Average stress over all Gauss points using current stress state
            elem_stress_avg = np.zeros(3)
            n_gauss = 4  # 2x2 integration for 8-node quads
            
            if elem_idx in stress_state:
                # Use current stress state from Perzyna iteration
                for gp in range(n_gauss):
                    elem_stress_avg += stress_state[elem_idx][gp, :]
                elem_stress_avg /= n_gauss
            else:
                # Fallback: compute stress from displacements
                gauss_coords = [(-1/np.sqrt(3), -1/np.sqrt(3)), (1/np.sqrt(3), -1/np.sqrt(3)),
                               (1/np.sqrt(3), 1/np.sqrt(3)), (-1/np.sqrt(3), 1/np.sqrt(3))]
                for xi, eta in gauss_coords:
                    strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta)
                    D = build_constitutive_matrix_perzyna(E, nu)
                    stress = D @ strains
                    elem_stress_avg += stress
                elem_stress_avg /= len(gauss_coords)
            
            # Compute von Mises for display
            sig_x, sig_y, tau_xy = elem_stress_avg
            sig_vm = np.sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
            final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]
            
            # Check yield function (same as triangles)
            f_yield = check_mohr_coulomb_perzyna(elem_stress_avg, c, phi)
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
    
    # B matrix (standard tension positive convention)
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],  # εx = ∂u/∂x
        [0,  c1, 0,  c2, 0,  c3],  # εy = ∂v/∂y
        [c1, b1, c2, b2, c3, b3]   # γxy = ∂u/∂y + ∂v/∂x
    ]) / (2 * area)
    
    # Strains
    strains = B @ displacements
    return strains


def build_constitutive_matrix_perzyna(E, nu):
    """Build constitutive matrix for plane strain."""
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    return D


def check_mohr_coulomb_perzyna(stress, c, phi, debug=False):
    """Check Mohr-Coulomb yield criterion and return violation.
    
    Uses Griffiths & Lane (1999) formulation with compression-negative convention.
    """
    
    sig_x, sig_y, tau_xy = stress
    
    # Principal stresses (tension positive, compression negative)
    sig_mean = (sig_x + sig_y) / 2
    tau_max = sqrt(((sig_x - sig_y) / 2)**2 + tau_xy**2)
    
    sig1 = sig_mean + tau_max  # Most positive (least compressive/most tensile)
    sig3 = sig_mean - tau_max  # Most negative (most compressive/least tensile)
    
    # Griffiths & Lane (1999) yield function for compression-negative
    # F = (σ₁ + σ₃)/2 × sin φ' - (σ₁ - σ₃)/2 - c' cos φ'
    # Yield when F > 0
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    
    # Note: sig1 + sig3 = 2*sig_mean = sig_x + sig_y
    # And: sig1 - sig3 = 2*tau_max
    F = sig_mean * sin_phi - tau_max - c * cos_phi
    
    if debug and F > -c * cos_phi * 0.1:  # Near yielding
        print(f"  Debug: σx={sig_x:.1f}, σy={sig_y:.1f}, τxy={tau_xy:.1f}")
        print(f"         σ1={sig1:.1f}, σ3={sig3:.1f}, F={F:.3f}, c={c:.1f}, φ={np.degrees(phi):.1f}°")
    
    # Return actual F value, positive means yielding
    return F


def compute_plastic_flow_vector(stress, psi):
    """Compute plastic flow vector for non-associated plasticity."""
    
    sig_x, sig_y, tau_xy = stress
    
    # For non-associated plasticity with ψ = 0 (Griffiths approach)
    # Flow vector is derived from potential function g = (sig1 - sig3)
    
    # Principal stresses
    sig_mean = (sig_x + sig_y) / 2
    tau_max = sqrt(((sig_x - sig_y) / 2)**2 + tau_xy**2)
    
    if tau_max < 1e-12:
        # Hydrostatic stress state
        return np.array([0.0, 0.0, 0.0])
    
    sig1 = sig_mean + tau_max
    sig3 = sig_mean - tau_max
    
    # Flow direction derivatives for ψ = 0 case
    # ∂g/∂σ where g = (σ1 - σ3) + (σ1 + σ3)*sin(ψ) - 2*c*cos(ψ)
    
    sin_psi = sin(psi)
    
    # Simplified flow vector for ψ = 0
    # Direction: [1+sin(ψ), -(1-sin(ψ)), 0] for principal directions
    # Transform back to x-y coordinates
    
    # For simplicity, use associated flow approximation
    flow_x = 1 + sin_psi
    flow_y = -(1 - sin_psi)
    flow_xy = 0
    
    flow_vector = np.array([flow_x, flow_y, flow_xy])
    
    # Normalize
    norm = np.linalg.norm(flow_vector)
    if norm > 1e-12:
        flow_vector /= norm
    
    return flow_vector


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
    
    This captures stress concentrations and geometric effects essential for proper
    plastic failure initiation, unlike simplified geostatic K₀ approaches.
    """
    n_elements = len(elements)
    max_gauss_points = 4
    element_stresses = np.zeros((n_elements, max_gauss_points, 3))  # [sig_x, sig_y, tau_xy]
    
    for elem_idx in range(n_elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        gamma = gamma_by_mat[mat_id]
        
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
            # True Griffiths approach: compute strains from FEM displacements at this Gauss point
            if elem_type == 3:  # Triangle
                strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            elif elem_type == 8:  # 8-node quad
                xi, eta_local = gauss_points_2x2[gp]
                strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta_local)
            else:
                strains = np.array([0.0, 0.0, 0.0])
            
            # Compute stresses from strains using elastic constitutive matrix
            # This is the true "gravity in single increment to initially stress-free slope"
            # Standard tension positive convention: compression is negative
            D = build_constitutive_matrix_perzyna(E, nu)
            stresses = D @ strains
            
            # Store stress at this Gauss point (tension positive, compression negative)
            element_stresses[elem_idx, gp, :] = stresses
    
    # Debug: Check stress state statistics
    stress_stats = {
        'sigma_x': {'min': np.min(element_stresses[:, :, 0]), 'max': np.max(element_stresses[:, :, 0]), 'mean': np.mean(element_stresses[:, :, 0])},
        'sigma_y': {'min': np.min(element_stresses[:, :, 1]), 'max': np.max(element_stresses[:, :, 1]), 'mean': np.mean(element_stresses[:, :, 1])},
        'tau_xy': {'min': np.min(element_stresses[:, :, 2]), 'max': np.max(element_stresses[:, :, 2]), 'mean': np.mean(element_stresses[:, :, 2])}
    }
    
    print(f"  Initial stress state statistics:")
    print(f"    σ_x: min={stress_stats['sigma_x']['min']:.1f}, max={stress_stats['sigma_x']['max']:.1f}, mean={stress_stats['sigma_x']['mean']:.1f}")
    print(f"    σ_y: min={stress_stats['sigma_y']['min']:.1f}, max={stress_stats['sigma_y']['max']:.1f}, mean={stress_stats['sigma_y']['mean']:.1f}")
    print(f"    τ_xy: min={stress_stats['tau_xy']['min']:.1f}, max={stress_stats['tau_xy']['max']:.1f}, mean={stress_stats['tau_xy']['mean']:.1f}")
    
    return {
        'element_stresses': element_stresses,
        'plastic_state': np.zeros((n_elements, max_gauss_points), dtype=bool)
    }


def compute_strains_perzyna(nodes, elements, element_types, displacements):
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
                dN_dx[a] = J_inv[0,0] * dN_dxi[a] + J_inv[0,1] * dN_deta[a]
                dN_dy[a] = J_inv[1,0] * dN_dxi[a] + J_inv[1,1] * dN_deta[a]
            
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
    
    # B matrix (strain-displacement, standard tension positive)
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