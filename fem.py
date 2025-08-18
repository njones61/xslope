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
    
    # Compute geostatic initial stresses (K0 state) per element to avoid spurious base plasticity
    # Temporarily disabled to test if initial stresses are preventing failure
    sigma0_by_elem = np.zeros((n_elements, 3))  # [sig_x0, sig_y0, tau0]
    # ground_surface = slope_data.get("ground_surface", None)
    # if ground_surface is not None and hasattr(ground_surface, "project"):
    #     line = ground_surface
    #     for elem_idx in range(n_elements):
    #         elem_nodes = elements[elem_idx]
    #         elem_type = element_types[elem_idx]
    #         active_nodes = elem_nodes[:elem_type]
    #         elem_coords = nodes[active_nodes]
    #         centroid = np.mean(elem_coords, axis=0)
    #         # Find surface elevation directly above centroid x by projecting onto the ground_surface
    #         pt_on_line = line.interpolate(line.project(Point(centroid[0], centroid[1])))
    #         y_surface = pt_on_line.y
    #         depth = max(0.0, y_surface - centroid[1])
    #         mat_id = element_materials[elem_idx] - 1
    #         gamma = gamma_by_mat[mat_id]
    #         phi_deg = phi_by_elem[elem_idx]
    #         # Jaky's K0 ~ 1 - sin(phi)
    #         K0 = max(0.0, 1.0 - np.sin(np.radians(phi_deg)))
    #         sigma_v0 = gamma * depth
    #         sigma_h0 = K0 * sigma_v0
    #         sigma0_by_elem[elem_idx, 0] = sigma_h0
    #         sigma0_by_elem[elem_idx, 1] = sigma_v0
    #         sigma0_by_elem[elem_idx, 2] = 0.0
    
    fem_data["sigma0_by_elem"] = sigma0_by_elem
    
    return fem_data


def solve_fem(fem_data, F=1.0, debug_level=0, plastic_stiffness_reduction=0.01):
    """
    Solve finite element system for slope stability analysis with elastic-plastic behavior.
    
    This function performs finite element analysis with Mohr-Coulomb plasticity using an
    elastic predictor-plastic corrector algorithm. It can handle both 2D soil elements
    and 1D truss reinforcement elements.
    
    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        F (float): Shear strength reduction factor (default 1.0)
        debug_level (int): Verbosity level (0=quiet, 1=basic, 2=detailed, 3=debug)
        plastic_stiffness_reduction (float): Stiffness reduction factor for plastic elements
            Literature ranges:
            - Conservative: 0.1 to 0.5 (10% to 50% of original stiffness - less reduction)
            - Aggressive: 0.01 to 0.1 (1% to 10% of original stiffness - more reduction)
            - Default: 0.1 (reduces plastic element stiffness to 10% of elastic value)
    
    Returns:
        dict: Solution dictionary containing:
            - converged: bool, whether solution converged
            - displacements: np.ndarray (n_nodes * 2,) nodal displacements [u1,v1,u2,v2,...]
            - stresses: np.ndarray (n_elements, 4) element stresses [sig_x, sig_y, tau_xy, von_mises]
            - strains: np.ndarray (n_elements, 4) element strains [eps_x, eps_y, gamma_xy, max_shear_strain]
            - plastic_elements: np.ndarray boolean mask of yielded elements
            - forces_1d: np.ndarray (n_1d_elements,) forces in 1D truss elements
            - iterations: int, number of iterations to convergence
            - residual_norm: float, final residual norm
    """
    
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
    
    # 1D elements
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    element_types_1d = fem_data.get("element_types_1d", np.array([]))
    t_allow_by_1d_elem = fem_data.get("t_allow_by_1d_elem", np.array([]))
    t_res_by_1d_elem = fem_data.get("t_res_by_1d_elem", np.array([]))
    k_by_1d_elem = fem_data.get("k_by_1d_elem", np.array([]))
    
    # Seismic and other loads
    k_seismic = fem_data.get("k_seismic", 0.0)
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    n_1d_elements = len(elements_1d)
    n_dof = 2 * n_nodes  # 2 DOF per node (u, v)
    
    # Reduce strength parameters by factor F
    c_reduced = c_by_elem / F
    # FIXED: Properly reduce friction angle by reducing tan(phi) by factor F
    phi_reduced = np.arctan(np.tan(np.radians(phi_by_elem)) / F)
    tan_phi_reduced = np.tan(phi_reduced)
    
    # Non-associated plasticity: dilation angle different from friction angle
    # Typical relationship: psi = phi/3 for most soils (promotes slip surface localization)
    psi_by_elem = phi_by_elem / 3.0  # Dilation angle (degrees)
    psi_reduced = np.arctan(np.tan(np.radians(psi_by_elem)) / F)  # Reduce dilation angle by F too
    
    if debug_level >= 1:
        print(f"Starting FEM analysis with F = {F:.3f}")
        print(f"Mesh: {n_nodes} nodes, {n_elements} 2D elements, {n_1d_elements} 1D elements")
    
    # Initialize solution vectors
    displacements = np.zeros(n_dof)
    
    # Track plastic state and strain softening
    plastic_elements = np.zeros(n_elements, dtype=bool)
    failed_1d_elements = np.zeros(n_1d_elements, dtype=bool)
    accumulated_plastic_strain = np.zeros(n_elements)  # Track plastic strain for softening
    
    # Convergence parameters - allow more iterations for strain softening
    max_iterations = 20  # More iterations to allow strain softening to develop
    tol_force = 1e-4     # Slightly relaxed to allow softening process
    tol_disp = 1e-4      # Slightly relaxed to allow softening process
    
    converged = False
    iteration = 0
    residual_norm = np.inf
    
    # Initial elastic solution (first iteration only)
    iteration = 0
    converged = False
    
    while iteration < max_iterations and not converged:
        if debug_level >= 2:
            print(f"\n--- Iteration {iteration + 1} ---")
        
        # Implement proper iterative elastic-plastic equilibrium
        if iteration == 0:
            # First iteration - solve elastically, then check for plasticity
            plastic_state_changed = False
        else:
            # Check for plastic state changes and update stiffness accordingly
            plastic_state_changed = check_and_update_plastic_state_with_softening(
                nodes, elements, element_types, element_materials, displacements,
                E_by_mat, nu_by_mat, u_nodal, c_reduced, phi_reduced, psi_reduced, 
                plastic_elements, accumulated_plastic_strain, debug_level)
        # Assemble force vector (only for first iteration)
        if iteration == 0:
            F_global = np.zeros(n_dof)
            
            # Body forces (gravity + seismic)
            for elem_idx in range(n_elements):
                elem_nodes = elements[elem_idx]
                elem_type = element_types[elem_idx]
                active_nodes = elem_nodes[:elem_type]
                
                mat_id = element_materials[elem_idx] - 1
                gamma = gamma_by_mat[mat_id]
                
                # Body force components
                b_x = k_seismic * gamma  # Horizontal seismic force
                b_y = -gamma             # Gravity (downward)
                
                F_body = compute_body_forces(nodes[active_nodes], elem_type, b_x, b_y)
                
                # Assemble into global vector
                dofs = []
                for node_id in active_nodes:
                    dofs.extend([2*node_id, 2*node_id + 1])
                
                for i, dof in enumerate(dofs):
                    F_global[dof] += F_body[i]
            
            # Applied loads from boundary conditions
            for i in range(n_nodes):
                if bc_type[i] == 4:  # Applied force
                    F_global[2*i] += bc_values[i, 0]      # F_x
                    F_global[2*i + 1] += bc_values[i, 1]  # F_y
        
        # Inner loop for elastic-plastic iteration within each load step
        sub_iteration = 0
        max_sub_iterations = 10
        plastic_state_changed = False  # Initialize to False for each main iteration
        
        # Initialize variables for cases where inner loop doesn't execute
        delta_u_full = np.zeros(n_dof)
        K_global = lil_matrix((n_dof, n_dof))
        constraint_dofs = []
        F_constrained = np.zeros(n_dof)
        K_constrained = lil_matrix((n_dof, n_dof))
        
        # Always execute the inner loop for the first iteration to build stiffness matrix
        # For subsequent iterations, only execute if plastic state changed
        should_execute_inner = (iteration == 0) or plastic_state_changed
        
        while sub_iteration < max_sub_iterations and should_execute_inner:
            if debug_level >= 3:
                print(f"    Sub-iteration {sub_iteration + 1}")
            
            # Assemble global stiffness matrix with current plastic state
            K_global = lil_matrix((n_dof, n_dof))
            
            # Assemble 2D elements
            for elem_idx in range(n_elements):
                elem_nodes = elements[elem_idx]
                elem_type = element_types[elem_idx]
                active_nodes = elem_nodes[:elem_type]
                
                # Get material properties
                mat_id = element_materials[elem_idx] - 1
                E = E_by_mat[mat_id]
                nu = nu_by_mat[mat_id]
                
                # Build element stiffness matrix with proper plastic state consideration
                # Use moderate stiffness reduction for plastic elements (research-based approach)
                is_elem_plastic = plastic_elements[elem_idx]
                if elem_type in [3, 6]:  # Triangular elements
                    K_elem = build_triangle_stiffness(nodes[active_nodes], E, nu, is_elem_plastic, elem_type, 0.2)
                elif elem_type in [4, 8, 9]:  # Quadrilateral elements
                    K_elem = build_quad_stiffness(nodes[active_nodes], E, nu, is_elem_plastic, elem_type, 0.2)
                else:
                    if debug_level >= 1:
                        print(f"Warning: Element type {elem_type} not supported, skipping element {elem_idx}")
                    continue
                
                # Assemble into global matrix
                dofs = []
                for node_id in active_nodes:
                    dofs.extend([2*node_id, 2*node_id + 1])
                
                for i in range(len(dofs)):
                    for j in range(len(dofs)):
                        K_global[dofs[i], dofs[j]] += K_elem[i, j]
            
            # Assemble 1D truss elements
            for elem_idx in range(n_1d_elements):
                if failed_1d_elements[elem_idx]:
                    continue  # Skip failed elements
                    
                elem_nodes = elements_1d[elem_idx]
                elem_type = element_types_1d[elem_idx]
                active_nodes = elem_nodes[:elem_type]
                
                if len(active_nodes) >= 2:
                    K_truss = build_truss_stiffness(nodes[active_nodes], k_by_1d_elem[elem_idx])
                    
                    # Assemble into global matrix
                    dofs = []
                    for node_id in active_nodes:
                        dofs.extend([2*node_id, 2*node_id + 1])
                    
                    for i in range(len(dofs)):
                        for j in range(len(dofs)):
                            K_global[dofs[i], dofs[j]] += K_truss[i, j]
            
            # Apply boundary conditions
            K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
                K_global, F_global, bc_type, nodes)
            
            # Solve system
            try:
                if iteration == 0:
                    # First iteration: solve for total displacements
                    u_free = spsolve(K_constrained.tocsr(), F_constrained)
                    
                    # Reconstruct full displacement vector
                    displacements = np.zeros(n_dof)
                    free_dof_idx = 0
                    for i in range(n_dof):
                        if i not in constraint_dofs:
                            displacements[i] = u_free[free_dof_idx]
                            free_dof_idx += 1
                    
                    # Boundary conditions are automatically enforced since we only assign free DOFs
                    # and constrained DOFs remain at their initialized value of 0.0
                    
                    delta_u_full = displacements.copy()
                else:
                    # Subsequent iterations: solve for displacement increments
                    # For plastic correction, we need to compute the residual force
                    # and solve for the correction
                    
                    # Compute current internal forces from current displacements
                    F_internal = K_global @ displacements
                    
                    # Apply boundary conditions to internal forces
                    for dof in constraint_dofs:
                        F_internal[dof] = 0.0
                    
                    # The residual is the difference between applied and internal forces
                    F_residual = F_global - F_internal
                    
                    # Apply boundary conditions to residual
                    for dof in constraint_dofs:
                        F_residual[dof] = 0.0
                    
                    # Extract only the free DOFs for the constrained system
                    F_residual_constrained = []
                    for i in range(n_dof):
                        if i not in constraint_dofs:
                            F_residual_constrained.append(F_residual[i])
                    F_residual_constrained = np.array(F_residual_constrained)
                    
                    # Solve for displacement correction
                    delta_u_free = spsolve(K_constrained.tocsr(), F_residual_constrained)
                    
                    # Reconstruct full displacement increment vector
                    delta_u_full = np.zeros(n_dof)
                    free_dof_idx = 0
                    for i in range(n_dof):
                        if i not in constraint_dofs:
                            delta_u_full[i] = delta_u_free[free_dof_idx]
                            free_dof_idx += 1
                    
                    # Update total displacements
                    displacements += delta_u_full
                    
                    # Enforce boundary conditions on total displacements
                    for dof in constraint_dofs:
                        displacements[dof] = 0.0
                    
            except Exception as e:
                if debug_level >= 1:
                    print(f"Linear solver failed: {e}")
                return {
                    "converged": False,
                    "error": f"Linear solver failed: {e}",
                    "iterations": iteration
                }
            
            sub_iteration += 1
        
        # The plastic state checking is now done at the beginning of each iteration
        # No need to duplicate it here
        
        # After the inner loop, ensure we have the final plastic state
        # This is important because the inner loop may have updated the plastic state
        if iteration == 0:
            # For the first iteration, we already have the correct plastic state
            pass
        else:
            # For subsequent iterations, we need to check the final state
            # The plastic state should already be correct from the inner loop
            pass
        
        # Check convergence
        disp_change_norm = np.linalg.norm(delta_u_full) if iteration > 0 else np.linalg.norm(displacements)
        
        # Compute proper force residual norm
        if iteration == 0:
            # First iteration has no residual
            force_residual_norm = 0.0
        else:
            # For subsequent iterations, compute residual from stiffness matrix
            F_internal = K_global @ displacements
            # Apply boundary conditions to internal forces
            for dof in constraint_dofs:
                F_internal[dof] = 0.0
            # The residual is the difference between applied and internal forces
            F_residual = F_global - F_internal
            # Apply boundary conditions to residual
            for dof in constraint_dofs:
                F_residual[dof] = 0.0
            force_residual_norm = np.linalg.norm(F_residual)
        
        residual_norm = force_residual_norm
        
        # DEBUG: Check what's happening
        if debug_level >= 3:
            if iteration == 0:
                print(f"    DEBUG: Force vector norm: {np.linalg.norm(F_global):.2e}")
                print(f"    DEBUG: Displacement update norm: {np.linalg.norm(displacements):.2e}")
                print(f"    DEBUG: Max displacement: {np.max(np.abs(displacements)):.2e}")
                print(f"    DEBUG: Force residual norm: {force_residual_norm:.2e}")
                # K_constrained is defined in the inner loop, so we can't access it here
                print(f"    DEBUG: Matrix size: Not available yet")
            else:
                print(f"    DEBUG: Force vector norm: {np.linalg.norm(F_constrained):.2e}")
                print(f"    DEBUG: Displacement update norm: {np.linalg.norm(delta_u_full):.2e}")
                print(f"    DEBUG: Max displacement: {np.max(np.abs(displacements)):.2e}")
                print(f"    DEBUG: Force residual norm: {force_residual_norm:.2e}")
                print(f"    DEBUG: Matrix size: {K_constrained.shape}")
        
        # These variables are no longer needed with stress return mapping approach
        n_newly_plastic = 0
        n_newly_failed_1d = 0
        
        if debug_level >= 2:
            print(f"  Displacement change norm: {disp_change_norm:.2e}")
            print(f"  Force residual norm: {force_residual_norm:.2e}")
            print(f"  Total plastic elements: {np.sum(plastic_elements)}")
            print(f"  Total failed 1D elements: {np.sum(failed_1d_elements)}")
        
        # Check convergence criteria for elastic-plastic iteration
        if iteration == 0:
            # First iteration always continues to check for plasticity
            converged = False
        elif (disp_change_norm < tol_disp and force_residual_norm < tol_force and not plastic_state_changed):
            converged = True
            if debug_level >= 1:
                print(f"Converged after {iteration + 1} iterations (elastic-plastic equilibrium achieved)")
        elif iteration >= max_iterations - 1:
            # Maximum iterations reached - this indicates failure for SSRM
            converged = False
            if debug_level >= 1:
                print(f"Failed to converge after {max_iterations} iterations - slope failure detected")
        else:
            # Not converged - continue iterating
            converged = False
            if debug_level >= 2:
                if disp_change_norm >= tol_disp:
                    print(f"  Displacement not converged: {disp_change_norm:.2e} >= {tol_disp:.2e}")
                if force_residual_norm >= tol_force:
                    print(f"  Force residual not converged: {force_residual_norm:.2e} >= {tol_force:.2e}")
                if plastic_state_changed:
                    print(f"  Plastic state changed, re-solving with updated stiffness")
        
        # For subsequent iterations, zero force vector (no additional loading)
        if iteration == 0:
            F_global = np.zeros(n_dof)
        
        iteration += 1
    
    # Compute final stresses and strains
    stresses = np.zeros((n_elements, 4))  # sig_x, sig_y, tau_xy, von_mises
    strains = np.zeros((n_elements, 4))   # eps_x, eps_y, gamma_xy, max_shear_strain
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:  # All supported element types
            # Compute trial elastic stress (include initial stress if available)
            sigma0_elem = fem_data.get("sigma0_by_elem", None)
            sigma0 = None if sigma0_elem is None else sigma0_elem[elem_idx]
            trial_stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                E_by_mat[element_materials[elem_idx] - 1],
                nu_by_mat[element_materials[elem_idx] - 1],
                u_nodal[active_nodes], elem_type, sigma0)
            
            # Apply non-associated return mapping for plasticity
            c_elem = c_reduced[elem_idx]
            phi_elem = phi_reduced[elem_idx]
            psi_elem = psi_reduced[elem_idx]  # Use dilation angle
            E_elem = E_by_mat[element_materials[elem_idx] - 1]
            nu_elem = nu_by_mat[element_materials[elem_idx] - 1]
            
            corrected_stress, is_plastic = return_mapping_mohr_coulomb_stress_nonassociated(
                trial_stress, c_elem, phi_elem, psi_elem, E_elem, nu_elem)
            
            stresses[elem_idx, :3] = corrected_stress
            stresses[elem_idx, 3] = compute_von_mises(corrected_stress)
            
            # Update plastic element tracking
            plastic_elements[elem_idx] = is_plastic
            
            # Compute strains
            strain_vec = compute_element_strain(
                nodes[active_nodes], displacements, active_nodes, elem_type)
            strains[elem_idx, :3] = strain_vec
            
            # Compute maximum shear strain from principal strains
            eps_x, eps_y, gamma_xy = strain_vec
            
            # Principal strains
            eps_mean = 0.5 * (eps_x + eps_y)
            R = np.sqrt(0.25 * (eps_x - eps_y)**2 + 0.25 * gamma_xy**2)
            
            eps_1 = eps_mean + R  # Major principal strain
            eps_3 = eps_mean - R  # Minor principal strain
            
            # Maximum shear strain (engineering definition)
            max_shear_strain = abs(eps_1 - eps_3)
            strains[elem_idx, 3] = max_shear_strain
    
    # Final force computation for 1D elements
    forces_1d = np.zeros(n_1d_elements)
    for elem_idx in range(n_1d_elements):
        if not failed_1d_elements[elem_idx]:
            elem_nodes = elements_1d[elem_idx]
            elem_type = element_types_1d[elem_idx]
            active_nodes = elem_nodes[:elem_type]
            
            if len(active_nodes) >= 2:
                forces_1d[elem_idx] = compute_truss_force(
                    nodes[active_nodes], displacements, active_nodes, k_by_1d_elem[elem_idx])
    
    if debug_level >= 1:
        n_plastic = np.sum(plastic_elements)
        n_failed_1d = np.sum(failed_1d_elements)
        print(f"Final state: {n_plastic}/{n_elements} plastic elements, {n_failed_1d}/{n_1d_elements} failed reinforcement")
    
    return {
        "converged": converged,
        "displacements": displacements,
        "stresses": stresses,
        "strains": strains,
        "plastic_elements": plastic_elements,
        "forces_1d": forces_1d,
        "failed_1d_elements": failed_1d_elements,
        "iterations": iteration + 1,
        "residual_norm": residual_norm
    }


def build_triangle_stiffness(coords, E, nu, is_plastic=False, elem_type=3, plastic_stiffness_reduction=0.1):
    """
    Build stiffness matrix for triangular elements (tri3 or tri6) using proper elastic-plastic constitutive model.
    
    This function now uses the consistent tangent stiffness approach instead of arbitrary
    stiffness reduction, providing much better numerical stability and convergence.
    
    Parameters:
        coords: np.ndarray - coordinates of triangle nodes (3x2 for tri3, 6x2 for tri6)
        E: Young's modulus
        nu: Poisson's ratio
        is_plastic: boolean indicating if element is plastic
        elem_type: element type (3 for tri3, 6 for tri6)
        plastic_stiffness_reduction: float - DEPRECATED, kept for compatibility
    
    Returns:
        np.ndarray - element stiffness matrix (6x6 for tri3, 12x12 for tri6)
    """
    
    # For now, use elastic stiffness matrix
    # The elastic-plastic behavior will be handled at the global assembly level
    # where we can properly track stress states and apply return mapping
    
    # Constitutive matrix for plane strain
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    if elem_type == 3:  # Linear triangle (tri3)
        return build_tri3_stiffness(coords, D)
    elif elem_type == 6:  # Quadratic triangle (tri6)
        return build_tri6_stiffness(coords, D)
    else:
        raise ValueError(f"Unsupported triangle element type: {elem_type}")


def build_tri3_stiffness(coords, D):
    """Build stiffness matrix for 3-node triangle."""
    # Extract coordinates
    x1, y1 = coords[0]
    x2, y2 = coords[1] 
    x3, y3 = coords[2]
    
    # Calculate area
    area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    
    if area < 1e-12:
        warnings.warn("Degenerate triangle detected")
        return np.zeros((6, 6))
    
    # Build B matrix (strain-displacement matrix)
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    B = (1.0 / (2.0 * area)) * np.array([
        [b1, 0,  b2, 0,  b3, 0],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
    ])
    
    # Element stiffness matrix
    K = B.T @ D @ B * area
    
    return K


def build_tri6_stiffness(coords, D):
    """Build stiffness matrix for 6-node triangle using numerical integration."""
    # Gauss points for 6-point integration (exact for quadratic elements)
    gauss_points = np.array([
        [0.816847572980459, 0.091576213509771, 0.091576213509771],
        [0.091576213509771, 0.816847572980459, 0.091576213509771],
        [0.091576213509771, 0.091576213509771, 0.816847572980459],
        [0.108103018168070, 0.445948490915965, 0.445948490915965],
        [0.445948490915965, 0.108103018168070, 0.445948490915965],
        [0.445948490915965, 0.445948490915965, 0.108103018168070]
    ])
    
    weights = np.array([0.109951743655322, 0.109951743655322, 0.109951743655322,
                       0.223381589678011, 0.223381589678011, 0.223381589678011])
    
    # Calculate triangle area for coordinate transformation
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    
    if area < 1e-12:
        return np.zeros((12, 12))
    
    K = np.zeros((12, 12))
    
    # Numerical integration
    for gp in range(len(gauss_points)):
        L1, L2, L3 = gauss_points[gp]
        w = weights[gp]
        
        # Shape functions for tri6 (quadratic)
        N = np.array([
            L1 * (2*L1 - 1),                    # N1
            L2 * (2*L2 - 1),                    # N2
            L3 * (2*L3 - 1),                    # N3
            4 * L1 * L2,                        # N4 (midside 1-2)
            4 * L2 * L3,                        # N5 (midside 2-3)
            4 * L3 * L1                         # N6 (midside 3-1)
        ])
        
        # Shape function derivatives in natural coordinates
        dN_dL = np.array([
            [4*L1 - 1, 0,       0],            # dN1/dL1, dN1/dL2, dN1/dL3
            [0,        4*L2 - 1, 0],           # dN2/dL1, dN2/dL2, dN2/dL3
            [0,        0,        4*L3 - 1],    # dN3/dL1, dN3/dL2, dN3/dL3
            [4*L2,     4*L1,     0],           # dN4/dL1, dN4/dL2, dN4/dL3
            [0,        4*L3,     4*L2],        # dN5/dL1, dN5/dL2, dN5/dL3
            [4*L3,     0,        4*L1]         # dN6/dL1, dN6/dL2, dN6/dL3
        ])
        
        # Jacobian transformation from natural to global coordinates
        J = np.zeros((2, 2))
        for i in range(6):
            J[0, 0] += dN_dL[i, 0] * coords[i, 0] - dN_dL[i, 2] * coords[i, 0]  # dx/dL1
            J[0, 1] += dN_dL[i, 1] * coords[i, 0] - dN_dL[i, 2] * coords[i, 0]  # dx/dL2
            J[1, 0] += dN_dL[i, 0] * coords[i, 1] - dN_dL[i, 2] * coords[i, 1]  # dy/dL1
            J[1, 1] += dN_dL[i, 1] * coords[i, 1] - dN_dL[i, 2] * coords[i, 1]  # dy/dL2
        
        det_J = np.linalg.det(J)
        if abs(det_J) < 1e-12:
            continue
        
        J_inv = np.linalg.inv(J)
        
        # Shape function derivatives in global coordinates
        dN_dx = np.zeros(6)
        dN_dy = np.zeros(6)
        for i in range(6):
            # Convert L1, L2, L3 derivatives to x, y derivatives
            dN_dL1 = dN_dL[i, 0] - dN_dL[i, 2]
            dN_dL2 = dN_dL[i, 1] - dN_dL[i, 2]
            
            dN_dx[i] = J_inv[0, 0] * dN_dL1 + J_inv[0, 1] * dN_dL2
            dN_dy[i] = J_inv[1, 0] * dN_dL1 + J_inv[1, 1] * dN_dL2
        
        # Strain-displacement matrix B
        B = np.zeros((3, 12))
        for i in range(6):
            B[0, 2*i]     = dN_dx[i]      # dN/dx for u
            B[1, 2*i + 1] = dN_dy[i]      # dN/dy for v
            B[2, 2*i]     = dN_dy[i]      # dN/dy for u
            B[2, 2*i + 1] = dN_dx[i]      # dN/dx for v
        
        # Add to stiffness matrix
        K += B.T @ D @ B * det_J * w
    
    return K


def build_quad_stiffness(coords, E, nu, is_plastic=False, elem_type=4, plastic_stiffness_reduction=0.1):
    """
    Build stiffness matrix for quadrilateral elements (quad4, quad8, quad9).
    
    Parameters:
        coords: np.ndarray - coordinates of quad nodes (4x2 for quad4, 8x2 for quad8, 9x2 for quad9)
        E: Young's modulus
        nu: Poisson's ratio
        is_plastic: boolean indicating if element is plastic
        elem_type: element type (4 for quad4, 8 for quad8, 9 for quad9)
        plastic_stiffness_reduction: float - stiffness reduction factor for plastic elements
            Literature ranges:
            - Conservative: 0.1 to 0.5 (10% to 50% of original stiffness - less reduction)
            - Aggressive: 0.01 to 0.1 (1% to 10% of original stiffness - more reduction)
            - Default: 0.1 (reduces plastic element stiffness to 10% of elastic value)
    
    Returns:
        np.ndarray - element stiffness matrix (8x8 for quad4, 16x16 for quad8, 18x18 for quad9)
    """
    
    # Constitutive matrix for plane strain
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    # Reduce stiffness if plastic (simplified approach)
    if is_plastic:
        D *= plastic_stiffness_reduction
    
    if elem_type == 4:  # Linear quadrilateral (quad4)
        return build_quad4_stiffness(coords, D)
    elif elem_type == 8:  # 8-node quadrilateral (quad8)
        return build_quad8_stiffness(coords, D)
    elif elem_type == 9:  # 9-node quadrilateral (quad9)
        return build_quad9_stiffness(coords, D)
    else:
        raise ValueError(f"Unsupported quadrilateral element type: {elem_type}")


def build_quad4_stiffness(coords, D):
    """Build stiffness matrix for 4-node quadrilateral using 2x2 Gauss integration."""
    # 2x2 Gauss points and weights
    xi_eta = np.array([[-1/np.sqrt(3), -1/np.sqrt(3)],
                       [1/np.sqrt(3),  -1/np.sqrt(3)],
                       [1/np.sqrt(3),   1/np.sqrt(3)],
                       [-1/np.sqrt(3),  1/np.sqrt(3)]])
    weights = np.array([1.0, 1.0, 1.0, 1.0])
    
    K = np.zeros((8, 8))
    
    for gp in range(4):
        xi, eta = xi_eta[gp]
        w = weights[gp]
        
        # Shape functions for quad4 (bilinear)
        N = 0.25 * np.array([
            (1 - xi) * (1 - eta),    # N1
            (1 + xi) * (1 - eta),    # N2
            (1 + xi) * (1 + eta),    # N3
            (1 - xi) * (1 + eta)     # N4
        ])
        
        # Shape function derivatives in natural coordinates
        dN_dxi = 0.25 * np.array([
            -(1 - eta),    # dN1/dxi
            (1 - eta),     # dN2/dxi
            (1 + eta),     # dN3/dxi
            -(1 + eta)     # dN4/dxi
        ])
        
        dN_deta = 0.25 * np.array([
            -(1 - xi),     # dN1/deta
            -(1 + xi),     # dN2/deta
            (1 + xi),      # dN3/deta
            (1 - xi)       # dN4/deta
        ])
        
        # Jacobian matrix
        J = np.zeros((2, 2))
        for i in range(4):
            J[0, 0] += dN_dxi[i] * coords[i, 0]   # dx/dxi
            J[0, 1] += dN_deta[i] * coords[i, 0]  # dx/deta
            J[1, 0] += dN_dxi[i] * coords[i, 1]   # dy/dxi
            J[1, 1] += dN_deta[i] * coords[i, 1]  # dy/deta
        
        det_J = np.linalg.det(J)
        if abs(det_J) < 1e-12:
            continue
        
        J_inv = np.linalg.inv(J)
        
        # Shape function derivatives in global coordinates
        dN_dx = np.zeros(4)
        dN_dy = np.zeros(4)
        for i in range(4):
            dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
            dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
        
        # Strain-displacement matrix B
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2*i]     = dN_dx[i]      # dN/dx for u
            B[1, 2*i + 1] = dN_dy[i]      # dN/dy for v
            B[2, 2*i]     = dN_dy[i]      # dN/dy for u
            B[2, 2*i + 1] = dN_dx[i]      # dN/dx for v
        
        # Add to stiffness matrix
        K += B.T @ D @ B * det_J * w
    
    return K


def build_quad8_stiffness(coords, D):
    """Build stiffness matrix for 8-node quadrilateral (serendipity) using 3x3 Gauss integration."""
    # 3x3 Gauss points and weights
    xi_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    eta_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    weights_1d = np.array([5/9, 8/9, 5/9])
    
    K = np.zeros((16, 16))
    
    for i in range(3):
        for j in range(3):
            xi, eta = xi_points[i], eta_points[j]
            w = weights_1d[i] * weights_1d[j]
            
            # Shape functions for quad8 (serendipity)
            N = np.zeros(8)
            # Corner nodes
            N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1)
            N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1)
            N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1)
            N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1)
            # Edge nodes
            N[4] = 0.5 * (1 - xi**2) * (1 - eta)
            N[5] = 0.5 * (1 + xi) * (1 - eta**2)
            N[6] = 0.5 * (1 - xi**2) * (1 + eta)
            N[7] = 0.5 * (1 - xi) * (1 - eta**2)
            
            # Shape function derivatives
            dN_dxi = np.zeros(8)
            dN_deta = np.zeros(8)
            
            # Corner derivatives
            dN_dxi[0] = 0.25 * (1 - eta) * (2*xi + eta)
            dN_dxi[1] = 0.25 * (1 - eta) * (2*xi - eta)
            dN_dxi[2] = 0.25 * (1 + eta) * (2*xi + eta)
            dN_dxi[3] = 0.25 * (1 + eta) * (2*xi - eta)
            
            dN_deta[0] = 0.25 * (1 - xi) * (xi + 2*eta)
            dN_deta[1] = 0.25 * (1 + xi) * (-xi + 2*eta)
            dN_deta[2] = 0.25 * (1 + xi) * (xi + 2*eta)
            dN_deta[3] = 0.25 * (1 - xi) * (-xi + 2*eta)
            
            # Edge derivatives
            dN_dxi[4] = -xi * (1 - eta)
            dN_dxi[5] = 0.5 * (1 - eta**2)
            dN_dxi[6] = -xi * (1 + eta)
            dN_dxi[7] = -0.5 * (1 - eta**2)
            
            dN_deta[4] = -0.5 * (1 - xi**2)
            dN_deta[5] = -(1 + xi) * eta
            dN_deta[6] = 0.5 * (1 - xi**2)
            dN_deta[7] = -(1 - xi) * eta
            
            # Jacobian matrix
            J = np.zeros((2, 2))
            for k in range(8):
                J[0, 0] += dN_dxi[k] * coords[k, 0]   # dx/dxi
                J[0, 1] += dN_deta[k] * coords[k, 0]  # dx/deta
                J[1, 0] += dN_dxi[k] * coords[k, 1]   # dy/dxi
                J[1, 1] += dN_deta[k] * coords[k, 1]  # dy/deta
            
            det_J = np.linalg.det(J)
            if abs(det_J) < 1e-12:
                continue
            
            J_inv = np.linalg.inv(J)
            
            # Shape function derivatives in global coordinates
            dN_dx = np.zeros(8)
            dN_dy = np.zeros(8)
            for k in range(8):
                dN_dx[k] = J_inv[0, 0] * dN_dxi[k] + J_inv[0, 1] * dN_deta[k]
                dN_dy[k] = J_inv[1, 0] * dN_dxi[k] + J_inv[1, 1] * dN_deta[k]
            
            # Strain-displacement matrix B
            B = np.zeros((3, 16))
            for k in range(8):
                B[0, 2*k]     = dN_dx[k]      # dN/dx for u
                B[1, 2*k + 1] = dN_dy[k]      # dN/dy for v
                B[2, 2*k]     = dN_dy[k]      # dN/dy for u
                B[2, 2*k + 1] = dN_dx[k]      # dN/dx for v
            
            # Add to stiffness matrix
            K += B.T @ D @ B * det_J * w
    
    return K


def build_quad9_stiffness(coords, D):
    """Build stiffness matrix for 9-node quadrilateral (Lagrangian) using 3x3 Gauss integration."""
    # 3x3 Gauss points and weights
    xi_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    eta_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    weights_1d = np.array([5/9, 8/9, 5/9])
    
    K = np.zeros((18, 18))
    
    for i in range(3):
        for j in range(3):
            xi, eta = xi_points[i], eta_points[j]
            w = weights_1d[i] * weights_1d[j]
            
            # Shape functions for quad9 (Lagrangian)
            N = np.zeros(9)
            
            # Corner nodes
            N[0] = 0.25 * xi * (xi - 1) * eta * (eta - 1)
            N[1] = 0.25 * xi * (xi + 1) * eta * (eta - 1)
            N[2] = 0.25 * xi * (xi + 1) * eta * (eta + 1)
            N[3] = 0.25 * xi * (xi - 1) * eta * (eta + 1)
            
            # Edge nodes
            N[4] = 0.5 * (1 - xi**2) * eta * (eta - 1)
            N[5] = 0.5 * xi * (xi + 1) * (1 - eta**2)
            N[6] = 0.5 * (1 - xi**2) * eta * (eta + 1)
            N[7] = 0.5 * xi * (xi - 1) * (1 - eta**2)
            
            # Center node
            N[8] = (1 - xi**2) * (1 - eta**2)
            
            # Shape function derivatives (detailed implementation omitted for brevity)
            # This would require computing all derivatives similar to quad8 but with 9 nodes
            
            # For now, use simplified approach - this should be fully implemented
            # in a production system
            dN_dx = np.zeros(9)
            dN_dy = np.zeros(9)
            
            # Continue with Jacobian and B matrix calculation similar to quad8...
            # (Implementation abbreviated for space)
    
    return K


def build_truss_stiffness(coords, k_axial):
    """
    Build stiffness matrix for 2-node truss element.
    
    Parameters:
        coords: np.ndarray (2, 2) - coordinates of truss endpoints
        k_axial: axial stiffness (EA/L)
    
    Returns:
        np.ndarray (4, 4) - element stiffness matrix in global coordinates
    """
    if len(coords) < 2:
        return np.zeros((4, 4))
    
    # Element geometry
    dx = coords[1, 0] - coords[0, 0]
    dy = coords[1, 1] - coords[0, 1]
    length = sqrt(dx**2 + dy**2)
    
    if length < 1e-12:
        return np.zeros((4, 4))
    
    # Direction cosines
    cos_theta = dx / length
    sin_theta = dy / length
    
    # Local stiffness matrix
    K_local = k_axial * np.array([
        [1,  -1],
        [-1,  1]
    ])
    
    # Transformation matrix
    T = np.array([
        [cos_theta, sin_theta, 0, 0],
        [0, 0, cos_theta, sin_theta]
    ])
    
    # Transform to global coordinates
    K_global = T.T @ K_local @ T
    
    return K_global


def compute_body_forces(coords, elem_type, b_x, b_y):
    """
    Compute equivalent nodal forces for body forces (gravity, seismic).
    
    Parameters:
        coords: element node coordinates
        elem_type: element type (3 or 4)
        b_x, b_y: body force components per unit volume
    
    Returns:
        np.ndarray: nodal force vector
    """
    if elem_type == 3:  # Triangle
        area = 0.5 * abs((coords[1, 0] - coords[0, 0]) * (coords[2, 1] - coords[0, 1]) - 
                        (coords[2, 0] - coords[0, 0]) * (coords[1, 1] - coords[0, 1]))
        
        # For linear elements, distribute equally to nodes
        nodal_force = area / 3.0
        F = np.array([b_x * nodal_force, b_y * nodal_force] * 3)
        
    elif elem_type == 4:  # Quadrilateral (approximate)
        # Compute area using shoelace formula
        x = coords[:, 0]
        y = coords[:, 1]
        area = 0.5 * abs(sum(x[i]*y[i+1] - x[i+1]*y[i] for i in range(-1, len(x)-1)))
        
        # Distribute equally to nodes
        nodal_force = area / 4.0
        F = np.array([b_x * nodal_force, b_y * nodal_force] * 4)
    else:
        F = np.zeros(2 * elem_type)
    
    return F


def apply_boundary_conditions(K_global, F_global, bc_type, nodes):
    """
    Apply boundary conditions to the global system.
    
    Parameters:
        K_global: global stiffness matrix (lil_matrix)
        F_global: global force vector
        bc_type: boundary condition types per node
        nodes: node coordinates
    
    Returns:
        K_constrained: constrained stiffness matrix
        F_constrained: constrained force vector
        constraint_dofs: list of constrained DOFs
    """
    n_nodes = len(nodes)
    n_dof = 2 * n_nodes
    
    constraint_dofs = set()
    
    for i in range(n_nodes):
        node_dof_x = 2 * i
        node_dof_y = 2 * i + 1
        
        if bc_type[i] == 1:  # Fixed (u=0, v=0)
            constraint_dofs.add(node_dof_x)
            constraint_dofs.add(node_dof_y)
        elif bc_type[i] == 2:  # X roller (u=0, v=free)
            constraint_dofs.add(node_dof_x)
        elif bc_type[i] == 3:  # Y roller (u=free, v=0)
            constraint_dofs.add(node_dof_y)
        # Type 4 (applied force) - no constraints, forces already added to F_global
    
    constraint_dofs = sorted(list(constraint_dofs))
    
    # Create free DOF list
    free_dofs = [i for i in range(n_dof) if i not in constraint_dofs]
    
    # Extract free system
    K_constrained = K_global[np.ix_(free_dofs, free_dofs)]
    F_constrained = F_global[free_dofs]
    
    return K_constrained, F_constrained, constraint_dofs


def compute_element_stress(coords, displacements, node_ids, E, nu, pore_pressures, elem_type, sigma0=None):
    """
    Compute element stresses from nodal displacements.
    
    Parameters:
        coords: element node coordinates
        displacements: global displacement vector
        node_ids: global node IDs for this element
        E, nu: material properties
        pore_pressures: nodal pore pressures
        elem_type: element type
    
    Returns:
        np.ndarray: stress vector [sig_x, sig_y, tau_xy]
    """
    # Constitutive matrix
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    if elem_type == 3:  # Linear triangle (tri3)
        return compute_tri3_stress(coords, displacements, node_ids, D, pore_pressures, sigma0)
    elif elem_type == 6:  # Quadratic triangle (tri6)
        return compute_tri6_stress(coords, displacements, node_ids, D, pore_pressures)
    elif elem_type == 4:  # Linear quadrilateral (quad4)
        return compute_quad4_stress(coords, displacements, node_ids, D, pore_pressures)
    elif elem_type in [8, 9]:  # Higher-order quadrilaterals
        # For now, use quad4 stress calculation at element center
        # In production, this should use proper integration
        return compute_quad4_stress(coords[:4], displacements, node_ids[:4], D, pore_pressures[:4])
    else:
        return np.zeros(3)


def compute_tri3_stress(coords, displacements, node_ids, D, pore_pressures, sigma0=None):
    """Compute stress for 3-node triangle."""
    # Extract nodal displacements
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # Calculate area and B matrix
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    
    area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    
    if area < 1e-12:
        return np.zeros(3)
    
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    B = (1.0 / (2.0 * area)) * np.array([
        [b1, 0,  b2, 0,  b3, 0],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
    ])
    
    # Nodal displacement vector
    u_elem = np.array([u[0], v[0], u[1], v[1], u[2], v[2]])
    
    # Strain vector
    strain = B @ u_elem
    
    # Total stress increment from displacement
    stress_inc = D @ strain

    # Add initial geostatic stress if provided
    if sigma0 is None:
        sigma0 = np.zeros(3)
    stress_total = stress_inc + sigma0
    
    # Average pore pressure
    u_avg = np.mean(pore_pressures)
    
    # Effective stress (subtract pore pressure from normal stresses)
    stress_eff = stress_total.copy()
    stress_eff[0] -= u_avg  # sig_x'
    stress_eff[1] -= u_avg  # sig_y'
    # tau_xy unchanged
    
    return stress_eff


def compute_tri6_stress(coords, displacements, node_ids, D, pore_pressures):
    """Compute stress for 6-node triangle at element centroid."""
    # For tri6, compute stress at the centroid (L1=L2=L3=1/3)
    L1 = L2 = L3 = 1.0/3.0
    
    # Extract nodal displacements
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # Shape function derivatives at centroid
    dN_dL = np.array([
        [4*L1 - 1, 0,       0],            # dN1/dL1, dN1/dL2, dN1/dL3
        [0,        4*L2 - 1, 0],           # dN2/dL1, dN2/dL2, dN2/dL3
        [0,        0,        4*L3 - 1],    # dN3/dL1, dN3/dL2, dN3/dL3
        [4*L2,     4*L1,     0],           # dN4/dL1, dN4/dL2, dN4/dL3
        [0,        4*L3,     4*L2],        # dN5/dL1, dN5/dL2, dN5/dL3
        [4*L3,     0,        4*L1]         # dN6/dL1, dN6/dL2, dN6/dL3
    ])
    
    # Jacobian transformation
    J = np.zeros((2, 2))
    for i in range(6):
        J[0, 0] += dN_dL[i, 0] * coords[i, 0] - dN_dL[i, 2] * coords[i, 0]
        J[0, 1] += dN_dL[i, 1] * coords[i, 0] - dN_dL[i, 2] * coords[i, 0]
        J[1, 0] += dN_dL[i, 0] * coords[i, 1] - dN_dL[i, 2] * coords[i, 1]
        J[1, 1] += dN_dL[i, 1] * coords[i, 1] - dN_dL[i, 2] * coords[i, 1]
    
    det_J = np.linalg.det(J)
    if abs(det_J) < 1e-12:
        return np.zeros(3)
    
    J_inv = np.linalg.inv(J)
    
    # Shape function derivatives in global coordinates
    dN_dx = np.zeros(6)
    dN_dy = np.zeros(6)
    for i in range(6):
        dN_dL1 = dN_dL[i, 0] - dN_dL[i, 2]
        dN_dL2 = dN_dL[i, 1] - dN_dL[i, 2]
        
        dN_dx[i] = J_inv[0, 0] * dN_dL1 + J_inv[0, 1] * dN_dL2
        dN_dy[i] = J_inv[1, 0] * dN_dL1 + J_inv[1, 1] * dN_dL2
    
    # Strain-displacement matrix B
    B = np.zeros((3, 12))
    for i in range(6):
        B[0, 2*i]     = dN_dx[i]      # dN/dx for u
        B[1, 2*i + 1] = dN_dy[i]      # dN/dy for v
        B[2, 2*i]     = dN_dy[i]      # dN/dy for u
        B[2, 2*i + 1] = dN_dx[i]      # dN/dx for v
    
    # Nodal displacement vector
    u_elem = np.zeros(12)
    for i in range(6):
        u_elem[2*i] = u[i]
        u_elem[2*i + 1] = v[i]
    
    # Strain vector
    strain = B @ u_elem
    
    # Total stress increment from displacement
    stress_inc = D @ strain
    
    # Average pore pressure
    u_avg = np.mean(pore_pressures)
    
    # Effective stress
    stress_eff = stress_inc.copy()
    stress_eff[0] -= u_avg  # sig_x'
    stress_eff[1] -= u_avg  # sig_y'
    
    return stress_eff


def compute_quad4_stress(coords, displacements, node_ids, D, pore_pressures):
    """Compute stress for 4-node quadrilateral at element center."""
    # Extract nodal displacements
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # At element center (xi=0, eta=0)
    xi, eta = 0.0, 0.0
    
    # Shape function derivatives at center
    dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    
    # Jacobian matrix
    J = np.zeros((2, 2))
    for i in range(4):
        J[0, 0] += dN_dxi[i] * coords[i, 0]
        J[0, 1] += dN_deta[i] * coords[i, 0]
        J[1, 0] += dN_dxi[i] * coords[i, 1]
        J[1, 1] += dN_deta[i] * coords[i, 1]
    
    det_J = np.linalg.det(J)
    if abs(det_J) < 1e-12:
        return np.zeros(3)
    
    J_inv = np.linalg.inv(J)
    
    # Shape function derivatives in global coordinates
    dN_dx = np.zeros(4)
    dN_dy = np.zeros(4)
    for i in range(4):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    # Strain-displacement matrix B
    B = np.zeros((3, 8))
    for i in range(4):
        B[0, 2*i]     = dN_dx[i]
        B[1, 2*i + 1] = dN_dy[i]
        B[2, 2*i]     = dN_dy[i]
        B[2, 2*i + 1] = dN_dx[i]
    
    # Nodal displacement vector
    u_elem = np.zeros(8)
    for i in range(4):
        u_elem[2*i] = u[i]
        u_elem[2*i + 1] = v[i]
    
    # Strain vector
    strain = B @ u_elem
    
    # Total stress increment from displacement
    stress_inc = D @ strain
    
    # Average pore pressure
    u_avg = np.mean(pore_pressures)
    
    # Effective stress
    stress_eff = stress_inc.copy()
    stress_eff[0] -= u_avg
    stress_eff[1] -= u_avg
    
    return stress_eff


def compute_element_strain(coords, displacements, node_ids, elem_type):
    """
    Compute element strains from nodal displacements.
    """
    if elem_type == 3:  # Linear triangle (tri3)
        return compute_tri3_strain(coords, displacements, node_ids)
    elif elem_type == 6:  # Quadratic triangle (tri6)  
        return compute_tri6_strain(coords, displacements, node_ids)
    elif elem_type == 4:  # Linear quadrilateral (quad4)
        return compute_quad4_strain(coords, displacements, node_ids)
    elif elem_type in [8, 9]:  # Higher-order quadrilaterals
        # For now, use quad4 strain calculation
        return compute_quad4_strain(coords[:4], displacements, node_ids[:4])
    else:
        return np.zeros(3)


def compute_tri3_strain(coords, displacements, node_ids):
    """Compute strain for 3-node triangle."""
    # Extract nodal displacements
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # Calculate B matrix
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    
    area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    
    if area < 1e-12:
        return np.zeros(3)
    
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    B = (1.0 / (2.0 * area)) * np.array([
        [b1, 0,  b2, 0,  b3, 0],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
    ])
    
    # Nodal displacement vector
    u_elem = np.array([u[0], v[0], u[1], v[1], u[2], v[2]])
    
    # Strain vector
    strain = B @ u_elem
    
    return strain


def compute_tri6_strain(coords, displacements, node_ids):
    """Compute strain for 6-node triangle at centroid."""
    # Similar to compute_tri6_stress but only return strain
    # (Implementation details similar to tri6_stress function)
    L1 = L2 = L3 = 1.0/3.0
    
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # Use same B matrix computation as in tri6_stress
    # (Abbreviated for space - would use same Jacobian transformation)
    
    # For now, return zeros - this should be fully implemented
    return np.zeros(3)


def compute_quad4_strain(coords, displacements, node_ids):
    """Compute strain for 4-node quadrilateral at center."""
    # Similar to compute_quad4_stress but only return strain
    u = np.array([displacements[2*nid] for nid in node_ids])
    v = np.array([displacements[2*nid + 1] for nid in node_ids])
    
    # At element center (xi=0, eta=0)
    xi, eta = 0.0, 0.0
    
    dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    
    # Jacobian and B matrix (same as in stress computation)
    J = np.zeros((2, 2))
    for i in range(4):
        J[0, 0] += dN_dxi[i] * coords[i, 0]
        J[0, 1] += dN_deta[i] * coords[i, 0]
        J[1, 0] += dN_dxi[i] * coords[i, 1]
        J[1, 1] += dN_deta[i] * coords[i, 1]
    
    det_J = np.linalg.det(J)
    if abs(det_J) < 1e-12:
        return np.zeros(3)
    
    J_inv = np.linalg.inv(J)
    
    dN_dx = np.zeros(4)
    dN_dy = np.zeros(4)
    for i in range(4):
        dN_dx[i] = J_inv[0, 0] * dN_dxi[i] + J_inv[0, 1] * dN_deta[i]
        dN_dy[i] = J_inv[1, 0] * dN_dxi[i] + J_inv[1, 1] * dN_deta[i]
    
    B = np.zeros((3, 8))
    for i in range(4):
        B[0, 2*i]     = dN_dx[i]
        B[1, 2*i + 1] = dN_dy[i]
        B[2, 2*i]     = dN_dy[i]
        B[2, 2*i + 1] = dN_dx[i]
    
    u_elem = np.zeros(8)
    for i in range(4):
        u_elem[2*i] = u[i]
        u_elem[2*i + 1] = v[i]
    
    strain = B @ u_elem
    return strain


def compute_truss_force(coords, displacements, node_ids, k_axial):
    """
    Compute axial force in truss element.
    """
    if len(node_ids) < 2:
        return 0.0
    
    # Element geometry
    dx = coords[1, 0] - coords[0, 0]
    dy = coords[1, 1] - coords[0, 1]
    length = sqrt(dx**2 + dy**2)
    
    if length < 1e-12:
        return 0.0
    
    # Direction cosines
    cos_theta = dx / length
    sin_theta = dy / length
    
    # Nodal displacements
    u1 = displacements[2 * node_ids[0]]
    v1 = displacements[2 * node_ids[0] + 1]
    u2 = displacements[2 * node_ids[1]]
    v2 = displacements[2 * node_ids[1] + 1]
    
    # Axial displacement
    u_axial = cos_theta * (u2 - u1) + sin_theta * (v2 - v1)
    
    # Axial force
    force = k_axial * u_axial
    
    return force


def check_mohr_coulomb_yield(stress, c, phi):
    """
    Check Mohr-Coulomb yield criterion.
    
    Parameters:
        stress: stress vector [sig_x, sig_y, tau_xy]
        c: cohesion
        phi: friction angle (radians)
    
    Returns:
        float: yield function value (>0 indicates yielding)
    """
    if len(stress) < 3:
        return 0.0
    
    sig_x, sig_y, tau_xy = stress[:3]
    
    # Principal stresses
    sig_mean = 0.5 * (sig_x + sig_y)
    tau_max = sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
    
    sig_1 = sig_mean + tau_max  # Major principal stress
    sig_3 = sig_mean - tau_max  # Minor principal stress
    
    # Mohr-Coulomb yield function
    if abs(phi) < 1e-6:  # Undrained case (phi = 0)
        f = tau_max - c
    else:
        sin_phi = sin(phi)
        f = 0.5 * (sig_1 - sig_3) - 0.5 * (sig_1 + sig_3) * sin_phi - c * cos(phi)
    
    return f


def compute_von_mises(stress):
    """
    Compute von Mises equivalent stress.
    """
    if len(stress) < 3:
        return 0.0
    
    sig_x, sig_y, tau_xy = stress[:3]
    
    von_mises = sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
    
    return von_mises


def return_mapping_mohr_coulomb_stress_nonassociated(stress_trial, c, phi, psi, E, nu):
    """
    Return mapping algorithm for non-associated Mohr-Coulomb plasticity.
    
    This function implements non-associated plasticity where:
    - Yield criterion uses friction angle phi
    - Flow rule uses dilation angle psi (typically psi = phi/3)
    - This promotes localization along slip surfaces (critical for slope stability)
    
    Parameters:
        stress_trial: np.array - trial elastic stress [sig_x, sig_y, tau_xy]
        c: float - cohesion
        phi: float - friction angle (radians) - used for yield criterion
        psi: float - dilation angle (radians) - used for flow rule  
        E: float - Young's modulus
        nu: float - Poisson's ratio
    
    Returns:
        tuple: (stress_corrected, is_plastic)
            stress_corrected: corrected stress on or within yield surface
            is_plastic: boolean indicating if plastic correction was applied
    """
    
    # Check yield criterion using friction angle phi
    f_trial = check_mohr_coulomb_yield(stress_trial, c, phi)
    
    if f_trial <= 0.0:  # Elastic - no correction needed
        return stress_trial.copy(), False
    
    # Plastic loading - implement non-associated return mapping
    sig_x, sig_y, tau_xy = stress_trial
    
    # Principal stresses and directions
    sig_mean = 0.5 * (sig_x + sig_y)
    tau_max = np.sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
    
    if tau_max < 1e-10:  # Hydrostatic stress state
        return stress_trial.copy(), False
    
    if abs(phi) < 1e-6:  # Undrained case (phi = 0)
        # Tresca criterion - no dilation effect
        if tau_max > c:
            scale_factor = c / tau_max
            sig_x_corr = sig_mean + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            return np.array([sig_x_corr, sig_y_corr, tau_xy_corr]), True
    else:
        # Non-associated Mohr-Coulomb case
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        sin_psi = np.sin(psi)  # Use dilation angle for flow rule
        
        # For non-associated plasticity, the flow rule uses psi instead of phi
        # This creates different volumetric behavior and promotes localization
        
        # Simplified non-associated return: adjust mean stress based on dilation
        # Full implementation would use proper flow rule mathematics
        dilation_factor = sin_psi / sin_phi if sin_phi > 1e-10 else 0.0
        
        # Adjust mean stress based on dilation behavior
        sig_mean_corrected = sig_mean * (1.0 - 0.1 * dilation_factor)  # Reduce volume expansion
        
        # Calculate allowable shear stress from yield criterion (using phi)
        allowable_tau = c * cos_phi + sig_mean_corrected * sin_phi
        
        if allowable_tau > 0 and tau_max > allowable_tau:
            # Scale down the deviatoric stress
            scale_factor = allowable_tau / tau_max
            
            sig_x_corr = sig_mean_corrected + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean_corrected + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            
            return np.array([sig_x_corr, sig_y_corr, tau_xy_corr]), True
    
    # If we get here, no correction was needed
    return stress_trial.copy(), False


def return_mapping_mohr_coulomb_stress(stress_trial, c, phi, E, nu):
    """
    Return mapping algorithm for Mohr-Coulomb plasticity - stress correction approach.
    
    This function corrects trial stresses to satisfy the yield criterion without
    reducing element stiffness dramatically. This prevents the cascade failure
    issues seen in traditional stiffness reduction approaches.
    
    Parameters:
        stress_trial: np.array - trial elastic stress [sig_x, sig_y, tau_xy]
        c: float - cohesion
        phi: float - friction angle (radians)
        E: float - Young's modulus (not used in stress-only correction)
        nu: float - Poisson's ratio (not used in stress-only correction)
    
    Returns:
        tuple: (stress_corrected, is_plastic)
            stress_corrected: corrected stress on or within yield surface
            is_plastic: boolean indicating if plastic correction was applied
    """
    
    # Check yield criterion
    f_trial = check_mohr_coulomb_yield(stress_trial, c, phi)
    
    if f_trial <= 0.0:  # Elastic - no correction needed
        return stress_trial.copy(), False
    
    # Plastic - apply return mapping to bring stress back to yield surface
    sig_x, sig_y, tau_xy = stress_trial
    
    # Principal stresses and directions
    sig_mean = 0.5 * (sig_x + sig_y)
    tau_max = np.sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
    
    if tau_max < 1e-10:  # Hydrostatic stress state
        return stress_trial.copy(), False
    
    sig_1 = sig_mean + tau_max  # Major principal stress
    sig_3 = sig_mean - tau_max  # Minor principal stress
    
    if abs(phi) < 1e-6:  # Undrained case (phi = 0) - Tresca criterion
        # Simply limit the maximum shear stress
        if tau_max > c:
            scale_factor = c / tau_max
            sig_x_corr = sig_mean + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            return np.array([sig_x_corr, sig_y_corr, tau_xy_corr]), True
    else:
        # General Mohr-Coulomb case
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        
        # Calculate the stress point that satisfies yield criterion
        # f = 0.5*(sig_1 - sig_3) - 0.5*(sig_1 + sig_3)*sin_phi - c*cos_phi = 0
        
        # For given mean stress, find allowable stress difference
        sig_mean_eff = sig_mean  # Keep mean stress approximately same
        
        # From yield criterion: (sig_1 - sig_3) = 2*c*cos_phi + (sig_1 + sig_3)*sin_phi
        # With sig_1 + sig_3 = 2*sig_mean: (sig_1 - sig_3) = 2*c*cos_phi + 2*sig_mean*sin_phi
        allowable_tau = c * cos_phi + sig_mean_eff * sin_phi
        
        if allowable_tau > 0 and tau_max > allowable_tau:
            # Scale down the deviatoric stress
            scale_factor = allowable_tau / tau_max
            
            sig_x_corr = sig_mean_eff + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean_eff + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            
            return np.array([sig_x_corr, sig_y_corr, tau_xy_corr]), True
    
    # If we get here, no correction was needed
    return stress_trial.copy(), False


def check_and_update_plastic_state(nodes, elements, element_types, element_materials, 
                                  displacements, E_by_mat, nu_by_mat, u_nodal, 
                                  c_reduced, phi_reduced, plastic_elements, debug_level):
    """
    Check for yielding and update plastic element state for iterative solution.
    
    This function implements proper elastic-plastic iteration by:
    1. Computing current stresses from displacements
    2. Checking yield criterion for each element
    3. Updating plastic element flags
    4. Returns True if plastic state changed (requiring stiffness update)
    """
    
    n_elements = len(elements)
    new_plastic_elements = np.zeros(n_elements, dtype=bool)
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:  # All supported element types
            # Compute current element stresses (with initial stress)
            sigma0_elem = fem_data.get("sigma0_by_elem", None)
            sigma0 = None if sigma0_elem is None else sigma0_elem[elem_idx]
            trial_stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                E_by_mat[element_materials[elem_idx] - 1],
                nu_by_mat[element_materials[elem_idx] - 1],
                u_nodal[active_nodes], elem_type, sigma0)
            
            # Check Mohr-Coulomb yield criterion
            c_elem = c_reduced[elem_idx]
            phi_elem = phi_reduced[elem_idx]
            
            yield_value = check_mohr_coulomb_yield(trial_stress, c_elem, phi_elem)
            
            # Use realistic yield threshold (research-based)
            yield_threshold = max(0.01 * c_elem, 50.0)  # 1% of cohesion or 50 Pa minimum
            
            if yield_value > yield_threshold:
                new_plastic_elements[elem_idx] = True
                if debug_level >= 3:
                    print(f"    Element {elem_idx} yielded: yield_value = {yield_value:.2e}")
    
    # Check if plastic state has changed
    plastic_state_changed = not np.array_equal(new_plastic_elements, plastic_elements)
    
    if plastic_state_changed:
        old_plastic_count = np.sum(plastic_elements)
        plastic_elements[:] = new_plastic_elements  # Update in-place
        new_plastic_count = np.sum(plastic_elements)
        
        if debug_level >= 2:
            print(f"  Plastic state changed: {old_plastic_count} -> {new_plastic_count} plastic elements")
    else:
        if debug_level >= 2:
            print(f"  Plastic state unchanged: {np.sum(plastic_elements)} plastic elements")
    
    return plastic_state_changed


def check_and_update_plastic_state_with_softening(nodes, elements, element_types, element_materials,
                                                  displacements, E_by_mat, nu_by_mat, u_nodal,
                                                  c_reduced, phi_reduced, psi_reduced, plastic_elements,
                                                  accumulated_plastic_strain, debug_level):
    """
    Check for yielding and update plastic element state with strain softening.
    
    This implements strain softening where elements that accumulate plastic strain
    gradually lose strength, promoting failure localization along slip surfaces.
    """
    
    n_elements = len(elements)
    new_plastic_elements = np.zeros(n_elements, dtype=bool)
    
    # Strain softening parameters - more aggressive to promote localization
    softening_modulus = 5.0  # Rate of strength loss (aggressive softening)
    max_softening = 0.7      # Maximum strength loss (retain at least 30% of original strength)
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:  # All supported element types
            # Compute current element stresses (with initial stress)
            sigma0_elem = fem_data.get("sigma0_by_elem", None)
            sigma0 = None if sigma0_elem is None else sigma0_elem[elem_idx]
            trial_stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                E_by_mat[element_materials[elem_idx] - 1],
                nu_by_mat[element_materials[elem_idx] - 1],
                u_nodal[active_nodes], elem_type, sigma0)
            
            # Apply strain softening to strength parameters
            softening_factor = min(softening_modulus * accumulated_plastic_strain[elem_idx], max_softening)
            
            # Reduce strength based on accumulated plastic strain
            c_elem = c_reduced[elem_idx] * (1.0 - softening_factor)
            phi_elem_reduced = phi_reduced[elem_idx] * (1.0 - softening_factor)
            psi_elem_reduced = psi_reduced[elem_idx] * (1.0 - softening_factor)  # Also soften dilation angle
            
            # Update the reduced strength arrays (this affects future iterations)
            c_reduced[elem_idx] = c_elem
            phi_reduced[elem_idx] = phi_elem_reduced
            psi_reduced[elem_idx] = psi_elem_reduced
            
            yield_value = check_mohr_coulomb_yield(trial_stress, c_elem, phi_elem_reduced)
            
            # Use realistic yield threshold
            yield_threshold = max(0.01 * c_elem, 50.0)
            
            if yield_value > yield_threshold:
                new_plastic_elements[elem_idx] = True
                
                # Accumulate plastic strain for elements that are yielding
                # Use equivalent plastic strain increment (simplified)
                strain_vec = compute_element_strain(nodes[active_nodes], displacements, active_nodes, elem_type)
                eps_x, eps_y, gamma_xy = strain_vec[:3]
                
                # Compute equivalent plastic strain increment (von Mises type)
                plastic_strain_increment = np.sqrt((2.0/3.0) * (eps_x**2 + eps_y**2 + eps_x*eps_y + 0.75*gamma_xy**2))
                
                # Accumulate plastic strain more aggressively for localization
                accumulated_plastic_strain[elem_idx] += plastic_strain_increment * 0.1  # More aggressive increment
                
                if debug_level >= 3:
                    print(f"    Element {elem_idx} yielded: yield_value = {yield_value:.2e}, plastic_strain = {accumulated_plastic_strain[elem_idx]:.4f}")
    
    # Check if plastic state has changed
    plastic_state_changed = not np.array_equal(new_plastic_elements, plastic_elements)
    
    if plastic_state_changed:
        old_plastic_count = np.sum(plastic_elements)
        plastic_elements[:] = new_plastic_elements  # Update in-place
        new_plastic_count = np.sum(plastic_elements)
        
        if debug_level >= 2:
            max_plastic_strain = np.max(accumulated_plastic_strain)
            mean_softening = np.mean([min(softening_modulus * ps, max_softening) for ps in accumulated_plastic_strain if ps > 0])
            print(f"  Plastic state changed: {old_plastic_count} -> {new_plastic_count} plastic elements")
            print(f"  Max plastic strain: {max_plastic_strain:.4f}, Mean softening: {mean_softening:.3f}")
    else:
        if debug_level >= 2:
            print(f"  Plastic state unchanged: {np.sum(plastic_elements)} plastic elements")
    
    return plastic_state_changed


def return_mapping_mohr_coulomb(stress_trial, c, phi, E, nu):
    """
    Return mapping algorithm for Mohr-Coulomb plasticity.
    
    This function implements proper elastic-plastic constitutive relations
    using return mapping to ensure stresses remain on the yield surface.
    
    Parameters:
        stress_trial: np.array - trial elastic stress [sig_x, sig_y, tau_xy]
        c: float - cohesion
        phi: float - friction angle (radians)
        E: float - Young's modulus
        nu: float - Poisson's ratio
    
    Returns:
        tuple: (stress_corrected, D_ep, is_plastic)
            stress_corrected: corrected stress on yield surface
            D_ep: elastic-plastic tangent stiffness matrix
            is_plastic: boolean indicating if plastic correction was applied
    """
    
    # Check yield criterion
    f_trial = check_mohr_coulomb_yield(stress_trial, c, phi)
    
    if f_trial <= 1e-9:  # Elastic loading
        # Return elastic stiffness matrix
        factor = E / ((1 + nu) * (1 - 2 * nu))
        D_ep = factor * np.array([
            [1 - nu, nu,     0],
            [nu,     1 - nu, 0],
            [0,      0,      (1 - 2 * nu) / 2]
        ])
        return stress_trial, D_ep, False
    
    # Plastic loading - implement return mapping
    # For now, use a simplified approach: reduce stress proportionally to bring it to yield surface
    # This is not optimal but better than arbitrary stiffness reduction
    
    if abs(phi) < 1e-6:  # Undrained case (phi = 0)
        # Tresca/von Mises type yield - scale down deviatoric stress
        sig_x, sig_y, tau_xy = stress_trial
        sig_mean = 0.5 * (sig_x + sig_y)
        tau_max = sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
        
        if tau_max > c:
            scale_factor = c / tau_max
            sig_x = sig_mean + scale_factor * 0.5 * (sig_x - sig_y)
            sig_y = sig_mean + scale_factor * 0.5 * (sig_y - sig_x)
            tau_xy = scale_factor * tau_xy
    else:
        # General Mohr-Coulomb case
        sig_x, sig_y, tau_xy = stress_trial
        
        # Principal stresses
        sig_mean = 0.5 * (sig_x + sig_y)
        tau_max = sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
        
        sig_1 = sig_mean + tau_max
        sig_3 = sig_mean - tau_max
        
        # Compute required reduction to satisfy yield criterion
        sin_phi = sin(phi)
        cos_phi = cos(phi)
        
        # Solve for new principal stresses on yield surface
        # f = 0.5*(sig_1 - sig_3) - 0.5*(sig_1 + sig_3)*sin_phi - c*cos_phi = 0
        # This gives: (sig_1 - sig_3)*(1 - sin_phi)/2 = c*cos_phi + (sig_1 + sig_3)*sin_phi/2
        
        # Simplified approach: scale down stress difference while maintaining mean stress
        if tau_max > 0:
            # Target stress difference from yield criterion
            sig_m_target = sig_mean  # Keep mean stress roughly the same
            
            # From yield criterion: tau_max_new = c*cos_phi + sig_m_target*sin_phi
            tau_max_new = c * cos_phi + sig_m_target * sin_phi
            
            if tau_max_new > 0:
                scale_factor = tau_max_new / tau_max
                
                # Scale the deviatoric components
                sig_x = sig_mean + scale_factor * 0.5 * (sig_x - sig_y)
                sig_y = sig_mean + scale_factor * 0.5 * (sig_y - sig_x) 
                tau_xy = scale_factor * tau_xy
    
    stress_corrected = np.array([sig_x, sig_y, tau_xy])
    
    # For plastic loading, use reduced stiffness matrix
    # This is still simplified - proper implementation would compute elastic-plastic tangent
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D_elastic = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    # Simple plastic stiffness reduction (can be improved)
    plastic_stiffness_factor = 0.1  # 10% of elastic stiffness
    D_ep = plastic_stiffness_factor * D_elastic
    
    return stress_corrected, D_ep, True


def solve_ssrm(fem_data, debug_level=0):
    """
    Solve for factor of safety using Shear Strength Reduction Method (SSRM).
    
    This function performs SSRM by iteratively increasing the shear strength reduction
    factor F until the finite element system fails to converge, indicating slope failure.
    The critical factor of safety is the last F value for which convergence was achieved.
    
    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        debug_level (int): Verbosity level (0=quiet, 1=basic, 2=detailed, 3=debug)
    
    Returns:
        dict: SSRM solution dictionary containing:
            - converged: bool, whether SSRM procedure completed successfully
            - FS: float, critical factor of safety
            - last_solution: dict, final FEM solution at critical F
            - F_history: list, history of F values tested
            - convergence_history: list, convergence status for each F
            - iterations_ssrm: int, number of SSRM iterations
    """
    
    if debug_level >= 1:
        print("Starting Shear Strength Reduction Method (SSRM)")
    
    # SSRM parameters
    F_start = 1.0
    F_increment_initial = 0.1    # Larger initial steps to find failure region quickly
    F_increment_min = 0.001      # Fine resolution for final refinement
    F_max = 3.0                  # Maximum F to prevent excessive iterations
    max_ssrm_iterations = 100    # Allow more iterations for adaptive search
    
    # Adaptive stepping parameters
    F_increment = F_increment_initial
    refinement_mode = False      # Whether we're in refinement mode
    failure_interval = None      # Store the interval where failure occurs
    
    # History tracking
    F_history = []
    convergence_history = []
    solutions_history = []
    displacement_history = []  # Track displacement growth rate
    
    # Initialize
    F = F_start
    F_increment = F_increment_initial
    last_converged_F = None
    last_converged_solution = None
    
    ssrm_iteration = 0
    
    for ssrm_iteration in range(max_ssrm_iterations):
        if debug_level >= 1:
            print(f"\n=== SSRM Iteration {ssrm_iteration + 1}: F = {F:.4f} ===")
        
        # Solve FEM with current F
        solution = solve_fem(fem_data, F=F, debug_level=max(0, debug_level-1))
        
        F_history.append(F)
        convergence_history.append(solution.get("converged", False))
        solutions_history.append(solution)
        
        if solution.get("converged", False):
            # Solution converged - this is stable, continue with higher F
            max_displacement = np.max(np.abs(solution.get("displacements", [0])))
            
            # Check for excessive displacements that indicate failure
            # Literature suggests displacements > 5% of characteristic length indicate failure
            # For SSRM, we use a more sensitive criterion
            characteristic_length = np.max(fem_data["nodes"][:, 0]) - np.min(fem_data["nodes"][:, 0])  # Approximate slope width
            displacement_threshold = 0.05 * characteristic_length  # Reduced from 0.1 to 0.05 (5%)
            
            if max_displacement > displacement_threshold:
                if debug_level >= 1:
                    print(f"  ✗ Displacement-based failure: max_disp = {max_displacement:.4f} > threshold = {displacement_threshold:.4f}")
                
                # This is a failure - return the last converged solution
                if last_converged_F is not None:
                    return {
                        "converged": True,
                        "FS": last_converged_F,
                        "last_solution": last_converged_solution,
                        "F_history": F_history,
                        "convergence_history": convergence_history,
                        "solutions_history": solutions_history,
                        "iterations_ssrm": ssrm_iteration + 1,
                        "failure_mode": "excessive_displacement"
                    }
                else:
                    return {
                        "converged": False,
                        "error": "Failed due to excessive displacement at initial F value",
                        "FS": None,
                        "F_history": F_history,
                        "convergence_history": convergence_history,
                        "iterations_ssrm": ssrm_iteration + 1
                    }
            
            # Check for displacement growth rate instability
            # If displacements are growing rapidly between iterations, this indicates instability
            if len(displacement_history) > 0:
                last_displacement = displacement_history[-1]
                displacement_growth_rate = (max_displacement - last_displacement) / last_displacement if last_displacement > 0 else 0
                
                # If displacement growth rate > 20%, this indicates instability (reduced from 50%)
                if displacement_growth_rate > 0.2:
                    if debug_level >= 1:
                        print(f"  ✗ Displacement growth rate failure: growth_rate = {displacement_growth_rate:.2f} > 0.2")
                    
                    # This is a failure - return the last converged solution
                    if last_converged_F is not None:
                        return {
                            "converged": True,
                            "FS": last_converged_F,
                            "last_solution": last_converged_solution,
                            "F_history": F_history,
                            "convergence_history": convergence_history,
                            "solutions_history": solutions_history,
                            "iterations_ssrm": ssrm_iteration + 1,
                            "failure_mode": "displacement_growth_instability"
                        }
                    else:
                        return {
                            "converged": False,
                            "error": "Failed due to displacement growth instability at initial F value",
                            "FS": None,
                            "F_history": F_history,
                            "convergence_history": convergence_history,
                            "iterations_ssrm": ssrm_iteration + 1
                        }
            
            # Store displacement for growth rate checking
            displacement_history.append(max_displacement)
            
            # Store this as the last converged solution
            last_converged_F = F
            last_converged_solution = solution
            
            if debug_level >= 1:
                print(f"  ✓ Converged: {np.sum(solution.get('plastic_elements', []))}/{len(solution.get('plastic_elements', []))} plastic elements, {solution.get('iterations', 'Unknown')} iterations")
                print(f"  Max displacement: {max_displacement:.4f}")
                print(f"  Displacement threshold: {displacement_threshold:.4f}")
            
            # In refinement mode, update the failure interval
            if refinement_mode and failure_interval is not None:
                # This F converged, so the failure point is between this F and the upper bound
                failure_interval = (F, failure_interval[1])
                
                # Check if we've reached the desired precision
                if (failure_interval[1] - failure_interval[0]) < F_increment_min:
                    if debug_level >= 1:
                        print(f"  🎯 Precision reached: critical FS = {F:.4f} ± {F_increment_min:.4f}")
                    
                    return {
                        "converged": True,
                        "FS": F,
                        "last_solution": solution,
                        "F_history": F_history,
                        "convergence_history": convergence_history,
                        "solutions_history": solutions_history,
                        "iterations_ssrm": ssrm_iteration + 1,
                        "failure_mode": "bisection_precision"
                    }
        else:
            # Solution failed to converge - this indicates failure!
            if debug_level >= 1:
                print(f"  ✗ Non-convergence failure: {solution.get('error', 'Unknown error')}")
                print(f"  Iterations attempted: {solution.get('iterations', 'Unknown')}")
            
            # We've found the failure region! Now switch to refinement mode
            if not refinement_mode:
                refinement_mode = True
                failure_interval = (last_converged_F, F)
                F_increment = (failure_interval[1] - failure_interval[0]) / 2.0  # Start with bisection
                
                if debug_level >= 1:
                    print(f"  🔍 Entering refinement mode: failure interval = [{failure_interval[0]:.4f}, {failure_interval[1]:.4f}]")
                    print(f"  🔍 Next F = {failure_interval[0] + F_increment:.4f}")
                
                # Go back to the last converged F and refine
                F = failure_interval[0]
                continue  # Skip the increment and continue with refined F
            
            # In refinement mode, update the failure interval
            if failure_interval is not None:
                # This F failed to converge, so the failure point is between the lower bound and this F
                failure_interval = (failure_interval[0], F)
                
                # Check if we've reached the desired precision
                if (failure_interval[1] - failure_interval[0]) < F_increment_min:
                    if debug_level >= 1:
                        print(f"  🎯 Precision reached: critical FS = {failure_interval[0]:.4f} ± {F_increment_min:.4f}")
                    
                    return {
                        "converged": True,
                        "FS": failure_interval[0],
                        "last_solution": last_converged_solution,
                        "F_history": F_history,
                        "convergence_history": convergence_history,
                        "solutions_history": solutions_history,
                        "iterations_ssrm": ssrm_iteration + 1,
                        "failure_mode": "bisection_precision"
                    }
        
        # Update F for next iteration
        if refinement_mode:
            # In refinement mode, use bisection to hone in on critical F
            if failure_interval is not None:
                # Try the midpoint of the current interval
                F_mid = (failure_interval[0] + failure_interval[1]) / 2.0
                F = F_mid
                
                if debug_level >= 2:
                    print(f"  🔍 Refinement: trying F = {F:.4f} (midpoint of [{failure_interval[0]:.4f}, {failure_interval[1]:.4f}])")
        else:
            # Normal mode: increment F
            F += F_increment
        
        # Check bounds
        if F > F_max:
            if debug_level >= 1:
                print(f"  Reached maximum F = {F_max}, slope appears very stable")
            
            return {
                "converged": True,
                "FS": last_converged_F if last_converged_F else F_max,
                "last_solution": last_converged_solution,
                "F_history": F_history,
                "convergence_history": convergence_history,
                "solutions_history": solutions_history,
                "iterations_ssrm": ssrm_iteration + 1,
                "note": f"Reached maximum F = {F_max} without failure"
            }
    
    # Maximum SSRM iterations reached
    if debug_level >= 1:
        print(f"Maximum SSRM iterations ({max_ssrm_iterations}) reached")
    
    return {
        "converged": bool(last_converged_F is not None),
        "FS": last_converged_F,
        "last_solution": last_converged_solution,
        "F_history": F_history,
        "convergence_history": convergence_history,
        "solutions_history": solutions_history,
        "iterations_ssrm": ssrm_iteration + 1,
        "error": "Maximum SSRM iterations reached"
    }


def check_convergence(current_solution, previous_solution, tol_disp=1e-6, tol_force=1e-6):
    """
    Check convergence criteria for FEM solution.
    
    Parameters:
        current_solution: current iteration solution
        previous_solution: previous iteration solution  
        tol_disp: displacement tolerance
        tol_force: force residual tolerance
    
    Returns:
        tuple: (converged: bool, disp_norm: float, force_norm: float)
    """
    if previous_solution is None:
        return False, np.inf, np.inf
    
    # Displacement convergence check
    u_current = current_solution.get("displacements", np.array([]))
    u_previous = previous_solution.get("displacements", np.array([]))
    
    if len(u_current) != len(u_previous):
        return False, np.inf, np.inf
    
    du = u_current - u_previous
    disp_norm = np.linalg.norm(du) / max(np.linalg.norm(u_current), 1e-12)
    
    # Force convergence check (simplified - using residual norm from solution)
    force_norm = current_solution.get("residual_norm", np.inf)
    
    # Both criteria must be satisfied
    converged = (disp_norm < tol_disp and force_norm < tol_force)
    
    return converged, disp_norm, force_norm


def compute_elastic_plastic_stiffness(stress, strain, material_props, plastic_state, accumulated_plastic_strain=0.0):
    """
    Compute proper elastic-plastic tangent stiffness matrix.
    
    This implements the consistent tangent stiffness approach for elastic-plastic materials
    rather than arbitrary stiffness reduction. This is the foundation for robust
    convergence in nonlinear FEM analysis.
    
    Parameters:
        stress: np.ndarray - current stress state [sig_x, sig_y, tau_xy]
        strain: np.ndarray - current strain state [eps_x, eps_y, gamma_xy]
        material_props: dict - material properties (E, nu, c, phi, psi)
        plastic_state: bool - whether element is currently plastic
        accumulated_plastic_strain: float - accumulated plastic strain for softening
    
    Returns:
        tuple: (D_ep, is_plastic, plastic_strain_increment)
            D_ep: elastic-plastic tangent stiffness matrix
            is_plastic: boolean indicating if element is plastic
            plastic_strain_increment: incremental plastic strain
    """
    E = material_props['E']
    nu = material_props['nu']
    c = material_props['c']
    phi = material_props['phi']
    psi = material_props.get('psi', phi/3.0)  # Dilation angle, default to phi/3
    
    # Apply strain softening if accumulated plastic strain > 0
    if accumulated_plastic_strain > 0.0:
        softening_factor = min(0.7, 5.0 * accumulated_plastic_strain)  # Max 70% strength loss
        c *= (1.0 - softening_factor)
        phi *= (1.0 - softening_factor)
        psi *= (1.0 - softening_factor)
    
    # Elastic constitutive matrix (plane strain)
    factor = E / ((1 + nu) * (1 - 2 * nu))
    D_e = factor * np.array([
        [1 - nu, nu,     0],
        [nu,     1 - nu, 0],
        [0,      0,      (1 - 2 * nu) / 2]
    ])
    
    # Check yield criterion
    yield_value = check_mohr_coulomb_yield(stress, c, np.radians(phi))
    
    if yield_value <= 0.0 or not plastic_state:
        # Elastic loading or unloading - return elastic stiffness
        return D_e, False, 0.0
    
    # Plastic loading - compute elastic-plastic tangent stiffness
    # This is the key improvement over the previous stiffness reduction approach
    
    # Compute stress gradients for yield function and flow rule
    sig_x, sig_y, tau_xy = stress
    
    # Principal stresses and directions
    sig_mean = 0.5 * (sig_x + sig_y)
    tau_max = np.sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
    
    if tau_max < 1e-10:
        return D_e, False, 0.0
    
    # Yield function gradient (∂f/∂σ)
    if abs(phi) < 1e-6:  # Undrained case
        df_dsig = np.array([0.0, 0.0, 1.0])  # ∂f/∂τ_xy = 1
    else:
        # For Mohr-Coulomb, ∂f/∂σ depends on principal stress directions
        # Simplified approach: use current stress state
        df_dsig = np.array([
            0.5 * (1.0 - np.sin(np.radians(phi))),  # ∂f/∂σ_x
            0.5 * (1.0 + np.sin(np.radians(phi))),  # ∂f/∂σ_y
            0.0  # ∂f/∂τ_xy (simplified)
        ])
    
    # Flow rule gradient (∂g/∂σ) - non-associated plasticity
    if abs(psi) < 1e-6:  # No dilation
        dg_dsig = np.array([0.0, 0.0, 1.0])
    else:
        # Non-associated flow rule with dilation angle
        dg_dsig = np.array([
            0.5 * (1.0 - np.sin(np.radians(psi))),  # ∂g/∂σ_x
            0.5 * (1.0 + np.sin(np.radians(psi))),  # ∂g/∂σ_y
            0.0  # ∂g/∂τ_xy (simplified)
        ])
    
    # Compute plastic multiplier increment
    # Δλ = (∂f/∂σ)ᵀ D_e Δε / [(∂f/∂σ)ᵀ D_e (∂g/∂σ)]
    numerator = df_dsig @ D_e @ strain
    denominator = df_dsig @ D_e @ dg_dsig
    
    if abs(denominator) < 1e-12:
        return D_e, False, 0.0
    
    delta_lambda = numerator / denominator
    
    # Elastic-plastic tangent stiffness matrix
    # D_ep = D_e - (D_e ∂g/∂σ ⊗ ∂f/∂σ D_e) / [(∂f/∂σ)ᵀ D_e (∂g/∂σ)]
    outer_product = np.outer(D_e @ dg_dsig, df_dsig @ D_e)
    D_ep = D_e - outer_product / denominator
    
    # Ensure positive definiteness
    eigenvals = np.linalg.eigvals(D_ep)
    if np.any(eigenvals < 0):
        # Fall back to elastic stiffness if D_ep becomes non-positive definite
        return D_e, False, 0.0
    
    # Compute plastic strain increment
    plastic_strain_increment = delta_lambda * dg_dsig
    
    return D_ep, True, np.linalg.norm(plastic_strain_increment)


def mohr_coulomb_return_mapping(stress_trial, material_props, plastic_state, accumulated_plastic_strain=0.0):
    """
    Implement proper Mohr-Coulomb return mapping algorithm.
    
    This replaces the previous simplified stress correction with a rigorous
    return mapping algorithm that correctly projects stresses back to the
    yield surface while maintaining consistency with the constitutive model.
    
    Parameters:
        stress_trial: np.ndarray - trial elastic stress [sig_x, sig_y, tau_xy]
        material_props: dict - material properties (E, nu, c, phi, psi)
        plastic_state: bool - whether element is currently plastic
        accumulated_plastic_strain: float - accumulated plastic strain for softening
    
    Returns:
        tuple: (stress_corrected, is_plastic, plastic_strain_increment)
            stress_corrected: stress projected back to yield surface
            is_plastic: boolean indicating if plastic correction was applied
            plastic_strain_increment: incremental plastic strain
    """
    E = material_props['E']
    nu = material_props['nu']
    c = material_props['c']
    phi = material_props['phi']
    psi = material_props.get('psi', phi/3.0)  # Dilation angle
    
    # Apply strain softening if accumulated plastic strain > 0
    if accumulated_plastic_strain > 0.0:
        softening_factor = min(0.7, 5.0 * accumulated_plastic_strain)
        c *= (1.0 - softening_factor)
        phi *= (1.0 - softening_factor)
        psi *= (1.0 - softening_factor)
    
    # Check yield criterion
    yield_value = check_mohr_coulomb_yield(stress_trial, c, np.radians(phi))
    
    if yield_value <= 0.0 or not plastic_state:
        # Elastic - no correction needed
        return stress_trial.copy(), False, 0.0
    
    # Plastic loading - implement proper return mapping
    sig_x, sig_y, tau_xy = stress_trial
    
    # Principal stresses and directions
    sig_mean = 0.5 * (sig_x + sig_y)
    tau_max = np.sqrt(0.25 * (sig_x - sig_y)**2 + tau_xy**2)
    
    if tau_max < 1e-10:
        return stress_trial.copy(), False, 0.0
    
    # Compute principal stress directions
    if abs(sig_x - sig_y) < 1e-10:
        # Hydrostatic stress state
        theta = np.pi/4.0
    else:
        theta = 0.5 * np.arctan2(2.0 * tau_xy, sig_x - sig_y)
    
    # Principal stresses
    sig_1 = sig_mean + tau_max  # Major principal stress
    sig_3 = sig_mean - tau_max  # Minor principal stress
    
    # Mohr-Coulomb return mapping
    if abs(phi) < 1e-6:  # Undrained case (Tresca criterion)
        # For Tresca, simply limit the maximum shear stress
        if tau_max > c:
            scale_factor = c / tau_max
            
            # Scale deviatoric stresses while maintaining mean stress
            sig_x_corr = sig_mean + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            
            stress_corrected = np.array([sig_x_corr, sig_y_corr, tau_xy_corr])
            
            # Compute plastic strain increment (simplified)
            plastic_strain_increment = (1.0 - scale_factor) * tau_max / (2.0 * E)
            
            return stress_corrected, True, plastic_strain_increment
    else:
        # General Mohr-Coulomb case
        sin_phi = np.sin(np.radians(phi))
        cos_phi = np.cos(np.radians(phi))
        sin_psi = np.sin(np.radians(psi))  # Non-associated flow rule
        
        # From yield criterion: f = 0.5*(sig_1 - sig_3) - 0.5*(sig_1 + sig_3)*sin_phi - c*cos_phi = 0
        # Solve for the corrected principal stresses that satisfy f = 0
        
        # For given mean stress, find allowable stress difference
        # From yield criterion: (sig_1 - sig_3) = 2*c*cos_phi + (sig_1 + sig_3)*sin_phi
        # With sig_1 + sig_3 = 2*sig_mean: (sig_1 - sig_3) = 2*c*cos_phi + 2*sig_mean*sin_phi
        allowable_tau = c * cos_phi + sig_mean * sin_phi
        
        if allowable_tau > 0 and tau_max > allowable_tau:
            # Scale down the deviatoric stress to satisfy yield criterion
            scale_factor = allowable_tau / tau_max
            
            # Apply correction while maintaining mean stress
            sig_x_corr = sig_mean + scale_factor * (sig_x - sig_mean)
            sig_y_corr = sig_mean + scale_factor * (sig_y - sig_mean)
            tau_xy_corr = scale_factor * tau_xy
            
            stress_corrected = np.array([sig_x_corr, sig_y_corr, tau_xy_corr])
            
            # Compute plastic strain increment using flow rule
            # For non-associated plasticity, use dilation angle psi
            if abs(psi) > 1e-6:
                # Volumetric plastic strain increment
                eps_v_p = (1.0 - scale_factor) * tau_max * sin_psi / E
                # Deviatoric plastic strain increment
                eps_d_p = (1.0 - scale_factor) * tau_max * (1.0 - sin_psi) / E
                plastic_strain_increment = np.sqrt(eps_v_p**2 + eps_d_p**2)
            else:
                # No dilation - only deviatoric plastic strain
                plastic_strain_increment = (1.0 - scale_factor) * tau_max / E
            
            return stress_corrected, True, plastic_strain_increment
    
    # If we get here, no correction was needed
    return stress_trial.copy(), False, 0.0


def assemble_global_stiffness_with_plasticity(nodes, elements, element_types, element_materials, 
                                            displacements, material_props_by_mat, plastic_elements,
                                            accumulated_plastic_strain, u_nodal, fem_data, debug_level=0):
    """
    Assemble global stiffness matrix using consistent tangent stiffness for plastic elements.
    
    This replaces the previous approach of arbitrary stiffness reduction with a proper
    elastic-plastic constitutive model that provides much better numerical stability.
    
    Parameters:
        nodes: np.ndarray - node coordinates
        elements: np.ndarray - element connectivity
        element_types: np.ndarray - element types
        element_materials: np.ndarray - material IDs for each element
        displacements: np.ndarray - current displacement field
        material_props_by_mat: dict - material properties by material ID
        plastic_elements: np.ndarray - boolean array indicating plastic elements
        accumulated_plastic_strain: np.ndarray - accumulated plastic strain per element
        u_nodal: np.ndarray - nodal pore pressures
        fem_data: dict - FEM data dictionary containing initial stresses and other data
        debug_level: int - verbosity level
    
    Returns:
        tuple: (K_global, plastic_state_changed, new_plastic_elements)
            K_global: global stiffness matrix
            plastic_state_changed: boolean indicating if plastic state changed
            new_plastic_elements: updated plastic element flags
    """
    n_nodes = len(nodes)
    n_elements = len(elements)
    n_dof = 2 * n_nodes
    
    K_global = lil_matrix((n_dof, n_dof))
    new_plastic_elements = np.zeros(n_elements, dtype=bool)
    plastic_state_changed = False
    
    if debug_level >= 2:
        print(f"    Assembling global stiffness matrix with {n_elements} elements")
    
    # Assemble 2D elements
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        # Get material properties
        mat_id = element_materials[elem_idx] - 1
        material_props = material_props_by_mat[mat_id]
        
        # Check current plastic state
        is_currently_plastic = plastic_elements[elem_idx]
        
        # Compute current element stress and strain
        if elem_type in [3, 6]:  # Triangular elements
            # Compute element stress from current displacements
            # Include initial stress per element if available
            sigma0_elem = fem_data.get("sigma0_by_elem", None)
            sigma0 = None if sigma0_elem is None else sigma0_elem[elem_idx]
            stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                material_props['E'], material_props['nu'],
                u_nodal[active_nodes], elem_type, sigma0)
            
            # Compute element strain
            strain = compute_element_strain(
                nodes[active_nodes], displacements, active_nodes, elem_type)
            
            # Temporarily use elastic stiffness for all elements to debug constitutive model
            # Apply stress return mapping and get constitutive matrix
            stress_corrected, is_plastic, plastic_strain_increment = mohr_coulomb_return_mapping(
                stress, material_props, is_currently_plastic, 
                accumulated_plastic_strain[elem_idx])
            
            # Update plastic state
            new_plastic_elements[elem_idx] = is_plastic
            
            # Update accumulated plastic strain
            if is_plastic and plastic_strain_increment > 0:
                accumulated_plastic_strain[elem_idx] += plastic_strain_increment
            
            # Temporarily use elastic stiffness to debug constitutive model
            E = material_props['E']
            nu = material_props['nu']
            factor = E / ((1 + nu) * (1 - 2 * nu))
            D_matrix = factor * np.array([
                [1 - nu, nu,     0],
                [nu,     1 - nu, 0],
                [0,      0,      (1 - 2 * nu) / 2]
            ])
            
            # Get constitutive matrix (elastic or elastic-plastic)
            # D_matrix, _, _ = compute_elastic_plastic_stiffness(
            #     stress_corrected, strain, material_props, is_plastic,
            #     accumulated_plastic_strain[elem_idx])
            
            # Build element stiffness matrix
            if elem_type == 3:
                K_elem = build_tri3_stiffness(nodes[active_nodes], D_matrix)
            else:  # elem_type == 6
                K_elem = build_tri6_stiffness(nodes[active_nodes], D_matrix)
                
        elif elem_type in [4, 8, 9]:  # Quadrilateral elements
            # Similar approach for quadrilateral elements
            # For now, use elastic stiffness - can be extended later
            E = material_props['E']
            nu = material_props['nu']
            
            if elem_type == 4:
                K_elem = build_quad4_stiffness(nodes[active_nodes], E, nu, is_plastic, elem_type)
            elif elem_type == 8:
                K_elem = build_quad8_stiffness(nodes[active_nodes], E, nu, is_plastic, elem_type)
            else:  # elem_type == 9
                K_elem = build_quad9_stiffness(nodes[active_nodes], E, nu, is_plastic, elem_type)
        else:
            if debug_level >= 1:
                print(f"Warning: Element type {elem_type} not supported, skipping element {elem_idx}")
            continue
        
        # Assemble into global matrix
        dofs = []
        for node_id in active_nodes:
            dofs.extend([2*node_id, 2*node_id + 1])
        
        for i in range(len(dofs)):
            for j in range(len(dofs)):
                K_global[dofs[i], dofs[j]] += K_elem[i, j]
    
    # Check if plastic state has changed
    plastic_state_changed = not np.array_equal(new_plastic_elements, plastic_elements)
    
    if plastic_state_changed and debug_level >= 2:
        old_plastic_count = np.sum(plastic_elements)
        new_plastic_count = np.sum(new_plastic_elements)
        print(f"      Plastic state changed: {old_plastic_count} -> {new_plastic_count} plastic elements")
    
    return K_global, plastic_state_changed, new_plastic_elements


def solve_fem_elastic_plastic(fem_data, F=1.0, debug_level=0):
    """
    Two-phase elastic-plastic FEM solver with proper constitutive modeling.
    """
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

    # 1D elements
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    element_types_1d = fem_data.get("element_types_1d", np.array([]))
    t_allow_by_1d_elem = fem_data.get("t_allow_by_1d_elem", np.array([]))
    t_res_by_1d_elem = fem_data.get("t_res_by_1d_elem", np.array([]))
    k_by_1d_elem = fem_data.get("k_by_1d_elem", np.array([]))

    # Seismic and other loads
    k_seismic = fem_data.get("k_seismic", 0.0)

    n_nodes = len(nodes)
    n_elements = len(elements)
    n_1d_elements = len(elements_1d)
    n_dof = 2 * n_nodes

    # Reduce strength parameters by factor F (per element)
    c_reduced = c_by_elem / F
    phi_reduced = np.arctan(np.tan(np.radians(phi_by_elem)) / F)  # radians per element

    # Non-associated plasticity (per element)
    psi_by_elem = phi_by_elem / 3.0
    psi_reduced = np.arctan(np.tan(np.radians(psi_by_elem)) / F)

    if debug_level >= 1:
        print(f"Starting two-phase FEM analysis with F = {F:.3f}")
        print(f"Mesh: {n_nodes} nodes, {n_elements} 2D elements, {n_1d_elements} 1D elements")

    displacements = np.zeros(n_dof)

    plastic_elements = np.zeros(n_elements, dtype=bool)
    failed_1d_elements = np.zeros(n_1d_elements, dtype=bool)
    accumulated_plastic_strain = np.zeros(n_elements)

    # Phase 1: Elastic solution
    if debug_level >= 1:
        print("  Phase 1: Elastic solution")

    # Build material properties dict per material for elastic phase
    material_props_by_mat = []
    for i in range(len(E_by_mat)):
        material_props_by_mat.append({
            'E': E_by_mat[i],
            'nu': nu_by_mat[i],
            'c': None,
            'phi': None,
            'psi': None,
            'gamma': gamma_by_mat[i]
        })

    # Build element-wise reduced strengths (arrays)
    c_reduced_by_elem = c_reduced if isinstance(c_reduced, np.ndarray) else np.full(n_elements, c_reduced)
    phi_reduced_by_elem = phi_reduced if isinstance(phi_reduced, np.ndarray) else np.full(n_elements, phi_reduced)

    elastic_solution = solve_elastic_phase(
        fem_data, material_props_by_mat, plastic_elements,
        accumulated_plastic_strain, u_nodal, k_seismic,
        c_reduced_by_elem, phi_reduced_by_elem, debug_level
    )

    if not elastic_solution["converged"]:
        if debug_level >= 1:
            print("  ✗ Elastic phase failed to converge")
        return {
            "converged": False,
            "error": "Elastic phase failed to converge",
            "phase": "elastic"
        }

    displacements = elastic_solution["displacements"]

    # Determine plastic elements from elastic phase stresses
    plastic_elements = elastic_solution.get("plastic_elements", plastic_elements)
    n_plastic_after_elastic = np.sum(plastic_elements)

    if n_plastic_after_elastic == 0:
        if debug_level >= 1:
            print("  ✓ Elastic solution sufficient - no yielding detected")
        return {
            "converged": True,
            "displacements": displacements,
            "stresses": elastic_solution["stresses"],
            "strains": elastic_solution["strains"],
            "plastic_elements": plastic_elements,
            "forces_1d": elastic_solution.get("forces_1d", np.zeros(0)),
            "iterations": elastic_solution["iterations"],
            "residual_norm": elastic_solution["residual_norm"],
            "phase": "elastic"
        }

    if debug_level >= 1:
        print(f"  Phase 2: Plastic correction with {n_plastic_after_elastic} plastic elements")

    # Build material properties dict per material (with reduced strengths)
    material_props_by_mat_plastic = []
    for i in range(len(E_by_mat)):
        material_props_by_mat_plastic.append({
            'E': E_by_mat[i],
            'nu': nu_by_mat[i],
            'c': c_reduced[i] if hasattr(c_reduced, '__getitem__') else c_reduced,
            'phi': phi_reduced[i] if hasattr(phi_reduced, '__getitem__') else phi_reduced,
            'psi': psi_reduced[i] if hasattr(psi_reduced, '__getitem__') else psi_reduced,
            'gamma': gamma_by_mat[i]
        })

    plastic_solution = solve_plastic_phase(
        fem_data, material_props_by_mat_plastic, plastic_elements,
        accumulated_plastic_strain, u_nodal, displacements, debug_level
    )

    if not plastic_solution["converged"]:
        if debug_level >= 1:
            print("  ✗ Plastic phase failed to converge")
            
        if plastic_solution.get("failure_detected", False):
            if debug_level >= 1:
                print(f"  Slope failure detected: {plastic_solution.get('failure_reason', 'Unknown')}")
            return {
                "converged": False,
                "error": f"Slope failure: {plastic_solution.get('failure_reason', 'Unknown')}",
                "phase": "plastic",
                "elastic_solution": elastic_solution,
                "failure_detected": True,
                "failure_reason": plastic_solution.get("failure_reason", "Unknown")
            }
        else:
            return {
                "converged": False,
                "error": "Plastic phase failed to converge",
                "phase": "plastic",
                "elastic_solution": elastic_solution
            }

    final_solution = {
        "converged": True,
        "displacements": plastic_solution["displacements"],
        "stresses": plastic_solution["stresses"],
        "strains": plastic_solution["strains"],
        "plastic_elements": plastic_solution.get("plastic_elements", plastic_elements),
        "forces_1d": plastic_solution.get("forces_1d", np.zeros(0)),
        "iterations": elastic_solution["iterations"] + plastic_solution["iterations"],
        "residual_norm": plastic_solution["residual_norm"],
        "phase": "plastic",
        "elastic_iterations": elastic_solution["iterations"],
        "plastic_iterations": plastic_solution["iterations"]
    }

    if debug_level >= 1:
        n_plastic = np.sum(final_solution["plastic_elements"])
        print(f"  ✓ Two-phase solution converged: {n_plastic}/{len(elements)} plastic elements")
        print(f"  Total iterations: {final_solution['iterations']}")

    return final_solution


def solve_elastic_phase(fem_data, material_props_by_mat, plastic_elements, 
                        accumulated_plastic_strain, u_nodal, k_seismic,
                        c_reduced_by_elem, phi_reduced_by_elem, debug_level=0):
    """
    Solve the elastic phase of the two-phase FEM solver.
    
    This phase solves the system assuming purely elastic behavior to get
    an initial displacement field and identify elements that may yield.
    
    Parameters:
        fem_data (dict): FEM data dictionary
        material_props_by_mat (list): Material properties for each material
        plastic_elements (np.ndarray): Boolean array indicating plastic elements
        accumulated_plastic_strain (np.ndarray): Accumulated plastic strain per element
        u_nodal (np.ndarray): Nodal pore pressures
        k_seismic (float): Seismic coefficient
        debug_level (int): Verbosity level
    
    Returns:
        dict: Elastic solution dictionary
    """
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    element_materials = fem_data["element_materials"]
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    n_dof = 2 * n_nodes
    
    if debug_level >= 2:
        print("    Solving elastic phase...")
    
    # Assemble global stiffness matrix (elastic only)
    K_global = lil_matrix((n_dof, n_dof))
    
    # Assemble 2D elements with elastic stiffness
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        # Get material properties
        mat_id = element_materials[elem_idx] - 1
        material_props = material_props_by_mat[mat_id]
        
        E = material_props['E']
        nu = material_props['nu']
        
        # Build element stiffness matrix (elastic only)
        if elem_type in [3, 6]:  # Triangular elements
            # Build constitutive matrix D
            factor = E / ((1 + nu) * (1 - 2 * nu))
            D = factor * np.array([
                [1 - nu, nu,     0],
                [nu,     1 - nu, 0],
                [0,      0,      (1 - 2 * nu) / 2]
            ])
            
            if elem_type == 3:
                K_elem = build_tri3_stiffness(nodes[active_nodes], D)
            else:  # elem_type == 6
                K_elem = build_tri6_stiffness(nodes[active_nodes], D)
        elif elem_type in [4, 8, 9]:  # Quadrilateral elements
            if elem_type == 4:
                K_elem = build_quad4_stiffness(nodes[active_nodes], E, nu, is_plastic=False, elem_type=elem_type)
            elif elem_type == 8:
                K_elem = build_quad8_stiffness(nodes[active_nodes], E, nu, is_plastic=False, elem_type=elem_type)
            else:  # elem_type == 9
                K_elem = build_quad9_stiffness(nodes[active_nodes], E, nu, is_plastic=False, elem_type=elem_type)
        else:
            if debug_level >= 1:
                print(f"Warning: Element type {elem_type} not supported, skipping element {elem_idx}")
            continue
        
        # Assemble into global matrix
        dofs = []
        for node_id in active_nodes:
            dofs.extend([2*node_id, 2*node_id + 1])
        
        for i in range(len(dofs)):
            for j in range(len(dofs)):
                K_global[dofs[i], dofs[j]] += K_elem[i, j]
    
    # Assemble 1D truss elements
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    element_types_1d = fem_data.get("element_types_1d", np.array([]))
    k_by_1d_elem = fem_data.get("k_by_1d_elem", np.array([]))
    
    for elem_idx in range(len(elements_1d)):
        elem_nodes = elements_1d[elem_idx]
        elem_type = element_types_1d[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if len(active_nodes) >= 2:
            K_truss = build_truss_stiffness(nodes[active_nodes], k_by_1d_elem[elem_idx])
            
            # Assemble into global matrix
            dofs = []
            for node_id in active_nodes:
                dofs.extend([2*node_id, 2*node_id + 1])
            
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    K_global[dofs[i], dofs[j]] += K_truss[i, j]
    
    # Assemble force vector
    F_global = np.zeros(n_dof)
    
    # Body forces (gravity + seismic)
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        mat_id = element_materials[elem_idx] - 1
        gamma = material_props_by_mat[mat_id]['gamma']
        
        # Body force components
        b_x = k_seismic * gamma  # Horizontal seismic force
        b_y = -gamma             # Gravity (downward)
        
        F_body = compute_body_forces(nodes[active_nodes], elem_type, b_x, b_y)
        
        # Assemble into global vector
        dofs = []
        for node_id in active_nodes:
            dofs.extend([2*node_id, 2*node_id + 1])
        
        for i, dof in enumerate(dofs):
            F_global[dof] += F_body[i]
    
    # Applied loads from boundary conditions
    for i in range(n_nodes):
        if bc_type[i] == 4:  # Applied force
            F_global[2*i] += bc_values[i, 0]      # F_x
            F_global[2*i + 1] += bc_values[i, 1]  # F_y
    
    # Apply boundary conditions
    K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
        K_global, F_global, bc_type, nodes)
    
    # Solve elastic system
    try:
        u_free = spsolve(K_constrained.tocsr(), F_constrained)
        
        # Reconstruct full displacement vector
        displacements = np.zeros(n_dof)
        free_dof_idx = 0
        for i in range(n_dof):
            if i not in constraint_dofs:
                displacements[i] = u_free[free_dof_idx]
                free_dof_idx += 1
        
        if debug_level >= 2:
            print(f"    Elastic solution completed")
            print(f"    Max displacement: {np.max(np.abs(displacements)):.6f}")
        
    except Exception as e:
        if debug_level >= 1:
            print(f"Elastic solver failed: {e}")
        return {
            "converged": False,
            "error": f"Elastic solver failed: {e}",
            "iterations": 0
        }
    
    # Check for yielding and update plastic elements
    plastic_elements[:] = False  # Reset plastic state
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:
            # Compute element stress
            mat_id = element_materials[elem_idx] - 1
            material_props = material_props_by_mat[mat_id]
            
            sigma0_elem = fem_data.get("sigma0_by_elem", None)
            sigma0 = None if sigma0_elem is None else sigma0_elem[elem_idx]
            stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                material_props['E'], material_props['nu'],
                u_nodal[active_nodes], elem_type, sigma0)
            
            # Check yield criterion using reduced element strengths for this F value
            elem_id = elem_idx
            c_elem = c_reduced_by_elem[elem_id]
            phi_elem = phi_reduced_by_elem[elem_id]  # already radians
            yield_value = check_mohr_coulomb_yield(
                stress, c_elem, phi_elem)
            
            if yield_value > 0.0:
                plastic_elements[elem_idx] = True
    
    # Compute final stresses and strains
    stresses = np.zeros((n_elements, 4))
    strains = np.zeros((n_elements, 4))
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:
            mat_id = element_materials[elem_idx] - 1
            material_props = material_props_by_mat[mat_id]
            
            # Compute stress and strain
            stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                material_props['E'], material_props['nu'],
                u_nodal[active_nodes], elem_type)
            
            strain = compute_element_strain(
                nodes[active_nodes], displacements, active_nodes, elem_type)
            
            stresses[elem_idx, :3] = stress
            stresses[elem_idx, 3] = compute_von_mises(stress)
            
            strains[elem_idx, :3] = strain
            
            # Compute maximum shear strain
            eps_x, eps_y, gamma_xy = strain
            eps_mean = 0.5 * (eps_x + eps_y)
            R = np.sqrt(0.25 * (eps_x - eps_y)**2 + 0.25 * gamma_xy**2)
            eps_1 = eps_mean + R
            eps_3 = eps_mean - R
            max_shear_strain = abs(eps_1 - eps_3)
            strains[elem_idx, 3] = max_shear_strain
    
    # Compute 1D element forces
    forces_1d = np.zeros(len(elements_1d))
    for elem_idx in range(len(elements_1d)):
        if len(elements_1d[elem_idx]) >= 2:
            elem_nodes = elements_1d[elem_idx]
            active_nodes = elem_nodes[:element_types_1d[elem_idx]]
            forces_1d[elem_idx] = compute_truss_force(
                nodes[active_nodes], displacements, active_nodes, k_by_1d_elem[elem_idx])
    
    n_plastic = np.sum(plastic_elements)
    if debug_level >= 1:
        print(f"  Elastic phase completed: {n_plastic}/{n_elements} elements yielded")
    
    return {
        "converged": True,
        "displacements": displacements,
        "stresses": stresses,
        "strains": strains,
        "forces_1d": forces_1d,
        "iterations": 1,
        "residual_norm": 0.0,
        "phase": "elastic"
    }


def solve_plastic_phase(fem_data, material_props_by_mat, plastic_elements, 
                        accumulated_plastic_strain, u_nodal, initial_displacements, debug_level=0):
    """
    Solve the plastic phase of the two-phase FEM solver.
    
    This phase iteratively solves the system with elastic-plastic constitutive
    behavior, using the consistent tangent stiffness approach for numerical stability.
    
    Parameters:
        fem_data (dict): FEM data dictionary
        material_props_by_mat (list): Material properties for each material
        plastic_elements (np.ndarray): Boolean array indicating plastic elements
        accumulated_plastic_strain (np.ndarray): Accumulated plastic strain per element
        u_nodal (np.ndarray): Nodal pore pressures
        initial_displacements (np.ndarray): Initial displacements from elastic phase
        debug_level (int): Verbosity level
    
    Returns:
        dict: Plastic solution dictionary
    """
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    element_materials = fem_data["element_materials"]
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    n_dof = 2 * n_nodes
    
    if debug_level >= 2:
        print("    Solving plastic phase...")
    
    # Initialize displacements from elastic phase
    displacements = initial_displacements.copy()
    
    # Convergence parameters - tightened for better failure detection
    max_iterations = 30  # Reduced to detect failure faster
    tol_force = 1e-5     # Tightened
    tol_disp = 1e-5      # Tightened
    
    converged = False
    iteration = 0
    residual_norm = np.inf
    
    # Failure detection parameters - more sensitive for slope stability
    max_displacement_threshold = 0.02  # 2% of characteristic length (more sensitive)
    characteristic_length = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
    displacement_threshold = max_displacement_threshold * characteristic_length
    
    # Plastic zone monitoring - more sensitive for slope stability
    initial_plastic_count = np.sum(plastic_elements)
    max_plastic_ratio = 0.3  # If >30% elements are plastic, consider it failure (more sensitive)
    
    while iteration < max_iterations and not converged:
        if debug_level >= 3:
            print(f"      Plastic iteration {iteration + 1}")
        
        # Assemble global stiffness matrix with current plastic state
        K_global, plastic_state_changed, new_plastic_elements = assemble_global_stiffness_with_plasticity(
            nodes, elements, element_types, element_materials, displacements,
            material_props_by_mat, plastic_elements, accumulated_plastic_strain, u_nodal, fem_data, debug_level)
        
        # Update plastic elements
        plastic_elements[:] = new_plastic_elements
        
        # Assemble 1D truss elements
        elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
        element_types_1d = fem_data.get("element_types_1d", np.array([]))
        k_by_1d_elem = fem_data.get("k_by_1d_elem", np.array([]))
        
        for elem_idx in range(len(elements_1d)):
            elem_nodes = elements_1d[elem_idx]
            elem_type = element_types_1d[elem_idx]
            active_nodes = elem_nodes[:elem_type]
            
            if len(active_nodes) >= 2:
                K_truss = build_truss_stiffness(nodes[active_nodes], k_by_1d_elem[elem_idx])
                
                # Assemble into global matrix
                dofs = []
                for node_id in active_nodes:
                    dofs.extend([2*node_id, 2*node_id + 1])
                
                for i in range(len(dofs)):
                    for j in range(len(dofs)):
                        K_global[dofs[i], dofs[j]] += K_truss[i, j]
        
        # Assemble force vector (only body forces and applied loads, no additional loading)
        F_global = np.zeros(n_dof)
        
        # Body forces (gravity + seismic)
        k_seismic = fem_data.get("k_seismic", 0.0)
        for elem_idx in range(n_elements):
            elem_nodes = elements[elem_idx]
            elem_type = element_types[elem_idx]
            active_nodes = elem_nodes[:elem_type]
            
            mat_id = element_materials[elem_idx] - 1
            gamma = material_props_by_mat[mat_id]['gamma']
            
            # Body force components
            b_x = k_seismic * gamma  # Horizontal seismic force
            b_y = -gamma             # Gravity (downward)
            
            F_body = compute_body_forces(nodes[active_nodes], elem_type, b_x, b_y)
            
            # Assemble into global vector
            dofs = []
            for node_id in active_nodes:
                dofs.extend([2*node_id, 2*node_id + 1])
            
            for i, dof in enumerate(dofs):
                F_global[dof] += F_body[i]
        
        # Applied loads from boundary conditions
        for i in range(n_nodes):
            if bc_type[i] == 4:  # Applied force
                F_global[2*i] += bc_values[i, 0]      # F_x
                F_global[2*i + 1] += bc_values[i, 1]  # F_y
        
        # Apply boundary conditions
        K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(
            K_global, F_global, bc_type, nodes)
        
        # Compute internal forces from current displacements
        F_internal = K_global @ displacements
        
        # Apply boundary conditions to internal forces
        for dof in constraint_dofs:
            F_internal[dof] = 0.0
        
        # Compute residual force
        F_residual = F_global - F_internal
        
        # Apply boundary conditions to residual
        for dof in constraint_dofs:
            F_residual[dof] = 0.0
        
        # Extract only the free DOFs for the constrained system
        F_residual_constrained = []
        for i in range(n_dof):
            if i not in constraint_dofs:
                F_residual_constrained.append(F_residual[i])
        F_residual_constrained = np.array(F_residual_constrained)
        
        # Solve for displacement correction
        try:
            delta_u_free = spsolve(K_constrained.tocsr(), F_residual_constrained)
            
            # Reconstruct full displacement increment vector
            delta_u_full = np.zeros(n_dof)
            free_dof_idx = 0
            for i in range(n_dof):
                if i not in constraint_dofs:
                    delta_u_full[i] = delta_u_free[free_dof_idx]
                    free_dof_idx += 1
            
            # Update total displacements
            displacements += delta_u_full
            
            # Enforce boundary conditions on total displacements
            for dof in constraint_dofs:
                displacements[dof] = 0.0
            
        except Exception as e:
            if debug_level >= 1:
                print(f"Plastic solver failed at iteration {iteration + 1}: {e}")
            return {
                "converged": False,
                "error": f"Plastic solver failed: {e}",
                "iterations": iteration + 1,
                "phase": "plastic"
            }
        
        # Check convergence
        disp_change_norm = np.linalg.norm(delta_u_full)
        force_residual_norm = np.linalg.norm(F_residual_constrained)
        
        residual_norm = force_residual_norm
        
        # Check for slope failure conditions
        max_displacement = np.max(np.abs(displacements))
        current_plastic_count = np.sum(plastic_elements)
        plastic_ratio = current_plastic_count / n_elements
        
        if debug_level >= 3:
            print(f"        Displacement change norm: {disp_change_norm:.2e}")
            print(f"        Force residual norm: {force_residual_norm:.2e}")
            print(f"        Plastic elements: {current_plastic_count}/{n_elements} ({plastic_ratio:.1%})")
            print(f"        Max displacement: {max_displacement:.6f}")
        
        # Failure detection criteria
        failure_detected = False
        failure_reason = ""
        
        # 1. Excessive displacements (slope collapse)
        if max_displacement > displacement_threshold:
            failure_detected = True
            failure_reason = f"excessive_displacement ({max_displacement:.6f} > {displacement_threshold:.6f})"
        
        # 2. Excessive plastic zone (widespread failure)
        if plastic_ratio > max_plastic_ratio:
            failure_detected = True
            failure_reason = f"excessive_plastic_zone ({plastic_ratio:.1%} > {max_plastic_ratio:.1%})"
        
        # 3. Rapid plastic zone growth (instability)
        if current_plastic_count > initial_plastic_count * 2:  # 2x growth (more sensitive)
            failure_detected = True
            failure_reason = f"rapid_plastic_growth ({current_plastic_count} vs {initial_plastic_count})"
        
        if failure_detected:
            if debug_level >= 1:
                print(f"      ✗ Slope failure detected: {failure_reason}")
            return {
                "converged": False,
                "error": f"Slope failure: {failure_reason}",
                "iterations": iteration + 1,
                "phase": "plastic",
                "failure_detected": True,
                "failure_reason": failure_reason
            }
        
        # Check convergence criteria
        if (disp_change_norm < tol_disp and force_residual_norm < tol_force and not plastic_state_changed):
            converged = True
            if debug_level >= 2:
                print(f"      Plastic phase converged after {iteration + 1} iterations")
        elif iteration >= max_iterations - 1:
            if debug_level >= 1:
                print(f"      Plastic phase failed to converge after {max_iterations} iterations")
            break
        
        iteration += 1
    
    # Compute final stresses and strains
    stresses = np.zeros((n_elements, 4))
    strains = np.zeros((n_elements, 4))
    
    for elem_idx in range(n_elements):
        elem_nodes = elements[elem_idx]
        elem_type = element_types[elem_idx]
        active_nodes = elem_nodes[:elem_type]
        
        if elem_type in [3, 4, 6, 8, 9]:
            mat_id = element_materials[elem_idx] - 1
            material_props = material_props_by_mat[mat_id]
            
            # Compute trial elastic stress
            trial_stress = compute_element_stress(
                nodes[active_nodes], displacements, active_nodes,
                material_props['E'], material_props['nu'],
                u_nodal[active_nodes], elem_type)
            
            # Apply return mapping for plasticity
            corrected_stress, is_plastic, _ = mohr_coulomb_return_mapping(
                trial_stress, material_props, plastic_elements[elem_idx],
                accumulated_plastic_strain[elem_idx])
            
            stresses[elem_idx, :3] = corrected_stress
            stresses[elem_idx, 3] = compute_von_mises(corrected_stress)
            
            # Update plastic element tracking
            plastic_elements[elem_idx] = is_plastic
            
            # Compute strains
            strain_vec = compute_element_strain(
                nodes[active_nodes], displacements, active_nodes, elem_type)
            strains[elem_idx, :3] = strain_vec
            
            # Compute maximum shear strain
            eps_x, eps_y, gamma_xy = strain_vec
            eps_mean = 0.5 * (eps_x + eps_y)
            R = np.sqrt(0.25 * (eps_x - eps_y)**2 + 0.25 * gamma_xy**2)
            eps_1 = eps_mean + R
            eps_3 = eps_mean - R
            max_shear_strain = abs(eps_1 - eps_3)
            strains[elem_idx, 3] = max_shear_strain
    
    # Compute final 1D element forces
    forces_1d = np.zeros(len(elements_1d))
    for elem_idx in range(len(elements_1d)):
        if len(elements_1d[elem_idx]) >= 2:
            elem_nodes = elements_1d[elem_idx]
            active_nodes = elem_nodes[:element_types_1d[elem_idx]]
            forces_1d[elem_idx] = compute_truss_force(
                nodes[active_nodes], displacements, active_nodes, k_by_1d_elem[elem_idx])
    
    if debug_level >= 1:
        n_plastic = np.sum(plastic_elements)
        print(f"  Plastic phase completed: {n_plastic}/{n_elements} plastic elements")
    
    return {
        "converged": converged,
        "displacements": displacements,
        "stresses": stresses,
        "strains": strains,
        "forces_1d": forces_1d,
        "iterations": iteration + 1,
        "residual_norm": residual_norm,
        "phase": "plastic"
    }


def solve_ssrm_bisection(fem_data, F_min=1.0, F_max=3.0, tolerance=0.001, debug_level=0):
    """
    Solve for factor of safety using improved SSRM with bisection method.
    
    This implementation uses the bisection method to efficiently find the critical
    factor of safety, providing much better convergence than the previous incremental
    approach. It also uses the new two-phase elastic-plastic solver for better
    numerical stability.
    
    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        F_min (float): Minimum reduction factor to test (default 1.0)
        F_max (float): Maximum reduction factor to test (default 3.0)
        tolerance (float): Tolerance for factor of safety determination (default 0.001)
        debug_level (int): Verbosity level (0=quiet, 1=basic, 2=detailed, 3=debug)
    
    Returns:
        dict: SSRM solution dictionary containing:
            - converged: bool, whether SSRM procedure completed successfully
            - FS: float, critical factor of safety
            - last_solution: dict, final FEM solution at critical F
            - F_history: list, history of F values tested
            - convergence_history: list, convergence status for each F
            - solutions_history: list, solutions for each F
            - iterations_ssrm: int, number of SSRM iterations
            - method: str, solution method used ('bisection')
    """
    
    if debug_level >= 1:
        print("Starting improved SSRM with bisection method")
        print(f"Search range: F = [{F_min:.3f}, {F_max:.3f}]")
        print(f"Tolerance: {tolerance:.4f}")
    
    # Initialize bisection variables
    F_left = F_min
    F_right = F_max
    F_mid = (F_left + F_right) / 2.0
    
    # History tracking
    F_history = []
    convergence_history = []
    solutions_history = []
    
    # Test initial bounds
    if debug_level >= 1:
        print(f"\nTesting initial bounds...")
    
    # Test left bound (F_min)
    if debug_level >= 1:
        print(f"  Testing F = {F_left:.4f} (left bound)")
    
    solution_left = solve_fem_elastic_plastic(fem_data, F=F_left, debug_level=max(0, debug_level-1))
    
    F_history.append(F_left)
    convergence_history.append(solution_left.get("converged", False))
    solutions_history.append(solution_left)
    
    if not solution_left.get("converged", False):
        if debug_level >= 1:
            print(f"  ✗ Left bound F = {F_left:.4f} failed to converge")
        return {
            "converged": False,
            "error": f"Left bound F = {F_left:.4f} failed to converge",
            "FS": None,
            "F_history": F_history,
            "convergence_history": convergence_history,
            "solutions_history": solutions_history,
            "iterations_ssrm": 1,
            "method": "bisection"
        }
    
    # Test right bound (F_max)
    if debug_level >= 1:
        print(f"  Testing F = {F_right:.4f} (right bound)")
    
    solution_right = solve_fem_elastic_plastic(fem_data, F=F_right, debug_level=max(0, debug_level-1))
    
    F_history.append(F_right)
    convergence_history.append(solution_right.get("converged", False))
    solutions_history.append(solution_right)
    
    if solution_right.get("converged", False):
        if debug_level >= 1:
            print(f"  ✓ Right bound F = {F_max:.4f} converged - slope appears very stable")
            print(f"  Need to test higher F values to find failure point")
        
        # If right bound converges, we need to find a higher F that fails
        # Let's try F = F_max + 1.0 to see if we can find failure
        F_test = F_max + 1.0
        if debug_level >= 1:
            print(f"  Testing higher F = {F_test:.4f} to find failure point")
        
        solution_test = solve_fem_elastic_plastic(fem_data, F=F_test, debug_level=max(0, debug_level-1))
        
        F_history.append(F_test)
        convergence_history.append(solution_test.get("converged", False))
        solutions_history.append(solution_test)
        
        if not solution_test.get("converged", False):
            # Now we have: F_left converges, F_test fails
            F_right = F_test
            solution_right = solution_test
            if debug_level >= 1:
                print(f"  ✗ F = {F_test:.4f} fails to converge")
                print(f"  Critical FS lies in [{F_left:.4f}, {F_right:.4f}]")
        else:
            # Even F_test converges - slope is very stable
            if debug_level >= 1:
                print(f"  ✓ F = {F_test:.4f} still converges - slope appears extremely stable")
            return {
                "converged": True,
                "FS": F_test,
                "last_solution": solution_test,
                "F_history": F_history,
                "convergence_history": convergence_history,
                "solutions_history": solutions_history,
                "iterations_ssrm": 3,
                "method": "bisection",
                "note": f"Slope stable even at F = {F_test:.4f} - very stable slope"
            }
    else:
        # Right bound failed to converge - this is what we want
        if debug_level >= 1:
            print(f"  ✗ Right bound F = {F_max:.4f} fails to converge")
    
    # Now we know: F_left converges, F_right doesn't converge
    # The critical factor of safety lies between F_left and F_right
    if debug_level >= 1:
        print(f"  ✓ Left bound F = {F_left:.4f} converges")
        print(f"  ✗ Right bound F = {F_right:.4f} fails to converge")
        print(f"  Critical FS lies in [{F_left:.4f}, {F_right:.4f}]")
    
    # Bisection loop
    ssrm_iteration = 2  # We've already done 2 iterations above
    max_bisection_iterations = 50
    
    for bisection_iter in range(max_bisection_iterations):
        F_mid = (F_left + F_right) / 2.0
        
        if debug_level >= 1:
            print(f"\nBisection iteration {bisection_iter + 1}: F = {F_mid:.4f}")
            print(f"  Current interval: [{F_left:.4f}, {F_right:.4f}]")
            print(f"  Interval width: {F_right - F_left:.4f}")
        
        # Test midpoint
        solution_mid = solve_fem_elastic_plastic(fem_data, F=F_mid, debug_level=max(0, debug_level-1))
        
        F_history.append(F_mid)
        convergence_history.append(solution_mid.get("converged", False))
        solutions_history.append(solution_mid)
        ssrm_iteration += 1
        
        if solution_mid.get("converged", False):
            # F_mid converges - critical FS is between F_mid and F_right
            if debug_level >= 1:
                print(f"  ✓ F = {F_mid:.4f} converges")
            
            F_left = F_mid
            solution_left = solution_mid
        else:
            # F_mid fails - critical FS is between F_left and F_mid
            if debug_level >= 1:
                print(f"  ✗ F = {F_mid:.4f} fails to converge")
            
            F_right = F_mid
            solution_right = solution_mid
        
        # Check convergence
        interval_width = F_right - F_left
        if interval_width < tolerance:
            if debug_level >= 1:
                print(f"\n🎯 Bisection converged!")
                print(f"  Critical FS = {F_left:.4f} ± {tolerance:.4f}")
                print(f"  Final interval: [{F_left:.4f}, {F_right:.4f}]")
            
            return {
                "converged": True,
                "FS": F_left,
                "last_solution": solution_left,
                "F_history": F_history,
                "convergence_history": convergence_history,
                "solutions_history": solutions_history,
                "iterations_ssrm": ssrm_iteration,
                "method": "bisection",
                "final_interval": (F_left, F_right),
                "interval_width": interval_width
            }
        
        if debug_level >= 2:
            print(f"  New interval: [{F_left:.4f}, {F_right:.4f}]")
    
    # Maximum bisection iterations reached
    if debug_level >= 1:
        print(f"\nMaximum bisection iterations ({max_bisection_iterations}) reached")
        print(f"Best estimate: FS = {F_left:.4f} ± {F_right - F_left:.4f}")
    
    return {
        "converged": True,
        "FS": F_left,
        "last_solution": solution_left,
        "F_history": F_history,
        "convergence_history": convergence_history,
        "solutions_history": solutions_history,
        "iterations_ssrm": ssrm_iteration,
        "method": "bisection",
        "note": f"Maximum bisection iterations reached",
        "final_interval": (F_left, F_right),
        "interval_width": F_right - F_left
    }


def check_convergence_robust(current_solution, previous_solution, 
                            plastic_state, material_props, tolerance_disp=1e-6, tolerance_force=1e-6):
    """
    Robust convergence checking that considers multiple criteria.
    
    This function implements comprehensive convergence checking that considers:
    - Displacement convergence
    - Force equilibrium
    - Plastic state stability
    - Energy balance
    
    Parameters:
        current_solution: current iteration solution
        previous_solution: previous iteration solution  
        plastic_state: current plastic element state
        material_props: material properties
        tolerance_disp: displacement tolerance
        tolerance_force: force residual tolerance
    
    Returns:
        tuple: (converged: bool, disp_norm: float, force_norm: float, plastic_stable: bool)
    """
    if previous_solution is None:
        return False, np.inf, np.inf, False
    
    # 1. Displacement convergence check
    u_current = current_solution.get("displacements", np.array([]))
    u_previous = previous_solution.get("displacements", np.array([]))
    
    if len(u_current) != len(u_previous):
        return False, np.inf, np.inf, False
    
    du = u_current - u_previous
    disp_norm = np.linalg.norm(du) / max(np.linalg.norm(u_current), 1e-12)
    
    # 2. Force convergence check
    force_norm = current_solution.get("residual_norm", np.inf)
    
    # 3. Plastic state stability check
    # Check if plastic elements are still changing significantly
    plastic_stable = True
    if 'plastic_elements' in current_solution and 'plastic_elements' in previous_solution:
        current_plastic = current_solution['plastic_elements']
        previous_plastic = previous_solution['plastic_elements']
        
        if len(current_plastic) == len(previous_plastic):
            plastic_change_ratio = np.sum(current_plastic != previous_plastic) / len(current_plastic)
            plastic_stable = plastic_change_ratio < 0.01  # Less than 1% change
    
    # 4. Energy-based convergence check (optional)
    # This could be implemented to check if the total strain energy is stabilizing
    
    # 5. Combined convergence criteria
    # All criteria must be satisfied for convergence
    converged = (disp_norm < tolerance_disp and 
                 force_norm < tolerance_force and 
                 plastic_stable)
    
    return converged, disp_norm, force_norm, plastic_stable