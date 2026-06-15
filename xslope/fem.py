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
import warnings
from math import degrees, sin, cos, sqrt, asin, tan


import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import splu
from shapely.geometry import LineString, Point


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


def export_fem_solution(fem_data, solution, output_stem):
    """Export FEM nodal and element results to CSV files using a common stem."""
    import pandas as pd
    from pathlib import Path

    output_stem = Path(output_stem)
    nodes_file = output_stem.parent / f"{output_stem.name}_fem_nodes.csv"
    elements_file = output_stem.parent / f"{output_stem.name}_fem_elements.csv"

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

    with open(nodes_file, "w") as f:
        node_df.to_csv(f, index=False)

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

    with open(elements_file, "w") as f:
        element_df.to_csv(f, index=False)

    print(f"Exported FEM nodal results to {nodes_file}")
    print(f"Exported FEM element results to {elements_file}")


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
            - cos_theta_1d: np.ndarray (n_1d_elements,) of direction cosines (x) for each 1D element
            - sin_theta_1d: np.ndarray (n_1d_elements,) of direction cosines (y) for each 1D element
            - dof_indices_1d: np.ndarray (n_1d_elements, 4) of global DOF indices using dof_offset
            - K_global_1d_elems: list of np.ndarray (4, 4) global stiffness matrices for each 1D element
            - dof_offset: np.ndarray (n_nodes+1,) cumulative DOF count; pile nodes get 3 DOFs, others get 2
            - is_pile_node: np.ndarray (n_nodes,) boolean, True for nodes belonging to pile elements
            - n_dof_total: int, total number of DOFs (dof_offset[n_nodes])
            - dof_indices_pile: np.ndarray (n_pile_elements, 6) of global DOF indices for 6-DOF beam elements
            - K_global_pile_elems: list of np.ndarray (6, 6) global stiffness matrices for pile beam elements
            - EI_by_pile_elem: np.ndarray (n_pile_elements,) of flexural rigidity per unit width
            - EA_by_pile_elem: np.ndarray (n_pile_elements,) of axial rigidity per unit width
            - pile_head_nodes: np.ndarray of node indices for each pile line's top node
            - pile_head_fixed: np.ndarray of booleans for fixity of each pile head
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
        else:
            c_by_elem[elem_idx] = c_by_mat[mat_id]
            phi_by_elem[elem_idx] = phi_by_mat[mat_id]
    
    # Process pore pressures
    u = np.zeros(n_nodes)
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
        if "seep_u" in slope_data:
            seep_u = slope_data["seep_u"]
            if isinstance(seep_u, np.ndarray) and len(seep_u) == n_nodes:
                u = np.maximum(0.0, seep_u)
            else:
                n_seep = len(seep_u) if hasattr(seep_u, "__len__") else "?"
                warnings.warn(
                    f"Seepage pore pressures NOT applied: the stored seep solution has {n_seep} "
                    f"values but the FEM mesh has {n_nodes} nodes. Pore pressures default to 0, "
                    "which over-predicts the factor of safety — run the FEM on the same mesh as "
                    "the seepage solution.")
    
    # Process 1D reinforcement elements
    elements_1d = np.array([]).reshape(0, 3) if 'elements_1d' not in mesh else mesh['elements_1d']
    element_types_1d = np.array([]) if 'element_types_1d' not in mesh else mesh['element_types_1d'] 
    element_materials_1d = np.array([]) if 'element_materials_1d' not in mesh else mesh['element_materials_1d']
    
    n_1d_elements = len(elements_1d)
    
    t_allow_by_1d_elem = np.zeros(n_1d_elements)
    t_res_by_1d_elem = np.zeros(n_1d_elements)
    k_by_1d_elem = np.zeros(n_1d_elements)
    cos_theta_1d = np.zeros(n_1d_elements)
    sin_theta_1d = np.zeros(n_1d_elements)
    dof_indices_1d = np.zeros((n_1d_elements, 4), dtype=int)
    K_global_1d_elems = []

    if n_1d_elements > 0 and "reinforcement_lines" in slope_data:
        reinforcement_lines = slope_data["reinforcement_lines"]

        for elem_idx in range(n_1d_elements):
            line_id = element_materials_1d[elem_idx] - 1  # Convert to 0-based

            if line_id < len(reinforcement_lines):
                line_data = reinforcement_lines[line_id]

                # Get element geometry — use only end nodes [0] and [1],
                # ignore mid-node for quadratic elements
                elem_nodes_1d = elements_1d[elem_idx]
                node_0 = elem_nodes_1d[0]
                node_1 = elem_nodes_1d[1]
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

                    # DOF indices for end nodes (4 DOFs total)
                    dof_indices_1d[elem_idx] = [2*node_0, 2*node_0+1, 2*node_1, 2*node_1+1]

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
                    k_val = E * A / elem_length
                    k_by_1d_elem[elem_idx] = k_val

                    # Build 4x4 truss element stiffness matrix in global coordinates
                    # T = [[cos, sin, 0, 0], [0, 0, cos, sin]]  (2x4)
                    # K_local = k * [[1, -1], [-1, 1]]           (2x2)
                    # K_global_elem = T^T @ K_local @ T          (4x4)
                    c2 = cos_t * cos_t
                    cs = cos_t * sin_t
                    s2 = sin_t * sin_t
                    K_global_elem = k_val * np.array([
                        [ c2,  cs, -c2, -cs],
                        [ cs,  s2, -cs, -s2],
                        [-c2, -cs,  c2,  cs],
                        [-cs, -s2,  cs,  s2]
                    ])
                    K_global_1d_elems.append(K_global_elem)
                else:
                    K_global_1d_elems.append(np.zeros((4, 4)))
            else:
                K_global_1d_elems.append(np.zeros((4, 4)))

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
    pile_elem_indices = []  # maps pile element index to global 1D element index
    V_cap_by_pile_elem = []
    M_cap_by_pile_elem = []
    elem_length_by_pile_elem = []
    S_by_pile_elem = []
    EI_by_pile_elem = []
    EA_by_pile_elem = []
    pile_node_pairs = []  # (node_0, node_1) for each pile element
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

                # Build full 6x6 Euler-Bernoulli beam stiffness in local coords
                # DOFs: [u1, v1, theta1, u2, v2, theta2]
                # u = axial, v = transverse, theta = rotation
                L2 = L * L
                L3 = L2 * L
                K_local = np.array([
                    [ EA/L,   0.0,          0.0,        -EA/L,  0.0,          0.0       ],
                    [ 0.0,    12*EI/L3,     6*EI/L2,     0.0,  -12*EI/L3,    6*EI/L2   ],
                    [ 0.0,    6*EI/L2,      4*EI/L,      0.0,  -6*EI/L2,     2*EI/L    ],
                    [-EA/L,   0.0,          0.0,         EA/L,   0.0,         0.0       ],
                    [ 0.0,   -12*EI/L3,    -6*EI/L2,     0.0,   12*EI/L3,   -6*EI/L2   ],
                    [ 0.0,    6*EI/L2,      2*EI/L,      0.0,  -6*EI/L2,     4*EI/L    ],
                ])

                # 6x6 rotation matrix T (local -> global)
                # local x along element, local y perpendicular
                c = cos_t
                s = sin_t
                T = np.array([
                    [ c,  s, 0, 0, 0, 0],
                    [-s,  c, 0, 0, 0, 0],
                    [ 0,  0, 1, 0, 0, 0],
                    [ 0,  0, 0, c, s, 0],
                    [ 0,  0, 0,-s, c, 0],
                    [ 0,  0, 0, 0, 0, 1],
                ])
                K_beam = T.T @ K_local @ T

                cos_theta_pile.append(cos_t)
                sin_theta_pile.append(sin_t)
                pile_node_pairs.append((node_0, node_1))
                K_global_pile_elems.append(K_beam)
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
    is_pile_node = np.zeros(n_nodes, dtype=bool)
    for p_idx in range(n_pile_elements):
        n0, n1 = pile_node_pairs[p_idx]
        is_pile_node[n0] = True
        is_pile_node[n1] = True

    dof_offset = np.zeros(n_nodes + 1, dtype=int)
    for i in range(n_nodes):
        dof_offset[i + 1] = dof_offset[i] + (3 if is_pile_node[i] else 2)
    n_dof_total = int(dof_offset[n_nodes])

    # Build 6-element DOF indices for pile elements (using dof_offset)
    dof_indices_pile = np.zeros((n_pile_elements, 6), dtype=int) if n_pile_elements > 0 else np.zeros((0, 6), dtype=int)
    for p_idx in range(n_pile_elements):
        n0, n1 = pile_node_pairs[p_idx]
        dof_indices_pile[p_idx] = [
            dof_offset[n0], dof_offset[n0] + 1, dof_offset[n0] + 2,
            dof_offset[n1], dof_offset[n1] + 1, dof_offset[n1] + 2,
        ]

    # Rebuild 1D truss DOF indices using dof_offset (in case any share nodes with piles)
    for elem_idx in range(n_1d_elements):
        if pile_elem_mask[elem_idx]:
            continue
        elem_nodes_1d = elements_1d[elem_idx]
        node_0 = elem_nodes_1d[0]
        node_1 = elem_nodes_1d[1]
        dof_indices_1d[elem_idx] = [dof_offset[node_0], dof_offset[node_0] + 1,
                                     dof_offset[node_1], dof_offset[node_1] + 1]

    # Identify pile head nodes and their fixity for boundary conditions
    # The pile head is the top node (highest y) of each pile line.
    pile_head_nodes = []
    pile_head_fixed = []
    for pl_idx in range(n_pile_lines):
        pile_data = pile_lines[pl_idx]
        fixity = pile_data.get("fixity", "free")

        # Collect all nodes belonging to this pile line
        pile_nodes_for_line = set()
        for p_idx in range(n_pile_elements):
            if pile_line_idx_by_pile_elem[p_idx] == pl_idx:
                n0, n1 = pile_node_pairs[p_idx]
                pile_nodes_for_line.add(n0)
                pile_nodes_for_line.add(n1)

        if pile_nodes_for_line:
            # Top node = highest y coordinate
            top_node = max(pile_nodes_for_line, key=lambda nd: nodes[nd, 1])
            pile_head_nodes.append(top_node)
            pile_head_fixed.append(fixity == "fixed")

    pile_head_nodes = np.array(pile_head_nodes, dtype=int)
    pile_head_fixed = np.array(pile_head_fixed, dtype=bool)

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
    
    # Save displacement constraints before force BCs can overwrite them.
    # Nodes on boundary faces that also receive distributed loads need both
    # their displacement constraint (roller/fixed) AND the applied force.
    fixed_nodes = set(np.where(bc_type == 1)[0])
    roller_x_nodes = set(np.where(bc_type == 2)[0])
    roller_y_nodes = set(np.where(bc_type == 3)[0])

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
        tolerance = 1e-1  # Tolerance for finding nodes on load lines

        for load_idx, load_line in enumerate(distributed_loads):
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
            _seen_edges = set()
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
                            continue
                        _seen_edges.add(_key)
                        _load_edges.append(_gn)

            if _load_edges:
                nodal_fx = {}
                nodal_fy = {}
                _g = 1.0 / np.sqrt(3.0)
                for _gn in _load_edges:
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
                print(f"  Distributed load {load_idx}: {len(_load_edges)} edges / "
                      f"{len(load_nodes)} nodes (consistent), "
                      f"total force = {total_force:.1f}, expected ~{expected_force:.1f}")
                continue

            # Pass 2b (fallback when no boundary edges found on the line):
            # tributary lumping per node.
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
            print(f"  Distributed load {load_idx}: {n_load_nodes} nodes, "
                  f"total force = {total_force:.1f}, expected ~{expected_force:.1f}, "
                  f"sum(trib) = {total_trib:.2f}")
            
    
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
        "fixed_nodes": fixed_nodes,
        "roller_x_nodes": roller_x_nodes,
        "roller_y_nodes": roller_y_nodes,
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
        "cos_theta_1d": cos_theta_1d,
        "sin_theta_1d": sin_theta_1d,
        "dof_indices_1d": dof_indices_1d,
        "K_global_1d_elems": K_global_1d_elems,
        "unit_weight": unit_weight,
        "k_seismic": k_seismic,
        "pp_option": pp_option,
        "piezo_line_coords": piezo_line_coords,
        "gamma_water": slope_data.get("gamma_water", 9.81),
        # DOF offset map (pile nodes get 3 DOFs, others get 2)
        "dof_offset": dof_offset,
        "is_pile_node": is_pile_node,
        "n_dof_total": n_dof_total,
        # Pile beam elements (6-DOF Euler-Bernoulli)
        "n_pile_elements": n_pile_elements,
        "pile_elem_mask": pile_elem_mask,
        "pile_elem_indices": pile_elem_indices,
        "cos_theta_pile": cos_theta_pile,
        "sin_theta_pile": sin_theta_pile,
        "dof_indices_pile": dof_indices_pile,
        "K_global_pile_elems": K_global_pile_elems,
        "V_cap_by_pile_elem": V_cap_by_pile_elem,
        "M_cap_by_pile_elem": M_cap_by_pile_elem,
        "elem_length_by_pile_elem": elem_length_by_pile_elem,
        "S_by_pile_elem": S_by_pile_elem,
        "EI_by_pile_elem": EI_by_pile_elem,
        "EA_by_pile_elem": EA_by_pile_elem,
        "pile_node_pairs": pile_node_pairs,
        "pile_line_idx_by_pile_elem": pile_line_idx_by_pile_elem,
        "pile_head_nodes": pile_head_nodes,
        "pile_head_fixed": pile_head_fixed,
    }


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

def solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=3000, tolerance=1e-3,
              max_disp_factor=0.1, tension_cutoff=False, staged=False, dt_scale=1.0,
              pp_formulation='effective',
              early_exit=True):
    """
    Solve FEM using Griffiths & Lane (1999) viscoplastic algorithm.

    Implements the exact algorithm from the 1999 Geotechnique paper:
    - 8-node quadrilateral elements with reduced integration (4 Gauss points)
    - Viscoplastic stress redistribution with accumulated plastic strains
    - Non-convergence failure criterion with displacement limit safeguard
    - Pre-factored elastic stiffness matrix for efficiency
    - No damping (stability from dt parameter)
    - Direct solve each iteration (not residual-based)

    Parameters:
        fem_data (dict): FEM data dictionary from build_fem_data
        F (float): Shear strength reduction factor (c/F, tan(phi)/F)
        debug_level (int): 0=silent, 1=summary, 2=per-iteration
        max_iterations (int): Maximum viscoplastic iterations (default 3000)
        tolerance (float): Convergence tolerance ||du|| / ||u|| (default 1e-3).
            Normalized by the current displacement (Smith & Griffiths CHECON-style),
            so steady benign viscoplastic creep is accepted as converged; false
            convergence at true failure is guarded by max_disp_factor below.
        max_disp_factor (float): Displacement limit as fraction of mesh height (default 0.1).
            If max VP displacement (total - elastic) exceeds this fraction of the mesh height,
            the slope is declared as failed regardless of convergence. This prevents false
            convergence when large displacements make the relative change appear small.
            Set to None to disable.
        pp_formulation (str): How pore pressures enter the analysis.
            'effective' (default): u is moved into the load vector
            (K du = F_ext + int B^T m u dV) and the elastic stresses are
            EFFECTIVE stresses directly. This is the standard effective-stress
            formulation (Potts & Zdravkovic); under a flooded boundary it
            yields sigma'_v = gamma'*z with compressive lateral stresses.
            'total': legacy recipe — solve the total-stress elastic problem
            and subtract u pointwise at Gauss points before the yield check.
            With nu < 0.5 this produces a spurious effective-tension zone of
            magnitude ~(1-2*nu)/(1-nu)*u at submerged boundaries (the lateral
            elastic response to the water load is nu/(1-nu) of it while u
            subtracts all of it), which yields and creeps indefinitely.
        tension_cutoff (bool): If True (default False), states with effective mean
            TENSION are relaxed volumetrically as a second viscoplastic mechanism
            (damped, ~10% per iteration). The psi=0 Mohr-Coulomb flow is purely
            deviatoric and cannot return near-apex tensile states to the envelope.
            Equivalent in spirit to the tension cutoff used by commercial SSRM
            codes; Griffiths & Lane (1999) include no tension treatment. Rarely
            needed now that dt uses the stable p61 value (the apparent tension
            stall was a dt-induced limit cycle), but retained as an option.
        staged (bool): If True and the model has water loads (applied boundary
            forces and/or pore pressures), solve in two stages: stage 1 applies
            gravity only (dry), stage 2 adds the water loads and pore pressures,
            continuing from the converged stage-1 state with accumulated
            viscoplastic strains (construction history: built, then filled).
        dt_scale (float): Multiplier on the viscoplastic pseudo-timestep
            (research/diagnostic knob; default 1.0).

    Returns:
        dict: Solution dictionary with keys:
            - converged (bool): Whether iterations converged
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

    # Material properties
    c_by_elem = fem_data.get("c_by_elem", fem_data["c_by_mat"][element_materials - 1])
    phi_by_elem = fem_data.get("phi_by_elem", fem_data["phi_by_mat"][element_materials - 1])
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

    # Apply strength reduction (Griffiths & Lane 1999): c_r = c/F, phi_r = atan(tan(phi)/F)
    # Note: Only soil strength (c, phi) is reduced by F. Reinforcement properties
    # (T_allow, T_res, EA/L) are NOT reduced — they are structural capacities.
    c_reduced = c_by_elem / F
    tan_phi_reduced = np.tan(np.radians(phi_by_elem)) / F
    phi_reduced = np.arctan(tan_phi_reduced)  # radians

    if debug_level >= 1:
        print(f"  c: {c_by_elem[0]:.1f} -> {c_reduced[0]:.1f}")
        print(f"  phi: {phi_by_elem[0]:.1f} -> {np.degrees(phi_reduced[0]):.1f}")

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

        # Tracking arrays for 1D element status
        forces_1d = np.zeros(n_1d_elements)
        failed_1d = np.zeros(n_1d_elements, dtype=bool)

        if debug_level >= 1:
            print(f"  1D truss elements: {n_1d_elements}")

    # Extract pile beam element data (6-DOF Euler-Bernoulli)
    n_pile_elements = fem_data.get("n_pile_elements", 0)
    has_pile_elements = n_pile_elements > 0
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(n_1d_elements, dtype=bool))

    if has_pile_elements:
        cos_theta_pile = fem_data["cos_theta_pile"]
        sin_theta_pile = fem_data["sin_theta_pile"]
        dof_indices_pile = fem_data["dof_indices_pile"]
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

        if debug_level >= 1:
            print(f"  Pile beam elements: {n_pile_elements} (6-DOF Euler-Bernoulli)")

    # ---- Step 1: Build K_global (elastic, constant — includes 2D soil + 1D truss + pile beam) ----
    K_global = build_global_stiffness(nodes, elements, element_types,
                                      element_materials, E_by_mat, nu_by_mat,
                                      fem_data=fem_data)

    # ---- Step 2: Build gravity load vector ----
    F_gravity = build_gravity_loads(nodes, elements, element_types,
                                    element_materials, gamma_by_mat, k_seismic,
                                    fem_data=fem_data)

    # Gravity-only copy for staged loading (water/applied loads enter at stage 2)
    F_grav_pure = F_gravity.copy()

    # Add boundary condition forces (using dof_offset)
    for i in range(n_nodes):
        if bc_type[i] == 4:  # Force boundary condition
            dof_x = dof_offset[i] if dof_offset is not None else 2 * i
            F_gravity[dof_x] += bc_values[i, 0]
            F_gravity[dof_x + 1] += bc_values[i, 1]

    # ---- Step 3: Identify free/constrained DOFs ONCE ----
    # Use saved constraint sets so that nodes with both a displacement constraint
    # (roller/fixed) and an applied force (bc_type overwritten to 4) are still constrained.
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

    # Add rotation constraints for fixed pile heads
    if has_pile_elements:
        for ph_idx in range(len(pile_head_nodes)):
            if pile_head_fixed[ph_idx]:
                ph_node = pile_head_nodes[ph_idx]
                rot_dof = dof_offset[ph_node] + 2 if dof_offset is not None else None
                if rot_dof is not None:
                    constraint_dofs.append(rot_dof)

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


    if debug_level >= 1:
        print(f"  DOFs: {n_dof} total, {n_free} free, {len(constraint_dofs)} constrained")
        print(f"  K factorized (reused for all iterations)")

    # ---- Step 5: Compute dt from material properties (Smith & Griffiths) ----
    # Viscoplastic pseudo-timestep dt = 4*(1+nu)/(3*E) (S&G p61 form). The
    # Mohr-Coulomb variant 4*(1+nu)*(1-2nu)/(E*(1-2nu+sin^2(phi_r))) is ~2.6x
    # larger and was found to drive a limit cycle at Gauss points in mild
    # effective tension under reservoir loading (the redistribution overshoots
    # and never settles); the p61 dt is in the stable regime. Note that the
    # per-iteration displacement increment scales with dt, so the convergence
    # tolerance and the SSRM failure criterion are calibrated to this dt (see
    # docs/fem/overview.md).
    dt = 1.0e15
    dt_t = np.zeros(n_elements)   # tension-cutoff (volumetric) timestep per element
    for elem_idx in range(n_elements):
        mat_id = element_materials[elem_idx] - 1
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        ddt = 4.0 * (1.0 + nu) / (3.0 * E)
        if ddt < dt:
            dt = ddt
        # volumetric stiffness K = E/(3(1-2nu)). Damped relaxation: relax ~10%
        # of the mean tension per iteration (dt_t*K = 0.1). A full single-step
        # return (dt_t*K = 1) pumps against the elastic re-solve at water-
        # pressure boundaries and never settles; gentle relaxation lets the
        # surrounding field re-equilibrate and the zone converge.
        dt_t[elem_idx] = 0.3 * (1.0 - 2.0 * nu) / E

    dt *= dt_scale

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
        D4 = build_constitutive_matrix_4(E, nu)   # 4-comp (with sigma_z) for the VP loop

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

    # ---- Step 6b: Precompute pore pressure at each Gauss point ----
    u_nodes = fem_data.get("u", np.zeros(n_nodes))
    pp_option = fem_data.get("pp_option", "none")

    u_gp = []
    if pp_option == "none":
        for elem_idx in range(n_elements):
            u_gp.append([0.0] * len(elem_gp_data[elem_idx]))
    elif pp_option == "piezo":
        piezo_line_coords = fem_data.get("piezo_line_coords", None)
        gamma_water = fem_data.get("gamma_water", 9.81)
        if piezo_line_coords:
            piezo_line = LineString(piezo_line_coords)
            for elem_idx in range(n_elements):
                elem_type = element_types[elem_idx]
                elem_nodes_idx = elements[elem_idx][:elem_type]
                elem_coords = nodes[elem_nodes_idx]
                gp_u_list = []
                for gp_data in elem_gp_data[elem_idx]:
                    N = gp_data['N']
                    x_gp = N @ elem_coords[:, 0]
                    y_gp = N @ elem_coords[:, 1]
                    gp_point = Point(x_gp, y_gp)
                    closest_pt = piezo_line.interpolate(piezo_line.project(gp_point))
                    piezo_elev = closest_pt.y
                    gp_u_list.append(max(0.0, gamma_water * (piezo_elev - y_gp)))
                u_gp.append(gp_u_list)
        else:
            for elem_idx in range(n_elements):
                u_gp.append([0.0] * len(elem_gp_data[elem_idx]))
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
    else:
        for elem_idx in range(n_elements):
            u_gp.append([0.0] * len(elem_gp_data[elem_idx]))

    if debug_level >= 1 and pp_option != "none":
        all_u = [u_gp[e][g] for e in range(n_elements) for g in range(len(u_gp[e]))]
        max_u = max(all_u) if all_u else 0.0
        n_nonzero = sum(1 for v in all_u if v > 0.0)
        print(f"  Pore pressure ({pp_option}): max u_gp = {max_u:.3f}, {n_nonzero}/{len(all_u)} GPs with u > 0")

    # ---- Step 6c: Flatten Gauss-point data into per-element-type groups for
    # vectorized iteration (the per-GP Python loop dominated runtime; batched
    # numpy is 10-50x faster). Groups are needed because the B-matrix width
    # differs by element type. ----
    gp_groups = []
    _by_ndof = {}
    for _e in range(n_elements):
        for _g, _gpd in enumerate(elem_gp_data[_e]):
            _nd = len(_gpd['dof_indices'])
            _by_ndof.setdefault(_nd, []).append((_e, _g))
    for _nd, _pairs in _by_ndof.items():
        G = len(_pairs)
        grp = {
            'pairs': _pairs,
            'B': np.array([elem_gp_data[e][g]['B'] for e, g in _pairs]),
            'D4': np.array([elem_gp_data[e][g]['D4'] for e, g in _pairs]),
            'w': np.array([elem_gp_data[e][g]['weight'] for e, g in _pairs]),
            'dof': np.array([elem_gp_data[e][g]['dof_indices'] for e, g in _pairs], dtype=int),
            'c_r': np.array([c_reduced[e] for e, g in _pairs]),
            'phi_r': np.array([phi_reduced[e] for e, g in _pairs]),
            'dt_t': np.array([dt_t[e] for e, g in _pairs]),
            'evp': np.zeros((G, 4)),
        }
        grp['snph'] = np.sin(grp['phi_r'])
        grp['csph'] = np.cos(grp['phi_r'])
        gp_groups.append(grp)
    n_total_gp = sum(len(g['pairs']) for g in gp_groups)


    # ---- Step 7: Displacement-limit setup ----
    mesh_height = float(np.max(nodes[:, 1]) - np.min(nodes[:, 1]))
    if max_disp_factor is not None and mesh_height > 0:
        vp_disp_limit = max_disp_factor * mesh_height
    else:
        vp_disp_limit = None
    if debug_level >= 1 and vp_disp_limit is not None:
        print(f"  VP displacement limit: {vp_disp_limit:.2f} ({max_disp_factor:.0%} of mesh height {mesh_height:.1f})")

    # ---- Step 8: Initialize viscoplastic strains (zero) ----
    # evp[elem_idx][gp_idx] = array of shape (4,): [ex, ey, gxy, ez]
    # (4-component plane strain after Smith & Griffiths: plastic eps_z is
    # tracked so sigma_z can relax; total eps_z = 0 so elastic eps_z = -evp_z)
    evp = []
    for elem_idx in range(n_elements):
        n_gp = len(elem_gp_data[elem_idx])
        evp.append([np.zeros(4) for _ in range(n_gp)])


    # ---- Step 9: Viscoplastic iteration loop (optionally staged) ----
    # Staged loading: stage 1 applies gravity only (dry); stage 2 adds the
    # water loads (applied boundary forces, e.g. reservoir pressure) and the
    # pore-pressure field, continuing from the converged stage-1 state with
    # accumulated viscoplastic strains. This follows construction history
    # (dam built, then reservoir filled) and avoids the spurious effective-
    # tension zone produced by one-shot gravity+water elastic loading.
    water_present = bool(np.any(bc_type == 4)) or pp_option != "none"
    if staged and water_present:
        u_gp_dry = [[0.0] * len(g) for g in u_gp]
        stage_list = [(F_grav_pure, u_gp_dry, 'stage 1: gravity (dry)'),
                      (F_gravity, u_gp, 'stage 2: + water loads and pore pressures')]
    else:
        stage_list = [(F_gravity, u_gp, None)]

    total_iterations = 0
    u = np.zeros(n_dof)
    converged = False
    iteration = 0
    unbalanced_force_ratio = 0.0

    for stage_idx, (base_loads, u_gp_active, stage_label) in enumerate(stage_list):
        if debug_level >= 1 and stage_label is not None:
            print(f"  {stage_label}")

        # flatten this stage's per-GP pore pressures into the groups
        for grp in gp_groups:
            grp['u_gp'] = np.array([u_gp_active[e][g] for e, g in grp['pairs']])

        if pp_formulation == 'effective':
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

        # per-stage elastic reference (pure elastic response to this stage's loads)
        u_e_free = K_factor.solve(base_loads[free_dofs])
        u_elastic = np.zeros(n_dof)
        u_elastic[free_dofs] = u_e_free
        if stage_idx == 0:
            u = u_elastic.copy()   # start from the elastic solution
        if debug_level >= 1 and stage_idx == 0:
            print(f"  Initial elastic: max|u| = {np.max(np.abs(u)):.6f}")

        norm_F_gravity = np.linalg.norm(base_loads)
        converged = False
        unbalanced_force_ratio = 0.0
        ufr_prev = 0.0       # previous-iteration UFR
        ufr_rate_peak = 0.0  # per-solve peak of dUFR (reference for settled test)
        ufr_rate_best = float('inf')   # lowest dUFR since the peak
        last_progress_iter = 0         # iteration of last meaningful dUFR improvement

        for iteration in range(max_iterations):
            # Build body load correction from accumulated viscoplastic strains
            loads = base_loads.copy()

            n_yielding = 0

            for grp in gp_groups:
                Bg, D4g, wg = grp['B'], grp['D4'], grp['w']
                dofg, evpg = grp['dof'], grp['evp']
                u_e = u[dofg]                                   # (G, ndof)
                eps = np.einsum('gij,gj->gi', Bg, u_e)          # (G, 3)
                eps4 = np.empty((len(wg), 4))
                eps4[:, :3] = eps - evpg[:, :3]
                eps4[:, 3] = -evpg[:, 3]
                sig4 = np.einsum('gij,gj->gi', D4g, eps4)       # (G, 4) tension-positive
                sig_eff = sig4.copy()
                sig_eff[:, [0, 1, 3]] += grp['u_gp'][:, None]

                sx, sy, txy, sz = sig_eff.T
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
                sq3 = np.sqrt(3.0)
                f = (sigm * grp['snph']
                     + dsbar * (csth / sq3 - snth * grp['snph'] / 3.0)
                     - grp['c_r'] * grp['csph'])

                m = (f > 0) & (dsbar > 1e-20)
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

                if tension_cutoff:
                    tm = sigm > 0.0
                    if np.any(tm):
                        evpg[tm] += (sigm[tm] * grp['dt_t'][tm])[:, None] * \
                            np.array([1.0/3.0, 1.0/3.0, 0.0, 1.0/3.0])

                # body-load correction: B^T (D4 evp)[:3] * w, scattered to dofs
                s4 = np.einsum('gij,gj->gi', D4g, evpg)
                contrib = np.einsum('gij,gi->gj', Bg, s4[:, :3]) * wg[:, None]
                np.add.at(loads, dofg.ravel(), contrib.ravel())

            # ---- 1D Truss element body-force corrections ----
            if has_1d_elements:
                n_1d_compression = 0
                n_1d_exceeded = 0

                for elem_idx_1d in range(n_1d_elements):
                    if pile_elem_mask[elem_idx_1d]:
                        continue  # pile elements handled separately below
                    dof_idx = dof_indices_1d[elem_idx_1d]
                    k = k_by_1d_elem[elem_idx_1d]
                    cos_t = cos_theta_1d[elem_idx_1d]
                    sin_t = sin_theta_1d[elem_idx_1d]

                    # Relative displacement projected along element axis
                    u_elem = u[dof_idx]  # [u_x0, u_y0, u_x1, u_y1]
                    du_x = u_elem[2] - u_elem[0]
                    du_y = u_elem[3] - u_elem[1]
                    delta = du_x * cos_t + du_y * sin_t

                    # Axial force: T = k * delta (positive = tension)
                    T = k * delta
                    forces_1d[elem_idx_1d] = T

                    # Determine effective capacity
                    if failed_1d[elem_idx_1d]:
                        T_cap = t_res_by_1d_elem[elem_idx_1d]
                    else:
                        T_cap = t_allow_by_1d_elem[elem_idx_1d]

                    # Check violations and compute correction
                    correction_T = 0.0

                    if T < 0:
                        # Compression: cancel it entirely
                        correction_T = -T
                        n_1d_compression += 1
                    elif T > T_cap:
                        if not failed_1d[elem_idx_1d]:
                            # First time exceeding T_allow: mark as failed
                            failed_1d[elem_idx_1d] = True
                            T_cap = t_res_by_1d_elem[elem_idx_1d]
                        # Reduce force to T_cap
                        correction_T = T_cap - T
                        n_1d_exceeded += 1

                    if abs(correction_T) > 1e-30:
                        # Convert axial correction to nodal forces:
                        # Internal force pattern for tension T: [-cos, -sin, +cos, +sin]
                        # correction_T acts as additional axial force along element
                        loads[dof_idx[0]] += correction_T * (-cos_t)
                        loads[dof_idx[1]] += correction_T * (-sin_t)
                        loads[dof_idx[2]] += correction_T * cos_t
                        loads[dof_idx[3]] += correction_T * sin_t

                if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                    print(f"    1D elements: {n_1d_compression} in compression, "
                          f"{n_1d_exceeded} exceeded capacity, "
                          f"{np.sum(failed_1d)} total failed")

            # ---- Pile beam element force computation and capacity checks (6-DOF) ----
            if has_pile_elements:
                n_pile_yielded_V = 0
                n_pile_yielded_M = 0
                for p_idx in range(n_pile_elements):
                    dof_idx = dof_indices_pile[p_idx]
                    cos_t = cos_theta_pile[p_idx]
                    sin_t = sin_theta_pile[p_idx]
                    L = L_pile_elem[p_idx]
                    EI_val = EI_pile[p_idx]
                    EA_val = EA_pile[p_idx]

                    # Extract 6 global DOFs: [ux1, uy1, theta1, ux2, uy2, theta2]
                    u_elem = u[dof_idx]

                    # Transform to local coordinates using rotation matrix T
                    c = cos_t
                    s = sin_t
                    T = np.array([
                        [ c,  s, 0, 0, 0, 0],
                        [-s,  c, 0, 0, 0, 0],
                        [ 0,  0, 1, 0, 0, 0],
                        [ 0,  0, 0, c, s, 0],
                        [ 0,  0, 0,-s, c, 0],
                        [ 0,  0, 0, 0, 0, 1],
                    ])
                    u_local = T @ u_elem  # [u1_axial, v1_trans, theta1, u2_axial, v2_trans, theta2]

                    # Axial force: T = EA/L * (u2_axial - u1_axial)
                    T_force = EA_val / L * (u_local[3] - u_local[0])
                    forces_pile_axial[p_idx] = T_force

                    # Shear force: V = dM/dx, from beam theory
                    # V = 12*EI/L^3 * (v1 - v2) + 6*EI/L^2 * (theta1 + theta2)
                    L2 = L * L
                    L3 = L2 * L
                    V = 12*EI_val/L3 * (u_local[1] - u_local[4]) + 6*EI_val/L2 * (u_local[2] + u_local[5])

                    # Bending moments at node 1 and node 2 from K_local rows 2 and 5
                    M1 = EI_val * (6.0/L2 * u_local[1] + 4.0/L * u_local[2] - 6.0/L2 * u_local[4] + 2.0/L * u_local[5])
                    M2 = EI_val * (6.0/L2 * u_local[1] + 2.0/L * u_local[2] - 6.0/L2 * u_local[4] + 4.0/L * u_local[5])

                    forces_pile_moment[p_idx] = [M1, M2]

                    # --- V_cap check ---
                    V_limit = V_cap_pile[p_idx]
                    correction_V = 0.0
                    if abs(V) > V_limit:
                        correction_V = np.sign(V) * V_limit - V
                        V = np.sign(V) * V_limit
                        yielded_pile_V[p_idx] = True
                        n_pile_yielded_V += 1

                    forces_pile_lateral[p_idx] = V

                    if abs(correction_V) > 1e-30:
                        # Convert lateral correction to global nodal forces
                        # In local coords, shear internal force pattern: [0, 1, 0, 0, -1, 0]
                        # Transform to global: correction_V * T^T @ [0, 1, 0, 0, -1, 0]
                        f_local = np.array([0.0, correction_V, 0.0, 0.0, -correction_V, 0.0])
                        f_global = T.T @ f_local
                        for k in range(6):
                            loads[dof_idx[k]] += f_global[k]

                    # --- M_cap check (plastic hinge at each node) ---
                    M_cap_uw = M_cap_pile[p_idx]
                    if M_cap_uw < float('inf'):
                        # Node 1 moment check
                        if abs(M1) > M_cap_uw:
                            correction_M1 = np.sign(M1) * M_cap_uw - M1
                            rot_dof_1 = dof_idx[2]  # theta1 DOF
                            loads[rot_dof_1] += correction_M1
                            yielded_pile_M[p_idx] = True
                            n_pile_yielded_M += 1

                        # Node 2 moment check
                        if abs(M2) > M_cap_uw:
                            correction_M2 = np.sign(M2) * M_cap_uw - M2
                            rot_dof_2 = dof_idx[5]  # theta2 DOF
                            loads[rot_dof_2] += correction_M2
                            yielded_pile_M[p_idx] = True

                if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                    if n_pile_yielded_V > 0 or n_pile_yielded_M > 0:
                        print(f"    Pile elements: {n_pile_yielded_V} V-yielded, {n_pile_yielded_M} M-yielded")

            # Compute unbalanced force ratio: ||VP corrections|| / ||F_gravity||
            # This measures how far the system is from force equilibrium.
            if norm_F_gravity > 1e-30:
                unbalanced_force_ratio = np.linalg.norm(loads - base_loads) / norm_F_gravity
            else:
                unbalanced_force_ratio = np.linalg.norm(loads - base_loads)

            # Solve K * u_new = loads
            loads_free = loads[free_dofs]
            u_free_new = K_factor.solve(loads_free)

            u_new = np.zeros(n_dof)
            u_new[free_dofs] = u_free_new

            # Convergence check: max|du| / max|u| < tolerance — Smith & Griffiths'
            # CHECON test (infinity norms), exactly as in p62.f90. The max-norm
            # matters: a small cluster of Gauss points sustaining benign localized
            # creep (possible under pp_formulation='total', whose subtract-at-GP
            # recipe leaves mild effective tension at submerged boundaries)
            # produces a bounded per-iteration du measured against the GLOBAL
            # maximum displacement, so the ratio decays and the stable state is
            # accepted within the iteration ceiling. A true failure mechanism
            # keeps feeding the global du and stays above tolerance. False
            # convergence from large failure
            # displacements is guarded by the max_disp_factor limit below.
            norm_diff = np.max(np.abs(u_new - u))
            norm_u_new = np.max(np.abs(u_new))

            if norm_u_new > 1e-30:
                relative_change = norm_diff / norm_u_new
            else:
                relative_change = norm_diff

            # Plastic-settled condition: the rate of change of the accumulated
            # body-load vector (dUFR) decays to zero when viscoplastic flow has
            # genuinely ceased, while a failing (permanently creeping) state
            # plateaus at a constant rate. The threshold is PEAK-RELATIVE — each
            # solve provides its own reference — so the test is dimensionless
            # and scale-free (independent of units, domain size, and yielding-
            # zone fraction). Measured separation between settled and creeping
            # states is ~700x; the 1%-of-peak threshold has wide margin.
            ufr_rate = abs(unbalanced_force_ratio - ufr_prev)
            ufr_prev = unbalanced_force_ratio
            if ufr_rate > ufr_rate_peak:
                ufr_rate_peak = ufr_rate
            plastic_settled = ufr_rate <= max(0.01 * ufr_rate_peak, 1e-14)

            # No-progress early exit: a settling state's dUFR keeps decaying
            # toward the threshold; a failing state's dUFR plateaus. If 500
            # iterations pass with no meaningful improvement (>1%) of the best
            # dUFR seen, the trial cannot settle - declare failure early
            # rather than burning the ceiling. Disabled in displacement-limit
            # mode, where failure is defined by the limit trip itself.
            if ufr_rate < 0.99 * ufr_rate_best:
                ufr_rate_best = ufr_rate
                last_progress_iter = iteration
            # Window calibration: genuinely settling states can stall for
            # >500 iterations mid-decay (reinforced slope at F=1.6 settles at
            # ~2900 iters with a ~1000-iter dUFR plateau on the way), so the
            # window must be generous; post-vectorization the extra iterations
            # cost seconds.
            if (early_exit and not plastic_settled
                    and iteration - last_progress_iter > 1500):
                converged = False
                u = u_new
                if debug_level >= 1:
                    print(f"  Early exit at iteration {iteration+1}: no dUFR "
                          f"progress for 1500 iterations (plateau "
                          f"{ufr_rate:.2e} vs settled target "
                          f"{0.01*ufr_rate_peak:.2e}) - declared FAILED")
                break

            if debug_level >= 2 and (iteration % 10 == 0 or iteration < 5):
                print(f"  Iter {iteration+1:4d}: max|du|/max|u| = {relative_change:.3e}, "
                      f"UFR = {unbalanced_force_ratio:.3e}, dUFR = {ufr_rate:.2e}, "
                      f"yielding = {n_yielding}/{n_total_gp}, max|u| = {np.max(np.abs(u_new)):.6f}")

            # Displacement limit check: detect false convergence from unbounded plastic flow
            # When VP displacements exceed a fraction of mesh height, the slope has physically
            # failed even if the relative convergence criterion is satisfied (see FLAC manual;
            # Griffiths & Lane 1999 displacement-vs-F plots).
            if vp_disp_limit is not None:
                u_vp = u_new - u_elastic
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
                    u = u_new
                    if debug_level >= 1:
                        print(f"  Displacement limit exceeded at iteration {iteration+1}: "
                              f"max VP disp = {max_vp_disp:.2f} > limit {vp_disp_limit:.2f}")
                    break

            if relative_change < tolerance and plastic_settled:
                converged = True
                u = u_new
                if debug_level >= 1:
                    print(f"  Converged after {iteration+1} iterations "
                          f"(max|du|/max|u| = {relative_change:.3e}, "
                          f"dUFR = {ufr_rate:.2e} <= 1% of peak {ufr_rate_peak:.2e})")
                break

            u = u_new


        total_iterations += iteration + 1
        if not converged:
            break   # stage failed -> overall failure

    if not converged and debug_level >= 1:
        print(f"  Did NOT converge after {max_iterations} iterations (max|du|/max|u| = {relative_change:.3e})")

    # Copy grouped viscoplastic strains back into the per-element list used by
    # the post-processing blocks below.
    for grp in gp_groups:
        for k, (e, g) in enumerate(grp['pairs']):
            evp[e][g][:] = grp['evp'][k]

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

        stress_avg_tp4 /= n_gp
        u_elem_avg = sum(u_gp[elem_idx]) / len(u_gp[elem_idx]) if u_gp[elem_idx] else 0.0
        if pp_formulation == 'effective':
            # D*B*u is effective; report total stresses for output (legacy
            # convention) and use the stresses as-is for the yield check.
            stress_total_tp4 = stress_avg_tp4 - np.array(
                [u_elem_avg, u_elem_avg, 0.0, u_elem_avg])
            sig_eff4 = stress_avg_tp4
        else:
            stress_total_tp4 = stress_avg_tp4
            sig_eff4 = stress_avg_tp4 + np.array(
                [u_elem_avg, u_elem_avg, 0.0, u_elem_avg])
        # compression-positive in-plane components for output; sig_vm = dsbar
        sig_x, sig_y = -stress_total_tp4[0], -stress_total_tp4[1]
        tau_xy = stress_total_tp4[2]
        _, sig_vm, _ = stress_invariants(stress_total_tp4)
        final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, sig_vm]

        sigm, dsbar, theta = stress_invariants(sig_eff4)
        f_yield = mc_yield_invariants(sigm, dsbar, theta,
                                      c_reduced[elem_idx], phi_reduced[elem_idx])
        yield_function_out[elem_idx] = f_yield
        plastic_elements[elem_idx] = f_yield > 1e-8

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
            dof_idx = dof_indices_1d[elem_idx_1d]
            k = k_by_1d_elem[elem_idx_1d]
            cos_t = cos_theta_1d[elem_idx_1d]
            sin_t = sin_theta_1d[elem_idx_1d]

            u_elem = u[dof_idx]
            du_x = u_elem[2] - u_elem[0]
            du_y = u_elem[3] - u_elem[1]
            delta = du_x * cos_t + du_y * sin_t
            T = k * delta

            # Apply constraints to final reported force
            if failed_1d[elem_idx_1d]:
                T = min(T, t_res_by_1d_elem[elem_idx_1d])
            else:
                T = min(T, t_allow_by_1d_elem[elem_idx_1d])
            T = max(T, 0.0)  # No compression

            forces_1d[elem_idx_1d] = T

    # ---- Step 10c: Compute final pile beam element forces (capped at capacity) ----
    if has_pile_elements:
        for p_idx in range(n_pile_elements):
            dof_idx = dof_indices_pile[p_idx]
            cos_t = cos_theta_pile[p_idx]
            sin_t = sin_theta_pile[p_idx]
            L = L_pile_elem[p_idx]
            EI_val = EI_pile[p_idx]
            EA_val = EA_pile[p_idx]

            u_elem = u[dof_idx]

            c = cos_t
            s = sin_t
            T = np.array([
                [ c,  s, 0, 0, 0, 0],
                [-s,  c, 0, 0, 0, 0],
                [ 0,  0, 1, 0, 0, 0],
                [ 0,  0, 0, c, s, 0],
                [ 0,  0, 0,-s, c, 0],
                [ 0,  0, 0, 0, 0, 1],
            ])
            u_local = T @ u_elem

            # Axial force
            T_force = EA_val / L * (u_local[3] - u_local[0])
            forces_pile_axial[p_idx] = T_force

            # Shear force
            L2 = L * L
            L3 = L2 * L
            V = 12*EI_val/L3 * (u_local[1] - u_local[4]) + 6*EI_val/L2 * (u_local[2] + u_local[5])

            # Bending moments
            M1 = EI_val * (6.0/L2 * u_local[1] + 4.0/L * u_local[2] - 6.0/L2 * u_local[4] + 2.0/L * u_local[5])
            M2 = EI_val * (6.0/L2 * u_local[1] + 2.0/L * u_local[2] - 6.0/L2 * u_local[4] + 4.0/L * u_local[5])
            forces_pile_moment[p_idx] = [M1, M2]

            # Cap shear at V_cap
            V_limit = V_cap_pile[p_idx]
            if abs(V) > V_limit:
                V = np.sign(V) * V_limit
                yielded_pile_V[p_idx] = True
            forces_pile_lateral[p_idx] = V

            # Cap moments at M_cap
            M_cap_uw = M_cap_pile[p_idx]
            if M_cap_uw < float('inf'):
                if abs(M1) > M_cap_uw or abs(M2) > M_cap_uw:
                    yielded_pile_M[p_idx] = True
                M1_capped = np.sign(M1) * min(abs(M1), M_cap_uw) if M_cap_uw < float('inf') else M1
                M2_capped = np.sign(M2) * min(abs(M2), M_cap_uw) if M_cap_uw < float('inf') else M2
                forces_pile_moment[p_idx] = [M1_capped, M2_capped]

    n_plastic = np.sum(plastic_elements)
    if debug_level >= 1:
        print(f"  Plastic elements: {n_plastic}/{n_elements}")
        print(f"  Max displacement: {np.max(np.abs(u)):.6f}")
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
        "iterations": total_iterations,
        "displacements": u,
        "displacements_elastic": u_elastic,
        "stresses": final_stresses,
        "strains": strains,
        "vp_shear_strain": vp_shear_strain,
        "plastic_elements": plastic_elements,
        "yield_function": yield_function_out,
        "max_displacement": np.max(np.abs(u)),
        "plastic_strains": {i: np.array(evp[i]) for i in range(n_elements)},
        "algorithm": "Griffiths & Lane (1999) Viscoplastic",
        "F": F,
        "residual": relative_change if 'relative_change' in locals() else 0.0,
        "unbalanced_force_ratio": unbalanced_force_ratio,
        "plastic_fraction": n_plastic / n_elements if n_elements > 0 else 0.0,
        "forces_1d": forces_1d if has_1d_elements else np.array([]),
        "failed_1d_elements": failed_1d if has_1d_elements else np.array([], dtype=bool),
        "forces_pile_axial": forces_pile_axial if has_pile_elements else np.array([]),
        "forces_pile_lateral": forces_pile_lateral if has_pile_elements else np.array([]),
        "forces_pile_moment": forces_pile_moment if has_pile_elements else np.zeros((0, 2)),
        "yielded_pile_V": yielded_pile_V if has_pile_elements else np.array([], dtype=bool),
        "yielded_pile_M": yielded_pile_M if has_pile_elements else np.array([], dtype=bool),
        "yielded_pile": (yielded_pile_V | yielded_pile_M) if has_pile_elements else np.array([], dtype=bool),
    }


def print_reinforcement_summary(fem_data, solution):
    """
    Print a summary table of reinforcement line results.

    Groups 1D elements by reinforcement line and reports per-line statistics
    including element counts, force ranges, and failure modes.
    """
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

    # Filter out pile elements — reinforcement reported separately
    reinf_mask = ~pile_elem_mask
    if not np.any(reinf_mask):
        return

    # Get per-line Tmax and Tres from slope_data stored in fem_data
    # element_materials_1d is 1-based line ID
    line_ids = np.unique(element_materials_1d[reinf_mask])

    print("\n=== Reinforcement Summary ===")
    print(f"{'Line':>4}  {'Elems':>5}  {'Max T':>8}  {'Avg T':>8}  "
          f"{'Tension':>7}  {'In Lp':>5}  {'At Tres':>7}  {'Broken':>6}  {'Status'}")
    print("-" * 80)

    statuses_seen = set()
    for line_id in sorted(line_ids):
        mask = element_materials_1d == line_id
        n_elem = int(mask.sum())
        line_forces = forces[mask]
        line_failed = failed[mask]
        line_t_allow = t_allow_by_elem[mask]
        line_t_res = t_res_by_elem[mask]

        # Max Tallow for this line (the full-capacity elements)
        t_max_line = line_t_allow.max() if n_elem > 0 else 0.0

        # Active: elements carrying tension
        n_active = int((line_forces > 0).sum())

        # Pullout zone: elements where T_allow < T_max (reduced by proximity to end)
        n_pullout = int(((line_t_allow < t_max_line - 1e-6) & (line_t_allow > 1e-6)).sum())

        # At Tres: failed elements still carrying residual force
        n_at_tres = int((line_failed & (line_t_res > 1e-6)).sum())

        # Broken: failed elements with zero residual (complete failure)
        broken_mask = line_failed & (line_t_res < 1e-6)
        n_broken = int(broken_mask.sum())

        # Determine if broken elements are in Lp zone or outside
        in_lp_mask = (line_t_allow < t_max_line - 1e-6) & (line_t_allow > 1e-6)
        # Elements at the very end (t_allow ~ 0) are also in Lp zone
        at_end_mask = line_t_allow < 1e-6
        lp_zone_mask = in_lp_mask | at_end_mask
        n_broken_in_lp = int((broken_mask & lp_zone_mask).sum())
        n_broken_outside_lp = n_broken - n_broken_in_lp

        max_t = line_forces.max() if n_elem > 0 else 0.0
        active_forces = line_forces[line_forces > 0]
        avg_t = active_forces.mean() if len(active_forces) > 0 else 0.0

        # Status
        if n_broken == n_elem:
            status = "RUPTURED"
        elif n_broken_outside_lp > 0:
            status = "RUPTURED"
        elif n_at_tres > 0:
            status = "YIELDED"
        elif n_broken_in_lp > 0:
            status = "PULLOUT"
        elif max_t > 0.95 * t_max_line and t_max_line > 0:
            status = "NEAR CAPACITY"
        elif n_active > 0:
            status = "OK"
        else:
            status = "INACTIVE"

        statuses_seen.add(status)

        print(f"{line_id:>4}  {n_elem:>5}  {max_t:>8.1f}  {avg_t:>8.1f}  "
              f"{n_active:>7}  {n_pullout:>5}  {n_at_tres:>7}  {n_broken:>6}  {status}")

    print("-" * 80)

    # Print notes for statuses that appeared
    status_notes = {
        "OK": "OK: All elements within allowable capacity, no failures.",
        "NEAR CAPACITY": "NEAR CAPACITY: Maximum force exceeds 95% of Tmax. Close to yielding.",
        "PULLOUT": "PULLOUT: Elements near the reinforcement ends (within Lp) have failed due to insufficient embedment length. Interior elements are intact.",
        "YIELDED": "YIELDED: One or more elements have exceeded Tallow and dropped to residual capacity Tres. The line is still carrying load at reduced strength.",
        "RUPTURED": "RUPTURED: One or more elements outside the pullout zone have broken (T exceeded Tallow with Tres=0). The reinforcement line has lost structural continuity.",
        "INACTIVE": "INACTIVE: No elements are carrying tension. The reinforcement is not engaged.",
    }
    notes = [status_notes[s] for s in ["OK", "NEAR CAPACITY", "PULLOUT", "YIELDED", "RUPTURED", "INACTIVE"] if s in statuses_seen]
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

            if is_failed and t_res < 1e-6:
                # Distinguish pullout (near ends, reduced t_allow) from rupture (full capacity exceeded)
                # Get max t_allow for this line to determine if element is in pullout zone
                line_mask = element_materials_1d == line_id
                t_max_line = t_allow_by_elem[line_mask].max()
                if t_allow < t_max_line - 1e-6:
                    status = "PULLOUT"
                else:
                    status = "BROKEN"
            elif is_failed:
                status = "AT TRES"
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


def solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.01, debug_level=0,
               max_iterations=3000, convergence_tol=1e-3, max_disp_factor=0.1,
               failure_criterion="non_convergence", n_sweep=10,
               staged=False, tension_cutoff=False, char_point=None,
               pp_formulation='effective', dt_scale=1.0):
    """
    Shear Strength Reduction Method using bisection on solve_fem convergence.

    Finds the critical strength reduction factor F where the viscoplastic
    algorithm transitions from converging (stable) to non-converging (failure).

    Parameters:
        fem_data (dict): FEM data from build_fem_data
        F_min (float): Lower bound for F (must converge). Default 1.0.
        F_max (float): Upper bound for F (should not converge). Default 2.0.
        tolerance (float): Bisection stops when F_right - F_left < tolerance. Default 0.01.
            The reported FS is the midpoint of the final bracket (+/- tolerance/2);
            the bracket itself is returned in 'final_interval'.
        debug_level (int): Verbosity (0=silent, 1=summary, 2=detailed)
        max_iterations (int): Max viscoplastic iterations passed to solve_fem
        convergence_tol (float): Convergence tolerance passed to solve_fem
        max_disp_factor (float): Displacement limit (fraction of mesh height) used as a
            backstop/early-termination cap inside solve_fem trials (default 0.1).
        failure_criterion (str): How to determine failure.
            "non_convergence" (default) - Bisection on TRUE viscoplastic
                equilibrium: a trial converges only if both the CHECON
                displacement test and the plastic-flow-stopped test (dUFR
                below 1% of its per-solve peak) are satisfied. Scale-free and
                insensitive to dt/tolerance/ceiling. Use for problems without
                reservoir loading.
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

    Returns:
        dict: Result with keys FS, converged, last_solution, final_interval, etc.
    """

    t_start = time.perf_counter()

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

    if failure_criterion == "non_convergence":
        result = _ssrm_displacement_limit(
            fem_data, F_min=F_min, F_max=F_max, tolerance=tolerance,
            debug_level=debug_level, max_iterations=max_iterations,
            convergence_tol=convergence_tol, max_disp_factor=None, staged=staged,
            tension_cutoff=tension_cutoff, pp_formulation=pp_formulation, dt_scale=dt_scale)
    elif failure_criterion == "displacement_limit":
        result = _ssrm_displacement_limit(
            fem_data, F_min=F_min, F_max=F_max, tolerance=tolerance,
            debug_level=debug_level, max_iterations=max_iterations,
            convergence_tol=convergence_tol, max_disp_factor=max_disp_factor,
            staged=staged, tension_cutoff=tension_cutoff, pp_formulation=pp_formulation, dt_scale=dt_scale)
    elif failure_criterion == "displacement_increase":
        result = _ssrm_displacement_increase(
            fem_data, F_min=F_min, F_max=F_max, tolerance=tolerance,
            debug_level=debug_level, max_iterations=max_iterations,
            convergence_tol=convergence_tol, n_sweep=n_sweep,
            tension_cutoff=tension_cutoff, char_point=char_point, pp_formulation=pp_formulation, dt_scale=dt_scale)
    else:
        raise ValueError(
            f"Unknown failure_criterion '{failure_criterion}'. Supported: "
            "'non_convergence' (default; bisection on true viscoplastic "
            "equilibrium), 'displacement_limit' (displacement-budget "
            "backstop), and "
            "'displacement_increase' (displacement-catastrophe sweep) - see "
            "docs/fem/overview.md, 'Choosing a Failure Criterion'.")

    elapsed = time.perf_counter() - t_start
    result["elapsed_time"] = elapsed
    if debug_level >= 1:
        print(f"  SSRM completed in {elapsed:.1f} seconds")

    return result


def _ssrm_displacement_limit(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05,
                              debug_level=0, max_iterations=500,
                              convergence_tol=1e-3, max_disp_factor=0.1,
                              staged=False, tension_cutoff=False,
                 pp_formulation='effective',
                 dt_scale=1.0):
    """SSRM using fixed VP displacement limit as failure criterion."""

    if debug_level >= 1:
        label = "Non-Convergence" if max_disp_factor is None else "Displacement Limit"
        print(f"=== SSRM Analysis ({label} Method) ===")
        print(f"  Bisection range: [{F_min:.2f}, {F_max:.2f}], tolerance: {tolerance}")
        if max_disp_factor is not None:
            print(f"  Displacement limit: {max_disp_factor:.0%} of mesh height")

    F_left = F_min
    F_right = F_max

    # Verify lower bound converges
    if debug_level >= 1:
        print(f"  Verifying lower bound F={F_min:.2f} converges...")
    solution_min = solve_fem(fem_data, F=F_min, debug_level=max(0, debug_level-1), dt_scale=dt_scale, pp_formulation=pp_formulation,
                             max_iterations=max_iterations, tolerance=convergence_tol,
                             max_disp_factor=max_disp_factor, staged=staged,
                             tension_cutoff=tension_cutoff,
                             early_exit=(max_disp_factor is None))
    if not solution_min["converged"]:
        msg = (f"SSRM: the slope does not reach equilibrium even at F = {F_min:.2f} — it is "
               f"unstable at or below this strength-reduction factor (FS < {F_min:.2f}). "
               "Lower F_min to bracket the factor of safety.")
        print(f"\n{msg}")
        return {
            "converged": False,
            "error": msg,
            "FS": None
        }
    if debug_level >= 1:
        print(f"    -> Converged in {solution_min['iterations']} iters")

    # Verify upper bound does not converge
    if debug_level >= 1:
        print(f"  Verifying upper bound F={F_max:.2f} does not converge...")
    solution_max = solve_fem(fem_data, F=F_max, debug_level=max(0, debug_level-1), dt_scale=dt_scale, pp_formulation=pp_formulation,
                             max_iterations=max_iterations, tolerance=convergence_tol,
                             max_disp_factor=max_disp_factor, staged=staged,
                             tension_cutoff=tension_cutoff,
                             early_exit=(max_disp_factor is None))
    if solution_max["converged"]:
        msg = (f"SSRM: the slope still reaches equilibrium at F = {F_max:.2f}, so the factor "
               f"of safety exceeds the search ceiling (FS > {F_max:.2f}). No failure was found "
               f"in [{F_min:.2f}, {F_max:.2f}] — increase F_max to bracket it. This is expected "
               "for a heavily reinforced or very stable slope that deforms ductilely without a "
               "displacement catastrophe.")
        print(f"\n{msg}")
        return {
            "converged": False,
            "FS": None,
            "last_solution": solution_max,
            "error": msg
        }
    if debug_level >= 1:
        print(f"    -> Did NOT converge ({solution_max['iterations']} iters)")

    last_converged_solution = solution_min
    iteration = 0

    while (F_right - F_left) > tolerance and iteration < 50:
        F_mid = (F_left + F_right) / 2.0

        if debug_level >= 1:
            print(f"\n  SSRM step {iteration+1}: F = {F_mid:.4f}  [{F_left:.4f}, {F_right:.4f}]")

        solution = solve_fem(fem_data, F=F_mid, debug_level=max(0, debug_level-1), dt_scale=dt_scale, pp_formulation=pp_formulation,
                             max_iterations=max_iterations, tolerance=convergence_tol,
                             max_disp_factor=max_disp_factor, staged=staged,
                             tension_cutoff=tension_cutoff,
                             early_exit=(max_disp_factor is None))

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

    # Report the midpoint of the final bracket (unbiased, +/- tolerance/2);
    # the full bracket is returned in 'final_interval'.
    critical_FS = 0.5 * (F_left + F_right)

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
        "method": "SSRM — Non-Convergence (Griffiths & Lane 1999)" if max_disp_factor is None else "SSRM — Displacement Limit"
    }


def _ssrm_displacement_increase(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05,
                                 debug_level=0, max_iterations=500,
                                 convergence_tol=1e-3, n_sweep=10,
                                 tension_cutoff=False, char_point=None,
                 pp_formulation='effective',
                 dt_scale=1.0):
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

    def _get_max_vp_disp(F_val):
        """Run solve_fem; return the displacement measure (characteristic-point
        VP displacement, or global max before the point is selected)."""
        sol = solve_fem(fem_data, F=F_val, debug_level=max(0, debug_level-1), dt_scale=dt_scale, pp_formulation=pp_formulation,
                        max_iterations=max_iterations, tolerance=convergence_tol,
                        max_disp_factor=early_term_factor,
                        tension_cutoff=tension_cutoff)
        # Use VP displacement (total - elastic) to isolate plastic deformation.
        # The elastic component is roughly constant regardless of F and masks
        # the catastrophic growth in plastic displacement at failure.
        return _measure(sol), sol

    # Phase 1: Coarse sweep
    F_values = np.linspace(F_min, F_max, n_sweep)
    displacements = []
    solutions = []

    if debug_level >= 1:
        print(f"\n  Phase 1: Coarse sweep ({n_sweep} points)")

    for i, F_val in enumerate(F_values):
        max_vp, sol = _get_max_vp_disp(F_val)
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
        F_mid = (F_left + F_right) / 2.0
        d_mid, sol_mid = _get_max_vp_disp(F_mid)

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


def build_global_stiffness(nodes, elements, element_types, element_materials, E_by_mat, nu_by_mat, fem_data=None):
    """
    Build global stiffness matrix from 2D soil elements and (optionally) 1D truss + pile beam elements.

    Uses dof_offset from fem_data to support mixed DOF systems (pile nodes have 3 DOFs).
    """
    n_nodes = len(nodes)

    # Get DOF offset map from fem_data if available
    dof_offset = fem_data.get("dof_offset", None) if fem_data is not None else None
    if dof_offset is not None:
        n_dof = int(dof_offset[n_nodes])
    else:
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

        # Assemble into global matrix using dof_offset for global DOF indices
        for i in range(elem_type):
            for j in range(elem_type):
                node_i = elem_nodes[i]
                node_j = elem_nodes[j]

                if dof_offset is not None:
                    base_i = dof_offset[node_i]
                    base_j = dof_offset[node_j]
                else:
                    base_i = 2 * node_i
                    base_j = 2 * node_j

                for di in range(2):
                    for dj in range(2):
                        global_i = base_i + di
                        global_j = base_j + dj
                        local_i = 2 * i + di
                        local_j = 2 * j + dj

                        if local_i < K_elem.shape[0] and local_j < K_elem.shape[1]:
                            K_global[global_i, global_j] += K_elem[local_i, local_j]

    # Assemble 1D truss element stiffness matrices (reinforcement only — skip pile elements)
    if fem_data is not None:
        K_global_1d_elems = fem_data.get("K_global_1d_elems", [])
        dof_indices_1d = fem_data.get("dof_indices_1d", np.zeros((0, 4), dtype=int))
        pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(K_global_1d_elems), dtype=bool))

        for elem_idx_1d in range(len(K_global_1d_elems)):
            if pile_elem_mask[elem_idx_1d]:
                continue  # pile elements use beam stiffness, assembled below
            K_elem_1d = K_global_1d_elems[elem_idx_1d]
            dof_idx = dof_indices_1d[elem_idx_1d]

            for i in range(4):
                for j in range(4):
                    K_global[dof_idx[i], dof_idx[j]] += K_elem_1d[i, j]

        # Assemble pile beam element stiffness matrices (6x6 Euler-Bernoulli)
        K_global_pile_elems = fem_data.get("K_global_pile_elems", [])
        dof_indices_pile = fem_data.get("dof_indices_pile", np.zeros((0, 6), dtype=int))

        for p_idx in range(len(K_global_pile_elems)):
            K_beam = K_global_pile_elems[p_idx]
            dof_idx = dof_indices_pile[p_idx]

            for i in range(6):
                for j in range(6):
                    K_global[dof_idx[i], dof_idx[j]] += K_beam[i, j]

    return K_global.tocsr()



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
    # Add numerical stability check for near-incompressible materials
    if nu >= 0.45:
        print(f"Warning: Poisson's ratio {nu:.3f} is close to incompressible limit (0.5)")
        print("Consider using nu <= 0.4 for better numerical stability")

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


def mc_flow_vector_4(stress4_tp, psi=0.0):
    """Plastic-potential gradient dQ/dsigma for the 4-component plane-strain
    stress (tension-positive), after Smith & Griffiths' MOCOUQ + FORMM.

    Q has the same invariant form as the yield function with phi replaced by
    the dilation angle psi. Near the Lode-angle corners (|sin(theta)| > 0.49,
    i.e. theta within ~0.7 deg of +/-30 deg) the theta-dependence is frozen at
    the corner value and the J3 term dropped — S&G's corner treatment, which
    keeps the flow direction finite where tan(3*theta) blows up.

    Returns the 4-component flow vector [dQ/dsx, dQ/dsy, dQ/dtxy, dQ/dsz]
    (engineering shear).
    """
    sx, sy, txy, sz = stress4_tp
    sigm, dsbar, theta = stress_invariants(stress4_tp)
    if dsbar < 1e-20:
        return np.zeros(4)

    snps = sin(psi)
    sq3 = sqrt(3.0)
    snth = sin(theta)

    # deviator components and J2
    dx = (2.0 * sx - sy - sz) / 3.0
    dy = (2.0 * sy - sz - sx) / 3.0
    dz = (2.0 * sz - sx - sy) / 3.0
    xj2 = dsbar**2 / 3.0

    # dQ/dsigma = dq1 * d(sigm)/dsig + C2 * d(dsbar)/dsig + C3 * d(J3)/dsig
    a1 = np.array([1.0/3.0, 1.0/3.0, 0.0, 1.0/3.0])
    a2 = (3.0 / (2.0 * dsbar)) * np.array([dx, dy, 2.0 * txy, dz])
    a3 = np.array([
        dx*dx + txy*txy - (2.0/3.0) * xj2,
        dy*dy + txy*txy - (2.0/3.0) * xj2,
        -2.0 * dz * txy,
        dz*dz - (2.0/3.0) * xj2,
    ])

    dq1 = snps
    if abs(snth) > 0.49:
        # corner: freeze K(theta) at theta = +/-30 deg, drop J3 term
        c1 = 1.0 if snth >= 0.0 else -1.0
        C2 = 0.5 - c1 * snps / 6.0          # K(+/-30 deg) = 1/2 -/+ snps/6
        C3 = 0.0
    else:
        csth = cos(theta)
        cs3th = cos(3.0 * theta)
        tn3th = tan(3.0 * theta)
        # K(theta) = cos(theta)/sqrt(3) - sin(theta)*snps/3
        # K'(theta) = -sin(theta)/sqrt(3) - cos(theta)*snps/3
        # From sin(3*theta) = -13.5*J3/dsbar^3:
        #   d(theta) = -4.5/(cos3th*dsbar^3) dJ3 - (tan3th/dsbar) d(dsbar)
        # so C2 = K - K'*tan(3*theta);  C3 = -4.5*K'/(cos(3*theta)*dsbar^2)
        K = csth / sq3 - snth * snps / 3.0
        Kp = -snth / sq3 - csth * snps / 3.0
        C2 = K - Kp * tn3th
        C3 = -4.5 * Kp / (cs3th * dsbar**2)

    return dq1 * a1 + C2 * a2 + C3 * a3



def check_mohr_coulomb_cp(stress_cp, c, phi, u=0.0):
    """Mohr-Coulomb yield function for compression-positive stresses with pore pressure.

    For compression-positive convention (compression > 0, tension < 0):
    F = tau_max - sigma'_mean * sin(phi) - c * cos(phi)

    Where:
    - tau_max = maximum shear stress = sqrt((sig_x - sig_y)^2/4 + tau_xy^2)
    - sigma'_mean = effective mean normal stress = (sig_x + sig_y)/2 - u
    - u = pore pressure (positive value reduces effective compressive stress)
    - Positive F indicates yielding

    Args:
        stress_cp: Array [sig_x, sig_y, tau_xy] in compression-positive convention
        c: Cohesion
        phi: Friction angle in radians
        u: Pore pressure (default 0.0). Positive u reduces effective stress.

    Returns:
        F: Yield function value (F > 0 means yielding)
    """
    sig_x, sig_y, tau_xy = stress_cp
    sig_mean = (sig_x + sig_y) / 2.0 - u  # effective mean stress
    tau_max = sqrt(((sig_x - sig_y) / 2.0)**2 + tau_xy**2)
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    F = tau_max - sig_mean * sin_phi - c * cos_phi
    return F



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
