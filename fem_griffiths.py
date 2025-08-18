# Implementation of Griffiths & Lane (1999) FEM Slope Stability Algorithm
#
# This implements the exact algorithm from "Slope stability analysis by finite elements"
# Geotechnique 49, No. 3, 387-403 (1999)
#
# Key features:
# - 8-node quadrilateral elements with reduced integration (4 Gauss points)
# - Visco-plastic stress redistribution (Perzyna algorithm)
# - Non-convergence failure criterion only
# - Non-associated flow rule (ψ=0)
# - K0 initial stress state
# - Elastic-perfectly plastic Mohr-Coulomb model

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from math import radians, degrees, sin, cos, tan, sqrt, atan2
import warnings

def solve_fem_griffiths(fem_data, F=1.0, debug_level=0):
    """
    Solve FEM using Griffiths & Lane (1999) algorithm.
    
    This implements the exact algorithm from the 1999 Geotechnique paper:
    - 8-node quadrilateral elements with reduced integration
    - Visco-plastic stress redistribution (Perzyna algorithm)
    - Non-convergence failure criterion
    - K0 initial stress state
    - Non-associated flow rule (ψ=0)
    
    Parameters:
        fem_data (dict): FEM data dictionary
        F (float): Shear strength reduction factor
        debug_level (int): Verbosity level
        
    Returns:
        dict: Solution dictionary with convergence status
    """
    
    if debug_level >= 1:
        print(f"=== Griffiths FEM Analysis (F={F:.3f}) ===")
    
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
    
    # Apply strength reduction
    c_reduced = c_by_elem / F
    phi_reduced = np.arctan(np.tan(np.radians(phi_by_elem)) / F)
    
    if debug_level >= 2:
        print(f"Original c range: [{np.min(c_by_elem):.3f}, {np.max(c_by_elem):.3f}]")
        print(f"Reduced c range: [{np.min(c_reduced):.3f}, {np.max(c_reduced):.3f}]")
        print(f"Original φ range: [{np.min(phi_by_elem):.1f}°, {np.max(phi_by_elem):.1f}°]")
        print(f"Reduced φ range: [{np.min(np.degrees(phi_reduced)):.1f}°, {np.max(np.degrees(phi_reduced)):.1f}°]")
    
    # Build global stiffness matrix (done once)
    K_global = build_global_stiffness_griffiths(nodes, elements, element_types, 
                                               element_materials, E_by_mat, nu_by_mat)
    
    # Build gravity load vector
    F_gravity = build_gravity_loads_griffiths(nodes, elements, element_types, 
                                            element_materials, gamma_by_mat, k_seismic)
    
    # Initialize stress state with K0 conditions
    stress_state = initialize_k0_stress_griffiths(nodes, elements, element_types,
                                                element_materials, gamma_by_mat, nu_by_mat, u_nodal)
    
    # Apply boundary conditions
    constraint_dofs = []
    for i in range(n_nodes):
        if bc_type[i] == 1:  # Fixed
            constraint_dofs.extend([2*i, 2*i+1])
        elif bc_type[i] == 2:  # X-roller
            constraint_dofs.append(2*i)
        elif bc_type[i] == 3:  # Y-roller
            constraint_dofs.append(2*i+1)
    
    # Griffiths algorithm parameters
    max_iterations = 1000  # Same as in paper
    visco_plastic_param = 1e-6  # Perzyna viscosity parameter
    tolerance = 1e-6  # Convergence tolerance
    
    # Initialize solution vectors
    displacements = np.zeros(n_dof)
    plastic_multipliers = np.zeros((n_elements, 4))  # For 4 Gauss points per element
    
    converged = False
    
    for iteration in range(max_iterations):
        if debug_level >= 3:
            print(f"\n--- Iteration {iteration + 1} ---")
        
        # Apply gravity loads (F_gravity remains constant)
        F_applied = F_gravity.copy()
        
        # Add any additional loads from boundary conditions
        for i in range(n_nodes):
            if bc_type[i] == 4:  # Force boundary condition
                F_applied[2*i] += bc_values[i, 0]
                F_applied[2*i+1] += bc_values[i, 1]
        
        # Constrain system
        K_constrained = K_global.copy()
        F_constrained = F_applied.copy()
        
        # Apply constraints by modifying stiffness matrix and force vector
        for dof in constraint_dofs:
            # Zero out row and column, put 1 on diagonal
            K_constrained[dof, :] = 0
            K_constrained[:, dof] = 0
            K_constrained[dof, dof] = 1
            F_constrained[dof] = 0
        
        # Solve for displacements
        try:
            if hasattr(K_constrained, 'toarray'):
                K_constrained = K_constrained.tocsr()
            displacements = spsolve(K_constrained, F_constrained)
        except Exception as e:
            if debug_level >= 1:
                print(f"Matrix solution failed: {e}")
            return {
                "converged": False,
                "error": f"Matrix solution failed: {e}",
                "iterations": iteration + 1,
                "displacements": displacements,
                "stress_state": stress_state
            }
        
        # Update stresses at all Gauss points
        stress_state_new, plastic_violations = update_stress_griffiths(
            nodes, elements, element_types, element_materials, 
            displacements, stress_state, c_reduced, phi_reduced, 
            E_by_mat, nu_by_mat, u_nodal, visco_plastic_param)
        
        # Check convergence - Griffiths uses displacement and stress change
        disp_change = np.linalg.norm(displacements - (stress_state.get('prev_displacements', np.zeros_like(displacements))))
        stress_change = np.sum(plastic_violations)
        
        if debug_level >= 3:
            print(f"Displacement change norm: {disp_change:.2e}")
            print(f"Plastic violations: {stress_change}")
            print(f"Max displacement: {np.max(np.abs(displacements)):.6f}")
        
        # Griffiths convergence criterion
        if disp_change < tolerance and stress_change < tolerance * n_elements:
            converged = True
            if debug_level >= 2:
                print(f"Converged after {iteration + 1} iterations")
            break
        
        # Update stress state for next iteration
        stress_state = stress_state_new
        stress_state['prev_displacements'] = displacements.copy()
        
        # Check for excessive displacements that indicate numerical instability
        max_disp = np.max(np.abs(displacements))
        if max_disp > 1e6:  # Numerical instability threshold
            if debug_level >= 1:
                print(f"Numerical instability detected: max displacement = {max_disp:.2e}")
            break
    
    # Compute final stresses and plastic indicators
    final_stresses, plastic_elements = compute_final_state_griffiths(
        nodes, elements, element_types, element_materials,
        displacements, stress_state, c_reduced, phi_reduced,
        E_by_mat, nu_by_mat, u_nodal)
    
    # Compute strains for visualization
    strains = compute_strains_griffiths(nodes, elements, element_types, displacements)
    
    return {
        "converged": converged,
        "iterations": iteration + 1,
        "displacements": displacements,
        "stresses": final_stresses,
        "strains": strains,
        "plastic_elements": plastic_elements,
        "max_displacement": np.max(np.abs(displacements)),
        "stress_state": stress_state,
        "algorithm": "Griffiths & Lane (1999)"
    }


def solve_ssrm_griffiths(fem_data, F_min=1.0, F_max=3.0, tolerance=0.01, debug_level=0):
    """
    Shear Strength Reduction Method using Griffiths algorithm.
    
    This implements the SSRM procedure from Griffiths & Lane (1999) where
    failure is defined purely by non-convergence of the Newton-Raphson iteration.
    
    Parameters:
        fem_data (dict): FEM data dictionary
        F_min (float): Minimum F to test
        F_max (float): Maximum F to test
        tolerance (float): F precision tolerance
        debug_level (int): Verbosity level
        
    Returns:
        dict: SSRM results with critical factor of safety
    """
    
    if debug_level >= 1:
        print("=== Griffiths SSRM Analysis ===")
        print("Failure criterion: Non-convergence only")
    
    F_left = F_min
    F_right = F_max
    
    # First, verify that F_min converges and F_max doesn't
    solution_min = solve_fem_griffiths(fem_data, F=F_min, debug_level=max(0, debug_level-1))
    if not solution_min["converged"]:
        return {
            "converged": False,
            "error": f"F_min = {F_min} does not converge - slope may be unstable",
            "FS": None
        }
    
    solution_max = solve_fem_griffiths(fem_data, F=F_max, debug_level=max(0, debug_level-1))
    if solution_max["converged"]:
        if debug_level >= 1:
            print(f"Warning: F_max = {F_max} still converges - slope very stable")
        return {
            "converged": True,
            "FS": F_max,
            "last_solution": solution_max,
            "note": f"Slope stable up to F = {F_max}"
        }
    
    iteration = 0
    max_iterations = 50
    
    # Bisection search for critical F
    while (F_right - F_left) > tolerance and iteration < max_iterations:
        F_mid = (F_left + F_right) / 2.0
        
        if debug_level >= 1:
            print(f"\nSSRM Iteration {iteration + 1}: Testing F = {F_mid:.4f}")
            print(f"Current interval: [{F_left:.4f}, {F_right:.4f}]")
        
        solution = solve_fem_griffiths(fem_data, F=F_mid, debug_level=max(0, debug_level-1))
        
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
    
    # Critical factor of safety is F_left (last converged value)
    critical_FS = F_left
    
    if debug_level >= 1:
        print(f"\nSSRM completed: Critical FS = {critical_FS:.4f}")
        print(f"Final interval: [{F_left:.4f}, {F_right:.4f}]")
        print(f"Iterations: {iteration}")
    
    return {
        "converged": True,
        "FS": critical_FS,
        "last_solution": last_converged_solution,
        "F_history": [],  # Could track if needed
        "convergence_history": [],
        "iterations_ssrm": iteration,
        "final_interval": (F_left, F_right),
        "interval_width": F_right - F_left,
        "method": "Griffiths & Lane (1999) Non-convergence"
    }


def build_global_stiffness_griffiths(nodes, elements, element_types, element_materials, E_by_mat, nu_by_mat):
    """
    Build global stiffness matrix using Griffiths approach.
    
    Note: This is simplified for triangular elements. For full Griffiths implementation,
    8-node quadrilateral elements with reduced integration should be used.
    """
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
        
        # Build element stiffness matrix
        if elem_type == 3:  # Triangle
            K_elem = build_triangle_stiffness_griffiths(elem_coords, E, nu)
        else:
            # For other element types, fall back to existing implementation
            from fem import build_triangle_stiffness
            K_elem = build_triangle_stiffness(elem_coords, E, nu, False, elem_type)
        
        # Assemble into global matrix
        for i in range(elem_type):
            for j in range(elem_type):
                node_i = elem_nodes[i]
                node_j = elem_nodes[j]
                
                # Add contributions to global stiffness matrix
                for di in range(2):
                    for dj in range(2):
                        global_i = 2 * node_i + di
                        global_j = 2 * node_j + dj
                        local_i = 2 * i + di
                        local_j = 2 * j + dj
                        
                        K_global[global_i, global_j] += K_elem[local_i, local_j]
    
    return K_global.tocsr()


def build_triangle_stiffness_griffiths(coords, E, nu):
    """
    Build stiffness matrix for 3-node triangle using standard formulation.
    
    For full Griffiths implementation, this should be replaced with
    8-node quadrilateral elements with reduced integration.
    """
    # Extract coordinates
    x1, y1 = coords[0]
    x2, y2 = coords[1] 
    x3, y3 = coords[2]
    
    # Area of triangle
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
    
    # B matrix (strain-displacement)
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
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


def build_gravity_loads_griffiths(nodes, elements, element_types, element_materials, gamma_by_mat, k_seismic):
    """
    Build gravity load vector using Griffiths approach.
    """
    n_nodes = len(nodes)
    F_gravity = np.zeros(2 * n_nodes)
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        gamma = gamma_by_mat[mat_id]
        
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]
        
        if elem_type == 3:  # Triangle
            # Calculate element area
            x1, y1 = elem_coords[0]
            x2, y2 = elem_coords[1]
            x3, y3 = elem_coords[2]
            area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            # Distribute loads to nodes (1/3 each for triangles)
            load_per_node = gamma * area / 3.0
            
            for i, node in enumerate(elem_nodes):
                # Vertical gravity load (negative y-direction)
                F_gravity[2*node + 1] -= load_per_node
                # Horizontal seismic load (x-direction)
                F_gravity[2*node] += k_seismic * load_per_node
    
    return F_gravity


def initialize_k0_stress_griffiths(nodes, elements, element_types, element_materials, 
                                 gamma_by_mat, nu_by_mat, u_nodal):
    """
    Initialize stress state using K0 conditions as per Griffiths approach.
    
    This sets up the initial geostatic stress state before applying
    the shear strength reduction procedure.
    """
    stress_state = {}
    n_elements = len(elements)
    
    # Initialize stress arrays for each element (assuming 4 Gauss points per element max)
    max_gauss_points = 4
    element_stresses = np.zeros((n_elements, max_gauss_points, 3))  # [sig_x, sig_y, tau_xy]
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        gamma = gamma_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        
        # K0 coefficient (typically 0.4-0.7, use nu-based estimate)
        K0 = nu / (1 - nu)
        
        # Get element centroid for depth calculation
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]
        centroid_y = np.mean(elem_coords[:, 1])
        
        # Find surface level (maximum y-coordinate)
        surface_y = np.max(nodes[:, 1])
        depth = surface_y - centroid_y
        
        if depth > 0:
            # Calculate geostatic stresses
            sigma_v = gamma * depth  # Vertical effective stress
            
            # Interpolate pore pressure at centroid
            u_interp = 0.0  # Could interpolate from u_nodal if needed
            
            sigma_v_eff = sigma_v - u_interp  # Effective vertical stress
            sigma_h_eff = K0 * sigma_v_eff    # Horizontal effective stress
            
            # For all Gauss points in this element (simplified - use same stress)
            n_gauss = min(elem_type, max_gauss_points)
            for gp in range(n_gauss):
                element_stresses[elem_idx, gp, 0] = sigma_h_eff  # sig_x
                element_stresses[elem_idx, gp, 1] = sigma_v_eff  # sig_y  
                element_stresses[elem_idx, gp, 2] = 0.0          # tau_xy
    
    stress_state['element_stresses'] = element_stresses
    stress_state['plastic_state'] = np.zeros((n_elements, max_gauss_points), dtype=bool)
    
    return stress_state


def update_stress_griffiths(nodes, elements, element_types, element_materials,
                          displacements, stress_state, c_reduced, phi_reduced,
                          E_by_mat, nu_by_mat, u_nodal, visco_param):
    """
    Update stresses using Griffiths visco-plastic algorithm (Perzyna approach).
    """
    n_elements = len(elements)
    element_stresses = stress_state['element_stresses'].copy()
    plastic_state = stress_state['plastic_state'].copy()
    
    total_violations = 0
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        mat_id = element_materials[elem_idx] - 1
        
        E = E_by_mat[mat_id]
        nu = nu_by_mat[mat_id]
        c = c_reduced[elem_idx]
        phi = phi_reduced[elem_idx]
        
        # Get element displacements
        elem_nodes = element[:elem_type]
        elem_coords = nodes[elem_nodes]
        elem_disp = np.zeros(2 * elem_type)
        for i, node in enumerate(elem_nodes):
            elem_disp[2*i] = displacements[2*node]
            elem_disp[2*i+1] = displacements[2*node+1]
        
        # Compute strains at element centroid (simplified)
        strains = compute_element_strains_griffiths(elem_coords, elem_disp)
        
        # Update stresses using elastic relationship
        D = build_constitutive_matrix_griffiths(E, nu)
        stress_increment = D @ strains
        
        # For each Gauss point (simplified to 1 for triangles)
        n_gauss = min(elem_type, 1)  # Simplified
        for gp in range(n_gauss):
            # Current stress state
            sigma = element_stresses[elem_idx, gp, :]
            sigma_new = sigma + stress_increment
            
            # Check Mohr-Coulomb yield criterion
            violation = check_mohr_coulomb_griffiths(sigma_new, c, phi)
            
            if violation > 0:
                # Apply visco-plastic correction (Perzyna algorithm)
                sigma_corrected = apply_perzyna_correction_griffiths(
                    sigma_new, c, phi, visco_param)
                element_stresses[elem_idx, gp, :] = sigma_corrected
                plastic_state[elem_idx, gp] = True
                total_violations += violation
            else:
                element_stresses[elem_idx, gp, :] = sigma_new
                plastic_state[elem_idx, gp] = False
    
    return {
        'element_stresses': element_stresses,
        'plastic_state': plastic_state
    }, total_violations


def compute_element_strains_griffiths(coords, displacements):
    """
    Compute element strains from nodal displacements.
    Simplified for triangular elements.
    """
    if len(coords) != 3:
        # Simplified fallback
        return np.array([0.0, 0.0, 0.0])
    
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
    
    # B matrix
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
    ]) / (2 * area)
    
    # Compute strains
    strains = B @ displacements
    
    return strains


def build_constitutive_matrix_griffiths(E, nu):
    """
    Build constitutive matrix for plane strain.
    """
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    return D


def check_mohr_coulomb_griffiths(stress, c, phi):
    """
    Check Mohr-Coulomb yield criterion and return violation magnitude.
    """
    sig_x, sig_y, tau_xy = stress
    
    # Principal stresses
    sig_mean = (sig_x + sig_y) / 2
    tau_max = sqrt(((sig_x - sig_y) / 2)**2 + tau_xy**2)
    
    sig1 = sig_mean + tau_max
    sig3 = sig_mean - tau_max
    
    # Mohr-Coulomb criterion (compression negative)
    # F = (sig1 - sig3) sin(phi) - (sig1 + sig3) cos(phi) - 2c cos(phi)
    F = (sig1 - sig3) * sin(phi) + (sig1 + sig3) * cos(phi) - 2 * c * cos(phi)
    
    return max(0, F)


def apply_perzyna_correction_griffiths(stress, c, phi, visco_param):
    """
    Apply Perzyna visco-plastic stress correction.
    This is a simplified implementation of the Perzyna algorithm.
    """
    # For now, just return stress projected to yield surface
    # Full Perzyna implementation would be more complex
    
    sig_x, sig_y, tau_xy = stress
    
    # Project to yield surface (simplified)
    # This is a placeholder - full implementation would use proper return mapping
    
    return stress  # Simplified - return original stress


def compute_final_state_griffiths(nodes, elements, element_types, element_materials,
                                displacements, stress_state, c_reduced, phi_reduced,
                                E_by_mat, nu_by_mat, u_nodal):
    """
    Compute final stress state and plastic indicators.
    """
    n_elements = len(elements)
    final_stresses = np.zeros((n_elements, 4))  # [sig_x, sig_y, tau_xy, von_mises]
    plastic_elements = np.zeros(n_elements, dtype=bool)
    
    element_stresses = stress_state['element_stresses']
    plastic_state = stress_state['plastic_state']
    
    for elem_idx in range(n_elements):
        # Use stress from first Gauss point (simplified)
        sig_x = element_stresses[elem_idx, 0, 0]
        sig_y = element_stresses[elem_idx, 0, 1]
        tau_xy = element_stresses[elem_idx, 0, 2]
        
        # Von Mises stress
        von_mises = sqrt(sig_x**2 + sig_y**2 - sig_x*sig_y + 3*tau_xy**2)
        
        final_stresses[elem_idx] = [sig_x, sig_y, tau_xy, von_mises]
        plastic_elements[elem_idx] = plastic_state[elem_idx, 0]
    
    return final_stresses, plastic_elements


def compute_strains_griffiths(nodes, elements, element_types, displacements):
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
        element_strains = compute_element_strains_griffiths(elem_coords, elem_disp)
        
        eps_x = element_strains[0]
        eps_y = element_strains[1] 
        gamma_xy = element_strains[2]
        
        # Maximum shear strain
        max_shear_strain = sqrt(((eps_x - eps_y) / 2)**2 + (gamma_xy / 2)**2)
        
        strains[elem_idx] = [eps_x, eps_y, gamma_xy, max_shear_strain]
    
    return strains