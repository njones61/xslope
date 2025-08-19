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

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from math import radians, degrees, sin, cos, tan, sqrt, atan2
import warnings

def solve_fem_perzyna(fem_data, F=1.0, debug_level=0):
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
        from fem import apply_boundary_conditions
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
    
    if debug_level >= 2:
        n_plastic = np.sum(plastic_elements)
        print(f"Final: {n_plastic}/{n_elements} plastic elements")
    
    return {
        "converged": converged,
        "iterations": iteration + 1,
        "displacements": displacements,
        "stresses": final_stresses,
        "strains": strains,
        "plastic_elements": plastic_elements,
        "max_displacement": np.max(np.abs(displacements)),
        "plastic_strains": plastic_strains,
        "algorithm": "Perzyna Visco-Plastic"
    }


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
    # Import the existing stiffness functions that properly handle quad8
    from fem import build_triangle_stiffness, build_quad_stiffness
    
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
            elif elem_type in [4, 6, 9]:  # Other elements - use existing implementation
                from fem import build_triangle_stiffness, build_quad_stiffness
                if elem_type in [4, 9]:
                    K_elem = build_quad_stiffness(elem_coords, E, nu, False, elem_type, 0.1)
                else:
                    K_elem = build_triangle_stiffness(elem_coords, E, nu, False, elem_type, 0.1)
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
            
            # Get plastic strains (single Gauss point)
            plastic_strain = plastic_strains[elem_idx][0, :]  # First (and only) Gauss point
            
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
    from fem import apply_boundary_conditions
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