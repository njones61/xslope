#!/usr/bin/env python3
"""
Debug the displacement solution to understand why strains/stresses are wrong.
"""

import numpy as np
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons  
from fem import build_fem_data
from fem_perzyna import build_global_stiffness_perzyna, build_gravity_loads_perzyna
from fem import apply_boundary_conditions
from scipy.sparse.linalg import spsolve

def debug_gravity_solution():
    """Debug the basic gravity loading solution."""
    print("=== Debugging Gravity Loading Solution ===")
    
    # Load the same data as the failing test
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    mesh = build_mesh_from_polygons(polygons, 5, 'quad8')
    fem_data = build_fem_data(slope_data, mesh)
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"] 
    element_materials = fem_data["element_materials"]
    bc_type = fem_data["bc_type"]
    E_by_mat = fem_data["E_by_mat"]
    nu_by_mat = fem_data["nu_by_mat"]
    gamma_by_mat = fem_data["gamma_by_mat"]
    
    print(f"Mesh info:")
    print(f"  Nodes: {len(nodes)}")
    print(f"  Elements: {len(elements)}")
    print(f"  Element types: {np.unique(element_types)}")
    
    print(f"Material properties:")
    print(f"  E = {E_by_mat[0]:,.0f}")
    print(f"  nu = {nu_by_mat[0]}")
    print(f"  gamma = {gamma_by_mat[0]}")
    
    print(f"Mesh extent:")
    print(f"  X: [{np.min(nodes[:, 0]):.1f}, {np.max(nodes[:, 0]):.1f}]")
    print(f"  Y: [{np.min(nodes[:, 1]):.1f}, {np.max(nodes[:, 1]):.1f}]")
    
    # Expected vertical stress at bottom of slope
    height = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    expected_stress = gamma_by_mat[0] * height
    print(f"  Height = {height:.1f}")
    print(f"  Expected stress at bottom = {expected_stress:.0f}")
    
    print(f"Boundary conditions:")
    bc_counts = np.bincount(bc_type)
    for i, count in enumerate(bc_counts):
        if count > 0:
            print(f"  Type {i}: {count} nodes")
    
    # Build stiffness and loads
    print(f"\\nBuilding system matrices...")
    K_global = build_global_stiffness_perzyna(nodes, elements, element_types, 
                                             element_materials, E_by_mat, nu_by_mat)
    F_gravity = build_gravity_loads_perzyna(nodes, elements, element_types, 
                                           element_materials, gamma_by_mat, k_seismic=0.0)
    
    print(f"Stiffness matrix: {K_global.shape}, nnz = {K_global.nnz}")
    print(f"Force vector: {F_gravity.shape}, max = {np.max(np.abs(F_gravity)):.1f}")
    
    # Apply boundary conditions and solve
    print(f"\\nApplying boundary conditions and solving...")
    K_constrained, F_constrained, constraint_dofs = apply_boundary_conditions(K_global, F_gravity, bc_type, nodes)
    
    print(f"Constrained system: {K_constrained.shape}, constrained DOFs = {len(constraint_dofs)}")
    
    # Solve
    u_free = spsolve(K_constrained.tocsr(), F_constrained)
    
    # Reconstruct full displacement vector
    displacements = np.zeros(2 * len(nodes))
    free_dof_idx = 0
    for i in range(2 * len(nodes)):
        if i not in constraint_dofs:
            displacements[i] = u_free[free_dof_idx]
            free_dof_idx += 1
    
    print(f"\\nDisplacement solution:")
    print(f"  Max |u_x| = {np.max(np.abs(displacements[0::2])):.6f}")
    print(f"  Max |u_y| = {np.max(np.abs(displacements[1::2])):.6f}")
    print(f"  Max total = {np.max(np.abs(displacements)):.6f}")
    
    # Check displacement reasonableness
    # For typical slopes, displacements should be small (mm to cm scale)
    # With mesh ~100 units, displacements should be ~0.001 to 0.1
    max_disp = np.max(np.abs(displacements))
    mesh_size = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
    disp_ratio = max_disp / mesh_size
    
    print(f"  Displacement/mesh ratio = {disp_ratio:.6f}")
    if disp_ratio > 0.01:
        print(f"  WARNING: Displacements seem too large relative to mesh size")
    if disp_ratio < 1e-6:
        print(f"  WARNING: Displacements seem too small")
    
    # Check for problematic nodes
    disp_magnitudes = np.sqrt(displacements[0::2]**2 + displacements[1::2]**2)
    max_disp_node = np.argmax(disp_magnitudes)
    
    print(f"  Max displacement at node {max_disp_node}:")
    print(f"    Location: ({nodes[max_disp_node, 0]:.1f}, {nodes[max_disp_node, 1]:.1f})")
    print(f"    Displacement: ({displacements[2*max_disp_node]:.6f}, {displacements[2*max_disp_node+1]:.6f})")
    print(f"    BC type: {bc_type[max_disp_node]}")
    
    return displacements, nodes, elements, element_types, fem_data

def debug_simple_element(displacements, nodes, elements, element_types, fem_data):
    """Debug strain calculation for a single element."""
    print(f"\\n=== Debugging Single Element ===")
    
    # Pick first 8-node quad
    quad8_indices = np.where(element_types == 8)[0]
    if len(quad8_indices) == 0:
        print("No quad8 elements found")
        return
    
    elem_idx = quad8_indices[0]
    element = elements[elem_idx]
    elem_nodes = element[:8]
    elem_coords = nodes[elem_nodes]
    
    print(f"Element {elem_idx} (8-node quad):")
    print(f"  Nodes: {elem_nodes}")
    print(f"  Coordinates:")
    for i, (x, y) in enumerate(elem_coords):
        print(f"    Node {elem_nodes[i]}: ({x:.2f}, {y:.2f})")
    
    # Element displacements
    elem_disp = np.zeros(16)  # 8 nodes x 2 DOF
    for i, node in enumerate(elem_nodes):
        elem_disp[2*i] = displacements[2*node]
        elem_disp[2*i+1] = displacements[2*node+1]
    
    print(f"  Element displacements:")
    for i in range(8):
        print(f"    Node {elem_nodes[i]}: ({elem_disp[2*i]:.8f}, {elem_disp[2*i+1]:.8f})")
    
    max_elem_disp = np.max(np.abs(elem_disp))
    print(f"  Max element displacement: {max_elem_disp:.8f}")
    
    # Compute strains manually with debug output
    from fem_perzyna import compute_quad8_strains_at_xi_eta
    
    # Test at element center
    xi, eta = 0.0, 0.0
    try:
        strains = compute_quad8_strains_at_xi_eta(elem_coords, elem_disp, xi, eta)
        print(f"  Strains at center (xi=0, eta=0): [{strains[0]:.8f}, {strains[1]:.8f}, {strains[2]:.8f}]")
        print(f"  Strain magnitudes: [{abs(strains[0]):.8f}, {abs(strains[1]):.8f}, {abs(strains[2]):.8f}]")
        
        # Check if strains are reasonable
        # Typical engineering strains should be < 0.01 (1%) for small deformations
        max_strain = np.max(np.abs(strains))
        print(f"  Max strain magnitude: {max_strain:.8f} = {max_strain*100:.6f}%")
        
        if max_strain > 0.1:
            print(f"  ERROR: Strains are too large (>{max_strain*100:.2f}%)")
        elif max_strain < 1e-8:
            print(f"  WARNING: Strains are too small (<1e-6%)")
        else:
            print(f"  Strains appear reasonable")
            
        # Compute corresponding stress
        E = fem_data["E_by_mat"][0]
        nu = fem_data["nu_by_mat"][0] 
        from fem_perzyna import build_constitutive_matrix_perzyna
        D = build_constitutive_matrix_perzyna(E, nu)
        stresses = D @ strains
        
        print(f"  Resulting stresses: [{stresses[0]:.1f}, {stresses[1]:.1f}, {stresses[2]:.1f}]")
        
        # Expected stress for comparison
        gamma = fem_data["gamma_by_mat"][0]
        elem_y = np.mean(elem_coords[:, 1])
        max_y = np.max(nodes[:, 1])
        depth = max_y - elem_y
        expected_vertical_stress = gamma * depth if depth > 0 else 0
        
        print(f"  Element depth: {depth:.1f}")
        print(f"  Expected vertical stress: {expected_vertical_stress:.1f}")
        
        stress_ratio = abs(stresses[1]) / max(expected_vertical_stress, 1.0)
        print(f"  Actual/Expected stress ratio: {stress_ratio:.1f}")
        
    except Exception as e:
        print(f"  ERROR computing strains: {e}")

if __name__ == "__main__":
    displacements, nodes, elements, element_types, fem_data = debug_gravity_solution()
    debug_simple_element(displacements, nodes, elements, element_types, fem_data)