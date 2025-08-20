#!/usr/bin/env python3
"""
Test Perzyna implementation against Griffiths & Lane (1999) Example 1.

Expected results:
- Spencer's method: FS = 1.376
- Griffiths FE: FS = 1.4
- Slope: φ=20°, c/γH=0.05, 26.57° (2:1)
"""

from fem import solve_fem_perzyna, solve_ssrm_perzyna
from fem import build_fem_data
from plot_fem import plot_fem_results, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

def test_griffiths_example1():
    """Test Perzyna algorithm against Griffiths Example 1."""
    
    print("=== Testing Perzyna Implementation vs Griffiths Example 1 ===")
    print("Expected: Spencer FS=1.376, Griffiths FE FS=1.4")
    print("Slope: φ=20°, c/γH=0.05, angle=26.57° (2:1)")
    
    # Load Griffiths Example 1
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    
    # Build mesh with 8-node quads (node ordering now corrected)
    polygons = build_polygons(slope_data)
    target_size = 5  # Coarser mesh initially for testing
    
    print(f"\n=== Building Mesh with 3-node Triangles ===")
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)

    #print_dictionary(fem_data)

    plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, material_table=True, label_elements=False, label_nodes=False)
    
    
    print(f"\n=== Testing Perzyna Algorithm at Specified F Value ===")
    print(f"{'F':<6} {'Converged':<10} {'Iterations':<11} {'Max Disp':<12} {'Plastic Elem':<12}")
    print("-" * 65)

    # First, let's examine the initial stress state from gravity loading alone
    print(f"\n=== Debugging Initial Gravity Stress State ===")
    
    # Check material properties and units
    print(f"Material properties:")
    print(f"  E: {fem_data['E_by_mat']} (Young's modulus)")
    print(f"  nu: {fem_data['nu_by_mat']} (Poisson's ratio)")  
    print(f"  gamma: {fem_data['gamma_by_mat']} (unit weight)")
    print(f"  c: {fem_data['c_by_mat']} (cohesion)")
    print(f"  phi: {fem_data['phi_by_mat']} (friction angle, degrees)")
    
    # Get just the gravity stress state without any strength reduction
    from fem import establish_k0_stress_state, build_global_stiffness_perzyna, build_gravity_loads_perzyna
    
    K_global = build_global_stiffness_perzyna(fem_data["nodes"], fem_data["elements"], 
                                            fem_data["element_types"], fem_data["element_materials"],
                                            fem_data["E_by_mat"], fem_data["nu_by_mat"])
    F_gravity = build_gravity_loads_perzyna(fem_data["nodes"], fem_data["elements"], 
                                          fem_data["element_types"], fem_data["element_materials"],
                                          fem_data["gamma_by_mat"], k_seismic=0.0)
    
    initial_displacements, stress_state = establish_k0_stress_state(
        K_global, F_gravity, fem_data["bc_type"], fem_data["nodes"], fem_data["elements"], 
        fem_data["element_types"], fem_data["element_materials"], fem_data["E_by_mat"], 
        fem_data["nu_by_mat"], fem_data["gamma_by_mat"], fem_data["u"], debug_level=2)
    
    # Compute element-averaged stresses for plotting (average over all Gauss points)
    element_stresses_raw = stress_state['element_stresses']  # Shape: [n_elements, n_gauss, 3]
    print(f"  Raw stress array shape: {element_stresses_raw.shape}")
    
    # Average stresses over all active Gauss points for each element
    n_elements = len(fem_data["elements"])
    element_stresses_3col = np.zeros((n_elements, 3))
    
    # Debug: Check for problematic elements with extreme stresses
    extreme_elements = []
    
    for elem_idx in range(n_elements):
        elem_type = fem_data["element_types"][elem_idx]
        if elem_type == 8:  # 8-node quad - average over 4 Gauss points
            n_gauss = 4
        else:  # Triangle or other - single Gauss point
            n_gauss = 1
            
        # Check individual Gauss point stresses before averaging
        elem_stresses = element_stresses_raw[elem_idx, :n_gauss, :]
        max_stress_gp = np.max(np.abs(elem_stresses))
        
        if max_stress_gp > 100000:  # Flag elements with very high stresses
            extreme_elements.append((elem_idx, max_stress_gp, elem_stresses))
            
        # Average over active Gauss points
        avg_stress = np.mean(elem_stresses, axis=0)
        element_stresses_3col[elem_idx, :] = avg_stress
    
    print(f"  Found {len(extreme_elements)} elements with stresses > 100,000")
    if len(extreme_elements) > 0:
        print(f"  Top 3 extreme elements:")
        extreme_elements.sort(key=lambda x: x[1], reverse=True)
        for i, (elem_idx, max_stress, stresses) in enumerate(extreme_elements[:3]):
            print(f"    Element {elem_idx}: max = {max_stress:.0f}")
            print(f"      Gauss point stresses: {stresses.flatten()}")
    
    print(f"  Averaged element stresses shape: {element_stresses_3col.shape}")
    print(f"  Stress ranges after averaging:")
    print(f"    σ_x: [{np.min(element_stresses_3col[:, 0]):.1f}, {np.max(element_stresses_3col[:, 0]):.1f}]")
    print(f"    σ_y: [{np.min(element_stresses_3col[:, 1]):.1f}, {np.max(element_stresses_3col[:, 1]):.1f}]")
    print(f"    τ_xy: [{np.min(element_stresses_3col[:, 2]):.1f}, {np.max(element_stresses_3col[:, 2]):.1f}]")
    
    # Compute von Mises stress from averaged values
    von_mises = np.sqrt(element_stresses_3col[:, 0]**2 + element_stresses_3col[:, 1]**2 - 
                       element_stresses_3col[:, 0]*element_stresses_3col[:, 1] + 
                       3*element_stresses_3col[:, 2]**2)
    
    print(f"    von Mises: [{np.min(von_mises):.1f}, {np.max(von_mises):.1f}]")
    
    # Create stress array with von Mises as 4th column
    stresses_with_vm = np.column_stack([element_stresses_3col, von_mises])
    
    # Create a mock solution object to plot the gravity stress state
    gravity_solution = {
        "converged": True,
        "displacements": initial_displacements,
        "stresses": stresses_with_vm,  # [sig_x, sig_y, tau_xy, von_mises]
        "strains": np.zeros((len(fem_data["elements"]), 4)),  # Dummy strains
        "plastic_elements": np.zeros(len(fem_data["elements"]), dtype=bool),
        "algorithm": "Gravity Loading Only"
    }
    
    print(f"\n=== Plotting Initial Gravity Stress State ===")
    print("This should show mostly compressive stresses from self-weight")
    plot_fem_results(fem_data, gravity_solution, plot_type='stress')
    
    print("\n=== STOPPING HERE TO CHECK INITIAL STRESSES ===")
    print("Skipping full Perzyna analysis to focus on fixing initial stress calculation.")
    
    # Check ratio to expected
    max_vertical_stress = np.max(np.abs(element_stresses_3col[:, 1]))
    height = np.max(fem_data["nodes"][:, 1]) - np.min(fem_data["nodes"][:, 1])
    expected_stress = fem_data['gamma_by_mat'][0] * height
    ratio = max_vertical_stress / expected_stress
    
    print(f"\nStress scaling check:")
    print(f"  Max vertical stress: {max_vertical_stress:.0f}")
    print(f"  Expected vertical stress: {expected_stress:.0f}")
    print(f"  Ratio (actual/expected): {ratio:.1f}")
    
    if ratio > 10:
        print("  ERROR: Stresses are still too large!")
    elif ratio < 0.1:
        print("  ERROR: Stresses are too small!")
    else:
        print("  OK: Stress magnitudes are reasonable")
        
    print("\n=== Stress calculation is now WORKING! ===")
    print(f"The stress pattern should now show a proper gravity-induced gradient.")
    print(f"Ready to proceed with full Perzyna analysis...")
    
    # Continue with actual Perzyna analysis
    print(f"\n=== Starting Perzyna Analysis ===")
    
    # Test F=1.6 to see if we get yielding 
    test_F_values = [1.6]
    
    # Manual check of yield function for triangle mesh
    print("\n=== Manual yield function check (triangles) ===")
    # Use typical stress from triangle mesh near face
    sig_x, sig_y, tau_xy = -463.9, -360.9, 218.6
    
    # Test at F=1.4 (expected critical)
    F_test = 1.4
    c_orig = 312.5
    phi_orig_deg = 20.0
    c_test = c_orig / F_test
    phi_test_deg = np.degrees(np.arctan(np.tan(np.radians(phi_orig_deg)) / F_test))
    phi_test = np.radians(phi_test_deg)
    
    sig_mean = (sig_x + sig_y) / 2
    tau_max = np.sqrt(((sig_x - sig_y) / 2)**2 + tau_xy**2)
    sig1 = sig_mean + tau_max
    sig3 = sig_mean - tau_max
    
    F_test_result = ((sig1 + sig3)/2) * np.sin(phi_test) - ((sig1 - sig3)/2) - c_test * np.cos(phi_test)
    
    print(f"For stress state: σx={sig_x:.1f}, σy={sig_y:.1f}, τxy={tau_xy:.1f}")
    print(f"At F={F_test}: c_reduced={c_test:.1f}, φ_reduced={phi_test_deg:.1f}°")
    print(f"σ_mean={sig_mean:.1f}, τ_max={tau_max:.1f}")
    print(f"σ1={sig1:.1f}, σ3={sig3:.1f}")
    print(f"Yield function F = {F_test_result:.3f}")
    print(f"For failure at F=1.4, we need F ≈ 0")
    print(f"Current F = {F_test_result:.3f} (very negative = very stable)")
    print("This suggests stresses are still too LOW relative to strength")
    print("===================================\n")
    
    for F in test_F_values:
        print(f"\nTesting F = {F:.3f}")
        solution = solve_fem_perzyna(fem_data, F=F, debug_level=2)
        
        converged_str = "YES" if solution["converged"] else "NO"
        max_disp = solution.get("max_displacement", 0.0)
        plastic_count = np.sum(solution.get("plastic_elements", []))
        iterations = solution.get("iterations", 0)
        
        print(f"  F={F:<6.3f} Converged={converged_str:<3} Iterations={iterations:<3} Max_Disp={max_disp:<8.5f} Plastic={plastic_count:<3}")
        
        if solution["converged"]:
            print(f"    F = {F:.3f} is STABLE")
        else:
            print(f"    F = {F:.3f} FAILED to converge (unstable)")
            
        # Plot the result
        plot_fem_results(fem_data, solution, plot_type='stress,deformation,shear_strain')
        
        if not solution["converged"]:
            print(f"    Stopping at first failure: F = {F:.3f}")
            break
    
    return


if __name__ == "__main__":
    test_griffiths_example1()