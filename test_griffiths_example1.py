#!/usr/bin/env python3
"""
Test Perzyna implementation against Griffiths & Lane (1999) Example 1.

Expected results:
- Spencer's method: FS = 1.376
- Griffiths FE: FS = 1.4
- Slope: φ=20°, c/γH=0.05, 26.57° (2:1)
"""

from fem_perzyna import solve_fem_perzyna, solve_ssrm_perzyna
from fem import build_fem_data, solve_ssrm
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
    
    # Build mesh with 8-node quadrilaterals (key for Griffiths approach)
    polygons = build_polygons(slope_data)
    target_size = 5  # Coarser mesh initially for testing
    
    print(f"\n=== Building Mesh with 8-node Quadrilaterals ===")
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)

    #print_dictionary(fem_data)

    plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, material_table=True, label_elements=False, label_nodes=False)
    
    print(f"Mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    print(f"Element types: {np.unique(fem_data['element_types'])}")
    
    # Verify material properties
    print(f"\n=== Material Properties ===")
    print(f"c: {fem_data['c_by_mat'][0]:.1f} psf")
    print(f"φ: {fem_data['phi_by_mat'][0]:.1f}°")
    print(f"γ: {fem_data['gamma_by_mat'][0]:.1f} pcf")
    print(f"E: {fem_data['E_by_mat'][0]:.0f} psf")
    print(f"ν: {fem_data['nu_by_mat'][0]:.2f}")
    
    # Check c/γH ratio
    H = 50  # Height from geometry
    c_gamma_H = fem_data['c_by_mat'][0] / (fem_data['gamma_by_mat'][0] * H)
    print(f"c/γH: {c_gamma_H:.3f} (should be 0.05)")
    
    # Test single analysis near critical F (simplified for faster testing)
    test_F_values = [1.0, 1.3]
    
    print(f"\n=== Testing Perzyna Algorithm at Various F Values ===")
    print(f"{'F':<6} {'Converged':<10} {'Iterations':<11} {'Max Disp':<12} {'Plastic Elem':<12}")
    print("-" * 65)
    
    results = []
    
    for F in test_F_values:
        try:
            solution = solve_fem_perzyna(fem_data, F=F, debug_level=1)
            converged = solution.get("converged", False)
            iterations = solution.get("iterations", 0)
            max_disp = solution.get("max_displacement", 0)
            plastic_count = np.sum(solution.get("plastic_elements", []))
            
            print(f"{F:<6.2f} {str(converged):<10} {iterations:<11} {max_disp:<12.6f} {plastic_count:<12}")
            
            results.append((F, converged, iterations, max_disp, plastic_count, solution))
            
            if not converged:
                print(f"\n*** NON-CONVERGENCE at F = {F:.2f} ***")
                print("This indicates slope failure!")
                break
                
        except Exception as e:
            print(f"{F:<6.2f} ERROR: {str(e)[:30]}")
            results.append((F, False, 0, 0, 0, None))
    
    # Skip SSRM for now - focus on testing basic solver
    print(f"\n=== Skipping SSRM Analysis (focus on basic solver first) ===")
    
    # Plot results from individual tests
    critical_solution = None
    critical_F = None
    
    # Find any converged solution from individual tests
    for F, converged, iterations, max_disp, plastic_count, solution in results:
        if converged and solution:
            critical_solution = solution
            critical_F = F
            print(f"\n=== Plotting Individual Test Solution (F = {critical_F:.2f}) ===")
            break
    
    if critical_solution:
        try:
            plot_fem_results(fem_data, critical_solution, plot_type='deformation, stress, shear_strain')
            print("Check plots for rotational failure pattern!")
            print(f"This shows the slope state at critical F = {critical_F:.3f}")
        except Exception as e:
            print(f"Plotting error: {e}")
    else:
        print("\n=== No suitable solution found for plotting ===")
        print("All test cases failed to converge - indicates slope failure")
    
    print(f"\n=== Summary ===")
    print(f"Target: Reproduce Griffiths Example 1 rotational failure")
    print(f"Key test: Does failure surface show circular/rotational pattern?")
    print(f"Expected FS: ~1.4 (Griffiths FE) vs 1.376 (Spencer)")


if __name__ == "__main__":
    test_griffiths_example1()