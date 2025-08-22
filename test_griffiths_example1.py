#!/usr/bin/env python3
"""
Test Perzyna implementation against Griffiths & Lane (1999) Example 1.

Expected results:
- Spencer's method: FS = 1.376
- Griffiths FE: FS = 1.4
- Slope: φ=20°, c/γH=0.05, 26.57° (2:1)
"""

from fem import solve_fem
from fem import build_fem_data
from plot_fem import plot_fem_results, plot_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

def test_griffiths_example1():
    """Test Perzyna algorithm against Griffiths Example 1."""
    
    print("=== Testing Perzyna Implementation vs Griffiths Example 1 ===")
    print("Expected: Spencer FS=1.376, Griffiths FE FS=1.4")
    print("Slope: φ=20°, c/γH=0.05, angle=26.57° (2:1)")
    
    # Load Griffiths Example 1
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    
    # Build mesh with 3-node triangles
    polygons = build_polygons(slope_data)
    target_size = 5  # Coarser mesh initially for testing
    
    print(f"\n=== Building Mesh with 3-node Triangles ===")
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)

    # Plot the initial mesh
    plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, 
                  material_table=True, label_elements=False, label_nodes=False)
    
    F_test = 1.2
    
    print(f"\n=== Testing Initial Gravity Loading (abort_after=0) ===")
    # First check gravity loading only
    gravity_solution = solve_fem(fem_data, F=F_test, debug_level=2, abort_after=0)
    
    # Plot gravity stress state and yield function
    print(f"\n=== Plotting Initial Gravity Stress and Yield Function ===")
    plot_fem_results(fem_data, gravity_solution, 
                     plot_type='stress,yield', 
                     label_elements=True)  # This will show element IDs on stress plot and F values on yield plot
    
    # Check yield function statistics
    yield_fn = gravity_solution.get("yield_function", np.array([]))
    if len(yield_fn) > 0:
        n_yielding = np.sum(yield_fn > 0)
        max_yield = np.max(yield_fn)
        min_yield = np.min(yield_fn)
        print(f"\nYield function statistics after gravity loading at F={F_test:.3f}:")
        print(f"  Elements yielding: {n_yielding}/{len(yield_fn)}")
        print(f"  Max yield function: {max_yield:.3f}")
        print(f"  Min yield function: {min_yield:.3f}")
        
        if n_yielding == 0:
            print("  WARNING: No elements yielding at F=1.6 - may need higher F for failure")
    
    print(f"\n=== Testing One Perzyna Iteration (abort_after=1) ===")
    # Test one iteration to see plastic correction
    one_iter_solution = solve_fem(fem_data, F=F_test, debug_level=2, abort_after=1)
    
    # Plot after one iteration
    print(f"\n=== Plotting After One Perzyna Iteration ===")
    plot_fem_results(fem_data, one_iter_solution, 
                     plot_type='stress,yield,shear_strain', 
                     label_elements=False)
    
    print(f"\n=== Testing Full Perzyna Analysis at F = {F_test:.3f} ===")
    solution = solve_fem(fem_data, F=F_test, debug_level=2)
    
    converged_str = "YES" if solution["converged"] else "NO"
    max_disp = solution.get("max_displacement", 0.0)
    plastic_count = np.sum(solution.get("plastic_elements", []))
    iterations = solution.get("iterations", 0)
    
    print(f"\nResults:")
    print(f"  F={F_test:<6.3f}")
    print(f"  Converged={converged_str}")
    print(f"  Iterations={iterations}")
    print(f"  Max Displacement={max_disp:.5f}")
    print(f"  Plastic Elements={plastic_count}")
    
    if solution["converged"]:
        print(f"  Status: F = {F_test:.3f} is STABLE")
    else:
        print(f"  Status: F = {F_test:.3f} FAILED to converge (unstable)")
    
    # Plot the final state
    print(f"\n=== Plotting Final State at F = {F_test:.3f} ===")
    plot_fem_results(fem_data, solution, 
                   plot_type='stress,deformation,shear_strain',
                   label_elements=False)
    
    return


if __name__ == "__main__":
    test_griffiths_example1()