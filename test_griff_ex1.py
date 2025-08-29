#!/usr/bin/env python3
"""
Test implementation against Griffiths & Lane (1999) Example 1.

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
    
    # Build mesh with 8-node quadrilaterals
    polygons = build_polygons(slope_data)
    target_size = 5  # Coarser mesh initially for testing
    
    print(f"\n=== Building Mesh with 8-node Quads ===")
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)

    # Plot the initial mesh
    plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, 
                  material_table=True, label_elements=False, label_nodes=False)
    
    F_test = 1.0
    
    print(f"\n\n=== Testing Initial Gravity Loading (abort_after=0) ===")
    # First check gravity loading only
    gravity_solution = solve_fem(fem_data, F=F_test, debug_level=2, abort_after=0)
    
    # Plot gravity stress state and yield function
    print(f"\n\n=== Plotting Initial Gravity Stress and Yield Function ===")
    plot_fem_results(fem_data, gravity_solution, 
                     plot_type='stress,yield', 
                     label_elements=True)  # This will show element IDs on stress plot and F values on yield plot
    
    # Check yield function statistics
    yield_fn = gravity_solution.get("yield_function", np.array([]))
    if len(yield_fn) > 0:
        n_yielding = np.sum(yield_fn > 0)
        max_yield = np.max(yield_fn)
        min_yield = np.min(yield_fn)
        print(f"\n\nYield function statistics after gravity loading at F={F_test:.3f}:")
        print(f"  Elements with F>0 (meeting yield criterion): {n_yielding}/{len(yield_fn)}")
        print(f"  Max yield function: {max_yield:.3f}")
        print(f"  Min yield function: {min_yield:.3f}")
        print(f"  Note: F>0 means the element meets the yield criterion but hasn't developed plastic strains yet")
        
        if n_yielding == 0:
            print(f"  WARNING: No elements yielding at F={F_test:.3f} - may need higher F for failure")
    
    print(f"\n\n=== Testing One Perzyna Iteration (abort_after=1) ===")
    # Test one iteration to see plastic correction
    one_iter_solution = solve_fem(fem_data, F=F_test, debug_level=2, abort_after=1)
    
    # Plot after one iteration
    print(f"\n\n=== Plotting After One Perzyna Iteration ===")
    plot_fem_results(fem_data, one_iter_solution, 
                     plot_type='stress,yield,shear_strain', 
                     label_elements=False)
    
    print(f"\n\n=== Testing Full Perzyna Analysis at F = {F_test:.3f} ===")
    solution = solve_fem(fem_data, F=F_test, debug_level=2)
    
    converged_str = "YES" if solution["converged"] else "NO"
    max_disp = solution.get("max_displacement", 0.0)
    plastic_count = np.sum(solution.get("plastic_elements", []))
    yielding_count = np.sum(solution.get("yield_function", []) > 0)
    iterations = solution.get("iterations", 0)
    
    print(f"\n\nResults:")
    print(f"  F={F_test:<6.3f}")
    print(f"  Converged={converged_str}")
    print(f"  Iterations={iterations}")
    print(f"  Max Displacement={max_disp:.5f}")
    print(f"  Elements with F>0 (yielding)={yielding_count}")
    print(f"  Elements with plastic strains={plastic_count}")
    
    if solution["converged"]:
        print(f"  Status: F = {F_test:.3f} is STABLE")
    else:
        print(f"  Status: F = {F_test:.3f} FAILED to converge (unstable)")
    
    # Plot the final state
    print(f"\n\n=== Plotting Final State at F = {F_test:.3f} ===")
    plot_fem_results(fem_data, solution,
                  plot_type='stress,deformation,shear_strain',
                  label_elements=False)
    
    return


if __name__ == "__main__":
    test_griffiths_example1()