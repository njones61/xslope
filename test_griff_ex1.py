#!/usr/bin/env python3
"""
Test implementation against Griffiths & Lane (1999) Example 1.

Expected results:
- Spencer's method: FS = 1.376
- Griffiths FE: FS = 1.4
- Slope: φ=20°, c/γH=0.05, 26.57° (2:1)
"""

import numpy as np

from xslope.fem import solve_fem, build_fem_data
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons
from xslope.plot_fem import plot_fem_results, plot_fem_data

def test_griffiths_example1():
    """Test Perzyna algorithm against Griffiths Example 1."""
    
    print("=== Testing Perzyna Implementation vs Griffiths Example 1 ===")
    print("Expected: Spencer FS=1.376, Griffiths FE FS=1.4")
    print("Slope: φ=20°, c/γH=0.05, angle=26.57° (2:1)")
    
    # Load Griffiths Example 1
    slope_data = load_slope_data("docs/inputs/slope/input_template_griffiths1_6.xlsx")
    
    # Build mesh with 8-node quadrilaterals
    polygons = build_polygons(slope_data)
    target_size = 4  # Coarser mesh initially for testing
    
    print(f"\n=== Building Mesh with 8-node Quads ===")
    mesh = build_mesh_from_polygons(polygons, target_size, 'quad8')
    fem_data = build_fem_data(slope_data, mesh)

    # Plot the initial mesh
    plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, 
                   material_table=True, label_elements=False, label_nodes=False)
    
    F_test = 1.3
    
    
    print(f"\n\n=== Testing Full Perzyna Analysis at F = {F_test:.3f} ===")
    solution = solve_fem(fem_data, F=F_test, debug_level=3)
    
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
                  plot_type=['displace_vector', 'deformation', 'shear_strain'],
                  label_elements=False)
    
    return


if __name__ == "__main__":
    test_griffiths_example1()