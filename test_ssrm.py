#!/usr/bin/env python3

from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_ssrm_convergence
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

# Load data and build FEM model - same as main_fem2.py
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 2, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

print("=== Testing SSRM Implementation ===")
print("Expected factor of safety around 1.28 based on limit equilibrium analysis")

# Run SSRM analysis
ssrm_result = solve_ssrm(fem_data, debug_level=2)

# Print results
if ssrm_result.get("converged", False):
    print(f"\n✓ SSRM converged successfully")
    print(f"Factor of Safety: {ssrm_result['FS']:.4f}")
    print(f"SSRM iterations: {ssrm_result['iterations_ssrm']}")
    
    # Plot final solution with stress contours
    print("\n=== Plotting Final Solution ===")
    if 'last_solution' in ssrm_result:
        plot_fem_results(fem_data, ssrm_result['last_solution'], plot_type='stress')
        
    # Plot SSRM convergence history
    print("\n=== Plotting SSRM Convergence ===")
    plot_ssrm_convergence(ssrm_result)
        
else:
    print(f"\n✗ SSRM failed to converge")
    if 'error' in ssrm_result:
        print(f"Error: {ssrm_result['error']}")
    if 'FS' in ssrm_result and ssrm_result['FS'] is not None:
        print(f"Last converged FS: {ssrm_result['FS']:.4f}")