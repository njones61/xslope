
from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
from plot import plot_inputs, plot_mesh, plot_polygons
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")

#print_dictionary(slope_data)

#plot_inputs(slope_data)

polygons = build_polygons(slope_data)

#plot_polygons(polygons)

target_size = 8

mesh = build_mesh_from_polygons(polygons, target_size, 'quad8')

#plot_mesh(mesh, materials=slope_data['materials'])

fem_data = build_fem_data(slope_data, mesh)

#print_dictionary(fem_data)

plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, material_table=True, label_elements=False, label_nodes=False)

print("=== Testing New Improved FEM Solver ===")

# Test 1: Single FEM analysis with new two-phase solver
print("\n1. Testing new two-phase elastic-plastic solver...")
solution_new = solve_fem(fem_data, F=1.28, debug_level=2)

if solution_new.get("converged", False):
    print("✓ New solver converged successfully!")
    print(f"  Phase: {solution_new.get('phase', 'unknown')}")
    print(f"  Total iterations: {solution_new.get('iterations', 'Unknown')}")
    if 'elastic_iterations' in solution_new:
        print(f"  Elastic iterations: {solution_new['elastic_iterations']}")
        print(f"  Plastic iterations: {solution_new['plastic_iterations']}")
    print(f"  Max displacement: {np.max(np.abs(solution_new.get('displacements', [0]))):.6f}")
    print(f"  Plastic elements: {np.sum(solution_new.get('plastic_elements', []))}/{len(solution_new.get('plastic_elements', []))}")
else:
    print("✗ New solver failed to converge")
    print(f"  Error: {solution_new.get('error', 'Unknown error')}")

# Test 2: Compare with old solver (if it still works)
print("\n2. Testing old solver for comparison...")
try:
    solution_old = solve_fem(fem_data, F=1.28, debug_level=1)
    if solution_old.get("converged", False):
        print("✓ Old solver also converged")
        print(f"  Iterations: {solution_old.get('iterations', 'Unknown')}")
        print(f"  Max displacement: {np.max(np.abs(solution_old.get('displacements', [0]))):.6f}")
        print(f"  Plastic elements: {np.sum(solution_old.get('plastic_elements', []))}/{len(solution_old.get('plastic_elements', []))}")
    else:
        print("✗ Old solver failed to converge")
        print(f"  Error: {solution_old.get('error', 'Unknown error')}")
except Exception as e:
    print(f"✗ Old solver encountered an error: {e}")

# Test 3: Test new SSRM with bisection method
print("\n3. Testing new SSRM with bisection method...")
solution_ssrm_new = solve_ssrm(fem_data, F_min=1.0, F_max=5.0, tolerance=0.01, debug_level=1)

if solution_ssrm_new.get("converged", False):
    print("✓ SSRM bisection completed successfully!")
    print(f"  Factor of Safety: {solution_ssrm_new.get('FS', 'Unknown'):.4f}")
    print(f"  Method: {solution_ssrm_new.get('method', 'Unknown')}")
    print(f"  Total SSRM iterations: {solution_ssrm_new.get('iterations_ssrm', 'Unknown')}")
    if 'final_interval' in solution_ssrm_new:
        F_left, F_right = solution_ssrm_new['final_interval']
        print(f"  Final interval: [{F_left:.4f}, {F_right:.4f}]")
        print(f"  Interval width: {solution_ssrm_new.get('interval_width', 'Unknown'):.4f}")
else:
    print("✗ SSRM bisection failed")
    print(f"  Error: {solution_ssrm_new.get('error', 'Unknown error')}")

# Plot results from new solver
if solution_new.get("converged", False):
    print("\n=== FEM Results Visualization (New Solver) ===")
    plot_fem_results(fem_data, solution_new, plot_type='stress, shear_strain, deformation')
else:
    print("\nCannot plot results - solver did not converge")
