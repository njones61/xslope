
from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
from plot import plot_inputs, plot_mesh, plot_polygons
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")

#print_dictionary(slope_data)

#plot_inputs(slope_data)

polygons = build_polygons(slope_data)

#plot_polygons(polygons)

target_size = 2

mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')

#plot_mesh(mesh, materials=slope_data['materials'])

fem_data = build_fem_data(slope_data, mesh)

#print_dictionary(fem_data)

#plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True, material_table=True, label_elements=False, label_nodes=False)


#solution = solve_fem(fem_data, F=1.0, debug_level=2)

# Test single FEM analysis first
# print("Testing single FEM analysis with F=1.0...")
# solution = solve_fem(fem_data, F=1.0, debug_level=3)

# if solution.get("converged", False):
#     print("✓ FEM analysis converged")
#     print(f"Max displacement: {np.max(np.abs(solution.get('displacements', [0]))):.6f}")
#     print(f"Plastic elements: {np.sum(solution.get('plastic_elements', []))}/{len(solution.get('plastic_elements', []))}")
#     print(f"Plastic elements array: {solution.get('plastic_elements', [])[:20]}...")  # Show first 20 elements
# else:
#     print("✗ FEM analysis failed to converge")
#     print(f"Error: {solution.get('error', 'Unknown error')}")
#     print(f"Iterations: {solution.get('iterations', 'Unknown')}")

solution = solve_ssrm(fem_data, debug_level=2)

# Plot multiple result types
print("\n=== FEM Results Visualization ===")
plot_fem_results(fem_data, solution, plot_type='stress, deformation')
