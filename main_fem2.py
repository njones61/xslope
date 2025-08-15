
from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np
from plot import plot_inputs, plot_mesh, plot_polygons

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

print("=== Testing F=1.0 ===")
solution1 = solve_fem(fem_data, F=1.0, debug_level=1)
print(f"F=1.0: Max displacement = {np.max(np.abs(solution1['displacements'])):.6f}")
print(f"F=1.0: Plastic elements = {np.sum(solution1.get('plastic_elements', []))}")

print("\n=== Testing F=0.05 ===")
solution2 = solve_fem(fem_data, F=0.05, debug_level=1)
print(f"F=0.05: Max displacement = {np.max(np.abs(solution2['displacements'])):.6f}")
print(f"F=0.05: Plastic elements = {np.sum(solution2.get('plastic_elements', []))}")

print(f"\nDisplacement ratio (F=1.0 / F=0.05) = {np.max(np.abs(solution1['displacements'])) / np.max(np.abs(solution2['displacements'])):.2f}")

solution = solution1  # Use F=1.0 for plotting

print(f"Solution converged: {solution['converged']}")
print(f"Max displacement: {np.max(np.abs(solution['displacements']))}")
print(f"Min displacement: {np.min(solution['displacements'])}")
print(f"NaN values: {np.sum(np.isnan(solution['displacements']))}")
print(f"Inf values: {np.sum(np.isinf(solution['displacements']))}")

# Debug: Check if solution dictionaries are actually different
print("\n=== Debugging solution data ===")
print(f"solution1 ID: {id(solution1)}")
print(f"solution2 ID: {id(solution2)}")
print(f"solution1['displacements'] ID: {id(solution1['displacements'])}")
print(f"solution2['displacements'] ID: {id(solution2['displacements'])}")
print(f"Are displacement arrays equal? {np.array_equal(solution1['displacements'], solution2['displacements'])}")
print(f"First 10 displacements F=1.0: {solution1['displacements'][:10]}")
print(f"First 10 displacements F=0.05: {solution2['displacements'][:10]}")

# Plot both solutions to compare
print("\n=== Plotting F=1.0 solution ===")
if np.all(np.isfinite(solution1['displacements'])):
    print(f"Plotting F=1.0: Max displacement = {np.max(np.abs(solution1['displacements'])):.6f}")
    plot_fem_results(fem_data, solution1, plot_type='displacement')
else:
    print("Skipping F=1.0 plot due to non-finite displacement values")

print("\n=== Plotting F=0.05 solution ===")
if np.all(np.isfinite(solution2['displacements'])):
    print(f"Plotting F=0.05: Max displacement = {np.max(np.abs(solution2['displacements'])):.6f}")
    plot_fem_results(fem_data, solution2, plot_type='displacement')
else:
    print("Skipping F=0.05 plot due to non-finite displacement values")