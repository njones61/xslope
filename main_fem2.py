
from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
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


solution = solve_fem(fem_data, F=1.0, debug_level=2)

# Plot multiple result types
print("\n=== FEM Results Visualization ===")
plot_fem_results(fem_data, solution, plot_type='stress, deformation')
