from fileio import load_globals

from mesh import build_polygons, build_mesh_from_polygons, plot_mesh, plot_polygons, plot_polygons_separately
from mesh import save_mesh_to_json, load_mesh_from_json
from plot import plot_inputs
from seep import build_seep_data, run_seepage_analysis, save_seep_data_to_json
from plot_seep import plot_seep_data, plot_seep_solution
import numpy as np

data = load_globals("inputs/slopes/input_template_lface3.xlsx")

# plot_inputs(data)

polygons = build_polygons(data['profile_lines'], max_depth=data['max_depth'])

# plot_polygons_separately(polygons)

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in data['ground_surface'].coords), max(x for x, _ in data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 150

# target_size = 10

# Build quadrilateral mesh
mesh = build_mesh_from_polygons(polygons, target_size, 'tri')

plot_mesh(mesh, materials=data['materials'])

seep_data = build_seep_data(mesh, data)

plot_seep_data(seep_data, show_nodes=True, show_bc=True, material_table=True, label_elements=False)

solution = run_seepage_analysis(seep_data)

plot_seep_solution(seep_data, solution, levels=30, base_mat=2, fill_contours=False, phreatic=True)
