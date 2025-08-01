from fileio import load_slope_data

from mesh import build_polygons, build_mesh_from_polygons
from mesh import export_mesh_to_json, import_mesh_from_json
from plot import plot_inputs, plot_mesh, plot_polygons, plot_polygons_separately
from seep import build_seep_data, run_seepage_analysis, save_seep_data_to_json, export_seep_solution
from plot_seep import plot_seep_data, plot_seep_solution
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_lface5.xlsx")

plot_inputs(slope_data)

polygons = build_polygons(slope_data)

# plot_polygons_separately(polygons)

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 100

target_size = 5

# Build quadrilateral mesh
mesh = build_mesh_from_polygons(polygons, target_size, 'tri6')

plot_mesh(mesh, materials=slope_data['materials'])

seep_data = build_seep_data(mesh, slope_data)

plot_seep_data(seep_data, show_nodes=True, show_bc=True, material_table=True, label_elements=False, label_nodes=False)

solution = run_seepage_analysis(seep_data)

plot_seep_solution(seep_data, solution, levels=30, base_mat=2, fill_contours=False, phreatic=True)

export_mesh_to_json(mesh, "inputs/slope/seep_mesh_lface5.json")
export_seep_solution(seep_data, solution, "inputs/slope/seep_solution_lface5.csv")