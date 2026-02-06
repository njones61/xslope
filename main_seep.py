import numpy as np

from xslope.fileio import load_slope_data, print_dictionary
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json, import_mesh_from_json
from xslope.plot import plot_inputs, plot_mesh, plot_polygons, plot_polygons_separately
from xslope.plot_seep import plot_seep_data, plot_seep_solution
from xslope.seep import build_seep_data, run_seepage_analysis, save_seep_data_to_json, export_seep_solution

path = "docs/seep/files/xslope_levee_full.xlsx"

path_prefix = path.replace(".xlsx", "")  # make a copy of the path minus the extension for later use

slope_data = load_slope_data(path)

plot_inputs(slope_data, figsize=(12, 6), mode='seep', mat_table=False, tab_loc='top', save_png=True)

polygons = build_polygons(slope_data, debug=True)

plot_polygons(polygons, figsize=(12, 6), materials=slope_data['materials'], title="Material Zones", nodes=False, legend=False, save_png=True)
# plot_polygons_separately(polygons, materials=slope_data['materials'], save_png=True)

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 100

# Build quadrilateral mesh
mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')

# plot_mesh(mesh, materials=slope_data['materials'])

seep_data = build_seep_data(mesh, slope_data)

# print_dictionary(seep_data)

plot_seep_data(seep_data, figsize=(12, 6), show_nodes=True, show_bc=True, label_elements=False, label_nodes=False)

solution = run_seepage_analysis(seep_data, tol=1e-4)

plot_seep_solution(seep_data, solution, figsize=(12, 6), variable="head", vectors=False, flowlines=True, 
                          mesh=False, levels=20, base_mat=2, fill_contours=False, phreatic=True, save_png=True)

# export_mesh_to_json(mesh, path_prefix + "_mesh.json")
# export_seep_solution(seep_data, solution, path_prefix + "_seep.csv")