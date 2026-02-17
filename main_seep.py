import numpy as np
from pathlib import Path

from xslope.fileio import load_slope_data, print_dictionary
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json
from xslope.plot import plot_inputs, plot_mesh, plot_polygons, plot_polygons_separately
from xslope.plot_seep import plot_seep_data, plot_seep_solution
from xslope.seep import build_seep_data, run_seepage_analysis, save_seep_data_to_json, export_seep_solution

input_file = "docs/seep/files/xslope_johnson_res.xlsx"
input_path = Path(input_file)

slope_data = load_slope_data(input_file)

plot_inputs(slope_data, figsize=(12, 6), mode='seep', mat_table=False, tab_loc='top', save_png=True)

element_type = 'tri6'
re_mesh = True

# Use existing mesh from slope_data if available, otherwise build a new one
if slope_data.get("mesh") is not None and not re_mesh:
    print("Using existing mesh file.")
    mesh = slope_data["mesh"]
else:
    print("No existing mesh found in slope_data, building new mesh from profile line data.")
    polygons = build_polygons(slope_data, debug=True)

    # find the x-range of the ground_surface and use it to set the target size
    x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
    target_size = (x_range[1] - x_range[0]) / 80

    mesh = build_mesh_from_polygons(polygons, target_size, element_type)
    mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_file)

seep_data = build_seep_data(mesh, slope_data)

plot_seep_data(seep_data, figsize=(12, 6), show_nodes=True, show_bc=True, label_elements=False, label_nodes=False)

solution = run_seepage_analysis(seep_data, tol=1e-4)

plot_seep_solution(seep_data, solution, figsize=(12, 6), variable="head", vectors=False, flowlines=True, 
                          mesh=False, levels=20, base_mat=3, fill_contours=False, phreatic=True, save_png=True)

# Save seep solution to CSV
seep_file = input_path.parent / f"{input_path.stem}_seep.csv"
export_seep_solution(seep_data, solution, seep_file)