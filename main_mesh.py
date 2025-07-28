from fileio import load_slope_data

from mesh import build_polygons, build_mesh_from_polygons
from mesh import export_mesh_to_json, import_mesh_from_json
from plot import plot_inputs, plot_polygons, plot_polygons_separately, plot_mesh
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_lface5.xlsx")

# plot_inputs(slope_data)

polygons = build_polygons(slope_data['profile_lines'], max_depth=slope_data['max_depth'])

# plot_polygons_separately(polygons)

# build a list of region ids
region_ids = [i for i in range(len(polygons))]

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 150

target_size = 10

# Build triangular mesh
print("Building triangular mesh...")
mesh_tri = build_mesh_from_polygons(polygons, target_size, 'tri6')

export_mesh_to_json(mesh_tri, "mesh_tri.json")

plot_mesh(mesh_tri, materials=slope_data['materials'])

# Build quadrilateral mesh
print("\nBuilding quadrilateral mesh...")
mesh_quad = build_mesh_from_polygons(polygons, target_size, 'quad8')

export_mesh_to_json(mesh_quad, "mesh_quad.json")

plot_mesh(mesh_quad, materials=slope_data['materials'])
