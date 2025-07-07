from fileio import load_globals

from mesh import build_polygons, build_tri_mesh_with_regions, plot_mesh_with_materials, plot_polygons
from plot import plot_inputs
import numpy as np

data = load_globals("inputs/slopes/input_template_dam2.xlsx")

plot_inputs(data)

polygons = build_polygons(data['profile_lines'], max_depth=data['max_depth'])

plot_polygons(polygons)

# build a list of region ids
region_ids = [i for i in range(len(polygons))]

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in data['ground_surface'].coords), max(x for x, _ in data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 150

nodes, elements, mat_ids = build_tri_mesh_with_regions(polygons, region_ids, target_size)

plot_mesh_with_materials(nodes, elements, mat_ids)
