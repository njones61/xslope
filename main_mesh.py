from fileio import load_globals

from mesh import build_polygons, build_mesh_with_regions, plot_mesh_with_materials, plot_polygons, plot_polygons_separately, reorder_mesh
from plot import plot_inputs
import numpy as np

data = load_globals("inputs/slopes/input_template_lface3.xlsx")

# plot_inputs(data)

polygons = build_polygons(data['profile_lines'], max_depth=data['max_depth'])

# plot_polygons_separately(polygons)

# build a list of region ids
region_ids = [i for i in range(len(polygons))]

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in data['ground_surface'].coords), max(x for x, _ in data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 150

# Build triangular mesh
print("Building triangular mesh...")
nodes_tri, elements_tri, mat_ids_tri = build_mesh_with_regions(polygons, region_ids, target_size, 'tri')

# reorder the mesh
nodes_tri, elements_tri = reorder_mesh(nodes_tri, elements_tri)

plot_mesh_with_materials(nodes_tri, elements_tri, mat_ids_tri, materials=data['materials'])

# Build quadrilateral mesh
print("\nBuilding quadrilateral mesh...")
nodes_quad, elements_quad, mat_ids_quad = build_mesh_with_regions(polygons, region_ids, target_size, 'quad')


# reorder the mesh
nodes_quad, elements_quad = reorder_mesh(nodes_quad, elements_quad)

plot_mesh_with_materials(nodes_quad, elements_quad, mat_ids_quad, materials=data['materials'])
