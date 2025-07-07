from mesh import build_tri_mesh_with_regions, plot_mesh, plot_mesh_with_materials
import numpy as np

# Define two adjacent rectangles
left_rect = [(0, 0), (1, 0), (1, 1), (0, 1)]
right_rect = [(1, 0), (2, 0), (2, 1), (1, 1)]

polygons = [left_rect, right_rect]
region_ids = [1, 2]  # material IDs
target_size = 0.1

nodes, elements, mat_ids = build_tri_mesh_with_regions(polygons, region_ids, target_size)

plot_mesh_with_materials(nodes, elements, mat_ids)