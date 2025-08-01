from fileio import load_slope_data

from mesh import build_polygons, build_mesh_from_polygons, get_quad_mesh_presets
from mesh import export_mesh_to_json, import_mesh_from_json, test_1d_element_alignment
from mesh import add_intersection_points_to_polygons, extract_reinforcement_line_geometry
from plot import plot_inputs, plot_polygons, plot_polygons_separately, plot_mesh
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_reinf5.xlsx")

# Extract reinforcement lines in the correct format
test_lines = extract_reinforcement_line_geometry(slope_data)
print(f"Extracted {len(test_lines)} reinforcement lines")

# Build polygons with reinforcement line intersections
polygons = build_polygons(slope_data, reinf_lines=test_lines, debug=True)

# plot_polygons_separately(polygons)

# Use fixed target size that was working before
target_size = 2

print(f"\nGenerating triangular mesh with {len(test_lines)} reinforcement lines...")
print(f"Target mesh size: {target_size:.2f}")

# Build triangular mesh with 1D reinforcement elements
mesh_tri = build_mesh_from_polygons(
    polygons, 
    target_size, 
    element_type='tri3', 
    lines=test_lines,
    debug=True
)

print(f"Generated mesh with {len(mesh_tri['nodes'])} nodes")
if 'elements_1d' in mesh_tri:
    print(f"Generated {len(mesh_tri['elements_1d'])} 1D elements")
    
    # Test reinforcement line alignment
    test_success = test_1d_element_alignment(mesh_tri, test_lines, debug=True)
    print(f"1D element alignment test: {'PASSED' if test_success else 'FAILED'}")
else:
    print("No 1D elements were generated")

# Export and plot the final mesh
export_mesh_to_json(mesh_tri, "mesh_tri.json")
plot_mesh(mesh_tri, materials=slope_data['materials'])