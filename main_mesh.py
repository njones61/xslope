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
target_size = 5

# Test the new linear-first approach with quadratic elements
print("Testing  with reinforcement lines:")
mesh = build_mesh_from_polygons(
    polygons, 
    target_size, 
    element_type='tri6', 
    lines=None,
    debug=True
)
print(f"Generated mesh with {len(mesh['nodes'])} nodes")
if 'elements_1d' in mesh:
    print(f"Generated {len(mesh['elements_1d'])} 1D elements")
print()


print(f"Generated mesh with {len(mesh['nodes'])} nodes")
if 'elements_1d' in mesh:
    print(f"Generated {len(mesh['elements_1d'])} 1D elements")
    
    # Test reinforcement line alignment
    test_success = test_1d_element_alignment(mesh, test_lines, debug=True)
    print(f"1D element alignment test: {'PASSED' if test_success else 'FAILED'}")
else:
    print("No 1D elements were generated")

# Export and plot the final mesh
export_mesh_to_json(mesh, "mesh.json")
plot_mesh(mesh, materials=slope_data['materials'])