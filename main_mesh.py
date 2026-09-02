import os
import tempfile

from xslope.fileio import load_slope_data
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         export_mesh_to_json, test_1d_element_alignment,
                         extract_reinforcement_line_geometry)
from xslope.plot import plot_mesh

# A reinforced model: the mesh has to place 1D elements along every soil-nail
# line, which is what the alignment test below checks.
slope_data = load_slope_data("docs/inputs/slope/xslope_reinf.xlsx")

# Extract reinforcement lines in the format the mesher wants
test_lines = extract_reinforcement_line_geometry(slope_data)
print(f"Extracted {len(test_lines)} reinforcement lines")

# Build polygons with reinforcement line intersections
polygons = get_material_polygons(slope_data, reinf_lines=test_lines)

target_size = 5

print("Building a tri6 mesh with reinforcement lines:")
mesh = build_mesh_from_polygons(
    polygons,
    target_size,
    element_type='tri6',
    lines=test_lines,
    debug=False
)
print(f"Generated mesh with {len(mesh['nodes'])} nodes")

if 'elements_1d' in mesh:
    print(f"Generated {len(mesh['elements_1d'])} 1D elements")
    test_success = test_1d_element_alignment(mesh, test_lines, debug=True)
    print(f"1D element alignment test: {'PASSED' if test_success else 'FAILED'}")
else:
    print("No 1D elements were generated")

# Export and plot the final mesh
mesh_path = os.path.join(tempfile.gettempdir(), "xslope_reinf_mesh.json")
export_mesh_to_json(mesh, mesh_path)
print(f"Mesh written to {mesh_path}")
plot_mesh(mesh, materials=slope_data['materials'])
