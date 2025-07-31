from fileio import load_slope_data

from mesh import build_polygons, build_mesh_from_polygons, get_quad_mesh_presets
from mesh import export_mesh_to_json, import_mesh_from_json, test_1d_element_alignment
from plot import plot_inputs, plot_polygons, plot_polygons_separately, plot_mesh
import numpy as np

slope_data = load_slope_data("inputs/slope/input_template_reinf5.xlsx")

plot_inputs(slope_data)

# Extract reinforcement lines in the correct format
test_lines = []
if 'reinforce_lines' in slope_data and slope_data['reinforce_lines']:
    for line in slope_data['reinforce_lines']:
        # Convert from dict format to tuple format
        line_coords = [(point['X'], point['Y']) for point in line]
        test_lines.append(line_coords)
    print(f"Extracted {len(test_lines)} reinforcement lines")
    for i, line in enumerate(test_lines):
        print(f"Line {i}: {line}")

# Note: Distributed loads are not part of the mesh - they are only for load application
print(f"Total reinforcement lines for meshing: {len(test_lines)}")

# plot_inputs(slope_data)

polygons = build_polygons(slope_data)

print(f"\nGenerating mesh with {len(test_lines)} 1D lines...")
target_size = 10

# Test with just the first reinforcement line to see what happens
test_single_line = [test_lines[0]]  # Just the first line: [(0.0, 0.0), (4.0, 0.0), (16.0, 0.0), (20.0, 0.0)]
print(f"Testing with single line: {test_single_line[0]}")

# Test 1D line meshing in isolation (without polygons)
print("\\nTesting 1D line in isolation...")
mesh_1d_only = build_mesh_from_polygons(
    [],  # No polygons - just the 1D line
    target_size, 
    element_type='tri3', 
    lines=test_single_line,
    debug=True
)

print(f"\\n1D-only mesh: {len(mesh_1d_only.get('nodes', []))} nodes")
if 'elements_1d' in mesh_1d_only:
    print(f"1D-only mesh: {len(mesh_1d_only['elements_1d'])} 1D elements")
    # Show all node coordinates
    print("\\nAll node coordinates:")
    nodes = mesh_1d_only['nodes']
    for i, node in enumerate(nodes):
        print(f"  Node {i}: {node}")
    
    # Show actual coordinates of 1D elements
    print("\\n1D-only element coordinates:")
    elements_1d = mesh_1d_only['elements_1d']
    for i, element in enumerate(elements_1d):
        # Check if element has valid node indices (not zero-padded)
        if len(element) >= 2 and element[1] != 0:  # Skip only if second index is 0 (padding)
            coord1 = nodes[element[0]]
            coord2 = nodes[element[1]]
            print(f"  Element {i}: {coord1} -> {coord2}")
            print(f"    Node indices: {element[0]} -> {element[1]}")
else:
    print("  No 1D elements found")

# Test mesh generation with 1D lines and polygons
print("\\nTesting 1D line with polygons...")
mesh_with_1d = build_mesh_from_polygons(
    polygons, 
    target_size, 
    element_type='tri3', 
    lines=test_single_line,
    debug=True
)

print(f"Generated mesh with {len(mesh_with_1d['nodes'])} nodes")
if 'elements_1d' in mesh_with_1d:
    print(f"Generated {len(mesh_with_1d['elements_1d'])} 1D elements")
    print(f"1D element types: {mesh_with_1d['element_types_1d']}")
    print(f"1D material IDs: {mesh_with_1d['element_materials_1d']}")
    
    # Show actual coordinates of 1D elements
    print("\nActual 1D element coordinates:")
    nodes = mesh_with_1d['nodes']
    elements_1d = mesh_with_1d['elements_1d']
    for i, element in enumerate(elements_1d[:10]):  # Show first 10 elements
        # Check if element has valid node indices (not zero-padded)
        if len(element) >= 2 and element[1] != 0:  # Skip only if second index is 0 (padding)
            coord1 = nodes[element[0]]
            coord2 = nodes[element[1]]
            print(f"  Element {i}: {coord1} -> {coord2}")
    
    print(f"\nExpected line path: {test_single_line[0]}")
    
    # Test 1D element alignment
    test_success = test_1d_element_alignment(mesh_with_1d, test_single_line, debug=True)
    print(f"1D element alignment test: {'PASSED' if test_success else 'FAILED'}")
    
    # Plot the mesh with 1D elements
    plot_mesh(mesh_with_1d, materials=slope_data['materials'], label_elements=True, label_nodes=True)
else:
    print("No 1D elements were generated")

# plot_polygons_separately(polygons)

# build a list of region ids
region_ids = [i for i in range(len(polygons))]

# find the x-range of the ground_surface and use it to set the target size
x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
target_size = (x_range[1] - x_range[0]) / 150

# Check if we have reinforcement lines to test intersection preprocessing
reinforcement_lines = []
if 'reinforce_lines' in slope_data and slope_data['reinforce_lines']:
    for reinforcement in slope_data['reinforce_lines']:
        # Convert from dict format to tuple format
        line_coords = [(point['X'], point['Y']) for point in reinforcement]
        if len(line_coords) >= 2:
            reinforcement_lines.append(line_coords)
    print(f"Found {len(reinforcement_lines)} reinforcement lines for testing")
else:
    print("No reinforcement lines found in slope_data")

# Use only reinforcement lines for mesh generation (distributed loads are not part of mesh)
all_test_lines = reinforcement_lines
print(f"Total reinforcement lines for mesh generation: {len(all_test_lines)}")
    
# Print some slope_data keys to see what's available
print("Available slope_data keys:", list(slope_data.keys()))

target_size = 10

# Test mesh generation with 1D lines if available
if all_test_lines:
    print(f"\nTesting mesh generation with {len(all_test_lines)} 1D lines...")
    mesh_with_lines = build_mesh_from_polygons(
        polygons, 
        target_size, 
        element_type='tri3', 
        lines=all_test_lines,
        debug=True
    )
    print(f"Generated mesh with {len(mesh_with_lines['nodes'])} nodes")
    if 'elements_1d' in mesh_with_lines:
        print(f"Generated {len(mesh_with_lines['elements_1d'])} 1D elements")
        
        # Test all reinforcement lines
        print(f"\nTesting all {len(all_test_lines)} reinforcement lines...")
        
        # Debug: Show all 1D elements found
        print(f"\nAll 1D elements found ({len(mesh_with_lines['elements_1d'])} total):")
        nodes = mesh_with_lines['nodes']
        elements_1d = mesh_with_lines['elements_1d']
        for i, element in enumerate(elements_1d):
            if len(element) >= 2 and element[1] != 0:
                coord1 = nodes[element[0]]
                coord2 = nodes[element[1]]
                print(f"  Element {i}: {coord1} -> {coord2}")
        
        all_lines_success = test_1d_element_alignment(mesh_with_lines, all_test_lines, debug=True)
        print(f"All reinforcement lines alignment test: {'PASSED' if all_lines_success else 'FAILED'}")
        
        plot_mesh(mesh_with_lines, materials=slope_data['materials'], label_elements=True)
    else:
        print("No 1D elements were generated")
else:
    print("No test lines available for intersection testing")


# Build triangular mesh
print("Building triangular mesh...")
mesh_tri = build_mesh_from_polygons(polygons, target_size, element_type='tri3', lines=all_test_lines, debug=True)

export_mesh_to_json(mesh_tri, "mesh_tri.json")

plot_mesh(mesh_tri, materials=slope_data['materials'])

# Build quadrilateral mesh
print("\nBuilding quadrilateral mesh...")

# Try different meshing presets
presets = get_quad_mesh_presets()

mesh_quad = build_mesh_from_polygons(polygons, target_size, element_type='quad4', mesh_params=presets['default'], lines=all_test_lines ,debug=True)

export_mesh_to_json(mesh_quad, "mesh_quad.json")

plot_mesh(mesh_quad, materials=slope_data['materials'])
