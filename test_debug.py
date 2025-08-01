#!/usr/bin/env python3

from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons, extract_reinforcement_line_geometry

def test_element_type(element_type, lines):
    """Test a specific element type with timeout-like behavior"""
    print(f"\n=== Testing {element_type} ===")
    try:
        mesh = build_mesh_from_polygons(
            polygons, 
            2, 
            element_type=element_type, 
            lines=lines,
            debug=False
        )
        print(f"SUCCESS: {element_type} generated mesh with {len(mesh['nodes'])} nodes")
        if 'elements_1d' in mesh:
            print(f"Generated {len(mesh['elements_1d'])} 1D elements")
        return True
    except Exception as e:
        print(f"ERROR: {element_type} failed with: {e}")
        return False

if __name__ == "__main__":
    print("Loading data...")
    slope_data = load_slope_data("inputs/slope/input_template_reinf5.xlsx")
    test_lines = extract_reinforcement_line_geometry(slope_data)
    print(f"Extracted {len(test_lines)} reinforcement lines")
    
    polygons = build_polygons(slope_data, reinf_lines=test_lines, debug=False)
    print("Built polygons with intersection points")
    
    # Test different element types
    element_types = ['tri3', 'quad4', 'tri6', 'quad8']
    
    for element_type in element_types:
        success = test_element_type(element_type, test_lines)
        if not success:
            print(f"STOPPING: {element_type} failed")
            break