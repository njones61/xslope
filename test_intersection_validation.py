#!/usr/bin/env python3
"""
Test script to verify intersection validation in generate_slices function.
"""

import numpy as np
from shapely.geometry import LineString
import slice

def create_test_data_with_circle_above_ground():
    """Create test data where the circle is above the ground surface (no intersection)."""
    
    # Ground surface that's below the circle
    x_ground = np.linspace(0, 100, 50)
    y_ground = 30 + 5 * np.sin(x_ground * 0.1)  # Low ground surface
    ground_surface = LineString(list(zip(x_ground, y_ground)))
    
    # Profile lines
    y_profile = y_ground - 5
    profile_lines = [list(zip(x_ground, y_profile))]
    
    # Materials
    materials = [{
        'gamma': 120.0,
        'option': 'mc',
        'c': 500.0,
        'phi': 30.0,
        'd': 0.0,
        'psi': 0.0
    }]
    
    # Circle that's above the ground surface (no intersection)
    circle = {
        'Xo': 50.0,
        'Yo': 60.0,  # Circle center is above ground
        'Depth': 20.0,
        'R': 10.0    # Small radius so it doesn't reach ground
    }
    
    # Data dictionary
    data = {
        'profile_lines': profile_lines,
        'ground_surface': ground_surface,
        'materials': materials,
        'piezo_line': [],
        'piezo_line2': [],
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'dloads': [],
        'dloads2': [],
        'max_depth': 50.0,
        'reinforce_lines': []
    }
    
    return data, circle

def create_test_data_with_circle_tangent():
    """Create test data where the circle is tangent to the ground surface (1 intersection)."""
    
    # Ground surface
    x_ground = np.linspace(0, 100, 50)
    y_ground = 50 + 10 * np.sin(x_ground * 0.1)
    ground_surface = LineString(list(zip(x_ground, y_ground)))
    
    # Profile lines
    y_profile = y_ground - 5
    profile_lines = [list(zip(x_ground, y_profile))]
    
    # Materials
    materials = [{
        'gamma': 120.0,
        'option': 'mc',
        'c': 500.0,
        'phi': 30.0,
        'd': 0.0,
        'psi': 0.0
    }]
    
    # Circle that's tangent to the ground surface
    circle = {
        'Xo': 50.0,
        'Yo': 65.0,  # Circle center positioned so it's tangent
        'Depth': 20.0,
        'R': 15.0    # Radius that makes it tangent
    }
    
    # Data dictionary
    data = {
        'profile_lines': profile_lines,
        'ground_surface': ground_surface,
        'materials': materials,
        'piezo_line': [],
        'piezo_line2': [],
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'dloads': [],
        'dloads2': [],
        'max_depth': 50.0,
        'reinforce_lines': []
    }
    
    return data, circle

def test_intersection_validation():
    """Test that generate_slices properly validates intersection points."""
    
    print("=== TESTING INTERSECTION VALIDATION ===\n")
    
    # Test 1: Circle above ground (no intersection)
    print("Test 1: Circle above ground surface (should fail)")
    data, circle = create_test_data_with_circle_above_ground()
    
    # Check intersection points directly
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    intersection_points = slice.circle_polyline_intersections(Xo, Yo, R, data['ground_surface'])
    print(f"  Found {len(intersection_points)} intersection points")
    
    # Test generate_slices
    success, result = slice.generate_slices(data, circle=circle, num_slices=20, debug=False)
    if success:
        print("  ❌ ERROR: Should have failed but succeeded")
    else:
        print(f"  ✅ Correctly failed: {result}")
    
    print()
    
    # Test 2: Circle tangent to ground (1 intersection)
    print("Test 2: Circle tangent to ground surface (should fail)")
    data, circle = create_test_data_with_circle_tangent()
    
    # Check intersection points directly
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    intersection_points = slice.circle_polyline_intersections(Xo, Yo, R, data['ground_surface'])
    print(f"  Found {len(intersection_points)} intersection points")
    
    # Test generate_slices
    success, result = slice.generate_slices(data, circle=circle, num_slices=20, debug=False)
    if success:
        print("  ❌ ERROR: Should have failed but succeeded")
    else:
        print(f"  ✅ Correctly failed: {result}")
    
    print()
    
    # Test 3: Valid case (2 intersections)
    print("Test 3: Valid circle intersecting ground surface (should succeed)")
    data, circle = create_test_data_with_circle_above_ground()
    # Modify circle to intersect properly
    circle['R'] = 35.0  # Larger radius to intersect ground
    
    # Check intersection points directly
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    intersection_points = slice.circle_polyline_intersections(Xo, Yo, R, data['ground_surface'])
    print(f"  Found {len(intersection_points)} intersection points")
    
    # Test generate_slices
    success, result = slice.generate_slices(data, circle=circle, num_slices=20, debug=False)
    if success:
        print("  ✅ Correctly succeeded")
        df, clipped_surface = result
        print(f"  Generated {len(df)} slices")
    else:
        print(f"  ❌ ERROR: Should have succeeded but failed: {result}")

if __name__ == "__main__":
    test_intersection_validation() 