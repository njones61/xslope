#!/usr/bin/env python3
"""
Test script to verify that generate_slices handles multiple intersections correctly.
"""

import numpy as np
from shapely.geometry import LineString
import slice

def create_test_data_with_multiple_intersections():
    """Create test data where a circle intersects the ground surface multiple times."""
    
    # Create a ground surface with a crest and flat areas
    x_ground = np.linspace(0, 100, 100)
    y_ground = np.where(x_ground < 30, 40,  # Left flat area
                np.where(x_ground < 70, 60 + 10 * np.sin((x_ground - 30) * 0.1),  # Crest area
                40))  # Right flat area
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
    
    # Circle that intersects multiple times
    circle = {
        'Xo': 50.0,
        'Yo': 70.0,  # Circle center above the crest
        'Depth': 30.0,
        'R': 35.0    # Large radius to intersect multiple times
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

def test_multiple_intersections():
    """Test that generate_slices handles multiple intersections correctly."""
    
    print("=== TESTING MULTIPLE INTERSECTIONS ===\n")
    
    data, circle = create_test_data_with_multiple_intersections()
    
    # Check intersection points directly
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    intersection_points = slice.circle_polyline_intersections(Xo, Yo, R, data['ground_surface'])
    
    print(f"Found {len(intersection_points)} intersection points:")
    for i, pt in enumerate(intersection_points):
        print(f"  Point {i+1}: ({pt.x:.3f}, {pt.y:.3f})")
    
    print()
    
    # Test generate_slices
    print("Testing generate_slices with multiple intersections...")
    success, result = slice.generate_slices(data, circle=circle, num_slices=20, debug=False)
    
    if success:
        df, clipped_surface = result
        print(f"✅ Successfully generated {len(df)} slices")
        print(f"✅ Failure surface type: {type(clipped_surface)}")
        
        # Check that we're using the correct intersection points
        print(f"✅ First slice x_l: {df.iloc[0]['x_l']:.3f}")
        print(f"✅ Last slice x_r: {df.iloc[-1]['x_r']:.3f}")
        
        # Verify that the intersection points are reasonable
        if len(intersection_points) >= 2:
            sorted_points = sorted(intersection_points, key=lambda p: p.x)
            y_first, y_last = sorted_points[0].y, sorted_points[-1].y
            
            if y_first > y_last:
                # Right-facing slope: should use first two points
                expected_left = sorted_points[0].x
                expected_right = sorted_points[1].x
                print(f"✅ Right-facing slope detected (y_first={y_first:.3f} > y_last={y_last:.3f})")
            else:
                # Left-facing slope: should use last two points
                expected_left = sorted_points[-2].x
                expected_right = sorted_points[-1].x
                print(f"✅ Left-facing slope detected (y_first={y_first:.3f} <= y_last={y_last:.3f})")
            
            print(f"✅ Expected x range: {expected_left:.3f} to {expected_right:.3f}")
            print(f"✅ Actual x range: {df.iloc[0]['x_l']:.3f} to {df.iloc[-1]['x_r']:.3f}")
            
            # Check if the actual range matches the expected range (within tolerance)
            tol = 1e-3
            if (abs(df.iloc[0]['x_l'] - expected_left) < tol and 
                abs(df.iloc[-1]['x_r'] - expected_right) < tol):
                print("✅ Intersection points correctly selected!")
            else:
                print("❌ Intersection points not correctly selected!")
        else:
            print("❌ Not enough intersection points to test pruning logic")
    else:
        print(f"❌ Failed: {result}")

if __name__ == "__main__":
    test_multiple_intersections() 