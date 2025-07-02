#!/usr/bin/env python3
"""
Debug script to investigate thin slices at the edges.
"""

import numpy as np
from shapely.geometry import LineString
import slice

def create_simple_test_data():
    """Create a simple test case to debug thin slices."""
    
    # Simple ground surface
    x_ground = np.linspace(0, 100, 50)
    y_ground = 50 + 10 * np.sin(x_ground * 0.1)
    ground_surface = LineString(list(zip(x_ground, y_ground)))
    
    # Simple profile line
    y_profile = y_ground - 5
    profile_lines = [list(zip(x_ground, y_profile))]
    
    # Simple material
    materials = [{
        'gamma': 120.0,
        'option': 'mc',
        'c': 500.0,
        'phi': 30.0,
        'd': 0.0,
        'psi': 0.0
    }]
    
    # Circular failure surface
    circle = {
        'Xo': 50.0,
        'Yo': 60.0,
        'Depth': 20.0,
        'R': 25.0
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

def debug_slice_boundaries():
    """Debug the slice boundary generation process."""
    
    print("=== DEBUGGING THIN SLICES ===\n")
    
    data, circle = create_simple_test_data()
    
    # Get the intersection points
    Xo, Yo, R = circle['Xo'], circle['Yo'], circle['R']
    ground_surface = data['ground_surface']
    
    print("Circle parameters:")
    print(f"  Center: ({Xo}, {Yo})")
    print(f"  Radius: {R}")
    print()
    
    # Find intersection points using analytic method
    intersection_points = slice.circle_polyline_intersections(Xo, Yo, R, ground_surface)
    
    print(f"Found {len(intersection_points)} intersection points:")
    for i, pt in enumerate(intersection_points):
        print(f"  Point {i+1}: ({pt.x:.6f}, {pt.y:.6f})")
    
    print()
    
    # Sort by x-coordinate
    sorted_points = sorted(intersection_points, key=lambda p: p.x)
    x_min, x_max = sorted_points[0].x, sorted_points[-1].x
    y_left, y_right = sorted_points[0].y, sorted_points[-1].y
    
    print(f"Boundary points:")
    print(f"  Left:  ({x_min:.6f}, {y_left:.6f})")
    print(f"  Right: ({x_max:.6f}, {y_right:.6f})")
    print(f"  Width: {x_max - x_min:.6f}")
    print()
    
    # Now let's trace through the slice boundary generation step by step
    print("=== TRACING SLICE BOUNDARY GENERATION ===\n")
    
    # Step 1: Start with intersection points
    fixed_xs = set()
    fixed_xs.update([x_min, x_max])
    print(f"Step 1 - Added intersection points: {sorted(fixed_xs)}")
    
    # Step 2: Add profile line points above failure surface
    profile_lines = data["profile_lines"]
    for line_idx, line in enumerate(profile_lines):
        line_coords = np.array(line)
        x_coords = line_coords[:, 0]
        y_coords = line_coords[:, 1]
        
        # Filter points within x-range
        mask = (x_coords >= x_min) & (x_coords <= x_max)
        x_filtered = x_coords[mask]
        y_filtered = y_coords[mask]
        
        if len(x_filtered) == 0:
            continue
            
        # Check if points are above the failure surface
        failure_y = slice.get_circular_y_coordinates(x_filtered, Xo, Yo, R)
        above_mask = y_filtered > failure_y
        
        # Add points that are above the failure surface
        new_points = x_filtered[above_mask]
        fixed_xs.update(new_points)
        if len(new_points) > 0:
            print(f"Step 2 - Added {len(new_points)} profile points from layer {line_idx}: {sorted(new_points)}")
    
    print(f"After Step 2: {sorted(fixed_xs)}")
    
    # Step 3: Add dload points (none in this test)
    dloads = data["dloads"]
    if dloads:
        dload_points = set()
        for line in dloads:
            for pt in line:
                if x_min <= pt['X'] <= x_max:
                    dload_points.add(pt['X'])
        fixed_xs.update(dload_points)
        if dload_points:
            print(f"Step 3 - Added dload points: {sorted(dload_points)}")
    
    # Step 4: Add non_circ points (none in this test)
    non_circ = None
    if non_circ:
        non_circ_points = set()
        for pt in non_circ:
            if x_min <= pt['X'] <= x_max:
                non_circ_points.add(pt['X'])
        fixed_xs.update(non_circ_points)
        if non_circ_points:
            print(f"Step 4 - Added non_circ points: {sorted(non_circ_points)}")
    
    # Step 5: Find intersections with profile lines and failure surface
    print("\nStep 5 - Finding profile/failure surface intersections...")
    for line_idx, line in enumerate(profile_lines):
        line_geom = LineString(line)
        # Create a dense circle representation for intersection
        theta_range = np.linspace(np.pi, 2 * np.pi, 200)
        circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        circle_line = LineString(circle_coords)
        
        intersection = line_geom.intersection(circle_line)
        if not intersection.is_empty:
            intersection_points = []
            if hasattr(intersection, 'x'):
                intersection_points.append(intersection.x)
            elif hasattr(intersection, 'geoms'):
                for geom in intersection.geoms:
                    if hasattr(geom, 'x'):
                        intersection_points.append(geom.x)
            
            if intersection_points:
                fixed_xs.update(intersection_points)
                print(f"  Layer {line_idx}: Added intersection points {intersection_points}")
    
    print(f"After Step 5: {sorted(fixed_xs)}")
    
    # Step 6: Find intersections with piezometric lines (none in this test)
    piezo_line = data["piezo_line"]
    if piezo_line:
        print("\nStep 6 - Finding piezometric/failure surface intersections...")
        piezo_geom = LineString(piezo_line)
        theta_range = np.linspace(np.pi, 2 * np.pi, 200)
        circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        circle_line = LineString(circle_coords)
        intersection = piezo_geom.intersection(circle_line)
        if not intersection.is_empty:
            intersection_points = []
            if hasattr(intersection, 'x'):
                if x_min <= intersection.x <= x_max:
                    intersection_points.append(intersection.x)
            elif hasattr(intersection, 'geoms'):
                for geom in intersection.geoms:
                    if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                        intersection_points.append(geom.x)
            
            if intersection_points:
                fixed_xs.update(intersection_points)
                print(f"  Added piezometric intersection points: {intersection_points}")
    
    # Step 7: Clean duplicates
    tolerance = 1e-6
    cleaned_xs = []
    for x in sorted(fixed_xs):
        if not cleaned_xs or abs(x - cleaned_xs[-1]) > tolerance:
            cleaned_xs.append(x)
    
    fixed_xs = cleaned_xs
    print(f"\nStep 7 - After cleaning duplicates: {fixed_xs}")
    
    # Step 8: Generate slice boundaries
    print("\nStep 8 - Generating slice boundaries...")
    segment_lengths = [fixed_xs[i + 1] - fixed_xs[i] for i in range(len(fixed_xs) - 1)]
    total_length = sum(segment_lengths)
    print(f"  Segment lengths: {[f'{l:.6f}' for l in segment_lengths]}")
    print(f"  Total length: {total_length:.6f}")
    
    all_xs = [fixed_xs[0]]
    for i in range(len(fixed_xs) - 1):
        x_start = fixed_xs[i]
        x_end = fixed_xs[i + 1]
        segment_length = x_end - x_start
        n_subdiv = max(1, int(round((segment_length / total_length) * 20)))  # num_slices=20
        xs = np.linspace(x_start, x_end, n_subdiv + 1).tolist()
        all_xs.extend(xs[1:])
        print(f"  Segment {i+1}: {x_start:.6f} to {x_end:.6f}, length={segment_length:.6f}, subdivisions={n_subdiv}")
        if n_subdiv > 1:
            print(f"    Added points: {xs[1:]}")
    
    print(f"\nStep 8 - All boundaries before cleaning: {[f'{x:.6f}' for x in all_xs]}")
    
    # Step 9: Remove thin slices
    print("\nStep 9 - Removing thin slices...")
    min_width = 1e-2
    cleaned_xs = [all_xs[0]]
    for x in all_xs[1:-1]:
        if abs(x - cleaned_xs[-1]) >= min_width:
            cleaned_xs.append(x)
        else:
            print(f"  Removed thin slice boundary: {x:.6f} (too close to {cleaned_xs[-1]:.6f})")
    cleaned_xs.append(all_xs[-1])
    all_xs = cleaned_xs
    
    print(f"Step 9 - Final boundaries: {[f'{x:.6f}' for x in all_xs]}")
    
    # Show the final slice widths
    print("\nFinal slice widths:")
    for i in range(len(all_xs) - 1):
        width = all_xs[i+1] - all_xs[i]
        print(f"  Slice {i+1}: {all_xs[i]:.6f} to {all_xs[i+1]:.6f}, width={width:.6f}")
        if width < min_width:
            print(f"    ⚠️  WARNING: Thin slice detected!")

if __name__ == "__main__":
    debug_slice_boundaries() 