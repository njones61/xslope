#!/usr/bin/env python3
"""
Debug boundary conditions for Griffiths Example 1.
"""

from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

def debug_boundary_conditions():
    """Debug the boundary condition application."""
    
    print("=== Debugging Boundary Conditions ===")
    
    # Load case
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    
    # Check ground surface
    ground_surface = slope_data['ground_surface']
    coords = list(ground_surface.coords)
    print(f"Ground surface coordinates: {coords}")
    
    x_left = coords[0][0]  # Should be 0.0
    x_right = coords[-1][0]  # Should be 160.0
    print(f"Expected x_left: {x_left}, x_right: {x_right}")
    
    # Build a simple mesh
    polygons = build_polygons(slope_data)
    target_size = 20  # Coarse for debugging
    mesh = build_mesh_from_polygons(polygons, target_size, 'quad8')
    
    nodes = mesh['nodes']
    print(f"\nMesh has {len(nodes)} nodes")
    print(f"X coordinate range: [{np.min(nodes[:, 0]):.3f}, {np.max(nodes[:, 0]):.3f}]")
    print(f"Y coordinate range: [{np.min(nodes[:, 1]):.3f}, {np.max(nodes[:, 1]):.3f}]")
    
    # Check for nodes at boundaries
    tolerance = 1e-6
    left_nodes = np.abs(nodes[:, 0] - x_left) < tolerance
    right_nodes = np.abs(nodes[:, 0] - x_right) < tolerance
    
    print(f"\nNodes at x_left ({x_left}):")
    left_indices = np.where(left_nodes)[0]
    print(f"  Count: {len(left_indices)}")
    if len(left_indices) > 0:
        for i in left_indices[:5]:  # Show first 5
            print(f"    Node {i}: ({nodes[i, 0]:.6f}, {nodes[i, 1]:.6f})")
    
    print(f"\nNodes at x_right ({x_right}):")
    right_indices = np.where(right_nodes)[0]
    print(f"  Count: {len(right_indices)}")
    if len(right_indices) > 0:
        for i in right_indices[:5]:  # Show first 5
            print(f"    Node {i}: ({nodes[i, 0]:.6f}, {nodes[i, 1]:.6f})")
    
    # If no exact matches, find closest nodes
    if len(left_indices) == 0:
        closest_left = np.argmin(np.abs(nodes[:, 0] - x_left))
        print(f"\nClosest to x_left: Node {closest_left} at ({nodes[closest_left, 0]:.6f}, {nodes[closest_left, 1]:.6f})")
        print(f"  Distance: {abs(nodes[closest_left, 0] - x_left):.6f}")
    
    if len(right_indices) == 0:
        closest_right = np.argmin(np.abs(nodes[:, 0] - x_right))
        print(f"\nClosest to x_right: Node {closest_right} at ({nodes[closest_right, 0]:.6f}, {nodes[closest_right, 1]:.6f})")
        print(f"  Distance: {abs(nodes[closest_right, 0] - x_right):.6f}")
    
    # Now test the actual boundary condition application
    print(f"\n=== Testing Boundary Condition Application ===")
    fem_data = build_fem_data(slope_data, mesh)
    
    bc_type = fem_data['bc_type']
    bc_summary = np.bincount(bc_type, minlength=5)
    print(f"\nBoundary condition summary:")
    print(f"  Type 0 (free): {bc_summary[0]} nodes")
    print(f"  Type 1 (fixed): {bc_summary[1]} nodes") 
    print(f"  Type 2 (x-roller): {bc_summary[2]} nodes")
    print(f"  Type 3 (y-roller): {bc_summary[3]} nodes")
    print(f"  Type 4 (force): {bc_summary[4]} nodes")

if __name__ == "__main__":
    debug_boundary_conditions()