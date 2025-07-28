#!/usr/bin/env python3
"""
Test consistency between tri3 and tri6 stiffness matrices.
For the same triangle geometry, the tri6 matrix should reduce to tri3 behavior
when only corner nodes are used.
"""

import numpy as np
from seep import tri3_stiffness_matrix, tri6_stiffness_matrix

def test_triangle_consistency():
    """Test that tri6 reduces to tri3 for corner nodes only."""
    
    # Define a simple triangle
    nodes_tri3 = np.array([
        [0.0, 0.0],  # Node 0
        [1.0, 0.0],  # Node 1  
        [0.0, 1.0]   # Node 2
    ])
    
    # For tri6, add midpoint nodes
    nodes_tri6 = np.array([
        [0.0, 0.0],  # Node 0 (corner)
        [1.0, 0.0],  # Node 1 (corner)
        [0.0, 1.0],  # Node 2 (corner)
        [0.5, 0.0],  # Node 3 (midpoint 0-1)
        [0.5, 0.5],  # Node 4 (midpoint 1-2)
        [0.0, 0.5]   # Node 5 (midpoint 2-0)
    ])
    
    # Isotropic conductivity matrix
    Kmat = np.eye(2)
    
    # Compute stiffness matrices
    ke_tri3 = tri3_stiffness_matrix(nodes_tri3, Kmat)
    ke_tri6 = tri6_stiffness_matrix(nodes_tri6, Kmat)
    
    print("TRI3 stiffness matrix:")
    print(ke_tri3)
    print(f"TRI3 diagonal sum: {np.trace(ke_tri3):.6f}")
    
    print("\nTRI6 stiffness matrix:")
    print(ke_tri6)
    print(f"TRI6 diagonal sum: {np.trace(ke_tri6):.6f}")
    
    # Extract corner-corner part of tri6 matrix
    corner_part = ke_tri6[:3, :3]
    print("\nTRI6 corner-corner part:")
    print(corner_part)
    print(f"Corner diagonal sum: {np.trace(corner_part):.6f}")
    
    # Check if corner part is consistent
    diff = np.abs(ke_tri3 - corner_part)
    max_diff = np.max(diff)
    print(f"\nMax difference between TRI3 and TRI6 corner part: {max_diff:.2e}")
    
    # Test conservation: row sums should be zero for both
    print("\nRow sums (should be ~0 for Laplace operator):")
    print(f"TRI3 row sums: {np.sum(ke_tri3, axis=1)}")
    print(f"TRI6 row sums: {np.sum(ke_tri6, axis=1)}")
    
    return max_diff < 1e-10

def test_clockwise_vs_counterclockwise():
    """Test tri6 behavior with different node orderings."""
    
    # Counterclockwise triangle
    nodes_ccw = np.array([
        [0.0, 0.0],  # Node 0
        [1.0, 0.0],  # Node 1  
        [0.0, 1.0],  # Node 2
        [0.5, 0.0],  # Node 3 (midpoint 0-1)
        [0.5, 0.5],  # Node 4 (midpoint 1-2)
        [0.0, 0.5]   # Node 5 (midpoint 2-0)
    ])
    
    # Clockwise triangle (swap nodes 1 and 2)
    nodes_cw = np.array([
        [0.0, 0.0],  # Node 0
        [0.0, 1.0],  # Node 1 (was 2)
        [1.0, 0.0],  # Node 2 (was 1)
        [0.0, 0.5],  # Node 3 (midpoint 0-1)
        [0.5, 0.5],  # Node 4 (midpoint 1-2)
        [0.5, 0.0]   # Node 5 (midpoint 2-0)
    ])
    
    Kmat = np.eye(2)
    
    ke_ccw = tri6_stiffness_matrix(nodes_ccw, Kmat)
    ke_cw = tri6_stiffness_matrix(nodes_cw, Kmat)
    
    print("\nCounterclockwise triangle diagonal sum:", np.trace(ke_ccw))
    print("Clockwise triangle diagonal sum:", np.trace(ke_cw))
    
    # Both should have same magnitude (differ by sign in det, but abs is used)
    return abs(np.trace(ke_ccw) - np.trace(ke_cw)) < 1e-10

if __name__ == "__main__":
    print("Testing tri3 vs tri6 consistency...")
    consistent = test_triangle_consistency()
    print(f"\nConsistency test: {'PASSED' if consistent else 'FAILED'}")
    
    print("\n" + "="*50)
    print("Testing clockwise vs counterclockwise...")
    orientation_ok = test_clockwise_vs_counterclockwise()
    print(f"\nOrientation test: {'PASSED' if orientation_ok else 'FAILED'}")