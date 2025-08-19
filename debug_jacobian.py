#!/usr/bin/env python3
"""
Debug the Jacobian calculation and coordinate transformation for 8-node quads.

This tests the fundamental coordinate transformation to identify why all 
Gauss points give identical stress values.
"""

import numpy as np
from fem_perzyna import compute_quad8_shape_derivatives

def test_simple_rectangle():
    """Test coordinate transformation on a simple rectangular element."""
    print("=== Testing Simple Rectangle Element ===")
    
    # Define a simple 2x1 rectangular element centered at origin
    # Node ordering: 0,4,1,5,2,6,3,7 (CCW from bottom-left)
    coords = np.array([
        [-1.0, -0.5],  # 0: bottom-left
        [ 1.0, -0.5],  # 1: bottom-right  
        [ 1.0,  0.5],  # 2: top-right
        [-1.0,  0.5],  # 3: top-left
        [ 0.0, -0.5],  # 4: bottom-mid
        [ 1.0,  0.0],  # 5: right-mid
        [ 0.0,  0.5],  # 6: top-mid
        [-1.0,  0.0],  # 7: left-mid
    ])
    
    print(f"Element coordinates:")
    for i, (x, y) in enumerate(coords):
        print(f"  Node {i}: ({x:4.1f}, {y:4.1f})")
    
    # Test at 4 Gauss points for 2x2 integration
    gauss_coords = [
        (-1/np.sqrt(3), -1/np.sqrt(3)),  # GP 0
        ( 1/np.sqrt(3), -1/np.sqrt(3)),  # GP 1  
        ( 1/np.sqrt(3),  1/np.sqrt(3)),  # GP 2
        (-1/np.sqrt(3),  1/np.sqrt(3)),  # GP 3
    ]
    
    print(f"\nTesting coordinate transformation at Gauss points:")
    
    for gp, (xi, eta) in enumerate(gauss_coords):
        print(f"\n--- Gauss Point {gp}: (ξ={xi:6.3f}, η={eta:6.3f}) ---")
        
        # Get shape function derivatives
        dN_dxi, dN_deta = compute_quad8_shape_derivatives(xi, eta)
        
        print(f"Shape function derivatives:")
        print(f"  dN/dξ:  {dN_dxi}")
        print(f"  dN/dη:  {dN_deta}")
        
        # Check derivatives sum to zero (partition of unity)
        sum_dxi = np.sum(dN_dxi)
        sum_deta = np.sum(dN_deta)
        print(f"  Sum checks: dN/dξ={sum_dxi:.6f}, dN/dη={sum_deta:.6f} (should be ~0)")
        
        # Compute Jacobian matrix
        J = np.zeros((2, 2))
        for i in range(8):
            x, y = coords[i]
            J[0, 0] += dN_dxi[i] * x   # dx/dξ
            J[0, 1] += dN_dxi[i] * y   # dy/dξ
            J[1, 0] += dN_deta[i] * x  # dx/dη
            J[1, 1] += dN_deta[i] * y  # dy/dη
        
        det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
        
        print(f"Jacobian matrix:")
        print(f"  J = [[{J[0,0]:7.4f}, {J[0,1]:7.4f}],")
        print(f"       [{J[1,0]:7.4f}, {J[1,1]:7.4f}]]")
        print(f"  det(J) = {det_J:.6f}")
        
        # For a rectangular element, Jacobian should be constant
        # dx/dξ = width/2 = 1.0, dy/dη = height/2 = 0.5
        # dy/dξ = dx/dη = 0 (for aligned rectangle)
        expected_J = np.array([[1.0, 0.0], [0.0, 0.5]])
        expected_det = 0.5
        
        print(f"Expected J for rectangle:")
        print(f"  J = [[{expected_J[0,0]:7.4f}, {expected_J[0,1]:7.4f}],")
        print(f"       [{expected_J[1,0]:7.4f}, {expected_J[1,1]:7.4f}]]")
        print(f"  det(J) = {expected_det:.6f}")
        
        # Check if Jacobian is correct
        J_error = np.max(np.abs(J - expected_J))
        det_error = abs(det_J - expected_det)
        
        print(f"Errors:")
        print(f"  Max J error: {J_error:.6f}")
        print(f"  det(J) error: {det_error:.6f}")
        
        # Compute physical coordinates at this Gauss point
        # x = Σ N_i * x_i, y = Σ N_i * y_i
        N = compute_quad8_shape_functions(xi, eta)
        x_phys = np.sum(N * coords[:, 0])
        y_phys = np.sum(N * coords[:, 1])
        
        print(f"Physical coordinates: ({x_phys:.6f}, {y_phys:.6f})")
        
        # For reference rectangle, should be: x = xi, y = 0.5*eta
        expected_x = xi
        expected_y = 0.5 * eta
        print(f"Expected coordinates: ({expected_x:.6f}, {expected_y:.6f})")
        
        coord_error = np.sqrt((x_phys - expected_x)**2 + (y_phys - expected_y)**2)
        print(f"Coordinate error: {coord_error:.6f}")
        
        # If Jacobian varies between Gauss points, that's the problem!
        if gp == 0:
            reference_J = J.copy()
            reference_det = det_J
        else:
            variation = np.max(np.abs(J - reference_J))
            det_variation = abs(det_J - reference_det)
            print(f"Variation from GP 0:")
            print(f"  Max J variation: {variation:.6f}")
            print(f"  det(J) variation: {det_variation:.6f}")
            
            if variation > 1e-6:
                print(f"  ERROR: Jacobian should be constant for rectangular element!")


def compute_quad8_shape_functions(xi, eta):
    """Compute 8-node quad shape functions at (xi, eta) using serendipity formulation."""
    N = np.zeros(8)
    
    # Corner nodes (serendipity formulation)
    N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1)  # Node 0: (-1,-1)
    N[1] = 0.25 * (1 + xi) * (1 - eta) * ( xi - eta - 1)  # Node 1: ( 1,-1)
    N[2] = 0.25 * (1 + xi) * (1 + eta) * ( xi + eta - 1)  # Node 2: ( 1, 1)
    N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1)  # Node 3: (-1, 1)
    
    # Edge nodes
    N[4] = 0.5 * (1 - xi*xi) * (1 - eta)  # Node 4: ( 0,-1) bottom edge
    N[5] = 0.5 * (1 + xi) * (1 - eta*eta) # Node 5: ( 1, 0) right edge
    N[6] = 0.5 * (1 - xi*xi) * (1 + eta)  # Node 6: ( 0, 1) top edge
    N[7] = 0.5 * (1 - xi) * (1 - eta*eta) # Node 7: (-1, 0) left edge
    
    return N


def test_stress_displacement_chain():
    """Test the complete displacement -> strain -> stress chain."""
    print("\n=== Testing Displacement -> Strain -> Stress Chain ===")
    
    # Use same rectangle
    coords = np.array([
        [-1.0, -0.5],  # 0: bottom-left
        [ 1.0, -0.5],  # 1: bottom-right  
        [ 1.0,  0.5],  # 2: top-right
        [-1.0,  0.5],  # 3: top-left
        [ 0.0, -0.5],  # 4: bottom-mid
        [ 1.0,  0.0],  # 5: right-mid
        [ 0.0,  0.5],  # 6: top-mid
        [-1.0,  0.0],  # 7: left-mid
    ])
    
    # Create a known displacement field: linear vertical variation
    # For a uniform strain field, we need a linear displacement variation
    # u_x = 0 everywhere, u_y = constant * y (linear variation with height)
    # For ε_y = -0.01, we need du_y/dy = -0.01
    # Since element height = 1.0, we need u_y = -0.01 * y
    
    displacements = np.zeros(16)  # 8 nodes × 2 DOF
    for i in range(8):
        y = coords[i, 1]  # y-coordinate of node i
        displacements[2*i]     = 0.0           # u_x = 0
        displacements[2*i + 1] = -0.01 * y     # u_y = -0.01 * y (linear variation)
    
    print(f"Applied displacement field:")
    print(f"  u_x = 0 everywhere")
    print(f"  u_y = -0.01 * y (linear variation with height)")
    print(f"Expected strains: ε_x=0, ε_y=-0.01, γ_xy=0")
    
    print(f"Node displacements:")
    for i in range(8):
        print(f"  Node {i}: u=({displacements[2*i]:7.4f}, {displacements[2*i+1]:7.4f})")
    
    # Test strain calculation at each Gauss point
    from fem_perzyna import compute_quad8_strains_at_xi_eta
    
    gauss_coords = [
        (-1/np.sqrt(3), -1/np.sqrt(3)),
        ( 1/np.sqrt(3), -1/np.sqrt(3)),
        ( 1/np.sqrt(3),  1/np.sqrt(3)),
        (-1/np.sqrt(3),  1/np.sqrt(3)),
    ]
    
    for gp, (xi, eta) in enumerate(gauss_coords):
        strains = compute_quad8_strains_at_xi_eta(coords, displacements, xi, eta)
        print(f"GP {gp}: strains = [{strains[0]:8.5f}, {strains[1]:8.5f}, {strains[2]:8.5f}]")
        
        # Check if strains are correct
        expected_strains = np.array([0.0, -0.01, 0.0])
        strain_error = np.max(np.abs(strains - expected_strains))
        print(f"       strain error = {strain_error:.6f}")
        
        if strain_error > 1e-6:
            print(f"       ERROR: Strains incorrect!")


if __name__ == "__main__":
    test_simple_rectangle()
    test_stress_displacement_chain()