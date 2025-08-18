#!/usr/bin/env python3
"""
Debug plasticity in Perzyna implementation.
"""

import numpy as np
from math import sin, cos, sqrt, radians
from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

def debug_element_stresses():
    """Debug stress calculation at element level."""
    
    print("=== Debugging Element-Level Plasticity ===")
    
    # Load simple case
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 15  # Very coarse for debugging
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Debug mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    
    # Material properties
    c_val = fem_data['c_by_mat'][0]
    phi_val = fem_data['phi_by_mat'][0]
    
    print(f"Material: c={c_val:.1f}, φ={phi_val:.1f}°")
    
    # Test at F = 1.4 (should be near failure)
    F = 1.4
    c_reduced = c_val / F
    phi_reduced = radians(phi_val) / F
    
    print(f"Reduced at F={F}: c={c_reduced:.1f}, φ={phi_reduced*180/np.pi:.1f}°")
    
    # Compute theoretical yield stress for this material
    print(f"\n=== Theoretical Yield Analysis ===")
    
    # For φ=20°, c=312.5/1.4=223.2
    # Maximum principal stress difference for failure:
    # (σ1 - σ3) = 2c*cos(φ) + (σ1 + σ3)*sin(φ)
    
    # Simple case: uniaxial compression (σ3=0)
    # σ1 = 2c*cos(φ)/(1-sin(φ))
    
    cos_phi = cos(phi_reduced)
    sin_phi = sin(phi_reduced)
    
    if sin_phi < 0.999:  # Avoid division by zero
        sigma1_yield = 2 * c_reduced * cos_phi / (1 - sin_phi)
        print(f"Uniaxial yield stress: σ1 = {sigma1_yield:.1f}")
        print(f"Expected plastic zones where stress > {sigma1_yield:.1f}")
    
    # Test manual stress calculation
    print(f"\n=== Manual Stress Check ===")
    test_manual_yield_criterion(c_reduced, phi_reduced)

def test_manual_yield_criterion(c, phi):
    """Test yield criterion with known stress states."""
    
    print(f"Testing Mohr-Coulomb with c={c:.1f}, φ={phi*180/np.pi:.1f}°")
    
    # Test cases: [sig_x, sig_y, tau_xy]
    test_stresses = [
        [0, 0, 0],           # No stress - should be elastic
        [-100, -100, 0],     # Hydrostatic compression - should be elastic  
        [-500, 0, 0],        # High compression - should yield
        [-1000, -100, 50],   # Complex stress - check yield
    ]
    
    for i, stress in enumerate(test_stresses):
        sig_x, sig_y, tau_xy = stress
        violation = check_mohr_coulomb_manual(stress, c, phi)
        status = "PLASTIC" if violation > 1e-6 else "ELASTIC"
        print(f"Test {i+1}: σ=[{sig_x:4.0f}, {sig_y:4.0f}, {tau_xy:4.0f}] → F={violation:.3f} ({status})")

def check_mohr_coulomb_manual(stress, c, phi):
    """Manual implementation of Mohr-Coulomb yield criterion."""
    
    sig_x, sig_y, tau_xy = stress
    
    # Principal stresses
    sig_mean = (sig_x + sig_y) / 2
    tau_max = sqrt(((sig_x - sig_y) / 2)**2 + tau_xy**2)
    
    sig1 = sig_mean + tau_max  # Major principal (less negative)
    sig3 = sig_mean - tau_max  # Minor principal (more negative)
    
    # Mohr-Coulomb failure criterion 
    # F = (sig1 - sig3) - (sig1 + sig3)*sin(phi) - 2*c*cos(phi)
    # Note: Using compression as negative convention
    
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    
    F = (sig1 - sig3) + (sig1 + sig3) * sin_phi - 2 * c * cos_phi
    
    return max(0, F)

def debug_fem_stress_calculation():
    """Debug actual FEM stress calculation."""
    
    print(f"\n=== Debugging FEM Stress Calculation ===")
    
    # Load and solve at F=1.4
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 15
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    # Get a solution
    solution = solve_fem_perzyna(fem_data, F=1.4, debug_level=1)
    
    if solution.get("converged", False):
        displacements = solution["displacements"]
        max_disp = np.max(np.abs(displacements))
        print(f"✓ Converged with max displacement: {max_disp:.6f}")
        
        # Check if displacement is reasonable
        if max_disp < 1e-6:
            print("⚠ Displacements very small - may indicate rigid behavior")
        elif max_disp > 1.0:
            print("⚠ Displacements very large - may indicate instability")
        else:
            print("✓ Displacements in reasonable range")
            
        # Manual stress calculation for one element
        elements = fem_data["elements"]
        nodes = fem_data["nodes"]
        
        if len(elements) > 0:
            elem_idx = 0  # First element
            elem = elements[elem_idx]
            elem_nodes = elem[:3]  # Triangle
            elem_coords = nodes[elem_nodes]
            
            print(f"\nElement {elem_idx} analysis:")
            print(f"  Nodes: {elem_nodes}")
            print(f"  Coordinates: {elem_coords}")
            
            # Get element displacements
            elem_disp = np.zeros(6)
            for i, node in enumerate(elem_nodes):
                elem_disp[2*i] = displacements[2*node]
                elem_disp[2*i+1] = displacements[2*node+1]
            
            print(f"  Displacements: {elem_disp}")
            
            # Compute element strains manually
            strains = compute_triangle_strains_manual(elem_coords, elem_disp)
            print(f"  Strains: {strains}")
            
            # Compute stresses  
            E = fem_data['E_by_mat'][0]
            nu = fem_data['nu_by_mat'][0]
            stresses = compute_triangle_stresses_manual(strains, E, nu)
            print(f"  Stresses: {stresses}")
            
            # Check yield
            c = fem_data['c_by_mat'][0] / 1.4
            phi = radians(fem_data['phi_by_mat'][0]) / 1.4
            yield_val = check_mohr_coulomb_manual(stresses, c, phi)
            print(f"  Yield function: {yield_val:.6f}")
            
            if yield_val > 1e-6:
                print("  → Should be PLASTIC")
            else:
                print("  → Should be ELASTIC")
    else:
        print("✗ Solution failed to converge")

def compute_triangle_strains_manual(coords, displacements):
    """Manually compute triangle strains."""
    
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    if area < 1e-12:
        return np.array([0.0, 0.0, 0.0])
    
    # Shape function derivatives
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2
    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1
    
    # B matrix
    B = np.array([
        [b1, 0,  b2, 0,  b3, 0 ],
        [0,  c1, 0,  c2, 0,  c3],
        [c1, b1, c2, b2, c3, b3]
    ]) / (2 * area)
    
    # Strains
    strains = B @ displacements
    return strains

def compute_triangle_stresses_manual(strains, E, nu):
    """Manually compute triangle stresses from strains."""
    
    # Constitutive matrix (plane strain)
    factor = E / ((1 + nu) * (1 - 2*nu))
    D = factor * np.array([
        [1-nu, nu,   0        ],
        [nu,   1-nu, 0        ],
        [0,    0,    (1-2*nu)/2]
    ])
    
    stresses = D @ strains
    return stresses

if __name__ == "__main__":
    debug_element_stresses()
    debug_fem_stress_calculation()