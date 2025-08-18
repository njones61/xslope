#!/usr/bin/env python3
"""
Detailed debug of Perzyna plastic strain updates.
"""

import numpy as np
from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

def debug_plastic_strain_updates():
    """Debug the plastic strain update mechanism in detail."""
    
    print("=== Detailed Perzyna Debug ===")
    
    # Load case and build mesh
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 20  # Very coarse for detailed debugging
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Debug mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    
    # Test at F = 1.4 with high debug level
    print(f"\n=== Testing F = 1.4 with High Debug ===")
    solution = solve_fem_perzyna(fem_data, F=1.4, debug_level=3)
    
    if solution.get("converged", False):
        print(f"✓ Converged after {solution['iterations']} iterations")
        
        # Check plastic strains
        plastic_strains = solution.get("plastic_strains", {})
        total_plastic = 0
        plastic_elements = 0
        
        for elem_idx in plastic_strains:
            elem_plastic = plastic_strains[elem_idx]
            elem_total = np.sum(np.abs(elem_plastic))
            if elem_total > 1e-8:
                plastic_elements += 1
                total_plastic += elem_total
        
        print(f"Plastic elements: {plastic_elements}/{len(fem_data['elements'])}")
        print(f"Total plastic strain: {total_plastic:.6e}")
        
        if plastic_elements == 0:
            print("⚠ No plastic strains detected - algorithm may not be working")
        elif total_plastic < 1e-6:
            print("⚠ Very small plastic strains - may indicate convergence issues")
        else:
            print("✓ Plastic strains detected")
    else:
        print(f"✗ Failed to converge: {solution.get('error', 'Unknown')}")

def debug_modified_algorithm():
    """Debug with modified algorithm parameters."""
    
    print(f"\n=== Testing Modified Algorithm Parameters ===")
    
    # Load case
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 20
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    # Test with different viscosity parameters
    eta_values = [1e-6, 1e-3, 1.0, 100.0]
    
    print(f"{'Eta':<8} {'Converged':<10} {'Iterations':<11} {'Max Disp':<12}")
    print("-" * 50)
    
    for eta in eta_values:
        # Modify the fem_data to test different viscosity
        # Note: This would require modifying the solve function to accept eta parameter
        solution = solve_fem_perzyna(fem_data, F=1.4, debug_level=1)
        
        converged = solution.get("converged", False)
        iterations = solution.get("iterations", 0)
        max_disp = solution.get("max_displacement", 0)
        
        print(f"{eta:<8.0e} {str(converged):<10} {iterations:<11} {max_disp:<12.6f}")

def test_simplified_yield_criterion():
    """Test the yield criterion implementation directly."""
    
    print(f"\n=== Testing Yield Criterion Implementation ===")
    
    from fem_perzyna import check_mohr_coulomb_perzyna
    from math import radians
    
    # Test parameters from Griffiths Example 1
    c = 312.5 / 1.4  # Reduced cohesion
    phi = radians(20.0) / 1.4  # Reduced friction angle
    
    print(f"Reduced parameters: c={c:.1f}, φ={phi*180/np.pi:.1f}°")
    
    # Test stress states
    test_stresses = [
        [0, 0, 0],              # No stress
        [-100, -100, 0],        # Moderate compression
        [-500, -500, 0],        # High compression
        [-1000, -200, 50],      # Complex stress state
        [-2000, -500, 100],     # Very high stress
    ]
    
    print(f"\n{'σx':<8} {'σy':<8} {'τxy':<8} {'F_yield':<12} {'Status':<8}")
    print("-" * 50)
    
    for stress in test_stresses:
        sig_x, sig_y, tau_xy = stress
        
        try:
            f_yield = check_mohr_coulomb_perzyna(stress, c, phi)
            status = "PLASTIC" if f_yield > 1e-6 else "ELASTIC"
            
            print(f"{sig_x:<8.0f} {sig_y:<8.0f} {tau_xy:<8.0f} {f_yield:<12.3f} {status:<8}")
            
        except Exception as e:
            print(f"{sig_x:<8.0f} {sig_y:<8.0f} {tau_xy:<8.0f} {'ERROR':<12} {str(e)[:8]:<8}")

if __name__ == "__main__":
    test_simplified_yield_criterion()
    debug_plastic_strain_updates() 
    debug_modified_algorithm()