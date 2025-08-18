#!/usr/bin/env python3
"""
Debug yield criterion and stress states in Perzyna implementation.
"""

from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

def debug_yield_criterion():
    """Debug the yield criterion detection."""
    
    print("=== Debugging Yield Criterion ===")
    
    # Load simple case
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 20  # Coarse mesh
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Material properties:")
    print(f"  c = {fem_data['c_by_mat'][0]:.1f} psf")
    print(f"  φ = {fem_data['phi_by_mat'][0]:.1f}°")
    print(f"  γ = {fem_data['gamma_by_mat'][0]:.1f} pcf")
    
    # Test at F=1.5 where we expect plastic behavior
    print(f"\n=== Testing F=1.5 with detailed debugging ===")
    solution = solve_fem_perzyna(fem_data, F=1.5, debug_level=3)
    
    print(f"\nResults:")
    print(f"  Converged: {solution.get('converged', False)}")
    print(f"  Iterations: {solution.get('iterations', 0)}")
    print(f"  Max displacement: {solution.get('max_displacement', 0):.6f}")
    
    # Check for plastic elements
    if 'plastic_elements' in solution:
        plastic_count = np.sum(solution['plastic_elements'])
        print(f"  Plastic elements: {plastic_count}")
        
        if plastic_count > 0:
            plastic_indices = np.where(solution['plastic_elements'])[0]
            print(f"  Plastic element indices: {plastic_indices[:10]}")  # Show first 10
        else:
            print(f"  ⚠ No plastic elements detected - yield criterion may not be working")
    
    # Check stress magnitudes
    if 'stresses' in solution:
        stresses = solution['stresses']
        print(f"\nStress analysis:")
        print(f"  Number of elements: {len(stresses)}")
        
        # Find element with maximum stress
        max_stress_elem = -1
        max_stress_mag = 0
        
        for elem_idx, stress in stresses.items():
            # stress = [sig_x, sig_y, tau_xy]
            stress_mag = np.linalg.norm(stress)
            if stress_mag > max_stress_mag:
                max_stress_mag = stress_mag
                max_stress_elem = elem_idx
        
        print(f"  Max stress magnitude: {max_stress_mag:.1f} psf in element {max_stress_elem}")
        
        if max_stress_elem >= 0:
            stress = stresses[max_stress_elem]
            print(f"  Max stress element stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
            
            # Check yield function manually
            from fem_perzyna import check_mohr_coulomb_perzyna
            from math import radians
            
            c = fem_data['c_by_mat'][0] / 1.5  # Reduced strength
            phi_rad = radians(fem_data['phi_by_mat'][0])
            
            f_yield = check_mohr_coulomb_perzyna(stress, c, phi_rad)
            print(f"  Yield function value: {f_yield:.6f}")
            print(f"  c_reduced: {c:.1f} psf")
            
            if f_yield > 1e-8:
                print(f"  ✓ Element should be plastic (f > 1e-8)")
            else:
                print(f"  ⚠ Element should be elastic (f ≤ 1e-8)")

if __name__ == "__main__":
    debug_yield_criterion()