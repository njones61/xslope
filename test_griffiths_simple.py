#!/usr/bin/env python3
"""
Simple test of Griffiths Example 1 to debug issues.
"""

from fem_perzyna import solve_fem_perzyna, solve_ssrm_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
from solve import spencer
from slice import generate_slices
import numpy as np

def test_griffiths_example1_simple():
    """Test Griffiths Example 1 with debugging."""
    
    print("=== Griffiths Example 1 Debug Test ===")
    
    # Load Griffiths Example 1
    try:
        slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
        print("✓ Loaded Griffiths Example 1 data")
    except Exception as e:
        print(f"✗ Failed to load data: {e}")
        return
    
    # Check material properties
    print(f"\n=== Material Properties ===")
    materials = slope_data['materials']
    print(f"Materials type: {type(materials)}")
    
    # Handle both dict and list cases
    if isinstance(materials, list):
        mat = materials[0]  # First material
    else:
        mat = materials  # Single material dict
    
    print(f"c: {mat['c']:.1f}")
    print(f"φ: {mat['phi']:.1f}°")
    print(f"γ: {mat['gamma']:.1f}")
    
    print(f"φ: {mat['phi']:.1f}° (Griffiths Example 1)")
    print(f"Spencer FS: 1.376 (user confirmed)")
    print(f"Target FE FS: 1.400 (Griffiths paper)")
    
    # Build FE mesh with triangles first
    print(f"\n=== Building FE Mesh ===")
    try:
        polygons = build_polygons(slope_data)
        target_size = 8  # Moderate mesh size
        mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
        fem_data = build_fem_data(slope_data, mesh)
        
        print(f"Mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
        print(f"FE material properties:")
        print(f"  c: {fem_data['c_by_mat'][0]:.1f}")
        print(f"  φ: {fem_data['phi_by_mat'][0]:.1f}°")
        print(f"  γ: {fem_data['gamma_by_mat'][0]:.1f}")
        print(f"  E: {fem_data['E_by_mat'][0]:.0f}")
        print(f"  ν: {fem_data['nu_by_mat'][0]:.2f}")
        
    except Exception as e:
        print(f"✗ Mesh build failed: {e}")
        return
    
    # Test Perzyna at key F values
    test_F_values = [1.0, 1.2, 1.3, 1.35, 1.4, 1.45]
    
    print(f"\n=== Perzyna Analysis ===")
    print(f"{'F':<6} {'Status':<12} {'Iter':<6} {'Max Disp':<10}")
    print("-" * 40)
    
    for F in test_F_values:
        try:
            solution = solve_fem_perzyna(fem_data, F=F, debug_level=1)
            
            if solution.get("converged", False):
                status = "CONVERGED"
                iterations = solution.get("iterations", 0)
                max_disp = solution.get("max_displacement", 0)
                print(f"{F:<6.2f} {status:<12} {iterations:<6} {max_disp:<10.6f}")
            else:
                status = "FAILED"
                iterations = solution.get("iterations", 0)
                print(f"{F:<6.2f} {status:<12} {iterations:<6} {'---':<10}")
                print(f"  → First failure at F = {F:.2f}")
                break
                
        except Exception as e:
            print(f"{F:<6.2f} ERROR: {str(e)[:20]}")
            break
    
    # Try SSRM if basic tests work
    print(f"\n=== Perzyna SSRM ===")
    try:
        ssrm_result = solve_ssrm_perzyna(fem_data, F_min=1.0, F_max=2.0, tolerance=0.02, debug_level=1)
        
        if ssrm_result.get("converged", False):
            fs_fe = ssrm_result.get("FS", None)
            print(f"Perzyna SSRM: FS = {fs_fe:.3f}")
            print(f"Target:       FS = 1.400 (Griffiths)")
            print(f"Spencer:      FS = 1.376 (expected)")
            
            if fs_fe:
                error_vs_target = abs(fs_fe - 1.4)
                print(f"Error vs Griffiths: {error_vs_target:.3f}")
                
                if error_vs_target < 0.05:
                    print("✓ Good agreement with Griffiths!")
                elif abs(fs_fe - 1.376) < 0.05:
                    print("✓ Good agreement with Spencer!")
                else:
                    print("⚠ Different result - needs investigation")
        else:
            print(f"SSRM failed: {ssrm_result.get('error', 'Unknown')}")
            
    except Exception as e:
        print(f"SSRM error: {e}")
    
    print(f"\n=== Summary ===")
    print(f"Goal: Reproduce Griffiths Example 1 with FS ≈ 1.4")
    print(f"Key: Look for rotational failure pattern, not horizontal base failure")

if __name__ == "__main__":
    test_griffiths_example1_simple()