#!/usr/bin/env python3
"""
Test SSRM with improved Perzyna implementation.
"""

from fem_perzyna import solve_fem_perzyna, solve_ssrm_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

def test_ssrm_focused():
    """Test SSRM with focus on finding failure."""
    
    print("=== SSRM Test with Improved Perzyna ===")
    
    # Load Griffiths Example 1
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 8  # Finer mesh for better accuracy
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    print(f"Material: c={fem_data['c_by_mat'][0]:.1f}, φ={fem_data['phi_by_mat'][0]:.1f}°")
    print(f"Target: FS ≈ 1.4 (Griffiths), Spencer = 1.376")
    
    # Test specific F values to find failure range
    print(f"\n=== Finding Failure Range ===")
    test_F_values = [1.0, 1.5, 2.0, 2.5, 3.0]
    
    print(f"{'F':<6} {'Status':<12} {'Iterations':<11}")
    print("-" * 35)
    
    failure_range = None
    
    for F in test_F_values:
        print(f"{F:<6.1f} ", end="")
        
        try:
            solution = solve_fem_perzyna(fem_data, F=F, debug_level=0)
            converged = solution.get("converged", False)
            iterations = solution.get("iterations", 0)
            
            if converged:
                print(f"{'CONVERGED':<12} {iterations:<11}")
            else:
                print(f"{'FAILED':<12} {iterations:<11}")
                if not failure_range:
                    failure_range = (test_F_values[test_F_values.index(F)-1] if F != test_F_values[0] else F, F)
                    
        except Exception as e:
            print(f"{'ERROR':<12} {str(e)[:9]:<11}")
    
    # Run SSRM based on findings
    print(f"\n=== SSRM Analysis ===")
    
    if failure_range:
        F_min, F_max = failure_range
        print(f"Detected failure between F={F_min:.1f} and F={F_max:.1f}")
        F_max = F_max + 0.5  # Extend range slightly
    else:
        print("No failure detected in test range, using wider SSRM range")
        F_min, F_max = 1.0, 4.0
    
    print(f"Running SSRM with F_min={F_min:.1f}, F_max={F_max:.1f}")
    
    try:
        ssrm_result = solve_ssrm_perzyna(fem_data, F_min=F_min, F_max=F_max, 
                                        tolerance=0.05, debug_level=2)
        
        if ssrm_result.get("converged", False):
            fs = ssrm_result.get("FS", None)
            print(f"\n✓ SSRM Success: FS = {fs:.3f}")
            print(f"Target Griffiths: FS = 1.400")
            print(f"Spencer:          FS = 1.376")
            
            if fs:
                error_griffiths = abs(fs - 1.4)
                error_spencer = abs(fs - 1.376)
                print(f"Error vs Griffiths: {error_griffiths:.3f}")
                print(f"Error vs Spencer:   {error_spencer:.3f}")
                
                if error_griffiths < 0.1:
                    print("🎯 EXCELLENT: Close to Griffiths target!")
                elif error_spencer < 0.1:
                    print("✅ GOOD: Close to Spencer result!")
                elif fs > 1.2 and fs < 1.6:
                    print("✓ Reasonable: In expected range")
                else:
                    print("⚠ Different: Needs investigation")
        else:
            print(f"\n✗ SSRM Failed: {ssrm_result.get('error', 'Unknown error')}")
            
    except Exception as e:
        print(f"\n✗ SSRM Error: {e}")
    
    print(f"\n=== Summary ===")
    print(f"Goal: Verify Perzyna implementation can reproduce Griffiths FS ≈ 1.4")
    print(f"Method: Non-convergence failure criterion (1000 iterations)")

if __name__ == "__main__":
    test_ssrm_focused()