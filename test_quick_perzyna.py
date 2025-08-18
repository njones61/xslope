#!/usr/bin/env python3
"""
Quick test to identify specific issues in Perzyna implementation.
"""

from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

def test_quick():
    """Quick test with minimal output."""
    
    print("=== Quick Perzyna Test ===")
    
    # Load case
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 30  # Very coarse for speed
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    
    # Test at F=1.5 (should fail)
    try:
        solution = solve_fem_perzyna(fem_data, F=1.5, debug_level=0)  # Minimal debug
        
        print(f"Converged: {solution.get('converged', False)}")
        print(f"Iterations: {solution.get('iterations', 0)}")
        print(f"Plastic elements: {np.sum(solution.get('plastic_elements', []))}")
        
        if solution.get('converged', False):
            print("⚠ F=1.5 converged - unexpected")
        else:
            print("✓ F=1.5 failed to converge - expected")
            
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    import numpy as np
    test_quick()