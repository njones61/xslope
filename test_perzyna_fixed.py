#!/usr/bin/env python3
"""
Test fixed Perzyna implementation.
"""

from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

def test_fixed_perzyna():
    """Test the fixed Perzyna implementation."""
    
    print("=== Testing Fixed Perzyna Implementation ===")
    
    # Load Griffiths Example 1
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 20  # Coarse mesh for quick testing
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    
    # Test at F=1.0 (should converge easily)
    print(f"\n=== Testing F=1.0 ===")
    solution = solve_fem_perzyna(fem_data, F=1.0, debug_level=2)
    
    if solution.get("converged", False):
        print(f"✓ F=1.0 converged in {solution['iterations']} iterations")
        print(f"  Max displacement: {solution['max_displacement']:.6f}")
    else:
        print(f"✗ F=1.0 failed: {solution.get('error', 'Unknown')}")
    
    # Test at F=1.5 (should have more plastic behavior)
    print(f"\n=== Testing F=1.5 ===")
    solution = solve_fem_perzyna(fem_data, F=1.5, debug_level=2)
    
    if solution.get("converged", False):
        print(f"✓ F=1.5 converged in {solution['iterations']} iterations")
        print(f"  Max displacement: {solution['max_displacement']:.6f}")
    else:
        print(f"✗ F=1.5 failed to converge in {solution['iterations']} iterations")
        print(f"  This indicates slope failure at F=1.5")
    
    # Test at F=2.0 (should definitely fail)
    print(f"\n=== Testing F=2.0 ===")
    solution = solve_fem_perzyna(fem_data, F=2.0, debug_level=1)
    
    if solution.get("converged", False):
        print(f"✓ F=2.0 converged in {solution['iterations']} iterations")
        print(f"  Max displacement: {solution['max_displacement']:.6f}")
        print(f"  ⚠ Unexpected - slope should fail at F=2.0")
    else:
        print(f"✓ F=2.0 failed to converge in {solution['iterations']} iterations")
        print(f"  Expected failure behavior")

if __name__ == "__main__":
    test_fixed_perzyna()