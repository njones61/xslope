#!/usr/bin/env python3
"""
Debug Perzyna implementation step by step.
"""

from fem_perzyna import solve_fem_perzyna
from fem import build_fem_data
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

def debug_perzyna_simple():
    """Debug basic Perzyna functionality."""
    
    print("=== Debugging Perzyna Implementation ===")
    
    # Load test case
    try:
        slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
        print("✓ Successfully loaded slope data")
    except Exception as e:
        print(f"✗ Failed to load slope data: {e}")
        return
    
    # Build mesh with triangles first (simpler)
    try:
        polygons = build_polygons(slope_data)
        target_size = 10  # Large elements for quick testing
        mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
        fem_data = build_fem_data(slope_data, mesh)
        print(f"✓ Built mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    except Exception as e:
        print(f"✗ Failed to build mesh: {e}")
        return
    
    # Check material properties
    print("\n=== Material Properties ===")
    print(f"c: {fem_data['c_by_mat'][0]:.1f}")
    print(f"φ: {fem_data['phi_by_mat'][0]:.1f}°")
    print(f"γ: {fem_data['gamma_by_mat'][0]:.1f}")
    print(f"E: {fem_data['E_by_mat'][0]:.0f}")
    print(f"ν: {fem_data['nu_by_mat'][0]:.2f}")
    
    # Test at F=1.0 (should definitely converge)
    print(f"\n=== Testing Perzyna at F=1.0 ===")
    try:
        solution = solve_fem_perzyna(fem_data, F=1.0, debug_level=2)
        
        if solution.get("converged", False):
            print("✓ F=1.0 converged successfully")
            print(f"  Iterations: {solution.get('iterations', 'Unknown')}")
            print(f"  Max displacement: {solution.get('max_displacement', 0):.6f}")
        else:
            print("✗ F=1.0 failed to converge")
            print(f"  Error: {solution.get('error', 'Unknown')}")
            print(f"  Iterations: {solution.get('iterations', 'Unknown')}")
            
    except Exception as e:
        print(f"✗ Error in Perzyna analysis: {e}")
        import traceback
        traceback.print_exc()

def debug_matrix_issues():
    """Debug matrix singularity issues."""
    
    print("\n=== Debugging Matrix Issues ===")
    
    # Load minimal test case
    slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 20  # Very coarse mesh
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Minimal mesh: {len(fem_data['nodes'])} nodes, {len(fem_data['elements'])} elements")
    
    # Check boundary conditions
    bc_type = fem_data["bc_type"]
    n_fixed = np.sum(bc_type == 1)
    n_xroller = np.sum(bc_type == 2) 
    n_yroller = np.sum(bc_type == 3)
    n_free = np.sum(bc_type == 0)
    
    print(f"Boundary conditions:")
    print(f"  Fixed: {n_fixed}")
    print(f"  X-roller: {n_xroller}")
    print(f"  Y-roller: {n_yroller}")
    print(f"  Free: {n_free}")
    
    # Check for adequate constraints
    n_dof = 2 * len(fem_data['nodes'])
    n_constrained = 2*n_fixed + n_xroller + n_yroller
    
    print(f"DOF: {n_dof}, Constrained: {n_constrained}")
    
    if n_constrained < 3:
        print("⚠ WARNING: Insufficient constraints for rigid body motion")
    else:
        print("✓ Adequate constraints")

if __name__ == "__main__":
    debug_matrix_issues()
    debug_perzyna_simple()