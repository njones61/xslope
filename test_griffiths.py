#!/usr/bin/env python3
"""
Test script to compare Griffiths algorithm with current implementation.

This script tests the new Griffiths & Lane (1999) implementation against
the existing FEM solver to see if we can achieve circular failure surfaces.
"""

from fem import build_fem_data, solve_fem, solve_ssrm
from fem_griffiths import solve_fem_griffiths, solve_ssrm_griffiths
from plot_fem import plot_fem_results, plot_fem_data
from fileio import load_slope_data, print_dictionary
from mesh import build_polygons, build_mesh_from_polygons
from plot import plot_inputs, plot_mesh, plot_polygons
import numpy as np
import matplotlib.pyplot as plt

def compare_algorithms():
    """Compare Griffiths algorithm with current implementation."""
    
    print("=== Comparing Griffiths vs Current FEM Implementation ===\n")
    
    # Load test case
    slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
    
    # Build mesh
    polygons = build_polygons(slope_data)
    target_size = 2
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"Mesh info:")
    print(f"  Nodes: {len(fem_data['nodes'])}")
    print(f"  Elements: {len(fem_data['elements'])}")
    print(f"  Element types: {np.unique(fem_data['element_types'])}")
    
    # Test factor of safety values
    test_F_values = [1.0, 1.2, 1.4, 1.6]
    
    print(f"\n{'F Value':<8} {'Current':<12} {'Griffiths':<12} {'Max Disp Curr':<15} {'Max Disp Grif':<15}")
    print("-" * 70)
    
    results_current = []
    results_griffiths = []
    
    for F in test_F_values:
        print(f"{F:<8.1f}", end=" ")
        
        # Test current implementation
        try:
            sol_current = solve_fem(fem_data, F=F, debug_level=0, plastic_stiffness_reduction=0.01)
            if sol_current.get("converged", False):
                max_disp_curr = np.max(np.abs(sol_current.get('displacements', [0])))
                status_curr = "CONVERGED"
                results_current.append((F, True, max_disp_curr, sol_current))
            else:
                max_disp_curr = float('inf')
                status_curr = "FAILED"
                results_current.append((F, False, max_disp_curr, sol_current))
        except Exception as e:
            max_disp_curr = float('inf')
            status_curr = f"ERROR: {str(e)[:20]}"
            results_current.append((F, False, max_disp_curr, None))
        
        print(f"{status_curr:<12}", end=" ")
        
        # Test Griffiths implementation
        try:
            sol_griffiths = solve_fem_griffiths(fem_data, F=F, debug_level=0)
            if sol_griffiths.get("converged", False):
                max_disp_grif = np.max(np.abs(sol_griffiths.get('displacements', [0])))
                status_grif = "CONVERGED"
                results_griffiths.append((F, True, max_disp_grif, sol_griffiths))
            else:
                max_disp_grif = float('inf')
                status_grif = "FAILED"
                results_griffiths.append((F, False, max_disp_grif, sol_griffiths))
        except Exception as e:
            max_disp_grif = float('inf')
            status_grif = f"ERROR: {str(e)[:20]}"
            results_griffiths.append((F, False, max_disp_grif, None))
        
        print(f"{status_grif:<12}", end=" ")
        
        # Print displacement magnitudes
        if max_disp_curr != float('inf'):
            print(f"{max_disp_curr:<15.6f}", end=" ")
        else:
            print(f"{'N/A':<15}", end=" ")
            
        if max_disp_grif != float('inf'):
            print(f"{max_disp_grif:<15.6f}")
        else:
            print(f"{'N/A':<15}")
    
    print("\n" + "="*70)
    
    # Find critical F values using SSRM
    print("\n=== SSRM Comparison ===")
    
    # Current SSRM
    print("\nTesting current SSRM...")
    try:
        ssrm_current = solve_ssrm(fem_data, debug_level=1)
        if ssrm_current.get("converged", False):
            fs_current = ssrm_current.get("FS", None)
            print(f"Current SSRM: FS = {fs_current:.4f}")
        else:
            print(f"Current SSRM failed: {ssrm_current.get('error', 'Unknown error')}")
            fs_current = None
    except Exception as e:
        print(f"Current SSRM error: {e}")
        fs_current = None
    
    # Griffiths SSRM
    print("\nTesting Griffiths SSRM...")
    try:
        ssrm_griffiths = solve_ssrm_griffiths(fem_data, F_min=1.0, F_max=3.0, tolerance=0.01, debug_level=1)
        if ssrm_griffiths.get("converged", False):
            fs_griffiths = ssrm_griffiths.get("FS", None)
            print(f"Griffiths SSRM: FS = {fs_griffiths:.4f}")
        else:
            print(f"Griffiths SSRM failed: {ssrm_griffiths.get('error', 'Unknown error')}")
            fs_griffiths = None
    except Exception as e:
        print(f"Griffiths SSRM error: {e}")
        fs_griffiths = None
    
    # Plot comparisons for interesting cases
    print("\n=== Visualization Comparison ===")
    
    # Find a case where both converged for comparison
    comparison_F = None
    for F, curr_conv, _, curr_sol in results_current:
        for F2, grif_conv, _, grif_sol in results_griffiths:
            if F == F2 and curr_conv and grif_conv and curr_sol and grif_sol:
                comparison_F = F
                comparison_curr = curr_sol
                comparison_grif = grif_sol
                break
        if comparison_F:
            break
    
    if comparison_F:
        print(f"\nPlotting comparison for F = {comparison_F}")
        
        # Plot current results
        print("Plotting current implementation results...")
        plot_fem_results(fem_data, comparison_curr, plot_type='stress, shear_strain, deformation')
        
        # Plot Griffiths results  
        print("Plotting Griffiths implementation results...")
        plot_fem_results(fem_data, comparison_grif, plot_type='stress, shear_strain, deformation')
    else:
        print("No suitable case found for visualization comparison")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY:")
    print("="*70)
    
    if fs_current and fs_griffiths:
        print(f"Critical Factor of Safety:")
        print(f"  Current Implementation: {fs_current:.4f}")
        print(f"  Griffiths Implementation: {fs_griffiths:.4f}")
        print(f"  Difference: {abs(fs_current - fs_griffiths):.4f}")
        
        if abs(fs_current - fs_griffiths) < 0.1:
            print("  → Similar results - algorithms agree")
        else:
            print("  → Different results - investigate further")
    
    print(f"\nKey Differences Expected:")
    print(f"  • Griffiths: Pure non-convergence failure criterion")
    print(f"  • Current: Multiple failure criteria (displacement + non-convergence)")
    print(f"  • Griffiths: K₀ initial stress state")
    print(f"  • Current: Zero initial stress state")
    print(f"  • Griffiths: Visco-plastic stress redistribution")
    print(f"  • Current: Elastic-plastic with stiffness reduction")


def test_single_griffiths():
    """Test single Griffiths analysis with detailed output."""
    
    print("=== Detailed Griffiths Test ===\n")
    
    # Load test case
    slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
    polygons = build_polygons(slope_data)
    target_size = 2
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    fem_data = build_fem_data(slope_data, mesh)
    
    # Test at F = 1.3 (should be near critical)
    F_test = 1.3
    
    print(f"Testing Griffiths algorithm at F = {F_test}")
    solution = solve_fem_griffiths(fem_data, F=F_test, debug_level=3)
    
    if solution.get("converged", False):
        print(f"\n✓ Griffiths algorithm CONVERGED")
        print(f"  Iterations: {solution.get('iterations', 'Unknown')}")
        print(f"  Max displacement: {solution.get('max_displacement', 0):.6f}")
        print(f"  Plastic elements: {np.sum(solution.get('plastic_elements', []))}")
        
        # Plot results
        plot_fem_results(fem_data, solution, plot_type='stress, shear_strain, deformation')
    else:
        print(f"\n✗ Griffiths algorithm FAILED TO CONVERGE")
        print(f"  Error: {solution.get('error', 'Unknown error')}")
        print(f"  Iterations: {solution.get('iterations', 'Unknown')}")


if __name__ == "__main__":
    # Run tests
    print("Starting Griffiths vs Current Implementation Comparison\n")
    
    # Test single case first
    test_single_griffiths()
    print("\n" + "="*80 + "\n")
    
    # Full comparison
    compare_algorithms()
    
    print("\nTesting complete!")