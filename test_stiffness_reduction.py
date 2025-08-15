#!/usr/bin/env python3

from fem import build_fem_data, solve_fem
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

# Load data and build FEM model
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
target_size = 2
mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

print("=== Testing Different Plastic Stiffness Reduction Factors ===")
print("F=1.0 (near-failure condition)\n")

# Test different plastic stiffness reduction factors
reduction_factors = [0.01, 0.05, 0.1, 0.2, 0.5]

for reduction in reduction_factors:
    print(f"Plastic stiffness reduction = {reduction}")
    
    solution = solve_fem(fem_data, F=1.0, debug_level=0, plastic_stiffness_reduction=reduction)
    
    # Extract results
    displacements = solution['displacements']
    u = displacements[0::2]  # x-displacements
    v = displacements[1::2]  # y-displacements
    max_disp = np.max(np.sqrt(u**2 + v**2))
    plastic_count = np.sum(solution['plastic_elements'])
    total_elements = len(fem_data['elements'])
    converged = solution['converged']
    iterations = solution['iterations']
    
    print(f"  Max displacement: {max_disp:.4f} ft")
    print(f"  Plastic elements: {plastic_count}/{total_elements} ({100*plastic_count/total_elements:.1f}%)")
    print(f"  Converged: {converged} in {iterations} iterations")
    print(f"  U range: [{np.min(u):.4f}, {np.max(u):.4f}] ft")
    print(f"  V range: [{np.min(v):.4f}, {np.max(v):.4f}] ft")
    print()

print("=== Analysis ===")
print("Conservative values (0.2-0.5): Less stiffness reduction, smaller displacements")
print("Moderate values (0.1-0.2): Balanced behavior, typical engineering practice")  
print("Aggressive values (0.01-0.05): More stiffness reduction, larger displacements")
print("Too aggressive (<0.01): Risk of excessive/unrealistic deformations")