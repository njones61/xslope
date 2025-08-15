#!/usr/bin/env python3

from fem import build_fem_data, solve_fem
from plot_fem import plot_fem_results
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons

# Load data and build FEM model
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 2, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

# Run FEM analysis
solution = solve_fem(fem_data, F=1.0, debug_level=0)

print("=== Testing Single Deformation Plot ===")
print("This should auto-scale properly as a standalone plot")

# Test single deformation plot (should still auto-scale properly)
plot_fem_results(fem_data, solution, plot_type='deformation')