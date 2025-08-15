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

print("=== Testing Multiple Plot Types ===\n")

print("1. Single stress plot:")
plot_fem_results(fem_data, solution, plot_type='stress')

print("\n2. Single deformation plot:")
plot_fem_results(fem_data, solution, plot_type='deformation')

print("\n3. Combined stress and deformation:")
plot_fem_results(fem_data, solution, plot_type='stress,deformation')

print("\n4. All three plot types:")
plot_fem_results(fem_data, solution, plot_type='displacement,stress,deformation')

print("\n5. Custom order:")
plot_fem_results(fem_data, solution, plot_type='deformation,displacement')