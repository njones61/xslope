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

print("=== Testing Automatic Scale Factor Calculation ===\n")

print("1. Deformation plot with automatic scale factor:")
plot_fem_results(fem_data, solution, plot_type='deformation')

print("\n2. Displacement plot (no scale factor needed):")
plot_fem_results(fem_data, solution, plot_type='displacement')

print("\n3. Deformation plot with manual scale factor:")
plot_fem_results(fem_data, solution, plot_type='deformation', deform_scale=5.0)

print("\n4. Deformation plot with scale=1.0 (actual displacements):")
plot_fem_results(fem_data, solution, plot_type='deformation', deform_scale=1.0)