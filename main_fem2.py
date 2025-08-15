
from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np
from plot import plot_inputs, plot_mesh

slope_data = load_slope_data("inputs/slope/input_template_lface6.xlsx")

plot_inputs(slope_data)

polygons = build_polygons(slope_data)

target_size = 5

mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')

plot_mesh(mesh, materials=slope_data['materials'])