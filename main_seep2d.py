from seep import import_seep2d, run_seepage_analysis, print_seep_data_diagnostics, export_solution_csv
from plot_seep import plot_seep_mesh, plot_seep_solution
import numpy as np


# Load input
seep_data = import_seep2d("inputs/seep/lface.s2d")

# Print diagnostics
print_seep_data_diagnostics(seep_data)

# Plot mesh
plot_seep_mesh(seep_data, show_nodes=True, show_bc=True)

# Run analysis
solution = run_seepage_analysis(seep_data)

# # Export solution to CSV
# export_solution_csv("seep_solution.csv", seep_data, solution)

# # Plot solution
plot_seep_solution(seep_data, solution, base_mat=2, levels=30, fill_contours=True, phreatic=True) 