#!/usr/bin/env python3

from fem import build_fem_data, solve_fem
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import signal
import sys
import time

# Timeout handler
def timeout_handler(signum, frame):
    print("\nTimeout! FEM solver is hanging.")
    sys.exit(1)

# Set 30 second timeout
signal.signal(signal.SIGALRM, timeout_handler)

# Load minimal data
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)

# Use coarser mesh to speed up testing
mesh = build_mesh_from_polygons(polygons, 4, 'tri3')  # target_size=4 instead of 2
fem_data = build_fem_data(slope_data, mesh)

print(f"Testing coarser mesh: {len(fem_data['nodes'])} nodes")

# Start timeout
signal.alarm(30)

print("Testing F=1.0...")
start_time = time.time()
solution = solve_fem(fem_data, F=1.0, debug_level=2)
elapsed = time.time() - start_time

# Cancel timeout
signal.alarm(0)

print(f"\nSolution completed in {elapsed:.1f} seconds")
print(f"Converged: {solution.get('converged', False)}")
print(f"Iterations: {solution.get('iterations', 'unknown')}")

if solution.get('converged', False):
    max_disp = max(abs(solution['displacements']))
    print(f"Max displacement: {max_disp:.3f}")
else:
    print(f"Error: {solution.get('error', 'unknown')}")