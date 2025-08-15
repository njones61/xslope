#!/usr/bin/env python3

from fem import build_fem_data, solve_fem, solve_ssrm
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import sys

# Load data and build FEM model
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 2, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

print("Testing single FEM solution first...")
solution = solve_fem(fem_data, F=1.0, debug_level=1)
print(f"FEM converged: {solution.get('converged', False)}")
if not solution.get('converged', False):
    print("FEM not converging - exiting")
    sys.exit(1)

print("\nTesting higher F values...")
for F in [1.1, 1.2, 1.3, 1.4, 1.5]:
    solution = solve_fem(fem_data, F=F, debug_level=1)
    print(f"F={F:.1f}: converged={solution.get('converged', False)}")
    if not solution.get('converged', False):
        print(f"First failure at F={F:.1f}")
        break

print("\nTesting SSRM with max 5 iterations...")
# Run SSRM but add some safeguards
try:
    ssrm_result = solve_ssrm(fem_data, debug_level=2)
    print(f"SSRM result: {ssrm_result}")
except Exception as e:
    print(f"SSRM failed with error: {e}")