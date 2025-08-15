#!/usr/bin/env python3

from fem import build_fem_data, solve_fem, solve_ssrm
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import signal
import sys

# Timeout handler
def timeout_handler(signum, frame):
    print("\nTimeout!")
    sys.exit(1)

signal.signal(signal.SIGALRM, timeout_handler)

# Load data with coarser mesh
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 4, 'tri3')  # Coarser mesh
fem_data = build_fem_data(slope_data, mesh)

print(f"Testing final SSRM with aggressive plastic reduction (0.001) and 5% displacement threshold")
print(f"Mesh: {len(fem_data['nodes'])} nodes")

# Set timeout
signal.alarm(120)  # 2 minutes

print("\nRunning SSRM...")
ssrm_result = solve_ssrm(fem_data, debug_level=2)

signal.alarm(0)

# Print results
if ssrm_result.get("converged", False):
    print(f"\n✓ SSRM converged successfully")
    print(f"Factor of Safety: {ssrm_result['FS']:.4f}")
    print(f"Expected: ~1.28 from limit equilibrium")
    print(f"SSRM iterations: {ssrm_result['iterations_ssrm']}")
    
    if 'failure_mode' in ssrm_result:
        print(f"Failure mode: {ssrm_result['failure_mode']}")
        
else:
    print(f"\n✗ SSRM failed to converge")
    if 'error' in ssrm_result:
        print(f"Error: {ssrm_result['error']}")
    if 'FS' in ssrm_result and ssrm_result['FS'] is not None:
        print(f"Last converged FS: {ssrm_result['FS']:.4f}")