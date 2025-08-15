#!/usr/bin/env python3

from fem import build_fem_data, solve_fem, solve_ssrm
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import signal
import sys
import time

# Timeout handler
def timeout_handler(signum, frame):
    print("\nTimeout! SSRM is hanging.")
    sys.exit(1)

# Set 2 minute timeout
signal.signal(signal.SIGALRM, timeout_handler)

# Load data with coarser mesh
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 4, 'tri3')  # Coarser mesh
fem_data = build_fem_data(slope_data, mesh)

print(f"Testing SSRM with coarser mesh: {len(fem_data['nodes'])} nodes")

# Test individual F values first
print("\nTesting individual F values:")
for F in [1.0, 1.1, 1.2, 1.25, 1.3]:
    start_time = time.time()
    solution = solve_fem(fem_data, F=F, debug_level=1)
    elapsed = time.time() - start_time
    converged = solution.get('converged', False)
    iterations = solution.get('iterations', 'unknown')
    print(f"F={F:.2f}: converged={converged}, iterations={iterations}, time={elapsed:.1f}s")
    
    if not converged:
        print(f"  First failure at F={F:.2f}")
        break

print("\nNow testing full SSRM...")

# Start timeout
signal.alarm(120)  # 2 minutes

start_time = time.time()
ssrm_result = solve_ssrm(fem_data, debug_level=2)
elapsed = time.time() - start_time

# Cancel timeout
signal.alarm(0)

print(f"\nSSRM completed in {elapsed:.1f} seconds")

if ssrm_result.get("converged", False):
    print(f"✓ SSRM converged successfully")
    print(f"Factor of Safety: {ssrm_result['FS']:.4f}")
    print(f"Expected: ~1.28 from limit equilibrium")
    print(f"SSRM iterations: {ssrm_result['iterations_ssrm']}")
else:
    print(f"✗ SSRM failed to converge")
    if 'error' in ssrm_result:
        print(f"Error: {ssrm_result['error']}")
    if 'FS' in ssrm_result and ssrm_result['FS'] is not None:
        print(f"Last converged FS: {ssrm_result['FS']:.4f}")