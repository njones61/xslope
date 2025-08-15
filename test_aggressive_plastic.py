#!/usr/bin/env python3

from fem import build_fem_data, solve_fem
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

print(f"Testing with aggressive plastic stiffness reduction (0.01)")
print(f"Mesh: {len(fem_data['nodes'])} nodes")

# Test individual F values with new default
signal.alarm(120)  # 2 minute timeout

print("\nTesting individual F values:")
for F in [1.0, 1.1, 1.2, 1.25, 1.3, 1.35, 1.4]:
    solution = solve_fem(fem_data, F=F, debug_level=1)
    converged = solution.get('converged', False)
    iterations = solution.get('iterations', 'unknown')
    max_disp = max(abs(solution.get('displacements', [0])))
    print(f"F={F:.2f}: converged={converged}, iterations={iterations}, max_disp={max_disp:.4f}")
    
    if not converged:
        print(f"  First failure at F={F:.2f}")
        break

signal.alarm(0)