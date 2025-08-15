#!/usr/bin/env python3

from fem import build_fem_data, solve_fem
from fileio import load_slope_data
from mesh import build_polygons, build_mesh_from_polygons
import numpy as np

# Load data and build FEM model
slope_data = load_slope_data("inputs/slope/input_template_simple1_6.xlsx")
polygons = build_polygons(slope_data)
target_size = 2
mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

# Run FEM analysis
solution = solve_fem(fem_data, F=1.0, debug_level=0)

# Extract data
nodes = fem_data['nodes']
bc_type = fem_data['bc_type']
displacements = solution['displacements']

# Extract displacements
u = displacements[0::2]  # x-displacements
v = displacements[1::2]  # y-displacements

print("=== Boundary Condition Verification ===")

# Check left boundary (x=0) nodes
left_boundary_indices = []
tolerance = 1e-6
for i, node in enumerate(nodes):
    if abs(node[0] - 0.0) < tolerance:  # x = 0
        left_boundary_indices.append(i)

print(f"Found {len(left_boundary_indices)} nodes on left boundary (x=0)")

# Check their boundary conditions and displacements
for i in left_boundary_indices:
    x, y = nodes[i]
    bc = bc_type[i]
    u_val = u[i]
    v_val = v[i]
    print(f"Node {i}: ({x:.6f}, {y:.6f}) BC_type={bc} u={u_val:.6f} v={v_val:.6f}")

# Check right boundary (x=100) nodes  
right_boundary_indices = []
for i, node in enumerate(nodes):
    if abs(node[0] - 100.0) < tolerance:  # x = 100
        right_boundary_indices.append(i)

print(f"\nFound {len(right_boundary_indices)} nodes on right boundary (x=100)")

# Check their boundary conditions and displacements
for i in right_boundary_indices:
    x, y = nodes[i]
    bc = bc_type[i]
    u_val = u[i]
    v_val = v[i]
    print(f"Node {i}: ({x:.6f}, {y:.6f}) BC_type={bc} u={u_val:.6f} v={v_val:.6f}")

# Summary
left_u_values = [u[i] for i in left_boundary_indices]
right_u_values = [u[i] for i in right_boundary_indices]

print(f"\n=== Summary ===")
print(f"Left boundary X-displacements: min={min(left_u_values):.6f}, max={max(left_u_values):.6f}")
print(f"Right boundary X-displacements: min={min(right_u_values):.6f}, max={max(right_u_values):.6f}")

# Verify boundary conditions are enforced
max_left_u = max(abs(val) for val in left_u_values)
max_right_u = max(abs(val) for val in right_u_values)

if max_left_u < 1e-10 and max_right_u < 1e-10:
    print("✓ PASS: X-roller boundary conditions correctly enforced")
else:
    print("✗ FAIL: X-roller boundary conditions NOT properly enforced")
    print(f"  Max left boundary u: {max_left_u:.2e}")
    print(f"  Max right boundary u: {max_right_u:.2e}")