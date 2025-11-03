#!/usr/bin/env python3
"""Diagnostic to understand premature yielding."""

from xslope.fem import solve_fem, build_fem_data
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons
import numpy as np

# Load slope
slope_data = load_slope_data('inputs/slope/input_template_griffiths1_6.xlsx')
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, 5, 'tri3')
fem_data = build_fem_data(slope_data, mesh)

# Get material properties
mat = slope_data['materials'][0]
c_orig = mat['c']
phi_orig = mat['phi']
gamma = mat['gamma']

print(f"Original material properties:")
print(f"  c = {c_orig:.1f} kPa")
print(f"  φ = {phi_orig:.1f}°") 
print(f"  γ = {gamma:.1f} kN/m³")

# Run with F=1.0
print(f"\n=== Testing with F=1.0 (no reduction) ===")
solution = solve_fem(fem_data, F=1.0, debug_level=1, abort_after=0)

# Check a few elements
yield_fn = solution.get("yield_function", [])
stresses = solution.get("stresses", [])

print(f"\n=== Sample element analysis ===")
for i in [0, 50, 100, 150, 200]:
    if i < len(yield_fn):
        F_val = yield_fn[i]
        sig_x, sig_y, tau_xy = stresses[i][:3]
        sig_mean = (sig_x + sig_y) / 2
        tau_max = np.sqrt(((sig_x - sig_y)/2)**2 + tau_xy**2)
        
        print(f"\nElement {i}:")
        print(f"  Stresses (compression+): σx={sig_x:.1f}, σy={sig_y:.1f}, τxy={tau_xy:.1f}")
        print(f"  σ_mean={sig_mean:.1f}, τ_max={tau_max:.1f}")
        print(f"  Yield function F={F_val:.1f} {'(YIELDING!)' if F_val > 0 else '(safe)'}")
        
        # Manual check
        phi_rad = np.radians(phi_orig)
        F_manual = tau_max - sig_mean * np.sin(phi_rad) - c_orig * np.cos(phi_rad)
        print(f"  Manual F check: {F_manual:.1f}")
