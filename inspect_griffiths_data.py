#!/usr/bin/env python3
"""
Inspect Griffiths Example 1 data structure.
"""

from fileio import load_slope_data, print_dictionary
import numpy as np

def inspect_griffiths_data():
    """Inspect the loaded data structure."""
    
    print("=== Inspecting Griffiths Example 1 Data ===")
    
    # Load data
    slope_data = load_slope_data("inputs/slope/input_template_griffiths1_6.xlsx")
    
    # Print full structure
    print("\n=== Full Data Structure ===")
    print_dictionary(slope_data)
    
    # Check materials specifically
    print(f"\n=== Materials ===")
    materials = slope_data.get('materials', [])
    print(f"Number of materials: {len(materials)}")
    for i, mat in enumerate(materials):
        print(f"Material {i}: {mat}")
    
    # Check ground surface
    print(f"\n=== Ground Surface ===")
    ground_surface = slope_data.get('ground_surface', [])
    print(f"Ground surface points: {len(ground_surface)}")
    if len(ground_surface) > 0:
        print(f"First few points: {ground_surface[:3]}")
        y_values = ground_surface[:, 1] if hasattr(ground_surface, 'shape') else [p[1] for p in ground_surface]
        print(f"Y range: {min(y_values):.1f} to {max(y_values):.1f}")
    
    # Check circles
    print(f"\n=== Circles ===")
    circles = slope_data.get('circles', [])
    print(f"Number of circles: {len(circles)}")
    for i, circle in enumerate(circles):
        print(f"Circle {i}: {circle}")

if __name__ == "__main__":
    inspect_griffiths_data()