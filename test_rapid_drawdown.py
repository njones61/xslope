#!/usr/bin/env python3
"""
Test script for the rapid_drawdown function
"""

import pandas as pd
import numpy as np
from solve import rapid_drawdown, oms

def test_rapid_drawdown():
    """Test the rapid_drawdown function with a simple example"""
    
    # Create a simple test dataframe with 3 slices
    data = {
        'slice #': [1, 2, 3],
        'alpha': [10, 15, 20],  # degrees
        'phi': [30, 25, 35],    # degrees
        'c': [500, 400, 600],   # psf
        'w': [1000, 1200, 800], # lbs/ft
        'u': [100, 150, 200],   # psf
        'dl': [10, 12, 8],      # ft
        'dload': [0, 0, 0],     # lbs/ft
        'd_x': [0, 0, 0],       # ft
        'd_y': [0, 0, 0],       # ft
        'beta': [0, 0, 0],      # degrees
        'kw': [0, 0, 0],        # lbs
        't': [0, 0, 0],         # lbs
        'y_t': [0, 0, 0],       # ft
        'p': [0, 0, 0],         # lbs
        'x_c': [5, 15, 25],     # ft
        'y_cg': [10, 12, 8],    # ft
        'r': [50, 50, 50],      # ft
        'xo': [25, 25, 25],     # ft
        'yo': [30, 30, 30],     # ft
        # Rapid drawdown specific data
        'c1': [500, 400, 600],  # Original c values
        'phi1': [30, 25, 35],   # Original phi values
        'd': [0, 200, 0],       # d values (only slice 2 has low-K)
        'psi': [0, 15, 0],      # psi values (only slice 2 has low-K)
        'u2': [50, 75, 100],    # Stage 2 pore pressures
        'dload2': [0, 0, 0],    # Stage 2 distributed loads
        'd_x2': [0, 0, 0],      # Stage 2 d_x
        'd_y2': [0, 0, 0],      # Stage 2 d_y
    }
    
    df = pd.DataFrame(data)
    
    # Ensure numeric columns have float dtype to avoid warnings
    numeric_columns = ['c', 'phi', 'c1', 'phi1', 'd', 'psi', 'u', 'u2', 'dload', 'dload2', 'd_x', 'd_x2', 'd_y', 'd_y2']
    for col in numeric_columns:
        if col in df.columns:
            df[col] = df[col].astype(float)
    
    print("Testing rapid_drawdown function...")
    print("Input data:")
    print(df[['slice #', 'c', 'phi', 'd', 'psi', 'u', 'u2']])
    
    # Test with debug_level=1 to see the stages
    success, result = rapid_drawdown(df, oms, debug_level=1)
    
    if success:
        print("\nRapid drawdown analysis successful!")
        print(f"Final FS: {result['FS']:.4f}")
        print(f"Stage 1 FS: {result['stage1_FS']:.4f}")
        print(f"Stage 2 FS: {result['stage2_FS']:.4f}")
        print(f"Stage 3 FS: {result['stage3_FS']:.4f}")
    else:
        print(f"\nRapid drawdown analysis failed: {result}")
    
    return success, result

if __name__ == "__main__":
    test_rapid_drawdown() 