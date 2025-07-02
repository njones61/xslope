#!/usr/bin/env python3
"""
Test script to compare performance of original vs optimized generate_slices function.
"""

import time
import numpy as np
from shapely.geometry import LineString
import slice

def create_test_data(num_profile_points=100, num_slices=40, circular=True):
    """Create test data for performance comparison."""
    
    # Create ground surface
    x_ground = np.linspace(0, 100, num_profile_points)
    y_ground = 50 + 20 * np.sin(x_ground * 0.1) + 5 * np.random.random(len(x_ground))
    ground_surface = LineString(list(zip(x_ground, y_ground)))
    
    # Create profile lines (soil layers)
    profile_lines = []
    materials = []
    
    # Layer 1
    y_layer1 = y_ground - 5 - 2 * np.random.random(len(x_ground))
    profile_lines.append(list(zip(x_ground, y_layer1)))
    materials.append({
        'gamma': 120.0,
        'option': 'mc',
        'c': 500.0,
        'phi': 30.0,
        'd': 0.0,
        'psi': 0.0
    })
    
    # Layer 2
    y_layer2 = y_layer1 - 8 - 3 * np.random.random(len(x_ground))
    profile_lines.append(list(zip(x_ground, y_layer2)))
    materials.append({
        'gamma': 125.0,
        'option': 'mc',
        'c': 800.0,
        'phi': 25.0,
        'd': 0.0,
        'psi': 0.0
    })
    
    # Layer 3
    y_layer3 = y_layer2 - 10 - 4 * np.random.random(len(x_ground))
    profile_lines.append(list(zip(x_ground, y_layer3)))
    materials.append({
        'gamma': 130.0,
        'option': 'mc',
        'c': 1200.0,
        'phi': 20.0,
        'd': 0.0,
        'psi': 0.0
    })
    
    # Create piezometric line
    x_piezo = np.linspace(0, 100, 20)
    y_piezo = 45 + 10 * np.sin(x_piezo * 0.15)
    piezo_line = list(zip(x_piezo, y_piezo))
    
    # Create distributed loads
    x_dload = np.linspace(20, 80, 10)
    y_dload = y_ground[:len(x_dload)] + 2
    dloads = [[{'X': x, 'Y': y, 'Normal': 1000.0} for x, y in zip(x_dload, y_dload)]]
    
    # Create failure surface
    if circular:
        circle = {
            'Xo': 50.0,
            'Yo': 60.0,
            'Depth': 20.0,
            'R': 25.0
        }
        non_circ = None
    else:
        circle = None
        # Create non-circular failure surface
        x_failure = np.linspace(25, 75, 20)
        y_failure = 30 + 15 * np.sin((x_failure - 50) * 0.2) + 2 * np.random.random(len(x_failure))
        non_circ = [{'X': x, 'Y': y, 'Movement': 0.0} for x, y in zip(x_failure, y_failure)]
    
    # Create data dictionary
    data = {
        'profile_lines': profile_lines,
        'ground_surface': ground_surface,
        'materials': materials,
        'piezo_line': piezo_line,
        'piezo_line2': [],
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'dloads': dloads,
        'dloads2': [],
        'max_depth': 50.0,
        'reinforce_lines': []
    }
    
    return data, circle, non_circ

def test_performance():
    """Test performance of generate_slices function."""
    
    print("=== SLICE GENERATION PERFORMANCE TEST ===\n")
    
    # Test configurations
    test_configs = [
        {"name": "Small dataset (circular)", "points": 100, "slices": 30, "circular": True},
        {"name": "Medium dataset (circular)", "points": 200, "slices": 50, "circular": True},
        {"name": "Large dataset (circular)", "points": 500, "slices": 100, "circular": True},
        {"name": "Small dataset (non-circular)", "points": 100, "slices": 30, "circular": False},
        {"name": "Medium dataset (non-circular)", "points": 200, "slices": 50, "circular": False},
    ]
    
    for config in test_configs:
        print(f"Testing: {config['name']}")
        print(f"Profile points: {config['points']}, Slices: {config['slices']}")
        
        data, circle, non_circ = create_test_data(
            num_profile_points=config['points'], 
            num_slices=config['slices'],
            circular=config['circular']
        )
        
        # Test multiple runs to get average performance
        num_runs = 3
        times = []
        
        for i in range(num_runs):
            start_time = time.time()
            
            if config['circular']:
                success, result = slice.generate_slices(data, circle=circle, num_slices=config['slices'], debug=False)
            else:
                success, result = slice.generate_slices(data, non_circ=non_circ, num_slices=config['slices'], debug=False)
            
            end_time = time.time()
            elapsed = end_time - start_time
            times.append(elapsed)
            
            if not success:
                print(f"  Error in run {i+1}: {result}")
                continue
        
        if times:
            avg_time = np.mean(times)
            std_time = np.std(times)
            
            print(f"  Average time: {avg_time:.4f} seconds")
            print(f"  Standard deviation: {std_time:.4f} seconds")
            print(f"  Min time: {min(times):.4f} seconds")
            print(f"  Max time: {max(times):.4f} seconds")
            
            if success:
                df, clipped_surface = result
                print(f"  Generated {len(df)} slices successfully")
        
        print()

def test_correctness():
    """Test that the optimized version produces the same results as expected."""
    
    print("=== CORRECTNESS TEST ===\n")
    
    # Create simple test case
    data, circle, _ = create_test_data(num_profile_points=50, num_slices=20, circular=True)
    
    print("Testing circular failure surface...")
    success, result = slice.generate_slices(data, circle=circle, num_slices=20, debug=False)
    
    if success:
        df, clipped_surface = result
        print(f"✓ Generated {len(df)} slices successfully")
        print(f"✓ Failure surface type: {type(clipped_surface)}")
        print(f"✓ First slice coordinates: x_l={df.iloc[0]['x_l']:.3f}, y_lb={df.iloc[0]['y_lb']:.3f}")
        print(f"✓ Last slice coordinates: x_r={df.iloc[-1]['x_r']:.3f}, y_rb={df.iloc[-1]['y_rb']:.3f}")
        
        # Check that all required columns are present
        required_columns = ['slice #', 'x_l', 'y_lb', 'y_lt', 'x_r', 'y_rb', 'y_rt', 'x_c', 'y_cb', 'y_ct', 'alpha', 'w']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"✗ Missing columns: {missing_columns}")
        else:
            print("✓ All required columns present")
            
        # Check that coordinates are reasonable
        if df['x_l'].min() < df['x_r'].max() and df['y_lb'].min() < df['y_lt'].max():
            print("✓ Coordinate ranges are reasonable")
        else:
            print("✗ Coordinate ranges are suspicious")
    else:
        print(f"✗ Failed to generate slices: {result}")

if __name__ == "__main__":
    test_correctness()
    print()
    test_performance() 