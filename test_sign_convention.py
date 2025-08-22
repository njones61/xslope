#!/usr/bin/env python3
"""
Test that the sign convention changes are working correctly.
In geotechnical engineering: compression positive, tension negative.
"""

import numpy as np
from fem import check_mohr_coulomb
from math import radians, degrees

def test_mohr_coulomb_compression():
    """Test Mohr-Coulomb criterion with compression-positive convention."""
    
    print("="*60)
    print("Testing Mohr-Coulomb yield function with compression-positive")
    print("Using formula: f = tau_max - (sig_mean * sin(φ) + c * cos(φ))")
    print("="*60)
    
    # Test case 1: Pure compression (should be stable for reasonable values)
    print("\nTest 1: Pure hydrostatic compression")
    stress = np.array([100.0, 100.0, 0.0])  # Hydrostatic compression (positive)
    c = 10.0  # cohesion
    phi = radians(30)  # friction angle
    
    F = check_mohr_coulomb(stress, c, phi)
    print(f"  Stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
    print(f"  c={c:.1f}, φ={degrees(phi):.1f}°")
    print(f"  Principal stresses: σ1=100.0, σ3=100.0")
    print(f"  tau_max=0.0, sig_mean=100.0")
    print(f"  F={F:.3f} (negative = stable, positive = yielding)")
    
    # Test case 2: Compression with shear
    print("\nTest 2: Compression with shear")
    stress = np.array([100.0, 50.0, 20.0])  # Compression with shear
    F = check_mohr_coulomb(stress, c, phi, debug=True)
    print(f"  Stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
    print(f"  F={F:.3f}")
    
    # Test case 3: Near failure
    print("\nTest 3: Near failure state")
    stress = np.array([30.0, 10.0, 15.0])  # Lower compression, higher shear
    F = check_mohr_coulomb(stress, c, phi, debug=True)
    print(f"  Stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
    print(f"  F={F:.3f}")
    
    # Test case 4: Tension (should yield)
    print("\nTest 4: Tension state (should yield)")
    stress = np.array([-10.0, -10.0, 0.0])  # Tension (negative values)
    F = check_mohr_coulomb(stress, c, phi, debug=True)
    print(f"  Stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
    print(f"  F={F:.3f}")
    
    # Test case 5: Failure case (should be positive F)
    print("\nTest 5: Known failure state")
    stress = np.array([20.0, 5.0, 25.0])  # High shear relative to normal stress
    F = check_mohr_coulomb(stress, c, phi, debug=True)
    print(f"  Stress: σx={stress[0]:.1f}, σy={stress[1]:.1f}, τxy={stress[2]:.1f}")
    print(f"  F={F:.3f}")
    
    # Verify sign convention
    print("\n" + "="*60)
    print("Sign Convention Summary:")
    print("  Compression: positive values")
    print("  Tension: negative values")
    print("  F > 0: Material is yielding (failure)")
    print("  F < 0: Material is stable")
    print("="*60)

if __name__ == "__main__":
    test_mohr_coulomb_compression()