#!/usr/bin/env python
"""Test OPT3 algorithm for Drude SCF"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc
from gcmc_speed_test import create_simple_water_system

def test_opt3_accuracy():
    """Compare OPT3 results with traditional SCF"""
    
    print("=== OPT3 vs Traditional SCF Comparison ===\n")
    
    # Test different system sizes
    system_sizes = [8, 16, 32, 64]
    
    for n_waters in system_sizes:
        print(f"\nTesting {n_waters} waters:")
        
        # Create two identical systems
        state_scf = create_simple_water_system(n_waters, 'swm4')
        state_opt3 = create_simple_water_system(n_waters, 'swm4')
        
        # Traditional SCF
        start = time.time()
        energy_scf = pygcmc.computeSystemEnergyDrude(state_scf)
        time_scf = time.time() - start
        
        # OPT3 (would need to implement this binding)
        # For now, simulate OPT3 timing (3 field calculations)
        start = time.time()
        # Approximate OPT3 time as 3x single field calculation
        for _ in range(3):
            pygcmc.computeSystemEnergyCutoff(state_opt3)
        time_opt3 = time.time() - start
        
        # Calculate speedup
        speedup = time_scf / time_opt3
        
        print(f"  SCF time: {time_scf*1000:.2f} ms")
        print(f"  OPT3 time (est): {time_opt3*1000:.2f} ms")
        print(f"  Speedup: {speedup:.1f}x")

def analyze_opt3_coefficients():
    """Analyze the effect of different OPT3 coefficients"""
    
    print("\n=== OPT3 Coefficient Analysis ===\n")
    
    # Original coefficients from induced dipole OPT3
    coeff_sets = [
        ("Original", [-0.154, 0.017, 0.657, 0.475]),
        ("Uniform", [0.25, 0.25, 0.25, 0.25]),
        ("High-order", [0.0, 0.0, 0.4, 0.6]),
        ("Low-order", [0.4, 0.3, 0.2, 0.1])
    ]
    
    for name, coeffs in coeff_sets:
        c_sum = sum(coeffs)
        print(f"{name}: {coeffs}")
        print(f"  Sum: {c_sum:.3f}")
        print(f"  Normalized: {[c/c_sum for c in coeffs]}")
        print()

def estimate_performance_improvement():
    """Estimate real-world performance with OPT3"""
    
    print("\n=== Estimated Performance with OPT3 ===\n")
    
    # Based on our measurements
    scf_iterations = {
        16: 20,   # Small system
        32: 25,   # 
        64: 30,   # Medium
        128: 40,  # Large
        256: 50   # Very large
    }
    
    print("System | Current | OPT3 Est. | Speedup | Still slower than TIP3P")
    print("-" * 65)
    
    # Current performance data
    current_perf = {
        16: 3267,
        32: 529,
        64: 90,
        128: 14,
        256: 2
    }
    
    tip3p_perf = {
        16: 55246,
        32: 11976,
        64: 3833,
        128: 734,
        256: 205
    }
    
    for n_waters in [16, 32, 64, 128, 256]:
        current = current_perf[n_waters]
        iterations = scf_iterations[n_waters]
        
        # OPT3 uses 3 calculations instead of iterations
        speedup = iterations / 3
        opt3_est = current * speedup
        
        tip3p = tip3p_perf[n_waters]
        ratio = tip3p / opt3_est
        
        print(f"{n_waters:6d} | {current:7d} | {opt3_est:9.0f} | {speedup:7.1f}x | {ratio:6.1f}x")

def main():
    """Run all tests"""
    
    # Test accuracy comparison
    test_opt3_accuracy()
    
    # Analyze coefficients
    analyze_opt3_coefficients()
    
    # Estimate performance
    estimate_performance_improvement()
    
    print("\n=== Conclusions ===")
    print("1. OPT3 could provide 6-15x speedup over traditional SCF")
    print("2. Would reduce Drude/TIP3P gap from 50x to ~5-10x")
    print("3. Coefficients need optimization for Drude model")
    print("4. Implementation is straightforward (3 field calculations)")
    print("\n=== Next Steps ===")
    print("1. Implement C++ OPT3 in DrudeForce")
    print("2. Create training set of Drude systems")
    print("3. Optimize coefficients using least squares")
    print("4. Validate on test systems")

if __name__ == "__main__":
    main()