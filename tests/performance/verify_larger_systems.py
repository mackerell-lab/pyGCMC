#!/usr/bin/env python
"""Verify performance trends with larger systems"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')
import pygcmc

# Import the create function from existing script
from gcmc_speed_test import create_simple_water_system

def measure_single_system(n_waters, n_iterations=50):
    """Measure performance for a single system size with more iterations"""
    
    results = {}
    
    for model in ['tip3p', 'swm4']:
        # Create system
        state = create_simple_water_system(n_waters, model)
        
        # Warm up
        if model == 'tip3p':
            pygcmc.computeSystemEnergyCutoff(state)
        else:
            pygcmc.computeSystemEnergyDrude(state)
        
        # Measure with more iterations for accuracy
        times = []
        for _ in range(5):  # Multiple measurements
            start = time.time()
            for _ in range(n_iterations):
                if model == 'tip3p':
                    pygcmc.computeSystemEnergyCutoff(state)
                else:
                    pygcmc.computeSystemEnergyDrude(state)
            elapsed = time.time() - start
            times.append(elapsed)
        
        # Use median to reduce outlier effect
        median_time = np.median(times)
        calcs_per_sec = n_iterations / median_time
        gcmc_steps_per_sec = calcs_per_sec / 2  # Assume 2 calcs per move
        
        results[model] = {
            'calcs_per_sec': calcs_per_sec,
            'gcmc_steps_per_sec': gcmc_steps_per_sec,
            'time_per_calc_ms': median_time / n_iterations * 1000
        }
    
    return results

def main():
    """Test larger systems and verify trends"""
    
    print("=== Verifying Performance Trends ===\n")
    
    # Test systems: skip 8 (too small), use larger sizes
    system_sizes = [16, 32, 64, 128, 256]
    
    all_results = []
    
    for n_waters in system_sizes:
        print(f"\nTesting {n_waters} waters...")
        
        # Skip very large if too slow
        if n_waters > 128:
            n_iter = 10
        elif n_waters > 64:
            n_iter = 25
        else:
            n_iter = 50
        
        results = measure_single_system(n_waters, n_iter)
        
        tip3p_speed = results['tip3p']['gcmc_steps_per_sec']
        swm4_speed = results['swm4']['gcmc_steps_per_sec']
        ratio = tip3p_speed / swm4_speed
        
        all_results.append({
            'n_waters': n_waters,
            'tip3p': tip3p_speed,
            'swm4': swm4_speed,
            'ratio': ratio
        })
        
        print(f"  TIP3P: {tip3p_speed:8.0f} steps/s")
        print(f"  SWM4:  {swm4_speed:8.0f} steps/s")
        print(f"  Ratio: {ratio:6.1f}x")
    
    # Analyze monotonicity
    print("\n=== Trend Analysis ===")
    print("N_waters | TIP3P steps/s | SWM4 steps/s | Speed Ratio")
    print("-" * 55)
    
    for r in all_results:
        print(f"{r['n_waters']:8d} | {r['tip3p']:13.0f} | {r['swm4']:12.0f} | {r['ratio']:11.1f}x")
    
    # Check if ratios are monotonically increasing
    ratios = [r['ratio'] for r in all_results]
    is_monotonic = all(ratios[i] <= ratios[i+1] for i in range(len(ratios)-1))
    
    print(f"\nRatios monotonically increasing: {is_monotonic}")
    
    # Calculate scaling exponents
    if len(all_results) >= 3:
        n_vals = np.array([r['n_waters'] for r in all_results])
        tip3p_vals = np.array([r['tip3p'] for r in all_results])
        swm4_vals = np.array([r['swm4'] for r in all_results])
        
        # Log-log fit for scaling
        p_tip3p = np.polyfit(np.log(n_vals), np.log(tip3p_vals), 1)
        p_swm4 = np.polyfit(np.log(n_vals), np.log(swm4_vals), 1)
        
        print(f"\nScaling exponents:")
        print(f"TIP3P: N^{-p_tip3p[0]:.2f}")
        print(f"SWM4:  N^{-p_swm4[0]:.2f}")
        print(f"Difference: {abs(-p_swm4[0] - (-p_tip3p[0])):.2f}")
    
    # Generate updated data for PPT
    print("\n=== Updated PPT Data ===")
    print("# Use these values for the presentation:")
    for r in all_results[:4]:  # First 4 for PPT
        print(f"{r['n_waters']}, {r['tip3p']:.0f}, {r['swm4']:.0f}")

if __name__ == "__main__":
    main()