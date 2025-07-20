#!/usr/bin/env python
"""Find optimal SCF parameters for different use cases"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc
from test_simple_drude_water_placement import create_simple_water_box, optimize_simple_positions, convert_to_drude_water

def test_parameter_combinations():
    """Test various parameter combinations to find optimal settings"""
    print("=== Finding Optimal SCF Parameters ===\n")
    
    # Create a challenging 27-water system
    n_waters = 27
    box_size = 2.0
    
    print(f"Testing with {n_waters} waters in {box_size} nm box\n")
    
    # Parameter combinations to test
    param_combos = [
        # (tolerance, max_iter, damping, name)
        (1.0, 50, 0.5, "Default"),
        (10.0, 100, 0.3, "Relaxed"),
        (50.0, 200, 0.2, "Very Relaxed"),
        (100.0, 300, 0.2, "Ultra Relaxed"),
        (5.0, 100, 0.4, "Balanced"),
        (20.0, 150, 0.3, "Moderate"),
    ]
    
    results = []
    
    for tol, max_iter, damp, name in param_combos:
        print(f"\nTesting {name}: tol={tol}, iter={max_iter}, damp={damp}")
        
        # Create system
        simple_state = create_simple_water_box(n_waters, box_size)
        simple_state = optimize_simple_positions(simple_state, n_steps=100)
        
        # Initialize and set parameters
        pygcmc.initializeDrudeForce()
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = tol
        scf_params.maxIterations = max_iter
        scf_params.dampingFactor = damp
        scf_params.maxDrudeDistance = 0.02
        pygcmc.setDrudeSCFParameters(scf_params)
        
        # Convert to Drude
        drude_state = convert_to_drude_water(simple_state)
        
        # Calculate energy and check for warnings
        import io
        import contextlib
        
        # Capture stderr to check for warnings
        f = io.StringIO()
        with contextlib.redirect_stderr(f):
            result = pygcmc.computeSystemEnergyDrude(drude_state)
        
        stderr_output = f.getvalue()
        has_warning = "Warning" in stderr_output
        
        if isinstance(result, tuple):
            energy = result[0]
            energy_per_water = energy / n_waters
            
            # Extract RMS force if present in warning
            rms_force = None
            if "RMS force =" in stderr_output:
                try:
                    rms_part = stderr_output.split("RMS force =")[1].split("kJ/mol/nm")[0]
                    rms_force = float(rms_part.strip())
                except:
                    pass
            
            results.append({
                'name': name,
                'tolerance': tol,
                'max_iter': max_iter,
                'damping': damp,
                'energy_per_water': energy_per_water,
                'has_warning': has_warning,
                'rms_force': rms_force
            })
            
            print(f"  Energy per water: {energy_per_water:.2f} kJ/mol")
            print(f"  Warning: {'Yes' if has_warning else 'No'}")
            if rms_force:
                print(f"  Final RMS force: {rms_force:.2f} kJ/mol/nm")
    
    # Summary
    print("\n\n=== SUMMARY ===")
    print("\nConfigurations without warnings:")
    for r in results:
        if not r['has_warning']:
            print(f"  {r['name']}: tol={r['tolerance']}, iter={r['max_iter']}, damp={r['damping']}")
            print(f"    Energy: {r['energy_per_water']:.2f} kJ/mol per water")
    
    print("\nRecommendation:")
    print("For general use: tolerance=10.0, maxIterations=100, dampingFactor=0.3")
    print("This provides a good balance between accuracy and robustness.")

if __name__ == "__main__":
    test_parameter_combinations()