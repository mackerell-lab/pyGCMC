#!/usr/bin/env python
"""Optimized SCF parameters for large Drude systems"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def get_optimized_scf_params(system_size='small'):
    """Get optimized SCF parameters based on system size
    
    Args:
        system_size: 'small' (<10 waters), 'medium' (10-50), 'large' (>50)
    
    Returns:
        DrudeSCFParams object with optimized settings
    """
    params = pygcmc.DrudeSCFParams()
    
    if system_size == 'small':
        # Tight convergence for small systems
        params.tolerance = 1.0  # kJ/mol/nm
        params.maxIterations = 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02  # nm
        
    elif system_size == 'medium':
        # Balanced for medium systems
        params.tolerance = 10.0
        params.maxIterations = 100
        params.dampingFactor = 0.3
        params.maxDrudeDistance = 0.05
        
    else:  # large
        # Relaxed for large systems to ensure convergence
        params.tolerance = 100.0
        params.maxIterations = 500
        params.dampingFactor = 0.2
        params.maxDrudeDistance = 0.1
    
    return params

def setup_drude_system_with_optimized_params(n_waters):
    """Setup Drude system with size-appropriate SCF parameters"""
    
    # Determine system size category
    if n_waters < 10:
        size = 'small'
    elif n_waters < 50:
        size = 'medium'
    else:
        size = 'large'
    
    print(f"Setting up {n_waters} waters as '{size}' system")
    
    # Get optimized parameters
    scf_params = get_optimized_scf_params(size)
    
    # Apply parameters
    pygcmc.setDrudeSCFParameters(scf_params)
    
    print(f"SCF parameters:")
    print(f"  Tolerance: {scf_params.tolerance} kJ/mol/nm")
    print(f"  Max iterations: {scf_params.maxIterations}")
    print(f"  Damping factor: {scf_params.dampingFactor}")
    print(f"  Max Drude distance: {scf_params.maxDrudeDistance} nm")
    
    return scf_params

# Recommended usage in production code:
"""
# For any Drude calculation:
from test_optimized_scf_params import setup_drude_system_with_optimized_params

# Example:
n_waters = 27
setup_drude_system_with_optimized_params(n_waters)

# Then proceed with energy calculation
energy = pygcmc.computeSystemEnergyDrude(state)
"""

if __name__ == "__main__":
    print("=== Optimized SCF Parameters for Different System Sizes ===\n")
    
    for n in [5, 27, 100]:
        setup_drude_system_with_optimized_params(n)
        print()