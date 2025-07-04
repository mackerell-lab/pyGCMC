# tests/simulation/energyPGP/test_wrapper.py
"""
Test if wrapper is working.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc

print("\n=== Testing Wrapper ===")
print(f"computeSystemEnergyPGP function: {pygcmc.computeSystemEnergyPGP}")
print(f"Is it our wrapper? {hasattr(pygcmc.computeSystemEnergyPGP, '__name__')}")
if hasattr(pygcmc.computeSystemEnergyPGP, '__name__'):
    print(f"Function name: {pygcmc.computeSystemEnergyPGP.__name__}")
    
print(f"\n_orig_computeSystemEnergyPGP: {pygcmc._orig_computeSystemEnergyPGP if hasattr(pygcmc, '_orig_computeSystemEnergyPGP') else 'Not found'}")