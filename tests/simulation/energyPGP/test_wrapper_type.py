#!/usr/bin/env python3
"""Check if Python wrapper is active."""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')

print("Testing wrapper type...")

# First import just the module
import pygcmc

print(f"\npygcmc.computeSystemEnergyPGP: {pygcmc.computeSystemEnergyPGP}")
print(f"Type: {type(pygcmc.computeSystemEnergyPGP)}")
print(f"Module: {getattr(pygcmc.computeSystemEnergyPGP, '__module__', 'No module')}")
print(f"Doc: {pygcmc.computeSystemEnergyPGP.__doc__}")

# Check if it's defined in Python
print(f"\nIs Python function? {type(pygcmc.computeSystemEnergyPGP).__name__ == 'function'}")
print(f"Has __code__? {hasattr(pygcmc.computeSystemEnergyPGP, '__code__')}")

# Check for _orig
print(f"\n_orig_computeSystemEnergyPGP exists? {hasattr(pygcmc, '_orig_computeSystemEnergyPGP')}")
if hasattr(pygcmc, '_orig_computeSystemEnergyPGP'):
    print(f"Type of _orig: {type(pygcmc._orig_computeSystemEnergyPGP)}")