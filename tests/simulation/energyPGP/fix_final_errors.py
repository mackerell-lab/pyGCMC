#!/usr/bin/env python3
"""Fix final errors in PGP test files"""

import os
import re
from pathlib import Path

def fix_real_space_errors():
    """Fix real_space variable not defined errors"""
    
    # Fix pgp_real_space.py
    file_path = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP/pgp_real_space.py")
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Add real_space extraction after computeSystemEnergyPGP
    content = re.sub(
        r'(computeSystemEnergyPGP\(state\))\n(\s+)reciprocal = state\.ewald_energy\.get\(\'reciprocal\'\)',
        r'\1\n\2real_space = state.ewald_energy.get("real_space", 0.0)\n\2reciprocal = state.ewald_energy.get("reciprocal")',
        content
    )
    
    # Fix line 226 specifically
    content = re.sub(
        r'(precomputeGridPotential\(state\))\n(\s+)grid_energy = state\.ewald_energy\.get\(\'reciprocal\'\)',
        r'\1\n\2computeSystemEnergyPGP(state)\n\2real_space = state.ewald_energy.get("real_space", 0.0)\n\2grid_energy = state.ewald_energy.get("reciprocal")',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)
    print("Fixed pgp_real_space.py")
    
    # Fix pgp_real_space_debug.py (water molecule test)
    file_path = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP/pgp_real_space_debug.py")
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Add real_space extraction
    content = re.sub(
        r'(precomputeGridPotential\(state\))\n(\s+)self_energy = state\.ewald_energy\.get\(\'self\'\)',
        r'\1\n\2computeSystemEnergyPGP(state)\n\2real_space = state.ewald_energy.get("real_space", 0.0)\n\2self_energy = state.ewald_energy.get("self")',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)
    print("Fixed pgp_real_space_debug.py")

def fix_grid_size_errors():
    """Fix grid_size not defined errors"""
    
    file_path = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP/pgp_initialization_order.py")
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Find test_pgp_initialization_methods and add grid_size
    content = re.sub(
        r'(# Calculate PGP energy)\n(\s+computeSystemEnergyPGP\(state\))',
        r'# Calculate PGP energy\n        grid_size = 32  # Define grid_size\n\2',
        content
    )
    
    # Find test_pgp_with_correct_initialization and add grid_size
    content = re.sub(
        r'(# Grid size tests)\n(\s+grids = \[16, 32, 64\])',
        r'\1\n        grid_size = 32  # Default grid size\n\2',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)
    print("Fixed pgp_initialization_order.py")

def main():
    print("Fixing final errors in PGP test files...")
    fix_real_space_errors()
    fix_grid_size_errors()
    print("\nDone!")

if __name__ == "__main__":
    main()