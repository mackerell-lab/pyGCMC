#!/usr/bin/env python3
"""Fix remaining errors in PGP test files"""

import os
import re
from pathlib import Path

def fix_pgp_grid_convergence(file_path):
    """Fix specific issues in pgp_grid_convergence.py"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Fix state_pgp -> state
    content = re.sub(r'state_pgp\.info\.cutoff', 'state.info.cutoff', content)
    
    # Remove duplicate mesh_size lines
    content = re.sub(
        r'mesh_size = \[grid_size, grid_size, grid_size\]\n\s*mesh_size = \[grid_size, grid_size, grid_size\]',
        'mesh_size = [grid_size, grid_size, grid_size]',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)

def fix_pgp_cutoff_continuity(file_path):
    """Fix specific issues in pgp_cutoff_continuity.py"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Fix undefined grid_size
    content = re.sub(
        r'mesh_size = \[grid_size, grid_size, grid_size\]',
        '# mesh_size already defined above',
        content
    )
    
    # Fix state_pgp in setPGPParameters where it should be state
    content = re.sub(
        r'setPGPParameters\(alpha, mesh_size, state_pgp\.info\.cutoff, mesh_size, 4, 1e-6\)',
        'setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)

def fix_pgp_real_space_debug(file_path):
    """Fix pgp_real_space_debug.py"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Add pgp_real extraction
    if 'print(f"PGP Real-space: {pgp_real' in content and 'pgp_real = ' not in content:
        content = re.sub(
            r'(precomputeGridPotential\(state\))\n(\s+pgp_recip = state\.ewald_energy\.get\(\'reciprocal\'\))',
            r'\1\n        computeSystemEnergyPGP(state)\n        pgp_real = state.ewald_energy.get("real_space", 0.0)\n\2',
            content
        )
    
    with open(file_path, 'w') as f:
        f.write(content)

def fix_pgp_real_minimal(file_path):
    """Fix pgp_real_minimal.py"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Add real_space extraction
    if 'print(f"Real-space energy: {real_space' in content and 'real_space = ' not in content:
        content = re.sub(
            r'(precomputeGridPotential\(state\))\n(\s+reciprocal = state\.ewald_energy\.get\(\'reciprocal\'\))',
            r'\1\n        computeSystemEnergyPGP(state)\n        real_space = state.ewald_energy.get("real_space", 0.0)\n\2',
            content
        )
    
    with open(file_path, 'w') as f:
        f.write(content)

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Fix specific files
    print("Fixing pgp_grid_convergence.py...")
    fix_pgp_grid_convergence(test_dir / "pgp_grid_convergence.py")
    
    print("Fixing pgp_cutoff_continuity.py...")
    fix_pgp_cutoff_continuity(test_dir / "pgp_cutoff_continuity.py")
    
    print("Fixing pgp_real_space_debug.py...")
    fix_pgp_real_space_debug(test_dir / "pgp_real_space_debug.py")
    
    print("Fixing pgp_real_minimal.py...")
    fix_pgp_real_minimal(test_dir / "pgp_real_minimal.py")
    
    print("\nDone!")

if __name__ == "__main__":
    main()