#!/usr/bin/env python3
"""Fix common errors in PGP test files"""

import os
import re
from pathlib import Path

def fix_file(file_path):
    """Fix common errors in a single file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Fix 1: create_nacl_crystal missing n_cells argument
    content = re.sub(
        r'system = create_nacl_crystal\(box_size\)',
        'system = create_nacl_crystal(box_size, n_cells)',
        content
    )
    
    # Fix 2: mesh_size not defined in test_pgp_grid_convergence
    content = re.sub(
        r'setPGPParameters\(alpha, mesh_size, state\.info\.cutoff, mesh_size, 4, 1e-6\)',
        'mesh_size = [grid_size, grid_size, grid_size]\n    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)',
        content
    )
    
    # Fix 3: state should be state_pgp in pgp_cutoff_continuity.py
    content = re.sub(
        r'setPGPParameters\(alpha, mesh_size, state\.info\.cutoff, mesh_size, 4, 1e-6\)',
        'setPGPParameters(alpha, mesh_size, state_pgp.info.cutoff, mesh_size, 4, 1e-6)',
        content
    )
    
    # Fix 4: relative_error not defined in pgp_lj_minimum_energy.py
    if 'relative_error' in content and 'relative_error =' not in content:
        # Find the pattern and add the calculation
        pattern = r'(\s+# Note: Energy is stored in both residues\) / expected_lj\))\n(\s+print\(f"Relative error:)'
        replacement = r'\1\n\2'
        if re.search(pattern, content):
            content = re.sub(
                r'# Note: Energy is stored in both residues\) / expected_lj\)',
                '# Note: Energy is stored in both residues\n        lj_total = state.residues[0].energy_vdw + state.residues[1].energy_vdw\n        expected_lj = -eps  # At minimum\n        relative_error = abs((lj_total - expected_lj) / expected_lj) if expected_lj != 0 else abs(lj_total)',
                content
            )
    
    # Fix 5: real_space not defined - add after computeSystemEnergyPGP
    if 'print(f"  Real-space energy: {real_space' in content and 'real_space = ' not in content:
        # Add real_space extraction after computeSystemEnergyPGP
        content = re.sub(
            r'(computeSystemEnergyPGP\(state[^)]*\))\n',
            r'\1\n        real_space = state.ewald_energy.get("real_space", 0.0)\n',
            content
        )
    
    # Fix 6: pgp_real not defined
    if 'print(f"PGP Real-space: {pgp_real' in content and 'pgp_real = ' not in content:
        content = re.sub(
            r'(computeSystemEnergyPGP\(state[^)]*\))\n',
            r'\1\n        pgp_real = state.ewald_energy.get("real_space", 0.0)\n',
            content
        )
    
    # Fix 7: Fix initializePMEParameters missing arguments
    content = re.sub(
        r'initializePMEParameters\(state2\.info\.cutoff\)',
        'initializePMEParameters(state2.info.cutoff, state2.info.box, alpha)',
        content
    )
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Files to fix
    files_to_fix = [
        "basic_operations.py",
        "pgp_lj_minimum_energy.py",
        "pgp_real_space.py",
        "pgp_real_space_debug.py",
        "pgp_real_minimal.py",
        "pgp_cutoff_continuity.py",
        "pgp_grid_convergence.py",
        "pgp_initialization_order.py"
    ]
    
    fixed_count = 0
    for file_name in files_to_fix:
        file_path = test_dir / file_name
        if file_path.exists():
            print(f"Fixing {file_name}...")
            if fix_file(file_path):
                print(f"  ✓ Fixed")
                fixed_count += 1
            else:
                print(f"  - No changes needed")
        else:
            print(f"  ! File not found: {file_name}")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()