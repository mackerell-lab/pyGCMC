#!/usr/bin/env python3
"""Fix all import issues for test_energy_PGP.py tests."""

import os
import re
from pathlib import Path

def fix_imports_in_file(file_path):
    """Fix imports in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Fix imports from pygcmc to pgp_wrapper
    pgp_functions = [
        'setPGPParameters', 'initializePMEParameters', 'setPMEParameters',
        'precomputeGridPotential', 'computeSystemEnergyPGP', 'computeMovementEnergyPGP',
        'calculateMoleculeEnergy', 'computeSystemEnergyPME', 'computeMovementEnergyPME',
        'computeSystemEnergyPMEComplete', 'computeMovementEnergyPMEFixed',
        'computeSystemEnergyPMEFixed', 'computeSystemEnergyEwald',
        'computeSystemVdwEnergyCutoff', 'interpolateMoleculeEnergy'
    ]
    
    # Replace direct imports from pygcmc
    for func in pgp_functions:
        # Pattern 1: from pygcmc import ..., func, ...
        pattern1 = rf'(from\s+pygcmc\s+import\s+[^;\n]*?)(\b{func}\b)'
        content = re.sub(pattern1, lambda m: m.group(0).replace(m.group(0), ''), content)
        
        # Pattern 2: from pygcmc import func
        pattern2 = rf'^from\s+pygcmc\s+import\s+{func}$'
        content = re.sub(pattern2, f'from .pgp_wrapper import {func}', content, flags=re.MULTILINE)
    
    # Add pgp_wrapper import if needed
    if 'from .pgp_wrapper import' in content and 'from . import pgp_wrapper' not in content:
        # Add after other imports
        content = re.sub(r'(import pygcmc\n)', r'\1from . import pgp_wrapper\n', content)
    
    # Clean up multiple imports on same line
    content = re.sub(r'from pygcmc import\s*,', 'from pygcmc import', content)
    content = re.sub(r',\s*,', ',', content)
    content = re.sub(r'from pygcmc import\s*\n', '', content)
    
    # Add missing imports to pgp_wrapper
    for func in pgp_functions:
        if func in content and f'from .pgp_wrapper import.*{func}' not in content:
            # Check if we already have some pgp_wrapper imports
            if 'from .pgp_wrapper import' in content:
                # Add to existing import
                content = re.sub(r'(from \.pgp_wrapper import [^\n]+)', rf'\1, {func}', content, count=1)
            else:
                # Add new import after pygcmc imports
                content = re.sub(r'(from pygcmc import [^\n]+\n)', rf'\1from .pgp_wrapper import {func}\n', content, count=1)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Files that are imported by test_energy_PGP.py
    target_files = [
        'vspme_two_atom.py',
        'method_comparison.py',
        'complex_systems.py',
        'planar_systems.py',
        'asymmetric_water_complete.py',
        'asymmetric_nacl.py',
        'pgp_real_space.py',
        'pgp_real_space_debug.py',
        'pgp_real_minimal.py',
        'pgp_real_fixed.py',
        'pgp_initialization_order.py',
        'pgp_pme_cutoff.py',
        'pgp_cutoff_issue.py',
        'pgp_erfc_debug.py',
        'pgp_pme_debug_electrostatic.py',
        'pgp_pme_electrostatic_validation.py',
        'debug_movement_residues.py',
        'debug_vdw_movement.py',
        'pgp_complete.py'
    ]
    
    fixed_count = 0
    
    for file_name in target_files:
        file_path = test_dir / file_name
        if file_path.exists():
            if fix_imports_in_file(file_path):
                print(f"Fixed imports in: {file_name}")
                fixed_count += 1
        else:
            print(f"Warning: File not found: {file_name}")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()
