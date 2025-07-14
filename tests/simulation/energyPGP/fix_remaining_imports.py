#!/usr/bin/env python3
"""
Fix remaining files that still import PGP functions directly from pygcmc.
"""

import re
import os
from pathlib import Path

def fix_imports_in_file(file_path):
    """Fix imports in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Patterns to replace
    replacements = [
        # Replace direct imports from pygcmc
        (r'from pygcmc import (.*?)setPGPParameters(.*?, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)$', 
         lambda m: m.group(0).replace('from pygcmc import', 'from . import pgp_wrapper\nfrom pgp_wrapper import'))
        
        (r'from pygcmc import (.*?)setPMEParameters(.*?)$'
         lambda m: m.group(0).replace('from pygcmc import', 'from . import pgp_wrapper\nfrom pgp_wrapper import'))
         
        (r'from pygcmc import (.*?)initializePMEParameters(.*?)$'
         lambda m: m.group(0).replace('from pygcmc import', 'from . import pgp_wrapper\nfrom pgp_wrapper import'))
         
        (r'from pygcmc import (.*?)precomputeGridPotential(.*?)$'
         lambda m: m.group(0).replace('from pygcmc import', 'from . import pgp_wrapper\nfrom pgp_wrapper import'))
         
        (r'from pygcmc import (.*?)computeSystemEnergyPGP(.*?)$'
         lambda m: m.group(0).replace('from pygcmc import'))
    ]
    
    # Apply replacements
    for pattern)
    
    # Also handle cases where multiple functions are imported on one line
    # Replace any remaining direct PGP function imports
    pgp_funcs = [
        'setPGPParameters', 'setPMEParameters', 'initializePMEParameters',
        'precomputeGridPotential', 'computeSystemEnergyPGP', 'computeMovementEnergyPGP',
        'calculateMoleculeEnergy', 'interpolateMoleculeEnergy'
    ]
    
    for func in pgp_funcs:
        # Handle: from pygcmc import func1, func2, setPGPParameters, func3
        pattern = rf'from pygcmc import ([^\\n]*){func}'
        if re.search(pattern, content):
            # First add wrapper import if not present
            if 'from . import pgp_wrapper' not in content:
                content = re.sub(r'(import pygcmc\n)', r'\1from . import pgp_wrapper\n', content)
            # Then add specific import
            if f'from pgp_wrapper import.*{func}' not in content:
                content = re.sub(r'(from . import pgp_wrapper\n)', 
                               rf'\1from pgp_wrapper import {func}\n', content)
            # Remove from pygcmc import
            content = re.sub(rf'from pygcmc import ([^\\n]*?)[\s,]*{func}[\s,]*', 
                           r'from pygcmc import \1', content)
    
    # Clean up empty imports
    content = re.sub(r'from pygcmc import\s*\n', '', content)
    content = re.sub(r'from pgp_wrapper import\s*\n', '', content)
    
    # Write back if changed
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    """Fix all files with direct PGP imports."""
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Files that still have direct imports
    files_to_fix = [
        'pgp_absolute_validation.py',
        'pgp_consistency_tests.py',
        'pgp_cutoff_continuity.py',
        'pgp_cutoff_issue.py',
        'pgp_erfc_debug.py',
        'pgp_grid_convergence.py',
        'pgp_initialization_order.py',
        'pgp_lj_boundary_tests.py',
        'pgp_lj_diagnosis.py',
        'pgp_lj_extreme_distances.py',
        'pgp_lj_minimum_energy.py',
        'pgp_lj_neutral_systems.py',
        'pgp_mixed_interactions.py',
        'pgp_pme_complete_validation.py',
        'pgp_pme_cutoff.py',
        'pgp_pme_debug_electrostatic.py',
        'pgp_pme_electrostatic_validation.py',
        'pgp_real_fixed.py',
        'pgp_real_minimal.py',
        'pgp_real_space.py',
        'pgp_real_space_debug.py',
        'vspme_combined.py',
        'vspme_delta_energies.py',
        'vspme_direct_lj.py',
        'vspme_lj_energy.py',
        'vspme_two_atom.py',
    ]
    
    fixed_count = 0
    for file_name in files_to_fix:
        file_path = test_dir / file_name
        if file_path.exists():
            if fix_imports_in_file(file_path):
                print(f"Fixed: {file_name}")
                fixed_count += 1
        else:
            print(f"File not found: {file_name}")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()
