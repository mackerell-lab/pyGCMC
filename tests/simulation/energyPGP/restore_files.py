#!/usr/bin/env python3
"""
Restore files from git to fix broken syntax.
"""

import subprocess
import sys
from pathlib import Path

def restore_file(file_path):
    """Restore a file from git."""
    try:
        subprocess.run(['git', 'checkout', file_path], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False

def main():
    # Files with syntax errors that need to be restored
    broken_files = [
        'asymmetric_nacl_helpers.py',
        'asymmetric_water_movement.py',
        'debug_movement_residues.py',
        'debug_vdw_movement.py',
        'find_2_source.py',
        'helpers.py',
        'pgp_absolute_validation.py',
        'pgp_complete.py',
        'pgp_complete_validation.py',
        'pgp_consistency_tests.py',
        'pgp_cutoff_continuity.py',
        'pgp_cutoff_issue.py',
        'pgp_erfc_debug.py',
        'pgp_grid_convergence.py',
        'pgp_grid_convergence_fixed.py',
        'pgp_initialization_order.py',
        'pgp_lj_boundary_tests.py',
        'pgp_lj_diagnosis.py',
        'pgp_lj_extreme_distances.py',
        'pgp_lj_helpers.py',
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
        'pgp_simple_test.py',
        'run_isolated.py',
        'test_pgp_context.py',
        'test_pgp_context_migration.py',
        'test_pgp_grid_convergence_independent.py',
        'test_pgp_independent.py',
        'verify_erfc_init.py',
        'vspme_direct_helpers.py',
        'vspme_energy_diagnostics.py',
        'vspme_lj_helpers.py',
        'vspme_two_atom.py',
    ]
    
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    restored = 0
    for file_name in broken_files:
        file_path = test_dir / file_name
        if restore_file(file_path):
            print(f"Restored: {file_name}")
            restored += 1
        else:
            print(f"Failed to restore: {file_name}")
    
    print(f"\nRestored {restored} files")
    
    # Now apply the wrapper import fix again
    print("\nApplying wrapper import fixes...")
    subprocess.run([sys.executable, 'fix_pgp_to_context.py'], cwd=test_dir)

if __name__ == "__main__":
    main()
