#!/usr/bin/env python
"""
Simple runner to execute PGP Complete tests in isolation.

Run from build directory:
    cd /home/zhaomt/gcmc/test107/pygcmc_dev/build
    python ../tests/simulation/energyPGP/run_isolated.py
"""

import subprocess
import sys
import os
from pathlib import Path

def run_test_isolated(test_name, workdir=None):
    """Run a single test in an isolated process"""
    workdir = Path(workdir or Path.cwd())
    cmd = [
        sys.executable, "-m", "pytest",
        f"../tests/simulation/test_energy_PGP.py::{test_name}",
        "-v", "-s", "--tb=short"
    ]
    
    env = os.environ.copy()
    env['PYTHONPATH'] = f"{os.getcwd()}/modules/bindings:{env.get('PYTHONPATH', '')}"
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    
    print(f"\n{'='*70}")
    print(f"Running {test_name} in isolated process...")
    print('='*70)
    
    result = subprocess.run(cmd, env=env, cwd=str(workdir))
    return result.returncode == 0

def main():
    """Run all PGP Complete tests"""
    tests = [
        "test_pgp_complete_vs_pme_complete",
        "test_pgp_complete_movement_energy",
        "test_pgp_complete_pure_lj", 
        "test_pgp_complete_multi_atom_residue",
        "test_pgp_complete_extreme_distances",
        "test_pgp_complete_direct_movement_test"
    ]
    
    print("PGP Complete Tests - Running in Isolated Processes")
    print("="*70)
    print("This avoids the global state bug in the C++ PGP implementation")
    print("="*70)
    
    passed = 0
    failed = 0
    
    workdir = Path.cwd()
    for test in tests:
        if run_test_isolated(test, workdir=workdir):
            passed += 1
            print(f"✓ {test}")
        else:
            failed += 1
            print(f"✗ {test}")
    
    print("\n" + "="*70)
    print(f"Summary: {passed} passed, {failed} failed")
    
    if failed > 0:
        print("\nNote: Failures here indicate real test failures, not crashes.")
        print("The isolation prevents the memory corruption bug.")
    
    return 0 if failed == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
