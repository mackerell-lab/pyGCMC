#!/usr/bin/env python
"""
独立运行有问题的测试，避免内存崩溃
"""

import subprocess
import sys
import os

def run_test_in_subprocess(test_path):
    """在子进程中运行测试"""
    env = os.environ.copy()
    env['PYTHONPATH'] = env.get('PYTHONPATH', '') + ':./modules/bindings'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MALLOC_CHECK_'] = '0'
    
    cmd = [sys.executable, '-m', 'pytest', test_path, '-v', '-s', '--tb=short']
    
    try:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=60)
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"Test {test_path} timed out")
        return False
    except Exception as e:
        print(f"Error running test: {e}")
        return False

if __name__ == "__main__":
    os.chdir('/home/zhaomt/gcmc/test107/pygcmc_dev/build')
    
    tests = [
        '../tests/simulation/energyPGP/pgp_pme_debug_electrostatic.py::test_pgp_self_consistency',
        '../tests/simulation/energyPGP/debug_movement_residues.py::test_movement_residues',
        '../tests/simulation/energyPGP/debug_vdw_movement.py::test_vdw_movement_debug'
    ]
    
    for test in tests:
        print(f"\n{'='*60}")
        print(f"Running: {test}")
        print('='*60)
        
        success = run_test_in_subprocess(test)
        
        if not success:
            print(f"\n⚠️  Test {test} failed or crashed")
        else:
            print(f"\n✓ Test {test} passed")
    
    print("\n" + "="*60)
    print("Note: These tests have memory issues and should be run in isolation")