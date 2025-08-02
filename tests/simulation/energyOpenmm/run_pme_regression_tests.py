#!/usr/bin/env python
"""
Independent PME regression test script

Since PME uses global state, these tests cannot run in parallel test environments.
Use this script to separately verify the correctness of PME Total energy calculation.
"""

import sys
import os

# Add path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Ensure pygcmc can be imported
if 'PYTHONPATH' not in os.environ:
    os.environ['PYTHONPATH'] = ''
if 'modules/bindings' not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = os.path.abspath('modules/bindings') + ':' + os.environ['PYTHONPATH']
    
# Reset sys.path
bindings_path = os.path.abspath("modules/bindings")
if bindings_path not in sys.path:
    sys.path.insert(0, bindings_path)

# Import test functions (remove @pytest.mark.skip decorator effects)
from pme_residue_offset_pattern import test_residue_offset_pattern
from analyze_pme_total_bug import test_pme_total_configurations

def run_tests():
    """Run PME regression tests"""
    print("=" * 80)
    print("Running PME Total energy calculation regression tests")
    print("=" * 80)
    
    # Test 1: Residue offset pattern test
    print("\nTest 1: test_residue_offset_pattern")
    print("-" * 40)
    try:
        # Get original function (skip decorators)
        if hasattr(test_residue_offset_pattern, '__wrapped__'):
            test_func = test_residue_offset_pattern.__wrapped__
        else:
            test_func = test_residue_offset_pattern
        
        test_func()
        print("✅ test_residue_offset_pattern passed")
    except Exception as e:
        print(f"❌ test_residue_offset_pattern failed: {e}")
        return False
    
    # Test 2: PME Total configuration test
    print("\nTest 2: test_pme_total_configurations")
    print("-" * 40)
    try:
        # Get original function (skip decorators)
        if hasattr(test_pme_total_configurations, '__wrapped__'):
            test_func = test_pme_total_configurations.__wrapped__
        else:
            test_func = test_pme_total_configurations
            
        test_func()
        print("✅ test_pme_total_configurations passed")
    except Exception as e:
        print(f"❌ test_pme_total_configurations failed: {e}")
        return False
    
    print("\n" + "=" * 80)
    print("✅ All PME regression tests passed!")
    print("PME Total energy calculation formula has been correctly fixed: total = real_space + reciprocal + self")
    print("=" * 80)
    return True

if __name__ == "__main__":
    # Run tests
    success = run_tests()
    sys.exit(0 if success else 1)