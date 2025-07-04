# tests/simulation/energyPGP/test_erfc_approx.py
"""
Test what erfcApprox is actually returning.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import setPMEParameters, initializePMEParameters
import math


def test_erfc_approx():
    """Test erfcApprox function."""
    print("\n=== Testing erfcApprox ===")
    
    # Initialize PME parameters
    alpha = 2.5
    mesh_size = [32, 32, 32]
    cutoff = 1.2
    box = [10.0, 10.0, 10.0]
    
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(cutoff, box, alpha)
    
    # Test at r = 1.0
    r = 1.0
    expected_erfc = math.erfc(alpha * r)
    
    print(f"Testing at r = {r} nm:")
    print(f"  Alpha = {alpha}")
    print(f"  Expected erfc(αr) = erfc({alpha * r}) = {expected_erfc:.6f}")
    
    # Unfortunately we can't directly call erfcApprox from Python
    # But we can infer what it's returning from the energy calculation
    
    # From our previous tests:
    # - Residue energy (no COULOMB) = -1.0 for qi=1, qj=-1, r=1
    # - Formula: energy = qi * qj * erfcApprox(r) / r
    # - So: -1.0 = (-1) * erfcApprox(1) / 1
    # - Therefore: erfcApprox(1) = 1.0
    
    inferred_erfc = 1.0  # Based on the residue energy we're seeing
    
    print(f"\nInferred erfcApprox(1.0) = {inferred_erfc}")
    print(f"Expected erfcApprox(1.0) = {expected_erfc:.6f}")
    print(f"Ratio: {inferred_erfc / expected_erfc:.1f}")
    
    # This ratio is 2457.3, which is 1/erfc(2.5)
    # This strongly suggests erfcApprox is returning 1.0 instead of erfc(αr)
    
    print(f"\nPossible explanations:")
    print(f"1. erfcApprox is returning 1.0 (constant)")
    print(f"2. erfcApprox is returning 1/r (would give 1.0 at r=1)")
    print(f"3. erfcApprox is returning r/erfc(αr) (would give 2457.3 at r=1)")
    
    # Let's check what happens if erfcApprox returns different values
    print(f"\nIf erfcApprox returns:")
    print(f"  1.0: energy = -1 * 1.0 / 1 = -1.0 ✓ (matches what we see)")
    print(f"  1/r: energy = -1 * 1.0 / 1 = -1.0 ✓ (also matches)")
    print(f"  erfc(αr): energy = -1 * {expected_erfc:.6f} / 1 = {-expected_erfc:.6f} ✗")
    
    print(f"\n⚠️ CONCLUSION: erfcApprox is NOT returning erfc(αr)!")
    print(f"   It appears to be returning a constant 1.0 or 1/r")


if __name__ == "__main__":
    test_erfc_approx()