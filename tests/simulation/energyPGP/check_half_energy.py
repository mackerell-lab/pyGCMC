# tests/simulation/energyPGP/check_half_energy.py
"""
Check if residues are getting half or full energy.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import math


def check_half_energy():
    """Analyze if residues get half or full energy."""
    print("\n=== Checking Half vs Full Energy ===")
    
    COULOMB = 138.935456
    observed = -19303.060547
    
    # If the issue is:
    # 1. erfcApprox = 1.0
    # 2. pair_energy = -1.0
    # 3. Each residue gets FULL energy (not half)
    # 4. COULOMB applied once: -138.935456
    # 5. COULOMB applied again: -19303.06
    
    print("Scenario 1: Residues get half energy (normal)")
    half_energy = -1.0 / 2.0
    after_coulomb1 = half_energy * COULOMB
    after_coulomb2 = after_coulomb1 * COULOMB
    print(f"  pair_energy = -1.0")
    print(f"  per residue = {half_energy}")
    print(f"  after COULOMB = {after_coulomb1}")
    print(f"  after COULOMB² = {after_coulomb2}")
    
    print("\nScenario 2: Residues get full energy (GCMC)")
    full_energy = -1.0
    after_coulomb1 = full_energy * COULOMB
    after_coulomb2 = after_coulomb1 * COULOMB
    print(f"  pair_energy = -1.0")
    print(f"  per residue = {full_energy}")
    print(f"  after COULOMB = {after_coulomb1}")
    print(f"  after COULOMB² = {after_coulomb2}")
    print(f"  This matches observed: {observed}")
    
    print("\n\nConclusion:")
    print("1. Residues ARE getting full pair energy (GCMC pattern)")
    print("2. COULOMB is being applied twice")
    print("3. erfcApprox returns 1.0 instead of 0.000407")


if __name__ == "__main__":
    check_half_energy()