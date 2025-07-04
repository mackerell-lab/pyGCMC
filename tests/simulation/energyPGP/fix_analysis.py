# tests/simulation/energyPGP/fix_analysis.py
"""
Correct analysis of the issue.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import math


def fix_analysis():
    """Correct analysis."""
    print("\n=== Correct Analysis ===")
    
    COULOMB = 138.935456
    observed = -19303.060547
    
    print(f"Observed residue energy: {observed:.6f}")
    print(f"COULOMB = {COULOMB:.6f}")
    print(f"COULOMB² = {COULOMB**2:.6f}")
    
    # Check what multiple of COULOMB this is
    ratio = observed / (-COULOMB)
    print(f"\nObserved / (-COULOMB) = {ratio:.6f}")
    print(f"This is very close to COULOMB = {COULOMB:.6f}")
    
    # So observed ≈ -COULOMB²
    print(f"\nSo the observed value is -COULOMB²")
    
    # This means:
    # 1. pair_energy starts as -1.0 (because erfcApprox returns 1.0)
    # 2. Each residue gets FULL pair energy: -1.0
    # 3. PMEComposite.cpp line 78 applies COULOMB: -138.935456
    # 4. Somewhere else COULOMB is applied again: -19303.06
    
    print("\n\nThe flow must be:")
    print("1. erfcApprox returns 1.0 instead of 0.000407")
    print("2. pair_energy = -1.0 * 1.0 / 1.0 = -1.0")
    print("3. Each residue gets FULL energy = -1.0 (GCMC pattern)")
    print("4. PMEComposite multiplies by COULOMB: -138.935456")
    print("5. Something else multiplies by COULOMB again: -19303.06")
    
    print("\n\nTo fix:")
    print("1. Fix erfcApprox to return correct values")
    print("2. Check if residues should get half or full energy")
    print("3. Find where COULOMB is applied twice")


if __name__ == "__main__":
    fix_analysis()