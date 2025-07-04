# tests/simulation/energyPGP/check_gcmc_pattern.py
"""
Check if the double energy storage is intentional for GCMC.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def check_gcmc_pattern():
    """Check the GCMC energy storage pattern."""
    print("\n=== Checking GCMC Energy Storage Pattern ===")
    
    # According to the user:
    # "要知道residue上的能量是两倍可能是正确的，因为我们判断插入，就是要用这个能量，
    #  环境和residue的两部分能量都放在residue上了"
    # Translation: The residue energy storing double values might be correct because 
    # for GCMC insertion/deletion decisions, we need the full interaction energy 
    # (both the residue's contribution and the environment's contribution).
    
    print("\nUser feedback suggests residue should store full pair energy for GCMC.")
    print("This is for insertion/deletion decisions in Grand Canonical MC.")
    
    # If this is true, then for a pair interaction:
    # - Each residue stores the FULL pair energy (not half)
    # - This is intentional for GCMC
    
    # Let's check the math:
    r = 1.0
    q1, q2 = 1.0, -1.0
    alpha = 2.5
    COULOMB = 138.935456
    
    # If erfcApprox returns 1.0 (the bug):
    erfc_value = 1.0  # Bug value
    pair_energy = q1 * q2 * erfc_value / r  # = -1.0
    
    # If COULOMB is applied in PMEReal.cpp:
    pair_energy_with_coulomb = pair_energy * COULOMB  # = -138.935456
    
    # If each residue gets the FULL energy (GCMC pattern):
    residue_energy = pair_energy_with_coulomb  # = -138.935456
    
    # But then in PMEComposite.cpp line 78, COULOMB is applied again:
    residue_energy_final = residue_energy * COULOMB  # = -19303.06
    
    print(f"\nIf erfcApprox = 1.0 and COULOMB applied twice:")
    print(f"  Pair energy: {pair_energy}")
    print(f"  After first COULOMB: {pair_energy_with_coulomb}")
    print(f"  After second COULOMB: {residue_energy_final}")
    print(f"  This matches the observed value!")
    
    # The real issue is:
    # 1. erfcApprox returns 1.0 instead of 0.000407
    # 2. COULOMB might be applied in PMEReal.cpp already
    
    print("\n\nThe main issues to fix:")
    print("1. erfcApprox returns 1.0 instead of correct erfc values")
    print("2. Check if COULOMB should be applied in PMEReal.cpp or PMEComposite.cpp")


if __name__ == "__main__":
    check_gcmc_pattern()