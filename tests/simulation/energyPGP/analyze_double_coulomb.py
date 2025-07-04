# tests/simulation/energyPGP/analyze_double_coulomb.py
"""
Analyze where COULOMB is being applied twice.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def analyze_double_coulomb():
    """Analyze the double COULOMB application."""
    print("\n=== Analyzing Double COULOMB Application ===")
    
    COULOMB = 138.935456
    
    # The observed value is -19303.06 for each residue
    # This is exactly COULOMB² / 2
    print(f"COULOMB = {COULOMB}")
    print(f"COULOMB² = {COULOMB**2:.2f}")
    print(f"COULOMB² / 2 = {COULOMB**2 / 2:.2f}")
    print(f"Observed residue energy: -19303.06")
    
    # Working backwards:
    # If final = -19303.06 = -COULOMB²/2
    # And we know COULOMB is applied in PMEComposite.cpp line 78
    # Then before that multiplication: -19303.06 / COULOMB = -138.935456
    
    print("\n\nWorking backwards:")
    before_composite = -19303.06 / COULOMB
    print(f"Before PMEComposite COULOMB: {before_composite:.6f}")
    
    # This suggests pair_energy in PMEReal.cpp is already -138.935456
    # Which means pair_energy = q1 * q2 * erfc * COULOMB / r
    # With q1*q2 = -1, erfc = 1.0, r = 1.0:
    # pair_energy = -1 * 1.0 * COULOMB / 1.0 = -COULOMB
    
    print("\nThis suggests PMEReal.cpp is calculating:")
    print("  pair_energy = q1 * q2 * erfcApprox(r) * COULOMB / r")
    print("  Instead of: q1 * q2 * erfcApprox(r) / r")
    
    print("\n\nOr possibly:")
    print("1. erfcApprox returns 1.0")
    print("2. pair_energy = -1.0 * 1.0 / 1.0 = -1.0")
    print("3. Each residue gets full energy: -1.0 (not half)")
    print("4. PMEComposite applies COULOMB: -138.935456")
    print("5. Some other place applies COULOMB again: -19303.06")


if __name__ == "__main__":
    analyze_double_coulomb()