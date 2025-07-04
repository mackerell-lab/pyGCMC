# tests/simulation/energyPGP/debug_table_lookup.py
"""
Debug the table lookup issue.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def debug_table_lookup():
    """Debug why erfcApprox returns 1.0."""
    print("\n=== Debugging Table Lookup ===")
    
    # The issue might be:
    # 1. erfcDXInv is 0, making x = r * 0 = 0, so index = 0
    # 2. erfcTable[0] = erfc(0) = 1.0
    
    print("\nExpected erfc values:")
    alpha = 2.5
    for r in [0.0, 0.5, 1.0, 1.5]:
        print(f"  erfc({alpha}*{r}) = {math.erfc(alpha*r):.6f}")
    
    print("\nIf erfcApprox always returns 1.0:")
    print("  - The table lookup always returns index 0")
    print("  - erfcTable[0] = erfc(0) = 1.0")
    print("  - This happens if erfcDXInv = 0 or very small")
    
    # erfcDXInv = (NUM_TABLE_POINTS - 1) / tableRange
    # tableRange = cutoff
    # So erfcDXInv = 19999 / 1.2 = 16665.83
    
    print("\nExpected erfcDXInv calculation:")
    cutoff = 1.2
    NUM_TABLE_POINTS = 20000
    expected_erfcDXInv = (NUM_TABLE_POINTS - 1) / cutoff
    print(f"  erfcDXInv = ({NUM_TABLE_POINTS}-1) / {cutoff} = {expected_erfcDXInv:.2f}")
    
    print("\nFor r=1.0:")
    print(f"  x = r * erfcDXInv = 1.0 * {expected_erfcDXInv:.2f} = {expected_erfcDXInv:.2f}")
    print(f"  index = int(x) = {int(expected_erfcDXInv)}")
    print(f"  This should access erfcTable[{int(expected_erfcDXInv)}], not erfcTable[0]")
    
    print("\nPossible causes:")
    print("1. erfcDXInv is not initialized (= 0)")
    print("2. erfcTable is not properly populated")
    print("3. The table lookup calculation has a bug")
    print("4. We're calling a different erfcApprox function")
    
    # Let me check if it's an initialization order issue
    print("\n\nTesting initialization order:")
    
    # Create state
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    for i in range(2):
        atom = MCAtom()
        atom.x = 5.0 + i * 0.5  # 0.5 nm apart
        atom.y = 5.0
        atom.z = 5.0
        atom.charge = 1.0 if i == 0 else -1.0
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Initialize and compute
    setPMEParameters(2.5, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    
    # At 0.5 nm distance
    res_energy = state.residues[0].energy_elec
    # res_energy = -erfcApprox(0.5) / 0.5 * COULOMB²
    inferred_erfc = -res_energy * 0.5 / (138.935456**2)
    
    print(f"\nAt r=0.5 nm:")
    print(f"  Residue energy: {res_energy:.2f}")
    print(f"  Inferred erfcApprox(0.5): {inferred_erfc:.6f}")
    print(f"  Expected erfc(2.5*0.5): {math.erfc(2.5*0.5):.6f}")
    
    if abs(inferred_erfc - 1.0) < 0.01:
        print("\n⚠️ erfcApprox still returns 1.0 even for r=0.5!")
        print("This confirms the table lookup is broken.")


if __name__ == "__main__":
    debug_table_lookup()