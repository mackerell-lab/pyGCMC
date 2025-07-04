# tests/simulation/energyPGP/debug_erfc_table.py
"""
Debug the erfc table issue.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def debug_erfc_table():
    """Debug erfc table initialization."""
    print("\n=== Debugging erfc Table ===")
    
    # Create minimal system with varying distances
    for test_r in [0.5, 0.75, 1.0, 1.1]:
        print(f"\n--- Testing at r = {test_r} nm ---")
        
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.cutoff = 1.2
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two particles at distance test_r
        atoms = []
        
        atom1 = MCAtom()
        atom1.x = 5.0
        atom1.y = 5.0
        atom1.z = 5.0
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 5.0 + test_r
        atom2.y = 5.0
        atom2.z = 5.0
        atom2.charge = -1.0
        atom2.type = 0
        atoms.append(atom2)
        
        state.atoms = atoms
        state.activeAtomCount = 2
        
        # Two residues
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
        
        # Initialize PME
        alpha = 2.5
        mesh_size = [32, 32, 32]
        setPMEParameters(alpha, mesh_size, 4, 1e-6)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        
        # Compute
        computeSystemEnergyPME(state)
        
        # Get residue energy (before COULOMB²)
        res_energy = state.residues[0].energy_elec
        
        # Infer what erfcApprox returned
        # energy = qi * qj * erfcApprox(r) / r
        # res_energy = energy/2 * COULOMB²
        # So: erfcApprox(r) = res_energy * 2 / (qi * qj * COULOMB²) * r
        
        kC = 138.935456
        inferred_erfc = res_energy * 2 / (-1.0 * kC * kC) * test_r
        expected_erfc = math.erfc(alpha * test_r)
        
        print(f"  Expected erfc({alpha}*{test_r}) = {expected_erfc:.6f}")
        print(f"  Inferred erfcApprox({test_r}) = {inferred_erfc:.6f}")
        print(f"  Ratio: {inferred_erfc / expected_erfc:.1f}")
        print(f"  Residue energy: {res_energy:.2f}")
        
        # Check if it's returning a constant
        if abs(inferred_erfc - 1.0) < 0.01:
            print(f"  ⚠️ erfcApprox appears to return constant 1.0")
        elif abs(inferred_erfc / test_r - 1.0) < 0.01:
            print(f"  ⚠️ erfcApprox appears to return 1/r")
    
    print("\n\nCONCLUSION:")
    print("The erfcApprox function is returning 1.0 for all distances!")
    print("This suggests the erfc table is not being properly initialized or accessed.")


if __name__ == "__main__":
    debug_erfc_table()