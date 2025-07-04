# tests/simulation/energyPGP/diagnose_erfc_issue.py
"""
Diagnose the erfc issue by process of elimination.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def diagnose_erfc_issue():
    """Diagnose erfc issue."""
    print("\n=== Diagnosing erfc Issue ===")
    
    # Test at different distances to see pattern
    distances = [0.5, 0.8, 1.0, 1.2]
    kC = 138.935456
    alpha = 2.5
    
    print("\nExpected vs Actual erfcApprox values:")
    print("Distance | Expected erfc | Inferred erfc | Ratio")
    print("-" * 50)
    
    for r in distances:
        # Create system
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.cutoff = 1.5  # Larger than all test distances
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        atoms = []
        atom1 = MCAtom()
        atom1.x = 5.0
        atom1.y = 5.0
        atom1.z = 5.0
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 5.0 + r
        atom2.y = 5.0
        atom2.z = 5.0
        atom2.charge = -1.0
        atom2.type = 0
        atoms.append(atom2)
        
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
        mesh_size = [32, 32, 32]
        setPMEParameters(alpha, mesh_size, 4, 1e-6)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        computeSystemEnergyPME(state)
        
        # Get residue energy
        res_energy = state.residues[0].energy_elec
        
        # Infer erfcApprox value
        # res_energy = 0.5 * (-1) * erfcApprox(r) / r * COULOMB²
        # So: erfcApprox(r) = -2 * res_energy * r / COULOMB²
        inferred_erfc = -2 * res_energy * r / (kC * kC)
        expected_erfc = math.erfc(alpha * r)
        
        if expected_erfc > 0:
            ratio = inferred_erfc / expected_erfc
        else:
            ratio = float('inf')
        
        print(f"{r:8.1f} | {expected_erfc:13.6f} | {inferred_erfc:13.6f} | {ratio:8.1f}")
    
    print("\n\nAnalysis:")
    print("If erfcApprox always returns the same value regardless of r,")
    print("then the inferred values should scale with r.")
    print("If the inferred values are constant, then erfcApprox(r) = constant/r")
    
    # Check if it's returning 1/r
    print("\nChecking if erfcApprox(r) = 1/r:")
    for r in distances:
        expected = 1.0 / r
        print(f"  1/{r} = {expected:.6f}")


if __name__ == "__main__":
    diagnose_erfc_issue()