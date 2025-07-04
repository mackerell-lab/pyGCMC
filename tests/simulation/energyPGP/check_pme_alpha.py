# tests/simulation/energyPGP/check_pme_alpha.py
"""
Check if PME alpha is being set correctly.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def check_pme_alpha():
    """Check PME alpha parameter."""
    print("\n=== Checking PME Alpha ===")
    
    # Test with different alpha values
    for alpha in [1.0, 2.0, 2.5, 3.0]:
        print(f"\n--- Testing with alpha = {alpha} ---")
        
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.cutoff = 1.2
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two particles at 1 nm
        atoms = []
        
        atom1 = MCAtom()
        atom1.x = 5.0
        atom1.y = 5.0
        atom1.z = 5.0
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 6.0
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
        
        # Initialize PME with this alpha
        mesh_size = [32, 32, 32]
        setPMEParameters(alpha, mesh_size, 4, 1e-6)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        
        # Compute
        computeSystemEnergyPME(state)
        
        # Check self energy
        self_energy = state.ewald_energy.get('self', 0.0)
        expected_self = -2 * alpha / math.sqrt(math.pi) * 138.935456  # 2 particles
        
        print(f"  Self energy: {self_energy:.2f}")
        print(f"  Expected: {expected_self:.2f}")
        print(f"  Ratio: {self_energy / expected_self:.3f}")
        
        # The self energy should scale with alpha
        # If the ratio is always 1.0, then alpha is being set correctly
        # If the ratio is constant but not 1.0, then there's a scaling issue
    
    # Now check if the erfc table issue is alpha-related
    print("\n\n--- Checking if erfcApprox scales with alpha ---")
    
    r = 1.0
    for alpha in [1.0, 2.0, 3.0]:
        expected_erfc = math.erfc(alpha * r)
        
        # From our previous tests, erfcApprox returns 2.0
        # If it's always 2.0 regardless of alpha, that's the bug
        # If it scales with alpha, then it might be returning 2*alpha/sqrt(pi) or similar
        
        print(f"\nAlpha = {alpha}:")
        print(f"  Expected erfc({alpha}) = {expected_erfc:.6f}")
        print(f"  Special value 2α/√π = {2*alpha/math.sqrt(math.pi):.6f}")


if __name__ == "__main__":
    check_pme_alpha()