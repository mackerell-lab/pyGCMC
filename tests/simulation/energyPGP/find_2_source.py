# tests/simulation/energyPGP/find_2_source.py
"""
Find where the value 2.0 is coming from.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def find_2_source():
    """Test different scenarios to find where 2.0 comes from."""
    print("\n=== Finding Source of 2.0 ===")
    
    # Test 1: Check if it happens without initialization
    print("\nTest 1: No initialization")
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
        atom.x = 5.0 + i
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
    
    # Try without proper initialization
    try:
        computeSystemEnergyPME(state)
        print("  Energy without init: SUCCEEDED (shouldn't happen!)")
        print(f"  Residue energy: {state.residues[0].energy_elec:.2f}")
    except:
        print("  Energy without init: FAILED (expected)")
    
    # Test 2: Initialize with different alpha values
    print("\nTest 2: Different alpha values")
    for alpha in [1.0, 2.0, 2.5, 3.0]:
        state_copy = MCState()
        state_copy.info = state.info
        state_copy.forcefield = state.forcefield
        state_copy.atoms = state.atoms
        state_copy.activeAtomCount = state.activeAtomCount
        state_copy.residues = list(state.residues)  # Copy
        state_copy.activeResidueCount = state.activeResidueCount
        
        setPMEParameters(alpha, [32, 32, 32], 4, 1e-6)
        initializePMEParameters(state_copy.info.cutoff, state_copy.info.box, alpha)
        computeSystemEnergyPME(state_copy)
        
        res_energy = state_copy.residues[0].energy_elec
        # Infer erfcApprox value
        # res_energy = pair_energy * COULOMB² where pair_energy = -1 * erfcApprox(1) / 1
        # So: erfcApprox(1) = -res_energy / COULOMB²
        inferred_erfc = -res_energy / (138.935456**2)
        
        print(f"  Alpha={alpha}: Residue energy={res_energy:.2f}, Inferred erfcApprox={inferred_erfc:.6f}")
    
    # Test 3: The special case
    print("\nTest 3: Check if 2.0 comes from a specific formula")
    alpha = 2.5
    print(f"  2*alpha/sqrt(pi) = {2*alpha/math.sqrt(math.pi):.6f}")
    print(f"  2*1/sqrt(pi) = {2*1.0/math.sqrt(math.pi):.6f}")
    print(f"  2.0 = 2.0")
    
    print("\nConclusion: erfcApprox is returning exactly 2.0 for all inputs")


if __name__ == "__main__":
    find_2_source()