# tests/simulation/energyPGP/test_pgp_pme_cutoff.py
"""
Test to check if pme_params.cutoff is set when using PGP.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


def test_pme_cutoff_in_pgp():
    """Test if pme_params.cutoff is properly set when using PGP."""
    
    print("\n=== Test: PME Cutoff in PGP ===")
    
    # Create minimal system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Two atoms
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 2.8
    atom2.y = 2.5
    atom2.z = 2.5
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
    
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    print("Test 1: Only setPGPParameters (no initializePMEParameters)")
    print(f"System cutoff: {state.info.cutoff} nm")
    
    # Only set PGP parameters
    setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=state.info.cutoff,
        potentialGridSize=mesh_size,
        splineOrder=4,
        tolerance=1e-5
    )
    
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate
    computeSystemEnergyPGP(state)
    real_space = state.ewald_energy.get('real_space', 0.0)
    
    print(f"PGP real-space: {real_space:.6f} kJ/mol")
    
    if abs(real_space) < 0.1:
        print("❌ Zero - pme_params.cutoff likely not set")
        print("This happens because setPGPParameters calls setPMEParameters")
        print("but setPMEParameters doesn't set cutoff - only initializePMETables does!")
    else:
        print("✅ Non-zero")
    
    print("\n\nTest 2: Call initializePMEParameters after setPGPParameters")
    
    # Now initialize PME parameters
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Calculate again
    computeSystemEnergyPGP(state)
    real_space = state.ewald_energy.get('real_space', 0.0)
    
    print(f"PGP real-space: {real_space:.6f} kJ/mol")
    
    if abs(real_space) < 0.1:
        print("❌ Still zero - the issue is deeper than cutoff")
    else:
        print("✅ Non-zero - fixed by initializing PME")


if __name__ == "__main__":
    test_pme_cutoff_in_pgp()