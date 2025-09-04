# tests/simulation/energyPGP/test_pgp_initialization_order.py
"""
Test different initialization orders to fix PGP real-space calculation.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters
from pygcmc import precomputeGridPotential, computeSystemEnergyPGP
import os


def test_pgp_with_correct_initialization():
    """Test PGP with correct initialization order."""
    
    print("\n=== Test: Correct PGP Initialization Order ===")
    
    # Enable debug
    os.environ['PYGCMC_DEBUG'] = '1'
    
    # Create system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Two atoms in different residues
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 2.8  # 0.3 nm away
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
    
    print("System: +1 and -1 charges at 0.3 nm distance")
    
    # METHOD 1: Initialize PME then PGP
    print("\nMETHOD 1: Initialize PME then PGP")
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    # Initialize PME first
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Then set PGP parameters
    setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=state.info.cutoff,
        potentialGridSize=mesh_size,
        splineOrder=4,
        tolerance=1e-5
    )
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy
    computeSystemEnergyPGP(state)
    
    # Get results
    real_space = state.ewald_energy.get('real_space', 0.0)
    reciprocal = state.ewald_energy.get('reciprocal', 0.0)
    self_energy = state.ewald_energy.get('self', 0.0)
    total = state.ewald_energy.get('total', 0.0)
    
    print(f"\nEnergy components:")
    print(f"  Real-space: {real_space:.6f} kJ/mol")
    print(f"  Reciprocal: {reciprocal:.6f} kJ/mol")
    print(f"  Self: {self_energy:.6f} kJ/mol")
    print(f"  Total: {total:.6f} kJ/mol")
    
    # Expected real-space
    r = 0.3
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected = -erfc_val / r * kC
    
    print(f"\nExpected real-space: {expected:.6f} kJ/mol")
    
    if abs(real_space) < 0.1:
        print("\n❌ Still zero with initializePGPParameters!")
    else:
        print(f"\n✅ Real-space is non-zero: {real_space:.6f} kJ/mol")
        rel_error = abs((real_space - expected) / expected)
        if rel_error < 0.01:
            print("✅ Matches expected value!")
        else:
            print(f"⚠️  Differs from expected by {rel_error*100:.1f}%")


def test_pgp_initialization_methods():
    """Compare different initialization methods."""
    
    print("\n=== Compare PGP Initialization Methods ===")
    
    # Create identical systems
    def create_system():
        state = MCState()
        state.info.box = [5.0, 5.0, 5.0]
        state.info.cutoff = 2.0
        
        state.forcefield.numTotalTypes = 1
        state.forcefield.numMovementTypes = 1
        state.forcefield.ljEps = [0.0]
        state.forcefield.ljSigma = [0.3]
        
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
        
        return state
    
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    # Test 1: setPGPParameters THEN initializePMEParameters
    print("\nTest 1: setPGPParameters -> initializePMEParameters")
    state1 = create_system()
    setPGPParameters(alpha, mesh_size, state1.info.cutoff, mesh_size, 4, 1e-5)
    initializePMEParameters(state1.info.cutoff, state1.info.box, alpha)
    precomputeGridPotential(state1, fixed_only=True)
    computeSystemEnergyPGP(state1)
    print(f"Real-space energy: {state1.ewald_energy.get('real_space', 0.0):.6f} kJ/mol")
    
    # Test 2: initializePMEParameters THEN setPGPParameters
    print("\nTest 2: initializePMEParameters -> setPGPParameters")
    state2 = create_system()
    initializePMEParameters(state2.info.cutoff, state2.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state2.info.cutoff, mesh_size, 4, 1e-5)
    precomputeGridPotential(state2, fixed_only=True)
    computeSystemEnergyPGP(state2)
    print(f"Real-space energy: {state2.ewald_energy.get('real_space', 0.0):.6f} kJ/mol")
    
    # Test 3: Skip the non-existent function
    print("\nTest 3: Skipped (initializePGPParameters not available in Python bindings)")


if __name__ == "__main__":
    test_pgp_with_correct_initialization()
    test_pgp_initialization_methods()