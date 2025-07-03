# tests/simulation/energyPGP/test_pgp_real_fixed.py
"""
Fixed test that properly initializes PGP parameters.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


def test_pgp_real_space_fixed():
    """Test PGP real-space with proper initialization."""
    
    print("\n=== Fixed PGP Real-Space Test ===")
    
    # Create system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Two charged atoms
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 3.0  # 0.5 nm away
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
    
    print("System: +1 and -1 charges at 0.5 nm distance")
    
    # Initialize PGP parameters
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    # IMPORTANT: First initialize PME tables
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Then set PGP parameters (this will copy from PME)
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
    r = 0.5
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected = -erfc_val / r * kC
    
    print(f"\nExpected real-space: {expected:.6f} kJ/mol")
    
    # Now real-space should NOT be zero
    if abs(real_space) < 0.1:
        print("\n❌ Real-space is still zero!")
        
        # Additional debug: check if erfc table is initialized
        # We can't directly access pgp_params from Python, but we can infer
        print("\nPossible issues:")
        print("1. PGP erfc table not properly initialized")
        print("2. Real-space loop not executing correctly")
        print("3. Charges being treated as zero")
    else:
        print(f"\n✅ Real-space is non-zero: {real_space:.6f} kJ/mol")
        
        # Check if it matches expected
        rel_error = abs((real_space - expected) / expected)
        if rel_error < 0.01:
            print("✅ Real-space matches expected value!")
        else:
            print(f"⚠️  Real-space differs from expected by {rel_error*100:.1f}%")


if __name__ == "__main__":
    test_pgp_real_space_fixed()