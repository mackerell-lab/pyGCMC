# tests/simulation/energyPGP/test_pgp_erfc_debug.py
"""
Debug test to check if PGP erfc table is properly initialized.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP, computeSystemEnergyPME


def test_pgp_erfc_table():
    """Test if PGP erfc table is properly initialized."""
    
    print("\n=== PGP erfc Table Debug Test ===")
    
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
    
    # Test different initialization orders
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    print("\n--- Test 1: Initialize PME first, then PGP ---")
    
    # Initialize PME first (this creates the erfc tables)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Then set PGP parameters (should copy the tables)
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
    
    # Calculate with both methods
    print("\nPME calculation:")
    computeSystemEnergyPME(state)
    pme_real = state.ewald_energy.get('real_space', 0.0)
    print(f"PME real-space: {pme_real:.6f} kJ/mol")
    
    print("\nPGP calculation:")
    computeSystemEnergyPGP(state)
    pgp_real = state.ewald_energy.get('real_space', 0.0)
    print(f"PGP real-space: {pgp_real:.6f} kJ/mol")
    
    # Manual calculation
    r = 0.3
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected = -erfc_val / r * kC
    print(f"\nExpected real-space: {expected:.6f} kJ/mol")
    print(f"erfc({alpha}*{r}) = {erfc_val:.6f}")
    
    if abs(pgp_real) < 0.1:
        print("\n❌ PGP real-space is still zero!")
        print("This suggests pgp_params.erfcTable is not properly initialized")
    else:
        print("\n✅ PGP real-space is non-zero")


def test_compare_initialization_sequences():
    """Compare different ways to initialize PGP."""
    
    print("\n=== Compare PGP Initialization Sequences ===")
    
    # Helper to create identical systems
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
    
    # Try different orders
    sequences = [
        ("PME->PGP", lambda s: [
            initializePMEParameters(s.info.cutoff, s.info.box, alpha),
            setPGPParameters(alpha, mesh_size, s.info.cutoff, mesh_size, 4, 1e-5)
        ]),
        ("PGP->PME", lambda s: [
            setPGPParameters(alpha, mesh_size, s.info.cutoff, mesh_size, 4, 1e-5),
            initializePMEParameters(s.info.cutoff, s.info.box, alpha)
        ])
    ]
    
    for name, init_func in sequences:
        print(f"\n--- {name} ---")
        state = create_system()
        
        # Run initialization sequence
        init_func(state)
        precomputeGridPotential(state, fixed_only=True)
        
        # Calculate energy
        computeSystemEnergyPGP(state)
        real_space = state.ewald_energy.get('real_space', 0.0)
        
        print(f"Real-space energy: {real_space:.6f} kJ/mol")
        
        if abs(real_space) < 0.1:
            print("❌ Zero - initialization failed")
        else:
            print("✅ Non-zero - initialization succeeded")


if __name__ == "__main__":
    test_pgp_erfc_table()
    test_compare_initialization_sequences()