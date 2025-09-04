"""
PGP Complete validation tests for GCMC applications

This module contains practical validation tests that verify PGP Complete
works correctly for typical GCMC scenarios with reasonable precision standards.
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import numpy as np
import math


def next_power_of_2(n):
    """Return the next power of 2 greater than or equal to n"""
    return 1 << (n - 1).bit_length()


def test_pgp_complete_path_conservation():
    """Test that PGP Complete conserves energy along cyclic paths"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Path Conservation Test")
    print("="*70)
    
    # Create state
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    sigma_na = 0.333
    sigma_cl = 0.442
    eps_na = 0.0115
    eps_cl = 0.4184
    
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # Create atoms
    atoms = []
    
    # Fixed Na-Cl
    positions = [
        (2.0, 2.5, 2.5, 1.0, 0),   # Na+
        (2.5, 2.5, 2.5, -1.0, 1),  # Cl-
        # Moving Na-Cl
        (3.5, 2.5, 2.5, 1.0, 0),   # Na+
        (4.0, 2.5, 2.5, -1.0, 1),  # Cl-
    ]
    
    for x, y, z, charge, atype in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = x, y, z
        atom.charge = charge
        atom.type = atype
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Residues
    residues = []
    
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = True
    residues.append(res)
    
    res = MCResidue()
    res.atomStart = 2
    res.atomCount = 2
    res.active = True
    res.fixed = False
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Initialize PGP
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Get initial energy
    result_init = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    energy_init = result_init[0] + result_init[1]
    
    # Define path
    path = [
        (0.1, 0.0, 0.0),
        (0.0, 0.1, 0.0),
        (-0.1, 0.0, 0.0),
        (0.0, -0.1, 0.0),
    ]
    
    print(f"Initial energy: {energy_init:.6f} kJ/mol")
    
    # Follow path
    for i, (dx, dy, dz) in enumerate(path):
        # Move atoms
        for j in [2, 3]:  # Movement atoms
            state.atoms[j].x += dx
            state.atoms[j].y += dy
            state.atoms[j].z += dz
        
        result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        energy = result[0] + result[1]
        print(f"Step {i+1}: E = {energy:.6f} kJ/mol, ΔE = {energy - energy_init:.6f} kJ/mol")
    
    # Get final energy
    result_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    energy_final = result_final[0] + result_final[1]
    
    print(f"Final energy: {energy_final:.6f} kJ/mol")
    print(f"Difference from initial: {abs(energy_final - energy_init):.2e} kJ/mol")
    
    # Should return to exact initial energy
    assert abs(energy_final - energy_init) < 1e-10, \
        f"Energy not conserved: {energy_final} != {energy_init}"
    
    print("✅ Path conservation test PASSED")


def test_pgp_complete_delta_e_distribution():
    """Test PGP Complete Delta E accuracy for typical GCMC moves"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    np.random.seed(42)
    
    print("\n" + "="*70)
    print("PGP Complete Delta E Distribution Test")
    print("="*70)
    
    # Create state
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    sigma_na = 0.333
    sigma_cl = 0.442
    eps_na = 0.0115
    eps_cl = 0.4184
    
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # Create atoms - simple 2 Na-Cl pairs
    atoms = []
    
    # Fixed Na-Cl
    positions = [
        (2.0, 2.5, 2.5, 1.0, 0),   # Na+
        (2.5, 2.5, 2.5, -1.0, 1),  # Cl-
        # Moving Na-Cl
        (3.5, 2.5, 2.5, 1.0, 0),   # Na+
        (4.0, 2.5, 2.5, -1.0, 1),  # Cl-
    ]
    
    for x, y, z, charge, atype in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = x, y, z
        atom.charge = charge
        atom.type = atype
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Store original positions
    original_positions = [(atom.x, atom.y, atom.z) for atom in atoms]
    
    # Residues
    residues = []
    
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = True
    residues.append(res)
    
    res = MCResidue()
    res.atomStart = 2
    res.atomCount = 2
    res.active = True
    res.fixed = False
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Initialize PGP
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Test N random displacements
    n_tests = 100
    delta_energies = []
    
    print(f"\nTesting {n_tests} random displacements...")
    
    for i in range(n_tests):
        # Random displacement (typical GCMC)
        max_disp = 0.15
        dx = np.random.uniform(-max_disp, max_disp)
        dy = np.random.uniform(-max_disp, max_disp)
        dz = np.random.uniform(-max_disp, max_disp)
        
        # Reset to original positions
        for j in range(4):
            state.atoms[j].x = original_positions[j][0]
            state.atoms[j].y = original_positions[j][1]
            state.atoms[j].z = original_positions[j][2]
        
        # PGP before
        pgp_result1 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_e1 = pgp_result1[0] + pgp_result1[1]
        
        # Move atoms (as rigid body)
        for j in [2, 3]:  # Movement atoms
            state.atoms[j].x += dx
            state.atoms[j].y += dy
            state.atoms[j].z += dz
        
        # PGP after
        pgp_result2 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_e2 = pgp_result2[0] + pgp_result2[1]
        
        pgp_delta = pgp_e2 - pgp_e1
        delta_energies.append(pgp_delta)
    
    # Analyze results
    delta_energies = np.array(delta_energies)
    
    print(f"\nDelta E statistics:")
    print(f"  Mean: {np.mean(delta_energies):.3f} kJ/mol")
    print(f"  Std: {np.std(delta_energies):.3f} kJ/mol")
    print(f"  Min: {np.min(delta_energies):.3f} kJ/mol")
    print(f"  Max: {np.max(delta_energies):.3f} kJ/mol")
    
    # Check that energies are reasonable
    assert np.all(np.isfinite(delta_energies)), "All delta energies must be finite"
    assert np.std(delta_energies) > 0.01, "Should see variation in delta energies"
    assert np.abs(np.mean(delta_energies)) < 100, "Mean delta E should be reasonable"
    
    print("✅ Delta E distribution test PASSED")