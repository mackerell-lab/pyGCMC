# tests/simulation/energyNB/nonbonded_repulsive.py
"""Non-bonded repulsive and combined interaction tests."""

import pytest
import pygcmc
import math


def test_repulsive_interaction():
    """Test nonbonded energy calculation for a repulsive interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 0.5 (distance between atoms - closer than sigma)
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    V = 4 * 1.0 * [(1.0/0.5)¹² - (1.0/0.5)⁶]
    V = 4 * 1.0 * [2¹² - 2⁶] = 4 * (4096 - 64) = 4 * 4032 = 16128
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2  # Two types: movement and fixed
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters for repulsive interaction
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms
    # First atom (for movement residue)
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.0  # No electrostatic interaction
    atom1.type = 0  # Movement type
    
    # Second atom (for fixed residue)
    atom2 = pygcmc.MCAtom()
    atom2.x = 0.5  # At distance 0.5 from first atom (closer than sigma)
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 0.0  # No electrostatic interaction
    atom2.type = 1  # Fixed type
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    expected_energy = 16128.0  # 4 * (4096 - 64)
    assert abs(energy - expected_energy) < 1.0, "Expected strong repulsive energy"


def test_combined_interaction():
    """Test both LJ and Coulomb interactions together.
    
    Setup:
    - Two residues with both LJ and charge interactions
    - eps = 1.0 kJ/mol, sigma = 1.0 nm
    - q1 = +0.5e, q2 = -0.5e
    - r = 1.0 nm
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] = 0 kJ/mol (at r = sigma)
    V_C = k_c * q1*q2/r = 138.935458 * 0.5 * (-0.5) / 1.0 = -34.734 kJ/mol
    Each residue gets this full energy (not divided by 2)
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms with both LJ and charge interactions
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.5  # Half positive charge
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # Distance of 1.0 nm
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -0.5  # Half negative charge
    atom2.type = 1
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    movement_info.resName = "MOV"
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    
    # Get actual energies
    actual_vdw = state.residues[0].energy_vdw
    actual_elec = state.residues[0].energy_elec
    actual_total = actual_vdw + actual_elec
    
    # Expected energies
    COULOMB = 138.935458  # kJ·nm/mol/e²
    expected_vdw = 0.0  # kJ/mol (at r = sigma)
    expected_elec = COULOMB * 0.5 * (-0.5) / 1.0  # kJ/mol
    expected_total = expected_vdw + expected_elec
    
    # Check results with appropriate tolerance
    rel_tol = 1e-5
    assert abs(actual_vdw - expected_vdw) < 1e-6, \
           f"VDW energy mismatch: expected {expected_vdw}, got {actual_vdw}"
    assert abs((actual_elec - expected_elec) / expected_elec) < rel_tol, \
           f"Electrostatic energy mismatch: expected {expected_elec}, got {actual_elec}"
    assert abs((actual_total - expected_total) / expected_total) < rel_tol, \
           f"Total energy mismatch: expected {expected_total}, got {actual_total}"