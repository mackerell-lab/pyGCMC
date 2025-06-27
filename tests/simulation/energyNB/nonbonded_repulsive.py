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
    
    # Expected energy based on LJ formula
    expected_energy = 4.0 * 1.0 * ((1.0/0.5)**12 - (1.0/0.5)**6)
    expected_energy = 4.0 * (2**12 - 2**6)  # 4 * (4096 - 64) = 4 * 4032 = 16128
    
    # Check result - allow for numerical precision
    assert abs(energy - expected_energy) < 1e-6, f"Expected {expected_energy}, got {energy}"


def test_combined_interaction():
    """Test nonbonded energy calculation for combined LJ and electrostatic interactions.
    
    Setup:
    - Two residues: one movement (+0.5 charge) and one fixed (-0.5 charge)
    - LJ parameters: eps = 1.0, sigma = 1.0
    - r = 1.122462 (near LJ minimum)
    
    Expected:
    - Both LJ and electrostatic contributions
    - Net attractive interaction due to opposite charges and LJ minimum
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2  # Two types: movement and fixed
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Set both LJ and charge interactions
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # LJ interaction
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # LJ size parameter
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms with charges and positioned near LJ minimum
    # First atom (for movement residue)
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.5  # Small positive charge
    atom1.type = 0  # Movement type
    
    # Second atom (for fixed residue)
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.122462  # Distance at LJ minimum (2^(1/6) * sigma)
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -0.5  # Small negative charge
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
    
    # Check that we get an attractive (negative) total interaction
    assert energy < 0, "Expected negative (attractive) total energy for combined interaction at LJ minimum"