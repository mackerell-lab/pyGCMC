# tests/simulation/energyNB/nonbonded_attractive.py
"""Non-bonded attractive and electrostatic interaction tests."""

import pytest
import pygcmc
import math


def test_attractive_interaction():
    """Test nonbonded energy calculation for an attractive interaction between movement and fixed residues.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with zero charge (only vdw interaction)
    - Both residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (distance between atoms)
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    V = 4 * 1.0 * [(1.0/1.0)¹² - (1.0/1.0)⁶]
    V = 4 * 1.0 * (1 - 1) = 0
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2  # Two types: movement and fixed
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters for attractive interaction
    # For type 0:
    #   - interaction with type 0: index = 0 * 2 + 0 = 0
    #   - interaction with type 1: index = 0 * 2 + 1 = 1
    # For type 1:
    #   - interaction with type 0: index = 1 * 2 + 0 = 2
    #   - interaction with type 1: index = 1 * 2 + 1 = 3
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
    atom2.x = 1.0  # At distance 1.0 from first atom
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
    assert abs(energy) < 1e-6, "Expected zero interaction energy at r = sigma"


def test_electrostatic_interaction():
    """Test nonbonded energy calculation for electrostatic interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with opposite charges (+1 and -1)
    - Both residues are active
    - Force field parameters:
        eps = 0.0 (no vdw interaction)
        r = 1.0 nm
    
    Expected:
    V = k_c * q1*q2/r where:
    - k_c = 138.935458 kJ·nm/mol/e^2 (Coulomb constant)
    - q1 = +1e, q2 = -1e
    - r = 1.0 nm
    Therefore:
    V = 138.935458 * (+1) * (-1) / 1.0 = -138.935458 kJ/mol
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field (no vdw interaction)
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0, 0.0, 0.0, 0.0]    # Complete 2x2 matrix with eps = 0.0
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix with sigma = 1.0
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms with opposite charges
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0  # Positive charge
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # Distance of 1.0 nm
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -1.0  # Negative charge
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
    movement_info.resName = "MOV"  # Add residue name
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Expected energy with Coulomb constant
    COULOMB = 138.935458  # kJ·nm/mol/e^2
    expected_energy = -COULOMB  # k_c * (+1) * (-1) / 1.0
    
    # Check result with appropriate tolerance for single-precision float
    assert abs(energy - expected_energy) < 1e-5, f"Expected energy {expected_energy}, got {energy}"