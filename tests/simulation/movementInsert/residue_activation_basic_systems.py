# tests/simulation/movementInsert/residue_activation_systems.py
"""
Helper functions for creating test systems for residue activation tests
"""

import pytest
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²


def create_two_residue_system():
    """Create a system with 2 residues"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    type_idx = state.atomTypes.get_or_add_type("ION")
    
    # Set up force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]  # epsilon = 1.0 kJ/mol
    state.forcefield.ljSigma = [0.3]  # sigma = 0.3 nm
    
    # Create 2 atoms
    atom0 = pygcmc.MCAtom()
    atom0.x = 1.0
    atom0.y = 1.0
    atom0.z = 1.0
    atom0.charge = 1.0
    atom0.type = type_idx
    
    atom1 = pygcmc.MCAtom()
    atom1.x = 2.0
    atom1.y = 1.0
    atom1.z = 1.0
    atom1.charge = -1.0
    atom1.type = type_idx
    
    state.atoms = [atom0, atom1]
    state.activeAtomCount = 2
    
    # Create residues
    res0 = pygcmc.MCResidue()
    res0.active = True
    res0.atomStart = 0
    res0.atomCount = 1
    res0.type = 0
    
    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.atomStart = 1
    res1.atomCount = 1
    res1.type = 0
    
    state.residues = [res0, res1]
    state.activeResidueCount = 2
    
    return state


def create_three_residue_system():
    """Create a system with 3 residues (same as 2-residue system + 1 new)"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    type_idx = state.atomTypes.get_or_add_type("ION")
    
    # Set up force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [0.3]
    
    # Create 3 atoms (first 2 same as before, plus 1 new)
    atom0 = pygcmc.MCAtom()
    atom0.x = 1.0
    atom0.y = 1.0
    atom0.z = 1.0
    atom0.charge = 1.0
    atom0.type = type_idx
    
    atom1 = pygcmc.MCAtom()
    atom1.x = 2.0
    atom1.y = 1.0
    atom1.z = 1.0
    atom1.charge = -1.0
    atom1.type = type_idx
    
    # New atom
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.5
    atom2.y = 1.5
    atom2.z = 1.0
    atom2.charge = 0.5
    atom2.type = type_idx
    
    state.atoms = [atom0, atom1, atom2]
    state.activeAtomCount = 3
    
    # Create residues
    res0 = pygcmc.MCResidue()
    res0.active = True
    res0.atomStart = 0
    res0.atomCount = 1
    res0.type = 0
    
    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.atomStart = 1
    res1.atomCount = 1
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.active = True
    res2.atomStart = 2
    res2.atomCount = 1
    res2.type = 0
    
    state.residues = [res0, res1, res2]
    state.activeResidueCount = 3
    
    return state


def create_mixed_residue_system():
    """Create a system with movement and fixed residues"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types (0 = movement, 1 = fixed)
    movement_type = state.atomTypes.get_or_add_type("GUEST")
    fixed_type = state.atomTypes.get_or_add_type("FRAMEWORK")
    
    # Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1  # Only first type is movement
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]
    state.forcefield.ljSigma = [0.3, 0.3, 0.3, 0.3]
    
    # Create atoms
    # Fixed atom (framework)
    atom_fixed = pygcmc.MCAtom()
    atom_fixed.x = 1.0
    atom_fixed.y = 1.0
    atom_fixed.z = 1.0
    atom_fixed.charge = 0.0
    atom_fixed.type = fixed_type
    
    # Movement atom (guest)
    atom_movement = pygcmc.MCAtom()
    atom_movement.x = 1.5
    atom_movement.y = 1.0
    atom_movement.z = 1.0
    atom_movement.charge = 1.0
    atom_movement.type = movement_type
    
    state.atoms = [atom_fixed, atom_movement]
    state.activeAtomCount = 2
    
    # Create residues
    res_fixed = pygcmc.MCResidue()
    res_fixed.active = True
    res_fixed.atomStart = 0
    res_fixed.atomCount = 1
    res_fixed.type = 1  # Fixed type
    
    res_movement = pygcmc.MCResidue()
    res_movement.active = True
    res_movement.atomStart = 1
    res_movement.atomCount = 1
    res_movement.type = 0  # Movement type
    
    state.residues = [res_fixed, res_movement]
    state.activeResidueCount = 2
    
    return state


def create_system_with_inactive_residue():
    """Create a system with 3 residues, one inactive"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    type_idx = state.atomTypes.get_or_add_type("ION")
    
    # Set up force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [0.3]
    
    # Create 3 atoms
    atom0 = pygcmc.MCAtom()
    atom0.x = 1.0
    atom0.y = 1.0
    atom0.z = 1.0
    atom0.charge = 1.0
    atom0.type = type_idx
    
    atom1 = pygcmc.MCAtom()
    atom1.x = 2.0
    atom1.y = 1.0
    atom1.z = 1.0
    atom1.charge = -1.0
    atom1.type = type_idx
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.5
    atom2.y = 1.5
    atom2.z = 1.0
    atom2.charge = 0.5
    atom2.type = type_idx
    
    state.atoms = [atom0, atom1, atom2]
    state.activeAtomCount = 2  # Only 2 active
    
    # Create residues with third one inactive
    res0 = pygcmc.MCResidue()
    res0.active = True
    res0.atomStart = 0
    res0.atomCount = 1
    res0.type = 0
    
    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.atomStart = 1
    res1.atomCount = 1
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.active = False  # Inactive
    res2.atomStart = 2
    res2.atomCount = 1
    res2.type = 0
    
    state.residues = [res0, res1, res2]
    state.activeResidueCount = 2
    
    return state


