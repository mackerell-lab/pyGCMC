# tests/simulation/movementInsert/insertion_deletion_helpers.py
"""
Helper functions for GCMC insertion/deletion tests
"""

import pytest
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²
kB = 0.008314463  # Boltzmann constant in kJ/mol/K


def create_empty_system():
    """Create an empty system with no molecules"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom type but no atoms
    state.atomTypes.get_or_add_type("ION")
    
    # Set up force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [0.3]
    
    state.atoms = []
    state.residues = []
    state.activeAtomCount = 0
    state.activeResidueCount = 0
    
    return state


def create_system_with_one_molecule():
    """Create system with a single molecule"""
    state = create_empty_system()
    
    # Add one atom/molecule
    atom = pygcmc.MCAtom()
    atom.x = 2.5
    atom.y = 2.5
    atom.z = 2.5
    atom.charge = 1.0
    atom.type = 0
    
    state.atoms = [atom]
    state.activeAtomCount = 1
    
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = 0
    res.atomCount = 1
    res.type = 0
    
    state.residues = [res]
    state.activeResidueCount = 1
    
    return state


def create_system_with_molecules(n):
    """Create system with n molecules in a regular arrangement"""
    state = create_empty_system()
    
    if n == 0:
        return state
    
    # Arrange molecules in a line with 0.5 nm spacing
    atoms = []
    residues = []
    
    for i in range(n):
        atom = pygcmc.MCAtom()
        atom.x = 1.0 + i * 0.5
        atom.y = 2.5
        atom.z = 2.5
        atom.charge = 1.0 if i % 2 == 0 else -1.0  # Alternate charges
        atom.type = 0
        atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n
    state.activeResidueCount = n
    
    return state


def create_framework_system():
    """Create a system with framework atoms forming a cavity"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    framework_type = state.atomTypes.get_or_add_type("FRAMEWORK")
    guest_type = state.atomTypes.get_or_add_type("GUEST")
    
    # Set up force field (2 types)
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1  # Only guest moves
    state.forcefield.ljEps = [1.0, 0.5, 0.5, 1.0]  # Guest-guest, guest-framework, framework-framework
    state.forcefield.ljSigma = [0.3, 0.35, 0.35, 0.3]
    
    # Create framework atoms in a cubic arrangement (cavity in center)
    atoms = []
    residues = []
    atom_idx = 0
    
    # Framework atoms at corners of a cube
    corners = [
        (1.0, 1.0, 1.0), (1.0, 1.0, 4.0),
        (1.0, 4.0, 1.0), (1.0, 4.0, 4.0),
        (4.0, 1.0, 1.0), (4.0, 1.0, 4.0),
        (4.0, 4.0, 1.0), (4.0, 4.0, 4.0)
    ]
    
    for x, y, z in corners:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = -0.1  # Slightly negative to attract guest
        atom.type = framework_type
        atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 1  # Framework type
        residues.append(res)
        atom_idx += 1
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state


def create_guest_molecule(x, y, z, charge=0.5, atom_type=0):
    """Create a guest molecule at specified position"""
    atom = pygcmc.MCAtom()
    atom.x = x
    atom.y = y
    atom.z = z
    atom.charge = charge
    atom.type = atom_type
    return atom


def add_molecule_to_state(state, molecule, movement_type=False):
    """Add a molecule to an existing state"""
    # Create a copy of the state
    new_state = pygcmc.MCState()
    # Deep copy the info structure to avoid shared pointers
    new_state.info.box = list(state.info.box)
    new_state.info.cutoff = state.info.cutoff
    new_state.atomTypes = state.atomTypes
    new_state.forcefield = state.forcefield
    
    # Copy existing atoms and residues
    atoms = list(state.atoms)
    residues = list(state.residues)
    
    # Add new molecule
    atom_idx = len(atoms)
    atoms.append(molecule)
    
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = atom_idx
    res.atomCount = 1
    res.type = 0 if movement_type else 1
    residues.append(res)
    
    new_state.atoms = atoms
    new_state.residues = residues
    new_state.activeAtomCount = len(atoms)
    new_state.activeResidueCount = len(residues)
    
    return new_state


def get_system_energy(state):
    """Get total system energy, correcting for double counting"""
    total = 0.0
    for res in state.residues:
        if res.active:
            total += res.energy_vdw + res.energy_elec
    return total / 2.0  # Divide by 2 to correct for double counting
