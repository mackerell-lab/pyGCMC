# tests/simulation/movementInsert/advanced_insertion_molecules.py
"""
Molecule creation and manipulation functions for advanced insertion tests
"""

import pytest
import random
import math
import pygcmc

def create_guest_molecule(x, y, z, charge=0.0):
    """Create a guest atom/molecule"""
    atom = pygcmc.MCAtom()
    atom.x = x
    atom.y = y
    atom.z = z
    atom.charge = charge
    atom.type = 1  # GUEST
    return atom


def create_water_molecule(x, y, z):
    """Create water molecule at position"""
    atoms = []
    
    # Oxygen
    o_atom = pygcmc.MCAtom()
    o_atom.x = x
    o_atom.y = y
    o_atom.z = z
    o_atom.charge = -0.834
    o_atom.type = 2  # OT
    atoms.append(o_atom)
    
    # Hydrogens
    h1_atom = pygcmc.MCAtom()
    h1_atom.x = x + 0.0957
    h1_atom.y = y
    h1_atom.z = z
    h1_atom.charge = 0.417
    h1_atom.type = 3  # HT
    atoms.append(h1_atom)
    
    h2_atom = pygcmc.MCAtom()
    h2_atom.x = x - 0.0757
    h2_atom.y = y + 0.0587
    h2_atom.z = z
    h2_atom.charge = 0.417
    h2_atom.type = 3  # HT
    atoms.append(h2_atom)
    
    return atoms


def insert_single_atom(system, atom):
    """Insert single atom as guest molecule"""
    new_system = pygcmc.MCState()
    new_system.info.box = system.info.box
    new_system.info.cutoff = system.info.cutoff
    new_system.atomTypes = system.atomTypes
    new_system.forcefield = system.forcefield
    
    # Copy existing atoms
    new_system.atoms = []
    for old_atom in system.atoms:
        new_atom = pygcmc.MCAtom()
        new_atom.x = old_atom.x
        new_atom.y = old_atom.y
        new_atom.z = old_atom.z
        new_atom.charge = old_atom.charge
        new_atom.type = old_atom.type
        new_system.atoms.append(new_atom)
    
    # Add new atom
    new_system.atoms.append(atom)
    
    # Copy existing residues
    new_system.residues = []
    for res in system.residues:
        new_res = pygcmc.MCResidue()
        new_res.active = res.active
        new_res.atomStart = res.atomStart
        new_res.atomCount = res.atomCount
        new_res.type = res.type
        new_system.residues.append(new_res)
    
    # Add new residue
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = len(system.atoms)
    res.atomCount = 1
    res.type = 0  # Guest type
    new_system.residues.append(res)
    
    new_system.activeAtomCount = len(new_system.atoms)
    new_system.activeResidueCount = len(new_system.residues)
    
    return new_system


def insert_molecule(system, atoms):
    """Insert molecule (list of atoms)"""
    new_system = pygcmc.MCState()
    new_system.info.box = system.info.box
    new_system.info.cutoff = system.info.cutoff
    new_system.atomTypes = system.atomTypes
    new_system.forcefield = system.forcefield
    
    # Copy existing atoms
    new_system.atoms = []
    for old_atom in system.atoms:
        new_atom = pygcmc.MCAtom()
        new_atom.x = old_atom.x
        new_atom.y = old_atom.y
        new_atom.z = old_atom.z
        new_atom.charge = old_atom.charge
        new_atom.type = old_atom.type
        new_system.atoms.append(new_atom)
    
    atom_start = len(new_system.atoms)
    
    # Add new atoms
    for atom in atoms:
        new_system.atoms.append(atom)
    
    # Copy existing residues
    new_system.residues = []
    for res in system.residues:
        new_res = pygcmc.MCResidue()
        new_res.active = res.active
        new_res.atomStart = res.atomStart
        new_res.atomCount = res.atomCount
        new_res.type = res.type
        new_system.residues.append(new_res)
    
    # Add new residue
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = atom_start
    res.atomCount = len(atoms)
    res.type = 0  # Guest type
    new_system.residues.append(res)
    
    new_system.activeAtomCount = len(new_system.atoms)
    new_system.activeResidueCount = len(new_system.residues)
    
    return new_system


def calculate_system_energy(state):
    """Calculate total system energy"""
    total = 0.0
    for res in state.residues:
        if res.active:
            total += res.energy_vdw + res.energy_elec
    return total / 2.0 if len(state.residues) > 1 else total
