# tests/simulation/movementInsert/basic_insertion_helpers.py
"""
Helper functions for basic insertion tests
"""

import pytest
import random
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²


def create_empty_system():
    """Create an empty system ready for insertions"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Pre-define common atom types
    state.atomTypes.get_or_add_type("OT")   # Water oxygen
    state.atomTypes.get_or_add_type("HT")   # Water hydrogen
    state.atomTypes.get_or_add_type("SOD")  # Sodium
    state.atomTypes.get_or_add_type("CLA")  # Chloride
    
    # Set up force field BEFORE setting atoms
    # 4 atom types: OT, HT, SOD, CLA
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4  # All can move
    
    # LJ parameters - need full 4x4 matrix = 16 values
    # Order: OT-OT, OT-HT, OT-SOD, OT-CLA, HT-OT, HT-HT, HT-SOD, HT-CLA, etc.
    state.forcefield.ljEps = [
        0.6364, 0.0, 0.5, 0.7,      # OT with OT, HT, SOD, CLA
        0.0, 0.0, 0.0, 0.0,         # HT with OT, HT, SOD, CLA
        0.5, 0.0, 0.4, 0.6,         # SOD with OT, HT, SOD, CLA
        0.7, 0.0, 0.6, 0.8          # CLA with OT, HT, SOD, CLA
    ]
    
    state.forcefield.ljSigma = [
        0.3166, 0.0, 0.28, 0.35,    # OT interactions
        0.0, 0.0, 0.0, 0.0,         # HT interactions (no LJ)
        0.28, 0.0, 0.24, 0.32,      # SOD interactions
        0.35, 0.0, 0.32, 0.40       # CLA interactions
    ]
    
    # Now we can set atoms
    state.atoms = []
    state.residues = []
    state.activeAtomCount = 0
    state.activeResidueCount = 0
    
    return state


def create_water_ion_forcefield():
    """Create force field for water and ions"""
    ff = pygcmc.MCForceField()
    
    # 4 atom types: OT, HT, SOD, CLA
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4  # All can move
    
    # LJ parameters (epsilon in kJ/mol, sigma in nm)
    # Order: OT-OT, OT-HT, OT-SOD, OT-CLA, HT-HT, HT-SOD, HT-CLA, SOD-SOD, SOD-CLA, CLA-CLA
    ff.ljEps = [
        0.6364, 0.0, 0.5, 0.7,      # OT interactions
        0.0, 0.0, 0.0,              # HT interactions
        0.4, 0.6,                   # SOD interactions
        0.8                         # CLA interactions
    ]
    
    ff.ljSigma = [
        0.3166, 0.0, 0.28, 0.35,    # OT interactions
        0.0, 0.0, 0.0,              # HT interactions  
        0.24, 0.32,                 # SOD interactions
        0.40                        # CLA interactions
    ]
    
    return ff


def create_ion_forcefield():
    """Create force field for ions only"""
    ff = pygcmc.MCForceField()
    
    # 4 atom types to match atomTypes indices: OT(0), HT(1), SOD(2), CLA(3)
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    
    # LJ parameters - full 4x4 matrix = 16 values
    ff.ljEps = [
        0.0, 0.0, 0.0, 0.0,      # OT interactions (not used)
        0.0, 0.0, 0.0, 0.0,      # HT interactions (not used)
        0.0, 0.0, 0.4, 0.6,      # SOD with OT, HT, SOD, CLA
        0.0, 0.0, 0.6, 0.8       # CLA with OT, HT, SOD, CLA
    ]
    
    ff.ljSigma = [
        0.0, 0.0, 0.0, 0.0,      # OT interactions (not used)
        0.0, 0.0, 0.0, 0.0,      # HT interactions (not used)
        0.0, 0.0, 0.24, 0.32,    # SOD interactions
        0.0, 0.0, 0.32, 0.40     # CLA interactions
    ]
    
    return ff


def create_water_molecule(x, y, z):
    """Create a water molecule (TIP3P-like) at given position"""
    atoms = []
    
    # Oxygen
    o_atom = pygcmc.MCAtom()
    o_atom.x = x
    o_atom.y = y
    o_atom.z = z
    o_atom.charge = -0.834
    o_atom.type = 0  # OT
    atoms.append(o_atom)
    
    # Hydrogen 1 (along x-axis from O)
    h1_atom = pygcmc.MCAtom()
    h1_atom.x = x + 0.0957
    h1_atom.y = y
    h1_atom.z = z
    h1_atom.charge = 0.417
    h1_atom.type = 1  # HT
    atoms.append(h1_atom)
    
    # Hydrogen 2 (104.5 degree angle)
    angle = math.radians(104.5)
    h2_atom = pygcmc.MCAtom()
    h2_atom.x = x + 0.0957 * math.cos(angle)
    h2_atom.y = y + 0.0957 * math.sin(angle)
    h2_atom.z = z
    h2_atom.charge = 0.417
    h2_atom.type = 1  # HT
    atoms.append(h2_atom)
    
    return atoms


def insert_molecule(system, molecule_atoms):
    """Insert a molecule (list of atoms) into the system"""
    # Create a new system to avoid modifying the original
    new_system = pygcmc.MCState()
    new_system.info = system.info
    new_system.atomTypes = system.atomTypes
    new_system.forcefield = system.forcefield
    
    atom_start = len(system.atoms)
    
    # Build new atoms list
    new_atoms = []
    
    # Copy existing atoms
    for atom in system.atoms:
        new_atoms.append(atom)
    
    # Add new molecule atoms
    for atom in molecule_atoms:
        new_atoms.append(atom)
    
    # Build new residues list
    new_residues = []
    
    # Copy existing residues
    for res in system.residues:
        new_residues.append(res)
    
    # Add new residue
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = atom_start
    res.atomCount = len(molecule_atoms)
    res.type = 0  # Movement type
    new_residues.append(res)
    
    # Update new system
    new_system.atoms = new_atoms
    new_system.residues = new_residues
    new_system.activeAtomCount = len(new_atoms)
    new_system.activeResidueCount = len(new_residues)
    
    return new_system


def insert_single_atom_molecule(system, atom):
    """Insert a single atom as a molecule"""
    return insert_molecule(system, [atom])


def create_water_system(n_molecules):
    """Create a system with n water molecules"""
    system = create_empty_system()
    
    # Place molecules in a grid pattern
    grid_size = int(math.ceil(n_molecules ** (1/3)))
    spacing = 4.0 / grid_size
    
    count = 0
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                if count >= n_molecules:
                    break
                
                x = 0.5 + i * spacing
                y = 0.5 + j * spacing  
                z = 0.5 + k * spacing
                
                molecule = create_water_molecule(x, y, z)
                system = insert_molecule(system, molecule)
                count += 1
    
    return system


def calculate_system_energy(state):
    """Calculate total system energy, correcting for double counting"""
    total = 0.0
    for res in state.residues:
        if res.active:
            total += res.energy_vdw + res.energy_elec
    # Divide by 2 to correct for double counting if multiple residues
    return total / 2.0 if len(state.residues) > 0 else 0.0
