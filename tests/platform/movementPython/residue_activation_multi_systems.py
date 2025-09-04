# tests/simulation/movementInsert/residue_activation_multi_systems.py
"""
Helper functions for creating multi-type test systems
"""

import pytest
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²


def create_multi_type_system():
    """Create a system with multiple residue types"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    type0 = state.atomTypes.get_or_add_type("TYPE0")
    type1 = state.atomTypes.get_or_add_type("TYPE1")
    
    # Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1  # Only type 0 moves
    state.forcefield.ljEps = [1.0, 0.8, 0.8, 0.6]
    state.forcefield.ljSigma = [0.3, 0.35, 0.35, 0.4]
    
    # Create atoms
    atoms = []
    
    # Type 0 atoms (movement)
    for i in range(2):
        atom = pygcmc.MCAtom()
        atom.x = 1.0 + i * 0.5
        atom.y = 1.0
        atom.z = 1.0
        atom.charge = 0.5 if i == 0 else -0.5
        atom.type = type0
        atoms.append(atom)
    
    # Type 1 atom (fixed)
    atom = pygcmc.MCAtom()
    atom.x = 2.5
    atom.y = 1.0
    atom.z = 1.0
    atom.charge = 0.0
    atom.type = type1
    atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # Create residues
    residues = []
    
    # Type 0 residues
    for i in range(2):
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Type 1 residue
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = 2
    res.atomCount = 1
    res.type = 1
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    return state


def calculate_total_energy(state):
    """Calculate total system energy, correcting for double counting"""
    total = 0.0
    for res in state.residues:
        if res.active:
            energy = res.energy_vdw + res.energy_elec
            total += energy
    # Divide by 2 to correct for double counting in system energy calculation
    return total / 2.0
