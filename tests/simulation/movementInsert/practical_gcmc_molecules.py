# tests/simulation/movementInsert/practical_gcmc_molecules.py
"""
Molecule creation and utility functions for practical GCMC tests
"""

import pytest
import math
import random
import pygcmc
import os

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²

def create_water_box_system(n_waters, box_size):
    """Create water box system"""
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.0
    
    # Water parameters
    o_type = state.atomTypes.get_or_add_type("OT")
    h_type = state.atomTypes.get_or_add_type("HT")
    
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [0.6364, 0.0, 0.0, 0.0]
    state.forcefield.ljSigma = [0.3166, 0.0, 0.0, 0.0]
    
    state.atoms = []
    state.residues = []
    
    # Add waters in grid
    spacing = box_size / (int(n_waters**(1/3)) + 1)
    water_count = 0
    
    for i in range(int(n_waters**(1/3)) + 1):
        for j in range(int(n_waters**(1/3)) + 1):
            for k in range(int(n_waters**(1/3)) + 1):
                if water_count >= n_waters:
                    break
                
                x = spacing * (i + 0.5)
                y = spacing * (j + 0.5)
                z = spacing * (k + 0.5)
                
                # Add water
                atom_start = len(state.atoms)
                
                # Oxygen
                o_atom = pygcmc.MCAtom()
                o_atom.x = x
                o_atom.y = y
                o_atom.z = z
                o_atom.charge = -0.834
                o_atom.type = o_type
                state.atoms.append(o_atom)
                
                # Hydrogens
                for h in range(2):
                    h_atom = pygcmc.MCAtom()
                    h_atom.x = x + 0.1 * (h - 0.5)
                    h_atom.y = y + 0.05
                    h_atom.z = z
                    h_atom.charge = 0.417
                    h_atom.type = h_type
                    state.atoms.append(h_atom)
                
                # Residue
                res = pygcmc.MCResidue()
                res.active = True
                res.atomStart = atom_start
                res.atomCount = 3
                res.type = 0
                state.residues.append(res)
                
                water_count += 1
    
    state.atoms = state.atoms
    state.residues = state.residues
    state.activeAtomCount = len(state.atoms)
    state.activeResidueCount = len(state.residues)
    
    return state


def find_water_insertion_site(system):
    """Find suitable site for water insertion"""
    for _ in range(20):
        x = random.uniform(0.3, system.info.box[0] - 0.3)
        y = random.uniform(0.3, system.info.box[1] - 0.3)
        z = random.uniform(0.3, system.info.box[2] - 0.3)
        
        min_dist = calculate_minimum_distance(system, x, y, z)
        if min_dist > 0.25:  # 2.5 Å minimum
            return (x, y, z)
    
    return None


def create_water_molecule(x, y, z):
    """Create water molecule at position"""
    o_atom = pygcmc.MCAtom()
    o_atom.x = x
    o_atom.y = y
    o_atom.z = z
    o_atom.charge = -0.834
    o_atom.type = 0  # OT
    
    h1_atom = pygcmc.MCAtom()
    h1_atom.x = x + 0.0957
    h1_atom.y = y
    h1_atom.z = z
    h1_atom.charge = 0.417
    h1_atom.type = 1  # HT
    
    h2_atom = pygcmc.MCAtom()
    h2_atom.x = x - 0.0757
    h2_atom.y = y + 0.0587
    h2_atom.z = z
    h2_atom.charge = 0.417
    h2_atom.type = 1  # HT
    
    return [o_atom, h1_atom, h2_atom]


def calculate_minimum_distance(system, x, y, z):
    """Calculate minimum distance to existing atoms"""
    min_dist = float('inf')
    for atom in system.atoms:
        dist = math.sqrt((x - atom.x)**2 + (y - atom.y)**2 + (z - atom.z)**2)
        min_dist = min(min_dist, dist)
    return min_dist


def count_water_residues(system):
    """Count water residues"""
    return sum(1 for res in system.residues if res.atomCount == 3)


def create_system_with_density(density):
    """Create system with specified occupancy density"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Add atoms to achieve density
    n_atoms = int(1000 * density)  # Total possible sites = 1000
    
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [0.3]
    
    state.atoms = []
    for _ in range(n_atoms):
        atom = pygcmc.MCAtom()
        atom.x = random.uniform(0.5, 4.5)
        atom.y = random.uniform(0.5, 4.5)
        atom.z = random.uniform(0.5, 4.5)
        atom.charge = 0.0
        atom.type = 0
        state.atoms.append(atom)
    
    state.atoms = state.atoms
    state.activeAtomCount = n_atoms
    
    return state


if __name__ == "__main__":
    # Run tests
    print("=== Practical GCMC Insertion Tests ===\n")
    
    print("1. Testing benzene insertion...")
    test_practical_benzene_insertion()
    
    print("\n2. Testing water insertion...")
    test_water_insertion_with_real_energy()
    
    print("\n3. Testing cavity bias effect...")
