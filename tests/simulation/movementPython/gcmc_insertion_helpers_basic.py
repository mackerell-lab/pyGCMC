# tests/simulation/movementInsert/gcmc_insertion_helpers.py
"""Helper functions for GCMC insertion tests"""

import pytest, math, random, os
import pygcmc

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²

# Get data directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(TEST_DIR))), "data")

def create_protein_system(pdb_file):
    """Create MCState from protein PDB file"""
    # This is a placeholder - would need actual PDB parsing
    state = pygcmc.MCState()
    state.info.box = [13.2126, 13.2213, 12.2070]  # From 4wp7 CRYST1
    state.info.cutoff = 1.2
    
    # Set up basic forcefield
    state.forcefield.numTotalTypes = 20  # Approximate for protein
    state.forcefield.numMovementTypes = 1  # Only guest molecules move
    
    # Placeholder for protein atoms (would be loaded from PDB)
    # For now, just create a few atoms to represent protein
    for i in range(100):
        atom = pygcmc.MCAtom()
        atom.x = random.uniform(2.0, 11.0)
        atom.y = random.uniform(2.0, 11.0)
        atom.z = random.uniform(2.0, 10.0)
        atom.charge = random.choice([-0.5, 0.0, 0.5])
        atom.type = random.randint(1, 19)  # Protein atom types
    
    return state

def load_benzene_molecule(pdb_file):
    """Load benzene molecule structure"""
    # Benzene coordinates from PDB (simplified)
    atoms = []
    coords = [
        (-0.013, -0.023, 0.000, -0.115),  # CG
        (1.388, -0.023, 0.000, -0.115),   # CD1
        (-0.714, 1.191, 0.000, -0.115),   # CD2
        (2.089, 1.191, 0.000, -0.115),    # CE1
        (-0.013, 2.405, 0.000, -0.115),   # CE2
        (1.388, 2.405, 0.000, -0.115),    # CZ
    ]
    
    for x, y, z, charge in coords:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charge
        atom.type = 0  # Guest molecule type
        atoms.append(atom)
    
    # Add hydrogens (charge = 0.115)
    for i in range(6):
        atom = pygcmc.MCAtom()
        atom.charge = 0.115
        atom.type = 0
        atoms.append(atom)
    
    return atoms

def find_cavity_sites(system, probe_radius=0.3, grid_spacing=0.05):
    """Find cavity sites in protein using grid search"""
    cavities = []
    
    # Simple grid search (would use more sophisticated method in practice)
    nx = int(system.info.box[0] / grid_spacing)
    ny = int(system.info.box[1] / grid_spacing)
    nz = int(system.info.box[2] / grid_spacing)
    
    # Sample a subset of grid points for efficiency
    n_samples = min(1000, nx * ny * nz // 100)
    
    for _ in range(n_samples):
        x = random.uniform(probe_radius, system.info.box[0] - probe_radius)
        y = random.uniform(probe_radius, system.info.box[1] - probe_radius)
        z = random.uniform(probe_radius, system.info.box[2] - probe_radius)
        
        # Check if point is in cavity (simplified)
        min_dist = float('inf')
        for atom in system.atoms:
            dist = math.sqrt((x - atom.x)**2 + (y - atom.y)**2 + (z - atom.z)**2)
            min_dist = min(min_dist, dist)
        
        if min_dist > probe_radius * 2:  # Cavity criterion
            cavities.append((x, y, z))
    
    return cavities


def insert_benzene_at_position(system, benzene, x, y, z, theta, phi, psi):
    """Insert benzene molecule at specified position and orientation"""
    # This would create a new system with benzene inserted
    # For now, return the original system
    return system


def calculate_insertion_energy(trial_system, original_system):
    """Calculate energy difference for insertion"""
    # Simplified - would calculate actual VDW and electrostatic energies
    return random.uniform(-20.0, 50.0)  # kJ/mol


def count_benzene_molecules(system):
    """Count number of benzene molecules in system"""
    # Simplified - would count actual benzene residues
    return 0


def create_water_box(n_waters, box_size):
    """Create a box with water molecules"""
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.2
    
    # Water forcefield parameters
    state.forcefield.numTotalTypes = 2  # O and H
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [0.6364, 0.0, 0.0, 0.0]  # TIP3P-like
    state.forcefield.ljSigma = [0.3166, 0.0, 0.0, 0.0]
    
    # Add water molecules
    state.atoms = []
    state.residues = []
    
    # Place waters in a grid
    n_per_dim = int(math.ceil(n_waters ** (1/3)))
    spacing = box_size / (n_per_dim + 1)
    
    water_count = 0
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_count >= n_waters:
                    break
                
                x = (i + 1) * spacing
                y = (j + 1) * spacing
                z = (k + 1) * spacing
                
                # Add water molecule
                atom_start = len(state.atoms)
                
