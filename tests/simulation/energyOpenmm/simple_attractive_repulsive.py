# tests/simulation/energyOpenmm/simple_attractive_repulsive.py

import pytest
from .simple_helpers import *

def test_attractive_interaction():
    """
    Test nonbonded energy calculation for an attractive interaction.
    
    Tests the attractive regime of the Lennard-Jones and Coulomb potentials:
    1. LJ attraction: r > r_min where r_min = 2^(1/6)σ
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb attraction: opposite charges
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    System configuration:
    - Benzene (+0.1e per C) interacting with water (-0.834e on O, +0.417e on H)
    - Default separation should result in net attractive force
    
    Verification:
    - Ensures total energy is negative (attractive)
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be negative (attractive) due to opposite charges
    assert energy.value_in_unit(kilojoules_per_mole) < 0, \
           f"Expected attractive interaction, got {energy}"
    print(f"Attractive interaction energy: {energy}")

def test_repulsive_interaction():
    """
    Test nonbonded energy calculation for a repulsive interaction.
    
    Tests the repulsive regime of the Lennard-Jones potential:
    1. LJ repulsion: r < r_min where r_min = 2^(1/6)σ
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    System configuration:
    - Moves water molecule very close to benzene (0.1 nm)
    - At this distance, LJ repulsion dominates over electrostatic attraction
    
    Verification:
    - Ensures total energy is positive (repulsive)
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule very close to benzene to create repulsion
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be positive (repulsive) due to close distance
    assert energy.value_in_unit(kilojoules_per_mole) > 0, \
           f"Expected repulsive interaction, got {energy}"
    print(f"Repulsive interaction energy: {energy}")

