# tests/simulation/energyOpenmm/intermediate_energy_symmetry.py

import pytest
from .intermediate_helpers import *

def test_energy_symmetry():
    """
    Test that energy calculation is symmetric (A->B equals B->A).
    
    Tests the fundamental physical principle that nonbonded interactions are symmetric:
    1. LJ interaction symmetry:
       E_LJ(1,2) = E_LJ(2,1)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb interaction symmetry:
       E_coul(1,2) = E_coul(2,1)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    Verification:
    - Calculates energy both ways (group1->group2 and group2->group1)
    - Ensures absolute difference is negligible
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define groups
    group1 = set(range(6))  # Benzene
    group2 = set(range(6, 9))  # Water
    
    # Calculate energy both ways
    energy1 = calculate_nonbonded_energy(system, positions, group1, group2)
    energy2 = calculate_nonbonded_energy(system, positions, group2, group1)
    
    # Energies should be equal
    assert abs(energy1.value_in_unit(kilojoules_per_mole) - 
              energy2.value_in_unit(kilojoules_per_mole)) < 1e-6, \
           f"Energy calculation should be symmetric: {energy1} != {energy2}"
    print(f"Forward energy: {energy1}")
    print(f"Reverse energy: {energy2}")

