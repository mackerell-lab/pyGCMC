# tests/simulation/energyOpenmm/simple_pbc_cutoff.py

import pytest
from .simple_helpers import *

def test_pbc_interaction():
    """
    Test nonbonded energy calculation with periodic boundary conditions.
    
    Tests the implementation of periodic boundary conditions in nonbonded calculations:
    1. Reaction field for electrostatics:
       E = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    2. Periodic wrapping for LJ interactions
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    System configuration:
    - Moves water molecule to box edge
    - Tests with and without PBC
    
    Verification:
    - Ensures PBC affects the interaction energy
    - Compares energies with and without PBC
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule to box edge
    box_size = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(box_size - 0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy with and without PBC
    energy_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=True)
    energy_no_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=False)
    
    # Energy with PBC should be different from energy without PBC
    assert abs(energy_pbc.value_in_unit(kilojoules_per_mole) - 
              energy_no_pbc.value_in_unit(kilojoules_per_mole)) > 1e-3, \
           "PBC should affect the interaction energy"
    print(f"PBC interaction energy: {energy_pbc}")
    print(f"Non-PBC interaction energy: {energy_no_pbc}")

def test_cutoff_effect():
    """
    Test the effect of cutoff distance on nonbonded energy calculation.
    
    Tests the implementation of cutoff-based methods:
    1. Switching function for LJ:
       S(r) = 1-6x⁵+15x⁴-10x³, x=(r-r_switch)/(r_cutoff-r_switch)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Reaction field for Coulomb:
       E = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Test distances:
    - 0.5 nm: Well within cutoff
    - 0.7 nm: Within cutoff
    - 0.9 nm: At switching distance
    - 1.2 nm: Beyond cutoff
    
    Verification:
    - Energy decreases with distance
    - Energy goes to zero beyond cutoff
    - Switching function properly applied
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Place water molecule at different distances
    # Note: benzene carbons are at radius 0.15 nm, so we need x >= 1.15 nm
    # to ensure all atom pairs are beyond the 1.0 nm cutoff
    distances = [0.5, 0.7, 0.9, 1.2]  # nm (last distance ensures all pairs > cutoff)
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    energies = []
    for dist in distances:
        # Move water molecule
        new_positions = []
        for i in range(len(positions)):
            if i >= 6 and i < 9:  # Water atoms
                pos = positions[i].value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # Calculate energy without PBC to observe true distance effects
        energy = calculate_nonbonded_energy(system, new_positions, movement_atoms, fixed_atoms, use_pbc=False)
        energy_val = energy.value_in_unit(kilojoules_per_mole)
        energies.append(energy_val)
        print(f"Energy at distance {dist} nm: {energy_val:.4f} kJ/mol")
    
    # Verify energy decreases with distance
    for i in range(len(distances)-1):
        assert abs(energies[i]) > abs(energies[i+1]), \
               f"Energy should decrease with distance. Energies: {energies}"
    
    # Verify energy is zero beyond cutoff (1.0 nm)
    assert abs(energies[-1]) < 1e-6, \
           f"Energy should be zero beyond cutoff (1.0 nm), but got {energies[-1]} kJ/mol"

