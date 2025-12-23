# tests/simulation/movementInsert/advanced_insertion_sequential_tests.py
"""
Sequential molecule insertion tests

This module tests sequential cavity filling and water cluster formation.
"""

import pytest
import random
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²
kB = 0.008314463  # Boltzmann constant in kJ/mol/K



from .advanced_insertion_helpers import calculate_insertion_energy
from .advanced_insertion_systems import *
from .advanced_insertion_molecules import *

def test_sequential_cavity_filling():
    """Test filling cavities sequentially based on favorability"""
    # Create system with ranked cavities
    system = create_ranked_cavity_system()
    
    # Define guest molecules to insert
    n_guests = 4
    filled_system = system
    insertion_energies = []
    
    for i in range(n_guests):
        # Find best insertion position by sampling
        best_pos = None
        best_energy = float('inf')
        
        # Sample positions around known cavity regions
        cavity_regions = [
            (1.5, 2.5, 2.5),  # Best cavity
            (3.5, 2.5, 2.5),  # Second best
            (2.5, 1.5, 2.5),  # Third
            (2.5, 3.5, 2.5),  # Fourth
        ]
        
        for base_x, base_y, base_z in cavity_regions:
            # Try positions near cavity center
            for dx in [-0.2, 0, 0.2]:
                for dy in [-0.2, 0, 0.2]:
                    x = base_x + dx
                    y = base_y + dy
                    z = base_z
                    
                    # Test insertion
                    test_guest = create_guest_molecule(x, y, z, charge=0.5)
                    test_system = insert_single_atom(filled_system, test_guest)
                    
                    # Calculate insertion energy
                    energy = calculate_insertion_energy(filled_system, test_system)
                    if energy < best_energy:
                        best_energy = energy
                        best_pos = (x, y, z)
        
        # Insert at best position
        if best_pos:
            x, y, z = best_pos
            guest = create_guest_molecule(x, y, z, charge=0.5)
            filled_system = insert_single_atom(filled_system, guest)
            insertion_energies.append(best_energy)
    
    # Debug output
    print(f"Insertion energies: {insertion_energies}")
    
    # Require non-zero energy signal (no silent pass).
    assert any(abs(e) > 1e-6 for e in insertion_energies), (
        "All insertion energies are ~0; expected non-zero values. "
        f"Energies={insertion_energies}"
    )
    
    # First insertions should be more favorable (or at least not worse)
    if len(insertion_energies) > 1:
        assert insertion_energies[0] <= insertion_energies[-1] + 0.1, \
            f"First cavity ({insertion_energies[0]:.4f}) should be no worse than last ({insertion_energies[-1]:.4f})"
    
    # Energy should generally become less favorable
    assert insertion_energies[0] < 0, "First insertion should be favorable"
    
    # Verify all molecules were inserted
    n_guests_inserted = sum(1 for res in filled_system.residues if res.type == 0)
    assert n_guests_inserted == n_guests, f"Should have inserted {n_guests} guests"


def test_water_cluster_formation():
    """Test water molecule insertion leading to cluster formation"""
    # Start with one water molecule
    system = create_empty_system()
    first_water = create_water_molecule(2.5, 2.5, 2.5)
    system = insert_molecule(system, first_water)
    
    # Insert additional water molecules nearby
    n_waters = 5
    cluster_positions = [
        (2.8, 2.5, 2.5),  # Hydrogen bonding distance
        (2.2, 2.5, 2.5),  # Other side
        (2.5, 2.8, 2.5),  # Above
        (2.5, 2.2, 2.5),  # Below
        (2.5, 2.5, 2.8),  # Front
    ]
    
    total_energies = []
    
    # Calculate initial energy
    pygcmc.computeSystemEnergyCutoff(system)
    total_energies.append(calculate_system_energy(system))
    
    for x, y, z in cluster_positions:
        water = create_water_molecule(x, y, z)
        system = insert_molecule(system, water)
        
        # Calculate new total energy
        pygcmc.computeSystemEnergyCutoff(system)
        total_energies.append(calculate_system_energy(system))
    
    # Debug output
    print(f"Water cluster formation energies: {[f'{e:.4f}' for e in total_energies]}")
    
    # Require non-zero energy signal (no silent pass).
    assert any(abs(e) > 1e-6 for e in total_energies), (
        "All water cluster energies are ~0; expected non-zero values. "
        f"Energies={total_energies}"
    )
    
    # Adding waters should change the total energy magnitude (direction can vary by FF).
    abs_values = [abs(e) for e in total_energies]
    assert max(abs_values) > 1.0, (
        "Cluster formation should produce a noticeable energy change. "
        f"Energies={total_energies}"
    )
    rounded = {round(e, 6) for e in total_energies}
    assert len(rounded) > 1, "Total energy should vary across insertions"
