# tests/simulation/movementInsert/advanced_insertion_biased_tests.py
"""
Biased molecule insertion tests

This module tests biased cavity insertion and energy-guided insertion strategies.
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

def test_biased_cavity_insertion():
    """Test biased insertion into detected cavities"""
    # Create system with multiple cavities
    system = create_multi_cavity_system()
    
    # Define cavity centers
    cavity_centers = [
        (1.5, 1.5, 2.5),  # Small cavity
        (3.5, 3.5, 2.5),  # Large cavity
        (2.5, 2.5, 1.0),  # Edge cavity
    ]
    
    # Test guest molecule insertion in each cavity
    cavity_energies = []
    for center in cavity_centers:
        x, y, z = center
        guest = create_guest_molecule(x, y, z, charge=0.5)
        test_system = insert_single_atom(system, guest)
        
        # Calculate insertion energy
        insertion_energy = calculate_insertion_energy(system, test_system)
        cavity_energies.append(insertion_energy)
    
    # Debug: print cavity energies
    print(f"Cavity energies: Small={cavity_energies[0]:.4f}, Large={cavity_energies[1]:.4f}, Edge={cavity_energies[2]:.4f}")
    
    # Require energy discrimination across cavities (no silent pass).
    spread = max(cavity_energies) - min(cavity_energies)
    assert spread > 1e-6, (
        "All cavities have the same energy; expected discrimination. "
        f"Energies={cavity_energies}"
    )
    
    # Large cavity should be most favorable
    assert cavity_energies[1] < cavity_energies[0], \
        f"Large cavity ({cavity_energies[1]:.4f}) should be more favorable than small cavity ({cavity_energies[0]:.4f})"
    
    # Edge cavity might have different energy due to boundary
    assert len(set(cavity_energies)) > 1, \
        "Different cavities should have different energies"


def test_energy_guided_insertion():
    """Test insertion guided by energy calculations"""
    # Create a system with existing molecules
    system = create_mixed_system()
    
    # Try multiple insertion attempts
    n_attempts = 20
    random.seed(42)
    
    insertion_attempts = []
    for i in range(n_attempts):
        # Generate random position
        x = random.uniform(0.5, 4.5)
        y = random.uniform(0.5, 4.5)
        z = random.uniform(0.5, 4.5)
        
        # Create test molecule
        test_molecule = create_guest_molecule(x, y, z, charge=1.0)
        test_system = insert_single_atom(system, test_molecule)
        
        # Calculate insertion energy
        insertion_energy = calculate_insertion_energy(system, test_system)
        
        insertion_attempts.append({
            'position': (x, y, z),
            'energy': insertion_energy
        })
    
    # Sort by energy
    insertion_attempts.sort(key=lambda x: x['energy'])
    
    # Best positions should have negative energy (favorable)
    best_energy = insertion_attempts[0]['energy']
    worst_energy = insertion_attempts[-1]['energy']
    
    print(f"Energy range: best={best_energy:.4f}, worst={worst_energy:.4f}")
    
    # Require non-zero energy signal to avoid silent pass.
    assert any(abs(att["energy"]) > 1e-6 for att in insertion_attempts), (
        "All insertion energies are ~0; expected non-zero values. "
        f"Energies={[att['energy'] for att in insertion_attempts]}"
    )
    
    # If we have non-zero energies, check the pattern
    assert best_energy <= 0, f"Best insertion should be favorable or neutral, got {best_energy:.2f}"
    assert worst_energy >= best_energy, "Energy range should span from favorable to unfavorable"
    
    # Calculate acceptance probabilities at 300K
    T = 300.0
    beta = 1.0 / (kB * T)
    
    # Best position should have high acceptance
    p_best = min(1.0, math.exp(-beta * best_energy))
    assert p_best > 0.9, f"Best position should have high acceptance, got {p_best:.3f}"
    
    # Worst position should have low acceptance
    p_worst = min(1.0, math.exp(-beta * worst_energy))
    assert p_worst < 0.1, f"Worst position should have low acceptance, got {p_worst:.3f}"

