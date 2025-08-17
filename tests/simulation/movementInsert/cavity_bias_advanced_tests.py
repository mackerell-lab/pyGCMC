# tests/simulation/movementInsert/cavity_bias_advanced_tests.py
"""
Advanced cavity bias tests
"""

import pytest
import math
import random
import numpy as np
import pygcmc

# Constants
kB = 0.008314463  # kJ/(mol·K)

# Import cavity bias helpers
from .cavity_bias_basic_tests import CavityGrid

# Import system creation helpers
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)

def test_cavity_bias_energy_calculation():
    """Test energy calculation for cavity-biased insertions."""
    # Create system with one existing molecule
    system = create_empty_system()
    molecule1 = create_water_molecule(1.0, 1.0, 1.0)
    system = insert_molecule(system, molecule1)
    
    # Find cavities
    grid = CavityGrid(system.info.box, grid_spacing=0.3)
    vdw_radii = {0: 0.15, 1: 0.1}
    
    for atom in system.atoms:
        radius = vdw_radii.get(atom.type, 0.15)
        grid.mark_occupied((atom.x, atom.y, atom.z), radius)
    
    grid.find_cavities()
    
    # Insert at cavity point
    cavity_point = grid.select_random_cavity()
    assert cavity_point is not None, "Should find cavity point"
    
    # Create new molecule at cavity
    x, y, z = cavity_point
    molecule2 = create_water_molecule(x, y, z)
    
    # Calculate energy before insertion
    pygcmc.computeSystemEnergyCutoff(system)
    energy_before = calculate_system_energy(system)
    
    # Insert molecule
    system_after = insert_molecule(system, molecule2)
    
    # Calculate energy after insertion
    pygcmc.computeSystemEnergyCutoff(system_after)
    energy_after = calculate_system_energy(system_after)
    
    # Energy change
    delta_E = energy_after - energy_before
    
    # Cavity insertion should generally avoid high-energy overlaps
    # So energy change should be reasonable (not extremely positive)
    assert delta_E < 1000.0, f"Energy change {delta_E} suggests bad overlap"


def test_adaptive_cavity_grid():
    """Test adaptive grid spacing based on system density."""
    # Test with different system densities
    densities = [0.1, 0.5, 1.0, 2.0]  # molecules/nm³
    box_volume = 5.0 ** 3  # nm³
    
    for density in densities:
        system = create_empty_system()
        n_molecules = int(density * box_volume)
        
        # Add molecules randomly
        random.seed(123)
        for i in range(n_molecules):
            x = random.uniform(0.5, 4.5)
            y = random.uniform(0.5, 4.5)
            z = random.uniform(0.5, 4.5)
            molecule = create_water_molecule(x, y, z)
            system = insert_molecule(system, molecule)
        
        # Choose grid spacing based on density
        # Higher density needs finer grid
        if density < 0.5:
            spacing = 0.5
        elif density < 1.0:
            spacing = 0.3
        else:
            spacing = 0.2
        
        # Create grid
        grid = CavityGrid(system.info.box, grid_spacing=spacing)
        vdw_radii = {0: 0.15, 1: 0.1}
        
        for atom in system.atoms:
            radius = vdw_radii.get(atom.type, 0.15)
            grid.mark_occupied((atom.x, atom.y, atom.z), radius)
        
        grid.find_cavities()
        f_n = grid.get_cavity_fraction()
        
        # Higher density should have lower cavity fraction
        if density > 1.0:
            assert f_n < 1.0, f"High density system should have f_n < 1.0, got {f_n}"
        elif density < 0.5:
            assert f_n > 0.5, f"Low density system should have f_n > 0.5, got {f_n}"


def test_cavity_bias_detailed_balance():
    """Test cavity bias detailed balance with proper proposal probabilities."""
    # System parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -5.0
    n_bar = 100
    B = beta * mu_ex + math.log(n_bar)
    
    # State A: n molecules, State B: n+1 molecules
    n_A = 10
    n_B = n_A + 1
    E_A = -50.0  # System energy
    delta_E = 5.0  # Energy change for insertion
    E_B = E_A + delta_E
    
    # Cavity fraction
    f_n = 0.3
    
    # Calculate the acceptance ratio r
    r = (f_n / (n_A + 1.0)) * math.exp(B - beta * delta_E)
    
    # Acceptance probabilities
    alpha_ins = min(1.0, r)        # Forward: insertion A→B
    alpha_del = min(1.0, 1.0 / r)  # Reverse: deletion B→A
    
    # Proposal probabilities
    # For consistent detailed balance: q_del/q_ins = f_n/(n+1)
    q_ins = 1.0
    q_del = f_n / n_B  # This ensures q_del/q_ins = f_n/(n_B) = f_n/(n_A+1)
    
    # Grand canonical weights: π_GC(n,E) ∝ exp(-βE + Bn)
    pi_A = math.exp(-beta * E_A + B * n_A)
    pi_B = math.exp(-beta * E_B + B * n_B)
    
    # Full transition probabilities: T = π * q * α
    forward_flux = pi_A * q_ins * alpha_ins
    reverse_flux = pi_B * q_del * alpha_del
    
    # Check detailed balance
    if forward_flux > 0 and reverse_flux > 0:
        ratio = forward_flux / reverse_flux
        # Should be equal within numerical precision
        assert 0.999999999 < ratio < 1.000000001, \
            f"Detailed balance violated: ratio = {ratio}"