# tests/simulation/movementInsert/practical_gcmc_water_test.py
"""
Practical GCMC water insertion test

This module tests GCMC insertion of water molecules with real energy calculations.
"""

import pytest
import math
import random
import pygcmc
import os

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²


from .practical_gcmc_systems import *
from .practical_gcmc_molecules import *

def test_water_insertion_with_real_energy():
    """Test water insertion with real PyGCMC energy calculations"""
    
    # GCMC parameters for water
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -6.3  # kJ/mol for water at 300K
    n_bar = 55.5  # mol/L water concentration
    B = beta * mu_ex + math.log(n_bar)
    
    print("\n=== Water Insertion with Real Energy ===")
    
    # Create initial water box
    system = create_water_box_system(n_waters=10, box_size=2.5)
    
    # Calculate initial energy
    pygcmc.computeSystemEnergyCutoff(system)
    initial_energy = calculate_total_energy(system)
    initial_n = count_water_residues(system)
    print(f"Initial energy ({initial_n} waters): {initial_energy:.2f} kJ/mol")
    
    # Perform insertions
    n_attempts = 10
    accepted = 0
    
    for attempt in range(n_attempts):
        # Find suitable position
        pos = find_water_insertion_site(system)
        if pos is None:
            print(f"Attempt {attempt+1}: No suitable site found")
            continue
        
        x, y, z = pos
        
        # Create water molecule
        water_atoms = create_water_molecule(x, y, z)
        
        # Calculate insertion energy
        trial_system = create_system_with_guest(system, water_atoms)
        pygcmc.computeSystemEnergyCutoff(trial_system)
        trial_energy = calculate_total_energy(trial_system)
        delta_e = trial_energy - initial_energy
        
        # Cavity bias (based on minimum distance)
        min_dist = calculate_minimum_distance(system, x, y, z)
        f_n = 0.8 if min_dist > 0.35 else 0.2  # Simple cavity factor
        
        # GCMC acceptance
        n_current = count_water_residues(system)
        acc_arg = (f_n / (n_current + 1)) * math.exp(B - beta * delta_e)
        acc_prob = min(1.0, acc_arg)
        
        if random.random() < acc_prob:
            accepted += 1
            system = trial_system
            initial_energy = trial_energy
            print(f"Attempt {attempt+1}: ACCEPTED - ΔE = {delta_e:6.2f} kJ/mol, "
                  f"n = {n_current+1}")
        else:
            print(f"Attempt {attempt+1}: REJECTED - ΔE = {delta_e:6.2f} kJ/mol, "
                  f"P = {acc_prob:.4f}")
    
    acceptance_rate = accepted / n_attempts
    print(f"\nWater acceptance rate: {acceptance_rate:.2%}")
    final_n = count_water_residues(system)
    print(f"Final number of waters: {final_n}")
    
    # Should have reasonable acceptance
    assert accepted > 0, "No water insertions accepted"
    # With our test parameters, we might not add many waters
    assert final_n >= initial_n, f"Should not lose waters (started with {initial_n}, ended with {final_n})"


