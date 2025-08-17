# tests/simulation/movementInsert/advanced_insertion_helpers.py
"""
Helper functions for advanced molecule insertion tests

This module provides energy calculation utilities.
"""

import pytest
import random
import math
import pygcmc

# Constants
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²
kB = 0.008314463  # Boltzmann constant in kJ/mol/K

def calculate_insertion_energy(original_system, new_system):
    """Calculate energy change for insertion using system energies"""
    # Calculate original system energy
    pygcmc.computeSystemEnergyCutoff(original_system)
    orig_energy = sum(res.energy_vdw + res.energy_elec for res in original_system.residues) / 2.0
    
    # Calculate new system energy
    pygcmc.computeSystemEnergyCutoff(new_system)
    new_energy = sum(res.energy_vdw + res.energy_elec for res in new_system.residues) / 2.0
    
    return new_energy - orig_energy

