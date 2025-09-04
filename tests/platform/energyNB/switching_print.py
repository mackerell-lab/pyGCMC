# tests/simulation/energyNB/switching_print.py
"""Switching function print and display tests."""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCInfo, MCMovementResidueInfo
import os
import math
import sys

# Set log level to INFO to view debug output
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)


def create_test_system(box_size=4.0):
    """
    Create a simple test system with two atoms
    
    Args:
        box_size: box size (nm)
    Returns:
        state: MC state object with two atoms
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff (12 Å) - Using CHARMM's ctofnb value as cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Type 0 and Type 1
    
    # LJ parameters
    sigma = 0.3  # nm
    eps = 0.5    # kJ/mol
    
    # Set LJ parameter matrix (2x2 matrix flattened to array)
    ff.ljSigma = [
        sigma, sigma,  # type 0 to types 0, 1
        sigma, sigma   # type 1 to types 0, 1
    ]
    
    ff.ljEps = [
        eps, eps,  # type 0 to types 0, 1
        eps, eps   # type 1 to types 0, 1
    ]
    
    state.forcefield = ff
    
    # First atom
    atom1 = MCAtom()
    atom1.type = 0
    atom1.charge = 0.0  # Neutral
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    
    # Second atom
    atom2 = MCAtom()
    atom2.type = 1
    atom2.charge = 0.0  # Neutral
    atom2.x = 1.0  # Will be adjusted in tests
    atom2.y = 0.0
    atom2.z = 0.0
    
    # Add atoms to state
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # Create two residues and atoms
    res1 = MCResidue()
    res1.active = True
    res1.atomStart = 0
    res1.atomCount = 1
    res1.type = 0
    
    res2 = MCResidue()
    res2.active = True
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 1
    
    # Add residues to state
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    return state


def test_print_switching_values():
    """
    Print switching function values and energy data for debugging and documentation purposes
    """
    # According to CHARMM convention:
    # ctonnb = 1.0 - inner cutoff, where smooth decay begins (10 Å)
    # ctofnb = 1.2 - outer cutoff, where potential finally decays to 0 (12 Å)
    r_on = 1.0  # 10 Å
    r_off = 1.2  # 12 Å
    
    # Create MonteCarloSystem object
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set switching function parameters - Only set mc_system parameters for calculate_switching_function call
    mc_system.set_switching_function(True, r_on, r_off)
    
    # Calculate switching function values
    key_distances = [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
    switch_values = [mc_system.calculate_switching_function(r) for r in key_distances]
    
    # Print header
    print("\nCHARMM Switching Function Values:")
    print(f"{'Distance (nm)':15s} | {'Switch Value':15s}")
    print("-" * 33)
    
    # Print switching function values
    for r, s in zip(key_distances, switch_values):
        print(f"{r:15.3f} | {s:15.4f}")
    
    # Disable switching function
    mc_system.set_switching_function(False)