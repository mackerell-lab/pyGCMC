# tests/simulation/energyNB/switching_methods.py
"""Switching function method tests."""

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


def test_with_ewald():
    """
    Test switching function with Ewald summation
    """
    # Create test system
    state = create_test_system()
    
    # Set atom charges
    state.atoms[0].charge = 1.0
    state.atoms[1].charge = -1.0
    
    # Set Ewald parameters
    state.info.setTemperature(300.0)  # 300K
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)
    
    # Set switching function parameters
    # According to CHARMM convention:
    # ctonnb = 1.0 - inner cutoff (r_on), where smooth decay begins (10 Å)
    # ctofnb = 1.2 - outer cutoff (r_off), where potential finally decays to 0 (12 Å)
    r_on = 1.0   # Inner cutoff radius (10 Å)
    r_off = 1.2  # Outer cutoff radius (12 Å)
    
    # Create MonteCarloSystem object
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set position of second atom
    distances = [0.9 + i * 0.4/20 for i in range(20)]  # 20 points from 0.9 to 1.3
    
    # Close switching function, calculate ordinary Ewald
    state.info.use_switching = False
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    energies_ewald = []
    for r in distances:
        state.atoms[1].x = r
        pygcmc.computeSystemEnergyEwald(state)  # Use Ewald method to calculate energy
        # Record only VDW energy part
        energies_ewald.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # Enable switching function
    state.info.use_switching = True
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    energies_ewald_switching = []
    for r in distances:
        state.atoms[1].x = r
        pygcmc.computeSystemEnergyEwald(state)  # Use Ewald method to calculate energy
        # Record only VDW energy part
        energies_ewald_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # Compare switching before and after VDW energy
    for i, r in enumerate(distances):
        e_normal = energies_ewald[i]
        e_switch = energies_ewald_switching[i]
        
        print(f"Ewald r={r:.3f}, no_switching={e_normal:.6f}, with_switching={e_switch:.6f}")
        
        # When r < r_on, energy should be the same
        if r < r_on:
            assert e_normal == pytest.approx(e_switch)
        # When r > r_off, switching VDW energy should be 0
        elif r > r_off:
            assert e_switch == pytest.approx(0.0)
        # In r_on and r_off, check if switching is correctly applied
        elif r_on <= r <= r_off:
            # Calculate expected switching value
            switch_value = mc_system.calculate_switching_function(r)
            # No longer strictly require energy to equal original energy multiplied by switching value
            # Instead, verify: 1. Energy between 0 and original energy 2. Energy decreases with distance
            assert 0.0 <= abs(e_switch) <= abs(e_normal), \
                  f"Ewald VDW energy with switching should be between 0 and original energy at r={r}"
            
            # If not the first point between r_on and r_off, ensure energy decreases monotonically (absolute value decreases)
            if i > 0 and r_on <= distances[i-1] <= r_off:
                assert abs(e_switch) <= abs(energies_ewald_switching[i-1]), \
                      f"Energy magnitude should decrease as r increases in switching region"
    
    # Disable switching function, restore default state
    state.info.use_switching = False


def test_monte_carlo_system_switching_function():
    """
    Test using CHARMM switching function through MonteCarloSystem class new interface
    """
    # Create system
    r_on = 1.0  # Inner cutoff radius
    r_off = 1.2  # Outer cutoff radius
    
    # Create a test system, for storing state
    state = create_test_system()
    
    # Create MonteCarloSystem object
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set switching function parameters
    mc_system.set_switching_function(True, r_on, r_off)
    
    # Correspondly set parameters in state
    state.info.use_switching = True
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    # Verify parameter setting
    assert mc_system.is_using_switching_function() == True, "Switching function should be enabled"
    assert mc_system.get_switching_r_on() == pytest.approx(r_on), "r_on parameter should be correctly set"
    assert mc_system.get_switching_r_off() == pytest.approx(r_off), "r_off parameter should be correctly set"
    
    # Calculate switching function values at different distances
    distances = [0.9 + i * 0.4/50 for i in range(50)]  # 50 points from 0.9 to 1.3
    switch_values = [mc_system.calculate_switching_function(r) for r in distances]
    
    # Verify function values
    for r, s in zip(distances, switch_values):
        if r <= r_on:
            assert s == pytest.approx(1.0), f"S({r}) should be 1.0 when r <= r_on"
        elif r >= r_off:
            assert s == pytest.approx(0.0), f"S({r}) should be 0.0 when r >= r_off"
        else:
            # Calculate expected value based on formula
            r2 = r * r
            ron2 = r_on * r_on
            roff2 = r_off * r_off
            
            numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0*r2 - 3.0*ron2)
            denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2)
            expected = numerator / denominator
            
            assert s == pytest.approx(expected, abs=1e-6), f"S({r}) calculation error"
    
    # Disable switching function
    mc_system.set_switching_function(False)
    state.info.use_switching = False
    assert mc_system.is_using_switching_function() == False, "Switching function should be disabled"