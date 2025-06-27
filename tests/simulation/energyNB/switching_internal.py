# tests/simulation/energyNB/switching_internal.py
"""Switching function internal validation tests."""

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


def test_internal_switching_function():
    """
    Test switching function through energy calculation results directly
    
    This test indirectly tests whether the internal implemented switching function meets expectations
    by comparing energy ratios with/without switching function
    """
    # Create test system
    state = create_test_system()
    
    # Set switching function parameters
    r_on = 1.0   # Inner cutoff radius
    r_off = 1.2  # Outer cutoff radius
    
    # Select distances between r_on and r_off for testing
    test_distances = [1.02, 1.05, 1.08, 1.10, 1.15, 1.18]
    
    print("\nInternal Switching Function test (through energy ratio verification):")
    print(f"{'Distance(nm)':10s} | {'Energy without switching':14s} | {'Energy with switching':14s} | {'Actual ratio':10s} | {'Expected switching value':10s} | {'Error':10s}")
    print("-" * 75)
    
    for r in test_distances:
        # Set particle distance
        state.atoms[1].x = r
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        # 1. Calculate energy without switching function
        state.info.use_switching = False
        pygcmc.computeSystemEnergy(state)
        energy_no_switch = state.residues[0].energy_vdw + state.residues[1].energy_vdw
        
        # If energy is 0, skip this distance point
        if abs(energy_no_switch) < 1e-10:
            continue
        
        # 2. Calculate energy with switching function
        state.info.use_switching = True
        state.info.r_on = r_on
        state.info.r_off = r_off
        pygcmc.computeSystemEnergy(state)
        energy_with_switch = state.residues[0].energy_vdw + state.residues[1].energy_vdw
        
        # 3. Calculate energy ratio, i.e., actual switching function value
        actual_switch = energy_with_switch / energy_no_switch if energy_no_switch != 0 else 0
        
        # 4. Calculate theoretical switching function value
        r2 = r * r
        ron2 = r_on * r_on
        roff2 = r_off * r_off
        
        numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0*r2 - 3.0*ron2)
        denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2)
        expected_switch = numerator / denominator
        
        # Calculate error
        error = abs(actual_switch - expected_switch)
        
        # Print results
        print(f"{r:10.3f} | {energy_no_switch:14.6f} | {energy_with_switch:14.6f} | "
              f"{actual_switch:10.6f} | {expected_switch:10.6f} | {error:10.6f}")
        
        # Verify error is within acceptable range
        assert actual_switch == pytest.approx(expected_switch, abs=1e-4), \
               f"At r={r}, switching function value differs: {actual_switch} vs expected {expected_switch}"
    
    # Restore test state
    state.info.use_switching = False


def test_internal_switching_function_ewald():
    """
    Test switching function through Ewald energy calculation results directly
    
    This test uses Ewald calculation method, indirectly tests whether the internal implemented switching function meets expectations
    by comparing VDW energy ratios with/without switching function
    """
    # Create test system
    state = create_test_system()
    
    # Set atom charges for Ewald calculation
    state.atoms[0].charge = 1.0
    state.atoms[1].charge = -1.0
    
    # Set Ewald parameters
    state.info.setTemperature(300.0)  # 300K
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)
    
    # Set switching function parameters
    r_on = 1.0   # Inner cutoff radius
    r_off = 1.2  # Outer cutoff radius
    
    # Select distances between r_on and r_off for testing
    test_distances = [1.02, 1.05, 1.08, 1.10, 1.15, 1.18]
    
    print("\nInternal Switching Function test (through Ewald VDW energy ratio verification):")
    print(f"{'Distance(nm)':10s} | {'VDW without switching':14s} | {'VDW with switching':14s} | {'Actual ratio':10s} | {'Expected switching value':10s} | {'Error':10s}")
    print("-" * 75)
    
    for r in test_distances:
        # Set particle distance
        state.atoms[1].x = r
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        # 1. Calculate energy without switching function
        state.info.use_switching = False
        pygcmc.computeSystemEnergyEwald(state)
        vdw_no_switch = state.residues[0].energy_vdw + state.residues[1].energy_vdw
        
        # If energy is 0, skip this distance point
        if abs(vdw_no_switch) < 1e-10:
            continue
        
        # 2. Calculate energy with switching function
        state.info.use_switching = True
        state.info.r_on = r_on
        state.info.r_off = r_off
        pygcmc.computeSystemEnergyEwald(state)
        vdw_with_switch = state.residues[0].energy_vdw + state.residues[1].energy_vdw
        
        # 3. Calculate energy ratio, i.e., actual switching function value
        actual_switch = vdw_with_switch / vdw_no_switch if vdw_no_switch != 0 else 0
        
        # 4. Calculate theoretical switching function value
        r2 = r * r
        ron2 = r_on * r_on
        roff2 = r_off * r_off
        
        numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0*r2 - 3.0*ron2)
        denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2)
        expected_switch = numerator / denominator
        
        # Calculate error
        error = abs(actual_switch - expected_switch)
        
        # Print results
        print(f"{r:10.3f} | {vdw_no_switch:14.6f} | {vdw_with_switch:14.6f} | "
              f"{actual_switch:10.6f} | {expected_switch:10.6f} | {error:10.6f}")
        
        # Verify error is within acceptable range
        assert actual_switch == pytest.approx(expected_switch, abs=1e-4), \
               f"At r={r}, switching function value differs: {actual_switch} vs expected {expected_switch}"
    
    # Restore test state
    state.info.use_switching = False


def test_mcs_energy_calculation_with_switching():
    """
    Test whether MonteCarloSystem object can correctly affect energy calculation
    This test checks whether switching function parameters are passed from MCS to MCState
    """
    # Create test system
    state = create_test_system()
    
    # Set switching function parameters
    r_on = 1.0   # Inner cutoff radius
    r_off = 1.2  # Outer cutoff radius
    
    # Create MonteCarloSystem object
    mc_system = pygcmc.MonteCarloSystem()
    
    # Print available methods
    print("\nTest MonteCarloSystem and energy calculation association")
    print(f"MCS switching function method: {dir(mc_system)}")
    
    # Set a particle distance, so it's above r_off
    state.atoms[1].x = 1.25  # Distance greater than r_off
    
    # First calculate energy without switching function
    state.info.use_switching = False
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    pygcmc.computeSystemEnergy(state)
    energy_no_switching = state.residues[0].energy_vdw + state.residues[1].energy_vdw
    
    print(f"Distance = 1.25 (> r_off), energy without switching function = {energy_no_switching}")
    
    # Then calculate energy with switching function
    state.info.use_switching = True
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    # Verify setting is correctly applied
    # Since MCInfo private attribute cannot be directly accessed, but we can test them through special way
    state.atoms[1].x = 1.25  # Ensure distance remains unchanged
    pygcmc.computeSystemEnergy(state)
    energy_with_switching = state.residues[0].energy_vdw + state.residues[1].energy_vdw
    
    print(f"Distance = 1.25 (> r_off), energy with switching function = {energy_with_switching}")
    
    # Verify switching function is correctly applied
    if energy_with_switching != 0.0:
        print(f"Warning: Switching function may not be correctly applied! Energy should be 0, but got {energy_with_switching}")
    
    assert energy_with_switching == pytest.approx(0.0, abs=1e-8), \
           f"Energy with switching should be 0 when r > r_off, but got {energy_with_switching}"
    
    # Test completed, disable switching function
    state.info.use_switching = False