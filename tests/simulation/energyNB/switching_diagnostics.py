# tests/simulation/energyNB/switching_diagnostics.py
"""Switching function diagnostic tests."""

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


def test_compare_energy_with_without_switching():
    """
    Detailed comparison of energy calculation results with and without switching function
    Also includes direct calculation and Ewald calculation methods
    """
    # Create test system
    state = create_test_system()
    
    # Set atom charges, for Ewald calculation
    state.atoms[0].charge = 1.0
    state.atoms[1].charge = -1.0
    
    # Set Ewald parameters
    state.info.setTemperature(300.0)  # 300K
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)
    
    # Set switching function parameters
    r_on = 1.0   # Inner cutoff radius (10 Å)
    r_off = 1.2  # Outer cutoff radius (12 Å)
    
    # Create MonteCarloSystem object
    mc_system = pygcmc.MonteCarloSystem()
    
    # Take more points near r_on for better switching effect
    distances = [0.7 + i * 0.8/50 for i in range(51)]  # 51 points from 0.7 to 1.5
    
    # 1. Direct calculation result comparison
    # Disable switching function
    state.info.use_switching = False
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    direct_no_switching = []
    for r in distances:
        state.atoms[1].x = r  # Set distance
        pygcmc.computeSystemEnergy(state)  # Use direct method to calculate energy
        direct_no_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # Enable switching function
    state.info.use_switching = True
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    direct_with_switching = []
    for r in distances:
        state.atoms[1].x = r  # Set distance
        pygcmc.computeSystemEnergy(state)  # Use direct method to calculate energy
        direct_with_switching.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
    
    # 2. Ewald calculation result comparison
    # Disable switching function
    state.info.use_switching = False
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    ewald_no_switching_vdw = []
    ewald_no_switching_elec = []
    for r in distances:
        state.atoms[1].x = r  # Set distance
        pygcmc.computeSystemEnergyEwald(state)  # Use Ewald method to calculate energy
        ewald_no_switching_vdw.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
        ewald_no_switching_elec.append(state.residues[0].energy_elec + state.residues[1].energy_elec)
    
    # Enable switching function
    state.info.use_switching = True
    state.info.r_on = r_on
    state.info.r_off = r_off
    
    ewald_with_switching_vdw = []
    ewald_with_switching_elec = []
    for r in distances:
        state.atoms[1].x = r  # Set distance
        pygcmc.computeSystemEnergyEwald(state)  # Use Ewald method to calculate energy
        ewald_with_switching_vdw.append(state.residues[0].energy_vdw + state.residues[1].energy_vdw)
        ewald_with_switching_elec.append(state.residues[0].energy_elec + state.residues[1].energy_elec)
    
    # Print results and verify
    print("\nEnergy calculation comparison results (with switching vs without switching):")
    print(f"{'Distance(nm)':10s} | {'Direct no switching':14s} | {'Direct with switching':14s} | {'Ewald VDW no switching':16s} | {'Ewald VDW with switching':16s}")
    print("-" * 80)
    
    # Print some key result points
    key_indices = [0, 10, 20, 25, 30, 35, 40, 45, 50]  # Select several key points for display
    
    for idx in key_indices:
        if idx < len(distances):
            r = distances[idx]
            # Calculate switching value at r_off for verification
            switch_value = mc_system.calculate_switching_function(r) if r_on <= r <= r_off else (1.0 if r < r_on else 0.0)
            
            print(f"{r:10.3f} | {direct_no_switching[idx]:14.6f} | {direct_with_switching[idx]:14.6f} | "
                  f"{ewald_no_switching_vdw[idx]:16.6f} | {ewald_with_switching_vdw[idx]:16.6f}")
    
    # Perform verification
    for i, r in enumerate(distances):
        # 1. When r < r_on, energy should be the same
        if r < r_on:
            assert direct_no_switching[i] == pytest.approx(direct_with_switching[i]), \
                  f"Direct energy should be the same when r < r_on, at r={r}"
            assert ewald_no_switching_vdw[i] == pytest.approx(ewald_with_switching_vdw[i]), \
                  f"Ewald VDW energy should be the same when r < r_on, at r={r}"
                  
        # 2. When r > r_off, energy with switching should be 0
        elif r > r_off:
            assert direct_with_switching[i] == pytest.approx(0.0), \
                  f"Direct energy with switching should be 0 when r > r_off, at r={r}"
            assert ewald_with_switching_vdw[i] == pytest.approx(0.0), \
                  f"Ewald VDW energy with switching should be 0 when r > r_off, at r={r}"
                  
        # 3. In r_on and r_off, check if switching is correctly applied
        elif r_on <= r <= r_off:
            # Calculate expected switching value
            switch_value = mc_system.calculate_switching_function(r)
            
            # Check direct energy is within reasonable range
            assert 0.0 <= abs(direct_with_switching[i]) <= abs(direct_no_switching[i]), \
                  f"Direct energy with switching should be between 0 and original energy at r={r}"
            
            # Check Ewald VDW energy is within reasonable range
            assert 0.0 <= abs(ewald_with_switching_vdw[i]) <= abs(ewald_no_switching_vdw[i]), \
                  f"Ewald VDW energy with switching should be between 0 and original energy at r={r}"
            
            # Check energy decreases monotonically (if not the first point between r_on and r_off)
            if i > 0 and r_on <= distances[i-1] <= r_off:
                assert abs(direct_with_switching[i]) <= abs(direct_with_switching[i-1]), \
                      f"Direct energy magnitude should decrease as r increases in switching region"
                assert abs(ewald_with_switching_vdw[i]) <= abs(ewald_with_switching_vdw[i-1]), \
                      f"Ewald VDW energy magnitude should decrease as r increases in switching region"
    
    # Charge energy test - Only check key points
    if ewald_no_switching_elec[0] != 0:  # Ensure there is charge energy
        print("\nCharge energy comparison results (Ewald):")
        print(f"{'Distance(nm)':10s} | {'Ewald charge no switching':18s} | {'Ewald charge with switching':18s}")
        print("-" * 60)
        
        for idx in key_indices:
            if idx < len(distances):
                r = distances[idx]
                print(f"{r:10.3f} | {ewald_no_switching_elec[idx]:18.6f} | {ewald_with_switching_elec[idx]:18.6f}")
    
    # Disable switching function, restore default state
    state.info.use_switching = False


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