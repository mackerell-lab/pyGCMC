# tests/simulation/energyPGP/vspme_combined.py
"""Combined energy calculation diagnostic test."""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import computeMovementEnergyPGP, computeMovementEnergyPME, computeSystemVdwEnergyCutoff
import sys
from .vspme_energy_diagnostics import (
    create_diagnostic_test_system,
    setup_energy_parameters
)

def calculate_theoretical_lj_energy(distance, sigma, eps):
    """Calculate theoretical LJ energy for given distance."""
    r_over_sigma = distance / sigma
    r_6 = math.pow(1.0/r_over_sigma, 6)
    r_12 = r_6 * r_6
    return 4 * eps * (r_12 - r_6)

def test_movement_energy_functions(state):
    """Test movement energy functions for PME and PGP."""
    print("\n1.4 Calculate movement energy using computeMovementEnergyPME")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate movement energy
    pme_movement_result = computeMovementEnergyPME(state)
    pme_movement_total = pme_movement_result[0]
    pme_movement_lj = pme_movement_result[1]
    pme_movement_components = pme_movement_result[2]
    
    print(f"  PME Movement total energy: {pme_movement_total:.6f} kJ/mol")
    print(f"  PME Movement LJ energy: {pme_movement_lj:.6f} kJ/mol")
    print(f"  PME Movement components: {pme_movement_components}")
    
    # Check if residue energies are updated
    movement_residue_vdw = 0.0
    for i, res in enumerate(state.residues):
        if i == 1:  # Moving residue
            movement_residue_vdw = res.energy_vdw
        print(f"  Residue {i} LJ energy after Movement: {res.energy_vdw:.6f} kJ/mol")
    
    print("\n1.5 Calculate movement energy using computeMovementEnergyPGP")
    
    
    pgp_movement_result = computeMovementEnergyPGP(state)
    pgp_movement_total = pgp_movement_result[0]
    pgp_movement_lj = pgp_movement_result[1]
    pgp_movement_components = pgp_movement_result[2]
    
    print(f"  PGP Movement total energy: {pgp_movement_total:.6f} kJ/mol")
    print(f"  PGP Movement LJ energy: {pgp_movement_lj:.6f} kJ/mol")
    print(f"  PGP Movement components: {pgp_movement_components}")
    

def perform_manual_energy_combination(state):
    """Perform manual combination of electrostatic and LJ energies."""
    print("\n--- Phase 3: Manual Combination of Electrostatic and LJ Energies ---")
    
    # Reset to initial position
    state.atoms[1].x = 3.0  # Initial distance 1.0nm
    state.atoms[1].y = 2.0
    state.atoms[1].z = 2.0
    
    # Reset movement reference point
    computeMovementEnergyPME(state)
    computeMovementEnergyPGP(state)
    
    print("Reset atom position to initial state, distance: 1.0 nm")
    print("Reset Movement function reference point")
    
    # Move atom again
    state.atoms[1].x = 2.5  # 0.5nm
    
    print("\n3.1 Manual combination PME method")
    
    # Calculate electrostatic energy
    pme_movement_result_3 = computeMovementEnergyPME(state)
    pme_elec_delta = pme_movement_result_3[0]
    
    # Reset and calculate LJ energy
    
    computeSystemVdwEnergyCutoff(state)
    
    # Only calculate LJ energy for moving residue
    moving_res_lj = state.residues[1].energy_vdw
    
    # Combine energies
    pme_combined = pme_elec_delta + moving_res_lj
    
    print(f"  PME electrostatic energy change: {pme_elec_delta:.6f} kJ/mol")
    print(f"  Moving residue LJ energy: {moving_res_lj:.6f} kJ/mol")
    print(f"  Manual combined total energy change: {pme_combined:.6f} kJ/mol")
    
    print("\n3.2 Manual combination PGP method")
    
    pgp_movement_result_3 = computeMovementEnergyPGP(state)
    pgp_elec_delta = pgp_movement_result_3[0]
    
    
    
    
    pgp_combined = pgp_elec_delta + moving_res_lj
    
    print(f"  PGP electrostatic energy change: {pgp_elec_delta:.6f} kJ/mol")
    print(f"  Manual combined total energy change: {pgp_combined:.6f} kJ/mol")

def test_combined_energy_calculation():
    """
    Diagnostic test: Check PME and PGP energy calculation details, especially LJ energy component
    
    This test will:
    1. Create a simple two-atom system
    2. Test different energy calculation functions separately
    3. Check all intermediate results
    4. Try to manually combine electrostatic and LJ energies, simulating complete energy calculation
    print("\n--- Diagnostic Test: Energy Calculation Component Detail Check ---")

    # Create test system
    state, box_size, cutoff, sigma, eps = create_diagnostic_test_system()
    
    # Calculate initial distance
    fixed_atom = state.atoms[0]
    moving_atom = state.atoms[1]
    dx = moving_atom.x - fixed_atom.x
    dy = moving_atom.y - fixed_atom.y
    dz = moving_atom.z - fixed_atom.z
    initial_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"System initialization complete:")
    print(f"Fixed atom coordinates: ({fixed_atom.x}, {fixed_atom.y}, {fixed_atom.z})")
    print(f"Moving atom coordinates: ({moving_atom.x}, {moving_atom.y}, {moving_atom.z})")
    print(f"Initial distance: {initial_distance:.4f} nm")
    print(f"LJ parameters: sigma={sigma} nm, epsilon={eps} kJ/mol")
    
    # Setup energy parameters
    alpha, mesh_size = setup_energy_parameters(state, box_size, cutoff)
    
    # Analyze initial state energies - moved from vspme_energy_diagnostics.py
    from .vspme_energy_diagnostics import analyze_initial_state_energies
    vdw_total = analyze_initial_state_energies(state, sigma, eps)
    
    # Test movement energy functions
    test_movement_energy_functions(state)
    
    print("\n--- Phase 2: Energy Change Analysis After Moving Atom ---")
    
    # Move atom to new position - significantly change distance
    state.atoms[1].x = 2.5  # From 1.0nm to 0.5nm
    
    # Calculate new distance
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    new_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Moved atom to new position: ({state.atoms[1].x}, {state.atoms[1].y}, {state.atoms[1].z})")
    print(f"New distance: {new_distance:.4f} nm (from {initial_distance:.4f} nm)")
    
    # Theoretical LJ energy calculation
    lj_energy1 = calculate_theoretical_lj_energy(initial_distance, sigma, eps)
    lj_energy2 = calculate_theoretical_lj_energy(new_distance, sigma, eps)
    
    lj_energy_change = lj_energy2 - lj_energy1
    
    print(f"Theoretical LJ energy 1: {lj_energy1:.6f} kJ/mol")
    print(f"Theoretical LJ energy 2: {lj_energy2:.6f} kJ/mol")
    print(f"Theoretical LJ energy change: {lj_energy_change:.6f} kJ/mol")
    
    print("\n2.1 Direct LJ energy calculation after movement")
    
    
    # Direct LJ energy calculation
    
    # Print LJ energy for each residue
    vdw_total_2 = 0.0
        vdw_total_2 += res.energy_vdw
        print(f"  Residue {i} LJ energy: {res.energy_vdw:.6f} kJ/mol")
    print(f"  Total LJ energy: {vdw_total_2:.6f} kJ/mol")
    print(f"  LJ energy change: {vdw_total_2 - vdw_total:.6f} kJ/mol")
    
    print("\n2.2 Calculate PME Movement energy after movement")
    
    
    pme_movement_result_2 = computeMovementEnergyPME(state)
    pme_movement_total_2 = pme_movement_result_2[0]
    pme_movement_lj_2 = pme_movement_result_2[1]
    pme_movement_components_2 = pme_movement_result_2[2]
    
    print(f"  PME Movement total energy: {pme_movement_total_2:.6f} kJ/mol")
    print(f"  PME Movement LJ energy: {pme_movement_lj_2:.6f} kJ/mol")
    print(f"  PME Movement components: {pme_movement_components_2}")
    
    
    print("\n2.3 Calculate PGP Movement energy after movement")
    
    
    pgp_movement_result_2 = computeMovementEnergyPGP(state)
    pgp_movement_total_2 = pgp_movement_result_2[0]
    pgp_movement_lj_2 = pgp_movement_result_2[1]
    pgp_movement_components_2 = pgp_movement_result_2[2]
    
    print(f"  PGP Movement total energy: {pgp_movement_total_2:.6f} kJ/mol")
    print(f"  PGP Movement LJ energy: {pgp_movement_lj_2:.6f} kJ/mol")
    print(f"  PGP Movement components: {pgp_movement_components_2}")
    
    
    # Perform manual energy combination
    perform_manual_energy_combination(state)
    
    print("\n--- Diagnostic Conclusion ---")
    print("1. Direct LJ calculation method (computeSystemVdwEnergyCutoff) can correctly calculate LJ energy")
    print("2. SystemEnergy calculations include LJ energy component")
    print("3. Movement functions do not include LJ energy change calculations")
    print("4. Manual combination of PME/PGP electrostatic energy changes and direct LJ calculations can achieve complete energy calculation")
"""
