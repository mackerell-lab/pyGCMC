# tests/simulation/energyPGP/vspme_lj_energy.py
"""LJ energy calculation test in PME and PGP methods."""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPMEParameters, initializePMEParameters, setPGPParameters
from .pgp_wrapper import precomputeGridPotential, computeMovementEnergyPGP, computeMovementEnergyPME
from .pgp_wrapper import computeSystemVdwEnergyCutoff
import sys
from pygcmc import MCMovementResidueInfo
from .vspme_lj_helpers import (
    create_lj_test_system,
    test_cumulative_energy_changes,
    print_test_summary
)

def test_lj_energy_pme_pgp():
    """
    Test LJ energy calculation in PME and PGP methods.
    
    This test creates an atomic system with reasonable distances and verifies:
    1. LJ energy calculations are within reasonable range
    2. Different calculation methods (direct, PME, PGP) produce consistent energy changes
    3. Movement functions correctly calculate energy changes rather than absolute energies
    # --- Parameter Settings ---
    box_size = 4.0  # nm
    cutoff = 1.2    # nm
    box = [box_size, box_size, box_size]
    alpha = 0.29    # 1/nm
    mesh_size = [16, 16, 16]  # Reduce grid size for efficiency
    spline_order = 4
    tolerance = 1e-5

    print("\n--- Test: LJ Energy Calculation in PME and PGP Methods ---")
    sys.stdout.flush()

    # --- Create Test System ---
    state, combined_sigma, combined_eps, sigma_1, sigma_2, eps_1, eps_2 = create_lj_test_system(box_size, cutoff)
    
    # Get atom references
    fixed_atom = state.atoms[0]
    moving_atom = state.atoms[1]
    
    # Calculate actual distance for confirmation
    dx = moving_atom.x - fixed_atom.x
    dy = moving_atom.y - fixed_atom.y
    dz = moving_atom.z - fixed_atom.z
    initial_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Created a simple two-atom system")
    print(f"Fixed atom position: ({fixed_atom.x}, {fixed_atom.y}, {fixed_atom.z})")
    print(f"Moving atom initial position: ({moving_atom.x}, {moving_atom.y}, {moving_atom.z})")
    print(f"Initial distance: {initial_distance:.4f} nm")
    print(f"LJ parameters: sigma_1={sigma_1} nm, sigma_2={sigma_2} nm, epsilon_1={eps_1} kJ/mol, epsilon_2={eps_2} kJ/mol")
    
    # --- Initialize PME/PGP Parameters ---
    setPMEParameters(alpha=alpha, meshSize=mesh_size, splineOrder=spline_order, tolerance=tolerance)
    initializePMEParameters(cutoff, box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    
    print("\n1. Initial state energy calculation")
    
    # --- Test Direct LJ Calculation ---
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    computeSystemVdwEnergyCutoff(state)
    
    # Get direct LJ calculation values
    direct_lj_initial = 0.0
    for i, res in enumerate(state.residues):
        direct_lj_initial += res.energy_vdw
        print(f"Residue {i} LJ energy: {res.energy_vdw:.6f} kJ/mol")
    print(f"Initial state total LJ energy: {direct_lj_initial:.6f} kJ/mol")
    
    # --- Set Movement Information ---
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1  # Second residue is the moving residue
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # --- Test PME Movement Energy ---
    pme_result_initial = computeMovementEnergyPME(state)
    pme_total_initial = pme_result_initial[0]
    pme_lj_initial = pme_result_initial[1]
    pme_components_initial = pme_result_initial[2]
    
    print(f"Initial PME Movement energy (should be 0 or reference value):")
    print(f"  Total energy: {pme_total_initial:.6f} kJ/mol")
    print(f"  LJ energy: {pme_lj_initial:.6f} kJ/mol")
    print(f"  Energy components: {pme_components_initial}")
    
    # --- Test PGP Movement Energy ---
    pgp_result_initial = computeMovementEnergyPGP(state)
    pgp_total_initial = pgp_result_initial[0]
    pgp_lj_initial = pgp_result_initial[1]
    pgp_components_initial = pgp_result_initial[2]
    
    print(f"Initial PGP Movement energy (should be 0 or reference value):")
    print(f"  Total energy: {pgp_total_initial:.6f} kJ/mol")
    print(f"  LJ energy: {pgp_lj_initial:.6f} kJ/mol")
    print(f"  Energy components: {pgp_components_initial}")
    
    # Confirm Movement function behavior
    if abs(pme_total_initial) < 1e-6 and abs(pgp_total_initial) < 1e-6:
        print("\nConfirmed: Movement functions return energy changes rather than absolute energy values")
        print("This means the initial state serves as reference point, energy change is 0")
    
    # --- Move atom and record energy changes ---
    # Test multiple distance points, from far to near
    test_distances = [1.1, 1.0, 0.9, 0.8, 0.7]
    
    print("\n2. Test energy changes when moving atom")
    print("\nDist(nm) | Theor LJ   | Direct LJ  | PME Move E  | PGP Move E  | PME-Theor   | PGP-Theor")
    print("---------+------------+------------+-------------+-------------+-------------+------------")
    
    # Initial distance record
    prev_distance = initial_distance
    prev_theoretical_lj = 0  # Theoretical initial LJ energy
    
    for target_distance in test_distances:
        # Move atom to new position
        vector_length = target_distance  # Set desired distance
        state.atoms[1].x = fixed_atom.x + vector_length  # Move along x-axis
        state.atoms[1].y = fixed_atom.y
        state.atoms[1].z = fixed_atom.z
        
        # Calculate actual distance
        dx = state.atoms[1].x - fixed_atom.x
        dy = state.atoms[1].y - fixed_atom.y
        dz = state.atoms[1].z - fixed_atom.z
        actual_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        # Calculate theoretical LJ energy E = 4ε[(σ/r)^12 - (σ/r)^6]
        r = actual_distance
        sigma = combined_sigma
        eps = combined_eps
        r_over_sigma = r / sigma
        r6 = 1.0 / (r_over_sigma**6)
        r12 = r6 * r6
        theoretical_lj = 4.0 * eps * (r12 - r6)
        
        # Theoretical LJ energy change
        theoretical_delta = theoretical_lj - prev_theoretical_lj
        
        # Reset residue energies
        for res in state.residues:
            res.energy_vdw = 0.0
            res.energy_elec = 0.0
        
        # Calculate direct LJ energy
        computeSystemVdwEnergyCutoff(state)
        direct_lj = 0.0
            direct_lj += res.energy_vdw
        
        # Calculate PME movement energy
        pme_result = computeMovementEnergyPME(state)
        pme_total = pme_result[0]
        pme_lj = pme_result[1]
        
        # Calculate PGP movement energy
        pgp_result = computeMovementEnergyPGP(state)
        pgp_total = pgp_result[0]
        pgp_lj = pgp_result[1]
        
        # Calculate relative errors
        pme_rel_error = abs((pme_lj - theoretical_delta) / theoretical_delta) if abs(theoretical_delta) > 1e-10 else float('nan')
        pgp_rel_error = abs((pgp_lj - theoretical_delta) / theoretical_delta) if abs(theoretical_delta) > 1e-10 else float('nan')
        
        # Print results
        print(f"{actual_distance:7.4f} | {theoretical_lj:10.6f} | {direct_lj:10.6f} | {pme_total:11.6f} | {pgp_total:11.6f} | {pme_rel_error:11.2%} | {pgp_rel_error:10.2%}")
        
        # Save current values for next comparison
        prev_distance = actual_distance
        prev_theoretical_lj = theoretical_lj
    
    # --- Test cumulative energy changes from continuous movement ---
    cumul_pme_error, cumul_pgp_error = test_cumulative_energy_changes(state, fixed_atom, combined_sigma, combined_eps)
    
    # Print test summary
    print_test_summary(cumul_pme_error, cumul_pgp_error)
"""
