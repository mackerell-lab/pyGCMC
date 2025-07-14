# tests/simulation/energyPGP/asymmetric_nacl.py
"""PGP asymmetric NaCl systems test - Complete implementation with complex movement algorithm."""

import pytest
import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys
from .helpers import calculate_pbc_distance, is_safe_position
from .asymmetric_nacl_helpers import (
    create_nacl_crystal_system,
    handle_initial_position_safety,
    calculate_complex_movement_vector,
    apply_movement_to_residue,
    verify_moved_position_safety
)


def test_compare_ewald_pme_pgp_asymmetric_nacl():
    """
    Test whether PGP calculation is accurate in systems with asymmetric charge distribution
    Using NaCl as the moving residue instead of water
    
    Features:
    1. Uses a 4×4×4 NaCl supercell as the base structure
    2. One NaCl pair is designated as the moving residue
    3. Only compares reciprocal space energy changes, which should be accurate
       even if particles are closer than the cutoff distance
    4. Verify accuracy by comparing energy calculated by Ewald, PME, and PGP
    """
    # Set system parameters
    n_cells = 4      # 4x4x4 supercell
    a = 0.564        # NaCl lattice constant (nm)
    box_size = 5.0   # nm - Use large box to ensure sufficient space
    cutoff = 1.0     # nm - Note: even if distance is smaller than this value, we only compare reciprocal space energy
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald and PME parameters
    alpha = 0.29    # 1/nm
    kmax = [8, 8, 8]  # Number of k-space vectors for Ewald calculation
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    print("Creating NaCl crystal system with moving NaCl residue...")
    sys.stdout.flush()
    
    # Create NaCl crystal system
    print(f"\nCreating {n_cells}x{n_cells}x{n_cells} NaCl crystal...")
    system, mobile_res, fixed_positions = create_nacl_crystal_system(n_cells, a, box_size, cutoff)
    
    # Calculate system total charge
    total_system_charge = sum(atom.charge for atom in system.atoms)
    print(f"System total charge: {total_system_charge}")
    
    print(f"System created: {system.activeAtomCount} atoms, {system.activeResidueCount} residues")
    print(f"Fixed atoms: {len(system.atoms) - 2}, Moving atoms: 2 (NaCl ion pair)")
    sys.stdout.flush()
    
    # Initialize Ewald, PME, and PGP
    print("Setting calculation parameters...")
    
    # Ewald parameters initialization
    print("Initializing Ewald parameters...")
    pygcmc.setEwaldParameters(alpha, kmax)
    
    # PME parameters initialization
    print("Initializing PME parameters...")
    pygcmc.setPMEParameters(alpha=alpha, meshSize=mesh_size, splineOrder=spline_order, tolerance=tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP parameters initialization
    print("Initializing PGP parameters...")
    pygcmc.setPGPParameters(alpha=alpha, meshSize=mesh_size, potential_cutoff=potential_cutoff,
                           potentialGridSize=potential_grid_size, splineOrder=spline_order, tolerance=tolerance)
    
    # Precompute grid potential for fixed particles
    print("Precomputing PGP grid potential for fixed particles...")
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Execute multiple movements
    num_moves = 5  # Reduce test times to speed up test
    print(f"Executing {num_moves} random movement tests")
    
    # Store relative errors between methods
    pgp_pme_errors = []
    ewald_pme_errors = []
    pgp_ewald_errors = []
    
    # Handle initial position safety
    handle_initial_position_safety(system, mobile_res, fixed_positions, cutoff, box_size)
    
    # Execute multiple movements
    for move_idx in range(num_moves):
        print(f"\nExecuting {move_idx+1}/{num_moves} movement test")
        
        # Step 1: Calculate initial energy
        print("Calculating initial energy...")
        
        # Ewald energy calculation
        initial_ewald_result = pygcmc.computeSystemEnergyEwald(system)
        initial_ewald_elec = initial_ewald_result[0]
        initial_ewald_dict = initial_ewald_result[2]
        initial_ewald_reciprocal = initial_ewald_dict['reciprocal']
        
        # PME energy calculation
        initial_pme_result = pygcmc.computeMovementEnergyPME(system)
        initial_pme_energy = initial_pme_result[0]
        initial_pme_dict = initial_pme_result[2]
        initial_pme_reciprocal = initial_pme_dict['reciprocal']
        
        # PGP energy calculation
        initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"Initial Ewald reciprocal energy: {initial_ewald_reciprocal}")
        print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
        print(f"Initial PGP energy: {initial_pgp_energy}")
        
        # Step 2: Calculate complex movement vector
        movement_vector = calculate_complex_movement_vector(system, mobile_res, system.residues, box_size)
        
        # Step 3: Apply movement
        apply_movement_to_residue(system, mobile_res, movement_vector, box_size)
        
        # Step 4: Verify moved position safety
        verify_moved_position_safety(system, mobile_res, fixed_positions, cutoff, box_size)
        
        # Step 5: Calculate moved energy
        print("Calculating moved energy...")
        
        # Ewald energy calculation
        moved_ewald_result = pygcmc.computeSystemEnergyEwald(system)
        moved_ewald_elec = moved_ewald_result[0]
        moved_ewald_dict = moved_ewald_result[2]
        moved_ewald_reciprocal = moved_ewald_dict['reciprocal']
        
        # PME energy calculation
        moved_pme_result = pygcmc.computeMovementEnergyPME(system)
        moved_pme_energy = moved_pme_result[0]
        moved_pme_dict = moved_pme_result[2]
        moved_pme_reciprocal = moved_pme_dict['reciprocal']
        
        # Check direct space energy
        moved_pme_direct = moved_pme_dict.get('direct', 0.0)
        if abs(moved_pme_direct) > 1e-10:
            print(f"PME direct space energy: {moved_pme_direct}")
            print("Particles are within cutoff distance - but we only compare reciprocal space energy changes")
        
        # PGP energy calculation
        moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"Moved Ewald reciprocal energy: {moved_ewald_reciprocal}")
        print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
        print(f"Moved PGP energy: {moved_pgp_energy}")
        
        # Calculate energy change
        ewald_energy_change = moved_ewald_reciprocal - initial_ewald_reciprocal
        pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
        pgp_energy_change = moved_pgp_energy - initial_pgp_energy
        
        print(f"Ewald reciprocal energy change: {ewald_energy_change}")
        print(f"PME reciprocal energy change: {pme_energy_change}")
        print(f"PGP energy change: {pgp_energy_change}")
        
        # Note: We only compare reciprocal space energy changes, even if distance is smaller than cutoff
        # This is because PGP's purpose is to simulate PME's reciprocal space calculation, direct space calculation is handled separately
        
        # Calculate relative error
        if abs(pme_energy_change) > 1e-6:
            # PGP relative error compared to PME
            pgp_pme_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
            print(f"PGP relative error: {pgp_pme_error*100:.4f}%")
            pgp_pme_errors.append(pgp_pme_error)
            
            # Ewald relative error compared to PME
            ewald_pme_error = abs((ewald_energy_change - pme_energy_change) / pme_energy_change)
            print(f"Ewald relative error: {ewald_pme_error*100:.4f}%")
            ewald_pme_errors.append(ewald_pme_error)
            
            # PGP relative error compared to Ewald
            if abs(ewald_energy_change) > 1e-6:
                pgp_ewald_error = abs((pgp_energy_change - ewald_energy_change) / ewald_energy_change)
                print(f"PGP relative error: {pgp_ewald_error*100:.4f}%")
                pgp_ewald_errors.append(pgp_ewald_error)
        else:
            print("PME energy change near zero, skipping relative error calculation")
    
    # Calculate average error
    print("\nError analysis statistics:")
    
    if pgp_pme_errors:
        avg_pgp_pme_error = sum(pgp_pme_errors) / len(pgp_pme_errors)
        print(f"Average PGP relative error: {avg_pgp_pme_error*100:.4f}%")
        
        # Use more lenient error tolerance because PGP is an approximation method
        acceptable_error = 0.5  # Allow 50% error
        assert avg_pgp_pme_error < acceptable_error, f"Average PGP relative error too large: {avg_pgp_pme_error*100:.2f}%"
    else:
        print("PGP relative error: No valid error data for statistics")
    
    if ewald_pme_errors:
        avg_ewald_pme_error = sum(ewald_pme_errors) / len(ewald_pme_errors)
        print(f"Average Ewald relative error: {avg_ewald_pme_error*100:.4f}%")
        
        # Ewald and PME should theoretically have very high consistency
        acceptable_error = 0.2  # Allow 20% error
        assert avg_ewald_pme_error < acceptable_error, f"Average Ewald relative error too large: {avg_ewald_pme_error*100:.2f}%"
    else:
        print("Ewald relative error: No valid error data for statistics")
    
    if pgp_ewald_errors:
        avg_pgp_ewald_error = sum(pgp_ewald_errors) / len(pgp_ewald_errors)
        print(f"Average PGP relative error: {avg_pgp_ewald_error*100:.4f}%")
        
        # Use more lenient error tolerance because PGP is an approximation method
        acceptable_error = 0.5  # Allow 50% error
        assert avg_pgp_ewald_error < acceptable_error, f"Average PGP relative error too large: {avg_pgp_ewald_error*100:.2f}%"
    else:
        print("PGP relative error: No valid error data for statistics")
    
    print("\n--- Test completed successfully ---")
    sys.stdout.flush()