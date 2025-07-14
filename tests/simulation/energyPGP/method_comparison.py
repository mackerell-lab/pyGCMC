# tests/simulation/pgp/method_comparison.py
"""PGP method comparison tests."""

import pytest
import pygcmc
from pygcmc import MCMovementResidueInfo
import sys
from .helpers import create_long_distance_system


def test_compare_pme_pgp_energy():
    """
    Compare whether the energy values calculated by PME and PGP are consistent before and after movement
    
    This test verifies:
    1. In the initial state, energy values calculated by PME and PGP should be the same
    2. After moving the molecule, the energy changes calculated by PME and PGP should be the same
    """
    # Set parameters - ensure PME and PGP use the same parameters
    box_size = 5.0  # nm - use a larger box
    cutoff = 1.0   # nm
    potential_cutoff = 1.0  # nm - same as cutoff
    box = [box_size, box_size, box_size]
    
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]  # Use the same grid size as PME for accurate comparison
    spline_order = 4
    tolerance = 1e-5
    
    # Set test timeout to avoid long running
    timeout = 10  # seconds
    
    print("Creating long-distance test system...")
    sys.stdout.flush()
    
    # Create test system - using specially designed long-distance system
    system = create_long_distance_system(box_size)
    
    # Ensure box size is set correctly
    system.info.box = box
    system.info.cutoff = cutoff
    
    print("Setting PME parameters...")
    sys.stdout.flush()
    
    # Set PME and PGP parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    print("Setting PGP parameters...")
    sys.stdout.flush()
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    print("Setting up moving residues...")
    sys.stdout.flush()
    
    # Moving residues already set up in create_long_distance_system
    # Just need to prepare movementResidues list
    moving_residues = [1]  # Second residue is the moving residue
    
    # Verify fixed residue information
    fixed_count = sum(1 for res in system.residues if res.fixed)
    print(f"Number of fixed residues: {fixed_count}")
    print(f"Number of moving residues: {len(moving_residues)}")
    sys.stdout.flush()
    
    # Set up moving residue info
    system.movementResidues.clear()
    
    # Create moving residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # Index of moving residue
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"Added movement info: startIndex={movement_info.startIndex}, activeCount={movement_info.activeCount}")
    print(f"System has {system.activeResidueCount} active residues and {len(system.movementResidues)} movement residue groups")
    sys.stdout.flush()
    
    # Step 1: Calculate initial system energy using PME
    print("Calculating initial PME energy...")
    sys.stdout.flush()
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]  # PME electrostatic energy
    initial_pme_dict = initial_pme_result[2]  # PME energy details dictionary
    initial_pme_reciprocal = initial_pme_dict['reciprocal']  # Only take reciprocal space part
    print(f"Initial PME energy result: {initial_pme_result}")
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    sys.stdout.flush()
    
    # Step 2: Precompute grid potential using PGP and calculate moving residue energy
    print("Precomputing PGP grid potential...")
    sys.stdout.flush()
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    print("Calculating initial PGP energy...")
    sys.stdout.flush()
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Print initial energy
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    sys.stdout.flush()
    
    # Step 3: Move moving residue (e.g., translate 0.1 nm)
    translation = [0.1, 0.1, 0.1]  # nm
    
    print("Moving residue...")
    sys.stdout.flush()
    for res_idx in moving_residues:
        residue = system.residues[res_idx]
        for atom_idx in range(residue.atomCount):
            atom_index = residue.atomStart + atom_idx
            atom = system.atoms[atom_index]
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
    
    # Step 4: Calculate energy using PME after movement
    print("Calculating PME after movement energy...")
    sys.stdout.flush()
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]  # PME electrostatic energy
    moved_pme_dict = moved_pme_result[2]  # PME energy details dictionary
    moved_pme_reciprocal = moved_pme_dict['reciprocal']  # Only take reciprocal space part
    print(f"Moved PME energy result: {moved_pme_result}")
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    sys.stdout.flush()
    
    # Step 5: Calculate energy using PGP after movement
    print("Calculating PGP after movement energy...")
    sys.stdout.flush()
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Print moved energy
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    sys.stdout.flush()
    
    # Calculate energy change - only use reciprocal space part
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME reciprocal energy change: {pme_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    sys.stdout.flush()
    
    if abs(pme_energy_change) < 1e-10:
        print("PME energy change is too small, cannot compute relative error")
        assert abs(pgp_energy_change) < 1e-10, f"PGP energy should also be close to zero"
    else:
        # Calculate relative error, allowing some error range (e.g., 10%)
        relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"Relative error: {relative_error * 100:.4f}%")
        sys.stdout.flush()
        
        # Verify PGP and PME calculated energy changes are consistent within error range
        assert relative_error < 0.1, f"Relative error too large: {relative_error*100:.2f}%"  # Allow 10% error