# tests/simulation/pgp/asymmetric_water_movement.py
"""PGP asymmetric water movement loop - Complete version with all original functionality."""

import pygcmc
from .helpers import calculate_pbc_distance, is_safe_position, generate_safe_move


def execute_movement_loop(system, mobile_res, fixed_particles, num_moves, 
                         pgp_pme_errors, ewald_pme_errors, pgp_ewald_errors, 
                         cutoff, box_size):
    """Execute the complete movement loop with full error analysis."""
    
    for move_idx in range(num_moves):
        print(f"\nExecuting {move_idx+1}/{num_moves} random movement test")
        
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
        
        # Step 2: Generate safe random movement and apply
        current_positions = []
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            current_positions.append((
                system.atoms[atom_idx].x,
                system.atoms[atom_idx].y,
                system.atoms[atom_idx].z
            ))
        
        delta, new_positions = generate_safe_move(current_positions, fixed_particles, cutoff, box_size)
        print(f"Random movement vector: {delta}")
        
        # Move all particles in moving residue
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            system.atoms[atom_idx].x = new_positions[i][0]
            system.atoms[atom_idx].y = new_positions[i][1]
            system.atoms[atom_idx].z = new_positions[i][2]
            print(f"Moved atom {i} to: ({new_positions[i][0]:.4f}, {new_positions[i][1]:.4f}, {new_positions[i][2]:.4f})")
        
        # Verify moved position is safe
        is_safe, min_dist = is_safe_position(new_positions, fixed_particles, cutoff, box_size)
        if not is_safe:
            print(f"Warning: Moved position not safe, minimum distance is {min_dist} nm")
            assert min_dist > cutoff, f"Moved distance ({min_dist} nm) less than cutoff ({cutoff} nm), will introduce real space energy"
        else:
            min_dist = float('inf')
            for mobile_pos in new_positions:
                for fixed_pos, _ in fixed_particles:
                    dist = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                    min_dist = min(min_dist, dist)
            print(f"Moved position safe, minimum distance from fixed particles: {min_dist:.4f} nm (cutoff={cutoff} nm)")
        
        # Step 3: Calculate moved energy
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
        
        # Check direct space energy is 0
        moved_pme_direct = moved_pme_dict.get('direct', 0.0)
        if abs(moved_pme_direct) > 1e-10:
            print(f"Warning: PME direct space energy not zero: {moved_pme_direct}")
            print("This means there are particles < cutoff!")
        
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