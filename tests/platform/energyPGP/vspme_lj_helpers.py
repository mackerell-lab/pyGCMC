# tests/simulation/energyPGP/vspme_lj_helpers.py
"""Helper functions for LJ energy calculation tests."""

import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import computeSystemVdwEnergyCutoff, computeMovementEnergyPME, computeMovementEnergyPGP
from pygcmc import setPMEParameters, setPGPParameters, initializePMEParameters, precomputeGridPotential


def create_lj_test_system(box_size, cutoff):
    """Create a simple two-atom system for LJ energy testing."""
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = cutoff

    # Set more reasonable LJ parameters to avoid extreme energy values
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    # Use more conservative parameters: larger sigma, smaller epsilon
    sigma_1 = 0.6  # nm - increase sigma
    sigma_2 = 0.6  # nm - use same value to simplify
    eps_1 = 0.001  # kJ/mol - significantly reduce epsilon
    eps_2 = 0.001  # kJ/mol - use same value to simplify
    
    # Set LJ parameter matrix
    combined_sigma = (sigma_1 + sigma_2) / 2.0
    combined_eps = math.sqrt(eps_1 * eps_2)
    
    ff.ljSigma = [
        sigma_1, combined_sigma,
        combined_sigma, sigma_2
    ]
    ff.ljEps = [
        eps_1, combined_eps,
        combined_eps, eps_2
    ]
    state.forcefield = ff

    # Create atoms
    atoms = []
    
    # Fixed central atom
    fixed_atom = MCAtom()
    fixed_atom.x = 2.0  # Place at box center
    fixed_atom.y = 2.0
    fixed_atom.z = 2.0
    fixed_atom.charge = 0.0  # No charge, focus on LJ interactions
    fixed_atom.type = 0
    atoms.append(fixed_atom)
    
    # Create fixed residue
    residues = []
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 1
    fixed_res.active = True
    fixed_res.fixed = True
    fixed_res.energy_vdw = 0.0
    fixed_res.energy_elec = 0.0
    residues.append(fixed_res)
    
    # Create a moving atom at a safe distance
    moving_atom = MCAtom()
    moving_atom.x = 2.0 + 1.2  # Initial distance 1.2nm, away from repulsion zone
    moving_atom.y = 2.0
    moving_atom.z = 2.0
    moving_atom.charge = 0.0  # No charge, focus on LJ interactions
    moving_atom.type = 1
    atoms.append(moving_atom)
    
    # Create moving residue
    moving_res = MCResidue()
    moving_res.atomStart = 1
    moving_res.atomCount = 1
    moving_res.active = True
    moving_res.fixed = False
    moving_res.energy_vdw = 0.0
    moving_res.energy_elec = 0.0
    residues.append(moving_res)
    
    # Set system state
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state, combined_sigma, combined_eps, sigma_1, sigma_2, eps_1, eps_2


def test_cumulative_energy_changes(state, fixed_atom, combined_sigma, combined_eps):
    """Test cumulative energy changes from continuous movement."""
    print("\n3. Test cumulative energy changes from continuous movement")
    
    # Reset atom position to initial state
    state.atoms[1].x = 2.0 + 1.2  # Initial distance 1.2nm
    state.atoms[1].y = 2.0
    state.atoms[1].z = 2.0
    
    # Ensure reset of movement energy reference
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Calculate initial PME/PGP energies to set reference point
    computeMovementEnergyPME(state)
    computeMovementEnergyPGP(state)
    
    # Record cumulative energy changes
    cumulative_theoretical = 0.0
    cumulative_pme = 0.0
    cumulative_pgp = 0.0
    
    print("\nDist(nm) | Step(nm) | Theor ΔLJ  | PME ΔLJ    | PGP ΔLJ    | Cumul Th   | Cumul PME  | Cumul PGP  | PME Rel Err | PGP Rel Err")
    print("---------+----------+------------+------------+------------+------------+------------+------------+-------------+------------")
    
    # Initial distance
    dx = state.atoms[1].x - fixed_atom.x
    dy = state.atoms[1].y - fixed_atom.y
    dz = state.atoms[1].z - fixed_atom.z
    current_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    prev_distance = current_distance
    
    # Gradually decrease distance
    step_sizes = [0.1, 0.1, 0.1, 0.1]  # Move 0.1nm each time
    
    for step_size in step_sizes:
        # Move atom
        state.atoms[1].x -= step_size  # Move towards fixed atom
        
        # Calculate new distance
        dx = state.atoms[1].x - fixed_atom.x
        dy = state.atoms[1].y - fixed_atom.y
        dz = state.atoms[1].z - fixed_atom.z
        current_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        # Calculate theoretical LJ energy (current and previous state)
        r = current_distance
        sigma = combined_sigma
        eps = combined_eps
        r_over_sigma = r / sigma
        r6 = 1.0 / (r_over_sigma**6)
        r12 = r6 * r6
        current_theoretical_lj = 4.0 * eps * (r12 - r6)
        
        r_prev = prev_distance
        r_over_sigma_prev = r_prev / sigma
        r6_prev = 1.0 / (r_over_sigma_prev**6)
        r12_prev = r6_prev * r6_prev
        prev_theoretical_lj = 4.0 * eps * (r12_prev - r6_prev)
        
        # Theoretical LJ energy change
        theoretical_delta = current_theoretical_lj - prev_theoretical_lj
        
        # Calculate PME energy change
        pme_result = computeMovementEnergyPME(state)
        pme_total = pme_result[0]
        pme_lj = pme_result[1]
        
        # Calculate PGP energy change
        pgp_result = computeMovementEnergyPGP(state)
        pgp_total = pgp_result[0]
        pgp_lj = pgp_result[1]
        
        # Cumulative energy changes
        cumulative_theoretical += theoretical_delta
        cumulative_pme += pme_lj
        cumulative_pgp += pgp_lj
        
        # Calculate relative errors
        step_pme_error = abs((pme_lj - theoretical_delta) / theoretical_delta) if abs(theoretical_delta) > 1e-10 else float('nan')
        step_pgp_error = abs((pgp_lj - theoretical_delta) / theoretical_delta) if abs(theoretical_delta) > 1e-10 else float('nan')
        
        cumul_pme_error = abs((cumulative_pme - cumulative_theoretical) / cumulative_theoretical) if abs(cumulative_theoretical) > 1e-10 else float('nan')
        cumul_pgp_error = abs((cumulative_pgp - cumulative_theoretical) / cumulative_theoretical) if abs(cumulative_theoretical) > 1e-10 else float('nan')
        
        # Print results
        print(f"{current_distance:7.4f} | {step_size:8.4f} | {theoretical_delta:10.6f} | {pme_lj:10.6f} | {pgp_lj:10.6f} | {cumulative_theoretical:10.6f} | {cumulative_pme:10.6f} | {cumulative_pgp:10.6f} | {cumul_pme_error:11.2%} | {cumul_pgp_error:10.2%}")
        
        # Update previous distance
        prev_distance = current_distance
    
    return cumul_pme_error, cumul_pgp_error


def print_test_summary(cumul_pme_error, cumul_pgp_error):
    """Print test summary and conclusions."""
    print("\n--- LJ Energy Calculation Test Summary ---")
    print("1. Movement functions (PME and PGP) calculate energy changes rather than absolute energies")
    print("2. Initial state is set as energy reference point (zero point)")
    print("3. After each movement, Movement functions return energy change relative to previous state")
    print("4. Cumulative energy changes from continuous movements can be compared with theoretical predictions")
    
    # If there are serious errors, print warning
    if cumul_pme_error > 0.1 or cumul_pgp_error > 0.1:
        print("\nWarning: PME or PGP LJ energy calculations have significant differences from theoretical values")
        print("Possible causes:")
        print("- Movement distance steps too large causing numerical issues")
        print("- Cumulative errors in Movement functions")
        print("- LJ potential calculations sensitive to distance changes")
    else:
        print("\nValidation passed: PME and PGP Movement functions can correctly calculate LJ energy changes")
    
    print("--- Test Complete ---")
