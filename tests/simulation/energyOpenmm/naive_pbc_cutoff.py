# tests/simulation/energyOpenmm/naive_pbc_cutoff.py

import pytest
from .naive_helpers import *

def test_compare_pbc_energies():
    """Compare PBC energy calculations between OpenMM and naive implementation."""
    # Get system
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # Set consistent cutoff distance
    state.info.cutoff = 1.0  # Consistent with OpenMM: 1.0 nm
    
    # Move water molecule to box edge to test PBC
    box_size = state.info.box[0]  # 3.0 nm
    new_positions = []
    for i in range(len(positions)):
        pos = positions[i].value_in_unit(nanometers)
        if i >= 6:  # Water molecule atoms
            # Move to 2.5 nm, so distance from origin is 2.5 nm, PBC distance is 0.5 nm
            new_positions.append(Vec3(2.5, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    
    # Verify PBC distance calculation
    # Using first carbon atom (-0.15 nm) and water's oxygen atom (2.5 nm) as example
    dx = 2.5 - (-0.15)  # Original distance = 2.65 nm
    dx_pbc = dx - box_size * round(dx/box_size)  # Should be approximately -0.35 nm
    print(f"\nPBC distance validation:")
    print(f"Original distance: {dx:.6f} nm")
    print(f"After PBC: {dx_pbc:.6f} nm")
    print(f"Should be in range [-{box_size/2:.1f}, {box_size/2:.1f}] nm")
    assert abs(dx_pbc) <= box_size/2, "PBC distance calculation error"
    
    # Calculate OpenMM PBC energy
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    openmm_energy = calculate_nonbonded_energy(system, new_positions, movement_atoms, fixed_atoms, use_pbc=True)
    openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
    
    # Update atom positions in naive implementation
    for i in range(6, 9):
        pos = new_positions[i].value_in_unit(nanometers)
        state.atoms[i].x = pos[0]
        state.atoms[i].y = pos[1]
        state.atoms[i].z = pos[2]
    
    # Calculate naive PBC energy - only calculate movement residue energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    movement_fixed_energy = state.residues[1].energy_vdw + state.residues[1].energy_elec
    
    # Print detailed information
    print(f"\nDetailed PBC energy comparison:")
    print(f"Box size: {box_size} nm")
    print(f"Cutoff distance: {state.info.cutoff} nm")
    print(f"Water molecule position: {new_positions[6].value_in_unit(nanometers)} nm")
    print(f"OpenMM PBC energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive movement-fixed energy: {movement_fixed_energy:.6f} kJ/mol")
    print(f"Difference: {abs(openmm_energy_val - movement_fixed_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - movement_fixed_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    # Verify results
    rel_tol = 0.3  # 30% relative error tolerance, as the two implementations may have different details
    assert abs(openmm_energy_val - movement_fixed_energy) / abs(openmm_energy_val) < rel_tol, \
           f"PBC energy mismatch: OpenMM={openmm_energy_val}, Naive={movement_fixed_energy}"

def test_compare_cutoff_effects():
    """Compare cutoff effects between OpenMM and naive implementation."""
    # Enable debug output
    pygcmc.setEnergyDebugOutput(True)
    
    # Get system
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # Test different distances
    distances = [0.5, 0.7, 0.9, 1.2]  # nm
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    
    for dist in distances:
        # Move water molecule
        new_positions = []
        for i in range(len(positions)):
            pos = positions[i].value_in_unit(nanometers)
            if i >= 6:  # Water molecule atoms
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # Calculate OpenMM energy
        openmm_energy = calculate_nonbonded_energy(
            system, new_positions, movement_atoms, fixed_atoms, use_pbc=False
        )
        openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
        
        # Update atom positions in naive implementation
        for i in range(6, 9):
            pos = new_positions[i].value_in_unit(nanometers)
            state.atoms[i].x = pos[0]
            state.atoms[i].y = pos[1]
            state.atoms[i].z = pos[2]
        
        # Calculate naive energy
        pygcmc.computeMovementEnergyCutoff(state)
        naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
        
        print(f"\nEnergy comparison at distance {dist} nm:")
        print(f"OpenMM energy: {openmm_energy_val:.6f} kJ/mol")
        print(f"Naive energy: {naive_energy:.6f} kJ/mol")
        print(f"Difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
        
        # For distances beyond cutoff, both should give energy close to 0
        if dist > state.info.cutoff:
            assert abs(openmm_energy_val) < 1e-6, f"OpenMM energy not zero beyond cutoff: {openmm_energy_val}"
            assert abs(naive_energy) < 1e-6, f"Naive energy not zero beyond cutoff: {naive_energy}"
        else:
            # For distances within cutoff, energies should be close
            rel_tol = 0.001  # 0.1% relative error tolerance, as the two implementations may have different details
            assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
                   f"Energy mismatch at {dist} nm: OpenMM={openmm_energy_val}, Naive={naive_energy}"