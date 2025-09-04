"""
Test PGP performance under strict tolerance standards

Recreate original PGP test scenario but with stricter tolerance standards
"""

import sys
import os
import random
import statistics

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import (setPGPParameters, initializePMEParameters, 
                    precomputeGridPotential, computeSystemEnergyPGP,
                    calculateMoleculeEnergy, computeMovementEnergyPME,
                    computeSystemEnergyPME, computeSystemEnergyEwald)


def test_pgp_with_strict_tolerance():
    """Test PGP performance under different tolerance standards"""
    
    print("\n" + "="*80)
    print("PGP Strict Tolerance Test")
    print("="*80)
    
    # Create a simple test system
    box_size = 5.0  # nm
    cutoff = 1.8    # nm
    
    # Fixed molecules
    fixed_positions = [
        [1.0, 1.0, 2.5],
        [4.0, 1.0, 2.5],
        [4.0, 4.0, 2.5],
        [1.0, 4.0, 2.5]
    ]
    fixed_charges = [1.0, -1.0, 1.0, -1.0]
    
    # Moving molecule (initial position)
    moving_position_initial = [2.5, 2.5, 2.5]
    moving_charge = 0.5
    
    # Create system
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]    # No LJ, only test electrostatics
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Add atoms
    all_positions_initial = fixed_positions + [moving_position_initial]
    all_charges = fixed_charges + [moving_charge]
    
    atoms = []
    for i in range(5):
        atom = MCAtom()
        atom.x, atom.y, atom.z = all_positions_initial[i]
        atom.charge = all_charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # Add residues
    residues = []
    for i in range(5):
        res = MCResidue()
        res.active = True
        res.fixed = (i < 4)  # First 4 are fixed
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 5
    
    # Initialize parameters
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    initializePMEParameters(cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, cutoff, mesh_size, 4, 1e-6)
    
    # Initialize Ewald parameters
    pygcmc.initializeEwaldParameters(cutoff, state.info.box, alpha)
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Set movement residues
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    print("\nTesting multiple random movements...")
    print("-" * 60)
    
    errors_pgp_pme = []
    errors_pgp_ewald = []
    n_moves = 10
    
    for move_idx in range(n_moves):
        # Calculate initial energy
        initial_pgp = calculateMoleculeEnergy(state)
        initial_pme_result = computeMovementEnergyPME(state)
        initial_pme_recip = initial_pme_result[2].get('reciprocal', 0.0)
        initial_ewald_result = computeSystemEnergyEwald(state)
        initial_ewald_recip = initial_ewald_result[2].get('reciprocal', 0.0)
        
        # Random movement (small displacement)
        dx = random.uniform(-0.3, 0.3)
        dy = random.uniform(-0.3, 0.3)
        dz = random.uniform(-0.3, 0.3)
        
        # Move atom
        moving_atom = state.atoms[4]
        moving_atom.x = (moving_position_initial[0] + dx) % box_size
        moving_atom.y = (moving_position_initial[1] + dy) % box_size
        moving_atom.z = (moving_position_initial[2] + dz) % box_size
        
        # Calculate energy after movement
        moved_pgp = calculateMoleculeEnergy(state)
        moved_pme_result = computeMovementEnergyPME(state)
        moved_pme_recip = moved_pme_result[2].get('reciprocal', 0.0)
        moved_ewald_result = computeSystemEnergyEwald(state)
        moved_ewald_recip = moved_ewald_result[2].get('reciprocal', 0.0)
        
        # Calculate energy change
        delta_pgp = moved_pgp - initial_pgp
        delta_pme = moved_pme_recip - initial_pme_recip
        delta_ewald = moved_ewald_recip - initial_ewald_recip
        
        print(f"\nMove {move_idx + 1}: ({dx:.3f}, {dy:.3f}, {dz:.3f})")
        print(f"ΔE_PGP:   {delta_pgp:10.6f} kJ/mol")
        print(f"ΔE_PME:   {delta_pme:10.6f} kJ/mol")
        print(f"ΔE_Ewald: {delta_ewald:10.6f} kJ/mol")
        
        # Calculate relative error
        if abs(delta_pme) > 1e-6:
            error_pgp_pme = abs((delta_pgp - delta_pme) / delta_pme)
            errors_pgp_pme.append(error_pgp_pme)
            print(f"PGP vs PME error: {error_pgp_pme*100:.2f}%")
        
        if abs(delta_ewald) > 1e-6:
            error_pgp_ewald = abs((delta_pgp - delta_ewald) / delta_ewald)
            errors_pgp_ewald.append(error_pgp_ewald)
            print(f"PGP vs Ewald error: {error_pgp_ewald*100:.2f}%")
        
        # Restore original position
        moving_atom.x = moving_position_initial[0]
        moving_atom.y = moving_position_initial[1]
        moving_atom.z = moving_position_initial[2]
    
    # Analyze results
    print("\n" + "="*80)
    print("Error Analysis")
    print("="*80)
    
    if errors_pgp_pme:
        avg_error = statistics.mean(errors_pgp_pme) * 100
        max_error = max(errors_pgp_pme) * 100
        min_error = min(errors_pgp_pme) * 100
        
        print(f"\nPGP vs PME error statistics:")
        print(f"Average error: {avg_error:.2f}%")
        print(f"Maximum error: {max_error:.2f}%")
        print(f"Minimum error: {min_error:.2f}%")
        
        # Test different tolerance standards
        tolerances = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
        print(f"\nTest results under different tolerance standards:")
        print("-" * 40)
        
        for tol in tolerances:
            passed = avg_error/100 < tol
            status = "✓ Pass" if passed else "✗ Fail"
            print(f"Tolerance {tol*100:5.1f}%: {status}")
        
        # Check how many individual tests would fail
        print(f"\nIndividual test pass rates:")
        for tol in tolerances:
            n_passed = sum(1 for e in errors_pgp_pme if e < tol)
            pass_rate = n_passed / len(errors_pgp_pme) * 100
            print(f"Tolerance {tol*100:5.1f}%: {n_passed}/{len(errors_pgp_pme)} ({pass_rate:.1f}%)")
    
    # Test system total energy
    print("\n" + "="*80)
    print("System Total Energy Test")
    print("="*80)
    
    # Calculate system total energy
    computeSystemEnergyPGP(state)
    pgp_total = state.ewald_energy.get('total', 0.0)
    pgp_recip = state.ewald_energy.get('reciprocal', 0.0)
    
    computeSystemEnergyPME(state)
    pme_total = state.ewald_energy.get('total', 0.0)
    pme_recip = state.ewald_energy.get('reciprocal', 0.0)
    
    print(f"\nPME reciprocal space: {pme_recip:.2f} kJ/mol")
    print(f"PGP reciprocal space: {pgp_recip:.2f} kJ/mol")
    print(f"Ratio PGP/PME: {pgp_recip/pme_recip if pme_recip != 0 else 0:.3f}")
    
    if abs(pgp_recip/pme_recip - 2.0) < 0.01:
        print("\n⚠️  Found PGP reciprocal space energy is 2x PME!")
    
    print("\n" + "="*80)
    print("Conclusion")
    print("="*80)
    
    if errors_pgp_pme and avg_error > 10:
        print("✗ PGP fails test under strict tolerance standards")
        print(f"  Average error {avg_error:.1f}% far exceeds reasonable range")
    else:
        print("? Need more test data")


if __name__ == "__main__":
    test_pgp_with_strict_tolerance()