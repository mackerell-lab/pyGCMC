"""
Validate PGP design principle: PGP should accurately calculate the
interaction energy between moving atoms and fixed atoms.

Key insight: PGP's "reciprocal" energy is NOT the same as PME reciprocal.
PGP calculates the energy of moving atoms in the potential field created
by fixed atoms, which is exactly what GCMC needs.
"""

import pytest
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

# Only import OpenMM if available
try:
    import openmm as mm
    from openmm import app
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    pytest.skip("OpenMM not available", allow_module_level=True)


def test_pgp_moving_fixed_interaction():
    """Test that PGP accurately calculates moving-fixed interactions"""
    
    print("\n" + "="*60)
    print("Testing PGP Moving-Fixed Interaction Accuracy")
    print("="*60)
    
    # Create system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+, Cl-
    ff.numMovementTypes = 2
    eps_na = 0.5216
    eps_cl = 0.5216
    sig_na = 0.2160
    sig_cl = 0.4830
    
    ff.ljEps = [
        eps_na,
        (eps_na * eps_cl)**0.5,
        (eps_na * eps_cl)**0.5,
        eps_cl
    ]
    ff.ljSigma = [
        sig_na,
        (sig_na + sig_cl) / 2.0,
        (sig_na + sig_cl) / 2.0,
        sig_cl
    ]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed molecules (4 atoms)
    fixed_positions = [
        ([1.0, 1.0, 2.5], 1.0, 0),   # Na+
        ([1.5, 1.0, 2.5], -1.0, 1),  # Cl-
        ([3.5, 3.5, 2.5], 1.0, 0),   # Na+
        ([4.0, 3.5, 2.5], -1.0, 1),  # Cl-
    ]
    
    for i, (pos, charge, atom_type) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = i
        res.atomCount = 1
        res.active = True
        res.fixed = True
        res.type = atom_type
        residues.append(res)
    
    # Moving molecule (2 atoms)
    moving_positions = [
        ([2.5, 2.5, 2.5], 1.0, 0),   # Na+
        ([2.7, 2.5, 2.5], -1.0, 1),  # Cl-
    ]
    
    for i, (pos, charge, atom_type) in enumerate(moving_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = len(fixed_positions) + i
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = atom_type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Setup movement residues
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4
    movement_info.activeCount = 2
    state.movementResidues.append(movement_info)
    
    # Initialize PME/PGP
    alpha = 2.84
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    
    # Precompute grid for fixed atoms
    print("\nPrecomputing PGP grid for fixed atoms...")
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Test multiple configurations
    print("\nTesting different moving atom positions...")
    print("-" * 50)
    
    test_positions = [
        ([2.5, 2.5, 2.5], [2.7, 2.5, 2.5]),  # Original
        ([2.0, 2.0, 2.0], [2.2, 2.0, 2.0]),  # Closer to fixed
        ([3.0, 3.0, 3.0], [3.2, 3.0, 3.0]),  # Different location
        ([1.5, 2.5, 3.5], [1.7, 2.5, 3.5]),  # Another location
    ]
    
    for pos_idx, (pos1, pos2) in enumerate(test_positions):
        # Set positions
        state.atoms[4].x, state.atoms[4].y, state.atoms[4].z = pos1
        state.atoms[5].x, state.atoms[5].y, state.atoms[5].z = pos2
        
        # Calculate moving-fixed interaction using PME
        # First, get total system energy
        pme_total = pygcmc.computeSystemEnergyPMEComplete(state)
        
        # Then get fixed-only energy by temporarily deactivating moving atoms
        for i in range(4, 6):
            state.residues[i].active = False
        pme_fixed = pygcmc.computeSystemEnergyPMEComplete(state)
        
        # Reactivate moving atoms
        for i in range(4, 6):
            state.residues[i].active = True
        
        # Now get moving-only energy
        for i in range(4):
            state.residues[i].active = False
        pme_moving = pygcmc.computeSystemEnergyPMEComplete(state)
        
        # Reactivate all
        for i in range(6):
            state.residues[i].active = True
        
        # Moving-fixed interaction = Total - Fixed - Moving
        pme_mov_fix_interaction = pme_total[0] - pme_fixed[0] - pme_moving[0]
        
        # Calculate using PGP interpolation
        pgp_energy = pygcmc.calculateMoleculeEnergy(state)
        
        # The PGP energy should match the moving-fixed PME interaction
        print(f"\nPosition {pos_idx + 1}:")
        print(f"  PME moving-fixed interaction: {pme_mov_fix_interaction:.6f} kJ/mol")
        print(f"  PGP interpolation energy:      {pgp_energy:.6f} kJ/mol")
        
        if abs(pme_mov_fix_interaction) > 0.1:
            error = abs((pgp_energy - pme_mov_fix_interaction) / pme_mov_fix_interaction)
            print(f"  Relative error: {error:.2%}")
            
            # PGP should approximate this interaction reasonably well
            assert error < 0.15, f"Error {error:.2%} exceeds 15% threshold"
    
    print("\n✓ PGP correctly approximates moving-fixed interactions!")
    
    # Now test ΔE accuracy
    print("\n\nTesting ΔE accuracy for GCMC moves...")
    print("-" * 50)
    
    # Reset to first position
    state.atoms[4].x, state.atoms[4].y, state.atoms[4].z = test_positions[0][0]
    state.atoms[5].x, state.atoms[5].y, state.atoms[5].z = test_positions[0][1]
    
    # Initial energy
    pgp_initial = pygcmc.calculateMoleculeEnergy(state)
    
    # Calculate initial PME moving-fixed interaction
    pme_total_init = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(4, 6):
        state.residues[i].active = False
    pme_fixed_init = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(4, 6):
        state.residues[i].active = True
    for i in range(4):
        state.residues[i].active = False
    pme_moving_init = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(6):
        state.residues[i].active = True
    pme_mov_fix_init = pme_total_init[0] - pme_fixed_init[0] - pme_moving_init[0]
    
    # Move to second position
    state.atoms[4].x, state.atoms[4].y, state.atoms[4].z = test_positions[1][0]
    state.atoms[5].x, state.atoms[5].y, state.atoms[5].z = test_positions[1][1]
    
    # Final energy
    pgp_final = pygcmc.calculateMoleculeEnergy(state)
    
    # Calculate final PME moving-fixed interaction
    pme_total_final = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(4, 6):
        state.residues[i].active = False
    pme_fixed_final = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(4, 6):
        state.residues[i].active = True
    for i in range(4):
        state.residues[i].active = False
    pme_moving_final = pygcmc.computeSystemEnergyPMEComplete(state)
    for i in range(6):
        state.residues[i].active = True
    pme_mov_fix_final = pme_total_final[0] - pme_fixed_final[0] - pme_moving_final[0]
    
    # Compare ΔE
    pgp_delta = pgp_final - pgp_initial
    pme_delta = pme_mov_fix_final - pme_mov_fix_init
    
    print(f"\nΔE comparison:")
    print(f"  PME moving-fixed ΔE: {pme_delta:.6f} kJ/mol")
    print(f"  PGP ΔE:              {pgp_delta:.6f} kJ/mol")
    
    if abs(pme_delta) > 0.01:
        delta_error = abs((pgp_delta - pme_delta) / pme_delta)
        print(f"  ΔE relative error: {delta_error:.2%}")
        
        # For GCMC, ΔE accuracy is most important
        assert delta_error < 0.10, f"ΔE error {delta_error:.2%} exceeds 10% threshold"
    
    print("\n✓ PGP provides accurate ΔE for GCMC moves!")
    
    print("\n" + "="*60)
    print("PGP validation complete: Design goals achieved!")
    print("="*60)


if __name__ == "__main__":
    test_pgp_moving_fixed_interaction()