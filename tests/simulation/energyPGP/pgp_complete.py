"""
Test PGP Complete implementation
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import math


def create_test_system():
    """Create a test system with charges and LJ"""
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field with both charges and LJ
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create a system with fixed and moveable particles
    particle_data = [
        # (position, charge, fixed)
        ([2.0, 2.5, 2.5], 1.0, True),    # Fixed positive
        ([3.0, 2.5, 2.5], -1.0, True),   # Fixed negative
        ([2.5, 2.5, 2.5], 0.5, False),   # Moveable
    ]
    
    for i, (pos, charge, fixed) in enumerate(particle_data):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = fixed
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    # Set up movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2  # The moveable particle
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    return state


def pgp_complete_vs_pme_complete():
    """Test that PGP Complete gives similar results to PME Complete"""
    
    print("\n" + "="*70)
    print("PGP Complete vs PME Complete Test")
    print("="*70)
    
    state = create_test_system()
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    # Initialize PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Initialize PGP
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    
    # Precompute PGP grid
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate with PME Complete
    elec_pme, vdw_pme, total_pme = pygcmc.computeSystemEnergyPMEComplete(state)
    print(f"\nPME Complete:")
    print(f"  Electrostatic: {elec_pme:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pme:.6f} kJ/mol")
    print(f"  Total:        {total_pme:.6f} kJ/mol")
    
    # Calculate with PGP Complete
    elec_pgp, vdw_pgp, total_pgp = pygcmc.computeSystemEnergyPGPComplete(state)
    print(f"\nPGP Complete:")
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pgp:.6f} kJ/mol")
    print(f"  Total:        {total_pgp:.6f} kJ/mol")
    
    # Compare
    elec_diff = abs(elec_pgp - elec_pme)
    vdw_diff = abs(vdw_pgp - vdw_pme)
    total_diff = abs(total_pgp - total_pme)
    
    print(f"\nDifferences:")
    print(f"  Electrostatic: {elec_diff:.6f} kJ/mol")
    print(f"  VdW:          {vdw_diff:.6f} kJ/mol")
    print(f"  Total:        {total_diff:.6f} kJ/mol")
    
    # VdW should be identical
    assert vdw_diff < 1e-6, f"VdW energies should be identical, got difference {vdw_diff}"
    
    # Electrostatic should be very close
    rel_elec_error = elec_diff / abs(elec_pme) * 100 if elec_pme != 0 else 0
    print(f"\nRelative electrostatic error: {rel_elec_error:.4f}%")
    assert rel_elec_error < 0.1, f"Electrostatic error too large: {rel_elec_error:.4f}%"


def pgp_complete_movement_energy():
    """Test PGP Complete movement energy calculation"""
    
    print("\n" + "="*70)
    print("PGP Complete Movement Energy Test")
    print("="*70)
    
    state = create_test_system()
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    # Initialize
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Use system energy to test movement residue changes
    # Get initial system energies
    elec_pme_init, vdw_pme_init, total_pme_init = pygcmc.computeSystemEnergyPMEComplete(state)
    elec_pgp_init, vdw_pgp_init, total_pgp_init = pygcmc.computeSystemEnergyPGPComplete(state)
    
    print(f"\nInitial system energies:")
    print(f"  PME Complete: {total_pme_init:.6f} kJ/mol")
    print(f"  PGP Complete: {total_pgp_init:.6f} kJ/mol")
    
    # Move the particle
    state.atoms[2].x += 0.2
    
    # Get final system energies
    elec_pme_final, vdw_pme_final, total_pme_final = pygcmc.computeSystemEnergyPMEComplete(state)
    elec_pgp_final, vdw_pgp_final, total_pgp_final = pygcmc.computeSystemEnergyPGPComplete(state)
    
    print(f"\nFinal system energies:")
    print(f"  PME Complete: {total_pme_final:.6f} kJ/mol")
    print(f"  PGP Complete: {total_pgp_final:.6f} kJ/mol")
    
    # Calculate energy changes
    pme_delta = total_pme_final - total_pme_init
    pgp_delta = total_pgp_final - total_pgp_init
    
    print(f"\nEnergy changes:")
    print(f"  PME Complete: {pme_delta:.6f} kJ/mol")
    print(f"  PGP Complete: {pgp_delta:.6f} kJ/mol")
    print(f"  Difference:   {abs(pme_delta - pgp_delta):.6f} kJ/mol")
    
    # The energy changes should be very similar
    rel_error = abs(pme_delta - pgp_delta) / abs(pme_delta) * 100 if pme_delta != 0 else 0
    print(f"  Relative error: {rel_error:.4f}%")
    
    # Since PGP precomputes fixed potential, the difference might be larger but should still be small
    assert rel_error < 5.0, f"Energy change error too large: {rel_error:.4f}%"


def pgp_complete_pure_lj():
    """Test PGP Complete with pure LJ system"""
    
    print("\n" + "="*70)
    print("PGP Complete Pure LJ Test")
    print("="*70)
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field with LJ only
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.4]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Two particles, no charges
    positions = [[2.0, 2.5, 2.5], [2.5, 2.5, 2.5]]
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # No charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = (i == 0)
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energies
    elec_pme, vdw_pme, total_pme = pygcmc.computeSystemEnergyPMEComplete(state)
    elec_pgp, vdw_pgp, total_pgp = pygcmc.computeSystemEnergyPGPComplete(state)
    
    print(f"\nPME Complete:")
    print(f"  Electrostatic: {elec_pme:.6f} kJ/mol (should be 0)")
    print(f"  VdW:          {vdw_pme:.6f} kJ/mol")
    
    print(f"\nPGP Complete:")
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol (should be 0)")
    print(f"  VdW:          {vdw_pgp:.6f} kJ/mol")
    
    # Electrostatic should be zero
    assert abs(elec_pme) < 1e-6, f"PME electrostatic should be zero, got {elec_pme}"
    assert abs(elec_pgp) < 1e-6, f"PGP electrostatic should be zero, got {elec_pgp}"
    
    # VdW should match
    vdw_diff = abs(vdw_pgp - vdw_pme)
    print(f"\nVdW difference: {vdw_diff:.6f} kJ/mol")
    assert vdw_diff < 1e-6, f"VdW energies should match, got difference {vdw_diff}"


if __name__ == "__main__":
    pgp_complete_vs_pme_complete()
    pgp_complete_movement_energy()
    pgp_complete_pure_lj()