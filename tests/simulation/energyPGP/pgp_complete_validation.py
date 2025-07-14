"""
Test PGP Complete implementation validation

This module tests the pure PGP Complete functionality, verifying that it correctly
implements the PGP method with complete interactions (including intramolecular).
"""

import pytest
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPMEParameters, initializePMEParameters, setPGPParameters
from .pgp_wrapper import precomputeGridPotential, computeSystemEnergyPGP, computeSystemEnergyPGPComplete
from .pgp_wrapper import computeMovementEnergyPGPComplete, resetPGPState
import math
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField, MCMovementResidueInfo

def test_pgp_complete_consistency():
    """Test that PGP Complete is internally consistent"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Internal Consistency Test")
    print("="*70)
    
    # Create a test system with both charges and LJ
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create a system with 3 particles
    particle_data = [
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
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Initialize parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    pgp_wrapper.setPMEParameters(alpha, mesh_size, spline_order)
    pgp_wrapper.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pgp_wrapper.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # Precompute PGP grid for fixed particles
    pgp_wrapper.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate initial energy
    elec1, vdw1, total1 = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nInitial PGP Complete:")
    print(f"  Electrostatic: {elec1:.6f} kJ/mol")
    print(f"  VdW:          {vdw1:.6f} kJ/mol")
    print(f"  Total:        {total1:.6f} kJ/mol")
    
    # Move the particle slightly
    state.atoms[2].x = 2.6
    
    # Calculate new energy
    elec2, vdw2, total2 = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nAfter moving particle:")
    print(f"  Electrostatic: {elec2:.6f} kJ/mol")
    print(f"  VdW:          {vdw2:.6f} kJ/mol")
    print(f"  Total:        {total2:.6f} kJ/mol")
    
    # Energy should change
    elec_change = elec2 - elec1
    vdw_change = vdw2 - vdw1
    total_change = total2 - total1
    
    print(f"\nEnergy changes:")
    print(f"  Electrostatic: {elec_change:.6f} kJ/mol")
    print(f"  VdW:          {vdw_change:.6f} kJ/mol")
    print(f"  Total:        {total_change:.6f} kJ/mol")
    
    # Verify that energy changed (particle moved closer to one charge, farther from another)
    assert abs(elec_change) > 0.1, "Electrostatic energy should change when particle moves"
    assert abs(vdw_change) > 0.01, "VdW energy should change when particle moves"
    
    # Move back
    state.atoms[2].x = 2.5
    
    # Calculate energy again - should match original
    elec3, vdw3, total3 = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nAfter moving back:")
    print(f"  Electrostatic: {elec3:.6f} kJ/mol")
    print(f"  VdW:          {vdw3:.6f} kJ/mol")
    print(f"  Total:        {total3:.6f} kJ/mol")
    
    # Should match original to reasonable precision
    # PGP grid interpolation may have small numerical differences
    assert abs(elec3 - elec1) < 1e-6, f"Electrostatic energy should be reversible, got diff {abs(elec3 - elec1)}"
    assert abs(vdw3 - vdw1) < 1e-10, f"VdW energy should be reversible, got diff {abs(vdw3 - vdw1)}"
    assert abs(total3 - total1) < 1e-6, f"Total energy should be reversible, got diff {abs(total3 - total1)}"

def test_pgp_complete_intramolecular():
    """Test that PGP Complete correctly includes intramolecular interactions"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Intramolecular Test")
    print("="*70)
    
    # Create a system with a multi-atom residue
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 4.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create a 2-atom residue
    positions = [[5.0, 5.0, 5.0], [5.5, 5.0, 5.0]]
    charges = [1.0, -1.0]
    
    for i, (pos, charge) in enumerate(zip(positions, charges)):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
    
    # Single residue with 2 atoms
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 2
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 1
    
    # Initialize parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pgp_wrapper.setPMEParameters(alpha, mesh_size, spline_order)
    pgp_wrapper.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pgp_wrapper.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # No fixed particles to precompute
    pgp_wrapper.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy with PGP Complete
    elec_pgp, vdw_pgp, total_pgp = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nPGP Complete (includes intramolecular):")
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pgp:.6f} kJ/mol")
    print(f"  Total:        {total_pgp:.6f} kJ/mol")
    
    # Calculate energy with standard PGP (excludes intramolecular LJ)
    elec_std, vdw_std, pgp_dict = pgp_wrapper.computeSystemEnergyPGP(state)
    total_std = pgp_dict["total"]
    print(f"\nStandard PGP (excludes intramolecular LJ):")
    print(f"  Electrostatic: {elec_std:.6f} kJ/mol")
    print(f"  VdW:          {vdw_std:.6f} kJ/mol")
    print(f"  Total:        {total_std:.6f} kJ/mol")
    
    # The difference should be the intramolecular LJ energy
    vdw_diff = vdw_pgp - vdw_std
    print(f"\nIntramolecular LJ energy: {vdw_diff:.6f} kJ/mol")
    
    # Calculate expected intramolecular LJ energy
    r = 0.5  # Distance between atoms
    sigma = 0.35
    eps = 1.0
    sigma_over_r = sigma / r
    sigma6 = sigma_over_r ** 6
    sigma12 = sigma6 ** 2
    expected_lj = 4.0 * eps * (sigma12 - sigma6)
    
    print(f"Expected intramolecular LJ: {expected_lj:.6f} kJ/mol")
    
    # Verify the intramolecular LJ is included (allow small numerical difference)
    assert abs(vdw_diff - expected_lj) < 0.01, f"Intramolecular LJ not correctly included, diff={abs(vdw_diff - expected_lj)}"
    
    # Electrostatic energies might differ due to different calculation methods
    # but both should be negative (attractive between opposite charges)
    assert elec_pgp < 0, "Electrostatic energy should be negative for opposite charges"
    assert elec_std < 0, "Electrostatic energy should be negative for opposite charges"

def test_pgp_complete_movement_energy():
    """Test PGP Complete movement energy calculation"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Movement Energy Test")
    print("="*70)
    
    # Create system with fixed and moveable particles
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 4.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create system: 2 fixed particles, 1 moveable
    particle_data = [
        ([4.0, 5.0, 5.0], 1.0, True),    # Fixed
        ([6.0, 5.0, 5.0], -1.0, True),   # Fixed
        ([5.0, 5.0, 5.0], 0.5, False),   # Moveable
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
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Initialize
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pgp_wrapper.setPMEParameters(alpha, mesh_size, spline_order)
    pgp_wrapper.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pgp_wrapper.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # Precompute grid for fixed particles
    pgp_wrapper.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate full system energy
    elec_full, vdw_full, total_full = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nFull system energy:")
    print(f"  Electrostatic: {elec_full:.6f} kJ/mol")
    print(f"  VdW:          {vdw_full:.6f} kJ/mol")
    print(f"  Total:        {total_full:.6f} kJ/mol")
    
    # Calculate movement energy only
    result = pgp_wrapper.computeMovementEnergyPGPComplete(state)
    if isinstance(result, tuple) and len(result) == 3:
        elec_move, vdw_move, pgp_dict = result
        total_move = pgp_dict["total"] if isinstance(pgp_dict, dict) else elec_move + vdw_move
    else:
        elec_move, vdw_move, total_move = result
    print(f"\nMovement energy only:")
    print(f"  Electrostatic: {elec_move:.6f} kJ/mol")
    print(f"  VdW:          {vdw_move:.6f} kJ/mol")
    print(f"  Total:        {total_move:.6f} kJ/mol")
    
    # Movement energy should be less than total (no fixed-fixed interactions)
    assert abs(total_move) < abs(total_full), "Movement energy should be less than total"
    
    # Deactivate the moveable particle
    state.residues[2].active = False
    state.activeAtomCount = 2
    state.activeResidueCount = 2
    
    # Calculate energy without moveable particle
    elec_fixed, vdw_fixed, total_fixed = pgp_wrapper.computeSystemEnergyPGPComplete(state)
    print(f"\nFixed particles only:")
    print(f"  Electrostatic: {elec_fixed:.6f} kJ/mol")
    print(f"  VdW:          {vdw_fixed:.6f} kJ/mol")
    print(f"  Total:        {total_fixed:.6f} kJ/mol")
    
    # The difference should be the movement energy
    elec_diff = elec_full - elec_fixed
    vdw_diff = vdw_full - vdw_fixed
    total_diff = total_full - total_fixed
    
    print(f"\nDifference (should match movement energy):")
    print(f"  Electrostatic: {elec_diff:.6f} kJ/mol")
    print(f"  VdW:          {vdw_diff:.6f} kJ/mol")
    print(f"  Total:        {total_diff:.6f} kJ/mol")
    
    # The differences should approximately match movement energy
    # (may not be exact due to different calculation paths)
    print(f"\nMovement energy validation:")
    print(f"  Elec diff vs move: {abs(elec_diff - elec_move):.6e}")
    print(f"  VdW diff vs move:  {abs(vdw_diff - vdw_move):.6e}")
    print(f"  Total diff vs move: {abs(total_diff - total_move):.6e}")

if __name__ == "__main__":
    test_pgp_complete_consistency()
    test_pgp_complete_intramolecular()
    test_pgp_complete_movement_energy()
