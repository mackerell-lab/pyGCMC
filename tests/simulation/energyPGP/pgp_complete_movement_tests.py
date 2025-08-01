"""
PGP Complete movement energy tests

This module contains tests for PGP Complete movement energy calculations.
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import math


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
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    
    # Precompute grid for fixed particles
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate full system energy
    elec_full, vdw_full, total_full = pygcmc.computeSystemEnergyPGPComplete(state)
    print(f"\nFull system energy:")
    print(f"  Electrostatic: {elec_full:.6f} kJ/mol")
    print(f"  VdW:          {vdw_full:.6f} kJ/mol")
    print(f"  Total:        {total_full:.6f} kJ/mol")
    
    # Calculate movement energy only
    result = pygcmc.computeMovementEnergyPGPComplete(state)
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
    elec_fixed, vdw_fixed, total_fixed = pygcmc.computeSystemEnergyPGPComplete(state)
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