"""
PGP Complete pure LJ system tests

This module tests PGP Complete functionality for pure LJ systems (no electrostatics)
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import math


def test_pgp_complete_pure_lj():
    """Test PGP Complete with pure LJ system"""
    
    # Reset PGP state to avoid memory corruption from previous tests
    pygcmc.resetPGPState()
    
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
    
    # Electrostatic should be exactly zero (machine precision)
    assert abs(elec_pme) < 1e-14, f"PME electrostatic must be zero to machine precision, got {elec_pme:.2e}"
    assert abs(elec_pgp) < 1e-14, f"PGP electrostatic must be zero to machine precision, got {elec_pgp:.2e}"
    
    # VdW should match to machine precision
    vdw_diff = abs(vdw_pgp - vdw_pme)
    rel_vdw_diff = vdw_diff / abs(vdw_pme) * 100 if vdw_pme != 0 else 0
    print(f"\nVdW absolute difference: {vdw_diff:.2e} kJ/mol")
    print(f"VdW relative difference: {rel_vdw_diff:.2e}%")
    
    assert vdw_diff < 1e-7, f"VdW energies must match to high precision, got difference {vdw_diff:.2e}"
    assert rel_vdw_diff < 1e-5, f"VdW relative error too large: {rel_vdw_diff:.2e}%"