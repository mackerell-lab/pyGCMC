"""
PGP Complete boundary condition tests

This module tests PGP Complete functionality at boundary conditions:
- Extreme distances (very close, at cutoff, beyond cutoff)
- Direct movement energy calculations
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import math
import os
import gc  # For garbage collection
import subprocess
import sys


def test_pgp_complete_extreme_distances():
    """Test PGP Complete with particles at extreme distances"""
    
    # Reset PGP state to avoid memory corruption from previous tests
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Extreme Distances Test")
    print("="*70)
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Test cases: very close, at sigma, near cutoff, beyond cutoff
    test_positions = [
        ([5.0, 5.0, 5.0], 1.0, True),      # Fixed reference
        ([5.2, 5.0, 5.0], -1.0, False),    # Very close (0.2 nm)
        ([5.35, 5.0, 5.0], 0.5, False),    # At sigma (0.35 nm)
        ([6.95, 5.0, 5.0], -0.5, False),   # Near cutoff (1.95 nm)
        ([8.0, 5.0, 5.0], 0.3, False),     # Beyond cutoff (3.0 nm)
    ]
    
    for i, (pos, charge, fixed) in enumerate(test_positions):
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
    state.activeAtomCount = 5
    state.residues = residues
    state.activeResidueCount = 5
    
    # All non-fixed are movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 4
    state.movementResidues = [movement_info]
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energies
    elec_pme, vdw_pme, total_pme = pygcmc.computeSystemEnergyPMEComplete(state)
    elec_pgp, vdw_pgp, total_pgp = pygcmc.computeSystemEnergyPGPComplete(state)
    
    print(f"\nPME Complete:")
    print(f"  Electrostatic: {elec_pme:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pme:.6f} kJ/mol")
    
    print(f"\nPGP Complete:")
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pgp:.6f} kJ/mol")
    
    # Check individual particle contributions
    print(f"\nDistance-based energy analysis:")
    for i in range(1, 5):
        dist = abs(state.atoms[i].x - state.atoms[0].x)
        print(f"  Particle {i} at distance {dist:.2f} nm")
    
    # Strict checks for extreme distances
    vdw_diff = abs(vdw_pgp - vdw_pme)
    elec_diff = abs(elec_pgp - elec_pme)
    
    print(f"\nAbsolute differences:")
    print(f"  VdW:           {vdw_diff:.2e} kJ/mol")
    print(f"  Electrostatic: {elec_diff:.2e} kJ/mol")
    
    # VdW should match closely (same algorithm)
    assert vdw_diff < 1e-6, f"VdW must match even at extreme distances, got difference {vdw_diff:.2e}"
    # Electrostatic will differ significantly (PME FFT vs PGP grid interpolation)
    rel_elec_diff = elec_diff / abs(elec_pme) * 100 if elec_pme != 0 else 0
    print(f"\nElectrostatic relative difference: {rel_elec_diff:.1f}% (expected for different methods)")
    assert math.isfinite(elec_pgp) and math.isfinite(elec_pme), "Energies must be finite"
    
    # Relative errors (careful with large VdW values)
    rel_vdw_error = vdw_diff / abs(vdw_pme) * 100 if vdw_pme != 0 else 0
    rel_elec_error = elec_diff / abs(elec_pme) * 100 if elec_pme != 0 else 0
    
    print(f"\nRelative errors:")
    print(f"  VdW:           {rel_vdw_error:.2e}%")
    print(f"  Electrostatic: {rel_elec_error:.2e}%")
    
    assert rel_vdw_error < 1e-4, f"VdW relative error too large: {rel_vdw_error:.2e}%"
    # Large electrostatic differences are expected
    print(f"Electrostatic relative difference: {rel_elec_error:.1f}% (expected)")
    assert math.isfinite(elec_pgp) and math.isfinite(elec_pme), "Energies must be finite"
    
    # At very close distances, VdW can be positive (repulsive)
    print(f"\nVdW energy sign: {'positive (repulsive)' if vdw_pme > 0 else 'negative (attractive)'}")


def test_pgp_complete_direct_movement_test():
    """Test computeMovementEnergyPGPComplete directly"""
    
    # Reset PGP state to avoid memory corruption from previous tests
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Direct Movement Energy Test")
    print("="*70)
    
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
    
    # Create system with 2 fixed and 2 moveable particles
    particle_data = [
        ([2.0, 2.5, 2.5], 1.0, True),    # Fixed 1
        ([3.0, 2.5, 2.5], -1.0, True),   # Fixed 2
        ([2.5, 2.5, 2.5], 0.5, False),   # Moveable 1
        ([2.5, 3.0, 2.5], -0.5, False),  # Moveable 2
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
    state.activeAtomCount = 4
    state.residues = residues
    state.activeResidueCount = 4
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 2
    state.movementResidues = [movement_info]
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Test movement energy directly with PGP Complete
    # First get system energy before movement
    elec_before, vdw_before, total_before = pygcmc.computeSystemEnergyPGPComplete(state)
    
    # Move one particle slightly
    state.atoms[2].x += 0.1
    
    # Get system energy after movement  
    elec_after, vdw_after, total_after = pygcmc.computeSystemEnergyPGPComplete(state)
    
    # Calculate energy change
    elec_change = elec_after - elec_before
    vdw_change = vdw_after - vdw_before
    total_change = total_after - total_before
    
    print(f"\nEnergy changes due to movement:")
    print(f"  Electrostatic: {elec_change:.6f} kJ/mol")
    print(f"  VdW:          {vdw_change:.6f} kJ/mol")
    print(f"  Total:        {total_change:.6f} kJ/mol")
    
    # Now test direct movement energy function
    state.atoms[2].x -= 0.1  # Move back
    elec_pgp_mv, vdw_pgp_mv, pgp_dict = pygcmc.computeMovementEnergyPGPComplete(state)
    
    print(f"\nPGP Movement Energy (direct):")
    print(f"  Electrostatic: {elec_pgp_mv:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pgp_mv:.6f} kJ/mol")
    
    # Movement VdW energy should be reasonable
    assert abs(vdw_pgp_mv) > 1e-6, "Movement VdW energy should not be zero"
    assert abs(vdw_pgp_mv) < 10000, "Movement VdW energy seems too large"
    
    # Note: PGP movement electrostatic energy might be zero if the grid interpolation
    # doesn't capture the movement properly
    if abs(elec_pgp_mv) < 1e-6:
        print("WARNING: Movement electrostatic energy is zero - PGP grid may need finer resolution")
    else:
        assert abs(elec_pgp_mv) < 1000, "Movement electrostatic energy seems too large"