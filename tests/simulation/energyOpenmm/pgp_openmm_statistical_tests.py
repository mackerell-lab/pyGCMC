"""
PGP Complete vs OpenMM statistical tests

This module contains statistical tests comparing PGP Complete with OpenMM PME
for multiple random displacements.
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import numpy as np
import math
import openmm
import openmm.unit as unit


def test_pgp_vs_openmm_multiple_displacements():
    """Test PGP Complete vs OpenMM for multiple random displacements"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    np.random.seed(42)
    
    print("\n" + "="*70)
    print("PGP Complete vs OpenMM Multiple Displacements Test")
    print("="*70)
    
    # Create state
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    sigma_na = 0.333
    sigma_cl = 0.442
    eps_na = 0.0115
    eps_cl = 0.4184
    
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # Create atoms - simple 2 Na-Cl pairs
    atoms = []
    
    # Fixed Na-Cl
    positions = [
        (2.0, 2.5, 2.5, 1.0, 0),   # Na+
        (2.5, 2.5, 2.5, -1.0, 1),  # Cl-
        # Moving Na-Cl
        (3.5, 2.5, 2.5, 1.0, 0),   # Na+
        (4.0, 2.5, 2.5, -1.0, 1),  # Cl-
    ]
    
    for x, y, z, charge, atype in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = x, y, z
        atom.charge = charge
        atom.type = atype
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Store original positions
    original_positions = [(atom.x, atom.y, atom.z) for atom in atoms]
    
    # Residues
    residues = []
    
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = True
    residues.append(res)
    
    res = MCResidue()
    res.atomStart = 2
    res.atomCount = 2
    res.active = True
    res.fixed = False
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Initialize PGP
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Setup OpenMM
    system = openmm.System()
    for _ in range(4):
        system.addParticle(22.99 if _ % 2 == 0 else 35.45)
    
    system.setDefaultPeriodicBoxVectors(
        [5.0, 0, 0] * unit.nanometer,
        [0, 5.0, 0] * unit.nanometer,
        [0, 0, 5.0] * unit.nanometer
    )
    
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(1.2 * unit.nanometer)
    nonbonded.setEwaldErrorTolerance(1e-5)
    
    for atom in atoms:
        charge = atom.charge * unit.elementary_charge
        if atom.type == 0:  # Na
            sigma = 0.333 * unit.nanometer
            epsilon = 0.0115 * unit.kilojoule_per_mole
        else:  # Cl
            sigma = 0.442 * unit.nanometer
            epsilon = 0.4184 * unit.kilojoule_per_mole
        nonbonded.addParticle(charge, sigma, epsilon)
    
    system.addForce(nonbonded)
    
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    
    # Test N random displacements
    n_tests = 20
    errors = []
    rel_errors = []
    
    print(f"\nTesting {n_tests} random displacements...")
    
    for i in range(n_tests):
        # Random displacement
        max_disp = 0.15
        dx = np.random.uniform(-max_disp, max_disp)
        dy = np.random.uniform(-max_disp, max_disp)
        dz = np.random.uniform(-max_disp, max_disp)
        
        # Reset to original positions
        for j in range(4):
            state.atoms[j].x = original_positions[j][0]
            state.atoms[j].y = original_positions[j][1]
            state.atoms[j].z = original_positions[j][2]
        
        # PGP before
        pgp_result1 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_e1 = pgp_result1[0] + pgp_result1[1]
        
        # Move atoms (as rigid body)
        for j in [2, 3]:  # Movement atoms
            state.atoms[j].x += dx
            state.atoms[j].y += dy
            state.atoms[j].z += dz
        
        # PGP after
        pgp_result2 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_e2 = pgp_result2[0] + pgp_result2[1]
        
        pgp_delta = pgp_e2 - pgp_e1
        
        # OpenMM calculation
        # Before positions
        positions_before = [
            [original_positions[j][0], original_positions[j][1], original_positions[j][2]] * unit.nanometer
            for j in range(4)
        ]
        context.setPositions(positions_before)
        
        state_before = context.getState(getEnergy=True)
        total_before = state_before.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        # Zero movement atoms to get fixed energy
        saved_params = []
        for idx in [2, 3]:
            params = nonbonded.getParticleParameters(idx)
            saved_params.append((idx, params))
            nonbonded.setParticleParameters(idx, 0.0, 1.0*unit.nanometer, 0.0*unit.kilojoule_per_mole)
        nonbonded.updateParametersInContext(context)
        
        state_fixed = context.getState(getEnergy=True)
        fixed_energy = state_fixed.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        # Restore parameters
        for idx, params in saved_params:
            nonbonded.setParticleParameters(idx, *params)
        nonbonded.updateParametersInContext(context)
        
        omm_movement_before = total_before - fixed_energy
        
        # After positions
        positions_after = []
        for j in range(4):
            if j >= 2:  # Movement atoms
                positions_after.append([
                    original_positions[j][0] + dx,
                    original_positions[j][1] + dy,
                    original_positions[j][2] + dz
                ] * unit.nanometer)
            else:
                positions_after.append(positions_before[j])
        
        context.setPositions(positions_after)
        
        state_after = context.getState(getEnergy=True)
        total_after = state_after.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        omm_movement_after = total_after - fixed_energy
        omm_delta = omm_movement_after - omm_movement_before
        
        # Record errors
        error = pgp_delta - omm_delta
        errors.append(error)
        if abs(omm_delta) > 0.1:  # Only calculate relative error for significant energy changes
            rel_error = abs(error) / abs(omm_delta) * 100
            rel_errors.append(rel_error)
    
    # Analyze results
    errors = np.array(errors)
    rel_errors = np.array(rel_errors)
    
    print(f"\nAbsolute error statistics:")
    print(f"  Mean: {np.mean(np.abs(errors)):.2e} kJ/mol")
    print(f"  Max: {np.max(np.abs(errors)):.2e} kJ/mol")
    print(f"  95th percentile: {np.percentile(np.abs(errors), 95):.2e} kJ/mol")
    
    if len(rel_errors) > 0:
        print(f"\nRelative error statistics (for ΔE > 0.1 kJ/mol):")
        print(f"  Mean: {np.mean(rel_errors):.2f}%")
        print(f"  Max: {np.max(rel_errors):.2f}%")
        print(f"  95th percentile: {np.percentile(rel_errors, 95):.2f}%")
    
    # Check reasonable accuracy
    assert np.percentile(np.abs(errors), 95) < 0.1, "95% of errors should be < 0.1 kJ/mol"
    
    print("\n✅ Multiple displacements test PASSED")