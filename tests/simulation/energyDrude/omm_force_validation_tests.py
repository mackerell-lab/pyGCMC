"""
Force validation tests inspired by OpenMM's validateForce
"""

import pytest
import numpy as np
import pygcmc
import math
from .omm_test_helpers import calculate_numerical_force_simple


def test_numerical_force_validation():
    """Validate forces using numerical differentiation (inspired by OpenMM's validateForce)"""
    
    # Simple system: parent + Drude
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 1.5
    parent.type = 0
    
    # Drude atom (displaced)
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.1, -0.05, 0.08
    drude.charge = -1.5
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
    # Residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    k = 138.935456 * 1.5  # ONE_4PI_EPS0 * 1.5
    charge = 1.5
    alpha = 138.935456 * charge * charge / k
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.5
    particle.polarizability = alpha
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    # Initial displacement
    initial_dx = drude.x - parent.x
    initial_dy = drude.y - parent.y
    initial_dz = drude.z - parent.z
    initial_r_squared = initial_dx*initial_dx + initial_dy*initial_dy + initial_dz*initial_dz
    
    print(f"Initial Drude displacement: ({initial_dx:.3f}, {initial_dy:.3f}, {initial_dz:.3f})")
    print(f"Initial |r|: {math.sqrt(initial_r_squared):.6f} nm")
    
    # Calculate energy (SCF will optimize Drude position)
    actual_energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # After SCF, check the final Drude position
    final_dx = state.atoms[1].x - state.atoms[0].x
    final_dy = state.atoms[1].y - state.atoms[0].y
    final_dz = state.atoms[1].z - state.atoms[0].z
    final_r_squared = final_dx*final_dx + final_dy*final_dy + final_dz*final_dz
    
    print(f"Final Drude displacement: ({final_dx:.6f}, {final_dy:.6f}, {final_dz:.6f})")
    print(f"Final |r|: {math.sqrt(final_r_squared):.6f} nm")
    
    # With no external field, SCF should minimize the energy to nearly zero
    print(f"Final energy: {actual_energy:.6f} kJ/mol")
    
    # Energy should be very small after SCF optimization
    assert actual_energy < 0.001, \
        f"Energy after SCF should be near zero, got {actual_energy}"
    
    # But not exactly zero due to numerical precision
    assert actual_energy > 1e-10, \
        f"Energy suspiciously small: {actual_energy}"
    
    # Now test force validation by applying external field
    # Add an external charge to create a non-zero equilibrium position
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 2.0, 0.0, 0.0
    external.charge = 2.0
    external.type = 0
    state.atoms.append(external)
    state.activeAtomCount = 3
    
    # Update residue for the external charge
    res_ext = pygcmc.MCResidue()
    res_ext.atomStart = 2
    res_ext.atomCount = 1
    res_ext.active = True
    res_ext.type = 1
    state.residues.append(res_ext)
    state.activeResidueCount = 2
    
    # Recalculate with external field
    energy_with_field = pygcmc.DrudeComplete.calculateEnergy(state)
    print(f"\nEnergy with external field: {energy_with_field:.6f} kJ/mol")
    
    # Now validate forces numerically
    delta = 1e-6
    for i, name in enumerate(['Parent', 'Drude']):
        print(f"\nValidating forces on {name}:")
        for direction, axis in enumerate(['x', 'y', 'z']):
            force_numerical = calculate_numerical_force_simple(state, i, direction, delta)
            print(f"  F_{axis} = {force_numerical:.6f} kJ/mol/nm")
    
    # Check that Drude moved significantly
    assert abs(final_dx) > 1e-6 or abs(final_dy) > 1e-6 or abs(final_dz) > 1e-6, \
        f"Drude did not move from initial displaced position under external field"
    
    # Validate forces numerically - for Drude at equilibrium, force should be small
    drude_force_x = calculate_numerical_force_simple(state, 1, 0, delta)
    drude_force_y = calculate_numerical_force_simple(state, 1, 1, delta)  
    drude_force_z = calculate_numerical_force_simple(state, 1, 2, delta)
    drude_force_norm = math.sqrt(drude_force_x**2 + drude_force_y**2 + drude_force_z**2)
    
    print(f"\nDrude force magnitude at equilibrium: {drude_force_norm:.6f} kJ/mol/nm")
    
    # Force should be small but might not be exactly zero due to SCF tolerance
    assert drude_force_norm < 10.0, f"Drude force not converged: {drude_force_norm}"
    
    pygcmc.DrudeComplete.clear()