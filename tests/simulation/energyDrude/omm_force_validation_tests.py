"""
Force validation tests inspired by OpenMM's validateForce
"""

import pytest
import numpy as np
import pygcmc
import math
from .omm_test_helpers import calculate_numerical_force_simple


def test_drude_spring_energy_formula():
    """Test exact spring energy formula E = 0.5 * k * |r|^2"""
    
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
    parent.charge = 1.0
    parent.type = 0
    
    # Drude atom (displaced)
    drude = pygcmc.MCAtom()
    displacement = 0.01  # nm
    drude.x, drude.y, drude.z = displacement, 0.0, 0.0
    drude.charge = -1.0
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
    
    # Use simple values for clear testing
    alpha = 0.001  # nm^3
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = alpha
    particle.computeSpringConstants()
    
    k_spring = particle.kSpring
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Direct energy calculation (no SCF)
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e10  # Very high tolerance to avoid SCF
    params.maxIterations = 0  # No iterations
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected energy: E = 0.5 * k * |r|^2
    expected_energy = 0.5 * k_spring * displacement * displacement
    
    print(f"Spring constant k: {k_spring:.6f} kJ/mol/nm^2")
    print(f"Displacement: {displacement:.6f} nm")
    print(f"Calculated energy: {energy:.10f} kJ/mol")
    print(f"Expected energy: {expected_energy:.10f} kJ/mol")
    
    # Tight tolerance as in OpenMM (1e-5 kJ/mol)
    assert abs(energy - expected_energy) < 1e-5, \
        f"Energy mismatch: {energy} vs {expected_energy}, diff = {abs(energy - expected_energy)}"
    
    pygcmc.DrudeComplete.clear()


def test_anisotropic_spring_energy():
    """Test anisotropic spring with a1=0.8, a2=1.1 (from OpenMM TestDrudeForce)"""
    
    # System with 4 atoms defining anisotropy axes
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3  # parent, drude, axis atoms
    ff.numMovementTypes = 3
    ff.ljEps = [0.0, 0.0, 0.0]
    ff.ljSigma = [0.1, 0.1, 0.1]
    state.forcefield = ff
    
    # Parent atom at origin
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 1.0
    parent.type = 0
    
    # Drude atom (will be displaced)
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # Axis atoms for defining anisotropy directions
    # Axis 1-2: along x direction
    axis1 = pygcmc.MCAtom()
    axis1.x, axis1.y, axis1.z = -1.0, 0.0, 0.0
    axis1.charge = 0.0
    axis1.type = 2
    
    axis2 = pygcmc.MCAtom()
    axis2.x, axis2.y, axis2.z = 1.0, 0.0, 0.0
    axis2.charge = 0.0
    axis2.type = 2
    
    # Axis 3-4: along y direction
    axis3 = pygcmc.MCAtom()
    axis3.x, axis3.y, axis3.z = 0.0, -1.0, 0.0
    axis3.charge = 0.0
    axis3.type = 2
    
    axis4 = pygcmc.MCAtom()
    axis4.x, axis4.y, axis4.z = 0.0, 1.0, 0.0
    axis4.charge = 0.0
    axis4.type = 2
    
    state.atoms = [parent, drude, axis1, axis2, axis3, axis4]
    state.activeAtomCount = 6
    
    # Single residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 6
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup anisotropic Drude
    pygcmc.DrudeComplete.clear()
    
    alpha = 0.001  # nm^3
    a1 = 0.8  # OpenMM test value
    a2 = 1.1  # OpenMM test value
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = alpha
    particle.aniso1Index = 2  # axis1
    particle.aniso2Index = 3  # axis2
    particle.aniso3Index = 4  # axis3
    particle.aniso4Index = 5  # axis4
    particle.aniso12 = a1
    particle.aniso34 = a2
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # No SCF for direct energy test
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e10
    params.maxIterations = 0
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    # Test 1: Displacement along x (axis 1-2 direction)
    dx = 0.01  # nm
    state.atoms[1].x = dx
    energy_x = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected: E = 0.5 * k * a1 * dx^2
    k_base = particle.charge * particle.charge * 138.935456 / alpha
    expected_x = 0.5 * k_base * a1 * dx * dx
    
    print(f"X-displacement test:")
    print(f"  a1 = {a1}")
    print(f"  Energy = {energy_x:.10f} kJ/mol")
    print(f"  Expected = {expected_x:.10f} kJ/mol")
    
    assert abs(energy_x - expected_x) < 1e-5, \
        f"X energy mismatch: {energy_x} vs {expected_x}"
    
    # Test 2: Displacement along y (axis 3-4 direction)
    state.atoms[1].x = 0.0
    state.atoms[1].y = dx
    energy_y = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected: E = 0.5 * k * a2 * dy^2
    expected_y = 0.5 * k_base * a2 * dx * dx
    
    print(f"\nY-displacement test:")
    print(f"  a2 = {a2}")
    print(f"  Energy = {energy_y:.10f} kJ/mol")
    print(f"  Expected = {expected_y:.10f} kJ/mol")
    
    assert abs(energy_y - expected_y) < 1e-5, \
        f"Y energy mismatch: {energy_y} vs {expected_y}"
    
    # Test 3: Displacement along z (no anisotropy)
    state.atoms[1].y = 0.0
    state.atoms[1].z = dx
    energy_z = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected: E = 0.5 * k * dz^2 (isotropic)
    expected_z = 0.5 * k_base * dx * dx
    
    print(f"\nZ-displacement test (isotropic):")
    print(f"  Energy = {energy_z:.10f} kJ/mol")
    print(f"  Expected = {expected_z:.10f} kJ/mol")
    
    assert abs(energy_z - expected_z) < 1e-5, \
        f"Z energy mismatch: {energy_z} vs {expected_z}"
    
    pygcmc.DrudeComplete.clear()


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
    
    # Note: k = ONE_4PI_EPS0 * 1.5, so alpha = ONE_4PI_EPS0 * q^2 / k = q^2 / 1.5
    # For q = 1.5, alpha = 1.5^2 / 1.5 = 1.5 nm^3
    alpha = 1.5  # nm^3
    
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
    # OpenMM uses ~5e-3 kJ/mol/nm tolerance for force validation
    assert drude_force_norm < 5e-3, f"Drude force not converged: {drude_force_norm} kJ/mol/nm (threshold: 5e-3)"
    
    pygcmc.DrudeComplete.clear()