#!/usr/bin/env python3
"""
Improved test documenting the optimization difference between PyGCMC and OpenMM.
Uses dynamic calculations instead of hardcoded values.
"""

import pytest
import numpy as np
import pygcmc
import logging

# Configure logging
logger = logging.getLogger(__name__)

# Try importing OpenMM
try:
    import openmm
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False


def calculate_analytical_drude_displacement(q_drude, q_external, r0, alpha):
    """
    Calculate analytical displacement for a Drude oscillator in external field.
    
    This is the exact solution for force equilibrium (F=0).
    """
    k_coulomb = 138.935456  # kJ·nm/mol/e²
    k_spring = k_coulomb * q_drude**2 / alpha
    
    # Electric field at origin from external charge
    E_field = k_coulomb * q_external / r0**2
    
    # Force equilibrium: k_spring * x = q_drude * E_field
    displacement = q_drude * E_field / k_spring
    
    return displacement


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_drude_optimization_difference_dynamic():
    """
    Test and document the difference between PyGCMC and OpenMM optimization.
    Uses dynamic calculations instead of hardcoded values.
    """
    logger.info("=== Test: PyGCMC vs OpenMM Optimization Difference ===")
    
    # Physical parameters
    q_external = 1.0    # e
    q_drude = -1.0      # e
    alpha = 0.001       # nm³
    r0 = 0.5            # nm
    
    # Calculate theoretical displacement
    theoretical_disp = calculate_analytical_drude_displacement(
        q_drude, q_external, r0, alpha
    )
    logger.info(f"Theoretical displacement (F=0): {theoretical_disp:.8f} nm")
    
    # 1. PyGCMC calculation
    state = create_drude_system(q_drude, q_external, r0, alpha)
    energy_pygcmc = pygcmc.DrudeComplete.calculateEnergy(state)
    x_pygcmc = state.atoms[1].x
    
    logger.info(f"PyGCMC displacement: {x_pygcmc:.8f} nm")
    logger.info(f"PyGCMC vs theory: {100*(x_pygcmc/theoretical_disp - 1):.2f}%")
    
    # 2. OpenMM calculation
    x_openmm, energy_openmm = run_openmm_calculation(
        q_drude, q_external, r0, alpha
    )
    
    logger.info(f"OpenMM displacement: {x_openmm:.8f} nm")
    logger.info(f"OpenMM vs theory: {100*(x_openmm/theoretical_disp - 1):.2f}%")
    
    # 3. Analyze the difference
    rel_diff = abs(x_openmm / x_pygcmc - 1)
    logger.info(f"Relative difference: {100*rel_diff:.2f}%")
    
    # Dynamic assertions based on physics
    # PyGCMC should be very close to theoretical (F=0)
    assert abs(x_pygcmc / theoretical_disp - 1) < 0.001, \
        "PyGCMC should match theoretical F=0 solution"
    
    # OpenMM and PyGCMC should differ by less than 1%
    assert rel_diff < 0.01, \
        f"PyGCMC and OpenMM differ by {100*rel_diff:.1f}%, expected < 1%"
    
    # Calculate residual force at OpenMM position in PyGCMC framework
    residual_force = calculate_residual_force(
        state, x_openmm, q_drude, alpha
    )
    logger.info(f"Residual force at OpenMM position: {residual_force:.6f} kJ/mol/nm")
    
    # The residual force should be small but non-zero
    assert abs(residual_force) > 1e-3, \
        "OpenMM position should have non-zero force (confirms different objectives)"
    assert abs(residual_force) < 10.0, \
        "Residual force should still be reasonably small"


def create_drude_system(q_drude, q_external, r0, alpha):
    """Create a simple Drude system for testing"""
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = q_drude
    drude.type = 1
    
    # External charge
    external = pygcmc.MCAtom()
    external.x = r0
    external.y = external.z = 0.0
    external.charge = q_external
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    
    # Setup residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = q_drude
    particle.polarizability = alpha
    particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
    particle.aniso12 = particle.aniso34 = 1.0
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 1000
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    return state


def run_openmm_calculation(q_drude, q_external, r0, alpha):
    """Run OpenMM calculation and return displacement and energy"""
    # Create OpenMM system
    system = openmm.System()
    system.addParticle(1.0 * unit.dalton)  # Parent
    system.addParticle(0.4 * unit.dalton)  # Drude
    system.addParticle(1.0 * unit.dalton)  # External
    
    # Drude force
    drude_force = openmm.DrudeForce()
    drude_force.addParticle(
        1, 0, -1, -1, -1,
        q_drude,
        alpha,
        1.0, 1.0
    )
    system.addForce(drude_force)
    
    # Nonbonded force
    nonbonded = openmm.NonbondedForce()
    nonbonded.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
    nonbonded.addParticle(q_drude * unit.elementary_charge, 0.1 * unit.nanometer, 0.0)
    nonbonded.addParticle(q_external * unit.elementary_charge, 0.1 * unit.nanometer, 0.0)
    
    # Exclude parent-drude interaction
    nonbonded.addException(0, 1, 0.0, 1.0, 0.0)
    system.addForce(nonbonded)
    
    # Create context
    integrator = openmm.DrudeSCFIntegrator(0.001 * unit.picoseconds)
    integrator.setMinimizationErrorTolerance(1e-8)
    
    context = openmm.Context(system, integrator)
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [r0, 0.0, 0.0]
    ] * unit.nanometer
    context.setPositions(positions)
    
    # Minimize
    integrator.step(1)
    
    # Get results
    state = context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    x_openmm = pos[1][0] - pos[0][0]
    
    return x_openmm, energy


def calculate_residual_force(state, x_test, q_drude, alpha):
    """Calculate the residual force at a test position"""
    k_coulomb = 138.935456
    k_spring = k_coulomb * q_drude**2 / alpha
    
    # Set Drude to test position
    state.atoms[1].x = x_test
    
    # Calculate forces
    # Spring force
    f_spring = -k_spring * x_test
    
    # Coulomb force from external charge
    r = state.atoms[2].x - x_test
    f_coulomb = k_coulomb * q_drude * state.atoms[2].charge / r**2
    
    # Total force
    f_total = f_spring + f_coulomb
    
    return f_total


@pytest.mark.parametrize("field_strength", [0.1, 1.0, 10.0])
def test_optimization_difference_vs_field_strength(field_strength):
    """
    Test how the PyGCMC-OpenMM difference scales with field strength.
    This verifies the difference is systematic, not random.
    """
    if not HAS_OPENMM:
        pytest.skip("OpenMM not available")
    
    # Base parameters
    q_drude = -1.0
    alpha = 0.001
    r0 = 0.5
    
    # Vary external charge to change field strength
    q_external = field_strength
    
    # Calculate displacements
    state = create_drude_system(q_drude, q_external, r0, alpha)
    energy_pygcmc = pygcmc.DrudeComplete.calculateEnergy(state)
    x_pygcmc = state.atoms[1].x
    
    x_openmm, energy_openmm = run_openmm_calculation(
        q_drude, q_external, r0, alpha
    )
    
    # The relative difference should be roughly constant
    rel_diff = abs(x_openmm / x_pygcmc - 1)
    
    logger.info(f"Field strength {field_strength}: rel diff = {100*rel_diff:.2f}%")
    
    # All field strengths should show similar relative difference
    assert 0.003 < rel_diff < 0.007, \
        f"Relative difference {100*rel_diff:.2f}% outside expected range"


if __name__ == "__main__":
    test_drude_optimization_difference_dynamic()
    if HAS_OPENMM:
        for strength in [0.1, 1.0, 10.0]:
            test_optimization_difference_vs_field_strength(strength)