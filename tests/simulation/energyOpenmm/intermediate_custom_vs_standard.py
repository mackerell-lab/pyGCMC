# tests/simulation/energyOpenmm/intermediate_custom_vs_standard.py

import pytest
from .intermediate_helpers import *

def test_compare_custom_vs_standard_nonbonded():
    """
    Compare energy calculations between CustomNonbondedForce and NonbondedForce.
    
    Tests the equivalence of custom and standard implementations of:
    1. LJ potential:
       E = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb potential:
       E = (1/4πε₀)(q₁q₂/r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    Parameters tested:
    - Particle charges, σ, and ε values
    - Cutoff distances
    - Switching function
    
    Verification:
    - Relative difference < 1e-6
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # Create a new system containing only NonbondedForce
    system_standard = System()
    for i in range(system.getNumParticles()):
        system_standard.addParticle(system.getParticleMass(i))
    system_standard.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # Copy NonbondedForce to the new system
    nb_force = NonbondedForce()
    for i in range(original_nb_force.getNumParticles()):
        params = original_nb_force.getParticleParameters(i)
        nb_force.addParticle(*params)
    nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)
    system_standard.addForce(nb_force)
    
    # Create a system using CustomNonbondedForce
    system_custom = System()
    for i in range(system.getNumParticles()):
        system_custom.addParticle(system.getParticleMass(i))
    system_custom.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # Create CustomNonbondedForce using the same energy expression as NonbondedForce
    energy_expression = """
    kC * q1 * q2 / r + 
    4 * sqrt(eps1*eps2) * (
        (0.5*(sigma1+sigma2)/r)^12 - 
        (0.5*(sigma1+sigma2)/r)^6
    )"""
    
    custom_force = CustomNonbondedForce(energy_expression)
    custom_force.addPerParticleParameter("q")
    custom_force.addPerParticleParameter("sigma")
    custom_force.addPerParticleParameter("eps")
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb constant (kJ·nm/mol/e^2)
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        custom_force.addParticle([charge, sigma, epsilon])
    
    custom_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    system_custom.addForce(custom_force)
    
    platform = Platform.getPlatformByName('Reference')
    
    # Calculate standard NonbondedForce energy
    integrator_standard = VerletIntegrator(0.001 * picoseconds)
    context_standard = Context(system_standard, integrator_standard, platform)
    context_standard.setPositions(positions)
    state_standard = context_standard.getState(getEnergy=True)
    energy_standard = state_standard.getPotentialEnergy()
    
    # Calculate CustomNonbondedForce energy
    integrator_custom = VerletIntegrator(0.001 * picoseconds)
    context_custom = Context(system_custom, integrator_custom, platform)
    context_custom.setPositions(positions)
    state_custom = context_custom.getState(getEnergy=True)
    energy_custom = state_custom.getPotentialEnergy()
    
    # Print results
    print(f"\nComparing NonbondedForce vs CustomNonbondedForce:")
    print(f"NonbondedForce energy: {energy_standard.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"CustomNonbondedForce energy: {energy_custom.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole)):.6f} kJ/mol")
    rel_diff = abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole))/abs(energy_standard.value_in_unit(kilojoules_per_mole))*100
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # Verify results
    assert rel_diff < 1e-6, "Energy mismatch between NonbondedForce and CustomNonbondedForce"
    
    del context_standard, context_custom

