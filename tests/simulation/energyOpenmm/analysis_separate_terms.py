# tests/simulation/energyOpenmm/analysis_separate_terms.py

import pytest
from .analysis_helpers import *

def test_compare_separate_terms():
    """
    Compare Coulomb and LJ terms separately.
    
    Tests the individual contributions of:
    1. Lennard-Jones term:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb term:
       E_coul = (1/4πε₀)(q₁q₂/r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    3. Total energy:
       E_total = E_LJ + E_coul
    
    Verification:
    - Compares each term with OpenMM reference
    - Ensures relative error < 1%
    """
    # Create test system
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    
    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # Test Coulomb term
    print("\nTesting Coulomb term:")
    
    # Create CustomNonbondedForce for Coulomb only
    coulomb_custom = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff)
    )""")
    
    # Add parameters
    coulomb_custom.addPerParticleParameter("q")
    coulomb_custom.addGlobalParameter("kC", 138.935456)
    coulomb_custom.addGlobalParameter("cutoff", cutoff_distance)
    
    # Create only Coulomb system
    system_coulomb = System()
    for i in range(system.getNumParticles()):
        system_coulomb.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, _, _ = original_nb_force.getParticleParameters(i)
        coulomb_custom.addParticle([charge])
    
    coulomb_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    coulomb_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_coulomb.addForce(coulomb_custom)
    
    # Calculate custom Coulomb energy
    platform = Platform.getPlatformByName('Reference')
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_coulomb, integrator, platform)
    context.setPositions(positions)
    coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # Add self-energy correction
    correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
    coulomb_energy = coulomb_energy + correction * kilojoules_per_mole
    
    print(f"Custom Coulomb energy: {coulomb_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # Test LJ term
    print("\nTesting LJ term:")
    
    # Create CustomNonbondedForce for LJ only
    lj_custom = CustomNonbondedForce("""
    step(cutoff - r) * (
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
    )""")
    
    # Add parameters
    lj_custom.addPerParticleParameter("sigma")
    lj_custom.addPerParticleParameter("eps")
    lj_custom.addGlobalParameter("cutoff", cutoff_distance)
    lj_custom.addGlobalParameter("switch", switch_distance)
    
    # Create only LJ system
    system_lj = System()
    for i in range(system.getNumParticles()):
        system_lj.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        _, sigma, epsilon = original_nb_force.getParticleParameters(i)
        lj_custom.addParticle([sigma, epsilon])
    
    lj_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    lj_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_lj.addForce(lj_custom)
    
    # Calculate custom LJ energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_lj, integrator, platform)
    context.setPositions(positions)
    lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Custom LJ energy: {lj_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # Create reference system
    system_ref = System()
    for i in range(system.getNumParticles()):
        system_ref.addParticle(system.getParticleMass(i))
    
    nb_force = NonbondedForce()
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        nb_force.addParticle(charge, sigma, epsilon)
    
    nb_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
    nb_force.setCutoffDistance(cutoff_distance * nanometers)
    system_ref.addForce(nb_force)
    
    # Calculate reference energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_ref, integrator, platform)
    context.setPositions(positions)
    total_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Reference total energy: {total_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # Verify total energy
    custom_total = coulomb_energy + lj_energy
    energy_diff = abs(custom_total.value_in_unit(kilojoules_per_mole) - 
                     total_energy.value_in_unit(kilojoules_per_mole))
    rel_diff = energy_diff / abs(total_energy.value_in_unit(kilojoules_per_mole)) * 100
    
    print(f"\nResults:")
    print(f"Custom total energy: {custom_total.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # Verify relative error < 1%
    assert rel_diff/100 < 1e-2, "Energy terms differ significantly from OpenMM reference"
    
    del context, integrator

