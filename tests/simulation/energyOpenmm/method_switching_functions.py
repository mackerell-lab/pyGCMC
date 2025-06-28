# tests/simulation/energyOpenmm/method_switching_functions.py

import pytest
from .method_helpers import *

def test_compare_switching_functions():
    """
    Test nonbonded interactions with switching function.
    
    Tests switching function implementation:
    1. LJ switching function:
       S(r) = 1-6x⁵+15x⁴-10x³, x=(r-r_switch)/(r_cutoff-r_switch)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    Compares:
    1. With switching function
    2. Without switching function
    
    Verification:
    - Switching function properly modifies energy
    - Energy smoothly approaches zero at cutoff
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
    
    # Test with switching function
    print("\nTesting with switching function:")
    
    # Create custom force with switching function
    switched_force = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
    )""")
    
    # Create unswitched force
    unswitched_force = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )""")
    
    # Add parameters
    switched_force.addPerParticleParameter("q")
    switched_force.addPerParticleParameter("sigma")
    switched_force.addPerParticleParameter("eps")
    switched_force.addGlobalParameter("kC", 138.935456)
    switched_force.addGlobalParameter("cutoff", cutoff_distance)
    switched_force.addGlobalParameter("switch", switch_distance)
    
    unswitched_force.addPerParticleParameter("q")
    unswitched_force.addPerParticleParameter("sigma")
    unswitched_force.addPerParticleParameter("eps")
    unswitched_force.addGlobalParameter("kC", 138.935456)
    unswitched_force.addGlobalParameter("cutoff", cutoff_distance)
    
    # Create system
    system_switched = System()
    for i in range(system.getNumParticles()):
        system_switched.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        switched_force.addParticle([charge, sigma, epsilon])
        unswitched_force.addParticle([charge, sigma, epsilon])
    
    switched_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    switched_force.setCutoffDistance(cutoff_distance * nanometers)
    system_switched.addForce(switched_force)
    
    # Calculate switched energy
    platform = Platform.getPlatformByName('Reference')
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_switched, integrator, platform)
    context.setPositions(positions)
    switched_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # Add self-energy correction
    correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
    switched_energy = switched_energy + correction * kilojoules_per_mole
    
    print(f"Energy with switching: {switched_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # Test without switching
    print("\nTesting without switching:")
    
    # Create unswitched force
    unswitched_force = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )""")
    
    # Add parameters
    unswitched_force.addPerParticleParameter("q")
    unswitched_force.addPerParticleParameter("sigma")
    unswitched_force.addPerParticleParameter("eps")
    unswitched_force.addGlobalParameter("kC", 138.935456)
    unswitched_force.addGlobalParameter("cutoff", cutoff_distance)
    
    # Create system
    system_unswitched = System()
    for i in range(system.getNumParticles()):
        system_unswitched.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        unswitched_force.addParticle([charge, sigma, epsilon])
    
    unswitched_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    unswitched_force.setCutoffDistance(cutoff_distance * nanometers)
    system_unswitched.addForce(unswitched_force)
    
    # Calculate unswitched energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_unswitched, integrator, platform)
    context.setPositions(positions)
    unswitched_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # Add self-energy correction
    unswitched_energy = unswitched_energy + correction * kilojoules_per_mole
    
    print(f"Energy without switching: {unswitched_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
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
    nb_force.setUseSwitchingFunction(True)
    nb_force.setSwitchingDistance(switch_distance * nanometers)
    system_ref.addForce(nb_force)
    
    # Calculate reference energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_ref, integrator, platform)
    context.setPositions(positions)
    ref_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Reference energy: {ref_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # Verify results
    print("\nResults:")
    
    # Verify switched energy vs reference
    energy_diff = abs(switched_energy.value_in_unit(kilojoules_per_mole) - 
                     ref_energy.value_in_unit(kilojoules_per_mole))
    rel_diff = energy_diff / abs(ref_energy.value_in_unit(kilojoules_per_mole)) * 100
    
    print(f"Difference from reference (with switching):")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # Verify relative error < 0.1%
    assert rel_diff/100 < 1e-3, "Energy with switching differs significantly from OpenMM reference"
    
    # Verify switching function effect
    energy_diff = abs(switched_energy.value_in_unit(kilojoules_per_mole) - 
                     unswitched_energy.value_in_unit(kilojoules_per_mole))
    print(f"\nDifference between switched and unswitched:")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    
    # Verify switching function affects energy
    assert energy_diff > 0, "Switching function should affect the energy"
    
    del context, integrator

