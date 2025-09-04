# tests/simulation/energyOpenmm/method_combined.py

import pytest
from .method_helpers import *

def test_compare_nonbonded_methods():
    """
    Compare different nonbonded methods between CustomNonbondedForce and NonbondedForce.
    
    Tests different calculation methods:
    1. No cutoff:
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    2. Cutoff with reaction field:
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Tolerances:
    - NoCutoff: 1e-6 (higher precision)
    - CutoffNonPeriodic: 5e-4 (allows for switching function effects)
    
    Verification:
    - Compares energies between custom and standard implementations
    - Ensures differences are within method-specific tolerances
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Test different nonbonded methods
    methods = [
        (NonbondedForce.NoCutoff, CustomNonbondedForce.NoCutoff, "NoCutoff"),
        (NonbondedForce.CutoffNonPeriodic, CustomNonbondedForce.CutoffNonPeriodic, "CutoffNonPeriodic")
    ]
    
    # Define different tolerances for different methods
    tolerances = {
        "NoCutoff": 1e-6,           # Higher precision required for no cutoff
        "CutoffNonPeriodic": 5e-4   # Allow 0.05% error when using cutoff
    }
    
    platform = Platform.getPlatformByName('Reference')
    
    for nb_method, custom_method, method_name in methods:
        print(f"\nTesting {method_name}:")
        
        # Create standard NonbondedForce system
        system_standard = System()
        for i in range(system.getNumParticles()):
            system_standard.addParticle(system.getParticleMass(i))
        system_standard.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        
        nb_force = NonbondedForce()
        original_nb_force = None
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                original_nb_force = force
                break
        
        for i in range(original_nb_force.getNumParticles()):
            params = original_nb_force.getParticleParameters(i)
            nb_force.addParticle(*params)
        
        nb_force.setNonbondedMethod(nb_method)
        cutoff_distance = 1.0  # nm
        switch_distance = 0.9  # nm, distance where switching function starts
        if nb_method != NonbondedForce.NoCutoff:
            nb_force.setCutoffDistance(cutoff_distance * nanometers)
            nb_force.setUseSwitchingFunction(True)
            nb_force.setSwitchingDistance(switch_distance * nanometers)
        system_standard.addForce(nb_force)
        
        # Create CustomNonbondedForce system
        system_custom = System()
        for i in range(system.getNumParticles()):
            system_custom.addParticle(system.getParticleMass(i))
        system_custom.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        
        # Modify energy expression
        if custom_method == CustomNonbondedForce.NoCutoff:
            energy_expression = """
            kC * q1 * q2 / r + 
            4 * sqrt(eps1*eps2) * (
                (0.5*(sigma1+sigma2)/r)^12 - 
                (0.5*(sigma1+sigma2)/r)^6
            )"""
        else:
            energy_expression = """
            step(cutoff - r) * (
                kC * q1 * q2 * (1/r - 1/cutoff) + 
                4 * sqrt(eps1*eps2) * (
                    (0.5*(sigma1+sigma2)/r)^12 - 
                    (0.5*(sigma1+sigma2)/r)^6
                ) * (
                    step(switch - r) +
                    step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
                )
            )"""
        
        custom_force = CustomNonbondedForce(energy_expression)
        custom_force.addPerParticleParameter("q")
        custom_force.addPerParticleParameter("sigma")
        custom_force.addPerParticleParameter("eps")
        custom_force.addGlobalParameter("kC", 138.935456)
        custom_force.addGlobalParameter("cutoff", cutoff_distance)
        custom_force.addGlobalParameter("switch", switch_distance)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            custom_force.addParticle([charge, sigma, epsilon])
        
        custom_force.setNonbondedMethod(custom_method)
        if custom_method != CustomNonbondedForce.NoCutoff:
            custom_force.setCutoffDistance(cutoff_distance * nanometers)
        system_custom.addForce(custom_force)
        
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
        
        # If using cutoff, subtract self-energy correction
        if custom_method != CustomNonbondedForce.NoCutoff:
            correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
            energy_custom = energy_custom - correction * kilojoules_per_mole
        
        # Calculate difference from reference energy
        energy_diff = abs(energy_standard.value_in_unit(kilojoules_per_mole) - 
                         energy_custom.value_in_unit(kilojoules_per_mole))
        rel_diff = energy_diff / abs(energy_standard.value_in_unit(kilojoules_per_mole)) * 100
        
        print(f"NonbondedForce energy: {energy_standard.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
        print(f"CustomNonbondedForce energy: {energy_custom.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
        print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
        print(f"Relative difference: {rel_diff:.6f}%")
        
        # Use corresponding tolerance for verification
        rel_tol = tolerances[method_name]
        print(f"Using relative tolerance: {rel_tol:.6e}")
        
        assert rel_diff/100 < rel_tol, \
               f"Energy mismatch for {method_name}"
        
        del context_standard, context_custom