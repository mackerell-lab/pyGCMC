# tests/simulation/energyOpenmm/periodic_self_energy.py

import pytest
from .periodic_helpers import *

def test_compare_all_methods_with_self_energy():
    """Compare three methods (simple formula, custom cutoff, and standard OpenMM) with self-energy correction.
    
    All methods use the same interaction groups (only calculating interactions between benzene and water molecules).
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # Define interaction groups
    movement_atoms = set(range(6))  # benzene
    fixed_atoms = set(range(6, 9))  # water
    
    # Calculate self-energy correction
    correction = calculate_self_energy_correction(original_nb_force, 1.0)
    print(f"\nSelf-energy correction: {correction:.6f} kJ/mol")
    print("(Note: Standard OpenMM already handles self-energy correction internally)")
    
    # Define test distances
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]
    
    print("\nDetailed Energy Comparison:")
    print("Distance (nm) | Simple (kJ/mol) | Custom (kJ/mol) | Standard (kJ/mol) | Max Diff (kJ/mol)")
    print("-" * 100)
    
    platform = Platform.getPlatformByName('Reference')
    
    for dist in distances:
        # Generate new positions
        new_positions = []
        for i, pos in enumerate(positions):
            if 6 <= i < 9:  # water molecule
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)
        
        # 1. Simple formula (hard cutoff)
        sys_naive = System()
        for i in range(system.getNumParticles()):
            sys_naive.addParticle(system.getParticleMass(i))
        
        naive_expression = """
        step(cutoff - r) * (
            kC * q1 * q2 / r +
            4 * sqrt(eps1*eps2) * (
                (0.5*(sigma1+sigma2)/r)^12 -
                (0.5*(sigma1+sigma2)/r)^6
            )
        )"""
        
        naive_force = CustomNonbondedForce(naive_expression)
        naive_force.addPerParticleParameter("q")
        naive_force.addPerParticleParameter("sigma")
        naive_force.addPerParticleParameter("eps")
        naive_force.addGlobalParameter("kC", 138.935456)
        naive_force.addGlobalParameter("cutoff", 1.0)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            naive_force.addParticle([charge, sigma, epsilon])
        
        naive_force.addInteractionGroup(movement_atoms, fixed_atoms)
        naive_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        naive_force.setCutoffDistance(1.0 * nanometers)
        sys_naive.addForce(naive_force)
        
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_naive, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        simple_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # 2. Custom cutoff formula (with shifting and switching)
        sys_custom = System()
        for i in range(system.getNumParticles()):
            sys_custom.addParticle(system.getParticleMass(i))
        
        cutoff_expression = """
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
        
        custom_force = CustomNonbondedForce(cutoff_expression)
        custom_force.addPerParticleParameter("q")
        custom_force.addPerParticleParameter("sigma")
        custom_force.addPerParticleParameter("eps")
        custom_force.addGlobalParameter("kC", 138.935456)
        custom_force.addGlobalParameter("cutoff", 1.0)
        custom_force.addGlobalParameter("switch", 0.9)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            custom_force.addParticle([charge, sigma, epsilon])
        
        custom_force.addInteractionGroup(movement_atoms, fixed_atoms)
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        custom_force.setCutoffDistance(1.0 * nanometers)
        sys_custom.addForce(custom_force)
        
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_custom, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        custom_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        # No longer apply self-energy correction since standard OpenMM handles it internally
        custom_energy_corrected = custom_energy
        del context, integrator
        
        # 3. Standard OpenMM NonbondedForce
        sys_standard = System()
        for i in range(system.getNumParticles()):
            sys_standard.addParticle(system.getParticleMass(i))
        
        nb_force = NonbondedForce()
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nb_force.addParticle(charge, sigma, epsilon)
        
        # Set interaction groups by adding exceptions
        # Set all non-required interactions to zero
        for i in range(original_nb_force.getNumParticles()):
            for j in range(i+1, original_nb_force.getNumParticles()):
                # If not one atom in movement_atoms and one in fixed_atoms, set as exception
                if not ((i in movement_atoms and j in fixed_atoms) or 
                       (i in fixed_atoms and j in movement_atoms)):
                    nb_force.addException(i, j, 0.0, 1.0, 0.0)
        
        nb_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
        nb_force.setCutoffDistance(1.0 * nanometers)
        nb_force.setUseSwitchingFunction(True)
        nb_force.setSwitchingDistance(0.9 * nanometers)
        sys_standard.addForce(nb_force)
        
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_standard, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        standard_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # Calculate maximum difference
        energies = [simple_energy, custom_energy_corrected, standard_energy]
        max_diff = max([abs(e1 - e2) for e1 in energies for e2 in energies])
        
        print(f"{dist:11.2f} | {simple_energy:13.6f} | {custom_energy_corrected:17.6f} | {standard_energy:15.6f} | {max_diff:16.6f}")
        
        # If difference is too large, output detailed information
        if max_diff > 1.0:  # Output detailed info when difference > 1 kJ/mol
            print(f"  Detailed differences at {dist} nm:")
            print(f"  Custom-Simple: {abs(custom_energy_corrected - simple_energy):.6f} kJ/mol")
            print(f"  Standard-Simple: {abs(standard_energy - simple_energy):.6f} kJ/mol")
            print(f"  Standard-Custom: {abs(standard_energy - custom_energy_corrected):.6f} kJ/mol")

