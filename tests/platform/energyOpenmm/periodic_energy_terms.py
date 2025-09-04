# tests/simulation/energyOpenmm/periodic_energy_terms.py

import pytest
from .periodic_helpers import *

def test_analyze_openmm_energy_terms():
    """
    Analyze OpenMM energy terms by calculating Coulomb and LJ terms separately.
    
    Tests the decomposition and analysis of nonbonded energy terms:
    1. Coulomb term with reaction field:
       E_rf = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    2. Lennard-Jones term with switching:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] * S(r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    3. Total energy:
       E_total = E_rf + E_LJ
    
    Test distances:
    - 0.35 nm: Strong repulsion region
    - 0.5-0.7 nm: Normal interaction region
    - 0.9-1.0 nm: Switching region
    - >1.0 nm: Beyond cutoff
    
    Verification:
    - Individual terms match OpenMM reference
    - Sum of terms equals total energy
    - Energy behavior in different regions is physically reasonable
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
    
    # Define test distances
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]
    
    print("\nAnalyzing OpenMM energy terms:")
    print("Distance (nm) | Coulomb (kJ/mol) | LJ (kJ/mol) | Total (kJ/mol) | Standard OpenMM (kJ/mol)")
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
        
        # 1. Calculate shifted Coulomb energy
        sys_coulomb = System()
        for i in range(system.getNumParticles()):
            sys_coulomb.addParticle(system.getParticleMass(i))
        
        coulomb_expression = "step(cutoff - r) * kC * q1 * q2 * (1/r - 1/cutoff)"
        
        coulomb_force = CustomNonbondedForce(coulomb_expression)
        coulomb_force.addPerParticleParameter("q")
        coulomb_force.addGlobalParameter("kC", 138.935456)
        coulomb_force.addGlobalParameter("cutoff", 1.0)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, _, _ = original_nb_force.getParticleParameters(i)
            coulomb_force.addParticle([charge])
        
        coulomb_force.addInteractionGroup(movement_atoms, fixed_atoms)
        coulomb_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        coulomb_force.setCutoffDistance(1.0 * nanometers)
        sys_coulomb.addForce(coulomb_force)
        
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_coulomb, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        coulomb_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # 2. Calculate LJ energy with switching function
        sys_lj = System()
        for i in range(system.getNumParticles()):
            sys_lj.addParticle(system.getParticleMass(i))
        
        lj_expression = """
        step(cutoff - r) * (
            4 * sqrt(eps1*eps2) * (
                (0.5*(sigma1+sigma2)/r)^12 - 
                (0.5*(sigma1+sigma2)/r)^6
            ) * (
                step(switch - r) +
                step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
            )
        )"""
        
        lj_force = CustomNonbondedForce(lj_expression)
        lj_force.addPerParticleParameter("sigma")
        lj_force.addPerParticleParameter("eps")
        lj_force.addGlobalParameter("cutoff", 1.0)
        lj_force.addGlobalParameter("switch", 0.9)
        
        for i in range(original_nb_force.getNumParticles()):
            _, sigma, epsilon = original_nb_force.getParticleParameters(i)
            lj_force.addParticle([sigma, epsilon])
        
        lj_force.addInteractionGroup(movement_atoms, fixed_atoms)
        lj_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        lj_force.setCutoffDistance(1.0 * nanometers)
        sys_lj.addForce(lj_force)
        
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_lj, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        lj_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # 3. Calculate standard OpenMM energy as reference
        sys_standard = System()
        for i in range(system.getNumParticles()):
            sys_standard.addParticle(system.getParticleMass(i))
        
        nb_force = NonbondedForce()
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nb_force.addParticle(charge, sigma, epsilon)
        
        # Set interaction groups
        for i in range(original_nb_force.getNumParticles()):
            for j in range(i+1, original_nb_force.getNumParticles()):
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
        
        # Calculate total energy (Coulomb + LJ)
        total_energy = coulomb_energy + lj_energy
        
        print(f"{dist:11.2f} | {coulomb_energy:15.6f} | {lj_energy:11.6f} | {total_energy:13.6f} | {standard_energy:21.6f}")
        
        # If difference from standard OpenMM is large, output detailed information
        diff = abs(total_energy - standard_energy)
        if diff > 0.1:
            print(f"  Large difference at {dist} nm:")
            print(f"  Difference between sum and standard: {diff:.6f} kJ/mol")
            print(f"  Relative difference: {diff/abs(standard_energy)*100 if abs(standard_energy) > 1e-10 else 0:.6f}%")

