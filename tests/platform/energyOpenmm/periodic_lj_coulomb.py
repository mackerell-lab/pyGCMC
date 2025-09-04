# tests/simulation/energyOpenmm/periodic_lj_coulomb.py

import pytest
from .periodic_helpers import *

def test_separate_lj_coulomb_periodic():
    """
    Compare LJ and reaction field Coulomb terms separately in periodic boundary conditions.
    
    Tests individual terms:
    1. LJ term with switching:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] * S(r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Reaction field Coulomb:
       E_rf = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Parameters:
    - Cutoff: 1.0 nm
    - Switching: 0.9 nm
    - Reaction field dielectric: 78.5
    
    Verification:
    - Compares each term with OpenMM reference
    - Analyzes relative contributions
    """
    # Create test system
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    epsilon_rf = 78.5  # relative dielectric constant of water

    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("No NonbondedForce found in system")

    # Define LJ energy expression (with switching function)
    lj_expression = """
    4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * sw;
    epsilon = sqrt(epsilon1*epsilon2);
    sigma = 0.5*(sigma1+sigma2);
    sw = step(cutoff - r) * (step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3) + step(switch - r));
    """

    # Define reaction field Coulomb energy expression
    coulomb_expression = """
    kC * q1 * q2 * (1/r + krf * r^2 - crf);
    krf = (epsilon_rf - 1) / (2*epsilon_rf + 1) / cutoff^3;
    crf = (3*epsilon_rf) / (2*epsilon_rf + 1) / cutoff;
    """

    # Test different distances
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.1]
    print("\nTesting LJ and Coulomb terms separately:")
    print("\nDistance(nm) |  Custom LJ  |  Custom Coulomb  |  Total Custom  |  OpenMM Total  |  Diff(%)")
    print("-" * 85)

    platform = Platform.getPlatformByName('Reference')

    for dist in distances:
        # Move water molecule to new position
        new_positions = []
        for i, pos in enumerate(positions):
            if i >= 6 and i < 9:  # water molecule
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)

        # 1. Calculate LJ energy
        lj_force = CustomNonbondedForce(lj_expression)
        lj_force.addPerParticleParameter("sigma")
        lj_force.addPerParticleParameter("epsilon")
        lj_force.addGlobalParameter("cutoff", cutoff_distance)
        lj_force.addGlobalParameter("switch", switch_distance)

        # Add particle parameters (LJ parameters only)
        for i in range(original_nb_force.getNumParticles()):
            _, sigma, epsilon = original_nb_force.getParticleParameters(i)
            lj_force.addParticle([sigma, epsilon])

        lj_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        lj_force.setCutoffDistance(cutoff_distance * nanometers)

        # Create LJ system
        lj_system = System()
        for i in range(system.getNumParticles()):
            lj_system.addParticle(system.getParticleMass(i))
        lj_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        lj_system.addForce(lj_force)

        # Calculate LJ energy
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(lj_system, integrator, platform)
        context.setPositions(new_positions)
        lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 2. Calculate Coulomb energy
        coulomb_force = CustomNonbondedForce(coulomb_expression)
        coulomb_force.addPerParticleParameter("q")
        coulomb_force.addGlobalParameter("kC", 138.935456)
        coulomb_force.addGlobalParameter("cutoff", cutoff_distance)
        coulomb_force.addGlobalParameter("epsilon_rf", epsilon_rf)

        # Add particle parameters (charges only)
        for i in range(original_nb_force.getNumParticles()):
            charge, _, _ = original_nb_force.getParticleParameters(i)
            coulomb_force.addParticle([charge])

        coulomb_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        coulomb_force.setCutoffDistance(cutoff_distance * nanometers)

        # Create Coulomb system
        coulomb_system = System()
        for i in range(system.getNumParticles()):
            coulomb_system.addParticle(system.getParticleMass(i))
        coulomb_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        coulomb_system.addForce(coulomb_force)

        # Calculate Coulomb energy
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(coulomb_system, integrator, platform)
        context.setPositions(new_positions)
        coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 3. Calculate standard OpenMM energy as reference
        ref_force = NonbondedForce()
        ref_force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        ref_force.setCutoffDistance(cutoff_distance * nanometers)
        ref_force.setUseSwitchingFunction(True)
        ref_force.setSwitchingDistance(switch_distance * nanometers)
        ref_force.setReactionFieldDielectric(epsilon_rf)

        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            ref_force.addParticle(charge, sigma, epsilon)

        ref_system = System()
        for i in range(system.getNumParticles()):
            ref_system.addParticle(system.getParticleMass(i))
        ref_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        ref_system.addForce(ref_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(ref_system, integrator, platform)
        context.setPositions(new_positions)
        ref_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # Convert to kJ/mol and calculate total energy
        lj_val = lj_energy.value_in_unit(kilojoules_per_mole)
        coulomb_val = coulomb_energy.value_in_unit(kilojoules_per_mole)
        custom_total = lj_val + coulomb_val
        ref_val = ref_energy.value_in_unit(kilojoules_per_mole)

        # Calculate relative difference
        abs_diff = abs(custom_total - ref_val)
        rel_diff = abs_diff / abs(ref_val) * 100 if abs(ref_val) > 1e-6 else abs_diff

        print(f"{dist:8.2f} | {lj_val:10.4f} | {coulomb_val:14.4f} | {custom_total:12.4f} | {ref_val:13.4f} | {rel_diff:8.4f}")

        # If difference is large, output detailed information
        if rel_diff > 0.05:  # Output detailed info when difference > 0.05%
            print(f"\n  Detailed information at {dist} nm:")
            print(f"    LJ energy:        {lj_val:.6f} kJ/mol")
            print(f"    Coulomb energy:   {coulomb_val:.6f} kJ/mol")
            print(f"    Custom total:     {custom_total:.6f} kJ/mol")
            print(f"    OpenMM energy:    {ref_val:.6f} kJ/mol")
            print(f"    Absolute diff:    {abs_diff:.6f} kJ/mol")
            print(f"    Relative diff:    {rel_diff:.6f}%")

        # Verify results: use different tolerances based on distance
        if dist <= switch_distance:
            # Use stricter tolerance within switching distance
            assert rel_diff < 0.1, f"Energy difference too large at {dist} nm: {rel_diff:.6f}% > 0.1%"
        elif dist < cutoff_distance:
            # Use looser tolerance in switching region
            assert rel_diff < 0.5, f"Energy difference too large in switching region at {dist} nm: {rel_diff:.6f}% > 0.5%"
        else:
            # Energy should be close to zero beyond cutoff
            assert abs_diff < 2e-1, f"Energy should be close to zero beyond cutoff at {dist} nm, but difference is {abs_diff:.6f} kJ/mol"

