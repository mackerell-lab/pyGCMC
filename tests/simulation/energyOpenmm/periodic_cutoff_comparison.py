# tests/simulation/energyOpenmm/periodic_cutoff_comparison.py

import pytest
from .periodic_helpers import *

def test_cutoff_periodic_comparison():
    """
    Compare energy between CustomNonbondedForce and NonbondedForce in CutoffPeriodic mode.
    
    Tests periodic boundary implementations with:
    1. Reaction-field electrostatics:
       E = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       where:
       k_rf = (εₛ-1)/(2εₛ+1)/r_c³
       c_rf = (3εₛ)/(2εₛ+1)/r_c
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    2. LJ with switching function:
       S(r) = 1-6x⁵+15x⁴-10x³
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    Parameters:
    - Cutoff: 1.0 nm
    - Switching: 0.9 nm
    - Reaction field dielectric: 78.5
    
    Verification:
    - Compares energies at various distances
    - Analyzes switching function behavior
    """
    # Create test system
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    epsilon_rf = 78.5  # relative dielectric constant of water

    # Output system configuration
    print("\n=== System Configuration ===")
    print(f"Total atoms: {system.getNumParticles()}")
    print(f"Periodic box vectors:")
    vectors = system.getDefaultPeriodicBoxVectors()
    for i, v in enumerate(vectors):
        print(f"  Vector{i+1}: ({v[0]}, {v[1]}, {v[2]}) nm")
    print(f"Cutoff distance: {cutoff_distance} nm")
    print(f"Switching distance: {switch_distance} nm")
    print(f"Reaction field dielectric constant: {epsilon_rf}")

    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("No NonbondedForce found in system")

    # Output particle parameters
    print("\n=== Particle Parameters ===")
    print("Index  Charge(e)  Sigma(nm)  Epsilon(kJ/mol)")
    print("-" * 45)
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        # Convert to dimensionless values
        charge_val = charge.value_in_unit(elementary_charge)
        sigma_val = sigma.value_in_unit(nanometers)
        epsilon_val = epsilon.value_in_unit(kilojoules_per_mole)
        print(f"{i:3d}  {charge_val:8.3f}  {sigma_val:9.3f}  {epsilon_val:13.3f}")

    # Define energy expression
    custom_energy_expression = """
    U_LJ + U_Coulomb;
    U_LJ = 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * sw;
    U_Coulomb = kC * q1 * q2 * (1/r + krf * r^2 - crf);
    sw = step(cutoff - r) * (step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3) + step(switch - r));
    epsilon = sqrt(epsilon1*epsilon2);
    sigma = 0.5*(sigma1+sigma2);
    krf = (epsilon_rf - 1) / (2*epsilon_rf + 1) / cutoff^3;
    crf = (3*epsilon_rf) / (2*epsilon_rf + 1) / cutoff;
    """

    print("\n=== Energy Expression ===")
    print(custom_energy_expression)

    # Calculate and output reaction field parameters
    krf = (epsilon_rf - 1) / (2 * epsilon_rf + 1) / cutoff_distance**3
    crf = (3 * epsilon_rf) / (2 * epsilon_rf + 1) / cutoff_distance
    print("\n=== Reaction Field Parameters ===")
    print(f"krf = {krf:.6f} nm^-3")
    print(f"crf = {crf:.6f} nm^-1")

    # Test different distances
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.1]
    print("\n=== Energy Calculation Results ===")
    print("\nDistance(nm)  Switch Value  CustomNonbondedForce  NonbondedForce    Diff(kJ/mol)  Diff(%)")
    print("-" * 85)

    # Define switching function
    def calc_switching_function(r):
        if r >= cutoff_distance:
            return 0.0
        elif r <= switch_distance:
            return 1.0
        else:
            x = (cutoff_distance - r)**2 * (cutoff_distance + 2*r - 3*switch_distance)
            x /= (cutoff_distance - switch_distance)**3
            return x

    platform = Platform.getPlatformByName('Reference')

    for dist in distances:
        # Calculate switching function value
        switch_value = calc_switching_function(dist)

        # Move water molecule to new position
        new_positions = []
        for i, pos in enumerate(positions):
            if i >= 6 and i < 9:  # water molecule
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)

        # Create new CustomNonbondedForce
        custom_force = CustomNonbondedForce(custom_energy_expression)
        custom_force.addPerParticleParameter("q")
        custom_force.addPerParticleParameter("sigma")
        custom_force.addPerParticleParameter("epsilon")
        custom_force.addGlobalParameter("kC", 138.935456)
        custom_force.addGlobalParameter("cutoff", cutoff_distance)
        custom_force.addGlobalParameter("switch", switch_distance)
        custom_force.addGlobalParameter("epsilon_rf", epsilon_rf)

        # Copy particle parameters from NonbondedForce
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            custom_force.addParticle([charge, sigma, epsilon])

        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        custom_force.setCutoffDistance(cutoff_distance * nanometers)

        # Create CustomNonbondedForce system
        custom_system = System()
        for i in range(system.getNumParticles()):
            custom_system.addParticle(system.getParticleMass(i))
        custom_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        custom_system.addForce(custom_force)

        # Calculate CustomNonbondedForce energy
        custom_integrator = VerletIntegrator(0.001 * picoseconds)
        custom_context = Context(custom_system, custom_integrator, platform)
        custom_context.setPositions(new_positions)
        custom_energy = custom_context.getState(getEnergy=True).getPotentialEnergy()
        
        # Clean up CustomNonbondedForce resources
        del custom_context, custom_integrator, custom_system

        # Create new NonbondedForce
        ref_force = NonbondedForce()
        ref_force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        ref_force.setCutoffDistance(cutoff_distance * nanometers)
        ref_force.setUseSwitchingFunction(True)
        ref_force.setSwitchingDistance(switch_distance * nanometers)
        ref_force.setReactionFieldDielectric(epsilon_rf)

        # Copy particle parameters
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            ref_force.addParticle(charge, sigma, epsilon)

        # Create reference system
        ref_system = System()
        for i in range(system.getNumParticles()):
            ref_system.addParticle(system.getParticleMass(i))
        ref_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        ref_system.addForce(ref_force)

        # Calculate reference energy
        ref_integrator = VerletIntegrator(0.001 * picoseconds)
        ref_context = Context(ref_system, ref_integrator, platform)
        ref_context.setPositions(new_positions)
        ref_energy = ref_context.getState(getEnergy=True).getPotentialEnergy()

        # Clean up reference system resources
        del ref_context, ref_integrator, ref_system

        # Calculate differences
        custom_val = custom_energy.value_in_unit(kilojoules_per_mole)
        ref_val = ref_energy.value_in_unit(kilojoules_per_mole)
        abs_diff = abs(custom_val - ref_val)
        rel_diff = abs_diff / abs(ref_val) * 100 if abs(ref_val) > 1e-6 else abs_diff

        print(f"{dist:6.2f}    {switch_value:10.4f}  {custom_val:16.6f}    {ref_val:12.6f}    {abs_diff:12.6f}  {rel_diff:8.4f}")

        # If difference is large, output detailed information
        if rel_diff > 0.01:  # Lower threshold to 0.01% to get more detailed information
            print(f"\n=== Detailed Analysis at {dist} nm ===")
            print(f"1. Switching Function Analysis:")
            print(f"   - Switch value: {switch_value:.6f}")
            print(f"   - Relative to switch distance: {'Inside' if dist <= switch_distance else 'Switching region' if dist < cutoff_distance else 'Beyond cutoff'}")
            
            print(f"\n2. Energy Analysis:")
            print(f"   - CustomNonbondedForce: {custom_val:.6f} kJ/mol")
            print(f"   - NonbondedForce:       {ref_val:.6f} kJ/mol")
            print(f"   - Absolute difference:   {abs_diff:.6f} kJ/mol")
            print(f"   - Relative difference:   {rel_diff:.6f}%")

            print(f"\n3. Position Analysis:")
            print(f"   - Distance: {dist:.6f} nm")
            print(f"   - Relative to cutoff: {dist/cutoff_distance:.2%}")
            if switch_distance < dist < cutoff_distance:
                print(f"   - Switching progress: {((dist-switch_distance)/(cutoff_distance-switch_distance)):.2%}")

        # Verify results: use different tolerances based on distance
        if dist <= switch_distance:
            # Use stricter tolerance within switching distance
            assert rel_diff < 0.1, f"Energy difference too large at {dist} nm: {rel_diff:.6f}% > 0.1%"
        elif dist < cutoff_distance:
            # Use looser tolerance in switching region
            assert rel_diff < 0.5, f"Energy difference too large in switching region at {dist} nm: {rel_diff:.6f}% > 0.5%"
        else:
            # Energy should be very close to zero beyond cutoff
            assert abs_diff < 2e-1, f"Energy should be close to zero beyond cutoff at {dist} nm, but difference is {abs_diff:.6f} kJ/mol"

    print("\n=== Test Summary ===")
    print("Energy calculations at all distances are within acceptable error ranges")
    print("- Within switching distance (r ≤ 0.9 nm): relative error < 0.1%")
    print("- Switching region (0.9 nm < r < 1.0 nm): relative error < 0.5%")
    print("- Beyond cutoff (r ≥ 1.0 nm): absolute error < 0.2 kJ/mol")

