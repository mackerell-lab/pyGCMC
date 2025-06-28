# tests/simulation/energyOpenmm/periodic_force_parameters.py

import pytest
from .periodic_helpers import *

def test_compare_force_parameters_and_energies():
    """
    Compare parameter settings and energy calculations between NonbondedForce and CustomNonbondedForce.
    
    Tests parameter consistency:
    1. NonbondedForce parameters:
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#nonbondedforce
       - Cutoff distance
       - Switching function
       - Reaction field dielectric
    
    2. CustomNonbondedForce implementation:
       - Matching parameters
       - Equivalent energy expressions
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Verification:
    - Parameter values match
    - Energy calculations agree within tolerances:
      * Within switching: < 0.1%
      * Switching region: < 0.5%
      * Beyond cutoff: < 0.2 kJ/mol
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Get original NonbondedForce
    ref_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            ref_force = force
            break
    if ref_force is None:
        raise ValueError("No NonbondedForce found in system")
    
    # Set NonbondedMethod to CutoffPeriodic
    ref_force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    
    # Get current parameter values
    cutoff = ref_force.getCutoffDistance().value_in_unit(nanometers)
    use_switching = ref_force.getUseSwitchingFunction()
    switching_distance = ref_force.getSwitchingDistance().value_in_unit(nanometers)
    rf_dielectric = ref_force.getReactionFieldDielectric()
    
    print("\n=== NonbondedForce Parameter Settings ===")
    print(f"Cutoff distance: {cutoff} nm")
    print(f"Using switching function: {use_switching}")
    print(f"Switching distance: {switching_distance} nm")
    print(f"Reaction field dielectric: {rf_dielectric}")
    
    # Create CustomNonbondedForce
    custom_expression = """
    U_LJ + U_Coulomb;
    U_LJ = 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * sw;
    U_Coulomb = kC * q1 * q2 * (1/r + krf * r^2 - crf);
    epsilon = sqrt(epsilon1*epsilon2);
    sigma = 0.5*(sigma1+sigma2);
    sw = step(cutoff - r) * (step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3) + step(switch - r));
    krf = (epsilon_rf - 1) / (2*epsilon_rf + 1) / cutoff^3;
    crf = (3*epsilon_rf) / (2*epsilon_rf + 1) / cutoff;
    """
    
    custom_force = CustomNonbondedForce(custom_expression)
    custom_force.addPerParticleParameter("q")
    custom_force.addPerParticleParameter("sigma")
    custom_force.addPerParticleParameter("epsilon")
    custom_force.addGlobalParameter("kC", 138.935456)
    custom_force.addGlobalParameter("cutoff", cutoff)
    custom_force.addGlobalParameter("switch", switching_distance)
    custom_force.addGlobalParameter("epsilon_rf", rf_dielectric)
    
    # Add particle parameters
    for i in range(ref_force.getNumParticles()):
        charge, sigma, epsilon = ref_force.getParticleParameters(i)
        custom_force.addParticle([charge, sigma, epsilon])
    
    custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_force.setCutoffDistance(cutoff * nanometers)
    
    # Create system and add force
    custom_system = System()
    
    for i in range(system.getNumParticles()):
        custom_system.addParticle(system.getParticleMass(i))
    
    custom_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    custom_system.addForce(custom_force)
    
    # Define test distances
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.1]
    
    print("\n=== Energy Comparison ===")
    print("Distance(nm)  NonbondedForce(kJ/mol)  CustomNonbondedForce(kJ/mol)  Relative Diff(%)")
    print("-" * 75)
    
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
        
        # Calculate NonbondedForce energy
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(system, integrator, platform)
        context.setPositions(new_positions)
        ref_energy = context.getState(getEnergy=True).getPotentialEnergy()
        ref_val = ref_energy.value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # Calculate CustomNonbondedForce energy
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(custom_system, integrator, platform)
        context.setPositions(new_positions)
        custom_energy = context.getState(getEnergy=True).getPotentialEnergy()
        custom_val = custom_energy.value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # Calculate relative difference
        abs_diff = abs(custom_val - ref_val)
        rel_diff = abs_diff / abs(ref_val) * 100 if abs(ref_val) > 1e-6 else abs_diff
        
        print(f"{dist:7.2f}  {ref_val:20.6f}  {custom_val:25.6f}  {rel_diff:12.6f}")
        
        # If difference is large, output detailed information
        if rel_diff > 0.05:  # Output detailed info when difference > 0.05%
            print(f"\n  Detailed information at {dist} nm:")
            print(f"    NonbondedForce energy:     {ref_val:.6f} kJ/mol")
            print(f"    CustomNonbondedForce energy: {custom_val:.6f} kJ/mol")
            print(f"    Absolute difference:        {abs_diff:.6f} kJ/mol")
            print(f"    Relative difference:        {rel_diff:.6f}%")
        
        # Verify results: use different tolerances based on distance
        if dist <= switching_distance:
            # Use stricter tolerance within switching distance
            assert rel_diff < 0.1, f"Energy difference too large at {dist} nm: {rel_diff:.6f}% > 0.1%"
        elif dist < cutoff:
            # Use looser tolerance in switching region
            assert rel_diff < 0.5, f"Energy difference too large in switching region at {dist} nm: {rel_diff:.6f}% > 0.5%"
        else:
            # Energy should be close to zero beyond cutoff
            assert abs_diff < 2e-1, f"Energy should be close to zero beyond cutoff at {dist} nm, but difference is {abs_diff:.6f} kJ/mol"
    
    print("\n=== Test Summary ===")
    print("Energy calculations at all distances are within acceptable error ranges")
    print(f"- Within switching distance (r ≤ {switching_distance} nm): relative error < 0.1%")
    print(f"- Switching region ({switching_distance} nm < r < {cutoff} nm): relative error < 0.5%")
    print(f"- Beyond cutoff (r ≥ {cutoff} nm): absolute error < 0.2 kJ/mol")
