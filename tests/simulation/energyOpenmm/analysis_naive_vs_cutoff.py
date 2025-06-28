# tests/simulation/energyOpenmm/analysis_naive_vs_cutoff.py

import pytest
from .analysis_helpers import *

def test_compare_naive_vs_cutoff_energy():
    """
    Compare energy differences between simple formula and cutoff formula.
    
    Tests two implementations:
    1. Simple formula (hard cutoff):
       E = step(cutoff - r) * (kC * q₁q₂/r + 4ε[(σ/r)¹² - (σ/r)⁶])
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    2. Cutoff formula:
       - Shifted Coulomb: E_coul = kC * q₁q₂ * (1/r - 1/r_cutoff)
       - LJ with switching: E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] * S(r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Test distances:
    - Inside switching region (0.9-1.0 nm)
    - Beyond cutoff (>1.0 nm)
    
    Verification:
    - Energy properly goes to zero at cutoff
    - Switching function properly applied
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define test distances
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.2]  # nm
    movement_atoms = set(range(6))  # Benzene carbons
    fixed_atoms = set(range(6, 9))  # Water atoms
    
    print("\nComparing simple formula and cutoff formula energies:")
    print("Distance(nm)  Simple(kJ/mol)  Cutoff(kJ/mol)  Difference(%)")
    print("-" * 60)
    
    # Add debug function
    def analyze_switching_function(r):
        """Analyze switching function value at given distance r"""
        cutoff = 1.0
        switch = 0.9
        if r <= switch:
            return 1.0
        elif r >= cutoff:
            return 0.0
        else:
            x = (cutoff - r)**2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)**3)
            return x
    
    def debug_energy_components(r, nb_force, movement_atoms, fixed_atoms, platform):
        """Analyze energy components at distance r"""
        # 1. Calculate Coulomb term only
        coulomb_expression = """
        step(cutoff - r) * kC * q1 * q2 * (1/r - 1/cutoff)
        """
        coulomb_force = CustomNonbondedForce(coulomb_expression)
        coulomb_force.addPerParticleParameter("q")
        coulomb_force.addGlobalParameter("kC", 138.935456)
        coulomb_force.addGlobalParameter("cutoff", 1.0)
        
        # 2. Calculate LJ term (without switching function)
        lj_expression = """
        step(cutoff - r) * 4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
        """
        lj_force = CustomNonbondedForce(lj_expression)
        lj_force.addPerParticleParameter("sigma")
        lj_force.addPerParticleParameter("eps")
        lj_force.addGlobalParameter("cutoff", 1.0)
        
        # 3. Calculate LJ term (with switching function)
        lj_switched_expression = """
        step(cutoff - r) * 4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
        """
        lj_switched_force = CustomNonbondedForce(lj_switched_expression)
        lj_switched_force.addPerParticleParameter("sigma")
        lj_switched_force.addPerParticleParameter("eps")
        lj_switched_force.addGlobalParameter("cutoff", 1.0)
        lj_switched_force.addGlobalParameter("switch", 0.9)
        
        # Add particle parameters to all forces
        for i in range(nb_force.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            coulomb_force.addParticle([charge])
            lj_force.addParticle([sigma, epsilon])
            lj_switched_force.addParticle([sigma, epsilon])
        
        # Set interaction groups and cutoff method
        for force in [coulomb_force, lj_force, lj_switched_force]:
            force.addInteractionGroup(movement_atoms, fixed_atoms)
            force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
            force.setCutoffDistance(1.0 * nanometers)
        
        # Create system and calculate energy
        def calc_energy(force):
            sys = System()
            for i in range(nb_force.getNumParticles()):
                sys.addParticle(system.getParticleMass(i))
            sys.addForce(force)
            integrator = VerletIntegrator(0.001 * picoseconds)
            context = Context(sys, integrator, platform)
            context.setPositions(new_positions)
            energy = context.getState(getEnergy=True).getPotentialEnergy()
            del context, integrator
            return energy.value_in_unit(kilojoules_per_mole)
        
        # Calculate component energies
        coulomb_energy = calc_energy(coulomb_force)
        lj_energy = calc_energy(lj_force)
        lj_switched_energy = calc_energy(lj_switched_force)
        
        # Calculate switching function value
        switch_value = analyze_switching_function(r)
        
        print(f"\n=== Energy Analysis (r = {r:.3f} nm) ===")
        print(f"Switching function value: {switch_value:.6f}")
        print(f"Coulomb energy: {coulomb_energy:.6f} kJ/mol")
        print(f"LJ energy (no switching): {lj_energy:.6f} kJ/mol")
        print(f"LJ energy (with switching): {lj_switched_energy:.6f} kJ/mol")
        print(f"LJ energy ratio (switched/unswitched): {lj_switched_energy/lj_energy if abs(lj_energy) > 1e-10 else 0:.6f}")
        
        return coulomb_energy, lj_energy, lj_switched_energy, switch_value
    
    for dist in distances:
        # Move water molecule to specified distance
        new_positions = []
        for i in range(len(positions)):
            if i >= 6 and i < 9:  # Water atoms
                pos = positions[i].value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # Get original NonbondedForce
        nb_force = None
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                nb_force = force
                break
        
        # Analyze energy components
        platform = Platform.getPlatformByName('Reference')
        coulomb_energy, lj_energy, lj_switched_energy, switch_value = debug_energy_components(
            dist, nb_force, movement_atoms, fixed_atoms, platform
        )
        
        # Calculate energy using simple formula
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
        
        # Add particle parameters
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            naive_force.addParticle([charge, sigma, epsilon])
        
        naive_force.addInteractionGroup(movement_atoms, fixed_atoms)
        naive_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        naive_force.setCutoffDistance(1.0 * nanometers)
        
        # Create system and calculate energy
        naive_system = System()
        for i in range(system.getNumParticles()):
            naive_system.addParticle(system.getParticleMass(i))
        naive_system.addForce(naive_force)
        
        # Calculate simple formula energy
        integrator_naive = VerletIntegrator(0.001 * picoseconds)
        context = Context(naive_system, integrator_naive, platform)
        context.setPositions(new_positions)
        naive_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator_naive
        
        # Calculate energy using cutoff formula
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
        
        cutoff_force = CustomNonbondedForce(cutoff_expression)
        cutoff_force.addPerParticleParameter("q")
        cutoff_force.addPerParticleParameter("sigma")
        cutoff_force.addPerParticleParameter("eps")
        cutoff_force.addGlobalParameter("kC", 138.935456)
        cutoff_force.addGlobalParameter("cutoff", 1.0)
        cutoff_force.addGlobalParameter("switch", 0.9)
        
        # Add particle parameters
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            cutoff_force.addParticle([charge, sigma, epsilon])
        
        cutoff_force.addInteractionGroup(movement_atoms, fixed_atoms)
        cutoff_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        cutoff_force.setCutoffDistance(1.0 * nanometers)
        
        # Create system and calculate energy
        cutoff_system = System()
        for i in range(system.getNumParticles()):
            cutoff_system.addParticle(system.getParticleMass(i))
        cutoff_system.addForce(cutoff_force)
        
        # Calculate energy
        integrator_cutoff = VerletIntegrator(0.001 * picoseconds)
        context = Context(cutoff_system, integrator_cutoff, platform)
        context.setPositions(new_positions)
        cutoff_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator_cutoff
        
        # Calculate differences
        naive_val = naive_energy.value_in_unit(kilojoules_per_mole)
        cutoff_val = cutoff_energy.value_in_unit(kilojoules_per_mole)
        
        # Calculate relative difference (use absolute difference if energy is close to zero)
        if abs(naive_val) < 1e-6:
            diff_percent = abs(cutoff_val - naive_val)
        else:
            diff_percent = abs(cutoff_val - naive_val) / abs(naive_val) * 100
        
        print(f"{dist:6.2f}  {naive_val:14.6f}  {cutoff_val:16.6f}  {diff_percent:8.2f}")
        
        # For distances beyond cutoff, cutoff formula should give 0 energy
        if dist > 1.0:  # cutoff distance
            assert abs(cutoff_val) < 1e-6, f"Energy should be zero beyond cutoff, got {cutoff_val}"
        
        # For distances close to cutoff, cutoff formula should give smaller energy
        if 0.9 < dist < 1.0:  # switching region
            # Use relative tolerance for comparison
            rel_tol = 1e-10  # relative tolerance: 1e-10
            abs_tol = 1e-10  # absolute tolerance: 1e-10 kJ/mol
            
            # If energy is small, use absolute tolerance; otherwise use relative tolerance
            if abs(naive_val) < 1e-6:
                assert abs(cutoff_val) <= abs_tol, \
                       f"Energy with switching should be near zero at {dist} nm, got {cutoff_val}"
            else:
                # Check if the energy with switching is smaller than or equal to (considering tolerance) the simple formula energy
                assert abs(cutoff_val) <= abs(naive_val) * (1 + rel_tol) + abs_tol, \
                       f"Energy with switching ({cutoff_val}) should be smaller than or equal to naive ({naive_val}) at {dist} nm"
                
                # Output detailed comparison information
                print(f"\nEnergy comparison details (r = {dist} nm):")
                print(f"Simple formula energy: {naive_val:.15f} kJ/mol")
                print(f"Cutoff formula energy: {cutoff_val:.15f} kJ/mol")
                print(f"Relative difference: {abs(cutoff_val - naive_val)/abs(naive_val)*100:.15f}%")
                print(f"Absolute difference: {abs(cutoff_val - naive_val):.15e} kJ/mol")
        
        del naive_system, cutoff_system

