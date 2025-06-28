# tests/simulation/energyOpenmm/analysis_detailed_comparison.py

import pytest
from .analysis_helpers import *

def test_detailed_energy_comparison_simple_vs_openmm_cutoff():
    """
    Detailed comparison of energy between simple formula and OpenMM cutoff formula.
    
    Tests energy calculations at multiple distances with:
    1. Simple cutoff formula:
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    2. OpenMM cutoff formula with:
       - Reaction field electrostatics
       - LJ switching function
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Output table includes:
    - simple_energy: Energy from simple cutoff formula
    - cutoff_energy: Energy from OpenMM cutoff formula
    - abs_diff: Absolute difference
    - rel_diff: Relative difference
    
    System configuration:
    - Uses test system from create_test_system()
    - Water molecule position varied along x-axis
    """
    # Get initial system, topology and positions
    system, topology, positions = create_test_system()
    # Save original parameters from NonbondedForce (for extracting particle parameters and calculating self-energy correction)
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("No NonbondedForce found in system")
    
    # Define interacting atom groups: benzene molecule atoms 0-5 and water molecule atoms 6-8
    movement_atoms = set(range(6))  # benzene
    fixed_atoms = set(range(6, 9))    # water

    # Define energy expressions
    # (1) Simple formula: hard cutoff, no shift or switching
    naive_expression = """
    step(cutoff - r) * (
        kC * q1 * q2 / r +
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 -
            (0.5*(sigma1+sigma2)/r)^6
        )
    )"""
    # (2) OpenMM cutoff formula: shifted Coulomb (1/r - 1/cutoff) and LJ with switching function
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
    # Helper function for creating new CustomNonbondedForce instances
    def create_naive_force():
        force = CustomNonbondedForce(naive_expression)
        force.addPerParticleParameter("q")
        force.addPerParticleParameter("sigma")
        force.addPerParticleParameter("eps")
        force.addGlobalParameter("kC", 138.935456)
        force.addGlobalParameter("cutoff", 1.0)
        # Extract all particle parameters from original_nb_force
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            force.addParticle([charge, sigma, epsilon])
        force.addInteractionGroup(movement_atoms, fixed_atoms)
        force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        force.setCutoffDistance(1.0 * nanometers)
        return force

    def create_cutoff_force():
        force = CustomNonbondedForce(cutoff_expression)
        force.addPerParticleParameter("q")
        force.addPerParticleParameter("sigma")
        force.addPerParticleParameter("eps")
        force.addGlobalParameter("kC", 138.935456)
        force.addGlobalParameter("cutoff", 1.0)
        force.addGlobalParameter("switch", 0.9)
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            force.addParticle([charge, sigma, epsilon])
        force.addInteractionGroup(movement_atoms, fixed_atoms)
        force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        force.setCutoffDistance(1.0 * nanometers)
        return force

    # Define test distances (nm); includes equilibrium region, switching region, and beyond cutoff
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]

    print("\nDetailed Energy Comparison: Simple vs OpenMM Cutoff")
    print("Distance (nm) | Simple Energy (kJ/mol) | OpenMM Cutoff Energy (kJ/mol) | Abs Diff (kJ/mol) | Rel Diff (%)")
    print("-" * 100)

    platform = Platform.getPlatformByName('Reference')
    # Calculate self-energy correction (Note: using cutoff = 1.0 nm)
    correction = calculate_self_energy_correction(original_nb_force, 1.0)
    print(f"\nSelf-energy correction: {correction:.6f} kJ/mol")
    print("(Note: Standard OpenMM already handles self-energy correction internally)")
    
    for dist in distances:
        # Generate new positions: set water molecule (atoms 6-8) x-coordinate to dist, keep others unchanged
        new_positions = []
        for i, pos in enumerate(positions):
            if 6 <= i < 9:  # water molecule
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)
        
        # --- Calculate simple formula energy ---
        sys_naive = System()
        for i in range(system.getNumParticles()):
            sys_naive.addParticle(system.getParticleMass(i))
        naive_force = create_naive_force()
        sys_naive.addForce(naive_force)
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_naive, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        simple_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # --- Calculate OpenMM cutoff formula energy ---
        sys_cutoff = System()
        for i in range(system.getNumParticles()):
            sys_cutoff.addParticle(system.getParticleMass(i))
        cutoff_force = create_cutoff_force()
        sys_cutoff.addForce(cutoff_force)
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(sys_cutoff, integrator, platform)
        context.setPositions(new_positions)
        state = context.getState(getEnergy=True)
        cutoff_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        del context, integrator
        
        # For OpenMM cutoff formula, add self-energy correction (Note: original test used subtraction)
        cutoff_energy_corrected = cutoff_energy - correction
        
        # Calculate absolute and relative differences (when simple_energy is close to 0, only compare absolute difference)
        abs_diff = abs(cutoff_energy_corrected - simple_energy)
        if abs(simple_energy) > 1e-6:
            rel_diff = abs_diff / abs(simple_energy) * 100
        else:
            rel_diff = 0.0
        
        
        print(f"{dist:13.2f} | {simple_energy:22.6f} | {cutoff_energy_corrected:29.6f} | {abs_diff:16.6f} | {rel_diff:11.6f}")

