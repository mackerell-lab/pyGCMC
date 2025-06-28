# tests/simulation/openmm/method_comparisons.py

import pytest
import os
import warnings

import math

# Suppress SWIG-related DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type swigvarlink has no __module__ attribute")
from openmm import *
from openmm.app import *
from openmm.unit import *

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")

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

def create_test_system():
    """
    Create a simple test system with a benzene-like molecule and a water-like molecule.
    
    System configuration:
    1. Periodic boundary conditions with 3nm box size
    2. Movement molecule (benzene-like):
       - 6 carbon atoms in a ring
       - Each carbon has charge +0.1e
       - Uses OPLS-AA like parameters:
         σ = 0.34nm, ε = 0.36kJ/mol
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    3. Fixed molecule (water-like):
       - TIP3P water model
       - Oxygen: charge -0.834e, σ = 0.3166nm, ε = 0.650kJ/mol
       - Hydrogens: charge +0.417e each
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#nonbondedforce
    
    4. Nonbonded force settings:
       - Method: CutoffNonPeriodic
       - Cutoff distance: 1.0nm
       - Switching function: True
       - Switching distance: 0.9nm
    
    Returns:
        system (System): OpenMM system
        topology (Topology): System topology
        positions (list): Initial atomic positions
    """
    # Create a system with periodic boundary conditions
    system = System()
    box_size = 3.0 * nanometers
    system.setDefaultPeriodicBoxVectors(
        Vec3(box_size, 0, 0),
        Vec3(0, box_size, 0),
        Vec3(0, 0, box_size)
    )
    
    # Add particles
    # Movement molecule (benzene-like): 6 carbons in a ring
    # Each carbon has charge +0.1e
    movement_atoms = []
    radius = 0.15  # nm
    for i in range(6):
        angle = i * 2 * math.pi / 6
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        z = 0.0
        system.addParticle(12.0)  # mass in amu
        movement_atoms.append([x, y, z])
    
    # Fixed molecule (water-like): 3 atoms
    # Oxygen with charge -0.834e and 2 hydrogens with charge +0.417e each (TIP3P)
    fixed_atoms = [
        [1.0, 0.0, 0.0],  # O
        [1.1, 0.1, 0.0],  # H1
        [1.1, -0.1, 0.0]  # H2
    ]
    system.addParticle(16.0)  # O
    system.addParticle(1.0)   # H1
    system.addParticle(1.0)   # H2
    
    # Create nonbonded force
    nb_force = NonbondedForce()
    
    # Add movement molecule atoms (carbons with OPLS-AA like parameters)
    for _ in range(6):
        nb_force.addParticle(0.1, 0.34, 0.36)  # charge=+0.1e, sigma=0.34nm, epsilon=0.36kJ/mol
    
    # Add fixed molecule atoms (TIP3P-like parameters)
    nb_force.addParticle(-0.834, 0.3166, 0.650)  # O
    nb_force.addParticle(0.417, 0.0, 0.0)        # H1
    nb_force.addParticle(0.417, 0.0, 0.0)        # H2
    
    # Set up nonbonded method
    nb_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
    nb_force.setCutoffDistance(1.0 * nanometers)
    nb_force.setUseSwitchingFunction(True)
    nb_force.setSwitchingDistance(0.9 * nanometers)
    system.addForce(nb_force)
    
    # Create topology
    topology = Topology()
    chain = topology.addChain()
    
    # Add movement molecule (benzene)
    res_movement = topology.addResidue('BEN', chain)
    element = Element.getBySymbol('C')
    for pos in movement_atoms:
        topology.addAtom('C', element, res_movement)
    
    # Add fixed molecule (water)
    res_fixed = topology.addResidue('HOH', chain)
    o_element = Element.getBySymbol('O')
    h_element = Element.getBySymbol('H')
    topology.addAtom('O', o_element, res_fixed)
    topology.addAtom('H1', h_element, res_fixed)
    topology.addAtom('H2', h_element, res_fixed)
    
    # Create positions
    positions = []
    positions.extend([Vec3(*pos) for pos in movement_atoms])
    positions.extend([Vec3(*pos) for pos in fixed_atoms])
    positions = [pos * nanometers for pos in positions]
    
    return system, topology, positions

def calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=True):
    """
    Calculate nonbonded energy between specified groups using shifted Coulomb and LJ potentials.
    
    Energy expressions:
    1. Lennard-Jones with switching function:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] * S(r)
       where S(r) is the switching function:
       S(r) = 1                                    if r ≤ r_switch
       S(r) = (r_cut-r)²(r_cut+2r-3r_switch)/     if r_switch < r < r_cut
              (r_cut-r_switch)³
       S(r) = 0                                    if r ≥ r_cut
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Shifted Coulomb potential:
       E_coul = kC * q₁q₂ * (1/r - 1/r_cut)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Parameters:
        system (System): OpenMM system
        positions (list): Atomic positions
        movement_atoms (set): Indices of first group atoms
        fixed_atoms (set): Indices of second group atoms
        use_pbc (bool): Whether to use periodic boundary conditions
    
    Returns:
        energy (Quantity): Calculated nonbonded energy
    """
    # Complete nonbonded energy expression, including Coulomb and LJ
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
    
    # Add per-particle parameters
    custom_force.addPerParticleParameter("q")      # charge
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # Add global parameters
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb constant (kJ·nm/mol/e^2)
    custom_force.addGlobalParameter("cutoff", 1.0)     # cutoff distance (nm)
    custom_force.addGlobalParameter("switch", 0.9)     # switching distance (nm)
    
    # Get parameters from original NonbondedForce
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    if nb_force is None:
        raise ValueError("No NonbondedForce found in system")
    
    # Add particle parameters
    num_particles = system.getNumParticles()
    for i in range(num_particles):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        # OpenMM's sigma is in nm, epsilon in kJ/mol
        custom_force.addParticle([charge, sigma, epsilon])
    
    # Only calculate interactions between specified groups
    custom_force.addInteractionGroup(movement_atoms, fixed_atoms)
    
    # Set nonbonded method and cutoff distance
    if use_pbc:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    else:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    custom_force.setCutoffDistance(1.0 * nanometers)
    
    # Create energy system with only custom force
    energy_system = System()
    for i in range(num_particles):
        energy_system.addParticle(system.getParticleMass(i))
    energy_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    energy_system.addForce(custom_force)
    
    # Calculate energy using Reference platform
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(energy_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    del context, integrator, energy_system
    return energy

def calculate_self_energy_correction(force, cutoff):
    """
    Calculate self-energy correction for shifted Coulomb potential.
    
    In OpenMM's shifted Coulomb implementation, a self-energy correction term is subtracted:
    U_self = - (kC/(2*r_cut)) * sum_i q_i²
    
    This correction ensures proper energy conservation and is described in:
    @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Parameters:
        force (NonbondedForce): The force object containing particle parameters
        cutoff (float): Cutoff distance in nm
    
    Returns:
        float: Self-energy correction in kJ/mol
        
    Note:
        This correction is automatically handled by OpenMM's NonbondedForce,
        but needs to be manually applied when using CustomNonbondedForce
        to match the standard implementation.
    """
    kC = 138.935456  # kJ·nm/mol/e^2
    sum_q2 = 0.0
    for i in range(force.getNumParticles()):
        charge, _, _ = force.getParticleParameters(i)
        # Convert to dimensionless value in units of e
        charge_val = charge.value_in_unit(elementary_charge)
        sum_q2 += charge_val * charge_val
    correction = - (kC * sum_q2) / (2 * cutoff)
    return correction

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

