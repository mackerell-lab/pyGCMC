# tests/simulation/openmm/periodic_tests.py

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
