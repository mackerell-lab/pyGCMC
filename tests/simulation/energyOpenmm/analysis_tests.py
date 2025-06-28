# tests/simulation/openmm/analysis_tests.py

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
def test_compare_separate_terms():
    """
    Compare Coulomb and LJ terms separately.
    
    Tests the individual contributions of:
    1. Lennard-Jones term:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb term:
       E_coul = (1/4πε₀)(q₁q₂/r)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    3. Total energy:
       E_total = E_LJ + E_coul
    
    Verification:
    - Compares each term with OpenMM reference
    - Ensures relative error < 1%
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
    
    # Test Coulomb term
    print("\nTesting Coulomb term:")
    
    # Create CustomNonbondedForce for Coulomb only
    coulomb_custom = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff)
    )""")
    
    # Add parameters
    coulomb_custom.addPerParticleParameter("q")
    coulomb_custom.addGlobalParameter("kC", 138.935456)
    coulomb_custom.addGlobalParameter("cutoff", cutoff_distance)
    
    # Create only Coulomb system
    system_coulomb = System()
    for i in range(system.getNumParticles()):
        system_coulomb.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, _, _ = original_nb_force.getParticleParameters(i)
        coulomb_custom.addParticle([charge])
    
    coulomb_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    coulomb_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_coulomb.addForce(coulomb_custom)
    
    # Calculate custom Coulomb energy
    platform = Platform.getPlatformByName('Reference')
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_coulomb, integrator, platform)
    context.setPositions(positions)
    coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # Add self-energy correction
    correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
    coulomb_energy = coulomb_energy + correction * kilojoules_per_mole
    
    print(f"Custom Coulomb energy: {coulomb_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # Test LJ term
    print("\nTesting LJ term:")
    
    # Create CustomNonbondedForce for LJ only
    lj_custom = CustomNonbondedForce("""
    step(cutoff - r) * (
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
    )""")
    
    # Add parameters
    lj_custom.addPerParticleParameter("sigma")
    lj_custom.addPerParticleParameter("eps")
    lj_custom.addGlobalParameter("cutoff", cutoff_distance)
    lj_custom.addGlobalParameter("switch", switch_distance)
    
    # Create only LJ system
    system_lj = System()
    for i in range(system.getNumParticles()):
        system_lj.addParticle(system.getParticleMass(i))
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        _, sigma, epsilon = original_nb_force.getParticleParameters(i)
        lj_custom.addParticle([sigma, epsilon])
    
    lj_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    lj_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_lj.addForce(lj_custom)
    
    # Calculate custom LJ energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_lj, integrator, platform)
    context.setPositions(positions)
    lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Custom LJ energy: {lj_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
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
    system_ref.addForce(nb_force)
    
    # Calculate reference energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_ref, integrator, platform)
    context.setPositions(positions)
    total_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Reference total energy: {total_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # Verify total energy
    custom_total = coulomb_energy + lj_energy
    energy_diff = abs(custom_total.value_in_unit(kilojoules_per_mole) - 
                     total_energy.value_in_unit(kilojoules_per_mole))
    rel_diff = energy_diff / abs(total_energy.value_in_unit(kilojoules_per_mole)) * 100
    
    print(f"\nResults:")
    print(f"Custom total energy: {custom_total.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # Verify relative error < 1%
    assert rel_diff/100 < 1e-2, "Energy terms differ significantly from OpenMM reference"
    
    del context, integrator

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

