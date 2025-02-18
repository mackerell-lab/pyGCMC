# /tests/simulation/test_openmm_nonbonded.py

import pytest
import os
import numpy as np
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
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def create_test_system():
    """Create a simple test system with a few molecules."""
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
    """Calculate nonbonded energy between specified groups using shifted Coulomb and LJ potentials.
    
    The energy expression includes both electrostatic and van der Waals terms:
    E = step(cutoff - r) * [
        kC * q1 * q2 * (1/r - 1/cutoff) +  # shifted Coulomb
        4 * sqrt(eps1*eps2) * ((sigma/r)^12 - (sigma/r)^6)  # LJ
    ]
    where:
    - step(x) is the Heaviside step function (0 for x < 0, 1 for x >= 0)
    - sigma = (sigma1 + sigma2)/2  # Lorentz-Berthelot mixing rule for sigma
    - eps = sqrt(eps1*eps2)        # Lorentz-Berthelot mixing rule for epsilon
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

    The standard NonbondedForce subtracts a self-energy term when using shifted Coulomb:
      U_self = - (kC/(2*r_cut)) * sum_i q_i^2.
    To make the custom pairwise energy sum match the standard implementation,
    we need to add this correction term (note: correction is negative).
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

def test_attractive_interaction():
    """Test nonbonded energy calculation for an attractive interaction."""
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be negative (attractive) due to opposite charges
    assert energy.value_in_unit(kilojoules_per_mole) < 0, \
           f"Expected attractive interaction, got {energy}"
    print(f"Attractive interaction energy: {energy}")

def test_repulsive_interaction():
    """Test nonbonded energy calculation for a repulsive interaction."""
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule very close to benzene to create repulsion
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be positive (repulsive) due to close distance
    assert energy.value_in_unit(kilojoules_per_mole) > 0, \
           f"Expected repulsive interaction, got {energy}"
    print(f"Repulsive interaction energy: {energy}")

def test_pbc_interaction():
    """Test nonbonded energy calculation with periodic boundary conditions."""
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule to box edge
    box_size = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(box_size - 0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy with and without PBC
    energy_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=True)
    energy_no_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=False)
    
    # Energy with PBC should be different from energy without PBC
    assert abs(energy_pbc.value_in_unit(kilojoules_per_mole) - 
              energy_no_pbc.value_in_unit(kilojoules_per_mole)) > 1e-3, \
           "PBC should affect the interaction energy"
    print(f"PBC interaction energy: {energy_pbc}")
    print(f"Non-PBC interaction energy: {energy_no_pbc}")

def test_cutoff_effect():
    """Test the effect of cutoff distance on nonbonded energy calculation.
    
    This test verifies that the nonbonded energy decreases with distance and
    becomes zero beyond the cutoff distance. We disable PBC to ensure we're
    measuring true distance-dependent effects.
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Place water molecule at different distances
    # Note: benzene carbons are at radius 0.15 nm, so we need x >= 1.15 nm
    # to ensure all atom pairs are beyond the 1.0 nm cutoff
    distances = [0.5, 0.7, 0.9, 1.2]  # nm (last distance ensures all pairs > cutoff)
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    energies = []
    for dist in distances:
        # Move water molecule
        new_positions = []
        for i in range(len(positions)):
            if i >= 6 and i < 9:  # Water atoms
                pos = positions[i].value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # Calculate energy without PBC to observe true distance effects
        energy = calculate_nonbonded_energy(system, new_positions, movement_atoms, fixed_atoms, use_pbc=False)
        energy_val = energy.value_in_unit(kilojoules_per_mole)
        energies.append(energy_val)
        print(f"Energy at distance {dist} nm: {energy_val:.4f} kJ/mol")
    
    # Verify energy decreases with distance
    for i in range(len(distances)-1):
        assert abs(energies[i]) > abs(energies[i+1]), \
               f"Energy should decrease with distance. Energies: {energies}"
    
    # Verify energy is zero beyond cutoff (1.0 nm)
    assert abs(energies[-1]) < 1e-6, \
           f"Energy should be zero beyond cutoff (1.0 nm), but got {energies[-1]} kJ/mol"

def test_energy_symmetry():
    """Test that energy calculation is symmetric (A->B equals B->A)."""
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define groups
    group1 = set(range(6))  # Benzene
    group2 = set(range(6, 9))  # Water
    
    # Calculate energy both ways
    energy1 = calculate_nonbonded_energy(system, positions, group1, group2)
    energy2 = calculate_nonbonded_energy(system, positions, group2, group1)
    
    # Energies should be equal
    assert abs(energy1.value_in_unit(kilojoules_per_mole) - 
              energy2.value_in_unit(kilojoules_per_mole)) < 1e-6, \
           f"Energy calculation should be symmetric: {energy1} != {energy2}"
    print(f"Forward energy: {energy1}")
    print(f"Reverse energy: {energy2}")

def test_compare_custom_vs_standard_nonbonded():
    """Compare energy calculations between CustomNonbondedForce and NonbondedForce."""
    # Create test system
    system, topology, positions = create_test_system()
    
    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # Create a new system containing only NonbondedForce
    system_standard = System()
    for i in range(system.getNumParticles()):
        system_standard.addParticle(system.getParticleMass(i))
    system_standard.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # Copy NonbondedForce to the new system
    nb_force = NonbondedForce()
    for i in range(original_nb_force.getNumParticles()):
        params = original_nb_force.getParticleParameters(i)
        nb_force.addParticle(*params)
    nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)
    system_standard.addForce(nb_force)
    
    # Create a system using CustomNonbondedForce
    system_custom = System()
    for i in range(system.getNumParticles()):
        system_custom.addParticle(system.getParticleMass(i))
    system_custom.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # Create CustomNonbondedForce using the same energy expression as NonbondedForce
    energy_expression = """
    kC * q1 * q2 / r + 
    4 * sqrt(eps1*eps2) * (
        (0.5*(sigma1+sigma2)/r)^12 - 
        (0.5*(sigma1+sigma2)/r)^6
    )"""
    
    custom_force = CustomNonbondedForce(energy_expression)
    custom_force.addPerParticleParameter("q")
    custom_force.addPerParticleParameter("sigma")
    custom_force.addPerParticleParameter("eps")
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb constant (kJ·nm/mol/e^2)
    
    # Add particle parameters
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        custom_force.addParticle([charge, sigma, epsilon])
    
    custom_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    system_custom.addForce(custom_force)
    
    platform = Platform.getPlatformByName('Reference')
    
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
    
    # Print results
    print(f"\nComparing NonbondedForce vs CustomNonbondedForce:")
    print(f"NonbondedForce energy: {energy_standard.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"CustomNonbondedForce energy: {energy_custom.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole)):.6f} kJ/mol")
    rel_diff = abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole))/abs(energy_standard.value_in_unit(kilojoules_per_mole))*100
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # Verify results
    assert rel_diff < 1e-6, "Energy mismatch between NonbondedForce and CustomNonbondedForce"
    
    del context_standard, context_custom

def test_compare_nonbonded_methods():
    """Compare different nonbonded methods between CustomNonbondedForce and NonbondedForce."""
    # 创建测试系统
    system, topology, positions = create_test_system()
    
    # 测试不同的非键方法
    methods = [
        (NonbondedForce.NoCutoff, CustomNonbondedForce.NoCutoff, "NoCutoff"),
        (NonbondedForce.CutoffNonPeriodic, CustomNonbondedForce.CutoffNonPeriodic, "CutoffNonPeriodic")
    ]
    
    # 为不同方法定义不同的容差
    tolerances = {
        "NoCutoff": 1e-6,           # 无截断时要求更高的精度
        "CutoffNonPeriodic": 5e-4   # 使用截断时允许 0.05% 的误差
    }
    
    platform = Platform.getPlatformByName('Reference')
    
    for nb_method, custom_method, method_name in methods:
        print(f"\nTesting {method_name}:")
        
        # 创建标准NonbondedForce系统
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
        switch_distance = 0.9  # nm, 切换函数开始的距离
        if nb_method != NonbondedForce.NoCutoff:
            nb_force.setCutoffDistance(cutoff_distance * nanometers)
            nb_force.setUseSwitchingFunction(True)
            nb_force.setSwitchingDistance(switch_distance * nanometers)
        system_standard.addForce(nb_force)
        
        # 创建CustomNonbondedForce系统
        system_custom = System()
        for i in range(system.getNumParticles()):
            system_custom.addParticle(system.getParticleMass(i))
        system_custom.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        
        # 修改能量表达式
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
        
        # 计算标准NonbondedForce的能量
        integrator_standard = VerletIntegrator(0.001 * picoseconds)
        context_standard = Context(system_standard, integrator_standard, platform)
        context_standard.setPositions(positions)
        state_standard = context_standard.getState(getEnergy=True)
        energy_standard = state_standard.getPotentialEnergy()
        
        # 计算CustomNonbondedForce的能量
        integrator_custom = VerletIntegrator(0.001 * picoseconds)
        context_custom = Context(system_custom, integrator_custom, platform)
        context_custom.setPositions(positions)
        state_custom = context_custom.getState(getEnergy=True)
        energy_custom = state_custom.getPotentialEnergy()
        
        # 如果使用截断，减去自能补偿
        if custom_method != CustomNonbondedForce.NoCutoff:
            correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
            energy_custom = energy_custom - correction * kilojoules_per_mole
        
        # 计算与参考能量的差异
        energy_diff = abs(energy_standard.value_in_unit(kilojoules_per_mole) - 
                         energy_custom.value_in_unit(kilojoules_per_mole))
        rel_diff = energy_diff / abs(energy_standard.value_in_unit(kilojoules_per_mole)) * 100
        
        print(f"NonbondedForce energy: {energy_standard.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
        print(f"CustomNonbondedForce energy: {energy_custom.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
        print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
        print(f"Relative difference: {rel_diff:.6f}%")
        
        # 使用对应方法的容差进行验证
        rel_tol = tolerances[method_name]
        print(f"Using relative tolerance: {rel_tol:.6e}")
        
        assert rel_diff/100 < rel_tol, \
               f"Energy mismatch for {method_name}"
        
        del context_standard, context_custom

def test_compare_switching_functions():
    """Test nonbonded interactions with switching function.
    
    This test compares the energy calculations between:
    1. With switching function
    2. Without switching function
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

def test_compare_separate_terms():
    """Compare Coulomb and LJ terms separately.
    
    This test verifies that the energy calculated by CustomNonbondedForce matches
    the results from standard OpenMM NonbondedForce by comparing:
    1. Coulomb term
    2. LJ term
    3. Total energy
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
    """Compare energy differences between simple formula and cutoff formula.
    
    Simple formula (with hard cutoff):
    E = step(cutoff - r) * (
        kC * q1 * q2 / r + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )
        
    Cutoff formula (using shifted Coulomb potential and LJ switching function):
    E = step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
    )
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
        # 移动水分子到指定距离
        new_positions = []
        for i in range(len(positions)):
            if i >= 6 and i < 9:  # Water atoms
                pos = positions[i].value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # 获取原始NonbondedForce
        nb_force = None
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                nb_force = force
                break
        
        # 分析能量组分
        platform = Platform.getPlatformByName('Reference')
        coulomb_energy, lj_energy, lj_switched_energy, switch_value = debug_energy_components(
            dist, nb_force, movement_atoms, fixed_atoms, platform
        )
        
        # 计算简单公式的能量
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
        
        # 添加粒子参数
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            naive_force.addParticle([charge, sigma, epsilon])
        
        naive_force.addInteractionGroup(movement_atoms, fixed_atoms)
        naive_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        naive_force.setCutoffDistance(1.0 * nanometers)
        
        # 创建系统并计算能量
        naive_system = System()
        for i in range(system.getNumParticles()):
            naive_system.addParticle(system.getParticleMass(i))
        naive_system.addForce(naive_force)
        
        # 计算简单公式能量
        integrator_naive = VerletIntegrator(0.001 * picoseconds)
        context = Context(naive_system, integrator_naive, platform)
        context.setPositions(new_positions)
        naive_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator_naive
        
        # 计算带截断公式能量
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
        
        # 添加粒子参数
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            cutoff_force.addParticle([charge, sigma, epsilon])
        
        cutoff_force.addInteractionGroup(movement_atoms, fixed_atoms)
        cutoff_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        cutoff_force.setCutoffDistance(1.0 * nanometers)
        
        # 创建系统并计算能量
        cutoff_system = System()
        for i in range(system.getNumParticles()):
            cutoff_system.addParticle(system.getParticleMass(i))
        cutoff_system.addForce(cutoff_force)
        
        # 计算能量
        integrator_cutoff = VerletIntegrator(0.001 * picoseconds)
        context = Context(cutoff_system, integrator_cutoff, platform)
        context.setPositions(new_positions)
        cutoff_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator_cutoff
        
        # 计算差异
        naive_val = naive_energy.value_in_unit(kilojoules_per_mole)
        cutoff_val = cutoff_energy.value_in_unit(kilojoules_per_mole)
        
        # 计算相对差异（如果能量接近0，使用绝对差异）
        if abs(naive_val) < 1e-6:
            diff_percent = abs(cutoff_val - naive_val)
        else:
            diff_percent = abs(cutoff_val - naive_val) / abs(naive_val) * 100
        
        print(f"{dist:6.2f}  {naive_val:14.6f}  {cutoff_val:16.6f}  {diff_percent:8.2f}")
        
        # 对于超出截断距离的情况，带截断公式应该给出0能量
        if dist > 1.0:  # cutoff distance
            assert abs(cutoff_val) < 1e-6, f"Energy should be zero beyond cutoff, got {cutoff_val}"
        
        # 对于接近截断距离的情况，带截断公式应该给出较小的能量
        if 0.9 < dist < 1.0:  # switching region
            # 使用相对容差进行比较
            rel_tol = 1e-10  # 相对容差：1e-10
            abs_tol = 1e-10  # 绝对容差：1e-10 kJ/mol
            
            # 如果能量很小，使用绝对容差；否则使用相对容差
            if abs(naive_val) < 1e-6:
                assert abs(cutoff_val) <= abs_tol, \
                       f"Energy with switching should be near zero at {dist} nm, got {cutoff_val}"
            else:
                # 检查带切换的能量是否小于或等于（考虑容差）简单公式的能量
                assert abs(cutoff_val) <= abs(naive_val) * (1 + rel_tol) + abs_tol, \
                       f"Energy with switching ({cutoff_val}) should be smaller than or equal to naive ({naive_val}) at {dist} nm"
                
                # 输出详细的比较信息
                print(f"\n能量比较详情 (r = {dist} nm):")
                print(f"简单公式能量: {naive_val:.15f} kJ/mol")
                print(f"带切换能量: {cutoff_val:.15f} kJ/mol")
                print(f"相对差异: {abs(cutoff_val - naive_val)/abs(naive_val)*100:.15f}%")
                print(f"绝对差异: {abs(cutoff_val - naive_val):.15e} kJ/mol")
        
        del naive_system, cutoff_system

def test_detailed_energy_comparison_simple_vs_openmm_cutoff():
    """
    在多个距离下详细比较简单公式（硬截断）和 OpenMM 截断公式（移位库伦和 LJ 切换函数）
    计算的能量差异。输出的表格中包括每个距离下：
      - simple_energy：使用简单截断公式计算的能量（kJ/mol）
      - cutoff_energy：使用 OpenMM 截断公式（移位/切换函数）计算的能量，经自能补偿后的值（kJ/mol）
      - abs_diff：两者的绝对差异（kJ/mol）
      - rel_diff：两者的相对差异（百分比）
    注意：本测试中使用的系统与粒子参数与 create_test_system() 中定义的保持一致，
    而 water 分子（原系统中最后 3 个原子）的 x 坐标会被设置为指定的测试距离（nm）。
    """
    # 获取初始系统、拓扑结构和初始位置
    system, topology, positions = create_test_system()
    # 保存 NonbondedForce 中的原始参数（用于提取粒子参数和计算自能补偿）
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("系统中未找到 NonbondedForce")
    
    # 定义参与相互作用的原子组：benzene 分子的原子编号 0-5 与 water 分子的编号 6-8
    movement_atoms = set(range(6))  # benzene
    fixed_atoms = set(range(6, 9))    # water

    # 定义能量表达式
    # （1）简单公式：硬截断，不含移位或切换
    naive_expression = """
    step(cutoff - r) * (
        kC * q1 * q2 / r +
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 -
            (0.5*(sigma1+sigma2)/r)^6
        )
    )"""
    # （2）OpenMM 截断公式：移位库伦（1/r - 1/cutoff）和 LJ 带切换函数
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
    # 用于每次新建 CustomNonbondedForce 的辅助函数
    def create_naive_force():
        force = CustomNonbondedForce(naive_expression)
        force.addPerParticleParameter("q")
        force.addPerParticleParameter("sigma")
        force.addPerParticleParameter("eps")
        force.addGlobalParameter("kC", 138.935456)
        force.addGlobalParameter("cutoff", 1.0)
        # 将所有粒子的参数从 original_nb_force 中提取出来
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

    # 定义测试距离（nm）；其中包括平衡区、切换区和超过截断距离的情况
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]

    print("\nDetailed Energy Comparison: Simple vs OpenMM Cutoff")
    print("Distance (nm) | Simple Energy (kJ/mol) | OpenMM Cutoff Energy (kJ/mol) | Abs Diff (kJ/mol) | Rel Diff (%)")
    print("-" * 100)

    platform = Platform.getPlatformByName('Reference')
    # 计算自能补偿项（注意：计算时 cutoff 使用 1.0 nm）
    correction = calculate_self_energy_correction(original_nb_force, 1.0)
    print(f"\nSelf-energy correction: {correction:.6f} kJ/mol")
    print("(Note: Standard OpenMM already handles self-energy correction internally)")
    
    for dist in distances:
        # 生成新的位置：将 water 分子（原子编号 6-8）的 x 坐标设置为 dist，其余保持不变
        new_positions = []
        for i, pos in enumerate(positions):
            if 6 <= i < 9:
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)
        
        # --- 计算简单公式能量 ---
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
        
        # --- 计算 OpenMM 截断公式能量 ---
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
        
        # 对于 OpenMM 截断公式，需要加上自能补偿（注意：原测试中使用减去 correction）
        cutoff_energy_corrected = cutoff_energy - correction
        
        # 计算绝对和相对差异（当 simple_energy 接近 0 时，仅比较绝对差异）
        abs_diff = abs(cutoff_energy_corrected - simple_energy)
        if abs(simple_energy) > 1e-6:
            rel_diff = abs_diff / abs(simple_energy) * 100
        else:
            rel_diff = 0.0
        
        print(f"{dist:13.2f} | {simple_energy:22.6f} | {cutoff_energy_corrected:29.6f} | {abs_diff:16.6f} | {rel_diff:11.6f}")

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
    """Analyze OpenMM energy terms by calculating Coulomb and LJ terms separately."""
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
    """Compare energy between CustomNonbondedForce and NonbondedForce in CutoffPeriodic mode.
    
    This test verifies that the energy calculated by CustomNonbondedForce using reaction-field
    Coulomb potential and LJ switching function matches the results from standard OpenMM
    NonbondedForce in periodic boundary conditions.
    
    Test includes:
    1. Reaction-field electrostatics
    2. LJ potential with switching function
    3. Periodic boundary conditions
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
    """Compare LJ and reaction field Coulomb terms separately in periodic boundary conditions.
    
    This test decomposes the CustomNonbondedForce energy into:
    1. LJ term: Lennard-Jones potential with switching function
    2. Coulomb term: Reaction field electrostatic potential
    
    Each term is compared with standard OpenMM NonbondedForce results to better understand
    the contributions and sources of error from different terms.
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
    """Compare parameter settings and energy calculations between NonbondedForce and CustomNonbondedForce.
    
    This test will:
    1. Get parameter settings from standard NonbondedForce
    2. Apply these parameters to CustomNonbondedForce
    3. Compare energy results from both methods
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
