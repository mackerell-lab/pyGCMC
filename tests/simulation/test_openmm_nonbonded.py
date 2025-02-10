# /tests/simulation/test_openmm_nonbonded.py

import pytest
import os
import numpy as np
import warnings

import math

# 提前屏蔽 SWIG 相关的 DeprecationWarning
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
    
    # Set up periodic boundary conditions
    nb_force.setNonbondedMethod(NonbondedForce.PME)
    nb_force.setCutoffDistance(1.0 * nanometers)
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
    # 完整的nonbonded能量表达式，包括Coulomb和LJ
    energy_expression = """
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    );
    """
    custom_force = CustomNonbondedForce(energy_expression)
    
    # 添加每个粒子的参数
    custom_force.addPerParticleParameter("q")      # 电荷
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # 添加全局参数
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb常数 (kJ·nm/mol/e^2)
    custom_force.addGlobalParameter("cutoff", 1.0)     # 截断距离 (nm)
    
    # 从原始NonbondedForce中获取参数
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    if nb_force is None:
        raise ValueError("No NonbondedForce found in system")
    
    # 添加粒子参数
    num_particles = system.getNumParticles()
    for i in range(num_particles):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        # OpenMM的sigma单位是nm，epsilon单位是kJ/mol
        custom_force.addParticle([charge, sigma, epsilon])
    
    # 只计算指定组之间的相互作用
    custom_force.addInteractionGroup(movement_atoms, fixed_atoms)
    
    # 设置非键方法和截断距离
    if use_pbc:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    else:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    custom_force.setCutoffDistance(1.0 * nanometers)
    
    # 创建只包含custom force的能量系统
    energy_system = System()
    for i in range(num_particles):
        energy_system.addParticle(system.getParticleMass(i))
    energy_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    energy_system.addForce(custom_force)
    
    # 使用Reference平台计算能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(energy_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    del context, integrator, energy_system
    return energy

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

