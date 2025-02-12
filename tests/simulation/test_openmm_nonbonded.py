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
    # 完整的nonbonded能量表达式，包括Coulomb和LJ
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
    
    # 添加每个粒子的参数
    custom_force.addPerParticleParameter("q")      # 电荷
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # 添加全局参数
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb常数 (kJ·nm/mol/e^2)
    custom_force.addGlobalParameter("cutoff", 1.0)     # 截断距离 (nm)
    custom_force.addGlobalParameter("switch", 0.9)     # 切换距离 (nm)
    
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

def calculate_self_energy_correction(force, cutoff):
    """
    计算移位库伦势下的自能补偿项。

    标准 NonbondedForce 在使用移位库伦势时扣除的自能为
      U_self = - (kC/(2*r_cut)) * sum_i q_i^2.
    为了使对成对相互作用求和的自定义能量与标准实现一致，
    需要加上该补偿项（注意：补偿项为负）。
    """
    kC = 138.935456  # kJ·nm/mol/e^2
    sum_q2 = 0.0
    for i in range(force.getNumParticles()):
        charge, _, _ = force.getParticleParameters(i)
        # 转换为以e为单位的无量纲数值
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
    # 创建测试系统
    system, topology, positions = create_test_system()
    
    # 获取原始的NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 创建一个新的系统，只包含NonbondedForce
    system_standard = System()
    for i in range(system.getNumParticles()):
        system_standard.addParticle(system.getParticleMass(i))
    system_standard.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # 复制NonbondedForce到新系统
    nb_force = NonbondedForce()
    for i in range(original_nb_force.getNumParticles()):
        params = original_nb_force.getParticleParameters(i)
        nb_force.addParticle(*params)
    nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)
    system_standard.addForce(nb_force)
    
    # 创建一个使用CustomNonbondedForce的系统
    system_custom = System()
    for i in range(system.getNumParticles()):
        system_custom.addParticle(system.getParticleMass(i))
    system_custom.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    
    # 创建CustomNonbondedForce，使用与NonbondedForce相同的能量表达式
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
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb常数 (kJ·nm/mol/e^2)
    
    # 添加粒子参数
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        custom_force.addParticle([charge, sigma, epsilon])
    
    custom_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    system_custom.addForce(custom_force)
    
    platform = Platform.getPlatformByName('Reference')
    
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
    
    # 打印结果
    print(f"\nComparing NonbondedForce vs CustomNonbondedForce:")
    print(f"NonbondedForce energy: {energy_standard.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"CustomNonbondedForce energy: {energy_custom.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole)):.6f} kJ/mol")
    rel_diff = abs(energy_standard.value_in_unit(kilojoules_per_mole) - energy_custom.value_in_unit(kilojoules_per_mole))/abs(energy_standard.value_in_unit(kilojoules_per_mole))*100
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # 验证结果
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
    """测试带切换函数的非键相互作用。"""
    # 创建测试系统
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 测试带切换函数的情况
    print("\nTesting with switching function:")
    
    # 创建带切换函数的自定义力
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
    
    # 创建不带切换函数的自定义力
    unswitched_force = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )""")
    
    # 添加参数
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
    
    # 创建系统
    system_switched = System()
    for i in range(system.getNumParticles()):
        system_switched.addParticle(system.getParticleMass(i))
    
    # 添加粒子参数
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        switched_force.addParticle([charge, sigma, epsilon])
        unswitched_force.addParticle([charge, sigma, epsilon])
    
    switched_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    switched_force.setCutoffDistance(cutoff_distance * nanometers)
    system_switched.addForce(switched_force)
    
    # 计算带切换函数的能量
    platform = Platform.getPlatformByName('Reference')
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_switched, integrator, platform)
    context.setPositions(positions)
    switched_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # 加上自能补偿项
    correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
    switched_energy = switched_energy + correction * kilojoules_per_mole
    
    print(f"Energy with switching: {switched_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # 测试不带切换函数的情况
    print("\nTesting without switching:")
    
    # 创建不带切换函数的自定义力
    unswitched_force = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )""")
    
    # 添加参数
    unswitched_force.addPerParticleParameter("q")
    unswitched_force.addPerParticleParameter("sigma")
    unswitched_force.addPerParticleParameter("eps")
    unswitched_force.addGlobalParameter("kC", 138.935456)
    unswitched_force.addGlobalParameter("cutoff", cutoff_distance)
    
    # 创建系统
    system_unswitched = System()
    for i in range(system.getNumParticles()):
        system_unswitched.addParticle(system.getParticleMass(i))
    
    # 添加粒子参数
    for i in range(original_nb_force.getNumParticles()):
        charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
        unswitched_force.addParticle([charge, sigma, epsilon])
    
    unswitched_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    unswitched_force.setCutoffDistance(cutoff_distance * nanometers)
    system_unswitched.addForce(unswitched_force)
    
    # 计算不带切换函数的能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_unswitched, integrator, platform)
    context.setPositions(positions)
    unswitched_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # 加上自能补偿项
    unswitched_energy = unswitched_energy + correction * kilojoules_per_mole
    
    print(f"Energy without switching: {unswitched_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # 创建标准参考系统
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
    
    # 计算参考能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_ref, integrator, platform)
    context.setPositions(positions)
    ref_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Reference energy: {ref_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # 验证结果
    print("\nResults:")
    
    # 验证带切换函数的能量与参考值的差异
    energy_diff = abs(switched_energy.value_in_unit(kilojoules_per_mole) - 
                     ref_energy.value_in_unit(kilojoules_per_mole))
    rel_diff = energy_diff / abs(ref_energy.value_in_unit(kilojoules_per_mole)) * 100
    
    print(f"Difference from reference (with switching):")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # 验证相对误差小于0.1%
    assert rel_diff/100 < 1e-3, "Energy with switching differs significantly from OpenMM reference"
    
    # 验证切换函数的效果
    energy_diff = abs(switched_energy.value_in_unit(kilojoules_per_mole) - 
                     unswitched_energy.value_in_unit(kilojoules_per_mole))
    print(f"\nDifference between switched and unswitched:")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    
    # 验证切换函数确实影响了能量
    assert energy_diff > 0, "Switching function should affect the energy"
    
    del context, integrator

def test_compare_separate_terms():
    """分别比较库伦项和LJ项的能量计算。"""
    # 创建测试系统
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 测试库伦项
    print("\nTesting Coulomb term:")
    
    # 创建只有库伦项的CustomNonbondedForce
    coulomb_custom = CustomNonbondedForce("""
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff)
    )""")
    
    # 添加参数
    coulomb_custom.addPerParticleParameter("q")
    coulomb_custom.addGlobalParameter("kC", 138.935456)
    coulomb_custom.addGlobalParameter("cutoff", cutoff_distance)
    
    # 创建只有库伦项的系统
    system_coulomb = System()
    for i in range(system.getNumParticles()):
        system_coulomb.addParticle(system.getParticleMass(i))
    
    # 添加粒子参数
    for i in range(original_nb_force.getNumParticles()):
        charge, _, _ = original_nb_force.getParticleParameters(i)
        coulomb_custom.addParticle([charge])
    
    coulomb_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    coulomb_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_coulomb.addForce(coulomb_custom)
    
    # 计算自定义库伦能量
    platform = Platform.getPlatformByName('Reference')
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_coulomb, integrator, platform)
    context.setPositions(positions)
    coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
    
    # 加上自能补偿项
    correction = calculate_self_energy_correction(original_nb_force, cutoff_distance)
    coulomb_energy = coulomb_energy + correction * kilojoules_per_mole
    
    print(f"Custom Coulomb energy: {coulomb_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # 测试LJ项
    print("\nTesting LJ term:")
    
    # 创建只有LJ项的CustomNonbondedForce
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
    
    # 添加参数
    lj_custom.addPerParticleParameter("sigma")
    lj_custom.addPerParticleParameter("eps")
    lj_custom.addGlobalParameter("cutoff", cutoff_distance)
    lj_custom.addGlobalParameter("switch", switch_distance)
    
    # 创建只有LJ项的系统
    system_lj = System()
    for i in range(system.getNumParticles()):
        system_lj.addParticle(system.getParticleMass(i))
    
    # 添加粒子参数
    for i in range(original_nb_force.getNumParticles()):
        _, sigma, epsilon = original_nb_force.getParticleParameters(i)
        lj_custom.addParticle([sigma, epsilon])
    
    lj_custom.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    lj_custom.setCutoffDistance(cutoff_distance * nanometers)
    system_lj.addForce(lj_custom)
    
    # 计算自定义LJ能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_lj, integrator, platform)
    context.setPositions(positions)
    lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Custom LJ energy: {lj_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    del context, integrator
    
    # 创建标准参考系统
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
    
    # 计算参考能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = Context(system_ref, integrator, platform)
    context.setPositions(positions)
    total_energy = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Reference total energy: {total_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # 验证总能量
    custom_total = coulomb_energy + lj_energy
    energy_diff = abs(custom_total.value_in_unit(kilojoules_per_mole) - 
                     total_energy.value_in_unit(kilojoules_per_mole))
    rel_diff = energy_diff / abs(total_energy.value_in_unit(kilojoules_per_mole)) * 100
    
    print(f"\nResults:")
    print(f"Custom total energy: {custom_total.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Absolute difference: {energy_diff:.6f} kJ/mol")
    print(f"Relative difference: {rel_diff:.6f}%")
    
    # 验证相对误差小于1%
    assert rel_diff/100 < 1e-2, "Energy terms differ significantly from OpenMM reference"
    
    del context, integrator

