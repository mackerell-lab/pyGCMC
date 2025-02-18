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

def test_compare_naive_vs_cutoff_energy():
    """比较简单公式和带截断公式的能量差别。
    
    简单公式（带硬截断）：
    E = step(cutoff - r) * (
        kC * q1 * q2 / r + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )
    )
        
    带截断公式（使用移位库伦势和LJ切换函数）：
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
    # 创建测试系统
    system, topology, positions = create_test_system()
    
    # 定义测试距离
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.2]  # nm
    movement_atoms = set(range(6))  # Benzene carbons
    fixed_atoms = set(range(6, 9))  # Water atoms
    
    print("\n比较简单公式和带截断公式的能量差别：")
    print("距离(nm)  简单公式(kJ/mol)  带截断公式(kJ/mol)  差异(%)")
    print("-" * 60)
    
    # 添加调试函数
    def analyze_switching_function(r):
        """分析在给定距离r处的切换函数值"""
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
        """分析在距离r处的能量组分"""
        # 1. 只计算库伦项
        coulomb_expression = """
        step(cutoff - r) * kC * q1 * q2 * (1/r - 1/cutoff)
        """
        coulomb_force = CustomNonbondedForce(coulomb_expression)
        coulomb_force.addPerParticleParameter("q")
        coulomb_force.addGlobalParameter("kC", 138.935456)
        coulomb_force.addGlobalParameter("cutoff", 1.0)
        
        # 2. 只计算LJ项（不带切换函数）
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
        
        # 3. 只计算LJ项（带切换函数）
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
        
        # 添加粒子参数到所有力场
        for i in range(nb_force.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            coulomb_force.addParticle([charge])
            lj_force.addParticle([sigma, epsilon])
            lj_switched_force.addParticle([sigma, epsilon])
        
        # 设置相互作用组和截断方法
        for force in [coulomb_force, lj_force, lj_switched_force]:
            force.addInteractionGroup(movement_atoms, fixed_atoms)
            force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
            force.setCutoffDistance(1.0 * nanometers)
        
        # 创建系统并计算能量
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
        
        # 计算各组分能量
        coulomb_energy = calc_energy(coulomb_force)
        lj_energy = calc_energy(lj_force)
        lj_switched_energy = calc_energy(lj_switched_force)
        
        # 计算切换函数值
        switch_value = analyze_switching_function(r)
        
        print(f"\n=== 能量分析 (r = {r:.3f} nm) ===")
        print(f"切换函数值: {switch_value:.6f}")
        print(f"库伦能量: {coulomb_energy:.6f} kJ/mol")
        print(f"LJ能量 (无切换): {lj_energy:.6f} kJ/mol")
        print(f"LJ能量 (带切换): {lj_switched_energy:.6f} kJ/mol")
        print(f"LJ能量比例 (带切换/无切换): {lj_switched_energy/lj_energy if abs(lj_energy) > 1e-10 else 0:.6f}")
        
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
    """比较三种方法（简单公式、自定义截断和标准OpenMM）在考虑自能补偿后的能量差异。
    
    所有方法都使用相同的相互作用组（只计算benzene和water分子之间的相互作用）。
    """
    # 创建测试系统
    system, topology, positions = create_test_system()
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 定义相互作用组
    movement_atoms = set(range(6))  # benzene
    fixed_atoms = set(range(6, 9))  # water
    
    # 计算自能补偿项
    correction = calculate_self_energy_correction(original_nb_force, 1.0)
    print(f"\nSelf-energy correction: {correction:.6f} kJ/mol")
    print("(Note: Standard OpenMM already handles self-energy correction internally)")
    
    # 定义测试距离
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]
    
    print("\nDetailed Energy Comparison:")
    print("Distance (nm) | Simple (kJ/mol) | Custom (kJ/mol) | Standard (kJ/mol) | Max Diff (kJ/mol)")
    print("-" * 100)
    
    platform = Platform.getPlatformByName('Reference')
    
    for dist in distances:
        # 生成新的位置
        new_positions = []
        for i, pos in enumerate(positions):
            if 6 <= i < 9:  # water分子
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)
        
        # 1. 简单公式（硬截断）
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
        
        # 2. 自定义截断公式（带移位和切换）
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
        # 不再应用自能补偿，因为标准OpenMM已经在内部处理了
        custom_energy_corrected = custom_energy
        del context, integrator
        
        # 3. 标准OpenMM NonbondedForce
        sys_standard = System()
        for i in range(system.getNumParticles()):
            sys_standard.addParticle(system.getParticleMass(i))
        
        nb_force = NonbondedForce()
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nb_force.addParticle(charge, sigma, epsilon)
        
        # 设置相互作用组：通过添加例外（exceptions）来实现
        # 将所有不需要计算的相互作用设置为0
        for i in range(original_nb_force.getNumParticles()):
            for j in range(i+1, original_nb_force.getNumParticles()):
                # 如果两个原子不是一个在movement_atoms一个在fixed_atoms，就设置为例外
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
        
        # 计算最大差异
        energies = [simple_energy, custom_energy_corrected, standard_energy]
        max_diff = max([abs(e1 - e2) for e1 in energies for e2 in energies])
        
        print(f"{dist:11.2f} | {simple_energy:13.6f} | {custom_energy_corrected:17.6f} | {standard_energy:15.6f} | {max_diff:16.6f}")
        
        # 如果差异太大，输出详细信息
        if max_diff > 1.0:  # 差异大于1 kJ/mol时输出详细信息
            print(f"  Detailed differences at {dist} nm:")
            print(f"  Custom-Simple: {abs(custom_energy_corrected - simple_energy):.6f} kJ/mol")
            print(f"  Standard-Simple: {abs(standard_energy - simple_energy):.6f} kJ/mol")
            print(f"  Standard-Custom: {abs(standard_energy - custom_energy_corrected):.6f} kJ/mol")

def test_analyze_openmm_energy_terms():
    """分析OpenMM的能量计算公式，分别计算库伦项和LJ项。"""
    # 创建测试系统
    system, topology, positions = create_test_system()
    
    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # 定义相互作用组
    movement_atoms = set(range(6))  # benzene
    fixed_atoms = set(range(6, 9))  # water
    
    # 定义测试距离
    distances = [0.35, 0.5, 0.7, 0.9, 0.95, 1.0, 1.1, 1.2]
    
    print("\nAnalyzing OpenMM energy terms:")
    print("Distance (nm) | Coulomb (kJ/mol) | LJ (kJ/mol) | Total (kJ/mol) | Standard OpenMM (kJ/mol)")
    print("-" * 100)
    
    platform = Platform.getPlatformByName('Reference')
    
    for dist in distances:
        # 生成新的位置
        new_positions = []
        for i, pos in enumerate(positions):
            if 6 <= i < 9:  # water分子
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)
        
        # 1. 计算移位库伦能
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
        
        # 2. 计算带切换函数的LJ能
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
        
        # 3. 计算标准OpenMM能量作为参考
        sys_standard = System()
        for i in range(system.getNumParticles()):
            sys_standard.addParticle(system.getParticleMass(i))
        
        nb_force = NonbondedForce()
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nb_force.addParticle(charge, sigma, epsilon)
        
        # 设置相互作用组
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
        
        # 计算总能量（库伦 + LJ）
        total_energy = coulomb_energy + lj_energy
        
        print(f"{dist:11.2f} | {coulomb_energy:15.6f} | {lj_energy:11.6f} | {total_energy:13.6f} | {standard_energy:21.6f}")
        
        # 如果与标准OpenMM结果差异较大，输出详细信息
        diff = abs(total_energy - standard_energy)
        if diff > 0.1:
            print(f"  Large difference at {dist} nm:")
            print(f"  Difference between sum and standard: {diff:.6f} kJ/mol")
            print(f"  Relative difference: {diff/abs(standard_energy)*100 if abs(standard_energy) > 1e-10 else 0:.6f}%")

def test_cutoff_periodic_comparison():
    """比较 CustomNonbondedForce 与 NonbondedForce 在 CutoffPeriodic 模式下的能量。
    
    本测试验证在周期性边界条件下，CustomNonbondedForce 使用反应场库仑势和 LJ 切换函数
    计算的能量是否与 OpenMM 标准 NonbondedForce 的结果一致。
    
    测试内容包括：
    1. 使用反应场库仑势（reaction-field electrostatics）
    2. 带切换函数的 LJ 势
    3. 周期性边界条件
    """
    # 创建测试系统
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    epsilon_rf = 78.5  # 水的相对介电常数

    # 获取原始 NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("系统中未找到 NonbondedForce")

    # 定义能量表达式
    # 注意：反应场势已经包含了长程效应的校正，不需要额外的自能校正
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

    # 测试不同距离
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.1]
    print("\n测试不同位置的能量:")
    print("\n距离(nm)  CustomNonbondedForce  NonbondedForce    差异(%)")
    print("-" * 60)

    # 使用 Reference 平台
    platform = Platform.getPlatformByName('Reference')

    for dist in distances:
        # 移动水分子到新位置
        new_positions = []
        for i, pos in enumerate(positions):
            if i >= 6 and i < 9:  # water分子
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)

        # 创建新的 CustomNonbondedForce
        custom_force = CustomNonbondedForce(custom_energy_expression)
        custom_force.addPerParticleParameter("q")
        custom_force.addPerParticleParameter("sigma")
        custom_force.addPerParticleParameter("epsilon")
        custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb 常数
        custom_force.addGlobalParameter("cutoff", cutoff_distance)
        custom_force.addGlobalParameter("switch", switch_distance)
        custom_force.addGlobalParameter("epsilon_rf", epsilon_rf)

        # 从 NonbondedForce 复制粒子参数
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            custom_force.addParticle([charge, sigma, epsilon])

        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        custom_force.setCutoffDistance(cutoff_distance * nanometers)

        # 创建 CustomNonbondedForce 系统
        custom_system = System()
        for i in range(system.getNumParticles()):
            custom_system.addParticle(system.getParticleMass(i))
        custom_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        custom_system.addForce(custom_force)

        # 计算 CustomNonbondedForce 能量
        custom_integrator = VerletIntegrator(0.001 * picoseconds)
        custom_context = Context(custom_system, custom_integrator, platform)
        custom_context.setPositions(new_positions)
        custom_energy = custom_context.getState(getEnergy=True).getPotentialEnergy()
        
        # 注意：不再添加移位库仑势的自能校正，因为反应场方法已经包含了长程效应的校正

        # 清理 CustomNonbondedForce 资源
        del custom_context, custom_integrator, custom_system

        # 创建新的 NonbondedForce
        ref_force = NonbondedForce()
        ref_force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        ref_force.setCutoffDistance(cutoff_distance * nanometers)
        ref_force.setUseSwitchingFunction(True)
        ref_force.setSwitchingDistance(switch_distance * nanometers)
        ref_force.setReactionFieldDielectric(epsilon_rf)

        # 复制粒子参数
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            ref_force.addParticle(charge, sigma, epsilon)

        # 创建参考系统
        ref_system = System()
        for i in range(system.getNumParticles()):
            ref_system.addParticle(system.getParticleMass(i))
        ref_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        ref_system.addForce(ref_force)

        # 计算参考能量
        ref_integrator = VerletIntegrator(0.001 * picoseconds)
        ref_context = Context(ref_system, ref_integrator, platform)
        ref_context.setPositions(new_positions)
        ref_energy = ref_context.getState(getEnergy=True).getPotentialEnergy()

        # 清理参考系统资源
        del ref_context, ref_integrator, ref_system

        # 计算差异
        custom_val = custom_energy.value_in_unit(kilojoules_per_mole)
        ref_val = ref_energy.value_in_unit(kilojoules_per_mole)
        abs_diff = abs(custom_val - ref_val)

        # 使用更稳健的相对差异计算方法
        if abs(ref_val) > 1e-6:
            rel_diff = abs_diff / abs(ref_val) * 100
        else:
            rel_diff = abs_diff * 100 if abs_diff > 1e-6 else 0.0

        print(f"{dist:6.2f}    {custom_val:16.6f}    {ref_val:12.6f}    {rel_diff:8.4f}")

        # 如果差异较大，输出详细信息
        if rel_diff > 0.05:  # 当差异超过0.05%时输出详细信息
            print(f"\n  距离 {dist} nm 处的详细信息:")
            print(f"    CustomNonbondedForce: {custom_val:.6f} kJ/mol")
            print(f"    NonbondedForce:       {ref_val:.6f} kJ/mol")
            print(f"    绝对差异:             {abs_diff:.6f} kJ/mol")
            print(f"    相对差异:             {rel_diff:.6f}%")

        # 验证结果：根据距离使用不同的容差
        if dist <= switch_distance:
            # 在切换距离内使用较严格的容差
            assert rel_diff < 0.1, f"在距离 {dist} nm 处能量差异过大: {rel_diff:.6f}% > 0.1%"
        elif dist < cutoff_distance:
            # 在切换区域使用较宽松的容差
            assert rel_diff < 0.5, f"在切换区域 {dist} nm 处能量差异过大: {rel_diff:.6f}% > 0.5%"
        else:
            # 在截断距离之外，能量应该非常接近零
            assert abs_diff < 2e-1, f"在截断距离外 {dist} nm 处能量应该接近零，但差异为 {abs_diff:.6f} kJ/mol"

def test_separate_lj_coulomb_periodic():
    """分别比较周期性边界条件下的LJ项和反应场库仑项。
    
    本测试将CustomNonbondedForce的能量分解为：
    1. LJ项：带切换函数的Lennard-Jones势
    2. Coulomb项：反应场静电势
    
    分别与标准OpenMM NonbondedForce的结果进行比较，以便更好地理解
    不同项的贡献和误差来源。
    """
    # 创建测试系统
    system, topology, positions = create_test_system()
    cutoff_distance = 1.0  # nm
    switch_distance = 0.9  # nm
    epsilon_rf = 78.5  # 水的相对介电常数

    # 获取原始NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    if original_nb_force is None:
        raise ValueError("系统中未找到NonbondedForce")

    # 定义LJ能量表达式（带切换函数）
    lj_expression = """
    4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * sw;
    epsilon = sqrt(epsilon1*epsilon2);
    sigma = 0.5*(sigma1+sigma2);
    sw = step(cutoff - r) * (step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3) + step(switch - r));
    """

    # 定义反应场库仑能量表达式
    coulomb_expression = """
    kC * q1 * q2 * (1/r + krf * r^2 - crf);
    krf = (epsilon_rf - 1) / (2*epsilon_rf + 1) / cutoff^3;
    crf = (3*epsilon_rf) / (2*epsilon_rf + 1) / cutoff;
    """

    # 测试不同距离
    distances = [0.5, 0.7, 0.9, 0.95, 1.0, 1.1]
    print("\n分别测试LJ项和Coulomb项的能量:")
    print("\n距离(nm)  |  Custom LJ  |  Custom Coulomb  |  Total Custom  |  OpenMM Total  |  差异(%)")
    print("-" * 85)

    platform = Platform.getPlatformByName('Reference')

    for dist in distances:
        # 移动water分子到新位置
        new_positions = []
        for i, pos in enumerate(positions):
            if i >= 6 and i < 9:  # water分子
                pos_val = pos.value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos_val[1], pos_val[2]) * nanometers)
            else:
                new_positions.append(pos)

        # 1. 计算LJ能量
        lj_force = CustomNonbondedForce(lj_expression)
        lj_force.addPerParticleParameter("sigma")
        lj_force.addPerParticleParameter("epsilon")
        lj_force.addGlobalParameter("cutoff", cutoff_distance)
        lj_force.addGlobalParameter("switch", switch_distance)

        # 添加粒子参数（只需要LJ参数）
        for i in range(original_nb_force.getNumParticles()):
            _, sigma, epsilon = original_nb_force.getParticleParameters(i)
            lj_force.addParticle([sigma, epsilon])

        lj_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        lj_force.setCutoffDistance(cutoff_distance * nanometers)

        # 创建LJ系统
        lj_system = System()
        for i in range(system.getNumParticles()):
            lj_system.addParticle(system.getParticleMass(i))
        lj_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        lj_system.addForce(lj_force)

        # 计算LJ能量
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(lj_system, integrator, platform)
        context.setPositions(new_positions)
        lj_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 2. 计算Coulomb能量
        coulomb_force = CustomNonbondedForce(coulomb_expression)
        coulomb_force.addPerParticleParameter("q")
        coulomb_force.addGlobalParameter("kC", 138.935456)
        coulomb_force.addGlobalParameter("cutoff", cutoff_distance)
        coulomb_force.addGlobalParameter("epsilon_rf", epsilon_rf)

        # 添加粒子参数（只需要电荷）
        for i in range(original_nb_force.getNumParticles()):
            charge, _, _ = original_nb_force.getParticleParameters(i)
            coulomb_force.addParticle([charge])

        coulomb_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        coulomb_force.setCutoffDistance(cutoff_distance * nanometers)

        # 创建Coulomb系统
        coulomb_system = System()
        for i in range(system.getNumParticles()):
            coulomb_system.addParticle(system.getParticleMass(i))
        coulomb_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        coulomb_system.addForce(coulomb_force)

        # 计算Coulomb能量
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(coulomb_system, integrator, platform)
        context.setPositions(new_positions)
        coulomb_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # 3. 计算标准OpenMM能量作为参考
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

        # 转换为kJ/mol并计算总能量
        lj_val = lj_energy.value_in_unit(kilojoules_per_mole)
        coulomb_val = coulomb_energy.value_in_unit(kilojoules_per_mole)
        custom_total = lj_val + coulomb_val
        ref_val = ref_energy.value_in_unit(kilojoules_per_mole)

        # 计算相对差异
        abs_diff = abs(custom_total - ref_val)
        rel_diff = abs_diff / abs(ref_val) * 100 if abs(ref_val) > 1e-6 else abs_diff

        print(f"{dist:8.2f} | {lj_val:10.4f} | {coulomb_val:14.4f} | {custom_total:12.4f} | {ref_val:13.4f} | {rel_diff:8.4f}")

        # 如果差异较大，输出详细信息
        if rel_diff > 0.05:  # 当差异超过0.05%时输出详细信息
            print(f"\n  距离 {dist} nm 处的详细信息:")
            print(f"    LJ能量:        {lj_val:.6f} kJ/mol")
            print(f"    Coulomb能量:   {coulomb_val:.6f} kJ/mol")
            print(f"    Custom总能量:  {custom_total:.6f} kJ/mol")
            print(f"    OpenMM能量:    {ref_val:.6f} kJ/mol")
            print(f"    绝对差异:      {abs_diff:.6f} kJ/mol")
            print(f"    相对差异:      {rel_diff:.6f}%")

        # 验证结果：根据距离使用不同的容差
        if dist <= switch_distance:
            # 在切换距离内使用较严格的容差
            assert rel_diff < 0.1, f"在距离 {dist} nm 处能量差异过大: {rel_diff:.6f}% > 0.1%"
        elif dist < cutoff_distance:
            # 在切换区域使用较宽松的容差
            assert rel_diff < 0.5, f"在切换区域 {dist} nm 处能量差异过大: {rel_diff:.6f}% > 0.5%"
        else:
            # 在截断距离之外，能量应该接近零
            assert abs_diff < 2e-1, f"在截断距离外 {dist} nm 处能量应该接近零，但差异为 {abs_diff:.6f} kJ/mol"
