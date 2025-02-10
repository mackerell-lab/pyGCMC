# /tests/simulation/test_openmm_naive_nonbonded.py

import pytest
import numpy as np
import math
import pygcmc
import os
import warnings

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
    
    # 将system中的约束(constraints)复制到energy_system中
    for i in range(system.getNumConstraints()):
        energy_system.addConstraint(*system.getConstraintParameters(i))
        
    energy_system.addForce(custom_force)
    
    # 使用Reference平台计算能量
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(energy_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    
    # 打印详细的能量信息用于调试
    print(f"\nDetailed energy calculation:")
    print(f"Number of movement atoms: {len(movement_atoms)}")
    print(f"Number of fixed atoms: {len(fixed_atoms)}")
    print(f"PBC: {use_pbc}")
    print(f"Energy: {energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    del context, integrator, energy_system
    return energy

def convert_openmm_state_to_mcstate():
    """Convert OpenMM test system to MCState for naive implementation."""
    # 创建OpenMM测试系统
    system, topology, positions = create_test_system()
    
    # 创建MCState
    state = pygcmc.MCState()
    
    # 1. 设置力场
    state.forcefield.numTotalTypes = 2  # 两种类型：C和O (H的LJ参数为0)
    state.forcefield.numMovementTypes = 1  # C是movement类型
    
    # 从OpenMM系统中获取力场参数
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    
    # 获取参数
    c_params = nb_force.getParticleParameters(0)  # Carbon parameters
    o_params = nb_force.getParticleParameters(6)  # Oxygen parameters
    h_params = nb_force.getParticleParameters(7)  # Hydrogen parameters
    
    # 注意：在naive实现中，epsilon已经包含了factor 4，所以不需要除以4
    # 因为在naive实现的C++代码中已经包含了这个factor
    c_c_eps = c_params[2].value_in_unit(kilojoules_per_mole)
    c_c_sigma = c_params[1].value_in_unit(nanometers)
    o_o_eps = o_params[2].value_in_unit(kilojoules_per_mole)
    o_o_sigma = o_params[1].value_in_unit(nanometers)
    
    # C-O参数使用Lorentz-Berthelot混合规则
    c_o_sigma = (c_c_sigma + o_o_sigma) / 2
    c_o_eps = math.sqrt(c_c_eps * o_o_eps)
    
    # 设置力场参数
    # 注意：在naive实现中，参数矩阵是按照movement type和total type组织的
    # 对于每个movement type，需要它与所有total type的相互作用参数
    # 这里只有一个movement type (C)，它需要与两个total types (C和O)的相互作用参数
    state.forcefield.ljEps = [c_c_eps, c_o_eps]  # [C-C, C-O]
    state.forcefield.ljSigma = [c_c_sigma, c_o_sigma]  # [C-C, C-O]
    
    # 设置movement类型
    state.movementAtomTypes = [0]  # Type 0 (C) is movement type
    state.numMovementAtomTypes = 1
    
    # 2. 设置原子
    atoms = []
    # 添加苯环原子
    for i in range(6):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        atom.charge = c_params[0].value_in_unit(elementary_charge)
        atom.type = 0  # Carbon type
        atoms.append(atom)
    
    # 添加水分子原子
    for i in range(6, 9):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        if i == 6:  # Oxygen
            atom.charge = o_params[0].value_in_unit(elementary_charge)
            atom.type = 1  # Oxygen type
        else:  # Hydrogen
            atom.charge = h_params[0].value_in_unit(elementary_charge)
            atom.type = 1  # Same type as oxygen for simplicity
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # 3. 设置residues
    # Movement residue (benzene)
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    movement_res.atomStart = 0
    movement_res.atomCount = 6
    
    # Fixed residue (water)
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    fixed_res.atomStart = 6
    fixed_res.atomCount = 3
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. 设置movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # 5. 设置盒子和截断
    box_vectors = system.getDefaultPeriodicBoxVectors()
    state.info.box = [
        box_vectors[0][0].value_in_unit(nanometers),
        box_vectors[1][1].value_in_unit(nanometers),
        box_vectors[2][2].value_in_unit(nanometers)
    ]
    state.info.cutoff = nb_force.getCutoffDistance().value_in_unit(nanometers)
    
    return state, system, positions

def test_compare_openmm_naive_nonbonded():
    """Compare nonbonded energy calculations between OpenMM and naive implementation."""
    # 获取两种实现的系统
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 计算OpenMM能量
    movement_atoms = set(range(6))  # Benzene carbons
    fixed_atoms = set(range(6, 9))  # Water atoms
    openmm_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
    
    # 计算naive实现的能量
    pygcmc.computeMovementResiduesEnergy(state)
    naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # 打印详细信息用于调试
    print(f"\nDetailed energy comparison:")
    print(f"OpenMM energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive energy: {naive_energy:.6f} kJ/mol")
    print(f"Difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - naive_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    # 验证结果
    # 由于数值精度和实现细节的差异，我们使用相对宽松的容差
    rel_tol = 0.05  # 5%的相对误差容忍度
    assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
           f"Energy mismatch: OpenMM={openmm_energy_val}, Naive={naive_energy}"

def test_compare_pbc_energies():
    """Compare PBC energy calculations between OpenMM and naive implementation."""
    # 获取系统
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 移动水分子到盒子边缘以测试PBC
    box_size = state.info.box[0]
    new_positions = []
    for i in range(len(positions)):
        pos = positions[i].value_in_unit(nanometers)
        if i >= 6:  # Water molecule atoms
            new_positions.append(Vec3(box_size - 0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    
    # 计算OpenMM PBC能量
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    openmm_energy = calculate_nonbonded_energy(system, new_positions, movement_atoms, fixed_atoms, use_pbc=True)
    openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
    
    # 更新naive实现中的原子位置
    for i in range(6, 9):
        pos = new_positions[i].value_in_unit(nanometers)
        state.atoms[i].x = pos[0]
        state.atoms[i].y = pos[1]
        state.atoms[i].z = pos[2]
    
    # 计算naive PBC能量
    pygcmc.computeMovementResiduesEnergy(state)
    naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # 打印详细信息
    print(f"\nDetailed PBC energy comparison:")
    print(f"OpenMM PBC energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive PBC energy: {naive_energy:.6f} kJ/mol")
    print(f"Difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - naive_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    # 验证结果
    rel_tol = 0.05  # 5%的相对误差容忍度
    assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
           f"PBC energy mismatch: OpenMM={openmm_energy_val}, Naive={naive_energy}"

def test_compare_cutoff_effects():
    """Compare cutoff effects between OpenMM and naive implementation."""
    # 获取系统
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 测试不同距离
    distances = [0.5, 0.7, 0.9, 1.2]  # nm
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    
    for dist in distances:
        # 移动水分子
        new_positions = []
        for i in range(len(positions)):
            pos = positions[i].value_in_unit(nanometers)
            if i >= 6:  # Water molecule atoms
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # 计算OpenMM能量
        openmm_energy = calculate_nonbonded_energy(
            system, new_positions, movement_atoms, fixed_atoms, use_pbc=False
        )
        openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
        
        # 更新naive实现中的原子位置
        for i in range(6, 9):
            pos = new_positions[i].value_in_unit(nanometers)
            state.atoms[i].x = pos[0]
            state.atoms[i].y = pos[1]
            state.atoms[i].z = pos[2]
        
        # 计算naive能量
        pygcmc.computeMovementResiduesEnergy(state)
        naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
        
        print(f"\nEnergy comparison at distance {dist} nm:")
        print(f"OpenMM energy: {openmm_energy_val:.6f} kJ/mol")
        print(f"Naive energy: {naive_energy:.6f} kJ/mol")
        print(f"Difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
        
        # 对于超出截断的距离，两者都应该给出接近0的能量
        if dist > state.info.cutoff:
            assert abs(openmm_energy_val) < 1e-6, f"OpenMM energy not zero beyond cutoff: {openmm_energy_val}"
            assert abs(naive_energy) < 1e-6, f"Naive energy not zero beyond cutoff: {naive_energy}"
        else:
            # 对于截断内的距离，能量应该接近
            rel_tol = 0.05  # 5%的相对误差容忍度
            assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
                   f"Energy mismatch at {dist} nm: OpenMM={openmm_energy_val}, Naive={naive_energy}"

