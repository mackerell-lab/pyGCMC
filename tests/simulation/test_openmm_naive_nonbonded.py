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
    # Each carbon has charge 0 (for testing LJ only)
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
    # All charges set to 0 (for testing LJ only)
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
        nb_force.addParticle(-0.1, 0.34, 0.36)  # charge=-0.1e, sigma=0.34nm, epsilon=0.36kJ/mol
    
    # Add fixed molecule atoms (TIP3P-like parameters)
    nb_force.addParticle(-0.834, 0.3166, 0.650)  # O: charge=-0.834e
    nb_force.addParticle(0.417, 0.0, 0.0)       # H1: charge=0.417e
    nb_force.addParticle(0.417, 0.0, 0.0)       # H2: charge=0.417e
    
    # Set up periodic boundary conditions
    nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)  # 改为NoCutoff，与naive实现一致
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

def calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=False, custom_expression=None):
    """Calculate nonbonded energy between specified groups using direct Coulomb and LJ potentials.
    
    The energy expression includes both electrostatic and van der Waals terms:
    E = kC * q1 * q2 / r +  # direct Coulomb
        4 * sqrt(eps1*eps2) * ((sigma/r)^12 - (sigma/r)^6)  # LJ
    
    Args:
        system: OpenMM System object
        positions: List of Vec3 positions
        movement_atoms: Set of atom indices for movement group
        fixed_atoms: Set of atom indices for fixed group
        use_pbc: Whether to use periodic boundary conditions
        custom_expression: Optional custom energy expression to use instead of default
    """
    # 默认的nonbonded能量表达式，包括Coulomb和LJ
    if custom_expression is None:
        energy_expression = """
        kC * q1 * q2 / r + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )"""
    else:
        energy_expression = custom_expression
    
    custom_force = CustomNonbondedForce(energy_expression)
    
    # 添加每个粒子的参数
    custom_force.addPerParticleParameter("q")      # 电荷
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # 添加全局参数
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb常数 (kJ·nm/mol/e^2)
    
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
        custom_force.setCutoffDistance(1.0 * nanometers)  # 设置与naive实现相同的截断距离
    else:
        custom_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    
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
    state.forcefield.numTotalTypes = 3  # 三种类型：C、O和H
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
    
    # 注意：OpenMM的epsilon中已经包含了4的因子，所以不需要除以4
    c_c_eps = c_params[2].value_in_unit(kilojoules_per_mole)
    c_c_sigma = c_params[1].value_in_unit(nanometers)
    o_o_eps = o_params[2].value_in_unit(kilojoules_per_mole)
    o_o_sigma = o_params[1].value_in_unit(nanometers)
    h_h_eps = h_params[2].value_in_unit(kilojoules_per_mole)  # 应该是0
    h_h_sigma = h_params[1].value_in_unit(nanometers)         # 应该是0
    
    # 使用OpenMM的混合规则
    def mix_params(eps1, sigma1, eps2, sigma2):
        if eps1 == 0 or eps2 == 0 or sigma1 == 0 or sigma2 == 0:
            return 0.0, 0.0
        # OpenMM的混合规则：
        # - sigma: 算术平均 0.5*(sigma1 + sigma2)
        # - epsilon: 几何平均 sqrt(eps1 * eps2)
        mixed_sigma = 0.5 * (sigma1 + sigma2)  # 修改这里，使用0.5*(sigma1 + sigma2)
        mixed_eps = math.sqrt(eps1 * eps2)
        return mixed_eps, mixed_sigma
    
    # 计算混合参数
    c_o_eps, c_o_sigma = mix_params(c_c_eps, c_c_sigma, o_o_eps, o_o_sigma)
    c_h_eps, c_h_sigma = mix_params(c_c_eps, c_c_sigma, h_h_eps, h_h_sigma)
    o_h_eps, o_h_sigma = mix_params(o_o_eps, o_o_sigma, h_h_eps, h_h_sigma)
    
    # 设置力场参数矩阵 (numTotalTypes * numTotalTypes = 3 * 3)
    # 完整的交互矩阵:
    # [C-C, C-O, C-H]
    # [O-C, O-O, O-H]
    # [H-C, H-O, H-H]
    state.forcefield.ljEps = [
        c_c_eps, c_o_eps, c_h_eps,    # C与(C,O,H)的相互作用
        c_o_eps, o_o_eps, o_h_eps,    # O与(C,O,H)的相互作用
        c_h_eps, o_h_eps, h_h_eps     # H与(C,O,H)的相互作用
    ]
    state.forcefield.ljSigma = [
        c_c_sigma, c_o_sigma, c_h_sigma,    # C与(C,O,H)的相互作用
        c_o_sigma, o_o_sigma, o_h_sigma,    # O与(C,O,H)的相互作用
        c_h_sigma, o_h_sigma, h_h_sigma     # H与(C,O,H)的相互作用
    ]
    
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
            atom.type = 2  # Hydrogen type (separate from oxygen)
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

def test_openmm_energy_components():
    """Test OpenMM energy components separately."""
    state, system, positions = convert_openmm_state_to_mcstate()
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    
    # 计算总能量
    total_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # 只计算静电能
    elec_expression = """
    kC * q1 * q2 / r;
    """
    elec_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, 
                                           custom_expression=elec_expression)
    
    # 只计算LJ能
    lj_expression = """
    4 * sqrt(eps1*eps2) * (
        (0.5*(sigma1+sigma2)/r)^12 - 
        (0.5*(sigma1+sigma2)/r)^6
    );
    """
    lj_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms,
                                         custom_expression=lj_expression)
    
    print(f"\nOpenMM energy components:")
    print(f"Total energy: {total_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Electrostatic: {elec_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"LJ: {lj_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # 验证总能量约等于分量之和
    total_val = total_energy.value_in_unit(kilojoules_per_mole)
    components_sum = (elec_energy.value_in_unit(kilojoules_per_mole) + 
                     lj_energy.value_in_unit(kilojoules_per_mole))
    assert abs(total_val - components_sum) < 1e-6, \
           f"Energy components don't sum to total: {total_val} != {components_sum}"

def test_naive_energy_components():
    """Test naive implementation energy components."""
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 计算能量
    pygcmc.computeMovementResiduesEnergy(state)
    
    print(f"\nNaive implementation energy components:")
    print(f"VDW energy: {state.residues[0].energy_vdw:.6f} kJ/mol")
    print(f"Elec energy: {state.residues[0].energy_elec:.6f} kJ/mol")
    print(f"Total energy: {(state.residues[0].energy_vdw + state.residues[0].energy_elec):.6f} kJ/mol")

def print_force_field_params():
    """Print force field parameters for both implementations."""
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 打印OpenMM参数
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
            
    print("\nOpenMM parameters:")
    for i in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        print(f"Atom {i}: q={charge.value_in_unit(elementary_charge):.3f}e, "
              f"sigma={sigma.value_in_unit(nanometers):.3f}nm, "
              f"epsilon={epsilon.value_in_unit(kilojoules_per_mole):.3f}kJ/mol")
        
    print("\nNaive implementation parameters:")
    print("LJ Epsilon matrix [kJ/mol]:")
    n = int(math.sqrt(len(state.forcefield.ljEps)))
    for i in range(n):
        row = state.forcefield.ljEps[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")
    
    print("\nLJ Sigma matrix [nm]:")
    for i in range(n):
        row = state.forcefield.ljSigma[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")

def test_compare_openmm_naive_nonbonded():
    """Compare nonbonded energy calculations between OpenMM and naive implementation."""
    # 首先打印所有调试信息
    print("\n=== Force Field Parameters ===")
    print_force_field_params()
    
    print("\n=== OpenMM Energy Components ===")
    test_openmm_energy_components()
    
    print("\n=== Naive Implementation Energy Components ===")
    test_naive_energy_components()
    
    # 原始的比较测试代码
    state, system, positions = convert_openmm_state_to_mcstate()
    
    movement_atoms = set(range(6))  # Benzene carbons
    fixed_atoms = set(range(6, 9))  # Water atoms
    openmm_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
    
    pygcmc.computeMovementResiduesEnergy(state)
    naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    print(f"\n=== Final Energy Comparison ===")
    print(f"OpenMM energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive energy: {naive_energy:.6f} kJ/mol")
    print(f"Absolute difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - naive_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    rel_tol = 0.001  # 0.%的相对误差容忍度，因为两种实现的细节可能有所不同
    assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
           f"Energy mismatch: OpenMM={openmm_energy_val}, Naive={naive_energy}"

def test_compare_pbc_energies():
    """Compare PBC energy calculations between OpenMM and naive implementation."""
    # 获取系统
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # 设置一致的截断距离
    state.info.cutoff = 1.0  # 与OpenMM保持一致：1.0 nm
    
    # 移动水分子到盒子边缘以测试PBC
    box_size = state.info.box[0]  # 3.0 nm
    new_positions = []
    for i in range(len(positions)):
        pos = positions[i].value_in_unit(nanometers)
        if i >= 6:  # Water molecule atoms
            # 移动到2.5 nm处，这样与原点的距离为2.5 nm，通过PBC的距离为0.5 nm
            new_positions.append(Vec3(2.5, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    
    # 验证PBC距离计算
    # 以第一个碳原子（约-0.15 nm）和水的氧原子（2.5 nm）为例
    dx = 2.5 - (-0.15)  # 原始距离 = 2.65 nm
    dx_pbc = dx - box_size * round(dx/box_size)  # 应该约为 -0.35 nm
    print(f"\nPBC distance validation:")
    print(f"Original distance: {dx:.6f} nm")
    print(f"After PBC: {dx_pbc:.6f} nm")
    print(f"Should be in range [-{box_size/2:.1f}, {box_size/2:.1f}] nm")
    assert abs(dx_pbc) <= box_size/2, "PBC distance calculation error"
    
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
    
    # 计算naive PBC能量 - 只计算movement residue的能量
    pygcmc.computeFullSystemCutoffPBCEnergy(state)
    movement_fixed_energy = state.residues[1].energy_vdw + state.residues[1].energy_elec
    
    # 打印详细信息
    print(f"\nDetailed PBC energy comparison:")
    print(f"Box size: {box_size} nm")
    print(f"Cutoff distance: {state.info.cutoff} nm")
    print(f"Water molecule position: {new_positions[6].value_in_unit(nanometers)} nm")
    print(f"OpenMM PBC energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive movement-fixed energy: {movement_fixed_energy:.6f} kJ/mol")
    print(f"Difference: {abs(openmm_energy_val - movement_fixed_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - movement_fixed_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    # 验证结果
    rel_tol = 0.3  # 30%的相对误差容忍度，因为两种实现的细节可能有所不同
    assert abs(openmm_energy_val - movement_fixed_energy) / abs(openmm_energy_val) < rel_tol, \
           f"PBC energy mismatch: OpenMM={openmm_energy_val}, Naive={movement_fixed_energy}"

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
            rel_tol = 0.3  # 30%的相对误差容忍度，因为两种实现的细节可能有所不同
            assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
                   f"Energy mismatch at {dist} nm: OpenMM={openmm_energy_val}, Naive={naive_energy}"

