#!/usr/bin/env python3
"""
使用OpenMM生成优化好的SWM4-NDP水体系用于测试
"""

import numpy as np
import os
import pickle

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    print("警告: OpenMM未安装，将使用简化的系统生成")
    HAS_OPENMM = False

def create_water_box_openmm(n_waters, output_prefix="water_system"):
    """
    使用OpenMM创建并优化SWM4-NDP水盒子
    """
    if not HAS_OPENMM:
        print("OpenMM未安装，无法生成优化的水体系")
        return None
    
    print(f"\n创建 {n_waters} 个水分子的系统...")
    
    # 创建拓扑
    topology = app.Topology()
    positions = []
    
    # 使用初始密度稍低于水的密度，让NPT自然平衡
    # 初始密度约 0.9 g/mL = 0.9 g/cm³
    initial_density = 0.9 * unit.gram / unit.centimeter**3
    mass_per_water = 18.015 * unit.gram / unit.mole
    total_mass = n_waters * mass_per_water / unit.AVOGADRO_CONSTANT_NA
    volume = total_mass / initial_density
    box_length = (volume ** (1.0/3.0)).in_units_of(unit.nanometer)
    
    # 添加周期性盒子
    topology.setPeriodicBoxVectors([
        [box_length, 0, 0],
        [0, box_length, 0],
        [0, 0, box_length]
    ] * unit.nanometer)
    
    # 创建链
    chain = topology.addChain()
    
    # 添加水分子
    spacing = box_length / (n_waters ** (1.0/3.0))
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 添加残基
                residue = topology.addResidue(f'WAT', chain)
                
                # 原子位置（初始随机扰动）
                x = (i + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                y = (j + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                z = (k + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                
                # 添加原子 (O, D, H1, H2, M for SWM4-NDP)
                o_atom = topology.addAtom('O', mm.Element.getBySymbol('O'), residue)
                d_atom = topology.addAtom('D', mm.Element.getBySymbol('H'), residue)  # Drude
                h1_atom = topology.addAtom('H1', mm.Element.getBySymbol('H'), residue)
                h2_atom = topology.addAtom('H2', mm.Element.getBySymbol('H'), residue)
                m_atom = topology.addAtom('M', None, residue)  # Virtual site
                
                # 位置（简化的TIP4P几何）
                positions.extend([
                    [x, y, z],  # O
                    [x, y, z],  # D (初始与O重合)
                    [x + 0.09572, y, z + 0.03],  # H1
                    [x - 0.04786, y + 0.08288, z + 0.03],  # H2
                    [x, y - 0.024034, z]  # M
                ])
                
                water_count += 1
    
    positions = positions * unit.nanometer
    
    # 创建系统
    print("创建OpenMM系统...")
    forcefield = app.ForceField()
    
    # 创建自定义的SWM4-NDP力场
    system = mm.System()
    
    # 添加粒子
    for i in range(topology.getNumAtoms()):
        if i % 5 == 0:  # O
            system.addParticle(15.999 * unit.amu)
        elif i % 5 == 1:  # D
            system.addParticle(0.4 * unit.amu)
        elif i % 5 in [2, 3]:  # H
            system.addParticle(1.008 * unit.amu)
        else:  # M
            system.addParticle(0.0 * unit.amu)
    
    # 添加Drude力
    drude_force = mm.DrudeForce()
    
    # SWM4-NDP参数
    charge_O = 1.71636 * unit.elementary_charge
    charge_D = -1.71636 * unit.elementary_charge
    charge_H = 0.55733 * unit.elementary_charge
    charge_M = -1.11466 * unit.elementary_charge
    
    k_spring = 418400.0 * unit.kilojoule_per_mole / unit.nanometer**2
    polarizability = 0.0009782237 * unit.nanometer**3
    
    # 添加Drude粒子
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        drude_force.addParticle(d_idx, o_idx, -1, -1, -1, 
                               charge_D, polarizability, 1.0, 1.0)
    
    # 添加Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drude_force.addScreenedPair(i, j, 1.3)
    
    system.addForce(drude_force)
    
    # 添加非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(1.2 * unit.nanometer)
    
    # 设置电荷和LJ参数
    for i in range(topology.getNumAtoms()):
        if i % 5 == 0:  # O
            nonbonded.addParticle(charge_O, 0.318395 * unit.nanometer, 0.88257 * unit.kilojoule_per_mole)
        elif i % 5 == 1:  # D
            nonbonded.addParticle(charge_D, 0.0 * unit.nanometer, 0.0 * unit.kilojoule_per_mole)
        elif i % 5 in [2, 3]:  # H
            nonbonded.addParticle(charge_H, 0.0 * unit.nanometer, 0.0 * unit.kilojoule_per_mole)
        else:  # M
            nonbonded.addParticle(charge_M, 0.0 * unit.nanometer, 0.0 * unit.kilojoule_per_mole)
    
    # 排除分子内相互作用
    for i in range(n_waters):
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base + j, base + k, 0, 0, 0)
    
    system.addForce(nonbonded)
    
    # 创建积分器（NPT）
    print("设置NPT模拟...")
    temperature = 300 * unit.kelvin
    pressure = 1 * unit.bar
    
    integrator = mm.DrudeLangevinIntegrator(
        temperature,
        1.0 / unit.picosecond,  # friction coefficient
        1.0 * unit.kelvin,      # Drude temperature
        20.0 / unit.picosecond, # Drude friction
        0.001 * unit.picosecond # timestep
    )
    
    barostat = mm.MonteCarloBarostat(pressure, temperature)
    system.addForce(barostat)
    
    # 创建模拟
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("运行能量最小化...")
    simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole)
    
    # NPT平衡 - 更长时间让盒子达到平衡
    print("运行NPT平衡...")
    
    # 先快速平衡
    print("  快速平衡 (5000步)...")
    simulation.step(5000)
    
    # 再精细平衡
    print("  精细平衡 (10000步)...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(10000)
    
    # 生产运行收集平均盒子大小
    print("  生产运行 (5000步)...")
    simulation.step(5000)
    
    # 获取最终位置和盒子
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    final_positions = state.getPositions()
    final_box = state.getPeriodicBoxVectors()
    
    # 保存结果
    result = {
        'n_waters': n_waters,
        'positions': final_positions,
        'box_vectors': final_box,
        'topology': topology
    }
    
    # 保存为pickle文件
    output_file = f"{output_prefix}_{n_waters}.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump(result, f)
    
    print(f"系统已保存到: {output_file}")
    
    # 同时保存为PDB
    pdb_file = f"{output_prefix}_{n_waters}.pdb"
    app.PDBFile.writeFile(topology, final_positions, open(pdb_file, 'w'))
    print(f"PDB文件已保存到: {pdb_file}")
    
    return result

def convert_openmm_to_pygcmc(openmm_data):
    """
    将OpenMM数据转换为pygcmc格式
    """
    import pygcmc
    
    n_waters = openmm_data['n_waters']
    positions = openmm_data['positions']
    box_vectors = openmm_data['box_vectors']
    
    # 创建pygcmc状态
    atoms = []
    residues = []
    
    # SWM4-NDP电荷
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters):
        for j in range(5):
            atom = pygcmc.MCAtom()
            idx = i * 5 + j
            atom.x = positions[idx][0].value_in_unit(unit.nanometer)
            atom.y = positions[idx][1].value_in_unit(unit.nanometer)
            atom.z = positions[idx][2].value_in_unit(unit.nanometer)
            atom.charge = charges[j]
            atom.type = types[j]
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    # 设置盒子
    box_x = box_vectors[0][0].value_in_unit(unit.nanometer)
    box_y = box_vectors[1][1].value_in_unit(unit.nanometer)
    box_z = box_vectors[2][2].value_in_unit(unit.nanometer)
    state.info.box = np.array([box_x, box_y, box_z])
    state.info.cutoff = min(1.2, min(box_x, box_y, box_z) / 2 - 0.01)
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def generate_test_systems():
    """
    生成一系列测试系统
    """
    # 2的幂次，方便测试扩展性
    # 先测试小系统
    system_sizes = [2, 4, 8, 16, 32]  # 暂时不生成大系统
    
    if not HAS_OPENMM:
        print("OpenMM未安装，使用简化方法生成系统")
        return generate_simple_systems(system_sizes)
    
    print("使用OpenMM生成优化的水体系")
    print("="*60)
    
    # 创建输出目录
    os.makedirs("optimized_systems", exist_ok=True)
    os.chdir("optimized_systems")
    
    for n_waters in system_sizes:
        try:
            create_water_box_openmm(n_waters)
        except Exception as e:
            print(f"生成 {n_waters} 水系统时出错: {e}")
    
    print("\n所有系统生成完成!")

def generate_simple_systems(system_sizes):
    """
    不使用OpenMM的简化系统生成
    """
    import pygcmc
    
    print("生成简化的水体系（未优化）")
    print("="*60)
    
    os.makedirs("simple_systems", exist_ok=True)
    
    for n_waters in system_sizes:
        print(f"\n创建 {n_waters} 个水分子的系统...")
        
        # 基于密度估算盒子大小
        density = 1000  # kg/m^3
        mass_per_water = 18.015e-3  # kg/mol
        avogadro = 6.022e23
        volume_per_water = mass_per_water / (density * avogadro) * 1e27  # nm^3
        total_volume = n_waters * volume_per_water
        box_length = total_volume ** (1.0/3.0)
        
        # 创建系统
        spacing = box_length / (n_waters ** (1.0/3.0))
        n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
        
        atoms = []
        residues = []
        
        water_count = 0
        for i in range(n_per_side):
            for j in range(n_per_side):
                for k in range(n_per_side):
                    if water_count >= n_waters:
                        break
                    
                    # 添加随机扰动
                    x = (i + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                    y = (j + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                    z = (k + 0.5 + 0.1 * (np.random.rand() - 0.5)) * spacing
                    
                    # SWM4-NDP水模型
                    positions = [
                        (x, y, z, 1.71636, 0),   # O
                        (x, y, z, -1.71636, 1),  # D
                        (x + 0.09572, y, z + 0.03, 0.55733, 2),  # H1
                        (x - 0.04786, y + 0.08288, z + 0.03, 0.55733, 2),  # H2
                        (x, y - 0.024034, z, -1.11466, 3)  # M
                    ]
                    
                    for px, py, pz, charge, typ in positions:
                        atom = pygcmc.MCAtom()
                        atom.x = px
                        atom.y = py
                        atom.z = pz
                        atom.charge = charge
                        atom.type = typ
                        atoms.append(atom)
                    
                    res = pygcmc.MCResidue()
                    res.atomStart = 5 * water_count
                    res.atomCount = 5
                    res.active = True
                    res.type = 0
                    residues.append(res)
                    
                    water_count += 1
        
        state = pygcmc.MCState()
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = len(atoms)
        state.activeResidueCount = n_waters
        
        state.info.box = np.array([box_length, box_length, box_length])
        state.info.cutoff = min(1.2, box_length / 2 - 0.01)
        
        # SWM4-NDP力场参数
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 保存
        result = {
            'n_waters': n_waters,
            'state': state,
            'box_length': box_length
        }
        
        output_file = f"simple_systems/water_{n_waters}.pkl"
        with open(output_file, 'wb') as f:
            pickle.dump(result, f)
        
        print(f"系统已保存到: {output_file}")
    
    print("\n注意: 这些系统未经优化，可能有较高的初始能量")

if __name__ == "__main__":
    generate_test_systems()