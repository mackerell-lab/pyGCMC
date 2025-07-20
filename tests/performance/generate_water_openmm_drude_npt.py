#!/usr/bin/env python3
"""
使用OpenMM生成优化的SWM4-NDP水分子体系
通过NPT模拟优化密度和结构
"""

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import pickle
import os

def create_swm4ndp_system(n_waters):
    """
    创建SWM4-NDP水分子体系
    """
    print(f"\n创建 {n_waters} 个SWM4-NDP水分子...")
    
    # 创建拓扑
    topology = app.Topology()
    chain = topology.addChain()
    
    # 添加水分子
    positions = []
    
    # 估算盒子大小 (初始密度约0.95 g/cm³)
    volume_per_water = 30.0  # Å³
    total_volume = n_waters * volume_per_water
    box_length = (total_volume ** (1.0/3.0)) * unit.angstrom
    
    # 在格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                # 添加残基
                residue = topology.addResidue('HOH', chain)
                
                # 水分子中心位置
                x = i * spacing
                y = j * spacing  
                z = k * spacing
                
                # 添加原子和位置
                # O原子
                o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
                positions.append([x, y, z])
                
                # H原子 (简单的几何结构)
                h1_atom = topology.addAtom('H1', app.Element.getBySymbol('H'), residue)
                positions.append([x + 0.0756*unit.nanometer, y + 0.0586*unit.nanometer, z])
                
                h2_atom = topology.addAtom('H2', app.Element.getBySymbol('H'), residue)
                positions.append([x - 0.0756*unit.nanometer, y + 0.0586*unit.nanometer, z])
                
                # 添加键
                topology.addBond(o_atom, h1_atom)
                topology.addBond(o_atom, h2_atom)
                
                water_count += 1
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    # 设置周期性盒子
    topology.setPeriodicBoxVectors([
        [box_length.value_in_unit(unit.nanometer), 0, 0],
        [0, box_length.value_in_unit(unit.nanometer), 0],
        [0, 0, box_length.value_in_unit(unit.nanometer)]
    ] * unit.nanometer)
    
    positions = positions * unit.nanometer
    
    return topology, positions, box_length

def setup_drude_system(topology, positions):
    """
    设置Drude极化力场
    """
    # 创建力场
    forcefield = app.ForceField()
    
    # 创建系统
    system = forcefield.createSystem(topology)
    
    # 添加Drude力
    drude = mm.DrudeForce()
    
    # SWM4-NDP参数
    # 每个水分子有5个粒子: O, H1, H2, M (虚拟位点), D (Drude)
    # 这里简化处理，只考虑O原子的Drude粒子
    
    # 添加非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(1.2*unit.nanometer)
    
    # 参数来自SWM4-NDP模型
    # O: sigma=0.318395 nm, epsilon=0.88257 kJ/mol, charge=1.71636 e
    # H: charge=0.55733 e
    # M: charge=-1.11466 e
    # D: charge=-1.71636 e, polarizability=0.0009782237 nm³
    
    # 为每个原子添加参数
    n_atoms = topology.getNumAtoms()
    atom_index = 0
    
    for residue in topology.residues():
        # O原子
        nonbonded.addParticle(1.71636*unit.elementary_charge, 
                            0.318395*unit.nanometer, 
                            0.88257*unit.kilojoule_per_mole)
        
        # 添加Drude粒子到O原子
        k_spring = 418400.0 * unit.kilojoule_per_mole / unit.nanometer**2
        drude_charge = -1.71636 * unit.elementary_charge
        drude_index = system.addParticle(0.4 * unit.amu)  # Drude粒子质量
        drude.addParticle(drude_index, atom_index, -1, -1, -1, -1,
                         drude_charge, 0.0009782237*unit.nanometer**3, 1.0, 1.0)
        
        # H原子
        nonbonded.addParticle(0.55733*unit.elementary_charge, 
                            0.0*unit.nanometer, 
                            0.0*unit.kilojoule_per_mole)
        nonbonded.addParticle(0.55733*unit.elementary_charge, 
                            0.0*unit.nanometer, 
                            0.0*unit.kilojoule_per_mole)
        
        atom_index += 3
    
    # 添加Thole屏蔽
    thole_cutoff = 0.8 * unit.nanometer
    residues = list(topology.residues())
    for i in range(len(residues)):
        for j in range(i+1, len(residues)):
            # 计算O-O距离
            o1_idx = i * 3
            o2_idx = j * 3
            
            dx = positions[o1_idx][0] - positions[o2_idx][0]
            dy = positions[o1_idx][1] - positions[o2_idx][1]
            dz = positions[o1_idx][2] - positions[o2_idx][2]
            
            # 应用PBC
            box = topology.getPeriodicBoxVectors()[0][0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = (dx*dx + dy*dy + dz*dz)**0.5
            
            if dist < thole_cutoff:
                drude.addScreenedPair(i, j, 1.3)
    
    system.addForce(nonbonded)
    system.addForce(drude)
    
    # 添加键约束
    bonds = mm.HarmonicBondForce()
    # O-H键: 0.09572 nm
    for bond in topology.bonds():
        bonds.addBond(bond[0].index, bond[1].index, 
                     0.09572*unit.nanometer, 
                     502416.0*unit.kilojoule_per_mole/unit.nanometer**2)
    system.addForce(bonds)
    
    # 添加角约束
    angles = mm.HarmonicAngleForce()
    # H-O-H角: 104.52度
    for residue in topology.residues():
        atoms = list(residue.atoms())
        if len(atoms) >= 3:
            angles.addAngle(atoms[1].index, atoms[0].index, atoms[2].index,
                          104.52*unit.degree, 
                          628.02*unit.kilojoule_per_mole/unit.radian**2)
    system.addForce(angles)
    
    return system

def run_npt_equilibration(topology, system, positions, n_waters, output_dir):
    """
    运行NPT平衡模拟
    """
    print(f"\n运行NPT平衡模拟...")
    
    # 创建积分器 - 使用Drude Langevin积分器
    temperature = 298.0 * unit.kelvin
    drude_temperature = 1.0 * unit.kelvin
    friction = 1.0 / unit.picosecond
    drude_friction = 50.0 / unit.picosecond
    timestep = 1.0 * unit.femtosecond
    
    integrator = mm.DrudeLangevinIntegrator(
        temperature, friction, drude_temperature, drude_friction, timestep
    )
    integrator.setMaxDrudeDistance(0.02*unit.nanometer)  # 硬墙约束
    
    # 添加压力控制
    pressure = 1.0 * unit.atmosphere
    barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
    system.addForce(barostat)
    
    # 创建模拟
    platform = mm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("  能量最小化...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # 平衡运行
    print("  NPT平衡 (100 ps)...")
    simulation.context.setVelocitiesToTemperature(temperature)
    
    # 记录轨迹
    simulation.reporters.append(
        app.StateDataReporter(f'{output_dir}/equilibration_{n_waters}.log', 
                            1000, step=True, time=True, 
                            potentialEnergy=True, temperature=True, 
                            density=True, volume=True)
    )
    
    # 运行平衡
    simulation.step(100000)  # 100 ps
    
    # 获取最终状态
    state = simulation.context.getState(getPositions=True, 
                                      enforcePeriodicBox=True)
    final_positions = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()
    
    # 生产运行
    print("  生产运行 (500 ps)...")
    simulation.reporters.clear()
    simulation.reporters.append(
        app.StateDataReporter(f'{output_dir}/production_{n_waters}.log', 
                            5000, step=True, time=True, 
                            potentialEnergy=True, temperature=True, 
                            density=True, volume=True)
    )
    
    # PDB轨迹输出
    simulation.reporters.append(
        app.PDBReporter(f'{output_dir}/trajectory_{n_waters}.pdb', 10000)
    )
    
    simulation.step(500000)  # 500 ps
    
    # 获取最终优化的结构
    state = simulation.context.getState(getPositions=True, 
                                      enforcePeriodicBox=True)
    final_positions = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()
    
    return final_positions, box_vectors

def save_optimized_system(n_waters, positions, box_vectors, output_dir):
    """
    保存优化后的系统
    """
    # 提取位置数据
    positions_nm = []
    for pos in positions:
        positions_nm.append([
            pos[0].value_in_unit(unit.nanometer),
            pos[1].value_in_unit(unit.nanometer),
            pos[2].value_in_unit(unit.nanometer)
        ])
    
    # 盒子长度
    box_length = box_vectors[0][0].value_in_unit(unit.nanometer)
    
    # 计算密度
    volume_nm3 = box_length ** 3
    mass_g = n_waters * 18.015 / 6.022e23
    volume_cm3 = volume_nm3 * 1e-21
    density = mass_g / volume_cm3
    
    print(f"\n最终系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {density:.3f} g/cm³")
    
    # 保存pickle文件
    data = {
        'n_waters': n_waters,
        'positions': positions_nm,
        'box_length': box_length,
        'density': density
    }
    
    pickle_file = f'{output_dir}/water_{n_waters}_optimized.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    print(f"  保存到: {pickle_file}")
    
    # 保存PDB文件
    pdb_file = f'{output_dir}/water_{n_waters}_optimized.pdb'
    with open(pdb_file, 'w') as f:
        f.write(f"CRYST1{box_length*10:9.3f}{box_length*10:9.3f}{box_length*10:9.3f}  90.00  90.00  90.00 P 1           1\n")
        
        atom_idx = 1
        for i in range(n_waters):
            # 只保存O和H原子 (跳过Drude粒子)
            base_idx = i * 3
            
            # O原子
            f.write(f"ATOM  {atom_idx:5d}  O   HOH A{i+1:4d}    ")
            f.write(f"{positions_nm[base_idx][0]*10:8.3f}")
            f.write(f"{positions_nm[base_idx][1]*10:8.3f}")
            f.write(f"{positions_nm[base_idx][2]*10:8.3f}")
            f.write(f"  1.00  0.00           O\n")
            atom_idx += 1
            
            # H1原子
            f.write(f"ATOM  {atom_idx:5d}  H1  HOH A{i+1:4d}    ")
            f.write(f"{positions_nm[base_idx+1][0]*10:8.3f}")
            f.write(f"{positions_nm[base_idx+1][1]*10:8.3f}")
            f.write(f"{positions_nm[base_idx+1][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # H2原子
            f.write(f"ATOM  {atom_idx:5d}  H2  HOH A{i+1:4d}    ")
            f.write(f"{positions_nm[base_idx+2][0]*10:8.3f}")
            f.write(f"{positions_nm[base_idx+2][1]*10:8.3f}")
            f.write(f"{positions_nm[base_idx+2][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
        
        f.write("END\n")
    print(f"  PDB文件: {pdb_file}")

def main():
    """
    主函数
    """
    print("OpenMM SWM4-NDP水模型NPT优化")
    print("="*60)
    
    # 创建输出目录
    output_dir = '../tests/performance/optimized_water_openmm'
    os.makedirs(output_dir, exist_ok=True)
    
    # 要生成的系统大小
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    
    for n_waters in system_sizes:
        try:
            print(f"\n{'='*60}")
            print(f"处理 {n_waters} 水分子系统")
            print(f"{'='*60}")
            
            # 创建初始系统
            topology, positions, box_length = create_swm4ndp_system(n_waters)
            
            # 设置Drude力场
            system = setup_drude_system(topology, positions)
            
            # 运行NPT平衡
            final_positions, box_vectors = run_npt_equilibration(
                topology, system, positions, n_waters, output_dir
            )
            
            # 保存优化的系统
            save_optimized_system(n_waters, final_positions, box_vectors, output_dir)
            
        except Exception as e:
            print(f"\n错误: 处理 {n_waters} 水分子失败")
            print(f"原因: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n\n完成所有系统的生成和优化!")

if __name__ == "__main__":
    main()