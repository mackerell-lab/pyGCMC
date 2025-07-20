#!/usr/bin/env python3
"""
使用OpenMM对SWM4-NDP水体系进行NPT优化
"""

import numpy as np
import os
import pickle
import time

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    print("错误: 需要安装OpenMM")
    print("conda install -c conda-forge openmm")
    exit(1)

def create_swm4ndp_system(topology, nonbondedMethod=mm.NonbondedForce.PME):
    """
    创建SWM4-NDP水模型的OpenMM系统
    """
    system = mm.System()
    
    # 添加粒子质量
    n_waters = topology.getNumResidues()
    for i in range(n_waters):
        # O, D, H1, H2, M
        system.addParticle(15.999 * unit.amu)  # O
        system.addParticle(0.4 * unit.amu)     # D (Drude)
        system.addParticle(1.008 * unit.amu)   # H1
        system.addParticle(1.008 * unit.amu)   # H2
        system.addParticle(0.0 * unit.amu)     # M (virtual site)
    
    # 创建Drude力
    drudeForce = mm.DrudeForce()
    
    # SWM4-NDP参数
    charge_O = 1.71636
    charge_D = -1.71636
    charge_H = 0.55733
    charge_M = -1.11466
    
    k_spring = 1005.0 * unit.kilocalorie_per_mole / unit.angstrom**2  # 418400 kJ/mol/nm^2
    polarizability = 0.97822e-3 * unit.nanometer**3
    thole = 1.3
    
    # 添加Drude粒子
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        # addParticle(drude, parent, aniso12, aniso34, aniso56, charge, polarizability, aniso12scale, aniso34scale)
        drudeForce.addParticle(d_idx, o_idx, -1, -1, -1, 
                              charge_D * unit.elementary_charge, 
                              polarizability, 1.0, 1.0)
    
    # 添加Thole屏蔽对
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drudeForce.addScreenedPair(i, j, thole)
    
    system.addForce(drudeForce)
    
    # 创建非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(nonbondedMethod)
    nonbonded.setCutoffDistance(1.2 * unit.nanometer)
    nonbonded.setUseDispersionCorrection(True)
    
    # LJ参数
    sigma_O = 3.18395 * unit.angstrom
    epsilon_O = 0.21094 * unit.kilocalorie_per_mole
    
    # 添加粒子
    for i in range(n_waters):
        # O
        nonbonded.addParticle(charge_O * unit.elementary_charge, sigma_O, epsilon_O)
        # D
        nonbonded.addParticle(charge_D * unit.elementary_charge, 0.0, 0.0)
        # H1
        nonbonded.addParticle(charge_H * unit.elementary_charge, 0.0, 0.0)
        # H2
        nonbonded.addParticle(charge_H * unit.elementary_charge, 0.0, 0.0)
        # M
        nonbonded.addParticle(charge_M * unit.elementary_charge, 0.0, 0.0)
    
    # 排除分子内相互作用
    for i in range(n_waters):
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base + j, base + k, 0.0, 1.0, 0.0)
    
    system.addForce(nonbonded)
    
    # 添加虚拟位点
    for i in range(n_waters):
        o_idx = i * 5
        h1_idx = i * 5 + 2
        h2_idx = i * 5 + 3
        m_idx = i * 5 + 4
        
        # M位点在O-H平分线上，距离O 0.024034 nm
        # 权重: wO = 0.786646558, wH = 0.106676721
        w_O = 0.786646558
        w_H = 0.106676721
        
        virtual_site = mm.ThreeParticleAverageSite(o_idx, h1_idx, h2_idx, w_O, w_H, w_H)
        system.setVirtualSite(m_idx, virtual_site)
    
    return system

def create_initial_positions(n_waters):
    """
    创建初始位置和拓扑
    """
    # 基于0.9 g/cm³的初始密度
    density = 0.9 * unit.gram / unit.centimeter**3
    mass_per_water = 18.015 * unit.gram / unit.mole
    total_mass = n_waters * mass_per_water / unit.AVOGADRO_CONSTANT_NA
    volume = total_mass / density
    box_length = (volume ** (1.0/3.0)).in_units_of(unit.nanometer)
    
    # 创建拓扑
    topology = app.Topology()
    chain = topology.addChain()
    
    # 计算格子排列
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    
    # 水分子几何
    oh_bond = 0.09572 * unit.nanometer
    hoh_angle = 104.52 * unit.degree
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 添加残基
                residue = topology.addResidue('WAT', chain)
                
                # O位置（添加随机扰动）
                x = (i + 0.5 + 0.2*(np.random.rand()-0.5)) * spacing
                y = (j + 0.5 + 0.2*(np.random.rand()-0.5)) * spacing
                z = (k + 0.5 + 0.2*(np.random.rand()-0.5)) * spacing
                
                # 随机旋转水分子
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                psi = np.random.rand() * 2 * np.pi
                
                # 构建旋转矩阵
                R = rotation_matrix(theta, phi, psi)
                
                # 初始H位置（相对于O）
                angle_rad = hoh_angle.value_in_unit(unit.radian)
                h1_rel = np.array([
                    oh_bond.value_in_unit(unit.nanometer) * np.sin(angle_rad/2),
                    0,
                    oh_bond.value_in_unit(unit.nanometer) * np.cos(angle_rad/2)
                ])
                h2_rel = np.array([
                    -oh_bond.value_in_unit(unit.nanometer) * np.sin(angle_rad/2),
                    0,
                    oh_bond.value_in_unit(unit.nanometer) * np.cos(angle_rad/2)
                ])
                
                # 应用旋转
                h1_rel = R @ h1_rel
                h2_rel = R @ h2_rel
                
                # M位置（沿平分线）
                m_direction = (h1_rel + h2_rel) / 2
                m_direction = m_direction / np.linalg.norm(m_direction)
                m_rel = m_direction * 0.024034  # nm
                
                # 添加原子到拓扑
                o_atom = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
                d_atom = topology.addAtom('D', app.Element.getBySymbol('H'), residue)
                h1_atom = topology.addAtom('H1', app.Element.getBySymbol('H'), residue)
                h2_atom = topology.addAtom('H2', app.Element.getBySymbol('H'), residue)
                m_atom = topology.addAtom('M', None, residue)
                
                # 添加位置
                o_pos = np.array([x.value_in_unit(unit.nanometer), 
                                 y.value_in_unit(unit.nanometer), 
                                 z.value_in_unit(unit.nanometer)])
                
                positions.extend([
                    o_pos,                    # O
                    o_pos,                    # D (初始与O重合)
                    o_pos + h1_rel,          # H1
                    o_pos + h2_rel,          # H2
                    o_pos + m_rel            # M
                ])
                
                water_count += 1
    
    # 设置周期性盒子
    topology.setPeriodicBoxVectors([
        [box_length, 0*unit.nanometer, 0*unit.nanometer],
        [0*unit.nanometer, box_length, 0*unit.nanometer],
        [0*unit.nanometer, 0*unit.nanometer, box_length]
    ])
    
    return topology, positions * unit.nanometer

def rotation_matrix(theta, phi, psi):
    """
    创建3D旋转矩阵
    """
    # Z-Y-Z欧拉角
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)
    
    R = np.array([
        [cos_theta*cos_psi - cos_phi*sin_theta*sin_psi, -cos_theta*sin_psi - cos_phi*sin_theta*cos_psi, sin_phi*sin_theta],
        [sin_theta*cos_psi + cos_phi*cos_theta*sin_psi, -sin_theta*sin_psi + cos_phi*cos_theta*cos_psi, -sin_phi*cos_theta],
        [sin_phi*sin_psi, sin_phi*cos_psi, cos_phi]
    ])
    
    return R

def optimize_water_system(n_waters, output_prefix="optimized_water"):
    """
    使用OpenMM NPT优化水体系
    """
    print(f"\n优化 {n_waters} 个水分子系统...")
    
    # 创建初始结构
    print("  创建初始结构...")
    topology, positions = create_initial_positions(n_waters)
    
    # 创建系统
    print("  创建OpenMM系统...")
    system = create_swm4ndp_system(topology)
    
    # 创建积分器
    temperature = 300 * unit.kelvin
    pressure = 1 * unit.bar
    timestep = 0.5 * unit.femtosecond
    
    # Drude积分器
    integrator = mm.DrudeLangevinIntegrator(
        temperature,
        1.0 / unit.picosecond,      # 摩擦系数
        1.0 * unit.kelvin,          # Drude温度
        20.0 / unit.picosecond,     # Drude摩擦系数
        timestep
    )
    integrator.setMaxDrudeDistance(0.02 * unit.nanometer)  # 硬墙约束
    
    # 添加压力控制
    barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
    system.addForce(barostat)
    
    # 创建模拟
    platform = mm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("  能量最小化...")
    print(f"    初始能量: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole, maxIterations=1000)
    print(f"    最小化后: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    
    # NPT平衡
    print("  NPT平衡...")
    
    # 阶段1：快速平衡（高温）
    print("    阶段1: 快速平衡 (400K, 2 ps)")
    simulation.context.setVelocitiesToTemperature(400*unit.kelvin)
    simulation.step(4000)  # 2 ps
    
    # 阶段2：降温
    print("    阶段2: 降温到300K (3 ps)")
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(6000)  # 3 ps
    
    # 阶段3：平衡密度
    print("    阶段3: 平衡密度 (10 ps)")
    simulation.reporters.append(app.StateDataReporter(
        f'{output_prefix}_{n_waters}_equilibration.log', 
        1000, step=True, potentialEnergy=True, temperature=True, 
        density=True, volume=True
    ))
    simulation.step(20000)  # 10 ps
    
    # 阶段4：生产运行
    print("    阶段4: 生产运行 (5 ps)")
    simulation.step(10000)  # 5 ps
    
    # 获取最终状态
    state = simulation.context.getState(getPositions=True, getVelocities=True, 
                                      getEnergy=True, enforcePeriodicBox=True)
    final_positions = state.getPositions()
    final_box = state.getPeriodicBoxVectors()
    final_energy = state.getPotentialEnergy()
    
    # 计算最终密度
    box_volume = (final_box[0][0] * final_box[1][1] * final_box[2][2]).in_units_of(unit.nanometer**3)
    mass = n_waters * 18.015 * unit.gram / unit.mole / unit.AVOGADRO_CONSTANT_NA
    density = (mass / box_volume).in_units_of(unit.gram/unit.centimeter**3)
    
    print(f"\n  最终结果:")
    print(f"    盒子: {final_box[0][0]:.3f} x {final_box[1][1]:.3f} x {final_box[2][2]:.3f}")
    print(f"    密度: {density:.3f}")
    print(f"    能量: {final_energy}")
    print(f"    能量/水: {final_energy/n_waters}")
    
    # 保存结果
    result = {
        'n_waters': n_waters,
        'topology': topology,
        'positions': final_positions,
        'box_vectors': final_box,
        'energy': final_energy,
        'density': density,
        'temperature': temperature,
        'pressure': pressure
    }
    
    # 保存pickle
    os.makedirs('optimized_systems', exist_ok=True)
    pickle_file = f'optimized_systems/{output_prefix}_{n_waters}.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(result, f)
    print(f"    保存到: {pickle_file}")
    
    # 保存PDB
    pdb_file = f'optimized_systems/{output_prefix}_{n_waters}.pdb'
    with open(pdb_file, 'w') as f:
        app.PDBFile.writeFile(topology, final_positions, f)
    print(f"    PDB文件: {pdb_file}")
    
    return result

def main():
    """
    主函数：生成一系列优化的水体系
    """
    print("使用OpenMM NPT优化SWM4-NDP水体系")
    print("="*80)
    
    # 要生成的系统
    system_sizes = [2, 4, 8, 16, 32]  # 先测试小系统
    
    # 测试小系统
    if True:  # 设为False来跳过测试
        print("\n先测试小系统...")
        result = optimize_water_system(4)
        print("\n测试成功！")
    
    # 生成所有系统
    print("\n开始生成所有系统...")
    for n_waters in system_sizes:
        try:
            start_time = time.time()
            optimize_water_system(n_waters)
            elapsed = time.time() - start_time
            print(f"  用时: {elapsed:.1f} 秒\n")
        except Exception as e:
            print(f"\n错误: 优化 {n_waters} 水系统失败")
            print(f"  原因: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n完成！优化的系统保存在 optimized_systems/ 目录")

if __name__ == "__main__":
    main()