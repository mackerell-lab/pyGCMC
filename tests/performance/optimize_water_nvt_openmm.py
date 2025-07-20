#!/usr/bin/env python3
"""
使用OpenMM对水系统进行NVT优化
保持密度1.0 g/cm³不变
"""

import numpy as np
import pickle
import os
import time

try:
    import openmm as mm
    import openmm.app as app
    from openmm import unit
    print("OpenMM导入成功")
except ImportError:
    print("错误：需要安装OpenMM")
    print("请运行: conda install -c conda-forge openmm")
    exit(1)

def create_openmm_system(positions, box_length, n_waters):
    """
    创建OpenMM系统用于SWM4-NDP水模型
    """
    print("\n创建OpenMM系统...")
    
    # 创建系统
    system = mm.System()
    
    # SWM4-NDP参数
    # 质量 (g/mol -> amu)
    mass_O = 15.99943  # 氧
    mass_H = 1.007947  # 氢
    mass_D = 0.4       # Drude粒子质量
    mass_M = 0.0       # 虚拟位点
    
    # 添加粒子
    for i in range(n_waters):
        # O
        system.addParticle(mass_O * unit.amu)
        # D (Drude)
        system.addParticle(mass_D * unit.amu)
        # H1
        system.addParticle(mass_H * unit.amu)
        # H2
        system.addParticle(mass_H * unit.amu)
        # M (虚拟位点)
        system.addParticle(mass_M * unit.amu)
    
    # 设置盒子
    system.setDefaultPeriodicBoxVectors(
        [box_length, 0, 0] * unit.nanometer,
        [0, box_length, 0] * unit.nanometer,
        [0, 0, box_length] * unit.nanometer
    )
    
    # 1. 非键相互作用
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    # 截断距离必须小于半盒子
    cutoff = min(0.9, box_length/2 - 0.01) * unit.nanometer
    nonbonded.setCutoffDistance(cutoff)
    nonbonded.setEwaldErrorTolerance(0.0005)
    print(f"  使用截断距离: {cutoff}")
    
    # SWM4-NDP电荷和LJ参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]  # e
    sigma_O = 0.318395 * unit.nanometer
    epsilon_O = 0.88257 * unit.kilojoule_per_mole
    
    for i in range(n_waters):
        # O - 有LJ参数
        nonbonded.addParticle(
            charges[0] * unit.elementary_charge,
            sigma_O,
            epsilon_O
        )
        # D - 无LJ
        nonbonded.addParticle(
            charges[1] * unit.elementary_charge,
            1.0 * unit.nanometer,  # sigma=1表示无LJ
            0.0 * unit.kilojoule_per_mole
        )
        # H1 - 无LJ
        nonbonded.addParticle(
            charges[2] * unit.elementary_charge,
            1.0 * unit.nanometer,
            0.0 * unit.kilojoule_per_mole
        )
        # H2 - 无LJ
        nonbonded.addParticle(
            charges[3] * unit.elementary_charge,
            1.0 * unit.nanometer,
            0.0 * unit.kilojoule_per_mole
        )
        # M - 无LJ
        nonbonded.addParticle(
            charges[4] * unit.elementary_charge,
            1.0 * unit.nanometer,
            0.0 * unit.kilojoule_per_mole
        )
        
        # 添加分子内排除
        base = i * 5
        for j in range(5):
            for k in range(j+1, 5):
                nonbonded.addException(base+j, base+k, 0, 1, 0)
    
    system.addForce(nonbonded)
    
    # 2. Drude力
    drudeForce = mm.DrudeForce()
    
    # Drude参数
    k_drude = 418400.0 * unit.kilojoule_per_mole / unit.nanometer**2
    drude_charge = -1.71636 * unit.elementary_charge
    polarizability = 0.00097825258 * unit.nanometer**3
    
    for i in range(n_waters):
        parent_idx = i * 5      # O
        drude_idx = i * 5 + 1   # D
        
        # 添加Drude粒子对
        drudeForce.addParticle(
            drude_idx,      # particle (Drude)
            parent_idx,     # particle1 (parent)
            -1,             # particle2 (aniso)
            -1,             # particle3 (aniso)
            -1,             # particle4 (aniso)
            drude_charge,   # charge
            polarizability, # polarizability
            1.0,            # aniso12
            1.0             # aniso34
        )
    
    # 3. Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drudeForce.addScreenedPair(i, j, 1.3)  # Thole参数 = 1.3
    
    system.addForce(drudeForce)
    
    # 4. 约束（O-H和H-H距离）
    # 注意：OpenMM的DrudeForce自动处理Drude-parent约束
    
    print(f"  系统创建完成：")
    print(f"  - {system.getNumParticles()} 个粒子")
    print(f"  - {system.getNumForces()} 个力场")
    print(f"  - {n_waters} 个水分子")
    print(f"  - {n_waters * (n_waters-1) // 2} 个Thole对")
    
    return system

def run_nvt_equilibration(system, positions, box_length, temperature=300*unit.kelvin, 
                         steps=10000, dt=0.5*unit.femtosecond):
    """
    运行NVT平衡
    """
    print(f"\n运行NVT平衡...")
    print(f"  温度: {temperature}")
    print(f"  步数: {steps}")
    print(f"  时间步长: {dt}")
    
    # 创建积分器 - 使用DrudeLangevinIntegrator
    integrator = mm.DrudeLangevinIntegrator(
        temperature,      # 温度
        1/unit.picosecond,  # 摩擦系数
        1*unit.kelvin,    # Drude温度（低温）
        10/unit.picosecond,  # Drude摩擦系数
        dt                # 时间步长
    )
    integrator.setMaxDrudeDistance(0.02 * unit.nanometer)  # 硬墙约束
    
    # 创建拓扑（简单版本）
    topology = app.Topology()
    chain = topology.addChain()
    
    n_waters = len(positions) // 5
    for i in range(n_waters):
        residue = topology.addResidue('HOH', chain)
        # 简化：只添加可见原子
        o_atom = topology.addAtom('O', app.element.oxygen, residue)
        h1_atom = topology.addAtom('H1', app.element.hydrogen, residue)
        h2_atom = topology.addAtom('H2', app.element.hydrogen, residue)
    
    # 设置周期性盒子
    topology.setPeriodicBoxVectors([
        [box_length, 0, 0],
        [0, box_length, 0],
        [0, 0, box_length]
    ] * unit.nanometer)
    
    # 创建模拟
    simulation = app.Simulation(topology, system, integrator)
    
    # 设置初始位置
    simulation.context.setPositions(positions * unit.nanometer)
    
    # 能量最小化
    print("\n  能量最小化...")
    initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"    初始能量: {initial_energy}")
    
    simulation.minimizeEnergy(maxIterations=1000)
    
    minimized_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"    最小化后能量: {minimized_energy}")
    
    # NVT平衡
    print(f"\n  开始NVT模拟...")
    
    # 记录能量
    energies = []
    report_interval = 1000
    
    start_time = time.time()
    
    for i in range(0, steps, report_interval):
        simulation.step(report_interval)
        
        # 获取能量
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        energies.append(energy)
        
        # 进度报告
        if i % (report_interval * 10) == 0:
            elapsed = time.time() - start_time
            progress = (i + report_interval) / steps * 100
            print(f"    进度: {progress:.0f}%, 能量: {energy:.0f} kJ/mol, 时间: {elapsed:.1f}s")
    
    # 获取最终位置
    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    final_positions = final_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    final_energy = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    print(f"\n  NVT平衡完成!")
    print(f"  最终能量: {final_energy:.0f} kJ/mol")
    print(f"  能量/水: {final_energy/n_waters:.1f} kJ/mol")
    
    return final_positions, energies

def save_optimized_system(n_waters, positions, box_length, output_dir):
    """
    保存优化后的系统
    """
    data = {
        'n_waters': n_waters,
        'positions': positions,
        'box_length': box_length,
        'density': 1.0,  # g/cm³
        'charges': [1.71636, -1.71636, 0.55733, 0.55733, -1.11466],
        'atom_names': ['O', 'D', 'H1', 'H2', 'M'],
        'optimized': True,
        'method': 'OpenMM_NVT'
    }
    
    pickle_file = os.path.join(output_dir, f'water_{n_waters}_nvt.pkl')
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"\n保存优化后的系统到: {pickle_file}")

def main():
    """
    主函数
    """
    print("使用OpenMM进行NVT优化")
    print("保持密度1.0 g/cm³")
    print("="*70)
    
    # 选择要优化的系统
    n_waters = 256  # 从256水开始
    
    # 加载初始结构
    input_file = f'large_water_systems/water_{n_waters}.pkl'
    print(f"\n加载初始结构: {input_file}")
    
    with open(input_file, 'rb') as f:
        data = pickle.load(f)
    
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {data['density']:.3f} g/cm³")
    
    # 创建OpenMM系统
    system = create_openmm_system(positions, box_length, n_waters)
    
    # 运行NVT平衡
    try:
        final_positions, energies = run_nvt_equilibration(
            system, positions, box_length,
            temperature=300*unit.kelvin,
            steps=10000,  # 5 ps
            dt=0.5*unit.femtosecond
        )
        
        # 保存结果
        output_dir = 'optimized_water_systems'
        os.makedirs(output_dir, exist_ok=True)
        
        save_optimized_system(n_waters, final_positions, box_length, output_dir)
        
        # 分析能量收敛
        print("\n能量收敛分析:")
        print(f"  初始能量/水: {energies[0]/n_waters:.1f} kJ/mol")
        print(f"  最终能量/水: {energies[-1]/n_waters:.1f} kJ/mol")
        print(f"  能量变化: {(energies[-1]-energies[0])/n_waters:.1f} kJ/mol/水")
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()