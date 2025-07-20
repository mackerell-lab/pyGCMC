#!/usr/bin/env python3
"""
使用OpenMM对水系统进行简单的NVT优化
不使用Drude力场，只使用标准TIP3P水模型
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

def create_simple_water_system(positions, box_length, n_waters):
    """
    创建简单的TIP3P水系统（只有O和H原子）
    """
    print("\n创建OpenMM系统...")
    
    # 创建拓扑
    topology = app.Topology()
    chain = topology.addChain()
    
    # 从5位点水模型中提取3位点坐标
    new_positions = []
    
    for i in range(n_waters):
        residue = topology.addResidue('HOH', chain)
        
        # O原子
        o_pos = positions[i*5]
        o_atom = topology.addAtom('O', app.element.oxygen, residue)
        new_positions.append(o_pos)
        
        # H1原子
        h1_pos = positions[i*5+2]
        h1_atom = topology.addAtom('H1', app.element.hydrogen, residue)
        new_positions.append(h1_pos)
        
        # H2原子
        h2_pos = positions[i*5+3]
        h2_atom = topology.addAtom('H2', app.element.hydrogen, residue)
        new_positions.append(h2_pos)
        
        # 添加键
        topology.addBond(o_atom, h1_atom)
        topology.addBond(o_atom, h2_atom)
    
    new_positions = np.array(new_positions) * unit.nanometer
    
    # 设置周期性盒子
    topology.setPeriodicBoxVectors([
        [box_length, 0, 0],
        [0, box_length, 0],
        [0, 0, box_length]
    ] * unit.nanometer)
    
    # 创建力场
    forcefield = app.ForceField('tip3p.xml')
    
    # 创建系统
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=min(0.9, box_length/2 - 0.01)*unit.nanometer,
        constraints=app.HBonds
    )
    
    print(f"  系统创建完成：")
    print(f"  - {system.getNumParticles()} 个粒子")
    print(f"  - {n_waters} 个水分子")
    print(f"  - 使用TIP3P力场")
    
    return system, topology, new_positions

def run_nvt_equilibration(system, topology, positions, temperature=300*unit.kelvin, 
                         steps=10000, dt=2*unit.femtosecond):
    """
    运行NVT平衡
    """
    print(f"\n运行NVT平衡...")
    print(f"  温度: {temperature}")
    print(f"  步数: {steps}")
    print(f"  时间步长: {dt}")
    
    # 创建积分器
    integrator = mm.LangevinIntegrator(
        temperature,
        1/unit.picosecond,
        dt
    )
    
    # 创建模拟
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
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
    n_residues = topology.getNumResidues()
    print(f"  能量/水: {final_energy/n_residues:.1f} kJ/mol")
    
    return final_positions, energies

def convert_back_to_5site(positions_3site, n_waters, original_positions):
    """
    将3位点坐标转换回5位点格式
    保持D和M位点在原始位置
    """
    positions_5site = []
    
    for i in range(n_waters):
        # O原子（从优化后的位置）
        positions_5site.append(positions_3site[i*3])
        
        # D粒子（保持原始相对位置）
        o_new = positions_3site[i*3]
        o_old = original_positions[i*5]
        d_old = original_positions[i*5+1]
        d_new = d_old - o_old + o_new  # 平移D粒子
        positions_5site.append(d_new)
        
        # H1原子（从优化后的位置）
        positions_5site.append(positions_3site[i*3+1])
        
        # H2原子（从优化后的位置）
        positions_5site.append(positions_3site[i*3+2])
        
        # M位点（保持原始相对位置）
        m_old = original_positions[i*5+4]
        m_new = m_old - o_old + o_new  # 平移M位点
        positions_5site.append(m_new)
    
    return np.array(positions_5site)

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
        'method': 'OpenMM_NVT_TIP3P'
    }
    
    pickle_file = os.path.join(output_dir, f'water_{n_waters}_nvt_simple.pkl')
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"\n保存优化后的系统到: {pickle_file}")

def main():
    """
    主函数
    """
    print("使用OpenMM进行简单NVT优化")
    print("使用TIP3P水模型")
    print("="*70)
    
    # 选择要优化的系统
    import sys
    n_waters = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    
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
    system, topology, positions_3site = create_simple_water_system(positions, box_length, n_waters)
    
    # 运行NVT平衡
    try:
        final_positions_3site, energies = run_nvt_equilibration(
            system, topology, positions_3site,
            temperature=300*unit.kelvin,
            steps=20000,  # 40 ps
            dt=2*unit.femtosecond
        )
        
        # 转换回5位点格式
        final_positions_5site = convert_back_to_5site(
            final_positions_3site, n_waters, positions
        )
        
        # 保存结果
        output_dir = 'optimized_water_systems'
        os.makedirs(output_dir, exist_ok=True)
        
        save_optimized_system(n_waters, final_positions_5site, box_length, output_dir)
        
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