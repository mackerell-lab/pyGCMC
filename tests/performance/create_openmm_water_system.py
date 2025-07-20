#!/usr/bin/env python3
"""
使用OpenMM创建并优化128水分子系统
然后用于测试我们的FBP算法
"""

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    print("OpenMM未安装，尝试简化版本...")

import numpy as np
import pickle

def create_water_box_openmm(n_waters=128):
    """使用OpenMM创建水盒子"""
    if not HAS_OPENMM:
        print("需要安装OpenMM: conda install -c conda-forge openmm")
        return None
    
    print(f"使用OpenMM创建{n_waters}水分子系统...")
    
    # 创建拓扑
    topology = app.Topology()
    positions = []
    
    # 添加水分子
    chain = topology.addChain()
    
    # 计算盒子大小 - 水密度1g/cm³，但需要更大的盒子以满足截断要求
    volume_per_water = 30.0  # Å³
    total_volume = n_waters * volume_per_water
    box_length = np.cbrt(total_volume) * 1.5  # Å，增大1.5倍确保截断距离小于盒子一半
    
    # 在立方体中均匀分布
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_length / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                residue = topology.addResidue('HOH', chain)
                
                # 添加原子
                o = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
                h1 = topology.addAtom('H1', app.Element.getBySymbol('H'), residue)
                h2 = topology.addAtom('H2', app.Element.getBySymbol('H'), residue)
                
                # 添加键
                topology.addBond(o, h1)
                topology.addBond(o, h2)
                
                # 位置 (转换为纳米)
                x = (i + 0.5) * spacing * 0.1  # Å to nm
                y = (j + 0.5) * spacing * 0.1
                z = (k + 0.5) * spacing * 0.1
                
                # 标准水分子几何
                positions.append([x, y, z])  # O
                positions.append([x + 0.09572, y, z])  # H1
                positions.append([x - 0.04786, y + 0.08288, z])  # H2
                
                water_count += 1
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    positions = positions * unit.nanometer
    
    # 设置周期性边界
    topology.setPeriodicBoxVectors([
        [box_length * 0.1, 0, 0],
        [0, box_length * 0.1, 0],
        [0, 0, box_length * 0.1]
    ] * unit.nanometer)
    
    print(f"创建完成: {n_waters}水, 盒子尺寸 {box_length * 0.1:.2f} nm")
    
    return topology, positions

def minimize_with_openmm(topology, positions):
    """使用OpenMM能量最小化"""
    if not HAS_OPENMM:
        return None
    
    print("\n使用OpenMM进行能量最小化...")
    
    # 创建系统 - 使用TIP3P水模型（更简单）
    forcefield = app.ForceField('tip3p.xml')
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds
    )
    
    # 创建积分器和模拟
    integrator = mm.LangevinIntegrator(
        300*unit.kelvin,
        1/unit.picosecond,
        0.002*unit.picoseconds
    )
    
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("开始最小化...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # 获取最小化后的位置
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    positions = state.getPositions()
    energy = state.getPotentialEnergy()
    
    print(f"最小化完成，能量: {energy}")
    
    return positions

def save_water_configuration(positions, n_waters, filename):
    """保存水分子配置"""
    print(f"\n保存配置到 {filename}...")
    
    # 转换为简单格式
    config = {
        'n_waters': n_waters,
        'positions': []
    }
    
    # 假设每3个原子是一个水分子(O, H1, H2)
    for i in range(n_waters):
        water = {
            'O': [positions[3*i].x, positions[3*i].y, positions[3*i].z],
            'H1': [positions[3*i+1].x, positions[3*i+1].y, positions[3*i+1].z],
            'H2': [positions[3*i+2].x, positions[3*i+2].y, positions[3*i+2].z]
        }
        config['positions'].append(water)
    
    with open(filename, 'wb') as f:
        pickle.dump(config, f)
    
    print("保存完成！")

def create_simple_water_config(n_waters=128):
    """不使用OpenMM，创建简单的预配置水系统"""
    print(f"创建简单的{n_waters}水分子配置...")
    
    # 使用较大间距确保稳定
    volume_per_water = 40.0  # Å³ (比实际密度稍大)
    total_volume = n_waters * volume_per_water
    box_length = np.cbrt(total_volume)  # Å
    
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_length / n_per_side
    
    config = {
        'n_waters': n_waters,
        'box_length': box_length * 0.1,  # nm
        'positions': []
    }
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 中心位置
                x = (i + 0.5) * spacing * 0.1  # nm
                y = (j + 0.5) * spacing * 0.1
                z = (k + 0.5) * spacing * 0.1
                
                # 添加随机旋转
                angle = np.random.rand() * 2 * np.pi
                
                # 水分子几何 (相对坐标)
                o_pos = [x, y, z]
                h1_rel = [0.09572, 0, 0]
                h2_rel = [-0.04786, 0.08288, 0]
                
                # 旋转H原子
                cos_a = np.cos(angle)
                sin_a = np.sin(angle)
                
                h1_x = h1_rel[0] * cos_a - h1_rel[1] * sin_a
                h1_y = h1_rel[0] * sin_a + h1_rel[1] * cos_a
                
                h2_x = h2_rel[0] * cos_a - h2_rel[1] * sin_a
                h2_y = h2_rel[0] * sin_a + h2_rel[1] * cos_a
                
                water = {
                    'O': o_pos,
                    'H1': [x + h1_x, y + h1_y, z],
                    'H2': [x + h2_x, y + h2_y, z]
                }
                
                config['positions'].append(water)
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return config

def save_simple_config(config, filename):
    """保存简单配置"""
    with open(filename, 'wb') as f:
        pickle.dump(config, f)
    print(f"配置已保存到: {filename}")

def main():
    """主函数"""
    n_waters = 128
    
    if HAS_OPENMM:
        # 使用OpenMM
        print("检测到OpenMM，使用OpenMM创建优化的水系统")
        
        # 创建系统
        topology, positions = create_water_box_openmm(n_waters)
        
        if topology and positions:
            # 最小化
            optimized_positions = minimize_with_openmm(topology, positions)
            
            if optimized_positions:
                # 保存
                save_water_configuration(
                    optimized_positions, 
                    n_waters, 
                    'water128_openmm_optimized.pkl'
                )
    else:
        # 不使用OpenMM
        print("OpenMM未安装，创建简单的预配置系统")
        config = create_simple_water_config(n_waters)
        save_simple_config(config, 'water128_simple.pkl')
        
        print("\n配置信息:")
        print(f"  水分子数: {config['n_waters']}")
        print(f"  盒子尺寸: {config['box_length']:.2f} nm")
        print(f"  平均间距: {config['box_length']/n_waters**(1/3):.3f} nm")

if __name__ == "__main__":
    main()