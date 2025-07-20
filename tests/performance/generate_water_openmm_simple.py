#!/usr/bin/env python3
"""
使用OpenMM生成NPT优化的水分子体系
使用内置的水模型进行测试
"""

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    print("错误: 需要安装OpenMM")
    print("请运行: conda install -c conda-forge openmm")
    exit(1)

import numpy as np
import pickle
import os

def create_water_box(n_waters):
    """
    创建水盒子
    """
    print(f"\n创建 {n_waters} 个水分子...")
    
    # 创建PDB对象
    pdb = app.PDBFile('../tests/data/single_water.pdb')
    
    # 创建Modeller
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # 删除原有分子
    modeller.delete(modeller.topology.atoms())
    
    # 计算盒子大小 (目标密度 ~1.0 g/cm³)
    # 确保盒子足够大以容纳截断距离
    min_box_size = 2.5  # nm (截断距离1.0nm的2.5倍)
    volume_per_water = 30.0e-3  # nm³
    total_volume = n_waters * volume_per_water
    box_size = max((total_volume ** (1.0/3.0)) * 1.2, min_box_size)
    
    # 使用TIP3P力场添加水
    forcefield = app.ForceField('tip3p.xml')
    modeller.addSolvent(forcefield, 
                       boxSize=(box_size, box_size, box_size))
    
    return modeller.topology, modeller.positions

def run_npt_simulation(topology, positions, n_waters):
    """
    运行NPT模拟优化结构
    """
    print("设置NPT模拟...")
    
    # 创建系统
    forcefield = app.ForceField('tip3p.xml')
    system = forcefield.createSystem(topology,
                                   nonbondedMethod=app.PME,
                                   nonbondedCutoff=1.0*unit.nanometer,
                                   constraints=app.HBonds)
    
    # NPT积分器
    temperature = 300*unit.kelvin
    pressure = 1*unit.bar
    integrator = mm.LangevinMiddleIntegrator(temperature, 
                                            1/unit.picosecond, 
                                            2*unit.femtoseconds)
    
    # 添加压力控制
    barostat = mm.MonteCarloBarostat(pressure, temperature)
    system.addForce(barostat)
    
    # 创建模拟
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
    # 能量最小化
    print("  能量最小化...")
    simulation.minimizeEnergy()
    
    # 平衡
    print("  NPT平衡 (10 ps)...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(5000)  # 10 ps
    
    # 生产运行
    print("  生产运行 (20 ps)...")
    simulation.step(10000)  # 20 ps
    
    # 获取最终状态
    state = simulation.context.getState(getPositions=True, 
                                      enforcePeriodicBox=True)
    final_positions = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()
    
    # 计算密度
    box_length = box_vectors[0][0].value_in_unit(unit.nanometer)
    volume_nm3 = box_length ** 3
    mass_g = n_waters * 18.015 / 6.022e23
    volume_cm3 = volume_nm3 * 1e-21
    density = mass_g / volume_cm3
    
    print(f"  最终密度: {density:.3f} g/cm³")
    print(f"  盒子大小: {box_length:.3f} nm")
    
    return final_positions, box_length

def convert_to_swm4ndp_format(positions, n_waters, box_length):
    """
    转换为SWM4-NDP格式（5个位点）
    """
    swm4_positions = []
    
    for i in range(n_waters):
        # 获取O和H原子位置
        o_idx = i * 3
        h1_idx = i * 3 + 1
        h2_idx = i * 3 + 2
        
        o_pos = [positions[o_idx][j].value_in_unit(unit.nanometer) for j in range(3)]
        h1_pos = [positions[h1_idx][j].value_in_unit(unit.nanometer) for j in range(3)]
        h2_pos = [positions[h2_idx][j].value_in_unit(unit.nanometer) for j in range(3)]
        
        # O原子
        swm4_positions.append(o_pos)
        
        # Drude粒子 (初始与O重合)
        swm4_positions.append(o_pos)
        
        # H1原子
        swm4_positions.append(h1_pos)
        
        # H2原子
        swm4_positions.append(h2_pos)
        
        # M位点 (质心修正)
        # 简化计算：M在O的相反方向
        h_center = [(h1_pos[j] + h2_pos[j])/2 for j in range(3)]
        m_vec = [o_pos[j] - h_center[j] for j in range(3)]
        m_norm = np.sqrt(sum(x*x for x in m_vec))
        if m_norm > 0:
            m_vec = [x/m_norm * 0.024034 for x in m_vec]
        m_pos = [o_pos[j] + m_vec[j] for j in range(3)]
        swm4_positions.append(m_pos)
    
    return swm4_positions

def save_system(n_waters, positions, box_length, output_dir):
    """
    保存系统
    """
    # 计算密度
    volume_nm3 = box_length ** 3
    mass_g = n_waters * 18.015 / 6.022e23
    volume_cm3 = volume_nm3 * 1e-21
    density = mass_g / volume_cm3
    
    data = {
        'n_waters': n_waters,
        'positions': positions,
        'box_length': box_length,
        'density': density
    }
    
    # 保存pickle
    pickle_file = f'{output_dir}/water_{n_waters}.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"\n保存系统:")
    print(f"  文件: {pickle_file}")
    print(f"  水分子: {n_waters}")
    print(f"  盒子: {box_length:.3f} nm")
    print(f"  密度: {density:.3f} g/cm³")

def main():
    """
    主函数
    """
    print("使用OpenMM生成NPT优化的水体系")
    print("="*60)
    
    # 检查单个水分子模板
    if not os.path.exists('../tests/data/single_water.pdb'):
        print("创建单个水分子模板...")
        with open('../tests/data/single_water.pdb', 'w') as f:
            f.write("HETATM    1  O   HOH     1       0.000   0.000   0.000  1.00  0.00           O\n")
            f.write("HETATM    2  H1  HOH     1       0.757   0.586   0.000  1.00  0.00           H\n")
            f.write("HETATM    3  H2  HOH     1      -0.757   0.586   0.000  1.00  0.00           H\n")
            f.write("END\n")
    
    # 创建输出目录
    output_dir = '../tests/performance/water_openmm'
    os.makedirs(output_dir, exist_ok=True)
    
    # 系统大小
    system_sizes = [2, 4, 8, 16, 32]  # 先测试小系统
    
    for n_waters in system_sizes:
        try:
            print(f"\n{'='*60}")
            print(f"生成 {n_waters} 水分子系统")
            
            # 创建水盒子
            topology, positions = create_water_box(n_waters)
            
            # NPT优化
            final_positions, box_length = run_npt_simulation(topology, positions, n_waters)
            
            # 转换为SWM4-NDP格式
            swm4_positions = convert_to_swm4ndp_format(final_positions, n_waters, box_length)
            
            # 保存
            save_system(n_waters, swm4_positions, box_length, output_dir)
            
        except Exception as e:
            print(f"\n错误: {n_waters} 水分子失败")
            print(f"原因: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n\n完成!")

if __name__ == "__main__":
    main()