#!/usr/bin/env python3
"""
使用OpenMM生成SWM4-NDP水分子体系 - 简化版本
专注于正确的水模型实现
"""

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import pickle
import os

def create_water_box_with_modeller(n_waters):
    """
    使用OpenMM Modeller创建SWM4-NDP水盒子
    """
    print(f"\n创建 {n_waters} 个水分子的盒子...")
    
    # 创建单个水分子作为模板
    pdb = app.PDBFile('../tests/data/single_water.pdb')
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # 删除所有原子，从空系统开始
    modeller.delete(modeller.topology.atoms())
    
    # 计算合适的盒子大小 (目标密度 ~1.0 g/cm³)
    # 每个水分子体积约30 Å³
    volume_per_water = 30.0e-3  # nm³
    total_volume = n_waters * volume_per_water
    box_size = (total_volume ** (1.0/3.0)) * 1.1  # 稍微大一点的初始盒子
    
    # 添加溶剂
    forcefield = app.ForceField('amber14/tip3p.xml')  # 先用TIP3P创建拓扑
    modeller.addSolvent(forcefield, model='tip3p', 
                       numAdded=n_waters,
                       boxSize=mm.Vec3(box_size, box_size, box_size)*unit.nanometer)
    
    # 现在需要为SWM4-NDP添加额外粒子
    # SWM4-NDP每个水分子有5个位点: O, H1, H2, M(虚拟位点), D(Drude)
    
    print(f"  初始盒子大小: {box_size:.3f} nm")
    print(f"  水分子数: {modeller.topology.getNumResidues()}")
    
    return modeller

def setup_swm4ndp_forcefield():
    """
    创建SWM4-NDP力场参数
    """
    # SWM4-NDP参数
    params = {
        'sigma_O': 0.318395 * unit.nanometer,
        'epsilon_O': 0.88257 * unit.kilojoule_per_mole,
        'charge_O': 1.71636 * unit.elementary_charge,
        'charge_H': 0.55733 * unit.elementary_charge,
        'charge_M': -1.11466 * unit.elementary_charge,
        'charge_D': -1.71636 * unit.elementary_charge,
        'polarizability': 0.0009782237 * unit.nanometer**3,
        'k_spring': 418400.0 * unit.kilojoule_per_mole / unit.nanometer**2,
        'thole': 1.3,
        'drude_mass': 0.4 * unit.amu
    }
    
    return params

def create_simple_water_system(n_waters):
    """
    创建简单的水分子系统用于测试
    """
    print(f"\n创建简化的 {n_waters} 水分子系统...")
    
    # 计算盒子大小
    volume_per_water = 30.0e-3  # nm³ (对应密度 ~1.0 g/cm³)
    total_volume = n_waters * volume_per_water
    box_length = (total_volume ** (1.0/3.0))
    
    # 在立方格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    water_count = 0
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 水分子中心
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                # 添加随机旋转
                theta = np.random.random() * 2 * np.pi
                phi = np.random.random() * np.pi
                
                # O-H键长 0.09572 nm, H-O-H角 104.52°
                r_oh = 0.09572
                angle_hoh = 104.52 * np.pi / 180.0
                
                # 计算H原子位置
                h1_x = r_oh * np.sin(angle_hoh/2) * np.cos(theta)
                h1_y = r_oh * np.sin(angle_hoh/2) * np.sin(theta)
                h1_z = r_oh * np.cos(angle_hoh/2)
                
                h2_x = r_oh * np.sin(angle_hoh/2) * np.cos(theta + np.pi)
                h2_y = r_oh * np.sin(angle_hoh/2) * np.sin(theta + np.pi)
                h2_z = r_oh * np.cos(angle_hoh/2)
                
                # O原子位置
                positions.append([x, y, z])
                
                # Drude粒子位置 (初始时与O重合)
                positions.append([x, y, z])
                
                # H1原子位置
                positions.append([x + h1_x, y + h1_y, z + h1_z])
                
                # H2原子位置
                positions.append([x + h2_x, y + h2_y, z + h2_z])
                
                # M位点位置 (质心位置的修正)
                m_x = x - 0.024034 * (h1_x + h2_x)
                m_y = y - 0.024034 * (h1_y + h2_y)
                m_z = z - 0.024034 * (h1_z + h2_z)
                positions.append([m_x, m_y, m_z])
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return positions, box_length

def optimize_with_gradient_descent(positions, n_waters, box_length, max_steps=1000):
    """
    使用简单的梯度下降优化水分子位置
    """
    print(f"\n优化水分子构型...")
    
    # SWM4-NDP参数
    sigma_O = 0.318395  # nm
    epsilon_O = 0.88257  # kJ/mol
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]  # O, D, H1, H2, M
    
    positions = np.array(positions)
    n_particles = len(positions)
    
    step_size = 0.001  # nm
    
    for step in range(max_steps):
        forces = np.zeros_like(positions)
        total_energy = 0.0
        
        # 计算LJ相互作用 (只有O-O)
        for i in range(n_waters):
            o1_idx = i * 5  # O原子索引
            
            for j in range(i+1, n_waters):
                o2_idx = j * 5
                
                # 计算距离
                dr = positions[o2_idx] - positions[o1_idx]
                
                # 应用PBC
                dr = dr - box_length * np.round(dr / box_length)
                
                r = np.linalg.norm(dr)
                
                if r < 1.2:  # 截断距离
                    # LJ势能和力
                    r6 = (sigma_O / r) ** 6
                    r12 = r6 * r6
                    
                    energy = 4 * epsilon_O * (r12 - r6)
                    force_mag = 24 * epsilon_O / r * (2 * r12 - r6)
                    
                    force_vec = force_mag * dr / r
                    forces[o1_idx] += force_vec
                    forces[o2_idx] -= force_vec
                    
                    total_energy += energy
        
        # 计算静电相互作用
        for i in range(n_particles):
            if i % 5 == 1:  # 跳过Drude粒子
                continue
                
            for j in range(i+1, n_particles):
                if j % 5 == 1:  # 跳过Drude粒子
                    continue
                
                # 不同水分子间的相互作用
                if i // 5 != j // 5:
                    dr = positions[j] - positions[i]
                    dr = dr - box_length * np.round(dr / box_length)
                    r = np.linalg.norm(dr)
                    
                    if r < 1.2 and r > 0.1:
                        q_i = charges[i % 5]
                        q_j = charges[j % 5]
                        
                        # 库仑能量和力
                        k_e = 138.935  # kJ/mol·nm/e²
                        energy = k_e * q_i * q_j / r
                        force_mag = k_e * q_i * q_j / (r * r)
                        
                        force_vec = force_mag * dr / r
                        forces[i] -= force_vec
                        forces[j] += force_vec
                        
                        total_energy += energy
        
        # 保持分子内结构
        for i in range(n_waters):
            base = i * 5
            # 保持H和M相对于O的位置
            forces[base+2] = forces[base]  # H1跟随O
            forces[base+3] = forces[base]  # H2跟随O
            forces[base+4] = forces[base]  # M跟随O
        
        # 更新位置
        max_force = np.max(np.abs(forces))
        if max_force > 100:
            forces = forces * 100 / max_force
        
        positions -= step_size * forces
        
        if step % 100 == 0:
            print(f"  步骤 {step}: 能量 = {total_energy:.2f} kJ/mol")
        
        # 收敛检查
        if np.max(np.abs(forces)) < 0.1:
            print(f"  收敛于步骤 {step}")
            break
    
    return positions.tolist()

def save_water_system(n_waters, positions, box_length, output_dir):
    """
    保存水分子系统
    """
    # 计算密度
    volume_nm3 = box_length ** 3
    mass_g = n_waters * 18.015 / 6.022e23
    volume_cm3 = volume_nm3 * 1e-21
    density = mass_g / volume_cm3
    
    print(f"\n系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {density:.3f} g/cm³")
    
    # 保存pickle文件 (用于pygcmc测试)
    data = {
        'n_waters': n_waters,
        'positions': positions,
        'box_length': box_length,
        'density': density
    }
    
    pickle_file = f'{output_dir}/water_{n_waters}.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    print(f"  保存到: {pickle_file}")
    
    # 保存PDB文件
    pdb_file = f'{output_dir}/water_{n_waters}.pdb'
    with open(pdb_file, 'w') as f:
        f.write(f"CRYST1{box_length*10:9.3f}{box_length*10:9.3f}{box_length*10:9.3f}  90.00  90.00  90.00 P 1           1\n")
        
        atom_idx = 1
        for i in range(n_waters):
            base = i * 5
            
            # O原子
            f.write(f"ATOM  {atom_idx:5d}  O   HOH A{i+1:4d}    ")
            f.write(f"{positions[base][0]*10:8.3f}")
            f.write(f"{positions[base][1]*10:8.3f}")
            f.write(f"{positions[base][2]*10:8.3f}")
            f.write(f"  1.00  0.00           O\n")
            atom_idx += 1
            
            # Drude粒子 (作为注释)
            f.write(f"REMARK  Drude at ({positions[base+1][0]:.3f}, {positions[base+1][1]:.3f}, {positions[base+1][2]:.3f})\n")
            
            # H1原子
            f.write(f"ATOM  {atom_idx:5d}  H1  HOH A{i+1:4d}    ")
            f.write(f"{positions[base+2][0]*10:8.3f}")
            f.write(f"{positions[base+2][1]*10:8.3f}")
            f.write(f"{positions[base+2][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # H2原子
            f.write(f"ATOM  {atom_idx:5d}  H2  HOH A{i+1:4d}    ")
            f.write(f"{positions[base+3][0]*10:8.3f}")
            f.write(f"{positions[base+3][1]*10:8.3f}")
            f.write(f"{positions[base+3][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # M位点 (作为注释)
            f.write(f"REMARK  M-site at ({positions[base+4][0]:.3f}, {positions[base+4][1]:.3f}, {positions[base+4][2]:.3f})\n")
        
        f.write("END\n")
    
    print(f"  PDB文件: {pdb_file}")

def main():
    """
    主函数
    """
    print("生成SWM4-NDP水分子测试系统")
    print("="*60)
    
    # 创建输出目录
    output_dir = '../tests/performance/water_systems'
    os.makedirs(output_dir, exist_ok=True)
    
    # 要生成的系统大小
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256]
    
    for n_waters in system_sizes:
        try:
            print(f"\n{'='*60}")
            print(f"生成 {n_waters} 水分子系统")
            print(f"{'='*60}")
            
            # 创建初始水分子位置
            positions, box_length = create_simple_water_system(n_waters)
            
            # 简单优化
            if n_waters <= 64:
                positions = optimize_with_gradient_descent(
                    positions, n_waters, box_length, max_steps=500
                )
            else:
                print("  跳过优化 (系统太大)")
            
            # 保存系统
            save_water_system(n_waters, positions, box_length, output_dir)
            
        except Exception as e:
            print(f"\n错误: 生成 {n_waters} 水分子失败")
            print(f"原因: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n\n完成!")
    print("现在可以运行 test_cg_final_comparison.py 进行测试")

if __name__ == "__main__":
    main()