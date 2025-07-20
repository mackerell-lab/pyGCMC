#!/usr/bin/env python3
"""
生成密度恰好为1.0 g/cm³的水分子系统
"""

import numpy as np
import pickle
import os

def calculate_box_length_for_density(n_waters, target_density=1.0):
    """
    计算给定水分子数和目标密度下的盒子长度
    
    target_density: g/cm³
    返回: 盒子长度 (nm)
    """
    # 水的摩尔质量
    molar_mass_water = 18.01528  # g/mol
    avogadro = 6.02214076e23     # mol⁻¹
    
    # 总质量
    total_mass_g = n_waters * molar_mass_water / avogadro  # g
    
    # 目标体积
    target_volume_cm3 = total_mass_g / target_density  # cm³
    target_volume_nm3 = target_volume_cm3 * 1e21      # nm³
    
    # 立方体盒子边长
    box_length = target_volume_nm3 ** (1.0/3.0)  # nm
    
    return box_length

def verify_density(n_waters, box_length):
    """
    验证实际密度
    """
    molar_mass_water = 18.01528  # g/mol
    avogadro = 6.02214076e23     # mol⁻¹
    
    total_mass_g = n_waters * molar_mass_water / avogadro
    volume_nm3 = box_length ** 3
    volume_cm3 = volume_nm3 * 1e-21
    actual_density = total_mass_g / volume_cm3
    
    return actual_density

def create_water_system_exact_density(n_waters, target_density=1.0):
    """
    创建精确密度的水分子系统
    """
    # 计算精确的盒子长度
    box_length = calculate_box_length_for_density(n_waters, target_density)
    
    # 验证密度
    actual_density = verify_density(n_waters, box_length)
    
    print(f"\n{n_waters} 水分子系统:")
    print(f"  目标密度: {target_density:.6f} g/cm³")
    print(f"  实际密度: {actual_density:.6f} g/cm³")
    print(f"  盒子长度: {box_length:.6f} nm")
    print(f"  误差: {abs(actual_density - target_density):.2e} g/cm³")
    
    # 在立方格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    water_count = 0
    
    # 设置随机种子以确保可重复性
    np.random.seed(42)
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础位置（格子点）
                base_x = (i + 0.5) * spacing
                base_y = (j + 0.5) * spacing
                base_z = (k + 0.5) * spacing
                
                # 添加小的随机扰动避免完全对称
                # 扰动幅度不超过spacing的10%
                dx = (np.random.random() - 0.5) * spacing * 0.1
                dy = (np.random.random() - 0.5) * spacing * 0.1
                dz = (np.random.random() - 0.5) * spacing * 0.1
                
                x = base_x + dx
                y = base_y + dy
                z = base_z + dz
                
                # 确保在盒子内（应用周期性边界）
                x = x % box_length
                y = y % box_length
                z = z % box_length
                
                # 生成随机旋转的水分子
                # 使用欧拉角
                alpha = np.random.random() * 2 * np.pi
                beta = np.random.random() * np.pi
                gamma = np.random.random() * 2 * np.pi
                
                # 创建旋转矩阵
                R = euler_to_rotation_matrix(alpha, beta, gamma)
                
                # SWM4-NDP水分子几何参数
                r_oh = 0.09572  # O-H键长 (nm)
                angle_hoh = 104.52 * np.pi / 180.0  # H-O-H键角
                
                # 标准水分子坐标（O在原点）
                h1_std = np.array([
                    r_oh * np.sin(angle_hoh/2),
                    0.0,
                    r_oh * np.cos(angle_hoh/2)
                ])
                
                h2_std = np.array([
                    -r_oh * np.sin(angle_hoh/2),
                    0.0,
                    r_oh * np.cos(angle_hoh/2)
                ])
                
                # 应用旋转
                h1_rot = R @ h1_std
                h2_rot = R @ h2_std
                
                # SWM4-NDP的5个位点
                # 1. O原子
                positions.append([x, y, z])
                
                # 2. Drude粒子（初始时与O重合）
                positions.append([x, y, z])
                
                # 3. H1原子
                positions.append([
                    x + h1_rot[0],
                    y + h1_rot[1],
                    z + h1_rot[2]
                ])
                
                # 4. H2原子
                positions.append([
                    x + h2_rot[0],
                    y + h2_rot[1],
                    z + h2_rot[2]
                ])
                
                # 5. M位点（虚拟位点）
                # M位于O沿HOH角平分线反方向0.024034 nm处
                bisector = -(h1_rot + h2_rot)
                bisector_norm = bisector / np.linalg.norm(bisector)
                m_pos = bisector_norm * 0.024034
                
                positions.append([
                    x + m_pos[0],
                    y + m_pos[1],
                    z + m_pos[2]
                ])
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return positions, box_length, actual_density

def euler_to_rotation_matrix(alpha, beta, gamma):
    """
    将欧拉角转换为旋转矩阵
    使用Z-Y-Z约定
    """
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)
    cg = np.cos(gamma)
    sg = np.sin(gamma)
    
    R = np.array([
        [ca*cb*cg - sa*sg, -ca*cb*sg - sa*cg, ca*sb],
        [sa*cb*cg + ca*sg, -sa*cb*sg + ca*cg, sa*sb],
        [-sb*cg, sb*sg, cb]
    ])
    
    return R

def save_system(n_waters, positions, box_length, density, output_dir):
    """
    保存系统到文件
    """
    data = {
        'n_waters': n_waters,
        'positions': positions,
        'box_length': box_length,
        'density': density
    }
    
    # 保存pickle文件
    pickle_file = os.path.join(output_dir, f'water_{n_waters}.pkl')
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    
    # 保存PDB文件
    pdb_file = os.path.join(output_dir, f'water_{n_waters}.pdb')
    with open(pdb_file, 'w') as f:
        # PDB头部，包含晶胞信息
        f.write(f"CRYST1{box_length*10:9.3f}{box_length*10:9.3f}{box_length*10:9.3f}  90.00  90.00  90.00 P 1           1\n")
        f.write(f"REMARK   Density = {density:.6f} g/cm³\n")
        f.write(f"REMARK   SWM4-NDP water model with 5 sites per molecule\n")
        
        atom_idx = 1
        for i in range(n_waters):
            base_idx = i * 5
            
            # O原子
            f.write(f"ATOM  {atom_idx:5d}  O   WAT A{i+1:4d}    ")
            f.write(f"{positions[base_idx][0]*10:8.3f}")
            f.write(f"{positions[base_idx][1]*10:8.3f}")
            f.write(f"{positions[base_idx][2]*10:8.3f}")
            f.write(f"  1.00  0.00           O\n")
            atom_idx += 1
            
            # Drude粒子（注释形式）
            f.write(f"REMARK  Drude particle {i+1} at ")
            f.write(f"({positions[base_idx+1][0]*10:.3f}, ")
            f.write(f"{positions[base_idx+1][1]*10:.3f}, ")
            f.write(f"{positions[base_idx+1][2]*10:.3f})\n")
            
            # H1原子
            f.write(f"ATOM  {atom_idx:5d}  H1  WAT A{i+1:4d}    ")
            f.write(f"{positions[base_idx+2][0]*10:8.3f}")
            f.write(f"{positions[base_idx+2][1]*10:8.3f}")
            f.write(f"{positions[base_idx+2][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # H2原子
            f.write(f"ATOM  {atom_idx:5d}  H2  WAT A{i+1:4d}    ")
            f.write(f"{positions[base_idx+3][0]*10:8.3f}")
            f.write(f"{positions[base_idx+3][1]*10:8.3f}")
            f.write(f"{positions[base_idx+3][2]*10:8.3f}")
            f.write(f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # M位点（注释形式）
            f.write(f"REMARK  M-site {i+1} at ")
            f.write(f"({positions[base_idx+4][0]*10:.3f}, ")
            f.write(f"{positions[base_idx+4][1]*10:.3f}, ")
            f.write(f"{positions[base_idx+4][2]*10:.3f})\n")
        
        f.write("END\n")
    
    print(f"  保存到: {pickle_file}")
    print(f"  PDB文件: {pdb_file}")

def main():
    """
    主函数
    """
    print("生成密度精确为1.0 g/cm³的SWM4-NDP水分子系统")
    print("="*70)
    
    # 创建输出目录
    output_dir = '../tests/performance/water_density_1.0'
    os.makedirs(output_dir, exist_ok=True)
    
    # 要生成的系统大小
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    
    # 汇总表头
    print(f"\n{'水分子数':>8} {'盒子长度(nm)':>12} {'密度(g/cm³)':>12} {'密度误差':>12}")
    print("-"*50)
    
    for n_waters in system_sizes:
        # 生成系统
        positions, box_length, density = create_water_system_exact_density(n_waters, 1.0)
        
        # 保存系统
        save_system(n_waters, positions, box_length, density, output_dir)
        
        # 汇总输出
        print(f"{n_waters:8d} {box_length:12.6f} {density:12.6f} {abs(density-1.0):12.2e}")
    
    print("\n完成！所有系统的密度都精确为1.0 g/cm³")
    print("\n物理常数说明:")
    print(f"  水的摩尔质量: 18.01528 g/mol")
    print(f"  阿伏伽德罗常数: 6.02214076e23 mol⁻¹")
    print(f"  密度计算: ρ = (n_waters × M_water) / (N_A × V)")
    
    # 额外验证
    print("\n密度验证（以256水为例）:")
    n = 256
    box = calculate_box_length_for_density(n, 1.0)
    mass = n * 18.01528 / 6.02214076e23  # g
    volume = box**3 * 1e-21  # cm³
    print(f"  质量: {mass:.6e} g")
    print(f"  体积: {volume:.6e} cm³")
    print(f"  密度: {mass/volume:.6f} g/cm³")

if __name__ == "__main__":
    main()