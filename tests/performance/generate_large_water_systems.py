#!/usr/bin/env python3
"""
生成大型水系统用于Drude模型测试
密度：1.0 g/cm³
系统大小：256, 512, 1024, 2048, 4096 水分子
"""

import numpy as np
import pickle
import os
import time

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

def euler_to_rotation_matrix(alpha, beta, gamma):
    """
    将欧拉角转换为旋转矩阵
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

def create_large_water_system(n_waters, target_density=1.0):
    """
    创建大型水系统
    使用更高效的方法处理大系统
    """
    print(f"\n创建 {n_waters} 水分子系统...")
    start_time = time.time()
    
    # 计算盒子大小
    box_length = calculate_box_length_for_density(n_waters, target_density)
    
    # 验证密度
    molar_mass_water = 18.01528
    avogadro = 6.02214076e23
    total_mass_g = n_waters * molar_mass_water / avogadro
    volume_nm3 = box_length ** 3
    volume_cm3 = volume_nm3 * 1e-21
    actual_density = total_mass_g / volume_cm3
    
    print(f"  目标密度: {target_density:.6f} g/cm³")
    print(f"  实际密度: {actual_density:.6f} g/cm³")
    print(f"  盒子长度: {box_length:.6f} nm")
    print(f"  半盒子: {box_length/2:.6f} nm")
    print(f"  Thole(0.8nm)/半盒子: {0.8/(box_length/2):.3f}")
    
    # 在立方格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    print(f"  格子: {n_per_side}×{n_per_side}×{n_per_side}")
    print(f"  间距: {spacing:.3f} nm")
    
    # 预分配数组（5个原子/水）
    positions = np.zeros((n_waters * 5, 3))
    
    # 设置随机种子
    np.random.seed(42)
    
    # SWM4-NDP水分子几何参数
    r_oh = 0.09572  # O-H键长 (nm)
    angle_hoh = 104.52 * np.pi / 180.0  # H-O-H键角
    r_om = 0.024034  # O-M距离 (nm)
    
    water_count = 0
    
    # 生成水分子位置
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础位置（格子点）+ 小扰动
                base_pos = np.array([
                    (ix + 0.5) * spacing,
                    (iy + 0.5) * spacing,
                    (iz + 0.5) * spacing
                ])
                
                # 添加小的随机扰动（±5%）
                perturbation = (np.random.random(3) - 0.5) * spacing * 0.1
                o_pos = base_pos + perturbation
                
                # 确保在盒子内
                o_pos = o_pos % box_length
                
                # 随机旋转
                euler_angles = np.random.random(3) * 2 * np.pi
                R = euler_to_rotation_matrix(*euler_angles)
                
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
                
                # M位点（沿角平分线反方向）
                bisector = -(h1_std + h2_std)
                bisector_norm = bisector / np.linalg.norm(bisector)
                m_std = bisector_norm * r_om
                
                # 应用旋转
                h1_rot = R @ h1_std
                h2_rot = R @ h2_std
                m_rot = R @ m_std
                
                # 存储5个原子的位置
                idx = water_count * 5
                positions[idx] = o_pos           # O
                positions[idx+1] = o_pos         # D (初始与O重合)
                positions[idx+2] = o_pos + h1_rot  # H1
                positions[idx+3] = o_pos + h2_rot  # H2
                positions[idx+4] = o_pos + m_rot   # M
                
                water_count += 1
                
                # 进度报告
                if water_count % 1000 == 0:
                    print(f"    已生成 {water_count}/{n_waters} 水分子...")
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    elapsed = time.time() - start_time
    print(f"  生成完成，用时 {elapsed:.1f} 秒")
    
    return positions, box_length, actual_density

def save_system(n_waters, positions, box_length, density, output_dir):
    """
    保存系统到文件
    """
    print(f"\n保存 {n_waters} 水系统...")
    
    data = {
        'n_waters': n_waters,
        'positions': positions,
        'box_length': box_length,
        'density': density,
        'charges': [1.71636, -1.71636, 0.55733, 0.55733, -1.11466],  # O, D, H1, H2, M
        'atom_names': ['O', 'D', 'H1', 'H2', 'M']
    }
    
    # 保存pickle文件
    pickle_file = os.path.join(output_dir, f'water_{n_waters}.pkl')
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    print(f"  保存到: {pickle_file}")
    
    # 对于大系统，PDB文件可能太大，只为较小的系统生成
    if n_waters <= 512:
        pdb_file = os.path.join(output_dir, f'water_{n_waters}.pdb')
        print(f"  生成PDB文件: {pdb_file}")
        
        with open(pdb_file, 'w') as f:
            # PDB头部
            f.write(f"CRYST1{box_length*10:9.3f}{box_length*10:9.3f}{box_length*10:9.3f}")
            f.write(f"  90.00  90.00  90.00 P 1           1\n")
            f.write(f"REMARK   Density = {density:.6f} g/cm³\n")
            f.write(f"REMARK   SWM4-NDP water model with 5 sites per molecule\n")
            f.write(f"REMARK   System size: {n_waters} water molecules\n")
            
            atom_idx = 1
            for i in range(n_waters):
                base_idx = i * 5
                res_num = (i % 9999) + 1  # PDB限制残基号为4位
                
                # O原子
                f.write(f"ATOM  {atom_idx:5d}  O   WAT A{res_num:4d}    ")
                f.write(f"{positions[base_idx][0]*10:8.3f}")
                f.write(f"{positions[base_idx][1]*10:8.3f}")
                f.write(f"{positions[base_idx][2]*10:8.3f}")
                f.write(f"  1.00  0.00           O\n")
                atom_idx += 1
                
                # H1原子
                f.write(f"ATOM  {atom_idx:5d}  H1  WAT A{res_num:4d}    ")
                f.write(f"{positions[base_idx+2][0]*10:8.3f}")
                f.write(f"{positions[base_idx+2][1]*10:8.3f}")
                f.write(f"{positions[base_idx+2][2]*10:8.3f}")
                f.write(f"  1.00  0.00           H\n")
                atom_idx += 1
                
                # H2原子
                f.write(f"ATOM  {atom_idx:5d}  H2  WAT A{res_num:4d}    ")
                f.write(f"{positions[base_idx+3][0]*10:8.3f}")
                f.write(f"{positions[base_idx+3][1]*10:8.3f}")
                f.write(f"{positions[base_idx+3][2]*10:8.3f}")
                f.write(f"  1.00  0.00           H\n")
                atom_idx += 1
                
                # 每1000个水分子报告进度
                if (i + 1) % 1000 == 0:
                    print(f"    写入PDB: {i+1}/{n_waters}")
            
            f.write("END\n")
    else:
        print(f"  跳过PDB文件（系统太大）")

def main():
    """
    主函数
    """
    print("生成大型水系统用于Drude模型测试")
    print("="*70)
    
    # 创建输出目录
    output_dir = '../tests/performance/large_water_systems'
    os.makedirs(output_dir, exist_ok=True)
    print(f"输出目录: {output_dir}")
    
    # 要生成的系统大小
    system_sizes = [256, 512, 1024, 2048, 4096]
    
    # 汇总信息
    print(f"\n计划生成的系统:")
    print(f"{'水分子数':>8} {'盒子(nm)':>10} {'原子总数':>10} {'文件大小(MB)':>14}")
    print("-"*50)
    
    for n_waters in system_sizes:
        box = calculate_box_length_for_density(n_waters, 1.0)
        n_atoms = n_waters * 5
        # 估算pickle文件大小（每个位置3个float64 = 24字节）
        size_mb = n_atoms * 3 * 8 / 1024 / 1024
        print(f"{n_waters:8d} {box:10.3f} {n_atoms:10d} {size_mb:14.1f}")
    
    # 生成系统
    total_start = time.time()
    
    for n_waters in system_sizes:
        print(f"\n{'='*70}")
        positions, box_length, density = create_large_water_system(n_waters, 1.0)
        save_system(n_waters, positions, box_length, density, output_dir)
    
    total_elapsed = time.time() - total_start
    print(f"\n\n所有系统生成完成！")
    print(f"总用时: {total_elapsed:.1f} 秒 ({total_elapsed/60:.1f} 分钟)")
    
    # 统计信息
    print(f"\n生成的文件:")
    for n_waters in system_sizes:
        pickle_file = os.path.join(output_dir, f'water_{n_waters}.pkl')
        if os.path.exists(pickle_file):
            size_mb = os.path.getsize(pickle_file) / 1024 / 1024
            print(f"  water_{n_waters}.pkl: {size_mb:.1f} MB")

if __name__ == "__main__":
    main()