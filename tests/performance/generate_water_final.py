#!/usr/bin/env python3
"""
最终版本：生成合理密度的水分子系统用于CG vs SCF测试
"""

import numpy as np
import pickle
import os

def create_water_system_with_density(n_waters, target_density=1.0):
    """
    创建指定密度的水分子系统
    target_density: g/cm³
    """
    # 计算盒子大小
    mass_per_water = 18.015 / 6.022e23  # g
    total_mass = n_waters * mass_per_water  # g
    volume_cm3 = total_mass / target_density  # cm³
    volume_nm3 = volume_cm3 * 1e21  # nm³
    box_length = (volume_nm3 ** (1.0/3.0))  # nm
    
    # 在立方格子上放置水分子，加入随机扰动
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    water_count = 0
    
    np.random.seed(42)  # 可重复性
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础位置 + 随机扰动
                x = (i + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.2
                y = (j + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.2
                z = (k + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.2
                
                # 确保在盒子内
                x = x % box_length
                y = y % box_length
                z = z % box_length
                
                # 随机旋转水分子
                theta = np.random.random() * 2 * np.pi
                phi = np.random.random() * np.pi
                psi = np.random.random() * 2 * np.pi
                
                # 创建旋转矩阵
                R = create_rotation_matrix(theta, phi, psi)
                
                # SWM4-NDP几何参数
                r_oh = 0.09572  # nm
                angle_hoh = 104.52 * np.pi / 180.0
                
                # 原始H原子位置（相对于O）
                h1_local = np.array([
                    r_oh * np.sin(angle_hoh/2),
                    0,
                    r_oh * np.cos(angle_hoh/2)
                ])
                h2_local = np.array([
                    -r_oh * np.sin(angle_hoh/2),
                    0,
                    r_oh * np.cos(angle_hoh/2)
                ])
                
                # 旋转
                h1_rotated = R @ h1_local
                h2_rotated = R @ h2_local
                
                # SWM4-NDP的5个位点
                # 1. O原子
                positions.append([x, y, z])
                
                # 2. Drude粒子（初始与O重合）
                positions.append([x, y, z])
                
                # 3. H1原子
                positions.append([
                    x + h1_rotated[0],
                    y + h1_rotated[1],
                    z + h1_rotated[2]
                ])
                
                # 4. H2原子
                positions.append([
                    x + h2_rotated[0],
                    y + h2_rotated[1],
                    z + h2_rotated[2]
                ])
                
                # 5. M位点（虚拟位点）
                # M位于O原子沿着HOH角平分线反方向0.024034 nm处
                bisector = -(h1_rotated + h2_rotated)
                bisector = bisector / np.linalg.norm(bisector) * 0.024034
                positions.append([
                    x + bisector[0],
                    y + bisector[1],
                    z + bisector[2]
                ])
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return positions, box_length

def create_rotation_matrix(theta, phi, psi):
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
        [cos_psi*cos_theta - sin_psi*cos_phi*sin_theta, 
         -sin_psi*cos_theta - cos_psi*cos_phi*sin_theta,
         sin_phi*sin_theta],
        [cos_psi*sin_theta + sin_psi*cos_phi*cos_theta,
         -sin_psi*sin_theta + cos_psi*cos_phi*cos_theta,
         -sin_phi*cos_theta],
        [sin_psi*sin_phi, 
         cos_psi*sin_phi, 
         cos_phi]
    ])
    
    return R

def check_minimum_distances(positions, n_waters, box_length):
    """
    检查最小原子间距离
    """
    min_dist = float('inf')
    
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            # O-O距离
            o1_idx = i * 5
            o2_idx = j * 5
            
            dr = np.array(positions[o2_idx]) - np.array(positions[o1_idx])
            # PBC
            dr = dr - box_length * np.round(dr / box_length)
            dist = np.linalg.norm(dr)
            
            if dist < min_dist:
                min_dist = dist
    
    return min_dist

def save_system(n_waters, positions, box_length, output_dir):
    """
    保存系统
    """
    # 验证密度
    volume_nm3 = box_length ** 3
    mass_g = n_waters * 18.015 / 6.022e23
    volume_cm3 = volume_nm3 * 1e-21
    density = mass_g / volume_cm3
    
    # 检查最小距离
    min_dist = check_minimum_distances(positions, n_waters, box_length)
    
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
    
    print(f"{n_waters:4d} 水: 盒子 {box_length:6.3f} nm, "
          f"密度 {density:5.3f} g/cm³, 最小O-O距离 {min_dist:5.3f} nm")

def main():
    """
    主函数
    """
    print("生成合理密度的SWM4-NDP水分子系统")
    print("="*70)
    
    # 创建输出目录
    output_dir = '../tests/performance/water_systems_final'
    os.makedirs(output_dir, exist_ok=True)
    
    # 系统大小
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256]
    
    print(f"\n{'水数':>6} {'盒子(nm)':>10} {'密度(g/cm³)':>12} {'最小O-O距离(nm)':>15}")
    print("-"*50)
    
    for n_waters in system_sizes:
        # 根据系统大小调整密度，避免原子重叠
        if n_waters <= 8:
            target_density = 0.8  # 小系统用较低密度
        else:
            target_density = 1.0  # 大系统用标准密度
        
        positions, box_length = create_water_system_with_density(n_waters, target_density)
        save_system(n_waters, positions, box_length, output_dir)
    
    print("\n完成! 生成的系统具有合理的密度和原子间距。")
    print("\n使用方法:")
    print("1. 更新test_cg_final_comparison.py中的目录为'water_systems_final'")
    print("2. 运行: cd build && PYTHONPATH=$PYTHONPATH:./modules/bindings python ../tests/performance/test_cg_final_comparison.py")

if __name__ == "__main__":
    main()