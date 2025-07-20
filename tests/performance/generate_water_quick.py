#!/usr/bin/env python3
"""
快速生成SWM4-NDP水分子测试系统 (不做优化)
"""

import numpy as np
import pickle
import os

def create_water_system(n_waters):
    """
    创建水分子系统，密度接近1.0 g/cm³
    """
    # 目标密度 1.0 g/cm³
    # 每个水分子质量 18.015 g/mol / 6.022e23 = 2.992e-23 g
    # 每个水分子体积 2.992e-23 g / 1.0 g/cm³ = 2.992e-23 cm³ = 29.92 Å³ = 0.02992 nm³
    volume_per_water = 0.02992  # nm³
    total_volume = n_waters * volume_per_water
    box_length = (total_volume ** (1.0/3.0))
    
    # 在立方格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    water_count = 0
    
    # 添加一些随机性避免完全规则排列
    np.random.seed(42)
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 水分子中心 (加入小的随机偏移)
                x = (i + 0.5) * spacing + (np.random.random() - 0.5) * 0.02
                y = (j + 0.5) * spacing + (np.random.random() - 0.5) * 0.02
                z = (k + 0.5) * spacing + (np.random.random() - 0.5) * 0.02
                
                # 随机旋转
                theta = np.random.random() * 2 * np.pi
                phi = np.random.random() * np.pi
                psi = np.random.random() * 2 * np.pi
                
                # 旋转矩阵
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
                    [sin_psi*sin_phi, cos_psi*sin_phi, cos_phi]
                ])
                
                # SWM4-NDP几何参数
                r_oh = 0.09572  # nm
                angle_hoh = 104.52 * np.pi / 180.0
                
                # 原始坐标 (分子坐标系)
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
                
                # 旋转到实验室坐标系
                h1_lab = R @ h1_local
                h2_lab = R @ h2_local
                
                # O原子位置
                positions.append([x, y, z])
                
                # Drude粒子位置 (初始时与O重合)
                positions.append([x, y, z])
                
                # H1原子位置
                positions.append([x + h1_lab[0], y + h1_lab[1], z + h1_lab[2]])
                
                # H2原子位置
                positions.append([x + h2_lab[0], y + h2_lab[1], z + h2_lab[2]])
                
                # M位点位置 (在O-H键的延长线上)
                # SWM4-NDP: M位点距O原子0.024034 nm
                m_vec = -(h1_lab + h2_lab) / np.linalg.norm(h1_lab + h2_lab) * 0.024034
                positions.append([x + m_vec[0], y + m_vec[1], z + m_vec[2]])
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return positions, box_length

def save_system(n_waters, positions, box_length, output_dir):
    """
    保存系统数据
    """
    # 计算实际密度
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
    
    print(f"{n_waters:4d} 水: 盒子 {box_length:6.3f} nm, 密度 {density:5.3f} g/cm³")
    
    return data

def main():
    """
    主函数
    """
    print("快速生成SWM4-NDP水分子系统")
    print("="*50)
    
    # 创建输出目录
    output_dir = '../tests/performance/water_systems'
    os.makedirs(output_dir, exist_ok=True)
    
    # 系统大小
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256]
    
    print(f"\n{'水数':>6} {'盒子(nm)':>10} {'密度(g/cm³)':>12}")
    print("-"*30)
    
    for n_waters in system_sizes:
        positions, box_length = create_water_system(n_waters)
        save_system(n_waters, positions, box_length, output_dir)
    
    print("\n完成! 现在运行测试:")
    print("cd build && python ../tests/performance/test_cg_final_comparison.py")

if __name__ == "__main__":
    main()