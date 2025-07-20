#!/usr/bin/env python3
"""
生成优化的水体系PDB文件
"""

import numpy as np
import os

def write_pdb(filename, n_waters, positions, box_vectors):
    """
    写入PDB文件
    """
    with open(filename, 'w') as f:
        # 写入标题
        f.write(f"TITLE     SWM4-NDP Water System with {n_waters} molecules\n")
        f.write(f"REMARK    Generated for Drude model testing\n")
        
        # 写入盒子信息
        if box_vectors is not None:
            a = box_vectors[0] * 10  # nm to Å
            b = box_vectors[1] * 10
            c = box_vectors[2] * 10
            f.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}  90.00  90.00  90.00 P 1           1\n")
        
        # 写入原子
        atom_idx = 1
        for i in range(n_waters):
            res_idx = i + 1
            base = i * 5
            
            # O atom
            f.write(f"ATOM  {atom_idx:5d}  O   WAT {res_idx:5d}    "
                   f"{positions[base][0]:8.3f}{positions[base][1]:8.3f}{positions[base][2]:8.3f}"
                   f"  1.00  0.00           O\n")
            atom_idx += 1
            
            # D atom (Drude)
            f.write(f"ATOM  {atom_idx:5d}  D   WAT {res_idx:5d}    "
                   f"{positions[base+1][0]:8.3f}{positions[base+1][1]:8.3f}{positions[base+1][2]:8.3f}"
                   f"  1.00  0.00           D\n")
            atom_idx += 1
            
            # H1 atom
            f.write(f"ATOM  {atom_idx:5d}  H1  WAT {res_idx:5d}    "
                   f"{positions[base+2][0]:8.3f}{positions[base+2][1]:8.3f}{positions[base+2][2]:8.3f}"
                   f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # H2 atom
            f.write(f"ATOM  {atom_idx:5d}  H2  WAT {res_idx:5d}    "
                   f"{positions[base+3][0]:8.3f}{positions[base+3][1]:8.3f}{positions[base+3][2]:8.3f}"
                   f"  1.00  0.00           H\n")
            atom_idx += 1
            
            # M atom (virtual site)
            f.write(f"ATOM  {atom_idx:5d}  M   WAT {res_idx:5d}    "
                   f"{positions[base+4][0]:8.3f}{positions[base+4][1]:8.3f}{positions[base+4][2]:8.3f}"
                   f"  1.00  0.00           M\n")
            atom_idx += 1
        
        f.write("END\n")

def generate_optimized_water_pdb(n_waters, output_dir="water_pdbs"):
    """
    生成优化的水体系PDB文件
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 基于密度计算盒子大小
    density = 1.0  # g/cm³
    mass_per_water = 18.015  # g/mol
    avogadro = 6.022e23
    volume_per_water = mass_per_water / (density * avogadro) * 1e24  # nm³
    total_volume = n_waters * volume_per_water
    box_length = (total_volume) ** (1.0/3.0)
    
    # 生成初始位置
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    
    # 水分子的标准几何
    oh_bond = 0.9572  # Å
    angle = 104.52 * np.pi / 180
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 氧原子位置（Å）
                x = (i + 0.5) * spacing * 10  # nm to Å
                y = (j + 0.5) * spacing * 10
                z = (k + 0.5) * spacing * 10
                
                # 添加随机扰动
                x += (np.random.rand() - 0.5) * 1.0
                y += (np.random.rand() - 0.5) * 1.0
                z += (np.random.rand() - 0.5) * 1.0
                
                # 随机旋转
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                
                # 计算H原子位置
                h1_x = oh_bond * np.sin(angle/2) * np.cos(theta)
                h1_y = oh_bond * np.sin(angle/2) * np.sin(theta)
                h1_z = oh_bond * np.cos(angle/2)
                
                h2_x = oh_bond * np.sin(angle/2) * np.cos(theta + np.pi)
                h2_y = oh_bond * np.sin(angle/2) * np.sin(theta + np.pi)
                h2_z = oh_bond * np.cos(angle/2)
                
                # M site (along bisector, 0.24034 Å from O)
                m_vec = np.array([h1_x + h2_x, h1_y + h2_y, h1_z + h2_z])
                m_vec = m_vec / np.linalg.norm(m_vec) * 0.24034
                
                # 添加5个原子的位置
                positions.append([x, y, z])  # O
                positions.append([x, y, z])  # D (初始与O重合)
                positions.append([x + h1_x, y + h1_y, z + h1_z])  # H1
                positions.append([x + h2_x, y + h2_y, z + h2_z])  # H2
                positions.append([x + m_vec[0], y + m_vec[1], z + m_vec[2]])  # M
                
                water_count += 1
    
    # 写入PDB文件
    filename = os.path.join(output_dir, f"water_{n_waters}.pdb")
    write_pdb(filename, n_waters, positions, [box_length, box_length, box_length])
    
    print(f"生成了 {filename}")
    print(f"  盒子大小: {box_length:.3f} nm")
    print(f"  密度: {density:.3f} g/cm³")
    
    return filename

def generate_test_systems():
    """
    生成一系列测试系统
    """
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256]
    
    print("生成优化的水体系PDB文件")
    print("="*60)
    
    for n_waters in system_sizes:
        generate_optimized_water_pdb(n_waters)
    
    print("\n完成！PDB文件保存在 water_pdbs/ 目录下")
    print("可以使用VMD或PyMOL查看这些文件")

if __name__ == "__main__":
    generate_test_systems()