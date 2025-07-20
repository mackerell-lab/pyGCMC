#!/usr/bin/env python3
"""
简化版：生成预优化的水体系用于测试
"""

import numpy as np
import os
import pickle

def generate_water_system(n_waters, optimize=True):
    """
    生成水体系，可选择是否优化
    """
    # 基于水的密度计算盒子
    density = 1.0  # g/cm³
    mass_per_water = 18.015  # g/mol
    avogadro = 6.022e23
    volume_per_water = mass_per_water / (density * avogadro) * 1e24  # nm³
    
    # 稍微大一点的初始体积（密度约0.95）
    total_volume = n_waters * volume_per_water * 1.05
    box_length = (total_volume) ** (1.0/3.0)
    
    # 确保盒子足够大以满足截断距离要求
    min_box = 2.5  # nm，最小盒子大小
    if box_length < min_box:
        box_length = min_box
    
    print(f"\n生成 {n_waters} 水分子系统:")
    print(f"  盒子大小: {box_length:.3f} nm")
    print(f"  密度: {n_waters * mass_per_water / (box_length**3 * avogadro * 1e-24):.3f} g/cm³")
    
    # 创建格子排列
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    positions = []
    
    # 水分子的标准几何（SWM4-NDP）
    oh_bond = 0.09572  # nm
    hoh_angle = 104.52 * np.pi / 180  # rad
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # O原子位置
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                if optimize:
                    # 添加随机扰动避免完美晶格
                    x += (np.random.rand() - 0.5) * spacing * 0.3
                    y += (np.random.rand() - 0.5) * spacing * 0.3
                    z += (np.random.rand() - 0.5) * spacing * 0.3
                    
                    # 随机旋转水分子
                    theta = np.random.rand() * 2 * np.pi
                    phi = np.random.rand() * np.pi
                    psi = np.random.rand() * 2 * np.pi
                else:
                    # 固定取向
                    theta = 0
                    phi = 0
                    psi = 0
                
                # 计算H原子相对位置
                # H1在xz平面
                h1_x = oh_bond * np.sin(hoh_angle/2) * np.cos(theta)
                h1_y = oh_bond * np.sin(hoh_angle/2) * np.sin(theta)
                h1_z = oh_bond * np.cos(hoh_angle/2)
                
                # H2在xz平面，对称
                h2_x = -oh_bond * np.sin(hoh_angle/2) * np.cos(theta)
                h2_y = -oh_bond * np.sin(hoh_angle/2) * np.sin(theta)
                h2_z = oh_bond * np.cos(hoh_angle/2)
                
                # 应用额外旋转
                if optimize and (phi != 0 or psi != 0):
                    # 简化：只绕z轴旋转
                    cos_psi = np.cos(psi)
                    sin_psi = np.sin(psi)
                    
                    h1_x_new = h1_x * cos_psi - h1_y * sin_psi
                    h1_y_new = h1_x * sin_psi + h1_y * cos_psi
                    h1_x, h1_y = h1_x_new, h1_y_new
                    
                    h2_x_new = h2_x * cos_psi - h2_y * sin_psi
                    h2_y_new = h2_x * sin_psi + h2_y * cos_psi
                    h2_x, h2_y = h2_x_new, h2_y_new
                
                # M位点：沿O-H平分线，距离0.024034 nm
                m_vec_x = (h1_x + h2_x) / 2
                m_vec_y = (h1_y + h2_y) / 2
                m_vec_z = (h1_z + h2_z) / 2
                m_norm = np.sqrt(m_vec_x**2 + m_vec_y**2 + m_vec_z**2)
                if m_norm > 0:
                    m_vec_x = m_vec_x / m_norm * 0.024034
                    m_vec_y = m_vec_y / m_norm * 0.024034
                    m_vec_z = m_vec_z / m_norm * 0.024034
                
                # 添加5个原子的位置
                positions.extend([
                    [x, y, z],  # O
                    [x, y, z],  # D (初始与O重合)
                    [x + h1_x, y + h1_y, z + h1_z],  # H1
                    [x + h2_x, y + h2_y, z + h2_z],  # H2
                    [x + m_vec_x, y + m_vec_y, z + m_vec_z]  # M
                ])
                
                water_count += 1
    
    # 转换为numpy数组
    positions = np.array(positions)
    
    # 应用周期性边界条件
    positions = positions % box_length
    
    return positions, box_length

def save_water_system(n_waters, positions, box_length, output_dir="optimized_water_systems"):
    """
    保存水体系
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存为简单格式
    data = {
        'n_waters': n_waters,
        'positions': positions,  # nm单位
        'box_length': box_length,  # nm
        'charges': [1.71636, -1.71636, 0.55733, 0.55733, -1.11466],  # 每个水分子的电荷
        'atom_types': [0, 1, 2, 2, 3],  # 原子类型
        'description': 'SWM4-NDP water system'
    }
    
    # 保存pickle
    pickle_file = os.path.join(output_dir, f'water_{n_waters}.pkl')
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    print(f"  保存到: {pickle_file}")
    
    # 保存PDB
    pdb_file = os.path.join(output_dir, f'water_{n_waters}.pdb')
    write_pdb(pdb_file, n_waters, positions, box_length)
    print(f"  PDB文件: {pdb_file}")
    
    return data

def write_pdb(filename, n_waters, positions, box_length):
    """
    写入PDB文件
    """
    with open(filename, 'w') as f:
        f.write(f"TITLE     SWM4-NDP Water System with {n_waters} molecules\n")
        f.write(f"REMARK    Optimized structure for Drude testing\n")
        
        # 盒子信息
        box_angstrom = box_length * 10  # nm to Å
        f.write(f"CRYST1{box_angstrom:9.3f}{box_angstrom:9.3f}{box_angstrom:9.3f}"
                f"  90.00  90.00  90.00 P 1           1\n")
        
        # 原子
        atom_idx = 1
        atom_names = ['O', 'D', 'H1', 'H2', 'M']
        elements = ['O', 'D', 'H', 'H', 'M']
        
        for i in range(n_waters):
            res_idx = i + 1
            for j in range(5):
                idx = i * 5 + j
                x, y, z = positions[idx] * 10  # nm to Å
                
                f.write(f"ATOM  {atom_idx:5d}  {atom_names[j]:<3s} WAT {res_idx:5d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}"
                       f"  1.00  0.00           {elements[j]}\n")
                atom_idx += 1
        
        f.write("END\n")

def convert_to_pygcmc(data):
    """
    将数据转换为pygcmc格式
    """
    import pygcmc
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    charges = data['charges']
    atom_types = data['atom_types']
    
    # 创建状态
    atoms = []
    residues = []
    
    for i in range(n_waters):
        for j in range(5):
            atom = pygcmc.MCAtom()
            idx = i * 5 + j
            atom.x = positions[idx][0]
            atom.y = positions[idx][1]
            atom.z = positions[idx][2]
            atom.charge = charges[j]
            atom.type = atom_types[j]
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)
    
    # SWM4-NDP力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def main():
    """
    生成一系列优化的水体系
    """
    print("生成优化的SWM4-NDP水体系")
    print("="*60)
    
    system_sizes = [2, 4, 8, 16, 32, 64, 128, 256]
    
    for n_waters in system_sizes:
        positions, box_length = generate_water_system(n_waters, optimize=True)
        save_water_system(n_waters, positions, box_length)
    
    print("\n完成！系统保存在 optimized_water_systems/ 目录")
    
    # 测试转换
    print("\n测试pygcmc转换...")
    test_file = "optimized_water_systems/water_4.pkl"
    if os.path.exists(test_file):
        with open(test_file, 'rb') as f:
            data = pickle.load(f)
        try:
            state = convert_to_pygcmc(data)
            print(f"  成功转换 {data['n_waters']} 水系统")
            print(f"  原子数: {state.activeAtomCount}")
            print(f"  盒子: {state.info.box}")
        except Exception as e:
            print(f"  转换失败: {e}")

if __name__ == "__main__":
    main()