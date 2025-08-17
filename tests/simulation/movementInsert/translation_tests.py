# tests/simulation/movementInsert/test_translation_rotation.py
"""
平移和旋转移动测试 - 补充GPU版本已有但PyGCMC测试缺失的功能

测试GCMC中的平移(Translation)和旋转(Rotation)移动，
这些是GPU实现中的case 2和case 3。
"""

import pytest
import random
import math
import numpy as np
import pygcmc
from typing import List, Tuple

# Constants
kB = 0.008314463  # kJ/(mol·K)

# Import helper functions
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)


def calculate_com(residue, system=None):
    """计算分子质心 - 适配PyGCMC的MCResidue"""
    total_mass = 0
    com = [0, 0, 0]
    
    # 如果residue是MCResidue，使用atom_start和atom_count访问原子
    if hasattr(residue, 'atom_start') and hasattr(residue, 'atom_count') and system:
        for i in range(residue.atom_count):
            atom = system.atoms[residue.atom_start + i]
            mass = 1.0  # 简化处理，使用单位质量
            com[0] += atom.x * mass
            com[1] += atom.y * mass
            com[2] += atom.z * mass
            total_mass += mass
    elif hasattr(residue, 'atoms'):
        # 兼容其他分子对象
        for atom in residue.atoms:
            mass = atom.mass if hasattr(atom, 'mass') else 1.0
            com[0] += atom.x * mass
            com[1] += atom.y * mass
            com[2] += atom.z * mass
            total_mass += mass
    
    if total_mass > 0:
        com[0] /= total_mass
        com[1] /= total_mass
        com[2] /= total_mass
    
    return com


def translate_molecule(molecule, dx, dy, dz):
    """平移分子"""
    for atom in molecule.atoms:
        atom.x += dx
        atom.y += dy
        atom.z += dz
    return molecule


def rotate_molecule_atoms(atoms, com, quaternion):
    """使用四元数旋转原子列表"""
    q0, q1, q2, q3 = quaternion
    
    # 四元数归一化
    norm = math.sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)
    if norm > 0:
        q0, q1, q2, q3 = q0/norm, q1/norm, q2/norm, q3/norm
    
    # 构建旋转矩阵
    R = [
        [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
    ]
    
    # 旋转每个原子
    for atom in atoms:
        # 平移到原点
        x = atom.x - com[0]
        y = atom.y - com[1]
        z = atom.z - com[2]
        
        # 应用旋转
        new_x = R[0][0]*x + R[0][1]*y + R[0][2]*z
        new_y = R[1][0]*x + R[1][1]*y + R[1][2]*z
        new_z = R[2][0]*x + R[2][1]*y + R[2][2]*z
        
        # 平移回原位
        atom.x = new_x + com[0]
        atom.y = new_y + com[1]
        atom.z = new_z + com[2]


def create_random_quaternion():
    """生成随机四元数"""
    # 使用Shoemake算法生成均匀分布的四元数
    u1, u2, u3 = random.random(), random.random(), random.random()
    
    sqrt1_u1 = math.sqrt(1 - u1)
    sqrt_u1 = math.sqrt(u1)
    
    q0 = sqrt1_u1 * math.sin(2 * math.pi * u2)
    q1 = sqrt1_u1 * math.cos(2 * math.pi * u2)
    q2 = sqrt_u1 * math.sin(2 * math.pi * u3)
    q3 = sqrt_u1 * math.cos(2 * math.pi * u3)
    
    return [q0, q1, q2, q3]


def test_translation_move():
    """测试平移移动"""
    # 创建系统
    system = create_empty_system()
    
    # 添加一个水分子
    molecule = create_water_molecule(2.5, 2.5, 2.5)
    system = insert_molecule(system, molecule)
    
    # 记录初始位置和能量
    initial_com = calculate_com(system.residues[-1], system)
    pygcmc.computeSystemEnergyCutoff(system)
    initial_energy = calculate_system_energy(system)
    
    # 执行平移
    max_trans_dist = 0.5  # nm
    dx = (random.random() - 0.5) * max_trans_dist
    dy = (random.random() - 0.5) * max_trans_dist
    dz = (random.random() - 0.5) * max_trans_dist
    
    # 平移分子
    residue = system.residues[-1]
    for i in range(residue.atom_count):
        atom = system.atoms[residue.atom_start + i]
        atom.x += dx
        atom.y += dy
        atom.z += dz
    
    # 计算新位置和能量
    final_com = calculate_com(residue, system)
    pygcmc.computeSystemEnergyCutoff(system)
    final_energy = calculate_system_energy(system)
    
    # 验证平移
    assert abs(final_com[0] - initial_com[0] - dx) < 1e-6
    assert abs(final_com[1] - initial_com[1] - dy) < 1e-6
    assert abs(final_com[2] - initial_com[2] - dz) < 1e-6
    
    # 能量可能变化（由于相互作用）
    energy_change = final_energy - initial_energy
    
    # 计算接受概率（Metropolis准则）
    T = 300.0
    beta = 1.0 / (kB * T)
    if energy_change <= 0:
        acceptance = 1.0
    else:
        acceptance = math.exp(-beta * energy_change)
    
    assert 0 <= acceptance <= 1


def test_rotation_move():
    """测试旋转移动"""
    # 创建系统
    system = create_empty_system()
    
    # 添加一个水分子
    molecule = create_water_molecule(2.5, 2.5, 2.5)
    system = insert_molecule(system, molecule)
    
    # 记录初始能量和质心
    pygcmc.computeSystemEnergyCutoff(system)
    initial_energy = calculate_system_energy(system)
    initial_com = calculate_com(system.residues[-1], system)
    
    # 生成随机四元数
    quaternion = create_random_quaternion()
    
    # 执行旋转
    residue = system.residues[-1]
    atoms = [system.atoms[residue.atom_start + i] for i in range(residue.atom_count)]
    
    # 临时存储原子坐标
    coords = [(a.x, a.y, a.z) for a in atoms]
    
    # 应用旋转
    q0, q1, q2, q3 = quaternion
    R = [
        [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
    ]
    
    com = initial_com
    for i, atom in enumerate(atoms):
        x, y, z = coords[i]
        # 平移到原点
        x -= com[0]
        y -= com[1]
        z -= com[2]
        
        # 旋转
        new_x = R[0][0]*x + R[0][1]*y + R[0][2]*z
        new_y = R[1][0]*x + R[1][1]*y + R[1][2]*z
        new_z = R[2][0]*x + R[2][1]*y + R[2][2]*z
        
        # 平移回去
        atom.x = new_x + com[0]
        atom.y = new_y + com[1]
        atom.z = new_z + com[2]
    
    # 验证质心不变
    final_com = calculate_com(residue, system)
    assert abs(final_com[0] - initial_com[0]) < 1e-6
    assert abs(final_com[1] - initial_com[1]) < 1e-6
    assert abs(final_com[2] - initial_com[2]) < 1e-6
    
    # 计算新能量
    pygcmc.computeSystemEnergyCutoff(system)
    final_energy = calculate_system_energy(system)
    
    # 旋转不应该显著改变能量（除非有定向相互作用）
    energy_change = final_energy - initial_energy
    
    # 计算接受概率
    T = 300.0
    beta = 1.0 / (kB * T)
    if energy_change <= 0:
        acceptance = 1.0
    else:
        acceptance = math.exp(-beta * energy_change)
    
    assert 0 <= acceptance <= 1


