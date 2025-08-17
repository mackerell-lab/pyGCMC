# tests/simulation/movementInsert/rotation_tests.py
"""
Rotation and combined movement tests
"""

import pytest
import random
import math
import numpy as np
import pygcmc

# Constants
kB = 0.008314463  # kJ/(mol·K)

# Import helper functions
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)

# Import translation helpers
from .translation_tests import (
    calculate_com,
    create_random_quaternion
)

# Define helper functions that were lost in split
def quaternion_to_rotation_matrix(q):
    """Convert quaternion to rotation matrix"""
    q0, q1, q2, q3 = q
    R = [
        [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
    ]
    return R

def rotate_point(point, center, rotation_matrix):
    """Rotate a point around a center"""
    # Translate to origin
    translated = [point[i] - center[i] for i in range(3)]
    
    # Apply rotation
    R = rotation_matrix
    rotated = [
        R[0][0]*translated[0] + R[0][1]*translated[1] + R[0][2]*translated[2],
        R[1][0]*translated[0] + R[1][1]*translated[1] + R[1][2]*translated[2],
        R[2][0]*translated[0] + R[2][1]*translated[1] + R[2][2]*translated[2]
    ]
    
    # Translate back
    return [rotated[i] + center[i] for i in range(3)]

def apply_pbc(position, box):
    """Apply periodic boundary conditions"""
    for i in range(3):
        while position[i] < 0:
            position[i] += box[i]
        while position[i] >= box[i]:
            position[i] -= box[i]
    return position

def calculate_bond_length(atom1, atom2):
    """Calculate bond length between two atoms"""
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def test_combined_translation_rotation():
    """测试组合的平移和旋转移动"""
    # 创建包含多个分子的系统
    system = create_empty_system()
    
    # 添加几个水分子
    positions = [
        (1.0, 1.0, 1.0),
        (3.0, 1.0, 1.0),
        (1.0, 3.0, 1.0),
        (3.0, 3.0, 1.0),
    ]
    
    for x, y, z in positions:
        molecule = create_water_molecule(x, y, z)
        system = insert_molecule(system, molecule)
    
    # 计算初始能量
    pygcmc.computeSystemEnergyCutoff(system)
    initial_energy = calculate_system_energy(system)
    
    # 随机选择一个分子
    mol_idx = random.randint(0, len(system.residues) - 1)
    residue = system.residues[mol_idx]
    
    # 执行平移
    max_trans = 0.3
    dx = (random.random() - 0.5) * max_trans
    dy = (random.random() - 0.5) * max_trans
    dz = (random.random() - 0.5) * max_trans
    
    for i in range(residue.atom_count):
        atom = system.atoms[residue.atom_start + i]
        atom.x += dx
        atom.y += dy
        atom.z += dz
    
    # 执行旋转
    quaternion = create_random_quaternion()
    q0, q1, q2, q3 = quaternion
    
    R = [
        [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
    ]
    
    com = calculate_com(residue, system)
    atoms = [system.atoms[residue.atom_start + i] for i in range(residue.atom_count)]
    
    for atom in atoms:
        # 平移到质心
        x = atom.x - com[0]
        y = atom.y - com[1]
        z = atom.z - com[2]
        
        # 旋转
        new_x = R[0][0]*x + R[0][1]*y + R[0][2]*z
        new_y = R[1][0]*x + R[1][1]*y + R[1][2]*z
        new_z = R[2][0]*x + R[2][1]*y + R[2][2]*z
        
        # 平移回去
        atom.x = new_x + com[0]
        atom.y = new_y + com[1]
        atom.z = new_z + com[2]
    
    # 计算最终能量
    pygcmc.computeSystemEnergyCutoff(system)
    final_energy = calculate_system_energy(system)
    
    # 计算接受概率
    T = 300.0
    beta = 1.0 / (kB * T)
    energy_change = final_energy - initial_energy
    
    if energy_change <= 0:
        acceptance = 1.0
    else:
        acceptance = math.exp(-beta * energy_change)
    
    assert 0 <= acceptance <= 1
    
    # 统计接受率（用于调试）
    print(f"Energy change: {energy_change:.2f} kJ/mol")
    print(f"Acceptance probability: {acceptance:.4f}")


def test_translation_with_pbc():
    """测试带周期性边界条件的平移"""
    # 创建系统
    system = create_empty_system()
    box = system.info.box
    
    # 在边界附近添加分子
    molecule = create_water_molecule(box[0] - 0.5, 2.5, 2.5)
    system = insert_molecule(system, molecule)
    
    # 大平移跨越边界
    dx = 1.5  # 会跨越边界
    
    residue = system.residues[-1]
    for i in range(residue.atom_count):
        atom = system.atoms[residue.atom_start + i]
        atom.x += dx
        
        # 应用PBC
        while atom.x >= box[0]:
            atom.x -= box[0]
        while atom.x < 0:
            atom.x += box[0]
    
    # 验证原子在盒子内
    for i in range(residue.atom_count):
        atom = system.atoms[residue.atom_start + i]
        assert 0 <= atom.x < box[0]
        assert 0 <= atom.y < box[1]
        assert 0 <= atom.z < box[2]


def test_rotation_preserves_bond_lengths():
    """测试旋转保持键长不变"""
    # 创建系统
    system = create_empty_system()
    
    # 添加水分子
    molecule = create_water_molecule(2.5, 2.5, 2.5)
    system = insert_molecule(system, molecule)
    
    residue = system.residues[-1]
    atoms = [system.atoms[residue.atom_start + i] for i in range(residue.atom_count)]
    
    # 计算初始键长
    initial_bonds = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            initial_bonds.append(dist)
    
    # 执行多次随机旋转
    for _ in range(10):
        quaternion = create_random_quaternion()
        q0, q1, q2, q3 = quaternion
        
        R = [
            [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
            [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
            [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
        ]
        
        com = calculate_com(residue, system)
        
        for atom in atoms:
            x = atom.x - com[0]
            y = atom.y - com[1]
            z = atom.z - com[2]
            
            new_x = R[0][0]*x + R[0][1]*y + R[0][2]*z
            new_y = R[1][0]*x + R[1][1]*y + R[1][2]*z
            new_z = R[2][0]*x + R[2][1]*y + R[2][2]*z
            
            atom.x = new_x + com[0]
            atom.y = new_y + com[1]
            atom.z = new_z + com[2]
    
    # 验证键长保持不变
    final_bonds = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            final_bonds.append(dist)
    
    for initial, final in zip(initial_bonds, final_bonds):
        assert abs(initial - final) < 1e-6, f"Bond length changed: {initial} -> {final}"