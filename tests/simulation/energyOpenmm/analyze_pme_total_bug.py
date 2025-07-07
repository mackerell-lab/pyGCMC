"""
深入分析 PME Total 字段的 bug

找出偏移量的规律和来源
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

def test_pme_total_configurations():
    """测试 PME Total 在不同配置下的正确性（回归测试）
    
    验证修复后的 PME total 在各种原子/残基配置下都正确计算
    为每个测试案例创建新的 MCState 实例以避免状态污染
    """
    
    # 固定参数
    box_size = 5.0
    cutoff = 2.0
    alpha = 2.5
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # ------------------------------------------------------------
    # 只在函数开头初始化一次 PME（设好 α / 网格 / 样条阶数）
    # 后续针对不同 MCState 直接调用 computeSystemEnergyPME 即可复用同一组
    # 全局参数，避免重复释放/重建网格带来的内存问题
    # ------------------------------------------------------------
    try:
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(
            cutoff,
            [box_size, box_size, box_size],
            alpha,
            mesh_size,
            spline_order
        )
    except Exception as e:
        # 如果初始化失败，可能是因为已经初始化过
        print(f"PME initialization warning: {e}")
        pass
    
    # 测试案例1：2个原子（电中性对）
    print("测试案例1：2个原子（电中性对）")
    
    # 创建新的 state
    state1 = MCState()
    state1.info.box = [box_size, box_size, box_size]
    state1.info.cutoff = cutoff
    
    ff1 = MCForceField()
    ff1.numTotalTypes = 1
    ff1.numMovementTypes = 1
    ff1.ljEps = [0.0]
    ff1.ljSigma = [0.3]
    state1.forcefield = ff1
    
    # 创建2个原子
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 2.0, 2.0, 2.5
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 3.0, 2.0, 2.5
    atom2.charge = -1.0
    atom2.type = 0
    
    state1.atoms = [atom1, atom2]
    state1.activeAtomCount = 2
    
    # 2个残基
    res1_1 = MCResidue()
    res1_1.active = True
    res1_1.fixed = False
    res1_1.atomStart = 0
    res1_1.atomCount = 1
    res1_1.type = 0
    
    res1_2 = MCResidue()
    res1_2.active = True
    res1_2.fixed = False
    res1_2.atomStart = 1
    res1_2.atomCount = 1
    res1_2.type = 0
    
    state1.residues = [res1_1, res1_2]
    state1.activeResidueCount = 2
    
    # 直接计算能量（PME 已在函数开头初始化）
    pygcmc.computeSystemEnergyPME(state1)
    
    pme_total1 = state1.ewald_energy.get('total', 0.0)
    pme_sum1 = (state1.ewald_energy.get('real_space', 0.0) + 
                state1.ewald_energy.get('reciprocal', 0.0) + 
                state1.ewald_energy.get('self', 0.0))
    offset1 = pme_total1 - pme_sum1
    
    print(f"  Total: {pme_total1:.6f}, Sum: {pme_sum1:.6f}, Offset: {offset1:.6f}")
    assert abs(offset1) < 1e-6, f"2原子系统的偏移量应该为 0，实际为 {offset1}"
    
    # 测试案例2：3个原子
    print("\n测试案例2：3个原子")
    
    # 创建新的 state
    state2 = MCState()
    state2.info.box = [box_size, box_size, box_size]
    state2.info.cutoff = cutoff
    
    ff2 = MCForceField()
    ff2.numTotalTypes = 1
    ff2.numMovementTypes = 1
    ff2.ljEps = [0.0]
    ff2.ljSigma = [0.3]
    state2.forcefield = ff2
    
    # 创建3个原子
    positions2 = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [2.5, 3.0, 2.5]]
    charges2 = [1.0, -1.0, 0.5]
    
    atoms2 = []
    for i in range(3):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions2[i]
        atom.charge = charges2[i]
        atom.type = 0
        atoms2.append(atom)
    
    state2.atoms = atoms2
    state2.activeAtomCount = 3
    
    # 3个残基
    residues2 = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues2.append(res)
    
    state2.residues = residues2
    state2.activeResidueCount = 3
    
    # Calculate energy directly, no need to reinitialize
    pygcmc.computeSystemEnergyPME(state2)
    
    pme_total2 = state2.ewald_energy.get('total', 0.0)
    pme_sum2 = (state2.ewald_energy.get('real_space', 0.0) + 
                state2.ewald_energy.get('reciprocal', 0.0) + 
                state2.ewald_energy.get('self', 0.0))
    offset2 = pme_total2 - pme_sum2
    
    print(f"  Total: {pme_total2:.6f}, Sum: {pme_sum2:.6f}, Offset: {offset2:.6f}")
    assert abs(offset2) < 1e-6, f"3原子系统的偏移量应该为 0，实际为 {offset2}"
    
    # 测试案例3：4原子系统（1个4原子残基）
    print("\n测试案例3：4原子系统（1个4原子残基）")
    
    # 创建新的 state
    state3 = MCState()
    state3.info.box = [box_size, box_size, box_size]
    state3.info.cutoff = cutoff
    
    ff3 = MCForceField()
    ff3.numTotalTypes = 1
    ff3.numMovementTypes = 1
    ff3.ljEps = [0.0]
    ff3.ljSigma = [0.3]
    state3.forcefield = ff3
    
    # 创建4个原子
    positions3 = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]]
    charges3 = [1.0, -1.0, 1.0, -1.0]
    
    atoms3 = []
    for i in range(4):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions3[i]
        atom.charge = charges3[i]
        atom.type = 0
        atoms3.append(atom)
    
    state3.atoms = atoms3
    state3.activeAtomCount = 4
    
    # 1个4原子残基
    res3 = MCResidue()
    res3.active = True
    res3.fixed = False
    res3.atomStart = 0
    res3.atomCount = 4
    res3.type = 0
    
    state3.residues = [res3]
    state3.activeResidueCount = 1
    
    # Calculate energy directly, no need to reinitialize
    pygcmc.computeSystemEnergyPME(state3)
    
    pme_total3 = state3.ewald_energy.get('total', 0.0)
    pme_sum3 = (state3.ewald_energy.get('real_space', 0.0) + 
                state3.ewald_energy.get('reciprocal', 0.0) + 
                state3.ewald_energy.get('self', 0.0))
    offset3 = pme_total3 - pme_sum3
    
    print(f"  Total: {pme_total3:.6f}, Sum: {pme_sum3:.6f}, Offset: {offset3:.6f}")
    assert abs(offset3) < 1e-6, f"1个4原子残基的偏移量应该为 0，实际为 {offset3}"
    
    # 验证所有偏移量都为 0
    print(f"\n所有测试通过！偏移量: {offset1:.6f}, {offset2:.6f}, {offset3:.6f}")
    assert abs(offset1) < 1e-6 and abs(offset2) < 1e-6 and abs(offset3) < 1e-6, \
        "PME Total 偏移量应该在所有配置下都为 0"


if __name__ == "__main__":
    test_pme_total_configurations()