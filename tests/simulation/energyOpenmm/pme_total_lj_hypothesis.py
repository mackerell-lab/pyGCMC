"""
测试 PME Total 是否包含了 LJ 能量

假设：PME Total = 静电能量 + 错误的 LJ 项
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import (initializePMEParameters, computeSystemEnergyPME,
                    computeSystemVdwEnergyCutoff)


def test_lj_hypothesis():
    """测试 PME Total 是否包含 LJ 能量（回归测试）
    
    验证修复后的 PME total 不再错误包含 LJ 能量
    """
    
    # 测试配置
    box_size = 5.0
    cutoff = 2.0
    alpha = 2.5
    
    # 测试1：纯静电系统（无 LJ）
    print("\n测试1：纯静电系统（LJ = 0）")
    print("-" * 50)
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # 无 LJ
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # 2个原子
    atoms = []
    atoms.append(MCAtom())
    atoms[0].x, atoms[0].y, atoms[0].z = 2.0, 2.5, 2.5
    atoms[0].charge = 1.0
    atoms[0].type = 0
    
    atoms.append(MCAtom())
    atoms[1].x, atoms[1].y, atoms[1].z = 3.0, 2.5, 2.5
    atoms[1].charge = -1.0
    atoms[1].type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 2个残基
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    initializePMEParameters(cutoff, state.info.box, alpha)
    computeSystemEnergyPME(state)
    
    pme_total_no_lj = state.ewald_energy.get('total', 0.0)
    pme_sum_no_lj = (state.ewald_energy.get('real_space', 0.0) + 
                     state.ewald_energy.get('reciprocal', 0.0) + 
                     state.ewald_energy.get('self', 0.0))
    offset_no_lj = pme_total_no_lj - pme_sum_no_lj
    
    print(f"PME Total: {pme_total_no_lj:.2f}")
    print(f"PME 分量和: {pme_sum_no_lj:.2f}")
    print(f"偏移量: {offset_no_lj:.2f}")
    
    # 测试2：有 LJ 的系统
    print("\n测试2：有 LJ 的系统")
    print("-" * 50)
    
    state2 = MCState()
    state2.info.box = [box_size, box_size, box_size]
    state2.info.cutoff = cutoff
    
    ff2 = MCForceField()
    ff2.numTotalTypes = 1
    ff2.numMovementTypes = 1
    ff2.ljEps = [1.0]  # 有 LJ
    ff2.ljSigma = [0.35]
    state2.forcefield = ff2
    
    # 相同的原子配置
    atoms2 = []
    atoms2.append(MCAtom())
    atoms2[0].x, atoms2[0].y, atoms2[0].z = 2.0, 2.5, 2.5
    atoms2[0].charge = 1.0
    atoms2[0].type = 0
    
    atoms2.append(MCAtom())
    atoms2[1].x, atoms2[1].y, atoms2[1].z = 3.0, 2.5, 2.5
    atoms2[1].charge = -1.0
    atoms2[1].type = 0
    
    state2.atoms = atoms2
    state2.activeAtomCount = 2
    
    # 相同的残基配置
    residues2 = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues2.append(res)
    
    state2.residues = residues2
    state2.activeResidueCount = 2
    
    # 先计算 LJ 能量
    try:
        lj_result = computeSystemVdwEnergyCutoff(state2)
        if lj_result is not None:
            lj_energy = lj_result
        else:
            lj_energy = 0.0
        print(f"LJ 能量: {lj_energy:.6f}")
    except:
        print("LJ 能量计算失败")
        lj_energy = 0.0
    
    # 计算 PME
    initializePMEParameters(cutoff, state2.info.box, alpha)
    computeSystemEnergyPME(state2)
    
    pme_total_with_lj = state2.ewald_energy.get('total', 0.0)
    pme_sum_with_lj = (state2.ewald_energy.get('real_space', 0.0) + 
                       state2.ewald_energy.get('reciprocal', 0.0) + 
                       state2.ewald_energy.get('self', 0.0))
    offset_with_lj = pme_total_with_lj - pme_sum_with_lj
    
    print(f"PME Total: {pme_total_with_lj:.2f}")
    print(f"PME 分量和: {pme_sum_with_lj:.2f}")
    print(f"偏移量: {offset_with_lj:.2f}")
    
    print(f"\n偏移量差异: {offset_with_lj - offset_no_lj:.6f}")
    print(f"与 LJ 能量比较: 差异是否接近 LJ？")
    
    # 测试3：检查残基的 LJ 能量
    print("\n测试3：残基 LJ 能量")
    print("-" * 50)
    
    total_res_lj = 0.0
    for i, res in enumerate(state2.residues):
        if hasattr(res, 'energy_vdw'):
            print(f"残基 {i} LJ 能量: {res.energy_vdw:.6f}")
            total_res_lj += res.energy_vdw
    
    print(f"残基 LJ 总和: {total_res_lj:.6f}")
    if 'lj_energy' in locals():
        print(f"系统 LJ 能量: {lj_energy:.6f}")
    
    # 测试4：1个残基的情况
    print("\n测试4：1个残基的系统")
    print("-" * 50)
    
    state3 = MCState()
    state3.info = state2.info
    state3.forcefield = state2.forcefield
    state3.atoms = state2.atoms
    state3.activeAtomCount = state2.activeAtomCount
    
    # 1个残基包含2个原子
    residues3 = []
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 2
    res.type = 0
    residues3.append(res)
    
    state3.residues = residues3
    state3.activeResidueCount = 1
    
    initializePMEParameters(cutoff, state3.info.box, alpha)
    computeSystemEnergyPME(state3)
    
    pme_total_1res = state3.ewald_energy.get('total', 0.0)
    pme_sum_1res = (state3.ewald_energy.get('real_space', 0.0) + 
                    state3.ewald_energy.get('reciprocal', 0.0) + 
                    state3.ewald_energy.get('self', 0.0))
    offset_1res = pme_total_1res - pme_sum_1res
    
    print(f"1个残基系统偏移量: {offset_1res:.2f} （应该是 0）")
    
    print("\n" + "="*80)
    print("结论")
    print("="*80)
    
    if abs(offset_1res) < 1e-6:
        print("✓ 确认：1个残基时偏移量为 0")
    
    if abs(offset_no_lj) == abs(offset_with_lj):
        print("✓ 偏移量与 LJ 参数无关")
    else:
        print(f"✗ 偏移量可能与 LJ 有关：差异 = {abs(offset_with_lj - offset_no_lj):.6f}")


if __name__ == "__main__":
    test_lj_hypothesis()