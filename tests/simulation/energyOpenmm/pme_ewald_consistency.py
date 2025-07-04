"""
测试 PME 和 Ewald 的一致性

理论上 PME 是 Ewald 的快速实现，结果应该高度一致
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import (initializePMEParameters, initializeEwaldParameters,
                    computeSystemEnergyPME, computeSystemEnergyEwald)


def test_pme_ewald_consistency():
    """测试 PME 和 Ewald 的一致性"""
    
    print("\n" + "="*80)
    print("PME vs Ewald 一致性测试")
    print("="*80)
    
    # 测试不同的系统配置
    test_cases = [
        {
            'name': '两个带电粒子',
            'positions': [[2.0, 2.0, 2.0], [3.0, 2.0, 2.0]],
            'charges': [1.0, -1.0],
            'box_size': 5.0,
            'cutoff': 1.8
        },
        {
            'name': '四个带电粒子（正方形）',
            'positions': [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], 
                         [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]],
            'charges': [1.0, -1.0, 1.0, -1.0],
            'box_size': 5.0,
            'cutoff': 1.8
        },
        {
            'name': '六个带电粒子（随机）',
            'positions': [[1.5, 1.5, 2.5], [3.5, 1.5, 2.5],
                         [3.5, 3.5, 2.5], [1.5, 3.5, 2.5],
                         [2.5, 2.5, 1.5], [2.5, 2.5, 3.5]],
            'charges': [1.0, -1.0, 1.0, -1.0, 0.5, -0.5],
            'box_size': 5.0,
            'cutoff': 2.0
        }
    ]
    
    for test_case in test_cases:
        print(f"\n测试案例：{test_case['name']}")
        print("-" * 60)
        
        # 创建系统
        state = MCState()
        box_size = test_case['box_size']
        cutoff = test_case['cutoff']
        
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = cutoff
        
        # 设置力场（无 LJ）
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # 添加原子
        positions = test_case['positions']
        charges = test_case['charges']
        n_atoms = len(positions)
        
        atoms = []
        for i in range(n_atoms):
            atom = MCAtom()
            atom.x, atom.y, atom.z = positions[i]
            atom.charge = charges[i]
            atom.type = 0
            atoms.append(atom)
        
        state.atoms = atoms
        state.activeAtomCount = n_atoms
        
        # 添加残基
        residues = []
        for i in range(n_atoms):
            res = MCResidue()
            res.active = True
            res.fixed = False
            res.atomStart = i
            res.atomCount = 1
            res.type = 0
            residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = n_atoms
        
        # 测试不同的 alpha 值
        alpha_values = [2.0, 2.5, 3.0]
        
        for alpha in alpha_values:
            print(f"\n  Alpha = {alpha:.1f}")
            
            # 初始化 PME
            mesh_size = [32, 32, 32]
            initializePMEParameters(cutoff, state.info.box, alpha)
            
            # 初始化 Ewald
            initializeEwaldParameters(cutoff, state.info.box, alpha)
            
            # 计算 PME 能量
            computeSystemEnergyPME(state)
            pme_total = state.ewald_energy.get('total', 0.0)
            pme_real = state.ewald_energy.get('real_space', 0.0)
            pme_recip = state.ewald_energy.get('reciprocal', 0.0)
            pme_self = state.ewald_energy.get('self', 0.0)
            
            # 计算 Ewald 能量
            ewald_result = computeSystemEnergyEwald(state)
            ewald_total = ewald_result[0]
            ewald_dict = ewald_result[2]
            ewald_real = ewald_dict.get('real_space', 0.0)
            ewald_recip = ewald_dict.get('reciprocal', 0.0)
            ewald_self = ewald_dict.get('self', 0.0)
            
            # 比较结果
            print(f"    PME:   Real={pme_real:10.4f}, Recip={pme_recip:10.4f}, "
                  f"Self={pme_self:10.4f}, Total={pme_total:10.4f}")
            print(f"    Ewald: Real={ewald_real:10.4f}, Recip={ewald_recip:10.4f}, "
                  f"Self={ewald_self:10.4f}, Total={ewald_total:10.4f}")
            
            # 计算差异
            diff_real = abs(pme_real - ewald_real)
            diff_recip = abs(pme_recip - ewald_recip)
            diff_self = abs(pme_self - ewald_self)
            diff_total = abs(pme_total - ewald_total)
            
            # 计算相对误差
            if abs(ewald_real) > 1e-6:
                rel_err_real = diff_real / abs(ewald_real) * 100
            else:
                rel_err_real = 0
                
            if abs(ewald_recip) > 1e-6:
                rel_err_recip = diff_recip / abs(ewald_recip) * 100
            else:
                rel_err_recip = 0
                
            if abs(ewald_self) > 1e-6:
                rel_err_self = diff_self / abs(ewald_self) * 100
            else:
                rel_err_self = 0
                
            if abs(ewald_total) > 1e-6:
                rel_err_total = diff_total / abs(ewald_total) * 100
            else:
                rel_err_total = 0
            
            print(f"    差异:  Real={rel_err_real:6.2f}%, Recip={rel_err_recip:6.2f}%, "
                  f"Self={rel_err_self:6.2f}%, Total={rel_err_total:6.2f}%")
            
            # 检查 PME total 是否等于分量之和
            pme_sum = pme_real + pme_recip + pme_self
            print(f"    PME 分量之和: {pme_sum:10.4f} (应该等于 Total)")
            if abs(pme_total - pme_sum) > 1e-6:
                print(f"    ⚠️  PME Total 与分量之和不一致！差异: {abs(pme_total - pme_sum):.2f}")
    
    # 测试能量变化的一致性
    print("\n" + "="*80)
    print("能量变化（ΔE）一致性测试")
    print("="*80)
    
    # 使用第二个测试案例
    test_case = test_cases[1]
    
    # 创建系统
    state = MCState()
    box_size = test_case['box_size']
    cutoff = test_case['cutoff']
    
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    positions = np.array(test_case['positions'])
    charges = test_case['charges']
    n_atoms = len(positions)
    
    atoms = []
    for i in range(n_atoms):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = n_atoms
    
    residues = []
    for i in range(n_atoms):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = n_atoms
    
    # 初始化
    alpha = 2.5
    initializePMEParameters(cutoff, state.info.box, alpha)
    initializeEwaldParameters(cutoff, state.info.box, alpha)
    
    # 计算初始能量
    computeSystemEnergyPME(state)
    initial_pme = state.ewald_energy.get('reciprocal', 0.0)
    
    ewald_result = computeSystemEnergyEwald(state)
    initial_ewald = ewald_result[2].get('reciprocal', 0.0)
    
    # 移动一个原子
    state.atoms[0].x += 0.2
    state.atoms[0].y += 0.1
    
    # 计算移动后能量
    computeSystemEnergyPME(state)
    moved_pme = state.ewald_energy.get('reciprocal', 0.0)
    
    ewald_result = computeSystemEnergyEwald(state)
    moved_ewald = ewald_result[2].get('reciprocal', 0.0)
    
    # 计算 ΔE
    delta_pme = moved_pme - initial_pme
    delta_ewald = moved_ewald - initial_ewald
    
    print(f"\n初始倒空间能量:")
    print(f"  PME:   {initial_pme:10.4f} kJ/mol")
    print(f"  Ewald: {initial_ewald:10.4f} kJ/mol")
    
    print(f"\n移动后倒空间能量:")
    print(f"  PME:   {moved_pme:10.4f} kJ/mol")
    print(f"  Ewald: {moved_ewald:10.4f} kJ/mol")
    
    print(f"\n能量变化 ΔE:")
    print(f"  PME:   {delta_pme:10.4f} kJ/mol")
    print(f"  Ewald: {delta_ewald:10.4f} kJ/mol")
    print(f"  差异:  {abs(delta_pme - delta_ewald):10.6f} kJ/mol")
    
    if abs(delta_ewald) > 1e-6:
        rel_err = abs(delta_pme - delta_ewald) / abs(delta_ewald) * 100
        print(f"  相对误差: {rel_err:.2f}%")
    
    print("\n" + "="*80)
    print("结论")
    print("="*80)
    
    if abs(delta_pme - delta_ewald) / abs(delta_ewald) > 0.05:
        print("⚠️  PME 和 Ewald 的能量变化差异超过 5%")
        print("可能存在 PME 实现问题")
    else:
        print("✓ PME 和 Ewald 的能量变化一致性良好")


if __name__ == "__main__":
    test_pme_ewald_consistency()