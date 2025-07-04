"""
验证 PME Total 值的来源

检查 state.ewald_energy 和函数返回值的差异
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME


def test_pme_total_sources():
    """测试 PME Total 值的一致性（回归测试）
    
    验证 state.ewald_energy 和函数返回值是否一致
    """
    
    # 创建简单系统
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # 2个原子
    positions = [[2.0, 2.5, 2.5], [3.0, 2.5, 2.5]]
    charges = [1.0, -1.0]
    
    atoms = []
    for i in range(2):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 测试不同的残基配置
    configs = [
        (1, "1个残基"),
        (2, "2个残基")
    ]
    
    for n_residues, desc in configs:
        print(f"\n测试：{desc}")
        print("-" * 50)
        
        # 设置残基
        residues = []
        if n_residues == 1:
            res = MCResidue()
            res.active = True
            res.fixed = False
            res.atomStart = 0
            res.atomCount = 2
            res.type = 0
            residues.append(res)
        else:  # 2 residues
            for i in range(2):
                res = MCResidue()
                res.active = True
                res.fixed = False
                res.atomStart = i
                res.atomCount = 1
                res.type = 0
                residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = n_residues
        
        # 初始化 PME
        initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
        
        # 调用 computeSystemEnergyPME 并获取返回值
        result = computeSystemEnergyPME(state)
        
        # result 是一个元组: (electrostatic_total, vdw, pme_dict)
        if isinstance(result, tuple) and len(result) == 3:
            elec_total, vdw, pme_dict = result
            
            print(f"\n从函数返回值获取：")
            print(f"  静电总能量: {elec_total:.2f}")
            print(f"  VDW 能量: {vdw:.6f}")
            print(f"  字典中的 total: {pme_dict.get('total', 'NOT SET')}")
            print(f"  字典中的 real_space: {pme_dict.get('real_space', 'NOT SET')}")
            print(f"  字典中的 reciprocal: {pme_dict.get('reciprocal', 'NOT SET')}")
            print(f"  字典中的 self: {pme_dict.get('self', 'NOT SET')}")
        
        # 从 state.ewald_energy 获取
        print(f"\n从 state.ewald_energy 获取：")
        print(f"  total: {state.ewald_energy.get('total', 'NOT SET')}")
        print(f"  real_space: {state.ewald_energy.get('real_space', 'NOT SET')}")
        print(f"  reciprocal: {state.ewald_energy.get('reciprocal', 'NOT SET')}")
        print(f"  self: {state.ewald_energy.get('self', 'NOT SET')}")
        
        # 手动计算
        manual_total = (state.ewald_energy.get('real_space', 0) + 
                       state.ewald_energy.get('reciprocal', 0) + 
                       state.ewald_energy.get('self', 0))
        
        print(f"\n手动计算的静电总能量: {manual_total:.2f}")
        
        # 检查差异
        if isinstance(result, tuple) and len(result) == 3:
            _, _, pme_dict = result
            dict_total = pme_dict.get('total', 0)
            state_total = state.ewald_energy.get('total', 0)
            
            print(f"\n差异分析：")
            print(f"  返回字典的 total: {dict_total:.2f}")
            print(f"  state 的 total: {state_total:.2f}")
            print(f"  差异: {abs(dict_total - state_total):.2f}")
            
            if abs(dict_total - state_total) > 1e-6:
                print("  ⚠️  返回值和 state 中的 total 不一致！")
    
    print("\n" + "="*80)


if __name__ == "__main__":
    verify_pme_total_source()