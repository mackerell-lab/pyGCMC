"""
精确测试 PME Total 偏移量与残基的关系

假设：偏移量 = f(残基数量)
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME


@pytest.mark.skip(reason="PME global state causes segfault with multiple MCState instances - needs C++ fix")
def test_residue_offset_pattern():
    """测试 PME Total 不应该有与残基数量相关的偏移（回归测试）
    
    验证修复后的 PME total 在不同残基配置下都正确
    """
    
    # 固定系统参数
    box_size = 5.0
    cutoff = 2.0
    alpha = 2.5
    
    # 只初始化一次 PME 参数
    initializePMEParameters(cutoff, [box_size, box_size, box_size], alpha)
    
    # 使用4个原子的系统
    positions = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]]
    charges = [1.0, -1.0, 1.0, -1.0]
    n_atoms = 4
    
    # 测试案例1：1个4原子残基
    state1 = MCState()
    state1.info.box = [box_size, box_size, box_size]
    state1.info.cutoff = cutoff
    
    ff1 = MCForceField()
    ff1.numTotalTypes = 1
    ff1.numMovementTypes = 1
    ff1.ljEps = [0.0]
    ff1.ljSigma = [0.3]
    state1.forcefield = ff1
    
    atoms1 = []
    for i in range(n_atoms):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms1.append(atom)
    
    state1.atoms = atoms1
    state1.activeAtomCount = n_atoms
    
    res1 = MCResidue()
    res1.active = True
    res1.fixed = False
    res1.atomStart = 0
    res1.atomCount = 4
    res1.type = 0
    
    state1.residues = [res1]
    state1.activeResidueCount = 1
    
    computeSystemEnergyPME(state1)
    
    total1 = state1.ewald_energy.get('total', 0.0)
    sum1 = (state1.ewald_energy.get('real_space', 0.0) + 
            state1.ewald_energy.get('reciprocal', 0.0) + 
            state1.ewald_energy.get('self', 0.0))
    offset1 = total1 - sum1
    
    assert abs(offset1) < 1e-6, f"1个残基配置的偏移应该为0，实际为 {offset1}"
    
    # 测试案例2：2个2原子残基
    state2 = MCState()
    state2.info.box = [box_size, box_size, box_size]
    state2.info.cutoff = cutoff
    
    ff2 = MCForceField()
    ff2.numTotalTypes = 1
    ff2.numMovementTypes = 1
    ff2.ljEps = [0.0]
    ff2.ljSigma = [0.3]
    state2.forcefield = ff2
    
    atoms2 = []
    for i in range(n_atoms):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms2.append(atom)
    
    state2.atoms = atoms2
    state2.activeAtomCount = n_atoms
    
    residues2 = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i * 2
        res.atomCount = 2
        res.type = 0
        residues2.append(res)
    
    state2.residues = residues2
    state2.activeResidueCount = 2
    
    computeSystemEnergyPME(state2)
    
    total2 = state2.ewald_energy.get('total', 0.0)
    sum2 = (state2.ewald_energy.get('real_space', 0.0) + 
            state2.ewald_energy.get('reciprocal', 0.0) + 
            state2.ewald_energy.get('self', 0.0))
    offset2 = total2 - sum2
    
    assert abs(offset2) < 1e-6, f"2个残基配置的偏移应该为0，实际为 {offset2}"
    
    # 测试案例3：4个1原子残基
    state3 = MCState()
    state3.info.box = [box_size, box_size, box_size]
    state3.info.cutoff = cutoff
    
    ff3 = MCForceField()
    ff3.numTotalTypes = 1
    ff3.numMovementTypes = 1
    ff3.ljEps = [0.0]
    ff3.ljSigma = [0.3]
    state3.forcefield = ff3
    
    atoms3 = []
    for i in range(n_atoms):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms3.append(atom)
    
    state3.atoms = atoms3
    state3.activeAtomCount = n_atoms
    
    residues3 = []
    for i in range(4):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues3.append(res)
    
    state3.residues = residues3
    state3.activeResidueCount = 4
    
    computeSystemEnergyPME(state3)
    
    total3 = state3.ewald_energy.get('total', 0.0)
    sum3 = (state3.ewald_energy.get('real_space', 0.0) + 
            state3.ewald_energy.get('reciprocal', 0.0) + 
            state3.ewald_energy.get('self', 0.0))
    offset3 = total3 - sum3
    
    assert abs(offset3) < 1e-6, f"4个残基配置的偏移应该为0，实际为 {offset3}"
    
    # 验证所有偏移都为 0
    assert abs(offset1) < 1e-6 and abs(offset2) < 1e-6 and abs(offset3) < 1e-6, \
        "PME Total 修复后不应该有与残基数量相关的偏移"


if __name__ == "__main__":
    test_residue_offset_pattern()