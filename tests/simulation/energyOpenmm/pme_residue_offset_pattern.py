"""
精确测试 PME Total 偏移量与残基的关系

假设：偏移量 = f(残基数量)
"""

import numpy as np
import sys
import os
import subprocess

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def compute_pme_with_state_isolation(n_residues, residue_config):
    """在隔离的子进程中计算PME能量，避免全局状态污染"""
    
    # 构建Python脚本字符串，直接返回结果
    script = f"""
import sys
sys.path.insert(0, '{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}')

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

# 固定系统参数
box_size = 5.0
cutoff = 2.0
alpha = 2.5
mesh_size = [32, 32, 32]
spline_order = 4

# 创建新的 state
state = MCState()
state.info.box = [box_size, box_size, box_size]
state.info.cutoff = cutoff

ff = MCForceField()
ff.numTotalTypes = 1
ff.numMovementTypes = 1
ff.ljEps = [0.0]
ff.ljSigma = [0.3]
state.forcefield = ff

# 共享的原子位置和电荷
positions = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]]
charges = [1.0, -1.0, 1.0, -1.0]
n_atoms = 4

# 创建原子
atoms = []
for i in range(n_atoms):
    atom = MCAtom()
    atom.x, atom.y, atom.z = positions[i]
    atom.charge = charges[i]
    atom.type = 0
    atoms.append(atom)

state.atoms = atoms
state.activeAtomCount = n_atoms

# 创建残基配置
residues = []
residue_config = {residue_config}
for res_info in residue_config:
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = res_info['atomStart']
    res.atomCount = res_info['atomCount']
    res.type = 0
    residues.append(res)

state.residues = residues
state.activeResidueCount = len(residues)

# 初始化PME
pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
pygcmc.initializePMEParameters(
    cutoff,
    [box_size, box_size, box_size],
    alpha,
    mesh_size,
    spline_order
)

# 计算能量
pygcmc.computeSystemEnergyPME(state)

# 直接打印结果字典
print({{
    'total': state.ewald_energy.get('total', 0.0),
    'real_space': state.ewald_energy.get('real_space', 0.0),
    'reciprocal': state.ewald_energy.get('reciprocal', 0.0),
    'self': state.ewald_energy.get('self', 0.0)
}})
"""

    # 在子进程中运行
    try:
        result = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, 'PYTHONPATH': os.environ.get('PYTHONPATH', '')}
        )
        
        # 使用eval解析字典（安全因为我们控制输出）
        return eval(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e.stderr}")
        raise


def test_residue_offset_pattern():
    """测试 PME Total 不应该有与残基数量相关的偏移（回归测试）
    
    验证修复后的 PME total 在不同残基配置下都正确
    为每个测试案例创建新的 MCState 实例以避免状态污染
    """
    
    # 测试案例1：1个4原子残基
    print("测试案例1：1个4原子残基")
    residue_config1 = [{'atomStart': 0, 'atomCount': 4}]
    result1 = compute_pme_with_state_isolation(1, residue_config1)
    
    total1 = result1['total']
    sum1 = result1['real_space'] + result1['reciprocal'] + result1['self']
    offset1 = total1 - sum1
    
    print(f"  Total: {total1:.6f}, Sum: {sum1:.6f}, Offset: {offset1:.6f}")
    assert abs(offset1) < 1e-6, f"1个残基配置的偏移应该为0，实际为 {offset1}"
    
    # 测试案例2：2个2原子残基
    print("\n测试案例2：2个2原子残基")
    residue_config2 = [
        {'atomStart': 0, 'atomCount': 2},
        {'atomStart': 2, 'atomCount': 2}
    ]
    result2 = compute_pme_with_state_isolation(2, residue_config2)
    
    total2 = result2['total']
    sum2 = result2['real_space'] + result2['reciprocal'] + result2['self']
    offset2 = total2 - sum2
    
    print(f"  Total: {total2:.6f}, Sum: {sum2:.6f}, Offset: {offset2:.6f}")
    assert abs(offset2) < 1e-6, f"2个残基配置的偏移应该为0，实际为 {offset2}"
    
    # 测试案例3：4个1原子残基
    print("\n测试案例3：4个1原子残基")
    residue_config3 = [
        {'atomStart': 0, 'atomCount': 1},
        {'atomStart': 1, 'atomCount': 1},
        {'atomStart': 2, 'atomCount': 1},
        {'atomStart': 3, 'atomCount': 1}
    ]
    result3 = compute_pme_with_state_isolation(4, residue_config3)
    
    total3 = result3['total']
    sum3 = result3['real_space'] + result3['reciprocal'] + result3['self']
    offset3 = total3 - sum3
    
    print(f"  Total: {total3:.6f}, Sum: {sum3:.6f}, Offset: {offset3:.6f}")
    assert abs(offset3) < 1e-6, f"4个残基配置的偏移应该为0，实际为 {offset3}"
    
    # 验证所有偏移都为 0
    print(f"\n所有测试通过！偏移量: {offset1:.6f}, {offset2:.6f}, {offset3:.6f}")
    assert abs(offset1) < 1e-6 and abs(offset2) < 1e-6 and abs(offset3) < 1e-6, \
        "PME Total 修复后不应该有与残基数量相关的偏移"


if __name__ == "__main__":
    test_residue_offset_pattern()