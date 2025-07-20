#!/usr/bin/env python3
"""
测试SCF中各种力的贡献
"""

import pygcmc
import numpy as np

def test_force_contributions():
    """测试各种力的贡献"""
    print("SCF力贡献分析")
    print("="*60)
    
    # 创建最简单的系统
    atoms = []
    
    # Parent
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 0.0  # 不带电
    parent.type = 0
    atoms.append(parent)
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    atoms.append(drude)
    
    # External
    external = pygcmc.MCAtom()
    external.x = 1.0
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.0
    external.type = 2
    atoms.append(external)
    
    # 残基
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 1
    res2.active = True
    res2.type = 1
    
    residues = [res1, res2]
    
    # 状态
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 3
    state.activeResidueCount = 2
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    state.forcefield.numTotalTypes = 3
    state.forcefield.numMovementTypes = 3
    state.forcefield.ljSigma = [0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0]
    
    # DrudeForce
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.0,
        polarizability=0.001,
        aniso12=1.0, aniso34=1.0
    )
    
    # 手动计算预期的力
    print("\n手动计算预期的力：")
    print("-"*40)
    
    k_spring = 138935.0
    k_elec = 138.935
    
    # 初始状态：Drude在原点
    disp = 0.0
    
    # 1. 弹簧力
    F_spring = -k_spring * disp  # 0
    print(f"弹簧力: F_spring = -k*d = {F_spring}")
    
    # 2. 库仑力（Drude-External）
    r_DE = 1.0  # 距离
    # 按照我们的代码逻辑
    delta_x = 1.0 - 0.0  # external.x - drude.x
    forceMag = k_elec * (-1.0) * 1.0 / r_DE**2  # -138.935
    F_coulomb_code = delta_x * (forceMag / r_DE)  # 1.0 * (-138.935) = -138.935
    print(f"库仑力（代码）: F = {F_coulomb_code}")
    
    # 物理上正确的库仑力
    # 负电荷被正电荷吸引，力应该向右
    F_coulomb_physics = 138.935  # 向右
    print(f"库仑力（物理）: F = {F_coulomb_physics}")
    
    print(f"\n总力：")
    print(f"代码预测: F_total = {F_spring + F_coulomb_code}")
    print(f"物理预期: F_total = {F_spring + F_coulomb_physics}")
    
    # 测试不同的算法
    print("\n\n实际测试结果：")
    print("="*60)
    
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}算法：")
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001
        params.maxIterations = 10  # 限制迭代次数以便观察
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        # 计算
        energy = force.calculateEnergySCF(state)
        
        disp = state.atoms[1].x
        print(f"  最终位移: {disp*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")
        
        # 分析更新方向
        if disp < 0:
            print(f"  Drude向左移动 ← （远离正电荷）")
        else:
            print(f"  Drude向右移动 → （靠近正电荷）")

def test_with_parent_charge():
    """测试Parent带电的情况"""
    print("\n\n测试Parent带电的情况")
    print("="*60)
    
    # 创建系统（Parent也带电）
    atoms = []
    
    # Parent（带正电）
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 1.71636  # 水模型的实际电荷
    parent.type = 0
    atoms.append(parent)
    
    # Drude（带负电）
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    # External（带正电）
    external = pygcmc.MCAtom()
    external.x = 1.0
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.71636
    external.type = 2
    atoms.append(external)
    
    # 残基设置...（省略重复代码）
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 1
    res2.active = True
    res2.type = 1
    
    residues = [res1, res2]
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 3
    state.activeResidueCount = 2
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    state.forcefield.numTotalTypes = 3
    state.forcefield.numMovementTypes = 3
    state.forcefield.ljSigma = [0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0]
    
    # Force
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.71636,
        polarizability=0.000978253,
        aniso12=1.0, aniso34=1.0
    )
    
    print("\n注意：Parent-Drude相互作用应该被排除（分子内）")
    print("只有External-Drude和External-Parent的相互作用")
    
    # 测试
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.001
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}结果：")
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置
        state.atoms[1].x = 0.0
        
        # 计算
        energy = force.calculateEnergySCF(state)
        
        disp = state.atoms[1].x
        print(f"  位移: {disp*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")

if __name__ == "__main__":
    test_force_contributions()
    test_with_parent_charge()