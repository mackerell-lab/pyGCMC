#!/usr/bin/env python3
"""
简化测试：理解SCF和FBP的差异
"""

import pygcmc
import numpy as np

def test_simple_system():
    """测试最简单的系统"""
    print("简化测试：SCF vs FBP")
    print("="*60)
    
    # 创建极简系统
    atoms = []
    
    # 情况1：Parent不带电
    print("\n情况1：Parent不带电，只有Drude和外部电荷相互作用")
    print("-"*40)
    
    # Parent (charge = 0)
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 0.0  # 不带电！
    parent.type = 0
    atoms.append(parent)
    
    # Drude (charge = -1)
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    atoms.append(drude)
    
    # External positive charge
    external = pygcmc.MCAtom()
    external.x = 1.0
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.0
    external.type = 2
    atoms.append(external)
    
    # 残基设置
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
    
    # 创建状态
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
    
    # 设置参数
    k_spring = 138935.0  # 简单值
    polarizability = 0.001
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.0,
        polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # 理论预测
    print("理论分析:")
    print(f"  弹簧常数 k = {k_spring} kJ/(mol·nm²)")
    print(f"  外部电场 E = 138.935 kJ/(mol·nm·e) (向左)")
    print(f"  Drude受力 F = -1 * E = -138.935 kJ/(mol·nm) (向左)")
    print(f"  平衡位移 d = F/k = -0.001 nm = -1.0 pm")
    print(f"  结论：Drude应该向左移动1.0 pm")
    
    # 测试算法
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}结果:")
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001
        params.maxIterations = 100
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
        
        disp = state.atoms[1].x - state.atoms[0].x
        print(f"  位移: {disp*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")
        
        # 验证力
        F_spring = -k_spring * disp
        r = 1.0 - disp  # 到外部电荷的距离
        F_electric = -138.935 / (r*r) * (-1.0)  # 负电荷受到的力
        print(f"  弹簧力: {F_spring:.1f} kJ/(mol·nm)")
        print(f"  电力: {F_electric:.1f} kJ/(mol·nm)")
        print(f"  净力: {F_spring + F_electric:.1f} kJ/(mol·nm)")
    
    # 情况2：真实水模型参数
    print("\n\n情况2：真实水模型参数")
    print("-"*40)
    
    # 更新电荷
    state.atoms[0].charge = 1.71636   # Parent带正电
    state.atoms[1].charge = -1.71636  # Drude带负电
    state.atoms[2].charge = 1.71636   # External带正电
    
    # 更新力常数
    k_spring_real = 418400.0
    polarizability_real = 0.000978253
    
    force2 = pygcmc.DrudeForce()
    force2.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.71636,
        polarizability=polarizability_real,
        aniso12=1.0, aniso34=1.0
    )
    
    print("系统参数:")
    print(f"  Parent: +1.71636 at (0,0,0)")
    print(f"  Drude: -1.71636 at (0,0,0)")
    print(f"  External: +1.71636 at (1,0,0)")
    print(f"  k = {k_spring_real} kJ/(mol·nm²)")
    
    for algo in algorithms:
        print(f"\n{algo}结果:")
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force2.setSCFParameters(params)
        
        if algo == "SCF":
            force2.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force2.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        # 计算
        energy = force2.calculateEnergySCF(state)
        
        disp = state.atoms[1].x - state.atoms[0].x
        print(f"  位移: {disp*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")

if __name__ == "__main__":
    test_simple_system()