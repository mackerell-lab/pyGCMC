#!/usr/bin/env python3
"""
调试SCF算法的问题
"""

import pygcmc
import numpy as np

def analyze_scf_logic():
    """分析SCF的更新逻辑"""
    print("SCF算法逻辑分析")
    print("="*60)
    
    # 重现SCF的关键逻辑
    print("\nSCF更新公式：")
    print("-"*40)
    print("1. 计算所有粒子受到的力 F")
    print("2. 更新位置: Δr = damping * F / k")
    print("3. 新位置: r_new = r_old + Δr")
    
    # 简单例子
    print("\n数值例子：")
    print("-"*40)
    
    # 初始状态：Drude在原点
    drude_pos = 0.0
    parent_pos = 0.0
    external_pos = 1.0
    
    k_spring = 138935.0
    k_elec = 138.935
    q_drude = -1.0
    q_external = 1.0
    
    print(f"初始配置：")
    print(f"  Drude: x={drude_pos}, q={q_drude}")
    print(f"  Parent: x={parent_pos}")
    print(f"  External: x={external_pos}, q={q_external}")
    print(f"  弹簧常数: k={k_spring}")
    
    # 迭代几步
    damping = 0.5
    
    for iter in range(5):
        print(f"\n迭代 {iter+1}:")
        
        # 计算力
        # 1. 弹簧力
        displacement = drude_pos - parent_pos
        F_spring = -k_spring * displacement
        
        # 2. 电场力
        r_to_external = external_pos - drude_pos
        if abs(r_to_external) > 1e-10:
            F_electric = k_elec * q_drude * q_external / (r_to_external**2)
            # 注意：这里的力已经考虑了方向
            # 如果r_to_external > 0 (external在右边)，负电荷受到向右的力(正)
            if r_to_external > 0:
                F_electric = -F_electric  # 负电荷被吸引，力为正
        else:
            F_electric = 0.0
        
        # 总力
        F_total = F_spring + F_electric
        
        print(f"  位移: {displacement:.6f} nm = {displacement*1000:.3f} pm")
        print(f"  弹簧力: {F_spring:.1f} kJ/(mol·nm)")
        print(f"  电场力: {F_electric:.1f} kJ/(mol·nm)")
        print(f"  总力: {F_total:.1f} kJ/(mol·nm)")
        
        # SCF更新
        delta = damping * F_total / k_spring
        drude_pos_new = drude_pos + delta
        
        print(f"  更新量: Δr = {damping} * {F_total:.1f} / {k_spring} = {delta:.6f} nm")
        print(f"  新位置: {drude_pos_new:.6f} nm = {drude_pos_new*1000:.3f} pm")
        
        drude_pos = drude_pos_new
        
        # 检查收敛
        if abs(F_total) < 0.1:
            print(f"\n收敛！最终位置: {drude_pos*1000:.3f} pm")
            break
    
    # 分析问题
    print("\n\n问题分析：")
    print("="*60)
    print("SCF的核心问题在于力的计算和组合方式！")
    print("")
    print("观察上面的迭代过程，我们发现：")
    print("1. 初始时，Drude在原点")
    print("2. 电场力计算可能有符号问题")
    print("3. 更新方向可能与预期相反")

def test_scf_vs_fbp_forces():
    """对比SCF和FBP的力计算"""
    print("\n\nSCF vs FBP 力计算对比")
    print("="*60)
    
    # 创建测试系统
    atoms = []
    
    # Parent
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 0.0
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
    
    # Force对象
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
    
    # 测试不同位移下的行为
    print("\n不同位移下的力分析：")
    print(f"{'位移(pm)':<10} {'理论弹簧力':<15} {'理论电力':<15} {'理论净力':<15}")
    print("-"*55)
    
    k_spring = 138935.0
    displacements = [-2, -1, 0, 1, 2]  # pm
    
    for disp_pm in displacements:
        disp = disp_pm / 1000.0  # 转换为nm
        
        # 理论计算
        F_spring = -k_spring * disp
        r_to_ext = 1.0 - disp
        F_electric = -138.935 / (r_to_ext**2)  # 负电荷受到的力
        F_net = F_spring + F_electric
        
        print(f"{disp_pm:<10} {F_spring:>15.1f} {F_electric:>15.1f} {F_net:>15.1f}")
    
    print("\n关键洞察：")
    print("- 在位移=0时，净力应该是负的（向左），因为负电荷被右侧正电荷吸引")
    print("- 平衡点应该在正位移处（约+1 pm）")
    print("- SCF如果更新方向错误，会向相反方向移动")

def analyze_scf_code_issue():
    """分析SCF代码的具体问题"""
    print("\n\nSCF代码问题分析")
    print("="*60)
    
    print("查看SCF的更新逻辑（DrudeForce.cpp第216-217行）：")
    print("```cpp")
    print("double factor = damping / drude.kIsotropic;")
    print("delta = force * factor;")
    print("```")
    print("")
    print("这里直接使用了力来更新位置，但问题可能在于：")
    print("")
    print("1. **力的计算来源**")
    print("   - calculateForces()计算所有力")
    print("   - 包括calculateHarmonicEnergy（弹簧力）")
    print("   - 包括calculateCoulombEnergy（静电力）")
    print("")
    print("2. **可能的问题**")
    print("   - 弹簧力的符号可能有误")
    print("   - 静电力的计算可能包含了不该包含的项")
    print("   - 力的组合方式可能有问题")
    print("")
    print("3. **FBP的不同之处**")
    print("   - FBP直接求解平衡方程：F_spring + F_external = 0")
    print("   - 得到：r_d = r_p + F_external/k")
    print("   - 避免了迭代中的误差累积")

if __name__ == "__main__":
    analyze_scf_logic()
    test_scf_vs_fbp_forces()
    analyze_scf_code_issue()