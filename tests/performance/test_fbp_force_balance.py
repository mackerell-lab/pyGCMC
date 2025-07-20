#!/usr/bin/env python3
"""
测试FBP的力平衡逻辑
"""

import pygcmc
import numpy as np

def test_force_balance_logic():
    """测试力平衡的基本逻辑"""
    print("测试力平衡基本逻辑")
    print("="*60)
    
    # 创建简单系统：一个原子带正电，一个Drude带负电
    atoms = []
    
    # Parent原子
    atom = pygcmc.MCAtom()
    atom.x = 0.0
    atom.y = 0.0
    atom.z = 0.0
    atom.charge = 0.0
    atom.type = 0
    atoms.append(atom)
    
    # Drude粒子
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    atoms.append(drude)
    
    # 外部正电荷
    external = pygcmc.MCAtom()
    external.x = 1.0
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.0
    external.type = 2
    atoms.append(external)
    
    # 创建残基
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
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    state.forcefield.numTotalTypes = 3
    state.forcefield.numMovementTypes = 3
    state.forcefield.ljSigma = [0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0]
    
    # 理论分析
    print("\n理论分析:")
    print("-"*40)
    
    # 参数
    q_drude = -1.0
    q_external = 1.0
    r = 1.0  # 距离
    k_elec = 138.935  # 静电常数
    polarizability = 0.001  # nm^3
    k_spring = q_drude**2 * k_elec / polarizability  # 弹簧常数
    
    print(f"Drude电荷: {q_drude} e")
    print(f"外部电荷: {q_external} e")
    print(f"距离: {r} nm")
    print(f"弹簧常数: k = {k_spring:.1f} kJ/(mol·nm²)")
    
    # 力平衡分析
    print("\n力平衡分析:")
    print("-"*40)
    
    # 电场（在Drude位置，由外部电荷产生）
    # E = k*q/r² 指向远离正电荷
    E = k_elec * q_external / r**2
    print(f"电场强度: E = {E:.1f} kJ/(mol·nm·e)")
    
    # Drude受到的电力
    F_electric = q_drude * E
    print(f"电力: F_electric = q*E = {F_electric:.1f} kJ/(mol·nm)")
    print(f"  方向: {'向左(远离正电荷)' if F_electric < 0 else '向右(朝向正电荷)'}")
    
    # 平衡位置
    # 力平衡: F_spring + F_electric = 0
    # -k*Δx + F_electric = 0
    # Δx = F_electric/k
    displacement = F_electric / k_spring
    print(f"\n平衡位移: Δx = F_electric/k = {displacement:.6f} nm = {displacement*1000:.3f} pm")
    print(f"  方向: {'向左' if displacement < 0 else '向右'}")
    
    # 验证力平衡
    F_spring = -k_spring * displacement
    print(f"\n验证力平衡:")
    print(f"  弹簧力: F_spring = -k*Δx = {F_spring:.1f} kJ/(mol·nm)")
    print(f"  电力: F_electric = {F_electric:.1f} kJ/(mol·nm)")
    print(f"  总力: F_total = {F_spring + F_electric:.3f} kJ/(mol·nm)")
    
    # 能量计算
    print("\n能量分析:")
    print("-"*40)
    
    # 谐振子能量
    E_harmonic = 0.5 * k_spring * displacement**2
    print(f"谐振子能量: E_harmonic = 0.5*k*Δx² = {E_harmonic:.6f} kJ/mol")
    
    # 静电能量变化
    # 初始位置(0,0,0)到最终位置(Δx,0,0)
    r_initial = r
    r_final = r - displacement  # Drude向左移动，距离减小
    E_coulomb_initial = k_elec * q_drude * q_external / r_initial
    E_coulomb_final = k_elec * q_drude * q_external / r_final
    ΔE_coulomb = E_coulomb_final - E_coulomb_initial
    
    print(f"初始静电能: {E_coulomb_initial:.6f} kJ/mol")
    print(f"最终静电能: {E_coulomb_final:.6f} kJ/mol")
    print(f"静电能变化: ΔE_coulomb = {ΔE_coulomb:.6f} kJ/mol")
    
    print(f"\n总能量变化: ΔE_total = {E_harmonic + ΔE_coulomb:.6f} kJ/mol")
    
    # 用FBP和SCF测试
    print("\n\n实际计算测试:")
    print("="*60)
    
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=q_drude,
        polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}算法:")
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001 if algo == "SCF" else 0.1
        params.maxIterations = 500 if algo == "SCF" else 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.1
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude位置
        state.atoms[1].x = state.atoms[0].x
        state.atoms[1].y = state.atoms[0].y
        state.atoms[1].z = state.atoms[0].z
        
        # 计算
        energy = force.calculateEnergySCF(state)
        
        # 获取位移
        actual_disp = state.atoms[1].x - state.atoms[0].x
        
        print(f"  计算位移: {actual_disp:.6f} nm = {actual_disp*1000:.3f} pm")
        print(f"  理论位移: {displacement:.6f} nm = {displacement*1000:.3f} pm")
        print(f"  差异: {abs(actual_disp - displacement)*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")
        
        # 验证力平衡
        F_spring_actual = -k_spring * actual_disp
        r_actual = r - actual_disp
        F_electric_actual = k_elec * q_drude * q_external / r_actual**2
        
        print(f"  实际弹簧力: {F_spring_actual:.3f} kJ/(mol·nm)")
        print(f"  实际电力: {F_electric_actual:.3f} kJ/(mol·nm)")
        print(f"  实际总力: {F_spring_actual + F_electric_actual:.3f} kJ/(mol·nm)")

if __name__ == "__main__":
    test_force_balance_logic()