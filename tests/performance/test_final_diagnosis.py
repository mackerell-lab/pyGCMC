#!/usr/bin/env python3
"""
最终诊断：理解SCF和FBP的本质差异
"""

import pygcmc
import numpy as np

def test_energy_landscape():
    """测试能量景观"""
    print("能量景观分析")
    print("="*60)
    
    # 创建简单系统
    atoms = []
    
    # Parent (不带电)
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
    
    # 参数
    k_spring = 138935.0
    k_elec = 138.935
    
    # 创建force
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
    
    # 扫描不同的Drude位置
    print("扫描Drude位置的能量景观:")
    print(f"{'位移(pm)':<10} {'谐振子能':<15} {'静电能':<15} {'总能量':<15} {'弹簧力':<15} {'电力':<15} {'净力':<15}")
    print("-"*110)
    
    displacements = np.linspace(-5, 5, 21) / 1000  # -5到+5 pm，单位nm
    
    min_energy = float('inf')
    min_disp = 0
    zero_force_disp = None
    
    for disp in displacements:
        # 设置Drude位置
        state.atoms[1].x = disp
        
        # 计算能量组分
        # 1. 谐振子能量
        E_harmonic = 0.5 * k_spring * disp**2
        
        # 2. 静电能量
        r_DE = abs(disp - 1.0)  # Drude到External的距离
        E_coulomb = k_elec * (-1.0) * 1.0 / r_DE
        
        # 总能量
        E_total = E_harmonic + E_coulomb
        
        # 力分析
        F_spring = -k_spring * disp  # 弹簧力
        # 电场力（注意方向）
        if disp < 1.0:  # Drude在External左边
            F_electric = k_elec * 1.0 / r_DE**2  # 吸引力，向右
        else:  # Drude在External右边
            F_electric = -k_elec * 1.0 / r_DE**2  # 吸引力，向左
        
        F_net = F_spring + F_electric
        
        print(f"{disp*1000:>10.2f} {E_harmonic:>15.6f} {E_coulomb:>15.6f} {E_total:>15.6f} "
              f"{F_spring:>15.1f} {F_electric:>15.1f} {F_net:>15.1f}")
        
        if E_total < min_energy:
            min_energy = E_total
            min_disp = disp
        
        if zero_force_disp is None and abs(F_net) < 1.0:
            zero_force_disp = disp
    
    print("\n分析结果:")
    print(f"  最低能量: {min_energy:.6f} kJ/mol at {min_disp*1000:.3f} pm")
    print(f"  力平衡点: {zero_force_disp*1000:.3f} pm" if zero_force_disp else "  力平衡点: 未找到")
    
    # 测试算法
    print("\n\n算法结果:")
    print("-"*40)
    
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001
        params.maxIterations = 500
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置
        state.atoms[1].x = 0.0
        
        # 计算
        energy = force.calculateEnergySCF(state)
        
        disp = state.atoms[1].x
        print(f"\n{algo}:")
        print(f"  位移: {disp*1000:.3f} pm")
        print(f"  能量: {energy:.6f} kJ/mol")
        
        # 验证
        r_DE = abs(disp - 1.0)
        E_manual = 0.5 * k_spring * disp**2 + k_elec * (-1.0) * 1.0 / r_DE
        print(f"  手动计算能量: {E_manual:.6f} kJ/mol")
        
        # 力
        F_spring = -k_spring * disp
        if disp < 1.0:
            F_electric = k_elec * 1.0 / r_DE**2
        else:
            F_electric = -k_elec * 1.0 / r_DE**2
        print(f"  弹簧力: {F_spring:.1f} kJ/(mol·nm)")
        print(f"  电力: {F_electric:.1f} kJ/(mol·nm)")
        print(f"  净力: {F_spring + F_electric:.1f} kJ/(mol·nm)")
    
    print("\n\n最终诊断:")
    print("="*60)
    print("1. 能量最小点在正位移处（Drude被吸引向右）")
    print("2. FBP找到了力平衡点（正位移）")
    print("3. SCF给出了负位移，这不是能量最小点")
    print("4. 问题：SCF可能在优化不同的目标函数，或有其他约束")

if __name__ == "__main__":
    test_energy_landscape()