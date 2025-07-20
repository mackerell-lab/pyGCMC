#!/usr/bin/env python3
"""
分析PyGCMC Drude的力平衡
理解为什么位移是预期的一半
"""

import numpy as np
import pygcmc

def analyze_drude_forces():
    """
    详细分析Drude粒子的受力情况
    """
    print("分析Drude粒子的力平衡")
    print("="*70)
    
    # 创建简单系统：Drude粒子 + 外部点电荷
    state = pygcmc.MCState()
    
    box_size = 3.0  # 大盒子
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 1.4  # 足够大的截断
    
    atoms = []
    
    # O原子（Drude parent）
    o_atom = pygcmc.MCAtom()
    o_atom.x = 1.5
    o_atom.y = 1.5
    o_atom.z = 1.5
    o_atom.charge = 1.71636  # SWM4-NDP氧电荷
    o_atom.type = 0
    atoms.append(o_atom)
    
    # Drude粒子
    d_atom = pygcmc.MCAtom()
    d_atom.x = 1.5  # 初始在O位置
    d_atom.y = 1.5
    d_atom.z = 1.5
    d_atom.charge = -1.71636  # Drude电荷
    d_atom.type = 1
    atoms.append(d_atom)
    
    # 外部点电荷
    ext_atom = pygcmc.MCAtom()
    ext_atom.x = 2.0  # 距离0.5 nm
    ext_atom.y = 1.5
    ext_atom.z = 1.5
    ext_atom.charge = 1.0  # 单位正电荷
    ext_atom.type = 0
    atoms.append(ext_atom)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 3
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 3
    state.activeResidueCount = 1
    
    # 力场
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.0] * 4
    state.forcefield.ljEps = [0.0] * 4
    
    # 理论计算
    ONE_4PI_EPS0 = 138.935456
    q_drude = -1.71636
    q_ext = 1.0
    r_ext = 0.5  # nm
    alpha = 0.0009782237  # nm³
    
    print("系统参数：")
    print(f"  Drude电荷: {q_drude} e")
    print(f"  极化率: {alpha} nm³")
    print(f"  外部电荷: {q_ext} e，距离 {r_ext} nm")
    
    # 1. 外部电场
    E_ext = ONE_4PI_EPS0 * q_ext / (r_ext * r_ext)
    print(f"\n外部电场：")
    print(f"  E = {E_ext:.1f} kJ/(mol·nm·e)")
    
    # 2. Drude受到的电场力
    F_elec = q_drude * E_ext
    print(f"\nDrude受到的电场力：")
    print(f"  F_elec = q_drude * E = {F_elec:.1f} kJ/(mol·nm)")
    
    # 3. 弹簧常数
    k_spring = q_drude * q_drude * ONE_4PI_EPS0 / (2 * alpha)
    print(f"\n弹簧常数：")
    print(f"  k = q²/(2α) = {k_spring:.1f} kJ/(mol·nm²)")
    
    # 4. 平衡位置（F_spring + F_elec = 0）
    # F_spring = -k * x
    # F_elec = q * E
    # 平衡：-k * x + q * E = 0
    # x = q * E / k
    x_eq = q_drude * E_ext / k_spring
    print(f"\n平衡位置：")
    print(f"  x = F_elec / k = {x_eq*1000:.2f} pm")
    
    # 但这里有个问题：Drude和parent之间也有库仑相互作用！
    print("\n\n考虑Drude-parent库仑相互作用：")
    print("当Drude位移x时，它与parent的库仑相互作用能为：")
    print("  U_self = k_e * q_drude * q_parent / x")
    print("  但当x很小时，这个能量会发散！")
    
    # 5. 运行PyGCMC SCF
    print("\n\n运行PyGCMC SCF：")
    
    force = pygcmc.DrudeForce()
    force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=-1.71636,
        polarizability=0.0009782237,
        aniso12=1.0,
        aniso34=1.0
    )
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 1000
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    state_opt = state.copy()
    
    try:
        energy = force.calculateEnergySCF(state_opt)
        
        dx = state_opt.atoms[1].x - state_opt.atoms[0].x
        disp = dx * 1000  # pm
        
        print(f"  SCF位移: {disp:.2f} pm")
        print(f"  预期位移: {x_eq*1000:.2f} pm")
        print(f"  比率: {disp/(x_eq*1000):.2f}")
        
        # 分析可能的原因
        print("\n可能的原因：")
        print("1. Drude-parent之间的库仑排斥被排除了")
        print("2. 实际的有效弹簧常数可能不同")
        print("3. SCF算法的实现细节")
        
        # 测试不同的极化率
        print("\n\n测试双倍极化率：")
        force2 = pygcmc.DrudeForce()
        force2.addParticle(
            drudeIndex=1,
            parentIndex=0,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.0009782237 * 2,  # 双倍极化率
            aniso12=1.0,
            aniso34=1.0
        )
        force2.setSCFParameters(params)
        force2.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        state_opt2 = state.copy()
        energy2 = force2.calculateEnergySCF(state_opt2)
        
        dx2 = state_opt2.atoms[1].x - state_opt2.atoms[0].x
        disp2 = dx2 * 1000
        
        print(f"  双倍极化率位移: {disp2:.2f} pm")
        print(f"  位移比: {disp2/disp:.2f} (应该是2.0)")
        
    except Exception as e:
        print(f"  失败: {e}")

def test_thole_mechanism():
    """
    测试Thole机制如何产生有效电场
    """
    print("\n\n测试Thole机制")
    print("="*70)
    
    # 创建两个Drude粒子系统
    state = pygcmc.MCState()
    
    box_size = 3.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 1.4
    
    atoms = []
    
    # 第一个Drude系统
    # O1
    o1 = pygcmc.MCAtom()
    o1.x = 1.0
    o1.y = 1.5
    o1.z = 1.5
    o1.charge = 1.71636
    o1.type = 0
    atoms.append(o1)
    
    # D1
    d1 = pygcmc.MCAtom()
    d1.x = 1.0
    d1.y = 1.5
    d1.z = 1.5
    d1.charge = -1.71636
    d1.type = 1
    atoms.append(d1)
    
    # 第二个Drude系统
    # O2
    o2 = pygcmc.MCAtom()
    o2.x = 1.5  # 距离0.5 nm
    o2.y = 1.5
    o2.z = 1.5
    o2.charge = 1.71636
    o2.type = 0
    atoms.append(o2)
    
    # D2
    d2 = pygcmc.MCAtom()
    d2.x = 1.5
    d2.y = 1.5
    d2.z = 1.5
    d2.charge = -1.71636
    d2.type = 1
    atoms.append(d2)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 4
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 4
    state.activeResidueCount = 1
    
    # 力场
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.0] * 4
    state.forcefield.ljEps = [0.0] * 4
    
    print("两个Drude系统，相距0.5 nm")
    
    # 1. 无Thole屏蔽
    print("\n1. 无Thole屏蔽：")
    force1 = pygcmc.DrudeForce()
    
    for i in range(2):
        force1.addParticle(
            drudeIndex=2*i+1,
            parentIndex=2*i,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.0009782237,
            aniso12=1.0,
            aniso34=1.0
        )
    
    # 不添加Thole对
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 1000
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1
    force1.setSCFParameters(params)
    force1.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    state_opt1 = state.copy()
    energy1 = force1.calculateEnergySCF(state_opt1)
    
    dx1_1 = state_opt1.atoms[1].x - state_opt1.atoms[0].x
    dx2_1 = state_opt1.atoms[3].x - state_opt1.atoms[2].x
    
    print(f"  Drude1位移: {dx1_1*1000:.2f} pm")
    print(f"  Drude2位移: {dx2_1*1000:.2f} pm")
    
    # 2. 有Thole屏蔽
    print("\n2. 有Thole屏蔽：")
    force2 = pygcmc.DrudeForce()
    
    for i in range(2):
        force2.addParticle(
            drudeIndex=2*i+1,
            parentIndex=2*i,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.0009782237,
            aniso12=1.0,
            aniso34=1.0
        )
    
    # 添加Thole对
    force2.addScreenedPair(0, 1, 1.3)
    
    force2.setSCFParameters(params)
    force2.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    state_opt2 = state.copy()
    energy2 = force2.calculateEnergySCF(state_opt2)
    
    dx1_2 = state_opt2.atoms[1].x - state_opt2.atoms[0].x
    dx2_2 = state_opt2.atoms[3].x - state_opt2.atoms[2].x
    
    print(f"  Drude1位移: {dx1_2*1000:.2f} pm")
    print(f"  Drude2位移: {dx2_2*1000:.2f} pm")
    print(f"\n  能量差异: {energy2 - energy1:.4f} kJ/mol")
    
    print("\n结论：Thole屏蔽改变了Drude之间的相互作用")

def main():
    """
    主函数
    """
    analyze_drude_forces()
    test_thole_mechanism()
    
    print("\n\n最终结论：")
    print("="*70)
    print("1. PyGCMC的Drude确实考虑了静电相互作用")
    print("2. 实际位移约为理论值的50%，可能因为：")
    print("   - Drude-parent库仑相互作用的处理方式")
    print("   - 有效弹簧常数的定义")
    print("3. Thole屏蔽提供了额外的相互作用机制")

if __name__ == "__main__":
    main()