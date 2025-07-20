#!/usr/bin/env python3
"""
测试FBP算法的单位一致性问题
"""

import pygcmc
import numpy as np

def test_single_drude_in_field():
    """测试单个Drude在外电场中的行为"""
    print("测试：单个Drude在均匀外电场中")
    print("="*60)
    
    # 创建一个简单系统：一个原子和一个Drude
    atoms = []
    
    # Parent原子（不带电）
    atom = pygcmc.MCAtom()
    atom.x = 0.0
    atom.y = 0.0
    atom.z = 0.0
    atom.charge = 0.0  # 母原子不带电
    atom.type = 0
    atoms.append(atom)
    
    # Drude粒子
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0  # Drude带-1电荷
    drude.type = 1
    atoms.append(drude)
    
    # 外部电荷（产生电场）
    external = pygcmc.MCAtom()
    external.x = 1.0  # 1 nm away
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.0  # +1电荷
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
    
    # 理论计算
    # 电场 E = k*q/r^2 = 138.935 * 1.0 / 1.0^2 = 138.935 kJ/(mol·nm·e)
    # 力 F = q*E = -1.0 * 138.935 = -138.935 kJ/(mol·nm)
    # 极化率 α 和力常数 k 的关系：k = q^2/(4πε0*α)
    # 设置一个简单的极化率
    polarizability = 0.001  # nm^3
    k = 1.0 * 1.0 * 138.935 / polarizability  # = 138935 kJ/(mol·nm^2)
    
    print(f"理论计算：")
    print(f"  电场 E = 138.935 kJ/(mol·nm·e)")
    print(f"  力 F = q*E = -138.935 kJ/(mol·nm)")
    print(f"  力常数 k = {k:.1f} kJ/(mol·nm^2)")
    print(f"  平衡位移 r = F/k = {-138.935/k:.6f} nm = {-138935/k:.3f} pm")
    
    # 创建Drude力
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
    
    # 测试不同算法
    algorithms = ["SCF", "FBP"]
    
    for algo in algorithms:
        print(f"\n{algo}算法:")
        
        # 设置参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.001 if algo == "SCF" else 0.1
        params.maxIterations = 500 if algo == "SCF" else 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.1  # 较大的限制
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
        disp_x = state.atoms[1].x - state.atoms[0].x
        disp_y = state.atoms[1].y - state.atoms[0].y
        disp_z = state.atoms[1].z - state.atoms[0].z
        disp = np.sqrt(disp_x**2 + disp_y**2 + disp_z**2)
        
        print(f"  能量: {energy:.6f} kJ/mol")
        print(f"  位移: ({disp_x:.6f}, {disp_y:.6f}, {disp_z:.6f}) nm")
        print(f"  位移大小: {disp*1000:.3f} pm")
        
        # 计算能量组分
        # 谐振子能量
        E_harmonic = 0.5 * k * disp**2
        # 静电能量（Drude与外电荷）
        r_drude_external = np.sqrt((state.atoms[1].x - 1.0)**2 + 
                                  state.atoms[1].y**2 + 
                                  state.atoms[1].z**2)
        E_coulomb = 138.935 * (-1.0) * 1.0 / r_drude_external
        
        print(f"  谐振子能量: {E_harmonic:.6f} kJ/mol")
        print(f"  静电能量: {E_coulomb:.6f} kJ/mol")
        print(f"  计算的总能量: {E_harmonic + E_coulomb:.6f} kJ/mol")
        
        # 验证力平衡
        F_spring = -k * disp_x  # 弹簧力（只考虑x方向）
        F_electric = -138.935 / r_drude_external**2 * (state.atoms[1].x - 1.0) / r_drude_external
        print(f"  弹簧力(x): {F_spring:.3f} kJ/(mol·nm)")
        print(f"  电场力(x): {F_electric:.3f} kJ/(mol·nm)")
        print(f"  净力(x): {F_spring + F_electric:.3f} kJ/(mol·nm)")

def test_force_calculation_consistency():
    """测试力计算的一致性"""
    print("\n\n测试：能量-力一致性")
    print("="*60)
    
    # 创建两个水分子
    atoms = []
    residues = []
    
    for i in range(2):
        base_x = i * 0.5
        positions = [
            (base_x, 0.0, 0.0, 1.71636, 0),   # O
            (base_x, 0.0, 0.0, -1.71636, 1),  # D
            (base_x + 0.09572, 0.0, 0.0, 0.55733, 2),  # H1
            (base_x - 0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
            (base_x, -0.024034, 0.0, -1.11466, 3)  # M
        ]
        
        for x, y, z, charge, typ in positions:
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charge
            atom.type = typ
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = 2
    
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建force
    force = pygcmc.DrudeForce()
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    force.addScreenedPair(0, 1, 1.3)
    
    # 数值计算力（通过能量的有限差分）
    print("通过有限差分计算力：")
    
    # 设置一个小的Drude位移
    delta = 0.0001  # nm
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1000.0  # 不优化
    params.maxIterations = 1
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 计算不同位移下的能量
    print(f"{'Displacement (pm)':<20} {'Energy (kJ/mol)':<20} {'Force (kJ/mol/nm)':<20}")
    print("-"*60)
    
    displacements = [-2*delta, -delta, 0, delta, 2*delta]
    energies = []
    
    for disp in displacements:
        # 设置Drude位置
        state.atoms[1].x = state.atoms[0].x + disp
        state.atoms[1].y = state.atoms[0].y
        state.atoms[1].z = state.atoms[0].z
        
        state.atoms[6].x = state.atoms[5].x + disp
        state.atoms[6].y = state.atoms[5].y
        state.atoms[6].z = state.atoms[5].z
        
        energy = force.calculateEnergySCF(state)
        energies.append(energy)
        
        # 计算力（中心差分）
        if len(energies) >= 3:
            idx = len(energies) - 2
            if idx > 0 and idx < len(displacements) - 1:
                f = -(energies[idx+1] - energies[idx-1]) / (2 * delta)
                print(f"{disp*1000:<20.1f} {energy:<20.6f} {f:<20.3f}")
        else:
            print(f"{disp*1000:<20.1f} {energy:<20.6f}")

def main():
    """主函数"""
    print("FBP单位一致性测试")
    print("="*80)
    
    test_single_drude_in_field()
    test_force_calculation_consistency()
    
    print("\n\n分析：")
    print("="*60)
    print("1. 检查电场到力的转换是否正确")
    print("2. 检查力到位移的转换是否正确")
    print("3. 验证力平衡时能量是否最小")

if __name__ == "__main__":
    main()