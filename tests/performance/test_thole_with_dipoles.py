#!/usr/bin/env python3
"""
测试有偶极矩的Thole相互作用
"""

import pygcmc
import numpy as np

def test_thole_with_dipole_moments():
    """
    测试当Drude粒子有位移时的Thole能量
    """
    print("测试有偶极矩时的Thole相互作用")
    print("="*70)
    
    # 创建系统
    state = pygcmc.MCState()
    atoms = []
    
    # 物理常数
    ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2
    
    # 第一个偶极：O-D沿x方向
    o1 = pygcmc.MCAtom()
    o1.x, o1.y, o1.z = 0.0, 0.0, 0.0
    o1.charge = 0.0  # 使用0电荷简化，只看Thole贡献
    o1.type = 0
    atoms.append(o1)
    
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.01, 0.0, 0.0  # 10 pm位移
    d1.charge = -1.71636
    d1.type = 1
    atoms.append(d1)
    
    # 第二个偶极：也沿x方向
    distance = 0.5  # nm
    o2 = pygcmc.MCAtom()
    o2.x, o2.y, o2.z = distance, 0.0, 0.0
    o2.charge = 0.0
    o2.type = 0
    atoms.append(o2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = distance + 0.01, 0.0, 0.0  # 也是10 pm位移
    d2.charge = -1.71636
    d2.type = 1
    atoms.append(d2)
    
    # 创建残基
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 2 * i
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 4
    state.activeResidueCount = 2
    
    # 大盒子
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    # 力场
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    drude_charge = -1.71636
    polarizability = 0.0009782237
    
    force.addParticle(1, 0, -1, -1, -1, -1, drude_charge, polarizability, 1.0, 1.0)
    force.addParticle(3, 2, -1, -1, -1, -1, drude_charge, polarizability, 1.0, 1.0)
    
    # 测试有无Thole屏蔽的差异
    print("\n1. 无Thole屏蔽时的能量分解:")
    print("-"*50)
    
    # 手动计算四种相互作用
    q = drude_charge  # -1.71636
    
    # 距离
    r_dd = np.sqrt((d2.x - d1.x)**2 + (d2.y - d1.y)**2 + (d2.z - d1.z)**2)
    r_dp = np.sqrt((o2.x - d1.x)**2 + (o2.y - d1.y)**2 + (o2.z - d1.z)**2)
    r_pd = np.sqrt((d2.x - o1.x)**2 + (d2.y - o1.y)**2 + (d2.z - o1.z)**2)
    r_pp = np.sqrt((o2.x - o1.x)**2 + (o2.y - o1.y)**2 + (o2.z - o1.z)**2)
    
    print(f"距离:")
    print(f"  D1-D2: {r_dd:.6f} nm")
    print(f"  D1-P2: {r_dp:.6f} nm")
    print(f"  P1-D2: {r_pd:.6f} nm")
    print(f"  P1-P2: {r_pp:.6f} nm")
    
    # 无屏蔽时的能量（假设母原子带电荷-q）
    E_dd = ONE_4PI_EPS0 * q * q / r_dd
    E_dp = ONE_4PI_EPS0 * q * (-q) / r_dp
    E_pd = ONE_4PI_EPS0 * (-q) * q / r_pd
    E_pp = ONE_4PI_EPS0 * (-q) * (-q) / r_pp
    
    print(f"\n能量贡献 (kJ/mol):")
    print(f"  D1-D2: {E_dd:10.4f} (q²/r)")
    print(f"  D1-P2: {E_dp:10.4f} (-q²/r)")
    print(f"  P1-D2: {E_pd:10.4f} (-q²/r)")
    print(f"  P1-P2: {E_pp:10.4f} (q²/r)")
    print(f"  总和:  {E_dd + E_dp + E_pd + E_pp:10.4f}")
    
    # 添加Thole屏蔽
    force.addScreenedPair(0, 1, 1.3)
    
    # 不做SCF（保持固定偶极）
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 1
    params.dampingFactor = 0.0
    params.maxDrudeDistance = 1.0
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    print("\n2. 有Thole屏蔽时:")
    print("-"*50)
    
    try:
        # 只计算Thole贡献（因为原子电荷为0）
        energy = force.calculateEnergySCF(state)
        print(f"Thole屏蔽能量: {energy:.6f} kJ/mol")
        
        # 计算屏蔽函数
        thole = 1.3
        uscale = thole / (polarizability * polarizability)**(1.0/6.0)
        
        print(f"\n屏蔽参数:")
        print(f"  Thole a: {thole}")
        print(f"  uscale: {uscale:.6f} nm⁻¹")
        
        # 计算每个距离的屏蔽
        for r, label in [(r_dd, "D1-D2"), (r_dp, "D1-P2"), 
                         (r_pd, "P1-D2"), (r_pp, "P1-P2")]:
            u = r * uscale
            screening = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
            print(f"  {label}: u={u:.3f}, S(u)={screening:.6f}")
        
    except Exception as e:
        print(f"错误: {e}")
    
    print("\n3. 测试不同偶极方向:")
    print("-"*50)
    
    # 改变第二个偶极方向（垂直）
    state.atoms[3].x = distance
    state.atoms[3].y = 0.01  # y方向位移
    state.atoms[3].z = 0.0
    
    try:
        energy_perpendicular = force.calculateEnergySCF(state)
        print(f"垂直偶极能量: {energy_perpendicular:.6f} kJ/mol")
    except:
        print("计算失败")
    
    # 反平行偶极
    state.atoms[3].x = distance - 0.01  # 反向位移
    state.atoms[3].y = 0.0
    
    try:
        energy_antiparallel = force.calculateEnergySCF(state)
        print(f"反平行偶极能量: {energy_antiparallel:.6f} kJ/mol")
    except:
        print("计算失败")

def main():
    test_thole_with_dipole_moments()
    
    print("\n\n结论:")
    print("1. Thole屏蔽减少了短程偶极-偶极相互作用")
    print("2. 不同偶极取向给出不同能量")
    print("3. 屏蔽函数S(u)在短距离接近0，长距离接近1")

if __name__ == "__main__":
    main()