#!/usr/bin/env python3
"""
调试为什么只有O-D对时Drude不移动
"""

import numpy as np
import pygcmc

def test_gradual_system_building():
    """
    逐步构建系统，看哪一步开始有Drude位移
    """
    print("逐步构建系统测试")
    print("="*70)
    
    # 基础设置
    box_size = 2.0
    
    # 测试1：最简单的O-D对
    print("\n测试1：只有O-D对")
    print("-"*60)
    
    state1 = pygcmc.MCState()
    state1.info.box = np.array([box_size, box_size, box_size])
    state1.info.cutoff = 0.9
    
    atoms1 = []
    
    # O原子
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 1.0, 1.0, 1.0
    o.charge = 1.71636
    o.type = 0
    atoms1.append(o)
    
    # D原子
    d = pygcmc.MCAtom()
    d.x, d.y, d.z = 1.0, 1.0, 1.0
    d.charge = -1.71636
    d.type = 1
    atoms1.append(d)
    
    # 外部正电荷
    ext = pygcmc.MCAtom()
    ext.x, ext.y, ext.z = 1.5, 1.0, 1.0
    ext.charge = 5.0
    ext.type = 0
    atoms1.append(ext)
    
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 3
    res1.active = True
    res1.type = 0
    
    state1.atoms = atoms1
    state1.residues = [res1]
    state1.activeAtomCount = 3
    state1.activeResidueCount = 1
    
    state1.forcefield.numTotalTypes = 2
    state1.forcefield.numMovementTypes = 2
    state1.forcefield.ljSigma = [0.0] * 4
    state1.forcefield.ljEps = [0.0] * 4
    
    test_scf(state1, "O-D对 + 外部电荷")
    
    # 测试2：O-D对作为独立残基
    print("\n\n测试2：O-D在残基1，外部电荷在残基2")
    print("-"*60)
    
    state2 = state1.copy()
    
    # 修改残基定义
    res2_1 = pygcmc.MCResidue()
    res2_1.atomStart = 0
    res2_1.atomCount = 2  # 只包含O和D
    res2_1.active = True
    res2_1.type = 0
    
    res2_2 = pygcmc.MCResidue()
    res2_2.atomStart = 2
    res2_2.atomCount = 1  # 外部电荷
    res2_2.active = True
    res2_2.type = 1
    
    state2.residues = [res2_1, res2_2]
    state2.activeResidueCount = 2
    
    test_scf(state2, "分离的残基")
    
    # 测试3：添加一个氢原子
    print("\n\n测试3：O-D-H系统")
    print("-"*60)
    
    state3 = pygcmc.MCState()
    state3.info.box = np.array([box_size, box_size, box_size])
    state3.info.cutoff = 0.9
    
    atoms3 = []
    
    # O原子
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 1.0, 1.0, 1.0
    o.charge = 1.71636
    o.type = 0
    atoms3.append(o)
    
    # D原子
    d = pygcmc.MCAtom()
    d.x, d.y, d.z = 1.0, 1.0, 1.0
    d.charge = -1.71636
    d.type = 1
    atoms3.append(d)
    
    # H原子
    h = pygcmc.MCAtom()
    h.x, h.y, h.z = 1.096, 1.0, 1.0
    h.charge = 0.55733
    h.type = 2
    atoms3.append(h)
    
    # 外部电荷
    ext = pygcmc.MCAtom()
    ext.x, ext.y, ext.z = 1.5, 1.0, 1.0
    ext.charge = 5.0
    ext.type = 0
    atoms3.append(ext)
    
    res3_1 = pygcmc.MCResidue()
    res3_1.atomStart = 0
    res3_1.atomCount = 3
    res3_1.active = True
    res3_1.type = 0
    
    res3_2 = pygcmc.MCResidue()
    res3_2.atomStart = 3
    res3_2.atomCount = 1
    res3_2.active = True
    res3_2.type = 1
    
    state3.atoms = atoms3
    state3.residues = [res3_1, res3_2]
    state3.activeAtomCount = 4
    state3.activeResidueCount = 2
    
    state3.forcefield.numTotalTypes = 3
    state3.forcefield.numMovementTypes = 3
    state3.forcefield.ljSigma = [0.0] * 9
    state3.forcefield.ljEps = [0.0] * 9
    
    test_scf(state3, "O-D-H系统")
    
    # 测试4：完整水分子
    print("\n\n测试4：完整水分子（O-D-H-H-M）")
    print("-"*60)
    
    state4 = pygcmc.MCState()
    state4.info.box = np.array([box_size, box_size, box_size])
    state4.info.cutoff = 0.9
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    positions = [
        [1.0, 1.0, 1.0],      # O
        [1.0, 1.0, 1.0],      # D
        [1.096, 1.0, 1.0],    # H1
        [0.952, 1.077, 1.0],  # H2
        [1.0, 1.0, 1.0]       # M
    ]
    
    atoms4 = []
    
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = i if i < 4 else 3
        atoms4.append(atom)
    
    # 外部电荷
    ext = pygcmc.MCAtom()
    ext.x, ext.y, ext.z = 1.5, 1.0, 1.0
    ext.charge = 5.0
    ext.type = 0
    atoms4.append(ext)
    
    res4_1 = pygcmc.MCResidue()
    res4_1.atomStart = 0
    res4_1.atomCount = 5
    res4_1.active = True
    res4_1.type = 0
    
    res4_2 = pygcmc.MCResidue()
    res4_2.atomStart = 5
    res4_2.atomCount = 1
    res4_2.active = True
    res4_2.type = 1
    
    state4.atoms = atoms4
    state4.residues = [res4_1, res4_2]
    state4.activeAtomCount = 6
    state4.activeResidueCount = 2
    
    state4.forcefield.numTotalTypes = 5
    state4.forcefield.numMovementTypes = 5
    ljSigma = []
    ljEps = []
    for i in range(5):
        for j in range(5):
            if i == 0 and j == 0:
                ljSigma.append(0.318395)
                ljEps.append(0.88257)
            else:
                ljSigma.append(0.0)
                ljEps.append(0.0)
    
    state4.forcefield.ljSigma = ljSigma
    state4.forcefield.ljEps = ljEps
    
    test_scf(state4, "完整水分子")

def test_scf(state, description):
    """
    测试SCF优化
    """
    print(f"{description}:")
    
    # 计算外部电场（理论值）
    ONE_4PI_EPS0 = 138.935456
    r = 0.5  # nm
    q_ext = 5.0
    E_field = ONE_4PI_EPS0 * q_ext / (r * r)
    
    # 预期位移
    q_drude = -1.71636
    alpha = 0.0009782237
    k_spring = q_drude**2 * ONE_4PI_EPS0 / (2 * alpha)
    expected_disp = abs(q_drude) * E_field / k_spring * 1000  # pm
    
    print(f"  理论电场: {E_field:.1f} kJ/(mol·nm·e)")
    print(f"  预期位移: {expected_disp:.2f} pm")
    
    # 创建DrudeForce
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
    params.tolerance = 0.1
    params.maxIterations = 500
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    state_opt = state.copy()
    
    try:
        energy = force.calculateEnergySCF(state_opt)
        
        dx = state_opt.atoms[1].x - state_opt.atoms[0].x
        dy = state_opt.atoms[1].y - state_opt.atoms[0].y
        dz = state_opt.atoms[1].z - state_opt.atoms[0].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        
        print(f"  实际位移: {disp:.2f} pm")
        print(f"  能量: {energy:.4f} kJ/mol")
        
        if disp > 1.0:
            print(f"  ✓ Drude响应了电场")
        else:
            print(f"  ✗ Drude没有响应")
            
    except Exception as e:
        print(f"  失败: {e}")

def main():
    """
    主函数
    """
    test_gradual_system_building()
    
    print("\n\n结论：")
    print("="*70)
    print("PyGCMC的DrudeForce可能需要：")
    print("1. 完整的分子结构")
    print("2. 正确的残基定义")
    print("3. 分子内的其他原子来建立正确的电场环境")
    print("\n这可能是为了确保分子内相互作用被正确处理")

if __name__ == "__main__":
    main()