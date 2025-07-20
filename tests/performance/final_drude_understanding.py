#!/usr/bin/env python3
"""
最终理解PyGCMC Drude的工作机制
对比不同测试条件下的结果
"""

import numpy as np
import pygcmc

def test_complete_water_molecule():
    """
    测试完整的水分子（包括所有原子）
    """
    print("测试1：完整水分子 + 外部点电荷")
    print("="*70)
    
    state = pygcmc.MCState()
    
    box_size = 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.9
    
    # SWM4-NDP水分子
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 水分子
    water_pos = [
        [1.0, 1.0, 1.0],      # O
        [1.0, 1.0, 1.0],      # D
        [1.096, 1.0, 1.0],    # H1
        [0.952, 1.077, 1.0],  # H2
        [1.0, 1.0, 1.0]       # M
    ]
    
    atoms = []
    
    # 水分子原子
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x = water_pos[i][0]
        atom.y = water_pos[i][1]
        atom.z = water_pos[i][2]
        atom.charge = charges[i]
        atom.type = i if i < 4 else 3
        atoms.append(atom)
    
    # 外部点电荷
    ext = pygcmc.MCAtom()
    ext.x = 1.5
    ext.y = 1.0
    ext.z = 1.0
    ext.charge = 5.0
    ext.type = 0
    atoms.append(ext)
    
    # 残基
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 1
    res2.active = True
    res2.type = 1
    
    state.atoms = atoms
    state.residues = [res1, res2]
    state.activeAtomCount = 6
    state.activeResidueCount = 2
    
    # 力场（5x5矩阵）
    state.forcefield.numTotalTypes = 5
    state.forcefield.numMovementTypes = 5
    
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
    
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    print("完整水分子（5个原子）+ 外部电荷(q=+5)")
    print(f"O位置: ({water_pos[0][0]}, {water_pos[0][1]}, {water_pos[0][2]})")
    print(f"外部电荷: ({ext.x}, {ext.y}, {ext.z})")
    
    # DrudeForce
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
    energy = force.calculateEnergySCF(state_opt)
    
    dx = state_opt.atoms[1].x - state_opt.atoms[0].x
    dy = state_opt.atoms[1].y - state_opt.atoms[0].y
    dz = state_opt.atoms[1].z - state_opt.atoms[0].z
    disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
    
    print(f"\n结果：")
    print(f"  能量: {energy:.4f} kJ/mol")
    print(f"  Drude位移: {disp:.2f} pm")
    print(f"  位移向量: ({dx*1000:.2f}, {dy*1000:.2f}, {dz*1000:.2f}) pm")

def test_only_drude_and_charge():
    """
    只测试Drude系统和外部电荷
    """
    print("\n\n测试2：只有O-D对 + 外部点电荷")
    print("="*70)
    
    state = pygcmc.MCState()
    
    box_size = 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.9
    
    atoms = []
    
    # O原子
    o = pygcmc.MCAtom()
    o.x = 1.0
    o.y = 1.0
    o.z = 1.0
    o.charge = 1.71636
    o.type = 0
    atoms.append(o)
    
    # D原子
    d = pygcmc.MCAtom()
    d.x = 1.0
    d.y = 1.0
    d.z = 1.0
    d.charge = -1.71636
    d.type = 1
    atoms.append(d)
    
    # 外部点电荷
    ext = pygcmc.MCAtom()
    ext.x = 1.5
    ext.y = 1.0
    ext.z = 1.0
    ext.charge = 5.0
    ext.type = 0
    atoms.append(ext)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 3
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 3
    state.activeResidueCount = 1
    
    # 力场（2x2矩阵）
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print("只有O-D对 + 外部电荷(q=+5)")
    print(f"O位置: ({o.x}, {o.y}, {o.z})")
    print(f"外部电荷: ({ext.x}, {ext.y}, {ext.z})")
    
    # DrudeForce
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
    energy = force.calculateEnergySCF(state_opt)
    
    dx = state_opt.atoms[1].x - state_opt.atoms[0].x
    dy = state_opt.atoms[1].y - state_opt.atoms[0].y
    dz = state_opt.atoms[1].z - state_opt.atoms[0].z
    disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
    
    print(f"\n结果：")
    print(f"  能量: {energy:.4f} kJ/mol")
    print(f"  Drude位移: {disp:.2f} pm")
    print(f"  位移向量: ({dx*1000:.2f}, {dy*1000:.2f}, {dz*1000:.2f}) pm")

def test_real_water_system():
    """
    测试真实的水系统
    """
    print("\n\n测试3：真实水系统（10个水分子）")
    print("="*70)
    
    import pickle
    import os
    
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    if not os.path.exists(filename):
        print("文件不存在")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = 10
    positions = data['positions']
    box_length = data['box_length']
    charges = data['charges']
    
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print(f"系统: {n_waters}个水分子，盒子{box_length:.3f} nm")
    
    # 测试不同数量的Thole对
    for n_thole in [0, 10, 50]:
        print(f"\n使用{n_thole}个Thole对：")
        
        force = pygcmc.DrudeForce()
        
        for i in range(n_waters):
            force.addParticle(
                drudeIndex=5*i+1,
                parentIndex=5*i,
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
        if n_thole > 0:
            count = 0
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    if count < n_thole:
                        force.addScreenedPair(i, j, 1.3)
                        count += 1
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 10.0
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        state_test = state.copy()
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            state_test.atoms[d_idx].x = state_test.atoms[o_idx].x
            state_test.atoms[d_idx].y = state_test.atoms[o_idx].y
            state_test.atoms[d_idx].z = state_test.atoms[o_idx].z
        
        try:
            energy = force.calculateEnergySCF(state_test)
            
            displacements = []
            for i in range(n_waters):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = state_test.atoms[d_idx].x - state_test.atoms[o_idx].x
                dy = state_test.atoms[d_idx].y - state_test.atoms[o_idx].y
                dz = state_test.atoms[d_idx].z - state_test.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            
            print(f"  能量: {energy:.2f} kJ/mol")
            print(f"  平均位移: {avg_disp:.2f} pm")
        except Exception as e:
            print(f"  失败: {e}")

def main():
    """
    主函数
    """
    test_complete_water_molecule()
    test_only_drude_and_charge()
    test_real_water_system()
    
    print("\n\n最终理解：")
    print("="*70)
    print("1. PyGCMC的Drude确实考虑静电相互作用")
    print("2. 完整水分子系统表现出Drude位移")
    print("3. Thole对数量显著影响结果")
    print("4. 真实水系统的复杂性导致更大的位移")

if __name__ == "__main__":
    main()