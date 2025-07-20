#!/usr/bin/env python3
"""
测试OpenMM风格的Thole实现
确保电荷模型正确：母原子电荷需要调整
"""

import pygcmc
import numpy as np

def create_test_system_openmm_style():
    """
    创建测试系统，遵循OpenMM的电荷模型
    """
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    # 原始电荷（未加Drude之前）
    original_charges = [1.71636, 0.0, 0.55733, 0.55733, -1.11466]  # O, D位置, H1, H2, M
    drude_charge = -1.71636  # Drude粒子电荷
    
    # 调整后的电荷（OpenMM风格）
    # O原子电荷 = 原始电荷 - Drude电荷 = 1.71636 - (-1.71636) = 3.43272
    adjusted_charges = [
        original_charges[0] - drude_charge,  # O: 3.43272
        drude_charge,                         # D: -1.71636
        original_charges[2],                  # H1: 0.55733
        original_charges[3],                  # H2: 0.55733
        original_charges[4]                   # M: -1.11466
    ]
    
    atom_types = [0, 1, 2, 2, 3]
    
    # 创建两个水分子，距离0.5 nm
    positions = [
        # 第一个水分子
        [0.0, 0.0, 0.0],      # O
        [0.0, 0.0, 0.0],      # D (初始与O重合)
        [0.09572, 0.0, 0.0],  # H1
        [-0.09572, 0.0, 0.0], # H2
        [0.0, 0.024034, 0.0], # M
        # 第二个水分子
        [0.5, 0.0, 0.0],      # O
        [0.5, 0.0, 0.0],      # D (初始与O重合)
        [0.59572, 0.0, 0.0],  # H1
        [0.40428, 0.0, 0.0],  # H2
        [0.5, 0.024034, 0.0]  # M
    ]
    
    # 创建原子
    for i, (pos, charge, atype) in enumerate(zip(positions, 
                                                  adjusted_charges + adjusted_charges, 
                                                  atom_types + atom_types)):
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        atom.charge = charge
        atom.type = atype
        atoms.append(atom)
    
    # 创建残基
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = 2
    
    # 大盒子避免PBC
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, drude_charge

def test_charge_setup():
    """
    测试电荷设置
    """
    print("测试OpenMM风格的电荷设置")
    print("="*60)
    
    state, drude_charge = create_test_system_openmm_style()
    
    print("\n电荷验证:")
    print("-"*40)
    
    for i in range(2):
        print(f"\n水分子 {i+1}:")
        o_charge = state.atoms[i*5].charge
        d_charge = state.atoms[i*5+1].charge
        h1_charge = state.atoms[i*5+2].charge
        h2_charge = state.atoms[i*5+3].charge
        m_charge = state.atoms[i*5+4].charge
        
        print(f"  O原子电荷: {o_charge:.5f} (应该是 3.43272)")
        print(f"  D粒子电荷: {d_charge:.5f} (应该是 -1.71636)")
        print(f"  H1原子电荷: {h1_charge:.5f}")
        print(f"  H2原子电荷: {h2_charge:.5f}")
        print(f"  M位点电荷: {m_charge:.5f}")
        
        total_charge = o_charge + d_charge + h1_charge + h2_charge + m_charge
        print(f"  总电荷: {total_charge:.5f} (应该是 0)")
        
        # 验证偶极
        print(f"  O-D偶极电荷: {o_charge + d_charge:.5f} (应该是原始O电荷 1.71636)")

def test_thole_energy():
    """
    测试Thole能量计算
    """
    print("\n\n测试Thole屏蔽能量")
    print("="*60)
    
    state, drude_charge_value = create_test_system_openmm_style()
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    polarizability = 0.0009782237  # nm³
    
    for i in range(2):
        force.addParticle(
            drudeIndex=5*i+1,
            parentIndex=5*i,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=drude_charge_value,  # 使用定义的Drude电荷
            polarizability=polarizability,
            aniso12=1.0,
            aniso34=1.0
        )
    
    # 添加Thole屏蔽对
    force.addScreenedPair(0, 1, 1.3)
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 计算能量
    print("\n不同距离下的能量:")
    print("-"*40)
    
    distances = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    
    for dist in distances:
        # 移动第二个水分子
        for j in range(5):
            state.atoms[5+j].x = dist + (state.atoms[5+j].x - 0.5)
        
        try:
            energy = force.calculateEnergySCF(state)
            
            # 计算Drude位移
            displacements = []
            for i in range(2):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = state.atoms[d_idx].x - state.atoms[o_idx].x
                dy = state.atoms[d_idx].y - state.atoms[o_idx].y
                dz = state.atoms[d_idx].z - state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
                displacements.append(disp)
            
            print(f"  距离 {dist:.1f} nm: 能量 = {energy:10.4f} kJ/mol, "
                  f"Drude位移 = [{displacements[0]:.1f}, {displacements[1]:.1f}] pm")
            
        except Exception as e:
            print(f"  距离 {dist:.1f} nm: 错误 - {str(e)}")
        
        # 恢复原始位置
        state.atoms[5].x = 0.5
        state.atoms[6].x = 0.5
        state.atoms[7].x = 0.59572
        state.atoms[8].x = 0.40428
        state.atoms[9].x = 0.5

def main():
    """
    主函数
    """
    test_charge_setup()
    test_thole_energy()
    
    print("\n\n关键点:")
    print("1. OpenMM风格要求母原子电荷 = 原始电荷 - Drude电荷")
    print("2. 这确保了偶极矩正确：μ = q_drude * d")
    print("3. Thole屏蔽修正了偶极-偶极相互作用")
    print("4. 如果能量随距离合理变化，说明实现正确")

if __name__ == "__main__":
    main()