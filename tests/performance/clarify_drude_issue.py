#!/usr/bin/env python3
"""
理清PyGCMC Drude实现的问题
"""

import numpy as np
import pygcmc

def test_what_pygcmc_drude_does():
    """
    搞清楚PyGCMC的Drude到底在做什么
    """
    print("理清PyGCMC Drude实现")
    print("="*70)
    
    # 1. 我们之前的测试显示了什么？
    print("\n1. 之前测试的关键结果：")
    print("-"*60)
    print("a) test_scf_energy_change.py:")
    print("   - 初始能量（Drude在parent）: -97.42 kJ/mol/水")
    print("   - SCF优化后: -165.14 kJ/mol/水")
    print("   - 能量降低: 67.72 kJ/mol/水")
    print("   - Drude位移: 平均11.24 pm, 最大19.39 pm")
    print("   => 说明SCF确实在工作！")
    
    print("\nb) compare_openmm_pygcmc_scf.py:")
    print("   - OpenMM: -1.88 kJ/mol/水, 位移0.48 pm")
    print("   - PyGCMC: -26.31 kJ/mol/水, 位移0.63 pm")
    print("   => 位移相似，但能量差异巨大")
    
    print("\nc) 刚才的测试:")
    print("   - 两个水分子系统，PyGCMC显示0位移")
    print("   - computeSystemEnergyDrude只返回Drude内部能量")
    print("   => 这与之前的结果矛盾！")
    
    # 2. 重新测试256水系统
    print("\n\n2. 重新测试之前成功的256水系统")
    print("-"*60)
    
    import pickle
    import os
    
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    if not os.path.exists(filename):
        print("文件不存在")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    # 使用前10个水分子
    n_waters = 10
    positions = data['positions']
    box_length = data['box_length']
    charges = data['charges']
    
    print(f"系统: {n_waters}个水分子, 盒子{box_length:.3f} nm")
    
    # 创建状态
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
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
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
    
    # 添加一些Thole对
    n_pairs = 0
    for i in range(min(5, n_waters)):
        for j in range(i+1, min(i+5, n_waters)):
            force.addScreenedPair(i, j, 1.3)
            n_pairs += 1
    
    print(f"添加了{n_pairs}个Thole对")
    
    # SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 确保Drude在parent位置
    state_test = state.copy()
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state_test.atoms[d_idx].x = state_test.atoms[o_idx].x
        state_test.atoms[d_idx].y = state_test.atoms[o_idx].y
        state_test.atoms[d_idx].z = state_test.atoms[o_idx].z
    
    # 运行SCF
    print("\n运行SCF...")
    try:
        energy = force.calculateEnergySCF(state_test)
        print(f"SCF能量: {energy:.2f} kJ/mol ({energy/n_waters:.2f} kJ/mol/水)")
        
        # 检查位移
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
        max_disp = np.max(displacements)
        
        print(f"平均Drude位移: {avg_disp:.2f} pm")
        print(f"最大Drude位移: {max_disp:.2f} pm")
        
        if avg_disp > 1.0:
            print("\n✓ Drude确实在移动！")
        else:
            print("\n✗ Drude几乎没有移动")
            
    except Exception as e:
        print(f"SCF失败: {e}")
    
    # 3. 分析差异
    print("\n\n3. 为什么简单系统和复杂系统表现不同？")
    print("-"*60)
    print("可能的原因：")
    print("a) 密度差异：")
    print("   - 简单2水系统：密度很低")
    print("   - 256水系统：密度1.0 g/cm³")
    print("   => 高密度下水分子更近，电场更强")
    
    print("\nb) Thole对数量：")
    print("   - 简单系统：1个Thole对")
    print("   - 复杂系统：很多Thole对")
    print("   => Thole屏蔽可能影响电场传递")
    
    print("\nc) 系统配置：")
    print("   - 简单系统：人工摆放")
    print("   - 256水系统：经过NVT优化")
    print("   => 优化的结构可能有更强的相互作用")

def test_drude_with_external_field():
    """
    测试外加电场下的Drude响应
    """
    print("\n\n4. 测试外加电场")
    print("="*70)
    
    # 创建单个水分子
    state = pygcmc.MCState()
    
    box_size = 1.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    positions = [
        [0.5, 0.5, 0.5],      # O
        [0.5, 0.5, 0.5],      # D
        [0.596, 0.5, 0.5],    # H1
        [0.452, 0.577, 0.5],  # H2
        [0.5, 0.5, 0.5]       # M
    ]
    
    atoms = []
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i]
        atom.type = i if i < 4 else 3
        atoms.append(atom)
    
    # 添加一个外部点电荷来产生电场
    external_charge = pygcmc.MCAtom()
    external_charge.x = 0.7
    external_charge.y = 0.5
    external_charge.z = 0.5
    external_charge.charge = 10.0  # 大电荷产生强电场
    external_charge.type = 0
    atoms.append(external_charge)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 6  # 包括外部电荷
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 6
    state.activeResidueCount = 1
    
    # 力场
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
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
    params.tolerance = 100.0
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    print("单个水分子 + 外部点电荷(q=+10)")
    
    try:
        energy = force.calculateEnergySCF(state)
        print(f"能量: {energy:.2f} kJ/mol")
        
        dx = state.atoms[1].x - state.atoms[0].x
        dy = state.atoms[1].y - state.atoms[0].y
        dz = state.atoms[1].z - state.atoms[0].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        
        print(f"Drude位移: {disp:.2f} pm")
        print(f"位移方向: ({dx*1000:.2f}, {dy*1000:.2f}, {dz*1000:.2f}) pm")
        
        if disp > 1.0:
            print("\n✓ Drude响应外部电场！")
        else:
            print("\n✗ Drude不响应外部电场")
            
    except Exception as e:
        print(f"失败: {e}")

def main():
    """
    主函数
    """
    test_what_pygcmc_drude_does()
    test_drude_with_external_field()
    
    print("\n\n最终结论：")
    print("="*70)
    print("需要确认：")
    print("1. PyGCMC的SCF是否考虑了系统中其他原子的电场？")
    print("2. 为什么256水系统有位移，而简单系统没有？")
    print("3. calculateEnergySCF到底包含哪些能量项？")

if __name__ == "__main__":
    main()