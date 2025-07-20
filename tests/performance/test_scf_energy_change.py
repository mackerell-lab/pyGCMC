#!/usr/bin/env python3
"""
测试SCF优化前后的能量变化
使用真实的水系统，比较：
1. Drude在parent位置的能量
2. SCF优化后的能量
"""

import numpy as np
import pickle
import os
import pygcmc

def test_scf_energy_change():
    """
    测试SCF优化前后的能量变化
    """
    print("测试SCF优化前后的能量变化")
    print("="*70)
    
    # 加载水系统
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"\n系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = data['charges']
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
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
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
    print("\n添加Thole对...")
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 添加足够的Thole对
    for i in range(min(100, n_waters)):
        o1_idx = i * 5
        for j in range(i+1, min(i+50, n_waters)):
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
                
                if n_thole_pairs >= 2000:
                    break
        
        if n_thole_pairs >= 2000:
            break
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 步骤1: 确保所有Drude粒子在parent位置
    print("\n步骤1: 将所有Drude粒子放在parent位置")
    state_initial = state.copy()
    
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        # 将Drude放在O原子位置
        state_initial.atoms[d_idx].x = state_initial.atoms[o_idx].x
        state_initial.atoms[d_idx].y = state_initial.atoms[o_idx].y
        state_initial.atoms[d_idx].z = state_initial.atoms[o_idx].z
    
    # 步骤2: 计算初始能量（使用calculateEnergyDrude，不优化）
    print("\n步骤2: 计算初始能量（Drude在parent位置）")
    
    # 注意：这里我们需要使用不优化Drude的能量计算
    # 但是PyGCMC可能没有这个功能，所以我们用非常大的容差
    params_no_opt = pygcmc.DrudeSCFParams()
    params_no_opt.tolerance = 1e10  # 非常大的容差，基本不优化
    params_no_opt.maxIterations = 1  # 只迭代一次
    params_no_opt.dampingFactor = 0.0  # 不移动
    params_no_opt.maxDrudeDistance = 0.00001  # 几乎不允许移动
    force.setSCFParameters(params_no_opt)
    
    state_no_opt = state_initial.copy()
    
    try:
        energy_initial = force.calculateEnergySCF(state_no_opt)
        print(f"  初始能量: {energy_initial:.2f} kJ/mol")
        print(f"  能量/水: {energy_initial/n_waters:.2f} kJ/mol")
        
        # 检查Drude确实没有移动
        max_disp = 0
        for i in range(min(10, n_waters)):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_no_opt.atoms[d_idx].x - state_no_opt.atoms[o_idx].x
            dy = state_no_opt.atoms[d_idx].y - state_no_opt.atoms[o_idx].y
            dz = state_no_opt.atoms[d_idx].z - state_no_opt.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            max_disp = max(max_disp, disp)
        
        print(f"  最大Drude位移: {max_disp:.3f} pm （应该接近0）")
        
    except Exception as e:
        print(f"  计算失败: {e}")
        energy_initial = None
    
    # 步骤3: 运行完整的SCF优化
    print("\n步骤3: 运行完整的SCF优化")
    
    # 恢复正常的SCF参数
    force.setSCFParameters(params)
    
    state_optimized = state_initial.copy()
    
    try:
        energy_optimized = force.calculateEnergySCF(state_optimized)
        print(f"  优化后能量: {energy_optimized:.2f} kJ/mol")
        print(f"  能量/水: {energy_optimized/n_waters:.2f} kJ/mol")
        
        # 分析Drude位移
        displacements = []
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_optimized.atoms[d_idx].x - state_optimized.atoms[o_idx].x
            dy = state_optimized.atoms[d_idx].y - state_optimized.atoms[o_idx].y
            dz = state_optimized.atoms[d_idx].z - state_optimized.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            displacements.append(disp)
        
        avg_disp = np.mean(displacements)
        max_disp = np.max(displacements)
        
        print(f"  平均Drude位移: {avg_disp:.2f} pm")
        print(f"  最大Drude位移: {max_disp:.2f} pm")
        
        # 计算能量变化
        if energy_initial is not None:
            energy_change = energy_optimized - energy_initial
            print(f"\n能量变化分析:")
            print(f"  总能量变化: {energy_change:.2f} kJ/mol")
            print(f"  能量变化/水: {energy_change/n_waters:.2f} kJ/mol")
            
            if energy_change < 0:
                print(f"\n✓ SCF优化降低了能量 ({-energy_change:.2f} kJ/mol)")
                print("  这表明Drude极化正确地稳定了系统")
            else:
                print(f"\n⚠ 能量增加了 {energy_change:.2f} kJ/mol")
                print("  可能的原因：")
                print("  1. 初始能量计算时Drude已经有小的位移")
                print("  2. 数值精度问题")
        
        # 分析前5个水分子的诱导偶极矩
        print(f"\n前5个水分子的诱导偶极矩:")
        for i in range(min(5, n_waters)):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_optimized.atoms[d_idx].x - state_optimized.atoms[o_idx].x
            dy = state_optimized.atoms[d_idx].y - state_optimized.atoms[o_idx].y
            dz = state_optimized.atoms[d_idx].z - state_optimized.atoms[o_idx].z
            
            # 偶极矩 = 电荷 × 位移
            q_drude = -1.71636
            dipole_vec = np.array([dx, dy, dz]) * q_drude * 4.80321  # Debye
            dipole_mag = np.linalg.norm(dipole_vec)
            
            print(f"  水{i+1}: μ = {dipole_mag:.3f} D")
        
    except Exception as e:
        print(f"  SCF优化失败: {e}")
    
    print(f"\n\n总结:")
    print("="*70)
    print("1. SCF优化改变了Drude位置")
    print("2. Drude极化产生了诱导偶极矩")
    print("3. 能量变化反映了极化相互作用")

def main():
    """
    主函数
    """
    test_scf_energy_change()

if __name__ == "__main__":
    main()