#!/usr/bin/env python3
"""
简化测试：检查SCF优化后Drude粒子是否达到平衡
通过检查优化前后的位移变化来判断
"""

import pygcmc
import numpy as np
import pickle
import os

def check_drude_equilibrium(filename, description):
    """
    检查Drude粒子平衡状态
    """
    print(f"\n{'='*70}")
    print(f"{description}")
    print(f"{'='*70}")
    
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return None
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"系统信息:")
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
    
    # 添加少量Thole对进行快速测试
    n_thole_pairs = 0
    for i in range(min(50, n_waters)):
        for j in range(i+1, min(i+20, n_waters)):
            force.addScreenedPair(i, j, 1.3)
            n_thole_pairs += 1
            if n_thole_pairs >= 500:
                break
        if n_thole_pairs >= 500:
            break
    
    print(f"  Thole对数: {n_thole_pairs}")
    
    # 记录初始Drude位置
    initial_drude_positions = []
    for i in range(n_waters):
        d_idx = i * 5 + 1
        initial_drude_positions.append([
            state.atoms[d_idx].x,
            state.atoms[d_idx].y,
            state.atoms[d_idx].z
        ])
    
    # 设置SCF参数 - 使用不同的容差测试
    tolerances = [500.0, 200.0, 100.0, 50.0]
    
    print(f"\n测试不同容差下的平衡状态:")
    print(f"{'容差(kJ/mol/nm)':>18} {'收敛':>8} {'能量/水':>12} {'平均位移(pm)':>15} {'最大位移(pm)':>15}")
    print("-"*80)
    
    for tolerance in tolerances:
        # 复制原始状态
        test_state = state.copy()
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tolerance
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        try:
            energy = force.calculateEnergySCF(test_state)
            
            # 计算位移
            displacements = []
            for i in range(n_waters):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
                dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
                dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            max_disp = np.max(displacements)
            
            print(f"{tolerance:18.1f} {'是':>8} {energy/n_waters:12.1f} {avg_disp:15.2f} {max_disp:15.2f}")
            
        except Exception as e:
            print(f"{tolerance:18.1f} {'否':>8} {'N/A':>12} {'N/A':>15} {'N/A':>15}")
    
    # 测试力平衡的一个简单方法：检查两次SCF的结果是否一致
    print(f"\n\n一致性测试（运行SCF两次）:")
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 100.0
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 第一次SCF
    state1 = state.copy()
    try:
        energy1 = force.calculateEnergySCF(state1)
        
        # 记录Drude位置
        drude_pos_1 = []
        for i in range(min(20, n_waters)):  # 只检查前20个
            d_idx = i * 5 + 1
            drude_pos_1.append([
                state1.atoms[d_idx].x,
                state1.atoms[d_idx].y,
                state1.atoms[d_idx].z
            ])
        
        # 第二次SCF（从第一次的结果开始）
        state2 = state1.copy()
        energy2 = force.calculateEnergySCF(state2)
        
        # 检查Drude位置变化
        max_change = 0.0
        for i in range(min(20, n_waters)):
            d_idx = i * 5 + 1
            dx = state2.atoms[d_idx].x - drude_pos_1[i][0]
            dy = state2.atoms[d_idx].y - drude_pos_1[i][1]
            dz = state2.atoms[d_idx].z - drude_pos_1[i][2]
            change = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            max_change = max(max_change, change)
        
        print(f"  第一次SCF能量: {energy1/n_waters:.1f} kJ/mol/水")
        print(f"  第二次SCF能量: {energy2/n_waters:.1f} kJ/mol/水")
        print(f"  能量差: {abs(energy2-energy1)/n_waters:.3f} kJ/mol/水")
        print(f"  最大Drude位移变化: {max_change:.3f} pm")
        
        if max_change < 0.1:
            print(f"\n✓ Drude粒子已达到平衡（位置变化 < 0.1 pm）")
        elif max_change < 1.0:
            print(f"\n✓ Drude粒子基本平衡（位置变化 < 1 pm）")
        else:
            print(f"\n⚠ Drude粒子可能未完全平衡（位置变化 = {max_change:.1f} pm）")
            
    except Exception as e:
        print(f"  SCF未收敛: {e}")
    
    # 分析Drude位移分布
    print(f"\n\nDrude位移分布分析:")
    
    # 使用最后一次成功的状态
    if 'state1' in locals():
        test_state = state1
    else:
        test_state = state.copy()
        params.tolerance = 500.0  # 非常宽松
        try:
            force.calculateEnergySCF(test_state)
        except:
            pass
    
    displacements = []
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
        dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
        dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        displacements.append(disp)
    
    # 统计
    bins = [0, 2, 5, 10, 15, 20, 30, 50]
    hist, _ = np.histogram(displacements, bins=bins)
    
    print(f"  位移范围(pm)    数量    百分比")
    print("-"*40)
    for i in range(len(bins)-1):
        percentage = hist[i] / n_waters * 100
        print(f"  {bins[i]:2d}-{bins[i+1]:2d}         {hist[i]:4d}    {percentage:5.1f}%")
    
    # 最终判断
    avg_disp = np.mean(displacements)
    max_disp = np.max(displacements)
    
    print(f"\n总体统计:")
    print(f"  平均位移: {avg_disp:.2f} pm")
    print(f"  最大位移: {max_disp:.2f} pm")
    print(f"  标准差: {np.std(displacements):.2f} pm")
    
    if avg_disp < 10 and max_disp < 20:
        print(f"\n✓ 结论：Drude粒子整体达到良好平衡")
    elif avg_disp < 15 and max_disp < 30:
        print(f"\n✓ 结论：Drude粒子基本平衡")
    else:
        print(f"\n⚠ 结论：部分Drude粒子可能未完全平衡")

def main():
    """
    主函数
    """
    print("Drude粒子平衡状态检查")
    print("="*70)
    print("通过SCF优化和一致性测试检查力平衡")
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
    ]
    
    for filename, description in test_systems:
        check_drude_equilibrium(filename, description)

if __name__ == "__main__":
    main()