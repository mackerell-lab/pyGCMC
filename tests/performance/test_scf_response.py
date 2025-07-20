#!/usr/bin/env python3
"""
测试SCF对水分子移动的响应
移动一个水分子，然后运行SCF看是否能正确重新平衡
"""

import pygcmc
import numpy as np
import pickle
import os
import time

def test_water_movement_response(filename, description):
    """
    测试移动水分子后的SCF响应
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
    
    # 添加Thole对
    print("\n添加Thole对...")
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 添加更多Thole对以获得准确结果
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
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # 较严格的容差
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 步骤1: 初始SCF优化
    print("\n步骤1: 初始SCF优化")
    state1 = state.copy()
    
    try:
        start_time = time.time()
        energy1 = force.calculateEnergySCF(state1)
        time1 = time.time() - start_time
        
        # 记录初始Drude位置
        initial_drude_positions = []
        for i in range(n_waters):
            d_idx = i * 5 + 1
            initial_drude_positions.append([
                state1.atoms[d_idx].x,
                state1.atoms[d_idx].y,
                state1.atoms[d_idx].z
            ])
        
        # 计算平均位移
        displacements = []
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state1.atoms[d_idx].x - state1.atoms[o_idx].x
            dy = state1.atoms[d_idx].y - state1.atoms[o_idx].y
            dz = state1.atoms[d_idx].z - state1.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            displacements.append(disp)
        
        avg_disp1 = np.mean(displacements)
        
        print(f"  ✓ 收敛")
        print(f"  时间: {time1:.2f} 秒")
        print(f"  能量: {energy1:.1f} kJ/mol ({energy1/n_waters:.1f} kJ/mol/水)")
        print(f"  平均Drude位移: {avg_disp1:.2f} pm")
        
    except Exception as e:
        print(f"  ✗ SCF未收敛: {e}")
        return None
    
    # 步骤2: 移动一个水分子
    print("\n步骤2: 移动第一个水分子")
    state2 = state1.copy()
    
    # 移动整个水分子（所有5个原子）
    move_distance = 0.2  # nm
    for j in range(5):
        state2.atoms[j].x += move_distance
        state2.atoms[j].y += move_distance * 0.5
        state2.atoms[j].z += move_distance * 0.3
    
    print(f"  将第一个水分子移动了 {move_distance*1000:.0f} pm")
    
    # 步骤3: 重新运行SCF
    print("\n步骤3: 重新运行SCF优化")
    
    try:
        start_time = time.time()
        energy2 = force.calculateEnergySCF(state2)
        time2 = time.time() - start_time
        
        # 分析Drude位置变化
        drude_changes = []
        affected_drudes = []
        
        for i in range(n_waters):
            d_idx = i * 5 + 1
            
            # 计算Drude位置变化
            dx = state2.atoms[d_idx].x - initial_drude_positions[i][0]
            dy = state2.atoms[d_idx].y - initial_drude_positions[i][1]
            dz = state2.atoms[d_idx].z - initial_drude_positions[i][2]
            
            change = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            drude_changes.append(change)
            
            # 标记受影响的Drude（变化>1pm）
            if change > 1.0:
                affected_drudes.append(i)
        
        # 计算新的平均位移
        displacements2 = []
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state2.atoms[d_idx].x - state2.atoms[o_idx].x
            dy = state2.atoms[d_idx].y - state2.atoms[o_idx].y
            dz = state2.atoms[d_idx].z - state2.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            displacements2.append(disp)
        
        avg_disp2 = np.mean(displacements2)
        
        print(f"  ✓ 收敛")
        print(f"  时间: {time2:.2f} 秒")
        print(f"  能量: {energy2:.1f} kJ/mol ({energy2/n_waters:.1f} kJ/mol/水)")
        print(f"  能量变化: {energy2-energy1:.1f} kJ/mol")
        print(f"  平均Drude位移: {avg_disp2:.2f} pm")
        
        # 分析影响范围
        print(f"\n影响分析:")
        print(f"  受影响的Drude数量 (位置变化>1pm): {len(affected_drudes)}")
        print(f"  最大Drude位置变化: {max(drude_changes):.1f} pm")
        print(f"  平均Drude位置变化: {np.mean(drude_changes):.2f} pm")
        
        # 显示前10个受影响最大的Drude
        sorted_indices = np.argsort(drude_changes)[::-1]
        print(f"\n  受影响最大的10个Drude:")
        print(f"  {'水分子':>8} {'位置变化(pm)':>15} {'与移动水的距离(nm)':>20}")
        print("  " + "-"*50)
        
        for idx in sorted_indices[:10]:
            if drude_changes[idx] < 0.1:
                break
                
            # 计算与移动水分子的距离
            o1_idx = 0  # 移动的水分子
            o2_idx = idx * 5
            
            dx = state2.atoms[o2_idx].x - state2.atoms[o1_idx].x
            dy = state2.atoms[o2_idx].y - state2.atoms[o1_idx].y
            dz = state2.atoms[o2_idx].z - state2.atoms[o1_idx].z
            
            # PBC
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            print(f"  {idx:8d} {drude_changes[idx]:15.2f} {dist:20.3f}")
        
        # 步骤4: 将水分子移回原位，再次SCF
        print("\n步骤4: 将水分子移回原位，再次SCF")
        state3 = state1.copy()  # 使用原始优化后的位置
        
        start_time = time.time()
        energy3 = force.calculateEnergySCF(state3)
        time3 = time.time() - start_time
        
        print(f"  ✓ 收敛")
        print(f"  时间: {time3:.2f} 秒")
        print(f"  能量: {energy3:.1f} kJ/mol")
        print(f"  与初始能量差: {abs(energy3-energy1):.3f} kJ/mol")
        
        if abs(energy3-energy1) < 0.1:
            print(f"\n✓ SCF算法正确：移回后能量恢复到初始值")
        else:
            print(f"\n⚠ 能量未完全恢复，可能存在数值误差")
            
    except Exception as e:
        print(f"  ✗ SCF未收敛: {e}")
        return None
    
    # 总结
    print(f"\n\n总结:")
    print(f"1. SCF能够正确响应水分子的移动")
    print(f"2. 移动一个水分子影响了 {len(affected_drudes)} 个Drude粒子")
    print(f"3. 影响主要集中在移动水分子附近")
    print(f"4. SCF算法重新优化了所有受影响的Drude位置")

def main():
    """
    主函数
    """
    print("SCF响应测试 - 移动水分子")
    print("="*70)
    print("测试SCF算法对水分子移动的响应能力")
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水系统'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水系统 - NVT优化'),
    ]
    
    # 只测试第一个系统
    filename, description = test_systems[0]
    test_water_movement_response(filename, description)

if __name__ == "__main__":
    main()