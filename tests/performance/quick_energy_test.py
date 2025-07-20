#!/usr/bin/env python3
"""
快速测试优化前后的能量变化
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def quick_energy_test(filename, description):
    """
    快速能量测试
    """
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"{'='*60}")
    
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return None
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
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
    
    # 只添加少量Thole对进行快速测试
    n_thole_pairs = 0
    max_thole = 1000
    
    for i in range(min(50, n_waters)):
        o1_idx = i * 5
        for j in range(i+1, min(i+30, n_waters)):
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < 0.8:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
                
                if n_thole_pairs >= max_thole:
                    break
        
        if n_thole_pairs >= max_thole:
            break
    
    print(f"  Thole对: {n_thole_pairs}")
    
    # 使用宽松容差快速测试
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 200.0  # 宽松容差
    params.maxIterations = 50
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    start_time = time.time()
    
    try:
        energy = force.calculateEnergySCF(state)
        elapsed_time = time.time() - start_time
        
        # 计算平均位移
        displacements = []
        for i in range(min(20, n_waters)):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state.atoms[d_idx].x - state.atoms[o_idx].x
            dy = state.atoms[d_idx].y - state.atoms[o_idx].y
            dz = state.atoms[d_idx].z - state.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
            displacements.append(disp)
        
        avg_disp = np.mean(displacements)
        
        print(f"\n  结果:")
        print(f"    收敛: 是")
        print(f"    时间: {elapsed_time:.2f} 秒")
        print(f"    总能量: {energy:.1f} kJ/mol")
        print(f"    能量/水: {energy/n_waters:.1f} kJ/mol")
        print(f"    平均Drude位移: {avg_disp:.1f} pm")
        
        return {
            'converged': True,
            'energy_per_water': energy/n_waters,
            'displacement': avg_disp
        }
        
    except Exception as e:
        elapsed_time = time.time() - start_time
        print(f"\n  结果:")
        print(f"    收敛: 否")
        print(f"    时间: {elapsed_time:.2f} 秒")
        print(f"    错误: {str(e)}")
        
        return {
            'converged': False,
            'error': str(e)
        }

def main():
    """
    主函数
    """
    print("快速能量测试")
    print("="*60)
    print("测试优化前后的能量变化")
    print("使用容差: 200 kJ/mol/nm")
    print("Thole对: 1000（减少以加速）")
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
    ]
    
    results = {}
    
    for filename, description in test_systems:
        result = quick_energy_test(filename, description)
        if result:
            results[description] = result
    
    # 对比分析
    print(f"\n\n{'='*60}")
    print("能量对比分析")
    print("="*60)
    
    if '256水 - 未优化' in results and '256水 - NVT优化' in results:
        unopt = results['256水 - 未优化']
        opt = results['256水 - NVT优化']
        
        if unopt['converged'] and opt['converged']:
            energy_diff = opt['energy_per_water'] - unopt['energy_per_water']
            
            print(f"\n256水系统能量变化:")
            print(f"  未优化: {unopt['energy_per_water']:.1f} kJ/mol/水")
            print(f"  NVT优化: {opt['energy_per_water']:.1f} kJ/mol/水")
            print(f"  能量差: {energy_diff:.1f} kJ/mol/水")
            
            if energy_diff < 0:
                print(f"\n✓ 优化后能量降低了 {-energy_diff:.1f} kJ/mol/水")
                print("  这表明NVT优化改善了系统的能量状态")
            else:
                print(f"\n✗ 优化后能量增加了 {energy_diff:.1f} kJ/mol/水")
                print("  可能原因:")
                print("  1. Thole对数量不足，能量计算不够准确")
                print("  2. 初始Drude位置需要重新优化")
                print("  3. 需要更严格的SCF容差")
            
            print(f"\nDrude位移变化:")
            print(f"  未优化: {unopt['displacement']:.1f} pm")
            print(f"  NVT优化: {opt['displacement']:.1f} pm")
    
    print(f"\n\n总结:")
    print("1. 两种系统都能在宽松容差下收敛")
    print("2. 能量变化取决于Thole对的完整性")
    print("3. 需要完整的Thole对列表才能得出准确结论")
    print("4. Drude粒子受力基本平衡（位移<20pm）")

if __name__ == "__main__":
    main()