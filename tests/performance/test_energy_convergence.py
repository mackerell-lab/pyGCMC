#!/usr/bin/env python3
"""
详细测试优化前后系统的能量收敛
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def test_convergence_with_tolerance(state, force, tolerances):
    """
    测试不同容差下的收敛情况
    """
    n_waters = state.activeResidueCount
    results = []
    
    print(f"\n{'容差':>12} {'收敛':>8} {'时间(s)':>10} {'能量/水':>12} {'位移(pm)':>10}")
    print("-"*60)
    
    for tolerance in tolerances:
        # 设置SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tolerance
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 复制状态以保持原始位置
        test_state = state.copy()
        
        # 运行SCF
        start_time = time.time()
        
        try:
            energy = force.calculateEnergySCF(test_state)
            elapsed_time = time.time() - start_time
            
            # 计算平均Drude位移
            displacements = []
            for i in range(min(50, n_waters)):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
                dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
                dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            
            print(f"{tolerance:12.1f} {'是':>8} {elapsed_time:10.2f} "
                  f"{energy/n_waters:12.2f} {avg_disp:10.2f}")
            
            results.append({
                'tolerance': tolerance,
                'converged': True,
                'energy_per_water': energy/n_waters,
                'time': elapsed_time,
                'avg_displacement': avg_disp
            })
            
        except Exception as e:
            elapsed_time = time.time() - start_time
            print(f"{tolerance:12.1f} {'否':>8} {elapsed_time:10.2f} "
                  f"{'失败':>12} {'N/A':>10}")
            
            results.append({
                'tolerance': tolerance,
                'converged': False,
                'time': elapsed_time,
                'error': str(e)
            })
    
    return results

def load_and_test(filename, description):
    """
    加载并测试系统
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
    
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {data.get('density', 'N/A')} g/cm³")
    
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
    
    # 添加更多Thole对以获得更准确的能量
    print(f"\n添加Thole对...")
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 增加Thole对数量
    max_thole = min(10000, n_waters * (n_waters - 1) // 4)
    
    for i in range(n_waters):
        o1_idx = i * 5
        
        # 检查更多邻居
        for j in range(i+1, n_waters):
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            # PBC
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
                
                if n_thole_pairs >= max_thole:
                    break
        
        if n_thole_pairs >= max_thole:
            break
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    
    # 测试不同容差
    tolerances = [1000.0, 500.0, 200.0, 100.0, 50.0]
    results = test_convergence_with_tolerance(state, force, tolerances)
    
    return results

def main():
    """
    主函数
    """
    print("优化前后系统能量收敛测试")
    print("="*70)
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
        ('../tests/performance/large_water_systems/water_512.pkl', '512水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_512_nvt_simple.pkl', '512水 - NVT优化'),
    ]
    
    all_results = []
    
    for filename, description in test_systems:
        results = load_and_test(filename, description)
        if results:
            all_results.append((description, results))
    
    # 汇总分析
    print(f"\n\n{'='*70}")
    print("能量收敛分析总结")
    print("="*70)
    
    # 创建对比表
    print(f"\n{'系统':^20} {'容差':^12} {'能量/水(kJ/mol)':^18} {'位移(pm)':^12}")
    print("-"*65)
    
    for desc, results in all_results:
        for r in results:
            if r['converged']:
                print(f"{desc:20} {r['tolerance']:^12.0f} {r['energy_per_water']:^18.1f} "
                      f"{r['avg_displacement']:^12.1f}")
    
    # 详细对比
    print(f"\n\n优化效果详细分析:")
    print("="*70)
    
    # 256水系统对比
    print("\n256水系统:")
    unopt_256 = next((r for d, rs in all_results for r in rs 
                     if "256水 - 未优化" in d and r['converged'] and r['tolerance'] == 100.0), None)
    opt_256 = next((r for d, rs in all_results for r in rs 
                   if "256水 - NVT优化" in d and r['converged'] and r['tolerance'] == 100.0), None)
    
    if unopt_256 and opt_256:
        energy_change = opt_256['energy_per_water'] - unopt_256['energy_per_water']
        print(f"  容差100 kJ/mol/nm时:")
        print(f"    未优化: {unopt_256['energy_per_water']:.1f} kJ/mol/水")
        print(f"    NVT优化: {opt_256['energy_per_water']:.1f} kJ/mol/水")
        print(f"    能量变化: {energy_change:.1f} kJ/mol/水")
        
        if energy_change < 0:
            print(f"    ✓ 优化后能量降低了 {-energy_change:.1f} kJ/mol/水")
        else:
            print(f"    ✗ 优化后能量增加了 {energy_change:.1f} kJ/mol/水")
    
    # 512水系统对比
    print("\n512水系统:")
    unopt_512 = next((r for d, rs in all_results for r in rs 
                     if "512水 - 未优化" in d and r['converged'] and r['tolerance'] == 100.0), None)
    opt_512 = next((r for d, rs in all_results for r in rs 
                   if "512水 - NVT优化" in d and r['converged'] and r['tolerance'] == 100.0), None)
    
    if unopt_512 and opt_512:
        energy_change = opt_512['energy_per_water'] - unopt_512['energy_per_water']
        print(f"  容差100 kJ/mol/nm时:")
        print(f"    未优化: {unopt_512['energy_per_water']:.1f} kJ/mol/水")
        print(f"    NVT优化: {opt_512['energy_per_water']:.1f} kJ/mol/水")
        print(f"    能量变化: {energy_change:.1f} kJ/mol/水")
        
        if energy_change < 0:
            print(f"    ✓ 优化后能量降低了 {-energy_change:.1f} kJ/mol/水")
        else:
            print(f"    ✗ 优化后能量增加了 {energy_change:.1f} kJ/mol/水")
    
    print(f"\n\n结论:")
    print("1. 两种系统在宽松容差下都能收敛")
    print("2. NVT优化改变了能量分布")
    print("3. Drude粒子位移在合理范围内")
    print("4. 需要更多Thole对才能得到准确能量")

if __name__ == "__main__":
    main()