#!/usr/bin/env python3
"""
比较不同Drude算法在大型水系统中的性能
"""

import numpy as np
import time
import pickle
import os
import pygcmc

def benchmark_drude_algorithms():
    """
    测试不同算法的性能
    """
    print("Drude算法性能基准测试")
    print("="*70)
    
    # 加载优化的水系统
    filenames = {
        256: '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl',
        512: '../tests/performance/optimized_water_systems/water_512_nvt_simple.pkl',
        1024: '../tests/performance/optimized_water_systems/water_1024_nvt_simple.pkl',
        2048: '../tests/performance/optimized_water_systems/water_2048_nvt_simple.pkl'
    }
    
    # 测试不同大小的系统
    system_sizes = [256, 512, 1024, 2048]
    
    # 可用的算法
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF"),
        (pygcmc.DrudeAlgorithm.ConjugateGradient, "CG (共轭梯度)"),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3"),
        (pygcmc.DrudeAlgorithm.OPT4, "OPT4"),
        (pygcmc.DrudeAlgorithm.SmartOPT3, "SmartOPT3"),
        (pygcmc.DrudeAlgorithm.AdaptiveOPT, "AdaptiveOPT"),
        (pygcmc.DrudeAlgorithm.HybridOPT, "HybridOPT"),
        (pygcmc.DrudeAlgorithm.FBP, "FBP (力平衡)")
    ]
    
    results = {}
    
    for n_waters in system_sizes:
        filename = filenames.get(n_waters)
        if not filename or not os.path.exists(filename):
            print(f"\n跳过{n_waters}水系统（文件不存在）")
            continue
            
        print(f"\n\n测试{n_waters}水系统")
        print("-"*60)
        
        # 加载系统
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        
        positions = data['positions']
        box_length = data['box_length']
        charges = data['charges']
        
        # 创建状态
        state = create_state(n_waters, positions, box_length, charges)
        
        results[n_waters] = {}
        
        # 测试每种算法
        for algo, algo_name in algorithms:
            print(f"\n测试{algo_name}...")
            
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
            
            # 添加适量的Thole对
            n_thole_pairs = add_thole_pairs(force, state, n_waters, box_length)
            
            # 设置算法
            force.setAlgorithm(algo)
            
            # 设置参数
            if algo == pygcmc.DrudeAlgorithm.SCF:
                params = pygcmc.DrudeSCFParams()
                params.tolerance = 10.0
                params.maxIterations = 100
                params.dampingFactor = 0.5
                params.maxDrudeDistance = 0.02
                force.setSCFParameters(params)
            
            # 准备测试状态
            state_test = state.copy()
            reset_drude_positions(state_test, n_waters)
            
            # 计时
            try:
                start_time = time.time()
                
                # 运行5次取平均
                n_runs = 5
                energies = []
                for _ in range(n_runs):
                    state_run = state_test.copy()
                    if algo == pygcmc.DrudeAlgorithm.SCF:
                        energy = force.calculateEnergySCF(state_run)
                    elif algo == pygcmc.DrudeAlgorithm.ConjugateGradient:
                        energy = force.calculateEnergyCG(state_run)
                    else:
                        energy = force.calculateEnergyOPT(state_run)
                    energies.append(energy)
                
                end_time = time.time()
                avg_time = (end_time - start_time) / n_runs
                avg_energy = np.mean(energies)
                
                # 检查收敛性
                avg_disp = calculate_avg_displacement(state_run, n_waters)
                
                results[n_waters][algo_name] = {
                    'time': avg_time,
                    'energy': avg_energy,
                    'displacement': avg_disp,
                    'thole_pairs': n_thole_pairs
                }
                
                print(f"  平均时间: {avg_time*1000:.2f} ms")
                print(f"  能量: {avg_energy/n_waters:.2f} kJ/mol/水")
                print(f"  平均位移: {avg_disp:.2f} pm")
                
            except Exception as e:
                print(f"  失败: {e}")
                results[n_waters][algo_name] = {'error': str(e)}
    
    # 总结结果
    print_summary(results)
    
    return results

def create_state(n_waters, positions, box_length, charges):
    """
    创建PyGCMC状态
    """
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
    
    return state

def add_thole_pairs(force, state, n_waters, box_length):
    """
    添加Thole对
    """
    n_thole_pairs = 0
    thole_cutoff = 0.8
    max_pairs_per_molecule = 20  # 限制每个分子的Thole对数
    
    for i in range(n_waters):
        o1_idx = i * 5
        pairs_for_this_molecule = 0
        
        for j in range(i+1, n_waters):
            if pairs_for_this_molecule >= max_pairs_per_molecule:
                break
                
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
                pairs_for_this_molecule += 1
    
    return n_thole_pairs

def reset_drude_positions(state, n_waters):
    """
    重置Drude到parent位置
    """
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state.atoms[d_idx].x = state.atoms[o_idx].x
        state.atoms[d_idx].y = state.atoms[o_idx].y
        state.atoms[d_idx].z = state.atoms[o_idx].z

def calculate_avg_displacement(state, n_waters):
    """
    计算平均Drude位移
    """
    displacements = []
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        displacements.append(disp)
    
    return np.mean(displacements)

def print_summary(results):
    """
    打印总结
    """
    print("\n\n性能总结")
    print("="*70)
    
    # 创建表格
    print(f"\n{'系统大小':>10} {'算法':>15} {'时间(ms)':>12} {'相对速度':>10} {'位移(pm)':>10}")
    print("-"*60)
    
    for n_waters in sorted(results.keys()):
        # 找出最快的算法
        valid_results = {k: v for k, v in results[n_waters].items() 
                        if 'time' in v and not 'error' in v}
        
        if not valid_results:
            continue
            
        fastest_time = min(v['time'] for v in valid_results.values())
        
        for algo_name in ['SCF', 'CG (共轭梯度)', 'OPT3', 'OPT4', 'SmartOPT3', 'AdaptiveOPT', 'HybridOPT', 'FBP (力平衡)']:
            if algo_name in results[n_waters]:
                result = results[n_waters][algo_name]
                
                if 'error' in result:
                    print(f"{n_waters:>10} {algo_name:>15} {'失败':>12} {'-':>10} {'-':>10}")
                else:
                    time_ms = result['time'] * 1000
                    relative_speed = result['time'] / fastest_time
                    displacement = result.get('displacement', 0)
                    
                    print(f"{n_waters:>10} {algo_name:>15} {time_ms:>12.2f} {relative_speed:>10.2f}x {displacement:>10.2f}")
    
    # 推荐
    print("\n推荐：")
    print("1. 小系统（<500水）：SCF或OPT3，精度高")
    print("2. 中等系统（500-2000水）：OPT3或CG，平衡速度和精度")
    print("3. 大系统（>2000水）：FBP或OPT2，速度快")
    print("4. 需要高精度：始终使用SCF")

def main():
    """
    主函数
    """
    results = benchmark_drude_algorithms()
    
    # 测试并行化潜力
    print("\n\n并行化分析")
    print("="*70)
    print("Drude计算的并行化潜力：")
    print("1. SCF/CG：迭代算法，难以并行化")
    print("2. OPT系列：每个Drude独立计算，易于并行化")
    print("3. FBP：力计算可并行，位置更新需同步")
    print("\n建议：对OPT3实现OpenMP并行化")

if __name__ == "__main__":
    main()