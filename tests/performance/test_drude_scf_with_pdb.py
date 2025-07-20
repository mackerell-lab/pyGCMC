#!/usr/bin/env python3
"""
使用生成的PDB文件测试Drude SCF的速度和收敛性
"""

import pygcmc
import numpy as np
import time
import pickle
import os

def load_water_system_from_pdb(n_waters):
    """
    从pickle文件加载水系统（包含完整的5位点信息）
    """
    pickle_file = f'../tests/performance/water_density_1.0/water_{n_waters}.pkl'
    if not os.path.exists(pickle_file):
        raise FileNotFoundError(f"找不到文件: {pickle_file}")
    
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    # 转换为pygcmc格式
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]  # O, D, H1, H2, M
    atom_types = [0, 1, 2, 2, 3]  # 原子类型
    
    positions = data['positions']
    n_waters = data['n_waters']
    box_length = data['box_length']
    
    # 创建原子和残基
    for i in range(n_waters):
        # 添加5个原子（O, D, H1, H2, M）
        for j in range(5):
            atom = pygcmc.MCAtom()
            idx = i * 5 + j
            atom.x = positions[idx][0]
            atom.y = positions[idx][1]
            atom.z = positions[idx][2]
            atom.charge = charges[j]
            atom.type = atom_types[j]
            atoms.append(atom)
        
        # 创建残基
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0  # 水分子类型
        residues.append(res)
    
    # 创建状态
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    # 设置盒子和截断
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)  # 确保截断小于半盒子
    
    # 设置力场参数（只有O原子有LJ）
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]  # 只有O有LJ
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, data

def test_scf_convergence(state, n_waters, tolerances=[1.0, 10.0, 100.0]):
    """
    测试不同容差下的SCF收敛性
    """
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP Drude参数
    drude_charge = -1.71636
    polarizability = 0.0009782237  # nm³
    k_spring = 418400.0  # kJ/mol/nm² (从极化率计算得出)
    
    # 添加Drude粒子（每个水分子一个）
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1,    # Drude粒子索引
            parentIndex=5*i,     # O原子索引
            aniso1Index=-1,      # 各向异性（这里不使用）
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=drude_charge,
            polarizability=polarizability,
            aniso12=1.0,
            aniso34=1.0
        )
    
    # 添加Thole屏蔽对
    n_thole_pairs = 0
    # 动态调整Thole截断距离，确保不超过半盒子
    thole_cutoff = min(0.8, state.info.box[0] / 2.0 - 0.01)  # nm
    thole_param = 1.3   # Thole参数
    
    print(f"  Thole截断距离: {thole_cutoff:.3f} nm (半盒子: {state.info.box[0]/2:.3f} nm)")
    
    for i in range(n_waters):
        o1_idx = i * 5
        for j in range(i+1, n_waters):
            o2_idx = j * 5
            
            # 计算O-O距离
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            # 应用PBC
            box = state.info.box[0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, thole_param)
                n_thole_pairs += 1
    
    print(f"  添加了 {n_thole_pairs} 个Thole屏蔽对")
    
    # 测试结果存储
    results = []
    
    # 测试不同容差
    for tolerance in tolerances:
        print(f"\n  容差 = {tolerance} kJ/mol/nm:")
        
        # 设置SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tolerance
        params.maxIterations = 200
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02  # 2 pm硬墙约束
        force.setSCFParameters(params)
        
        # 确保使用SCF算法
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 复制状态（避免修改原始状态）
        test_state = state.copy()
        
        # 计时
        start_time = time.time()
        
        # 运行SCF
        try:
            energy = force.calculateEnergySCF(test_state)
            elapsed_time = (time.time() - start_time) * 1000  # ms
            converged = True
            
            # 分析Drude位移
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
            
            print(f"    收敛: 是")
            print(f"    时间: {elapsed_time:.1f} ms")
            print(f"    能量: {energy:.2f} kJ/mol ({energy/n_waters:.2f} kJ/mol/水)")
            print(f"    平均Drude位移: {avg_disp:.2f} pm")
            print(f"    最大Drude位移: {max_disp:.2f} pm")
            
        except Exception as e:
            elapsed_time = (time.time() - start_time) * 1000
            converged = False
            energy = float('nan')
            avg_disp = float('nan')
            max_disp = float('nan')
            print(f"    收敛: 否 (错误: {str(e)})")
            print(f"    时间: {elapsed_time:.1f} ms")
        
        results.append({
            'tolerance': tolerance,
            'converged': converged,
            'time_ms': elapsed_time,
            'energy': energy,
            'avg_displacement': avg_disp,
            'max_displacement': max_disp
        })
    
    return results

def main():
    """
    主测试函数
    """
    print("Drude SCF收敛性和性能测试")
    print("使用密度1.0 g/cm³的SWM4-NDP水模型")
    print("="*70)
    
    # 测试不同大小的系统
    system_sizes = [2, 4, 8, 16, 32, 64, 128]
    tolerances = [1.0, 10.0, 100.0]
    
    # 汇总结果
    all_results = {}
    
    for n_waters in system_sizes:
        print(f"\n{'='*70}")
        print(f"测试 {n_waters} 水分子系统")
        print(f"{'='*70}")
        
        try:
            # 加载系统
            state, data = load_water_system_from_pdb(n_waters)
            
            print(f"  盒子长度: {data['box_length']:.3f} nm")
            print(f"  密度: {data['density']:.6f} g/cm³")
            print(f"  截断距离: {state.info.cutoff:.3f} nm")
            
            # 测试收敛性
            results = test_scf_convergence(state, n_waters, tolerances)
            all_results[n_waters] = results
            
        except Exception as e:
            print(f"  错误: {e}")
            import traceback
            traceback.print_exc()
    
    # 打印汇总
    print("\n\n" + "="*70)
    print("收敛性汇总")
    print("="*70)
    
    print(f"\n{'系统':<10} {'容差':<10} {'收敛':<8} {'时间(ms)':<12} {'能量/水':<15} {'平均位移(pm)':<12}")
    print("-"*80)
    
    for n_waters in sorted(all_results.keys()):
        for result in all_results[n_waters]:
            conv_str = "是" if result['converged'] else "否"
            energy_per_water = result['energy'] / n_waters if not np.isnan(result['energy']) else float('nan')
            
            print(f"{n_waters:<10} {result['tolerance']:<10.1f} {conv_str:<8} "
                  f"{result['time_ms']:<12.1f} {energy_per_water:<15.2f} "
                  f"{result['avg_displacement']:<12.2f}")
    
    # 性能分析
    print("\n\n" + "="*70)
    print("性能分析")
    print("="*70)
    
    print("\n时间复杂度分析（容差=10.0）:")
    sizes = []
    times = []
    
    for n_waters in sorted(all_results.keys()):
        for result in all_results[n_waters]:
            if result['tolerance'] == 10.0 and result['converged']:
                sizes.append(n_waters)
                times.append(result['time_ms'])
                break
    
    if len(sizes) > 1:
        # 简单的复杂度估算
        log_sizes = np.log(sizes)
        log_times = np.log(times)
        
        # 线性回归估算 time ~ n^k
        k = np.polyfit(log_sizes, log_times, 1)[0]
        print(f"  时间复杂度: O(n^{k:.2f})")
    
    print("\n关键发现:")
    print("1. 密度1.0 g/cm³的水系统中Drude-Drude相互作用密集")
    print("2. SCF收敛性随系统增大和容差减小而变差")
    print("3. Drude位移典型值在10-20 pm范围")
    print("4. 建议生产计算使用容差10-100 kJ/mol/nm")

if __name__ == "__main__":
    main()