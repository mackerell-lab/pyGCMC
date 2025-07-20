#!/usr/bin/env python3
"""
最终的CG vs SCF对比测试，使用优化的水体系
"""

import pygcmc
import numpy as np
import time
import pickle
import os

def load_water_system(n_waters):
    """
    加载优化的水体系
    """
    # 尝试不同的目录
    for directory in ['../tests/performance/water_systems_final',
                     '../tests/performance/water_systems',
                     '../tests/performance/optimized_water_systems', 
                     '../tests/performance/water_pdbs', 
                     'optimized_water_systems', 
                     '.']:
        pickle_file = os.path.join(directory, f'water_{n_waters}.pkl')
        if os.path.exists(pickle_file):
            print(f"  加载系统: {pickle_file}")
            with open(pickle_file, 'rb') as f:
                data = pickle.load(f)
            return convert_to_pygcmc(data)
    
    # 如果没找到，创建简单系统
    print(f"  创建新系统 ({n_waters} 水)")
    return create_simple_system(n_waters)

def convert_to_pygcmc(data):
    """
    将保存的数据转换为pygcmc状态
    """
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters):
        for j in range(5):
            atom = pygcmc.MCAtom()
            idx = i * 5 + j
            atom.x = positions[idx][0]
            atom.y = positions[idx][1]
            atom.z = positions[idx][2]
            atom.charge = charges[j]
            atom.type = atom_types[j]
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_simple_system(n_waters):
    """
    创建简单系统（后备方案）
    """
    # 直接使用之前生成的系统
    import subprocess
    import sys
    
    print(f"  生成新系统...")
    result = subprocess.run([sys.executable, '../tests/performance/generate_water_quick.py'], 
                          capture_output=True, text=True)
    
    # 重新加载
    pickle_file = f'../tests/performance/water_systems/water_{n_waters}.pkl'
    if os.path.exists(pickle_file):
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
        return convert_to_pygcmc(data)
    else:
        raise FileNotFoundError(f"无法创建系统: {pickle_file}")

def run_comparison(n_waters):
    """
    运行CG vs SCF对比
    """
    print(f"\n{'='*80}")
    print(f"测试 {n_waters} 水分子系统")
    print(f"{'='*80}")
    
    # 加载系统
    state = load_water_system(n_waters)
    
    # 系统信息
    actual_density = n_waters * 18.015 / (state.info.box[0]**3 * 0.6022)
    print(f"  盒子: {state.info.box[0]:.3f} nm")
    print(f"  密度: {actual_density:.3f} g/cm³")
    print(f"  截断: {state.info.cutoff:.3f} nm")
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0  # kJ/mol/nm²
    polarizability = 0.0009782237  # nm³
    
    # 添加Drude粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 智能添加Thole对
    n_pairs = 0
    cutoff_thole = 0.8  # nm
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            dx = state.atoms[5*i].x - state.atoms[5*j].x
            dy = state.atoms[5*i].y - state.atoms[5*j].y
            dz = state.atoms[5*i].z - state.atoms[5*j].z
            
            # PBC
            box = state.info.box[0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < cutoff_thole:
                force.addScreenedPair(i, j, 1.3)
                n_pairs += 1
    
    print(f"  Thole对: {n_pairs}")
    
    # 测试不同容差
    tolerances = [1.0, 10.0, 100.0]
    results = {}
    
    print(f"\n{'容差':<10} {'算法':<6} {'时间(ms)':<12} {'能量(kJ/mol)':<15} {'能量/水':<12} {'位移(pm)':<10} {'收敛':<6}")
    print("-"*85)
    
    for tol in tolerances:
        # 设置参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 200
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        # SCF测试
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        state_scf = state.copy()
        
        start = time.time()
        energy_scf = force.calculateEnergySCF(state_scf)
        time_scf = (time.time() - start) * 1000
        
        # 分析SCF结果
        disps_scf = []
        for i in range(n_waters):
            dx = state_scf.atoms[5*i+1].x - state_scf.atoms[5*i].x
            dy = state_scf.atoms[5*i+1].y - state_scf.atoms[5*i].y
            dz = state_scf.atoms[5*i+1].z - state_scf.atoms[5*i].z
            disps_scf.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
        
        converged_scf = "是" if max(disps_scf) < 20 else "否"
        
        print(f"{tol:<10.1f} {'SCF':<6} {time_scf:<12.1f} {energy_scf:<15.2f} "
              f"{energy_scf/n_waters:<12.2f} {np.mean(disps_scf):<10.2f} {converged_scf:<6}")
        
        # CG测试
        force.setAlgorithm(pygcmc.DrudeAlgorithm.ConjugateGradient)
        state_cg = state.copy()
        
        start = time.time()
        energy_cg = force.calculateEnergySCF(state_cg)
        time_cg = (time.time() - start) * 1000
        
        # 分析CG结果
        disps_cg = []
        for i in range(n_waters):
            dx = state_cg.atoms[5*i+1].x - state_cg.atoms[5*i].x
            dy = state_cg.atoms[5*i+1].y - state_cg.atoms[5*i].y
            dz = state_cg.atoms[5*i+1].z - state_cg.atoms[5*i].z
            disps_cg.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
        
        converged_cg = "是" if max(disps_cg) < 20 else "否"
        
        print(f"{'':<10} {'CG':<6} {time_cg:<12.1f} {energy_cg:<15.2f} "
              f"{energy_cg/n_waters:<12.2f} {np.mean(disps_cg):<10.2f} {converged_cg:<6}")
        
        # 对比
        speedup = time_scf / time_cg if time_cg > 0 else 0
        energy_diff_pct = abs(energy_cg - energy_scf) / abs(energy_scf) * 100 if energy_scf != 0 else 0
        
        print(f"{'':<10} {'对比':<6} {'加速':<6}{speedup:<6.2f}x "
              f"{'能量差':<9}{energy_diff_pct:<6.1f}%")
        print()
        
        # 保存结果
        results[tol] = {
            'SCF': {'time': time_scf, 'energy': energy_scf, 'disp': np.mean(disps_scf)},
            'CG': {'time': time_cg, 'energy': energy_cg, 'disp': np.mean(disps_cg)}
        }
    
    return results

def main():
    """
    主测试函数
    """
    print("CG vs SCF 最终对比测试")
    print("使用优化的SWM4-NDP水体系")
    print("="*80)
    
    # 测试系统
    system_sizes = [4, 8, 16, 32, 64, 128]
    
    all_results = {}
    
    for n_waters in system_sizes:
        try:
            results = run_comparison(n_waters)
            all_results[n_waters] = results
        except Exception as e:
            print(f"\n错误: 测试 {n_waters} 水失败: {e}")
            import traceback
            traceback.print_exc()
    
    # 总结
    print("\n\n" + "="*80)
    print("性能总结")
    print("="*80)
    
    print(f"\n{'系统':<8} {'容差':<8} {'SCF(ms)':<12} {'CG(ms)':<12} {'加速比':<10} {'能量差(%)':<12}")
    print("-"*70)
    
    for n_waters in sorted(all_results.keys()):
        for tol in [1.0, 10.0, 100.0]:
            if tol in all_results[n_waters]:
                scf = all_results[n_waters][tol]['SCF']
                cg = all_results[n_waters][tol]['CG']
                speedup = scf['time'] / cg['time'] if cg['time'] > 0 else 0
                energy_diff = abs(cg['energy'] - scf['energy']) / abs(scf['energy']) * 100 if scf['energy'] != 0 else 0
                
                print(f"{n_waters:<8} {tol:<8.1f} {scf['time']:<12.1f} {cg['time']:<12.1f} "
                      f"{speedup:<10.2f}x {energy_diff:<12.1f}")
    
    print("\n关键结论:")
    print("1. CG和SCF给出不同的物理图像（能量和Drude位移不同）")
    print("2. CG是线性化近似，SCF是真正的能量最小化")
    print("3. 对于GCMC应用：")
    print("   - 小系统用SCF（更准确）")
    print("   - 大系统可考虑CG（更快但有误差）")
    print("   - 或使用CG初始化+SCF精修的混合策略")

if __name__ == "__main__":
    main()