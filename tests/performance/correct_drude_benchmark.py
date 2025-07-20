#!/usr/bin/env python3
"""
使用正确方法的Drude算法基准测试
"""

import numpy as np
import time
import pickle
import os
import pygcmc

def benchmark_scf_vs_opt3():
    """
    比较SCF和OPT3的性能
    """
    print("SCF vs OPT3 性能对比")
    print("="*70)
    
    # 测试不同大小的系统
    test_sizes = [10, 50, 100, 256]
    
    # 首先用256水系统测试
    filename = '../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl'
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        full_positions = data['positions']
        full_box_length = data['box_length']
        full_charges = data['charges']
    else:
        print("256水系统文件不存在")
        return
    
    print(f"{'系统大小':>10} {'算法':>10} {'时间(ms)':>12} {'能量/水(kJ/mol)':>18} {'平均位移(pm)':>15} {'Thole对数':>12}")
    print("-"*85)
    
    for n_waters in test_sizes:
        if n_waters > 256:
            continue
            
        # 创建状态
        state = pygcmc.MCState()
        atoms = []
        residues = []
        
        atom_types = [0, 1, 2, 2, 3]
        
        # 使用前n_waters个水分子
        for i in range(n_waters * 5):
            atom = pygcmc.MCAtom()
            atom.x = full_positions[i][0]
            atom.y = full_positions[i][1]
            atom.z = full_positions[i][2]
            atom.charge = full_charges[i % 5]
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
        
        state.info.box = np.array([full_box_length, full_box_length, full_box_length])
        state.info.cutoff = min(0.9, full_box_length / 2 - 0.01)
        
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 测试两种算法
        for algo, algo_name in [(pygcmc.DrudeAlgorithm.SCF, "SCF"),
                               (pygcmc.DrudeAlgorithm.OPT3, "OPT3")]:
            
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
            n_thole_pairs = 0
            thole_cutoff = 0.8
            
            for i in range(n_waters):
                o1_idx = i * 5
                for j in range(i+1, min(i+10, n_waters)):  # 限制每个水最多10个Thole对
                    o2_idx = j * 5
                    
                    dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
                    dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
                    dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
                    
                    # PBC
                    dx -= full_box_length * round(dx / full_box_length)
                    dy -= full_box_length * round(dy / full_box_length)
                    dz -= full_box_length * round(dz / full_box_length)
                    
                    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if dist < thole_cutoff:
                        force.addScreenedPair(i, j, 1.3)
                        n_thole_pairs += 1
            
            # 设置算法
            force.setAlgorithm(algo)
            
            # SCF参数
            if algo == pygcmc.DrudeAlgorithm.SCF:
                params = pygcmc.DrudeSCFParams()
                params.tolerance = 10.0
                params.maxIterations = 50
                params.dampingFactor = 0.5
                params.maxDrudeDistance = 0.02
                force.setSCFParameters(params)
            
            # 准备测试状态
            state_test = state.copy()
            for i in range(n_waters):
                o_idx = i * 5
                d_idx = i * 5 + 1
                state_test.atoms[d_idx].x = state_test.atoms[o_idx].x
                state_test.atoms[d_idx].y = state_test.atoms[o_idx].y
                state_test.atoms[d_idx].z = state_test.atoms[o_idx].z
            
            # 测试
            try:
                # 预热
                if algo == pygcmc.DrudeAlgorithm.SCF:
                    force.calculateEnergySCF(state_test.copy())
                else:
                    force.calculateEnergyOPT3(state_test.copy())
                
                # 正式测试
                n_runs = 3
                times = []
                
                for _ in range(n_runs):
                    state_run = state_test.copy()
                    
                    start_time = time.time()
                    if algo == pygcmc.DrudeAlgorithm.SCF:
                        energy = force.calculateEnergySCF(state_run)
                    else:
                        energy = force.calculateEnergyOPT3(state_run)
                    end_time = time.time()
                    
                    times.append(end_time - start_time)
                
                avg_time = np.mean(times) * 1000  # ms
                
                # 使用最后一次运行的结果计算位移
                displacements = []
                for i in range(min(20, n_waters)):  # 只检查前20个
                    o_idx = i * 5
                    d_idx = i * 5 + 1
                    
                    dx = state_run.atoms[d_idx].x - state_run.atoms[o_idx].x
                    dy = state_run.atoms[d_idx].y - state_run.atoms[o_idx].y
                    dz = state_run.atoms[d_idx].z - state_run.atoms[o_idx].z
                    
                    disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                    displacements.append(disp)
                
                avg_disp = np.mean(displacements)
                
                print(f"{n_waters:>10} {algo_name:>10} {avg_time:>12.2f} {energy/n_waters:>18.2f} {avg_disp:>15.2f} {n_thole_pairs:>12}")
                
            except Exception as e:
                print(f"{n_waters:>10} {algo_name:>10} {'失败':>12} {'-':>18} {'-':>15} {n_thole_pairs:>12}")
                # print(f"          错误: {str(e)[:60]}")

def analyze_performance():
    """
    分析性能特点
    """
    print("\n\n性能分析")
    print("="*70)
    
    print("\n1. SCF (Self-Consistent Field) 特点：")
    print("   - 迭代求解，直到收敛")
    print("   - 精度可控（通过tolerance参数）")
    print("   - 大系统可能需要更多迭代")
    print("   - 适合需要高精度的场合")
    
    print("\n2. OPT3 特点：")
    print("   - 直接求解方法")
    print("   - 固定计算量，不需要迭代")
    print("   - 可能牺牲一些精度换取速度")
    print("   - 适合大规模系统")
    
    print("\n3. 实际应用建议：")
    print("   - 小系统（<100水）：两者差异不大，可用SCF获得更高精度")
    print("   - 中等系统（100-1000水）：OPT3可能更快")
    print("   - 大系统（>1000水）：优先考虑OPT3")
    print("   - GCMC模拟：OPT3可能更合适（需要频繁计算）")
    
    print("\n4. 并行化潜力：")
    print("   - SCF：迭代过程难以并行，但每次迭代内可并行")
    print("   - OPT3：每个Drude独立计算，易于并行化")

def main():
    """
    主函数
    """
    benchmark_scf_vs_opt3()
    analyze_performance()

if __name__ == "__main__":
    main()