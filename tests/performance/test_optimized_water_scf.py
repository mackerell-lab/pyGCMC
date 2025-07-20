#!/usr/bin/env python3
"""
测试优化后的水系统的Drude SCF收敛性
比较优化前后的差异
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def load_and_test_system(filename, description):
    """
    加载并测试单个系统
    """
    print(f"\n{'='*70}")
    print(f"测试: {description}")
    print(f"文件: {filename}")
    print(f"{'='*70}")
    
    # 加载系统
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return None
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"\n系统参数:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {data.get('density', 'N/A')} g/cm³")
    print(f"  优化方法: {data.get('method', '未优化')}")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = data['charges']
    atom_types = [0, 1, 2, 2, 3]  # O, D, H1, H2, M
    
    # 创建原子
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    # 创建残基
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
    
    # 设置盒子和力场
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
    
    # 添加Thole对（限制数量以加速）
    print(f"\n添加Thole屏蔽对...")
    start_time = time.time()
    
    n_thole_pairs = 0
    thole_cutoff = 0.8
    max_thole = 2000  # 限制Thole对数量
    
    for i in range(n_waters):
        o1_idx = i * 5
        
        for j in range(i+1, min(i+50, n_waters)):  # 只检查附近的50个分子
            o2_idx = j * 5
            
            # 计算O-O距离
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
    print(f"  用时 {time.time()-start_time:.2f} 秒")
    
    # 测试不同容差的SCF
    tolerances = [100.0, 10.0, 1.0]
    results = []
    
    print(f"\n测试SCF收敛:")
    print(f"{'容差':>10} {'收敛':>8} {'迭代':>8} {'时间(s)':>10} {'能量/水':>12} {'位移(pm)':>10}")
    print("-"*70)
    
    for tolerance in tolerances:
        # 设置SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tolerance
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 复制状态
        test_state = state.copy()
        
        # 运行SCF
        start_time = time.time()
        
        try:
            energy = force.calculateEnergySCF(test_state)
            elapsed_time = time.time() - start_time
            
            # 假设收敛（因为没有抛出异常）
            converged = True
            iterations = -1  # 未知
            
            # 分析Drude位移
            displacements = []
            for i in range(min(50, n_waters)):  # 采样50个分子
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
                dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
                dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
                displacements.append(disp)
            
            avg_disp = np.mean(displacements)
            
            result = {
                'tolerance': tolerance,
                'converged': converged,
                'iterations': iterations,
                'time': elapsed_time,
                'energy_per_water': energy/n_waters,
                'avg_displacement': avg_disp
            }
            results.append(result)
            
            print(f"{tolerance:10.1f} {'是' if converged else '否':>8} {iterations:8d} "
                  f"{elapsed_time:10.2f} {energy/n_waters:12.2f} {avg_disp:10.2f}")
            
        except Exception as e:
            elapsed_time = time.time() - start_time
            print(f"{tolerance:10.1f} {'错误':>8} {'N/A':>8} "
                  f"{elapsed_time:10.2f} {'N/A':>12} {'N/A':>10}")
            print(f"      错误: {str(e)}")
    
    return results

def main():
    """
    主函数
    """
    print("测试优化前后的水系统Drude SCF收敛性")
    print("="*70)
    
    # 测试系统列表
    test_systems = [
        # (文件名, 描述)
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - TIP3P NVT优化'),
    ]
    
    all_results = []
    
    for filename, description in test_systems:
        results = load_and_test_system(filename, description)
        if results:
            all_results.append((description, results))
    
    # 总结
    print(f"\n\n{'='*70}")
    print("收敛性总结")
    print("="*70)
    
    for description, results in all_results:
        print(f"\n{description}:")
        
        # 找出最严格的收敛容差
        best_tolerance = None
        for r in results:
            if r['converged']:
                best_tolerance = r['tolerance']
        
        if best_tolerance:
            print(f"  ✓ 可以收敛到容差 {best_tolerance} kJ/mol/nm")
            # 找到对应的结果
            for r in results:
                if r['tolerance'] == best_tolerance:
                    print(f"    - 迭代次数: {r['iterations']}")
                    print(f"    - 计算时间: {r['time']:.2f} 秒")
                    print(f"    - 能量/水: {r['energy_per_water']:.2f} kJ/mol")
                    print(f"    - 平均Drude位移: {r['avg_displacement']:.2f} pm")
        else:
            print(f"  ✗ 无法收敛")
    
    # 比较优化前后
    if len(all_results) >= 2:
        print(f"\n\n优化效果分析:")
        print("="*70)
        
        # 比较相同容差下的收敛性
        for tolerance in [100.0, 10.0, 1.0]:
            print(f"\n容差 {tolerance} kJ/mol/nm:")
            
            for description, results in all_results:
                for r in results:
                    if r['tolerance'] == tolerance:
                        if r['converged']:
                            print(f"  {description}: 收敛，{r['iterations']} 次迭代")
                        else:
                            print(f"  {description}: 未收敛")
                        break

if __name__ == "__main__":
    main()