#!/usr/bin/env python3
"""
测试大型水系统的Drude SCF收敛性和性能
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def load_water_system(n_waters):
    """
    加载预生成的水系统
    """
    pickle_file = f'../tests/performance/large_water_systems/water_{n_waters}.pkl'
    if not os.path.exists(pickle_file):
        raise FileNotFoundError(f"找不到文件: {pickle_file}")
    
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    # 转换为pygcmc格式
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = data['charges']  # [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]  # O, D, H1, H2, M
    
    positions = data['positions']
    n_waters = data['n_waters']
    box_length = data['box_length']
    
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
        res.type = 0  # 水分子类型
        residues.append(res)
    
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

def test_scf_performance(n_waters, max_thole_pairs=None):
    """
    测试SCF性能
    max_thole_pairs: 限制Thole对数量（用于大系统）
    """
    print(f"\n{'='*70}")
    print(f"测试 {n_waters} 水分子系统")
    print(f"{'='*70}")
    
    # 加载系统
    state, data = load_water_system(n_waters)
    
    print(f"系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  原子总数: {state.activeAtomCount}")
    print(f"  盒子长度: {data['box_length']:.3f} nm")
    print(f"  密度: {data['density']:.6f} g/cm³")
    print(f"  截断距离: {state.info.cutoff:.3f} nm")
    print(f"  Thole(0.8)/半盒子: {0.8/(data['box_length']/2):.3f}")
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP Drude参数
    drude_charge = -1.71636
    polarizability = 0.0009782237  # nm³
    
    # 添加Drude粒子
    print(f"\n添加Drude粒子...")
    start = time.time()
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1,    # D粒子索引
            parentIndex=5*i,     # O原子索引
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=drude_charge,
            polarizability=polarizability,
            aniso12=1.0,
            aniso34=1.0
        )
    print(f"  完成，用时 {time.time()-start:.2f} 秒")
    
    # 添加Thole屏蔽对
    print(f"\n添加Thole屏蔽对...")
    start = time.time()
    n_thole_pairs = 0
    thole_cutoff = 0.8  # nm
    thole_param = 1.3
    
    # 对于大系统，可以限制Thole对数量
    if max_thole_pairs and n_waters > 512:
        print(f"  注意：限制Thole对数量为 {max_thole_pairs}")
    
    # 使用简单的距离筛选
    # 为了效率，只检查一部分对
    skip_factor = 1
    if n_waters >= 2048:
        skip_factor = 2  # 对于大系统，跳过一些检查
        print(f"  使用跳跃因子 {skip_factor} 加速")
    
    checked_pairs = 0
    for i in range(0, n_waters, skip_factor):
        o1_idx = i * 5
        
        # 只检查附近的分子
        check_range = min(100, n_waters - i - 1)  # 限制检查范围
        
        for j in range(i+1, min(i+1+check_range, n_waters)):
            o2_idx = j * 5
            
            # 计算O-O距离
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            # 应用PBC
            box = data['box_length']
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, thole_param)
                n_thole_pairs += 1
                
                if max_thole_pairs and n_thole_pairs >= max_thole_pairs:
                    break
            
            checked_pairs += 1
            
        if max_thole_pairs and n_thole_pairs >= max_thole_pairs:
            print(f"  达到最大Thole对限制")
            break
    
    print(f"  检查了 {checked_pairs} 对")
    print(f"  添加了 {n_thole_pairs} 个Thole屏蔽对")
    print(f"  完成，用时 {time.time()-start:.2f} 秒")
    
    # 测试SCF收敛
    print(f"\n测试SCF收敛性:")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # kJ/mol/nm
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02  # 2 pm硬墙约束
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    print(f"  运行SCF (容差={params.tolerance} kJ/mol/nm)...")
    start_time = time.time()
    
    try:
        energy = force.calculateEnergySCF(state)
        elapsed_time = time.time() - start_time
        
        print(f"  ✓ SCF收敛!")
        print(f"  时间: {elapsed_time:.2f} 秒")
        print(f"  总能量: {energy:.2f} kJ/mol")
        print(f"  能量/水: {energy/n_waters:.2f} kJ/mol")
        
        # 分析Drude位移
        displacements = []
        sample_size = min(100, n_waters)  # 采样分析
        for i in range(0, n_waters, max(1, n_waters//sample_size)):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state.atoms[d_idx].x - state.atoms[o_idx].x
            dy = state.atoms[d_idx].y - state.atoms[o_idx].y
            dz = state.atoms[d_idx].z - state.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            displacements.append(disp)
        
        avg_disp = np.mean(displacements)
        max_disp = np.max(displacements)
        
        print(f"  平均Drude位移: {avg_disp:.2f} pm")
        print(f"  最大Drude位移: {max_disp:.2f} pm")
        
        return {
            'n_waters': n_waters,
            'converged': True,
            'time': elapsed_time,
            'energy': energy,
            'energy_per_water': energy/n_waters,
            'avg_displacement': avg_disp,
            'n_thole_pairs': n_thole_pairs
        }
        
    except Exception as e:
        elapsed_time = time.time() - start_time
        print(f"  ✗ SCF未收敛")
        print(f"  时间: {elapsed_time:.2f} 秒")
        print(f"  错误: {str(e)}")
        
        return {
            'n_waters': n_waters,
            'converged': False,
            'time': elapsed_time,
            'error': str(e)
        }

def main():
    """
    主测试函数
    """
    print("大型水系统Drude SCF测试")
    print("="*70)
    
    # 测试系统
    test_systems = [
        (256, None),      # 完整Thole
        (512, None),      # 完整Thole
        (1024, 5000),     # 限制Thole对
        (2048, 5000),     # 限制Thole对
        (4096, 5000)      # 限制Thole对
    ]
    
    results = []
    
    for n_waters, max_thole in test_systems:
        try:
            result = test_scf_performance(n_waters, max_thole)
            results.append(result)
        except Exception as e:
            print(f"\n错误测试 {n_waters} 水: {e}")
            import traceback
            traceback.print_exc()
    
    # 汇总结果
    print(f"\n\n{'='*70}")
    print("测试结果汇总")
    print("="*70)
    
    print(f"\n{'系统大小':>8} {'收敛':>6} {'时间(s)':>10} {'能量/水':>12} {'位移(pm)':>10} {'Thole对':>10}")
    print("-"*70)
    
    for result in results:
        if result.get('converged'):
            print(f"{result['n_waters']:8d} {'是':>6} {result['time']:10.2f} "
                  f"{result['energy_per_water']:12.2f} {result['avg_displacement']:10.2f} "
                  f"{result['n_thole_pairs']:10d}")
        else:
            print(f"{result['n_waters']:8d} {'否':>6} {result['time']:10.2f} "
                  f"{'N/A':>12} {'N/A':>10} {'N/A':>10}")
    
    # 性能分析
    if len([r for r in results if r.get('converged')]) > 1:
        print(f"\n\n性能分析:")
        sizes = []
        times = []
        
        for r in results:
            if r.get('converged'):
                sizes.append(r['n_waters'])
                times.append(r['time'])
        
        if len(sizes) > 1:
            # 简单的复杂度估算
            log_sizes = np.log(sizes)
            log_times = np.log(times)
            
            # 线性回归估算 time ~ n^k
            k = np.polyfit(log_sizes, log_times, 1)[0]
            print(f"  时间复杂度: O(n^{k:.2f})")

if __name__ == "__main__":
    main()