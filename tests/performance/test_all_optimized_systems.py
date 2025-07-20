#!/usr/bin/env python3
"""
测试所有优化系统的Drude SCF收敛性
"""

import pygcmc
import numpy as np
import pickle
import time
import os

def test_system_scf(filename, description, max_thole_pairs=None):
    """
    测试单个系统的SCF收敛性
    """
    print(f"\n{'='*70}")
    print(f"{description}")
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
    
    print(f"系统参数:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {data.get('density', 'N/A')} g/cm³")
    print(f"  优化方法: {data.get('method', '未优化')}")
    print(f"  Thole(0.8)/半盒子: {0.8/(box_length/2):.3f}")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = data['charges']
    atom_types = [0, 1, 2, 2, 3]
    
    # 创建原子和残基
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
    
    # 添加Thole对
    print(f"\n添加Thole屏蔽对...")
    start_time = time.time()
    
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 根据系统大小设置Thole对限制
    if max_thole_pairs is None:
        if n_waters <= 256:
            max_thole_pairs = 10000
        elif n_waters <= 512:
            max_thole_pairs = 20000
        elif n_waters <= 1024:
            max_thole_pairs = 30000
        else:
            max_thole_pairs = 40000
    
    # 采样方式添加Thole对
    skip_factor = max(1, n_waters // 512)
    check_range = min(100, n_waters // skip_factor)
    
    for i in range(0, n_waters, skip_factor):
        o1_idx = i * 5
        
        for j in range(i+1, min(i+check_range, n_waters)):
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
                
                if n_thole_pairs >= max_thole_pairs:
                    break
        
        if n_thole_pairs >= max_thole_pairs:
            break
    
    print(f"  添加了 {n_thole_pairs} 个Thole对 (限制: {max_thole_pairs})")
    print(f"  用时 {time.time()-start_time:.2f} 秒")
    
    # 测试SCF收敛
    print(f"\n测试SCF收敛:")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # kJ/mol/nm
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    start_time = time.time()
    
    try:
        energy = force.calculateEnergySCF(state)
        elapsed_time = time.time() - start_time
        
        # 分析Drude位移
        displacements = []
        sample_size = min(100, n_waters)
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
        
        print(f"  ✓ SCF收敛!")
        print(f"  时间: {elapsed_time:.2f} 秒")
        print(f"  总能量: {energy:.2f} kJ/mol")
        print(f"  能量/水: {energy/n_waters:.2f} kJ/mol")
        print(f"  平均Drude位移: {avg_disp:.2f} pm")
        print(f"  最大Drude位移: {max_disp:.2f} pm")
        
        return {
            'n_waters': n_waters,
            'converged': True,
            'time': elapsed_time,
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
    主函数
    """
    print("测试所有优化系统的Drude SCF收敛性")
    print("="*70)
    
    # 测试系统列表
    test_systems = [
        # 未优化系统
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/large_water_systems/water_512.pkl', '512水 - 未优化'),
        ('../tests/performance/large_water_systems/water_1024.pkl', '1024水 - 未优化'),
        
        # NVT优化系统
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
        ('../tests/performance/optimized_water_systems/water_512_nvt_simple.pkl', '512水 - NVT优化'),
        ('../tests/performance/optimized_water_systems/water_1024_nvt_simple.pkl', '1024水 - NVT优化'),
        ('../tests/performance/optimized_water_systems/water_2048_nvt_simple.pkl', '2048水 - NVT优化'),
        ('../tests/performance/optimized_water_systems/water_4096_nvt_simple.pkl', '4096水 - NVT优化'),
    ]
    
    results = []
    
    for filename, description in test_systems:
        result = test_system_scf(filename, description)
        if result:
            results.append((description, result))
    
    # 汇总结果
    print(f"\n\n{'='*70}")
    print("测试结果汇总")
    print("="*70)
    
    print(f"\n{'系统':^20} {'收敛':^6} {'时间(s)':^10} {'能量/水':^12} {'位移(pm)':^10}")
    print("-"*60)
    
    for description, result in results:
        if result['converged']:
            print(f"{description:20} {'是':^6} {result['time']:^10.2f} "
                  f"{result['energy_per_water']:^12.2f} {result['avg_displacement']:^10.2f}")
        else:
            print(f"{description:20} {'否':^6} {result['time']:^10.2f} "
                  f"{'N/A':^12} {'N/A':^10}")
    
    # 分析优化效果
    print(f"\n\n优化效果分析:")
    print("="*70)
    
    # 比较相同大小系统的优化前后
    sizes = [256, 512, 1024]
    for size in sizes:
        unopt = None
        opt = None
        
        for desc, res in results:
            if f"{size}水 - 未优化" in desc:
                unopt = res
            elif f"{size}水 - NVT优化" in desc:
                opt = res
        
        if unopt and opt:
            print(f"\n{size}水系统:")
            print(f"  未优化: {'收敛' if unopt['converged'] else '未收敛'}")
            print(f"  NVT优化: {'收敛' if opt['converged'] else '未收敛'}")
            
            if unopt['converged'] and opt['converged']:
                time_improve = (unopt['time'] - opt['time']) / unopt['time'] * 100
                print(f"  时间改善: {time_improve:.1f}%")
    
    # 大系统性能
    print(f"\n\n大系统性能:")
    print("="*70)
    
    large_systems = [opt for opt, res in results if res['n_waters'] >= 2048 and res['converged']]
    if large_systems:
        for desc, res in results:
            if res['n_waters'] >= 2048 and res['converged']:
                print(f"{desc}: {res['time']:.2f}秒, {res['energy_per_water']:.2f} kJ/mol/水")

if __name__ == "__main__":
    main()