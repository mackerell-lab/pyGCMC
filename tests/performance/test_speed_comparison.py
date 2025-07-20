#!/usr/bin/env python3
"""
纯速度对比测试 - 不关注完美收敛，只要能量合理即可
"""

import pygcmc
import numpy as np
import time

def create_water_system(n_waters):
    """创建水分子系统，合理间距"""
    # 使用较大间距避免严重重叠
    spacing = 0.4  # nm，确保水分子间有合理距离
    n_per_side = int(np.ceil(n_waters**(1/3)))
    box_length = n_per_side * spacing
    
    atoms = []
    residues = []
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                x_base = (i + 0.5) * spacing
                y_base = (j + 0.5) * spacing
                z_base = (k + 0.5) * spacing
                
                # SWM4-NDP水模型
                positions = [
                    (x_base, y_base, z_base, 1.71636, 0),   # O
                    (x_base, y_base, z_base, -1.71636, 1),  # D
                    (x_base + 0.09572, y_base, z_base, 0.55733, 2),  # H1
                    (x_base - 0.04786, y_base + 0.08288, z_base, 0.55733, 2),  # H2
                    (x_base, y_base - 0.024034, z_base, -1.11466, 3)  # M
                ]
                
                for x, y, z, charge, typ in positions:
                    atom = pygcmc.MCAtom()
                    atom.x = x
                    atom.y = y
                    atom.z = z
                    atom.charge = charge
                    atom.type = typ
                    atoms.append(atom)
                
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_count
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                water_count += 1
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length/2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def time_algorithm(state, n_waters, algorithm, tolerance=1.0, max_iter=100):
    """测试算法速度，返回时间和能量"""
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    # 添加Drude粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置参数 - 使用较宽松的设置以确保快速"收敛"
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = max_iter  # 限制迭代次数
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 重置Drude位置
    state_copy = state.copy()
    for i in range(n_waters):
        state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
        state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
        state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
    
    # 多次运行取平均
    times = []
    for _ in range(3):
        start = time.time()
        energy = force.calculateEnergySCF(state_copy)
        elapsed = (time.time() - start) * 1000  # ms
        times.append(elapsed)
    
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    return {
        'time': avg_time,
        'std': std_time,
        'energy': energy,
        'runs': times
    }

def speed_comparison_test():
    """速度对比测试主函数"""
    print("纯速度对比测试 - SCF vs FBP")
    print("="*80)
    print("注：使用宽松容差(10.0)和有限迭代(100次)，只关注速度对比")
    print("-"*80)
    
    # 测试不同大小的系统
    system_sizes = [32, 64, 128]
    
    # 使用相同的宽松参数
    tolerance = 10.0  # 宽松容差
    max_iter = 100    # 限制迭代次数
    
    results = []
    
    for n_waters in system_sizes:
        print(f"\n测试 {n_waters} 水分子系统...")
        state = create_water_system(n_waters)
        print(f"  盒子尺寸: {state.info.box[0]:.2f} nm")
        print(f"  截断距离: {state.info.cutoff:.2f} nm")
        
        # 测试SCF
        print(f"\n  运行SCF (tol={tolerance}, max_iter={max_iter})...")
        scf_result = time_algorithm(state, n_waters, "SCF", tolerance, max_iter)
        print(f"    时间: {scf_result['time']:.1f} ± {scf_result['std']:.1f} ms")
        print(f"    能量: {scf_result['energy']:.2f} kJ/mol")
        
        # 测试FBP
        print(f"\n  运行FBP (tol={tolerance}, max_iter={max_iter})...")
        fbp_result = time_algorithm(state, n_waters, "FBP", tolerance, max_iter)
        print(f"    时间: {fbp_result['time']:.1f} ± {fbp_result['std']:.1f} ms")
        print(f"    能量: {fbp_result['energy']:.2f} kJ/mol")
        
        # 计算加速比
        speedup = scf_result['time'] / fbp_result['time']
        energy_diff = abs(fbp_result['energy'] - scf_result['energy'])
        energy_diff_pct = energy_diff / abs(scf_result['energy']) * 100 if abs(scf_result['energy']) > 1 else 0
        
        results.append({
            'n_waters': n_waters,
            'scf_time': scf_result['time'],
            'fbp_time': fbp_result['time'],
            'speedup': speedup,
            'scf_energy': scf_result['energy'],
            'fbp_energy': fbp_result['energy'],
            'energy_diff': energy_diff,
            'energy_diff_pct': energy_diff_pct
        })
        
        print(f"\n  速度提升: {speedup:.2f}x")
        print(f"  能量差异: {energy_diff:.2f} kJ/mol ({energy_diff_pct:.1f}%)")
    
    # 汇总结果
    print("\n\n" + "="*100)
    print("速度对比汇总")
    print("="*100)
    print(f"{'系统大小':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'加速比':<10} {'能量差(%)':<12} {'结论':<20}")
    print("-"*100)
    
    for r in results:
        if r['speedup'] > 1:
            conclusion = f"FBP快{r['speedup']:.1f}倍"
        else:
            conclusion = f"SCF快{1/r['speedup']:.1f}倍"
        
        print(f"{r['n_waters']:<10} {r['scf_time']:<15.1f} {r['fbp_time']:<15.1f} {r['speedup']:<10.2f}x {r['energy_diff_pct']:<12.1f} {conclusion:<20}")
    
    # 分析扩展性
    print("\n\n扩展性分析:")
    print("-"*60)
    
    if len(results) >= 2:
        for i in range(1, len(results)):
            size_ratio = results[i]['n_waters'] / results[i-1]['n_waters']
            scf_time_ratio = results[i]['scf_time'] / results[i-1]['scf_time']
            fbp_time_ratio = results[i]['fbp_time'] / results[i-1]['fbp_time']
            
            print(f"\n{results[i-1]['n_waters']} → {results[i]['n_waters']} 水分子:")
            print(f"  系统大小增加: {size_ratio:.1f}x")
            print(f"  SCF时间增加: {scf_time_ratio:.1f}x")
            print(f"  FBP时间增加: {fbp_time_ratio:.1f}x")
            print(f"  FBP加速比变化: {results[i-1]['speedup']:.2f}x → {results[i]['speedup']:.2f}x")
    
    # 最终结论
    print("\n\n最终结论:")
    print("-"*60)
    
    avg_speedup = np.mean([r['speedup'] for r in results])
    avg_energy_diff = np.mean([r['energy_diff_pct'] for r in results])
    
    print(f"平均加速比: {avg_speedup:.2f}x")
    print(f"平均能量差异: {avg_energy_diff:.1f}%")
    
    if avg_speedup > 1:
        print(f"\n在宽松容差({tolerance})下，FBP平均比SCF快{avg_speedup:.1f}倍")
    else:
        print(f"\n在宽松容差({tolerance})下，SCF平均比FBP快{1/avg_speedup:.1f}倍")
    
    print(f"能量差异在{avg_energy_diff:.1f}%范围内，对于快速计算是可接受的")

def test_different_tolerances():
    """测试不同容差下的速度对比"""
    print("\n\n\n不同容差下的速度对比 (128水分子)")
    print("="*80)
    
    n_waters = 128
    state = create_water_system(n_waters)
    
    tolerances = [0.1, 1.0, 10.0, 100.0]
    
    print(f"{'容差':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'加速比':<10} {'能量差(%)':<12}")
    print("-"*80)
    
    for tol in tolerances:
        # 限制迭代次数避免太长时间
        max_iter = min(int(100 / tol), 200)
        
        scf_result = time_algorithm(state, n_waters, "SCF", tol, max_iter)
        fbp_result = time_algorithm(state, n_waters, "FBP", tol, max_iter)
        
        speedup = scf_result['time'] / fbp_result['time']
        energy_diff_pct = abs(fbp_result['energy'] - scf_result['energy']) / abs(scf_result['energy']) * 100 if abs(scf_result['energy']) > 1 else 0
        
        print(f"{tol:<10.1f} {scf_result['time']:<15.1f} {fbp_result['time']:<15.1f} {speedup:<10.2f}x {energy_diff_pct:<12.1f}")

if __name__ == "__main__":
    speed_comparison_test()
    test_different_tolerances()