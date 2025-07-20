#!/usr/bin/env python3
"""
FBP算法最终性能测试
使用合理的水密度和配置
"""

import pygcmc
import numpy as np
import time

def create_dense_water_system(n_waters=128):
    """
    创建密集水系统，接近真实密度
    """
    print(f"\n创建{n_waters}个水分子系统...")
    
    # 真实水密度：1 g/cm³
    # 1个水分子体积 ≈ 30 Å³ = 0.03 nm³
    volume_per_water = 0.03  # nm³
    total_volume = n_waters * volume_per_water
    box_length = np.cbrt(total_volume)
    
    # 稍微增大一点避免严重重叠
    box_length *= 1.1
    
    print(f"盒子尺寸: {box_length:.2f} nm (接近真实密度)")
    
    atoms = []
    residues = []
    
    # 网格放置
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_length / n_per_side
    
    print(f"网格: {n_per_side}³, 间距: {spacing:.3f} nm")
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础位置
                x_base = (i + 0.5) * spacing
                y_base = (j + 0.5) * spacing
                z_base = (k + 0.5) * spacing
                
                # 微小扰动避免完美对称
                x_base += 0.05 * (np.random.rand() - 0.5)
                y_base += 0.05 * (np.random.rand() - 0.5)
                z_base += 0.05 * (np.random.rand() - 0.5)
                
                # SWM4-NDP水分子
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
    
    print(f"截断距离: {state.info.cutoff:.2f} nm")
    
    return state

def test_single_configuration(state, n_waters, algorithm, tolerance, max_iter):
    """测试单个配置"""
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP
    charge = -1.71636
    k_spring = 418400.0
    polarizability = 1.71636**2 * 138.935456 / k_spring
    
    # 添加Drude
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # Thole屏蔽
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = max_iter
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    else:
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 重置Drude
    state_copy = state.copy()
    for i in range(n_waters):
        state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
        state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
        state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
    
    # 多次运行
    times = []
    energy = None
    
    n_runs = 3
    for run in range(n_runs):
        start = time.time()
        energy = force.calculateEnergySCF(state_copy)
        elapsed = (time.time() - start) * 1000
        times.append(elapsed)
    
    avg_time = np.mean(times)
    
    # 位移统计
    disps = []
    for i in range(n_waters):
        dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
        dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
        dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        disps.append(disp)
    
    return {
        'energy': energy,
        'time': avg_time,
        'times': times,
        'avg_disp': np.mean(disps),
        'max_disp': np.max(disps)
    }

def main_benchmark():
    """主要基准测试"""
    print("="*80)
    print("FBP vs SCF 性能基准测试 (128水分子)")
    print("="*80)
    
    n_waters = 128
    state = create_dense_water_system(n_waters)
    
    # 测试配置
    tests = [
        ("SCF高精度", "SCF", 0.01, 1000),
        ("SCF中精度", "SCF", 0.1, 500),
        ("SCF低精度", "SCF", 1.0, 500),
        ("FBP高精度", "FBP", 0.1, 200),
        ("FBP标准", "FBP", 1.0, 200),
        ("FBP快速", "FBP", 10.0, 200),
    ]
    
    results = []
    ref_energy = None
    
    print("\n开始测试...")
    for name, algo, tol, max_iter in tests:
        print(f"\n{name} (tol={tol}):")
        try:
            result = test_single_configuration(state, n_waters, algo, tol, max_iter)
            result['name'] = name
            results.append(result)
            
            if ref_energy is None and algo == "SCF":
                ref_energy = result['energy']
            
            print(f"  能量: {result['energy']:.2f} kJ/mol")
            print(f"  时间: {result['time']:.1f} ms (运行{len(result['times'])}次)")
            print(f"  位移: 平均{result['avg_disp']:.2f} pm, 最大{result['max_disp']:.2f} pm")
            
        except Exception as e:
            print(f"  失败: {e}")
    
    # 结果汇总
    print("\n\n" + "="*100)
    print("最终结果汇总")
    print("="*100)
    print(f"{'算法':<15} {'能量(kJ/mol)':<15} {'误差(kJ/mol)':<12} {'误差(%)':<10} {'时间(ms)':<12} {'速度提升':<10} {'平均位移(pm)':<12}")
    print("-"*100)
    
    # 找SCF中精度作为基准
    scf_ref = next((r for r in results if r['name'] == "SCF中精度"), None)
    if scf_ref:
        ref_time = scf_ref['time']
        ref_energy = scf_ref['energy']
        
        for r in results:
            error = abs(r['energy'] - ref_energy)
            error_pct = error / abs(ref_energy) * 100 if abs(ref_energy) > 1.0 else 0
            speedup = ref_time / r['time']
            
            print(f"{r['name']:<15} {r['energy']:<15.2f} {error:<12.2f} {error_pct:<10.3f} {r['time']:<12.1f} {speedup:<10.2f}x {r['avg_disp']:<12.2f}")
    
    # 核心对比
    print("\n\n核心性能对比:")
    print("-"*60)
    
    scf_mid = next((r for r in results if r['name'] == "SCF中精度"), None)
    fbp_std = next((r for r in results if r['name'] == "FBP标准"), None)
    
    if scf_mid and fbp_std:
        speedup = scf_mid['time'] / fbp_std['time']
        error = abs(fbp_std['energy'] - scf_mid['energy'])
        error_pct = error / abs(scf_mid['energy']) * 100 if abs(scf_mid['energy']) > 1.0 else 0
        
        print(f"SCF中精度 vs FBP标准:")
        print(f"  FBP速度提升: {speedup:.2f}倍")
        print(f"  FBP能量误差: {error:.2f} kJ/mol ({error_pct:.2f}%)")
        print(f"  FBP计算时间: {fbp_std['time']:.1f} ms")
        print(f"  SCF计算时间: {scf_mid['time']:.1f} ms")
        
        if speedup > 1 and error_pct < 1:
            print(f"\n结论: FBP在128水分子系统中提供{speedup:.1f}倍加速，误差仅{error_pct:.2f}%")
        elif speedup < 1:
            print(f"\n注意: FBP在此系统中比SCF慢")

def test_scaling():
    """测试不同大小系统的扩展性"""
    print("\n\n\n扩展性测试")
    print("="*80)
    
    sizes = [32, 64, 128]
    
    for n in sizes:
        print(f"\n{n}水分子系统:")
        state = create_dense_water_system(n)
        
        # 只测试关键算法
        scf = test_single_configuration(state, n, "SCF", 0.1, 500)
        fbp = test_single_configuration(state, n, "FBP", 1.0, 200)
        
        speedup = scf['time'] / fbp['time']
        error = abs(fbp['energy'] - scf['energy'])
        error_pct = error / abs(scf['energy']) * 100 if abs(scf['energy']) > 1.0 else 0
        
        print(f"  SCF: {scf['time']:.1f} ms, 能量 {scf['energy']:.1f} kJ/mol")
        print(f"  FBP: {fbp['time']:.1f} ms, 能量 {fbp['energy']:.1f} kJ/mol")
        print(f"  加速: {speedup:.2f}x, 误差: {error_pct:.2f}%")

if __name__ == "__main__":
    main_benchmark()
    test_scaling()