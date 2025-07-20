#!/usr/bin/env python3
"""
测试FBP在大体系（128水分子）中的性能和精度 - 改进版
使用更合理的初始配置避免原子重叠
"""

import pygcmc
import numpy as np
import time

def create_water_box_improved(n_waters=128):
    """
    创建水分子盒子，确保合理的初始间距
    """
    print(f"创建{n_waters}个水分子系统...")
    
    # 使用更大的初始间距避免重叠
    min_spacing = 0.35  # nm，水分子间最小距离
    
    # 计算需要的盒子大小
    n_per_side = int(np.ceil(n_waters**(1/3)))
    box_length = n_per_side * min_spacing * 1.2  # 留出一些空间
    
    print(f"盒子尺寸: {box_length:.2f} nm")
    print(f"网格: {n_per_side}x{n_per_side}x{n_per_side}")
    
    atoms = []
    residues = []
    
    # 在立方体网格中放置水分子
    spacing = box_length / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 网格中心位置
                x_base = (i + 0.5) * spacing
                y_base = (j + 0.5) * spacing
                z_base = (k + 0.5) * spacing
                
                # SWM4-NDP水模型原子位置
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
    
    # 设置周期性边界条件
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(box_length/2 - 0.1, 1.2)  # 典型的截断距离是1.2 nm
    
    # 力场参数（SWM4-NDP）
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]  # 只有O有LJ
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print(f"系统创建完成: {n_waters}个水分子, {len(atoms)}个原子")
    print(f"截断距离: {state.info.cutoff:.2f} nm")
    
    return state

def test_single_algorithm(state, n_waters, algorithm, tolerance, description=""):
    """测试单个算法"""
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0  # kJ/mol/nm²
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
    
    # 添加Thole屏蔽对
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = 500 if algorithm == "SCF" else 100  # FBP通常收敛更快
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 复制状态并重置Drude位置
    state_copy = state.copy()
    for i in range(n_waters):
        state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
        state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
        state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
    
    # 多次运行取平均时间
    n_runs = 3
    times = []
    energy = None
    
    print(f"  运行{n_runs}次...")
    for run in range(n_runs):
        start_time = time.time()
        energy = force.calculateEnergySCF(state_copy)
        end_time = time.time()
        times.append((end_time - start_time) * 1000)  # ms
        print(f"    第{run+1}次: {times[-1]:.1f} ms")
    
    avg_time = np.mean(times)
    
    # 计算Drude位移统计
    displacements = []
    for i in range(n_waters):
        dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
        dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
        dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # nm to pm
        displacements.append(disp)
    
    return {
        'energy': energy,
        'time': avg_time,
        'avg_displacement': np.mean(displacements),
        'max_displacement': np.max(displacements),
        'std_displacement': np.std(displacements)
    }

def test_128_water_system():
    """专门测试128水分子系统"""
    print("FBP在128水分子系统中的性能测试")
    print("="*80)
    
    # 创建系统
    n_waters = 128
    state = create_water_box_improved(n_waters)
    
    print("\n测试算法性能...")
    
    # 先运行高精度SCF作为参考
    print("\n1. SCF高精度 (tol=0.01) - 参考标准")
    scf_ref = test_single_algorithm(state, n_waters, "SCF", 0.01, "参考")
    
    # 运行中精度SCF
    print("\n2. SCF中精度 (tol=0.1)")
    scf_mid = test_single_algorithm(state, n_waters, "SCF", 0.1, "中精度")
    
    # 运行FBP标准
    print("\n3. FBP标准 (tol=1.0)")
    fbp_std = test_single_algorithm(state, n_waters, "FBP", 1.0, "标准")
    
    # 运行FBP高精度
    print("\n4. FBP高精度 (tol=0.1)")
    fbp_high = test_single_algorithm(state, n_waters, "FBP", 0.1, "高精度")
    
    # 结果汇总
    print("\n\n结果汇总 (128水分子系统):")
    print("="*90)
    print(f"{'算法':<20} {'能量(kJ/mol)':<15} {'能量误差':<12} {'误差(%)':<10} {'时间(ms)':<12} {'速度提升':<10}")
    print("-"*90)
    
    ref_energy = scf_ref['energy']
    ref_time = scf_ref['time']
    
    results = [
        ("SCF高精度(0.01)", scf_ref),
        ("SCF中精度(0.1)", scf_mid),
        ("FBP标准(1.0)", fbp_std),
        ("FBP高精度(0.1)", fbp_high)
    ]
    
    for name, result in results:
        energy_error = abs(result['energy'] - ref_energy)
        error_percent = (energy_error / abs(ref_energy) * 100) if abs(ref_energy) > 0.01 else 0
        speedup = ref_time / result['time']
        
        print(f"{name:<20} {result['energy']:<15.6f} {energy_error:<12.6f} {error_percent:<10.3f} {result['time']:<12.1f} {speedup:<10.1f}x")
    
    # Drude位移分析
    print("\n\nDrude位移分析:")
    print("-"*70)
    print(f"{'算法':<20} {'平均(pm)':<12} {'最大(pm)':<12} {'标准差(pm)':<12}")
    print("-"*70)
    
    for name, result in results:
        print(f"{name:<20} {result['avg_displacement']:<12.3f} {result['max_displacement']:<12.3f} {result['std_displacement']:<12.3f}")
    
    # 性能分析
    print("\n\n性能分析:")
    print("-"*60)
    
    fbp_speedup_std = ref_time / fbp_std['time']
    fbp_speedup_high = ref_time / fbp_high['time']
    fbp_error_std = (abs(fbp_std['energy'] - ref_energy) / abs(ref_energy) * 100) if abs(ref_energy) > 0.01 else 0
    fbp_error_high = (abs(fbp_high['energy'] - ref_energy) / abs(ref_energy) * 100) if abs(ref_energy) > 0.01 else 0
    
    print(f"FBP标准模式:")
    print(f"  - 速度提升: {fbp_speedup_std:.1f}x")
    print(f"  - 能量误差: {fbp_error_std:.3f}%")
    print(f"  - 效率指标: {fbp_speedup_std/max(fbp_error_std, 0.001):.1f} (速度/误差)")
    
    print(f"\nFBP高精度模式:")
    print(f"  - 速度提升: {fbp_speedup_high:.1f}x")
    print(f"  - 能量误差: {fbp_error_high:.3f}%")
    print(f"  - 效率指标: {fbp_speedup_high/max(fbp_error_high, 0.001):.1f} (速度/误差)")
    
    # 绝对能量分析
    print("\n\n绝对能量分析:")
    print("-"*60)
    if ref_energy > 0:
        print("注意：正能量表示系统处于排斥状态")
        print("这可能是由于初始配置中分子间距较小")
    else:
        print("负能量表示系统处于吸引状态")

def test_scaling_behavior():
    """测试算法的扩展行为"""
    print("\n\n\n算法扩展性测试")
    print("="*80)
    
    sizes = [32, 64, 128]
    scf_times = []
    fbp_times = []
    scf_energies = []
    fbp_energies = []
    
    for n in sizes:
        print(f"\n测试 {n} 水分子系统...")
        state = create_water_box_improved(n)
        
        # SCF
        print("  运行SCF...")
        scf_result = test_single_algorithm(state, n, "SCF", 0.1, "")
        scf_times.append(scf_result['time'])
        scf_energies.append(scf_result['energy'])
        
        # FBP
        print("  运行FBP...")
        fbp_result = test_single_algorithm(state, n, "FBP", 1.0, "")
        fbp_times.append(fbp_result['time'])
        fbp_energies.append(fbp_result['energy'])
    
    print("\n\n扩展性结果:")
    print("-"*80)
    print(f"{'系统大小':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'加速比':<10} {'FBP误差(%)':<12}")
    print("-"*80)
    
    for i, n in enumerate(sizes):
        speedup = scf_times[i] / fbp_times[i]
        error = abs(fbp_energies[i] - scf_energies[i]) / abs(scf_energies[i]) * 100 if abs(scf_energies[i]) > 0.01 else 0
        print(f"{n:<10} {scf_times[i]:<15.1f} {fbp_times[i]:<15.1f} {speedup:<10.1f}x {error:<12.3f}")
    
    # 分析时间复杂度
    print("\n\n时间复杂度分析:")
    print("-"*60)
    print("理论：O(N²) 由于所有对相互作用")
    
    for i in range(1, len(sizes)):
        size_ratio = sizes[i] / sizes[i-1]
        scf_time_ratio = scf_times[i] / scf_times[i-1]
        fbp_time_ratio = fbp_times[i] / fbp_times[i-1]
        
        print(f"\n{sizes[i-1]} → {sizes[i]} 水分子:")
        print(f"  大小比: {size_ratio:.1f}x")
        print(f"  理论时间比: {size_ratio**2:.1f}x")
        print(f"  SCF实际时间比: {scf_time_ratio:.1f}x")
        print(f"  FBP实际时间比: {fbp_time_ratio:.1f}x")

if __name__ == "__main__":
    test_128_water_system()
    test_scaling_behavior()