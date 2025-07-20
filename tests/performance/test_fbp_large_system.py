#!/usr/bin/env python3
"""
测试FBP在大体系（128水分子）中的性能和精度
"""

import pygcmc
import numpy as np
import time

def create_water_box(n_waters=128, density=1.0):
    """
    创建水分子盒子
    density: g/cm³，水的密度
    """
    print(f"创建{n_waters}个水分子系统...")
    
    # 计算盒子大小
    # 水的分子量 = 18.015 g/mol
    # 1个水分子的质量 = 18.015 / 6.022e23 g
    # n个水分子的质量 = n * 18.015 / 6.022e23 g
    # 体积 = 质量 / 密度
    # 体积 = n * 18.015 / (6.022e23 * density) cm³
    # 转换为nm³: 1 cm = 1e7 nm
    mass_per_molecule = 18.015 / 6.022e23  # g
    total_mass = n_waters * mass_per_molecule  # g
    volume_cm3 = total_mass / density  # cm³
    volume_nm3 = volume_cm3 * 1e21  # nm³
    box_length = np.cbrt(volume_nm3)  # nm
    
    print(f"盒子尺寸: {box_length:.2f} nm")
    
    atoms = []
    residues = []
    
    # 在立方体中随机放置水分子
    np.random.seed(42)  # 可重复性
    
    # 生成初始网格位置
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_length / n_per_side
    
    water_count = 0
    positions_used = []
    
    # 网格位置 + 小扰动
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础网格位置
                x_base = (i + 0.5) * spacing
                y_base = (j + 0.5) * spacing
                z_base = (k + 0.5) * spacing
                
                # 添加小的随机扰动
                x_base += (np.random.rand() - 0.5) * spacing * 0.2
                y_base += (np.random.rand() - 0.5) * spacing * 0.2
                z_base += (np.random.rand() - 0.5) * spacing * 0.2
                
                # 随机旋转
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                
                # SWM4-NDP水模型原子位置（相对坐标）
                rel_positions = [
                    (0.0, 0.0, 0.0, 1.71636, 0),   # O
                    (0.0, 0.0, 0.0, -1.71636, 1),  # D
                    (0.09572, 0.0, 0.0, 0.55733, 2),  # H1
                    (-0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
                    (0.0, -0.024034, 0.0, -1.11466, 3)  # M
                ]
                
                # 应用旋转并添加原子
                for dx, dy, dz, charge, typ in rel_positions:
                    # 简单旋转（绕z轴）
                    x_rot = dx * np.cos(theta) - dy * np.sin(theta)
                    y_rot = dx * np.sin(theta) + dy * np.cos(theta)
                    z_rot = dz
                    
                    atom = pygcmc.MCAtom()
                    atom.x = x_base + x_rot
                    atom.y = y_base + y_rot
                    atom.z = z_base + z_rot
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
    state.info.cutoff = min(box_length/2 - 0.1, 1.2)  # 典型截断距离
    
    # 力场参数（SWM4-NDP）
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]  # 只有O有LJ
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print(f"系统创建完成: {n_waters}个水分子, {len(atoms)}个原子")
    print(f"截断距离: {state.info.cutoff:.2f} nm")
    
    return state

def test_algorithm_performance(state, n_waters, algorithm, tolerance, max_iterations=1000):
    """测试单个算法的性能"""
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
    print(f"添加Thole屏蔽对...")
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = max_iterations
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 重置Drude位置到parent位置
    state_copy = state.copy()
    for i in range(n_waters):
        state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
        state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
        state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
    
    # 计时开始
    start_time = time.time()
    
    # 计算能量
    energy = force.calculateEnergySCF(state_copy)
    
    # 计时结束
    end_time = time.time()
    elapsed_time = (end_time - start_time) * 1000  # ms
    
    # 计算Drude位移统计
    displacements = []
    for i in range(n_waters):
        dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
        dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
        dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz)
        displacements.append(disp)
    
    displacements = np.array(displacements) * 1000  # nm to pm
    
    return {
        'energy': energy,
        'time': elapsed_time,
        'avg_displacement': np.mean(displacements),
        'max_displacement': np.max(displacements),
        'min_displacement': np.min(displacements),
        'std_displacement': np.std(displacements)
    }

def run_large_system_test():
    """运行大体系测试"""
    print("FBP大体系性能测试")
    print("="*80)
    
    # 测试不同大小的系统
    n_waters_list = [32, 64, 128]
    
    for n_waters in n_waters_list:
        print(f"\n\n{'='*80}")
        print(f"测试 {n_waters} 水分子系统")
        print(f"{'='*80}")
        
        # 创建系统
        state = create_water_box(n_waters)
        
        # 测试配置
        test_configs = [
            ("SCF (高精度)", "SCF", 0.01),
            ("SCF (中精度)", "SCF", 0.1),
            ("SCF (低精度)", "SCF", 1.0),
            ("FBP (标准)", "FBP", 1.0),
            ("FBP (宽松)", "FBP", 10.0)
        ]
        
        results = []
        ref_energy = None
        
        print(f"\n运行测试...")
        for name, algo, tol in test_configs:
            print(f"\n测试 {name}...")
            try:
                result = test_algorithm_performance(state, n_waters, algo, tol)
                result['name'] = name
                results.append(result)
                
                # 保存第一个SCF结果作为参考
                if ref_energy is None and algo == "SCF":
                    ref_energy = result['energy']
                    
            except Exception as e:
                print(f"  错误: {e}")
                results.append({
                    'name': name,
                    'energy': float('nan'),
                    'time': float('nan'),
                    'avg_displacement': float('nan')
                })
        
        # 打印结果表格
        print(f"\n\n结果汇总 ({n_waters}水分子):")
        print("-"*100)
        print(f"{'算法':<20} {'能量(kJ/mol)':<15} {'能量误差':<12} {'误差(%)':<10} {'时间(ms)':<12} {'速度提升':<10} {'平均位移(pm)':<12}")
        print("-"*100)
        
        # 找到SCF时间作为基准
        scf_time = None
        for r in results:
            if r['name'].startswith("SCF") and not np.isnan(r['time']):
                scf_time = r['time']
                break
        
        for r in results:
            if not np.isnan(r['energy']):
                energy_error = abs(r['energy'] - ref_energy) if ref_energy else 0
                error_percent = (energy_error / abs(ref_energy) * 100) if ref_energy and abs(ref_energy) > 0.01 else 0
                speedup = scf_time / r['time'] if scf_time and r['time'] > 0 else 1.0
                
                print(f"{r['name']:<20} {r['energy']:<15.6f} {energy_error:<12.6f} {error_percent:<10.3f} {r['time']:<12.1f} {speedup:<10.1f}x {r['avg_displacement']:<12.3f}")
            else:
                print(f"{r['name']:<20} {'失败':<15} {'N/A':<12} {'N/A':<10} {'N/A':<12} {'N/A':<10} {'N/A':<12}")
        
        # 详细的位移统计
        print(f"\n\nDrude位移统计 ({n_waters}水分子):")
        print("-"*80)
        print(f"{'算法':<20} {'平均(pm)':<12} {'最大(pm)':<12} {'最小(pm)':<12} {'标准差(pm)':<12}")
        print("-"*80)
        
        for r in results:
            if not np.isnan(r['energy']):
                print(f"{r['name']:<20} {r['avg_displacement']:<12.3f} {r.get('max_displacement', 0):<12.3f} {r.get('min_displacement', 0):<12.3f} {r.get('std_displacement', 0):<12.3f}")

def analyze_scaling():
    """分析算法的扩展性"""
    print("\n\n\n算法扩展性分析")
    print("="*80)
    
    sizes = [16, 32, 64, 128]
    scf_times = []
    fbp_times = []
    
    for n in sizes:
        print(f"\n测试 {n} 水分子...")
        state = create_water_box(n)
        
        # SCF
        scf_result = test_algorithm_performance(state, n, "SCF", 0.1)
        scf_times.append(scf_result['time'])
        
        # FBP
        fbp_result = test_algorithm_performance(state, n, "FBP", 1.0)
        fbp_times.append(fbp_result['time'])
    
    print("\n\n扩展性结果:")
    print("-"*60)
    print(f"{'系统大小':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'加速比':<10}")
    print("-"*60)
    
    for i, n in enumerate(sizes):
        speedup = scf_times[i] / fbp_times[i] if fbp_times[i] > 0 else 0
        print(f"{n:<10} {scf_times[i]:<15.1f} {fbp_times[i]:<15.1f} {speedup:<10.1f}x")
    
    # 计算扩展性
    print("\n\n扩展性分析:")
    print("-"*60)
    print("理论上，计算时间应该与N²成正比（N是水分子数）")
    
    # 计算时间比例
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
    run_large_system_test()
    analyze_scaling()