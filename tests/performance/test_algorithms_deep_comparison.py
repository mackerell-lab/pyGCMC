#!/usr/bin/env python3
"""
深入测试和比较 SCF, OPT3, OPT3+SCF, FBP 四种算法
"""

import pygcmc
import numpy as np
import time
from collections import defaultdict

def create_water_system(n_waters):
    """创建水分子系统"""
    atoms = []
    residues = []
    
    # 网格放置
    spacing = 0.35  # nm
    n_per_side = int(np.ceil(n_waters**(1/3)))
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                base_x = i * spacing + 1.0
                base_y = j * spacing + 1.0
                base_z = k * spacing + 1.0
                
                # 添加随机扰动使系统更真实
                base_x += np.random.uniform(-0.05, 0.05)
                base_y += np.random.uniform(-0.05, 0.05)
                base_z += np.random.uniform(-0.05, 0.05)
                
                # 水分子原子 (O, D, H, H, M)
                positions = [
                    (base_x, base_y, base_z, 1.71636, 0),   # O
                    (base_x, base_y, base_z, -1.71636, 1),  # D
                    (base_x + 0.09572, base_y, base_z, 0.55733, 2),  # H1
                    (base_x - 0.04786, base_y + 0.08288, base_z, 0.55733, 2),  # H2
                    (base_x, base_y - 0.024034, base_z, -1.11466, 3)  # M
                ]
                
                for x, y, z, charge, typ in positions:
                    atom = pygcmc.MCAtom()
                    atom.x = x
                    atom.y = y
                    atom.z = z
                    atom.charge = charge
                    atom.type = typ
                    atoms.append(atom)
                
                # 残基
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
    
    # 创建状态
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # 盒子信息
    box_size = (n_per_side + 1) * spacing + 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(4.5, box_size/2.0 - 0.1)
    
    # 力场
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_drude_force(n_waters):
    """创建Drude力对象"""
    force = pygcmc.DrudeForce()
    
    # Drude参数
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    # 添加粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge,
            polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # 添加屏蔽对
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force.addScreenedPair(i, j, 1.3)
    
    return force

def test_algorithm(state, n_waters, algorithm, tolerance, max_iter=50):
    """测试单个算法"""
    force = create_drude_force(n_waters)
    
    # 设置参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tolerance
    params.maxIterations = max_iter
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    
    # 设置算法
    if algorithm == "SCF":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    elif algorithm == "OPT3":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    elif algorithm == "OPT3+SCF":
        # 这需要特殊处理
        force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
    
    # 重置Drude位置
    for i in range(n_waters):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    # 计算能量
    start = time.time()
    energy = force.calculateEnergySCF(state)
    elapsed = time.time() - start
    
    # 保存Drude位置用于分析
    drude_positions = []
    for i in range(n_waters):
        idx = 5*i + 1
        drude_positions.append((
            state.atoms[idx].x,
            state.atoms[idx].y,
            state.atoms[idx].z
        ))
    
    return energy, elapsed, drude_positions

def calculate_position_rmsd(pos1, pos2):
    """计算两组位置的RMSD"""
    rmsd = 0.0
    for p1, p2 in zip(pos1, pos2):
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dz = p1[2] - p2[2]
        rmsd += dx*dx + dy*dy + dz*dz
    return np.sqrt(rmsd / len(pos1))

def main():
    """主测试函数"""
    print("深入比较 SCF, OPT3, OPT3+SCF, FBP 算法")
    print("="*80)
    
    # 设置随机种子
    np.random.seed(42)
    
    # 测试不同大小的系统
    system_sizes = [5, 10, 20, 40]
    
    # 不同的容差设置
    tolerances = [10.0, 5.0, 2.0, 1.0, 0.5, 0.1, 0.01]
    
    results = defaultdict(lambda: defaultdict(list))
    
    for n_waters in system_sizes:
        print(f"\n测试 {n_waters} 水分子系统")
        print("-"*60)
        
        # 创建系统
        state = create_water_system(n_waters)
        
        # 获取高精度参考能量 (SCF with very tight tolerance)
        print("计算参考能量 (SCF, tol=0.001)...")
        ref_energy, ref_time, ref_positions = test_algorithm(
            state, n_waters, "SCF", 0.001, max_iter=1000
        )
        print(f"参考能量: {ref_energy:.6f} kJ/mol (耗时: {ref_time*1000:.2f} ms)")
        
        # 测试每种算法在不同容差下的表现
        print("\n不同容差下的算法比较:")
        print(f"{'Tolerance':<10} {'Algorithm':<12} {'Energy':<15} {'Error':<12} {'RMSD (pm)':<12} {'Time (ms)':<10}")
        print("-"*80)
        
        for tol in tolerances:
            for algo in ["SCF", "OPT3", "FBP"]:  # OPT3+SCF暂时跳过
                energy, elapsed, positions = test_algorithm(state, n_waters, algo, tol)
                
                # 计算误差
                energy_error = abs(energy - ref_energy)
                position_rmsd = calculate_position_rmsd(positions, ref_positions) * 1000  # nm to pm
                
                # 保存结果
                results[n_waters][algo].append({
                    'tolerance': tol,
                    'energy': energy,
                    'error': energy_error,
                    'rmsd': position_rmsd,
                    'time': elapsed * 1000
                })
                
                print(f"{tol:<10.2f} {algo:<12} {energy:<15.6f} {energy_error:<12.6f} {position_rmsd:<12.3f} {elapsed*1000:<10.2f}")
    
    # 汇总分析
    print("\n\n" + "="*80)
    print("汇总分析")
    print("="*80)
    
    # 找出每种算法达到特定精度所需的最小容差
    print("\n达到不同精度所需的最小容差:")
    print(f"{'System':<10} {'Algorithm':<12} {'E<0.1 kJ/mol':<15} {'E<0.01 kJ/mol':<15} {'RMSD<1pm':<15}")
    print("-"*70)
    
    for n_waters in system_sizes:
        for algo in ["SCF", "OPT3", "FBP"]:
            tol_01 = None
            tol_001 = None
            tol_1pm = None
            
            for result in results[n_waters][algo]:
                if result['error'] < 0.1 and tol_01 is None:
                    tol_01 = result['tolerance']
                if result['error'] < 0.01 and tol_001 is None:
                    tol_001 = result['tolerance']
                if result['rmsd'] < 1.0 and tol_1pm is None:
                    tol_1pm = result['tolerance']
            
            tol_01_str = f"{tol_01:.2f}" if tol_01 else "Not achieved"
            tol_001_str = f"{tol_001:.2f}" if tol_001 else "Not achieved"
            tol_1pm_str = f"{tol_1pm:.2f}" if tol_1pm else "Not achieved"
            
            print(f"{n_waters:<10} {algo:<12} {tol_01_str:<15} {tol_001_str:<15} {tol_1pm_str:<15}")
    
    # 速度比较
    print("\n\n速度比较 (相对于SCF):")
    print(f"{'System':<10} {'Tolerance':<12} {'OPT3 Speedup':<15} {'FBP Speedup':<15}")
    print("-"*55)
    
    for n_waters in system_sizes:
        for i, tol in enumerate([1.0, 0.5, 0.1]):
            if i < len(results[n_waters]["SCF"]):
                scf_time = results[n_waters]["SCF"][i]['time']
                opt3_time = results[n_waters]["OPT3"][i]['time']
                fbp_time = results[n_waters]["FBP"][i]['time']
                
                opt3_speedup = scf_time / opt3_time
                fbp_speedup = scf_time / fbp_time
                
                print(f"{n_waters:<10} {tol:<12.1f} {opt3_speedup:<15.1f} {fbp_speedup:<15.1f}")
    
    # 精度-速度权衡分析
    print("\n\n精度-速度权衡分析 (20水系统):")
    print(f"{'Algorithm':<12} {'Best Time@0.1kJ':<20} {'Best Time@0.01kJ':<20}")
    print("-"*55)
    
    for algo in ["SCF", "OPT3", "FBP"]:
        time_01 = None
        time_001 = None
        
        for result in results[20][algo]:
            if result['error'] < 0.1 and (time_01 is None or result['time'] < time_01):
                time_01 = result['time']
            if result['error'] < 0.01 and (time_001 is None or result['time'] < time_001):
                time_001 = result['time']
        
        time_01_str = f"{time_01:.2f} ms" if time_01 else "Not achieved"
        time_001_str = f"{time_001:.2f} ms" if time_001 else "Not achieved"
        
        print(f"{algo:<12} {time_01_str:<20} {time_001_str:<20}")

if __name__ == "__main__":
    main()