#!/usr/bin/env python3
"""
比较不同算法的收敛性和能量精度
"""

import pygcmc
import numpy as np
import time

def create_simple_water_system(n_waters=5):
    """创建简单的水系统"""
    atoms = []
    residues = []
    
    # 线性排列，避免复杂相互作用
    spacing = 0.5  # nm
    
    for i in range(n_waters):
        base_x = i * spacing
        base_y = 0.0
        base_z = 0.0
        
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
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    # 创建状态
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # 盒子信息
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    # 力场
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_single_algorithm(state, n_waters, algorithm, tolerance, max_iter=50):
    """测试单个算法"""
    # 创建力对象
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
    elif algorithm == "Smart OPT3":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SmartOPT3)
    
    # 重置Drude位置
    for i in range(n_waters):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    # 计算能量
    start = time.time()
    energy = force.calculateEnergySCF(state)
    elapsed = time.time() - start
    
    # 计算Drude位移
    displacements = []
    for i in range(n_waters):
        drude_idx = 5*i + 1
        parent_idx = 5*i
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz)
        displacements.append(disp)
    
    avg_disp = np.mean(displacements) * 1000  # nm to pm
    max_disp = np.max(displacements) * 1000
    
    return energy, elapsed, avg_disp, max_disp

def main():
    """主测试函数"""
    print("算法收敛性比较")
    print("="*100)
    
    # 测试5个水分子
    n_waters = 5
    state = create_simple_water_system(n_waters)
    
    # 首先获取高精度参考
    print(f"\n{n_waters} 水分子系统:")
    print("-"*100)
    
    print("计算高精度参考 (SCF, tol=0.0001, maxIter=2000)...")
    ref_energy, ref_time, ref_avg_disp, ref_max_disp = test_single_algorithm(
        state, n_waters, "SCF", 0.0001, max_iter=2000
    )
    print(f"参考能量: {ref_energy:.8f} kJ/mol")
    print(f"平均Drude位移: {ref_avg_disp:.3f} pm")
    print(f"最大Drude位移: {ref_max_disp:.3f} pm")
    print(f"计算时间: {ref_time*1000:.2f} ms")
    
    # 测试不同算法和容差
    print("\n不同算法和容差的比较:")
    print(f"{'Algorithm':<12} {'Tolerance':<10} {'Energy':<15} {'Error':<12} {'Avg Disp (pm)':<15} {'Time (ms)':<10} {'Speedup':<10}")
    print("-"*100)
    
    algorithms = ["SCF", "OPT3", "Smart OPT3", "FBP"]
    tolerances = [10.0, 5.0, 2.0, 1.0, 0.5, 0.1]
    
    results = {}
    
    for algo in algorithms:
        results[algo] = []
        for tol in tolerances:
            energy, elapsed, avg_disp, max_disp = test_single_algorithm(
                state, n_waters, algo, tol
            )
            
            error = abs(energy - ref_energy)
            speedup = ref_time / elapsed
            
            results[algo].append({
                'tolerance': tol,
                'energy': energy,
                'error': error,
                'avg_disp': avg_disp,
                'time': elapsed * 1000,
                'speedup': speedup
            })
            
            print(f"{algo:<12} {tol:<10.2f} {energy:<15.8f} {error:<12.8f} {avg_disp:<15.3f} {elapsed*1000:<10.2f} {speedup:<10.1f}")
    
    # 分析能量精度
    print("\n\n能量精度分析:")
    print("="*80)
    print(f"{'Algorithm':<12} {'Min Error':<15} {'at Tolerance':<15} {'Max Error':<15}")
    print("-"*80)
    
    for algo in algorithms:
        errors = [r['error'] for r in results[algo]]
        min_error = min(errors)
        min_idx = errors.index(min_error)
        min_tol = results[algo][min_idx]['tolerance']
        max_error = max(errors)
        
        print(f"{algo:<12} {min_error:<15.8f} {min_tol:<15.2f} {max_error:<15.8f}")
    
    # 收敛性分析
    print("\n\n收敛性分析 (能量误差 < 0.01 kJ/mol):")
    print("="*60)
    print(f"{'Algorithm':<12} {'Min Tolerance':<15} {'Time (ms)':<15}")
    print("-"*60)
    
    for algo in algorithms:
        converged_tol = None
        converged_time = None
        
        for r in results[algo]:
            if r['error'] < 0.01:
                converged_tol = r['tolerance']
                converged_time = r['time']
                break
        
        if converged_tol:
            print(f"{algo:<12} {converged_tol:<15.2f} {converged_time:<15.2f}")
        else:
            print(f"{algo:<12} {'Not achieved':<15} {'-':<15}")
    
    # 测试OPT3+SCF混合方法
    print("\n\n测试OPT3+SCF混合方法:")
    print("="*60)
    
    # 先用OPT3预测
    force_opt3 = pygcmc.DrudeForce()
    for i in range(n_waters):
        force_opt3.addParticle(
            drudeIndex=5*i + 1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=1.71636**2 * 138.935456 / 418400.0,
            aniso12=1.0, aniso34=1.0
        )
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force_opt3.addScreenedPair(i, j, 1.3)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0
    params.maxIterations = 1
    force_opt3.setSCFParameters(params)
    force_opt3.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    
    # 重置位置
    for i in range(n_waters):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    start = time.time()
    energy_opt3 = force_opt3.calculateEnergySCF(state)
    time_opt3 = time.time() - start
    
    # 然后用SCF精修
    force_scf = pygcmc.DrudeForce()
    for i in range(n_waters):
        force_scf.addParticle(
            drudeIndex=5*i + 1, parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=1.71636**2 * 138.935456 / 418400.0,
            aniso12=1.0, aniso34=1.0
        )
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force_scf.addScreenedPair(i, j, 1.3)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1
    params.maxIterations = 50
    force_scf.setSCFParameters(params)
    force_scf.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 从OPT3的位置开始（不重置）
    start = time.time()
    energy_final = force_scf.calculateEnergySCF(state)
    time_scf = time.time() - start
    
    total_time = (time_opt3 + time_scf) * 1000
    error = abs(energy_final - ref_energy)
    speedup = ref_time / (time_opt3 + time_scf)
    
    print(f"OPT3 预测能量: {energy_opt3:.8f} kJ/mol")
    print(f"最终能量: {energy_final:.8f} kJ/mol")
    print(f"能量误差: {error:.8f} kJ/mol")
    print(f"总时间: {total_time:.2f} ms")
    print(f"加速比: {speedup:.1f}x")

if __name__ == "__main__":
    main()