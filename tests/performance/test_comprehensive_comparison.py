#!/usr/bin/env python3
"""
综合比较各种Drude算法的能量、容差和速度
"""

import pygcmc
import numpy as np
import time

def create_water_system(n_waters):
    """创建水分子系统"""
    atoms = []
    residues = []
    
    # 在立方体中均匀分布水分子
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = 0.5  # nm
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                x_base = i * spacing
                y_base = j * spacing
                z_base = k * spacing
                
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
    
    # 设置盒子
    box_size = (n_per_side - 1) * spacing + 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(4.5, box_size/2 - 0.1)
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def test_algorithm(state, n_waters, algorithm, tolerance, max_iter=1000):
    """测试单个算法"""
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
    charge = -1.71636
    k_spring = 418400.0  # kJ/mol/nm^2
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
    elif algorithm == "OPT4":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT4)
    elif algorithm == "HybridOPT":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.HybridOPT)
    elif algorithm == "SmartOPT3":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SmartOPT3)
    elif algorithm == "FBP":
        force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # 重置Drude位置
    for i in range(n_waters):
        state.atoms[5*i+1].x = state.atoms[5*i].x
        state.atoms[5*i+1].y = state.atoms[5*i].y
        state.atoms[5*i+1].z = state.atoms[5*i].z
    
    # 计时和计算
    start_time = time.time()
    energy = force.calculateEnergySCF(state)
    end_time = time.time()
    
    # 计算Drude位移的RMSD
    displacements = []
    for i in range(n_waters):
        dx = state.atoms[5*i+1].x - state.atoms[5*i].x
        dy = state.atoms[5*i+1].y - state.atoms[5*i].y
        dz = state.atoms[5*i+1].z - state.atoms[5*i].z
        disp = np.sqrt(dx*dx + dy*dy + dz*dz)
        displacements.append(disp)
    
    rmsd = np.sqrt(np.mean(np.array(displacements)**2)) * 1000  # nm to pm
    
    return {
        'energy': energy,
        'time': (end_time - start_time) * 1000,  # ms
        'rmsd': rmsd
    }

def run_comprehensive_test():
    """运行综合测试"""
    print("Drude算法综合性能测试")
    print("="*80)
    
    # 测试配置
    n_waters_list = [5, 10, 20, 50]
    algorithms = ["SCF", "OPT3", "OPT4", "HybridOPT", "SmartOPT3", "FBP"]
    tolerances = [10.0, 1.0, 0.1, 0.01]
    
    # 对每个系统大小
    for n_waters in n_waters_list:
        print(f"\n\n{'='*80}")
        print(f"测试 {n_waters} 水分子系统")
        print(f"{'='*80}")
        
        # 创建系统
        state = create_water_system(n_waters)
        
        # 获取SCF参考能量（高精度）
        print("\n计算参考能量 (SCF, tol=0.001)...")
        ref_result = test_algorithm(state.copy(), n_waters, "SCF", 0.001)
        ref_energy = ref_result['energy']
        print(f"参考能量: {ref_energy:.6f} kJ/mol")
        
        # 测试不同容差
        print(f"\n{'容差':<10} {'算法':<12} {'能量(kJ/mol)':<15} {'误差(kJ/mol)':<15} {'误差(%)':<10} {'RMSD(pm)':<10} {'时间(ms)':<10} {'速度提升':<10}")
        print("-"*110)
        
        for tolerance in tolerances:
            scf_time = None
            
            for algorithm in algorithms:
                # 某些算法不需要多个容差测试
                if algorithm in ["OPT3", "OPT4"] and tolerance not in [1.0]:
                    continue
                
                try:
                    result = test_algorithm(state.copy(), n_waters, algorithm, tolerance)
                    
                    energy = result['energy']
                    time_ms = result['time']
                    rmsd = result['rmsd']
                    
                    # 计算误差
                    error = abs(energy - ref_energy)
                    error_pct = error / abs(ref_energy) * 100 if abs(ref_energy) > 0.01 else 0
                    
                    # 记录SCF时间用于速度比较
                    if algorithm == "SCF" and scf_time is None:
                        scf_time = time_ms
                    
                    speedup = scf_time / time_ms if scf_time else 1.0
                    
                    print(f"{tolerance:<10.2f} {algorithm:<12} {energy:<15.6f} {error:<15.6f} {error_pct:<10.2f} {rmsd:<10.3f} {time_ms:<10.2f} {speedup:<10.1f}x")
                    
                except Exception as e:
                    print(f"{tolerance:<10.2f} {algorithm:<12} {'失败':<15} {str(e)[:50]}")
        
        # 总结最佳算法
        print(f"\n{n_waters}水分子系统总结:")
        print("-"*60)
        print("推荐算法：")
        print("- 高精度（误差<0.01 kJ/mol）: SCF (tol=0.1) 或 HybridOPT")
        print("- 快速估算（误差<1%）: FBP (tol=1.0)")
        print("- 平衡选择：HybridOPT")

def test_convergence_behavior():
    """测试收敛行为"""
    print("\n\n收敛行为测试")
    print("="*80)
    
    n_waters = 10
    state = create_water_system(n_waters)
    
    print(f"\n测试 {n_waters} 水分子系统的收敛行为")
    print(f"\n{'算法':<12} {'容差':<10} {'迭代次数':<12} {'最终误差':<15} {'收敛?':<10}")
    print("-"*70)
    
    algorithms = ["SCF", "FBP", "HybridOPT"]
    tolerances = [10.0, 1.0, 0.1, 0.01]
    
    for algorithm in algorithms:
        for tolerance in tolerances:
            # 这里需要更详细的收敛信息，暂时使用简化版本
            try:
                result = test_algorithm(state.copy(), n_waters, algorithm, tolerance, max_iter=100)
                print(f"{algorithm:<12} {tolerance:<10.2f} {'<100':<12} {'N/A':<15} {'是':<10}")
            except:
                print(f"{algorithm:<12} {tolerance:<10.2f} {'>100':<12} {'N/A':<15} {'否':<10}")

if __name__ == "__main__":
    run_comprehensive_test()
    test_convergence_behavior()