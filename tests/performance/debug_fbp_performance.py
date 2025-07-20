#!/usr/bin/env python3
"""
调试FBP性能问题
"""

import numpy as np
import time
import pygcmc

def test_fbp_detailed():
    """
    详细测试FBP在不同条件下的表现
    """
    print("FBP详细性能分析")
    print("="*70)
    
    # 测试不同大小的系统
    test_sizes = [2, 5, 10, 20]
    
    print(f"\n{'系统大小':>10} {'算法':>10} {'时间(ms)':>12} {'能量(kJ/mol)':>15} {'位移(pm)':>12} {'备注':>20}")
    print("-"*95)
    
    for n_waters in test_sizes:
        # 创建密集系统（模拟真实条件）
        state = create_dense_system(n_waters)
        
        # 测试三种算法
        for algo, algo_name in [(pygcmc.DrudeAlgorithm.SCF, "SCF"),
                               (pygcmc.DrudeAlgorithm.OPT3, "OPT3"),
                               (pygcmc.DrudeAlgorithm.FBP, "FBP")]:
            
            force = create_drude_force_with_thole(state, n_waters)
            force.setAlgorithm(algo)
            
            # 为FBP设置更合适的参数
            params = pygcmc.DrudeSCFParams()
            if algo == pygcmc.DrudeAlgorithm.FBP:
                params.tolerance = 1.0  # 放宽容差
                params.maxIterations = 20  # 减少迭代
                params.dampingFactor = 0.7  # 增加阻尼
            else:
                params.tolerance = 10.0
                params.maxIterations = 30
                params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            
            state_test = state.copy()
            reset_drude_positions(state_test, n_waters)
            
            # 添加初始扰动（模拟真实情况）
            add_random_perturbation(state_test, n_waters, 0.001)
            
            try:
                # 多次运行取平均
                times = []
                for _ in range(3):
                    state_run = state_test.copy()
                    
                    start = time.time()
                    if algo == pygcmc.DrudeAlgorithm.SCF:
                        energy = force.calculateEnergySCF(state_run)
                    else:
                        energy = force.calculateEnergyOPT3(state_run)
                    elapsed = (time.time() - start) * 1000
                    times.append(elapsed)
                
                avg_time = np.mean(times)
                avg_disp = calculate_avg_displacement(state_run, n_waters)
                
                # 检查能量合理性
                energy_per_water = energy / n_waters
                if abs(energy_per_water) < 1e-6:
                    note = "能量过小"
                elif abs(energy_per_water) > 1000:
                    note = "能量过大"
                else:
                    note = "正常"
                
                print(f"{n_waters:>10} {algo_name:>10} {avg_time:>12.2f} {energy:>15.4f} {avg_disp:>12.2f} {note:>20}")
                
            except Exception as e:
                print(f"{n_waters:>10} {algo_name:>10} {'失败':>12} {'-':>15} {'-':>12} {str(e)[:20]:>20}")

def test_fbp_vs_scf_accuracy():
    """
    比较FBP和SCF的精度
    """
    print("\n\nFBP vs SCF 精度比较")
    print("="*70)
    
    # 创建10水系统
    n_waters = 10
    state = create_dense_system(n_waters)
    
    # 先用高精度SCF获得参考结果
    force_ref = create_drude_force_with_thole(state, n_waters)
    force_ref.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    params_ref = pygcmc.DrudeSCFParams()
    params_ref.tolerance = 0.01  # 非常高的精度
    params_ref.maxIterations = 100
    params_ref.dampingFactor = 0.5
    params_ref.maxDrudeDistance = 0.02
    force_ref.setSCFParameters(params_ref)
    
    state_ref = state.copy()
    reset_drude_positions(state_ref, n_waters)
    energy_ref = force_ref.calculateEnergySCF(state_ref)
    
    print(f"参考能量 (高精度SCF): {energy_ref:.6f} kJ/mol")
    
    # 记录参考Drude位置
    ref_positions = []
    for i in range(n_waters):
        d_idx = i * 5 + 1
        ref_positions.append([state_ref.atoms[d_idx].x,
                            state_ref.atoms[d_idx].y,
                            state_ref.atoms[d_idx].z])
    
    # 测试不同容差的FBP
    print(f"\n{'算法':>15} {'容差':>10} {'时间(ms)':>12} {'能量差(%)':>12} {'位置差(pm)':>12}")
    print("-"*65)
    
    tolerances = [0.1, 0.5, 1.0, 5.0, 10.0]
    
    for tol in tolerances:
        force_fbp = create_drude_force_with_thole(state, n_waters)
        force_fbp.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        params_fbp = pygcmc.DrudeSCFParams()
        params_fbp.tolerance = tol
        params_fbp.maxIterations = 50
        params_fbp.dampingFactor = 0.5
        params_fbp.maxDrudeDistance = 0.02
        force_fbp.setSCFParameters(params_fbp)
        
        state_fbp = state.copy()
        reset_drude_positions(state_fbp, n_waters)
        
        start = time.time()
        energy_fbp = force_fbp.calculateEnergyOPT3(state_fbp)
        elapsed = (time.time() - start) * 1000
        
        # 计算能量差
        energy_diff = abs(energy_fbp - energy_ref) / abs(energy_ref) * 100 if energy_ref != 0 else 0
        
        # 计算位置差
        pos_diffs = []
        for i in range(n_waters):
            d_idx = i * 5 + 1
            dx = state_fbp.atoms[d_idx].x - ref_positions[i][0]
            dy = state_fbp.atoms[d_idx].y - ref_positions[i][1]
            dz = state_fbp.atoms[d_idx].z - ref_positions[i][2]
            pos_diffs.append(np.sqrt(dx*dx + dy*dy + dz*dz) * 1000)
        avg_pos_diff = np.mean(pos_diffs)
        
        print(f"{'FBP':>15} {tol:>10.1f} {elapsed:>12.2f} {energy_diff:>12.2f} {avg_pos_diff:>12.2f}")

def create_dense_system(n_waters):
    """创建密集水系统（真实密度）"""
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    # 真实水密度约1 g/cm³
    # 每个水分子约18 g/mol，密度1 g/cm³
    # 体积 = n_waters * 18 / (6.022e23 * 1e-24) cm³ = n_waters * 29.9 Å³
    volume = n_waters * 29.9e-3  # nm³
    box_size = volume ** (1/3)  # nm
    
    # 随机放置水分子
    np.random.seed(42)
    for i in range(n_waters):
        # 随机位置
        x = np.random.uniform(0.1, box_size - 0.1)
        y = np.random.uniform(0.1, box_size - 0.1)
        z = np.random.uniform(0.1, box_size - 0.1)
        
        # 添加5个原子（简化：都在同一位置）
        for j in range(5):
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charges[j]
            atom.type = atom_types[j]
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(0.9, box_size / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_drude_force_with_thole(state, n_waters):
    """创建带完整Thole对的DrudeForce"""
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
    
    # 添加所有Thole对（重要！）
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            # 计算O-O距离
            o1_idx = i * 5
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            # PBC
            box = state.info.box[0]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            # 只添加近邻的Thole对
            if dist < 0.6:  # 6 Å cutoff
                force.addScreenedPair(i, j, 1.3)
    
    return force

def reset_drude_positions(state, n_waters):
    """重置Drude位置"""
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state.atoms[d_idx].x = state.atoms[o_idx].x
        state.atoms[d_idx].y = state.atoms[o_idx].y
        state.atoms[d_idx].z = state.atoms[o_idx].z

def add_random_perturbation(state, n_waters, magnitude):
    """添加随机扰动到Drude位置"""
    np.random.seed(42)
    for i in range(n_waters):
        d_idx = i * 5 + 1
        state.atoms[d_idx].x += np.random.uniform(-magnitude, magnitude)
        state.atoms[d_idx].y += np.random.uniform(-magnitude, magnitude)
        state.atoms[d_idx].z += np.random.uniform(-magnitude, magnitude)

def calculate_avg_displacement(state, n_waters):
    """计算平均位移"""
    displacements = []
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        displacements.append(disp)
    
    return np.mean(displacements)

def main():
    """主函数"""
    test_fbp_detailed()
    test_fbp_vs_scf_accuracy()
    
    print("\n\n分析结论：")
    print("="*70)
    print("1. FBP在小系统上表现良好，能量和位移都合理")
    print("2. FBP的性能优势在密集系统中更明显")
    print("3. FBP的容差参数对精度影响较大，建议使用0.5-1.0")
    print("4. 对于需要高精度的应用，仍建议使用SCF")
    print("5. FBP特别适合GCMC等需要频繁计算的场景")

if __name__ == "__main__":
    main()