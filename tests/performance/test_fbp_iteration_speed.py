#!/usr/bin/env python3
"""
测试FBP和SCF每次迭代的纯速度
不考虑收敛，只运行固定次数迭代
"""

import pygcmc
import numpy as np
import time

def create_test_system(n_waters):
    """创建测试系统"""
    # 使用较大间距
    spacing = 0.6  # nm
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
                
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                # SWM4-NDP
                positions = [
                    (x, y, z, 1.71636, 0),   # O
                    (x, y, z, -1.71636, 1),  # D
                    (x + 0.09572, y, z, 0.55733, 2),  # H1
                    (x - 0.04786, y + 0.08288, z, 0.55733, 2),  # H2
                    (x, y - 0.024034, z, -1.11466, 3)  # M
                ]
                
                for px, py, pz, charge, typ in positions:
                    atom = pygcmc.MCAtom()
                    atom.x = px
                    atom.y = py
                    atom.z = pz
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

def measure_iteration_speed(state, n_waters, algorithm, n_iterations):
    """测量固定迭代次数的速度"""
    force = pygcmc.DrudeForce()
    
    # SWM4-NDP参数
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
    
    # 设置参数 - 使用极大容差确保不会提前收敛
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e10  # 永不收敛
    params.maxIterations = n_iterations
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
    
    # 运行多次
    times = []
    for _ in range(5):
        start = time.time()
        force.calculateEnergySCF(state_copy)
        elapsed = (time.time() - start) * 1000
        times.append(elapsed)
    
    avg_time = np.mean(times)
    std_time = np.std(times)
    per_iter = avg_time / n_iterations
    
    return avg_time, std_time, per_iter

def main():
    """主测试函数"""
    print("FBP vs SCF 纯迭代速度测试")
    print("="*80)
    print("测试每次迭代的计算时间，不考虑收敛性")
    print("-"*80)
    
    # 测试不同系统大小
    system_sizes = [32, 64, 128]
    iteration_counts = [10, 50, 100]
    
    for n_waters in system_sizes:
        print(f"\n\n{n_waters}水分子系统")
        print("="*60)
        
        state = create_test_system(n_waters)
        print(f"盒子: {state.info.box[0]:.2f} nm, 截断: {state.info.cutoff:.2f} nm")
        
        print(f"\n{'迭代次数':<10} {'SCF总时间(ms)':<15} {'FBP总时间(ms)':<15} {'SCF每迭代(ms)':<15} {'FBP每迭代(ms)':<15} {'FBP/SCF比':<10}")
        print("-"*90)
        
        for n_iter in iteration_counts:
            scf_time, scf_std, scf_per = measure_iteration_speed(state, n_waters, "SCF", n_iter)
            fbp_time, fbp_std, fbp_per = measure_iteration_speed(state, n_waters, "FBP", n_iter)
            
            ratio = fbp_per / scf_per
            
            print(f"{n_iter:<10} {scf_time:<15.1f} {fbp_time:<15.1f} {scf_per:<15.3f} {fbp_per:<15.3f} {ratio:<10.2f}")
    
    # 分析计算复杂度
    print("\n\n计算复杂度分析")
    print("="*80)
    
    print("\n测试O(N²)扩展性（100次迭代）:")
    print(f"{'系统大小':<10} {'SCF时间(ms)':<15} {'FBP时间(ms)':<15} {'相对32水系统':<20}")
    print("-"*60)
    
    times_scf = {}
    times_fbp = {}
    
    for n_waters in [32, 64, 128]:
        state = create_test_system(n_waters)
        scf_time, _, _ = measure_iteration_speed(state, n_waters, "SCF", 100)
        fbp_time, _, _ = measure_iteration_speed(state, n_waters, "FBP", 100)
        
        times_scf[n_waters] = scf_time
        times_fbp[n_waters] = fbp_time
        
        if n_waters == 32:
            ref_scf = scf_time
            ref_fbp = fbp_time
            relative = "1.0x / 1.0x"
        else:
            relative = f"{scf_time/ref_scf:.1f}x / {fbp_time/ref_fbp:.1f}x"
        
        print(f"{n_waters:<10} {scf_time:<15.1f} {fbp_time:<15.1f} {relative:<20}")
    
    print("\n结论:")
    print("-"*60)
    
    # 计算平均每迭代时间比
    avg_ratio = 0
    count = 0
    for n_waters in system_sizes:
        state = create_test_system(n_waters)
        _, _, scf_per = measure_iteration_speed(state, n_waters, "SCF", 50)
        _, _, fbp_per = measure_iteration_speed(state, n_waters, "FBP", 50)
        ratio = fbp_per / scf_per
        avg_ratio += ratio
        count += 1
    
    avg_ratio /= count
    
    if avg_ratio > 1:
        print(f"FBP每次迭代比SCF慢约{avg_ratio:.1f}倍")
        print("这解释了为什么FBP在相同迭代次数下更慢")
    else:
        print(f"FBP每次迭代比SCF快约{1/avg_ratio:.1f}倍")
    
    print("\nFBP慢的原因:")
    print("1. 每次迭代需要计算完整的力（不只是电场）")
    print("2. 需要额外的力平衡计算和收敛检查")
    print("3. 在密集系统中收敛性较差，需要更多迭代")

if __name__ == "__main__":
    main()